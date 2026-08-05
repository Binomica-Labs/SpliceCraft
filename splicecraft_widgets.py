"""Self-contained Textual widget primitives (layer 3).

`_InstantPressButton` (fires Pressed on mouse-down, used by StrandPickerModal)
and the xterm-256 color grid (`_XtermColorGrid` + its `_xterm_index_to_hex` /
`_ANSI16_HEX` palette helpers, used by ColorPickerModal). Pure Textual/Rich
widgets -- no dependency on hub state, persistence, or any other splicecraft
module.
"""
from __future__ import annotations

import functools as _functools
import re

from rich.text import Text
from textual.coordinate import Coordinate as _Coordinate
from textual.events import MouseDown
from textual.widgets import Button, DataTable, DirectoryTree, Input, Static

from splicecraft_logging import _log
from splicecraft_util import (
    _DEFAULT_TYPE_COLORS, _FASTA_EXTS, _is_fasta_path, _is_seq_zip_path, _natural_sort_key,
)
from splicecraft_dataaccess import _load_feature_colors


class _InstantPressButton(Button):
    """Button that fires `Pressed` on mouse-DOWN rather than
    waiting for the Click cycle to complete.

    Sweep #31 (2026-05-26): works around the Textual real-terminal
    bug where a click on a non-focused widget gets eaten by
    focus-transition, and `Button.Pressed` only fires on the
    SECOND press (the first lands during focus-shift, the second
    actually triggers Click). Pilot's synthetic clicks bypass
    that gate so the bug never surfaces in tests — but it does
    on a real WSL2 / xterm session. Posting `Pressed` on
    mouse-down + stopping the event means a single physical
    click always wins, regardless of which widget held focus
    beforehand.

    Used by `StrandPickerModal` so picking ▶ / ◀ / ▒ / ↔ is
    one-and-done. Drop-in replacement for `Button`; @on(Button
    .Pressed, ...) handlers in the owning modal still see the
    same `Pressed` message, so dispatch wiring is unchanged."""

    async def _on_mouse_down(self, event: MouseDown) -> None:
        if event.button != 1:
            return
        if self.disabled or not self.display:
            return
        event.stop()
        self.post_message(self.Pressed(self))


_ANSI16_HEX: list[str] = [
    "#000000", "#800000", "#008000", "#808000",
    "#000080", "#800080", "#008080", "#C0C0C0",
    "#808080", "#FF0000", "#00FF00", "#FFFF00",
    "#0000FF", "#FF00FF", "#00FFFF", "#FFFFFF",
]


@_functools.lru_cache(maxsize=256)
def _xterm_index_to_hex(idx: int) -> str:
    """Convert an xterm-256 color index (0..255) to the closest 24-bit RGB
    hex. Matches the xterm default palette — terminals may remap these but
    the vast majority follow the spec. Cube levels use the canonical
    ``[0, 95, 135, 175, 215, 255]`` ramp; grayscale uses
    ``8 + 10 * k`` for k in 0..23.

    LRU-cached at maxsize=256 (entire palette) — `_XtermColorGrid.render`
    calls this 256× per mount and the output is deterministic."""
    idx = max(0, min(255, int(idx)))
    if idx < 16:
        return _ANSI16_HEX[idx]
    if idx < 232:
        n = idx - 16
        levels = (0, 95, 135, 175, 215, 255)
        r = levels[(n // 36) % 6]
        g = levels[(n // 6)  % 6]
        b = levels[ n        % 6]
        return f"#{r:02X}{g:02X}{b:02X}"
    v = 8 + 10 * (idx - 232)
    return f"#{v:02X}{v:02X}{v:02X}"


class _XtermColorGrid(Static):
    """Single-Static replacement for the 256 individual color-cell
    Buttons that used to make up the xterm grid in ColorPickerModal.

    Profiling on a T480s baseline showed `ColorPickerModal.push_screen`
    settling in ~2 s, dominated by mounting 256 Buttons + iterating
    them in `on_mount` to set per-cell `styles.background`. Each
    Button is a full Textual widget with its own CSS pipeline; 256
    of them blow past the cost of any color picker should ever
    impose. This widget paints the entire grid as one Rich Text
    canvas (each cell = 3 spaces with a coloured background) and
    hit-tests clicks to a cell index via integer arithmetic — three
    orders of magnitude fewer widgets, ~30× faster modal mount.

    Layout (each cell 3 chars wide × 1 row tall, matches the
    legacy `.colorpick-xterm-cell` CSS rule):

      row 0       — 16 ANSI cells          (cells 0-15)
      rows 1..6   — 216-color cube         (cells 16-231, 36/row)
      row 7       — 24 grayscale cells     (cells 232-255)

    `cell_at(x, y)` is the inverse mapping used by
    `ColorPickerModal._cell_index_at` for click + drag.
    """

    DEFAULT_CSS = """
    _XtermColorGrid {
        height: 8;
        width: auto;
    }
    """

    _CELL_W = 3   # widget-relative pixel width per cell (matches old CSS)

    def render(self) -> Text:
        t = Text(no_wrap=True, overflow="crop")
        cell = " " * self._CELL_W
        # Row 0: 16 ANSI
        for i in range(16):
            t.append(cell, style=f"on {_xterm_index_to_hex(i)}")
        t.append("\n")
        # Rows 1..6: 216-cube
        for row in range(6):
            for col in range(36):
                idx = 16 + row * 36 + col
                t.append(cell, style=f"on {_xterm_index_to_hex(idx)}")
            t.append("\n")
        # Row 7: 24 grayscale
        for i in range(232, 256):
            t.append(cell, style=f"on {_xterm_index_to_hex(i)}")
        return t

    def cell_at(self, x: int, y: int) -> "int | None":
        """Convert widget-relative `(x, y)` to xterm cell index or
        `None` if the click landed outside any cell. `x` is in
        cells (Textual character columns); divide by `_CELL_W`."""
        col = x // self._CELL_W
        if col < 0:
            return None
        if y == 0:
            return col if 0 <= col < 16 else None
        if 1 <= y <= 6:
            return 16 + (y - 1) * 36 + col if 0 <= col < 36 else None
        if y == 7:
            return 232 + col if 0 <= col < 24 else None
        return None


# ── File-picker / search widget styling (moved from hub, Phase D) ───────────
_PICKER_PLASMID_STYLE = "bold #BFFF00"   # lime green — .gb / .gbk / .genbank
_PICKER_OTHER_STYLE   = "#FFFFFF"        # plain white for everything else
_FASTA_PICKER_FASTA_STYLE = "bold #FF69B4"
_FASTA_PICKER_OTHER_STYLE = "#FFFFFF"
_SEQ_ZIP_HIGHLIGHT_STYLE = "bold #BFFF00"
_SEQ_ZIP_OTHER_STYLE     = "#FFFFFF"


# ── Custom file-picker / search / table widgets (moved from hub, Phase D) ────
class _SearchInput(Input):
    """Reusable search-bar Input.

    Three reusable behaviours that every search bar in the app wants
    consistently:

    1. **Prefill / placeholder hybrid.** ``PREFILL`` is shown when the
       field is unfocused + empty; clicking into the field clears it
       so the user types into a fresh cursor. Blurring an empty field
       restores the prefill. The PREFILL string also counts as
       "empty" in `query()` so a user mashing Enter without focusing
       first doesn't search for the literal word "Search".

    2. **Optional live-filter debouncing.** When the parent passes
       ``debounce_s=<seconds>`` + ``on_filter=<callable>``, every
       keystroke schedules a single timer; rapid typing coalesces
       into one call. With ``debounce_s=None`` (default) the widget
       behaves identically to plain ``Input`` and the parent listens
       to ``Input.Submitted`` for Enter-to-apply UX. The debounce
       timer is cancelled on unmount so a queued tick can't fire
       against a disposed widget tree.

    3. **Length cap.** Defaults to ``DEFAULT_MAX_LEN = 200``. Longer
       pastes (e.g. an accidental clipboard dump) are silently
       truncated. The fuzzy matcher is O(query × name); without a
       cap a 100k-char paste against a 5k-plasmid library would lock
       the UI for seconds before any debounce fires.

    Pre-2026-05-10 this class had only the focus-clear behaviour and
    every site that wanted debounce / cleanup / cap rolled its own.
    The four search bars in the app (LibraryPanel, LibrarySearchModal,
    LoadPartSourceModal, SpeciesPickerModal) used 4 slightly-different
    patterns; the consolidated widget is the one place to fix
    behaviour bugs going forward.
    """
    PREFILL = "Search"
    DEFAULT_MAX_LEN = 200

    def __init__(
        self,
        *args,
        prefill: "str | None" = None,
        debounce_s: "float | None" = None,
        max_len: int = DEFAULT_MAX_LEN,
        on_filter=None,
        **kwargs,
    ):
        # Default prefill to the class constant; allow opt-out by
        # passing ``prefill=""`` so a caller that wants placeholder-
        # only behaviour gets a plain Input visually.
        if prefill is None:
            prefill = self.PREFILL
        # Use the prefill as the initial value if the caller didn't
        # pass one explicitly. ``None`` is the sentinel for "use the
        # prefill"; ``""`` is the sentinel for "start empty".
        if "value" not in kwargs:
            kwargs["value"] = prefill
        super().__init__(*args, **kwargs)
        self._prefill = prefill
        self._debounce_s = debounce_s
        self._max_len = max(1, int(max_len))
        self._on_filter = on_filter
        self._filter_timer = None
        # Latched by `on_unmount`. Cancelling the pending timer there is NOT
        # enough on its own: tearing the widget down changes `value` (focus
        # is dropped, the prefill is restored), which posts a fresh
        # `Input.Changed` — and if that lands AFTER `on_unmount`, it arms a
        # brand-new timer that nothing will ever cancel. Observed as a
        # load-dependent flake where the filter fired with an EMPTY query
        # (`['']`, not the typed `['queued']`) a fifth of a second after the
        # app had gone away. Once latched, this widget never schedules or
        # runs a filter again.
        self._filter_disabled = False

    def on_focus(self, _event) -> None:
        # Always blank the field on focus gain — matches the spec
        # "clicking into … the textbox clears and a cursor appears".
        if self.value == self._prefill:
            self.value = ""

    def on_blur(self, _event) -> None:
        # Restoring the prefill on blur (empty field only) keeps the
        # idle UI readable — without this, a user who clicked away
        # without typing would see an empty field with no affordance
        # cue. Whitespace-only counts as empty.
        if not self.value.strip():
            try:
                self.value = self._prefill
            except Exception:
                # Defensive — Input.value setter could fail mid-
                # teardown; the next on_focus / clear() will recover.
                pass

    def on_input_changed(self, event) -> None:
        # Length cap is a hard truncate. Setting `self.value` here
        # re-triggers Input.Changed; the second pass reads the
        # truncated value, hits the equality check, and skips the
        # truncate path. So we don't need an explicit re-entrancy
        # guard — but we DO need to bail before scheduling a debounce
        # for the truncate-pass, since we'll be back here in 1 µs.
        if len(self.value) > self._max_len:
            self.value = self.value[: self._max_len]
            return
        # No callback / no debounce → behave like a plain Input.
        if self._on_filter is None or self._debounce_s is None:
            return
        # Teardown re-entry: the value changes as the widget is disposed,
        # and arming a timer here outlives the tree it would fire against.
        if self._filter_disabled:
            return
        # Cancel any pending tick before scheduling a fresh one;
        # without this a burst of N keystrokes spawns N timers and
        # the filter callback fires N times after debounce_s.
        if self._filter_timer is not None:
            try:
                self._filter_timer.stop()
            except Exception:
                pass
        try:
            self._filter_timer = self.set_timer(
                self._debounce_s, self._fire_filter,
            )
        except Exception:
            # Timer setup can race with unmount; fall back to a
            # synchronous fire so the user's keystroke isn't lost.
            self._filter_timer = None
            self._fire_filter()

    def _fire_filter(self) -> None:
        """Invoke the ``on_filter`` callback with the sanitised query.
        Wrapped in try/except so a callback that raises (e.g. a stale
        widget query during unmount) doesn't poison the timer."""
        if self._on_filter is None:
            return
        # Backstop: a timer already DISPATCHED when unmount landed can no
        # longer be stopped, so the gate has to sit here too. Callbacks that
        # raise are merely logged; one that quietly MUTATES state would act
        # on a dead screen, which is the failure this exists to prevent.
        if self._filter_disabled:
            return
        try:
            self._on_filter(self.current_query())
        except Exception:
            _log.exception(
                "_SearchInput: on_filter callback raised on query %r",
                self.value,
            )

    def current_query(self) -> str:
        """Return the sanitised query string. Strips outer whitespace
        and treats the prefill as empty so a user who didn't focus
        the field doesn't search for the placeholder text.

        Named ``current_query`` (not ``query``) to avoid colliding
        with ``Widget.query(selector)`` from Textual's DOMNode API,
        which would otherwise be shadowed by the same-name override.
        """
        v = (self.value or "").strip()
        if v == self._prefill:
            return ""
        return v

    def clear(self) -> None:
        """Reset the input to prefill state and cancel any pending
        debounce. Idempotent — safe to call from a parent's
        view-mode-switch path that also fires `on_unmount` shortly
        afterward."""
        try:
            self.value = self._prefill
        except Exception:
            pass
        if self._filter_timer is not None:
            try:
                self._filter_timer.stop()
            except Exception:
                pass
            self._filter_timer = None

    def on_unmount(self) -> None:
        # Latch FIRST, then cancel. Setting the flag before stopping the
        # timer closes the window where a teardown-driven `Input.Changed`
        # re-arms after the cancel — the load-dependent path that let a
        # filter fire against a disposed tree.
        self._filter_disabled = True
        # Cancel the pending debounce so a queued tick doesn't fire
        # against a disposed widget tree — the callback typically
        # `query_one`s into the parent's table, which would raise
        # NoMatches (caught) but log a noisy warning. Wiping the
        # slot keeps the tear-down clean.
        if self._filter_timer is not None:
            try:
                self._filter_timer.stop()
            except Exception:
                pass
            self._filter_timer = None


class _SingleClickDataTable(DataTable):
    """Sweep #30 (2026-05-26): DataTable variant that fires
    `CellSelected` on the FIRST click of any data cell, instead
    of requiring two clicks (cursor-move + select). Textual's
    stock `DataTable._on_click` moves the cursor on first click
    and only posts `CellSelected` on the second click of an
    already-highlighted cell — that's confusing in the members
    table where the strand / color / label cells are meant to be
    direct one-click affordances. User-reported: "still no
    change to the segment list after trying to change arrow"
    because the picker never opened on single click.

    Header / row-label clicks fall through to the parent
    implementation so sort + label-click semantics are preserved.
    """

    async def _on_click(self, event) -> None:
        meta = event.style.meta
        if "row" not in meta or "column" not in meta:
            return
        if (self.cursor_type != "row"
                and meta.get("out_of_bounds", False)):
            return
        row_index = meta["row"]
        column_index = meta["column"]
        is_header_click = (self.show_header and row_index == -1)
        is_row_label_click = (self.show_row_labels
                              and column_index == -1)
        if is_header_click or is_row_label_click:
            # Defer to the stock handler for header / label
            # clicks so HeaderSelected / RowLabelSelected fire
            # exactly as Textual ships them.
            await super()._on_click(event)
            return
        if self.show_cursor and self.cursor_type != "none":
            new_coordinate = _Coordinate(row_index, column_index)
            self.cursor_coordinate = new_coordinate
            # The key change vs stock Textual: always post the
            # selected message, regardless of whether the cell
            # was already the cursor's location. Single-click
            # opens the picker.
            self._post_selected_message()
            self._scroll_cursor_into_view(animate=True)
            event.stop()


class _ExtensionAwareDirectoryTree(DirectoryTree):
    """DirectoryTree that highlights files per-extension via a colour map.

    Construction options (mutually compatible — `highlight_map` wins
    when both are supplied):

      `highlight_map: dict[str, str]` — explicit ``{".ext": "style"}``.
        Lets a single tree colour different formats with different
        styles (e.g. pink for FASTA, orange for .dna, green for .gb).
      `highlight_exts: frozenset[str]` + `highlight_style: str` — legacy
        single-colour API. Every extension in the set renders in
        `highlight_style`; everything else in `other_style`.

    Files that don't match any rule render in `other_style` (white by
    default). Directories keep Textual's default styling so folder
    navigation cues stay obvious.

    Replaces the per-modal `_FastaAwareDirectoryTree` /
    `_ZipAwareDirectoryTree` triplet — each modal just constructs this
    with the highlight rules it cares about.
    """

    def __init__(self, path, *,
                  highlight_map: "dict[str, str] | None" = None,
                  highlight_exts: "frozenset[str] | None" = None,
                  highlight_style: str = _PICKER_PLASMID_STYLE,
                  other_style:     str = _PICKER_OTHER_STYLE,
                  **kwargs) -> None:
        super().__init__(path, **kwargs)
        # Lower-case once on construction so the per-render check is a
        # single dict lookup; lets callers write the literal {".gb": ...}
        # without worrying about case variants on the input side.
        if highlight_map is not None:
            self._highlight_map: dict[str, str] = {
                k.lower(): v for k, v in highlight_map.items()
            }
        elif highlight_exts is not None:
            self._highlight_map = {
                e.lower(): highlight_style for e in highlight_exts
            }
        else:
            self._highlight_map = {}
        self._other_style = other_style

    def render_label(self, node, base_style, style):
        label = super().render_label(node, base_style, style)
        data = node.data
        if data is None:
            return label
        p = getattr(data, "path", None)
        if p is None:
            return label
        try:
            if not p.is_file():
                return label
        except OSError:
            return label
        styled = label.copy()
        suffix = (getattr(p, "suffix", "") or "").lower()
        styled.stylize(self._highlight_map.get(suffix, self._other_style))
        return styled

    def _populate_node(self, node, content):
        """Re-sort directory contents with `_natural_sort_key` before
        the base populates the tree, so ``FFE 2`` sorts before
        ``FFE 10`` (numerical-aware), not after (lexicographic).

        The base ``DirectoryTree._load_directory`` sorts by
        ``(not is_dir, name.lower())`` — directories first, files
        alphabetically by case-insensitive name. We keep the
        directories-first half but swap the file-name comparator
        to the same natural-sort key the plasmid library panel uses,
        so a folder full of `FFE 2 ENTRY A1.dna` / `FFE 10 …` sorts
        as a human would expect.
        """
        sorted_content = sorted(
            content,
            key=lambda p: (
                not self._safe_is_dir(p),
                _natural_sort_key(p.name),
            ),
        )
        super()._populate_node(node, sorted_content)


class _FastaAwareDirectoryTree(DirectoryTree):
    """DirectoryTree variant that colours FASTA files lime green and
    every other file white. Directories are left alone so Textual's
    default folder styling still applies. Sorts via the natural-key
    comparator so ``FFE 2`` precedes ``FFE 10`` — same convention
    the plasmid library panel uses."""

    def render_label(self, node, base_style, style):
        label = super().render_label(node, base_style, style)
        data = node.data
        if data is None:
            return label
        p = getattr(data, "path", None)
        if p is None:
            return label
        try:
            if not p.is_file():
                return label
        except OSError:
            return label
        styled = label.copy()
        if _is_fasta_path(p):
            styled.stylize(_FASTA_PICKER_FASTA_STYLE)
        else:
            styled.stylize(_FASTA_PICKER_OTHER_STYLE)
        return styled

    def _populate_node(self, node, content):
        sorted_content = sorted(
            content,
            key=lambda p: (
                not self._safe_is_dir(p),
                _natural_sort_key(p.name),
            ),
        )
        super()._populate_node(node, sorted_content)


class _ZipAwareDirectoryTree(DirectoryTree):
    """DirectoryTree variant that colours .zip archives lime green so
    the user can scan a downloads folder for the Plasmidsaurus run.
    Mirrors `_FastaAwareDirectoryTree`'s contract — including the
    natural-key sort so `run_2.zip` precedes `run_10.zip`."""

    def render_label(self, node, base_style, style):
        label = super().render_label(node, base_style, style)
        data = node.data
        if data is None:
            return label
        p = getattr(data, "path", None)
        if p is None:
            return label
        try:
            if not p.is_file():
                return label
        except OSError:
            return label
        styled = label.copy()
        if _is_seq_zip_path(p):
            styled.stylize(_SEQ_ZIP_HIGHLIGHT_STYLE)
        else:
            styled.stylize(_SEQ_ZIP_OTHER_STYLE)
        return styled

    def _populate_node(self, node, content):
        sorted_content = sorted(
            content,
            key=lambda p: (
                not self._safe_is_dir(p),
                _natural_sort_key(p.name),
            ),
        )
        super()._populate_node(node, sorted_content)


# ── Feature/Rich color helpers (moved from hub, Phase D) ────────────────────
_FEATURE_PALETTE: list[str] = [
    "color(39)",   "color(118)",  "color(208)",  "color(213)",  "color(51)",
    "color(220)",  "color(196)",  "color(46)",   "color(201)",  "color(129)",
    "color(166)",  "color(33)",   "color(226)",  "color(160)",  "color(87)",
    "color(105)",  "color(154)",  "color(203)",  "color(81)",   "color(185)",
]



_HEX3_RE = re.compile(r"^#[0-9A-Fa-f]{3}$")


_HEX6_RE = re.compile(r"^#[0-9A-Fa-f]{6}$")


_XTERM_RE = re.compile(r"^(?:color\()?(\d{1,3})\)?$")


def _resolve_feature_color(entry: dict) -> str:
    """Effective render color for a library entry. Never returns an empty
    string — falls back through entry override → user default → built-in
    default → palette[0]. Always produces something Rich will parse.

    Palette values that use ``color(N)`` syntax are safe inside Rich
    Style.parse but blow up Rich's markup lexer when rewrapped as
    ``[color(N)]...[/]``. Normalise them to their hex equivalent before
    returning so every downstream caller can use plain markup templating.
    """
    col = entry.get("color") if isinstance(entry, dict) else None
    if isinstance(col, str) and col:
        return _markup_safe_color(col)
    ft = entry.get("feature_type", "") if isinstance(entry, dict) else ""
    user_defaults = _load_feature_colors()
    if ft in user_defaults:
        return _markup_safe_color(user_defaults[ft])
    if ft in _DEFAULT_TYPE_COLORS:
        return _markup_safe_color(_DEFAULT_TYPE_COLORS[ft])
    return _markup_safe_color(_FEATURE_PALETTE[0])


def _markup_safe_color(col: str) -> str:
    """Return ``col`` unchanged if it's already safe for Rich markup (hex,
    named color, or ``rgb(r,g,b)``). Palette-style ``color(N)`` values
    are converted to hex via :func:`_normalise_color_input` because the
    parens collide with Rich's ``[...]`` lexer. Falls back to returning
    the raw string if normalisation fails — better to render something
    than to crash the preview.

    Sweep #34 (2026-05-26): also strip any Rich markup
    metacharacters (`[`, `]`) from the input. Pre-fix a
    `.gb` file with a hand-crafted `color="[red]X[/red]"`
    qualifier would round-trip through this helper unchanged
    and break out of `f"[bold {color}]..."` interpolation,
    corrupting the rendered label. No legitimate colour value
    contains brackets; rejecting them stops the markup
    injection cleanly.
    """
    if not isinstance(col, str) or not col:
        return col
    if "[" in col or "]" in col:
        # Drop the suspect value; the caller falls back to the
        # palette default via `_resolve_feature_color`'s chain.
        return _FEATURE_PALETTE[0]
    if col.startswith("color("):
        canon = _normalise_color_input(col)
        if canon:
            return canon
    return col


def _normalise_color_input(raw: str) -> "str | None":
    """Accept any of ``#RGB`` / ``#RRGGBB`` / ``0..255`` / ``color(N)`` and
    return a canonical hex string. Returns ``None`` if the input isn't
    parseable — callers surface a validation error.

    Hex is preferred over ``color(N)`` on save because hex survives round-
    trip through JSON unambiguously and renders identically on every
    terminal that supports truecolor. xterm indices are converted to their
    24-bit RGB equivalent via the standard 6×6×6 cube / grayscale ramp.
    """
    if not isinstance(raw, str):
        return None
    raw = raw.strip()
    if not raw:
        return None
    if _HEX6_RE.match(raw):
        return "#" + raw[1:].upper()
    if _HEX3_RE.match(raw):
        r, g, b = raw[1], raw[2], raw[3]
        return f"#{r}{r}{g}{g}{b}{b}".upper()
    m = _XTERM_RE.match(raw)
    if m:
        try:
            idx = int(m.group(1))
        except ValueError:
            return None
        if 0 <= idx <= 255:
            return _xterm_index_to_hex(idx)
    return None


# ── Plasmid-picker highlight / status / ext maps (moved from hub, Phase D) ──
_PLASMID_TEXT_EXTS: frozenset[str] = frozenset(
    {".gb", ".gbk", ".genbank"}
)


_PICKER_DNA_STYLE     = "bold #FFA500"   # orange    — .dna (binary plasmid)


_PICKER_FASTA_STYLE   = "bold #FF69B4"   # hot pink  — .fa / .fasta / .fna / …


_PLASMID_PICKER_HIGHLIGHT_MAP: dict[str, str] = {
    **{e: _PICKER_PLASMID_STYLE for e in _PLASMID_TEXT_EXTS},
    ".dna": _PICKER_DNA_STYLE,
    **{e: _PICKER_FASTA_STYLE for e in _FASTA_EXTS},
}


_PLASMID_STATUS_COLORS: dict[str, str] = {
    "DESIGNING":  "#B975FF",   # purple
    "CLONING":    "#FFA62B",   # orange (matches the app's accent)
    "SEQUENCING": "#5FB3FF",   # blue
    "VERIFIED":   "#4EBF71",   # green (matches success accent)
    "ERROR":      "#FF5C5C",   # red  (matches error/divergent accent)
}
