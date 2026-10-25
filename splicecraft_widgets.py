"""Self-contained Textual widget primitives (layer 3).

`_InstantPressButton` (fires Pressed on mouse-down, used by StrandPickerModal)
and the xterm-256 color grid (`_XtermColorGrid`, used by ColorPickerModal --
its `_xterm_index_to_hex` / `_ANSI16_HEX` palette helpers live in util L0 now
and are re-exported here). Pure Textual/Rich widgets -- no dependency on hub
state or persistence.
"""
from __future__ import annotations

import re

from rich.text import Text
from textual.coordinate import Coordinate as _Coordinate
from textual.events import MouseDown
from textual import work as _textual_work
from textual.widgets import Button, DataTable, DirectoryTree, Input, Static

from splicecraft_logging import _log
from splicecraft_util import (
    _ANSI16_HEX as _ANSI16_HEX,
    _DEFAULT_TYPE_COLORS, _FASTA_EXTS, _is_fasta_path, _is_seq_zip_path, _natural_sort_key,
    _xterm_index_to_hex as _xterm_index_to_hex,
)
from splicecraft_dataaccess import _load_feature_colors


def _safe_work(method=None, /, **kwargs):
    """Textual's ``@work`` with ``exit_on_error=False`` as the DEFAULT.

    Textual's own default ENDS THE APP when an exception escapes a worker
    (`app._handle_exception(WorkerFailed)`), taking unsaved edits with it;
    `PlasmidApp._on_worker_state_changed` — the log + toast safety net
    written for exactly that case — only ever ran after the shutdown had
    begun (audit 2026-10-23: an error raised in the autosave worker's
    deepcopy exited the app with return code 1). With this default the
    worker ends in its ERROR state, the net reports it, and the session
    lives on. A worker that genuinely should take the app down can still
    pass ``exit_on_error=True``. Hub and modals bind this as ``work``."""
    kwargs.setdefault("exit_on_error", False)
    if method is None:
        return _textual_work(**kwargs)
    return _textual_work(method, **kwargs)


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


# ── Virtualised row table ────────────────────────────────────────────────────

from rich.segment import Segment as _Segment          # noqa: E402
from textual import events as _events                 # noqa: E402
from textual.binding import Binding as _Binding       # noqa: E402
from textual.geometry import Size as _Size            # noqa: E402
from textual.message import Message as _Message       # noqa: E402
from textual.scroll_view import ScrollView as _ScrollView  # noqa: E402
from textual.strip import Strip as _Strip             # noqa: E402

# C0 / DEL / C1 control characters. A cell is written to the terminal
# verbatim, so an escape byte in a feature label (a file is untrusted input)
# would reach the terminal as a live control sequence.
_VRT_CONTROL_RE = re.compile(r"[\x00-\x1f\x7f-\x9f]")


class _VirtualColumn:
    """One column of a `_VirtualRowTable`: a header label and the width the
    cells are drawn at. `width` IS the rendered width (cells are truncated
    with an ellipsis or padded to it)."""

    __slots__ = ("key", "label", "width")

    def __init__(self, key: str, label: str, width: int) -> None:
        self.key = key
        self.label = label
        self.width = max(0, int(width))


class _VirtualRowTable(_ScrollView, can_focus=True):
    """A row-cursor table that only ever draws the rows on screen.

    Textual's `DataTable` materialises every row: `add_row` builds per-row
    bookkeeping, and the next idle pass renders AND measures every cell of
    every new row (`_update_dimensions`). That is roughly 200 µs a row, so a
    50,000-feature chromosome took seconds to list — and the Features sidebar
    re-lists on every map rotation, because rotating re-sorts it.

    This table keeps no per-row objects. The owner supplies a row count and a
    `row_cells(i)` callback returning the cells of display row `i`, and gives
    every column an explicit width; `render_line` draws only the visible
    lines. Replacing or re-ordering the rows is therefore O(1) here — the
    owner's own sort is the only O(n) work left.

    Messages carry the click's modifiers, because the owner needs them and a
    parent cannot see the `Click` itself (the table consumes it):

    * `RowHighlighted(table, cursor_row, shift=)` — the cursor moved to a
      different row (keyboard, click, or `move_cursor`).
    * `RowSelected(table, cursor_row, shift=, chain=)` — a click on the row
      the cursor is already on (`chain` >= 2 for a double-click), or Enter.

    Header geometry matches `DataTable`: the header is screen line 0 and data
    row `i` is drawn at line `i + 1 - scroll_y`.
    """

    DEFAULT_CSS = """
    _VirtualRowTable {
        background: $surface;
        color: $foreground;
        height: auto;
        max-height: 100%;
    }
    _VirtualRowTable:focus {
        background-tint: $foreground 5%;
    }
    _VirtualRowTable > .vrt--header {
        text-style: bold;
        background: $panel;
        color: $foreground;
    }
    _VirtualRowTable > .vrt--even-row {
        background: $surface-lighten-1 50%;
    }
    _VirtualRowTable:dark > .vrt--even-row {
        background: $surface-darken-1 40%;
    }
    _VirtualRowTable > .vrt--odd-row {
    }
    _VirtualRowTable > .vrt--cursor {
        background: $block-cursor-blurred-background;
        color: $block-cursor-blurred-foreground;
        text-style: $block-cursor-blurred-text-style;
    }
    _VirtualRowTable:focus > .vrt--cursor {
        background: $block-cursor-background;
        color: $block-cursor-foreground;
        text-style: $block-cursor-text-style;
    }
    """

    COMPONENT_CLASSES = {
        "vrt--header", "vrt--cursor", "vrt--even-row", "vrt--odd-row",
    }

    BINDINGS = [
        _Binding("enter", "select_cursor", "Select", show=False),
        _Binding("up", "cursor_up", "Cursor up", show=False),
        _Binding("down", "cursor_down", "Cursor down", show=False),
        _Binding("pageup", "page_up", "Page up", show=False),
        _Binding("pagedown", "page_down", "Page down", show=False),
        _Binding("ctrl+home", "cursor_top", "Top", show=False),
        _Binding("ctrl+end", "cursor_bottom", "Bottom", show=False),
        _Binding("home", "cursor_top", "Top", show=False),
        _Binding("end", "cursor_bottom", "Bottom", show=False),
        _Binding("left", "scroll_left", "Scroll left", show=False),
        _Binding("right", "scroll_right", "Scroll right", show=False),
    ]

    CELL_PADDING = 1

    class RowHighlighted(_Message):
        def __init__(self, table: "_VirtualRowTable", cursor_row: int,
                     *, shift: bool = False) -> None:
            self.table = table
            self.cursor_row = cursor_row
            self.shift = bool(shift)
            super().__init__()

        @property
        def control(self) -> "_VirtualRowTable":
            return self.table

    class RowSelected(_Message):
        def __init__(self, table: "_VirtualRowTable", cursor_row: int,
                     *, shift: bool = False, chain: int = 1) -> None:
            self.table = table
            self.cursor_row = cursor_row
            self.shift = bool(shift)
            self.chain = int(chain)
            super().__init__()

        @property
        def control(self) -> "_VirtualRowTable":
            return self.table

    def __init__(self, *, zebra_stripes: bool = True, show_header: bool = True,
                 name: "str | None" = None, id: "str | None" = None,
                 classes: "str | None" = None,
                 disabled: bool = False) -> None:
        super().__init__(name=name, id=id, classes=classes, disabled=disabled)
        self.zebra_stripes = bool(zebra_stripes)
        self.show_header = bool(show_header)
        self.columns: "dict[str, _VirtualColumn]" = {}
        self._row_count = 0
        self._row_cells = None
        self._cursor_row = 0

    # ── data ────────────────────────────────────────────────────────────────

    def add_columns(self, *labels: str) -> "list[str]":
        """Add columns sized to their header label; set real widths with
        `set_column_widths`."""
        keys: "list[str]" = []
        for label in labels:
            key = f"col{len(self.columns)}"
            self.columns[key] = _VirtualColumn(key, str(label), len(str(label)))
            keys.append(key)
        self._update_virtual_size()
        return keys

    def set_column_widths(self, widths) -> None:
        for col, w in zip(self.columns.values(), widths):
            col.width = max(0, int(w))
        self._update_virtual_size()
        self.refresh()

    def set_rows(self, row_count: int, row_cells) -> None:
        """Replace every row. `row_cells(i)` returns the cells of display row
        `i` (a `Text` or anything `str()`-able per column); it is only called
        for rows being drawn. The cursor returns to row 0 without posting
        `RowHighlighted` — like `DataTable.clear()` + re-add, a repopulate is
        not a user move."""
        self._row_count = max(0, int(row_count))
        self._row_cells = row_cells
        self._cursor_row = 0
        self._update_virtual_size()
        self.scroll_to(x=0, y=0, animate=False, force=True)
        self.refresh()

    def clear(self) -> "_VirtualRowTable":
        self.set_rows(0, None)
        return self

    @property
    def row_count(self) -> int:
        return self._row_count

    @property
    def cursor_row(self) -> int:
        return self._cursor_row

    def _total_width(self) -> int:
        pad = 2 * self.CELL_PADDING
        return sum(c.width + pad for c in self.columns.values())

    def _header_height(self) -> int:
        return 1 if self.show_header else 0

    def _update_virtual_size(self) -> None:
        self.virtual_size = _Size(self._total_width(),
                                  self._row_count + self._header_height())

    def _body_height(self) -> int:
        return max(1, self.scrollable_content_region.height
                   - self._header_height())

    # ── cursor ──────────────────────────────────────────────────────────────

    def move_cursor(self, *, row: "int | None" = None, animate: bool = False,
                    scroll: bool = True) -> None:
        """Programmatic cursor move. Posts `RowHighlighted` when the row
        changes, exactly like `DataTable.move_cursor` (owners that drive it
        from elsewhere suppress the echo themselves)."""
        if row is None:
            return
        self._set_cursor(row, scroll=scroll)

    def _set_cursor(self, row: int, *, shift: bool = False,
                    scroll: bool = True) -> None:
        if self._row_count <= 0:
            self._cursor_row = 0
            return
        row = max(0, min(self._row_count - 1, int(row)))
        old = self._cursor_row
        self._cursor_row = row
        if scroll:
            self._scroll_cursor_into_view()
        if row != old:
            self.refresh()
            self.post_message(self.RowHighlighted(self, row, shift=shift))

    def _scroll_cursor_into_view(self) -> None:
        top = int(self.scroll_offset.y)
        body = self._body_height()
        row = self._cursor_row
        if row < top:
            self.scroll_to(y=row, animate=False, force=True)
        elif row >= top + body:
            self.scroll_to(y=row - body + 1, animate=False, force=True)

    def action_cursor_up(self) -> None:
        self._set_cursor(self._cursor_row - 1)

    def action_cursor_down(self) -> None:
        self._set_cursor(self._cursor_row + 1)

    def action_page_up(self) -> None:
        self._set_cursor(self._cursor_row - self._body_height())

    def action_page_down(self) -> None:
        self._set_cursor(self._cursor_row + self._body_height())

    def action_cursor_top(self) -> None:
        self._set_cursor(0)

    def action_cursor_bottom(self) -> None:
        self._set_cursor(self._row_count - 1)

    def action_select_cursor(self) -> None:
        if 0 <= self._cursor_row < self._row_count:
            self.post_message(self.RowSelected(self, self._cursor_row))

    async def _on_click(self, event: _events.Click) -> None:
        offset = event.get_content_offset(self)
        if offset is None:
            return
        header = self._header_height()
        if offset.y < header:
            return
        row = int(self.scroll_offset.y) + offset.y - header
        if not 0 <= row < self._row_count:
            return
        event.stop()
        # Ctrl is a Shift synonym: many terminals keep Shift+click for their
        # own text selection.
        shift = bool(getattr(event, "shift", False)
                     or getattr(event, "ctrl", False))
        chain = int(getattr(event, "chain", 1) or 1)
        if row == self._cursor_row:
            self.post_message(self.RowSelected(self, row, shift=shift,
                                               chain=chain))
        else:
            self._set_cursor(row, shift=shift)

    # ── rendering ───────────────────────────────────────────────────────────

    @staticmethod
    def _cell_text(cell) -> Text:
        if isinstance(cell, Text):
            text = cell.copy()
            if not text.plain.isprintable():
                text = Text(_VRT_CONTROL_RE.sub("?", text.plain),
                            style=text.style)
            return text
        s = "" if cell is None else str(cell)
        if not s.isprintable():
            s = _VRT_CONTROL_RE.sub("?", s)
        return Text(s)

    def _line_text(self, cells, header: bool = False) -> Text:
        pad = " " * self.CELL_PADDING
        line = Text(no_wrap=True, overflow="crop", end="")
        cols = list(self.columns.values())
        if header:
            cells = [c.label for c in cols]
        for col, cell in zip(cols, cells):
            t = self._cell_text(cell)
            t.truncate(col.width, overflow="ellipsis", pad=True)
            line.append(pad)
            line.append_text(t)
            line.append(pad)
        return line

    def render_line(self, y: int) -> _Strip:
        width = self.size.width
        scroll_x, scroll_y = self.scroll_offset
        base = self.rich_style
        header = self._header_height()
        if header and y == 0:
            style = base + self.get_component_rich_style("vrt--header")
            text = self._line_text([], header=True)
            post = None
        else:
            row = scroll_y + y - header
            if not (0 <= row < self._row_count) or self._row_cells is None:
                return _Strip.blank(width, base)
            try:
                cells = self._row_cells(row)
            except Exception:
                _log.exception("_VirtualRowTable: row_cells(%d) failed", row)
                cells = ()
            text = self._line_text(cells)
            if row == self._cursor_row:
                style = base
                post = self.get_component_rich_style("vrt--cursor")
            else:
                stripe = ("vrt--even-row" if row % 2 == 0 else "vrt--odd-row")
                style = (base + self.get_component_rich_style(stripe)
                         if self.zebra_stripes else base)
                post = None
        segments = list(_Segment.apply_style(
            text.render(self.app.console), style, post))
        fill = style + post if post is not None else style
        strip = _Strip(segments).crop_extend(scroll_x, scroll_x + width, fill)
        return strip.apply_offsets(scroll_x, y)
