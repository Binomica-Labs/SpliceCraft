"""splicecraft_modals — extracted ModalScreen / Screen dialog classes (Phase D).

The hub's modal/screen dialog classes, relocated here so context-limited models
(and humans) can load the UI dialog layer as a bounded unit instead of scrolling
the 100k-line hub. The hub re-exports every class
(`from splicecraft_modals import X as X`) so `sc.<Class>` and every
`push_screen(...)` / `isinstance(...)` site resolves unchanged.

Layer L4 (see tests/test_import_layers.py::_LAYER_RULES): imports textual + the
lower siblings (logging L0 / widgets L3 / dataaccess L1 / render L1) ONLY. Modals
reach the running app through Textual's `self.app` (a runtime attribute) and pop
results via `self.dismiss(...)`; any `PlasmidApp` parameter/return type hints are
quoted strings that are never evaluated (`from __future__ import annotations`), so
this module does NOT import the hub — that would be an import cycle.

Started with the dependency-closure-clean leaf dialogs (confirm / prompt / picker
modals that reference only textual + `self.app` + each other). Modals that still
call bare hub helpers stay hub-side until those helpers move to siblings.
"""
from __future__ import annotations

import re

from datetime import date as _date
from pathlib import Path
from rich.text import Text
from textual import on
from textual.app import ComposeResult
from textual.await_complete import AwaitComplete
from textual.binding import Binding
from textual.containers import Horizontal, ScrollableContainer, Vertical, VerticalScroll
from textual.coordinate import Coordinate as _Coordinate
from textual.css.query import NoMatches
from textual.events import Click, MouseDown, MouseMove, MouseUp
from textual.screen import ModalScreen, Screen
from textual.widgets import Button, DataTable, DirectoryTree, Input, Label, ListItem, ListView, RadioButton, RadioSet, Select, Static, TextArea, Tree

import splicecraft_state as _state
from splicecraft_cloning import _simulate_cloned_plasmid, _simulate_primed_amplicon
from splicecraft_dataaccess import _BUILTIN_GRAMMARS, _all_grammars, _collection_name_taken, _find_collection, _find_hmm_db_entry, _get_active_collection_name, _grammar_dropdown_options, _hmm_db_name_taken, _iter_collections_readonly, _iter_library_readonly, _load_collections, _load_feature_colors, _load_library, _load_primer_collections, _normalise_hmm_db_entry, _sanitize_hmm_db_id, _sanitize_hmm_db_url, _save_collections, _search_collections_library
from splicecraft_history import _CommercialSaaSHistoryNode, _history_consistency_summary, _history_detail_lines, _history_former_name, _history_node_warnings, _history_populate_tree, _history_protocol_renderable, _history_tree_label
from splicecraft_logging import _log, _log_event
from splicecraft_util import _CONTROL_CHARS_RE, _PLASMID_STATUS_VALUES, _cursor_row_key, _natural_sort_key, _normalize_collection_name, _notify_save_failure, _primer_tm_safe, _sanitize_label, _sanitize_plasmid_name, _sanitize_plasmid_status, _scrub_path, _validate_group_members
from splicecraft_widgets import _DEFAULT_TYPE_COLORS, _ExtensionAwareDirectoryTree, _FastaAwareDirectoryTree, _HEX6_RE, _InstantPressButton, _PICKER_PLASMID_STYLE, _PLASMID_STATUS_COLORS, _SearchInput, _XtermColorGrid, _ZipAwareDirectoryTree, _markup_safe_color, _normalise_color_input, _xterm_index_to_hex



# ── Unsaved-changes quit dialog ────────────────────────────────────────────────

class UnsavedQuitModal(ModalScreen):
    """Shown when the user tries to quit with unsaved edits."""

    _blocks_undo: bool = True   # [INV-50] destructive confirm — Ctrl+Z above

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next button", show=False),
    ]

    def __init__(self) -> None:
        super().__init__()
        self._dismissed: bool = False

    def _dismiss_once(self, payload) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="quit-dlg"):
            yield Static(" Unsaved Changes ", id="quit-title")
            yield Static(
                "  You have unsaved edits. What would you like to do?",
                id="quit-msg",
            )
            with Horizontal(id="quit-btns"):
                yield Button("Save & Quit",      id="btn-save-quit", variant="primary")
                yield Button("Abandon Changes",  id="btn-abandon",   variant="error")
                yield Button("Cancel",           id="btn-cancel-quit")

    def on_mount(self) -> None:
        # Default focus on Cancel — match the "default No / safe" pattern
        # of every other confirm modal so a hammered Enter can't quit.
        self.query_one("#btn-cancel-quit", Button).focus()

    @on(Button.Pressed, "#btn-save-quit")
    def _save_quit(self, _): self._dismiss_once("save")

    @on(Button.Pressed, "#btn-abandon")
    def _abandon(self, _):   self._dismiss_once("abandon")

    @on(Button.Pressed, "#btn-cancel-quit")
    def _cancel_btn(self, _): self._dismiss_once(None)

    def action_cancel(self): self._dismiss_once(None)


class QuitConfirmModal(ModalScreen):
    """Confirm-quit modal for the no-unsaved-changes case. The unsaved
    branch goes through `UnsavedQuitModal` (with Save / Abandon / Cancel)
    instead. Default focus on `No`.
    """

    _blocks_undo: bool = True   # [INV-50] confirm modal — Ctrl+Z above

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self) -> None:
        super().__init__()
        self._dismissed: bool = False

    def _dismiss_once(self, payload: bool) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="quitcon-dlg"):
            yield Static(" Quit SpliceCraft? ", id="quitcon-title")
            yield Static(
                "  Are you sure you want to quit?",
                id="quitcon-msg",
            )
            with Horizontal(id="quitcon-btns"):
                yield Button("No",  id="btn-quitcon-no",  variant="default")
                yield Button("Yes", id="btn-quitcon-yes", variant="error")

    def on_mount(self) -> None:
        self.query_one("#btn-quitcon-no", Button).focus()

    @on(Button.Pressed, "#btn-quitcon-no")
    def _no(self, _): self._dismiss_once(False)

    @on(Button.Pressed, "#btn-quitcon-yes")
    def _yes(self, _): self._dismiss_once(True)

    def action_cancel(self): self._dismiss_once(False)


class _OneShotDismissScreen(Screen):
    """Mixin that makes ``dismiss()`` idempotent.

    Real terminals can deliver TWO ``Button.Pressed`` events for one
    physical click (the focus-transition + click-cycle interaction
    documented across sweeps #36 / #38 / #39). Textual's
    ``Screen.dismiss()`` (8.x) fires the result callback AND calls
    ``app.pop_screen()`` with no internal guard, so a double-fire would
    (a) re-run the caller's callback — e.g. a duplicate library save —
    and (b) pop a SECOND screen off the stack (the parent) or raise
    ``ScreenStackError``. Mixing this in BEFORE ``ModalScreen`` /
    ``Screen`` makes the first dismiss flip a one-shot flag and delegate
    to ``super().dismiss()``; every later call is an inert no-op.

    Self-contained (a class-attribute flag, no ``__init__``) so any
    modal can opt in with a one-line base change and nothing else. The
    return type matches ``Screen.dismiss`` (``AwaitComplete``) so the
    two-base MRO stays type-compatible; the no-op path returns ``None``
    at runtime (never awaited from a handler) to avoid building an
    ``AwaitComplete`` outside a running loop.
    """

    _dismissed_once: bool = False

    def dismiss(self, result=None) -> AwaitComplete:
        if self._dismissed_once:
            return None  # type: ignore[return-value]
        self._dismissed_once = True
        return super().dismiss(result)


class EditSeqDialog(ModalScreen):
    """Unified sequence-edit modal — insert (left/right), replace, or
    delete.

    Sweep #39 (2026-05-27): refactored from a single-mode dialog into
    a context-aware operations modal. User picks the mode in-modal
    via a RadioSet; the caller seeds the initial mode based on
    whether a selection or cursor is present.

    Modes (internal ``_mode`` names) — user request 2026-06-01 split
    Insert into LEFT / RIGHT and prefilled the textbox with the region.
    Ctrl+E (`action_edit_seq`) opens with a NON-destructive default —
    Replace for a selection, Insert-left for a bare cursor (never Delete,
    a hand-slip risk); the Delete KEY with no feature selected opens this
    same dialog forced to Delete mode (`_open_seq_edit_dialog`):
      * ``"insert_left"`` — insert ``new_seq`` at ``start`` (before the
        region / cursor bp). Dismisses as ``(new_seq, "insert", start,
        start)``.
      * ``"insert_right"`` — insert ``new_seq`` at ``end`` (after the
        region / cursor bp). Dismisses as ``(new_seq, "insert", end,
        end)``. A bare cursor opens as a 1-bp region ``[pos, pos+1)`` so
        left vs right is meaningful.
      * ``"replace"`` — replace ``seq[start:end]`` with ``new_seq``.
      * ``"delete"`` — drop ``seq[start:end]``. Returned as
        ``("", "replace", start, end)`` so the existing dispatch
        in ``_edit_dialog_result`` handles it as a replace-with-
        empty. ``"insert"`` is kept as a back-compat alias → insert_left.
        (Copy was dropped 2026-06-01 — five radios overflowed the dialog
        and the seq panel's Ctrl+C already copies a selection.)

    A colour-coded live WARNING line ("You are about to <action> <N> bp
    <position>") tracks the focused mode + textbox as the user types.

    Input handling:
      * ``TextArea`` instead of ``Input`` — multi-line, so paste
        of a multi-line FASTA / long fragment / spread-out paste
        from a browser works without each newline tripping
        validation.
      * Sanitization on submit (``_sanitize_pasted_sequence``):
        drop FASTA header lines (``>``), strip all whitespace,
        upper-case, U→T (RNA paste compatibility).
      * Validates the cleaned string against IUPAC alphabet.

    Returns via dismiss:
      * insert_left:  ``(new_seq, "insert", start, start)``
      * insert_right: ``(new_seq, "insert", end, end)``
      * replace:      ``(new_seq, "replace", start, end)``
      * delete:       ``("", "replace", start, end)``
      * cancel:       ``None``
    """

    _blocks_undo: bool = True   # sweep #10: app-level Ctrl+Z would race the insert

    _VALID = frozenset("ATCGNRYSWKMBDHV")   # IUPAC DNA codes
    # Insert is split into LEFT (insert before the region / cursor bp)
    # and RIGHT (insert after it) — user request 2026-06-01. The bare
    # ``"insert"`` name is kept as a back-compat alias that ``__init__``
    # maps to insert-left (old callers + boundary tests pass "insert").
    _MODE_INSERT_LEFT  = "insert_left"
    _MODE_INSERT_RIGHT = "insert_right"
    _MODE_INSERT       = "insert"          # legacy alias → insert_left
    _MODE_REPLACE = "replace"
    _MODE_DELETE  = "delete"
    # Copy was dropped 2026-06-01 — five horizontal radios overflowed the
    # dialog (user: "ui crowded and lil messy"), and copy is redundant
    # with the seq panel's Ctrl+C. The four real edit ops fit cleanly.
    _MODES        = (_MODE_INSERT_LEFT, _MODE_INSERT_RIGHT,
                     _MODE_REPLACE, _MODE_DELETE)
    _INSERT_MODES = (_MODE_INSERT_LEFT, _MODE_INSERT_RIGHT)
    _MODES_NEEDING_INPUT = (_MODE_INSERT_LEFT, _MODE_INSERT_RIGHT,
                            _MODE_REPLACE)
    _MODES_NEEDING_REGION = (_MODE_REPLACE, _MODE_DELETE)
    # Sweep #39 hardening (2026-05-27): cap pasted input at 200 kbp
    # (matches `_MAX_FEATURE_SEQ_LEN` precedent — anything longer
    # belongs as a record in the plasmid library, not as a one-shot
    # paste into the seq editor). The cap stops a paste-attack /
    # accidental megabyte-paste from freezing the TextArea + the
    # sanitiser. Counted after sanitisation so FASTA-decorated
    # pastes (headers + whitespace) get their real bp count
    # measured.
    _MAX_INPUT_BP = 200_000
    # Cap on how many bp we prefill into the textbox (2026-06-01). A
    # huge selection (e.g. the whole plasmid → Ctrl+E) would be slow to
    # render in the TextArea AND — for Replace, whose submit reads the
    # textbox — a truncated prefill would silently commit truncated
    # bases. Above the cap the box opens EMPTY: the context line + the
    # warning still convey the region size, and Delete / Copy act on the
    # region itself (never the textbox), so nothing is lost.
    _PREFILL_MAX_BP = 20_000

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("ctrl+s", "submit", "Submit"),
        Binding("ctrl+enter", "submit", "Submit"),
    ]

    DEFAULT_CSS = """
    #edit-dlg {
        width: 84; max-width: 96%; min-width: 76; height: auto;
        background: $surface; border: solid $primary;
        padding: 1 2;
    }
    #edit-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; height: 1; text-align: center;
    }
    #edit-context {
        height: auto; color: $text-muted;
        margin: 1 0 0 0; padding: 0 1;
    }
    #edit-mode-row {
        height: 3; margin: 1 0 0 0;
        align: left middle;
    }
    #edit-mode-row RadioSet {
        layout: horizontal; height: 3; width: 1fr;
        border: none; padding: 0;
    }
    #edit-mode-row RadioButton {
        width: auto; padding: 0 1; margin-right: 2;
    }
    #edit-preview {
        height: 1; margin: 0 0 1 0;
        color: $accent; padding: 0 1;
    }
    #edit-label {
        margin: 0 0 0 0; padding: 0 1;
    }
    #edit-input-area {
        height: 6; margin: 0 0 1 0;
        border: solid $primary-darken-2;
    }
    #edit-err { height: 1; padding: 0 1; }
    #edit-btns { height: 3; margin-top: 1; align: right middle; }
    #edit-btns Button { margin-right: 1; }
    """

    def __init__(self, mode: str, existing: str = "",
                 start: int = 0, end: int = 0,
                 total: "int | None" = None):
        super().__init__()
        # Legacy "insert" alias → insert-left.
        if mode == self._MODE_INSERT:
            mode = self._MODE_INSERT_LEFT
        if mode not in self._MODES:
            mode = self._MODE_INSERT_LEFT
        # If the caller picked a region-required mode but didn't
        # actually supply a region, demote to insert-left so the modal
        # doesn't open in a state that can't be confirmed.
        if mode in self._MODES_NEEDING_REGION and end <= start:
            mode = self._MODE_INSERT_LEFT
        self._mode       = mode
        self._existing   = existing
        self._start      = start
        self._end        = end
        self._total      = (int(total) if total is not None
                            else (start + len(existing)))
        self._has_region = end > start
        # Sweep #39 hardening: one-shot dismiss guard — a double-
        # click on OK / Esc / button-mouse-up-after-mouse-down posts
        # two Pressed events in real terminals (see sweep #36 +
        # InstantPressButton). Without the guard the second dispatch
        # tries to dismiss an already-popped screen and could re-fire
        # the caller's _edit_dialog_result twice (which would commit
        # the same edit twice + leave two undo snapshots).
        self._dismissed: bool = False

    # ── Layout ─────────────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        with Vertical(id="edit-dlg"):
            yield Static(" Edit sequence ", id="edit-title")
            yield Static(self._build_context_line(), id="edit-context",
                          markup=True)
            with Horizontal(id="edit-mode-row"):
                with RadioSet(id="edit-mode"):
                    yield RadioButton(
                        "Insert left", id="edit-mode-insert_left",
                        value=(self._mode == self._MODE_INSERT_LEFT),
                    )
                    yield RadioButton(
                        "Insert right", id="edit-mode-insert_right",
                        value=(self._mode == self._MODE_INSERT_RIGHT),
                    )
                    yield RadioButton(
                        "Replace", id="edit-mode-replace",
                        value=(self._mode == self._MODE_REPLACE),
                    )
                    yield RadioButton(
                        "Delete", id="edit-mode-delete",
                        value=(self._mode == self._MODE_DELETE),
                    )
            yield Static(self._build_preview_line(), id="edit-preview",
                          markup=True)
            yield Label(
                "New sequence  (IUPAC: A T C G N R Y W S K M B D H V — "
                "paste OK; FASTA headers + whitespace + U auto-stripped):",
                id="edit-label",
            )
            yield TextArea(id="edit-input-area")
            yield Static("", id="edit-err", markup=True)
            with Horizontal(id="edit-btns"):
                # Sweep #39 hardening: InstantPressButton so a
                # single physical click registers on first try in
                # real terminals (same focus-transition workaround
                # sweep #36 applied to ColorPickerModal). Paired
                # with the dismiss-once guard so the click-cycle
                # double-fire doesn't double-dispatch the result.
                yield _InstantPressButton(
                    "OK  [Ctrl+S]", id="btn-ok",
                    variant="primary",
                )
                yield _InstantPressButton(
                    "Cancel  [Esc]", id="btn-cancel",
                )

    def _prefill_text(self) -> str:
        """Bases to seed the textbox with — the region's existing bases,
        but only when small enough that rendering is snappy AND (for
        Replace, whose submit reads the box) it round-trips losslessly.
        A larger region yields ``""`` (see ``_PREFILL_MAX_BP``)."""
        return (self._existing
                if len(self._existing) <= self._PREFILL_MAX_BP else "")

    def on_mount(self) -> None:
        # Prefill the textbox with the region's existing bases so the
        # user sees exactly what's selected / under the cursor. For
        # delete / copy it's a read-only preview of what's affected; for
        # insert / replace it's an editable starting point. Capped — a
        # very large selection opens with an empty box (see
        # `_prefill_text`).
        try:
            # Insert modes open with an EMPTY box (you're typing NEW
            # bases — a prefilled region would insert a duplicate);
            # replace / delete / copy prefill with the region.
            self.query_one("#edit-input-area", TextArea).text = (
                "" if self._mode in self._INSERT_MODES
                else self._prefill_text())
        except NoMatches:
            pass
        self._apply_mode_visibility()
        # Focus the input for Insert/Replace; focus OK for the no-input
        # modes (Delete is the default on open — a single Enter confirms
        # the prefilled deletion).
        try:
            if self._mode in self._MODES_NEEDING_INPUT:
                self.query_one("#edit-input-area", TextArea).focus()
            else:
                self.query_one("#btn-ok", Button).focus()
        except NoMatches:
            pass
        # Sync the RadioSet's keyboard-highlight cursor onto the CHECKED
        # button. Textual's `RadioSet._on_mount` runs `action_next_button`
        # which highlights index 0 regardless of which radio is
        # `value=True` — so a non-first default (Replace, or Delete via
        # the Delete key) leaves the "-selected" highlight box on Insert
        # left while a different radio is filled (user report 2026-06-01:
        # "the delete radio is selected but the first radio is
        # highlighted"). Deferred so it runs AFTER RadioSet._on_mount.
        def _sync_radio_highlight() -> None:
            try:
                rs = self.query_one("#edit-mode", RadioSet)
                idx = rs.pressed_index
                if idx >= 0:
                    rs._selected = idx
            except Exception:
                # `_selected` is private — if a future Textual drops it,
                # the highlight just stays put (purely cosmetic).
                _log.debug(
                    "EditSeqDialog: radio highlight sync skipped",
                    exc_info=True,
                )
        self.call_after_refresh(_sync_radio_highlight)

    # ── Mode + preview state ───────────────────────────────────────────────────

    def _build_context_line(self) -> str:
        """One-line summary of what the modal is operating on."""
        if self._has_region:
            length = self._end - self._start
            excerpt = self._existing[:32]
            if len(self._existing) > 32:
                excerpt += "…"
            excerpt_safe = excerpt.replace("[", r"\[")
            return (
                f"Selection: [b]{length:,} bp[/]  "
                f"({self._start + 1:,}‥{self._end:,} of "
                f"{self._total:,})   "
                f"[dim]{excerpt_safe}[/]"
            )
        return (
            f"Cursor at position [b]{self._start + 1:,}[/] of "
            f"{self._total:,} bp  [dim](no selection — pick "
            f"bases first to enable Replace / Delete / Copy)[/]"
        )

    def _build_preview_line(self) -> str:
        """The live WARNING line: "You are about to <action> <N> bp
        <position>", with the action verb, the bp count, and the
        position phrase colour-coded by action (green insert / red
        delete / yellow replace / cyan copy). Rebuilt whenever the mode
        changes or the textarea contents change (user spec 2026-06-01).

        ``#bp`` is the bases the action WRITES for the input modes
        (insert / replace — so it updates live as the user types) and
        the region size for delete / copy.
        """
        n = len(self._sanitize_pasted_sequence(self._current_input_text()))
        region = self._end - self._start
        if self._mode == self._MODE_INSERT_LEFT:
            verb, count, where, color = (
                "insert", n, "to the left of cursor", "green")
        elif self._mode == self._MODE_INSERT_RIGHT:
            verb, count, where, color = (
                "insert", n, "to the right of cursor", "green")
        elif self._mode == self._MODE_REPLACE:
            if not self._has_region:
                return ("[yellow]Replace needs a selection — switch to "
                        "an Insert mode or pick bases first.[/]")
            verb, count, where, color = "replace", n, "at cursor", "yellow"
        elif self._mode == self._MODE_DELETE:
            if not self._has_region:
                return ("[yellow]Delete needs a selection — pick "
                        "bases first.[/]")
            verb, count, where, color = "delete", region, "at cursor", "red"
        else:
            return ""
        # Colour-coded action verb + bp count + position phrase.
        return (
            f"You are about to [b {color}]{verb}[/] "
            f"[b {color}]{count:,} bp[/] [{color}]{where}[/]"
        )

    def _apply_mode_visibility(self) -> None:
        """Keep the textbox visible in EVERY mode (user request
        2026-06-01 — the modal always shows the bases being acted on).
        For insert / replace it's the editable input; for delete / copy
        it's a READ-ONLY preview of the region the action affects, so the
        label adapts to never read "New sequence" over a preview."""
        needs_input = self._mode in self._MODES_NEEDING_INPUT
        try:
            ta = self.query_one("#edit-input-area", TextArea)
            ta.read_only = not needs_input
        except NoMatches:
            pass
        try:
            lbl = self.query_one("#edit-label", Label)
            note = ("too large to preview — see the count above"
                    if len(self._existing) > self._PREFILL_MAX_BP
                    else "read-only preview")
            if self._mode == self._MODE_DELETE:
                lbl.update(f"Bases that will be deleted ({note}):")
            else:
                lbl.update(
                    "New sequence  (IUPAC: A T C G N R Y W S K M B D H V — "
                    "paste OK; FASTA headers + whitespace + U auto-stripped):"
                )
        except NoMatches:
            pass
        self._refresh_preview()

    def _current_input_text(self) -> str:
        try:
            return self.query_one("#edit-input-area", TextArea).text
        except NoMatches:
            return ""

    def _refresh_preview(self) -> None:
        try:
            self.query_one(
                "#edit-preview", Static,
            ).update(self._build_preview_line())
        except NoMatches:
            pass

    # ── Sanitization ───────────────────────────────────────────────────────────

    # Sweep #39 hardening (2026-05-27): pre-compile the sanitiser
    # regexes so a fast-typist's keystroke doesn't trigger a
    # re-compile per change. The leading position-number stripper
    # is line-anchored (only eats digits at the START of a line)
    # so a future IUPAC extension that allows in-sequence digits
    # can't be silently corrupted — and a paste like "1234ATCG"
    # (no whitespace) won't lose the leading "1234". NCBI's
    # numbered FASTA always puts the position number BEFORE
    # whitespace, so anchoring to "^\d+\s+" matches the real
    # use case without over-stripping.
    _SANITIZE_WS_RE = re.compile(r"\s+")
    _SANITIZE_LEADING_POS_RE = re.compile(r"^\d+\s+")

    def _sanitize_pasted_sequence(self, raw: str) -> str:
        """Clean up pasted text:
          * Drop FASTA header lines (lines starting with ``>``).
          * Drop NCBI-style line-leading position numbers
            (``"  60 ATCG"`` → ``"ATCG"``). Anchored to start-of-
            line so in-sequence digits are left alone for
            validation to catch.
          * Strip whitespace + newlines.
          * Upper-case, U→T (RNA→DNA).
        Returns the cleaned sequence — no validation.
        """
        if not raw:
            return ""
        out_parts: list[str] = []
        for line in raw.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith(">"):
                continue
            # Strip NCBI-style leading position number ("60 ATCG…").
            stripped = self._SANITIZE_LEADING_POS_RE.sub("", stripped)
            out_parts.append(stripped)
        joined = "".join(out_parts)
        # Strip any remaining whitespace inside lines (spaces
        # between blocks in a 60-char-wide FASTA, tabs, etc.).
        joined = self._SANITIZE_WS_RE.sub("", joined)
        return joined.upper().replace("U", "T")

    def dismiss(self, result=None):  # type: ignore[override]
        """Sweep #39 hardening: one-shot dismiss. The Pressed
        message pump can deliver two Pressed events for a single
        physical button click in real terminals (sweep #36
        documents the focus-transition + click-cycle interaction);
        without this guard the second event would try to dismiss
        an already-popped screen and could fire the caller's
        ``_edit_dialog_result`` twice — committing the same edit
        twice and stacking two undo snapshots.

        Liskov note: base ``Screen.dismiss`` returns ``AwaitComplete``;
        we widen to ``AwaitComplete | None`` because the guard MUST
        NOT call ``super().dismiss`` again on the second event (the
        whole point of the one-shot). ``# type: ignore[override]``
        documents this is intentional — no caller awaits this
        modal's dismiss return value."""
        if self._dismissed:
            return None
        self._dismissed = True
        return super().dismiss(result)

    # ── Event wiring ──────────────────────────────────────────────────────────

    @on(RadioSet.Changed, "#edit-mode")
    def _on_mode_changed(self, event: RadioSet.Changed) -> None:
        btn = event.pressed
        if btn is None:
            return
        chosen = (btn.id or "").rsplit("-", 1)[-1]
        if chosen not in self._MODES or chosen == self._mode:
            return
        prev = self._mode
        self._mode = chosen
        # Textbox content follows the mode:
        #   * delete / copy → show the region as a read-only preview;
        #   * insert        → clear an UNTOUCHED prefill (or a preview
        #                     carried over from delete/copy) so the box
        #                     starts blank for fresh bases, but keep a
        #                     real edit the user already typed;
        #   * replace       → seed with the region if the box is empty.
        try:
            ta = self.query_one("#edit-input-area", TextArea)
            pf = self._prefill_text()
            if chosen == self._MODE_DELETE:
                ta.text = pf
            elif chosen in self._INSERT_MODES:
                if ta.text == pf or prev == self._MODE_DELETE:
                    ta.text = ""
            elif chosen == self._MODE_REPLACE:
                # Seed with the region when arriving from an insert with
                # an empty box, OR from a delete read-only preview (whose
                # text — possibly an empty huge-region placeholder — must
                # not become the replacement payload).
                if not ta.text or prev == self._MODE_DELETE:
                    ta.text = pf
        except NoMatches:
            pass
        self._apply_mode_visibility()
        # Refocus the relevant widget for fast keyboard-driven flow.
        try:
            if chosen in self._MODES_NEEDING_INPUT:
                self.query_one("#edit-input-area", TextArea).focus()
            else:
                self.query_one("#btn-ok", Button).focus()
        except NoMatches:
            pass

    @on(TextArea.Changed, "#edit-input-area")
    def _on_input_changed(self, _event: TextArea.Changed) -> None:
        raw = self._current_input_text()
        cleaned = self._sanitize_pasted_sequence(raw)
        try:
            err = self.query_one("#edit-err", Static)
        except NoMatches:
            return
        if not cleaned:
            err.update("")
        elif len(cleaned) > self._MAX_INPUT_BP:
            # Sweep #39 hardening: surface the cap LIVE so the user
            # sees the refusal as they paste, before clicking OK.
            err.update(Text(
                f"{len(cleaned):,} bp — exceeds the "
                f"{self._MAX_INPUT_BP:,} bp cap for a single edit. "
                f"Trim the paste or save as a library entry.",
                style="bold red",
            ))
        else:
            bad = sorted(
                set(c for c in cleaned if c not in self._VALID),
            )
            if bad:
                err.update(Text(
                    f"Invalid: {'  '.join(repr(c) for c in bad)} — "
                    f"only IUPAC DNA codes allowed",
                    style="bold red",
                ))
            else:
                err.update(Text(
                    f"{len(cleaned):,} bp (after stripping "
                    f"whitespace / FASTA headers)",
                    style="dim green",
                ))
        self._refresh_preview()

    # ── Submit / cancel ───────────────────────────────────────────────────────

    def action_submit(self) -> None:
        self._try_submit()

    def action_cancel(self) -> None:
        self.dismiss(None)

    @on(Button.Pressed, "#btn-ok")
    def _ok(self, _) -> None:
        self._try_submit()

    @on(Button.Pressed, "#btn-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def _set_err(self, msg: str, *, severity: str = "red") -> None:
        try:
            self.query_one("#edit-err", Static).update(
                Text(msg, style=f"bold {severity}"),
            )
        except NoMatches:
            pass

    def _try_submit(self) -> None:
        if self._mode in (self._MODE_REPLACE, self._MODE_DELETE) \
                and not self._has_region:
            self._set_err(
                "This mode needs a selection — "
                "switch to an Insert mode or pick bases on the canvas first.",
            )
            return
        if self._mode == self._MODE_DELETE:
            # Empty bases through the replace dispatch.
            self.dismiss(("", self._MODE_REPLACE,
                           self._start, self._end))
            return
        # Insert / Replace: sanitize + validate input.
        cleaned = self._sanitize_pasted_sequence(
            self._current_input_text(),
        )
        if not cleaned:
            self._set_err(
                "Please enter or paste a sequence.",
            )
            return
        # Sweep #39 hardening: cap pasted size BEFORE
        # validation so a 50 MB paste that's all valid IUPAC
        # still gets refused instead of locking up the
        # canvas-rebuild + the restriction-scan worker on
        # commit. 200 kbp matches `_MAX_FEATURE_SEQ_LEN` and is
        # well above realistic plasmid-edit sizes.
        if len(cleaned) > self._MAX_INPUT_BP:
            self._set_err(
                f"Input is {len(cleaned):,} bp — exceeds the "
                f"{self._MAX_INPUT_BP:,} bp cap for a single "
                f"sequence edit. Save the fragment as a feature "
                f"or plasmid library entry instead.",
            )
            return
        bad = [c for c in cleaned if c not in self._VALID]
        if bad:
            return   # live validation already shows the error
        # Map the internal mode → the (new_bases, mode, s, e) payload
        # `_edit_dialog_result` expects ("insert" at a position, or
        # "replace" of [s, e)). Insert-left inserts at the region start,
        # insert-right at the region end (= after the selection / cursor
        # bp); both reuse the App's existing "insert" dispatch.
        if self._mode == self._MODE_INSERT_LEFT:
            self.dismiss((cleaned, "insert", self._start, self._start))
        elif self._mode == self._MODE_INSERT_RIGHT:
            self.dismiss((cleaned, "insert", self._end, self._end))
        else:  # replace
            self.dismiss((cleaned, "replace", self._start, self._end))


class MigrateDataModal(_OneShotDismissScreen, ModalScreen):
    """Migrate Data entry point — export ALL user data to one portable
    file, or import such a file into this install. Dismisses with
    ``"export"`` / ``"import"`` / ``None``."""

    _blocks_undo: bool = True
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    DEFAULT_CSS = """
    #migrate-box { width: 86; max-width: 95%; height: auto; max-height: 90%;
        background: $surface; border: solid $primary; padding: 1 2; }
    #migrate-title { background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1; text-align: center; }
    #migrate-blurb { color: $text-muted; margin-bottom: 1; }
    #migrate-btns { align: right middle;  height: 3; margin-top: 1; }
    #migrate-btns Button { margin-right: 1; }
    """

    def compose(self) -> ComposeResult:
        with Vertical(id="migrate-box"):
            yield Static(" Migrate Data ", id="migrate-title")
            yield Static(
                "Move ALL your SpliceCraft data between machines as one "
                "portable, compressed file.\n\n"
                "• Export — package your entire library, collections, parts "
                "bin, primers, features, grammars, codon tables, entry "
                "vectors, settings, lab notebook, and full construction "
                "history into a single .zip.\n"
                "• Import — load such a .zip into this install. Your current "
                "data is backed up automatically first, then replaced — so a "
                "fresh install picks up exactly where you left off.",
                id="migrate-blurb", markup=False,
            )
            with Horizontal(id="migrate-btns"):
                yield Button("Export all my data", id="btn-migrate-export",
                             variant="primary")
                yield Button("Import a data file", id="btn-migrate-import",
                             variant="primary")
                yield Button("Cancel", id="btn-migrate-cancel")

    @on(Button.Pressed, "#btn-migrate-export")
    def _exp(self, _) -> None:
        self.dismiss("export")

    @on(Button.Pressed, "#btn-migrate-import")
    def _imp(self, _) -> None:
        self.dismiss("import")

    @on(Button.Pressed, "#btn-migrate-cancel")
    def _cxl(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class DropdownScreen(_OneShotDismissScreen, ModalScreen):
    """Lightweight overlay showing a positioned dropdown menu.

    Uses a near-transparent backdrop so the main app stays visible — the
    dropdown looks like a real popup anchored to the menu bar, not a
    separate "screen". Click outside the box dismisses.
    """

    # Sweep #26: opening a menu shouldn't let an absent-minded Ctrl+Z
    # while the dropdown floats roll back the canvas underneath.
    _blocks_undo: bool = True

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    DEFAULT_CSS = """
    DropdownScreen {
        background: rgba(0, 0, 0, 0.15);
    }
    """

    def __init__(self, items: list, x: int, y: int) -> None:
        super().__init__()
        self._items = items   # (label, action_str | None)
        self._x = x
        self._y = y
        self._cursor = next(
            (i for i, (_, a) in enumerate(items) if a is not None), 0
        )

    def compose(self) -> ComposeResult:
        yield Static(
            self._build_dropdown_text(),
            id="dropdown-box",
        )

    def on_mount(self) -> None:
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        box_h   = len(self._items) + 2
        box = self.query_one("#dropdown-box", Static)
        box.styles.offset = (self._x, self._y)
        box.styles.width  = inner_w
        box.styles.height = box_h
        box.styles.border = ("solid", "#555555")
        box.styles.background = "#1e1e1e"

    def _build_dropdown_text(self) -> Text:
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        sep_line = "\u2500" * (inner_w - 2)
        result = Text()
        for i, (label, action) in enumerate(self._items):
            is_sep      = (label == "---")
            is_selected = (i == self._cursor and not is_sep and action is not None)
            is_disabled = (action is None and not is_sep)

            if is_sep:
                line = Text(sep_line + "\n", style="white")
            else:
                padded = f" {label:<{inner_w - 3}}"
                if is_selected:
                    line = Text(padded + "\n", style="reverse white")
                elif is_disabled:
                    line = Text(padded + "\n", style="dim white")
                else:
                    line = Text(padded + "\n", style="white")
            result.append_text(line)
        return result

    def _refresh_box(self) -> None:
        box = self.query_one("#dropdown-box", Static)
        box.update(self._build_dropdown_text())

    def on_key(self, event) -> None:
        items = self._items
        if event.key == "up":
            pos = self._cursor - 1
            while pos >= 0 and (items[pos][0] == "---" or items[pos][1] is None):
                pos -= 1
            if pos >= 0:
                self._cursor = pos
                self._refresh_box()
            event.stop()
        elif event.key == "down":
            pos = self._cursor + 1
            while pos < len(items) and (items[pos][0] == "---" or items[pos][1] is None):
                pos += 1
            if pos < len(items):
                self._cursor = pos
                self._refresh_box()
            event.stop()
        elif event.key == "enter":
            label, action = items[self._cursor]
            if action is not None:
                self.dismiss(action)
            event.stop()

    def on_click(self, event: Click) -> None:
        bx = self._x
        by = self._y
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        bh = len(self._items) + 2
        cx, cy = event.screen_x, event.screen_y
        if bx <= cx < bx + inner_w and by <= cy < by + bh:
            row_in_box = cy - by - 1  # -1 for top border
            if 0 <= row_in_box < len(self._items):
                label, action = self._items[row_in_box]
                if label == "---" or action is None:
                    event.stop()
                    return
                self.dismiss(action)
        else:
            self.dismiss(None)
        event.stop()

    def action_cancel(self) -> None:
        self.dismiss(None)


class PlasmidFeaturePickerModal(_OneShotDismissScreen, ModalScreen):
    """Scrollable list of non-source features from a specific library entry.

    Dismisses with a feature-library-style entry dict
    ``{name, feature_type, sequence, strand, qualifiers, description}``,
    or None on cancel. No persistence side effects — the caller decides
    whether to save the picked entry or just use it to prefill a form.
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, entries: list[dict], plasmid_name: str = ""):
        super().__init__()
        self._entries = list(entries)
        self._plasmid_name = plasmid_name or "plasmid"

    def compose(self) -> ComposeResult:
        with Vertical(id="featpick-dlg"):
            yield Static(f" Feature from [{self._plasmid_name}] ",
                         id="featpick-title")
            yield DataTable(id="featpick-table", cursor_type="row",
                            zebra_stripes=True)
            with Horizontal(id="featpick-btns"):
                yield Button("Select", id="btn-featpick-ok", variant="primary")
                yield Button("Cancel", id="btn-featpick-cancel")

    def on_mount(self) -> None:
        t = self.query_one("#featpick-table", DataTable)
        t.add_columns("Name", "Type", "Strand", "Length")
        for i, e in enumerate(self._entries):
            strand_str = "+" if e.get("strand", 1) == 1 else "−"
            t.add_row(
                e.get("name", "?"),
                e.get("feature_type", "?"),
                strand_str,
                f"{len(e.get('sequence', ''))} bp",
                key=str(i),
            )
        if self._entries:
            t.move_cursor(row=0)
            t.focus()

    @on(Button.Pressed, "#btn-featpick-ok")
    def _select(self, _):
        self._dismiss_cursor()

    @on(DataTable.RowSelected, "#featpick-table")
    def _row_selected(self, event):
        if event.row_key and event.row_key.value is not None:
            try:
                idx = int(event.row_key.value)
            except (TypeError, ValueError):
                return
            if 0 <= idx < len(self._entries):
                self.dismiss(dict(self._entries[idx]))

    def _dismiss_cursor(self) -> None:
        t = self.query_one("#featpick-table", DataTable)
        if t.row_count == 0:
            self.dismiss(None)
            return
        row_keys = list(t.rows.keys())
        if 0 <= t.cursor_row < len(row_keys):
            key = row_keys[t.cursor_row].value
            if key is None:
                self.dismiss(None)
                return
            try:
                idx = int(key)
            except (TypeError, ValueError):
                self.dismiss(None)
                return
            if 0 <= idx < len(self._entries):
                self.dismiss(dict(self._entries[idx]))
                return
        self.dismiss(None)

    @on(Button.Pressed, "#btn-featpick-cancel")
    def _cancel_btn(self, _):
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class UngroupScopeModal(ModalScreen):
    """Sweep #29: tiny confirm modal that asks whether `Ungroup`
    should drop the qualifier from JUST this feature or from
    every member of the group.

    Returns `{"scope": "this" | "whole"}` on Save, `None` on
    Cancel."""

    _blocks_undo: bool = True

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    DEFAULT_CSS = """
    #ungroup-dlg {
        width: 70; height: auto;
        background: $surface; border: solid $primary;
        padding: 1 2;
    }
    #ungroup-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1;
        text-align: center; text-style: bold;
    }
    #ungroup-msg { color: $text-muted; height: auto;
                   margin-bottom: 1; padding: 0 1; }
    #ungroup-btns { height: 3; margin-top: 1;
                    align: right middle; }
    #ungroup-btns Button { margin-right: 1; min-width: 18; }
    """

    def __init__(self, group_id: str) -> None:
        super().__init__()
        self._group_id = str(group_id or "")
        self._dismissed: bool = False

    def _dismiss_once(self, result) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(result)

    def compose(self) -> ComposeResult:
        with Vertical(id="ungroup-dlg"):
            yield Static(" Ungroup feature ", id="ungroup-title")
            yield Static(
                f"This feature is part of group "
                f"[b]…{self._group_id[-6:]}[/b]. Choose the scope "
                f"of the ungroup operation:\n\n"
                f"  • [b]Just this[/b] — drop the qualifier from "
                f"this feature only. Other group members keep "
                f"their shared id.\n"
                f"  • [b]Whole group[/b] — dissolve the group: "
                f"every member loses the qualifier.",
                id="ungroup-msg", markup=True,
            )
            with Horizontal(id="ungroup-btns"):
                yield Button("Just this", id="btn-ungroup-this",
                             variant="primary")
                yield Button("Whole group",
                             id="btn-ungroup-whole",
                             variant="warning")
                yield Button("Cancel", id="btn-ungroup-cancel")

    @on(Button.Pressed, "#btn-ungroup-this")
    def _this(self, _) -> None:
        self._dismiss_once({"scope": "this"})

    @on(Button.Pressed, "#btn-ungroup-whole")
    def _whole(self, _) -> None:
        self._dismiss_once({"scope": "whole"})

    @on(Button.Pressed, "#btn-ungroup-cancel")
    def _cancel(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)


class SplitPositionPromptModal(ModalScreen):
    """Sweep #30: tiny modal for "split a member at a position".

    Used by the unified AddFeatureModal / FeatureEditModal — the
    user selects a row in the members table and clicks "Split".
    This pops a single-number input pre-filled with the row's
    midpoint; on Save it returns `{"pos": <int>}`.

    Validates that the user-typed position falls strictly inside
    `(rs, re)` — equal-to either endpoint would produce a zero-
    width head or tail, which `_split_member` rejects. Status
    line surfaces the validation error rather than dismissing on
    a bad value so the user can correct it inline.
    """

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("enter",  "submit", "Split"),
    ]

    DEFAULT_CSS = """
    #splitpos-dlg {
        width: 60; height: auto;
        background: $surface; border: solid $primary;
        padding: 1 2;
    }
    #splitpos-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1;
        text-align: center; text-style: bold;
    }
    #splitpos-dlg Label { color: $text-muted; margin-top: 1; }
    #splitpos-input  { margin-top: 1; }
    #splitpos-status { height: 2; margin-top: 1; }
    #splitpos-btns   { height: 3; margin-top: 1;
                       align: right middle; }
    #splitpos-btns Button { margin-right: 1; min-width: 12; }
    """

    def __init__(self, rs: int, re_: int) -> None:
        super().__init__()
        if not isinstance(rs, int) or not isinstance(re_, int):
            raise ValueError("rs and re_ must be ints")
        if rs + 1 >= re_:
            # No interior position exists — a row of width 1 can't
            # be split. Caller (the modal's Split button handler)
            # is expected to surface this BEFORE opening this
            # prompt, but we guard defensively here too.
            raise ValueError(
                f"row [{rs}, {re_}) has no interior position to "
                f"split at (need width ≥ 2)"
            )
        self._rs = rs
        self._re = re_
        self._dismissed = False

    def _dismiss_once(self, result) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(result)

    def compose(self) -> ComposeResult:
        with Vertical(id="splitpos-dlg"):
            yield Static(" Split row at position ",
                          id="splitpos-title")
            yield Label(
                f"Position must satisfy "
                f"[b]{self._rs} < pos < {self._re}[/b]. "
                f"Default = midpoint.",
                markup=True,
            )
            default = (self._rs + self._re) // 2
            yield Input(
                value=str(default),
                placeholder=f"e.g. {default}",
                id="splitpos-input",
            )
            yield Static("", id="splitpos-status", markup=True)
            with Horizontal(id="splitpos-btns"):
                yield Button("Split", id="btn-splitpos-save",
                             variant="primary")
                yield Button("Cancel", id="btn-splitpos-cancel")

    def on_mount(self) -> None:
        try:
            inp = self.query_one("#splitpos-input", Input)
            inp.focus()
            # Select all so the user can immediately type to
            # overwrite the default.
            try:
                inp.action_select_all()
            except (AttributeError, NoMatches):
                pass
        except NoMatches:
            pass

    def action_submit(self) -> None:
        self._save(None)

    @on(Button.Pressed, "#btn-splitpos-save")
    def _save(self, _) -> None:
        try:
            raw = self.query_one("#splitpos-input", Input).value
        except NoMatches:
            return
        raw = (raw or "").strip()
        if not raw:
            self._set_status(
                "[red]Enter a position between "
                f"{self._rs} and {self._re}.[/red]"
            )
            return
        try:
            pos = int(raw)
        except ValueError:
            self._set_status("[red]Position must be an integer.[/red]")
            return
        if not (self._rs < pos < self._re):
            self._set_status(
                f"[red]Position must satisfy "
                f"{self._rs} < pos < {self._re} "
                f"(got {pos}).[/red]"
            )
            return
        self._dismiss_once({"pos": pos})

    @on(Button.Pressed, "#btn-splitpos-cancel")
    def _cancel(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)

    def _set_status(self, msg: str) -> None:
        try:
            self.query_one("#splitpos-status", Static).update(msg)
        except NoMatches:
            pass


class StrandPickerModal(ModalScreen):
    """Sweep #30 (2026-05-26): per-row strand picker for the
    members tables. Click the strand cell of a row in the
    AddFeatureModal / FeatureEditModal members table → this modal
    opens with the row's current strand pre-selected. User picks
    one of the four valid strands; modal dismisses with
    `{"strand": <int>}` (Forward=1 / Reverse=-1 / Arrowless=0 /
    Double=2) or `None` on Cancel.

    Hardening: the only way OUT with a non-Cancel result is via
    the four hard-coded button handlers, so the returned int is
    always in `{-1, 0, 1, 2}`. No free-text input → no
    sanitization needed. Modal is `_blocks_undo=True` so a
    pending strand change can't be wiped by Ctrl+Z while the
    picker is open."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    DEFAULT_CSS = """
    /* 2026-05-26: widened from 60 → 56 vertical stack so the
       4 strand buttons no longer overflow the dialog horizontally.
       The previous horizontal Horizontal layout pushed the last
       button (`↔ Both strands`) past the dialog border, where
       the clipped click hit-zone swallowed the first press and
       the user had to click a second time. Vertical stack keeps
       each button full-width inside the dialog so a single click
       lands cleanly. */
    #strand-pick-dlg {
        width: 56; height: auto; max-height: 90%;
        background: $surface; border: solid $primary;
        padding: 1 2;
    }
    #strand-pick-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1;
        text-align: center; text-style: bold;
    }
    #strand-pick-dlg Label { color: $text-muted; margin-top: 1; }
    #strand-pick-btns {
        height: auto; margin-top: 1;
        align-horizontal: center;
    }
    #strand-pick-btns Button {
        width: 100%;
        margin: 0 0 1 0;
    }
    #strand-pick-cancel-row {
        height: 3; margin-top: 1;
        align: right middle;
    }
    """

    def __init__(self, current_strand: int = 0,
                 *,
                 row_idx: "int | None" = None) -> None:
        super().__init__()
        # Sanitise the current strand on the way in: any value
        # outside {-1, 0, 1, 2} falls back to arrowless (0). This
        # defends the modal against a malformed caller without
        # crashing.
        try:
            current_strand = int(current_strand)
        except (TypeError, ValueError):
            current_strand = 0
        if current_strand not in (-1, 0, 1, 2):
            current_strand = 0
        self._current  = current_strand
        self._row_idx  = row_idx
        self._dismissed: bool = False

    def _dismiss_once(self, result) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(result)

    def compose(self) -> ComposeResult:
        title_suffix = (f" — row {self._row_idx + 1}"
                        if isinstance(self._row_idx, int) else "")
        cur = self._current
        with Vertical(id="strand-pick-dlg"):
            yield Static(
                f" Pick strand{title_suffix} ",
                id="strand-pick-title",
            )
            yield Label(
                "Choose the strand this segment lives on. "
                "Forward = top (5'→3' on the sense strand); "
                "Reverse = bottom (5'→3' on the antisense / "
                "complementary strand).",
            )
            # 2026-05-26: stack vertically so each button has the
            # full dialog width — the horizontal row clipped the
            # last two buttons and made single-click unreliable.
            # Buttons use `_InstantPressButton` so a single
            # mouse-down fires `Pressed` immediately, bypassing
            # the Textual focus-transition gate that swallowed
            # the first click in a real terminal.
            with Vertical(id="strand-pick-btns"):
                yield _InstantPressButton(
                    "▶ Top (forward)",  id="btn-strand-fwd",
                    variant=("primary" if cur == 1 else "default"),
                )
                yield _InstantPressButton(
                    "◀ Bottom (reverse)", id="btn-strand-rev",
                    variant=("primary" if cur == -1 else "default"),
                )
                yield _InstantPressButton(
                    "▒ None (arrowless)", id="btn-strand-none",
                    variant=("primary" if cur == 0 else "default"),
                )
                yield _InstantPressButton(
                    "↔ Both strands", id="btn-strand-both",
                    variant=("primary" if cur == 2 else "default"),
                )
            with Horizontal(id="strand-pick-cancel-row"):
                yield _InstantPressButton(
                    "Cancel", id="btn-strand-cancel",
                )

    def on_mount(self) -> None:
        """Focus the currently-selected strand button so the user
        can hit Enter / Space to confirm without a mouse round-
        trip. Falls back to forward when no match (defensive — the
        sanitiser in __init__ already clamps invalid input)."""
        sel = {
             1: "#btn-strand-fwd",
            -1: "#btn-strand-rev",
             0: "#btn-strand-none",
             2: "#btn-strand-both",
        }.get(self._current, "#btn-strand-fwd")
        try:
            self.query_one(sel, Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-strand-fwd")
    def _fwd(self, _) -> None:
        self._dismiss_once({"strand": 1})

    @on(Button.Pressed, "#btn-strand-rev")
    def _rev(self, _) -> None:
        self._dismiss_once({"strand": -1})

    @on(Button.Pressed, "#btn-strand-none")
    def _none(self, _) -> None:
        self._dismiss_once({"strand": 0})

    @on(Button.Pressed, "#btn-strand-both")
    def _both(self, _) -> None:
        self._dismiss_once({"strand": 2})

    @on(Button.Pressed, "#btn-strand-cancel")
    def _cancel(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)


class PrimerDuplicatesModal(_OneShotDismissScreen, ModalScreen):
    """Surfaced at startup when `primers.json` carries duplicate-
    sequence entries OR name-collision groups (same name, different
    sequences — common in `.dna` imports that replay primer entries
    across multiple plasmids with different binding-region / tail
    variants).

    Defaults to KEEP — user must explicitly click Delete to remove
    duplicates. The Keep button has initial focus so a stray Enter on
    splash dismiss doesn't silently delete primer data.

    Two cleanup passes on confirm:
      * Sequence collisions: first occurrence wins (sacred policy).
      * Name collisions: longest sequence wins (typically the variant
        carrying the full cloning tail, not a truncated binding only).

    Dismiss payload:
      True  — user confirmed deletion
      False / None — user kept everything (cancel / Escape / Keep)
    """

    _blocks_undo: bool = True   # primer-library write, no Ctrl+Z race

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #pdup-dlg {
        width: 70;
        height: auto;
        background: $surface;
        border: solid $warning;
        padding: 1 2;
    }
    #pdup-title {
        text-align: center;
        background: $warning-darken-2;
        color: $text;
        padding: 0 1;
    }
    #pdup-body {
        margin: 1 0;
        padding: 0 1;
    }
    #pdup-btns {
        height: 3;
        align: right middle;
        padding-top: 1;
    }
    #pdup-btns Button { margin-left: 1; }
    """

    def __init__(self, total: int, seq_duplicates: int,
                 name_collisions: int, final_kept: int):
        super().__init__()
        self._total = total
        self._seq_duplicates = seq_duplicates
        self._name_collisions = name_collisions
        self._final_kept = final_kept

    def compose(self) -> ComposeResult:
        to_remove = self._seq_duplicates + self._name_collisions
        lines: list[str] = []
        if self._seq_duplicates:
            # `_seq_duplicates` is the number of entries that WOULD
            # BE REMOVED, not the number that share — for a group of
            # 3 primers with the same sequence we'd remove 2 and keep
            # 1. Pre-fix this said "{2} entries share a DNA sequence"
            # which under-counts the actual sharers.
            lines.append(
                f"• {self._seq_duplicates} entr"
                f"{'ies' if self._seq_duplicates != 1 else 'y'} "
                f"would be removed for sharing a DNA sequence "
                f"(longest preserved per group)"
            )
        if self._name_collisions:
            lines.append(
                f"• {self._name_collisions} entr"
                f"{'ies' if self._name_collisions != 1 else 'y'} "
                f"would be removed for sharing a name "
                f"(longest sequence preserved per name)"
            )
        breakdown = "\n".join(lines)
        with Vertical(id="pdup-dlg"):
            yield Static(" Duplicate primers detected ", id="pdup-title")
            yield Static(
                f"Your primer library has {self._total} entries.\n"
                f"{breakdown}\n\n"
                f"Delete duplicates? Sequence collisions keep the "
                f"first occurrence; name collisions keep the entry "
                f"with the longest sequence (typically the full "
                f"primer including its cloning tail). Library would "
                f"end up with {self._final_kept} entries.",
                id="pdup-body",
            )
            with Horizontal(id="pdup-btns"):
                yield Button("Keep all", id="btn-pdup-keep",
                              variant="default")
                yield Button(
                    f"Delete {to_remove} entr"
                    f"{'ies' if to_remove != 1 else 'y'}",
                    id="btn-pdup-delete", variant="warning",
                )

    def on_mount(self) -> None:
        # Default focus on Keep so a stray Enter on splash dismiss
        # doesn't trigger a delete.
        try:
            self.query_one("#btn-pdup-keep", Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-pdup-keep")
    def _keep(self, _) -> None:
        self.dismiss(False)

    @on(Button.Pressed, "#btn-pdup-delete")
    def _delete(self, _) -> None:
        self.dismiss(True)

    def action_cancel(self) -> None:
        self.dismiss(False)


class EditGrammarConfirmModal(ModalScreen):
    """Warning before opening `GrammarEditorModal` on a custom
    grammar. Surfaces the fact that edits propagate to every
    TU/MOD already saved under this grammar — overhang positions
    can shift, enzyme assignments change which fragments digest
    cleanly, and the persisted parts may stop chaining at the
    next cycle. Default focus on `Cancel`."""

    # Sweep #26: [INV-50] destructive confirm — Ctrl+Z above
    # the modal could race the user's grammar-edit launch.
    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, grammar_name: str, n_dependents: int = 0) -> None:
        super().__init__()
        self._grammar_name = grammar_name
        self._n_dependents = n_dependents
        self._dismissed: bool = False   # [INV-50]

    def _dismiss_once(self, payload) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="gec-dlg"):
            yield Static(" Edit Grammar? ", id="gec-title")
            from rich.markup import escape as _esc
            nm = _esc(self._grammar_name)
            if self._n_dependents > 0:
                msg = (
                    f"  You're about to edit grammar [bold]{nm}[/bold].\n"
                    f"  [yellow]{self._n_dependents}[/yellow] saved "
                    f"part{'s' if self._n_dependents != 1 else ''} "
                    f"reference{'' if self._n_dependents != 1 else 's'} "
                    f"this grammar.\n"
                    f"  Changes to enzyme, recognition, or overhang "
                    f"positions can break those parts' chaining."
                )
            else:
                msg = (
                    f"  You're about to edit grammar [bold]{nm}[/bold].\n"
                    f"  Changes to enzyme, recognition site, or overhang "
                    f"positions can\n"
                    f"  invalidate any TUs / MODs you save under this "
                    f"grammar later."
                )
            yield Static(msg, id="gec-msg", markup=True)
            with Horizontal(id="gec-btns"):
                yield Button("Cancel", id="btn-gec-no", variant="default")
                yield Button("Edit",   id="btn-gec-yes", variant="warning")

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-gec-no", Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-gec-no")
    def _no(self, _): self._dismiss_once(False)

    @on(Button.Pressed, "#btn-gec-yes")
    def _yes(self, _): self._dismiss_once(True)

    def action_cancel(self): self._dismiss_once(False)


class _ConfirmDeleteGrammarModal(ModalScreen):
    """Confirmation prompt before `GrammarManagerModal._delete`
    actually wipes a custom grammar + its entry-vector bindings.
    Same shape as `EditGrammarConfirmModal` but with a red Delete
    button and a body that quantifies the blast radius (how many
    parts reference this grammar, what cleanup will run). Default
    focus on Cancel — protects against an accidental Enter on the
    bare Delete button in the manager."""

    # Sweep #26: [INV-50] destructive confirm — Ctrl+Z above
    # could race the impending delete cascade.
    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, grammar_name: str, body_markup: str) -> None:
        super().__init__()
        self._grammar_name = grammar_name
        self._body = body_markup
        self._dismissed: bool = False   # [INV-50]

    def _dismiss_once(self, payload) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="gec-dlg"):
            yield Static(" Delete Grammar? ", id="gec-title")
            yield Static(self._body, id="gec-msg", markup=True)
            with Horizontal(id="gec-btns"):
                yield Button("Cancel", id="btn-gec-no",
                             variant="default")
                yield Button("Delete", id="btn-gec-yes",
                             variant="error")

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-gec-no", Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-gec-no")
    def _no(self, _): self._dismiss_once(False)

    @on(Button.Pressed, "#btn-gec-yes")
    def _yes(self, _): self._dismiss_once(True)

    def action_cancel(self): self._dismiss_once(False)


class ExperimentDeleteConfirmModal(ModalScreen):
    """Yes / No confirmation for entry deletion. Defaults focus on No
    so a stray Enter can't delete an entry. Esc → No."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next button", show=False),
    ]

    def __init__(self, title: str, body: str) -> None:
        super().__init__()
        self._title = title
        self._body  = body
        self._dismissed: bool = False   # [INV-50]

    def _dismiss_once(self, payload: bool) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="expdel-dlg"):
            yield Static(f" {self._title} ", id="expdel-title")
            yield Static(self._body, id="expdel-msg", markup=True)
            with Horizontal(id="expdel-btns"):
                yield Button("No", id="btn-expdel-no", variant="default")
                yield Button("Yes, delete", id="btn-expdel-yes",
                              variant="error")

    DEFAULT_CSS = """
    #expdel-dlg {
        width: 60; height: auto;
        background: $surface; border: solid $primary-darken-2;
        padding: 1 2; layout: vertical;
    }
    #expdel-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; height: 1; text-align: center;
    }
    #expdel-msg { height: auto; margin-top: 1; color: $text; }
    #expdel-btns { height: 3; margin-top: 1; align: center middle; }
    #expdel-btns Button { margin: 0 1; min-width: 14; }
    """

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-expdel-no", Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-expdel-no")
    def _no(self, _ev) -> None:
        self._dismiss_once(False)

    @on(Button.Pressed, "#btn-expdel-yes")
    def _yes(self, _ev) -> None:
        self._dismiss_once(True)

    def action_cancel(self) -> None:
        self._dismiss_once(False)


class ExperimentUnsavedChangesModal(ModalScreen):
    """Three-way prompt shown when the user tries to leave the
    `ExperimentsScreen` with a dirty compose buffer.

    Dismiss payload:
      ``"save"``    — save then exit
      ``"abandon"`` — exit without saving
      ``"cancel"``  — stay on screen (default; Esc + initial focus)
    """

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #expunsaved-dlg {
        width: 64; height: auto;
        background: $surface; border: solid $warning;
        padding: 1 2; layout: vertical;
    }
    #expunsaved-title {
        background: $warning; color: $text;
        padding: 0 1; height: 1; text-style: bold; text-align: center;
    }
    #expunsaved-msg { height: auto; margin-top: 1; color: $text; }
    #expunsaved-btns { height: 3; margin-top: 1; align: center middle; }
    #expunsaved-btns Button { margin: 0 1; min-width: 18; }
    """

    def __init__(self) -> None:
        super().__init__()
        self._dismissed: bool = False   # [INV-50]

    def _dismiss_once(self, payload: str) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="expunsaved-dlg"):
            yield Static(" Unsaved Changes ", id="expunsaved-title")
            yield Static(
                "The current entry has [b]unsaved changes[/b]. "
                "What would you like to do?",
                id="expunsaved-msg", markup=True,
            )
            with Horizontal(id="expunsaved-btns"):
                yield Button("Close", id="btn-expunsaved-cancel",
                             variant="default")
                yield Button("Save changes",
                             id="btn-expunsaved-save",
                             variant="primary")
                yield Button("Abandon and exit",
                             id="btn-expunsaved-abandon",
                             variant="error")

    def on_mount(self) -> None:
        # Default focus on Cancel/Close so a stray Enter can't
        # discard work or auto-save unintended changes.
        try:
            self.query_one(
                "#btn-expunsaved-cancel", Button,
            ).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-expunsaved-cancel")
    def _cancel_btn(self, _ev) -> None:
        self._dismiss_once("cancel")

    @on(Button.Pressed, "#btn-expunsaved-save")
    def _save_btn(self, _ev) -> None:
        self._dismiss_once("save")

    @on(Button.Pressed, "#btn-expunsaved-abandon")
    def _abandon_btn(self, _ev) -> None:
        self._dismiss_once("abandon")

    def action_cancel(self) -> None:
        self._dismiss_once("cancel")


class SpellcheckModal(_OneShotDismissScreen, ModalScreen):
    """Spellcheck sidebar for the active entry's body. Lists
    misspellings + suggestions; per-row buttons replace, add-to-
    dictionary, or skip.

    Dismiss payload: ``(replacements, dict_additions)`` where
    ``replacements`` is a ``{original_word: chosen_replacement}`` map
    and ``dict_additions`` is a list of words to append to the user's
    custom dictionary. The caller applies both — this modal is purely
    UX, no persistence side-effects.
    """

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "done",            "Done"),
        Binding("tab",    "app.focus_next",  "Next", show=False),
    ]

    DEFAULT_CSS = """
    #spell-dlg {
        width: 100; height: 36; max-width: 95%; max-height: 90%;
        background: $surface; border: solid $primary-darken-2;
        padding: 1 2; layout: vertical;
    }
    #spell-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; height: 1; text-align: center;
    }
    #spell-hint   { height: 1; color: $text-muted; margin-top: 1; }
    #spell-table  {
        height: 1fr; min-height: 8; margin-top: 1;
        border: solid $primary-darken-2;
    }
    #spell-suggestions { height: 1; color: $accent; margin-top: 1; }
    #spell-btns { height: 3; margin-top: 1; align: right middle; }
    #spell-btns Button { margin-right: 1; min-width: 14; }
    """

    def __init__(self,
                  misspellings: "list[tuple[str, list[str]]]") -> None:
        super().__init__()
        self._misspellings = list(misspellings)
        self._replacements: "dict[str, str]" = {}
        self._dict_additions: "list[str]" = []

    def compose(self) -> ComposeResult:
        with Vertical(id="spell-dlg"):
            yield Static(
                f" Spellcheck — {len(self._misspellings)} "
                f"unknown words ",
                id="spell-title",
            )
            yield Static(
                "[dim]Click a row to see suggestions. Replace applies "
                "the first suggestion; Add-to-dictionary keeps the word "
                "for future runs.[/]",
                id="spell-hint", markup=True,
            )
            yield DataTable(id="spell-table",
                              cursor_type="row",
                              zebra_stripes=True)
            yield Static(
                "[dim]Pick a row to see suggestions.[/]",
                id="spell-suggestions", markup=True,
            )
            with Horizontal(id="spell-btns"):
                yield Button("Replace", id="btn-spell-replace",
                              variant="primary", disabled=True)
                yield Button("Add to dictionary",
                              id="btn-spell-add", disabled=True)
                yield Button("Skip", id="btn-spell-skip", disabled=True)
                yield Button("Done", id="btn-spell-done")

    def on_mount(self) -> None:
        try:
            t = self.query_one("#spell-table", DataTable)
            t.add_columns("Word", "Suggestions", "Status")
        except NoMatches:
            return
        for word, sugs in self._misspellings:
            sug_str = ", ".join(sugs[:3]) if sugs else "—"
            t.add_row(
                word,
                Text(sug_str, no_wrap=True, overflow="ellipsis"),
                "",
                key=word,
            )
        if t.row_count > 0:
            try:
                t.focus()
            except Exception:
                pass

    @on(DataTable.RowSelected, "#spell-table")
    def _on_row(self, event) -> None:
        if event.row_key is None:
            return
        word = event.row_key.value or ""
        sugs = next(
            (s for w, s in self._misspellings if w == word), [],
        )
        try:
            sug_widget = self.query_one("#spell-suggestions", Static)
        except NoMatches:
            return
        if sugs:
            sug_widget.update(
                f"[accent]Suggestions for '{word}':[/] "
                f"{', '.join(sugs[:5])}",
            )
        else:
            sug_widget.update(
                f"[dim]No suggestions for '{word}'.[/]",
            )
        for bid in ("btn-spell-replace", "btn-spell-add",
                    "btn-spell-skip"):
            try:
                self.query_one(f"#{bid}", Button).disabled = False
            except NoMatches:
                pass

    @on(Button.Pressed, "#btn-spell-replace")
    def _replace(self, _ev) -> None:
        word, sugs = self._get_selected()
        if not word:
            return
        if not sugs:
            self.app.notify(
                "No suggestion to replace with.", severity="warning",
            )
            return
        self._replacements[word] = sugs[0]
        self._mark_done(word, f"→ {sugs[0]}")

    @on(Button.Pressed, "#btn-spell-add")
    def _add(self, _ev) -> None:
        word, _sugs = self._get_selected()
        if not word:
            return
        self._dict_additions.append(word.lower())
        self._mark_done(word, "+ dict")

    @on(Button.Pressed, "#btn-spell-skip")
    def _skip(self, _ev) -> None:
        word, _ = self._get_selected()
        if not word:
            return
        self._mark_done(word, "skipped")

    @on(Button.Pressed, "#btn-spell-done")
    def _done(self, _ev) -> None:
        self.action_done()

    def action_done(self) -> None:
        self.dismiss((self._replacements, self._dict_additions))

    def _get_selected(self) -> "tuple[str, list[str]]":
        try:
            t = self.query_one("#spell-table", DataTable)
        except NoMatches:
            return ("", [])
        if t.cursor_row < 0 or t.row_count == 0:
            return ("", [])
        try:
            row_key = t.coordinate_to_cell_key(
                _Coordinate(t.cursor_row, 0),
            ).row_key
            word = row_key.value or ""
        except Exception:
            return ("", [])
        sugs = next(
            (s for w, s in self._misspellings if w == word), [],
        )
        return (word, sugs)

    def _mark_done(self, word: str, status: str) -> None:
        try:
            t = self.query_one("#spell-table", DataTable)
            for i in range(t.row_count):
                try:
                    row_key = t.coordinate_to_cell_key(
                        _Coordinate(i, 0),
                    ).row_key
                except Exception:
                    continue
                if row_key.value == word:
                    try:
                        t.update_cell_at(_Coordinate(i, 2), status)
                    except Exception:
                        pass
                    if i + 1 < t.row_count:
                        try:
                            t.move_cursor(row=i + 1)
                        except Exception:
                            pass
                    break
        except NoMatches:
            pass


class SynthesisUnsavedChangesModal(ModalScreen):
    """Save / Abandon / Cancel prompt when leaving a dirty synthesis
    buffer. Mirrors the ExperimentUnsavedChangesModal pattern."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    DEFAULT_CSS = """
    #suc-dlg {
        width: 60; height: 12;
        background: $surface; padding: 1 2;
        border: solid $primary-darken-2;
    }
    #suc-title {
        background: $warning; color: black;
        padding: 0 1; height: 1; text-align: center;
    }
    #suc-body { height: 1fr; margin-top: 1; }
    #suc-btns { height: 3; align: right middle; }
    #suc-btns Button { margin-left: 1; min-width: 12; }
    """

    def __init__(self) -> None:
        super().__init__()
        self._dismissed: bool = False   # [INV-50]

    def _dismiss_once(self, payload) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="suc-dlg"):
            yield Static(" Unsaved fragment ", id="suc-title")
            yield Static(
                "The current fragment has unsaved edits. "
                "Save before continuing?",
                id="suc-body",
            )
            with Horizontal(id="suc-btns"):
                yield Button("Save", id="btn-suc-save", variant="primary")
                yield Button("Abandon", id="btn-suc-abandon",
                              variant="error")
                yield Button("Cancel", id="btn-suc-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-suc-cancel", Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-suc-save")
    def _save(self, _) -> None:
        self._dismiss_once("save")

    @on(Button.Pressed, "#btn-suc-abandon")
    def _abandon(self, _) -> None:
        self._dismiss_once("abandon")

    @on(Button.Pressed, "#btn-suc-cancel")
    def _cancel_btn(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)


class SynthesisReplaceDnaConfirmModal(ModalScreen):
    """Confirm before an 'Optimize → DNA' hand-off overwrites unsaved
    DNA-tab edits. Dismisses ``True`` (replace) or ``False`` (cancel)."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    DEFAULT_CSS = """
    #srd-dlg {
        width: 62; height: auto;
        background: $surface; padding: 1 2;
        border: solid $warning;
    }
    #srd-title {
        background: $warning-darken-2; color: $text; padding: 0 1; height: 1; text-align: center;
    }
    #srd-body { margin: 1 0; }
    #srd-btns { height: 3; align: right middle; }
    #srd-btns Button { margin-left: 1; min-width: 12; }
    """

    def __init__(self) -> None:
        super().__init__()
        self._dismissed: bool = False   # [INV-50]

    def _dismiss_once(self, payload) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="srd-dlg"):
            yield Static(" Replace DNA fragment? ", id="srd-title")
            yield Static(
                "The DNA tab has unsaved edits. Optimizing this protein "
                "will replace them with the new CDS.",
                id="srd-body",
            )
            with Horizontal(id="srd-btns"):
                yield Button("Replace", id="btn-srd-replace",
                              variant="warning")
                yield Button("Cancel", id="btn-srd-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-srd-cancel", Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-srd-replace")
    def _replace(self, _) -> None:
        self._dismiss_once(True)

    @on(Button.Pressed, "#btn-srd-cancel")
    def _cancel_btn(self, _) -> None:
        self._dismiss_once(False)

    def action_cancel(self) -> None:
        self._dismiss_once(False)


class NewMotifModal(_OneShotDismissScreen, ModalScreen):
    """Add a new custom protein motif. Dismisses with
    ``{name, sequence, feature_type, color, description}`` on Save (the
    caller persists via `_protein_motif_upsert`), or None on Cancel."""

    _blocks_undo: bool = True
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    DEFAULT_CSS = """
    #nm-dlg { width: 72; height: auto; max-height: 90%;
              background: $surface; border: solid $primary; padding: 1 2; }
    #nm-title { background: $primary-darken-2; color: $text; padding: 0 1;
                margin-bottom: 1; text-align: center; text-style: bold; }
    #nm-dlg Static.nm-label { color: $text-muted; margin-top: 1; }
    #nm-desc { height: 5; border: solid $primary-darken-2; }
    #nm-status { height: auto; margin-top: 1; color: $text-muted; }
    #nm-btns { height: 3; margin-top: 1; align: right middle; }
    #nm-btns Button { margin-right: 1; min-width: 12; }
    """

    def compose(self) -> ComposeResult:
        with Vertical(id="nm-dlg"):
            yield Static(" New protein motif ", id="nm-title")
            yield Static("Name", classes="nm-label")
            yield Input(placeholder="e.g. My tag", id="nm-name")
            yield Static("Amino-acid sequence", classes="nm-label")
            yield Input(placeholder="e.g. HHHHHH", id="nm-seq")
            yield Static("Type", classes="nm-label")
            yield Input(value="Motif", id="nm-type")
            yield Static("Color (hex, optional)", classes="nm-label")
            yield Input(placeholder="#1E40AF", id="nm-color")
            yield Static("Description (optional)", classes="nm-label")
            yield TextArea(id="nm-desc")
            yield Static("", id="nm-status", markup=True)
            with Horizontal(id="nm-btns"):
                yield Button("Save", id="btn-nm-save", variant="primary")
                yield Button("Cancel", id="btn-nm-cancel")

    def _status(self, msg: str) -> None:
        try:
            self.query_one("#nm-status", Static).update(msg)
        except NoMatches:
            pass

    def _val(self, wid: str) -> str:
        try:
            return self.query_one(wid, Input).value.strip()
        except NoMatches:
            return ""

    @on(Button.Pressed, "#btn-nm-save")
    def _save(self, _) -> None:
        name = self._val("#nm-name")
        seq = self._val("#nm-seq").upper()
        if not name:
            self._status("[red]Name is required.[/red]")
            return
        if not seq:
            self._status("[red]Amino-acid sequence is required.[/red]")
            return
        bad = sorted({c for c in seq if c not in "ACDEFGHIKLMNPQRSTVWY*"})
        if bad:
            self._status(
                f"[red]Non-canonical amino acids: {''.join(bad)}[/red]")
            return
        try:
            desc = self.query_one("#nm-desc", TextArea).text.strip()
        except NoMatches:
            desc = ""
        self.dismiss({
            "name":         name,
            "sequence":     seq,
            "feature_type": self._val("#nm-type") or "Motif",
            "color":        self._val("#nm-color"),
            "description":  desc,
        })

    @on(Button.Pressed, "#btn-nm-cancel")
    def _cancel(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class AminoAcidPickerModal(_OneShotDismissScreen, ModalScreen):
    """Tiny picker shown when the user clicks an AA in the Mutagenize
    preview. Returns the selected one-letter AA on dismiss, or None on
    cancel. The WT amino at the clicked position is filtered out so
    the user can't pick a no-op mutation."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    # 20 standard amino acids + stop. Ordered alphabetically by one-letter.
    _AA_CATALOG: list[tuple[str, str, str]] = [
        ("A", "Ala", "Alanine"),       ("C", "Cys", "Cysteine"),
        ("D", "Asp", "Aspartate"),     ("E", "Glu", "Glutamate"),
        ("F", "Phe", "Phenylalanine"), ("G", "Gly", "Glycine"),
        ("H", "His", "Histidine"),     ("I", "Ile", "Isoleucine"),
        ("K", "Lys", "Lysine"),        ("L", "Leu", "Leucine"),
        ("M", "Met", "Methionine"),    ("N", "Asn", "Asparagine"),
        ("P", "Pro", "Proline"),       ("Q", "Gln", "Glutamine"),
        ("R", "Arg", "Arginine"),      ("S", "Ser", "Serine"),
        ("T", "Thr", "Threonine"),     ("V", "Val", "Valine"),
        ("W", "Trp", "Tryptophan"),    ("Y", "Tyr", "Tyrosine"),
        ("*", "***", "Stop codon"),
    ]

    def __init__(self, position: int, wt_aa: str) -> None:
        super().__init__()
        self._position = position
        self._wt_aa    = (wt_aa or "").upper()
        self._choices: list[str] = [
            a for (a, _, _) in self._AA_CATALOG if a != self._wt_aa
        ]

    def compose(self) -> ComposeResult:
        with Vertical(id="aa-pick-box"):
            yield Static(f" Mutate {self._wt_aa}{self._position}  →  ? ",
                         id="aa-pick-title")
            yield Label("[dim]Pick the replacement amino acid. "
                        "Esc to cancel.[/dim]", markup=True)
            items: list = []
            for (a, tl, fn) in self._AA_CATALOG:
                if a == self._wt_aa:
                    continue
                items.append(ListItem(Label(
                    f"[bold]{a}[/bold]   {tl}   [dim]{fn}[/dim]",
                    markup=True,
                )))
            yield ListView(*items, id="aa-pick-list")
            with Horizontal(id="aa-pick-btns"):
                yield Button("Cancel  [Esc]", id="btn-aa-pick-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#aa-pick-list", ListView).focus()
        except NoMatches:
            pass

    @on(ListView.Selected, "#aa-pick-list")
    def _selected(self, _event) -> None:
        lv = self.query_one("#aa-pick-list", ListView)
        if lv.index is None or lv.index >= len(self._choices):
            return
        self.dismiss(self._choices[lv.index])

    @on(Button.Pressed, "#btn-aa-pick-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class PrimerSaveModal(ModalScreen):
    """Prompt for per-oligo names + target primer collection after a
    Primer3 design run produces results.

    2026-05-21: replaces the two-step "fill name inputs then click
    Save to Library" flow. Pressing Design now auto-opens this modal
    so naming + collection routing happens in one place.

    Each oligo in ``oligos`` gets its own name `Input` (scrollable
    when N > visible rows, so a 6-primer mutagenesis set still fits).
    The collection `Select` lists every existing primer collection
    plus a "+ New collection…" sentinel that prompts the user for a
    fresh name. Defaults to the currently-active primer collection.

    Sacred — preserves user-typed names verbatim (no underscore-for-
    space substitution). See feedback_no_underscores_in_names.

    Dismiss payload:
      ``dict`` with keys ``names: list[str]`` and ``collection: str`` on Save.
      ``None`` on Cancel.
    """

    _blocks_undo: bool = True   # carries Inputs; pitfall #41 / [INV-41]

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    _NEW_COLLECTION_SENTINEL = "__new__"

    def __init__(self, oligos: "list[dict]",
                 *, default_collection: str = "") -> None:
        super().__init__()
        # Each oligo dict carries at minimum: default_name (str), label
        # (e.g. "Forward", "Reverse", "Inner Fwd"). Optional: sequence,
        # tm — shown in the right-hand hint column so the user knows
        # which oligo they're naming.
        self._oligos = oligos
        self._default_collection = default_collection or ""
        self._dismissed = False

    def compose(self) -> ComposeResult:
        with Vertical(id="primer-save-dlg"):
            yield Static(" Save primers to library ",
                         id="primer-save-title")
            yield Label(
                f"Name each of the {len(self._oligos)} oligo"
                f"{'s' if len(self._oligos) != 1 else ''} below, "
                "then pick the primer collection. Spaces in names "
                "are preserved.",
                id="primer-save-help",
            )
            with VerticalScroll(id="primer-save-list"):
                for i, ol in enumerate(self._oligos):
                    label = str(ol.get("label") or f"Oligo {i+1}")
                    default = str(ol.get("default_name") or label)
                    tm = ol.get("tm")
                    seq = ol.get("sequence") or ""
                    hint_parts = []
                    if isinstance(tm, (int, float)):
                        hint_parts.append(f"Tm {tm:.1f}°C")
                    if seq:
                        hint_parts.append(f"{len(seq)} nt")
                    hint = " · ".join(hint_parts)
                    with Vertical(classes="primer-save-row"):
                        yield Label(f"[bold]{label}[/]"
                                    + (f"  [dim]{hint}[/]" if hint else ""),
                                    classes="primer-save-row-label",
                                    markup=True)
                        yield Input(
                            value=default,
                            placeholder=label,
                            id=f"primer-save-name-{i}",
                            classes="primer-save-name-input",
                        )
            yield Label("Save to primer collection:",
                        id="primer-save-coll-label")
            yield Select(
                options=self._collection_options(),
                value=self._default_collection
                      or Select.BLANK,
                id="primer-save-collection",
                allow_blank=True,
            )
            with Horizontal(id="primer-save-buttons"):
                yield Button("Save", id="primer-save-ok",
                             variant="primary")
                yield Button("Cancel", id="primer-save-cancel",
                             variant="default")

    def _collection_options(self) -> "list[tuple[str, str]]":
        """Build Select options: every existing primer collection +
        a "+ New collection…" sentinel that triggers a fresh-name
        sub-prompt on Save."""
        opts: list[tuple[str, str]] = []
        for c in _load_primer_collections():
            name = c.get("name")
            if isinstance(name, str) and name:
                opts.append((name, name))
        opts.append(("+ New collection…", self._NEW_COLLECTION_SENTINEL))
        return opts

    def _gather_names(self) -> "list[str] | None":
        names: list[str] = []
        for i in range(len(self._oligos)):
            try:
                inp = self.query_one(f"#primer-save-name-{i}", Input)
            except NoMatches:
                return None
            # SACRED: do NOT strip() or replace internal whitespace —
            # the user gets exactly what they typed. Leading/trailing
            # spaces are still trimmed via .strip() at the boundary
            # because they're almost always typos, not intent.
            val = inp.value.strip()
            if not val:
                self.app.notify(
                    f"Oligo {i+1}: name cannot be empty.",
                    severity="warning",
                )
                inp.focus()
                return None
            names.append(val)
        # Duplicate-name check within THIS batch — saving "F" + "F"
        # would collide downstream; surface here so the user can fix
        # it before the collision resolver fires.
        seen: set[str] = set()
        for n in names:
            if n in seen:
                self.app.notify(
                    f"Duplicate name '{n}' in this batch — "
                    "give each oligo a unique name.",
                    severity="warning",
                )
                return None
            seen.add(n)
        return names

    def _resolve_collection(self) -> "str | None":
        """Return the chosen primer-collection name, or None if the
        user cancelled the new-collection sub-prompt. Sentinel value
        triggers an inline NewCollectionNameModal."""
        try:
            sel = self.query_one("#primer-save-collection", Select)
        except NoMatches:
            return None
        val = sel.value
        if val == Select.BLANK or val is None:
            self.app.notify("Pick a primer collection.", severity="warning")
            return None
        if isinstance(val, str):
            return val
        return None

    @on(Button.Pressed, "#primer-save-ok")
    def _save_btn(self, _) -> None:
        names = self._gather_names()
        if names is None:
            return
        coll = self._resolve_collection()
        if coll is None:
            return
        if coll == self._NEW_COLLECTION_SENTINEL:
            # Inline sub-prompt: ask for a new collection name, then
            # commit. The new collection is created at commit time
            # (not now) so a cancel here doesn't litter empty
            # collections.
            def _on_new_name(new_name) -> None:
                if not isinstance(new_name, str) or not new_name.strip():
                    return
                self._commit(names, new_name.strip(), create=True)
            self.app.push_screen(
                _PrimerCollectionNameModal(taken=self._existing_coll_names()),
                callback=_on_new_name,
            )
            return
        self._commit(names, coll, create=False)

    def _existing_coll_names(self) -> "set[str]":
        out: set[str] = set()
        for c in _load_primer_collections():
            name = c.get("name")
            if isinstance(name, str):
                out.add(name)
        return out

    def _commit(self, names: "list[str]", collection: str,
                *, create: bool) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss({
            "names":      names,
            "collection": collection,
            "create":     create,
        })

    @on(Button.Pressed, "#primer-save-cancel")
    def _cancel_btn(self, _) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(None)

    def action_cancel(self) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(None)


class _PrimerCollectionNameModal(ModalScreen):
    """Sub-prompt: ask the user for a new primer-collection name when
    they pick "+ New collection…" in the save modal. Light-weight
    inline prompt — full collection-management lives elsewhere."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    DEFAULT_CSS = """
    _PrimerCollectionNameModal { align: center middle; }
    #primer-coll-name-dlg { width: 60; height: auto; background: $surface;
                            border: solid $primary; padding: 1 2; }
    #primer-coll-name-title { background: $primary-darken-2; color: $text;
                              padding: 0 1; margin-bottom: 1; text-align: center; }
    #primer-coll-name-input { margin-bottom: 1; }
    #primer-coll-name-status { min-height: 1; margin-bottom: 1; }
    #primer-coll-name-buttons { height: 3; align: right middle; }
    #primer-coll-name-buttons Button { margin-left: 2; }
    """

    def __init__(self, *, taken: "set[str]") -> None:
        super().__init__()
        self._taken = taken
        self._dismissed = False

    def compose(self) -> ComposeResult:
        with Vertical(id="primer-coll-name-dlg"):
            yield Static(" New primer collection ",
                         id="primer-coll-name-title")
            yield Label("Name (spaces allowed):")
            yield Input(placeholder="e.g. 'PCR validation primers'",
                        id="primer-coll-name-input")
            yield Static("", id="primer-coll-name-status", markup=True)
            with Horizontal(id="primer-coll-name-buttons"):
                yield Button("Create", id="primer-coll-name-ok",
                             variant="primary")
                yield Button("Cancel", id="primer-coll-name-cancel",
                             variant="default")

    @on(Input.Changed, "#primer-coll-name-input")
    def _on_changed(self, event: Input.Changed) -> None:
        val = (event.value or "").strip()
        try:
            status = self.query_one("#primer-coll-name-status", Static)
        except NoMatches:
            return
        if not val:
            status.update("")
        elif val in self._taken:
            status.update("[yellow]Name already in use.[/]")
        else:
            status.update("[green]OK[/]")

    @on(Button.Pressed, "#primer-coll-name-ok")
    def _ok(self, _) -> None:
        if self._dismissed:
            return
        try:
            val = self.query_one("#primer-coll-name-input", Input).value.strip()
        except NoMatches:
            val = ""
        if not val:
            self.app.notify("Enter a name.", severity="warning")
            return
        if val in self._taken:
            self.app.notify("Name already in use.", severity="warning")
            return
        self._dismissed = True
        self.dismiss(val)

    @on(Button.Pressed, "#primer-coll-name-cancel")
    def _cancel_btn(self, _) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(None)

    def action_cancel(self) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(None)


class PrimerTypePickerModal(ModalScreen):
    """Pick which mode the Primer Designer should open in.

    Pushed from `action_open_primer_design`; on dismiss with one
    of `{"detection", "cloning", "goldenbraid", "generic"}`, the
    designer opens with that tab active. None on Cancel.
    """

    _blocks_undo: bool = True

    DEFAULT_CSS = """
    #primer-type-dlg {
        width: 64;
        height: auto;
        max-height: 90%;
        background: $surface;
        border: solid $accent;
        padding: 1 2;
    }
    #primer-type-title {
        width: 100%;
        height: 1;
        background: $accent;
        color: $surface;
        text-style: bold;
        content-align: center middle;
        padding: 0 1;
    }
    #primer-type-help {
        width: 100%;
        margin: 1 0;
        color: $text-muted;
        padding: 0 1;
    }
    #primer-type-buttons {
        width: 100%;
        height: auto;
        align: center middle;
    }
    #primer-type-buttons Button {
        width: 100%;
        margin: 0 0 1 0;
    }
    /* Mode-specific accents — green = create / start (Detection),
       cyan = cloning, orange = Golden Braid (matches the GB
       overhang block in the results legend), gray = Generic /
       Cancel = panel-lighten so it sits below the primaries. */
    #primer-type-detection { background: $success; }
    #primer-type-detection:hover { background: $success-lighten-1; }
    #primer-type-cloning { background: #4CC4FF; color: $surface; }
    #primer-type-cloning:hover { background: #80D8FF; }
    #primer-type-goldenbraid { background: #FFB347; color: $surface; }
    #primer-type-goldenbraid:hover { background: #FFCB80; }
    #primer-type-generic { background: $primary; }
    #primer-type-generic:hover { background: $primary-lighten-1; }
    #primer-type-cancel {
        background: $panel-lighten-2;
        margin-top: 1;
    }
    #primer-type-cancel:hover { background: $panel-lighten-3; }
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    def __init__(self) -> None:
        super().__init__()
        self._dismissed = False

    def compose(self) -> ComposeResult:
        with Vertical(id="primer-type-dlg"):
            yield Static(" What kind of primers? ",
                         id="primer-type-title")
            yield Label(
                "Pick the workflow. The Primer Designer opens "
                "on the matching tab; your selection (if any) "
                "auto-fills the Start/End coordinates — and if "
                "it spans a full feature, the Feature dropdown "
                "is selected so you can press Design directly.",
                id="primer-type-help",
            )
            with Vertical(id="primer-type-buttons"):
                yield Button("Detection  —  diagnostic PCR",
                             id="primer-type-detection",
                             variant="primary")
                yield Button("Cloning  —  RE tails + GCGC flap",
                             id="primer-type-cloning")
                yield Button("Golden Braid  —  L0 domestication",
                             id="primer-type-goldenbraid")
                yield Button("Generic  —  binding only",
                             id="primer-type-generic")
                yield Button("Cancel", id="primer-type-cancel",
                             variant="default")

    @on(Button.Pressed, "#primer-type-detection")
    def _pick_det(self, _) -> None:
        self._commit("detection")

    @on(Button.Pressed, "#primer-type-cloning")
    def _pick_clo(self, _) -> None:
        self._commit("cloning")

    @on(Button.Pressed, "#primer-type-goldenbraid")
    def _pick_gb(self, _) -> None:
        self._commit("goldenbraid")

    @on(Button.Pressed, "#primer-type-generic")
    def _pick_gen(self, _) -> None:
        self._commit("generic")

    @on(Button.Pressed, "#primer-type-cancel")
    def _cancel_btn(self, _) -> None:
        self.action_cancel()

    def _commit(self, mode: str) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(mode)

    def action_cancel(self) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(None)


class _PrimerCollectionDeleteConfirmModal(_OneShotDismissScreen, ModalScreen):
    """Confirm deleting a PRIMER collection. Default focus on [No] (handslip
    guard, like every destructive confirm in the app). Dismisses True (delete)
    or False (keep)."""

    _blocks_undo: bool = True   # [INV-50] destructive confirm
    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab", "app.focus_next", "Next", show=False),
    ]
    DEFAULT_CSS = """
    _PrimerCollectionDeleteConfirmModal { align: center middle; }
    #pcolldel-dlg { width: 66; height: auto; background: $surface;
                    border: solid $error; padding: 1 2; }
    #pcolldel-title { background: $error; color: $text; padding: 0 1;
                      margin-bottom: 1; text-align: center; }
    #pcolldel-btns { height: 3; align: right middle; margin-top: 1; }
    #pcolldel-btns Button { margin-left: 2; }
    """

    def __init__(self, name: str, n_primers: int) -> None:
        super().__init__()
        self._name = name
        self._n = n_primers

    def compose(self) -> ComposeResult:
        from rich.markup import escape as _esc
        plural = "" if self._n == 1 else "s"
        with Vertical(id="pcolldel-dlg"):
            yield Static(" Delete primer collection ", id="pcolldel-title")
            yield Static(
                f"  Delete primer collection [bold]{_esc(self._name)}[/bold]?\n"
                f"  ({self._n} primer{plural})\n\n"
                f"  [dim]A backup is written to\n"
                f"  primer_collections.json.bak before the change.[/dim]",
                id="pcolldel-msg", markup=True)
            with Horizontal(id="pcolldel-btns"):
                yield Button("No", id="btn-pcolldel-no", variant="default")
                yield Button("Yes, delete", id="btn-pcolldel-yes",
                             variant="error")

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-pcolldel-no", Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-pcolldel-no")
    def _no(self, _) -> None:
        self.dismiss(False)

    @on(Button.Pressed, "#btn-pcolldel-yes")
    def _yes(self, _) -> None:
        self.dismiss(True)

    def action_cancel(self) -> None:
        self.dismiss(False)


class _PrimerMoveCopyModal(_OneShotDismissScreen, ModalScreen):
    """Pick a destination container to move / copy items into. Dismisses
    with the chosen container name (str) or None on cancel.

    Originally primer-collection-only; now a generic destination picker
    reused for parts→bin and notebook-entry→project moves via the
    ``item_singular`` / ``item_plural`` / ``dest_label`` params (all
    default to the primer-collection wording, so existing callers are
    unchanged). The destination rows are passed in by the caller already
    filtered (e.g. every container except the active one).

    Layout: the destination table and the Cancel button sit in a padded
    box with a 1-char void on every side; the row picks on Enter or click,
    so the picker stays a single uncluttered column."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab", "app.focus_next", "Next", show=False),
    ]
    DEFAULT_CSS = """
    _PrimerMoveCopyModal { align: center middle; }
    #pmove-box { width: 60; height: auto; max-height: 80%;
                 background: $surface; border: solid $accent; padding: 1 2; }
    #pmove-title { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; text-align: center; }
    #pmove-table { height: 14; margin-bottom: 1; }
    #pmove-btns { height: 3; align: right middle; }
    #pmove-btns Button { margin: 0 1; min-width: 12; }
    """

    def __init__(self, mode: str, collections: "list[str]",
                 count: int, *, item_singular: str = "primer",
                 item_plural: str = "primers",
                 dest_label: str = "collection") -> None:
        super().__init__()
        self._mode = mode
        self._collections = list(collections or [])
        self._count = count
        self._item_singular = item_singular
        self._item_plural = item_plural
        self._dest_label = dest_label

    def compose(self) -> ComposeResult:
        verb = "Copy" if self._mode == "copy" else "Move"
        noun = self._item_singular if self._count == 1 else self._item_plural
        with Vertical(id="pmove-box"):
            yield Static(f" {verb} {self._count} {noun} to… ",
                         id="pmove-title")
            yield DataTable(id="pmove-table", cursor_type="row")
            with Horizontal(id="pmove-btns"):
                yield Button("Cancel", id="btn-pmove-cancel")

    def on_mount(self) -> None:
        try:
            t = self.query_one("#pmove-table", DataTable)
        except NoMatches:
            return
        t.add_column(f"Destination {self._dest_label}")
        for name in self._collections:
            t.add_row(name, key=name)
        t.focus()

    @on(DataTable.RowSelected, "#pmove-table")
    def _row(self, event: DataTable.RowSelected) -> None:
        rk = event.row_key
        name = rk.value if rk else None
        if isinstance(name, str) and name:
            self.dismiss(name)

    @on(Button.Pressed, "#btn-pmove-cancel")
    def _cancel_btn(self) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class AmpliconSaveModal(ModalScreen):
    """Name a PCR amplicon and pick the collection to save the linear
    fragment into. Dismisses with ``{"name": str, "collection": str}``
    on confirm, or ``None`` on cancel.

    Hardened:
      * ``_blocks_undo`` — the name Input holds focus; a stray Ctrl+Z
        under it would otherwise revert the canvas.
      * One-shot ``_dismissed`` guard — a double-fire ``Pressed`` (real
        terminals can emit two for one physical click) can't dismiss
        twice / double-save.
      * Empty / whitespace name → refused with a status line, never a
        silent save under a blank name.
      * The collection list is guaranteed non-empty: the caller passes
        the active collection, which is de-duped to the front so the
        ``Select`` (``allow_blank=False``) always mounts with a valid
        value even if `collections.json` is empty / unreadable.
    """

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #ampsave-dlg {
        width: 70; height: auto; max-height: 90%;
        background: $surface; border: solid $primary; padding: 1 2;
    }
    #ampsave-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1; text-align: center;
    }
    #ampsave-dlg Label   { color: $text-muted; margin-top: 1; }
    #ampsave-name        { margin-top: 1; }
    #ampsave-collection  { margin-top: 1; }
    #ampsave-status      { height: 1; margin-top: 1; }
    #ampsave-btns        { align: right middle;  height: 3; margin-top: 1; }
    #ampsave-btns Button { margin-right: 1; min-width: 10; }
    """

    def __init__(self, *, default_name: str,
                 collections: "list[str]",
                 active_collection: "str | None",
                 title: str = "Save amplicon to library",
                 name_label: str = "Amplicon name",
                 parts_bins: "list[str] | None" = None,
                 active_bin: "str | None" = None,
                 fragment_default_name: "str | None" = None,
                 fragment_collections: "list[str] | None" = None,
                 fragment_active_collection: "str | None" = None,
                 fragment_label: str = "Fragment name") -> None:
        super().__init__()
        # Title + name-label are parametrised (back-compat defaults) so the
        # same hardened name+collection dialog serves the PCR amplicon save
        # AND the Domesticator cloned-plasmid save (dual-save flow).
        self._title = (title or "Save amplicon to library").strip() \
            or "Save amplicon to library"
        self._name_label = (name_label or "Amplicon name").strip() \
            or "Amplicon name"
        self._default_name = (default_name or "").strip() or "PCR amplicon"
        # Guarantee the active collection is present + first so the
        # Select always has a valid value; de-dup while preserving
        # order. Empty / blank names dropped. Final fallback: "Default".
        names: list[str] = []
        for raw in [active_collection, *collections]:
            nm = (raw or "").strip()
            if nm and nm not in names:
                names.append(nm)
        if not names:
            names = ["Default"]
        self._collections = names
        self._active = names[0]
        # Optional parts-bin picker (dual-save: choose the part's bin in
        # the SAME dialog). None → no bin field, so the amplicon path is
        # byte-unchanged. Normalised like collections: active first,
        # de-duped, never empty.
        if parts_bins is None:
            self._parts_bins: "list[str] | None" = None
            self._active_bin: "str | None" = None
        else:
            bins: list[str] = []
            for raw in [active_bin, *parts_bins]:
                nm = (raw or "").strip()
                if nm and nm not in bins:
                    bins.append(nm)
            self._parts_bins = bins or ["Main Parts Bin"]
            self._active_bin = self._parts_bins[0]
        # Optional primed-fragment fields (dual-save, 2026-06-02): name +
        # place the LINEAR primed fragment independently of the clone.
        # None → no fragment row, so the amplicon / clone-only paths are
        # byte-unchanged. Collections normalised like the clone's: active
        # first, de-duped, never empty.
        self._fragment_label = (fragment_label or "Fragment name").strip() \
            or "Fragment name"
        if fragment_default_name is None:
            self._fragment_default_name: "str | None" = None
            self._fragment_collections: "list[str]" = []
            self._fragment_active: "str | None" = None
        else:
            self._fragment_default_name = (
                (fragment_default_name or "").strip() or "FRAG-fragment")
            f_src = (fragment_collections
                     if fragment_collections is not None else collections)
            f_names: list[str] = []
            for raw in [fragment_active_collection, *(f_src or [])]:
                nm = (raw or "").strip()
                if nm and nm not in f_names:
                    f_names.append(nm)
            self._fragment_collections = f_names or list(self._collections)
            self._fragment_active = self._fragment_collections[0]
        self._dismissed = False

    def _dismiss_once(self, result) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(result)

    def compose(self) -> ComposeResult:
        with Vertical(id="ampsave-dlg"):
            yield Static(f" {self._title} ", id="ampsave-title")
            yield Label(f"{self._name_label}:")
            yield Input(value=self._default_name, id="ampsave-name",
                          placeholder=self._name_label)
            yield Label("Save to collection:")
            yield Select(
                [(n, n) for n in self._collections],
                value=self._active,
                allow_blank=False,
                id="ampsave-collection",
            )
            if self._fragment_default_name is not None:
                yield Label(f"{self._fragment_label}:")
                yield Input(value=self._fragment_default_name,
                            id="ampsave-frag-name",
                            placeholder=self._fragment_label)
                yield Label("Save fragment to collection:")
                yield Select(
                    [(n, n) for n in self._fragment_collections],
                    value=self._fragment_active,
                    allow_blank=False,
                    id="ampsave-frag-collection",
                )
            if self._parts_bins is not None:
                yield Label("Store part in parts bin:")
                yield Select(
                    [(n, n) for n in self._parts_bins],
                    value=self._active_bin,
                    allow_blank=False,
                    id="ampsave-bin",
                )
            yield Static("", id="ampsave-status", markup=True)
            with Horizontal(id="ampsave-btns"):
                yield Button("Save", id="btn-ampsave-ok", variant="primary")
                yield Button("Cancel", id="btn-ampsave-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#ampsave-name", Input).focus()
        except NoMatches:
            pass

    def _set_status(self, msg: str) -> None:
        try:
            self.query_one("#ampsave-status", Static).update(msg)
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-ampsave-ok")
    def _ok(self, _) -> None:
        self._submit()

    @on(Input.Submitted, "#ampsave-name")
    def _submitted(self, _) -> None:
        self._submit()

    def _submit(self) -> None:
        try:
            name = self.query_one("#ampsave-name", Input).value.strip()
            coll_w = self.query_one("#ampsave-collection", Select)
        except NoMatches:
            return
        if not name:
            self._set_status("[red]Enter an amplicon name.[/red]")
            return
        coll = coll_w.value
        if not isinstance(coll, str) or not coll.strip():
            coll = self._active
        result = {"name": name, "collection": coll}
        # Optional primed-fragment fields (dual-save). When the fragment
        # row is present its name must be non-empty too — never save a
        # fragment under a blank name (parity with the clone-name check).
        if self._fragment_default_name is not None:
            try:
                frag_name = self.query_one(
                    "#ampsave-frag-name", Input).value.strip()
            except NoMatches:
                frag_name = ""
            if not frag_name:
                self._set_status("[red]Enter a fragment name.[/red]")
                return
            frag_coll = self._fragment_active
            try:
                fv = self.query_one("#ampsave-frag-collection", Select).value
                if isinstance(fv, str) and fv.strip():
                    frag_coll = fv
            except NoMatches:
                pass
            result["frag_name"] = frag_name
            result["frag_collection"] = frag_coll or ""
        # Optional parts-bin choice (dual-save). Absent → caller keeps
        # the active bin.
        if self._parts_bins is not None:
            chosen_bin = self._active_bin
            try:
                bv = self.query_one("#ampsave-bin", Select).value
                if isinstance(bv, str) and bv.strip():
                    chosen_bin = bv
            except NoMatches:
                pass
            result["bin"] = chosen_bin or ""
        self._dismiss_once(result)

    @on(Button.Pressed, "#btn-ampsave-cancel")
    def _cancel_btn(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)


class UnsavedNavigateModal(ModalScreen):
    """Shown when the user tries to navigate (e.g. Back to Collections)
    with unsaved edits. Sibling of `UnsavedQuitModal` — kept separate
    because the button labels and verb differ ("go back" vs "quit"),
    and the wording matters for users to understand the consequence.

    Dismisses with ``"save"`` (caller saves then proceeds), ``"discard"``
    (caller reverts the in-memory record from the library copy then
    proceeds), or ``None`` (cancel — stay).
    """

    _blocks_undo: bool = True   # [INV-50] destructive confirm — Ctrl+Z above

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "app.focus_next", "Next button", show=False),
    ]

    def __init__(self, action_phrase: str = "leave"):
        super().__init__()
        self._action_phrase = action_phrase
        self._dismissed: bool = False   # [INV-50]

    def _dismiss_once(self, payload) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        with Vertical(id="navunsv-dlg"):
            yield Static(" Unsaved Changes ", id="navunsv-title")
            yield Static(
                f"  The loaded plasmid has unsaved edits.\n"
                f"  Save before you {self._action_phrase}?",
                id="navunsv-msg",
            )
            with Horizontal(id="navunsv-btns"):
                yield Button("Save",            id="btn-navunsv-save",
                             variant="primary")
                yield Button("Discard Changes", id="btn-navunsv-discard",
                             variant="error")
                yield Button("Cancel",          id="btn-navunsv-cancel")

    @on(Button.Pressed, "#btn-navunsv-save")
    def _save(self, _):     self._dismiss_once("save")

    @on(Button.Pressed, "#btn-navunsv-discard")
    def _discard(self, _):  self._dismiss_once("discard")

    @on(Button.Pressed, "#btn-navunsv-cancel")
    def _cancel_btn(self, _): self._dismiss_once(None)

    def action_cancel(self): self._dismiss_once(None)


class AnnotationTransferModal(_OneShotDismissScreen, ModalScreen):
    """Preview + confirm annotation transfer from a source plasmid
    onto the loaded record.

    Built from a list of transfer dicts produced by
    `_find_annotation_transfers` — the table shows feature label /
    target coords / strand / length so the user can sanity-check
    before any features land on their construct. "Apply all" adds
    every listed feature to the current record (the matcher already
    skipped exact-coord duplicates, so re-running on an already-
    annotated target is a no-op for visible duplicates).

    Dismiss payload:
      None       — cancelled
      list[dict] — the transfers the user accepted (same shape the
                   modal received). Caller turns them into SeqFeatures
                   and `_apply_record`s the result.
    """

    # Sweep #26: Apply-all triggers `_apply_record` (canvas mutation
    # + history push). Block app-level Ctrl+Z so a stray keystroke
    # over the preview doesn't roll back the canvas.
    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next",   show=False),
    ]

    def __init__(self, source_label: str, target_label: str,
                 transfers: list[dict]) -> None:
        super().__init__()
        self._source_label = source_label
        self._target_label = target_label
        self._transfers    = list(transfers)

    def compose(self) -> ComposeResult:
        with Vertical(id="annot-box"):
            yield Static(" Transfer annotations ", id="annot-title")
            yield Label(
                f"From: {self._source_label}    →    "
                f"Into: {self._target_label}"
            )
            yield Label(
                f"{len(self._transfers)} feature(s) matched by sequence"
                if self._transfers
                else "No features matched. Sequences too divergent?"
            )
            yield DataTable(id="annot-table", cursor_type="row",
                            zebra_stripes=True)
            with Horizontal(id="annot-btns"):
                yield Button("Apply all", id="btn-annot-apply",
                             variant="primary",
                             disabled=not self._transfers)
                yield Button("Cancel", id="btn-annot-cancel")

    def on_mount(self) -> None:
        try:
            t = self.query_one("#annot-table", DataTable)
        except NoMatches:
            return
        t.add_columns("Label", "Type", "Target start",
                      "Target end", "Strand", "Length (bp)")
        for tr in self._transfers:
            # Display is GenBank-style 1-based inclusive. Internal
            # coords are 0-based half-open: start += 1, end stays
            # (exclusive end == inclusive end numerically). Wrap
            # features (target_end < target_start) get a "(wrap)"
            # tag so the user can tell origin-spanning at a glance.
            wrap = (tr["target_end"] < tr["target_start"])
            start_disp = f"{tr['target_start'] + 1:,}"
            end_disp   = (f"{tr['target_end']:,} (wrap)"
                          if wrap else f"{tr['target_end']:,}")
            t.add_row(
                Text(tr["label"], no_wrap=True, overflow="ellipsis"),
                tr["type"],
                start_disp,
                end_disp,
                "+" if tr["target_strand"] == 1 else "-",
                f"{tr['length']:,}",
            )

    @on(Button.Pressed, "#btn-annot-apply")
    def _apply(self, _) -> None:
        self.dismiss(self._transfers)

    @on(Button.Pressed, "#btn-annot-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class LargeFileConfirmModal(ModalScreen):
    """Confirm before loading a large record (file or NCBI fetch).

    Shown whenever the on-disk size or the parsed sequence length
    crosses the "this is going to take a noticeable while" threshold.
    Default focus is on **No** so a tap of Enter (or Esc) bails — the
    user has to deliberately reach for Yes to commit. Dismisses True
    (load) or False (cancel).

    `description` is a short context line that explains what's about
    to load (file path, accession, member name) so the user can tell
    "yes, this is the chromosome I asked for" from "wait, that's not
    what I expected".
    """

    # Sweep #26: loading a 50+ MB file is destructive to the canvas
    # state — block app-level Ctrl+Z while the user reads the prompt.
    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #lfc-dlg {
        width: 72; height: auto; max-height: 60%;
        background: $surface; border: solid $warning; padding: 1 2;
    }
    #lfc-title { background: $warning-darken-2; color: $text;
                  padding: 0 1; margin-bottom: 1; text-align: center; }
    #lfc-msg   { margin-bottom: 1; }
    #lfc-btns  { align: center middle;  height: 3; margin-top: 1; }
    #lfc-btns Button { margin-right: 1; min-width: 14; }
    """

    def __init__(self, description: str,
                 size_text: str = "",
                 *, threshold_text: str = "") -> None:
        super().__init__()
        self.description    = description
        self.size_text      = size_text
        self.threshold_text = threshold_text
        self._dismissed: bool = False   # [INV-50]

    def _dismiss_once(self, payload: bool) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        from rich.markup import escape as _esc
        with Vertical(id="lfc-dlg"):
            yield Static(" Large record — confirm load ", id="lfc-title")
            body  = f"  About to load:\n  [bold]{_esc(self.description)}[/]\n"
            if self.size_text:
                body += f"\n  Size: {_esc(self.size_text)}\n"
            if self.threshold_text:
                body += f"  [dim]{_esc(self.threshold_text)}[/]\n"
            body += (
                "\n  Loading large records can take several seconds and use "
                "significant memory. Continue?"
            )
            yield Static(body, id="lfc-msg", markup=True)
            with Horizontal(id="lfc-btns"):
                yield Button("No (default)", id="btn-lfc-no",
                              variant="default")
                yield Button("Yes, load",    id="btn-lfc-yes",
                              variant="warning")

    def on_mount(self) -> None:
        # Default focus is No so a stray Enter bails out rather than
        # commits to a multi-second load.
        self.query_one("#btn-lfc-no", Button).focus()

    @on(Button.Pressed, "#btn-lfc-no")
    def _no(self, _) -> None:
        self._dismiss_once(False)

    @on(Button.Pressed, "#btn-lfc-yes")
    def _yes(self, _) -> None:
        self._dismiss_once(True)

    def action_cancel(self) -> None:
        self._dismiss_once(False)


class CollectionDeleteConfirmModal(ModalScreen):
    """Confirm-on-delete modal for collections — different copy from
    LibraryDeleteConfirmModal (which talks about library entries).

    Default focus on [No] to protect against handslip-deletes.
    Dismisses True (delete) or False (keep)."""

    _blocks_undo: bool = True   # [INV-50] destructive confirm — Ctrl+Z above

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, name: str, n_plasmids: int) -> None:
        super().__init__()
        self.coll_name = name
        self.n_plas = n_plasmids
        # [INV-50] dismiss-once guard against double-click / double-Enter
        # races landing two `dismiss()` calls on a destructive modal.
        self._dismissed: bool = False

    def compose(self) -> ComposeResult:
        plural = "" if self.n_plas == 1 else "s"
        with Vertical(id="colldel-dlg"):
            yield Static(" Delete collection ", id="colldel-title")
            yield Static(
                f"  Delete collection [bold]{self.coll_name}[/bold]?\n"
                f"  ({self.n_plas} plasmid{plural})\n\n"
                f"  [dim]A backup is written to\n"
                f"  collections.json.bak before the change.[/dim]",
                id="colldel-msg", markup=True,
            )
            with Horizontal(id="colldel-btns"):
                yield Button("No",          id="btn-colldel-no",  variant="default")
                yield Button("Yes, delete", id="btn-colldel-yes", variant="error")

    def on_mount(self) -> None:
        self.query_one("#btn-colldel-no", Button).focus()

    def _dismiss_once(self, payload: bool) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    @on(Button.Pressed, "#btn-colldel-no")
    def _no(self, _) -> None:
        self._dismiss_once(False)

    @on(Button.Pressed, "#btn-colldel-yes")
    def _yes(self, _) -> None:
        self._dismiss_once(True)

    def action_cancel(self) -> None:
        self._dismiss_once(False)


class ScaryDeleteConfirmModal(ModalScreen):
    """Second-stage confirmation for collection delete — deliberately
    visually loud (red border + warning banner + emphatic copy) to make
    the user pause. Default focus on [No] like every confirm modal in
    the app. Dismisses True (delete) or False (keep)."""

    _blocks_undo: bool = True   # [INV-50] destructive confirm — Ctrl+Z above

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, name: str, n_plasmids: int) -> None:
        super().__init__()
        self.coll_name = name
        self.n_plas = n_plasmids
        # [INV-50] dismiss-once guard
        self._dismissed: bool = False

    def _dismiss_once(self, payload: bool) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        plural = "" if self.n_plas == 1 else "s"
        with Vertical(id="scarydel-dlg"):
            yield Static(
                "  ⚠   ARE YOU ABSOLUTELY SURE?   ⚠  ",
                id="scarydel-title", markup=False,
            )
            yield Static(
                f"\n  This will [bold red]permanently delete[/bold red] the "
                f"collection\n"
                f"  [bold]{self.coll_name}[/bold] and its "
                f"[bold red]{self.n_plas} plasmid{plural}[/bold red].\n\n"
                f"  [yellow]The plasmids inside will also be removed from\n"
                f"  the library mirror.[/yellow]\n\n"
                f"  A backup of [italic]collections.json[/italic] is written "
                f"to\n"
                f"  [italic]collections.json.bak[/italic] in your data "
                f"directory —\n"
                f"  recover from it manually if you change your mind.\n",
                id="scarydel-msg", markup=True,
            )
            with Horizontal(id="scarydel-btns"):
                yield Button("No, keep it", id="btn-scarydel-no",
                             variant="default")
                yield Button("Yes, delete forever", id="btn-scarydel-yes",
                             variant="error")

    def on_mount(self) -> None:
        self.query_one("#btn-scarydel-no", Button).focus()

    @on(Button.Pressed, "#btn-scarydel-no")
    def _no(self, _) -> None:
        self._dismiss_once(False)

    @on(Button.Pressed, "#btn-scarydel-yes")
    def _yes(self, _) -> None:
        self._dismiss_once(True)

    def action_cancel(self) -> None:
        self._dismiss_once(False)


class MasterDeleteModal(_OneShotDismissScreen, ModalScreen):
    """Stage 1 of the Master Delete flow.

    Gates the destructive button behind a typed challenge: the user
    must type exactly ``YES`` (case-sensitive, no whitespace, no
    normalisation) into the Input field before the Delete button
    activates. Cancel is default-focused; Esc → cancel. Dismisses
    ``True`` (proceed to confirm modal) or ``False`` (abort).

    Design contract — this is the LEAST destructive of the two
    Master Delete modals. Even hitting "Delete" here only opens
    Stage 2; nothing on disk is touched until the Stage 2 button
    fires AFTER its 3 s cool-down expires.

    Sacred — must satisfy ALL of:
      * input matches `"YES"` byte-for-byte (case-sensitive); "yes",
        " YES", "YES ", "YESS", "" all keep Delete disabled.
      * Delete button starts disabled — no race where the user can
        click it before the first `Input.Changed` lands.
      * Cancel is the default-focused button so a stray Enter at the
        modal level cancels rather than commits.
      * `_blocks_undo = True` so app-level Ctrl+Z under the modal
        doesn't fire on whatever canvas is loaded underneath.
      * No keyboard binding for "delete" — only the visible button
        with an explicit click commits the user to Stage 2.
    """

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #md-dlg {
        width: 78; height: auto; max-height: 90%;
        background: $surface; border: solid $error; padding: 1 2;
    }
    #md-title {
        background: $error-darken-1; color: $text;
        padding: 0 1; margin-bottom: 1; text-align: center;
    }
    #md-warn { margin-bottom: 1; }
    #md-scope-label { margin-top: 1; color: $warning; }
    #md-scope { color: $text; margin-bottom: 1; }
    #md-prompt { margin-top: 1; color: $warning; }
    #md-input { margin-bottom: 1; }
    #md-status { margin-bottom: 1; }
    #md-btns { align: center middle;  height: 3; margin-top: 1; }
    #md-btns Button { margin-right: 1; min-width: 22; }
    """

    # Exact required input. Sacred — DO NOT loosen to a casefold or
    # `.strip()` compare. A user typing "yes" hasn't deliberately
    # spelled out the affirmative; they've muscle-memoried it.
    _REQUIRED_INPUT = "YES"

    def __init__(self, *,
                  files_count: int,
                  dirs_count: int,
                  pre_update_present: bool) -> None:
        super().__init__()
        # Pre-computed scope counts. Caller (the action handler)
        # runs the enumeration once on the main thread before
        # pushing the modal so the user sees the same numbers in
        # both stages.
        self._files_count = int(files_count)
        self._dirs_count = int(dirs_count)
        self._pre_update_present = bool(pre_update_present)

    def compose(self) -> ComposeResult:
        with Vertical(id="md-dlg"):
            yield Static(
                "  ⚠   MASTER DELETE — WIPE ALL USER DATA   ⚠  ",
                id="md-title", markup=False,
            )
            yield Static(
                "  This will [bold red]permanently and "
                "irreversibly[/bold red] delete every plasmid,\n"
                "  collection, experiment, gel, primer, part, "
                "feature, grammar,\n"
                "  codon table, custom setting and saved "
                "preference in this\n"
                "  SpliceCraft data directory.\n\n"
                "  There is [bold red]no built-in restore[/bold red]. "
                "Pre-update snapshots,\n"
                "  daily backups, .bak files, and lost-entries "
                "spillover\n"
                "  will [bold red]all be wiped[/bold red] too — true "
                "clean slate.",
                id="md-warn", markup=True,
            )
            scope_lines = [
                f"  • {self._files_count} user-data file(s) (plus "
                "every `.bak` sibling)",
                f"  • {self._dirs_count} user-data directory tree(s)",
            ]
            if self._pre_update_present:
                scope_lines.append(
                    "  • Sibling [italic]pre-update-backups[/italic] "
                    "directory (your last recovery copy)"
                )
            scope_lines.append(
                "  • Rotated log backups (active log stays open)"
            )
            yield Static(
                "  Scope of this deletion:", id="md-scope-label",
            )
            yield Static(
                "\n".join(scope_lines),
                id="md-scope", markup=True,
            )
            yield Static(
                f"  To enable the Delete button, type exactly  "
                f"[bold]{self._REQUIRED_INPUT}[/bold]  (case-sensitive)\n"
                f"  in the field below. Anything else keeps the "
                f"button disabled.",
                id="md-prompt", markup=True,
            )
            yield Input(
                value="",
                placeholder='type "YES" (case-sensitive, no spaces) to enable',
                id="md-input",
            )
            yield Static(
                "[dim]Delete button stays disabled until input "
                "matches.[/dim]",
                id="md-status", markup=True,
            )
            with Horizontal(id="md-btns"):
                yield Button(
                    "Cancel (default)", id="btn-md-cancel",
                    variant="default",
                )
                yield Button(
                    "Delete (disabled)", id="btn-md-delete",
                    variant="error", disabled=True,
                )

    def on_mount(self) -> None:
        # Cancel is the default focus so a stray Enter at the modal
        # level fires Cancel, NOT Delete. The Input below grabs
        # focus only after we explicitly Tab into it — see
        # invariant on Cancel-default-focus across destructive modals.
        try:
            self.query_one("#btn-md-cancel", Button).focus()
        except NoMatches:
            pass

    @on(Input.Changed, "#md-input")
    def _on_input_changed(self, event: Input.Changed) -> None:
        # Strict byte-for-byte compare. No `.strip()`, no `.upper()`,
        # no `.casefold()` — the user must type the exact 3 chars.
        # If we ever loosened this it would defeat the whole point
        # of the typed challenge.
        try:
            btn = self.query_one("#btn-md-delete", Button)
            status = self.query_one("#md-status", Static)
        except NoMatches:
            return
        if event.value == self._REQUIRED_INPUT:
            btn.disabled = False
            btn.label = "Delete"
            status.update(
                "[bold green]✓ Match — Delete button enabled.[/bold green] "
                "[dim](next stage has a 3 s cool-down)[/dim]"
            )
        else:
            btn.disabled = True
            btn.label = "Delete (disabled)"
            if event.value == "":
                status.update(
                    "[dim]Delete button stays disabled until "
                    "input matches.[/dim]"
                )
            else:
                # Helpful but not too helpful — tell the user the
                # value didn't match without revealing the exact
                # mismatch (case, whitespace, extra chars).
                status.update(
                    "[red]✗ Doesn't match. Type the exact word "
                    f"[/red][bold]{self._REQUIRED_INPUT}[/bold]"
                    "[red] — case-sensitive, no spaces.[/red]"
                )

    @on(Button.Pressed, "#btn-md-cancel")
    def _cancel_btn(self, _) -> None:
        _log_event("masterdelete.stage1.cancel")
        self.dismiss(False)

    @on(Button.Pressed, "#btn-md-delete")
    def _delete_btn(self, _) -> None:
        # Final guard — even if a stale `Input.Changed` missed a
        # disable, re-check the value here so the only path that
        # dismisses True has a verified match. Pulls the live
        # widget value (not a stashed one) so a fast retype-and-
        # click sequence still validates the visible string.
        try:
            inp = self.query_one("#md-input", Input)
        except NoMatches:
            self.dismiss(False)
            return
        if inp.value != self._REQUIRED_INPUT:
            # Re-disable + status nudge.
            try:
                btn = self.query_one("#btn-md-delete", Button)
                btn.disabled = True
                btn.label = "Delete (disabled)"
            except NoMatches:
                pass
            try:
                self.query_one("#md-status", Static).update(
                    "[red]✗ Input no longer matches — type "
                    f"[/red][bold]{self._REQUIRED_INPUT}[/bold]"
                    "[red] exactly.[/red]"
                )
            except NoMatches:
                pass
            return
        _log_event("masterdelete.stage1.confirmed")
        self.dismiss(True)

    def action_cancel(self) -> None:
        _log_event("masterdelete.stage1.cancel", reason="esc")
        self.dismiss(False)


class MasterDeleteResultModal(_OneShotDismissScreen, ModalScreen):
    """Stage 3 — post-wipe summary + restart nudge. One button
    ("OK, I'll restart") so the user can't accidentally dismiss the
    confirmation while they're still reading. Dismisses ``None``.
    """

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "ok", "OK"),
        Binding("enter",  "ok", "OK", show=False),
    ]

    DEFAULT_CSS = """
    #mdr-dlg {
        width: 70; height: auto; max-height: 70%;
        background: $surface; border: solid $accent; padding: 1 2;
    }
    #mdr-title {
        background: $accent-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1; text-align: center;
    }
    #mdr-msg { margin-bottom: 1; }
    #mdr-btns { align: center middle;  height: 3; margin-top: 1; }
    #mdr-btns Button { margin-right: 1; min-width: 22; }
    """

    def __init__(self, summary: "dict[str, int | bool]") -> None:
        super().__init__()
        self._summary = dict(summary or {})

    def compose(self) -> ComposeResult:
        s = self._summary
        files = int(s.get("files_removed", 0) or 0)
        dirs = int(s.get("dirs_removed", 0) or 0)
        logs = int(s.get("log_files_removed", 0) or 0)
        residual_f = int(s.get("residual_files", 0) or 0)
        residual_d = int(s.get("residual_dirs", 0) or 0)
        errors = int(s.get("errors", 0) or 0)
        pre = bool(s.get("pre_update_removed", False))
        residual_line = (
            f"  [dim]+ {residual_f} extra file(s) and "
            f"{residual_d} extra dir(s) caught by the "
            "residual-sweep pass.[/dim]\n"
            if (residual_f or residual_d) else ""
        )
        err_line = (
            f"\n  [yellow]⚠ {errors} path(s) could not be removed "
            "— see the log file for details.[/yellow]"
            if errors else ""
        )
        with Vertical(id="mdr-dlg"):
            yield Static(
                "  Master Delete complete  ",
                id="mdr-title", markup=False,
            )
            yield Static(
                f"  [bold]{files}[/bold] user-data file(s) removed.\n"
                f"  [bold]{dirs}[/bold] directory tree(s) removed.\n"
                f"  [bold]{logs}[/bold] rotated log backup(s) "
                "removed.\n"
                f"  Pre-update backups: "
                + ("[bold]wiped[/bold]" if pre else "[dim]none "
                   "present[/dim]") +
                ".\n"
                f"{residual_line}"
                f"{err_line}\n\n"
                "  Some panels may still show stale state in this\n"
                "  session — quit and restart SpliceCraft for a\n"
                "  fully clean canvas.",
                id="mdr-msg", markup=True,
            )
            with Horizontal(id="mdr-btns"):
                yield Button(
                    "OK", id="btn-mdr-ok", variant="primary",
                )

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-mdr-ok", Button).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-mdr-ok")
    def _ok(self, _) -> None:
        self.dismiss(None)

    def action_ok(self) -> None:
        self.dismiss(None)


class RenamePlasmidModal(_OneShotDismissScreen, ModalScreen):
    """Prompt for a new name for a library entry.

    Tab cycles between the Input and Save/Cancel buttons.
    Dismisses with the new name (a non-empty string) or None on cancel.
    Input validation (non-empty, trimmed, collision check) lives in the
    app-side handler — the modal just collects a value.
    """

    _blocks_undo: bool = True   # Input editing; result mutates library entry

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, current_name: str, entry_id: str):
        super().__init__()
        self.current_name = current_name
        self.entry_id     = entry_id

    def compose(self) -> ComposeResult:
        from rich.markup import escape as _md_escape
        with Vertical(id="rename-dlg"):
            yield Static(" Rename plasmid ", id="rename-title")
            # Escape — `current_name` is user-controlled and Label
            # interprets Rich markup by default, so a name like
            # "TU [draft]" would otherwise eat the trailing bracket
            # as a malformed tag. Same hygiene as NamePlasmidModal's
            # status line + the History viewer's parent-name render.
            yield Label(
                f"Current name:  {_md_escape(self.current_name)}",
            )
            yield Label("New name:")
            yield Input(
                value=self.current_name,
                placeholder="enter a new name",
                id="rename-input",
            )
            yield Static("", id="rename-status", markup=True)
            with Horizontal(id="rename-btns"):
                yield Button("Save",   id="btn-rename-save",   variant="primary")
                yield Button("Cancel", id="btn-rename-cancel")

    def on_mount(self) -> None:
        # Default focus on the Input, text pre-selected via select_on_focus
        # (Textual Input defaults to selecting all when focused, which is
        # what you want for a rename — typing replaces the old name).
        inp = self.query_one("#rename-input", Input)
        inp.focus()

    @on(Button.Pressed, "#btn-rename-save")
    def _save(self, _):
        self._try_submit()

    @on(Input.Submitted, "#rename-input")
    def _submitted(self, _):
        self._try_submit()

    def _try_submit(self) -> None:
        new_name = self.query_one("#rename-input", Input).value.strip()
        status   = self.query_one("#rename-status", Static)
        if not new_name:
            status.update("[red]Name cannot be empty.[/red]")
            return
        if new_name == self.current_name:
            # No-op rename — treat as cancel so the app doesn't bother writing.
            self.dismiss(None)
            return
        self.dismiss(new_name)

    @on(Button.Pressed, "#btn-rename-cancel")
    def _cancel_btn(self, _):
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class MinPrimerBindingModal(_OneShotDismissScreen, ModalScreen):
    """Prompt for the `min_primer_binding` threshold (integer bp).

    Replaces the old preset-cycle action with a free-typed integer
    input. Validates against `_settings_validator_int_range(1, 60)` —
    the same range the agent endpoint enforces — so the modal can't
    persist a value the loader would reject on next launch. Dismisses
    with the chosen integer (1–60) on submit, or `None` on cancel /
    no-change.
    """

    _blocks_undo: bool = True   # Input editing; result becomes settings value

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next",   show=False),
    ]

    DEFAULT_CSS = """
    #mpb-dlg {
        width: 56; height: auto;
        background: $surface; border: solid $primary; padding: 1 2;
    }
    #mpb-title  { background: $primary-darken-2; color: $text;
                  padding: 0 1; margin-bottom: 1; text-align: center; }
    #mpb-help   { height: auto; color: $text-muted; margin-bottom: 1; }
    #mpb-input  { margin-top: 0; margin-bottom: 1; }
    #mpb-status { height: 1; color: $text-muted; }
    #mpb-btns   { align: center middle;  height: 3; margin-top: 1; }
    #mpb-btns Button { margin-right: 1; }
    """

    _MIN_BP = 1
    _MAX_BP = 60

    def __init__(self, current_value: int) -> None:
        super().__init__()
        self.current_value = int(current_value)

    def compose(self) -> ComposeResult:
        with Vertical(id="mpb-dlg"):
            yield Static(" Minimum primer binding length ", id="mpb-title")
            yield Static(
                "Primers whose bound region is shorter than this "
                "threshold are flagged with a yellow ⚠ on the "
                "sequence panel. Allowed range: "
                f"{self._MIN_BP}–{self._MAX_BP} bp.",
                id="mpb-help",
            )
            yield Label("Threshold (bp):")
            yield Input(
                value=str(self.current_value),
                placeholder=f"{self._MIN_BP}–{self._MAX_BP}",
                id="mpb-input",
            )
            yield Static("", id="mpb-status", markup=True)
            with Horizontal(id="mpb-btns"):
                yield Button("Apply",  id="btn-mpb-apply",  variant="primary")
                yield Button("Cancel", id="btn-mpb-cancel")

    def on_mount(self) -> None:
        # Default focus on the input so the user can immediately type a
        # new value; the existing text is preselected (Textual default)
        # so a single keystroke replaces the displayed number.
        self.query_one("#mpb-input", Input).focus()

    @on(Button.Pressed, "#btn-mpb-apply")
    def _apply_btn(self, _) -> None:
        self._try_submit()

    @on(Input.Submitted, "#mpb-input")
    def _submitted(self, _) -> None:
        self._try_submit()

    def _try_submit(self) -> None:
        raw    = self.query_one("#mpb-input", Input).value.strip()
        status = self.query_one("#mpb-status", Static)
        if not raw:
            status.update("[red]Enter a value.[/red]")
            return
        try:
            value = int(raw)
        except ValueError:
            status.update("[red]Must be an integer.[/red]")
            return
        if not (self._MIN_BP <= value <= self._MAX_BP):
            status.update(
                f"[red]Out of range — pick "
                f"{self._MIN_BP}–{self._MAX_BP}.[/red]"
            )
            return
        if value == self.current_value:
            # No-op — treat as cancel so we don't bother re-stamping
            # primers + writing settings for an unchanged value.
            self.dismiss(None)
            return
        self.dismiss(value)

    @on(Button.Pressed, "#btn-mpb-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class PartsBinDeleteConfirmModal(_OneShotDismissScreen, ModalScreen):
    """Confirm-on-delete modal for the parts bin. Adapts its message,
    title, and primary-button label to the count of selected parts so
    the user always knows whether they're about to commit a single
    delete or a bulk delete.

    For single delete: shows the part name in bold.
    For multi delete: shows the count + a name preview (first 3, with
    "(+N more)" tail when over 3).

    Default focus on [No] so a stray Enter bails out — same handslip
    protection as `LibraryDeleteConfirmModal`. Dismisses True on
    confirm, False on cancel / Escape.
    """

    _blocks_undo: bool = True   # [INV-50] destructive confirm — Ctrl+Z above

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next button", show=False),
    ]

    DEFAULT_CSS = """
    #partsdel-dlg {
        width: 76; max-width: 95%; min-width: 56;
        height: auto; max-height: 90%;
        background: $surface; border: solid $error; padding: 1 2;
    }
    #partsdel-title { background: $error-darken-2; color: $text;
                      padding: 0 1; margin-bottom: 1; text-align: center; }
    #partsdel-msg   { height: auto; margin-bottom: 1; }
    #partsdel-btns  { align: center middle;  height: 3; margin-top: 1; }
    #partsdel-btns Button { margin-right: 1; min-width: 16; }
    """

    def __init__(self, names: list[str]) -> None:
        super().__init__()
        self._names = list(names)

    def compose(self) -> ComposeResult:
        from rich.markup import escape as _esc
        n = len(self._names)
        # Preview: first three names with markup escaped — a part
        # called `[red]boom[/]` shouldn't reformat the dialog.
        preview = ", ".join(_esc(name) for name in self._names[:3])
        if n > 3:
            preview += f" (+{n - 3} more)"
        if n == 1:
            title = " Remove part from bin "
            body = (
                f"  Remove [bold]{_esc(self._names[0])}[/bold] "
                f"from the parts bin?\n\n"
                f"  [dim]This cannot be undone from within the app. "
                f"A backup (.bak) of the parts-bin file is kept.[/dim]"
            )
            yes_label = "Yes, remove"
        else:
            # Bulk-delete: render the count in red bold so it's
            # impossible to miss before clicking Yes. The title +
            # button label echo the count too as belt + braces.
            title = f" Remove {n} parts from bin "
            body = (
                f"  Remove [bold red]{n}[/bold red] parts from "
                f"the parts bin?\n"
                f"  [dim]({preview})[/dim]\n\n"
                f"  [dim]This cannot be undone from within the app. "
                f"A backup (.bak) of the parts-bin file is kept.[/dim]"
            )
            yes_label = f"Yes, remove all {n}"
        with Vertical(id="partsdel-dlg"):
            yield Static(title, id="partsdel-title")
            yield Static(body, id="partsdel-msg", markup=True)
            with Horizontal(id="partsdel-btns"):
                yield Button("No",   id="btn-partsdel-no",
                              variant="default")
                yield Button(yes_label, id="btn-partsdel-yes",
                              variant="error")

    def on_mount(self) -> None:
        self.query_one("#btn-partsdel-no", Button).focus()

    @on(Button.Pressed, "#btn-partsdel-no")
    def _no(self, _) -> None:
        self.dismiss(False)

    @on(Button.Pressed, "#btn-partsdel-yes")
    def _yes(self, _) -> None:
        self.dismiss(True)

    def action_cancel(self) -> None:
        self.dismiss(False)


class ExactCopyConfirmModal(ModalScreen):
    """Prompt the user when a load batch contains exact duplicates of
    existing entries (same name AND same content). Two-way choice:

      * "Skip duplicates" — drop the duplicates from the load; existing
        entries unchanged. Default focus.
      * "Keep as COPY"    — append " COPY" (or " COPY 2", etc.) to each
        duplicate's name so it lands alongside the existing entry.

    Default focus on Skip so a stray Enter on splash dismiss doesn't
    silently double-add data. Generic across subsystems — used by parts
    bin, plasmid library, primers, experiments, collections.

    Dismiss payload:
      True  — Keep as COPY (caller renames + saves)
      False — Skip duplicates (caller drops them from the batch)
      None  — Escape (treated as skip, same as False)

    See ``_classify_collisions`` for the upstream pure classifier and
    ``_ensure_unique_copy_name`` for the rename helper.
    """

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #copydlg-dlg {
        width: 80; max-width: 95%; min-width: 60;
        height: auto; max-height: 90%;
        background: $surface; border: solid $warning; padding: 1 2;
    }
    #copydlg-title { background: $warning-darken-2; color: $text;
                     padding: 0 1; margin-bottom: 1; text-align: center; }
    #copydlg-msg   { height: auto; margin-bottom: 1; }
    #copydlg-btns  { height: 3; margin-top: 1; align: right middle; }
    #copydlg-btns Button { margin-left: 1; min-width: 18; }
    """

    def __init__(self, entity: str, names: "list[str]") -> None:
        """``entity`` — singular noun for the modal copy ("part",
        "plasmid", "primer", "experiment", "collection"). ``names`` —
        list of duplicate entry names; first three are previewed.
        """
        super().__init__()
        self._entity = entity
        self._names  = list(names)
        # Sweep #10 (2026-05-20): rapid double-click between
        # button-press and the async pop_screen completing can fire
        # the same handler twice with a stale `dismiss` call slipping
        # through. `_dismissed` gates every exit path so the modal
        # callback is invoked exactly once.
        self._dismissed: bool = False

    def compose(self) -> ComposeResult:
        from rich.markup import escape as _esc
        n = len(self._names)
        preview = ", ".join(_esc(name) for name in self._names[:3])
        if n > 3:
            preview += f" (+{n - 3} more)"
        ent = _esc(self._entity)
        plural = "" if n == 1 else "s"
        body = (
            f"  [bold red]{n}[/bold red] {ent}{plural} in this load "
            f"match {'an existing entry' if n == 1 else 'existing entries'} "
            f"by both name AND content (exact duplicate"
            f"{plural}).\n"
            f"  [dim]({preview})[/dim]\n\n"
            f"  • [b]Skip[/b] — drop the duplicate"
            f"{plural} from the load (default).\n"
            f"  • [b]Keep as COPY[/b] — append [b]“COPY”[/b] "
            f"to the name so {'it' if n == 1 else 'they'} can coexist "
            f"alongside the original{plural}.\n"
        )
        with Vertical(id="copydlg-dlg"):
            yield Static(
                f" Exact duplicate{plural} detected ", id="copydlg-title",
            )
            yield Static(body, id="copydlg-msg", markup=True)
            with Horizontal(id="copydlg-btns"):
                yield Button(
                    f"Skip duplicate{plural}",
                    id="btn-copydlg-skip", variant="default",
                )
                yield Button(
                    "Keep as COPY",
                    id="btn-copydlg-keep", variant="warning",
                )

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-copydlg-skip", Button).focus()
        except NoMatches:
            pass

    def _dismiss_once(self, payload) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    @on(Button.Pressed, "#btn-copydlg-skip")
    def _skip(self, _) -> None:
        _log_event(
            "collision.exact_copy.dismiss",
            entity=self._entity, n=len(self._names), choice="skip",
        )
        self._dismiss_once(False)

    @on(Button.Pressed, "#btn-copydlg-keep")
    def _keep(self, _) -> None:
        _log_event(
            "collision.exact_copy.dismiss",
            entity=self._entity, n=len(self._names), choice="keep",
        )
        self._dismiss_once(True)

    def action_cancel(self) -> None:
        _log_event(
            "collision.exact_copy.dismiss",
            entity=self._entity, n=len(self._names), choice="esc",
        )
        self._dismiss_once(False)


class NameCollisionModal(ModalScreen):
    """Prompt the user when a load batch contains entries whose name
    matches an existing entry but whose CONTENT differs. Three-way
    choice:

      * "Keep original"  — drop the new entries; existing entries
        unchanged. Default focus (safest — preserves user data).
      * "Overwrite"      — replace existing entries with the new content
        (matched by name).
      * "Cancel load"    — abort the entire load; nothing saved on
        either side.

    Used across subsystems (parts bin, plasmid library, primers,
    experiments, collections, gels) so a "load with a familiar name but
    different bytes" never silently clobbers or silently drops data.

    Dismiss payload (str):
      "keep"      — keep originals; new entries dropped
      "overwrite" — replace existing with new
      "cancel"    — abort the whole load (caller does NOTHING)
      None        — Escape (treated as cancel)

    See ``_classify_collisions`` for the upstream pure classifier.
    """

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #namecoll-dlg {
        width: 84; max-width: 95%; min-width: 64;
        height: auto; max-height: 90%;
        background: $surface; border: solid $warning; padding: 1 2;
    }
    #namecoll-title { background: $warning-darken-2; color: $text;
                      padding: 0 1; margin-bottom: 1; text-align: center; }
    #namecoll-msg   { height: auto; margin-bottom: 1; }
    #namecoll-btns  { height: 3; margin-top: 1; align: right middle; }
    #namecoll-btns Button { margin-left: 1; min-width: 18; }
    """

    def __init__(self, entity: str, names: "list[str]") -> None:
        """``entity`` — singular noun for the modal copy ("part",
        "plasmid", "primer", "experiment", "collection"). ``names`` —
        list of colliding entry names; first three are previewed.
        """
        super().__init__()
        self._entity = entity
        self._names  = list(names)
        # Sweep #10 (2026-05-20): double-fire guard — see
        # ExactCopyConfirmModal for rationale.
        self._dismissed: bool = False

    def compose(self) -> ComposeResult:
        from rich.markup import escape as _esc
        n = len(self._names)
        preview = ", ".join(_esc(name) for name in self._names[:3])
        if n > 3:
            preview += f" (+{n - 3} more)"
        ent = _esc(self._entity)
        plural = "" if n == 1 else "s"
        body = (
            f"  [bold red]{n}[/bold red] {ent}{plural} in this load "
            f"share a name with {'an existing entry' if n == 1 else 'existing entries'} "
            f"but the [b]content is different[/b].\n"
            f"  [dim]({preview})[/dim]\n\n"
            f"  • [b]Keep original[/b] — drop the new "
            f"{'entry' if n == 1 else 'entries'}; existing data unchanged "
            f"(default).\n"
            f"  • [b]Overwrite[/b] — replace the existing "
            f"{'entry' if n == 1 else 'entries'} with the new content.\n"
            f"  • [b]Cancel[/b] — abort the entire load; nothing saved.\n"
        )
        with Vertical(id="namecoll-dlg"):
            yield Static(
                f" Name collision{plural} (different content) ",
                id="namecoll-title",
            )
            yield Static(body, id="namecoll-msg", markup=True)
            with Horizontal(id="namecoll-btns"):
                yield Button(
                    "Keep original",
                    id="btn-namecoll-keep", variant="default",
                )
                yield Button(
                    "Overwrite",
                    id="btn-namecoll-overwrite", variant="warning",
                )
                yield Button(
                    "Cancel load",
                    id="btn-namecoll-cancel", variant="error",
                )

    def on_mount(self) -> None:
        try:
            self.query_one("#btn-namecoll-keep", Button).focus()
        except NoMatches:
            pass

    def _dismiss_once(self, payload: str) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    @on(Button.Pressed, "#btn-namecoll-keep")
    def _keep(self, _) -> None:
        _log_event(
            "collision.name_collision.dismiss",
            entity=self._entity, n=len(self._names), choice="keep",
        )
        self._dismiss_once("keep")

    @on(Button.Pressed, "#btn-namecoll-overwrite")
    def _overwrite(self, _) -> None:
        _log_event(
            "collision.name_collision.dismiss",
            entity=self._entity, n=len(self._names),
            choice="overwrite",
        )
        self._dismiss_once("overwrite")

    @on(Button.Pressed, "#btn-namecoll-cancel")
    def _cancel_btn(self, _) -> None:
        _log_event(
            "collision.name_collision.dismiss",
            entity=self._entity, n=len(self._names), choice="cancel",
        )
        self._dismiss_once("cancel")

    def action_cancel(self) -> None:
        _log_event(
            "collision.name_collision.dismiss",
            entity=self._entity, n=len(self._names), choice="esc",
        )
        self._dismiss_once("cancel")


class LibraryDeleteConfirmModal(ModalScreen):
    """Generic delete-confirmation modal. Used by the plasmid library,
    primer library, and any future list that needs handslip protection.

    Default focus is on [No]. Tab cycles between [No] and [Yes, remove].
    Escape dismisses as False (cancel).
    """

    _blocks_undo: bool = True   # [INV-50] destructive confirm — Ctrl+Z above

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next button", show=False),
    ]

    def __init__(self, name: str, size: int, entry_id: str,
                 dup_count: "int | None" = None,
                 bulk_names: "list[str] | None" = None,
                 bulk_last_copies: "int | None" = None):
        """`dup_count` = how many OTHER entries carry this entry's exact
        sequence (by `_entry_content_key` — `gb_ref` when dehydrated,
        else a hash of inline `gb_text`). None = not computed (primer
        library and other non-plasmid callers), keeping the original
        wording rather than inventing a reassurance or an alarm.

        Drives the "is this the last copy?" line ([INV-167]). The
        2026-07-25 incident: entries were deleted as presumed copies when
        they were distinct sequences that merely shared a name, so the
        modal now answers that question instead of leaving it to the name.

        `bulk_names` turns this into the BULK confirm: the dialog leads with
        how many plasmids are marked and names them (capped) instead of
        naming one. `bulk_last_copies` is how many of them are the only
        entry holding their sequence — the bulk analogue of the
        `dup_count == 0` alarm, and the number that actually decides whether
        this delete loses data ([INV-167]). Both None = the original
        single-entry dialog, unchanged.
        """
        super().__init__()
        self.entry_name = name
        self.entry_size = size
        self.entry_id   = entry_id
        self.dup_count  = dup_count
        self.bulk_names = [str(n) for n in (bulk_names or [])]
        self.bulk_last_copies = bulk_last_copies
        # [INV-50] dismiss-once guard
        self._dismissed: bool = False

    def _dismiss_once(self, payload: bool) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(payload)

    def compose(self) -> ComposeResult:
        # Escape the name: this Static renders with markup=True, so a
        # plasmid called `pUC19 [old]` would have the bracketed part
        # SILENTLY SWALLOWED as a style tag and the dialog would ask you
        # to confirm deleting "pUC19" — the wrong name, in the one
        # dialog where that matters most ([INV-167]). Bracketed text is
        # routine in NCBI-derived names (`[Candida] glabrata`).
        from rich.markup import escape as _esc
        size_str = f" ({self.entry_size:,} bp)" if self.entry_size > 0 else ""
        # Answer "am I about to lose this sequence?" up front — a name
        # alone never settles it ([INV-167]).
        if self.dup_count is None:
            dup_line = ""
        elif self.dup_count > 0:
            other = "entry" if self.dup_count == 1 else "entries"
            dup_line = (
                f"  [green]An identical sequence is kept by "
                f"{self.dup_count} other {other}[/green] — "
                f"removing this one does not lose it.\n\n"
            )
        else:
            dup_line = (
                "  [bold red]This is the ONLY entry holding this "
                "sequence.[/bold red]\n"
                "  [red]No other entry has the same bases — removing it "
                "loses it.[/red]\n\n"
            )
        if self.bulk_names:
            n = len(self.bulk_names)
            shown = [_esc(s) for s in self.bulk_names[:8]]
            listing = "\n".join(f"    • {s}" for s in shown)
            if n > len(shown):
                listing += f"\n    [dim]… and {n - len(shown)} more[/dim]"
            # Lead with the count — "how many am I about to lose" is the
            # question a bulk delete has to answer before anything else.
            head = (f"  Remove [bold]{n} marked plasmid"
                    f"{'s' if n != 1 else ''}[/bold] from the library?\n\n"
                    f"{listing}\n\n")
            if self.bulk_last_copies is None:
                dup_line = ""
            elif self.bulk_last_copies > 0:
                # Counted as SEQUENCES, not entries: two marked rows sharing
                # one sequence lose one sequence between them, and saying
                # "2 entries are the only holder" would be false.
                k = self.bulk_last_copies
                dup_line = (
                    f"  [bold red]{k} sequence{'s' if k != 1 else ''} "
                    f"{'are' if k != 1 else 'is'} held by no other "
                    f"entry.[/bold red]\n"
                    f"  [red]Removing these loses "
                    f"{'them' if k != 1 else 'it'} — nothing else has the "
                    f"same bases.[/red]\n\n")
            else:
                dup_line = (
                    "  [green]Every sequence here is also kept by an entry "
                    "you are not deleting[/green] — nothing is lost.\n\n")
        else:
            head = (f"  Remove [bold]{_esc(str(self.entry_name))}[/bold]"
                    f"{size_str} from the library?\n\n")
        with Vertical(id="libdel-dlg"):
            yield Static(" Remove from library ", id="libdel-title")
            yield Static(
                f"{head}"
                f"{dup_line}"
                f"  [dim]Undo with Ctrl+Z (or u) while the library is\n"
                f"  focused. A backup (.bak) of the library file is kept,\n"
                f"  and today's snapshot survives a restart.[/dim]",
                id="libdel-msg",
                markup=True,
            )
            with Horizontal(id="libdel-btns"):
                yield Button("No",           id="btn-libdel-no",  variant="default")
                yield Button("Yes, remove",  id="btn-libdel-yes", variant="error")

    def on_mount(self) -> None:
        self.query_one("#btn-libdel-no", Button).focus()

    @on(Button.Pressed, "#btn-libdel-no")
    def _no(self, _):
        self._dismiss_once(False)

    @on(Button.Pressed, "#btn-libdel-yes")
    def _yes(self, _):
        self._dismiss_once(True)

    def action_cancel(self) -> None:
        self._dismiss_once(False)


class FastaFilePickerModal(_OneShotDismissScreen, ModalScreen):
    """Modal file browser that returns the path to a selected FASTA file.

    Dismisses with ``str`` (absolute path) on Open, or ``None`` on Cancel /
    Escape. FASTA files are painted lime green in the tree; other files
    are white so the user can scan a mixed directory quickly. The tree
    starts in ``start_path`` when given (and readable), else ``$HOME``."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, start_path: "str | None" = None) -> None:
        super().__init__()
        start = Path(start_path).expanduser() if start_path else Path.home()
        try:
            if not start.is_dir():
                start = Path.home()
        except OSError:
            start = Path.home()
        self._start = str(start)
        self._selected: "str | None" = None

    def compose(self) -> ComposeResult:
        with Vertical(id="fasta-box"):
            yield Static(" Open FASTA File ", id="fasta-title")
            yield Static(
                f"[dim]{self._start}[/dim]", id="fasta-header", markup=True
            )
            yield _FastaAwareDirectoryTree(self._start, id="fasta-tree")
            yield Static(
                "[dim]FASTA files are highlighted in lime green. "
                "Click a file, then Open.[/dim]",
                id="fasta-hint", markup=True,
            )
            yield Static("", id="fasta-status", markup=True)
            with Horizontal(id="fasta-btns"):
                yield Button("Open", id="btn-fasta-open",
                             variant="primary", disabled=True)
                yield Button("Cancel", id="btn-fasta-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#fasta-tree", _FastaAwareDirectoryTree).focus()
        except NoMatches:
            pass

    @on(DirectoryTree.FileSelected)
    def _on_file_selected(self, event) -> None:
        self._selected = str(event.path)
        try:
            self.query_one("#fasta-header", Static).update(
                f"[dim]{self._selected}[/dim]"
            )
            self.query_one("#btn-fasta-open", Button).disabled = False
            self.query_one("#fasta-status", Static).update("")
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-fasta-open")
    def _open(self) -> None:
        if self._selected:
            self.dismiss(self._selected)
            return
        try:
            self.query_one("#fasta-status", Static).update(
                "[red]Pick a file first.[/red]"
            )
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-fasta-cancel")
    def _cancel_btn(self) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class GroupNamePromptModal(ModalScreen):
    """Sweep #29: simple name prompt for "Save group as library
    entry". Returns `{"name": <str>}` on Save, `None` on Cancel.

    Sweep #30 (2026-05-26): generalized — accepts custom title,
    prompt label, placeholder, and an `allow_empty` flag so the
    same modal serves both the original library-entry-naming
    use case AND the per-segment rename flow (where empty label
    is a legitimate "clear the label" intent)."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("enter",  "submit", "Save"),
    ]

    DEFAULT_CSS = """
    #gname-dlg {
        width: 70; height: auto;
        background: $surface; border: solid $primary;
        padding: 1 2;
    }
    #gname-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1;
        text-align: center; text-style: bold;
    }
    #gname-dlg Label { color: $text-muted; margin-top: 1; }
    #gname-input  { margin-top: 1; }
    #gname-status { height: 1; margin-top: 1; color: $text-muted; }
    #gname-btns   { height: 3; margin-top: 1;
                    align: right middle; }
    #gname-btns Button { margin-right: 1; min-width: 12; }
    """

    def __init__(
        self,
        default_name: str = "",
        *,
        title: str = " Save group as library entry ",
        prompt: str = "Name for the new feature-library group entry:",
        placeholder: str = "e.g. Esp3I → AATG adapter",
        allow_empty: bool = False,
    ) -> None:
        super().__init__()
        self._default     = str(default_name or "")
        self._title       = str(title)
        self._prompt      = str(prompt)
        self._placeholder = str(placeholder)
        self._allow_empty = bool(allow_empty)
        self._dismissed: bool = False

    def _dismiss_once(self, result) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(result)

    def compose(self) -> ComposeResult:
        with Vertical(id="gname-dlg"):
            yield Static(self._title, id="gname-title")
            yield Label(self._prompt)
            yield Input(value=self._default,
                         placeholder=self._placeholder,
                         id="gname-input")
            yield Static("", id="gname-status", markup=True)
            with Horizontal(id="gname-btns"):
                yield Button("Save", id="btn-gname-save",
                             variant="primary")
                yield Button("Cancel", id="btn-gname-cancel")

    def on_mount(self) -> None:
        try:
            inp = self.query_one("#gname-input", Input)
            inp.focus()
            # Select-all so the user can immediately type to
            # overwrite the prefilled default — matches the
            # SplitPositionPromptModal UX.
            try:
                inp.action_select_all()
            except (AttributeError, NoMatches):
                pass
        except NoMatches:
            pass

    def action_submit(self) -> None:
        self._save(None)

    @on(Button.Pressed, "#btn-gname-save")
    def _save(self, _) -> None:
        try:
            raw = self.query_one("#gname-input", Input).value
        except NoMatches:
            return
        # Sweep #29 hardening (2026-05-26): scrub control chars +
        # cap length so a pasted name with ANSI escapes / newlines
        # / null bytes can't corrupt the library JSON or downstream
        # Rich Text rendering. Same `_sanitize_label` used elsewhere
        # in the codebase for feature labels.
        name = _sanitize_label(raw, max_len=200)
        if not name and not self._allow_empty:
            try:
                self.query_one("#gname-status", Static).update(
                    "[red]Name cannot be empty (after stripping "
                    "control characters).[/red]"
                )
            except NoMatches:
                pass
            return
        self._dismiss_once({"name": name})

    @on(Button.Pressed, "#btn-gname-cancel")
    def _cancel(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)


class MigrateImportPickerModal(_OneShotDismissScreen, ModalScreen):
    """Browse for a SpliceCraft data ``.zip`` to import. ``.zip`` files
    are highlighted lime-green. Dismisses with the path or ``None``."""

    _blocks_undo: bool = True
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    DEFAULT_CSS = """
    #migimp-box { width: 96; max-width: 95%; height: auto; max-height: 92%;
        background: $surface; border: solid $primary; padding: 1 2; }
    #migimp-title { background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1; text-align: center; }
    #migimp-tree { height: 15; border: solid $primary-darken-2; }
    #migimp-status { height: 1; margin-top: 1; }
    #migimp-btns { height: 3; margin-top: 1; align: right middle; }
    #migimp-btns Button { margin-left: 1; }
    """

    def __init__(self, start_path: "str | None" = None) -> None:
        super().__init__()
        start = Path(start_path).expanduser() if start_path else Path.home()
        try:
            if not start.is_dir():
                start = Path.home()
        except OSError:
            start = Path.home()
        self._start = str(start)
        self._selected: "str | None" = None

    def compose(self) -> ComposeResult:
        with Vertical(id="migimp-box"):
            yield Static(" Import a data file (.zip) ", id="migimp-title")
            yield Static(f"[dim]{self._start}[/dim]", id="migimp-header",
                         markup=True)
            yield _ZipAwareDirectoryTree(self._start, id="migimp-tree")
            yield Static(
                "[dim].zip files are highlighted. Pick your exported data "
                "file, then Import.[/dim]", id="migimp-hint", markup=True)
            yield Static("", id="migimp-status", markup=True)
            with Horizontal(id="migimp-btns"):
                yield Button("Import", id="btn-migimp", variant="primary",
                             disabled=True)
                yield Button("Cancel", id="btn-migimp-cancel")

    @on(DirectoryTree.FileSelected)
    def _sel(self, event) -> None:
        self._selected = str(event.path)
        try:
            self.query_one("#migimp-header", Static).update(
                f"[dim]{self._selected}[/dim]")
            self.query_one("#btn-migimp", Button).disabled = False
            self.query_one("#migimp-status", Static).update("")
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-migimp")
    def _go(self, _) -> None:
        if self._selected:
            self.dismiss(self._selected)
            return
        try:
            self.query_one("#migimp-status", Static).update(
                "[red]Pick a .zip first.[/red]")
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-migimp-cancel")
    def _cxl(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class MultiAlignPickerModal(_OneShotDismissScreen, ModalScreen):
    """Modal: pick multiple library plasmids to align against the
    currently-loaded record. Each pick becomes one row on the linear-
    map alignment overlay.

    Layout: leading checkbox column + standard library columns. Space
    toggles the cursor row's selection; Align runs all picks. Filters
    the current plasmid out of the list (no self-self alignment).

    Dismiss payload:
      ``None``         — cancelled (Esc / Cancel)
      ``list[str]``    — selected entry IDs (may be empty if user
                          hit Align with nothing selected)
    """

    # Sweep #26: dispatches `_align_worker` on dismiss which mutates
    # per-target library mirrors. Block app-level Ctrl+Z while the
    # picker is open so the user can't accidentally roll back the
    # canvas between picking + running.
    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("space",  "toggle_selection",         "Toggle"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #mam-dlg {
        width: 96;
        /* Was rigid `height: 36;`; switched to flex height + cap so
           a short library list doesn't leave half a screen of dead
           space below the table (2026-05-20 UX audit). */
        height: 90%; max-height: 36; min-height: 18;
        background: $surface;
        border: solid $primary;
        padding: 1 2;
    }
    #mam-title {
        text-align: center;
        background: $primary-darken-2;
        color: $text;
        padding: 0 1;
    }
    #mam-help { color: $text-muted; padding: 0 1; margin-top: 1; }
    #mam-table { height: 1fr; margin-top: 1; }
    #mam-status { color: $text-muted; padding: 0 1; }
    /* `margin-top` adds breathing room ABOVE the container so the
       inner h=3 budget stays available for the buttons themselves.
       The old `padding-top: 1` ate one of those three rows, leaving
       only 2 for the h=3 button widget — its bottom border bled into
       the dialog's padding-bottom and got clipped (user report
       2026-05-23: "buttons cut off at the bottom"). */
    #mam-btns { height: 3; align: right middle; margin-top: 1; }
    """

    _CHECK_ON  = "[bold green]✓[/]"
    _CHECK_OFF = "[dim]·[/]"

    # Hard cap so a user with a 500-entry library can't kick off a
    # 500-target pairwise sweep in one click. The notify still
    # surfaces a count.
    _MAX_TARGETS = 20

    def __init__(self, current_id: "str | None" = None):
        super().__init__()
        self._current_id = current_id
        self._selected_ids: set[str] = set()

    def compose(self) -> ComposeResult:
        with Vertical(id="mam-dlg"):
            yield Static(" Align with library plasmids ", id="mam-title")
            yield Static(
                "Space toggles selection · Align runs pairwise "
                "alignments against the current plasmid · Esc cancels",
                id="mam-help", markup=False,
            )
            yield DataTable(id="mam-table", cursor_type="row",
                            zebra_stripes=True)
            yield Static("0 selected", id="mam-status")
            with Horizontal(id="mam-btns"):
                yield Button("Align", id="btn-mam-ok", variant="primary")
                yield Button("Cancel", id="btn-mam-cancel")

    def on_mount(self) -> None:
        t = self.query_one("#mam-table", DataTable)
        t.add_columns("", "Name", "ID", "Size", "Features")
        # Natural-sort by display name; current plasmid (if any)
        # filtered out so the picker never offers self-self alignment.
        entries = sorted(
            (
                e for e in _iter_library_readonly()
                if e.get("id") and e.get("id") != self._current_id
            ),
            key=lambda e: _natural_sort_key(
                e.get("name") or e.get("id") or ""
            ),
        )
        for e in entries:
            t.add_row(
                Text.from_markup(self._CHECK_OFF),
                Text(e.get("name", "?"), style="bold"),
                e.get("id", "?"),
                f"{e.get('size', 0):,} bp",
                f"{e.get('n_feats', 0)}",
                key=e.get("id"),
            )
        if entries:
            t.move_cursor(row=0)
            t.focus()
        else:
            try:
                self.query_one("#mam-status", Static).update(
                    "No other plasmids in the active collection."
                )
            except NoMatches:
                pass

    def action_toggle_selection(self) -> None:
        """Flip the cursor row's selection state. Updates the column-0
        marker and the running count line. Named distinctly from
        Textual's base ``DOMNode.action_toggle(attribute_name)`` (which
        expects an arg) so the override doesn't trip
        ``reportIncompatibleMethodOverride`` — the binding above routes
        to this method by name."""
        t = self.query_one("#mam-table", DataTable)
        key = _cursor_row_key(t)
        if not key:
            return
        if key in self._selected_ids:
            self._selected_ids.remove(key)
            mark = self._CHECK_OFF
        else:
            if len(self._selected_ids) >= self._MAX_TARGETS:
                self.app.notify(
                    f"At most {self._MAX_TARGETS} targets per batch — "
                    "deselect one first.",
                    severity="warning", timeout=4,
                )
                return
            self._selected_ids.add(key)
            mark = self._CHECK_ON
        # Update the marker cell at the cursor row.
        try:
            from textual.coordinate import Coordinate
            row_idx = t.cursor_row
            t.update_cell_at(
                Coordinate(row_idx, 0), Text.from_markup(mark),
            )
        except Exception:
            _log.exception("MultiAlignPickerModal: cell update failed")
        # Refresh the running count.
        try:
            self.query_one("#mam-status", Static).update(
                f"{len(self._selected_ids)} selected"
            )
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-mam-ok")
    def _ok(self, _):
        # Empty-selection guard: closing with [] would silently no-op
        # the action handler. Surface a notify and keep the modal open
        # so the user can fix their pick.
        if not self._selected_ids:
            try:
                self.app.notify(
                    "Select at least one plasmid (space toggles the "
                    "cursor row), then press Align.",
                    severity="warning", timeout=4,
                )
            except Exception:
                pass
            return
        # Stale-id filter: if a library mutation happened while the
        # modal was open (agent endpoint deletion, external file edit,
        # collection switch from a sibling pane), `_selected_ids` may
        # carry entry ids that no longer resolve. Drop them silently
        # rather than dismissing with ghost ids that downstream
        # `_action_open_align_picker` would just `continue` past with
        # a "not found" warning per target. Use the readonly-iter
        # helper to skip the per-call deepcopy.
        try:
            live_ids = {
                e.get("id") for e in _iter_library_readonly()
                if isinstance(e, dict) and e.get("id")
            }
        except Exception:
            _log.exception(
                "MultiAlignPickerModal: live-id filter failed; "
                "dismissing with the full selection",
            )
            live_ids = None
        if live_ids is not None:
            picked = [i for i in self._selected_ids if i in live_ids]
            n_dropped = len(self._selected_ids) - len(picked)
            if n_dropped:
                try:
                    self.app.notify(
                        f"{n_dropped} selected target(s) were removed "
                        "from the library while this picker was open — "
                        "aligning the remaining {n}.".format(n=len(picked))
                        if picked else
                        f"All {n_dropped} selected target(s) were "
                        "removed from the library while this picker "
                        "was open — nothing to align.",
                        severity="warning", timeout=5,
                    )
                except Exception:
                    pass
                if not picked:
                    return
            self.dismiss(picked)
            return
        self.dismiss(list(self._selected_ids))

    @on(Button.Pressed, "#btn-mam-cancel")
    def _cancel_btn(self, _):
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class PlasmidPickerModal(_OneShotDismissScreen, ModalScreen):
    """Scrollable plasmid-picker modal. Shows all entries from the library.
    Dismisses with the selected entry's id, or None on cancel.
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, current_id: "str | None" = None):
        super().__init__()
        self._current_id = current_id

    def compose(self) -> ComposeResult:
        with Vertical(id="pick-dlg"):
            yield Static(" Select plasmid from library ", id="pick-title")
            yield DataTable(id="pick-table", cursor_type="row",
                            zebra_stripes=True)
            with Horizontal(id="pick-btns"):
                yield Button("Select",  id="btn-pick-ok",     variant="primary")
                yield Button("Cancel",  id="btn-pick-cancel")

    def on_mount(self) -> None:
        t = self.query_one("#pick-table", DataTable)
        t.add_columns("Name", "ID", "Size", "Features")
        cursor = 0
        # Natural-sort by display name so `pBin2` lands before
        # `pBin10` rather than the lexicographic disk order.
        entries = sorted(
            _load_library(),
            key=lambda e: _natural_sort_key(
                e.get("name") or e.get("id") or ""
            ),
        )
        for i, e in enumerate(entries):
            t.add_row(
                Text(e.get("name", "?"), style="bold"),
                e.get("id", "?"),
                f"{e.get('size', 0):,} bp",
                f"{e.get('n_feats', 0)}",
                key=e.get("id"),
            )
            if self._current_id and e.get("id") == self._current_id:
                cursor = i
        if entries:
            t.move_cursor(row=cursor)
            t.focus()

    @on(Button.Pressed, "#btn-pick-ok")
    def _select(self, _):
        self.dismiss(_cursor_row_key(self.query_one("#pick-table", DataTable)))

    @on(DataTable.RowSelected, "#pick-table")
    def _row_selected(self, event):
        # Enter-key selection = same as clicking Select
        if event.row_key and event.row_key.value:
            self.dismiss(event.row_key.value)

    @on(Button.Pressed, "#btn-pick-cancel")
    def _cancel_btn(self, _):
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class PrimerCsvImportModal(_OneShotDismissScreen, ModalScreen):
    """File browser → returns the path of a primer-order CSV to import.
    Dismisses with the selected ``str`` path or ``None`` on cancel; the caller
    runs `_import_primers_from_csv` (which validates + reports skips)."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab", "app.focus_next", "Next", show=False),
    ]
    DEFAULT_CSS = """
    #pcsvimp-box { width: 86; height: auto; max-height: 90%;
                   background: $surface; border: solid $accent; padding: 1 2; }
    #pcsvimp-title { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; text-align: center; }
    #pcsvimp-header { color: $text-muted; }
    #pcsvimp-tree { height: 16; border: solid $primary-darken-2; margin: 1 0; }
    #pcsvimp-btns { height: 3; align: right middle; margin-top: 1; }
    #pcsvimp-btns Button { margin-left: 2; }
    """

    def __init__(self, start_path: "str | None" = None) -> None:
        super().__init__()
        start = Path(start_path).expanduser() if start_path else Path.home()
        try:
            if not start.is_dir():
                start = Path.home()
        except OSError:
            start = Path.home()
        self._start = str(start)
        self._selected: "str | None" = None

    def compose(self) -> ComposeResult:
        with Vertical(id="pcsvimp-box"):
            yield Static(" Import primers from CSV ", id="pcsvimp-title")
            yield Static(f"[dim]{self._start}[/dim]",
                         id="pcsvimp-header", markup=True)
            yield _ExtensionAwareDirectoryTree(
                self._start,
                highlight_map={".csv": _PICKER_PLASMID_STYLE},
                id="pcsvimp-tree")
            yield Static(
                "[dim]Pick a CSV (Name, Sequence[, Tm]), then Open. "
                "Invalid oligos are skipped + reported.[/dim]",
                id="pcsvimp-hint", markup=True)
            yield Static("", id="pcsvimp-status", markup=True)
            with Horizontal(id="pcsvimp-btns"):
                yield Button("Open", id="btn-pcsvimp-open",
                             variant="primary", disabled=True)
                yield Button("Cancel", id="btn-pcsvimp-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#pcsvimp-tree",
                           _ExtensionAwareDirectoryTree).focus()
        except NoMatches:
            pass

    @on(DirectoryTree.FileSelected, "#pcsvimp-tree")
    def _file_sel(self, event) -> None:
        self._selected = str(event.path)
        try:
            self.query_one("#pcsvimp-header", Static).update(
                f"[dim]{self._selected}[/dim]")
            self.query_one("#btn-pcsvimp-open", Button).disabled = False
            self.query_one("#pcsvimp-status", Static).update("")
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-pcsvimp-open")
    def _open(self) -> None:
        if self._selected:
            self.dismiss(self._selected)
            return
        try:
            self.query_one("#pcsvimp-status", Static).update(
                "[red]Pick a file first.[/red]")
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-pcsvimp-cancel")
    def _cancel_btn(self) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class PrimerPlasmidsModal(_OneShotDismissScreen, ModalScreen):
    """Surfaces every plasmid (across every collection) that carries a
    `primer_bind` feature matching a given primer-library entry's
    sequence. Lets the user pick one to jump to — the dismiss handler
    on `PrimerDesignScreen` then closes the primer-design screen and
    `PlasmidApp._goto_primer_in_plasmid` navigates to the chosen
    plasmid + scrolls the seq-panel cursor to the primer's binding
    region, mirroring the click-a-feature UX.

    Dismiss payload:
      ``None`` — cancelled (Escape / Close button)
      ``dict`` — chosen usage entry with keys
        ``collection``, ``plasmid_id``, ``plasmid_name``,
        ``start``, ``end``, ``strand``.
    """

    _blocks_undo: bool = True   # caller may swap `_current_record`

    BINDINGS = [
        Binding("escape",     "cancel",             "Cancel"),
        Binding("tab",        "app.focus_next",     "Next",  show=False),
        Binding("shift+tab",  "app.focus_previous", "Prev",  show=False),
    ]

    DEFAULT_CSS = """
    #pmp-dlg {
        width: 92;
        height: auto; max-height: 40;
        background: $surface;
        border: solid $primary;
        padding: 1 2;
    }
    #pmp-title {
        text-align: center;
        background: $primary-darken-2;
        color: $text;
        padding: 0 1;
        height: 1;
    }
    #pmp-info {
        margin: 1 0;
        padding: 0 1;
        color: $text;
        height: auto;
    }
    #pmp-table { height: 1fr; margin-top: 1; }
    /* `min-height: 3` + `height: auto` so the buttons get the height
       they declare even on terminals where the dialog box is the
       minimum size — previous fixed `height: 3` could clip on tiny
       terminals where the dialog box itself wasn't quite tall enough.
       The Vertical above (#pmp-dlg) is 40 rows fixed, so the table
       absorbs the slack. */
    #pmp-btns {
        height: auto;
        min-height: 3;
        align: right middle;
        padding-top: 1;
    }
    #pmp-btns Button { margin-left: 1; }
    """

    def __init__(self, primer_entry: dict, usages: list[dict]):
        super().__init__()
        self._primer = primer_entry
        # Sort once at construction (natural sort by collection, then
        # plasmid name) so cursor-row → usage mapping is stable.
        self._sorted_usages = sorted(
            usages,
            key=lambda u: (
                _natural_sort_key(u.get("collection") or ""),
                _natural_sort_key(u.get("plasmid_name") or ""),
            ),
        )

    def compose(self) -> ComposeResult:
        seq = (self._primer.get("sequence") or "")
        seq_preview = seq[:80] + ("…" if len(seq) > 80 else "")
        tm = self._primer.get("tm")
        tm_str = (f"{float(tm):.1f}°C"
                  if isinstance(tm, (int, float)) else "—")
        name = self._primer.get("name", "?")
        n_usages = len(self._sorted_usages)
        with Vertical(id="pmp-dlg"):
            yield Static(f" Plasmids using primer '{name}' ", id="pmp-title")
            yield Static(
                f"Sequence (5'→3'):  {seq_preview}\n"
                f"Length: {len(seq)} nt    Tm: {tm_str}    "
                f"Found in {n_usages} plasmid"
                f"{'s' if n_usages != 1 else ''}",
                id="pmp-info",
                markup=False,
            )
            yield DataTable(id="pmp-table", cursor_type="row",
                              zebra_stripes=True)
            with Horizontal(id="pmp-btns"):
                yield Button("Open plasmid", id="btn-pmp-open",
                              variant="primary")
                yield Button("Cancel", id="btn-pmp-cancel")

    def on_mount(self) -> None:
        t = self.query_one("#pmp-table", DataTable)
        t.add_columns("Collection", "Plasmid", "Position", "Strand")
        for u in self._sorted_usages:
            start = int(u.get("start") or 0)
            end = int(u.get("end") or 0)
            # 1-based inclusive for display (matches GenBank
            # convention); wrap-aware features (end < start) render
            # as "S..0..E" to match `_feat_span_label` style.
            if end < start:
                pos_str = f"{start + 1}..0..{end}"
            else:
                pos_str = f"{start + 1}–{end}"
            strand = int(u.get("strand") or 0)
            strand_str = ("+" if strand == 1
                            else "-" if strand == -1 else ".")
            t.add_row(
                u.get("collection") or "",
                u.get("plasmid_name") or "?",
                pos_str,
                strand_str,
            )
        t.focus()

    @on(DataTable.RowSelected, "#pmp-table")
    def _row_selected(self, event: DataTable.RowSelected) -> None:
        row = event.cursor_row
        if 0 <= row < len(self._sorted_usages):
            self.dismiss(self._sorted_usages[row])

    @on(Button.Pressed, "#btn-pmp-open")
    def _btn_open(self, _) -> None:
        try:
            t = self.query_one("#pmp-table", DataTable)
        except NoMatches:
            return
        row = t.cursor_row
        if 0 <= row < len(self._sorted_usages):
            self.dismiss(self._sorted_usages[row])

    @on(Button.Pressed, "#btn-pmp-cancel")
    def _btn_cancel(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class SynthesisLoadModal(_OneShotDismissScreen, ModalScreen):
    """Picker for loading a linear-topology library entry into the
    synthesis editor. Filters the active library to entries whose
    parsed topology is ``linear`` (circular plasmids don't belong in
    the synthesis editor — they have an origin invariant the linear
    editor can't represent).

    Dismiss payload:
      * ``str`` — the entry id to load.
      * ``None`` — user cancelled."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("enter",  "pick",   "Load"),
    ]

    DEFAULT_CSS = """
    #sl-dlg {
        width: 70%; max-width: 90; height: 80%; max-height: 38;
        background: $surface; padding: 1 2; border: solid $primary-darken-2;
    }
    #sl-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; height: 1; text-align: center;
    }
    #sl-search { height: 3; margin-top: 1; }
    #sl-table  { height: 1fr; border: solid $primary-darken-2;
                 margin-top: 1; }
    #sl-hint   { height: 1; color: $text-muted; margin-top: 1; }
    #sl-btns   { height: 3; align: right middle; margin-top: 1; }
    #sl-btns Button { margin-left: 1; min-width: 10; }
    """

    def __init__(self) -> None:
        super().__init__()
        self._rows: list[tuple[str, str, int]] = []  # (id, name, size)

    def compose(self) -> ComposeResult:
        with Vertical(id="sl-dlg"):
            yield Static(" Load fragment from library ",
                          id="sl-title")
            yield Input(placeholder="filter by name or id",
                          id="sl-search")
            yield DataTable(id="sl-table",
                              cursor_type="row",
                              zebra_stripes=True)
            yield Static(
                "[dim]Only linear-topology entries shown. "
                "Saving back overwrites the same entry (document model).[/]",
                id="sl-hint", markup=True,
            )
            with Horizontal(id="sl-btns"):
                yield Button("Load", id="btn-sl-load", variant="primary")
                yield Button("Cancel", id="btn-sl-cancel")

    def on_mount(self) -> None:
        try:
            t = self.query_one("#sl-table", DataTable)
            t.add_columns("Name", "ID", "bp")
        except NoMatches:
            return
        # Filter library to linear entries via a cheap LOCUS-line peek
        # (avoids parsing every gb_text just to read topology).
        # Sweep #26: readonly iter — loop body only `.get()`s fields
        # then appends to a fresh `rows` list.
        rows: list[tuple[str, str, int]] = []
        for e in _iter_library_readonly():
            if not isinstance(e, dict):
                continue
            gb_text = e.get("gb_text", "") or ""
            # LOCUS line in GenBank carries the topology word; cheap
            # substring test beats a full BioPython parse.
            first_line = gb_text.split("\n", 1)[0] if gb_text else ""
            if "linear" not in first_line.lower():
                continue
            rows.append((
                e.get("id", "") or "",
                e.get("name", "") or e.get("id", "") or "(unnamed)",
                int(e.get("size", 0) or 0),
            ))
        rows.sort(key=lambda r: _natural_sort_key(r[1]))
        self._rows = rows
        self._repopulate("")
        try:
            self.query_one("#sl-search", Input).focus()
        except NoMatches:
            pass

    def _repopulate(self, needle: str) -> None:
        try:
            t = self.query_one("#sl-table", DataTable)
        except NoMatches:
            return
        t.clear()
        needle_lo = needle.strip().lower()
        for eid, name, size in self._rows:
            if needle_lo and needle_lo not in name.lower() \
                    and needle_lo not in eid.lower():
                continue
            t.add_row(name, eid, f"{size:,}", key=eid)

    @on(Input.Changed, "#sl-search")
    def _on_search(self, event: Input.Changed) -> None:
        self._repopulate(event.value)

    @on(Input.Submitted, "#sl-search")
    def _on_search_submit(self, _) -> None:
        self.action_pick()

    @on(Button.Pressed, "#btn-sl-load")
    def _btn_load(self, _) -> None:
        self.action_pick()

    @on(Button.Pressed, "#btn-sl-cancel")
    def _btn_cancel(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)

    def action_pick(self) -> None:
        try:
            t = self.query_one("#sl-table", DataTable)
        except NoMatches:
            self.dismiss(None)
            return
        if t.cursor_row < 0 or t.row_count == 0:
            return
        try:
            row_key = t.coordinate_to_cell_key(
                _Coordinate(t.cursor_row, 0)
            ).row_key
        except (LookupError, AttributeError):
            self.dismiss(None)
            return
        eid = row_key.value if hasattr(row_key, "value") else str(row_key)
        self.dismiss(eid)


class CloneMethodChooserModal(_OneShotDismissScreen, ModalScreen):
    """Shown the moment the Synthesis screen's "Clone Fragment" button is
    pressed. The user picks HOW to clone the composed fragment — no save
    and no name prompt happen here (naming is deferred to the
    destination's own save step, where the user names the primed
    fragment + the cloned plasmid independently).

    Choices:
      * a modular-cloning GRAMMAR (Golden Braid L0, MoClo Plant, or any
        custom grammar from `_grammar_dropdown_options`) → routes to the
        Parts Domesticator, prefilled with the fragment, which designs
        domestication primers and produces the two L0 deliverables
        (primed linear fragment + cloned plasmid in the entry vector);
      * **Gibson assembly** → opens the Constructor's Gibson tab with the
        fragment pre-pasted into the lane's paste box;
      * **Traditional (restriction / ligation)** → opens the
        Constructor's Traditional tab, fragment pre-pasted.

    Dismisses with ``{"method": "grammar"|"gibson"|"traditional",
    "grammar_id": <id or "">}`` or ``None`` on cancel."""

    _blocks_undo: bool = True
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    DEFAULT_CSS = """
    #cmc-dlg {
        width: 76; height: auto; max-height: 90%;
        background: $surface; border: solid $primary; padding: 1 2;
    }
    #cmc-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1; text-align: center;
    }
    #cmc-hint { color: $text-muted; margin-bottom: 1; }
    #cmc-list { height: auto; max-height: 20; }
    #cmc-list Button { width: 100%; margin-bottom: 1; }
    #cmc-cancel-row { height: 3; margin-top: 1; align: right middle; }
    """

    def __init__(self, *, title: str = "", hint: str = "") -> None:
        super().__init__()
        # Snapshot the grammar list once so the button-id → grammar-id
        # map stays stable for the modal's lifetime (a grammar added in
        # another screen mid-modal can't shift the indices).
        self._grammars: "list[tuple[str, str]]" = _grammar_dropdown_options()
        # Optional caller wording so the same chooser serves both the Synthesis
        # "Clone Fragment" handoff (default) and the Alt+Shift+P selection→
        # pipeline hub ("Send selection to…").
        self._title: str = title or " Clone Fragment — choose a method "
        self._hint: str = hint or (
            "How should this fragment be cloned? Modular grammars "
            "domesticate it into an L0 part (primed fragment + cloned "
            "plasmid). Gibson / Traditional open the Constructor with "
            "the fragment pre-filled.")

    def compose(self) -> ComposeResult:
        with Vertical(id="cmc-dlg"):
            yield Static(self._title, id="cmc-title")
            yield Static(self._hint, id="cmc-hint", markup=False)
            with VerticalScroll(id="cmc-list"):
                for i, (label, _gid) in enumerate(self._grammars):
                    yield Button(label, id=f"cmc-g{i}", variant="primary",
                                 tooltip="Build an assembly with this grammar")
                yield Button("Gibson assembly", id="cmc-gibson")
                yield Button(
                    "Traditional (restriction / ligation)",
                    id="cmc-traditional",
                )
            with Horizontal(id="cmc-cancel-row"):
                yield Button("Cancel", id="cmc-cancel")

    @on(Button.Pressed)
    def _on_button(self, event: Button.Pressed) -> None:
        bid = event.button.id or ""
        if bid == "cmc-cancel":
            self.dismiss(None)
        elif bid == "cmc-gibson":
            self.dismiss({"method": "gibson", "grammar_id": ""})
        elif bid == "cmc-traditional":
            self.dismiss({"method": "traditional", "grammar_id": ""})
        elif bid.startswith("cmc-g"):
            try:
                idx = int(bid[len("cmc-g"):])
            except ValueError:
                return
            if 0 <= idx < len(self._grammars):
                self.dismiss({
                    "method": "grammar",
                    "grammar_id": self._grammars[idx][1],
                })

    def action_cancel(self) -> None:
        self.dismiss(None)


class SynthL0FragmentModal(_OneShotDismissScreen, ModalScreen):
    """Synthesis-tab "L0 Fragment": pick a grammar POSITION (Promoter / CDS /
    Terminator …) or enter CUSTOM overhangs, and the composed sequence is
    wrapped in the correct nested overhangs so the SYNTHESISED fragment becomes
    a proper Level-0 part once cloned into the configured UPD/pUPD entry vector
    with the grammar's Type IIS enzyme.

    Two-tier aware: when the entry vector's EXTERNAL acceptor overhangs (shown
    in the header) differ from the position's category overhangs, the category
    pair nests INSIDE the external pair — the same layout `_design_gb_primers`
    builds, so the fragment both enters the vector (external, e.g. CTCG/TGAG)
    and later releases the part on its category overhangs (e.g. AATG/GCTT).

    Inputs (passed by the caller so the modal stays data-access-free):
      * ``positions`` — ``[{name, type, oh5, oh3}, …]`` from the active grammar.
      * ``ext_overhangs`` — the vector's external ``(oh5, oh3)`` or ``None``.
      * ``ev_name`` / ``enzyme`` — for the header line.

    Dismisses ``{oh5, oh3, part_type, label}`` (category overhangs + the part
    type used for the ATG-fusion rule) or ``None`` on cancel."""

    _blocks_undo: bool = True
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    DEFAULT_CSS = """
    #l0-dlg {
        width: 84; height: auto; max-height: 90%;
        background: $surface; border: solid $primary; padding: 1 2;
    }
    #l0-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1; text-align: center;
    }
    #l0-info { color: $text-muted; margin-bottom: 1; }
    #l0-pos { height: auto; max-height: 16; margin-bottom: 1; }
    #l0-custom-row { height: 3; margin-bottom: 1; }
    #l0-custom-row Label { padding: 1 1 0 0; }
    #l0-custom-row Input { width: 14; margin-right: 2; }
    #l0-status { height: auto; margin-bottom: 1; }
    #l0-btns { height: 3; align: right middle; }
    #l0-btns Button { margin-left: 2; }
    """

    def __init__(self, positions, ext_overhangs, ev_name="", enzyme="") -> None:
        super().__init__()
        self._positions = list(positions or [])
        self._ext = tuple(ext_overhangs) if ext_overhangs else None
        self._ev_name = ev_name or ""
        self._enzyme = enzyme or "Esp3I"

    def compose(self) -> ComposeResult:
        with Vertical(id="l0-dlg"):
            yield Static(" L0 Fragment — pick the part position ", id="l0-title")
            if self._ext:
                info = (f"Entry vector [b]{self._ev_name or '—'}[/b]: a "
                        f"{self._enzyme} cut exposes [b]{self._ext[0]}[/b] / "
                        f"[b]{self._ext[1]}[/b] (external). The position's "
                        f"overhangs nest INSIDE these so it clones as a real L0 "
                        f"part.")
            else:
                info = ("No two-tier entry vector configured — the fragment "
                        "will carry the position's own overhangs directly "
                        "(one-tier). Bind a UPD vector in Settings to nest.")
            yield Static(info, id="l0-info")
            with RadioSet(id="l0-pos"):
                for i, p in enumerate(self._positions):
                    yield RadioButton(
                        f"{p.get('name','?')}  "
                        f"({p.get('oh5','')}→{p.get('oh3','')})  "
                        f"{p.get('type','')}",
                        value=(i == 0),
                    )
                yield RadioButton("Custom overhangs (enter below)")
            with Horizontal(id="l0-custom-row"):
                yield Label("Custom 5'")
                yield Input(placeholder="AATG", id="l0-oh5", max_length=8)
                yield Label("3'")
                yield Input(placeholder="GCTT", id="l0-oh3", max_length=8)
            yield Static("", id="l0-status", markup=True)
            with Horizontal(id="l0-btns"):
                yield Button("Make Fragment", id="l0-ok", variant="success")
                yield Button("Cancel", id="l0-cancel")

    @staticmethod
    def _clean_oh(val: str) -> str:
        return "".join(c for c in (val or "").upper() if c in "ACGT")

    @on(Button.Pressed, "#l0-cancel")
    def _cancel(self, _) -> None:
        self.dismiss(None)

    @on(Button.Pressed, "#l0-ok")
    def _ok(self, _) -> None:
        status = self.query_one("#l0-status", Static)
        try:
            idx = self.query_one("#l0-pos", RadioSet).pressed_index
        except NoMatches:
            idx = -1
        n = len(self._positions)
        if 0 <= idx < n:
            p = self._positions[idx]
            self.dismiss({"oh5": p.get("oh5", ""), "oh3": p.get("oh3", ""),
                          "part_type": p.get("type", ""),
                          "label": p.get("name", "")})
            return
        # Custom-overhang row (the entry AFTER every grammar position).
        oh5 = self._clean_oh(self.query_one("#l0-oh5", Input).value)
        oh3 = self._clean_oh(self.query_one("#l0-oh3", Input).value)
        if len(oh5) != 4 or len(oh3) != 4:
            status.update("[red]Custom overhangs must each be exactly 4 ACGT "
                          "bases (the Type IIS overhang width).[/red]")
            return
        self.dismiss({"oh5": oh5, "oh3": oh3, "part_type": "", "label": "Custom"})

    def action_cancel(self) -> None:
        self.dismiss(None)


class NamePlasmidModal(_OneShotDismissScreen, ModalScreen):
    """Prompt the user to name a freshly-assembled plasmid before
    it lands in the library.

    Default value is the auto-generated ``vector · part1+part2…``
    string from `ConstructorModal._compose_assembly_name`. The user
    can edit it freely; the dismiss flow re-runs
    ``_sanitize_plasmid_name`` so even a hand-pasted weird character
    can't reach the library.

    Shows a reference table of every plasmid already in the active
    collection so the user can pick a non-colliding name at a glance.
    Live duplicate-name check on `Input.Changed`: when the sanitised
    name matches an existing entry's name (case-insensitive) the Save
    button is disabled and the status line flags the collision.

    Dismiss payload:
      ``str``  — the sanitised name (always non-empty, non-duplicate).
      ``None`` — user cancelled; caller should NOT save.
    """

    _blocks_undo: bool = True   # Input editing; result becomes library entry name

    DEFAULT_CSS = """
    /* Tidy spacing: a min width so neither the inputs nor the optional
       "Primer Family Name" box clip, and one blank line between each
       stacked element so the labels + textboxes read cleanly. */
    #nameplasmid-dlg { min-width: 62; }
    #nameplasmid-input, #nameplasmid-collection,
    #nameplasmid-primer-family, #nameplasmid-status {
        margin-bottom: 1;
    }
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, default_name: str,
                 *, target_label: str = "plasmid",
                 default_collection: "str | None" = None,
                 primer_family: "str | None" = None,
                 exempt_id: "str | None" = None) -> None:
        super().__init__()
        # Opt-in: the id of an EXISTING library entry that this save will edit
        # in place (Alt+K re-saving the loaded plasmid). Its own name/id must
        # NOT count as a duplicate, else the modal blocks the user from
        # re-saving their own plasmid. Callers naming a FRESH artifact
        # (Constructor / Gibson / PCR / protein) pass nothing → no exemption.
        self._exempt_id = (exempt_id or "").strip()
        self._default_name = _sanitize_plasmid_name(
            default_name or "", fallback="assembly",
        )
        # ``target_label`` ("TU", "MOD", or just "plasmid") shows up
        # in the title so the user knows which level they're naming.
        self._target_label = target_label or "plasmid"
        # Opt-in "Primer Family Name" box: when given, the modal shows a
        # single labelled textbox and the dismiss dict carries
        # ``"primer_family"``. The operon save uses it to name the SOE
        # primer set ``{family}-DOM-{#}-{F|R}``. None → no box (every other
        # caller unaffected).
        self._primer_family_mode = primer_family is not None
        self._default_primer_family = (
            _sanitize_plasmid_name(primer_family or "", fallback="primers")
            if primer_family is not None else ""
        )
        # Opt-in collection picker: when ``default_collection`` is given
        # the modal shows a "Save to collection" Select and dismisses
        # with ``{"name": str, "collection": str}``; when it's None the
        # modal keeps its legacy contract (dismiss with a bare ``str``)
        # so existing callers are unaffected. The universal save flow
        # (Constructor / Gibson / Traditional / amplicon) passes the
        # active collection here so the user can redirect the save.
        self._collection_mode = default_collection is not None
        self._collections: list[str] = []
        self._default_collection = ""
        if self._collection_mode:
            names: list[str] = []
            try:
                for c in _iter_collections_readonly():
                    if not isinstance(c, dict):
                        continue
                    nm = (c.get("name") or "").strip()
                    if nm and nm not in names:
                        names.append(nm)
            except Exception:
                _log.exception("NamePlasmidModal: collection enum failed")
            active = (default_collection or "").strip()
            # Active collection first + guaranteed present so the Select
            # (allow_blank=False) always mounts with a valid value.
            ordered = ([active] if active else []) + [
                n for n in names if n != active
            ]
            if not ordered:
                ordered = [active or "Default"]
            self._collections = ordered
            self._default_collection = ordered[0]
        # Snapshot existing names / ids at modal-construct time —
        # the library doesn't mutate while the modal is open so a
        # one-shot read is enough and keeps the dup-check O(1) on
        # every keystroke. Names are case-folded so a user typing
        # "Demo 26" still flags the existing "DEMO 26".
        self._existing_names: dict[str, str] = {}
        self._existing_ids:   dict[str, str] = {}
        seen_name_keys: set[str] = set()
        dup_name_log: list[str] = []
        # Sweep #26: readonly walk — building maps is read-only over the
        # cache view, no deepcopy needed.
        for e in _iter_library_readonly():
            if not isinstance(e, dict):
                continue
            nm = (e.get("name") or "").strip()
            eid = (e.get("id") or "").strip()
            # The record being re-saved isn't a collision with ITSELF: drop its
            # own entry from the dup-scope so re-saving under the same name
            # reads as the intended edit-in-place (add_entry keys that on a
            # matching id), not a blocked "DUPLICATE". No exempt_id → no skip.
            if self._exempt_id and eid == self._exempt_id:
                continue
            if nm:
                key = nm.casefold()
                if key in seen_name_keys:
                    dup_name_log.append(nm)
                seen_name_keys.add(key)
                self._existing_names[key] = nm
            if eid:
                # Map id → DISPLAY name (fall back to id when name is
                # empty) so the dup-warning surfaces what the user sees
                # in the library row, not the bare stale id. After a
                # rename `e["id"]` is the OLD sanitised name (immutable
                # by design) while `e["name"]` is the user's new label
                # — showing the id confused users into thinking the
                # warning referenced a phantom old plasmid.
                self._existing_ids[eid.casefold()] = nm or eid
        # Surface data-integrity oddities: two library entries that
        # case-fold to the same display name are likely the
        # downstream of an unintended duplicate-save. The modal
        # itself prevents new duplicates; this log lets a user
        # diagnose existing ones via the diagnostic bundle.
        if dup_name_log:
            _log.warning(
                "NamePlasmidModal: library contains %d case-fold "
                "duplicate name(s): %s — first/last wins for dup-check",
                len(dup_name_log), dup_name_log[:5],
            )

    def compose(self) -> ComposeResult:
        with Vertical(id="nameplasmid-dlg"):
            yield Static(
                f" Name your {self._target_label} ",
                id="nameplasmid-title",
            )
            yield Label(
                "This name lands on the SeqRecord, the library row, "
                "and the Parts Bin entry. You can rename later.",
            )
            yield Input(
                value=self._default_name,
                placeholder="enter a name (default shown)",
                id="nameplasmid-input",
            )
            yield Static("", id="nameplasmid-status", markup=True)
            if self._collection_mode:
                yield Label("Save to collection:")
                yield Select(
                    [(c, c) for c in self._collections],
                    value=self._default_collection,
                    allow_blank=False,
                    id="nameplasmid-collection",
                )
            if self._primer_family_mode:
                yield Label("Primer Family Name:")
                yield Input(
                    value=self._default_primer_family,
                    placeholder="primer family (e.g. VhLux)",
                    id="nameplasmid-primer-family",
                )
            # Reference table of existing plasmids in the active
            # collection — read-only, sorted alphabetically (natural-
            # sort) so the user can scan for collisions at a glance.
            active_coll = _get_active_collection_name() or "library"
            yield Label(
                f"Existing plasmids in '{active_coll}' "
                f"({len(self._existing_names)}):",
                id="nameplasmid-list-label",
            )
            yield DataTable(
                id="nameplasmid-list",
                cursor_type="row",
                zebra_stripes=True,
                show_header=False,
            )
            with Horizontal(id="nameplasmid-btns"):
                yield Button("Save",   id="btn-nameplasmid-save",
                             variant="primary")
                yield Button("Cancel", id="btn-nameplasmid-cancel")

    def on_mount(self) -> None:
        # Populate the reference table — natural-sort by name so
        # `DEMO 2` lands before `DEMO 10`. Use the underlying display
        # names (not the case-folded keys) for visual fidelity.
        try:
            t = self.query_one("#nameplasmid-list", DataTable)
        except NoMatches:
            t = None
        if t is not None:
            t.add_columns("Name")
            display_names = sorted(
                set(self._existing_names.values()),
                key=_natural_sort_key,
            )
            if display_names:
                for nm in display_names:
                    t.add_row(
                        Text(nm, no_wrap=True, overflow="ellipsis"),
                    )
            else:
                # Empty-collection placeholder — better than a bare
                # empty DataTable which reads as "loading" or "broken".
                t.add_row(Text(
                    "(no plasmids yet in this collection)",
                    style="dim italic",
                ))
        try:
            inp = self.query_one("#nameplasmid-input", Input)
        except NoMatches:
            return
        inp.focus()
        # Run the dup check once on the default value so the user
        # sees the warning + the Save button reflects state without
        # having to type a character first.
        self._refresh_dup_state(inp.value)

    @on(Input.Changed, "#nameplasmid-input")
    def _on_input_changed(self, event: Input.Changed) -> None:
        # Live dup-check on every keystroke. `event.value` is the
        # post-keystroke string. Status-line update covers exact
        # match (red) AND substring near-miss (yellow) — see
        # `_refresh_dup_state` for the policy. The reference table
        # below the Input stays unfiltered so the user can scan all
        # existing names regardless of what they've typed.
        _log.debug(
            "NamePlasmidModal: Input.Changed value=%r len(value)=%d",
            event.value, len(event.value),
        )
        self._refresh_dup_state(event.value)

    def _refresh_dup_state(self, raw_value: str) -> "str | None":
        """Run the sanitise + duplicate check against the current
        Input value. Updates the status line + Save button. Returns
        the cleaned name when valid + non-duplicate, else ``None``.
        Centralised so `on_mount`, `Input.Changed`, and `_try_submit`
        all see consistent state without re-implementing the logic.

        Three severity levels:
          * **Exact dup** (case-folded name OR sanitised id matches an
            existing entry) → red status, Save disabled.
          * **Soft warning** (typed string is a substring of an
            existing name, OR an existing name is a substring of the
            typed string) → yellow status, Save enabled. User has to
            confirm by pressing Save anyway — soft warnings catch
            near-misses without blocking legitimate distinct names
            that happen to share a prefix.
          * **Available** → green status, Save enabled.
        """
        try:
            status = self.query_one("#nameplasmid-status", Static)
            save_btn = self.query_one(
                "#btn-nameplasmid-save", Button,
            )
        except NoMatches:
            return None
        cleaned = _sanitize_plasmid_name(
            raw_value, fallback=self._default_name,
        )
        if not cleaned:
            status.update("[bold red]Name cannot be empty.[/bold red]")
            save_btn.disabled = True
            return None
        # Markup-escape every user-controlled string that interpolates
        # into a `markup=True` Static. Without this, a saved entry
        # named "TU [draft]" would render the trailing "[draft]" as a
        # Rich-markup tag — visually broken at best, malformed-string
        # exception at worst. Sacred hygiene rule the History viewer
        # already follows (CLAUDE.md invariant #11). `cleaned` itself
        # comes from `_sanitize_plasmid_name` which strips paths and
        # control chars BUT preserves `[`, `]`, `<`, `>`, `&` so the
        # display side has to escape.
        from rich.markup import escape as _md_escape
        cleaned_safe = _md_escape(cleaned)
        # Detect leading / trailing whitespace separately from other
        # cleaning (illegal chars, length cap). Trailing whitespace
        # is the silent-cascade-breaker case (`'DEMO 27 ….dna'` → name
        # with trailing space → delete-cascade `==` miss against the
        # parts_bin row). Surface it as a distinct warning chip so the
        # user notices BEFORE saving, even though save will strip it
        # for them anyway.
        has_lead_ws = raw_value != raw_value.lstrip()
        has_trail_ws = raw_value != raw_value.rstrip()
        whitespace_warning = ""
        if has_lead_ws or has_trail_ws:
            sides = []
            if has_lead_ws:
                sides.append("leading")
            if has_trail_ws:
                sides.append("trailing")
            whitespace_warning = (
                f"[yellow]⚠ {' + '.join(sides)} whitespace will be "
                f"stripped[/yellow] — saved as [b]{cleaned_safe}[/b]"
            )
        # Cleaning-only mismatch (whitespace / illegal chars stripped).
        # Not an error — surface as info so the user can confirm.
        # When the only diff is leading/trailing whitespace, the
        # standalone `whitespace_warning` above already covers it;
        # otherwise (illegal chars, length cap) the generic hint
        # appended after status messages tells the user the final form.
        cleaning_hint = ""
        if cleaned != raw_value.strip():
            cleaning_hint = (
                f" [yellow](will save as[/yellow] [b]{cleaned_safe}[/b]"
                f"[yellow])[/yellow]"
            )
        # Case-fold for matching — typing "Demo 26" still flags the
        # existing "DEMO 26". Check both name and sanitised id space
        # since the library disambiguates by id.
        cleaned_cf = cleaned.casefold()
        cleaned_id = re.sub(r"[^A-Za-z0-9_]+", "_", cleaned)
        cleaned_id_cf = cleaned_id.casefold()
        if cleaned_cf in self._existing_names:
            actual = _md_escape(self._existing_names[cleaned_cf])
            status.update(
                f"[bold red]✗ DUPLICATE — already in use:[/bold red] "
                f"[b]{actual}[/b]{cleaning_hint}"
            )
            save_btn.disabled = True
            return None
        if cleaned_id_cf in self._existing_ids:
            actual = _md_escape(self._existing_ids[cleaned_id_cf])
            sanitised_safe = _md_escape(cleaned_id)
            status.update(
                f"[bold red]✗ Sanitised id[/bold red] [b]{sanitised_safe}[/b] "
                f"[bold red]clashes with[/bold red] [b]{actual}[/b] "
                f"(would auto-rename on save){cleaning_hint}"
            )
            save_btn.disabled = True
            return None
        # Soft warning: substring overlap with an existing name.
        # Either direction matters — a user typing "DEMO 32" sees the
        # existing "DEMO 32 VARA" as a near-match, and a user typing
        # "DEMO 32 NEW CDS" sees the existing "DEMO 32" as a near-match.
        # Save stays enabled — the names are distinct so the user can
        # legitimately commit. Just keeps near-misses visible.
        soft_hits: list[str] = []
        for existing_name in self._existing_names.values():
            ecf = existing_name.casefold()
            if (cleaned_cf in ecf or ecf in cleaned_cf) \
                    and cleaned_cf != ecf:
                soft_hits.append(existing_name)
                if len(soft_hits) >= 3:
                    break
        if soft_hits:
            preview = ", ".join(
                f"[b]{_md_escape(n)}[/b]" for n in soft_hits[:3]
            )
            status.update(
                f"[yellow]⚠ similar to:[/yellow] {preview}"
                f"{cleaning_hint}"
            )
            save_btn.disabled = False
            return cleaned
        if whitespace_warning:
            # Leading/trailing-space-only case gets the dedicated
            # warning chip (more prominent than the generic
            # "Will save as" hint).
            status.update(whitespace_warning)
        elif cleaning_hint:
            status.update(
                f"[yellow]Will save as[/yellow] [b]{cleaned_safe}[/b]"
            )
        else:
            status.update("[bold green]✓ Name available.[/bold green]")
        save_btn.disabled = False
        return cleaned

    @on(Button.Pressed, "#btn-nameplasmid-save")
    def _save(self, _) -> None:
        self._try_submit()

    @on(Input.Submitted, "#nameplasmid-input")
    def _submitted(self, _) -> None:
        self._try_submit()

    def _try_submit(self) -> None:
        try:
            inp = self.query_one("#nameplasmid-input", Input)
        except NoMatches:
            return
        # Final guard — `_refresh_dup_state` already disables Save on
        # dup / empty cases, but pressing Enter on the Input bypasses
        # the button so we re-validate here. Same dup-check the
        # button uses; returns None when not OK, in which case the
        # status line already shows why.
        cleaned = self._refresh_dup_state(inp.value)
        if cleaned is None:
            return
        # If sanitisation reshaped the user's input, reflect that
        # back into the field before dismiss so the next caller (rare
        # but possible — e.g. an agent inspecting the modal's last
        # value) sees the canonical form.
        try:
            if cleaned != inp.value:
                inp.value = cleaned
        except Exception:
            pass
        if self._collection_mode or self._primer_family_mode:
            payload: dict = {"name": cleaned}
            if self._collection_mode:
                payload["collection"] = self._selected_collection()
            if self._primer_family_mode:
                payload["primer_family"] = self._selected_primer_family()
            self.dismiss(payload)
        else:
            self.dismiss(cleaned)

    def _selected_primer_family(self) -> str:
        """Currently-entered primer family (primer-family mode only),
        sanitised; falls back to the default when blank / missing."""
        try:
            v = self.query_one("#nameplasmid-primer-family", Input).value
        except NoMatches:
            return self._default_primer_family or "primers"
        v = _sanitize_plasmid_name(v or "", fallback="primers")
        return v or self._default_primer_family or "primers"

    def _selected_collection(self) -> str:
        """Currently-picked collection (collection mode only). Falls back
        to the default if the Select is missing / blank."""
        try:
            v = self.query_one("#nameplasmid-collection", Select).value
        except NoMatches:
            return self._default_collection
        return v if isinstance(v, str) and v.strip() else \
            self._default_collection

    @on(Button.Pressed, "#btn-nameplasmid-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class RestrictionInsertModal(_OneShotDismissScreen, ModalScreen):
    """Searchable picker for inserting a restriction-enzyme recognition
    site at the synthesis cursor.

    Dismiss payload:
      * ``dict`` — ``{"enzyme": <name>, "strand": 1 | -1}``. The caller
        looks up the recognition site via ``_site_for_enzyme``; for
        ``strand == -1`` it inserts the reverse-complement (so a
        directional Type IIS enzyme cuts the other way) and draws the
        feature arrow ◀. ``strand == 1`` inserts the site as shown, arrow ▶.
      * ``None`` — user cancelled."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("enter",  "pick",   "Insert"),
    ]

    DEFAULT_CSS = """
    #ri-dlg {
        width: 60%; max-width: 80; height: 80%; max-height: 40;
        background: $surface; padding: 1 2; border: solid $primary-darken-2;
    }
    #ri-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; height: 1; text-align: center;
    }
    #ri-search { height: 3; margin-top: 1; }
    #ri-table  { height: 1fr; border: solid $primary-darken-2;
                 margin-top: 1; }
    #ri-dir-row { height: 3; margin-top: 1; }
    #ri-dir-label { width: auto; padding: 1 1 0 0; color: $text-muted; }
    #ri-direction { width: 1fr; }
    #ri-btns   { height: 3; align: right middle; margin-top: 1; }
    #ri-btns Button { margin-left: 1; min-width: 10; }
    """

    # Direction options for the Select — value is the feature strand.
    _DIR_OPTIONS = [
        ("Forward  ▶  (site as shown, cuts →)", "1"),
        ("Reverse  ◀  (reverse-complement, cuts ←)", "-1"),
    ]

    def __init__(self) -> None:
        super().__init__()
        # Pre-build the full enzyme list once so search re-filtering
        # is cheap. (name, site, n_bp).
        self._all_rows: list[tuple[str, str, int]] = sorted(
            (
                (name, info[0], len(info[0].replace("N", "")))
                for name, info in _state._all_enzymes_hook().items()
                if isinstance(info, tuple) and info
            ),
            key=lambda r: _natural_sort_key(r[0]),
        )

    def compose(self) -> ComposeResult:
        with Vertical(id="ri-dlg"):
            yield Static(" Insert restriction site at cursor ",
                          id="ri-title")
            yield Input(placeholder="filter by name or site "
                         "(e.g. EcoRI, GAATTC, BsaI)",
                          id="ri-search")
            yield DataTable(id="ri-table",
                              cursor_type="row",
                              zebra_stripes=True)
            with Horizontal(id="ri-dir-row"):
                yield Static("Direction:", id="ri-dir-label")
                yield Select(self._DIR_OPTIONS, value="1",
                              id="ri-direction", allow_blank=False)
            with Horizontal(id="ri-btns"):
                yield Button("Insert", id="btn-ri-insert", variant="primary")
                yield Button("Cancel", id="btn-ri-cancel")

    def on_mount(self) -> None:
        try:
            t = self.query_one("#ri-table", DataTable)
            t.add_columns("Enzyme", "Recognition", "bp")
        except NoMatches:
            return
        self._repopulate("")
        try:
            self.query_one("#ri-search", Input).focus()
        except NoMatches:
            pass

    def _repopulate(self, needle: str) -> None:
        try:
            t = self.query_one("#ri-table", DataTable)
        except NoMatches:
            return
        t.clear()
        needle_up = needle.strip().upper()
        for name, site, n_bp in self._all_rows:
            if needle_up and needle_up not in name.upper() \
                    and needle_up not in site.upper():
                continue
            t.add_row(name, site, str(n_bp), key=name)

    @on(Input.Changed, "#ri-search")
    def _on_search(self, event: Input.Changed) -> None:
        self._repopulate(event.value)

    @on(Input.Submitted, "#ri-search")
    def _on_search_submit(self, _) -> None:
        self.action_pick()

    @on(Button.Pressed, "#btn-ri-insert")
    def _btn_insert(self, _) -> None:
        self.action_pick()

    @on(Button.Pressed, "#btn-ri-cancel")
    def _btn_cancel(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)

    def action_pick(self) -> None:
        try:
            t = self.query_one("#ri-table", DataTable)
        except NoMatches:
            self.dismiss(None)
            return
        if t.cursor_row < 0 or t.row_count == 0:
            return
        try:
            row_key = t.coordinate_to_cell_key(
                _Coordinate(t.cursor_row, 0)
            ).row_key
        except (LookupError, AttributeError):
            self.dismiss(None)
            return
        enzyme = row_key.value if hasattr(row_key, "value") else str(row_key)
        # Read the direction Select — value is the feature strand as a
        # string ("1" / "-1"). Defensive: fall back to forward if the
        # Select is missing or holds an unexpected value.
        strand = 1
        try:
            raw = self.query_one("#ri-direction", Select).value
            strand = -1 if str(raw) == "-1" else 1
        except NoMatches:
            pass
        self.dismiss({"enzyme": enzyme, "strand": strand})


class CollectionNameModal(_OneShotDismissScreen, ModalScreen):
    """Tiny prompt modal for creating or renaming a collection.

    Dismisses with the trimmed name string, or None on cancel.
    Caller is responsible for collision-checking before persisting.
    """

    _blocks_undo: bool = True   # Input editing; result mutates collections

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    def __init__(self, title: str, current: str = "",
                 placeholder: str = "Collection name") -> None:
        super().__init__()
        self.title_text = title
        self.current = current
        self.placeholder_text = placeholder

    def compose(self) -> ComposeResult:
        with Vertical(id="collname-dlg"):
            yield Static(f" {self.title_text} ", id="collname-title")
            yield Label("Name:")
            yield Input(value=self.current,
                        placeholder=self.placeholder_text,
                        id="collname-input")
            yield Static("", id="collname-status", markup=True)
            with Horizontal(id="collname-btns"):
                yield Button("OK",     id="btn-collname-ok",     variant="primary")
                yield Button("Cancel", id="btn-collname-cancel")

    def on_mount(self) -> None:
        self.query_one("#collname-input", Input).focus()

    @on(Button.Pressed, "#btn-collname-ok")
    def _ok(self, _) -> None:
        self._submit()

    @on(Input.Submitted, "#collname-input")
    def _submitted(self, _) -> None:
        self._submit()

    def _submit(self) -> None:
        raw = self.query_one("#collname-input", Input).value
        # Same normaliser the agent-API uses: strip control chars,
        # trim, cap length. Stops a hand-typed `name\x00\n` from
        # breaking the panel header / collection-list rendering.
        name = _normalize_collection_name(raw)
        if name is None:
            self.query_one("#collname-status", Static).update(
                "[red]Name cannot be empty.[/red]"
            )
            return
        self.dismiss(name)

    @on(Button.Pressed, "#btn-collname-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class MultiRecordFastaModal(_OneShotDismissScreen, ModalScreen):
    """Prompt to import a multi-record FASTA into a new collection.

    Triggered by `OpenFileModal._do_load` when a picked FASTA contains
    more than one record. Shows the record count + a name input + a
    brief description of what will happen (each record becomes a
    separate library entry; topology is auto-detected from each
    record's description and otherwise defaults to linear).

    Pure prompt — doesn't touch disk. Dismisses with the validated
    collection name on submit, or `None` on cancel. The caller
    (`OpenFileModal`) builds the SeqRecords and calls
    `_create_fasta_collection` once the modal returns a name.
    """

    _blocks_undo: bool = True   # Input editing; Ctrl+Z = input undo, not app

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next",   show=False),
    ]

    DEFAULT_CSS = """
    #multi-fasta-dlg {
        width: 72; max-width: 95%; min-width: 56;
        height: auto; max-height: 90%;
        background: $surface; border: solid $primary; padding: 1 2;
    }
    #multi-fasta-title  { background: $primary-darken-2; color: $text;
                          padding: 0 1; margin-bottom: 1; text-align: center; }
    #multi-fasta-info   { height: auto; color: $text-muted;
                          margin-bottom: 1; }
    #multi-fasta-dlg Label { color: $text-muted; margin-top: 1; }
    #multi-fasta-input  { margin-top: 0; margin-bottom: 1; }
    #multi-fasta-status { height: 2; color: $text-muted; }
    #multi-fasta-btns   { align: right middle;  height: 3; margin-top: 1; }
    #multi-fasta-btns Button { margin-right: 1; }
    """

    def __init__(self, count: int, default_name: str = "") -> None:
        super().__init__()
        self.count = int(count)
        self.default_name = default_name

    def compose(self) -> ComposeResult:
        with Vertical(id="multi-fasta-dlg"):
            yield Static(" Multi-record FASTA ", id="multi-fasta-title")
            yield Static(
                f"This FASTA contains [b]{self.count}[/b] records.\n\n"
                f"Load all into a new collection? Each record becomes "
                f"a separate plasmid in the collection. Topology is "
                f"auto-detected per record from its description "
                f"(records mentioning [b]circular[/b] or [b]plasmid[/b] "
                f"load as circular; the rest load as linear).",
                id="multi-fasta-info", markup=True,
            )
            yield Label("Collection name:")
            yield Input(value=self.default_name,
                          placeholder="Collection name",
                          id="multi-fasta-input")
            yield Static("", id="multi-fasta-status", markup=True)
            with Horizontal(id="multi-fasta-btns"):
                yield Button("Create collection",
                              id="btn-multi-fasta-ok", variant="primary")
                yield Button("Cancel", id="btn-multi-fasta-cancel")

    def on_mount(self) -> None:
        self.query_one("#multi-fasta-input", Input).focus()

    @on(Button.Pressed, "#btn-multi-fasta-ok")
    def _ok(self, _) -> None:
        self._submit()

    @on(Input.Submitted, "#multi-fasta-input")
    def _submitted(self, _) -> None:
        self._submit()

    def _submit(self) -> None:
        # Re-use the agent-API normaliser: strips control chars + caps
        # length so a hand-typed `name\x00\n` can't break the panel
        # header / collection-list rendering.
        raw = self.query_one("#multi-fasta-input", Input).value
        name = _normalize_collection_name(raw)
        status = self.query_one("#multi-fasta-status", Static)
        if name is None:
            status.update("[red]Name cannot be empty.[/red]")
            return
        if _collection_name_taken(name):
            status.update(
                f"[red]A collection named '{name}' already exists. "
                f"Pick a different name.[/red]"
            )
            return
        self.dismiss(name)

    @on(Button.Pressed, "#btn-multi-fasta-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class ColorPickerModal(_OneShotDismissScreen, ModalScreen):
    """Pick a display color for a feature-library entry.

    Returns via ``dismiss``:
      * ``None`` — user cancelled.
      * ``{"color": "#RRGGBB", "set_default": False}`` — set entry color.
      * ``{"color": "#RRGGBB", "set_default": True}`` — also save as the
        default for this feature type.
      * ``{"color": None, "set_default": False}`` — clear the override and
        fall back to the type default.

    Three ways to pick a color:

      1. **Curated quick-picks** — the 20-color ``_SWATCHES`` that reuse
         the main map palette.
      2. **xterm 256-color grid** — full 256-cell grid (16 ANSI + 216 cube
         + 24 grayscale). Renders as tiny colored buttons; click one to
         load into the preview.
      3. **Custom hex / index input** — free-form ``#RGB`` / ``#RRGGBB`` /
         ``0..255`` / ``color(N)``. Validated via ``_normalise_color_input``
         which also converts xterm indices to their RGB equivalent so the
         stored value is always a canonical uppercase hex string.

    If the terminal only supports 8/16 colors (``console.color_system`` is
    ``"standard"`` or ``None``), a yellow warning explains that truecolor
    choices will be approximated. The picker still works — you just can't
    visually distinguish similar hex colors on that terminal.
    """

    _blocks_undo: bool = True   # Hex Input editing; Ctrl+Z = input undo

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    # Curated quick-picks — hex-encoded so they round-trip through JSON.
    _SWATCHES: list[str] = [
        "#FF6347", "#FFA500", "#FFD700", "#FFFF00", "#ADFF2F",
        "#7CFC00", "#00FF7F", "#00CED1", "#1E90FF", "#4169E1",
        "#9370DB", "#BA55D3", "#FF69B4", "#FF1493", "#DC143C",
        "#A0522D", "#CD853F", "#20B2AA", "#708090", "#2F4F4F",
    ]

    def __init__(self, feature_type: str, current_color: "str | None") -> None:
        super().__init__()
        self._feature_type = feature_type
        self._current      = current_color or ""
        self._pending:  "str | None" = current_color or None
        self._drag_active: bool = False

    def compose(self) -> ComposeResult:
        type_default = _DEFAULT_TYPE_COLORS.get(self._feature_type, "")
        user_default = _load_feature_colors().get(self._feature_type, "")
        effective_default = _markup_safe_color(
            user_default or type_default or "#808080")

        with Vertical(id="colorpick-dlg"):
            yield Static(f" Pick color for {self._feature_type} ",
                         id="colorpick-title")
            yield Label(
                f"Current: {self._current or '[dim](auto — using type default)[/]'}   "
                f"Type default: [{effective_default}]███[/]",
                markup=True, id="colorpick-current",
            )
            yield Static("", id="colorpick-capability", markup=True)

            with Horizontal(id="colorpick-preview-row"):
                yield Static("", id="colorpick-preview-swatch")
                yield Static("", id="colorpick-preview-label", markup=True)

            with ScrollableContainer(id="colorpick-scroll"):
                yield Label("Curated", classes="colorpick-section-hdr")
                with Horizontal(id="colorpick-row"):
                    # Sweep #36 (2026-05-27): swatches + action
                    # buttons use `_InstantPressButton` so a single
                    # mouse-down fires `Pressed` immediately,
                    # bypassing the Textual focus-transition gate
                    # that swallowed the first click in a real
                    # terminal. User-reported: picking a swatch +
                    # clicking "Auto (clear override)" both took
                    # two physical clicks because the modal wasn't
                    # focused on open. Drop-in replacement — the
                    # `@on(Button.Pressed, ...)` handlers below see
                    # the same `Pressed` message, so dispatch
                    # wiring is unchanged. Same fix that sweep #31
                    # applied to `StrandPickerModal`.
                    for i, hex_col in enumerate(self._SWATCHES):
                        yield _InstantPressButton(
                            "  ", id=f"colorpick-swatch-{i}",
                            classes="colorpick-swatch",
                        )

                yield Label("xterm 256  (click + drag any cell)",
                            classes="colorpick-section-hdr")
                # Single Static replaces the 256 individual Button
                # widgets that used to make up this grid (~2 s mount
                # time → ~70 ms after the swap). Click + drag still
                # work through `on_mouse_down` / `on_mouse_move` →
                # `_cell_index_at`, which now does integer math
                # against the grid's region instead of a widget-tree
                # `get_widget_at` lookup.
                yield _XtermColorGrid(id="colorpick-xterm-grid")

                yield Label("Custom", classes="colorpick-section-hdr")
                with Horizontal(id="colorpick-custom-row"):
                    yield Label("Hex / xterm idx:",
                                id="colorpick-custom-label")
                    yield Input(
                        value=(self._current if _HEX6_RE.match(self._current)
                               else ""),
                        placeholder="#FF6347, F63, 208, or color(208)",
                        id="colorpick-hex-input",
                    )
                    yield _InstantPressButton(
                        "Apply", id="btn-colorpick-apply",
                        variant="primary",
                    )

            yield Static("", id="colorpick-status", markup=True)
            with Horizontal(id="colorpick-btns"):
                yield _InstantPressButton(
                    "Auto (clear override)",
                    id="btn-colorpick-auto",
                )
                yield _InstantPressButton(
                    "Save",
                    id="btn-colorpick-save",
                    variant="primary",
                )
                yield _InstantPressButton(
                    "Save + set as type default",
                    id="btn-colorpick-default",
                    variant="success",
                )
                yield _InstantPressButton(
                    "Cancel", id="btn-colorpick-cancel",
                )

    def on_mount(self) -> None:
        # Paint curated swatches with their own background so the
        # palette is visible at a glance. The xterm grid is now a
        # single `_XtermColorGrid` Static that paints itself in its
        # own `render()` — no per-cell button-style iteration here.
        for i, hex_col in enumerate(self._SWATCHES):
            try:
                btn = self.query_one(f"#colorpick-swatch-{i}", Button)
            except NoMatches:
                continue
            btn.styles.background = hex_col
        self._refresh_capability_warning()
        self._refresh_status()

    def _refresh_capability_warning(self) -> None:
        """Warn the user if the terminal is 8/16-color — they can still
        pick truecolor hex, it'll just be approximated to the nearest
        ANSI when rendered."""
        try:
            cap = self.query_one("#colorpick-capability", Static)
        except NoMatches:
            return
        try:
            sys_name = (self.app.console.color_system or "").lower()
        except Exception:       # noqa: BLE001
            sys_name = ""
        if sys_name in ("truecolor", "256"):
            cap.update(
                f"[dim]Terminal palette: {sys_name} — full range available.[/]"
            )
        else:
            label = sys_name or "unknown / 8-color"
            cap.update(
                f"[yellow]Terminal palette: {label}. Truecolor choices "
                f"will be approximated to the nearest ANSI color.[/]"
            )

    def _refresh_status(self) -> None:
        """Repaint the big preview swatch + hex label. Called whenever
        ``self._pending`` changes — including during a live drag across
        xterm cells."""
        try:
            swatch = self.query_one("#colorpick-preview-swatch", Static)
            label  = self.query_one("#colorpick-preview-label",  Static)
        except NoMatches:
            return
        if self._pending:
            swatch.styles.background = self._pending
            t = Text()
            t.append("Selected: ", style="bold")
            t.append(self._pending, style=self._pending)
            label.update(t)
        else:
            swatch.styles.background = "transparent"
            label.update("[dim]Selected: Auto (use type default)[/]")

    def _set_pending(self, value: "str | None") -> None:
        """Central entry point for any source that changes the pending
        color — keeps the preview swatch in lock-step and clears any
        stale error message in #colorpick-status."""
        if value == self._pending:
            return
        self._pending = value
        self._refresh_status()
        try:
            self.query_one("#colorpick-status", Static).update("")
        except NoMatches:
            pass

    def _cell_index_at(self, sx: int, sy: int) -> "int | None":
        """Hit-test screen coords against the xterm grid. Returns the
        xterm index (0..255) under the cursor, or ``None`` if the point
        is outside any cell. Used by the live-drag preview so the user
        can sweep across the grid with a held click. Integer math
        against the grid's screen region — no widget-tree lookup,
        unlike the legacy 256-button implementation that called
        `self.get_widget_at` on every mouse-move event."""
        try:
            grid = self.query_one("#colorpick-xterm-grid", _XtermColorGrid)
        except NoMatches:
            return None
        region = grid.region
        if region.width == 0 or region.height == 0:
            return None
        if not (region.x <= sx < region.x + region.width
                and region.y <= sy < region.y + region.height):
            return None
        return grid.cell_at(sx - region.x, sy - region.y)

    def on_mouse_down(self, event: MouseDown) -> None:
        """Entering a drag: if the mouse-down lands on an xterm cell,
        arm drag-mode and load that cell's color into the preview
        immediately. Other targets are left alone so regular button
        clicks (Save, Cancel, Apply) still work."""
        if event.button != 1:
            return
        idx = self._cell_index_at(event.screen_x, event.screen_y)
        if idx is not None:
            self._drag_active = True
            self._set_pending(_xterm_index_to_hex(idx))

    def on_mouse_move(self, event: MouseMove) -> None:
        """During a drag, follow the cursor across xterm cells and keep
        the preview in sync. No effect outside of drag-mode."""
        if not self._drag_active:
            return
        idx = self._cell_index_at(event.screen_x, event.screen_y)
        if idx is not None:
            self._set_pending(_xterm_index_to_hex(idx))

    def on_mouse_up(self, event: MouseUp) -> None:
        self._drag_active = False

    @on(Button.Pressed, ".colorpick-swatch")
    def _swatch(self, event: Button.Pressed) -> None:
        btn_id = event.button.id or ""
        if not btn_id.startswith("colorpick-swatch-"):
            return
        try:
            idx = int(btn_id.rsplit("-", 1)[1])
        except ValueError:
            return
        if 0 <= idx < len(self._SWATCHES):
            self._set_pending(self._SWATCHES[idx])

    @on(Click, "#colorpick-xterm-grid")
    def _xterm_grid_click(self, event: Click) -> None:
        """Click on the new single-Static grid — convert to a cell
        index and set the pending color. Mirrors the legacy
        per-button `Button.Pressed` handler that the old design used
        before the 256-button grid was replaced. `on_mouse_down` also
        handles initial click-down for drag continuity, so this
        handler is the keyboard-friendly path for terminals that don't
        deliver a clean MouseDown chain.
        """
        idx = self._cell_index_at(event.screen_x, event.screen_y)
        if idx is not None:
            self._set_pending(_xterm_index_to_hex(idx))

    @on(Button.Pressed, "#btn-colorpick-apply")
    def _apply_custom(self, _) -> None:
        try:
            inp = self.query_one("#colorpick-hex-input", Input)
            status = self.query_one("#colorpick-status", Static)
        except NoMatches:
            return
        raw = inp.value.strip()
        canonical = _normalise_color_input(raw)
        if canonical is None:
            status.update(
                f"[red]Invalid color '{raw}' — use #RGB, #RRGGBB, or 0..255.[/]"
            )
            return
        self._set_pending(canonical)

    @on(Input.Submitted, "#colorpick-hex-input")
    def _hex_submitted(self, _) -> None:
        self._apply_custom(None)

    @on(Button.Pressed, "#btn-colorpick-auto")
    def _auto(self, _) -> None:
        self._set_pending(None)

    @on(Button.Pressed, "#btn-colorpick-save")
    def _save(self, _) -> None:
        self.dismiss({"color": self._pending, "set_default": False})

    @on(Button.Pressed, "#btn-colorpick-default")
    def _save_default(self, _) -> None:
        if not self._pending:
            self.query_one("#colorpick-status", Static).update(
                "[yellow]Pick a specific color before setting as default.[/]"
            )
            return
        self.dismiss({"color": self._pending, "set_default": True})

    @on(Button.Pressed, "#btn-colorpick-cancel")
    def _cancel(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class HmmDbAddEditModal(ModalScreen):
    """Add / edit a custom HMM database catalog entry. Three Inputs:
    display name, source URL, optional version-check URL. On submit,
    dismisses with a normalised entry dict (or None on cancel).

    Mode is `("add", "")` for a new entry or `("edit", entry_id)` to
    pre-fill the form with an existing entry. Builtin entries can
    only have their URL overridden (id + name preserved + builtin
    flag stays True); user-added entries are fully editable.
    """

    # Sweep #28: TextArea/Input inside → Ctrl+Z would leak to canvas.
    _blocks_undo: bool = True

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    DEFAULT_CSS = """
    HmmDbAddEditModal { align: center middle; }
    #hmmdb-form-dlg { width: 84; height: auto; max-height: 90%;
                      background: $surface; border: solid $primary; padding: 1 2; }
    #hmmdb-form-title { background: $primary-darken-2; color: $text;
                        padding: 0 1; margin-bottom: 1; text-align: center; }
    #hmmdb-form-dlg Input { margin-bottom: 1; }
    #hmmdb-form-status { min-height: 1; margin-bottom: 1; }
    #hmmdb-form-btns { height: 3; align: right middle; }
    #hmmdb-form-btns Button { margin-left: 2; }
    """

    def __init__(self, *, mode: str = "add",
                 entry_id: "str | None" = None,
                 **kwargs) -> None:
        super().__init__(**kwargs)
        self._mode = mode
        self._entry_id = entry_id or ""
        self._existing: "dict | None" = None
        if mode == "edit" and entry_id:
            self._existing = _find_hmm_db_entry(entry_id)
        # Sweep #50 double-fire guard.
        self._dismissed: bool = False

    def _dismiss_once(self, result) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(result)

    def compose(self) -> ComposeResult:
        title = ("Add HMM database" if self._mode == "add"
                 else "Edit HMM database")
        with Vertical(id="hmmdb-form-dlg"):
            yield Static(f" {title} ", id="hmmdb-form-title")
            yield Label("Display name (used in the BLAST modal picker):")
            yield Input(
                value=str(self._existing.get("name", "")
                          if self._existing else ""),
                placeholder="e.g. Dfam, my-organism HMMs",
                id="hmmdb-form-name",
                disabled=bool(self._existing
                              and self._existing.get("builtin")),
            )
            yield Label("Source URL (must point at a gzipped HMMER3 "
                        "`.hmm.gz` file):")
            yield Input(
                value=str(self._existing.get("url", "")
                          if self._existing else ""),
                placeholder="https://…/database.hmm.gz",
                id="hmmdb-form-url",
            )
            yield Label("Optional: version-check URL (small file, "
                        "polled for update-detection):")
            yield Input(
                value=str(self._existing.get("version_url", "")
                          if self._existing else ""),
                placeholder="https://…/version (leave blank for "
                            "HEAD-on-source fallback)",
                id="hmmdb-form-version-url",
            )
            yield Label("Optional: short description:")
            yield Input(
                value=str(self._existing.get("description", "")
                          if self._existing else ""),
                placeholder="What's in this database?",
                id="hmmdb-form-desc",
            )
            yield Static("", id="hmmdb-form-status", markup=True)
            with Horizontal(id="hmmdb-form-btns"):
                yield Button("Save", id="btn-hmmdb-form-save",
                             variant="primary")
                yield Button("Cancel", id="btn-hmmdb-form-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#hmmdb-form-name", Input).focus()
        except NoMatches:
            pass

    def _set_status(self, msg: str) -> None:
        try:
            self.query_one("#hmmdb-form-status", Static).update(msg)
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-hmmdb-form-save")
    def _save_btn(self, _) -> None:
        try:
            name = self.query_one("#hmmdb-form-name", Input).value.strip()
            url  = self.query_one("#hmmdb-form-url",  Input).value.strip()
            vurl = self.query_one("#hmmdb-form-version-url",
                                    Input).value.strip()
            desc = self.query_one("#hmmdb-form-desc",
                                    Input).value.strip()
        except NoMatches:
            return
        # For builtin edits, name + id stay; only URL overrides.
        if self._existing and self._existing.get("builtin"):
            entry_id = self._existing["id"]
            display_name = self._existing.get("name") or entry_id
            builtin_flag = True
        else:
            # User add: id = sanitise(name). Reject if id collides.
            if not name:
                self._set_status("[red]Display name is required.[/red]")
                return
            entry_id = _sanitize_hmm_db_id(
                name.replace(" ", "_").lower()
            )
            if entry_id is None:
                self._set_status(
                    "[red]Name can't be reduced to a valid id "
                    "(letters, digits, `_`, `-` only).[/red]"
                )
                return
            display_name = name[:200]
            builtin_flag = False
            if self._mode == "add" and _find_hmm_db_entry(entry_id):
                self._set_status(
                    f"[red]An entry with id [b]{entry_id}[/b] already "
                    f"exists. Pick a different name.[/red]"
                )
                return
            if _hmm_db_name_taken(
                display_name,
                exclude_id=(self._existing or {}).get("id"),
            ):
                self._set_status(
                    f"[red]Display name [b]{display_name}[/b] is "
                    f"taken.[/red]"
                )
                return
        if _sanitize_hmm_db_url(url) is None:
            self._set_status(
                "[red]URL must be a well-formed http(s):// link "
                "with no whitespace.[/red]"
            )
            return
        if vurl and _sanitize_hmm_db_url(vurl) is None:
            self._set_status(
                "[red]Version URL invalid — must be http(s):// or "
                "blank.[/red]"
            )
            return
        normalised = _normalise_hmm_db_entry({
            "id":          entry_id,
            "name":        display_name,
            "url":         url,
            "version_url": vurl,
            "format":      "hmm-gz",
            "builtin":     builtin_flag,
            "description": desc,
        })
        if normalised is None:
            self._set_status(
                "[red]Could not save — entry failed validation.[/red]"
            )
            return
        self._dismiss_once(normalised)

    @on(Button.Pressed, "#btn-hmmdb-form-cancel")
    def _cancel_btn(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)


class MoveCopyToCollectionModal(ModalScreen):
    """Sweep #28: target-collection picker for bulk-move / bulk-copy.

    For mode=="move", the source collection is hidden from the list
    (moving into itself is a no-op). For mode=="copy", the source
    appears in the list labelled as "(duplicate here)" so the user
    can clone marked plasmids in-place inside the active collection —
    each landing gets a " COPY" / " COPY 2" suffix via the same
    collision-rename path used for cross-collection copies.

    Dismisses with the picked target name (str) on confirm, or None
    on cancel. The caller (PlasmidApp's
    `on_library_panel_move_copy_requested` handler) runs the actual
    transactional commit.

    Hardening:
      * `_blocks_undo = True` — modal has Input focus via the search
        field; a stray Ctrl+Z under it would revert the canvas.
      * `_dismiss_once` flag — double-click on Confirm + a worker
        mid-flight can't fire the dismiss callback twice.
      * Default focus on the target-table (cursor highlights the
        first row) so a deliberate Enter confirms the pick. Cancel
        is still one Esc / Tab+Enter away. Pre-sweep the focus
        landed on [Cancel] to defend against a stray Enter on the
        only choice; with same-collection copy now in the list, the
        first row is no longer always the "right" answer, so we
        defer to the user picking.
      * Refuses if the source collection no longer exists (race
        between Space-marking and 'm'/'y' press).
      * Works even when no target collection exists yet — the
        [New collection] button creates an empty one inline (copies
        CollectionsModal's create-and-persist path), then selects it so
        the user can Confirm straight into it. For copy, the source is
        always offered as "duplicate here".
    """

    _blocks_undo: bool = True

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, *, source_collection: str,
                 entry_ids: "list[str]", mode: str,
                 **kwargs) -> None:
        super().__init__(**kwargs)
        self._source = source_collection
        self._entry_ids = list(entry_ids)
        self._mode = mode if mode in ("move", "copy") else "move"
        self._dismissed: bool = False
        self._row_to_name: list[str] = []
        # Re-entrancy guard for the [New collection] button: True while
        # the CollectionNameModal prompt is open so a double-click can't
        # stack two name prompts.
        self._naming: bool = False

    def _dismiss_once(self, result) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(result)

    def compose(self) -> ComposeResult:
        n = len(self._entry_ids)
        verb = "Move" if self._mode == "move" else "Copy"
        plural = "s" if n != 1 else ""
        title = f" {verb} {n} plasmid{plural} to… "
        help_text = (
            "Pick a destination (Enter to confirm, Esc to cancel). "
            "The active collection is offered as "
            "[b]duplicate here[/b] — each duplicate gets a "
            "[b]COPY[/b] suffix."
            if self._mode == "copy"
            else "Pick a target collection (Enter to confirm, "
                 "Esc to cancel):"
        )
        with Vertical(id="movecopy-dlg"):
            yield Static(title, id="movecopy-title", markup=True)
            yield Label(
                f"From collection: [b]{self._source}[/b]",
                id="movecopy-source",
                markup=True,
            )
            yield Label(help_text, id="movecopy-help", markup=True)
            yield DataTable(id="movecopy-table",
                              cursor_type="row",
                              show_cursor=True,
                              zebra_stripes=True)
            yield Static("", id="movecopy-status", markup=True)
            with Horizontal(id="movecopy-btns"):
                yield Button("Confirm", id="btn-movecopy-go",
                             variant="primary")
                yield Button("New collection", id="btn-movecopy-newcoll")
                yield Button("Cancel", id="btn-movecopy-cancel")

    def on_mount(self) -> None:
        self._repopulate()
        try:
            self.query_one("#movecopy-table", DataTable).focus()
        except NoMatches:
            pass

    def _set_status(self, msg: str) -> None:
        try:
            self.query_one("#movecopy-status", Static).update(msg)
        except NoMatches:
            pass

    def _repopulate(self) -> None:
        try:
            table = self.query_one("#movecopy-table", DataTable)
        except NoMatches:
            return
        table.clear(columns=True)
        table.add_columns("Collection", "Plasmids")
        self._row_to_name = []
        try:
            colls = list(_iter_collections_readonly())
        except (OSError, RuntimeError) as exc:
            self._set_status(
                f"[red]Couldn't load collections: {_scrub_path(str(exc))}"
                f"[/red]"
            )
            return
        seen = 0
        # For copy mode, the source collection is offered as
        # "(duplicate here)" so the user can clone marked entries
        # in-place — same collision-rename path as cross-collection
        # copy, no special-case at the commit layer beyond allowing
        # source==target for copy.
        allow_self = (self._mode == "copy")
        for c in colls:
            name = (c.get("name") or "").strip()
            if not name:
                continue
            is_source = (name == self._source)
            if is_source and not allow_self:
                continue
            n_plas = len(c.get("plasmids") or [])
            label = (f"{name}  [dim](duplicate here)[/dim]"
                     if is_source else name)
            table.add_row(label, str(n_plas), key=name)
            self._row_to_name.append(name)
            seen += 1
        if seen == 0:
            self._set_status(
                "[yellow]No target collections yet — press "
                "[b]New collection[/b] to make one.[/yellow]"
            )
        else:
            verb = "target" if self._mode == "move" else "destination"
            self._set_status(
                f"[dim]{seen} {verb}"
                f"{'s' if seen != 1 else ''} available.[/dim]"
            )

    def _selected_name(self) -> "str | None":
        try:
            table = self.query_one("#movecopy-table", DataTable)
        except NoMatches:
            return None
        row = table.cursor_row
        if row is None or row < 0 or row >= len(self._row_to_name):
            return None
        return self._row_to_name[row]

    @on(Button.Pressed, "#btn-movecopy-go")
    def _go_btn(self, _) -> None:
        target = self._selected_name()
        if not target:
            self._set_status("[red]Pick a target collection first.[/red]")
            return
        # Re-check the source still exists (race between marking and
        # confirming — a parallel agent / second-window collection
        # delete could have removed it).
        if _find_collection(self._source) is None:
            self._set_status(
                f"[red]Source collection {self._source!r} no longer "
                f"exists — the operation is aborted.[/red]"
            )
            return
        if _find_collection(target) is None:
            self._set_status(
                f"[red]Target collection {target!r} disappeared mid-"
                f"flight. Refresh and retry.[/red]"
            )
            self._repopulate()
            return
        self._dismiss_once(target)

    @on(DataTable.RowSelected, "#movecopy-table")
    def _row_selected(self, _event) -> None:
        # Enter on a row = Confirm. Same UX as the other Picker modals.
        self._go_btn(None)

    @on(Button.Pressed, "#btn-movecopy-newcoll")
    def _new_collection_btn(self, _) -> None:
        """Create a brand-new empty target collection without leaving the
        picker. Copies CollectionsModal._save's create-and-persist path
        (normalise via CollectionNameModal → collision check → append an
        empty collection → `_save_collections` → repopulate). On success
        the new collection is selected so Enter / Confirm lands on it."""
        if self._naming:                       # double-click guard
            return

        def _named(name: "str | None") -> None:
            self._naming = False
            if not name:                        # cancelled / empty
                return
            # Collision check mirrors CollectionsModal._save. `name` is
            # already normalised by CollectionNameModal.
            if _collection_name_taken(name):
                self._set_status(
                    f"[red]A collection named '{name}' already exists — "
                    f"pick it from the list above instead.[/red]"
                )
                return
            try:
                existing = _load_collections()
                existing.append({
                    "name":        name,
                    "description": "",
                    "plasmids":    [],
                    "saved":       _date.today().isoformat(),
                })
                _save_collections(existing)
            except (OSError, RuntimeError) as exc:
                _notify_save_failure(self.app, "Collections", exc)
                return
            self._repopulate()
            self._select_row_by_name(name)
            verb = "move" if self._mode == "move" else "copy"
            self._set_status(
                f"[green]Created '{name}'. Confirm to {verb} the marked "
                f"plasmid(s) here.[/green]"
            )

        self._naming = True
        try:
            self.app.push_screen(
                CollectionNameModal("New collection"), callback=_named,
            )
        except Exception:
            # Push should never fail for a mounted modal, but if it does
            # don't leave the button permanently dead-locked.
            self._naming = False
            _log.debug("movecopy: new-collection prompt push failed",
                       exc_info=True)

    def _select_row_by_name(self, name: str) -> None:
        """Move the table cursor onto the row whose key is `name` and
        focus the table, so a fresh Enter confirms it. No-op if the row
        isn't present (e.g. it was filtered out as the move source)."""
        try:
            table = self.query_one("#movecopy-table", DataTable)
        except NoMatches:
            return
        try:
            idx = self._row_to_name.index(name)
        except ValueError:
            return
        try:
            table.move_cursor(row=idx)
            table.focus()
        except Exception:
            _log.debug("movecopy: select-row-by-name failed", exc_info=True)

    @on(Button.Pressed, "#btn-movecopy-cancel")
    def _cancel_btn(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)


class SplitFeatureModal(ModalScreen):
    """Sweep #29: divide one feature into N sub-features that
    share a `feature_group=<uuid>` qualifier.

    The user types a per-line breakdown into a multi-line input:
      ``<rel_start>-<rel_end> <label> <color> <strand>``

    e.g.:
      ``0-4 GCGC #888888 0``
      ``4-10 Esp3I #FF3333 1``
      ``10-11 N #666666 0``
      ``11-15 AATG #00CC66 1``

    Preset templates (button row above the input) pre-fill the
    text area with common Golden Gate adapter shapes so users
    don't have to remember the rel coords. On Save, the lines
    are parsed + validated; failures surface inline status text
    and DON'T dismiss — the user keeps editing.

    Returns `{"members": [<validated member dicts>]}` on Save,
    `None` on Cancel."""

    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #split-dlg {
        width: 100; max-width: 95%; height: auto; max-height: 90%;
        background: $surface; border: solid $primary;
        padding: 1 2;
    }
    #split-title {
        background: $primary-darken-2; color: $text;
        padding: 0 1; margin-bottom: 1;
        text-align: center; text-style: bold;
    }
    #split-dlg Label { color: $text-muted; margin-top: 1; }
    #split-presets { height: 3; margin-top: 0; }
    #split-presets Button { margin-right: 1; min-width: 16; }
    #split-rows { height: 12; min-height: 6; margin-top: 0; }
    #split-help { color: $text-muted; padding: 0 1;
                  height: auto; margin-top: 1; }
    #split-status { height: auto; margin-top: 1; padding: 0 1; }
    #split-btns { height: 3; margin-top: 1; align: right middle; }
    #split-btns Button { margin-right: 1; min-width: 12; }
    """

    def __init__(self, feat_label: str, span: int,
                 *, base_color: "str | None" = None,
                 base_strand: int = 1,
                 base_type: str = "misc_feature") -> None:
        super().__init__()
        self._feat_label  = str(feat_label or "")
        self._span        = int(span)
        self._base_color  = base_color
        self._base_strand = base_strand
        self._base_type   = base_type
        self._dismissed: bool = False

    def _dismiss_once(self, result) -> None:
        if self._dismissed:
            return
        self._dismissed = True
        self.dismiss(result)

    def compose(self) -> ComposeResult:
        with Vertical(id="split-dlg"):
            yield Static(
                f" Split '{self._feat_label or '(unnamed)'}' "
                f"into sub-features ({self._span} bp total) ",
                id="split-title",
            )
            yield Label("Presets:")
            with Horizontal(id="split-presets"):
                yield Button(
                    "Equal thirds", id="btn-split-preset-thirds",
                    tooltip=(
                        "Pre-fill the rows with 3 equal-sized "
                        "sub-features sharing the parent's type."
                    ),
                )
                yield Button(
                    "5'-Esp3I-AATG adapter",
                    id="btn-split-preset-esp3i",
                    tooltip=(
                        "Pre-fill with GCGC pad + Esp3I + N + "
                        "AATG overhang (Golden Gate 5'-adapter "
                        "shape). Requires the feature to be "
                        "exactly 15 bp."
                    ),
                )
                yield Button(
                    "Clear", id="btn-split-preset-clear",
                    tooltip="Erase the rows.",
                )
            yield Label(
                "Sub-features  (one per line: "
                "`<rel_start>-<rel_end> <label> <color> <strand>`):"
            )
            yield TextArea("", id="split-rows")
            yield Static(
                "[dim]Strand: `1` = forward, `-1` = reverse, "
                "`0` = arrowless, `2` = double. Color: `#RRGGBB` "
                "or leave empty for the type default. "
                "Rel coords are within `[0, span)` half-open. "
                "Members CAN overlap and CAN leave gaps "
                "(unannotated bases stay unannotated).[/dim]",
                id="split-help", markup=True,
            )
            yield Static("", id="split-status", markup=True)
            with Horizontal(id="split-btns"):
                yield Button("Save", id="btn-split-save",
                             variant="primary")
                yield Button("Cancel", id="btn-split-cancel")

    @on(Button.Pressed, "#btn-split-preset-thirds")
    def _preset_thirds(self, _) -> None:
        n = self._span
        if n < 3:
            self._set_status(
                "[red]Span too short for 3 equal thirds "
                "(need ≥ 3 bp).[/red]"
            )
            return
        third = n // 3
        rows = [
            f"0-{third} part_1 #888888 {self._base_strand}",
            f"{third}-{2 * third} part_2 #BBBBBB "
            f"{self._base_strand}",
            f"{2 * third}-{n} part_3 #EEEEEE {self._base_strand}",
        ]
        self._set_rows("\n".join(rows))

    @on(Button.Pressed, "#btn-split-preset-esp3i")
    def _preset_esp3i(self, _) -> None:
        if self._span != 15:
            self._set_status(
                f"[red]Esp3I→AATG preset wants exactly 15 bp; "
                f"this feature is {self._span} bp.[/red]"
            )
            return
        rows = [
            "0-4 GCGC #888888 0",
            "4-10 Esp3I #FF3333 1",
            "10-11 N #666666 0",
            "11-15 AATG #00CC66 1",
        ]
        self._set_rows("\n".join(rows))

    @on(Button.Pressed, "#btn-split-preset-clear")
    def _preset_clear(self, _) -> None:
        self._set_rows("")
        self._set_status("")

    def _set_rows(self, text: str) -> None:
        try:
            ta = self.query_one("#split-rows", TextArea)
        except NoMatches:
            return
        ta.text = text

    def _set_status(self, msg: str) -> None:
        try:
            self.query_one("#split-status", Static).update(msg)
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-split-save")
    def _save(self, _) -> None:
        try:
            text = self.query_one("#split-rows", TextArea).text
        except NoMatches:
            return
        # Sweep #29 hardening (2026-05-26): copy-paste of a huge
        # blob into the TextArea (binary file, syslog, etc.) is a
        # real risk. Cap input size before splitting on newlines
        # so the parser can't be DoS'd; cap line count downstream
        # via `_MAX_GROUP_MEMBERS`.
        if len(text) > 32_768:
            self._set_status(
                "[red]Input too long (>32 KB) — paste a smaller "
                "breakdown or split into multiple groups.[/red]"
            )
            return
        # Parse each non-empty line into a member dict. Format:
        # `<rel_start>-<rel_end> <label> [color] [strand]`. Color
        # and strand are optional; missing color → palette
        # default (None), missing strand → base_strand.
        members: list[dict] = []
        for lineno, raw in enumerate(text.splitlines(), 1):
            # Strip control chars on the way in so a paste with
            # ANSI escape sequences / tabs / nulls can't smuggle
            # into the parser. Tabs collapse to spaces below.
            cleaned = _CONTROL_CHARS_RE.sub("", raw)
            line = cleaned.strip()
            if not line or line.startswith("#"):
                continue
            # Cap individual line length too — a line of
            # 100 KB of "label" text shouldn't be tolerated.
            if len(line) > 1024:
                self._set_status(
                    f"[red]Line {lineno} too long "
                    f"({len(line)} chars > 1024) — wrap or "
                    f"shorten.[/red]"
                )
                return
            # Split on whitespace; first token is the rel range.
            parts = line.split()
            if not parts:
                continue
            rng = parts[0]
            if "-" not in rng:
                self._set_status(
                    f"[red]Line {lineno}: expected "
                    f"`<rel_start>-<rel_end>`, got '{rng}'.[/red]"
                )
                return
            try:
                rs_str, re_str = rng.split("-", 1)
                rs = int(rs_str)
                re_ = int(re_str)
            except ValueError:
                self._set_status(
                    f"[red]Line {lineno}: rel coords must be "
                    f"ints (got '{rng}').[/red]"
                )
                return
            if not (0 <= rs < re_ <= self._span):
                self._set_status(
                    f"[red]Line {lineno}: range {rs}-{re_} "
                    f"out of `[0, {self._span}]`.[/red]"
                )
                return
            # Remaining tokens: label, color, strand (in order).
            tail = parts[1:]
            label = tail[0] if len(tail) >= 1 else ""
            color = tail[1] if len(tail) >= 2 else ""
            strand_str = tail[2] if len(tail) >= 3 else str(
                self._base_strand,
            )
            try:
                strand = int(strand_str)
            except ValueError:
                strand = self._base_strand
            if strand not in (-1, 0, 1, 2):
                strand = self._base_strand
            members.append({
                "rel_start": rs,
                "rel_end":   re_,
                "feature_type": self._base_type,
                "label":     label,
                "color":     color if color else None,
                "strand":    strand,
                "qualifiers": {},
                "description": "",
            })
        if not members:
            self._set_status(
                "[red]No sub-features defined — type at least "
                "one row first (or use a preset).[/red]"
            )
            return
        # Re-validate via the shared helper so the same rules
        # apply as for library entries. Catches any edge case
        # the line-by-line parse missed.
        try:
            members = _validate_group_members(members, self._span)
        except ValueError as exc:
            self._set_status(
                f"[red]Validation failed: {exc}[/red]"
            )
            return
        self._dismiss_once({"members": members})

    @on(Button.Pressed, "#btn-split-cancel")
    def _cancel(self, _) -> None:
        self._dismiss_once(None)

    def action_cancel(self) -> None:
        self._dismiss_once(None)


class PlasmidStatusPickerModal(_OneShotDismissScreen, ModalScreen):
    """Tiny modal to set the workflow status on a library entry.

    Five RadioButtons (the four canonical statuses + "no status").
    Dismisses with the chosen status string (one of
    `_PLASMID_STATUS_VALUES` or "" for cleared) or None on cancel.

    The picker doesn't persist on its own — caller is responsible
    for routing the result through `_save_library` so the on-disk
    entry + active-collection mirror stay in sync.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    DEFAULT_CSS = """
    #status-dlg {
        width: 56; height: auto; max-height: 22;
        background: $surface;
        border: solid $primary;
        padding: 1 2;
    }
    #status-title {
        height: 1;
        text-style: bold;
        color: $accent;
        margin-bottom: 1; text-align: center;
    }
    #status-current { color: $text-muted; margin-bottom: 1; }
    #status-radios { height: auto; padding: 0 1; }
    #status-btns { height: 3; margin-top: 1; align: right middle; }
    #status-btns Button { margin-right: 1; }
    """

    def __init__(self, plasmid_name: str, current_status: str = "") -> None:
        super().__init__()
        self._plasmid_name = str(plasmid_name or "")
        self._current      = _sanitize_plasmid_status(current_status)

    def compose(self) -> ComposeResult:
        with Vertical(id="status-dlg"):
            yield Static(
                f"Set status: {self._plasmid_name or '(unnamed)'}",
                id="status-title",
            )
            yield Static(
                f"Current: {self._current or '(none)'}",
                id="status-current",
            )
            with RadioSet(id="status-radios"):
                # Order: workflow-natural (designing → verified) plus
                # an explicit "(none)" sentinel so the user can clear
                # a previous assignment without leaving the modal.
                yield RadioButton(
                    "(none)", id="status-radio-none",
                    value=(self._current == ""),
                )
                for s in _PLASMID_STATUS_VALUES:
                    color = _PLASMID_STATUS_COLORS.get(s, "white")
                    # Embed the colour swatch in the label so the
                    # radio reads as the colour the user will see
                    # in the library table.
                    yield RadioButton(
                        f"[{color}]●[/]  {s}",
                        id=f"status-radio-{s.lower()}",
                        value=(self._current == s),
                    )
            with Horizontal(id="status-btns"):
                yield Button("Save",   id="btn-status-save",
                             variant="primary")
                yield Button("Cancel", id="btn-status-cancel")

    def on_mount(self) -> None:
        # Focus the radio set so up/down navigates immediately.
        try:
            self.query_one("#status-radios", RadioSet).focus()
        except NoMatches:
            pass

    def _read_status(self) -> str:
        for s in ("",) + _PLASMID_STATUS_VALUES:
            rid = ("status-radio-none" if s == ""
                    else f"status-radio-{s.lower()}")
            try:
                if self.query_one(f"#{rid}", RadioButton).value:
                    return s
            except NoMatches:
                continue
        return self._current

    @on(Button.Pressed, "#btn-status-save")
    def _save(self, _) -> None:
        chosen = _sanitize_plasmid_status(self._read_status())
        if chosen == self._current:
            # No-op — same as cancel so the caller doesn't bother
            # re-writing the library.
            self.dismiss(None)
            return
        self.dismiss(chosen)

    @on(Button.Pressed, "#btn-status-cancel")
    def _cancel(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class PartEditModal(_OneShotDismissScreen, ModalScreen):
    """View / edit a parts-bin row.

    Read-only by default — `Edit` flips the form editable, `Save`
    commits and dismisses with ``{"idx": <row idx>, "entry": <new
    dict>}``, `Cancel` dismisses with ``None``.

    Grammar is locked: changing it invalidates the type / position /
    overhang semantics. Users who need to migrate a part to a
    different grammar should delete + re-create through the
    Domesticator. Type changes auto-fill position + overhangs from
    the matching grammar position so the common "re-tag" edit
    doesn't require manual oh re-entry.

    Sequence + overhang edits trigger re-derivation of ``primed_seq``
    / ``cloned_seq`` so the Copy Primed / Copy Cloned actions on the
    Parts Bin keep serving the right amplicon. Primer edits re-run
    the primer3 Tm calc (or fall back to 0.0 when primer3 is
    unavailable).
    """

    _blocks_undo: bool = True   # Input editing; Save mutates parts bin

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "app.focus_next", "Next", show=False),
    ]

    # IUPAC alphabet for DNA. Pre-built once at class scope so every
    # save doesn't rebuild the lookup set. `frozenset` makes accidental
    # mutation a TypeError.
    _VALID_IUPAC = frozenset("ACGTRYWSMKBDHVN")

    DEFAULT_CSS = """
    #partedit-dlg {
        width: 108; max-width: 95%; height: auto; max-height: 36;
        background: $surface;
        border: solid $accent;
        padding: 1 2;
    }
    #partedit-title {
        height: 1;
        text-style: bold;
        background: $accent;
        color: $text;
        padding: 0 1;
        margin-bottom: 1; text-align: center;
    }
    /* Body sizes to content (auto), with overflow scroll as a safety
       net for sub-baseline terminals. */
    #partedit-body { height: auto; padding: 0 1; overflow-y: auto; }
    #partedit-body Label { color: $text-muted; margin: 0; height: 1; }
    /* Each label+input pair is 4 rows: 1 label + 3 input/select.
       No margin-top between rows — borders give enough separation. */
    #partedit-row1, #partedit-row2,
    #partedit-row3, #partedit-row4 { height: 4; }
    /* Row 1: Name (3fr) | Grammar (2fr) | Type (2fr) — Name leads because it's
       the field users actually read + edit, and a long part name was cramped
       at 2fr (too narrow to see the words); grammar is a real Select so a part
       can be re-tagged to a different cloning grammar from this modal. */
    #partedit-row1 #partedit-name-col    { width: 3fr; padding-right: 1; }
    #partedit-row1 #partedit-grammar-col { width: 2fr; padding-right: 1; }
    #partedit-row1 #partedit-type-col    { width: 2fr; }
    /* Row 2: 5'OH | 3'OH | Position. Position widens to 2fr because
       its values ("Pos 3-4 (CDS)") are wider than 4-bp overhangs. */
    #partedit-row2 #partedit-oh5-col      { width: 1fr; padding-right: 1; }
    #partedit-row2 #partedit-oh3-col      { width: 1fr; padding-right: 1; }
    #partedit-row2 #partedit-position-col { width: 2fr; }
    /* Rows 3 / 4: pairs of equal-width columns. */
    #partedit-row3 Vertical,
    #partedit-row4 Vertical { width: 1fr; padding-right: 1; }
    #partedit-seq { height: 6; border: solid $primary-darken-1; margin-top: 1; }
    #partedit-status { height: 1; padding: 0 1; }
    #partedit-btns { align: right middle;  height: 3; margin-top: 1; }
    #partedit-btns Button { margin-right: 1; }
    """

    def __init__(self, idx: int, part: dict) -> None:
        super().__init__()
        self._idx = idx
        # Defensive copy so partial edits don't leak back if the user
        # cancels — the committed dict is rebuilt from the form widgets
        # at Save time, so this snapshot only matters for the read-only
        # initial render.
        self._part = dict(part)
        self._editing = False
        # Grammar id is the source of truth — `_grammar` and
        # `_position_index` are derived caches refreshed by
        # `_apply_grammar` whenever the user picks a different
        # grammar from the dropdown. Storing the id (not the dict)
        # keeps a stale grammar from leaking after a re-pick.
        self._grammar_id: str = self._part.get("grammar", "gb_l0") or "gb_l0"
        self._grammar: dict = {}
        self._position_index: dict[str, dict] = {}
        self._apply_grammar(self._grammar_id, refresh_widgets=False)

    def _apply_grammar(self, gid: str, *, refresh_widgets: bool) -> None:
        """Snap `_grammar` + `_position_index` to the grammar with id
        ``gid``. Falls back to ``gb_l0`` for an unknown id (e.g., a
        custom grammar that was deleted while the modal was open).

        ``refresh_widgets=True`` rebuilds the Type select and refreshes
        the position / overhang inputs from the new grammar's matching
        position. Pass False from `__init__` (widgets don't exist
        yet) and True from the grammar-change handler.
        """
        all_g = _all_grammars()
        new_grammar = all_g.get(gid) or _BUILTIN_GRAMMARS["gb_l0"]
        self._grammar_id = gid if gid in all_g else "gb_l0"
        self._grammar = new_grammar
        self._position_index = {}
        for pos in new_grammar.get("positions", []):
            ptype = pos.get("type")
            if (isinstance(ptype, str) and ptype
                    and ptype not in self._position_index):
                self._position_index[ptype] = pos
        if not refresh_widgets:
            return

        # Preserve the user's current type if it survives the grammar
        # change; otherwise default to the first option.
        try:
            type_sel = self.query_one("#partedit-type", Select)
            cur_val = type_sel.value
        except NoMatches:
            return
        cur_type = (cur_val if isinstance(cur_val, str)
                    and cur_val != Select.BLANK else "")
        opts, _ = self._type_options(cur_type)
        new_type = (cur_type if cur_type and cur_type in self._position_index
                    else (opts[0][1] if opts else ""))
        type_sel.set_options(opts)
        if new_type:
            type_sel.value = new_type
        # Refresh position / overhang inputs from the new grammar's
        # matching position so the form is internally consistent
        # immediately after the grammar pick.
        pos = self._position_index.get(new_type)
        if pos is not None:
            for sel, key, upper in (
                ("#partedit-position", "name", False),
                ("#partedit-oh5",       "oh5",  True),
                ("#partedit-oh3",       "oh3",  True),
            ):
                try:
                    val = str(pos.get(key, "") or "")
                    self.query_one(sel, Input).value = (
                        val.upper() if upper else val
                    )
                except NoMatches:
                    pass

    def _type_options(self, current_type: "str | None" = None) -> tuple[list[tuple[str, str]], str]:
        """Build ``(label, ptype)`` options for the Type select keyed
        off the active grammar's positions. ``current_type`` (if
        non-empty + not in the grammar) is added with a ``(legacy)``
        suffix so a Save round-trips it without silently rewriting
        the field. A grammar with no positions yields a placeholder
        so the Select widget composes (Select with ``allow_blank=False``
        and an empty list raises at mount)."""
        opts: list[tuple[str, str]] = []
        for ptype, pos in self._position_index.items():
            label = (
                f"{ptype}  ({pos.get('name','?')}: "
                f"{pos.get('oh5','')}→{pos.get('oh3','')})"
            )
            opts.append((label, ptype))
        current = current_type if current_type is not None else (
            self._part.get("type", "") or ""
        )
        if current and current not in self._position_index:
            opts.insert(0, (f"{current} (legacy)", current))
        if not opts:
            # Pathological grammar (no positions, no part type). Surface
            # an inert placeholder rather than crashing the Select.
            opts = [("(no types defined)", "")]
        default = (current
                   if current and any(v == current for _, v in opts)
                   else opts[0][1])
        return opts, default

    def _validate_iupac_chars(self, label: str, value: str,
                               status_widget) -> bool:
        """Render a red status + return False if ``value`` contains
        any non-IUPAC bases. Empty string passes. Used by every DNA
        field on save (sequence, overhangs, primers) so the error
        format stays consistent and the validation block doesn't
        repeat 4 times in `_on_save`."""
        if not value:
            return True
        bad = [c for c in value if c not in self._VALID_IUPAC]
        if not bad:
            return True
        if status_widget is not None:
            status_widget.update(
                f"[red]{label} has invalid bases: "
                f"{''.join(sorted(set(bad)))[:10]}[/red]"
            )
        return False

    def _primer_label(self, base: str, tm: object) -> str:
        """Format a primer field label including its current Tm so
        the user knows what Tm they're editing against. ``tm`` is
        coerced to float when possible; non-numeric / zero values
        render as a plain label without a Tm suffix."""
        try:
            tm_f = float(tm)  # type: ignore[arg-type] # handles int, str-of-float, np.float64
        except (TypeError, ValueError):
            tm_f = 0.0
        if tm_f > 0:
            return f"{base} (Tm {tm_f:.1f}°C):"
        return f"{base}:"

    def compose(self) -> ComposeResult:
        type_options, default_type = self._type_options()
        # Grammar dropdown — the canonical option list shared with the
        # Domesticator + future grammar pickers. The part's stored
        # grammar is selected by default; if it's not in the registry
        # (e.g., a custom grammar deleted since the part was saved)
        # we splice it in with a "(missing)" suffix so a Save still
        # round-trips the value rather than silently rewriting it.
        grammar_options = list(_grammar_dropdown_options())
        if not any(gid == self._grammar_id for _, gid in grammar_options):
            grammar_options.insert(
                0, (f"{self._grammar_id} (missing)", self._grammar_id),
            )
        p = self._part
        # Title interpolates the part name verbatim. Static defaults to
        # markup=True, which would interpret a name like "[red]X" as
        # Rich markup. markup=False renders the name literally.
        with Vertical(id="partedit-dlg"):
            yield Static(f" Part: {p.get('name', '?')} ",
                         id="partedit-title", markup=False)
            with ScrollableContainer(id="partedit-body"):
                # Row 1 — identity: Name | Grammar | Type
                with Horizontal(id="partedit-row1"):
                    with Vertical(id="partedit-name-col"):
                        yield Label("Name:")
                        yield Input(value=p.get("name", ""),
                                     id="partedit-name", disabled=True)
                    with Vertical(id="partedit-grammar-col"):
                        yield Label("Cloning grammar:")
                        yield Select(grammar_options,
                                      value=self._grammar_id,
                                      id="partedit-grammar",
                                      allow_blank=False, disabled=True)
                    with Vertical(id="partedit-type-col"):
                        yield Label("Type:")
                        yield Select(type_options, value=default_type,
                                     id="partedit-type",
                                     allow_blank=False, disabled=True)
                # Row 2 — cloning context: 5'OH | 3'OH | Position
                with Horizontal(id="partedit-row2"):
                    with Vertical(id="partedit-oh5-col"):
                        yield Label("5' overhang:")
                        yield Input(value=p.get("oh5", ""),
                                     id="partedit-oh5", disabled=True)
                    with Vertical(id="partedit-oh3-col"):
                        yield Label("3' overhang:")
                        yield Input(value=p.get("oh3", ""),
                                     id="partedit-oh3", disabled=True)
                    with Vertical(id="partedit-position-col"):
                        yield Label("Position:")
                        yield Input(value=p.get("position", ""),
                                     id="partedit-position", disabled=True)
                # Row 3 — vector: Backbone | Selection marker
                with Horizontal(id="partedit-row3"):
                    with Vertical():
                        yield Label("Backbone:")
                        yield Input(value=p.get("backbone", ""),
                                     id="partedit-backbone", disabled=True)
                    with Vertical():
                        yield Label("Selection marker:")
                        yield Input(value=p.get("marker", ""),
                                     id="partedit-marker", disabled=True)
                # Row 4 — primers (Tms shown in labels)
                with Horizontal(id="partedit-row4"):
                    with Vertical():
                        yield Label(self._primer_label("Forward primer",
                                                         p.get("fwd_tm")))
                        yield Input(value=p.get("fwd_primer", ""),
                                     id="partedit-fwd", disabled=True)
                    with Vertical():
                        yield Label(self._primer_label("Reverse primer",
                                                         p.get("rev_tm")))
                        yield Input(value=p.get("rev_primer", ""),
                                     id="partedit-rev", disabled=True)
                seq = p.get("sequence", "")
                seq_ta = TextArea(seq, id="partedit-seq",
                                   read_only=True, soft_wrap=True)
                seq_ta.border_title = (
                    f"Insert sequence  (5'→3', {len(seq):,} bp)"
                )
                yield seq_ta
            yield Static("", id="partedit-status", markup=True)
            with Horizontal(id="partedit-btns"):
                yield Button("Edit",   id="btn-partedit-edit",
                             variant="primary")
                yield Button("Save",   id="btn-partedit-save",
                             variant="success", disabled=True)
                yield Button("Cancel", id="btn-partedit-cancel")

    def _set_editing(self, on: bool) -> None:
        """Toggle every form field between read-only and editable.
        TextArea uses `read_only`; everything else uses `disabled` —
        Inputs in the disabled state still show their value clearly,
        which is the right read-mode look."""
        self._editing = on
        for sel in (
            "#partedit-name", "#partedit-grammar", "#partedit-type",
            "#partedit-oh5", "#partedit-oh3",
            "#partedit-position",
            "#partedit-backbone", "#partedit-marker",
            "#partedit-fwd", "#partedit-rev",
        ):
            try:
                self.query_one(sel).disabled = not on
            except NoMatches:
                pass
        try:
            self.query_one("#partedit-seq", TextArea).read_only = not on
        except NoMatches:
            pass
        try:
            self.query_one("#btn-partedit-edit", Button).display = not on
            self.query_one("#btn-partedit-save", Button).disabled = not on
        except NoMatches:
            pass
        if on:
            try:
                self.query_one("#partedit-name", Input).focus()
            except NoMatches:
                pass

    @on(Button.Pressed, "#btn-partedit-edit")
    def _on_edit(self) -> None:
        self._set_editing(True)
        try:
            self.query_one("#partedit-status", Static).update(
                "[dim]Edit mode — make changes and press Save.[/dim]"
            )
        except NoMatches:
            pass

    @on(Select.Changed, "#partedit-grammar")
    def _on_grammar_changed(self, event: Select.Changed) -> None:
        """Grammar change → swap the active grammar, rebuild the Type
        select, and refresh position / overhangs from the new
        grammar's matching position. No-op in read-only mode and on
        the initial Select.Changed that fires during compose."""
        if not self._editing:
            return
        new_gid = (event.value
                    if isinstance(event.value, str) else "")
        if not new_gid or new_gid == self._grammar_id:
            return
        self._apply_grammar(new_gid, refresh_widgets=True)
        try:
            self.query_one("#partedit-status", Static).update(
                "[dim]Grammar updated — overhangs refreshed for the "
                "new grammar.[/dim]"
            )
        except NoMatches:
            pass

    @on(Select.Changed, "#partedit-type")
    def _on_type_changed(self, event: Select.Changed) -> None:
        """Type change → prefill position + overhangs from the
        grammar's matching position so re-tagging a part doesn't
        require manual oh re-entry. No-op while in read-only mode
        (Select.Changed fires on initial render too)."""
        if not self._editing:
            return
        ptype = event.value if isinstance(event.value, str) else ""
        pos = self._position_index.get(ptype) if ptype else None
        if not pos:
            return
        try:
            pos_inp = self.query_one("#partedit-position", Input)
            oh5_inp = self.query_one("#partedit-oh5",      Input)
            oh3_inp = self.query_one("#partedit-oh3",      Input)
        except NoMatches:
            return
        pos_inp.value = pos.get("name", pos_inp.value)
        oh5_inp.value = (pos.get("oh5", "") or "").upper()
        oh3_inp.value = (pos.get("oh3", "") or "").upper()

    @on(Button.Pressed, "#btn-partedit-cancel")
    def _on_cancel(self) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        # Required so the `escape` binding actually closes the modal —
        # `Binding("escape", "cancel", …)` dispatches to `action_cancel`,
        # which `ModalScreen` doesn't supply by default. Without this,
        # Escape would be silently swallowed.
        self.dismiss(None)

    @on(Button.Pressed, "#btn-partedit-save")
    def _on_save(self) -> None:
        if not self._editing:
            return
        try:
            name     = self.query_one("#partedit-name",     Input).value
            ptype    = self.query_one("#partedit-type",     Select).value
            oh5      = self.query_one("#partedit-oh5",      Input).value
            oh3      = self.query_one("#partedit-oh3",      Input).value
            position = self.query_one("#partedit-position", Input).value
            backbone = self.query_one("#partedit-backbone", Input).value
            marker   = self.query_one("#partedit-marker",   Input).value
            fwd      = self.query_one("#partedit-fwd",      Input).value
            rev      = self.query_one("#partedit-rev",      Input).value
            seq_text = self.query_one("#partedit-seq",      TextArea).text
        except NoMatches:
            return
        try:
            status = self.query_one("#partedit-status", Static)
        except NoMatches:
            status = None

        # ── Sanitise form values ────────────────────────────────────
        clean_name = _sanitize_label(name, max_len=200)
        if not clean_name:
            if status:
                status.update("[red]Name cannot be empty.[/red]")
            return
        clean_ptype = _sanitize_label(str(ptype) if ptype else "", max_len=64)
        if not clean_ptype:
            if status:
                status.update("[red]Type cannot be empty.[/red]")
            return
        clean_seq = "".join(str(seq_text or "").split()).upper()
        clean_oh5 = "".join(str(oh5 or "").split()).upper()
        clean_oh3 = "".join(str(oh3 or "").split()).upper()
        clean_fwd = "".join(str(fwd or "").split()).upper()
        clean_rev = "".join(str(rev or "").split()).upper()

        # ── DNA validation (one helper, four call sites) ───────────
        for label_, val in (
            ("Sequence",       clean_seq),
            ("5' OH",          clean_oh5),
            ("3' OH",          clean_oh3),
            ("Forward primer", clean_fwd),
            ("Reverse primer", clean_rev),
        ):
            if not self._validate_iupac_chars(label_, val, status):
                return

        # ── Build updated entry ─────────────────────────────────────
        # Preserve any unrelated fields the user didn't see (schema-
        # version stamps, legacy qualifiers) so partial schemas survive.
        out = dict(self._part)
        out["name"]       = clean_name
        out["type"]       = clean_ptype
        out["position"]   = _sanitize_label(position, max_len=64)
        out["oh5"]        = clean_oh5
        out["oh3"]        = clean_oh3
        out["backbone"]   = _sanitize_label(backbone, max_len=120)
        out["marker"]     = _sanitize_label(marker,   max_len=120)
        out["sequence"]   = clean_seq
        out["fwd_primer"] = clean_fwd
        out["rev_primer"] = clean_rev
        out["grammar"]    = self._grammar_id
        # Re-derive Tms whenever the primer text changed (including
        # the empty → empty no-op via the equality short-circuit, so
        # we don't burn a primer3 call on a value the user didn't
        # touch). Round to 0.1 °C for stable JSON round-trips.
        if clean_fwd != self._part.get("fwd_primer", ""):
            tm = _primer_tm_safe(clean_fwd)
            out["fwd_tm"] = round(float(tm), 1) if tm is not None else 0.0
        if clean_rev != self._part.get("rev_primer", ""):
            tm = _primer_tm_safe(clean_rev)
            out["rev_tm"] = round(float(tm), 1) if tm is not None else 0.0
        # Re-derive simulator outputs when sequence, overhangs, OR
        # grammar changed (different grammar → different enzyme tail
        # in `primed_seq`). Empty sequence drops the derived fields so
        # a user who hand-clears the seq doesn't ship a stale primed
        # amplicon to Copy Primed.
        seq_or_oh_changed = (
            clean_seq != self._part.get("sequence", "")
            or clean_oh5 != self._part.get("oh5", "")
            or clean_oh3 != self._part.get("oh3", "")
            or self._grammar_id != (self._part.get("grammar") or "gb_l0")
        )
        if seq_or_oh_changed:
            if clean_seq:
                # Preserve any two-tier Golden-Braid nesting stored on the part
                # (external UPD entry overhangs wrapping the category pair) so an
                # edit to the body/overhangs doesn't silently collapse the primed
                # amplicon back to a one-tier layout that won't clone.
                _e5 = str(self._part.get("entry_oh5") or "")
                _e3 = str(self._part.get("entry_oh3") or "")
                _entry_overhangs = (_e5, _e3) if _e5 and _e3 else None
                out["primed_seq"] = _simulate_primed_amplicon(
                    clean_seq, clean_oh5, clean_oh3, grammar=self._grammar,
                    part_type=clean_ptype, entry_overhangs=_entry_overhangs,
                )
                out["cloned_seq"] = _simulate_cloned_plasmid(
                    clean_seq, clean_oh5, clean_oh3, clean_ptype,
                )
            else:
                out.pop("primed_seq", None)
                out.pop("cloned_seq", None)

        # No-op detection: if every comparable field round-trips
        # unchanged, dismiss without writing the file. Avoids burning
        # a JSON write + UI repopulate when the user clicked Edit
        # → Save with no actual modifications. Derived fields (Tm /
        # primed_seq / cloned_seq) follow from the exposed set, so
        # only comparing exposed fields is sufficient.
        exposed = ("name", "type", "position", "oh5", "oh3",
                   "backbone", "marker", "sequence",
                   "fwd_primer", "rev_primer", "grammar")
        original_grammar = self._part.get("grammar") or "gb_l0"
        baseline = {**self._part, "grammar": original_grammar}
        if all(out.get(k, "") == baseline.get(k, "") for k in exposed):
            self.dismiss(None)
            return
        self.dismiss({"idx": self._idx, "entry": out})


class LibrarySearchModal(_OneShotDismissScreen, ModalScreen):
    """Cross-collection plasmid search.

    Lists every plasmid across every collection on disk, fuzzy-
    filtered as the user types in the input box. The current
    library panel only filters within the active collection — this
    modal is the "where did I save that pUC19 variant" affordance.

    Dismiss payload:
      None                    — cancelled
      (collection, entry_id)  — user picked a row; caller switches
                                to that collection and loads the
                                plasmid via the existing
                                `_apply_record` flow.
    """

    # Sweep #26: search-query Input editing — block app-level Ctrl+Z.
    _blocks_undo: bool = True

    BINDINGS = [
        Binding("escape", "cancel",         "Cancel"),
        Binding("tab",    "app.focus_next", "Next",   show=False),
    ]

    # Debounce window for live filter. 150 ms feels instant on modern
    # terminals but coalesces a "type 5 chars in 200 ms" burst into
    # 1 search instead of 5 — matters when collections.json contains
    # thousands of plasmids and `_search_collections_library` is the
    # hot path.
    _LIVE_FILTER_DEBOUNCE_S = 0.15

    def __init__(self, *, initial_query: str = "") -> None:
        super().__init__()
        self._initial_query = initial_query
        self._matches: list[dict] = []
        self._filter_timer = None

    def compose(self) -> ComposeResult:
        with Vertical(id="libsearch-box"):
            yield Static(" Find plasmid (across all collections) ",
                         id="libsearch-title")
            yield Input(
                value=self._initial_query,
                placeholder="type to filter (fuzzy match)…",
                id="libsearch-input",
            )
            yield DataTable(id="libsearch-table",
                            cursor_type="row",
                            zebra_stripes=True)
            yield Static("", id="libsearch-status", markup=True)
            with Horizontal(id="libsearch-btns"):
                yield Button("Open", id="btn-libsearch-ok",
                             variant="primary")
                yield Button("Close", id="btn-libsearch-close")

    def on_mount(self) -> None:
        try:
            t = self.query_one("#libsearch-table", DataTable)
        except NoMatches:
            return
        t.add_columns("Plasmid", "Collection", "Status", "bp")
        # Initial population — user sees something immediately even
        # before they start typing. Cap at 200 (matches the helper
        # default; large libraries get a hint to narrow the query).
        self._refresh()
        try:
            self.query_one("#libsearch-input", Input).focus()
        except NoMatches:
            pass

    def _refresh(self) -> None:
        try:
            inp = self.query_one("#libsearch-input", Input)
            t   = self.query_one("#libsearch-table", DataTable)
        except NoMatches:
            return
        query = (inp.value or "").strip()
        # Strip the prefill text so a user who hasn't focused the
        # input yet sees all plasmids instead of "search for the
        # word `Search`".
        if query == _SearchInput.PREFILL:
            query = ""
        self._matches = _search_collections_library(query, limit=300)
        t.clear()
        for m in self._matches:
            status = m.get("status") or ""
            color  = _PLASMID_STATUS_COLORS.get(status)
            status_cell = (Text(status, style=f"{color} bold")
                           if status and color is not None
                           else Text("—", style="dim"))
            t.add_row(
                Text(m["name"], no_wrap=True, overflow="ellipsis"),
                Text(m["collection"], no_wrap=True, overflow="ellipsis"),
                status_cell,
                f"{m.get('size', 0):,}",
                key=f"{m['collection']}\x00{m['id']}",
            )
        try:
            self.query_one("#libsearch-status", Static).update(
                f"[dim]{len(self._matches)} match(es)"
                + (" — refine the query to see more"
                   if len(self._matches) >= 300 else "")
                + "[/dim]"
            )
        except NoMatches:
            pass

    @on(Input.Changed, "#libsearch-input")
    def _on_query_changed(self, _event: Input.Changed) -> None:
        # Live filter as the user types — debounced via `set_timer` so
        # a burst of keystrokes coalesces into a single
        # `_search_collections_library` call. Per-keystroke firing was
        # noticeably laggy on 5 k+ plasmid libraries (the matcher walks
        # every entry across every collection); the debounce keeps
        # typing snappy without sacrificing live-update behaviour.
        if self._filter_timer is not None:
            try:
                self._filter_timer.stop()
            except Exception:
                pass
        self._filter_timer = self.set_timer(
            self._LIVE_FILTER_DEBOUNCE_S, self._refresh,
        )

    @on(Input.Submitted, "#libsearch-input")
    def _on_query_submitted(self, _event: Input.Submitted) -> None:
        # Enter in the input == press Open if there's a match.
        if not self._matches:
            # Sweep #9 (2026-05-19): notify on empty result so the
            # user gets feedback instead of silent no-op.
            try:
                self.app.notify(
                    "No matches — try a different query.",
                    severity="warning", timeout=4,
                )
            except Exception:
                pass
            return
        try:
            t = self.query_one("#libsearch-table", DataTable)
        except NoMatches:
            return
        idx = t.cursor_row if t.cursor_row is not None else 0
        if 0 <= idx < len(self._matches):
            m = self._matches[idx]
            self.dismiss((m["collection"], m["id"]))

    @on(Button.Pressed, "#btn-libsearch-ok")
    def _ok_btn(self, _) -> None:
        self._on_query_submitted(None)  # type: ignore[arg-type]

    @on(DataTable.RowSelected, "#libsearch-table")
    def _row_selected(self, event: DataTable.RowSelected) -> None:
        idx = event.cursor_row
        if 0 <= idx < len(self._matches):
            m = self._matches[idx]
            self.dismiss((m["collection"], m["id"]))

    @on(Button.Pressed, "#btn-libsearch-close")
    def _close_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class HistoryViewerModal(_OneShotDismissScreen, ModalScreen):
    """Construction-history viewer for a library plasmid.

    Renders the parsed `<HistoryTree>` from
    `_parse_commercialsaas_history` into a Textual Tree widget — each
    history step is one tree node, expandable to show its parent
    fragments. Used by the LibraryPanel's `h` key binding; the App
    handler builds the tree object and pushes this modal.

    The display is read-only — modifying history is the
    responsibility of the construction action that produced it
    (Traditional cloning Save, future: Mutagenesis, Gibson, etc.).
    Closing returns to the library panel; the tree is rebuilt
    fresh each open so library edits in another window are
    reflected on the next view.
    """

    BINDINGS = [
        Binding("escape", "dismiss_history", "Close"),
        Binding("q",      "dismiss_history", "Close"),
    ]

    DEFAULT_CSS = """
    #hist-box {
        /* Was `height: 36;` (rigid) — switched to flex height with a
           cap so the dialog shrinks to content on a sparse history
           and still leaves Tree room to scroll on a deep one
           (2026-05-20 UX audit). */
        width: 110; max-width: 95%; height: 90%; max-height: 36; min-height: 18;
        background: $surface; border: solid $accent;
        padding: 1 2;
    }
    #hist-title {
        background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; text-align: center;
    }
    #hist-formerly { color: $text-muted; margin-bottom: 1; text-align: center; }
    #hist-warn { color: $warning; margin-bottom: 1; }
    #hist-proto-label { color: $accent; text-style: bold; }
    #hist-proto {
        height: auto; max-height: 9; margin-bottom: 1;
        border: solid $primary-darken-2; padding: 0 1; overflow-y: auto;
    }
    #hist-proto Static { padding: 0 1; }
    #hist-tree { height: 1fr; }
    #hist-detail {
        /* Was 8. The detail block now carries conditions, primer
           binding sites, end chemistry and the unmodelled-field
           backstop, so 8 rows meant scrolling past the Properties
           header to reach anything. Still scrollable — this just sets
           how much is visible without it. */
        height: 11; border: solid $primary-darken-2;
        padding: 0 1; margin-top: 1; overflow-y: auto;
    }
    #hist-detail Static { padding: 0 1; }
    #hist-btns { height: 3; align: right middle; margin-top: 1; }
    #hist-btns Button { min-width: 10; }
    """

    def __init__(self, title: str,
                  root_node: "_CommercialSaaSHistoryNode",
                  seq: str = "") -> None:
        super().__init__()
        self._title = title
        self._root_node = root_node
        # The entry's bases, used ONLY to check the ROOT node's claims
        # against the molecule — that's the one node we can identify with
        # certainty. Read-only: nothing here writes a sequence or a
        # history record.
        self._seq = str(seq or "").upper()
        # Map Textual tree-node-id → CommercialSaaS history node, populated
        # in `on_mount`. Lets the Selected handler look up the
        # backing history node without re-parsing the XML.
        self._node_by_id: "dict[int, _CommercialSaaSHistoryNode]" = {}

    def compose(self) -> ComposeResult:
        from rich.markup import escape as _esc
        with Vertical(id="hist-box"):
            yield Static(f" Construction history — {_esc(self._title)} ",
                          id="hist-title")
            # A plasmid renamed after it was built keeps the build-time
            # name on its history root, so the tree opens on a name the
            # user may not recognise. Say so instead of leaving them to
            # guess (177 of the reference library's 1,184 plasmids).
            was = _history_former_name(self._root_node, self._title)
            if was:
                yield Static(f"[dim]built as[/dim] {_esc(was)}",
                              id="hist-formerly", markup=True)
            # Flag up front when the record claims something the molecule
            # doesn't bear out. Read-only — nothing is corrected for them.
            warn = _history_consistency_summary(_history_node_warnings(
                self._root_node, self._seq,
                enzymes=_state._all_enzymes_hook()))
            if warn:
                yield Static(f" {warn}", id="hist-warn", markup=True)
            yield Static("Protocol", id="hist-proto-label")
            with VerticalScroll(id="hist-proto"):
                yield Static("", id="hist-proto-text", markup=True)
            yield Tree("History", id="hist-tree")
            with Vertical(id="hist-detail"):
                yield Static(
                    "[dim]Pick a node to see its details.[/]",
                    id="hist-detail-text", markup=True,
                )
            with Horizontal(id="hist-btns"):
                yield Button("Close", id="btn-hist-close")

    def on_mount(self) -> None:
        tree = self.query_one("#hist-tree", Tree)
        tree.show_root = False
        tree.guide_depth = 4
        # De-noised, iterative populate — shared with `HistoryScreen` so
        # both viewers render identically (only the product auto-expands;
        # repeated ancestors collapse to references). Replaces the old
        # recursive `_add`, which both auto-expanded everything AND could
        # blow the recursion limit on a deeply-nested hostile `.dna`.
        if _history_populate_tree(tree, self._root_node, self._node_by_id):
            _log.warning("history: modal tree render truncated by cap "
                         "for %r", self._title)
        try:
            proto = self.query_one("#hist-proto-text", Static)
            proto.update(_history_protocol_renderable(self._root_node))
        except NoMatches:
            pass
        # Auto-select the root so the detail pane has something
        # interesting on open.
        try:
            top = tree.root.children[0]
            tree.select_node(top)
        except IndexError:
            pass

    @staticmethod
    def _tree_label_for(node: "_CommercialSaaSHistoryNode") -> str:
        """One-line tree label — thin wrapper around the module-level
        `_history_tree_label` helper so the legacy modal and the
        fullscreen `HistoryScreen` render identical rows. Kept as a
        staticmethod for backwards compatibility with code that calls
        `HistoryViewerModal._tree_label_for(node)`."""
        return _history_tree_label(node)

    @on(Button.Pressed, "#btn-hist-close")
    def _on_close(self, _: Button.Pressed) -> None:
        self.dismiss(None)

    def action_dismiss_history(self) -> None:
        """Close the history viewer modal. Named distinctly from the
        Textual base `Screen.action_dismiss` (which is async + takes a
        result kwarg) so the override doesn't trip
        `reportIncompatibleMethodOverride` — the binding above routes
        to this method by name."""
        self.dismiss(None)

    @on(Tree.NodeSelected, "#hist-tree")
    def _on_node_selected(self, event) -> None:
        """When the user picks a node, render its full details into the
        bottom detail panel via the shared `_history_detail_lines` (same
        block the full-screen viewer shows; every dynamic field is
        Rich-escaped). Falls back gracefully if the lookup fails —
        defence-in-depth against the dict going stale."""
        hist = self._node_by_id.get(event.node.id)
        if hist is None:
            return
        try:
            detail = self.query_one("#hist-detail-text", Static)
        except NoMatches:
            return
        detail.update("\n".join(_history_detail_lines(
            hist,
            # Root only — an ancestor is just a name and a length here, and
            # resolving that to a library entry can pick the wrong molecule.
            self._seq if hist.element is self._root_node.element else "",
            enzymes=_state._all_enzymes_hook())))
