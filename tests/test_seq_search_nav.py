"""
test_seq_search_nav — Ctrl+F sequence search + alignment-click navigation.

Two user features landed 2026-06-01:

  1. **Alignment-lane click → seq-panel navigation.** Clicking an
     alignment row on the linear map no longer opens the AlignmentScreen
     detail modal; it jumps the per-base sequence panel to the clicked
     region (centered + highlighted) so misaligned / to-be-re-edited
     bases are one click away. Covered: the column→bp inverse
     (`_linear_x_to_bp` / `_linear_click_bp_span`), the `on_click`
     dispatch (posts `AlignmentLaneClicked`, does NOT push a screen), and
     the App handler that drives the seq panel.

  2. **Ctrl+F sequence search.** Find a DNA subsequence — fuzzy
     (allowable mismatches, Hamming), both strands, circular-aware — and
     step through hits with `n` / `N`. Covered: the pure matcher
     `_search_subsequence`, the shared `SequencePanel.focus_span`
     primitive, the find/next/prev actions (including wrap + the
     stale-sequence guard), and the keybinding remap (Ctrl+F = find,
     Alt+Shift+F = add-feature — NOT Alt+F, which Textual decodes as
     ctrl+right).
"""
from __future__ import annotations

import types

import pytest

import splicecraft as sc


TERMINAL_SIZE = (160, 48)


def _fake_click(x: int, y: int, *, shift: bool = False, ctrl: bool = False):
    """Minimal stand-in for a Textual Click event carrying just the
    attributes `PlasmidMap.on_click` reads."""
    return types.SimpleNamespace(
        x=x, y=y, shift=shift, ctrl=ctrl, button=1, chain=1,
        stop=lambda: None,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Pure matcher: _search_subsequence
# ═══════════════════════════════════════════════════════════════════════════════

class TestSearchSubsequence:

    def test_exact_forward_hits(self):
        hits = sc._search_subsequence(
            "AACCGGTTAACC", "AACC", max_mismatches=0,
            circular=False, both_strands=False,
        )
        assert [(h["start"], h["end"], h["strand"]) for h in hits] == [
            (0, 4, "+"), (8, 12, "+"),
        ]
        assert all(h["mismatches"] == 0 for h in hits)

    def test_reverse_strand_hit_uses_forward_coords(self):
        # RC(AACC) == GGTT, which sits at [4, 8) on the top strand.
        hits = sc._search_subsequence(
            "AACCGGTTAACC", "AACC", max_mismatches=0, both_strands=True,
        )
        minus = [h for h in hits if h["strand"] == "-"]
        assert minus == [{"start": 4, "end": 8, "strand": "-",
                          "mismatches": 0}]

    def test_both_strands_false_omits_reverse(self):
        hits = sc._search_subsequence(
            "AACCGGTTAACC", "AACC", both_strands=False,
        )
        assert all(h["strand"] == "+" for h in hits)

    def test_mismatch_tolerance_finds_near_matches(self):
        hits = sc._search_subsequence(
            "AAACAAAA", "AAAA", max_mismatches=1, both_strands=False,
        )
        # Windows 0..3 each carry one mismatch; window 4 is exact.
        starts = {(h["start"], h["mismatches"]) for h in hits}
        assert (4, 0) in starts
        assert (0, 1) in starts and (1, 1) in starts

    def test_k_zero_is_exact_only(self):
        hits = sc._search_subsequence(
            "AAACAAAA", "AAAA", max_mismatches=0, both_strands=False,
        )
        assert [h["start"] for h in hits] == [4]

    def test_iupac_n_in_query_matches_any(self):
        hits = sc._search_subsequence(
            "ATGCATGC", "ANGC", max_mismatches=0, both_strands=False,
        )
        assert [h["start"] for h in hits] == [0, 4]
        assert all(h["mismatches"] == 0 for h in hits)

    def test_iupac_in_template_matches(self):
        # Template R = {A,G}; query A is compatible at that position.
        hits = sc._search_subsequence(
            "RTGC", "ATGC", max_mismatches=0, both_strands=False,
        )
        assert [h["start"] for h in hits] == [0]

    def test_circular_wrap_hit(self):
        hits = sc._search_subsequence(
            "TTAAACGG", "GGTT", max_mismatches=0,
            circular=True, both_strands=False,
        )
        # GG at [6,8) + TT at [0,2) → wraps; start 6, end 10 (> n).
        assert hits == [{"start": 6, "end": 10, "strand": "+",
                         "mismatches": 0}]

    def test_linear_does_not_wrap(self):
        hits = sc._search_subsequence(
            "TTAAACGG", "GGTT", max_mismatches=0,
            circular=False, both_strands=False,
        )
        assert hits == []

    def test_palindrome_dedups_to_single_plus(self):
        # GAATTC is its own reverse-complement: the '-' hit collapses
        # into the '+' so find-next never visits the bp twice.
        hits = sc._search_subsequence(
            "GGGAATTCGG", "GAATTC", both_strands=True,
        )
        assert hits == [{"start": 2, "end": 8, "strand": "+",
                         "mismatches": 0}]

    def test_empty_query(self):
        assert sc._search_subsequence("ACGT", "") == []

    def test_query_longer_than_seq(self):
        assert sc._search_subsequence("ACGT", "ACGTACGT") == []

    def test_negative_mismatches(self):
        assert sc._search_subsequence("ACGT", "ACGT", max_mismatches=-1) == []

    def test_foreign_char_raises(self):
        with pytest.raises(ValueError):
            sc._search_subsequence("ACGTACGT", "ACXT")

    def test_rna_query_mapped_to_dna(self):
        # U in the query is normalised to T before matching.
        hits = sc._search_subsequence(
            "ACGTACGT", "ACGU", both_strands=False,
        )
        assert [h["start"] for h in hits] == [0, 4]

    def test_whitespace_and_fasta_stripped_from_query(self):
        hits = sc._search_subsequence(
            "ACGTACGT", ">primer\nac gt\n", both_strands=False,
        )
        assert [h["start"] for h in hits] == [0, 4]

    def test_hits_sorted_by_start(self):
        hits = sc._search_subsequence(
            "ACGTACGTACGT", "ACGT", both_strands=True,
        )
        starts = [h["start"] for h in hits]
        assert starts == sorted(starts)

    def test_max_hits_cap(self):
        hits = sc._search_subsequence(
            "A" * 100, "A", max_mismatches=0, both_strands=False,
            max_hits=10,
        )
        assert len(hits) == 10


# ═══════════════════════════════════════════════════════════════════════════════
# SequencePanel.focus_span — the shared jump-to-bp primitive
# ═══════════════════════════════════════════════════════════════════════════════

class TestFocusSpan:

    async def test_sets_cursor_and_selection(self, tiny_record,
                                             isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp.focus_span(10, 20, center=True, select=True)
            assert sp._cursor_pos == 10
            assert sp._user_sel == (10, 20)
            assert sp._sel_anchor == 10

    async def test_single_bp_span_is_cursor_only(self, tiny_record,
                                                 isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp.focus_span(7, 8, center=True, select=True)
            assert sp._cursor_pos == 7
            # A 1-bp span never masquerades as a region selection.
            assert sp._user_sel is None

    async def test_select_false_leaves_no_selection(self, tiny_record,
                                                    isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp.focus_span(10, 20, center=False, select=False)
            assert sp._cursor_pos == 10
            assert sp._user_sel is None

    async def test_clamps_out_of_range(self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            n = len(sp._seq)
            sp.focus_span(10_000, 20_000, center=True, select=True)
            # Cursor clamped into [0, n-1]; never out of bounds.
            assert 0 <= sp._cursor_pos <= n - 1

    async def test_empty_seq_is_noop(self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._seq = ""
            # Must not raise on an empty sequence.
            sp.focus_span(5, 10)

    async def test_wrap_hit_sets_two_piece_selection(self, tiny_record,
                                                     isolated_library):
        # [INV-96 / L2a] An origin-crossing search hit arrives with
        # end > n; focus_span must encode it as a WRAP selection
        # (start, head_end) with start > head_end, so the FULL match —
        # tail [start, n) + head [0, head_end) — highlights / copies /
        # annotates, not just the tail half (Ctrl+C + Alt+Shift+F + the
        # DNA-row `in_usr` wrap branch all read this encoding).
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            n = len(sp._seq)
            sp.focus_span(n - 3, n + 2, center=True, select=True)   # [n-3, n+2)
            assert sp._cursor_pos == n - 3
            assert sp._user_sel == (n - 3, 2), (
                "wrap hit must store (start, head_end) with start > head_end"
            )
            # The in_usr wrap branch must render without error.
            sp._refresh_view()
            await pilot.pause()

    async def test_hit_ending_exactly_at_origin_is_not_a_wrap(
            self, tiny_record, isolated_library):
        # Boundary: end == n is a normal tail selection, NOT a wrap
        # (head_end == 0 fails the `1 <= head_end < start` guard).
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            n = len(sp._seq)
            sp.focus_span(n - 5, n, center=True, select=True)
            assert sp._user_sel == (n - 5, n)


# ═══════════════════════════════════════════════════════════════════════════════
# Alignment-lane click → seq-panel navigation
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlignmentClickNavigation:

    @staticmethod
    def _linear_map(app):
        pm = app.query_one("#plasmid-map", sc.PlasmidMap)
        pm._map_mode = "linear"
        pm._linear_zoom = 1.0
        pm._linear_offset_bp = 0
        return pm

    async def test_click_bp_span_within_band(self, tiny_record,
                                             isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pm = self._linear_map(app)
            total = pm._total
            w = pm.size.width
            margin_l, margin_r = pm._LINEAR_MARGIN_L, pm._LINEAR_MARGIN_R
            usable_w = w - margin_l - margin_r
            # A click at the left edge of the band → bp 0; mid → ~total/2.
            assert pm._linear_x_to_bp(margin_l) == 0
            mid_bp = pm._linear_x_to_bp(margin_l + usable_w // 2)
            assert 0 < mid_bp < total
            lo, hi = pm._linear_click_bp_span(margin_l + usable_w // 2)
            assert 0 <= lo < hi <= total

    async def test_click_bp_span_rejects_margin_and_circular(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pm = self._linear_map(app)
            # Left margin (x < margin_l) → no bp.
            assert pm._linear_click_bp_span(0) == (-1, -1)
            # Circular mode paints no band → no navigation.
            pm._map_mode = "circular"
            assert pm._linear_click_bp_span(pm._LINEAR_MARGIN_L + 3) == (-1, -1)

    async def test_handler_moves_seq_panel(self, tiny_record,
                                           isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            msg = sc.PlasmidMap.AlignmentLaneClicked(0, 30, 42)
            app._map_alignment_lane_clicked(msg)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            assert sp._cursor_pos == 30
            assert sp._user_sel == (30, 42)

    async def test_handler_survives_missing_seq_panel(self, tiny_record,
                                                      isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            msg = sc.PlasmidMap.AlignmentLaneClicked(0, 30, 42)
            # Monkeypatch query_one to raise NoMatches for the seq panel
            # — the handler must swallow it, not crash the click.
            from textual.css.query import NoMatches
            orig = app.query_one

            def _boom(sel, *a, **k):
                if sel == "#seq-panel":
                    raise NoMatches("seq panel gone")
                return orig(sel, *a, **k)

            app.query_one = _boom  # type: ignore[assignment]
            try:
                app._map_alignment_lane_clicked(msg)  # must not raise
            finally:
                app.query_one = orig  # type: ignore[assignment]

    async def test_on_click_navigates_not_drills_in(self, tiny_record,
                                                    isolated_library):
        """The core behaviour change end-to-end: an alignment-lane click
        jumps the seq panel to the clicked bp and does NOT push the
        AlignmentScreen detail modal. Asserted via real machinery (the
        posted message is processed by the App handler) rather than by
        monkeypatching core Textual methods — patching `post_message` /
        `push_screen` on live objects deadlocks run_test teardown."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            pm = self._linear_map(app)
            # Seed an alignment + a click bbox by hand so the hit-test
            # resolves without depending on exact render geometry.
            pm._alignments = [{
                "name": "read1", "query_label": "q", "target_label": "t",
                "target_record": tiny_record,
                "result": {"aligned_q": "A", "aligned_t": "A"},
                "aligned_q": "A", "aligned_t": "A",
                "segments": [(0, 1, "match")], "t_lo": 0, "t_hi": 1,
                "axis": "target", "letters": None,
            }]
            row = max(4, pm.size.height // 2) + 2
            x0 = pm._LINEAR_MARGIN_L + 2
            pm._align_bboxes = [(x0, x0 + 6, row, 0)]
            expected_bp = pm._linear_x_to_bp(x0 + 1)

            n_screens = len(app.screen_stack)
            pm.on_click(_fake_click(x0 + 1, row))
            await pilot.pause()        # let the posted message be handled

            sp = app.query_one("#seq-panel", sc.SequencePanel)
            assert sp._cursor_pos == expected_bp      # jumped to the click
            assert len(app.screen_stack) == n_screens  # no drill-in modal


# ═══════════════════════════════════════════════════════════════════════════════
# Ctrl+F find flow: search → next/prev → wrap → stale guard
# ═══════════════════════════════════════════════════════════════════════════════

class TestFindSequenceFlow:

    @staticmethod
    def _result(hits, **kw):
        base = {"query": "AAA", "mismatches": 0, "both_strands": True,
                "hits": hits}
        base.update(kw)
        return base

    async def test_on_search_stores_and_jumps(self, tiny_record,
                                              isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 0
            hits = [
                {"start": 5, "end": 8, "strand": "+", "mismatches": 0},
                {"start": 40, "end": 43, "strand": "+", "mismatches": 0},
            ]
            app._on_seq_search(self._result(hits))
            assert app._seq_search is not None
            assert app._seq_search["idx"] == 0
            assert sp._cursor_pos == 5

    async def test_next_and_prev_wrap(self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            hits = [
                {"start": 5, "end": 8, "strand": "+", "mismatches": 0},
                {"start": 40, "end": 43, "strand": "+", "mismatches": 0},
            ]
            app._on_seq_search(self._result(hits))
            assert app._seq_search["idx"] == 0
            app.action_find_next()
            assert app._seq_search["idx"] == 1
            assert sp._cursor_pos == 40
            app.action_find_next()           # wraps back to 0
            assert app._seq_search["idx"] == 0
            assert sp._cursor_pos == 5
            app.action_find_prev()           # wraps to last
            assert app._seq_search["idx"] == 1
            assert sp._cursor_pos == 40

    async def test_next_with_no_search_is_safe(self, tiny_record,
                                               isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app._seq_search = None
            app.action_find_next()          # notify, no crash
            app.action_find_prev()
            assert app._seq_search is None

    async def test_stale_sequence_invalidates_search(self, tiny_record,
                                                     isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            hits = [{"start": 5, "end": 8, "strand": "+", "mismatches": 0}]
            app._on_seq_search(self._result(hits))
            assert app._seq_search is not None
            # Simulate an edit: the recorded seq_len no longer matches.
            app._seq_search["seq_len"] = 999999
            app._goto_search_hit()
            assert app._seq_search is None   # dropped rather than mis-jump

    async def test_start_idx_picks_hit_after_cursor(self, tiny_record,
                                                    isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 30           # between the two hits
            hits = [
                {"start": 5, "end": 8, "strand": "+", "mismatches": 0},
                {"start": 40, "end": 43, "strand": "+", "mismatches": 0},
            ]
            app._on_seq_search(self._result(hits))
            # Search continues from where the user was looking → hit @40.
            assert app._seq_search["idx"] == 1
            assert sp._cursor_pos == 40

    async def test_n_and_N_keystrokes_fire(self, tiny_record,
                                           isolated_library):
        """The real keystrokes — not just the binding registration —
        drive next / previous. Catches the Shift+n → "N" delivery quirk
        that a binding-map check alone would miss."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            # Focus the map so the plain letters bubble to the App
            # bindings (an Input would otherwise swallow them as typing).
            app.set_focus(app.query_one("#plasmid-map", sc.PlasmidMap))
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            hits = [
                {"start": 5, "end": 8, "strand": "+", "mismatches": 0},
                {"start": 40, "end": 43, "strand": "+", "mismatches": 0},
            ]
            app._on_seq_search(self._result(hits))
            assert app._seq_search["idx"] == 0
            await pilot.press("n")
            await pilot.pause()
            assert app._seq_search["idx"] == 1
            assert sp._cursor_pos == 40
            await pilot.press("N")          # standard-terminal Shift+n
            await pilot.pause()
            assert app._seq_search["idx"] == 0
            assert sp._cursor_pos == 5

    async def test_record_swap_clears_search(self, tiny_record,
                                             isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            hits = [{"start": 5, "end": 8, "strand": "+", "mismatches": 0}]
            app._on_seq_search(self._result(hits))
            assert app._seq_search is not None
            # Loading a fresh plasmid must invalidate the stale search.
            app._apply_record(tiny_record)
            assert app._seq_search is None

    async def test_ctrl_f_does_not_stack_modals(self, tiny_record,
                                                isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.action_find_sequence()
            await pilot.pause()
            depth = len(app.screen_stack)
            assert isinstance(app.screen, sc.SeqSearchModal)
            # Re-pressing Ctrl+F while the box is open must NOT stack.
            app.action_find_sequence()
            await pilot.pause()
            assert len(app.screen_stack) == depth


# ═══════════════════════════════════════════════════════════════════════════════
# Keybinding remap: Ctrl+F = find, Alt+Shift+F = add-feature (not Alt+F)
# ═══════════════════════════════════════════════════════════════════════════════

class TestAddFeatureKeystroke:
    """The actual `alt+shift+f` keystroke opens the Add Feature modal —
    the fix for the user bug where Alt+F just nudged the selection
    (Textual decodes ESC-f as ctrl+right, so the old alt+f binding could
    never fire). `pilot.press` uses Textual's key name, so this proves
    the BINDING fires; the ESC-f → ctrl+right and ESC-F → alt+shift+f
    decoding is verified separately at the parser level in the work log."""

    async def test_alt_shift_f_opens_add_feature_modal(self, tiny_record,
                                                       isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.set_focus(app.query_one("#plasmid-map", sc.PlasmidMap))
            n0 = len(app.screen_stack)
            await pilot.press("alt+shift+f")
            await pilot.pause()
            assert len(app.screen_stack) > n0
            assert type(app.screen).__name__ == "AddFeatureModal"


class TestKeybindingRemap:

    @staticmethod
    def _binding_map(screen_cls):
        """{(key): action} for a screen/app's BINDINGS (Binding objects
        or (key, action, desc) tuples)."""
        out = {}
        for b in getattr(screen_cls, "BINDINGS", []):
            key = getattr(b, "key", None)
            action = getattr(b, "action", None)
            if key is None and isinstance(b, (tuple, list)):
                key, action = b[0], b[1]
            out.setdefault(key, action)
        return out

    def test_app_ctrl_f_is_find_sequence(self):
        m = self._binding_map(sc.PlasmidApp)
        assert m.get("ctrl+f") == "find_sequence"

    def test_app_add_feature_on_alt_shift_f_not_alt_f(self):
        # add_feature lives on Alt+Shift+F. It must NOT be on Alt+F:
        # Textual decodes ESC-f (Alt+F) as `ctrl+right`, so an alt+f
        # binding can never fire and would just nudge the cursor.
        m = self._binding_map(sc.PlasmidApp)
        assert m.get("alt+shift+f") == "add_feature"
        assert m.get("alt+f") != "add_feature"

    def test_app_n_keys_step_search(self):
        m = self._binding_map(sc.PlasmidApp)
        assert m.get("n") == "find_next"
        # find_prev must be bound on BOTH the bare uppercase "N" (what a
        # standard terminal sends for Shift+n) and "shift+n" (kitty /
        # enhanced protocols) so it fires everywhere.
        assert m.get("N") == "find_prev"
        assert m.get("shift+n") == "find_prev"

    def test_nothing_bound_to_undeliverable_alt_f_or_alt_b(self):
        # Alt+F (ESC-f → ctrl+right) and Alt+B (ESC-b → ctrl+left) can't
        # be delivered, so the App must bind NOTHING to them — no
        # add_feature / File menu on alt+f, no BLAST menu on alt+b
        # (the dead alt+b binding was removed 2026-06-01).
        m = self._binding_map(sc.PlasmidApp)
        assert m.get("alt+f") != "add_feature"
        assert "open_named_menu" not in str(m.get("alt+f") or "")
        assert m.get("alt+b") is None

    def test_synthesis_screen_add_feature_on_alt_shift_f(self):
        m = self._binding_map(sc.SynthesisScreen)
        assert m.get("alt+shift+f") == "add_feature"
        # Neither Ctrl+F nor the undeliverable Alt+F adds a feature here.
        assert m.get("ctrl+f") != "add_feature"
        assert m.get("alt+f") != "add_feature"

    def test_actions_exist_on_app(self):
        for name in ("action_find_sequence", "action_find_next",
                     "action_find_prev", "action_add_feature"):
            assert callable(getattr(sc.PlasmidApp, name, None)), name
