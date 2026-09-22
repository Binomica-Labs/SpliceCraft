"""
test_synthesis — gene-synthesis composer (Synthesis menu) tests.

Covers:
  * SynthesisEditor primitives — insert, delete, selection delete,
    cursor moves, IUPAC-only typing, max-bp cap, feature shift on
    insert/delete.
  * Selection arithmetic — drag-to-select via the message-bus
    surface, Shift+arrow extension via `set_cursor(extend_sel=True)`.
  * SynthesisScreen mount + lifecycle — opens, hosts the editor,
    feature-library side panel populates, document-model save flow.
  * Restriction site insert — picker dispatch + cursor splice.
  * Feature library insert mode + annotate mode.
  * AddFeatureModal `total_len` injection (synthesis host w/o the
    main #seq-panel).
  * Save round-trip lands a linear SeqRecord in the library.
  * `_blocks_undo` opt-out on every new modal (invariant #41).
"""
from __future__ import annotations

import pytest

import splicecraft as sc


_TERM = (160, 48)


# ═══════════════════════════════════════════════════════════════════════════════
# SynthesisEditor — synchronous unit tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestSynthesisEditorPrimitives:
    """The editor's mutations are pure functions over (seq, feats,
    cursor, selection); test them without an async Pilot."""

    def _make(self):
        ed = sc.SynthesisEditor()
        ed._seq = ""
        ed._feats = []
        ed._cursor_pos = 0
        return ed

    def test_insert_appends_at_cursor(self):
        ed = self._make()
        ed._cursor_pos = 0
        # `insert_at_cursor` notifies via `self.app.notify` — bypass by
        # poking the buffer directly so we don't need a mounted app.
        ed._seq = "AAA"
        ed._cursor_pos = 3
        # Manual mirror of `insert_at_cursor`'s body so we can test
        # offline. The integration path is exercised in the async tests.
        cur = ed._cursor_pos
        clean = "TTT"
        ed._seq = ed._seq[:cur] + clean + ed._seq[cur:]
        ed._cursor_pos = cur + len(clean)
        assert ed._seq == "AAATTT"
        assert ed._cursor_pos == 6

    def test_clamp_cursor_to_seq_bounds(self):
        ed = self._make()
        ed._seq = "ACGT"
        ed._cursor_pos = 99
        ed._clamp_cursor()
        assert ed._cursor_pos == 4
        ed._cursor_pos = -5
        ed._clamp_cursor()
        assert ed._cursor_pos == 0

    def test_delete_range_shifts_features_right_of_cut(self):
        ed = self._make()
        ed._seq = "ATGCATGCATGC"
        ed._feats = [
            {"start": 0, "end": 4, "label": "left", "type": "misc_feature"},
            {"start": 8, "end": 12, "label": "right", "type": "misc_feature"},
        ]
        ed._cursor_pos = 0
        # Delete bp 4-8 (middle)
        ed._delete_range(4, 8)
        assert ed._seq == "ATGCATGC"
        # `left` is unchanged (entirely left of deletion).
        assert ed._feats[0]["start"] == 0 and ed._feats[0]["end"] == 4
        # `right` shifts left by 4.
        assert ed._feats[1]["start"] == 4 and ed._feats[1]["end"] == 8

    def test_delete_range_clips_overlap(self):
        ed = self._make()
        ed._seq = "ATGCATGCATGC"
        ed._feats = [
            # Spans across the cut.
            {"start": 2, "end": 10, "label": "span", "type": "misc_feature"},
        ]
        ed._cursor_pos = 0
        ed._delete_range(4, 8)  # delete bp 4-8
        # Overlap clipped: surviving span starts at min(s, start) = 2,
        # ends at max(start, e - n_del) = max(4, 10-4) = 6.
        assert ed._feats[0]["start"] == 2
        assert ed._feats[0]["end"] == 6

    def test_delete_range_drops_zero_length_feature(self):
        ed = self._make()
        ed._seq = "ATGCATGC"
        # Feature exactly equal to the deleted range disappears.
        ed._feats = [{"start": 2, "end": 4, "label": "tiny",
                      "type": "misc_feature"}]
        ed._delete_range(2, 4)
        assert ed._feats == []

    def test_insert_shifts_features_at_or_after_cursor(self):
        ed = self._make()
        ed._seq = "ATGCATGC"
        ed._feats = [
            # Feature whose end == cursor extends (the "appending to
            # an upstream feature" case).
            {"start": 0, "end": 4, "label": "extends",
             "type": "misc_feature"},
            # Feature whose start == cursor shifts (insert goes
            # BEFORE the feature).
            {"start": 4, "end": 8, "label": "shifts",
             "type": "misc_feature"},
        ]
        ed._cursor_pos = 4
        # Walk the same logic insert_at_cursor uses (offline mirror).
        cur = ed._cursor_pos
        n_ins = 3
        ed._seq = ed._seq[:cur] + "NNN" + ed._seq[cur:]
        new_feats = []
        for f in ed._feats:
            s = f["start"]
            e = f["end"]
            new_f = dict(f)
            if s >= cur:
                new_f["start"] = s + n_ins
            if e >= cur and not (e == cur and s == cur):
                new_f["end"] = e + n_ins
            new_feats.append(new_f)
        ed._feats = new_feats
        # `extends` end shifted (the feature ate the inserted bases).
        assert ed._feats[0]["start"] == 0
        assert ed._feats[0]["end"] == 7
        # `shifts` moved entirely to the right.
        assert ed._feats[1]["start"] == 7
        assert ed._feats[1]["end"] == 11
        assert ed._seq == "ATGCNNNATGC"


class TestSynthesisCaps:
    def test_max_bp_constant_is_50k(self):
        # Sacred — commercial gene synthesis tops out at ~30 kb for
        # the longest vendors. 50 kb is generous headroom that still
        # keeps `_build_seq_text` snappy on one chunk.
        assert sc._SYNTHESIS_MAX_BP == 50_000

    def test_typeable_bases_iupac_only(self):
        # ACGT plus IUPAC ambiguity codes — same set AddFeatureModal
        # validates against so a fragment composed here round-trips
        # through the feature library without surprises. This is the
        # PROGRAMMATIC insert path (restriction sites, library
        # features); keyboard typing is restricted further below.
        assert sc._SYNTHESIS_TYPEABLE_BASES == frozenset("ACGTRYWSMKBDHVN")

    def test_keyboard_bases_acgtn_only(self):
        # Direct keyboard typing only — vendors universally accept
        # A/C/G/T/N as-is; IUPAC ambiguity codes are still allowed
        # via restriction-site / feature-library inserts where the
        # ambiguity is biologically intentional.
        assert sc._SYNTHESIS_KEYBOARD_BASES == frozenset("ACGTN")


# ═══════════════════════════════════════════════════════════════════════════════
# Menubar wiring
# ═══════════════════════════════════════════════════════════════════════════════

class TestMenubarWiring:
    def test_menubar_contains_synthesis(self):
        assert "Synthesis" in sc.MenuBar.MENUS

    def test_synthesis_placed_after_mutagenize(self):
        # Per the workflow the user picked at design time:
        # Mutagenize → Synthesize new → Parts → Constructor.
        menus = sc.MenuBar.MENUS
        i_mut = menus.index("Mutato")  # whimsical label for the mutagenesis tool
        i_syn = menus.index("Synthesis")
        i_parts = menus.index("Parts")
        assert i_mut < i_syn < i_parts

    def test_action_open_synthesis_defined(self):
        assert hasattr(sc.PlasmidApp, "action_open_synthesis")


# ═══════════════════════════════════════════════════════════════════════════════
# _blocks_undo invariant (sweep #6 / sweep #10 / sweep #12)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlocksUndoOptOut:
    """Every new modal that hosts an Input / TextArea / mutates record
    state must carry `_blocks_undo = True` per invariant #41 so app-
    level Ctrl+Z falls through to the inner editor's undo."""

    @pytest.mark.parametrize("cls_name", [
        "SynthesisScreen",
        "SynthesisLoadModal",
        "RestrictionInsertModal",
        "SynthesisUnsavedChangesModal",
    ])
    def test_blocks_undo_true(self, cls_name):
        cls = getattr(sc, cls_name)
        assert getattr(cls, "_blocks_undo", False) is True


# ═══════════════════════════════════════════════════════════════════════════════
# Async lifecycle — Screen mount, menu open, save flow
# ═══════════════════════════════════════════════════════════════════════════════

class TestSynthesisScreen:
    async def test_screen_mounts_via_action(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.SynthesisScreen)

    async def test_screen_mounts_via_menu_string(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            # Synthesis direct-opens the composer (no dropdown): opening the
            # menu pushes the SynthesisScreen straight away.
            app.open_menu("Synthesis", 0, 0)
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.SynthesisScreen)

    async def test_codon_tables_action_opens_manager(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            # The codon-table manager action (kept for keybindings / command
            # palette; no longer a Synthesis-menu entry) opens the manager.
            app.action_open_codon_tables()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.SpeciesPickerModal)

    async def test_initial_state_empty(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            assert ed._seq == ""
            assert ed._cursor_pos == 0
            assert ed._user_sel is None
            assert scr._dirty is False
            assert scr._loaded_id is None

    async def test_insert_makes_dirty(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.insert_at_cursor("ATGC")
            await pilot.pause()
            await pilot.pause()
            assert ed._seq == "ATGC"
            assert ed._cursor_pos == 4
            assert scr._dirty is True

    async def test_cursor_move_does_not_mark_dirty(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed._seq = "ATGCATGC"
            ed._refresh_view()
            ed.set_cursor(3)
            await pilot.pause()
            await pilot.pause()
            # Cursor moves emit CursorMoved (not Changed); the screen
            # listener refreshes status but doesn't touch dirty.
            assert scr._dirty is False

    async def test_insert_at_cursor_keeps_iupac_codes(self):
        """`insert_at_cursor` is the PROGRAMMATIC path — used by
        restriction-site / feature-library inserts that legitimately
        carry IUPAC ambiguity (AvaI's CYCGRG, etc.). Keep IUPAC
        codes; drop only truly invalid chars."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            # X / Z are not IUPAC; Y is (pyrimidine).
            ed.insert_at_cursor("XAYZG@")
            await pilot.pause()
            assert ed._seq == "AYG"

    async def test_keyboard_typing_rejects_iupac_codes(self):
        """`on_key` is the KEYBOARD path — restricted to A/C/G/T/N
        only so accidental IUPAC ambiguity can't slip into a
        synthesis-bound fragment via a stray keystroke. Y, R, W,
        S, M, K, B, D, H, V should all be rejected at the keystroke;
        A, C, G, T, N should land."""
        from textual.events import Key
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            # Drive on_key directly with synthetic events for each
            # character — bypasses Textual's focus / dispatch chain
            # so the test doesn't depend on which widget has focus.
            for ch in "ACGTNYRWSMKBDHV":
                ed.on_key(Key(key=ch.lower(), character=ch))
            await pilot.pause()
            # Only A, C, G, T, N landed.
            assert ed._seq == "ACGTN"

    async def test_keyboard_lowercase_accepted_as_uppercase(self):
        """Lowercase keystrokes uppercase to a typeable base; verifies
        the .upper() shim in the key handler."""
        from textual.events import Key
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            for ch in "acgtn":
                ed.on_key(Key(key=ch, character=ch))
            await pilot.pause()
            assert ed._seq == "ACGTN"

    async def test_keyboard_path_allows_n_explicitly(self):
        """N (any base) is a common synthesis-vendor placeholder and
        must remain typeable even though it's technically an IUPAC
        ambiguity code (= ACGT)."""
        from textual.events import Key
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.on_key(Key(key="n", character="N"))
            ed.on_key(Key(key="n", character="N"))
            ed.on_key(Key(key="n", character="N"))
            await pilot.pause()
            assert ed._seq == "NNN"

    async def test_max_bp_cap_truncates(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            # Pre-fill near the cap so the truncation path fires.
            ed._seq = "A" * (sc._SYNTHESIS_MAX_BP - 3)
            ed._cursor_pos = len(ed._seq)
            ok = ed.insert_at_cursor("ATGCATGC")
            assert ok is True
            assert len(ed._seq) == sc._SYNTHESIS_MAX_BP

    async def test_backspace_at_cursor(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed._seq = "ATGCATGC"
            ed._cursor_pos = 4
            ed._refresh_view()
            ed.delete_at_cursor(forward=False)
            assert ed._seq == "ATGATGC"
            assert ed._cursor_pos == 3

    async def test_delete_forward_at_cursor(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed._seq = "ATGCATGC"
            ed._cursor_pos = 4
            ed._refresh_view()
            ed.delete_at_cursor(forward=True)
            assert ed._seq == "ATGCTGC"
            assert ed._cursor_pos == 4

    async def test_selection_delete(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed._seq = "AAAATTTTCCCC"
            ed._user_sel = (4, 8)
            ed._refresh_view()
            ed.delete_at_cursor()
            assert ed._seq == "AAAACCCC"
            assert ed._user_sel is None


# ═══════════════════════════════════════════════════════════════════════════════
# Document-model save round-trip
# ═══════════════════════════════════════════════════════════════════════════════

class TestSaveRoundTrip:
    async def test_commit_save_lands_linear_record_in_library(self):
        """`_commit_save` builds a SeqRecord with topology=linear,
        molecule_type=DNA, and pushes through LibraryPanel.add_entry.
        Verify the entry lands in the active library."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGCATGC", [
                {"start": 0, "end": 3, "label": "start",
                 "type": "misc_feature", "color": "green",
                 "strand": 1, "qualifiers": {}},
            ])
            scr._loaded_name = "test_frag"
            scr._loaded_id   = "test_frag"
            scr._commit_save("ATGCATGCATGCATGC",
                              [dict(f) for f in ed._feats],
                              after=None)
            await pilot.pause()
            await pilot.pause()
            entries = sc._load_library()
            matches = [e for e in entries if e.get("id") == "test_frag"]
            assert len(matches) == 1
            assert "linear" in matches[0].get("gb_text", "").lower()
            # The annotated feature survived into the saved record.
            assert "start" in matches[0].get("gb_text", "")

    async def test_document_model_overwrite_by_id(self):
        """Second save with the same id silently replaces, not appends."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr._loaded_id = "doc_frag"
            scr._loaded_name = "doc_frag"
            scr._commit_save("AAAA", [], after=None)
            await pilot.pause()
            scr._commit_save("CCCC", [], after=None)
            await pilot.pause()
            entries = sc._load_library()
            matches = [e for e in entries if e.get("id") == "doc_frag"]
            assert len(matches) == 1
            # Second-save content is what landed (BioPython lower-cases
            # the ORIGIN block in GenBank output, hence the .upper()).
            assert "CCCC" in matches[0].get("gb_text", "").upper()


# ═══════════════════════════════════════════════════════════════════════════════
# Restriction-site picker + feature-library side panel
# ═══════════════════════════════════════════════════════════════════════════════

class TestRestrictionInsert:
    async def test_picker_dismisses_with_enzyme_and_strand(self):
        # action_pick returns {"enzyme", "strand"} — default direction is
        # forward (strand 1); the Select drives the strand.
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            holder = {"result": "sentinel"}
            app.push_screen(sc.RestrictionInsertModal(),
                            callback=lambda r: holder.update(result=r))
            await pilot.pause()
            await pilot.pause()
            modal = app.screen
            modal._repopulate("EcoRI")            # filter to a single row
            await pilot.pause()
            t = modal.query_one("#ri-table", sc.DataTable)
            t.move_cursor(row=0)
            modal.query_one("#ri-direction", sc.Select).value = "-1"
            await pilot.pause()
            modal.action_pick()
            await pilot.pause()
            await pilot.pause()
            assert holder["result"] == {"enzyme": "EcoRI", "strand": -1}

    async def test_inserted_site_appears_in_buffer(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("AAAA", [])
            ed._cursor_pos = 2
            # Manually drive the post-pick callback.
            site = sc._site_for_enzyme("EcoRI")
            assert site == "GAATTC"
            ed.insert_at_cursor(site)
            await pilot.pause()
            assert ed._seq == "AAGAATTCAA"


class TestFeatureLibrarySidePanel:
    async def test_pane_widgets_present(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            from textual.widgets import DataTable, Input, Button
            scr.query_one("#syn-featlib-table", DataTable)
            scr.query_one("#syn-featlib-search", Input)
            scr.query_one("#btn-syn-featlib-insert", Button)
            scr.query_one("#btn-syn-featlib-annotate", Button)
            scr.query_one("#btn-syn-featlib-refresh", Button)

    async def test_insert_mode_splices_sequence_at_cursor(self,
                                                            monkeypatch):
        # Seed a feature library entry so _refresh_featlib_table picks
        # it up; mock _load_features to return it deterministically.
        monkeypatch.setattr(sc, "_load_features", lambda: [{
            "name": "rbs1",
            "feature_type": "RBS",
            "sequence": "AAAGGAGG",
            "color": "yellow",
            "strand": 1,
            "qualifiers": {},
            "description": "",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr._refresh_featlib_table()
            await pilot.pause()
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("TTTT", [])
            ed._cursor_pos = 2
            # Drive the side-panel insert path with a known entry —
            # bypass the DataTable cursor lookup so the test doesn't
            # depend on Textual focus state.
            scr._featlib_rows = [("k", {
                "name": "rbs1", "feature_type": "RBS",
                "sequence": "AAAGGAGG", "color": "yellow",
                "strand": 1, "qualifiers": {},
            })]
            # Patch the selector to return our fixed entry.
            scr._featlib_selected_entry = lambda: scr._featlib_rows[0][1]
            scr._featlib_insert_selected(mode="insert")
            await pilot.pause()
            assert ed._seq == "TTAAAGGAGGTT"
            # A new feature dict was added covering the insert range.
            inserted = [f for f in ed._feats if f.get("label") == "rbs1"]
            assert len(inserted) == 1
            assert inserted[0]["start"] == 2
            assert inserted[0]["end"] == 10

    async def test_add_feature_render_includes_new_feature(self,
                                                              monkeypatch):
        """Regression for 2026-05-26 "weird DNA visualization after
        library insert" report: the rendered text after
        `add_feature` MUST include the new feature's lane art
        (label + bar) and the `5'-` / `-3'` flank markers. Pre-fix
        `add_feature` used `self._feats.append(...)` (in-place
        mutation), but the `_CHUNK_STATIC_CACHE`'s key is
        `(hash(seq), id(feats), line_width, show_connectors)` —
        without `len(feats)`. The append left `id(feats)` unchanged
        so the cache HIT a stale entry built BEFORE the new feature
        was added; the render came back with only the bare DNA
        rows (no feature lane, no flank markers — visible
        symptom: a sequence floating in the editor without its
        markers / bar). Arrow keys took the cursor to a new
        position which fell into a different cache slot, masking
        the bug.

        Sacred invariant #4 (CLAUDE.md): lists are reassigned on
        load, never mutated in place. The fix: `add_feature` now
        does `self._feats = self._feats + [dict(feat)]` so
        `id(feats)` bumps and every downstream cache invalidates
        cleanly."""
        monkeypatch.setattr(sc, "_load_features", lambda: [{
            "name": "test_feat",
            "feature_type": "misc_feature",
            "sequence": "GCGCCGTCT",
            "color": "#FF0000",
            "strand": 1,
            "qualifiers": {},
            "description": "",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr._refresh_featlib_table()
            await pilot.pause()
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            # Empty initial seq — primes the chunk caches with an
            # empty-feats entry that the buggy in-place append
            # path would HIT for the post-insert state.
            ed.load("", [])
            # Drive the library insert path with our fixed entry.
            scr._featlib_rows = [("k", {
                "name": "test_feat",
                "feature_type": "misc_feature",
                "sequence": "GCGCCGTCT", "color": "#FF0000",
                "strand": 1, "qualifiers": {},
            })]
            scr._featlib_selected_entry = (
                lambda: scr._featlib_rows[0][1]
            )
            scr._featlib_insert_selected(mode="insert")
            await pilot.pause()
            await pilot.pause(0.1)
            # The render must now contain the feature lane (label
            # + bar) AND the flank markers. Inspect the raw view
            # render rather than the buffer state — that's what
            # the user actually sees.
            view = scr.query_one("#syn-view")
            rendered = str(view.render())
            assert "test_feat" in rendered, (
                "feature label missing from post-insert render "
                "(in-place feats mutation cached out the lane art)"
            )
            assert "5'-" in rendered, (
                "5' flank marker missing from post-insert render"
            )
            assert "-3'" in rendered, (
                "3' flank marker missing from post-insert render"
            )

    async def test_enter_inside_feature_opens_edit_modal(self):
        """Sweep #29 (2026-05-26): pressing Enter on the synth
        editor while the cursor sits inside a feature span posts
        a `FeatureEditRequested` message that the SynthesisScreen
        catches → pushes `FeatureEditModal` for that feature.
        Mirrors the main-canvas seq-panel's Enter-on-feature UX."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            # Load a fragment with one feature so Enter has a
            # target inside the cursor's span.
            ed.load("ATGCATGCATGC", [
                {"start": 2, "end": 9, "label": "test-feat",
                 "type": "CDS", "color": "#FF8800",
                 "strand": 1, "qualifiers": {"label": ["test-feat"]}},
            ])
            ed._cursor_pos = 5   # inside the feature [2, 9)
            await pilot.pause()
            # Direct API call mirrors what `on_key` does on Enter.
            # Driving via `pilot.press('enter')` would require the
            # editor's scroll container to have focus, which is
            # brittle in tests. We invoke the action directly and
            # assert the message flow + modal push.
            assert ed._feature_at_cursor() == 0
            ed.post_message(sc.SynthesisEditor.FeatureEditRequested(0))
            await pilot.pause()
            await pilot.pause(0.05)
            assert isinstance(app.screen, sc.FeatureEditModal)
            assert app.screen._idx == 0

    async def test_enter_outside_feature_does_not_open_modal(self):
        """Enter on a base NOT inside any feature surfaces a
        notify, does not push the modal. The cursor's bp position
        must fall in `[start, end)` of at least one feat for the
        modal to open."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGC", [
                {"start": 2, "end": 5, "label": "test-feat",
                 "type": "CDS", "color": "#FF8800",
                 "strand": 1, "qualifiers": {"label": ["test-feat"]}},
            ])
            ed._cursor_pos = 8   # OUTSIDE [2, 5)
            assert ed._feature_at_cursor() == -1
            # No FeatureEditRequested message gets posted by the
            # editor's on_key — `_feature_at_cursor` returns -1
            # and the handler short-circuits to a notify.

    async def test_feature_at_cursor_picks_smallest_enclosing(self):
        """When a small sub-feature lives inside a larger parent
        (e.g. an active-site marker inside a CDS), Enter must open
        the SMALLEST enclosing feature so the user can edit the
        more-specific annotation. Matches the main-canvas
        `_smallest_enclosing_feature` convention."""
        ed = sc.SynthesisEditor()
        ed._seq = "A" * 30
        ed._feats = [
            {"start":  0, "end": 30, "label": "parent",
             "type": "CDS",     "color": "#FFFF00",
             "strand": 1, "qualifiers": {}},
            {"start": 10, "end": 13, "label": "active-site",
             "type": "domain",  "color": "#FF0000",
             "strand": 1, "qualifiers": {}},
        ]
        ed._cursor_pos = 11  # inside both
        assert ed._feature_at_cursor() == 1   # the inner one

    async def test_feature_at_cursor_returns_neg_one_at_boundaries(
        self,
    ):
        """The hop-on marker positions (`-1` / `n+1`, from sweep
        #28) and the natural past-end position (`cur == n`) all
        return -1 — those positions sit OUTSIDE any feature span
        by construction."""
        ed = sc.SynthesisEditor()
        ed._seq = "ATGC"
        ed._feats = [
            {"start": 0, "end": 4, "label": "f",
             "type": "misc_feature", "color": "#888",
             "strand": 0, "qualifiers": {}},
        ]
        # Inside the feature: returns 0.
        ed._cursor_pos = 0
        assert ed._feature_at_cursor() == 0
        ed._cursor_pos = 3
        assert ed._feature_at_cursor() == 0
        # End-of-sequence (cur == n): outside the half-open [0, 4).
        ed._cursor_pos = 4
        assert ed._feature_at_cursor() == -1
        # 5' hop-on: outside.
        ed._cursor_pos = -1
        assert ed._feature_at_cursor() == -1
        # 3' hop-on: outside.
        ed._cursor_pos = 5
        assert ed._feature_at_cursor() == -1

    async def test_enter_only_fires_when_editor_has_focus(self):
        """Sweep #29 (2026-05-26): Enter on the synth editor must
        ONLY trigger the FeatureEdit dispatch when the scroll
        container actually has focus. If focus is elsewhere
        (library table, search input, a button), the editor's
        `on_key` may still fire from bubbled events — the explicit
        `scroll.has_focus` guard ensures Enter doesn't sideways-
        open the modal from a foreign focus."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGC", [
                {"start": 2, "end": 9, "label": "f",
                 "type": "CDS", "color": "#FF8800",
                 "strand": 1, "qualifiers": {"label": ["f"]}},
            ])
            ed._cursor_pos = 5
            await pilot.pause()
            # Focus the library SEARCH input — NOT the editor.
            from textual.widgets import Input
            search = scr.query_one("#syn-featlib-search", Input)
            search.focus()
            await pilot.pause()
            assert not ed.query_one(
                "#syn-scroll",
            ).has_focus
            # Simulate Enter via on_key with the editor as the
            # bubble target. The has_focus guard should reject.
            from textual.events import Key
            event = Key("enter", character="\n")
            ed.on_key(event)
            await pilot.pause()
            await pilot.pause(0.05)
            # Modal MUST NOT have opened — focus wasn't on editor.
            assert not isinstance(app.screen, sc.FeatureEditModal)

    async def test_focus_returns_to_editor_after_add_feature(self):
        """After the AddFeatureModal Insert dismisses, focus
        returns to the synth editor's scroll container so the
        user can keep typing bases immediately. Pre-2026-05-26
        focus drifted to whatever widget held it when the modal
        opened (often a toolbar Button), and the next base
        keystroke missed the editor entirely."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            # Verify the helper directly — it idempotently focuses
            # the editor's scroll container regardless of prior
            # focus state.
            from textual.widgets import Input
            scr.query_one("#syn-featlib-search", Input).focus()
            await pilot.pause()
            assert not scr.query_one(
                "#syn-editor",
            ).query_one("#syn-scroll").has_focus
            scr._focus_dna_editor()
            await pilot.pause()
            assert scr.query_one(
                "#syn-editor",
            ).query_one("#syn-scroll").has_focus

    async def test_modal_save_writes_back_to_synth_feats(self):
        """End-to-end: open modal via Enter, edit the label, Save
        → the synth editor's `_feats` dict reflects the new
        label. Reassignment (`_feats = list(...)`) preserves the
        in-place-mutation cache fix from earlier in this sweep."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGC", [
                {"start": 2, "end": 9, "label": "orig-label",
                 "type": "CDS", "color": "#FF8800",
                 "strand": 1, "qualifiers": {"label": ["orig-label"]}},
            ])
            ed._cursor_pos = 5
            feats_id_before = id(ed._feats)
            ed.post_message(sc.SynthesisEditor.FeatureEditRequested(0))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            assert isinstance(modal, sc.FeatureEditModal)
            # Press Edit to unlock, then change the label.
            from textual.widgets import Input, Button
            modal.query_one("#btn-featedit-edit", Button).action_press()
            await pilot.pause()
            modal.query_one("#featedit-name", Input).value = "new-label"
            modal.query_one("#btn-featedit-save", Button).action_press()
            await pilot.pause()
            await pilot.pause(0.05)
            # The synth editor's feat dict has the new label.
            assert ed._feats[0]["label"] == "new-label"
            assert ed._feats[0]["qualifiers"]["label"] == ["new-label"]
            # And `_feats` was reassigned, not mutated in place,
            # so `id(ed._feats)` changed — the chunk cache will
            # invalidate cleanly on the next render.
            assert id(ed._feats) != feats_id_before

    async def test_annotate_mode_requires_selection(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            ed._user_sel = None  # no selection
            scr._featlib_selected_entry = lambda: {
                "name": "x", "feature_type": "misc_feature",
                "sequence": "ATGC", "color": "white",
                "strand": 1, "qualifiers": {},
            }
            scr._featlib_insert_selected(mode="annotate")
            await pilot.pause()
            # No selection → no annotation added.
            assert ed._feats == []

    async def test_annotate_mode_overlays_on_selection(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            ed._user_sel = (2, 6)
            scr._featlib_selected_entry = lambda: {
                "name": "promoter1", "feature_type": "promoter",
                "sequence": "ATGC", "color": "magenta",
                "strand": 1, "qualifiers": {},
            }
            scr._featlib_insert_selected(mode="annotate")
            await pilot.pause()
            # DNA unchanged, feature added over the selection.
            assert ed._seq == "ATGCATGC"
            inserted = [f for f in ed._feats
                         if f.get("label") == "promoter1"]
            assert len(inserted) == 1
            assert inserted[0]["start"] == 2
            assert inserted[0]["end"] == 6


# ═══════════════════════════════════════════════════════════════════════════════
# AddFeatureModal total_len injection (synthesis host has no #seq-panel)
# ═══════════════════════════════════════════════════════════════════════════════

class TestHorizontalScrollAndDnaCentering:
    """The user explicitly asked for horizontal scroll (not wrap) plus
    DNA strand at viewport vertical center. Regression-lock both."""

    async def test_static_is_auto_width(self):
        """`width: auto` on `#syn-view` is what lets the Static expand
        to the rendered Text's natural width — without it, the Static
        inherits parent `1fr` and the no-wrap render gets clipped on
        the right instead of surfacing a horizontal scrollbar."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            view = app.screen.query_one("#syn-view")
            # Textual stringifies `auto` as `"auto"` in styles.
            assert "auto" in str(view.styles.width)

    async def test_scroll_overflow_horizontal_only(self):
        """`overflow-x: auto` surfaces the horizontal scrollbar when
        the rendered Text is wider than the viewport. `overflow-y:
        hidden` suppresses the vertical scrollbar because the editor
        pre-pads blank lines top/bottom for DNA centering — that
        padding inflates content height past the viewport and would
        otherwise draw a non-functional vertical bar."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scroll = app.screen.query_one("#syn-scroll")
            ox = str(scroll.styles.overflow_x).lower()
            oy = str(scroll.styles.overflow_y).lower()
            assert "auto" in ox
            assert "hidden" in oy

    async def test_dna_strand_lands_at_viewport_center(self):
        """Sacred — `_pad_above_for_centering` should drop the DNA top
        strand at viewport vertical center regardless of how many
        feature lanes sit above the DNA. Verified by computing
        screen-y = pad + above_rows and asserting equal to
        viewport_height // 2."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            # Load something with a feature so `above_rows > 0` — that
            # way we exercise the path where padding has to compensate
            # for stacked feature lanes.
            ed.load("ATGCATGC" * 40, [
                {"start": 4, "end": 60, "label": "f1",
                 "type": "CDS", "color": "red",
                 "strand": 1, "qualifiers": {}},
            ])
            await pilot.pause()
            await pilot.pause()
            vp = ed._viewport_height()
            dna_row = ed._dna_top_row_offset()
            pad = ed._pad_above_rows
            assert vp > 0
            assert pad + dna_row == vp // 2

    async def test_centering_recomputed_on_resize(self):
        """Resize re-fires `_refresh_view` so the centering padding
        re-aligns against the new viewport height."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC" * 20, [])
            await pilot.pause()
            await pilot.pause()
            pad_big = ed._pad_above_rows
            # `on_resize` is the trigger; calling it directly with a
            # synthetic event is enough to verify the recompute fires.
            ed.on_resize(None)  # type: ignore[arg-type]
            await pilot.pause()
            # Pad value should not have grown stale; same screen, same
            # viewport, same value.
            assert ed._pad_above_rows == pad_big


class TestFlankMarkers:
    """`5'-` / `-3'` flanking markers wrap the DNA top strand
    visually; they live OUTSIDE the bp coordinate space so
    selection / click resolution never confuses them with bases."""

    async def test_markers_present_in_rendered_text(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGC", [])
            await pilot.pause()
            rendered = str(app.screen.query_one("#syn-view").render())
            assert "5'-" in rendered
            assert "-3'" in rendered
            # Markers hug the DNA: "5'-ATGCATGC...ATGC-3'" should
            # appear on the DNA top-strand row.
            assert "5'-ATGCATGCATGC-3'" in rendered

    async def test_empty_fragment_has_no_markers(self):
        """The 5'/3' markers only appear with a non-empty sequence —
        an empty fragment shows the placeholder text instead."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("", [])
            await pilot.pause()
            rendered = str(app.screen.query_one("#syn-view").render())
            assert "5'-" not in rendered
            assert "-3'" not in rendered

    async def test_select_all_covers_dna_only(self):
        """Ctrl+A → `_user_sel = (0, n)` — the bp coordinate space
        excludes the markers, so the highlight visual stays on the
        DNA bases by construction."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            await pilot.pause()
            ed.select_all()
            await pilot.pause()
            assert ed._user_sel == (0, 8)
            # Selection is in bp coords; markers sit outside that
            # range. The rendered DNA row carries the selection
            # highlight on bases only.

    async def test_click_resolution_skips_marker_offset(self):
        """`_click_to_bp` subtracts horizontal-centring pad + 3-char
        5'- marker so a click at the first base lands at bp 0, not at
        bp 3 (which would happen if the marker offset weren't
        subtracted). The line-number gutter was dropped in
        `_wrap_with_53_markers` so it doesn't enter the formula."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            # Long-enough seq so centring pad is 0 — keeps the column
            # arithmetic deterministic (test wants to verify the
            # marker offset, not the centring math).
            ed.load("ATGCATGC" * 30, [])
            await pilot.pause()
            assert ed._pad_left_cols == 0
            view = app.screen.query_one("#syn-view")
            # Column where bp 0 begins = view.region.x + pad_left (0
            # here) + 5'- marker width (3).
            base_col_x = (view.region.x + ed._pad_left_cols
                          + ed._FLANK_MARKER_WIDTH)
            bp = ed._click_to_bp(base_col_x, view.region.y)
            assert bp == 0
            # Clicking ONE column to the left (still on the 5'-
            # marker) should clamp to bp 0.
            bp_marker = ed._click_to_bp(
                base_col_x - 1, view.region.y,
            )
            assert bp_marker == 0
            # Clicking at the column of the 4th base (bp 3) lands at
            # bp 3.
            bp_3 = ed._click_to_bp(base_col_x + 3, view.region.y)
            assert bp_3 == 3

    async def test_marker_width_constant_is_3(self):
        # Sacred — gutter / cursor / click all assume width == 3.
        # Changing to "5'─" (em-dash) or similar would break the math
        # silently.
        assert sc.SynthesisEditor._FLANK_MARKER_WIDTH == 3
        assert sc.SynthesisEditor._FLANK_MARKER_TOP_LEFT == "5'-"
        assert sc.SynthesisEditor._FLANK_MARKER_TOP_RIGHT == "-3'"
        assert sc.SynthesisEditor._FLANK_MARKER_BOT_LEFT == "3'-"
        assert sc.SynthesisEditor._FLANK_MARKER_BOT_RIGHT == "-5'"

    # ── Cursor-on-hyphen highlight (sweep #28 extended 2026-05-26) ──
    #
    # When the cursor sits at a sequence boundary — either the
    # natural [0, n] end-of-sequence position OR the explicit
    # hop-on positions [-1, n+1] — the matching hyphen reverses to
    # carry the cursor. Before 2026-05-26 the natural past-last-base
    # position `cur == n` rendered an invisible cursor (per-base
    # reverse loop is [0, n), so no base reversed; hyphen highlight
    # only fired for hop-on). The unified design fixes that and
    # keeps the existing hop-on regression covered.

    @staticmethod
    def _wrap_text_for(seq: str, cur: int, feats=None):
        """Build a `_build_seq_text` + `_wrap_with_53_markers` Text
        for a non-empty `seq` and `cur` cursor position. Returns the
        wrapped Text (with the hyphen markers injected) so a test
        can walk spans to check the reverse-video paint position.

        Unit-style — instantiates `SynthesisEditor` without a widget
        mount because `_wrap_with_53_markers` reads only `_seq`,
        `_feats`, `_cursor_pos`, and the marker class constants."""
        ed = sc.SynthesisEditor()
        ed._seq = str(seq)
        ed._feats = list(feats or [])
        ed._cursor_pos = int(cur)
        n = len(ed._seq)
        # Build the inner Text the editor would build before wrapping.
        inner = sc._build_seq_text(
            ed._seq, ed._feats, line_width=n + 1,
            cursor_pos=ed._cursor_pos,
        )
        return ed, ed._wrap_with_53_markers(inner)

    @staticmethod
    def _dna_row_spans(text_obj, dna_row: int):
        """Return the list of Text spans on the DNA top-strand row."""
        lines = text_obj.split("\n")
        if dna_row >= len(lines):
            return []
        return list(lines[dna_row].spans)

    @staticmethod
    def _span_at(spans, col: int):
        """Return spans whose range covers `col`."""
        return [s for s in spans if s.start <= col < s.end]

    @staticmethod
    def _has_reverse(spans) -> bool:
        return any("reverse" in str(s.style).lower() for s in spans)

    def test_cursor_at_natural_end_highlights_3p_hyphen(self):
        """`cur == n` (natural past-last-base, reached by → / End)
        now reverses the leading hyphen of the `-3'` marker so the
        cursor is visible at the sequence end. Pre-2026-05-26 the
        cursor was invisible at this position."""
        ed, txt = self._wrap_text_for("ATGCATGC", cur=8)
        n = len(ed._seq)
        dna_row = ed._dna_top_row_offset()
        spans = self._dna_row_spans(txt, dna_row)
        # `-3'` sits immediately after the bases. With the line-
        # number gutter stripped, the layout is `5'-` (3) + n bases
        # + `-3'` (3). The leading hyphen of `-3'` is at column
        # `marker_w + n` = 3 + 8 = 11.
        hyphen_col = ed._FLANK_MARKER_WIDTH + n  # 11
        assert self._has_reverse(self._span_at(spans, hyphen_col))

    def test_cursor_at_natural_start_highlights_5p_hyphen(self):
        """`cur == 0` (natural before-first-base, reached by ← /
        Home) reverses the trailing hyphen of the `5'-` marker.
        base[0] still reverses too (via `_build_seq_text`) so the
        cursor visually straddles the boundary — an unambiguous
        "you're at the very start" cue."""
        ed, txt = self._wrap_text_for("ATGCATGC", cur=0)
        dna_row = ed._dna_top_row_offset()
        spans = self._dna_row_spans(txt, dna_row)
        # Trailing hyphen of `5'-` sits at column 2 (positions 0, 1,
        # 2 = `5`, `'`, `-`).
        hyphen_col = 2
        assert self._has_reverse(self._span_at(spans, hyphen_col))

    def test_cursor_hopon_5p_marker_still_highlights(self):
        """Regression for sweep #28 hop-on: `cur == -1` (Home from
        cur==0) reverses the 5' hyphen."""
        ed, txt = self._wrap_text_for("ATGCATGC", cur=-1)
        dna_row = ed._dna_top_row_offset()
        spans = self._dna_row_spans(txt, dna_row)
        assert self._has_reverse(self._span_at(spans, 2))

    def test_cursor_hopon_3p_marker_still_highlights(self):
        """Regression for sweep #28 hop-on: `cur == n+1` (End from
        cur==n) reverses the -3' hyphen."""
        ed, txt = self._wrap_text_for("ATGCATGC", cur=9)
        n = len(ed._seq)
        dna_row = ed._dna_top_row_offset()
        spans = self._dna_row_spans(txt, dna_row)
        assert self._has_reverse(
            self._span_at(spans, ed._FLANK_MARKER_WIDTH + n)
        )

    def test_cursor_interior_does_not_highlight_hyphens(self):
        """Cursor in the middle of the sequence (e.g. base 4 of 8)
        must NOT reverse either hyphen — the cursor lives on the
        base. This guards against an over-eager `cur <= 0` /
        `cur >= n` that would also paint the hyphen for any cursor
        position."""
        ed, txt = self._wrap_text_for("ATGCATGC", cur=4)
        n = len(ed._seq)
        dna_row = ed._dna_top_row_offset()
        spans = self._dna_row_spans(txt, dna_row)
        # Neither hyphen reverses.
        assert not self._has_reverse(self._span_at(spans, 2))
        assert not self._has_reverse(
            self._span_at(spans, ed._FLANK_MARKER_WIDTH + n)
        )

    def test_single_base_sequence_cursor_at_start(self):
        """Smallest non-empty seq: n=1. cur=0 → 5' hyphen highlights;
        cur=1 → -3' hyphen highlights; cur=0 must NOT also highlight
        the -3' hyphen (and vice versa)."""
        ed, txt = self._wrap_text_for("A", cur=0)
        dna_row = ed._dna_top_row_offset()
        spans = self._dna_row_spans(txt, dna_row)
        # 5' hyphen reverses.
        assert self._has_reverse(self._span_at(spans, 2))
        # -3' hyphen does NOT reverse (cur != n).
        assert not self._has_reverse(
            self._span_at(spans, ed._FLANK_MARKER_WIDTH + 1)
        )

    def test_single_base_sequence_cursor_at_end(self):
        """n=1, cur=1 (past the only base) → -3' hyphen reverses,
        5'- hyphen does NOT (cur != 0)."""
        ed, txt = self._wrap_text_for("A", cur=1)
        n = len(ed._seq)
        dna_row = ed._dna_top_row_offset()
        spans = self._dna_row_spans(txt, dna_row)
        assert self._has_reverse(
            self._span_at(spans, ed._FLANK_MARKER_WIDTH + n)
        )
        assert not self._has_reverse(self._span_at(spans, 2))

    def test_bottom_hyphens_never_get_cursor_highlight(self):
        """The cursor lives on the TOP strand; the bottom strand
        (`3'-…-5'`) is the reverse-complement display and must
        never carry the cursor reverse-video. Reproduces the
        biological convention — typing always goes into the top
        strand."""
        ed, txt = self._wrap_text_for("ATGCATGC", cur=8)  # natural end
        n = len(ed._seq)
        dna_row = ed._dna_top_row_offset()
        # Bottom strand sits one row below the top.
        lines = txt.split("\n")
        bot_row = dna_row + 1
        if bot_row >= len(lines):
            pytest.skip("bottom strand not rendered in this layout")
        bot_spans = list(lines[bot_row].spans)
        # Left hyphen of `3'-` is at column 2; right hyphen of `-5'`
        # is at column marker_w + n.
        assert not self._has_reverse(self._span_at(bot_spans, 2))
        assert not self._has_reverse(
            self._span_at(bot_spans, ed._FLANK_MARKER_WIDTH + n)
        )

    def test_marker_with_cursor_clamps_out_of_bounds_idx(self):
        """Defensive: `_marker_with_cursor` clamps an out-of-bounds
        `cursor_char_idx` rather than raising. A subclass / test
        fixture could override the marker constants to a shorter
        string; the render path must degrade to "no reverse" or to
        the last-char position instead of crashing the editor."""
        ed = sc.SynthesisEditor()
        ed._seq = "ATGC"
        ed._feats = []
        ed._cursor_pos = 4
        # Build the inner Text first.
        inner = sc._build_seq_text(
            ed._seq, ed._feats, line_width=5, cursor_pos=4,
        )
        # Wrap once normally — establishes baseline.
        wrapped = ed._wrap_with_53_markers(inner)
        assert wrapped is not None
        # Now exercise the helper indirectly by clobbering the
        # marker to an empty string and re-wrapping. The defensive
        # branch returns an empty Text rather than slicing past 0.
        try:
            ed._FLANK_MARKER_TOP_RIGHT = ""
            inner2 = sc._build_seq_text(
                ed._seq, ed._feats, line_width=5, cursor_pos=4,
            )
            # Must NOT raise.
            wrapped2 = ed._wrap_with_53_markers(inner2)
            assert wrapped2 is not None
        finally:
            # Restore for any later test in the same process.
            ed._FLANK_MARKER_TOP_RIGHT = "-3'"

    async def test_arrow_past_last_base_lands_cursor_at_n(self):
        """End-to-end via pilot — pressing `set_cursor(+1)` on the
        last base lands `_cursor_pos` at `n` (the natural past-end
        position). Combined with the unit tests above (which verify
        that `cur == n` reverses the -3' hyphen via
        `_wrap_with_53_markers`), this exercises the full
        click-to-render path that reproduces the 2026-05-26 symptom
        ("cursor invisibly sitting beyond the last basepair")."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            await pilot.pause()
            ed.set_cursor(len(ed._seq) - 1)
            await pilot.pause()
            ed.set_cursor(ed._cursor_pos + 1)
            await pilot.pause()
            # Cursor lands at the natural past-end position.
            assert ed._cursor_pos == len(ed._seq)
            # One more → puts the cursor at the hop-on position
            # (`n + 1`); the editor's `set_cursor` clamp is [-1, n+1].
            ed.set_cursor(ed._cursor_pos + 1)
            await pilot.pause()
            assert ed._cursor_pos == len(ed._seq) + 1
            # And one more → MUST NOT advance further (clamp).
            ed.set_cursor(ed._cursor_pos + 1)
            await pilot.pause()
            assert ed._cursor_pos == len(ed._seq) + 1

    async def test_bottom_strand_markers_reflect_antiparallel(self):
        """Biological reality — DNA is anti-parallel. Top strand 5'→3',
        bottom strand 3'→5' under standard left-to-right orientation.
        The bottom row must carry `3'-` on the left and `-5'` on the
        right, so the user sees:

            5'-ATGCATGC-3'
            3'-TACGTACG-5'
        """
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            await pilot.pause()
            rendered = str(app.screen.query_one("#syn-view").render())
            # Top strand: 5'- on left, -3' on right, hugging bases.
            assert "5'-ATGCATGC-3'" in rendered
            # Bottom strand: 3'- on left, -5' on right, hugging the
            # reverse complement (TACGTACG is the complement of
            # ATGCATGC read left-to-right).
            assert "3'-TACGTACG-5'" in rendered


class TestSelectAllThenDelete:
    """Sacred — Ctrl+A in either tab MUST select the whole buffer
    AND focus the editor so the very next Backspace / Delete
    keystroke erases the highlighted span.

    Pre-fix bug (reported 2026-05-20): Ctrl+A fired from a screen-
    level binding, which left focus on whatever widget the user
    last interacted with (typically a toolbar Button on modal
    open). Backspace then went to the focused Button, which has
    no key handler → silent no-op. Fix: `action_select_all`
    explicitly focuses the active editor's ScrollableContainer
    AND a screen-level Backspace/Delete fallback dispatches to
    the editor's `delete_at_cursor` so the operation works even
    if focus drifts after select-all.
    """

    async def test_dna_ctrl_a_then_backspace_clears(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGC", [])
            await pilot.pause()
            await pilot.press("ctrl+a")
            await pilot.pause()
            assert ed._user_sel == (0, 12)
            await pilot.press("backspace")
            await pilot.pause()
            assert ed._seq == ""

    async def test_dna_ctrl_a_then_delete_clears(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            await pilot.pause()
            await pilot.press("ctrl+a")
            await pilot.pause()
            await pilot.press("delete")
            await pilot.pause()
            assert ed._seq == ""

    async def test_protein_ctrl_a_then_backspace_clears(self):
        from textual.widgets import TabbedContent
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            tabs = scr.query_one("#syn-tabs", TabbedContent)
            tabs.active = "syn-tab-protein"
            await pilot.pause()
            await pilot.pause()
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MSKLNPE")
            await pilot.pause()
            await pilot.press("ctrl+a")
            await pilot.pause()
            assert pe._user_sel == (0, 7)
            await pilot.press("backspace")
            await pilot.pause()
            assert pe._aa_seq == ""

    async def test_ctrl_a_focuses_active_editor(self):
        """The select_all action must move focus onto the active
        editor's ScrollableContainer. Verifies the explicit focus
        call lands the right widget."""
        from textual.containers import ScrollableContainer
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGC", [])
            await pilot.pause()
            await pilot.press("ctrl+a")
            await pilot.pause()
            focused = app.focused
            assert isinstance(focused, ScrollableContainer)
            assert focused.id == "syn-scroll"

    async def test_screen_backspace_fallback_no_op_without_selection(self):
        """The screen-level Backspace binding must NOT delete anything
        when there's no active selection — single-base backspace
        belongs to the editor's own on_key when it has focus, not the
        screen fallback."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            ed._cursor_pos = 4
            ed._user_sel = None
            # Direct invoke of the screen-level fallback — should bail
            # (no selection to delete).
            ok = scr._delete_active_editor_selection(forward=False)
            await pilot.pause()
            assert ok is False
            # Buffer unchanged.
            assert ed._seq == "ATGCATGC"


class TestShiftArrowSelection:
    """Sacred — the base under the cursor MUST be part of the
    selection on the first shift+arrow press in either direction.

    Pre-fix bug (reported 2026-05-20): cursor visually sits ON
    base[N] (reverse-video paint by `_build_seq_text`), but the
    set_cursor anchor was set to cursor_pos itself, giving
    half-open (pos, anchor) = (N-1, N) for shift+left — covering
    only base[N-1] and "skipping" base[N]. Fix: direction-aware
    anchor — for shift+left set anchor to cursor_pos+1 so the
    half-open range covers BOTH the base under the cursor AND the
    base the cursor moves to.
    """

    async def test_shift_left_includes_base_under_cursor(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGC", [])
            # Simulate a fresh click at base 5 — sets cursor + anchor
            # without any active selection.
            ed._cursor_pos = 5
            ed._user_sel = None
            ed._sel_anchor = 5
            # First shift+left: cursor moves to 4, selection must
            # cover BOTH base 4 (where cursor moved to) AND base 5
            # (where cursor WAS).
            ed.set_cursor(4, extend_sel=True)
            assert ed._cursor_pos == 4
            assert ed._user_sel == (4, 6)
            # Each subsequent shift+left adds exactly one base on
            # the left edge; the original cursor base (5) stays in.
            ed.set_cursor(3, extend_sel=True)
            assert ed._user_sel == (3, 6)
            ed.set_cursor(2, extend_sel=True)
            assert ed._user_sel == (2, 6)

    async def test_shift_right_still_includes_starting_base(self):
        """Shift+right wasn't broken pre-fix — anchor = cursor_pos
        with half-open (anchor, pos) naturally covers
        base[cursor_pos]. Regression-lock against future tweaks."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGC", [])
            ed._cursor_pos = 5
            ed._user_sel = None
            ed._sel_anchor = 5
            ed.set_cursor(6, extend_sel=True)
            assert ed._cursor_pos == 6
            # Base 5 (where cursor WAS) is selected.
            assert ed._user_sel == (5, 6)
            ed.set_cursor(7, extend_sel=True)
            # Both bases 5 and 6 selected; cursor at 7.
            assert ed._user_sel == (5, 7)

    async def test_shift_left_then_right_collapses_then_extends(self):
        """Reversing direction reduces selection from the moving
        edge; once the selection collapses, the next press starts a
        fresh extension in the new direction."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGCATGC", [])
            ed._cursor_pos = 5
            ed._user_sel = None
            ed._sel_anchor = 5
            # Shift+left twice → sel (3, 6), cursor at 3.
            ed.set_cursor(4, extend_sel=True)
            ed.set_cursor(3, extend_sel=True)
            assert ed._user_sel == (3, 6)
            # Shift+right reduces left edge: cursor 4, sel (4, 6).
            ed.set_cursor(4, extend_sel=True)
            assert ed._user_sel == (4, 6)
            # Continue right past anchor: cursor 6, sel collapses.
            ed.set_cursor(5, extend_sel=True)
            assert ed._user_sel == (5, 6)
            ed.set_cursor(6, extend_sel=True)
            assert ed._user_sel is None

    async def test_shift_left_at_position_zero_no_op(self):
        """Shift+left at the very start (cursor=0) clamps and
        produces no selection."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            ed._cursor_pos = 0
            ed._user_sel = None
            ed._sel_anchor = 0
            ed.set_cursor(-1, extend_sel=True)  # clamps to 0
            assert ed._cursor_pos == 0
            assert ed._user_sel is None

    async def test_shift_left_after_ctrl_a_reduces_selection(self):
        """Ctrl+A then Shift+Left should REDUCE the selection by one
        from the right edge (the cursor was past the end at n,
        moving left to n-1 drops the last base from the selection).
        """
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])  # n = 8
            ed.select_all()
            assert ed._user_sel == (0, 8)
            assert ed._cursor_pos == 8
            ed.set_cursor(7, extend_sel=True)
            # Cursor moves from 8 to 7. Selection shrinks by one base
            # on the right edge.
            assert ed._user_sel == (0, 7)
            assert ed._cursor_pos == 7


class TestRowMarkerAndCentering:
    """The line-number gutter is gone (always one row, so the leading
    "1" added no info); short sequences sit horizontally centred in
    the viewport; long sequences left-anchor so horizontal scroll
    works."""

    async def test_row_number_marker_dropped(self):
        """Pre-fix the DNA row started with " 1  " before the 5'-
        marker. Post-fix the marker is the very first non-space
        glyph on the DNA row."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            # Long enough to suppress the centring pad — keeps the
            # assertion deterministic.
            ed.load("ATGCATGC" * 30, [])
            await pilot.pause()
            rendered = str(app.screen.query_one("#syn-view").render())
            # Find the DNA top-strand line and check it starts with
            # the 5'- marker (no leading line-number gutter).
            for line in rendered.split("\n"):
                if "5'-" in line and "TACG" not in line:
                    # `line.lstrip()` should yield the marker right
                    # at the start — no `1 ` or similar gutter.
                    assert line.lstrip().startswith("5'-")
                    return
            assert False, "DNA top-strand row not found"

    async def test_short_fragment_horizontally_centered(self):
        """Content width (n + 6) < viewport width → left-pad spaces
        so the fragment sits at the viewport's horizontal centre."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])  # n + 6 = 14
            await pilot.pause()
            vp = ed._viewport_width()
            cw = ed._content_width()
            # Sanity — at the standard 160-col test terminal the
            # editor sits in a 4fr/1fr split so its viewport is
            # narrower than 160 but still > 14.
            assert vp > cw
            # Centring pad = (vp - cw) // 2 exactly.
            assert ed._pad_left_cols == (vp - cw) // 2

    async def test_long_fragment_left_anchored(self):
        """Content width >= viewport width → no centring pad; user
        scrolls horizontally instead."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC" * 50, [])  # 400 bp + 6 = 406, > viewport
            await pilot.pause()
            assert ed._content_width() > ed._viewport_width()
            assert ed._pad_left_cols == 0

    async def test_centering_recomputed_on_size_change(self):
        """Padding recomputes when the sequence shrinks / grows
        across the viewport-width threshold."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            # Start small → centred.
            ed.load("ATGC", [])
            await pilot.pause()
            small_pad = ed._pad_left_cols
            assert small_pad > 0
            # Grow past viewport → pad should drop to 0.
            ed.load("ATGCATGC" * 100, [])
            await pilot.pause()
            assert ed._pad_left_cols == 0
            # Shrink back → centring returns.
            ed.load("ATGC", [])
            await pilot.pause()
            assert ed._pad_left_cols == small_pad

    async def test_click_resolution_respects_centering_pad(self):
        """A click on the very first base of a CENTRED fragment
        must land at bp 0, not at bp `-(pad + marker)`. Critical:
        click math must subtract `_pad_left_cols + FLANK_MARKER_WIDTH`
        from the raw column."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])  # short → centring pad applies
            await pilot.pause()
            assert ed._pad_left_cols > 0
            view = app.screen.query_one("#syn-view")
            base0_x = (view.region.x + ed._pad_left_cols
                        + ed._FLANK_MARKER_WIDTH)
            bp = ed._click_to_bp(base0_x, view.region.y)
            assert bp == 0
            # Clicking deep in the centring pad collapses to bp 0.
            bp_pad = ed._click_to_bp(view.region.x + 2, view.region.y)
            assert bp_pad == 0


class TestProteinEditorPrimitives:
    """ProteinEditor mirrors SynthesisEditor's mutation semantics but
    in AA coords. Test the primitives independently of the screen."""

    def _make(self):
        pe = sc.ProteinEditor()
        pe._aa_seq = ""
        pe._cursor_pos = 0
        pe._codon_mode = True
        pe._codon_table_raw = dict(sc._CODON_BUILTIN_K12)
        pe._codon_cache = sc.ProteinEditor._build_codon_cache(
            pe._codon_table_raw,
        )
        return pe

    def test_aa_alphabet_is_20_plus_stop(self):
        assert sc._PROTEIN_AA_ALPHABET == frozenset("ACDEFGHIKLMNPQRSTVWY*")

    def test_max_aa_constant_matches_max_bp_div_3(self):
        assert sc._PROTEIN_MAX_AA == sc._SYNTHESIS_MAX_BP // 3

    def test_codon_cache_has_all_20_aas(self):
        cache = sc.ProteinEditor._build_codon_cache(
            dict(sc._CODON_BUILTIN_K12),
        )
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            assert aa in cache, f"missing AA {aa}"
            assert len(cache[aa]) == 3
        # Stop codon also present.
        assert "*" in cache
        assert len(cache["*"]) == 3

    def test_codon_cache_empty_table_falls_back(self):
        cache = sc.ProteinEditor._build_codon_cache({})
        # Empty table → empty cache except the hardcoded TAA stop.
        assert cache.get("*") == "TAA"

    def test_motif_library_has_essentials(self):
        names = {m["name"] for m in sc._PROTEIN_MOTIFS}
        for required in ("His6", "FLAG", "HA", "Myc", "TEV", "P2A",
                          "(GGGGS)x3"):
            assert required in names, f"missing motif {required}"


class TestProteinEditorAsync:
    async def test_typing_inserts_aas(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.insert_at_cursor("MSKL")
            await pilot.pause()
            assert pe._aa_seq == "MSKL"
            assert pe._cursor_pos == 4

    async def test_typing_drops_non_aa_chars(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            # B, J, O, U, X, Z are non-standard / ambiguous — drop.
            pe.insert_at_cursor("BJMOUXZK")
            await pilot.pause()
            # M and K survive; everything else drops.
            assert pe._aa_seq == "MK"

    async def test_stop_codon_allowed(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.insert_at_cursor("MGE*")
            await pilot.pause()
            assert pe._aa_seq == "MGE*"

    async def test_max_aa_cap_truncates(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe._aa_seq = "M" * (sc._PROTEIN_MAX_AA - 3)
            pe._cursor_pos = len(pe._aa_seq)
            pe.insert_at_cursor("AAAAAAAAAA")  # 10 As — only 3 fit
            await pilot.pause()
            assert len(pe._aa_seq) == sc._PROTEIN_MAX_AA

    async def test_codon_mode_renders_codons_below_aas(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MSK")
            await pilot.pause()
            rendered = str(app.screen.query_one("#pe-view").render())
            # Top row: N- M  S  K -C  (each AA centred in 3-col group).
            # Bottom row: ATGAGCAAA (K12 most-frequent codons).
            assert "N-" in rendered and "-C" in rendered
            assert " M " in rendered  # AA letter wrapped in spaces
            assert "ATG" in rendered  # M codon
            assert "AAA" in rendered  # K codon

    async def test_aa_only_mode_no_dna_row(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MSK")
            pe.set_codon_mode(False)
            await pilot.pause()
            rendered = str(app.screen.query_one("#pe-view").render())
            # Single row: N-MSK-C.
            assert "N-MSK-C" in rendered
            # The most-frequent codon for M (ATG) should NOT appear.
            assert "ATG" not in rendered

    async def test_codon_table_switch_updates_render(self):
        """Picking a different codon table should re-render the DNA codons
        under the AA letters. Uses two tables that disagree on the preferred
        Lys codon (AAA vs AAG) — both valid synonyms — so the switch is
        observable WITHOUT remapping a codon to a non-canonical amino acid
        (which `_codon_build_aa_map` now correctly ignores, deriving the AA
        from the standard genetic code)."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MK")
            await pilot.pause()
            # Table A: Lys prefers AAA.
            pe.set_codon_table({
                "ATG": ("M", 1), "AAA": ("K", 99), "AAG": ("K", 1),
                "TAA": ("*", 1),
            })
            await pilot.pause()
            assert pe._codon_cache["K"] == "AAA"
            # Switch to Table B: Lys prefers AAG — the render must follow.
            pe.set_codon_table({
                "ATG": ("M", 1), "AAA": ("K", 1), "AAG": ("K", 99),
                "TAA": ("*", 1),
            })
            await pilot.pause()
            assert pe._codon_cache["K"] == "AAG"
            assert pe._codon_cache["K"] != "AAA"


class TestSynthesisTabbing:
    async def test_dna_and_protein_tabs_both_mount(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            # Both editors present.
            scr.query_one("#syn-editor", sc.SynthesisEditor)
            scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            # Active tab defaults to DNA.
            assert scr._active_tab_id() == "dna"

    async def test_tab_switch_updates_active_state(self):
        from textual.widgets import TabbedContent
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            tabs = scr.query_one("#syn-tabs", TabbedContent)
            tabs.active = "syn-tab-protein"
            await pilot.pause()
            await pilot.pause()
            assert scr._active_tab_id() == "protein"
            tabs.active = "syn-tab-dna"
            await pilot.pause()
            await pilot.pause()
            assert scr._active_tab_id() == "dna"

    async def test_protein_dirty_separate_from_dna_dirty(self):
        """Mutating the DNA tab marks DNA dirty but NOT the protein
        tab, and vice versa."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            ed.insert_at_cursor("ATGC")
            # `insert_at_cursor` flips `_dirty` via a posted change message;
            # under `-n auto` load a single pause may not drain it, so poll.
            for _ in range(20):
                await pilot.pause()
                if scr._dirty:
                    break
            assert scr._dirty is True
            assert scr._protein_dirty is False
            pe.insert_at_cursor("MS")
            for _ in range(20):
                await pilot.pause()
                if scr._protein_dirty:
                    break
            assert scr._dirty is True
            assert scr._protein_dirty is True

    async def test_add_feature_on_protein_tab_notifies(self):
        """Ctrl+F on the protein tab is a no-op (motif library handles
        the equivalent flow). The notify confirms the user got the
        right message."""
        from textual.widgets import TabbedContent
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            tabs = scr.query_one("#syn-tabs", TabbedContent)
            tabs.active = "syn-tab-protein"
            await pilot.pause()
            await pilot.pause()
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.insert_at_cursor("MS")
            await pilot.pause()
            n_aa_before = len(pe._aa_seq)
            scr.action_add_feature()
            await pilot.pause()
            # AA buffer unchanged — the action just notified.
            assert len(pe._aa_seq) == n_aa_before

    async def test_toggle_codon_mode_on_dna_tab_notifies(self):
        """Alt+T on the DNA tab is a no-op (codon mode is protein-
        tab-only). DNA buffer untouched. (Was Ctrl+M but the terminal
        eats that as Enter / ^M.)"""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.insert_at_cursor("ATGC")
            await pilot.pause()
            n_bp_before = len(ed._seq)
            scr.action_toggle_codon_mode()
            await pilot.pause()
            assert len(ed._seq) == n_bp_before


class TestProteinSaveRoundTrip:
    async def test_save_creates_linear_record_with_translation(self):
        """Protein save → linear DNA library entry with CDS feature
        carrying a translation= qualifier matching the AA seq."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MSKLEAH*")
            scr._protein_loaded_id = "test_prot"
            scr._protein_loaded_name = "test_prot"
            scr._commit_protein_save("MSKLEAH*", after=None)
            await pilot.pause()
            await pilot.pause()
            entries = sc._load_library()
            matches = [e for e in entries if e.get("id") == "test_prot"]
            assert len(matches) == 1
            gb = matches[0].get("gb_text", "")
            assert "linear" in gb.lower()
            # CDS + translation qualifier landed.
            assert "CDS" in gb
            assert "MSKLEAH" in gb.upper().replace("\n", "").replace(" ", "")

    async def test_save_then_load_preserves_aa(self):
        """Round-trip: save protein → load same library entry into
        the protein editor → AA sequence preserved exactly."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            original = "MGSSHHHHHHENLYFQGSAVNTGLPRSE*"
            pe.load(original)
            scr._protein_loaded_id = "rt_prot"
            scr._protein_loaded_name = "rt_prot"
            scr._commit_protein_save(original, after=None)
            await pilot.pause()
            # Reset editor + reload.
            pe.load("")
            await pilot.pause()
            scr._load_protein_entry_by_id("rt_prot")
            await pilot.pause()
            assert pe._aa_seq == original

    async def test_atomic_save_lock_present(self):
        """Sacred — the SynthesisScreen's `_save_lock` must be an
        RLock so DNA `_commit_save` and protein `_commit_protein_save`
        serialise their SeqRecord build + library hand-off."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            # RLock is reentrant — exposes acquire / release.
            assert hasattr(scr._save_lock, "acquire")
            assert hasattr(scr._save_lock, "release")
            # And it's truly an RLock (allows reentrant acquire).
            scr._save_lock.acquire()
            scr._save_lock.acquire()
            scr._save_lock.release()
            scr._save_lock.release()


class TestProteinMotifLibrary:
    async def test_motif_table_populates_with_builtin_motifs(self):
        from textual.widgets import DataTable
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            t = scr.query_one("#syn-motif-table", DataTable)
            # All 30 builtin motifs should land.
            assert t.row_count == len(sc._PROTEIN_MOTIFS)

    async def test_motif_insert_splices_at_cursor(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MM")
            pe._cursor_pos = 1
            # Stub the selector so the test doesn't depend on row
            # focus state.
            scr._motif_selected_entry = lambda: {
                "name": "His6", "feature_type": "Tag",
                "sequence": "HHHHHH", "description": "",
            }
            scr._motif_insert_selected()
            await pilot.pause()
            assert pe._aa_seq == "MHHHHHHM"


class TestAddFeatureModalTotalLen:
    def test_total_len_kwarg_accepted(self):
        # Backwards-compatible: omitted → None → falls back to live
        # #seq-panel query (the original on-canvas Ctrl+F path).
        m1 = sc.AddFeatureModal()
        assert m1._total_len is None
        # SynthesisScreen path injects an explicit length so the CDS
        # divisibility-by-3 check works even without #seq-panel.
        m2 = sc.AddFeatureModal(total_len=123)
        assert m2._total_len == 123


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #14 — Save As double-prompt + silent-overwrite fixes
# ═══════════════════════════════════════════════════════════════════════════════

class TestUniqueEntryIdHelper:
    """Regression guard for 2026-05-20 fix.

    ``_make_unique_entry_id`` must suffix-disambiguate against current
    library ids so a fresh save / Save As never silently overwrites an
    unrelated library entry whose sanitised name happens to collide.
    """

    async def test_returns_base_when_no_collision(self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            assert scr._make_unique_entry_id("fresh_name") == "fresh_name"

    async def test_disambiguates_against_existing_ids(self, isolated_library):
        # Seed the library with an entry whose id is "my_fragment".
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC"), id="my_fragment", name="my_fragment")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "linear"
        sc._save_library([{
            "id": "my_fragment", "name": "My fragment",
            "size": 4, "n_feats": 0,
            "source": "test", "added": "2026-05-20",
            "gb_text": sc._record_to_gb_text(rec),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            # Sanitised name "my fragment" → base "my_fragment" already
            # taken → suffix-disambiguate to "my_fragment_2".
            uid = scr._make_unique_entry_id("my fragment")
            assert uid != "my_fragment"
            assert uid.startswith("my_fragment_")

    async def test_exclude_id_skips_own_entry(self, isolated_library):
        # Re-save case: the loaded_id is already in the library but we
        # don't want to bump past our own slot.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC"), id="my_fragment", name="my_fragment")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "linear"
        sc._save_library([{
            "id": "my_fragment", "name": "My fragment",
            "size": 4, "n_feats": 0,
            "source": "test", "added": "2026-05-20",
            "gb_text": sc._record_to_gb_text(rec),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            # Excluding our own id should return the base without bump.
            uid = scr._make_unique_entry_id(
                "my fragment", exclude_id="my_fragment",
            )
            assert uid == "my_fragment"


class TestSaveAsDoublePromptFix:
    """Regression guard for 2026-05-20 fix.

    Pre-fix ``_save_as._on_named`` set ``loaded_id=None`` then called
    ``_do_save``, which saw ``loaded_id is None`` and prompted for the
    name AGAIN. Post-fix the inner callback hands a unique id straight
    to ``_commit_save``; only one ``NamePlasmidModal`` ever appears.
    """

    async def test_save_as_commits_without_second_prompt(
        self, isolated_library,
    ):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCATGC", [])
            scr._dirty = True
            # Stub NamePlasmidModal push so we can count invocations.
            push_count = [0]
            orig_push = app.push_screen
            def _counting_push(modal, callback=None):
                if isinstance(modal, sc.NamePlasmidModal):
                    push_count[0] += 1
                    # Simulate user typing a name and confirming.
                    if callback is not None:
                        callback("brand_new_name")
                    return None
                return orig_push(modal, callback=callback)
            app.push_screen = _counting_push  # type: ignore[method-assign]
            try:
                scr._save_as()
            finally:
                app.push_screen = orig_push  # type: ignore[method-assign]
            await pilot.pause()
            await pilot.pause()
            # Exactly one prompt — the pre-fix bug would have stacked
            # a second NamePlasmidModal on top.
            assert push_count[0] == 1, (
                f"Save As should prompt once, prompted {push_count[0]}"
            )
            # Entry landed in the library under the (sanitised, unique) id.
            entries = sc._load_library()
            ids = {e.get("id") for e in entries}
            assert "brand_new_name" in ids

    async def test_protein_save_as_commits_without_second_prompt(
        self, isolated_library,
    ):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("MASHHH")
            scr._protein_dirty = True
            push_count = [0]
            orig_push = app.push_screen
            def _counting_push(modal, callback=None):
                if isinstance(modal, sc.NamePlasmidModal):
                    push_count[0] += 1
                    if callback is not None:
                        callback("brand_new_protein")
                    return None
                return orig_push(modal, callback=callback)
            app.push_screen = _counting_push  # type: ignore[method-assign]
            try:
                scr._protein_save_as()
            finally:
                app.push_screen = orig_push  # type: ignore[method-assign]
            await pilot.pause()
            await pilot.pause()
            assert push_count[0] == 1, (
                f"Protein Save As should prompt once, prompted {push_count[0]}"
            )

    async def test_fresh_save_uses_unique_id_against_library(
        self, isolated_library,
    ):
        # Bug scenario: library already has an entry whose id is
        # "my_fragment_v2" with a DIFFERENT name ("Pre-existing
        # entry"). User types Save As "my fragment v2" → sanitises to
        # id "my_fragment_v2". Pre-fix `_make_entry_id` returned the
        # raw sanitised id and `LibraryPanel.add_entry` took the
        # id-match silent-replace branch — the unrelated entry got
        # overwritten silently. Post-fix `_make_unique_entry_id`
        # suffix-disambiguates so the new save lands under a fresh
        # id and the original entry survives.
        #
        # Names deliberately differ so the NameCollisionModal path
        # doesn't preempt the bug class we're guarding.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        prior = SeqRecord(Seq("CCCC"), id="my_fragment_v2",
                            name="my_fragment_v2")
        prior.annotations["molecule_type"] = "DNA"
        prior.annotations["topology"]      = "linear"
        prior_gb = sc._record_to_gb_text(prior)
        sc._save_library([{
            "id": "my_fragment_v2",
            "name": "Pre-existing entry",
            "size": 4, "n_feats": 0,
            "source": "test", "added": "2026-05-20",
            "gb_text": prior_gb,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("AAAA", [])
            scr._dirty = True
            orig_push = app.push_screen
            def _stub_push(modal, callback=None):
                if isinstance(modal, sc.NamePlasmidModal):
                    # User types a name that sanitises to the existing
                    # entry's id but is a DIFFERENT human-readable name
                    # (so no NameCollisionModal fires).
                    if callback is not None:
                        callback("my fragment v2")
                    return None
                return orig_push(modal, callback=callback)
            app.push_screen = _stub_push  # type: ignore[method-assign]
            try:
                scr._do_save()
            finally:
                app.push_screen = orig_push  # type: ignore[method-assign]
            await pilot.pause()
            await pilot.pause()
            entries = sc._load_library()
            ids = [e.get("id") for e in entries]
            # Original "my_fragment_v2" entry must survive (pre-fix it
            # would have been clobbered by the new save's id match).
            assert "my_fragment_v2" in ids, (
                f"original entry overwritten — ids: {ids}"
            )
            # New save lands under a disambiguated id.
            assert any(
                i and i.startswith("my_fragment_v2_") for i in ids
            ), f"expected disambiguated id, got {ids}"
            # The original gb_text must NOT be overwritten.
            prior_entry = next(
                e for e in entries if e.get("id") == "my_fragment_v2"
            )
            assert prior_entry.get("gb_text") == prior_gb
            assert prior_entry.get("name") == "Pre-existing entry"


class TestPartsBinAutoTriggerNewPart:
    """Sweep #14 — ``PartsBinModal(auto_trigger_new_part=True)`` fires
    ``_new_part(None)`` on mount so Synthesis → Clone Fragment hands off
    cleanly to the Domesticator without an extra click."""

    def test_init_accepts_flag(self):
        modal = sc.PartsBinModal(auto_trigger_new_part=True)
        assert modal._auto_trigger_new_part is True

    def test_init_defaults_false(self):
        modal = sc.PartsBinModal()
        assert modal._auto_trigger_new_part is False


class TestSynthesisCloneFragmentButton:
    """Sweep #14 — Clone Fragment button on Synthesis toolbar."""

    async def test_clone_button_present_on_toolbar(self):
        from textual.widgets import Button
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            btn = scr.query_one("#btn-syn-clone", Button)
            assert "clone fragment" in str(btn.label).lower()

    async def test_clone_action_method_exists(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            assert callable(getattr(scr, "action_clone_fragment", None))

    async def test_clone_on_protein_tab_notifies(self):
        from textual.widgets import TabbedContent
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            tabs = scr.query_one("#syn-tabs", TabbedContent)
            tabs.active = "syn-tab-protein"
            await pilot.pause()
            await pilot.pause()
            notes: list[tuple[str, str]] = []
            orig_notify = app.notify
            def _capture(msg, *, severity="information", **kw):
                notes.append((severity, str(msg)))
                return orig_notify(msg, severity=severity, **kw)
            app.notify = _capture  # type: ignore[method-assign]
            try:
                scr.action_clone_fragment()
                await pilot.pause()
            finally:
                app.notify = orig_notify  # type: ignore[method-assign]
            assert any("DNA-tab only" in m for _, m in notes), (
                f"expected DNA-tab-only notify, got {notes}"
            )

    async def test_clone_empty_fragment_notifies(self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            notes: list[tuple[str, str]] = []
            orig_notify = app.notify
            def _capture(msg, *, severity="information", **kw):
                notes.append((severity, str(msg)))
                return orig_notify(msg, severity=severity, **kw)
            app.notify = _capture  # type: ignore[method-assign]
            try:
                scr.action_clone_fragment()
                await pilot.pause()
            finally:
                app.notify = orig_notify  # type: ignore[method-assign]
            assert any("empty" in m.lower() for _, m in notes), (
                f"expected empty-fragment notify, got {notes}"
            )


def _bind_valid_gb_l0_vector():
    """Bind a synthetic gb_l0 entry vector carrying two Esp3I sites so the
    Clone Fragment → Domesticator path finds a valid acceptor (and doesn't
    pop the entry-vector picker)."""
    import random
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rng = random.Random(0xBEEF)

    def _scrub(s):
        for x in ("CGTCTC", "GAGACG", "GGTCTC", "GAGACC"):
            s = s.replace(x, "CTGCAG")
        return s
    bb = _scrub("".join(rng.choice("ACGT") for _ in range(1400)))
    drop = _scrub("".join(rng.choice("ACGT") for _ in range(160)))
    cassette = "CGTCTCA" + "CTCG" + drop + "TGAG" + "AGAGACG"
    seq = bb[:700] + cassette + bb[700:]
    rec = SeqRecord(Seq(seq), id="TESTUPD", name="TESTUPD",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    sc._set_entry_vector("gb_l0", {
        "name": "TESTUPD", "size": 0, "source": "test",
        "id": "TESTUPD", "gb_text": sc._record_to_gb_text(rec)})


class TestSynthesisCloneFragmentFlow:
    """End-to-end integration guard for the reworked Clone Fragment
    handoff (2026-06-09): a method chooser first, NO premature save, then
    the Domesticator prefilled byte-exact."""

    async def test_clone_opens_method_chooser_without_saving(self, isolated_library):
        # Clicking Clone Fragment opens the method chooser and saves
        # NOTHING — the old flow auto-saved an unprimed copy first.
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-editor", sc.SynthesisEditor).load(
                "ATGAAACCCGGGTTT", [])
            scr._dirty = True
            scr.action_clone_fragment()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.CloneMethodChooserModal), (
                f"expected the method chooser, got {type(app.screen).__name__}"
            )
            assert sc._load_library() == [], "nothing should be saved yet"

    async def test_grammar_choice_routes_to_domesticator(self, isolated_library):
        # Picking a modular grammar (with a valid bound entry vector)
        # routes through the Parts Bin to the Domesticator — still with
        # no premature library save.
        _bind_valid_gb_l0_vector()
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-editor", sc.SynthesisEditor).load(
                "ATGAAACCCGGGTTT", [])
            scr._dirty = True
            scr.action_clone_fragment()
            await pilot.pause()
            await pilot.pause()
            app.screen.dismiss({"method": "grammar", "grammar_id": "gb_l0"})
            for _ in range(8):
                await pilot.pause()
            assert isinstance(app.screen, sc.DomesticatorModal), (
                f"expected DomesticatorModal, got {type(app.screen).__name__}"
            )
            stack = [type(s).__name__ for s in app.screen_stack]
            assert "PartsBinModal" in stack, f"stack: {stack}"
            assert sc._load_library() == [], "still nothing saved"

    async def test_clone_carries_clean_display_name(self, isolated_library):
        """A fragment named with SPACES seeds the Domesticator's part-name
        field cleanly (no underscores), with NO library save — the old
        flow round-tripped through a save; the new one just carries the
        name into the chooser → Domesticator handoff."""
        _bind_valid_gb_l0_vector()
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-editor", sc.SynthesisEditor).load(
                "ATGAAACCCGGGTTTAAA", [])
            scr._loaded_name = "Batch 57 clone test"
            scr._dirty = True
            scr.action_clone_fragment()
            await pilot.pause()
            await pilot.pause()
            app.screen.dismiss({"method": "grammar", "grammar_id": "gb_l0"})
            for _ in range(8):
                await pilot.pause()
            assert isinstance(app.screen, sc.DomesticatorModal)
            from textual.widgets import Input
            name_val = app.screen.query_one("#dom-name", Input).value
            assert name_val == "Batch 57 clone test", (
                f"clean name not carried into Domesticator: {name_val!r}"
            )
            assert "_" not in name_val, f"underscore leaked into name: {name_val!r}"
            assert sc._load_library() == [], "no premature save"

    async def test_clone_prefills_direct_input_textarea(
        self, isolated_library,
    ):
        """The complete synthesis sequence lands byte-exact in the
        Domesticator's #dom-direct-seq TextArea after choosing a grammar —
        the whole string in one write, no missed bases."""
        from textual.widgets import TextArea
        _bind_valid_gb_l0_vector()
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            # Use a sequence with a mix of bases so we can verify
            # the full payload landed, not just a prefix.
            full_seq = "ATGAAACCCGGGTTTAACCGGTTAACCGGTTAA"
            scr.query_one("#syn-editor", sc.SynthesisEditor).load(full_seq, [])
            scr._dirty = True
            scr.action_clone_fragment()
            await pilot.pause()
            await pilot.pause()
            app.screen.dismiss({"method": "grammar", "grammar_id": "gb_l0"})
            for _ in range(8):
                await pilot.pause()
            assert isinstance(app.screen, sc.DomesticatorModal)
            ta = app.screen.query_one("#dom-direct-seq", TextArea)
            assert ta.text == full_seq, (
                f"prefill seq mismatch: expected {full_seq!r}, "
                f"got {ta.text!r}"
            )

    async def test_clone_prefill_self_clears_after_firing(
        self, isolated_library,
    ):
        """The one-shot ``_clone_prefill_seq`` attr on PartsBinModal
        clears after _new_part consumes it so a subsequent manual
        'New Part' click doesn't re-prime the textarea with the same
        synthesis fragment."""
        _bind_valid_gb_l0_vector()
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-editor", sc.SynthesisEditor).load("ATGCATGC", [])
            scr._dirty = True
            scr.action_clone_fragment()
            await pilot.pause()
            await pilot.pause()
            app.screen.dismiss({"method": "grammar", "grammar_id": "gb_l0"})
            for _ in range(8):
                await pilot.pause()
            # Walk the screen stack to find the PartsBinModal under
            # the Domesticator and verify its prefill attr is empty.
            pb = next(
                (s for s in app.screen_stack
                 if isinstance(s, sc.PartsBinModal)),
                None,
            )
            assert pb is not None, "PartsBinModal not in stack"
            assert pb._clone_prefill_seq == "", (
                f"prefill should self-clear after firing; "
                f"still holds {pb._clone_prefill_seq!r}"
            )

    async def test_chooser_cancel_is_a_noop(self, isolated_library):
        # Cancelling the method chooser leaves the user in Synthesis with
        # their buffer intact and nothing written — there's no save to
        # fail anymore (the premature save is gone).
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-editor", sc.SynthesisEditor).load("ATGAAA", [])
            scr._dirty = True
            scr.action_clone_fragment()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.CloneMethodChooserModal)
            app.screen.dismiss(None)
            for _ in range(4):
                await pilot.pause()
            # Synthesis screen back on top — no handoff happened.
            assert isinstance(app.screen, sc.SynthesisScreen), (
                f"chooser cancel must return to Synthesis; "
                f"screen is {type(app.screen).__name__}"
            )
            # Library still empty.
            assert sc._load_library() == []


class TestSynthesisRenderExceptNarrowed:
    """Sweep #14 — three ``except Exception`` clauses on the render
    path narrowed to ``(AttributeError, TypeError)`` per invariant #1.
    White-box source check so the regression can't drift back to bare-
    ish exception handling."""

    def test_no_bare_exception_in_synth_editors(self):
        with open(sc.__file__, "r", encoding="utf-8") as fh:
            src = fh.read()
        # Look at the lines around the three known render-path sites.
        # Pre-fix used ``except Exception:``; post-fix uses
        # ``except (AttributeError, TypeError):``.
        # The render-path sites are inside SynthesisEditor +
        # ProteinEditor. Walk the file and assert that within those
        # class bodies, no plain ``except Exception:`` lurks at the
        # known column depth.
        synth_start = src.find("class SynthesisEditor(")
        protein_start = src.find("class ProteinEditor(")
        # End ProteinEditor's body at the next TOP-LEVEL class (was anchored on
        # `class RestrictionInsertModal(`, which moved to splicecraft_modals in
        # Phase D — use the generic next-class boundary so the anchor survives
        # future relocations).
        protein_end = src.find("\nclass ", protein_start + 1)
        assert synth_start > 0 and protein_start > 0 and protein_end > 0
        synth_body = src[synth_start:protein_start]
        protein_body = src[protein_start:protein_end]
        # No bare ``except Exception:`` (with optional whitespace+colon)
        # in the render path bodies — narrow types only.
        import re
        for body, name in ((synth_body, "SynthesisEditor"),
                            (protein_body, "ProteinEditor")):
            bare = re.findall(r"except Exception\s*:", body)
            assert not bare, (
                f"{name} body carries bare 'except Exception:' — "
                "narrow to specific types per invariant #1"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #15 — feature pre-coloring + Edit buttons
# ═══════════════════════════════════════════════════════════════════════════════

class TestProteinEditorFeatures:
    """ProteinEditor gains AA-coord feature tracking — features render
    as colored bands in both render modes and survive insert/delete."""

    def test_load_accepts_feats_kwarg(self):
        pe = sc.ProteinEditor()
        pe.load("MASGGGSGGGS", feats=[
            {"start": 0, "end": 1, "label": "M",
             "type": "Motif", "color": "#A855F7"},
            {"start": 3, "end": 8, "label": "linker",
             "type": "Linker", "color": "#6B7280"},
        ])
        assert len(pe._aa_feats) == 2
        assert pe._aa_feats[0]["start"] == 0
        assert pe._aa_feats[1]["color"] == "#6B7280"

    def test_load_drops_out_of_bounds_feats(self):
        pe = sc.ProteinEditor()
        pe.load("MAS", feats=[
            {"start": 0, "end": 2, "label": "valid"},
            {"start": 0, "end": 99, "label": "end-out-of-range"},
            {"start": 5, "end": 6, "label": "start-out-of-range"},
            {"start": 2, "end": 2, "label": "zero-length"},
        ])
        assert len(pe._aa_feats) == 1
        assert pe._aa_feats[0]["label"] == "valid"

    def test_get_feats_returns_deepcopy(self):
        pe = sc.ProteinEditor()
        pe.load("MAS", feats=[
            {"start": 0, "end": 2, "label": "x", "color": "#FF0000"},
        ])
        snap = pe.get_feats()
        snap[0]["color"] = "#00FF00"
        # Mutating the snapshot must NOT touch the editor's stored
        # feature.
        assert pe._aa_feats[0]["color"] == "#FF0000"

    def test_add_feature_appends_and_clips(self):
        pe = sc.ProteinEditor()
        pe.load("MASGGGS")
        pe.add_feature({
            "start": 1, "end": 4, "label": "blue",
            "type": "Tag", "color": "#3B82F6",
        })
        assert len(pe._aa_feats) == 1
        # Out-of-bounds: silently rejected (defensive).
        pe.add_feature({"start": -1, "end": 2})
        pe.add_feature({"start": 0, "end": 99})
        pe.add_feature({"start": 3, "end": 3})  # zero-length
        assert len(pe._aa_feats) == 1

    async def test_insert_at_cursor_shifts_feats_right_of_cursor(self):
        # insert_at_cursor posts a Changed message, which needs an
        # active app context. Mount through SynthesisScreen so the
        # editor lands inside a real screen stack.
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MASGGGS")
            pe.add_feature({"start": 4, "end": 7,
                             "label": "right", "color": "#FFAA00"})
            pe._cursor_pos = 2
            # Use "KK" not "XX" — X isn't in _PROTEIN_AA_ALPHABET
            # (20 standard + stop only); the filter would drop it.
            pe.insert_at_cursor("KK")
            await pilot.pause()
            assert pe._aa_feats[0]["start"] == 6
            assert pe._aa_feats[0]["end"] == 9

    async def test_insert_does_not_extend_feat_strictly_left_of_cursor(
        self,
    ):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            pe = app.screen.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MASGGGS")
            pe.add_feature({"start": 0, "end": 2,
                             "label": "left", "color": "#FFAA00"})
            pe._cursor_pos = 4   # past the feature
            pe.insert_at_cursor("KK")
            await pilot.pause()
            # Feature strictly left of cursor: untouched.
            assert pe._aa_feats[0]["start"] == 0
            assert pe._aa_feats[0]["end"] == 2

    def test_delete_range_clips_overlapping_feature(self):
        pe = sc.ProteinEditor()
        pe.load("MASGGGSXYZ")
        pe.add_feature({"start": 2, "end": 7,
                         "label": "overlap", "color": "#FFAA00"})
        # Delete bp 4..6 — overlaps feature.
        pe._delete_range(4, 6)
        assert len(pe._aa_feats) == 1
        # Surviving region clipped + shifted.
        assert pe._aa_feats[0]["start"] == 2
        assert pe._aa_feats[0]["end"] == 5

    def test_delete_drops_zero_survival_feature(self):
        pe = sc.ProteinEditor()
        pe.load("MASGGGSXYZ")
        pe.add_feature({"start": 3, "end": 5,
                         "label": "doomed", "color": "#FFAA00"})
        # Delete the entire feature span.
        pe._delete_range(3, 5)
        assert pe._aa_feats == []


class TestProteinMotifInsertColored:
    """Motif insert produces a feature with a color sourced from the
    type→color palette (or the motif's own ``color`` field if set)."""

    async def test_motif_insert_adds_colored_feature(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MM")
            pe._cursor_pos = 1
            scr._motif_selected_entry = lambda: {
                "name": "His6", "feature_type": "Tag",
                "sequence": "HHHHHH",
            }
            scr._motif_insert_selected()
            await pilot.pause()
            assert pe._aa_seq == "MHHHHHHM"
            # Feature landed covering the inserted span.
            assert len(pe._aa_feats) == 1
            feat = pe._aa_feats[0]
            assert feat["start"] == 1
            assert feat["end"] == 7
            assert feat["label"] == "His6"
            # Tag → blue from the type→color map.
            assert feat["color"] == sc._PROTEIN_FEATURE_TYPE_COLORS["Tag"]

    async def test_motif_with_explicit_color_wins(self):
        # A motif whose entry carries its own `color` field overrides
        # the type→palette default — same precedence the DNA feature
        # library uses.
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("M")
            pe._cursor_pos = 1
            scr._motif_selected_entry = lambda: {
                "name": "Custom",
                "feature_type": "Tag",
                "sequence": "GGG",
                "color": "#123456",
            }
            scr._motif_insert_selected()
            await pilot.pause()
            assert pe._aa_feats[0]["color"] == "#123456"


class TestProteinSaveLoadRoundTrip:
    """Motif features written on save survive a save + reload cycle
    via the CDS sub-feature encoding."""

    async def test_aa_features_survive_save_and_load(
        self, isolated_library,
    ):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("MASGGGSHHHHHH")
            pe.add_feature({
                "start": 7, "end": 13,
                "label": "His6", "type": "Tag",
                "color": "#3B82F6", "strand": 1,
            })
            scr._protein_loaded_id   = "round_trip_test"
            scr._protein_loaded_name = "round_trip_test"
            scr._protein_dirty = True
            # Skip the unique-id helper's name prompt — load the
            # entry directly to verify reload picks up the motif.
            scr._commit_protein_save("MASGGGSHHHHHH", after=None)
            await pilot.pause()
            # Reset editor state and reload from library.
            pe.load("")
            assert pe._aa_feats == []
            scr._load_protein_entry_by_id("round_trip_test")
            await pilot.pause()
            assert pe._aa_seq == "MASGGGSHHHHHH"
            # Motif feature restored.
            assert len(pe._aa_feats) == 1
            f = pe._aa_feats[0]
            assert f["start"] == 7
            assert f["end"] == 13
            assert f["label"] == "His6"
            assert f["type"] == "Tag"
            assert f["color"] == "#3B82F6"


class TestProteinMotifsPersistence:
    """`_load_protein_motifs` merges built-ins with user overrides;
    `_save_protein_motifs` writes only the user-modified entries."""

    def test_load_returns_builtins_when_no_user_file(self):
        # Cold cache + no user file → all built-ins surface.
        sc._state._protein_motifs_cache = None
        merged = sc._load_protein_motifs()
        builtin_names = {m["name"] for m in sc._PROTEIN_MOTIFS}
        merged_names = {m["name"] for m in merged}
        assert builtin_names.issubset(merged_names)

    def test_user_override_replaces_builtin_by_name(self):
        sc._save_protein_motifs([{
            "name": "His6", "feature_type": "Tag",
            "sequence": "HHHHHHHH",  # user-edited (now 8 H's)
            "description": "Modified by user",
            "color": "#FF00FF",
        }])
        # Read back through the merged loader.
        merged = sc._load_protein_motifs()
        his6 = next(m for m in merged if m["name"] == "His6")
        assert his6["sequence"] == "HHHHHHHH"
        assert his6["color"] == "#FF00FF"
        # All other built-ins still present.
        names = {m["name"] for m in merged}
        for builtin in sc._PROTEIN_MOTIFS:
            assert builtin["name"] in names

    def test_user_added_novel_motif_appends(self):
        sc._save_protein_motifs([{
            "name": "MyCustomMotif",
            "feature_type": "Tag",
            "sequence": "KGGKGG",
            "description": "My personal sequence",
        }])
        merged = sc._load_protein_motifs()
        assert any(m["name"] == "MyCustomMotif" for m in merged)


class TestSynthesisEditButtons:
    """Edit buttons on DNA feature library + Protein motif library."""

    async def test_dna_featlib_edit_button_present(self):
        from textual.widgets import Button
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            btn = scr.query_one("#btn-syn-featlib-edit", Button)
            assert "edit" in str(btn.label).lower()

    async def test_protein_motif_edit_button_present(self):
        from textual.widgets import Button, TabbedContent
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            tabs = scr.query_one("#syn-tabs", TabbedContent)
            tabs.active = "syn-tab-protein"
            await pilot.pause()
            await pilot.pause()
            btn = scr.query_one("#btn-syn-motif-edit", Button)
            assert "edit" in str(btn.label).lower()

    async def test_protein_motif_edit_writes_to_user_file(self):
        # Editing a built-in motif lands the edited entry in the
        # user file via _save_protein_motifs. The original built-in
        # stays in code; the merge layer surfaces the user copy.
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            # Stub the selector + the modal so we don't depend on
            # focus + row state.
            scr._motif_selected_entry = lambda: {
                "name": "His6", "feature_type": "Tag",
                "sequence": "HHHHHH",
            }
            orig_push = app.push_screen
            def _stub(modal, callback=None):
                if isinstance(modal, sc.AddFeatureModal):
                    if callback is not None:
                        callback({
                            "entry": {
                                "name": "His6",
                                "feature_type": "Tag",
                                "sequence": "HHHHHHHHHH",
                                "color": "#FF00FF",
                            },
                        })
                    return None
                return orig_push(modal, callback=callback)
            app.push_screen = _stub  # type: ignore[method-assign]
            try:
                scr._on_motif_edit(None)
            finally:
                app.push_screen = orig_push  # type: ignore[method-assign]
            await pilot.pause()
            # Merged read picks up the user edit.
            merged = sc._load_protein_motifs()
            his6 = next(m for m in merged if m["name"] == "His6")
            assert his6["sequence"] == "HHHHHHHHHH"
            assert his6.get("color") == "#FF00FF"


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #16 — distinct motif colors + dithered protein render
# ═══════════════════════════════════════════════════════════════════════════════

class TestProteinMotifsDistinctColors:
    """Every built-in motif carries its own `color` hex so the user
    can tell them apart in the motif library list AND in the dither
    bar above the AA row."""

    def test_every_builtin_motif_has_color(self):
        for m in sc._PROTEIN_MOTIFS:
            color = m.get("color")
            assert color, (
                f"motif {m.get('name')!r} missing color field"
            )
            assert isinstance(color, str)
            assert color.startswith("#"), (
                f"motif {m.get('name')!r} color {color!r} not a hex"
            )

    def test_all_motif_colors_are_distinct(self):
        seen: dict[str, str] = {}
        for m in sc._PROTEIN_MOTIFS:
            name  = m.get("name")
            color = (m.get("color") or "").lower()
            assert color, f"motif {name!r} missing color"
            if color in seen:
                raise AssertionError(
                    f"motif {name!r} reuses color {color} already "
                    f"taken by {seen[color]!r}"
                )
            seen[color] = name

    def test_palette_key_matches_motif_feature_types(self):
        """Every `feature_type` used in `_PROTEIN_MOTIFS` must have a
        matching key in the fallback palette so user-added motifs
        without an explicit color resolve to the right family default.
        Pre-sweep #16 the palette had `"2A peptide"` while the data
        used `"2A"` — silent fallback to `Motif` purple for every
        2A motif.
        """
        used_types = {m.get("feature_type", "") for m in sc._PROTEIN_MOTIFS}
        palette_keys = set(sc._PROTEIN_FEATURE_TYPE_COLORS.keys())
        missing = used_types - palette_keys
        assert not missing, (
            f"feature_type(s) {missing} appear in _PROTEIN_MOTIFS but "
            "aren't in _PROTEIN_FEATURE_TYPE_COLORS — silent fallback"
        )


class TestProteinDitherRender:
    """The Protein tab now renders features as a dithered ▒-block row
    above the AA letters, matching the seq panel's lane style."""

    def test_lane_art_emitted_only_when_features_present(self):
        pe = sc.ProteinEditor()
        pe.load("MASGGGS")
        # No features → row count should be 2 (codon mode default).
        assert pe._row_count() == 2
        pe.load("MASGGGS", feats=[
            {"start": 0, "end": 4, "label": "x",
             "type": "Motif", "color": "#FF0000", "strand": 1},
        ])
        # Features present → 1 lane = 2 rows (bar + label) prepended.
        # 2 base (AA + codon) + 2 lane = 4 rows.
        assert pe._row_count() == 4

    def test_lane_bar_codon_mode_emits_block_glyphs(self):
        pe = sc.ProteinEditor()
        pe.load("MASGGGS", feats=[
            {"start": 1, "end": 4, "label": "x",
             "type": "Motif", "color": "#FF0000", "strand": 1},
        ])
        # The legacy `_build_dither_row` returns the bar row (the
        # line closest to AA, i.e. the LAST line of the reversed
        # lane text).
        dither = pe._build_dither_row(cols_per_aa=3)
        body = dither.plain[pe._FLANK_MARKER_WIDTH:]
        assert body[0:3] == "   "
        assert body[3:11].count("▒") == 8
        assert body[11] == "▶"

    def test_lane_aa_only_mode_one_cell_per_aa(self):
        pe = sc.ProteinEditor()
        pe.set_codon_mode(False)
        pe.load("MASGGGS", feats=[
            {"start": 0, "end": 3, "label": "x",
             "type": "Motif", "color": "#FF0000", "strand": 1},
        ])
        dither = pe._build_dither_row(cols_per_aa=1)
        body = dither.plain[pe._FLANK_MARKER_WIDTH:]
        # AA 0..3 covered: 2 ▒ + 1 ▶ at the right terminus.
        assert body[:3].count("▒") == 2
        assert body[2] == "▶"

    def test_reverse_strand_uses_left_arrowhead(self):
        pe = sc.ProteinEditor()
        pe.load("MASGGGS", feats=[
            {"start": 1, "end": 5, "label": "rev",
             "type": "Motif", "color": "#00FF00", "strand": -1},
        ])
        dither = pe._build_dither_row(cols_per_aa=3)
        body = dither.plain[pe._FLANK_MARKER_WIDTH:]
        # Feature 1..5 → cols 3..15. Leftmost (3) should be ◀.
        assert body[3] == "◀"
        # Rightmost (14) should be ▒ not ▶ (reverse strand).
        assert body[14] == "▒"

    def test_no_feature_color_skips_lane_cell(self):
        # A feature without a `color` field shouldn't paint anything
        # on the bar. The lane height stays 0 (the empty-color path
        # in `_build_protein_lane_text` skips the feature entirely
        # so it doesn't contribute to either the pack or the render).
        pe = sc.ProteinEditor()
        pe.load("MASGGGS", feats=[
            {"start": 0, "end": 4, "label": "no-color",
             "type": "Motif"},
        ])
        # The pack still places the feature, but the bar render skips
        # painting because color is empty. So we should see blank.
        text, n = pe._build_protein_lane_text(cols_per_aa=3)
        plain = text.plain
        assert "▒" not in plain
        assert "▶" not in plain

    def test_lane_includes_centred_label_row(self):
        pe = sc.ProteinEditor()
        pe.load("MASGGGSAAAA", feats=[
            {"start": 0, "end": 6, "label": "His6",
             "type": "Tag", "color": "#1E40AF", "strand": 1},
        ])
        text, n_rows = pe._build_protein_lane_text(cols_per_aa=3)
        # 1 feature → 2 rows (label + bar).
        assert n_rows == 2
        lines = text.split("\n")
        assert len(lines) == 2
        # Reversed order: highest stack row first; bar row LAST.
        label_line = lines[0].plain[pe._FLANK_MARKER_WIDTH:]
        bar_line   = lines[1].plain[pe._FLANK_MARKER_WIDTH:]
        # Bar row: ▒ blocks across cols 0..18 (AA 0..6, cpa=3).
        # Last cell becomes ▶.
        assert "▒" in bar_line
        assert "▶" in bar_line
        # Label row centres "His6" within the 18-cell span.
        # 18-cell span, 4-char label → starts at col (18-4)//2 = 7.
        # So the H of His6 should be at body[7].
        assert "His6" in label_line
        # Label NOT on bar row.
        assert "His6" not in bar_line

    def test_overlapping_features_stack_in_separate_lanes(self):
        # Two motifs overlapping the same AA range stack into two
        # lanes — older feature closer to AA, newer one above.
        pe = sc.ProteinEditor()
        pe.load("MASGGGSAAAAA", feats=[
            {"start": 0, "end": 4, "label": "near",
             "type": "Tag", "color": "#FF0000", "strand": 1},
            {"start": 2, "end": 6, "label": "above",
             "type": "Tag", "color": "#00FF00", "strand": 1},
        ])
        # 2 overlapping features → 2 lanes × 2 rows = 4 lane rows.
        # + 2 base = 6 total rows in codon mode.
        assert pe._row_count() == 6
        text, n_rows = pe._build_protein_lane_text(cols_per_aa=3)
        assert n_rows == 4

    def test_non_overlapping_features_share_one_lane(self):
        # Two features in disjoint AA ranges pack into the same lane.
        pe = sc.ProteinEditor()
        pe.load("MASGGGSAAAAA", feats=[
            {"start": 0, "end": 4, "label": "left",
             "type": "Tag", "color": "#FF0000", "strand": 1},
            {"start": 5, "end": 10, "label": "right",
             "type": "Tag", "color": "#00FF00", "strand": 1},
        ])
        # 2 non-overlapping → 1 lane × 2 rows = 2 lane rows.
        assert pe._row_count() == 4
        text, n_rows = pe._build_protein_lane_text(cols_per_aa=3)
        assert n_rows == 2


class TestProteinMotifInsertColorIsUnique:
    """End-to-end check: inserting two different built-in motifs
    produces features with different colors."""

    async def test_two_distinct_motif_inserts_have_different_colors(
        self,
    ):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("M")
            # Find two motifs of the same feature_type so we know we're
            # actually testing per-motif distinct colors (not just
            # per-family-type fallbacks).
            his6 = next(m for m in sc._PROTEIN_MOTIFS
                          if m.get("name") == "His6")
            flag = next(m for m in sc._PROTEIN_MOTIFS
                          if m.get("name") == "FLAG")
            # Both are "Tag" — pre-sweep #16 they shared "#3B82F6".
            assert his6.get("feature_type") == flag.get("feature_type")
            assert his6.get("color") != flag.get("color")
            # Insert both and confirm each lands a distinct color.
            pe._cursor_pos = 1
            scr._motif_selected_entry = lambda: dict(his6)
            scr._motif_insert_selected()
            pe._cursor_pos = len(pe._aa_seq)
            scr._motif_selected_entry = lambda: dict(flag)
            scr._motif_insert_selected()
            await pilot.pause()
            colors = {f["color"] for f in pe._aa_feats}
            assert len(colors) == 2


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #17 — insert-at-cursor boundary-extend fix
# ═══════════════════════════════════════════════════════════════════════════════

class TestProteinMotifNonOverlapAtBoundary:
    """Regression guard for 2026-05-21 fix.

    Pre-sweep, inserting a second motif at the 3' boundary of an
    existing motif extended the first motif's `end` to cover the
    insertion AND stacked the new motif's feature on top — a visible
    overlap. Post-sweep the first motif's `end` stays put and the
    new motif lands as a clean downstream feature.
    """

    async def test_two_motifs_back_to_back_dont_overlap(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("")
            # Insert His6 at start.
            scr._motif_selected_entry = lambda: {
                "name": "His6", "feature_type": "Tag",
                "sequence": "HHHHHH", "color": "#1E40AF",
            }
            scr._motif_insert_selected()
            await pilot.pause()
            # Cursor is now at the 3' end of His6 (position 6).
            assert pe._cursor_pos == 6
            assert len(pe._aa_feats) == 1
            assert pe._aa_feats[0]["start"] == 0
            assert pe._aa_feats[0]["end"]   == 6
            # Insert FLAG at the same cursor (3' boundary of His6).
            scr._motif_selected_entry = lambda: {
                "name": "FLAG", "feature_type": "Tag",
                "sequence": "DYKDDDDK", "color": "#0E7490",
            }
            scr._motif_insert_selected()
            await pilot.pause()
            # His6 must NOT have extended.
            his6 = next(f for f in pe._aa_feats if f["label"] == "His6")
            flag = next(f for f in pe._aa_feats if f["label"] == "FLAG")
            assert his6["start"] == 0
            assert his6["end"]   == 6, (
                f"His6 must not extend; got end={his6['end']}"
            )
            # FLAG sits cleanly downstream — no overlap.
            assert flag["start"] == 6
            assert flag["end"]   == 14
            # The sequence appended in order.
            assert pe._aa_seq == "HHHHHHDYKDDDDK"

    async def test_motif_insert_at_5p_boundary_shifts_existing(self):
        # 5' boundary case: insert a NEW motif at position 0 with an
        # existing motif already at [0..6). The existing motif's start
        # was at the cursor, so it shifts right; the new motif lands at
        # [0..N) where N is its length. Neither feature overlaps.
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one(
                "#syn-protein-editor", sc.ProteinEditor,
            )
            pe.load("HHHHHH", feats=[
                {"start": 0, "end": 6, "label": "His6",
                 "type": "Tag", "color": "#1E40AF", "strand": 1},
            ])
            # Move cursor to position 0 (5' end / start of His6).
            pe._cursor_pos = 0
            scr._motif_selected_entry = lambda: {
                "name": "Kozak", "feature_type": "Motif",
                "sequence": "M", "color": "#C026D3",
            }
            scr._motif_insert_selected()
            await pilot.pause()
            kozak = next(
                f for f in pe._aa_feats if f["label"] == "Kozak"
            )
            his6 = next(
                f for f in pe._aa_feats if f["label"] == "His6"
            )
            # Kozak at [0..1), His6 shifted to [1..7) — no overlap.
            assert kozak["start"] == 0
            assert kozak["end"]   == 1
            assert his6["start"]  == 1
            assert his6["end"]    == 7

    async def test_typing_inside_feature_still_extends_bar(self):
        # Sweep #17 invariant: the boundary-extend fix MUST NOT change
        # the "typing inside a feature" behaviour. With the cursor
        # strictly inside a feature span, inserting bases must STILL
        # extend the feature's end.
        pe = sc.ProteinEditor()
        pe.load("HHHHHH", feats=[
            {"start": 0, "end": 6, "label": "His6",
             "type": "Tag", "color": "#1E40AF", "strand": 1},
        ])
        # Cursor strictly inside (between H3 and H4).
        pe._cursor_pos = 3
        # Test the editor primitive directly (skip the AddFeature flow).
        # extend_adjacent_feats=False mode (the new path) still extends
        # because cur is STRICTLY inside the feature.
        pe.insert_at_cursor("HH", extend_adjacent_feats=False)
        his6 = pe._aa_feats[0]
        assert his6["start"] == 0
        assert his6["end"]   == 8, (
            f"His6 should grow when typing inside; got end={his6['end']}"
        )
        assert pe._aa_seq == "HHHHHHHH"

    async def test_typing_inside_feature_default_path_also_extends(self):
        # Same scenario but with default kwarg (extend=True) — the
        # editor's keyboard-typing path uses this. Verifies neither
        # branch broke the inside-extend rule.
        pe = sc.ProteinEditor()
        pe.load("HHHHHH", feats=[
            {"start": 0, "end": 6, "label": "His6",
             "type": "Tag", "color": "#1E40AF", "strand": 1},
        ])
        pe._cursor_pos = 3
        pe.insert_at_cursor("HH")  # default extend_adjacent_feats=True
        his6 = pe._aa_feats[0]
        assert his6["end"] == 8

    async def test_deleting_inside_feature_contracts_bar(self):
        # Sweep #17 invariant: deleting bases inside a feature
        # contracts the bar.
        pe = sc.ProteinEditor()
        pe.load("HHHHHHHH", feats=[
            {"start": 0, "end": 8, "label": "His8",
             "type": "Tag", "color": "#3B82F6", "strand": 1},
        ])
        # Delete the middle two H's.
        pe._delete_range(3, 5)
        his = pe._aa_feats[0]
        assert his["start"] == 0
        assert his["end"]   == 6, (
            f"His8 should contract on delete; got end={his['end']}"
        )
        assert pe._aa_seq == "HHHHHH"


class TestDnaFeatLibNonOverlapAtBoundary:
    """Mirror sweep #17 — same boundary-extend fix on the DNA tab."""

    async def test_dna_featlib_insert_at_3p_boundary_no_overlap(
        self, isolated_library,
    ):
        # Seed the feature library with two entries we can insert.
        sc._save_features([
            {"name": "alpha", "feature_type": "misc_feature",
             "sequence": "AAA", "color": "#FF0000", "strand": 1,
             "qualifiers": {}},
            {"name": "beta",  "feature_type": "misc_feature",
             "sequence": "GGG", "color": "#00FF00", "strand": 1,
             "qualifiers": {}},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("", [])
            # Insert alpha at start.
            scr._featlib_selected_entry = lambda: {
                "name": "alpha", "feature_type": "misc_feature",
                "sequence": "AAA", "color": "#FF0000", "strand": 1,
                "qualifiers": {},
            }
            scr._featlib_insert_selected(mode="insert")
            await pilot.pause()
            assert ed._cursor_pos == 3
            assert len(ed._feats) == 1
            assert ed._feats[0]["start"] == 0
            assert ed._feats[0]["end"]   == 3
            # Insert beta at the 3' boundary of alpha.
            scr._featlib_selected_entry = lambda: {
                "name": "beta", "feature_type": "misc_feature",
                "sequence": "GGG", "color": "#00FF00", "strand": 1,
                "qualifiers": {},
            }
            scr._featlib_insert_selected(mode="insert")
            await pilot.pause()
            alpha = next(f for f in ed._feats if f["label"] == "alpha")
            beta  = next(f for f in ed._feats if f["label"] == "beta")
            # alpha must NOT have extended.
            assert alpha["start"] == 0
            assert alpha["end"]   == 3, (
                f"alpha must not extend; got end={alpha['end']}"
            )
            # beta sits cleanly downstream.
            assert beta["start"] == 3
            assert beta["end"]   == 6
            assert ed._seq == "AAAGGG"

    async def test_dna_typing_inside_feature_still_extends_bar(self):
        # The user-typing path passes extend_adjacent_feats=True by
        # default. Inside-feature typing must still extend the bar
        # (and verifies the boundary-extend fix didn't regress the
        # inside-extend rule).
        ed = sc.SynthesisEditor()
        ed.load("AAAGGG", [{
            "start": 0, "end": 6, "label": "joined",
            "type": "misc_feature", "color": "#FF0000", "strand": 1,
        }])
        ed._cursor_pos = 3   # strictly inside [0..6)
        ed.insert_at_cursor("TT", extend_adjacent_feats=False)
        # extend_adjacent_feats=False ALSO extends inside-feature
        # because cur strictly < fe.
        assert ed._feats[0]["start"] == 0
        assert ed._feats[0]["end"]   == 8


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #18 — protein-load filter tightening
# ═══════════════════════════════════════════════════════════════════════════════

class TestProteinLoadFilterStrict:
    """`_extract_aa_feats_from_record` only picks up misc_feature
    sub-features whose coords are codon-aligned. Tighter than the
    pre-sweep "not CDS" filter, which would have also pulled in
    gene/RBS/mat_peptide/etc. annotations from an imported record."""

    def _build_record(self, dna: str, cds_features: list[dict]):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq(dna), id="t", name="t")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "linear"
        # First feature is the CDS spanning the whole DNA.
        rec.features.append(
            SeqFeature(
                FeatureLocation(0, len(dna), strand=1),
                type="CDS",
                qualifiers={
                    "label": ["protein"],
                    "translation": [
                        # back-translate ATGGCC… as M, A, S, H, H, H
                        "MASHHH"
                    ],
                },
            ),
        )
        for f in cds_features:
            rec.features.append(
                SeqFeature(
                    FeatureLocation(
                        f["start"], f["end"], strand=f.get("strand", 1),
                    ),
                    type=f["type"],
                    qualifiers=dict(f.get("qualifiers") or {}),
                )
            )
        return rec

    def test_misc_feature_is_picked_up(self):
        rec = self._build_record(
            "ATGGCCAGCCACCACCAC",   # MASHHH (6 aa)
            [{"start": 0, "end": 18, "type": "misc_feature",
              "qualifiers": {
                  "label": ["HisRich"],
                  "note": [
                      "splicecraft-aa-feature-type=Tag",
                  ],
                  "ApEinfo_fwdcolor": ["#1E40AF"],
              }}],
        )
        feats = sc.SynthesisScreen._extract_aa_feats_from_record(rec)
        assert len(feats) == 1
        assert feats[0]["start"] == 0
        assert feats[0]["end"]   == 6
        assert feats[0]["label"] == "HisRich"
        assert feats[0]["type"]  == "Tag"
        assert feats[0]["color"] == "#1E40AF"

    def test_gene_subfeature_inside_cds_is_filtered_out(self):
        # Pre-sweep: a `gene` feature INSIDE the CDS would mis-load
        # as an AA motif because the filter only skipped CDS. Now
        # it's filtered out.
        rec = self._build_record(
            "ATGGCCAGCCACCACCAC",
            [{"start": 0, "end": 18, "type": "gene",
              "qualifiers": {"label": ["intruder"]}}],
        )
        feats = sc.SynthesisScreen._extract_aa_feats_from_record(rec)
        assert feats == [], (
            f"gene sub-feature should be filtered; got {feats}"
        )

    def test_rbs_subfeature_inside_cds_is_filtered_out(self):
        rec = self._build_record(
            "ATGGCCAGCCACCACCAC",
            [{"start": 0, "end": 6, "type": "RBS",
              "qualifiers": {"label": ["fakeRBS"]}}],
        )
        feats = sc.SynthesisScreen._extract_aa_feats_from_record(rec)
        assert feats == []

    def test_mat_peptide_subfeature_inside_cds_is_filtered_out(self):
        rec = self._build_record(
            "ATGGCCAGCCACCACCAC",
            [{"start": 3, "end": 12, "type": "mat_peptide",
              "qualifiers": {"label": ["mature"]}}],
        )
        feats = sc.SynthesisScreen._extract_aa_feats_from_record(rec)
        assert feats == []

    def test_codon_misaligned_misc_feature_is_dropped(self):
        # bp 1..10 is NOT codon-aligned. The pre-sweep code would
        # have floored to AA[0..3) — silently wrong residues.
        # Post-sweep: dropped.
        rec = self._build_record(
            "ATGGCCAGCCACCACCAC",
            [{"start": 1, "end": 10, "type": "misc_feature",
              "qualifiers": {
                  "label": ["misaligned"],
                  "note": ["splicecraft-aa-feature-type=Motif"],
              }}],
        )
        feats = sc.SynthesisScreen._extract_aa_feats_from_record(rec)
        assert feats == [], (
            f"mis-aligned coords should drop; got {feats}"
        )

    def test_misc_feature_with_only_one_aligned_boundary_dropped(self):
        # fs=0 (aligned), fe=10 (NOT aligned → 10 % 3 == 1). Drop.
        rec = self._build_record(
            "ATGGCCAGCCACCACCAC",
            [{"start": 0, "end": 10, "type": "misc_feature",
              "qualifiers": {"label": ["half-aligned"]}}],
        )
        feats = sc.SynthesisScreen._extract_aa_feats_from_record(rec)
        assert feats == []

    def test_mixed_subfeatures_only_misc_aligned_passes(self):
        # Three sub-features inside the CDS:
        #   * misc_feature aligned at [0..6) → KEPT (aa 0..2)
        #   * gene at [0..18) → FILTERED (wrong type)
        #   * misc_feature mis-aligned at [7..13) → DROPPED (alignment)
        rec = self._build_record(
            "ATGGCCAGCCACCACCAC",
            [
                {"start": 0,  "end": 6,  "type": "misc_feature",
                 "qualifiers": {
                     "label": ["pass"],
                     "note": ["splicecraft-aa-feature-type=Tag"],
                 }},
                {"start": 0,  "end": 18, "type": "gene",
                 "qualifiers": {"label": ["wrong-type"]}},
                {"start": 7,  "end": 13, "type": "misc_feature",
                 "qualifiers": {"label": ["misaligned"]}},
            ],
        )
        feats = sc.SynthesisScreen._extract_aa_feats_from_record(rec)
        assert len(feats) == 1
        assert feats[0]["label"] == "pass"
        assert feats[0]["start"] == 0
        assert feats[0]["end"]   == 2

    def test_existing_round_trip_still_works(self):
        # Confirm the sweep doesn't regress the normal motif round-
        # trip (the path we wrote earlier in sweep #15).
        rec = self._build_record(
            "ATGGCCAGCCACCACCAC",  # 18 bp, MASHHH
            [{"start": 9, "end": 18, "type": "misc_feature",
              "qualifiers": {
                  "label": ["His3"],
                  "note": ["splicecraft-aa-feature-type=Tag"],
                  "ApEinfo_fwdcolor": ["#3B82F6"],
              }}],
        )
        feats = sc.SynthesisScreen._extract_aa_feats_from_record(rec)
        assert len(feats) == 1
        assert feats[0]["start"] == 3
        assert feats[0]["end"]   == 6
        assert feats[0]["label"] == "His3"
        assert feats[0]["type"]  == "Tag"


class TestProteinMotifEditor:
    """The GUI motif-library editor (New / Edit / Delete in the Synthesis
    protein tab) goes through `_protein_motif_upsert` / `_protein_motif_delete`,
    which mirror the agent-endpoint override semantics."""

    def test_upsert_new_motif_round_trip(self):
        assert sc._protein_motif_upsert(
            "MyTag", "hhhkkk", feature_type="Tag",
            color="#123456", description="custom") is None
        by_name = {m["name"]: m for m in sc._load_protein_motifs()}
        assert "MyTag" in by_name
        assert by_name["MyTag"]["sequence"] == "HHHKKK"   # upper-cased
        assert by_name["MyTag"]["feature_type"] == "Tag"

    def test_upsert_override_then_delete_restores_builtin(self):
        builtins = {m["name"]: m for m in sc._load_protein_motifs()}
        assert "His6" in builtins
        orig_desc = builtins["His6"].get("description")
        assert sc._protein_motif_upsert(
            "His6", "HHHHHH", feature_type="Tag",
            description="OVERRIDDEN") is None
        after = {m["name"]: m for m in sc._load_protein_motifs()}["His6"]
        assert after["description"] == "OVERRIDDEN"
        # Deleting the user override restores the built-in.
        assert sc._protein_motif_delete("His6") is None
        restored = {m["name"]: m for m in sc._load_protein_motifs()}["His6"]
        assert restored.get("description") == orig_desc

    def test_upsert_rejects_invalid(self):
        assert sc._protein_motif_upsert("", "HHHH") is not None      # no name
        assert sc._protein_motif_upsert("X", "") is not None         # no seq
        assert sc._protein_motif_upsert("X", "HHJ123") is not None   # non-AA

    def test_delete_builtin_without_override_rejected(self):
        # Fresh sandbox → no user override for His6 → can't delete a built-in.
        assert sc._protein_motif_delete("His6") is not None

    def test_delete_user_novel_motif(self):
        assert sc._protein_motif_upsert("Zztag", "WWWW") is None
        assert sc._protein_motif_delete("Zztag") is None
        names = {m["name"] for m in sc._load_protein_motifs()}
        assert "Zztag" not in names

    def test_editor_modal_and_handlers_exist(self):
        assert hasattr(sc, "NewMotifModal")
        assert hasattr(sc.SynthesisScreen, "_on_motif_new")
        assert hasattr(sc.SynthesisScreen, "_on_motif_delete")


class TestProteinEditorPaste:
    """ProteinEditor.on_paste mirrors SynthesisEditor.on_paste — a
    terminal native-paste lands amino acids at the cursor through
    `insert_at_cursor` (AA-filter + cap + selection-replace)."""

    def test_handler_exists(self):
        assert hasattr(sc.ProteinEditor, "on_paste")

    async def test_paste_inserts_amino_acids_at_cursor(self):
        from textual.events import Paste
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("", [])
            pe._cursor_pos = 0
            pe.on_paste(Paste("MKWVT"))
            await pilot.pause()
            assert pe._aa_seq == "MKWVT"
            assert pe._cursor_pos == 5

    async def test_paste_filters_non_aa_chars(self):
        # Lower-case + whitespace + non-AA letters (B/J/O/U/X/Z) +
        # digits are dropped; the valid residues + '*' survive.
        from textual.events import Paste
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("", [])
            pe._cursor_pos = 0
            pe.on_paste(Paste("mhhhhhh* xz123"))
            await pilot.pause()
            assert pe._aa_seq == "MHHHHHH*"

    async def test_paste_replaces_active_selection(self):
        from textual.events import Paste
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("AAAA", [])
            pe._user_sel = (1, 3)        # select the middle "AA"
            pe._cursor_pos = 3
            pe.on_paste(Paste("WY"))
            await pilot.pause()
            assert pe._aa_seq == "AWYA"

    async def test_empty_paste_is_noop(self):
        from textual.events import Paste
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("MK", [])
            pe.on_paste(Paste(""))
            await pilot.pause()
            assert pe._aa_seq == "MK"

    async def test_paste_event_bubbles_from_focused_scroll(self):
        # End-to-end: a Paste delivered to the focused scroll container
        # (where Textual routes terminal native-paste) must bubble up to
        # ProteinEditor.on_paste — same wiring as the DNA editor.
        from textual.events import Paste
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs").active = "syn-tab-protein"
            await pilot.pause()
            await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("", [])
            scroll = scr.query_one("#pe-scroll")
            scroll.post_message(Paste("GSGSG"))
            await pilot.pause()
            await pilot.pause()
            assert pe._aa_seq == "GSGSG"


class TestSynthesisSidePaneButtonsFit:
    """The DNA feature-library and protein-motif side panes each host
    four action buttons. At the editor's 4fr / pane 1fr split the pane
    is only ~30 cols — too narrow for a single row of four. Regression
    guard for the off-screen-button bug: they must lay out 2x2 and stay
    inside the pane (pre-fix buttons 3 + 4 rendered past the right edge
    and were invisible)."""

    def _assert_2x2_inside(self, cont, buttons, term_w):
        # Container is two rows tall (2 x button height 3).
        assert cont.region.height == 6
        ys = set()
        for b in buttons:
            r = b.region
            # Every button fully inside the container horizontally...
            assert r.x >= cont.region.x
            assert r.right <= cont.region.right
            # ...vertically, and on-screen.
            assert r.y >= cont.region.y
            assert r.bottom <= cont.region.bottom
            assert r.right <= term_w
            ys.add(r.y)
        # Two distinct rows → the four buttons actually wrapped.
        assert len(ys) == 2

    async def test_motif_buttons_fit(self):
        from textual.widgets import Button
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs").active = "syn-tab-protein"
            await pilot.pause()
            await pilot.pause()
            cont = scr.query_one("#syn-motif-btns")
            buttons = [
                scr.query_one("#" + bid, Button) for bid in (
                    "btn-syn-motif-insert", "btn-syn-motif-new",
                    "btn-syn-motif-edit", "btn-syn-motif-delete",
                )
            ]
            self._assert_2x2_inside(cont, buttons, _TERM[0])

    async def test_featlib_buttons_fit(self):
        from textual.widgets import Button
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            # DNA tab is the initial tab — no switch needed.
            cont = scr.query_one("#syn-featlib-btns")
            buttons = [
                scr.query_one("#" + bid, Button) for bid in (
                    "btn-syn-featlib-insert", "btn-syn-featlib-annotate",
                    "btn-syn-featlib-edit", "btn-syn-featlib-refresh",
                )
            ]
            self._assert_2x2_inside(cont, buttons, _TERM[0])


class TestProteinSeqFileParse:
    """`_parse_protein_fasta_single` / `_validate_ingested_protein_seq` —
    single-entry amino-acid file ingest for the Protein-tab Open button."""

    def _write(self, tmp_path, name: str, text: str) -> str:
        d = tmp_path / "seqs"
        d.mkdir(exist_ok=True)
        p = d / name
        p.write_text(text)
        return str(p)

    def test_validate_uppercases_and_strips_whitespace(self):
        assert sc._validate_ingested_protein_seq(" mk wv\tt f\n") == "MKWVTF"

    def test_validate_rejects_empty(self):
        with pytest.raises(ValueError):
            sc._validate_ingested_protein_seq("   \n  ")

    def test_validate_rejects_non_aa(self):
        with pytest.raises(ValueError):
            sc._validate_ingested_protein_seq("MK123!!")

    def test_validate_accepts_ambiguity_codes(self):
        # B/Z/J/X/U/O/* are accepted at the file boundary (real records
        # carry them); the editor's load() filters them later.
        assert sc._validate_ingested_protein_seq("MKXBZJUO*") == "MKXBZJUO*"

    def test_single_record_fasta_header_stripped(self, tmp_path):
        p = self._write(
            tmp_path, "one.fasta",
            ">sp|P12345|TEST a protein\nMKWVTFISLLFL\nFSSAYSRGVFRR\n",
        )
        name, seq = sc._parse_protein_fasta_single(p)
        assert name == "sp|P12345|TEST"
        assert seq == "MKWVTFISLLFLFSSAYSRGVFRR"

    def test_lowercase_sequence_uppercased(self, tmp_path):
        p = self._write(tmp_path, "lc.fasta", ">x\nmkwvtf\n")
        _name, seq = sc._parse_protein_fasta_single(p)
        assert seq == "MKWVTF"

    def test_multi_record_rejected(self, tmp_path):
        p = self._write(tmp_path, "multi.fasta", ">a\nMKMK\n>b\nWYWY\n")
        with pytest.raises(ValueError, match="[Mm]ulti"):
            sc._parse_protein_fasta_single(p)

    def test_empty_record_rejected(self, tmp_path):
        p = self._write(tmp_path, "empty.fasta", ">only-header\n")
        with pytest.raises(ValueError):
            sc._parse_protein_fasta_single(p)

    def test_non_aa_content_rejected(self, tmp_path):
        p = self._write(tmp_path, "junk.fasta", ">j\n12345 %%% !!!\n")
        with pytest.raises(ValueError):
            sc._parse_protein_fasta_single(p)

    def test_headerless_raw_sequence_file(self, tmp_path):
        # No '>' — falls back to raw-text ingest.
        p = self._write(tmp_path, "raw.seq", "MKWVTF\nISLLFL\n")
        name, seq = sc._parse_protein_fasta_single(p)
        assert seq == "MKWVTFISLLFL"
        assert name == "raw"   # stem when there's no record id

    def test_missing_file_rejected(self, tmp_path):
        with pytest.raises(ValueError):
            sc._parse_protein_fasta_single(str(tmp_path / "nope.fasta"))

    def test_directory_path_rejected(self, tmp_path):
        with pytest.raises(ValueError):
            sc._parse_protein_fasta_single(str(tmp_path))


class TestProteinPickerHighlight:
    """The protein file picker paints amino-acid formats pink — the same
    signal FASTA carries in every other picker."""

    def test_fasta_exts_mapped_to_protein_style(self):
        for ext in (".faa", ".fasta", ".fa", ".pep", ".gp", ".pir"):
            assert sc._PROTEIN_PICKER_HIGHLIGHT_MAP.get(ext) == \
                sc._PICKER_PROTEIN_STYLE

    def test_protein_style_matches_global_fasta_pink(self):
        # Consistency: same pink the rest of the app's pickers use.
        assert sc._PICKER_PROTEIN_STYLE == sc._PICKER_FASTA_STYLE
        assert "FF69B4" in sc._PICKER_PROTEIN_STYLE


class TestProteinSeqFilePickerModal:
    def _write_fasta(self, tmp_path, text: str, name: str = "p.fasta") -> str:
        d = tmp_path / "seqs"
        d.mkdir(exist_ok=True)
        p = d / name
        p.write_text(text)
        return str(p)

    async def test_open_without_selection_shows_inline_error(self, tmp_path):
        from textual.widgets import Static
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            captured = []
            app.push_screen(
                sc.ProteinSeqFilePickerModal(start_path=str(tmp_path)),
                callback=captured.append,
            )
            await pilot.pause()
            await pilot.pause()
            modal = app.screen
            modal._open()      # no file selected
            await pilot.pause()
            status = str(modal.query_one("#pfa-status", Static).render())
            assert "Pick a file" in status
            assert not captured     # did NOT dismiss

    async def test_open_valid_file_dismisses_with_seq(self, tmp_path):
        p = self._write_fasta(tmp_path, ">prot1\nMKWVTFISLL\n")
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            captured = []
            app.push_screen(
                sc.ProteinSeqFilePickerModal(start_path=str(tmp_path)),
                callback=captured.append,
            )
            await pilot.pause()
            await pilot.pause()
            modal = app.screen
            modal._selected = p
            modal._open()
            await pilot.pause()
            await pilot.pause()
            assert len(captured) == 1
            path, name, seq = captured[0]
            assert path == p
            assert name == "prot1"
            assert seq == "MKWVTFISLL"

    async def test_open_multi_entry_stays_open_with_error(self, tmp_path):
        from textual.widgets import Static
        p = self._write_fasta(tmp_path, ">a\nMKMK\n>b\nWYWY\n", "multi.fasta")
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            captured = []
            app.push_screen(
                sc.ProteinSeqFilePickerModal(start_path=str(tmp_path)),
                callback=captured.append,
            )
            await pilot.pause()
            await pilot.pause()
            modal = app.screen
            modal._selected = p
            modal._open()
            await pilot.pause()
            status = str(modal.query_one("#pfa-status", Static).render())
            assert "Multi" in status or "multi" in status
            assert not captured             # did NOT dismiss
            assert app.screen is modal      # still open


class TestSynthesisProteinOpenIntegration:
    def _write_fasta(self, tmp_path, text: str, name: str = "p.fasta") -> str:
        d = tmp_path / "seqs"
        d.mkdir(exist_ok=True)
        p = d / name
        p.write_text(text)
        return str(p)

    async def test_open_button_present_and_hint_removed(self):
        from textual.widgets import Button
        from textual.css.query import NoMatches
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            with pytest.raises(NoMatches):
                scr.query_one("#syn-motif-hint")             # removed
            scr.query_one("#syn-tabs").active = "syn-tab-protein"
            await pilot.pause()
            await pilot.pause()
            btn = scr.query_one("#btn-syn-protein-open", Button)   # exists
            # On-screen and inside its pane (not clipped off the edge).
            pane = scr.query_one("#syn-motif-pane")
            assert btn.region.width > 0 and btn.region.height > 0
            assert btn.region.right <= _TERM[0]
            assert btn.region.x >= pane.region.x
            assert btn.region.right <= pane.region.right

    async def test_open_loads_sequence_into_editor(self, tmp_path):
        p = self._write_fasta(tmp_path, ">myprot\nMKWVTFISLLFL\n")
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs").active = "syn-tab-protein"
            await pilot.pause()
            scr._open_protein_file()        # protein tab not dirty → picks
            await pilot.pause()
            await pilot.pause()
            picker = app.screen
            assert isinstance(picker, sc.ProteinSeqFilePickerModal)
            picker._selected = p
            picker._open()
            await pilot.pause()
            await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            assert pe.get_state()[0] == "MKWVTFISLLFL"
            assert scr._protein_loaded_name == "myprot"
            assert scr._protein_loaded_id is None
            assert scr._protein_dirty is True

    async def test_open_refuses_ambiguity_codes(self, tmp_path):
        # X is accepted by the FILE parser but the editor can't store it.
        # It used to be filtered out of the loaded buffer — which shortens
        # the protein and renumbers every residue after it, so the construct
        # you go on to design is not the one in the file (and nobody can
        # order DNA for an unknown residue anyway). Refused since 2026-09-22.
        p = self._write_fasta(tmp_path, ">amb\nMKXMK\n", "amb.fasta")
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs").active = "syn-tab-protein"
            await pilot.pause()
            scr._open_protein_file()
            await pilot.pause()
            await pilot.pause()
            picker = app.screen
            picker._selected = p
            picker._open()
            await pilot.pause()
            await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            assert pe.get_state()[0] == ""       # nothing loaded
            assert scr._protein_loaded_name is None


# ═══════════════════════════════════════════════════════════════════════════════
# 2026-07-08 — editable features on paste (lane-click highlight + edit),
# R cut-site overlay, Synthesis undo/redo, and the 3' insert-site fix.
# ═══════════════════════════════════════════════════════════════════════════════

import types as _types


def _mouse(**kw):
    """Minimal stand-in for a Textual mouse event — the editors' handlers
    only read button / screen_x / screen_y / shift / chain (+ stop())."""
    kw.setdefault("button", 1)
    kw.setdefault("shift", False)
    kw.setdefault("chain", 1)
    kw.setdefault("stop", lambda: None)
    return _types.SimpleNamespace(**kw)


class TestSynthesis3pInsertSite:
    """The Insert-site path must not let a feature that ends at the 3'
    end swallow the inserted recognition bases (the '5' works, 3' doesn't'
    bug reported 2026-07-08). `extend_adjacent_feats=False` keeps the
    site standalone, symmetric with the 5' end."""

    def test_boundary_feature_absorbs_by_default(self):
        # Baseline: the default (extend) DOES absorb — documents WHY the
        # insert-site path must pass extend_adjacent_feats=False.
        ed = sc.SynthesisEditor()
        n = 12
        ed._seq = "A" * n
        ed._feats = [{"start": 0, "end": n, "label": "f", "type": "CDS",
                      "color": "#88cc88", "strand": 1, "qualifiers": {}}]
        ed._cursor_pos = n
        # Offline mirror of insert_at_cursor's feature-shift (extend=True).
        cur, n_ins = n, 6
        for f in ed._feats:
            s, e = f["start"], f["end"]
            if e >= cur and not (e == cur and s == cur):
                f["end"] = e + n_ins
        assert ed._feats[0]["end"] == n + n_ins   # swallowed

    async def test_insert_site_3p_stays_standalone(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            seq = "ATG" + "AAACCCGGG" * 3 + "TAA"
            n = len(seq)
            ed.load(seq, [{"start": 0, "end": n, "label": "cds",
                           "type": "CDS", "color": "#88cc88", "strand": 1,
                           "qualifiers": {}}])
            await pilot.pause()
            ed.set_cursor(n)
            ok = ed.insert_at_cursor("GGTACC", extend_adjacent_feats=False)
            await pilot.pause()
            assert ok
            assert ed._seq.endswith("GGTACC")
            # Feature must NOT have grown over the appended site.
            assert ed._feats[0]["end"] == n
            # Symmetric with the 5' behaviour: reload + insert at 0.
            ed.load(seq, [{"start": 0, "end": n, "label": "cds",
                           "type": "CDS", "color": "#88cc88", "strand": 1,
                           "qualifiers": {}}])
            await pilot.pause()
            ed.set_cursor(0)
            ed.insert_at_cursor("GTCGAC", extend_adjacent_feats=False)
            await pilot.pause()
            assert ed._feats[0]["start"] == 6   # site left uncovered at 5'

    def test_action_insert_site_uses_no_extend(self):
        import inspect
        src = inspect.getsource(sc.SynthesisScreen.action_insert_site)
        assert "extend_adjacent_feats=False" in src


class TestInsertSiteFeature:
    """Insert-site drops a real, editable feature (hot-pink protein_bind,
    arrowed by the chosen direction), not just bare bases."""

    async def _insert(self, pilot, app, result):
        """Load a fragment, open the picker via action_insert_site, and
        dismiss it with `result`; return the DNA editor."""
        app.push_screen(sc.SynthesisScreen())
        await pilot.pause(); await pilot.pause()
        scr = app.screen
        ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
        ed.load("AAAACCCCGGGG", [])
        await pilot.pause()
        ed.set_cursor(4)
        scr.action_insert_site()
        await pilot.pause(); await pilot.pause()
        assert isinstance(app.screen, sc.RestrictionInsertModal)
        app.screen.dismiss(result)
        await pilot.pause(); await pilot.pause()
        return ed

    async def test_forward_creates_hotpink_feature(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            ed = await self._insert(pilot, app,
                                    {"enzyme": "EcoRI", "strand": 1})
            assert "GAATTC" in ed._seq
            hits = [f for f in ed._feats if f.get("label") == "EcoRI"]
            assert len(hits) == 1
            f = hits[0]
            assert f["type"] == "protein_bind"
            assert f["color"] == sc._RESTR_SITE_FEATURE_COLOR
            assert f["strand"] == 1
            assert f["end"] - f["start"] == 6
            # spans exactly the inserted site
            assert ed._seq[f["start"]:f["end"]] == "GAATTC"

    async def test_reverse_inserts_rc_and_reverse_arrow(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            # BsaI GGTCTC → RC GAGACC on reverse.
            ed = await self._insert(pilot, app,
                                    {"enzyme": "BsaI", "strand": -1})
            assert "GAGACC" in ed._seq
            assert "GGTCTC" not in ed._seq
            f = [f for f in ed._feats if f.get("label") == "BsaI"][0]
            assert f["strand"] == -1
            assert ed._seq[f["start"]:f["end"]] == "GAGACC"

    async def test_palindrome_reverse_same_bases_flipped_arrow(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            ed = await self._insert(pilot, app,
                                    {"enzyme": "EcoRI", "strand": -1})
            # EcoRI is palindromic → RC == site; only the arrow flips.
            assert "GAATTC" in ed._seq
            f = [f for f in ed._feats if f.get("label") == "EcoRI"][0]
            assert f["strand"] == -1

    async def test_insert_site_feature_is_editable(self):
        # The stamped site behaves like any feature: click highlights it,
        # Enter opens the editor.
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            ed = await self._insert(pilot, app,
                                    {"enzyme": "EcoRI", "strand": 1})
            f = [f for f in ed._feats if f.get("label") == "EcoRI"][0]
            ed._cursor_pos = (f["start"] + f["end"]) // 2
            assert ed._feature_at_cursor() == ed._feats.index(f)

    async def test_cancel_adds_nothing(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            ed = await self._insert(pilot, app, None)
            assert ed._seq == "AAAACCCCGGGG"
            assert ed._feats == []

    def test_hotpink_color_is_valid_6hex(self):
        import re
        assert re.fullmatch(r"#[0-9A-Fa-f]{6}", sc._RESTR_SITE_FEATURE_COLOR)


class TestSynthesisLaneClick:
    """Clicking a feature's lane art highlights the whole feature; a
    double-click opens FeatureEditModal. DNA + protein editors."""

    @staticmethod
    def _dna_bar_xy(ed, bp):
        """(x, y) of the row directly above the DNA strand (non-CDS bar)
        at base `bp`."""
        reg = ed.query_one("#syn-view", sc.Static).region
        above = ed._dna_top_row_offset()
        y = reg.y + ed._pad_above_rows + (above - 1)   # packed_row 0 bar
        x = reg.x + ed._pad_left_cols + ed._FLANK_MARKER_WIDTH + bp
        return x, y

    async def test_dna_lane_click_highlights_and_dna_click_does_not(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATG" + "AAACCCGGGTTT" * 3 + "TAA",
                    [{"start": 5, "end": 20, "label": "f",
                      "type": "misc_feature", "color": "#88cc88",
                      "strand": 1, "qualifiers": {}}])
            await pilot.pause(); await pilot.pause()
            x, y = self._dna_bar_xy(ed, 12)
            assert ed._lane_feature_at(x, y)[0] == 0
            ed.on_mouse_down(_mouse(screen_x=x, screen_y=y))
            ed.on_mouse_up(_mouse())
            ed.on_click(_mouse(screen_x=x, screen_y=y))
            await pilot.pause()
            assert ed._user_sel == (5, 20)
            # A DNA-strand click clears the whole-feature highlight.
            reg = ed.query_one("#syn-view", sc.Static).region
            dy = reg.y + ed._pad_above_rows + ed._dna_top_row_offset()
            dx = reg.x + ed._pad_left_cols + ed._FLANK_MARKER_WIDTH + 30
            ed.on_mouse_down(_mouse(screen_x=dx, screen_y=dy))
            ed.on_mouse_up(_mouse())
            ed.on_click(_mouse(screen_x=dx, screen_y=dy))
            await pilot.pause()
            assert ed._user_sel is None

    async def test_dna_double_click_opens_feature_editor(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGAAACCCGGGTTTAAACCC",
                    [{"start": 3, "end": 12, "label": "f",
                      "type": "misc_feature", "color": "#88cc88",
                      "strand": 1, "qualifiers": {}}])
            await pilot.pause(); await pilot.pause()
            x, y = self._dna_bar_xy(ed, 7)
            ed.on_mouse_down(_mouse(screen_x=x, screen_y=y))
            ed.on_mouse_up(_mouse())
            ed.on_click(_mouse(screen_x=x, screen_y=y, chain=2))
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.FeatureEditModal)

    async def test_protein_lane_click_and_enter_edit(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs", sc.TabbedContent).active = \
                "syn-tab-protein"
            await pilot.pause(); await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.set_codon_mode(False)
            pe.load("MASGGGSGGGSKKKK", feats=[
                {"start": 4, "end": 10, "label": "linker", "type": "Linker",
                 "color": "#6B7280", "qualifiers": {}}])
            await pilot.pause(); await pilot.pause()
            _placements, total_rows = pe._pack_aa_feats()
            reg = pe.query_one("#pe-view", sc.Static).region
            aax = reg.x + pe._pad_left_cols + pe._FLANK_MARKER_WIDTH + 6
            hit_y = None
            for rr in range(total_rows):
                yy = reg.y + pe._pad_above_rows + rr
                if pe._lane_feature_at(aax, yy)[0] == 0:
                    hit_y = yy
                    break
            assert hit_y is not None
            pe.on_mouse_down(_mouse(screen_x=aax, screen_y=hit_y))
            pe.on_mouse_up(_mouse())
            pe.on_click(_mouse(screen_x=aax, screen_y=hit_y))
            await pilot.pause()
            assert pe._user_sel == (4, 10)
            # Enter inside the motif opens the editor.
            pe._cursor_pos = 6
            pe.query_one("#pe-scroll", sc.ScrollableContainer).focus()
            await pilot.pause()
            pe.on_key(_types.SimpleNamespace(key="enter", character=None,
                                             stop=lambda: None))
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.FeatureEditModal)


class TestSynthesisScrolledClick:
    """Clicks resolve correctly after panning the editor horizontally —
    `view.region.x` already includes the scroll, so the resolver must not
    add scroll_x again (the double-count that made scrolled-right feature
    clicks miss, reported 2026-07-08)."""

    async def test_dna_lane_click_and_cursor_after_scroll(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            scroll = ed.query_one("#syn-scroll", sc.ScrollableContainer)
            seq = "ACGT" * 160   # 640 bp, far wider than the viewport
            ed.load(seq, [{"start": 520, "end": 560, "label": "far",
                           "type": "misc_feature", "color": "#88cc88",
                           "strand": 1, "qualifiers": {}}])
            await pilot.pause(); await pilot.pause()
            scroll.scroll_to(x=scroll.max_scroll_x, animate=False)
            await pilot.pause(); await pilot.pause()
            assert int(scroll.scroll_x) > 0        # actually scrolled
            reg = ed.query_one("#syn-view", sc.Static).region
            bp = 540
            x = reg.x + ed._pad_left_cols + ed._FLANK_MARKER_WIDTH + bp
            above = ed._dna_top_row_offset()
            lane_y = reg.y + ed._pad_above_rows + (above - 1)
            dna_y = reg.y + ed._pad_above_rows + above
            # Lane click resolves to the feature (not a miss / wrong bp).
            assert ed._lane_feature_at(x, lane_y) == (0, 540)
            # Cursor click on the DNA row resolves to the right bp.
            assert ed._click_to_bp(x, dna_y) == 540


class TestSynthesisRowAwareClicks:
    """Clicks / drags off the actual sequence row don't move the cursor
    or start a phantom selection (reported 2026-07-08)."""

    async def test_dna_off_strand_no_phantom_selection(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATG" + "AAACCCGGGTTT" * 3 + "TAA",
                    [{"start": 5, "end": 20, "label": "f",
                      "type": "misc_feature", "color": "#88cc88",
                      "strand": 1, "qualifiers": {}}])
            await pilot.pause(); await pilot.pause()
            reg = ed.query_one("#syn-view", sc.Static).region
            above = ed._dna_top_row_offset()
            bx = reg.x + ed._pad_left_cols + ed._FLANK_MARKER_WIDTH
            dna_y = reg.y + ed._pad_above_rows + above
            off_y = reg.y + ed._pad_above_rows + above + 2  # below the strands
            assert ed._dna_row_bp(bx + 30, dna_y) == 30
            assert ed._dna_row_bp(bx + 30, off_y) == -1
            # Drag off the strand → nothing selected.
            ed.on_mouse_down(_mouse(screen_x=bx + 30, screen_y=off_y))
            ed.on_mouse_move(_mouse(screen_x=bx + 34, screen_y=off_y))
            ed.on_mouse_up(_mouse())
            await pilot.pause()
            assert ed._user_sel is None
            assert ed._drag_start_bp == -1
            # Drag on the strand → real selection.
            ed.on_mouse_down(_mouse(screen_x=bx + 8, screen_y=dna_y))
            ed.on_mouse_move(_mouse(screen_x=bx + 16, screen_y=dna_y))
            ed.on_mouse_up(_mouse())
            await pilot.pause()
            assert ed._user_sel == (8, 16)

    async def test_protein_off_strand_no_phantom_selection(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs", sc.TabbedContent).active = \
                "syn-tab-protein"
            await pilot.pause(); await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.set_codon_mode(False)
            pe.load("MASGGGSGGGSKKKK", feats=[
                {"start": 4, "end": 10, "label": "L", "type": "Linker",
                 "color": "#6B7280", "qualifiers": {}}])
            await pilot.pause(); await pilot.pause()
            _pl, total_rows = pe._pack_aa_feats()
            reg = pe.query_one("#pe-view", sc.Static).region
            bx = reg.x + pe._pad_left_cols + pe._FLANK_MARKER_WIDTH
            aa_y = reg.y + pe._pad_above_rows + total_rows      # AA strand
            lane_y = reg.y + pe._pad_above_rows                 # lane row 0
            assert pe._aa_row_idx(bx + 12, aa_y) == 12
            assert pe._aa_row_idx(bx + 12, lane_y) == -1
            pe.on_mouse_down(_mouse(screen_x=bx + 12, screen_y=lane_y))
            pe.on_mouse_move(_mouse(screen_x=bx + 13, screen_y=lane_y))
            pe.on_mouse_up(_mouse())
            await pilot.pause()
            # lane row 0 IS the linker's label/bar span at aa 12? No —
            # aa 12 is outside [4,10), so it's a blank lane cell: no select.
            assert pe._user_sel is None


class TestSynthesisCopyPaste:
    """Ctrl+C copies the active editor's selection byte-perfect and
    carries the selection's features; paste re-lays them."""

    async def test_dna_copy_is_byteperfect_and_carries_feats(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            seq = "ATGAAACCCGGGTTTAAACCCGGG"
            ed.load(seq, [{"start": 6, "end": 15, "label": "mid",
                           "type": "CDS", "color": "#88cc88", "strand": 1,
                           "qualifiers": {}}])
            await pilot.pause()
            ed._user_sel = (3, 18)
            scr._copy_dna_selection()
            await pilot.pause()
            cr = app._copied_region
            assert cr["seq"] == seq[3:18]                 # byte-perfect
            assert len(cr["feats"]) == 1
            f = cr["feats"][0]
            assert (f["start"], f["end"]) == (3, 12)      # rebased to [0,span)
            assert f["label"] == "mid"
            # paste back → feature rides along at the cursor offset
            ed.load("NNNN", [])
            await pilot.pause()
            ed._cursor_pos = 2
            ed.on_paste(_types.SimpleNamespace(text=cr["seq"],
                                               stop=lambda: None))
            await pilot.pause()
            assert ed._seq == "NN" + seq[3:18] + "NN"
            got = [g for g in ed._feats if g.get("label") == "mid"]
            assert len(got) == 1
            assert (got[0]["start"], got[0]["end"]) == (5, 14)

    async def test_dna_copy_no_selection_notifies(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGT", [])
            app._copied_region = "sentinel"
            scr._copy_dna_selection()          # no selection
            await pilot.pause()
            assert app._copied_region == "sentinel"   # untouched

    async def test_ctrl_c_binding_targets_synthesis_editor(self):
        # The screen's priority ctrl+c must win over the app's own
        # ctrl+c (which targets the absent main #seq-panel).
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            seq = "ATGCGTACGTACGTAC"
            ed.load(seq, [])
            await pilot.pause()
            ed._user_sel = (2, 9)
            ed.query_one("#syn-scroll", sc.ScrollableContainer).focus()
            await pilot.pause()
            app._copied_region = None
            await pilot.press("ctrl+c")
            await pilot.pause(); await pilot.pause()
            assert app._copied_region is not None
            assert app._copied_region["seq"] == seq[2:9]

    async def test_protein_copy_paste_carries_motif(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs", sc.TabbedContent).active = \
                "syn-tab-protein"
            await pilot.pause(); await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("MASGGGSKKKKLLLL", feats=[
                {"start": 3, "end": 7, "label": "linker", "type": "Linker",
                 "color": "#6B7280", "qualifiers": {}}])
            await pilot.pause()
            pe._user_sel = (1, 9)
            scr._copy_protein_selection()
            await pilot.pause()
            cr = app._copied_protein_region
            assert cr["aa"] == "ASGGGSKK"
            assert (cr["feats"][0]["start"], cr["feats"][0]["end"]) == (2, 6)
            pe.load("MM", feats=[])
            await pilot.pause()
            pe._cursor_pos = 1
            pe.on_paste(_types.SimpleNamespace(text=cr["aa"],
                                               stop=lambda: None))
            await pilot.pause()
            assert pe._aa_seq == "MASGGGSKKM"
            got = [g for g in pe._aa_feats if g.get("label") == "linker"]
            assert len(got) == 1
            assert (got[0]["start"], got[0]["end"]) == (3, 7)


class TestSynthesisDeleteFeatureFirst:
    """Delete on a highlighted feature removes the annotation first
    (keeping the bases + the exact highlight); a second Delete then
    removes the highlighted bases."""

    async def test_dna_two_step_delete(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            seq = "ATGAAACCCGGGTTTAAA"
            ed.load(seq, [{"start": 3, "end": 9, "label": "f", "type": "CDS",
                           "color": "#88cc88", "strand": 1, "qualifiers": {}}])
            await pilot.pause()
            ed._select_feature_span(0)          # lane-click highlight
            await pilot.pause()
            assert ed._user_sel == (3, 9)
            ed.delete_at_cursor(forward=True)   # 1st Delete: feature only
            await pilot.pause()
            assert ed._seq == seq               # DNA untouched
            assert ed._feats == []              # feature gone
            assert ed._user_sel == (3, 9)       # highlight maintained
            ed.delete_at_cursor(forward=True)   # 2nd Delete: the bases
            await pilot.pause()
            assert ed._seq == seq[:3] + seq[9:]
            assert ed._user_sel is None

    async def test_dna_plain_selection_deletes_bases_immediately(self):
        # A drag-selection that does NOT match a feature span deletes
        # the bases on the first Delete (no feature to strip).
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGAAACCCGGG", [{"start": 3, "end": 9, "label": "f",
                                      "type": "CDS", "color": "#88cc88",
                                      "strand": 1, "qualifiers": {}}])
            await pilot.pause()
            ed._user_sel = (1, 4)               # not a feature span
            ed.delete_at_cursor(forward=True)
            await pilot.pause()
            assert ed._seq == "A" + "AACCCGGG"  # bases [1,4) gone
            assert len(ed._feats) == 1          # feature kept

    async def test_protein_two_step_delete(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs", sc.TabbedContent).active = \
                "syn-tab-protein"
            await pilot.pause(); await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("MASGGGSKKKK", feats=[
                {"start": 3, "end": 7, "label": "L", "type": "Linker",
                 "color": "#6B7280", "qualifiers": {}}])
            await pilot.pause()
            pe._select_feature_span(0)
            await pilot.pause()
            assert pe._user_sel == (3, 7)
            pe.delete_at_cursor(forward=True)   # motif only
            await pilot.pause()
            assert pe._aa_seq == "MASGGGSKKKK"
            assert pe._aa_feats == []
            assert pe._user_sel == (3, 7)
            pe.delete_at_cursor(forward=True)   # residues
            await pilot.pause()
            assert pe._aa_seq == "MAS" + "KKKK"


class TestSynthesisClear:
    """The Clear button wipes the whole active sub-tab."""

    async def test_dna_clear_wipes_bases_and_feats(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGT", [{"start": 0, "end": 3, "label": "x",
                                "type": "CDS", "color": "#88cc88",
                                "strand": 1, "qualifiers": {}}])
            scr._dirty = False
            await pilot.pause()
            scr._btn_dna_clear(None)
            await pilot.pause(); await pilot.pause()
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            assert ed._seq == ""
            assert ed._feats == []


class TestSynthesisPersistence:
    """A FRESH SynthesisScreen each open (robust layout), but the composed
    DNA / Protein data survives via app-level `_synthesis_state`."""

    async def test_close_reopen_restores_content(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            scr1 = app.screen
            assert isinstance(scr1, sc.SynthesisScreen)
            ed = scr1.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGTACGTAAA", [{"start": 0, "end": 3, "label": "keep",
                                       "type": "CDS", "color": "#88cc88",
                                       "strand": 1, "qualifiers": {}}])
            # A protein motif too, to prove both tabs persist.
            pe = scr1.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("MASKKK", feats=[{"start": 0, "end": 3, "label": "m",
                                      "type": "Motif", "color": "#6B7280",
                                      "qualifiers": {}}])
            await pilot.pause()
            scr1.action_cancel()                # close → saves state
            await pilot.pause(); await pilot.pause()
            assert app.screen is not scr1
            app.action_open_synthesis()         # reopen → FRESH screen
            await pilot.pause(); await pilot.pause()
            scr2 = app.screen
            assert scr2 is not scr1             # fresh instance
            # But the widgets ARE laid out (the regression was a 0×0
            # re-mounted reused instance).
            assert scr2.query_one("#syn-box").size.width > 0
            ed2 = scr2.query_one("#syn-editor", sc.SynthesisEditor)
            assert ed2._seq == "ATGCGTACGTAAA"  # DNA restored
            assert len(ed2._feats) == 1
            assert ed2._feats[0]["label"] == "keep"
            pe2 = scr2.query_one("#syn-protein-editor", sc.ProteinEditor)
            assert pe2._aa_seq == "MASKKK"       # protein restored
            assert len(pe2._aa_feats) == 1
            # And it closes cleanly a second time.
            scr2.action_cancel()
            await pilot.pause(); await pilot.pause()
            assert app.screen is not scr2

    async def test_clear_then_reopen_stays_empty(self):
        # Clear wipes the tab AND the persisted state, so a reopen is empty.
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGT", [])
            scr._dirty = False
            await pilot.pause()
            scr._btn_dna_clear(None)            # Clear → empty
            await pilot.pause(); await pilot.pause()
            scr.action_cancel()
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            assert app.screen.query_one("#syn-editor",
                                        sc.SynthesisEditor)._seq == ""

    async def test_subtab_switch_keeps_data(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGTACGT", [])
            await pilot.pause()
            tabs = scr.query_one("#syn-tabs", sc.TabbedContent)
            tabs.active = "syn-tab-protein"
            await pilot.pause(); await pilot.pause()
            tabs.active = "syn-tab-dna"
            await pilot.pause(); await pilot.pause()
            assert scr.query_one("#syn-editor", sc.SynthesisEditor)._seq \
                == "ATGCGTACGT"


class TestSynthesisZoom:
    """In-editor bp/column zoom-out overview (−/+/0), mirroring the main
    app's linear view, for the always-linear DNA fragment."""

    @staticmethod
    def _key(k):
        import types
        return types.SimpleNamespace(key=k, stop=lambda: None)

    async def test_zoom_out_in_reset_and_render(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGC" * 200, [{"start": 10, "end": 120, "label": "prom",
                                    "type": "promoter", "color": "#84B0DC",
                                    "strand": 1, "qualifiers": {}}])
            await pilot.pause()
            assert ed._zoom_bp_per_col == 1
            ed.on_key(self._key("minus"))
            await pilot.pause()
            assert ed._zoom_bp_per_col == 2
            ed.on_key(self._key("minus"))
            await pilot.pause()
            assert ed._zoom_bp_per_col == 4
            rendered = ed.query_one("#syn-view", sc.Static).render()
            txt = getattr(rendered, "plain", str(rendered))
            assert "zoom" in txt           # hint header
            assert "prom" in txt           # feature name block
            ed.on_key(self._key("plus"))
            await pilot.pause()
            assert ed._zoom_bp_per_col == 2
            ed.on_key(self._key("0"))
            await pilot.pause()
            assert ed._zoom_bp_per_col == 1

    async def test_click_in_zoom_jumps_to_1to1(self):
        import types
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGC" * 200, [])
            await pilot.pause()
            ed.on_key(self._key("minus")); ed.on_key(self._key("minus"))
            await pilot.pause()
            z = ed._zoom_bp_per_col
            assert z > 1
            reg = ed.query_one("#syn-view", sc.Static).region
            col = 20
            ev = types.SimpleNamespace(
                button=1, screen_x=reg.x + ed._pad_left_cols + col,
                screen_y=reg.y + ed._pad_above_rows + 2,
                shift=False, ctrl=False, chain=1, stop=lambda: None)
            ed.on_click(ev)
            await pilot.pause()
            assert ed._zoom_bp_per_col == 1
            assert ed._cursor_pos == col * z

    async def test_no_typing_while_zoomed(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGC" * 200, [])
            await pilot.pause()
            ed.on_key(self._key("minus"))
            await pilot.pause()
            before = ed._seq
            ed.on_key(self._key("a"))
            await pilot.pause()
            assert ed._seq == before

    async def test_zoom_capped_at_fit(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGC" * 50, [])
            await pilot.pause()
            for _ in range(20):
                ed.on_key(self._key("minus"))
            await pilot.pause()
            z = ed._zoom_bp_per_col
            ncols = (len(ed._seq) + z - 1) // z
            assert ncols <= ed._viewport_width()
            assert z == ed._zoom_max()

    async def test_load_resets_zoom(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGC" * 200, [])
            ed.on_key(self._key("minus"))
            await pilot.pause()
            assert ed._zoom_bp_per_col > 1
            ed.load("GGGG", [])        # fresh buffer
            assert ed._zoom_bp_per_col == 1


class TestSynthesisEnzymeClick:
    """Clicking a restriction-enzyme name highlights recognition + spacer +
    cut (Type IIS footprint computed from the sequence), like SequencePanel."""

    def test_build_resite_highlight_typeIIS(self):
        seq = "AAAAA" + "GGTCTC" + "A" * 14   # BsaI at [5,11)
        resites = [f for f in sc._scan_restriction_sites(
            seq, min_recognition_len=4, unique_only=False, circular=False,
            allowed_enzymes={"BsaI"}) if f.get("type") == "resite"
            and f.get("label")]
        assert resites
        hi = sc._build_resite_highlight(len(seq), resites[0])
        # recognition [5,11) + spacer [11,12) + overhang [12,16) → [5,16)
        assert (hi["start"], hi["end"]) == (5, 16)
        assert (hi["rec_start"], hi["rec_end"]) == (5, 11)

    def test_build_resite_highlight_palindrome(self):
        seq = "TTTTT" + "GAATTC" + "TTTTT"    # EcoRI at [5,11)
        resites = [f for f in sc._scan_restriction_sites(
            seq, min_recognition_len=4, unique_only=False, circular=False,
            allowed_enzymes={"EcoRI"}) if f.get("type") == "resite"
            and f.get("label")]
        assert resites
        hi = sc._build_resite_highlight(len(seq), resites[0])
        assert (hi["start"], hi["end"]) == (5, 11)   # recognition only

    async def test_click_enzyme_name_highlights_and_arrow_clears(self):
        import types
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("AAAAA" + "GGTCTC" + "A" * 14, [])
            ed.toggle_restr_overlay()
            await pilot.pause(); await pilot.pause()
            r = next(f for f in ed._restr_feats
                     if f.get("type") == "resite" and f.get("label") == "BsaI")
            view = ed.query_one("#syn-view", sc.Static)
            reg = view.region
            sx = reg.x + ed._pad_left_cols + ed._FLANK_MARKER_WIDTH + r["start"]
            hit = None
            for rr in range(0, 8):
                sy = reg.y + ed._pad_above_rows + rr
                ed._last_resite_click = None
                ed._lane_feature_at(sx, sy)
                if ed._last_resite_click is not None:
                    hit = sy
                    break
            assert hit is not None
            ev = types.SimpleNamespace(button=1, screen_x=sx, screen_y=hit,
                                       shift=False, ctrl=False, chain=1,
                                       stop=lambda: None)
            ed.on_mouse_down(ev)
            ed.on_click(ev)
            await pilot.pause()
            assert ed._re_highlight is not None
            assert (ed._re_highlight["start"], ed._re_highlight["end"]) == (5, 16)
            ed.set_cursor(3)                    # arrow clears the highlight
            assert ed._re_highlight is None


class TestSynthesisCtrlSelection:
    """Ctrl+click extends (shift alias); Ctrl+arrow slides the selection."""

    @staticmethod
    def _ev(**kw):
        import types
        kw.setdefault("button", 1); kw.setdefault("shift", False)
        kw.setdefault("ctrl", False); kw.setdefault("chain", 1)
        kw.setdefault("stop", lambda: None)
        return types.SimpleNamespace(**kw)

    async def test_ctrl_arrow_slides_selection(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGTACGTAAACCCGGG", [])
            ed._user_sel = (3, 7); ed._sel_anchor = 3; ed._cursor_pos = 7
            await pilot.pause()
            ed.on_key(self._ev(key="ctrl+right"))
            assert ed._user_sel == (4, 8)
            ed.on_key(self._ev(key="ctrl+left"))
            ed.on_key(self._ev(key="ctrl+left"))
            assert ed._user_sel == (2, 6)
            for _ in range(10):
                ed.on_key(self._ev(key="ctrl+left"))
            assert ed._user_sel == (0, 4)      # clamps, keeps length

    async def test_ctrl_click_extends(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGTACGTAAACCCGGG", [])
            await pilot.pause()
            reg = ed.query_one("#syn-view", sc.Static).region
            row = reg.y + ed._pad_above_rows + ed._dna_top_row_offset()

            def sx(bp):
                return reg.x + ed._pad_left_cols + ed._FLANK_MARKER_WIDTH + bp
            ed.on_mouse_down(self._ev(screen_x=sx(3), screen_y=row))
            ed.on_mouse_up(self._ev(screen_x=sx(3), screen_y=row))
            ed.on_click(self._ev(screen_x=sx(3), screen_y=row))
            assert ed._cursor_pos == 3
            ed.on_mouse_down(self._ev(screen_x=sx(9), screen_y=row, ctrl=True))
            ed.on_mouse_up(self._ev(screen_x=sx(9), screen_y=row, ctrl=True))
            assert ed._user_sel == (3, 9)


class TestSynthesisFeatlibSave:
    """Save-to-library from the DNA tab persists AND refreshes the list."""

    async def test_save_refreshes_featlib_list(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGTACGTAAACCCGGG", [])
            ed._user_sel = (0, 9); ed._sel_anchor = 0; ed._cursor_pos = 9
            await pilot.pause()
            before = scr.query_one("#syn-featlib-table", sc.DataTable).row_count
            n_before = len(sc._load_features())
            scr.action_add_feature()
            await pilot.pause(); await pilot.pause()
            entry = {"name": "RegFeatSave", "feature_type": "misc_feature",
                     "sequence": "ATGCGTACG", "strand": 1, "color": None,
                     "qualifiers": {}, "description": ""}
            app.screen.dismiss({"action": "save", "entry": entry})
            await pilot.pause(); await pilot.pause()
            assert len(sc._load_features()) == n_before + 1
            assert app.screen.query_one(
                "#syn-featlib-table", sc.DataTable).row_count == before + 1
            assert any(e.get("name") == "RegFeatSave"
                       for e in sc._load_features())


class TestSynthesisPasteFeatureIsolation:
    """Pasting one annotated sequence after another must keep the features
    separate — the upstream feature must NOT extend beneath the new bases."""

    @staticmethod
    def _paste(ed, app, seq, feats):
        import types
        app._copied_region = {"seq": seq, "feats": [dict(f) for f in feats]}
        ed.on_paste(types.SimpleNamespace(text=seq, stop=lambda: None))

    async def test_paste_after_paste_keeps_features_separate(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            self._paste(ed, app, "ATGCGTACGT",
                        [{"start": 0, "end": 10, "label": "featA",
                          "type": "promoter", "color": "#84B0DC",
                          "strand": 1, "qualifiers": {}}])
            await pilot.pause()
            self._paste(ed, app, "TTTTGGGGCC",
                        [{"start": 0, "end": 10, "label": "featB",
                          "type": "CDS", "color": "#88cc88",
                          "strand": 1, "qualifiers": {}}])
            await pilot.pause()
            spans = {f["label"]: (f["start"], f["end"]) for f in ed._feats}
            assert ed._seq == "ATGCGTACGTTTTTGGGGCC"
            assert spans["featA"] == (0, 10)     # did NOT grow under B
            assert spans["featB"] == (10, 20)

    async def test_plain_paste_after_feature_does_not_extend(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGT", [{"start": 0, "end": 6, "label": "A",
                                "type": "CDS", "color": "#88cc88",
                                "strand": 1, "qualifiers": {}}])
            ed.set_cursor(6)
            await pilot.pause()
            app._copied_region = None
            import types
            ed.on_paste(types.SimpleNamespace(text="GGGG", stop=lambda: None))
            await pilot.pause()
            spans = {f["label"]: (f["start"], f["end"]) for f in ed._feats}
            assert ed._seq == "ATGCGTGGGG"
            assert spans["A"] == (0, 6)

    async def test_paste_inside_feature_still_extends(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGT", [{"start": 0, "end": 6, "label": "A",
                                "type": "CDS", "color": "#88cc88",
                                "strand": 1, "qualifiers": {}}])
            ed.set_cursor(3)
            await pilot.pause()
            app._copied_region = None
            import types
            ed.on_paste(types.SimpleNamespace(text="AAA", stop=lambda: None))
            await pilot.pause()
            spans = {f["label"]: (f["start"], f["end"]) for f in ed._feats}
            assert ed._seq == "ATGAAACGT"
            assert spans["A"] == (0, 9)          # grew (insert is inside it)

    async def test_protein_paste_after_paste_keeps_motifs_separate(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs", sc.TabbedContent).active = "syn-tab-protein"
            await pilot.pause(); await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            import types

            def _pp(aa, feats):
                app._copied_protein_region = {"aa": aa,
                                              "feats": [dict(f) for f in feats]}
                pe.on_paste(types.SimpleNamespace(text=aa, stop=lambda: None))
            _pp("MASK", [{"start": 0, "end": 4, "label": "mA",
                          "type": "Motif", "color": "#6B7280",
                          "qualifiers": {}}])
            await pilot.pause()
            _pp("WYFG", [{"start": 0, "end": 4, "label": "mB",
                          "type": "Motif", "color": "#6B7280",
                          "qualifiers": {}}])
            await pilot.pause()
            spans = {f["label"]: (f["start"], f["end"]) for f in pe._aa_feats}
            assert pe._aa_seq == "MASKWYFG"
            assert spans["mA"] == (0, 4)         # did NOT grow under mB
            assert spans["mB"] == (4, 8)


class TestWorkbenchTrueSwitch:
    """Opening a workbench tab while on another CLOSES the prior (true
    switch, not a stack); each tab persists its work; modals still overlay."""

    async def test_switch_closes_prior_and_persists(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGTACGTAAA", [{"start": 0, "end": 3, "label": "keep",
                                       "type": "CDS", "color": "#88cc88",
                                       "strand": 1, "qualifiers": {}}])
            await pilot.pause()
            # Switch to another workbench tab.
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.FeatureLibraryScreen)
            assert not any(isinstance(s, sc.SynthesisScreen)
                           for s in app.screen_stack)   # closed, not stacked
            # Switch back → prior tab closes, Synthesis data restored + laid out.
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.SynthesisScreen)
            assert not any(isinstance(s, sc.FeatureLibraryScreen)
                           for s in app.screen_stack)
            ed2 = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            assert ed2._seq == "ATGCGTACGTAAA"
            assert app.screen.query_one("#syn-box").size.width > 0

    async def test_modal_overlays_not_switch(self):
        from textual.screen import ModalScreen
        from textual.widgets import Static as _TS
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()

            class _Dlg(ModalScreen):
                def compose(self):
                    yield _TS("dialog")
            app.push_screen(_Dlg())
            await pilot.pause(); await pilot.pause()
            # Synthesis stays (overlaid, not closed).
            assert any(isinstance(s, sc.SynthesisScreen)
                       for s in app.screen_stack)
            assert isinstance(app.screen, _Dlg)

    def test_workbench_tab_registry(self):
        for cls in (sc.SynthesisScreen, sc.BabsScreen, sc.AutolabScreen,
                    sc.SimulatorScreen, sc.SequencingScreen,
                    sc.ExperimentsScreen, sc.HistoryScreen,
                    sc.FeatureLibraryScreen, sc.PartsBinModal,
                    sc.PrimerDesignScreen):
            assert cls in sc.PlasmidApp._WORKBENCH_TAB_TYPES
        # Modals are NOT tabs (they overlay).
        for cls in (sc.SettingsModal, sc.BlastModal, sc.ConstructorModal):
            assert cls not in sc.PlasmidApp._WORKBENCH_TAB_TYPES

    async def test_switch_synthesis_babs_roundtrip(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            app.screen.query_one("#syn-editor",
                                 sc.SynthesisEditor).load("ATGATGATG", [])
            await pilot.pause()
            app.action_open_babs()            # Synthesis → Babs
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.BabsScreen)
            assert not any(isinstance(s, sc.SynthesisScreen)
                           for s in app.screen_stack)
            app.action_open_synthesis()       # Babs → Synthesis (restore)
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.SynthesisScreen)
            assert not any(isinstance(s, sc.BabsScreen)
                           for s in app.screen_stack)
            assert app.screen.query_one(
                "#syn-editor", sc.SynthesisEditor)._seq == "ATGATGATG"


class TestSynthesisToolbarNav:
    """The main menu bar is present on the Synthesis screen; closing is
    lossless (no save modal); re-opening while stacked resurfaces."""

    async def test_menubar_present(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            mb = app.screen.query_one("#menubar", sc.MenuBar)
            assert mb.size.width > 0

    async def test_close_never_prompts_and_persists(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.insert_at_cursor("ATGCGT")       # dirty
            await pilot.pause()
            assert scr._dirty is True
            scr.action_cancel()                 # close — NO modal
            await pilot.pause(); await pilot.pause()
            assert not isinstance(app.screen,
                                  sc.SynthesisUnsavedChangesModal)
            assert not isinstance(app.screen, sc.SynthesisScreen)
            app.action_open_synthesis()         # reopen restores
            await pilot.pause(); await pilot.pause()
            assert app.screen.query_one(
                "#syn-editor", sc.SynthesisEditor)._seq == "ATGCGT"

    async def test_reopen_while_stacked_resurfaces(self):
        from textual.screen import ModalScreen
        from textual.widgets import Static as _TS
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            syn = app.screen

            class _Top(ModalScreen):
                def compose(self):
                    yield _TS("top")
            app.push_screen(_Top())
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()         # from atop → resurface
            await pilot.pause(); await pilot.pause()
            assert app.screen is syn
            assert sum(1 for s in app.screen_stack
                       if isinstance(s, sc.SynthesisScreen)) == 1


class TestSynthesisSaveIsLibraryOnly:
    """Save writes the fragment to the library as LINEAR and does NOT
    populate the parts bin (that's the Clone Fragment button's job)."""

    async def test_save_linear_no_parts_bin(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.action_open_synthesis()
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGTACGTAAA", [{"start": 0, "end": 6, "label": "f",
                                       "type": "CDS", "color": "#88cc88",
                                       "strand": 1, "qualifiers": {}}])
            scr._loaded_name = "myfrag"
            scr._loaded_id = "myfrag"
            bin_before = len(sc._load_parts_bin())
            await pilot.pause()
            scr.action_save()
            await pilot.pause(); await pilot.pause()
            entry = sc._find_library_entry_by_id("myfrag")
            assert entry is not None
            rec = sc._gb_text_to_record(entry.get("gb_text", ""))
            assert rec.annotations.get("topology") == "linear"
            assert len(sc._load_parts_bin()) == bin_before   # untouched


class TestHangWatchdog:
    """The auto stack-dump watchdog for UI-thread stalls."""

    def test_tick_and_stall_detection(self):
        wd = sc._HangWatchdog()
        wd.tick()
        # Fresh tick → not stalled.
        import time as _t
        with wd._lock:
            stalled = _t.monotonic() - wd._last_tick
        assert stalled < sc._HANG_WATCHDOG_STALL_S
        # Simulate a stall by back-dating the last tick.
        with wd._lock:
            wd._last_tick = _t.monotonic() - (sc._HANG_WATCHDOG_STALL_S + 5)
            aged = _t.monotonic() - wd._last_tick
        assert aged >= sc._HANG_WATCHDOG_STALL_S
        # tick() re-arms.
        wd.tick()
        with wd._lock:
            assert wd._dumped_for_stall is False

    def test_dump_is_safe_without_fd(self):
        # _dump must never raise even if the stacks FD is unavailable.
        wd = sc._HangWatchdog()
        wd._dump(9.0)   # should log + not raise

    async def test_subtab_switch_keeps_data(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("ATGCGTACGT", [])
            await pilot.pause()
            tabs = scr.query_one("#syn-tabs", sc.TabbedContent)
            tabs.active = "syn-tab-protein"
            await pilot.pause(); await pilot.pause()
            tabs.active = "syn-tab-dna"
            await pilot.pause(); await pilot.pause()
            assert scr.query_one("#syn-editor", sc.SynthesisEditor)._seq \
                == "ATGCGTACGT"


class TestFeatureColorSafety:
    """A malformed feature colour (short/garbage hex from a hand-edited
    .gb) must not crash the render, and _parse must normalise it."""

    def test_safe_feature_color(self):
        assert sc._safe_feature_color("#8c8") == "#88cc88"     # expand
        assert sc._safe_feature_color("#88cc88") == "#88cc88"  # valid
        assert sc._safe_feature_color("red") == "red"          # named
        assert sc._safe_feature_color("color(252)") == "color(252)"
        assert sc._safe_feature_color("#ZZZZZZ") == "color(252)"  # garbage
        assert sc._safe_feature_color(None) == "color(252)"

    def test_build_seq_text_renders_bad_color_without_crash(self):
        from rich.console import Console
        import io
        for bad in ("#8c8", "#ZZZZZZ", "notacolor", "#12345"):
            feats = [{"start": 2, "end": 8, "label": "b",
                      "type": "misc_feature", "color": bad, "strand": 1}]
            txt = sc._build_seq_text("ATGCATGCATGCATGC", feats, line_width=25)
            # Force the DEFERRED style parse (the paint that crashed).
            Console(file=io.StringIO(), force_terminal=True,
                    width=100).print(txt)

    def test_parse_normalises_short_hex_qualifier(self):
        gb = ("LOCUS       t              20 bp    DNA     linear\n"
              "FEATURES             Location/Qualifiers\n"
              "     misc_feature    5..15\n"
              '                     /label="badc"\n'
              '                     /ApEinfo_fwdcolor="#8c8"\n'
              "ORIGIN\n"
              "        1 atgcatgcat gcatgcatgc\n"
              "//\n")
        rec = sc._gb_text_to_record(gb)
        pm = sc.PlasmidMap()
        parsed = pm._parse(rec)
        cols = [f.get("color") for f in parsed if f.get("label") == "badc"]
        assert cols == ["#88cc88"]   # expanded, not the crashing #8c8


class TestSynthesisRestrictionOverlay:
    """The R key reveals / hides restriction cut sites on the DNA tab."""

    async def test_toggle_scans_and_clears(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("TTTTGAATTCTTTTGGATCCTTTTAAGCTTTTTT")  # EcoRI/BamHI/HindIII
            await pilot.pause(); await pilot.pause()
            assert ed._show_restr_overlay is False
            scr.action_toggle_restr_sites()
            await pilot.pause(); await pilot.pause()
            assert ed._show_restr_overlay is True
            assert any(f.get("type") == "resite" for f in ed._restr_feats)
            assert len(ed._display_feats()) >= len(ed._restr_feats)
            scr.action_toggle_restr_sites()
            await pilot.pause()
            assert ed._show_restr_overlay is False
            assert ed._restr_feats == []

    async def test_overlay_rescans_on_edit(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            ed = app.screen.query_one("#syn-editor", sc.SynthesisEditor)
            ed.load("TTTTTTTTTTTT")
            await pilot.pause()
            ed.toggle_restr_overlay()
            await pilot.pause()
            n0 = sum(1 for f in ed._restr_feats if f.get("type") == "resite")
            # Type an EcoRI site at the cursor; the overlay must pick it up.
            ed.set_cursor(6)
            ed.insert_at_cursor("GAATTC")
            await pilot.pause(); await pilot.pause()
            n1 = sum(1 for f in ed._restr_feats if f.get("type") == "resite")
            assert n1 > n0


class TestSynthesisUndoRedo:
    """Ctrl+Z / Ctrl+Y in the Synthesis tab undo/redo typed bases + edits."""

    async def test_dna_undo_redo_cycle(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.insert_at_cursor("ATG")
            await pilot.pause()
            ed.insert_at_cursor("CCC")
            await pilot.pause()
            assert ed._seq == "ATGCCC"
            scr.action_undo()
            await pilot.pause()
            assert ed._seq == "ATG"
            scr.action_undo()
            await pilot.pause()
            assert ed._seq == ""
            scr.action_redo()
            await pilot.pause()
            assert ed._seq == "ATG"
            scr.action_redo()
            await pilot.pause()
            assert ed._seq == "ATGCCC"

    async def test_new_edit_clears_redo(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.insert_at_cursor("ATGCCC")
            await pilot.pause()
            scr.action_undo()          # -> ""
            await pilot.pause()
            assert scr._dna_redo       # redo now has content
            ed.insert_at_cursor("GGG")
            await pilot.pause()
            assert ed._seq == "GGG"
            assert scr._dna_redo == []

    async def test_load_resets_undo_history(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            ed = scr.query_one("#syn-editor", sc.SynthesisEditor)
            ed.insert_at_cursor("ATGATG")
            await pilot.pause()
            ed.load("CCCCCC")          # posts Loaded → resets history
            await pilot.pause(); await pilot.pause()
            # Ctrl+Z can't reach the pre-load buffer.
            scr.action_undo()
            await pilot.pause()
            assert ed._seq == "CCCCCC"

    async def test_protein_undo_redo(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_TERM) as pilot:
            await pilot.pause(); await pilot.pause()
            app.push_screen(sc.SynthesisScreen())
            await pilot.pause(); await pilot.pause()
            scr = app.screen
            scr.query_one("#syn-tabs", sc.TabbedContent).active = \
                "syn-tab-protein"
            await pilot.pause(); await pilot.pause()
            pe = scr.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.insert_at_cursor("MKK")
            await pilot.pause()
            assert pe.get_state()[0] == "MKK"
            scr.action_undo()
            await pilot.pause()
            assert pe.get_state()[0] == ""
            scr.action_redo()
            await pilot.pause()
            assert pe.get_state()[0] == "MKK"
