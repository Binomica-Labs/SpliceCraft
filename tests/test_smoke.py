"""
test_smoke — Textual TUI bootstrap smoke tests.

These are the minimum checks that a human running `python3 splicecraft.py` with
a pre-loaded GenBank file will not hit a Python error during mount, compose, or
the first render pass. They are NOT pixel-level rendering tests.

All tests run with `asyncio_mode = "auto"` (see pyproject.toml) so async test
functions are picked up without a `@pytest.mark.asyncio` decorator.

Each test starts the app with a synthetic SeqRecord via `_preload_record` so
NO network (NCBI) access is required, and isolates the library JSON with the
`isolated_library` fixture so the real library file is never touched.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


TERMINAL_SIZE = (160, 48)   # wide enough for the three-pane layout


def _build_app(tiny_record, isolated_library) -> sc.PlasmidApp:
    """Build a PlasmidApp with a pre-loaded record. `isolated_library` is
    required as a parameter even though we don't touch it here — it's a
    fixture side-effect that monkeypatches `_LIBRARY_FILE`."""
    app = sc.PlasmidApp()
    app._preload_record = tiny_record
    return app


# ═══════════════════════════════════════════════════════════════════════════════
# Bootstrap
# ═══════════════════════════════════════════════════════════════════════════════

class TestAppBootstrap:
    async def test_app_mounts_with_preloaded_record(self, tiny_record,
                                                     isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            # Let the on_mount _apply_record call_after_refresh run
            await pilot.pause()
            await pilot.pause(0.05)
            assert app._current_record is not None
            assert app._current_record.id == tiny_record.id

    async def test_all_panels_present(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Every one of these must exist; query_one raises if not.
            app.query_one("#plasmid-map", sc.PlasmidMap)
            app.query_one("#sidebar", sc.FeatureSidebar)
            app.query_one("#seq-panel", sc.SequencePanel)
            app.query_one("#library", sc.LibraryPanel)

    async def test_features_loaded_into_map(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            # The fixture has at least 2 features (CDS + misc_feature); the
            # load path may add a 'source' record. Assert non-empty.
            assert len(pm._feats) >= 2

    async def test_sequence_panel_has_sequence(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            assert sp._seq == str(tiny_record.seq)

    async def test_restriction_scan_ran_on_load(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # _restr_cache should be populated (tiny_record contains EcoRI
            # sites; depending on unique_only filter it may yield 0 or more
            # hits — here we just check the field was set to a list).
            assert isinstance(app._restr_cache, list)

    async def test_empty_app_mounts_without_preload(self, isolated_library):
        """App must also mount cleanly with no preloaded record. Pre-populate
        the library with a dummy entry (using the correct `size` field schema
        — see LibraryPanel._repopulate line ~2010) so the on_mount seeder's
        `not _load_library()` guard is False and no network fetch is attempted.
        """
        app = sc.PlasmidApp()
        sc._save_library([{
            "name":    "dummy",
            "id":      "DUMMY",
            "size":    1,
            "n_feats": 0,
            "source":  "test",
            "added":   "2026-04-11",
            "gb_text": "",
        }])
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            assert app._current_record is None


# ═══════════════════════════════════════════════════════════════════════════════
# Basic interactions
# ═══════════════════════════════════════════════════════════════════════════════

class TestBasicKeybindings:
    async def test_rotation_keys_change_origin(self, tiny_record,
                                                isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            origin_before = pm.origin_bp
            await pilot.press("[")
            await pilot.pause(0.1)
            assert pm.origin_bp != origin_before

    async def test_view_toggle_key(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            view_before = getattr(pm, "_view_mode", None) or \
                          getattr(pm, "view_mode", None)
            await pilot.press("v")
            await pilot.pause(0.1)
            view_after = getattr(pm, "_view_mode", None) or \
                         getattr(pm, "view_mode", None)
            # If the widget uses a private attr we may not find it — soft check
            if view_before is not None:
                assert view_before != view_after

    async def test_restr_toggle_changes_state(self, tiny_record,
                                               isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            before = app._show_restr
            await pilot.press("r")
            await pilot.pause(0.1)
            assert app._show_restr != before


# ═══════════════════════════════════════════════════════════════════════════════
# No network / no library pollution guard
# ═══════════════════════════════════════════════════════════════════════════════

class TestNoNetworkAccess:
    async def test_mount_does_not_call_fetch_genbank(self, tiny_record,
                                                      isolated_library,
                                                      monkeypatch):
        """If _preload_record is set, the app must never fall through to
        _seed_default_library, which would call fetch_genbank and try NCBI."""
        calls = []

        def _fake_fetch(*args, **kwargs):
            calls.append((args, kwargs))
            raise RuntimeError("fetch_genbank should not be called in tests")

        monkeypatch.setattr(sc, "fetch_genbank", _fake_fetch)
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            assert not calls, f"fetch_genbank was called {len(calls)} time(s)"


# ═══════════════════════════════════════════════════════════════════════════════
# pLannotate UI entry points (button + shortcut, pLannotate itself mocked)
# ═══════════════════════════════════════════════════════════════════════════════

class TestDeleteFocusRouting:
    """Delete key must be focus-aware: pressing Delete with library focus
    should offer to delete the library entry (with a confirmation defaulting
    to No), NOT silently delete a feature the user forgot they had selected
    in the map. Pressing Delete elsewhere still deletes the selected feature."""

    async def test_focus_is_in_library_helper_true_for_library(
        self, tiny_record, isolated_library
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Focus the library's DataTable
            lib_table = app.query_one("#lib-table")
            lib_table.focus()
            await pilot.pause(0.05)
            assert app._focus_is_in_library() is True

    async def test_focus_is_in_library_helper_false_for_map(
        self, tiny_record, isolated_library
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm.focus()
            await pilot.pause(0.05)
            assert app._focus_is_in_library() is False

    async def test_library_focus_clears_feature_selection(
        self, tiny_record, isolated_library
    ):
        """When focus moves INTO the library from elsewhere, any currently-
        selected feature in the map must be deselected. Mount auto-focuses
        the library table on first load, so we explicitly move focus to the
        map first to create a real transition."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            if len(pm._feats) == 0:
                pytest.skip("fixture has no features")
            # Move focus OUT of the library first (mount auto-focused it)
            pm.focus()
            await pilot.pause()
            await pilot.pause(0.05)
            pm.select_feature(0)
            assert pm.selected_idx == 0
            # Now focus the library's DataTable — this is the real transition.
            # GainedFocus dispatch is async; pause twice to let the message
            # be posted, routed, and the handler run.
            app.query_one("#lib-table").focus()
            await pilot.pause()
            await pilot.pause(0.1)
            assert pm.selected_idx == -1, (
                "feature selection should clear when library gains focus"
            )

    async def test_delete_with_library_focus_opens_confirm_modal(
        self, tiny_record, isolated_library
    ):
        """Delete key with library focused must push the confirmation modal,
        NOT silently delete a feature."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Focus the library DataTable and make sure cursor is on a row
            lib_table = app.query_one("#lib-table")
            lib_table.focus()
            await pilot.pause(0.05)
            # Invoke the action directly to avoid focus/key-routing races
            app.action_delete_feature()
            await pilot.pause(0.05)
            # The modal should now be on top of the screen stack
            from splicecraft import LibraryDeleteConfirmModal
            top = app.screen
            assert isinstance(top, LibraryDeleteConfirmModal), (
                f"expected LibraryDeleteConfirmModal on top, got {type(top).__name__}"
            )

    async def test_confirm_modal_default_focus_is_no(
        self, tiny_record, isolated_library
    ):
        """Modal mounts → the [No] button must be focused. This is the whole
        point of the dialog — Enter should be a safe no-op."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.query_one("#lib-table").focus()
            await pilot.pause(0.05)
            app.action_delete_feature()
            await pilot.pause(0.05)
            from splicecraft import LibraryDeleteConfirmModal
            modal = app.screen
            assert isinstance(modal, LibraryDeleteConfirmModal)
            no_btn = modal.query_one("#btn-libdel-no", sc.Button)
            # Either app.focused IS the No button, or the No button has
            # `has_focus` set
            assert app.focused is no_btn or no_btn.has_focus, (
                f"expected [No] focused; got {app.focused!r}"
            )

    async def test_confirm_no_keeps_entry_in_library(
        self, tiny_record, isolated_library
    ):
        """Pressing No in the dialog must leave the library unchanged."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            before_ids = [e["id"] for e in sc._load_library()]
            app.query_one("#lib-table").focus()
            await pilot.pause(0.05)
            app.action_delete_feature()
            await pilot.pause(0.05)
            from splicecraft import LibraryDeleteConfirmModal
            modal = app.screen
            assert isinstance(modal, LibraryDeleteConfirmModal)
            modal.dismiss(False)
            await pilot.pause(0.05)
            after_ids = [e["id"] for e in sc._load_library()]
            assert after_ids == before_ids

    async def test_confirm_yes_removes_entry_from_library(
        self, tiny_record, isolated_library
    ):
        """Pressing Yes in the dialog must delete the highlighted entry."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # tiny_record was auto-persisted during mount; it should be in lib
            before_ids = [e["id"] for e in sc._load_library()]
            assert tiny_record.id in before_ids
            app.query_one("#lib-table").focus()
            await pilot.pause(0.05)
            # Move DataTable cursor to the tiny_record row (should already be
            # there since it's the only entry)
            app.action_delete_feature()
            await pilot.pause(0.05)
            from splicecraft import LibraryDeleteConfirmModal
            modal = app.screen
            assert isinstance(modal, LibraryDeleteConfirmModal)
            modal.dismiss(True)
            await pilot.pause(0.05)
            after_ids = [e["id"] for e in sc._load_library()]
            assert tiny_record.id not in after_ids, (
                f"expected {tiny_record.id} removed; library now: {after_ids}"
            )

    async def test_delete_with_map_focus_still_deletes_feature(
        self, tiny_record, isolated_library
    ):
        """Classic feature-delete path must still work when the library does
        NOT have focus. Guards against over-broad routing."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            if len(pm._feats) == 0:
                pytest.skip("fixture has no features")
            n_feats_before = len(pm._feats)
            # Focus the map and select a feature
            pm.focus()
            pm.select_feature(0)
            await pilot.pause(0.05)
            assert not app._focus_is_in_library()
            app.action_delete_feature()
            await pilot.pause(0.05)
            # Feature should be gone
            pm_after = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert len(pm_after._feats) == n_feats_before - 1


class TestImportAutoPersist:
    """Every 'user imports a plasmid' entry point should auto-save the
    record to the library. Library loads, pLannotate merges, and undo/redo
    should NOT re-save."""

    async def test_preload_record_is_auto_added_to_library(
        self, tiny_record, isolated_library
    ):
        """A CLI-preloaded record (python3 splicecraft.py myplasmid.gb)
        should appear in the library JSON after mount."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Library JSON on disk should contain the record
            lib_entries = sc._load_library()
            ids = [e["id"] for e in lib_entries]
            assert tiny_record.id in ids, (
                f"preloaded record {tiny_record.id} not saved to library; "
                f"library contains {ids}"
            )

    async def test_library_load_does_not_duplicate(
        self, tiny_record, isolated_library
    ):
        """Clicking a library row fires _library_load → _apply_record (NOT
        _import_and_persist), so the same record must not be added twice."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            before = len(sc._load_library())
            # Simulate clicking the library row by sending the same message
            # the DataTable would post
            app.post_message(
                sc.LibraryPanel.PlasmidLoad(sc._load_library()[0])
            )
            await pilot.pause(0.05)
            after = len(sc._load_library())
            assert after == before, (
                f"library_load should not add entries: {before} → {after}"
            )

    async def test_fetch_callback_adds_to_library(
        self, tiny_record, isolated_library, monkeypatch
    ):
        """When FetchModal dismisses with a record, the app callback
        (_import_and_persist) should save it to the library."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Call the helper directly — the modal → callback route is
            # equivalent to this once the modal dismisses.
            app._import_and_persist(tiny_record)
            await pilot.pause(0.05)
            entries = sc._load_library()
            assert any(e["id"] == tiny_record.id for e in entries), (
                f"fetched record not persisted; library: "
                f"{[e['id'] for e in entries]}"
            )

    async def test_import_of_duplicate_id_updates_in_place(
        self, tiny_record, isolated_library, monkeypatch
    ):
        """Re-importing a record with the same id should update the existing
        entry rather than create a duplicate (the add_entry dedup contract)."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._import_and_persist(tiny_record)
            await pilot.pause(0.05)
            n_first  = len(sc._load_library())
            app._import_and_persist(tiny_record)
            await pilot.pause(0.05)
            n_second = len(sc._load_library())
            assert n_first == n_second, (
                f"re-import duplicated the entry: {n_first} → {n_second}"
            )
            # And the record is present exactly once
            ids = [e["id"] for e in sc._load_library()]
            assert ids.count(tiny_record.id) == 1

    async def test_import_none_is_noop(self, isolated_library):
        """Cancelled fetch/open modals dismiss with None — the helper must
        handle it silently."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            before = len(sc._load_library())
            app._import_and_persist(None)
            await pilot.pause(0.05)
            assert len(sc._load_library()) == before


class TestPlannotateUIEntryPoints:
    async def test_annotate_button_exists_in_library(self, tiny_record,
                                                      isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            btn = app.query_one("#btn-lib-annot", sc.Button)
            assert btn is not None

    async def test_shift_a_binding_registered(self):
        keys = [b.key for b in sc.PlasmidApp.BINDINGS]
        assert "A" in keys, "shift+A (key='A') binding is missing"
        # And it's distinct from lowercase a
        assert "a" in keys

    async def test_shift_a_with_no_record_notifies_not_crashes(
        self, isolated_library, monkeypatch
    ):
        """With no record loaded, Shift+A must notify a warning and return
        without touching pLannotate."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # No record loaded; invoke the action directly (avoids key
            # routing which may target a different widget)
            app.action_annotate_plasmid()
            await pilot.pause(0.05)
            # The app should still be alive — assertion is "didn't crash"

    async def test_annotate_action_with_plannotate_missing_notifies_install(
        self, tiny_record, isolated_library, monkeypatch
    ):
        """With pLannotate absent from PATH, the action notifies instead of
        attempting to run anything. Verify by counting subprocess.run calls."""
        import shutil, subprocess
        monkeypatch.setattr(sc, "_PLANNOTATE_CHECK_CACHE", None)
        monkeypatch.setattr(shutil, "which", lambda *a, **k: None)
        calls = []
        monkeypatch.setattr(
            subprocess, "run",
            lambda *a, **k: calls.append((a, k)) or None,
        )
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.action_annotate_plasmid()
            await pilot.pause(0.05)
            assert not calls, (
                "subprocess.run should NOT be invoked when pLannotate "
                "is not on PATH"
            )
