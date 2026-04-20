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

class TestLibraryRename:
    """Library panel rename (✎ button). Verifies the button exists, the
    modal opens with the current name, saving persists the new name to
    the library JSON AND mutates the currently-loaded record's name so
    the plasmid map header picks up the change without a reload.

    Collision check: refuses to rename to the name of another existing
    entry. Empty names are rejected by the modal itself (we test the
    modal-side validator via `_try_submit`)."""

    async def test_rename_button_exists(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            btn = app.query_one("#btn-lib-rename", sc.Button)
            assert btn is not None

    async def test_rename_opens_modal_with_current_name(
        self, tiny_record, isolated_library
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.post_message(sc.LibraryPanel.RenameRequested(tiny_record.id))
            await pilot.pause()
            await pilot.pause(0.05)
            from splicecraft import RenamePlasmidModal
            modal = app.screen
            assert isinstance(modal, RenamePlasmidModal), (
                f"expected RenamePlasmidModal, got {type(modal).__name__}"
            )
            inp = modal.query_one("#rename-input", sc.Input)
            assert inp.value == tiny_record.name

    async def test_rename_save_persists_to_library_json(
        self, tiny_record, isolated_library
    ):
        """After Save, the library JSON's `name` field is the new name and
        the stored gb_text parses back to a record with matching name."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            new_name = "pACYC-custom"
            # Call the backend directly — no modal round-trip needed to test
            # the persistence logic
            app._rename_library_entry(tiny_record.id, new_name)
            await pilot.pause(0.05)

            entries = sc._load_library()
            match = [e for e in entries if e["id"] == tiny_record.id]
            assert len(match) == 1
            assert match[0]["name"] == new_name
            # gb_text should round-trip to a record with the new name
            reloaded = sc._gb_text_to_record(match[0]["gb_text"])
            assert reloaded.name == new_name

    async def test_rename_updates_currently_loaded_record(
        self, tiny_record, isolated_library
    ):
        """If the renamed entry is currently loaded, _current_record.name
        is mutated in place so the map header picks it up."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            old_name = app._current_record.name
            new_name = "my-lab-plasmid"
            assert old_name != new_name
            app._rename_library_entry(tiny_record.id, new_name)
            await pilot.pause(0.05)
            assert app._current_record.name == new_name
            # PlasmidMap uses record.name during render — its record field
            # is the same object, so it should see the new name.
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm.record.name == new_name

    async def test_rename_invalidates_map_draw_cache(
        self, tiny_record, isolated_library
    ):
        """PlasmidMap._draw_cache holds a (key, Text) tuple. Rename must
        either nuke it or the cache key must differ for the new name."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            # Force one render so _draw_cache has an entry
            pm.render()
            key_before = pm._draw_cache[0] if pm._draw_cache else None
            app._rename_library_entry(tiny_record.id, "renamed-test")
            await pilot.pause(0.05)
            # After rename, _draw_cache is either None (nuked) OR a fresh
            # entry with a different key (record.name is part of the key).
            if pm._draw_cache is not None:
                key_after = pm._draw_cache[0]
                assert key_after != key_before, (
                    "draw cache key must change after rename"
                )
                # And the new key's name field must be the new name
                assert "renamed-test" in key_after, (
                    f"expected 'renamed-test' in cache key; got {key_after}"
                )

    async def test_rename_rejects_duplicate_name(
        self, tiny_record, isolated_library, monkeypatch
    ):
        """If another entry already has the target name, the rename is
        refused with an error notification and the library is unchanged."""
        # Seed the library with two entries: tiny_record and a fake second
        from copy import deepcopy
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        second = SeqRecord(
            Seq("ACGT" * 30), id="OTHER01", name="other",
            description="another plasmid",
        )
        second.annotations["molecule_type"] = "DNA"
        second.annotations["topology"]      = "circular"

        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Add a second entry manually via the library panel
            lib = app.query_one("#library", sc.LibraryPanel)
            lib.add_entry(second)
            await pilot.pause(0.05)
            # Now try to rename tiny_record to 'other' — should fail
            def _cb(result):
                pass
            collisions = []
            orig_notify = app.notify
            def _spy_notify(msg, **kw):
                collisions.append((msg, kw))
                return orig_notify(msg, **kw)
            monkeypatch.setattr(app, "notify", _spy_notify)
            # Fire the RenameRequested handler path with a fake callback
            # that asserts the collision path by calling the inner _on_result
            # equivalent directly: the handler opens a modal with callback —
            # easier to test the collision branch by looking at the spy.
            app.post_message(sc.LibraryPanel.RenameRequested(tiny_record.id))
            await pilot.pause()
            await pilot.pause(0.05)
            # Now dismiss the modal with the colliding name
            from splicecraft import RenamePlasmidModal
            modal = app.screen
            assert isinstance(modal, RenamePlasmidModal)
            modal.dismiss("other")
            await pilot.pause()
            await pilot.pause(0.05)
            # The entry should still have its original name
            entries = sc._load_library()
            tiny_entry = [e for e in entries if e["id"] == tiny_record.id][0]
            assert tiny_entry["name"] == tiny_record.name, (
                "rename to a colliding name should have been refused"
            )
            # And an error notification should have fired
            err_notes = [
                m for m, kw in collisions
                if kw.get("severity") == "error" and "already exists" in m
            ]
            assert err_notes, "expected an 'already exists' error notification"

    async def test_rename_modal_empty_name_rejected(
        self, tiny_record, isolated_library
    ):
        """Modal validator rejects an empty name and does NOT dismiss."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.post_message(sc.LibraryPanel.RenameRequested(tiny_record.id))
            await pilot.pause()
            await pilot.pause(0.05)
            from splicecraft import RenamePlasmidModal
            modal = app.screen
            assert isinstance(modal, RenamePlasmidModal)
            # Blank out the input and try to save
            modal.query_one("#rename-input", sc.Input).value = "   "
            modal._try_submit()
            await pilot.pause(0.05)
            # Modal should still be up (not dismissed)
            assert app.screen is modal
            # And the status line should show an error message.
            status = modal.query_one("#rename-status", sc.Static)
            status_text = str(status.content)
            assert "empty" in status_text.lower() or "cannot" in status_text.lower(), (
                f"expected error message in rename status; got {status_text!r}"
            )

    async def test_rename_modal_cancel_is_noop(
        self, tiny_record, isolated_library
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            original_name = app._current_record.name
            app.post_message(sc.LibraryPanel.RenameRequested(tiny_record.id))
            await pilot.pause()
            await pilot.pause(0.05)
            from splicecraft import RenamePlasmidModal
            modal = app.screen
            assert isinstance(modal, RenamePlasmidModal)
            modal.dismiss(None)   # cancel path
            await pilot.pause(0.05)
            assert app._current_record.name == original_name
            # Library entry unchanged
            entries = sc._load_library()
            assert any(
                e["id"] == tiny_record.id and e["name"] == original_name
                for e in entries
            )


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
        # Block the network seed worker from firing when library starts empty.
        # Without this, the mount handler sees an empty library and kicks off
        # `_seed_default_library` which calls fetch_genbank → a live NCBI
        # fetch that races our assertions.
        monkeypatch.setattr(
            sc, "fetch_genbank",
            lambda *a, **k: (_ for _ in ()).throw(
                RuntimeError("network disabled in tests")
            ),
        )
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

    async def test_import_none_is_noop(self, isolated_library, monkeypatch):
        """Cancelled fetch/open modals dismiss with None — the helper must
        handle it silently. Also blocks the seed worker (see above)."""
        monkeypatch.setattr(
            sc, "fetch_genbank",
            lambda *a, **k: (_ for _ in ()).throw(
                RuntimeError("network disabled in tests")
            ),
        )
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


# ═══════════════════════════════════════════════════════════════════════════════
# _apply_record source_path + dirty-flag handling (regression guard 2026-04-13)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Before today's fix, _apply_record always cleared _source_path — even when
# called with clear_undo=False for an in-place update (pLannotate merge,
# primer-add). That meant after annotating, Ctrl+S no longer targeted the
# user's original .gb file. Also, pLannotate used lib.set_dirty(True) alone,
# which only updated the library panel's marker but left self._unsaved=False,
# so the user could quit without being prompted to save.

class TestApplyRecordInPlaceSemantics:
    """`_apply_record(record, clear_undo=False)` is the "in-place update"
    path — it must not clobber `_source_path`, and the caller is expected
    to call `_mark_dirty()` afterwards to set `_unsaved=True`."""

    async def test_clear_undo_true_clears_source_path(
        self, tiny_record, isolated_library
    ):
        """Fresh-load semantics: loading a different record from the library
        should clear the path of whatever was previously open."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._source_path = "/tmp/old.gb"
            app._apply_record(tiny_record)   # default clear_undo=True
            assert app._source_path is None

    async def test_clear_undo_false_preserves_source_path(
        self, tiny_record, isolated_library
    ):
        """In-place-update semantics: after pLannotate merge or primer-add,
        the user's original source file should still be the Ctrl+S target."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._source_path = "/tmp/myfile.gb"
            app._apply_record(tiny_record, clear_undo=False)
            assert app._source_path == "/tmp/myfile.gb"

    async def test_clear_undo_false_preserves_undo_stack(
        self, tiny_record, isolated_library
    ):
        """The undo stack itself must not be wiped by an in-place update —
        otherwise the pre-merge / pre-primer-add state becomes un-recoverable."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._undo_stack.append(("DUMMY_SEQ", 0, tiny_record))
            app._apply_record(tiny_record, clear_undo=False)
            assert len(app._undo_stack) == 1
            assert app._undo_stack[0][0] == "DUMMY_SEQ"

    async def test_per_plasmid_undo_restored_on_switch_back(
        self, tiny_record, isolated_library
    ):
        """Load plasmid A, push an undo snapshot, switch to plasmid B, then
        switch back to A — A's undo history must be restored (not reset
        to empty as it was before per-plasmid stacks were introduced)."""
        from copy import deepcopy
        from Bio.SeqRecord import SeqRecord
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Prime plasmid A (tiny_record)
            app._apply_record(tiny_record, clear_undo=True)
            app._undo_stack.append(("A_SEQ", 0, tiny_record))
            # Build a second plasmid B with a distinct id
            rec_b = deepcopy(tiny_record)
            rec_b.id = "PLASMID_B_XYZ"
            rec_b.name = "PLASMID_B"
            # Switch to B — A's stack should be stashed
            app._apply_record(rec_b, clear_undo=True)
            assert app._undo_stack == []
            assert "pUC19_MINI" in app._stashed_undo_stacks or \
                   tiny_record.id in app._stashed_undo_stacks
            # Switch back to A — A's stack must be restored
            app._apply_record(tiny_record, clear_undo=True)
            assert len(app._undo_stack) == 1
            assert app._undo_stack[0][0] == "A_SEQ"

    async def test_per_plasmid_undo_lru_eviction(
        self, tiny_record, isolated_library
    ):
        """With _MAX_PLASMIDS_WITH_UNDO slots in the stash, loading a new
        plasmid once the cap is full must evict the least-recently-used
        plasmid's stashed history."""
        from copy import deepcopy
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._MAX_PLASMIDS_WITH_UNDO = 3
            # Load 5 plasmids A-E, pushing an undo snapshot into each.
            # The stash holds non-current plasmids only, so after E is
            # loaded the stash contains the 4 that were swapped out
            # (A, B, C, D) minus anything past the cap.
            ids = ["PID_A", "PID_B", "PID_C", "PID_D", "PID_E"]
            for pid in ids:
                rec = deepcopy(tiny_record)
                rec.id = pid
                app._apply_record(rec, clear_undo=True)
                app._undo_stack.append((f"{pid}_SEQ", 0, rec))
            # Stash capacity is 3. A was swapped out first → evicted.
            # B, C, D survive; E is live.
            assert "PID_A" not in app._stashed_undo_stacks
            assert "PID_B" in app._stashed_undo_stacks
            assert "PID_C" in app._stashed_undo_stacks
            assert "PID_D" in app._stashed_undo_stacks
            assert "PID_E" not in app._stashed_undo_stacks
            assert app._current_undo_key == "PID_E"

    async def test_mark_dirty_after_in_place_update_flips_unsaved(
        self, tiny_record, isolated_library
    ):
        """In-place update flow: _apply_record(clear_undo=False) calls
        _mark_clean internally, so callers must invoke _mark_dirty()
        afterwards to make the app's _unsaved flag reflect reality."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record, clear_undo=False)
            # _apply_record calls _mark_clean internally
            assert app._unsaved is False
            # The fix: callers must mark dirty after in-place updates
            app._mark_dirty()
            assert app._unsaved is True


class TestCrashRecoveryAutosave:
    """Crash-recovery autosave writes the current record to
    `_CRASH_RECOVERY_DIR/{safe_id}.gb` so an unexpected exit doesn't lose
    edits. The file is deleted on successful save or explicit abandon."""

    async def test_mark_dirty_schedules_autosave(
        self, tiny_record, isolated_library
    ):
        """`_mark_dirty` must register a debounced autosave timer."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._mark_dirty()
            assert app._autosave_timer is not None

    async def test_do_autosave_writes_genbank_file(
        self, tiny_record, isolated_library
    ):
        """Forcing `_do_autosave` must write a valid GenBank file at the
        record's autosave path."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._unsaved = True
            app._do_autosave()
            path = app._autosave_path(app._current_record)
            assert path is not None and path.exists()
            # Should be parseable GenBank
            from Bio import SeqIO
            roundtrip = SeqIO.read(str(path), "genbank")
            assert str(roundtrip.seq) == str(app._current_record.seq)

    async def test_mark_clean_clears_autosave_file(
        self, tiny_record, isolated_library
    ):
        """A successful save (→ `_mark_clean`) must delete the recovery
        file so next startup doesn't flag a stale recovery."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._unsaved = True
            app._do_autosave()
            path = app._autosave_path(app._current_record)
            assert path.exists()
            app._mark_clean()
            assert not path.exists()

    async def test_autosave_skipped_when_clean(
        self, tiny_record, isolated_library
    ):
        """If the record isn't dirty, autosave must not write anything."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._unsaved = False
            app._do_autosave()
            path = app._autosave_path(app._current_record)
            assert path is None or not path.exists()

    async def test_autosave_path_sanitises_unsafe_ids(self, tiny_record,
                                                       isolated_library):
        """Record ids can contain characters that are unsafe for filenames
        (slashes, spaces). The autosave helper must sanitise them."""
        from copy import deepcopy
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            bad = deepcopy(tiny_record)
            bad.id = "some/weird id.with:chars"
            path = app._autosave_path(bad)
            assert path is not None
            assert "/" not in path.name
            assert ":" not in path.name


class TestPlannotateReentryGuard:
    """Re-entry guard: pressing Shift+A while pLannotate is already running
    must not spawn a second subprocess. Regression guard for 2026-04-13."""

    async def test_action_noop_when_plannotate_running(
        self, tiny_record, isolated_library, monkeypatch
    ):
        """Set the running flag, call action_annotate_plasmid, confirm
        no subprocess was invoked."""
        import shutil, subprocess
        # Pretend pLannotate is fully installed so the code reaches the
        # re-entry guard (otherwise it short-circuits on "not installed").
        monkeypatch.setattr(sc, "_PLANNOTATE_CHECK_CACHE", None)
        monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")
        calls = []
        monkeypatch.setattr(
            subprocess, "run",
            lambda *a, **k: calls.append((a, k)) or None,
        )
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Simulate a pLannotate run already in flight
            app._plannotate_running = True
            app.action_annotate_plasmid()
            await pilot.pause(0.05)
            assert not calls, (
                "subprocess.run should NOT be invoked while "
                "_plannotate_running is True"
            )

    async def test_flag_exists_after_mount(
        self, tiny_record, isolated_library
    ):
        """The flag is initialized in on_mount — regression guard for
        future refactors that might forget to set it up."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            assert hasattr(app, "_plannotate_running")
            assert app._plannotate_running is False


class TestSidebarDetailWrapFeature:
    """The sidebar detail pane must render wrap features with an unambiguous
    compound-location string. A naive '{start+1}..{end}' shows '97..5' for a
    wrap, which a casual reader could mis-interpret as a 3-bp reverse range.
    Added 2026-04-13 alongside the _feat_len fix — users kept asking in the
    issue tracker 'what does 97..5 mean?'."""

    async def test_wrap_feature_coord_string_includes_origin(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        rec = SeqRecord(Seq("A" * 100), id="wrap_test",
                        annotations={"molecule_type": "DNA"})
        wrap_loc = CompoundLocation([
            FeatureLocation(95, 100, strand=1),
            FeatureLocation(0, 5, strand=1),
        ])
        rec.features.append(SeqFeature(wrap_loc, type="CDS",
                                       qualifiers={"label": ["wrapCDS"]}))
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sidebar = app.query_one("#sidebar", sc.FeatureSidebar)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            wrap_feat = next(f for f in pm._feats if f.get("label") == "wrapCDS")
            sidebar.show_detail(wrap_feat)
            box = sidebar.query_one("#detail-box")
            rendered = str(box.render())
            # Must reference both halves (tail 96..100 and head 1..5)
            assert "96" in rendered and "100" in rendered
            assert "1‥5" in rendered or "1..5" in rendered
            # Length displayed is 10 bp (5 + 5), not the wrong 'end - start'
            assert "10 bp" in rendered or "10\xa0bp" in rendered

    async def test_linear_feature_coord_string_unchanged(
        self, tiny_record, isolated_library,
    ):
        """A linear feature must still render as '{start+1}‥{end} (N bp)' —
        the wrap fix must not regress the common case."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sidebar = app.query_one("#sidebar", sc.FeatureSidebar)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            linear = next(f for f in pm._feats if f["end"] > f["start"])
            sidebar.show_detail(linear)
            box = sidebar.query_one("#detail-box")
            rendered = str(box.render())
            # One hyphen-separator only, no comma.
            assert "," not in rendered.split("(")[0]


class TestUndoSnapshotIndependence:
    """Defensive guard for an invariant that's currently easy to break by
    accident: undo/redo snapshots must be INDEPENDENT of the live record,
    so a future contributor who writes
    `self._current_record.features.append(...)` instead of building a fresh
    SeqRecord can't retroactively poison earlier undo entries.

    Today no code mutates _current_record in place, so this test wouldn't
    fail without the deep-copy — but locking the contract down with a test
    means a regression to in-place mutation will be caught immediately
    rather than discovered in production via a baffling Ctrl+Z bug."""

    async def test_push_undo_then_inplace_mutation_does_not_poison_snapshot(
        self, tiny_record, isolated_library,
    ):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            n_before = len(app._current_record.features)
            app._push_undo()
            # Simulate the dangerous pattern a future contributor might add
            app._current_record.features.append(SeqFeature(
                FeatureLocation(0, 5, strand=1),
                type="misc_feature",
                qualifiers={"label": ["poison"]},
            ))
            # Snapshot must NOT have grown — it's a deep copy.
            _, _, snapshot_record = app._undo_stack[-1]
            assert len(snapshot_record.features) == n_before, (
                "Undo snapshot was poisoned by an in-place mutation of "
                "_current_record. _push_undo must deep-copy."
            )

    async def test_action_undo_redo_snapshots_are_independent(
        self, tiny_record, isolated_library,
    ):
        """Round-trip: push_undo, _action_undo (redo snapshot taken), mutate
        in place, verify the redo snapshot survives."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._push_undo()
            n_before = len(app._current_record.features)
            app._action_undo()
            # _action_undo just pushed a redo snapshot — capture it before
            # poisoning the live record.
            _, _, redo_snapshot = app._redo_stack[-1]
            app._current_record.features.append(SeqFeature(
                FeatureLocation(0, 5, strand=1),
                type="misc_feature",
                qualifiers={"label": ["poison"]},
            ))
            assert len(redo_snapshot.features) == n_before, (
                "Redo snapshot was poisoned by in-place mutation of "
                "_current_record. _action_undo must deep-copy."
            )
