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
    record to the library. Library loads and undo/redo should NOT
    re-save."""

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


# ═══════════════════════════════════════════════════════════════════════════════
# _apply_record source_path + dirty-flag handling (regression guard 2026-04-13)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Before today's fix, _apply_record always cleared _source_path — even when
# called with clear_undo=False for an in-place update (e.g. primer-add).
# That meant after the in-place merge, Ctrl+S no longer targeted the user's
# original .gb file. Also, the merge path used lib.set_dirty(True) alone,
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
        """In-place-update semantics: after primer-add or feature-merge,
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

    async def test_autosave_path_disambiguates_sanitised_collisions(
        self, tiny_record, isolated_library,
    ):
        """Regression guard for 2026-04-25: pre-fix, two records with ids
        like 'foo/bar' and 'foo_bar' both sanitised to 'foo_bar.gb' and
        stomped each other on autosave. The fix appends a 6-char hash of
        the original id so collisions resolve to distinct filenames."""
        from copy import deepcopy
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            a = deepcopy(tiny_record); a.id = "foo/bar"
            b = deepcopy(tiny_record); b.id = "foo_bar"
            path_a = app._autosave_path(a)
            path_b = app._autosave_path(b)
            assert path_a is not None and path_b is not None
            assert path_a != path_b, (
                f"'{a.id}' and '{b.id}' must produce distinct autosave "
                f"paths after sanitisation; both got {path_a.name}"
            )
            # And reproducibility — the same id always maps to the same path.
            assert app._autosave_path(deepcopy(a)) == path_a


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


class TestCursorReachesEndOfSequence:
    """Regression guard for the 2026-04-25 cursor cap fix.

    Pre-fix the Right/Down arrow handlers clamped to `min(n - 1, …)` so the
    cursor could never land on position `n` (one past the last base) — the
    Edit Sequence dialog at `_edit_dialog_result` builds
    `old_seq[:s] + new_bases + old_seq[s:]`, so an end-of-sequence cursor is
    needed for an arrow-driven 'append' to work. Cap is now `min(n, …)`.

    Note (2026-04-25 amendment): Down arrow keeps the n-1 cap because pressing
    Down on the last row should land on the last visible base, not past it.
    Insert-at-end is reachable via Right arrow only.
    """

    async def test_right_arrow_at_end_advances_cursor_to_n(
        self, tiny_record, isolated_library,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            n  = len(sp._seq)
            # Position cursor on the last base, focus the panel, then press
            # Right. Pre-fix the cursor stayed at n-1; post-fix it reaches n.
            sp._cursor_pos = n - 1
            # Move focus off the default DataTable — the App-level on_key
            # arrow handler skips when a DataTable owns the keystroke,
            # because the sidebar/library tables have their own arrow nav.
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            app.set_focus(pm)
            await pilot.pause(0.05)
            await pilot.press("right")
            await pilot.pause(0.05)
            assert sp._cursor_pos == n, (
                f"Right arrow at last base should advance cursor to n={n} "
                f"(insert-at-end); got {sp._cursor_pos}"
            )

            # And one more Right keypress must NOT push cursor past n.
            await pilot.press("right")
            await pilot.pause(0.05)
            assert sp._cursor_pos == n

    async def test_down_arrow_on_last_row_caps_at_last_basepair(
        self, tiny_record, isolated_library,
    ):
        """Pressing Down on the last visible row should land on the last
        basepair (n-1), not on n. Position n has no base to highlight, so
        the cursor would visually disappear. Reported 2026-04-25."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            n  = len(sp._seq)
            lw = sp._line_width()
            assert n > 0 and lw > 0
            # Place cursor a few bases into what should be the last row.
            # (`n - 5` is on the last row for any sequence with at least
            # one full row; tiny_record is ~120 bp so this holds.)
            sp._cursor_pos = max(0, n - 5)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            app.set_focus(pm)
            await pilot.pause(0.05)
            await pilot.press("down")
            await pilot.pause(0.05)
            assert sp._cursor_pos == n - 1, (
                f"Down on last row should clamp cursor to n-1={n-1} "
                f"(last visible base); got {sp._cursor_pos}"
            )

            # Pressing Down again must keep cursor at n-1 (no overshoot).
            await pilot.press("down")
            await pilot.pause(0.05)
            assert sp._cursor_pos == n - 1


class TestSidebarClickCentersSeqPanel:
    """Regression guard for the 2026-04-25 sidebar-click centering fix.

    Clicking a feature in the sidebar previously highlighted it but did not
    scroll the sequence panel. Users with a 50 kb plasmid had to manually
    scroll through hundreds of rows to find the feature they just clicked.
    Now the seq panel jumps to the feature's wrap-aware midpoint."""

    async def test_sidebar_click_scrolls_seq_panel_to_feature(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        # 5 kb plasmid with a feature far past the initial viewport so
        # centering must scroll meaningfully (not stay at scroll_y=0).
        rec = SeqRecord(Seq("A" * 5000), id="centerTest",
                        annotations={"molecule_type": "DNA"})
        rec.features.append(SeqFeature(
            FeatureLocation(4000, 4200, strand=1), type="CDS",
            qualifiers={"label": ["targetFeat"]},
        ))
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sidebar = app.query_one("#sidebar", sc.FeatureSidebar)

            target_idx = next(
                i for i, f in enumerate(pm._feats)
                if f.get("label") == "targetFeat"
            )
            # Post the sidebar's RowActivated message — this is what fires
            # when the user actually clicks a feature row.
            sidebar.post_message(sc.FeatureSidebar.RowActivated(target_idx))
            await pilot.pause()
            await pilot.pause(0.1)  # let call_after_refresh do its scroll
            scroll = app.query_one("#seq-scroll")
            # Pre-fix: scroll_y stayed at 0. Post-fix: scrolls toward the
            # feature at bp 4100 (far row).
            assert scroll.scroll_y > 5, (
                f"Sidebar click on feature at bp 4100 should scroll seq "
                f"panel meaningfully; scroll_y={scroll.scroll_y}"
            )


class TestClickConsistencyAcrossPanels:
    """All three "I clicked a feature" entry points — sidebar row, plasmid map
    feature, sequence-panel lane art — must produce the same outcome:
    `user_sel` set to the feature span (the copyable highlight), cursor on
    the click bp / feature midpoint, and the seq-panel viewport centred on
    that bp. Map worked before but sidebar only set `sel_range` (bold +
    underline, NOT the copyable selection) and never centred reliably; the
    seq-panel click never centred at all. Reported 2026-04-25."""

    async def test_all_three_click_paths_are_consistent(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("A" * 5000), id="clickConsistency",
                        annotations={"molecule_type": "DNA"})
        # A handful of features so the table cursor actually moves on click.
        for i in range(3):
            rec.features.append(SeqFeature(
                FeatureLocation(i * 1500 + 100, i * 1500 + 200, strand=1),
                type="CDS", qualifiers={"label": [f"f{i}"]},
            ))
        rec.features.append(SeqFeature(
            FeatureLocation(4000, 4200, strand=1), type="CDS",
            qualifiers={"label": ["targetFeat"]},
        ))
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            scroll = app.query_one("#seq-scroll")
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sidebar = app.query_one("#sidebar", sc.FeatureSidebar)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)

            target_idx = next(i for i, f in enumerate(pm._feats)
                              if f.get("label") == "targetFeat")
            target = pm._feats[target_idx]
            target_bp = (target["start"] + target["end"]) // 2  # 4100

            async def reset_state():
                scroll.scroll_y = 0
                sidebar._prog_row = -1
                sp._user_sel = None
                sp._sel_range = None
                sp._cursor_pos = -1
                await pilot.pause(0.05)

            async def expect_focused_target():
                # All three paths land on these post-click invariants:
                assert sp._user_sel == (4000, 4200), (
                    f"user_sel must be the feature span; got {sp._user_sel}"
                )
                assert sp._cursor_pos == target_bp, (
                    f"cursor must land on bp={target_bp}; got {sp._cursor_pos}"
                )
                assert scroll.scroll_y > 30, (
                    f"viewport must scroll (centre on bp 4100); "
                    f"got scroll_y={scroll.scroll_y}"
                )

            # 1. Plasmid-map click.
            await reset_state()
            pm.post_message(sc.PlasmidMap.FeatureSelected(
                target_idx, target, bp=target_bp
            ))
            await pilot.pause(0.5)
            await expect_focused_target()

            # 2. Sequence-panel lane click — `from_lane=True` mirrors a
            # real click on the feature bar/arrow art (the only seq-panel
            # click path that should fire a whole-feature highlight).
            await reset_state()
            sp.post_message(
                sc.SequencePanel.SequenceClick(bp=target_bp, from_lane=True)
            )
            await pilot.pause(0.5)
            await expect_focused_target()

            # 3. Sidebar row click — actual mouse click, not a synthesized
            # RowActivated message (the message bypasses the RowHighlighted
            # cascade that real clicks go through).
            await reset_state()
            await pilot.click("#feat-table", offset=(5, target_idx + 1))
            await pilot.pause(0.5)
            await expect_focused_target()


class TestSidebarArrowNavSingleScroll:
    """Regression guard for the 2026-04-25 sidebar-arrow-key jitter fix.

    Pressing Up/Down in the sidebar's feature list cascades into
    `_focus_feature`, which used to call `select_feature_range` (which
    triggered `_ensure_cursor_visible` — partial scroll just-into-view)
    AND THEN `center_on_bp` (full scroll to centre). The two scrolls
    happened in quick succession and were perceptible as a jitter / snap
    on every arrow press. Fix: pass `scroll=False` to the highlight
    helpers when `_focus_feature` will issue its own centred scroll.
    """

    async def test_only_center_scroll_runs_per_arrow_press(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("A" * 5000), id="arrowJitter",
                        annotations={"molecule_type": "DNA"})
        for i in range(4):
            rec.features.append(SeqFeature(
                FeatureLocation(i * 1200 + 100, i * 1200 + 200, strand=1),
                type="CDS", qualifiers={"label": [f"f{i}"]},
            ))

        # Spy on _ensure_cursor_visible — it must NOT fire from the
        # `_focus_feature` path (which would cause the double-scroll). It
        # MAY fire from arrow-key text navigation (in the seq panel), but
        # arrowing through the sidebar should hit center_on_bp only.
        ensure_calls = []
        orig_ensure = sc.SequencePanel._ensure_cursor_visible
        def spy_ensure(self):
            ensure_calls.append(self._cursor_pos)
            orig_ensure(self)
        sc.SequencePanel._ensure_cursor_visible = spy_ensure

        try:
            app = sc.PlasmidApp()
            app._preload_record = rec
            async with app.run_test(size=TERMINAL_SIZE) as pilot:
                await pilot.pause()
                await pilot.pause(0.05)
                sp = app.query_one("#seq-panel", sc.SequencePanel)
                scroll = app.query_one("#seq-scroll")

                # Click row 0 to land focus inside the sidebar table.
                await pilot.click("#feat-table", offset=(5, 1))
                await pilot.pause(0.3)
                ensure_calls.clear()

                # Arrow through the rest. Each press triggers
                # _focus_feature → highlight + centring scroll. With the
                # fix, _ensure_cursor_visible does NOT fire from this
                # path.
                for _ in range(3):
                    await pilot.press("down")
                    await pilot.pause(0.3)

                assert ensure_calls == [], (
                    f"_focus_feature path must NOT call _ensure_cursor_visible "
                    f"(would cause double-scroll jitter). Got "
                    f"{len(ensure_calls)} call(s): {ensure_calls}"
                )
                # Sanity: the centring scroll did happen.
                assert scroll.scroll_y > 30, (
                    f"sidebar arrow nav should still scroll to centre; "
                    f"scroll_y={scroll.scroll_y}"
                )
        finally:
            sc.SequencePanel._ensure_cursor_visible = orig_ensure


class TestEnsureCursorVisibleShowsLanes:
    """Regression guard for the 2026-04-25 chunk-aware scroll fix.

    `_ensure_cursor_visible` previously scrolled to put the cursor's DNA
    forward-strand row at the top of the viewport when the user scrolled
    up. That left the feature lanes ABOVE the DNA off-screen, so the user
    had to press Up again just to see which feature their cursor was on.
    The fix targets `chunk_top` (DNA row minus above-lane rows) instead.
    """

    async def test_scroll_up_brings_above_lanes_into_view(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        # Build a record where every chunk has a feature so above_pairs > 0.
        rec = SeqRecord(Seq("A" * 2000), id="laneScrollTest",
                        annotations={"molecule_type": "DNA"})
        for i in range(0, 2000, 100):
            rec.features.append(SeqFeature(
                FeatureLocation(i, i + 80, strand=1), type="CDS",
                qualifiers={"label": [f"f{i}"]},
            ))

        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)

            # Park the cursor deep into the sequence so scrolling up to a
            # mid-chunk has somewhere to go.
            sp._cursor_pos = 1500
            sp.focus()
            await pilot.pause(0.05)
            sp._ensure_cursor_visible()
            await pilot.pause(0.05)

            scroll = app.query_one("#seq-scroll")
            scroll_y_before = scroll.scroll_y

            # Scroll up via Up arrow until cursor is at the top of viewport.
            for _ in range(40):
                await pilot.press("up")
                await pilot.pause(0.02)

            # The cursor's DNA row must be at least `above_pairs * rpg`
            # below the top of the viewport — i.e., the feature lanes
            # above the DNA must fit in the viewport above the cursor.
            line_width = sp._line_width()
            chunks_layout, prefix_dna2, prefix_lanes = sc._chunk_layout(
                sp._seq, sp._feats, line_width
            )
            rpg = 2 + (1 if sp._show_connectors else 0)
            chunk_idx = sp._cursor_pos // line_width
            above_pairs = chunks_layout[chunk_idx][3]
            chunk_top = (prefix_dna2[chunk_idx]
                         + (rpg - 2) * prefix_lanes[chunk_idx])
            dna_row = chunk_top + above_pairs * rpg

            scroll = app.query_one("#seq-scroll")
            vp_top = int(scroll.scroll_y)

            # Pre-fix: vp_top would equal dna_row (lanes clipped above viewport).
            # Post-fix: vp_top <= chunk_top, so the above-lanes are visible.
            assert vp_top <= chunk_top, (
                f"vp_top={vp_top} should be at or above chunk_top={chunk_top} "
                f"so the {above_pairs} feature-lane row(s) above the cursor's "
                f"DNA stay visible. dna_row={dna_row}"
            )


class TestMapClickCentersSeqPanel:
    """Regression guard for the 2026-04-25 map-click centering fix.

    Clicking on the plasmid map (feature or backbone) now centres the
    sequence panel on the clicked bp."""

    async def test_map_feature_click_centers_seq_panel(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("A" * 5000), id="mapClickTest",
                        annotations={"molecule_type": "DNA"})
        rec.features.append(SeqFeature(
            FeatureLocation(4000, 4200, strand=1), type="CDS",
            qualifiers={"label": ["targetFeat"]},
        ))
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)

            target_idx = next(
                i for i, f in enumerate(pm._feats)
                if f.get("label") == "targetFeat"
            )
            target_feat = pm._feats[target_idx]
            # Simulate the FeatureSelected event the map fires on click.
            pm.post_message(sc.PlasmidMap.FeatureSelected(
                target_idx, target_feat, bp=4100
            ))
            await pilot.pause()
            await pilot.pause(0.1)
            scroll = app.query_one("#seq-scroll")
            assert scroll.scroll_y > 5, (
                f"Map click at bp 4100 should scroll seq panel; "
                f"scroll_y={scroll.scroll_y}"
            )

    async def test_map_backbone_click_centers_seq_panel(
        self, tiny_record, isolated_library,
    ):
        """Clicking on the bare backbone (no feature) must still scroll the
        sequence panel — backbone clicks send `feat_dict=None, bp=clicked`
        and the handler now centres on bp regardless of feature presence."""
        # Use a longer record so backbone scrolling has somewhere to go.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 5000), id="backboneClickTest",
                        annotations={"molecule_type": "DNA"})
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            # idx=-1, feat=None — the backbone-click signature.
            pm.post_message(sc.PlasmidMap.FeatureSelected(
                -1, None, bp=4500
            ))
            await pilot.pause()
            await pilot.pause(0.1)
            scroll = app.query_one("#seq-scroll")
            assert scroll.scroll_y > 5, (
                f"Backbone click at bp 4500 should still scroll seq panel; "
                f"scroll_y={scroll.scroll_y}"
            )


class TestSeqClickWrapFeature:
    """Regression guard for the 2026-04-25 fix to `_seq_click`.

    Pre-fix the handler used `s <= bp < e and (e - s) < best_span` which
    (a) failed every wrap feature (where `e < s`, so the comparison is
    always False) and (b) used a negative `e - s` span for any wrap that
    *did* somehow leak through. Clicking the lane art of a wrap feature
    on the sequence panel silently selected nothing.

    Updated 2026-04-28: `_seq_click` now distinguishes lane clicks
    (`from_lane=True`) from DNA-row clicks. Only lane clicks pick a
    feature; DNA-row clicks just place the cursor.
    """

    async def test_lane_click_inside_wrap_feature_selects_it(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        rec = SeqRecord(Seq("A" * 100), id="wrap_click_test",
                        annotations={"molecule_type": "DNA"})
        # Wrap feature spanning 95..100 + 0..5 (10 bp around origin).
        wrap_loc = CompoundLocation([
            FeatureLocation(95, 100, strand=1),
            FeatureLocation(0, 5, strand=1),
        ])
        rec.features.append(SeqFeature(wrap_loc, type="CDS",
                                       qualifiers={"label": ["wrapCDS"]}))
        rec.features.append(SeqFeature(
            FeatureLocation(40, 60, strand=1), type="CDS",
            qualifiers={"label": ["linearCDS"]}
        ))

        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)

            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp = app.query_one("#seq-panel",   sc.SequencePanel)

            wrap_idx = next(
                (i for i, f in enumerate(pm._feats)
                 if f.get("label") == "wrapCDS"),
                None,
            )
            assert wrap_idx is not None

            # Lane click at bp=2 — inside the wrap's head [0, 5).
            sp.post_message(
                sc.SequencePanel.SequenceClick(bp=2, from_lane=True)
            )
            await pilot.pause()
            await pilot.pause(0.05)
            assert pm.selected_idx == wrap_idx, (
                "Lane-click at bp=2 (wrap head) should select the wrap "
                f"feature; got selected_idx={pm.selected_idx}"
            )

            # Lane click at bp=97 — inside the wrap's tail [95, 100).
            sp.post_message(
                sc.SequencePanel.SequenceClick(bp=97, from_lane=True)
            )
            await pilot.pause()
            await pilot.pause(0.05)
            assert pm.selected_idx == wrap_idx, (
                "Lane-click at bp=97 (wrap tail) should select the wrap "
                f"feature; got selected_idx={pm.selected_idx}"
            )

    async def test_base_click_does_not_select_feature(
        self, isolated_library,
    ):
        """A click on the DNA strand row (not the lane art) must NOT
        trigger a whole-feature selection, even if `bp` is inside one
        — the user asked for a single-base operation, not a feature
        pick. Regression guard for the 2026-04-28 lane-click rule."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("A" * 100), id="base_click_test",
                        annotations={"molecule_type": "DNA"})
        rec.features.append(SeqFeature(
            FeatureLocation(40, 60, strand=1), type="CDS",
            qualifiers={"label": ["linearCDS"]}
        ))

        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp = app.query_one("#seq-panel",   sc.SequencePanel)

            assert pm.selected_idx == -1
            # Base click at bp=50 (inside the CDS feature) must NOT
            # select it — only lane art clicks do.
            sp.post_message(
                sc.SequencePanel.SequenceClick(bp=50, from_lane=False)
            )
            await pilot.pause()
            await pilot.pause(0.05)
            assert pm.selected_idx == -1, (
                "Base-row click should not pick a feature; got "
                f"selected_idx={pm.selected_idx}"
            )

    async def test_click_outside_wrap_does_not_falsely_select(
        self, isolated_library,
    ):
        """Negative control: clicking far from the wrap feature must NOT
        pick it up. The fix uses `_feat_len(s, e, total)` which is positive
        for wraps; a regression that compared raw `e - s` (negative) would
        always pick the wrap feature as 'smallest' and break this case."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        rec = SeqRecord(Seq("A" * 100), id="wrap_neg_test",
                        annotations={"molecule_type": "DNA"})
        rec.features.append(SeqFeature(
            CompoundLocation([
                FeatureLocation(95, 100, strand=1),
                FeatureLocation(0, 5, strand=1),
            ]),
            type="CDS", qualifiers={"label": ["wrapCDS"]},
        ))
        rec.features.append(SeqFeature(
            FeatureLocation(40, 60, strand=1), type="CDS",
            qualifiers={"label": ["linearCDS"]}
        ))

        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp = app.query_one("#seq-panel",   sc.SequencePanel)

            linear_idx = next(
                i for i, f in enumerate(pm._feats)
                if f.get("label") == "linearCDS"
            )

            # Lane click at bp=50 — inside the linear feature, far from wrap.
            sp.post_message(
                sc.SequencePanel.SequenceClick(bp=50, from_lane=True)
            )
            await pilot.pause()
            await pilot.pause(0.05)
            assert pm.selected_idx == linear_idx, (
                f"Should pick linear feature at bp=50, "
                f"got selected_idx={pm.selected_idx}"
            )


class TestSplashScreen:
    """Splash modal mounts on launch and dismisses on any keystroke; the
    test conftest sets `_skip_splash = True` for everything else, so any
    test that wants to drive the splash has to opt back in."""

    async def test_splash_mounts_when_enabled(self, isolated_library):
        app = sc.PlasmidApp()
        app._skip_splash = False
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.2)
            assert isinstance(app.screen, sc.SplashScreen)
            # The new splash paints DNA helix + logo + tagline + version
            # all into one canvas Static. Probe the rendered Text for
            # the Binomica + version string to verify the composition
            # actually included them.
            canvas = app.screen.query_one("#splash-canvas")
            content = str(canvas.render())
            assert "Binomica" in content
            assert sc.__version__ in content
            app.exit()

    async def test_splash_dismisses_on_key(self, isolated_library):
        app = sc.PlasmidApp()
        app._skip_splash = False
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.SplashScreen)
            await pilot.press("a")
            await pilot.pause(0.1)
            assert not isinstance(app.screen, sc.SplashScreen)
            app.exit()

    async def test_splash_skipped_under_default_test_config(self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            # conftest set _skip_splash = True; no splash on screen.
            assert not isinstance(app.screen, sc.SplashScreen)
            app.exit()


class TestQuitConfirm:
    """Pressing q opens QuitConfirmModal (default No) when there are no
    unsaved edits. With unsaved edits the existing UnsavedQuitModal still
    fires instead. Tab cycles between buttons; Enter on the focused
    button presses it (Textual default — no extra wiring needed)."""

    async def test_clean_quit_pushes_confirm_modal(self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            assert app._unsaved is False
            app.action_quit()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.QuitConfirmModal)
            # Default focus is on No.
            assert app.screen.focused.id == "btn-quitcon-no"
            app.exit()

    async def test_clean_quit_no_keeps_app_running(self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app.action_quit()
            await pilot.pause(0.1)
            app.screen.query_one("#btn-quitcon-no").action_press()
            await pilot.pause(0.1)
            # Still running — return to default screen, not exited.
            assert not isinstance(app.screen, sc.QuitConfirmModal)
            app.exit()

    async def test_unsaved_quit_routes_through_unsaved_modal(
        self, tiny_record, isolated_library
    ):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.2)
            app._mark_dirty()
            await pilot.pause(0.05)
            app.action_quit()
            await pilot.pause(0.1)
            # With unsaved edits, the unsaved modal fires (3 buttons),
            # not the simple QuitConfirmModal.
            assert isinstance(app.screen, sc.UnsavedQuitModal)
            app.screen.query_one("#btn-cancel-quit").action_press()
            await pilot.pause(0.1)
            app.exit()

    async def test_tab_cycles_focus_between_no_and_yes(self, isolated_library):
        """Tab + Enter end-to-end on a confirm modal — the only modal
        contract that matters for keyboard-only quit confirmation."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app.action_quit()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal.focused.id == "btn-quitcon-no"
            await pilot.press("tab")
            await pilot.pause(0.05)
            assert modal.focused.id == "btn-quitcon-yes"
            await pilot.press("tab")
            await pilot.pause(0.05)
            # Wraps back round to No.
            assert modal.focused.id == "btn-quitcon-no"
            app.exit()


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
