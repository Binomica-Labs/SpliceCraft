"""
test_collections — plasmid collection persistence + CollectionsModal flows.

Covers:
  * `_load_collections` / `_save_collections` envelope round-trip + legacy
    bare-list back-compat (matches sacred invariant #7).
  * CollectionsModal: save current library as a new collection, load
    (replace), delete, name-collision rejection.
  * Auto-load of the first library entry on PlasmidApp startup.
"""
from __future__ import annotations

import json

import pytest
import splicecraft as sc


# ── Persistence round-trip ─────────────────────────────────────────────────────

class TestCollectionsPersistence:
    def test_empty_load(self):
        # Fresh tmp file → no collections yet.
        assert sc._load_collections() == []

    def test_save_load_round_trip(self):
        sample = [
            {"name": "yeast", "description": "S. cerevisiae plasmids",
             "plasmids": [{"id": "P1", "name": "yp1", "size": 100,
                           "gb_text": "FAKE"}]},
            {"name": "ecoli", "description": "E. coli toolkit",
             "plasmids": []},
        ]
        sc._save_collections(sample)
        sc._collections_cache = None  # force a cold reload from disk
        out = sc._load_collections()
        assert len(out) == 2
        assert out[0]["name"] == "yeast"
        assert out[0]["plasmids"][0]["id"] == "P1"
        assert out[1]["name"] == "ecoli"

    def test_envelope_schema_on_disk(self):
        sc._save_collections([{"name": "t", "plasmids": []}])
        raw = json.loads(sc._COLLECTIONS_FILE.read_text())
        assert raw.get("_schema_version") == 1
        assert isinstance(raw.get("entries"), list)

    def test_legacy_bare_list_loads(self):
        """Pre-envelope save format must still load (sacred invariant #7)."""
        sc._COLLECTIONS_FILE.write_text(
            json.dumps([{"name": "legacy", "plasmids": []}])
        )
        sc._collections_cache = None
        out = sc._load_collections()
        assert len(out) == 1
        assert out[0]["name"] == "legacy"

    def test_non_dict_entries_dropped(self):
        sc._COLLECTIONS_FILE.write_text(json.dumps({
            "_schema_version": 1,
            "entries": [{"name": "ok", "plasmids": []}, "garbage", 42, None],
        }))
        sc._collections_cache = None
        out = sc._load_collections()
        assert len(out) == 1
        assert out[0]["name"] == "ok"


# ── Modal flows ────────────────────────────────────────────────────────────────

class TestCollectionsModalFlows:
    async def test_save_current_library_as_collection(self):
        # Seed a non-empty library. App startup auto-creates a "Main
        # Collection" wrapping it, then the modal Save adds "myset"
        # alongside.
        sc._save_library([{"id": "X", "name": "X", "size": 10,
                           "gb_text": "GB"}])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.CollectionsModal())
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#coll-save-name").value = "myset"
            modal.query_one("#btn-coll-save").action_press()
            await pilot.pause(0.2)
            colls = sc._load_collections()
            myset = [c for c in colls if c.get("name") == "myset"]
            assert len(myset) == 1
            assert len(myset[0]["plasmids"]) == 1
            assert myset[0]["plasmids"][0]["id"] == "X"
            # Main Collection is also present — auto-created on app start.
            assert any(c.get("name") == sc._DEFAULT_COLLECTION_NAME for c in colls)
            app.exit()

    async def test_save_rejects_duplicate_name(self):
        sc._save_collections([{"name": "dup", "plasmids": []}])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.CollectionsModal())
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#coll-save-name").value = "dup"
            modal.query_one("#btn-coll-save").action_press()
            await pilot.pause(0.2)
            # Still only the original collection — no new row.
            assert len(sc._load_collections()) == 1
            app.exit()

    async def test_save_rejects_blank_name(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.CollectionsModal())
            await pilot.pause(0.3)
            modal = app.screen
            before = sc._load_collections()
            modal.query_one("#coll-save-name").value = "   "
            modal.query_one("#btn-coll-save").action_press()
            await pilot.pause(0.2)
            # Same set as before — blank name was rejected. Main Collection
            # may have been auto-created on startup but no new entry should
            # have been added.
            after = sc._load_collections()
            assert [c.get("name") for c in after] == [c.get("name") for c in before]
            app.exit()

    async def test_load_replaces_library(self):
        # Library starts with one record; collection holds two different
        # records — Load must end with the library == collection's plasmids.
        sc._save_library([{"id": "OLD", "name": "old", "size": 1,
                           "gb_text": "GB"}])
        sc._save_collections([{
            "name": "newset",
            "plasmids": [
                {"id": "A", "name": "a", "size": 100, "gb_text": "GB"},
                {"id": "B", "name": "b", "size": 200, "gb_text": "GB"},
            ],
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.CollectionsModal())
            await pilot.pause(0.3)
            modal = app.screen
            t = modal.query_one("#coll-table")
            t.move_cursor(row=0)
            modal.query_one("#btn-coll-load").action_press()
            await pilot.pause(0.3)
            lib = sc._load_library()
            ids = {e.get("id") for e in lib}
            assert ids == {"A", "B"}
            app.exit()

    async def test_delete_collection(self):
        sc._save_collections([
            {"name": "keep", "plasmids": []},
            {"name": "drop", "plasmids": []},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.CollectionsModal())
            await pilot.pause(0.3)
            modal = app.screen
            t = modal.query_one("#coll-table")
            # Cursor row 1 = "drop" (collection order on save).
            t.move_cursor(row=1)
            modal.query_one("#btn-coll-del").action_press()
            await pilot.pause(0.2)
            names = {c["name"] for c in sc._load_collections()}
            assert names == {"keep"}
            app.exit()


# ── Auto-load first library entry on startup ──────────────────────────────────

class TestStartupAutoLoad:
    """If the library is non-empty, PlasmidApp.on_mount must call
    _apply_record on the first entry so the canvas isn't blank."""

    @pytest.fixture(autouse=True)
    def _no_seed(self, monkeypatch):
        # Suppress the first-run seed worker so it can't mask our assertions.
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

    async def test_first_entry_loaded(self, tiny_record):
        from io import StringIO
        from Bio import SeqIO
        buf = StringIO()
        SeqIO.write(tiny_record, buf, "genbank")
        sc._save_library([{
            "id":   tiny_record.id,
            "name": tiny_record.name,
            "size": len(tiny_record.seq),
            "gb_text": buf.getvalue(),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await pilot.pause(0.3)
            assert app._current_record is not None
            assert app._current_record.id == tiny_record.id
            app.exit()

    async def test_empty_library_no_record(self):
        # Library empty + seed-worker patched out → canvas stays blank.
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await pilot.pause(0.3)
            assert app._current_record is None
            app.exit()

    async def test_unparseable_first_entry_does_not_crash(self):
        sc._save_library([{
            "id":   "BROKEN",
            "name": "broken",
            "size": 0,
            "gb_text": "this is not GenBank text",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await pilot.pause(0.3)
            # App still up; current record falls back to None.
            assert app._current_record is None
            app.exit()


# ── Main Collection migration ─────────────────────────────────────────────────

class TestMainCollectionMigration:
    """First-run behaviour: an existing library gets wrapped in
    'Main Collection' the first time the new code sees it. Active
    collection is auto-pointed at it so library writes start mirroring
    immediately."""

    @pytest.fixture(autouse=True)
    def _no_seed(self, monkeypatch):
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

    def test_ensure_default_creates_main_with_existing_library(self):
        # Library has content, collections.json doesn't exist yet.
        sc._save_library([{"id": "A", "name": "a", "size": 1, "gb_text": "GB"}])
        assert sc._load_collections() == []
        sc._ensure_default_collection()
        colls = sc._load_collections()
        assert len(colls) == 1
        assert colls[0]["name"] == sc._DEFAULT_COLLECTION_NAME
        assert colls[0]["plasmids"][0]["id"] == "A"
        assert sc._get_active_collection_name() == sc._DEFAULT_COLLECTION_NAME

    def test_ensure_default_idempotent(self):
        sc._save_collections([{"name": "Existing", "plasmids": []}])
        sc._set_active_collection_name("Existing")
        sc._ensure_default_collection()
        # Did not add Main Collection; active unchanged.
        assert [c["name"] for c in sc._load_collections()] == ["Existing"]
        assert sc._get_active_collection_name() == "Existing"

    def test_ensure_default_sets_active_when_missing(self):
        # Collections present but no active set → first one becomes active.
        sc._save_collections([{"name": "First", "plasmids": []},
                              {"name": "Second", "plasmids": []}])
        # Explicitly clear any active that may have been set by save_collections
        sc._set_active_collection_name(None)
        assert sc._get_active_collection_name() is None
        sc._ensure_default_collection()
        assert sc._get_active_collection_name() == "First"

    def test_app_startup_creates_main(self):
        # Empty everything → app boot should leave a Main Collection on disk.
        assert sc._load_collections() == []
        app = sc.PlasmidApp()
        async def go():
            async with app.run_test(size=(140, 50)) as pilot:
                await pilot.pause()
                await pilot.pause(0.1)
                colls = sc._load_collections()
                assert any(c.get("name") == sc._DEFAULT_COLLECTION_NAME
                           for c in colls)
                assert sc._get_active_collection_name() == sc._DEFAULT_COLLECTION_NAME
                app.exit()
        import asyncio
        asyncio.run(go())


# ── Library/collection sync ───────────────────────────────────────────────────

class TestLibraryCollectionSync:
    """Every library mutation must mirror into the active collection so
    the two on-disk files don't drift."""

    def test_save_library_mirrors_to_active(self):
        sc._save_collections([{"name": "Active", "plasmids": []}])
        sc._set_active_collection_name("Active")
        new_entries = [{"id": "P1", "name": "p1", "size": 100, "gb_text": "GB"}]
        sc._save_library(new_entries)
        # Active collection on disk now holds the same entries.
        coll = next(c for c in sc._load_collections()
                    if c.get("name") == "Active")
        assert len(coll["plasmids"]) == 1
        assert coll["plasmids"][0]["id"] == "P1"

    def test_save_library_silent_with_no_active(self):
        # No active collection → library writes don't crash; collection
        # file stays empty.
        sc._set_active_collection_name(None)
        sc._save_library([{"id": "X", "name": "x", "size": 1, "gb_text": "GB"}])
        # No collection should have been touched.
        assert sc._load_collections() == []

    def test_save_library_with_deleted_active_no_crash(self):
        # Active points at a name that no longer exists — should silently
        # skip the mirror, never raise.
        sc._save_collections([{"name": "Stale", "plasmids": []}])
        sc._set_active_collection_name("Stale")
        sc._save_collections([])  # delete it directly
        sc._save_library([{"id": "P", "name": "p", "size": 1, "gb_text": "GB"}])
        assert sc._load_collections() == []  # still empty


# ── LibraryPanel two-mode flows ──────────────────────────────────────────────

class TestLibraryPanelModes:
    """The redesigned panel toggles between collections and plasmids
    views. Click a collection → plasmids view; ← → collections view."""

    @pytest.fixture(autouse=True)
    def _no_seed(self, monkeypatch):
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

    async def test_starts_in_plasmids_view_when_active_set(self):
        sc._save_collections([{"name": "Pre", "plasmids": []}])
        sc._set_active_collection_name("Pre")
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            assert lib._view_mode == "plasmids"
            assert lib.query_one("#lib-table").display is True
            assert lib.query_one("#lib-coll-table").display is False
            app.exit()

    async def test_back_button_returns_to_collections_view(self):
        sc._save_collections([{"name": "X", "plasmids": []}])
        sc._set_active_collection_name("X")
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            assert lib._view_mode == "plasmids"
            lib.query_one("#btn-lib-back").action_press()
            await pilot.pause(0.1)
            assert lib._view_mode == "collections"
            assert lib.query_one("#lib-coll-table").display is True
            assert lib.query_one("#lib-table").display is False
            app.exit()

    async def test_clicking_collection_requires_double_activation(self):
        """Loading a collection swaps the entire library, so a stray
        RowSelected must NOT fire the load. The first activation arms
        the collection and a second activation on the same row commits.
        """
        sc._save_collections([
            {"name": "A", "plasmids": [
                {"id": "p1", "name": "p1", "size": 1, "gb_text": "GB"}]},
            {"name": "B", "plasmids": []},
        ])
        sc._set_active_collection_name(None)
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.move_cursor(row=0)  # "A"
            t.action_select_cursor()  # arm
            await pilot.pause(0.1)
            # Still in collections view — first click only armed.
            assert lib._view_mode == "collections"
            assert lib._coll_armed_name == "A"
            # Second activation on the same row commits.
            t.action_select_cursor()
            await pilot.pause(0.2)
            assert lib._view_mode == "plasmids"
            assert sc._get_active_collection_name() == "A"
            assert [e["id"] for e in sc._load_library()] == ["p1"]
            app.exit()

    async def test_switching_row_disarms_previous_collection(self):
        """Arming row A then activating row B should NOT load either —
        the arm transfers to B and one more activation on B is needed."""
        sc._save_collections([
            {"name": "A", "plasmids": [
                {"id": "pA", "name": "pA", "size": 1, "gb_text": "GB"}]},
            {"name": "B", "plasmids": [
                {"id": "pB", "name": "pB", "size": 1, "gb_text": "GB"}]},
        ])
        sc._set_active_collection_name(None)
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.move_cursor(row=0)
            t.action_select_cursor()  # arm A
            await pilot.pause(0.1)
            t.move_cursor(row=1)
            t.action_select_cursor()  # arm B (disarms A)
            await pilot.pause(0.1)
            # Still in collections view; B is now armed.
            assert lib._view_mode == "collections"
            assert lib._coll_armed_name == "B"
            t.action_select_cursor()  # commit B
            await pilot.pause(0.2)
            assert sc._get_active_collection_name() == "B"
            app.exit()


class TestPanelCollectionCRUD:
    """Add / remove / rename a collection from the panel."""

    @pytest.fixture(autouse=True)
    def _no_seed(self, monkeypatch):
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

    async def test_add_collection_via_panel(self):
        sc._save_collections([{"name": "Existing", "plasmids": []}])
        sc._set_active_collection_name(None)  # start in collections view
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            lib.query_one("#btn-coll-add").action_press()
            await pilot.pause(0.2)
            modal = app.screen
            assert isinstance(modal, sc.CollectionNameModal)
            modal.query_one("#collname-input").value = "MyNew"
            modal.query_one("#btn-collname-ok").action_press()
            await pilot.pause(0.2)
            names = [c["name"] for c in sc._load_collections()]
            assert "MyNew" in names
            app.exit()

    async def test_add_collection_rejects_duplicate(self):
        sc._save_collections([{"name": "Same", "plasmids": []}])
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib.query_one("#btn-coll-add").action_press()
            await pilot.pause(0.2)
            modal = app.screen
            modal.query_one("#collname-input").value = "Same"
            modal.query_one("#btn-collname-ok").action_press()
            await pilot.pause(0.2)
            # Still only one — duplicate rejected.
            same = [c for c in sc._load_collections() if c.get("name") == "Same"]
            assert len(same) == 1
            app.exit()

    async def test_rename_collection_via_panel(self):
        sc._save_collections([{"name": "OldName", "plasmids": []}])
        sc._set_active_collection_name("OldName")
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.move_cursor(row=0)
            lib.query_one("#btn-coll-rename").action_press()
            await pilot.pause(0.2)
            modal = app.screen
            assert isinstance(modal, sc.CollectionNameModal)
            modal.query_one("#collname-input").value = "NewName"
            modal.query_one("#btn-collname-ok").action_press()
            await pilot.pause(0.2)
            names = [c["name"] for c in sc._load_collections()]
            assert names == ["NewName"]
            # Active pointer follows the rename.
            assert sc._get_active_collection_name() == "NewName"
            app.exit()

    async def test_delete_collection_via_panel(self):
        sc._save_collections([
            {"name": "Keep", "plasmids": []},
            {"name": "Drop", "plasmids": []},
        ])
        sc._set_active_collection_name("Keep")
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.move_cursor(row=1)  # "Drop"
            lib.query_one("#btn-coll-del").action_press()
            await pilot.pause(0.2)
            # Stage 1: friendly confirm
            assert isinstance(app.screen, sc.CollectionDeleteConfirmModal)
            app.screen.query_one("#btn-colldel-yes").action_press()
            await pilot.pause(0.2)
            # Stage 2: scary red second confirm
            assert isinstance(app.screen, sc.ScaryDeleteConfirmModal)
            app.screen.query_one("#btn-scarydel-yes").action_press()
            await pilot.pause(0.2)
            names = [c["name"] for c in sc._load_collections()]
            assert names == ["Keep"]
            app.exit()

    async def test_collection_delete_first_no_keeps_collection(self):
        """Cancel at the friendly stage — collection stays."""
        sc._save_collections([{"name": "Stay", "plasmids": []}])
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.move_cursor(row=0)
            lib.query_one("#btn-coll-del").action_press()
            await pilot.pause(0.2)
            assert isinstance(app.screen, sc.CollectionDeleteConfirmModal)
            app.screen.query_one("#btn-colldel-no").action_press()
            await pilot.pause(0.2)
            assert [c["name"] for c in sc._load_collections()] == ["Stay"]
            app.exit()

    async def test_collection_delete_second_no_keeps_collection(self):
        """Yes through the friendly stage but No on the scary stage —
        the collection still survives. Belt + suspenders confirm pattern."""
        sc._save_collections([{"name": "Saved", "plasmids": []}])
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.move_cursor(row=0)
            lib.query_one("#btn-coll-del").action_press()
            await pilot.pause(0.2)
            app.screen.query_one("#btn-colldel-yes").action_press()
            await pilot.pause(0.2)
            assert isinstance(app.screen, sc.ScaryDeleteConfirmModal)
            app.screen.query_one("#btn-scarydel-no").action_press()
            await pilot.pause(0.2)
            assert [c["name"] for c in sc._load_collections()] == ["Saved"]
            app.exit()

    async def test_delete_active_collection_clears_active(self):
        sc._save_collections([{"name": "Going", "plasmids": []}])
        sc._set_active_collection_name("Going")
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.move_cursor(row=0)
            lib.query_one("#btn-coll-del").action_press()
            await pilot.pause(0.2)
            app.screen.query_one("#btn-colldel-yes").action_press()
            await pilot.pause(0.2)
            app.screen.query_one("#btn-scarydel-yes").action_press()
            await pilot.pause(0.2)
            assert sc._load_collections() == []
            assert sc._get_active_collection_name() is None
            app.exit()


class TestBackButtonUnsavedPrompt:
    """Pressing the Back button must check for unsaved edits and prompt
    the user. Save → save then leave; Discard → revert + leave; Cancel
    → stay in plasmids view."""

    @pytest.fixture(autouse=True)
    def _no_seed(self, monkeypatch):
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

    async def test_back_button_visible_in_plasmids_view(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            assert lib._view_mode == "plasmids"
            btn = lib.query_one("#btn-lib-back")
            # The back button lives in the bottom button row (it took
            # over the slot the old ◈ annotate button used). Label is
            # the universal "←" back arrow.
            assert lib.query_one("#lib-btns").display is True
            assert "←" in str(btn.label)
            app.exit()

    async def test_back_button_hidden_in_collections_view(self):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name(None)  # force collections view
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            # In collections view the whole plasmid-button row hides,
            # which carries the back button along with it.
            assert lib.query_one("#lib-btns").display is False
            app.exit()

    async def test_back_with_no_unsaved_goes_back_directly(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            assert app._unsaved is False
            lib = app.query_one("#library", sc.LibraryPanel)
            lib.query_one("#btn-lib-back").action_press()
            await pilot.pause(0.1)
            # No modal — went straight to collections view.
            assert lib._view_mode == "collections"
            assert type(app.screen).__name__ != "UnsavedNavigateModal"
            app.exit()

    async def test_back_with_unsaved_pushes_modal(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._mark_dirty()  # simulate an unsaved edit
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            assert lib._view_mode == "plasmids"
            lib.query_one("#btn-lib-back").action_press()
            await pilot.pause(0.2)
            # Modal pushed; still in plasmids view.
            assert isinstance(app.screen, sc.UnsavedNavigateModal)
            assert lib._view_mode == "plasmids"
            app.exit()

    async def test_modal_cancel_keeps_plasmids_view(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._mark_dirty()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib.query_one("#btn-lib-back").action_press()
            await pilot.pause(0.2)
            assert isinstance(app.screen, sc.UnsavedNavigateModal)
            app.screen.query_one("#btn-navunsv-cancel").action_press()
            await pilot.pause(0.1)
            assert lib._view_mode == "plasmids"
            assert app._unsaved is True  # state preserved
            app.exit()

    async def test_modal_discard_reverts_and_goes_back(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.2)  # let preload save to library
            app._mark_dirty()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib.query_one("#btn-lib-back").action_press()
            await pilot.pause(0.2)
            app.screen.query_one("#btn-navunsv-discard").action_press()
            await pilot.pause(0.2)
            # Switched view, dirty cleared.
            assert lib._view_mode == "collections"
            assert app._unsaved is False
            app.exit()

    async def test_modal_save_persists_and_goes_back(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.2)
            app._mark_dirty()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib.query_one("#btn-lib-back").action_press()
            await pilot.pause(0.2)
            app.screen.query_one("#btn-navunsv-save").action_press()
            await pilot.pause(0.3)
            # Saved → marked clean → switched to collections view.
            assert app._unsaved is False
            assert lib._view_mode == "collections"
            app.exit()


class TestPlasmidDeleteConfirmFlow:
    """The plasmid `−` button and the Delete key both go through
    LibraryDeleteConfirmModal (default focus on No) so a stray click
    or keypress can't silently nuke a saved plasmid."""

    @pytest.fixture(autouse=True)
    def _no_seed(self, monkeypatch):
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

    async def test_button_pushes_confirm_modal(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.2)  # let preload save
            lib = app.query_one("#library", sc.LibraryPanel)
            t = lib.query_one("#lib-table")
            t.move_cursor(row=0)
            lib.query_one("#btn-lib-del").action_press()
            await pilot.pause(0.2)
            assert isinstance(app.screen, sc.LibraryDeleteConfirmModal)
            # Default focus is on No; cancel dismisses without deleting.
            app.screen.query_one("#btn-libdel-no").action_press()
            await pilot.pause(0.2)
            assert any(e.get("id") == tiny_record.id
                       for e in sc._load_library())
            app.exit()

    async def test_button_yes_deletes(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.2)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib.query_one("#lib-table").move_cursor(row=0)
            lib.query_one("#btn-lib-del").action_press()
            await pilot.pause(0.2)
            app.screen.query_one("#btn-libdel-yes").action_press()
            await pilot.pause(0.2)
            assert all(e.get("id") != tiny_record.id
                       for e in sc._load_library())
            app.exit()

    async def test_delete_key_routes_through_same_confirm(self, tiny_record):
        sc._save_collections([{"name": "C", "plasmids": []}])
        sc._set_active_collection_name("C")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.2)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib.query_one("#lib-table").focus()
            lib.query_one("#lib-table").move_cursor(row=0)
            await pilot.pause(0.05)
            app.action_delete_feature()
            await pilot.pause(0.2)
            assert isinstance(app.screen, sc.LibraryDeleteConfirmModal)
            app.exit()


class TestCollectionDeleteKeyFlow:
    """The Delete key in collections view triggers the same two-stage
    confirm as the `−` button — friendly modal then loud red modal."""

    @pytest.fixture(autouse=True)
    def _no_seed(self, monkeypatch):
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

    async def test_delete_key_in_collections_view_pushes_first_confirm(self):
        sc._save_collections([{"name": "Doomed", "plasmids": []}])
        sc._set_active_collection_name(None)
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.focus()
            t.move_cursor(row=0)
            await pilot.pause(0.05)
            app.action_delete_feature()
            await pilot.pause(0.2)
            assert isinstance(app.screen, sc.CollectionDeleteConfirmModal)
            app.exit()

    async def test_delete_key_full_two_stage_flow(self):
        sc._save_collections([{"name": "Doomed", "plasmids": []}])
        sc._set_active_collection_name(None)
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "collections"
            lib._apply_view_mode()
            lib._repopulate_collections()
            t = lib.query_one("#lib-coll-table")
            t.focus()
            t.move_cursor(row=0)
            await pilot.pause(0.05)
            app.action_delete_feature()
            await pilot.pause(0.2)
            app.screen.query_one("#btn-colldel-yes").action_press()
            await pilot.pause(0.2)
            assert isinstance(app.screen, sc.ScaryDeleteConfirmModal)
            app.screen.query_one("#btn-scarydel-yes").action_press()
            await pilot.pause(0.2)
            assert sc._load_collections() == []
            app.exit()


class TestSaveLoadedPlasmidToCollection:
    """The '+' button in plasmids view writes the current record into the
    active collection (via _save_library → _sync_active_collection_plasmids)."""

    @pytest.fixture(autouse=True)
    def _no_seed(self, monkeypatch):
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

    async def test_add_record_lands_in_active_collection(self, tiny_record):
        sc._save_collections([{"name": "MyCol", "plasmids": []}])
        sc._set_active_collection_name("MyCol")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.2)  # let preload + add_entry settle
            coll = next(c for c in sc._load_collections()
                        if c.get("name") == "MyCol")
            ids = [p.get("id") for p in coll["plasmids"]]
            assert tiny_record.id in ids
            app.exit()
