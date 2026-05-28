"""
test_universal_save — the universal "name + collection" save flow.

Covers the shared pieces every library-save path now routes through:
  * `_name_modal_result` — normalises a NamePlasmidModal payload
    (dict OR bare str OR None) into (name, collection) | None.
  * `_commit_library_entry_to_collection` — appends an entry to a chosen
    collection, collision-renaming name (" COPY") + id ("_N"), mirroring
    the active library only when the target IS active.
  * `NamePlasmidModal` collection mode — mounts the picker + dismisses
    with {"name", "collection"}.

The per-pane wiring (Constructor / Gibson / Traditional) is exercised by
their own suites; here we lock the shared contract + a non-active route.
"""
from __future__ import annotations

import splicecraft as sc

TERMINAL_SIZE = (160, 48)


# ═══════════════════════════════════════════════════════════════════════
# _name_modal_result
# ═══════════════════════════════════════════════════════════════════════

class TestNameModalResult:
    def test_dict_with_collection(self):
        assert sc._name_modal_result(
            {"name": "My Plasmid", "collection": "Bin A"}, "Default"
        ) == ("My Plasmid", "Bin A")

    def test_dict_without_collection_uses_default(self):
        assert sc._name_modal_result(
            {"name": "My Plasmid"}, "Default"
        ) == ("My Plasmid", "Default")

    def test_bare_string_uses_default_collection(self):
        assert sc._name_modal_result("Legacy Name", "Active") == (
            "Legacy Name", "Active")

    def test_blank_and_none_are_cancel(self):
        assert sc._name_modal_result(None, "Default") is None
        assert sc._name_modal_result("", "Default") is None
        assert sc._name_modal_result("   ", "Default") is None
        assert sc._name_modal_result({"name": "  "}, "Default") is None
        assert sc._name_modal_result({}, "Default") is None


# ═══════════════════════════════════════════════════════════════════════
# _commit_library_entry_to_collection
# ═══════════════════════════════════════════════════════════════════════

def _entry(name, *, eid=None, seq="ACGTACGTACGT"):
    return {
        "id": eid or name.replace(" ", "_"),
        "name": name,
        "size": len(seq),
        "n_feats": 0,
        "source": "test",
        "added": "2026-05-28",
        "gb_text": (f"LOCUS       {(eid or 'x')[:16]:<16} {len(seq)} bp "
                    "DNA     linear\nORIGIN\n//\n"),
    }


class TestCommitToCollection:
    def test_creates_collection_and_lands_entry(self, isolated_library):
        final = sc._commit_library_entry_to_collection(
            _entry("Widget"), "Fresh Collection")
        assert final == "Widget"
        colls = sc._load_collections()
        target = next(c for c in colls if c.get("name") == "Fresh Collection")
        names = [e.get("name") for e in target["plasmids"]]
        assert "Widget" in names

    def test_name_collision_gets_copy_suffix(self, isolated_library):
        sc._commit_library_entry_to_collection(_entry("Dup"), "C1")
        final = sc._commit_library_entry_to_collection(_entry("Dup"), "C1")
        assert final == "Dup COPY"
        final3 = sc._commit_library_entry_to_collection(_entry("Dup"), "C1")
        assert final3 == "Dup COPY 2"

    def test_collision_suffix_uses_spaces_not_underscores(self,
                                                          isolated_library):
        sc._commit_library_entry_to_collection(_entry("My Plasmid"), "C2")
        final = sc._commit_library_entry_to_collection(
            _entry("My Plasmid"), "C2")
        assert "_" not in final  # display name stays underscore-free

    def test_id_collision_gets_underscore_suffix(self, isolated_library):
        sc._commit_library_entry_to_collection(
            _entry("A", eid="dup_id"), "C3")
        sc._commit_library_entry_to_collection(
            _entry("B", eid="dup_id"), "C3")
        colls = sc._load_collections()
        target = next(c for c in colls if c.get("name") == "C3")
        ids = [e.get("id") for e in target["plasmids"]]
        assert "dup_id" in ids and "dup_id_2" in ids

    def test_non_active_collection_not_mirrored_to_library(self,
                                                           isolated_library):
        before = {e.get("id") for e in sc._load_library()}
        sc._commit_library_entry_to_collection(
            _entry("OffStage", eid="offstage"), "Some Other Collection")
        after = {e.get("id") for e in sc._load_library()}
        # The active library mirror is untouched for a non-active target.
        assert "offstage" not in (after - before)


# ═══════════════════════════════════════════════════════════════════════
# NamePlasmidModal — collection mode
# ═══════════════════════════════════════════════════════════════════════

class TestNamePlasmidCollectionMode:
    async def test_legacy_mode_has_no_collection_picker(self, tiny_record,
                                                        isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("thing"))
            await pilot.pause()
            modal = app.screen
            from textual.css.query import NoMatches
            try:
                modal.query_one("#nameplasmid-collection")
                found = True
            except NoMatches:
                found = False
            assert found is False

    async def test_collection_mode_mounts_picker_and_returns_dict(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            captured = {}

            def _cb(result):
                captured["result"] = result

            app.push_screen(
                sc.NamePlasmidModal("widget", default_collection="Default"),
                callback=_cb,
            )
            await pilot.pause()
            modal = app.screen
            modal.query_one("#nameplasmid-collection", sc.Select)
            modal.query_one("#nameplasmid-input", sc.Input).value = "Widget 2"
            await pilot.pause()
            modal._try_submit()
            await pilot.pause()
            assert isinstance(captured["result"], dict)
            assert captured["result"]["name"] == "Widget 2"
            assert isinstance(captured["result"]["collection"], str)
