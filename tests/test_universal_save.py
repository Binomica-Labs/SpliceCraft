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

    # ── replace_id (document re-save, 2026-05-30) ──────────────────────
    def test_replace_id_replaces_in_place(self, isolated_library):
        """A re-save with `replace_id` overwrites the existing slot —
        no ` COPY`, id preserved, count unchanged."""
        sc._commit_library_entry_to_collection(
            _entry("Doc v1", eid="doc"), "Edit Coll")
        final = sc._commit_library_entry_to_collection(
            _entry("Doc v2", eid="doc", seq="TTTTGGGGCCCC"), "Edit Coll",
            replace_id="doc")
        target = next(c for c in sc._load_collections()
                      if c.get("name") == "Edit Coll")
        docs = [e for e in target["plasmids"] if e.get("id") == "doc"]
        assert final == "Doc v2"
        assert len(docs) == 1 and docs[0]["name"] == "Doc v2"
        assert len(target["plasmids"]) == 1

    def test_replace_id_not_found_appends(self, isolated_library):
        sc._commit_library_entry_to_collection(_entry("X", eid="x"), "C9")
        sc._commit_library_entry_to_collection(
            _entry("Y", eid="y"), "C9", replace_id="ghost")
        target = next(c for c in sc._load_collections()
                      if c.get("name") == "C9")
        assert {e.get("id") for e in target["plasmids"]} == {"x", "y"}

    def test_replace_id_keeps_own_name_no_copy(self, isolated_library):
        """Replacing must not ` COPY`-rename the entry against itself."""
        sc._commit_library_entry_to_collection(_entry("Keeper", eid="k"), "C10")
        sc._commit_library_entry_to_collection(_entry("Other", eid="o"), "C10")
        final = sc._commit_library_entry_to_collection(
            _entry("Keeper", eid="k", seq="AAAACCCCGGGG"), "C10",
            replace_id="k")
        assert final == "Keeper"


# ═══════════════════════════════════════════════════════════════════════
# Library-entry kind classification (Kind column, 2026-05-30)
# ═══════════════════════════════════════════════════════════════════════

class TestEntryKind:
    _CIRC = "LOCUS       x   100 bp ds-DNA circular SYN 01-JAN-2026\n//\n"
    _LIN = "LOCUS       x   100 bp ds-DNA linear SYN 01-JAN-2026\n//\n"

    def test_derive_plasmid_from_circular(self):
        assert sc._derive_entry_kind(
            {"source": "constructor:gb", "gb_text": self._CIRC}) == "plasmid"

    def test_derive_fragment_from_linear(self):
        assert sc._derive_entry_kind(
            {"source": "file:x.gb", "gb_text": self._LIN}) == "fragment"

    def test_derive_protein_from_source(self):
        # Protein source beats topology (a protein CDS is linear DNA).
        assert sc._derive_entry_kind(
            {"source": "synthesis-protein:p1", "gb_text": self._LIN}) == "protein"

    def test_derive_amplicon_from_source(self):
        assert sc._derive_entry_kind(
            {"source": "simulator:pcr", "gb_text": self._LIN}) == "amplicon"

    def test_derive_topology_field_beats_gb_text(self):
        assert sc._derive_entry_kind(
            {"source": "x", "topology": "circular"}) == "plasmid"

    def test_entry_kind_uses_explicit_field(self):
        # Explicit kind overrides what derive would infer.
        assert sc._entry_kind(
            {"kind": "amplicon", "source": "constructor:gb",
             "gb_text": self._CIRC}) == "amplicon"

    def test_entry_kind_invalid_field_falls_back_to_derive(self):
        assert sc._entry_kind(
            {"kind": "bogus", "gb_text": self._CIRC}) == "plasmid"

    def test_entry_kind_missing_field_derives(self):
        assert sc._entry_kind({"gb_text": self._LIN}) == "fragment"


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


# ═══════════════════════════════════════════════════════════════════════
# add_entry(target_collection=…) routing — the Synthesis save-to-
# collection integration point (2026-05-30)
# ═══════════════════════════════════════════════════════════════════════

class TestAddEntryCollectionRouting:
    async def test_routes_to_non_active_collection_and_replaces(
            self, tiny_record, isolated_library):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            lib = app.query_one("#library", sc.LibraryPanel)
            active = sc._get_active_collection_name() or ""
            rec = SeqRecord(Seq("ACGT" * 20), id="frag_rt", name="frag_rt")
            rec.annotations["topology"] = "linear"
            rec.annotations["molecule_type"] = "DNA"
            rec._tui_kind = "fragment"   # synthesis stamps this
            assert lib.add_entry(rec, target_collection="Side Bin") is True
            await pilot.pause()
            side = next(c for c in sc._load_collections()
                        if c.get("name") == "Side Bin")
            landed = [e for e in side["plasmids"] if e.get("id") == "frag_rt"]
            assert len(landed) == 1
            assert landed[0].get("kind") == "fragment"   # #3 stamp carried
            if active != "Side Bin":
                # Non-active target → active library mirror untouched.
                assert not any(e.get("id") == "frag_rt"
                               for e in sc._load_library())
            # Re-save in place (document model) → no duplicate.
            rec2 = SeqRecord(Seq("TTTT" * 20), id="frag_rt", name="frag_rt")
            rec2.annotations["topology"] = "linear"
            rec2.annotations["molecule_type"] = "DNA"
            assert lib.add_entry(rec2, target_collection="Side Bin",
                                 replace_in_target=True) is True
            await pilot.pause()
            side = next(c for c in sc._load_collections()
                        if c.get("name") == "Side Bin")
            assert sum(1 for e in side["plasmids"]
                       if e.get("id") == "frag_rt") == 1
