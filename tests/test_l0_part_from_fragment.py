"""L0 part from a synthesised fragment.

Turning a saved synthetic fragment (FRAG) into a Level-0 part by cloning it
into the grammar's entry vector — the headless core shared by the Constructor
modal, the agent endpoint and BABS.

The interesting property is the TWO-TIER layout ([INV-136] nesting): a nested
fragment carries entry overhangs wrapping the category pair, so anything that
recovers overhangs by digesting reads back the WRONG pair. These tests pin
that the part's overhangs come from the grammar's position table and that a
mismatched part type is refused rather than stored.
"""
from __future__ import annotations

import pytest

import splicecraft as sc
import splicecraft_cloning as cl


BODY = "ATG" + "GCTAGCTAGCTAGCATCGATCGGATCC" * 3 + "TAA"


@pytest.fixture
def grammar():
    return dict(sc._BUILTIN_GRAMMARS["gb_l0"])


def _gb(seq: str, name: str, topology: str = "linear") -> str:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    r = SeqRecord(Seq(seq), id=name[:16], name=name[:16])
    r.annotations["molecule_type"] = "DNA"
    r.annotations["topology"] = topology
    return sc._record_to_gb_text(r)


def _position(grammar: dict, want: str) -> dict:
    return next(p for p in grammar["positions"]
                if (p.get("type") or "").upper() == want.upper())


def _entry_vector(grammar: dict, e5: str, e3: str) -> str:
    """A pUPD2-shaped acceptor that releases `e5`/`e3` when cut."""
    pad = grammar.get("pad")
    site = grammar.get("site")
    spacer = grammar.get("spacer")
    return ("GGCGCGTTAACCGGTTAACCGG" * 4 + pad + site + spacer + e5
            + "TTTTGGGGCCCCAAAA" * 2 + e3 + spacer + sc._rc(site) + pad
            + "ACGTACGTTGCA" * 4)


def _seed_frag_and_vector(grammar: dict, *, nested: bool = False,
                          part_type: str = "CDS", entry_id: str = "frag1",
                          name: str = "MyFrag"):
    """A FRAG in the library + the matching entry vector configured for
    gb_l0 — the state a user is in after Synthesis → L0 Fragment → Save."""
    built = _fragment(grammar, part_type, nested=nested)
    vec_seq = _entry_vector(grammar, built["entry_oh5"], built["entry_oh3"])
    lib = sc._load_library()
    lib.append({"id": entry_id, "name": name, "kind": "fragment",
                "size": len(built["fragment"]),
                "gb_text": _gb(built["fragment"], name)})
    sc._save_library(lib)
    sc._set_entry_vector("gb_l0", {"name": "pUPD2-test",
                                   "gb_text": _gb(vec_seq, "pUPD2",
                                                  "circular")})
    return built


def _fragment(grammar: dict, part_type: str, *, nested: bool):
    pos = _position(grammar, part_type)
    ext = ("CTCG", "TGAG") if nested else None
    return cl._build_synthesis_l0_fragment(
        BODY, pos["oh5"], pos["oh3"], grammar=grammar,
        part_type=part_type, entry_overhangs=ext)


class TestReleasedInsert:
    def test_plain_fragment_is_refused_by_name(self, grammar):
        with pytest.raises(ValueError, match="plain fragment"):
            cl._released_insert_from_fragment("ATGC" * 20,
                                              grammar["enzyme"])

    def test_one_sided_fragment_is_refused(self, grammar):
        half = (grammar["pad"] + grammar["site"] + grammar["spacer"]
                + "AATG" + BODY)
        with pytest.raises(ValueError, match="only one"):
            cl._released_insert_from_fragment(half, grammar["enzyme"])

    def test_missing_enzyme_is_refused(self):
        with pytest.raises(ValueError, match="no Type IIS enzyme"):
            cl._released_insert_from_fragment("ATGC" * 20, "")

    def test_enzyme_not_in_the_catalog_blames_the_grammar(self, grammar):
        """A custom grammar naming an unknown enzyme finds no sites — which
        without this reads as "your fragment isn't wrapped", sending the user
        to rebuild a fragment that was fine."""
        frag = _fragment(grammar, "CDS", nested=False)["fragment"]
        with pytest.raises(ValueError, match="isn't in the enzyme catalog"):
            cl._released_insert_from_fragment(frag, "NotAnEnzyme")


class TestPartFromFragment:
    @pytest.mark.parametrize("nested", [False, True])
    def test_overhangs_come_from_the_grammar_not_the_digest(self, grammar,
                                                            nested):
        # THE two-tier property: a nested fragment's DIGEST yields the entry
        # overhangs (CTCG/TGAG), but the part must carry the CATEGORY pair.
        built = _fragment(grammar, "CDS", nested=nested)
        vec = _entry_vector(grammar, built["entry_oh5"], built["entry_oh3"])
        part = cl._l0_part_from_syn_fragment(
            built["fragment"], vec, grammar=grammar, part_type="CDS",
            name="MyPart")
        cds = _position(grammar, "CDS")
        assert (part["oh5"], part["oh3"]) == (cds["oh5"], cds["oh3"])
        assert part["level"] == 0
        assert part["grammar"] == grammar["id"]
        assert part["nested"] is nested

    @pytest.mark.parametrize("nested", [False, True])
    def test_stored_sequence_is_the_body_between_the_overhangs(self, grammar,
                                                                nested):
        """The L0 storage contract: `sequence` carries NEITHER overhang —
        `oh5`/`oh3` are held apart and the assembly chain re-adds them
        (`_assembly_fragment_strip_oh5`). Storing the flanked region instead
        duplicates 8 bp per part (16 for a two-tier fragment, whose entry pair
        would ride along too) — a frameshift in a CDS, invisible in the row."""
        built = _fragment(grammar, "CDS", nested=nested)
        vec = _entry_vector(grammar, built["entry_oh5"], built["entry_oh3"])
        part = cl._l0_part_from_syn_fragment(
            built["fragment"], vec, grammar=grammar, part_type="CDS")
        seq = part["sequence"]
        assert not seq.startswith(part["oh5"])
        assert not seq.endswith(part["oh3"])
        assert "CTCG" not in seq and "TGAG" not in seq   # entry pair, if any
        # BODY minus the ATG the AATG overhang absorbed — the same body either
        # tier, since how a part ENTERED the vector isn't part of the part.
        assert seq == BODY[3:]

    def test_both_tiers_file_the_identical_part(self, grammar):
        made = []
        for nested in (False, True):
            built = _fragment(grammar, "CDS", nested=nested)
            vec = _entry_vector(grammar, built["entry_oh5"],
                                built["entry_oh3"])
            made.append(cl._l0_part_from_syn_fragment(
                built["fragment"], vec, grammar=grammar, part_type="CDS"))
        flat, nest = made
        assert flat["sequence"] == nest["sequence"]
        assert (flat["oh5"], flat["oh3"]) == (nest["oh5"], nest["oh3"])
        assert (flat["nested"], nest["nested"]) == (False, True)

    @pytest.mark.parametrize("nested", [False, True])
    def test_the_filed_part_reclones_into_the_same_plasmid(self, grammar,
                                                            nested):
        """The storage contract, proved end to end: hand the part we FILE to
        the pre-existing `_clone_part_into_entry_vector` machinery and it must
        rebuild the very circle this feature simulated. When `sequence` still
        carried its overhangs inline, the round trip came back 8 bp longer
        (16 nested) — the duplication that would have shipped."""
        built = _fragment(grammar, "CDS", nested=nested)
        vec_seq = _entry_vector(grammar, built["entry_oh5"],
                                built["entry_oh3"])
        part = cl._l0_part_from_syn_fragment(
            built["fragment"], vec_seq, grammar=grammar, part_type="CDS",
            name="P")
        rec = sc._clone_part_into_entry_vector(
            {"name": "P", "type": part["type"], "oh5": part["oh5"],
             "oh3": part["oh3"], "sequence": part["sequence"]},
            {"name": "EV", "gb_text": _gb(vec_seq, "EV", "circular")},
            grammar)
        assert rec is not None, "the filed part no longer clones at all"
        rt = str(rec.seq).upper()
        sim = part["cloned_seq"].upper()
        # Same circle, allowing for where the record starts.
        assert len(rt) == len(sim) and rt in (sim + sim)

    def test_the_product_is_a_real_cloned_plasmid(self, grammar):
        built = _fragment(grammar, "CDS", nested=False)
        vec = _entry_vector(grammar, built["entry_oh5"], built["entry_oh3"])
        part = cl._l0_part_from_syn_fragment(
            built["fragment"], vec, grammar=grammar, part_type="CDS")
        # Simulated, not hand-built ([INV-127]): the product carries backbone
        # AND insert, so it is longer than either alone.
        assert len(part["cloned_seq"]) > len(built["fragment"])
        assert BODY in part["cloned_seq"] + part["cloned_seq"]

    def test_part_type_absent_from_the_grammar_lists_what_it_has(self,
                                                                 grammar):
        built = _fragment(grammar, "CDS", nested=False)
        vec = _entry_vector(grammar, built["entry_oh5"], built["entry_oh3"])
        with pytest.raises(ValueError) as ei:
            cl._l0_part_from_syn_fragment(
                built["fragment"], vec, grammar=grammar,
                part_type="not_a_position")
        msg = str(ei.value)
        assert "no 'not_a_position' position" in msg
        assert "Promoter" in msg           # names the real options

    @pytest.mark.parametrize("wrong", ["Promoter", "Terminator"])
    def test_wrong_part_type_for_this_fragment_is_refused(self, grammar,
                                                          wrong):
        # A CDS fragment filed as a promoter would store GGAG/AATG against a
        # body that actually carries AATG/GCTT — a part that cannot assemble.
        built = _fragment(grammar, "CDS", nested=False)
        vec = _entry_vector(grammar, built["entry_oh5"], built["entry_oh3"])
        with pytest.raises(ValueError, match="built as a different part type"):
            cl._l0_part_from_syn_fragment(
                built["fragment"], vec, grammar=grammar, part_type=wrong)

    def test_wrong_entry_vector_is_refused(self, grammar):
        built = _fragment(grammar, "CDS", nested=False)
        # A vector that releases overhangs this insert doesn't present.
        vec = _entry_vector(grammar, "TTTT", "AAAA")
        with pytest.raises(ValueError, match="match nothing the vector"):
            cl._l0_part_from_syn_fragment(
                built["fragment"], vec, grammar=grammar, part_type="CDS")

    def test_overhang_check_is_positional_not_a_substring_search(self,
                                                                 grammar):
        """A body that merely CONTAINS the other type's overhang near its
        start must not satisfy that type — the ends are what matter."""
        cds = _position(grammar, "CDS")
        promoter = _position(grammar, "Promoter")
        # A CDS fragment whose body opens with the promoter's 5' overhang.
        body = "ATG" + promoter["oh5"] + "GGCATCGATCGGATCCTAA" * 3
        built = cl._build_synthesis_l0_fragment(
            body, cds["oh5"], cds["oh3"], grammar=grammar, part_type="CDS")
        vec = _entry_vector(grammar, built["entry_oh5"], built["entry_oh3"])
        with pytest.raises(ValueError, match="built as a different part type"):
            cl._l0_part_from_syn_fragment(
                built["fragment"], vec, grammar=grammar, part_type="Promoter")

    def test_vector_leaving_one_overhang_at_both_cuts_is_refused(self,
                                                                  grammar):
        """A vector whose two cuts leave the SAME overhang has no piece left
        over as a backbone — and an insert with matching ends would happily
        close on ITSELF, reporting a backbone-free circle as a clean clone."""
        frag = cl._build_synthesis_l0_fragment(
            BODY, "AATG", "AATG", grammar=grammar, part_type="CDS",
            entry_overhangs=None)["fragment"]
        with pytest.raises(ValueError, match="backbone"):
            cl._clone_syn_fragment_into_entry_vector(
                frag, _entry_vector(grammar, "AATG", "AATG"), grammar=grammar)

    def test_a_fragment_with_no_insert_is_refused(self, grammar):
        """Type IIS ends back to back: the overhangs are both there, so the
        positional guard passes, and the part would be a pair of overhangs
        around nothing."""
        cds = _position(grammar, "CDS")
        empty = (grammar["pad"] + grammar["site"] + grammar["spacer"]
                 + cds["oh5"] + cds["oh3"] + grammar["spacer"]
                 + sc._rc(grammar["site"]) + grammar["pad"])
        vec = _entry_vector(grammar, cds["oh5"], cds["oh3"])
        with pytest.raises(ValueError, match="carries no insert"):
            cl._l0_part_from_syn_fragment(empty, vec, grammar=grammar,
                                          part_type="CDS")

    def test_iupac_bases_survive_the_trip(self, grammar):
        """Ns and degenerate codes are legal in an ordered fragment; they must
        not be silently dropped out of the part's body."""
        cds = _position(grammar, "CDS")
        body = "ATGNNNRYSWGCATCGATCGGATCCTAA" * 2
        built = cl._build_synthesis_l0_fragment(
            body, cds["oh5"], cds["oh3"], grammar=grammar, part_type="CDS")
        vec = _entry_vector(grammar, built["entry_oh5"], built["entry_oh3"])
        part = cl._l0_part_from_syn_fragment(built["fragment"], vec,
                                             grammar=grammar,
                                             part_type="CDS")
        assert part["sequence"] == body[3:]

    def test_empty_vector_is_refused(self, grammar):
        built = _fragment(grammar, "CDS", nested=False)
        with pytest.raises(ValueError, match="no sequence"):
            cl._l0_part_from_syn_fragment(
                built["fragment"], "", grammar=grammar, part_type="CDS")


class TestAgentEndpoint:
    """`make-l0-part-from-fragment` — the surface agents and BABS drive."""

    def _handler(self):
        return sc._state._AGENT_HANDLERS["make-l0-part-from-fragment"][0]

    def _seed(self, grammar, *, nested=False, part_type="CDS"):
        return _seed_frag_and_vector(grammar, nested=nested,
                                     part_type=part_type)

    def test_registered_as_a_write_endpoint(self):
        H = sc._state._AGENT_HANDLERS
        assert H["make-l0-part-from-fragment"][1] is True

    def test_builds_and_files_the_part(self, grammar):
        self._seed(grammar)
        res = self._handler()(None, {"fragment": "frag1", "part_type": "CDS"})
        assert not isinstance(res, tuple), res
        cds = _position(grammar, "CDS")
        assert (res["oh5"], res["oh3"]) == (cds["oh5"], cds["oh3"])
        assert res["plasmid_len"] > 0 and res["entry_vector"] == "pUPD2-test"
        filed = [p for p in sc._load_parts_bin() if p.get("name") == "MyFrag"]
        assert len(filed) == 1 and filed[0]["level"] == 0

    def test_nested_fragment_files_the_category_overhangs(self, grammar):
        # The two-tier case, end to end through the endpoint.
        self._seed(grammar, nested=True)
        res = self._handler()(None, {"fragment": "frag1", "part_type": "CDS"})
        assert not isinstance(res, tuple), res
        cds = _position(grammar, "CDS")
        assert (res["oh5"], res["oh3"]) == (cds["oh5"], cds["oh3"])

    def test_unknown_fragment_is_404(self, grammar):
        payload = {"fragment": "nope", "part_type": "CDS"}
        assert self._handler()(None, payload)[1] == 404

    def test_missing_part_type_is_400(self, grammar):
        self._seed(grammar)
        assert self._handler()(None, {"fragment": "frag1"})[1] == 400

    def test_grammar_without_an_entry_vector_is_refused(self, grammar):
        self._seed(grammar)
        sc._set_entry_vector("gb_l0", None)
        payload, status = self._handler()(None, {"fragment": "frag1",
                                                 "part_type": "CDS"})
        assert status == 400 and "entry vector" in payload["error"]

    def test_wrong_part_type_is_refused_with_the_reason(self, grammar):
        self._seed(grammar)
        payload, status = self._handler()(None, {"fragment": "frag1",
                                                 "part_type": "Promoter"})
        assert status == 400
        assert "different part type" in payload["error"]

    def test_the_filed_sequence_is_the_body_not_the_flanked_region(self,
                                                                    grammar):
        """What lands in the bin obeys the same contract as `create-part`."""
        self._seed(grammar, nested=True)
        res = self._handler()(None, {"fragment": "frag1", "part_type": "CDS"})
        assert not isinstance(res, tuple), res
        filed = next(p for p in sc._load_parts_bin()
                     if p.get("name") == "MyFrag")
        assert filed["sequence"] == BODY[3:]
        assert res["part_len"] == len(BODY) - 3 and res["nested"] is True

    @pytest.mark.parametrize("nested", [False, True])
    def test_the_l0_plasmid_is_saved_to_the_library(self, grammar, nested):
        """Both artifacts, or the feature misses its point: the part goes to
        the bin AND the circular L0 plasmid to the library."""
        self._seed(grammar, nested=nested)
        before = {e.get("id") for e in sc._load_library()}
        res = self._handler()(None, {"fragment": "frag1", "part_type": "CDS"})
        assert not isinstance(res, tuple), res
        # Named apart from the FRAG it was built from — which is sitting in
        # the library as "MyFrag" — so neither gets a bare dedup suffix.
        assert res["saved_id"] and res["saved_name"] == "MyFrag (L0)"
        new = [e for e in sc._load_library() if e.get("id") not in before]
        assert len(new) == 1, "expected exactly one new library entry"
        rec = sc._gb_text_to_record(new[0]["gb_text"])
        assert (rec.annotations.get("topology") or "").lower() == "circular"
        # It IS the simulated construct — the part's body rode in with it.
        seq = str(rec.seq).upper()
        assert BODY[3:] in seq + seq
        assert res["plasmid_len"] == len(seq) == new[0]["size"]

    def test_the_saved_plasmid_carries_its_lineage(self, grammar):
        """History must show where it came from — the insert and the entry
        vector — the same tree `domesticate-part` builds."""
        self._seed(grammar)
        self._handler()(None, {"fragment": "frag1", "part_type": "CDS"})
        entry = next(e for e in sc._load_library()
                     if e.get("name") == "MyFrag (L0)")
        assert entry.get("history_xml"), "no construction lineage attached"

    def test_save_plasmid_false_files_only_the_part(self, grammar):
        self._seed(grammar)
        before = len(sc._load_library())
        res = self._handler()(None, {"fragment": "frag1", "part_type": "CDS",
                                     "save_plasmid": False})
        assert not isinstance(res, tuple), res
        assert res["saved_id"] is None
        assert len(sc._load_library()) == before
        assert any(p.get("name") == "MyFrag" for p in sc._load_parts_bin())

    def test_unknown_collection_fails_before_the_part_is_filed(self, grammar):
        """Ordering guarantee: a bad destination must not leave a part filed
        with no plasmid beside it."""
        self._seed(grammar)
        payload, status = self._handler()(None, {
            "fragment": "frag1", "part_type": "CDS", "collection": "nope"})
        assert status == 404 and "collection" in payload["error"]
        assert not any(p.get("name") == "MyFrag" for p in sc._load_parts_bin())

    def test_a_circular_entry_is_refused_by_name(self, grammar):
        """Digesting a plasmid as if it were linear isn't an error — it
        "releases" whatever sits between two sites, so aiming this at pUPD2
        would file its stuffer as a part."""
        self._seed(grammar)
        lib = sc._load_library()
        lib.append({"id": "circ", "name": "SomePlasmid", "kind": "plasmid",
                    "size": 40,
                    "gb_text": _gb("ATGC" * 40, "SomePlasmid", "circular")})
        sc._save_library(lib)
        payload, status = self._handler()(None, {"fragment": "circ",
                                                 "part_type": "CDS"})
        assert status == 400 and "circular plasmid" in payload["error"]
        assert not sc._load_parts_bin()

    def test_a_mistyped_routing_param_is_refused_not_ignored(self, grammar):
        """`parts_bin` isn't this endpoint's spelling of `bin` — routing a
        part somewhere else while answering ok is the SC-D footgun."""
        self._seed(grammar)
        payload, status = self._handler()(None, {
            "fragment": "frag1", "part_type": "CDS", "parts_bin": "Other"})
        assert status == 400 and "routing" in payload["error"]
        assert not sc._load_parts_bin()

    def test_product_name_names_the_plasmid_only(self, grammar):
        self._seed(grammar)
        res = self._handler()(None, {"fragment": "frag1", "part_type": "CDS",
                                     "name": "ThePart",
                                     "product_name": "ThePlasmid"})
        assert not isinstance(res, tuple), res
        assert res["saved_name"] == "ThePlasmid" and res["name"] == "ThePart"
        assert any(p.get("name") == "ThePart" for p in sc._load_parts_bin())
        assert any(e.get("name") == "ThePlasmid" for e in sc._load_library())

    def test_log_label_says_ui_when_the_button_is_the_caller(self, grammar,
                                                             monkeypatch):
        """The Constructor routes through this handler, so a button press
        must not be logged as a script's doing. (Moved handlers resolve
        `_log_event` in the SIBLING's namespace — patch there, not on `sc`.)"""
        import splicecraft_agent as ag
        seen = []
        monkeypatch.setattr(ag, "_log_event",
                            lambda ev, **kw: seen.append((ev, kw)))
        self._seed(grammar)
        self._handler()(None, {"fragment": "frag1", "part_type": "CDS",
                               "via": "ui"})
        assert seen and seen[-1][1]["via"] == "ui"
        seen.clear()
        self._handler()(None, {"fragment": "frag1", "part_type": "CDS",
                               "name": "Second"})
        assert seen and seen[-1][1]["via"] == "agent"

    def test_duplicate_part_is_409(self, grammar):
        self._seed(grammar)
        assert not isinstance(
            self._handler()(None, {"fragment": "frag1", "part_type": "CDS"}),
            tuple)
        again = self._handler()(None, {"fragment": "frag1",
                                       "part_type": "CDS"})
        assert isinstance(again, tuple) and again[1] == 409


class TestPickerModal:
    """`PartFromSynFragModal` — the Constructor's picker."""

    def _seed_library(self, grammar):
        built = _seed_frag_and_vector(grammar, entry_id="f1", name="AFrag")
        lib = sc._load_library()
        lib.append({"id": "p1", "name": "APlasmid", "kind": "plasmid",
                    "size": 100, "gb_text": _gb("ATGC" * 40, "APlasmid")})
        sc._save_library(lib)
        return built

    async def test_lists_only_fragments_not_plasmids(self, grammar):
        from textual.widgets import DataTable
        self._seed_library(grammar)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PartFromSynFragModal(grammar))
            await pilot.pause()
            modal = app.screen
            names = [str(e.get("name")) for e in modal._frags]
            assert names == ["AFrag"], "a non-fragment entry was offered"
            assert modal.query_one("#pfsf-table", DataTable).row_count == 1
            app.exit()

    async def test_part_types_come_from_the_grammar(self, grammar):
        from textual.widgets import Select
        self._seed_library(grammar)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PartFromSynFragModal(grammar))
            await pilot.pause()
            sel = app.screen.query_one("#pfsf-type", Select)
            # `_options` carries the prompt's Select.NULL sentinel, which is
            # truthy — filter on type, not truthiness.
            offered = {v for _lbl, v in sel._options
                       if isinstance(v, str) and v}
            grammar_types = {str(p.get("type")) for p in grammar["positions"]
                             if p.get("oh5") and p.get("oh3")}
            # Only this grammar's own positions — that's the guard against
            # pairing a part type with a grammar that doesn't define it.
            assert offered and offered <= grammar_types
            app.exit()

    async def test_missing_entry_vector_is_flagged_up_front(self, grammar):
        from textual.widgets import Static
        self._seed_library(grammar)
        sc._set_entry_vector("gb_l0", None)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PartFromSynFragModal(grammar))
            await pilot.pause()
            txt = str(app.screen.query_one("#pfsf-ev", Static).render())
            assert "No entry vector" in txt
            app.exit()

    async def test_constructor_button_opens_the_picker(self, grammar):
        """The button is routed by ConstructorModal's SINGLE `_on_button`
        dispatcher; a second catch-all `@on(Button.Pressed)` would
        double-dispatch every button in the modal, so the wiring is worth
        pinning."""
        from textual.widgets import Button
        self._seed_library(grammar)
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 55)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            ctor = app.screen
            btn = ctor.query_one("#btn-ctor-from-frag-gb_l0", Button)
            btn.press()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.PartFromSynFragModal)
            app.exit()

    async def test_the_whole_button_flow_files_part_and_saves_plasmid(
            self, grammar):
        """End to end from the button: pick the fragment, choose the type,
        press Create Part — BOTH artifacts must land.

        This one runs the REAL threaded worker, which is the point: the
        endpoint's library save finishes by refreshing the panel through
        `app.call_from_thread`, and that RAISES if it is already on the app's
        own thread. Calling the handler straight from the modal callback (the
        obvious way to write this) blows up here, after the part is filed."""
        from textual.widgets import Button, DataTable, Select
        self._seed_library(grammar)
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 55)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            app.screen.query_one("#btn-ctor-from-frag-gb_l0", Button).press()
            await pilot.pause()
            picker = app.screen
            assert isinstance(picker, sc.PartFromSynFragModal)
            picker.query_one("#pfsf-table", DataTable).move_cursor(row=0)
            picker.query_one("#pfsf-type", Select).value = "CDS"
            picker.query_one("#btn-pfsf-go", Button).press()
            await pilot.pause()
            await app.workers.wait_for_complete()
            await pilot.pause()
            parts = [p for p in sc._load_parts_bin()
                     if p.get("name") == "AFrag"]
            assert len(parts) == 1, "the part was not filed"
            assert parts[0]["sequence"] == BODY[3:]
            saved = [e for e in sc._load_library()
                     if e.get("name") == "AFrag (L0)"]
            assert len(saved) == 1, "the L0 plasmid was not saved"
            app.exit()

    async def test_nothing_picked_is_refused_not_crashed(self, grammar):
        from textual.widgets import Static
        self._seed_library(grammar)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.PartFromSynFragModal(grammar))
            await pilot.pause()
            modal = app.screen
            modal._go(None)                      # no part type selected
            await pilot.pause()
            txt = str(modal.query_one("#pfsf-status", Static).render())
            assert "part type" in txt.lower()
            app.exit()


class TestGrammarPositionLookup:
    def test_case_tolerant_but_delegates_to_the_exact_matcher(self, grammar):
        assert cl._grammar_position_for_type(grammar, "CDS") is not None
        assert cl._grammar_position_for_type(grammar, "cds") is not None
        # Same row the pre-existing exact matcher finds — one definition of
        # what a position is, not two.
        assert (cl._grammar_position_for_type(grammar, "CDS")
                is cl._grammar_position_by_type(grammar, "CDS"))

    def test_position_names_are_not_matched(self, grammar):
        # `Pos 3-4` is the CDS row's NAME; matching it could alias another
        # row's type.
        assert cl._grammar_position_for_type(grammar, "Pos 3-4") is None

    def test_unknown_and_blank_return_none(self, grammar):
        assert cl._grammar_position_for_type(grammar, "nope") is None
        assert cl._grammar_position_for_type(grammar, "") is None



# A real ORF: ATG, all-Ala codons, one terminator, and no Type IIS site on
# either strand. The module-level `BODY` carries internal `TAG`s (harmless
# for the storage-contract tests above, useless for a frame test).
ORF_BODY = "ATG" + "GCTGCAGCAGCTGCAGCA" * 4 + "TAA"


class TestClonedCdsFrame:
    """The CDS annotation on the saved L0 plasmid must sit ON the ORF.

    2026-07-28 (user report: "why is there a broken frame"): Constructor →
    new part from syn fragment saved a plasmid whose CDS spanned the WHOLE
    released insert. A two-tier acceptor releases `[CTCG][AATG]<body>[GCTT]`,
    so that annotation began on the vector-side external overhang — 5 bases
    before the ATG, out of frame. The map + seq panel translated garbage and
    stamped the premature-stop ⚠ over a perfectly good ORF.

    The DNA was always right; only the annotation moved. Pinned on BOTH tiers
    because they are different branches: the one-tier
    `_clone_part_synthesise_insert` had the mirror-image bug — it started at
    codon 2 and dropped the ATG that the AATG overhang supplies.
    """

    @staticmethod
    def _cloned(grammar, *, nested):
        pos = _position(grammar, "CDS")
        built = cl._build_synthesis_l0_fragment(
            ORF_BODY, pos["oh5"], pos["oh3"], grammar=grammar,
            part_type="CDS",
            entry_overhangs=("CTCG", "TGAG") if nested else None)
        vec_seq = _entry_vector(grammar, built["entry_oh5"],
                                built["entry_oh3"])
        part = cl._l0_part_from_syn_fragment(
            built["fragment"], vec_seq, grammar=grammar,
            part_type="CDS", name="MyPart")
        if nested:
            # Exactly what `domesticate-part` stamps when the configured
            # vector's external pair differs from the category pair.
            part["entry_oh5"] = built["entry_oh5"]
            part["entry_oh3"] = built["entry_oh3"]
        rec = sc._clone_part_into_entry_vector(
            part, {"gb_text": _gb(vec_seq, "pUPD2", "circular")}, grammar)
        assert rec is not None, "the clone simulation bailed"
        cds = [f for f in rec.features
               if f.type == "CDS"
               and (f.qualifiers.get("label") or [""])[0] == "MyPart"]
        assert len(cds) == 1
        return str(rec.seq).upper(), cds[0]

    @pytest.mark.parametrize("nested", [False, True])
    def test_cds_starts_on_the_atg_the_overhang_supplies(self, grammar,
                                                         nested):
        seq, cds = self._cloned(grammar, nested=nested)
        start = int(cds.location.start)
        assert seq[start:start + 3] == "ATG", (
            f"CDS starts on {seq[start:start + 3]!r}, not the start codon "
            f"(context {seq[max(0, start - 6):start + 6]!r})")

    @pytest.mark.parametrize("nested", [False, True])
    def test_cds_is_in_frame_and_has_no_premature_stop(self, grammar,
                                                       nested):
        from Bio.Seq import Seq
        seq, cds = self._cloned(grammar, nested=nested)
        start, end = int(cds.location.start), int(cds.location.end)
        assert (end - start) % 3 == 0, "CDS isn't a whole number of codons"
        aa = str(Seq(seq[start:end]).translate())
        # The rule `PlasmidMap._parse` applies before stamping ⚠: a trailing
        # run of stops collapses to one terminator, anything else is premature.
        assert aa.startswith("M")
        assert aa.rstrip("*").count("*") == 0, f"premature stop in {aa!r}"

    def test_the_external_overhang_stays_outside_the_cds(self, grammar):
        """The regression proper: annotating the released insert from index 0
        swallows the vector-side `CTCG` and shifts the frame by 5."""
        seq, cds = self._cloned(grammar, nested=True)
        start, end = int(cds.location.start), int(cds.location.end)
        assert start != 0
        assert "CTCG" not in seq[start:end]
        # …and the body it DOES cover is the part's own sequence.
        assert ORF_BODY[3:] in seq[start:end]

    @pytest.mark.parametrize("nested", [False, True])
    def test_a_non_coding_part_is_not_extended_into_its_overhang(
            self, grammar, nested):
        """`_atg_offset_for_part` is coding-only: a Promoter's GGAG 5'
        overhang carries no start codon, so its annotation must stop at the
        body — otherwise every non-coding part grows 3 phantom bases."""
        pos = _position(grammar, "Promoter")
        built = cl._build_synthesis_l0_fragment(
            ORF_BODY, pos["oh5"], pos["oh3"], grammar=grammar,
            part_type="Promoter",
            entry_overhangs=("CTCG", "TGAG") if nested else None)
        vec_seq = _entry_vector(grammar, built["entry_oh5"],
                                built["entry_oh3"])
        part = cl._l0_part_from_syn_fragment(
            built["fragment"], vec_seq, grammar=grammar,
            part_type="Promoter", name="MyProm")
        if nested:
            part["entry_oh5"], part["entry_oh3"] = ("CTCG", "TGAG")
        rec = sc._clone_part_into_entry_vector(
            part, {"gb_text": _gb(vec_seq, "pUPD2", "circular")}, grammar)
        assert rec is not None
        f = next(f for f in rec.features
                 if (f.qualifiers.get("label") or [""])[0] == "MyProm")
        seq = str(rec.seq).upper()
        body = seq[int(f.location.start):int(f.location.end)]
        assert body == part["sequence"]
        assert not body.startswith(part["oh5"])


class TestPartFeaturePlacement:
    """Unit-level hardening of the placer behind [INV-173].

    `TestClonedCdsFrame` drives it through the whole clone; these pin the
    layouts and the degenerate inputs a clone won't reach on its own.
    """

    BODY = "GTGTCCAGCGGGGAGGACATCTTCTCGGGCCTCGTGCCC"      # codon 2 on, 39 nt
    TPL = {"start": 0, "end": 0, "strand": 1, "type": "CDS",
           "label": "MyPart", "color": "white"}

    def _place(self, top, insert, oh5="AATG", oh3="GCTT", ptype="CDS"):
        f = sc._clone_part_place_part_feature(
            dict(self.TPL), top, insert, oh5, oh3, ptype)
        return f["start"], f["end"], f.get("strand", 1)

    def test_two_tier_release_skips_the_vector_side_overhang(self):
        top = "CTCG" + "AATG" + self.BODY + "GCTT"
        assert self._place(top, self.BODY) == (5, 4 + 4 + len(self.BODY), 1)
        assert top[5:8] == "ATG"

    def test_one_tier_release_is_flush_with_the_body(self):
        # No 3' overhang on this fragment — it lives at the START of the
        # next one (`_fragments_from_cuts`), so the body runs to the edge.
        top = "AATG" + self.BODY
        assert self._place(top, self.BODY) == (1, len(top), 1)

    def test_a_duplicated_overhang_picks_the_inner_one(self):
        """`_clone_part_synthesise_insert` builds `left_oh + oh5 + body`, so
        on a one-tier acceptor the overhang appears TWICE. The annotation
        belongs on the second."""
        top = "AATG" + "AATG" + self.BODY + "GCTT"
        assert self._place(top, self.BODY)[0] == 5
        assert top[5:8] == "ATG"

    def test_reverse_complemented_fragment_flips_the_strand(self):
        """`_clone_part_try_primer_path` RCs the candidate when the forward
        orientation finds no matching dropout."""
        top = sc._rc("CTCG" + "AATG" + self.BODY + "GCTT")
        start, end, strand = self._place(top, self.BODY)
        assert strand == -1
        # The ATG sits at the far end, reverse-complemented.
        assert top[end - 3:end] == "CAT"
        assert (end - start) % 3 == 0

    def test_a_body_repeating_with_the_overhang_period_is_not_mismatched(self):
        """The substring-vs-positional trap ([INV-172]'s guard, again): a bare
        leftmost `find` locks onto the FIRST copy and shifts the frame."""
        top = "AATG" + self.BODY + "AATG" + self.BODY
        start, _end, _s = self._place(top, self.BODY)
        assert start == 4 + len(self.BODY) + 1     # the second copy
        assert top[start:start + 3] == "ATG"

    def test_a_non_coding_part_covers_the_body_only(self):
        top = "CTCG" + "GGAG" + self.BODY + "AATG"
        start, end, _ = self._place(top, self.BODY, oh5="GGAG", oh3="AATG",
                                    ptype="Promoter")
        assert top[start:end] == self.BODY

    def test_mixed_case_still_matches(self):
        """The digest yields upper case, but a hand-built part dict need not
        — and a silent no-match falls back to the pre-fix whole span."""
        top = "CTCG" + "AATG" + self.BODY + "GCTT"
        assert self._place(top, self.BODY.lower()) == self._place(top, self.BODY)

    def test_an_unlocatable_body_falls_back_and_says_so(self):
        """Fail loud, not silent: the whole-fragment span IS the pre-[INV-173]
        placement, so reverting to it without a word is how the frameshift
        survived a release.

        NOT `caplog`: `_log` sets `propagate=False`, so pytest's handler on
        the root logger never sees these records — the assertion would pass
        vacuously against an empty list.
        """
        import logging
        seen: list[str] = []

        class _Capture(logging.Handler):
            def emit(self, record):
                seen.append(record.getMessage())

        h = _Capture()
        sc._log.addHandler(h)
        try:
            span = self._place("TTTT" * 20, self.BODY)
        finally:
            sc._log.removeHandler(h)
        assert span == (0, 80, 1)
        assert any("could not locate the body" in m for m in seen), seen

    def test_an_empty_body_does_not_invent_a_feature(self):
        """`_clone_part_into_entry_vector` bails before this, but the placer
        must not annotate a bare start codon as if it were the part."""
        top = "CTCG" + "AATG" + "GCTT"
        assert self._place(top, "") == (0, len(top), 1)
