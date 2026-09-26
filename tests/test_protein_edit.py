"""Protein-level editing — residue -> genomic bp, and the edit itself.

The mapping is the exact INVERSE of `_translate_cds`, so the central test is a
cross-check against it: rebuild the coding sequence from the position list and
assert the two translations agree, over every geometry that matters (both
strands x /codon_start 1-3 x spliced x origin wrap). Reading the mapping back to
yourself proves nothing; agreeing with the function that has shipped for months
does.

Then the edit: exactly one residue changes, to the intended amino acid, and the
DNA changes ONLY inside the three planned positions.
"""

from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

import splicecraft as sc


def _dna(n=600, seed=7):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def _rec(seq, feat, circular=True):
    r = SeqRecord(Seq(seq), id="p", name="p",
                  annotations={"molecule_type": "DNA",
                               "topology": "circular" if circular else "linear"})
    r.features = [feat]
    return r


def _cds(start, end, strand, *, codon_start=None):
    q = {}
    if codon_start is not None:
        q["codon_start"] = [str(codon_start)]
    return SeqFeature(FeatureLocation(start, end, strand), type="CDS",
                      qualifiers=q)


def _spliced(spans, strand):
    return SeqFeature(
        CompoundLocation([FeatureLocation(a, b, strand) for a, b in spans]),
        type="CDS")


def _wrap(tail, head, strand, n):
    return SeqFeature(
        CompoundLocation([FeatureLocation(tail, n, strand),
                          FeatureLocation(0, head, strand)]), type="CDS")


# (label, feature factory, exons for _translate_cds)
_GEOMETRIES = [
    ("plain +",        lambda n: (_cds(30, 330, 1), None)),
    ("plain -",        lambda n: (_cds(30, 330, -1), None)),
    ("codon_start2 +", lambda n: (_cds(30, 330, 1, codon_start=2), None)),
    ("codon_start2 -", lambda n: (_cds(30, 330, -1, codon_start=2), None)),
    ("codon_start3 +", lambda n: (_cds(30, 330, 1, codon_start=3), None)),
    ("codon_start3 -", lambda n: (_cds(30, 330, -1, codon_start=3), None)),
    ("spliced +",      lambda n: (_spliced([(30, 180), (300, 450)], 1),
                                  [(30, 180), (300, 450)])),
    # Minus-strand parts in READING order (descending) — what Biopython
    # stores for `complement(join(31..180,301..450))`. The ascending list
    # this fixture used is a DIFFERENT molecule to Biopython's `extract()`
    # (rc(A)+rc(B), an origin-crossing reading); SpliceCraft now reads parts
    # the way every exporter and external tool does (audit 2026-09-22).
    ("spliced -",      lambda n: (_spliced([(300, 450), (30, 180)], -1),
                                  [(30, 180), (300, 450)])),
    ("wrap +",         lambda n: (_wrap(500, 200, 1, n), None)),
    ("wrap -",         lambda n: (_wrap(500, 200, -1, n), None)),
]


def _translate_via_feature(seq, feat, exons):
    st, en, strand = sc._feat_bounds(feat, len(seq))
    return sc._translate_cds(seq, st, en, strand, exons=exons,
                             codon_start=sc._cds_codon_start(feat))


def _translate_via_positions(seq, feat):
    strand = sc._feat_bounds(feat, len(seq))[2]
    coding = "".join(
        sc._comp_base(seq[p]) if strand == -1 else seq[p].upper()
        for p in sc._cds_coding_positions(len(seq), feat))
    tbl = sc._codon_table_for(1)
    return "".join(tbl.get(coding[i:i + 3], "?")
                   for i in range(0, len(coding) - 2, 3))


class TestResidueMappingAgreesWithTranslation:

    @pytest.mark.parametrize("label,factory", _GEOMETRIES,
                             ids=[g[0] for g in _GEOMETRIES])
    def test_rebuilt_coding_sequence_matches(self, label, factory):
        seq = _dna()
        feat, exons = factory(len(seq))
        assert _translate_via_positions(seq, feat) == \
            _translate_via_feature(seq, feat, exons)

    @pytest.mark.parametrize("label,factory", _GEOMETRIES,
                             ids=[g[0] for g in _GEOMETRIES])
    def test_every_residue_codon_translates_to_that_residue(self, label,
                                                            factory):
        seq = _dna()
        feat, _exons = factory(len(seq))
        aa = _translate_via_positions(seq, feat)
        tbl = sc._codon_table_for(1)
        strand = sc._feat_bounds(feat, len(seq))[2]
        for r in range(1, len(aa) + 1):
            p3 = sc._residue_codon_positions(len(seq), feat, r)
            assert p3 is not None and len(p3) == 3, f"residue {r}"
            codon = "".join(
                sc._comp_base(seq[x]) if strand == -1 else seq[x].upper()
                for x in p3)
            assert tbl.get(codon) == aa[r - 1], f"residue {r}"

    def test_a_wrap_is_read_tail_first(self):
        """A wrap and a spliced CDS are BOTH two-part CompoundLocations, so the
        part ORDER must not decide which this is — trusting declared order put
        every GenBank-parsed wrap's arcs backwards once already."""
        seq = _dna()
        n = len(seq)
        feat = _wrap(500, 200, 1, n)
        pos = sc._cds_coding_positions(n, feat)
        assert pos[0] == 500, "must start at the tail, not bp 0"
        assert 500 <= pos[99] < n or pos[99] < 200
        assert len(pos) == 300

    def test_out_of_range_residue_is_none(self):
        seq = _dna()
        feat = _cds(30, 330, 1)
        assert sc._residue_codon_positions(len(seq), feat, 0) is None
        assert sc._residue_codon_positions(len(seq), feat, 101) is None
        assert sc._residue_codon_positions(len(seq), feat, 100) is not None

    def test_incomplete_trailing_codon_has_no_residue(self):
        """A 2-base remainder is not a residue; returning two positions would
        let a caller write a 2-base 'codon'."""
        seq = _dna()
        feat = _cds(30, 331, 1)          # 301 bases = 100 codons + 1
        assert sc._residue_codon_positions(len(seq), feat, 101) is None

    @pytest.mark.parametrize("cs,expect", [(1, 1), (2, 2), (3, 3),
                                           (0, 1), (9, 3), ("2", 2),
                                           (None, 1), ("junk", 1)])
    def test_codon_start_is_clamped(self, cs, expect):
        q = {} if cs is None else {"codon_start": [cs]}
        f = SeqFeature(FeatureLocation(0, 30, 1), type="CDS", qualifiers=q)
        assert sc._cds_codon_start(f) == expect

    def test_degenerate_inputs(self):
        assert sc._cds_coding_positions(0, _cds(0, 3, 1)) == []
        assert sc._cds_coding_positions(100, None) == []


class TestCodonChoice:

    def test_orders_by_usage(self):
        alts = sc._codon_choices_for_aa("A")
        assert alts and all(len(c) == 3 for c, _f in alts)
        fracs = [f for _c, f in alts]
        assert fracs == sorted(fracs, reverse=True)

    def test_every_residue_has_a_codon(self):
        for aa in sc._PROTEIN_AA_LETTERS + "*":
            assert sc._codon_choices_for_aa(aa), aa

    def test_unknown_letter_has_none(self):
        assert sc._codon_choices_for_aa("Z") == []
        assert sc._codon_choices_for_aa("") == []

    def test_comp_base_uses_the_iupac_table(self):
        """Sacred #3: one complement table, not a second hand-written map."""
        assert sc._comp_base("A") == "T"
        assert sc._comp_base("g") == "C"
        assert sc._comp_base("R") == "Y"      # IUPAC, not just ACGT
        assert sc._comp_base("N") == "N"


class TestEditPlan:

    def test_plan_reports_the_resolved_codon(self):
        seq = _dna()
        feat = _cds(30, 330, 1)
        rec = _rec(seq, feat)
        plan = sc._protein_edit_plan(rec, feat, 50, "W")
        assert plan["new_aa"] == "W" and plan["new_codon"] == "TGG"
        assert len(plan["positions"]) == 3
        assert len(plan["writes"]) == 3
        assert plan["strand"] == 1

    def test_silent_change_is_flagged_not_refused(self):
        """Swapping to a synonymous codon to kill a restriction site is a real
        thing people do on purpose."""
        seq = _dna()
        feat = _cds(30, 330, 1)
        rec = _rec(seq, feat)
        wt = sc._protein_edit_plan(rec, feat, 50, "W")["wt_aa"]
        plan = sc._protein_edit_plan(rec, feat, 50, wt)
        assert plan["silent"] is True

    def test_stop_creation_and_removal_are_flagged(self):
        seq = _dna()
        feat = _cds(30, 330, 1)
        rec = _rec(seq, feat)
        plan = sc._protein_edit_plan(rec, feat, 50, "*")
        assert plan["creates_stop"] is (plan["wt_aa"] != "*")

    def test_explicit_codon_is_honoured(self):
        seq = _dna()
        feat = _cds(30, 330, 1)
        plan = sc._protein_edit_plan(_rec(seq, feat), feat, 50, "A",
                                     codon="GCG")
        assert plan["new_codon"] == "GCG"

    def test_explicit_codon_must_encode_the_target(self):
        seq = _dna()
        feat = _cds(30, 330, 1)
        with pytest.raises(ValueError, match="codes"):
            sc._protein_edit_plan(_rec(seq, feat), feat, 50, "A",
                                  codon="TGG")

    @pytest.mark.parametrize("bad", ["GC", "GCGA", "GCN", "xyz"])
    def test_malformed_codon_is_refused(self, bad):
        seq = _dna()
        feat = _cds(30, 330, 1)
        with pytest.raises(ValueError, match="3 unambiguous"):
            sc._protein_edit_plan(_rec(seq, feat), feat, 50, "A", codon=bad)

    @pytest.mark.parametrize("bad", ["Z", "", "1", "XX"])
    def test_bad_target_residue_is_refused(self, bad):
        seq = _dna()
        feat = _cds(30, 330, 1)
        with pytest.raises(ValueError, match="amino-acid"):
            sc._protein_edit_plan(_rec(seq, feat), feat, 50, bad)

    def test_out_of_range_residue_names_the_protein_length(self):
        seq = _dna()
        feat = _cds(30, 330, 1)
        with pytest.raises(ValueError, match="100 complete"):
            sc._protein_edit_plan(_rec(seq, feat), feat, 500, "A")


class TestApplyEdit:

    @pytest.mark.parametrize("label,factory", _GEOMETRIES,
                             ids=[g[0] for g in _GEOMETRIES])
    @pytest.mark.parametrize("target", ["M", "W", "K", "*"])
    def test_exactly_one_residue_changes(self, label, factory, target):
        seq = _dna()
        feat, exons = factory(len(seq))
        rec = _rec(seq, feat)
        before = _translate_via_feature(seq, feat, exons)
        residue = max(1, len(before) // 2)
        plan = sc._protein_edit_plan(rec, feat, residue, target)
        new = sc._apply_protein_edit(rec, plan)
        new.features = [feat]
        after = _translate_via_feature(str(new.seq), feat, exons)

        assert len(after) == len(before)
        assert after[residue - 1] == target
        changed = [i for i, (a, b) in enumerate(zip(before, after)) if a != b]
        expect = [residue - 1] if before[residue - 1] != target else []
        assert changed == expect

    @pytest.mark.parametrize("label,factory", _GEOMETRIES,
                             ids=[g[0] for g in _GEOMETRIES])
    def test_dna_changes_only_inside_the_planned_codon(self, label, factory):
        seq = _dna()
        feat, _exons = factory(len(seq))
        rec = _rec(seq, feat)
        plan = sc._protein_edit_plan(rec, feat, 20, "W")
        new = sc._apply_protein_edit(rec, plan)
        diffs = {i for i, (a, b) in enumerate(zip(seq, str(new.seq))) if a != b}
        assert diffs <= {p for p, _b in plan["writes"]}

    def test_length_never_changes(self):
        seq = _dna()
        feat = _cds(30, 330, 1)
        rec = _rec(seq, feat)
        plan = sc._protein_edit_plan(rec, feat, 10, "W")
        assert len(sc._apply_protein_edit(rec, plan).seq) == len(seq)

    def test_caller_record_is_not_mutated(self):
        seq = _dna()
        feat = _cds(30, 330, 1)
        rec = _rec(seq, feat)
        plan = sc._protein_edit_plan(rec, feat, 10, "W")
        sc._apply_protein_edit(rec, plan)
        assert str(rec.seq) == seq

    def test_out_of_range_write_is_refused(self):
        seq = _dna(60)
        feat = _cds(0, 60, 1)
        rec = _rec(seq, feat)
        plan = sc._protein_edit_plan(rec, feat, 1, "W")
        plan["writes"] = [(999, "A")]
        with pytest.raises(ValueError, match="outside"):
            sc._apply_protein_edit(rec, plan)


class TestAgentEndpoints:
    """Agent parity for the protein-edit path."""

    class _App:
        _current_record = None
        _unsaved = False
        def _push_undo(self): self._undo = True
        def _apply_record(self, rec, clear_undo=True): self._current_record = rec
        def _mark_dirty(self): self._unsaved = True

    def _call(self, name, payload, rec=None):
        import splicecraft_state as st
        app = self._App()
        app._current_record = rec
        fn, _w = st._AGENT_HANDLERS[name]
        return (fn(app, payload), app)

    def _rec_two_cds(self):
        seq = _dna()
        f1 = _cds(30, 330, 1)
        f1.qualifiers["label"] = ["ampR"]
        f2 = _cds(400, 550, -1)
        f2.qualifiers["label"] = ["kanR"]
        r = SeqRecord(Seq(seq), id="p", name="p",
                      annotations={"molecule_type": "DNA",
                                   "topology": "circular"})
        r.features = [f1, f2]
        return r

    def test_endpoints_registered(self):
        import splicecraft_state as st
        assert "plan-residue-edit" in st._AGENT_HANDLERS
        assert st._AGENT_HANDLERS["edit-residue"][1] is True   # write=True

    def test_plan_by_cds_name(self):
        out, _app = self._call("plan-residue-edit",
                               {"cds": "ampR", "residue": 50, "to": "W"},
                               rec=self._rec_two_cds())
        assert out["new_codon"] == "TGG" and out["cds_label"] == "ampR"

    def test_plan_by_index(self):
        out, _app = self._call("plan-residue-edit",
                               {"cds_index": 1, "residue": 10, "to": "A"},
                               rec=self._rec_two_cds())
        assert out["strand"] == -1

    def test_ambiguous_cds_lists_the_choices(self):
        (body, code), _app = self._call(
            "plan-residue-edit", {"residue": 5, "to": "A"},
            rec=self._rec_two_cds())
        assert code == 400 and "choices" in body

    def test_unknown_cds_is_404(self):
        (body, code), _app = self._call(
            "plan-residue-edit", {"cds": "nope", "residue": 5, "to": "A"},
            rec=self._rec_two_cds())
        assert code == 404 and "available" in body["error"]

    def test_no_cds_features_is_422(self):
        r = SeqRecord(Seq(_dna()), id="p", name="p",
                      annotations={"molecule_type": "DNA"})
        (body, code), _app = self._call("plan-residue-edit",
                                        {"residue": 1, "to": "A"}, rec=r)
        assert code == 422 and "no CDS" in body["error"]

    def test_missing_to_is_400(self):
        (body, code), _app = self._call(
            "plan-residue-edit", {"cds": "ampR", "residue": 5},
            rec=self._rec_two_cds())
        assert code == 400 and "'to'" in body["error"]

    def test_edit_applies_and_reports(self):
        rec = self._rec_two_cds()
        before = str(rec.seq)
        out, app = self._call("edit-residue",
                              {"cds": "ampR", "residue": 50, "to": "W"},
                              rec=rec)
        assert out["applied"] is True
        assert out["bp_changed"] >= 1
        assert str(app._current_record.seq) != before
        assert len(str(app._current_record.seq)) == len(before)
        assert app._unsaved is True

    def test_a_load_right_after_the_read_is_not_overwritten(self):
        """The handler runs on the agent's thread. It read the record FIRST and
        captured the load counter after, so a GUI load landing between the two
        was invisible to the stale check and the edit — computed from the
        plasmid that had just been closed — was installed over the new one
        (audit 2026-09-22)."""
        import splicecraft_state as st
        rec, other = self._rec_two_cds(), self._rec_two_cds()

        class _RacingApp(self._App):
            _record_load_counter = 0
            installed = None
            _reads = 0

            @property
            def _current_record(self):
                self._reads += 1
                if self._reads == 1:           # the GUI loads right after
                    self._record_load_counter += 1
                    return rec
                return other

            def _apply_record(self, rec, clear_undo=True):
                self.installed = rec
        app = _RacingApp()
        fn, _w = st._AGENT_HANDLERS["edit-residue"]
        body, code = fn(app, {"cds": "ampR", "residue": 50, "to": "W"})
        assert code == 409 and app.installed is None

    def test_edit_refuses_a_bad_residue_without_touching_the_record(self):
        rec = self._rec_two_cds()
        before = str(rec.seq)
        (body, code), app = self._call(
            "edit-residue", {"cds": "ampR", "residue": 9999, "to": "W"},
            rec=rec)
        assert code == 400
        assert str(rec.seq) == before
