"""Regression tests for the 2026-10-24 assembly-design audit.

A read-only auditor drove the segmenter, the linter, the IIS cloner and the
synthesis-fragment designer against hand-derived ligations and brute-force
oracles. Each test below fails on the tree as it stood before this pass
(v1.3.11).

The recurring shape: ONE FUNCTION TREATING THE SAME QUESTION TWO WAYS. The
segment planner's junction-uniqueness check is reverse-complement-aware while
its off-target scan — three lines later, answering the same "does this
homology occur anywhere else" question — was forward-only, linear-only and
non-overlapping. `design-synthesis-fragment` collapses a duplicated start
codon on its typed path and not on its custom-overhang path. `homopolymer_max`
flagged a run of exactly max while `gc_window_lo`/`hi` beside it are
inclusive bounds.
"""
from __future__ import annotations

import random

import pytest

import splicecraft as sc
import splicecraft_biology as _bio
import splicecraft_seqanalysis as _sa

GB = sc._BUILTIN_GRAMMARS["gb_l0"]
_SITES = ("CGTCTC", "GAGACG", "GGTCTC", "GAGACC", "GAAGAC", "GTCTTC")


def _scrub(s: str) -> str:
    for st in _SITES:
        s = s.replace(st, "CTGCAG")
    return s


def _rseq(n: int, seed: int) -> str:
    r = random.Random(seed)
    return _scrub("".join(r.choice("ACGT") for _ in range(n)))


def _handler(name):
    return sc._state._AGENT_HANDLERS[name][0]


# ══════════════════════════════════════════════════════════════════════════
# Off-target homology: the inconsistent half of the same function
# ══════════════════════════════════════════════════════════════════════════

class TestOffTargetHomologyIsDoubleStrandedAndCircular:
    """`extra = s.count(j["overlap"]) - 1` is forward-only, blind to the
    origin, and non-overlapping, while the junction-vs-junction uniqueness
    check in the SAME function already tests `_rc(...)`. Gibson and yeast-HR
    homology is double-stranded, and the docstring promises "whether the
    overlap occurs ANYWHERE ELSE in the target"."""

    def _plan(self, seq, **kw):
        body = {"sequence": seq, "segment_bp": 3_000, "strategy": "gibson"}
        body.update(kw)
        r = _handler("plan-chromosome-segments")(None, body)
        assert not isinstance(r, tuple), r
        return r

    def test_a_reverse_complement_copy_is_reported(self):
        seq = _rseq(9_000, 3)
        plan = self._plan(seq, circular=False)
        ov = plan["junctions"][0]["overlap"]
        # Plant the RC of junction 0's homology well away from any junction.
        at = 7_000
        planted = seq[:at] + sc._rc(ov) + seq[at + len(ov):]
        got = self._plan(planted, circular=False)
        j0 = got["junctions"][0]
        assert j0.get("offtarget_hits"), j0
        assert j0.get("offtarget_hits_reverse") == 1, j0
        assert any("INVERTED" in w for w in j0["warnings"]), j0["warnings"]

    def test_a_copy_across_the_origin_is_reported(self):
        seq = _rseq(9_000, 5)
        plan = self._plan(seq, circular=True)
        ov = plan["junctions"][0]["overlap"]
        half = len(ov) // 2
        # A second copy straddling bp 0: invisible to a scan of the
        # linearised string (sacred invariant #6).
        planted = (ov[half:] + seq[len(ov) - half:-half] + ov[:half])
        assert len(planted) == len(seq)
        got = self._plan(planted, circular=True)
        hits = [j for j in got["junctions"] if j.get("offtarget_hits")]
        assert hits, [j["overlap"] for j in got["junctions"]]

    def test_a_plain_forward_copy_is_still_reported(self):
        """The control: if this breaks, the fix broke the working half."""
        seq = _rseq(9_000, 7)
        plan = self._plan(seq, circular=False)
        ov = plan["junctions"][0]["overlap"]
        at = 7_000
        planted = seq[:at] + ov + seq[at + len(ov):]
        got = self._plan(planted, circular=False)
        j0 = got["junctions"][0]
        assert j0.get("offtarget_hits_forward") == 1, j0

    def test_overlapping_occurrences_are_all_counted(self):
        """`"ATATATATAT".count("ATAT")` is 2 where 4 occur — which matters at
        Golden Gate's 4 nt."""
        import re
        pat = re.compile(re.escape("ATAT"))
        assert len(list(_bio._iter_match_starts(pat, "ATATATATAT"))) == 4


class TestThePlannerSaysWhatItDroppedFromTheTarget:
    """Non-ACGT bases were removed before planning and nothing said so: every
    coordinate then indexed the filtered sequence while the caller reads them
    against the molecule they passed, and the segments together spelled a
    molecule missing those bases — a plan no vendor can fill and that would
    not reassemble into the input."""

    def test_removed_bases_are_reported(self):
        seq = _rseq(6_000, 11)
        with_n = seq[:1_000] + "N" * 25 + seq[1_025:]
        r = _handler("plan-chromosome-segments")(
            None, {"sequence": with_n, "segment_bp": 2_000,
                   "strategy": "gibson", "circular": False})
        assert not isinstance(r, tuple), r
        assert r["n_degenerate"] == 25
        assert r["input_bp"] == len(with_n)
        assert r["total_bp"] == len(with_n) - 25
        assert any("REMOVED" in w for w in r["warnings"]), r["warnings"]

    def test_a_clean_target_says_nothing(self):
        r = _handler("plan-chromosome-segments")(
            None, {"sequence": _rseq(6_000, 13), "segment_bp": 2_000,
                   "strategy": "gibson", "circular": False})
        assert not isinstance(r, tuple), r
        assert r["n_degenerate"] == 0
        assert not [w for w in r["warnings"] if "REMOVED" in w]


# ══════════════════════════════════════════════════════════════════════════
# Threshold naming
# ══════════════════════════════════════════════════════════════════════════

class TestHomopolymerMaxIsTheLongestRunAllowed:
    """`homopolymer_max=9` flagged a run of exactly 9, reading the name as
    "minimum to flag", while `gc_window_lo` / `gc_window_hi` right beside it
    are inclusive bounds (exactly 20% / 80% is not flagged). One of the two
    names had to be wrong."""

    def _runs(self, run_len, **kw):
        body = _rseq(400, 17).replace("A", "C")
        raw = body[:200] + "A" * run_len + body[200:]
        r = _handler("lint-synthesis")(
            None, {"sequence": raw, "circular": False, **kw})
        assert not isinstance(r, tuple), r
        return [w for w in r["warnings"] if w["kind"] == "homopolymer"]

    def test_a_run_of_exactly_max_is_allowed(self):
        assert self._runs(9, homopolymer_max=9) == []

    def test_a_run_one_longer_is_flagged(self):
        runs = self._runs(10, homopolymer_max=9)
        assert len(runs) == 1 and "10×A" in runs[0]["message"], runs

    def test_the_default_still_flags_a_real_hazard(self):
        runs = self._runs(14)
        assert len(runs) == 1 and "14×A" in runs[0]["message"], runs

    def test_the_gc_window_bounds_stay_inclusive(self):
        """The parameter this one was made consistent WITH."""
        import inspect
        sig = inspect.signature(_sa._synthesis_lint)
        assert sig.parameters["gc_window_lo"].default == 0.20
        assert sig.parameters["gc_window_hi"].default == 0.80
        assert sig.parameters["homopolymer_max"].default == 9


# ══════════════════════════════════════════════════════════════════════════
# design-synthesis-fragment: the two paths disagreed about the start codon
# ══════════════════════════════════════════════════════════════════════════

class TestACustomOverhangSaysWhenItDoublesTheStartCodon:
    """A 4+ nt 5' overhang ending in ATG CARRIES the start codon (Golden
    Braid's AATG = A spacer + ATG). The typed path collapses a body that also
    opens on ATG; the custom-overhang path has no part type, so the same
    insert with the same overhangs came out 3 bp longer reading
    `CTCG·AATG·ATG·…` — an extra Met, [INV-136] through another door. Custom
    overhangs are verbatim by design, so this is SAID rather than silently
    fixed."""

    BODY = "ATG" + _rseq(150, 19) + "TAA"

    def _design(self, **kw):
        body = {"sequence": self.BODY, "grammar": "gb_l0"}
        body.update(kw)
        r = _handler("design-synthesis-fragment")(None, body)
        assert not isinstance(r, tuple), r
        return r

    def test_the_custom_path_warns(self):
        r = self._design(oh5="AATG", oh3="GCTT")
        warns = r["result"].get("warnings") or []
        assert any("two (M-M)" in w for w in warns), warns

    def test_the_typed_path_does_not_warn(self):
        """It collapses the duplicate instead, so there is nothing to say."""
        r = self._design(part_type="CDS")
        warns = r["result"].get("warnings") or []
        assert not [w for w in warns if "M-M" in w], warns

    def test_a_body_not_opening_on_atg_does_not_warn(self):
        r = _handler("design-synthesis-fragment")(
            None, {"sequence": _rseq(120, 23), "grammar": "gb_l0",
                   "oh5": "AATG", "oh3": "GCTT"})
        assert not isinstance(r, tuple), r
        warns = r["result"].get("warnings") or []
        assert not [w for w in warns if "M-M" in w], warns

    def test_a_non_atg_overhang_does_not_warn(self):
        r = self._design(oh5="AGGT", oh3="GCTT")
        warns = r["result"].get("warnings") or []
        assert not [w for w in warns if "M-M" in w], warns


# ══════════════════════════════════════════════════════════════════════════
# The chain cloner resolved each fragment's part type by NAME
# ══════════════════════════════════════════════════════════════════════════

def _mk_vector(name, site, oh5, oh3, seed=1):
    """Circular acceptor: backbone + SITE+A + oh5 + dropout + oh3 + rc(SITE+A),
    carrying a marker CDS so backbone-marker exclusion finds the stuffer."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    bb = _rseq(1_400, seed)
    dropout = _rseq(220, seed + 77)
    left = site + "A"
    cassette = left + oh5 + dropout + oh3 + sc._rc(left)
    seq = bb[:700] + cassette + bb[700:]
    rec = SeqRecord(Seq(seq), id=name, name=name,
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    rec.features.append(SeqFeature(
        FeatureLocation(0, 600, 1), type="CDS",
        qualifiers={"label": ["AmpR"], "gene": ["bla"]}))
    rec.features.append(SeqFeature(
        FeatureLocation(len(bb[:700]) + len(cassette) + 10,
                        len(bb[:700]) + len(cassette) + 300, 1),
        type="rep_origin", qualifiers={"label": ["ori"]}))
    return {"name": name, "id": name, "size": len(seq),
            "gb_text": sc._record_to_gb_text(rec)}


class TestTheChainClonerResolvesPartTypeByPosition:
    """`fragments` is built in `sources` order, one per source, so the index
    IS the mapping — but the type was looked up by NAME. A nameless source
    becomes "src", matched nothing, left the type empty and silently skipped
    the AATG start-codon collapse: an assembled CDS with two start codons,
    [INV-136] again."""

    P_BODY = _rseq(90, 31)
    C_BODY = "ATG" + _rseq(150, 32) + "TAA"
    T_BODY = _rseq(120, 33)

    def _parts(self, cds_name="C1"):
        def l0(name, ptype, oh5, oh3, body):
            return {"name": name, "type": ptype, "oh5": oh5, "oh3": oh3,
                    "sequence": body, "level": 0, "grammar": "gb_l0"}
        return [l0("P1", "Promoter", "GGAG", "AATG", self.P_BODY),
                l0(cds_name, "CDS", "AATG", "GCTT", self.C_BODY),
                l0("T1", "Terminator", "GCTT", "CGCT", self.T_BODY)]

    @pytest.mark.parametrize("cds_name", ["C1", ""])
    def test_the_start_codon_collapses_even_for_a_nameless_source(
            self, cds_name):
        alpha = _mk_vector("ALPHA1", "CGTCTC", "GGAG", "CGCT", seed=1)
        res = sc._clone_assembly_into_entry_vector(
            self._parts(cds_name), alpha, GB, source_level=0, name="TU1")
        assert res is not None
        prod = str(res.seq).upper()
        circ = prod + prod
        collapsed = ("GGAG" + self.P_BODY + "AATG" + self.C_BODY[3:]
                     + "GCTT" + self.T_BODY + "CGCT")
        doubled = ("GGAG" + self.P_BODY + "AATG" + self.C_BODY
                   + "GCTT" + self.T_BODY + "CGCT")
        assert collapsed in circ, (cds_name, "collapse did not run")
        assert doubled not in circ, (cds_name, "two start codons")


class TestAnInternalStepEnzymeSiteIsNeverSilent:
    """A part whose body carries the STEP's own enzyme site is cut inside the
    one-pot reaction, so the product that assembles is not the one drawn. The
    assembler never consulted the grammar's own `forbidden_sites` machinery.
    Warned rather than refused (the user's call for this class), carried on
    the record so the Constructor toast and the agent response say the same
    thing."""

    def test_the_chain_path_warns_with_the_position(self):
        def l0(name, ptype, oh5, oh3, body):
            return {"name": name, "type": ptype, "oh5": oh5, "oh3": oh3,
                    "sequence": body, "level": 0, "grammar": "gb_l0"}
        dirty = l0("Cdirty", "CDS", "AATG", "GCTT",
                   "ATG" + _rseq(60, 41) + "CGTCTC" + _rseq(60, 42) + "TAA")
        parts = [l0("P1", "Promoter", "GGAG", "AATG", _rseq(90, 43)),
                 dirty,
                 l0("T1", "Terminator", "GCTT", "CGCT", _rseq(120, 44))]
        alpha = _mk_vector("ALPHA1", "CGTCTC", "GGAG", "CGCT", seed=2)
        res = sc._clone_assembly_into_entry_vector(
            parts, alpha, GB, source_level=0, name="TU1")
        assert res is not None
        warns = list(getattr(res, "_sc_clone_warnings", []) or [])
        assert any("internal Esp3I site" in w for w in warns), warns
        assert any("cutting at +" in w for w in warns), warns
