"""Bacterial regulatory scanners + transcription mapping.

`splicecraft_regulatory` (sigma-70 promoters, intrinsic terminators),
`splicecraft_seqanalysis._map_transcription`, and the five agent endpoints
added with them.

Three things are checked against something OUTSIDE the implementation, because
reading a scanner back to yourself agrees with whatever it already does:

  * REAL sequences — the T7 Tphi terminator and this repo's own verified
    preset catalogue (extracted from public GenBank records), not hand-typed
    motifs.
  * ROTATION and STRAND invariance — a circular molecule has no origin, so the
    same construct must answer identically however it is drawn. This is what
    caught the terminator scanner clamping its upstream search at index 0 and
    losing every hairpin that straddled bp 0.
  * NON-REDUNDANCY — one hairpin is one record. The field report that prompted
    this work had a hand-rolled scan report 22 terminators where there was
    one, because each overlapping stem/loop register scored separately.
"""

from __future__ import annotations

import random

import pytest

import splicecraft as sc
import splicecraft_regulatory as reg
import splicecraft_presets as presets
import splicecraft_seqanalysis as sa


_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _rand(n: int, seed: int, alphabet: str = "ACGT") -> str:
    r = random.Random(seed)
    return "".join(r.choice(alphabet) for _ in range(n))


# A textbook-perfect sigma-70 promoter: consensus -35, 17 bp spacer,
# consensus -10.
_PERFECT_PROMOTER = "TTGACA" + _rand(17, 1) + "TATAAT"
# A designed intrinsic terminator: 10 bp stem, 4 nt loop, 8 T tract.
_STEM = "GCCGCCAGTT"
_TERMINATOR = _STEM + "CGAA" + _rc(_STEM) + "TTTTTTTT"


# ══════════════════════════════════════════════════════════════════════
#  Promoters
# ══════════════════════════════════════════════════════════════════════

class TestPromoterScan:

    def test_perfect_consensus_scores_one(self):
        seq = _rand(200, 2) + _PERFECT_PROMOTER + _rand(200, 3)
        top = reg._scan_promoters(seq, circular=False)[0]
        assert top["start"] == 200
        assert top["minus35"] == "TTGACA"
        assert top["minus10"] == "TATAAT"
        assert top["spacer"] == 17
        assert top["score"] == pytest.approx(1.0)
        assert top["strength"] == "strong"

    def test_tss_sits_downstream_of_the_minus10(self):
        seq = _rand(200, 2) + _PERFECT_PROMOTER + _rand(200, 3)
        top = reg._scan_promoters(seq, circular=False)[0]
        # The -10 hexamer sits at -12..-7 and the start is +1, so the TSS is 12
        # bases past the hexamer's FIRST base — i.e. 7 past its last, which is
        # the spacing the literature quotes. Written out from the textbook
        # numbering rather than from the module's constant on purpose.
        minus10_start = 200 + 6 + 17
        assert top["minus10_start"] == minus10_start
        assert top["tss"] == minus10_start + 12
        assert top["tss"] == (minus10_start + 5) + 7   # same, from the 3' end

    def test_a_degraded_hexamer_scores_below_a_perfect_one(self):
        good = _rand(120, 4) + _PERFECT_PROMOTER + _rand(120, 5)
        bad_el = "TTGACA" + _rand(17, 1) + "TACAAG"      # -10 broken
        bad = _rand(120, 4) + bad_el + _rand(120, 5)
        g = reg._scan_promoters(good, circular=False, min_score=0.0)[0]
        b = reg._scan_promoters(bad, circular=False, min_score=0.0)[0]
        assert b["score"] < g["score"]

    def test_components_are_reported_so_the_score_shows_its_working(self):
        seq = _rand(100, 6) + _PERFECT_PROMOTER + _rand(100, 7)
        top = reg._scan_promoters(seq, circular=False)[0]
        c = top["components"]
        assert c["minus35_match"] == pytest.approx(1.0)
        assert c["minus10_match"] == pytest.approx(1.0)
        assert c["spacer_penalty"] == pytest.approx(0.0)

    def test_the_score_says_it_is_not_calibrated(self):
        """A number that looks trained and is not would be worse than none —
        the same discipline `_score_guide` and the splice model apply."""
        seq = _rand(100, 8) + _PERFECT_PROMOTER + _rand(100, 9)
        top = reg._scan_promoters(seq, circular=False)[0]
        assert "not calibrated" in top["model"]

    def test_minus_strand_hit_reports_forward_coordinates(self):
        """Sacred invariant #2: a reverse hit is stored at its FORWARD
        coordinate, not mirrored."""
        seq = _rand(200, 10) + _rc(_PERFECT_PROMOTER) + _rand(200, 11)
        hits = [h for h in reg._scan_promoters(seq, circular=False)
                if h["strand"] == -1]
        assert hits and hits[0]["score"] == pytest.approx(1.0)
        assert hits[0]["start"] == 200
        # A minus-strand promoter fires LEFTWARDS, so its TSS is at a lower
        # coordinate than the element itself.
        assert hits[0]["tss"] < hits[0]["start"]

    def test_both_strands_false_finds_only_forward(self):
        seq = _rand(150, 12) + _rc(_PERFECT_PROMOTER) + _rand(150, 13)
        assert not [h for h in reg._scan_promoters(
            seq, circular=False, both_strands=False) if h["strand"] == -1]

    @pytest.mark.parametrize("rot", [1, 17, 150, 299])
    def test_rotation_invariant(self, rot):
        circ = _PERFECT_PROMOTER + _rand(300, 14)
        n = len(circ)
        base = sorted((h["start"], h["score"], h["strand"])
                      for h in reg._scan_promoters(circ, circular=True))
        r = circ[rot:] + circ[:rot]
        got = sorted(((h["start"] + rot) % n, h["score"], h["strand"])
                     for h in reg._scan_promoters(r, circular=True))
        assert got == base

    def test_reverse_complementing_the_molecule_finds_the_same_promoter(self):
        seq = _rand(200, 15) + _PERFECT_PROMOTER + _rand(200, 16)
        fwd = max(h["score"] for h in reg._scan_promoters(seq,
                                                          circular=False))
        rev = max(h["score"] for h in reg._scan_promoters(_rc(seq),
                                                          circular=False))
        assert fwd == pytest.approx(rev)

    def test_real_preset_promoters_score_above_the_default_floor(self):
        """Calibration anchor. The default `min_score` is set from these
        numbers; if the scoring drifts, this is where it shows."""
        wanted = {"trc promoter", "J23119 promoter", "tac promoter"}
        seen = {}
        for f in presets._preset_features():
            if f.get("feature_type") != "promoter" or f["name"] not in wanted:
                continue
            hits = reg._scan_promoters(f["sequence"], circular=False,
                                       min_score=0.0)
            seen[f["name"]] = hits[0]["score"] if hits else 0.0
        assert seen, "no promoter presets found — catalogue changed?"
        for name, score in seen.items():
            assert score >= 0.90, f"{name} scored {score}"

    def test_default_floor_keeps_random_dna_quiet(self):
        """The floor is calibrated, not arbitrary: ~1 chance hit per kb per
        strand. A regression that made it permissive again would bury real
        hits in noise, which is what 245 hits per 2.7 kb looked like."""
        seq = _rand(20_000, 17)
        hits = reg._scan_promoters(seq, circular=False, both_strands=False)
        per_kb = len(hits) / 20.0
        assert per_kb < 3.0, f"{per_kb:.2f} hits/kb is too noisy"

    def test_empty_and_short_inputs(self):
        assert reg._scan_promoters("") == []
        assert reg._scan_promoters("ACGT") == []

    def test_ambiguity_codes_do_not_match(self):
        """An N is not evidence. A promoter call resting on a base the record
        does not know would be a guess presented as a finding."""
        masked = "NNNNNN" + _rand(17, 1) + "NNNNNN"
        seq = _rand(100, 18) + masked + _rand(100, 19)
        assert not [h for h in reg._scan_promoters(seq, circular=False)
                    if h["start"] == 100]


# ══════════════════════════════════════════════════════════════════════
#  Terminators
# ══════════════════════════════════════════════════════════════════════

class TestTerminatorScan:

    def test_designed_hairpin_geometry_is_exact(self):
        seq = _rand(150, 20) + _TERMINATOR + _rand(150, 21)
        hit = [h for h in reg._scan_terminators(seq, circular=False)
               if h["strand"] == 1 and h["start"] == 150]
        assert len(hit) == 1
        h = hit[0]
        assert h["stem_len"] == 10
        assert h["loop_len"] == 4
        assert h["u_tract"] >= 8
        assert h["dg"] is not None and h["dg"] < 0
        assert h["efficiency"] == "strong"

    def test_real_t7_terminator_is_found_once(self):
        """A real third-party sequence, not a hand-built motif."""
        t7 = "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG"
        seq = _rand(120, 22) + t7 + _rand(120, 23)
        hits = [h for h in reg._scan_terminators(seq, circular=False)
                if h["strand"] == 1 and 120 <= h["start"] < 120 + len(t7)]
        assert len(hits) == 1, f"T7 terminator reported {len(hits)} times"
        assert hits[0]["u_tract"] >= 6
        assert hits[0]["dg"] < -15
        assert hits[0]["efficiency"] == "strong"

    def test_real_preset_terminators_are_all_found_and_called_strong(self):
        """Calibration anchor for the `strong` tier — the tier that makes
        `map-transcription` call a read-through BLOCKED."""
        wanted = {"rrnB T1 terminator", "T7 terminator",
                  "lambda t0 terminator", "T7Te terminator",
                  "BBa_B1008 terminator"}
        found = {}
        for f in presets._preset_features():
            if f.get("feature_type") != "terminator" or f["name"] not in wanted:
                continue
            hits = reg._scan_terminators(f["sequence"], circular=False,
                                         both_strands=False)
            found[f["name"]] = hits
        assert found, "no terminator presets found — catalogue changed?"
        for name, hits in found.items():
            assert hits, f"{name} was missed entirely"
            best = max(hits, key=lambda h: h["stem_len"])
            assert best["efficiency"] == "strong", f"{name}: {best}"

    def test_eukaryotic_terminators_are_not_called(self):
        """CYC1 / NOS / ADH1 are not Rho-independent terminators. Calling them
        would mean the scanner is matching AT-richness, not hairpins."""
        for f in presets._preset_features():
            if f.get("feature_type") != "terminator":
                continue
            if f["name"] not in ("CYC1 terminator", "NOS terminator",
                                 "ADH1 terminator"):
                continue
            hits = reg._scan_terminators(f["sequence"], circular=False,
                                         both_strands=False)
            assert not hits, f"{f['name']} was called as intrinsic: {hits}"

    def test_a_non_folding_stem_is_not_reported(self):
        """Bases can be written as complementary and still not form. A
        candidate that folds to ~0 kcal/mol is not a terminator."""
        for h in reg._scan_terminators(_rand(5000, 24), circular=False):
            assert h["dg"] is None or h["dg"] < 0

    def test_minus_strand_hit_reports_forward_coordinates(self):
        """Sacred invariant #2. On the forward strand a minus-strand
        terminator reads A-tract first, then the hairpin — so its forward
        `start` is the tract, at the LOW coordinate. The tract may run a base
        or two further left if the flank happens to end in A, so the test
        pins the hairpin geometry and bounds the start rather than demanding
        an exact offset the flanking sequence can legitimately move."""
        seq = _rand(150, 25) + _rc(_TERMINATOR) + _rand(150, 26)
        hits = [h for h in reg._scan_terminators(seq, circular=False)
                if h["strand"] == -1 and h["stem_len"] == 10]
        assert len(hits) == 1
        assert hits[0]["u_tract"] >= 8
        assert hits[0]["start"] <= 150
        assert 150 < hits[0]["end"] <= 150 + len(_TERMINATOR)

    @pytest.mark.parametrize("rot", [1, 5, 60, 351])
    def test_rotation_invariant_including_an_origin_spanning_hairpin(self, rot):
        """THE regression test. The terminator search runs BACKWARDS from the
        U-tract, so a suffix-only wrap window clamped at index 0 and either
        truncated the stem or lost the terminator outright — invisible until
        the molecule was rotated."""
        circ = _TERMINATOR + _rand(320, 27)
        n = len(circ)
        base = sorted((h["start"], h["strand"], h["stem_len"], h["u_tract"])
                      for h in reg._scan_terminators(circ, circular=True))
        r = circ[rot:] + circ[:rot]
        got = sorted(((h["start"] + rot) % n, h["strand"], h["stem_len"],
                      h["u_tract"])
                     for h in reg._scan_terminators(r, circular=True))
        assert got == base

    def test_origin_spanning_terminator_keeps_its_full_stem(self):
        circ = _TERMINATOR + _rand(320, 28)
        n = len(circ)
        rotated = circ[1:] + circ[:1]
        hits = [h for h in reg._scan_terminators(rotated, circular=True)
                if h["stem_len"] == 10]
        assert hits, "the origin-spanning hairpin lost its stem"

    def test_a_t_free_strand_has_no_terminator(self):
        clean = _rand(400, 29, alphabet="ACG")
        assert reg._scan_terminators(clean, circular=False,
                                     both_strands=False) == []

    def test_empty_input(self):
        assert reg._scan_terminators("") == []


class TestMinusStrandAgainstTheBottomStrand:
    """An INDEPENDENT check, written from the biology rather than from the
    scanner: whatever span is reported for a minus-strand hit, reading those
    bases on the BOTTOM strand must produce the element itself.

    Reading a coordinate convention back to yourself agrees with whatever the
    code already does. This does not — it reconstructs the element from the
    sequence the span points at, so an off-by-one or a mirrored coordinate
    (sacred invariant #2) fails here even when the scanner is self-consistent.
    """

    @pytest.mark.parametrize("trial", range(12))
    def test_a_minus_strand_promoter_span_reads_back_as_the_promoter(
            self, trial):
        r = random.Random(1000 + trial)
        spacer = "".join(r.choice("ACGT") for _ in range(17))
        el = "TTGACA" + spacer + "TATAAT"
        left = "".join(r.choice("ACGT") for _ in range(r.randint(40, 200)))
        right = "".join(r.choice("ACGT") for _ in range(r.randint(40, 200)))
        seq = left + _rc(el) + right
        hits = [h for h in reg._scan_promoters(seq, circular=False)
                if h["strand"] == -1 and h["score"] > 0.99]
        assert hits, "the planted minus-strand promoter was not found"
        for h in hits:
            span = seq[h["start"]:h["end"]]
            bottom = _rc(span)
            assert bottom.startswith(h["minus35"])
            assert bottom.endswith(h["minus10"])

    @pytest.mark.parametrize("trial", range(12))
    def test_a_minus_strand_terminator_span_ends_in_its_u_tract(self, trial):
        r = random.Random(2000 + trial)
        # A GC-RICH stem, as every real intrinsic terminator has. A stem drawn
        # from all four bases is often AT-rich enough to fold above the -8
        # kcal/mol reporting floor, and is then correctly NOT called — which
        # would make this test fail for a reason that has nothing to do with
        # the coordinates it exists to check.
        stem = "".join(r.choice("GCGCAT") for _ in range(10))
        el = stem + "CGAA" + _rc(stem) + "T" * 8
        left = "".join(r.choice("ACGT") for _ in range(r.randint(60, 200)))
        right = "".join(r.choice("ACGT") for _ in range(r.randint(60, 200)))
        seq = left + _rc(el) + right
        planted = len(left)
        # Selected by POSITION, not by stem length. When `rc(stem)` happens to
        # end in T the scanner reads a 9 bp stem with a 9 T tract rather than
        # 10 with 8 — both are valid readings of the same molecule and it
        # cannot know which base the designer meant as stem. The coordinate is
        # what this test is about.
        hits = [h for h in reg._scan_terminators(seq, circular=False)
                if h["strand"] == -1 and abs(h["start"] - planted) <= 2]
        assert hits, "the planted minus-strand terminator was not found"
        for h in hits:
            bottom = _rc(seq[h["start"]:h["end"]])
            # A terminator's U-tract lies 3' of its hairpin, so on the strand
            # the element is actually read from, the span must END in it.
            assert bottom.endswith("T" * h["u_tract"])


class TestNonRedundancy:
    """One hairpin is ONE record. The bug this module exists to not have."""

    def test_one_designed_hairpin_yields_one_record(self):
        seq = _rand(200, 30) + _TERMINATOR + _rand(200, 31)
        overlapping = [h for h in reg._scan_terminators(seq, circular=False)
                       if h["strand"] == 1
                       and 200 <= h["start"] < 200 + len(_TERMINATOR)]
        assert len(overlapping) == 1

    def test_no_two_same_strand_hits_overlap(self):
        """The general property, over a molecule with many candidates."""
        seq = _rand(8000, 32)
        for strand in (1, -1):
            hits = sorted((h for h in reg._scan_terminators(seq,
                                                            circular=False)
                           if h["strand"] == strand),
                          key=lambda h: h["start"])
            for a, b in zip(hits, hits[1:]):
                assert b["start"] >= a["end"], f"{a} overlaps {b}"

    def test_a_longer_u_tract_does_not_multiply_the_record(self):
        """Every offset into a T-run is a candidate tract start; only the run's
        real start may anchor one."""
        long_tract = _STEM + "CGAA" + _rc(_STEM) + "T" * 9
        seq = _rand(150, 33) + long_tract + _rand(150, 34)
        hits = [h for h in reg._scan_terminators(seq, circular=False)
                if h["strand"] == 1 and 150 <= h["start"] < 150 + len(long_tract)]
        assert len(hits) == 1


# ══════════════════════════════════════════════════════════════════════
#  map-transcription
# ══════════════════════════════════════════════════════════════════════

def _readthrough_construct(with_insulator: bool = False):
    """The reporter's topology: marker promoter -> marker -> ori -> phage
    promoter -> cargo, with no designed insulator anywhere."""
    kan_p = _PERFECT_PROMOTER
    kan_cds = _rand(800, 40)
    insulator = _TERMINATOR if with_insulator else ""
    ori = _rand(600, 41)
    t7_p = "TAATACGACTCACTATAG"
    cargo = _rand(700, 42)
    seq = kan_p + kan_cds + insulator + ori + t7_p + cargo + _rand(200, 43)
    a = len(kan_p)
    b = a + len(kan_cds)
    c = b + len(insulator)
    d = c + len(ori)
    e = d + len(t7_p)
    feats = [
        {"start": 0, "end": a, "strand": 1, "type": "promoter",
         "label": "KanR promoter"},
        {"start": a, "end": b, "strand": 1, "type": "CDS", "label": "KanR"},
        {"start": d, "end": e, "strand": 1, "type": "promoter",
         "label": "T7 promoter"},
        {"start": e, "end": e + len(cargo), "strand": 1, "type": "CDS",
         "label": "GeneX"},
    ]
    if with_insulator:
        feats.insert(2, {"start": b, "end": c, "strand": 1,
                         "type": "terminator", "label": "T1"})
    return seq, feats


def _gene(res, label):
    return next(g for g in res["genes"] if g["cds"]["label"] == label)


class TestMapTranscription:

    def test_a_marker_promoter_reaches_the_cargo_with_nothing_between(self):
        """The finding the whole endpoint exists for: checking terminators
        DOWNSTREAM of a cassette cannot see a read-through arriving from
        behind, around the circle."""
        seq, feats = _readthrough_construct()
        res = sa._map_transcription(seq, feats, circular=True)
        cargo_gene = _gene(res, "GeneX")
        hit = [r for r in cargo_gene["reachable_from"]
               if r["promoter"]["label"] == "KanR promoter"]
        assert hit, "the read-through was missed"
        assert hit[0]["readthrough"] == "clear"
        assert hit[0]["terminators_between"] == []
        assert not cargo_gene["silent"]

    def test_an_intervening_strong_terminator_blocks_it(self):
        seq, feats = _readthrough_construct(with_insulator=True)
        res = sa._map_transcription(seq, feats, circular=True)
        cargo_gene = _gene(res, "GeneX")
        hit = [r for r in cargo_gene["reachable_from"]
               if r["promoter"]["label"] == "KanR promoter"]
        assert hit and hit[0]["readthrough"] == "blocked"
        assert hit[0]["terminators_between"]

    def test_a_phage_promoter_is_not_host_recognised(self):
        """`silent` and `silent_in_host` are different claims, and the gap
        between them is the point."""
        seq, feats = _readthrough_construct()
        res = sa._map_transcription(seq, feats, circular=True)
        by_label = {p["label"]: p for p in res["promoters"]
                    if p["source"] == "annotated"}
        assert by_label["T7 promoter"]["host_recognised"] is False
        assert "T7" in by_label["T7 promoter"]["host_basis"]
        assert by_label["KanR promoter"]["host_recognised"] is True

    def test_the_open_construct_is_not_silent_in_a_host_without_t7(self):
        seq, feats = _readthrough_construct()
        res = sa._map_transcription(seq, feats, circular=True)
        cargo_gene = _gene(res, "GeneX")
        assert cargo_gene["silent_in_host"] is False
        assert cargo_gene["n_host_sources"] >= 1

    def test_strand_is_respected(self):
        """A + promoter cannot drive a - CDS."""
        seq = _PERFECT_PROMOTER + _rand(500, 44)
        feats = [
            {"start": 0, "end": len(_PERFECT_PROMOTER), "strand": 1,
             "type": "promoter", "label": "P"},
            {"start": 100, "end": 400, "strand": -1, "type": "CDS",
             "label": "rev gene"},
        ]
        res = sa._map_transcription(seq, feats, circular=True,
                                    include_predicted=False)
        g = _gene(res, "rev gene")
        assert all(r["promoter"]["strand"] == -1
                   for r in g["reachable_from"])

    def test_a_linear_molecule_does_not_wrap(self):
        """On a linear molecule a promoter downstream of the gene simply
        cannot reach it — there is no way round."""
        seq = _rand(300, 45) + _PERFECT_PROMOTER
        feats = [
            {"start": 300, "end": 300 + len(_PERFECT_PROMOTER), "strand": 1,
             "type": "promoter", "label": "P"},
            {"start": 10, "end": 200, "strand": 1, "type": "CDS",
             "label": "upstream gene"},
        ]
        res = sa._map_transcription(seq, feats, circular=False,
                                    include_predicted=False)
        g = _gene(res, "upstream gene")
        assert g["reachable_from"] == []
        assert g["silent"] is True

    def test_the_same_construct_answers_identically_when_rotated(self):
        seq, feats = _readthrough_construct()
        n = len(seq)
        base = sa._map_transcription(seq, feats, circular=True)
        rot = 777
        rseq = seq[rot:] + seq[:rot]
        rfeats = [{**f, "start": (f["start"] - rot) % n,
                   "end": (f["end"] - rot) % n} for f in feats]
        got = sa._map_transcription(rseq, rfeats, circular=True)
        b = _gene(base, "GeneX")
        g = _gene(got, "GeneX")
        assert b["silent"] == g["silent"]
        assert b["silent_in_host"] == g["silent_in_host"]
        assert sorted(r["distance_bp"] for r in b["reachable_from"]) == \
            sorted(r["distance_bp"] for r in g["reachable_from"])

    def test_no_cds_is_a_warning_not_a_crash(self):
        res = sa._map_transcription(_rand(500, 46), [], circular=True)
        assert res["ok"] and res["genes"] == []
        assert any("no CDS" in w for w in res["warnings"])

    def test_empty_sequence(self):
        res = sa._map_transcription("", [], circular=True)
        assert res["ok"] is False

    def test_coordinates_past_the_end_are_wrapped_but_reported(self):
        """HARDENING REGRESSION. Every position is taken `% total`, so a CDS
        annotated at bp 50,000 on a 1,000 bp plasmid became one at bp 0 — a
        confident answer about a gene that is not where the caller thinks it
        is. The modulo is what makes the wrap-aware arithmetic work; doing it
        silently is the part that was wrong."""
        seq = _rand(1000, 90)
        feats = [
            {"start": 0, "end": 30, "strand": 1, "type": "promoter",
             "label": "P"},
            {"start": 50_000, "end": 50_300, "strand": 1, "type": "CDS",
             "label": "ghost"},
        ]
        res = sa._map_transcription(seq, feats, circular=True,
                                    include_predicted=False)
        assert res["ok"]
        assert any("outside 0..1000" in w for w in res["warnings"])

    def test_in_range_coordinates_raise_no_such_warning(self):
        seq, feats = _readthrough_construct()
        res = sa._map_transcription(seq, feats, circular=True,
                                    include_predicted=False)
        assert not any("outside 0.." in w for w in res["warnings"])

    def test_a_fuzzy_coordinate_is_dropped_with_a_warning(self):
        """One unusable feature must not sink the whole map."""
        seq, feats = _readthrough_construct()
        feats = feats + [{"start": "<1", "end": 500, "strand": 1,
                          "type": "CDS", "label": "fuzzy"}]
        res = sa._map_transcription(seq, feats, circular=True)
        assert res["ok"]
        assert any("aren't whole numbers" in w for w in res["warnings"])


# ══════════════════════════════════════════════════════════════════════
#  Agent endpoints
# ══════════════════════════════════════════════════════════════════════

def _call(name, payload, app=None):
    return sc._state._AGENT_HANDLERS[name][0](app, payload)


class TestAgentEndpoints:

    def test_all_five_are_registered_read_only(self):
        """Read-only matters beyond tidiness: a `write=True` endpoint is
        refused outright by a `--read-only` attach."""
        for name in ("scan-promoters", "scan-terminators", "map-transcription",
                     "analyse-cds", "analyse-read-heterogeneity"):
            entry = sc._state._AGENT_HANDLERS.get(name)
            assert entry is not None, f"{name} not registered"
            assert entry[1] is False, f"{name} is registered as a write"

    def test_scan_promoters_over_a_bare_sequence(self):
        seq = _rand(200, 50) + _PERFECT_PROMOTER + _rand(200, 51)
        res = _call("scan-promoters", {"sequence": seq, "circular": False})
        assert res["ok"] and res["n_hits"] >= 1
        assert res["hits"][0]["score"] == pytest.approx(1.0)
        assert res["length"] == len(seq)

    def test_scan_promoters_limit_flags_truncation(self):
        res = _call("scan-promoters",
                    {"sequence": _rand(20_000, 52), "min_score": 0.0,
                     "limit": 5, "circular": False})
        assert len(res["hits"]) == 5
        assert res["truncated"] is True

    def test_scan_terminators_over_a_bare_sequence(self):
        seq = _rand(150, 53) + _TERMINATOR + _rand(150, 54)
        res = _call("scan-terminators", {"sequence": seq, "circular": False})
        assert res["ok"]
        assert any(h["start"] == 150 and h["stem_len"] == 10
                   for h in res["hits"])

    @pytest.mark.parametrize("body,code", [
        ({"sequence": "ACGT" * 50, "min_score": 5}, 400),
        ({"sequence": "ACGT" * 50, "min_score": "x"}, 400),
        ({"sequence": "ACGT" * 50, "limit": 0}, 400),
        ({}, 422),
    ])
    def test_scan_promoters_rejects_bad_input(self, body, code):
        res = _call("scan-promoters", body)
        assert isinstance(res, tuple) and res[1] == code

    @pytest.mark.parametrize("body,code", [
        ({"sequence": "ACGT" * 50, "min_stem": 1}, 400),
        ({"sequence": "ACGT" * 50, "min_u_tract": 0}, 400),
        ({"sequence": "ACGT" * 50, "min_u_tract": 99}, 400),
    ])
    def test_scan_terminators_rejects_bad_input(self, body, code):
        res = _call("scan-terminators", body)
        assert isinstance(res, tuple) and res[1] == code

    @pytest.mark.parametrize("key,value", [
        ("min_stem", reg._TERM_MAX_STEM + 1),
        ("min_stem", 30),
        ("min_u_tract", reg._TERM_U_TRACT_WINDOW + 1),
    ])
    def test_a_filter_that_can_never_match_is_refused_not_answered_empty(
            self, key, value):
        """HARDENING REGRESSION. `min_stem: 25` was accepted and returned
        `n_hits: 0` for every sequence, because the scanner only searches
        stems up to `_TERM_MAX_STEM`. "0 terminators" reads as "this construct
        has none", not as "that filter excluded every candidate there is" —
        the vacuous-pass shape `list-restriction-sites` already refuses for an
        over-long `min_length`."""
        res = _call("scan-terminators",
                    {"sequence": "ACGT" * 200, key: value, "circular": False})
        assert isinstance(res, tuple) and res[1] == 400
        assert "zero terminators" in res[0]["error"]

    @pytest.mark.parametrize("key,value", [
        ("min_stem", reg._TERM_MAX_STEM),
        ("min_u_tract", reg._TERM_U_TRACT_WINDOW),
    ])
    def test_the_boundary_value_itself_is_still_accepted(self, key, value):
        res = _call("scan-terminators",
                    {"sequence": "ACGT" * 200, key: value, "circular": False})
        assert not isinstance(res, tuple), res

    def test_map_transcription_refuses_a_bare_sequence(self):
        """It needs annotations. Silently scanning with no features would
        report every gene as absent rather than saying the input was wrong."""
        res = _call("map-transcription", {"sequence": "ACGT" * 100})
        assert isinstance(res, tuple) and res[1] == 400
        assert "annotations" in res[0]["error"]

    def test_map_transcription_with_no_plasmid_loaded(self):
        res = _call("map-transcription", {})
        assert isinstance(res, tuple) and res[1] == 422

    def test_unknown_library_name_is_404(self):
        for name in ("scan-promoters", "scan-terminators",
                     "map-transcription"):
            res = _call(name, {"name": "no such plasmid anywhere"})
            assert isinstance(res, tuple) and res[1] == 404, name


# ══════════════════════════════════════════════════════════════════════
#  Strand mirror (audit 2026-10-23, A1 / A5 / A6)
# ══════════════════════════════════════════════════════════════════════
#
# The reverse complement of a construct is the same molecule read from the
# other strand, so every answer about it must come back identical. Each check
# below failed on exactly one strand before the fix — a terminator touching
# the gene or the promoter was seen on one strand only, the minus-strand box
# starts pointed at each box's LAST base, and a predicted minus-strand
# promoter reported its leader one base short.

def _mirror(seq, feats):
    """The same construct drawn from the other strand."""
    n = len(seq)
    out = []
    for f in feats:
        g = dict(f)
        g["start"], g["end"] = n - f["end"], n - f["start"]
        g["strand"] = -f["strand"]
        out.append(g)
    return _rc(seq), out


def _p_to_g(res):
    """The P -> G reach, reduced to what must not depend on the strand."""
    gene = _gene(res, "G")
    hits = [r for r in gene["reachable_from"]
            if r["promoter"]["label"] == "P"]
    assert len(hits) == 1, hits
    r = hits[0]
    return (r["distance_bp"], r["readthrough"],
            [t["label"] for t in r["terminators_between"]])


# (terminator start, terminator end) relative to a layout of
#   pad 60 | promoter P [60, 90) | 120 bp | gene G [210, 510) | pad 60
_TERM_PLACEMENTS = {
    "touches the gene's 5' end": ((170, 210), ["T"]),
    "touches the promoter's 3' end": ((90, 130), ["T"]),
    "straddles the gene start": ((190, 230), ["T"]),
    "strictly between": ((130, 170), ["T"]),
    "inside the gene": ((300, 340), []),
    "straddles the promoter's fire point": ((80, 100), []),
}


class TestTranscriptionMapIsStrandSymmetric:

    @pytest.mark.parametrize("circular", [False, True])
    @pytest.mark.parametrize("where", sorted(_TERM_PLACEMENTS))
    def test_a_terminator_reads_the_same_from_either_strand(
            self, where, circular):
        (ts, te), expect = _TERM_PLACEMENTS[where]
        seq = _rand(570, 5100)
        feats = [
            {"start": 60, "end": 90, "strand": 1, "type": "promoter",
             "label": "P"},
            {"start": ts, "end": te, "strand": 1, "type": "terminator",
             "label": "T"},
            {"start": 210, "end": 510, "strand": 1, "type": "CDS",
             "label": "G"},
        ]
        fwd = _p_to_g(sa._map_transcription(seq, feats, circular=circular,
                                            include_predicted=False))
        mseq, mfeats = _mirror(seq, feats)
        rev = _p_to_g(sa._map_transcription(mseq, mfeats, circular=circular,
                                            include_predicted=False))
        assert fwd == rev, (where, fwd, rev)
        assert fwd[0] == 120 and fwd[2] == expect, (where, fwd)

    @pytest.mark.parametrize("trial", range(6))
    def test_a_predicted_promoter_reports_its_true_leader_on_both_strands(
            self, trial):
        """The leader is the transcribed stretch from the start site to the
        CDS: `gene_start - tss` bases on the plus strand. Its mirror was one
        base short — the minus-strand TSS is a base, and was used as the
        boundary on its far side."""
        r = random.Random(5200 + trial)
        leader = r.randint(10, 60)
        left = _rand(r.randint(80, 200), 5300 + trial)
        head = left + _PERFECT_PROMOTER
        tss = len(left) + len(_PERFECT_PROMOTER) - 6 + reg._TSS_OFFSET_FROM_MINUS10
        gap = _rand(tss + leader - len(head), 5400 + trial)
        seq = head + gap + _rand(300, 5500 + trial) + _rand(80, 5600 + trial)
        g0 = tss + leader
        feats = [{"start": g0, "end": g0 + 300, "strand": 1, "type": "CDS",
                  "label": "G"}]

        def leader_of(s, fs, strand):
            res = sa._map_transcription(s, fs, circular=False,
                                        include_predicted=True)
            planted = [x for x in _gene(res, "G")["reachable_from"]
                       if x["promoter"]["source"] == "predicted"
                       and x["promoter"]["strand"] == strand
                       and x["promoter"].get("score", 0) > 0.99]
            assert planted, "the planted promoter did not reach the gene"
            return min(x["distance_bp"] for x in planted)

        assert leader_of(seq, feats, 1) == leader
        mseq, mfeats = _mirror(seq, feats)
        assert leader_of(mseq, mfeats, -1) == leader


class TestPromoterBoxesReadBackOnTheirOwnStrand:
    """`minus35_start` / `minus10_start` are forward-strand coordinates like
    every other field: the six bases there, read on the hit's own strand, ARE
    the box. For a minus-strand hit they pointed at the box's last base."""

    @pytest.mark.parametrize("circular", [False, True])
    @pytest.mark.parametrize("trial", range(8))
    def test_box_starts_point_at_the_boxes(self, trial, circular):
        r = random.Random(5700 + trial)
        el = "TTGACA" + _rand(17, 5800 + trial) + "TATAAT"
        left = _rand(r.randint(40, 200), 5900 + trial)
        right = _rand(r.randint(40, 200), 6000 + trial)
        for strand, seq in ((1, left + el + right),
                            (-1, left + _rc(el) + right)):
            n = len(seq)
            hits = [h for h in reg._scan_promoters(seq, circular=circular)
                    if h["strand"] == strand and h["score"] > 0.99]
            assert hits, (strand, "the planted promoter was not found")
            for h in hits:
                for box in ("minus35", "minus10"):
                    p = h[box + "_start"]
                    six = (seq + seq)[p:p + 6] if circular else seq[p:p + 6]
                    read = six if strand == 1 else _rc(six)
                    assert read == h[box], (strand, box, p, n, read, h[box])


# ══════════════════════════════════════════════════════════════════════
#  Speed without a change of answer (audit 2026-10-23, A9)
# ══════════════════════════════════════════════════════════════════════
#
# The promoter scan rescored every window and built a components dict per
# candidate (32 s/Mbp); the terminator scan tried every stem/loop register with
# a fresh slice (44 s/Mbp) — `map-transcription` with predicted elements was a
# multi-minute call on a genome. Both are now exact rewrites. The digest below
# was computed with the PRE-rewrite scanners over a fixed corpus; a deliberate
# change to what the scanners report updates this ONE constant, an accidental
# one fails here.

_SCANNER_OUTPUT_DIGEST = (
    "8da0bbff1d8a04af7b6f84c5774ef9905b7dd16c249b91ed844110eef7814289")


def _scanner_corpus():
    """Ten molecules with planted promoters and terminators. A 64-bit LCG,
    not `random`, so the corpus is the same bytes on every Python."""
    prom = "TTGACA" + "ACGTGCATGCATGCATG" + "TATAAT"
    stem = "GCCGCCAGTT"
    term = stem + "CGAA" + _rc(stem) + "TTTTTTTT"
    x = 20261024

    def nxt(m):
        nonlocal x
        x = (x * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        return (x >> 33) % m
    out = []
    for k in range(10):
        n = 800 + nxt(3200)
        alph = ("ACGT", "GCGCAT", "ATATGC", "ACGTN")[k % 4]
        s = [alph[nxt(len(alph))] for _ in range(n)]
        for _ in range(4):
            el = (prom, term)[nxt(2)]
            el = el if nxt(2) else _rc(el)
            p = nxt(n - len(el))
            s[p:p + len(el)] = el
        out.append(("".join(s), k % 2 == 0))
    return out


class TestScannersAreFastAndUnchanged:

    def test_the_answers_are_the_pre_rewrite_answers(self):
        import hashlib
        import json
        h = hashlib.sha256()
        for seq, circ in _scanner_corpus():
            for res in (reg._scan_promoters(seq, circular=circ),
                        reg._scan_terminators(seq, circular=circ),
                        reg._scan_promoters(seq, circular=circ, min_score=0.5),
                        reg._scan_terminators(seq, circular=circ, min_stem=4,
                                              min_u_tract=3)):
                h.update(json.dumps(res, sort_keys=True, default=str).encode())
        assert h.hexdigest() == _SCANNER_OUTPUT_DIGEST

    def test_each_promoter_box_is_scored_once_per_hexamer(self, monkeypatch):
        calls = [0]
        real = reg._hexamer_score

        def counting(*a, **k):
            calls[0] += 1
            return real(*a, **k)
        monkeypatch.setattr(reg, "_hexamer_score", counting)
        reg._scan_promoters(_rand(20_000, 3), circular=True)
        # <= 4,096 distinct hexamers per box per strand; was 320,320.
        assert calls[0] <= 20_000, calls[0]

    def test_the_hairpin_search_does_not_retry_every_register(
            self, monkeypatch):
        class Counting(frozenset):
            n = 0

            def __contains__(self, item):
                Counting.n += 1
                return frozenset.__contains__(self, item)
        monkeypatch.setattr(reg, "_COMPLEMENT_PAIRS",
                            Counting(reg._COMPLEMENT_PAIRS))
        reg._scan_terminators(_rand(20_000, 3), circular=True)
        # 285k pair tests now; 959k when every register was sliced afresh.
        assert Counting.n <= 500_000, Counting.n
