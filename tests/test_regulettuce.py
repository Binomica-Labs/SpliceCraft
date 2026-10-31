"""Plant cis-regulatory annotation — the Regulettuce engine.

`splicecraft_regulettuce`: the curated plant cis-element catalogue, anchored
TATA boxes, and 3'-end-like regions.

Checked against things OUTSIDE the implementation, because reading a scanner
back to itself agrees with whatever it already does:

  * **REAL sequences** from this repo's own GenBank-verified preset catalogue.
    The CaMV 35S promoter's TATA box must land in the measured -40..-31 peak;
    the NOS terminator must flag as a plant 3'-end region; a 720 bp EGFP CDS
    must flag as nothing. That last pair is the whole discrimination claim.
  * **A BRUTE-FORCE reference** for motif finding — an independent nested-loop
    IUPAC matcher, so the regex/cache path cannot be wrong in agreement with
    itself.
  * **ROTATION and STRAND invariance** — a circular molecule has no origin, so
    the same construct must answer identically however it is drawn. This is
    what the `_signed_offset` wrap rule exists for: without it the one promoter
    straddling bp 0 is the one whose TATA box is missed.
  * **NON-REDUNDANCY** — one stretch of DNA is one record. A palindromic
    consensus, or a degenerate one matching on both strands, must not be
    emitted twice (`strand: 0` says so), the way `splicecraft_regulatory`'s
    terminator scan must not report one hairpin 22 times.
  * **THE REFUSALS** — the module's three measured non-claims are contract.
    `_scan_3prime_risk` must carry no cleavage-site field and
    `_scan_tata_boxes` must refuse to run unanchored; a later "improvement"
    that adds either is a regression, not a feature.
"""
from __future__ import annotations

import hashlib
import random

import pytest

import splicecraft as sc
import splicecraft_presets as presets
import splicecraft_regulettuce as rl


# Aggregate digest over the catalogue: name, consensus and both measured rates
# of every entry, in declaration order. A deliberate recalibration updates this
# constant in the same commit; an accidental edit to a consensus or a rate
# fails loudly here. `responds_to` prose and aliases are OUT of the digest.
_CATALOGUE_DIGEST = "5943c5c5d8ba0883c6714e20304bc15851b3d46a3f457cbf61d7fc77c7047853"
_CATALOGUE_COUNT = 34

_COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _catalogue_digest() -> str:
    h = hashlib.sha256()
    for e in rl._CIS_ELEMENTS:
        h.update(str(e["name"]).encode())
        h.update(b"\x00")
        h.update(str(e["consensus"]).encode())
        h.update(b"\x00")
        h.update(f"{float(e['promoter_rate']):.4f}".encode())
        h.update(b"\x00")
        h.update(f"{float(e['null_rate']):.4f}".encode())
        h.update(b"\x00")
    return h.hexdigest()


_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT",
          "S": "GC", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
          "H": "ACT", "V": "ACG", "N": "ACGT"}


def _brute_find(seq: str, consensus: str, circular: bool) -> "list[int]":
    """Every start index where `consensus` matches, by nested loops only.

    Independent of `_iupac_pattern`, the regex cache and `_iter_match_starts`,
    so a shared bug cannot make both agree. Overlapping matches included."""
    n = len(seq)
    width = len(consensus)
    hay = seq + (seq[:width - 1] if circular and width > 1 else "")
    out = []
    last = n if circular else n - width + 1
    for i in range(max(0, last)):
        for j, ch in enumerate(consensus):
            if i + j >= len(hay) or hay[i + j] not in _IUPAC[ch]:
                break
        else:
            out.append(i)
    return out


def _preset(name: str) -> dict:
    for f in presets._preset_features():
        if f["name"] == name:
            return f
    raise AssertionError(f"preset {name!r} missing from the catalogue")


def _rand_seq(n: int, seed: int, gc: float = 0.5) -> str:
    r = random.Random(seed)
    return "".join(r.choices("ACGT", weights=[(1 - gc) / 2, gc / 2,
                                              gc / 2, (1 - gc) / 2], k=n))


# ──────────────────────────────────────────────────────────────────────
class TestCatalogueIntegrity:
    def test_count_and_digest(self):
        assert len(rl._CIS_ELEMENTS) == _CATALOGUE_COUNT
        assert _catalogue_digest() == _CATALOGUE_DIGEST, (
            "cis-element catalogue changed. If deliberate, update "
            "_CATALOGUE_DIGEST in this file in the same commit.")

    def test_every_entry_has_the_required_fields(self):
        for row in rl._cis_element_catalogue():
            assert row["name"] and row["consensus"]
            assert row["responds_to"], f"{row['name']} has no documented role"
            assert 0.0 <= row["promoter_rate"] <= 1.0
            assert 0.0 < row["null_rate"] <= 1.0
            assert row["length"] == len(row["consensus"])

    def test_consensus_is_valid_iupac(self):
        for row in rl._cis_element_catalogue():
            assert set(row["consensus"]) <= set(_IUPAC), row["name"]

    def test_informativeness_is_derived_not_stored(self):
        """A stored label could disagree with the numbers beside it."""
        for e in rl._CIS_ELEMENTS:
            assert "informativeness" not in e, (
                f"{e['name']} stores a label; it must be derived from the "
                "measured rates so the two cannot drift apart")
        for row in rl._cis_element_catalogue():
            assert row["informativeness"] == rl._cis_informativeness(
                row["promoter_rate"], row["null_rate"])

    def test_palindromic_flag_matches_reverse_complement(self):
        for row in rl._cis_element_catalogue():
            assert row["palindromic"] == (_rc(row["consensus"])
                                          == row["consensus"]), row["name"]

    def test_ubiquitous_elements_are_flagged_not_hidden(self):
        """A motif in ~95% of random DNA must be KEPT and labelled.

        Dropping it would be worse than reporting it: people come looking for
        `GT1-motif` because other tools report it, and the useful answer is
        the null rate, not silence."""
        ubi = [r for r in rl._cis_element_catalogue()
               if r["informativeness"] == "ubiquitous"]
        assert ubi, "expected at least one ubiquitous element in the catalogue"
        for row in ubi:
            assert row["null_rate"] >= rl._CIS_UBIQUITOUS_NULL

    def test_aliases_resolve_and_are_not_duplicate_entries(self):
        """as-1 / TGA-element / TGACG-motif / CGTCA-motif are ONE site."""
        names = {r["name"] for r in rl._cis_element_catalogue()}
        for spelling in ("as-1", "TGA-element", "TGACG-motif", "CGTCA-motif"):
            row = rl._find_cis_element(spelling)
            assert row is not None, spelling
            assert row["name"] == "as-1"
        for dup in ("TGACG-motif", "CGTCA-motif"):
            assert dup not in names, (
                f"{dup} is a second name for as-1 (TGACG reverse-complements "
                "to CGTCA); a separate entry reports one site twice")

    def test_lookup_is_case_insensitive_and_rejects_unknown(self):
        assert rl._find_cis_element("g-BOX")["name"] == "G-box"
        assert rl._find_cis_element("not-an-element") is None
        assert rl._find_cis_element("") is None
        with pytest.raises(ValueError):
            rl._scan_cis_elements("ACGT" * 20, names=["not-an-element"])


# ──────────────────────────────────────────────────────────────────────
class TestAgainstBruteForce:
    @pytest.mark.parametrize("circular", [False, True])
    def test_every_catalogued_motif_matches_a_nested_loop_reference(self,
                                                                    circular):
        seq = _rand_seq(3000, seed=4)
        for row in rl._cis_element_catalogue():
            hits = rl._scan_cis_elements(seq, circular=circular,
                                         both_strands=False,
                                         names=[row["name"]],
                                         max_hits_per_element=10 ** 9)
            got = sorted(h["start"] for h in hits)
            want = sorted(_brute_find(seq, row["consensus"], circular))
            assert got == want, row["name"]

    def test_overlapping_tandem_matches_are_all_found(self):
        """`re.finditer` would drop the second of two overlapping hits."""
        # ATTAAT (Box 4) is palindromic; ATTAATTAAT holds two overlapping
        # occurrences, at 0 and 4.
        hits = rl._scan_cis_elements("GGG" + "ATTAATTAAT" + "GGG",
                                     circular=False, names=["Box 4"])
        assert sorted(h["start"] for h in hits) == [3, 7]


# ──────────────────────────────────────────────────────────────────────
class TestNonRedundancy:
    def test_palindromic_consensus_is_one_record_on_strand_zero(self):
        row = next(r for r in rl._cis_element_catalogue() if r["palindromic"])
        seq = "TTTT" + row["consensus"].replace("N", "A").replace(
            "W", "A").replace("Y", "C").replace("R", "A") + "TTTT"
        hits = [h for h in rl._scan_cis_elements(seq, circular=False,
                                                 names=[row["name"]])]
        assert len(hits) == 1, f"{row['name']} emitted {len(hits)} records"
        assert hits[0]["strand"] == 0, (
            "a palindromic consensus reads as the motif on both strands; "
            "strand 0 says one duplex, one record")

    def test_same_bases_never_emitted_twice(self):
        seq = _rand_seq(5000, seed=9)
        for name in ("as-1", "G-box", "W-box", "ABRE"):
            hits = rl._scan_cis_elements(seq, circular=True, names=[name],
                                         max_hits_per_element=10 ** 9)
            spans = [(h["start"], h["length"]) for h in hits]
            assert len(spans) == len(set(spans)), name


# ──────────────────────────────────────────────────────────────────────
class TestRotationAndStrandInvariance:
    def test_cis_elements_survive_every_rotation(self):
        seq = _rand_seq(900, seed=12)
        n = len(seq)
        base = rl._scan_cis_elements(seq, circular=True)
        want = sorted((h["name"], h["start"], h["strand"]) for h in base)
        assert want, "expected some hits to compare"
        for shift in (1, 7, 123, 449, 899):
            rot = seq[shift:] + seq[:shift]
            got = sorted(((h["name"], (h["start"] + shift) % n, h["strand"])
                          for h in rl._scan_cis_elements(rot, circular=True)))
            assert got == want, f"rotation by {shift} changed the hit set"

    def test_three_prime_regions_survive_rotation(self):
        # An AT-rich island inside a GC-rich molecule, placed across bp 0.
        seq = _rand_seq(600, seed=15, gc=0.72)
        island = _rand_seq(300, seed=16, gc=0.22)
        mol = island[150:] + seq + island[:150]
        n = len(mol)
        base = rl._scan_3prime_risk(mol, circular=True)
        assert base, "expected the AT-rich island to be flagged"
        assert len(base) == 1
        want_start, want_len = base[0]["start"], base[0]["length"]
        step = rl._3P_STEP
        for shift in (1, 50, 371, 740):
            rot = mol[shift:] + mol[:shift]
            got = rl._scan_3prime_risk(rot, circular=True)
            assert len(got) == 1, f"rotation by {shift}: {got}"
            start = (got[0]["start"] + shift) % n
            # Boundaries are quantised to the step grid, which is anchored on
            # bp 0 and therefore moves with the rotation. Detection and extent
            # must not.
            assert min((start - want_start) % n,
                       (want_start - start) % n) < step, shift
            assert abs(got[0]["length"] - want_len) <= step, shift

    def test_step_one_makes_boundaries_rotation_exact(self):
        """The quantisation above is the grid, not a coordinate bug.

        At step 1 there is no grid, and the region must land on exactly the
        same bases however the molecule is drawn."""
        seq = _rand_seq(400, seed=15, gc=0.75)
        island = _rand_seq(260, seed=16, gc=0.18)
        mol = island[130:] + seq + island[:130]
        n = len(mol)
        base = rl._scan_3prime_risk(mol, circular=True, step=1)
        assert base
        want = sorted((h["start"], h["length"]) for h in base)
        for shift in (1, 37, 200, 519):
            rot = mol[shift:] + mol[:shift]
            got = sorted(((h["start"] + shift) % n, h["length"])
                         for h in rl._scan_3prime_risk(rot, circular=True,
                                                       step=1))
            assert got == want, f"rotation by {shift} moved a step-1 region"

    def test_wrapping_region_is_one_record_not_two_stubs(self):
        seq = _rand_seq(600, seed=18, gc=0.72)
        island = _rand_seq(320, seed=19, gc=0.20)
        mol = island[160:] + seq + island[:160]
        hits = rl._scan_3prime_risk(mol, circular=True)
        wrapping = [h for h in hits if h["end"] < h["start"]]
        assert len(wrapping) == 1, (
            f"an origin-spanning AT-rich region must be ONE arc, got {hits}")

    def test_reverse_strand_hit_reports_forward_coordinates(self):
        """Sacred #2: a reverse hit at p is stored at p, not mirrored."""
        motif = "TAACAAA"                      # GARE, not palindromic
        seq = "GC" * 30 + _rc(motif) + "GC" * 30
        hits = rl._scan_cis_elements(seq, circular=False, names=["GARE"])
        assert len(hits) == 1
        h = hits[0]
        assert h["strand"] == -1
        assert h["start"] == 60
        assert seq[h["start"]:h["start"] + len(motif)] == _rc(motif)


# ──────────────────────────────────────────────────────────────────────
class TestRealSequences:
    def test_camv_35s_tata_box_lands_in_the_measured_peak(self):
        """The 35S promoter's TATA box is real and documented.

        Anchored at the promoter's 3' end (where its TSS sits by construction),
        it must fall inside -40..-31 — the peak bin measured over 3,089
        verified plant TSSs. This ties the scanner, the anchor convention and
        the calibration together on a sequence extracted from a public GenBank
        record rather than typed here."""
        seq = _preset("CaMV 35S promoter")["sequence"]
        hits = rl._scan_tata_boxes(
            seq, anchors=[{"position": len(seq) - 1, "strand": 1,
                           "label": "35S TSS"}], circular=False)
        assert hits, "no TATA box found in the CaMV 35S promoter"
        peak = [h for h in hits if h["tier"] == "peak"]
        assert peak, f"35S TATA not in the measured peak: {hits}"
        assert rl._TATA_PEAK[0] <= peak[0]["offset"] <= rl._TATA_PEAK[1]

    def test_nos_terminator_reads_as_a_plant_three_prime_region(self):
        seq = _preset("NOS terminator")["sequence"]
        hits = rl._scan_3prime_risk(seq, circular=False)
        assert hits, "the NOS terminator is a real plant 3' element"
        assert hits[0]["at_fraction"] > 0.6

    def test_a_coding_sequence_is_not_flagged(self):
        """The discrimination claim, on a real CDS.

        AUC 0.895 against coding sequence is only meaningful if an actual
        720 bp reporter CDS comes back clean at the shipped threshold."""
        seq = _preset("EGFP")["sequence"]
        assert rl._scan_3prime_risk(seq, circular=False) == []

    def test_bacterial_terminator_presets_do_not_dominate(self):
        """A sanity check that the 3' scanner is not just "AT-rich".

        The rrnB T1 terminator is a GC-stem bacterial element; it should not
        outscore the NOS terminator, a genuine plant 3' end."""
        nos = rl._scan_3prime_risk(_preset("NOS terminator")["sequence"],
                                   circular=False)
        rrnb = rl._scan_3prime_risk(_preset("rrnB T1 terminator")["sequence"],
                                    circular=False)
        best_nos = max((h["score"] for h in nos), default=0.0)
        best_rrnb = max((h["score"] for h in rrnb), default=0.0)
        assert best_nos > best_rrnb


# ──────────────────────────────────────────────────────────────────────
class TestTataAnchoring:
    def test_refuses_to_run_without_an_anchor(self):
        """An unanchored TATAWA scan is chance hits, one per 429 bp."""
        with pytest.raises(ValueError, match="anchor"):
            rl._scan_tata_boxes(_rand_seq(5000, seed=21), anchors=[])
        with pytest.raises(ValueError, match="anchor"):
            rl._scan_tata_boxes(_rand_seq(5000, seed=21), anchors=None)

    def test_reverse_anchor_searches_higher_coordinates(self):
        """A reverse-strand gene's promoter lies at HIGHER forward coords.

        Defaulting it to forward would search the transcript body."""
        tata = "TATAAA"
        # Gene on the reverse strand with its +1 at index 100; its TATA box
        # 33 bp upstream sits at forward 133..138, reading on the complement.
        seq = list(_rand_seq(300, seed=24, gc=0.7))
        seq[133:139] = list(_rc(tata))
        seq = "".join(seq)
        hits = rl._scan_tata_boxes(
            seq, anchors=[{"position": 100, "strand": -1}], circular=False)
        offs = [h["offset"] for h in hits]
        assert -38 in offs or -33 in offs, (
            f"expected the upstream reverse TATA, got {hits}")
        got = next(h for h in hits if h["start"] == 133)
        assert got["strand"] == -1
        assert got["offset"] == -(138 - 100)

    def test_offset_takes_the_short_way_round_the_origin(self):
        """A promoter straddling bp 0 must not report a huge positive offset."""
        tata = "TATAAA"
        n = 400
        seq = list(_rand_seq(n, seed=27, gc=0.7))
        seq[n - 36:n - 30] = list(tata)          # 6 nt ending at n-31
        seq = "".join(seq)
        # Start site just past the origin, so the TATA is upstream ACROSS bp 0.
        hits = rl._scan_tata_boxes(seq, anchors=[5], circular=True)
        offs = [h["offset"] for h in hits]
        assert any(-45 <= o <= -30 for o in offs), (
            f"origin-spanning TATA lost; offsets were {offs}")
        assert all(o < 0 for o in offs), offs

    def test_rejects_a_bad_anchor_or_window(self):
        with pytest.raises(ValueError):
            rl._scan_tata_boxes("ACGT" * 50, anchors=[{"strand": 1}])
        with pytest.raises(ValueError):
            rl._scan_tata_boxes("ACGT" * 50, anchors=[{"position": 1,
                                                       "strand": 7}])
        with pytest.raises(ValueError):
            rl._scan_tata_boxes("ACGT" * 50, anchors=[10], window=(-5, -60))


# ──────────────────────────────────────────────────────────────────────
class TestChanceExpectation:
    def test_expected_count_tracks_a_real_count_on_random_dna(self):
        """The null must be believable on the sequence actually queried.

        A uniform 25% null would make every AT-rich motif look enriched in
        AT-rich DNA, which is most plant DNA and most backbones."""
        seq = _rand_seq(40000, seed=30, gc=0.38)
        rows = {r["name"]: r for r in
                rl._cis_element_summary(seq, circular=True)}
        for name in ("G-box", "Box 4", "GARE", "as-1", "DRE / CRT"):
            r = rows[name]
            assert r["expected"] > 0
            assert 0.5 <= r["ratio"] <= 2.0, (
                f"{name}: observed {r['observed']} vs expected "
                f"{r['expected']} on random DNA -- the null is off")

    def test_palindromic_motif_is_not_double_counted_in_the_null(self):
        seq = _rand_seq(20000, seed=33)
        rows = {r["name"]: r for r in rl._cis_element_summary(seq)}
        pal = next(r for r in rl._cis_element_catalogue() if r["palindromic"])
        r = rows[pal["name"]]
        assert 0.5 <= r["ratio"] <= 2.0, (
            f"{pal['name']} is palindromic; counting its reverse complement "
            f"separately would double the expectation: {r}")

    def test_summary_counts_every_hit_even_past_the_hit_cap(self):
        """The hit list is capped for display; the arithmetic must not be.

        A capped `observed` would understate the count and make a ubiquitous
        element look enriched."""
        seq = _rand_seq(30000, seed=36)
        ubi = next(r for r in rl._cis_element_catalogue()
                   if r["informativeness"] == "ubiquitous")
        capped = rl._scan_cis_elements(seq, circular=True,
                                       names=[ubi["name"]],
                                       max_hits_per_element=5)
        rows = {r["name"]: r for r in rl._cis_element_summary(seq)}
        observed = rows[ubi["name"]]["observed"]
        assert len(capped) <= 10, "cap not applied to the hit list"
        assert observed > 50, f"expected many hits, got {observed}"
        assert observed > len(capped) * 4, (
            f"summary observed ({observed}) looks capped like the hit list "
            f"({len(capped)})")


# ──────────────────────────────────────────────────────────────────────
class TestTheRefusalsAreContract:
    """The module's measured non-claims. Adding either back is a regression."""

    def test_three_prime_never_reports_a_cleavage_site(self):
        seq = _preset("NOS terminator")["sequence"]
        for h in rl._scan_3prime_risk(seq, circular=False):
            assert "cleavage_site" not in h, (
                "a position-specific model located real plant poly(A) sites "
                "BELOW chance (27.2% within +-20 nt vs 29.1% expected). A "
                "null field would still invite a feature arrow.")
            assert "polya_site" not in h and "position" not in h
            assert h["length"] >= rl._3P_WINDOW

    def test_no_promoter_caller_or_y_patch_scanner_exists(self):
        for banned in ("_scan_plant_promoters", "_score_plant_promoter",
                       "_scan_y_patch", "_scan_polya_sites",
                       "_predict_polya_site"):
            assert not hasattr(rl, banned), (
                f"{banned} was measured and refused -- see the module "
                "docstring for the number that killed it")

    def test_calibration_publishes_the_refusals_and_their_numbers(self):
        cal = rl._regulettuce_calibration()
        assert cal["three_prime"]["predicts_cleavage_site"] is False
        assert cal["tata"]["requires_anchor"] is True
        whats = [d["what"] for d in cal["not_modelled"]]
        assert len(whats) == 3
        for d in cal["not_modelled"]:
            assert d["why"].strip(), d["what"]
            assert any(ch.isdigit() for ch in d["why"]), (
                f"{d['what']}: a refusal must cite the measurement that "
                "produced it, not an opinion")

    def test_calibration_reports_its_reference_set_sizes(self):
        cal = rl._regulettuce_calibration()
        assert cal["reference_sets"]["plant_promoters"]["n"] == 3089
        assert cal["reference_sets"]["polya_sites"]["n"] == 600
        assert "5'UTR" in cal["reference_sets"]["plant_promoters"]["source"]
        assert cal["cis_elements"]["total"] == _CATALOGUE_COUNT


# ──────────────────────────────────────────────────────────────────────
class TestEdgeCases:
    @pytest.mark.parametrize("seq", ["", "N" * 50, "ACGT", "   "])
    def test_degenerate_input_returns_empty_not_raises(self, seq):
        assert rl._scan_cis_elements(seq) == []
        assert rl._scan_3prime_risk(seq) == []

    def test_ambiguity_codes_do_not_match(self):
        """An element called on a base the record does not know is a guess."""
        assert rl._scan_cis_elements("TTT" + "NNNNNN" + "TTT",
                                     circular=False) == []

    def test_sequence_shorter_than_the_window_is_not_flagged(self):
        assert rl._scan_3prime_risk("AT" * 40, circular=False) == []

    def test_lowercase_and_whitespace_are_handled(self):
        a = rl._scan_cis_elements("gc" * 20 + "cacgtg" + "gc" * 20,
                                  circular=False, names=["G-box"])
        assert len(a) == 1 and a[0]["strand"] == 0

    def test_hit_cap_is_per_element(self):
        seq = _rand_seq(20000, seed=39)
        hits = rl._scan_cis_elements(seq, circular=True,
                                     max_hits_per_element=5)
        per = {}
        for h in hits:
            per[h["name"]] = per.get(h["name"], 0) + 1
        assert per, "expected hits"
        # Both strands of a non-palindromic motif each get the cap.
        assert max(per.values()) <= 10
        assert len(per) > 5, "one ubiquitous element must not starve the rest"


# ──────────────────────────────────────────────────────────────────────
class TestHubReExport:
    def test_names_resolve_through_the_hub(self):
        for name in ("_scan_cis_elements", "_cis_element_summary",
                     "_cis_element_catalogue", "_find_cis_element",
                     "_scan_tata_boxes", "_scan_3prime_risk",
                     "_regulettuce_calibration", "_TATA_CONSENSUS",
                     "_3P_MIN_SCORE"):
            assert hasattr(sc, name), f"sc.{name} does not resolve"
        assert sc._scan_cis_elements is rl._scan_cis_elements


# ══════════════════════════════════════════════════════════════════════
#  Agent endpoints
# ══════════════════════════════════════════════════════════════════════

def _call(name, payload, app=None):
    return sc._state._AGENT_HANDLERS[name][0](app, payload)


class TestAgentEndpoints:

    def test_all_five_are_registered_read_only(self):
        """A `write=True` endpoint is refused outright by a `--read-only`
        attach, which is how the GUI-is-open case is served."""
        for name in ("list-cis-elements", "scan-cis-elements",
                     "scan-tata-boxes", "scan-3prime-risk",
                     "regulettuce-calibration"):
            entry = sc._state._AGENT_HANDLERS.get(name)
            assert entry is not None, f"{name} not registered"
            assert entry[1] is False, f"{name} is registered as a write"

    def test_list_cis_elements(self):
        res = _call("list-cis-elements", {})
        assert res["ok"] and res["n_elements"] == _CATALOGUE_COUNT
        assert res["by_informativeness"]["ubiquitous"] >= 1
        row = next(r for r in res["elements"] if r["name"] == "GT1-motif")
        assert row["informativeness"] == "ubiquitous"
        assert row["null_rate"] > 0.9, (
            "the point of the field: a caller must be able to see that this "
            "element is in most random DNA before reporting a hit")

    def test_scan_cis_elements_returns_hits_and_the_chance_summary(self):
        seq = _preset("CaMV 35S promoter")["sequence"]
        res = _call("scan-cis-elements", {"sequence": seq, "circular": False})
        assert res["ok"]
        assert res["hits"] and res["summary"]
        names = {h["name"] for h in res["hits"]}
        assert "ABRE" in names
        s = {r["name"]: r for r in res["summary"]}
        assert s["ABRE"]["observed"] >= 1
        assert s["ABRE"]["expected"] > 0
        assert s["ABRE"]["verdict"]

    def test_summary_is_not_capped_by_limit(self):
        """`limit` trims the display list; the arithmetic must be complete."""
        seq = _rand_seq(20000, seed=42)
        res = _call("scan-cis-elements",
                    {"sequence": seq, "circular": True, "limit": 3})
        assert len(res["hits"]) == 3 and res["truncated"] is True
        total = sum(r["observed"] for r in res["summary"])
        assert total == res["n_hits"] > 3

    def test_scan_cis_elements_name_filter_and_unknown_name(self):
        seq = _preset("CaMV 35S promoter")["sequence"]
        ok = _call("scan-cis-elements",
                   {"sequence": seq, "circular": False, "names": ["ABRE"]})
        assert ok["ok"] and {h["name"] for h in ok["hits"]} == {"ABRE"}
        bad = _call("scan-cis-elements",
                    {"sequence": seq, "names": ["nope"]})
        assert bad[1] == 400 and "list-cis-elements" in bad[0]["error"]

    def test_scan_tata_boxes_requires_anchors(self):
        seq = _rand_seq(4000, seed=45)
        res = _call("scan-tata-boxes", {"sequence": seq, "circular": True})
        assert res[1] == 400
        assert "429" in res[0]["error"], (
            "the refusal must carry the number that justifies it")

    def test_scan_tata_boxes_accepts_int_and_object_anchors(self):
        seq = _preset("CaMV 35S promoter")["sequence"]
        a = _call("scan-tata-boxes", {"sequence": seq, "circular": False,
                                      "anchors": [len(seq) - 1]})
        b = _call("scan-tata-boxes", {"sequence": seq, "circular": False,
                                      "anchors": {"position": len(seq) - 1,
                                                  "strand": 1,
                                                  "label": "35S TSS"}})
        assert a["ok"] and b["ok"]
        assert [h["offset"] for h in a["hits"]] == [h["offset"]
                                                    for h in b["hits"]]
        assert b["hits"][0]["anchor_label"] == "35S TSS"
        assert any(h["tier"] == "peak" for h in b["hits"])

    def test_scan_3prime_risk_declares_its_non_claim_in_the_payload(self):
        seq = _preset("NOS terminator")["sequence"]
        res = _call("scan-3prime-risk", {"sequence": seq, "circular": False})
        assert res["ok"] and res["hits"]
        assert res["predicts_cleavage_site"] is False, (
            "the refusal belongs in the response, not only in the docs -- an "
            "agent reads the payload")
        for h in res["hits"]:
            assert "cleavage_site" not in h

    def test_scan_3prime_risk_refuses_an_unbounded_scan(self):
        res = _call("scan-3prime-risk",
                    {"sequence": "AT" * 300_000, "step": 1, "circular": True})
        assert res[1] == 400
        assert "50,000,000" in res[0]["error"]

    def test_regulettuce_calibration_publishes_the_refusals(self):
        res = _call("regulettuce-calibration", {})
        assert res["ok"]
        assert len(res["not_modelled"]) == 3
        assert res["three_prime"]["predicts_cleavage_site"] is False
        assert res["reference_sets"]["plant_promoters"]["n"] == 3089

    @pytest.mark.parametrize("name,body,code", [
        ("scan-cis-elements", {"sequence": "ACGT" * 50, "limit": 0}, 400),
        ("scan-cis-elements", {"sequence": "ACGT" * 50,
                               "min_enrichment": -1}, 400),
        ("scan-cis-elements", {"sequence": "ACGT" * 50, "names": []}, 400),
        ("scan-cis-elements", {}, 422),
        ("scan-tata-boxes", {"sequence": "ACGT" * 50, "anchors": []}, 400),
        ("scan-tata-boxes", {"sequence": "ACGT" * 50, "anchors": [1],
                             "window": [-5, -60]}, 400),
        ("scan-tata-boxes", {"sequence": "ACGT" * 50,
                             "anchors": [{"strand": 1}]}, 400),
        ("scan-tata-boxes", {"sequence": "ACGT" * 50,
                             "anchors": [{"position": 1,
                                          "strand": 9}]}, 400),
        ("scan-tata-boxes", {}, 422),
        ("scan-3prime-risk", {"sequence": "ACGT" * 50, "window": 5}, 400),
        ("scan-3prime-risk", {"sequence": "ACGT" * 50, "step": 0}, 400),
        ("scan-3prime-risk", {"sequence": "ACGT" * 50, "min_score": 7}, 400),
        ("scan-3prime-risk", {}, 422),
    ])
    def test_rejects_bad_input(self, name, body, code):
        res = _call(name, body)
        assert isinstance(res, tuple), f"{name} accepted {body!r}"
        assert res[1] == code, f"{name} {body!r} -> {res}"
        assert res[0]["error"]

    def test_unknown_keys_are_reported_not_silently_ignored(self):
        res = _call("scan-cis-elements",
                    {"sequence": "ACGT" * 60, "circular": False,
                     "nonsense": 1})
        assert "nonsense" in res["ignored"]


# ══════════════════════════════════════════════════════════════════════
#  UI: the modal, the anchor reader, and the annotate path
# ══════════════════════════════════════════════════════════════════════

def _plant_record():
    """A small circular record with a real plant promoter and a CDS."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    prom = _preset("CaMV 35S promoter")["sequence"]
    gfp = _preset("EGFP")["sequence"]
    term = _preset("NOS terminator")["sequence"]
    seq = prom + gfp + term
    rec = SeqRecord(Seq(seq), id="pPlant", name="pPlant")
    rec.annotations["topology"] = "circular"
    rec.features = [
        SeqFeature(FeatureLocation(0, len(prom), strand=1), type="promoter",
                   qualifiers={"label": ["35S"]}),
        SeqFeature(FeatureLocation(len(prom), len(prom) + len(gfp), strand=1),
                   type="CDS", qualifiers={"label": ["EGFP"]}),
    ]
    return rec


class TestAnchorReader:
    """`_regulettuce_anchors` runs on the real record, not a stub.

    It is the one place a typo'd helper name would only surface at runtime —
    which is exactly what happened while this was being built."""

    def test_reads_cds_starts_off_a_real_record(self):
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        rec = _plant_record()
        anchors = sc.PlasmidApp._regulettuce_anchors(app, rec)
        assert anchors, "a record with a CDS must yield an anchor"
        cds = next(a for a in anchors if "EGFP" in a["label"])
        assert cds["strand"] == 1
        assert cds["position"] == len(_preset("CaMV 35S promoter")["sequence"])

    def test_reverse_strand_feature_anchors_at_its_highest_coordinate(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("ACGT" * 200), id="r", name="r")
        rec.features = [SeqFeature(FeatureLocation(100, 400, strand=-1),
                                   type="CDS",
                                   qualifiers={"label": ["rev"]})]
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        a = sc.PlasmidApp._regulettuce_anchors(app, rec)[0]
        assert a["strand"] == -1
        assert a["position"] == 399, (
            "a reverse gene's start is its HIGHEST forward coordinate; "
            "anchoring at the low end searches the transcript body")

    def test_empty_record_yields_no_anchors(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        rec = SeqRecord(Seq(""), id="e", name="e")
        assert sc.PlasmidApp._regulettuce_anchors(app, rec) == []


async def _settle(pilot):
    """Wait for the scan worker, then let the UI apply its result.

    The scan runs on a thread (a 2 Mbp record is 6 s of work and the app opens
    chromosome-scale records), so a bare `pause()` is a race dressed as a
    test — it would pass on a fast machine and flake on a loaded one."""
    await pilot.pause()
    await pilot.app.workers.wait_for_complete()
    await pilot.pause()


class TestModal:
    @pytest.mark.asyncio
    async def test_scans_on_mount_and_switches_mode(self):
        """Drives the real modal: mount, read the table, change mode."""
        rec = _plant_record()
        modal = sc.RegulettuceModal(
            str(rec.seq), circular=True, label="pPlant",
            anchors=[{"position": 346, "strand": 1, "label": "EGFP"}])

        class _Host(sc.App):
            def on_mount(self) -> None:
                self.push_screen(modal)

        async with _Host().run_test() as pilot:
            await _settle(pilot)
            from textual.widgets import DataTable, Static
            table = modal.query_one("#reg-table", DataTable)
            assert modal._mode == "cis"
            assert table.row_count == len(modal._hits) > 0
            status = modal.query_one("#reg-status", Static)
            assert "chance" in str(status.render()).lower(), (
                "the chance count must be on screen beside the hits")
            # 3'-end mode on a record ending in the NOS terminator.
            modal._mode = "3prime"
            modal._rebuild_columns()
            modal._run_scan()
            await _settle(pilot)
            assert modal._hits, "the NOS terminator should flag"
            assert all("cleavage_site" not in h for h in modal._hits)
            assert table.row_count == len(modal._hits)

    @pytest.mark.asyncio
    async def test_tata_mode_without_anchors_says_so_rather_than_showing_none(
            self):
        """An empty table reads as "no TATA boxes"; the truth is "not asked"."""
        modal = sc.RegulettuceModal(_rand_seq(3000, seed=51),
                                    circular=True, label="x", anchors=[])

        class _Host(sc.App):
            def on_mount(self) -> None:
                self.push_screen(modal)

        async with _Host().run_test() as pilot:
            await _settle(pilot)
            modal._mode = "tata"
            modal._rebuild_columns()
            modal._run_scan()
            await _settle(pilot)
            from textual.widgets import Static
            text = str(modal.query_one("#reg-status", Static).render())
            assert modal._hits == []
            assert "429" in text, (
                "the refusal must explain itself with the measured rate")

    def test_dismiss_payload_shape(self):
        modal = sc.RegulettuceModal("ACGT" * 100, circular=False, label="x")
        modal._hits = [{"start": 1, "end": 7, "strand": 1, "name": "G-box"}]
        sent = {}
        modal.dismiss = lambda v: sent.update(v or {})
        modal._on_annotate_all(None)
        assert sent["action"] == "annotate"
        assert sent["kind"] == "cis"
        assert len(sent["hits"]) == 1


class TestAnnotateResult:
    def _app_with(self, rec):
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        app._current_record = rec
        calls = {"transfers": None, "dirty": 0, "notes": []}

        def _apply(transfers, target):
            calls["transfers"] = transfers
            # The real `_apply_annotation_transfers` marks the record dirty
            # and raises its own success toast; the stub does both so a
            # DUPLICATE from `_regulettuce_result` is visible as a second
            # count rather than hidden by the stub swallowing the first.
            calls["dirty"] += 1
            calls["notes"].append(("applied",))
        # Stub only I/O: the transfer builder under test stays real.
        app._apply_annotation_transfers = _apply
        app._mark_dirty = lambda: calls.__setitem__("dirty",
                                                    calls["dirty"] + 1)
        app._notify_success = lambda *a, **k: calls["notes"].append(a)
        app.notify = lambda *a, **k: calls["notes"].append(a)
        return app, calls

    def test_builds_one_transfer_per_hit_in_a_single_call(self):
        """One call to `_apply_annotation_transfers` = one undo for all N."""
        rec = _plant_record()
        app, calls = self._app_with(rec)
        hits = sc._scan_cis_elements(str(rec.seq), circular=True,
                                     include_uninformative=False)[:5]
        assert len(hits) >= 2
        sc.PlasmidApp._regulettuce_result(
            app, {"action": "annotate", "kind": "cis", "hits": hits})
        assert calls["transfers"] is not None
        assert len(calls["transfers"]) == len(hits)
        assert calls["dirty"] == 1, (
            "exactly one dirty mark -- `_apply_annotation_transfers` already "
            "makes it, so a second one here means duplicated bookkeeping")
        assert len(calls["notes"]) == 2, (
            "one toast from the transfer helper plus the one this flow adds "
            "about the single Ctrl+Z; a third means a duplicate notification")

    def test_note_carries_the_measurement_not_just_the_name(self):
        rec = _plant_record()
        app, calls = self._app_with(rec)
        hits = [h for h in sc._scan_cis_elements(str(rec.seq), circular=True)
                if h["informativeness"] == "enriched"][:1]
        sc.PlasmidApp._regulettuce_result(
            app, {"action": "annotate", "kind": "cis", "hits": hits})
        note = calls["transfers"][0]["qualifiers"]["note"]
        assert "enriched" in note
        assert "not measured occupancy" in note, (
            "an element arrow with no caveat in the note is the overclaim "
            "this subsystem exists to avoid -- and the note is what survives "
            "into the exported GenBank file")

    def test_three_prime_feature_key_does_not_claim_a_polya_signal(self):
        rec = _plant_record()
        app, calls = self._app_with(rec)
        hits = sc._scan_3prime_risk(str(rec.seq), circular=True)[:1]
        assert hits
        sc.PlasmidApp._regulettuce_result(
            app, {"action": "annotate", "kind": "3prime", "hits": hits})
        tr = calls["transfers"][0]
        assert tr["type"] == "misc_feature", (
            "polyA_signal would assert a site the scanner measured itself "
            "unable to call")
        assert "NO cleavage site is predicted" in tr["qualifiers"]["note"]

    def test_palindromic_strand_zero_is_recorded_on_the_plus_strand(self):
        """GenBank has no "both strands"; 0 must not become strand 0."""
        rec = _plant_record()
        app, calls = self._app_with(rec)
        sc.PlasmidApp._regulettuce_result(app, {
            "action": "annotate", "kind": "cis",
            "hits": [{"start": 10, "end": 16, "strand": 0, "name": "G-box",
                      "consensus": "CACGTG", "responds_to": "light",
                      "enrichment": 1.9, "informativeness": "enriched"}]})
        assert calls["transfers"][0]["target_strand"] == 1

    def test_ignores_a_cancel_or_malformed_dismissal(self):
        rec = _plant_record()
        app, calls = self._app_with(rec)
        for bad in (None, {}, {"action": "nope"}, "x", 7):
            sc.PlasmidApp._regulettuce_result(app, bad)
        assert calls["transfers"] is None and calls["dirty"] == 0

    def test_empty_hit_list_notifies_instead_of_annotating_nothing(self):
        rec = _plant_record()
        app, calls = self._app_with(rec)
        sc.PlasmidApp._regulettuce_result(
            app, {"action": "annotate", "kind": "cis", "hits": []})
        assert calls["transfers"] is None
        assert calls["notes"], "a no-op must say so"


class TestBothHalvesAgree:
    """The scan and the summary must answer the same way about the same input.

    They did not: the summary resolved an unknown element name to ``None``,
    built a filter set containing it, and returned an EMPTY summary while the
    scan raised ``ValueError``. Pyright found it; this pins it."""

    def test_unknown_name_raises_from_both(self):
        seq = "ACGT" * 100
        with pytest.raises(ValueError, match="unknown cis-element"):
            rl._scan_cis_elements(seq, names=["nope"])
        with pytest.raises(ValueError, match="unknown cis-element"):
            rl._cis_element_summary(seq, names=["nope"])

    def test_alias_resolves_the_same_in_both(self):
        seq = _rand_seq(4000, seed=57)
        hits = rl._scan_cis_elements(seq, names=["TGACG-motif"])
        rows = rl._cis_element_summary(seq, names=["TGACG-motif"])
        assert {h["name"] for h in hits} <= {"as-1"}
        assert [r["name"] for r in rows] == ["as-1"]
        assert rows[0]["observed"] == len(hits)

    def test_name_filter_agrees_on_counts(self):
        seq = _rand_seq(6000, seed=58)
        for name in ("G-box", "Box 4", "as-1", "GT1-motif"):
            hits = rl._scan_cis_elements(seq, names=[name],
                                         max_hits_per_element=10 ** 9)
            rows = rl._cis_element_summary(seq, names=[name])
            assert len(rows) == 1
            assert rows[0]["observed"] == len(hits), name


class TestAnchorEdgeCases:
    """`_feat_bounds` has two traps the anchor reader must not fall into."""

    def test_spliced_cds_anchors_at_its_first_exon_on_a_linear_record(self):
        """On a LINEAR record a two-part location is SPLICED, not a wrap.

        Letting `_feat_bounds` assume the wrap re-encodes it as
        ``(exon2_start, exon1_end)``, which would anchor the gene at its last
        exon and measure every TATA offset from the wrong place."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import (SeqFeature, FeatureLocation,
                                    CompoundLocation)
        rec = SeqRecord(Seq("ACGT" * 100), id="lin", name="lin")
        rec.annotations["topology"] = "linear"
        rec.features = [SeqFeature(
            CompoundLocation([FeatureLocation(40, 90, strand=1),
                              FeatureLocation(200, 300, strand=1)]),
            type="CDS", qualifiers={"label": ["spliced"]})]
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        a = sc.PlasmidApp._regulettuce_anchors(app, rec)[0]
        assert a["position"] == 40, (
            "a spliced CDS on a linear record starts at its FIRST exon")

    def test_feature_with_no_integer_start_is_skipped_not_fatal(self):
        """`_feat_bounds` returns None for UnknownPosition — documented, not
        an error. The reader must skip it and keep going."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import (SeqFeature, FeatureLocation,
                                    UnknownPosition)
        rec = SeqRecord(Seq("ACGT" * 100), id="u", name="u")
        rec.features = [
            SeqFeature(FeatureLocation(UnknownPosition(), 50, strand=1),
                       type="CDS", qualifiers={"label": ["fuzzy"]}),
            SeqFeature(FeatureLocation(100, 200, strand=1),
                       type="CDS", qualifiers={"label": ["solid"]}),
        ]
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        got = sc.PlasmidApp._regulettuce_anchors(app, rec)
        assert [a["label"] for a in got] == ["solid"]

    def test_anchor_count_is_capped(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("ACGT" * 2000), id="many", name="many")
        rec.features = [
            SeqFeature(FeatureLocation(i * 10, i * 10 + 9, strand=1),
                       type="CDS", qualifiers={"label": [f"c{i}"]})
            for i in range(700)]
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        assert len(sc.PlasmidApp._regulettuce_anchors(app, rec)) == 500


class TestScanWorker:
    """The scan runs off the UI thread, and a stale result must not land.

    Measured on the UI thread: 1.1 s for the cis scan and 6.0 s for the 3'-end
    scan on a 2 Mbp record. SpliceCraft opens chromosome-scale records
    (`[INV-223]`), so this cannot be synchronous — and once it is a worker, a
    slow scan can finish *after* the user has switched modes."""

    def test_run_scan_does_not_compute_on_the_calling_thread(self):
        """`_run_scan` dispatches; it must not have the answer when it returns."""
        modal = sc.RegulettuceModal(_rand_seq(4000, seed=61), circular=True,
                                    label="x")
        sent = {}
        modal._scan_worker = lambda *a: sent.setdefault("args", a)
        modal._set_scan_busy = lambda busy: sent.setdefault("busy", busy)
        modal._set_status = lambda markup: sent.setdefault("status", markup)
        modal.query_one = lambda *a, **k: (_ for _ in ()).throw(
            sc.NoMatches("no widgets mounted"))
        modal._run_scan()
        assert sent["args"][0] == "cis"
        assert sent["busy"] is True, "the Scan button must disable while busy"
        assert modal._hits == [], "the scan must not have run inline"

    def test_a_stale_scan_result_is_dropped(self):
        """A scan that finishes after a newer one must not overwrite it."""
        modal = sc.RegulettuceModal("ACGT" * 100, circular=False, label="x")
        modal._fill_table = lambda: None
        modal._set_status = lambda markup: None
        modal._set_scan_busy = lambda busy: None
        modal._scan_gen = 7
        stale = [{"start": 1, "end": 7, "strand": 1, "name": "stale"}]
        modal._scan_done("cis", True, stale, [], [], None, 3)
        assert modal._hits == [], "a gen-3 result landed over gen 7"
        fresh = [{"start": 2, "end": 8, "strand": 1, "name": "fresh"}]
        modal._scan_done("cis", True, fresh, [], [], None, 7)
        assert modal._hits == fresh

    def test_every_mode_switch_bumps_the_generation(self):
        modal = sc.RegulettuceModal(_rand_seq(2000, seed=62), circular=True,
                                    label="x")
        modal._scan_worker = lambda *a: None
        modal._set_scan_busy = lambda busy: None
        modal._set_status = lambda markup: None
        modal.query_one = lambda *a, **k: (_ for _ in ()).throw(
            sc.NoMatches("no widgets mounted"))
        seen = set()
        for mode in ("cis", "3prime", "cis", "3prime"):
            modal._mode = mode
            modal._run_scan()
            seen.add(modal._scan_gen)
        assert len(seen) == 4, f"generations repeated: {seen}"

    def test_worker_error_clears_the_table_and_shows_the_reason(self):
        modal = sc.RegulettuceModal("ACGT" * 100, circular=False, label="x")
        modal._hits = [{"start": 0, "end": 6, "strand": 1, "name": "old"}]
        shown = {}
        modal._fill_table = lambda: None
        modal._set_status = lambda markup: shown.setdefault("s", markup)
        modal._set_scan_busy = lambda busy: None
        modal._scan_done("cis", True, [], None, [], "boom", modal._scan_gen)
        assert modal._hits == []
        assert "boom" in shown["s"]

    @pytest.mark.asyncio
    async def test_tata_without_anchors_answers_without_dispatching(self):
        """No anchor means no scan at all — not a worker that returns []."""
        modal = sc.RegulettuceModal(_rand_seq(3000, seed=63), circular=True,
                                    label="x", anchors=[])

        class _Host(sc.App):
            def on_mount(self) -> None:
                self.push_screen(modal)

        async with _Host().run_test() as pilot:
            await _settle(pilot)
            dispatched = []
            modal._scan_worker = lambda *a: dispatched.append(a)
            modal._mode = "tata"
            modal._rebuild_columns()
            modal._run_scan()
            await _settle(pilot)
            assert dispatched == [], (
                "an unanchored TATA scan has nothing to compute; dispatching "
                "a worker to return [] would read as 'none found'")
            from textual.widgets import Static
            assert "429" in str(modal.query_one("#reg-status", Static).render())


class TestRefusalsCannotBeBypassed:
    """A contract that a degenerate input can slip past is not a contract.

    `_scan_tata_boxes`'s anchors rule is this module's loudest refusal — an
    unanchored TATAWA scan is ~47 chance hits on a 10 kb plasmid. It used to
    check the SEQUENCE first and return `[]` for an empty one, so the one call
    shape that should have been refused outright got "no TATA boxes" instead.
    """

    def test_empty_sequence_still_demands_an_anchor(self):
        with pytest.raises(ValueError, match="requires at least one anchor"):
            rl._scan_tata_boxes("", anchors=[])
        with pytest.raises(ValueError, match="requires at least one anchor"):
            rl._scan_tata_boxes("", anchors=None)

    def test_all_ambiguous_sequence_still_demands_an_anchor(self):
        """`_clean` turns every non-ACGT base into N, so this is non-empty and
        reaches the scan proper — but a caller who passed no anchor still gets
        the refusal, not an empty list."""
        with pytest.raises(ValueError, match="requires at least one anchor"):
            rl._scan_tata_boxes("N" * 500, anchors=[])

    # `None` is NOT in this list: it means "use the measured default".
    @pytest.mark.parametrize("bad", [(1,), (), "xx", 7, object()])
    def test_malformed_window_is_a_valueerror(self, bad):
        """Every caller above translates ValueError into a 400 and lets
        anything else reach the 500 path, so a short tuple must not come back
        as IndexError and a stray string must not come back as a bare int()
        message."""
        with pytest.raises(ValueError, match="low, high"):
            rl._scan_tata_boxes("ACGT" * 50, anchors=[10], window=bad)


class TestAnchorsAreNotWrappedOnALinearMolecule:
    """``% total`` is the coordinate system of a CIRCLE and nothing else.

    On a linear molecule bp -40 and bp ``total - 40`` are different places,
    and one of them does not exist. An anchor typed as 99999 on a 1,206 bp
    linear sequence used to come back as bp 640 and report a TATA box 40 bp
    upstream of a base the molecule has not got."""

    SEQ = ("GC" * 300) + "TATAAA" + ("GC" * 300)

    def test_linear_anchor_past_the_end_raises(self):
        with pytest.raises(ValueError, match="outside a linear"):
            rl._scan_tata_boxes(self.SEQ, anchors=[640 + len(self.SEQ)],
                                circular=False)

    def test_linear_negative_anchor_raises(self):
        with pytest.raises(ValueError, match="outside a linear"):
            rl._scan_tata_boxes(self.SEQ, anchors=[-566], circular=False)

    def test_circular_anchor_still_wraps(self):
        """The wrap is CORRECT on a circle and must survive the fix: bp 640
        and bp 640 + total are the same base there, and the hit at -40 is the
        one `_signed_offset` exists to keep."""
        inside = rl._scan_tata_boxes(self.SEQ, anchors=[640], circular=True)
        wrapped = rl._scan_tata_boxes(self.SEQ,
                                      anchors=[640 + len(self.SEQ)],
                                      circular=True)
        assert [h["offset"] for h in inside] == [-40]
        assert wrapped == inside

    def test_a_linear_anchor_inside_the_molecule_is_untouched(self):
        hits = rl._scan_tata_boxes(self.SEQ, anchors=[640], circular=False)
        assert [h["offset"] for h in hits] == [-40]
        assert [h["anchor"] for h in hits] == [640]


class TestModalHardening:
    """Three dead ends the dialog could reach, all of them silent."""

    def test_typed_anchor_off_the_molecule_is_refused_with_the_range(self):
        """Taking it would wrap on a circle (bp 99999 on a 5 kb plasmid
        becomes bp 4999) and raise on a linear one. Falling back to the
        record's own starts is just as wrong — the user typed a number and the
        table would answer about somewhere else."""
        modal = sc.RegulettuceModal(
            "ACGT" * 250, circular=True, label="x",
            anchors=[{"position": 10, "strand": 1, "label": "cds"}])
        modal.query_one = lambda *a, **k: _FakeInput("99999")
        assert modal._anchor_list() == []
        assert "99,999" in modal._anchor_error
        assert "0 and 999" in modal._anchor_error, (
            "the refusal has to name the range that would work")

    def test_typed_anchor_inside_the_molecule_is_used(self):
        modal = sc.RegulettuceModal("ACGT" * 250, circular=True, label="x",
                                    anchors=[{"position": 10, "strand": 1}])
        modal.query_one = lambda *a, **k: _FakeInput("640")
        assert modal._anchor_list() == [
            {"position": 640, "strand": 1, "label": "typed"}]
        assert modal._anchor_error == ""

    def test_annotate_all_refuses_past_the_cap_rather_than_truncating(self):
        """A 6 bp motif turns up by chance every ~2 kb per strand, so the
        whole catalogue on a 50 kb record is ~2,100 hits and feature COUNT is
        what makes a record slow to draw (`[INV-223]`). Annotating the first
        N silently would leave the user believing the map shows them all."""
        modal = sc.RegulettuceModal("ACGT" * 100, circular=False, label="x")
        modal._hits = [{"start": i, "end": i + 6, "strand": 1, "name": "G-box"}
                       for i in range(modal._MAX_ANNOTATE_ALL + 1)]
        sent, status = {}, {}
        modal.dismiss = lambda v: sent.update({"v": v})
        modal._set_status = lambda m: status.update({"m": m})
        modal._on_annotate_all(None)
        assert sent == {}, "it must not dismiss with an unbounded payload"
        assert "too many" in status["m"]
        # Exactly at the cap still goes through.
        modal._hits = modal._hits[:modal._MAX_ANNOTATE_ALL]
        modal._on_annotate_all(None)
        assert len(sent["v"]["hits"]) == modal._MAX_ANNOTATE_ALL

    @pytest.mark.asyncio
    async def test_unexpected_scanner_error_reaches_the_dialog(self):
        """The worker caught only ValueError. Anything else killed the worker
        with `_set_scan_busy(True)` still in force, so the Scan button stayed
        disabled and the status stayed on "Scanning…" for the rest of the
        session — Escape was the only way out."""
        modal = sc.RegulettuceModal(_rand_seq(2000, seed=77), circular=True,
                                    label="x")

        class _Host(sc.App):
            def on_mount(self) -> None:
                self.push_screen(modal)

        async with _Host().run_test() as pilot:
            await _settle(pilot)
            import splicecraft_modals as mods

            def _boom(*a, **k):
                raise RuntimeError("scanner exploded")

            from textual.widgets import Button, Static
            orig = mods._scan_3prime_risk
            mods._scan_3prime_risk = _boom
            try:
                modal._mode = "3prime"
                modal._rebuild_columns()
                modal._run_scan()
                await _settle(pilot)
            finally:
                mods._scan_3prime_risk = orig
            text = str(modal.query_one("#reg-status", Static).render())
            assert "scanner exploded" in text
            assert modal.query_one("#reg-btn-scan", Button).disabled is False, (
                "a stuck Scan button is a dead end, not an error message")


class _FakeInput:
    """Minimal stand-in for the anchor `Input` — `_anchor_list` reads only
    `.value`, and stubbing the widget keeps this a unit test."""

    def __init__(self, value: str) -> None:
        self.value = value
