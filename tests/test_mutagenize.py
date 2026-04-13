"""
test_mutagenize — SOE-PCR mutagenesis primer design.

Covers the module-level helpers ported from mutagenesis_primers.py:
_mut_parse, _mut_translate, _mut_revcomp, _mut_design_outer,
_mut_design_inner, _mut_extract_cds, and the edge-case branch for
mutations within 60 nt of either CDS end.
"""
import pytest

import splicecraft as sc


# ── Hand-verifiable fixtures ──────────────────────────────────────────────────

# 246-nt realistic CDS (~50 % GC, varied codons) — 81 aa protein + stop.
# Hand-verified aa positions used in tests: 2=A, 3=E, 40=V, 78=A.
_CDS_LONG = (
    "ATG"
    "GCTGAAGTTCAGGATAACCTGGCGAAAGTTCAGGAAGCGGTTGATACCCTGAAACGTGGT"
    "CTGGAAGCGGCGAAAGCGACCCTGGAAAAAGCGGGTGAAGATATCGCGAAAGCGGTTGAT"
    "GGTAAACGTAAAGGCGATCTGGAAAAACTGGCGGAAGCGCTGCAGAAAGTTGAAGCGGAT"
    "ATCGCGAAAGCGGTTGATGGTAAACGTAAAGGCGATCTGGAAAAACTGGCGGAAGCGCTG"
    "TAA"
)


# ── _mut_parse ────────────────────────────────────────────────────────────────

class TestParseMutation:
    def test_basic(self):
        assert sc._mut_parse("W140F") == ("W", 140, "F")

    def test_lowercase_accepted(self):
        assert sc._mut_parse("w140f") == ("W", 140, "F")

    def test_stop_codon(self):
        assert sc._mut_parse("W140*") == ("W", 140, "*")

    def test_rejects_bad_format(self):
        with pytest.raises(ValueError):
            sc._mut_parse("W140")
        with pytest.raises(ValueError):
            sc._mut_parse("140F")
        with pytest.raises(ValueError):
            sc._mut_parse("not a mutation")


# ── _mut_revcomp ──────────────────────────────────────────────────────────────

class TestRevcomp:
    def test_simple(self):
        assert sc._mut_revcomp("ATGC") == "GCAT"

    def test_involutive(self):
        seq = "ATGCATGCGGTTAA"
        assert sc._mut_revcomp(sc._mut_revcomp(seq)) == seq


# ── _mut_translate ────────────────────────────────────────────────────────────

class TestTranslate:
    def test_stops_at_stop_codon(self):
        # "ATG AAA TAA GGG" — stop at codon 3, GGG never translated.
        assert sc._mut_translate("ATGAAATAAGGG") == "MK"

    def test_long_cds(self):
        protein = sc._mut_translate(_CDS_LONG)
        assert len(protein) == 81
        assert protein[0]  == "M"
        assert protein[1]  == "A"   # GCT
        assert protein[2]  == "E"   # GAA


# ── Outer primer design — BsaI tails (regression guard) ───────────────────────

class TestDesignOuter:
    """The outer primers are constant per CDS and must carry the BsaI-AATG
    (FWD) / BsaI-AACG (REV) tails that produce GB B3 / B5 overhangs after
    digestion. Changes here will break every Golden Braid L0 assembly."""

    def test_fwd_carries_bsai_aatg_tail(self):
        outer = sc._mut_design_outer(_CDS_LONG)
        assert outer["fwd"]["full"].startswith("CCCCGGTCTCAAATG")
        assert outer["b3_overhang"] == "AATG"

    def test_rev_carries_bsai_aacg_tail(self):
        outer = sc._mut_design_outer(_CDS_LONG)
        assert outer["rev"]["full"].startswith("CCCCGGTCTCAAACG")
        # The vector-side overhang name is CGTT (= revcomp of AACG on insert)
        assert outer["b5_overhang"] == "CGTT"

    def test_fwd_anneal_starts_after_atg(self):
        """FWD_outer anneal must begin at nt 4 (index 3) of the CDS so the
        AATG overhang reconstitutes the start codon in the assembled part."""
        outer = sc._mut_design_outer(_CDS_LONG)
        assert outer["fwd_anneal_start"] == 3
        # Anneal region is taken from _CDS_LONG[3:3+length]
        anneal = outer["fwd"]["anneal"]
        assert _CDS_LONG[3:3 + len(anneal)] == anneal

    def test_rev_anneal_is_revcomp_of_cds_end(self):
        outer = sc._mut_design_outer(_CDS_LONG)
        anneal = outer["rev"]["anneal"]
        end_rc = sc._mut_revcomp(_CDS_LONG)
        assert end_rc.startswith(anneal)


# ── Inner pair — revcomp invariant, WT codon check ────────────────────────────

class TestDesignInner:
    def test_rev_is_revcomp_of_fwd(self):
        """Inner REV must be the exact revcomp of inner FWD — this is the
        whole point of the SOE joint primer pair.
        Signature: _mut_design_inner(dna, mut_pos_1, mut_aa, wt_aa)."""
        inner = sc._mut_design_inner(_CDS_LONG, 40, "F", "V")  # V40F
        best = inner["candidates"][0]
        assert sc._mut_revcomp(best["fwd"]) == best["rev"]

    def test_wt_codon_mismatch_raises(self):
        """If the caller says WT='W' but the DNA codon at that position
        doesn't actually encode W, we must error rather than produce a
        nonsense mutation primer."""
        # Position 2 is A (GCT); caller claims WT='W' → error.
        with pytest.raises(ValueError, match="mutation says WT='W'"):
            sc._mut_design_inner(_CDS_LONG, 2, "F", "W")

    def test_mut_codon_differs_from_wt(self):
        """mut_codon must encode the requested mutant aa and differ from the
        wt codon so the DNA actually changes."""
        inner = sc._mut_design_inner(_CDS_LONG, 40, "F", "V")  # V40F
        assert inner["wt_codon"] == "GTT"                 # codon 40 of CDS
        assert sc._MUT_CODON_TO_AA[inner["mut_codon"]] == "F"
        assert inner["mut_codon"] != inner["wt_codon"]

    def test_mutation_string_format(self):
        """Mutation string format is WT_AA + pos + MUT_AA."""
        inner = sc._mut_design_inner(_CDS_LONG, 40, "F", "V")  # mut=F, wt=V
        assert inner["mutation"] == "V40F"

    def test_no_alt_codon_for_single_codon_aa(self):
        """Met has only one codon (ATG). Asking to mutate an interior Met
        back to Met must error — there is no alternative codon."""
        cds = "ATG" + ("GCT" * 30) + "ATG" + ("GCT" * 30) + "TAA"
        assert sc._mut_translate(cds)[31] == "M"
        with pytest.raises(ValueError):
            sc._mut_design_inner(cds, 32, "M", "M")


# ── Edge-case branch ──────────────────────────────────────────────────────────

class TestEdgeCase:
    """Mutations within _MUT_MIN_SOE_FRAG (60 nt) of either CDS end must
    trigger the modified-outer branch and skip the inner pair."""

    def test_near_start_triggers_modified_fwd(self):
        # Position 3 is E (codon 3 = GAA). Fragment A ≈ 9 nt → far below 60.
        inner = sc._mut_design_inner(_CDS_LONG, 3, "F", "E")   # E3F
        ec = inner["edge_case"]
        assert ec is not None
        assert ec["near_start"] is True
        assert ec["near_end"] is False
        assert ec["modified_outer"]["label"] == "modified_FWD_outer"
        # The modified FWD carries the BsaI-AATG tail like the normal FWD.
        assert ec["modified_outer"]["full"].startswith("CCCCGGTCTCAAATG")

    def test_near_end_triggers_modified_rev(self):
        # Position 78 is A (codon 78 = GCG). Fragment B ≈ 12 nt → below 60.
        inner = sc._mut_design_inner(_CDS_LONG, 78, "F", "A")  # A78F
        ec = inner["edge_case"]
        assert ec is not None
        assert ec["near_end"] is True
        assert ec["modified_outer"]["label"] == "modified_REV_outer"
        assert ec["modified_outer"]["full"].startswith("CCCCGGTCTCAAACG")

    def test_middle_mutation_no_edge_case(self):
        # Position 40 is V, well away from both ends (~120 nt from either).
        inner = sc._mut_design_inner(_CDS_LONG, 40, "F", "V")  # V40F
        assert inner["edge_case"] is None


# ── CDS extraction — strand and wrap handling ─────────────────────────────────

class TestExtractCds:
    """_mut_extract_cds must return the CDS in its biological 5'→3'
    orientation regardless of strand or origin-wrap."""

    def test_forward_strand_simple(self):
        seq = "AAAA" + _CDS_LONG + "TTTT"
        cds = sc._mut_extract_cds(seq, 4, 4 + len(_CDS_LONG), 1)
        assert cds == _CDS_LONG

    def test_reverse_strand_is_revcomp(self):
        """A CDS on the reverse strand at plasmid[a:b] should come back as
        revcomp(plasmid[a:b]) so the first codon is ATG."""
        rc = sc._mut_revcomp(_CDS_LONG)
        seq = "AAAA" + rc + "TTTT"
        cds = sc._mut_extract_cds(seq, 4, 4 + len(rc), -1)
        assert cds == _CDS_LONG
        assert cds.startswith("ATG")

    def test_wrap_around_origin(self):
        """A feature with end < start spans the origin. The extracted CDS
        must be tail + head, in order."""
        # Build a "plasmid" where the CDS starts near the end and wraps:
        # CDS = ATG + GCT*5 + TAA = 21 nt. Place first 15 nt at the end of
        # the plasmid and the last 6 nt at the start.
        cds = "ATG" + ("GCT" * 5) + "TAA"
        assert len(cds) == 21
        padding = "N" * 30
        # plasmid layout: [last 6 nt of cds][padding][first 15 nt of cds]
        plasmid = cds[15:] + padding + cds[:15]
        start = len(plasmid) - 15  # where the CDS head lives
        end   = 6                  # where the CDS tail ends (wrapped)
        extracted = sc._mut_extract_cds(plasmid, start, end, 1)
        assert extracted == cds
