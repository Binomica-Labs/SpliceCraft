"""
test_dna_sanity — biological correctness tests for SpliceCraft.

If any of these fail, DO NOT SHIP. These tests protect the sacred invariants
listed in CLAUDE.md. Every assertion is either:

  (a) cross-validated against Biopython's reference implementation, or
  (b) verifiable by hand on a short sequence you can stare at, or
  (c) a regression guard for a specific past bug (comment cites the commit
      or the dated entry in CLAUDE.md).

Organised by subsystem:

  TestReverseComplement      — `_rc` IUPAC-aware RC, involution
  TestIUPACPattern           — `_iupac_pattern` degenerate regex + cache
  TestCodonTable             — SpliceCraft's hand-rolled codon table vs Biopython
  TestNEBEnzymesCatalog      — `_NEB_ENZYMES` dict schema and key uniqueness
  TestRestrictionScan        — `_scan_restriction_sites` palindrome-aware scanning
  TestTranslateCds           — `_translate_cds` forward/reverse strands
"""
from __future__ import annotations

import random
import re

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Reverse complement — sacred invariant #3 (CLAUDE.md)
# ═══════════════════════════════════════════════════════════════════════════════

class TestReverseComplement:
    """`_rc()` must handle full IUPAC, not just ACGT. Regression guard for the
    2026-03-30 fix where ambiguity codes were passing through unchanged."""

    def test_rc_acgt_ground_truth(self):
        # Hand-verifiable: ATGC → GCAT
        assert sc._rc("ATGC") == "GCAT"
        # Longer: AAATTT → AAATTT (palindrome)
        assert sc._rc("AAATTT") == "AAATTT"
        # EcoRI recognition site — palindromic
        assert sc._rc("GAATTC") == "GAATTC"
        # BsaI recognition — NOT palindromic
        assert sc._rc("GGTCTC") == "GAGACC"

    def test_rc_matches_biopython(self):
        """Cross-check against Biopython on random ACGT sequences."""
        from Bio.Seq import Seq
        rng = random.Random(0xC0DE)
        for _ in range(100):
            n = rng.randint(1, 200)
            seq = "".join(rng.choice("ACGT") for _ in range(n))
            assert sc._rc(seq) == str(Seq(seq).reverse_complement())

    @pytest.mark.parametrize("code,comp", [
        # IUPAC complement pairs from the table in splicecraft.py:433
        ("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"),
        ("R", "Y"), ("Y", "R"),   # R={A,G} ↔ Y={C,T}
        ("W", "W"), ("S", "S"),   # W={A,T}, S={C,G} — self-complementary
        ("M", "K"), ("K", "M"),   # M={A,C} ↔ K={G,T}
        ("B", "V"), ("V", "B"),   # B={!A} ↔ V={!T}
        ("D", "H"), ("H", "D"),   # D={!C} ↔ H={!G}
        ("N", "N"),               # any ↔ any
    ])
    def test_rc_handles_each_iupac_code(self, code, comp):
        assert sc._rc(code) == comp

    def test_rc_preserves_length(self):
        for n in [0, 1, 2, 6, 20, 100]:
            seq = "ACGTRY" * (n // 6 + 1)
            seq = seq[:n]
            assert len(sc._rc(seq)) == n

    def test_rc_is_involutive(self):
        """rc(rc(x)) == x for all IUPAC-containing strings — 2026-03-30 regression."""
        rng = random.Random(0xBEEF)
        alphabet = "ACGTRYWSMKBDHVN"
        for _ in range(200):
            n = rng.randint(0, 100)
            seq = "".join(rng.choice(alphabet) for _ in range(n))
            assert sc._rc(sc._rc(seq)) == seq

    def test_rc_uppercases(self):
        """Lowercase input must be folded to uppercase — all other code paths
        (restriction scan, pattern matching) assume uppercase."""
        assert sc._rc("acgt") == "ACGT"
        assert sc._rc("gaattc") == "GAATTC"


# ═══════════════════════════════════════════════════════════════════════════════
# IUPAC regex patterns — sacred invariant #4
# ═══════════════════════════════════════════════════════════════════════════════

class TestIUPACPattern:
    def test_plain_acgt_pattern_matches_literal(self):
        p = sc._iupac_pattern("GAATTC")
        assert p.search("TTTGAATTCAAA")
        assert not p.search("TTTGAATTTAAA")

    def test_degenerate_pattern_expansions(self):
        # AvaI = CYCGRG : Y={C,T}, R={A,G} → 4 canonical sequences
        p = sc._iupac_pattern("CYCGRG")
        matches = {"CCCGAG", "CCCGGG", "CTCGAG", "CTCGGG"}
        nonmatches = {"CACGAG", "CCCGCC", "CTCGCC"}
        for m in matches:
            assert p.search(m), f"CYCGRG must match {m}"
        for m in nonmatches:
            assert not p.search(m), f"CYCGRG must NOT match {m}"

    def test_n_is_any_base(self):
        p = sc._iupac_pattern("GGTNACC")    # BstEII
        for b in "ACGT":
            assert p.search(f"GGT{b}ACC"), f"N must match {b}"

    def test_pattern_cache_is_populated(self):
        """The cache (sacred invariant #4) must hold the compiled pattern."""
        sc._PATTERN_CACHE.clear()
        sc._iupac_pattern("ATGCAT")
        assert "ATGCAT" in sc._PATTERN_CACHE

    def test_pattern_cache_returns_same_object(self):
        """Second call must return the SAME compiled re.Pattern, not a fresh one."""
        sc._PATTERN_CACHE.clear()
        p1 = sc._iupac_pattern("GCCGGC")
        p2 = sc._iupac_pattern("GCCGGC")
        assert p1 is p2


# ═══════════════════════════════════════════════════════════════════════════════
# Codon table
# ═══════════════════════════════════════════════════════════════════════════════

class TestCodonTable:
    """SpliceCraft's hand-rolled `_CODON_TABLE` must equal Biopython's standard
    table. If Biopython is ever updated to a different table, this test will
    fail noisily rather than silently mistranslate."""

    def test_codon_table_is_complete(self):
        """Must have exactly 64 entries — all possible trinucleotides."""
        assert len(sc._CODON_TABLE) == 64
        # Every entry is a 3-letter ACGT codon
        for codon in sc._CODON_TABLE:
            assert len(codon) == 3
            assert set(codon) <= set("ACGT")

    def test_stop_codons_are_canonical(self):
        stops = {c for c, aa in sc._CODON_TABLE.items() if aa == "*"}
        assert stops == {"TAA", "TAG", "TGA"}

    def test_atg_is_methionine(self):
        assert sc._CODON_TABLE["ATG"] == "M"

    def test_matches_biopython_standard_table(self):
        """Translate each of the 64 codons with Biopython and compare."""
        from Bio.Seq import Seq
        from Bio.Data.CodonTable import standard_dna_table
        for codon, sc_aa in sc._CODON_TABLE.items():
            if sc_aa == "*":
                assert codon in standard_dna_table.stop_codons
            else:
                bp_aa = standard_dna_table.forward_table[codon]
                assert sc_aa == bp_aa, (
                    f"codon {codon}: splicecraft={sc_aa} biopython={bp_aa}"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# NEB enzyme catalog
# ═══════════════════════════════════════════════════════════════════════════════

class TestNEBEnzymesCatalog:
    """Guard against duplicate keys (the 2026-03-30 SbfI/NspI incident) and
    malformed entries (wrong tuple length, out-of-range cut positions)."""

    def test_at_least_180_enzymes(self):
        # We advertise ~200; a drop below 180 means something was deleted.
        assert len(sc._NEB_ENZYMES) >= 180

    def test_no_duplicate_enzyme_names(self):
        # dict literal auto-dedupes, but this still catches a hypothetical
        # refactor to a list-of-tuples or a merge from another source.
        names = list(sc._NEB_ENZYMES.keys())
        assert len(names) == len(set(names))

    def test_each_entry_is_well_formed(self):
        for name, entry in sc._NEB_ENZYMES.items():
            assert isinstance(entry, tuple), f"{name}: not a tuple"
            assert len(entry) == 3, f"{name}: expected (site, fwd, rev)"
            site, fwd, rev = entry
            assert isinstance(site, str) and site
            # Every character must be a valid IUPAC code
            assert set(site) <= set("ACGTRYWSMKBDHVN"), (
                f"{name}: unknown base in '{site}'"
            )
            assert isinstance(fwd, int)
            assert isinstance(rev, int)

    def test_every_site_has_a_color(self):
        """`_RESTR_COLOR` is indexed by enzyme name in the scan path — every
        key in `_NEB_ENZYMES` must also be present in `_RESTR_COLOR`."""
        for name in sc._NEB_ENZYMES:
            assert name in sc._RESTR_COLOR, f"{name}: missing from _RESTR_COLOR"

    def test_every_site_compiles_as_iupac_pattern(self):
        """Every recognition sequence must be a valid regex after IUPAC expansion."""
        for name, (site, _, _) in sc._NEB_ENZYMES.items():
            try:
                pat = sc._iupac_pattern(site)
                assert isinstance(pat, re.Pattern)
            except re.error as e:
                pytest.fail(f"{name}: site {site!r} does not compile: {e}")

    def test_common_enzymes_present(self):
        for must in ["EcoRI", "BamHI", "HindIII", "NcoI", "XbaI", "PstI",
                     "BsaI", "BsmBI", "BbsI", "SapI"]:
            assert must in sc._NEB_ENZYMES, f"{must} missing from catalog"

    def test_ecori_site_is_canonical(self):
        site, fwd, rev = sc._NEB_ENZYMES["EcoRI"]
        assert site == "GAATTC"
        # G^AATTC : cut between position 0 and 1 of recognition → fwd_cut = 1
        assert fwd == 1

    def test_bsai_is_type_iis(self):
        """BsaI = GGTCTC(1/5). The cut must fall OUTSIDE the 6-bp recognition
        site (fwd_cut > site_len). This is the defining Type IIS property."""
        site, fwd, rev = sc._NEB_ENZYMES["BsaI"]
        assert site == "GGTCTC"
        assert fwd > len(site), (
            f"BsaI fwd_cut={fwd} should be > {len(site)} (Type IIS)"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Restriction site scanning — sacred invariants #1 and #2
# ═══════════════════════════════════════════════════════════════════════════════

class TestRestrictionScan:
    """These tests protect the two hardest-won bugfixes from 2026-03-30:
      1. Palindromic enzymes must not be double-counted (forward + reverse scan
         was producing 2 resites per physical site).
      2. Reverse-strand hits for NON-palindromic enzymes must use the forward
         coordinate `p`, not `n - p - site_len`.
    """

    @staticmethod
    def _resites(feats, enzyme=None):
        out = [f for f in feats if f["type"] == "resite"]
        if enzyme is not None:
            out = [f for f in out if f["label"] == enzyme]
        return out

    @staticmethod
    def _recuts(feats, enzyme=None):
        out = [f for f in feats if f["type"] == "recut"]
        if enzyme is not None:
            out = [f for f in out if f["label"] == enzyme]
        return out

    def test_ecori_single_site_not_double_counted(self):
        """Sacred invariant #1: palindromic EcoRI must give exactly 1 resite
        even though the scanner looks at both strands."""
        seq = "AAA" + "GAATTC" + "AAA"   # 12 bp, one EcoRI site at position 3
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True)
        eco_resites = self._resites(feats, "EcoRI")
        assert len(eco_resites) == 1, (
            f"EcoRI should have 1 resite, got {len(eco_resites)}"
        )
        r = eco_resites[0]
        assert r["start"] == 3 and r["end"] == 9
        assert r["strand"] == 1   # palindrome → reported on forward strand

    def test_ecori_three_sites(self):
        seq = "AA" + "GAATTC" + "AAA" + "GAATTC" + "AAA" + "GAATTC" + "AA"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=False)
        eco = self._resites(feats, "EcoRI")
        assert len(eco) == 3

    def test_palindromes_produce_one_recut_per_site(self):
        """Each palindromic cut should emit exactly 1 `recut` entry
        (bottom-strand recut was part of the 2026-03-30 fix)."""
        seq = "AAA" + "GAATTC" + "AAA"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True)
        cuts = self._recuts(feats, "EcoRI")
        assert len(cuts) == 1
        # EcoRI cuts G^AATTC → position 3 (start of recognition) + 1 = 4
        assert cuts[0]["start"] == 4

    def test_non_palindrome_on_forward_strand(self):
        """BsaI forward site: recognition at p, cut OUTSIDE recognition site."""
        prefix = "AAAAAA"                         # 6 bp
        site = "GGTCTC"                           # BsaI recognition, 6 bp
        spacer_then_cut = "NNNNNNNNNN"            # 10 bp after, so cut (at +7) is valid
        seq = prefix + site + spacer_then_cut
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True)
        resites = self._resites(feats, "BsaI")
        assert len(resites) == 1
        r = resites[0]
        assert r["strand"] == 1
        assert r["start"] == 6 and r["end"] == 12

        # The associated recut lies OUTSIDE the recognition span (Type IIS)
        cuts = self._recuts(feats, "BsaI")
        assert len(cuts) == 1
        cut_bp = cuts[0]["start"]
        # BsaI fwd_cut=7 → recognition starts at 6, cut at 6+7=13
        assert cut_bp == 13
        assert cut_bp >= r["end"]   # cut is downstream of recognition

    def test_non_palindrome_on_reverse_strand_uses_forward_coordinate(self):
        """Sacred invariant #2: when BsaI is found via its RC on the forward
        strand (i.e. the enzyme binds the reverse strand), the recorded
        `start` must be the forward-strand coordinate `p`, NOT the pre-2026-
        03-30 buggy value `n - p - site_len`.
        """
        # Put GAGACC (= rc of GGTCTC) at a known forward position.
        prefix = "AAAAAAAAA"        # 9 bp
        rc_site = "GAGACC"          # reverse-strand BsaI binding
        suffix = "AAAAAAAAA"        # 9 bp  (total 24)
        seq = prefix + rc_site + suffix
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True)
        resites = self._resites(feats, "BsaI")
        assert len(resites) == 1
        r = resites[0]
        assert r["strand"] == -1
        # Critical: start is 9 (position of GAGACC on forward strand),
        # not len(seq) - 9 - 6 = 9 — in this symmetric test they happen to
        # coincide, so use an asymmetric one to prove the point.
        assert r["start"] == 9
        assert r["end"] == 15

    def test_non_palindrome_reverse_strand_asymmetric(self):
        """Asymmetric positioning: GAGACC near the left edge so the buggy
        `n - p - site_len` coordinate would clearly differ from the correct
        `p`."""
        seq = "AA" + "GAGACC" + "AAAAAAAAAAAAAAAAAAAA"   # 28 bp
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True)
        resites = self._resites(feats, "BsaI")
        assert len(resites) == 1
        r = resites[0]
        assert r["strand"] == -1
        assert r["start"] == 2            # correct: position of GAGACC
        # Buggy formula would give len(seq) - 2 - 6 = 20 — explicitly reject it:
        assert r["start"] != len(seq) - 2 - 6

    def test_unique_only_filter_excludes_multi_cutters(self):
        # Two EcoRI sites: unique_only=True should drop EcoRI entirely.
        seq = "AA" + "GAATTC" + "AAA" + "GAATTC" + "AA"
        feats_unique = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                                  unique_only=True)
        feats_all = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                               unique_only=False)
        assert not self._resites(feats_unique, "EcoRI")
        assert len(self._resites(feats_all, "EcoRI")) == 2

    def test_min_length_filter_excludes_4_cutters(self):
        """HaeIII = GGCC (4 bp). Should be skipped when min_recognition_len=6."""
        seq = "AAA" + "GGCC" + "AAA"    # HaeIII site
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=False)
        assert not self._resites(feats, "HaeIII")
        # And included when min_recognition_len=4
        feats4 = sc._scan_restriction_sites(seq, min_recognition_len=4,
                                            unique_only=False)
        assert self._resites(feats4, "HaeIII")

    def test_empty_sequence_yields_empty_list(self):
        assert sc._scan_restriction_sites("") == []
        assert sc._scan_restriction_sites("AAAA") == []

    def test_degenerate_site_matches_all_expansions(self):
        """BstEII = GGTNACC — N is any base. Seed one of each and verify the
        scan integrates _iupac_pattern correctly. BstEII is chosen over AvaI
        because AvaI's canonical expansions (CTCGAG, CCCGGG) overlap with
        XhoI and SmaI, which the scan dedups as isoschizomers."""
        seq = (
            "AA" + "GGTAACC" + "AA"   # N = A
            + "GGTCACC" + "AA"         # N = C
            + "GGTGACC" + "AA"         # N = G
            + "GGTTACC" + "AA"         # N = T
        )
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=False)
        bst = self._resites(feats, "BstEII")
        assert len(bst) == 4


# ═══════════════════════════════════════════════════════════════════════════════
# CDS translation
# ═══════════════════════════════════════════════════════════════════════════════

class TestTranslateCds:
    def test_forward_strand_simple(self):
        # ATG AAA TAG → M K *
        seq = "ATGAAATAG"
        aa = sc._translate_cds(seq, 0, 9, strand=1)
        assert aa == "MK*"

    def test_forward_strand_with_context(self):
        # Stop is added if missing
        pre  = "GGGG"
        cds  = "ATGAAAAAA"   # M K K (no stop in window)
        post = "CCCC"
        full = pre + cds + post
        aa = sc._translate_cds(full, 4, 13, strand=1)
        # _translate_cds appends a trailing * if the last codon isn't a stop
        assert aa == "MKK*"

    def test_reverse_strand_matches_forward_after_rc(self):
        """The same CDS, placed on the reverse strand, must translate to the
        same protein when called with strand=-1."""
        from Bio.Seq import Seq
        fwd_cds = "ATGAAATAG"   # M K *
        rc_cds  = str(Seq(fwd_cds).reverse_complement())
        # Build a genome where the RC version appears on the forward strand;
        # calling _translate_cds with strand=-1 should RC it back and translate.
        full = "NN" + rc_cds + "NN"
        aa = sc._translate_cds(full, 2, 11, strand=-1)
        assert aa == "MK*"

    def test_all_stops_terminate(self):
        for stop in ("TAA", "TAG", "TGA"):
            seq = "ATGAAA" + stop
            aa = sc._translate_cds(seq, 0, 9, strand=1)
            assert aa.endswith("*")
            assert len(aa) == 3   # M K *

    def test_partial_codon_at_end_is_dropped(self):
        """Translation stops at the last complete codon; a trailing 1–2 bp
        overhang must not produce a '?' or crash."""
        seq = "ATGAAA" + "T"        # 7 bp; last codon is incomplete
        aa = sc._translate_cds(seq, 0, 7, strand=1)
        # Exactly 2 amino acids should translate; trailing stop appended.
        assert aa.replace("*", "") == "MK"

    def test_unknown_codon_becomes_question_mark(self):
        # Inserting 'N' creates a codon not in the table
        seq = "ATG" + "ANA" + "TAG"
        aa = sc._translate_cds(seq, 0, 9, strand=1)
        assert "?" in aa
