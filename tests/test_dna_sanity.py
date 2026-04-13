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

    def test_circular_wraparound_ecori_found(self):
        """A palindromic site that spans the origin is found when circular=True."""
        # seq[-3:] + seq[:3] == 'GAATTC' (EcoRI)
        seq = "TTC" + "ACGTACGTACGTACGTACGT" + "GAA"  # 27 bp
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True, circular=True)
        eco = self._resites(feats, "EcoRI")
        # Wrap site emits two pieces; the labeled tail is what counts as a site.
        assert len(eco) == 1, f"Expected 1 labeled wrap-around EcoRI, got {len(eco)}"
        r = eco[0]
        assert r["start"] == len(seq) - 3 and r["end"] == len(seq), r
        # And there is an unlabeled continuation piece on the head.
        all_eco_resites = [f for f in feats
                           if f.get("type") == "resite" and f.get("color") == r["color"]]
        heads = [h for h in all_eco_resites if h.get("start") == 0]
        assert len(heads) == 1 and heads[0]["end"] == 3

    def test_circular_wraparound_not_found_when_linear(self):
        """Same wrap site does NOT appear with circular=False."""
        seq = "TTC" + "ACGTACGTACGTACGTACGT" + "GAA"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True, circular=False)
        assert not self._resites(feats, "EcoRI")

    def test_circular_wraparound_recut_position(self):
        """Wrap-around EcoRI (cuts after G^AATTC). Site spans n-3..n+2, so the
        cut is at (n-3)+1 = n-2 (just after the G, still in the tail)."""
        seq = "TTC" + "ACGTACGTACGTACGTACGT" + "GAA"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True, circular=True)
        cuts = self._recuts(feats, "EcoRI")
        assert len(cuts) == 1
        assert cuts[0]["start"] == len(seq) - 2

    def test_circular_wraparound_unique_filter(self):
        """A single wrap-around site must pass unique_only=True even though
        it produces two resite pieces — only the labeled one counts."""
        seq = "TTC" + "ACGTACGTACGTACGTACGT" + "GAA"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True, circular=True)
        labeled = [f for f in feats if f.get("label") == "EcoRI" and f.get("type") == "resite"]
        assert len(labeled) == 1

    def test_circular_wrap_type_iis_ext_cut_bp_preserved(self):
        """Regression guard for 2026-04-13 fix. A Type IIS enzyme (BsaI:
        recognition + external cut 7 bp downstream) that wraps the origin was
        losing `ext_cut_bp` on both wrap pieces — so the cut-arrow glyph in
        the sequence panel disappeared. `ext_cut_bp` must now be preserved on
        both tail and head pieces regardless of where the absolute cut lands."""
        # Place GGTCTC wrapping origin: last 3 bases 'GGT' + first 3 bases 'CTC'
        n = 30
        seq = "CTC" + "X" * (n - 6) + "GGT"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True, circular=True)
        resites = [f for f in feats
                   if f.get("label") == "BsaI" and f.get("type") == "resite"]
        assert len(resites) == 1
        tail = resites[0]
        assert tail["start"] == 27 and tail["end"] == 30
        # fwd_cut=7 → abs cut position = (27 + 7) % 30 = 4
        assert tail["ext_cut_bp"] == 4, (
            f"ext_cut_bp lost on tail of wrap Type IIS; got {tail['ext_cut_bp']}"
        )
        # Head piece (unlabeled) also carries ext_cut_bp — the chunk-range
        # check in _build_seq_text makes the double-attach idempotent.
        heads = [f for f in feats
                 if f.get("type") == "resite" and f.get("color") == tail["color"]
                 and f["start"] == 0]
        assert len(heads) == 1
        assert heads[0]["ext_cut_bp"] == 4

    def test_circular_wrap_type_iis_cut_arrow_renders(self):
        """End-to-end render check: the cut-arrow glyph (↑ or ↓) must appear
        in the rendered sequence panel output for a wrap Type IIS site.
        Without the 2026-04-13 ext_cut_bp fix, the arrow was silently dropped."""
        n = 30
        seq = "CTC" + "X" * (n - 6) + "GGT"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True, circular=True)
        rendered = sc._build_seq_text(seq, feats, line_width=60).plain
        assert ("↓" in rendered or "↑" in rendered), (
            f"cut arrow missing from rendered wrap Type IIS site:\n{rendered}"
        )

    def test_circular_wrap_reverse_strand_non_palindrome(self):
        """A non-palindromic enzyme whose RC binding sequence spans the origin
        on the forward strand (enzyme binds the reverse strand across the
        origin). The resite must be found with strand=-1 and correct coords."""
        # Place GAGACC (= rc of GGTCTC) across origin: head 'ACC', tail 'GAG'
        n = 30
        seq = "ACC" + "X" * (n - 6) + "GAG"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True, circular=True)
        bsai = [f for f in feats if f.get("label") == "BsaI" and f.get("type") == "resite"]
        assert len(bsai) == 1
        r = bsai[0]
        assert r["strand"] == -1
        # Wrap: tail [n-3, n) labeled; head [0, 3) unlabeled.
        assert r["start"] == 27 and r["end"] == 30

    def test_circular_wrap_recut_per_strand(self):
        """Non-palindromic enzyme with BOTH forward + reverse-strand hits
        emits one `recut` per strand (bottom-strand cut mirrors top-strand cut)."""
        # Put GGTCTC at p=0 forward and GAGACC at p=20 (reverse binding)
        seq = "GGTCTC" + "AAAA" + "GAGACC" + "AAAAAAAAAAAAA"
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=False, circular=False)
        cuts = [f for f in feats if f.get("label") == "BsaI" and f.get("type") == "recut"]
        fwd_cuts = [c for c in cuts if c["strand"] == 1]
        rev_cuts = [c for c in cuts if c["strand"] == -1]
        assert len(fwd_cuts) == 1, fwd_cuts
        assert len(rev_cuts) == 1, rev_cuts

    def test_plasmid_shorter_than_longest_enzyme(self):
        """A plasmid shorter than the scan catalog's max_site_len must not
        crash the augmented-sequence wrap logic. Just verify it runs and
        returns a sensible (possibly empty) list of hits."""
        # Several tiny plasmids — all should scan cleanly without crashing.
        for short_seq in ["", "A", "GAATT", "GAATTC", "GAATTCA"]:
            out = sc._scan_restriction_sites(short_seq, min_recognition_len=6,
                                             unique_only=False, circular=True)
            assert isinstance(out, list)

    def test_circular_multiple_origin_spanning_sites_same_enzyme(self):
        """Two distinct palindromic EcoRI sites, both spanning the origin
        (possible on a very short plasmid) — unique_only=True must drop this
        enzyme since there are 2 labeled resites."""
        # Tricky to set up: we need a palindromic site ≥ 6 bp that can appear
        # twice in a plasmid, each straddling the origin. Instead of forcing
        # multiple wrap sites, confirm the simpler case: 1 wrap + 1 linear.
        # Plasmid: "TTC" + ... + "GAATTC" + ... + "GAA" → 1 wrap at origin +
        # 1 plain linear EcoRI in the middle → unique_only drops EcoRI.
        seq = "TTC" + "ACGT" + "GAATTC" + "ACGTACGT" + "GAA"
        feats_unique = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                                  unique_only=True, circular=True)
        feats_all = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                               unique_only=False, circular=True)
        labeled_unique = [f for f in feats_unique
                          if f.get("label") == "EcoRI" and f.get("type") == "resite"]
        labeled_all = [f for f in feats_all
                       if f.get("label") == "EcoRI" and f.get("type") == "resite"]
        assert labeled_unique == []      # dropped because 2 sites
        assert len(labeled_all) == 2     # 1 linear + 1 wrap

    def test_no_duplicate_match_at_wrap_boundary(self):
        """A recognition that starts exactly at p=n-site_len+1 is fully within
        the augmented sequence at that position AND has a duplicate at p=1 of
        the augmented-only region. The `if p >= n: continue` guard must drop
        the duplicate. Verified by checking resite count matches a hand count."""
        # EcoRI 'GAATTC' at exactly p=n-6: fits without wrap AND doesn't wrap.
        # Also ensure no phantom second match from augmented tail.
        seq = "AAAAAAAA" + "GAATTC"   # 14 bp, EcoRI at p=8 (n-6)
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=True, circular=True)
        eco = [f for f in feats if f.get("label") == "EcoRI" and f.get("type") == "resite"]
        assert len(eco) == 1
        assert eco[0]["start"] == 8 and eco[0]["end"] == 14


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

    def test_reverse_strand_uses_full_iupac_rc(self):
        """Regression for the 2026-04-12 fix: _translate_cds used a bare
        ACGT→TGCA maketrans for reverse-strand RC, which left IUPAC codes
        (R, Y, W, S, M, K, B, D, H, V) unchanged instead of complementing
        them. For example, R (A/G) should complement to Y (C/T), but the
        old code left it as R — producing the wrong codon on the reverse
        strand and a silent mistranslation.

        This test builds a CDS with R (purine = A or G) and verifies that
        the reverse-strand translation matches the forward-strand translation.
        With the old buggy code, the RC'd codon would be ARG→AYG, which
        maps to a different (or unknown) amino acid."""
        # Forward CDS: ATG ARG TAG → has R at position 4
        fwd_cds = "ATGARGTAG"
        from Bio.Seq import Seq
        rc_cds = str(Seq(fwd_cds).reverse_complement())
        # Build a genome with the RC version on the forward strand
        full = "AA" + rc_cds + "AA"
        # Forward translation for reference
        fwd_aa = sc._translate_cds(fwd_cds, 0, 9, strand=1)
        # Reverse-strand translation should produce the SAME protein
        rev_aa = sc._translate_cds(full, 2, 11, strand=-1)
        assert rev_aa == fwd_aa, (
            f"reverse-strand translation should match forward: "
            f"fwd={fwd_aa!r} rev={rev_aa!r} "
            f"(old bug: R not complemented to Y in RC step)"
        )

    def test_wrapped_cds_forward_strand(self):
        """Regression guard for 2026-04-13 fix. An origin-spanning CDS stored
        as `end < start` must concatenate tail + head before translating.
        Pre-fix behaviour: `full_seq[start:end]` returns "" for end<start, so
        a wrapped CDS silently translated to empty — users saw (0 aa) for a
        real protein."""
        # Plasmid of length 11; CDS = tail 'ATG' at bp 8..11 + head 'AAATAG' at bp 0..6
        # → forward reads 'ATG' + 'AAATAG' = 'ATGAAATAG' = MK*
        seq = "AAATAG" + "XX" + "ATG"
        assert len(seq) == 11
        aa = sc._translate_cds(seq, start=8, end=6, strand=1)
        assert aa == "MK*", f"wrapped CDS forward translation wrong: {aa!r}"

    def test_wrapped_cds_reverse_strand(self):
        """As above but the CDS is on the reverse strand. After concatenating
        tail + head we must RC before translating."""
        # RC of 'ATGAAATAG' = 'CTATTTCAT'. Put CTATTT at tail and CAT at head:
        # tail [5,11) = 'CTATTT', head [0,3) = 'CAT'. Reading forward tail+head
        # = 'CTATTTCAT', RC = 'ATGAAATAG' → MK*
        seq = "CAT" + "YY" + "CTATTT"
        assert len(seq) == 11
        aa = sc._translate_cds(seq, start=5, end=3, strand=-1)
        assert aa == "MK*", f"wrapped CDS reverse translation wrong: {aa!r}"

    def test_wrapped_cds_preserves_iupac_across_join(self):
        """IUPAC codes that happen to span the tail/head boundary must still
        translate via `?` rather than silently becoming a letter. Tests that
        concatenation doesn't destroy IUPAC semantics."""
        # CDS 'ATG' + 'NNA TAG' wraps origin: tail 'ATG' + head 'NNATAG'
        # 'ATGNNATAG' → ATG NNA TAG → M ? *
        seq = "NNATAG" + "XX" + "ATG"
        aa = sc._translate_cds(seq, start=8, end=6, strand=1)
        assert aa == "M?*"

    def test_wrapped_cds_length_not_multiple_of_three(self):
        """Wrapped CDS with trailing partial codon: the extra 1-2 nt must be
        dropped identically to the non-wrap partial-codon case."""
        # tail 'ATG' + head 'AAAT' = 'ATGAAAT' (7 nt) → MK + 1 leftover → 'MK*'
        seq = "AAAT" + "XXXX" + "ATG"
        assert len(seq) == 11
        aa = sc._translate_cds(seq, start=8, end=4, strand=1)
        assert aa == "MK*"

    def test_empty_cds_returns_empty(self):
        """start == end is treated as empty (consistent with `_bp_in` which
        returns False for zero-width features)."""
        assert sc._translate_cds("ATGAAATAG", 5, 5, strand=1) == ""
        assert sc._translate_cds("ATGAAATAG", 5, 5, strand=-1) == ""

    def test_lowercase_input_is_accepted(self):
        """Lowercase DNA must be folded to uppercase before translation."""
        assert sc._translate_cds("atgaaatag", 0, 9, strand=1) == "MK*"
        # Reverse strand with lowercase
        from Bio.Seq import Seq
        full = "nn" + str(Seq("atgaaatag").reverse_complement()) + "nn"
        assert sc._translate_cds(full, 2, 11, strand=-1) == "MK*"

    def test_single_base_and_two_base_cds_return_empty(self):
        """1 or 2 bp can't form a codon → empty string, no crash."""
        assert sc._translate_cds("ATGAAA", 0, 1, strand=1) == ""
        assert sc._translate_cds("ATGAAA", 0, 2, strand=1) == ""

    def test_wrapped_cds_codon_spans_origin(self):
        """A codon whose three bases are split across the tail/head boundary
        must still decode correctly — the wrap fix concatenates tail+head so
        codon boundaries are a natural consequence of the joined string."""
        # tail 'AT' + head 'GAAATAG' = 'ATGAAATAG' = MK*
        # Here the first codon 'ATG' crosses the origin: A,T are at tail
        # positions [n-2, n) and G is at head position 0.
        seq = "GAAATAG" + "XX" + "AT"
        assert len(seq) == 11
        aa = sc._translate_cds(seq, start=9, end=7, strand=1)
        assert aa == "MK*"

    def test_wrapped_cds_matches_biopython(self):
        """Cross-check wrap translation against Biopython on a bigger CDS."""
        from Bio.Seq import Seq
        # 30-bp CDS: ATG + 27 codons of random sense + TAA = 10 aa + stop
        cds = "ATG" + "AAACCCGGGTTTAAACCCGGGTTT" + "TAA"   # 30 bp
        assert len(cds) == 30
        expected = str(Seq(cds).translate())
        # Embed with a wrap: plasmid = cds[15:] + 'XX' + cds[:15]
        # So CDS spans tail [17, 32) = cds[15:30] and head [0, 15) = cds[0:15]
        plasmid = cds[15:] + "XX" + cds[:15]
        assert len(plasmid) == 32
        aa = sc._translate_cds(plasmid, start=17, end=15, strand=1)
        # _translate_cds force-appends '*' if absent; expected already has one
        if not expected.endswith("*"):
            expected += "*"
        assert aa == expected
