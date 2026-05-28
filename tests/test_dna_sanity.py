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


class TestNEBEnzymesAgreeWithBioPython:
    """Cross-validate `_NEB_ENZYMES` tuples against `Bio.Restriction` for every
    enzyme name shared between the two catalogs. Regression guard added
    2026-05-11 after an audit caught ~25 wrong tuples (wrong overhang length
    on BsmI/BstEII/SfiI/Eco53kI/BmtI/BstAPI/BseYI/BsrI/BtsCI/EarI/BssSI; wrong
    cut location on Type IIS far-cutters HphI/MboII/MmeI/NmeAIII/BceAI/BmrI/
    BsbI/BseMII/BseRI). Convention is `(fst5, size + fst3)` where `fst5` /
    `fst3` are BioPython's cut offsets — fst5 = top-strand cut from 5' of
    recognition, fst3 = bottom-strand cut from 3' end of recognition. This
    test would have caught every entry the audit found.

    HF / v2 variants are validated separately (they share a recognition and
    cut pattern with their parent, but BioPython doesn't always carry them)."""

    @pytest.fixture(scope="class")
    def biopy(self):
        import Bio.Restriction as R
        return R

    def test_every_overlapping_enzyme_matches_biopython(self, biopy):
        """For every name in BOTH catalogs, the SpliceCraft tuple's top and
        bottom cut positions must equal `(fst5, size + fst3)` from BioPython."""
        mismatches = []
        for name, (_site, sc_top, sc_bot) in sc._NEB_ENZYMES.items():
            if not hasattr(biopy, name):
                continue
            enz = getattr(biopy, name)
            expected = (enz.fst5, enz.size + enz.fst3)
            if (sc_top, sc_bot) != expected:
                mismatches.append((name, (sc_top, sc_bot), expected,
                                   enz.elucidate()))
        if mismatches:
            lines = [f"  {n:<12} SC={sc} BioPy={exp}  {elu}"
                     for n, sc, exp, elu in mismatches]
            pytest.fail(
                f"{len(mismatches)} enzymes disagree with BioPython:\n"
                + "\n".join(lines)
            )

    def test_every_enzyme_exists_in_rebase(self):
        """Regression guard for issue #14 (cory-mozza, 2026-05-11):
        SpliceCraft was listing enzymes that don't exist (BspLU11III) or
        have no commercial supplier (BsbI). The first crashes any wet-lab
        attempt; the second silently suggests a digest the user can never
        order. Every catalog entry must either appear in REBASE under
        that exact name, or be an HF/v2 variant of a REBASE-listed parent
        — with at least one supplier."""
        import Bio.Restriction.Restriction_Dictionary as RD
        rd = RD.rest_dict
        hallucinated = []
        no_supplier  = []
        for name in sc._NEB_ENZYMES:
            base = name
            for suffix in ("-HF", "-v2"):
                if name.endswith(suffix):
                    base = name[: -len(suffix)]
                    break
            if base not in rd:
                hallucinated.append(name)
                continue
            if not rd[base].get("suppl"):
                no_supplier.append(name)
        assert not hallucinated, (
            f"{len(hallucinated)} hallucinated / non-REBASE enzymes "
            f"in catalog: {hallucinated}"
        )
        assert not no_supplier, (
            f"{len(no_supplier)} enzymes with no commercial supplier "
            f"(user can't buy them) in catalog: {no_supplier}"
        )

    def test_recognition_sites_match_rebase(self):
        """Catalog recognition sequence must match REBASE for every
        non-variant enzyme. Caught AccI (was GTYRAC, real is GTMKAC —
        different IUPAC codes match different sequence sets) and BstXI
        (was 5 Ns between CCA/TGG, real is 6 Ns) on 2026-05-11."""
        import Bio.Restriction.Restriction_Dictionary as RD
        rd = RD.rest_dict
        mismatches = []
        for name, (sc_site, _, _) in sc._NEB_ENZYMES.items():
            if name.endswith("-HF") or name.endswith("-v2"):
                continue  # commercial variant; recognition equals parent
            if name not in rd:
                continue  # caught by `test_every_enzyme_exists_in_rebase`
            rb_site = rd[name]["site"]
            if sc_site.upper() != rb_site.upper():
                mismatches.append((name, sc_site, rb_site))
        if mismatches:
            lines = [f"  {n}: SC={s!r} REBASE={r!r}"
                     for n, s, r in mismatches]
            pytest.fail(
                f"{len(mismatches)} recognition sites disagree with REBASE:\n"
                + "\n".join(lines)
            )

    def test_hf_and_v2_variants_match_parent(self):
        """`*-HF` and `*-v2` variants are enzyme-prep improvements (better
        buffer compatibility, faster cleavage) — they share recognition and
        cleavage with their parent. Diverging values silently break cloning
        sims using the variant name."""
        for name, (site, top, bot) in sc._NEB_ENZYMES.items():
            for suffix in ("-HF", "-v2"):
                if not name.endswith(suffix):
                    continue
                parent = name[: -len(suffix)]
                if parent not in sc._NEB_ENZYMES:
                    continue
                p_site, p_top, p_bot = sc._NEB_ENZYMES[parent]
                assert (site, top, bot) == (p_site, p_top, p_bot), (
                    f"{name} {(site, top, bot)} disagrees with parent "
                    f"{parent} {(p_site, p_top, p_bot)}"
                )

    @pytest.mark.parametrize("name,padding", [
        ("EcoRI",   "AAAAA"),
        ("BamHI",   "AAAAA"),
        ("BsaI",    "AAAAAAAAAA"),
        ("BsmBI",   "AAAAAAAAAA"),
        ("BbsI",    "AAAAAAAAAA"),
        ("SapI",    "AAAAAAAAAA"),
        ("BsmI",    "AAAAAAAAAA"),
        ("BsrI",    "AAAAAAAAAA"),
        ("SfiI",    "AAAAAAAAAA"),
        ("HphI",    "AAAAAAAAAAAAAAAAAA"),
        ("MmeI",    "AAAAAAAAAAAAAAAAAAAAAAAAAA"),
        ("BssSI",   "AAAAAAAAAA"),
        ("Eco53kI", "AAAAAAAAAA"),
        ("BmtI",    "AAAAAAAAAA"),
        ("BstEII",  "AAAAAAAAAA"),
        ("BseYI",   "AAAAAAAAAA"),
    ])
    def test_synthetic_template_top_cut_matches_biopython(self, biopy, name, padding):
        """End-to-end: run `_enzyme_cuts` on a synthetic `pad + site + pad`
        template and confirm the top-strand cut position equals BioPython's
        cut for the same enzyme on the same template."""
        from Bio.Seq import Seq
        enz = getattr(biopy, name)
        # Resolve any IUPAC bases in the recognition site to a single
        # canonical match so BioPython's strict-match search will hit.
        site = (enz.site.replace("N", "A").replace("R", "A")
                          .replace("Y", "C").replace("W", "A")
                          .replace("S", "C").replace("M", "A")
                          .replace("K", "G").replace("B", "C")
                          .replace("D", "A").replace("H", "A")
                          .replace("V", "A"))
        seq = padding + site + padding
        bp_cut_1based = enz.search(Seq(seq), linear=True)
        assert bp_cut_1based, f"BioPython didn't find {name} site in synthetic template"
        bp_cut_0based = bp_cut_1based[0] - 1
        sc_hits = sc._enzyme_cuts(seq, [name], circular=False)
        sc_top_cuts = [h["top"] for h in sc_hits if h.get("top") is not None]
        assert sc_top_cuts, f"SpliceCraft didn't return a cut for {name}"
        assert bp_cut_0based in sc_top_cuts, (
            f"{name}: BioPython top cut at {bp_cut_0based} (0-based) "
            f"not in SpliceCraft top cuts {sc_top_cuts}"
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

    def test_non_palindrome_reverse_strand_recut_position(self):
        """Pre-2026-04-25, line 922's `(p + site_len - 1 - fwd_cut) % n`
        formula was off by one: for BsaI (GGTCTC, fwd_cut=7, rev_cut=11) the
        bottom-strand cut on a reverse-bound site landed one base too far to
        the LEFT.

        Derivation: BsaI's enzyme-local inter-base cut is at position 7 (=
        fwd_cut). On a reverse binding the enzyme's local 0 sits at top
        position p+site_len-1, so local-7 maps to top p+site_len-1-7 and the
        cut splits between local-6 (top p+site_len-1-6) and local-7. In top
        5'→3' order the base immediately right of the cut is at position
        `p + site_len - fwd_cut`, NOT `p + site_len - 1 - fwd_cut`.
        """
        # Place GAGACC (= rc of GGTCTC) at position 10. Sequence padded so
        # all rev-strand cuts land in-bounds for a clean linear test.
        seq = ("A" * 10) + "GAGACC" + ("A" * 10)        # 26 bp
        feats = sc._scan_restriction_sites(seq, min_recognition_len=6,
                                           unique_only=False, circular=False)
        cuts = [f for f in feats
                if f.get("label") == "BsaI" and f.get("type") == "recut"]
        assert len(cuts) == 1
        c = cuts[0]
        assert c["strand"] == -1
        # Correct: p + site_len - fwd_cut = 10 + 6 - 7 = 9.
        # Buggy:   p + site_len - 1 - fwd_cut = 10 + 6 - 1 - 7 = 8.
        assert c["start"] == 9, (
            f"Reverse-strand BsaI bottom-strand cut landed at "
            f"{c['start']} (expected 9 — off-by-one regression?)"
        )

        # And the top-strand `ext_cut_bp` (rev_cut path in
        # `_scan_restriction_sites`) should land at
        # p + site_len - rev_cut = 10 + 6 - 11 = 5, NOT 4 (buggy).
        resites = [f for f in feats
                   if f.get("label") == "BsaI" and f.get("type") == "resite"]
        assert len(resites) == 1
        assert resites[0]["ext_cut_bp"] == 5, (
            f"Reverse-strand BsaI top-strand ext_cut_bp at "
            f"{resites[0]['ext_cut_bp']} (expected 5 — off-by-one regression?)"
        )

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


class TestResiteLabelWidth:
    """Resite parens widen to fit the full enzyme short-form name.

    Per user spec 2026-05-08: "make sure the full enzyme name fits
    the enzyme label i.e EcoRI not (EcoR). Do not scale the
    parentheses width to the recognition site width. Keep it just
    wide enough to fit the full shorthand enzyme name."
    """

    @staticmethod
    def _paint(label: str, rec_len: int, *, rec_start: int = 10):
        # Build a synthetic resite. `rec_start` defaults to 10 so a
        # long name's parens can extend a few cols left of the
        # recognition without falling off the chunk's left edge.
        f = {
            "type":   "resite",
            "start":  rec_start,
            "end":    rec_start + rec_len,
            "strand": 1,
            "color":  "white",
            "label":  label,
            "cut_col": None,
            "ext_cut_bp": None,
            "top_cut_bp": -1,
            "bottom_cut_bp": -1,
        }
        chunk_w = 60
        arr = [(" ", "")] * chunk_w
        sc._paint_feature_label(arr, f, 0, chunk_w)
        return "".join(ch for ch, _ in arr).rstrip()

    def test_5char_name_on_6bp_recognition_fits_full_name(self):
        # EcoRI = 5 chars on GAATTC = 6 bp. Pre-fix the parens
        # bracketed the 6-bp recognition exactly, leaving 4 cols
        # interior — too narrow for "EcoRI" (truncated to "EcoR").
        # Now paren width = max(rec_len=6, len(name)+2=7) = 7.
        rendered = self._paint("EcoRI", rec_len=6).lstrip()
        assert "EcoRI" in rendered, (
            f"full name must fit the label; got {rendered!r}"
        )
        assert rendered.startswith("("), rendered
        assert ")" in rendered, rendered

    def test_4char_name_on_6bp_recognition_uses_recognition_width(self):
        # BsaI = 4 chars on GGTCTC = 6 bp. Paren width =
        # max(6, 4+2=6) = 6 — same as recognition. Label centered.
        rendered = self._paint("BsaI", rec_len=6).lstrip()
        assert "BsaI" in rendered, rendered
        assert rendered.startswith("("), rendered
        assert ")" in rendered, rendered

    def test_paren_width_does_not_scale_with_recognition(self):
        # 4-char name on a HUGE recognition (10 bp) — paren width
        # SHOULD scale to recognition (so the parens still bracket
        # the recognition span as a position cue), label centered.
        rendered = self._paint("BsaI", rec_len=10).lstrip()
        assert "BsaI" in rendered, rendered
        # Per spec: paren width = max(rec_len=10, 4+2=6) = 10.
        assert rendered.startswith("(")
        assert ")" in rendered

    def test_long_name_widens_paren_past_recognition(self):
        # Hypothetical 8-char name on a 4-bp recognition. Paren
        # width = max(4, 8+2=10) = 10. The parens extend past the
        # recognition's right edge so the full name fits.
        rendered = self._paint("MyLongNm", rec_len=4).lstrip()
        assert "MyLongNm" in rendered, rendered


class TestResiteHighlightWrap:
    """Regression guards for the 2026-05-11 wrap-cut highlight fix.

    A Type IIS forward-strand cut whose recognition sits at the end of
    a circular sequence (raw `p + offset` exceeds `n`) lands at a small
    bp value post-modulo. The old `_resite_highlight_dict` then dragged
    `hi_start` down to the wrapped cut, producing a highlight that spanned
    most of the plasmid the wrong way. The new code wrap-encodes the
    highlight (`hi_end < hi_start` signals origin-wrap) and the renderer
    `_render_chunk` tests `i >= reh_s or i < reh_e` in that case.
    """

    @staticmethod
    def _highlight(seq, resite):
        class _MockMap:
            pass
        m = _MockMap()
        m._seq = seq
        return sc.SequencePanel._resite_highlight_dict(m, resite)

    def test_forward_type_iis_at_origin_wraps(self):
        """BsaI forward at p=24 on n=30 (cut wraps origin). The highlight
        must encode the wrap as hi_end < hi_start, NOT drag hi_start back."""
        resite = {
            "start": 24, "end": 30, "strand": 1, "color": "red",
            "label": "BsaI",
            "top_cut_bp": 1, "bottom_cut_bp": 5, "ext_cut_bp": 1,
        }
        hi = self._highlight("X" * 30, resite)
        assert hi["start"] == 24, (
            f"hi_start dragged to wrapped cut: got {hi['start']} "
            "(should stay at recognition start 24)"
        )
        # Wrap-encoded: hi_end < hi_start. Furthest-forward wrapped cut is
        # bot=5, so hi_end = max(wrap_cuts) = 5.
        assert hi["end"] == 5, f"expected wrap-encoded hi_end=5, got {hi['end']}"
        assert hi["end"] < hi["start"], "highlight must be wrap-encoded"

    def test_linear_type_iis_in_middle_extends_normally(self):
        """BsaI forward at p=10 on n=50 (cuts in middle, no wrap). The
        highlight must extend hi_end past the recognition to enclose
        the Type IIS cuts as before — wrap detection should NOT fire."""
        resite = {
            "start": 10, "end": 16, "strand": 1, "color": "red",
            "label": "BsaI",
            "top_cut_bp": 17, "bottom_cut_bp": 21, "ext_cut_bp": 17,
        }
        hi = self._highlight("X" * 50, resite)
        assert hi["start"] == 10
        assert hi["end"] == 21
        assert hi["end"] >= hi["start"], "must NOT wrap-encode normal case"

    def test_reverse_type_iis_genuine_upstream_extends_hi_start(self):
        """Reverse-strand Type IIS (BsaI bound on the bottom strand) cuts
        upstream of recognition in forward coords. With recognition in the
        middle of a linear seq, the cut is genuinely upstream (cut <
        hi_start) — hi_start must extend back, NOT wrap-encode."""
        resite = {
            "start": 20, "end": 26, "strand": -1, "color": "red",
            "label": "BsaI",
            "top_cut_bp": 15, "bottom_cut_bp": 19, "ext_cut_bp": 15,
        }
        hi = self._highlight("X" * 50, resite)
        assert hi["start"] == 15
        assert hi["end"] == 26
        assert hi["end"] >= hi["start"]

    def test_palindromic_cut_inside_recognition_unchanged(self):
        """EcoRI palindrome cuts inside recognition. Both cuts at bp 5
        for a site at [4, 10) — inside hi_start/hi_end so neither
        ext_cuts nor wrap_cuts gets populated. Highlight matches recog."""
        resite = {
            "start": 4, "end": 10, "strand": 1, "color": "red",
            "label": "EcoRI",
            "top_cut_bp": 5, "bottom_cut_bp": 9,
        }
        hi = self._highlight("X" * 30, resite)
        assert hi["start"] == 4
        assert hi["end"] == 10


class TestRestrictionScanLinearVsCircular:
    """Linear records must NOT scan past their end. Pre-2026-05-08
    every caller of `_scan_restriction_sites` defaulted to
    ``circular=True``, so an EcoRI site that only exists by joining
    a linear record's tail + head would still appear in the panel —
    biologically impossible because the linear ends don't ligate.
    """

    def test_circular_finds_origin_spanning_site(self):
        # Recognition GAATTC straddles the origin: 'GAA' at the
        # tail, 'TTC' at the head. Only valid on circular records.
        seq = "TTC" + "AAAAAA" + "GAA"   # 12 bp
        circular = sc._scan_restriction_sites(
            seq, min_recognition_len=6,
            unique_only=True, circular=True,
        )
        eco = [f for f in circular
               if f["type"] == "resite" and f["label"] == "EcoRI"]
        assert len(eco) == 1, (
            f"circular EcoRI scan should find the origin-spanning "
            f"site; got {len(eco)}"
        )
        # Wrap site emits two pieces — one tail, one head.
        assert eco[0]["start"] == 9 or eco[0]["end"] == 12

    def test_linear_does_not_find_origin_spanning_site(self):
        # Same sequence, but as a linear record — the recognition
        # only exists by joining bp 11 → bp 0, which a linear
        # molecule can't do. Scanner must return zero EcoRI hits.
        seq = "TTC" + "AAAAAA" + "GAA"   # 12 bp
        linear = sc._scan_restriction_sites(
            seq, min_recognition_len=6,
            unique_only=True, circular=False,
        )
        eco = [f for f in linear
               if f["type"] == "resite" and f["label"] == "EcoRI"]
        assert eco == [], (
            f"linear EcoRI scan must NOT report the origin-spanning "
            f"site; got {len(eco)} resites: {eco}"
        )

    def test_linear_in_body_site_still_found(self):
        # Sanity: a normal mid-record EcoRI site is still found in
        # linear mode — only the wrap-spanning one is suppressed.
        seq = "AAA" + "GAATTC" + "AAAA"  # 13 bp, EcoRI at pos 3
        linear = sc._scan_restriction_sites(
            seq, min_recognition_len=6,
            unique_only=True, circular=False,
        )
        eco = [f for f in linear
               if f["type"] == "resite" and f["label"] == "EcoRI"]
        assert len(eco) == 1
        assert eco[0]["start"] == 3 and eco[0]["end"] == 9


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

    def test_codon_start_2_skips_first_base(self):
        """Regression guard for 2026-05-11. `/codon_start=2` means the first
        base of the CDS is the last base of an incomplete leading codon —
        the real reading frame starts at offset 1. Pre-fix the qualifier
        was silently ignored, frame-shifting every AA past position 0."""
        # First base 'X' is the partial leading codon; real CDS is ATGGCATAG → M A *
        seq = "X" + "ATGGCATAG"
        aa = sc._translate_cds(seq, 0, len(seq), strand=1, codon_start=2)
        # First complete codon = ATG (M), then GCA (A), then TAG (*).
        # Note `?` may appear from the leading partial codon's bases
        # depending on the slice; check the in-frame portion is right.
        assert "MA*" in aa, f"expected MA* in {aa!r}"

    def test_codon_start_3_skips_first_two_bases(self):
        seq = "XX" + "ATGGCATAG"     # /codon_start=3
        aa = sc._translate_cds(seq, 0, len(seq), strand=1, codon_start=3)
        assert "MA*" in aa, f"expected MA* in {aa!r}"

    def test_codon_start_matches_biopython(self):
        """Cross-check against Biopython's feature.extract().translate(cds=False)
        which DOES honour /codon_start when reading from a SeqRecord. We
        replicate via Seq.translate on a manually-offset slice."""
        from Bio.Seq import Seq
        for cs in (1, 2, 3):
            seq = "X" * (cs - 1) + "ATGAAATTTGGGCCCTAG"
            sc_aa = sc._translate_cds(seq, 0, len(seq), strand=1, codon_start=cs)
            bp_aa = str(Seq(seq[cs - 1:]).translate())
            assert sc_aa.rstrip("*") == bp_aa.rstrip("*"), (
                f"codon_start={cs}: SC={sc_aa!r} vs BioPython={bp_aa!r}"
            )

    def test_codon_start_reverse_strand(self):
        """Reverse-strand CDS with /codon_start=2: the offset is applied
        after RC, so the trailing base of the RC'd sequence is the partial
        leading codon and the in-frame translation starts one base in."""
        from Bio.Seq import Seq
        fwd = "ATGAAATAG"    # M K *
        rc = str(Seq(fwd).reverse_complement())
        # Reverse-strand CDS spanning the RC'd region + 1 leading 'X' base on
        # the forward template (which becomes a TRAILING base after RC) —
        # codon_start=1 sees a frame-shifted CDS, codon_start=2 corrects it.
        full = rc + "X"
        # With codon_start=1 the trailing X drops into the frame and breaks
        # the protein; with codon_start=2 the X is skipped and we get back MK*.
        cs2 = sc._translate_cds(full, 0, len(full), strand=-1, codon_start=2)
        assert "MK*" in cs2, f"expected MK* in {cs2!r}"

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


class TestWrapCDSInlineTranslation:
    """Regression guard for 2026-04-30: the inline AA-translation row
    drawn into the seq panel's lane art must place AA letters at the
    correct codon midpoints AND with the correct one-letter codes for
    wrap-CDS features. Pre-fix `_paint_cds_aa` ran codon math against
    the SPLIT virtual half (start=0 for the head), which both painted
    at the wrong bp and translated the wrong reading frame.

    The fix: `_feats_in_chunk` stamps `_orig_start` / `_orig_end` on
    each split half, and `_paint_cds_aa` / `_cds_aa_list` use those
    for the codon-midpoint formula and translation source so every
    half of a wrap-CDS shares one canonical translation.
    """

    def test_wrap_cds_paints_correct_letters_and_positions(self):
        import random
        random.seed(42)
        n = 1000
        seq = "".join(random.choice("ACGT") for _ in range(n))
        # Wrap CDS: tail [900, 1000) + head [0, 30); 130 bp = 43 codons.
        feats = [{"start": 900, "end": 30, "type": "CDS", "strand": 1,
                  "label": "wrapCDS", "color": "cyan"}]

        # Head chunk [0, 60) — should host codons 33..52 (codon 33's
        # midpoint at bp 0, codon 52's at bp 57).
        chunk_start, chunk_end = 0, 60
        in_chunk = sc._feats_in_chunk(feats, chunk_start, chunk_end, n)
        # Head half should carry the original wrap coords.
        head = next(f for f in in_chunk if f.get("type") == "CDS")
        assert head["_orig_start"] == 900
        assert head["_orig_end"] == 30

        arr = [(" ", "")] * (chunk_end - chunk_start)
        sc._paint_cds_aa(arr, head, chunk_start, chunk_end,
                          seq.upper(), None)
        painted = [(i, ch) for i, (ch, _sty) in enumerate(arr) if ch != " "]

        # Compute the canonical answer from the joined CDS sequence.
        cds_seq = (seq[900:] + seq[:30]).upper()
        expected = []
        for i in range(len(cds_seq) // 3):
            mid = (900 + 3 * i + 1) % n
            if 0 <= mid < 60:
                aa = sc._CODON_TABLE.get(cds_seq[3 * i:3 * i + 3], "?")
                expected.append((mid, aa))
        assert painted == expected, (
            f"Wrap-CDS AA letters mis-painted in head half.\n"
            f"painted={painted}\nexpected={expected}"
        )

    def test_wrap_cds_aa_cache_keys_on_orig_coords(self):
        """Both halves of a wrap-CDS should share one cached translation
        — keyed on `(_orig_start, _orig_end)` — so we don't translate
        the same protein twice per render. Pre-fix the cache keyed on
        the half's local `(start, end)` and produced two stale entries
        per wrap-CDS."""
        seq = "ATGAAATGCAAAAAAAAAACCCTAA" + "G" * 75   # 100 bp
        n = len(seq)
        # Wrap CDS [80, 25); covers tail [80, 100) + head [0, 25).
        feats = [{"start": 80, "end": 25, "type": "CDS", "strand": 1,
                  "label": "f", "color": "white"}]
        head = sc._feats_in_chunk(feats, 0, 30, n)[0]
        tail = sc._feats_in_chunk(feats, 80, 100, n)[0]
        sc._CDS_AA_CACHE.clear()
        aa_head, _, _ = sc._cds_aa_list(seq, head)
        aa_tail, _, _ = sc._cds_aa_list(seq, tail)
        assert aa_head is aa_tail   # same cached list
        assert len(sc._CDS_AA_CACHE) == 1

    def test_non_wrap_cds_unchanged(self):
        """Sanity: non-wrap CDS features (the common case) must paint
        identically with or without the `_orig_*` keys — i.e. the
        fallback `f.get("_orig_start", f["start"])` path is correct."""
        seq = "ATGAAATGCCCCAAATGCAAA" + "G" * 79   # 100 bp
        feats = [{"start": 0, "end": 21, "type": "CDS", "strand": 1,
                  "label": "f", "color": "white"}]
        n = len(seq)
        chunk_start, chunk_end = 0, 30
        in_chunk = sc._feats_in_chunk(feats, chunk_start, chunk_end, n)
        cds = next(f for f in in_chunk if f.get("type") == "CDS")
        # Non-wrap: no `_orig_*` keys stamped.
        assert "_orig_start" not in cds

        arr = [(" ", "")] * (chunk_end - chunk_start)
        sc._paint_cds_aa(arr, cds, chunk_start, chunk_end,
                          seq.upper(), None)
        painted = [(i, ch) for i, (ch, _sty) in enumerate(arr) if ch != " "]
        # Codon midpoints at bp 1, 4, 7, ..., 19 (7 codons in 21 bp).
        cds_seq = seq[0:21].upper()
        expected = [
            (3 * i + 1, sc._CODON_TABLE.get(cds_seq[3 * i:3 * i + 3], "?"))
            for i in range(7)
        ]
        assert painted == expected


# ═══════════════════════════════════════════════════════════════════════════════
# ORF finder — six-frame open-reading-frame scan, wrap-aware
# ═══════════════════════════════════════════════════════════════════════════════

class TestFindOrfs:
    """`_find_orfs` six-frame scan. Standard table 1, ATG-only by
    default. Wrap-aware on circular plasmids: ORFs that cross the
    origin are reported with `end < start` matching the wrap-feature
    convention used elsewhere."""

    def test_single_forward_orf(self):
        """ATG…TAA with 30 GCC codons in between → 1 ORF (+) strand.
        `length_aa` excludes the stop, so M + 30 GCC = 31 residues."""
        seq = "AA" + "ATG" + "GCC" * 30 + "TAA" + "AA"
        orfs = sc._find_orfs(seq, min_aa=30, circular=False)
        assert len(orfs) == 1
        o = orfs[0]
        assert o["strand"] == 1
        assert o["length_aa"] == 31
        assert o["start"] == 2
        assert o["end"] == 2 + 3 + 30 * 3 + 3
        # AA seq starts with M and ends with *
        assert o["aa_seq"].startswith("M")
        assert o["aa_seq"].endswith("*")

    def test_min_length_filter(self):
        """An ORF below `min_aa` must be filtered out."""
        # M + 19 GCC = 20 coded residues (excluding stop).
        seq = "ATG" + "GCC" * 19 + "TAA"
        assert sc._find_orfs(seq, min_aa=30, circular=False) == []
        assert len(sc._find_orfs(seq, min_aa=20, circular=False)) == 1

    def test_reverse_strand_orf(self):
        """An ORF on the (-) strand: forward seq carries the RC of the
        ATG…TAA pattern. The reported forward-strand coords cover the
        recognition span on the top strand even though the ORF itself
        reads on the bottom."""
        body = "ATG" + "GCC" * 30 + "TAA"
        rc_body = sc._rc(body)
        seq = "AAA" + rc_body + "AAA"
        orfs = sc._find_orfs(seq, min_aa=30, circular=False)
        assert len(orfs) == 1
        o = orfs[0]
        assert o["strand"] == -1
        assert o["length_aa"] == 31  # M + 30 GCC, excluding stop
        # Forward coords cover [3, 3 + len(body)) since rc_body starts at 3.
        assert o["start"] == 3
        assert o["end"] == 3 + len(body)

    def test_alt_starts_off_by_default(self):
        """GTG / TTG must NOT start an ORF unless `include_alt_starts=True`."""
        seq = "GTG" + "GCC" * 30 + "TAA"
        assert sc._find_orfs(seq, min_aa=30, circular=False) == []
        orfs = sc._find_orfs(seq, min_aa=30, circular=False,
                              include_alt_starts=True)
        assert len(orfs) == 1
        # The AA seq from a GTG start translates as Val (or Met by the
        # alt-start convention; we use the canonical table so it's V).
        assert orfs[0]["aa_seq"][0] in {"V", "M"}

    def test_circular_wrap_orf(self):
        """An ORF that crosses the origin on a circular plasmid is
        reported with `end < start`. Reproduce the case by placing
        the start codon near the right edge so the body wraps."""
        # n = 120; place ATG at bp 110, body+stop wraps to bp 110 + 3 + 90 + 3 - 120 = 86.
        body = "ATG" + "GCC" * 30 + "TAA"   # 99 bp
        n = 120
        prefix_len = 110
        wrapped = body[(n - prefix_len):]   # the part that goes back to start
        leading = body[:n - prefix_len]
        seq = "A" * prefix_len + leading + "G" * (n - prefix_len - len(wrapped))
        # Hand-build a clean circular layout instead — the splicing above
        # is fiddly. Use seq = filler[0:110] + "ATG..." with the body
        # split at bp 120.
        filler = "A" * 110
        seq = filler + body[:10]            # bp 110..120 holds first 10 bp of body
        seq = seq + body[10:]               # then the rest hangs off the end
        # That's 110 + len(body) = 209 bp. Now make the SEQUENCE itself n=120
        # by taking the first 120 chars and stuffing the tail at the front.
        full_with_tail = filler + body
        n = 120
        head = full_with_tail[n:]           # bp the wrap should cover at the start
        seq_circ = head + filler[len(head):] + body[:n - 110]
        # `seq_circ` is 120 bp with ATG at bp 110 and the body wrapping to bp `len(head)`.
        assert len(seq_circ) == 120
        orfs = sc._find_orfs(seq_circ, min_aa=30, circular=True)
        # At least one wrap ORF expected; some non-wrap junk ORFs may
        # also appear depending on filler — pick the wrap one and
        # verify its shape.
        wrap_orfs = [o for o in orfs if o["end"] < o["start"]]
        assert wrap_orfs, "expected at least one wrap ORF"
        o = wrap_orfs[0]
        assert o["strand"] == 1
        assert o["length_aa"] == 31  # M + 30 GCC, excluding stop

    def test_too_short_sequence(self):
        """A 5 bp seq has no room for any ORF."""
        assert sc._find_orfs("ACGTA", circular=False) == []

    def test_dedupe_no_duplicate_in_doubled_scan(self):
        """A circular ORF that does NOT wrap should not be reported twice
        even though the doubled-scan visits the same region twice."""
        body = "ATG" + "GCC" * 30 + "TAA"     # 99 bp
        seq = body + "G" * 21                  # 120 bp circular plasmid
        orfs = sc._find_orfs(seq, min_aa=30, circular=True)
        # The ATG at position 0 should map to exactly one ORF.
        zero_starts = [o for o in orfs if o["start"] == 0 and o["strand"] == 1]
        assert len(zero_starts) == 1


# ═══════════════════════════════════════════════════════════════════════════════
# Annotation transfer — exact-match feature propagation between plasmids
# ═══════════════════════════════════════════════════════════════════════════════

class TestAnnotationTransfer:
    """`_find_annotation_transfers` matches features by sequence
    identity and reports their target-strand coords. Exact-match for
    v1.0; wrap-aware on circular targets."""

    @staticmethod
    def _rec(seq: str, feats: list[tuple] = (), *,  # type: ignore[assignment]
              circular: bool = True):
        """Tiny SeqRecord builder. `feats` is `[(start, end, strand,
        type, label), ...]`."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq      import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        rec = SeqRecord(Seq(seq), id="t", name="t")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular" if circular else "linear"
        for s, e, strand, ftype, label in feats:
            n = len(seq)
            if e < s:
                # Wrap-feature CompoundLocation: [s, n) + [0, e).
                loc = CompoundLocation([
                    FeatureLocation(s, n, strand=strand),
                    FeatureLocation(0, e, strand=strand),
                ])
            else:
                loc = FeatureLocation(s, e, strand=strand)
            rec.features.append(
                SeqFeature(loc, type=ftype,
                            qualifiers={"label": [label]})
            )
        return rec

    def test_exact_forward_match(self):
        """A 60 bp CDS in source A should transfer onto an identical
        substring in target B."""
        body = "ATG" + "GCC" * 19 + "TAA"   # 60 bp
        prefix = "AAAAAAAAAA"                # 10 bp
        suffix = "AAAAAAAAAA"                # 10 bp
        src = self._rec(prefix + body + suffix,
                         feats=[(10, 70, 1, "CDS", "myCDS")])
        tgt = self._rec("TTTTT" + body + "TTTTT")
        out = sc._find_annotation_transfers(src, tgt, min_len=30)
        assert len(out) == 1
        t = out[0]
        assert t["label"] == "myCDS"
        assert t["type"] == "CDS"
        assert t["target_start"] == 5
        assert t["target_end"] == 5 + 60
        assert t["target_strand"] == 1
        assert t["length"] == 60

    def test_reverse_complement_match(self):
        """A feature whose RC appears in the target is reported as
        a reverse-strand match. Target carries the RC of the
        feature's coding-strand bases; the transfer flips strand and
        reports the forward-coord span where the RC sits."""
        body = "ATG" + "GCC" * 19 + "TAA"   # 63 bp
        src = self._rec("AA" + body + "AA",
                         feats=[(2, 62, 1, "CDS", "fwdCDS")])
        # Feature bases = body[0:60] (the slice [2, 62) of "AA"+body+"AA").
        feat_bases = body[0:60]
        # Target seq = "TT" + rc(feat_bases) + "TT" so the only hit is
        # an RC match. n_tgt = 64; rc(feat_bases) sits at [2, 62) on
        # the forward strand. The bottom-strand feature reads 5'→3'
        # going right-to-left, so the forward span is [n - 62, n - 2)
        # = [2, 62).
        tgt = self._rec("TT" + sc._rc(feat_bases) + "TT")
        out = sc._find_annotation_transfers(src, tgt, min_len=30)
        assert len(out) == 1
        t = out[0]
        assert t["target_strand"] == -1
        assert t["target_start"] == 2
        assert t["target_end"] == 62

    def test_skip_features_below_min_len(self):
        """Short features (e.g. primer-binding sites) generate noise;
        scan must skip them at `min_len`."""
        src = self._rec("AAAAATGCATG" + "AAA" * 30,
                         feats=[(5, 11, 1, "primer_bind", "tiny")])
        tgt = self._rec("AAAAATGCATG" + "AAA" * 30)
        out = sc._find_annotation_transfers(src, tgt, min_len=30)
        assert out == []

    def test_circular_wrap_target(self):
        """A feature whose bases cross the origin in the target must
        be reported with `end < start` (wrap convention)."""
        body = "ATG" + "GCC" * 19 + "TAA"   # 63 bp
        # Source feature spans [2, 62) of "AA"+body+"AA"; feat_bases is
        # body[0:60].
        src = self._rec("AA" + body + "AA",
                         feats=[(2, 62, 1, "CDS", "spanCDS")])
        feat_bases = body[0:60]
        # Build a 80-bp circular target where feat_bases occupies
        # forward-strand positions [70, 80) ∪ [0, 50). target_start
        # should be 70, target_end 50 (wrap).
        tgt_seq = feat_bases[10:60] + ("G" * 20) + feat_bases[0:10]
        assert len(tgt_seq) == 80
        tgt = self._rec(tgt_seq, circular=True)
        out = sc._find_annotation_transfers(src, tgt, min_len=30)
        assert len(out) == 1
        t = out[0]
        assert t["target_strand"] == 1
        assert t["target_start"] == 70
        assert t["target_end"] == 50

    def test_no_match_returns_empty(self):
        body = "ATG" + "GCC" * 19 + "TAA"
        src = self._rec("AA" + body + "AA",
                         feats=[(2, 62, 1, "CDS", "alpha")])
        tgt = self._rec("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
        assert sc._find_annotation_transfers(src, tgt, min_len=30) == []

    def test_skip_source_feature(self):
        """`source` features (the GenBank metadata feature spanning the
        whole record) must not be transferred."""
        body = "AAAAAAAA" * 10
        src = self._rec(body,
                         feats=[(0, len(body), 1, "source", "src")])
        tgt = self._rec(body)
        out = sc._find_annotation_transfers(src, tgt, min_len=30)
        assert out == []

    def test_whole_plasmid_match_emits_one_transfer(self):
        """Regression guard for 2026-05-06: when `feat_len == n_tgt`
        the wrap-fold collapsed `t_e` to `t_s` and the dedupe key
        aliased every full match. Now special-cased to a single
        full-record transfer at `[0, n_tgt)`."""
        body = "ATG" + "GCC" * 30 + "TAA"   # 99 bp
        src = self._rec(body, feats=[(0, len(body), 1, "CDS", "full")])
        tgt = self._rec(body)
        out = sc._find_annotation_transfers(src, tgt, min_len=30)
        # At least one match (forward); shape should be a non-wrap
        # full-record span so the apply path emits a normal
        # FeatureLocation rather than a degenerate CompoundLocation.
        full_matches = [t for t in out
                        if t["target_start"] == 0
                        and t["target_end"] == len(body)
                        and t["target_strand"] == 1]
        assert len(full_matches) == 1


class TestCustomEnzymeListFilter:
    """`allowed_enzymes` parameter on `_scan_restriction_sites` (GH #13,
    Cory Tobin 2026-05-14). When supplied, the scan restricts to JUST
    those enzymes — overrides `min_recognition_len` and `unique_only`
    so a hand-picked list always shows in full, regardless of cut
    count or recognition length.
    """

    def test_allow_list_overrides_min_len(self):
        """A 4-cutter like Sau3AI (GATC) wouldn't appear with the
        default `min_recognition_len=6`, but should when the user
        explicitly picks it."""
        seq = "AAAGATCAAAATCGAAAAAGATCAAA"
        # Default scan with min_recognition_len=6: Sau3AI excluded
        default = sc._scan_restriction_sites(
            seq, min_recognition_len=6, unique_only=False, circular=False,
        )
        assert not any(h.get("label") == "Sau3AI" for h in default)
        # Allow-list with Sau3AI: hits surface
        allowed = sc._scan_restriction_sites(
            seq, circular=False,
            allowed_enzymes=frozenset({"Sau3AI"}),
        )
        sau_hits = [h for h in allowed
                    if h.get("label") == "Sau3AI"
                    and h.get("type") == "resite"]
        assert len(sau_hits) >= 1

    def test_allow_list_overrides_unique_only(self):
        """A repeat-laden sequence that breaks the unique-cutter
        filter should still show every hit when the user explicitly
        picks the enzyme. Pre-fix the `unique_only=True` filter
        applied even with `allowed_enzymes` set, hiding multi-cutters
        from the user's hand-picked list — the opposite of what the
        user actually wants."""
        # GAATTC × 100 — way too many EcoRI sites to be "unique"
        seq = "GAATTCAAAA" * 100
        default = sc._scan_restriction_sites(
            seq, unique_only=True, circular=False,
        )
        assert not any(h.get("label") == "EcoRI" for h in default)
        allowed = sc._scan_restriction_sites(
            seq, unique_only=True, circular=False,
            allowed_enzymes=frozenset({"EcoRI"}),
        )
        ecori_hits = [h for h in allowed
                       if h.get("label") == "EcoRI"
                       and h.get("type") == "resite"]
        assert len(ecori_hits) == 100

    def test_allow_list_excludes_other_enzymes(self):
        """Only enzymes in the allow-list should appear. A site like
        BamHI shouldn't surface if only EcoRI is allowed, even when
        both recognition sequences are present."""
        seq = "GAATTCAAAAGGATCCAAAAGAATTCAAAA"
        allowed = sc._scan_restriction_sites(
            seq, circular=False,
            allowed_enzymes=frozenset({"EcoRI"}),
        )
        labels = {h.get("label") for h in allowed if h.get("label")}
        assert "EcoRI" in labels
        assert "BamHI" not in labels

    def test_unknown_enzyme_in_allow_list_silently_dropped(self):
        """A typo or HF-variant rename in the allow-list shouldn't
        crash the scan — unknown names just don't match anything."""
        seq = "GAATTCAAAA"
        result = sc._scan_restriction_sites(
            seq, circular=False,
            allowed_enzymes=frozenset({"EcoRI", "BogusEnzyme123"}),
        )
        labels = {h.get("label") for h in result if h.get("label")}
        assert labels == {"EcoRI"}

    def test_empty_allow_list_returns_nothing(self):
        """An empty frozenset means "no enzymes" — not "use defaults".
        The caller is responsible for passing None when they mean
        defaults; the dispatch in PlasmidApp handles that."""
        seq = "GAATTCAAAAGGATCC"
        result = sc._scan_restriction_sites(
            seq, circular=False, allowed_enzymes=frozenset(),
        )
        assert result == []

    def test_settings_validator_canonicalises_csv(self):
        """`_settings_validator_custom_enzymes_csv` drops unknown
        names, dedupes, sorts, and produces a canonical CSV that
        survives a settings.json round-trip."""
        v_fn = sc._settings_validator_custom_enzymes_csv
        result, err = v_fn("BamHI, EcoRI, BamHI, BogusEnzyme, BsaI")
        assert err is None
        assert result == "BamHI,BsaI,EcoRI"

    def test_settings_validator_handles_empty_string(self):
        v_fn = sc._settings_validator_custom_enzymes_csv
        result, err = v_fn("")
        assert err is None
        assert result == ""

    def test_settings_validator_rejects_non_string(self):
        v_fn = sc._settings_validator_custom_enzymes_csv
        result, err = v_fn(42)
        assert result is None
        assert err is not None


class TestTranslTableNonStandard:
    """Sweep #30 (2026-05-28): a CDS carrying /transl_table must translate
    with the named NCBI genetic code, not the hardcoded standard table —
    otherwise a mito / Mycoplasma CDS renders the wrong protein AND a
    reassigned stop (TGA→Trp) trips a false premature-stop ⚠."""

    def test_standard_table_is_canonical_object(self):
        assert sc._codon_table_for(1) is sc._CODON_TABLE
        assert sc._codon_table_for(None) is sc._CODON_TABLE
        assert sc._codon_table_for(0) is sc._CODON_TABLE   # falsy → standard

    def test_table4_reassigns_tga_to_trp(self):
        m = sc._codon_table_for(4)   # mold/protozoan mito + Mycoplasma
        assert m["TGA"] == "W"       # the reassignment that matters
        assert m["TAA"] == "*"       # other stops unchanged
        assert m["TAG"] == "*"
        assert m["ATG"] == "M"

    def test_vertebrate_mito_table2(self):
        m = sc._codon_table_for(2)
        assert m["TGA"] == "W"       # Trp
        assert m["AGA"] == "*"       # AGR = stop in vertebrate mito
        assert m["AGG"] == "*"
        assert m["ATA"] == "M"       # Met (sense in table 2)

    def test_unknown_table_falls_back_to_standard(self):
        # A hand-edited /transl_table=99 must not crash — fall back to std.
        m = sc._codon_table_for(99)
        assert m["TGA"] == "*"

    def test_translate_cds_honours_table4(self):
        seq = "ATGTGATAA"            # ATG TGA TAA
        assert sc._translate_cds(seq, 0, 9, 1) == "M**"
        assert sc._translate_cds(seq, 0, 9, 1, transl_table=4) == "MW*"

    def test_cds_aa_list_internal_tga_not_stop_under_table4(self):
        # ATG TGA AAA TAA — the internal TGA is a real premature stop under
        # the standard code but a Trp under table 4 (so no false ⚠).
        seq = "ATGTGAAAATAA"
        aa_std, _len, _ve = sc._cds_aa_list(
            seq, {"start": 0, "end": 12, "strand": 1})
        assert aa_std == ["M", "*", "K", "*"]
        aa_t4, _len, _ve = sc._cds_aa_list(
            seq, {"start": 0, "end": 12, "strand": 1, "transl_table": 4})
        assert aa_t4 == ["M", "W", "K", "*"]
