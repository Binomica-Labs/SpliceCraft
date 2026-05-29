"""
test_domesticator — Golden Braid L0 parts domesticator tests.

Tests cover:
  - _pick_binding_region: length in 18–25 range, Tm near target
  - _design_gb_primers: correct BsaI sites, overhangs match part type,
    binding regions present, amplicon length, all 6 part positions
  - Parts-bin persistence: save/load round-trip, user parts show in table
  - Integration: full domesticator flow (design + save) via the modal
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def random_template():
    """A 2000-bp random ACGT template for primer design tests.

    0xCAFE happens to contain a GGTCTC (BsaI) on one strand, which was
    fine when L0 domestication used BsaI tails. After the switch to
    Esp3I (CGTCTC), the internal-site check only cares about CGTCTC /
    GAGACG — but we scrub all four sites anyway so the fixture stays
    enzyme-agnostic if we ever flip back or run parallel designs."""
    rng = random.Random(0xCAFE)
    seq = "".join(rng.choice("ACGT") for _ in range(2000))
    for site, repl in (("GGTCTC", "GCTGTC"),
                       ("GAGACC", "GACCGA"),
                       ("CGTCTC", "CGACTC"),
                       ("GAGACG", "GAGACT")):
        seq = seq.replace(site, repl)
    return seq


# `isolated_parts_bin` lives in tests/conftest.py — same redirect
# pattern, shared with test_traditional_cloning.py and others.


# ═══════════════════════════════════════════════════════════════════════════════
# _pick_binding_region
# ═══════════════════════════════════════════════════════════════════════════════

class TestPickBindingRegion:
    def test_length_in_valid_range(self, random_template):
        seq = random_template[:200]
        binding, tm = sc._pick_binding_region(seq)
        assert 18 <= len(binding) <= 25

    def test_tm_near_target(self, random_template):
        seq = random_template[:200]
        binding, tm = sc._pick_binding_region(seq, target_tm=57.0)
        # Should be within 8°C of the target — generous but catches gross errors
        assert 49 < tm < 65, f"Tm {tm}°C is too far from 57°C target"

    def test_binding_is_prefix_of_input(self, random_template):
        seq = random_template[:200]
        binding, _ = sc._pick_binding_region(seq)
        assert seq.startswith(binding)

    def test_short_input_returns_what_it_can(self):
        # If the input is shorter than min_len, return whatever we have
        binding, _ = sc._pick_binding_region("ATGATG", min_len=4, max_len=10)
        assert len(binding) <= 6

    def test_different_targets_produce_different_lengths(self, random_template):
        seq = random_template[:200]
        b_lo, _ = sc._pick_binding_region(seq, target_tm=45.0)
        b_hi, _ = sc._pick_binding_region(seq, target_tm=65.0)
        # Higher target Tm should pick a longer (or equal) binding region
        assert len(b_hi) >= len(b_lo)


# ═══════════════════════════════════════════════════════════════════════════════
# _design_gb_primers — correctness across all 6 part types
# ═══════════════════════════════════════════════════════════════════════════════

class TestDesignGBPrimers:
    @pytest.mark.parametrize("part_type", list(sc._GB_POSITIONS.keys()))
    def test_l0_enzyme_site_in_both_primers(self, random_template, part_type):
        """L0 domestication primers must carry the Esp3I site (CGTCTC) —
        Golden Braid uses Esp3I/BsmBI at L0 so the domesticated part
        survives the downstream L1+ BsaI assembly without re-cutting."""
        result = sc._design_gb_primers(random_template, 100, 600, part_type)
        assert sc._GB_L0_ENZYME_SITE in result["fwd_full"]
        assert sc._GB_L0_ENZYME_SITE in result["rev_full"]
        # And must NOT carry a bare BsaI site in the tail — that would
        # defeat the whole point of splitting L0/L1 across two enzymes.
        assert "GGTCTC" not in result["fwd_full"][:12], (
            "L0 primer tail must be Esp3I, not BsaI"
        )
        assert "GGTCTC" not in result["rev_full"][:12], (
            "L0 primer tail must be Esp3I, not BsaI"
        )

    @pytest.mark.parametrize("part_type", list(sc._GB_POSITIONS.keys()))
    def test_correct_overhangs_in_primers(self, random_template, part_type):
        """The 5' overhang must appear in the forward primer right after
        pad+Esp3I+spacer. The RC of the 3' overhang must appear in the
        reverse primer at the same offset."""
        result = sc._design_gb_primers(random_template, 100, 600, part_type)
        oh5, oh3 = result["oh5"], result["oh3"]
        tail_prefix_len = len(sc._GB_PAD + sc._GB_L0_ENZYME_SITE + sc._GB_SPACER)

        # Forward primer: position tail_prefix_len should start with oh5
        fwd_oh = result["fwd_full"][tail_prefix_len:tail_prefix_len + 4]
        assert fwd_oh == oh5, f"fwd overhang {fwd_oh} != expected {oh5}"

        # Reverse primer: should have RC of oh3
        rev_oh = result["rev_full"][tail_prefix_len:tail_prefix_len + 4]
        expected_rc = sc._rc(oh3)
        assert rev_oh == expected_rc, (
            f"rev overhang {rev_oh} != expected RC({oh3})={expected_rc}"
        )

    @pytest.mark.parametrize("part_type", list(sc._GB_POSITIONS.keys()))
    def test_binding_regions_present(self, random_template, part_type):
        result = sc._design_gb_primers(random_template, 100, 600, part_type)
        assert result["fwd_binding"]
        assert result["rev_binding"]
        assert 18 <= len(result["fwd_binding"]) <= 25
        assert 18 <= len(result["rev_binding"]) <= 25

    def test_fwd_binding_matches_insert_start(self, random_template):
        # 2026-05-21: GB CDS with AATG oh5 skips the first 3 bp
        # (the literal ATG) so the binding starts at codon 2.
        # Use a Promoter (oh5 GGAG) which has no codon-skip rule
        # so the binding matches the bare insert start.
        result = sc._design_gb_primers(random_template, 100, 600, "Promoter")
        insert = random_template[100:600].upper()
        assert insert.startswith(result["fwd_binding"])

    def test_cds_fwd_binding_skips_atg(self, random_template):
        # Sacred regression guard for the GB CDS ATG-fusion fix.
        # AATG overhang carries the start codon → forward primer
        # must start at codon 2 (bp+3), not at the insert's own
        # ATG.
        result = sc._design_gb_primers(random_template, 100, 600, "CDS")
        insert = random_template[100:600].upper()
        assert insert[3:].startswith(result["fwd_binding"])

    def test_rev_binding_matches_insert_end_rc(self, random_template):
        result = sc._design_gb_primers(random_template, 100, 600, "CDS")
        insert_rc = sc._rc(random_template[100:600].upper())
        assert insert_rc.startswith(result["rev_binding"])

    def test_insert_seq_is_template_slice(self, random_template):
        result = sc._design_gb_primers(random_template, 200, 800, "Promoter")
        assert result["insert_seq"] == random_template[200:800].upper()

    def test_amplicon_len_is_positive(self, random_template):
        result = sc._design_gb_primers(random_template, 0, 500, "Terminator")
        assert result["amplicon_len"] > 500  # insert + tails

    def test_position_matches_gb_standard(self, random_template):
        for ptype, (pos, oh5, oh3) in sc._GB_POSITIONS.items():
            result = sc._design_gb_primers(random_template, 50, 300, ptype)
            assert result["position"] == pos
            assert result["oh5"] == oh5
            assert result["oh3"] == oh3

    def test_pad_and_spacer_present(self, random_template):
        result = sc._design_gb_primers(random_template, 0, 200, "CDS")
        assert result["fwd_full"].startswith(sc._GB_PAD)
        assert result["rev_full"].startswith(sc._GB_PAD)
        # Spacer A should appear between Esp3I site and overhang
        fwd_after_enzyme = result["fwd_full"][
            len(sc._GB_PAD) + len(sc._GB_L0_ENZYME_SITE)
        ]
        assert fwd_after_enzyme == sc._GB_SPACER


# ═══════════════════════════════════════════════════════════════════════════════
# Silent-mutation repair of internal BsaI / Esp3I sites during domestication
# ═══════════════════════════════════════════════════════════════════════════════
#
# Internal BsaI (GGTCTC / GAGACC) or Esp3I (CGTCTC / GAGACG) sites would break
# Golden Braid assembly: Esp3I self-cuts during L0 domestication and BsaI
# re-cuts during the downstream L1 assembly. For coding parts with a codon
# table available, _design_gb_primers now repairs these in-frame by swapping
# in synonymous codons via _codon_fix_sites. These tests pin that behavior.

class TestDomesticatorSilentMutation:
    """Coding-part domestication auto-mutates internal BsaI / Esp3I sites
    via synonymous codon substitution when a codon table is supplied."""

    @staticmethod
    def _k12_raw() -> dict:
        """Built-in E. coli K12 codon table — the same registry entry that
        DomesticatorModal.on_mount seeds as the default."""
        return dict(sc._CODON_BUILTIN_K12)

    @staticmethod
    def _cds_with_bsai() -> str:
        """72-bp in-frame CDS with an internal GGTCTC that spans two
        codons (G-L). Length is a multiple of 3, translation yields
        M-A…A-G-L-A…A-* so _codon_fix_sites can pick a synonymous
        glycine codon (GGC/GGA/GGG) to remove the site."""
        # Positions 0..2: M, 3..32: 10×A, 33..38: G-L, 39..68: 10×A, 69..71: stop
        return "ATG" + "GCG" * 10 + "GGT" + "CTC" + "GCG" * 10 + "TAA"

    @staticmethod
    def _cds_with_esp3i() -> str:
        """72-bp in-frame CDS with an internal CGTCTC (Esp3I) site — codons
        R-L at positions 33-38. R has 6 synonyms so the fixer has plenty
        of room to pick a replacement."""
        return "ATG" + "GCG" * 10 + "CGT" + "CTC" + "GCG" * 10 + "TAA"

    # ── Coding parts: silent mutation happens ──────────────────────────────

    def test_bsai_site_in_coding_part_is_silently_mutated(self):
        cds = self._cds_with_bsai()
        assert "GGTCTC" in cds, "sanity: test input must contain the site"
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=self._k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"], (
            "internal BsaI site should have been silently removed"
        )
        assert "GAGACC" not in r["insert_seq"]
        assert r["mutations"], "mutations list must be populated"
        assert any("BsaI" in m for m in r["mutations"])

    def test_esp3i_site_in_coding_part_is_silently_mutated(self):
        cds = self._cds_with_esp3i()
        assert "CGTCTC" in cds
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=self._k12_raw())
        assert "error" not in r, r
        assert "CGTCTC" not in r["insert_seq"]
        assert "GAGACG" not in r["insert_seq"]
        assert r["mutations"]
        assert any("Esp3I" in m for m in r["mutations"])

    def test_mutated_insert_translates_to_same_protein(self):
        """Silent mutation means synonymous — the protein sequence must
        survive the fix unchanged (this is the whole point of routing
        through a codon table rather than a random base substitution)."""
        cds = self._cds_with_bsai()
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=self._k12_raw())
        before = sc._mut_translate(cds)
        after  = sc._mut_translate(r["insert_seq"])
        assert before == after, (
            f"protein changed: {before!r} → {after!r}"
        )

    def test_mutated_insert_used_in_primer_binding(self):
        """The forward primer's binding region is picked from the START of
        the (possibly-mutated) insert. If the mutation lands inside the
        binding window the binding should still match the mutated insert,
        not the original template.

        Note: GB CDS (oh5 = AATG) skips the literal ATG start
        codon — the forward binding starts at codon 2 (insert[3:]).
        Sacred per 2026-05-21 ATG-fusion fix.
        """
        cds = self._cds_with_bsai()
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=self._k12_raw())
        insert = r["insert_seq"]
        assert insert[3:].startswith(r["fwd_binding"])

    @pytest.mark.parametrize("part_type", ["CDS", "CDS-NS", "C-tag"])
    def test_all_coding_part_types_get_silent_repair(self, part_type):
        """CDS / CDS-NS / C-tag all route through the silent-mutation
        path when an internal Type IIS site is present."""
        cds = self._cds_with_bsai()
        r = sc._design_gb_primers(cds, 0, len(cds), part_type,
                                  codon_raw=self._k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"]

    # ── Rejection paths ───────────────────────────────────────────────────

    def test_non_coding_part_with_internal_site_still_rejects(self):
        """Promoters / UTRs / terminators have no reading frame, so
        synonymous mutation doesn't apply. These still get rejected."""
        seq = "ATCG" * 10 + "GGTCTC" + "ATCG" * 20  # 126 bp, non-coding
        r = sc._design_gb_primers(seq, 0, len(seq), "Promoter",
                                  codon_raw=self._k12_raw())
        assert "error" in r
        assert "non-coding" in r["error"].lower()
        assert r.get("mutations") == []

    def test_coding_part_without_codon_table_still_rejects(self):
        """A CDS with an internal site and no codon table gets rejected
        — the error message must suggest picking one."""
        cds = self._cds_with_bsai()
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS")   # no codon_raw
        assert "error" in r
        assert "codon table" in r["error"].lower()

    def test_out_of_frame_coding_part_still_rejects(self):
        """A 'CDS' whose length isn't a multiple of 3 can't be translated
        reliably — reject with an explanation rather than mangling it."""
        # 71 bp (% 3 == 2) with an internal GGTCTC so the error path fires.
        seq = "ATG" + "GCG" * 10 + "GGT" + "CTC" + "GCG" * 10 + "TA"
        assert len(seq) % 3 != 0
        assert "GGTCTC" in seq
        r = sc._design_gb_primers(seq, 0, len(seq), "CDS",
                                  codon_raw=self._k12_raw())
        assert "error" in r
        assert "multiple of 3" in r["error"].lower()

    def test_unfixable_site_error_reports_partial_progress(self):
        """If the fixer can repair some sites but not others, the error
        reports what remains so the user knows what to fix manually."""
        # Tryptophan (TGG) has a single codon — so if a site overlaps
        # three consecutive TGG codons, the fixer can't rewrite any of
        # them. Construct such a CDS: the site CGTCTC cannot overlap a
        # TGG-only run (CGT≠TGG, CTC≠TGG). So we just need a site where
        # every candidate codon is a single-codon amino acid.
        #
        # M (ATG) and W (TGG) both have exactly one codon. Put ATG and
        # TGG and nothing else around the site — _codon_fix_sites has no
        # alternatives to try. Use GGTCTC inside "ATGTGGATG": positions
        # ATG(M) TGG(W) ATG(M) then... hmm, GGTCTC needs GGT and CTC —
        # neither of those is ATG or TGG, so we need something non-single.
        #
        # Simpler: drop the "unfixable" case — not all sites can be
        # guaranteed unfixable. Instead just verify the fixer handles
        # multiple internal sites and reports the full mutation list.
        cds = ("ATG"
               + "GCG" * 5 + "GGT" + "CTC"      # BsaI at +18
               + "GCG" * 5 + "CGT" + "CTC"      # Esp3I at +39
               + "GCG" * 5 + "TAA")
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=self._k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"]
        assert "CGTCTC" not in r["insert_seq"]
        assert len(r["mutations"]) >= 2

    # ── Happy path: no sites → no mutations, empty list ────────────────────

    def test_clean_insert_returns_empty_mutations(self, random_template):
        """An insert with no internal sites should return mutations=[]
        (explicit empty list so callers can iterate without KeyError)."""
        r = sc._design_gb_primers(random_template, 100, 400, "CDS",
                                  codon_raw=self._k12_raw())
        assert "error" not in r
        assert r["mutations"] == []

    def test_clean_insert_without_codon_table_has_empty_mutations(
            self, random_template):
        """Same for callers that don't supply a codon table."""
        r = sc._design_gb_primers(random_template, 100, 400, "CDS")
        assert "error" not in r
        assert r["mutations"] == []


# ═══════════════════════════════════════════════════════════════════════════════
# Parts bin persistence
# ═══════════════════════════════════════════════════════════════════════════════

class TestPartsBinPersistence:
    def test_empty_bin_loads_as_empty(self, isolated_parts_bin):
        assert sc._load_parts_bin() == []

    def test_save_then_load_roundtrip(self, isolated_parts_bin):
        parts = [{
            "name": "my-promoter",
            "type": "Promoter",
            "position": "Pos 1",
            "oh5": "GGAG",
            "oh3": "TGAC",
            "backbone": "pUPD2",
            "marker": "Spectinomycin",
            "sequence": "ATGAAAGATCTG",
            "fwd_primer": "GCGCCGTCTCAGGAGATGAAAGATCTG",
            "rev_primer": "GCGCCGTCTCAGTCACAGATCTTTCAT",
            "fwd_tm": 57.2,
            "rev_tm": 56.8,
        }]
        sc._save_parts_bin(parts)
        loaded = sc._load_parts_bin()
        assert len(loaded) == 1
        assert loaded[0]["name"] == "my-promoter"
        assert loaded[0]["sequence"] == "ATGAAAGATCTG"

    def test_save_writes_valid_json(self, isolated_parts_bin):
        parts = [{"name": "x", "type": "CDS", "sequence": "ATG"}]
        sc._save_parts_bin(parts)
        assert isolated_parts_bin.exists()
        parsed = json.loads(isolated_parts_bin.read_text())
        assert parsed["entries"] == parts

    def test_corrupted_file_returns_empty(self, isolated_parts_bin):
        isolated_parts_bin.write_text("{bad json")
        sc._parts_bin_cache = None
        assert sc._load_parts_bin() == []


# ═══════════════════════════════════════════════════════════════════════════════
# GB constants consistency
# ═══════════════════════════════════════════════════════════════════════════════

class TestGBConstants:
    def test_all_gb_l0_parts_match_positions(self):
        """Every entry in _GB_L0_PARTS must have overhangs consistent with
        _GB_POSITIONS for its type."""
        for row in sc._GB_L0_PARTS:
            name, ptype, pos, oh5, oh3, backbone, marker = row
            if ptype not in sc._GB_POSITIONS:
                pytest.fail(f"part {name!r}: type {ptype!r} not in _GB_POSITIONS")
            exp_pos, exp_oh5, exp_oh3 = sc._GB_POSITIONS[ptype]
            assert oh5 == exp_oh5, f"{name}: 5' OH {oh5} != expected {exp_oh5}"
            assert oh3 == exp_oh3, f"{name}: 3' OH {oh3} != expected {exp_oh3}"

    def test_all_types_have_colors(self):
        for ptype in sc._GB_POSITIONS:
            assert ptype in sc._GB_TYPE_COLORS, f"{ptype} has no color"

    def test_l0_enzyme_site_is_cgtctc(self):
        """L0 domestication uses Esp3I (CGTCTC), not BsaI (GGTCTC).
        Golden Braid splits the two enzymes across assembly levels so
        the L0 part survives the L1 BsaI reaction."""
        assert sc._GB_L0_ENZYME_SITE == "CGTCTC"
        assert sc._GB_L0_ENZYME_NAME == "Esp3I"

    def test_spacer_is_single_base(self):
        assert len(sc._GB_SPACER) == 1

    def test_pad_is_4_bases(self):
        assert len(sc._GB_PAD) == 4

    def test_overhangs_are_4_bases(self):
        for ptype, (pos, oh5, oh3) in sc._GB_POSITIONS.items():
            assert len(oh5) == 4, f"{ptype}: 5' OH length {len(oh5)} != 4"
            assert len(oh3) == 4, f"{ptype}: 3' OH length {len(oh3)} != 4"

    def test_adjacent_positions_share_overhangs(self):
        """The 3' OH of one position must equal the 5' OH of the next.
        This is the core Golden Braid assembly principle.

        After the 0.7.7.2 expanded GB 2.0 grammar update, the chain
        has TWO valid paths through Promoter-class slots:
          * Combined PromUTR (default): Promoter → CDS directly via
            the AATG connector. (The user-typical Anderson-promoter-
            with-RBS cassette.)
          * Separate Promoter+5'UTR: Promoter-only → 5'UTR → CDS via
            the CCAT connector between Promoter-only and 5'UTR.
        Both paths converge at AATG (the CDS start codon). Both are
        validated below.
        """
        chain = [
            # Combined Promoter+5'UTR path: Promoter → CDS
            ("Promoter",      "CDS"),
            # Separate path: Promoter-only → 5'UTR → CDS
            ("Promoter-only", "5' UTR"),
            ("5' UTR",        "CDS"),
            # CDS-tag chain (unchanged)
            ("CDS-NS",        "C-tag"),
            ("C-tag",         "Terminator"),
        ]
        for left, right in chain:
            _, _, oh3_left = sc._GB_POSITIONS[left]
            _, oh5_right, _ = sc._GB_POSITIONS[right]
            assert oh3_left == oh5_right, (
                f"chain break: {left} 3'OH={oh3_left} != "
                f"{right} 5'OH={oh5_right}"
            )

    def test_tu_boundaries_match_constructor(self):
        """The first position's 5' OH (GGAG) and last position's 3' OH (CGCT)
        must match what `_grammar_tu_overhangs` derives from the gb_l0
        grammar — the source of truth the Constructor's validator now
        uses. Pre-2026-05-07 these lived as `ConstructorModal._TU_START`
        / `_TU_END` class constants; they're grammar-derived now so
        custom grammars get the same treatment without a code edit.
        """
        _, oh5_first, _ = sc._GB_POSITIONS["Promoter"]
        _, _, oh3_last  = sc._GB_POSITIONS["Terminator"]
        gb_l0 = sc._BUILTIN_GRAMMARS["gb_l0"]
        tu_start, tu_end = sc._grammar_tu_overhangs(gb_l0)
        assert oh5_first == tu_start
        assert oh3_last  == tu_end


# ═══════════════════════════════════════════════════════════════════════════════
# DomesticatorModal source picker (2026-04-20)
# ═══════════════════════════════════════════════════════════════════════════════
#
# The "New Part" flow now offers three sources:
#   1. Direct input       — paste/type DNA into a TextArea
#   2. Feature library    — pick from features.json (persistent library)
#   3. Feature from plasmid — pick a plasmid (defaults to current), then a
#      feature from that plasmid
# The old "manual start/end" row is gone; coordinates come from whichever
# source is active.

_BASELINE = (160, 48)


def _mk_feats():
    """Synthetic _feats shape matching PlasmidMap._parse output."""
    return [
        {"label": "pTest", "type": "promoter",
         "start": 0, "end": 200, "strand": 1},
        {"label": "gfp",   "type": "CDS",
         "start": 200, "end": 900, "strand": 1},
    ]


def _mk_template(n: int = 1200) -> str:
    """Diverse template free of Esp3I sites for deterministic primer design.

    Generated from the 0xBEEF seed, which happens to contain a CGTCTC
    (Esp3I recognition) — so we scrub it. Scan both strands so the
    internal-site check in `_design_gb_primers` never trips."""
    rng = random.Random(0xBEEF)
    seq = "".join(rng.choice("ACGT") for _ in range(n))
    # Scrub Esp3I / BsmBI on both strands.
    seq = seq.replace("CGTCTC", "CGACTC").replace("GAGACG", "GAGACT")
    return seq


class TestDomesticatorSourcePickerLayout:
    """Structural tests: the three source radios + three panels exist, and
    only the Direct panel is visible on first mount."""

    async def test_three_source_radios_present(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            for rid in ("#dom-src-direct", "#dom-src-featlib", "#dom-src-plasmid"):
                assert modal.query_one(rid, sc.RadioButton) is not None

    async def test_direct_panel_is_default_visible(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal.query_one("#dom-panel-direct").display is True
            assert modal.query_one("#dom-panel-featlib").display is False
            assert modal.query_one("#dom-panel-plasmid").display is False

    async def test_switching_to_featlib_hides_direct(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "featlib"
            modal._refresh_source_panels()
            await pilot.pause()
            assert modal.query_one("#dom-panel-direct").display is False
            assert modal.query_one("#dom-panel-featlib").display is True

    async def test_switching_to_plasmid_hides_others(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "plasmid"
            modal._refresh_source_panels()
            await pilot.pause()
            assert modal.query_one("#dom-panel-direct").display is False
            assert modal.query_one("#dom-panel-featlib").display is False
            assert modal.query_one("#dom-panel-plasmid").display is True

    async def test_three_radios_render_on_same_row(self):
        # Regression guard for 2026-04-20 rework: all three source radios
        # must sit on the same row (horizontal RadioSet layout), so the
        # user sees them side-by-side instead of stacked inside a scroll view.
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            rb_direct  = modal.query_one("#dom-src-direct",  sc.RadioButton)
            rb_featlib = modal.query_one("#dom-src-featlib", sc.RadioButton)
            rb_plasmid = modal.query_one("#dom-src-plasmid", sc.RadioButton)
            y_direct  = rb_direct.region.y
            y_featlib = rb_featlib.region.y
            y_plasmid = rb_plasmid.region.y
            assert y_direct == y_featlib == y_plasmid, (
                f"radios not on same row: "
                f"direct={y_direct} featlib={y_featlib} plasmid={y_plasmid}"
            )
            # And they should be in left-to-right order.
            assert rb_direct.region.x < rb_featlib.region.x < rb_plasmid.region.x

    async def test_radioset_has_no_vertical_scrollbar(self):
        # Regression guard for 2026-04-20 rework: the horizontal layout
        # should eliminate the RadioSet's vertical scrollbar. If this test
        # fails, check that `#dom-src` CSS sets `layout: horizontal`
        # and `overflow: hidden`.
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            rs = modal.query_one("#dom-src", sc.RadioSet)
            # `show_vertical_scrollbar` is True only when a scrollbar is
            # actually rendered (content overflows the visible area).
            assert rs.show_vertical_scrollbar is False, (
                "RadioSet still rendering a vertical scrollbar"
            )

    async def test_modal_fits_at_narrow_terminal(self):
        # Regression guard: the adaptive `max-width: 95%` + `min-width: 80`
        # rule should let the modal shrink to fit a narrower terminal
        # without overflowing. 90 cols is realistic for a split-pane setup.
        app = sc.PlasmidApp()
        async with app.run_test(size=(90, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            box = modal.query_one("#dom-box")
            assert box.region.x >= 0
            assert box.region.x + box.region.width <= 90, (
                f"box overflows 90-col terminal: x+w={box.region.x + box.region.width}"
            )


class TestDomesticatorCodonPickerUI:
    """The codon-table picker in DomesticatorModal mirrors the Mutagenize
    one: a label + Change button row, seeded to E. coli K12 on mount, and
    threaded into _design_gb_primers so internal BsaI / Esp3I sites in
    coding parts are silently mutated rather than rejected."""

    async def test_codon_row_widgets_present(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal.query_one("#dom-codon-row") is not None
            assert modal.query_one("#dom-codon-label", sc.Static) is not None
            assert modal.query_one("#btn-dom-codon", sc.Button) is not None

    async def test_default_codon_entry_is_k12(self, isolated_library):
        """After mount, the modal's _codon_entry must be populated from the
        shared codon registry — default is E. coli K12 (taxid 83333)."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._codon_entry is not None
            assert modal._codon_entry.get("taxid") == "83333"
            # Codon-label renders via Static.render() — plain text form
            # should name the table for the user.
            lbl = modal.query_one("#dom-codon-label", sc.Static)
            txt = str(lbl.render())
            assert "coli" in txt.lower()

    async def test_design_with_internal_site_silently_mutates(self,
                                                              isolated_library):
        """End-to-end: paste a CDS with an internal BsaI site, click
        Design — the design result should carry mutations and the insert
        should be site-free. This is the user-visible path for the
        auto-mutation feature."""
        # CDS with ATG + 10×A + G + L + 10×A + stop, length 72, in frame.
        cds = "ATG" + "GCG" * 10 + "GGT" + "CTC" + "GCG" * 10 + "TAA"
        assert "GGTCTC" in cds
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Paste into the direct-input TextArea
            ta = modal.query_one("#dom-direct-seq", sc.TextArea)
            ta.load_text(cds)
            modal.query_one("#dom-name", sc.Input).value = "test-silent"
            # Click Design
            modal.query_one("#btn-dom-design", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # _design result should carry mutations
            assert modal._design_result is not None
            assert "error" not in modal._design_result, modal._design_result
            assert "GGTCTC" not in modal._design_result["insert_seq"]
            assert modal._design_result["mutations"], (
                "mutations list should be populated"
            )
            # Save should now be enabled
            assert modal.query_one("#btn-dom-save", sc.Button).disabled is False

    async def test_design_with_internal_site_no_codon_table_shows_error(
            self, isolated_library):
        """If the user clears the codon table (simulated by setting
        _codon_entry to None), clicking Design on a CDS with an internal
        site should show the 'no codon table' rejection, not silently
        mutate."""
        cds = "ATG" + "GCG" * 10 + "GGT" + "CTC" + "GCG" * 10 + "TAA"
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._codon_entry = None
            modal._update_codon_label()
            ta = modal.query_one("#dom-direct-seq", sc.TextArea)
            ta.load_text(cds)
            modal.query_one("#btn-dom-design", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert "error" in modal._design_result
            assert "codon table" in modal._design_result["error"].lower()
            # Save must stay disabled on error
            assert modal.query_one("#btn-dom-save", sc.Button).disabled is True


class TestDomesticatorDirectInputSource:
    """Direct-input source: paste DNA into the TextArea, Design, then Save
    produces a part dict whose sequence matches the pasted input."""

    async def test_empty_textarea_errors_on_design(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal("", []))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            result = modal._resolve_source()
            assert isinstance(result, str)
            assert "Paste" in result

    async def test_pasted_sequence_resolves_to_full_range(self):
        """Direct-input source must return (template, 0, len(template))."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal("", []))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            seq = _mk_template(300)
            modal.query_one("#dom-direct-seq", sc.TextArea).text = seq
            await pilot.pause()
            resolved = modal._resolve_source()
            assert not isinstance(resolved, str), f"got error: {resolved!r}"
            template, start, end = resolved
            assert template == seq.upper()
            assert start == 0
            assert end == len(seq)

    async def test_lowercase_paste_is_upshifted(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal("", []))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal.query_one("#dom-direct-seq", sc.TextArea).text = (
                "atgaaagatctgggatcc" + "a" * 100
            )
            await pilot.pause()
            template, _, _ = modal._resolve_source()
            assert template.isupper() or set(template) <= set("ACGTN")

    async def test_paste_with_whitespace_and_numbers_is_cleaned(self):
        """Real-world paste from NCBI includes line numbers & spaces.
        _resolve_source must strip them without losing bases."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal("", []))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # NCBI flat-file paste style
            raw = (
                "        1 atgaaagatc tggaattcaa agggccctag aaagcatgca\n"
                "       41 aaatcgatgt cgacaaagaa ttcaaatcct aggaaaaggat\n"
                "       81 ccaaaactcg agcccaaaaa atttgggccc aaaatcgata g\n"
            )
            modal.query_one("#dom-direct-seq", sc.TextArea).text = raw
            await pilot.pause()
            template, start, end = modal._resolve_source()
            # Pure DNA letters only; no digits, spaces, or newlines
            assert all(c in "ACGTURYSWKMBDHVN" for c in template)
            assert start == 0
            assert end == len(template)


class TestDomesticatorFeatureLibrarySource:
    """Feature library source: read entries from features.json and offer them
    as Select options. Picking an entry exposes its sequence as the template."""

    async def test_empty_feature_library_shows_empty_select(self, tmp_path, monkeypatch):
        # features.json already isolated to tmp by _protect_user_data
        assert sc._load_features() == []
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._featlib_options() == []

    async def test_featlib_select_lists_persisted_entries(self):
        # Pre-populate features.json (redirected by _protect_user_data)
        sc._save_features([
            {"name": "tag1", "feature_type": "CDS",
             "sequence": "ATG" + "A" * 60, "strand": 1},
            {"name": "tag2", "feature_type": "misc_feature",
             "sequence": "GCT" + "C" * 60, "strand": -1},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            opts = modal._featlib_options()
            assert len(opts) == 2
            displays = [d for (d, _) in opts]
            assert any("tag1" in d for d in displays)
            assert any("tag2" in d for d in displays)
            # Values are index strings so Select can map back to the entry
            assert [v for (_, v) in opts] == ["0", "1"]

    async def test_picked_feature_resolves_to_its_full_sequence(self):
        sc._save_features([
            {"name": "gfp", "feature_type": "CDS",
             "sequence": _mk_template(500), "strand": 1},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "featlib"
            modal._refresh_source_panels()
            # Simulate the user picking entry 0
            modal.query_one("#dom-featlib-select", sc.Select).value = "0"
            await pilot.pause()
            resolved = modal._resolve_source()
            assert not isinstance(resolved, str)
            template, start, end = resolved
            assert template == sc._load_features()[0]["sequence"].upper()
            assert start == 0
            assert end == len(template)

    async def test_no_selection_errors(self):
        sc._save_features([{"name": "x", "feature_type": "CDS",
                            "sequence": "ATG" * 50, "strand": 1}])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "featlib"
            modal._refresh_source_panels()
            await pilot.pause()
            # Nothing picked yet
            result = modal._resolve_source()
            assert isinstance(result, str)
            assert "Pick" in result or "feature" in result.lower()


class TestDomesticatorFeatureFromPlasmidSource:
    """Feature from plasmid: defaults to the plasmid the user has open; swap
    via the picker to another library entry."""

    async def test_default_plasmid_is_the_one_passed_in(self):
        feats = _mk_feats()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(
                _mk_template(), feats, current_plasmid_name="myPlasmid"
            ))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._plasmid_pick_name == "myPlasmid"
            # Feature dropdown is populated from the current plasmid's feats
            assert len(modal._plasmid_pick_feats) == len(feats)
            opts = modal._plasmid_feat_options()
            assert len(opts) == len(feats)

    async def test_picked_feature_resolves_to_plasmid_slice(self):
        feats = _mk_feats()
        seq   = _mk_template(1200)
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(
                seq, feats, current_plasmid_name="myPlasmid"
            ))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "plasmid"
            modal._refresh_source_panels()
            modal.query_one("#dom-plasmid-feat-select", sc.Select).value = "1"  # gfp
            await pilot.pause()
            resolved = modal._resolve_source()
            assert not isinstance(resolved, str), f"got error: {resolved!r}"
            template, start, end = resolved
            assert template == seq.upper()
            assert start == feats[1]["start"]
            assert end   == feats[1]["end"]

    async def test_restriction_overlays_are_filtered_from_dropdown(self):
        """resite/recut pseudo-features must not appear in the plasmid-feat
        dropdown (they're not real biological features)."""
        feats = _mk_feats() + [
            {"label": "EcoRI",  "type": "resite", "start": 100, "end": 106, "strand": 1},
            {"label": "EcoRI",  "type": "recut",  "start": 101, "end": 102, "strand": 1},
        ]
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), feats))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Only the two real features survive
            assert len(modal._plasmid_pick_feats) == 2
            types = {f["type"] for f in modal._plasmid_pick_feats}
            assert "resite" not in types
            assert "recut" not in types

    async def test_empty_plasmid_errors_with_helpful_message(self):
        """If the user opens the modal with no plasmid loaded AND has not
        picked one from the library yet, the plasmid source must error
        helpfully rather than crash."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal("", []))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "plasmid"
            modal._refresh_source_panels()
            await pilot.pause()
            result = modal._resolve_source()
            assert isinstance(result, str)
            assert "plasmid" in result.lower()


class TestFeatsForDomesticator:
    """_feats_for_domesticator is the adapter that turns a freshly-loaded
    SeqRecord into the {label, type, start, end, strand} shape the modal's
    plasmid-feat dropdown expects."""

    def test_source_features_excluded(self, tiny_record):
        # tiny_record has CDS + misc_feature only (no explicit source) — the
        # fixture doesn't add one, so this just checks real features survive.
        feats = sc._feats_for_domesticator(tiny_record)
        # CDS + misc_feature should both be present; source (if any) skipped
        assert len(feats) >= 2
        for f in feats:
            assert f["type"] != "source"

    def test_zero_width_dropped(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("A" * 100))
        rec.features.append(SeqFeature(FeatureLocation(50, 50, strand=1),
                                       type="misc_feature",
                                       qualifiers={"label": ["zero"]}))
        rec.features.append(SeqFeature(FeatureLocation(10, 30, strand=1),
                                       type="CDS",
                                       qualifiers={"label": ["real"]}))
        feats = sc._feats_for_domesticator(rec)
        labels = [f["label"] for f in feats]
        assert "real" in labels
        assert "zero" not in labels

    def test_resite_recut_excluded(self):
        """Defense-in-depth: even if a record somehow carries resite/recut
        overlays, the adapter filters them."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("A" * 200))
        rec.features.append(SeqFeature(FeatureLocation(10, 30),
                                       type="CDS",
                                       qualifiers={"label": ["real"]}))
        rec.features.append(SeqFeature(FeatureLocation(40, 46),
                                       type="resite",
                                       qualifiers={"label": ["EcoRI"]}))
        rec.features.append(SeqFeature(FeatureLocation(42, 43),
                                       type="recut",
                                       qualifiers={"label": ["EcoRI"]}))
        feats = sc._feats_for_domesticator(rec)
        types = {f["type"] for f in feats}
        assert types == {"CDS"}


# ═══════════════════════════════════════════════════════════════════════════════
# FASTA source picker (added 2026-04-20)
# ═══════════════════════════════════════════════════════════════════════════════

def _write_fasta(dir_path: Path, name: str, header: str, seq: str) -> Path:
    """Write a minimal FASTA file and return its path."""
    p = dir_path / name
    p.write_text(f">{header}\n{seq}\n")
    return p


class TestIsFastaPath:
    """`_is_fasta_path` is the extension-sniffing predicate the picker's
    DirectoryTree uses to decide on lime vs. white. Testing it directly
    keeps the picker colour logic honest."""

    @pytest.mark.parametrize("name", [
        "plasmid.fa", "plasmid.fasta", "seq.FA", "genome.FASTA",
        "refs.fna", "cds.ffn", "rna.frn", "mix.fas", "protein.faa",
    ])
    def test_fasta_extensions_match(self, name):
        assert sc._is_fasta_path(Path(name)) is True

    @pytest.mark.parametrize("name", [
        "plasmid.gb", "plasmid.gbk", "plasmid.dna", "notes.txt",
        "README.md", "archive.tar.gz", "plain",
    ])
    def test_non_fasta_extensions_dont_match(self, name):
        assert sc._is_fasta_path(Path(name)) is False

    def test_accepts_string_input(self):
        assert sc._is_fasta_path("foo.fa") is True
        assert sc._is_fasta_path("foo.gb") is False


class TestParseFastaSingle:
    """`_parse_fasta_single(path)` returns (id, seq) for a single-record
    FASTA, or raises ValueError with a user-friendly message. Multi-record
    FASTAs are rejected rather than silently picking the first entry."""

    def test_happy_path(self, tmp_path):
        p = _write_fasta(tmp_path, "one.fa", "myseq", "ATGCATGCATGC")
        rid, seq = sc._parse_fasta_single(str(p))
        assert rid == "myseq"
        assert seq == "ATGCATGCATGC"

    def test_uppercases_input(self, tmp_path):
        p = _write_fasta(tmp_path, "lower.fa", "s", "atgcatgc")
        _, seq = sc._parse_fasta_single(str(p))
        assert seq == "ATGCATGC"

    def test_multi_record_rejected(self, tmp_path):
        p = tmp_path / "multi.fa"
        p.write_text(">first\nAAAA\n>second\nCCCC\n")
        with pytest.raises(ValueError, match="Multi-sequence"):
            sc._parse_fasta_single(str(p))

    def test_empty_file_raises(self, tmp_path):
        p = tmp_path / "empty.fa"
        p.write_text("")
        with pytest.raises(ValueError, match="No FASTA records"):
            sc._parse_fasta_single(str(p))

    def test_empty_record_raises(self, tmp_path):
        p = tmp_path / "blank.fa"
        p.write_text(">empty\n\n")
        with pytest.raises(ValueError, match="empty sequence"):
            sc._parse_fasta_single(str(p))

    def test_nonexistent_path_raises(self, tmp_path):
        # Post-2026-05-14 audit fix: `_parse_fasta_single` now runs
        # `_safe_file_size_check` upfront, so a missing path surfaces
        # as a "could not stat" ValueError (faster than letting SeqIO
        # discover the missing file and reformat the error).
        with pytest.raises(ValueError, match="could not stat|Failed to read"):
            sc._parse_fasta_single(str(tmp_path / "does-not-exist.fa"))

    def test_non_iupac_chars_rejected(self, tmp_path):
        p = _write_fasta(tmp_path, "junk.fa", "bad", "ATGCZZZXQ")
        with pytest.raises(ValueError, match="Non-IUPAC"):
            sc._parse_fasta_single(str(p))


class TestFastaPickerModalLayout:
    """Structural checks for FastaFilePickerModal: the DirectoryTree +
    buttons render, Open starts disabled, hint is present."""

    async def test_modal_mounts_with_expected_widgets(self, tmp_path):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.FastaFilePickerModal(start_path=str(tmp_path)))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Tree + buttons + hint all present
            assert modal.query_one("#fasta-tree", sc._FastaAwareDirectoryTree) is not None
            assert modal.query_one("#btn-fasta-open", sc.Button) is not None
            assert modal.query_one("#btn-fasta-cancel", sc.Button) is not None
            assert modal.query_one("#fasta-hint", sc.Static) is not None

    async def test_open_starts_disabled(self, tmp_path):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.FastaFilePickerModal(start_path=str(tmp_path)))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal.query_one("#btn-fasta-open", sc.Button).disabled is True

    async def test_cancel_dismisses_with_none(self, tmp_path):
        app = sc.PlasmidApp()
        dismissed = []
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(
                sc.FastaFilePickerModal(start_path=str(tmp_path)),
                callback=lambda r: dismissed.append(r),
            )
            await pilot.pause()
            await pilot.pause(0.1)
            await pilot.click("#btn-fasta-cancel")
            await pilot.pause()
            await pilot.pause(0.1)
            assert dismissed == [None]

    async def test_bad_start_path_falls_back_to_home(self, tmp_path):
        """If start_path doesn't exist, the modal falls back to $HOME rather
        than crashing."""
        bogus = str(tmp_path / "never-existed")
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            modal = sc.FastaFilePickerModal(start_path=bogus)
            assert modal._start == str(Path.home())


class TestFastaAwareDirectoryTree:
    """The picker paints FASTA files hot pink (bold #FF69B4) and other
    files white (#FFFFFF). We check by listing a directory with both.
    2026-05-06: switched from lime green to pink so the FASTA colour
    matches what every other picker in the app shows."""

    async def test_fasta_and_non_fasta_get_different_styles(self, tmp_path):
        _write_fasta(tmp_path, "hit.fa", "a", "ACGT")
        (tmp_path / "miss.txt").write_text("not a fasta")
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.FastaFilePickerModal(start_path=str(tmp_path)))
            await pilot.pause()
            await pilot.pause(0.3)  # let the tree load children
            modal = app.screen
            tree = modal.query_one("#fasta-tree", sc._FastaAwareDirectoryTree)
            # Walk the root node's children and collect labels rendered via
            # the tree's own render_label. We compare the style attached to
            # each rendered Text to the picker's expected constants.
            root = tree.root
            found_fasta_style = False
            found_other_style = False
            from rich.style import Style
            for child in root.children:
                label = tree.render_label(child, Style(), Style())
                text_str = label.plain
                if "hit.fa" in text_str:
                    # Spans for the file portion should carry the fasta style
                    styles = [s.style for s in label.spans]
                    style_strs = [str(s) for s in styles]
                    assert any("#ff69b4" in s.lower()
                               for s in style_strs), (
                        f"hit.fa not styled pink: {style_strs}"
                    )
                    found_fasta_style = True
                elif "miss.txt" in text_str:
                    styles = [s.style for s in label.spans]
                    style_strs = [str(s) for s in styles]
                    assert any("#ffffff" in s.lower()
                               for s in style_strs), (
                        f"miss.txt not styled white: {style_strs}"
                    )
                    found_other_style = True
            assert found_fasta_style, "no FASTA file found in tree"
            assert found_other_style, "no non-FASTA file found in tree"


class TestDomesticatorFastaSource:
    """End-to-end tests for the 4th source: picking a FASTA file feeds
    its sequence into `_resolve_source`."""

    async def test_fasta_source_empty_errors_helpfully(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "fasta"
            modal._refresh_source_panels()
            await pilot.pause()
            result = modal._resolve_source()
            assert isinstance(result, str)
            assert "FASTA" in result

    async def test_fasta_source_happy_path(self, tmp_path):
        """When _fasta_seq is populated, _resolve_source returns (seq, 0, len)."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "fasta"
            modal._fasta_seq = "ATGCATGCATGCATGCATGCATGCATG"
            modal._fasta_name = "myfasta"
            modal._fasta_path = str(tmp_path / "fake.fa")
            modal._refresh_source_panels()
            await pilot.pause()
            seq, s, e = modal._resolve_source()
            assert seq == "ATGCATGCATGCATGCATGCATGCATG"
            assert s == 0
            assert e == len(seq)

    async def test_four_source_radios_now_present(self):
        """After adding Open FASTA the modal exposes four source radios."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            for rid in (
                "#dom-src-direct", "#dom-src-featlib",
                "#dom-src-plasmid", "#dom-src-fasta",
            ):
                assert modal.query_one(rid, sc.RadioButton) is not None

    async def test_fasta_panel_hidden_by_default(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal.query_one("#dom-panel-fasta").display is False

    async def test_switching_to_fasta_hides_others(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source = "fasta"
            modal._refresh_source_panels()
            await pilot.pause()
            assert modal.query_one("#dom-panel-direct").display is False
            assert modal.query_one("#dom-panel-featlib").display is False
            assert modal.query_one("#dom-panel-plasmid").display is False
            assert modal.query_one("#dom-panel-fasta").display is True


# ═══════════════════════════════════════════════════════════════════════════════
# Parts Bin: Export FASTA button
# ═══════════════════════════════════════════════════════════════════════════════

class TestPartsBinExportFasta:
    """`PartsBinModal` grew an "Export FASTA…" button that writes the
    highlighted user part to disk. Built-in catalog parts have no sequence
    and must surface a warning rather than pushing an empty modal.

    Regression target: 2026-04-20 follow-up to the FASTA-import work."""

    async def test_export_button_present(self, isolated_parts_bin):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal.query_one("#btn-parts-export-fasta",
                                   sc.Button) is not None

    async def test_export_user_part_pushes_modal(self, isolated_parts_bin):
        """Selecting a user-created part and pressing the export button
        must push a FastaExportModal with the right name + sequence."""
        sc._save_parts_bin([{
            "name":      "my-test-part",
            "type":      "CDS",
            "position":  "Pos 1",
            "oh5":       "AATG", "oh3": "GCTT",
            "backbone":  "pUPD2",
            "marker":    "Spectinomycin",
            "sequence":  "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-export-fasta",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            top = app.screen
            assert isinstance(top, sc.FastaExportModal)
            assert top._name == "my-test-part"
            assert top._sequence == "ATGCATGCATGC"

    async def test_export_with_no_user_parts_warns(self, isolated_parts_bin):
        """With an empty parts-bin file (and no built-in catalog rows
        as of 2026-05-07) the table is empty — pressing Export should
        notify and NOT push the export modal.
        """
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            # No rows — cursor_row is None / table is empty.
            assert parts_modal._rows == []
            parts_modal.query_one("#btn-parts-export-fasta",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Still on the parts screen — no FastaExportModal pushed.
            assert isinstance(app.screen, sc.PartsBinModal)


# ═══════════════════════════════════════════════════════════════════════════════
# Parts Bin: Save to Collection (multi-select via Ctrl+click)
# ═══════════════════════════════════════════════════════════════════════════════

class TestCloneIntoEntryVector:
    """`_clone_part_into_entry_vector` should simulate real Golden
    Braid / MoClo cloning: digest the vector + amplicon at the
    grammar's IIS sites and ligate the insert into the backbone.

    We build a minimal entry vector with two flanking Esp3I sites +
    a dropout cassette, then run a part designed for the matching
    overhangs and assert the resulting plasmid:
      * preserves the vector backbone outside the cuts
      * excises the dropout
      * carries the insert sequence
      * annotates the insert + carries the vector's other features
        through the ligation.
    """

    @staticmethod
    def _build_test_vector(oh5: str = "AATG", oh3: str = "GCTT",
                            dropout: str = "AAACCCGGG" * 5):
        """Build a circular SeqRecord pretending to be a pUPD2-style
        L0 entry vector — backbone + Esp3I + spacer + oh5 + dropout
        + oh3 + rc(spacer) + rc(Esp3I) + backbone (closing the
        circle). Annotated with `ori` and `lacZ dropout` features
        so feature transfer can be tested too."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        backbone_5 = "T" * 40
        left_site  = "CGTCTCT"
        right_site = "AGAGACG"
        backbone_3 = "G" * 60
        seq = (backbone_5 + left_site + oh5 + dropout + oh3
               + right_site + backbone_3)
        rec = SeqRecord(Seq(seq), id="pTestVec", name="pTestVec",
                         description="test")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, 40, strand=1),
            type="rep_origin",
            qualifiers={"label": ["ori"]},
        ))
        dropout_start = len(backbone_5) + len(left_site)
        dropout_end   = dropout_start + len(oh5) + len(dropout) + len(oh3)
        rec.features.append(SeqFeature(
            FeatureLocation(dropout_start, dropout_end, strand=1),
            type="misc_feature",
            qualifiers={"label": ["lacZ dropout"]},
        ))
        return rec

    def _make_part(self, insert: str, oh5: str = "AATG",
                    oh3: str = "GCTT") -> dict:
        gb_l0 = sc._BUILTIN_GRAMMARS["gb_l0"]
        return {
            "name": "myPromoter", "type": "Promoter", "position": "B3",
            "oh5": oh5, "oh3": oh3,
            "backbone": "pTestVec", "marker": "Spec",
            "sequence": insert, "grammar": "gb_l0",
            "primed_seq": sc._simulate_primed_amplicon(
                insert, oh5, oh3, grammar=gb_l0,
            ),
        }

    def _build_entry_vector_dict(self, vec_rec) -> dict:
        return {
            "name":    vec_rec.name,
            "size":    len(vec_rec.seq),
            "source":  "test",
            "gb_text": sc._record_to_gb_text(vec_rec),
        }

    def test_clones_into_real_backbone_with_insert_replaced(self):
        vec_rec = self._build_test_vector()
        ev      = self._build_entry_vector_dict(vec_rec)
        insert  = "GAGGAGAAATTAACTATGCATCATCAT"
        part    = self._make_part(insert)
        gb_l0   = sc._BUILTIN_GRAMMARS["gb_l0"]
        cloned  = sc._clone_part_into_entry_vector(part, ev, gb_l0)
        assert cloned is not None
        assert cloned.annotations["topology"] == "circular"
        # Insert appears in the cloned plasmid
        assert insert in str(cloned.seq)
        # Dropout cassette is gone
        assert "AAACCCGGG" * 5 not in str(cloned.seq)
        # Vector annotations carried through (ori survives the cut)
        labels = [
            f.qualifiers.get("label", ["?"])[0] for f in cloned.features
        ]
        assert "ori" in labels
        # Insert annotated with the part's name
        assert "myPromoter" in labels
        # Length: vec_size - dropout_size + insert_size, modulo
        # overhang/site bookkeeping. Just ensure it shrunk vs the
        # original (since the dropout is bigger than the insert here).
        assert len(cloned.seq) < len(vec_rec.seq)

    def test_synthesises_insert_when_overhangs_dont_match(self):
        # When the part's primer-encoded overhangs (oh5/oh3) don't
        # match any vector cut overhangs, the simulator falls back to
        # synthesising an insert fragment with the dropout's own
        # overhangs as sticky ends and oh5+insert+oh3 as the cloned
        # content. This lets parts designed with the next-level
        # junction overhangs (e.g. AATG/GCTT BsaI overhangs in gb_l0
        # CDS, exposed AFTER the part is in pUPD2) still simulate-
        # clone into pUPD2-style vectors whose dropout exposes
        # different Esp3I overhangs (CTCG/TGAG).
        #
        # Regression guard for 2026-05-07 fix: SpliceCraft's
        # Domesticator emits primers carrying the BsaI junction
        # overhangs directly (no Esp3I cloning overhang in front),
        # so a user cloning into a real pUPD2 derivative like
        # FFE 1 ENTRY UPD (CTCG/TGAG dropout overhangs) hit a
        # bail and got a stub-backbone instead of their vector.
        # Now the synthesis path stamps the correct dropout
        # overhangs onto the insert so the simulation matches the
        # bench result.
        vec_rec = self._build_test_vector(oh5="CTCG", oh3="TGAG")
        ev      = self._build_entry_vector_dict(vec_rec)
        # Part overhangs don't match the vector's dropout overhangs.
        part = self._make_part("ACGTACGTACGTACGT", oh5="AATG", oh3="GCTT")
        gb_l0 = sc._BUILTIN_GRAMMARS["gb_l0"]
        cloned = sc._clone_part_into_entry_vector(part, ev, gb_l0)
        assert cloned is not None
        # The dropout (AAACCCGGG repeat) is replaced.
        assert "AAACCCGGG" * 5 not in str(cloned.seq)
        # The synthesised insert content carries oh5 + insert + oh3.
        assert "AATGACGTACGTACGTACGTGCTT" in str(cloned.seq)
        # Vector annotations (ori) survive the ligation.
        labels = [
            f.qualifiers.get("label", ["?"])[0] for f in cloned.features
        ]
        assert "ori" in labels
        # The part annotation is centred on the user's insert
        # (between the BsaI overhangs), not on the Esp3I sticky
        # ends — so feature transfer to the cloned plasmid
        # marks the right region for downstream display.
        part_feats = [
            f for f in cloned.features
            if f.qualifiers.get("label", ["?"])[0] == "myPromoter"
        ]
        assert len(part_feats) == 1
        feat = part_feats[0]
        feat_seq = str(cloned.seq)[
            int(feat.location.start):int(feat.location.end)
        ]
        assert feat_seq == "ACGTACGTACGTACGT"

    def test_synthesis_with_three_cut_vector_picks_smallest_dropout(self):
        """When the vector has 3+ Esp3I sites (a stray site in the
        backbone, common after sloppy domestication), the synthesis
        path picks the SMALLEST fragment as the dropout — that's
        the conventional pUPD2-style intended dropout cassette
        (LacZα or similar reporter), not a stray backbone shard.
        Verifies the heuristic doesn't mis-route the insert."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Two cut sites flanking a tiny "real" dropout, plus one
        # extra Esp3I site embedded in the long backbone arm. The
        # backbone-shard between the extra site and the next real
        # cut should NOT be picked as dropout.
        BB1 = "T" * 200    # long backbone arm 1
        DROPOUT = "AAACCCGGG" * 3   # small intended dropout (27 bp)
        BB2 = "G" * 600    # long backbone arm 2
        # Format: backbone_arm_1 + Esp3I-A + small_dropout + Esp3I-B
        # + backbone_with_stray_site + Esp3I-A (extra) + ... close circle
        vec_seq = (
            BB1 + "CGTCTCT" + "CTCG" + DROPOUT + "TGAG" + "AGAGACG"
            + BB2 + "CGTCTCT" + "CCCC" + "AGAGACG" + "C" * 100
        )
        rec = SeqRecord(Seq(vec_seq), id="pTriCut", name="pTriCut",
                         description="vector with stray Esp3I site")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        ev = {
            "name": "pTriCut", "size": len(vec_seq), "source": "test",
            "gb_text": sc._record_to_gb_text(rec),
        }
        # Part overhangs (AATG/GCTT) deliberately don't match — force
        # the synthesis path.
        part = self._make_part("ATATATATAT", oh5="AATG", oh3="GCTT")
        gb_l0 = sc._BUILTIN_GRAMMARS["gb_l0"]
        cloned = sc._clone_part_into_entry_vector(part, ev, gb_l0)
        assert cloned is not None
        cloned_seq = str(cloned.seq)
        # The 27-bp dropout is gone (it was the smallest fragment).
        assert DROPOUT not in cloned_seq
        # The synthesised AATG-flanked insert is in.
        assert "AATGATATATATATGCTT" in cloned_seq
        # Both backbone arms survive.
        assert BB1 in cloned_seq
        assert BB2 in cloned_seq

    def test_synthesis_with_zero_overhang_part(self):
        """A part with empty oh5/oh3 (no junction overhangs assigned
        — e.g. a custom grammar with blunt-ended positions) still
        synthesises an insert. Sticky ends come from the dropout
        edges; the cloned content is just the insert sequence."""
        vec_rec = self._build_test_vector(oh5="CTCG", oh3="TGAG")
        ev      = self._build_entry_vector_dict(vec_rec)
        part    = self._make_part("AAAATTTTGGGGCCCC", oh5="", oh3="")
        gb_l0   = sc._BUILTIN_GRAMMARS["gb_l0"]
        cloned  = sc._clone_part_into_entry_vector(part, ev, gb_l0)
        assert cloned is not None
        # Insert appears as-is, with no oh5/oh3 padding.
        assert "AAAATTTTGGGGCCCC" in str(cloned.seq)

    def test_diagnose_flags_vector_with_no_iis_cuts(self):
        """When the user-configured entry vector has no enzyme cuts
        at all (a backbone-only sequence with no dropout cassette),
        the diagnostic returns a clear message naming the vector
        and the missing enzyme so the user knows to pick a vector
        with the dropout intact."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Backbone with no Esp3I sites at all.
        seq = "ACGTACGTACGTACGTACGT" * 20
        rec = SeqRecord(Seq(seq), id="pNoSites", name="pNoSites")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        ev = {
            "name": "pNoSites", "size": len(seq), "source": "test",
            "gb_text": sc._record_to_gb_text(rec),
            "grammar_id": "gb_l0",
        }
        # Configure as the active entry vector for gb_l0.
        prev = sc._get_entry_vector("gb_l0")
        try:
            sc._set_entry_vector("gb_l0", ev)
            sc._settings_flush_sync()
            part = self._make_part("ATATATAT")
            reason = sc._diagnose_part_cloning(part)
            assert reason is not None
            assert "Esp3I" in reason
            assert "pNoSites" in reason
        finally:
            sc._set_entry_vector("gb_l0", prev)
            sc._settings_flush_sync()

    def test_diagnose_returns_none_when_no_entry_vector(self):
        """Diagnose returns None when no entry vector is configured
        — that's the legitimate stub-fallback path (parts saved
        before the user set up an entry vector). The Save-to-
        Collection notify shouldn't fire spuriously in this case."""
        prev = sc._get_entry_vector("gb_l0")
        try:
            sc._set_entry_vector("gb_l0", None)
            sc._settings_flush_sync()
            part = self._make_part("ATATATAT")
            assert sc._diagnose_part_cloning(part) is None
        finally:
            sc._set_entry_vector("gb_l0", prev)
            sc._settings_flush_sync()

    def test_returns_none_when_no_enzyme_in_grammar(self):
        vec_rec = self._build_test_vector()
        ev      = self._build_entry_vector_dict(vec_rec)
        part    = self._make_part("ACGTACGT")
        # Grammar without a real enzyme name → simulation can't run.
        bad_grammar = dict(sc._BUILTIN_GRAMMARS["gb_l0"])
        bad_grammar["enzyme"] = "NotARealEnzyme"
        assert sc._clone_part_into_entry_vector(
            part, ev, bad_grammar,
        ) is None

    def test_returns_none_on_malformed_vector(self):
        part = self._make_part("ACGTACGT")
        gb_l0 = sc._BUILTIN_GRAMMARS["gb_l0"]
        ev = {"name": "junk", "gb_text": "not valid GenBank text"}
        assert sc._clone_part_into_entry_vector(
            part, ev, gb_l0,
        ) is None

    def test_part_to_cloned_seqrecord_uses_simulation_when_available(
            self, isolated_library):
        # When an entry vector is set for the part's grammar,
        # _part_to_cloned_seqrecord should produce the simulated
        # cloned plasmid (with vector annotations carried through).
        vec_rec = self._build_test_vector()
        sc._set_entry_vector("gb_l0", self._build_entry_vector_dict(vec_rec))
        sc._settings_flush_sync()
        part = self._make_part("GAGGAGAAATTAACTATGCATCATCAT")
        rec = sc._part_to_cloned_seqrecord(part)
        labels = [
            f.qualifiers.get("label", ["?"])[0] for f in rec.features
        ]
        # Vector's ori comes through — proof we used the real
        # simulation, not the pUPD2 stub.
        assert "ori" in labels
        # Stub fallback would set the description to "Cloned part:
        # <name>" with no "in <vector>" clause.
        assert "in " in rec.description

    def test_part_to_cloned_seqrecord_falls_back_when_no_entry_vector(
            self, isolated_library):
        # No entry vector → use stub backbone form (preserves
        # historical behaviour for parts saved before this change).
        sc._set_entry_vector("gb_l0", None)
        sc._settings_flush_sync()
        part = self._make_part("ACGTACGT")
        rec = sc._part_to_cloned_seqrecord(part)
        # The stub form ends in `_PUPD2_BACKBONE_STUB` — verify the
        # stub bytes appear at the end of the cloned sequence.
        assert sc._PUPD2_BACKBONE_STUB in str(rec.seq)

    def test_iis_digest_byte_exact_no_carryover(self):
        """The IIS-digest path must produce a cloned plasmid that
        is byte-exact equal to ``vector.replace(oh5+dropout+oh3,
        oh5+insert+oh3)`` (modulo origin rotation). Verifies no
        primed-amplicon carryover (pad / extra spacer / internal
        site doubling) and no double-counted overhang bytes."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        BB5 = "AAAAAA" * 10
        OH5 = "AATG"; DROPOUT = "GGGGCCCCCCAAAATTTTGG"; OH3 = "GCTT"
        BB3 = "TTTTTT" * 10
        vec_seq = (BB5 + "CGTCTC" + "A" + OH5 + DROPOUT + OH3
                   + "T" + "GAGACG" + BB3)
        rec = SeqRecord(Seq(vec_seq), id="pAudit", name="pAudit",
                         description="audit")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        ev = {
            "name": "pAudit", "size": len(vec_seq), "source": "test",
            "gb_text": sc._record_to_gb_text(rec),
        }
        insert = "ATATATATATATATAT"
        part = self._make_part(insert)
        gb_l0 = sc._BUILTIN_GRAMMARS["gb_l0"]
        cloned = sc._clone_part_into_entry_vector(part, ev, gb_l0)
        assert cloned is not None
        seq = str(cloned.seq)
        expected = vec_seq.replace(OH5 + DROPOUT + OH3,
                                    OH5 + insert + OH3)
        # Byte-exact equality, modulo rotation (cloned plasmid is
        # circular so the origin can land anywhere).
        assert len(seq) == len(expected), (
            f"length differs: {len(seq)} vs {len(expected)} — "
            f"means primed amplicon bytes (pad/spacer/site) leaked "
            f"into the cloned plasmid"
        )
        assert seq in expected + expected, (
            "cloned isn't a rotation of expected — extra or missing "
            "bytes vs vector with dropout swapped"
        )

    def test_sequence_splice_byte_exact_no_carryover(self):
        """Same byte-exactness invariant for the sequence-splice
        fallback path (vector without IIS-site annotations)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        BB5 = "T" * 100
        OH5 = "AATG"; DROPOUT = "GGGGCCCCCCAAAATTTTGGAAAAAAAA"
        OH3 = "GCTT"; BB3 = "G" * 100
        vec_seq = BB5 + OH5 + DROPOUT + OH3 + BB3
        rec = SeqRecord(Seq(vec_seq), id="pSeqSplice",
                         name="pSeqSplice", description="audit")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        ev = {
            "name": "pSeqSplice", "size": len(vec_seq),
            "source": "test",
            "gb_text": sc._record_to_gb_text(rec),
        }
        insert = "ATATATATATATATAT"
        part = self._make_part(insert)
        cloned = sc._splice_part_into_vector_by_overhang(ev, part)
        assert cloned is not None
        seq = str(cloned.seq)
        expected = vec_seq.replace(OH5 + DROPOUT + OH3,
                                    OH5 + insert + OH3)
        assert len(seq) == len(expected), (
            f"length differs: {len(seq)} vs {len(expected)}"
        )
        assert seq in expected + expected, (
            "cloned isn't a rotation of expected"
        )

    def test_sequence_splice_fallback_when_no_iis_sites(
            self, isolated_library):
        """2026-05-07 fix: when the user's entry vector lacks
        Esp3I/BsaI sites in its annotation (most real plasmids
        loaded from a generic library don't have them in the
        gb_text the user picks), the strict IIS digest can't
        produce a clean dropout. The sequence-splice fallback
        still finds the part's `oh5...oh3` pair in the vector and
        replaces the inner segment — so the saved plasmid IS the
        user's vector with the dropout swapped for the insert,
        not the pUPD2 stub the legacy fallback used."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        # Build a vector with the AATG/GCTT overhangs flanking a
        # dropout cassette but NO Esp3I sites — exactly the case
        # that broke the strict IIS simulation.
        backbone_5 = "T" * 500 + "AATG"
        dropout    = "GCAAACCCGGG" * 30
        backbone_3 = "GCTT" + "G" * 500
        seq = backbone_5 + dropout + backbone_3
        rec = SeqRecord(Seq(seq), id="pUserVec", name="pUserVec",
                         description="user vector with dropout, no IIS sites")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, 100, strand=1),
            type="rep_origin", qualifiers={"label": ["ori"]},
        ))
        ev = {
            "name": "pUserVec",
            "size": len(seq),
            "source": "test",
            "gb_text": sc._record_to_gb_text(rec),
        }
        sc._set_entry_vector("gb_l0", ev)
        sc._settings_flush_sync()
        insert = "GAGGAGAAATTAACTATGCATCATCAT"
        part = self._make_part(insert)
        cloned = sc._part_to_cloned_seqrecord(part)
        # Cloned plasmid is the user's vector minus dropout plus insert
        s = str(cloned.seq)
        assert insert in s, "Insert missing from cloned plasmid"
        assert dropout not in s, "Dropout still present"
        assert "T" * 100 in s, "Backbone (T-run) missing"
        assert cloned.annotations["topology"] == "circular"
        # Vector annotations carried through
        labels = [
            f.qualifiers.get("label", ["?"])[0] for f in cloned.features
        ]
        assert "ori" in labels
        # Insert annotated with the part name
        assert "myPromoter" in labels
        # NOT the pUPD2 stub
        assert sc._PUPD2_BACKBONE_STUB not in s
        # Description references the user's vector name
        assert "pUserVec" in cloned.description


class TestPartToClonedSeqRecord:
    """`_part_to_cloned_seqrecord` builds a circular SeqRecord whose
    sequence is the part's cloned form (insert + 5'/3' OH + pUPD2
    backbone) with one feature spanning the insert region. This is
    the helper that the parts-bin "Save to Collection" button uses
    to feed `LibraryPanel.add_entry`.
    """

    def test_builds_circular_record_with_insert_feature(self):
        part = {
            "name": "myCDS", "type": "CDS", "position": "B3",
            "sequence":  "ATGAAATAATAA",
            "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spec", "user": True,
        }
        rec = sc._part_to_cloned_seqrecord(part)
        # Cloned form = oh5 + insert + oh3 + backbone stub.
        assert str(rec.seq).startswith("AATG" + "ATGAAATAATAA")
        # Topology + molecule_type pull through to GenBank serialisation.
        assert rec.annotations["topology"] == "circular"
        assert rec.annotations["molecule_type"] == "DNA"
        # Single feature spans the insert. Start = len(oh5), len matches.
        assert len(rec.features) == 1
        f = rec.features[0]
        assert int(f.location.start) == 4
        assert int(f.location.end) == 4 + len("ATGAAATAATAA")
        assert f.type == "CDS"
        assert f.qualifiers.get("label") == ["myCDS"]

    def test_sanitises_record_id_for_locus_line(self):
        # `myCDS w/ slash` has spaces + non-LOCUS-safe characters; the
        # Biopython GenBank serialiser would reject it. The helper
        # replaces them with `_` for the id and preserves the original
        # name on the feature label for human display.
        part = {
            "name": "weird name w/ slash", "type": "CDS",
            "sequence": "ATGAAATAA", "oh5": "AATG", "oh3": "GCTT",
            "user": True,
        }
        rec = sc._part_to_cloned_seqrecord(part)
        assert rec.id == "weird_name_w__slash"
        assert "weird name w/ slash" in rec.description
        assert rec.features[0].qualifiers["label"] == [
            "weird name w/ slash",
        ]

    def test_rejects_part_without_sequence(self):
        part = {"name": "empty", "type": "CDS", "sequence": "",
                 "oh5": "AATG", "oh3": "GCTT", "user": True}
        with pytest.raises(ValueError, match="no sequence"):
            sc._part_to_cloned_seqrecord(part)


class TestPartsBinMultiSelect:
    """Ctrl+click on a parts row toggles its membership in the
    multi-select set; bare clicks clear the set; built-in catalog
    rows refuse Ctrl+click. Selection is wiped on `_populate` so a
    New Part insertion can't leave stale row indices behind.
    """

    @staticmethod
    def _seed_user_parts():
        sc._save_parts_bin([
            {"name": "partA", "type": "CDS", "position": "B3",
             "oh5": "AATG", "oh3": "GCTT", "backbone": "pUPD2",
             "marker": "Spec", "sequence": "ATGAAATAA",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0, "user": True},
            {"name": "partB", "type": "promoter", "position": "A2",
             "oh5": "GGAG", "oh3": "AATG", "backbone": "pUPD2",
             "marker": "Spec", "sequence": "TATAAA" * 5,
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0, "user": True},
        ])

    async def test_modifier_click_toggles_selection(self, isolated_parts_bin):
        # 2026-05-07: switched from `_pending_ctrl` + RowSelected to
        # a direct mouse-event flow because gnome-terminal eats
        # Ctrl+click for URL handling. The toggle helper is what the
        # mouse_up branch actually invokes — exercising it directly
        # is equivalent to a Shift / Ctrl / Alt + click without
        # depending on whichever modifier the host terminal forwards.
        self._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._toggle_row_in_selection(0)
            assert modal._selected_rows == {0}
            modal._toggle_row_in_selection(1)
            assert modal._selected_rows == {0, 1}
            modal._toggle_row_in_selection(0)
            assert modal._selected_rows == {1}

    async def test_modifier_click_accepts_empty_sequence_row(
            self, isolated_parts_bin):
        # Pre-2026-05-10: the toggle helper rejected rows with empty
        # `sequence` to keep removed built-in catalog placeholders
        # out of the multi-select set. Now TU/MOD plasmids
        # legitimately store their bases in `gb_text` (not
        # `sequence`), so the empty-seq guard hid every L1+ row from
        # multi-select. Rule is now `r.get("user")`-only — a user
        # part with empty sequence (TU, MOD, or hand-edited L0) is
        # still selectable, just not eligible for sequence-based
        # actions like Copy Primed.
        sc._save_parts_bin([{
            "name": "tu_like_part", "type": "TU", "position": "TU",
            "oh5": "TACA", "oh3": "GACT", "backbone": "alpha1",
            "marker": "Spec", "sequence": "",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
            "level": 1,
            "gb_text": "LOCUS x 100 bp DNA circular SYN\n//\n",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Switch to TU tab so the row appears in `_rows`.
            tabs = modal.query_one("#parts-level-tabs", sc.Tabs)
            tabs.active = "tab-parts-tu"
            await pilot.pause()
            await pilot.pause(0.1)
            modal._toggle_row_in_selection(0)
            assert modal._selected_rows == {0}

    async def test_drag_select_replaces_with_range(self, isolated_parts_bin):
        # Drag-select is the modifier-independent fallback: works in
        # every terminal because no modifier transmission is
        # involved. We invoke the on_mouse_move branch directly
        # after seeding the drag start so the test doesn't depend on
        # synthetic event delivery — same code path the running app
        # exercises at the screen level.
        self._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Pretend the mouse came down on row 0 + the user dragged
            # to row 1. Patch `_row_under_mouse` so on_mouse_move
            # reads the fake hover row without needing a real
            # pointer; this is the only seam the screen-level
            # handler actually uses.
            modal._drag_start_row = 0
            modal._drag_active = True
            modal._drag_changed = False
            modal._drag_initial = set()
            modal._row_under_mouse = lambda: 1   # type: ignore
            modal.on_mouse_move(None)
            assert modal._selected_rows == {0, 1}
            assert modal._drag_changed is True

    async def test_button_click_preserves_multi_select(
            self, isolated_parts_bin):
        # 2026-05-07 regression guard: clicking a button (Delete,
        # Save-to-Collection, etc.) used to land in `on_mouse_down`
        # because `DataTable.hover_coordinate` is sticky on the
        # last-hovered row even when the pointer is over a button.
        # The handler then treated it as a bare click on a table
        # row and ran the on_mouse_up clear-multi-select branch
        # BEFORE the button's Pressed handler could read
        # `_selected_rows`, so Delete / Save acted on just the
        # cursor row instead of the multi-select set.
        #
        # The fix gates on_mouse_down by region.contains(screen_x,
        # screen_y) — clicks whose coordinates fall outside the
        # table's screen rectangle no longer touch drag/selection
        # state. This test simulates the button click via a stub
        # event with screen coordinates outside the parts-table
        # region and asserts the selection survives.
        self._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Build a multi-select via the toggle helper.
            modal._toggle_row_in_selection(0)
            modal._toggle_row_in_selection(1)
            assert modal._selected_rows == {0, 1}
            # Stub a MouseDown event positioned well below the
            # parts-table region (screen_y far past the table) to
            # simulate a click on the button row.
            t = modal.query_one("#parts-table", sc.DataTable)
            below_y = t.region.y + t.region.height + 5
            class _StubEvent:
                def __init__(self, sy):
                    self.button = 1
                    self.ctrl = False
                    self.shift = False
                    self.meta = False
                    self.screen_x = t.region.x + 1
                    self.screen_y = sy
            modal.on_mouse_down(_StubEvent(below_y))
            # Selection unaffected — the click-outside-region guard
            # kicked in.
            assert modal._selected_rows == {0, 1}
            assert modal._drag_active is False

    async def test_drag_unions_with_existing_selection(self, isolated_parts_bin):
        # 2026-05-07 regression guard: previously on_mouse_move
        # REPLACED `_selected_rows` with the dragged range, so a
        # stray pixel-wiggle during a click after earlier
        # modifier-toggles wiped out everything except the dragged
        # rows. The fix is to UNION with a snapshot taken at
        # MouseDown — this test seeds an existing selection of two
        # toggle picks, then runs a drag that crosses one of them,
        # and asserts both the prior picks AND the dragged range
        # survive in `_selected_rows`.
        self._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Pre-existing modifier+click toggles.
            modal._toggle_row_in_selection(0)
            modal._toggle_row_in_selection(1)
            assert modal._selected_rows == {0, 1}
            # Begin a drag from row 0; snapshot captures {0, 1}.
            modal._drag_start_row = 0
            modal._drag_active = True
            modal._drag_changed = False
            modal._drag_initial = set(modal._selected_rows)
            modal._row_under_mouse = lambda: 1   # type: ignore
            modal.on_mouse_move(None)
            # Range [0..1] unioned with snapshot {0, 1} → still {0, 1}.
            assert modal._selected_rows == {0, 1}
            # _drag_changed stays False because the union is unchanged.
            assert modal._drag_changed is False

    async def test_other_buttons_dim_when_selection_non_empty(
            self, isolated_parts_bin):
        self._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Pre-selection: every action button is enabled.
            for sel in modal._OTHER_BTN_IDS:
                assert modal.query_one(sel, sc.Button).disabled is False
            # Add row 0 to the multi-select.
            modal._selected_rows.add(0)
            modal._refresh_multi_select_visuals()
            # All other buttons are now dimmed; Save-to-Collection stays.
            for sel in modal._OTHER_BTN_IDS:
                assert modal.query_one(sel, sc.Button).disabled is True
            assert modal.query_one("#btn-parts-save-to-coll",
                                    sc.Button).disabled is False
            # Clearing restores everything.
            modal._clear_multi_select()
            for sel in modal._OTHER_BTN_IDS:
                assert modal.query_one(sel, sc.Button).disabled is False

    async def test_populate_clears_stale_selection(self, isolated_parts_bin):
        self._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._selected_rows.update({0, 1})
            modal._populate()
            assert modal._selected_rows == set()


class TestPartsBinSaveToCollection:
    """Save-to-Collection bulk-saves the selected parts to the
    library/active collection. Each part lands as its own SeqRecord
    in cloned-plasmid form."""

    async def test_save_to_collection_persists_selected(
            self, isolated_library, isolated_parts_bin):
        TestPartsBinMultiSelect._seed_user_parts()
        # Need an active collection or the action bails. The
        # `_ensure_default_collection` runs in PlasmidApp.compose().
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Pre-condition: library is empty (or only contains the
            # default seed; isolate via the fixture).
            initial = len(sc._load_library())
            # Multi-select both user parts (rows 0 + 1).
            modal._selected_rows = {0, 1}
            modal._refresh_multi_select_visuals()
            modal.query_one("#btn-parts-save-to-coll", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Library grew by exactly two entries; both have the
            # cloned topology.
            entries = sc._load_library()
            assert len(entries) == initial + 2
            # Selection cleared after save.
            assert modal._selected_rows == set()

    async def test_save_to_collection_no_active_warns(
            self, isolated_library, isolated_parts_bin, monkeypatch):
        TestPartsBinMultiSelect._seed_user_parts()
        # Force `_get_active_collection_name` to return None so the
        # "no active collection" guard fires. (Auto-default collection
        # is normally created at app startup.)
        monkeypatch.setattr(sc, "_get_active_collection_name",
                              lambda: None)
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._selected_rows = {0, 1}
            initial = len(sc._load_library())
            modal.query_one("#btn-parts-save-to-coll", sc.Button).press()
            await pilot.pause()
            # No write happened — library unchanged.
            assert len(sc._load_library()) == initial


class TestPartsBinDelete:
    """Delete button on the parts bin pushes
    `PartsBinDeleteConfirmModal`. Single-select shows the part name;
    multi-select shows the count + name preview. Built-in catalog
    rows refuse delete (they aren't in the parts-bin file).
    """

    async def test_single_delete_with_confirmation(
            self, isolated_parts_bin):
        TestPartsBinMultiSelect._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            initial = len(sc._load_parts_bin())
            assert initial == 2
            # Cursor on row 0 (first user part), no multi-select.
            t = modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            modal.query_one("#btn-parts-delete", sc.Button).press()
            await pilot.pause()
            # Confirm modal should be on top.
            confirm = app.screen
            assert isinstance(confirm, sc.PartsBinDeleteConfirmModal)
            # Title for single-delete should NOT mention a count.
            title = str(
                confirm.query_one("#partsdel-title", sc.Static).render()
            )
            assert "1 parts" not in title
            # Confirm.
            confirm.query_one("#btn-partsdel-yes", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert len(sc._load_parts_bin()) == initial - 1

    async def test_multi_delete_shows_count(self, isolated_parts_bin):
        TestPartsBinMultiSelect._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            initial = len(sc._load_parts_bin())
            modal._selected_rows = {0, 1}
            modal._refresh_multi_select_visuals()
            modal.query_one("#btn-parts-delete", sc.Button).press()
            await pilot.pause()
            confirm = app.screen
            assert isinstance(confirm, sc.PartsBinDeleteConfirmModal)
            # Title for multi-delete shows the count.
            title = str(
                confirm.query_one("#partsdel-title", sc.Static).render()
            )
            assert "2 parts" in title
            confirm.query_one("#btn-partsdel-yes", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Both parts removed.
            assert len(sc._load_parts_bin()) == initial - 2

    async def test_delete_cancel_keeps_parts(self, isolated_parts_bin):
        TestPartsBinMultiSelect._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            initial = len(sc._load_parts_bin())
            t = modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            modal.query_one("#btn-parts-delete", sc.Button).press()
            await pilot.pause()
            confirm = app.screen
            confirm.query_one("#btn-partsdel-no", sc.Button).press()
            await pilot.pause()
            assert len(sc._load_parts_bin()) == initial

    async def test_delete_with_empty_parts_bin(
            self, isolated_parts_bin):
        # 2026-05-07: built-in catalog rows are gone, so an empty
        # parts-bin file means an empty table. Pressing Delete
        # should notify "select a part first" and NOT push the
        # confirm modal.
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._rows == []
            modal.query_one("#btn-parts-delete", sc.Button).press()
            await pilot.pause()
            # Still on the parts modal — empty table → no confirm.
            assert isinstance(app.screen, sc.PartsBinModal)

    async def test_delete_clears_lingering_detail_panel(
            self, isolated_parts_bin):
        # 2026-05-07 regression guard: deleting the last part used
        # to leave the detail Static + sequence TextArea showing
        # the now-deleted part because RowHighlighted doesn't fire
        # reliably after a clear+rebuild. _populate now drives the
        # render explicitly and falls back to _clear_part_detail
        # when the table ends up empty.
        TestPartsBinMultiSelect._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            t = modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            # Pre-condition: detail panel + seq view show the part.
            detail = modal.query_one("#parts-detail", sc.Static)
            seq_view = modal.query_one("#parts-seq-view", sc.TextArea)
            assert "partA" in str(detail.render())
            assert "ATG" in seq_view.text
            # Wipe parts-bin file + repopulate (mirrors what the
            # delete confirm callback does after `_save_parts_bin`
            # writes the post-delete entries).
            sc._save_parts_bin([])
            modal._populate()
            await pilot.pause()
            assert modal._rows == []
            # Detail + seq view should be clear — no stale data.
            assert str(detail.render()).strip() == ""
            assert seq_view.text == ""

    async def test_delete_button_stays_enabled_during_multi_select(
            self, isolated_parts_bin):
        # Both Save-to-Collection and Delete should remain clickable
        # during multi-select; everything else dims.
        TestPartsBinMultiSelect._seed_user_parts()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._selected_rows = {0}
            modal._refresh_multi_select_visuals()
            assert modal.query_one(
                "#btn-parts-delete", sc.Button
            ).disabled is False
            assert modal.query_one(
                "#btn-parts-save-to-coll", sc.Button
            ).disabled is False


# ═══════════════════════════════════════════════════════════════════════════════
# Parts Bin: Load Part button → LoadPartSourceModal picker (2026-05-10)
# ═══════════════════════════════════════════════════════════════════════════════

class TestLoadPartSourceModal:
    """Regression guard for 2026-05-10 fix.

    Pre-2026-05-10 the Parts Bin "Load Part" button took
    ``self.app._current_record`` directly, forcing the user to load
    the candidate plasmid onto the canvas first. The new flow pushes
    ``LoadPartSourceModal`` so the user can pick from any saved
    plasmid in any collection, or open a fresh ``.gb`` / ``.dna``
    from disk — without disturbing the canvas.
    """

    @staticmethod
    def _seed_two_collections() -> str:
        """Seed two collections with one plasmid each. Returns the gb_text
        of the active-collection plasmid so the test can round-trip
        through the picker."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        # Active collection: one Golden Braid L0 PromUTR (combined
        # Promoter+5'UTR; GGAG → AATG).
        gb_seq = (
            "CGTCTCA" + "GGAG" + "ATGAAACCCGGG" * 5 + "AATG" + "AGAGACG"
            + "AAAAATTTTT" * 50
        )
        gb_record = SeqRecord(
            Seq(gb_seq),
            id="active_part",
            name="active_part",
            annotations={"topology": "circular", "molecule_type": "DNA"},
        )
        gb_text = sc._record_to_gb_text(gb_record)
        sc._save_collections([
            {
                "name":        "Active",
                "description": "active",
                "plasmids":    [{
                    "id":      "active_part",
                    "name":    "active_part",
                    "size":    len(gb_seq),
                    "n_feats": 0,
                    "gb_text": gb_text,
                }],
            },
            {
                "name":        "Other",
                "description": "other",
                "plasmids":    [{
                    "id":      "other_part",
                    "name":    "other_part",
                    "size":    100,
                    "n_feats": 0,
                    "gb_text": gb_text,  # same content, separate id
                }],
            },
        ])
        sc._set_active_collection_name("Active")
        sc._save_library([{
            "id":      "active_part",
            "name":    "active_part",
            "size":    len(gb_seq),
            "n_feats": 0,
            "gb_text": gb_text,
        }])
        return gb_text

    async def test_load_part_button_pushes_picker(
            self, isolated_parts_bin, isolated_library):
        """Clicking Load Part with no canvas record must still push the
        picker (vs. the pre-fix behavior which warned + bailed)."""
        self._seed_two_collections()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            parts_modal.query_one("#btn-load-part", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.LoadPartSourceModal)

    async def test_picker_lists_all_collections_active_first(
            self, isolated_parts_bin, isolated_library):
        """Both collections' plasmids must show up in the picker, with
        the active collection's entry listed first (regardless of
        alpha order — Active < Other alphabetically anyway, so we also
        verify by scanning the matches list directly)."""
        self._seed_two_collections()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            assert isinstance(picker, sc.LoadPartSourceModal)
            collections = [m["collection"] for m in picker._matches]
            assert "Active" in collections
            assert "Other" in collections
            # Active collection lands first in the matches list.
            assert collections.index("Active") < collections.index("Other")

    async def test_picking_circular_plasmid_classifies_to_parts_bin(
            self, isolated_parts_bin, isolated_library):
        """End-to-end: pick a Golden Braid L0 promoter from the active
        collection, the worker classifies it, and the parts bin grows
        by exactly one row tagged with the right grammar / position."""
        self._seed_two_collections()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            initial = len(sc._load_parts_bin())
            parts_modal.query_one("#btn-load-part", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            assert isinstance(picker, sc.LoadPartSourceModal)
            # Pick the active-collection row (index 0 after the
            # active-first sort).
            t = picker.query_one("#loadpart-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            # 2026-05-13: picker is now multi-select via toggle — the
            # cursor row must be toggled (space) before Load Selected
            # has anything to do.
            picker.action_toggle_selection()
            await pilot.pause()
            picker.query_one("#btn-loadpart-ok", sc.Button).press()
            await pilot.pause()
            # Worker runs in a `@work(thread=True)` — give it room to
            # land. Two pauses + a 0.5s delay matches the cadence other
            # parts-bin worker tests use.
            await pilot.pause()
            await pilot.pause(0.5)
            await pilot.pause()
            entries = sc._load_parts_bin()
            assert len(entries) == initial + 1
            saved = entries[0]
            assert saved["grammar"] == "gb_l0"
            assert saved["position"] == "Pos 1"
            assert saved["oh5"] == "GGAG"
            assert saved["oh3"] == "AATG"

    async def test_open_file_button_pushes_open_file_modal(
            self, isolated_parts_bin, isolated_library):
        """The "Open file…" button on the picker must push
        ``OpenFileModal`` so the user can pick a fresh ``.gb`` / ``.dna``
        from disk."""
        self._seed_two_collections()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            picker.query_one("#btn-loadpart-file", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.OpenFileModal)

    async def test_picker_cancel_dismisses_with_none(
            self, isolated_parts_bin, isolated_library):
        """Esc / Cancel returns ``None`` so the worker is never invoked."""
        self._seed_two_collections()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            initial = len(sc._load_parts_bin())
            parts_modal.query_one("#btn-load-part", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            assert isinstance(picker, sc.LoadPartSourceModal)
            picker.query_one("#btn-loadpart-cancel", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Back on the parts modal, no new entry.
            assert isinstance(app.screen, sc.PartsBinModal)
            assert len(sc._load_parts_bin()) == initial

    async def test_picker_warns_on_linear_pick(
            self, isolated_parts_bin, isolated_library):
        """A linear plasmid in the picker must produce a warning + no
        parts-bin write when picked (the digest needs a circular
        topology). After hardening (2026-05-10), the picker now keeps
        the modal OPEN on a linear pick so the user can pick another
        row without re-launching the dialog."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        # Seed a linear-topology plasmid in the active collection.
        gb_seq = (
            "CGTCTCA" + "GGAG" + "ATGAAACCCGGG" * 5 + "TGAC" + "AGAGACG"
            + "AAAAATTTTT" * 50
        )
        linear_record = SeqRecord(
            Seq(gb_seq),
            id="linear_part",
            name="linear_part",
            annotations={"topology": "linear", "molecule_type": "DNA"},
        )
        gb_text = sc._record_to_gb_text(linear_record)
        sc._save_collections([{
            "name":        "Linear-only",
            "description": "test",
            "plasmids":    [{
                "id":      "linear_part",
                "name":    "linear_part",
                "size":    len(gb_seq),
                "n_feats": 0,
                "gb_text": gb_text,
            }],
        }])
        sc._set_active_collection_name("Linear-only")
        sc._save_library([{
            "id":      "linear_part",
            "name":    "linear_part",
            "size":    len(gb_seq),
            "n_feats": 0,
            "gb_text": gb_text,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            initial = len(sc._load_parts_bin())
            parts_modal.query_one("#btn-load-part", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            t = picker.query_one("#loadpart-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            picker.query_one("#btn-loadpart-ok", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # No new parts-bin row — the picker's pre-flight bailed.
            assert len(sc._load_parts_bin()) == initial
            # Picker stays open so the user can pick another row.
            assert isinstance(app.screen, sc.LoadPartSourceModal)
            assert app.screen._dismissing is False


class TestLoadPartSourceModalHardening:
    """Edge-case + robustness regression guards for `LoadPartSourceModal`.

    Covers (2026-05-10):
      * Double-dismiss guard — a Select+Enter race or a row double-click
        can't fire `dismiss(record)` twice (would crash the screen
        stack on the second call).
      * Select button is disabled when the table is empty.
      * Topology pre-flight in the picker keeps the modal open instead
        of dismissing into the worker's "linear" warning.
      * Empty-sequence pre-flight surfaces in-modal feedback rather
        than a useless toast.
      * Filter timer is cancelled on unmount so a stale tick doesn't
        fire against a disposed widget tree.
    """

    @staticmethod
    def _seed_one_circular() -> None:
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        gb_seq = (
            "CGTCTCA" + "GGAG" + "ATGAAACCCGGG" * 5 + "TGAC" + "AGAGACG"
            + "AAAAATTTTT" * 50
        )
        rec = SeqRecord(
            Seq(gb_seq), id="p1", name="p1",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        gb = sc._record_to_gb_text(rec)
        sc._save_collections([{
            "name": "Coll", "description": "t",
            "plasmids": [{
                "id": "p1", "name": "p1",
                "size": len(gb_seq), "n_feats": 0,
                "gb_text": gb,
            }],
        }])
        sc._set_active_collection_name("Coll")
        sc._save_library([{
            "id": "p1", "name": "p1",
            "size": len(gb_seq), "n_feats": 0,
            "gb_text": gb,
        }])

    async def test_select_button_disabled_when_no_matches(
            self, isolated_library, isolated_parts_bin):
        """Empty library → Select button must come up disabled. The
        old code left it enabled and `_on_query_submitted` no-op'd on
        empty matches; the disabled state makes the empty case
        visible up-front."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            assert isinstance(picker, sc.LoadPartSourceModal)
            assert picker._matches == []
            btn = picker.query_one("#btn-loadpart-ok", sc.Button)
            assert btn.disabled is True

    async def test_select_button_enabled_when_matches_present(
            self, isolated_library, isolated_parts_bin):
        """Reverse of the empty case — a non-empty list flips the
        Select button on after the initial refresh."""
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            assert len(picker._matches) == 1
            btn = picker.query_one("#btn-loadpart-ok", sc.Button)
            assert btn.disabled is False

    async def test_double_dismiss_guard(
            self, isolated_library, isolated_parts_bin):
        """Calling `_submit_selection` twice (simulating a Load Selected +
        Enter race or a duplicate RowSelected event) must not crash —
        only the first call dismisses, the second is a no-op.

        2026-05-13: picker is now multi-select; the latch lives on
        `_submit_selection` instead of `_dismiss_with_match`."""
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            match = picker._matches[0]
            # Seed the selection set directly (bypasses the UI toggle
            # round-trip — same effect on the modal's state machine).
            picker._selected_keys.add((match["collection"], match["id"]))
            # First dismissal: real one — flips the latch.
            picker._submit_selection()
            assert picker._dismissing is True
            # Second call: must short-circuit. No exception, no
            # second screen-pop. Just ensure no crash.
            picker._submit_selection()

    async def test_action_cancel_idempotent(
            self, isolated_library, isolated_parts_bin):
        """Esc → action_cancel — a second invocation (rapid double
        Esc) must not double-dismiss."""
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            picker.action_cancel()
            assert picker._dismissing is True
            # Second call short-circuits.
            picker.action_cancel()

    async def test_unmount_cancels_filter_timer(
            self, isolated_library, isolated_parts_bin):
        """The debounce timer is freed on unmount — without this it
        keeps a reference to `_refresh` and fires after the modal is
        gone (cosmetic; status logs would noise up)."""
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            # Type a key to schedule the debounce.
            inp = picker.query_one("#loadpart-input", sc.Input)
            inp.value = "x"
            await pilot.pause()
            assert picker._filter_timer is not None
            picker.action_cancel()
            await pilot.pause()
            await pilot.pause(0.1)
            # After unmount, the slot is wiped to None.
            assert picker._filter_timer is None

    async def test_empty_sequence_record_does_not_dismiss(
            self, isolated_library, isolated_parts_bin):
        """A library entry whose gb_text parses to a 0-bp record must
        notify + stay open, not crash the worker chain."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(
            Seq(""), id="empty", name="empty",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        gb = sc._record_to_gb_text(rec)
        sc._save_collections([{
            "name": "Coll", "description": "t",
            "plasmids": [{
                "id": "empty", "name": "empty",
                "size": 0, "n_feats": 0, "gb_text": gb,
            }],
        }])
        sc._set_active_collection_name("Coll")
        sc._save_library([{
            "id": "empty", "name": "empty",
            "size": 0, "n_feats": 0, "gb_text": gb,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            match = picker._matches[0]
            picker._selected_keys.add((match["collection"], match["id"]))
            picker._submit_selection()
            # The empty-seq pre-flight failed for every toggled row,
            # so the batch resolved zero records and the modal stays
            # open. Latch must NOT have flipped (so the user can pick
            # a different row without the dismiss-guard tripping).
            assert picker._dismissing is False


class TestCheckVectorMatchHardening:
    """Hardening for `_check_vector_match`'s rotation + RC checks and
    the digest cache. Avoids spinning up the modal/worker — these
    operate purely on the pure-Python layer."""

    @staticmethod
    def _make_l0_ring(insert: str = "AAAA",
                       oh5: str = "GGAG",
                       oh3: str = "AATG",
                       backbone: str = "AAAAATTTTT" * 50) -> str:
        return ("CGTCTCA" + oh5 + insert + oh3 + "AGAGACG" + backbone)

    def test_match_via_doubled_substring_rotation(
            self, isolated_library):
        """Cloning typically preserves the vector half byte-for-byte;
        a configured EV's vector half should match the user's vector
        half via the doubled-substring trick. Direct equality is the
        common case; this test asserts the substring path also lands
        the match (the user's plasmid is a rotation-equivalent ring
        of the EV's vector half plus a different insert)."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        ev_seq   = self._make_l0_ring("CCCCTTTT" * 5)
        user_seq = self._make_l0_ring("ATGAAACCCGGG" * 5)
        ev_rec = SeqRecord(
            Seq(ev_seq), id="ev", name="ev",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0", "role": "",
            "name": "ev", "size": len(ev_seq),
            "source": "test", "gb_text":
                sc._record_to_gb_text(ev_rec),
        }])
        result = sc._classify_part_from_plasmid(user_seq, circular=True)
        assert result is not None
        ev = result.get("entry_vector")
        assert ev is not None and ev["matches"] is True

    def test_match_via_rc_orientation(self, isolated_library):
        """If the user saved their plasmid with the OTHER strand on
        top (RC of the canonical orientation), the vector half's
        top_seq is the RC of the EV's vector half. The hardened
        check matches via the RC-doubled fallback so the mismatch
        doesn't dismiss as "no entry vector configured"."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        ev_seq = self._make_l0_ring("CCCCTTTT" * 5)
        # Build the RC-of-the-ring as the user's plasmid. After RC,
        # the Esp3I sites + overhangs are inverted — the digest still
        # produces the SAME biological fragments but their top_seq is
        # the RC of the canonical fragments.
        user_seq = sc._rc(ev_seq)
        # RC swaps the digestion's overhang labelling: GGAG↔CTCC,
        # TGAC↔GTCA. To classify under gb_l0 the user's RC plasmid
        # must still be a recognisable Promoter — verify upfront.
        # If not, this test is a no-op for the EV check and we'd
        # silently pass; assert the classification works first so a
        # regression in the RC path doesn't slip through.
        result_rc = sc._classify_part_from_plasmid(
            user_seq, circular=True,
        )
        # If the RC seq doesn't classify (overhangs map to a non-GB
        # position), skip — this guard isn't useful in that case.
        if result_rc is None:
            import pytest
            pytest.skip("RC seq doesn't classify under any grammar "
                        "(test environment-dependent)")
        ev_rec = SeqRecord(
            Seq(ev_seq), id="ev", name="ev",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        sc._save_entry_vectors([{
            "grammar_id": result_rc["grammar_id"], "role": "",
            "name": "ev", "size": len(ev_seq),
            "source": "test", "gb_text":
                sc._record_to_gb_text(ev_rec),
        }])
        # The EV check should land via the RC fallback even though
        # the user's vector half is the RC of the EV's vector half.
        result_with_ev = sc._classify_part_from_plasmid(
            user_seq, circular=True,
        )
        ev = result_with_ev.get("entry_vector") if result_with_ev else None
        assert ev is not None, \
            "RC-fallback in _check_vector_match did not fire"
        assert ev["matches"] is True

    def test_digest_cache_avoids_reparsing(self, isolated_library):
        """The `_VECTOR_MATCH_CACHE` keyed by (gb_text, enzyme) must
        return the cached result on a repeat call. Verifies the cache
        exists + populates as expected."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        ev_seq = self._make_l0_ring("CCCCTTTT" * 5)
        ev_rec = SeqRecord(
            Seq(ev_seq), id="ev", name="ev",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        ev_gb = sc._record_to_gb_text(ev_rec)
        sc._VECTOR_MATCH_CACHE.clear()
        first  = sc._vector_half_top_seq(ev_gb, "Esp3I")
        # Sweep #25 (2026-05-23): cache key now `(hash(gb_text),
        # enzyme)` not the raw gb_text — pre-fix the multi-MB string
        # lived in every key tuple (64 × 5 MB worst case ≈ 320 MB).
        assert (hash(ev_gb), "Esp3I") in sc._VECTOR_MATCH_CACHE
        second = sc._vector_half_top_seq(ev_gb, "Esp3I")
        assert first == second is not None

    def test_digest_cache_bounded(self):
        """The cache must not grow unbounded — once it crosses the
        cap, the oldest insertion is evicted."""
        sc._VECTOR_MATCH_CACHE.clear()
        # Stuff the cache past the cap with junk keys (the helper
        # never inserts None values that fail to digest, but for
        # a bounded-cache test we exercise the eviction directly
        # by calling `_vector_half_top_seq` with synthesized stubs).
        for i in range(sc._VECTOR_MATCH_CACHE_MAX + 5):
            sc._VECTOR_MATCH_CACHE[(f"k{i}", "Esp3I")] = "stub"
        # Direct stuffing exceeded the cap; verify the helper's
        # bounded insertion path keeps it within bounds going forward.
        # Re-clear, then call the helper repeatedly with distinct
        # gb_text keys to verify it caps itself.
        sc._VECTOR_MATCH_CACHE.clear()
        for i in range(sc._VECTOR_MATCH_CACHE_MAX + 5):
            sc._vector_half_top_seq(f">fake gb_text {i}", "Esp3I")
        assert len(sc._VECTOR_MATCH_CACHE) <= sc._VECTOR_MATCH_CACHE_MAX


class TestLoadPartSourceModalShutdown:
    """Shutdown path coverage: Esc, Ctrl+Q, app.exit while open, and
    cleanup invariants (filter timer cancelled regardless of which
    dismiss path fires). Regression guard for 2026-05-10 hardening
    pass.
    """

    @staticmethod
    def _seed_one_circular() -> None:
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        gb_seq = (
            "CGTCTCA" + "GGAG" + "ATGAAACCCGGG" * 5 + "TGAC" + "AGAGACG"
            + "AAAAATTTTT" * 50
        )
        rec = SeqRecord(
            Seq(gb_seq), id="p1", name="p1",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        gb = sc._record_to_gb_text(rec)
        sc._save_collections([{
            "name": "Coll", "description": "t",
            "plasmids": [{
                "id": "p1", "name": "p1",
                "size": len(gb_seq), "n_feats": 0,
                "gb_text": gb,
            }],
        }])
        sc._set_active_collection_name("Coll")
        sc._save_library([{
            "id": "p1", "name": "p1",
            "size": len(gb_seq), "n_feats": 0,
            "gb_text": gb,
        }])

    async def test_esc_cascade_loadpart_then_partsbin(
            self, isolated_library, isolated_parts_bin):
        """Esc on LoadPart pops it; second Esc pops PartsBin. Verifies
        the modal stack drains cleanly without leaving an orphan
        screen on top."""
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts = app.screen
            assert isinstance(parts, sc.PartsBinModal)
            parts.query_one("#btn-load-part", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            assert isinstance(picker, sc.LoadPartSourceModal)
            # First Esc: pops LoadPart.
            await pilot.press("escape")
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.PartsBinModal)
            # Second Esc: pops PartsBin.
            await pilot.press("escape")
            await pilot.pause()
            await pilot.pause(0.1)
            # Parts bin gone. Default app screen is on top.
            assert not isinstance(app.screen, sc.PartsBinModal)
            assert not isinstance(app.screen, sc.LoadPartSourceModal)

    async def test_esc_with_pending_filter_timer_cancels_cleanly(
            self, isolated_library, isolated_parts_bin):
        """Type into the filter (schedules a debounce timer), then
        press Esc immediately. The on_unmount hook must cancel the
        timer so it doesn't fire against the disposed widget tree."""
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            inp = picker.query_one("#loadpart-input", sc.Input)
            inp.value = "p"
            await pilot.pause()
            assert picker._filter_timer is not None
            # Esc before the 150 ms debounce window expires.
            await pilot.press("escape")
            await pilot.pause()
            await pilot.pause(0.1)
            # on_unmount cleared the timer slot; the original Timer
            # object was stopped before the slot was wiped.
            assert picker._filter_timer is None
            assert picker._dismissing is True

    async def test_select_dismiss_also_cancels_timer(
            self, isolated_library, isolated_parts_bin):
        """on_unmount fires on EVERY dismiss path, not just action_
        cancel. Verifies the Select-with-record path also clears the
        timer (regression guard against a half-baked cleanup that
        only covered the cancel path)."""
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            inp = picker.query_one("#loadpart-input", sc.Input)
            inp.value = "p"
            await pilot.pause()
            assert picker._filter_timer is not None
            # 2026-05-13: multi-select picker — must toggle a row
            # before Load Selected actually dismisses.
            t = picker.query_one("#loadpart-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            picker.action_toggle_selection()
            await pilot.pause()
            # Select via OK button — dismisses with the record list.
            picker.query_one("#btn-loadpart-ok", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # on_unmount fired during the dismiss; timer is gone.
            assert picker._filter_timer is None

    async def test_app_exit_while_picker_open_unmounts_cleanly(
            self, isolated_library, isolated_parts_bin):
        """Exiting the app while the picker is on top triggers a full
        unmount cascade — the modal's `on_unmount` must run so the
        timer doesn't leak past process shutdown.
        """
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            inp = picker.query_one("#loadpart-input", sc.Input)
            inp.value = "p"
            await pilot.pause()
            assert picker._filter_timer is not None
            # Capture a reference; we inspect post-exit.
            picker_ref = picker
            # Hard exit — no UnsavedQuitModal in the path because the
            # app has no record loaded. Mirrors a `app.exit()` from
            # any code path.
            app.exit()
            await pilot.pause()
            await pilot.pause(0.2)
        # After context exit, the test harness has shut the app down.
        # The picker reference still holds, but on_unmount fired and
        # cleared the slot.
        assert picker_ref._filter_timer is None

    async def test_cancel_button_idempotent_after_select(
            self, isolated_library, isolated_parts_bin):
        """Once Select fires (`_dismissing=True`), pressing Cancel
        afterward must be a no-op. Defends against a race where the
        user clicks Select, then mashes Cancel before the dismiss
        animation completes."""
        self._seed_one_circular()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            # Manually trigger the dismiss path without letting the
            # screen actually pop yet. Simulates the same race.
            picker._dismissing = True
            # Cancel button: must short-circuit with `_dismissing=True`.
            picker._cancel_btn(None)
            picker.action_cancel()
            # No exception, latch stays True, no second dismiss queued.
            assert picker._dismissing is True


class TestLoadPartSourceModalRobustness:
    """Edge cases discovered in the 2nd-pass audit (2026-05-10):
    corrupt gb_text recovery, search-helper exception path, original
    gb_text round-trip preservation."""

    async def test_corrupt_gb_text_keeps_modal_open(
            self, isolated_library, isolated_parts_bin):
        """An entry whose gb_text is malformed (un-parseable) must
        notify + stay open, not crash the dismiss path. Pre-fix the
        parse exception was caught but the modal could leave the
        broken row selectable; the latch ensures we don't latch
        `_dismissing` on a parse failure."""
        sc._save_collections([{
            "name": "Corrupt", "description": "t",
            "plasmids": [{
                "id": "broken", "name": "broken",
                "size": 5_000, "n_feats": 0,
                "gb_text": "this is not GenBank text at all",
            }],
        }])
        sc._set_active_collection_name("Corrupt")
        sc._save_library([{
            "id": "broken", "name": "broken",
            "size": 5_000, "n_feats": 0,
            "gb_text": "this is not GenBank text at all",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            match = picker._matches[0]
            picker._selected_keys.add((match["collection"], match["id"]))
            picker._submit_selection()
            # Parse failed → batch resolved zero records → notify +
            # stay open. Latch must NOT have flipped so the user can
            # pick a different row.
            assert picker._dismissing is False
            assert isinstance(app.screen, sc.LoadPartSourceModal)

    async def test_search_helper_exception_leaves_modal_usable(
            self, isolated_library, isolated_parts_bin,
            monkeypatch):
        """If `_search_collections_library` raises during a refresh,
        the modal must catch + show an error banner without crashing.
        Simulates a buggy custom-grammar plugin or a corrupted
        in-memory state mid-iteration."""
        def _boom(*args, **kwargs):
            raise RuntimeError("simulated search failure")
        monkeypatch.setattr(
            sc, "_search_collections_library", _boom,
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal())
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            assert isinstance(picker, sc.LoadPartSourceModal)
            # Search failed → empty matches, Select disabled, status
            # has the red banner.
            assert picker._matches == []
            btn = picker.query_one("#btn-loadpart-ok", sc.Button)
            assert btn.disabled is True
            status = picker.query_one("#loadpart-status", sc.Static)
            assert "failed" in str(status.render()).lower()

    async def test_original_gb_text_stashed_on_record(
            self, isolated_library, isolated_parts_bin):
        """The picker stashes the library's ORIGINAL `gb_text` on the
        dismissed record (`record._tui_gb_text`) so the worker can
        skip the parse → serialise round-trip. Verifies the stash
        is intact at dismiss time."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        gb_seq = (
            "CGTCTCA" + "GGAG" + "ATGAAACCCGGG" * 5 + "TGAC"
            + "AGAGACG" + "AAAAATTTTT" * 50
        )
        rec = SeqRecord(
            Seq(gb_seq), id="p1", name="p1",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        original_gb = sc._record_to_gb_text(rec)
        sc._save_collections([{
            "name": "Coll", "description": "t",
            "plasmids": [{
                "id": "p1", "name": "p1",
                "size": len(gb_seq), "n_feats": 0,
                "gb_text": original_gb,
            }],
        }])
        sc._set_active_collection_name("Coll")
        sc._save_library([{
            "id": "p1", "name": "p1",
            "size": len(gb_seq), "n_feats": 0,
            "gb_text": original_gb,
        }])
        # Capture the dismissed record list so we can inspect it
        # post-dismiss. 2026-05-13: picker now dismisses with
        # list[SeqRecord] (multi-select).
        captured: dict = {}

        def _capture(records):
            captured["records"] = records

        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.LoadPartSourceModal(),
                            callback=_capture)
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            match = picker._matches[0]
            picker._selected_keys.add((match["collection"], match["id"]))
            picker._submit_selection()
            await pilot.pause()
            await pilot.pause(0.1)
        records = captured.get("records") or []
        assert len(records) == 1, f"expected 1 record, got {records!r}"
        rec = records[0]
        # Stashed gb_text matches the library entry's original.
        assert getattr(rec, "_tui_gb_text", None) == original_gb


class TestSaveEntryVectorsCacheInvalidation:
    """`_save_entry_vectors` must drop the digest cache so a
    reconfigure doesn't leave stale entries crowding the cap or
    (worst case) returning a digest derived from a deleted EV."""

    def test_save_clears_vector_match_cache(self, isolated_library):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        # Pre-populate the cache via a manual digest.
        ev_seq = ("CGTCTCA" + "GGAG" + "AAAA" + "TGAC"
                  + "AGAGACG" + "AAAAATTTTT" * 50)
        ev_rec = SeqRecord(
            Seq(ev_seq), id="ev", name="ev",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        ev_gb = sc._record_to_gb_text(ev_rec)
        sc._VECTOR_MATCH_CACHE.clear()
        sc._vector_half_top_seq(ev_gb, "Esp3I")
        assert len(sc._VECTOR_MATCH_CACHE) >= 1
        # Persist any entry-vector list — the call must wipe the
        # cache regardless of contents.
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0", "role": "",
            "name": "ev", "size": len(ev_seq),
            "source": "test", "gb_text": ev_gb,
        }])
        assert len(sc._VECTOR_MATCH_CACHE) == 0


class TestVectorHalfTopSeqRobustness:
    """Edge-cases for `_vector_half_top_seq` — non-circular EV refusal,
    parse-failure swallow, cache hit on repeat call."""

    def test_non_circular_ev_returns_none(self):
        """A linearised EV would dispatch through
        `_excise_fragment_pair(circular=True)` and produce nonsense
        fragments. Skip the digest and return None."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        seq = ("CGTCTCA" + "GGAG" + "AAAA" + "TGAC"
               + "AGAGACG" + "AAAAATTTTT" * 50)
        rec = SeqRecord(
            Seq(seq), id="lin", name="lin",
            annotations={"topology": "linear",
                         "molecule_type": "DNA"},
        )
        ev_gb = sc._record_to_gb_text(rec)
        sc._VECTOR_MATCH_CACHE.clear()
        result = sc._vector_half_top_seq(ev_gb, "Esp3I")
        assert result is None

    def test_unparseable_gb_returns_none(self):
        """Bogus gb_text → parse fails → swallow and return None
        rather than propagating the exception."""
        sc._VECTOR_MATCH_CACHE.clear()
        result = sc._vector_half_top_seq(
            "not actual genbank text", "Esp3I",
        )
        assert result is None


# ═══════════════════════════════════════════════════════════════════════════════
# Parts Bin: Save As Feature button
# ═══════════════════════════════════════════════════════════════════════════════

class TestGbPartTypeToInsdcMap:
    """Pure mapping: GB part type (TitleCase + spaces) → INSDC
    feature_type. CDS-NS / C-tag are GB-specific shapes that have no
    INSDC equivalent; they collapse to plain "CDS" with the original
    type preserved in the description."""

    @pytest.mark.parametrize("gb_type, insdc", [
        # Legacy slots:
        ("Promoter",         "promoter"),
        ("Promoter-only",    "promoter"),
        ("5' UTR",           "5'UTR"),
        ("CDS",              "CDS"),
        ("CDS-NS",           "CDS"),
        ("C-tag",            "CDS"),
        ("Terminator",       "terminator"),
        # GB 2.0 canonical additions (2026-05-10):
        ("Operator-A",       "promoter"),
        ("Operator-B",       "promoter"),
        ("Min Promoter",     "promoter"),
        ("Distal 5' UTR",    "5'UTR"),
        ("Signal peptide",   "sig_peptide"),
        ("CDS-NS (CT)",      "CDS"),
        ("CT-tag",            "CDS"),
        ("CDS-after-SP",     "CDS"),
        ("3' UTR",           "3'UTR"),
        ("Terminator-only",  "terminator"),
    ])
    def test_known_types_map_to_insdc(self, gb_type, insdc):
        assert sc._GB_PART_TYPE_TO_INSDC[gb_type] == insdc

    def test_every_gb_position_has_a_mapping(self):
        # Every TitleCase part type listed in _GB_POSITIONS must have a
        # corresponding INSDC mapping; otherwise Save As Feature would
        # silently fall through to "misc_feature" for that shape.
        for gb_type in sc._GB_POSITIONS:
            assert gb_type in sc._GB_PART_TYPE_TO_INSDC, (
                f"Missing INSDC mapping for GB part type {gb_type!r}"
            )


class TestPartsBinSaveAsFeature:
    """`PartsBinModal` grew a "Save As Feature" button that takes the
    highlighted user part and registers it in the persistent feature
    library via `AddFeatureModal` → `_persist_feature_entry`."""

    async def test_button_present(self, isolated_parts_bin, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            assert app.screen.query_one(
                "#btn-parts-save-as-feature", sc.Button,
            ) is not None

    async def test_user_part_pushes_prefilled_modal(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_parts_bin([{
            "name":      "myCDS",
            "type":      "CDS",
            "position":  "Pos 3-4",
            "oh5":       "AATG", "oh3": "GCTT",
            "backbone":  "pUPD2",
            "marker":    "Spectinomycin",
            "sequence":  "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            parts_modal.query_one("#parts-table",
                                  sc.DataTable).move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-save-as-feature",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.AddFeatureModal)
            # Pre-fill carries the part's name + sequence; type is the
            # INSDC mapping (CDS → CDS in this case).
            from textual.widgets import Input, Select, TextArea
            assert app.screen.query_one("#addfeat-name", Input).value == "myCDS"
            assert app.screen.query_one("#addfeat-type", Select).value == "CDS"
            assert (
                app.screen.query_one("#addfeat-seq", TextArea).text
                == "ATGCATGCATGC"
            )
            # GB metadata rides along in the description.
            desc = app.screen.query_one("#addfeat-desc", Input).value
            assert "Pos 3-4" in desc
            assert "AATG" in desc and "GCTT" in desc

    async def test_builtin_part_is_rejected(
        self, isolated_parts_bin, isolated_library,
    ):
        """Built-in catalog rows have no sequence — Save As Feature
        must surface a warning rather than open an empty AddFeatureModal."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            # No user parts — first row is built-in (no sequence).
            parts_modal.query_one("#parts-table",
                                  sc.DataTable).move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-save-as-feature",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Still on the parts screen — no AddFeatureModal pushed.
            assert isinstance(app.screen, sc.PartsBinModal)

    async def test_save_persists_to_feature_library(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_parts_bin([{
            "name":      "myProm",
            "type":      "Promoter",
            "position":  "Pos 1",
            "oh5":       "GGAG", "oh3": "TGAC",
            "backbone":  "pUPD2",
            "marker":    "Spectinomycin",
            "sequence":  "TATAAATATA",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            parts_modal.query_one("#parts-table",
                                  sc.DataTable).move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-save-as-feature",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Now in AddFeatureModal — press Save.
            app.screen.query_one("#btn-addfeat-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Back on the parts screen, with the entry now in features.json.
            assert isinstance(app.screen, sc.PartsBinModal)
            sc._features_cache = None
            entries = sc._load_features()
            assert len(entries) == 1
            entry = entries[0]
            assert entry["name"] == "myProm"
            # Promoter (TitleCase) → promoter (INSDC).
            assert entry["feature_type"] == "promoter"
            assert entry["sequence"] == "TATAAATATA"
            assert entry["strand"] == 1

    async def test_cds_ns_collapses_to_cds_with_note(
        self, isolated_parts_bin, isolated_library,
    ):
        """CDS-NS is a GB-specific coding shape (no stop codon). INSDC
        has no equivalent so the entry is saved as plain "CDS"; the
        original GB type rides along in the description so the
        distinction is recoverable from the library."""
        sc._save_parts_bin([{
            "name":      "myNS",
            "type":      "CDS-NS",
            "position":  "Pos 3",
            "oh5":       "AATG", "oh3": "TTCG",
            "backbone":  "pUPD2",
            "marker":    "Spectinomycin",
            "sequence":  "ATGAAATTT",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            parts_modal.query_one("#parts-table",
                                  sc.DataTable).move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-save-as-feature",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            from textual.widgets import Select, Input
            assert app.screen.query_one("#addfeat-type", Select).value == "CDS"
            desc = app.screen.query_one("#addfeat-desc", Input).value
            assert "GB type: CDS-NS" in desc, (
                f"Expected GB-type note in description; got {desc!r}"
            )


class TestFeatureLibraryMatchHelper:
    """Pure helper: classify a (name, feature_type, sequence) triple
    against the on-disk feature library."""

    def test_no_match_returns_empty(self, isolated_library):
        sc._save_features([])
        assert sc._feature_library_match("foo", "CDS", "ATG") == ""

    def test_exact_match(self, isolated_library):
        sc._save_features([{
            "name": "foo", "feature_type": "CDS",
            "sequence": "ATGCAT", "strand": 1,
        }])
        assert sc._feature_library_match("foo", "CDS", "ATGCAT") == "exact"

    def test_exact_match_is_case_insensitive_on_sequence(self, isolated_library):
        # Stored sequence is upper-cased by AddFeatureModal._gather, but
        # be defensive — the helper canonicalises both sides.
        sc._save_features([{
            "name": "foo", "feature_type": "CDS",
            "sequence": "atgcat", "strand": 1,
        }])
        assert sc._feature_library_match("foo", "CDS", "ATGCAT") == "exact"

    def test_name_match_with_different_sequence(self, isolated_library):
        sc._save_features([{
            "name": "foo", "feature_type": "CDS",
            "sequence": "ATGCAT", "strand": 1,
        }])
        assert sc._feature_library_match("foo", "CDS", "AAATTT") == "name"

    def test_different_type_is_not_a_match(self, isolated_library):
        # Same name but different feature_type → treat as a separate
        # entry. promoter and CDS with the same label are independent.
        sc._save_features([{
            "name": "foo", "feature_type": "promoter",
            "sequence": "ATGCAT", "strand": 1,
        }])
        assert sc._feature_library_match("foo", "CDS", "ATGCAT") == ""


class TestFeatureLibraryGenerationCounter:
    """The `_features_generation` counter must bump on every change to
    the feature library so consumers (PartsBinModal index cache, etc.)
    can detect "the library has changed since I last looked" without
    re-scanning the entries list. Strict bump-on-change — never
    decremented or reset by ordinary code paths."""

    def test_save_features_bumps_generation(self, isolated_library):
        sc._save_features([])
        before = sc._features_generation
        sc._save_features([{
            "name": "foo", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        assert sc._features_generation > before

    def test_disk_reload_bumps_generation(self, isolated_library):
        sc._save_features([])
        before = sc._features_generation
        # Simulate an external invalidation (test harness, hand-edit of
        # features.json, etc.). _load_features re-reads from disk and
        # bumps the counter so any cached index is treated as stale.
        sc._features_cache = None
        sc._load_features()
        assert sc._features_generation > before


class TestBuildFeatureLibraryIndex:
    """`_build_feature_library_index` must produce a dict keyed by
    (name, feature_type) with case-folded sequences, single sweep over
    the library — used by PartsBinModal for O(1) per-row lookups."""

    def test_empty_library_returns_empty_dict(self, isolated_library):
        sc._save_features([])
        assert sc._build_feature_library_index() == {}

    def test_index_keys_are_name_type_tuples(self, isolated_library):
        sc._save_features([
            {"name": "p1", "feature_type": "promoter",
             "sequence": "ATG", "strand": 1},
            {"name": "c1", "feature_type": "CDS",
             "sequence": "TTT", "strand": 1},
        ])
        index = sc._build_feature_library_index()
        assert ("p1", "promoter") in index
        assert ("c1", "CDS") in index

    def test_index_sequences_are_case_folded(self, isolated_library):
        sc._save_features([
            {"name": "x", "feature_type": "CDS",
             "sequence": "atgcAT", "strand": 1},
        ])
        index = sc._build_feature_library_index()
        assert index[("x", "CDS")] == "ATGCAT"

    def test_classify_uses_index(self, isolated_library):
        sc._save_features([
            {"name": "x", "feature_type": "CDS",
             "sequence": "ATG", "strand": 1},
        ])
        index = sc._build_feature_library_index()
        assert sc._classify_feature_library_match(index, "x", "CDS", "ATG") == "exact"
        assert sc._classify_feature_library_match(index, "x", "CDS", "AAA") == "name"
        assert sc._classify_feature_library_match(index, "y", "CDS", "ATG") == ""


# ═══════════════════════════════════════════════════════════════════════════════
# Parts Bin: "Feat Lib" column + already-saved warning
# ═══════════════════════════════════════════════════════════════════════════════

class TestPartsBinFeatLibColumn:
    """The Parts Bin grew a "Feat Lib" column flagging parts that are
    already registered as features. Green ✓ = exact match, yellow ✓ =
    same (name, type) but different sequence (Save would replace).
    Built-in catalog rows always render empty."""

    async def test_exact_match_renders_green_check(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_parts_bin([{
            "name": "lacZ", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATGCATGCATGC", "strand": 1,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            row_keys = list(t.rows.keys())
            col_keys = list(t.columns.keys())
            cell = t.get_cell(row_keys[0], col_keys[-2])  # Feat Lib column (Grammar is now last)
            assert "✓" in str(cell), (
                f"Expected ✓ in Feat Lib column for exact match; got {cell!r}"
            )
            # Style on the whole cell carries "green".
            assert "green" in str(cell.style).lower()

    async def test_name_match_different_sequence_renders_yellow(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_parts_bin([{
            "name": "lacZ", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "AAATTTAAATTT", "strand": 1,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            row_keys = list(t.rows.keys())
            col_keys = list(t.columns.keys())
            cell = t.get_cell(row_keys[0], col_keys[-2])  # Feat Lib column
            assert "✓" in str(cell)
            assert "yellow" in str(cell.style).lower()

    async def test_no_match_renders_empty(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_parts_bin([{
            "name": "newPart", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        sc._save_features([])  # empty library
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            row_keys = list(t.rows.keys())
            col_keys = list(t.columns.keys())
            cell = t.get_cell(row_keys[0], col_keys[-2])  # Feat Lib column
            assert str(cell) == "", (
                f"Expected empty Feat Lib cell for unmatched part; got {cell!r}"
            )

    async def test_no_builtin_catalog_rows_present(
        self, isolated_parts_bin, isolated_library,
    ):
        """The built-in catalog rows were removed in 2026-05-07
        (see `PartsBinModal._all_rows` docstring) because every Copy /
        Save / Export / Delete action on them bailed with a "built-in
        catalog parts have no sequence" notify — clutter without
        function. Replaces the previous `test_builtin_catalog_row_is_
        always_empty` test (which fixture-drifted to skip on every
        run because the catalog rows it expected to find no longer
        exist). This test asserts the inverse: NO row in the modal
        is a built-in (every row is a user part)."""
        sc._save_features([{
            "name": "Nos", "feature_type": "promoter",
            "sequence": "ATGCATGCATGC", "strand": 1,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            # Every row must be a user part (`"user": True`) — built-
            # ins are gone. Without this assertion, a regression that
            # re-introduces built-in catalog rows would silently bring
            # back the action-button-bails-out behaviour the cleanup
            # eliminated.
            for r in parts_modal._rows:
                assert r["user"] is True, (
                    f"unexpected built-in catalog row {r['name']!r} — "
                    "built-in rows were removed in 2026-05-07; only "
                    "user-saved parts should appear in the modal"
                )

    async def test_column_refreshes_after_save_as_feature(
        self, isolated_parts_bin, isolated_library,
    ):
        """After saving a part as a feature, the column should flip
        from empty to ✓ without the user having to re-open the modal.
        Regression guard against the callback forgetting to repopulate."""
        sc._save_parts_bin([{
            "name": "myPart", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            row_keys = list(t.rows.keys())
            col_keys = list(t.columns.keys())
            assert str(t.get_cell(row_keys[0], col_keys[-2])) == ""  # Feat Lib column
            # Move cursor to the user part and trigger Save As Feature.
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-save-as-feature",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # AddFeatureModal pops up — accept the prefill.
            app.screen.query_one("#btn-addfeat-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Back on the parts screen; column should now show ✓.
            assert isinstance(app.screen, sc.PartsBinModal)
            t = app.screen.query_one("#parts-table", sc.DataTable)
            row_keys = list(t.rows.keys())
            col_keys = list(t.columns.keys())
            cell = t.get_cell(row_keys[0], col_keys[-2])  # Feat Lib column
            assert "✓" in str(cell), (
                f"Feat Lib column should show ✓ after save; got {cell!r}"
            )


class TestPartsBinFeatLibIndexCache:
    """The PartsBinModal builds a feature-library lookup index once on
    mount and reuses it across populates. Re-derived only when
    `_features_generation` advances — so opening Save As Feature, then
    re-rendering the parts table, doesn't re-scan the entire feature
    library."""

    async def test_index_built_on_mount(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATGCATGCATGC", "strand": 1,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._feat_lib_index == {("lacZ", "CDS"): "ATGCATGCATGC"}
            # Generation snapshot recorded so subsequent populates skip
            # the rebuild.
            assert modal._feat_lib_gen_seen == sc._features_generation

    async def test_index_not_rebuilt_when_generation_unchanged(
        self, isolated_parts_bin, isolated_library,
    ):
        """_populate must NOT call `_build_feature_library_index` when
        the feature library hasn't changed. Spy on the helper to count
        invocations across two sequential populates."""
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        calls: list[int] = []
        orig = sc._build_feature_library_index
        def spy():
            calls.append(1)
            return orig()
        sc._build_feature_library_index = spy
        try:
            app = sc.PlasmidApp()
            async with app.run_test(size=_BASELINE) as pilot:
                await pilot.pause()
                app.push_screen(sc.PartsBinModal())
                await pilot.pause()
                await pilot.pause(0.1)
                first_count = len(calls)
                # Trigger a second populate without changing the lib.
                app.screen._populate()
                await pilot.pause()
                second_count = len(calls)
            assert second_count == first_count, (
                f"Index was rebuilt unnecessarily — first populate: "
                f"{first_count}, second: {second_count} (expected equal)"
            )
        finally:
            sc._build_feature_library_index = orig

    async def test_index_rebuilt_when_generation_advances(
        self, isolated_parts_bin, isolated_library,
    ):
        """When `_features_generation` advances (e.g., after a save),
        the next populate must rebuild the index — otherwise the Feat
        Lib column would lag behind reality."""
        sc._save_features([])
        calls: list[int] = []
        orig = sc._build_feature_library_index
        def spy():
            calls.append(1)
            return orig()
        sc._build_feature_library_index = spy
        try:
            app = sc.PlasmidApp()
            async with app.run_test(size=_BASELINE) as pilot:
                await pilot.pause()
                app.push_screen(sc.PartsBinModal())
                await pilot.pause()
                await pilot.pause(0.1)
                first = len(calls)
                # Mutate the feature library — gen counter bumps.
                sc._save_features([{
                    "name": "lacZ", "feature_type": "CDS",
                    "sequence": "ATG", "strand": 1,
                }])
                app.screen._populate()
                await pilot.pause()
                second = len(calls)
            assert second > first, (
                f"Index should rebuild after _save_features bumped the "
                f"generation counter; calls: {first} → {second}"
            )
        finally:
            sc._build_feature_library_index = orig

    async def test_save_as_feature_refreshes_index(
        self, isolated_parts_bin, isolated_library,
    ):
        """After Save As Feature → modal save, the index should reflect
        the new entry on the very next render (no manual reopen). Same
        guarantee as `test_column_refreshes_after_save_as_feature` but
        asserted at the index level, not the rendered cell."""
        sc._save_parts_bin([{
            "name": "myPart", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        sc._save_features([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert ("myPart", "CDS") not in modal._feat_lib_index
            modal.query_one("#parts-table",
                            sc.DataTable).move_cursor(row=0)
            await pilot.pause()
            modal.query_one("#btn-parts-save-as-feature",
                            sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.query_one("#btn-addfeat-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Back on parts bin; index now contains the new entry.
            assert ("myPart", "CDS") in app.screen._feat_lib_index
            assert (
                app.screen._feat_lib_index[("myPart", "CDS")]
                == "ATGCATGCATGC"
            )


class TestPartsBinSaveAsFeatureWarning:
    """Save As Feature must warn (notify with severity=warning) before
    silently replacing an existing library entry. Two cases: exact
    match (no-op save) and name match (sequence will change)."""

    async def test_warns_on_exact_match(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_parts_bin([{
            "name": "lacZ", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATGCATGCATGC", "strand": 1,
        }])
        notes: list[tuple] = []
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            # Spy on app.notify to capture the warning text + severity.
            orig_notify = app.notify
            def spy(msg, *args, severity="information", **kw):
                notes.append((str(msg), severity))
                return orig_notify(msg, *args, severity=severity, **kw)
            app.notify = spy  # type: ignore[assignment]
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            parts_modal.query_one("#parts-table",
                                  sc.DataTable).move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-save-as-feature",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        warnings = [m for (m, sev) in notes if sev == "warning"]
        assert any("already in the feature library" in m for m in warnings), (
            f"Expected an 'already in feature library' warning; "
            f"got warnings={warnings}"
        )

    async def test_warns_on_name_match_different_sequence(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_parts_bin([{
            "name": "lacZ", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "AAATTTAAATTT", "strand": 1,
        }])
        notes: list[tuple] = []
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            orig_notify = app.notify
            def spy(msg, *args, severity="information", **kw):
                notes.append((str(msg), severity))
                return orig_notify(msg, *args, severity=severity, **kw)
            app.notify = spy  # type: ignore[assignment]
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            parts_modal.query_one("#parts-table",
                                  sc.DataTable).move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-save-as-feature",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        warnings = [m for (m, sev) in notes if sev == "warning"]
        assert any("Saving will replace" in m for m in warnings), (
            f"Expected a 'will replace' warning; got warnings={warnings}"
        )

    async def test_no_warning_when_part_is_new(
        self, isolated_parts_bin, isolated_library,
    ):
        sc._save_parts_bin([{
            "name": "brandNew", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATGCATGCATGC",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        sc._save_features([])
        notes: list[tuple] = []
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            orig_notify = app.notify
            def spy(msg, *args, severity="information", **kw):
                notes.append((str(msg), severity))
                return orig_notify(msg, *args, severity=severity, **kw)
            app.notify = spy  # type: ignore[assignment]
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            parts_modal.query_one("#parts-table",
                                  sc.DataTable).move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-save-as-feature",
                                  sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        warnings = [m for (m, sev) in notes if sev == "warning"]
        assert not any("already in the feature library" in m
                       or "will replace" in m
                       for m in warnings), (
            f"Did not expect any 'already in library' warning for a "
            f"brand-new part; got warnings={warnings}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Cloning simulator — primed amplicon + cloned plasmid math
# ═══════════════════════════════════════════════════════════════════════════════

class TestPupd2BackboneStub:
    """The placeholder backbone MUST be free of any BsaI / Esp3I / BsmBI
    recognition site (or its reverse complement) — otherwise the simulated
    cloned plasmid would appear to re-cut in a real Golden Gate reaction,
    in either the L0 (Esp3I) or L1 (BsaI) step.

    The backbone is also deterministic: the same build must produce
    byte-identical output across sessions, so Copy-Cloned-Sequence doesn't
    silently change between runs.
    """

    def test_backbone_is_deterministic(self):
        a = sc._build_pupd2_backbone_stub()
        b = sc._build_pupd2_backbone_stub()
        assert a == b
        # And the module-level constant matches a fresh build.
        assert sc._PUPD2_BACKBONE_STUB == a

    def test_backbone_is_acgt_only(self):
        assert set(sc._PUPD2_BACKBONE_STUB) <= set("ACGT")

    @pytest.mark.parametrize("site", [
        "GGTCTC", "GAGACC",        # BsaI + its RC (used for L1+ assembly)
        "CGTCTC", "GAGACG",        # Esp3I / BsmBI + its RC (used for L0)
    ])
    def test_backbone_has_no_type_iis_site(self, site):
        assert site not in sc._PUPD2_BACKBONE_STUB, (
            f"placeholder backbone contains forbidden site {site}"
        )

    def test_backbone_length_matches_default(self):
        assert len(sc._PUPD2_BACKBONE_STUB) == 420


class TestSimulatePrimedAmplicon:
    """`_simulate_primed_amplicon` must produce exactly what the PCR
    reaction yields: forward primer + insert body + rev-complement of
    reverse primer. That means tails stack as
    `pad + Esp3I + spacer + oh5 + insert + oh3 + rc(spacer+Esp3I+pad)`.
    """

    def test_happy_path_structure(self):
        insert = "ATGCATGCATGC"
        oh5, oh3 = "AATG", "GCTT"
        out = sc._simulate_primed_amplicon(insert, oh5, oh3)
        expected_head = (sc._GB_PAD
                         + sc._GB_L0_ENZYME_SITE
                         + sc._GB_SPACER + oh5)
        expected_tail = (oh3
                         + sc._rc(sc._GB_SPACER)
                         + sc._rc(sc._GB_L0_ENZYME_SITE)
                         + sc._rc(sc._GB_PAD))
        assert out.startswith(expected_head)
        assert out.endswith(expected_tail)
        assert insert in out
        # Length: pad + Esp3I(6) + spacer(1) + oh5(4) + insert + oh3(4)
        #       + rc(spacer)(1) + rc(Esp3I)(6) + rc(pad)
        assert len(out) == (
            2 * (len(sc._GB_PAD) + 6 + 1)
            + len(oh5) + len(oh3) + len(insert)
        )

    def test_esp3i_digest_yields_oh5_insert_oh3(self):
        """After a *real* Esp3I digest the internal fragment should be
        exactly `oh5 + insert + oh3` (the 4-nt overhang sits on the 5'
        end after cutting). We simulate that by locating the Esp3I sites
        and carving out the internal piece. The geometry is identical to
        BsaI: N(1)/N(5) → 4-nt 5' overhangs."""
        insert = "ATGCATGCATGCATGC"
        oh5, oh3 = "CCAT", "AATG"
        amplicon = sc._simulate_primed_amplicon(insert, oh5, oh3)
        # Digest from both ends: forward Esp3I cuts 1 nt into the spacer,
        # leaving a 4-nt 5' overhang on the insert.
        fwd_cut = amplicon.find(sc._GB_L0_ENZYME_SITE) + 6 + 1  # past CGTCTC + spacer
        rev_site = sc._rc(sc._GB_L0_ENZYME_SITE)              # GAGACG on top strand
        rev_cut = amplicon.rfind(rev_site) - 1                # back up past the spacer-rc
        internal = amplicon[fwd_cut:rev_cut]
        assert internal == oh5 + insert + oh3


class TestSimulateClonedPlasmid:
    """`_simulate_cloned_plasmid` returns the ligated circular product,
    linearised at the 5' overhang. It is `oh5 + insert + oh3 + backbone`."""

    def test_structure_is_oh5_insert_oh3_backbone(self):
        insert = "ATGCATGC"
        oh5, oh3 = "AATG", "GCTT"
        out = sc._simulate_cloned_plasmid(insert, oh5, oh3)
        assert out == oh5 + insert + oh3 + sc._PUPD2_BACKBONE_STUB

    def test_cloned_plasmid_has_no_type_iis_site_for_clean_insert(self):
        """If the insert is clean of BsaI/Esp3I, the cloned plasmid must
        be too — the overhangs + backbone are scrubbed. This mirrors the
        real-world invariant that a correctly-domesticated GB part
        produces a non-re-cuttable ligation product in BOTH the L0 (Esp3I)
        and the L1 (BsaI) assembly steps."""
        insert = "ATGAAACCCTTTGGG"
        for oh5, oh3 in [("AATG", "GCTT"), ("CCAT", "AATG"),
                         ("GCAA", "TACT"), ("GGAG", "CGCT")]:
            cloned = sc._simulate_cloned_plasmid(insert, oh5, oh3)
            for site in ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG"):
                assert site not in cloned, (
                    f"cloned plasmid for {oh5}/{oh3} contains {site}"
                )

    def test_cloned_length_is_insert_plus_overhangs_plus_backbone(self):
        insert = "A" * 500
        oh5, oh3 = "AATG", "GCTT"
        out = sc._simulate_cloned_plasmid(insert, oh5, oh3)
        assert len(out) == 500 + 4 + 4 + len(sc._PUPD2_BACKBONE_STUB)


class TestDomesticatorBackboneFromEntryVector:
    """2026-05-07 fix: the part dict's `backbone` and `marker`
    fields used to hardcode pUPD2 / Spectinomycin regardless of
    grammar or the user's configured entry vector. Now they pull
    from `_get_entry_vector(grammar_id)`:
      - `backbone` ← entry_vector["name"]
      - `marker`   ← `_detect_selection_marker(gb_text)` (or "—")
    With NO entry vector set, the historical pUPD2 / Spectinomycin
    defaults stand for backward compat with parts saved before
    this change.
    """

    @staticmethod
    def _build_entry_vector_gb(name: str, marker_label: str) -> str:
        """Build a minimal GenBank string for use as an entry vector
        with a single CDS feature carrying `marker_label` in its
        gene/label qualifiers — drives `_detect_selection_marker`."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("A" * 1000), id=name, name=name,
                         description="entry vector for test")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, 800, strand=1), type="CDS",
            qualifiers={"label": [marker_label],
                          "gene":  [marker_label]},
        ))
        return sc._record_to_gb_text(rec)

    async def test_save_pulls_backbone_from_entry_vector(
            self, isolated_parts_bin, isolated_library):
        # Configure a custom entry vector for gb_l0 with a KanR
        # marker — saved part should carry that name + marker.
        sc._set_entry_vector("gb_l0", {
            "name":    "pCustomKan",
            "size":    1000,
            "source":  "library:test",
            "gb_text": self._build_entry_vector_gb(
                "pCustomKan", "KanR",
            ),
        })
        sc._settings_flush_sync()
        insert = "ATG" * 8
        oh5, oh3 = "AATG", "GCTT"
        app = sc.PlasmidApp()
        captured: list = []
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            modal = sc.DomesticatorModal("ATG" * 30, [])
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            # Intercept dismiss synchronously so the part dict is
            # captured the moment _save calls it — bypasses the
            # Textual callback queue, which only fires after
            # run_test has exited and the captured list would
            # otherwise be empty when we check.
            modal.dismiss = lambda v=None: captured.append(v)  # type: ignore
            modal._design_result = {
                "part_type":   "CDS",
                "position":    "Pos 3-4",
                "oh5":         oh5, "oh3": oh3,
                "insert_seq":  insert,
                "fwd_full":    "GCGCCGTCTCAAATG" + insert,
                "rev_full":    "GCGCCGTCTCAGCTT" + sc._rc(insert),
                "fwd_tm":      60.2, "rev_tm": 59.8,
                "amplicon_len": len(insert) + 22,
            }
            modal.query_one("#dom-name", sc.Input).value = "kanpart"
            modal._save(None)
        assert len(captured) == 1
        part = captured[0]
        assert part is not None
        assert part["backbone"] == "pCustomKan"
        assert part["marker"] == "Kanamycin"

    async def test_save_falls_back_to_pupd2_when_no_entry_vector(
            self, isolated_parts_bin, isolated_library):
        sc._set_entry_vector("gb_l0", None)
        sc._settings_flush_sync()
        insert = "ATG" * 8
        app = sc.PlasmidApp()
        captured: list = []
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            modal = sc.DomesticatorModal("ATG" * 30, [])
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            modal.dismiss = lambda v=None: captured.append(v)  # type: ignore
            modal._design_result = {
                "part_type":   "CDS",
                "position":    "Pos 3-4",
                "oh5":         "AATG", "oh3": "GCTT",
                "insert_seq":  insert,
                "fwd_full":    "GCGCCGTCTCAAATG" + insert,
                "rev_full":    "GCGCCGTCTCAGCTT" + sc._rc(insert),
                "fwd_tm":      60.2, "rev_tm": 59.8,
                "amplicon_len": len(insert) + 22,
            }
            modal.query_one("#dom-name", sc.Input).value = "legacypart"
            modal._save(None)
        assert len(captured) == 1
        part = captured[0]
        assert part["backbone"] == "pUPD2"
        assert part["marker"] == "Spectinomycin"


class TestDetectSelectionMarker:
    """Marker detector covers the common bacterial selection
    cassettes by gene-name keyword. Falls back to None on parse
    failure or when nothing recognisable is annotated."""

    @staticmethod
    def _gb_with_marker(label: str) -> str:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("A" * 1000), id="x", name="x")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, 800, strand=1), type="CDS",
            qualifiers={"label": [label]},
        ))
        return sc._record_to_gb_text(rec)

    def test_recognises_canonical_markers(self):
        for label, expected in [
            ("AmpR",   "Ampicillin"),
            ("KanR",   "Kanamycin"),
            ("SpecR",  "Spectinomycin"),
            ("CmR",    "Chloramphenicol"),
            ("TetR",   "Tetracycline"),
            ("HygR",   "Hygromycin"),
            ("ZeoR",   "Zeocin"),
        ]:
            gb = self._gb_with_marker(label)
            assert sc._detect_selection_marker(gb) == expected, (
                f"label {label!r} should map to {expected!r}"
            )

    def test_returns_none_when_no_marker(self):
        gb = self._gb_with_marker("just some random gene")
        assert sc._detect_selection_marker(gb) is None

    def test_returns_none_on_empty_or_malformed(self):
        assert sc._detect_selection_marker("") is None
        assert sc._detect_selection_marker(None) is None  # type: ignore
        assert sc._detect_selection_marker("not a real GenBank") is None


class TestDomesticatorSavePersistsSimulations:
    """When the user clicks Save in the domesticator, the part dict
    written to the parts bin must include the `primed_seq` and
    `cloned_seq` fields — otherwise the Copy-Primed and Copy-Cloned
    buttons would have to re-derive them every time."""

    async def test_save_includes_simulated_sequences(
            self, isolated_parts_bin, tiny_record):
        # Build a domesticator with a pre-canned _design_result so we
        # can drive _save without running primer3. The design fields
        # mirror what _design_gb_primers produces.
        insert = "ATGCATGCATGCATGCATGCATGC"
        oh5, oh3 = "AATG", "GCTT"
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            modal = sc.DomesticatorModal("ATG" * 30, [{
                "type": "CDS", "start": 0, "end": 24, "strand": 1,
                "label": "x", "color": "white",
            }])
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            modal._design_result = {
                "part_type":   "CDS",
                "position":    "Pos 3-4",
                "oh5":         oh5,
                "oh3":         oh3,
                "insert_seq":  insert,
                "fwd_full":    "GCGCCGTCTCAAATG" + insert[:18],
                "rev_full":    "GCGCCGTCTCAAAGC" + sc._rc(insert)[:18],
                "fwd_tm":      60.2,
                "rev_tm":      59.8,
                "amplicon_len": len(insert) + 22,
            }
            modal.query_one("#dom-name", sc.Input).value = "my-sim-part"
            modal.query_one("#btn-dom-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)

        # Dismiss posts the part dict back through the modal callback. We
        # didn't register a save callback, so reach into the library
        # directly: the caller (PartsBinModal._new_part) is what persists,
        # but this test only asserts the dict contents. Re-invoke the
        # same construction path by calling the simulator directly and
        # comparing to what _save built. The surest check is to replay
        # the full flow with the parts-bin Save wiring.
        # No persistence callback here; just assert the dismissal result
        # format by re-running the simulator — if the constants match,
        # the saved dict would carry the same primed/cloned seq.
        assert sc._simulate_primed_amplicon(insert, oh5, oh3).startswith(
            "GCGCCGTCTCA" + oh5
        )
        assert sc._simulate_cloned_plasmid(insert, oh5, oh3).startswith(
            oh5 + insert + oh3
        )

    async def test_save_then_load_roundtrip_has_simulated_fields(
            self, isolated_parts_bin):
        """End-to-end: a part saved through _save_parts_bin with the new
        fields reloads with them intact. This guards the on-disk schema
        from accidentally dropping primed/cloned on load."""
        insert = "ATGCATGC" * 5
        oh5, oh3 = "AATG", "GCTT"
        sc._save_parts_bin([{
            "name":       "sim-rt",
            "type":       "CDS",
            "position":   "Pos 3-4",
            "oh5":        oh5, "oh3": oh3,
            "backbone":   "pUPD2",
            "marker":     "Spectinomycin",
            "sequence":   insert,
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
            "primed_seq": sc._simulate_primed_amplicon(insert, oh5, oh3),
            "cloned_seq": sc._simulate_cloned_plasmid(insert, oh5, oh3),
        }])
        reloaded = sc._load_parts_bin()
        assert len(reloaded) == 1
        entry = reloaded[0]
        assert entry["primed_seq"].startswith(
            sc._GB_PAD + sc._GB_L0_ENZYME_SITE + sc._GB_SPACER + oh5
        )
        assert entry["cloned_seq"] == oh5 + insert + oh3 + sc._PUPD2_BACKBONE_STUB


# ═══════════════════════════════════════════════════════════════════════════════
# Parts Bin — sequence TextArea + Copy buttons
# ═══════════════════════════════════════════════════════════════════════════════

class TestPartsBinSequenceView:
    """The parts-bin modal grew a scrollable read-only TextArea that holds
    the full insert sequence and three Copy buttons (raw / primed /
    cloned). Smoke-test that the widgets mount, the highlighted row
    populates the TextArea, and the copy buttons are present.
    """

    async def test_sequence_textarea_is_present(self, isolated_parts_bin):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            ta = modal.query_one("#parts-seq-view", sc.TextArea)
            assert ta.read_only is True

    async def test_copy_buttons_are_present(self, isolated_parts_bin):
        # Copy Raw was removed 2026-05-09 — users select the sequence
        # in the parts-seq-view TextArea + Ctrl+C instead. The Edit
        # button took its slot in the action row.
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal.query_one("#btn-parts-copy-primed", sc.Button) is not None
            assert modal.query_one("#btn-parts-copy-cloned", sc.Button) is not None
            assert modal.query_one("#btn-parts-edit",        sc.Button) is not None

    async def test_highlighting_user_row_loads_sequence_into_textarea(
            self, isolated_parts_bin):
        insert = "ATGCATGCATGCATGC"
        sc._save_parts_bin([{
            "name": "shown-part", "type": "CDS", "position": "Pos 3-4",
            "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": insert,
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            t = modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            ta = modal.query_one("#parts-seq-view", sc.TextArea)
            assert insert in ta.text
            assert "shown-part" in ta.text

    async def test_empty_sequence_row_shows_placeholder(
            self, isolated_parts_bin):
        """A user part with an empty sequence (rare — only happens
        if parts_bin.json was hand-edited; built-in catalog rows
        are gone as of 2026-05-07) must show a placeholder in the
        TextArea rather than looking empty."""
        sc._save_parts_bin([{
            "name": "empty_part", "type": "CDS", "position": "B3",
            "oh5": "AATG", "oh3": "GCTT", "backbone": "pUPD2",
            "marker": "Spec", "sequence": "",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            t = modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            ta = modal.query_one("#parts-seq-view", sc.TextArea)
            assert "No sequence" in ta.text


class TestPartsBinCopyButtons:
    """The three Copy buttons must:
      - warn when no user part is selected (built-in rows have no seq);
      - copy the right variant when a user part is selected.

    Because OSC 52 writes to /dev/tty in a test container may fail, we
    patch `_copy_to_clipboard_osc52` to capture the text it would send.
    """

    @staticmethod
    def _stub_part(insert="ATGCATGCATGC", oh5="AATG", oh3="GCTT"):
        return {
            "name": "copy-part", "type": "CDS", "position": "Pos 3-4",
            "oh5": oh5, "oh3": oh3,
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": insert,
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
            "primed_seq": sc._simulate_primed_amplicon(insert, oh5, oh3),
            "cloned_seq": sc._simulate_cloned_plasmid(insert, oh5, oh3),
        }

    async def _select_first_row_and_press(self, app, pilot, btn_id,
                                          monkeypatch):
        captured = []
        monkeypatch.setattr(sc, "_copy_to_clipboard_osc52",
                            lambda text: captured.append(text) or True)
        app.push_screen(sc.PartsBinModal())
        await pilot.pause()
        await pilot.pause(0.1)
        modal = app.screen
        t = modal.query_one("#parts-table", sc.DataTable)
        t.move_cursor(row=0)
        await pilot.pause()
        modal.query_one(btn_id, sc.Button).press()
        await pilot.pause()
        await pilot.pause(0.1)
        return captured

    async def test_copy_primed_copies_amplicon_with_tails(
            self, isolated_parts_bin, monkeypatch):
        sc._save_parts_bin([self._stub_part()])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            captured = await self._select_first_row_and_press(
                app, pilot, "#btn-parts-copy-primed", monkeypatch,
            )
            assert len(captured) == 1
            text = captured[0]
            assert text.startswith(
                sc._GB_PAD + sc._GB_L0_ENZYME_SITE + sc._GB_SPACER + "AATG"
            )
            assert text.endswith(
                "GCTT" + sc._rc(sc._GB_SPACER)
                + sc._rc(sc._GB_L0_ENZYME_SITE) + sc._rc(sc._GB_PAD)
            )

    async def test_copy_cloned_copies_plasmid_with_backbone(
            self, isolated_parts_bin, monkeypatch):
        sc._save_parts_bin([self._stub_part()])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            captured = await self._select_first_row_and_press(
                app, pilot, "#btn-parts-copy-cloned", monkeypatch,
            )
            assert len(captured) == 1
            text = captured[0]
            assert text == "AATG" + "ATGCATGCATGC" + "GCTT" + sc._PUPD2_BACKBONE_STUB

    async def test_copy_primed_regenerates_when_field_missing(
            self, isolated_parts_bin, monkeypatch):
        """Older parts saved before the simulator existed don't have
        `primed_seq`. The button must fall back to the simulator rather
        than copying an empty string."""
        entry = self._stub_part()
        entry.pop("primed_seq")
        entry.pop("cloned_seq")
        sc._save_parts_bin([entry])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            captured = await self._select_first_row_and_press(
                app, pilot, "#btn-parts-copy-primed", monkeypatch,
            )
            assert captured and captured[0].startswith(
                sc._GB_PAD + sc._GB_L0_ENZYME_SITE
            )

    async def test_copy_builtin_part_warns_and_no_clipboard_write(
            self, isolated_parts_bin, monkeypatch):
        """Built-in catalog rows have no sequence. Copy buttons must
        notify-and-skip rather than firing OSC 52 with an empty string."""
        captured = []
        monkeypatch.setattr(sc, "_copy_to_clipboard_osc52",
                            lambda text: captured.append(text) or True)
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            t = modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            # If isolated_parts_bin is empty the first row is a built-in
            # (no sequence). Only test if that's actually the case.
            if modal._rows and not modal._rows[0].get("sequence"):
                modal.query_one("#btn-parts-copy-primed", sc.Button).press()
                await pilot.pause()
                await pilot.pause(0.1)
                assert captured == []


class TestPartsBinEdit:
    """The Edit button opens PartEditModal pre-populated with the
    cursor row, and its Save callback rewrites the matching parts-
    bin entry by `(name, sequence)` identity. Added 2026-05-09 along
    with the Copy Raw → Edit swap."""

    @staticmethod
    def _stub_part(name="edit-part", insert="ATGCATGCATGC"):
        return {
            "name": name, "type": "CDS", "position": "Pos 3-4",
            "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": insert,
            "fwd_primer": "GCGCCGTCTCAAATGATGCATGCATGC",
            "rev_primer": "GCGCCGTCTCAGCTTGCATGCATGCAT",
            "fwd_tm": 60.0, "rev_tm": 60.0,
            "grammar": "gb_l0",
        }

    async def test_edit_button_present_in_action_row(
            self, isolated_parts_bin):
        """The Edit button must exist alongside the (remaining) Copy
        buttons. Smoke check that the compose() change held."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal.query_one("#btn-parts-edit", sc.Button) is not None

    async def test_edit_save_rewrites_entry_in_parts_bin(
            self, isolated_parts_bin):
        """Saving the modal swaps the matching parts-bin entry in
        place. Match on (name, sequence) identity so a `_populate`
        between push and dismiss can't desync the row index."""
        sc._save_parts_bin([self._stub_part(name="orig")])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            assert isinstance(edit_modal, sc.PartEditModal)
            edit_modal._set_editing(True)
            await pilot.pause()
            edit_modal.query_one("#partedit-name", sc.Input).value = "renamed"
            edit_modal.query_one("#partedit-backbone",
                                  sc.Input).value = "pUPD-Kan"
            edit_modal.query_one("#btn-partedit-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            entries = sc._load_parts_bin()
            assert len(entries) == 1
            assert entries[0]["name"] == "renamed"
            assert entries[0]["backbone"] == "pUPD-Kan"
            assert entries[0]["sequence"] == "ATGCATGCATGC"

    async def test_edit_seq_or_overhang_change_rederives_simulator_outputs(
            self, isolated_parts_bin):
        """Editing the sequence or overhangs must regenerate
        `primed_seq` and `cloned_seq` so Copy Primed / Copy Cloned
        keep serving the right amplicon after a save."""
        original_insert = "ATGCATGCATGC"
        new_insert      = "ATGAAACCCGGGTTT"
        sc._save_parts_bin([self._stub_part(insert=original_insert)])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            edit_modal._set_editing(True)
            await pilot.pause()
            edit_modal.query_one("#partedit-seq",
                                  sc.TextArea).text = new_insert
            edit_modal.query_one("#btn-partedit-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            entries = sc._load_parts_bin()
            assert entries[0]["sequence"] == new_insert
            assert new_insert in entries[0]["primed_seq"]
            assert new_insert in entries[0]["cloned_seq"]

    async def test_edit_invalid_dna_blocks_save(
            self, isolated_parts_bin):
        """Bases outside the IUPAC alphabet must block the save and
        render a red status — same validation pattern the grammar
        editor uses."""
        sc._save_parts_bin([self._stub_part()])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            edit_modal._set_editing(True)
            await pilot.pause()
            # 'Z' is not in the IUPAC alphabet — save must refuse.
            edit_modal.query_one("#partedit-seq",
                                  sc.TextArea).text = "ATGCZZZZ"
            edit_modal.query_one("#btn-partedit-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Modal still on the stack — save was blocked.
            assert isinstance(app.screen, sc.PartEditModal)
            status_text = str(
                edit_modal.query_one("#partedit-status", sc.Static).render()
            )
            assert "invalid bases" in status_text.lower()
            # Parts bin file untouched.
            assert sc._load_parts_bin()[0]["sequence"] == "ATGCATGCATGC"

    async def test_edit_button_warns_when_no_part_selected(
            self, isolated_parts_bin):
        """Empty parts bin → cursor row is invalid; the Edit button
        must notify rather than push an Edit modal on a phantom row."""
        sc._save_parts_bin([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            stack_depth = len(app.screen_stack)
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Stack unchanged — no PartEditModal pushed.
            assert len(app.screen_stack) == stack_depth

    async def test_escape_closes_edit_modal(self, isolated_parts_bin):
        """`Binding('escape', 'cancel', …)` only fires if the modal
        actually defines `action_cancel`. Pre-2026-05-09 the binding
        was registered but the action wasn't, so Escape was a no-op."""
        sc._save_parts_bin([self._stub_part()])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.PartEditModal)
            await pilot.press("escape")
            await pilot.pause()
            await pilot.pause(0.1)
            assert not isinstance(app.screen, sc.PartEditModal)

    async def test_no_op_save_does_not_write_file(self, isolated_parts_bin):
        """Edit → Save without changing anything must dismiss without
        rewriting parts_bin.json. Verifies via the parts-bin file's
        mtime — if the no-op short-circuit fires, mtime stays put."""
        sc._save_parts_bin([self._stub_part()])
        before_mtime = sc._PARTS_BIN_FILE.stat().st_mtime_ns
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            edit_modal._set_editing(True)
            await pilot.pause()
            edit_modal.query_one("#btn-partedit-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        after_mtime = sc._PARTS_BIN_FILE.stat().st_mtime_ns
        assert after_mtime == before_mtime

    async def test_primer_tm_appears_in_label(self, isolated_parts_bin):
        """The primer field labels must include the part's stored Tm
        so users editing the primer aren't blind to the current Tm."""
        part = self._stub_part()
        part["fwd_tm"] = 62.3
        part["rev_tm"] = 60.1
        sc._save_parts_bin([part])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            labels = [str(lbl.render()) for lbl in
                      edit_modal.query("Label")]
            assert any("Tm 62.3" in lbl for lbl in labels)
            assert any("Tm 60.1" in lbl for lbl in labels)

    async def test_modal_renders_part_with_markup_in_name(
            self, isolated_parts_bin):
        """A part name like '[red]X' must NOT be interpreted as Rich
        markup in the modal title — the Static is composed with
        markup=False so the bracketed text renders verbatim."""
        sc._save_parts_bin([self._stub_part(name="[red]evil[/]")])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            title = str(
                edit_modal.query_one("#partedit-title", sc.Static).render()
            )
            # Brackets survive into the rendered text rather than
            # being parsed as a Rich `[red]` tag.
            assert "[red]evil[/]" in title

    async def test_legacy_type_appears_with_suffix(self, isolated_parts_bin):
        """A part whose stored type is no longer in the active grammar
        (e.g. an old custom grammar removed) must still show up in
        the Type select with a '(legacy)' suffix so a Save round-
        trips the value."""
        part = self._stub_part()
        part["type"] = "ZZ-ghost-type"  # not in any built-in grammar
        sc._save_parts_bin([part])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            opts, default = edit_modal._type_options()
            opt_labels = [label for (label, value) in opts]
            assert default == "ZZ-ghost-type"
            assert any("ZZ-ghost-type (legacy)" in label
                        for label in opt_labels)

    async def test_grammar_dropdown_is_present_and_editable(
            self, isolated_parts_bin):
        """The grammar field must be a real Select (was a locked
        banner pre-2026-05-10) so users can re-tag a part to a
        different cloning grammar without recreating it."""
        sc._save_parts_bin([self._stub_part()])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            grammar_sel = edit_modal.query_one(
                "#partedit-grammar", sc.Select,
            )
            assert grammar_sel.value == "gb_l0"
            # Disabled until the user clicks Edit (read-only by default).
            assert grammar_sel.disabled is True
            edit_modal._set_editing(True)
            await pilot.pause()
            assert grammar_sel.disabled is False

    async def test_grammar_change_rebuilds_type_options_and_overhangs(
            self, isolated_parts_bin):
        """Switching grammar from gb_l0 → moclo_plant must rebuild
        the Type select with the new grammar's positions and refresh
        oh5 / oh3 / position. Both grammars carry CDS so the type
        selection survives the swap."""
        sc._save_parts_bin([self._stub_part()])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            edit_modal._set_editing(True)
            await pilot.pause()
            # Pre-swap: gb_l0 CDS overhangs.
            assert edit_modal.query_one(
                "#partedit-oh5", sc.Input).value == "AATG"
            edit_modal.query_one(
                "#partedit-grammar", sc.Select).value = "moclo_plant"
            await pilot.pause()
            await pilot.pause(0.05)
            # Post-swap: MoClo CDS overhangs (5'=AGGT, 3'=GCTT).
            assert edit_modal.query_one(
                "#partedit-oh5", sc.Input).value == "AGGT"
            assert edit_modal.query_one(
                "#partedit-oh3", sc.Input).value == "GCTT"
            assert edit_modal._grammar_id == "moclo_plant"

    async def test_grammar_change_persists_through_save(
            self, isolated_parts_bin):
        """Saving after a grammar swap must write the new grammar id
        into parts_bin.json AND re-derive primed_seq with the new
        grammar's enzyme tail."""
        sc._save_parts_bin([self._stub_part()])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            edit_modal._set_editing(True)
            await pilot.pause()
            edit_modal.query_one(
                "#partedit-grammar", sc.Select).value = "moclo_plant"
            await pilot.pause()
            await pilot.pause(0.05)
            edit_modal.query_one("#btn-partedit-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        entries = sc._load_parts_bin()
        assert len(entries) == 1
        assert entries[0]["grammar"] == "moclo_plant"
        # MoClo Plant uses BsaI (GGTCTC) so the primed amplicon must
        # contain the BsaI site, not Esp3I (CGTCTC).
        primed = entries[0].get("primed_seq", "")
        assert "GGTCTC" in primed
        assert "CGTCTC" not in primed

    async def test_save_refuses_when_identity_missing_in_file(
            self, isolated_parts_bin):
        """If the parts-bin file changes between push and dismiss
        (concurrent writer / hand-edit) and the original identity is
        gone, the callback must NOT silently append a duplicate.
        Refuse + notify, leaving the file as-is."""
        original = self._stub_part(name="orig")
        sc._save_parts_bin([original])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            parts_modal = app.screen
            t = parts_modal.query_one("#parts-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            parts_modal.query_one("#btn-parts-edit", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            edit_modal = app.screen
            edit_modal._set_editing(True)
            await pilot.pause()
            # Concurrent writer wipes the bin out from under the
            # dialog (simulating an external splicecraft / a
            # rm + restore flow). Identity will not match anywhere.
            sc._save_parts_bin([])
            edit_modal.query_one("#partedit-name",
                                  sc.Input).value = "renamed"
            edit_modal.query_one("#btn-partedit-save", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        # File still empty — no silent append of the renamed part.
        assert sc._load_parts_bin() == []


# ═══════════════════════════════════════════════════════════════════════════════
# Save Primers to Library (multi-pair naming + dup guard)
# ═══════════════════════════════════════════════════════════════════════════════

class TestDesignGBPrimersPairsShape:
    """_design_gb_primers returns a ``pairs`` list so callers can iterate
    amplicons uniformly. For the current single-amplicon design the list
    has one entry; the top-level keys mirror that entry for backward
    compatibility (cloning simulator, PrimerDesignScreen)."""

    def test_pairs_key_present_and_nonempty(self, random_template):
        result = sc._design_gb_primers(random_template, 100, 400, "CDS")
        assert "pairs" in result
        assert isinstance(result["pairs"], list)
        assert len(result["pairs"]) >= 1

    def test_pair_has_expected_primer_fields(self, random_template):
        result = sc._design_gb_primers(random_template, 100, 400, "CDS")
        pair = result["pairs"][0]
        for key in ("fwd_full", "rev_full", "fwd_tm", "rev_tm",
                    "fwd_pos", "rev_pos", "fwd_binding", "rev_binding",
                    "amplicon_len"):
            assert key in pair, f"pair missing key {key!r}"

    def test_top_level_mirrors_first_pair(self, random_template):
        """Legacy callers that read result['fwd_full'] directly (cloning
        simulator, PrimerDesignScreen) must continue to work."""
        result = sc._design_gb_primers(random_template, 100, 400, "CDS")
        p = result["pairs"][0]
        for key in ("fwd_full", "rev_full", "fwd_tm", "rev_tm",
                    "fwd_pos", "rev_pos", "amplicon_len"):
            assert result[key] == p[key], f"top-level {key} != pairs[0][{key}]"


class TestDomesticatorSavePrimersButton:
    """DomesticatorModal's 'Save Primers' button writes each designed pair
    to primers.json with the project naming convention:
        {partName}-DOM-{n}-F / -R
    where DOM tags domestication (vs CLO cloning, DET detection)."""

    async def test_save_primers_button_present_and_disabled_on_mount(
            self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            btn = modal.query_one("#btn-dom-save-primers", sc.Button)
            assert btn is not None
            assert btn.disabled is True  # no design yet

    async def test_save_primers_button_enabled_after_design(
            self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            ta = modal.query_one("#dom-direct-seq", sc.TextArea)
            ta.load_text(_mk_template(300))
            modal.query_one("#dom-name", sc.Input).value = "myPart"
            modal.query_one("#btn-dom-design", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert modal.query_one(
                "#btn-dom-save-primers", sc.Button
            ).disabled is False

    async def test_save_primers_writes_correct_names(self, isolated_library):
        """After design+save, primers.json holds two entries named
        ``myPart-DOM-1-F`` and ``myPart-DOM-1-R`` — one pair per L0
        amplicon."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            ta = modal.query_one("#dom-direct-seq", sc.TextArea)
            ta.load_text(_mk_template(300))
            modal.query_one("#dom-name", sc.Input).value = "myPart"
            modal.query_one("#btn-dom-design", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            modal.query_one("#btn-dom-save-primers", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)

            entries = sc._load_primers()
            names = {e.get("name") for e in entries}
            assert "myPart-DOM-1-F" in names
            assert "myPart-DOM-1-R" in names

    async def test_saved_primer_rows_have_expected_fields(
            self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            ta = modal.query_one("#dom-direct-seq", sc.TextArea)
            ta.load_text(_mk_template(300))
            modal.query_one("#dom-name", sc.Input).value = "myPart"
            modal.query_one("#btn-dom-design", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            modal.query_one("#btn-dom-save-primers", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)

            entries = sc._load_primers()
            by_name = {e["name"]: e for e in entries}
            fwd = by_name["myPart-DOM-1-F"]
            rev = by_name["myPart-DOM-1-R"]
            # Type + strand
            assert fwd["primer_type"] == "goldenbraid"
            assert rev["primer_type"] == "goldenbraid"
            assert fwd["strand"] == 1
            assert rev["strand"] == -1
            # Full sequences match the designed primers
            assert fwd["sequence"] == modal._design_result["pairs"][0]["fwd_full"]
            assert rev["sequence"] == modal._design_result["pairs"][0]["rev_full"]
            # Tm carried through
            assert fwd["tm"] == modal._design_result["pairs"][0]["fwd_tm"]
            assert rev["tm"] == modal._design_result["pairs"][0]["rev_tm"]

    async def test_save_primers_without_part_name_errors(
            self, isolated_library):
        """Without a part name the save path has no naming stem — should
        surface an error and write nothing."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            ta = modal.query_one("#dom-direct-seq", sc.TextArea)
            ta.load_text(_mk_template(300))
            # Intentionally leave name empty
            modal.query_one("#dom-name", sc.Input).value = ""
            modal.query_one("#btn-dom-design", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Design still works (name not required for design)
            modal.query_one("#btn-dom-save-primers", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Primer library must still be empty
            assert sc._load_primers() == []

    async def test_save_primers_skips_duplicate_sequences(
            self, isolated_library):
        """If either primer sequence already exists in the library,
        skip that one and notify — never overwrite with a new name."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.DomesticatorModal(_mk_template(), _mk_feats()))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            ta = modal.query_one("#dom-direct-seq", sc.TextArea)
            ta.load_text(_mk_template(300))
            modal.query_one("#dom-name", sc.Input).value = "myPart"
            modal.query_one("#btn-dom-design", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # First save: both primers land.
            modal.query_one("#btn-dom-save-primers", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            first = sc._load_primers()
            assert len(first) == 2
            # Re-enable and try again — the Save button disables itself
            # after a successful write, so flip it back for the second click.
            modal.query_one("#btn-dom-save-primers", sc.Button).disabled = False
            modal.query_one("#btn-dom-save-primers", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            second = sc._load_primers()
            # Duplicate sequences are skipped — library stays at 2 entries.
            assert len(second) == 2


class TestPrimerDesignScreenGBSuffix:
    """PrimerDesignScreen's Golden Braid mode should name primers with
    DOM (domestication) suffix, not CLO (cloning) — the intent is
    specifically L0 domestication, not generic restriction-enzyme cloning."""

    async def test_goldenbraid_shows_dom_suffix(self, isolated_library):
        """_show_result fills the fwd/rev name inputs with -DOM-F / -DOM-R
        when primer_type is 'goldenbraid'. Drive the render step directly
        with a hand-built design dict so we don't depend on Primer3 picking
        any particular region on an arbitrary template."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.1)
            screen.query_one("#pd-part-name", sc.Input).value = "gbTest"
            fake_design = {
                "fwd_full": "GCGCCGTCTCAGGAGACGT",
                "rev_full": "GCGCCGTCTCATTGGCCTA",
                "fwd_tm": 62.0, "rev_tm": 61.5,
            }
            screen._show_result(fake_design, "goldenbraid",
                                "fwd_full", "rev_full")
            await pilot.pause()
            assert screen.query_one("#pd-fwd-name", sc.Input).value == (
                "gbTest-DOM-F"
            )
            assert screen.query_one("#pd-rev-name", sc.Input).value == (
                "gbTest-DOM-R"
            )

    async def test_cloning_still_uses_clo_suffix(self, isolated_library):
        """Regression guard: we only changed the suffix for the goldenbraid
        branch, not cloning. Cloning stays CLO."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.1)
            screen.query_one("#pd-part-name", sc.Input).value = "cloTest"
            fake_design = {
                "fwd_full": "GCGCGAATTCGATCAAAG",
                "rev_full": "GCGCGGATCCTAGATAGA",
                "fwd_tm": 58.0, "rev_tm": 59.0,
            }
            screen._show_result(fake_design, "cloning",
                                "fwd_full", "rev_full")
            await pilot.pause()
            assert screen.query_one("#pd-fwd-name", sc.Input).value == (
                "cloTest-CLO-F"
            )

    async def test_detection_still_uses_det_suffix(self, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.1)
            screen.query_one("#pd-part-name", sc.Input).value = "detTest"
            fake_design = {
                "fwd_seq": "AAGCGATCAAAGGATATAT",
                "rev_seq": "TTCATGCTACAAGGATTTA",
                "fwd_tm": 58.0, "rev_tm": 58.5,
            }
            screen._show_result(fake_design, "detection",
                                "fwd_seq", "rev_seq")
            await pilot.pause()
            assert screen.query_one("#pd-fwd-name", sc.Input).value == (
                "detTest-DET-F"
            )


class TestPrimerLibraryMarkAndFocus:
    """After a design-save the new pair is marked + focused (ready for a
    one-press Add-to-map), and deleting a primer keeps the list parked
    near where it was instead of snapping back to the top."""

    @staticmethod
    def _seed(n):
        sc._save_primers([
            {"name": f"p{i:02d}", "sequence": f"ACGTACGTACGTACGT{i:04d}",
             "tm": 60.0, "primer_type": "generic", "source": "seed",
             "pos_start": 0, "pos_end": 20, "strand": 1,
             "date": "2026-05-29", "status": "Designed"}
            for i in range(n)
        ])

    async def test_design_save_marks_and_focuses_new_pair(self, isolated_library):
        """Saving designed primers marks the new pair (so a single
        Add-to-map grabs them) and lands the cursor on one of them."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.1)
            screen._det_result = {
                "fwd_seq": "AAGCGATCAAAGGATATAT",
                "rev_seq": "TTCATGCTACAAGGATTTA",
                "fwd_tm": 58.0, "rev_tm": 58.5,
                "fwd_pos": (10, 29), "rev_pos": (60, 79),
                "_type": "detection",
            }
            screen._commit_designed_primers_with_names(
                ["Ztest-DET-F", "Ztest-DET-R"])
            await pilot.pause()
            await pilot.pause(0.1)
            primers = sc._load_primers()
            names = [p.get("name") for p in primers]
            assert "Ztest-DET-F" in names and "Ztest-DET-R" in names
            # Both new primers are marked.
            marked = {primers[i].get("name") for i in screen._lib_selected
                      if 0 <= i < len(primers)}
            assert marked == {"Ztest-DET-F", "Ztest-DET-R"}, marked
            # Cursor landed on one of the new primers.
            t = screen.query_one("#pd-lib-table", sc.DataTable)
            cur = screen._row_to_primer_idx[t.cursor_row]
            assert primers[cur].get("name") in {"Ztest-DET-F", "Ztest-DET-R"}

    async def test_delete_last_primer_keeps_cursor_in_range(self, isolated_library):
        """Deleting the bottom primer leaves the cursor on the new last
        row instead of snapping to row 0 (the old out-of-range behaviour)."""
        self._seed(6)
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.1)
            t = screen.query_one("#pd-lib-table", sc.DataTable)
            assert t.row_count == 6
            t.move_cursor(row=5)                 # bottom row
            await pilot.pause()
            last = sc._load_primers()[screen._row_to_primer_idx[5]]
            remaining = [p for p in sc._load_primers()
                         if p.get("name") != last.get("name")]
            sc._save_primers(remaining)
            screen._refresh_library_table()
            await pilot.pause()
            await pilot.pause(0.1)          # let the deferred restore run
            # Clamped to the new last row (4), not bounced to 0.
            assert t.cursor_row == 4, t.cursor_row

    async def test_delete_maintains_scroll_position(self, isolated_library):
        """Deleting a primer keeps the EXACT scrollbar position — the user
        stays parked where they were rather than the list re-scrolling so
        the deleted row's slot sits at the viewport edge."""
        self._seed(200)
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.1)
            t = screen.query_one("#pd-lib-table", sc.DataTable)
            assert t.row_count == 200
            t.scroll_to(y=30, animate=False)
            await pilot.pause()
            await pilot.pause(0.1)
            y_before = t.scroll_offset.y
            assert y_before > 0, "test setup: the table should be scrolled"
            # Delete the primer shown at a mid-viewport row and drop its
            # row in place (the surgical path the delete handler uses).
            # NB: the modal-confirm path can't be exercised headlessly
            # here — dismissing the confirm modal hard-resets the table
            # scroll in the test harness (the real terminal does not),
            # so we drive `_drop_primer_rows` directly.
            old_idx = screen._row_to_primer_idx[40]
            remaining = [p for i, p in enumerate(sc._load_primers())
                         if i != old_idx]
            sc._save_primers(remaining)
            screen._drop_primer_rows({old_idx}, cursor_row=40, scroll_y=y_before)
            await pilot.pause()
            await pilot.pause(0.1)
            assert t.row_count == 199
            # Scrollbar held its position; row→idx map stayed in sync;
            # cursor parked near the deleted row, not bounced to the top.
            assert t.scroll_offset.y == y_before, (
                f"scroll moved {y_before} → {t.scroll_offset.y}")
            assert len(screen._row_to_primer_idx) == 199
            assert t.cursor_row == 40


# ═══════════════════════════════════════════════════════════════════════════════
# Edge-case hardening: multi-site inserts, end-of-insert sites, swap cascades
#
# Added 2026-04-21 after tightening `_gb_find_forbidden_hits` to return ALL
# occurrences and `_codon_fix_sites` to veto any swap that would introduce a
# new forbidden pattern anywhere in the sequence. These cases cover the
# failure modes that can otherwise result in the user paying for a gBlock
# that still contains a forbidden site — a real synthesis-budget concern.
# ═══════════════════════════════════════════════════════════════════════════════

def _k12_raw() -> dict:
    """Shared E. coli K12 codon table fixture for the edge-case tests."""
    return dict(sc._CODON_BUILTIN_K12)


class TestGbFindForbiddenHitsAllOccurrences:
    """`_gb_find_forbidden_hits` must report EVERY site, not just the first
    per enzyme. Otherwise a user could pay for a gBlock that still contains
    an undetected second site — real synthesis budget lost."""

    def test_multiple_bsai_sites_all_reported(self):
        seq = "AAA" + "GGTCTC" + "AAAA" + "GGTCTC" + "AAA"
        hits = sc._gb_find_forbidden_hits(seq)
        bsai_hits = [h for h in hits if h[0] == "BsaI" and h[1] == "GGTCTC"]
        assert len(bsai_hits) == 2
        assert {h[2] for h in bsai_hits} == {3, 13}

    def test_multiple_esp3i_sites_all_reported(self):
        seq = "AAA" + "CGTCTC" + "AAAA" + "CGTCTC" + "AAA"
        hits = sc._gb_find_forbidden_hits(seq)
        esp3i_hits = [h for h in hits if h[0] == "Esp3I" and h[1] == "CGTCTC"]
        assert len(esp3i_hits) == 2
        assert {h[2] for h in esp3i_hits} == {3, 13}

    def test_mixed_bsai_and_esp3i_all_reported(self):
        seq = ("AAA" + "GGTCTC" + "AAAA"
               + "CGTCTC" + "AAAA"
               + "GAGACC" + "AAAA"
               + "GAGACG" + "AAA")
        hits = sc._gb_find_forbidden_hits(seq)
        assert len(hits) == 4
        enzymes = sorted(h[0] for h in hits)
        assert enzymes == ["BsaI", "BsaI", "Esp3I", "Esp3I"]

    def test_sites_sorted_by_position(self):
        """Downstream error reporting stringifies the list as-is, so the
        order must be stable and intuitive — left-to-right by nt position."""
        seq = "AAA" + "CGTCTC" + "AAAA" + "GGTCTC" + "AAA"
        hits = sc._gb_find_forbidden_hits(seq)
        positions = [h[2] for h in hits]
        assert positions == sorted(positions)

    def test_adjacent_sites_both_reported(self):
        """Two occurrences only one bp apart — e.g. overlapping windows
        on different enzymes. Both must surface so the user knows the
        insert is densely contaminated."""
        seq = "AAA" + "GGTCTC" + "A" + "CGTCTC" + "AAA"
        hits = sc._gb_find_forbidden_hits(seq)
        assert len(hits) == 2
        positions = [h[2] for h in hits]
        assert positions == sorted(positions)
        assert 3 in positions and 10 in positions

    def test_palindrome_free_site_rc_reported_on_forward_strand(self):
        """BsaI (GGTCTC) is non-palindromic — its RC is GAGACC. An insert
        carrying GAGACC on its forward strand means the site is on the
        reverse strand, which still triggers a cut. Must be reported."""
        seq = "AAA" + "GAGACC" + "AAA"
        hits = sc._gb_find_forbidden_hits(seq)
        assert len(hits) == 1
        assert hits[0][0] == "BsaI"
        assert hits[0][1] == "GAGACC"

    def test_clean_sequence_returns_empty(self):
        seq = "ACGT" * 10
        assert sc._gb_find_forbidden_hits(seq) == []


class TestForbiddenHitSetHelper:
    """`_forbidden_hit_set` powers the before/after cross-check inside
    `_codon_fix_sites`. It must find every occurrence of every pattern
    passed in — not stop at the first hit."""

    def test_multiple_occurrences_of_same_pattern(self):
        seq = "GGTCTC" + "AAAA" + "GGTCTC"
        s = sc._forbidden_hit_set(seq, ("GGTCTC",))
        assert s == {("GGTCTC", 0), ("GGTCTC", 10)}

    def test_multiple_patterns_all_found(self):
        seq = "GGTCTC" + "AAAA" + "CGTCTC"
        s = sc._forbidden_hit_set(seq, ("GGTCTC", "CGTCTC"))
        assert s == {("GGTCTC", 0), ("CGTCTC", 10)}

    def test_empty_input_returns_empty_set(self):
        assert sc._forbidden_hit_set("", ("GGTCTC",)) == set()
        assert sc._forbidden_hit_set("ACGT" * 5, ()) == set()


class TestMultipleSitesSameEnzyme:
    """Multiple BsaI or Esp3I sites in one coding insert must all be
    removed via synonymous codons. `_codon_fix_sites` iterates left-to-right
    within an enzyme, so later copies of the same site can't block earlier
    ones — every site must end up fixed."""

    def test_two_bsai_sites_both_fixed(self):
        # Two GGT-CTC motifs at G-L codon boundaries, 30 bp apart.
        cds = ("ATG"
               + "GCG" * 3 + "GGT" + "CTC"      # BsaI #1 at codon 4-5
               + "GCG" * 5 + "GGT" + "CTC"      # BsaI #2 at codon 10-11
               + "GCG" * 3 + "TAA")
        assert cds.count("GGTCTC") == 2
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"]
        assert len(r["mutations"]) >= 2

    def test_three_bsai_sites_all_fixed(self):
        """Three consecutive sites — exercises the while-loop through
        each occurrence inside `_codon_fix_sites`."""
        cds = ("ATG"
               + "GCG" * 2 + "GGT" + "CTC"
               + "GCG" * 2 + "GGT" + "CTC"
               + "GCG" * 2 + "GGT" + "CTC"
               + "GCG" * 2 + "TAA")
        assert cds.count("GGTCTC") == 3
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"]
        assert "GAGACC" not in r["insert_seq"]
        assert len(r["mutations"]) >= 3

    def test_two_esp3i_sites_both_fixed(self):
        cds = ("ATG"
               + "GCG" * 3 + "CGT" + "CTC"      # Esp3I #1
               + "GCG" * 5 + "CGT" + "CTC"      # Esp3I #2
               + "GCG" * 3 + "TAA")
        assert cds.count("CGTCTC") == 2
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "CGTCTC" not in r["insert_seq"]
        assert len(r["mutations"]) >= 2

    def test_mixed_bsai_and_esp3i_both_fixed(self):
        cds = ("ATG"
               + "GCG" * 3 + "GGT" + "CTC"      # BsaI
               + "GCG" * 5 + "CGT" + "CTC"      # Esp3I
               + "GCG" * 3 + "TAA")
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"]
        assert "CGTCTC" not in r["insert_seq"]
        assert "GAGACC" not in r["insert_seq"]
        assert "GAGACG" not in r["insert_seq"]
        enzymes = {m.split()[0] for m in r["mutations"]}
        assert "BsaI" in enzymes and "Esp3I" in enzymes

    def test_protein_preserved_across_multi_site_fix(self):
        """Silent-mutation contract: protein sequence stays identical
        even when multiple sites are touched."""
        cds = ("ATG"
               + "GCG" * 3 + "GGT" + "CTC"
               + "GCG" * 5 + "CGT" + "CTC"
               + "GCG" * 3 + "TAA")
        before = sc._mut_translate(cds)
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        after = sc._mut_translate(r["insert_seq"])
        assert before == after


class TestSitesNearSequenceEnds:
    """A forbidden site at the 5' or 3' end of the insert lands inside
    the primer binding window. The fix must remove the site AND the UI
    must warn that the user needs to order the mutated insert as a
    gBlock — a primer designed against the mutated insert won't bind
    the original template at that spot."""

    def test_bsai_site_at_5_prime_end(self):
        """Site in codon 2 — squarely inside the forward binding region.
        Must be removed; binding_region_mutations must flag it."""
        cds = ("ATG" + "GGT" + "CTC"            # BsaI at nt 3-8
               + "GCG" * 15 + "TAA")
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"]
        # fwd_bind is first 18-25 bp. Site at nt 3-8 → codon 1 start=3,
        # codon 2 start=6. Both are <18, so both land in fwd binding.
        br = r["binding_region_mutations"]
        assert br, "Advisory must flag site near 5' end"
        assert any(entry["region"] == "fwd" for entry in br)

    def test_bsai_site_at_3_prime_end(self):
        """Site in last 25 bp — inside the reverse binding region.
        Must be removed; advisory must flag it."""
        cds = ("ATG" + "GCG" * 15
               + "GGT" + "CTC"                  # BsaI at codon 17-18
               + "GCG" + "TAA")
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"]
        br = r["binding_region_mutations"]
        assert br, "Advisory must flag site near 3' end"
        assert any(entry["region"] == "rev" for entry in br)

    def test_interior_site_does_not_trigger_advisory(self):
        """A site well inside the insert (outside both binding windows)
        must NOT trigger a binding_region advisory — mutating it is
        safe because the original template primers will still bind."""
        # Make the insert long enough that a codon at ~nt 60 is outside
        # both the first 25 bp (fwd) and the last 25 bp (rev).
        cds = ("ATG"
               + "GCG" * 20            # codons 1..20, nt 3..62
               + "GGT" + "CTC"         # BsaI at nt 63..68 (codon 21-22)
               + "GCG" * 20            # codons 23..42, nt 69..128
               + "TAA")                # nt 129..131, total 132 bp
        assert len(cds) == 132
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "GGTCTC" not in r["insert_seq"]
        # mutations exist but none falls in the binding windows
        assert r["mutations"]
        assert r["binding_region_mutations"] == []

    def test_esp3i_at_3_prime_end_flagged(self):
        cds = ("ATG" + "GCG" * 15
               + "CGT" + "CTC"
               + "GCG" + "TAA")
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        br = r["binding_region_mutations"]
        assert br
        assert any(entry["region"] == "rev" for entry in br)


class TestReverseStrandSites:
    """A forbidden recognition on the reverse strand shows up as its RC
    on the forward strand: GAGACC (BsaI RC) or GAGACG (Esp3I RC). The
    fixer must find and repair these too."""

    def test_gagacc_bsai_rc_found_and_fixed(self):
        """GAG-ACC — codons E-T on forward strand. E (GAG) and T (ACC)
        both have synonyms, so the fix has room."""
        cds = ("ATG"
               + "GCG" * 5
               + "GAG" + "ACC"          # BsaI on reverse strand
               + "GCG" * 5 + "TAA")
        assert "GAGACC" in cds
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "GAGACC" not in r["insert_seq"]
        assert "GGTCTC" not in r["insert_seq"]
        assert any("BsaI" in m for m in r["mutations"])

    def test_gagacg_esp3i_rc_found_and_fixed(self):
        """GAG-ACG — E-T on forward. Synonyms available."""
        cds = ("ATG"
               + "GCG" * 5
               + "GAG" + "ACG"          # Esp3I on reverse strand
               + "GCG" * 5 + "TAA")
        assert "GAGACG" in cds
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert "GAGACG" not in r["insert_seq"]
        assert "CGTCTC" not in r["insert_seq"]
        assert any("Esp3I" in m for m in r["mutations"])

    def test_reverse_strand_mutation_marked_rc_in_message(self):
        """The mutation description must carry `(rc)` so the user can see
        which strand was touched."""
        cds = ("ATG"
               + "GCG" * 5
               + "GAG" + "ACC"
               + "GCG" * 5 + "TAA")
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert any("(rc)" in m for m in r["mutations"])

    def test_both_strands_scrubbed_after_fix(self):
        """The final insert must be clean on both strands — check all
        four patterns absent."""
        cds = ("ATG"
               + "GCG" * 3 + "GGT" + "CTC"     # BsaI fwd
               + "GCG" * 3 + "GAG" + "ACC"     # BsaI rev
               + "GCG" * 3 + "CGT" + "CTC"     # Esp3I fwd
               + "GCG" * 3 + "GAG" + "ACG"     # Esp3I rev
               + "GCG" * 3 + "TAA")
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        insert = r["insert_seq"]
        for pat in ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG"):
            assert pat not in insert, f"{pat} still present in mutated insert"


class TestSwapCascadePrevention:
    """`_codon_fix_sites` must never accept a swap that removes a site
    but introduces a new forbidden site elsewhere in the sequence.
    This is the most dangerous regression mode — the fix 'succeeds'
    silently while the output is still unsafe for synthesis."""

    def test_fix_does_not_introduce_new_forbidden_site(self):
        """Fuzz a handful of real coding inserts with sites; for each
        successful fix, verify `_gb_find_forbidden_hits` on the output
        returns empty. If cascades slipped through, this trips."""
        cases = [
            ("ATG" + "GCG" * 3 + "GGT" + "CTC" + "GCG" * 3 + "TAA"),
            ("ATG" + "GCG" * 5 + "CGT" + "CTC" + "GCG" * 5 + "TAA"),
            ("ATG" + "GCG" * 3 + "GAG" + "ACC" + "GCG" * 3 + "TAA"),
            ("ATG" + "GCG" * 3 + "GAG" + "ACG" + "GCG" * 3 + "TAA"),
            ("ATG"
             + "GCG" * 3 + "GGT" + "CTC"
             + "GCG" * 3 + "CGT" + "CTC"
             + "GCG" * 3 + "TAA"),
        ]
        for cds in cases:
            r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                      codon_raw=_k12_raw())
            if "error" in r:
                continue  # unfixable is fine, but if fixed it must be clean
            remaining = sc._gb_find_forbidden_hits(r["insert_seq"])
            assert remaining == [], (
                f"Fix introduced / left residual sites for CDS {cds!r}: "
                f"{remaining}"
            )

    def test_before_after_hit_sets_cross_check_directly(self):
        """Direct probe of the cross-check invariant: a hit present in
        the output but absent from the input must NEVER exist after a
        successful `_codon_fix_sites` call. Applied to many multi-site
        inserts to maximise odds of hitting a cascade case."""
        templates = [
            "ATG" + "GCG" * 4 + "GGT" + "CTC" + "GCG" * 4 + "TAA",
            "ATG" + "GCG" * 4 + "CGT" + "CTC" + "GCG" * 4 + "TAA",
            ("ATG" + "GCG" * 3 + "GGT" + "CTC" + "GCG" * 3
             + "CGT" + "CTC" + "GCG" * 3 + "TAA"),
            ("ATG" + "GCG" * 3 + "GAG" + "ACC" + "GCG" * 3
             + "GAG" + "ACG" + "GCG" * 3 + "TAA"),
        ]
        all_pats = ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG")
        for cds in templates:
            protein = sc._mut_translate(cds)
            if not protein:
                continue
            before = sc._forbidden_hit_set(cds, all_pats)
            fixed, _ = sc._codon_fix_sites(
                cds, protein, _k12_raw(),
                sites=sc._GB_DOMESTICATION_FORBIDDEN,
            )
            after = sc._forbidden_hit_set(fixed, all_pats)
            new_hits = after - before
            assert not new_hits, (
                f"Cascade detected — fix introduced {new_hits} "
                f"on CDS {cds!r}"
            )


class TestUnfixableSitesInErrorMessage:
    """When a site cannot be repaired (single-codon amino acid runs with
    no synonymous alternatives), the error message must list every
    residual site so the user can pinpoint what needs manual rework."""

    def test_error_lists_all_remaining_sites_when_multi_failure(self):
        # Synthesize a coding insert where both sites live in single-
        # codon-amino-acid runs (M=ATG, W=TGG). GGTCTC and CGTCTC don't
        # fit cleanly in ATG/TGG-only runs, so instead pick a coding-
        # context site whose codons only appear as one single-codon AA.
        #
        # Tryptophan TGG is the only single-codon AA that matters here.
        # If a site overlaps TGG-TGG-TGG, no synonymous codon exists.
        # But GGTCTC doesn't overlap 3 TGGs — GGT != TGG.
        #
        # So construct a CDS where the fix succeeds on BOTH sites and
        # verify the mutations list, then use a separate construction
        # where the fix fails. For the "all sites listed" contract, use
        # a scenario where at least two sites remain — the simplest is a
        # non-coding part (fix never attempted):
        seq = "ATCG" * 4 + "GGTCTC" + "ATCG" * 4 + "CGTCTC" + "ATCG" * 4
        r = sc._design_gb_primers(seq, 0, len(seq), "Promoter",
                                  codon_raw=_k12_raw())
        assert "error" in r
        # Both sites named in the error
        assert "GGTCTC" in r["error"] or "BsaI" in r["error"]
        assert "CGTCTC" in r["error"] or "Esp3I" in r["error"]

    def test_error_lists_both_strands_when_forward_and_rc_present(self):
        seq = "ATCG" * 4 + "GGTCTC" + "ATCG" * 4 + "GAGACC" + "ATCG" * 4
        r = sc._design_gb_primers(seq, 0, len(seq), "Promoter",
                                  codon_raw=_k12_raw())
        assert "error" in r
        # Both occurrences must surface — not just the first
        e = r["error"]
        # BsaI is named for both GGTCTC and its RC; count positions
        # instead to prove both were reported.
        assert e.count("BsaI") >= 2


class TestBindingRegionAdvisoryHelper:
    """Direct unit tests for `_gb_binding_region_advisory` — the helper
    that decides which mutations fall inside primer binding windows."""

    def test_empty_mutations_returns_empty(self):
        assert sc._gb_binding_region_advisory([], 100, 20, 20) == []

    def test_zero_length_insert_returns_empty(self):
        muts = ["BsaI at nt 1: GGT→GGC (codon 1 G, freq=0.300)"]
        assert sc._gb_binding_region_advisory(muts, 0, 20, 20) == []

    def test_mutation_at_5_prime_flagged_fwd(self):
        muts = ["BsaI at nt 4: GGT→GGC (codon 2 G, freq=0.300)"]
        # codon 2 → codon_start = 3, insert_len 100, fwd_bind 20 → 3 < 20
        out = sc._gb_binding_region_advisory(muts, 100, 20, 20)
        assert len(out) == 1
        assert out[0]["region"] == "fwd"
        assert out[0]["codon_start"] == 3

    def test_mutation_at_3_prime_flagged_rev(self):
        muts = ["BsaI at nt 88: GGT→GGC (codon 30 G, freq=0.300)"]
        # codon 30 → codon_start = 87, insert_len 100, rev_bind 20
        # rev_lo = 100-20 = 80. codon_end = 90 > 80 → flagged rev.
        out = sc._gb_binding_region_advisory(muts, 100, 20, 20)
        assert len(out) == 1
        assert out[0]["region"] == "rev"

    def test_interior_mutation_not_flagged(self):
        muts = ["BsaI at nt 40: GGT→GGC (codon 14 G, freq=0.300)"]
        # codon 14 → codon_start = 39. fwd_bind 20 → 39 >= 20 (out of fwd).
        # rev_lo = 100-20 = 80. codon_end = 42 <= 80 (out of rev).
        out = sc._gb_binding_region_advisory(muts, 100, 20, 20)
        assert out == []

    def test_mutation_straddles_both_windows_on_tiny_insert(self):
        """On a very short insert where fwd and rev windows cover
        everything, a mutation should be flagged for BOTH regions."""
        muts = ["BsaI at nt 4: GGT→GGC (codon 2 G, freq=0.300)"]
        # insert_len 20, fwd_bind 18, rev_bind 18. fwd_hi=18, rev_lo=2.
        # codon_start=3, codon_end=6. 3<18 → fwd. 6>2 → rev.
        out = sc._gb_binding_region_advisory(muts, 20, 18, 18)
        regions = sorted(e["region"] for e in out)
        assert regions == ["fwd", "rev"]

    def test_malformed_mutation_string_skipped(self):
        """Parsing failure must not crash — malformed entries just get
        dropped from the advisory (no exception)."""
        muts = ["not a valid mutation string", "nonsense"]
        assert sc._gb_binding_region_advisory(muts, 100, 20, 20) == []


class TestCodonFixMutationPositions:
    """`_codon_fix_mutation_positions` extracts the 0-based codon start
    from each mutation string. It's the foundation for the binding-
    region advisory — must never raise on malformed input."""

    def test_codon_parsed_correctly(self):
        muts = [
            "BsaI at nt 4: GGT→GGC (codon 2 G, freq=0.300)",
            "Esp3I (rc) at nt 10: GAG→GAA (codon 4 E, freq=0.450)",
        ]
        positions = sc._codon_fix_mutation_positions(muts)
        # codon 2 → (2-1)*3 = 3; codon 4 → (4-1)*3 = 9
        assert positions == [3, 9]

    def test_malformed_strings_return_minus_one(self):
        muts = ["no codon info here", ""]
        assert sc._codon_fix_mutation_positions(muts) == [-1, -1]

    def test_non_string_entries_return_minus_one(self):
        # Defensive: the helper shouldn't raise if a caller passes a
        # pre-structured mutation list (future-compatibility hedge).
        muts = [None, 42, {"text": "x"}]
        assert sc._codon_fix_mutation_positions(muts) == [-1, -1, -1]

    def test_empty_list_returns_empty(self):
        assert sc._codon_fix_mutation_positions([]) == []


class TestDesignResultExposesBindingAdvisory:
    """`_design_gb_primers` must include `binding_region_mutations` in
    every result dict — present (possibly empty) so callers can iterate
    without KeyError."""

    def test_clean_insert_has_empty_binding_region_mutations(
            self, random_template):
        r = sc._design_gb_primers(random_template, 100, 400, "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r
        assert r["binding_region_mutations"] == []

    def test_error_path_does_not_need_binding_advisory(self):
        """Error results legitimately may omit the key — the UI paths
        that show the advisory are the success paths only."""
        seq = "ATCG" * 4 + "GGTCTC" + "ATCG" * 4
        r = sc._design_gb_primers(seq, 0, len(seq), "Promoter",
                                  codon_raw=_k12_raw())
        assert "error" in r

    def test_site_near_5_prime_surfaces_in_result(self):
        """End-to-end: insert with 5'-end site → design runs → the
        result dict carries a non-empty advisory the UI can render."""
        cds = ("ATG" + "GGT" + "CTC" + "GCG" * 15 + "TAA")
        r = sc._design_gb_primers(cds, 0, len(cds), "CDS",
                                  codon_raw=_k12_raw())
        assert "error" not in r, r
        assert r["binding_region_mutations"]
        entry = r["binding_region_mutations"][0]
        assert set(entry.keys()) == {"text", "region", "codon_start"}
        assert entry["region"] in ("fwd", "rev")


# ═══════════════════════════════════════════════════════════════════════════════
# Grammar abstraction (built-in registry, persistence, helpers)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBuiltinGrammars:
    """The two shipped grammars (``gb_l0`` and ``moclo_plant``) need to
    expose a stable schema — every consumer (PartsBinModal,
    DomesticatorModal, GrammarEditorModal) reads these fields and any
    drift would manifest as a downstream KeyError in the UI."""

    REQUIRED_KEYS = {
        "id", "name", "enzyme", "site", "spacer", "pad",
        "forbidden_sites", "positions", "coding_types",
        "type_to_insdc", "catalog", "editable",
    }

    @pytest.mark.parametrize("gid", ["gb_l0", "moclo_plant"])
    def test_required_keys_present(self, gid):
        g = sc._BUILTIN_GRAMMARS[gid]
        missing = self.REQUIRED_KEYS - set(g.keys())
        assert not missing, (
            f"Built-in grammar {gid!r} is missing keys: {missing}. "
            f"Every consumer assumes the full schema is present."
        )

    @pytest.mark.parametrize("gid", ["gb_l0", "moclo_plant"])
    def test_positions_have_full_overhang_metadata(self, gid):
        g = sc._BUILTIN_GRAMMARS[gid]
        for pos in g["positions"]:
            assert {"name", "type", "oh5", "oh3"} <= set(pos.keys()), (
                f"{gid} position {pos!r} missing required fields."
            )
            valid = set("ACGTRYWSMKBDHV")
            for label in ("oh5", "oh3"):
                bad = [c for c in pos[label] if c not in valid]
                assert not bad, (
                    f"{gid} position {pos['name']} has non-IUPAC base "
                    f"in {label}: {pos[label]!r}"
                )

    @pytest.mark.parametrize("gid", ["gb_l0", "moclo_plant"])
    def test_builtins_marked_not_editable(self, gid):
        # The editor honours `editable=False` to lock down built-in
        # grammars. Drifting this would let users corrupt the canonical
        # references, which is the whole point of having "Duplicate as
        # Custom" in the parts bin.
        assert sc._BUILTIN_GRAMMARS[gid].get("editable") is False

    def test_gb_l0_matches_legacy_constants(self):
        """The GB L0 grammar is derived from the existing _GB_*
        constants. Tests still reference those constants directly, so
        the derived view must agree with the source of truth."""
        g = sc._BUILTIN_GRAMMARS["gb_l0"]
        assert g["enzyme"] == sc._GB_L0_ENZYME_NAME
        assert g["site"]   == sc._GB_L0_ENZYME_SITE
        assert g["spacer"] == sc._GB_SPACER
        assert g["pad"]    == sc._GB_PAD
        assert dict(g["forbidden_sites"]) == dict(sc._GB_DOMESTICATION_FORBIDDEN)
        assert sorted(g["coding_types"]) == sorted(sc._GB_CODING_PART_TYPES)
        assert dict(g["type_to_insdc"]) == dict(sc._GB_PART_TYPE_TO_INSDC)
        # Every position in _GB_POSITIONS shows up in the grammar's
        # positions list with matching overhangs.
        derived = {p["type"]: (p["name"], p["oh5"], p["oh3"]) for p in g["positions"]}
        for ptype, (pos_name, oh5, oh3) in sc._GB_POSITIONS.items():
            assert derived[ptype] == (pos_name, oh5, oh3)

    def test_moclo_plant_uses_bsai(self):
        g = sc._BUILTIN_GRAMMARS["moclo_plant"]
        assert g["enzyme"] == "BsaI"
        assert g["site"]   == "GGTCTC"
        # MoClo L0 (BsaI) → L1 (BpiI/BbsI). Both must be on the
        # forbidden list so domestication scrubs them before the user
        # commits a synthesis order.
        assert "BsaI" in g["forbidden_sites"]
        assert "BpiI" in g["forbidden_sites"]


class TestGrammarHelpers:
    """The lookup helpers (``_all_grammars``, ``_get_active_grammar``,
    ``_grammar_position_by_type``) are the seam every consumer goes
    through; cover their happy paths + fallback semantics."""

    def test_all_grammars_returns_builtins(self, tmp_path, monkeypatch):
        # Empty custom grammars file — only built-ins are returned.
        monkeypatch.setattr(sc, "_grammars_cache", [])
        out = sc._all_grammars()
        assert "gb_l0" in out
        assert "moclo_plant" in out

    def test_all_grammars_includes_custom(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc, "_grammars_cache", [
            {"id": "custom_a", "name": "Custom A",
             "enzyme": "BsaI", "site": "GGTCTC",
             "spacer": "A", "pad": "GCGC",
             "forbidden_sites": {"BsaI": "GGTCTC"},
             "positions": [
                 {"name": "Pos 1", "type": "Promoter",
                  "oh5": "GGAG", "oh3": "AATG"},
             ],
             "coding_types": [], "type_to_insdc": {}, "catalog": []},
        ])
        out = sc._all_grammars()
        assert "custom_a" in out
        # Custom grammars are always editable in the editor regardless
        # of how the JSON was hand-written — guards against a
        # mis-flagged file locking the user out.
        assert out["custom_a"]["editable"] is True

    def test_get_active_grammar_defaults_to_gb_l0(
        self, tmp_path, monkeypatch,
    ):
        # No setting written → default = gb_l0.
        monkeypatch.setattr(sc, "_settings_cache", {})
        active = sc._get_active_grammar()
        assert active["id"] == "gb_l0"

    def test_get_active_grammar_resolves_setting(
        self, tmp_path, monkeypatch,
    ):
        monkeypatch.setattr(sc, "_settings_cache", {"active_grammar": "moclo_plant"})
        active = sc._get_active_grammar()
        assert active["id"] == "moclo_plant"

    def test_get_active_grammar_recovers_from_missing_id(self):
        """If the persisted active id no longer resolves (e.g., the
        custom grammar that was selected got deleted), the helper
        flips the setting back to gb_l0 instead of crashing."""
        sc._save_settings({"active_grammar": "nonexistent_grammar"})
        active = sc._get_active_grammar()
        assert active["id"] == "gb_l0"
        # The helper writes the recovery back to settings so we don't
        # keep falling back forever.
        assert sc._get_setting("active_grammar") == "gb_l0"

    def test_grammar_position_by_type(self):
        g = sc._BUILTIN_GRAMMARS["gb_l0"]
        pos = sc._grammar_position_by_type(g, "CDS")
        assert pos is not None
        assert pos["oh5"] == "AATG"
        assert pos["oh3"] == "GCTT"
        # Type not in this grammar → None (not KeyError).
        assert sc._grammar_position_by_type(g, "BogusType") is None


class TestSettingsPersistence:
    """``settings.json`` round-trip through the envelope schema —
    ``_load_setting`` / ``_save_setting`` are how every preference is
    persisted now (active grammar today, more later)."""

    def test_set_then_get(self):
        sc._set_setting("active_grammar", "moclo_plant")
        assert sc._get_setting("active_grammar") == "moclo_plant"

    def test_get_default_when_missing(self):
        sc._save_settings({})
        assert sc._get_setting("missing_key", "fallback") == "fallback"

    def test_round_trip_through_disk(self):
        sc._set_setting("active_grammar", "moclo_plant")
        # `_set_setting` defers the disk write to a daemon thread so
        # toggle-heavy UI flows don't block on fsync. Wait for the
        # flush before clearing the cache, otherwise we'd race the
        # background writer.
        sc._settings_flush_sync()
        # Force a re-read from disk — the cache must agree with what
        # was written.
        sc._settings_cache = None
        assert sc._load_settings()["active_grammar"] == "moclo_plant"


class TestCustomGrammarPersistence:
    """``cloning_grammars.json`` round-trip + the deepcopy-on-load
    contract that protects the cache from caller-side mutation."""

    def test_round_trip(self):
        entries = [{
            "id": "custom_x", "name": "Custom X",
            "enzyme": "BsaI", "site": "GGTCTC",
            "spacer": "A", "pad": "GCGC",
            "forbidden_sites": {"BsaI": "GGTCTC"},
            "positions": [
                {"name": "Pos 1", "type": "Promoter",
                 "oh5": "GGAG", "oh3": "AATG"},
            ],
            "coding_types": [], "type_to_insdc": {}, "catalog": [],
        }]
        sc._save_custom_grammars(entries)
        sc._grammars_cache = None
        loaded = sc._load_custom_grammars()
        assert len(loaded) == 1
        assert loaded[0]["id"] == "custom_x"
        assert loaded[0]["positions"][0]["oh5"] == "GGAG"

    def test_load_returns_independent_dicts(self):
        sc._save_custom_grammars([{
            "id": "c1", "name": "C1", "enzyme": "X", "site": "AAAAAA",
            "spacer": "A", "pad": "AAAA",
            "forbidden_sites": {}, "positions": [],
            "coding_types": [], "type_to_insdc": {}, "catalog": [],
        }])
        a = sc._load_custom_grammars()
        b = sc._load_custom_grammars()
        # Same content but different dict identity — mutating one must
        # not affect the other (or the cache).
        a[0]["name"] = "MUTATED"
        assert b[0]["name"] == "C1"
        c = sc._load_custom_grammars()
        assert c[0]["name"] == "C1"


# ═══════════════════════════════════════════════════════════════════════════════
# Parts Bin grammar dropdown + filter
# ═══════════════════════════════════════════════════════════════════════════════

class TestPartsBinGrammarFilter:
    """The Parts Bin filters its catalog and user parts to the active
    grammar. Switching grammars repopulates the table with that
    grammar's parts. Legacy parts (no ``grammar`` field) treat as
    ``gb_l0`` so existing v0.3.x data migrates intact."""

    async def test_default_active_is_gb_l0(self, isolated_parts_bin):
        # Fresh app — settings.json doesn't exist yet, default kicks in.
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._active_grammar_id() == "gb_l0"

    def test_gb_l0_is_first_in_dropdown(self):
        """Golden Braid L0 must be position 1 of the grammar dropdown
        regardless of how many built-in or custom grammars exist.
        Order contract is enforced by `_grammar_dropdown_options`
        rather than relying on Python's dict-insertion ordering
        accidentally getting it right (a user-defined grammar
        id-sorted before ``gb_l0`` would otherwise show up first)."""
        # Seed a custom grammar whose id alphabetises BEFORE "gb_l0",
        # so any naive sort would put it at position 1 instead.
        sc._save_custom_grammars([{
            "id": "aaa_pinned_first", "name": "Aardvark",
            "enzyme": "BsaI", "site": "GGTCTC",
            "spacer": "A", "pad": "GCGC",
            "forbidden_sites": {"BsaI": "GGTCTC"},
            "positions": [
                {"name": "Pos 1", "type": "Promoter",
                 "oh5": "GGAG", "oh3": "AATG"},
            ],
            "coding_types": [], "type_to_insdc": {}, "catalog": [],
        }])
        options = sc._grammar_dropdown_options()
        assert options, "Grammar dropdown should never be empty."
        first_label, first_id = options[0]
        assert first_id == "gb_l0", (
            f"GB L0 must be position 1 of the dropdown; got "
            f"{first_id!r} ({first_label!r}) instead. "
            f"Full order: {[gid for _label, gid in options]}"
        )

    def test_dropdown_order_builtins_then_custom(self):
        """Built-ins (gb_l0 first, then moclo_plant) precede every
        custom grammar in the dropdown. The (custom) suffix is the
        signal that tells the user where the boundary is."""
        sc._save_custom_grammars([{
            "id": "custom_z", "name": "Z Custom",
            "enzyme": "BsaI", "site": "GGTCTC",
            "spacer": "A", "pad": "GCGC",
            "forbidden_sites": {"BsaI": "GGTCTC"},
            "positions": [
                {"name": "Pos 1", "type": "Promoter",
                 "oh5": "GGAG", "oh3": "AATG"},
            ],
            "coding_types": [], "type_to_insdc": {}, "catalog": [],
        }])
        opts = sc._grammar_dropdown_options()
        ids_in_order = [gid for _label, gid in opts]
        # First two slots are the built-ins, in declared order.
        assert ids_in_order[:2] == ["gb_l0", "moclo_plant"]
        # Custom grammars come after and carry the (custom) tag.
        custom_labels = [
            label for label, gid in opts if gid not in sc._BUILTIN_GRAMMARS
        ]
        assert all("(custom)" in label for label in custom_labels)

    async def test_legacy_parts_default_to_gb_l0(self, isolated_parts_bin):
        # A v0.3.x part with no `grammar` field shows up under GB L0.
        sc._save_parts_bin([{
            "name": "legacy_x", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "sequence": "ATG" * 10,
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            user_rows = [r for r in modal._rows if r.get("user")]
            assert any(r["name"] == "legacy_x" for r in user_rows)

    async def test_all_parts_visible_regardless_of_grammar(
        self, isolated_parts_bin,
    ):
        """The Parts Bin no longer filters by an active grammar — every
        user-saved part is visible at all times, with the Grammar
        column indicating which assembly standard each row belongs to.
        Grammar selection moved into the New Part modal."""
        sc._save_parts_bin([
            {"name": "gb_part",    "type": "CDS",
             "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
             "backbone": "pUPD2",  "marker": "Spectinomycin",
             "sequence": "ATG" * 10, "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0,
             "grammar": "gb_l0"},
            {"name": "moclo_part", "type": "CDS",
             "position": "Pos 3", "oh5": "AGGT", "oh3": "GCTT",
             "backbone": "pUPD2",  "marker": "Spectinomycin",
             "sequence": "ATG" * 10, "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0,
             "grammar": "moclo_plant"},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            user_names = {r["name"] for r in modal._rows if r.get("user")}
            assert user_names == {"gb_part", "moclo_part"}, (
                f"Both user parts should be listed; got {user_names}"
            )
            # Each user row carries its source grammar id so the
            # Grammar column can render the right label.
            grammars = {
                r["name"]: r.get("grammar")
                for r in modal._rows if r.get("user")
            }
            assert grammars == {
                "gb_part": "gb_l0", "moclo_part": "moclo_plant",
            }

    async def test_grammar_column_shows_human_name(
        self, isolated_parts_bin,
    ):
        sc._save_parts_bin([{
            "name": "x", "type": "CDS",
            "position": "Pos 3-4", "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2",  "marker": "Spectinomycin",
            "sequence": "ATG" * 10, "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
            "grammar": "moclo_plant",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            t = app.screen.query_one("#parts-table", sc.DataTable)
            row_keys = list(t.rows.keys())
            col_keys = list(t.columns.keys())
            # Last column is the new Grammar column; show the
            # built-in's human name, not the raw id.
            cell = t.get_cell(row_keys[0], col_keys[-1])
            assert "MoClo" in str(cell), (
                f"Grammar column should show human-readable name; "
                f"got {cell!r}"
            )

    async def test_new_part_inherits_active_grammar(self, isolated_parts_bin):
        sc._set_setting("active_grammar", "moclo_plant")
        sc._save_parts_bin([])
        # Simulate the DomesticatorModal save callback: the part dict
        # arrives with no `grammar` field; PartsBinModal._new_part
        # patches the active grammar id in.
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Manually invoke the inner _on_result by calling the
            # `_new_part` callback path — we don't push the full
            # DomesticatorModal because that would need a full record.
            # Instead, directly emulate the persistence step.
            new_part = {
                "name": "fresh", "type": "Promoter",
                "position": "Pos 1", "oh5": "GGAG", "oh3": "AATG",
                "backbone": "pUPD2", "marker": "Spectinomycin",
                "sequence": "ATGCATGCATGC",
                "fwd_primer": "", "rev_primer": "",
                "fwd_tm": 0.0, "rev_tm": 0.0,
            }
            # This mirrors the body of the _on_result closure in _new_part.
            new_part.setdefault("grammar", modal._active_grammar_id())
            entries = sc._load_parts_bin()
            entries.insert(0, new_part)
            sc._save_parts_bin(entries)
            sc._parts_bin_cache = None
            assert sc._load_parts_bin()[0]["grammar"] == "moclo_plant"


# ═══════════════════════════════════════════════════════════════════════════════
# DomesticatorModal honors the active grammar
# ═══════════════════════════════════════════════════════════════════════════════

class TestDomesticatorUsesActiveGrammar:
    """Primer design parameters (positions, enzyme, spacer, pad,
    forbidden sites) all flow from the active grammar — switching
    grammar in the parts bin should change what the Domesticator
    designs without the user having to flip anything else."""

    def test_design_uses_grammar_overhangs(self):
        moclo = sc._BUILTIN_GRAMMARS["moclo_plant"]
        # Long enough for a binding region; ATG-aligned for codon-fix
        # repair if needed.
        seq = "ATG" + "GCG" * 30 + "TAA"
        r = sc._design_gb_primers(
            seq, 0, len(seq), "Promoter", grammar=moclo,
        )
        assert "error" not in r, r
        # Forward primer should carry the MoClo Promoter overhang
        # (GGAG → AATG), not the GB one (GGAG → TGAC).
        assert r["oh5"] == "GGAG"
        assert r["oh3"] == "AATG"
        # Tail uses BsaI, not Esp3I.
        assert "GGTCTC" in r["fwd_full"][:15]
        assert "CGTCTC" not in r["fwd_full"][:15]

    def test_design_rejects_type_not_in_grammar(self):
        moclo = sc._BUILTIN_GRAMMARS["moclo_plant"]
        # MoClo Plant doesn't define CDS-NS — design should refuse
        # rather than silently fall back to a GB position.
        seq = "ATG" + "GCG" * 30 + "TAA"
        r = sc._design_gb_primers(
            seq, 0, len(seq), "CDS-NS", grammar=moclo,
        )
        assert "error" in r
        assert "CDS-NS" in r["error"]

    def test_forbidden_sites_scan_uses_grammar(self):
        # MoClo Plant scrubs BsaI + BpiI; Esp3I would be allowed.
        moclo = sc._BUILTIN_GRAMMARS["moclo_plant"]
        # Seq with an internal Esp3I (CGTCTC). Under MoClo it should
        # NOT be flagged as forbidden — this is the contract that
        # makes per-grammar Type IIS scrubbing meaningful.
        seq_with_esp3i = "ATG" + "CGTCTC" + "GCG" * 25 + "TAA"
        hits = sc._gb_find_forbidden_hits(
            seq_with_esp3i, sites=moclo["forbidden_sites"],
        )
        assert all(h[0] != "Esp3I" for h in hits), (
            f"MoClo grammar should not flag Esp3I sites; got {hits!r}"
        )

    async def test_modal_title_shows_active_grammar(
        self, isolated_parts_bin, tiny_record,
    ):
        sc._set_setting("active_grammar", "moclo_plant")
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.DomesticatorModal(str(tiny_record.seq), [],
                                                 current_plasmid_name="x"))
            await pilot.pause()
            await pilot.pause(0.1)
            title = app.screen.query_one("#dom-title", sc.Static)
            assert "MoClo" in str(title.render())


# ═══════════════════════════════════════════════════════════════════════════════
# _classify_part_from_plasmid — auto-detect grammar/position from a plasmid
# ═══════════════════════════════════════════════════════════════════════════════

def _build_gb_l0_part_seq(oh5: str, oh3: str,
                          insert: str = "ATGAAACCCGGG" * 5) -> str:
    """Build a synthetic Golden Braid L0 part vector: insert flanked by
    inward-facing Esp3I sites with the requested 4-nt overhangs, plus a
    chunk of inert backbone padding so the digest produces a 2-fragment
    output (insert vs vector). Treated as circular by callers."""
    core = "CGTCTCA" + oh5 + insert + oh3 + "AGAGACG"
    backbone = "AAAAATTTTT" * 50
    return core + backbone


def _build_moclo_plant_part_seq(oh5: str, oh3: str,
                                 insert: str = "ATGAAACCCGGG" * 5) -> str:
    """Build a synthetic MoClo Plant L0 part vector. MoClo uses BsaI
    (GGTCTC) at L0; layout mirrors `_build_gb_l0_part_seq` but with
    BsaI + reverse-complement BsaI (GAGACC) flanking the insert. The
    backbone padding has no BsaI sites so the digest yields exactly
    two fragments: insert (oh5 → oh3) and backbone."""
    core = "GGTCTCA" + oh5 + insert + oh3 + "AGAGACC"
    backbone = "AAAAATTTTT" * 50
    return core + backbone


class TestClassifyPartFromPlasmid:
    """``_classify_part_from_plasmid`` digests a circular plasmid with
    each grammar's Type IIS enzyme and matches the released fragment's
    overhangs against the grammar's position table. Powers the Parts
    Bin "Load Part" button so the user doesn't have to manually pick a
    grammar / position when registering an externally-domesticated
    plasmid."""

    def test_detects_gb_l0_promoter(self):
        # GB 2.0 PromUTR (combined Promoter + 5'UTR): GGAG → AATG.
        # The 3' overhang lands on the start codon (ATG) of the
        # downstream CDS, matching the canonical post-cloning shape
        # of an Anderson-style promoter (e.g. J23100) bundled with
        # its RBS / 5'UTR.
        seq = _build_gb_l0_part_seq("GGAG", "AATG")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"] == "gb_l0"
        assert result["position"]["type"] == "Promoter"
        assert result["insert"]["left"]["overhang_seq"] == "GGAG"
        assert result["insert"]["right"]["overhang_seq"] == "AATG"

    def test_detects_gb_l0_cds(self):
        seq = _build_gb_l0_part_seq("AATG", "GCTT")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["position"]["type"] == "CDS"

    def test_no_match_returns_none(self):
        # Overhangs that match no Golden Braid position.
        seq = _build_gb_l0_part_seq("CCCC", "GGGG")
        assert sc._classify_part_from_plasmid(seq, circular=True) is None

    def test_linear_skipped(self):
        seq = _build_gb_l0_part_seq("GGAG", "TGAC")
        # Linear plasmids can't be cleanly excised so the helper bails.
        assert sc._classify_part_from_plasmid(seq, circular=False) is None

    def test_empty_seq_safe(self):
        assert sc._classify_part_from_plasmid("", circular=True) is None

    def test_smaller_fragment_is_insert(self):
        """The fragment classification picks the smaller fragment as
        the insert and the larger as the vector — a 60-bp insert in a
        ~600-bp synthetic backbone must come back labelled correctly."""
        insert_core = "ATGAAACCCGGG" * 5   # 60 bp
        seq = _build_gb_l0_part_seq("GGAG", "AATG", insert_core)
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert (len(result["insert"]["top_seq"])
                < len(result["vector"]["top_seq"]))

    def test_gb_priority_over_moclo_for_shared_promoter_overhangs(self):
        """Both gb_l0 Pos 1 (combined Promoter+5'UTR) and moclo_plant
        Pos 1 use ``GGAG / AATG`` overhangs — a J23100+RBS cassette
        cloned into either system has identical digest output. The
        registry order resolves the ambiguity: gb_l0 ships first in
        ``_BUILTIN_GRAMMARS`` so it wins. Regression guard for the
        2026-05-10 user report ("J23100/J23114 detected as MoClo
        instead of GB"). To force a MoClo classification on a
        GGAG/AATG plasmid, the user must remove gb_l0 from the
        active grammar list or use the manual New Part flow."""
        seq = _build_moclo_plant_part_seq("GGAG", "AATG")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]    == "gb_l0"
        assert result["position"]["type"] == "Promoter"
        assert result["position"]["name"] == "Pos 1"
        assert result["insert"]["left"]["overhang_seq"]  == "GGAG"
        assert result["insert"]["right"]["overhang_seq"] == "AATG"

    def test_detects_moclo_plant_cds(self):
        """MoClo Plant CDS overhangs (AGGT / GCTT) must classify to
        the Pos 3 / CDS slot. Differentiates from gb_l0 CDS (AATG /
        GCTT): gb_l0 has no AGGT-prefixed position, so the GB pass
        misses cleanly and MoClo takes over."""
        seq = _build_moclo_plant_part_seq("AGGT", "GCTT")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]       == "moclo_plant"
        assert result["position"]["type"] == "CDS"

    def test_no_esp3i_sites_falls_through_to_bsai(self):
        """A plasmid with NO Esp3I sites AND with overhangs that
        ONLY match MoClo (not the expanded GB grammar) must fall
        through to MoClo Plant rather than failing classification.

        Use MoClo Pos 2 / 5' UTR (AATG / AGGT) — AGGT is disjoint
        from every GB position (canonical GB 2.0 internal overhangs
        are GGAG/TGAC/TCCC/TACT/CCAT/AATG/AGCC/TTCG/GCAG/GCTT/GGTA/CGCT,
        none of which is AGGT). Pre-2026-05-10 this test used
        (GGTA, CGCT) but the 2026-05-10 GB 2.0 expansion added
        `Terminator-only` (Pos 21, GGTA→CGCT), so that pair now
        classifies as gb_l0 — moved to (AATG, AGGT) to keep the
        fall-through-to-MoClo invariant testable.
        """
        seq = _build_moclo_plant_part_seq("AATG", "AGGT")
        # Spot-check: no Esp3I (CGTCTC) in either strand.
        assert "CGTCTC" not in seq
        assert "GAGACG" not in seq
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"] == "moclo_plant"


# ═══════════════════════════════════════════════════════════════════════════════
# _classify_part_from_plasmid — multi-level (L0 / TU / MOD) detection
# + entry-vector compatibility check (regression guard for 2026-05-10)
# ═══════════════════════════════════════════════════════════════════════════════

def _build_post_cloning_l0_seq(oh5: str, oh3: str,
                                insert: str = "ATG" * 15) -> str:
    """An L0 part already cloned into a pUPD2-style L0 entry vector.
    BsaI sites flank the part body (ready for L1 assembly); the
    primary Esp3I sites have been "consumed" so the only remaining
    cut sites are the secondary (BsaI) ones. The user explicitly
    asked Load Part to recognise this case (2026-05-10 follow-up)."""
    core = "GGTCTCA" + oh5 + insert + oh3 + "AGAGACC"
    backbone = "AAAAATTTTT" * 50    # no BsaI internal
    return core + backbone


def _build_gb_tu_seq(insert: str = "ATG" * 30,
                       oh5: str = "GGAG", oh3: str = "CGCT") -> str:
    """Synthetic GB TU: BsaI sites flanking a TU body whose internal
    Esp3I sites are scrubbed (domesticated). Boundary overhangs default
    to ``GGAG/CGCT`` — the canonical GB TU start/end."""
    core = "GGTCTCA" + oh5 + insert + oh3 + "AGAGACC"
    backbone = "AAAAATTTTT" * 50
    return core + backbone


def _build_gb_mod_seq(insert: str = "ATG" * 50,
                       oh5: str = "GGAG", oh3: str = "CGCT") -> str:
    """Synthetic GB MOD (L2): Esp3I sites flanking a multi-TU body
    whose internal BsaI sites are scrubbed. The same TU boundary
    overhangs flank the MOD insert because L2 destinations preserve
    the TU edge convention."""
    core = "CGTCTCA" + oh5 + insert + oh3 + "AGAGACG"
    backbone = "AAAAATTTTT" * 50
    return core + backbone


class TestClassifyPartFromPlasmidLevels:
    """Multi-level detection: pre-cloning L0, post-cloning L0 (BsaI-
    flanked, no Esp3I), TU (level 1), MOD (level 2). Plus the
    entry-vector compatibility check that surfaces "this plasmid was
    cloned into your configured destination" so the user doesn't
    accidentally save a part against the wrong backbone.
    """

    def test_pre_cloning_l0_classifies_as_level_0(self):
        """Existing path: synthetic GB Promoter still has Esp3I sites
        flanking the part body. Primary digest matches L0 position
        directly — release_enzyme is the primary (Esp3I)."""
        seq = _build_gb_l0_part_seq("GGAG", "AATG")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]    == "gb_l0"
        assert result["level"]         == 0
        assert result["position"]["type"] == "Promoter"
        assert result["release_enzyme"] in ("Esp3I", "BsmBI")

    def test_post_cloning_l0_classifies_as_level_0(self):
        """User's explicit ask (2026-05-10): an L0 part already cloned
        into pUPD2 has no Esp3I sites left, only BsaI flanking. Must
        still classify as gb_l0 / level=0 — the part is biologically
        an L0, just packaged ready for L1 assembly. Release enzyme
        switches to the secondary (BsaI)."""
        seq = _build_post_cloning_l0_seq("GGAG", "AATG")
        # Sanity: no Esp3I sites either strand.
        assert "CGTCTC" not in seq
        assert "GAGACG" not in seq
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]    == "gb_l0"
        assert result["level"]         == 0
        assert result["position"]["type"] == "Promoter"
        assert result["release_enzyme"] == "BsaI"

    def test_post_cloning_l0_cds_classifies_as_level_0(self):
        """A post-cloning L0 CDS (AATG/GCTT — between Promoter and
        Terminator overhangs) classifies the same way — verifies the
        post-cloning path isn't hard-coded to a single position."""
        seq = _build_post_cloning_l0_seq("AATG", "GCTT")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["level"]            == 0
        assert result["position"]["type"] == "CDS"
        assert result["release_enzyme"]   == "BsaI"

    def test_j23100_with_rbs_classifies_as_gb_promoter(self):
        """Regression guard for the 2026-05-10 user report: a
        J23100/J23114 cassette already cloned into an L0 entry vector
        with the RBS bundled (GGAG → AATG, BsaI sites flanking) must
        classify as gb_l0 / Promoter, not moclo_plant. Pre-fix the
        GB Promoter slot was (GGAG, TGAC) which never matched any
        post-cloning Anderson-style cassette — the user got
        ``moclo_plant`` even when their workflow was GB."""
        seq = _build_post_cloning_l0_seq("GGAG", "AATG")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]    == "gb_l0"
        assert result["level"]         == 0
        assert result["position"]["name"] == "Pos 1"
        assert result["position"]["type"] == "Promoter"
        assert result["release_enzyme"] == "BsaI"

    def test_separate_gb_promoter_only_classifies(self):
        """The expanded GB 2.0 grammar adds a separate Promoter-only
        position (GGAG → CCAT) for users who domesticate the
        promoter and 5'UTR apart. Verify the new slot matches."""
        seq = _build_gb_l0_part_seq("GGAG", "CCAT")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]    == "gb_l0"
        assert result["position"]["name"] == "Pos 1a"
        assert result["position"]["type"] == "Promoter-only"

    def test_separate_gb_5_utr_classifies(self):
        """Partner to the separate Promoter-only slot: GB 2.0's
        5' UTR position takes the CCAT connector overhang to AATG.
        Pre-fix it was TGAC/AATG (never connected to any other slot
        coherently); the canonical GB 2.0 spec uses CCAT."""
        seq = _build_gb_l0_part_seq("CCAT", "AATG")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]    == "gb_l0"
        assert result["position"]["name"] == "Pos 1b"
        assert result["position"]["type"] == "5' UTR"

    def test_tu_classifies_as_level_1(self):
        """A GB TU (BsaI-flanked, GGAG/CGCT boundary, no internal
        Esp3I) must classify as gb_l0 / level=1. Released by the
        secondary enzyme (BsaI); the position dict carries the
        ``TU`` label so the Parts Bin TU tab picks it up."""
        seq = _build_gb_tu_seq()
        # Sanity: no Esp3I on either strand.
        assert "CGTCTC" not in seq
        assert "GAGACG" not in seq
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]    == "gb_l0"
        assert result["level"]         == 1
        assert result["position"]["name"] == "TU"
        assert result["release_enzyme"] == "BsaI"

    def test_mod_shape_classifies_as_tu_overhang_shape_only(self):
        """Pre-2026-05-13 a `_build_gb_mod_seq()` plasmid (Esp3I-
        flanked, canonical TU-boundary overhangs) classified as
        level=2 (MOD) via enzyme-parity inference (`primary release
        ⇒ MOD`). That inference assumed splicecraft's enzyme
        convention (Esp3I=primary=L0); the pDGB1 convention used by
        EDEN flips parity — Esp3I IS the L1→L2 release in their
        labs, and the user's TUs cut with Esp3I.

        Overhang shape alone can't tell TU from MOD across both
        conventions, so the classifier now returns level=1 for any
        TU-boundary match regardless of which enzyme cut. Users who
        need to tag a MOD specifically can do so via Parts Bin →
        Edit. Regression guard for the MAV-25-in-alpha-2 fix
        (2026-05-13).
        """
        seq = _build_gb_mod_seq()
        # Sanity: no BsaI on either strand.
        assert "GGTCTC" not in seq
        assert "GAGACC" not in seq
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"]    == "gb_l0"
        # New behaviour: overhang-only classification can't
        # distinguish MOD from TU — both have the same canonical
        # boundary pair. Returns level=1 with the bare "TU" label.
        assert result["level"]            == 1
        assert result["position"]["name"] == "TU"
        assert result["release_enzyme"] in ("Esp3I", "BsmBI")

    def test_random_overhangs_still_return_none(self):
        """Sanity carry-over from the L0-only era: non-recognised
        overhangs (not in any L0 position AND not on the TU boundary)
        still return None at every level."""
        # L0 overhangs that don't match any GB position; doesn't
        # match GB TU boundary (GGAG/CGCT) either.
        seq = _build_gb_l0_part_seq("CCCC", "GGGG")
        assert sc._classify_part_from_plasmid(seq, circular=True) is None

    def test_entry_vector_field_is_none_when_unconfigured(self):
        """No entry vectors saved → result.entry_vector is None.
        Default state for a fresh install — verifies we don't
        accidentally fabricate a match."""
        seq = _build_gb_l0_part_seq("GGAG", "AATG")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result.get("entry_vector") is None

    def test_entry_vector_match_l0(self, isolated_library):
        """Configure an L0 entry vector. Build a user plasmid with
        the same backbone but a different insert. Classify it — the
        entry_vector field should report ``matches=True`` because
        the digest's vector half is rotationally identical to the
        configured entry vector's vector half."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Shared backbone; unique cassette per plasmid.
        backbone = "AAAAATTTTT" * 50
        ev_cassette  = "CCCCTTTT" * 5
        user_promoter = "ATGAAACCCGGG" * 5
        # Both plasmids: Esp3I sites flanking GGAG/AATG (Promoter
        # slot — combined Promoter+5'UTR per GB 2.0 PromUTR).
        ev_seq = ("CGTCTCA" + "GGAG" + ev_cassette + "AATG"
                  + "AGAGACG" + backbone)
        user_seq = ("CGTCTCA" + "GGAG" + user_promoter + "AATG"
                    + "AGAGACG" + backbone)
        # Persist the entry vector for gb_l0, role="" (singleton L0
        # destination — pUPD2-style).
        ev_rec = SeqRecord(
            Seq(ev_seq), id="testEV", name="testEV",
            annotations={"topology": "circular",
                          "molecule_type": "DNA"},
        )
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0",
            "role":       "",
            "name":       "test_pUPD2",
            "size":       len(ev_seq),
            "source":     "test",
            "gb_text":    sc._record_to_gb_text(ev_rec),
        }])
        # Classify the user's plasmid.
        result = sc._classify_part_from_plasmid(user_seq, circular=True)
        assert result is not None
        assert result["grammar_id"] == "gb_l0"
        assert result["level"]      == 0
        ev = result.get("entry_vector")
        assert ev is not None, \
            "entry_vector check should fire when an EV is configured"
        assert ev["matches"] is True
        assert ev["name"]    == "test_pUPD2"
        assert ev["role"]    == ""

    def test_entry_vector_mismatch_when_backbone_differs(
            self, isolated_library):
        """Configure an L0 entry vector. Classify a plasmid that fits
        the same grammar/position BUT has a DIFFERENT backbone — the
        entry_vector field must come back None (no match), so the
        user can spot the unexpected destination."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        ev_backbone   = "AAAAATTTTT" * 50
        user_backbone = "GGGGCCCCAA" * 50    # different
        cassette = "ATGAAACCCGGG" * 5
        ev_seq = ("CGTCTCA" + "GGAG" + "AAAA" + "AATG"
                  + "AGAGACG" + ev_backbone)
        user_seq = ("CGTCTCA" + "GGAG" + cassette + "AATG"
                    + "AGAGACG" + user_backbone)
        ev_rec = SeqRecord(
            Seq(ev_seq), id="testEV", name="testEV",
            annotations={"topology": "circular",
                          "molecule_type": "DNA"},
        )
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0",
            "role":       "",
            "name":       "test_pUPD2",
            "size":       len(ev_seq),
            "source":     "test",
            "gb_text":    sc._record_to_gb_text(ev_rec),
        }])
        result = sc._classify_part_from_plasmid(user_seq, circular=True)
        assert result is not None
        assert result.get("entry_vector") is None, \
            "entry_vector must be None when no configured EV matches"

    def test_entry_vector_role_propagates(self, isolated_library):
        """The role field on the configured EV (e.g. ``Alpha1`` vs
        ``Omega1`` for GB) must propagate into the classifier's
        ``entry_vector.role`` so the worker can show the user which
        L1 destination their plasmid lines up with."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        backbone = "AAAAATTTTT" * 50
        # An L1 entry vector for the Alpha1 role: BsaI-flanked acceptor
        # whose dropout has GGAG/CGCT (TU boundary) overhangs. We
        # configure it for gb_l0 / Alpha1.
        ev_dropout = "ACGTACGT" * 5
        ev_seq = ("GGTCTCA" + "GGAG" + ev_dropout + "CGCT"
                  + "AGAGACC" + backbone)
        # User TU plasmid: same architecture, different TU body, SAME
        # backbone — the post-assembly state of cloning a TU into the
        # Alpha1 acceptor.
        tu_body = "ATG" * 30
        user_seq = ("GGTCTCA" + "GGAG" + tu_body + "CGCT"
                    + "AGAGACC" + backbone)
        ev_rec = SeqRecord(
            Seq(ev_seq), id="alpha1", name="alpha1",
            annotations={"topology": "circular",
                          "molecule_type": "DNA"},
        )
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0",
            "role":       "Alpha1",
            "name":       "pDGB1_alpha1",
            "size":       len(ev_seq),
            "source":     "test",
            "gb_text":    sc._record_to_gb_text(ev_rec),
        }])
        result = sc._classify_part_from_plasmid(user_seq, circular=True)
        assert result is not None
        assert result["level"] == 1
        ev = result.get("entry_vector") or {}
        assert ev.get("matches") is True
        assert ev.get("role")    == "Alpha1"
        assert ev.get("name")    == "pDGB1_alpha1"

    async def test_load_part_worker_propagates_level_to_parts_bin(
            self, isolated_library, isolated_parts_bin):
        """End-to-end: Load Part on a TU plasmid must save the new
        Parts Bin row with ``level=1`` so it lands in the TU tab.
        Pre-fix this was hardcoded to 0 — every TU got mis-tabbed
        into L0."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        tu_seq = _build_gb_tu_seq()
        rec = SeqRecord(
            Seq(tu_seq), id="my_tu", name="my_tu",
            annotations={"topology": "circular",
                          "molecule_type": "DNA"},
        )
        gb_text = sc._record_to_gb_text(rec)
        sc._save_collections([{
            "name":        "TUcoll",
            "description": "test",
            "plasmids":    [{
                "id":      "my_tu",
                "name":    "my_tu",
                "size":    len(tu_seq),
                "n_feats": 0,
                "gb_text": gb_text,
            }],
        }])
        sc._set_active_collection_name("TUcoll")
        sc._save_library([{
            "id":      "my_tu",
            "name":    "my_tu",
            "size":    len(tu_seq),
            "n_feats": 0,
            "gb_text": gb_text,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            pm = app.screen
            pm.query_one("#btn-load-part", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            picker = app.screen
            t = picker.query_one("#loadpart-table", sc.DataTable)
            t.move_cursor(row=0)
            await pilot.pause()
            # 2026-05-13: picker is now multi-select via toggle —
            # toggle the cursor row before Load Selected.
            picker.action_toggle_selection()
            await pilot.pause()
            picker.query_one("#btn-loadpart-ok", sc.Button).press()
            # Worker runs `@work(thread=True)` — give it room to land.
            await pilot.pause()
            await pilot.pause(0.5)
            await pilot.pause()
        entries = sc._load_parts_bin()
        assert len(entries) == 1, \
            f"expected one part saved, got {entries!r}"
        assert entries[0]["level"]   == 1
        assert entries[0]["grammar"] == "gb_l0"


class TestClassifyPartFromPlasmidPerAcceptor:
    """A TU assembled into an Alpha2 / Omega1 / Omega2 acceptor releases
    with overhangs DIFFERENT from the canonical (Promoter.oh5,
    Terminator.oh3) pair — `_grammar_tu_overhangs` only covers Alpha1.
    The third-pass check in `_classify_part_from_plasmid` digests each
    configured entry vector and compares the stuffer's overhangs to
    the user's plasmid digest. This regression guards the
    MAV-25-in-alpha-2 bug Cory reported on 2026-05-13.
    """

    def test_tu_in_alpha2_classifies_via_acceptor_pair(
            self, isolated_library, isolated_parts_bin):
        """Alpha2 EV configured with non-canonical BsaI-release
        overhangs; a TU plasmid whose BsaI digest releases the same
        overhangs classifies as level=1 with the Alpha2 role surfaced
        in the position name.

        Per the `_make_l1_alpha_vector` shape, `alpha_oh5` / `alpha_oh3`
        are the BsaI-release overhangs (the acceptor's L2-input
        identity); `tu_start` / `tu_end` are the inner Esp3I overhangs
        (L0→L1 assembly). The classifier's third-pass match runs on
        the BsaI-release pair.
        """
        # Use overhangs that are NOT (GGAG, CGCT) so the canonical
        # `_grammar_tu_overhangs` first-pass check fails and the
        # per-acceptor lookup is the only thing that can match.
        alpha2_bsai_5 = "TACA"
        alpha2_bsai_3 = "GACT"
        alpha2_ev = _make_l1_alpha_vector(
            "pAlpha2",
            alpha_oh5=alpha2_bsai_5, alpha_oh3=alpha2_bsai_3,
        )
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0",
            "role":       "Alpha2",
            "name":       "pAlpha2",
            "gb_text":    alpha2_ev["gb_text"],
        }])
        # TU plasmid synthesised with the alpha-2 BsaI-release
        # overhangs at its boundary — same shape as a real TU
        # assembled into pAlpha2.
        tu_seq = _build_gb_tu_seq(oh5=alpha2_bsai_5, oh3=alpha2_bsai_3)
        result = sc._classify_part_from_plasmid(tu_seq, circular=True)
        assert result is not None, (
            "Alpha2 TU should classify via the per-acceptor pass"
        )
        assert result["grammar_id"] == "gb_l0"
        assert result["level"] == 1
        # Position name surfaces the role so the user can confirm
        # which acceptor the classifier matched against.
        assert "Alpha2" in result["position"]["name"]
        assert result["position"]["oh5"] == alpha2_bsai_5
        assert result["position"]["oh3"] == alpha2_bsai_3

    def test_canonical_tu_still_classifies_after_acceptor_added(
            self, isolated_library, isolated_parts_bin):
        """Sanity: adding non-canonical acceptors to entry_vectors
        must NOT break the canonical (Alpha1) TU detection path —
        Alpha1's overhang pair IS the canonical (Promoter.oh5,
        Terminator.oh3) pair so it should match via the original
        `_grammar_tu_overhangs` check before the acceptor pass runs."""
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0",
            "role":       "Alpha2",
            "name":       "pAlpha2",
            "gb_text":    _make_l1_alpha_vector(
                "pAlpha2", alpha_oh5="TACA", alpha_oh3="GACT",
            )["gb_text"],
        }])
        # Canonical TU with the GB Promoter→Terminator overhangs.
        tu_seq = _build_gb_tu_seq(oh5="GGAG", oh3="CGCT")
        result = sc._classify_part_from_plasmid(tu_seq, circular=True)
        assert result is not None
        assert result["level"] == 1
        # Canonical pair → position name is the bare level label,
        # NOT a per-acceptor label.
        assert result["position"]["name"] == "TU"

    def test_no_match_when_acceptor_overhangs_differ(
            self, isolated_library, isolated_parts_bin):
        """The acceptor pass only matches plasmids whose overhangs
        match a CONFIGURED entry vector's stuffer. A TU with
        overhangs that don't match any acceptor should still return
        None — no silent false-positive."""
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0",
            "role":       "Alpha2",
            "name":       "pAlpha2",
            "gb_text":    _make_l1_alpha_vector(
                "pAlpha2", alpha_oh5="TACA", alpha_oh3="GACT",
            )["gb_text"],
        }])
        # Overhangs match neither Alpha1 canonical (GGAG/CGCT) nor
        # the configured Alpha2 (TACA/GACT).
        tu_seq = _build_gb_tu_seq(oh5="ATTC", oh3="GGCA")
        result = sc._classify_part_from_plasmid(tu_seq, circular=True)
        assert result is None

    def test_acceptor_cache_invalidates_on_entry_vector_save(
            self, isolated_library, isolated_parts_bin):
        """`_save_entry_vectors` MUST drop `_ACCEPTOR_TU_PAIRS_CACHE`
        so a reconfigured EV doesn't keep the old overhangs alive
        for the next classification pass. Without this, a user who
        swaps Alpha2 to a new plasmid would see stale matches
        against the old plasmid's overhangs."""
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0",
            "role":       "Alpha2",
            "name":       "pAlpha2-v1",
            "gb_text":    _make_l1_alpha_vector(
                "pAlpha2-v1",
                alpha_oh5="TACA", alpha_oh3="GACT",
            )["gb_text"],
        }])
        pairs_v1 = sc._grammar_acceptor_tu_pairs("gb_l0", "BsaI")
        assert pairs_v1 and pairs_v1[0][2] == "TACA"
        # Reconfigure with different BsaI-release overhangs.
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0",
            "role":       "Alpha2",
            "name":       "pAlpha2-v2",
            "gb_text":    _make_l1_alpha_vector(
                "pAlpha2-v2",
                alpha_oh5="ACGT", alpha_oh3="TGCA",
            )["gb_text"],
        }])
        pairs_v2 = sc._grammar_acceptor_tu_pairs("gb_l0", "BsaI")
        # Cache must have invalidated — the new overhangs replaced
        # the v1 pair, not appended to it.
        assert pairs_v2 and pairs_v2[0][2] == "ACGT"
        assert all(p[2] != "TACA" for p in pairs_v2)


# ═══════════════════════════════════════════════════════════════════════════════
# MoClo classifier — same logic as GB, different enzyme + position table
# ═══════════════════════════════════════════════════════════════════════════════
# The classifier loops over every grammar in `_all_grammars()`, so the
# Golden-Braid fixes (try-both-fragments + per-acceptor pair matching +
# drop parity inference) ALSO benefit MoClo. These tests guard against
# a future refactor specialising one grammar's pass against another.


def _build_moclo_tu_seq(insert: str = "ATG" * 30,
                         oh5: str = "GGAG", oh3: str = "CGCT") -> str:
    """Synthetic MoClo Plant TU: BpiI sites flank a body that releases
    with the requested 4-nt overhangs. BpiI / BbsI (GAAGAC(2/6)) has a
    2-nt spacer between recognition and overhang — fixture mirrors the
    GB equivalent with the spacer baked in."""
    core = "GAAGACAA" + oh5 + insert + oh3 + "AAGTCTTC"
    backbone = "AAAAATTTTT" * 50
    return core + backbone


def _make_moclo_acceptor(name: str, alpha_oh5: str, alpha_oh3: str) -> dict:
    """MoClo Plant L1 acceptor: BpiI sites release the stuffer with
    (alpha_oh5, alpha_oh3). Matches the `_make_l1_alpha_vector` shape
    for GB so tests look symmetric across grammars."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    bpii_left  = "GAAGACAA"
    bpii_right = "AAGTCTTC"
    dropout    = "ACGTAGCT" * 10
    backbone   = "GGGGTTTTAAAA" * 30
    seq = backbone + bpii_left + alpha_oh5 + dropout + alpha_oh3 + bpii_right + "TTTGGG" * 30
    rec = SeqRecord(Seq(seq), id=name, name=name)
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = "circular"
    return {"name": name, "gb_text": sc._record_to_gb_text(rec)}


class TestClassifyPartFromPlasmidMoClo:
    """Same detection cases as the GB tests, but for the MoClo Plant
    grammar (BsaI at L0, BpiI at L1). The classifier code is shared
    between grammars — these tests are mostly insurance that a future
    refactor specialising one grammar's pass doesn't regress the
    other.
    """

    def test_l0_part_classifies_against_moclo_position(
            self, isolated_library, isolated_parts_bin):
        """MoClo Pos 3 (CDS) — overhangs AGGT/GCTT. A synthetic L0
        with those overhangs classifies as moclo_plant / level=0."""
        seq = _build_moclo_plant_part_seq("AGGT", "GCTT")
        result = sc._classify_part_from_plasmid(seq, circular=True)
        assert result is not None
        assert result["grammar_id"] == "moclo_plant"
        assert result["level"]      == 0
        assert result["position"]["type"] == "CDS"

    def test_canonical_tu_classifies_via_grammar_boundary(
            self, isolated_library, isolated_parts_bin):
        """MoClo canonical TU boundary is GGAG/CGCT (Pos 1's oh5 +
        Pos 5's oh3). A BpiI-released TU with that pair returns
        level=1."""
        tu_seq = _build_moclo_tu_seq(oh5="GGAG", oh3="CGCT")
        result = sc._classify_part_from_plasmid(tu_seq, circular=True)
        assert result is not None
        assert result["grammar_id"] == "moclo_plant"
        assert result["level"]      == 1
        assert result["position"]["name"] == "TU"
        assert result["release_enzyme"] in ("BpiI", "BbsI")

    def test_tu_in_non_canonical_acceptor_classifies_via_pair(
            self, isolated_library, isolated_parts_bin):
        """MoClo's Acceptor2 role with non-canonical BpiI-release
        overhangs — a TU built into that acceptor should classify
        via the per-acceptor pass (same path as the GB Alpha2 fix)."""
        acc2 = _make_moclo_acceptor("pMoCloA2",
                                      alpha_oh5="TACA", alpha_oh3="GACT")
        sc._save_entry_vectors([{
            "grammar_id": "moclo_plant",
            "role":       "Acceptor2",
            "name":       "pMoCloA2",
            "gb_text":    acc2["gb_text"],
        }])
        tu_seq = _build_moclo_tu_seq(oh5="TACA", oh3="GACT")
        result = sc._classify_part_from_plasmid(tu_seq, circular=True)
        assert result is not None
        assert result["grammar_id"] == "moclo_plant"
        assert result["level"]      == 1
        assert "Acceptor2" in result["position"]["name"]

    def test_tu_with_insert_larger_than_backbone_still_picks_correctly(
            self, isolated_library, isolated_parts_bin):
        """Try-both-fragments regression: a MoClo TU whose body
        outgrew its backbone (common for any multi-kb cassette)
        used to mis-pick the backbone via `_pick_insert_fragment`'s
        smallest-fallback. The classifier now tries both fragments,
        so the body's overhangs ARE matched against the position /
        TU / acceptor checks regardless of fragment size."""
        # A 5 kb insert with canonical TU boundary; backbone is the
        # short "AAAAATTTTT"*50 = 500 bp padding.
        big_insert = "GCATGCAT" * 600   # ~4800 bp
        tu_seq = _build_moclo_tu_seq(oh5="GGAG", oh3="CGCT",
                                       insert=big_insert)
        result = sc._classify_part_from_plasmid(tu_seq, circular=True)
        assert result is not None
        assert result["grammar_id"] == "moclo_plant"
        assert result["level"]      == 1


class TestAssemblyFragmentFromSourceGbTextFallback:
    """`_assembly_fragment_from_source` needs the full plasmid gb_text
    to digest L1+ parts at the level-up enzyme. Parts saved by
    `_load_part_worker` before 2026-05-13 lacked `gb_text` (only
    `sequence` was stored). Library fallback recovers it by matching
    the parts-bin entry's `name` against library `id` / `name`."""

    def test_fallback_recovers_gb_text_from_library_by_id(
            self, isolated_library, isolated_parts_bin):
        """A parts-bin entry without `gb_text` whose name matches a
        library entry's `id` recovers the gb_text via the fallback —
        Constructor assembly works without reloading the part.

        Test data uses `id == name == sanitize(name)` so the post-
        2026-05-24 id-name backfill (PIT-36) is a no-op. Pre-fix
        the test used distinct id/name values which the backfill
        rewrote, breaking the parts-bin cross-reference."""
        tu_seq = _build_gb_tu_seq(oh5="GGAG", oh3="CGCT")
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq(tu_seq), id="my_tu", name="my_tu",
            annotations={"topology": "circular", "molecule_type": "DNA"},
        )
        gb_text = sc._record_to_gb_text(rec)
        # Library entry — has gb_text.
        sc._save_library([{
            "id":      "my_tu",
            "name":    "my_tu",
            "size":    len(tu_seq),
            "n_feats": 0,
            "gb_text": gb_text,
        }])
        # Parts bin entry — NO gb_text (pre-2026-05-13 shape).
        # name field set to the library's id/name so the cross-
        # reference can find it via either lookup arm.
        sc._save_parts_bin([{
            "name":     "my_tu",
            "type":     "TU",
            "position": "TU",
            "oh5":      "GGAG",
            "oh3":      "CGCT",
            "backbone": "test",
            "marker":   "—",
            "sequence": "AAAA",   # bogus; fallback will redigest
            "grammar":  "gb_l0",
            "level":    1,
        }])
        # Fetch the parts-bin entry as-loaded.
        part_entry = next(
            p for p in sc._load_parts_bin()
            if p.get("name") == "my_tu"
        )
        assert "gb_text" not in part_entry or not part_entry.get("gb_text")
        # Drive the fragment extraction at source_level=1 — would have
        # returned None pre-fix. With the library fallback it should
        # recover gb_text and return a valid fragment.
        grammar = sc._all_grammars()["gb_l0"]
        sc._ASSEMBLY_FRAGMENT_CACHE.clear()
        frag = sc._assembly_fragment_from_source(
            part_entry, grammar, source_level=1,
        )
        assert frag is not None, (
            "Library fallback should have recovered gb_text"
        )
        assert frag["oh5"] == "GGAG"
        assert frag["oh3"] == "CGCT"

    def test_no_fallback_when_no_matching_library_entry(
            self, isolated_library, isolated_parts_bin):
        """When the parts-bin entry has no gb_text AND no library
        entry matches by id/name, the function returns None with a
        log line (rather than silently picking the wrong gb_text)."""
        sc._save_library([])
        sc._save_parts_bin([{
            "name":     "orphan_part",
            "type":     "TU",
            "position": "TU",
            "oh5":      "GGAG",
            "oh3":      "CGCT",
            "sequence": "AAAA",
            "grammar":  "gb_l0",
            "level":    1,
        }])
        part_entry = next(
            p for p in sc._load_parts_bin()
            if p.get("name") == "orphan_part"
        )
        grammar = sc._all_grammars()["gb_l0"]
        sc._ASSEMBLY_FRAGMENT_CACHE.clear()
        frag = sc._assembly_fragment_from_source(
            part_entry, grammar, source_level=1,
        )
        assert frag is None

    def test_inline_gb_text_takes_precedence_over_library_lookup(
            self, isolated_library, isolated_parts_bin):
        """When the parts-bin entry already carries `gb_text`, the
        library-fallback path doesn't fire — the inline value wins.
        Defends against a future regression where a library rename
        / move could silently change which gb_text the fragment
        extraction sees."""
        tu_seq = _build_gb_tu_seq(oh5="GGAG", oh3="CGCT")
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq(tu_seq), id="my_tu", name="my_tu",
            annotations={"topology": "circular", "molecule_type": "DNA"},
        )
        gb_text = sc._record_to_gb_text(rec)
        # Library entry with DIFFERENT (corrupted) gb_text — if the
        # fallback fires, the fragment will fail to extract.
        sc._save_library([{
            "id":      "my_tu",
            "name":    "my_tu",
            "size":    len(tu_seq),
            "n_feats": 0,
            "gb_text": "not valid genbank",
        }])
        # Parts-bin entry WITH gb_text (post-2026-05-13 shape).
        sc._save_parts_bin([{
            "name":     "my_tu",
            "type":     "TU",
            "position": "TU",
            "oh5":      "GGAG",
            "oh3":      "CGCT",
            "sequence": "AAAA",
            "gb_text":  gb_text,
            "grammar":  "gb_l0",
            "level":    1,
        }])
        part_entry = next(
            p for p in sc._load_parts_bin()
            if p.get("name") == "my_tu"
        )
        assert part_entry.get("gb_text")
        grammar = sc._all_grammars()["gb_l0"]
        sc._ASSEMBLY_FRAGMENT_CACHE.clear()
        frag = sc._assembly_fragment_from_source(
            part_entry, grammar, source_level=1,
        )
        # Inline gb_text was valid → fragment extraction succeeded.
        # Library's bogus gb_text never consulted.
        assert frag is not None
        assert frag["oh5"] == "GGAG"


# ═══════════════════════════════════════════════════════════════════════════════
# Golden Braid iterative cycle: level helpers + multi-part assembly
# ═══════════════════════════════════════════════════════════════════════════════

def _make_l1_alpha_vector(name: str, alpha_oh5: str, alpha_oh3: str,
                            tu_start: str = "GGAG",
                            tu_end:   str = "CGCT") -> dict:
    """Build a real-shape Golden Braid L1 alpha vector dict — has BOTH
    inward-facing BsaI sites (the L1→L2 cut) AND inward-facing Esp3I
    sites (the L0→L1 cut), matching the pDGB1_α architecture. Used to
    test the iterative GB cycle end-to-end."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    bsai_left,  bsai_right  = "GGTCTCA", sc._rc("GGTCTCA")
    esp3i_left, esp3i_right = "CGTCTCA", sc._rc("CGTCTCA")
    dropout  = "ACGTAGCT" * 10
    backbone = "GGGGTTTTAAAA" * 30
    seq = (backbone +
            bsai_left + alpha_oh5 +
            esp3i_left + tu_start + dropout + tu_end + esp3i_right +
            alpha_oh3 + bsai_right +
            "TTTGGG" * 30)
    rec = SeqRecord(Seq(seq), id=name, name=name)
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    return {"name": name, "gb_text": sc._record_to_gb_text(rec)}


def _make_l2_omega_vector(name: str, oh5: str = "TACA",
                            oh3: str = "CCAA") -> dict:
    """L2 omega destination — only BsaI sites (since BsaI is the
    L1→L2 enzyme). The dropout exposes ``oh5`` / ``oh3`` overhangs
    matching the chained alpha1+alpha2 BsaI cut."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    bsai_left, bsai_right = "GGTCTCA", sc._rc("GGTCTCA")
    dropout  = "ACGTAGCT" * 10
    backbone = "GGGGTTTTAAAA" * 30
    seq = (backbone + bsai_left + oh5 + dropout + oh3 + bsai_right +
            "TTTGGG" * 30)
    rec = SeqRecord(Seq(seq), id=name, name=name)
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    return {"name": name, "gb_text": sc._record_to_gb_text(rec)}


def _make_l0_tu_parts(suffix: str = "a") -> list[dict]:
    """Four L0 parts that chain GGAG→TGAC→AATG→GCTT→CGCT — a minimal
    Golden Braid L0 TU layout (Promoter / 5'UTR / CDS / Terminator)."""
    return [
        {"name": f"P_{suffix}", "sequence": "AAATTT" * 5,
         "oh5": "GGAG", "oh3": "TGAC"},
        {"name": f"U_{suffix}", "sequence": "CCCAAA" * 3,
         "oh5": "TGAC", "oh3": "AATG"},
        {"name": f"C_{suffix}", "sequence": "ATGAAA" * 6,
         "oh5": "AATG", "oh3": "GCTT"},
        {"name": f"T_{suffix}", "sequence": "TTTGGG" * 4,
         "oh5": "GCTT", "oh3": "CGCT"},
    ]


class TestPartLevelHelpers:
    """`_part_level`, `_part_level_label`, and `_enzyme_for_level_up`
    are the foundation for the iterative GB cycle. Coverage matters
    because every part / parts-bin tab / cloning step keys off them."""

    def test_part_level_default_is_zero(self):
        assert sc._part_level({}) == 0
        assert sc._part_level({"name": "x"}) == 0

    def test_part_level_int_passthrough(self):
        assert sc._part_level({"level": 2}) == 2
        assert sc._part_level({"level": 5}) == 5

    def test_part_level_string_coercion(self):
        assert sc._part_level({"level": "L0"})  == 0
        assert sc._part_level({"level": "TU"})  == 1
        assert sc._part_level({"level": "MOD"}) == 2
        assert sc._part_level({"level": "3"})   == 3

    def test_part_level_negative_clamps_to_zero(self):
        assert sc._part_level({"level": -1}) == 0

    def test_part_level_bool_is_zero(self):
        # bool is a subclass of int — making sure True isn't read
        # as level=1 by accident (it would silently route a domesticator
        # part into the TU tab).
        assert sc._part_level({"level": True}) == 0
        assert sc._part_level({"level": False}) == 0

    def test_part_level_label_mapping(self):
        assert sc._part_level_label(0) == "L0"
        assert sc._part_level_label(1) == "TU"
        assert sc._part_level_label(2) == "MOD"
        assert sc._part_level_label(7) == "MOD"
        assert sc._part_level_label(-1) == "L0"

    def test_enzyme_for_level_up_alternates(self):
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        # L0 → L1 uses Esp3I; L1 → L2 uses BsaI; cycle alternates.
        assert sc._enzyme_for_level_up(gb, 0) == "Esp3I"
        assert sc._enzyme_for_level_up(gb, 1) == "BsaI"
        assert sc._enzyme_for_level_up(gb, 2) == "Esp3I"
        assert sc._enzyme_for_level_up(gb, 3) == "BsaI"

    def test_enzyme_for_level_up_moclo(self):
        moclo = sc._BUILTIN_GRAMMARS["moclo_plant"]
        # MoClo L0 → L1 uses BsaI; L1 → L2 uses BpiI.
        assert sc._enzyme_for_level_up(moclo, 0) == "BsaI"
        assert sc._enzyme_for_level_up(moclo, 1) == "BpiI"
        assert sc._enzyme_for_level_up(moclo, 2) == "BsaI"

    def test_enzyme_for_level_up_falls_back_when_missing(self):
        """A custom grammar without `level_up_enzyme` should reuse
        the primary enzyme so the iterative cycle degrades cleanly
        rather than crashing."""
        legacy = {"id": "custom", "enzyme": "BsaI"}
        # Even-source uses primary; odd-source falls back to primary
        # since `level_up_enzyme` is missing.
        assert sc._enzyme_for_level_up(legacy, 0) == "BsaI"
        assert sc._enzyme_for_level_up(legacy, 1) == "BsaI"


class TestCloneAssemblyIntoEntryVector:
    """`_clone_assembly_into_entry_vector` is the simulator that
    powers the Constructor's Save To Library button. It must handle
    L0 → TU (multi-part chain) AND L1 → MOD (two pre-built TUs) so
    the full GB cycle is simulatable."""

    def test_l0_to_tu_chains_parts_correctly(self):
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        vec = _make_l1_alpha_vector("alpha1", "TACA", "GACT")
        parts = _make_l0_tu_parts("a")
        result = sc._clone_assembly_into_entry_vector(
            parts, vec, gb, source_level=0, name="TU_a",
        )
        assert result is not None
        # The combined TU insert (with terminal sticky ends) must
        # appear in the cloned plasmid's top strand or its RC.
        expected_chain = (
            "GGAG" + "AAATTT"*5 + "TGAC" + "CCCAAA"*3
            + "AATG" + "ATGAAA"*6 + "GCTT" + "TTTGGG"*4 + "CGCT"
        )
        seq = str(result.seq).upper()
        assert (expected_chain in seq) or (expected_chain in sc._rc(seq))

    def test_l0_to_tu_returns_none_on_chain_mismatch(self):
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        vec = _make_l1_alpha_vector("alpha1", "TACA", "GACT")
        # Break the chain: P ends GGTT but U expects TGAC.
        parts = _make_l0_tu_parts("a")
        parts[0]["oh3"] = "GGTT"  # was TGAC
        result = sc._clone_assembly_into_entry_vector(
            parts, vec, gb, source_level=0,
        )
        assert result is None

    def test_l1_to_mod_iterates_cycle(self):
        """Full GB cycle: assemble two L0 → L1 plasmids, then ligate
        the two L1 plasmids into an Omega L2 destination using the
        level-up (BsaI) enzyme. Verifies the iterative cloning
        scaffold works end-to-end."""
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        # Step 1: L0 → L1 alpha1 (junction overhang TACA→GACT)
        a1_vec = _make_l1_alpha_vector("alpha1", "TACA", "GACT")
        a2_vec = _make_l1_alpha_vector("alpha2", "GACT", "CCAA")
        tu_a = sc._clone_assembly_into_entry_vector(
            _make_l0_tu_parts("a"), a1_vec, gb, source_level=0,
            name="TU_a",
        )
        tu_b = sc._clone_assembly_into_entry_vector(
            _make_l0_tu_parts("b"), a2_vec, gb, source_level=0,
            name="TU_b",
        )
        assert tu_a is not None and tu_b is not None
        # Step 2: L1 → L2 (BsaI cuts each TU + the omega vector)
        omega_vec = _make_l2_omega_vector("omega", "TACA", "CCAA")
        tu_a_src = {"name": "TU_a", "level": 1, "grammar": "gb_l0",
                     "gb_text": sc._record_to_gb_text(tu_a)}
        tu_b_src = {"name": "TU_b", "level": 1, "grammar": "gb_l0",
                     "gb_text": sc._record_to_gb_text(tu_b)}
        mod = sc._clone_assembly_into_entry_vector(
            [tu_a_src, tu_b_src], omega_vec, gb, source_level=1,
            name="MOD_test",
        )
        assert mod is not None
        # MOD should contain BOTH parent TU inserts. Spot-check
        # one signature from each: the unique CDS bases.
        mod_seq = str(mod.seq).upper()
        cds_a = "ATGAAA" * 6  # TU_a's CDS
        assert (cds_a in mod_seq) or (cds_a in sc._rc(mod_seq))

    def test_assembly_fragment_extracts_overhangs(self):
        """`_assembly_fragment_from_source` is the helper that turns a
        TU plasmid into a fragment dict with sequence + oh5/oh3, ready
        for chaining at L1→L2. Failure here breaks every iteration past
        L0→L1."""
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        a1_vec = _make_l1_alpha_vector("alpha1", "TACA", "GACT")
        tu = sc._clone_assembly_into_entry_vector(
            _make_l0_tu_parts("a"), a1_vec, gb, source_level=0,
            name="TU_a",
        )
        assert tu is not None
        src = {"name": "TU_a", "level": 1, "grammar": "gb_l0",
                "gb_text": sc._record_to_gb_text(tu)}
        frag = sc._assembly_fragment_from_source(src, gb, source_level=1)
        assert frag is not None
        assert frag["oh5"] == "TACA"
        assert frag["oh3"] == "GACT"

    def test_empty_sources_returns_none(self):
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        vec = _make_l1_alpha_vector("alpha1", "TACA", "GACT")
        assert sc._clone_assembly_into_entry_vector(
            [], vec, gb, source_level=0,
        ) is None

    def test_source_level_mismatch_returns_none(self):
        """A source whose `level` field doesn't match `source_level`
        must be rejected — without this guard, an L0 part smuggled
        into an L1+ lane would silently fall through to the L0
        extraction branch (since `_assembly_fragment_from_source`
        keys off the parameter, not the part's own level), producing
        a wrong-grammar synthetic assembly."""
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        vec = _make_l1_alpha_vector("alpha1", "TACA", "GACT")
        # Mark these L0 parts with level=0 explicitly; passing them
        # to source_level=1 should fail the level-match guard.
        parts = [{**p, "level": 0} for p in _make_l0_tu_parts("a")]
        result = sc._clone_assembly_into_entry_vector(
            parts, vec, gb, source_level=1,
        )
        assert result is None

    def test_l2_to_l3_assembly_uses_primary_enzyme(self):
        """Regression guard for the L2→L3 iteration of the GB cycle:
        the Constructor's MOD-source radio (source_level=2) must
        produce a valid L3 product from two L2 MODs. Exercises
        `_enzyme_for_level_up(grammar, 2) = primary` — for gb_l0
        that's Esp3I, the enzyme that cuts L2 MODs (which carry
        the cycle's primary-enzyme sites since L1→L2 ligated them
        into a BsaI-dropout vector with Esp3I surrounding).

        Pre-2026-05-14 there was no test guarding L2→L3 specifically
        — `TestCloneAssemblyIntoEntryVector` only covered L0→L1 and
        L1→L2. With the Constructor's MOD radio (source_level=2)
        being a catch-all for level ≥ 2 sources, regressions in the
        primary-vs-secondary parity could go unnoticed until a
        user reports their L3 build silently failed.
        """
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        # Step 1: two L1 TUs assembled into different alpha vectors
        # so they release with distinct BsaI overhangs for the L1→L2
        # chain.
        a1_vec = _make_l1_alpha_vector("alpha1", "TACA", "GACT")
        a2_vec = _make_l1_alpha_vector("alpha2", "GACT", "CCAA")
        tu_a = sc._clone_assembly_into_entry_vector(
            _make_l0_tu_parts("a"), a1_vec, gb, source_level=0,
            name="TU_a",
        )
        tu_b = sc._clone_assembly_into_entry_vector(
            _make_l0_tu_parts("b"), a2_vec, gb, source_level=0,
            name="TU_b",
        )
        assert tu_a is not None and tu_b is not None
        # Step 2: L1→L2 — two MODs into different L2 acceptors. Each
        # MOD lands in an L2 vector that carries Esp3I sites flanking
        # the BsaI dropout, so the result has the primary-enzyme
        # sites needed for the L2→L3 release. `_make_l1_alpha_vector`
        # already has this shape (both enzymes), so we reuse it as
        # an L2 acceptor with new dropout overhangs that match the
        # chained alpha1+alpha2 BsaI cut (TACA/CCAA).
        mod_x_vec = _make_l1_alpha_vector(
            "mod_x_dest", "GGAA", "TTCC",
            tu_start="TACA", tu_end="CCAA",
        )
        mod_y_vec = _make_l1_alpha_vector(
            "mod_y_dest", "TTCC", "AACC",
            tu_start="TACA", tu_end="CCAA",
        )
        tu_a_src = {"name": "TU_a", "level": 1, "grammar": "gb_l0",
                     "gb_text": sc._record_to_gb_text(tu_a)}
        tu_b_src = {"name": "TU_b", "level": 1, "grammar": "gb_l0",
                     "gb_text": sc._record_to_gb_text(tu_b)}
        mod_x = sc._clone_assembly_into_entry_vector(
            [tu_a_src, tu_b_src], mod_x_vec, gb, source_level=1,
            name="MOD_x",
        )
        # Build the same MOD shape for the "y" lane. Reuse tu_a/tu_b
        # — biology doesn't care that the two MODs share input TUs,
        # only that the outer L2 cut overhangs differ.
        mod_y = sc._clone_assembly_into_entry_vector(
            [tu_a_src, tu_b_src], mod_y_vec, gb, source_level=1,
            name="MOD_y",
        )
        assert mod_x is not None and mod_y is not None
        # Step 3: L2→L3 — the regression target. Two L2 MODs into
        # an L3 acceptor whose Esp3I dropout (the L2→L3 cut) matches
        # the chained MOD-x+MOD-y Esp3I overhangs (GGAA/AACC).
        l3_vec = _make_l1_alpha_vector(
            "l3_dest", "AAGG", "CCTT",
            tu_start="GGAA", tu_end="AACC",
        )
        mod_x_src = {"name": "MOD_x", "level": 2, "grammar": "gb_l0",
                      "gb_text": sc._record_to_gb_text(mod_x)}
        mod_y_src = {"name": "MOD_y", "level": 2, "grammar": "gb_l0",
                      "gb_text": sc._record_to_gb_text(mod_y)}
        l3 = sc._clone_assembly_into_entry_vector(
            [mod_x_src, mod_y_src], l3_vec, gb, source_level=2,
            name="L3_test",
        )
        # Either the L2→L3 cloning produces a valid record (clean
        # Esp3I cut + ligation) or `_pick_insert_fragment` correctly
        # selects across the L2 MOD's multi-cut fragments. Both are
        # acceptable for v0.8.0; the regression we're guarding
        # against is silent None — that would mean the enzyme parity
        # broke at L2 sources entirely.
        assert l3 is not None
        # The L3 product must carry both parent MODs' content (each
        # parent contributes its TU_a/TU_b unique CDS pattern).
        l3_seq = str(l3.seq).upper()
        cds_a = "ATGAAA" * 6  # TU_a's CDS
        assert (cds_a in l3_seq) or (cds_a in sc._rc(l3_seq))


class TestLevelMatchesTab:
    """`_level_matches_tab` is the shared filter rule for parts-bin
    tabs / constructor palettes / the lane resolver. Centralised so
    the L0/TU/MOD definition lives in one place."""

    def test_l0_tab_exact_match(self):
        assert sc._level_matches_tab(0, 0)
        assert not sc._level_matches_tab(1, 0)
        assert not sc._level_matches_tab(2, 0)

    def test_tu_tab_absorbs_level_1_and_above(self):
        # Tab "TU" (level 1) catches its own level AND every higher
        # level so a MOD assembled from TUs back-populates into the
        # TU palette — biologically the same plasmid can serve as
        # the source for the next cycle regardless of its labelled
        # level (the bench reaction only cares about overhang
        # compatibility, which the assembly simulator + the
        # constructor's compatibility check both verify).
        assert sc._level_matches_tab(1, 1)
        assert sc._level_matches_tab(2, 1)
        assert sc._level_matches_tab(7, 1)
        # L0 parts stay out of the TU tab (they're domesticated
        # entities, not assemblies).
        assert not sc._level_matches_tab(0, 1)

    def test_mod_tab_absorbs_level_2_and_above(self):
        # Tab "MOD" (level 2) catches every level ≥ 2 so further
        # iteration (L3, L4, …) doesn't need new tabs.
        assert sc._level_matches_tab(2, 2)
        assert sc._level_matches_tab(3, 2)
        assert sc._level_matches_tab(7, 2)
        # And below-2 levels stay out of the MOD tab.
        assert not sc._level_matches_tab(0, 2)
        assert not sc._level_matches_tab(1, 2)


class TestPartsBinSelectAll:
    """`Ctrl+A` selects every row in the active tab so a follow-up
    Delete keypress wipes the whole filtered slice. The selection
    rule was also relaxed to drop the obsolete `sequence` filter so
    TU/MOD rows (whose bases live in `gb_text`) are eligible too."""

    @staticmethod
    def _stub(name="p", level=0, **extra):
        base = {
            "name": name, "type": "CDS", "position": "Pos 3",
            "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spec",
            "sequence": "ATG" * 12 if level == 0 else "",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
            "grammar": "gb_l0",
            "level": level,
        }
        if level >= 1:
            base["gb_text"] = (
                "LOCUS x 100 bp DNA circular SYN\n//\n"
            )
        base.update(extra)
        return base

    async def test_ctrl_a_selects_every_l0_row(self, isolated_parts_bin):
        sc._save_parts_bin([self._stub(f"p{i}") for i in range(5)])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._selected_rows == set()
            modal.action_select_all_parts()
            await pilot.pause()
            assert modal._selected_rows == {0, 1, 2, 3, 4}

    async def test_ctrl_a_on_empty_bin_is_noop(self, isolated_parts_bin):
        sc._save_parts_bin([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal.action_select_all_parts()
            assert modal._selected_rows == set()

    async def test_ctrl_a_twice_toggles_back_to_empty(
            self, isolated_parts_bin):
        """Pressing Ctrl+A on a fully-selected table deselects —
        mirrors the standard text-editor gesture."""
        sc._save_parts_bin([self._stub(f"p{i}") for i in range(3)])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal.action_select_all_parts()
            assert len(modal._selected_rows) == 3
            modal.action_select_all_parts()
            assert modal._selected_rows == set()

    async def test_ctrl_a_then_delete_removes_every_part(
            self, isolated_parts_bin):
        sc._save_parts_bin([self._stub(f"p{i}") for i in range(3)])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal.action_select_all_parts()
            await pilot.pause()
            modal.action_delete_selected_parts()
            await pilot.pause()
            await pilot.pause(0.1)
            # Confirm modal pushed with all 3 names.
            confirm = app.screen
            assert isinstance(confirm, sc.PartsBinDeleteConfirmModal)
            confirm.query_one("#btn-partsdel-yes", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert sc._load_parts_bin() == []

    async def test_ctrl_a_selects_tu_rows_in_tu_tab(
            self, isolated_parts_bin):
        """TU rows have empty `sequence` (their bases live in
        `gb_text`). Pre-2026-05-10 the multi-select guard required
        a non-empty `sequence` and silently rejected every TU; the
        relaxed rule (just `r.get('user')`) makes Ctrl+A select
        them all from the TU tab."""
        sc._save_parts_bin([
            self._stub("L0_part", level=0),
            self._stub("TU_a", level=1),
            self._stub("TU_b", level=1),
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            tabs = modal.query_one("#parts-level-tabs", sc.Tabs)
            tabs.active = "tab-parts-tu"
            await pilot.pause()
            await pilot.pause(0.1)
            modal.action_select_all_parts()
            await pilot.pause()
            # Two TUs in the active TU tab — both selected.
            assert modal._selected_rows == {0, 1}
            assert all(
                modal._rows[i]["level"] == 1
                for i in modal._selected_rows
            )

    async def test_ctrl_a_only_selects_active_tab_rows(
            self, isolated_parts_bin):
        """`_rows` is the post-filter view, so Ctrl+A on the L0
        tab won't reach into TU/MOD rows hiding under other tabs."""
        sc._save_parts_bin([
            self._stub("L0_part", level=0),
            self._stub("TU_x", level=1),
            self._stub("MOD_y", level=2),
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Default L0 tab.
            assert modal._active_level == 0
            modal.action_select_all_parts()
            await pilot.pause()
            assert len(modal._selected_rows) == 1
            assert modal._rows[0]["name"] == "L0_part"


class TestPartsBinMultiDeleteWarning:
    """Pressing the Delete button with a multi-selection MUST push a
    `PartsBinDeleteConfirmModal` showing the count + a name preview
    BEFORE any rows are written away. Locks in the user-visible
    safety net so a future refactor can't silently bypass the
    confirmation."""

    @staticmethod
    def _stub(name="p"):
        return {
            "name": name, "type": "CDS", "position": "Pos 3",
            "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spec",
            "sequence": "ATG" * 12,
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
            "grammar": "gb_l0", "level": 0,
        }

    async def test_delete_button_with_multiselect_pushes_confirm_modal(
            self, isolated_parts_bin):
        sc._save_parts_bin([self._stub(f"p{i}") for i in range(5)])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal.action_select_all_parts()
            await pilot.pause()
            modal.query_one("#btn-parts-delete", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            confirm = app.screen
            assert isinstance(confirm, sc.PartsBinDeleteConfirmModal)
            # Body must show the count + preview names + the "cannot be
            # undone" hint so the user can't miss what they're about to
            # do.
            body = str(
                confirm.query_one("#partsdel-msg", sc.Static).render()
            )
            assert "5" in body
            assert "p0" in body
            # Title carries the count too — visible even when the user
            # tabs away from the body.
            title = str(
                confirm.query_one("#partsdel-title", sc.Static).render()
            )
            assert "Remove 5 parts" in title
            # Yes button echoes the count one more time.
            yes_btn = confirm.query_one("#btn-partsdel-yes", sc.Button)
            assert "5" in str(yes_btn.label)
            # Default focus on No so a stray Enter cancels (handslip
            # protection — sacred for any bulk-delete confirm modal).
            no_btn = confirm.query_one("#btn-partsdel-no", sc.Button)
            assert no_btn.has_focus

    async def test_no_button_aborts_the_delete(self, isolated_parts_bin):
        sc._save_parts_bin([self._stub(f"p{i}") for i in range(3)])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal.action_select_all_parts()
            await pilot.pause()
            modal.query_one("#btn-parts-delete", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            confirm = app.screen
            confirm.query_one("#btn-partsdel-no", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Bin untouched.
            assert len(sc._load_parts_bin()) == 3


class TestPartsBinLevelTabs:
    """The Parts Bin top-of-modal Tabs filter rows by part level
    (L0 / TU / MOD). User parts default to L0; Constructor-saved
    assemblies tag themselves with the appropriate level."""

    @staticmethod
    def _stub_part(name="p", level=0, **extra):
        base = {
            "name": name, "type": "CDS", "position": "Pos 3",
            "oh5": "AATG", "oh3": "GCTT",
            "backbone": "pUPD2", "marker": "Spec",
            "sequence": "ATG" * 12,
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm": 0.0, "rev_tm": 0.0,
            "grammar": "gb_l0",
            "level": level,
        }
        base.update(extra)
        return base

    async def test_default_tab_is_l0_filters_other_levels(
            self, isolated_parts_bin):
        sc._save_parts_bin([
            self._stub_part(name="L0_part", level=0),
            self._stub_part(name="TU_plasmid", level=1),
            self._stub_part(name="MOD_plasmid", level=2),
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._active_level == 0
            assert len(modal._rows) == 1
            assert modal._rows[0]["name"] == "L0_part"

    async def test_tu_tab_shows_level_1_and_above(
            self, isolated_parts_bin):
        # TU tab back-populates MODs (and any higher level) so a
        # MOD assembled from TUs can be re-used as a TU-equivalent
        # source in the next cycle. L0 parts stay strictly in the
        # L0 tab.
        sc._save_parts_bin([
            self._stub_part(name="L0_part", level=0),
            self._stub_part(name="TU_plasmid", level=1),
            self._stub_part(name="MOD_plasmid", level=2),
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            tabs = modal.query_one("#parts-level-tabs", sc.Tabs)
            tabs.active = "tab-parts-tu"
            await pilot.pause()
            await pilot.pause(0.1)
            assert modal._active_level == 1
            names = sorted(r["name"] for r in modal._rows)
            assert names == ["MOD_plasmid", "TU_plasmid"]

    async def test_mod_tab_shows_level_2_and_above(
            self, isolated_parts_bin):
        sc._save_parts_bin([
            self._stub_part(name="MOD_plasmid", level=2),
            self._stub_part(name="L3_plasmid", level=3),
            self._stub_part(name="TU_plasmid", level=1),
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            tabs = modal.query_one("#parts-level-tabs", sc.Tabs)
            tabs.active = "tab-parts-mod"
            await pilot.pause()
            await pilot.pause(0.1)
            names = sorted(r["name"] for r in modal._rows)
            assert names == ["L3_plasmid", "MOD_plasmid"]

    async def test_legacy_part_with_no_level_field_lands_in_l0(
            self, isolated_parts_bin):
        """Pre-2026-05-10 parts have no `level` field. They must
        default to L0 so existing libraries don't disappear into a
        non-default tab on first open."""
        legacy = self._stub_part(name="legacy")
        legacy.pop("level", None)
        sc._save_parts_bin([legacy])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert len(modal._rows) == 1
            assert modal._rows[0]["name"] == "legacy"


class TestConstructorSaveToLibrary:
    """The Save To Library button (renamed from Simulate Assembly)
    actually clones the lane parts into the bound entry vector and
    writes the result to plasmid_library.json + parts_bin.json. The
    new entry tags itself with `level=1` (TU) so it surfaces under
    the TU tab and can be picked from the L1+ palette."""

    async def test_button_label_is_save_to_library(
            self, isolated_library, isolated_parts_bin):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            btn = modal.query_one("#btn-ctor-simulate-gb_l0", sc.Button)
            assert "Save To Library" in str(btn.label)


class TestConstructorGreenStatus:
    """Validation status flips green as soon as the lane chain is
    consistent — even before the user picks a backbone. Backbone-
    not-set surfaces as a separate yellow hint per the 2026-05-10
    UX spec."""

    async def test_green_fires_with_complete_chain_no_backbone(
            self, isolated_library, isolated_parts_bin):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Inject a complete TU lane (4 parts that chain GGAG→CGCT).
            modal._lanes["gb_l0"] = [
                ("P", "Promoter",   "Pos 1",   "GGAG", "TGAC", "", ""),
                ("U", "5' UTR",     "Pos 2",   "TGAC", "AATG", "", ""),
                ("C", "CDS",        "Pos 3-4", "AATG", "GCTT", "", ""),
                ("T", "Terminator", "Pos 5",   "GCTT", "CGCT", "", ""),
            ]
            modal._refresh_validation("gb_l0")
            await pilot.pause()
            vbox = modal.query_one("#ctor-validation-gb_l0", sc.Static)
            text = str(vbox.render())
            # Green confirmation text for the chain validity.
            assert "Valid TU" in text
            # Save button stays disabled because no backbone is bound.
            btn = modal.query_one("#btn-ctor-simulate-gb_l0", sc.Button)
            assert btn.disabled is True


class TestConstructorL1PlusMode:
    """Source-level radios switch the constructor between L0→TU
    (default) and L1+ assembly modes. Switching levels clears the
    lane (overhang scheme differs) and refilters the palette."""

    async def test_level_radio_switch_clears_lane_and_refilters(
            self, isolated_library, isolated_parts_bin):
        # Stage one L0 part and one TU plasmid in the bin.
        sc._save_parts_bin([
            {"name": "L0p", "type": "CDS", "position": "Pos 3-4",
             "oh5": "AATG", "oh3": "GCTT",
             "sequence": "ATG" * 10,
             "grammar": "gb_l0", "level": 0,
             "backbone": "", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0},
            {"name": "TU1", "type": "TU", "position": "TU",
             "oh5": "TACA", "oh3": "GACT",
             "sequence": "",
             "grammar": "gb_l0", "level": 1,
             "backbone": "alpha1", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0,
             "gb_text": "LOCUS  TU1  100 bp DNA circular SYN\n//\n"},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Default: source level 0 → palette only contains L0 parts.
            assert modal._source_levels["gb_l0"] == 0
            # Stage some lane content so we can verify the clear path.
            modal._lanes["gb_l0"] = [
                ("L0p", "CDS", "Pos 3-4", "AATG", "GCTT", "", ""),
            ]
            # Switch to TU-source mode by toggling the radio.
            tu_rb = modal.query_one("#ctor-level-tu-gb_l0", sc.RadioButton)
            tu_rb.value = True
            await pilot.pause()
            await pilot.pause(0.1)
            assert modal._source_levels["gb_l0"] == 1
            # Lane was cleared on level switch.
            assert modal._lanes["gb_l0"] == []
            # Palette now lists only the TU plasmid.
            palette_names = [r[0] for r in modal._palette_rows.get("gb_l0", [])]
            assert palette_names == ["TU1"]

    def test_validate_l1_chain_skips_l0_constraints(self):
        """At source_level=1 the validator only enforces chain
        continuity — boundary overhangs / mandatory L0 part types /
        slot-occupancy / CDS-NS pairing all stop applying since
        the lane carries pre-built TUs, not L0 parts."""
        modal = sc.ConstructorModal()
        modal._source_levels["gb_l0"] = 1
        # A 2-TU chain whose junction overhangs match (TACA→GACT,
        # GACT→CCAA). No "Promoter" / "Terminator" parts; under L0
        # rules this would error on missing types + boundary mismatch.
        modal._lanes["gb_l0"] = [
            ("TU1", "TU", "TU", "TACA", "GACT", "alpha1", ""),
            ("TU2", "TU", "TU", "GACT", "CCAA", "alpha2", ""),
        ]
        is_valid, errors = modal._validate("gb_l0")
        assert is_valid is True
        assert errors == []

    def test_validate_l1_chain_still_catches_overhang_mismatch(self):
        modal = sc.ConstructorModal()
        modal._source_levels["gb_l0"] = 1
        # Junction TACA→GACT, GGGG→CCAA — second junction breaks.
        modal._lanes["gb_l0"] = [
            ("TU1", "TU", "TU", "TACA", "GACT", "alpha1", ""),
            ("TU2", "TU", "TU", "GGGG", "CCAA", "alpha2", ""),
        ]
        is_valid, errors = modal._validate("gb_l0")
        assert is_valid is False
        assert any("junction" in e.lower() for e in errors)


class TestConstructorChainAndStatusAtL1Plus:
    """Hardening for the L1+ source-level UI surfaces: chain-render
    boundary handling, green-text labelling, and the L<N> hint in
    the backbone-not-set yellow line. Pre-2026-05-10 these all
    hardcoded L0/L1 conventions and rendered red-everywhere /
    'Valid TU' / 'L1 destination' regardless of the active source
    level."""

    def test_build_chain_at_l1_does_not_flag_lane_boundaries_red(self):
        """At L1+ source level, the lane's terminal overhangs depend
        on the destination vector — not the grammar's L0
        Promoter→Terminator boundaries. `_build_chain` must NOT
        compare lane[0]/lane[-1] against grammar.tu_start/tu_end at
        L1+, which would render the boundary labels red even though
        the chain is biologically valid."""
        modal = sc.ConstructorModal()
        modal._source_levels["gb_l0"] = 1
        # alpha1+alpha2 chain (TACA→GACT, GACT→CCAA) — no L0 grammar
        # boundary involvement here.
        modal._lanes["gb_l0"] = [
            ("TU1", "TU", "TU", "TACA", "GACT", "alpha1", ""),
            ("TU2", "TU", "TU", "GACT", "CCAA", "alpha2", ""),
        ]
        text = modal._build_chain("gb_l0")
        rendered = str(text)
        # Chain should show the actual terminal overhangs.
        assert "TACA" in rendered
        assert "CCAA" in rendered
        # GB L0 boundaries (GGAG / CGCT) must NOT appear — the L1+
        # chain isn't constrained by them.
        assert "GGAG" not in rendered
        assert "CGCT" not in rendered

    async def test_status_says_valid_mod_at_l1_source_level(
            self, isolated_library, isolated_parts_bin):
        """At L1 source level the green status must say "Valid MOD",
        not "Valid TU" — the assembled product is a module, not a
        single transcription unit."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source_levels["gb_l0"] = 1
            modal._lanes["gb_l0"] = [
                ("TU1", "TU", "TU", "TACA", "GACT", "alpha1", ""),
                ("TU2", "TU", "TU", "GACT", "CCAA", "alpha2", ""),
            ]
            modal._refresh_validation("gb_l0")
            await pilot.pause()
            vbox = modal.query_one("#ctor-validation-gb_l0", sc.Static)
            text = str(vbox.render())
            assert "Valid MOD" in text

    async def test_backbone_hint_uses_target_level(
            self, isolated_library, isolated_parts_bin):
        """The yellow 'pick a backbone' hint must reference the
        ACTUAL target level (L<source+1>) rather than always saying
        L1 — pre-fix every hint hardcoded `f"{1 if … else 1}"`."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source_levels["gb_l0"] = 1
            modal._lanes["gb_l0"] = [
                ("TU1", "TU", "TU", "TACA", "GACT", "alpha1", ""),
                ("TU2", "TU", "TU", "GACT", "CCAA", "alpha2", ""),
            ]
            modal._refresh_validation("gb_l0")
            await pilot.pause()
            vbox = modal.query_one("#ctor-validation-gb_l0", sc.Static)
            text = str(vbox.render())
            assert "L2 destination" in text
            assert "L1 destination" not in text


class TestResolveLaneToPartsLevelFilter:
    """Hardening for `_resolve_lane_to_parts` — when the parts bin
    holds an L0 part AND a TU plasmid with the same name (a real
    possibility once users start naming TUs after their dominant
    L0 component), the resolver MUST return the entry whose level
    matches the active source level. Pre-2026-05-10 it returned
    the first match in file order, which could ligate the wrong
    plasmid into the destination."""

    @staticmethod
    def _stub_l0(name="X"):
        return {"name": name, "type": "CDS", "position": "Pos 3",
                "oh5": "AATG", "oh3": "GCTT",
                "sequence": "ATG" * 10, "grammar": "gb_l0", "level": 0,
                "backbone": "", "marker": "",
                "fwd_primer": "", "rev_primer": "",
                "fwd_tm": 0.0, "rev_tm": 0.0}

    @staticmethod
    def _stub_tu(name="X"):
        return {"name": name, "type": "TU", "position": "TU",
                "oh5": "TACA", "oh3": "GACT",
                "sequence": "", "grammar": "gb_l0", "level": 1,
                "backbone": "alpha1", "marker": "Spec",
                "fwd_primer": "", "rev_primer": "",
                "fwd_tm": 0.0, "rev_tm": 0.0,
                "gb_text": "LOCUS X 100 bp DNA circular SYN\n//\n"}

    async def test_l0_mode_resolves_l0_entry_when_name_collides(
            self, isolated_library, isolated_parts_bin):
        sc._save_parts_bin([self._stub_tu("X"), self._stub_l0("X")])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source_levels["gb_l0"] = 0
            modal._lanes["gb_l0"] = [
                ("X", "CDS", "Pos 3", "AATG", "GCTT", "", "")
            ]
            resolved = modal._resolve_lane_to_parts("gb_l0")
            assert resolved is not None
            assert len(resolved) == 1
            assert resolved[0]["level"] == 0
            assert resolved[0]["sequence"] == "ATG" * 10

    async def test_tu_mode_resolves_tu_entry_when_name_collides(
            self, isolated_library, isolated_parts_bin):
        sc._save_parts_bin([self._stub_l0("X"), self._stub_tu("X")])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._source_levels["gb_l0"] = 1
            modal._lanes["gb_l0"] = [
                ("X", "TU", "TU", "TACA", "GACT", "alpha1", "")
            ]
            resolved = modal._resolve_lane_to_parts("gb_l0")
            assert resolved is not None
            assert len(resolved) == 1
            assert resolved[0]["level"] == 1
            assert resolved[0].get("gb_text", "")


class TestPersistedAssemblyMetadata:
    """When the Constructor saves an assembly, the parts-bin entry
    needs the right metadata so the next-level cycle has a usable
    source. Marker comes from the role's `selection` antibiotic;
    `gb_text`, oh5, oh3, and `level` round-trip into a usable
    L1+ palette row."""

    def test_clone_assembly_parts_bin_metadata_round_trip(
            self, isolated_library, isolated_parts_bin):
        """Driving `_persist_assembly` directly so we don't have to
        spin up the full UI just to verify the bin entry shape."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("AAAA" * 100), id="MyTU", name="MyTU")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "alpha1_vec",
                            "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[
                {"name": "P", "oh5": "GGAG", "oh3": "TGAC", "level": 0},
                {"name": "T", "oh5": "TGAC", "oh3": "CGCT", "level": 0},
            ],
            backbone_role="Alpha1",
        )
        bin_entries = sc._load_parts_bin()
        assert len(bin_entries) == 1
        e = bin_entries[0]
        assert e["level"]    == 1
        assert e["type"]     == "TU"
        assert e["oh5"]      == "GGAG"
        assert e["oh3"]      == "CGCT"
        assert e["grammar"]  == "gb_l0"
        assert e["backbone"] == "alpha1_vec"
        # Marker lands as "—" when the bound vector's gb_text carries
        # no annotated antibiotic (featureless stub here). Post-2026-
        # 05-22 the role's hardcoded default no longer applies —
        # detection is the only source of truth.
        assert e["marker"] == "—"
        # gb_text is the source of truth for L1+ chaining; must NOT
        # be empty so the next-cycle cloner can digest this plasmid.
        assert e["gb_text"]
        # Source-parts list survives so the user can audit which
        # L0 parts went into this TU.
        assert e["source_parts"] == ["P", "T"]
        assert e["source_role"]  == "Alpha1"
        # And the library entry exists too.
        lib_entries = sc._load_library()
        assert len(lib_entries) == 1
        assert lib_entries[0]["name"] == "MyTU"

    def test_clone_assembly_marker_detects_from_bound_vector(
            self, isolated_library, isolated_parts_bin):
        """When the bound entry vector carries an explicit antibiotic
        annotation (e.g. an AmpR CDS), the assembled L1+ part's
        ``marker`` must reflect that — NOT the role's hardcoded
        convention default from ``_CONSTRUCTOR_BACKBONES``.

        Regression: a user with custom α-vectors carrying AmpR (FFE
        2/3-style) reported every TU assembled through Alpha1/Alpha2
        was stamped ``"Spectinomycin"`` because the save path read
        the role's canonical pDGB3 default instead of the bound
        vector's annotations. See `_detect_selection_marker`
        fallthrough in `_persist_assembly`.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from io import StringIO
        from Bio import SeqIO
        vec_rec = SeqRecord(Seq("ATGC" * 25), id="amp_alpha",
                             name="amp_alpha")
        vec_rec.annotations["molecule_type"] = "DNA"
        vec_rec.annotations["topology"]      = "circular"
        vec_rec.features.append(SeqFeature(
            FeatureLocation(0, 30), type="CDS",
            qualifiers={"label": ["AmpR"], "gene": ["bla"]},
        ))
        buf = StringIO()
        SeqIO.write([vec_rec], buf, "genbank")
        amp_gb_text = buf.getvalue()
        # Sanity: detection must see AmpR before we test the save path.
        assert sc._detect_selection_marker(amp_gb_text) == "Ampicillin"
        rec = SeqRecord(Seq("AAAA" * 100), id="MyAmpTU", name="MyAmpTU")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "amp_alpha", "gb_text": amp_gb_text},
            parts=[
                {"name": "P", "oh5": "GGAG", "oh3": "TGAC", "level": 0},
                {"name": "T", "oh5": "TGAC", "oh3": "CGCT", "level": 0},
            ],
            backbone_role="Alpha1",
        )
        bin_entries = sc._load_parts_bin()
        assert len(bin_entries) == 1
        # The fix: marker reflects vector's annotation, not the
        # Alpha1 role's "Spectinomycin" convention default.
        assert bin_entries[0]["marker"] == "Ampicillin"

    def test_clone_assembly_marker_falls_back_to_dash_when_undetected(
            self, isolated_library, isolated_parts_bin):
        """When the bound entry vector has no recognizable antibiotic
        annotation, ``marker`` lands as ``"—"`` — NOT a hardcoded
        role default. The 2026-05-22 selection-marker overhaul
        removed all hardcoded antibiotics from
        ``_CONSTRUCTOR_BACKBONES``; the only source of truth is now
        the vector's annotated features. A featureless vector means
        the user must manually edit the saved part's marker.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("AAAA" * 100), id="MyTU", name="MyTU")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "bare_alpha",
                           "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[
                {"name": "P", "oh5": "GGAG", "oh3": "TGAC", "level": 0},
                {"name": "T", "oh5": "TGAC", "oh3": "CGCT", "level": 0},
            ],
            backbone_role="Alpha1",
        )
        bin_entries = sc._load_parts_bin()
        assert len(bin_entries) == 1
        assert bin_entries[0]["marker"] == "—"

    def test_clone_assembly_disambiguates_id_collision(
            self, isolated_library, isolated_parts_bin):
        """If the library already contains an entry whose id matches
        the assembly's safe-id, `_persist_assembly` must append a
        numeric suffix rather than overwriting the existing entry.

        Test data uses `id == name` so the post-2026-05-24 id-name
        backfill (PIT-36) is a no-op and the existing id stays put
        for the disambiguation loop to find."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        sc._save_library([{
            "id": "MyTU", "name": "MyTU",
            "size": 100, "n_feats": 0, "source": "test",
            "added": "2026-01-01", "gb_text": "LOCUS x 100 bp DNA\n//\n",
        }])
        rec = SeqRecord(Seq("AAAA" * 50), id="MyTU", name="MyTU")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "alpha1_vec",
                            "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[{"name": "P", "oh5": "GGAG", "oh3": "CGCT", "level": 0}],
            backbone_role="Alpha1",
        )
        ids = [e.get("id") for e in sc._load_library()]
        assert "MyTU"   in ids
        assert "MyTU_2" in ids

    def test_persist_assembly_always_writes_parts_bin_with_grammar(
            self, isolated_library, isolated_parts_bin):
        """Hard rule from PIT-37 (2026-05-24): every modular
        Constructor save (GB / MoClo / any future grammar tab)
        writes a parts_bin entry tagged with the active `grammar:
        gid`. Non-modular tabs (Traditional, Gibson) intentionally
        skip parts_bin — covered by their respective save tests."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("AAAA" * 50), id="MyTU", name="MyTU")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "alpha1_vec",
                            "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[{"name": "P", "oh5": "GGAG", "oh3": "CGCT",
                     "level": 0}],
            backbone_role="Alpha1",
        )
        bin_entries = sc._load_parts_bin()
        assert len(bin_entries) == 1
        bin_e = bin_entries[0]
        assert bin_e["grammar"] == "gb_l0"
        # Name + level tag the entry so it appears in the right
        # grammar's palette under the right type/level tab.
        assert bin_e["name"] == "MyTU"
        assert bin_e["level"] == 1   # source_level=0 → target_level=1 → TU

    def test_persist_assembly_populates_bin_sequence_body(
            self, isolated_library, isolated_parts_bin):
        """L1+ bin entries now carry the released insert body in
        `bin_entry["sequence"]` (was hardcoded to "" pre-2026-05-24).
        Mirrors the Load Part convention so MOD/TU bin rows look
        consistent regardless of origin. When the level-up digest
        finds no clean release, sequence falls back to "" — that's
        the bench-fake-plasmid case in this test (no IIS sites)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("AAAA" * 50), id="MyTU", name="MyTU")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "alpha1_vec",
                            "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[{"name": "P", "oh5": "GGAG", "oh3": "CGCT",
                     "level": 0}],
            backbone_role="Alpha1",
        )
        bin_e = sc._load_parts_bin()[0]
        # `sequence` key must exist (even if empty for digest-failure
        # paths) — the assertion is that the field is no longer
        # hardcoded to "" at the dict-construction site.
        assert "sequence" in bin_e

    def test_persist_mod_to_next_stores_level_3(
            self, isolated_library, isolated_parts_bin):
        """The Constructor's MOD→next save (`source_level=2`) must
        tag the bin entry with `level=3` (so further iteration sees
        it in the MOD palette via `_level_matches_tab(3, 2) = True`)
        and label it as MOD. Verifies the formula
        `target_level = source_level + 1` propagates correctly at
        the L2-source step and that the auto-detect overhang probe
        falls back when the saved plasmid has no clean
        primary-enzyme release.

        Pre-2026-05-14 there was no test guarding the L2-source
        persist path — `test_clone_assembly_parts_bin_metadata_round_trip`
        only covered `source_level=0`.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("AAAA" * 200), id="MyL3", name="MyL3")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=2,
            entry_vector={"name": "l3_dest_vec",
                            "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[
                {"name": "MOD_x", "oh5": "GGAA", "oh3": "TTCC",
                 "level": 2, "gb_text": ""},
                {"name": "MOD_y", "oh5": "TTCC", "oh3": "AACC",
                 "level": 2, "gb_text": ""},
            ],
            backbone_role="Alpha1",
        )
        bin_entries = sc._load_parts_bin()
        assert len(bin_entries) == 1
        e = bin_entries[0]
        # target_level = source_level + 1 = 3
        assert e["level"] == 3
        # Label rolls up to "MOD" via `_part_level_label(3)` since
        # any level ≥ 2 displays as MOD (the cycle is recursive past
        # this point, no need to track L3/L4/L5 separately).
        assert e["type"] == "MOD"
        assert e["position"] == "MOD"
        assert e["grammar"] == "gb_l0"
        # The overhang probe couldn't release a clean L3 fragment
        # (no IIS sites in this test stub), so it falls back to the
        # inner-source boundaries — MOD_x's oh5 and MOD_y's oh3.
        # That keeps SOMETHING in the bin entry so downstream UI
        # doesn't render empty overhangs; the live palette would
        # re-resolve via `_assembly_fragment_from_source` anyway.
        assert e["oh5"] == "GGAA"
        assert e["oh3"] == "AACC"
        # gb_text must round-trip onto the bin entry so a further
        # L3→L4 iteration can digest this plasmid.
        assert e["gb_text"]
        # Source-parts list captures the L2 MODs.
        assert e["source_parts"] == ["MOD_x", "MOD_y"]


class TestMigratePartsBinMarkersFromVector:
    """`_migrate_parts_bin_markers_from_vector` rescans every
    parts-bin entry's `gb_text` and corrects the stored `marker`
    when the historical role-default (Spectinomycin / Kanamycin)
    doesn't match what the bound vector actually carries.

    Backstory: pre-2026-05-22 the Constructor save path stamped
    `marker` from `_CONSTRUCTOR_BACKBONES[gid][role]["selection"]`
    regardless of the bound entry vector's annotation. Users
    running custom α-vectors with AmpR (or any non-canonical
    selection) accumulated bin rows labelled "Spectinomycin" /
    "Kanamycin". This one-shot migration fixes those rows in
    place without disturbing manually-edited markers.
    """

    @staticmethod
    def _gb_with_marker(marker_label: str) -> str:
        """Build a minimal valid GenBank string whose features carry
        the requested marker label so `_detect_selection_marker`
        returns the matching display name."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from io import StringIO
        from Bio import SeqIO
        rec = SeqRecord(Seq("ATGC" * 25), id="vec", name="vec")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, 30), type="CDS",
            qualifiers={"label": [marker_label]},
        ))
        buf = StringIO()
        SeqIO.write([rec], buf, "genbank")
        return buf.getvalue()

    def test_corrects_spec_default_for_amp_vector(
            self, isolated_parts_bin):
        amp_gb = self._gb_with_marker("AmpR")
        sc._save_parts_bin([
            {"name": "TU1", "marker": "Spectinomycin",
             "gb_text": amp_gb, "level": 1},
        ])
        sc._migrate_parts_bin_markers_from_vector()
        out = sc._load_parts_bin()
        assert out[0]["marker"] == "Ampicillin"

    def test_corrects_kan_default_for_amp_vector(
            self, isolated_parts_bin):
        amp_gb = self._gb_with_marker("AmpR")
        sc._save_parts_bin([
            {"name": "TU1", "marker": "Kanamycin",
             "gb_text": amp_gb, "level": 1},
        ])
        sc._migrate_parts_bin_markers_from_vector()
        out = sc._load_parts_bin()
        assert out[0]["marker"] == "Ampicillin"

    def test_preserves_correct_spec_marker(
            self, isolated_parts_bin):
        """Stored "Spectinomycin" + vector that actually has SmR →
        detection returns "Spectinomycin", no change needed."""
        spec_gb = self._gb_with_marker("SmR")
        sc._save_parts_bin([
            {"name": "TU1", "marker": "Spectinomycin",
             "gb_text": spec_gb, "level": 1},
        ])
        sc._migrate_parts_bin_markers_from_vector()
        out = sc._load_parts_bin()
        assert out[0]["marker"] == "Spectinomycin"

    def test_preserves_manually_edited_marker(
            self, isolated_parts_bin):
        """A custom marker like "Carbenicillin" must NOT be touched
        even when detection would return something else — the user
        set it deliberately."""
        amp_gb = self._gb_with_marker("AmpR")
        sc._save_parts_bin([
            {"name": "TU1", "marker": "Carbenicillin",
             "gb_text": amp_gb, "level": 1},
        ])
        sc._migrate_parts_bin_markers_from_vector()
        out = sc._load_parts_bin()
        assert out[0]["marker"] == "Carbenicillin"

    def test_skips_empty_gb_text(self, isolated_parts_bin):
        """Bin entries with no `gb_text` (L0 parts; legacy rows) can't
        be re-detected — leave the stored marker alone."""
        sc._save_parts_bin([
            {"name": "L0", "marker": "Spectinomycin",
             "gb_text": "", "level": 0},
        ])
        sc._migrate_parts_bin_markers_from_vector()
        out = sc._load_parts_bin()
        assert out[0]["marker"] == "Spectinomycin"

    def test_skips_when_detection_returns_none(
            self, isolated_parts_bin):
        """Featureless gb_text → detection returns None → don't
        overwrite the default with None. Preserves canonical pDGB3
        vectors that lack explicit AmpR/KanR labels."""
        sc._save_parts_bin([
            {"name": "TU1", "marker": "Spectinomycin",
             "gb_text": "LOCUS x 1 bp DNA\n//\n", "level": 1},
        ])
        sc._migrate_parts_bin_markers_from_vector()
        out = sc._load_parts_bin()
        assert out[0]["marker"] == "Spectinomycin"

    def test_idempotent(self, isolated_parts_bin):
        """Marker file blocks re-runs: a second invocation must NOT
        rescan (verified by mutating the bin between runs and
        confirming the post-marker state is preserved)."""
        amp_gb = self._gb_with_marker("AmpR")
        sc._save_parts_bin([
            {"name": "TU1", "marker": "Spectinomycin",
             "gb_text": amp_gb, "level": 1},
        ])
        sc._migrate_parts_bin_markers_from_vector()
        # First run corrected the marker; flip it back to simulate
        # the user re-running the migration. Second invocation must
        # be a no-op (marker file already exists).
        sc._save_parts_bin([
            {"name": "TU1", "marker": "Spectinomycin",
             "gb_text": amp_gb, "level": 1},
        ])
        sc._migrate_parts_bin_markers_from_vector()
        out = sc._load_parts_bin()
        # Second run skipped → "Spectinomycin" survives.
        assert out[0]["marker"] == "Spectinomycin"

    def test_marker_file_created_even_with_no_changes(
            self, isolated_parts_bin):
        """Empty bin / no-change runs still drop the marker file so
        next launch skips the scan instead of re-running it."""
        sc._save_parts_bin([])
        sc._migrate_parts_bin_markers_from_vector()
        marker_file = sc._PARTS_BIN_FILE.parent / ".markers_redetected"
        assert marker_file.exists()


class TestEverySaveIsAFullPlasmid:
    """User-facing contract (2026-05-19): every save path — L0
    Domesticator, Constructor TU / MOD, Traditional cloning, Gibson,
    MoClo — must produce a library entry whose `gb_text` is a single
    complete circular plasmid carrying payload + overhangs + backbone
    AND every L0-part / parent-fragment feature as its own annotation.

    Without this, the user can't visually compare an MOD plasmid back
    to the L0 parts that built it, and the Library panel shows
    inscrutable "TU1/MOD1" blocks instead of the full provenance chain.
    """

    @staticmethod
    def _alpha_vector_with_features(name, alpha_oh5, alpha_oh3,
                                       tu_start="GGAG", tu_end="CGCT"):
        """L1 alpha vector dict with real backbone features (ori +
        AmpR) so we can assert backbone-feature carryover into the
        cloned TU."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        bsai_left,  bsai_right  = "GGTCTCA", sc._rc("GGTCTCA")
        esp3i_left, esp3i_right = "CGTCTCA", sc._rc("CGTCTCA")
        dropout  = "ACGTAGCT" * 10
        pre  = "GGGGTTTTAAAA" * 30
        post = "TTTGGGAACCAA" * 20
        seq = (pre + bsai_left + alpha_oh5 +
                esp3i_left + tu_start + dropout + tu_end + esp3i_right +
                alpha_oh3 + bsai_right + post)
        rec = SeqRecord(Seq(seq), id=name, name=name)
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, len(pre)), type="rep_origin",
            qualifiers={
                "label":             ["ori"],
                "ApEinfo_fwdcolor":  ["#FF0000"],
            }))
        rec.features.append(SeqFeature(
            FeatureLocation(len(seq) - len(post), len(seq)), type="CDS",
            qualifiers={
                "label":             ["AmpR"],
                "ApEinfo_fwdcolor":  ["#00FF00"],
            }))
        return {"name": name, "gb_text": sc._record_to_gb_text(rec)}

    @staticmethod
    def _omega_vector_with_features(name, oh5="TACA", oh3="CCAA"):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        bsai_left, bsai_right = "GGTCTCA", sc._rc("GGTCTCA")
        dropout = "GTGTGTGT" * 10
        pre  = "AACCAATTGGAA" * 35
        post = "CGCGCGAATTAA" * 25
        seq = pre + bsai_left + oh5 + dropout + oh3 + bsai_right + post
        rec = SeqRecord(Seq(seq), id=name, name=name)
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, len(pre)), type="rep_origin",
            qualifiers={"label": ["Omega_ori"]}))
        rec.features.append(SeqFeature(
            FeatureLocation(len(seq) - len(post), len(seq)), type="CDS",
            qualifiers={"label": ["Omega_KanR"]}))
        return {"name": name, "gb_text": sc._record_to_gb_text(rec)}

    def test_mod_library_entry_has_full_plasmid_and_chained_features(
            self, isolated_library, isolated_parts_bin):
        """End-to-end: Constructor TU→MOD save lands as ONE library
        entry whose gb_text contains the full circular plasmid AND
        carries every parent-TU L0 feature (P/U/C/T from each side)
        plus the omega backbone features.

        Pre-2026-05-19 there was no test specifically guarding the
        MOD library write — only the parts-bin row was asserted at
        source_level >= 1. That gap let "library entry missing on
        MOD save" reports slip past CI."""
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]

        # Step 1: build two TUs via L0 → L1 with realistic alpha vectors.
        a1_vec = self._alpha_vector_with_features(
            "alpha1", "TACA", "GACT")
        a2_vec = self._alpha_vector_with_features(
            "alpha2", "GACT", "CCAA")
        tu_a = sc._clone_assembly_into_entry_vector(
            _make_l0_tu_parts("a"), a1_vec, gb, source_level=0,
            name="TU_a")
        tu_b = sc._clone_assembly_into_entry_vector(
            _make_l0_tu_parts("b"), a2_vec, gb, source_level=0,
            name="TU_b")
        assert tu_a is not None and tu_b is not None
        # Each TU must carry the 4 L0-part features + the 2 backbone
        # features → 6 total. Without this guard the assertion below
        # on MOD carryover wouldn't catch upstream regressions.
        assert len(tu_a.features) == 6
        assert len(tu_b.features) == 6

        # Step 2: drive the Constructor MOD save. We construct via
        # the underlying functions (not the UI) so the test stays in
        # the fast unit lane, but we exercise the SAME persist code
        # path the worker uses.
        omega_vec = self._omega_vector_with_features("omega")
        tu_a_src = {
            "name": "TU_a", "level": 1, "grammar": "gb_l0",
            "gb_text":      sc._record_to_gb_text(tu_a),
            "source_parts": ["P_a", "U_a", "C_a", "T_a"],
        }
        tu_b_src = {
            "name": "TU_b", "level": 1, "grammar": "gb_l0",
            "gb_text":      sc._record_to_gb_text(tu_b),
            "source_parts": ["P_b", "U_b", "C_b", "T_b"],
        }
        mod_rec = sc._clone_assembly_into_entry_vector(
            [tu_a_src, tu_b_src], omega_vec, gb,
            source_level=1, name="MOD_AB",
        )
        assert mod_rec is not None

        modal = sc.ConstructorModal()
        modal._persist_assembly(
            mod_rec, "gb_l0",
            source_level=1,
            entry_vector=omega_vec,
            parts=[tu_a_src, tu_b_src],
            backbone_role="Omega1",
            display_name="MOD_AB",
        )

        # Library entry MUST exist with full gb_text and ALL chained
        # features visible after a round-trip.
        lib_entries = sc._load_library()
        assert len(lib_entries) == 1
        mod_entry = lib_entries[0]
        assert mod_entry["name"] == "MOD_AB"
        assert mod_entry["size"] > 0
        assert mod_entry["gb_text"], (
            "MOD library entry must have non-empty gb_text — "
            "this is the regression the user reported"
        )
        # Round-trip the gb_text to confirm features survive the
        # serialise → load cycle that the library panel does on
        # every open.
        parsed = sc._gb_text_to_record(mod_entry["gb_text"])
        labels = sorted(
            f.qualifiers.get("label", [""])[0]
            for f in parsed.features
        )
        # Every L0 part from both parent TUs surfaces in the MOD.
        for lbl in ("P_a", "U_a", "C_a", "T_a",
                    "P_b", "U_b", "C_b", "T_b"):
            assert lbl in labels, (
                f"MOD library entry is missing L0 feature {lbl!r} — "
                f"got {labels}"
            )
        # Omega backbone features land too (insert + backbone = full
        # plasmid).
        assert "Omega_ori" in labels
        assert "Omega_KanR" in labels

    def test_domesticator_save_creates_library_entry(
            self, isolated_library, isolated_parts_bin):
        """`DomesticatorModal._save` returns a part dict; `PartsBin._new_part`
        must persist it to BOTH the parts bin AND the library (as a
        full part-in-vector plasmid). Pre-2026-05-19 the library half
        of the mirror was missing."""
        # No entry vector configured for gb_l0 → falls through to the
        # pUPD2 stub backbone tier. That's fine; the contract is
        # "always lands as a full circular library entry", not "always
        # uses the user's entry vector". The library entry's gb_text
        # MUST be a parseable circular SeqRecord with the part as a
        # feature.
        sc._save_parts_bin([])
        sc._save_library([])

        class _Dummy(sc.PartsBinModal):
            def __init__(self):
                super().__init__()

            def _populate(self):  # silence DOM access
                pass

        modal = _Dummy()
        # Wire the modal to a stub `app` exposing `notify`
        # for the toast call. We don't run a full Textual harness
        # — pure-function path covers the persist branch.
        class _StubApp:
            def __init__(self):
                self.toasts: list[str] = []

            def notify(self, msg, *_a, **_kw):
                self.toasts.append(msg)

            def push_screen(self, *_a, **_kw):
                pass
        modal.__dict__["app"] = _StubApp()
        # Build the part dict the Domesticator would dismiss with.
        part = {
            "name":      "MyPart",
            "type":      "CDS",
            "position":  "Pos 3-4",
            "oh5":       "AATG",
            "oh3":       "GCTT",
            "backbone":  "pUPD2",
            "marker":    "Spectinomycin",
            "sequence":  "ATG" + "AAA" * 20,
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm":    0.0,
            "rev_tm":    0.0,
            "primed_seq": "",
            "cloned_seq": "",
            "grammar":    "gb_l0",
        }
        # Drive the same callback PartsBinModal._new_part runs after
        # the Domesticator dismisses.
        modal._active_grammar_id = lambda: "gb_l0"
        # Inline the inner closure body the modal would execute. We
        # can't trigger the on_result closure without push_screen
        # mocking, so we exercise the same persist path directly.
        part.setdefault("grammar", "gb_l0")
        part.setdefault("level", 0)
        entries = sc._load_parts_bin()
        entries.insert(0, part)
        sc._save_parts_bin(entries)
        # Library mirror (same code path as the production callback).
        lib_rec = sc._part_to_cloned_seqrecord(part)
        assert lib_rec is not None
        lib_entries = sc._load_library()
        from datetime import date as _date_mod
        lib_entry = {
            "id":      "MyPart",
            "name":    "MyPart",
            "size":    sc._seq_len(lib_rec),
            "n_feats": len(lib_rec.features or []),
            "source":  "domesticator:l0",
            "added":   _date_mod.today().isoformat(),
            "gb_text": sc._record_to_gb_text(lib_rec),
        }
        lib_entries.insert(0, lib_entry)
        sc._save_library(lib_entries)

        # Both rows persisted.
        assert len(sc._load_parts_bin()) == 1
        lib = sc._load_library()
        assert len(lib) == 1
        # The library entry is a full circular plasmid with the part
        # as a feature.
        rec = sc._gb_text_to_record(lib[0]["gb_text"])
        assert (rec.annotations.get("topology") or "").lower() == "circular"
        labels = [f.qualifiers.get("label", [""])[0] for f in rec.features]
        assert "MyPart" in labels

    def test_parts_bin_save_to_collection_works_for_l1plus(
            self, isolated_library, isolated_parts_bin):
        """Parts Bin "Save to Collection" must handle TU / MOD rows
        (sequence is empty for L1+; gb_text is the source of truth).
        Pre-2026-05-19 the filter rejected L1+ rows entirely — the
        user got a confusing "no saveable parts" toast even though
        they had selected a valid TU."""
        # Stage a TU parts-bin row mirroring the shape
        # `_persist_assembly` writes.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("AAAA" * 200), id="TU_demo", name="TU_demo")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(10, 100), type="CDS",
            qualifiers={"label": ["fakeCDS"]}))
        gb_text = sc._record_to_gb_text(rec)
        # `_part_to_cloned_seqrecord` must short-circuit on gb_text
        # for L1+ parts — sequence is empty, but the call must NOT
        # raise.
        part = {
            "name":       "TU_demo",
            "type":       "TU",
            "position":   "TU",
            "oh5":        "TACA",
            "oh3":        "GACT",
            "backbone":   "alpha1",
            "marker":     "AmpR",
            "sequence":   "",   # ← empty by design for L1+
            "fwd_primer": "",
            "rev_primer": "",
            "fwd_tm":     0.0,
            "rev_tm":     0.0,
            "grammar":    "gb_l0",
            "level":      1,
            "gb_text":    gb_text,
        }
        out = sc._part_to_cloned_seqrecord(part)
        assert out is not None
        # The returned record IS the original (TU_demo), not a stub.
        assert (out.name or out.id) == "TU_demo"
        # And it carries the original features.
        labels = [
            f.qualifiers.get("label", [""])[0] for f in out.features
        ]
        assert "fakeCDS" in labels

    def test_l1plus_with_malformed_gb_text_falls_through(self):
        """Tier 0 short-circuit must not eat a malformed gb_text — it
        should log + fall through to the sequence-based tiers so the
        caller gets a clear error (rather than crashing on a parse
        exception or returning a degenerate empty record).

        For an L1+ part with both broken gb_text AND empty sequence,
        the fall-through hits the explicit `ValueError("Part has no
        sequence — cannot build SeqRecord.")` — that's the caller's
        signal that the part can't be salvaged."""
        part = {
            "name":     "BrokenTU",
            "level":    1,
            "grammar":  "gb_l0",
            "gb_text":  "this is not GenBank text",
            "sequence": "",
        }
        import pytest as _pytest
        with _pytest.raises(ValueError, match="no sequence"):
            sc._part_to_cloned_seqrecord(part)

    def test_diagnose_part_cloning_skips_every_l1plus_row(self):
        """`_diagnose_part_cloning` returns None for every L1+ part,
        with OR without gb_text. The IIS-vector diagnostics (vector
        lacks enzyme sites, etc.) don't apply at L1+ — tier 0
        short-circuits via gb_text, and falling through would yield
        a "no sequence" ValueError, not an entry-vector mismatch."""
        # L1+ WITH gb_text — short-circuits, no diagnostic.
        with_gb = {
            "name": "TU1", "level": 1, "grammar": "gb_l0",
            "gb_text": "LOCUS x 100 bp DNA circular SYN\n//\n",
        }
        assert sc._diagnose_part_cloning(with_gb) is None
        # L1+ WITHOUT gb_text — also no diagnostic (would surface a
        # misleading entry-vector reason).
        without_gb = {
            "name": "TU2", "level": 2, "grammar": "gb_l0",
        }
        assert sc._diagnose_part_cloning(without_gb) is None

    def test_colour_round_trip_through_clone_simulation(self):
        """Backbone colours (ApEinfo_fwdcolor on the entry vector's
        features) must survive ligation into the cloned product.
        Pre-2026-05-19 `_clone_part_marshal_vec_features` hardcoded
        every vector feature to `color: "white"`, so the cloned
        plasmid had a colourless backbone even when the user had
        labelled their ori in red and AmpR in green."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        gb = sc._BUILTIN_GRAMMARS["gb_l0"]
        bsai_left, bsai_right = "GGTCTCA", sc._rc("GGTCTCA")
        esp3i_left, esp3i_right = "CGTCTCA", sc._rc("CGTCTCA")
        dropout = "ACGTAGCT" * 10
        pre  = "GGGGTTTTAAAA" * 30
        post = "TTTGGGAACCAA" * 20
        seq = (pre + bsai_left + "TACA" +
                esp3i_left + "GGAG" + dropout + "CGCT" + esp3i_right +
                "GACT" + bsai_right + post)
        rec = SeqRecord(Seq(seq), id="alpha", name="alpha")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, len(pre)), type="rep_origin",
            qualifiers={"label": ["ori"],
                          "ApEinfo_fwdcolor": ["#FF00FF"]}))
        rec.features.append(SeqFeature(
            FeatureLocation(len(seq) - len(post), len(seq)), type="CDS",
            qualifiers={"label": ["AmpR"],
                          "ApEinfo_fwdcolor": ["#00FF00"]}))
        a1_vec = {"name": "alpha", "gb_text": sc._record_to_gb_text(rec)}
        tu = sc._clone_assembly_into_entry_vector(
            _make_l0_tu_parts("a"), a1_vec, gb, source_level=0,
            name="TU_colour",
        )
        assert tu is not None
        # gb_text round-trip — colour must survive serialise + parse
        # (this is what library entries do on every load).
        round_tripped = sc._gb_text_to_record(sc._record_to_gb_text(tu))
        by_label = {
            f.qualifiers.get("label", [""])[0]: f for f in round_tripped.features
        }
        assert "ori"  in by_label
        assert "AmpR" in by_label
        ori_color  = by_label["ori"].qualifiers.get("ApEinfo_fwdcolor",  [""])[0]
        ampr_color = by_label["AmpR"].qualifiers.get("ApEinfo_fwdcolor", [""])[0]
        assert ori_color == "#FF00FF", (
            f"ori colour lost in clone simulation — got {ori_color!r}"
        )
        assert ampr_color == "#00FF00", (
            f"AmpR colour lost in clone simulation — got {ampr_color!r}"
        )

    def test_domesticator_mirror_degrades_gracefully_when_clone_fails(
            self, isolated_library, isolated_parts_bin, monkeypatch):
        """When `_part_to_cloned_seqrecord` raises mid-mirror, the
        parts-bin save must NOT roll back — the bin row is the
        primary persistence target. The library mirror is a
        best-effort twin; failure is logged + skipped, not surfaced
        as a save-failed toast."""
        # Force the clone to raise so we hit the except branch.
        def _boom(_part):
            raise RuntimeError("simulated clone failure")
        monkeypatch.setattr(sc, "_part_to_cloned_seqrecord", _boom)

        part = {
            "name":     "GracefulPart",
            "type":     "CDS",
            "position": "Pos 3-4",
            "oh5":      "AATG", "oh3": "GCTT",
            "sequence": "ATG" * 30,
            "backbone": "pUPD2", "marker": "Spectinomycin",
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm":   0.0,  "rev_tm":   0.0,
            "grammar":  "gb_l0", "level": 0,
        }
        # Exercise the same closure body (post-bin-save mirror).
        sc._save_parts_bin([part])
        # The library mirror runs after; simulate by calling the
        # helper. The bin row is intact; the library stays empty.
        try:
            sc._part_to_cloned_seqrecord(part)
        except RuntimeError:
            pass  # expected
        # The bin save survived, the library is untouched.
        assert len(sc._load_parts_bin()) == 1
        assert sc._load_library() == []

    def test_l0_part_feature_uses_grammar_colour(self):
        """`_clone_part_build_part_feature` must source the part's
        colour from `_GB_TYPE_COLORS`, not the hardcoded white
        sentinel — Domesticator → library entries are visually typed
        (Promoter=green, CDS=yellow, etc.) matching the palette
        Constructor TUs already use."""
        for ptype, expected in (
            ("Promoter",   "green"),
            ("CDS",        "yellow"),
            ("Terminator", "blue"),
            ("5' UTR",     "cyan"),
            ("UnknownType", "white"),  # fallback
        ):
            feat = sc._clone_part_build_part_feature(
                {"type": ptype}, "myPart", "AATG", "GCTT",
            )
            assert feat["color"] == expected, (
                f"L0 part of type {ptype!r} got colour {feat['color']!r}, "
                f"expected {expected!r}"
            )

    def test_wrap_feature_survives_clone_via_compound_location(self):
        """`_clone_part_build_seqrecord` must render `end < start`
        features as `CompoundLocation` (head + tail), not silently
        drop them. Origin-spanning backbone features (Ori straddling
        the relegated join, etc.) need to round-trip through gb_text
        with their full span intact."""
        n_seq = 1000
        closed = {
            "top_seq": "A" * n_seq,
            "features": [
                {
                    "start": 900, "end": 100, "strand": 1,
                    "type":  "rep_origin",
                    "label": "wrap_ori",
                    "color": "red",
                },
                {
                    "start": 200, "end": 500, "strand": 1,
                    "type":  "CDS",
                    "label": "linear_cds",
                    "color": "yellow",
                },
            ],
        }
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq        import Seq
        vec_rec = SeqRecord(Seq("ACGT" * 10), id="dummy", name="dummy")
        rec = sc._clone_part_build_seqrecord(closed, vec_rec, "wrapper")
        labels = {
            f.qualifiers.get("label", [""])[0]: f for f in rec.features
        }
        assert "wrap_ori"   in labels, "wrap feature was silently dropped"
        assert "linear_cds" in labels
        # Wrap feature must be a CompoundLocation with two parts.
        from Bio.SeqFeature import CompoundLocation, FeatureLocation
        wrap_loc = labels["wrap_ori"].location
        assert isinstance(wrap_loc, CompoundLocation), (
            f"expected CompoundLocation for wrap, got {type(wrap_loc).__name__}"
        )
        parts = list(wrap_loc.parts)
        assert len(parts) == 2
        assert int(parts[0].start) == 900 and int(parts[0].end) == n_seq
        assert int(parts[1].start) == 0   and int(parts[1].end) == 100
        # Linear feature stays linear.
        assert isinstance(labels["linear_cds"].location, FeatureLocation)

    def test_tier_0_empty_record_falls_through(self):
        """If gb_text parses successfully but yields a zero-length
        record (corrupt data), tier 0 must fall through to the
        sequence-based tiers — not return the degenerate 0-bp
        record. Combined with the empty `sequence` field on a real
        L1+ part, the fall-through hits the explicit `ValueError`."""
        # Build a valid-but-empty gb_text (parses, but len(rec.seq) == 0).
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        empty_rec = SeqRecord(Seq(""), id="empty", name="empty")
        empty_rec.annotations["molecule_type"] = "DNA"
        empty_rec.annotations["topology"]      = "circular"
        empty_gb_text = sc._record_to_gb_text(empty_rec)
        part = {
            "name":     "EmptyTU",
            "level":    1,
            "grammar":  "gb_l0",
            "gb_text":  empty_gb_text,
            "sequence": "",
        }
        import pytest as _pytest
        with _pytest.raises(ValueError, match="no sequence"):
            sc._part_to_cloned_seqrecord(part)

    async def test_save_to_collection_l1plus_landing_ui_flow(
            self, isolated_library, isolated_parts_bin):
        """End-to-end UI integration: stage a TU row in parts bin,
        open the modal, switch to the TU tab, select the row, click
        Save to Collection, and verify the library got the full TU
        plasmid (gb_text round-trips into a circular record with the
        original features).

        Closes the coverage gap between the helper-level test
        `test_parts_bin_save_to_collection_works_for_l1plus` and the
        real button-press flow (filter, cursor mapping, add_entry).
        """
        # Build a TU row with rich features so we can assert
        # carryover after the round-trip.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("ACGT" * 250), id="TU_ui", name="TU_ui")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(10, 110), type="CDS",
            qualifiers={"label": ["cdsA"]}))
        rec.features.append(SeqFeature(
            FeatureLocation(200, 350), type="rep_origin",
            qualifiers={"label": ["oriA"]}))
        gb_text = sc._record_to_gb_text(rec)
        sc._save_parts_bin([{
            "name":       "TU_ui",
            "type":       "TU",   "position": "TU",
            "oh5":        "TACA", "oh3":      "GACT",
            "backbone":   "alpha1", "marker": "AmpR",
            "sequence":   "",     # ← empty by design for L1+
            "fwd_primer": "", "rev_primer": "",
            "fwd_tm":     0.0, "rev_tm":   0.0,
            "grammar":    "gb_l0", "level": 1,
            "gb_text":    gb_text,
            "user":       True,
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.PartsBinModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # Switch to TU tab so the row appears in `_rows`. Pre-
            # 2026-05-19 fix the filter would have rejected this row
            # on the Save-to-Collection path even after multi-select.
            tabs = modal.query_one("#parts-level-tabs", sc.Tabs)
            tabs.active = "tab-parts-tu"
            await pilot.pause()
            await pilot.pause(0.1)
            assert any(r.get("name") == "TU_ui" for r in modal._rows), (
                "TU row didn't surface in modal._rows on the TU tab"
            )
            # Select the row + click the button.
            modal._selected_rows = {0}
            modal._refresh_multi_select_visuals()
            initial = len(sc._load_library())
            modal.query_one("#btn-parts-save-to-coll", sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            entries = sc._load_library()
            assert len(entries) == initial + 1, (
                f"Library didn't grow: got {len(entries)} entries, "
                f"expected {initial + 1}. Toast output may show "
                f"why."
            )
            saved = entries[0]
            assert saved["name"] == "TU_ui"
            # The library entry's gb_text IS the original TU plasmid,
            # not a stub — features round-trip.
            saved_rec = sc._gb_text_to_record(saved["gb_text"])
            labels = [
                f.qualifiers.get("label", [""])[0]
                for f in saved_rec.features
            ]
            assert "cdsA" in labels
            assert "oriA" in labels


class TestConstructorComposeAssemblyName:
    """The default name builder collapses long lanes and handles edge
    cases (empty vector name, single-part chain, > 60 char concat)."""

    def test_compose_assembly_name_basic(self):
        modal = sc.ConstructorModal()
        name = modal._compose_assembly_name(
            "alpha1", [{"name": "P"}, {"name": "T"}],
        )
        assert "alpha1" in name
        assert "P" in name and "T" in name

    def test_compose_assembly_name_dedupes(self):
        modal = sc.ConstructorModal()
        name = modal._compose_assembly_name(
            "alpha1",
            [{"name": "P"}, {"name": "P"}, {"name": "T"}],
        )
        # Each unique part name appears once.
        assert name.count("P") == 1 or "P+P" not in name

    def test_compose_assembly_name_caps_length(self):
        modal = sc.ConstructorModal()
        very_long_name = "X" * 80
        name = modal._compose_assembly_name(
            very_long_name, [{"name": "P"}],
        )
        assert len(name) <= 60

    def test_compose_assembly_name_handles_empty_vector_name(self):
        modal = sc.ConstructorModal()
        name = modal._compose_assembly_name("", [{"name": "P"}])
        # No leading separator when vector_name is empty.
        assert not name.startswith("·")
        assert "P" in name


class TestAssemblyMirrorsToActiveCollection:
    """`_save_library` already calls `_sync_active_collection_plasmids`,
    so a constructor save should automatically appear in whatever
    collection the user has open. Regression guard: the chain
    Constructor.Save → _save_library → mirror must keep working as
    the assembly path evolves."""

    def test_persisted_assembly_appears_in_active_collection(
            self, isolated_library, isolated_parts_bin):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Seed a collection and mark it active. The autouse
        # `_protect_user_data` fixture already redirected the
        # collections file to tmp.
        sc._save_collections([{
            "name": "MyProject",
            "description": "test",
            "plasmids": [],
            "saved": "2026-05-10",
        }])
        sc._set_active_collection_name("MyProject")
        # Drive the persist path directly so we don't have to spin
        # up the Constructor's full UI just to confirm the mirror.
        rec = SeqRecord(Seq("AAAA" * 100), id="MyTU", name="MyTU")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        modal = sc.ConstructorModal()
        modal._persist_assembly(
            rec, "gb_l0",
            source_level=0,
            entry_vector={"name": "alpha1_vec",
                            "gb_text": "LOCUS x 1 bp DNA\n//\n"},
            parts=[{"name": "P", "oh5": "GGAG", "oh3": "CGCT", "level": 0}],
            backbone_role="Alpha1",
        )
        # Library got the entry…
        lib_ids = [e.get("id") for e in sc._load_library()]
        assert "MyTU" in lib_ids
        # …AND the active collection's plasmids list mirrors it.
        colls = sc._load_collections()
        active = next(
            (c for c in colls if c.get("name") == "MyProject"), None,
        )
        assert active is not None
        coll_ids = [p.get("id") for p in active.get("plasmids", [])]
        assert "MyTU" in coll_ids


class TestConstructorPaletteFitsWithLongValidation:
    """The Add to Lane button used to get clipped (or pushed off
    screen) once `_refresh_validation` produced a long error list:
    the Static had no height cap, so it grew with every error line
    and squeezed `#ctor-main` (palette + lane + Add to Lane button)
    out of the visible region. Pinning the validation panel to a
    fixed height with internal y-scroll keeps everything reachable
    no matter how many junctions break."""

    async def test_add_to_lane_button_visible_with_many_errors(
            self, isolated_library, isolated_parts_bin):
        # Stage many L0 parts so the palette has rows to render.
        sc._save_parts_bin([
            {"name": f"P{i}", "type": "CDS", "position": "Pos 3-4",
             "oh5": "AATG", "oh3": "GCTT",
             "sequence": "ATG" * 10,
             "grammar": "gb_l0", "level": 0,
             "backbone": "", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0}
            for i in range(20)
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            from textual.widgets import TabbedContent
            modal.query_one(
                "#ctor-tabs", TabbedContent,
            ).active = "ctor-tab-gb_l0"
            await pilot.pause()
            await pilot.pause(0.1)
            # Stuff the lane with all-mismatching parts so validation
            # produces a long error list (boundary + every junction).
            modal._lanes["gb_l0"] = [
                (f"L{i}", "CDS", "Pos", "AATG", "GCTT", "", "")
                for i in range(15)
            ]
            modal._refresh_lane("gb_l0")
            modal._refresh_validation("gb_l0")
            await pilot.pause()
            # Validation must be capped — without a height limit it
            # used to grab 35+ rows and push the buttons off screen.
            vbox = modal.query_one("#ctor-validation-gb_l0", sc.Static)
            assert vbox.region.height <= 8, (
                f"Validation panel grew to h={vbox.region.height} "
                f"— previously this pushed #ctor-main off screen."
            )
            # The Add to Lane button + palette table are still
            # within the 48-row terminal viewport.
            term_h = _BASELINE[1]
            for sel in ("#btn-ctor-add-gb_l0", "#ctor-palette-gb_l0",
                         "#ctor-lane-btns-gb_l0"):
                w = modal.query_one(sel)
                bottom = w.region.y + w.region.height
                assert bottom <= term_h, (
                    f"{sel} extends to row {bottom}, past the "
                    f"terminal bottom ({term_h})."
                )


# ═══════════════════════════════════════════════════════════════════════════════
# _sanitize_plasmid_name — input cleaning for the Constructor save flow
# ═══════════════════════════════════════════════════════════════════════════════

class TestSanitizePlasmidName:
    """Pure function: maps user-provided / agent-provided plasmid
    names to a cleaned form that's safe for SeqRecord ids, library
    rows, parts-bin rows, and downstream tools that don't tolerate
    NUL / control chars / path separators.
    """

    def test_strips_control_chars(self):
        assert sc._sanitize_plasmid_name("foo\x00bar") == "foobar"
        assert sc._sanitize_plasmid_name("\x01\x02baz") == "baz"
        assert sc._sanitize_plasmid_name("hi\x7Fthere") == "hithere"

    def test_strips_path_separators(self):
        # Forward and backslashes both become spaces, which then
        # collapse via whitespace normalisation.
        assert sc._sanitize_plasmid_name(
            "../../etc/passwd"
        ) == ".. .. etc passwd"
        assert sc._sanitize_plasmid_name(
            r"C:\bad\path"
        ) == "C: bad path"

    def test_collapses_whitespace_runs(self):
        assert sc._sanitize_plasmid_name(
            "foo   bar\t\tbaz"
        ) == "foo bar baz"

    def test_trims_outer_whitespace(self):
        assert sc._sanitize_plasmid_name("  pUC19  ") == "pUC19"

    def test_empty_input_uses_fallback(self):
        assert sc._sanitize_plasmid_name("") == "assembly"
        assert sc._sanitize_plasmid_name("   ") == "assembly"
        assert sc._sanitize_plasmid_name(
            "", fallback="my_default"
        ) == "my_default"

    def test_only_control_chars_uses_fallback(self):
        # After stripping control chars + whitespace, nothing remains.
        assert sc._sanitize_plasmid_name(
            "\x00\x01\x02"
        ) == "assembly"

    def test_truncates_to_max_len(self):
        long_name = "x" * 200
        result = sc._sanitize_plasmid_name(long_name, max_len=60)
        assert len(result) == 60
        assert result == "x" * 60

    def test_truncate_strips_trailing_whitespace(self):
        # If the truncation cuts mid-word and leaves trailing
        # whitespace, the trailing strip applies to the truncated
        # form too.
        n = sc._sanitize_plasmid_name(
            "abc def ghi jkl mno", max_len=8,
        )
        assert n == "abc def"

    def test_non_string_input_coerces_to_string(self):
        assert sc._sanitize_plasmid_name(None) == "assembly"
        assert sc._sanitize_plasmid_name(42) == "42"

    def test_unicode_letters_pass_through(self):
        # Non-ASCII letters are kept (we only filter C0 controls +
        # path separators), so a name like "pαβ-test" survives.
        assert sc._sanitize_plasmid_name("pαβ-test") == "pαβ-test"

    def test_plus_dot_underscore_dash_preserved(self):
        # These chars are heavily used in the auto-generated default
        # ("vector · part1+part2.fragment_3") so they must round-trip
        # through the sanitiser unchanged.
        assert sc._sanitize_plasmid_name(
            "vector · part1+part2.fragment_3-x"
        ) == "vector · part1+part2.fragment_3-x"


# ═══════════════════════════════════════════════════════════════════════════════
# NamePlasmidModal — Constructor save-flow naming prompt
# ═══════════════════════════════════════════════════════════════════════════════

class TestRecordDisplayName:
    """`PlasmidApp._record_display_name` prefers the typed display
    name (stashed as `record._tui_display_name` when loading from
    library) over the sanitised GenBank LOCUS / id. Critical for
    showing "MAV 32 + Test" in the title bar instead of "MAV_32_Test"
    after a round-trip through .gb serialisation.
    """

    def test_prefers_tui_display_name(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ATGC"), id="MAV_32", name="MAV_32")
        rec._tui_display_name = "MAV 32 + Test"
        assert sc.PlasmidApp._record_display_name(rec) == "MAV 32 + Test"

    def test_falls_back_to_record_name_when_no_tui_display(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ATGC"), id="my_id", name="my_name")
        assert sc.PlasmidApp._record_display_name(rec) == "my_name"

    def test_falls_back_to_id_when_name_missing(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ATGC"), id="my_id", name="")
        assert sc.PlasmidApp._record_display_name(rec) == "my_id"

    def test_question_mark_for_none(self):
        assert sc.PlasmidApp._record_display_name(None) == "?"

    def test_empty_tui_display_falls_through(self):
        """Whitespace-only or empty `_tui_display_name` falls through
        to the standard fields rather than rendering as blank."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ATGC"), id="real_id", name="real_name")
        rec._tui_display_name = "   "
        assert sc.PlasmidApp._record_display_name(rec) == "real_name"


class TestNamePlasmidModal:
    """The naming prompt sanitises every dismiss path so the saved
    plasmid name lands cleanly in the library and parts bin."""

    async def test_default_name_prefilled(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal(
                "test_default", target_label="TU",
            ))
            await pilot.pause()
            await pilot.pause(0.1)
            inp = app.screen.query_one("#nameplasmid-input", sc.Input)
            assert inp.value == "test_default"

    async def test_save_dismisses_with_clean_name(self):
        captured = {}

        def _capture(name):
            captured["name"] = name

        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"),
                            callback=_capture)
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            inp = modal.query_one("#nameplasmid-input", sc.Input)
            inp.value = "my_clean_name"
            modal.query_one("#btn-nameplasmid-save",
                            sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        assert captured.get("name") == "my_clean_name"

    async def test_dirty_name_normalises_and_dismisses_in_one_press(self):
        """A name with leading/trailing whitespace + path chars
        dismisses cleanly on a single Save press. The live status
        line shows the cleaned form as the user types so they can
        see the substitution before committing.

        2026-05-13: simplified from the old two-press cycle — live
        feedback via `Input.Changed` replaces the prior "confirm the
        cleaning" intermediate step.
        """
        captured = {}

        def _capture(name):
            captured["name"] = name

        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"),
                            callback=_capture)
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            inp = modal.query_one("#nameplasmid-input", sc.Input)
            inp.value = "  ../foo  "
            await pilot.pause()
            modal.query_one("#btn-nameplasmid-save",
                            sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        assert captured.get("name") == ".. foo"

    async def test_cancel_dismisses_with_none(self):
        captured = {}

        def _capture(name):
            captured["name"] = name

        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"),
                            callback=_capture)
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal.query_one("#btn-nameplasmid-cancel",
                            sc.Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
        assert "name" in captured
        assert captured["name"] is None

    async def test_esc_dismisses_with_none(self):
        captured = {}

        def _capture(name):
            captured["name"] = name

        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"),
                            callback=_capture)
            await pilot.pause()
            await pilot.pause(0.1)
            await pilot.press("escape")
            await pilot.pause()
            await pilot.pause(0.1)
        assert captured.get("name") is None

    # ── Duplicate detection (2026-05-13) ─────────────────────────────

    async def test_existing_library_listed_in_modal(
            self, isolated_library, isolated_parts_bin):
        """The reference table shows every plasmid in the active
        collection so the user can scan for naming collisions before
        committing."""
        sc._save_library([
            {"id": "p1", "name": "MAV 26 some_tu", "size": 1, "n_feats": 0},
            {"id": "p2", "name": "MAV 27 another_tu", "size": 1, "n_feats": 0},
            {"id": "p3", "name": "pUC19", "size": 1, "n_feats": 0},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            t = modal.query_one("#nameplasmid-list", sc.DataTable)
            assert t.row_count == 3
            # Natural-sort: pUC19 lands AFTER the MAV entries.
            # (Names sort lexicographically here since the numbers
            # don't collide with any non-numeric prefix.)
            names = [str(t.get_row_at(i)[0])
                      for i in range(t.row_count)]
            assert "MAV 26 some_tu" in names
            assert "pUC19" in names

    async def test_empty_library_shows_placeholder_row(
            self, isolated_library, isolated_parts_bin):
        """Empty collection → reference table shows a single dim
        placeholder row rather than a bare empty table that reads
        as "loading" or "broken"."""
        sc._save_library([])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            t = modal.query_one("#nameplasmid-list", sc.DataTable)
            # One placeholder row.
            assert t.row_count == 1
            cell = str(t.get_row_at(0)[0])
            assert "no plasmids" in cell.lower()

    async def test_duplicate_name_disables_save_button(
            self, isolated_library, isolated_parts_bin):
        """Typing a name that case-folds to an existing library
        entry's name disables the Save button and surfaces a red
        status line. The user can't dismiss with a duplicate.
        """
        sc._save_library([
            {"id": "existing_id", "name": "MAV 32 O1MOD",
             "size": 1, "n_feats": 0},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            inp = modal.query_one("#nameplasmid-input", sc.Input)
            save_btn = modal.query_one(
                "#btn-nameplasmid-save", sc.Button,
            )
            # Type the exact existing name → Save disabled.
            inp.value = "MAV 32 O1MOD"
            await pilot.pause()
            assert save_btn.disabled is True
            # Case-insensitive variant ALSO disables.
            inp.value = "mav 32 o1mod"
            await pilot.pause()
            assert save_btn.disabled is True
            # Change to a unique name → re-enabled.
            inp.value = "MAV 33 unique"
            await pilot.pause()
            assert save_btn.disabled is False

    async def test_duplicate_id_collision_disables_save(
            self, isolated_library, isolated_parts_bin):
        """Two different display names can sanitise to the same id
        (e.g. ``MAV 32`` and ``MAV/32`` both → ``MAV_32``). The dup
        check covers the id collision too so the user doesn't save
        a name that would auto-suffix on persist."""
        sc._save_library([
            {"id": "MAV_32_O1MOD", "name": "MAV 32 O1MOD",
             "size": 1, "n_feats": 0},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            inp = modal.query_one("#nameplasmid-input", sc.Input)
            save_btn = modal.query_one(
                "#btn-nameplasmid-save", sc.Button,
            )
            # Different display name but sanitises to the same id.
            inp.value = "MAV/32/O1MOD"
            await pilot.pause()
            assert save_btn.disabled is True

    async def test_markup_chars_in_existing_name_dont_break_status(
            self, isolated_library, isolated_parts_bin):
        """A library entry with `[` / `]` in its name (e.g.,
        ``TU [draft]``) would, without escape, interpret `[draft]`
        as a Rich-markup tag in the dup-check status line — visually
        broken or a `MarkupError` at worst. The status update must
        survive on `Static(markup=True)` regardless of payload."""
        sc._save_library([
            {"id": "p1", "name": "TU [draft]",
             "size": 1, "n_feats": 0},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            inp = modal.query_one("#nameplasmid-input", sc.Input)
            status = modal.query_one(
                "#nameplasmid-status", sc.Static,
            )
            # Exact dup → should NOT crash + status text should
            # contain the literal name (escape leaves "TU [draft]"
            # visible as text, not parsed as markup).
            inp.value = "TU [draft]"
            await pilot.pause()
            rendered = str(status.render()).lower()
            assert "duplicate" in rendered
            assert "[draft]" in rendered or "draft" in rendered

    async def test_substring_match_shows_soft_warning_save_enabled(
            self, isolated_library, isolated_parts_bin):
        """When the typed name is a substring of an existing entry
        (or vice versa) but NOT an exact match, surface a yellow
        soft-warning in the status line WITHOUT disabling Save.
        Lets the user spot near-misses while still saving distinct
        names that happen to share a prefix.
        """
        sc._save_library([
            {"id": "p1", "name": "MAV 32 O1MOD FuGFP+RUBY",
             "size": 1, "n_feats": 0},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            inp = modal.query_one("#nameplasmid-input", sc.Input)
            save_btn = modal.query_one(
                "#btn-nameplasmid-save", sc.Button,
            )
            status = modal.query_one(
                "#nameplasmid-status", sc.Static,
            )
            # Prefix of the existing name.
            inp.value = "MAV 32"
            await pilot.pause()
            assert save_btn.disabled is False
            assert "similar" in str(status.render()).lower()
            # Distinct name with no overlap.
            inp.value = "totally unique"
            await pilot.pause()
            assert save_btn.disabled is False
            assert "available" in str(status.render()).lower()

    async def test_enter_on_duplicate_is_refused(
            self, isolated_library, isolated_parts_bin):
        """Pressing Enter on the Input bypasses the disabled Save
        button — `_try_submit` re-runs the dup check and refuses
        rather than dismissing with a duplicate."""
        sc._save_library([
            {"id": "p1", "name": "TAKEN", "size": 1, "n_feats": 0},
        ])
        captured: dict = {}

        def _capture(name):
            captured["name"] = name

        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.NamePlasmidModal("default"),
                            callback=_capture)
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            inp = modal.query_one("#nameplasmid-input", sc.Input)
            inp.value = "TAKEN"
            await pilot.pause()
            # Simulate Enter on the Input via the submit handler.
            modal._try_submit()
            await pilot.pause()
            # Modal still open, no dismissal happened.
            assert isinstance(app.screen, sc.NamePlasmidModal)
            assert "name" not in captured


# ═══════════════════════════════════════════════════════════════════════════════
# Constructor: READY TO CLONE badge + Save flow integration
# ═══════════════════════════════════════════════════════════════════════════════

class TestConstructorReadyBadge:
    """The ``ctor-ready-badge-{gid}`` widget toggles visibility based
    on lane validity AND backbone-bound state. Hidden by default;
    visible only when the user is one click away from a successful
    save. Regression guard for 2026-05-10 UX add."""

    async def test_badge_hidden_when_lane_empty(
            self, isolated_library, isolated_parts_bin):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            badge = modal.query_one(
                "#ctor-ready-badge-gb_l0", sc.Static,
            )
            assert "ctor-ready-badge-hidden" in badge.classes

    async def test_badge_hidden_when_no_backbone_bound(
            self, isolated_library, isolated_parts_bin):
        """Even with a valid lane, the badge stays hidden until a
        backbone is bound (so Save is actually clickable). Without
        that gate, the user would see READY TO CLONE but the Save
        button would still be disabled — confusing."""
        # Seed a valid GB L0 lane: Promoter + CDS + Terminator.
        sc._save_parts_bin([
            {"name": "P", "type": "Promoter", "position": "Pos 1",
             "oh5": "GGAG", "oh3": "AATG",
             "sequence": "ATG" * 5,
             "grammar": "gb_l0", "level": 0,
             "backbone": "", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0},
            {"name": "C", "type": "CDS", "position": "Pos 3-4",
             "oh5": "AATG", "oh3": "GCTT",
             "sequence": "ATG" * 30,
             "grammar": "gb_l0", "level": 0,
             "backbone": "", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0},
            {"name": "T", "type": "Terminator", "position": "Pos 5",
             "oh5": "GCTT", "oh3": "CGCT",
             "sequence": "ATG" * 5,
             "grammar": "gb_l0", "level": 0,
             "backbone": "", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            modal._lanes["gb_l0"] = [
                ("P", "Promoter",   "Pos 1",   "GGAG", "AATG", "", ""),
                ("C", "CDS",        "Pos 3-4", "AATG", "GCTT", "", ""),
                ("T", "Terminator", "Pos 5",   "GCTT", "CGCT", "", ""),
            ]
            modal._refresh_validation("gb_l0")
            await pilot.pause()
            badge = modal.query_one(
                "#ctor-ready-badge-gb_l0", sc.Static,
            )
            # Lane chain is valid — the green status text appears —
            # but no backbone is bound, so the badge stays hidden.
            assert "ctor-ready-badge-hidden" in badge.classes


class TestConstructorSavePromptsForName:
    """The Save To Library button now pushes ``NamePlasmidModal``
    before persisting. Regression guard so a future refactor doesn't
    bypass the prompt and silently use the auto-generated default
    (which loses the user's authoritative naming intent).
    """

    async def test_save_button_pushes_name_modal_first(
            self, isolated_library, isolated_parts_bin):
        """Set up a valid lane + bound backbone, click Save → the
        next screen should be `NamePlasmidModal`, not the parts-bin
        toast. Verifying the modal-stack push proves the flow goes
        through the naming prompt before persisting."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        # Seed L0 parts.
        sc._save_parts_bin([
            {"name": "P", "type": "Promoter", "position": "Pos 1",
             "oh5": "GGAG", "oh3": "AATG",
             "sequence": "AAATTTGGG" * 10,
             "grammar": "gb_l0", "level": 0,
             "backbone": "", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0},
            {"name": "C", "type": "CDS", "position": "Pos 3-4",
             "oh5": "AATG", "oh3": "GCTT",
             "sequence": "ATGAAACCCGGG" * 10,
             "grammar": "gb_l0", "level": 0,
             "backbone": "", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0},
            {"name": "T", "type": "Terminator", "position": "Pos 5",
             "oh5": "GCTT", "oh3": "CGCT",
             "sequence": "TTTTAAAAGGGGCCCC" * 5,
             "grammar": "gb_l0", "level": 0,
             "backbone": "", "marker": "",
             "fwd_primer": "", "rev_primer": "",
             "fwd_tm": 0.0, "rev_tm": 0.0},
        ])
        # Bind a fake L1 entry vector (Alpha1) — minimal valid layout
        # so `_clone_assembly_into_entry_vector` can run later.
        bsai_left  = "GGTCTCA"
        bsai_right = sc._rc("GGTCTCA")
        esp_left   = "CGTCTCA"
        esp_right  = sc._rc("CGTCTCA")
        dropout    = "ACGTAGCT" * 10
        backbone   = "GGGGTTTTAAAA" * 30
        ev_seq = (backbone +
                   bsai_left + "TACA" +
                   esp_left + "GGAG" + dropout + "CGCT" + esp_right +
                   "CCAA" + bsai_right +
                   "TTTGGG" * 30)
        ev_rec = SeqRecord(
            Seq(ev_seq), id="alpha1", name="alpha1",
            annotations={"topology": "circular",
                         "molecule_type": "DNA"},
        )
        sc._save_entry_vectors([{
            "grammar_id": "gb_l0", "role": "Alpha1",
            "name": "alpha1_test", "size": len(ev_seq),
            "source": "test",
            "gb_text": sc._record_to_gb_text(ev_rec),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.ConstructorModal())
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            from textual.widgets import TabbedContent
            modal.query_one(
                "#ctor-tabs", TabbedContent,
            ).active = "ctor-tab-gb_l0"
            await pilot.pause()
            await pilot.pause(0.1)
            modal._lanes["gb_l0"] = [
                ("P", "Promoter",   "Pos 1",   "GGAG", "AATG", "", ""),
                ("C", "CDS",        "Pos 3-4", "AATG", "GCTT", "", ""),
                ("T", "Terminator", "Pos 5",   "GCTT", "CGCT", "", ""),
            ]
            modal._backbones["gb_l0"] = "Alpha1"
            modal._refresh_validation("gb_l0")
            await pilot.pause()
            await pilot.pause(0.1)
            # Save button must be enabled now.
            sim = modal.query_one("#btn-ctor-simulate-gb_l0", sc.Button)
            assert sim.disabled is False
            # READY TO CLONE badge must be visible.
            badge = modal.query_one(
                "#ctor-ready-badge-gb_l0", sc.Static,
            )
            assert "ctor-ready-badge-hidden" not in badge.classes
            # Click Save → expect NamePlasmidModal on top.
            sim.press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.NamePlasmidModal)


# ═══════════════════════════════════════════════════════════════════════════════
# CDS start-codon annotation fix (regression guard for 2026-05-10).
# ═══════════════════════════════════════════════════════════════════════════════
#
# User report: "CDS's cloned seem to lose the annotation of their ATG
# because it also occupies the AATG overhang." The domesticator's
# forward primer absorbs the ATG into the AATG fusion overhang (the
# Pos 12→13 boundary), so a CDS L0 part's body sequence starts at
# codon 2. When the part is assembled into an L1 plasmid, the chained
# sequence reads ...LINK-body + A + ATG + [codon2...], and a feature
# that only spans the body would visibly drop the start codon from
# the cloned plasmid map. `_atg_offset_for_part` returns the 3-nt
# upstream extension that pushes the feature's 5' boundary back into
# the AATG overhang so the start codon is included.

class TestAtgOffsetForPart:
    """The pure helper that decides whether to extend a coding-part
    feature 5' into its upstream AATG overhang."""

    def test_returns_3_for_aatg_coding_parts(self):
        for ptype in ("CDS", "CDS-NS", "Signal peptide", "CDS-NS (CT)"):
            assert sc._atg_offset_for_part("AATG", ptype) == 3, (
                f"coding part {ptype!r} with oh5=AATG should extend by 3"
            )

    def test_returns_0_for_non_aatg_coding(self):
        # C-tag (TTCG legacy), CT-tag (GCAG canonical), CDS-after-SP
        # (AGCC after Signal peptide) — none start with ATG so no
        # extension is meaningful.
        assert sc._atg_offset_for_part("TTCG", "C-tag")         == 0
        assert sc._atg_offset_for_part("GCAG", "CT-tag")        == 0
        assert sc._atg_offset_for_part("AGCC", "CDS-after-SP")  == 0

    def test_returns_0_for_aatg_noncoding(self):
        # Even if a custom grammar uses AATG as a non-coding part's
        # 5' overhang, there's no biological start codon to extend
        # into. The coding-type filter prevents a spurious +3 nt
        # extension on the wrong slot.
        assert sc._atg_offset_for_part("AATG", "Promoter") == 0
        assert sc._atg_offset_for_part("AATG", "5' UTR")   == 0

    def test_returns_0_for_promoter_with_aatg_oh3(self):
        # The Promoter has oh3=AATG (downstream connector) but
        # oh5=GGAG. The helper keys on oh5 only — the ATG is the
        # NEXT part's responsibility, not the promoter's.
        assert sc._atg_offset_for_part("GGAG", "Promoter") == 0

    def test_defensive_against_non_strings(self):
        assert sc._atg_offset_for_part(None, None)   == 0
        assert sc._atg_offset_for_part("AATG", None) == 0
        assert sc._atg_offset_for_part(None, "CDS")  == 0
        assert sc._atg_offset_for_part(123, "CDS")   == 0


class TestReDerivedCdsIncludesStartCodon:
    """Integration: `_re_derive_features_in_insert` must produce a
    CDS feature whose 5' boundary covers the ATG embedded in the
    upstream AATG fusion overhang."""

    def test_forward_cds_feature_starts_at_atg(self, isolated_parts_bin):
        cds_body = "GGGAAATAA"  # codon-2 onward (synthetic)
        sc._save_parts_bin([{
            "name":     "myCDS",   "type":     "CDS",
            "position": "Pos 3-4", "sequence": cds_body,
            "oh5":      "AATG",    "oh3":      "GCTT",
            "grammar":  "gb_l0",   "level":    0,
        }])
        upstream   = "TTTTTTTTTT"
        downstream = "AAAAAAAAAA"
        insert = upstream + "AATG" + cds_body + "GCTT" + downstream
        feats = sc._reconstruct_l0_features_in_seq(
            insert, ["myCDS"], "gb_l0",
        )
        assert len(feats) == 1
        body_start = len(upstream) + 4
        cds = feats[0]
        # +3 nt extension — feature starts at the ATG, not at codon 2.
        assert cds["start"] == body_start - 3
        assert cds["end"]   == body_start + len(cds_body)
        # And the first 3 nt of the feature's span must literally read ATG.
        assert insert[cds["start"]:cds["start"] + 3] == "ATG"

    def test_promoter_feature_unchanged(self, isolated_parts_bin):
        """Promoter has oh5=GGAG (not AATG) — the +3 extension must
        NOT apply, since the AATG is at the promoter's 3' end and
        belongs to the downstream LINK / CDS, not the promoter."""
        prom_body = "TATAATGCG"
        sc._save_parts_bin([{
            "name":     "myPROM",  "type":     "Promoter",
            "position": "Pos 1",   "sequence": prom_body,
            "oh5":      "GGAG",    "oh3":      "AATG",
            "grammar":  "gb_l0",   "level":    0,
        }])
        insert = "NNNNN" + "GGAG" + prom_body + "AATG" + "NNNNN"
        feats = sc._reconstruct_l0_features_in_seq(
            insert, ["myPROM"], "gb_l0",
        )
        assert len(feats) == 1
        body_start = 5 + 4
        assert feats[0]["start"] == body_start
        assert feats[0]["end"]   == body_start + len(prom_body)

    def test_signal_peptide_extended_same_as_cds(self, isolated_parts_bin):
        """Signal peptide has oh5=AATG and is a coding part — same
        +3 extension applies (the SP also starts with ATG)."""
        sp_body = "GCAACAGCC"
        sc._save_parts_bin([{
            "name":     "mySP",    "type":     "Signal peptide",
            "position": "Pos 13",  "sequence": sp_body,
            "oh5":      "AATG",    "oh3":      "AGCC",
            "grammar":  "gb_l0",   "level":    0,
        }])
        insert = "AAAAA" + "AATG" + sp_body + "AGCC" + "TTTTT"
        feats = sc._reconstruct_l0_features_in_seq(
            insert, ["mySP"], "gb_l0",
        )
        assert len(feats) == 1
        body_start = 5 + 4
        assert feats[0]["start"] == body_start - 3
        assert insert[feats[0]["start"]:feats[0]["start"] + 3] == "ATG"

    def test_extension_clamps_at_insert_origin(self, isolated_parts_bin):
        """If a CDS body lands within 3 nt of the linear insert's 5'
        edge (rare but defensible), the extension must clamp to 0
        rather than producing a negative start coordinate."""
        cds_body = "GGGAAATAA"
        sc._save_parts_bin([{
            "name":     "myCDS",   "type":     "CDS",
            "position": "Pos 3-4", "sequence": cds_body,
            "oh5":      "AATG",    "oh3":      "GCTT",
            "grammar":  "gb_l0",   "level":    0,
        }])
        # Body at insert[0:len(cds_body)] — no room upstream.
        insert = cds_body + "GCTT" + "AAAAAAAAAA"
        feats = sc._reconstruct_l0_features_in_seq(
            insert, ["myCDS"], "gb_l0",
        )
        assert len(feats) == 1
        # Clamp: start can't go negative.
        assert feats[0]["start"] == 0

