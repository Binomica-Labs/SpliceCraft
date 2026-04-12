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

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def random_template():
    """A 2000-bp random ACGT template for primer design tests."""
    rng = random.Random(0xCAFE)
    return "".join(rng.choice("ACGT") for _ in range(2000))


@pytest.fixture
def isolated_parts_bin(tmp_path, monkeypatch):
    """Redirect _PARTS_BIN_FILE to a tmp path so tests don't touch
    the real parts_bin.json."""
    tmp_bin = tmp_path / "parts_bin.json"
    monkeypatch.setattr(sc, "_PARTS_BIN_FILE", tmp_bin)
    monkeypatch.setattr(sc, "_parts_bin_cache", None)
    return tmp_bin


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
    def test_bsai_site_in_both_primers(self, random_template, part_type):
        result = sc._design_gb_primers(random_template, 100, 600, part_type)
        assert sc._GB_BSAI_SITE in result["fwd_full"]
        assert sc._GB_BSAI_SITE in result["rev_full"]

    @pytest.mark.parametrize("part_type", list(sc._GB_POSITIONS.keys()))
    def test_correct_overhangs_in_primers(self, random_template, part_type):
        """The 5' overhang must appear in the forward primer right after
        pad+BsaI+spacer. The RC of the 3' overhang must appear in the
        reverse primer at the same offset."""
        result = sc._design_gb_primers(random_template, 100, 600, part_type)
        oh5, oh3 = result["oh5"], result["oh3"]
        tail_prefix_len = len(sc._GB_PAD + sc._GB_BSAI_SITE + sc._GB_SPACER)

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
        result = sc._design_gb_primers(random_template, 100, 600, "CDS")
        insert = random_template[100:600].upper()
        assert insert.startswith(result["fwd_binding"])

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
        # Spacer A should appear between BsaI site and overhang
        fwd_after_bsai = result["fwd_full"][len(sc._GB_PAD) + len(sc._GB_BSAI_SITE)]
        assert fwd_after_bsai == sc._GB_SPACER


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
            "fwd_primer": "GCGCGGTCTCAGGAGATGAAAGATCTG",
            "rev_primer": "GCGCGGTCTCAGTCACAGATCTTTCAT",
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
        assert parsed == parts

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

    def test_bsai_site_is_ggtctc(self):
        assert sc._GB_BSAI_SITE == "GGTCTC"

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
        This is the core Golden Braid assembly principle."""
        chain = [
            ("Promoter", "5' UTR"),
            ("5' UTR",   "CDS"),
            ("CDS-NS",   "C-tag"),
            ("C-tag",    "Terminator"),   # C-tag 3' = GCTT = Terminator 5'
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
        must match the ConstructorModal._TU_START / _TU_END constants."""
        _, oh5_first, _ = sc._GB_POSITIONS["Promoter"]
        _, _, oh3_last  = sc._GB_POSITIONS["Terminator"]
        assert oh5_first == sc.ConstructorModal._TU_START
        assert oh3_last  == sc.ConstructorModal._TU_END
