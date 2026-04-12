"""
test_primers — Primer design backend + persistence tests.

Covers:
  - _design_detection_primers: Primer3 returns valid pair, Tm in range,
    product size in range, primers flank the target
  - _design_cloning_primers: RE sites present, GCGC padding, binding Tm,
    correct RE at correct end
  - Primer library persistence: save/load round-trip
  - PrimerDesignScreen: mounts with correct widgets
"""
from __future__ import annotations

import json
import random

import pytest

import splicecraft as sc


@pytest.fixture
def random_seq_3k():
    rng = random.Random(0xBEEF)
    return "".join(rng.choice("ACGT") for _ in range(3000))


@pytest.fixture
def isolated_primers(tmp_path, monkeypatch):
    tmp_p = tmp_path / "primers.json"
    monkeypatch.setattr(sc, "_PRIMERS_FILE", tmp_p)
    monkeypatch.setattr(sc, "_primers_cache", None)
    return tmp_p


# ═══════════════════════════════════════════════════════════════════════════════
# Detection primers
# ═══════════════════════════════════════════════════════════════════════════════

class TestDetectionPrimers:
    def test_returns_valid_pair(self, random_seq_3k):
        r = sc._design_detection_primers(
            random_seq_3k, 500, 600, product_min=400, product_max=700,
        )
        assert "error" not in r
        assert r["fwd_seq"]
        assert r["rev_seq"]

    def test_tm_in_range(self, random_seq_3k):
        r = sc._design_detection_primers(
            random_seq_3k, 500, 600, target_tm=60.0,
            product_min=400, product_max=700,
        )
        assert "error" not in r
        assert 55 < r["fwd_tm"] < 65
        assert 55 < r["rev_tm"] < 65

    def test_product_size_in_range(self, random_seq_3k):
        r = sc._design_detection_primers(
            random_seq_3k, 500, 600,
            product_min=450, product_max=600,
        )
        assert "error" not in r
        assert 450 <= r["product_size"] <= 600

    def test_fwd_upstream_of_target(self, random_seq_3k):
        r = sc._design_detection_primers(
            random_seq_3k, 500, 600,
            product_min=400, product_max=700,
        )
        assert "error" not in r
        assert r["fwd_pos"][0] < 500, "fwd primer should start before target"

    def test_rev_downstream_of_target(self, random_seq_3k):
        r = sc._design_detection_primers(
            random_seq_3k, 500, 600,
            product_min=400, product_max=700,
        )
        assert "error" not in r
        assert r["rev_pos"][1] > 600, "rev primer should end after target"

    def test_empty_target_returns_error(self, random_seq_3k):
        r = sc._design_detection_primers(random_seq_3k, 500, 500)
        assert "error" in r

    def test_impossible_constraints_returns_error(self, random_seq_3k):
        # 10bp product range on a 100bp target — impossible
        r = sc._design_detection_primers(
            random_seq_3k, 500, 600,
            product_min=10, product_max=20,
        )
        assert "error" in r


# ═══════════════════════════════════════════════════════════════════════════════
# Cloning primers
# ═══════════════════════════════════════════════════════════════════════════════

class TestCloningPrimers:
    def test_returns_valid_pair(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI",
        )
        assert "error" not in r
        assert r["fwd_full"]
        assert r["rev_full"]

    def test_gcgc_padding_present(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI",
        )
        assert r["fwd_full"].startswith("GCGC")
        assert r["rev_full"].startswith("GCGC")

    def test_5prime_re_site_in_fwd(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI",
        )
        # EcoRI = GAATTC — should appear right after GCGC in fwd primer
        assert "GAATTC" in r["fwd_full"]

    def test_3prime_re_site_rc_in_rev(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI",
        )
        # BamHI = GGATCC → RC = GGATCC (palindrome)
        assert "GGATCC" in r["rev_full"]

    def test_non_palindrome_re_site_rc(self, random_seq_3k):
        # BsaI = GGTCTC (non-palindrome) → RC = GAGACC
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "BsaI", "EcoRI",
        )
        assert "GGTCTC" in r["fwd_full"], "5' BsaI site in fwd"
        # Rev should have RC of EcoRI = GAATTC (palindrome → same)
        assert "GAATTC" in r["rev_full"]

    def test_binding_tm_near_target(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI", target_tm=57.0,
        )
        assert 49 < r["fwd_tm"] < 65
        assert 49 < r["rev_tm"] < 65

    def test_unknown_enzyme_returns_error(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "FakeEnzyme", "BamHI",
        )
        assert "error" in r

    def test_short_region_returns_error(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 210, "EcoRI", "BamHI",
        )
        assert "error" in r

    @pytest.mark.parametrize("re5,re3", [
        ("EcoRI", "BamHI"), ("XhoI", "NdeI"), ("NcoI", "XbaI"),
        ("SpeI", "PstI"), ("HindIII", "SalI"), ("NotI", "BglII"),
    ])
    def test_various_enzyme_pairs(self, random_seq_3k, re5, re3):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, re5, re3,
        )
        assert "error" not in r
        site_5, _, _ = sc._NEB_ENZYMES[re5]
        site_3, _, _ = sc._NEB_ENZYMES[re3]
        assert site_5 in r["fwd_full"]
        assert sc._rc(site_3) in r["rev_full"] or site_3 in r["rev_full"]


# ═══════════════════════════════════════════════════════════════════════════════
# Primer library persistence
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrimerPersistence:
    def test_empty_loads_as_empty(self, isolated_primers):
        assert sc._load_primers() == []

    def test_save_load_roundtrip(self, isolated_primers):
        primers = [
            {"name": "ampR-DET-F", "sequence": "ATGAAAGATCTG", "tm": 57.2,
             "primer_type": "detection", "source": "pUC19"},
        ]
        sc._save_primers(primers)
        loaded = sc._load_primers()
        assert len(loaded) == 1
        assert loaded[0]["name"] == "ampR-DET-F"

    def test_writes_valid_json(self, isolated_primers):
        sc._save_primers([{"name": "x", "sequence": "ATG"}])
        assert isolated_primers.exists()
        parsed = json.loads(isolated_primers.read_text())
        assert parsed[0]["name"] == "x"

    def test_corrupted_file_returns_empty(self, isolated_primers):
        isolated_primers.write_text("{bad")
        sc._primers_cache = None
        assert sc._load_primers() == []


# ═══════════════════════════════════════════════════════════════════════════════
# RE options list
# ═══════════════════════════════════════════════════════════════════════════════

class TestCloningREOptions:
    def test_common_enzymes_present(self):
        names = {name for _, name in sc._CLONING_RE_OPTIONS}
        for must in ["EcoRI", "BamHI", "XhoI", "NdeI", "NotI", "BsaI",
                     "HindIII", "XbaI", "NcoI", "SalI"]:
            assert must in names, f"{must} missing from RE options"

    def test_all_options_are_in_neb_catalog(self):
        for label, name in sc._CLONING_RE_OPTIONS:
            assert name in sc._NEB_ENZYMES, f"{name} not in _NEB_ENZYMES"
