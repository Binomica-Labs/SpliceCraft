"""
test_features_library — persistent feature library JSON round-trip + API.

The feature library stores user-saved GenBank features that can be inserted
into any plasmid. Entries are dicts like:

    {"name": "lacZ-alpha", "feature_type": "CDS", "sequence": "ATG...",
     "strand": 1, "qualifiers": {"gene": ["lacZ"]}, "description": ""}

These tests cover:
  - `_load_features` / `_save_features` JSON round-trip
  - Corruption recovery via the shared `_safe_save_json` .bak mechanism
  - Non-dict entries filtered on load (hand-edited file safety)
  - `_GENBANK_FEATURE_TYPES` contains the INSDC types SpliceCraft relies on
  - Cache invalidation after save
"""
from __future__ import annotations

import json

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Round-trip
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatureLibraryRoundtrip:
    """Save + reload must yield identical entries for valid inputs."""

    def test_save_creates_file(self):
        sc._save_features([{"name": "x", "feature_type": "CDS",
                            "sequence": "ATG"}])
        assert sc._FEATURES_FILE.exists()

    def test_roundtrip_preserves_entry(self):
        entries = [{"name": "lacZ-alpha", "feature_type": "CDS",
                    "sequence": "ATGACC", "strand": 1,
                    "qualifiers": {"gene": ["lacZ"]},
                    "description": ""}]
        sc._save_features(entries)
        # Bypass the cache to read raw JSON
        raw = json.loads(sc._FEATURES_FILE.read_text())
        assert raw["entries"] == entries

    def test_roundtrip_multiple_entries(self):
        entries = [
            {"name": "p1", "feature_type": "promoter", "sequence": "TATA",
             "strand": 1, "qualifiers": {}, "description": ""},
            {"name": "t1", "feature_type": "terminator", "sequence": "TTT",
             "strand": 1, "qualifiers": {"note": ["rho-independent"]},
             "description": ""},
        ]
        sc._save_features(entries)
        sc._features_cache = None  # force reload from disk
        reloaded = sc._load_features()
        assert len(reloaded) == 2
        assert reloaded[0]["name"] == "p1"
        assert reloaded[1]["qualifiers"]["note"] == ["rho-independent"]

    def test_envelope_schema_version(self):
        """Features file uses the shared schema envelope (sacred invariant #7)."""
        sc._save_features([{"name": "x", "feature_type": "CDS",
                            "sequence": "ATG"}])
        raw = json.loads(sc._FEATURES_FILE.read_text())
        assert raw["_schema_version"] == sc._CURRENT_SCHEMA_VERSION
        assert isinstance(raw["entries"], list)

    def test_save_creates_bak_on_overwrite(self):
        sc._save_features([{"name": "first", "feature_type": "CDS",
                            "sequence": "A"}])
        sc._save_features([{"name": "second", "feature_type": "CDS",
                            "sequence": "T"}])
        bak_path = sc._FEATURES_FILE.with_suffix(sc._FEATURES_FILE.suffix + ".bak")
        assert bak_path.exists()
        assert json.loads(bak_path.read_text())["entries"][0]["name"] == "first"


# ═══════════════════════════════════════════════════════════════════════════════
# Corruption recovery
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatureLibraryCorruptionRecovery:
    """Corrupt / missing / hand-edited files must not crash the loader."""

    def test_missing_file_returns_empty(self):
        """First run: file doesn't exist → empty list, no error."""
        sc._features_cache = None
        assert sc._load_features() == []

    def test_corrupt_json_returns_empty(self):
        sc._FEATURES_FILE.write_text("{bad json")
        sc._features_cache = None
        # No valid main or bak → empty list
        assert sc._load_features() == []

    def test_non_dict_entries_filtered(self):
        """A hand-edited file with garbage entries must not crash `.get()` callers."""
        sc._FEATURES_FILE.write_text(json.dumps({
            "_schema_version": 1,
            "entries": [
                {"name": "good", "feature_type": "CDS", "sequence": "ATG"},
                "not a dict",
                42,
                None,
                {"name": "also good", "feature_type": "gene", "sequence": "T"},
            ],
        }))
        sc._features_cache = None
        entries = sc._load_features()
        assert len(entries) == 2
        assert entries[0]["name"] == "good"
        assert entries[1]["name"] == "also good"

    def test_bak_restore_after_main_corruption(self):
        """If main is corrupt but a .bak exists, the .bak is restored."""
        # First, write a valid file (creates .bak on next write)
        sc._save_features([{"name": "first", "feature_type": "CDS",
                            "sequence": "ATG"}])
        sc._save_features([{"name": "second", "feature_type": "CDS",
                            "sequence": "TAA"}])
        # Now corrupt main; .bak holds the 'first' version
        sc._FEATURES_FILE.write_text("!!!corrupt!!!")
        sc._features_cache = None
        entries = sc._load_features()
        assert len(entries) == 1
        assert entries[0]["name"] == "first"


# ═══════════════════════════════════════════════════════════════════════════════
# Cache behaviour
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatureLibraryCache:
    """`_features_cache` must stay in sync with the on-disk state."""

    def test_save_updates_cache(self):
        sc._save_features([{"name": "x", "feature_type": "CDS",
                            "sequence": "ATG"}])
        # Next load should return the saved entries without hitting disk
        sc._FEATURES_FILE.unlink()  # disk gone, but cache populated
        assert sc._load_features() == [{"name": "x", "feature_type": "CDS",
                                        "sequence": "ATG"}]

    def test_load_returns_copy_not_reference(self):
        """Mutating the returned list must not poison the cache."""
        sc._save_features([{"name": "x", "feature_type": "CDS",
                            "sequence": "A"}])
        loaded = sc._load_features()
        loaded.append({"name": "SHOULD_NOT_PERSIST"})
        loaded2 = sc._load_features()
        assert len(loaded2) == 1


# ═══════════════════════════════════════════════════════════════════════════════
# Curated type list
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenbankFeatureTypes:
    """`_GENBANK_FEATURE_TYPES` is the dropdown source for the Add Feature
    modal. It must contain the INSDC types SpliceCraft relies on."""

    def test_contains_core_types(self):
        core = {"CDS", "gene", "promoter", "terminator", "RBS", "5'UTR",
                "3'UTR", "intron", "exon", "rep_origin", "misc_feature",
                "primer_bind"}
        assert core.issubset(set(sc._GENBANK_FEATURE_TYPES))

    def test_does_not_include_source(self):
        """`source` is excluded: each GenBank record already has exactly one
        `source` feature spanning the whole molecule. Adding another would
        be invalid. (Regression guard.)"""
        assert "source" not in sc._GENBANK_FEATURE_TYPES

    def test_all_types_are_strings(self):
        for t in sc._GENBANK_FEATURE_TYPES:
            assert isinstance(t, str)
            assert t  # non-empty

    def test_no_duplicates(self):
        assert len(set(sc._GENBANK_FEATURE_TYPES)) == len(sc._GENBANK_FEATURE_TYPES)
