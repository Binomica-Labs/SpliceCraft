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


# ═══════════════════════════════════════════════════════════════════════════════
# Color + strand extensions (v0.3.2)
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatureEntryColorField:
    """Entries gained an optional ``color`` field and accept ``strand=0``
    for arrowless features. Round-trip must preserve both verbatim."""

    def test_color_field_roundtrip(self):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1, "color": "#FF6347",
        }])
        sc._features_cache = None
        loaded = sc._load_features()
        assert loaded[0]["color"] == "#FF6347"

    def test_strand_zero_roundtrip(self):
        """strand=0 means "arrowless" — meaningful for rep_origin,
        misc_feature, and similar non-directional annotations."""
        sc._save_features([{
            "name": "pMB1", "feature_type": "rep_origin",
            "sequence": "GCA", "strand": 0,
        }])
        sc._features_cache = None
        loaded = sc._load_features()
        assert loaded[0]["strand"] == 0

    def test_missing_color_survives(self):
        """Legacy entries without a color field must still load."""
        sc._save_features([{
            "name": "legacy", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        sc._features_cache = None
        loaded = sc._load_features()
        assert "color" not in loaded[0] or loaded[0].get("color") is None


class TestFeatureColorsPersistence:
    """``_load_feature_colors`` / ``_save_feature_colors`` manage the user-
    editable type → default color map. Empty, missing, and corrupt files
    all degrade to {}."""

    def test_missing_file_returns_empty(self):
        sc._feature_colors_cache = None
        assert sc._load_feature_colors() == {}

    def test_save_and_reload(self):
        sc._save_feature_colors({"CDS": "#FF0000", "promoter": "#00FF00"})
        sc._feature_colors_cache = None
        loaded = sc._load_feature_colors()
        assert loaded == {"CDS": "#FF0000", "promoter": "#00FF00"}

    def test_envelope_schema_version(self):
        sc._save_feature_colors({"CDS": "#FF0000"})
        raw = json.loads(sc._FEATURE_COLORS_FILE.read_text())
        assert raw["_schema_version"] == sc._CURRENT_SCHEMA_VERSION
        assert isinstance(raw["entries"], list)
        assert raw["entries"][0] == {"feature_type": "CDS", "color": "#FF0000"}

    def test_non_dict_entries_filtered(self):
        sc._FEATURE_COLORS_FILE.write_text(json.dumps({
            "_schema_version": 1,
            "entries": [
                {"feature_type": "CDS",      "color": "#FF0000"},
                "garbage", 42, None,
                {"feature_type": "",         "color": "#FF0000"},  # empty key
                {"feature_type": "promoter", "color":  ""},        # empty color
                {"feature_type": "gene",     "color": "#00FF00"},
            ],
        }))
        sc._feature_colors_cache = None
        loaded = sc._load_feature_colors()
        assert loaded == {"CDS": "#FF0000", "gene": "#00FF00"}


class TestResolveFeatureColor:
    """``_resolve_feature_color`` enforces the precedence:
    entry color > user default > built-in default > palette fallback."""

    def test_entry_color_wins(self):
        sc._save_feature_colors({"CDS": "#000000"})
        sc._feature_colors_cache = None
        col = sc._resolve_feature_color(
            {"feature_type": "CDS", "color": "#AAAAAA"}
        )
        assert col == "#AAAAAA"

    def test_user_default_over_builtin(self):
        sc._save_feature_colors({"CDS": "#000000"})
        sc._feature_colors_cache = None
        col = sc._resolve_feature_color({"feature_type": "CDS"})
        assert col == "#000000"

    def test_builtin_default_when_no_user_override(self):
        sc._feature_colors_cache = None
        col = sc._resolve_feature_color({"feature_type": "CDS"})
        assert col == sc._DEFAULT_TYPE_COLORS["CDS"]

    def test_unknown_type_falls_back_to_palette(self):
        """Palette fallback returns the hex equivalent of
        ``_FEATURE_PALETTE[0]``. The 2026-04-20 ColorPicker rework
        normalises ``color(N)`` palette syntax to hex so downstream Rich
        markup never trips on the parens — the *color* is still the same
        palette entry, just expressed in a markup-safe form."""
        sc._feature_colors_cache = None
        col = sc._resolve_feature_color({"feature_type": "my_custom_type"})
        assert col == sc._normalise_color_input(sc._FEATURE_PALETTE[0])

    def test_empty_color_treated_as_missing(self):
        """Empty string color → treat as not set, fall through to type default."""
        sc._feature_colors_cache = None
        col = sc._resolve_feature_color(
            {"feature_type": "CDS", "color": ""}
        )
        assert col == sc._DEFAULT_TYPE_COLORS["CDS"]

    def test_all_default_types_have_hex_colors(self):
        """Every built-in default is a valid hex string so Rich can render it."""
        for ftype, col in sc._DEFAULT_TYPE_COLORS.items():
            assert isinstance(col, str)
            assert col.startswith("#")
            assert len(col) == 7  # "#RRGGBB"

    def test_all_insdc_types_have_builtin_color(self):
        """Every curated _GENBANK_FEATURE_TYPES entry has a built-in default
        so the UI never falls back to the palette for standard types."""
        missing = [t for t in sc._GENBANK_FEATURE_TYPES
                   if t not in sc._DEFAULT_TYPE_COLORS]
        assert not missing, f"No default color for: {missing}"
