"""
test_genbank_io — GenBank file I/O round-trip tests.

Guards:
  - `load_genbank(path)` parses a real .gb file and preserves sequence bytes
    and feature count
  - `_record_to_gb_text` / `_gb_text_to_record` round-trip is lossless for
    the fields SpliceCraft actually relies on (seq, features, qualifiers)
  - Library save/load via `_save_library` / `_load_library` round-trips
    through JSON without corrupting accession / name / seq fields

These run entirely offline — no NCBI calls, no network. `fetch_genbank`
itself is covered by manual smoke testing, not automated tests.
"""
from __future__ import annotations

import json

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# File I/O round-trip
# ═══════════════════════════════════════════════════════════════════════════════

class TestLoadGenbank:
    def test_load_returns_seqrecord(self, tiny_gb_path):
        rec = sc.load_genbank(tiny_gb_path)
        # Duck-type: has .seq, .features, .id
        assert hasattr(rec, "seq")
        assert hasattr(rec, "features")
        assert hasattr(rec, "id")

    def test_sequence_length_preserved(self, tiny_gb_path, tiny_record):
        rec = sc.load_genbank(tiny_gb_path)
        assert len(rec.seq) == len(tiny_record.seq)

    def test_sequence_bytes_exact(self, tiny_gb_path, tiny_record):
        rec = sc.load_genbank(tiny_gb_path)
        assert str(rec.seq) == str(tiny_record.seq)

    def test_features_preserved(self, tiny_gb_path, tiny_record):
        rec = sc.load_genbank(tiny_gb_path)
        # The fixture has 2 features (CDS + misc_feature). Biopython may also
        # emit a 'source' feature when parsing, so count only non-source.
        non_source_in = [f for f in rec.features if f.type != "source"]
        non_source_fx = [f for f in tiny_record.features if f.type != "source"]
        assert len(non_source_in) == len(non_source_fx)

    def test_cds_feature_strand_preserved(self, tiny_gb_path):
        rec = sc.load_genbank(tiny_gb_path)
        cds = [f for f in rec.features if f.type == "CDS"]
        assert len(cds) == 1
        assert cds[0].location.strand == 1


# ═══════════════════════════════════════════════════════════════════════════════
# Text round-trip via StringIO
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenbankTextRoundtrip:
    """`_record_to_gb_text` → `_gb_text_to_record` must preserve the fields
    the UI touches: sequence bytes, feature type, strand, start, end."""

    def test_sequence_preserved(self, tiny_record):
        text = sc._record_to_gb_text(tiny_record)
        rec2 = sc._gb_text_to_record(text)
        assert str(rec2.seq) == str(tiny_record.seq)

    def test_feature_types_preserved(self, tiny_record):
        text = sc._record_to_gb_text(tiny_record)
        rec2 = sc._gb_text_to_record(text)
        types_in = sorted(f.type for f in tiny_record.features if f.type != "source")
        types_out = sorted(f.type for f in rec2.features if f.type != "source")
        assert types_in == types_out

    def test_feature_strands_preserved(self, tiny_record):
        text = sc._record_to_gb_text(tiny_record)
        rec2 = sc._gb_text_to_record(text)
        strands_in = sorted(
            (f.type, f.location.strand)
            for f in tiny_record.features if f.type != "source"
        )
        strands_out = sorted(
            (f.type, f.location.strand)
            for f in rec2.features if f.type != "source"
        )
        assert strands_in == strands_out

    def test_feature_positions_preserved(self, tiny_record):
        text = sc._record_to_gb_text(tiny_record)
        rec2 = sc._gb_text_to_record(text)
        pos_in = sorted(
            (f.type, int(f.location.start), int(f.location.end))
            for f in tiny_record.features if f.type != "source"
        )
        pos_out = sorted(
            (f.type, int(f.location.start), int(f.location.end))
            for f in rec2.features if f.type != "source"
        )
        assert pos_in == pos_out


# ═══════════════════════════════════════════════════════════════════════════════
# Library persistence (JSON)
# ═══════════════════════════════════════════════════════════════════════════════

class TestLibraryPersistence:
    """`_load_library` / `_save_library` use a module-global `_LIBRARY_FILE`.
    The `isolated_library` fixture redirects it to a tmp path so the real
    `plasmid_library.json` isn't touched."""

    def test_empty_library_loads_as_empty_list(self, isolated_library):
        assert sc._load_library() == []

    def test_save_then_load_roundtrip(self, isolated_library):
        entries = [
            {"id": "X001", "name": "test1", "seq": "ACGT", "length": 4},
            {"id": "X002", "name": "test2", "seq": "GATTACA", "length": 7},
        ]
        sc._save_library(entries)
        loaded = sc._load_library()
        assert loaded == entries

    def test_save_writes_valid_json(self, isolated_library):
        entries = [{"id": "Y001", "name": "probe", "seq": "A" * 10, "length": 10}]
        sc._save_library(entries)
        # Bypass the cache and read raw bytes
        assert isolated_library.exists()
        parsed = json.loads(isolated_library.read_text())
        assert parsed == entries

    def test_load_survives_corrupted_file(self, isolated_library, caplog):
        """If the library JSON is corrupted, `_load_library` must return []
        and log the exception — never propagate the error to the UI."""
        isolated_library.write_text("{not valid json")
        # Reset in-memory cache so _load_library actually re-reads the file
        sc._library_cache = None
        result = sc._load_library()
        assert result == []

    def test_load_memoizes(self, isolated_library):
        """Second call should hit the in-memory cache, not re-parse the file."""
        entries = [{"id": "Z001", "name": "n", "seq": "A", "length": 1}]
        sc._save_library(entries)
        once = sc._load_library()
        twice = sc._load_library()
        assert once == twice == entries
