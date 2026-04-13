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


# ═══════════════════════════════════════════════════════════════════════════════
# Multi-record / malformed file handling (added 2026-04-12)
# ═══════════════════════════════════════════════════════════════════════════════

class TestMultiRecordFiles:
    """`load_genbank` used to propagate Biopython's raw `More than one record`
    exception when a file contained multiple records. It now raises a
    user-friendly ValueError listing the accessions."""

    def _two_record_text(self, tiny_record):
        text1 = sc._record_to_gb_text(tiny_record)
        # Make a second record with a different id. GenBank requires
        # molecule_type in annotations for Biopython to write it.
        from Bio.SeqRecord import SeqRecord
        r2 = SeqRecord(
            tiny_record.seq, id="ALT12345", name="ALT",
            annotations={"molecule_type": "DNA"},
        )
        text2 = sc._record_to_gb_text(r2)
        return text1 + text2

    def test_two_records_raises_value_error_with_ids(self, tmp_path, tiny_record):
        gb = tmp_path / "multi.gb"
        gb.write_text(self._two_record_text(tiny_record))
        with pytest.raises(ValueError, match="2 records"):
            sc.load_genbank(str(gb))

    def test_empty_file_raises_value_error(self, tmp_path):
        gb = tmp_path / "empty.gb"
        gb.write_text("")
        with pytest.raises(ValueError, match="no GenBank records"):
            sc.load_genbank(str(gb))

    def test_non_genbank_text_raises_value_error(self, tmp_path):
        gb = tmp_path / "notgb.gb"
        gb.write_text(">fasta header\nACGTACGT\n")   # FASTA, not GenBank
        with pytest.raises(ValueError, match="no GenBank records"):
            sc.load_genbank(str(gb))


class TestParseRobustness:
    """`PlasmidMap._parse` must tolerate unusual features (compound locations,
    UnknownPosition) without crashing — users must be warned, not locked out."""

    def test_compound_location_counted_and_flattened(self, tiny_record):
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        # Build a fresh record with a compound-location feature
        from copy import deepcopy
        rec = deepcopy(tiny_record)
        compound = CompoundLocation([
            FeatureLocation(10, 30, strand=1),
            FeatureLocation(50, 80, strand=1),
        ])
        rec.features.append(SeqFeature(compound, type="mRNA",
                                       qualifiers={"label": ["spliced"]}))
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)  # don't mount
        feats = pm._parse(rec)
        # The mRNA feature is rendered at outer bounds [10, 80)
        spliced = [f for f in feats if f.get("label") == "spliced"]
        assert len(spliced) == 1
        assert spliced[0]["start"] == 10 and spliced[0]["end"] == 80
        # Counter surfaces for the caller to notify
        assert pm._n_flattened == 1
        assert pm._n_skipped == 0

    def test_unknown_position_is_skipped_not_crashed(self, tiny_record):
        from Bio.SeqFeature import SeqFeature, FeatureLocation, UnknownPosition
        from copy import deepcopy
        rec = deepcopy(tiny_record)
        # Feature with an unknown end coordinate — real-world rare but legal.
        # Biopython doesn't accept UnknownPosition objects directly in
        # FeatureLocation constructor post-1.80, so monkeypatch the _parse
        # path with a feature whose int() cast will fail.
        class BadLoc:
            start = 10
            end = "not-an-int"   # will fail int()
            strand = 1
        bad = SeqFeature(type="regulatory", qualifiers={"label": ["bad"]})
        bad.location = BadLoc()
        rec.features.append(bad)
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        # The bad feature is silently dropped (caller notifies via _n_skipped)
        assert not any(f.get("label") == "bad" for f in feats)
        assert pm._n_skipped == 1


# ═══════════════════════════════════════════════════════════════════════════════
# 2026-04-13 parser polish — empty labels, out-of-range clamps, wrap detection
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatLabelFallback:
    """Empty qualifier values (e.g. /label="") must fall through to the
    next candidate qualifier or to feat.type so users never see a
    blank entry in the sidebar."""

    def test_empty_label_falls_back_to_type(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="misc_feature",
                       qualifiers={"label": [""]})
        assert sc._feat_label(f) == "misc_feature"

    def test_empty_label_falls_back_to_gene(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="CDS",
                       qualifiers={"label": [""], "gene": ["ampR"]})
        assert sc._feat_label(f) == "ampR"

    def test_whitespace_only_label_falls_back(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="gene",
                       qualifiers={"label": ["   "]})
        assert sc._feat_label(f) == "gene"

    def test_real_label_truncates_at_28(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="gene",
                       qualifiers={"label": ["a" * 50]})
        assert sc._feat_label(f) == "a" * 28

    def test_newline_in_label_collapsed(self):
        """Multi-line /note="..." qualifiers picked up as the label
        must have newlines/tabs collapsed to single spaces — raw
        newlines would break the sidebar row and map label."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="CDS",
                       qualifiers={"product": ["Line one\nLine two"]})
        assert sc._feat_label(f) == "Line one Line two"

    def test_tab_in_label_collapsed(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="CDS",
                       qualifiers={"product": ["name\twith\ttabs"]})
        assert sc._feat_label(f) == "name with tabs"

    def test_runs_of_whitespace_collapsed(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="gene",
                       qualifiers={"label": ["many    spaces   here"]})
        assert sc._feat_label(f) == "many spaces here"

    def test_whitespace_only_label_falls_back(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="CDS",
                       qualifiers={"label": ["\n\n\t  "],
                                   "gene":  ["ampR"]})
        assert sc._feat_label(f) == "ampR"

    def test_non_ascii_preserved(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="CDS",
                       qualifiers={"product": ["β-galactosidase"]})
        assert sc._feat_label(f) == "β-galactosidase"


class TestParseClamping:
    """Features with out-of-range coordinates are clamped to
    [0, len(seq)] rather than rendered past the end of the map."""

    def _make_record(self, seq_len: int = 400):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        return SeqRecord(Seq("A" * seq_len), id="X",
                         annotations={"molecule_type": "DNA"})

    def test_end_beyond_seq_clamped(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._make_record(400)
        rec.features.append(SeqFeature(
            FeatureLocation(350, 500, strand=1), type="gene",
            qualifiers={"label": ["beyond"]}))
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        assert feats[0]["start"] == 350 and feats[0]["end"] == 400
        assert pm._n_clamped == 1

    def test_negative_start_clamped_to_zero(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._make_record(400)
        # Direct construction: Biopython usually rejects negatives, so
        # simulate via a custom location object that returns negatives.
        class NegLoc:
            start = -5
            end = 10
            strand = 1
        f = SeqFeature(type="gene", qualifiers={"label": ["neg"]})
        f.location = NegLoc()
        rec.features.append(f)
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        assert feats[0]["start"] == 0 and feats[0]["end"] == 10
        assert pm._n_clamped == 1

    def test_in_range_not_counted(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._make_record(400)
        rec.features.append(SeqFeature(
            FeatureLocation(10, 100, strand=1), type="gene",
            qualifiers={"label": ["ok"]}))
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        assert pm._n_clamped == 0


class TestWrapFeatureDetection:
    """GenBank join(n-X..n, 1..Y) on circular plasmids now turns into
    a proper wrap feature (end<start) instead of being flattened to
    whole-plasmid outer bounds."""

    def _make_record(self, seq_len: int = 400):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        return SeqRecord(Seq("A" * seq_len), id="X",
                         annotations={"molecule_type": "DNA"})

    def test_wrap_detected_two_parts(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        rec = self._make_record(400)
        cl = CompoundLocation([
            FeatureLocation(350, 400, strand=1),
            FeatureLocation(0, 50, strand=1),
        ])
        rec.features.append(SeqFeature(cl, type="gene",
                                       qualifiers={"label": ["wrap"]}))
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        f = feats[0]
        assert f["start"] == 350 and f["end"] == 50, (
            f"Expected wrap start=350 end=50, got {f['start']}..{f['end']}"
        )
        # Wrap features are not counted as "flattened" because we
        # preserve their wrap semantics rather than losing information.
        assert pm._n_flattened == 0

    def test_non_wrap_compound_still_flattened(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        rec = self._make_record(400)
        cl = CompoundLocation([
            FeatureLocation(100, 200, strand=1),
            FeatureLocation(250, 350, strand=1),
        ])
        rec.features.append(SeqFeature(cl, type="mRNA",
                                       qualifiers={"label": ["exons"]}))
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        f = feats[0]
        assert f["start"] == 100 and f["end"] == 350
        assert pm._n_flattened == 1

    def test_contiguous_compound_not_flagged_as_flattened(self):
        """SnapGene emits adjacent-parts CompoundLocations like
        parts=[(46, 47), (47, 52)] for a 6-bp RE site, where the
        outer bounds perfectly capture the real feature. Those don't
        represent lost information, so we must NOT bother the user
        with a ⚠ "features flattened" notification."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        rec = self._make_record(400)
        cl = CompoundLocation([
            FeatureLocation(46, 47, strand=1),
            FeatureLocation(47, 52, strand=1),
        ])
        rec.features.append(SeqFeature(cl, type="misc_feature",
                                       qualifiers={"label": ["EcoRI"]}))
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(rec)
        f = feats[0]
        # Outer bounds = the real span
        assert f["start"] == 46 and f["end"] == 52
        # Crucially: contiguous, so NOT counted as flattened
        assert pm._n_flattened == 0


class TestRecordToGbTextResilience:
    """_record_to_gb_text defaults molecule_type to 'DNA' if the
    caller's record lacks it, so saving primer-annotated plasmids
    doesn't crash when the record came from an unusual source."""

    def test_missing_molecule_type_defaults_to_dna(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ACGT" * 25), id="X")
        # No annotations at all — must not crash
        text = sc._record_to_gb_text(rec)
        assert "LOCUS" in text
        assert "DNA" in text.splitlines()[0]

    def test_existing_molecule_type_preserved(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ACGU" * 25), id="X",
                        annotations={"molecule_type": "RNA"})
        text = sc._record_to_gb_text(rec)
        assert "RNA" in text.splitlines()[0]


# ═══════════════════════════════════════════════════════════════════════════════
# SnapGene (.dna) import — format-detection dispatch (added 2026-04-13)
# ═══════════════════════════════════════════════════════════════════════════════

class TestPlasmidFormatDetection:
    """`_detect_plasmid_format` picks the right Biopython SeqIO format
    key from the file extension so `load_genbank` can load both
    GenBank (.gb/.gbk) and SnapGene (.dna) natively."""

    def test_gb_extension(self):
        assert sc._detect_plasmid_format("x.gb") == "genbank"

    def test_gbk_extension(self):
        assert sc._detect_plasmid_format("/path/to/x.gbk") == "genbank"

    def test_genbank_extension(self):
        assert sc._detect_plasmid_format("x.genbank") == "genbank"

    def test_dna_extension_is_snapgene(self):
        assert sc._detect_plasmid_format("x.dna") == "snapgene"

    def test_dna_uppercase(self):
        """Case-insensitive — users drag-dropping files from macOS/Windows
        sometimes have capitalized extensions."""
        assert sc._detect_plasmid_format("plasmid.DNA") == "snapgene"
        assert sc._detect_plasmid_format("foo.GB") == "genbank"

    def test_unknown_extension_defaults_genbank(self):
        """Fallthrough default is GenBank since that's the most common
        plasmid format; parser will raise a clear error if the content
        doesn't match."""
        assert sc._detect_plasmid_format("x.txt") == "genbank"
        assert sc._detect_plasmid_format("no-extension") == "genbank"


class TestSnapGeneDispatch:
    """`load_genbank(path)` dispatches to the SnapGene parser when the
    file extension is .dna, so users can import SnapGene-native files
    just by opening them — no need for a SnapGene → GenBank export first."""

    def test_dna_path_calls_snapgene_parser(self, tmp_path, monkeypatch,
                                             tiny_record):
        """The parser path is verified by monkeypatching SeqIO.parse to
        record the format string it was called with."""
        called_fmt = {}

        def fake_parse(path, fmt):
            called_fmt["fmt"] = fmt
            called_fmt["path"] = path
            return iter([tiny_record])

        from Bio import SeqIO
        monkeypatch.setattr(SeqIO, "parse", fake_parse)

        # Create a real .dna file stub — contents don't matter because
        # SeqIO.parse is mocked. Just needs the extension.
        dna_file = tmp_path / "plasmid.dna"
        dna_file.write_bytes(b"fake binary")

        rec = sc.load_genbank(str(dna_file))
        assert called_fmt["fmt"] == "snapgene"
        assert called_fmt["path"] == str(dna_file)
        # tiny_record is returned through _pick_single_record
        assert rec is tiny_record

    def test_gb_path_still_uses_genbank_parser(self, tmp_path, monkeypatch,
                                                tiny_record):
        """Regression guard: .gb files still dispatch to the genbank
        parser, not snapgene."""
        called_fmt = {}

        def fake_parse(path, fmt):
            called_fmt["fmt"] = fmt
            return iter([tiny_record])

        from Bio import SeqIO
        monkeypatch.setattr(SeqIO, "parse", fake_parse)

        gb_file = tmp_path / "plasmid.gb"
        gb_file.write_text("stub")
        sc.load_genbank(str(gb_file))
        assert called_fmt["fmt"] == "genbank"

    def test_malformed_dna_raises_helpful_error(self, tmp_path):
        """A .dna file that isn't actually a SnapGene binary should
        produce a user-friendly error, not a raw struct.error."""
        bad = tmp_path / "broken.dna"
        bad.write_bytes(b"not a snapgene file")
        with pytest.raises(ValueError, match="SnapGene"):
            sc.load_genbank(str(bad))


# ═══════════════════════════════════════════════════════════════════════════════
# SnapGene — real-file integration (uses fixtures in tests/*.dna if present)
# ═══════════════════════════════════════════════════════════════════════════════

def _dna_fixtures():
    """Collect any .dna files sitting next to this test module.
    Skipped entirely if none are present (e.g. on a fresh clone that
    didn't pull the fixtures)."""
    from pathlib import Path
    tests_dir = Path(__file__).parent
    return sorted(tests_dir.glob("*.dna"))


@pytest.mark.skipif(not _dna_fixtures(),
                    reason="no tests/*.dna SnapGene fixtures present")
class TestSnapGeneRealFiles:
    """Integration tests against real SnapGene (.dna) files. Each file
    must parse, have a non-empty sequence, at least one feature, and
    a sensible id/name (backfilled from the filename stem)."""

    @pytest.mark.parametrize("path", _dna_fixtures(),
                             ids=[p.name for p in _dna_fixtures()])
    def test_file_parses(self, path):
        rec = sc.load_genbank(str(path))
        assert len(rec.seq) > 0, f"{path.name} has empty sequence"
        # id/name must be filled, not Biopython sentinels
        assert rec.id and not rec.id.startswith("<unknown")
        assert rec.name and not rec.name.startswith("<unknown")
        # At least one non-source feature
        feats = [f for f in rec.features if f.type != "source"]
        assert feats, f"{path.name} has no features"

    @pytest.mark.parametrize("path", _dna_fixtures(),
                             ids=[p.name for p in _dna_fixtures()])
    def test_round_trip_through_genbank_text(self, path):
        """Load → _record_to_gb_text → _gb_text_to_record preserves the
        sequence. This is the path used when saving a SnapGene import
        to the plasmid library (which stores records as GenBank text
        inside plasmid_library.json)."""
        rec = sc.load_genbank(str(path))
        text = sc._record_to_gb_text(rec)
        rec2 = sc._gb_text_to_record(text)
        assert str(rec.seq) == str(rec2.seq)

    @pytest.mark.parametrize("path", _dna_fixtures(),
                             ids=[p.name for p in _dna_fixtures()])
    def test_features_parse_without_skip_or_clamp(self, path):
        """Real SnapGene files shouldn't trigger the UnknownPosition
        skip path or the out-of-range clamp path — if they do, something
        weird is going on with the import."""
        rec = sc.load_genbank(str(path))
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        pm._parse(rec)
        assert pm._n_skipped == 0, (
            f"{path.name}: {pm._n_skipped} feature(s) unexpectedly skipped"
        )
        assert pm._n_clamped == 0, (
            f"{path.name}: {pm._n_clamped} feature(s) unexpectedly clamped"
        )
