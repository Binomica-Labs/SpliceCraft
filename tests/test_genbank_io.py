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

    def test_spaced_locus_name_does_not_break_save(self):
        """A spaced display name in `record.name` (e.g. leaked from
        `_tui_display_name`) must NOT raise — the GenBank LOCUS line forbids
        whitespace, and an unsanitised name broke EVERY save + autosave
        (user-reported: "Invalid whitespace in '…' for LOCUS line" after a
        scrub Add-to-Map). `_record_to_gb_text` sanitises the LOCUS on a copy."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("ACGT" * 20), id="FFE_6_ENTRY_pCambia2300_GREEN",
            name="FFE 6 ENTRY pCambia2300-GREEN",   # spaces — invalid LOCUS
            annotations={"molecule_type": "DNA", "topology": "circular"})
        gb = sc._record_to_gb_text(rec)              # must NOT raise
        assert gb and gb.splitlines()[0].startswith("LOCUS")
        # No whitespace-bearing display name leaked into the LOCUS token.
        assert "FFE 6 ENTRY" not in gb.splitlines()[0]
        # The caller's record is NOT mutated — the display name survives.
        assert rec.name == "FFE 6 ENTRY pCambia2300-GREEN"
        # Round-trips cleanly.
        assert str(sc._gb_text_to_record(gb).seq) == str(rec.seq)

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


class TestProvenanceStamp:
    """Every serialized GenBank carries a `Created by SpliceCraft v… on …`
    COMMENT (file-level provenance), stamped once and preserved across
    round-trips so the date reflects creation, not the latest re-save."""

    def test_fresh_record_gets_stamp(self, tiny_record):
        gb = sc._record_to_gb_text(tiny_record)
        assert "Created by SpliceCraft v" in gb
        assert sc.__version__ in gb

    def test_roundtrip_does_not_duplicate(self, tiny_record):
        gb1 = sc._record_to_gb_text(tiny_record)
        gb2 = sc._record_to_gb_text(sc._gb_text_to_record(gb1))
        assert gb2.count("Created by SpliceCraft v") == 1

    def test_existing_stamp_preserved(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(
            Seq("ACGTACGTACGT"), id="OLD", name="OLD",
            annotations={"molecule_type": "DNA",
                         "comment": "Created by SpliceCraft v0.9.0 on 2025-01-01"})
        gb = sc._record_to_gb_text(rec)
        assert "v0.9.0 on 2025-01-01" in gb
        assert gb.count("Created by SpliceCraft v") == 1


class TestArrowlessStrandRoundtrip:
    """Arrowless (strand 0) features must survive `_record_to_gb_text` →
    `_gb_text_to_record`. BioPython's GenBank writer emits a plain location
    for BOTH strand 0 and +1, so without the `SpliceCraft_strand=["none"]`
    marker an arrowless feature silently reloads as a forward arrow
    (user-reported 2026-06-05, FRAG-Demo311 "Primer Padding"). Reverse (-1)
    rides complement() natively; double (±2) rides the existing
    `["double"]` marker."""

    @staticmethod
    def _rec(features):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        r = SeqRecord(Seq("ACGT" * 40), id="t", name="t",
                      annotations={"molecule_type": "DNA", "topology": "linear"})
        r.features = list(features)
        return r

    @staticmethod
    def _roundtrip(rec):
        return sc._gb_text_to_record(sc._record_to_gb_text(rec), cache=False)

    def test_arrowless_survives_roundtrip(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._rec([SeqFeature(
            FeatureLocation(0, 10, strand=0), type="misc_feature",
            qualifiers={"label": ["pad"]})])
        f = self._roundtrip(rec).features[0]
        assert f.location.strand == 0                 # NOT coerced to +1
        # The marker is internal — stripped from the live record on read.
        assert "SpliceCraft_strand" not in f.qualifiers

    def test_forward_and_reverse_unchanged(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._rec([
            SeqFeature(FeatureLocation(0, 10, strand=1), type="misc_feature",
                       qualifiers={"label": ["f"]}),
            SeqFeature(FeatureLocation(20, 30, strand=-1), type="misc_feature",
                       qualifiers={"label": ["r"]}),
        ])
        got = {f.qualifiers["label"][0]: f.location.strand
               for f in self._roundtrip(rec).features}
        assert got == {"f": 1, "r": -1}

    def test_encode_does_not_mutate_original(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        feat = SeqFeature(FeatureLocation(0, 10, strand=0), type="misc_feature",
                          qualifiers={"label": ["pad"]})
        sc._record_to_gb_text(self._rec([feat]))
        assert feat.location.strand == 0
        assert "SpliceCraft_strand" not in feat.qualifiers

    def test_double_strand_marker_untouched(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._rec([SeqFeature(
            FeatureLocation(0, 10, strand=1), type="misc_feature",
            qualifiers={"label": ["d"], "SpliceCraft_strand": ["double"]})])
        rec2 = self._roundtrip(rec)
        # The arrowless encoder must not clobber an existing marker, and the
        # double bit must still read back as strand 2 in the map.
        assert rec2.features[0].qualifiers.get("SpliceCraft_strand") == ["double"]
        assert sc.PlasmidMap()._parse(rec2)[0]["strand"] == 2

    def test_wrap_arrowless_survives(self):
        from Bio.SeqFeature import (SeqFeature, FeatureLocation,
                                    CompoundLocation)
        wrap = SeqFeature(
            CompoundLocation([FeatureLocation(150, 160, strand=0),
                              FeatureLocation(0, 5, strand=0)]),
            type="misc_feature", qualifiers={"label": ["w"]})
        assert self._roundtrip(self._rec([wrap])).features[0].location.strand == 0

    def test_plasmidmap_dict_renders_arrowless(self):
        """The seq-panel / map arrowhead is drawn from PlasmidMap._parse's
        dict strand, which must be 0 after the round-trip (not 1)."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._rec([SeqFeature(
            FeatureLocation(0, 10, strand=0), type="misc_feature",
            qualifiers={"label": ["pad"]})])
        d = sc.PlasmidMap()._parse(self._roundtrip(rec))
        assert d[0]["strand"] == 0

    def test_source_feature_is_never_tagged(self):
        """The `source` feature (often strand 0 / None) is skipped — it's not
        drawn as an arrow, and tagging it would noise up every export."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        src = SeqFeature(FeatureLocation(0, 160, strand=0), type="source",
                         qualifiers={"organism": ["synthetic construct"]})
        assert "SpliceCraft_strand" not in sc._record_to_gb_text(self._rec([src]))

    def test_no_arrowless_is_zero_alloc_noop(self):
        """A record with no features, or only ±1-strand features, encodes to
        the SAME list object (no copy) and emits no marker — the common path
        must stay allocation-free and clean."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        assert "SpliceCraft_strand" not in sc._record_to_gb_text(self._rec([]))
        feats = [SeqFeature(FeatureLocation(0, 10, strand=1),
                            type="misc_feature", qualifiers={"label": ["f"]})]
        rec = self._rec(feats)
        assert sc._arrowless_encode_features(rec.features) is rec.features
        assert "SpliceCraft_strand" not in sc._record_to_gb_text(rec)


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
        assert parsed["_schema_version"] == sc._state._CURRENT_SCHEMA_VERSION
        assert parsed["entries"] == entries

    def test_load_survives_corrupted_file(self, isolated_library, caplog):
        """If the library JSON is corrupted, `_load_library` must return []
        and log the exception — never propagate the error to the UI."""
        isolated_library.write_text("{not valid json")
        # Reset in-memory cache so _load_library actually re-reads the file
        sc._state._library_cache = None
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
        from Bio.SeqFeature import SeqFeature
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

    def test_whitespace_only_label_falls_back_to_gene_qualifier(self):
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


class TestFeatLabelFull:
    """`_feat_label_full` is the CANONICAL (untruncated) feature label —
    the value `PlasmidMap._parse` stores and that cloning / export carry
    forward. The 28-char cap lives only in the `_feat_label` display
    wrapper. Regression: an ENA-style operon ("unnamed protein product;
    luxC" … "luxG") collapsed to six identical "unnamed protein product;
    lux" bars on the clone because the gene letter sat at char 29, past
    the cap, and the truncated label was what got persisted."""

    @staticmethod
    def _cds(note):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        return SeqFeature(FeatureLocation(0, 99), type="CDS",
                          qualifiers={"note": [note]})

    def test_full_label_not_truncated(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="gene",
                       qualifiers={"label": ["a" * 50]})
        assert sc._feat_label_full(f) == "a" * 50      # full extractor keeps it
        assert sc._feat_label(f) == "a" * 28           # display wrapper still caps

    def test_placeholder_gene_surfaced_full(self):
        # The reported operon: each gene's distinguishing letter survives.
        for gene in ("luxC", "luxD", "luxA", "luxB", "luxE", "luxG"):
            f = self._cds(f"unnamed protein product; {gene}")
            assert sc._feat_label_full(f) == gene
            assert sc._feat_label(f) == gene          # short → no truncation

    def test_sibling_operon_genes_stay_distinct(self):
        # The crux: six sibling CDS no longer collapse to one label.
        genes = ["luxC", "luxD", "luxA", "luxB", "luxE", "luxG"]
        labels = [sc._feat_label_full(self._cds(f"unnamed protein product; {g}"))
                  for g in genes]
        assert labels == genes
        assert len(set(labels)) == 6                  # all distinct

    def test_bare_placeholder_unchanged(self):
        f = self._cds("unnamed protein product")
        assert sc._feat_label_full(f) == "unnamed protein product"

    def test_placeholder_case_insensitive_and_separators(self):
        assert sc._surface_placeholder_gene(
            "Unnamed Protein Product: luxC") == "luxC"
        assert sc._surface_placeholder_gene(
            "unnamed protein product - luxD") == "luxD"
        assert sc._surface_placeholder_gene("luxC") == "luxC"      # no placeholder
        assert sc._surface_placeholder_gene(
            "ribosome binding site; weak") == "ribosome binding site; weak"

    def test_full_label_collapses_whitespace(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        f = SeqFeature(FeatureLocation(0, 10), type="CDS",
                       qualifiers={"product": ["Line one\nLine two"]})
        assert sc._feat_label_full(f) == "Line one Line two"

    def test_clone_lift_preserves_distinct_gene_labels(self):
        """End-to-end on the modal pipeline the bug travels through:
        `_stamp_insert_feats_on_part` (source feats → part insert_feats)
        then `_attach_insert_feats_to_record` (insert_feats → clone). The
        six lifted operon genes must land as six DISTINCT CDS labels, not
        one collapsed block."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        genes = ["luxC", "luxD", "luxA", "luxB", "luxE", "luxG"]
        # Source feats as `PlasmidMap._parse` would emit them (full labels).
        feats, pos = [], 0
        for g in genes:
            feats.append({"start": pos, "end": pos + 90, "strand": 1,
                          "type": "CDS", "label": g, "color": "white"})
            pos += 100
        insert = "ATG" + "GCT" * 199          # 600 bp, covers all 6 genes
        part = {"name": "pJE-Operon", "sequence": insert}
        n = sc._stamp_insert_feats_on_part(part, insert, feats)
        assert n == 6 and len(part["insert_feats"]) == 6
        rec = SeqRecord(Seq(insert), id="clone", name="clone",
                        annotations={"molecule_type": "DNA"})
        sc._attach_insert_feats_to_record(rec, insert, part["insert_feats"])
        got = [str(f.qualifiers.get("label", [""])[0])
               for f in rec.features if f.type == "CDS"]
        assert sorted(got) == sorted(genes)           # all six survive
        assert len(set(got)) == 6                     # and stay distinct


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
        from Bio.SeqFeature import SeqFeature
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
        pm._parse(rec)
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
        """CommercialSaaS emits adjacent-parts CompoundLocations like
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
# CommercialSaaS (.dna) import — format-detection dispatch (added 2026-04-13)
# ═══════════════════════════════════════════════════════════════════════════════

class TestPlasmidFormatDetection:
    """`_detect_plasmid_format` picks the right Biopython SeqIO format
    key from the file extension so `load_genbank` can load both
    GenBank (.gb/.gbk) and CommercialSaaS (.dna) natively."""

    def test_gb_extension(self):
        assert sc._detect_plasmid_format("x.gb") == "genbank"

    def test_gbk_extension(self):
        assert sc._detect_plasmid_format("/path/to/x.gbk") == "genbank"

    def test_genbank_extension(self):
        assert sc._detect_plasmid_format("x.genbank") == "genbank"

    def test_dna_extension_is_commercialsaas(self):
        assert sc._detect_plasmid_format("x.dna") == sc._BIOPYTHON_DNA_FMT
        assert sc._detect_plasmid_format("x.dna") != "genbank"

    def test_dna_uppercase(self):
        """Case-insensitive — users drag-dropping files from macOS/Windows
        sometimes have capitalized extensions."""
        assert sc._detect_plasmid_format("plasmid.DNA") == sc._BIOPYTHON_DNA_FMT
        assert sc._detect_plasmid_format("foo.GB") == "genbank"

    def test_unknown_extension_defaults_genbank(self):
        """Fallthrough default is GenBank since that's the most common
        plasmid format; parser will raise a clear error if the content
        doesn't match."""
        assert sc._detect_plasmid_format("x.txt") == "genbank"
        assert sc._detect_plasmid_format("no-extension") == "genbank"


class TestCommercialSaaSDispatch:
    """`load_genbank(path)` dispatches to the CommercialSaaS parser when the
    file extension is .dna, so users can import CommercialSaaS-native files
    just by opening them — no need for a CommercialSaaS → GenBank export first."""

    def test_dna_path_calls_commercialsaas_parser(self, tmp_path, monkeypatch,
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
        assert called_fmt["fmt"] == sc._BIOPYTHON_DNA_FMT
        assert called_fmt["path"] == str(dna_file)
        # tiny_record is returned through _pick_single_record
        assert rec is tiny_record

    def test_gb_path_still_uses_genbank_parser(self, tmp_path, monkeypatch,
                                                tiny_record):
        """Regression guard: .gb files still dispatch to the genbank
        parser, not commercialsaas."""
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
        """A .dna file that isn't actually a valid binary plasmid file
        should produce a user-friendly error, not a raw struct.error."""
        bad = tmp_path / "broken.dna"
        bad.write_bytes(b"not a valid binary file")
        with pytest.raises(ValueError, match=r"popular commercial plasmid editor"):
            sc.load_genbank(str(bad))


# ═══════════════════════════════════════════════════════════════════════════════
# CommercialSaaS — real-file integration (uses fixtures in tests/*.dna if present)
# ═══════════════════════════════════════════════════════════════════════════════

def _dna_fixtures():
    """Collect any .dna files sitting next to this test module.
    Skipped entirely if none are present (e.g. on a fresh clone that
    didn't pull the fixtures)."""
    from pathlib import Path
    tests_dir = Path(__file__).parent
    return sorted(tests_dir.glob("*.dna"))


@pytest.mark.skipif(not _dna_fixtures(),
                    reason="no tests/*.dna CommercialSaaS fixtures present")
class TestCommercialSaaSRealFiles:
    """Integration tests against real CommercialSaaS (.dna) files. Each file
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
        sequence. This is the path used when saving a CommercialSaaS import
        to the plasmid library (which stores records as GenBank text
        inside plasmid_library.json)."""
        rec = sc.load_genbank(str(path))
        text = sc._record_to_gb_text(rec)
        rec2 = sc._gb_text_to_record(text)
        assert str(rec.seq) == str(rec2.seq)

    @pytest.mark.parametrize("path", _dna_fixtures(),
                             ids=[p.name for p in _dna_fixtures()])
    def test_features_parse_without_skip_or_clamp(self, path):
        """Real CommercialSaaS files shouldn't trigger the UnknownPosition
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


# ═══════════════════════════════════════════════════════════════════════════════
# GenBank export — `_normalize_for_genbank` + `_export_genbank_to_path`
# ═══════════════════════════════════════════════════════════════════════════════

class TestNormalizeForGenbank:
    """`_normalize_for_genbank` fills in INSDC-mandated fields that Biopython's
    writer needs. Caller's record must never be mutated."""

    def test_fills_molecule_type_default(self, tiny_record):
        # Erase molecule_type to simulate an imported-from-somewhere record
        rec = tiny_record
        rec.annotations = {k: v for k, v in rec.annotations.items()
                           if k != "molecule_type"}
        normalized = sc._normalize_for_genbank(rec)
        assert normalized.annotations["molecule_type"] == "DNA"

    def test_fills_topology_linear_default(self, tiny_record):
        # 2026-05-27 (audit-3 H1): default changed from "circular"
        # to "linear" so a record imported from GFF3 / FASTA without
        # an explicit topology is NOT silently relabelled circular.
        # Topology is biologically load-bearing; circular-default
        # corrupted PCR sim + restriction wrap detection downstream.
        rec = tiny_record
        rec.annotations = {k: v for k, v in rec.annotations.items()
                           if k != "topology"}
        normalized = sc._normalize_for_genbank(rec)
        assert normalized.annotations["topology"] == "linear"

    def test_fills_division_syn_default(self, tiny_record):
        """Synthetic plasmids default to the SYN division code."""
        normalized = sc._normalize_for_genbank(tiny_record)
        assert normalized.annotations["data_file_division"] == "SYN"

    def test_fills_date_in_genbank_format(self, tiny_record):
        """Date must be DD-MMM-YYYY uppercase (e.g. '20-APR-2026')."""
        import re
        normalized = sc._normalize_for_genbank(tiny_record)
        assert re.match(r"^\d{2}-[A-Z]{3}-\d{4}$",
                        normalized.annotations["date"])

    def test_preserves_existing_topology(self, tiny_record):
        """If topology is already set to 'linear', don't silently flip it."""
        tiny_record.annotations["topology"] = "linear"
        normalized = sc._normalize_for_genbank(tiny_record)
        assert normalized.annotations["topology"] == "linear"

    def test_preserves_existing_molecule_type(self, tiny_record):
        tiny_record.annotations["molecule_type"] = "ss-DNA"
        normalized = sc._normalize_for_genbank(tiny_record)
        assert normalized.annotations["molecule_type"] == "ss-DNA"

    def test_caller_record_not_mutated(self, tiny_record):
        """Sacred: normalize is pure — caller's annotations dict is intact."""
        before_keys = set(tiny_record.annotations.keys())
        before_id = id(tiny_record.annotations)
        _ = sc._normalize_for_genbank(tiny_record)
        assert set(tiny_record.annotations.keys()) == before_keys
        assert id(tiny_record.annotations) == before_id

    def test_truncates_long_locus_name(self, tiny_record):
        """NCBI accepts LOCUS names up to 28 chars; longer must be truncated."""
        tiny_record.name = "A" * 50
        normalized = sc._normalize_for_genbank(tiny_record)
        assert len(normalized.name) == 28

    def test_fills_accessions_from_id(self, tiny_record):
        tiny_record.annotations.pop("accessions", None)
        normalized = sc._normalize_for_genbank(tiny_record)
        assert normalized.annotations["accessions"] == [tiny_record.id]

    def test_fills_organism_synthetic_construct(self, tiny_record):
        normalized = sc._normalize_for_genbank(tiny_record)
        assert normalized.annotations["organism"] == "synthetic construct"
        assert "artificial sequences" in normalized.annotations["taxonomy"]


class TestExportGenBankToPath:
    """`_export_genbank_to_path` writes atomically AND round-trip verifies
    before touching the filesystem. A failed round-trip must leave no file."""

    def test_writes_file(self, tiny_record, tmp_path):
        out = tmp_path / "out.gb"
        summary = sc._export_genbank_to_path(tiny_record, out)
        assert out.exists()
        assert summary["path"] == str(out)
        assert summary["bp"] == len(tiny_record.seq)

    def test_roundtrip_sequence_matches(self, tiny_record, tmp_path):
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(tiny_record, out)
        reloaded = sc.load_genbank(str(out))
        assert str(reloaded.seq).upper() == str(tiny_record.seq).upper()

    def test_roundtrip_feature_count_matches(self, tiny_record, tmp_path):
        """Non-source feature count must survive export + reload."""
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(tiny_record, out)
        reloaded = sc.load_genbank(str(out))
        orig_non_source = [f for f in tiny_record.features if f.type != "source"]
        reload_non_source = [f for f in reloaded.features if f.type != "source"]
        assert len(reload_non_source) == len(orig_non_source)

    def test_roundtrip_reverse_strand_preserved(self, tiny_record, tmp_path):
        """Strand=-1 features must round-trip as `complement(...)` locations."""
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(tiny_record, out)
        reloaded = sc.load_genbank(str(out))
        minus = [f for f in reloaded.features if f.location.strand == -1]
        assert len(minus) == 1
        assert minus[0].type == "misc_feature"

    def test_roundtrip_cds_qualifiers_preserved(self, tiny_record, tmp_path):
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(tiny_record, out)
        reloaded = sc.load_genbank(str(out))
        cds = [f for f in reloaded.features if f.type == "CDS"]
        assert len(cds) == 1
        assert cds[0].qualifiers.get("gene") == ["testA"]

    def test_roundtrip_topology_circular_preserved(self, tiny_record, tmp_path):
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(tiny_record, out)
        reloaded = sc.load_genbank(str(out))
        assert reloaded.annotations.get("topology") == "circular"

    def test_wrap_feature_roundtrip(self, tiny_record, tmp_path):
        """Compound (wrap) locations must survive export + reload as join(...)."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        # Build a wrap feature that spans the origin of the 120 bp plasmid
        total = len(tiny_record.seq)
        parts = [FeatureLocation(total - 10, total, strand=1),
                 FeatureLocation(0, 5, strand=1)]
        tiny_record.features.append(SeqFeature(
            CompoundLocation(parts),
            type="misc_feature",
            qualifiers={"label": ["wrap_region"]},
        ))
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(tiny_record, out)
        reloaded = sc.load_genbank(str(out))
        wraps = [f for f in reloaded.features
                 if "wrap_region" in f.qualifiers.get("label", [])]
        assert len(wraps) == 1
        # CompoundLocation has .parts; FeatureLocation does not
        assert hasattr(wraps[0].location, "parts")
        assert len(wraps[0].location.parts) == 2

    def test_atomic_no_leftover_tmp_on_write(self, tiny_record, tmp_path):
        """After a successful export, no `.tmp` hidden files must linger."""
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(tiny_record, out)
        leftovers = list(tmp_path.glob(".*.tmp"))
        assert leftovers == []

    def test_fills_missing_annotations_on_export(self, tmp_path):
        """A bare-bones SeqRecord with no annotations must still export."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        bare = SeqRecord(Seq("ATGCATGCATGC" * 10), id="BARE", name="BARE",
                         description="")
        out = tmp_path / "bare.gb"
        sc._export_genbank_to_path(bare, out)
        reloaded = sc.load_genbank(str(out))
        assert reloaded.annotations["molecule_type"] == "DNA"
        # 2026-05-27 (audit-3 H1): topology defaults to "linear"
        # for records with no explicit value.
        assert reloaded.annotations["topology"] == "linear"

    def test_export_parent_dir_created(self, tiny_record, tmp_path):
        """Intermediate directories are created as needed."""
        out = tmp_path / "nested" / "deep" / "out.gb"
        sc._export_genbank_to_path(tiny_record, out)
        assert out.exists()

    def test_target_path_contains_spaces(self, tiny_record, tmp_path):
        """Paths with spaces (common on macOS / Windows) must not break."""
        out = tmp_path / "my plasmids" / "some name.gb"
        sc._export_genbank_to_path(tiny_record, out)
        assert out.exists()


class TestExportRoundTripDiagnostics:
    """The round-trip guard blocks a lossy export — but it has to say WHAT
    to fix. It used to print raw signature tuples, which in the common case
    (a `strand=None` feature) differ only by a `(?)` buried inside two
    near-identical tuples."""

    @staticmethod
    def _rec(strands):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("ATGC" * 80), id="DIAG", name="DIAG",
                        description="round-trip diagnostics")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        rec.features = [
            SeqFeature(FeatureLocation(10 + 20 * i, 30 + 20 * i, strand=s),
                       type="CDS", qualifiers={"label": [f"feat{i}"]})
            for i, s in enumerate(strands)
        ]
        return rec

    def test_strandless_feature_is_named_with_the_fix(self, tmp_path):
        rec = self._rec([None])
        with pytest.raises(ValueError) as exc:
            sc._export_genbank_to_path(rec, tmp_path / "out.gb")
        msg = str(exc.value)
        assert "has no strand" in msg
        assert "'feat0'" in msg            # names the offending feature
        assert "strand=1" in msg           # and how to fix it
        assert "strand=-1" in msg

    def test_strandless_failure_leaves_no_file(self, tmp_path):
        """The round-trip runs BEFORE the target is touched."""
        out = tmp_path / "out.gb"
        with pytest.raises(ValueError):
            sc._export_genbank_to_path(self._rec([None]), out)
        assert not out.exists()

    def test_one_bad_feature_reports_one_divergence(self, tmp_path):
        """Multiset diff, not a zip of two sorted lists: a single changed
        signature re-sorts, and the old zip then flagged every feature
        after it too — burying the one that actually moved."""
        rec = self._rec([1, 1, None, -1, 1])
        with pytest.raises(ValueError) as exc:
            sc._export_genbank_to_path(rec, tmp_path / "out.gb")
        msg = str(exc.value)
        assert "1 of 5 features diverged" in msg
        assert "'feat2'" in msg            # the strandless one, not a neighbour

    def test_fully_stranded_record_exports_clean(self, tmp_path):
        """The guard must not fire on ordinary records."""
        out = tmp_path / "ok.gb"
        summary = sc._export_genbank_to_path(self._rec([1, -1, 1]), out)
        assert out.exists() and summary["features"] == 3


class TestExportFastaToPath:
    """`_export_fasta_to_path(name, sequence, path)` writes a minimal
    single-record FASTA atomically, validates its inputs, and cleans up
    after itself. Parts Bin + Feature Library both route through it."""

    def test_writes_expected_fasta_text(self, tmp_path):
        out = tmp_path / "out.fa"
        summary = sc._export_fasta_to_path("partA", "ATGCATGC", out)
        assert out.exists()
        assert summary == {"path": str(out), "bp": 8, "name": "partA"}
        assert out.read_text() == ">partA\nATGCATGC\n"

    def test_sequence_is_uppercased(self, tmp_path):
        out = tmp_path / "lower.fa"
        sc._export_fasta_to_path("seq", "atgcatgc", out)
        assert out.read_text().splitlines()[1] == "ATGCATGC"

    def test_roundtrips_through_biopython(self, tmp_path):
        """Anything we write must also be readable via `_parse_fasta_single`."""
        out = tmp_path / "rt.fa"
        sc._export_fasta_to_path("myfeat", "ACGTACGTACGT", out)
        rid, seq = sc._parse_fasta_single(str(out))
        assert rid == "myfeat"
        assert seq == "ACGTACGTACGT"

    def test_empty_name_rejected(self, tmp_path):
        out = tmp_path / "out.fa"
        with pytest.raises(ValueError, match="non-empty record name"):
            sc._export_fasta_to_path("   ", "ATGC", out)
        assert not out.exists()

    def test_empty_sequence_rejected(self, tmp_path):
        out = tmp_path / "out.fa"
        with pytest.raises(ValueError, match="non-empty sequence"):
            sc._export_fasta_to_path("name", "", out)
        assert not out.exists()

    def test_parent_dir_created(self, tmp_path):
        out = tmp_path / "nested" / "deeper" / "out.fa"
        sc._export_fasta_to_path("seq", "ATGC", out)
        assert out.exists()

    def test_atomic_no_leftover_tmp(self, tmp_path):
        """No hidden `.tmp` files must linger after a successful export."""
        out = tmp_path / "clean.fa"
        sc._export_fasta_to_path("seq", "ATGC", out)
        assert list(tmp_path.glob(".*.tmp")) == []

    def test_overwrites_existing_file(self, tmp_path):
        """Atomic replace — a second call clobbers the old content cleanly."""
        out = tmp_path / "out.fa"
        sc._export_fasta_to_path("first", "AAAA", out)
        sc._export_fasta_to_path("second", "CCCC", out)
        assert out.read_text() == ">second\nCCCC\n"

    def test_target_path_contains_spaces(self, tmp_path):
        out = tmp_path / "my parts" / "weird name.fa"
        sc._export_fasta_to_path("seq", "ATGC", out)
        assert out.exists()


class TestRecordToGff3:
    """`_record_to_gff3` produces spec-compliant GFF3 (1.26):
    `##gff-version 3` header, sequence-region pragma, region row
    with `Is_circular=true` for circular plasmids, one tab-separated
    feature row per location part, percent-encoded attribute values.
    Wrap features get two rows sharing the same `ID=...`."""

    def _build(self, seq: str, feats: list[tuple], *,
                circular: bool = True):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        rec = SeqRecord(Seq(seq), id="P", name="P")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular" if circular else "linear"
        n = len(seq)
        for s, e, strand, ftype, label in feats:
            if e < s:
                loc = CompoundLocation([
                    FeatureLocation(s, n, strand=strand),
                    FeatureLocation(0, e, strand=strand),
                ])
            else:
                loc = FeatureLocation(s, e, strand=strand)
            rec.features.append(SeqFeature(
                loc, type=ftype, qualifiers={"label": [label]}
            ))
        return rec

    def test_header_and_sequence_region(self):
        rec = self._build("ACGT" * 10, [])
        gff = sc._record_to_gff3(rec)
        lines = gff.splitlines()
        assert lines[0] == "##gff-version 3"
        assert lines[1] == f"##sequence-region P 1 {40}"

    def test_circular_record_carries_is_circular_flag(self):
        rec = self._build("AAAA", [], circular=True)
        gff = sc._record_to_gff3(rec)
        assert "Is_circular=true" in gff

    def test_linear_record_omits_is_circular_flag(self):
        rec = self._build("AAAA", [], circular=False)
        gff = sc._record_to_gff3(rec)
        assert "Is_circular" not in gff

    def test_feature_row_uses_one_based_inclusive_coords(self):
        # CDS at internal 0-based half-open [3, 9) → GFF 1-based
        # inclusive [4, 9].
        rec = self._build("AA" + "ATGGCCTAA" + "AA",
                            [(2, 11, 1, "CDS", "x")])
        gff = sc._record_to_gff3(rec)
        # CDS row has start=3, end=11 (2+1 .. 11). type column is "CDS".
        cds_lines = [l for l in gff.splitlines()
                     if l and not l.startswith("##") and "CDS" in l]
        assert len(cds_lines) == 1
        cols = cds_lines[0].split("\t")
        assert cols[2] == "CDS"
        assert cols[3] == "3"   # 1-based start = 0-based start + 1
        assert cols[4] == "11"

    def test_strand_column(self):
        rec = self._build("AAAAAAAAAA",
                            [(0, 5, 1, "CDS", "fwd"),
                             (5, 10, -1, "CDS", "rev")])
        gff = sc._record_to_gff3(rec)
        rows = [l.split("\t") for l in gff.splitlines()
                if l and not l.startswith("##") and l.split("\t")[2] == "CDS"]
        strands = [r[6] for r in rows]
        assert "+" in strands and "-" in strands

    def test_wrap_feature_emits_two_rows_same_id(self):
        # 12-bp circular plasmid with a CDS at [9, 3) (wraps origin).
        rec = self._build("AAAAAAAAAAAA",
                            [(9, 3, 1, "CDS", "wrap")])
        gff = sc._record_to_gff3(rec)
        rows = [l for l in gff.splitlines()
                if l and not l.startswith("##") and "CDS" in l]
        assert len(rows) == 2
        # Both rows share the same ID attribute.
        ids = []
        for r in rows:
            attrs = r.split("\t")[-1]
            for kv in attrs.split(";"):
                if kv.startswith("ID="):
                    ids.append(kv)
        assert len(set(ids)) == 1, f"split-feature rows must share ID: {ids}"

    def test_skip_source_feature(self):
        """The synthesised `region` row covers the whole record, so a
        GenBank source feature would double-list the span — filter it."""
        rec = self._build("AAAA", [(0, 4, 1, "source", "src")])
        gff = sc._record_to_gff3(rec)
        # Only one row of type "region"; no "source" rows.
        rows = [l for l in gff.splitlines()
                if l and not l.startswith("##")]
        types = {r.split("\t")[2] for r in rows}
        assert "source" not in types
        assert "region" in types

    def test_attribute_percent_encoding(self):
        """A label with `;` or `=` would break the GFF3 attributes
        column without percent-encoding."""
        rec = self._build("AAAAAAAA",
                            [(0, 8, 1, "misc_feature", "weird=label;here")])
        gff = sc._record_to_gff3(rec)
        # Raw `;` and `=` from the label must not appear unescaped in
        # the Name= value.
        rows = [l for l in gff.splitlines()
                if l and not l.startswith("##") and "misc_feature" in l]
        assert len(rows) == 1
        attrs = rows[0].split("\t")[-1]
        # Find the Name= attribute and verify the literal `;` / `=`
        # in the original label are encoded.
        for kv in attrs.split(";"):
            if kv.startswith("Name="):
                assert "weird%3Dlabel%3Bhere" in kv


class TestExportGffToPath:
    def test_writes_file(self, tmp_path):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("ACGT" * 10), id="t", name="t")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(0, 12, strand=1),
            type="CDS", qualifiers={"label": ["x"]}
        ))
        out = tmp_path / "x.gff3"
        summary = sc._export_gff_to_path(rec, out)
        assert out.exists()
        assert summary["bp"] == 40
        assert summary["features"] == 1
        text = out.read_text()
        assert text.startswith("##gff-version 3")


class TestActionExportFasta:
    """Whole-plasmid FASTA export wired in 0.6.0.0: File → Export as
    FASTA pushes `FastaExportModal` pre-populated with the loaded
    record's sequence, mirroring `action_export_genbank` for the .gb
    side. Bare app (no record) must surface a friendly notify rather
    than push a modal with empty content."""

    @pytest.mark.asyncio
    async def test_action_export_fasta_pushes_modal(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.action_export_fasta()
            await pilot.pause()
            top = app.screen
            assert isinstance(top, sc.FastaExportModal)
            assert top._name == tiny_record.name
            assert top._sequence == str(tiny_record.seq)
            # 2026-05-06: the typed-path Input was replaced with a
            # dir picker + filename Input. The default filename still
            # carries the .fa extension; the parent dir is split out
            # into `_start_dir`.
            assert top._default_filename.endswith(".fa")

    @pytest.mark.asyncio
    async def test_action_export_fasta_no_record_notifies(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.action_export_fasta()
            await pilot.pause()
            # No modal pushed — still on the bare app screen.
            assert not isinstance(app.screen, sc.FastaExportModal)


class TestFastaImportViaOpenModal:
    """OpenFileModal (Ctrl+O) supports FASTA in addition to .gb / .gbk /
    .dna / .genbank — picking a `.fa` file routes to the FASTA loader
    instead of `load_genbank` and yields a `topology="linear"` SeqRecord.
    Locks in the 2026-05-06 unification of plasmid + FASTA ingestion."""

    def test_fasta_path_to_record_builds_linear_record(self, tmp_path):
        fa = tmp_path / "demo.fa"
        fa.write_text(">my_seq\nATCGATCGATCG\n")
        rec = sc._fasta_path_to_record(str(fa))
        assert str(rec.seq) == "ATCGATCGATCG"
        assert rec.id == "my_seq"
        # FASTA carries no topology annotation — the loader defaults to
        # linear so a chromosome-chunk import doesn't display as a fake
        # circle. Users flip to circular via the map view-mode toggle.
        assert rec.annotations["topology"] == "linear"
        assert rec.annotations["molecule_type"] == "DNA"
        assert "Imported from FASTA" in rec.description

    def test_fasta_path_to_record_rejects_multi_record(self, tmp_path):
        fa = tmp_path / "multi.fa"
        fa.write_text(">a\nACGT\n>b\nGCAT\n")
        with pytest.raises(ValueError, match="Multi-sequence"):
            sc._fasta_path_to_record(str(fa))

    def test_open_modal_highlight_map_includes_fasta_and_dna(self):
        # Map covers plasmid + .dna + every FASTA extension; pink for
        # FASTA, orange for .dna, green for plasmid text formats.
        m = sc._PLASMID_PICKER_HIGHLIGHT_MAP
        assert m[".fa"] == sc._PICKER_FASTA_STYLE
        assert m[".fasta"] == sc._PICKER_FASTA_STYLE
        assert m[".dna"] == sc._PICKER_DNA_STYLE
        assert m[".gb"] == sc._PICKER_PLASMID_STYLE
        assert m[".gbk"] == sc._PICKER_PLASMID_STYLE
        # Hot pink + orange must remain distinct so the user can tell
        # the two binary/sequence formats apart at a glance.
        assert sc._PICKER_FASTA_STYLE != sc._PICKER_DNA_STYLE


class TestMultiRecordFastaImport:
    """Multi-record FASTA → new collection workflow.

    Picking a multi-record `.fa` in OpenFileModal pushes
    `MultiRecordFastaModal`. On confirm the helpers build SeqRecords
    (with topology auto-detected per-record) and create a collection
    via `_create_fasta_collection`.
    """

    def test_detect_fasta_topology_circular_hint(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec = SeqRecord(Seq("ACGT"), id="x",
                         description="x circular plasmid backbone")
        assert sc._detect_fasta_topology(rec) == "circular"

    def test_detect_fasta_topology_default_linear(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        # No "circular" / "plasmid" in description → linear default.
        rec = SeqRecord(Seq("ACGT"), id="x", description="x chromosome chunk")
        assert sc._detect_fasta_topology(rec) == "linear"
        # Empty description: also linear.
        rec2 = SeqRecord(Seq("ACGT"), id="x", description="")
        assert sc._detect_fasta_topology(rec2) == "linear"

    def test_fasta_records_to_seqrecords_dedupes_ids(self, tmp_path):
        fa = tmp_path / "dupes.fa"
        fa.write_text(
            ">seq\nACGT\n>seq\nGCAT\n>seq\nTTAA\n"
        )
        from Bio import SeqIO
        parsed = list(SeqIO.parse(str(fa), "fasta"))
        out = sc._fasta_records_to_seqrecords(parsed, str(fa))
        assert len(out) == 3
        ids = [r.id for r in out]
        # All three input records had id "seq"; the helper appends
        # `_2`, `_3` to keep them unique so library `add_entry`
        # doesn't dedupe-by-id and silently drop the duplicates.
        assert ids == ["seq", "seq_2", "seq_3"]

    def test_fasta_records_to_seqrecords_skips_non_iupac(self, tmp_path):
        fa = tmp_path / "mixed.fa"
        # Second record contains '@' which fails the IUPAC check.
        fa.write_text(">good\nACGT\n>bad\nAC@T\n>good2\nGGCC\n")
        from Bio import SeqIO
        parsed = list(SeqIO.parse(str(fa), "fasta"))
        out = sc._fasta_records_to_seqrecords(parsed, str(fa))
        # Bad record dropped; the other two come through.
        assert [r.id for r in out] == ["good", "good2"]

    def test_fasta_records_to_seqrecords_per_record_topology(self, tmp_path):
        fa = tmp_path / "mixed_topo.fa"
        fa.write_text(
            ">linear_one chromosome\nACGT\n"
            ">circular_two circular plasmid\nGCAT\n"
        )
        from Bio import SeqIO
        parsed = list(SeqIO.parse(str(fa), "fasta"))
        out = sc._fasta_records_to_seqrecords(parsed, str(fa))
        assert len(out) == 2
        assert out[0].annotations["topology"] == "linear"
        assert out[1].annotations["topology"] == "circular"

    def test_create_fasta_collection_round_trip(self, tmp_path,
                                                  isolated_library):
        fa = tmp_path / "set.fa"
        fa.write_text(
            ">a\nACGTACGT\n"
            ">b circular\nGCATGCAT\n"
            ">c\nTTAATTAA\n"
        )
        from Bio import SeqIO
        parsed = list(SeqIO.parse(str(fa), "fasta"))
        records = sc._fasta_records_to_seqrecords(parsed, str(fa))
        sc._create_fasta_collection("FASTA Set", records, str(fa))
        # Active collection switched to the new one.
        assert sc._get_active_collection_name() == "FASTA Set"
        # Every record present in the active collection's plasmid list.
        coll = sc._find_collection("FASTA Set")
        assert coll is not None
        assert len(coll["plasmids"]) == 3
        names = {p.get("id") for p in coll["plasmids"]}
        assert names == {"a", "b", "c"}
        # Library file mirrors the collection's plasmids.
        assert len(sc._load_library()) == 3

    def test_create_fasta_collection_rejects_duplicate_name(
            self, tmp_path, isolated_library):
        # Pre-create a collection so the next call collides.
        sc._save_collections([{
            "name": "Existing", "description": "",
            "plasmids": [], "saved": "2026-05-06",
        }])
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        records = [SeqRecord(Seq("ACGT"), id="x", name="x",
                              description="")]
        with pytest.raises(ValueError, match="already exists"):
            sc._create_fasta_collection(
                "Existing", records, str(tmp_path / "x.fa"),
            )

    @pytest.mark.asyncio
    async def test_multi_record_modal_validates_name(self,
                                                       isolated_library):
        # Pre-create a collision so we can also exercise the dupe path.
        sc._save_collections([{
            "name": "Already there", "description": "",
            "plasmids": [], "saved": "2026-05-06",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.MultiRecordFastaModal(3, default_name="")
            await app.push_screen(modal)
            await pilot.pause()
            inp = modal.query_one("#multi-fasta-input", sc.Input)
            status = modal.query_one("#multi-fasta-status", sc.Static)
            # Empty name → status error, modal stays.
            inp.value = ""
            modal._submit()
            assert app.screen_stack[-1] is modal
            assert "empty" in str(status.render()).lower()
            # Collision → status error, modal stays.
            inp.value = "Already there"
            modal._submit()
            assert app.screen_stack[-1] is modal
            assert "already exists" in str(status.render()).lower()
            # Valid → dismisses with the name.
            inp.value = "New collection"
            modal._submit()
            await pilot.pause()
            assert app.screen_stack[-1] is not modal


class TestFastaTopologyHeuristic:
    """FASTA carries no topology field, so import guesses from the
    description. Guessing CIRCULAR is the expensive direction: it fabricates
    a junction between two ends that are not joined, and the origin-wrap
    restriction scan (sacred invariant #6) then reports sites across it.

    Every "linear" case below is a real FASTA header shape that came back
    circular before 2026-08-18, because `"plasmid"`/`"circular"` were matched
    as bare substrings with nothing weighing the words next to them."""

    @staticmethod
    def _rec(desc):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        return SeqRecord(Seq("ACGT"), id="x", description=desc)

    @pytest.mark.parametrize("desc", [
        "pUC19 complete circular plasmid sequence",
        "pACYC184 plasmid",
        "circular chromosome, E. coli K-12",
        "x circular plasmid backbone",
        "circularized ligation product",
    ])
    def test_positive_signal_reads_circular(self, desc):
        assert sc._detect_fasta_topology(self._rec(desc)) == "circular"

    @pytest.mark.parametrize("desc", [
        "pUC19 plasmid, linearized",              # says both — linear wins
        "linear plasmid pBSSB1 Borrelia burgdorferi",
        "plasmid backbone fragment (PCR product)",
        "PCR amplicon from plasmid template",
        "contig_17 Streptomyces linear plasmid SCP1",
        "gBlock, non-circular",
        "synthetic plasmid, not circular",
        "cpGFP circularly permuted green fluorescent protein",
    ])
    def test_explicit_linearity_beats_the_hint(self, desc):
        assert sc._detect_fasta_topology(self._rec(desc)) == "linear", (
            f"{desc!r} was read as circular — SpliceCraft would close a "
            f"circle the description says is open"
        )

    @pytest.mark.parametrize("desc", [
        "", "chromosome 3 scaffold", "synthetic construct",
        "collinear repeat array",                 # "collinear" is not "linear"
    ])
    def test_no_signal_defaults_linear(self, desc):
        assert sc._detect_fasta_topology(self._rec(desc)) == "linear"


class TestTopologyFromLocusLine:
    """`_topology_from_gb_text` reads the LOCUS line's topology FIELD.

    Three call sites (library Kind column, agent `list-library`, the BABS
    index) used to ask `"linear" in gb_text[:200].lower()` instead, which
    also reads the entry's own NAME and its DEFINITION."""

    @staticmethod
    def _gb(locus_name, topo, definition="."):
        return (f"LOCUS       {locus_name:<24}2686 bp    DNA     {topo} "
                f"SYN 27-OCT-2024\nDEFINITION  {definition}\n"
                f"ACCESSION   x\n//\n")

    def test_reads_the_topology_field(self):
        assert sc._topology_from_gb_text(self._gb("pUC19", "circular")) \
            == "circular"
        assert sc._topology_from_gb_text(self._gb("pUC19cut", "linear")) \
            == "linear"

    def test_the_entrys_own_name_is_not_a_topology_claim(self):
        """An entry called `pLinear2` is not linear because of its name."""
        assert sc._topology_from_gb_text(
            self._gb("pLinear2", "circular")) == "circular"

    def test_the_definition_is_not_a_topology_claim(self):
        assert sc._topology_from_gb_text(
            self._gb("pACYC184", "circular",
                     "linearized derivative described in ref 3")) == "circular"
        assert sc._topology_from_gb_text(
            self._gb("frag1", "linear",
                     "circular vector backbone, cut with EcoRI")) == "linear"

    @pytest.mark.parametrize("text", [
        "", None, 42, "not a genbank file at all",
        "DEFINITION  linear thing\nLOCUS x 10 bp DNA linear\n",   # no LOCUS 1st
    ])
    def test_falls_back_to_the_default(self, text):
        assert sc._topology_from_gb_text(text) == "circular"
        assert sc._topology_from_gb_text(text, default="") == ""

    def test_locus_line_without_a_topology_token(self):
        assert sc._topology_from_gb_text(
            "LOCUS       x   100 bp ds-DNA SYN 01-JAN-2026\n//\n") == "circular"
        assert sc._topology_from_gb_text(
            "LOCUS       x   100 bp ds-DNA SYN 01-JAN-2026\n//\n",
            default="") == ""


class TestLoadFileHyphenDisplayName:
    """Loading a hyphenated file whose LOCUS got underscored (`FRAG-Ds3.gb`
    with `LOCUS FRAG_Ds3`, the petunia running-log) must surface the
    HYPHENATED filename as the display name, not the underscored LOCUS.

    The mechanism: the LOCUS-safe `safe_stem` folds spaces AND hyphens to
    `_` (matching the GenBank LOCUS charset), so `stem != safe_stem` fires
    the `_tui_display_name` override for hyphenated names. [INV-98]"""

    def _write_gb(self, path, locus, seq="ATGC" * 30):
        # Hand-write a minimal GenBank with an explicitly UNDERSCORED LOCUS,
        # reproducing how the petunia assembler emitted PCR-/FRAG- files.
        body = (
            f"LOCUS       {locus:<16} {len(seq)} bp    DNA     linear   UNK\n"
            "FEATURES             Location/Qualifiers\n"
            "ORIGIN\n"
            f"        1 {seq.lower()}\n"
            "//\n"
        )
        path.write_text(body, encoding="utf-8")

    def test_hyphenated_filename_keeps_hyphen(self, tmp_path):
        p = tmp_path / "FRAG-Ds3.gb"
        self._write_gb(p, "FRAG_Ds3")
        rec = sc.load_genbank(str(p))
        assert rec.name == "FRAG_Ds3"                     # LOCUS as written
        assert getattr(rec, "_tui_display_name", None) == "FRAG-Ds3"
        assert sc.PlasmidApp._record_display_name(rec) == "FRAG-Ds3"

    def test_pcr_prefix_name(self, tmp_path):
        p = tmp_path / "PCR-DODA.gb"
        self._write_gb(p, "PCR_DODA")
        rec = sc.load_genbank(str(p))
        assert sc.PlasmidApp._record_display_name(rec) == "PCR-DODA"

    def test_plain_name_no_spurious_override(self, tmp_path):
        # A filename that matches its LOCUS (no spaces/hyphens/truncation)
        # must NOT get a display-name override — the LOCUS shows as-is.
        p = tmp_path / "pUC19.gb"
        self._write_gb(p, "pUC19")
        rec = sc.load_genbank(str(p))
        assert getattr(rec, "_tui_display_name", None) is None
        assert sc.PlasmidApp._record_display_name(rec) == "pUC19"

    def test_meaningful_locus_generic_filename_unchanged(self, tmp_path):
        # Regression guard: a generic filename with a DIFFERENT meaningful
        # LOCUS (not a mangled form of the stem) keeps showing the LOCUS —
        # the fix must not hijack every load to the filename.
        p = tmp_path / "sequence.gb"
        self._write_gb(p, "NC_001416")
        rec = sc.load_genbank(str(p))
        assert getattr(rec, "_tui_display_name", None) is None
        assert sc.PlasmidApp._record_display_name(rec) == "NC_001416"


class TestDisplayNameCommentRoundTrip:
    """SC-E: the human display name survives a gb_text / FILE round-trip via a
    COMMENT marker, so an exported construct re-imports under its real name —
    no `rename-plasmid` after every save."""

    def _rec(self, name="FIRE 13 CAM RUBYFIRE-Binary", locus="FIRE_13_CAM_RUBY"):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        r = SeqRecord(Seq("ATGC" * 30), id=locus, name=locus,
                      annotations={"molecule_type": "DNA", "topology": "circular"})
        r._tui_display_name = name
        return r

    def test_gb_text_roundtrip_preserves_display_name(self):
        gb = sc._record_to_gb_text(self._rec())
        assert "SpliceCraft-name: FIRE 13 CAM RUBYFIRE-Binary" in gb
        r2 = sc._gb_text_to_record(gb)
        assert (sc.PlasmidApp._record_display_name(r2)
                == "FIRE 13 CAM RUBYFIRE-Binary")

    def test_stamp_is_idempotent(self):
        r = self._rec()
        once = sc._record_to_gb_text(r)
        twice = sc._record_to_gb_text(sc._gb_text_to_record(once))
        assert twice.count("SpliceCraft-name:") == 1

    def test_file_roundtrip_beats_mangled_filename(self, tmp_path):
        gb = sc._record_to_gb_text(self._rec())
        p = tmp_path / "FIRE_13_CAM_RUBY.gb"        # mangled-LOCUS filename
        p.write_text(gb, encoding="utf-8")
        rec = sc.load_genbank(str(p))
        # The COMMENT name wins over the underscored filename stem.
        assert (sc.PlasmidApp._record_display_name(rec)
                == "FIRE 13 CAM RUBYFIRE-Binary")

    def test_no_stamp_without_display_name(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        r = SeqRecord(Seq("ATGC" * 10), id="pUC19", name="pUC19",
                      annotations={"molecule_type": "DNA"})
        gb = sc._record_to_gb_text(r)
        assert "SpliceCraft-name:" not in gb
        assert getattr(sc._gb_text_to_record(gb), "_tui_display_name", None) is None

    def test_no_stamp_when_display_equals_locus(self):
        assert "SpliceCraft-name:" not in sc._record_to_gb_text(
            self._rec(name="pUC19", locus="pUC19"))
