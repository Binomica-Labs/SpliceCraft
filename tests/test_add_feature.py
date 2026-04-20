"""
test_add_feature — AddFeatureModal, PlasmidFeaturePickerModal, and the
Insert-at-cursor pipeline.

Covers:
  - `_parse_qualifier_string` round-trips via `_qualifiers_to_string`
  - `_extract_feature_entries_from_record` respects strand + wrap
  - Modal mount + form gather (save / insert / validation branches)
  - App-side insert shifts feature coords and appends a new SeqFeature
    via the `_rebuild_record_with_edit` pipeline (sacred invariant #9
    remains intact for other features; the new feature lands exactly
    at the cursor with the requested strand + qualifiers)
"""
from __future__ import annotations

import pytest

import splicecraft as sc


TERMINAL_SIZE = (160, 48)


def _build_app(tiny_record, isolated_library) -> sc.PlasmidApp:
    app = sc.PlasmidApp()
    app._preload_record = tiny_record
    return app


# ═══════════════════════════════════════════════════════════════════════════════
# Pure helpers
# ═══════════════════════════════════════════════════════════════════════════════

class TestQualifierParsing:

    def test_single_pair(self):
        assert sc._parse_qualifier_string("gene=lacZ") == {"gene": ["lacZ"]}

    def test_multiple_pairs_semicolon(self):
        got = sc._parse_qualifier_string("gene=lacZ; product=LacZ alpha")
        assert got == {"gene": ["lacZ"], "product": ["LacZ alpha"]}

    def test_whitespace_is_stripped(self):
        got = sc._parse_qualifier_string("  gene  =  lacZ  ;  note  =  test  ")
        assert got == {"gene": ["lacZ"], "note": ["test"]}

    def test_duplicate_keys_collapsed_into_list(self):
        got = sc._parse_qualifier_string("note=a; note=b; note=c")
        assert got == {"note": ["a", "b", "c"]}

    def test_missing_equals_is_ignored(self):
        got = sc._parse_qualifier_string("gene=lacZ; garbage; product=LacZ")
        assert got == {"gene": ["lacZ"], "product": ["LacZ"]}

    def test_empty_input(self):
        assert sc._parse_qualifier_string("") == {}
        assert sc._parse_qualifier_string("   ") == {}

    def test_roundtrip_via_to_string(self):
        original = {"gene": ["lacZ"], "product": ["LacZ alpha"]}
        rendered = sc._qualifiers_to_string(original)
        parsed   = sc._parse_qualifier_string(rendered)
        assert parsed == original


# ═══════════════════════════════════════════════════════════════════════════════
# Feature extraction from a record
# ═══════════════════════════════════════════════════════════════════════════════

class TestExtractFeatureEntries:

    def test_skips_source_feature(self, tiny_record):
        # tiny_record has CDS + misc_feature, no 'source'. Add one to test.
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        tiny_record.features.insert(0, SeqFeature(
            FeatureLocation(0, len(tiny_record.seq), strand=1),
            type="source",
            qualifiers={"organism": ["synthetic"]},
        ))
        entries = sc._extract_feature_entries_from_record(tiny_record)
        # All entries must be non-source
        assert all(e["feature_type"] != "source" for e in entries)
        # CDS + misc_feature → 2 entries
        assert len(entries) == 2

    def test_forward_strand_sequence_matches_slice(self, tiny_record):
        entries = sc._extract_feature_entries_from_record(tiny_record)
        # tiny_record[0] is CDS at [0, 27, +1)
        cds = next(e for e in entries if e["feature_type"] == "CDS")
        assert cds["strand"] == 1
        assert cds["sequence"] == str(tiny_record.seq[0:27]).upper()

    def test_reverse_strand_sequence_is_revcomp(self, tiny_record):
        entries = sc._extract_feature_entries_from_record(tiny_record)
        mf = next(e for e in entries if e["feature_type"] == "misc_feature")
        assert mf["strand"] == -1
        # Fixture: misc_feature at [50, 80, -1). Stored sequence must be the
        # revcomp of the genomic slice (5'→3' of the feature as read).
        genomic = str(tiny_record.seq[50:80]).upper()
        assert mf["sequence"] == sc._rc(genomic)

    def test_qualifiers_preserved(self, tiny_record):
        entries = sc._extract_feature_entries_from_record(tiny_record)
        cds = next(e for e in entries if e["feature_type"] == "CDS")
        assert cds["qualifiers"].get("gene") == ["testA"]


# ═══════════════════════════════════════════════════════════════════════════════
# App-side insert pipeline
# ═══════════════════════════════════════════════════════════════════════════════

class TestInsertFeatureAtCursor:
    """`_insert_feature_at_cursor` must splice in the DNA, shift existing
    feature coords via `_rebuild_record_with_edit`, and append a new
    SeqFeature at the right place."""

    async def test_insert_forward_feature_appends_to_record(self, tiny_record,
                                                             isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 10
            entry = {
                "name": "my-insert",
                "feature_type": "promoter",
                "sequence": "AAATTTGGG",
                "strand": 1,
                "qualifiers": {"note": ["user-added"]},
                "description": "",
            }
            orig_len = len(tiny_record.seq)
            app._insert_feature_at_cursor(entry)
            # Length grew by inserted size
            assert len(app._current_record.seq) == orig_len + 9
            # New feature is the last one
            last = app._current_record.features[-1]
            assert last.type == "promoter"
            assert int(last.location.start) == 10
            assert int(last.location.end) == 19
            assert last.location.strand == 1
            # Qualifiers include our note + auto-label
            assert last.qualifiers.get("note") == ["user-added"]
            assert last.qualifiers.get("label") == ["my-insert"]

    async def test_insert_reverse_splices_rc_into_sequence(self, tiny_record,
                                                            isolated_library):
        """Reverse-strand entries store the 5'→3' of the feature; the bases
        spliced into the genomic strand are the RC."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 30
            entry = {
                "name": "rev-feat",
                "feature_type": "CDS",
                "sequence": "ATGAAATAG",   # feature 5'→3'
                "strand": -1,
                "qualifiers": {},
                "description": "",
            }
            app._insert_feature_at_cursor(entry)
            inserted = str(app._current_record.seq[30:39])
            # Inserted genomic bases == RC of feature sequence
            assert inserted == sc._rc("ATGAAATAG")
            # New feature has strand=-1
            last = app._current_record.features[-1]
            assert last.location.strand == -1

    async def test_insert_shifts_downstream_features(self, tiny_record,
                                                      isolated_library):
        """Existing features after the insertion point must shift by len(insert)."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Capture pre-insert coords of the misc_feature at [50, 80)
            mf_pre = next(f for f in app._current_record.features
                          if f.type == "misc_feature")
            pre_start = int(mf_pre.location.start)
            pre_end   = int(mf_pre.location.end)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 5  # well before the misc_feature
            entry = {
                "name": "x", "feature_type": "misc_feature",
                "sequence": "GGGGGGGGGG",  # 10 bp
                "strand": 1, "qualifiers": {}, "description": "",
            }
            app._insert_feature_at_cursor(entry)
            mf_post = next(f for f in app._current_record.features
                           if f.type == "misc_feature")
            assert int(mf_post.location.start) == pre_start + 10
            assert int(mf_post.location.end)   == pre_end   + 10

    async def test_insert_without_cursor_raises(self, tiny_record,
                                                 isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = -1
            with pytest.raises(RuntimeError):
                app._insert_feature_at_cursor({
                    "name": "x", "feature_type": "CDS",
                    "sequence": "ATG", "strand": 1,
                    "qualifiers": {}, "description": "",
                })

    async def test_insert_marks_dirty(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._unsaved = False
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 0
            app._insert_feature_at_cursor({
                "name": "x", "feature_type": "CDS",
                "sequence": "ATGTAA", "strand": 1,
                "qualifiers": {}, "description": "",
            })
            assert app._unsaved is True


# ═══════════════════════════════════════════════════════════════════════════════
# Modal surface — mount + gather
# ═══════════════════════════════════════════════════════════════════════════════

class TestAddFeatureModal:

    async def test_modal_mounts(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal(have_cursor=False))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            assert isinstance(modal, sc.AddFeatureModal)
            # Every required widget present
            modal.query_one("#addfeat-name")
            modal.query_one("#addfeat-type")
            modal.query_one("#addfeat-seq")
            modal.query_one("#addfeat-quals")

    async def test_gather_rejects_empty_name(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = ""
            modal.query_one("#addfeat-seq").text = "ATG"
            assert modal._gather() is None

    async def test_gather_rejects_invalid_bases(self, tiny_record,
                                                  isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = "bad"
            modal.query_one("#addfeat-seq").text = "ATGXXZZ"
            assert modal._gather() is None

    async def test_gather_accepts_valid_entry(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = "lacZ"
            modal.query_one("#addfeat-seq").text = "atg aaa tag"   # spaces ok
            modal.query_one("#addfeat-quals").value = "gene=lacZ"
            entry = modal._gather()
            assert entry is not None
            assert entry["name"] == "lacZ"
            assert entry["sequence"] == "ATGAAATAG"
            assert entry["qualifiers"] == {"gene": ["lacZ"]}
            assert entry["strand"] == 1   # default = forward

    async def test_gather_iupac_bases_allowed(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = "degenerate"
            modal.query_one("#addfeat-seq").text = "RRYYWWN"
            entry = modal._gather()
            assert entry is not None
            assert entry["sequence"] == "RRYYWWN"


# ═══════════════════════════════════════════════════════════════════════════════
# Save-to-library flow via _add_feature_result
# ═══════════════════════════════════════════════════════════════════════════════

class TestSaveToLibraryFlow:
    """The modal dismisses with {"action": "save", "entry": ...}; the app's
    `_add_feature_result` must persist via _save_features."""

    async def test_save_appends_entry(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sc._features_cache = None
            assert sc._load_features() == []
            app._add_feature_result({
                "action": "save",
                "entry": {
                    "name": "lacZ-alpha",
                    "feature_type": "CDS",
                    "sequence": "ATGAAA",
                    "strand": 1,
                    "qualifiers": {"gene": ["lacZ"]},
                    "description": "",
                },
            })
            assert sc._FEATURES_FILE.exists()
            entries = sc._load_features()
            assert len(entries) == 1
            assert entries[0]["name"] == "lacZ-alpha"

    async def test_save_deduplicates_by_name_and_type(self, tiny_record,
                                                       isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sc._features_cache = None
            app._add_feature_result({"action": "save", "entry": {
                "name": "dup", "feature_type": "CDS",
                "sequence": "A", "strand": 1,
                "qualifiers": {}, "description": "",
            }})
            app._add_feature_result({"action": "save", "entry": {
                "name": "dup", "feature_type": "CDS",
                "sequence": "T", "strand": 1,
                "qualifiers": {}, "description": "",
            }})
            entries = sc._load_features()
            assert len(entries) == 1
            assert entries[0]["sequence"] == "T"   # latest wins


# ═══════════════════════════════════════════════════════════════════════════════
# Feature picker
# ═══════════════════════════════════════════════════════════════════════════════

class TestPlasmidFeaturePickerModal:

    async def test_picker_mounts_with_entries(self, tiny_record,
                                               isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            feats = sc._extract_feature_entries_from_record(tiny_record)
            app.push_screen(sc.PlasmidFeaturePickerModal(feats,
                                                         plasmid_name="tiny"))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            table = modal.query_one("#featpick-table")
            assert table.row_count == len(feats)
