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

class TestAnnotateWithFeature:
    """`_annotate_with_feature` adds a SeqFeature to the loaded record
    spanning the given range without modifying the underlying DNA. This
    is the shared backend for both the AddFeatureModal "Insert feature"
    button and the agent-API `add-feature` endpoint — single source of
    truth for "annotate existing bases".

    Pre-2026-04-30 the modal button instead spliced new DNA at the
    cursor (`_insert_feature_at_cursor`); that path was removed in
    favour of "select region → Ctrl+F → annotate range" which lets the
    user mark up an existing region without changing its length. New
    DNA insertion lives in Ctrl+E (EditSeqDialog) for users who need it.
    """

    async def test_forward_feature_appends_with_correct_coords(
        self, tiny_record, isolated_library,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            orig_len = len(tiny_record.seq)
            entry = {
                "name": "my-feat",
                "feature_type": "promoter",
                "strand": 1,
                "qualifiers": {"note": ["user-added"]},
            }
            app._annotate_with_feature(10, 25, entry)
            # Sequence is unchanged — the whole point of "annotate".
            assert len(app._current_record.seq) == orig_len
            # New feature is the last one.
            last = app._current_record.features[-1]
            assert last.type == "promoter"
            assert int(last.location.start) == 10
            assert int(last.location.end) == 25
            assert last.location.strand == 1
            # Qualifiers include the user's note + the auto-label.
            assert last.qualifiers.get("note") == ["user-added"]
            assert last.qualifiers.get("label") == ["my-feat"]

    async def test_reverse_strand_records_strand_minus_one(
        self, tiny_record, isolated_library,
    ):
        """Reverse-strand annotations don't touch the DNA — they only
        flag the SeqFeature's strand. (The displayed bases are still
        the underlying top-strand bases.)"""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            orig_seq = str(tiny_record.seq)
            entry = {
                "name": "rev-feat", "feature_type": "CDS",
                "strand": -1, "qualifiers": {},
            }
            app._annotate_with_feature(30, 39, entry)
            # Sequence unchanged.
            assert str(app._current_record.seq) == orig_seq
            last = app._current_record.features[-1]
            assert last.location.strand == -1
            assert int(last.location.start) == 30
            assert int(last.location.end) == 39

    async def test_other_features_keep_their_coords(
        self, tiny_record, isolated_library,
    ):
        """Annotating doesn't shift any existing features — opposite
        of the old splice behaviour."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            mf_pre = next(f for f in app._current_record.features
                          if f.type == "misc_feature")
            pre_start = int(mf_pre.location.start)
            pre_end   = int(mf_pre.location.end)
            entry = {
                "name": "x", "feature_type": "misc_feature",
                "strand": 1, "qualifiers": {},
            }
            # Annotate well before the existing misc_feature.
            app._annotate_with_feature(0, 5, entry)
            mf_post = next(f for f in app._current_record.features
                           if f.type == "misc_feature"
                           and f.qualifiers.get("label") != ["x"])
            assert int(mf_post.location.start) == pre_start
            assert int(mf_post.location.end)   == pre_end

    async def test_no_record_raises(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._current_record = None
            with pytest.raises(RuntimeError):
                app._annotate_with_feature(0, 5, {
                    "name": "x", "feature_type": "CDS",
                    "strand": 1, "qualifiers": {},
                })

    async def test_zero_length_range_raises(self, tiny_record,
                                              isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            with pytest.raises(ValueError):
                app._annotate_with_feature(5, 5, {
                    "name": "x", "feature_type": "CDS",
                    "strand": 1, "qualifiers": {},
                })

    async def test_out_of_range_raises(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            n = len(app._current_record.seq)
            with pytest.raises(ValueError):
                app._annotate_with_feature(n + 5, n + 10, {
                    "name": "x", "feature_type": "CDS",
                    "strand": 1, "qualifiers": {},
                })

    async def test_annotate_marks_dirty(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._unsaved = False
            app._annotate_with_feature(0, 6, {
                "name": "x", "feature_type": "CDS",
                "strand": 1, "qualifiers": {},
            })
            assert app._unsaved is True

    async def test_wrap_range_builds_compound_location(
        self, isolated_library,
    ):
        """end < start should produce a CompoundLocation with two
        FeatureLocation parts (tail [start, n) + head [0, end)) — the
        same wrap-aware shape every other code path expects (sacred
        invariant #9)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import CompoundLocation
        rec = SeqRecord(Seq("A" * 100), id="wrap_anno", name="wrap_anno",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(95, 5, {
                "name": "wrap-feat", "feature_type": "misc_feature",
                "strand": 1, "qualifiers": {},
            })
            last = app._current_record.features[-1]
            assert isinstance(last.location, CompoundLocation)
            parts = list(last.location.parts)
            assert int(parts[0].start) == 95
            assert int(parts[0].end)   == 100
            assert int(parts[1].start) == 0
            assert int(parts[1].end)   == 5

    async def test_strand_zero_accepted(self, tiny_record, isolated_library):
        """Arrowless / unknown-strand annotations (strand=0) should
        round-trip through `_annotate_with_feature` without crashing —
        BioPython needs `strand=None` for that case."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(0, 10, {
                "name": "ds", "feature_type": "misc_feature",
                "strand": 0, "qualifiers": {},
            })
            last = app._current_record.features[-1]
            # FeatureLocation(strand=None) → location.strand == None
            assert last.location.strand is None
            assert int(last.location.start) == 0
            assert int(last.location.end)   == 10

    async def test_color_in_entry_writes_apeinfo_qualifiers(
        self, tiny_record, isolated_library,
    ):
        """Regression for 2026-05-26 "choosing color does not seem to
        apply that color to the feature in the new feature modal"
        report: when the AddFeatureModal's color picker sets a
        custom color, the entry dict carries `color="#…"` to
        `_annotate_with_feature`. Pre-fix the implementation
        silently ignored `entry["color"]` — the SeqFeature went
        onto the record with no color qualifier, so the bar in the
        seq panel / map stayed on the type-default palette color.
        Fix writes both `ApEinfo_fwdcolor` and `ApEinfo_revcolor`
        (mirroring `_apply_feature_edit`), so create + edit are
        symmetric AND the choice round-trips through `.gb` export."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(0, 9, {
                "name": "colored", "feature_type": "misc_feature",
                "strand": 1, "qualifiers": {},
                "color": "#ff8800",
            })
            last = app._current_record.features[-1]
            assert last.qualifiers.get("ApEinfo_fwdcolor") == ["#ff8800"]
            assert last.qualifiers.get("ApEinfo_revcolor") == ["#ff8800"]

    async def test_no_color_in_entry_omits_apeinfo_qualifiers(
        self, tiny_record, isolated_library,
    ):
        """Auto-color (no `color` key, or `color=None`, or empty
        string) MUST NOT write the qualifiers — that would lock in
        whatever the type-default happened to be at annotate time
        instead of letting the renderer's palette logic pick a
        fresh one. Mirrors the `_apply_feature_edit` else branch."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            for color_val in (None, "", "  "):
                app._annotate_with_feature(0, 6, {
                    "name": f"auto-{color_val!r}",
                    "feature_type": "misc_feature",
                    "strand": 1, "qualifiers": {},
                    "color": color_val,
                })
                last = app._current_record.features[-1]
                assert "ApEinfo_fwdcolor" not in last.qualifiers
                assert "ApEinfo_revcolor" not in last.qualifiers

    async def test_non_string_color_in_entry_is_rejected(
        self, tiny_record, isolated_library,
    ):
        """Defensive: a programmatic / agent-API caller could hand
        a non-string color value (int from a misread JSON, list
        from a malformed agent payload). The isinstance guard must
        drop those rather than coercing them via `str()` into a
        qualifier value that downstream colour parsers can't make
        sense of (e.g. `str([])` → `'[]'`)."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            for color_val in (123, ["#ff0000"], {"hex": "#ff0000"}, True):
                app._annotate_with_feature(0, 6, {
                    "name": f"bad-{type(color_val).__name__}",
                    "feature_type": "misc_feature",
                    "strand": 1, "qualifiers": {},
                    "color": color_val,
                })
                last = app._current_record.features[-1]
                assert "ApEinfo_fwdcolor" not in last.qualifiers
                assert "ApEinfo_revcolor" not in last.qualifiers

    async def test_color_with_whitespace_is_stripped(
        self, tiny_record, isolated_library,
    ):
        """A color value with leading / trailing whitespace lands
        as the stripped form. Defensive (the modal doesn't
        introduce whitespace today, but a programmatic caller or a
        future Input widget could)."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(0, 6, {
                "name": "padded", "feature_type": "misc_feature",
                "strand": 1, "qualifiers": {},
                "color": "  #abc123  \n",
            })
            last = app._current_record.features[-1]
            assert last.qualifiers["ApEinfo_fwdcolor"] == ["#abc123"]
            assert last.qualifiers["ApEinfo_revcolor"] == ["#abc123"]

    async def test_double_strand_round_trips_via_qualifier(
        self, tiny_record, isolated_library,
    ):
        """Regression for 2026-05-26 "double arrow option" report:
        the AddFeatureModal has a `Double (◀▶)` radio (strand=2)
        but BioPython's FeatureLocation only encodes ±1 / 0 /
        None. Without a custom qualifier the save→reload cycle
        silently collapsed strand=2 → strand=0, so the user's
        double-arrow choice was lost the moment `_parse` re-ran.
        `_annotate_with_feature_impl` now writes
        `SpliceCraft_strand=["double"]` for strand=2;
        `PlasmidMap._parse` reads it back so the dict-side
        strand is 2 and the double-arrow rendering survives.
        Same for `_apply_feature_edit` (edit-path)."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(0, 12, {
                "name": "ds", "feature_type": "misc_feature",
                "strand": 2, "qualifiers": {},
            })
            last = app._current_record.features[-1]
            # BioPython strand is None (since 2 isn't representable),
            # but the SpliceCraft qualifier records the intent.
            assert last.location.strand is None
            assert last.qualifiers.get("SpliceCraft_strand") == [
                "double",
            ]
            # `_parse` reads the qualifier back as strand=2.
            feats = sc.PlasmidMap._parse(
                sc.PlasmidMap.__new__(sc.PlasmidMap),
                app._current_record,
            )
            ds = next(f for f in feats if f["label"] == "ds")
            assert ds["strand"] == 2

    async def test_double_strand_clears_qualifier_when_changed(
        self, tiny_record, isolated_library,
    ):
        """If the user edits a strand=2 feature down to strand=1
        (forward), the qualifier MUST be removed so the next
        reload doesn't re-promote it back to strand=2. Same for
        the reverse direction (1→2): the qualifier appears on
        write."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Start with a strand=2 feature.
            app._annotate_with_feature(0, 12, {
                "name": "ds", "feature_type": "misc_feature",
                "strand": 2, "qualifiers": {},
            })
            # Edit it down to strand=1 via `_apply_feature_edit`.
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            idx = next(i for i, f in enumerate(pm._feats)
                        if f.get("label") == "ds")
            app._apply_feature_edit({
                "idx": idx, "strand": 1,
            })
            last = next(f for f in app._current_record.features
                         if f.qualifiers.get("label") == ["ds"])
            assert "SpliceCraft_strand" not in last.qualifiers
            assert last.location.strand == 1
            # Now edit back up to strand=2 — qualifier reappears.
            app._apply_feature_edit({
                "idx": idx, "strand": 2,
            })
            last = next(f for f in app._current_record.features
                         if f.qualifiers.get("label") == ["ds"])
            assert last.qualifiers.get("SpliceCraft_strand") == [
                "double",
            ]
            assert last.location.strand is None

    async def test_color_does_not_clobber_user_qualifiers(
        self, tiny_record, isolated_library,
    ):
        """Writing the color qualifiers must not erase / overwrite
        unrelated qualifiers the user supplied via the Qualifiers
        Input (`key=value; key=value`). Regression guard for a
        future refactor that swaps the dict-update for a full
        replacement."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(0, 9, {
                "name": "tagged", "feature_type": "misc_feature",
                "strand": 1,
                "qualifiers": {"note": ["from-user"],
                                "db_xref": ["UniProt:Q12345"]},
                "color": "#abcdef",
            })
            last = app._current_record.features[-1]
            assert last.qualifiers["note"] == ["from-user"]
            assert last.qualifiers["db_xref"] == ["UniProt:Q12345"]
            assert last.qualifiers["ApEinfo_fwdcolor"] == ["#abcdef"]
            assert last.qualifiers["ApEinfo_revcolor"] == ["#abcdef"]
            # Auto-label still lands.
            assert last.qualifiers["label"] == ["tagged"]

    async def test_modal_dispatches_annotate_action(
        self, tiny_record, isolated_library,
    ):
        """End-to-end: setting a selection_range, opening the modal,
        and clicking "Insert feature" should fire `_add_feature_result`
        with `action="annotate"` and our captured range, which lands
        a feature at exactly that range."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            n_before = len(app._current_record.features)
            # Drive the result callback directly — the button is
            # wired to dismiss with this exact dict shape.
            app._add_feature_result({
                "action": "annotate",
                "range":  (12, 24),
                "entry":  {
                    "name": "marked",
                    "feature_type": "misc_feature",
                    "strand": 1,
                    "qualifiers": {},
                },
            })
            assert len(app._current_record.features) == n_before + 1
            new = app._current_record.features[-1]
            assert int(new.location.start) == 12
            assert int(new.location.end)   == 24
            assert new.qualifiers.get("label") == ["marked"]


class TestCDSDivisibleByThreeGate:
    """A CDS feature must span a whole number of codons. The modal
    blocks a non-divisible-by-3 selection inline (so the user fixes
    the highlight), and `_annotate_with_feature` repeats the check
    so direct callers (agent-API) get the same gate. The check uses
    the SELECTION SPAN, wrap-aware, not the typed-sequence length —
    the feature is anchored to the bp range, not whatever's in the
    Sequence textarea."""

    async def test_helper_rejects_non_divisible_cds(
        self, tiny_record, isolated_library,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            with pytest.raises(ValueError, match="multiple of 3|divisible by 3"):
                app._annotate_with_feature(0, 10, {   # 10 bp — not %3
                    "name": "x", "feature_type": "CDS",
                    "strand": 1, "qualifiers": {},
                })

    async def test_helper_accepts_divisible_cds(
        self, tiny_record, isolated_library,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(0, 9, {   # 9 bp — passes
                "name": "x", "feature_type": "CDS",
                "strand": 1, "qualifiers": {},
            })
            assert app._current_record.features[-1].type == "CDS"

    async def test_helper_rejects_wrap_cds_when_total_span_indivisible(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 100), id="wrap_cds", name="wrap_cds",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Wrap span = (100 - 95) + 5 = 10, not %3.
            with pytest.raises(ValueError, match="divisible by 3"):
                app._annotate_with_feature(95, 5, {
                    "name": "wcds", "feature_type": "CDS",
                    "strand": 1, "qualifiers": {},
                })

    async def test_helper_accepts_wrap_cds_when_total_span_divisible(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 100), id="wrap_ok", name="wrap_ok",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Wrap span = (100 - 94) + 3 = 9, %3 == 0 → accept.
            app._annotate_with_feature(94, 3, {
                "name": "wcds", "feature_type": "CDS",
                "strand": 1, "qualifiers": {},
            })
            assert app._current_record.features[-1].type == "CDS"

    async def test_non_cds_unchecked_even_when_indivisible(
        self, tiny_record, isolated_library,
    ):
        """The gate applies ONLY to CDS — promoter / misc_feature / etc.
        accept any span length, including non-multiples of 3."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(0, 10, {
                "name": "p", "feature_type": "promoter",
                "strand": 1, "qualifiers": {},
            })
            last = app._current_record.features[-1]
            assert last.type == "promoter"
            assert int(last.location.end) - int(last.location.start) == 10

    async def test_modal_blocks_indivisible_cds_inline(
        self, tiny_record, isolated_library,
    ):
        """Open the modal with a 10-bp selection, set type=CDS, click
        Insert. The button handler should NOT dismiss — the inline
        status box shows the divisible-by-3 error instead."""
        from textual.widgets import Static
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._user_sel = (0, 10)   # 10 bp — not %3
            app.action_add_feature()
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert isinstance(modal, sc.AddFeatureModal)
            modal.query_one("#addfeat-name").value = "bad-cds"
            # Type already defaults to "CDS"; force it explicitly.
            modal.query_one("#addfeat-type").value = "CDS"
            await pilot.pause(0.05)
            # Click Insert → still on the modal, status box flagged.
            await pilot.click("#btn-addfeat-insert")
            await pilot.pause(0.2)
            # Modal still on screen (not dismissed).
            assert isinstance(app.screen, sc.AddFeatureModal)
            status_text = str(modal.query_one("#addfeat-status", Static)
                                  .render())
            assert "multiple of 3" in status_text


class TestPackerNewOnTop:
    """Newly added features stack on top of older overlapping features.
    The packer iterates `feats` in insertion order; older features pack
    first and land at the bottom (closest to DNA), newer features get
    pushed to higher rows wherever their column range overlaps an
    older feature. Per-feature priority rotation will land in a future
    release; until then, recency is the rule."""

    def test_overlapping_new_feature_lands_above_old(self):
        old = {"start": 0, "end": 30, "type": "misc_feature",
                "label": "old", "strand": 1, "color": "white"}
        new = {"start": 10, "end": 20, "type": "misc_feature",
                "label": "new", "strand": 1, "color": "white"}
        placements = sc._pack_features_2d([old, new], 0, 30)
        rows = {p[0]["label"]: p[1] for p in placements}
        assert rows["old"] < rows["new"], (
            f"new feature should land above old; got rows={rows}"
        )

    def test_new_cds_pushed_above_existing_non_cds(self):
        """Pre-fix, CDS features pre-empted lane 0 over non-CDS.
        Post-fix, insertion order rules: an existing non-CDS keeps
        lane 0 and a newly added CDS lands above it."""
        old_promoter = {"start": 0, "end": 30, "type": "promoter",
                          "label": "p", "strand": 1, "color": "white"}
        new_cds = {"start": 10, "end": 22, "type": "CDS",
                     "label": "c", "strand": 1, "color": "white"}
        placements = sc._pack_features_2d(
            [old_promoter, new_cds], 0, 30,
        )
        rows = {p[0]["label"]: p[1] for p in placements}
        assert rows["p"] == 0
        assert rows["c"] > 0

    def test_non_overlapping_features_share_lane(self):
        """Features that don't overlap pack into the same row regardless
        of insertion order — recency only kicks in on collisions."""
        a = {"start": 0,  "end": 10, "type": "misc_feature",
              "label": "a", "strand": 1, "color": "white"}
        b = {"start": 20, "end": 30, "type": "misc_feature",
              "label": "b", "strand": 1, "color": "white"}
        placements = sc._pack_features_2d([a, b], 0, 30)
        rows = {p[0]["label"]: p[1] for p in placements}
        assert rows["a"] == 0
        assert rows["b"] == 0


class TestGlyphOwnerTracking:
    """`SequencePanel._chunk_glyph_owners` fills per-(packed_row, col)
    owners across the FULL `(footprint_rows × bp_range)` rectangle of
    each feature — clicks anywhere on the visible lane art (bar,
    label-text, label-padding, AA cells) all resolve to the right
    feature. `_check_packed` does the codon-vs-bar dispatch for CDS
    afterward."""

    def _build_panel(self, seq: str, feats: list) -> "sc.SequencePanel":
        sp = sc.SequencePanel.__new__(sc.SequencePanel)
        sp._seq = seq
        sp._feats = feats
        sp._chunks_owners = {}
        # `_chunk_glyph_owners` keys the per-chunk owner cache on
        # `_view_origin_bp` so rotated views don't reuse stale lookup
        # tables. The bare-attribute SequencePanel built here doesn't
        # run `__init__`, so wire the rotation default in by hand.
        sp._view_origin_bp = 0
        return sp

    def test_owners_fill_full_footprint_for_non_cds(self):
        """A non-CDS feature occupies 2 packed rows (bar + label).
        owners_above[0..1][bp_range] should all be the feature, even
        in label-padding cells (so clicking on padding still picks
        the feature). Outside the bp range, owner=None."""
        f = {"start": 5, "end": 15, "type": "misc_feature",
              "label": "f", "strand": 1, "color": "white"}
        above_p, below_p, above_rows, below_rows = sc._chunk_lane_groups(
            [f], 0, 20,
        )
        sp = self._build_panel("A" * 20, [f])
        result = sp._chunk_glyph_owners(
            0, 20, [f], above_p, below_p, above_rows, below_rows,
        )
        owners = result["owners_above"]
        assert len(owners) == above_rows == 2
        for r in (0, 1):
            for col in range(5, 15):
                assert owners[r][col] is f, (
                    f"row {r} col {col} should own f"
                )
            for col in (0, 4, 15, 19):
                assert owners[r][col] is None

    def test_nested_features_get_distinct_owners_per_cell(self):
        """Smaller feature inside a larger one — every cell in the
        smaller's rectangle should own the smaller, every cell in
        the larger's rectangle (outside the smaller) should own the
        larger. Greedy packing ensures the rectangles don't overlap
        per row, so each cell has one unambiguous owner."""
        outer = {"start": 0, "end": 30, "type": "misc_feature",
                  "label": "outer", "strand": 1, "color": "white"}
        inner = {"start": 12, "end": 18, "type": "misc_feature",
                  "label": "inner", "strand": 1, "color": "red"}
        above_p, below_p, above_rows, below_rows = sc._chunk_lane_groups(
            [outer, inner], 0, 30,
        )
        sp = self._build_panel("A" * 30, [outer, inner])
        result = sp._chunk_glyph_owners(
            0, 30, [outer, inner],
            above_p, below_p, above_rows, below_rows,
        )
        owners = result["owners_above"]
        # outer at rows 0-1, inner at rows 2-3 (greedy packing pushes
        # inner above outer where bp ranges overlap).
        assert owners[0][5]  is outer
        assert owners[0][15] is outer
        assert owners[0][25] is outer
        assert owners[1][5]  is outer
        # Inner's rectangle: rows 2-3, cols 12-17.
        for r in (2, 3):
            for col in range(12, 18):
                assert owners[r][col] is inner
            for col in (5, 11, 18, 25):
                assert owners[r][col] is None

    def test_below_strand_owners_filled(self):
        """Reverse-strand feature in the below-DNA stack — `_check_packed`
        for `is_below=True` looks up `owners_below` indexed by
        screen_row_idx_from_top (= packed_row directly). Owners must
        be filled across the bp range for that strand, mirroring above."""
        rev = {"start": 5, "end": 25, "type": "misc_feature",
                "label": "rev", "strand": -1, "color": "white"}
        above_p, below_p, above_rows, below_rows = sc._chunk_lane_groups(
            [rev], 0, 30,
        )
        # Reverse-strand goes to below_p, not above_p.
        assert any(p[0] is rev for p in below_p)
        sp = self._build_panel("A" * 30, [rev])
        result = sp._chunk_glyph_owners(
            0, 30, [rev], above_p, below_p, above_rows, below_rows,
        )
        owners = result["owners_below"]
        assert len(owners) == below_rows == 2
        for r in (0, 1):
            for col in range(5, 25):
                assert owners[r][col] is rev
            for col in (0, 4, 25, 29):
                assert owners[r][col] is None

    def test_owner_cache_invalidates_on_feats_change(self):
        """Cache key includes `id(self._feats)` — reassigning `_feats`
        forces a fresh compute. Without this, post-annotate clicks
        could lookup stale owners that don't include the new feature."""
        f1 = {"start": 0, "end": 10, "type": "misc_feature",
               "label": "f1", "strand": 1, "color": "white"}
        above_p, below_p, above_rows, below_rows = sc._chunk_lane_groups(
            [f1], 0, 30,
        )
        sp = self._build_panel("A" * 30, [f1])
        first = sp._chunk_glyph_owners(
            0, 30, [f1], above_p, below_p, above_rows, below_rows,
        )
        # Reassign feats — simulating annotate's update_seq path,
        # but skip the .clear() to verify cache key is the safety net.
        f2 = {"start": 20, "end": 28, "type": "misc_feature",
               "label": "f2", "strand": 1, "color": "red"}
        sp._feats = [f1, f2]
        above_p2, below_p2, above_rows2, below_rows2 = (
            sc._chunk_lane_groups([f1, f2], 0, 30)
        )
        second = sp._chunk_glyph_owners(
            0, 30, [f1, f2],
            above_p2, below_p2, above_rows2, below_rows2,
        )
        assert second is not first, (
            "id-based cache key must invalidate on feats reassignment"
        )

    def test_cds_owners_fill_full_footprint(self):
        """A CDS occupies 3 packed rows (AA + bar + label). All 3 rows
        own the CDS across its bp range — owner-fill is uniform across
        the rectangle. The codon-vs-bar dispatch happens in
        `_check_packed` based on `packed_row - bottom_row`, not in
        the owner data itself, so AA-row inter-letter cells still
        own the CDS (and `_check_packed` resolves them as bar clicks)."""
        cds = {"start": 0, "end": 30, "type": "CDS",
                "strand": 1, "color": "white", "label": "c"}
        above_p, below_p, above_rows, below_rows = sc._chunk_lane_groups(
            [cds], 0, 30,
        )
        sp = self._build_panel("A" * 30, [cds])
        result = sp._chunk_glyph_owners(
            0, 30, [cds], above_p, below_p, above_rows, below_rows,
        )
        owners = result["owners_above"]
        assert above_rows == 3
        for r in (0, 1, 2):
            for col in range(30):
                assert owners[r][col] is cds, (
                    f"row {r} col {col} should own CDS"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# Modal surface — mount + gather
# ═══════════════════════════════════════════════════════════════════════════════

class TestAddFeatureModal:

    async def test_modal_mounts(self, tiny_record, isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
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

    # ── Sweep #30: unified members table wiring ───────────────
    async def test_members_table_seeded_with_one_row_on_mount(
        self, tiny_record, isolated_library,
    ):
        """The modal always opens with at least one row in the
        members table — the auto-built solo row (rs=0,
        re=len(sequence))."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal(prefill={
                "name": "lacZ", "sequence": "ATGAAATAG",
                "feature_type": "CDS", "strand": 1,
            }))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            assert isinstance(modal, sc.AddFeatureModal)
            assert len(modal._members) == 1
            assert modal._members[0]["rel_start"] == 0
            assert modal._members[0]["rel_end"]   == 9
            assert modal._members[0]["feature_type"] == "CDS"

    async def test_gather_single_row_yields_solo_entry(
        self, tiny_record, isolated_library,
    ):
        """A 1-row table saves as a solo library entry — no
        `is_group`, no `members`. Back-compat preserved."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = "single"
            modal.query_one("#addfeat-seq").text = "ATGAAA"
            entry = modal._gather()
            assert entry is not None
            assert "is_group" not in entry
            assert "members" not in entry
            assert entry["sequence"] == "ATGAAA"

    async def test_gather_after_split_yields_group_entry(
        self, tiny_record, isolated_library,
    ):
        """Splitting the 1-row table via `_split_member` produces
        a 2-row table, which `_gather` saves as a group entry."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = "split-test"
            modal.query_one("#addfeat-seq").text = "ACGTACGT"
            await pilot.pause()
            # Manually split row 0 at position 4 — same data path
            # that the Split button + SplitPositionPromptModal
            # exercises.
            modal._members = sc._split_member(modal._members, 0, 4)
            modal._refresh_table()
            entry = modal._gather()
            assert entry is not None
            assert entry.get("is_group") is True
            assert len(entry["members"]) == 2
            assert entry["members"][0]["rel_end"] == 4
            assert entry["members"][1]["rel_start"] == 4

    async def test_seq_shrink_drops_invalid_rows(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 hardening: when the user types a long
        sequence, splits it, then deletes most bases, rows whose
        rel_start ends up past the new sequence length get
        dropped from `_members` — the in-memory state stays
        validator-clean."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-seq").text = "A" * 30
            await pilot.pause()
            # Split into 3 rows: [0,10], [10,20], [20,30]
            modal._members = sc._split_member(modal._members, 0, 10)
            modal._members = sc._split_member(modal._members, 1, 20)
            assert len(modal._members) == 3
            # Shrink seq to 12 bp — row 2 (rs=20) and row 3 part
            # of row 1 (rs=10, re=20→12) should clamp; row 2 drops
            # entirely (rs=20 > 12).
            modal.query_one("#addfeat-seq").text = "A" * 12
            await pilot.pause()
            await pilot.pause(0.05)
            # Row 2 had rs=20 — must be dropped.
            for m in modal._members:
                assert int(m["rel_start"]) < 12
                assert int(m["rel_end"])  <= 12
                assert int(m["rel_start"]) < int(m["rel_end"])

    async def test_seq_shrink_to_zero_synthesises_solo(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 hardening: if the user deletes the entire
        sequence after splitting into multiple rows, all rows
        would be invalid → we synthesise a 1-row solo from the
        first row's metadata so the modal isn't stuck at 0
        members."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-seq").text = "ATGCATGC"
            await pilot.pause()
            modal._members = sc._split_member(modal._members, 0, 4)
            assert len(modal._members) == 2
            modal.query_one("#addfeat-seq").text = ""
            await pilot.pause()
            await pilot.pause(0.05)
            assert len(modal._members) == 1
            # Solo row spans the (now-empty) sequence.
            assert modal._members[0]["rel_start"] == 0
            assert modal._members[0]["rel_end"]   == 0

    async def test_seq_change_auto_resyncs_solo_row(
        self, tiny_record, isolated_library,
    ):
        """Typing in the sequence box updates row 0's rel_end to
        match — the user never has to manually edit rs/re for
        the simple solo case."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-seq").text = "ATGCATGC"
            await pilot.pause()
            await pilot.pause(0.05)
            assert modal._members[0]["rel_end"] == 8
            modal.query_one("#addfeat-seq").text = "ATGCATGCATGCATGC"
            await pilot.pause()
            await pilot.pause(0.05)
            assert modal._members[0]["rel_end"] == 16

    async def test_add_row_midpoint_split_produces_adjacent_rows(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30: `+ Add row` splits the selected row at its
        midpoint to produce two CONTIGUOUS sub-features — head.end
        == tail.start, no overlap, no gap. (User-reported: the
        original tail-append behaviour created a 1-bp overlap
        when the table tiled the full sequence.)"""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-seq").text = "ATGCATGCATGC"   # 12 bp
            await pilot.pause()
            await pilot.pause(0.05)
            assert len(modal._members) == 1
            assert modal._members[0]["rel_end"] == 12
            modal._on_add_row(None)
            assert len(modal._members) == 2
            head, tail = modal._members
            # Contiguous: head.end == tail.start, no overlap.
            assert head["rel_end"] == tail["rel_start"]
            # Spans the full sequence with no gap.
            assert head["rel_start"] == 0
            assert tail["rel_end"]   == 12
            # Midpoint of [0, 12] is 6.
            assert head["rel_end"] == 6

    async def test_add_row_refuses_width_lt_2(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30: a row of width 1 has no interior position
        for the midpoint split; Add row surfaces a status-line
        message and leaves the members list unchanged."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-seq").text = "AC"   # 2 bp
            await pilot.pause()
            modal._members = [{
                "rel_start": 0, "rel_end": 1,
                "feature_type": "misc_feature", "label": "tiny",
                "color": None, "strand": 1, "qualifiers": {},
                "description": "",
            }]
            modal._refresh_table()
            n_before = len(modal._members)
            modal._on_add_row(None)
            assert len(modal._members) == n_before   # no split

    async def test_add_row_caps_at_max_members(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 hardening: clicking `+ Add row` past
        `_MAX_GROUP_MEMBERS` surfaces the cap message + refuses
        to grow the table — even if a malicious script keeps
        firing the button event."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            # Seed the modal at the cap.
            seq = "A" * 200
            modal.query_one("#addfeat-seq").text = seq
            await pilot.pause()
            modal._members = [
                {"rel_start": i, "rel_end": i + 1,
                 "feature_type": "misc_feature",
                 "label": f"m{i}", "color": None,
                 "strand": 1, "qualifiers": {}, "description": ""}
                for i in range(sc._MAX_GROUP_MEMBERS)
            ]
            modal._refresh_table()
            assert len(modal._members) == sc._MAX_GROUP_MEMBERS
            modal._on_add_row(None)
            assert len(modal._members) == sc._MAX_GROUP_MEMBERS

    async def test_split_with_width_lt_2_refused(
        self, tiny_record, isolated_library,
    ):
        """A row of width 1 has no interior position; splitting
        is refused with a status-line message."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-seq").text = "AC"
            await pilot.pause()
            modal._members = [{
                "rel_start": 0, "rel_end": 1,
                "feature_type": "misc_feature", "label": "tiny",
                "color": None, "strand": 1, "qualifiers": {},
                "description": "",
            }]
            modal._refresh_table()
            # Should be a no-op (no SplitPositionPromptModal
            # pushed); the members list stays unchanged.
            n_before = len(modal._members)
            modal._on_split_row(None)
            await pilot.pause()
            assert len(modal._members) == n_before

    async def test_remove_last_row_refused(
        self, tiny_record, isolated_library,
    ):
        """Calling `_on_remove_row` when only 1 row remains keeps
        the table intact — every entry needs at least one
        member."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-seq").text = "ATG"
            await pilot.pause()
            assert len(modal._members) == 1
            modal._on_remove_row(None)
            assert len(modal._members) == 1

    async def test_prefill_group_entry_seeds_multi_row_table(
        self, tiny_record, isolated_library,
    ):
        """Opening the modal with a prefilled group entry seeds
        the members table from the entry's `members` list."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            prefill = {
                "name": "BsaI-pad",
                "sequence": "ACGTACGTACGT",
                "is_group": True,
                "members": [
                    {"rel_start": 0, "rel_end": 4,
                     "feature_type": "misc_feature",
                     "label": "pad", "color": "#888888",
                     "strand": 0, "qualifiers": {}, "description": ""},
                    {"rel_start": 4, "rel_end": 10,
                     "feature_type": "protein_bind",
                     "label": "BsaI", "color": "#ff3333",
                     "strand": 1, "qualifiers": {}, "description": ""},
                    {"rel_start": 10, "rel_end": 12,
                     "feature_type": "misc_feature",
                     "label": "OH", "color": "#00cc00",
                     "strand": 1, "qualifiers": {}, "description": ""},
                ],
            }
            app.push_screen(sc.AddFeatureModal(prefill=prefill))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            assert len(modal._members) == 3
            labels = [m["label"] for m in modal._members]
            assert labels == ["pad", "BsaI", "OH"]

    async def test_prefill_malformed_group_falls_back_to_solo(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 hardening: a malformed group prefill (e.g.
        member with out-of-range coords) falls back to a 1-row
        solo synthesised from the entry's top-level fields rather
        than crashing the modal."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            prefill = {
                "name": "bad", "sequence": "ACGT",
                "feature_type": "CDS", "strand": 1,
                "is_group": True,
                "members": [{"rel_start": 0, "rel_end": 99,
                             "feature_type": "x", "label": "boom"}],
            }
            app.push_screen(sc.AddFeatureModal(prefill=prefill))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            # Fallback: 1-row solo
            assert len(modal._members) == 1
            assert modal._members[0]["feature_type"] == "CDS"

    async def test_gather_rejects_oversized_sequence(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 hardening: sequence > `_MAX_FEATURE_SEQ_LEN`
        is refused with a status-line message; no entry returned.
        Defends against a 1 GB paste of base content via the
        Sequence TextArea."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = "huge"
            # Just past the cap.
            modal.query_one("#addfeat-seq").text = (
                "A" * (sc._MAX_FEATURE_SEQ_LEN + 100)
            )
            await pilot.pause()
            assert modal._gather() is None

    async def test_gather_sanitises_name_control_chars(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 hardening: pasted control chars (ANSI escape,
        null, bell) in the Name field get scrubbed before the
        entry is built. Defends against terminal-escape smuggling
        through the library JSON / .gb export."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = (
                "good\x00\x1b[31mEVIL\x1b[0m\x07tail"
            )
            modal.query_one("#addfeat-seq").text = "ATGAAA"
            entry = modal._gather()
            assert entry is not None
            assert "\x00" not in entry["name"]
            assert "\x1b" not in entry["name"]
            assert "\x07" not in entry["name"]
            assert "good" in entry["name"]
            assert "tail" in entry["name"]

    async def test_gather_rejects_name_with_only_control_chars(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 hardening: Name of only control chars
        scrubs to empty string → reject with clear error
        message (not silently lose the entry)."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            # All-control-char Name (no printable chars at all)
            # scrubs to empty string → reject. The `[31m`
            # printable suffix from a typical ANSI sequence would
            # survive and produce a non-empty (junk) name; the
            # test focuses on the all-stripped case.
            modal.query_one("#addfeat-name").value = (
                "\x00\x1b\x07\t\r\n\x08"
            )
            modal.query_one("#addfeat-seq").text = "ATG"
            assert modal._gather() is None

    async def test_gather_sanitises_description(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30: Description field gets `_sanitize_note`
        scrub (strips C0 control chars except \\t / \\n which
        preserve formatting)."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal.query_one("#addfeat-name").value = "x"
            modal.query_one("#addfeat-seq").text = "ATG"
            modal.query_one("#addfeat-desc").value = (
                "intro\x00\x1b[31mEVIL\x1b[0m\x07"
            )
            entry = modal._gather()
            assert entry is not None
            assert "\x00" not in entry["description"]
            assert "\x1b" not in entry["description"]

    async def test_per_row_strand_picker_updates_only_target_row(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30: clicking the strand cell of a sub-feature
        row opens `StrandPickerModal` scoped to THAT row. The
        callback writes the picked strand to ONLY that row's
        strand field — other rows untouched. This verifies the
        full per-row arrow-change flow the user asked for."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Seed the modal with a 3-row group via prefill so
            # every row has its own strand to flip.
            prefill = {
                "name": "tri", "sequence": "ATGCATGCATGC",
                "is_group": True,
                "members": [
                    {"rel_start": 0,  "rel_end": 4,
                     "feature_type": "misc_feature",
                     "label": "a", "color": "#ff0000",
                     "strand": 1, "qualifiers": {},
                     "description": ""},
                    {"rel_start": 4,  "rel_end": 8,
                     "feature_type": "misc_feature",
                     "label": "b", "color": "#00ff00",
                     "strand": 1, "qualifiers": {},
                     "description": ""},
                    {"rel_start": 8,  "rel_end": 12,
                     "feature_type": "misc_feature",
                     "label": "c", "color": "#0000ff",
                     "strand": 1, "qualifiers": {},
                     "description": ""},
                ],
            }
            app.push_screen(sc.AddFeatureModal(prefill=prefill))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            assert len(modal._members) == 3
            # Capture other-row strands so we can assert they
            # stay put when we change row 1.
            strand_row0_before = modal._members[0]["strand"]
            strand_row2_before = modal._members[2]["strand"]
            # Drive the per-row strand picker callback directly
            # (the modal-push path uses Textual's screen stack
            # which doesn't expose pickable callback wiring in
            # the test harness — but the picker's only effect on
            # success is `self._members[row].strand = X`, which
            # we can verify by emulating the callback). The
            # picker's own validation tests live in
            # `TestStrandPicker` (button → strand mapping).
            # Open the picker for row 1, then dismiss with
            # strand=-1. We push the picker, wait for it to
            # mount, then call its `_rev` button handler
            # directly.
            modal._open_per_row_strand_picker(1)
            await pilot.pause()
            await pilot.pause(0.05)
            picker = app.screen
            assert isinstance(picker, sc.StrandPickerModal)
            # Click the Reverse button.
            picker._rev(None)
            await pilot.pause()
            await pilot.pause(0.05)
            # Picker dismissed → callback fired → row 1's strand
            # is -1, OTHER rows untouched.
            assert modal._members[1]["strand"] == -1
            assert modal._members[0]["strand"] == strand_row0_before
            assert modal._members[2]["strand"] == strand_row2_before

    async def test_per_row_color_picker_writes_to_row(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 regression: picking a color in the picker
        and clicking Save writes that color to ONLY the picked
        row. Defends against the 'save defaults to auto' bug."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            prefill = {
                "name": "x", "sequence": "ATGCATGCATGC",
                "is_group": True,
                "members": [
                    {"rel_start": 0, "rel_end": 4,
                     "feature_type": "misc_feature",
                     "label": "a", "color": None,
                     "strand": 1, "qualifiers": {},
                     "description": ""},
                    {"rel_start": 4, "rel_end": 8,
                     "feature_type": "misc_feature",
                     "label": "b", "color": "#888888",
                     "strand": 1, "qualifiers": {},
                     "description": ""},
                ],
            }
            app.push_screen(sc.AddFeatureModal(prefill=prefill))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            modal._open_per_row_color_picker(1)
            await pilot.pause()
            await pilot.pause(0.05)
            picker = app.screen
            assert isinstance(picker, sc.ColorPickerModal)
            # Simulate a swatch click — directly set pending.
            picker._set_pending("#FF0000")
            # Simulate clicking Save.
            picker._save(None)
            await pilot.pause()
            await pilot.pause(0.05)
            # Picker dismissed → callback ran → row 1's color = #FF0000.
            assert modal._members[1]["color"] == "#FF0000"
            # Row 0 unchanged.
            assert modal._members[0]["color"] is None

    async def test_per_row_strand_picker_arrowless(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30: the picker's Arrowless button maps to
        strand=0 on the target row. Validates that all four
        button paths land valid strand values."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal())
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            # The default solo row has strand=1 (Forward sticky
            # or fall-through default). Change it via the
            # picker.
            modal._open_per_row_strand_picker(0)
            await pilot.pause()
            await pilot.pause(0.05)
            picker = app.screen
            assert isinstance(picker, sc.StrandPickerModal)
            picker._none(None)
            await pilot.pause()
            await pilot.pause(0.05)
            assert modal._members[0]["strand"] == 0

    async def test_paste_bomb_in_prefill_members_refused(
        self, tiny_record, isolated_library,
    ):
        """Sweep #30 hardening: a prefill with 1000 members blows
        past `_MAX_GROUP_MEMBERS`; the modal falls back to a 1-row
        solo rather than rendering a huge table."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            prefill = {
                "name": "huge", "sequence": "A" * 200,
                "feature_type": "CDS", "strand": 1,
                "is_group": True,
                "members": [
                    {"rel_start": i, "rel_end": i + 1,
                     "feature_type": "x", "label": f"m{i}"}
                    for i in range(1000)
                ],
            }
            app.push_screen(sc.AddFeatureModal(prefill=prefill))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            assert len(modal._members) == 1   # fallback to solo


# ═══════════════════════════════════════════════════════════════════════════════
# Save-to-library flow via _add_feature_result
# ═══════════════════════════════════════════════════════════════════════════════

class TestAddFeatureSelectionPrefill:
    """Alt+Shift+F (`action_add_feature`) checks the seq panel for an
    active multi-bp selection and pre-fills the modal's Sequence body
    with those bases verbatim. Saves the typical "select region →
    Ctrl+C → paste into modal" round-trip when adding a feature for
    a region the user just highlighted. (Add-feature moved off Ctrl+F
    — now Find-sequence — to Alt+Shift+F on 2026-06-01; not Alt+F,
    which terminals send as a cursor-motion code.)

    Selection sources covered: drag/Shift-click (`_user_sel`) and
    feature picks (`_sel_range`). Single-bp selections are NOT
    pre-filled — a click that lands on one base shouldn't be treated
    as a selection. Wrap-around selections (end < start) splice
    tail+head correctly. Pre-existing 2026-04-30 add-feature path."""

    async def test_user_sel_prefills_sequence(self, tiny_record,
                                                isolated_library):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            # tiny_record sequence starts with "ATGAAAGATCTGGAATTC..."
            sp._user_sel = (0, 9)   # "ATGAAAGAT"
            app.action_add_feature()
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            assert isinstance(modal, sc.AddFeatureModal)
            seq_box = modal.query_one("#addfeat-seq")
            assert seq_box.text == "ATGAAAGAT"

    async def test_sel_range_prefills_sequence(self, tiny_record,
                                                 isolated_library):
        """`_sel_range` is set when a feature is highlighted via map /
        sidebar / lane click. The Ctrl+F flow should still pick it up."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._sel_range = (3, 12)   # "AAAGATCTG"
            app.action_add_feature()
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            assert isinstance(modal, sc.AddFeatureModal)
            seq_box = modal.query_one("#addfeat-seq")
            assert seq_box.text == "AAAGATCTG"

    async def test_user_sel_takes_precedence_over_sel_range(
        self, tiny_record, isolated_library,
    ):
        """If both selections happen to be set (e.g., feature picked
        then user dragged a different region), the user's drag wins —
        same precedence Ctrl+C uses."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._user_sel  = (0, 6)    # "ATGAAA"
            sp._sel_range = (10, 20)  # different region
            app.action_add_feature()
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            seq_box = modal.query_one("#addfeat-seq")
            assert seq_box.text == "ATGAAA"

    async def test_single_bp_selection_does_not_prefill(
        self, tiny_record, isolated_library,
    ):
        """A 1-bp 'selection' is what a plain click produces; treat it
        as no selection so the modal opens with an empty Sequence box."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._user_sel = (5, 6)   # 1 bp
            app.action_add_feature()
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            seq_box = modal.query_one("#addfeat-seq")
            assert seq_box.text == ""

    async def test_no_selection_opens_empty_modal(
        self, tiny_record, isolated_library,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._user_sel  = None
            sp._sel_range = None
            app.action_add_feature()
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            seq_box = modal.query_one("#addfeat-seq")
            assert seq_box.text == ""

    async def test_wrap_selection_splices_tail_plus_head(
        self, isolated_library,
    ):
        """Wrap-around selections (end < start) should splice the tail
        [start, n) + head [0, end) — same convention as `_user_sel` set
        by `select_feature_range` for an origin-spanning feature."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        seq = "AAAACCCCGGGGTTTT"   # n = 16
        rec = SeqRecord(Seq(seq), id="wrap_sel", name="wrap_sel",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._user_sel = (12, 4)   # tail "TTTT" + head "AAAA"
            app.action_add_feature()
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            seq_box = modal.query_one("#addfeat-seq")
            assert seq_box.text == "TTTTAAAA"

    async def test_prefill_uppercases_lowercase_input(
        self, isolated_library,
    ):
        """The seq panel can hold lowercase bases (some users edit
        input as lowercase to mark introns / annotation overlays).
        Ctrl+F should normalise to uppercase like Ctrl+C does, since
        the modal's downstream validator expects ACGT/IUPAC uppercase."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("aaaaccccGGGG"), id="case_test",
                        name="case_test",
                        annotations={"molecule_type": "DNA"})
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._user_sel = (0, 8)
            app.action_add_feature()
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            seq_box = modal.query_one("#addfeat-seq")
            assert seq_box.text == "AAAACCCC"


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


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #29 — group-style library entries: validator + annotate path
# ═══════════════════════════════════════════════════════════════════════════════


class TestGroupValidator:
    """`_validate_group_members` normalises + sanity-checks the
    members list of a group library entry. Raises `ValueError` on
    unrecoverable shape errors so the save / annotate paths surface
    a clean error rather than deep-tracing into BioPython."""

    def test_well_formed_members_pass_through(self):
        members = [
            {"rel_start": 0, "rel_end": 4,
             "feature_type": "misc_feature", "label": "GCGC pad",
             "color": "#888888", "strand": 0},
            {"rel_start": 4, "rel_end": 10,
             "feature_type": "protein_bind", "label": "Esp3I",
             "color": "#FF3333", "strand": 1},
        ]
        out = sc._validate_group_members(members, 15)
        assert len(out) == 2
        assert out[0]["label"] == "GCGC pad"
        assert out[0]["strand"] == 0
        assert out[1]["label"] == "Esp3I"
        assert out[1]["strand"] == 1

    def test_members_sorted_by_rel_start(self):
        members = [
            {"rel_start": 10, "rel_end": 12, "feature_type": "x",
             "label": "z"},
            {"rel_start":  0, "rel_end":  4, "feature_type": "x",
             "label": "a"},
            {"rel_start":  4, "rel_end": 10, "feature_type": "x",
             "label": "m"},
        ]
        out = sc._validate_group_members(members, 12)
        assert [m["label"] for m in out] == ["a", "m", "z"]

    def test_empty_members_rejected(self):
        with pytest.raises(ValueError, match="non-empty"):
            sc._validate_group_members([], 10)

    def test_non_list_members_rejected(self):
        with pytest.raises(ValueError):
            sc._validate_group_members(None, 10)

    def test_zero_length_sequence_rejected(self):
        with pytest.raises(ValueError, match="positive"):
            sc._validate_group_members(
                [{"rel_start": 0, "rel_end": 1, "feature_type": "x"}],
                0,
            )

    def test_rel_start_must_be_less_than_rel_end(self):
        with pytest.raises(ValueError, match="rel_start"):
            sc._validate_group_members(
                [{"rel_start": 5, "rel_end": 5,
                  "feature_type": "x"}],
                10,
            )

    def test_rel_end_beyond_sequence_rejected(self):
        with pytest.raises(ValueError, match="rel_start"):
            sc._validate_group_members(
                [{"rel_start": 0, "rel_end": 100,
                  "feature_type": "x"}],
                10,
            )

    def test_negative_rel_start_rejected(self):
        with pytest.raises(ValueError, match="rel_start"):
            sc._validate_group_members(
                [{"rel_start": -1, "rel_end": 4,
                  "feature_type": "x"}],
                10,
            )

    def test_unparseable_rel_coords_rejected(self):
        # Numeric strings ("0", "4") are accepted by `int(...)`
        # and behave fine; only truly unparseable values like
        # lists / None should raise.
        with pytest.raises(ValueError, match="ints"):
            sc._validate_group_members(
                [{"rel_start": [0, 1], "rel_end": 4,
                  "feature_type": "x"}],
                10,
            )
        with pytest.raises(ValueError, match="ints"):
            sc._validate_group_members(
                [{"rel_end": 4, "feature_type": "x"}],   # rel_start missing → None
                10,
            )

    def test_non_dict_member_rejected(self):
        with pytest.raises(ValueError, match="dict"):
            sc._validate_group_members(
                ["not a dict"], 10,
            )

    def test_missing_feature_type_defaults_to_misc_feature(self):
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4}], 4,
        )
        assert out[0]["feature_type"] == "misc_feature"

    def test_invalid_strand_clamped_to_one(self):
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x", "strand": 99}], 4,
        )
        assert out[0]["strand"] == 1

    def test_non_string_color_dropped_to_none(self):
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x", "color": 123}], 4,
        )
        assert out[0]["color"] is None

    def test_overlapping_members_allowed(self):
        # Overlap is rare but legitimate (parent CDS with an
        # active-site sub-feature inside). Validator must NOT reject.
        out = sc._validate_group_members(
            [
                {"rel_start": 0, "rel_end": 12,
                 "feature_type": "CDS",     "label": "parent"},
                {"rel_start": 3, "rel_end": 6,
                 "feature_type": "domain",  "label": "active_site"},
            ], 12,
        )
        assert len(out) == 2

    def test_gaps_between_members_allowed(self):
        # Gaps = unannotated bases inside the group's sequence.
        # Legitimate (e.g. spacer N kept unannotated).
        out = sc._validate_group_members(
            [
                {"rel_start": 0, "rel_end": 4,
                 "feature_type": "misc_feature", "label": "a"},
                {"rel_start": 10, "rel_end": 14,
                 "feature_type": "misc_feature", "label": "b"},
            ], 14,
        )
        assert len(out) == 2

    # ── Sweep #29 hardening: paste-attack defenses ───────────
    def test_too_many_members_rejected(self):
        # 65 > _MAX_GROUP_MEMBERS (64) → raises so a paste of a
        # giant JSON fragment can't smuggle thousands of features
        # into the library file.
        members = [
            {"rel_start": i, "rel_end": i + 1,
             "feature_type": "x", "label": f"m{i}"}
            for i in range(65)
        ]
        with pytest.raises(ValueError, match="too many"):
            sc._validate_group_members(members, 100)

    def test_at_member_cap_passes(self):
        # Exactly 64 (the cap) should still pass — fence-post check.
        members = [
            {"rel_start": i, "rel_end": i + 1,
             "feature_type": "x", "label": f"m{i}"}
            for i in range(64)
        ]
        out = sc._validate_group_members(members, 100)
        assert len(out) == 64

    def test_control_chars_stripped_from_label(self):
        # Null byte, ESC, bell + other C0 control chars MUST be
        # scrubbed so a paste of binary blob / ANSI-escape-laden
        # text can't smuggle escape sequences into Rich Text or
        # GenBank export.
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x",
              "label": "good\x00\x1b[31mEVIL\x1b[0m\x07tail"}],
            4,
        )
        # _sanitize_label strips C0; the visible chars survive.
        lbl = out[0]["label"]
        assert "\x00" not in lbl
        assert "\x1b" not in lbl
        assert "\x07" not in lbl
        assert "good" in lbl and "tail" in lbl

    def test_oversized_label_truncated(self):
        # _MAX_GROUP_LABEL_LEN = 200; over-length labels capped.
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x",
              "label": "X" * 5000}], 4,
        )
        assert len(out[0]["label"]) <= 200

    def test_control_chars_stripped_from_color(self):
        # A "color" string of pure C0 control chars → scrubbed
        # to empty → dropped to None.
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x",
              "color": "\x1b\x00 \x07\t"}], 4,
        )
        assert out[0]["color"] is None

    def test_control_chars_in_partial_color_scrubbed(self):
        # Mixed control-char + printable chars → control chars
        # stripped; whatever survives is preserved (caller's
        # responsibility to validate hex shape — Rich accepts
        # named colours too so we don't lock it down here).
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x",
              "color": "\x1b\x00#FF0000\x07"}], 4,
        )
        assert "\x1b" not in (out[0]["color"] or "")
        assert "\x00" not in (out[0]["color"] or "")
        assert "\x07" not in (out[0]["color"] or "")
        assert "#FF0000" in (out[0]["color"] or "")

    def test_oversized_color_dropped(self):
        # Suspiciously long "color" (>_MAX_GROUP_COLOR_LEN=32) →
        # drop to None rather than store nonsense.
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x",
              "color": "#" + "abcdef" * 10}], 4,
        )
        assert out[0]["color"] is None

    def test_non_string_qualifier_key_dropped(self):
        # Qualifier dict may come from JSON with weird keys; ints
        # / lists / None as keys silently dropped (preserves valid
        # entries, refuses bogus ones).
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x",
              "qualifiers": {
                  123: ["bad-key-type"],
                  "gene": ["good"],
              }}], 4,
        )
        assert "gene" in out[0]["qualifiers"]
        assert 123 not in out[0]["qualifiers"]

    def test_control_chars_stripped_from_qualifier_keys_and_values(
        self,
    ):
        # Both keys and values get C0 control char scrub before
        # round-tripping to .gb export.
        out = sc._validate_group_members(
            [{"rel_start": 0, "rel_end": 4,
              "feature_type": "x",
              "qualifiers": {
                  "ge\x00ne": ["va\x1blue", "ok"],
                  "label": "fine",
              }}], 4,
        )
        quals = out[0]["qualifiers"]
        for k in quals:
            assert "\x00" not in k and "\x1b" not in k
        for vlist in quals.values():
            if isinstance(vlist, list):
                for v in vlist:
                    assert "\x00" not in v and "\x1b" not in v


class TestAddFeatureModalGroupParser:
    """`AddFeatureModal._parse_group_members_text` parses the members
    TextArea content into a `_validate_group_members`-ready list.
    Sweep #29 (2026-05-26): runs after the validator's hardening
    layer so the worst a paste-attack can do is raise ValueError
    (caught by the save path and surfaced as a status-line error).

    Direct call (staticmethod) — no app harness needed."""

    def test_simple_parse(self):
        text = (
            "0-4  GCGC #888888 0\n"
            "4-10 Esp3I #FF3333 1\n"
            "10-14 AATG  #FFCC00 1\n"
        )
        out = sc.AddFeatureModal._parse_group_members_text(text, 14)
        assert [m["label"] for m in out] == ["GCGC", "Esp3I", "AATG"]
        assert out[0]["strand"] == 0
        assert out[1]["color"] == "#FF3333"

    def test_blank_and_comment_lines_ignored(self):
        text = (
            "# this is a comment\n"
            "\n"
            "0-4 alpha #ff0000 1\n"
            "   \n"
            "# another comment\n"
            "4-8 beta #00ff00 -1\n"
        )
        out = sc.AddFeatureModal._parse_group_members_text(text, 8)
        assert len(out) == 2
        assert out[0]["label"] == "alpha"
        assert out[1]["strand"] == -1

    def test_empty_text_returns_empty(self):
        out = sc.AddFeatureModal._parse_group_members_text("", 10)
        assert out == []

    def test_whitespace_only_text_returns_empty(self):
        out = sc.AddFeatureModal._parse_group_members_text(
            "  \n\n  \n# comment\n", 10,
        )
        assert out == []

    def test_non_string_text_returns_empty(self):
        # Defensive: caller shouldn't pass non-str but the parser
        # mustn't crash if it does.
        out = sc.AddFeatureModal._parse_group_members_text(
            None, 10,  # type: ignore[arg-type]
        )
        assert out == []

    def test_oversized_total_text_rejected(self):
        # > 32 KB total → ValueError so a paste-bomb of 5 MB JSON
        # blob can't make us try to parse line-by-line for ages.
        huge = "0-1 x #ff0000 1\n" * 5000   # ~80 KB
        with pytest.raises(ValueError, match="too long"):
            sc.AddFeatureModal._parse_group_members_text(huge, 10)

    def test_oversized_single_line_rejected(self):
        # Single line > 1024 chars → ValueError. Defends against a
        # paste of a giant single-line blob (no newlines = no line
        # parsing = wouldn't hit per-line caps without this check).
        long_label = "X" * 2000
        text = f"0-4 {long_label} #ff0000 1\n"
        with pytest.raises(ValueError, match="too long"):
            sc.AddFeatureModal._parse_group_members_text(text, 4)

    def test_malformed_first_token_rejected(self):
        # Missing "-" in coords → ValueError so we don't silently
        # drop a row the user thought was valid.
        with pytest.raises(ValueError, match="rel_start"):
            sc.AddFeatureModal._parse_group_members_text(
                "garbage_no_dash alpha\n", 10,
            )

    def test_non_int_coords_rejected(self):
        with pytest.raises(ValueError, match="ints"):
            sc.AddFeatureModal._parse_group_members_text(
                "abc-def alpha\n", 10,
            )

    def test_out_of_range_coords_rejected(self):
        # rel_end > sequence_len → ValueError
        with pytest.raises(ValueError, match="not in"):
            sc.AddFeatureModal._parse_group_members_text(
                "0-100 alpha\n", 10,
            )

    def test_inverted_coords_rejected(self):
        # rel_start >= rel_end → ValueError
        with pytest.raises(ValueError, match="not in"):
            sc.AddFeatureModal._parse_group_members_text(
                "5-3 alpha\n", 10,
            )

    def test_invalid_strand_clamped_to_one(self):
        # Unknown strand → falls back to 1 (matches SplitFeatureModal
        # parser behavior).
        out = sc.AddFeatureModal._parse_group_members_text(
            "0-4 alpha #ff0000 99\n", 4,
        )
        assert out[0]["strand"] == 1

    def test_non_int_strand_clamped_to_one(self):
        # Non-int strand → falls back to 1.
        out = sc.AddFeatureModal._parse_group_members_text(
            "0-4 alpha #ff0000 notanint\n", 4,
        )
        assert out[0]["strand"] == 1

    def test_color_optional(self):
        # No color column → color None (Auto fallback at render).
        out = sc.AddFeatureModal._parse_group_members_text(
            "0-4 alpha\n", 4,
        )
        assert out[0]["color"] is None
        assert out[0]["label"] == "alpha"

    def test_label_optional(self):
        # No label column at all → empty label OK.
        out = sc.AddFeatureModal._parse_group_members_text(
            "0-4\n", 4,
        )
        assert out[0]["label"] == ""

    def test_control_chars_in_line_stripped(self):
        # ANSI / null in line content gets scrubbed at parse-time
        # (pre-line-cap), then validator does another scrub for
        # defense-in-depth on the resulting dict.
        out = sc.AddFeatureModal._parse_group_members_text(
            "0-4 al\x1b[31mpha #ff0000 1\n", 4,
        )
        assert "\x1b" not in out[0]["label"]

    def test_too_many_members_rejected_at_parse_time(self):
        # 65 valid lines → validator's _MAX_GROUP_MEMBERS check
        # rejects after parse-shape passes. Caller gets a clean
        # ValueError, not a half-built list.
        lines = "\n".join(
            f"{i}-{i+1} m{i}" for i in range(65)
        )
        with pytest.raises(ValueError, match="too many"):
            sc.AddFeatureModal._parse_group_members_text(lines, 100)


class TestIsGroupEntry:
    """`_is_group_entry` is the back-compat detector — `is_group=True`
    AND a non-empty `members` list. Malformed flags / empty members
    fall through to the single-feature path."""

    def test_legacy_single_feature_entry_not_a_group(self):
        assert sc._is_group_entry({
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATGAAA", "strand": 1,
        }) is False

    def test_well_formed_group_entry_is_a_group(self):
        assert sc._is_group_entry({
            "name": "g", "sequence": "ATGCATGC",
            "is_group": True,
            "members": [{"rel_start": 0, "rel_end": 4,
                          "feature_type": "x"}],
        }) is True

    def test_is_group_true_but_no_members_is_not_a_group(self):
        # Malformed entry — fall through to single-feature path
        # so the user at least sees the entry behave SOMEHOW.
        assert sc._is_group_entry({
            "name": "g", "sequence": "ATGCATGC",
            "is_group": True, "members": [],
        }) is False

    def test_non_dict_not_a_group(self):
        assert sc._is_group_entry(None) is False
        assert sc._is_group_entry("string") is False
        assert sc._is_group_entry([]) is False


class TestAnnotateGroup:
    """Group-aware annotate path: a group library entry pasted at
    bp X creates N sub-features sharing a fresh `feature_group`
    qualifier. Sacred invariants #8 / #9 covered explicitly."""

    @pytest.fixture
    def adapter_group(self):
        """Mirrors the user's GCGC + Esp3I + N + AATG cassette."""
        return {
            "name":     "Esp3I→AATG 5' adapter",
            "sequence": "GCGCCGTCTCNAATG",  # 15 bp
            "is_group": True,
            "members":  [
                {"rel_start":  0, "rel_end":  4,
                 "feature_type": "misc_feature",
                 "label": "GCGC pad", "color": "#888888",
                 "strand": 0},
                {"rel_start":  4, "rel_end": 10,
                 "feature_type": "protein_bind",
                 "label": "Esp3I",    "color": "#FF3333",
                 "strand": 1},
                {"rel_start": 10, "rel_end": 11,
                 "feature_type": "misc_feature",
                 "label": "N",        "color": "#666666",
                 "strand": 0},
                {"rel_start": 11, "rel_end": 15,
                 "feature_type": "misc_feature",
                 "label": "AATG",     "color": "#00CC66",
                 "strand": 1},
            ],
        }

    async def test_paste_lands_all_members_sharing_group_id(
        self, tiny_record, isolated_library, adapter_group,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            n0 = len(app._current_record.features)
            # tiny_record is 50 bp, selection 5..20 (15 bp = group len).
            app._annotate_with_feature(5, 20, adapter_group)
            n1 = len(app._current_record.features)
            assert n1 - n0 == 4
            # All four new features share one `feature_group` value.
            new_feats = app._current_record.features[-4:]
            group_ids = {
                f.qualifiers["feature_group"][0]
                for f in new_feats
            }
            assert len(group_ids) == 1
            # And the id is a 12-char hex (uuid4 prefix).
            (gid,) = group_ids
            assert len(gid) == 12
            assert all(c in "0123456789abcdef" for c in gid)

    async def test_members_land_at_correct_absolute_coords(
        self, tiny_record, isolated_library, adapter_group,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(5, 20, adapter_group)
            new_feats = app._current_record.features[-4:]
            # Members are sorted by rel_start in the validator, so
            # the LANDED order matches the input order: GCGC, Esp3I,
            # N, AATG.
            labels_to_coords = {
                f.qualifiers["label"][0]:
                    (int(f.location.start), int(f.location.end))
                for f in new_feats
            }
            assert labels_to_coords["GCGC pad"] == (5, 9)
            assert labels_to_coords["Esp3I"]    == (9, 15)
            assert labels_to_coords["N"]        == (15, 16)
            assert labels_to_coords["AATG"]     == (16, 20)

    async def test_per_member_color_persists_in_apeinfo(
        self, tiny_record, isolated_library, adapter_group,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(5, 20, adapter_group)
            new_feats = app._current_record.features[-4:]
            by_label = {f.qualifiers["label"][0]: f
                         for f in new_feats}
            assert by_label["GCGC pad"].qualifiers[
                "ApEinfo_fwdcolor"] == ["#888888"]
            assert by_label["Esp3I"].qualifiers[
                "ApEinfo_revcolor"] == ["#FF3333"]
            assert by_label["AATG"].qualifiers[
                "ApEinfo_fwdcolor"] == ["#00CC66"]

    async def test_per_member_strand_persists(
        self, tiny_record, isolated_library, adapter_group,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(5, 20, adapter_group)
            new_feats = app._current_record.features[-4:]
            by_label = {f.qualifiers["label"][0]: f
                         for f in new_feats}
            # Arrowless members: BioPython strand None.
            assert by_label["GCGC pad"].location.strand is None
            assert by_label["N"].location.strand is None
            # Forward members: BioPython strand 1.
            assert by_label["Esp3I"].location.strand == 1
            assert by_label["AATG"].location.strand == 1

    async def test_selection_span_mismatch_rejected(
        self, tiny_record, isolated_library, adapter_group,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            n0 = len(app._current_record.features)
            # Selection is 10 bp; group is 15 bp. Mismatch → error.
            with pytest.raises(ValueError, match="does not match"):
                app._annotate_with_feature(5, 15, adapter_group)
            # No partial mutation — feature count unchanged.
            assert len(app._current_record.features) == n0

    async def test_empty_group_sequence_rejected(
        self, tiny_record, isolated_library,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            with pytest.raises(ValueError):
                app._annotate_with_feature(0, 4, {
                    "name": "empty", "sequence": "",
                    "is_group": True,
                    "members": [{"rel_start": 0, "rel_end": 4,
                                  "feature_type": "x"}],
                })

    async def test_linear_plasmid_refuses_wrap_member(
        self, tiny_record, isolated_library, adapter_group,
    ):
        # tiny_record is a linear record (no `topology: circular`
        # annotation). A group whose end exceeds n must refuse.
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            n = len(app._current_record.seq)
            # Place anchor so the last AATG member would land at
            # bp n+something — i.e. need start such that
            # start + 15 > n. start = n - 5 → end = n+10 → wrap.
            start = n - 10
            end = start + 15  # > n
            # _feat_len gives n - start + 0 — but for linear we
            # don't allow wrap. The validator catches this with
            # `end > n` check.
            # First, the basic range check refuses end > n.
            with pytest.raises(ValueError):
                app._annotate_with_feature(start, end, adapter_group)

    async def test_circular_plasmid_member_past_origin_no_straddle(
        self, isolated_library, adapter_group,
    ):
        """Edge case: a member whose `raw_start` is already past
        the origin (raw_s > n) lands as a plain FeatureLocation,
        not a CompoundLocation — the WHOLE member lives past the
        origin, so there's no straddle to model. With n=30,
        anchor=24, the AATG member's raw range is [35, 39); both
        ends past origin, abs range [5, 9), no straddle."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import FeatureLocation
        rec = SeqRecord(
            Seq("A" * 30), id="circ", name="circ",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Selection 24..9, wrap-aware = 15 bp.
            app._annotate_with_feature(24, 9, adapter_group)
            new_feats = app._current_record.features[-4:]
            by_label = {f.qualifiers["label"][0]: f
                         for f in new_feats}
            aatg = by_label["AATG"]
            assert isinstance(aatg.location, FeatureLocation)
            assert int(aatg.location.start) == 5
            assert int(aatg.location.end)   == 9

    async def test_circular_plasmid_member_straddles_origin(
        self, isolated_library, adapter_group,
    ):
        """Anchor placed so the middle Esp3I member straddles
        the origin: raw_start < n, raw_end > n. The renderer
        needs a CompoundLocation with two parts so `_feat_bounds`
        re-extracts the right `(tail_start, head_end)` shape on
        the next `_parse`."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import CompoundLocation
        rec = SeqRecord(
            Seq("A" * 30), id="circ", name="circ",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # anchor=22 → Esp3I lands at (22+4=26, 22+10=32) →
            # straddles origin (26<30, 32>30). Selection 22..7
            # wrap-aware: 15 bp.
            app._annotate_with_feature(22, 7, adapter_group)
            new_feats = app._current_record.features[-4:]
            by_label = {f.qualifiers["label"][0]: f
                         for f in new_feats}
            esp = by_label["Esp3I"]
            assert isinstance(esp.location, CompoundLocation)
            parts = list(esp.location.parts)
            assert len(parts) == 2
            # First part: tail from raw_s=26 to n=30.
            assert int(parts[0].start) == 26
            assert int(parts[0].end)   == 30
            # Second part: head from 0 to abs_e=2.
            assert int(parts[1].start) == 0
            assert int(parts[1].end)   == 2

    async def test_group_undo_is_single_step(
        self, tiny_record, isolated_library, adapter_group,
    ):
        """Pasting a 4-member group must be a SINGLE undo step —
        Ctrl+Z reverts all 4 sub-features at once, not one at a
        time. The path uses `_push_undo` once before the bulk
        mutation."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            n0 = len(app._current_record.features)
            app._annotate_with_feature(5, 20, adapter_group)
            assert len(app._current_record.features) == n0 + 4
            app._action_undo()
            await pilot.pause()
            assert len(app._current_record.features) == n0

    async def test_strand_2_member_uses_splicecraft_qualifier(
        self, tiny_record, isolated_library,
    ):
        """A group member with strand=2 (double) MUST persist via
        the `SpliceCraft_strand=["double"]` qualifier so the
        `_paint_feature_bar` `strand == 2` branch fires on the
        next reload — same convention as single-feature path."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._annotate_with_feature(0, 6, {
                "name": "ds-group", "sequence": "ATGCAT",
                "is_group": True,
                "members": [{"rel_start": 0, "rel_end": 6,
                              "feature_type": "misc_feature",
                              "label": "ds",
                              "strand": 2}],
            })
            last = app._current_record.features[-1]
            assert last.location.strand is None
            assert last.qualifiers.get("SpliceCraft_strand") == [
                "double",
            ]
            # Round-trip via _parse: dict-side strand becomes 2.
            feats = sc.PlasmidMap._parse(
                sc.PlasmidMap.__new__(sc.PlasmidMap),
                app._current_record,
            )
            ds = next(f for f in feats if f["label"] == "ds")
            assert ds["strand"] == 2

    async def test_cds_member_must_be_codon_multiple(
        self, tiny_record, isolated_library,
    ):
        """Per-member CDS divisibility check mirrors the single-
        feature path. A group with a CDS member spanning 5 bp
        (not a codon multiple) refuses; no partial mutation."""
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            n0 = len(app._current_record.features)
            with pytest.raises(ValueError, match="codon"):
                app._annotate_with_feature(0, 9, {
                    "name": "bad-cds-group",
                    "sequence": "ATGCATGCA",
                    "is_group": True,
                    "members": [
                        {"rel_start": 0, "rel_end": 5,
                         "feature_type": "CDS",
                         "label": "shortCDS", "strand": 1},
                    ],
                })
            assert len(app._current_record.features) == n0


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #29 — group ops on canvas features + save-as-library-entry
# ═══════════════════════════════════════════════════════════════════════════════


class TestGroupOpsOnCanvas:
    """`_set_feature_group` / `_clear_feature_group` /
    `_feature_group_id` / `_features_in_group` are the helper APIs
    behind the sidebar Group / Ungroup actions (UI binding follows
    in a separate sweep). They operate on `PlasmidMap._feats`
    indices and persist via the `feature_group=[<uuid>]` qualifier."""

    async def test_set_feature_group_stamps_shared_id(
        self, isolated_library,
    ):
        """Three independent features → group all → all carry the
        same `feature_group` qualifier."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(
            Seq("A" * 100), id="r", name="r",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        rec.features = [
            SeqFeature(FeatureLocation(10, 20, strand=1),
                        type="CDS",
                        qualifiers={"label": ["a"]}),
            SeqFeature(FeatureLocation(30, 40, strand=1),
                        type="CDS",
                        qualifiers={"label": ["b"]}),
            SeqFeature(FeatureLocation(50, 60, strand=1),
                        type="CDS",
                        qualifiers={"label": ["c"]}),
        ]
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            gid = app._set_feature_group([0, 1, 2])
            assert gid is not None and len(gid) == 12
            for f in app._current_record.features:
                if f.type == "source":
                    continue
                assert f.qualifiers.get("feature_group") == [gid]

    async def test_clear_feature_group_drops_qualifier(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(
            Seq("A" * 100), id="r", name="r",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        rec.features = [
            SeqFeature(FeatureLocation(10, 20, strand=1),
                        type="CDS",
                        qualifiers={"label": ["a"],
                                     "feature_group": ["abc123"]}),
            SeqFeature(FeatureLocation(30, 40, strand=1),
                        type="CDS",
                        qualifiers={"label": ["b"],
                                     "feature_group": ["abc123"]}),
        ]
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            cleared = app._clear_feature_group([0, 1])
            assert cleared == 2
            for f in app._current_record.features:
                if f.type == "source":
                    continue
                assert "feature_group" not in f.qualifiers

    async def test_features_in_group_returns_all_members(
        self, isolated_library,
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(
            Seq("A" * 100), id="r", name="r",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        rec.features = [
            SeqFeature(FeatureLocation(10, 20, strand=1),
                        type="CDS",
                        qualifiers={"label": ["a"],
                                     "feature_group": ["g1"]}),
            SeqFeature(FeatureLocation(30, 40, strand=1),
                        type="CDS",
                        qualifiers={"label": ["b"],
                                     "feature_group": ["g2"]}),
            SeqFeature(FeatureLocation(50, 60, strand=1),
                        type="CDS",
                        qualifiers={"label": ["c"],
                                     "feature_group": ["g1"]}),
            SeqFeature(FeatureLocation(70, 80, strand=1),
                        type="CDS",
                        qualifiers={"label": ["d"]}),    # ungrouped
        ]
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            assert app._features_in_group("g1") == [0, 2]
            assert app._features_in_group("g2") == [1]
            assert app._features_in_group("nonexistent") == []
            assert app._features_in_group("") == []

    async def test_save_features_as_group_entry_round_trips(
        self, isolated_library,
    ):
        """End-to-end: build a 3-feature group on the canvas →
        save as library entry → load the entry → annotate at a
        new bp position → the three features land correctly with
        per-member colors/strands preserved."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(
            Seq("ATGCATGCATGCATGCATGC" + "A" * 80),  # 100 bp
            id="r", name="r",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        rec.features = [
            SeqFeature(
                FeatureLocation(0, 4, strand=0),
                type="misc_feature",
                qualifiers={"label": ["GCGC pad"],
                             "ApEinfo_fwdcolor": ["#888888"],
                             "ApEinfo_revcolor": ["#888888"]},
            ),
            SeqFeature(
                FeatureLocation(4, 10, strand=1),
                type="protein_bind",
                qualifiers={"label": ["Esp3I"],
                             "ApEinfo_fwdcolor": ["#FF3333"],
                             "ApEinfo_revcolor": ["#FF3333"]},
            ),
            SeqFeature(
                FeatureLocation(11, 15, strand=1),
                type="misc_feature",
                qualifiers={"label": ["AATG"],
                             "ApEinfo_fwdcolor": ["#00CC66"],
                             "ApEinfo_revcolor": ["#00CC66"]},
            ),
        ]
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            ok = app._save_features_as_group_entry(
                "test-adapter", [0, 1, 2],
            )
            assert ok is True
            # Library now has the entry.
            entries = sc._load_features()
            entry = next(e for e in entries
                          if e.get("name") == "test-adapter")
            assert entry["is_group"] is True
            assert len(entry["members"]) == 3
            # Anchor was bp 0; members keep their rel coords.
            members = sorted(entry["members"],
                              key=lambda m: m["rel_start"])
            assert (members[0]["rel_start"], members[0]["rel_end"]) == (0, 4)
            assert members[0]["label"] == "GCGC pad"
            assert members[0]["color"] == "#888888"
            assert members[0]["strand"] == 0
            assert (members[1]["rel_start"], members[1]["rel_end"]) == (4, 10)
            assert members[1]["color"] == "#FF3333"
            assert members[1]["strand"] == 1
            # Gap between members 1 and 2: rel 10..11 stays
            # unannotated by design.
            assert (members[2]["rel_start"], members[2]["rel_end"]) == (11, 15)
            assert members[2]["label"] == "AATG"
            # Annotate at a new position (bp 50).
            n0 = len(app._current_record.features)
            app._annotate_with_feature(50, 65, entry)
            n1 = len(app._current_record.features)
            assert n1 - n0 == 3
            new_feats = app._current_record.features[-3:]
            by_label = {f.qualifiers["label"][0]: f
                         for f in new_feats}
            # Re-anchored: members at 50 + rel.
            assert (int(by_label["GCGC pad"].location.start),
                    int(by_label["GCGC pad"].location.end)) == (50, 54)
            assert (int(by_label["Esp3I"].location.start),
                    int(by_label["Esp3I"].location.end)) == (54, 60)
            assert (int(by_label["AATG"].location.start),
                    int(by_label["AATG"].location.end)) == (61, 65)
            # All three share a fresh group id.
            gids = {f.qualifiers["feature_group"][0]
                     for f in new_feats}
            assert len(gids) == 1

    async def test_save_features_as_group_entry_refuses_wrap(
        self, isolated_library,
    ):
        """Per-member CompoundLocation (wrap-spanning feature)
        refuses the save — the user is told to rotate the plasmid
        first so the rel-coord math is unambiguous."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import (
            SeqFeature, FeatureLocation, CompoundLocation,
        )
        rec = SeqRecord(
            Seq("A" * 50), id="r", name="r",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        # Wrap feature: [45, 50) + [0, 5).
        rec.features = [
            SeqFeature(
                CompoundLocation([
                    FeatureLocation(45, 50, strand=1),
                    FeatureLocation(0, 5, strand=1),
                ]),
                type="misc_feature",
                qualifiers={"label": ["wrapfeat"]},
            ),
        ]
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            with pytest.raises(ValueError,
                                match="spans the origin"):
                app._save_features_as_group_entry(
                    "wrap-group", [0],
                )

    async def test_save_features_as_group_entry_empty_idx_refused(
        self, tiny_record, isolated_library,
    ):
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            with pytest.raises(ValueError, match="no features"):
                app._save_features_as_group_entry("g", [])
            with pytest.raises(ValueError, match="non-empty"):
                app._save_features_as_group_entry("", [0])

    async def test_set_group_idempotent_same_id(
        self, isolated_library,
    ):
        """Stamping a feature with the SAME group id is a no-op
        write (qualifier already has the value). Stamping with a
        DIFFERENT id replaces — most-recent wins."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(
            Seq("A" * 100), id="r", name="r",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        rec.features = [
            SeqFeature(FeatureLocation(10, 20, strand=1),
                        type="CDS",
                        qualifiers={"label": ["a"]}),
        ]
        app = sc.PlasmidApp()
        app._preload_record = rec
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            gid1 = app._set_feature_group([0], group_id="aaa111")
            assert gid1 == "aaa111"
            # Same id again: idempotent.
            gid2 = app._set_feature_group([0], group_id="aaa111")
            assert gid2 == "aaa111"
            # Different id: replaces.
            gid3 = app._set_feature_group([0], group_id="bbb222")
            assert gid3 == "bbb222"
            target = next(f for f in app._current_record.features
                           if f.type != "source")
            assert target.qualifiers["feature_group"] == ["bbb222"]


class TestGroupLibraryEntryRoundTrip:
    """Sweep #29 — a group library entry persisted to features.json
    survives `_save_features` → cache invalidation → `_load_features`
    intact: schema, members, colors, strands, descriptions."""

    def test_group_entry_round_trip_via_json(self, isolated_library):
        entry = {
            "name":     "test-rt",
            "sequence": "ATGCATGCATGCATGC",
            "feature_type": "misc_feature",
            "strand":   1,
            "is_group": True,
            "members":  [
                {"rel_start":  0, "rel_end":  4,
                 "feature_type": "misc_feature",
                 "label": "alpha", "color": "#FF0000",
                 "strand": 1, "qualifiers": {}, "description": ""},
                {"rel_start":  4, "rel_end": 16,
                 "feature_type": "CDS",
                 "label": "beta",  "color": "#00FF00",
                 "strand": -1, "qualifiers": {"gene": ["betagene"]},
                 "description": "test member"},
            ],
            "qualifiers": {},
            "description": "round-trip test",
        }
        sc._save_features([entry])
        sc._features_cache = None
        entries = sc._load_features()
        loaded = next(e for e in entries
                       if e.get("name") == "test-rt")
        assert loaded["is_group"] is True
        assert loaded["sequence"] == "ATGCATGCATGCATGC"
        assert len(loaded["members"]) == 2
        m0, m1 = loaded["members"]
        assert m0["label"] == "alpha"
        assert m0["color"] == "#FF0000"
        assert m1["strand"] == -1
        assert m1["qualifiers"]["gene"] == ["betagene"]

    def test_legacy_single_feature_entries_still_load(
        self, isolated_library,
    ):
        """Back-compat: a legacy single-feature entry (no
        `is_group` / no `members`) loads as before. The detector
        `_is_group_entry` returns False for it."""
        legacy = {
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATGAAA",
            "strand": 1, "qualifiers": {},
            "description": "",
        }
        sc._save_features([legacy])
        sc._features_cache = None
        entries = sc._load_features()
        loaded = next(e for e in entries
                       if e.get("name") == "lacZ")
        assert sc._is_group_entry(loaded) is False
        assert loaded["feature_type"] == "CDS"


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #30: unified split/merge helpers for the refactored modal
# ═══════════════════════════════════════════════════════════════════════════════

class TestSplitMember:
    """`_split_member(members, idx, pos)` returns a NEW members list
    where the row at `idx` is split into two contiguous rows at
    `pos`. Inherits metadata onto head + tail; tail gets blank
    label so user re-labels."""

    def _row(self, rs=0, re_=10, **kw):
        base = {
            "rel_start": rs, "rel_end": re_,
            "feature_type": "CDS", "label": "lacZ",
            "color": "#ff0000", "strand": 1,
            "qualifiers": {"gene": ["lacZ"]},
            "description": "the gene",
        }
        base.update(kw)
        return base

    def test_simple_split(self):
        out = sc._split_member([self._row(0, 10)], 0, 5)
        assert len(out) == 2
        assert out[0]["rel_start"] == 0
        assert out[0]["rel_end"]   == 5
        assert out[1]["rel_start"] == 5
        assert out[1]["rel_end"]   == 10

    def test_metadata_inherited_on_head_and_tail(self):
        out = sc._split_member(
            [self._row(0, 10, label="origCDS",
                       color="#abcdef", strand=-1)],
            0, 4,
        )
        assert out[0]["color"] == "#abcdef"
        assert out[1]["color"] == "#abcdef"
        assert out[0]["strand"] == -1
        assert out[1]["strand"] == -1
        assert out[0]["feature_type"] == "CDS"
        assert out[1]["feature_type"] == "CDS"

    def test_label_cleared_on_tail(self):
        out = sc._split_member([self._row(0, 10, label="parent")], 0, 5)
        assert out[0]["label"] == "parent"
        assert out[1]["label"] == ""

    def test_split_pos_at_rs_rejected(self):
        # pos == rs would create [0,0]+[0,10] — zero-width head.
        with pytest.raises(ValueError, match="must satisfy"):
            sc._split_member([self._row(0, 10)], 0, 0)

    def test_split_pos_at_re_rejected(self):
        # pos == re would create [0,10]+[10,10] — zero-width tail.
        with pytest.raises(ValueError, match="must satisfy"):
            sc._split_member([self._row(0, 10)], 0, 10)

    def test_split_pos_outside_range_rejected(self):
        with pytest.raises(ValueError, match="must satisfy"):
            sc._split_member([self._row(0, 10)], 0, 15)
        with pytest.raises(ValueError, match="must satisfy"):
            sc._split_member([self._row(5, 10)], 0, 2)

    def test_invalid_idx_rejected(self):
        with pytest.raises(ValueError, match="out of range"):
            sc._split_member([self._row()], 5, 5)
        with pytest.raises(ValueError, match="out of range"):
            sc._split_member([self._row()], -1, 5)

    def test_empty_members_rejected(self):
        with pytest.raises(ValueError, match="empty"):
            sc._split_member([], 0, 5)

    def test_split_in_middle_of_list(self):
        # Splitting row 1 of a 3-row list inserts the tail right
        # after row 1, shifting the rest by 1.
        rows = [
            self._row(0, 4, label="a"),
            self._row(4, 12, label="b"),
            self._row(12, 16, label="c"),
        ]
        out = sc._split_member(rows, 1, 8)
        assert len(out) == 4
        assert [m["label"] for m in out] == ["a", "b", "", "c"]
        assert out[1]["rel_end"] == 8
        assert out[2]["rel_start"] == 8
        assert out[2]["rel_end"] == 12

    def test_split_at_member_cap_rejected(self):
        # Sweep #30 hardening: even if a programmatic caller
        # bypasses the modal's cap check, `_split_member` itself
        # refuses to land an N+1 list past `_MAX_GROUP_MEMBERS`.
        members = [
            self._row(rs=i, re_=i+2, label=f"m{i}")
            for i in range(sc._MAX_GROUP_MEMBERS)
        ]
        # Each is 2bp wide so split-at-midpoint IS possible by
        # range alone — only the cap should reject.
        with pytest.raises(ValueError, match="_MAX_GROUP_MEMBERS"):
            sc._split_member(members, 0, members[0]["rel_start"] + 1)

    def test_split_returns_new_list_not_mutation(self):
        # Sacred invariant #4 (CLAUDE.md): never mutate lists in
        # place — return a fresh one so caller's chunk caches
        # invalidate on reassignment.
        rows = [self._row(0, 10)]
        out = sc._split_member(rows, 0, 5)
        assert out is not rows
        assert len(rows) == 1   # original untouched
        # And the dicts inside aren't the same objects either.
        assert out[0] is not rows[0]


class TestMergeMembers:
    """`_merge_members(members, idxs)` collapses adjacent rows into
    one spanning [min(rs), max(re)], inheriting the first (lowest-
    rs) row's metadata. Non-adjacent selections are rejected."""

    def _rows(self):
        return [
            {"rel_start":  0, "rel_end":  4,
             "feature_type": "x", "label": "a", "color": "#aa0000",
             "strand": 1, "qualifiers": {}, "description": ""},
            {"rel_start":  4, "rel_end":  8,
             "feature_type": "y", "label": "b", "color": "#00aa00",
             "strand": -1, "qualifiers": {}, "description": ""},
            {"rel_start":  8, "rel_end": 12,
             "feature_type": "z", "label": "c", "color": "#0000aa",
             "strand": 1, "qualifiers": {}, "description": ""},
        ]

    def test_simple_two_row_merge(self):
        out = sc._merge_members(self._rows(), [0, 1])
        assert len(out) == 2
        merged = out[0]
        assert merged["rel_start"] == 0
        assert merged["rel_end"]   == 8
        # First-row metadata wins (the leftmost picked row).
        assert merged["label"] == "a"
        assert merged["color"] == "#aa0000"
        assert merged["strand"] == 1

    def test_three_row_full_merge(self):
        out = sc._merge_members(self._rows(), [0, 1, 2])
        assert len(out) == 1
        assert out[0]["rel_start"] == 0
        assert out[0]["rel_end"]   == 12
        assert out[0]["label"] == "a"

    def test_non_adjacent_selection_rejected(self):
        # idxs [0, 2] skip row 1 → not adjacent.
        with pytest.raises(ValueError, match="adjacent"):
            sc._merge_members(self._rows(), [0, 2])

    def test_less_than_2_rows_rejected(self):
        with pytest.raises(ValueError, match="at least 2"):
            sc._merge_members(self._rows(), [0])

    def test_duplicate_idx_rejected(self):
        with pytest.raises(ValueError, match="unique"):
            sc._merge_members(self._rows(), [0, 0])

    def test_invalid_idx_rejected(self):
        with pytest.raises(ValueError, match="out of range"):
            sc._merge_members(self._rows(), [0, 99])

    def test_overlapping_rows_can_merge(self):
        # Two rows that overlap (e.g. parent CDS + active site
        # member) are still "adjacent" in sorted order → merge
        # absorbs both spans into the union.
        rows = [
            {"rel_start": 0, "rel_end": 12, "feature_type": "CDS",
             "label": "parent", "color": "#FF0000", "strand": 1,
             "qualifiers": {}, "description": ""},
            {"rel_start": 3, "rel_end":  6, "feature_type": "site",
             "label": "active", "color": "#0000FF", "strand": 0,
             "qualifiers": {}, "description": ""},
        ]
        out = sc._merge_members(rows, [0, 1])
        assert len(out) == 1
        assert out[0]["rel_start"] == 0
        assert out[0]["rel_end"]   == 12

    def test_merge_preserves_unselected_rows(self):
        rows = self._rows() + [
            {"rel_start": 12, "rel_end": 16,
             "feature_type": "q", "label": "d", "color": "#999999",
             "strand": 1, "qualifiers": {}, "description": ""},
        ]
        out = sc._merge_members(rows, [1, 2])
        # rows[0] + merged(rows[1+2]) + rows[3]
        assert len(out) == 3
        labels = [m["label"] for m in out]
        assert labels == ["a", "b", "d"]

    def test_merge_returns_new_list(self):
        rows = self._rows()
        out = sc._merge_members(rows, [0, 1])
        assert out is not rows
        assert len(rows) == 3  # original untouched


class TestEntryFromMembers:
    """`_entry_from_members(name, sequence, members, quals, desc)`
    collapses 1-row → solo entry (no `is_group`), 2+ row → group
    entry. Top-level qualifiers / description come from the modal's
    qualifier line + notes box."""

    def test_single_row_yields_solo_entry(self):
        members = [{
            "rel_start": 0, "rel_end": 6,
            "feature_type": "CDS", "label": "lacZ",
            "color": "#ff0000", "strand": 1,
            "qualifiers": {}, "description": "",
        }]
        entry = sc._entry_from_members(
            "lacZ", "ATGAAA", members,
            qualifiers={"gene": ["lacZ"]},
            description="the gene",
        )
        assert entry["name"]         == "lacZ"
        assert entry["sequence"]     == "ATGAAA"
        assert entry["feature_type"] == "CDS"
        assert entry["strand"]       == 1
        assert entry["color"]        == "#ff0000"
        assert entry["qualifiers"]   == {"gene": ["lacZ"]}
        assert entry["description"]  == "the gene"
        assert "is_group" not in entry
        assert "members" not in entry

    def test_two_row_yields_group_entry(self):
        members = [
            {"rel_start": 0, "rel_end": 4,
             "feature_type": "misc_feature", "label": "pad",
             "color": "#888888", "strand": 0,
             "qualifiers": {}, "description": ""},
            {"rel_start": 4, "rel_end": 10,
             "feature_type": "protein_bind", "label": "Esp3I",
             "color": "#ff3333", "strand": 1,
             "qualifiers": {}, "description": ""},
        ]
        entry = sc._entry_from_members(
            "BsaI-pad", "ACGTACGTAC", members,
        )
        assert entry["is_group"] is True
        assert len(entry["members"]) == 2
        assert entry["sequence"] == "ACGTACGTAC"

    def test_solo_inherits_row_metadata_for_solo_path(self):
        # The single row's qualifiers + description merge with the
        # modal-level qualifier dict; modal-level wins on key
        # conflicts.
        members = [{
            "rel_start": 0, "rel_end": 4,
            "feature_type": "promoter", "label": "p1",
            "color": "#00ff00", "strand": -1,
            "qualifiers": {"note": ["row-note"], "row_only": ["yes"]},
            "description": "row-desc",
        }]
        entry = sc._entry_from_members(
            "P1", "ACGT", members,
            qualifiers={"note": ["modal-note"]},
            description="modal-desc",
        )
        assert entry["feature_type"] == "promoter"
        assert entry["strand"]       == -1
        assert entry["color"]        == "#00ff00"
        # Modal-level "note" overrides the row's "note"; row's
        # "row_only" survives as it has no modal counterpart.
        assert entry["qualifiers"]["note"]     == ["modal-note"]
        assert entry["qualifiers"]["row_only"] == ["yes"]
        assert entry["description"] == "modal-desc"

    def test_empty_members_raises(self):
        with pytest.raises(ValueError):
            sc._entry_from_members("x", "ATG", [], qualifiers={})

    def test_solo_row_round_trip_via_solo_row_helper(self):
        # `_solo_row_from_entry` builds the 1-row shape from a
        # solo entry; `_entry_from_members` collapses it back.
        # Round-trip must be lossless for the obvious fields.
        original = {
            "name": "lacZ", "sequence": "ATGAAATAG",
            "feature_type": "CDS", "strand": 1,
            "color": "#ff0000",
            "qualifiers": {"gene": ["lacZ"]},
            "description": "the gene",
        }
        row = sc._solo_row_from_entry(original)
        entry = sc._entry_from_members(
            original["name"], original["sequence"], [row],
            qualifiers=original["qualifiers"],
            description=original["description"],
        )
        assert entry["name"]         == original["name"]
        assert entry["sequence"]     == original["sequence"]
        assert entry["feature_type"] == original["feature_type"]
        assert entry["strand"]       == original["strand"]
        assert entry["color"]        == original["color"]
        assert entry["qualifiers"]   == original["qualifiers"]
        assert entry["description"]  == original["description"]


# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #30: post-Ctrl+E group reconciliation
# ═══════════════════════════════════════════════════════════════════════════════

class TestGroupMemberCounter:
    """`_count_group_members(record)` walks the canvas record's
    feature list and returns `{group_id: count}` — used by the
    post-seq-edit reconciliation in `_edit_dialog_result` to
    detect groups that shrunk or collapsed to singletons."""

    def _record_with(self, qualifier_groups: "list[list[str]]"):
        """Build a small SeqRecord with N features, each carrying
        the given list of group_ids (empty list = no group)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(
            Seq("A" * 100), id="r", name="r",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        feats = []
        for i, gids in enumerate(qualifier_groups):
            quals = {}
            if gids:
                quals["feature_group"] = list(gids)
            feats.append(SeqFeature(
                FeatureLocation(i * 5, i * 5 + 3, strand=1),
                type="CDS",
                qualifiers={"label": [f"f{i}"], **quals},
            ))
        rec.features = feats
        return rec

    def test_empty_record_returns_empty(self):
        assert sc._count_group_members(None) == {}
        rec = self._record_with([])
        assert sc._count_group_members(rec) == {}

    def test_no_groups_returns_empty(self):
        rec = self._record_with([[], [], []])
        assert sc._count_group_members(rec) == {}

    def test_single_group(self):
        rec = self._record_with([["g1"], ["g1"], ["g1"]])
        assert sc._count_group_members(rec) == {"g1": 3}

    def test_multiple_groups(self):
        rec = self._record_with([
            ["g1"], ["g1"], ["g2"], ["g2"], ["g2"], [],
        ])
        assert sc._count_group_members(rec) == {"g1": 2, "g2": 3}

    def test_source_pseudo_feature_skipped(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(
            Seq("A" * 50), id="r", name="r",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        rec.features = [
            SeqFeature(FeatureLocation(0, 50, strand=1),
                        type="source",
                        qualifiers={"feature_group": ["g1"]}),
            SeqFeature(FeatureLocation(5, 10, strand=1),
                        type="CDS",
                        qualifiers={"feature_group": ["g1"]}),
        ]
        # Source feature with a (spurious) group qualifier is
        # NOT counted.
        assert sc._count_group_members(rec) == {"g1": 1}

    def test_empty_string_group_qualifier_skipped(self):
        rec = self._record_with([[""], ["g1"]])
        assert sc._count_group_members(rec) == {"g1": 1}

    def test_malformed_record_does_not_raise(self):
        # Defensive: malformed records should yield empty rather
        # than crash the seq-edit pipeline.
        class _Bogus:
            features = None  # not iterable
        assert sc._count_group_members(_Bogus()) == {}


class TestFeatureEditModalGroupMembersSeed:
    """Sweep #30 (2026-05-26) regression: when the caller passes
    `group_members` (the rel-coord projection of a whole canvas
    group), `FeatureEditModal.__init__` must seed `self._members`
    from them — even if their rel coords exceed the CURSOR
    feature's own length (the group's anchor span is wider than
    any single member). Original bug: validator's sequence_len
    arg was `_members_span_int()` (single feat) → group members
    rejected → silent fallback to 1-row solo → "no grouping
    behavior nor arrow change per segment"."""

    def test_group_members_seed_picks_group_span(self):
        # The cursor feature has width 6 (e.g. a 6-bp sub-feature
        # in a group). The group spans 30 bp total. Group members
        # have rel_ends up to 30 — way beyond the cursor's width.
        # Modal must accept all 4 members, not fall back to solo.
        group_members = [
            {"rel_start": 0,  "rel_end": 6,
             "feature_type": "misc_feature", "label": "pad",
             "color": "#888888", "strand": 0,
             "qualifiers": {}, "description": ""},
            {"rel_start": 6,  "rel_end": 12,
             "feature_type": "misc_feature", "label": "site",
             "color": "#ff0000", "strand": 1,
             "qualifiers": {}, "description": ""},
            {"rel_start": 12, "rel_end": 24,
             "feature_type": "CDS", "label": "core",
             "color": "#00ff00", "strand": -1,
             "qualifiers": {}, "description": ""},
            {"rel_start": 24, "rel_end": 30,
             "feature_type": "misc_feature", "label": "tail",
             "color": "#0000ff", "strand": 0,
             "qualifiers": {}, "description": ""},
        ]
        # Cursor feature: rep matching the second member (width 6).
        cursor_feat = {
            "start": 106, "end": 112, "strand": 1,
            "type": "misc_feature", "label": "site",
            "color": "#ff0000", "feature_group": "abc12345",
        }
        # Build the modal — total record length irrelevant for
        # this test, just needs to be ≥ cursor_end.
        modal = sc.FeatureEditModal(
            idx=5, feat=cursor_feat, total=2000,
            sequence="ATGCAT",
            group_members=group_members,
        )
        assert len(modal._members) == 4
        assert [m["label"] for m in modal._members] == [
            "pad", "site", "core", "tail",
        ]
        # The cursor's matching row should be selected by default.
        # `_find_self_idx` matches on (label, strand, width).
        assert modal._selected_idx == 1   # "site" row

    def test_group_members_seed_solo_fallback_when_none(self):
        # Caller passes None / empty → seed a 1-row solo from the
        # cursor feature itself (back-compat path).
        cursor_feat = {
            "start": 100, "end": 109, "strand": 1,
            "type": "CDS", "label": "lacZ", "color": "#ff0000",
        }
        modal = sc.FeatureEditModal(
            idx=2, feat=cursor_feat, total=500,
            sequence="ATGAAATAG",
            group_members=None,
        )
        assert len(modal._members) == 1
        assert modal._members[0]["label"] == "lacZ"
        assert modal._members[0]["rel_end"] == 9

    def test_group_members_seed_malformed_falls_back_to_solo(self):
        # Out-of-range rel coords (`rel_end > derived span`) →
        # validator rejects → fallback to 1-row solo. Defensive
        # path; the caller's normal flow validates first, so
        # this is a defense-in-depth.
        bad_members = [
            {"rel_start": 0, "rel_end": 10,
             "feature_type": "misc_feature", "label": "ok",
             "color": None, "strand": 1,
             "qualifiers": {}, "description": ""},
            {"rel_start": 5, "rel_end": 3,   # rel_end < rel_start
             "feature_type": "misc_feature", "label": "bad",
             "color": None, "strand": 1,
             "qualifiers": {}, "description": ""},
        ]
        cursor_feat = {
            "start": 0, "end": 10, "strand": 1,
            "type": "CDS", "label": "ok", "color": None,
        }
        modal = sc.FeatureEditModal(
            idx=0, feat=cursor_feat, total=100,
            sequence="ATGCATGCAT",
            group_members=bad_members,
        )
        # Fallback: 1-row solo from cursor feat.
        assert len(modal._members) == 1
        assert modal._members[0]["label"] == "ok"


class TestFeatureModalSpaceMark:
    """Universal Space-to-mark: in the feature-merge modal, Space toggles the
    ★ mark on the highlighted member row — the same gesture as the plasmid +
    primer libraries (so marking is consistent app-wide)."""

    class _SpaceKey:
        key = "space"

        def stop(self):
            pass

    async def test_space_toggles_member_mark(self, isolated_library):
        import splicecraft as sc
        from textual.widgets import DataTable
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 40)) as pilot:
            await pilot.pause()
            modal = sc.AddFeatureModal(
                prefill={"name": "m1", "sequence": "ACGTACGTACGT",
                         "feature_type": "CDS", "strand": 1},
                total_len=120)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            t = modal.query_one("#addfeat-members-tbl", DataTable)
            assert t.row_count >= 1
            t.focus()
            await pilot.pause()
            t.move_cursor(row=0)
            await pilot.pause()
            assert len(modal._marked_ids) == 0
            modal.on_key(self._SpaceKey())          # Space marks
            await pilot.pause()
            assert len(modal._marked_ids) == 1
            modal.on_key(self._SpaceKey())          # Space again unmarks
            await pilot.pause()
            assert len(modal._marked_ids) == 0
