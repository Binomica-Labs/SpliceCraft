"""Whole-record orientation edits: flip (reverse-complement) + re-origin.

[INV-181]. Two record-level transforms the user drives from the canvas:

  * ``_reverse_complement_record`` — flip the molecule end-for-end.
  * ``_rotate_seq_record(..., keep_source=True)`` — re-cut a circle so a
    chosen base becomes base 1.

The load-bearing property for BOTH is that **a feature's DNA does not
change**. For a stranded (±1) feature that means ``extract()`` returns the
byte-identical string before and after; for a directionless one (arrowless
0 / double) it means the same physical bases, read from the other strand.
That single check is what proves the coordinate mirror AND the strand flip
are right together — either one alone can look plausible and be wrong, which
is exactly how the exon-order bug in the first cut of this code survived
reading and was caught by the assertion below.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


# ── fixtures ──────────────────────────────────────────────────────────────

# Deliberately carries IUPAC codes (R/Y/K/M) so a complement table that only
# knows ACGT corrupts them visibly (sacred invariant #3).
SEQ = "ATGGGCCCRYKMTAGCTAGCTAGCATCGATCGGGTTTAAACCC"


def _mk(seq: str = SEQ, *, circular: bool = True):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    n = len(seq)
    r = SeqRecord(Seq(seq), id="t", name="t", description="test",
                  annotations={"topology": "circular" if circular else "linear",
                               "molecule_type": "DNA"})
    if n < len(SEQ):
        # The feature set below is laid out against the full-length SEQ;
        # a short/empty record gets a bare one for the degenerate cases.
        return r
    r.features = [
        SeqFeature(FeatureLocation(0, n, strand=1), type="source",
                   qualifiers={"label": ["source"], "organism": ["E. coli"]}),
        SeqFeature(FeatureLocation(3, 15, strand=1), type="CDS",
                   qualifiers={"label": ["fwd"]}),
        SeqFeature(FeatureLocation(20, 30, strand=-1), type="CDS",
                   qualifiers={"label": ["rev"]}),
        SeqFeature(FeatureLocation(31, 36, strand=0), type="misc_feature",
                   qualifiers={"label": ["arrowless"],
                               "SpliceCraft_strand": ["none"]}),
        # origin-spanning, forward
        SeqFeature(CompoundLocation([FeatureLocation(38, n, strand=1),
                                     FeatureLocation(0, 3, strand=1)]),
                   type="misc_feature", qualifiers={"label": ["wrap"]}),
        # spliced on the reverse strand — parts in BioPython's 5'->3'
        # (decreasing) order, as `complement(join(...))` parses to
        SeqFeature(CompoundLocation([FeatureLocation(29, 40, strand=-1),
                                     FeatureLocation(9, 20, strand=-1)]),
                   type="CDS", qualifiers={"label": ["rev-spliced"]}),
    ]
    return r


def _by_label(rec):
    return {(f.qualifiers.get("label") or [""])[0]: f for f in rec.features}


def _sig(rec):
    """Location + strand signature, for identity comparisons."""
    return [(f.type, str(f.location), f.location.strand) for f in rec.features]


# ── flip: the DNA under every feature is preserved ────────────────────────


class TestFlipKeepsFeatureDna:
    def test_stranded_features_extract_byte_identical(self):
        """THE direction check. A ±1 feature flipped correctly reads the
        same 5'->3' sequence afterwards — coordinates mirrored AND strand
        reversed. Mirror without the strand flip (or vice versa) fails."""
        r = _mk()
        f = sc._reverse_complement_record(r)
        for a, b in zip(r.features, f.features):
            if a.location.strand not in (1, -1):
                continue
            label = (a.qualifiers.get("label") or ["?"])[0]
            assert str(a.extract(r.seq)) == str(b.extract(f.seq)), (
                f"{label}: direction wrong after flip")

    def test_reverse_spliced_exon_order_survives(self):
        """Re-reversing an already-5'->3' part list scrambles exon order and
        still *looks* right — this is the assertion that catches it."""
        r = _mk()
        f = sc._reverse_complement_record(r)
        a = _by_label(r)["rev-spliced"]
        b = _by_label(f)["rev-spliced"]
        assert str(a.extract(r.seq)) == str(b.extract(f.seq))
        assert len(b.location.parts) == 2

    def test_origin_spanning_feature_survives(self):
        r = _mk()
        f = sc._reverse_complement_record(r)
        a, b = _by_label(r)["wrap"], _by_label(f)["wrap"]
        assert str(a.extract(r.seq)) == str(b.extract(f.seq))
        assert len(b.location.parts) == 2

    def test_directionless_feature_covers_the_same_bases(self):
        """An arrowless feature has no direction to reverse, so its strand
        stays 0 and its bases come back as the reverse complement — the same
        physical DNA, read from the other strand."""
        r = _mk()
        f = sc._reverse_complement_record(r)
        a, b = _by_label(r)["arrowless"], _by_label(f)["arrowless"]
        assert b.location.strand == 0
        assert str(b.extract(f.seq)) == sc._rc(str(a.extract(r.seq)))

    def test_iupac_codes_complement_correctly(self):
        r = _mk()
        f = sc._reverse_complement_record(r)
        assert str(f.seq) == sc._rc(SEQ)
        assert set(str(f.seq)) <= set("ACGTRYKMSWBDHVN")


class TestFlipInvariants:
    def test_double_flip_is_the_identity(self):
        r = _mk()
        back = sc._reverse_complement_record(sc._reverse_complement_record(r))
        assert str(back.seq) == SEQ
        assert _sig(back) == _sig(r)

    def test_length_and_feature_count_unchanged(self):
        r = _mk()
        f = sc._reverse_complement_record(r)
        assert len(f.seq) == len(r.seq)
        assert len(f.features) == len(r.features)

    def test_source_and_qualifiers_are_preserved(self):
        r = _mk()
        f = sc._reverse_complement_record(r)
        src = _by_label(f)["source"]
        assert src.type == "source"
        assert src.qualifiers.get("organism") == ["E. coli"]
        assert _by_label(f)["arrowless"].qualifiers.get(
            "SpliceCraft_strand") == ["none"]

    def test_the_input_record_is_not_mutated(self):
        r = _mk()
        before = _sig(r)
        sc._reverse_complement_record(r)
        assert _sig(r) == before
        assert str(r.seq) == SEQ

    def test_qualifiers_are_deep_copied(self):
        """A shallow copy shares the qualifier LISTS, so editing the flipped
        record's label would reach back into the original."""
        r = _mk()
        f = sc._reverse_complement_record(r)
        _by_label(f)["fwd"].qualifiers["label"][0] = "changed"
        assert (_by_label(r)["fwd"].qualifiers["label"][0] == "fwd")

    def test_topology_and_identity_carry_over(self):
        r = _mk()
        f = sc._reverse_complement_record(r)
        assert f.annotations.get("topology") == "circular"
        assert f.id == "t" and f.name == "t"


class TestFlipEdgeCases:
    def test_empty_sequence_is_returned_unchanged(self):
        r = _mk("")
        assert sc._reverse_complement_record(r) is r

    def test_single_base(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        r = _mk("A")
        r.features = [SeqFeature(FeatureLocation(0, 1, strand=1), type="CDS",
                                 qualifiers={"label": ["one"]})]
        f = sc._reverse_complement_record(r)
        assert str(f.seq) == "T"
        assert (int(f.features[0].location.start),
                int(f.features[0].location.end)) == (0, 1)
        assert f.features[0].location.strand == -1

    def test_a_record_with_no_features_flips_cleanly(self):
        r = _mk()
        r.features = []
        f = sc._reverse_complement_record(r)
        assert str(f.seq) == sc._rc(SEQ) and f.features == []

    def test_out_of_range_feature_is_dropped_not_guessed(self, caplog):
        """A coordinate we cannot mirror must be dropped loudly. Placing it
        at a guessed position would be worse than losing it."""
        import logging
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        r = _mk()
        r.features.append(SeqFeature(
            FeatureLocation(0, len(SEQ) + 50, strand=1), type="misc_feature",
            qualifiers={"label": ["off-the-end"]}))
        # `_log` sets propagate=False, so caplog's root handler never sees it.
        logger = logging.getLogger(sc._log.name)
        prev = logger.level
        logger.addHandler(caplog.handler)
        logger.setLevel(logging.WARNING)
        try:
            f = sc._reverse_complement_record(r)
        finally:
            logger.removeHandler(caplog.handler)
            logger.setLevel(prev)
        assert "off-the-end" not in _by_label(f)
        assert "off-the-end" in caplog.text

    def test_unresolvable_position_does_not_crash(self):
        from Bio.SeqFeature import SeqFeature, SimpleLocation, UnknownPosition
        r = _mk()
        r.features.append(SeqFeature(
            SimpleLocation(UnknownPosition(), UnknownPosition(), strand=1),
            type="misc_feature", qualifiers={"label": ["unknown"]}))
        f = sc._reverse_complement_record(r)     # must not raise
        assert "unknown" not in _by_label(f)


class TestFlipStrandHelper:
    @pytest.mark.parametrize("given,want", [(1, -1), (-1, 1), (0, 0),
                                            (None, None)])
    def test_only_real_directions_flip(self, given, want):
        assert sc._flip_feature_strand(given) == want


# ── re-origin (rotate) ────────────────────────────────────────────────────


class TestReOrigin:
    def test_stranded_features_keep_their_dna(self):
        r = _mk()
        rot = sc._rotate_seq_record(r, 17, keep_source=True)
        for lbl in ("fwd", "rev", "wrap"):
            a, b = _by_label(r)[lbl], _by_label(rot)[lbl]
            assert str(a.extract(r.seq)) == str(b.extract(rot.seq)), lbl

    def test_rotating_all_the_way_round_is_the_identity(self):
        r = _mk()
        n = len(SEQ)
        back = sc._rotate_seq_record(
            sc._rotate_seq_record(r, 17, keep_source=True), n - 17,
            keep_source=True)
        assert str(back.seq) == SEQ

    def test_keep_source_preserves_the_source_feature(self):
        r = _mk()
        rot = sc._rotate_seq_record(r, 11, keep_source=True)
        assert "source" in _by_label(rot)

    def test_default_still_drops_source(self):
        """The alignment path relies on this; the new flag must not change
        the behaviour every existing caller already depends on."""
        r = _mk()
        rot = sc._rotate_seq_record(r, 11)
        assert "source" not in _by_label(rot)

    def test_a_full_length_feature_is_not_split_at_the_new_origin(self):
        """`source` spans the whole record, and a full-circle span has no
        meaningful break point — splitting it produces a degenerate
        two-part location that renders as a bogus wrap."""
        r = _mk()
        rot = sc._rotate_seq_record(r, 13, keep_source=True)
        src = _by_label(rot)["source"]
        assert len(src.location.parts) == 1
        assert (int(src.location.start), int(src.location.end)) == (0, len(SEQ))

    def test_zero_offset_returns_the_record_unchanged(self):
        r = _mk()
        assert sc._rotate_seq_record(r, 0, keep_source=True) is r


# ── the canvas actions (driven through a real app) ────────────────────────


class TestReOriginKeepsFeatureShape:
    """[INV-182] Re-origin used to run every feature through `_feat_bounds`,
    which flattens anything that isn't a plain two-part origin wrap to
    ``(min start, max end)``. A spliced CDS therefore came back as ONE
    contiguous block with its introns swallowed, and a minus-strand feature
    pushed across the new origin came back with its two halves in ascending
    order — right bases, wrong reading order, i.e. a cyclically rotated
    protein once exported. Both are invisible unless you compare `extract()`.
    """

    @pytest.mark.parametrize("offset", [1, 5, 11, 17, 23, 30, 37, 41])
    def test_every_stranded_feature_keeps_its_dna_at_every_offset(self, offset):
        r = _mk()
        rot = sc._rotate_seq_record(r, offset, keep_source=True)
        before, after = _by_label(r), _by_label(rot)
        for lbl in ("fwd", "rev", "wrap", "rev-spliced"):
            assert str(before[lbl].extract(r.seq)) == \
                str(after[lbl].extract(rot.seq)), f"{lbl} @ offset {offset}"

    @pytest.mark.parametrize("offset", [1, 5, 11, 17, 23, 30, 37, 41])
    def test_a_spliced_feature_keeps_its_exon_count(self, offset):
        """An intron is not an annotation detail — a CDS that absorbs one
        translates straight through it."""
        r = _mk()
        rot = sc._rotate_seq_record(r, offset, keep_source=True)
        got = _by_label(rot)["rev-spliced"]
        # 11 + 11 exonic bases. The pre-fix code emitted ONE part spanning
        # 9..40 = 31 bases: the 9 bp intron swallowed into the CDS. (An exon
        # that lands on the new origin legitimately splits into two pieces,
        # so the part COUNT varies with the offset — the exonic base count
        # does not.)
        assert sum(int(p.end) - int(p.start)
                   for p in got.location.parts) == 22, (
            f"offset {offset}: {got.location} — exonic length changed, so "
            f"the intron between the exons is now part of the CDS")
        assert len(got.location.parts) >= 2

    def test_a_minus_feature_pushed_across_the_origin_reads_head_first(self):
        r = _mk()
        n = len(SEQ)
        # "rev" sits at [20, 30); rotating by 25 puts the new origin inside it.
        rot = sc._rotate_seq_record(r, 25, keep_source=True)
        got = _by_label(rot)["rev"]
        parts = [(int(p.start), int(p.end)) for p in got.location.parts]
        assert len(parts) == 2 and parts[0][0] == 0 and parts[1][1] == n, parts
        assert str(_by_label(r)["rev"].extract(r.seq)) == \
            str(got.extract(rot.seq))

    def test_the_rotated_record_survives_a_genbank_round_trip(self):
        """The whole point of the part order: what a THIRD-PARTY tool reads
        back out of the file we write."""
        import io

        from Bio import SeqIO

        r = _mk()
        rot = sc._rotate_seq_record(r, 25, keep_source=True)
        back = SeqIO.read(io.StringIO(sc._record_to_gb_text(rot)), "genbank")
        want = {lbl: str(f.extract(r.seq)) for lbl, f in _by_label(r).items()}
        for lbl, f in _by_label(back).items():
            if lbl in ("source", "arrowless"):
                continue          # whole-record / directionless by design
            assert str(f.extract(back.seq)) == want[lbl], lbl


class TestFeatureLocationBuilder:
    """`_feature_location` is the one place that decides which half of an
    origin-wrap comes first. [INV-182]"""

    def test_plus_strand_wrap_is_tail_then_head(self):
        loc = sc._feature_location(90, 10, 100, 1)
        assert [(int(p.start), int(p.end)) for p in loc.parts] == \
            [(90, 100), (0, 10)]

    def test_minus_strand_wrap_is_head_then_tail(self):
        loc = sc._feature_location(90, 10, 100, -1)
        assert [(int(p.start), int(p.end)) for p in loc.parts] == \
            [(0, 10), (90, 100)]

    def test_it_round_trips_through_feat_bounds(self):
        from Bio.SeqFeature import SeqFeature
        for strand in (1, -1, 0):
            for start, end in ((10, 40), (90, 10), (0, 100)):
                loc = sc._feature_location(start, end, 100, strand or None)
                got = sc._feat_bounds(SeqFeature(loc, type="CDS"), 100)
                assert got == (start, end, strand), (start, end, strand, got)

    def test_a_wrap_ending_on_the_origin_collapses_to_one_part(self):
        loc = sc._feature_location(90, 0, 100, -1)
        assert (int(loc.start), int(loc.end)) == (90, 100)
        assert len(loc.parts) == 1

    def test_a_zero_length_span_is_refused(self):
        assert sc._feature_location(40, 40, 100, 1) is None

    def test_minus_strand_wrap_exports_the_protein_it_displays(self):
        """The bug in one assertion: SpliceCraft's own `_translate_cds` reads
        the dict model and was always right; the GenBank writer read the part
        ORDER and was not."""
        import io

        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature

        seq = "ATGGGCCCATTAGCTAGCTAGCATCGATCGGGTTTAAACCCGGGTTTAAA"
        n = len(seq)
        rec = SeqRecord(Seq(seq), id="w", name="w",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec.features.append(SeqFeature(
            sc._feature_location(44, 6, n, -1), type="CDS",
            qualifiers={"label": ["wrapCDS"]}))
        shown = sc._translate_cds(seq, 44, 6, -1)
        back = SeqIO.read(io.StringIO(sc._record_to_gb_text(rec)), "genbank")
        exported = str(Seq(str(back.features[0].extract(back.seq))).translate())
        assert exported == shown, (exported, shown)


@pytest.mark.asyncio
class TestOrientationActions:
    """The circular gate and the undo contract, exercised through the real
    `PlasmidApp` — the actions query live widgets, so a hand-rolled stub
    would prove nothing about the thing that actually runs."""

    @staticmethod
    async def _app(tiny_record, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        return _build_app(tiny_record, isolated_library), TERMINAL_SIZE

    async def test_flip_reverse_complements_the_canvas(
            self, tiny_record, isolated_library):
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            before = str(app._current_record.seq)
            app.action_flip_sequence()
            await pilot.pause()
            assert str(app._current_record.seq) == sc._rc(before)
            app.exit()

    async def test_flip_is_undoable(self, tiny_record, isolated_library):
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            before = str(app._current_record.seq)
            app.action_flip_sequence()
            await pilot.pause()
            app.action_undo()
            await pilot.pause()
            assert str(app._current_record.seq) == before
            app.exit()

    async def test_flip_keeps_every_stranded_feature_reading_the_same(
            self, tiny_record, isolated_library):
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            old = app._current_record
            old_dna = [(f.type, str(f.extract(old.seq)))
                       for f in old.features if f.location.strand in (1, -1)]
            app.action_flip_sequence()
            await pilot.pause()
            new = app._current_record
            new_dna = [(f.type, str(f.extract(new.seq)))
                       for f in new.features if f.location.strand in (1, -1)]
            assert old_dna == new_dna
            app.exit()

    async def test_set_origin_refuses_on_a_linear_record(
            self, tiny_record, isolated_library):
        """Rotating a linear molecule would join two ends that aren't
        joined — it must refuse, and leave the record untouched."""
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            app._current_record.annotations["topology"] = "linear"
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 10
            before = str(app._current_record.seq)
            app.action_set_origin_here()
            await pilot.pause()
            assert str(app._current_record.seq) == before
            app.exit()

    async def test_set_origin_rotates_a_circular_record(
            self, tiny_record, isolated_library):
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            app._current_record.annotations["topology"] = "circular"
            before = str(app._current_record.seq)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            sp._cursor_pos = 10
            app.action_set_origin_here()
            await pilot.pause()
            after = str(app._current_record.seq)
            assert after == before[10:] + before[:10]
            assert len(after) == len(before)
            app.exit()

    async def test_set_origin_works_from_LINEAR_VIEW_when_circular(
            self, tiny_record, isolated_library):
        """The whole point of the feature: the map is drawn as a line, the
        MOLECULE is still a circle, so re-origining is legitimate. The gate
        must read topology, never the view toggle."""
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            app._current_record.annotations["topology"] = "circular"
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            await pilot.pause()
            before = str(app._current_record.seq)
            app.query_one("#seq-panel", sc.SequencePanel)._cursor_pos = 7
            app.action_set_origin_here()
            await pilot.pause()
            assert str(app._current_record.seq) == before[7:] + before[:7]
            app.exit()

    async def test_a_linear_record_in_circular_VIEW_is_still_refused(
            self, tiny_record, isolated_library):
        """The mirror of the above — the view saying 'circular' must not be
        able to authorise rotating a linear molecule."""
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            app._current_record.annotations["topology"] = "linear"
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "circular"
            await pilot.pause()
            before = str(app._current_record.seq)
            app.query_one("#seq-panel", sc.SequencePanel)._cursor_pos = 7
            app.action_set_origin_here()
            await pilot.pause()
            assert str(app._current_record.seq) == before
            app.exit()

    async def test_set_origin_at_base_one_is_a_no_op(
            self, tiny_record, isolated_library):
        """No rotation to do — and it must not burn an undo step that would
        then appear to do nothing."""
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            app._current_record.annotations["topology"] = "circular"
            app.query_one("#seq-panel", sc.SequencePanel)._cursor_pos = 0
            before = str(app._current_record.seq)
            depth = len(app._undo_stack)
            app.action_set_origin_here()
            await pilot.pause()
            assert str(app._current_record.seq) == before
            assert len(app._undo_stack) == depth
            app.exit()

    async def test_neither_action_crashes_with_nothing_loaded(
            self, tiny_record, isolated_library):
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            app._current_record = None
            app.action_flip_sequence()        # must not raise
            app.action_set_origin_here()      # must not raise
            await pilot.pause()
            assert app._current_record is None
            app.exit()

    async def test_reframing_clears_stale_alignment_bands(
            self, tiny_record, isolated_library):
        """Every coordinate moves, so a band drawn against the old frame now
        points at the wrong bases."""
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            app._alignments.append({"id": "a1", "label": "x"})
            app.action_flip_sequence()
            await pilot.pause()
            assert not app._alignments
            app.exit()

    async def test_a_reframed_record_does_not_persist_stale_alignments(
            self, tiny_record, isolated_library):
        """`add_entry` carries a same-id entry's stored alignments forward.
        After a flip those coordinates are wrong, so they must NOT ride along
        into the saved entry — they'd hydrate back as a lie on next load."""
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            lib = app.query_one(sc.LibraryPanel)
            lib.add_entry(app._current_record)
            await pilot.pause(0.1)
            rid = app._current_record.id
            entries = sc._load_library()
            for e in entries:
                if e.get("id") == rid:
                    e["alignments"] = [{"id": "a1", "label": "stale"}]
            sc._save_library(entries)
            await pilot.pause(0.1)
            app.action_flip_sequence()
            await pilot.pause()
            lib.add_entry(app._current_record)
            await pilot.pause(0.1)
            saved = sc._find_library_entry_by_id(rid)
            assert saved is not None
            assert not saved.get("alignments")
            app.exit()

    async def test_an_unflipped_resave_still_keeps_its_alignments(
            self, tiny_record, isolated_library):
        """The guard must be scoped to re-framing — an ordinary re-save
        still preserves alignments (the v1.2.50 behaviour)."""
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            lib = app.query_one(sc.LibraryPanel)
            lib.add_entry(app._current_record)
            await pilot.pause(0.1)
            rid = app._current_record.id
            entries = sc._load_library()
            for e in entries:
                if e.get("id") == rid:
                    e["alignments"] = [{"id": "a1", "label": "keep"}]
            sc._save_library(entries)
            await pilot.pause(0.1)
            lib.add_entry(app._current_record)
            await pilot.pause(0.1)
            saved = sc._find_library_entry_by_id(rid)
            assert saved and saved.get("alignments")
            app.exit()

    async def test_a_transform_that_changed_the_length_is_aborted(
            self, tiny_record, isolated_library, monkeypatch):
        """The safety net: if the transform ever regressed and returned a
        different amount of DNA, the record must be left ALONE rather than
        the corrupt result installed — and the undo snapshot pushed a
        moment earlier must not be left behind as a dud Ctrl+Z."""
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            before = str(app._current_record.seq)
            depth = len(app._undo_stack)

            real = sc._reverse_complement_record      # capture before patching

            def _truncating(rec):
                bad = real(rec)
                bad.seq = bad.seq[:-5]          # simulate a regression
                return bad
            monkeypatch.setattr(sc, "_reverse_complement_record", _truncating)
            app.action_flip_sequence()
            await pilot.pause()

            assert str(app._current_record.seq) == before, "record was corrupted"
            assert len(app._undo_stack) == depth, "left a dud undo step"
            app.exit()

    async def test_pop_undo_silently_is_safe_on_an_empty_stack(
            self, tiny_record, isolated_library):
        app, size = await self._app(tiny_record, isolated_library)
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            app.undo._undo_stack.clear()
            app._pop_undo_silently()            # must not raise
            assert app._undo_stack == []
            app.exit()
