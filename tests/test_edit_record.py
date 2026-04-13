"""Tests for PlasmidApp._rebuild_record_with_edit — the function that shifts
feature coordinates after a user insert/replace so the SeqRecord stays in
sync with the edited sequence.

Two pre-fix bugs under regression-guard here (both found 2026-04-13):

  Bug A (wrap mangling) — `int(CompoundLocation.start)` returns min(parts),
  `int(.end)` returns max(parts), so reading `fs = int(feat.location.start)`
  on a wrap feature stored as `CompoundLocation([FL(tail..n), FL(0..head)])`
  gave `fs=0, fe=n`, i.e. the whole plasmid. The function then rebuilt
  every edit-affected wrap feature as a FeatureLocation covering the
  entire plasmid, silently mangling the wrap semantics.

  Bug B (1-bp ghost) — the clamp `new_fe = max(new_fs + 1, min(new_fe, new_len))`
  forced every post-edit feature to be at least 1 bp wide, so features
  fully consumed by a replace/delete survived as 1-bp stubs instead of
  being dropped.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

import splicecraft as sc


def _make_app(record):
    """Bare PlasmidApp shell that _rebuild_record_with_edit can run on
    without mounting the Textual runtime."""
    app = sc.PlasmidApp.__new__(sc.PlasmidApp)
    app._current_record = record
    return app


def _wrap_record(total: int = 100, tail_start: int = 95, head_end: int = 5,
                 strand: int = 1, label: str = "wrapCDS",
                 extra_feats=None) -> SeqRecord:
    """Build a circular plasmid with a single canonical-form wrap feature
    `join(tail_start+1..total, 1..head_end)`."""
    rec = SeqRecord(Seq("A" * total), id="T", annotations={"molecule_type": "DNA"})
    wrap = CompoundLocation([
        FeatureLocation(tail_start, total, strand=strand),
        FeatureLocation(0, head_end, strand=strand),
    ])
    rec.features.append(SeqFeature(wrap, type="CDS", qualifiers={"label": [label]}))
    for (s, e, typ, lab) in (extra_feats or []):
        rec.features.append(SeqFeature(
            FeatureLocation(s, e, strand=1),
            type=typ, qualifiers={"label": [lab]},
        ))
    return rec


def _first_by_label(record, label: str):
    for f in record.features:
        if f.qualifiers.get("label", [None])[0] == label:
            return f
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Wrap-feature preservation across edits (regression guard for bug A)
# ═══════════════════════════════════════════════════════════════════════════════

class TestWrapFeaturePreserved:
    """A wrap feature stored as CompoundLocation must survive an edit
    as a CompoundLocation — not get flattened into a whole-plasmid span."""

    def test_insert_far_from_wrap_preserves_structure(self):
        """Insert in the middle, nowhere near the wrap. Tail shifts by
        ins_len, head unchanged. Round-trips through _parse as wrap."""
        rec = _wrap_record()
        app = _make_app(rec)
        new_seq = "A" * 50 + "T" * 10 + "A" * 50  # insert 10 bp at pos 50
        new_rec = app._rebuild_record_with_edit(new_seq, "insert", 50, 50, "T" * 10)

        wrap = _first_by_label(new_rec, "wrapCDS")
        assert isinstance(wrap.location, CompoundLocation)
        parts = sorted([(int(p.start), int(p.end)) for p in wrap.location.parts])
        assert parts == [(0, 5), (105, 110)]

        # _parse should re-detect as a wrap with start=105, end=5
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(new_rec)
        wrap_feat = next(f for f in feats if f.get("label") == "wrapCDS")
        assert wrap_feat["start"] == 105 and wrap_feat["end"] == 5
        assert pm._n_flattened == 0

    def test_insert_inside_tail_stretches_tail(self):
        """Insert at position 97 falls inside the tail (95, 100). Tail
        spans the insert and stretches to (95, 100+ins_len); head (0, 5)
        unchanged."""
        rec = _wrap_record()
        app = _make_app(rec)
        new_seq = "A" * 97 + "T" * 10 + "A" * 3
        new_rec = app._rebuild_record_with_edit(new_seq, "insert", 97, 97, "T" * 10)

        wrap = _first_by_label(new_rec, "wrapCDS")
        assert isinstance(wrap.location, CompoundLocation)
        parts = sorted([(int(p.start), int(p.end)) for p in wrap.location.parts])
        assert parts == [(0, 5), (95, 110)]
        # new_len is 110, so parts[-1].end == new_len → still a wrap
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(new_rec)
        wrap_feat = next(f for f in feats if f.get("label") == "wrapCDS")
        assert wrap_feat["start"] == 95 and wrap_feat["end"] == 5

    def test_insert_inside_head_stretches_head(self):
        """Insert at position 3 falls inside the head (0, 5). Head spans
        insert → (0, 5+ins_len); tail shifts fully → (95+ins_len, 100+ins_len)."""
        rec = _wrap_record()
        app = _make_app(rec)
        new_seq = "A" * 3 + "T" * 10 + "A" * 97
        new_rec = app._rebuild_record_with_edit(new_seq, "insert", 3, 3, "T" * 10)

        wrap = _first_by_label(new_rec, "wrapCDS")
        assert isinstance(wrap.location, CompoundLocation)
        parts = sorted([(int(p.start), int(p.end)) for p in wrap.location.parts])
        assert parts == [(0, 15), (105, 110)]
        pm = sc.PlasmidMap.__new__(sc.PlasmidMap)
        feats = pm._parse(new_rec)
        wrap_feat = next(f for f in feats if f.get("label") == "wrapCDS")
        assert wrap_feat["start"] == 105 and wrap_feat["end"] == 15

    def test_replace_consumes_entire_tail_collapses_to_head(self):
        """Replacing [95, 100) with '' deletes the entire tail part. The
        remaining single part (head 0..5) should collapse from
        CompoundLocation to FeatureLocation."""
        rec = _wrap_record()
        app = _make_app(rec)
        new_seq = "A" * 95  # 95 bp now
        new_rec = app._rebuild_record_with_edit(new_seq, "replace", 95, 100, "")

        wrap = _first_by_label(new_rec, "wrapCDS")
        # Only 1 part survives — should no longer be a CompoundLocation
        assert not isinstance(wrap.location, CompoundLocation)
        assert int(wrap.location.start) == 0 and int(wrap.location.end) == 5

    def test_replace_consumes_entire_head_collapses_to_tail(self):
        """Inverse of above — delete [0, 5), head is consumed, only the
        tail remains."""
        rec = _wrap_record()
        app = _make_app(rec)
        new_seq = "A" * 95  # 95 bp now
        new_rec = app._rebuild_record_with_edit(new_seq, "replace", 0, 5, "")

        wrap = _first_by_label(new_rec, "wrapCDS")
        assert not isinstance(wrap.location, CompoundLocation)
        # Tail (95, 100) shifted by delta=-5 → (90, 95)
        assert int(wrap.location.start) == 90 and int(wrap.location.end) == 95

    def test_replace_consumes_entire_wrap_drops_feature(self):
        """Replace that swallows BOTH head and tail — realistically only
        possible via two separate edits, but the per-part shift should
        still drop the whole feature if every part is consumed."""
        rec = _wrap_record(total=100, tail_start=95, head_end=5,
                           extra_feats=[(20, 40, "gene", "survivor")])
        app = _make_app(rec)
        # Replace [0, 5) with '' first
        mid = app._rebuild_record_with_edit("A" * 95, "replace", 0, 5, "")
        # Then replace [90, 95) (new tail position) with ''
        app._current_record = mid
        end = app._rebuild_record_with_edit("A" * 90, "replace", 90, 95, "")

        # Wrap feature should be gone; survivor shifted by first delete
        # (replace [0, 5) with '' shifts everything after by -5 → (15, 35))
        assert _first_by_label(end, "wrapCDS") is None
        survivor = _first_by_label(end, "survivor")
        assert survivor is not None
        assert int(survivor.location.start) == 15 and int(survivor.location.end) == 35

    def test_reverse_strand_wrap_preserves_strand(self):
        """A wrap feature on the reverse strand must survive an edit with
        both parts still on strand=-1."""
        rec = _wrap_record(strand=-1)
        app = _make_app(rec)
        new_seq = "A" * 50 + "T" * 10 + "A" * 50
        new_rec = app._rebuild_record_with_edit(new_seq, "insert", 50, 50, "T" * 10)

        wrap = _first_by_label(new_rec, "wrapCDS")
        assert isinstance(wrap.location, CompoundLocation)
        for part in wrap.location.parts:
            assert part.strand == -1


# ═══════════════════════════════════════════════════════════════════════════════
# Ghost-drop (regression guard for bug B)
# ═══════════════════════════════════════════════════════════════════════════════

class TestFullyConsumedFeatureDropped:
    """Features fully contained in a replace region must be DROPPED, not
    left behind as 1-bp stubs."""

    def test_delete_whole_feature_drops_it(self):
        """Feature at (60, 65), replace [60, 65) with '' — feature is gone."""
        rec = _wrap_record(extra_feats=[(60, 65, "promoter", "doomed"),
                                        (20, 40, "gene", "survivor")])
        app = _make_app(rec)
        new_seq = "A" * 60 + "A" * 35  # 95 bp; bases 60..64 removed
        new_rec = app._rebuild_record_with_edit(new_seq, "replace", 60, 65, "")

        assert _first_by_label(new_rec, "doomed") is None
        survivor = _first_by_label(new_rec, "survivor")
        assert survivor is not None
        assert int(survivor.location.start) == 20 and int(survivor.location.end) == 40
        # Other non-wrap features — count check
        non_wrap = [f for f in new_rec.features
                    if not isinstance(f.location, CompoundLocation)]
        assert len(non_wrap) == 1  # just "survivor"

    def test_delete_larger_region_drops_contained_features(self):
        """Replace [50, 80) with '' — any feature entirely inside is dropped."""
        rec = _wrap_record(extra_feats=[(55, 60, "promoter", "doomed1"),
                                        (65, 75, "gene", "doomed2"),
                                        (10, 30, "gene", "survivor")])
        app = _make_app(rec)
        new_seq = "A" * 50 + "A" * 20  # 70 bp; bases 50..79 removed
        new_rec = app._rebuild_record_with_edit(new_seq, "replace", 50, 80, "")

        assert _first_by_label(new_rec, "doomed1") is None
        assert _first_by_label(new_rec, "doomed2") is None
        assert _first_by_label(new_rec, "survivor") is not None

    def test_replace_that_shrinks_feature_to_zero_drops_it(self):
        """Feature at (60, 65), replace [55, 70) with '' — also drops."""
        rec = _wrap_record(extra_feats=[(60, 65, "promoter", "doomed")])
        app = _make_app(rec)
        new_seq = "A" * 55 + "A" * 30  # 85 bp; bases 55..69 removed
        new_rec = app._rebuild_record_with_edit(new_seq, "replace", 55, 70, "")

        assert _first_by_label(new_rec, "doomed") is None

    def test_replace_with_nonempty_keeps_feature(self):
        """Feature at (60, 65), replace [60, 65) with 'TT' (2 bp) — the
        feature survives as (60, 62) since the region was replaced, not
        deleted."""
        rec = _wrap_record(extra_feats=[(60, 65, "promoter", "tiny")])
        app = _make_app(rec)
        new_seq = "A" * 60 + "TT" + "A" * 35
        new_rec = app._rebuild_record_with_edit(new_seq, "replace", 60, 65, "TT")

        tiny = _first_by_label(new_rec, "tiny")
        assert tiny is not None
        # fs=60 <= s=60, fe=65 >= e=65 → "spans replaced region" branch
        # → (fs, fe+delta) = (60, 65 + (2-5)) = (60, 62)
        assert int(tiny.location.start) == 60 and int(tiny.location.end) == 62


# ═══════════════════════════════════════════════════════════════════════════════
# Simple features (pre-existing behavior preserved)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSimpleFeatureShift:
    """The pre-existing simple-feature shift logic must still work — we
    only changed consumed-feature handling and wrap awareness."""

    def test_feature_before_insert_unchanged(self):
        rec = _wrap_record(extra_feats=[(10, 20, "gene", "before")])
        app = _make_app(rec)
        new_seq = "A" * 50 + "T" * 10 + "A" * 50
        new_rec = app._rebuild_record_with_edit(new_seq, "insert", 50, 50, "T" * 10)
        f = _first_by_label(new_rec, "before")
        assert int(f.location.start) == 10 and int(f.location.end) == 20

    def test_feature_after_insert_shifts(self):
        rec = _wrap_record(extra_feats=[(70, 80, "gene", "after")])
        app = _make_app(rec)
        new_seq = "A" * 50 + "T" * 10 + "A" * 50
        new_rec = app._rebuild_record_with_edit(new_seq, "insert", 50, 50, "T" * 10)
        f = _first_by_label(new_rec, "after")
        assert int(f.location.start) == 80 and int(f.location.end) == 90

    def test_feature_spans_insert_stretches(self):
        rec = _wrap_record(extra_feats=[(40, 70, "gene", "spans")])
        app = _make_app(rec)
        new_seq = "A" * 50 + "T" * 10 + "A" * 50
        new_rec = app._rebuild_record_with_edit(new_seq, "insert", 50, 50, "T" * 10)
        f = _first_by_label(new_rec, "spans")
        assert int(f.location.start) == 40 and int(f.location.end) == 80
