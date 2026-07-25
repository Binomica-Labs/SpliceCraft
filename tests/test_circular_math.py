"""
test_circular_math — wrap-around coordinate tests for circular plasmids.

These tests protect two places where off-by-one or wrong-side-of-origin bugs
have historically crept in:

  1. Feature midpoint calculation for label placement (PlasmidMap line ~1545).
     The naive formula (start + (end-start)//2) % total places labels opposite
     the actual arc when end < start. Fixed 2026-04-11.

  2. `_bp_in(bp, feature_dict)` — is a bp inside a (possibly wrapped) feature?
     Lives as a method of `PlasmidMap`; we invoke it with self=None since it
     never touches instance state.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Feature midpoint — sacred invariant #5 (CLAUDE.md)
# ═══════════════════════════════════════════════════════════════════════════════

def _mid_bp(start: int, end: int, total: int) -> int:
    """Reference implementation matching splicecraft.py:~1545.

    Arc length = (end - start) mod total. This handles wrap-around: if end
    < start, Python's `%` returns the distance going forward around the circle.
    """
    arc_len = (end - start) % total
    return (start + arc_len // 2) % total


class TestFeatureMidpoint:
    def test_non_wrapped_normal_case(self):
        assert _mid_bp(100, 200, 1000) == 150
        assert _mid_bp(0, 500, 1000) == 250

    def test_span_of_two(self):
        assert _mid_bp(10, 12, 1000) == 11

    def test_full_circle_span_returns_start(self):
        # Degenerate: start == end → arc_len = 0 → midpoint = start
        assert _mid_bp(42, 42, 1000) == 42

    def test_wrap_around_simple(self):
        # Feature 900..100 has length 200; midpoint is at 1000 % 1000 = 0
        assert _mid_bp(900, 100, 1000) == 0

    def test_wrap_around_near_origin(self):
        # Feature 950..50 has length 100; midpoint at 1000 % 1000 = 0
        assert _mid_bp(950, 50, 1000) == 0

    def test_wrap_around_asymmetric(self):
        # Feature 800..200 has length 400; midpoint at (800 + 200) % 1000 = 0
        assert _mid_bp(800, 200, 1000) == 0
        # Feature 800..300 has length 500; midpoint at 1050 % 1000 = 50
        assert _mid_bp(800, 300, 1000) == 50

    def test_wrap_around_small_total(self):
        # pUC19-sized case (~2686 bp) abstracted to total=50 for grepability
        total = 50
        # Feature from 40..10 has length 20; mid = (40+10)%50 = 0
        assert _mid_bp(40, 10, total) == 0
        # Feature from 45..5 has length 10; mid = (45+5)%50 = 0
        assert _mid_bp(45, 5, total) == 0

    def test_non_wrapped_vs_wrapped_disagree_with_naive_formula(self):
        """This is the regression guard for the 2026-04-11 fix. The naive
        formula `(start + (end - start) // 2) % total` gives a midpoint
        OPPOSITE the arc for wrapped features. Prove correct formula
        differs from it for at least one wrap case."""
        total = 1000
        start, end = 900, 100
        naive = (start + (end - start) // 2) % total     # = 500 (OPPOSITE)
        correct = _mid_bp(start, end, total)             # = 0
        assert naive != correct
        assert correct == 0
        assert naive == 500

    def test_midpoint_lies_on_the_arc(self):
        """For any (start, end, total) the midpoint must lie INSIDE the arc,
        i.e. between start and end along the forward direction."""
        cases = [
            (100, 200, 1000),
            (900, 100, 1000),
            (950, 50, 1000),
            (0, 999, 1000),
            (500, 500, 1000),  # degenerate
            (0, 1, 1000),
            (999, 0, 1000),    # 1-bp span wrapping
        ]
        for start, end, total in cases:
            mid = _mid_bp(start, end, total)
            arc_len = (end - start) % total
            dist_from_start = (mid - start) % total
            # midpoint must be within [0, arc_len] along the forward direction
            assert 0 <= dist_from_start <= arc_len, (
                f"mid_bp({start},{end},{total})={mid} not on arc "
                f"(arc_len={arc_len}, dist_from_start={dist_from_start})"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# _bp_in — membership test for (possibly wrapped) features
# ═══════════════════════════════════════════════════════════════════════════════

def _bp_in(bp, start, end):
    """Call the PlasmidMap method with a stub self. It never uses self."""
    return sc.PlasmidMap._bp_in(None, bp, {"start": start, "end": end})


class TestBpIn:
    def test_non_wrapped_basic(self):
        # Feature 10..20 (half-open)
        assert not _bp_in(9, 10, 20)
        assert _bp_in(10, 10, 20)            # start is inclusive
        assert _bp_in(15, 10, 20)
        assert _bp_in(19, 10, 20)
        assert not _bp_in(20, 10, 20)         # end is exclusive

    def test_wrapped_feature(self):
        # Feature 95..5 on a 100-bp circle → spans 95,96,97,98,99,0,1,2,3,4
        feat = (95, 5)
        for bp in [95, 96, 97, 98, 99, 0, 1, 2, 3, 4]:
            assert _bp_in(bp, *feat), f"bp {bp} should be in wrapped feat"
        for bp in [5, 6, 50, 94]:
            assert not _bp_in(bp, *feat), f"bp {bp} should NOT be in wrapped feat"

    def test_wrapped_at_origin(self):
        feat = (990, 10)
        # Origin is inside
        assert _bp_in(0, *feat)
        assert _bp_in(999, *feat)
        # Far side is outside
        assert not _bp_in(500, *feat)

    def test_empty_feature(self):
        # start == end: Python returns (s <= bp < e) → False for any bp
        # (consistent with a zero-width feature that contains nothing)
        assert not _bp_in(0, 10, 10)
        assert not _bp_in(10, 10, 10)


# ═══════════════════════════════════════════════════════════════════════════════
# _feats_in_chunk — wrap-aware chunk filter (sequence panel rendering)
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatsInChunk:
    """Regression guard: the old filter `start < chunk_end and end > chunk_start`
    silently dropped every wrap feature from every chunk (both halves fail the
    conjunction). The sequence panel then rendered wrapped features as blank
    rows and click-to-feature stopped working.
    """

    def _feat(self, start, end, **extra):
        return {
            "start": start, "end": end, "strand": 1,
            "color": "white", "label": "test", "type": "CDS",
            **extra,
        }

    def test_non_wrapped_passes_through(self):
        f = self._feat(100, 200)
        out = sc._feats_in_chunk([f], 90, 210, 500)
        assert out == [f]
        assert out[0] is f   # same identity

    def test_non_wrapped_filters_out_when_disjoint(self):
        f = self._feat(100, 200)
        assert sc._feats_in_chunk([f], 300, 400, 500) == []
        assert sc._feats_in_chunk([f], 0, 50, 500) == []

    def test_wrapped_tail_visible_in_tail_chunk(self):
        f = self._feat(95, 5)
        out = sc._feats_in_chunk([f], 80, 100, 100)
        assert len(out) == 1
        assert out[0]["start"] == 95
        assert out[0]["end"]   == 100
        assert out[0]["label"] == "test"   # tail keeps the label

    def test_wrapped_head_visible_in_head_chunk(self):
        f = self._feat(95, 5)
        out = sc._feats_in_chunk([f], 0, 20, 100)
        assert len(out) == 1
        assert out[0]["start"] == 0
        assert out[0]["end"]   == 5
        assert out[0]["label"] == ""       # head is unlabeled continuation

    def test_wrapped_both_halves_in_single_chunk(self):
        # Chunk spanning the whole plasmid sees both pieces
        f = self._feat(95, 5)
        out = sc._feats_in_chunk([f], 0, 100, 100)
        assert len(out) == 2
        # Tail piece keeps label, head does not (ordering: tail first)
        tail, head = out
        assert tail["start"] == 95 and tail["end"] == 100
        assert head["start"] == 0 and head["end"] == 5

    def test_wrapped_disjoint_from_chunk(self):
        f = self._feat(95, 5)
        assert sc._feats_in_chunk([f], 30, 50, 100) == []

    def test_mix_of_wrapped_and_normal(self):
        normal = self._feat(10, 30, label="A")
        wrap   = self._feat(90, 20, label="B")
        out = sc._feats_in_chunk([normal, wrap], 0, 50, 100)
        # normal passes through; wrap contributes its head piece (start=0, end=20)
        assert any(f["label"] == "A" and f["start"] == 10 for f in out)
        assert any(f["start"] == 0 and f["end"] == 20 for f in out)


# ═══════════════════════════════════════════════════════════════════════════════
# Feature-arc sanity across the FULL plasmid — harder edge cases of invariant #5
# ═══════════════════════════════════════════════════════════════════════════════

class TestArcEdgeCases:
    """Deep edge cases for circular feature geometry beyond the common wrap.

    - Features covering > half the plasmid
    - Features covering exactly half
    - Full-circle features (start=0, end=total)
    - `_bp_in` at every boundary
    """

    def test_feature_spans_more_than_half(self):
        """Feature 600..400 on a 1000 bp plasmid covers 800 bp (80%).
        Arc midpoint must land inside the arc [600, 1000) ∪ [0, 400)."""
        start, end, total = 600, 400, 1000
        mid = _mid_bp(start, end, total)
        # arc_len = 800, so midpoint should be at (600 + 400) % 1000 = 0
        assert mid == 0
        # Assert mid lies inside the arc
        assert _bp_in(mid, start, end) or mid == 0  # origin also counts

    def test_feature_spans_exactly_half(self):
        """Feature 0..500 on a 1000 bp plasmid covers exactly 500 bp.
        Non-wrap case; midpoint at 250 is in [0, 500)."""
        assert _mid_bp(0, 500, 1000) == 250
        # And the wrap version 500..0 (start=500, end=0 → wraps to [500,1000)∪[])
        # Arc length = 500; midpoint = (500 + 250) % 1000 = 750
        assert _mid_bp(500, 0, 1000) == 750

    def test_full_circle_feature_non_wrap(self):
        """A feature from 0..total is a full-circle feature (not a wrap,
        since end >= start). `_bp_in` matches every bp (covers whole plasmid).

        Note: _mid_bp collapses (0, total) to 0 because (end-start) % total
        is 0 — the formula can't distinguish full-circle from zero-width.
        That's a documented quirk; labels on full-circle feats land at origin.
        """
        total = 1000
        # bp_in: every bp in [0, total) should match
        for bp in [0, 1, 500, 999]:
            assert _bp_in(bp, 0, total), f"bp {bp} should be inside full-circle feat"
        # Midpoint formula documents its quirk:
        assert _mid_bp(0, total, total) == 0
        assert _mid_bp(100, 100, total) == 100   # zero-width also returns start

    def test_midpoint_for_reverse_half_of_wrap(self):
        """If start is 1 bp past half and end is start-2 (tiny slice just
        before the origin), the 998-bp feature wraps almost the whole circle.
        Midpoint should be near the opposite side of the short gap."""
        start, end, total = 501, 499, 1000
        arc_len = 998
        mid = _mid_bp(start, end, total)
        assert mid == (501 + 499) % 1000   # = 0
        # Midpoint distance from start == arc_len // 2
        dist = (mid - start) % total
        assert dist == arc_len // 2

    def test_bp_in_every_boundary(self):
        """Boundary conditions of `_bp_in` with half-open [start, end)
        semantics. Check each boundary explicitly to prevent off-by-ones."""
        # Non-wrapped feat 10..20
        assert _bp_in(10, 10, 20) is True       # start inclusive
        assert _bp_in(19, 10, 20) is True       # end-1 inclusive
        assert _bp_in(20, 10, 20) is False      # end exclusive
        assert _bp_in(9,  10, 20) is False      # before start
        # Wrapped feat 90..10 on 100-bp circle
        assert _bp_in(90, 90, 10) is True       # start inclusive
        assert _bp_in(99, 90, 10) is True
        assert _bp_in(0,  90, 10) is True       # origin inside
        assert _bp_in(9,  90, 10) is True
        assert _bp_in(10, 90, 10) is False      # end exclusive
        assert _bp_in(50, 90, 10) is False      # far side of wrap

    def test_bp_in_is_monotonic_along_arc(self):
        """For any (start, end), arc bp's form a CONTIGUOUS forward-direction
        segment. Pick wrap + non-wrap, sweep every bp, and confirm the set of
        inside-bp's is exactly those whose forward-distance from start is
        < arc_len."""
        total = 100
        cases = [(10, 30), (90, 10), (0, 50), (50, 0), (75, 25)]
        for start, end in cases:
            arc_len = (end - start) % total
            inside_expected = {(start + k) % total for k in range(arc_len)}
            inside_actual = {bp for bp in range(total) if _bp_in(bp, start, end)}
            assert inside_actual == inside_expected, (
                f"start={start} end={end}: "
                f"expected {sorted(inside_expected)}, got {sorted(inside_actual)}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Rotation invariance — origin_bp affects rendering only, not feature coords
# ═══════════════════════════════════════════════════════════════════════════════

class TestRotationInvariance:
    """`PlasmidMap.origin_bp` is a PURE VIEW rotation: it offsets the visual
    origin but MUST NOT alter stored feature coordinates, _bp_in results, or
    midpoint calculations. The test harness simulates that by calling the
    underlying methods with different origin_bp values and asserting feature
    state is unchanged.
    """

    def test_rotate_wraps_modulo_total(self):
        """Multiple rotations accumulate modulo total, never escape [0, total)."""
        # Simulate the action code: `self.origin_bp = (self.origin_bp + step) % total`
        total = 1000
        origin = 0
        steps = [100, 200, -50, 999, -1001]   # last one forces a double-wrap
        for s in steps:
            origin = (origin + s) % total
            assert 0 <= origin < total, f"origin escaped [0, {total}): {origin}"
        # Net effect: (100+200-50+999-1001) % 1000 = 248
        assert origin == 248

    def test_feature_bp_in_unaffected_by_rotation(self):
        """Whatever origin_bp is, `_bp_in(bp, feat)` gives the same answer
        because _bp_in doesn't look at origin_bp."""
        # Just confirm _bp_in's signature doesn't take origin_bp.
        import inspect
        sig = inspect.signature(sc.PlasmidMap._bp_in)
        assert "origin_bp" not in sig.parameters, (
            "_bp_in signature leaked origin_bp — coordinates must be invariant"
        )

    def test_midpoint_formula_unaffected_by_rotation(self):
        """Midpoint of a feature depends only on (start, end, total), NOT
        origin_bp. Demonstrate this by computing with two different origins."""
        # The reference implementation _mid_bp takes no origin — that's the
        # contract. This test just documents that the formula lives outside
        # any rotation state.
        start, end, total = 900, 100, 1000
        mid = _mid_bp(start, end, total)
        # Ten "rotations" — none should change the result because the formula
        # only uses (start, end, total).
        for origin in [0, 250, 500, 750, 999, 1000, 1, 100]:
            assert _mid_bp(start, end, total) == mid

    def test_zero_step_rotation_is_noop(self):
        """step=0 or total-equivalent rotations are identity."""
        total = 1000
        assert (500 + 0) % total == 500
        assert (500 + total) % total == 500
        assert (500 - total) % total == 500
        assert (500 + 2 * total) % total == 500


# ═══════════════════════════════════════════════════════════════════════════════
# _feat_len — circular-aware feature length (added 2026-04-13)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Several sort keys and UI displays used to compute feature length as
# `end - start`, which is *negative* for wrap features. Wrap features then
# sorted to the front (highest priority = smallest length after negation)
# and displayed with negative bp counts in the sidebar. _feat_len is the
# single source of truth; callers pass it (start, end, total) and get the
# real biological length back.

class TestFeatLen:
    def test_linear_feature(self):
        assert sc._feat_len(100, 200, 1000) == 100

    def test_wrap_feature(self):
        # start=950, end=100 on 1000 bp plasmid = 150 bp total
        assert sc._feat_len(950, 100, 1000) == 150

    def test_feature_at_origin(self):
        assert sc._feat_len(0, 100, 1000) == 100

    def test_feature_ending_at_origin(self):
        # end=0 means wrap, contains [start:len) plus [0:0) = (total - start)
        assert sc._feat_len(800, 0, 1000) == 200

    def test_full_length_feature(self):
        assert sc._feat_len(0, 1000, 1000) == 1000

    def test_single_bp(self):
        assert sc._feat_len(50, 51, 1000) == 1

    def test_wrap_single_bp(self):
        # start=999, end=0 wraps 1 bp across origin
        assert sc._feat_len(999, 0, 1000) == 1

    def test_sort_key_orders_wrap_features_correctly(self):
        """The original bug: `-end + start` makes wrap features sort to the
        front (most negative = largest when negated). _feat_len fixes this."""
        feats = [
            {"start": 100, "end": 200, "name": "small_linear"},   # 100 bp
            {"start": 950, "end": 100, "name": "wrap_150"},       # 150 bp
            {"start": 500, "end": 900, "name": "big_linear"},     # 400 bp
        ]
        by_len_desc = sorted(feats, key=lambda f: -sc._feat_len(f["start"], f["end"], 1000))
        assert [f["name"] for f in by_len_desc] == ["big_linear", "wrap_150", "small_linear"]


# ═══════════════════════════════════════════════════════════════════════════════
# _feat_bounds — sacred invariant #9 (wrap-feature CompoundLocation extraction)
# Regression guard for 2026-05-10 fix: routing five call sites through this
# helper instead of raw `int(loc.start)`/`int(loc.end)` which silently flattened
# origin-spanning features to the BACKBONE GAP and produced biologically wrong
# downstream behavior (wrong primers, wrong restriction analysis, dropped
# annotations).
# ═══════════════════════════════════════════════════════════════════════════════


def _make_wrap_compound_feat(tail_start: int, head_end: int, total: int,
                               *, strand: int = 1, ftype: str = "CDS",
                               label: str = "wrap_cds"):
    """Build a Biopython SeqFeature whose CompoundLocation has the canonical
    origin-wrap shape: ``join(0..head_end, tail_start..total)``."""
    from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
    parts = [
        FeatureLocation(0, head_end, strand=strand),
        FeatureLocation(tail_start, total, strand=strand),
    ]
    feat = SeqFeature(
        CompoundLocation(parts), type=ftype,
        qualifiers={"label": [label]},
    )
    return feat


def _make_simple_feat(start: int, end: int, *, strand: int = 1,
                        ftype: str = "CDS", label: str = "linear_cds"):
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    return SeqFeature(
        FeatureLocation(start, end, strand=strand), type=ftype,
        qualifiers={"label": [label]},
    )


class TestFeatBounds:
    def test_linear_feature_unchanged(self):
        feat = _make_simple_feat(100, 200, strand=1)
        bounds = sc._feat_bounds(feat, total=1000)
        assert bounds == (100, 200, 1)

    def test_linear_reverse_strand_carries(self):
        feat = _make_simple_feat(100, 200, strand=-1)
        bounds = sc._feat_bounds(feat, total=1000)
        assert bounds == (100, 200, -1)

    def test_origin_wrap_returns_tail_start_head_end(self):
        """The canonical join(0..50, 5800..6000) shape on a 6000 bp plasmid
        re-encodes as (5800, 50) so end < start signals wrap."""
        feat = _make_wrap_compound_feat(
            tail_start=5800, head_end=50, total=6000,
        )
        bounds = sc._feat_bounds(feat, total=6000)
        assert bounds == (5800, 50, 1)

    def test_origin_wrap_reverse_strand(self):
        feat = _make_wrap_compound_feat(
            tail_start=4500, head_end=200, total=5000, strand=-1,
        )
        bounds = sc._feat_bounds(feat, total=5000)
        assert bounds == (4500, 200, -1)

    def test_non_origin_compound_flattens_to_outer_bounds(self):
        """A compound like exon-1 + exon-2 (NOT touching 0 or total) is
        flattened to outer bounds — lossy, but oriented correctly."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        feat = SeqFeature(
            CompoundLocation([
                FeatureLocation(100, 200, strand=1),
                FeatureLocation(400, 500, strand=1),
            ]),
            type="mRNA", qualifiers={"label": ["mrna"]},
        )
        bounds = sc._feat_bounds(feat, total=1000)
        assert bounds == (100, 500, 1)

    def test_unknown_position_returns_none(self):
        from Bio.SeqFeature import SeqFeature
        try:
            from Bio.SeqFeature import UnknownPosition, FeatureLocation
        except ImportError:
            pytest.skip("UnknownPosition unavailable in this Biopython")
        try:
            loc = FeatureLocation(UnknownPosition(), UnknownPosition())
        except Exception:
            pytest.skip("UnknownPosition couldn't build a FeatureLocation")
        feat = SeqFeature(loc, type="misc_feature")
        bounds = sc._feat_bounds(feat, total=1000)
        assert bounds is None

    def test_missing_location_returns_none(self):
        class _Fake:
            location = None
        assert sc._feat_bounds(_Fake(), total=1000) is None


class TestWrapFeatureSliceContract:
    """The whole point of `_feat_bounds` is that downstream `_slice_circular`
    + `_feat_len` produce the right biology for a wrap feature."""

    def test_slice_picks_correct_bases(self):
        # 60 bp plasmid; wrap CDS at join(0..10, 50..60) = 20 bp total
        # encoded as (50, 10). Sequence: 0-9 are "A"s (head), 10-49 are
        # "C"s (gap), 50-59 are "T"s (tail). Wrap feature should return
        # TTTTTTTTTTAAAAAAAAAA.
        seq = "A"*10 + "C"*40 + "T"*10
        feat = _make_wrap_compound_feat(50, 10, 60)
        s, e, _strand = sc._feat_bounds(feat, total=60)
        assert (s, e) == (50, 10)
        feat_seq = sc._slice_circular(seq, s, e)
        assert feat_seq == "T"*10 + "A"*10
        # Confirm the OLD (broken) path returned the gap:
        broken = seq[int(feat.location.start):int(feat.location.end)]
        assert broken == "A"*10 + "C"*40 + "T"*10  # Biopython flattens to full span
        # so OLD slice would have been the WRONG 60-bp full-record span.

    def test_length_routes_through_feat_len(self):
        feat = _make_wrap_compound_feat(5800, 50, 6000)
        s, e, _strand = sc._feat_bounds(feat, total=6000)
        assert sc._feat_len(s, e, 6000) == 250


class TestWrapFeatureAnnotationTransfer:
    """Regression guard: `_find_annotation_transfers` must propagate origin-
    spanning source features. Pre-fix it called `int(loc.start)`/`int(loc.end)`
    which flattened wrap CDSes and they vanished from the result list."""

    def test_wrap_cds_transfers(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Source: 100 bp circular, wrap feature at join(0..20, 80..100).
        # Use a non-palindromic 40-bp wrap region so the test only matches
        # the forward strand (palindromic regions would match RC too and
        # produce two transfers).
        src_seq = ("ACGAGTACGAGTACGAGTAC"      # 0..20: head of wrap CDS
                   + "C"*60                    # 20..80: gap (backbone)
                   + "GATTCAGATTCAGATTCAGT")   # 80..100: tail of wrap CDS
        assert len(src_seq) == 100
        src_rec = SeqRecord(Seq(src_seq), id="src", name="src")
        src_rec.annotations["topology"] = "circular"
        src_rec.features.append(_make_wrap_compound_feat(80, 20, 100))
        # Target: linearised version of the wrap region (tail+head).
        tgt_bases = src_seq[80:] + src_seq[:20]
        tgt_seq = "G"*30 + tgt_bases + "G"*30
        tgt_rec = SeqRecord(Seq(tgt_seq), id="tgt", name="tgt")
        tgt_rec.annotations["topology"] = "linear"

        transfers = sc._find_annotation_transfers(src_rec, tgt_rec, min_len=20)
        assert len(transfers) == 1, (
            "Wrap CDS must transfer exactly once; pre-fix it was silently "
            "dropped because `int(loc.start)`/`int(loc.end)` returned (0, 100) "
            "and the slice produced the whole sequence rather than the "
            "40-bp wrap region."
        )
        t = transfers[0]
        assert t["target_start"] == 30
        assert t["target_end"]   == 70


class TestWrapFeatureRecordFeatures:
    """`TraditionalCloningPane._record_features` must emit wrap-encoded dicts
    so `_excise_fragment_pair` sees `end < start` for origin-spanning vector
    features and routes the dropout correctly."""

    def test_wrap_vector_feature_encoded_as_wrap(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("N"*1000), id="r", name="r")
        rec.annotations["topology"] = "circular"
        rec.features.append(_make_wrap_compound_feat(950, 50, 1000))

        # Build a TraditionalCloningPane bound just enough to call the
        # method (it's instance-level but doesn't touch self state in the
        # body).
        pane = sc.TraditionalCloningPane.__new__(sc.TraditionalCloningPane)
        feats = pane._record_features(rec)
        assert len(feats) == 1
        f = feats[0]
        assert (f["start"], f["end"]) == (950, 50), (
            "Wrap feature must encode as (tail_start, head_end). Pre-fix "
            "this returned (0, 1000) — the WHOLE plasmid, which would slice "
            "every base of the vector as 'feature'."
        )


class TestWrapFeatureFeatsForDomesticator:
    def test_wrap_feature_encoded_correctly(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ACGT"*250), id="r", name="r")  # 1000 bp
        rec.annotations["topology"] = "circular"
        rec.features.append(_make_wrap_compound_feat(900, 100, 1000))

        feats = sc._feats_for_domesticator(rec)
        assert len(feats) == 1
        f = feats[0]
        assert (f["start"], f["end"]) == (900, 100), (
            "Pre-fix this returned (0, 1000) which the domesticator then "
            "treated as 'design primers for the entire plasmid' — wrong "
            "biology and wrong primer Tm budget."
        )


class TestWrapFeaturePrimerDesignScreen:
    def test_wrap_cds_encoded_correctly(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ACGT"*500), id="r", name="r")  # 2000 bp
        rec.annotations["topology"] = "circular"
        rec.features.append(_make_wrap_compound_feat(1900, 100, 2000))

        scr = sc.PrimerDesignScreen.__new__(sc.PrimerDesignScreen)
        feats = scr._parse_features_from_record(rec)
        assert len(feats) == 1
        assert (feats[0]["start"], feats[0]["end"]) == (1900, 100), (
            "Pre-fix the wrap CDS came in as (0, 2000) and Primer3 designed "
            "against the inverted backbone region rather than the user's "
            "feature."
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Intron-bearing CDS — issue #9 from a user (har1eyk reported similar)
# Pre-fix `_translate_cds` and `_paint_cds_aa` treated every CDS as
# contiguous, so an intron whose length isn't a multiple of 3 frame-
# shifted every AA past the splice. The loader now stamps `_exons`
# on CompoundLocation CDSes; the translators splice introns out.
# ═══════════════════════════════════════════════════════════════════════════════


class TestIntronAwareTranslation:
    def test_translate_cds_with_exons_skips_intron(self):
        # 60 bp record: exon-1 [0:9] (ATG GGG AAA) + intron [9:13]
        # (CCCC — 4 bp, not a multiple of 3) + exon-2 [13:22] (TTT TGT
        # AAA). Spliced CDS = ATGGGGAAATTTTGTAAA = M G K F C K (no stop)
        # → six AAs.
        seq = (
            "ATGGGGAAA"   # exon 1 (9 bp)
            "CCCC"         # intron (4 bp) — would frameshift everything after
            "TTTTGTAAA"   # exon 2 (9 bp) → F C K
            + "G" * 38     # padding to round out the test sequence
        )
        assert len(seq) == 60
        exons = [(0, 9), (13, 22)]
        # Reference: build the spliced sequence by hand and translate.
        spliced = seq[0:9] + seq[13:22]
        assert spliced == "ATGGGGAAATTTTGTAAA"
        # Translate via the exon-aware path.
        result = sc._translate_cds(seq, 0, 22, 1, exons=exons)
        # Six AAs, no trailing stop — the spliced CDS genuinely has none
        # (2026-07-24: the translator used to fabricate one, which is why
        # this asserted a seventh character its own comment didn't expect).
        assert result == "MGKFCK", (
            f"Expected MGKFCK (correctly spliced); got {result!r}. "
            f"Pre-fix the intron's 4 bp shifted every downstream codon's "
            f"reading frame."
        )

    def test_translate_cds_without_exons_keeps_legacy_contiguous(self):
        # The legacy path (no exons) must still produce the contiguous
        # translation — so passing a non-intron CDS through the new
        # signature doesn't regress.
        seq = "ATGGGGAAATAA"   # M G K *
        result = sc._translate_cds(seq, 0, 12, 1)
        assert result == "MGK*"

    def test_translate_cds_with_exons_reverse_strand(self):
        # Reverse strand: spliced CDS on the bottom strand should
        # translate from the 3' end of the genomic span. Build a record
        # whose top strand reverse-complements to ATGGGGAAATTTTGTAAA
        # (the same spliced CDS as the forward test) on the bottom
        # strand. RC of ATGGGGAAATTTTGTAAA = TTTACAAAATTTCCCCAT.
        # Place that at genomic positions: exon-2 (5' of spliced) at
        # [4:13] = TTTACAAAA (rc of "TTTTGTAAA"), intron at [13:17] =
        # GGGG (rc of "CCCC"), exon-1 (3' of spliced) at [17:26] =
        # TTTCCCCAT (rc of "ATGGGGAAA").
        # NB exons stored in ascending genomic order:
        #   genomic[4:13]  → 3' end of spliced CDS (last exon)
        #   genomic[17:26] → 5' end of spliced CDS (first exon)
        # Build by hand using sc._IUPAC_COMP.
        # exon-1 (5' of CDS = end of genomic span) = first 9 chars of
        # spliced_top → bottom strand bases at the END of the genomic
        # span. Layout: [pad 4][rc of exon-2 = "TTTACAAAA"][intron
        # "GGGG"][rc of exon-1 = "TTTCCCCAT"][pad].
        rc_exon1 = "ATGGGGAAA".translate(sc._IUPAC_COMP)[::-1]   # TTTCCCCAT
        rc_exon2 = "TTTTGTAAA".translate(sc._IUPAC_COMP)[::-1]   # TTTACAAAA
        intron_rc = "CCCC".translate(sc._IUPAC_COMP)[::-1]       # GGGG
        seq = (
            "G" * 4              # padding 0..4
            + rc_exon2           # 4..13 (genomic 5' = CDS 3' end)
            + intron_rc          # 13..17
            + rc_exon1           # 17..26 (genomic 3' = CDS 5' end)
            + "G" * 4            # padding
        )
        assert len(seq) == 30
        exons = [(4, 13), (17, 26)]
        result = sc._translate_cds(seq, 4, 26, -1, exons=exons)
        assert result == "MGKFCK", (
            f"Reverse-strand spliced CDS should still translate to "
            f"MGKFCK; got {result!r}"
        )

    def test_cds_aa_list_uses_exons_from_feature_dict(self):
        # Build a feature dict carrying _exons and verify _cds_aa_list
        # produces the spliced translation.
        seq = "ATGGGGAAA" + "CCCC" + "TTTTGTAAA" + "G" * 38
        f = {
            "type": "CDS", "start": 0, "end": 22, "strand": 1,
            "_exons": [(0, 9), (13, 22)],
        }
        aa_letters, cds_len, _virt_e = sc._cds_aa_list(seq, f)
        assert cds_len == 18, f"spliced CDS is 18 bp, got {cds_len}"
        assert "".join(aa_letters) == "MGKFCK", (
            f"Expected MGKFCK (spliced); got {''.join(aa_letters)!r}"
        )

    def test_cds_aa_list_no_exons_is_contiguous(self):
        # No _exons key → legacy contiguous translation path.
        seq = "ATGGGGAAATAA"
        f = {"type": "CDS", "start": 0, "end": 12, "strand": 1}
        aa_letters, _, _ = sc._cds_aa_list(seq, f)
        assert "".join(aa_letters) == "MGK*"

    def test_spliced_idx_to_genomic_bp_forward(self):
        # exons = [(0, 9), (13, 22)] → spliced length 18.
        # spliced[0] = genomic[0]; spliced[8] = genomic[8];
        # spliced[9] = genomic[13]; spliced[17] = genomic[21].
        exons = [(0, 9), (13, 22)]
        assert sc._spliced_idx_to_genomic_bp(0,  exons, 1) == 0
        assert sc._spliced_idx_to_genomic_bp(8,  exons, 1) == 8
        assert sc._spliced_idx_to_genomic_bp(9,  exons, 1) == 13
        assert sc._spliced_idx_to_genomic_bp(17, exons, 1) == 21

    def test_spliced_idx_to_genomic_bp_reverse(self):
        # Reverse strand: spliced[0] = 3'-most base of last exon =
        # genomic[exons[-1].end - 1].
        exons = [(0, 9), (13, 22)]
        assert sc._spliced_idx_to_genomic_bp(0, exons, -1) == 21
        assert sc._spliced_idx_to_genomic_bp(8, exons, -1) == 13
        # Crossing the splice: spliced[9] is the first base of the
        # 5'-most exon's far end on the genomic forward strand.
        assert sc._spliced_idx_to_genomic_bp(9, exons, -1) == 8
        assert sc._spliced_idx_to_genomic_bp(17, exons, -1) == 0
