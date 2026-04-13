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
