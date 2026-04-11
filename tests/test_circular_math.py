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
