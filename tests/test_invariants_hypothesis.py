"""
test_invariants_hypothesis — property-based fuzzing of the sacred invariants.

The hand-written tests in `test_dna_sanity.py` and friends cover known edge
cases. These tests use `hypothesis` to generate thousands of random inputs
to catch edge cases we didn't think to hand-write. Each test checks a
property that must hold for ALL valid inputs, not just specific examples.

Properties covered:
  1. `_rc` is involutive: _rc(_rc(s)) == s for any IUPAC DNA string.
  2. `_rc` preserves length.
  3. `_rc` output contains only IUPAC characters if input does.
  4. `_feat_len` returns the same value as the non-wrap arithmetic when
     start <= end, and the wrap formula when end < start.
  5. `_feat_len` is always non-negative and <= total.
  6. `_bp_in` is consistent with `_feat_len`: for a non-wrap feature of
     length L, exactly L positions in [0,total) satisfy _bp_in.
  7. The wrap-aware midpoint formula always lies on the feature's arc.
"""
from __future__ import annotations

import pytest

from hypothesis import given, strategies as st, settings, assume

import splicecraft as sc


# ── Strategy helpers ──────────────────────────────────────────────────────────

# `_rc` uppercases its input, so involution only holds for uppercase strings.
# Mixed-case tests belong in test_dna_sanity.py.
IUPAC_CHARS = "ACGTRYWSMKBDHVN"

iupac_dna = st.text(alphabet=IUPAC_CHARS, min_size=1, max_size=500)
small_total = st.integers(min_value=4, max_value=10_000)


def _bp_in(bp: int, start: int, end: int, total: int) -> bool:
    """Standalone mirror of PlasmidMap._bp_in (which is a method) so tests
    can call it as a pure function. See splicecraft.py:2245."""
    return (start <= bp < end) if end >= start else (bp >= start or bp < end)


def _arc_positions(start: int, end: int, total: int) -> list[int]:
    """Return the list of bp positions covered by [start, end) on a
    circular plasmid of length total. Handles wrap (end < start)."""
    if end < start:
        return list(range(start, total)) + list(range(0, end))
    return list(range(start, end))


# ── Property: _rc ─────────────────────────────────────────────────────────────

class TestReverseComplementProperties:
    @given(seq=iupac_dna)
    @settings(max_examples=300, deadline=None)
    def test_rc_is_involutive(self, seq):
        """Sacred invariant #3: `_rc` must round-trip for any IUPAC seq.
        Double reverse-complement returns the original string unchanged."""
        assert sc._rc(sc._rc(seq)) == seq

    @given(seq=iupac_dna)
    @settings(max_examples=300, deadline=None)
    def test_rc_preserves_length(self, seq):
        assert len(sc._rc(seq)) == len(seq)

    @given(seq=iupac_dna)
    @settings(max_examples=300, deadline=None)
    def test_rc_output_in_iupac(self, seq):
        """Every character in the output must be a valid IUPAC code.
        Catches regressions where ambiguity codes silently pass through
        un-complemented (producing garbage like 'X' or case drift)."""
        rc = sc._rc(seq)
        for ch in rc:
            assert ch in IUPAC_CHARS, f"RC produced non-IUPAC char {ch!r}"

    @given(seq=st.text(alphabet="ACGT", min_size=1, max_size=500))
    @settings(max_examples=100, deadline=None)
    def test_rc_matches_biopython_for_acgt(self, seq):
        """Cross-validate against Biopython's authoritative reverse
        complement for ACGT-only inputs (Biopython doesn't case-preserve
        the same way for full IUPAC, so we stick to ACGT here)."""
        from Bio.Seq import Seq
        assert sc._rc(seq.upper()) == str(Seq(seq.upper()).reverse_complement())


# ── Property: _feat_len ────────────────────────────────────────────────────────

class TestFeatLenProperties:
    @given(start=st.integers(min_value=0, max_value=999),
           length=st.integers(min_value=1, max_value=999),
           total=small_total)
    @settings(max_examples=300, deadline=None)
    def test_feat_len_matches_linear_when_no_wrap(self, start, length, total):
        """For end >= start, `_feat_len` must equal `end - start`."""
        assume(start < total)
        end = start + length
        assume(end <= total)
        assert sc._feat_len(start, end, total) == length

    @given(total=small_total, data=st.data())
    @settings(max_examples=300, deadline=None)
    def test_feat_len_matches_wrap_formula(self, total, data):
        """For end < start (wrap), `_feat_len` must equal (total - start) + end.
        Generate end first (strictly below total), then start strictly above
        end — so every sample is a valid wrap feature (no filtering)."""
        end = data.draw(st.integers(min_value=0, max_value=total - 2))
        start = data.draw(st.integers(min_value=end + 1, max_value=total - 1))
        assert sc._feat_len(start, end, total) == (total - start) + end

    @given(start=st.integers(min_value=0, max_value=9999),
           end=st.integers(min_value=0, max_value=9999),
           total=small_total)
    @settings(max_examples=300, deadline=None)
    def test_feat_len_nonneg_and_bounded(self, start, end, total):
        """Sacred invariant #8: `_feat_len` must be non-negative and
        never exceed the total plasmid length. Breakages here corrupt
        sort orders and primer-design math."""
        assume(start < total and end < total)
        L = sc._feat_len(start, end, total)
        assert 0 <= L <= total


# ── Property: _bp_in ───────────────────────────────────────────────────────────

class TestBpInProperties:
    @given(total=st.integers(min_value=10, max_value=2000), data=st.data())
    @settings(max_examples=100, deadline=None)
    def test_bp_in_counts_match_feat_len_no_wrap(self, total, data):
        """For non-wrap features, the number of positions in [0,total) that
        return True from `_bp_in` must equal `_feat_len(start, end, total)`."""
        start = data.draw(st.integers(min_value=0, max_value=total - 1))
        end = data.draw(st.integers(min_value=start, max_value=total))
        count = sum(1 for i in range(total) if _bp_in(i, start, end, total))
        assert count == sc._feat_len(start, end, total)

    @given(total=st.integers(min_value=10, max_value=2000), data=st.data())
    @settings(max_examples=100, deadline=None)
    def test_bp_in_counts_match_feat_len_wrap(self, total, data):
        """Same as above, but for wrap features (end < start)."""
        end = data.draw(st.integers(min_value=0, max_value=total - 2))
        start = data.draw(st.integers(min_value=end + 1, max_value=total - 1))
        count = sum(1 for i in range(total) if _bp_in(i, start, end, total))
        assert count == sc._feat_len(start, end, total)


# ── Property: circular midpoint (sacred invariant #5) ─────────────────────────

class TestWrapMidpointProperties:
    @given(data=st.data(), total=small_total)
    @settings(max_examples=300, deadline=None)
    def test_midpoint_lies_on_arc(self, data, total):
        """Sacred invariant #5: label-placement midpoint must lie on the
        feature's arc. For wrap features, the naive `(start+end)//2` sits
        on the wrong side of the plasmid; the modular formula must not.
        Covers both wrap and non-wrap cases.

        Drawing start/end from `[0, total)` directly (rather than the
        old `assume(start < total and end < total)`) avoids hypothesis'
        `filter_too_much` health-check failures when `total` is small
        and ~99 % of fixed-range inputs would be rejected.
        """
        start = data.draw(st.integers(min_value=0, max_value=total - 1))
        end   = data.draw(st.integers(min_value=0, max_value=total - 1))
        assume(start != end)  # zero-width has no midpoint semantics
        arc_len = (end - start) % total
        mid = (start + arc_len // 2) % total
        valid_positions = set(_arc_positions(start, end, total))
        assert mid in valid_positions, (
            f"midpoint {mid} not on arc for start={start} end={end} "
            f"total={total} (arc_len={arc_len})"
        )

    @given(data=st.data(), total=small_total)
    @settings(max_examples=200, deadline=None)
    def test_midpoint_is_within_feat_len(self, data, total):
        """The distance from start to midpoint (along the arc, modular)
        must be less than `_feat_len`. Catches off-by-one mistakes that
        place the midpoint just outside the arc."""
        start = data.draw(st.integers(min_value=0, max_value=total - 1))
        end   = data.draw(st.integers(min_value=0, max_value=total - 1))
        assume(start != end)
        arc_len = (end - start) % total
        mid = (start + arc_len // 2) % total
        feat_len = sc._feat_len(start, end, total)
        dist_from_start = (mid - start) % total
        assert dist_from_start < feat_len
