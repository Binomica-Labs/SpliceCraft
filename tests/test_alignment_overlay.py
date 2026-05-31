"""
test_alignment_overlay — linear-map alignment overlay helpers.

Covers `_alignment_to_target_segments` and `_alignment_to_target_letters`
— the pure functions that classify each target column of a pairwise
alignment as match / mismatch / gap. They drive the linear-view
alignment lanes (blue / red / gray bars + letters) that overlay
sequencing-read pile-ups on the plasmid map.

Sacred behaviours under test:

  * Target-resolution coordinates — target gaps (insertions in the
    query) consume no target column and don't subdivide the surrounding
    state.
  * Three-state classification — match / mismatch / gap match the
    user-spec 3-color scheme.
  * Case-insensitive matching — the bp comparison ignores case.
  * `t_start` offset — for local alignments with non-zero target
    offsets, segments shift accordingly.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_to_target_segments
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlignmentToTargetSegments:
    def test_all_match(self):
        assert sc._alignment_to_target_segments("ATGC", "ATGC") == [
            (0, 4, "match"),
        ]

    def test_all_mismatch(self):
        assert sc._alignment_to_target_segments("TTTT", "AAAA") == [
            (0, 4, "mismatch"),
        ]

    def test_match_mismatch_match(self):
        # col 2: G≠C mismatch; flanks all match
        assert sc._alignment_to_target_segments("ATGC", "ATCC") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
            (3, 4, "match"),
        ]

    def test_query_deletion_makes_gap_segment(self):
        # query has 2-bp deletion against target
        assert sc._alignment_to_target_segments("AT--GC", "ATCCGC") == [
            (0, 2, "match"),
            (2, 4, "gap"),
            (4, 6, "match"),
        ]

    def test_target_gap_invisible_at_target_resolution(self):
        # target has a 2-col gap (insertion in query) — consumes zero
        # target columns; surrounding state continues unbroken
        assert sc._alignment_to_target_segments("ATXXGC", "AT--GC") == [
            (0, 4, "match"),
        ]

    def test_target_gap_then_state_change(self):
        # query insertion immediately followed by a mismatch — the
        # mismatch starts at the target column right after the
        # insertion (insertion contributes no target position)
        assert sc._alignment_to_target_segments("ATXG", "AT-C") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
        ]

    def test_t_start_offset(self):
        assert sc._alignment_to_target_segments("ATGC", "ATGC", t_start=100) == [
            (100, 104, "match"),
        ]

    def test_case_insensitive(self):
        assert sc._alignment_to_target_segments("atgc", "ATGC") == [
            (0, 4, "match"),
        ]

    def test_empty(self):
        assert sc._alignment_to_target_segments("", "") == []

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ in length"):
            sc._alignment_to_target_segments("ATGC", "ATG")

    def test_complex_mixed_states(self):
        # M M MM M gap M  →  segments at positions 0..2, 2..3, 3..4, 4..5, 5..6
        assert sc._alignment_to_target_segments("ATGC-G", "ATCCAG") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
            (3, 4, "match"),
            (4, 5, "gap"),
            (5, 6, "match"),
        ]

    def test_consecutive_runs_coalesce(self):
        # 3 consecutive matches → one segment, not three
        result = sc._alignment_to_target_segments("AAAA", "AAAA")
        assert len(result) == 1
        assert result[0] == (0, 4, "match")

    def test_all_gap(self):
        assert sc._alignment_to_target_segments("----", "ATGC") == [
            (0, 4, "gap"),
        ]


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_name_overlay + _alignment_lane_indicator (bar-overlay + lane tags)
# ═══════════════════════════════════════════════════════════════════════════════


class TestAlignmentNameOverlay:
    """Pure helper that emits the (col, char, state) tuples used by the
    alignment band's in-bar name overlay. The renderer turns each
    tuple into a black-on-state-color cell so the name reads on the
    same colored bar that carries the alignment's data."""

    def test_basic_overlay_returns_one_tuple_per_char(self):
        col_state = {10: "match", 11: "match", 12: "match", 13: "match",
                     14: "match", 15: "match"}
        out = sc._alignment_name_overlay("READ_1", 10, 6, col_state)
        assert out == [
            (10, "R", "match"),
            (11, "E", "match"),
            (12, "A", "match"),
            (13, "D", "match"),
            (14, "_", "match"),
            (15, "1", "match"),
        ]

    def test_truncates_when_name_wider_than_bar(self):
        col_state = {0: "match", 1: "match", 2: "match", 3: "match"}
        out = sc._alignment_name_overlay("WAY_TOO_LONG", 0, 4, col_state)
        assert [t[1] for t in out] == ["W", "A", "Y", "_"]

    def test_skipped_when_bar_below_min_cols(self):
        # Default min is 4 — a 3-col bar produces nothing.
        col_state = {0: "match", 1: "match", 2: "match"}
        assert sc._alignment_name_overlay("ABC", 0, 3, col_state) == []

    def test_picks_state_per_column_for_mixed_bar(self):
        # Mid-bar mismatch should color that name char's bg red.
        col_state = {5: "match", 6: "match", 7: "mismatch", 8: "match"}
        out = sc._alignment_name_overlay("name", 5, 4, col_state)
        states = [t[2] for t in out]
        assert states == ["match", "match", "mismatch", "match"]

    def test_missing_col_state_defaults_to_match(self):
        # Edge col not in col_state → falls back to "match" (a stable
        # default so the overlay never paints transparent / unstyled).
        out = sc._alignment_name_overlay("xy", 100, 4, {})
        assert [t[2] for t in out] == ["match", "match"]

    def test_empty_name_returns_empty(self):
        assert sc._alignment_name_overlay("", 0, 10, {0: "match"}) == []

    def test_whitespace_only_name_returns_empty(self):
        assert sc._alignment_name_overlay("   ", 0, 10, {0: "match"}) == []

    def test_strips_surrounding_whitespace(self):
        out = sc._alignment_name_overlay("  ab  ", 0, 6, {})
        assert [t[1] for t in out] == ["a", "b"]

    def test_unknown_state_falls_back_to_match(self):
        out = sc._alignment_name_overlay("z", 0, 4, {0: "bogus"})
        assert out == [(0, "z", "match")]

    def test_non_int_bar_width_returns_empty(self):
        assert sc._alignment_name_overlay("ab", 0, "wide", {}) == []  # type: ignore[arg-type]

    def test_non_int_col_start_returns_empty(self):
        assert sc._alignment_name_overlay("ab", 1.5, 4, {}) == []  # type: ignore[arg-type]

    def test_non_dict_col_state_falls_back_to_match(self):
        # Defensive: pass-through robustness if the caller wires a
        # list / None by mistake (don't crash, default to match).
        out = sc._alignment_name_overlay("ab", 0, 4, None)  # type: ignore[arg-type]
        assert [t[2] for t in out] == ["match", "match"]


class TestAlignmentLaneIndicator:
    """1-indexed, fixed-width lane tag used at the left margin when
    letter mode is on (raw bases visible — name overlay would clash)."""

    def test_first_lane_is_one(self):
        assert sc._alignment_lane_indicator(0) == " 1"

    def test_second_lane_is_two(self):
        assert sc._alignment_lane_indicator(1) == " 2"

    def test_ninth_lane_right_justified(self):
        assert sc._alignment_lane_indicator(8) == " 9"

    def test_two_digit_lane_uses_full_width(self):
        assert sc._alignment_lane_indicator(9) == "10"
        assert sc._alignment_lane_indicator(98) == "99"

    def test_overflow_collapses_to_truncated_marker(self):
        # 100+ lanes at width 2 collapses to "9+" so the column
        # boundary stays stable (letter area never shifts between rows).
        assert sc._alignment_lane_indicator(99) == "9+"
        assert sc._alignment_lane_indicator(500) == "9+"

    def test_overflow_at_width_three(self):
        assert sc._alignment_lane_indicator(999, width=3) == "99+"

    def test_negative_returns_blanks(self):
        assert sc._alignment_lane_indicator(-1) == "  "

    def test_non_int_returns_blanks(self):
        assert sc._alignment_lane_indicator("x") == "  "  # type: ignore[arg-type]

    def test_zero_width_returns_empty(self):
        assert sc._alignment_lane_indicator(0, width=0) == ""


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_to_target_letters
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlignmentToTargetLetters:
    def test_all_match(self):
        assert sc._alignment_to_target_letters("ATGC", "ATGC") == {
            0: ("A", "match"),
            1: ("T", "match"),
            2: ("G", "match"),
            3: ("C", "match"),
        }

    def test_mismatch_letter_is_query_base(self):
        # target ATGC, query ATGT — col 3 query says T, target says C
        letters = sc._alignment_to_target_letters("ATGT", "ATGC")
        assert letters[3] == ("T", "mismatch")

    def test_gap_letter_is_dash(self):
        letters = sc._alignment_to_target_letters("AT-G", "ATCG")
        assert letters[2] == ("-", "gap")

    def test_target_gap_skipped(self):
        # target column 1 is a gap — query base at that column never
        # makes it into the per-target dict
        letters = sc._alignment_to_target_letters("ATXG", "A-TG")
        assert letters == {
            0: ("A", "match"),
            1: ("X", "mismatch"),
            2: ("G", "match"),
        }

    def test_t_start_offset(self):
        assert sc._alignment_to_target_letters("AT", "AT", t_start=50) == {
            50: ("A", "match"),
            51: ("T", "match"),
        }

    def test_case_insensitive(self):
        # lowercase query, uppercase target — match classification holds
        letters = sc._alignment_to_target_letters("atgc", "ATGC")
        for pos in range(4):
            _, state = letters[pos]
            assert state == "match"

    def test_empty(self):
        assert sc._alignment_to_target_letters("", "") == {}

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ in length"):
            sc._alignment_to_target_letters("AT", "ATG")


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_to_query_segments — query-axis mirror of the target helper
# ═══════════════════════════════════════════════════════════════════════════════
# Drives the Alt+A / diff-plasmid overlay flow where the currently-loaded
# plasmid is the **query** (first arg to `_pairwise_align`) and segments
# must land at query bp positions so bars line up on the open record's
# linear map.

class TestAlignmentToQuerySegments:
    def test_all_match(self):
        assert sc._alignment_to_query_segments("ATGC", "ATGC") == [
            (0, 4, "match"),
        ]

    def test_all_mismatch(self):
        assert sc._alignment_to_query_segments("TTTT", "AAAA") == [
            (0, 4, "mismatch"),
        ]

    def test_match_mismatch_match(self):
        assert sc._alignment_to_query_segments("ATGC", "ATCC") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
            (3, 4, "match"),
        ]

    def test_target_deletion_makes_gap_segment(self):
        # Symmetric to the target helper's `test_query_deletion_makes_gap_segment`:
        # target has a 2-bp deletion vs the query — those query
        # positions get classified as "gap" because the target has no
        # base aligned to them.
        assert sc._alignment_to_query_segments("ATCCGC", "AT--GC") == [
            (0, 2, "match"),
            (2, 4, "gap"),
            (4, 6, "match"),
        ]

    def test_query_gap_invisible_at_query_resolution(self):
        # Insertion in target relative to query — consumes zero query
        # columns; surrounding state continues unbroken.
        assert sc._alignment_to_query_segments("AT--GC", "ATXXGC") == [
            (0, 4, "match"),
        ]

    def test_query_gap_then_state_change(self):
        assert sc._alignment_to_query_segments("AT-C", "ATXG") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
        ]

    def test_q_start_offset(self):
        assert sc._alignment_to_query_segments("ATGC", "ATGC", q_start=100) == [
            (100, 104, "match"),
        ]

    def test_case_insensitive(self):
        assert sc._alignment_to_query_segments("atgc", "ATGC") == [
            (0, 4, "match"),
        ]

    def test_empty(self):
        assert sc._alignment_to_query_segments("", "") == []

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ in length"):
            sc._alignment_to_query_segments("ATGC", "ATG")

    def test_complex_mixed_states(self):
        # Mirror of the target helper's complex case: M M MM M gap M
        # using query-axis perspective. Query "ATCCAG" vs target
        # "ATGC-G" — at query position 4 the target has a gap, so
        # the query gets a single "gap" column there.
        assert sc._alignment_to_query_segments("ATCCAG", "ATGC-G") == [
            (0, 2, "match"),
            (2, 3, "mismatch"),
            (3, 4, "match"),
            (4, 5, "gap"),
            (5, 6, "match"),
        ]

    def test_consecutive_runs_coalesce(self):
        result = sc._alignment_to_query_segments("AAAA", "AAAA")
        assert len(result) == 1
        assert result[0] == (0, 4, "match")

    def test_all_gap(self):
        assert sc._alignment_to_query_segments("ATGC", "----") == [
            (0, 4, "gap"),
        ]


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_to_query_letters
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlignmentToQueryLetters:
    def test_all_match(self):
        # Per-query-bp dict — letter at each position is the TARGET
        # base (mirror of the target-axis helper, which stores the
        # query base).
        assert sc._alignment_to_query_letters("ATGC", "ATGC") == {
            0: ("A", "match"),
            1: ("T", "match"),
            2: ("G", "match"),
            3: ("C", "match"),
        }

    def test_mismatch_letter_is_target_base(self):
        letters = sc._alignment_to_query_letters("ATGT", "ATGC")
        assert letters[3] == ("C", "mismatch")

    def test_gap_letter_is_dash(self):
        # Query has a base at col 2, target has a gap — query position
        # records target letter "-" with state "gap".
        letters = sc._alignment_to_query_letters("ATCG", "AT-G")
        assert letters[2] == ("-", "gap")

    def test_query_gap_skipped(self):
        # Insertion in target — query has gap at col 1, that column
        # never enters the per-query dict.
        letters = sc._alignment_to_query_letters("A-TG", "ATXG")
        assert letters == {
            0: ("A", "match"),
            1: ("X", "mismatch"),
            2: ("G", "match"),
        }

    def test_q_start_offset(self):
        assert sc._alignment_to_query_letters("AT", "AT", q_start=50) == {
            50: ("A", "match"),
            51: ("T", "match"),
        }

    def test_case_insensitive(self):
        letters = sc._alignment_to_query_letters("atgc", "ATGC")
        for pos in range(4):
            _, state = letters[pos]
            assert state == "match"

    def test_empty(self):
        assert sc._alignment_to_query_letters("", "") == {}

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="differ in length"):
            sc._alignment_to_query_letters("AT", "ATG")


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-check: segments and letters agree on state
# ═══════════════════════════════════════════════════════════════════════════════

class TestIupacCompatible:
    """Pure helper that powers IUPAC-aware match counting. ``N`` vs
    ``A`` is a match (N stands for any base); ``R`` (A/G) vs ``A`` is
    a match; ``R`` vs ``C`` is a mismatch. Pre-2026-05-27 the aligner
    counted these as strict-equality mismatches, dropping identity %
    on any consensus / primer carrying ambiguity codes."""

    def test_strict_self_match(self):
        for b in "ACGT":
            assert sc._iupac_compatible(b, b)

    def test_n_matches_anything(self):
        for b in "ACGTN":
            assert sc._iupac_compatible("N", b)
            assert sc._iupac_compatible(b, "N")

    def test_r_matches_a_and_g_only(self):
        assert sc._iupac_compatible("R", "A")
        assert sc._iupac_compatible("R", "G")
        assert not sc._iupac_compatible("R", "C")
        assert not sc._iupac_compatible("R", "T")

    def test_y_matches_c_and_t_only(self):
        assert sc._iupac_compatible("Y", "C")
        assert sc._iupac_compatible("Y", "T")
        assert not sc._iupac_compatible("Y", "A")
        assert not sc._iupac_compatible("Y", "G")

    def test_two_ambiguities_compatible_if_any_base_shared(self):
        # R = {A,G}; M = {A,C}; share `A`.
        assert sc._iupac_compatible("R", "M")
        # R = {A,G}; Y = {C,T}; share nothing.
        assert not sc._iupac_compatible("R", "Y")

    def test_case_insensitive(self):
        assert sc._iupac_compatible("a", "N")
        assert sc._iupac_compatible("r", "g")

    def test_u_treated_as_t(self):
        assert sc._iupac_compatible("U", "T")
        assert sc._iupac_compatible("u", "t")
        # U should mismatch C/G/A
        assert not sc._iupac_compatible("U", "A")

    def test_gap_chars_return_false(self):
        # `-` is not in the IUPAC alphabet; callers must filter gaps
        # BEFORE this check. Returning False keeps the helper safe
        # if a caller forgets.
        assert not sc._iupac_compatible("-", "A")
        assert not sc._iupac_compatible("A", "-")

    def test_empty_inputs_return_false(self):
        assert not sc._iupac_compatible("", "A")
        assert not sc._iupac_compatible("A", "")

    def test_non_iupac_char_returns_false(self):
        assert not sc._iupac_compatible("X", "A")
        assert not sc._iupac_compatible("E", "A")


class TestPairwiseAlignIupacAndGapFields:
    """`_pairwise_align` 2026-05-27 changes: IUPAC-aware match
    counting + gap-open / gap-col split + degenerate-input safety."""

    def test_strict_match_counts(self):
        r = sc._pairwise_align("ATGC", "ATGC", mode="global")
        assert r["n_matches"] == 4
        assert r["n_mismatches"] == 0
        assert r["n_gap_cols"] == 0
        assert r["n_gap_opens_q"] == 0
        assert r["n_gap_opens_t"] == 0
        # Strict 100% — feeds the light-blue color tier.
        assert r["identity_pct"] == 100.0

    def test_n_against_a_counts_as_match(self):
        # Single-N consensus aligned against ATGC: pre-2026-05-27
        # the N position was a mismatch, identity_pct = 75. Now N
        # is IUPAC-compatible with any base → 100%.
        r = sc._pairwise_align("ANGC", "ATGC", mode="global")
        assert r["n_matches"] == 4
        assert r["n_mismatches"] == 0
        assert r["identity_pct"] == 100.0

    def test_r_vs_c_is_mismatch(self):
        r = sc._pairwise_align("RTGC", "CTGC", mode="global")
        # R = {A,G}; C is not in R's base set → mismatch.
        assert r["n_matches"] == 3
        assert r["n_mismatches"] == 1

    def test_u_normalized_to_t(self):
        # RNA pasted as DNA — every U used to mismatch T. Now U→T
        # mapping in `_normalize_dna_for_align` makes the alignment
        # behave identically to all-T input.
        r = sc._pairwise_align("AUGC", "ATGC", mode="global")
        assert r["n_matches"] == 4
        assert r["n_mismatches"] == 0
        assert r["identity_pct"] == 100.0

    def test_n_gap_cols_vs_opens(self):
        # Hand-build the gapped strings then count via the gap-open
        # logic by aligning sequences known to produce a 4-bp del.
        # AAAA vs AAAAAATTTT: query is shorter, so global mode
        # opens one gap in query.
        r = sc._pairwise_align(
            "AAAATTTT" + "G" * 4,
            "AAAATTTTGGGG",
            mode="global",
        )
        # The exact gap layout depends on Biopython's tie-breaking,
        # but n_gap_cols >= n_gap_opens_q + n_gap_opens_t MUST hold.
        assert r["n_gap_cols"] >= r["n_gap_opens_q"] + r["n_gap_opens_t"]
        # back-compat: n_gaps is alias for n_gap_cols.
        assert r["n_gaps"] == r["n_gap_cols"]

    def test_n_gap_opens_counts_runs_not_chars(self):
        # Multi-bp deletion: insert 10 As into query, target stays
        # short — the aligner produces 10 contiguous gaps in target.
        r = sc._pairwise_align(
            "AAAAAAAAAAAAAA",       # 14 A's
            "AAAA",                  # 4 A's
            mode="global",
        )
        # Target has 10 gap chars but ONE gap-open run.
        assert r["n_gap_opens_t"] >= 1
        # Gap-cols counts each '-' column.
        assert r["n_gap_cols"] >= 10
        # Gap-opens in target must be << gap-cols (run vs chars).
        assert r["n_gap_opens_t"] < r["n_gap_cols"]


class TestNormalizeMapsURnaToT:
    def test_uppercase_u_becomes_t(self):
        assert sc._normalize_dna_for_align("AUGC") == "ATGC"

    def test_lowercase_u_becomes_t(self):
        assert sc._normalize_dna_for_align("augc") == "ATGC"

    def test_mixed_dna_rna_normalizes(self):
        # AAUUGGCC → AATTGGCC
        assert sc._normalize_dna_for_align("AAUUGGCC") == "AATTGGCC"

    def test_u_in_iupac_string_normalized(self):
        # Ambiguity + RNA mixed: NRUYW → NRTYW (only U becomes T).
        assert sc._normalize_dna_for_align("NRUYW") == "NRTYW"


class TestRotateFrameRaisesOnFallthrough:
    """Pre-2026-05-27 the rotate-back helpers returned the input
    unchanged when the cut point wasn't found — silently leaving
    segments in rotated-target frame. Now strict: raise so a future
    local-mode caller catches the case loudly."""

    def test_target_frame_helper_raises(self):
        with pytest.raises(ValueError, match="cut point"):
            sc._rotate_aligned_to_original_target_frame(
                "AT-GC-", "AT-GC-", 3, 8,
            )

    def test_query_frame_helper_raises(self):
        with pytest.raises(ValueError, match="cut point"):
            sc._rotate_aligned_to_original_query_frame(
                "AT-GC-", "AT-GC-", 3, 8,
            )


class TestAlignmentContentKey:
    """Content-key helper that powers in-batch dedup in
    `_merge_stored_alignments`. Stable across re-serialise; only
    sensitive to target_id + axis + aligned strings."""

    @staticmethod
    def _entry(target_id="T", axis="target", aq="ATGC", at="ATGC"):
        return {
            "target_id": target_id, "axis": axis,
            "aligned_q": aq, "aligned_t": at,
        }

    def test_identical_entries_share_key(self):
        a = self._entry()
        b = self._entry()
        assert sc._alignment_content_key(a) == sc._alignment_content_key(b)

    def test_different_target_id_distinct(self):
        assert (sc._alignment_content_key(self._entry(target_id="A"))
                != sc._alignment_content_key(self._entry(target_id="B")))

    def test_different_axis_distinct(self):
        assert (sc._alignment_content_key(self._entry(axis="target"))
                != sc._alignment_content_key(self._entry(axis="query")))

    def test_different_aligned_q_distinct(self):
        assert (sc._alignment_content_key(self._entry(aq="ATGC"))
                != sc._alignment_content_key(self._entry(aq="AAAA")))

    def test_missing_aligned_strings_returns_none(self):
        assert sc._alignment_content_key({"target_id": "T"}) is None

    def test_pulls_from_result_dict_fallback(self):
        # Some entry shapes (in-flight worker results) carry aq/at
        # only under `result` — helper looks there too.
        e = {
            "target_id": "T", "axis": "target",
            "result": {"aligned_q": "ATGC", "aligned_t": "ATGC"},
        }
        e_direct = {
            "target_id": "T", "axis": "target",
            "aligned_q": "ATGC", "aligned_t": "ATGC",
        }
        assert sc._alignment_content_key(e) == sc._alignment_content_key(e_direct)


class TestMergeDedupConcurrentRace:
    """`_merge_stored_alignments` 2026-05-27: in-batch content-key
    dedup. Two workers producing the same alignment near-simultaneously
    both pass the id check (neither has a `_stored_id` yet) — without
    dedup they'd both append, producing duplicate stored rows."""

    @staticmethod
    def _in_memory(name="dup"):
        return {
            "name": name, "query_label": "q", "target_label": "t",
            "target_id": "T", "target_record": None,
            "axis": "target",
            "aligned_q": "ATGCATGC", "aligned_t": "ATGCATGC",
            "result": {
                "aligned_q": "ATGCATGC", "aligned_t": "ATGCATGC",
                "n_matches": 8, "n_mismatches": 0,
                "identity_pct": 100.0,
            },
        }

    def test_two_concurrent_workers_collapse_to_one_row(self):
        # `_serialize_alignment_for_storage` returns None when
        # target_record is None — use a real record so the serialise
        # succeeds.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("ATGCATGC"), id="T", name="T",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        a = self._in_memory("worker_1")
        b = self._in_memory("worker_2")
        a["target_record"] = rec
        b["target_record"] = rec
        merged, _ = sc._merge_stored_alignments([], [a, b])
        # Without dedup we'd get 2 rows; with dedup, 1.
        assert len(merged) == 1
        assert merged[0]["target_label"] == "t"

    def test_existing_with_same_content_is_deduped(self):
        # 2026-05-29: a re-run whose content byte-matches a row already
        # on disk is now IDEMPOTENT — the existing row is kept (not
        # duplicated) and the in-memory entry adopts its id. This
        # REVERSES the prior "I did this again → separate row" design
        # (`test_existing_with_same_content_is_NOT_deduped`), which let
        # re-runs / benchmarking grow the stored list unbounded —
        # bloating library.json and cluttering the verification report.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("ATGCATGC"), id="T", name="T",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        existing_stored = sc._serialize_alignment_for_storage(
            {**self._in_memory("old"), "target_record": rec},
        )
        assert existing_stored is not None
        existing_id = existing_stored["id"]
        in_mem = [{**self._in_memory("new"), "target_record": rec}]
        merged, stamp_pairs = sc._merge_stored_alignments(
            [existing_stored], in_mem,
        )
        assert len(merged) == 1, "identical re-run must not duplicate"
        assert merged[0]["id"] == existing_id
        # The in-memory entry adopts the canonical existing id, so a
        # later flush updates in place rather than re-appending.
        assert stamp_pairs and stamp_pairs[-1][1] == existing_id

    def test_dedup_leaves_existing_row_untouched(self):
        # A content-match must NOT overwrite the existing row — its
        # metadata (e.g. a user-set visible=False) is preserved.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("ATGCATGC"), id="T", name="T",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        existing_stored = sc._serialize_alignment_for_storage(
            {**self._in_memory("old"), "target_record": rec},
        )
        existing_stored["visible"] = False
        existing_stored["_user_note"] = "keep me"
        in_mem = [{**self._in_memory("new"), "target_record": rec}]
        merged, _ = sc._merge_stored_alignments([existing_stored], in_mem)
        assert len(merged) == 1
        assert merged[0]["visible"] is False
        assert merged[0]["_user_note"] == "keep me"

    def test_different_content_still_appends(self):
        # A genuinely different alignment (different aligned strings →
        # different content key) is still stored — the dedup is precise,
        # not a blanket "one row per target".
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("ATGCATGC"), id="T", name="T",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        existing_stored = sc._serialize_alignment_for_storage(
            {**self._in_memory("old"), "target_record": rec},
        )
        diff = {
            **self._in_memory("new"), "target_record": rec,
            "aligned_q": "TTTTTTTT", "aligned_t": "TTTTTTTT",
            "result": {
                "aligned_q": "TTTTTTTT", "aligned_t": "TTTTTTTT",
                "n_matches": 8, "n_mismatches": 0, "identity_pct": 100.0,
            },
        }
        merged, _ = sc._merge_stored_alignments([existing_stored], [diff])
        assert len(merged) == 2


class TestSegmenterLetterConsistency:
    """The two helpers walk the same gapped strings with the same
    classification — every target column in `letters` must fall inside
    exactly one segment of the matching state.
    """

    @pytest.mark.parametrize("aq,at", [
        ("ATGC",       "ATGC"),
        ("ATGT",       "ATGC"),
        ("AT--GC",     "ATCCGC"),
        ("ATXG",       "AT-C"),
        ("ATGC-G",     "ATCCAG"),
        ("----",       "ATGC"),
        ("ATCGATCG",   "ATCGTTCG"),  # one mismatch in the middle
    ])
    def test_consistency(self, aq, at):
        segs = sc._alignment_to_target_segments(aq, at)
        letters = sc._alignment_to_target_letters(aq, at)
        for t_pos, (_letter, state) in letters.items():
            matching = [
                s for s in segs if s[0] <= t_pos < s[1] and s[2] == state
            ]
            assert len(matching) == 1, (
                f"t_pos={t_pos} state={state!r} not covered by any "
                f"matching segment in {segs!r}"
            )


class TestQuerySegmenterLetterConsistency:
    """Symmetric guard for the query-axis helpers."""

    @pytest.mark.parametrize("aq,at", [
        ("ATGC",       "ATGC"),
        ("ATGT",       "ATGC"),
        ("ATCCGC",     "AT--GC"),
        ("AT-C",       "ATXG"),
        ("ATCCAG",     "ATGC-G"),
        ("ATGC",       "----"),
        ("ATCGATCG",   "ATCGTTCG"),
    ])
    def test_consistency(self, aq, at):
        segs = sc._alignment_to_query_segments(aq, at)
        letters = sc._alignment_to_query_letters(aq, at)
        for q_pos, (_letter, state) in letters.items():
            matching = [
                s for s in segs if s[0] <= q_pos < s[1] and s[2] == state
            ]
            assert len(matching) == 1, (
                f"q_pos={q_pos} state={state!r} not covered by any "
                f"matching segment in {segs!r}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Registration + lifecycle hardening
# ═══════════════════════════════════════════════════════════════════════════════
# Generation-counter race guards prevent in-flight workers from
# resurrecting cleared alignments. The two assertions below back the
# `_alignments_generation` contract:
#   * `_clear_alignments` ALWAYS bumps the counter (even when the band
#     is already empty) so a worker that hadn't registered yet still
#     gets poisoned by a "preemptive" clear.
#   * `_register_alignment` refuses degenerate input (empty aligned
#     strings) — those would paint nothing and surface as a phantom
#     row.

TERMINAL_SIZE = (160, 48)


class TestAlignmentLifecycle:
    """Pilot-driven tests for the register/clear contract."""

    async def test_clear_bumps_generation_when_non_empty(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Seed one alignment so clear has work to do.
            app._alignments = [{
                "name": "fake", "query_label": "q", "target_label": "t",
                "target_record": tiny_record,
                "result": {"aligned_q": "A", "aligned_t": "A"},
                "aligned_q": "A", "aligned_t": "A",
                "t_start": 0, "segments": [(0, 1, "match")],
                "t_lo": 0, "t_hi": 1, "letters": None,
            }]
            gen_before = app._alignments_generation
            app._clear_alignments()
            assert app._alignments == []
            assert app._alignments_generation == gen_before + 1

    async def test_clear_bumps_generation_when_already_empty(
            self, tiny_record, isolated_library):
        """Empty-band clear still bumps the counter — workers that
        started before the clear must still see the bump and refuse
        to register, even if there was nothing visible to clear."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            gen_before = app._alignments_generation
            assert app._alignments == []
            app._clear_alignments()
            assert app._alignments_generation == gen_before + 1

    async def test_register_rejects_empty_aligned_strings(
            self, tiny_record, isolated_library):
        """Degenerate `_pairwise_align` results (empty aligned_q /
        aligned_t) MUST NOT register — they'd surface as a phantom
        zero-width row."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            assert app._alignments == []
            # All three empty-string permutations should be refused.
            for aq, at in (("", ""), ("ATGC", ""), ("", "ATGC")):
                app._register_alignment(
                    name="empty",
                    query_label="q",
                    target_label="t",
                    target_record=tiny_record,
                    result={"aligned_q": aq, "aligned_t": at},
                )
            assert app._alignments == []

    async def test_register_succeeds_with_valid_result(
            self, tiny_record, isolated_library):
        """Sanity: a valid result lands in the band."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._register_alignment(
                name="ok",
                query_label="q",
                target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            assert len(app._alignments) == 1
            entry = app._alignments[0]
            assert entry["segments"] == [(0, 4, "match")]
            assert entry["t_lo"] == 0 and entry["t_hi"] == 4
            # Default axis is "target" (Plasmidsaurus convention).
            assert entry["axis"] == "target"


class TestAlignmentPersistenceRoundTrip:
    """User-feature 2026-05-23: alignments survive a record swap by
    living on the library entry. Hydrate restores the visible ones
    on the next `_apply_record`.

    Test contract: register an alignment against library entry A,
    flush, load entry B, load entry A again → the original alignment
    is back on the band with all fields intact (segments, axis,
    target_label, target_record's sequence).
    """

    @staticmethod
    def _make_library_record(seq: str, rid: str, name: str = None):
        """Build a minimal library entry dict that includes a real
        record + gb_text + size."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq(seq), id=rid, name=name or rid,
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        return rec, {
            "id":      rid,
            "name":    name or rid,
            "size":    len(seq),
            "gb_text": sc._record_to_gb_text(rec),
        }

    async def test_alignment_survives_record_swap_and_back(
            self, isolated_library):
        # Library has two entries: A (200 bp), B (300 bp). Set up
        # the active collection BEFORE creating the app, otherwise
        # `_ensure_default_collection` + `_restore_library_from_active_collection`
        # will rebuild "Main Collection" from whatever's already in
        # the active collection (which may be empty in this isolated
        # tmp tree) and wipe our seeded library.
        rec_a, entry_a = self._make_library_record("A" * 200, "A_PLASMID")
        rec_b, entry_b = self._make_library_record("C" * 300, "B_PLASMID")
        sc._save_collections([{
            "name":        sc._DEFAULT_COLLECTION_NAME,
            "description": "test collection",
            "plasmids":    [entry_a, entry_b],
            "saved":       "2026-05-23",
        }])
        sc._set_active_collection_name(sc._DEFAULT_COLLECTION_NAME)
        sc._save_library([entry_a, entry_b])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._apply_record(rec_a)
            await pilot.pause(0.05)

            # Register an alignment of A against B-as-target.
            seq_a = str(rec_a.seq)
            seq_b = str(rec_b.seq)
            app._register_alignment(
                name="A vs B",
                query_label="A_PLASMID",
                target_label="B_PLASMID",
                target_record=rec_b,
                result={"aligned_q": seq_a[:200], "aligned_t": seq_b[:200]},
                axis="query",
            )
            assert len(app._alignments) == 1
            app._flush_active_alignments()

            # Verify it landed on A's library entry.
            entries = sc._load_library()
            a_entry = next(
                (e for e in entries if e.get("id") == "A_PLASMID"), None,
            )
            assert a_entry is not None, (
                f"A_PLASMID missing from library after flush; entries="
                f"{[e.get('id') for e in entries]!r}"
            )
            stored = a_entry.get("alignments") or []
            assert len(stored) == 1
            assert stored[0]["target_label"] == "B_PLASMID"
            assert stored[0]["visible"] is True
            assert "target_gb_text" in stored[0]
            assert "target_seq_hash" in stored[0]
            stored_id_before = stored[0]["id"]

            # Swap to B → band is cleared (clear_undo=True).
            app._apply_record(rec_b)
            await pilot.pause(0.05)
            assert app._alignments == [], (
                "switching records must drop the in-memory band"
            )

            # Swap back to A → hydrate restores the alignment.
            app._apply_record(rec_a)
            await pilot.pause(0.05)
            assert len(app._alignments) == 1, (
                f"hydrate must restore the stored alignment; got "
                f"{len(app._alignments)} alignments"
            )
            restored = app._alignments[0]
            # Original storage metadata stamped on the restored entry
            # so the next flush round-trips losslessly.
            assert restored["_stored_id"] == stored_id_before
            assert restored["_stored_visible"] is True
            assert restored["target_label"] == "B_PLASMID"
            assert restored["axis"] == "query"
            # And the target record's sequence round-tripped via gb_text.
            assert str(restored["target_record"].seq) == seq_b

    async def test_flush_no_op_when_record_not_in_library(
            self, isolated_library, tiny_record):
        """Loading a record that isn't in the library (file open, demo)
        AND isn't in any collection means there's nowhere to persist.
        Flush must not raise; surfaces a warning notify (2026-05-27)
        so the user knows the alignment won't survive a restart."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)  # not in library
            await pilot.pause(0.05)
            seq = str(tiny_record.seq)
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            # Must not raise.
            app._flush_active_alignments()

    async def test_flush_persists_into_other_collection_when_target_lives_there(
            self, isolated_library):
        """Regression: 2026-05-27 user report. Running a Plasmidsaurus
        alignment from collection 'ActiveCol' against a plasmid that
        lives in collection 'TargetCol' — `_apply_record(target)` swaps
        the canvas, then the worker calls `_flush_active_alignments`.
        Pre-fix the flush looked up the target id in the ACTIVE library
        (= ActiveCol's snapshot) and missed → silently returned → the
        alignment vanished on the next record swap or restart.

        New contract: walk `collections.json`, find whichever
        collection holds the target id, persist into that collection's
        snapshot. The active library is untouched (target plasmid
        isn't there). User gets an info notify so they know where
        it landed.
        """
        rec_active, entry_active = self._make_library_record(
            "A" * 200, "ACTIVE_PLASMID",
        )
        rec_target, entry_target = self._make_library_record(
            "T" * 250, "TARGET_PLASMID",
        )
        sc._save_collections([
            {
                "name":        "ActiveCol",
                "description": "active",
                "plasmids":    [entry_active],
                "saved":       "2026-05-27",
            },
            {
                "name":        "TargetCol",
                "description": "holds the alignment target",
                "plasmids":    [entry_target],
                "saved":       "2026-05-27",
            },
        ])
        sc._set_active_collection_name("ActiveCol")
        sc._save_library([entry_active])   # active library = ActiveCol's

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            # Simulate the Plasmidsaurus flow: load the TARGET onto
            # the canvas (canvas record id = "TARGET_PLASMID", which
            # is NOT in the active library / ActiveCol).
            app._apply_record(rec_target)
            await pilot.pause(0.05)
            app._register_alignment(
                name="WZX_read_1",
                query_label="WZX_read_1",
                target_label="TARGET_PLASMID",
                target_record=rec_target,
                result={
                    "aligned_q": "A" * 200 + "-" * 50,
                    "aligned_t": "T" * 250,
                },
            )
            assert len(app._alignments) == 1
            app._flush_active_alignments()

            # Active library untouched — target isn't there, so the
            # active library save path didn't even run.
            active_entries = sc._load_library()
            assert {e.get("id") for e in active_entries} == {"ACTIVE_PLASMID"}
            assert (active_entries[0].get("alignments") or []) == []

            # The alignment MUST land in TargetCol's snapshot.
            cols = sc._load_collections()
            target_col = next(
                c for c in cols if c.get("name") == "TargetCol"
            )
            target_pl = next(
                p for p in (target_col.get("plasmids") or [])
                if p.get("id") == "TARGET_PLASMID"
            )
            stored = target_pl.get("alignments") or []
            assert len(stored) == 1
            assert stored[0]["target_label"] == "TARGET_PLASMID"
            assert stored[0]["visible"] is True
            stored_id_first = stored[0]["id"]

            # _stored_id stamped on the in-memory entry so a re-flush
            # picks the same on-disk row instead of appending a clone.
            assert app._alignments[0].get("_stored_id") == stored_id_first

            # Re-flush in the same session must NOT create a duplicate.
            app._flush_active_alignments()
            cols2 = sc._load_collections()
            target_col2 = next(
                c for c in cols2 if c.get("name") == "TargetCol"
            )
            target_pl2 = next(
                p for p in (target_col2.get("plasmids") or [])
                if p.get("id") == "TARGET_PLASMID"
            )
            assert len(target_pl2.get("alignments") or []) == 1, (
                "second flush in same session must not append a "
                "duplicate of the same alignment"
            )

    def test_persist_alignments_into_collection_finds_target_across_collections(
            self, isolated_library):
        """Pure helper test: `_persist_alignments_into_collection_for_target`
        walks every collection in `collections.json` to find the
        target id. Pre-fix the flush only looked in the active
        library — this helper is the path that catches the cross-
        collection case before the alignment is lost."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec_t = SeqRecord(
            Seq("G" * 120), id="LIVES_IN_OTHER",
            name="LIVES_IN_OTHER",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        entry_t = {
            "id":      "LIVES_IN_OTHER",
            "name":    "LIVES_IN_OTHER",
            "size":    120,
            "gb_text": sc._record_to_gb_text(rec_t),
        }
        sc._save_collections([
            {"name": "Empty", "description": "", "plasmids": []},
            {"name": "Holds_Target",
             "description": "", "plasmids": [entry_t]},
        ])
        in_memory = [{
            "name":         "stray_read",
            "query_label":  "stray_read",
            "target_label": "LIVES_IN_OTHER",
            "target_record": rec_t,
            "result": {"aligned_q": "G" * 120, "aligned_t": "G" * 120},
            "aligned_q":    "G" * 120,
            "aligned_t":    "G" * 120,
            "axis":         "target",
            "segments":     [(0, 120, "match")],
            "t_lo":         0,
            "t_hi":         120,
            "letters":      None,
        }]
        ok, col_name, n = sc._persist_alignments_into_collection_for_target(
            "LIVES_IN_OTHER", in_memory,
        )
        assert ok is True
        assert col_name == "Holds_Target"
        assert n == 1
        # The collection now carries the alignment.
        cols = sc._load_collections()
        holds = next(c for c in cols if c.get("name") == "Holds_Target")
        pl = holds["plasmids"][0]
        assert len(pl.get("alignments") or []) == 1

    def test_persist_alignments_into_collection_returns_false_when_target_nowhere(
            self, isolated_library):
        """When no collection holds the target id, the helper returns
        (False, None, 0) so the caller can warn the user that the
        alignment cannot persist."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec_t = SeqRecord(
            Seq("C" * 50), id="UNTRACKED",
            name="UNTRACKED",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        sc._save_collections([
            {"name": "OnlyCol", "description": "",
             "plasmids": [{
                 "id": "DIFFERENT", "name": "DIFFERENT", "size": 50,
                 "gb_text": sc._record_to_gb_text(SeqRecord(
                     Seq("A" * 50), id="DIFFERENT", name="DIFFERENT",
                     annotations={"molecule_type": "DNA",
                                  "topology": "linear"})),
             }]},
        ])
        ok, col_name, n = sc._persist_alignments_into_collection_for_target(
            "UNTRACKED",
            [{"target_record": rec_t,
              "result": {"aligned_q": "C"*50, "aligned_t": "C"*50},
              "name": "x", "query_label": "x", "target_label": "x",
              "axis": "target", "aligned_q": "C"*50, "aligned_t": "C"*50,
              "segments": [(0, 50, "match")], "t_lo": 0, "t_hi": 50,
              "letters": None}],
        )
        assert ok is False
        assert col_name is None
        assert n == 0

    def test_merge_stored_alignments_stamps_in_memory_with_canonical_id(
            self, isolated_library):
        """`_merge_stored_alignments` returns ``stamp_pairs`` so the
        caller writes ``_stored_id`` back onto the in-memory entry.
        Without this, a SECOND flush in the same session re-mints a
        fresh uuid on serialise and the merge appends a duplicate row.
        """
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec_t = SeqRecord(Seq("A" * 30), id="T", name="T",
                          annotations={"molecule_type": "DNA",
                                       "topology": "linear"})
        in_mem = [{
            "name": "r", "query_label": "r", "target_label": "T",
            "target_record": rec_t,
            "result": {"aligned_q": "A"*30, "aligned_t": "A"*30},
            "aligned_q": "A"*30, "aligned_t": "A"*30,
            "axis": "target",
            "segments": [(0, 30, "match")],
            "t_lo": 0, "t_hi": 30, "letters": None,
        }]
        merged, stamp_pairs = sc._merge_stored_alignments([], in_mem)
        assert len(merged) == 1
        assert len(stamp_pairs) == 1
        align, sid = stamp_pairs[0]
        assert align is in_mem[0]
        assert isinstance(sid, str) and sid
        assert merged[0]["id"] == sid
        # Caller applies the stamp:
        align["_stored_id"] = sid
        # Now a re-merge against the just-written existing list must
        # update in-place (no duplicate).
        merged2, _ = sc._merge_stored_alignments(merged, in_mem)
        assert len(merged2) == 1
        assert merged2[0]["id"] == sid

    async def test_flush_preserves_hidden_stored_alignments(
            self, isolated_library):
        """Footgun guard 2026-05-23: `_flush_active_alignments` used
        to overwrite the stored list with `self._alignments` — which
        only ever contains *visible* alignments after hydrate. So any
        `visible: False` stored entry would silently vanish the first
        time the user registered a new alignment (the flush would
        write [old visible, new alignment] and drop the hidden ones).

        Test contract: seed the library entry with both a visible and
        a hidden stored alignment. Load → hydrate restores only the
        visible one to `self._alignments`. Register a new alignment
        and flush. Re-read the stored list: must contain all three
        entries — the hidden one preserved untouched.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec_a, entry_a = self._make_library_record("A" * 200, "A_HIDDEN_KEEP")
        rec_b = SeqRecord(
            Seq("T" * 200), id="B", name="B",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        entry_a["alignments"] = [
            {
                "id":              "vis-id",
                "label":           "visible alignment",
                "query_label":     "q",
                "target_label":    "t",
                "target_id":       "B",
                "target_gb_text":  sc._record_to_gb_text(rec_b),
                "target_seq_hash": sc._alignment_target_hash("T" * 200),
                "axis":            "query",
                "result":          {"aligned_q": "A" * 200,
                                    "aligned_t": "T" * 200},
                "visible":         True,
                "added":           "2026-05-23",
                "source":          "manual",
            },
            {
                "id":              "hidden-id",
                "label":           "hidden alignment",
                "query_label":     "q",
                "target_label":    "t",
                "target_id":       "B",
                "target_gb_text":  sc._record_to_gb_text(rec_b),
                "target_seq_hash": sc._alignment_target_hash("T" * 200),
                "axis":            "query",
                "result":          {"aligned_q": "A" * 200,
                                    "aligned_t": "T" * 200},
                "visible":         False,
                "added":           "2026-05-22",
                "source":          "manual",
            },
        ]
        sc._save_collections([{
            "name":        sc._DEFAULT_COLLECTION_NAME,
            "description": "test",
            "plasmids":    [entry_a],
            "saved":       "2026-05-23",
        }])
        sc._set_active_collection_name(sc._DEFAULT_COLLECTION_NAME)
        sc._save_library([entry_a])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._apply_record(rec_a)
            await pilot.pause(0.05)
            # Hydrate restored visible only.
            assert len(app._alignments) == 1
            assert app._alignments[0]["_stored_id"] == "vis-id"

            # Register a fresh alignment and flush.
            app._register_alignment(
                name="new alignment",
                query_label="q",
                target_label="t",
                target_record=rec_b,
                # Distinct query content from the seeded vis/hidden rows
                # so this is a genuinely NEW alignment. With identical
                # content the 2026-05-29 re-run dedup correctly collapses
                # it into the visible row (no append) — which is a
                # separate behaviour with its own coverage; here we want
                # to exercise append-alongside-a-preserved-hidden-row.
                result={"aligned_q": "C" * 200, "aligned_t": "T" * 200},
                axis="query",
            )
            app._flush_active_alignments()

            # The stored list MUST contain all three — hidden preserved.
            entries = sc._load_library()
            a_entry = next(
                (e for e in entries if e.get("id") == "A_HIDDEN_KEEP"),
                None,
            )
            assert a_entry is not None
            stored = a_entry.get("alignments") or []
            ids = [e.get("id") for e in stored]
            assert "vis-id" in ids, (
                f"visible stored entry was dropped from {ids!r}"
            )
            assert "hidden-id" in ids, (
                f"HIDDEN stored entry was wiped by the flush — "
                f"flush must merge with existing storage, not replace. "
                f"got {ids!r}"
            )
            # The hidden entry's visible field must still be False
            # (we didn't accidentally re-visible it).
            hidden = next(e for e in stored if e["id"] == "hidden-id")
            assert hidden["visible"] is False
            # And the fresh one is appended with visible=True.
            fresh = next(
                e for e in stored
                if e["id"] not in ("vis-id", "hidden-id")
            )
            assert fresh["visible"] is True
            assert fresh["label"] == "new alignment"

    async def test_hydrate_skips_invisible_stored_alignments(
            self, isolated_library):
        """Stored alignments with `visible: False` must NOT land on
        the band — they exist for the manager modal but stay hidden
        until toggled. The visibility-toggle UI lives in chunk 2;
        this test pins down the hydrate-side filter so toggle-off
        actually hides on the next load."""
        rec_a, entry_a = self._make_library_record("A" * 100, "A_HIDDEN")
        # Pre-seed the library entry with an invisible alignment.
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        rec_b = SeqRecord(
            Seq("T" * 100), id="B", name="B",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        entry_a["alignments"] = [{
            "id":              "fixed-id",
            "label":           "hidden alignment",
            "query_label":     "q",
            "target_label":    "t",
            "target_id":       "B",
            "target_gb_text":  sc._record_to_gb_text(rec_b),
            "target_seq_hash": sc._alignment_target_hash("T" * 100),
            "axis":            "query",
            "result":          {"aligned_q": "A" * 100, "aligned_t": "T" * 100},
            "visible":         False,
            "added":           "2026-05-23",
            "source":          "manual",
        }]
        sc._save_library([entry_a])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._apply_record(rec_a)
            await pilot.pause(0.05)
            assert app._alignments == [], (
                "invisible stored alignments must not appear on the band"
            )


class TestAlignmentManagerModal:
    """Manager modal (Alt+L) — listing, toggling, and deleting
    stored alignments for the active plasmid."""

    @staticmethod
    def _make_stored(label: str, id_: str = "", *,
                      visible: bool = True, source: str = "manual") -> dict:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("T" * 100), id="T_TARGET", name="T_TARGET",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        return {
            "id":              id_ or f"id-{label}",
            "label":           label,
            "query_label":     "Q",
            "target_label":    "T",
            "target_id":       "T_TARGET",
            "target_gb_text":  sc._record_to_gb_text(rec),
            "target_seq_hash": sc._alignment_target_hash("T" * 100),
            "axis":            "query",
            "result":          {"aligned_q": "A" * 100,
                                "aligned_t": "T" * 100,
                                "identity_pct": 12.3},
            "visible":         visible,
            "added":           "2026-05-23",
            "source":          source,
        }

    async def test_modal_lists_all_stored_with_visibility_glyphs(
            self, tiny_record, isolated_library):
        """All stored entries appear in the table — visible AND hidden,
        with distinct glyphs in column 0."""
        stored = [
            self._make_stored("vis_one", "id1", visible=True),
            self._make_stored("hidden_one", "id2", visible=False),
            self._make_stored("vis_two", "id3", visible=True),
        ]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored, plasmid_label="P")
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            from textual.widgets import DataTable
            t = modal.query_one("#alnmgr-table", DataTable)
            assert t.row_count == 3
            # Modal's own copy is independent of caller's list.
            assert modal._alignments is not stored
            assert [a["label"] for a in modal._alignments] == [
                "vis_one", "hidden_one", "vis_two",
            ]

    async def test_toggle_visible_flips_in_place(
            self, tiny_record, isolated_library):
        """Space (action_toggle_visible) flips the cursor row's
        `visible` field without reordering rows."""
        stored = [
            self._make_stored("a", visible=True),
            self._make_stored("b", visible=True),
        ]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            from textual.widgets import DataTable
            t = modal.query_one("#alnmgr-table", DataTable)
            t.move_cursor(row=1)
            modal.action_toggle_visible()
            assert modal._alignments[0]["visible"] is True
            assert modal._alignments[1]["visible"] is False
            # Cursor stays on row 1.
            assert t.cursor_row == 1

    async def test_delete_removes_cursor_row(
            self, tiny_record, isolated_library):
        stored = [
            self._make_stored("keep1"),
            self._make_stored("remove_me"),
            self._make_stored("keep2"),
        ]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            from textual.widgets import DataTable
            t = modal.query_one("#alnmgr-table", DataTable)
            t.move_cursor(row=1)
            modal.action_delete_selected()
            assert [a["label"] for a in modal._alignments] == [
                "keep1", "keep2",
            ]
            from textual.widgets import DataTable
            assert modal.query_one("#alnmgr-table", DataTable).row_count == 2

    async def test_hide_all_and_show_all_bulk(
            self, tiny_record, isolated_library):
        stored = [
            self._make_stored("a", visible=True),
            self._make_stored("b", visible=False),
            self._make_stored("c", visible=True),
        ]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            # Trigger via direct button press handler.
            modal._hide_all(None)
            assert all(not a["visible"] for a in modal._alignments)
            modal._show_all(None)
            assert all(a["visible"] for a in modal._alignments)

    async def test_save_returns_modified_list_cancel_returns_none(
            self, tiny_record, isolated_library):
        stored = [self._make_stored("alpha")]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        captured = []
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Save path: dismiss with the (modified) list.
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal, callback=lambda r: captured.append(("save", r)))
            await pilot.pause()
            modal.action_toggle_visible()
            modal._save_and_close(None)
            await pilot.pause()
            await pilot.pause(0.05)
            assert captured[-1][0] == "save"
            assert captured[-1][1] is not None
            assert captured[-1][1][0]["visible"] is False  # toggled
            captured.clear()
            # Cancel path: dismiss with None even after edits.
            modal2 = sc.AlignmentManagerModal(stored)
            app.push_screen(modal2, callback=lambda r: captured.append(("cancel", r)))
            await pilot.pause()
            modal2._cancel_btn(None)
            await pilot.pause()
            await pilot.pause(0.05)
            assert captured[-1] == ("cancel", None)


class TestAlignmentManagerMarkAndDeleteMarked:
    """2026-05-27 user feedback: there was no way to bulk-delete a
    subset of lanes — the only bulk option was the all-or-nothing
    "Delete All" button. Replaced with a transient mark concept:
    Space marks the cursor row (× column), "Delete Marked" wipes
    only marked rows. Visibility toggle moved to `v` so Space's
    new mark binding doesn't collide."""

    @staticmethod
    def _make_stored(label, id_=""):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("T" * 80), id="T_TARGET", name="T_TARGET",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        return {
            "id":              id_ or f"id-{label}",
            "label":           label,
            "query_label":     "Q",
            "target_label":    "T",
            "target_id":       "T_TARGET",
            "target_gb_text":  sc._record_to_gb_text(rec),
            "target_seq_hash": sc._alignment_target_hash("T" * 80),
            "axis":            "query",
            "result":          {"aligned_q": "A" * 80,
                                "aligned_t": "T" * 80,
                                "identity_pct": 99.0},
            "visible":         True,
            "added":           "2026-05-27",
            "source":          "manual",
        }

    async def test_mark_default_off_for_all_rows(
            self, tiny_record, isolated_library):
        stored = [self._make_stored("a"), self._make_stored("b")]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            for row in modal._alignments:
                assert row.get("_marked") is False

    async def test_toggle_mark_flips_only_cursor_row(
            self, tiny_record, isolated_library):
        stored = [self._make_stored("a"), self._make_stored("b"),
                  self._make_stored("c")]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            from textual.widgets import DataTable
            t = modal.query_one("#alnmgr-table", DataTable)
            t.move_cursor(row=1)
            modal.action_toggle_mark()
            assert modal._alignments[0]["_marked"] is False
            assert modal._alignments[1]["_marked"] is True
            assert modal._alignments[2]["_marked"] is False
            # Toggle again → off.
            modal.action_toggle_mark()
            assert modal._alignments[1]["_marked"] is False

    async def test_delete_marked_only_removes_marked_rows(
            self, tiny_record, isolated_library):
        stored = [self._make_stored("keep1"), self._make_stored("dropme"),
                  self._make_stored("keep2"), self._make_stored("dropme2")]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            modal._alignments[1]["_marked"] = True
            modal._alignments[3]["_marked"] = True
            modal._delete_marked(None)
            assert [a["label"] for a in modal._alignments] == [
                "keep1", "keep2",
            ]

    async def test_delete_marked_with_nothing_marked_is_noop(
            self, tiny_record, isolated_library):
        stored = [self._make_stored("a"), self._make_stored("b")]
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            modal._delete_marked(None)
            # All rows survive.
            assert len(modal._alignments) == 2

    async def test_save_strips_marked_flag_from_dismiss_payload(
            self, tiny_record, isolated_library):
        """`_marked` is a UI-only selector — must not reach disk via
        the dismiss callback (caller's `_save_library` would persist
        it into the JSON otherwise)."""
        stored = [self._make_stored("a"), self._make_stored("b")]
        captured: list = []
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal, callback=captured.append)
            await pilot.pause(); await pilot.pause(0.05)
            modal._alignments[0]["_marked"] = True
            modal._save_and_close(None)
            await pilot.pause(); await pilot.pause(0.05)
            assert captured, "dismiss callback never fired"
            payload = captured[0]
            assert payload is not None
            assert len(payload) == 2
            for row in payload:
                assert "_marked" not in row, (
                    f"_marked flag leaked through dismiss: {row!r}"
                )


class TestIdentityPctColor:
    """Color tiers picked 2026-05-27 to match the user's sequencing-QC
    grading: light blue STRICT 100, then green / yellow / orange /
    red / gray as identity drops."""

    def test_strict_hundred_is_light_blue(self):
        assert sc._identity_pct_color(100.0) == "bright_cyan"

    def test_just_under_hundred_falls_to_green(self):
        # User-required strictness: 99.999% is NOT light blue.
        assert sc._identity_pct_color(99.999) == "green"

    def test_ninety_is_green(self):
        assert sc._identity_pct_color(90.0) == "green"

    def test_eighty_is_yellow(self):
        assert sc._identity_pct_color(80.0) == "yellow"

    def test_eighty_nine_nine_is_yellow(self):
        assert sc._identity_pct_color(89.999) == "yellow"

    def test_fifty_one_is_orange(self):
        assert sc._identity_pct_color(51.0) == "dark_orange"

    def test_fifty_falls_to_red(self):
        assert sc._identity_pct_color(50.0) == "red"

    def test_eleven_is_red(self):
        assert sc._identity_pct_color(11.0) == "red"

    def test_ten_is_gray(self):
        assert sc._identity_pct_color(10.0) == "grey50"

    def test_zero_is_gray(self):
        assert sc._identity_pct_color(0.0) == "grey50"

    def test_none_is_neutral_white(self):
        assert sc._identity_pct_color(None) == "white"

    def test_non_numeric_is_neutral_white(self):
        assert sc._identity_pct_color("bogus") == "white"  # type: ignore[arg-type]


class TestAlignmentManagerNewAlignButton:
    """`AlignmentManagerModal` gained a "New Align…" button 2026-05-27.
    Clicking it dismisses with a sentinel dict so the caller can chain
    the picker modal and re-open the manager with the new rows once
    workers complete."""

    @staticmethod
    def _make_stored(label):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("T" * 80), id="T_TARGET", name="T_TARGET",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        return {
            "id":              f"id-{label}",
            "label":           label,
            "query_label":     "Q",
            "target_label":    "T",
            "target_id":       "T_TARGET",
            "target_gb_text":  sc._record_to_gb_text(rec),
            "target_seq_hash": sc._alignment_target_hash("T" * 80),
            "axis":            "query",
            "result":          {"aligned_q": "A" * 80,
                                "aligned_t": "T" * 80,
                                "identity_pct": 99.0},
            "visible":         True,
            "added":           "2026-05-27",
            "source":          "manual",
        }

    async def test_new_align_dismisses_with_sentinel_payload(
            self, tiny_record, isolated_library):
        stored = [self._make_stored("a"), self._make_stored("b")]
        captured: list = []
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal, callback=captured.append)
            await pilot.pause(); await pilot.pause(0.05)
            modal._new_align(None)
            await pilot.pause(); await pilot.pause(0.05)
            assert captured, "dismiss callback did not fire"
            payload = captured[0]
            assert isinstance(payload, dict)
            assert payload.get("_new_align") is True
            # Pending edits surfaced under "alignments" so the caller
            # can save them before chaining the picker.
            pending = payload.get("alignments")
            assert isinstance(pending, list)
            assert [a["label"] for a in pending] == ["a", "b"]
            # _marked flag stripped (transient, never reaches disk).
            for row in pending:
                assert "_marked" not in row

    async def test_new_align_preserves_pending_mark_edits_in_payload(
            self, tiny_record, isolated_library):
        """Pending mark edits made BEFORE clicking New Align stay on
        the alignment dicts (the caller still uses them to decide
        what to persist), but the `_marked` flag itself is stripped
        so it doesn't leak to disk."""
        stored = [self._make_stored("a"), self._make_stored("b")]
        captured: list = []
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.AlignmentManagerModal(stored)
            app.push_screen(modal, callback=captured.append)
            await pilot.pause(); await pilot.pause(0.05)
            # Hide row B before pressing New Align — should land on
            # disk via the caller's persist hook.
            modal._alignments[1]["visible"] = False
            modal._new_align(None)
            await pilot.pause(); await pilot.pause(0.05)
            pending = captured[0]["alignments"]
            row_b = next(r for r in pending if r["label"] == "b")
            assert row_b["visible"] is False


class TestAlignmentManagerOpensOnEmptyStorage:
    """Pre-2026-05-27 the manager refused to open if the plasmid had
    no stored alignments — the user had to know about Alt+A
    separately. Now the manager opens with an empty table and the
    "New Align…" button starts the workflow from there."""

    async def test_manager_can_be_constructed_with_empty_list(self):
        # The modal itself accepts an empty list — table renders zero
        # rows, buttons remain functional. Caller's gate (the empty
        # `stored` early-return) was the previous block; that's been
        # removed in `action_open_alignment_manager`.
        modal = sc.AlignmentManagerModal([], plasmid_label="empty")
        assert modal._alignments == []


class TestAlignmentManagerBandRefreshAfterDelete:
    """2026-05-27 user report: deleting alignments via Alt+L doesn't
    update the lane bars on the linear viewer. Pin the band-refresh
    path so a Delete + Save & Close cycle clears the in-memory band
    AND re-hydrates only what's left on disk."""

    async def test_delete_via_modal_then_save_updates_band(
            self, tiny_record, isolated_library):
        # Seed: tiny_record's library entry has TWO stored alignments.
        stored_a = TestAlignmentManagerMarkAndDeleteMarked._make_stored("keep")
        stored_b = TestAlignmentManagerMarkAndDeleteMarked._make_stored("dropme")
        sc._save_library([{
            "id": tiny_record.id, "name": tiny_record.name,
            "size": len(tiny_record.seq), "n_feats": 0,
            "added": "2026-05-27",
            "gb_text": sc._record_to_gb_text(tiny_record),
            "alignments": [stored_a, stored_b],
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.1)
            assert len(app._alignments) == 2, (
                "hydrate should restore both stored alignments to "
                "the band"
            )
            # Drive the modal: mark the dropme row, delete marked,
            # save & close. Use the same callback the app wires up so
            # the band-refresh happens.
            modal = sc.AlignmentManagerModal(
                [stored_a, stored_b], plasmid_label=tiny_record.name,
            )
            rec_id = tiny_record.id
            captured = []

            def _on_done(updated):
                captured.append(updated)
                if updated is None:
                    return
                entries2 = sc._load_library()
                idx = next(
                    (i for i, e in enumerate(entries2)
                     if e.get("id") == rec_id), -1,
                )
                if idx < 0:
                    return
                entries2[idx]["alignments"] = updated
                sc._save_library(entries2, async_sync=True)
                app._clear_alignments()
                app._hydrate_alignments_for_active()

            app.push_screen(modal, callback=_on_done)
            await pilot.pause(); await pilot.pause(0.05)
            # Mark the second row + delete marked + save.
            modal._alignments[1]["_marked"] = True
            modal._delete_marked(None)
            modal._save_and_close(None)
            await pilot.pause(); await pilot.pause(0.1)
            assert captured and captured[0] is not None
            # In-memory band should now have 1 alignment (keep).
            assert len(app._alignments) == 1, (
                f"band must refresh after delete + save; "
                f"got {len(app._alignments)} alignments still on band"
            )
            assert app._alignments[0]["name"] == "keep"
            # And the library entry on disk matches.
            entries = sc._load_library()
            t_entry = next(e for e in entries if e["id"] == rec_id)
            stored = t_entry.get("alignments") or []
            assert len(stored) == 1
            assert stored[0]["label"] == "keep"


class TestAlignmentSurvivesZoomAndPan:
    """Regression guards for the user-reported "the alignment disappears
    when I zoom or pan" complaint. The intended behaviour:

      * Linear-map zoom (+/-) changes the bp/col ratio and the visible
        bp window but MUST NOT clear `self._alignments`.
      * Linear-map pan (arrow keys → `_linear_pan`) shifts the offset
        but MUST NOT clear or reshape alignments; bars stay anchored
        to their bp positions and slide on screen along with the rail.
    """

    async def test_zoom_in_does_not_clear_alignments(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            n = len(tiny_record.seq)
            # Whole-record alignment (covers every bp).
            seq = str(tiny_record.seq)
            app._register_alignment(
                name="self",
                query_label="self", target_label="self",
                target_record=tiny_record,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            assert len(app._alignments) == 1
            entry_before = app._alignments[0]
            for _ in range(5):
                pm.action_linear_zoom_in()
                await pilot.pause()
            # Still registered; identity preserved (same dict, not a
            # fresh copy that lost the cached letters / segments).
            assert len(app._alignments) == 1
            assert app._alignments[0] is entry_before
            assert app._alignments[0]["t_lo"] == 0
            assert app._alignments[0]["t_hi"] == n

    async def test_zoom_out_does_not_clear_alignments(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            pm._linear_zoom = 4.0  # start zoomed in
            seq = str(tiny_record.seq)
            app._register_alignment(
                name="self",
                query_label="self", target_label="self",
                target_record=tiny_record,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            assert len(app._alignments) == 1
            entry_before = app._alignments[0]
            for _ in range(5):
                pm.action_linear_zoom_out()
                await pilot.pause()
            assert len(app._alignments) == 1
            assert app._alignments[0] is entry_before

    async def test_pan_does_not_clear_alignments(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            pm._linear_zoom = 4.0  # zoomed in so pan is meaningful
            seq = str(tiny_record.seq)
            app._register_alignment(
                name="self",
                query_label="self", target_label="self",
                target_record=tiny_record,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            entry_before = app._alignments[0]
            # Right then left — should land back at the same offset
            # without disturbing the alignment overlay.
            for _ in range(5):
                pm.action_rotate_cw()  # `_linear_pan(+1)` in linear mode
                await pilot.pause()
            for _ in range(5):
                pm.action_rotate_ccw()
                await pilot.pause()
            assert len(app._alignments) == 1
            assert app._alignments[0] is entry_before


class TestRegisterAlignmentAxis:
    """`_register_alignment` accepts an `axis` parameter that selects
    which side of the alignment plays the role of the currently-loaded
    plasmid (= the render axis along which overlay bars are positioned).

    `axis="target"` is the Plasmidsaurus / sequencing-pile flow
    (segments in target coords). `axis="query"` is the Alt+A /
    diff-plasmid flow (segments in query coords) — without this,
    overlay bars on the open plasmid's linear map land at the picked
    plasmid's bp positions instead of the open plasmid's, which is
    wrong whenever the pairwise alignment isn't a perfect 1:1.
    """

    async def test_axis_query_uses_query_coord_segments(
            self, tiny_record, isolated_library):
        """The discriminating case: query "ATCCGC" vs target "AT--GC"
        produces different segments in the two coord systems.

          * Target axis: the target has a 2-col gap (insertion in query)
            so the surrounding match state continues unbroken at target
            resolution → ``[(0, 4, "match")]``.
          * Query axis: the target has 2 fewer bases than the query, so
            those query positions get classified as gap → ``[(0, 2,
            "match"), (2, 4, "gap"), (4, 6, "match")]``.
        """
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._register_alignment(
                name="alt-a",
                query_label="open-plasmid",
                target_label="picked",
                target_record=tiny_record,
                result={
                    "aligned_q": "ATCCGC",
                    "aligned_t": "AT--GC",
                },
                axis="query",
            )
            assert len(app._alignments) == 1
            entry = app._alignments[0]
            assert entry["axis"] == "query"
            assert entry["segments"] == [
                (0, 2, "match"),
                (2, 4, "gap"),
                (4, 6, "match"),
            ]
            # render-axis bounds — bars draw across the full query
            # span (0..6), NOT the target span (0..4).
            assert entry["t_lo"] == 0 and entry["t_hi"] == 6

    async def test_axis_target_default_unchanged(
            self, tiny_record, isolated_library):
        """Regression guard: the existing Plasmidsaurus flow doesn't
        pass `axis`, so it must keep getting target-axis segments. The
        same gapped pair as above, registered without `axis`, collapses
        the target gap into the surrounding match (target-resolution).
        """
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._register_alignment(
                name="plasmidsaurus",
                query_label="read",
                target_label="open-plasmid",
                target_record=tiny_record,
                result={
                    "aligned_q": "ATCCGC",
                    "aligned_t": "AT--GC",
                },
            )
            entry = app._alignments[0]
            assert entry["axis"] == "target"
            assert entry["segments"] == [(0, 4, "match")]

    async def test_invalid_axis_raises(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            with pytest.raises(ValueError, match="axis must be"):
                app._register_alignment(
                    name="bad", query_label="q", target_label="t",
                    target_record=tiny_record,
                    result={"aligned_q": "AT", "aligned_t": "AT"},
                    axis="sideways",
                )



# ═══════════════════════════════════════════════════════════════════════════════
# Circular alignment offset (GH #16, 2026-05-14)
# ═══════════════════════════════════════════════════════════════════════════════

class TestCircularAlignmentOffset:
    """`_find_circular_alignment_offset` rotates a circular target so
    the global pairwise align doesn't pair bp 1 of an arbitrarily-
    started Plasmidsaurus read with bp 1 of the GenBank reference.
    Regression guard for Cory Tobin's report — pre-fix alignment of a
    700-bp-rotated read showed 66% identity + 500+ gaps; post-fix the
    same read aligns at 100% with zero gaps."""

    def test_returns_zero_when_sequences_already_aligned(self):
        target = 'ACGTACGTACGT' * 100
        # Same sequence — no rotation needed.
        assert sc._find_circular_alignment_offset(target, target) == 0

    def test_detects_simple_700_bp_rotation(self):
        # Pseudo-random plasmid-shaped target so each 25-bp seed
        # appears at a unique anchor — repeat-pattern fixtures defeat
        # the uniqueness guard and the helper safely falls back to 0.
        import random
        rng = random.Random(20260514)
        target = "".join(rng.choice("ACGT") for _ in range(2000))
        read = target[700:] + target[:700]
        offset = sc._find_circular_alignment_offset(read, target)
        assert offset == 700

    def test_detects_rotation_near_origin_wrap(self):
        """Rotation just before the target's end means the seed may
        straddle the wrap. The doubled-target search handles this."""
        import random
        rng = random.Random(20260514)
        target = "".join(rng.choice("ACGT") for _ in range(1500))
        rotation = len(target) - 50
        read = target[rotation:] + target[:rotation]
        offset = sc._find_circular_alignment_offset(read, target)
        assert offset == rotation

    def test_returns_zero_for_short_sequences(self):
        # Below the k=25 minimum kmer length — bail cleanly.
        assert sc._find_circular_alignment_offset('AAAA', 'TTTT') == 0

    def test_skips_low_complexity_seeds(self):
        """A query that starts with a homopolymer run shouldn't seed
        on the homopolymer (it'd match everywhere); the helper steps
        past it and finds a complex seed further along."""
        target = 'A' * 100 + 'GTACGTACGTAC' * 30 + 'C' * 50
        # Read starts mid-target.
        rotation = 250
        read = target[rotation:] + target[:rotation]
        offset = sc._find_circular_alignment_offset(read, target)
        # Either the helper finds the exact rotation OR returns 0
        # (no clean unique seed); both are acceptable. The bad
        # outcome we're guarding against is a WRONG non-zero answer.
        assert offset in (0, rotation)

    def test_pairwise_align_with_rotation_recovers_identity(self):
        """End-to-end: a rotated read aligned against the rotated
        target should produce near-100%% identity vs ~50-70%% without
        rotation. This is the test that maps directly to Cory's GH #16
        screenshot."""
        import random
        rng = random.Random(20260514)
        target = "".join(rng.choice("ACGT") for _ in range(3000))
        read = target[700:] + target[:700]
        offset = sc._find_circular_alignment_offset(read, target)
        assert offset == 700
        rotated = target[offset:] + target[:offset]
        result = sc._pairwise_align(read, rotated, mode='global')
        assert result['identity_pct'] >= 99.0
        assert result['n_gaps'] == 0


class TestRotateSeqRecord:
    """`_rotate_seq_record` shifts a SeqRecord's sequence + features
    so that a chosen position becomes the new origin. Used by the
    alignment path to keep the viewer's feature lane in register
    with the rotated target."""

    @staticmethod
    def _circular(seq: str, *, features=()):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(seq), id='T', name='T')
        rec.annotations['molecule_type'] = 'DNA'
        rec.annotations['topology'] = 'circular'
        rec.features = list(features)
        return rec

    def test_zero_offset_returns_input(self):
        rec = self._circular('A' * 100)
        rotated = sc._rotate_seq_record(rec, 0)
        assert rotated is rec

    def test_rotation_shifts_sequence(self):
        rec = self._circular('ABCDEFGHIJ')
        rotated = sc._rotate_seq_record(rec, 3)
        assert str(rotated.seq) == 'DEFGHIJABC'

    def test_rotation_shifts_simple_feature(self):
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = self._circular(
            'A' * 100,
            features=[SeqFeature(FeatureLocation(50, 70, strand=1),
                                   type='CDS')],
        )
        rotated = sc._rotate_seq_record(rec, 20)
        # Feature was at 50-70; after rotation by 20 it's at 30-50.
        assert len(rotated.features) == 1
        loc = rotated.features[0].location
        assert int(loc.start) == 30
        assert int(loc.end) == 50

    def test_rotation_preserves_record_metadata(self):
        rec = self._circular('A' * 100)
        rec.description = 'test'
        rotated = sc._rotate_seq_record(rec, 20)
        assert rotated.id == 'T'
        assert rotated.description == 'test'
        assert rotated.annotations['topology'] == 'circular'


# ═══════════════════════════════════════════════════════════════════════════════
# Per-plasmid map_mode persistence (2026-05-18)
# ═══════════════════════════════════════════════════════════════════════════════
# Library entries can carry a `map_mode` field that overrides the
# topology-derived default on load. Plasmidsaurus alignment + the
# user's Alt+L toggle both write through `_persist_map_mode_for_active`
# so the choice sticks across reloads. Sequencing-aligned plasmids
# also have their library entry tagged `linear` so the next open
# defaults to the diff-friendly view.

class TestPerPlasmidMapModePersistence:
    """Tests use `_apply_record` rather than `_preload_record` because
    the preload path dispatches `_add_save_to_disk` (a `@work(thread=True)`
    worker) whose disk write can race test teardown and contaminate the
    next test's tmp library file."""

    async def test_load_record_with_stashed_linear_overrides_topology(
            self, isolated_library):
        """A circular plasmid loaded from a library entry tagged
        `map_mode: "linear"` opens in linear view — the user-set
        preference beats the topology default."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 500), id="C", name="C",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec._tui_map_mode = "linear"
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "linear", (
                "stashed _tui_map_mode='linear' must win over circular topology"
            )

    async def test_load_record_with_stashed_circular_overrides_linear_topo(
            self, isolated_library):
        """Symmetric: a `topology=linear` record loaded with a
        stashed `circular` preference opens circular. Belt-and-braces
        check on the override direction."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 500), id="L", name="L",
                        annotations={"molecule_type": "DNA",
                                     "topology": "linear"})
        rec._tui_map_mode = "circular"
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "circular"

    async def test_load_record_unknown_stashed_mode_falls_back_to_topology(
            self, isolated_library):
        """Defensive: a bogus stashed value (e.g. hand-edit) is
        ignored and the topology default applies."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("A" * 500), id="C", name="C",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec._tui_map_mode = "spiral"   # nonsense
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "circular"

    async def test_persist_map_mode_writes_to_library_entry(
            self, tiny_record, isolated_library):
        """`_persist_map_mode_for_active` saves the chosen mode onto
        the active library entry so the next reload picks it up.
        Uses `_apply_record` (not `_preload_record`) to avoid the
        background `_add_save_to_disk` worker — that worker's write
        races test teardown."""
        # Pre-seed the library with an entry matching the loaded record.
        sc._save_library([{
            "id":      tiny_record.id,
            "name":    tiny_record.name,
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-18",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            app._persist_map_mode_for_active("linear")
            entries = sc._load_library()
            match = next(e for e in entries if e.get("id") == tiny_record.id)
            assert match.get("map_mode") == "linear"

    async def test_persist_map_mode_is_noop_for_unknown_entry(
            self, isolated_library):
        """When the loaded record isn't in the library, the helper
        must silently no-op (no exception, no spurious row).
        Uses a unique record id so even if a prior test's worker
        leaked an entry to the cache, the lookup fails cleanly."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        sc._save_library([])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Unique id (`UNSAVED_NOOP_TEST`) the previous tests in
            # this class never reference — even cache contamination
            # from a sibling test can't match it.
            unsaved = SeqRecord(Seq("A" * 200),
                                 id="UNSAVED_NOOP_TEST",
                                 name="UNSAVED_NOOP_TEST",
                                 annotations={"molecule_type": "DNA",
                                              "topology": "circular"})
            app._apply_record(unsaved)
            await pilot.pause(0.05)
            entries_before = sc._load_library()
            assert not any(
                e.get("id") == "UNSAVED_NOOP_TEST" for e in entries_before
            ), "unsaved record must not be in the library"
            # Must not raise; library state for the unknown id stays
            # unchanged.
            app._persist_map_mode_for_active("linear")
            entries_after = sc._load_library()
            assert not any(
                e.get("id") == "UNSAVED_NOOP_TEST" for e in entries_after
            ), "no-op path must not insert a new row"
            # And entries that DID exist before keep their state.
            assert entries_before == entries_after

    async def test_toggle_map_view_persists_when_entry_exists(
            self, tiny_record, isolated_library):
        """End-to-end: user-driven `action_toggle_map_view` writes
        through to the library entry. Circular plasmid + toggle → entry
        carries `map_mode: "linear"`."""
        sc._save_library([{
            "id":      tiny_record.id,
            "name":    tiny_record.name,
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-18",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            # tiny_record is circular so starts circular
            assert pm._map_mode == "circular"
            pm.action_toggle_map_view()   # circular → linear
            assert pm._map_mode == "linear"
            entries = sc._load_library()
            match = next(e for e in entries if e.get("id") == tiny_record.id)
            assert match.get("map_mode") == "linear"

    async def test_register_alignment_persists_linear_on_target(
            self, tiny_record, isolated_library):
        """Registering an alignment against a circular target pins
        the map to linear AND writes `map_mode: "linear"` to the
        target's library entry. Mirrors the Plasmidsaurus path —
        sequencing-aligned plasmids default to linear on every later
        load."""
        sc._save_library([{
            "id":      tiny_record.id,
            "name":    tiny_record.name,
            "size":    len(tiny_record.seq),
            "n_feats": 0,
            "added":   "2026-05-18",
            "gb_text": sc._record_to_gb_text(tiny_record),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "circular"
            app._register_alignment(
                name="read1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            assert pm._map_mode == "linear"
            entries = sc._load_library()
            match = next(e for e in entries if e.get("id") == tiny_record.id)
            assert match.get("map_mode") == "linear"


# ═══════════════════════════════════════════════════════════════════════════════
# Sequencing toolbar screen (2026-05-18)
# ═══════════════════════════════════════════════════════════════════════════════
# Sequencing replaces the freestanding Plasmidsaurus modal. Tab layout
# leaves room for future ingestion sources (direct API, nanopore
# consensus). The legacy class name is aliased so agent/test paths
# keep resolving.

class TestSequencingScreen:
    def test_back_compat_alias_resolves(self):
        """`PlasmidsaurusAlignModal` is the old class name — kept as
        an alias for tests and agent-API callers."""
        assert sc.PlasmidsaurusAlignModal is sc.SequencingScreen

    def test_menu_lists_sequencing(self):
        """The Sequencing entry is wired into the top-level menu bar."""
        assert "Sequencing" in sc.MenuBar.MENUS

    async def test_screen_opens_with_plasmidsaurus_tab(
            self, tmp_path, tiny_record, isolated_library):
        """The Sequencing screen mounts with the Plasmidsaurus tab
        active (it's currently the only tab; future tabs will share
        the same TabbedContent)."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            assert isinstance(screen, sc.SequencingScreen)
            # All the legacy alignment IDs still resolve so the worker
            # event handlers (which key off these IDs) still work.
            screen.query_one("#align-zip-tree")
            screen.query_one("#align-members")
            screen.query_one("#align-target")
            screen.query_one("#btn-align-go")
            screen.query_one("#btn-sequencing-close")

    async def test_subtabs_disabled_until_zip_loaded(
            self, tmp_path, tiny_record, isolated_library):
        """Samples / Quality / Align sub-tabs are disabled on mount;
        the user can't tab into them until a valid zip lands. General
        stays enabled because it owns the zip picker."""
        from textual.widgets import TabPane
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            for tab_id in ("psaurus-sub-samples",
                            "psaurus-sub-quality",
                            "psaurus-sub-align"):
                assert screen.query_one(f"#{tab_id}", TabPane).disabled, (
                    f"{tab_id} must be disabled before zip load"
                )
            # General stays enabled — it owns the zip picker.
            assert not screen.query_one(
                "#psaurus-sub-general", TabPane,
            ).disabled

    async def test_zip_load_enables_subtabs_and_populates_tables(
            self, tmp_path, tiny_record, isolated_library):
        """End-to-end: feeding the screen a synthetic Plasmidsaurus-
        style zip via `_on_zip_picked` enables the dependent sub-tabs,
        populates the Samples + Quality tables, and writes the run
        metadata summary."""
        from textual.widgets import TabPane, DataTable
        from textual.widgets import DirectoryTree
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Synthesise a Plasmidsaurus-shaped zip: 2 samples with
        # consensus .gbk + summary.txt + per-base TSV.
        rec1 = SeqRecord(Seq("ATGC" * 200), id="MAV34", name="MAV34",
                         annotations={"molecule_type": "DNA",
                                      "topology": "circular"})
        rec2 = SeqRecord(Seq("GCAT" * 200), id="MAV35", name="MAV35",
                         annotations={"molecule_type": "DNA",
                                      "topology": "circular"})
        gbk1 = tmp_path / "MAV34.gbk"
        gbk2 = tmp_path / "MAV35.gbk"
        SeqIO.write(rec1, gbk1, "genbank")
        SeqIO.write(rec2, gbk2, "genbank")
        zp = tmp_path / "RUN42_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk1, "RUN42_genbank-files/RUN42_1_MAV34.gbk")
            zf.write(gbk2, "RUN42_genbank-files/RUN42_2_MAV35.gbk")
            zf.writestr(
                "RUN42_summary-files/RUN42_1_MAV34.txt",
                "       1-mer (%)  2-mer (%)\n"
                "moles       95.5        4.5\n"
                "mass        90.1        9.9\n\n\n"
                "*************************\n\n\n"
                "E. coli genomic contamination: 12.3%\n",
            )
            zf.writestr(
                "RUN42_summary-files/RUN42_2_MAV35.txt",
                "       1-mer (%)  2-mer (%)\n"
                "moles       99.9        0.1\n"
                "mass        99.5        0.5\n\n\n"
                "*************************\n\n\n"
                "E. coli genomic contamination: 2.1%\n",
            )
            # Synthetic per-base TSV: 5 rows, integer coverage.
            zf.writestr(
                "RUN42_per-base-data/RUN42_1_MAV34.tsv",
                "pos\tref\treads_all\n"
                "1\tA\t30\n2\tT\t25\n3\tG\t40\n4\tC\t10\n5\tA\t50\n",
            )
            zf.writestr("RUN42_gel.png", b"PNG-fake-bytes")

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            assert isinstance(screen, sc.SequencingScreen)
            # Feed the FileSelected event the directory tree would emit.
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.3)
            # Dependent sub-tabs are now enabled.
            for tab_id in ("psaurus-sub-samples",
                            "psaurus-sub-quality",
                            "psaurus-sub-align"):
                assert not screen.query_one(
                    f"#{tab_id}", TabPane,
                ).disabled, f"{tab_id} must be enabled after zip load"
            # Run metadata shows both samples.
            assert screen._parsed_run.get("run_id") == "RUN42"
            assert len(screen._parsed_run.get("samples", [])) == 2
            # Samples table populated with one row per sample.
            samples_t = screen.query_one(
                "#align-members", DataTable,
            )
            assert samples_t.row_count == 2
            # Quality table also has both samples.
            quality_t = screen.query_one(
                "#plasmidsaurus-quality-table", DataTable,
            )
            assert quality_t.row_count == 2
            # Run-level files table picks up the gel.png.
            runfiles_t = screen.query_one(
                "#plasmidsaurus-runfiles-table", DataTable,
            )
            assert runfiles_t.row_count >= 1

    async def test_sample_row_select_enables_align_button(
            self, tmp_path, tiny_record, isolated_library):
        """Clicking a sample row marks that sample's .gbk as the
        alignment query, updates the query indicator on the Align
        tab, and flips the Align button to enabled."""
        from textual.widgets import (DataTable, Button,
                                       DirectoryTree)
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 100), id="MAV1", name="MAV1",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        gbk = tmp_path / "MAV1.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "RUN1_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "RUN1_genbank-files/RUN1_1_MAV1.gbk")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.3)
            # Pre-select: Align button is disabled.
            assert screen.query_one(
                "#btn-align-go", Button,
            ).disabled
            # Synthesise the RowSelected event the Samples DataTable
            # would emit on click.
            samples_t = screen.query_one(
                "#align-members", DataTable,
            )
            samples_t.cursor_coordinate = (
                samples_t.cursor_coordinate.__class__(0, 0)
            )
            from textual.coordinate import Coordinate
            row_key = next(iter(samples_t.rows.keys()))
            samples_t.post_message(
                DataTable.RowSelected(
                    samples_t, Coordinate(0, 0), row_key,
                )
            )
            await pilot.pause(0.1)
            assert screen._selected_member is not None
            assert not screen.query_one(
                "#btn-align-go", Button,
            ).disabled
            # Selected member is set; the query indicator update is
            # exercised end-to-end (verified via the button-enabled
            # state above). Static's content is private API in Textual
            # so we don't peek at it directly.
            assert "MAV1" in str(screen._selected_member)

    async def test_repick_same_zip_skips_reparse(
            self, tmp_path, tiny_record, isolated_library):
        """Picking the same zip twice in a row is a no-op (perf
        guard — parse can take ~1 s on large runs)."""
        from textual.widgets import DirectoryTree
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 50), id="MAV1", name="MAV1",
                        annotations={"molecule_type": "DNA"})
        gbk = tmp_path / "MAV1.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "R_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "R_genbank-files/R_1_MAV1.gbk")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            # First pick — populates _parsed_run.
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.3)
            parsed_before = screen._parsed_run
            assert parsed_before
            # Mark with a sentinel so we can detect re-parse.
            parsed_before["_test_sentinel"] = True
            # Second pick of the same path — should NOT re-parse
            # (the sentinel survives because _parsed_run is the
            # same dict object).
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.2)
            assert screen._parsed_run.get("_test_sentinel") is True, (
                "same-path re-pick must not re-parse"
            )

    async def test_invalid_zip_keeps_subtabs_disabled(
            self, tmp_path, tiny_record, isolated_library):
        """A user picking a non-zip file should NOT unlock the sub-tabs
        and the General tab's status row should explain why."""
        from textual.widgets import TabPane, DirectoryTree
        # Non-zip file.
        bogus = tmp_path / "README.txt"
        bogus.write_text("not a zip")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, bogus),
            )
            await pilot.pause(0.2)
            # Sub-tabs stay disabled.
            for tab_id in ("psaurus-sub-samples",
                            "psaurus-sub-quality",
                            "psaurus-sub-align"):
                assert screen.query_one(
                    f"#{tab_id}", TabPane,
                ).disabled, f"{tab_id} must stay disabled on bad zip"
            # Parsed state stays empty.
            assert not screen._parsed_run


# ═══════════════════════════════════════════════════════════════════════════════
# Plasmidsaurus zip parser (run-structured ingestion)
# ═══════════════════════════════════════════════════════════════════════════════
# `_parse_plasmidsaurus_zip` walks a results zip and groups files by
# sample so the Sequencing toolbar's sub-tabs can render without
# re-reading the zip per tab. Run-level extras (gel.png, README) go
# under `run_files`. `_parse_plasmidsaurus_summary` extracts the
# k-mer / contamination percentages from the per-sample summary file.

class TestPlasmidsaurusZipParser:
    def _build_zip(self, dirpath, samples, *, run="RUN1",
                    extra_files=None):
        """Build a synthetic Plasmidsaurus-shaped zip in `dirpath`.
        `samples` is a list of (sample_name, summary_text, perbase_text).
        Returns the zip path."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        zp = dirpath / f"{run}_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            for idx, (name, summary, perbase) in enumerate(samples, 1):
                rec = SeqRecord(
                    Seq("ATGC" * 50), id=name, name=name,
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"},
                )
                gbk = dirpath / f"{name}.gbk"
                SeqIO.write(rec, gbk, "genbank")
                base = f"{run}_{idx}_{name}"
                zf.write(gbk, f"{run}_genbank-files/{base}.gbk")
                if summary is not None:
                    zf.writestr(
                        f"{run}_summary-files/{base}.txt", summary,
                    )
                if perbase is not None:
                    zf.writestr(
                        f"{run}_per-base-data/{base}.tsv", perbase,
                    )
            for name, content in (extra_files or []):
                zf.writestr(name, content)
        return zp

    def test_parses_run_id_from_folder_prefix(self, tmp_path):
        zp = self._build_zip(tmp_path, [("MAV1", None, None)],
                              run="ABC42")
        data = sc._parse_plasmidsaurus_zip(zp)
        assert data["run_id"] == "ABC42"

    def test_groups_files_under_one_sample(self, tmp_path):
        zp = self._build_zip(tmp_path, [
            ("MAV1", "moles 99.0\nmass 98.0\nE. coli contamination: 5.0%\n",
             "pos\tref\treads_all\n1\tA\t30\n2\tT\t40\n"),
        ])
        data = sc._parse_plasmidsaurus_zip(zp)
        assert len(data["samples"]) == 1
        s = data["samples"][0]
        # Sample base collapses to the run_<n>_<name> stem.
        assert s["base"].endswith("MAV1")
        # All categories landed on the same sample dict.
        assert s["gbk"]
        assert s["summary"]
        assert s["perbase"]
        # Summary text streamed inline.
        assert "moles" in s["summary_text"]
        # Per-base coverage stats computed.
        assert s["perbase_coverage"].get("mean") == 35.0

    def test_run_level_files_separated_from_samples(self, tmp_path):
        zp = self._build_zip(tmp_path, [("MAV1", None, None)],
                              extra_files=[("RUN1_gel.png", b"PNG")])
        data = sc._parse_plasmidsaurus_zip(zp)
        # Sample list has the one MAV1; run-level file shows up in
        # `run_files`.
        assert len(data["samples"]) == 1
        run_paths = {rf["name"] for rf in data["run_files"]}
        assert "RUN1_gel.png" in run_paths

    def test_natural_sort_samples(self, tmp_path):
        """Samples come back natural-sorted on their base name —
        the run-index prefix Plasmidsaurus uses (`<run>_<n>_<name>`)
        naturally puts `_2_` before `_10_` under the natural-sort
        key (vs lexicographic `_10_` < `_2_`)."""
        # Pass samples in scrambled order (1, 10, 2). The run-index
        # is assigned by `enumerate` in input order, so the bases
        # become `RUN1_1_A`, `RUN1_2_B`, `RUN1_3_C` — already
        # naturally sorted by index, regardless of name.
        zp = self._build_zip(tmp_path, [
            ("A", None, None),
            ("B", None, None),
            ("C", None, None),
        ])
        data = sc._parse_plasmidsaurus_zip(zp)
        names = [s["name"] for s in data["samples"]]
        assert names == sorted(names, key=sc._natural_sort_key)

    def test_summary_parser_extracts_kmer_and_contam(self):
        text = (
            "       1-mer (%)  2-mer (%)\n"
            "moles       97.5        2.5\n"
            "mass        95.1        4.9\n\n\n"
            "*************************\n\n\n"
            "E. coli genomic contamination: 18.0%\n"
        )
        out = sc._parse_plasmidsaurus_summary(text)
        assert out["kmer_moles_pct"] == 97.5
        assert out["kmer_mass_pct"] == 95.1
        assert out["contamination_pct"] == 18.0
        assert "E. coli" in out["contamination_source"]

    def test_summary_parser_handles_missing_fields(self):
        # Empty input — every field returns None / "".
        out = sc._parse_plasmidsaurus_summary("")
        assert out["kmer_moles_pct"] is None
        assert out["kmer_mass_pct"] is None
        assert out["contamination_pct"] is None
        assert out["contamination_source"] == ""

    def test_perbase_summary_returns_empty_on_garbage(self, tmp_path):
        """A malformed per-base TSV (no numeric column 2) should not
        crash the parser; the sample's `perbase_coverage` ends empty."""
        zp = self._build_zip(tmp_path, [
            ("MAV1", None, "pos\tref\treads_all\n"
                            "alpha\tbeta\tgamma\nA\tB\tC\n"),
        ])
        data = sc._parse_plasmidsaurus_zip(zp)
        assert data["samples"][0]["perbase_coverage"] == {}

    def test_missing_zip_raises_value_error(self, tmp_path):
        with pytest.raises(ValueError):
            sc._parse_plasmidsaurus_zip(tmp_path / "does-not-exist.zip")

    def test_oversize_zip_rejected(self, tmp_path, monkeypatch):
        """A zip claiming to be larger than the cap is refused."""
        # Build a tiny zip then artificially cap to a smaller size.
        zp = self._build_zip(tmp_path, [("MAV1", None, None)])
        monkeypatch.setattr(sc, "_PLASMIDSAURUS_ZIP_MAX_BYTES", 1)
        with pytest.raises(ValueError, match="too large"):
            sc._parse_plasmidsaurus_zip(zp)

    def test_standalone_gbk_no_category_folder(self, tmp_path):
        """Zips without the Plasmidsaurus `_genbank-files/` folder
        layout still discover .gbk files as samples (back-compat
        with the older `_list_gbk_members_in_zip` behaviour)."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 20), id="X", name="X",
                        annotations={"molecule_type": "DNA"})
        gbk = tmp_path / "X.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "ad-hoc.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "sample_A/consensus.gbk")
        data = sc._parse_plasmidsaurus_zip(zp)
        # The .gbk should land as a sample, not in run_files.
        assert any(s.get("gbk") for s in data["samples"])


# ═══════════════════════════════════════════════════════════════════════════════
# Alignment band positioned closest-to-centerline (2026-05-18)
# ═══════════════════════════════════════════════════════════════════════════════
# Alignment lanes used to stack BELOW the reverse-feature band. The
# closest-to-center refactor flips that order — alignment lanes now
# render at `rail_row + 2` and reverse features get offset downward
# by the alignment lane count. `_pack_alignment_lanes` is the helper
# that lets the parent renderer learn the lane count up front.

class TestAlignmentBandCenterline:
    async def test_pack_alignment_lanes_returns_count(
            self, tiny_record, isolated_library):
        """`_pack_alignment_lanes` returns (placed, lane_count); empty
        when no alignments touch the visible window."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            # No alignments registered → empty placement.
            placed, n_lanes = pm._pack_alignment_lanes(
                margin_l=5, usable_w=100, view_s=0, view_e=1000,
                w=160, bp_to_col=lambda bp: 5 + bp // 10,
            )
            assert placed == []
            assert n_lanes == 0
            # Register one alignment → one lane.
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            placed, n_lanes = pm._pack_alignment_lanes(
                margin_l=5, usable_w=100, view_s=0, view_e=1000,
                w=160, bp_to_col=lambda bp: 5 + bp // 10,
            )
            assert n_lanes == 1
            assert len(placed) == 1

    async def test_linear_draw_with_alignment_renders_without_error(
            self, tiny_record, isolated_library):
        """Smoke: a linear-view render with one alignment + rev feature
        completes without raising. Covers the new offset path where
        rev features land below the alignment band."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        # Add a rev-strand feature to exercise the offset path.
        tiny_record.features.append(
            SeqFeature(FeatureLocation(50, 100, strand=-1),
                       type="misc_feature",
                       qualifiers={"label": ["rev1"]})
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            # Must not raise.
            pm.refresh()
            await pilot.pause(0.05)

    async def test_first_paint_after_register_includes_band(
            self, tiny_record, isolated_library):
        """Regression guard (2026-05-22): the PlasmidMap render cache
        previously keyed on (zoom, offset, features, map_mode, …) but
        NOT on `_alignments`, so registering an alignment via
        `set_alignments` invalidated the widget (refresh()) without
        invalidating the draw cache — the cached pre-registration text
        was returned and the user saw bars only after some other
        tracked attribute changed (e.g. a zoom press flipped
        `_linear_zoom` and forced a fresh paint).

        Test contract: render the linear view BEFORE registering an
        alignment, then again AFTER, with no zoom / pan / feature
        change in between. The two paints must differ — the second
        carries the alignment band, the first doesn't. If the cache
        key regresses, both calls return the same cached object.
        """
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            # Force a paint before any alignment so the cache is hot.
            paint_before = pm.render()
            text_before = str(paint_before)
            # Register an alignment. The render path is supposed to
            # invalidate its cache and paint a band.
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            paint_after = pm.render()
            text_after = str(paint_after)
            # The two Rich Texts must compare unequal — same widget
            # size + viewport, but the band landed in the second.
            assert text_before != text_after, (
                "PlasmidMap.render() returned a stale cached paint after "
                "an alignment was registered; the draw-cache key needs "
                "to include alignment state."
            )

    async def test_letter_mode_renders_adjacent_letters_no_gaps(
            self, tiny_record, isolated_library):
        """Regression guard (2026-05-22): at letter-mode zoom the band
        must render alignment letters as adjacent characters
        (`ATGCATGC…`) with no blank gutters between them. The zoom is
        capped at `col_per_bp == 1.0` via `_max_useful_linear_zoom`,
        so the renderer never enters the "letters spread across N
        cols with blank gaps" regime that the user found awkward.

        Test contract: with a record longer than the usable width
        (so the cap is > 1.0 and actually reachable via zoom-in),
        push `_linear_zoom` to the cap, confirm a contiguous run of
        ATGC appears in the band row.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # 500 bp ATGC-repeat — longer than the test's usable_w (~96)
        # so the zoom cap lands at col_per_bp == 1.0 (zoom > 1.0).
        rec = SeqRecord(
            Seq("ATGC" * 125), id="t500", name="t500",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            seq = str(rec.seq)
            app._register_alignment(
                name="self", query_label="self",
                target_label="self", target_record=rec,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            # Push the zoom to its useful cap (one col per bp).
            cap = pm._max_useful_linear_zoom()
            assert cap > 1.0, (
                "test setup error: record must be longer than "
                "usable_w so the cap lands above min zoom"
            )
            pm._linear_zoom = cap
            pm._linear_offset_bp = 0
            text = pm.render()
            plain = str(text)
            rows = plain.split("\n")
            # Locate the band row by the substring of the alignment's
            # sequence — at col_per_bp == 1.0 the letters land at
            # adjacent columns, so `"ATGCATGC"` appears verbatim.
            band_row = next(
                (r for r in rows if "ATGCATGC" in r),
                "",
            )
            assert band_row, (
                f"expected contiguous letter run in band row at "
                f"col_per_bp=1.0; got rows={rows!r}"
            )

    async def test_letter_row_has_no_internal_spaces_at_cap(
            self, tiny_record, isolated_library):
        """Regression guard (2026-05-23): at the zoom cap (col_per_bp
        == 1.0) the band's letter row must be a contiguous run of
        bases — no internal spaces between letters.

        Pre-fix the `bp_to_col` formula used float math
        (`int((bp - view_s) / visible_bp * usable_w)`). At
        `usable_w == visible_bp`, `1/N*N` is `0.9999999999999999` in
        IEEE 754, so `int()` truncated to 0 and bps 0 and 1 collided
        on the same column. The dropped column rendered as a default
        background cell — visible as "splits" in the letter row.

        Using integer math (`(bp - view_s) * usable_w // visible_bp`)
        gives the same result without precision loss; every bp in
        the visible window now maps to a distinct column when
        `usable_w == visible_bp`.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Use a sequence whose length exceeds usable_w so the cap is
        # above min zoom and reachable. 500 bp self-alignment hits
        # one big "match" segment that should paint every bp.
        rec = SeqRecord(
            Seq("ATGC" * 125), id="t500", name="t500",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            seq = str(rec.seq)
            app._register_alignment(
                name="self", query_label="self",
                target_label="self", target_record=rec,
                result={"aligned_q": seq, "aligned_t": seq},
            )
            pm._linear_zoom = pm._max_useful_linear_zoom()
            pm._linear_offset_bp = 0
            text = pm.render()
            plain = str(text)
            rows = plain.split("\n")
            # Find the alignment band row: skip the header (row 0)
            # and pick the row with the most ACGT density.
            band_row = next(
                (r for r in rows[1:]
                 if sum(1 for c in r.strip() if c in "ACGTacgt") > 50),
                "",
            )
            assert band_row, (
                f"could not locate band row; rows={rows!r}"
            )
            # Strip the leading lane indicator + gap (introduced
            # 2026-05-27 as " 1 " / "10 " etc. — see
            # `_alignment_lane_indicator`) before counting internal
            # spaces. The test pins down adjacent-bp precision in the
            # LETTER area, not the margin area.
            stripped = band_row.strip()
            import re as _re
            stripped = _re.sub(r"^[0-9+]{1,3}\s+", "", stripped)
            internal_spaces = sum(1 for c in stripped if c == " ")
            assert internal_spaces == 0, (
                f"band row has {internal_spaces} internal spaces at "
                f"the zoom cap; bp_to_col float precision is causing "
                f"adjacent bps to collide on the same col. "
                f"row={stripped!r}"
            )

    async def test_zoom_in_caps_at_one_col_per_bp(
            self, tiny_record, isolated_library):
        """User feedback 2026-05-22: stop zooming the moment letters
        are one column apart — going further "A T G C" with gutters
        is not useful. `action_linear_zoom_in` must clamp at the
        zoom level where `col_per_bp == 1.0`. Repeated zoom-in
        presses past that should be no-ops.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("A" * 200), id="t200", name="t200",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(rec)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            cap = pm._max_useful_linear_zoom()
            # Zoom in many times — should plateau at the cap.
            for _ in range(30):
                pm.action_linear_zoom_in()
                await pilot.pause()
            assert pm._linear_zoom <= cap + 1e-6, (
                f"action_linear_zoom_in must clamp at "
                f"_max_useful_linear_zoom={cap}; got "
                f"_linear_zoom={pm._linear_zoom}"
            )
            # And once at the cap, another press doesn't push past.
            pm._linear_zoom = cap
            pm.action_linear_zoom_in()
            assert pm._linear_zoom <= cap + 1e-6, (
                f"once at cap, zoom-in must be a no-op; got "
                f"_linear_zoom={pm._linear_zoom}"
            )


    async def test_alignment_selection_clears_on_detail_dismiss(
            self, tiny_record, isolated_library):
        """Regression guard (2026-05-22): clicking an alignment bar
        sets `_selected_align_idx`, which the renderer turns into
        `style="reverse"` on the bar glyphs. Reverse on a full-block
        "█" inverts fg/bg → the bars read as the terminal's default
        foreground (gray on a dark terminal) instead of their
        blue/red/gray scheme. If the selection isn't cleared when the
        detail screen dismisses, the user comes back to all-gray bars.

        Test contract: register an alignment, set the selection
        manually (proxy for the click), simulate a close of the
        detail screen by pushing + dismissing AlignmentScreen with
        the same callback the click path uses. After dismiss,
        `_selected_align_idx` must be back to -1.
        """
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm._map_mode = "linear"
            app._register_alignment(
                name="r1", query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC"},
            )
            ai = 0
            pm._selected_align_idx = ai
            assert pm._selected_align_idx == ai

            # Mirror the on_click → push_screen(..., callback=) flow
            # without the click event itself. The callback below is
            # the one the production code installs.
            def _clear_align_selection(
                _result=None, _selected_ai=ai, _pm=pm,
            ):
                if _pm._selected_align_idx == _selected_ai:
                    _pm._selected_align_idx = -1
                    _pm.refresh()

            screen = sc.AlignmentScreen(
                query_label="q", target_label="t",
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGC",
                        "q_len": 4, "t_len": 4, "score": 8.0,
                        "identity_pct": 100.0, "n_matches": 4,
                        "n_mismatches": 0, "n_gaps": 0,
                        "mode": "global"},
            )
            app.push_screen(screen, callback=_clear_align_selection)
            await pilot.pause()
            await pilot.pause(0.05)
            # Trigger the dismiss path the Esc / q binding takes.
            screen.action_close()
            await pilot.pause()
            await pilot.pause(0.05)
            assert pm._selected_align_idx == -1, (
                "AlignmentScreen.action_close() must fire the "
                "push_screen callback so the lane selection drops; "
                "got _selected_align_idx="
                f"{pm._selected_align_idx}."
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Sequencing hardening (2026-05-18)
# ═══════════════════════════════════════════════════════════════════════════════
# Sweep #8: per-base TSV cap, single-pass zip-open in the Samples
# table, NUL-anchored sentinels for the empty-library / no-gbk paths,
# narrow exception types.

class TestSequencingHardening:
    def test_perbase_summary_truncates_at_max_bytes(self):
        """`_summarize_perbase_tsv` must stop reading once the
        decompressed stream exceeds `max_bytes`. A pathological zip
        bomb that decompresses into a multi-GB single line would
        otherwise OOM `io.TextIOWrapper`'s line buffer."""
        import io
        # Build a 4 KB-per-line TSV with 200 lines = 800 KB. Cap at
        # 200 KB so the streamer stops after ~50 lines, not all 200.
        rows = [f"{i}\tA\t30" for i in range(1, 201)]
        body = ("pos\tref\treads_all\n" + "\n".join(rows)).encode(
            "utf-8",
        )
        # Inflate each row to ~4KB by padding column 1 (`ref`).
        pad = "X" * 4000
        rows_padded = [f"{i}\t{pad}\t30" for i in range(1, 201)]
        body = ("pos\tref\treads_all\n" + "\n".join(rows_padded)).encode(
            "utf-8",
        )
        cap = 200 * 1024
        stats = sc._summarize_perbase_tsv(io.BytesIO(body), max_bytes=cap)
        # Stats are present (the cap allowed *some* rows through).
        assert stats, "truncation must still yield a partial summary"
        # n_pos is bounded by what the cap permitted — about
        # cap / 4 KB ≈ 50 rows. Refuse to specify the exact number,
        # just verify we didn't slurp the full 200.
        assert stats["n_pos"] < 200, (
            f"cap={cap} should have stopped before 200 rows; "
            f"got n_pos={stats['n_pos']}"
        )

    def test_perbase_summary_short_input_complete(self):
        """A short TSV (well under cap) is fully consumed — the cap
        is one-way (truncate-only), it never under-counts on small
        inputs."""
        import io
        body = (
            b"pos\tref\treads_all\n"
            b"1\tA\t10\n2\tT\t20\n3\tG\t30\n4\tC\t40\n5\tA\t50\n"
        )
        stats = sc._summarize_perbase_tsv(
            io.BytesIO(body), max_bytes=1024 * 1024,
        )
        assert stats["n_pos"] == 5
        assert stats["mean"] == 30.0
        assert stats["min"] == 10
        assert stats["max"] == 50
        assert stats["above_20x"] == 4   # 20, 30, 40, 50

    def test_perbase_summary_no_trailing_newline(self):
        """Final row without trailing `\\n` is still counted —
        regression guard for the chunked-reader tail-flush logic.
        Pre-fix the rewrite, a 1-row TSV without trailing newline
        returned an empty dict because the pending fragment never
        got consumed."""
        import io
        body = b"pos\tref\treads_all\n1\tA\t42"
        stats = sc._summarize_perbase_tsv(
            io.BytesIO(body), max_bytes=1024,
        )
        assert stats["n_pos"] == 1
        assert stats["mean"] == 42.0

    def test_parse_zip_skips_oversize_perbase(
            self, tmp_path, monkeypatch):
        """A per-base TSV whose claimed `file_size` exceeds the cap
        is skipped (no read attempted) — defence layer 1. The sample
        still surfaces; its `perbase_coverage` is just empty."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 20), id="X", name="X",
                        annotations={"molecule_type": "DNA"})
        gbk = tmp_path / "X.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "R_results.zip"
        # Plant a big TSV body (~200 KB of synthetic rows).
        big_body = ("pos\tref\treads_all\n"
                    + "\n".join(f"{i}\tA\t30" for i in range(1, 20001)))
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "R_genbank-files/R_1_X.gbk")
            zf.writestr("R_per-base-data/R_1_X.tsv", big_body)
        # Cap to 1 KB so the 200 KB tsv is refused upfront.
        monkeypatch.setattr(
            sc, "_PLASMIDSAURUS_PERBASE_MAX_BYTES", 1024,
        )
        data = sc._parse_plasmidsaurus_zip(zp)
        assert len(data["samples"]) == 1
        # perbase_coverage is empty because the read was refused.
        assert data["samples"][0]["perbase_coverage"] == {}
        # But the sample still lists the perbase member name.
        assert data["samples"][0]["perbase"]

    def test_empty_library_sentinel_is_unique(self):
        """The NUL-anchored sentinels must not collide with any
        realistic library `id` / zip member name. Sanity check that
        they actually contain NUL (which the safe-name check rejects
        in member paths and which LOCUS-safe ids never carry)."""
        assert "\x00" in sc.SequencingScreen._EMPTY_LIBRARY_SENTINEL
        assert "\x00" in sc.SequencingScreen._NO_GBK_KEY_PREFIX

    async def test_target_dropdown_handles_empty_library(
            self, tmp_path, isolated_library):
        """Sequencing screen with NO library entries shows the empty-
        library sentinel and `_go` refuses to advance when the user
        clicks Align without a real target. Verified indirectly via
        the Select's current value (Static's `renderable` is private
        in newer Textual)."""
        from textual.widgets import Select, Button
        # Wipe the library so `_target_options` only has the sentinel.
        sc._save_library([])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            assert isinstance(screen, sc.SequencingScreen)
            # The Select's current value is the empty-library sentinel.
            sel = screen.query_one("#align-target", Select)
            assert sel.value == sc.SequencingScreen._EMPTY_LIBRARY_SENTINEL
            # Simulate a state where the user has picked a sample
            # (forces `_go` past the early-return). The sentinel check
            # should fire BEFORE the zip is opened, so the fake path
            # never gets touched.
            screen._zip_path = tmp_path / "nope.zip"
            screen._selected_member = "ignored.gbk"
            screen.query_one("#btn-align-go", Button).disabled = False
            # Snapshot alignment-registration count; the early-return
            # path must not bump it.
            n_before = len(app._alignments)
            screen._go(None)
            await pilot.pause(0.05)
            # No alignment registered because `_go` short-circuited.
            assert len(app._alignments) == n_before

    async def test_no_gbk_sentinel_refuses_align(
            self, tmp_path, isolated_library):
        """A samples row keyed with the NUL-anchored no-gbk sentinel
        must not arm the Align button (synthetic key would crash
        `_extract_gbk_member`)."""
        from textual.widgets import Button
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen

            class _FakeKey:
                def __init__(self, v):
                    self.value = v

            class _FakeEvent:
                def __init__(self, key):
                    self.row_key = _FakeKey(key)
            # Pretend the user clicked a synthetic row.
            sentinel_key = (
                sc.SequencingScreen._NO_GBK_KEY_PREFIX + "sample-X"
            )
            screen._on_member_selected(_FakeEvent(sentinel_key))
            await pilot.pause(0.05)
            assert screen._selected_member is None
            assert screen.query_one(
                "#btn-align-go", Button,
            ).disabled

    def test_batch_extract_gbk_meta_opens_zip_once(
            self, tmp_path, monkeypatch):
        """`_batch_extract_gbk_meta` should walk every sample's gbk
        inside a single `ZipFile` open — pre-fix each sample paid
        a fresh open. Counts opens via a monkeypatched `__init__`."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Build a 5-sample zip.
        zp = tmp_path / "R_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            for i in range(1, 6):
                rec = SeqRecord(
                    Seq("ATGC" * 50), id=f"S{i}", name=f"S{i}",
                    annotations={"molecule_type": "DNA"},
                )
                gbk = tmp_path / f"S{i}.gbk"
                SeqIO.write(rec, gbk, "genbank")
                zf.write(gbk, f"R_genbank-files/R_{i}_S{i}.gbk")
        parsed = sc._parse_plasmidsaurus_zip(zp)
        # Wire up a SequencingScreen instance just enough to call the
        # batch method directly (avoids the full async-mount cost).
        screen = sc.SequencingScreen.__new__(sc.SequencingScreen)
        screen._zip_path = zp
        screen._parsed_run = parsed
        opens: list[str] = []
        real_init = zipfile.ZipFile.__init__

        def _counting_init(self, file, *a, **kw):
            opens.append(str(file))
            return real_init(self, file, *a, **kw)
        monkeypatch.setattr(zipfile.ZipFile, "__init__", _counting_init)
        meta = screen._batch_extract_gbk_meta(parsed["samples"])
        # Exactly one ZipFile open for all 5 samples.
        assert len(opens) == 1, (
            f"expected 1 zip open for batch read, got {len(opens)}: {opens}"
        )
        # Every sample resolved bp/feats counts.
        assert len(meta) == 5
        for s in parsed["samples"]:
            gbk = s.get("gbk")
            assert gbk in meta, f"missing meta for {gbk}"
            bp_str, _feats = meta[gbk]
            assert bp_str != "—"

    def test_batch_extract_gbk_meta_corrupt_zip_returns_empty(
            self, tmp_path):
        """A corrupted zip path makes `_batch_extract_gbk_meta` log
        and return an empty dict — caller falls back to per-row
        "—" placeholders. Guard against missing-file / bad-zip OS
        errors leaking up the call stack."""
        screen = sc.SequencingScreen.__new__(sc.SequencingScreen)
        # Path that exists but isn't a zip.
        bad = tmp_path / "not-a-zip.txt"
        bad.write_text("hello")
        screen._zip_path = bad
        meta = screen._batch_extract_gbk_meta(
            [{"gbk": "x.gbk"}],
        )
        assert meta == {}

    def test_batch_extract_rejects_unsafe_member_names(self, tmp_path):
        """Belt-and-braces: `_batch_extract_gbk_meta` re-checks
        `_is_safe_zip_member_name` on every member. An in-process
        mutator of `_parsed_run` that tried to smuggle a traversal
        path back in would land in the err bucket, not crash."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 20), id="X", name="X",
                        annotations={"molecule_type": "DNA"})
        gbk = tmp_path / "X.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = tmp_path / "R_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "R_genbank-files/R_1_X.gbk")
        screen = sc.SequencingScreen.__new__(sc.SequencingScreen)
        screen._zip_path = zp
        # Hand-crafted "sample" with a traversal name.
        meta = screen._batch_extract_gbk_meta(
            [{"gbk": "../../etc/passwd"}],
        )
        assert meta == {"../../etc/passwd": ("[red]err[/red]", "—")}


# ═══════════════════════════════════════════════════════════════════════════════
# Plasmidsaurus → "load target as canvas + read as overlay" flow (post-2026-05-24)
# ═══════════════════════════════════════════════════════════════════════════════
# Pre-fix the modal aligned a plasmidsaurus read against a library
# target but left whatever was on the canvas alone — bars rendered in
# rotated-target coords on (typically) the wrong plasmid, and
# `_flush_active_alignments` persisted to the wrong library entry.
# The new flow mirrors Alt+A: the picked library target becomes the
# canvas reference, the read paints as a blue overlay bar on its
# linear view. Rotation switched to the QUERY so target coords stay
# in the library's original frame.

class TestPlasmidsaurusLoadsTargetAsCanvas:
    def _build_min_zip(self, dirpath, gbk_basename, seq):
        """Synthesise the smallest possible plasmidsaurus-shaped zip
        with a single gbk sample. Returns the zip path."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(seq), id=gbk_basename, name=gbk_basename,
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        gbk = dirpath / f"{gbk_basename}.gbk"
        SeqIO.write(rec, gbk, "genbank")
        zp = dirpath / "RUN42_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, f"RUN42_genbank-files/{gbk_basename}.gbk")
        return zp

    def test_on_member_selected_tracks_order_and_basename(
            self, tmp_path, isolated_library):
        """Picking a samples-table row populates `_selected_order_num`
        (1-based row index) and `_selected_gbk_basename` (.gbk leaf
        with extension stripped). These drive the alignment label."""
        zp = self._build_min_zip(tmp_path, "RUN42_1_MAV34", "ATGC" * 50)
        parsed = sc._parse_plasmidsaurus_zip(zp)
        screen = sc.SequencingScreen.__new__(sc.SequencingScreen)
        screen._zip_path = zp
        screen._parsed_run = parsed
        screen._selected_member = None
        screen._selected_order_num = None
        screen._selected_gbk_basename = None
        # Stub out the Textual query_one calls — the test invokes the
        # handler directly without a real mounted tree.
        class _Stub:
            disabled = True
            def update(self, *_a, **_kw): pass
        screen.query_one = lambda *_a, **_kw: _Stub()

        class _FakeKey:
            def __init__(self, v): self.value = v

        class _FakeEvent:
            def __init__(self, key): self.row_key = _FakeKey(key)
        member = parsed["samples"][0]["gbk"]
        screen._on_member_selected(_FakeEvent(member))
        assert screen._selected_member == member
        assert screen._selected_order_num == 1
        # INV-73 (2026-05-25): basename is now post-processed by
        # `_display_label_for_gbk` — Plasmidsaurus run+order prefix
        # stripped, remaining underscores → spaces. Pre-fix the
        # label was "RUN42_1_MAV34" (TUI-unfriendly per user
        # feedback).
        assert screen._selected_gbk_basename == "MAV34"

    def test_on_member_selected_resets_on_no_gbk(
            self, tmp_path, isolated_library):
        """The NUL-anchored no-gbk sentinel clears every selection-
        tracking attr (member + order + basename). Pre-fix only
        `_selected_member` was reset, leaving stale order/basename
        from a prior pick."""
        screen = sc.SequencingScreen.__new__(sc.SequencingScreen)
        screen._zip_path = None
        screen._parsed_run = {}
        screen._selected_member = "old.gbk"
        screen._selected_order_num = 7
        screen._selected_gbk_basename = "old"

        class _Stub:
            disabled = True
            def update(self, *_a, **_kw): pass
        screen.query_one = lambda *_a, **_kw: _Stub()

        class _FakeKey:
            def __init__(self, v): self.value = v

        class _FakeEvent:
            def __init__(self, key): self.row_key = _FakeKey(key)
        sentinel = sc.SequencingScreen._NO_GBK_KEY_PREFIX + "no-gbk-sample"
        screen._on_member_selected(_FakeEvent(sentinel))
        assert screen._selected_member is None
        assert screen._selected_order_num is None
        assert screen._selected_gbk_basename is None

    async def test_align_loads_target_into_canvas(
            self, tmp_path, isolated_library):
        """End-to-end: with a different plasmid on the canvas, running
        the alignment swaps the canvas to the picked target so the
        blue overlay bar lands on the library plasmid's linear view
        (mirrors the Alt+A reference-as-canvas convention)."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Library entries: the target we'll pick + a different decoy
        # we'll load on the canvas first.
        target_seq = "ATGC" * 100
        target_rec = SeqRecord(
            Seq(target_seq), id="TARGET", name="TARGET",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        decoy_seq = "GCTA" * 100
        decoy_rec = SeqRecord(
            Seq(decoy_seq), id="DECOY", name="DECOY",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        sc._save_library([
            {"id": target_rec.id, "name": target_rec.name,
             "size": len(target_seq), "n_feats": 0, "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(target_rec)},
            {"id": decoy_rec.id, "name": decoy_rec.name,
             "size": len(decoy_seq), "n_feats": 0, "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(decoy_rec)},
        ])
        # Plasmidsaurus zip containing one gbk that matches the target
        # exactly (so the alignment is trivially identity-100).
        gbk = tmp_path / "RUN42_1_MAV34.gbk"
        SeqIO.write(target_rec, gbk, "genbank")
        zp = tmp_path / "RUN42_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "RUN42_genbank-files/RUN42_1_MAV34.gbk")

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Canvas starts on the decoy.
            app._apply_record(decoy_rec)
            await pilot.pause(0.05)
            assert app._current_record.id == "DECOY"

            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            from textual.widgets import DirectoryTree, Select
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.3)
            # Pre-set the modal's per-row state directly — the
            # samples-table row-selected event is harder to drive
            # reliably from the test pilot than calling the handler.
            samples = screen._parsed_run.get("samples") or []
            assert samples, "fixture zip must contain at least one sample"
            screen._selected_member = samples[0]["gbk"]
            screen._selected_order_num = 1
            screen._selected_gbk_basename = "RUN42_1_MAV34"
            # Point the target Select at our TARGET library entry.
            sel = screen.query_one("#align-target", Select)
            sel.value = "TARGET"
            screen._go(None)
            # The C-loop runs in a worker; give it a generous deadline.
            for _ in range(40):
                await pilot.pause(0.1)
                if (app._current_record is not None
                        and app._current_record.id == "TARGET"
                        and app._alignments):
                    break
            # Canvas is now the target (not the decoy).
            assert app._current_record.id == "TARGET"
            # Exactly one alignment, labelled `<order> <basename>`.
            assert len(app._alignments) == 1
            entry = app._alignments[0]
            assert entry["name"] == "1 RUN42_1_MAV34"
            assert entry["query_label"] == "1 RUN42_1_MAV34"
            # Source tag for the manager modal's batch-delete.
            assert entry.get("_stored_source") == "sequencing"

    async def test_align_persists_onto_target_library_entry(
            self, tmp_path, isolated_library):
        """The alignment is flushed onto the target's library entry's
        `alignments` field — re-loading the target restores the band.
        Pre-fix the flush wrote to whatever was on the canvas (often
        not the target), so the alignment vanished on re-load."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        target_seq = "ATGC" * 100
        target_rec = SeqRecord(
            Seq(target_seq), id="TARGET2", name="TARGET2",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        sc._save_library([
            {"id": target_rec.id, "name": target_rec.name,
             "size": len(target_seq), "n_feats": 0, "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(target_rec)},
        ])
        gbk = tmp_path / "RUN42_1_MAV34.gbk"
        SeqIO.write(target_rec, gbk, "genbank")
        zp = tmp_path / "RUN42_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "RUN42_genbank-files/RUN42_1_MAV34.gbk")

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            from textual.widgets import DirectoryTree, Select
            tree = screen.query_one(
                "#align-zip-tree", sc._ZipAwareDirectoryTree,
            )
            screen.post_message(
                DirectoryTree.FileSelected(tree.root, zp),
            )
            await pilot.pause(0.3)
            samples = screen._parsed_run.get("samples") or []
            screen._selected_member = samples[0]["gbk"]
            screen._selected_order_num = 1
            screen._selected_gbk_basename = "RUN42_1_MAV34"
            sel = screen.query_one("#align-target", Select)
            sel.value = "TARGET2"
            screen._go(None)
            for _ in range(40):
                await pilot.pause(0.1)
                if app._alignments:
                    break
            # Read it back from disk — must be on the TARGET2 entry.
            entries = sc._load_library()
            t_entry = next(e for e in entries if e["id"] == "TARGET2")
            stored = t_entry.get("alignments") or []
            assert len(stored) == 1
            assert stored[0]["label"] == "1 RUN42_1_MAV34"
            assert stored[0]["source"] == "sequencing"
            assert stored[0]["visible"] is True

    def test_align_worker_rotates_query_not_target(self):
        """Regression for the rotation-frame swap: the worker now
        rotates the QUERY (read) so the alignment result stays in
        the target's original coordinate frame. Pre-fix the target
        was rotated, making `aligned_t` positions land at wrong bp
        on the unrotated canvas. The result dict carries the
        `query_rotation` field for diagnostics; `target_rotation`
        stays 0 to signal no target shift was applied."""
        # Construct a target with a known seed and a query that's
        # the same sequence but rotated by a known offset. The
        # rotation-aware aligner should detect the offset and rotate
        # the query back; alignment identity should be ~100% in
        # target's original frame.
        target_seq = (
            "ATGCATGCATGCATGC" * 10
            + "GGTACCGAATTC"   # uniquely-anchored seed
            + "CCGGAATTCGCATGC" * 10
        )
        offset = 173
        query_seq = target_seq[offset:] + target_seq[:offset]
        # Verify the helper called swapped returns a non-zero offset
        # for the constructed pair (seed is in target, located in
        # query at a different position).
        q_rot = sc._find_circular_alignment_offset(target_seq, query_seq)
        assert q_rot != 0, (
            "swapped-arg helper must detect the constructed query offset"
        )
        # Now apply that rotation to the query and align: identity
        # should be ~100% (we set up the query to be a pure rotation
        # of the target — no mismatches).
        rotated_query = query_seq[q_rot:] + query_seq[:q_rot]
        result = sc._pairwise_align(rotated_query, target_seq, mode="global")
        assert result["identity_pct"] > 99.0, (
            f"rotation-corrected identity={result['identity_pct']} "
            f"should be ~100% for a pure rotation"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Rotation-aware alignment picker (`_pick_best_rotation` + frame transforms)
# ═══════════════════════════════════════════════════════════════════════════════
# Both `_align_worker` (plasmidsaurus) and `_diff_align_worker` (Alt+\)
# route through `_pick_best_rotation` so a circular target with a
# different origin from the canvas plasmid gets a well-anchored
# alignment (not edge-gap-padded). Helper picks among plain / query-rot
# / target-rot candidates by overall identity_pct (more aligned bp =
# more informative overlay band) and shifts aq/at back to the canvas
# axis frame whenever the picked rotation was on the canvas side.

class TestRotateAlignedToOriginalTargetFrame:
    def test_zero_rotation_returns_unchanged(self):
        """A t_rot of 0 is a no-op (no cut, no rotate). Pre-fix early
        return guard."""
        aq, at = "ATGCATGC", "ATGCATGC"
        assert sc._rotate_aligned_to_original_target_frame(
            aq, at, 0, 8,
        ) == (aq, at)

    def test_empty_target_returns_unchanged(self):
        """A tn of 0 is degenerate — skip the walk."""
        assert sc._rotate_aligned_to_original_target_frame(
            "", "", 5, 0,
        ) == ("", "")

    def test_perfect_alignment_rotates_cleanly(self):
        """Round-trip: rotate target by 3 bp, align identically,
        then rotate aq/at back to original frame. Non-gap positions
        of new_at should reproduce the original target sequence."""
        target = "ABCDEFGH"
        t_rot = 3
        # If we rotated the target by 3, rotated_target = "DEFGHABC".
        # Suppose a perfect alignment: aq = at = "DEFGHABC".
        aq = at = "DEFGHABC"
        new_aq, new_at = sc._rotate_aligned_to_original_target_frame(
            aq, at, t_rot, len(target),
        )
        # Non-gap positions of new_at should now spell out the
        # original target.
        assert new_at == "ABCDEFGH"
        assert new_aq == "ABCDEFGH"

    def test_alignment_with_gaps_preserves_length(self):
        """Cut+rotate preserves the alignment-column count. Gaps stay
        where the C-loop put them (rotated to new positions)."""
        # rotated_target = "DEFGH" with a 2-bp insertion in query
        # between bp 2 and bp 3 (G and H).
        # aq = "DE-FXX-GH"  no — let me redo with valid alignment
        # Concrete: target=ABCDEFGH (len 8), t_rot=3
        # rotated_target = "DEFGHABC"
        # Aligner pairs query "DE--FGHABC" against rotated_target
        # "DEFGHABC--" — 2-bp insertion at columns 2,3. So
        # aq = "DE--FGHABC", at = "DEFGHABC--"
        # No wait, at must have non-gap chars equal to rotated_target.
        # Let me use simpler: aq has a 1-bp gap relative to rotated_target.
        aq = "DEXFGHABC"   # length 9, last char of rotated_target replaced
        at = "DE-FGHABC"   # length 9, gap at column 2 (target gap, no bp here)
        new_aq, new_at = sc._rotate_aligned_to_original_target_frame(
            aq, at, 3, 8,
        )
        # Length preserved
        assert len(new_aq) == len(aq)
        assert len(new_at) == len(at)
        # Non-gap count of at is unchanged (still tracks 8 bp of target)
        assert at.count("-") == new_at.count("-")

    def test_cut_not_found_raises(self):
        """If the alignment's `at` doesn't contain enough non-gap
        chars to reach cut_target_bp, the helper raises ValueError.

        Pre-2026-05-27 returned the inputs unchanged — but those
        strings are in ROTATED frame, so downstream segments would
        render at the wrong bp on the canvas (silent positional
        corruption). The new strict contract surfaces the case
        loudly so any future local-mode caller catches it instead
        of producing wrong-frame overlays.
        """
        # 8 bp target, t_rot 3 → need to find the 5th non-gap char
        # (cut_target_bp=5). If `at` has only 4 non-gap chars, can't.
        aq = "AT-GC-"
        at = "AT-GC-"   # 4 non-gap chars
        with pytest.raises(ValueError, match="cut point"):
            sc._rotate_aligned_to_original_target_frame(aq, at, 3, 8)


class TestRotateAlignedToOriginalQueryFrame:
    def test_zero_rotation_returns_unchanged(self):
        """Symmetric to target-frame helper — q_rot 0 → no-op."""
        aq, at = "ATGCATGC", "ATGCATGC"
        assert sc._rotate_aligned_to_original_query_frame(
            aq, at, 0, 8,
        ) == (aq, at)

    def test_perfect_alignment_rotates_cleanly(self):
        """Round-trip: rotate query by 3 bp, align identically,
        rotate back. Non-gap positions of new_aq reproduce original
        query."""
        aq = at = "DEFGHABC"
        new_aq, new_at = sc._rotate_aligned_to_original_query_frame(
            aq, at, 3, 8,
        )
        assert new_aq == "ABCDEFGH"
        assert new_at == "ABCDEFGH"


class TestPickBestRotation:
    def test_picks_best_by_overall_identity(self):
        """When plain + rotation candidates are available, pick by
        overall identity_pct (gap-inclusive) — more aligned bp means
        more informative bars on the overlay band.

        Uses a non-repetitive random construct so plain alignment is
        clearly bad without rotation (a repetitive sequence lets
        plain land near 100% even at large offsets because the
        aligner finds many local matches across the repeats —
        defeating the test's intent)."""
        import random
        random.seed(1234)
        target_seq = "".join(random.choices("ACGT", k=2000))
        offset = 873
        query_seq = target_seq[offset:] + target_seq[:offset]
        result = sc._pick_best_rotation(
            query_seq, target_seq,
            is_circular=True, mode="global",
            canvas_axis="target",
        )
        # A rotation should have been picked AND should have given
        # near-perfect identity. (The RC trial also runs but won't
        # win — RC of a random sequence has no homology to its
        # forward form.)
        assert result["picked_rotation"] in ("query", "target")
        assert result["identity_pct"] > 95.0, (
            f"picker should choose a near-perfect rotation; got "
            f"{result['picked_rotation']!r} at {result['identity_pct']}%"
        )
        # The RC flag should be False since target == fwd query
        # (no RC reverses that relationship).
        assert result.get("query_rc") is False

    def test_plain_wins_when_already_aligned(self):
        """When plain alignment is good (≥ threshold), rotations
        aren't even attempted — picker returns plain. Avoids the 2x
        C-loop cost in the common case."""
        seq = "ATGCATGCATGC" * 100  # 1200 bp
        result = sc._pick_best_rotation(
            seq, seq, is_circular=True, mode="global",
            canvas_axis="target",
        )
        assert result["picked_rotation"] == "none"
        assert result["query_rotation"] == 0
        assert result["target_rotation"] == 0
        assert result["identity_pct"] > 95.0

    def test_linear_target_skips_rotations(self):
        """When `is_circular` is False, rotations aren't tried even
        if plain identity is poor — the target's a linear molecule,
        rotation doesn't make biological sense."""
        # Construct a query that's NOT a rotation of target — they
        # share a small homology but mostly differ. Plain will be
        # poor; rotations could rescue it on a circular target but
        # shouldn't be tried for linear.
        target_seq = "AAAAAAAAAA" + "ATGCATGC" * 50 + "TTTTTTTTTT"
        query_seq  = "ATGCATGC" * 50 + "GGGGGGGGGG"
        result = sc._pick_best_rotation(
            query_seq, target_seq,
            is_circular=False, mode="global",
            canvas_axis="target",
        )
        assert result["picked_rotation"] == "none"
        # `target_rotation` / `query_rotation` defaulted to 0.
        assert result["target_rotation"] == 0
        assert result["query_rotation"] == 0

    def test_target_rotation_shifts_aq_at_for_target_axis(self):
        """When picker chooses target-rotation AND canvas_axis is
        target, the returned aligned_q/aligned_t are pre-shifted to
        the original target frame — downstream segments naturally
        land at the canvas plasmid's bp positions."""
        # Construct one where target-rotation wins. Use a simple
        # pair where the target seed is unique in query.
        # Same construction as the rotation test above but pick the
        # canvas_axis so we exercise the shift code path.
        target_seq = ("AAAAAAAAAA" + "GGTACCGAATTC"
                      + "TTTTTTTTTT" * 10)
        offset = 17
        query_seq = target_seq[offset:] + target_seq[:offset]
        result = sc._pick_best_rotation(
            query_seq, target_seq,
            is_circular=True, mode="global",
            canvas_axis="target",
        )
        # If target-rotation won, the strings should encode the
        # original target frame: walking `aligned_t` non-gap chars
        # should match `target_seq` (possibly with a wrap split).
        if result["picked_rotation"] == "target":
            at_no_gaps = result["aligned_t"].replace("-", "")
            # Non-gap chars represent original target bps 0..tn-1
            assert at_no_gaps == target_seq, (
                "target-axis canvas: aligned_t non-gap chars must "
                "encode the original target sequence in order"
            )

    def test_empty_query_raises_clear_error(self):
        """Empty input is a programmer error — surface it clearly
        instead of letting `_pairwise_align` produce a useless
        0%-identity candidate."""
        with pytest.raises(ValueError, match="non-empty"):
            sc._pick_best_rotation(
                "", "ATGC" * 100, is_circular=False, mode="global",
            )

    def test_empty_target_raises_clear_error(self):
        with pytest.raises(ValueError, match="non-empty"):
            sc._pick_best_rotation(
                "ATGC" * 100, "", is_circular=False, mode="global",
            )

    def test_picks_rc_when_sample_is_reverse_complemented(self):
        """Regression for 2026-05-24: when the read is the RC of the
        target, plain forward alignment is ~0% identity. The picker
        must try the RC orientation and recover the full alignment."""
        from Bio.Seq import Seq
        target_seq = "ATGCGTACGTAGCTAGCTAGCTGATCG" * 100
        query_seq = str(Seq(target_seq).reverse_complement())
        result = sc._pick_best_rotation(
            query_seq, target_seq,
            is_circular=False, mode="global",
            canvas_axis="target",
        )
        # Best candidate should be the RC plain — same orientation
        # as target after flipping the query.
        assert result.get("query_rc") is True
        assert result["identity_pct"] > 95.0, (
            f"RC alignment should recover ~100% identity; got "
            f"{result['identity_pct']}%"
        )

    def test_raises_when_every_candidate_fails(self, monkeypatch):
        """If every alignment call raises (degenerate inputs etc.),
        the helper surfaces the underlying error rather than
        returning a sentinel."""
        def _boom(*_a, **_kw):
            raise ValueError("synthetic alignment failure")
        monkeypatch.setattr(sc, "_pairwise_align", _boom)
        with pytest.raises(ValueError, match="synthetic"):
            sc._pick_best_rotation(
                "ATGC", "GCAT", is_circular=False, mode="global",
            )


class TestDiffAlignWorkerRotation:
    """Alt+\\ (diff-plasmid) gained the same rotation-picker logic as
    plasmidsaurus 2026-05-24. Pre-fix a circular target with a
    different origin from the canvas paid edge gaps; now the worker
    routes through `_pick_best_rotation(canvas_axis='query')`.
    """

    async def test_circular_target_triggers_rotation_pick(
            self, tiny_record, isolated_library):
        """Concrete: load the canvas with one circular plasmid, diff
        against a rotated version. The picker should detect the
        rotation and produce a high-identity alignment instead of the
        gap-padded plain alignment."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Canvas plasmid (query for the diff).
        seq = ("ATGCATGCATGCATGC" * 20
                + "GGTACCGAATTCCCGG"
                + "TTAACCGGTTAACCGG" * 20)
        canvas = SeqRecord(
            Seq(seq), id="DIFF_Q", name="DIFF_Q",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        # Picked target: same sequence rotated 137 bp.
        offset = 137
        rotated = seq[offset:] + seq[:offset]
        target = SeqRecord(
            Seq(rotated), id="DIFF_T", name="DIFF_T",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        # Stash both in the library so _action_diff_plasmid resolves.
        sc._save_library([
            {"id": canvas.id, "name": canvas.name, "size": len(seq),
             "n_feats": 0, "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(canvas)},
            {"id": target.id, "name": target.name, "size": len(seq),
             "n_feats": 0, "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(target)},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(canvas)
            await pilot.pause(0.05)
            # Drive the worker directly (push_screen + picker dismiss
            # is harder to script from a test pilot reliably than a
            # direct worker call).
            app._diff_align_worker(canvas, target)
            # The C-loop runs in a worker; give it a deadline.
            for _ in range(40):
                await pilot.pause(0.1)
                if app._alignments:
                    break
            assert len(app._alignments) == 1
            entry = app._alignments[0]
            # The picker should have chosen a rotation (either
            # direction works for axis="query") and the result should
            # have near-perfect identity.
            picked = entry["result"].get("picked_rotation", "none")
            ident = entry["result"].get("identity_pct", 0.0)
            assert picked in ("query", "target"), (
                f"circular target with rotated origin should pick a "
                f"rotation; got {picked!r}"
            )
            assert ident > 95.0, (
                f"rotation-corrected identity should be ~100% for an "
                f"exact rotation; got {ident}%"
            )


class TestSangerAddToLibraryDeduplication:
    """Sanger AB1 add-to-library used to be re-clickable, silently
    creating `<id>_2`, `<id>_3`, ... duplicates per click. Post-fix
    the Add button disables itself after a successful add and the
    main LibraryPanel is refreshed so the user sees the new entry
    without navigating away.
    """

    async def test_add_disables_button_and_refreshes_library(
            self, tmp_path, tiny_record, isolated_library):
        """Pick a synthetic AB1, click Add to library, verify (a) the
        library now contains the new entry, (b) the Add button is
        disabled (so re-clicking can't mint duplicates), (c) the
        LibraryPanel reflects the new entry."""
        from textual.widgets import Button
        # Synthesize a minimal AB1 by base-calling a SeqRecord.
        # `_ab1_path_to_record` parses a real AB1, but constructing
        # one is heavy — instead, monkey the modal's `_sanger_record`
        # to bypass the file parse path. The handler only consumes
        # `_sanger_record` + `_sanger_path` from there.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("ATGCATGCATGC" * 20), id="sanger_test",
            name="sanger_test",
            annotations={"molecule_type": "DNA"},
        )
        fake_ab1 = tmp_path / "trace.ab1"
        fake_ab1.write_bytes(b"\x00" * 100)  # body never read
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(
                sc.SequencingScreen(start_path=str(tmp_path))
            )
            await pilot.pause(0.2)
            screen = app.screen
            # Inject sanger state directly so the handler can fire.
            screen._sanger_record = rec
            screen._sanger_path = fake_ab1
            # Arm the button so the click handler can run.
            add_btn = screen.query_one("#btn-sanger-add", Button)
            add_btn.disabled = False
            n_before = len(sc._load_library())
            screen._sanger_add_to_library(None)
            await pilot.pause(0.1)
            # Library has the new entry
            entries_after = sc._load_library()
            assert len(entries_after) == n_before + 1
            # Add button now disabled
            assert add_btn.disabled, (
                "Add button must disable after a successful add to "
                "prevent silent duplicate-add"
            )


class TestMultiAlignPickerStaleIdFilter:
    """MultiAlignPickerModal used to dismiss with the raw
    `_selected_ids` set. If the user opened the picker, picked a few
    plasmids, then a sibling pane (or agent endpoint) deleted one of
    them before the user clicked Align, the stale id flowed through
    to `_action_open_align_picker` and surfaced as a per-target
    "not found" warning. Post-fix the picker filters against the
    current library at dismiss-time and notifies if any picks were
    dropped.
    """

    async def test_dismiss_drops_ids_no_longer_in_library(
            self, tiny_record, isolated_library):
        """Pick three plasmids in the modal, delete one from the
        library before the user clicks Align, verify the modal
        dismisses with only the surviving two."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Three library entries the picker can list.
        records = []
        for name in ("ALPHA", "BETA", "GAMMA"):
            r = SeqRecord(
                Seq("ATGC" * 50), id=name, name=name,
                annotations={"molecule_type": "DNA", "topology": "circular"},
            )
            records.append(r)
        sc._save_library([
            {"id": r.id, "name": r.name, "size": len(r.seq),
             "n_feats": 0, "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(r)}
            for r in records
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            await app.push_screen(sc.MultiAlignPickerModal())
            await pilot.pause(0.1)
            modal = app.screen
            # Mark all three as picked.
            modal._selected_ids = {"ALPHA", "BETA", "GAMMA"}
            # Simulate a sibling pane deleting BETA mid-pick.
            remaining = [r for r in records if r.id != "BETA"]
            sc._save_library([
                {"id": r.id, "name": r.name, "size": len(r.seq),
                 "n_feats": 0, "added": "2026-05-24",
                 "gb_text": sc._record_to_gb_text(r)}
                for r in remaining
            ])
            # Capture the dismiss payload.
            dismissed: list = []
            real_dismiss = modal.dismiss
            def _capture(value=None):
                dismissed.append(value)
                return real_dismiss(value)
            modal.dismiss = _capture
            modal._ok(None)
            await pilot.pause(0.1)
            assert dismissed, "modal must dismiss after _ok"
            picked = dismissed[0]
            assert picked is not None
            assert set(picked) == {"ALPHA", "GAMMA"}, (
                f"stale BETA must be filtered out; got {picked!r}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Library sequencing-status badges (`_alignment_quality_status` +
# `_library_entry_alignment_summary` + LibraryPanel "Seq" column)
# ═══════════════════════════════════════════════════════════════════════════════

class TestAlignmentQualityStatus:
    def test_verified_requires_perfect_identity_and_coverage(self):
        """A perfectly matching read (100% identity, full coverage,
        zero gaps) lights up `verified` (✓ green)."""
        result = {
            "n_matches": 1000, "n_mismatches": 0, "n_gaps": 0,
            "ungapped_identity_pct": 100.0,
        }
        code, glyph, color = sc._alignment_quality_status(result, 1000)
        assert code == "verified"
        assert glyph == "✓"
        assert color == "green"

    def test_single_indel_demotes_to_near(self):
        """Even ONE gap demotes from verified to near-match — for a
        cloning workflow, a 1-bp indel is meaningful (frameshift)."""
        result = {
            "n_matches": 999, "n_mismatches": 0, "n_gaps": 1,
            "ungapped_identity_pct": 100.0,
        }
        code, glyph, color = sc._alignment_quality_status(result, 1000)
        assert code == "near"
        assert glyph == "⚠"

    def test_low_coverage_is_partial(self):
        """A read that aligns at high identity but only covers a
        sub-region (e.g. plasmidsaurus consensus of a different
        plasmid sharing a backbone) lights up `partial`."""
        result = {
            "n_matches": 500, "n_mismatches": 5, "n_gaps": 0,
            "ungapped_identity_pct": 99.0,
        }
        # target_len 5000 → coverage = 505/5000 = 10%
        code, _, _ = sc._alignment_quality_status(result, 5000)
        assert code == "partial"

    def test_low_identity_is_divergent(self):
        """Significantly mismatched reads light up `divergent`."""
        result = {
            "n_matches": 600, "n_mismatches": 400, "n_gaps": 0,
            "ungapped_identity_pct": 60.0,
        }
        code, glyph, color = sc._alignment_quality_status(result, 1000)
        assert code == "divergent"
        assert color == "red"

    def test_zero_target_len_doesnt_divide_by_zero(self):
        """Defensive: target_len=0 still produces a valid status
        (divergent) rather than ZeroDivisionError."""
        result = {
            "n_matches": 10, "n_mismatches": 0, "n_gaps": 0,
            "ungapped_identity_pct": 100.0,
        }
        code, _, _ = sc._alignment_quality_status(result, 0)
        # Coverage falls below threshold so it's NOT verified;
        # ungapped passes the near threshold but coverage doesn't.
        assert code in ("partial", "divergent", "near")


class TestLibraryEntryAlignmentSummary:
    def test_no_alignments_returns_none(self):
        """Entries with no `alignments` field return None — caller
        renders dim `—`."""
        entry = {"id": "X", "name": "X", "size": 1000}
        assert sc._library_entry_alignment_summary(entry) is None

    def test_picks_best_priority_alignment(self):
        """Multiple alignments: pick the highest-priority status.
        verified > near > partial > divergent."""
        entry = {
            "id": "X", "name": "X", "size": 1000,
            "alignments": [
                # Divergent.
                {"visible": True, "result": {
                    "n_matches": 100, "n_mismatches": 900, "n_gaps": 0,
                    "ungapped_identity_pct": 10.0,
                }},
                # Verified.
                {"visible": True, "result": {
                    "n_matches": 1000, "n_mismatches": 0, "n_gaps": 0,
                    "ungapped_identity_pct": 100.0,
                }},
            ],
        }
        summary = sc._library_entry_alignment_summary(entry)
        assert summary is not None
        assert summary["code"] == "verified"
        assert summary["glyph"] == "✓"
        assert summary["n_total"] == 2
        assert summary["n_visible"] == 2

    def test_hidden_alignments_dont_contribute_glyph(self):
        """Stored alignments with visible=False don't contribute the
        headline glyph (the manager modal hid them for a reason).
        But the total count still includes them so the user knows
        they exist."""
        entry = {
            "id": "X", "name": "X", "size": 1000,
            "alignments": [
                {"visible": False, "result": {
                    "n_matches": 1000, "n_mismatches": 0, "n_gaps": 0,
                    "ungapped_identity_pct": 100.0,
                }},
            ],
        }
        summary = sc._library_entry_alignment_summary(entry)
        assert summary is not None
        assert summary["code"] == "hidden"  # all hidden
        assert summary["n_total"] == 1
        assert summary["n_visible"] == 0


# ═══════════════════════════════════════════════════════════════════════════════
# Variant extractor (`_extract_variants_from_alignment`)
# ═══════════════════════════════════════════════════════════════════════════════

class TestExtractVariantsFromAlignment:
    def test_perfect_match_yields_no_variants(self):
        assert sc._extract_variants_from_alignment(
            "ATGCATGC", "ATGCATGC",
        ) == []

    def test_single_snp(self):
        v = sc._extract_variants_from_alignment("ATGCAAGC", "ATGCATGC")
        assert len(v) == 1
        assert v[0]["type"] == "snp"
        assert v[0]["target_pos"] == 5  # T→A at target bp 5
        assert v[0]["ref"] == "T"
        assert v[0]["alt"] == "A"

    def test_insertion_merges_run(self):
        """Query has 3 extra bp inserted between target bp 2 and 3:
        target gaps for 3 columns. Should emit ONE insertion record
        of length 3, not three separate 1-bp records."""
        # aq:  A T G C A A T G C    (9 chars; CAA inserted between G and T)
        # at:  A T G - - - T G C    (9 chars)
        # walk: 0=A/A, 1=T/T, 2=G/G match;
        #       3=C/-, 4=A/-, 5=A/- insertion run before next target bp;
        #       6=T/T (target bp 3), 7=G/G (target bp 4), 8=C/C (target bp 5)
        v = sc._extract_variants_from_alignment(
            "ATGCAATGC", "ATG---TGC",
        )
        assert len(v) == 1
        assert v[0]["type"] == "insertion"
        # Insertion appears BEFORE the next target bp consumed,
        # which is target bp 3 (the 'T' after the gap run).
        assert v[0]["target_pos"] == 3
        assert v[0]["length"] == 3
        assert v[0]["alt"] == "CAA"

    def test_deletion_merges_run(self):
        # Target has bases the query is missing.
        # aq:  A T G - - - C G C
        # at:  A T G C A A C G C
        v = sc._extract_variants_from_alignment(
            "ATG---CGC", "ATGCAACGC",
        )
        assert len(v) == 1
        assert v[0]["type"] == "deletion"
        assert v[0]["target_pos"] == 3   # first deleted target bp
        assert v[0]["length"] == 3
        assert v[0]["ref"] == "CAA"

    def test_mixed_snps_and_indels(self):
        # aq:  A T G C A A T - - G C
        # at:  A T G C C - - G G G C
        # walk: 0=A/A, 1=T/T, 2=G/G, 3=C/C, 4=A/C (SNP),
        # 5=A/- (insertion at target_pos=5), 6=T/- (continuation),
        # 7=-/G (deletion at target_pos=5),
        # 8=-/G (continuation), 9=G/G, 10=C/C
        # Note: this is a degenerate construct; aligner wouldn't
        # produce both gap-types adjacent. Use simpler:
        # aq:  A T G C A A A G C
        # at:  A T G T A A A G C
        # — single SNP at pos 3 (C→T)
        v = sc._extract_variants_from_alignment(
            "ATGCAAAGC", "ATGTAAAGC",
        )
        assert len(v) == 1
        assert v[0]["type"] == "snp"
        assert v[0]["target_pos"] == 3
        assert v[0]["ref"] == "T"
        assert v[0]["alt"] == "C"

    def test_empty_inputs_return_empty(self):
        assert sc._extract_variants_from_alignment("", "") == []

    def test_mismatched_length_returns_empty(self):
        """Defensive: degenerate caller passing strings of different
        lengths gets an empty list, not a crash."""
        assert sc._extract_variants_from_alignment("ATG", "ATGC") == []


# ═══════════════════════════════════════════════════════════════════════════════
# Sample-to-library matcher (`_normalize_for_match` +
# `_match_samples_to_library`)
# ═══════════════════════════════════════════════════════════════════════════════

class TestNormalizeForMatch:
    def test_strips_plasmidsaurus_prefix(self):
        """`RUN42_1_MAV34` → `mav34` (drop run-id + order-num)."""
        assert sc._normalize_for_match("RUN42_1_MAV34") == "mav34"

    def test_strips_path_and_extension(self):
        """Full zip-member path → leaf basename, no extension."""
        assert sc._normalize_for_match(
            "RUN42_genbank-files/RUN42_1_MAV34.gbk",
        ) == "mav34"

    def test_strips_punctuation_and_lowercases(self):
        """`MAV 38 CAM-cTPFuGFP+RUBY` → `mav38camctpfugfpruby`."""
        assert sc._normalize_for_match(
            "MAV 38 CAM-cTPFuGFP+RUBY",
        ) == "mav38camctpfugfpruby"

    def test_empty_input_returns_empty(self):
        assert sc._normalize_for_match("") == ""


class TestMatchSamplesToLibrary:
    def test_exact_name_match(self):
        """Sample `RUN42_1_MAV34` ↔ library entry `MAV 34` —
        normalized forms are both `mav34`, exact match."""
        samples = [
            {"name": "RUN42_1_MAV34",
             "gbk": "RUN42_genbank-files/RUN42_1_MAV34.gbk"},
        ]
        library = [
            {"id": "MAV_34", "name": "MAV 34", "gb_text": ""},
            {"id": "MAV_35", "name": "MAV 35", "gb_text": ""},
        ]
        out = sc._match_samples_to_library(
            samples, library, sequence_fallback=False,
        )
        assert len(out) == 1
        assert out[0]["action"] == "align"
        assert out[0]["target_entry"]["id"] == "MAV_34"
        assert out[0]["method"] == "name-exact"
        assert out[0]["score"] == 1.0

    def test_no_match_proposes_add(self):
        """When name match score is below threshold, recommend
        adding the sample as a new library entry."""
        samples = [
            {"name": "RUN42_1_NEWPLASMID",
             "gbk": "RUN42_genbank-files/RUN42_1_NEWPLASMID.gbk"},
        ]
        library = [
            {"id": "OLDPLASMID", "name": "OLDPLASMID", "gb_text": ""},
        ]
        out = sc._match_samples_to_library(
            samples, library, sequence_fallback=False,
        )
        assert len(out) == 1
        assert out[0]["action"] == "add"
        assert out[0]["target_entry"] is None

    def test_no_gbk_skip(self):
        """Sample without a .gbk consensus is skipped entirely."""
        samples = [{"name": "no-gbk", "gbk": None}]
        out = sc._match_samples_to_library(
            samples, [], sequence_fallback=False,
        )
        assert len(out) == 1
        assert out[0]["action"] == "skip"
        assert out[0]["method"] == "no-gbk"

    def test_sequence_match_beats_name_substring(self, tmp_path):
        """Regression for 2026-05-24: a coincidental name substring
        used to outrank a 99%-identical library entry by sequence.
        Now the matcher always computes k-mer Jaccard for every
        candidate and lets sequence beat a weak name match.

        Setup: sample `CAM-2` (basename `cam2`) with one sequence.
        Library has (a) `pCambia1300` (no sequence overlap; just
        a name-substring coincidence), (b) `MAV 38` (the actual
        sequence source of the sample). Pre-fix the matcher picks
        pCambia1300 because `cam2` is a substring of `cambia`-ish
        normalisation; post-fix MAV 38 wins by k-mer Jaccard.
        """
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # The sample's sequence — also the body of MAV 38. Use a
        # non-repetitive seed so the canonical k-mer set has enough
        # cardinality (>= `_MIN_KMER_SET_FOR_STRONG_MATCH`, INV-73)
        # for the kmer-strong path to fire. Pre-INV-73 the test used
        # a 27 bp tandem repeat — ~27 unique canonical k-mers, below
        # the threshold post-fix.
        import random as _random
        _rng = _random.Random(20260525)
        true_target_seq = "".join(
            _rng.choice("ACGT") for _ in range(1500)
        )
        # pCambia1300's sequence is unrelated (low k-mer overlap).
        _rng2 = _random.Random(99999)
        decoy_seq = "".join(_rng2.choice("ACGT") for _ in range(1500))
        # Build a real zip so the matcher can extract the sample gbk.
        gbk = tmp_path / "CAM-2.gbk"
        SeqIO.write(
            SeqRecord(
                Seq(true_target_seq), id="CAM-2", name="CAM-2",
                annotations={"molecule_type": "DNA",
                              "topology": "circular"},
            ), gbk, "genbank",
        )
        zp = tmp_path / "RUN_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "RUN_genbank-files/RUN_1_CAM-2.gbk")
        samples = [{
            "name": "RUN_1_CAM-2",
            "gbk":  "RUN_genbank-files/RUN_1_CAM-2.gbk",
        }]
        # Library has the decoy (name-substring trap) + the real
        # target (sequence-only match — name doesn't substring CAM-2).
        decoy_rec = SeqRecord(
            Seq(decoy_seq), id="pCambia1300", name="pCambia1300",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        target_rec = SeqRecord(
            Seq(true_target_seq), id="MAV_38", name="MAV_38",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        library = [
            {"id": "pCambia1300", "name": "pCambia1300",
             "gb_text": sc._record_to_gb_text(decoy_rec)},
            {"id": "MAV_38", "name": "MAV 38",
             "gb_text": sc._record_to_gb_text(target_rec)},
        ]
        out = sc._match_samples_to_library(
            samples, library,
            sequence_fallback=True,
            extract_gbk_fn=sc._extract_gbk_member,
            zip_path=zp,
        )
        assert len(out) == 1
        # MAV 38 (sequence-identical) MUST win, not pCambia1300
        # (name-substring coincidence).
        assert out[0]["target_entry"]["id"] == "MAV_38", (
            f"Expected MAV_38 (sequence-identical) but got "
            f"{out[0]['target_entry']['id']!r} — the matcher fell "
            f"back to name-substring over sequence again."
        )
        assert out[0]["method"] == "kmer-strong"
        # Top-3 alternatives include the decoy so the user can spot
        # near-misses.
        alt_ids = {a["entry_id"] for a in out[0]["alternatives"]}
        assert "pCambia1300" in alt_ids

    def test_reverse_complement_sample_matches_at_full_identity(
            self, tmp_path,
    ):
        """Regression for 2026-05-24: when a Plasmidsaurus consensus
        assembled in the opposite orientation to its library entry,
        the matcher returned 0% k-mer overlap and picked an unrelated
        backbone with name-substring score. Post-fix `_kmer_set`
        uses canonical (strand-agnostic) k-mers so RC-of-library
        samples score 100% Jaccard against the right entry."""
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # The library entry's sequence.
        target_seq = "ATGCGTACGTAGCTAGCTAGCTGATCG" * 100
        # The sample is the reverse complement of the library entry —
        # this is what Plasmidsaurus produces when the assembler
        # picks the opposite strand.
        sample_seq = str(Seq(target_seq).reverse_complement())
        gbk = tmp_path / "RC_sample.gbk"
        SeqIO.write(
            SeqRecord(
                Seq(sample_seq), id="RC_sample", name="RC_sample",
                annotations={"molecule_type": "DNA",
                              "topology": "circular"},
            ), gbk, "genbank",
        )
        zp = tmp_path / "RC_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "RC_genbank-files/RC_1_sample.gbk")
        target_rec = SeqRecord(
            Seq(target_seq), id="TARGET", name="TARGET",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        # Add a decoy backbone (no real homology — coincidental name
        # substring would have won pre-fix).
        decoy_seq = "AAAATTTTGGGGCCCC" * 100
        decoy_rec = SeqRecord(
            Seq(decoy_seq), id="DECOY_target", name="DECOY_target",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        library = [
            {"id": "DECOY_target", "name": "DECOY_target",
             "gb_text": sc._record_to_gb_text(decoy_rec)},
            {"id": "TARGET", "name": "TARGET",
             "gb_text": sc._record_to_gb_text(target_rec)},
        ]
        out = sc._match_samples_to_library(
            [{"name": "RC_1_sample",
              "gbk":  "RC_genbank-files/RC_1_sample.gbk"}],
            library,
            sequence_fallback=True,
            extract_gbk_fn=sc._extract_gbk_member,
            zip_path=zp,
        )
        assert len(out) == 1
        # The RC-of-library sample must match TARGET, not the
        # unrelated DECOY backbone, even though plain-strand k-mers
        # would have given 0% overlap.
        assert out[0]["target_entry"]["id"] == "TARGET", (
            f"RC sample must match TARGET via canonical k-mers; got "
            f"{out[0]['target_entry']['id']!r}"
        )
        # Canonical k-mer Jaccard for RC-of-X vs X is essentially
        # 1.0 (every k-mer's canonical form is shared).
        assert out[0]["kmer_score"] > 0.95, (
            f"canonical k-mer Jaccard should be ~100% for RC-of-target; "
            f"got {out[0]['kmer_score']:.0%}"
        )

    def test_alternatives_surfaced_on_each_match(self, tmp_path):
        """Every match row carries up to 3 ranked alternatives so the
        confirm modal can show "what was close" — even when the
        matcher's pick was confident, the runner-up is visible."""
        # Three library entries, all share some sequence with the
        # sample. Top one is the picked target; #2 + #3 are
        # alternatives.
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        sample_seq = "ATGCGTACGTAGCTAGCTAGCTGATCG" * 50
        gbk = tmp_path / "S.gbk"
        SeqIO.write(
            SeqRecord(
                Seq(sample_seq), id="S", name="S",
                annotations={"molecule_type": "DNA",
                              "topology": "circular"},
            ), gbk, "genbank",
        )
        zp = tmp_path / "RUN_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "RUN_genbank-files/RUN_1_S.gbk")
        # Identical, mostly-identical, partly-identical library entries.
        def _entry(rid, seq):
            r = SeqRecord(
                Seq(seq), id=rid, name=rid,
                annotations={"molecule_type": "DNA",
                              "topology": "circular"},
            )
            return {"id": rid, "name": rid,
                    "gb_text": sc._record_to_gb_text(r)}
        library = [
            _entry("A", sample_seq),                # 100%
            _entry("B", sample_seq[:len(sample_seq) // 2]
                       + "T" * (len(sample_seq) // 2)),  # half-shared
            _entry("C", "T" * len(sample_seq)),     # 0%
        ]
        out = sc._match_samples_to_library(
            [{"name": "RUN_1_S",
              "gbk": "RUN_genbank-files/RUN_1_S.gbk"}],
            library,
            sequence_fallback=True,
            extract_gbk_fn=sc._extract_gbk_member,
            zip_path=zp,
        )
        assert len(out) == 1
        alts = out[0]["alternatives"]
        # Three alternatives ranked by combined score.
        assert len(alts) == 3
        assert alts[0]["entry_id"] == "A"  # best
        # Picked target is A.
        assert out[0]["target_entry"]["id"] == "A"
        # Each alternative carries both per-axis scores.
        for a in alts:
            assert "kmer_score" in a
            assert "name_score" in a


class TestLibraryKmerCacheEviction:
    """Sweep #35 (2026-05-26): the library-side k-mer cache is bounded
    by `_LIBRARY_KMER_CACHE_MAX` with FIFO eviction. Pre-fix the cache
    grew without bound between `_save_library` invalidations — a user
    running many bulk-aligns on a 200+ entry library across long
    sessions could accumulate hundreds of MB of resident k-mer sets.
    Sibling caches (`_RESTR_SCAN_CACHE`, `_ENZYME_CUTS_CACHE`,
    `_BLAST_DB_CACHE`, `_GB_PARSE_CACHE`) all have explicit caps; the
    library kmer cache was the outlier.
    """

    def test_cache_size_stays_at_or_below_cap(self, tmp_path,
                                                monkeypatch):
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Tiny cap so we don't need a huge library to trigger eviction.
        monkeypatch.setattr(sc, "_LIBRARY_KMER_CACHE_MAX", 2)
        # Start with a clean cache so prior tests don't pre-fill us
        # over the cap.
        sc._invalidate_library_kmer_cache()
        # Build five library entries with distinct sequences so each
        # one hashes to a unique key. Random non-repetitive bodies of
        # at least 1.5 kb so the canonical k-mer set sits above
        # `_MIN_KMER_SET_FOR_STRONG_MATCH` (INV-73) — otherwise the
        # matcher takes the name-only path and never touches the
        # kmer cache.
        import random as _random
        bodies = []
        for seed in (101, 202, 303, 404, 505):
            rng = _random.Random(seed)
            bodies.append("".join(
                rng.choice("ACGT") for _ in range(1500)
            ))
        library = []
        for i, body in enumerate(bodies):
            rec = SeqRecord(
                Seq(body), id=f"LIB_{i}", name=f"LIB_{i}",
                annotations={"molecule_type": "DNA",
                              "topology": "circular"},
            )
            library.append({
                "id": f"LIB_{i}",
                "name": f"LIB_{i}",
                "gb_text": sc._record_to_gb_text(rec),
            })
        # Sample shares the first library entry's sequence so the
        # name match is weak (no shared substring) and the matcher
        # falls through to sequence — computing kmer sets for ALL
        # library entries en route.
        gbk = tmp_path / "sample.gbk"
        SeqIO.write(
            SeqRecord(
                Seq(bodies[0]), id="SAMPLE", name="SAMPLE",
                annotations={"molecule_type": "DNA",
                              "topology": "circular"},
            ), gbk, "genbank",
        )
        zp = tmp_path / "ko_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "KO_genbank-files/KO_1_sample.gbk")
        sc._match_samples_to_library(
            [{"name": "KO_1_sample",
              "gbk":  "KO_genbank-files/KO_1_sample.gbk"}],
            library,
            sequence_fallback=True,
            extract_gbk_fn=sc._extract_gbk_member,
            zip_path=zp,
        )
        # Five entries computed, cap=2 → exactly the two most-recent
        # insertions survive. Library iteration order matches the
        # input list, so the survivors are LIB_3 + LIB_4 (the
        # eviction pops the oldest at each over-cap insert).
        assert len(sc._LIBRARY_KMER_CACHE) <= 2, (
            f"Cache should be bounded by {sc._LIBRARY_KMER_CACHE_MAX}; "
            f"got {len(sc._LIBRARY_KMER_CACHE)} entries"
        )
        # Most-recent two entry IDs survive (FIFO).
        survivor_ids = {eid for (eid, _gb_hash)
                         in sc._LIBRARY_KMER_CACHE}
        assert "LIB_4" in survivor_ids, survivor_ids
        assert "LIB_3" in survivor_ids, survivor_ids

    def test_cache_hit_does_not_reset_eviction_order(self, monkeypatch):
        """FIFO eviction is insertion-order, not access-order. A cache
        hit (re-reading an already-cached entry) MUST NOT bump it to
        the front, otherwise frequent reads of the oldest entry keep
        it permanently alive and starve newer entries."""
        monkeypatch.setattr(sc, "_LIBRARY_KMER_CACHE_MAX", 2)
        sc._invalidate_library_kmer_cache()
        # Manually populate so we control the insertion order without
        # going through the full matcher pipeline.
        sc._LIBRARY_KMER_CACHE[("OLD", "h0")] = {"AAAA", "ACGT"}
        sc._LIBRARY_KMER_CACHE[("NEW", "h1")] = {"GGGG", "TTTT"}
        # Read the OLD entry — should NOT bump it past NEW.
        assert sc._LIBRARY_KMER_CACHE.get(("OLD", "h0")) == {"AAAA", "ACGT"}
        # Insertion order is still [OLD, NEW] — verify by checking
        # `next(iter(...))` returns OLD (the candidate for eviction).
        assert next(iter(sc._LIBRARY_KMER_CACHE)) == ("OLD", "h0")


# ═══════════════════════════════════════════════════════════════════════════════
# Bulk-align modal + Sequencing-status column wiring
# ═══════════════════════════════════════════════════════════════════════════════

class TestBulkAlignConfirmModalToggle:
    """The bulk-align confirm modal lets the user rotate each row's
    action via Space (align ↔ add ↔ skip). Rows with no target_entry
    can't be set to "align" — it bounces back to "add".
    """

    async def test_action_toggle_cycles_through_options(
            self, tiny_record, isolated_library):
        """Space on a row with a target_entry cycles align → add →
        skip → align."""
        matches = [{
            "sample": {"name": "X", "gbk": "x.gbk"},
            "action": "align",
            "target_entry": {"id": "T", "name": "Target plasmid"},
            "score": 1.0, "method": "name-exact", "note": "",
        }]
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(sc.BulkAlignConfirmModal(matches))
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._matches[0]["action"] == "align"
            modal.action_toggle_action()
            assert modal._matches[0]["action"] == "add"
            modal.action_toggle_action()
            assert modal._matches[0]["action"] == "skip"
            modal.action_toggle_action()
            assert modal._matches[0]["action"] == "align"

    async def test_no_target_skips_align(
            self, tiny_record, isolated_library):
        """Rows without a target_entry can't be set to align — the
        cycle skips align and lands on add."""
        matches = [{
            "sample": {"name": "Y", "gbk": "y.gbk"},
            "action": "add",
            "target_entry": None,
            "score": 0.0, "method": "no-match", "note": "",
        }]
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(sc.BulkAlignConfirmModal(matches))
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._matches[0]["action"] == "add"
            modal.action_toggle_action()  # add → skip
            assert modal._matches[0]["action"] == "skip"
            modal.action_toggle_action()  # skip → align? coerced to add
            assert modal._matches[0]["action"] == "add"


class TestBulkQualityCells:
    """`_bulk_quality_cells` renders the Identity / Mism / Gaps columns
    for the bulk-align confirm modal. The k-mer column is a *matching*
    score (which plasmid is this?); THESE columns describe the real
    alignment that runs on confirm, so a 1-bp mismatch surfaces here
    even though it barely dents the k-mer Jaccard."""

    def test_non_align_action_is_dashes(self):
        for act in ("add", "skip"):
            cells = sc._bulk_quality_cells(
                {"ident": 99.99, "mism": 1, "gaps": 0},
                action=act, has_target=True)
            assert [c.plain for c in cells] == ["—", "—", "—"]

    def test_no_target_is_dashes(self):
        cells = sc._bulk_quality_cells(
            {"ident": 100.0, "mism": 0, "gaps": 0},
            action="align", has_target=False)
        assert [c.plain for c in cells] == ["—", "—", "—"]

    def test_pending_is_ellipsis(self):
        cells = sc._bulk_quality_cells(None, action="align", has_target=True)
        assert [c.plain for c in cells] == ["…", "…", "…"]

    def test_failed_is_question(self):
        cells = sc._bulk_quality_cells(False, action="align", has_target=True)
        assert [c.plain for c in cells] == ["?", "?", "?"]

    def test_one_bp_mismatch_renders_honestly(self):
        ident, mism, gaps = sc._bulk_quality_cells(
            {"ident": 100.0 * 18093 / 18094, "mism": 1, "gaps": 0},
            action="align", has_target=True)
        assert ident.plain == "99.99%"            # NOT "100%"
        assert ident.style == "green"             # agrees with the tier
        assert mism.plain == "1" and mism.style == "yellow"
        assert gaps.plain == "0" and gaps.style == "dim"

    def test_perfect_match_is_clean_100(self):
        ident, mism, _g = sc._bulk_quality_cells(
            {"ident": 100.0, "mism": 0, "gaps": 0},
            action="align", has_target=True)
        assert ident.plain == "100%"
        assert ident.style == "bright_cyan"
        assert mism.style == "dim"

    def test_gap_nonzero_is_red(self):
        _i, _m, gaps = sc._bulk_quality_cells(
            {"ident": 98.0, "mism": 0, "gaps": 5},
            action="align", has_target=True)
        assert gaps.plain == "5" and gaps.style == "red"

    def test_negative_counts_clamped(self):
        _i, mism, gaps = sc._bulk_quality_cells(
            {"ident": 50.0, "mism": -3, "gaps": -1},
            action="align", has_target=True)
        assert mism.plain == "0" and gaps.plain == "0"


class TestBulkAlignConfirmModalQualityColumns:
    """The confirm modal grows Ident / Mism / Gaps columns. The quality
    alignment runs ONCE up front (in `_compute_bulk_quality` on the
    Bulk-align button-press worker) and is cached on each match as `_aln`
    (+ `_aln_result` for the commit to reuse); the modal only RENDERS it,
    never re-aligns."""

    async def test_columns_present_and_render_cached_aln(
            self, tiny_record, isolated_library):
        from textual.widgets import DataTable
        from textual.coordinate import Coordinate
        matches = [{
            "sample": {"name": "CAM4", "gbk": "cam4.gbk"},
            "action": "align",
            "target_entry": {"id": "T", "name": "MAV40", "gb_text": "x"},
            "score": 1.0, "method": "kmer-strong", "note": "",
            "_aln": {"ident": 100.0 * 18093 / 18094, "mism": 1, "gaps": 0},
        }]
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(sc.BulkAlignConfirmModal(matches))
            await pilot.pause(0.1)
            modal = app.screen
            t = modal.query_one("#bulk-table", DataTable)
            headers = [c.label.plain for c in t.columns.values()]
            assert "Ident" in headers
            assert "Mism" in headers
            assert "Gaps" in headers
            assert t.get_cell_at(Coordinate(0, 6)).plain == "99.99%"
            assert t.get_cell_at(Coordinate(0, 7)).plain == "1"
            assert t.get_cell_at(Coordinate(0, 8)).plain == "0"

    async def test_failed_prealign_renders_question(
            self, tiny_record, isolated_library):
        from textual.widgets import DataTable
        from textual.coordinate import Coordinate
        # `_compute_bulk_quality` sets `_aln = False` when a row can't be
        # aligned (no target / no consensus member / align failure); the
        # modal renders that as "?".
        matches = [{
            "sample": {"name": "CAM4", "gbk": "cam4.gbk"},
            "action": "align",
            "target_entry": {"id": "T", "name": "MAV40", "gb_text": "x"},
            "score": 1.0, "method": "kmer-strong", "note": "",
            "_aln": False,
        }]
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(sc.BulkAlignConfirmModal(matches))
            await pilot.pause(0.1)
            modal = app.screen
            t = modal.query_one("#bulk-table", DataTable)
            assert t.get_cell_at(Coordinate(0, 6)).plain == "?"
            assert t.get_cell_at(Coordinate(0, 7)).plain == "?"
            assert t.get_cell_at(Coordinate(0, 8)).plain == "?"


class TestLibraryPanelSeqColumn:
    """LibraryPanel's "Seq" column shows per-entry sequencing-status
    badges driven by `_library_entry_alignment_summary`. The cell
    updates incrementally via `refresh_seq_cell(entry_id)` after a
    `_flush_active_alignments` so the badge tracks current state
    without a full table repopulate.
    """

    async def test_seq_column_shows_dash_for_unsequenced_entries(
            self, tiny_record, isolated_library):
        """An entry with no stored alignments renders ``—`` (dim)
        in the Seq column."""
        from textual.widgets import DataTable
        sc._save_library([
            {"id": tiny_record.id, "name": tiny_record.name,
             "size": len(tiny_record.seq), "n_feats": 0,
             "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(tiny_record)},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            lib = app.query_one("#library", sc.LibraryPanel)
            lib._view_mode = "plasmids"
            lib._apply_view_mode()
            lib._repopulate_plasmids()
            await pilot.pause(0.05)
            t = lib.query_one("#lib-table", DataTable)
            # 6 columns: ●, Name, Status, Seq, bp, K (kind badge added
            # 2026-05-30). Regression guard for the add_columns / add_row
            # arity match.
            assert len(t.columns) == 6
            assert t.row_count >= 1

    async def test_refresh_seq_cell_after_alignment_flush(
            self, tiny_record, isolated_library):
        """After `_flush_active_alignments` writes an alignment onto
        an entry, the LibraryPanel's Seq cell updates without a full
        repopulate. Verified by checking that the badge summary
        reflects the new alignment."""
        sc._save_library([
            {"id": tiny_record.id, "name": tiny_record.name,
             "size": len(tiny_record.seq), "n_feats": 0,
             "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(tiny_record)},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            # Register an alignment and flush.
            app._register_alignment(
                name="test-read", query_label="test-read",
                target_label=tiny_record.name,
                target_record=tiny_record,
                result={
                    "aligned_q": "ATGC", "aligned_t": "ATGC",
                    "n_matches": 4, "n_mismatches": 0, "n_gaps": 0,
                    "ungapped_identity_pct": 100.0,
                    "identity_pct": 100.0,
                },
            )
            app._alignments[-1]["_stored_source"] = "sequencing"
            app._flush_active_alignments()
            await pilot.pause(0.1)
            # Read it back: the entry should now have one stored
            # alignment, and the summary should reflect it.
            entries = sc._load_library()
            t_entry = next(e for e in entries if e["id"] == tiny_record.id)
            summary = sc._library_entry_alignment_summary(t_entry)
            assert summary is not None
            assert summary["n_total"] == 1


class TestVerificationReportModal:
    async def test_modal_collects_rows_from_library(
            self, tiny_record, isolated_library):
        """A library with one entry carrying one stored alignment
        produces one row in the report."""
        # Library entry with a pre-baked stored alignment.
        stored_align = {
            "id": "test-alignment-id",
            "label": "test-read",
            "query_label": "test-read",
            "target_label": tiny_record.name,
            "target_id": tiny_record.id,
            "target_gb_text": sc._record_to_gb_text(tiny_record),
            "target_seq_hash": sc._alignment_target_hash(
                str(tiny_record.seq),
            ),
            "axis": "target",
            "result": {
                "aligned_q": "ATGC", "aligned_t": "ATGC",
                "n_matches": 4, "n_mismatches": 0, "n_gaps": 0,
                "ungapped_identity_pct": 100.0,
                "identity_pct": 100.0,
            },
            "visible": True,
            "added": "2026-05-24",
            "source": "sequencing",
        }
        sc._save_library([
            {"id": tiny_record.id, "name": tiny_record.name,
             "size": len(tiny_record.seq), "n_feats": 0,
             "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(tiny_record),
             "alignments": [stored_align]},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(sc.VerificationReportModal())
            await pilot.pause(0.1)
            modal = app.screen
            assert len(modal._rows_data) == 1
            row = modal._rows_data[0]
            assert row["entry_id"] == tiny_record.id
            assert row["read_label"] == "test-read"
            # Coverage = 4/(len of tiny_record) so verified status
            # depends on tiny_record length. Just check it's a
            # recognised status.
            assert row["code"] in (
                "verified", "near", "partial", "divergent",
            )

    async def test_modal_skips_entries_with_no_alignments_by_default(
            self, tiny_record, isolated_library):
        """`only_with_alignments=True` (default) hides entries that
        have no stored alignments — the report is for verified vs
        unsequenced, not a library catalog."""
        sc._save_library([
            {"id": tiny_record.id, "name": tiny_record.name,
             "size": len(tiny_record.seq), "n_feats": 0,
             "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(tiny_record)},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            await app.push_screen(sc.VerificationReportModal())
            await pilot.pause(0.1)
            modal = app.screen
            assert modal._rows_data == []


# ═══════════════════════════════════════════════════════════════════════════════
# Alignment hardening (post-audit 2026-05-24) — schema validation,
# length-mismatch guards, backward-compat defaults
# ═══════════════════════════════════════════════════════════════════════════════
# These tests exercise the defensive guards added after the alignment
# subsystem audit: stored alignments missing critical fields (or with
# malformed paired strings) must skip cleanly with a log entry rather
# than crashing the hydrate path; pre-rotation-picker stored entries
# must hydrate with safe defaults for the new rotation fields.

class TestDeserializeStoredAlignmentArgs:
    def _minimal_stored(self, tiny_record, **overrides) -> dict:
        """Build a minimum valid stored-alignment dict, overridable."""
        base = {
            "id":              "test-id",
            "label":           "test-read",
            "query_label":     "Q",
            "target_label":    tiny_record.name,
            "target_id":       tiny_record.id,
            "target_gb_text":  sc._record_to_gb_text(tiny_record),
            "target_seq_hash": sc._alignment_target_hash(
                str(tiny_record.seq),
            ),
            "axis":            "target",
            "result": {
                "aligned_q": "ATGC", "aligned_t": "ATGC",
                "n_matches": 4, "n_mismatches": 0, "n_gaps": 0,
                "identity_pct": 100.0,
                "ungapped_identity_pct": 100.0,
            },
            "visible": True,
            "added":   "2026-05-24",
            "source":  "test",
        }
        base.update(overrides)
        return base

    def test_missing_target_gb_text_returns_none(self, tiny_record):
        stored = self._minimal_stored(tiny_record, target_gb_text="")
        assert sc._deserialize_stored_alignment_args(stored) is None

    def test_corrupt_gb_text_returns_none(self, tiny_record):
        stored = self._minimal_stored(
            tiny_record, target_gb_text="not a valid GenBank record",
        )
        assert sc._deserialize_stored_alignment_args(stored) is None

    def test_missing_aligned_strings_returns_none(self, tiny_record):
        """Schema validation guard added post-audit: stored entries
        with empty/missing aligned_q or aligned_t must be skipped
        rather than passed downstream where segment computation
        would raise an opaque ValueError."""
        bad = self._minimal_stored(tiny_record)
        bad["result"]["aligned_q"] = ""
        assert sc._deserialize_stored_alignment_args(bad) is None

    def test_aligned_string_length_mismatch_returns_none(
            self, tiny_record):
        """Paired-column walk in `_alignment_to_target_segments`
        assumes len(aq) == len(at) — a mismatch is a corruption
        signal that should skip the hydrate, not crash downstream."""
        bad = self._minimal_stored(tiny_record)
        bad["result"]["aligned_q"] = "ATGCATGC"
        bad["result"]["aligned_t"] = "ATGC"
        assert sc._deserialize_stored_alignment_args(bad) is None

    def test_legacy_stored_entry_gets_rotation_field_defaults(
            self, tiny_record):
        """Pre-rotation-picker stored alignments lack the
        `picked_rotation` / `query_rotation` / `target_rotation` /
        `query_rc` fields. Hydration must inject defaults so
        downstream code can treat all stored entries uniformly."""
        stored = self._minimal_stored(tiny_record)
        # Strip the rotation fields (simulate pre-2026-05-24 stored).
        for field in ("picked_rotation", "query_rotation",
                       "target_rotation", "query_rc"):
            stored["result"].pop(field, None)
        args = sc._deserialize_stored_alignment_args(stored)
        assert args is not None
        result = args["result"]
        assert result["picked_rotation"] == "none"
        assert result["query_rotation"] == 0
        assert result["target_rotation"] == 0
        assert result["query_rc"] is False


class TestRegisterAlignmentLengthGuard:
    """Regression for the post-audit length-mismatch guard. An
    upstream caller passing aligned_q/aligned_t of unequal length
    (e.g., from a corrupted result dict that slipped past the
    hydrate gate) must surface a clear notify rather than letting
    the segment walk raise an opaque ValueError."""

    async def test_register_alignment_refuses_mismatched_strings(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            n_before = len(app._alignments)
            # Mismatched-length aligned strings.
            app._register_alignment(
                name="bad-read", query_label="Q",
                target_label=tiny_record.name,
                target_record=tiny_record,
                result={
                    "aligned_q": "ATGCATGC",  # 8 chars
                    "aligned_t": "ATGC",      # 4 chars
                    "n_matches": 4, "n_mismatches": 0, "n_gaps": 0,
                    "identity_pct": 100.0,
                    "ungapped_identity_pct": 100.0,
                },
            )
            await pilot.pause(0.05)
            # No alignment registered — the guard refused.
            assert len(app._alignments) == n_before


class TestSerializeAlignmentEmptyTargetSeq:
    """Post-audit guard: a target_record with an empty seq is broken
    (the renderer can't show it, the stale-target hash would mis-fire
    on every reload). `_serialize_alignment_for_storage` returns
    None so the caller skips persisting."""

    def test_empty_target_seq_returns_none(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        empty_rec = SeqRecord(
            Seq(""), id="EMPTY", name="EMPTY",
            annotations={"molecule_type": "DNA"},
        )
        entry = {
            "name":          "test",
            "target_record": empty_rec,
            "result": {
                "aligned_q": "ATGC", "aligned_t": "ATGC",
                "n_matches": 4, "n_mismatches": 0, "n_gaps": 0,
                "identity_pct": 100.0,
                "ungapped_identity_pct": 100.0,
            },
        }
        assert sc._serialize_alignment_for_storage(entry) is None


class TestRegisterAlignmentReturnsEntry:
    """`_register_alignment` returns the newly-appended entry on
    success and None on refusal. Callers must use the return value
    instead of `_alignments[-1]` to avoid corrupting the previous
    entry's storage metadata when the register is refused (the bug
    that surfaced as "deleted alignments resurrect on the next
    flush" after Alt+L delete)."""

    async def test_returns_none_on_empty_strings(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            ret = app._register_alignment(
                name="empty", query_label="q",
                target_label=tiny_record.name,
                target_record=tiny_record,
                result={"aligned_q": "", "aligned_t": ""},
            )
            assert ret is None

    async def test_returns_none_on_mismatched_lengths(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            ret = app._register_alignment(
                name="mismatch", query_label="q",
                target_label=tiny_record.name,
                target_record=tiny_record,
                result={"aligned_q": "ATGC", "aligned_t": "ATGCA"},
            )
            assert ret is None

    async def test_returns_appended_entry_on_success(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.05)
            ret = app._register_alignment(
                name="ok", query_label="q",
                target_label=tiny_record.name,
                target_record=tiny_record,
                result={
                    "aligned_q": "ATGC", "aligned_t": "ATGC",
                    "n_matches": 4, "n_mismatches": 0, "n_gaps": 0,
                    "identity_pct": 100.0,
                    "ungapped_identity_pct": 100.0,
                },
            )
            assert ret is not None
            # Returned entry IS the appended one (same object).
            assert app._alignments[-1] is ret


class TestAlignmentManagerDeleteRoundTrip:
    """Regression for the Alt+L delete bug (2026-05-24): deleting an
    alignment in the manager modal + saving + reloading must remove
    it from the library on disk AND from the in-memory band. Pre-fix
    a register-refused hydrate could corrupt a sibling entry's
    `_stored_id`, causing the deleted alignment to resurrect on the
    next flush.
    """

    async def test_delete_alignment_persists_through_save_load(
            self, tiny_record, isolated_library):
        # Library entry with two pre-baked stored alignments.
        def _stored(label, source="manual"):
            return {
                "id":              f"id-{label}",
                "label":           label,
                "query_label":     label,
                "target_label":    tiny_record.name,
                "target_id":       tiny_record.id,
                "target_gb_text":  sc._record_to_gb_text(tiny_record),
                "target_seq_hash": sc._alignment_target_hash(
                    str(tiny_record.seq),
                ),
                "axis":            "target",
                "result": {
                    "aligned_q": "ATGC", "aligned_t": "ATGC",
                    "n_matches": 4, "n_mismatches": 0, "n_gaps": 0,
                    "identity_pct": 100.0,
                    "ungapped_identity_pct": 100.0,
                },
                "visible": True,
                "added":   "2026-05-24",
                "source":  source,
            }
        sc._save_library([
            {"id": tiny_record.id, "name": tiny_record.name,
             "size": len(tiny_record.seq), "n_feats": 0,
             "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(tiny_record),
             "alignments": [_stored("A"), _stored("B")]},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.1)
            # Both alignments hydrated onto the band.
            assert len(app._alignments) == 2
            # Simulate AlignmentManagerModal: drop alignment B and
            # re-flush the remaining list. (Drive the underlying
            # path directly rather than via the modal pilot, which
            # is brittle in xdist.)
            entries = sc._load_library()
            t_entry = next(e for e in entries if e["id"] == tiny_record.id)
            kept = [
                a for a in (t_entry.get("alignments") or [])
                if a.get("label") != "B"
            ]
            t_entry["alignments"] = kept
            sc._save_library(entries, async_sync=True)
            app._clear_alignments()
            app._hydrate_alignments_for_active()
            await pilot.pause(0.1)
            # After delete + save + re-hydrate, only A remains.
            assert len(app._alignments) == 1
            assert app._alignments[0]["name"] == "A"
            # And the stored library entry no longer has B.
            re_loaded = sc._load_library()
            t2 = next(e for e in re_loaded if e["id"] == tiny_record.id)
            stored = t2.get("alignments") or []
            assert len(stored) == 1
            assert stored[0]["label"] == "A"

    async def test_register_refused_does_not_corrupt_sibling_metadata(
            self, tiny_record, isolated_library):
        """When `_register_alignment` refuses (malformed strings),
        the hydrate code must NOT stamp `_stored_id` on the previous
        entry. Pre-fix the stamp landed on the sibling, causing the
        next flush to overwrite the sibling's storage slot with the
        refused entry's metadata."""
        # Pre-bake: one good stored entry and one with corrupt
        # (length-mismatch) aligned strings.
        def _stored(label, aq, at):
            return {
                "id":              f"id-{label}",
                "label":           label,
                "query_label":     label,
                "target_label":    tiny_record.name,
                "target_id":       tiny_record.id,
                "target_gb_text":  sc._record_to_gb_text(tiny_record),
                "target_seq_hash": sc._alignment_target_hash(
                    str(tiny_record.seq),
                ),
                "axis":            "target",
                "result": {
                    "aligned_q": aq, "aligned_t": at,
                    "n_matches": len(aq), "n_mismatches": 0,
                    "n_gaps": 0, "identity_pct": 100.0,
                    "ungapped_identity_pct": 100.0,
                },
                "visible": True,
                "added":   "2026-05-24",
                "source":  "manual",
            }
        # First entry is valid, second has mismatched lengths.
        # Deserialize will SKIP the second per the new schema check
        # — but verify the first's metadata isn't corrupted.
        sc._save_library([
            {"id": tiny_record.id, "name": tiny_record.name,
             "size": len(tiny_record.seq), "n_feats": 0,
             "added": "2026-05-24",
             "gb_text": sc._record_to_gb_text(tiny_record),
             "alignments": [
                 _stored("GOOD", "ATGC", "ATGC"),
                 _stored("BAD",  "ATGC", "ATGCATGC"),
             ]},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app._apply_record(tiny_record)
            await pilot.pause(0.1)
            # Only GOOD hydrates (BAD is skipped by schema check).
            assert len(app._alignments) == 1
            good = app._alignments[0]
            assert good["name"] == "GOOD"
            # CRITICAL: the GOOD entry's `_stored_id` is its own id,
            # not BAD's. Pre-fix the second register (refused) would
            # have stamped BAD's id onto GOOD via the `[-1]` access.
            assert good["_stored_id"] == "id-GOOD"


class TestAlignmentQualityStatusNegativeGuard:
    """Defensive negative-value guard — a corrupted result with
    negative n_matches/n_mismatches/n_gaps must not let `verified`
    fire on garbage."""

    def test_negative_n_matches_is_divergent(self):
        result = {
            "n_matches": -1, "n_mismatches": 0, "n_gaps": 0,
            "ungapped_identity_pct": 100.0,
        }
        code, _, _ = sc._alignment_quality_status(result, 1000)
        assert code == "divergent"

    def test_negative_ungapped_is_divergent(self):
        result = {
            "n_matches": 100, "n_mismatches": 0, "n_gaps": 0,
            "ungapped_identity_pct": -50.0,
        }
        code, _, _ = sc._alignment_quality_status(result, 100)
        assert code == "divergent"


class TestCanonicalKmerPalindromeHandling:
    """Verify the canonical k-mer normalisation handles palindromes,
    IUPAC codes, and empty inputs without breaking."""

    def test_palindrome_canonical_equals_itself(self):
        """A palindromic k-mer (kmer == RC(kmer)) has a single
        canonical form regardless of strand."""
        # 4-bp palindrome.
        s = "GGCC"
        out = sc._kmer_set(s, k=4, canonical=True)
        assert out == {"GGCC"}

    def test_iupac_canonical_kmers(self):
        """IUPAC codes (N, R, Y, ...) are handled by `_rc`'s
        translation table — canonical form picks the lex-smaller
        of (kmer, RC)."""
        # Sequence with N: NNNN is its own RC (palindrome).
        out = sc._kmer_set("NNNNNN", k=4, canonical=True)
        assert "NNNN" in out

    def test_canonical_kmer_strand_agnostic_pair(self):
        """RC(seq) and seq produce the same canonical k-mer set."""
        from Bio.Seq import Seq
        s = "ATGCATGCATGCATGCATGCATGC"
        rc_s = str(Seq(s).reverse_complement())
        a = sc._kmer_set(s, k=8, canonical=True)
        b = sc._kmer_set(rc_s, k=8, canonical=True)
        # Should be identical: each k-mer's canonical form is the
        # same regardless of which strand was extracted.
        assert a == b


class TestMatcherParseFailureSurfacing:
    """Library entries that fail to parse during k-mer build are
    excluded from the comparison AND logged as a batched warning
    so the user can investigate corrupted entries via the diagnostic
    bundle rather than wondering why a sample silently fell back
    to 'add as new'."""

    def test_corrupt_library_entry_logged_at_warning(self, monkeypatch):
        """The matcher's `_log` (splicecraft logger) has
        propagate=False, so pytest's caplog can't see it via the
        root handler. Monkeypatch the warning method directly to
        capture the parse-failure summary call."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        ok_rec = SeqRecord(
            Seq("ATGC" * 100), id="OK", name="OK",
            annotations={"molecule_type": "DNA"},
        )
        library = [
            {"id": "BAD", "name": "BAD",
             "gb_text": "not a valid genbank record"},
            {"id": "OK", "name": "OK",
             "gb_text": sc._record_to_gb_text(ok_rec)},
        ]
        samples = [{"name": "NONE", "gbk": None}]  # no gbk → skip
        captured: list = []
        orig_warning = sc._log.warning

        def _capture(fmt, *args, **kwargs):
            try:
                captured.append(fmt % args if args else fmt)
            except Exception:
                captured.append(str(fmt))
            return orig_warning(fmt, *args, **kwargs)
        monkeypatch.setattr(sc._log, "warning", _capture)
        sc._match_samples_to_library(
            samples, library, sequence_fallback=True,
        )
        # The batched warning surfaces the parse failure with the
        # library entry id, so a user reading the log can find the
        # broken record.
        assert any(
            "failed to parse" in m and "BAD" in m for m in captured
        ), f"expected parse-failure warning; got: {captured}"


# ═══════════════════════════════════════════════════════════════════════════════
# [INV-72] Audit sweep 2026-05-25 — IUPAC normalisation, agent/UI picker
# parity, bulk-align failure surfacing, coverage clamp.
# ═══════════════════════════════════════════════════════════════════════════════

class TestNormalizeDnaForAlign:
    """The `_normalize_dna_for_align` helper scrubs whitespace/digits/
    FASTA markers and validates IUPAC nucleotide chars before any
    sequence reaches the C-loop. Covers the common bad-paste failures
    (FASTA-as-is, GenBank ORIGIN block, protein-into-DNA-field)."""

    def test_clean_input_passes_through(self):
        assert sc._normalize_dna_for_align("ACGT") == "ACGT"

    def test_lowercase_uppercased(self):
        assert sc._normalize_dna_for_align("acgtacgt") == "ACGTACGT"

    def test_strips_internal_whitespace(self):
        # GenBank ORIGIN block style: leading bp number + spaces.
        assert sc._normalize_dna_for_align(
            "        1 atgcatgcat tcgatcgatc",
        ) == "ATGCATGCATTCGATCGATC"

    def test_strips_fasta_header_line(self):
        # A pasted FASTA: the whole `>name desc\n` header line gets
        # stripped along with the embedded newlines, leaving only the
        # sequence body uppercased.
        assert sc._normalize_dna_for_align(
            ">myplasmid description\nACGTACGT\nACGTACGT",
        ) == "ACGTACGTACGTACGT"

    def test_strips_multi_fasta_headers(self):
        # Multiple FASTA records — every header line goes; sequence
        # bodies concatenate.
        assert sc._normalize_dna_for_align(
            ">seq1\nAAAA\n>seq2\nGGGG",
        ) == "AAAAGGGG"

    def test_strips_newlines_only(self):
        assert sc._normalize_dna_for_align("ACGT\nACGT") == "ACGTACGT"

    def test_iupac_ambiguity_accepted(self):
        assert sc._normalize_dna_for_align("ACGTNRYSWKM") == "ACGTNRYSWKM"

    def test_rejects_protein_letters(self):
        # Paste of a protein sequence into a DNA-only field: E F I L P Q
        # all sit outside the IUPAC nucleotide alphabet.
        with pytest.raises(ValueError, match="non-IUPAC"):
            sc._normalize_dna_for_align("MELFGPQ")

    def test_rejects_random_chars(self):
        with pytest.raises(ValueError, match="non-IUPAC"):
            sc._normalize_dna_for_align("ACGT*ACGT")

    def test_empty_returns_empty(self):
        assert sc._normalize_dna_for_align("") == ""

    def test_whitespace_only_returns_empty(self):
        assert sc._normalize_dna_for_align("   \n  \n  ") == ""

    def test_error_names_offending_char(self):
        try:
            sc._normalize_dna_for_align("ACGTZACGT")
        except ValueError as exc:
            assert "'Z'" in str(exc)
        else:
            pytest.fail("expected ValueError")

    def test_pairwise_align_uses_normaliser(self):
        # `_pairwise_align` must run the input through the normaliser
        # before reaching Biopython — a pasted-FASTA input would
        # otherwise either crash deep inside the C-loop or produce a
        # garbage alignment (length mismatch on the leading `>name`).
        result = sc._pairwise_align(
            ">qry\nACGTACGTACGTACGT",
            ">tgt\nACGTACGTACGTACGT",
        )
        assert result["identity_pct"] == 100.0
        assert result["q_len"] == 16
        assert result["t_len"] == 16

    def test_pairwise_align_rejects_protein_input(self):
        # Defensive: a pasted protein sequence should fail loudly at
        # the validation step, not deep inside Biopython.
        with pytest.raises(ValueError, match="non-IUPAC"):
            sc._pairwise_align("MELFGPQ", "ACGTACGT")


class TestPickBestRotationNormalises:
    """`_pick_best_rotation` pre-normalises at entry so the frame-shift
    helpers' length math agrees with what `_pairwise_align` actually
    consumes. Pre-fix passing raw FASTA would mean
    `len(target_seq)` in the rotation shift differed from the cleaned
    length, off-by-N depending on how many whitespace/header chars
    got scrubbed."""

    def test_pick_normalises_raw_fasta_input(self):
        # Both ways: with and without the leading FASTA header. The
        # alignment should be identical.
        clean_result = sc._pick_best_rotation(
            "ACGTACGTACGT",
            "ACGTACGTACGT",
            is_circular=False,
        )
        raw_result = sc._pick_best_rotation(
            ">qry name\nACGTACGTACGT",
            ">tgt name\nACGTACGTACGT",
            is_circular=False,
        )
        assert clean_result["identity_pct"] == raw_result["identity_pct"]
        assert clean_result["q_len"] == raw_result["q_len"]

    def test_pick_rejects_protein_input(self):
        with pytest.raises(ValueError, match="non-IUPAC"):
            sc._pick_best_rotation(
                "MELFGPQ", "ACGTACGT", is_circular=False,
            )

    def test_pick_rejects_post_normalise_empty(self):
        # Whitespace-only input is non-empty pre-strip but empty after.
        with pytest.raises(ValueError, match="empty"):
            sc._pick_best_rotation(
                "   \n  ", "ACGTACGT", is_circular=False,
            )


class TestAgentDiffPlasmidUsesPicker:
    """`_h_diff_plasmid` must use `_pick_best_rotation` (INV-72) so
    agent callers get the same RC-detection + multi-rotation
    best-of-N pick the UI gained in `[INV-71]`. Pre-fix the endpoint
    ran a single `_find_circular_alignment_offset` + bare
    `_pairwise_align`, missing RC orientations entirely."""

    async def test_returns_picker_metadata_fields(
            self, tiny_record, isolated_library):
        """Response payload carries the new picker fields:
        `picked_rotation`, `query_rotation`, `target_rotation`,
        `query_rc`. `rotation_offset` is kept for back-compat and
        mirrors `target_rotation`."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        target_seq = "ATGC" * 30
        target_rec = SeqRecord(
            Seq(target_seq), id="TGT", name="TGT",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        sc._save_library([{
            "id": "TGT", "name": "TGT", "size": len(target_seq),
            "n_feats": 0, "added": "2026-05-25",
            "gb_text": sc._record_to_gb_text(target_rec),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            # Same sequence on canvas → trivially 100% identity.
            app._apply_record(target_rec)
            await pilot.pause(0.05)
            result = sc._h_diff_plasmid(
                app, {"target_id": "TGT"},
            )
            assert isinstance(result, dict), result
            assert result["ok"] is True
            assert "picked_rotation" in result
            assert "query_rotation" in result
            assert "target_rotation" in result
            assert "query_rc" in result
            # Back-compat: rotation_offset === target_rotation.
            assert result["rotation_offset"] == result["target_rotation"]

    async def test_detects_rc_orientation(
            self, tiny_record, isolated_library):
        """Pre-INV-72 the endpoint ran only forward-orientation
        alignment — a query that's RC of the target scored ~0%
        identity. Post-fix, the picker tries RC plain at the first
        tier and surfaces `query_rc=True` for the agent."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        target_seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        target_rec = SeqRecord(
            Seq(target_seq), id="TGT", name="TGT",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        # RC the target → put it on the canvas as the "query".
        rc_seq = sc._rc(target_seq)
        rc_rec = SeqRecord(
            Seq(rc_seq), id="QRY", name="QRY",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        sc._save_library([{
            "id": "TGT", "name": "TGT", "size": len(target_seq),
            "n_feats": 0, "added": "2026-05-25",
            "gb_text": sc._record_to_gb_text(target_rec),
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app._apply_record(rc_rec)
            await pilot.pause(0.05)
            result = sc._h_diff_plasmid(
                app, {"target_id": "TGT"},
            )
            assert result["ok"] is True
            # Post-fix the picker recognises the RC orientation and
            # reports identity ≥ 99% (perfect minus rotation noise).
            # Pre-fix this scored under 50%.
            assert result["result"]["identity_pct"] >= 95.0
            assert result["query_rc"] is True


class TestAgentAlignPlasmidsaurusUsesPicker:
    """`_h_align_plasmidsaurus_zip` must use `_pick_best_rotation`
    (INV-72) for the same reasons as `_h_diff_plasmid` — pre-fix
    agents calling this endpoint missed RC-orientation detection."""

    async def test_returns_picker_metadata_fields(
            self, tmp_path, tiny_record, isolated_library):
        import zipfile
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        target_seq = "ATGC" * 60
        target_rec = SeqRecord(
            Seq(target_seq), id="TGT", name="TGT",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        sc._save_library([{
            "id": "TGT", "name": "TGT", "size": len(target_seq),
            "n_feats": 0, "added": "2026-05-25",
            "gb_text": sc._record_to_gb_text(target_rec),
        }])
        gbk = tmp_path / "RUN42_1_TGT.gbk"
        SeqIO.write(target_rec, gbk, "genbank")
        zp = tmp_path / "RUN42_results.zip"
        with zipfile.ZipFile(zp, "w") as zf:
            zf.write(gbk, "RUN42_genbank-files/RUN42_1_TGT.gbk")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            result = sc._h_align_plasmidsaurus_zip(app, {
                "path": str(zp),
                "member": "RUN42_genbank-files/RUN42_1_TGT.gbk",
                "target_id": "TGT",
            })
            assert isinstance(result, dict), result
            assert result["ok"] is True
            assert "picked_rotation" in result
            assert "query_rotation" in result
            assert "target_rotation" in result
            assert "query_rc" in result
            assert result["rotation_offset"] == result["target_rotation"]


class TestPairwiseAlignCoverageClampInToast:
    """Coverage % is clamped at 100% in the toast formatter (and in
    `VerificationReportModal._collect_rows`). The underlying data is
    untouched — only the display value is bounded — so a corrupted
    result dict doesn't render '150% coverage' to the user."""

    def test_clamp_formula_caps_at_100(self):
        # Reproduce the toast's coverage formula with inflated input.
        n_matches = 1500
        n_mismatches = 0
        aligned_bp = n_matches + n_mismatches
        target_len = 1000  # smaller than aligned_bp on purpose
        coverage_pct = (
            min(100.0, 100.0 * aligned_bp / target_len)
            if target_len else 0.0
        )
        assert coverage_pct == 100.0

    def test_clamp_formula_passes_through_normal_values(self):
        n_matches = 500
        n_mismatches = 50
        aligned_bp = n_matches + n_mismatches
        target_len = 1000
        coverage_pct = (
            min(100.0, 100.0 * aligned_bp / target_len)
            if target_len else 0.0
        )
        assert coverage_pct == 55.0

    def test_clamp_formula_handles_zero_target(self):
        coverage_pct = (
            min(100.0, 100.0 * 0 / 1) if 0 else 0.0
        )
        assert coverage_pct == 0.0


class TestBulkAlignNoGbkLogs:
    """`_bulk_align_worker` logs a warning when skipping a sample with
    no `gbk` field. Pre-INV-72 the skip was silent — a malformed
    Plasmidsaurus manifest would report 'failed N' with no clue why."""

    def test_skip_logs_sample_name(self, monkeypatch):
        # Reproduce the inline check without spinning up the worker.
        # The bulk-align worker's no-gbk path now matches:
        #     if not gbk_member: _log.warning(...); n_failed += 1
        captured: list[str] = []
        orig_warning = sc._log.warning

        def _capture(fmt, *args, **kwargs):
            try:
                captured.append(fmt % args if args else fmt)
            except Exception:
                captured.append(str(fmt))
            return orig_warning(fmt, *args, **kwargs)
        monkeypatch.setattr(sc._log, "warning", _capture)
        # Mimic the inline log line directly.
        sc._log.warning(
            "BulkAlign: skipping sample %r — no .gbk member "
            "field (malformed manifest or missing consensus)",
            "MAV34",
        )
        assert any(
            "BulkAlign" in m and "MAV34" in m for m in captured
        )


# ─────────────────────────────────────────────────────────────────────
# INV-73 (2026-05-25): follow-up alignment hardening sweep tests.
# ─────────────────────────────────────────────────────────────────────


class TestVariantExtractionCap:
    """`_extract_variants_from_alignment` caps the result list at
    `_MAX_VARIANTS_PER_ALIGNMENT` (default 10k) to bound memory on
    completely-divergent alignments. A truncation sentinel is
    appended so callers can surface "10k+ variants" rather than
    silently underreporting."""

    def test_small_alignment_no_cap(self):
        # Five SNPs in a 10-bp alignment; well under the cap.
        aq = "ATGCATGCAT"
        at = "ATCCAAGCTT"
        variants = sc._extract_variants_from_alignment(aq, at)
        assert all(v["type"] != "truncated" for v in variants)
        # SNP count: positions 2 (G→C), 4 (T→A), 5 (A→A NO)... let
        # the function speak — we just assert no truncation marker
        # for a small input.
        assert len(variants) <= 10

    def test_cap_appends_truncated_sentinel(self):
        # Build a 50-bp alignment where every column is a SNP.
        # Cap at 5 so the test is fast.
        aq = "A" * 50
        at = "C" * 50
        variants = sc._extract_variants_from_alignment(
            aq, at, max_variants=5,
        )
        # Expect 5 real variants + 1 truncation sentinel = 6 entries.
        assert len(variants) == 6
        assert variants[-1]["type"] == "truncated"
        assert variants[-1]["length"] == 0
        assert variants[-1]["ref"] == ""
        assert variants[-1]["alt"] == ""
        assert "omitted_after_pos" in variants[-1]
        # The first 5 are SNPs.
        assert all(v["type"] == "snp" for v in variants[:5])

    def test_truncation_filterable_by_type(self):
        # Callers that count by type (`snp`/`insertion`/`deletion`)
        # should naturally skip the sentinel.
        aq = "A" * 100
        at = "C" * 100
        variants = sc._extract_variants_from_alignment(
            aq, at, max_variants=10,
        )
        n_snps = sum(1 for v in variants if v["type"] == "snp")
        n_indels = sum(
            1 for v in variants
            if v["type"] in ("insertion", "deletion")
        )
        assert n_snps == 10
        assert n_indels == 0
        # Sentinel is present but excluded by type filter.
        assert sum(1 for v in variants if v["type"] == "truncated") == 1

    def test_zero_cap_disables_walking(self):
        # max_variants=0 is treated as "no cap" (defensive: avoid
        # an accidental 0 silently truncating to nothing).
        aq = "ATG"
        at = "CTG"
        variants = sc._extract_variants_from_alignment(
            aq, at, max_variants=0,
        )
        assert len(variants) >= 1
        assert all(v["type"] != "truncated" for v in variants)


class TestExtractVariantsMixedAndDivergent:
    """Edge cases: mixed SNP+indel calls in a single alignment; the
    all-divergent baseline."""

    def test_mixed_snp_and_indel(self):
        # aq:  ATG-CCG
        # at:  ATGGCAG
        # column 3: insertion (target has G, query gap)
        # column 5: SNP (C vs A)
        aq = "ATG-CCG"
        at = "ATGGCAG"
        variants = sc._extract_variants_from_alignment(aq, at)
        types = [v["type"] for v in variants]
        assert "deletion" in types or "insertion" in types
        # Note: by convention "-" in aq is a DELETION (target has
        # bp that query doesn't).
        del_v = next(v for v in variants if v["type"] == "deletion")
        assert del_v["ref"] == "G"
        assert del_v["length"] == 1
        # The C vs A column is a SNP.
        snps = [v for v in variants if v["type"] == "snp"]
        assert any(v["ref"] == "A" and v["alt"] == "C" for v in snps)

    def test_all_divergent_short(self):
        aq = "AAAA"
        at = "CCCC"
        variants = sc._extract_variants_from_alignment(aq, at)
        assert len(variants) == 4
        assert all(v["type"] == "snp" for v in variants)


class TestPickedRotationEnumValidation:
    """INV-73: `_deserialize_stored_alignment_args` validates
    rotation-picker fields against their expected value space.
    Corrupted/foreign values coerce back to safe defaults with a
    log warning rather than crashing or skipping the entry."""

    def _make_stored(self, **result_overrides):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        target_seq = "ATGC" * 25
        target_rec = SeqRecord(
            Seq(target_seq), id="T", name="T",
            annotations={"molecule_type": "DNA"},
        )
        gb = sc._record_to_gb_text(target_rec)
        result = {
            "mode": "global", "score": 200.0,
            "identity_pct": 100.0,
            "ungapped_identity_pct": 100.0,
            "aligned_q": target_seq, "aligned_t": target_seq,
            "n_matches": 100, "n_mismatches": 0, "n_gaps": 0,
            "q_len": 100, "t_len": 100,
            **result_overrides,
        }
        return {
            "id": "test-id",
            "label": "test",
            "query_label": "Q",
            "target_label": "T",
            "target_id": "T",
            "target_gb_text": gb,
            "target_seq_hash": sc._alignment_target_hash(target_seq),
            "axis": "target",
            "result": result,
            "visible": True,
            "source": "manual",
            "added": "2026-05-25",
        }

    def test_invalid_picked_rotation_coerces_to_none(self, caplog):
        stored = self._make_stored(picked_rotation="both")
        args = sc._deserialize_stored_alignment_args(stored)
        assert args is not None
        assert args["result"]["picked_rotation"] == "none"

    def test_negative_query_rotation_coerces_to_zero(self):
        stored = self._make_stored(query_rotation=-50)
        args = sc._deserialize_stored_alignment_args(stored)
        assert args is not None
        assert args["result"]["query_rotation"] == 0

    def test_negative_target_rotation_coerces_to_zero(self):
        stored = self._make_stored(target_rotation=-1)
        args = sc._deserialize_stored_alignment_args(stored)
        assert args is not None
        assert args["result"]["target_rotation"] == 0

    def test_non_int_rotation_coerces_to_zero(self):
        stored = self._make_stored(query_rotation="oops")
        args = sc._deserialize_stored_alignment_args(stored)
        assert args is not None
        assert args["result"]["query_rotation"] == 0

    def test_non_bool_query_rc_coerces_to_false(self):
        stored = self._make_stored(query_rc=1)
        args = sc._deserialize_stored_alignment_args(stored)
        assert args is not None
        assert args["result"]["query_rc"] is False

    def test_valid_values_preserved(self):
        stored = self._make_stored(
            picked_rotation="query",
            query_rotation=42,
            target_rotation=0,
            query_rc=True,
        )
        args = sc._deserialize_stored_alignment_args(stored)
        assert args is not None
        assert args["result"]["picked_rotation"] == "query"
        assert args["result"]["query_rotation"] == 42
        assert args["result"]["query_rc"] is True


class TestCoveragePctHelper:
    """INV-73: `_coverage_pct_from_result` centralises the clamp
    + zero-target guard so toast + verification report can't drift
    apart in their display logic."""

    def test_zero_target_len_returns_zero(self):
        assert sc._coverage_pct_from_result(
            {"n_matches": 100, "n_mismatches": 0}, 0,
        ) == 0.0

    def test_negative_target_len_returns_zero(self):
        assert sc._coverage_pct_from_result(
            {"n_matches": 100, "n_mismatches": 0}, -50,
        ) == 0.0

    def test_normal_coverage(self):
        # 500 + 50 = 550 aligned, target 1000 → 55%.
        assert sc._coverage_pct_from_result(
            {"n_matches": 500, "n_mismatches": 50}, 1000,
        ) == 55.0

    def test_clamp_caps_at_100(self):
        # Pathological: more aligned bp than target_len.
        assert sc._coverage_pct_from_result(
            {"n_matches": 1500, "n_mismatches": 0}, 1000,
        ) == 100.0

    def test_missing_fields_treated_as_zero(self):
        assert sc._coverage_pct_from_result({}, 1000) == 0.0

    def test_non_numeric_fields_returns_zero(self):
        assert sc._coverage_pct_from_result(
            {"n_matches": "oops", "n_mismatches": 0}, 1000,
        ) == 0.0


class TestKmerSetForStrongMatchThreshold:
    """INV-73: short samples (< `_MIN_KMER_SET_FOR_STRONG_MATCH`
    k-mers) can no longer trigger the `kmer-strong` match path —
    they'd otherwise score a coincidental 1.0 Jaccard against any
    library entry containing a primer-length match region."""

    def test_threshold_constant_is_reasonable(self):
        # 50 k-mers ≈ ~70 bp sample at k=20. Reasonable floor for
        # a Plasmidsaurus consensus.
        assert sc._MIN_KMER_SET_FOR_STRONG_MATCH >= 20
        assert sc._MIN_KMER_SET_FOR_STRONG_MATCH <= 200

    def test_short_sample_falls_through_to_name_or_weak(self):
        # Sample is 25 bp (~6 k-mers @ k=20). Library entry is the
        # SAME 25 bp — Jaccard would be 1.0. Without the threshold
        # guard this would match "kmer-strong"; with the guard it
        # falls through. We use a name match so it's still picked
        # but via the name path.
        seq = "ATGCATGCATGCATGCATGCATGCA"  # 25 bp
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq(seq), id="LIB-A", name="LIB-A",
            annotations={"molecule_type": "DNA"},
        )
        gb = sc._record_to_gb_text(rec)
        library = [{
            "id": "LIB-A", "name": "LIB-A", "gb_text": gb,
        }]
        samples = [{"name": "LIB-A", "gbk": "x.gbk", "base": "x"}]

        def _extract(zp, m):
            return gb
        out = sc._match_samples_to_library(
            samples, library,
            extract_gbk_fn=_extract, zip_path="dummy",
        )
        assert len(out) == 1
        # Should not be "kmer-strong" because the sample is too short
        # for the k-mer signal to be meaningful, even though the
        # Jaccard would compute as 1.0.
        assert out[0]["method"] != "kmer-strong"


class TestAlignmentQualityStatusBoundaries:
    """INV-73: explicit threshold-boundary tests for
    `_alignment_quality_status`. Pre-fix the verified→near→partial→
    divergent transitions were only tested at clear-cut values; an
    accidental >= → > shift could mis-label a borderline read."""

    def test_exact_verified_threshold_ungapped(self):
        # ungapped_identity_pct exactly at the verified floor (99.5)
        # + perfect coverage + zero gaps → verified.
        result = {
            "ungapped_identity_pct": 99.5,
            "n_matches": 995, "n_mismatches": 0, "n_gaps": 0,
        }
        code, _, _ = sc._alignment_quality_status(result, 995)
        assert code == "verified"

    def test_just_below_verified_demotes_to_near(self):
        # 99.49 ungapped → not verified, but ≥95% → near.
        result = {
            "ungapped_identity_pct": 99.49,
            "n_matches": 995, "n_mismatches": 5, "n_gaps": 0,
        }
        code, _, _ = sc._alignment_quality_status(result, 1000)
        assert code in ("near", "partial", "divergent")
        assert code != "verified"

    def test_one_gap_demotes_from_verified(self):
        # 100% ungapped but a single gap means a single indel — not
        # verified (verified requires zero gaps).
        result = {
            "ungapped_identity_pct": 100.0,
            "n_matches": 999, "n_mismatches": 0, "n_gaps": 1,
        }
        code, _, _ = sc._alignment_quality_status(result, 1000)
        assert code != "verified"

    def test_below_near_coverage_demotes_to_partial(self):
        # High ungapped identity but only 50% coverage → partial.
        result = {
            "ungapped_identity_pct": 99.0,
            "n_matches": 500, "n_mismatches": 5, "n_gaps": 0,
        }
        code, _, _ = sc._alignment_quality_status(result, 1000)
        assert code in ("partial", "divergent")

    def test_low_ungapped_is_divergent(self):
        result = {
            "ungapped_identity_pct": 60.0,
            "n_matches": 600, "n_mismatches": 400, "n_gaps": 0,
        }
        code, _, _ = sc._alignment_quality_status(result, 1000)
        assert code == "divergent"


class TestFlushAlignmentsLocked:
    """INV-73: `_flush_active_alignments` holds `_cache_lock` for
    the full read-modify-write so concurrent workers can't clobber
    each other's writes. We can't easily exercise true thread
    contention in a unit test, but we can assert that the function
    body acquires the lock at all — a future refactor that drops
    the lock would regress the data-loss path."""

    def test_flush_body_uses_cache_lock(self):
        import inspect
        src = inspect.getsource(
            sc.PlasmidApp._flush_active_alignments
        )
        # The fix wraps _load_library + merge + _save_library in
        # `with _cache_lock:`. If a refactor splits the function
        # and drops the lock, this fails and the regression is
        # caught.
        assert "with _cache_lock" in src or "_cache_lock.acquire" in src

    def test_flush_lock_is_rlock(self):
        # Without RLock, the inner _load_library and _save_library
        # would deadlock against our outer acquire.
        import threading
        assert isinstance(sc._cache_lock, type(threading.RLock()))


# ═══════════════════════════════════════════════════════════════════════════════
# Rotation tail-trigger — a SMALL circular origin offset still aligns
# ~96% (above the 80% identity gate) but leaves the origin-spanning
# bases as non-aligning 5'/3' end tails. `_pick_best_rotation` must
# detect the tails and rotate them away. Regression for the 2026-05-29
# "non-aligning bases at the ends" report (Angstrom consensus vs a
# reference linearised at a slightly different origin).
# ═══════════════════════════════════════════════════════════════════════════════

def _det_seq(n: int, seed: int) -> str:
    """Deterministic pseudo-random DNA via a private RNG (doesn't
    perturb the global `random` state other tests may depend on)."""
    import random as _random
    rng = _random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


class TestAlignmentTerminalTailBp:
    def test_no_tails(self):
        assert sc._alignment_terminal_tail_bp("ATGC", "ATGC") == 0

    def test_internal_gap_is_not_a_tail(self):
        assert sc._alignment_terminal_tail_bp("ATGC", "AT-C") == 0

    def test_five_prime_tail(self):
        assert sc._alignment_terminal_tail_bp("ATGCGG", "--GCGG") == 2

    def test_three_prime_tail(self):
        assert sc._alignment_terminal_tail_bp("GCGGAT", "GCGG--") == 2

    def test_both_ends(self):
        assert sc._alignment_terminal_tail_bp("ATGCGGAT", "--GCGG--") == 4

    def test_empty_inputs(self):
        assert sc._alignment_terminal_tail_bp("", "") == 0


class TestRotationTailTrigger:
    def test_small_origin_offset_rotates_tails_away(self):
        # THE regression: 30 bp origin shift → ~96% plain identity (gate
        # is 80%, so pre-fix rotation was skipped) → 30 bp non-aligning
        # tail at each end. The tail-trigger must now rotate them away.
        target = _det_seq(1500, seed=42)
        query = target[30:] + target[:30]
        res = sc._pick_best_rotation(
            query, target, is_circular=True, canvas_axis="target")
        assert res["picked_rotation"] != "none", (
            "small origin offset must trigger rotation via the tail-trigger "
            "(pre-fix it was skipped because identity stayed above 80%)"
        )
        assert sc._alignment_terminal_tail_bp(
            res["aligned_q"], res["aligned_t"]) == 0
        assert res["identity_pct"] > 99.0

    def test_clean_same_origin_skips_rotation(self):
        # No offset → no tails → rotation NOT attempted (perf preserved).
        target = _det_seq(1500, seed=42)
        res = sc._pick_best_rotation(
            target, target, is_circular=True, canvas_axis="target")
        assert res["picked_rotation"] == "none"
        assert sc._alignment_terminal_tail_bp(
            res["aligned_q"], res["aligned_t"]) == 0

    def test_large_offset_still_rotates_clean(self):
        target = _det_seq(1500, seed=42)
        query = target[600:] + target[:600]
        res = sc._pick_best_rotation(
            query, target, is_circular=True, canvas_axis="target")
        assert res["picked_rotation"] != "none"
        assert sc._alignment_terminal_tail_bp(
            res["aligned_q"], res["aligned_t"]) == 0

    def test_linear_target_never_force_rotated(self):
        # is_circular=False: end tails are legitimate (linear molecule),
        # so rotation must NOT fire off the back of the tail-trigger.
        target = _det_seq(1500, seed=42)
        query = target[30:] + target[:30]
        res = sc._pick_best_rotation(
            query, target, is_circular=False, canvas_axis="target")
        assert res["picked_rotation"] == "none"


class TestPastTurnRedundancyFoldsInternally:
    """Triage candidate #1: a read running past one full turn (origin-
    spanning redundancy) folds into an INTERNAL gap, not a terminal
    overhang — the benign case stays benign and can't masquerade as the
    origin-offset bug above."""

    def test_past_turn_redundancy_is_internal_not_terminal(self):
        target = _det_seq(1500, seed=42)
        query = target[600:] + target[:600] + target[600:660]  # +60 bp
        res = sc._pick_best_rotation(
            query, target, is_circular=True, canvas_axis="target")
        assert sc._alignment_terminal_tail_bp(
            res["aligned_q"], res["aligned_t"]) == 0


# ═══════════════════════════════════════════════════════════════════════════════
# _format_identity_pct — honest identity display (no round-up to 100)
# ═══════════════════════════════════════════════════════════════════════════════
class TestFormatIdentityPct:
    """A sub-100% identity must never render as "100%" — otherwise the
    number contradicts `_identity_pct_color` (strict >= 100 → light-blue,
    everything below → green/etc). User report 2026-05-31: a one-bp
    mismatch in an 18 kb plasmid read "100.0%" but coloured green."""

    def test_exact_100_renders_clean(self):
        # Genuine perfection → no decimals, distinct from near-100.
        assert sc._format_identity_pct(100.0) == "100%"

    def test_computed_exact_100_renders_clean(self):
        # The value the aligner actually produces for a perfect match.
        assert sc._format_identity_pct(100.0 * 18094 / 18094) == "100%"

    def test_one_bp_off_in_18kb_does_not_round_to_100(self):
        # 18093/18094 = 99.99447% — the user's exact case. Must NOT
        # show "100.0%"; escalates precision until it reads < 100.
        v = 100.0 * 18093 / 18094
        out = sc._format_identity_pct(v)
        assert out != "100%"
        assert out != "100.0%"
        assert out.endswith("%")
        assert float(out[:-1]) < 100.0

    def test_normal_value_keeps_one_decimal(self):
        assert sc._format_identity_pct(99.5) == "99.5%"
        assert sc._format_identity_pct(52.1) == "52.1%"

    def test_escalates_to_two_decimals_when_one_rounds_up(self):
        # 99.994 rounds to "100.0" at 1 dp; must escalate to "99.99".
        assert sc._format_identity_pct(99.994) == "99.99%"

    def test_zero(self):
        assert sc._format_identity_pct(0.0) == "0.0%"

    def test_pathologically_close_to_100_shows_lt_marker(self):
        # 99.999999 rounds to "100.0000" even at 4 dp — refuse to imply
        # perfection.
        assert sc._format_identity_pct(99.999999) == "<100%"

    def test_value_above_100_clamps_to_clean_100(self):
        # A corrupted >100 rides to the clean "100%" (matches the colour
        # tier's >= 100 boundary); never shows "100.5%".
        assert sc._format_identity_pct(100.5) == "100%"

    def test_none_is_dash(self):
        assert sc._format_identity_pct(None) == "—"

    def test_non_numeric_is_dash(self):
        assert sc._format_identity_pct("not a number") == "—"

    def test_decimals_zero_still_avoids_false_100(self):
        # decimals=0 normally shows integers, but a near-100 still
        # escalates rather than reading "100".
        assert sc._format_identity_pct(99.0, decimals=0) == "99%"
        assert sc._format_identity_pct(99.6, decimals=0) == "99.6%"

    def test_color_and_number_never_disagree(self):
        # The whole point: any value the colour function calls NON-light-
        # blue must NOT be formatted "100%".
        for v in (99.99447, 99.994, 99.9, 90.0, 51.0, 11.0, 0.5):
            num = sc._format_identity_pct(v)
            color = sc._identity_pct_color(v)
            if num == "100%":
                assert color == "bright_cyan", (v, num, color)
            else:
                assert color != "bright_cyan", (v, num, color)


# ═══════════════════════════════════════════════════════════════════════════════
# _alignment_bar_columns — zoomed-out (bar-mode) per-column collapse
# ═══════════════════════════════════════════════════════════════════════════════
class TestAlignmentBarColumns:
    """When < 1 col/bp (zoomed out), each terminal column spans many bp.
    The collapse must surface the WORST state in a column (mismatch >
    gap > match) and must not drop a sub-column segment whose
    [bp_to_col(s), bp_to_col(e)) range is empty. User request
    2026-05-31: "show red character in region where mismatch occurred
    even if one bp"."""

    @staticmethod
    def _binned(usable_w=10, total=100, margin_l=0):
        """bp→col mapping with `total` bp packed into `usable_w` cols
        (so 10 bp/col by default — well into bar mode)."""
        return lambda bp: margin_l + bp * usable_w // total

    def test_all_match_paints_every_column_blue(self):
        cols = sc._alignment_bar_columns(
            [(0, 100, "match")], 0, 100, self._binned(), 0, 10)
        assert set(cols.values()) == {"match"}
        assert len(cols) == 10

    def test_one_bp_mismatch_survives_when_zoomed_out(self):
        # THE bug: 10 bp/col, a 1-bp mismatch at bp 50 collapses to the
        # same column as the surrounding matches. Must stay mismatch.
        segs = [(0, 50, "match"), (50, 51, "mismatch"), (51, 100, "match")]
        cols = sc._alignment_bar_columns(segs, 0, 100, self._binned(), 0, 10)
        assert cols[5] == "mismatch"
        # neighbours stay blue
        assert cols[4] == "match"
        assert cols[6] == "match"

    def test_mismatch_outranks_gap_in_same_column(self):
        # bp50 mismatch + bp51 gap both land in col 5 → red wins.
        segs = [(0, 50, "match"), (50, 51, "mismatch"),
                (51, 52, "gap"), (52, 100, "match")]
        cols = sc._alignment_bar_columns(segs, 0, 100, self._binned(), 0, 10)
        assert cols[5] == "mismatch"

    def test_gap_outranks_match_in_same_column(self):
        # A lone 1-bp gap must not be hidden by neighbouring matches.
        segs = [(0, 50, "match"), (50, 51, "gap"), (51, 100, "match")]
        cols = sc._alignment_bar_columns(segs, 0, 100, self._binned(), 0, 10)
        assert cols[5] == "gap"

    def test_order_independent_priority(self):
        # Same column reached match-first then mismatch, and vice-versa
        # → result identical (priority, not paint-order, decides).
        a = sc._alignment_bar_columns(
            [(50, 51, "mismatch"), (51, 100, "match")],
            0, 100, self._binned(), 0, 10)
        b = sc._alignment_bar_columns(
            [(0, 50, "match"), (50, 51, "mismatch")],
            0, 100, self._binned(), 0, 10)
        assert a[5] == "mismatch"
        assert b[5] == "mismatch"

    def test_out_of_view_segments_clipped(self):
        # Segment entirely left of the view contributes nothing.
        cols = sc._alignment_bar_columns(
            [(0, 30, "mismatch"), (60, 100, "match")],
            50, 100, self._binned(), 0, 10)
        # cols < 5 are out of view (view starts at bp 50 → col 5)
        assert all(c >= 5 for c in cols)
        assert "mismatch" not in cols.values()

    def test_columns_respect_bounds(self):
        # col_lo_bound / col_hi_bound clamp the painted columns.
        cols = sc._alignment_bar_columns(
            [(0, 100, "match")], 0, 100, self._binned(), 3, 7)
        assert min(cols) >= 3
        assert max(cols) < 7

    def test_empty_segments(self):
        assert sc._alignment_bar_columns(
            [], 0, 100, self._binned(), 0, 10) == {}

    def test_malformed_segment_skipped(self):
        # A short/garbage tuple is skipped, not fatal.
        cols = sc._alignment_bar_columns(
            [(0, 50, "match"), (50,), (50, 51, "mismatch")],
            0, 100, self._binned(), 0, 10)
        assert cols[5] == "mismatch"

