"""The answer for a read must not depend on which bases happen to sit beside
an event ([INV-213]).

The global engine is unit-cost: it charges a gap per column and nothing to
open one, so splitting one deletion into three is free, and an aligner that
can buy a chance match by splitting will. Every test here plants ONE known
event in many random plasmids and asserts the same answer comes back for all
of them. Before the fix the same read shape gave up to 95 different answers
over 120 plasmids, and a perfect read missing its first 60 bp of a linear
reference was a "real change" about one time in three.

Each property runs over `_SEEDS` plasmids: one lucky seed proves nothing (the
random-sequence test this replaces passed on its seed while real pUC19 failed).
"""
from __future__ import annotations

import random

import pytest

import splicecraft as sc

_SEEDS = range(24)
_N = 2000


def _ref(seed: int, n: int = _N) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mut(s: str, i: int) -> str:
    return s[:i] + {"A": "C", "C": "G", "G": "T", "T": "A"}[s[i]] + s[i + 1:]


def _verify(ref: str, read: str, circular: bool) -> dict:
    return sc._state._AGENT_HANDLERS["verify-against-reads"][0](
        None, {"reference": ref, "reads": [read], "circular": circular,
               "read_quality": [[40] * len(read)]})


def _calls(r: dict) -> "tuple[str, list[tuple[str, int]]]":
    vs = [v for v in r["consensus"]["variants"] if v["type"] != "truncated"]
    return (r["read_quality"][0]["verdict"],
            sorted((v["type"], v["length"]) for v in vs))


def _rebuild(ref: str, variants: list) -> str:
    """Apply reported variants to the reference. At one position the base
    change goes first and the insertion (which sits BEFORE that base) after."""
    seq = ref
    for v in sorted(variants, key=lambda v: (-v["target_pos"],
                                             v["type"] == "insertion")):
        p = v["target_pos"]
        if v["type"] == "snp":
            seq = seq[:p] + v["alt"] + seq[p + 1:]
        elif v["type"] == "deletion":
            seq = seq[:p] + seq[p + v["length"]:]
        elif v["type"] == "insertion":
            seq = seq[:p] + v["alt"] + seq[p:]
    return seq


# ══════════════════════════════════════════════════════════════════════════
# A read's ENDS
# ══════════════════════════════════════════════════════════════════════════

class TestLinearReadEnds:
    """On a linear reference a read may stop short of either end. The bases
    it never reached are not a deletion — and must not become one because the
    read's first bases happened to match the reference's first bases (pUC19
    begins TC; so does the read at bp 100)."""

    @pytest.mark.parametrize("shape", [
        "missing_start", "missing_end", "missing_both",
        "overhang_start", "overhang_end",
    ])
    def test_a_perfect_read_is_clean_on_every_plasmid(self, shape):
        X = "GATTACAGGT"
        make = {
            "missing_start": lambda s: s[60:],
            "missing_end": lambda s: s[:-60],
            "missing_both": lambda s: s[30:-30],
            "overhang_start": lambda s: X + s,
            "overhang_end": lambda s: s + X,
        }[shape]
        bad = []
        for seed in _SEEDS:
            ref = _ref(seed)
            got = _calls(_verify(ref, make(ref), circular=False))
            if got != ("clean", []):
                bad.append((seed, got))
        assert not bad, bad

    @pytest.mark.parametrize("where", ["start_first_base", "start_third_base",
                                       "end_last_base"])
    def test_an_error_beside_a_missing_end_is_one_snp(self, where):
        """The read's own base at its edge is wrong as well. Pinned to bp 0
        by a chance match it read as a 60 bp deletion, and the SNP was gone."""
        make = {
            "start_first_base": lambda s: _mut(s[60:], 0),
            "start_third_base": lambda s: _mut(s[60:], 2),
            "end_last_base": lambda s: _mut(s[:-60], _N - 61),
        }[where]
        bad = []
        for seed in _SEEDS:
            ref = _ref(seed)
            got = _calls(_verify(ref, make(ref), circular=False))
            if got != ("real_changes", [("snp", 1)]):
                bad.append((seed, got))
        assert not bad, bad

    def test_a_real_indel_near_a_read_start_is_still_called(self):
        """The other side of making the read's ends cheap: they must not be
        cheap enough to explain away a real 1 bp deletion a few bases in."""
        bad = []
        for seed in _SEEDS:
            ref = "GGGG" + "A" * 8 + "CCCC" + _ref(seed, 400)
            read = ref[:6] + ref[7:]
            got = _calls(_verify(ref, read, circular=False))
            if got != ("real_changes", [("deletion", 1)]):
                bad.append((seed, got))
        assert not bad, bad

    def test_an_snp_at_the_first_and_last_base_is_still_reported(self):
        for seed in _SEEDS:
            ref = _ref(seed)
            for read in (_mut(ref, 0), _mut(ref, _N - 1)):
                assert _calls(_verify(ref, read, circular=False)) == \
                    ("real_changes", [("snp", 1)]), seed


class TestCircularSeam:
    """A circle has no ends: a full-length read missing bases at the origin
    carries a deletion there, whatever it was rotated to."""

    @pytest.mark.parametrize("kind", ["deletion", "insertion"])
    def test_an_origin_indel_is_exactly_one_event(self, kind):
        bad = []
        for seed in _SEEDS:
            ref = _ref(seed)
            read = ref[60:] if kind == "deletion" else "GATTACAGGT" + ref
            want = [(kind, 60 if kind == "deletion" else 10)]
            got = _calls(_verify(ref, read, circular=True))
            if got != ("real_changes", want):
                bad.append((seed, got))
        assert not bad, bad


# ══════════════════════════════════════════════════════════════════════════
# Events with a neighbour
# ══════════════════════════════════════════════════════════════════════════

_NEIGHBOURS = {
    # name: (make read, planted number of events)
    "deletion_60_snp_3_after": (lambda s: s[:900] + _mut(s[960:], 2), 2),
    "deletion_60_snp_3_before": (lambda s: _mut(s[:900], 897) + s[960:], 2),
    "deletion_10_snp_2_after": (lambda s: s[:900] + _mut(s[910:], 1), 2),
    "insertion_10_snp_2_after": (
        lambda s: s[:900] + "GATTACAGGT" + _mut(s[900:], 1), 2),
    "two_deletions_5_apart": (lambda s: s[:900] + s[920:925] + s[945:], 2),
}


class TestEventsWithANeighbour:
    """An indel with an SNP or a second indel a few bp away. Unit-cost
    alignment split the indel to make the neighbour's base match by chance:
    a 60 bp deletion with an SNP 3 bp away came back as 78 different answers
    over 120 plasmids. Now every answer must be a valid explanation of the
    read (it rebuilds the read exactly) with no more events than planted —
    fewer is allowed, when one event genuinely explains the read."""

    @pytest.mark.parametrize("circular", [False, True])
    @pytest.mark.parametrize("name", sorted(_NEIGHBOURS))
    def test_the_calls_explain_the_read_with_no_extra_events(self, name,
                                                             circular):
        make, planted = _NEIGHBOURS[name]
        bad = []
        for seed in _SEEDS:
            ref = _ref(seed)
            read = make(ref)
            r = _verify(ref, read, circular=circular)
            vs = [v for v in r["consensus"]["variants"]
                  if v["type"] != "truncated"]
            if _rebuild(ref, vs) != read or len(vs) > planted:
                bad.append((seed, [(v["type"], v["target_pos"], v["length"])
                                   for v in vs]))
        assert not bad, bad

    def test_the_planted_answer_is_the_usual_answer(self):
        """Not just valid: for a deletion with an SNP 3 bp away the planted
        pair is what comes back almost every time (119/120 measured; the odd
        one out is a single deletion that explains the read exactly)."""
        make, _ = _NEIGHBOURS["deletion_60_snp_3_after"]
        hits = sum(
            _calls(_verify(_ref(s), make(_ref(s)), circular=False))
            == ("real_changes", [("deletion", 60), ("snp", 1)])
            for s in _SEEDS)
        assert hits >= len(_SEEDS) - 2, hits


# ══════════════════════════════════════════════════════════════════════════
# The pieces
# ══════════════════════════════════════════════════════════════════════════

class TestSlideGapsToEnds:
    def test_a_chance_anchor_is_moved_off_the_start(self):
        # read = ref[3:] where ref[3:5] == ref[0:2]: "TC" anchored at bp 0.
        aq = "TC---GGA"
        at = "TCATCGGA"
        assert sc._slide_gaps_to_ends(aq, at) == ("---TCGGA", "TCATCGGA")

    def test_a_gap_that_cannot_slide_without_a_mismatch_stays(self):
        aq = "TA---GGA"     # 'A' would land on 'C' — a mismatch: never worse only
        at = "TCATCGGA"
        assert sc._slide_gaps_to_ends(aq, at) == (aq, at)

    def test_the_end_side_mirrors_the_start(self):
        aq = "GGAT---C"
        at = "GGATCTAC"
        assert sc._slide_gaps_to_ends(aq, at) == ("GGATC---", at)

    def test_a_read_overhang_slides_too(self):
        # read = "TCA" + ref, with ref's "TC" anchored on the read's first two
        # bases: the overhang belongs before the reference, not inside it.
        aq = "TCATCGGAT"
        at = "TC---GGAT"
        assert sc._slide_gaps_to_ends(aq, at) == (aq, "---TCGGAT")

    def test_round_trip_and_degenerate_rows(self):
        assert sc._slide_gaps_to_ends("", "") == ("", "")
        assert sc._slide_gaps_to_ends("AC", "A") == ("AC", "A")
        rng = random.Random(3)
        for _ in range(300):
            a = "".join(rng.choice("ACGT-") for _ in range(rng.randint(1, 30)))
            b = "".join(rng.choice("ACGT") for _ in range(len(a)))
            q, t = sc._slide_gaps_to_ends(a, b)
            assert len(q) == len(t) == len(a)
            assert q.replace("-", "") == a.replace("-", "")
            assert t.replace("-", "") == b.replace("-", "")


class TestVariantGate:
    def test_linear_insertions_past_either_end_are_overhang(self):
        gate = sc._variant_gate_pos
        assert gate({"type": "insertion", "target_pos": 0}, 100, False) == -1
        assert gate({"type": "insertion", "target_pos": 100}, 100, False) == -1
        assert gate({"type": "insertion", "target_pos": 1}, 100, False) == 1

    def test_on_a_circle_after_the_last_base_is_bp_0(self):
        gate = sc._variant_gate_pos
        assert gate({"type": "insertion", "target_pos": 100}, 100, True) == 0
        assert gate({"type": "insertion", "target_pos": 0}, 100, True) == 0

    def test_other_variants_gate_on_their_position(self):
        gate = sc._variant_gate_pos
        assert gate({"type": "deletion", "target_pos": 0}, 100, False) == 0
        assert gate({"type": "snp", "target_pos": 99}, 100, False) == 99

    def test_every_gating_site_shares_the_one_gate(self):
        import splicecraft_agent as ag
        import splicecraft_biology as bio
        assert sc._variant_gate_pos is bio._variant_gate_pos
        assert ag._variant_gate_pos is bio._variant_gate_pos


class TestRefineIndelClusters:
    def test_a_clean_single_gap_is_left_exactly_where_it_was(self):
        ref = _ref(5, 400)
        aq = ref[:200] + "-" * 10 + ref[210:]
        assert sc._refine_indel_clusters(aq, ref) == (aq, ref)

    def test_a_split_deletion_is_rejoined(self):
        ref = _ref(6, 400)
        read = ref[:200] + _mut(ref[260:], 2)
        q, t = sc._myers_align_global(read, ref)
        q2, t2 = sc._refine_indel_clusters(q, t)
        assert q2.replace("-", "") == read and t2.replace("-", "") == ref
        assert q2.count("-") == 60 and "-" * 60 in q2

    def test_divergent_sequence_is_not_touched(self):
        # A flipped block reads as dense chance matches: not an event.
        ref = _ref(7, 3000)
        read = ref[:1000] + sc._rc(ref[1000:1800]) + ref[1800:]
        q, t = sc._myers_align_global(read, ref)
        q2, _ = sc._refine_indel_clusters(q, t)
        assert q2[1100:1700] == q[1100:1700]

    def test_unread_plasmid_beside_a_partial_read_is_never_entered(self):
        ref = _ref(8, 3000)
        read = ref[1500:2300]
        aq = "-" * 1500 + read + "-" * 700
        assert sc._refine_indel_clusters(aq, ref) == (aq, ref)


class TestReadRotation:
    def test_the_two_rotations_compose(self):
        f = sc._read_rotation_in_rows
        assert f({"query_rotation": 99, "query_frame_shift": 0}) == 99
        assert f({"query_rotation": 0, "query_frame_shift": 198}) == 198
        assert f({"query_rotation": float("inf")}) == 0
        assert f(None) == 0

    def test_the_extent_follows_a_query_rotation(self):
        """Rows as a query-rotation pick leaves them: a 10 bp read covering
        ref[15:20] + ref[0:5], rotated by 5 before aligning. Its own ends are
        at bp 15 and bp 4 — not at the first and last aligned columns, which
        would claim the 10 bp it never read."""
        ref = "ACGTTGCAAGCTTAGCCATG"
        aq = ref[0:5] + "-" * 10 + ref[15:20]
        res = {"aligned_q": aq, "aligned_t": ref, "query_rotation": 5,
               "query_frame_shift": 0}
        ext = sc._read_extent_for(res, len(ref), circular=True)
        assert sorted(ext) == [(0, 5), (15, 20)], ext

    def test_a_query_rotated_partial_read_keeps_its_own_ends(self):
        """A 500 bp read across the origin, picked with the READ rotated: its
        unread 999 bp arc came back as a deletion, because only a target
        pick's frame shift was used to find where the read begins."""
        core = _ref(9, 1493)
        mut = "AAA" + core + "AAA"
        read = (mut[1400:] + mut[:1400])[:500]
        r = sc._pick_best_rotation(read, "A" + mut, is_circular=True,
                                   mode="global")
        ext = sc._read_extent_for(r, len(mut) + 1, circular=True)
        assert ext is not None
        covered = sum(hi - lo for lo, hi in ext)
        assert covered < 600, (r["picked_rotation"], ext)
