"""
test_chromosome_segments — designing at chromosome scale, `[INV-223]`.

A 1 Gbp target cannot be ordered as one piece, so designing one means
deciding where to cut it and PROVING the pieces can only reassemble one way.

The headline test is reassembly: walking the plan must reproduce the original
sequence exactly. A segmenter whose overlaps don't actually match, or whose
segments leave a gap, produces a plan that looks fine in a table and fails at
the bench — and the table is all anyone reads.

The second is junction uniqueness. Two junctions sharing an overlap give an
assembly that is chemically happy and structurally wrong: segments swap or
insert inverted, and it surfaces after sequencing. That is an `error`, not a
warning.
"""
from __future__ import annotations

import random

import pytest

import splicecraft as sc
import splicecraft_seqanalysis as sa


def _seq(n, seed=5):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _plan(seq, **kw):
    kw.setdefault("strategy", "gibson")
    return sa._plan_chromosome_segments(seq, **kw)


# ═══════════════════════════════════════════════════════════════════════════
# Structural correctness: the plan must describe the actual molecule
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("strategy,segment_bp", [
    ("gibson", 15_000),
    ("gibson", 7_000),
    ("yeast_hr", 25_000),
    ("goldengate", 12_000),
])
def test_the_plan_reassembles_the_original_sequence(strategy, segment_bp):
    """THE test. Walk the segments, drop each one's leading overlap (it is the
    previous segment's tail), concatenate — and get the original back. A plan
    that cannot be reassembled is a table, not a design."""
    seq = _seq(100_000)
    p = _plan(seq, strategy=strategy, segment_bp=segment_bp, circular=True)
    assert p["ok"], p["errors"]
    ov = p["overlap_bp"]
    rebuilt = p["segments"][0]["seq"]
    for sg in p["segments"][1:]:
        assert sg["seq"][:ov] == rebuilt[-ov:], (
            f"segment {sg['index']}'s 5' overlap does not match the previous "
            f"segment's 3' end — the pieces would not join")
        rebuilt += sg["seq"][ov:]
    # Circular: the last segment's tail wraps onto the first segment's head,
    # so the walk ends one overlap past the origin.
    assert rebuilt[-ov:] == p["segments"][0]["seq"][:ov], \
        "the closing junction does not match"
    assert len(rebuilt) == len(seq) + ov, (
        f"the walk should end exactly one overlap past the origin: "
        f"{len(rebuilt)} vs {len(seq)}+{ov}")
    assert rebuilt[:len(seq)] == seq, "the plan does not reproduce the target"


def test_linear_plan_covers_the_whole_target_with_free_ends():
    seq = _seq(50_000)
    p = _plan(seq, segment_bp=12_000, circular=False)
    assert p["ok"], p["errors"]
    segs = p["segments"]
    assert segs[0]["overlap_5"] == "", "a linear build's first end is free"
    assert segs[-1]["overlap_3"] == "", "a linear build's last end is free"
    ov = p["overlap_bp"]
    rebuilt = segs[0]["seq"]
    for sg in segs[1:]:
        assert sg["seq"][:ov] == rebuilt[-ov:]
        rebuilt += sg["seq"][ov:]
    assert rebuilt == seq


def test_every_base_is_covered_by_some_segment():
    seq = _seq(40_000)
    p = _plan(seq, segment_bp=9_000, circular=True)
    covered = bytearray(len(seq))
    for sg in p["segments"]:
        for i in range(sg["length"]):
            covered[(sg["start"] + i) % len(seq)] = 1
    assert all(covered), f"{len(covered) - sum(covered)} bases uncovered"


@pytest.mark.parametrize("strategy,want", [
    ("gibson", 30), ("yeast_hr", 120), ("goldengate", 4),
])
def test_strategy_sets_the_overlap_size(strategy, want):
    """The overlap is the chemistry's, not a preference."""
    p = _plan(_seq(80_000), strategy=strategy, segment_bp=20_000)
    assert p["overlap_bp"] == want
    assert all(j["length"] == want for j in p["junctions"])


# ═══════════════════════════════════════════════════════════════════════════
# Junction uniqueness — the check that decides whether a build works
# ═══════════════════════════════════════════════════════════════════════════

def test_identical_junctions_are_an_error_when_boundaries_are_fixed():
    """A tandem repeat makes every boundary at the same phase identical. With
    `nudge_bp=0` the planner must take those cuts and REFUSE: the assembly is
    ambiguous and segments can swap."""
    unit = _seq(10_000, seed=3)
    p = _plan(unit * 6, segment_bp=10_030, circular=True, nudge_bp=0)
    assert p["ok"] is False
    assert p["n_ambiguous_junctions"] >= 1
    assert any("AMBIGUOUS" in e for e in p["errors"])
    assert any(j["unique"] is False for j in p["junctions"])


def test_nudging_resolves_a_repeat_that_fixed_boundaries_cannot():
    """And with nudging ON the planner SOLVES it — moving each cut to a
    distinct phase — then still says the homology is repeated elsewhere,
    which is the honest residual risk of building in a repeat."""
    unit = _seq(10_000, seed=3)
    p = _plan(unit * 6, segment_bp=10_030, circular=True)
    assert p["n_ambiguous_junctions"] == 0, p["errors"]
    assert any("moved" in w for w in p["warnings"]), p["warnings"]
    # Each junction's homology still occurs once per repeat copy.
    assert any(j.get("offtarget_hits") for j in p["junctions"]), \
        "building inside a tandem repeat must still warn about off-targets"


def test_a_reverse_complement_junction_is_caught():
    """One junction being another's reverse complement lets a segment insert
    INVERTED — chemically fine, structurally wrong.

    Constructed exactly: with `overlap_bp` 30 and `segment_bp` 1030 the
    stride is 1000, so junction i's window is `seq[(i+1)*1000 : +30]`. Plant
    a window at 1000 and its reverse complement at 2000."""
    OV, STRIDE = 30, 1_000
    base = list(_seq(5 * STRIDE, seed=8))
    w = _seq(OV, seed=77)
    base[STRIDE:STRIDE + OV] = list(w)
    base[2 * STRIDE:2 * STRIDE + OV] = list(sc._rc(w))
    target = "".join(base)
    p = _plan(target, segment_bp=STRIDE + OV, overlap_bp=OV, circular=True,
              nudge_bp=0)
    assert p["junctions"][0]["overlap"] == w, "fixture did not land the window"
    assert p["ok"] is False, "an invertible assembly must not pass"
    assert any("reverse complement" in e for e in p["errors"]), p["errors"]
    # Nudging should rescue the same target.
    fixed = _plan(target, segment_bp=STRIDE + OV, overlap_bp=OV,
                  circular=True)
    assert fixed["n_ambiguous_junctions"] == 0, fixed["errors"]


def test_unique_junctions_pass_cleanly():
    p = _plan(_seq(120_000, seed=21), segment_bp=20_000)
    assert p["ok"] is True, p["errors"]
    assert p["n_ambiguous_junctions"] == 0
    assert all(j["unique"] for j in p["junctions"])


def test_offtarget_homology_is_reported():
    """An overlap that also matches elsewhere can recombine there instead —
    which long homology arms make MORE likely, not less."""
    body = _seq(9_000, seed=4)
    dup = _seq(30, seed=99)
    # The junction window appears again in the middle of a later segment.
    target = dup + body + dup + body + body
    p = _plan(target, segment_bp=len(dup) + len(body), overlap_bp=30)
    assert p["offtarget_checked"] is True
    flagged = [j for j in p["junctions"] if j.get("offtarget_hits")]
    assert flagged or not p["ok"], (
        "a repeated junction window was neither flagged off-target nor "
        "reported ambiguous")


def test_offtarget_check_says_when_it_was_skipped(monkeypatch):
    """"None found" must never be confusable with "never looked"."""
    monkeypatch.setattr(sa, "_SEGMENT_OFFTARGET_MAX_BP", 1_000)
    p = _plan(_seq(60_000), segment_bp=12_000)
    assert p["offtarget_checked"] is False
    assert any("NOT checked" in w for w in p["warnings"])


def test_goldengate_flags_a_palindromic_overhang():
    """A palindromic 4 nt overhang ligates to itself."""
    body = _seq(5_000, seed=12)
    target = (body + "GGCC") * 4          # GGCC is its own reverse complement
    p = _plan(target, strategy="goldengate", segment_bp=len(body) + 4,
              nudge_bp=0)
    pal = [w for j in p["junctions"] for w in j["warnings"]
           if "palindromic" in w]
    assert pal or not p["ok"], "a self-ligating overhang went unreported"
    # With nudging the planner moves off the palindrome entirely.
    fixed = _plan(target, strategy="goldengate", segment_bp=len(body) + 4)
    assert not [w for j in fixed["junctions"] for w in j["warnings"]
                if "palindromic" in w]


# ═══════════════════════════════════════════════════════════════════════════
# Forbidden sites belong to the ASSEMBLY, not to DNA
# ═══════════════════════════════════════════════════════════════════════════

def _kinds(plan):
    from collections import Counter
    return Counter(w.get("kind") for sg in plan["segments"]
                   for w in sg["warnings"])


@pytest.mark.parametrize("strategy", ["gibson", "yeast_hr"])
def test_type_iis_sites_are_not_flagged_for_non_goldengate(strategy):
    """Lint-by-default reported 281 "errors" per 20 kb of ordinary sequence
    for a Gibson build, which neither chemistry consults. Noise that teaches
    the reader to skim the report it exists to focus them on."""
    p = _plan(_seq(100_000), strategy=strategy, segment_bp=20_000)
    assert _kinds(p).get("type_iis_site", 0) == 0


def test_goldengate_does_flag_type_iis_sites():
    p = _plan(_seq(100_000), strategy="goldengate", segment_bp=20_000)
    assert _kinds(p).get("type_iis_site", 0) > 0


def test_explicit_forbidden_enzymes_override_the_strategy():
    p = _plan(_seq(100_000), strategy="gibson", segment_bp=20_000,
              forbidden_enzymes=["BsaI"])
    assert _kinds(p).get("type_iis_site", 0) > 0


# ═══════════════════════════════════════════════════════════════════════════
# Refusals — each one actionable
# ═══════════════════════════════════════════════════════════════════════════

def test_unknown_strategy_lists_the_real_ones():
    p = _plan(_seq(10_000), strategy="nope")
    assert p["ok"] is False
    assert "gibson" in p["errors"][0]


def test_segment_smaller_than_its_overlap_is_refused():
    p = _plan(_seq(10_000), segment_bp=10)
    assert p["ok"] is False
    assert "overlap" in p["errors"][0]


def test_segment_bigger_than_the_target_is_refused():
    p = _plan(_seq(1_000), segment_bp=50_000)
    assert p["ok"] is False
    assert "does not need segmenting" in p["errors"][0]


def test_overlap_outside_the_chemistry_is_refused():
    p = _plan(_seq(50_000), strategy="gibson", overlap_bp=5, segment_bp=10_000)
    assert p["ok"] is False
    assert "18..80" in p["errors"][0]


def test_too_many_segments_refuses_with_a_usable_suggestion():
    """A 100,000-row plan is not executable; the refusal must say what to do
    instead rather than just declining."""
    p = _plan(_seq(200_000), segment_bp=300, max_segments=20)
    assert p["ok"] is False
    assert "segment_bp >=" in p["errors"][0]


def test_empty_sequence_is_refused():
    assert _plan("")["ok"] is False
    assert _plan("NNNN")["ok"] is False


# ═══════════════════════════════════════════════════════════════════════════
# The side door
# ═══════════════════════════════════════════════════════════════════════════

def _call(_ep, _app=None, **payload):
    entry = sc._state._AGENT_HANDLERS.get(_ep)
    assert entry is not None, f"{_ep} not registered"
    return entry[0](_app, payload)


def _status(res):
    return res[1] if isinstance(res, tuple) else 200


def _body(res):
    return res[0] if isinstance(res, tuple) else res


def test_strategies_are_discoverable():
    b = _body(_call("list-assembly-strategies"))
    names = {s["name"] for s in b["strategies"]}
    assert names == set(sa._SEGMENT_STRATEGIES)
    gib = next(s for s in b["strategies"] if s["name"] == "gibson")
    assert gib["junction_tm_range"] is not None
    yhr = next(s for s in b["strategies"] if s["name"] == "yeast_hr")
    assert yhr["overlap_max"] > gib["overlap_max"]


def test_endpoint_requires_an_explicit_strategy():
    """The chemistry sets the overlap size, so a default would silently pick
    someone's assembly method for them."""
    res = _call("plan-chromosome-segments", sequence=_seq(40_000))
    assert _status(res) == 400
    assert "strategy" in _body(res)["error"]


def test_endpoint_omits_sequences_by_default():
    """A 2,000-segment plan carrying every base is tens of MB."""
    b = _body(_call("plan-chromosome-segments", sequence=_seq(60_000),
                    strategy="gibson", segment_bp=12_000))
    assert b["sequences_omitted"] is True
    assert "seq" not in b["segments"][0]
    assert b["segments"][0]["overlap_3"], "the junctions must still be there"


def test_endpoint_can_return_sequences_on_request():
    b = _body(_call("plan-chromosome-segments", sequence=_seq(60_000),
                    strategy="gibson", segment_bp=12_000,
                    include_sequences=True))
    assert len(b["segments"][0]["seq"]) == 12_000


def test_endpoint_404s_on_an_unknown_entry():
    assert _status(_call("plan-chromosome-segments", name="ghost",
                         strategy="gibson")) == 404
