#!/usr/bin/env python3
"""rotation_ranker_bench.py — empirical comparison of alternative
rank keys for `_pick_best_rotation`.

Background (deferred audit item from 0.9.35):

  The current rank key in `_pick_best_rotation` is
  ``(n_matches, ungapped_identity_pct)``. The audit flagged that
  this can prefer over-padded rotations on length-mismatched pairs:
  a candidate with 1001 matches + 5000 gap columns wins over one
  with 1000 matches + 0 gaps even though the latter is biologically
  more meaningful.

  This script generates synthetic plasmid pairs with KNOWN ground
  truth (rotation offset, RC, length-mismatch), runs both the
  current ranker and an alternative ``score``-based ranker on each,
  and reports:

    * total pairs
    * disagreement rate
    * which ranker matches ground truth more often

  A score-based ranker uses the Biopython aligner's own score
  (already in the result dict) which naturally encodes gap penalties.

  If the alternative wins by a clear margin (>= 5 % more correct
  picks) AND doesn't regress on the "no rotation needed" case, the
  ranker change is worth landing. Otherwise current ranker stays
  documented and unchanged.

Usage:

    python3 scripts/rotation_ranker_bench.py [--n=1000]

Sandbox required (writes to splicecraft data dir if run unsandboxed):

    XDG_DATA_HOME=$(mktemp -d) python3 scripts/rotation_ranker_bench.py
"""
from __future__ import annotations

import argparse
import os
import random
import sys
import tempfile
from pathlib import Path as _Path

# Add repo root to sys.path so `import splicecraft` resolves whether
# the script is run from the repo root or `scripts/`.
_REPO_ROOT = _Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# Mandatory sandbox per CLAUDE.md SACRED rule — never run without an
# XDG_DATA_HOME override; the bench's import alone would compute
# `_DATA_DIR` against the user's real plasmid library.
if not os.environ.get("XDG_DATA_HOME"):
    os.environ["XDG_DATA_HOME"] = tempfile.mkdtemp(prefix="rotbench-")
os.environ.setdefault("SPLICECRAFT_SKIP_LOCK", "1")

import splicecraft as sc  # noqa: E402

_sandbox_root = tempfile.gettempdir()
assert str(sc._state._DATA_DIR).startswith(_sandbox_root), (
    f"refusing to run unsandboxed (data dir = {sc._state._DATA_DIR}; "
    f"expected under {_sandbox_root})"
)


_BASES = "ACGT"


def _rand_seq(length: int, rng: random.Random) -> str:
    return "".join(rng.choices(_BASES, k=length))


def _rc(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def _gen_scenario(rng: random.Random, length: int = 2000) -> dict:
    """Generate a pair + record the ground truth."""
    target = _rand_seq(length, rng)
    # Optionally inject up to 1% point mutations into the query so the
    # alignment doesn't always score 100 % (more representative of
    # real reads).
    mutations = rng.randrange(0, max(1, length // 100))
    truth = {
        "rotation": rng.choice([0, 50, 173, length // 3, length // 2]),
        "rc": rng.choice([False, True]),
        "length_delta": rng.choice([0, 0, 0, +200, -200]),
        "mutations": mutations,
    }
    q_base = target
    if truth["length_delta"]:
        ld = truth["length_delta"]
        if ld > 0:
            insert = _rand_seq(ld, rng)
            q_base = q_base + insert
        else:
            q_base = q_base[: max(1, len(q_base) + ld)]
    # Apply rotation.
    r = truth["rotation"] % max(1, len(q_base))
    q_rot = q_base[r:] + q_base[:r]
    # Apply RC if needed.
    q_final = _rc(q_rot) if truth["rc"] else q_rot
    # Apply mutations to the query (not the target).
    if mutations:
        q_list = list(q_final)
        for _ in range(mutations):
            pos = rng.randrange(0, len(q_list))
            q_list[pos] = rng.choice([b for b in _BASES if b != q_list[pos]])
        q_final = "".join(q_list)
    return {
        "query": q_final,
        "target": target,
        "truth": truth,
    }


def _picked_matches_truth(result: dict, truth: dict) -> bool:
    """Heuristic: does the picker's `picked_rotation` + `query_rc`
    match the truth's rotation/RC?"""
    picked = result.get("picked_rotation", "")
    picked_rc = bool(result.get("query_rc", False))
    if truth["rotation"] == 0 and not truth["rc"]:
        # Plain pair — picker should choose "none"/"plain" (no
        # rotation, no RC).
        return picked in ("none", "plain") and picked_rc is False
    if truth["rc"] and picked_rc:
        # RC detected. Rotation may or may not be applied; either
        # is acceptable for our match-truth check.
        return True
    if truth["rotation"] != 0 and picked in ("query", "target"):
        # Rotation detected on some axis. Either query- or target-
        # axis rotation can produce the correct alignment.
        return True
    return False


def _bench_alternative_score_ranker(query: str, target: str) -> dict:
    """Re-implement `_pick_best_rotation` with a score-based rank
    key instead of `(n_matches, ungapped_identity_pct)`. The score
    field is the Biopython aligner's own score which includes gap
    penalties — bias toward less-padded alignments."""
    # Reuse the picker's existing candidate generation by calling
    # it, then re-rank the result. Since the existing function
    # returns ONLY the winner, we'd need to expose candidates.
    # For this bench, just re-run `_pairwise_align` on each
    # transform and rank by score.
    import itertools

    qn = len(query)
    candidates = []
    for is_rc, rot_kind, offset in itertools.product(
            (False, True),
            ("none", "query"),
            (0, qn // 4, qn // 2, qn // 3),
    ):
        if rot_kind == "none" and offset != 0:
            continue
        eff_q = _rc(query) if is_rc else query
        if rot_kind == "query":
            eff_q = eff_q[offset:] + eff_q[:offset]
        try:
            r = sc._pairwise_align(eff_q, target, mode="global")
        except Exception:
            continue
        candidates.append(
            (r.get("score", 0.0), r.get("n_matches", 0), is_rc, rot_kind, offset, r),
        )
    if not candidates:
        return {"picked_rotation": "none", "query_rc": False,
                 "identity_pct": 0.0, "score": 0.0}
    # Rank by score primary, n_matches secondary.
    candidates.sort(key=lambda c: (c[0], c[1]), reverse=True)
    score, n_match, is_rc, rot_kind, offset, r = candidates[0]
    return {
        "picked_rotation": rot_kind,
        "query_rc": is_rc,
        "query_rotation": offset if rot_kind == "query" else 0,
        "identity_pct": r.get("identity_pct", 0.0),
        "score": score,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=200,
                          help="number of synthetic pairs to test")
    parser.add_argument("--length", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--report", default="rotation_ranker_bench.md",
                          help="output report path")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    current_correct = 0
    alternative_correct = 0
    disagreements = 0
    plain_pairs = 0
    rotation_pairs = 0
    rc_pairs = 0

    for i in range(args.n):
        scenario = _gen_scenario(rng, length=args.length)
        truth = scenario["truth"]
        if truth["rotation"] == 0 and not truth["rc"]:
            plain_pairs += 1
        if truth["rotation"] != 0:
            rotation_pairs += 1
        if truth["rc"]:
            rc_pairs += 1
        try:
            current = sc._pick_best_rotation(
                scenario["query"], scenario["target"],
                is_circular=True, mode="global",
            )
        except ValueError:
            continue
        try:
            alternative = _bench_alternative_score_ranker(
                scenario["query"], scenario["target"],
            )
        except Exception:
            continue
        cur_ok = _picked_matches_truth(current, truth)
        alt_ok = _picked_matches_truth(alternative, truth)
        if cur_ok:
            current_correct += 1
        if alt_ok:
            alternative_correct += 1
        if cur_ok != alt_ok:
            disagreements += 1
        if (i + 1) % 25 == 0:
            print(
                f"  [{i + 1}/{args.n}] cur={current_correct} "
                f"alt={alternative_correct} disagreements={disagreements}",
                file=sys.stderr,
            )

    print()
    print(f"Synthetic pair count:   {args.n}")
    print(f"  - plain (no rot/RC):  {plain_pairs}")
    print(f"  - with rotation:      {rotation_pairs}")
    print(f"  - with RC:            {rc_pairs}")
    print()
    print(f"Current ranker correct:     {current_correct} ({100*current_correct/args.n:.1f}%)")
    print(f"Alternative ranker correct: {alternative_correct} ({100*alternative_correct/args.n:.1f}%)")
    print(f"Disagreements:              {disagreements}")
    delta = alternative_correct - current_correct
    print()
    print(f"Δ (alt - current): {delta:+d}")
    if delta >= max(5, args.n // 20):
        verdict = "ALTERNATIVE RANKER WINS by a clear margin — consider landing."
    elif delta <= -max(5, args.n // 20):
        verdict = "CURRENT RANKER WINS — keep current behavior."
    else:
        verdict = (
            "TIE / WITHIN NOISE — current ranker is fine; revisit only "
            "if real-world data exposes a specific regression."
        )
    print(f"Verdict: {verdict}")
    # Write a brief markdown report so the result is reviewable later.
    with open(args.report, "w", encoding="utf-8") as f:
        f.write("# Rotation Ranker Synthetic-Pair Sweep\n\n")
        f.write(f"- Total pairs: **{args.n}**\n")
        f.write(f"- Plain (no rot/RC): {plain_pairs}\n")
        f.write(f"- With rotation: {rotation_pairs}\n")
        f.write(f"- With RC: {rc_pairs}\n\n")
        f.write("## Results\n\n")
        f.write("| Ranker | Correct | Pct |\n")
        f.write("|---|---|---|\n")
        f.write(f"| Current `(n_matches, ungapped_identity_pct)` | {current_correct} | {100*current_correct/args.n:.1f}% |\n")
        f.write(f"| Alternative `(score, n_matches)` | {alternative_correct} | {100*alternative_correct/args.n:.1f}% |\n")
        f.write("\n")
        f.write(f"- Disagreements: {disagreements}\n")
        f.write(f"- Δ: **{delta:+d}**\n\n")
        f.write(f"## Verdict\n\n{verdict}\n")
    print(f"\nReport written to {args.report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
