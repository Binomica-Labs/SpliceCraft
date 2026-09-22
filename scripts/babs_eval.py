#!/usr/bin/env python3
"""Babs live test bench — real agent tasks against your local Ollama models.

Runs a fixed set of everyday requests ("where is CmR?", "which enzymes cut
once?", "design primers for AmpR", "add a feature at 100–200", …) through the
REAL Babs agent loop — the same `BabsScreen._generate_agentic`, dispatch and
endpoints the app uses — against a real local model, and grades every answer
against exact ground truth. Use it to compare models, the native vs JSON tool
protocols, or the effect of a prompt change, instead of guessing.

    python scripts/babs_eval.py                                  # qwen2.5:7b, both protocols
    python scripts/babs_eval.py --model qwen2.5:7b --protocol json
    python scripts/babs_eval.py --tasks protein,primers --out results.json

SAFETY: everything runs in a throwaway sandbox (per
.claude/skills/verifier-splicecraft.md): the data dir AND $HOME point into a
temp dir before SpliceCraft is imported, so no real library, setting, log or
~/babs memory/corpus can be read or written — the run refuses to start if a
babs repo is still reachable. The plasmid is synthetic and built offline.

It is SLOW on a CPU-only machine (minutes per task): leave it running.
Nothing here is imported by the app; `tests/test_babs_eval_harness.py` keeps
the graders and the fixture honest.
"""
from __future__ import annotations

import argparse
import asyncio
import json
import os
import random
import re
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PLASMID_NAME = "pBabsEval"
PLASMID_LEN = 2400


# ── the fixture: a synthetic plasmid with realistic names ─────────────────────

def _cds(rng: random.Random, n_codons: int) -> str:
    sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
             if a + b + c not in ("TAA", "TAG", "TGA")]
    return "ATG" + "".join(rng.choice(sense) for _ in range(n_codons)) + "TAA"


def build_fixture(seed: int = 7):
    """(SeqRecord, layout). AmpR on the + strand, CmR on the − strand, plus an
    ori, a terminator, lacZα and an MCS. Deterministic for a given seed."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rng = random.Random(seed)
    seq = [rng.choice("ACGT") for _ in range(PLASMID_LEN)]

    def put(at: int, bases: str) -> None:
        seq[at:at + len(bases)] = list(bases)
    amp = _cds(rng, 180)                 # 546 bp  [300, 846)  +
    put(300, amp)
    cm = _cds(rng, 130)                  # 396 bp  [1100, 1496) −
    put(1100, str(Seq(cm).reverse_complement()))
    lacz = _cds(rng, 40)                 # 126 bp  [1700, 1826) +
    put(1700, lacz)
    rec = SeqRecord(Seq("".join(seq)), id=PLASMID_NAME, name=PLASMID_NAME,
                    description="Babs test bench plasmid")
    rec.annotations["topology"] = "circular"
    rec.annotations["molecule_type"] = "DNA"
    layout = {"AmpR": (300, 846, 1), "CmR": (1100, 1496, -1),
              "lacZα": (1700, 1826, 1)}
    rec.features = [
        SeqFeature(FeatureLocation(0, 180, strand=1), type="rep_origin",
                   qualifiers={"label": ["ori"]}),
        SeqFeature(FeatureLocation(300, 846, strand=1), type="CDS",
                   qualifiers={"label": ["AmpR"], "gene": ["bla"]}),
        SeqFeature(FeatureLocation(900, 980, strand=1), type="terminator",
                   qualifiers={"label": ["T1 terminator"]}),
        SeqFeature(FeatureLocation(1100, 1496, strand=-1), type="CDS",
                   qualifiers={"label": ["CmR"], "gene": ["cat"]}),
        SeqFeature(FeatureLocation(1650, 1690, strand=1), type="misc_feature",
                   qualifiers={"label": ["MCS"]}),
        SeqFeature(FeatureLocation(1700, 1826, strand=1), type="CDS",
                   qualifiers={"label": ["lacZα"]}),
    ]
    return rec, layout


# ── ground truth, from the plasmid as loaded ──────────────────────────────────

@dataclass
class Truth:
    name: str
    length: int
    rows: dict                        # label → features row
    unique: set                       # enzymes cutting exactly once
    multi: set                        # enzymes cutting 2+ times
    unique_in_amp: set                # unique cutters whose cut falls in AmpR
    amp_protein10: str
    sequence: str
    primer_tm: float
    digest_pair: "tuple[str, str]"    # two once-cutters well apart, for `digest`
    digest_sizes: "list[int]"         # the fragments they give, largest first


PRIMER_FOR_TM = "ACGTACGTTTGGCCAAGTGA"


def compute_truth(sc, app) -> Truth:
    rows = {r["label"]: r for r in sc._h_features(app, {})["features"]}
    scan = sc._h_list_restriction_sites(app, {"respect_active_collection": False})
    counts: dict = {}
    cuts: dict = {}
    for s in scan["sites"]:
        counts[s["enzyme"]] = counts.get(s["enzyme"], 0) + 1
        cuts[s["enzyme"]] = s["cut_bp"]
    unique = {e for e, n in counts.items() if n == 1}
    amp = rows["AmpR"]
    n = len(app._current_record.seq)
    in_amp = {e for e in unique
              if sc._bp_in_span(cuts[e], amp["start"], amp["end"], n)}
    aa = sc._h_find_feature(app, {"name": "AmpR"})["matches"][0]["translation"]
    seq = str(app._current_record.seq).upper()
    # Two once-cutters at least 300 bp apart both ways, picked by name so the
    # choice is stable; the sizes come from SpliceCraft's own digest.
    pair = next((a, b) for a in sorted(unique) for b in sorted(unique)
                if a < b and 300 <= (cuts[b] - cuts[a]) % n <= n - 300)
    digest = sc._AGENT_HANDLERS["digest"][0](
        app, {"sequence": seq, "enzymes": list(pair), "circular": True})
    sizes = sorted((int(f["length"]) for f in digest["fragments"]), reverse=True)
    return Truth(name=PLASMID_NAME, length=n, rows=rows, unique=unique,
                 multi={e for e, k in counts.items() if k > 1},
                 unique_in_amp=in_amp, amp_protein10=aa[:10], sequence=seq,
                 primer_tm=float(sc._primer_tm(PRIMER_FOR_TM)),
                 digest_pair=pair, digest_sizes=sizes)


# ── graders: (answer, record-after, truth) → (passed, detail) ─────────────────

def _ints(text: str) -> "set[int]":
    return {int(x.replace(",", "")) for x in re.findall(r"\d[\d,]*", text or "")}


def _named(text: str, enzymes) -> "set[str]":
    return {e for e in enzymes
            if re.search(rf"(?<![A-Za-z0-9]){re.escape(e)}(?![A-Za-z0-9])",
                         text or "", re.IGNORECASE)}


def grade_overview(answer, _rec, t: Truth):
    a = (answer or "").lower()
    ok_len = t.length in _ints(answer)
    ok = t.name.lower() in a and ok_len and "circular" in a
    return ok, f"name={t.name.lower() in a} length={ok_len} circular={'circular' in a}"


def grade_find_coords(answer, _rec, t: Truth):
    r = t.rows["CmR"]
    nums = _ints(answer)
    start_ok = r["start"] in nums or r["start"] + 1 in nums
    end_ok = r["end"] in nums
    minus = bool(re.search(r"minus|reverse|complement|negative|bottom|-1|−", answer or "",
                           re.IGNORECASE))
    return (start_ok and end_ok and minus,
            f"start={start_ok} end={end_ok} strand={minus} (truth {r['start']}..{r['end']} −)")


def grade_unique_cutters(answer, _rec, t: Truth):
    # The prompt asks for "at least three", so three right and none wrong passes.
    right, wrong = _named(answer, t.unique), _named(answer, t.multi)
    return (len(right) >= 3 and not wrong,
            f"{len(right)} correct, {len(wrong)} not unique {sorted(wrong)[:4]}")


def grade_cutters_in_gene(answer, _rec, t: Truth):
    # "Which enzymes …" asks for the set, so naming one of thirteen is not an
    # answer: pass needs at least half of them (all, if there are ≤ 2) and none
    # wrong. The first bench run passed 1-of-13 under a looser rule.
    right = _named(answer, t.unique_in_amp)
    wrong = _named(answer, (t.unique - t.unique_in_amp) | t.multi)
    need = len(t.unique_in_amp) if len(t.unique_in_amp) <= 2 else (len(t.unique_in_amp) + 1) // 2
    return (len(right) >= need and not wrong,
            f"{len(right)}/{len(t.unique_in_amp)} correct (need {need}), "
            f"{len(wrong)} wrong {sorted(wrong)[:4]}")


def grade_protein(answer, _rec, t: Truth):
    flat = re.sub(r"[^A-Za-z]", "", answer or "").upper()
    return t.amp_protein10 in flat, f"want {t.amp_protein10}"


def grade_primers(answer, _rec, t: Truth):
    amp = t.rows["AmpR"]
    s, e = amp["start"], amp["end"]
    top = t.sequence

    def rc(x):
        return x.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    oligos = re.findall(r"[ACGT]{16,}", (answer or "").upper())
    fwd = any(top[s:s + len(o)] == o for o in oligos)
    rev = any(rc(top[e - len(o):e]) == o for o in oligos)
    return fwd and rev, f"forward={fwd} reverse={rev} ({len(oligos)} oligo(s) quoted)"


def grade_add_feature(_answer, after, _t: Truth):
    hits = []
    for f in getattr(after.get("record"), "features", []) or []:
        lab = (f.qualifiers.get("label") or [""])[0] if f.qualifiers else ""
        if str(lab).strip().lower() == "test site":
            hits.append((int(f.location.start), int(f.location.end)))
    if (99, 200) in hits:
        return True, "added at [99, 200) = bases 100–200"
    return False, f"'test site' features: {hits or 'none'} (want [99, 200))"


def grade_primer_tm(answer, _rec, t: Truth):
    nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", answer or "")]
    ok = any(abs(x - t.primer_tm) <= 1.0 for x in nums)
    return ok, f"want {t.primer_tm:.1f} °C ±1"


def grade_digest(answer, _after, t: Truth):
    nums = _ints(answer)
    found = [s for s in t.digest_sizes if s in nums]
    return (len(found) == len(t.digest_sizes),
            f"sizes {t.digest_sizes} with {' + '.join(t.digest_pair)}; found {found}")


def grade_save_primers(_answer, after, t: Truth):
    def key(name) -> str:
        return re.sub(r"[^a-z0-9]", "", str(name or "").lower())
    saved = {key(p.get("name")): str(p.get("sequence") or "").upper()
             for p in after.get("primers") or [] if isinstance(p, dict)}
    fwd, rev = saved.get("cmrf"), saved.get("cmrr")
    if not fwd or not rev:
        return False, f"saved primers: {sorted(saved)[:6]} (want CmR-F and CmR-R)"
    cm = t.rows["CmR"]
    s, e = cm["start"], cm["end"]

    def rc(x):
        return x.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def left(o):                      # anneals at CmR's leftmost base (top strand)
        return len(o) >= 16 and t.sequence[s:s + len(o)] == o

    def right(o):                     # anneals at CmR's rightmost base (bottom strand)
        return len(o) >= 16 and rc(t.sequence[e - len(o):e]) == o
    # Which one a person calls "F" for a − strand gene varies; either way the
    # pair must be the two ends of CmR.
    ok = (left(fwd) and right(rev)) or (left(rev) and right(fwd))
    return ok, f"CmR-F {len(fwd)} nt, CmR-R {len(rev)} nt — bind CmR's two ends: {ok}"


@dataclass
class Task:
    key: str
    prompt: str
    autonomy: str
    grade: object


TASKS = [
    Task("overview", "What plasmid is open in SpliceCraft right now, how long is it, "
                     "and is it circular or linear?", "readonly", grade_overview),
    Task("find_coords", "Where exactly is the CmR gene on this plasmid? Give its start "
                        "and end positions and which strand it is on.", "readonly",
         grade_find_coords),
    Task("unique_cutters", "Which restriction enzymes cut this plasmid exactly once? "
                           "Name at least three.", "readonly", grade_unique_cutters),
    Task("cutters_in_gene", "Which enzymes cut this plasmid only once, with that single "
                            "cut inside the AmpR gene?", "readonly", grade_cutters_in_gene),
    Task("protein", "What are the first 10 amino acids of the protein encoded by AmpR?",
         "readonly", grade_protein),
    Task("primers", "Design PCR primers to amplify the whole AmpR gene. Give me both "
                    "primer sequences.", "readonly", grade_primers),
    Task("add_feature", "Add a misc_feature called 'test site' covering bases 100 to "
                        "200 of this plasmid.", "auto", grade_add_feature),
    Task("primer_tm", f"What is the melting temperature of the primer {PRIMER_FOR_TM}?",
         "readonly", grade_primer_tm),
    # Added after the first round: a 7B summed cut positions itself instead of
    # digesting, and designed primers it was asked to SAVE and stopped there.
    Task("digest", "If I cut this plasmid with {e1} and {e2} together, what fragment "
                   "sizes do I get?", "readonly", grade_digest),
    Task("save_primers", "Design PCR primers that amplify the CmR gene and save both "
                         "to my primer library, named CmR-F and CmR-R.", "auto",
         grade_save_primers),
]


@dataclass
class Result:
    task: str
    model: str
    protocol: str
    passed: bool
    detail: str
    seconds: float
    steps: "list[str]" = field(default_factory=list)
    answer: str = ""
    error: str = ""
    step_errors: "list[str]" = field(default_factory=list)
    calls: "list[str]" = field(default_factory=list)   # each call with its arguments
    notes: "list[str]" = field(default_factory=list)   # ⚠ warnings + announce nudges


# ── the runner ────────────────────────────────────────────────────────────────

async def run_task(sc, task: Task, model: str, protocol: str, *, timeout_s: float,
                   num_ctx: int) -> Result:
    from textual.widgets import Input
    rec, _layout = build_fixture()
    app = sc.PlasmidApp()
    steps: "list[str]" = []
    step_errors: "list[str]" = []
    calls: "list[str]" = []
    notes: "list[str]" = []
    async with app.run_test(size=(200, 50)) as pilot:
        for _ in range(6):
            await pilot.pause()
        while len(app.screen_stack) > 1:           # What's New on a fresh data dir
            app.pop_screen()
            await pilot.pause()
        app._apply_record(rec)
        for _ in range(6):
            await pilot.pause()
        # Every task starts from the same empty primer library. The sandbox is
        # shared by the whole run, and a save is refused (409) when the SAME
        # sequence is already filed — so without this, the second protocol's
        # save_primers found the first one's primers, "passed" on them, and
        # measured nothing.
        sc._save_primers([])
        assert sc._load_primers() == [], "could not reset the primer library"
        truth = compute_truth(sc, app)
        app.action_open_babs()
        for _ in range(6):
            await pilot.pause()
        scr = app.screen
        scr._agent_enabled = True
        scr._autonomy = task.autonomy
        scr._model = scr._agent_model = model
        scr._agent_protocol_pref = protocol
        scr._grounded = False
        scr._num_ctx = num_ctx
        real_run = scr._run_tool_call

        def spy(tc):                               # observe, never alter
            label = sc._babs_tool_call_label(tc)
            steps.append(label)
            calls.append(f"{sc._babs_tool_call_name(tc)} "
                         f"{json.dumps(sc._babs_tool_call_args(tc), default=str)[:240]}")
            result = real_run(tc)
            err = result.get("error") if isinstance(result, dict) else None
            if not err and isinstance(result, dict) and isinstance(result.get("result"), dict):
                err = result["result"].get("error")
            if err:
                step_errors.append(f"{label}: {str(err)[:120]}")
            return result
        scr._run_tool_call = spy
        real_note = scr._sys_note

        def note_spy(text):                        # what the loop told the user
            if "⚠" in str(text) or "go ahead" in str(text):   # warnings + the nudge
                notes.append(str(text)[:240])
            return real_note(text)
        scr._sys_note = note_spy
        errors: "list[str]" = []
        real_finalize = scr._finalize

        def finalize_spy(query, display, answer, err, stopped):
            if err:
                errors.append(str(err))
            return real_finalize(query, display, answer, err, stopped)
        scr._finalize = finalize_spy
        t0 = time.monotonic()
        scr.query_one("#babs-input", Input).value = task.prompt.format(
            e1=truth.digest_pair[0], e2=truth.digest_pair[1])
        scr._submit_current()
        while scr._generating and time.monotonic() - t0 < timeout_s:
            await asyncio.sleep(0.5)
            await pilot.pause()
        timed_out = scr._generating
        if timed_out:
            scr._stop_generation()
            for _ in range(60):
                await asyncio.sleep(0.5)
                await pilot.pause()
                if not scr._generating:
                    break
        seconds = time.monotonic() - t0
        answer = (scr._history[-1]["content"]
                  if scr._history and scr._history[-1].get("role") == "assistant"
                  else "")
        after = {"record": app._current_record, "primers": sc._load_primers()}
        passed, detail = task.grade(answer, after, truth)
        error = "timed out" if timed_out else "; ".join(errors)[:300]
        return Result(task.key, model, protocol, bool(passed), detail, round(seconds, 1),
                      steps, answer[:2000], error, step_errors, calls, notes)


def _num_ctx_for(sc, model: str) -> int:
    B = sc._babs
    try:
        real = B.context_length_from_show(B.show(model))
    except B.BabsError:
        real = None
    return max(B.DEFAULT_NUM_CTX, min(int(real), B.MAX_AUTO_NUM_CTX)) if real else B.DEFAULT_NUM_CTX


async def run_suite(sc, models, protocols, tasks, *, timeout_s: float, log=print):
    results = []
    for model in models:
        num_ctx = _num_ctx_for(sc, model)
        for protocol in protocols:
            for task in tasks:
                log(f"[{model} · {protocol}] {task.key} …")
                r = await run_task(sc, task, model, protocol, timeout_s=timeout_s,
                                   num_ctx=num_ctx)
                mark = "PASS" if r.passed else "FAIL"
                log(f"    {mark} in {r.seconds:.0f}s · {len(r.steps)} step(s) "
                    f"{r.steps[:6]} · {r.detail}{' · ' + r.error if r.error else ''}")
                for e in r.step_errors[:4]:
                    log(f"      step failed — {e}")
                for n in r.notes[:2]:
                    log(f"      noted — {n}")
                results.append(r)
    return results


def _summary(results) -> str:
    lines = []
    configs = sorted({(r.model, r.protocol) for r in results})
    for model, protocol in configs:
        rs = [r for r in results if r.model == model and r.protocol == protocol]
        ok = sum(r.passed for r in rs)
        secs = sum(r.seconds for r in rs)
        lines.append(f"{model} · {protocol}: {ok}/{len(rs)} passed · "
                     f"{secs / 60:.1f} min total · "
                     f"{sum(len(r.steps) for r in rs) / max(1, len(rs)):.1f} steps/task")
    return "\n".join(lines)


def _sandbox() -> str:
    """Point the data dir, $HOME and the log into a temp dir BEFORE SpliceCraft
    is imported (its data dir is fixed at import time)."""
    root = tempfile.mkdtemp(prefix="sc-verify-babs-eval-")
    home = Path(root) / "home"
    home.mkdir()
    os.environ["XDG_DATA_HOME"] = str(Path(root) / "data")
    os.environ["HOME"] = str(home)                 # ~/babs must not resolve
    for var in ("SPLICECRAFT_DATA_DIR", "SPLICECRAFT_BABS_HOME", "BABS_HOME"):
        os.environ.pop(var, None)
    os.environ["SPLICECRAFT_LOG"] = str(Path(root) / "splicecraft.log")
    os.environ.setdefault("SPLICECRAFT_SKIP_LOCK", "1")
    return root


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--model", action="append",
                    help="model to test (repeatable; default qwen2.5:7b)")
    ap.add_argument("--protocol", choices=("native", "json", "both"), default="both")
    ap.add_argument("--tasks", default="",
                    help="comma-separated task keys (default: all): "
                         + ",".join(t.key for t in TASKS))
    ap.add_argument("--timeout", type=float, default=1800,
                    help="seconds per task before it is stopped (default 1800)")
    ap.add_argument("--out", default="", help="write the results as JSON here")
    args = ap.parse_args(argv)
    root = _sandbox()
    sys.path.insert(0, str(ROOT))
    import splicecraft as sc
    assert "sc-verify-babs-eval-" in str(sc._state._DATA_DIR), sc._state._DATA_DIR
    sc._authorize_writes_for_sandbox(sc._state._DATA_DIR)
    if sc._learn_resolve_babs_home() is not None:
        print("refusing to run: a real babs repo is reachable from the sandbox")
        return 2
    for flag in ("_skip_seed", "_skip_snapshot", "_skip_update_check"):
        setattr(sc.PlasmidApp, flag, True)
    sc.PlasmidApp._preload_demo_record = None
    want = {k.strip() for k in args.tasks.split(",") if k.strip()}
    unknown = want - {t.key for t in TASKS}
    if unknown:
        print(f"unknown task(s): {', '.join(sorted(unknown))}")
        return 2
    tasks = [t for t in TASKS if not want or t.key in want]
    models = args.model or ["qwen2.5:7b"]
    protocols = ["native", "json"] if args.protocol == "both" else [args.protocol]
    # Preflight: without this, a stopped Ollama or an unpulled model makes
    # every task "fail" — which reads as a bad model, not a bad setup.
    B = sc._babs
    if not B.ping():
        print("Ollama isn't reachable — start it (`ollama serve`) and re-run.")
        return 2
    try:
        installed = B.list_installed()
    except B.BabsError as exc:
        print(f"couldn't list Ollama models: {exc}")
        return 2
    missing = [m for m in models if not B.model_is_installed(m, installed)]
    if missing:
        print(f"not pulled: {', '.join(missing)} — `ollama pull <name>` first.")
        return 2
    print(f"sandbox: {root}")
    results = asyncio.run(run_suite(sc, models, protocols, tasks,
                                    timeout_s=args.timeout))
    print("\n" + _summary(results))
    out = Path(args.out) if args.out else Path(root) / "babs-eval-results.json"
    out.write_text(json.dumps([r.__dict__ for r in results], indent=2,
                              ensure_ascii=False), encoding="utf-8")
    print(f"results: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
