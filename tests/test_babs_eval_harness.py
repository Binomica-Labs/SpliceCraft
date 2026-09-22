"""The Babs live test bench (`scripts/babs_eval.py`) must not rot.

The bench runs real models and is far too slow for the suite, but its parts
are deterministic: the fixture plasmid, the ground truth taken from it, and the
graders that score answers. A grader that passes a wrong answer (or fails a
right one) would make every bench run lie, so each is pinned here both ways,
and one task runs end to end through the REAL Babs loop with a scripted model.
"""
from __future__ import annotations

import asyncio
import importlib.util
import json
import sys
from pathlib import Path

import pytest

import splicecraft as sc
import splicecraft_babs as B

_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _ROOT / "scripts" / "babs_eval.py"
if not _SCRIPT.is_file():          # the sdist ships tests/ but not scripts/
    pytest.skip("scripts/babs_eval.py not present (source distribution)",
                allow_module_level=True)
_spec = importlib.util.spec_from_file_location("babs_eval", _SCRIPT)
assert _spec is not None and _spec.loader is not None
bench = importlib.util.module_from_spec(_spec)
sys.modules["babs_eval"] = bench     # dataclasses resolve their module by name
_spec.loader.exec_module(bench)

_TERM = (200, 50)


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.fixture
async def truth():
    rec, _layout = bench.build_fixture()
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        app._apply_record(rec)
        for _ in range(6):
            await pilot.pause()
        yield bench.compute_truth(sc, app)


async def test_fixture_supports_every_task(truth):
    assert truth.length == bench.PLASMID_LEN and truth.name == bench.PLASMID_NAME
    assert {"AmpR", "CmR", "lacZα"} <= set(truth.rows)
    assert truth.rows["CmR"]["strand"] == -1
    assert len(truth.unique) >= 3, "unique_cutters needs at least three"
    assert truth.unique_in_amp, "cutters_in_gene needs a unique cutter inside AmpR"
    assert truth.amp_protein10.startswith("M") and len(truth.amp_protein10) == 10
    assert 40.0 < truth.primer_tm < 80.0


async def test_graders_pass_right_and_fail_wrong(truth):
    t = truth
    cm = t.rows["CmR"]
    amp = t.rows["AmpR"]
    good_u, bad_u = sorted(t.unique)[:3], sorted(t.multi)[:1]
    inside = sorted(t.unique_in_amp)[:1]
    outside = sorted(t.unique - t.unique_in_amp)[:1]
    fwd = t.sequence[amp["start"]: amp["start"] + 22]
    rev = _rc(t.sequence[amp["end"] - 22: amp["end"]])
    cases = [
        (bench.grade_overview, f"It is {t.name}, {t.length:,} bp and circular.",
         f"It is {t.name}, {t.length + 1} bp and circular."),
        (bench.grade_find_coords, f"CmR runs {cm['start'] + 1}..{cm['end']} on the minus strand.",
         f"CmR runs {cm['start'] + 1}..{cm['end']} on the top strand."),
        (bench.grade_unique_cutters, "Single cutters: " + ", ".join(good_u),
         "Single cutters: " + ", ".join(good_u + bad_u)),
        (bench.grade_protein, "It starts " + "-".join(t.amp_protein10),
         "It starts " + t.amp_protein10[::-1]),
        (bench.grade_primers, f"Forward {fwd}\nReverse {rev}",
         f"Forward {fwd}\nReverse {_rc(rev)}"),
        (bench.grade_primer_tm, f"The Tm is about {t.primer_tm + 0.4:.1f} °C.",
         f"The Tm is about {t.primer_tm + 6:.1f} °C."),
    ]
    if inside:
        half = sorted(t.unique_in_amp)[: (len(t.unique_in_amp) + 1) // 2]
        cases.append((bench.grade_cutters_in_gene, "These do: " + ", ".join(half),
                      f"{inside[0]} and {outside[0] if outside else bad_u[0]}."))
        if len(t.unique_in_amp) > 2:        # one of many is not an answer
            assert not bench.grade_cutters_in_gene(f"Only {inside[0]}.", None, t)[0]
    for grade, right, wrong in cases:
        assert grade(right, None, t)[0] is True, (grade.__name__, right)
        assert grade(wrong, None, t)[0] is False, (grade.__name__, wrong)


def test_add_feature_grader_is_exact_about_the_one_based_range():
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    def rec_with(a, b):
        r = SeqRecord(Seq("A" * 300))
        r.features = [SeqFeature(FeatureLocation(a, b), type="misc_feature",
                                 qualifiers={"label": ["Test Site"]})]
        return r
    assert bench.grade_add_feature("", rec_with(99, 200), None)[0] is True
    for off_by_one in ((100, 200), (100, 201), (99, 201)):
        assert bench.grade_add_feature("", rec_with(*off_by_one), None)[0] is False


def test_task_keys_are_unique_and_graded():
    keys = [t.key for t in bench.TASKS]
    assert len(keys) == len(set(keys)) == 8
    assert all(callable(t.grade) and t.autonomy in ("readonly", "auto", "ask")
               for t in bench.TASKS)


async def test_one_task_end_to_end_with_a_scripted_model(monkeypatch):
    """The whole bench path — sandboxed app, the real Babs loop, the step spy,
    grading — with the model scripted: look AmpR up, then quote its protein."""
    monkeypatch.setattr(B, "list_installed", lambda **k: [{"name": "fake:1b"}])
    monkeypatch.setattr(B, "model_is_installed", lambda name, *a, **k: True)
    monkeypatch.setattr(B, "show", lambda *a, **k: {})
    monkeypatch.setattr(sc, "_babs_memory_index", lambda: "")
    calls = {"n": 0}

    def fake_chat(model, messages, **kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            yield {"content": "", "thinking": "", "done": True, "done_reason": "stop",
                   "error": None, "tool_calls": [{"function": {
                       "name": "splicecraft_find_feature", "arguments": {"name": "AmpR"}}}]}
            return
        found = json.loads(next(m for m in reversed(messages)
                                if m.get("role") == "tool")["content"])
        aa = found["result"]["matches"][0]["translation"][:10]   # dispatch wraps payloads
        yield {"content": f"The protein begins {aa}.", "thinking": "", "tool_calls": [],
               "done": True, "done_reason": "stop", "error": None}
    monkeypatch.setattr(B, "chat_stream", fake_chat)
    task = next(t for t in bench.TASKS if t.key == "protein")
    r = await bench.run_task(sc, task, "fake:1b", "native", timeout_s=60, num_ctx=16384)
    assert r.passed, r
    assert r.steps == ["find AmpR"] and r.error == ""
    assert calls["n"] == 2
    asyncio.get_running_loop()      # still inside the test's loop — nothing leaked
