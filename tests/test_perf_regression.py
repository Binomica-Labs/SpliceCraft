"""test_perf_regression — best-of-N regression gates against perf-baseline.json.

Complements tests/test_performance.py (which is a single-shot smoke). The
test here is `@pytest.mark.slow` so the default `pytest -n auto -q` run
INCLUDES it; CI's dedicated perf job runs it serially to avoid the
parallel-load flake the smoke perf tests have hit historically (see
the 2026-05-17 flake on test_build_large_under_200ms).

Why best-of-N rather than median/p99:
    Under pytest-xdist load + GC noise, even the median sample drifts
    by 2-3×. The MINIMUM sample is what the code is capable of when
    nothing else is competing; if best-of-N is over budget, the code
    has genuinely regressed. If best-of-N is under budget but variance
    is high, that's environmental noise — not a code issue.

Budget changes:
    Tighten freely when code gets faster.
    Loosen ONLY with a written rationale in the PR description.
"""
from __future__ import annotations

import json
import random
import time
from pathlib import Path

import pytest

import splicecraft as sc


_BASELINE_PATH = Path(__file__).parent / "perf-baseline.json"


@pytest.fixture(scope="module")
def baseline() -> dict:
    """Load perf-baseline.json once per module run.

    Schema-version check sits inside the fixture so a malformed /
    mismatched baseline fails every test in the module with a clear
    `pytest.fail` rather than `KeyError` somewhere later."""
    if not _BASELINE_PATH.exists():
        pytest.fail(
            f"Perf baseline missing at {_BASELINE_PATH}. "
            f"Restore from git or regenerate via the perf probe script."
        )
    try:
        data = json.loads(_BASELINE_PATH.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        pytest.fail(f"Perf baseline is not valid JSON: {exc}")
    if data.get("_schema_version") != 1:
        pytest.fail(
            f"Perf baseline schema version "
            f"{data.get('_schema_version')!r} unsupported "
            f"(this test expects 1)."
        )
    ops = data.get("operations")
    if not isinstance(ops, dict):
        pytest.fail(
            "Perf baseline missing 'operations' dict — "
            "every regression test needs an entry."
        )
    return data


def _measure(op_name: str, baseline_data: dict, work) -> float:
    """Run `work()` `samples` times after `warmup_iterations` warmups
    and return the best wall-clock duration in milliseconds.

    Raises pytest.fail if `op_name` is missing from the baseline so a
    new test is forced to land its baseline entry in the same PR.
    """
    ops = baseline_data["operations"]
    if op_name not in ops:
        pytest.fail(
            f"Operation {op_name!r} not declared in perf-baseline.json. "
            f"Add a `budget_ms` entry before the test can run."
        )
    spec = ops[op_name]
    samples       = int(spec.get("samples", 5))
    warmup_iters  = int(spec.get("warmup_iterations", 2))
    if samples < 1:
        pytest.fail(f"{op_name}: `samples` must be >= 1, got {samples}")
    for _ in range(warmup_iters):
        work()
    durations_ms = []
    for _ in range(samples):
        t0 = time.perf_counter()
        work()
        durations_ms.append((time.perf_counter() - t0) * 1000)
    return min(durations_ms)


def _assert_within_budget(op_name: str, baseline_data: dict,
                            best_ms: float) -> None:
    budget = float(baseline_data["operations"][op_name]["budget_ms"])
    desc   = baseline_data["operations"][op_name].get("description", "")
    assert best_ms < budget, (
        f"\n  Operation : {op_name}"
        f"\n  Detail    : {desc}"
        f"\n  Budget    : {budget:.2f} ms"
        f"\n  Best-of-N : {best_ms:.2f} ms  ← OVER BUDGET"
        f"\n  Action    : either fix the regression or update "
        f"perf-baseline.json with a written rationale."
    )


# ── Fixtures: deterministic synthetic sequences ────────────────────────────


@pytest.fixture(scope="module")
def random_puc19() -> str:
    """A pUC19-size (2686 bp) ACGT sequence with a fixed RNG seed so
    the runtime is reproducible across machines (no real biology — the
    test cares about scanning cost, not biological correctness)."""
    rng = random.Random(0xC0FFEE)
    return "".join(rng.choice("ACGT") for _ in range(2686))


@pytest.fixture(scope="module")
def random_10kb() -> str:
    rng = random.Random(0xBEEF)
    return "".join(rng.choice("ACGT") for _ in range(10_000))


def _synth_record_module(n_bp: int, n_feats: int):
    """Construct a minimal SeqRecord with random ACGT + `n_feats`
    spread-out misc_features. Shares the helper signature with
    test_performance.py for parity."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rng = random.Random(0xFEED)
    seq = "".join(rng.choice("ACGT") for _ in range(n_bp))
    feats = []
    if n_feats > 0:
        stride = max(1, n_bp // (n_feats + 1))
        feat_len = max(20, stride // 3)
        for i in range(n_feats):
            s = (i + 1) * stride
            e = min(n_bp, s + feat_len)
            if e <= s:
                continue
            feats.append(SeqFeature(
                FeatureLocation(s, e, strand=1),
                type="misc_feature",
                qualifiers={"label": [f"feat{i:02d}"]},
            ))
    rec = SeqRecord(Seq(seq), id="perf", name="perf",
                      description="synthetic perf record")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    rec.features = feats
    return rec


# ── Regression tests ──────────────────────────────────────────────────────────


@pytest.mark.slow
class TestScanRestrictionSitesPerf:
    """Restriction-site scanning is the most-touched hot path — every
    record load + every settings toggle re-runs it on the active seq."""

    def test_scan_puc19(self, baseline, random_puc19):
        def _work():
            sc._scan_restriction_sites(random_puc19)
        best = _measure("scan_restriction_sites_puc19", baseline, _work)
        _assert_within_budget("scan_restriction_sites_puc19",
                                baseline, best)

    def test_scan_10kb(self, baseline, random_10kb):
        def _work():
            sc._scan_restriction_sites(random_10kb)
        best = _measure("scan_restriction_sites_10kb", baseline, _work)
        _assert_within_budget("scan_restriction_sites_10kb",
                                baseline, best)


@pytest.mark.slow
class TestBiologyPrimitivesPerf:
    """The biology primitives below `_scan_restriction_sites`: every
    palindromic enzyme uses `_rc`; every IUPAC pattern is compiled
    once and cached."""

    def test_rc_10kb(self, baseline, random_10kb):
        def _work():
            sc._rc(random_10kb)
        best = _measure("rc_10kb", baseline, _work)
        _assert_within_budget("rc_10kb", baseline, best)

    def test_iupac_pattern_warm(self, baseline):
        """200 cache hits on _iupac_pattern. Cache miss is irrelevant —
        the cache is warm in steady state."""
        # Warm a fixed set of patterns first.
        patterns = ["GAATTC", "GGATCC", "AAGCTT", "GGTACC",
                    "GAGCTC", "GTCGAC", "CTGCAG", "GCATGC"]
        for p in patterns:
            sc._iupac_pattern(p)

        def _work():
            for _ in range(200):
                for p in patterns:
                    sc._iupac_pattern(p)

        best = _measure("iupac_pattern_warm_200", baseline, _work)
        _assert_within_budget("iupac_pattern_warm_200", baseline, best)


@pytest.mark.slow
class TestBuildSeqTextPerf:
    """_build_seq_text is the per-chunk seq-panel render. Hot path on
    every cursor move."""

    def test_puc19(self, baseline):
        rec = _synth_record_module(2686, 10)
        pm = sc.PlasmidMap()
        pm.load_record(rec)
        seq_str, feats = str(rec.seq), pm._feats

        def _work():
            sc._BUILD_SEQ_CACHE.clear()
            sc._build_seq_text(seq_str, feats, line_width=127)

        best = _measure("build_seq_text_puc19", baseline, _work)
        _assert_within_budget("build_seq_text_puc19", baseline, best)

    def test_20kb(self, baseline):
        rec = _synth_record_module(20_000, 80)
        pm = sc.PlasmidMap()
        pm.load_record(rec)
        seq_str, feats = str(rec.seq), pm._feats

        def _work():
            sc._BUILD_SEQ_CACHE.clear()
            sc._build_seq_text(seq_str, feats, line_width=127)

        best = _measure("build_seq_text_20kb", baseline, _work)
        _assert_within_budget("build_seq_text_20kb", baseline, best)


# ── Self-tests on the harness ──────────────────────────────────────────────


class TestPerfHarness:
    """Tests that the perf-baseline.json schema + the measurement
    helper itself do what they claim. These run without `slow` so a
    broken harness shows up in the fast CI matrix, not just in the
    perf job."""

    def test_baseline_schema_v1(self, baseline):
        """Schema-version check; ensures fixture loads + every op
        carries the minimal required keys."""
        assert baseline["_schema_version"] == 1
        for name, spec in baseline["operations"].items():
            assert "budget_ms" in spec, f"{name} missing budget_ms"
            assert spec["budget_ms"] > 0, (
                f"{name} budget_ms must be positive"
            )
            assert "description" in spec, f"{name} missing description"

    def test_baseline_unknown_op_fails_loudly(self, baseline):
        """An op not in the baseline must surface as a pytest.fail at
        measurement time — the new test author needs a forcing function
        to add their budget entry to perf-baseline.json."""
        try:
            _measure("not_a_real_op", baseline, lambda: None)
        except pytest.fail.Exception as exc:
            assert "not declared" in str(exc)
        else:
            pytest.fail(
                "Expected _measure to fail on unknown operation, "
                "but it returned."
            )

    def test_best_of_n_returns_minimum(self, baseline):
        """Sanity-check the harness: with a sleep-based work function,
        the returned best-of-N must equal (within tolerance) the
        sleeping duration, not the sum or the mean."""
        # Insert a temp operation into the in-memory baseline so the
        # harness has something to look up. We restore the original
        # ops dict in `finally` so concurrent tests aren't affected.
        original_ops = baseline["operations"]
        original_ops["_harness_self_test"] = {
            "description": "harness self-test (sleep 5 ms)",
            "budget_ms": 100,
            "samples": 5,
            "warmup_iterations": 0,
        }
        try:
            def _work():
                time.sleep(0.005)
            best = _measure("_harness_self_test", baseline, _work)
            # best-of-N must be close to 5 ms, NOT 25 ms (sum) or
            # ~5 ms (mean). Tolerance: 0 < best < 25 ms.
            assert 0 < best < 25, (
                f"best-of-N harness regressed: returned {best:.2f} ms "
                f"for 5x5ms sleeps — should be close to 5 ms"
            )
        finally:
            original_ops.pop("_harness_self_test", None)


# ── 2026-11-09: the gate covered 6 operations ([INV-234]) ───────────────────
# The digest, the codon optimizer, the scanners, the GenBank reader/writer and
# the data layer were all ungated, so any of them could regress silently. Each
# test below CLEARS the LRU it would otherwise hit: with a warm
# `_ENZYME_CUTS_CACHE` the digest measured 0.01 ms, i.e. it measured nothing.


@pytest.fixture(scope="module")
def random_5kb() -> str:
    rng = random.Random(0x5EED)
    return "".join(rng.choice("ACGT") for _ in range(5_000))


@pytest.fixture(scope="module")
def random_50kb() -> str:
    rng = random.Random(0x5A5A)
    return "".join(rng.choice("ACGT") for _ in range(50_000))


@pytest.fixture(scope="module")
def flat_codon_table() -> dict:
    """A full 64-codon usage table in the registry's own
    ``{codon: (aa, count)}`` shape. A bare ``{codon: count}`` unpacks wrong
    inside `_codon_build_aa_map` — the optimizer would raise instead of being
    timed."""
    import splicecraft_codon as _codon
    return {c: (aa, 100.0) for c, aa in _codon._CODON_GENETIC_CODE.items()}


@pytest.fixture(scope="module")
def protein_500() -> str:
    rng = random.Random(11)
    return "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY")
                         for _ in range(499))


@pytest.mark.slow
class TestEnginePerf:
    """The construction + serialization engines every save and every clone
    goes through."""

    def test_digest_with_enzymes_50kb(self, baseline, random_50kb):
        feats = [{"type": "CDS", "start": i * 400, "end": i * 400 + 300,
                  "strand": 1, "label": f"g{i}"} for i in range(100)]
        enzymes = ["EcoRI", "BamHI", "HindIII", "PstI", "SalI", "XbaI"]

        def _work():
            sc._state._ENZYME_CUTS_CACHE.clear()
            sc._digest_with_enzymes(random_50kb, enzymes, circular=True,
                                    features=list(feats))

        best = _measure("digest_with_enzymes_50kb", baseline, _work)
        _assert_within_budget("digest_with_enzymes_50kb", baseline, best)

    def test_record_to_gb_text_100kb(self, baseline):
        rec = _synth_record_module(100_000, 500)

        def _work():
            sc._record_to_gb_text(rec)

        best = _measure("record_to_gb_text_100kb", baseline, _work)
        _assert_within_budget("record_to_gb_text_100kb", baseline, best)

    def test_gb_text_to_record_100kb(self, baseline):
        import splicecraft_record as _rec
        gb = sc._record_to_gb_text(_synth_record_module(100_000, 500))

        def _work():
            _rec._GB_PARSE_CACHE.clear()
            sc._gb_text_to_record(gb)

        best = _measure("gb_text_to_record_100kb", baseline, _work)
        _assert_within_budget("gb_text_to_record_100kb", baseline, best)

    def test_insdc_parse_genbank_100kb(self, baseline):
        import splicecraft_insdc as _insdc
        gb = sc._record_to_gb_text(_synth_record_module(100_000, 500))

        def _work():
            _insdc.parse_genbank(gb)

        best = _measure("insdc_parse_genbank_100kb", baseline, _work)
        _assert_within_budget("insdc_parse_genbank_100kb", baseline, best)


@pytest.mark.slow
class TestDesignEnginePerf:
    """The codon optimizer is mission-critical and the splice scan runs
    unasked on every plant-host optimize, so both are hot paths."""

    def test_codon_optimize_500aa(self, baseline, protein_500,
                                   flat_codon_table):
        def _work():
            sc._codon_optimize(protein_500, flat_codon_table)

        best = _measure("codon_optimize_500aa", baseline, _work)
        _assert_within_budget("codon_optimize_500aa", baseline, best)

    def test_splice_scan_5kb(self, baseline, random_5kb):
        import splicecraft_splice as _splice

        def _work():
            _splice._splice_scan(random_5kb)

        best = _measure("splice_scan_5kb", baseline, _work)
        _assert_within_budget("splice_scan_5kb", baseline, best)

    def test_codon_fix_splice_sites_300aa(self, baseline, flat_codon_table):
        """The fixture genuinely needs TWO synonymous repairs, so this gates
        the repair loop rather than just the scan — asserted below, because a
        fixture that drifted to scanning clean would quietly turn this into a
        much weaker test."""
        import splicecraft_codon as _codon
        rng = random.Random(12)
        prot = "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY")
                             for _ in range(299))
        dna = _codon._codon_optimize(prot, flat_codon_table)
        _out, _fixes, _rem = _codon._codon_fix_splice_sites(
            dna, prot, flat_codon_table)
        assert len(_fixes) >= 1, (
            "fixture no longer needs a repair, so this would time the clean "
            "fast path instead of the repair loop")

        def _work():
            _codon._codon_fix_splice_sites(dna, prot, flat_codon_table)

        best = _measure("codon_fix_splice_sites_300aa", baseline, _work)
        _assert_within_budget("codon_fix_splice_sites_300aa", baseline, best)


@pytest.mark.slow
class TestScannerPerf:
    """The regulatory / repeat scanners. `scan-terminators` folds a hairpin
    per candidate, so it is the heaviest per base and is heavy-classified on
    the agent side for that reason."""

    def test_scan_promoters_10kb(self, baseline, random_10kb):
        import splicecraft_regulatory as _reg

        def _work():
            _reg._scan_promoters(random_10kb)

        best = _measure("scan_promoters_10kb", baseline, _work)
        _assert_within_budget("scan_promoters_10kb", baseline, best)

    def test_scan_terminators_5kb(self, baseline, random_5kb):
        import splicecraft_regulatory as _reg

        def _work():
            _reg._scan_terminators(random_5kb)

        best = _measure("scan_terminators_5kb", baseline, _work)
        _assert_within_budget("scan_terminators_5kb", baseline, best)

    def test_cassette_internal_repeats_10kb(self, baseline, random_10kb):
        import splicecraft_cassette as _cas

        def _work():
            _cas._cassette_internal_repeats(random_10kb, 100)

        best = _measure("cassette_internal_repeats_10kb", baseline, _work)
        _assert_within_budget("cassette_internal_repeats_10kb", baseline, best)

    def test_find_orfs_50kb(self, baseline, random_50kb):
        def _work():
            sc._find_orfs(random_50kb, circular=True)

        best = _measure("find_orfs_50kb", baseline, _work)
        _assert_within_budget("find_orfs_50kb", baseline, best)


@pytest.mark.slow
class TestDataLayerPerf:
    """`_typed_clone` is what sacred #17 costs on every `_load_X()`, and
    `PlasmidMap.load_record` is the repaint tail undo / redo / flip share —
    both scale with COUNT, not bp ([INV-223])."""

    def test_typed_clone_library_1000(self, baseline):
        lib = [{"id": f"id{i}", "name": f"plasmid_{i}", "size": 5000,
                "n_feats": 12, "gb_text": "", "source": "",
                "features": [{"type": "CDS", "start": 100 + j * 300,
                              "end": 400 + j * 300, "strand": 1,
                              "label": f"g{j}",
                              "qualifiers": {"product": [f"p{j}"]}}
                             for j in range(12)]}
               for i in range(1000)]

        def _work():
            sc._typed_clone(lib)

        best = _measure("typed_clone_library_1000", baseline, _work)
        _assert_within_budget("typed_clone_library_1000", baseline, best)

    def test_plasmid_map_load_record_2000_feats(self, baseline):
        rec = _synth_record_module(20_000, 2000)
        pm = sc.PlasmidMap()

        def _work():
            pm.load_record(rec)

        best = _measure("plasmid_map_load_record_2000_feats", baseline, _work)
        _assert_within_budget("plasmid_map_load_record_2000_feats",
                              baseline, best)
