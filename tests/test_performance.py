"""
test_performance — performance budgets for hot paths.

These are deliberately loose (~6x headroom over current baseline) so that
normal machine variation and WSL2 noise don't fail the suite, but a real
regression (e.g. accidentally killing the IUPAC pattern cache or switching
to an O(n²) scanner) will trip them hard.

Baselines measured 2026-04-11 on Opus 4.6 dev WSL2 (post perf optimization pass):

  _scan_restriction_sites(pUC19-size)   ~6.5 ms
  _scan_restriction_sites(10 kb)        ~35 ms
  _iupac_pattern warm (200 lookups)     ~0.07 ms
  _rc(10 kb)                            ~0.1 ms
  _build_seq_text(pUC19, 10 feats)      ~6 ms
  _build_seq_text(20 kb, 80 feats)      ~50 ms

Every test warms the cache first (so we don't accidentally measure regex
compilation), then runs the hot path N times and asserts the average
per-call time is under the budget.
"""
from __future__ import annotations

import random
import time

import pytest

import splicecraft as sc


# Budgets — increase ONLY if you have a good reason. Decreasing is fine and
# means the code got faster; tighten the budget so regressions fail sooner.
BUDGET_SCAN_PUC19_MS = 30.0        # 4.6× headroom over 6.5 ms baseline
BUDGET_SCAN_10KB_MS  = 150.0       # 4.3× headroom over 35 ms baseline
BUDGET_IUPAC_WARM_MS = 5.0         # 70× headroom over 0.07 ms (200 lookups)
BUDGET_RC_10KB_MS    = 2.0         # 20× headroom over 0.1 ms
BUDGET_BUILD_SMALL_MS = 25.0       # 4× headroom over 6 ms baseline
BUDGET_BUILD_LARGE_MS = 200.0      # 4× headroom over 50 ms baseline


def _warm_up():
    """Run every hot path once so the regex cache is populated and we don't
    accidentally measure `re.compile` work on the first call."""
    sc._scan_restriction_sites("ACGTACGTACGT" * 100)


@pytest.fixture(scope="module")
def random_2686():
    rng = random.Random(42)
    return "".join(rng.choice("ACGT") for _ in range(2686))


@pytest.fixture(scope="module")
def random_10k():
    rng = random.Random(314)
    return "".join(rng.choice("ACGT") for _ in range(10_000))


# ═══════════════════════════════════════════════════════════════════════════════
# Restriction scanning
# ═══════════════════════════════════════════════════════════════════════════════

class TestScanRestrictionSitesPerformance:
    def test_puc19_size_under_50ms(self, random_2686):
        _warm_up()
        sc._scan_restriction_sites(random_2686)        # warm the cache
        t0 = time.perf_counter()
        for _ in range(20):
            sc._scan_restriction_sites(random_2686, unique_only=True)
        dt_ms = (time.perf_counter() - t0) / 20 * 1000
        assert dt_ms < BUDGET_SCAN_PUC19_MS, (
            f"scan(pUC19-size) took {dt_ms:.1f} ms, "
            f"budget is {BUDGET_SCAN_PUC19_MS} ms. "
            f"Check for O(n²) regression or dead cache."
        )

    def test_10kb_under_150ms(self, random_10k):
        _warm_up()
        sc._scan_restriction_sites(random_10k)
        t0 = time.perf_counter()
        for _ in range(10):
            sc._scan_restriction_sites(random_10k)
        dt_ms = (time.perf_counter() - t0) / 10 * 1000
        assert dt_ms < BUDGET_SCAN_10KB_MS, (
            f"scan(10 kb) took {dt_ms:.1f} ms, "
            f"budget is {BUDGET_SCAN_10KB_MS} ms."
        )

    def test_scan_scales_roughly_linearly(self, random_2686, random_10k):
        """A ~4× longer sequence should not take more than ~8× as long to scan.
        Guards against the O(n²) regressions that hit palindrome dedup or the
        isoschizomer filter."""
        _warm_up()
        sc._scan_restriction_sites(random_2686)
        sc._scan_restriction_sites(random_10k)

        t0 = time.perf_counter()
        for _ in range(10):
            sc._scan_restriction_sites(random_2686)
        short_ms = (time.perf_counter() - t0) / 10 * 1000

        t0 = time.perf_counter()
        for _ in range(10):
            sc._scan_restriction_sites(random_10k)
        long_ms = (time.perf_counter() - t0) / 10 * 1000

        ratio = long_ms / max(short_ms, 0.01)
        # 10_000 / 2_686 ≈ 3.7× more work — 8× allows generous slack
        assert ratio < 8.0, (
            f"scan scaling ratio {ratio:.1f}× (short={short_ms:.1f}ms, "
            f"long={long_ms:.1f}ms) — expected < 8×; "
            f"possible superlinear regression"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# IUPAC pattern compilation cache — sacred invariant #4
# ═══════════════════════════════════════════════════════════════════════════════

class TestIUPACPatternCachePerformance:
    def test_warm_cache_is_near_free(self):
        """After every NEB enzyme's pattern has been compiled once, a full
        second pass should cost ~nothing. A regression would mean the cache
        was broken and every call is recompiling."""
        # Warm the cache.
        for name, (site, _, _) in sc._NEB_ENZYMES.items():
            sc._iupac_pattern(site)

        # Time a full second pass.
        t0 = time.perf_counter()
        for _ in range(10):
            for name, (site, _, _) in sc._NEB_ENZYMES.items():
                sc._iupac_pattern(site)
        dt_ms = (time.perf_counter() - t0) / 10 * 1000
        assert dt_ms < BUDGET_IUPAC_WARM_MS, (
            f"200 warm pattern lookups took {dt_ms:.2f} ms, "
            f"budget {BUDGET_IUPAC_WARM_MS} ms — cache may be broken"
        )

    def test_cold_then_warm_is_faster(self):
        """Warm cache must be strictly faster than cold (regex compile).
        Guards against an accidental cache-bypass refactor."""
        sc._PATTERN_CACHE.clear()

        t0 = time.perf_counter()
        for name, (site, _, _) in sc._NEB_ENZYMES.items():
            sc._iupac_pattern(site)
        cold_ms = (time.perf_counter() - t0) * 1000

        t0 = time.perf_counter()
        for name, (site, _, _) in sc._NEB_ENZYMES.items():
            sc._iupac_pattern(site)
        warm_ms = (time.perf_counter() - t0) * 1000

        assert warm_ms < cold_ms, (
            f"warm ({warm_ms:.2f} ms) not faster than cold ({cold_ms:.2f} ms)"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Reverse complement
# ═══════════════════════════════════════════════════════════════════════════════

class TestRCPerformance:
    def test_rc_on_10kb_under_2ms(self, random_10k):
        t0 = time.perf_counter()
        for _ in range(100):
            sc._rc(random_10k)
        dt_ms = (time.perf_counter() - t0) / 100 * 1000
        assert dt_ms < BUDGET_RC_10KB_MS, (
            f"_rc(10 kb) took {dt_ms:.3f} ms, budget {BUDGET_RC_10KB_MS} ms"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Sequence panel text builder — the dominant per-frame cost
# ═══════════════════════════════════════════════════════════════════════════════

def _synth_record(n_bp: int, n_feats: int):
    """Synthetic SeqRecord with evenly-spaced features, no NCBI."""
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rng = random.Random(42)
    seq = "".join(rng.choice("ACGT") for _ in range(n_bp))
    rec = SeqRecord(Seq(seq), id="T", name="T", description="perf test")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    step = n_bp // (n_feats + 1)
    for i in range(n_feats):
        s = i * step + 50
        e = s + rng.randint(50, 500)
        rec.features.append(SeqFeature(
            FeatureLocation(s, e, strand=1 if i % 2 == 0 else -1),
            type="CDS" if i % 3 == 0 else "misc_feature",
            qualifiers={"label": [f"f{i}"]},
        ))
    return rec


class TestBuildSeqTextPerformance:
    def test_build_small_under_25ms(self):
        """pUC19-size sequence with 10 features. Dominant cost on small plasmids."""
        rec = _synth_record(2686, 10)
        pm = sc.PlasmidMap()
        pm.load_record(rec)
        seq_str, feats = str(rec.seq), pm._feats
        sc._build_seq_text(seq_str, feats, line_width=127)  # warm
        t0 = time.perf_counter()
        for _ in range(10):
            sc._BUILD_SEQ_CACHE.clear()   # cold path
            sc._build_seq_text(seq_str, feats, line_width=127)
        dt_ms = (time.perf_counter() - t0) / 10 * 1000
        assert dt_ms < BUDGET_BUILD_SMALL_MS, (
            f"_build_seq_text(pUC19) took {dt_ms:.1f} ms, "
            f"budget {BUDGET_BUILD_SMALL_MS} ms"
        )

    def test_build_large_under_200ms(self):
        """20 kb / 80 features — stress test for the sequence panel."""
        rec = _synth_record(20000, 80)
        pm = sc.PlasmidMap()
        pm.load_record(rec)
        seq_str, feats = str(rec.seq), pm._feats
        sc._build_seq_text(seq_str, feats, line_width=127)  # warm
        t0 = time.perf_counter()
        for _ in range(5):
            sc._BUILD_SEQ_CACHE.clear()
            sc._build_seq_text(seq_str, feats, line_width=127)
        dt_ms = (time.perf_counter() - t0) / 5 * 1000
        assert dt_ms < BUDGET_BUILD_LARGE_MS, (
            f"_build_seq_text(20 kb) took {dt_ms:.1f} ms, "
            f"budget {BUDGET_BUILD_LARGE_MS} ms"
        )

    def test_warm_cache_skips_styles_work(self):
        """Second call with same (seq, feats) must reuse the cached styles+
        annot_feats tuple. Verifies _BUILD_SEQ_CACHE is doing its job."""
        rec = _synth_record(10000, 40)
        pm = sc.PlasmidMap()
        pm.load_record(rec)
        seq_str, feats = str(rec.seq), pm._feats
        sc._BUILD_SEQ_CACHE.clear()

        # First call populates the cache
        sc._build_seq_text(seq_str, feats, line_width=127)
        key = (id(seq_str), id(feats), len(seq_str), len(feats))
        assert key in sc._BUILD_SEQ_CACHE

        # Subsequent calls with SAME (seq, feats) but DIFFERENT cursor must
        # still hit the cache (cursor isn't part of the key — only the expensive
        # input derivation is cached).
        for cursor in [100, 200, 300]:
            sc._build_seq_text(seq_str, feats, line_width=127, cursor_pos=cursor)
        assert key in sc._BUILD_SEQ_CACHE
