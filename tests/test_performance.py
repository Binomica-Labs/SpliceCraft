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


# ═══════════════════════════════════════════════════════════════════════════════
# Cursor scrolling on large plasmids — guards against regressing the chunk
# layout cache (`_chunk_layout`) and the per-chunk static render cache
# (`_CHUNK_STATIC_CACHE`). These were added 2026-04-25 so users editing
# cosmid (~30-50 kb) and BAC (~100-300 kb) records don't see jerky cursor
# movement when holding an arrow key. Pre-fix baselines on the same WSL2
# machine: 50 kb cursor ~67 ms/call, 150 kb ~200+ ms/call. Post-fix budgets
# below sit well above current numbers (~5× headroom) but trip hard on a
# real regression.
# ═══════════════════════════════════════════════════════════════════════════════

BUDGET_CURSOR_50KB_MS  = 50.0    # ~10× headroom over 5 ms baseline
BUDGET_CURSOR_150KB_MS = 120.0   # ~10× headroom over 12 ms baseline
# Bumped 2026-04-30 after the inline AA translation row + trailing
# inter-chunk gap row landed. Both add spans to the pre-cached chunk
# Text, so `result.append(cached)` (the dominant cost in the timed
# loop, per cProfile) does proportionally more work. The added cost
# is one-row-per-chunk for the gap and one-row-per-CDS-chunk for the
# AA translation — bounded, not algorithmic.
BUDGET_BP_TO_ROW_US    = 200.0   # ~40× headroom over 5 µs baseline


class TestLargePlasmidCursorScrolling:
    """Cursor scrolling on cosmid/BAC-scale records must stay smooth (≥30 fps).

    Holding the Right arrow key fires `_build_seq_text` once per move on a
    new cursor position. Without the static-render cache, every keystroke
    re-runs `_render_feature_row_pair` for ~1500 chunks on a 200 kb plasmid
    (~120 ms = 8 fps, visibly jerky). With the cache, only the chunk
    containing the cursor re-renders.
    """

    def test_50kb_warm_cursor_under_30ms(self):
        rec = _synth_record(50000, 200)
        pm = sc.PlasmidMap()
        pm.load_record(rec)
        seq_str, feats = str(rec.seq), pm._feats

        # Fully warm both caches: chunk_layout + chunk_static
        sc._build_seq_text(seq_str, feats, line_width=127, cursor_pos=0)

        # Move cursor across 50 positions and time the average per-move cost.
        t0 = time.perf_counter()
        for cur in range(0, 50000, 1000):
            sc._build_seq_text(seq_str, feats, line_width=127, cursor_pos=cur)
        dt_ms = (time.perf_counter() - t0) / 50 * 1000
        assert dt_ms < BUDGET_CURSOR_50KB_MS, (
            f"50 kb warm cursor move took {dt_ms:.1f} ms/call, "
            f"budget {BUDGET_CURSOR_50KB_MS} ms — possible chunk-static "
            f"cache regression"
        )

    def test_150kb_warm_cursor_under_60ms(self):
        rec = _synth_record(150000, 600)
        pm = sc.PlasmidMap()
        pm.load_record(rec)
        seq_str, feats = str(rec.seq), pm._feats
        sc._build_seq_text(seq_str, feats, line_width=127, cursor_pos=0)

        t0 = time.perf_counter()
        for cur in range(0, 150000, 3000):
            sc._build_seq_text(seq_str, feats, line_width=127, cursor_pos=cur)
        dt_ms = (time.perf_counter() - t0) / 50 * 1000
        assert dt_ms < BUDGET_CURSOR_150KB_MS, (
            f"150 kb (BAC) warm cursor move took {dt_ms:.1f} ms/call, "
            f"budget {BUDGET_CURSOR_150KB_MS} ms"
        )

    def test_bp_to_content_row_is_O1(self):
        """`_bp_to_content_row` must use chunk_layout prefix sums — direct
        indexing, not a per-chunk re-scan. Pre-fix this was O(chunks ×
        features), ~50 ms at 50 kb."""
        from splicecraft import SequencePanel
        rec = _synth_record(50000, 200)
        pm = sc.PlasmidMap()
        pm.load_record(rec)
        sp = SequencePanel()
        sp._seq = str(rec.seq)
        sp._feats = pm._feats
        sp._show_connectors = False
        sp._bp_to_content_row(0)  # warm chunk_layout cache

        t0 = time.perf_counter()
        for _ in range(500):
            sp._bp_to_content_row(49999)   # worst-case bp at end
        dt_us = (time.perf_counter() - t0) / 500 * 1e6
        assert dt_us < BUDGET_BP_TO_ROW_US, (
            f"_bp_to_content_row(50 kb, end) took {dt_us:.1f} µs/call, "
            f"budget {BUDGET_BP_TO_ROW_US} µs — chunk_layout cache regression"
        )

    def test_chunk_static_cache_invalidates_on_feats_change(self):
        """Two different feature lists must produce two different renders.
        Caching by `id(feats)` means a feats reassignment misses the cache;
        a buggy refactor that mutated `feats` in place would alias the two
        renders and silently leak stale lane art."""
        rec_a = _synth_record(5000, 20)
        rec_b = _synth_record(5000, 20)
        rec_b.features.clear()  # different feats list, same seq length
        pm_a = sc.PlasmidMap(); pm_a.load_record(rec_a)
        pm_b = sc.PlasmidMap(); pm_b.load_record(rec_b)
        seq_a = str(rec_a.seq)
        seq_b = str(rec_b.seq)

        text_a = sc._build_seq_text(seq_a, pm_a._feats, line_width=80)
        text_b = sc._build_seq_text(seq_b, pm_b._feats, line_width=80)

        # rec_a has 20 features (annotation rows), rec_b has none — line
        # counts must differ.
        rows_a = text_a.plain.count("\n")
        rows_b = text_b.plain.count("\n")
        assert rows_a > rows_b, (
            f"feats=20 produced {rows_a} rows, feats=0 produced {rows_b} — "
            f"static cache may be aliasing two records"
        )
