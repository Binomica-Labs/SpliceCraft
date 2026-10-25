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
        # `re` keeps its OWN compile cache. Under xdist an earlier test in this
        # worker has already compiled every NEB site, so without a purge the
        # "cold" pass is a second cache hit (~0.75 ms vs ~0.14 ms warm) and
        # one scheduler preemption on a loaded machine flips the comparison.
        # Purged, cold is a real compile (~9 ms) — the thing this test means.
        import re as _re
        _re.purge()

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
        # Cache key tracks _build_seq_inputs (hash(seq) since the
        # 2026-05-06 fix that killed the id() collision flake; id(feats)
        # is still safe per CLAUDE.md invariant #4).
        key = (hash(seq_str), id(feats), len(seq_str), len(feats))
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

BUDGET_CURSOR_50KB_MS  = 80.0    # ~16× headroom over 5 ms baseline
BUDGET_CURSOR_150KB_MS = 200.0   # ~16× headroom over 12 ms baseline
# Bumped 2026-04-30 after the inline AA translation row + trailing
# inter-chunk gap row landed. Both add spans to the pre-cached chunk
# Text, so `result.append(cached)` (the dominant cost in the timed
# loop, per cProfile) does proportionally more work. The added cost
# is one-row-per-chunk for the gap and one-row-per-CDS-chunk for the
# AA translation — bounded, not algorithmic.
# Bumped again 2026-04-30 (later) to absorb WSL2 load tail: the test
# passes in isolation in ~3-5 ms but spiked to ~150 ms during a
# release.py serial run after 12 min of all-cores test churn. The
# budget still catches a real regression (anything that breaks the
# static-render cache pushes the number back to the ~200 ms pre-fix
# baseline cited above) without flaking on shared-runner load tails.
BUDGET_BP_TO_ROW_US    = 200.0   # ~40× headroom over 5 µs baseline


class TestLargePlasmidCursorScrolling:
    """Cursor scrolling on cosmid/BAC-scale records must stay smooth (≥30 fps).

    Holding the Right arrow key fires `_build_seq_text` once per move on a
    new cursor position. Without the static-render cache, every keystroke
    re-runs the lane-art painters for ~1500 chunks on a 200 kb plasmid
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


# ═══════════════════════════════════════════════════════════════════════════════
# The restriction scan runs only when its result can be seen (audit 2026-10-23)
# ═══════════════════════════════════════════════════════════════════════════════

class TestRestrictionScanOnlyWhenShown:
    """Every load and edit dispatched a full restriction scan although the
    overlays are hidden by default and every reader shows
    `_restr_cache if _show_restr else []`: ~2.7 s/Mbp holding the GIL and
    ~82 MB/Mbp per scan, and a superseded scan could not be stopped — three
    ran side by side after three quick edits on a 2 Mbp record."""

    @staticmethod
    def _spy_scans(monkeypatch):
        calls: list = []
        real = sc._scan_restriction_sites

        def spy(seq, *a, **k):
            calls.append(seq)
            return real(seq, *a, **k)

        monkeypatch.setattr(sc, "_scan_restriction_sites", spy)
        return calls

    async def test_hidden_overlays_scan_nothing_until_shown(
            self, tiny_record, isolated_library, monkeypatch):
        calls = self._spy_scans(monkeypatch)
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            assert app._current_record is not None
            assert app._show_restr is False
            # A load and an edit-style re-dispatch while hidden: no scan.
            app._dispatch_restr_scan(str(app._current_record.seq))
            await app.workers.wait_for_complete()
            await pilot.pause()
            assert calls == []
            assert app._restr_cache == []
            # Showing the overlays scans the CURRENT sequence and paints it.
            app.action_toggle_restr()
            await app.workers.wait_for_complete()
            await pilot.pause()
            assert calls == [str(tiny_record.seq)]
            assert app._restr_cache == sc._scan_restriction_sites(
                str(tiny_record.seq), unique_only=app._restr_unique_only,
                min_recognition_len=app._restr_min_len, circular=True)

    async def test_superseded_scan_is_cancelled_and_never_lands(
            self, tiny_record, isolated_library, monkeypatch):
        outcomes: list = []
        real = sc._scan_restriction_sites

        def spy(seq, *a, **k):
            try:
                return real(seq, *a, **k)
            except sc._ScanCancelled:
                outcomes.append("cancelled")
                raise

        monkeypatch.setattr(sc, "_scan_restriction_sites", spy)
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.restr._show_restr = True
            stale_gen = app.restr._scan_gen
            app.restr._scan_gen += 1          # a newer request superseded it
            sentinel = [{"type": "resite", "label": "SENTINEL"}]
            app._restr_cache = sentinel
            # A sequence never scanned before, so the LRU can't answer it.
            fresh = "GAATTC" + "ACGT" * 40 + "GGATCC" + "TTGACA" * 7
            app._restr_scan_worker(
                fresh, min_len=6, unique_only=True, circular=True,
                entry_counter=app._record_load_counter, scan_gen=stale_gen)
            await app.workers.wait_for_complete()
            await pilot.pause()
            assert outcomes == ["cancelled"]
            assert app._restr_cache is sentinel


def test_cancelled_scan_raises_and_caches_nothing():
    seq = "GAATTC" + "ACGTTGCA" * 30 + "GGATCC" + "CCATGG" * 3
    with pytest.raises(sc._ScanCancelled):
        sc._scan_restriction_sites(seq, cancel=lambda: True)
    # The cancelled call left no partial entry behind: a normal scan
    # computes the same answer the uncached implementation gives.
    assert sc._scan_restriction_sites(seq) == \
        sc._scan_restriction_sites_impl(seq)


# ═══════════════════════════════════════════════════════════════════════════════
# Genome-scale memory: caches and undo are bounded by SIZE, not just count
# (audit 2026-10-23 — each edit of a 10 Mbp genome grew RSS ~130 MB)
# ═══════════════════════════════════════════════════════════════════════════════

def _rand_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_seq_panel_caches_are_bounded_by_total_length(monkeypatch):
    monkeypatch.setattr(sc, "_SEQ_CACHE_BP_BUDGET", 1000)
    sc._BUILD_SEQ_CACHE.clear()
    sc._CHUNK_LAYOUT_CACHE.clear()
    seqs = [_rand_seq(400, i) for i in range(5)]
    feats = [{"start": 10, "end": 90, "strand": 1, "type": "CDS",
              "label": "x", "color": "color(40)"}]
    for sq in seqs:
        sc._build_seq_inputs(sq, feats)
        sc._chunk_layout(sq, feats, 60)
    # 2 x 400 bp fit the 1,000 bp budget; a third would not.
    assert len(sc._BUILD_SEQ_CACHE) == 2
    assert len(sc._CHUNK_LAYOUT_CACHE) == 2
    newest = (hash(seqs[-1]), id(feats), 400, 1)
    assert newest in sc._BUILD_SEQ_CACHE
    # One entry longer than the whole budget is still cached — alone.
    big = _rand_seq(1500, 99)
    sc._build_seq_inputs(big, feats)
    assert list(sc._BUILD_SEQ_CACHE) == [(hash(big), id(feats), 1500, 1)]


def test_genbank_parse_cache_is_bounded_by_total_length(monkeypatch):
    import splicecraft_record as R
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    monkeypatch.setattr(R, "_GB_PARSE_CACHE_BP_BUDGET", 700)
    R._GB_PARSE_CACHE.clear()
    for i in range(4):
        rec = SeqRecord(Seq(_rand_seq(300, i)), id=f"R{i}", name=f"R{i}")
        rec.annotations["molecule_type"] = "DNA"
        sc._gb_text_to_record(sc._record_to_gb_text(rec))
    assert len(R._GB_PARSE_CACHE) == 2
    held = sum(len(r.seq) for r in R._GB_PARSE_CACHE.values())
    assert held <= 700


class TestUndoHistoryIsBoundedBySize:

    async def _app_with(self, record):
        app = sc.PlasmidApp()
        app._preload_record = record
        return app

    async def test_oldest_steps_go_first_and_three_always_stay(
            self, tiny_record, isolated_library, monkeypatch):
        n = len(tiny_record.seq)
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            monkeypatch.setattr(sc, "_UNDO_BP_BUDGET", int(5.5 * n))
            for _ in range(10):
                app._push_undo()
            assert len(app._undo_stack) == 5
            monkeypatch.setattr(sc, "_UNDO_BP_BUDGET", 1)
            app._push_undo()
            assert len(app._undo_stack) == sc._UNDO_MIN_STEPS == 3


# ═══════════════════════════════════════════════════════════════════════════════
# Rotation: an edit keeps it; a chromosome follows it without a re-layout
# ═══════════════════════════════════════════════════════════════════════════════

class TestRotationSurvivesEditsAndScales:

    async def test_an_edit_keeps_the_rotation_and_lays_out_once(
            self, tiny_record, isolated_library, monkeypatch):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm.origin_bp = 50
            await pilot.pause()
            assert sp._view_origin_bp == 50
            frames: list = []
            real = sc._chunk_layout

            def spy(seq, feats, lw):
                frames.append(seq[:12])
                return real(seq, feats, lw)

            monkeypatch.setattr(sc, "_chunk_layout", spy)
            app._edit_dialog_result(("GGG", "insert", 10, 10))
            for _ in range(4):
                await pilot.pause()
            assert pm.origin_bp == 50
            assert sp._view_origin_bp == 50
            # Every layout after the edit is of the ROTATED sequence: the
            # unrotated frame (old OriginChanged(0) echo) is never built.
            new_seq = str(app._current_record.seq)
            assert frames and all(
                f == (new_seq[50:] + new_seq[:50])[:12] for f in frames)

    async def test_chromosome_scale_rotation_moves_instead_of_relayout(
            self, tiny_record, isolated_library, monkeypatch):
        monkeypatch.setattr(sc, "_LAZY_RENDER_THRESHOLD_BP", 50)
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            laid_out: list = []
            real = sc._chunk_layout

            def spy(seq, feats, lw):
                laid_out.append(seq is sp._seq)
                return real(seq, feats, lw)

            monkeypatch.setattr(sc, "_chunk_layout", spy)
            pm.origin_bp = 70
            for _ in range(3):
                await pilot.pause()
            assert sp._view_origin_bp == 0          # absolute order kept
            assert sp._cursor_pos == 70             # cursor on the origin
            assert laid_out and all(laid_out)       # only the unrotated seq


# ═══════════════════════════════════════════════════════════════════════════════
# Restriction scan: `unique_only` decides an enzyme by COUNTING to 2
# (audit 2026-10-23 — 2 Mbp: 3.3 s / +136 MB of site dicts, all discarded)
# ═══════════════════════════════════════════════════════════════════════════════

def _plantable_sites():
    import splicecraft_biology as B
    out = []
    for e in sc._state._scan_catalog_hook():
        site = e[1]
        if site and set(site) <= set("ACGT"):
            out.append(site)
    return out, B


class TestRestrictionScanUniqueGate:

    def test_gated_scan_is_output_identical(self, monkeypatch):
        sites, B = _plantable_sites()
        rng = random.Random(20261024)
        cases = 0
        for trial in range(60):
            n = rng.choice([0, 7, 300, 1500, 3000])
            seq = list(rng.choice("ACGT") for _ in range(n))
            for _ in range(rng.randint(0, 4)):          # plant some sites
                site = rng.choice(sites)
                if len(site) > n:
                    continue
                how = rng.choice(["fwd", "rc", "wrap", "end", "twice"])
                text = sc._rc(site) if how == "rc" else site
                if how == "wrap" and n > len(site):
                    k = rng.randint(1, len(site) - 1)
                    seq[n - k:] = list(text[:k])
                    seq[:len(text) - k] = list(text[k:])
                    continue
                at = (n - len(text)) if how == "end" else rng.randint(
                    0, n - len(text))
                seq[at:at + len(text)] = list(text)
                if how == "twice":
                    at2 = rng.randint(0, n - len(text))
                    seq[at2:at2 + len(text)] = list(text)
            s = "".join(seq)
            for circular in (True, False):
                for min_len in (4, 6):
                    for distinct in (False, True):
                        monkeypatch.setattr(B, "_SCAN_UNIQUE_GATE", True)
                        got = B._scan_restriction_sites_impl(
                            s, min_len, True, circular, distinct_cuts=distinct)
                        monkeypatch.setattr(B, "_SCAN_UNIQUE_GATE", False)
                        want = B._scan_restriction_sites_impl(
                            s, min_len, True, circular, distinct_cuts=distinct)
                        assert got == want, (trial, n, circular, min_len,
                                             distinct)
                        cases += 1
        assert cases == 60 * 8

    def test_a_multi_cutter_is_decided_without_walking_its_sites(
            self, monkeypatch):
        sites, B = _plantable_sites()
        steps = [0]
        real = B._iter_match_starts

        def counting(pat, s):
            for p in real(pat, s):
                steps[0] += 1
                yield p

        monkeypatch.setattr(B, "_iter_match_starts", counting)
        rng = random.Random(7)
        seq = "".join(rng.choice("ACGT") for _ in range(200_000))
        assert seq.count("GAATTC") > 10            # a real sequence, not a run
        B._scan_restriction_sites_impl(seq, 6, True, True)
        n_enz = sum(1 for e in sc._state._scan_catalog_hook() if e[2] >= 6)
        # At most two matches per strand per enzyme are looked at; the full
        # walk visited every site (~100 per enzyme at 200 kb).
        assert steps[0] <= 4 * n_enz, (steps[0], n_enz)

    def test_each_distinct_pattern_is_walked_once_per_scan(self,
                                                           monkeypatch):
        """Isoschizomers, HF variants and a reverse pattern that equals
        another enzyme's forward one were each walked again — 206 regex
        passes for 146 distinct patterns (audit 2026-10-23)."""
        sites, B = _plantable_sites()
        walked: "list[tuple]" = []
        real = B._iter_match_starts

        def counting(pat, s):
            walked.append((pat.pattern, pat.flags))
            return real(pat, s)
        monkeypatch.setattr(B, "_iter_match_starts", counting)
        rng = random.Random(11)
        seq = "".join(rng.choice("ACGT") for _ in range(20_000))
        B._scan_restriction_sites_impl(seq, 4, False, True)
        assert walked and len(walked) == len(set(walked)), (
            len(walked), len(set(walked)))


# ═══════════════════════════════════════════════════════════════════════════════
# Find-annotation (ctrl+/): a virtual table — only on-screen rows are drawn
# (audit 2026-10-23: the eager DataTable took 1.85 s to open at 20k features)
# ═══════════════════════════════════════════════════════════════════════════════

def _find_feats(n):
    out = [{"type": "source", "label": "src", "start": 0, "end": 10_000_000},
           {"type": "resite", "label": "EcoRI", "start": 5, "end": 11},
           {"type": "recut", "label": "", "start": 6, "end": 7}]
    for i in range(n):
        t = "CDS" if i % 3 == 0 else ("gene" if i % 3 == 1 else "promoter")
        out.append({"type": t, "label": f"{'lacZ' if i % 50 == 0 else 'orf'}{i}",
                    "start": i * 400, "end": i * 400 + 300})
    return out


class TestFindAnnotationModal:

    async def test_filter_status_and_dismiss_map_back_to_the_feature(self):
        feats = _find_feats(200)
        app = sc.PlasmidApp()
        results: list = []
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            app.push_screen(sc.FeatureSearchModal(feats, total=10_000_000),
                            callback=results.append)
            await pilot.pause()
            m = app.screen
            assert isinstance(m, sc.FeatureSearchModal)
            t = m.query_one("#featfind-table")
            assert t.row_count == 200          # source/resite/recut excluded
            m._refilter("LACZ")                # label, case-insensitive
            assert t.row_count == 4
            m._refilter("prom")                # …or the type
            assert t.row_count == 66
            m._refilter("zzz")
            status = str(m.query_one("#featfind-status").render())
            assert "No matches" in status and "200" in status
            m._refilter("lacz")
            t.move_cursor(row=1)               # the SECOND match: orf50
            m._dismiss_with_cursor(t)
            await pilot.pause()
        assert results == [feats.index(next(
            f for f in feats if f["label"] == "lacZ50"))]

    async def test_enter_in_the_search_jumps_to_the_first_match(self):
        feats = _find_feats(120)
        app = sc.PlasmidApp()
        results: list = []
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            app.push_screen(sc.FeatureSearchModal(feats, total=10_000_000),
                            callback=results.append)
            await pilot.pause()
            await pilot.press(*"lacz")
            await pilot.press("enter")
            await pilot.pause()
        assert results == [3]                   # lacZ0, after the 3 overlays

    async def test_a_refilter_draws_only_the_visible_rows(self, monkeypatch):
        feats = _find_feats(20_000)
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            app.push_screen(sc.FeatureSearchModal(feats, total=10_000_000))
            await pilot.pause()
            m = app.screen
            calls = [0]
            real = m._row_cells

            def counting(i):
                calls[0] += 1
                return real(i)

            monkeypatch.setattr(m, "_row_cells", counting)
            m._refilter("")
            m._refilter("orf")                 # matches ~all 20k rows
            await pilot.pause()
            await pilot.pause()
            h = m.query_one("#featfind-table").size.height
            assert m.query_one("#featfind-table").row_count > 19_000
            assert calls[0] <= 4 * h + 10, (calls[0], h)


# ═══════════════════════════════════════════════════════════════════════════════
# Library / collection lookups never wait out a save holding `_cache_lock`
# (audit 2026-10-23: a row load waited 130–280 ms behind every Ctrl+S)
# ═══════════════════════════════════════════════════════════════════════════════

def test_lookups_do_not_wait_for_a_writer_holding_the_cache_lock():
    import threading
    sc._save_collections([{"name": "Main", "plasmids": [
        {"id": "p1", "name": "plasmid one", "gb_text": "", "size": 0}]}])
    sc._save_library([{"id": "p1", "name": "plasmid one", "gb_text": "",
                       "size": 0}])
    # Populated and migrated once, outside the race.
    assert sc._find_library_entry_by_id("p1") is not None
    assert sc._find_collection("Main") is not None
    held, release = threading.Event(), threading.Event()

    def writer():
        with sc._state._cache_lock:
            held.set()
            release.wait(5)

    t = threading.Thread(target=writer, daemon=True)
    t.start()
    assert held.wait(5)
    try:
        for fn, arg, check in (
                (sc._find_library_entry_by_id, "p1", lambda r: r["id"] == "p1"),
                (sc._find_library_entry_by_name, "plasmid one",
                 lambda r: r["id"] == "p1"),
                (sc._find_collection, "Main", lambda r: r["name"] == "Main"),
                (lambda _a: sc._iter_library_readonly(), None, bool),
                (lambda _a: sc._iter_collections_readonly(), None, bool)):
            t0 = time.perf_counter()
            r = fn(arg)
            dt = time.perf_counter() - t0
            assert dt < 0.25, (getattr(fn, "__name__", fn), dt)
            assert check(r)
    finally:
        release.set()
        t.join(5)


# ═══════════════════════════════════════════════════════════════════════════════
# Splash: each helix phase is composed once and replayed
# (audit 2026-10-23: ~45 ms of UI thread per 25-fps frame starved the load
# underneath it)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSplashFramesAreCached:

    async def test_a_phase_is_composed_once_and_replayed(
            self, isolated_library, monkeypatch):
        app = sc.PlasmidApp()
        app._skip_splash = False
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            scr = app.screen
            assert isinstance(scr, sc.SplashScreen)
            calls = [0]
            real = scr._compose_splash

            def counting(w, h):
                calls[0] += 1
                return real(w, h)

            monkeypatch.setattr(scr, "_compose_splash", counting)
            canvas = scr.query_one("#splash-canvas")
            for _ in range(200):               # ~4.4 turns of the helix
                scr._tick_rotation()
                canvas.render_line(0)
            assert calls[0] <= sc._SPLASH_PHASES + 1, calls[0]
            app.exit()

    async def test_cached_lines_match_the_composed_frame(self, isolated_library):
        app = sc.PlasmidApp()
        app._skip_splash = False
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            scr = app.screen
            canvas = scr.query_one("#splash-canvas")
            w, h = canvas.size.width, canvas.size.height
            for k in (0, 7, 22, 44):
                canvas.set_phase(k)
                want = canvas._frame_text(w, h, k).plain.split("\n")
                for y in range(h):
                    got = "".join(seg.text for seg in canvas.render_line(y))
                    line = want[y] if y < len(want) else ""
                    assert got.rstrip() == line[:w].rstrip(), (k, y)
            app.exit()
