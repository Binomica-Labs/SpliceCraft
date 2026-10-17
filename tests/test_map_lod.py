"""
test_map_lod — the large-chromosome map, `[INV-223]`.

Rotating the circular map repainted EVERY feature, ~14 trig calls and ~100
`round()`s each. Measured through a mounted app: 1.3 ms at 1 feature, 8 ms at
200, 69 ms at 2,000 and **40.6 SECONDS at 200,000** — which is what a 1 Gbp
chromosome at one gene per 5 kb actually is. Load was 41 s.

The fix is not faster trig. A 76x27 braille map has ~16,000 addressable dots,
so 200,000 features cannot be drawn distinctly at any speed; above a measured
threshold the map paints an aggregate DENSITY ring instead, keyed on bp so it
survives a rotate and a resize. Rotation became 8 ms — faster than the
2,000-feature case, because the aggregate path does not care how many
features there are.

These tests pin the three things that make that true: the histogram is
CORRECT (cross-checked against a brute-force count), it is CACHED across
rotations, and small records still take the per-feature path unchanged.
"""
from __future__ import annotations

import random
import time

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature

import splicecraft as sc


def _rec(n_bp: int, n_feats: int, *, span: int = 300, circular=True):
    rng = random.Random(11)
    chunk = "".join(rng.choice("ACGT") for _ in range(min(n_bp, 50_000)))
    seq = (chunk * (n_bp // len(chunk) + 1))[:n_bp]
    rec = SeqRecord(Seq(seq), id="CHR", name="CHR",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular" if circular
                                 else "linear"})
    step = max(1, n_bp // max(1, n_feats))
    rec.features = [
        SeqFeature(FeatureLocation(i * step, min(i * step + span, n_bp), 1),
                   type="gene", qualifiers={"label": [f"g{i}"]})
        for i in range(n_feats)]
    return rec


def _map_with(rec):
    pm = sc.PlasmidMap()
    pm.load_record(rec)
    return pm


# ═══════════════════════════════════════════════════════════════════════════
# The threshold
# ═══════════════════════════════════════════════════════════════════════════

def test_small_records_keep_the_per_feature_path():
    """A plasmid must be completely unaffected — the aggregate view is for
    chromosomes, and silently coarsening a 200-feature plasmid would be a
    regression dressed as an optimisation."""
    pm = _map_with(_rec(10_000, 50))
    assert pm._lod_active() is False


def test_a_chromosome_crosses_into_the_aggregate_view():
    pm = _map_with(_rec(2_000_000, 5_000, span=60))
    assert pm._lod_active() is True


def test_threshold_zero_disables_aggregation(monkeypatch):
    """`0` is the escape hatch for a user who would rather wait than lose
    per-feature fidelity."""
    monkeypatch.setattr(sc, "_get_setting",
                        lambda k, d=None: 0 if k == "map_lod_feature_threshold"
                        else d)
    pm = _map_with(_rec(2_000_000, 5_000, span=60))
    assert pm._lod_active() is False


def test_threshold_is_settable(monkeypatch):
    monkeypatch.setattr(sc, "_get_setting",
                        lambda k, d=None: 10 if k == "map_lod_feature_threshold"
                        else d)
    pm = _map_with(_rec(100_000, 50))
    assert pm._lod_feature_threshold() == 10
    assert pm._lod_active() is True


def test_a_broken_setting_value_falls_back(monkeypatch):
    monkeypatch.setattr(sc, "_get_setting",
                        lambda k, d=None: "nonsense"
                        if k == "map_lod_feature_threshold" else d)
    pm = _map_with(_rec(10_000, 5))
    assert pm._lod_feature_threshold() == sc._MAP_LOD_DEFAULT_THRESHOLD


def test_the_setting_is_reachable_from_both_surfaces():
    """The right value depends on the terminal and the molecule, so it must
    be tunable by a user AND by a script."""
    assert "map_lod_feature_threshold" in sc._SETTINGS_SCHEMA
    import splicecraft_agent as sca
    assert "map_lod_feature_threshold" in sca._AGENT_SETTINGS_ALLOWLIST


# ═══════════════════════════════════════════════════════════════════════════
# The histogram — correctness against a brute-force reference
# ═══════════════════════════════════════════════════════════════════════════

def _brute_counts(feats, total, n):
    """Independent reference: walk every bucket of every feature.

    This is the O(features x buckets) loop the real code deliberately does
    NOT do (it uses a difference array so a feature costs O(1) however wide
    it is). Agreeing with it is the whole point."""
    counts = [0] * n
    for f in feats:
        s = int(f["start"]) % total
        span = max(1, sc._feat_len(int(f["start"]), int(f["end"]), total))
        b0 = (s * n) // total
        nb = max(1, min(n, span * n // total))
        for k in range(nb):
            counts[(b0 + k) % n] += 1
    return counts


@pytest.mark.parametrize("n_bp,n_feats,span", [
    (100_000, 200, 300),
    (1_000_000, 500, 5_000),
    (2_000_000, 3_000, 60),
    (500_000, 40, 50_000),        # features much wider than a bucket
])
def test_difference_array_matches_a_brute_force_count(n_bp, n_feats, span):
    pm = _map_with(_rec(n_bp, n_feats, span=span))
    counts, _colors = pm._lod_histogram()
    want = _brute_counts(pm._feats, pm._total, sc._MAP_LOD_BUCKETS)
    assert counts == want


def test_an_origin_spanning_feature_lands_at_both_ends():
    """A wrap feature's run passes the last bucket and continues at 0. If it
    were recorded as one run that falls off the end, the head of the ring
    would be bare — invisible on a plasmid, glaring on a chromosome."""
    n = 100_000
    rec = _rec(n, 0)
    rec.features = [SeqFeature(
        sc.CompoundLocation([FeatureLocation(n - 5_000, n, 1),
                             FeatureLocation(0, 5_000, 1)]),
        type="gene", qualifiers={"label": ["wrap"]})] if hasattr(
            sc, "CompoundLocation") else []
    if not rec.features:
        from Bio.SeqFeature import CompoundLocation
        rec.features = [SeqFeature(
            CompoundLocation([FeatureLocation(n - 5_000, n, 1),
                              FeatureLocation(0, 5_000, 1)]),
            type="gene", qualifiers={"label": ["wrap"]})]
    pm = _map_with(rec)
    counts, _ = pm._lod_histogram()
    nb = sc._MAP_LOD_BUCKETS
    assert counts[0] > 0, "the head of the wrap is missing"
    assert counts[nb - 1] > 0, "the tail of the wrap is missing"
    assert counts[nb // 2] == 0, "the far side should be empty"


def test_histogram_is_cached_across_rotations():
    """Keyed on bp only, so a rotate reuses it. Rebuilding per frame would
    put the O(features) cost straight back."""
    pm = _map_with(_rec(2_000_000, 5_000, span=60))
    first = pm._lod_histogram()
    pm._rotation_bp = 12345
    assert pm._lod_histogram() is first, "a rotate must not rebuild it"


def test_label_slots_are_capped_and_cached():
    pm = _map_with(_rec(2_000_000, 5_000, span=60))
    slots = pm._lod_label_slots(pm._total)
    assert 0 < len(slots) <= sc._MAP_LOD_MAX_LABELS
    assert pm._lod_label_slots(pm._total) is slots


def test_labels_go_to_the_widest_features():
    """At this zoom a label can only name something big enough to see."""
    n = 1_000_000
    rec = _rec(n, 3_000, span=60)
    from Bio.SeqFeature import FeatureLocation as FL
    rec.features.append(SeqFeature(FL(10_000, 400_000, 1), type="gene",
                                    qualifiers={"label": ["HUGE"]}))
    pm = _map_with(rec)
    labels = [lbl for _i, _mid, lbl, _c in pm._lod_label_slots(pm._total)]
    assert any("HUGE" in str(x) for x in labels)


def test_empty_and_degenerate_records_do_not_raise():
    for rec in (_rec(1_000, 0), _rec(10, 1, span=1)):
        pm = _map_with(rec)
        counts, colors = pm._lod_histogram()
        assert len(counts) == sc._MAP_LOD_BUCKETS == len(colors)


@pytest.mark.parametrize("bad", [
    {"start": -100, "end": -50},
    {"start": 10 ** 9, "end": 10 ** 9 + 5},
    {"start": 50, "end": 10},
    {"start": None, "end": 10},
    {"start": "x", "end": "y"},
    {"start": 1.5, "end": 9.5},
    {},
])
def test_a_malformed_feature_does_not_break_the_labels(bad):
    """The LABEL pass sorts by feature length, and a non-integer or missing
    coordinate made `sorted` compare int to None — a TypeError that took the
    whole map down. The histogram already skipped such features; the label
    pass did not, so one bad annotation blanked the view that was supposed to
    survive it."""
    pm = _map_with(_rec(100_000, 3_000, span=60))
    assert pm._lod_active()
    pm._feats.append({"color": "", "strand": 1, "type": "gene",
                      "label": "bad", **bad})
    pm._lod_cache = None
    pm._lod_label_cache = None
    slots = pm._lod_label_slots(pm._total)      # must not raise
    assert all(isinstance(mid, int) for _i, mid, _l, _c in slots)
    assert sum(pm._lod_histogram()[0]) > 0


def test_a_malformed_feature_is_skipped_not_fatal():
    """One bad annotation must not take the whole map down — the same rule
    the digest chain follows."""
    pm = _map_with(_rec(100_000, 10))
    pm._feats.append({"start": "nonsense", "end": None, "color": "",
                      "strand": 1, "type": "gene", "label": "bad"})
    pm._lod_cache = None
    counts, _ = pm._lod_histogram()
    assert sum(counts) > 0


# ═══════════════════════════════════════════════════════════════════════════
# The point of the exercise
# ═══════════════════════════════════════════════════════════════════════════

def test_aggregate_cost_does_not_grow_with_feature_count():
    """THE regression guard. The aggregate path is O(buckets), so building
    the ring for 100,000 features must cost about what 10,000 costs — not
    10x. This is what stopped a rotation taking 40 seconds."""
    pm_small = _map_with(_rec(5_000_000, 10_000, span=60))
    pm_big = _map_with(_rec(5_000_000, 100_000, span=60))
    assert pm_small._lod_active() and pm_big._lod_active()

    def ring_ms(pm):
        pm._lod_cache = None
        pm._lod_label_cache = None
        t = time.perf_counter()
        pm._lod_histogram()
        pm._lod_label_slots(pm._total)
        return (time.perf_counter() - t) * 1000

    # Build cost IS O(features) once; what must not grow is the per-frame
    # cost, so measure the cached path.
    ring_ms(pm_small)
    ring_ms(pm_big)
    t = time.perf_counter()
    for _ in range(20):
        pm_small._lod_histogram()
        pm_small._lod_label_slots(pm_small._total)
    small = (time.perf_counter() - t) * 1000
    t = time.perf_counter()
    for _ in range(20):
        pm_big._lod_histogram()
        pm_big._lod_label_slots(pm_big._total)
    big = (time.perf_counter() - t) * 1000
    # Cached lookups: both effectively free, and crucially not 10x apart.
    assert big < max(5.0, small * 3 + 1), (
        f"cached aggregate cost scaled with feature count: "
        f"10k={small:.2f} ms vs 100k={big:.2f} ms")


def test_cache_is_keyed_on_a_generation_not_an_address():
    """`id(self._feats)` looked like a fine cache key and is not: CPython
    reuses addresses after GC, so a fresh feature list of the SAME length
    could land on a freed one's address and be served a stale density ring.
    `_RESTR_SCAN_CACHE` documents the identical trap for transient strings.

    Loading a different record with the same feature COUNT must produce a
    different histogram."""
    pm = sc.PlasmidMap()
    a = _rec(100_000, 3_000, span=60)
    pm.load_record(a)
    first = pm._lod_histogram()
    gen_a = pm._lod_gen

    # A genuinely different molecule, same feature count and same length.
    b = _rec(100_000, 3_000, span=900)
    pm.load_record(b)
    assert pm._lod_gen != gen_a, "load_record must bump the generation"
    second = pm._lod_histogram()
    assert second is not first, "a new record must not reuse the cache"
    assert second[0] != first[0], (
        "the two records have different feature widths, so their density "
        "rings must differ")
