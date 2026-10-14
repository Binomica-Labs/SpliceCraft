#!/usr/bin/env python3
"""perf_probe — one-shot microbench for SpliceCraft hot paths.

Run from the repo root:  python3 scripts/perf_probe.py

Times the operations our perf audit flagged as suspect, and runs the
proposed alternatives side-by-side so we can pick fixes from data
rather than from agent reports. Not part of pytest. No assertions.

Every probe runs on synthetic data, and the data dir + log are sandboxed
BEFORE `import splicecraft` (CLAUDE.md's data-dir rules): the import alone
used to append to the real `logs/splicecraft.log`.
"""
from __future__ import annotations

import copy
import gc
import json
import os
import pickle
import random
import statistics
import sys
import tempfile
import time
from functools import lru_cache
from pathlib import Path

_SANDBOX = tempfile.mkdtemp(prefix="sc-perfprobe-")
os.environ["XDG_DATA_HOME"] = _SANDBOX
os.environ.pop("SPLICECRAFT_DATA_DIR", None)        # it outranks XDG_DATA_HOME
os.environ["SPLICECRAFT_LOG"] = os.path.join(_SANDBOX, "perf_probe.log")
os.environ.setdefault("SPLICECRAFT_SKIP_LOCK", "1")

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import splicecraft as sc  # noqa: E402

assert _SANDBOX in str(sc._state._DATA_DIR), f"unsandboxed: {sc._state._DATA_DIR}"


# ── timing harness ────────────────────────────────────────────────────────────

def bench(fn, *, iters: int = 5, warmup: int = 1) -> float:
    """Return median per-call ms across `iters` runs after `warmup`."""
    for _ in range(warmup):
        fn()
    samples = []
    for _ in range(iters):
        gc.collect()
        t0 = time.perf_counter()
        fn()
        samples.append((time.perf_counter() - t0) * 1000.0)
    return statistics.median(samples)


def header(title: str) -> None:
    print()
    print("═" * 78)
    print(f"  {title}")
    print("═" * 78)


def row(label: str, ms: float, *, per_call_us: float | None = None,
        note: str = "") -> None:
    bits = [f"  {label:<48}", f"{ms:>10.2f} ms"]
    if per_call_us is not None:
        bits.append(f"  ({per_call_us:>7.2f} µs/call)")
    if note:
        bits.append(f"  {note}")
    print("".join(bits))


# ── synthetic library builder ────────────────────────────────────────────────

def synth_gb_text(bp: int, n_features: int = 20, *, seed: int = 0) -> str:
    """A realistic GenBank-text-like blob: header + N features + sequence."""
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(bp))
    feats = "\n".join(
        f'     CDS             {1 + i * (bp // n_features)}..'
        f'{(i + 1) * (bp // n_features)}\n'
        f'                     /label="orf{i}"\n'
        f'                     /note="synthetic feature {i}"'
        for i in range(n_features)
    )
    body = "\n".join(
        f"{1 + i * 60:>9} " + " ".join(seq[i * 60:(i + 1) * 60][j:j + 10]
                                        for j in range(0, 60, 10))
        for i in range((bp + 59) // 60)
    )
    return (
        f"LOCUS       p{seed:04d}    {bp} bp ds-DNA circular SYN\n"
        f"FEATURES             Location/Qualifiers\n{feats}\n"
        f"ORIGIN\n{body}\n//\n"
    )


def synth_library(n: int, *, bp_per: int = 5_000) -> list[dict]:
    """Library of n entries each carrying gb_text of `bp_per` bp."""
    return [
        {
            "id": f"P{i:05d}",
            "name": f"plasmid_{i}",
            "gb_text": synth_gb_text(bp_per, seed=i),
            "notes": f"entry {i}",
            "topology": "circular",
            "length": bp_per,
            "primer_pairs": [],
        }
        for i in range(n)
    ]


# ── candidate alternatives to deepcopy ───────────────────────────────────────

_IMMUTABLE_TYPES = (str, int, float, bool, bytes, type(None))


def typed_clone(obj):
    """Deep-clone that shares immutables. Library entries hold only str /
    int / list / dict / tuple / None — strings (the bulk of gb_text
    payloads) don't need to be re-allocated since they're immutable.
    Tuples are recursively cloned only if they contain mutables."""
    t = type(obj)
    if t is dict:
        return {k: typed_clone(v) for k, v in obj.items()}
    if t is list:
        return [typed_clone(v) for v in obj]
    if t is tuple:
        return tuple(typed_clone(v) for v in obj)
    if isinstance(obj, _IMMUTABLE_TYPES):
        return obj
    return copy.deepcopy(obj)


def pickle_clone(obj):
    return pickle.loads(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))


def json_clone(obj):
    return json.loads(json.dumps(obj))


# ── lru-cached pure-function candidates ──────────────────────────────────────

@lru_cache(maxsize=8)
def _rc_cached(seq: str) -> str:
    return seq.upper().translate(sc._IUPAC_COMP)[::-1]


@lru_cache(maxsize=2048)
def _feat_len_cached(start: int, end: int, total: int) -> int:
    return (total - start) + end if end < start else end - start


# ── probes ───────────────────────────────────────────────────────────────────

def probe_library_clone() -> None:
    header("Library clone (`_load_library` returns deepcopy on every call)")
    print(f"  {'op':<48}{'median':>13}{'speedup':>20}")
    print("  " + "─" * 76)
    for n_entries, bp in [(10, 5_000), (100, 5_000), (500, 5_000),
                           (10, 100_000), (100, 100_000)]:
        lib = synth_library(n_entries, bp_per=bp)
        # current: deepcopy
        t_dc = bench(lambda: copy.deepcopy(lib), iters=5, warmup=1)
        # candidates
        t_tc = bench(lambda: typed_clone(lib), iters=5, warmup=1)
        t_pk = bench(lambda: pickle_clone(lib), iters=5, warmup=1)
        try:
            t_js = bench(lambda: json_clone(lib), iters=3, warmup=1)
        except (TypeError, ValueError):
            t_js = float("nan")

        tag = f"  n={n_entries}, gb_text={bp:>6} bp ({len(lib[0]['gb_text']):,} ch each)"
        print(f"\n{tag}")
        print(f"    deepcopy        {t_dc:>10.2f} ms      (baseline)")
        print(f"    typed_clone     {t_tc:>10.2f} ms      "
              f"{t_dc / max(t_tc, 1e-6):>5.1f}× faster")
        print(f"    pickle clone    {t_pk:>10.2f} ms      "
              f"{t_dc / max(t_pk, 1e-6):>5.1f}× faster")
        if t_js == t_js:  # not nan
            print(f"    json clone      {t_js:>10.2f} ms      "
                  f"{t_dc / max(t_js, 1e-6):>5.1f}× faster")


def probe_rc() -> None:
    header("`_rc(seq)` — current implementation vs lru_cache")
    for bp in [2_000, 10_000, 50_000, 200_000]:
        seq = "".join(random.choices("ACGT", k=bp))
        # cold: function call cost on a fresh string
        t_cold = bench(lambda: sc._rc(seq), iters=20, warmup=2)
        # warm: lru_cache hit
        _rc_cached.cache_clear()
        _rc_cached(seq)  # warm
        t_warm = bench(lambda: _rc_cached(seq), iters=200, warmup=10)
        print(f"  seq={bp:>6} bp  | call: {t_cold:>7.3f} ms  "
              f"| cached hit: {t_warm * 1000:>7.3f} µs  "
              f"| {t_cold / max(t_warm, 1e-9):>7.0f}× faster on hit")


def probe_feat_len() -> None:
    header("`_feat_len(start, end, total)` — uncached vs lru_cache")
    feats = [(random.randint(0, 100_000), random.randint(0, 100_000))
             for _ in range(2_000)]
    total = 100_000

    def call_native():
        for s, e in feats:
            sc._feat_len(s, e, total)

    def call_cached():
        for s, e in feats:
            _feat_len_cached(s, e, total)

    t_native = bench(call_native, iters=10, warmup=2)
    t_cached_cold = bench(call_cached, iters=1, warmup=0)
    t_cached_warm = bench(call_cached, iters=10, warmup=2)
    print(f"  2000 calls native       : {t_native:>7.3f} ms  "
          f"({t_native * 1000 / 2000:.2f} µs/call)")
    print(f"  2000 calls cached cold  : {t_cached_cold:>7.3f} ms  "
          f"({t_cached_cold * 1000 / 2000:.2f} µs/call)")
    print(f"  2000 calls cached warm  : {t_cached_warm:>7.3f} ms  "
          f"({t_cached_warm * 1000 / 2000:.2f} µs/call)")


def probe_scan_restriction() -> None:
    header("`_scan_restriction_sites(seq)` — current baseline")
    sc._scan_restriction_sites("ACGTACGT" * 100)  # warm regex cache
    for bp in [2_000, 10_000, 50_000, 200_000]:
        seq = "".join(random.choices("ACGT", k=bp))
        sc._scan_restriction_sites(seq)  # warm scan cache
        t = bench(lambda: sc._scan_restriction_sites(seq, unique_only=True),
                  iters=10, warmup=2)
        cache_size = len(sc._RESTR_SCAN_CACHE)
        print(f"  seq={bp:>6} bp  : {t:>8.2f} ms   "
              f"(scan cache holds {cache_size}/{sc._RESTR_SCAN_CACHE_MAX})")


def probe_safe_save_json() -> None:
    header("`_safe_save_json` — write latency at typical library sizes")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for n in [10, 100, 500]:
            lib = synth_library(n, bp_per=5_000)
            path = tmp / f"lib_{n}.json"
            t = bench(lambda: sc._safe_save_json(path, lib, "Probe library"),
                      iters=5, warmup=1)
            size_mb = path.stat().st_size / 1_000_000
            print(f"  n={n:>3} entries ({size_mb:>5.2f} MB)   : "
                  f"{t:>8.2f} ms")


def probe_sort_features() -> None:
    header("`PlasmidMap._draw` re-sorts features every paint — sort cost")
    for n in [20, 100, 500, 2_000]:
        feats = [
            {"start": random.randint(0, 100_000),
             "end": random.randint(0, 100_000),
             "label": f"feat_{i}"}
            for i in range(n)
        ]

        def do_sort():
            sorted(enumerate(feats), key=lambda iv: -sc._feat_len(
                iv[1]["start"], iv[1]["end"], 100_000))

        t = bench(do_sort, iters=20, warmup=2)
        print(f"  {n:>4} features  : {t:>7.3f} ms per sort  "
              f"(@60fps rotation = {t * 60:.1f} ms/sec budget burn)")


def probe_natural_sort() -> None:
    header("`_natural_sort_key` for library panel sort")
    for n in [50, 200, 1000]:
        names = [f"plasmid_{random.randint(1, 9999)}_v{random.randint(1, 9)}"
                 for _ in range(n)]

        def do_sort():
            sorted(names, key=sc._natural_sort_key)

        t = bench(do_sort, iters=20, warmup=2)
        print(f"  {n:>4} names     : {t:>7.3f} ms")


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> int:
    print(f"SpliceCraft perf_probe — {sc.__version__}")
    print(f"Python {sys.version.split()[0]} on {sys.platform}")

    probe_library_clone()
    probe_rc()
    probe_feat_len()
    probe_scan_restriction()
    probe_safe_save_json()
    probe_sort_features()
    probe_natural_sort()

    print()
    print("Done. Each row: median of 5–20 runs after warm-up. "
          "Lower is better.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
