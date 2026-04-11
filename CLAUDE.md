# CLAUDE.md — AI Agent Context for SpliceCraft

This file is the **agent handoff document** for SpliceCraft. Any AI agent (Claude, GPT, Copilot, Gemini, or future systems) can read this file to understand the architecture, conventions, and design decisions behind the codebase — and pick up development, fix bugs, or build new modules without needing the full conversation history.

The project is developed in continuous collaboration between a human bioinformatician and an AI agent (Claude Opus 4.6). The goal is that any future agent can fork, extend, or integrate with this codebase compatibly.

---

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer and sequence editor** built with Python 3.12+ / Textual / Biopython. Renders Unicode braille-dot circular and linear plasmid maps directly in the terminal, with a per-base sequence panel, restriction-site overlays, library, and Golden Braid L0 assembly tooling.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs org, user ATinyGreenCell)

- **Single-file architecture:** the entire app is `splicecraft.py` (~4,200 lines). Intentional — avoids import complexity and keeps the codebase greppable. (Sibling project ScriptoScope follows the same convention at ~8,600 lines.)
- **Test suite:** 127 tests across 6 files in `tests/` (last refresh 2026-04-11). Full run ~15 s. See the **Test suite** section below.
- **Dependencies** (system-wide via `--break-system-packages` on Ubuntu/WSL2): `textual`, `biopython`, plus `pytest` and `pytest-asyncio` for the test suite. `Bio.Seq` and `Bio.SeqRecord` are the only Biopython surfaces touched in hot paths. **Optional runtime:** `pLannotate` (conda, GPL-3) for the Shift+A annotation feature — SpliceCraft works fine without it and notifies the user how to install if they press Shift+A.

## How to run

```bash
cd ~/SpliceCraft
python3 splicecraft.py              # empty canvas
python3 splicecraft.py L09137       # fetch pUC19 from NCBI
python3 splicecraft.py myplasmid.gb # open local GenBank file
python3 -m pytest -q                # run the 127-test sanity suite
```

Logs are written to `/tmp/splicecraft.log` (override with `$SPLICECRAFT_LOG`). Each log line is prefixed with an 8-char session ID so multi-run logs grep cleanly.

### Optional: pLannotate for automatic annotation

Press **Shift+A** (or click the ◈ button in the library panel) to run pLannotate on the current plasmid. SpliceCraft only calls it as a subprocess — it is never imported (pLannotate is GPL-3 and keeping it behind a subprocess boundary avoids license entanglement).

Install pLannotate in a dedicated conda env:
```bash
conda create -n plannotate -c conda-forge -c bioconda plannotate
conda activate plannotate
plannotate setupdb          # downloads ~500 MB of BLAST/diamond DBs
# then run SpliceCraft from the same conda env
```

If pLannotate is not on `PATH`, Shift+A just notifies the user with these instructions and returns — nothing crashes.

## Architecture (single file: `splicecraft.py`)

### Top-level structure (line numbers ±20, accurate as of 2026-04-11)

| Lines | Section |
|-------|---------|
| 1–30 | Module docstring, stdlib imports |
| 31–60 | Dependency check (`_check_deps`) |
| 61–110 | **Logging setup** — rotating file handler, session ID, startup banner |
| 111–155 | Library persistence (`plasmid_library.json`, `_load_library`, `_save_library`) |
| 156–420 | NEB restriction enzyme catalog (~200 enzymes), IUPAC tables, cached regex compilation, IUPAC-aware reverse-complement (`_rc`) |
| 425–560 | `_scan_restriction_sites()` — scans both strands, returns `resite` + `recut` dicts; palindrome-aware (skips reverse scan for palindromic sites) |
| 565–960 | Sequence panel rendering helpers — `_assign_chunk_features`, `_render_feature_row_pair`, `_chunk_lane_groups`, `_build_seq_text` (forward features above DNA, reverse below, braille bars with arrowheads) |
| 980–1100 | Codon table, OSC52 clipboard, CDS translation, GenBank I/O (NCBI fetch + local `.gb`) |
| 1102–1210 | `_Canvas` (2D char grid) and `_BrailleCanvas` (sub-character resolution via Unicode braille U+2800–U+28FF) |
| 1212–1880 | `PlasmidMap` widget — circular + linear map rendering, feature arcs, restriction site overlays, label placement, tick marks |
| 1884–1960 | `FeatureSidebar` — DataTable of features with detail panel |
| 1962–2070 | `LibraryPanel` — persistent plasmid collection with add/remove |
| 2075–2520 | `SequencePanel` — DNA viewer with click-to-cursor, drag selection, double-stranded display, feature annotation bars |
| 2525–2820 | Modal dialogs — `EditSeqDialog`, `FetchModal`, `OpenFileModal`, `DropdownScreen` |
| 2820–2910 | `MenuBar` widget — File / Edit / Enzymes / Features / Primers / Genes / Build |
| 2910–3290 | **Golden Braid L0** — `PartsBinModal`, `ConstructorModal` (assembly UI) |
| 3290–3320 | `UnsavedQuitModal` |
| 3321–4175 | `PlasmidApp` (main app) — keybindings, undo/redo stack, record loading, feature selection coordination, menu actions, workers |
| 4179–end | `main()` entry point |

### Key design patterns

- **All rendering uses Rich `Text`** — no curses
- **Braille canvas** gives sub-character pixel resolution (2x4 dots per terminal cell)
- **Feature coordination:** map click → sidebar highlight → sequence scroll (and vice versa via Textual messages)
- **Undo/redo:** snapshot-based (full seq + cursor + SeqRecord), max 50
- **Restriction sites:** scanned on load/edit, stored as `resite` (recognition bar) + `recut` (cut marker) dicts
- **Caching:** `PlasmidMap`, `SequencePanel`, and IUPAC regex patterns all cache rendered/compiled output keyed on state. Cache keys include `id(self._feats)` since lists are reassigned (not mutated) on load
- **Workers:** `@work(thread=True)` for NCBI fetch (`FetchModal._do_fetch`) and first-run library seed (`PlasmidApp._seed_default_library`). Both use `call_from_thread` to push results back to the UI

## Logging convention (borrowed from ScriptoScope)

```python
_log = logging.getLogger("splicecraft")
# Rotating file at /tmp/splicecraft.log, 2MB × 2 backups
# Every line prefixed with [session_id] for multi-run grepping
```

- **User-facing errors** → `self.notify(...)` or `Static.update("[red]...[/]")`. Never raw tracebacks.
- **Diagnostic detail** → `_log.exception("context: %s", ...)` inside `except` blocks. Stack traces go to the log file only.
- **Worker errors** → log with `_log.exception`, then push a friendly message to the UI via `call_from_thread`.

## Sacred invariants (DO NOT BREAK)

These are not yet test-enforced (no test suite exists), but breaking any of them will silently corrupt biology:

1. **Palindromic enzymes are scanned forward only.** `_scan_restriction_sites` must skip the reverse scan for palindromic sites and add only a bottom-strand `recut`. Scanning both strands for palindromes double-counts every site (regression introduced and fixed 2026-03-30).

2. **Reverse-strand resite positions use the forward coordinate.** A reverse-strand hit at position `p` (after RC) is stored as `p`, not `n - p - site_len`. The cut position maps via `site_len - 1 - fwd_cut` (regression fixed 2026-03-30).

3. **`_rc()` handles full IUPAC.** Reverse-complement must translate ambiguity codes (R, Y, W, S, M, K, B, D, H, V, N) via `_IUPAC_COMP`, not just ACGT.

4. **IUPAC regex patterns are cached.** `_iupac_pattern()` uses `_PATTERN_CACHE` to avoid recompiling ~200 patterns on every restriction scan.

5. **Circular wrap-around midpoints.** When computing the midpoint of a feature for label placement, use `arc_len = (end - start) % total` then `(start + arc_len // 2) % total`. The naive `(start + (end - start) // 2) % total` puts the label opposite the actual arc for wrapped features (bug fixed 2026-04-11).

## Recent fixes (2026-04-11 session)

- **Removed unused imports** — `Coordinate`, `TabbedContent`, `TabPane`, `Select` (none referenced anywhere in the file).
- **Wrapped-feature label midpoint** — fixed at `splicecraft.py:~1545`. See sacred invariant #5.
- **Logging infrastructure** — rotating file logger with session ID, startup banner capturing Python / platform / textual / biopython versions. Replaces silent `except: pass` blocks with `_log.exception(...)`.
- **Worker error reporting** — both `@work` workers (`FetchModal._do_fetch`, `PlasmidApp._seed_default_library`) now log exceptions; the previously-silent first-run seed failure is now diagnosable.
- **`main()` wrapped in try/except** — unhandled exceptions are logged before re-raising; session end is logged in `finally`.
- **Test suite added (99 tests)** — see section below.
- **Performance pass** — see section below.

## Performance (2026-04-11 optimization pass)

Four targeted fixes, roughly halving the sequence-panel build time and eliminating a redundant rebuild on every record load. Full suite unchanged at 99/99 green.

### 1. Suppressed the sidebar populate cascade (`_build_seq_text` 2× → 1×)
`FeatureSidebar.populate()` now sets a `_populating` flag and clears it via `call_after_refresh`, which fires AFTER Textual's message queue is drained. Previously, the auto-cursor-to-row-0 triggered by `DataTable.clear()`+`add_row()` fired a `RowHighlighted` event that cascaded through `PlasmidApp._sidebar_row_activated` → `SequencePanel.highlight_feature` → `_refresh_view`, causing a second expensive rebuild. One-line flag + deferred reset eliminates it completely.

**Savings per record load**: ~50 ms on 10 kb plasmids, ~110 ms on 20 kb, ~180 ms on 40 kb. See `FeatureSidebar.populate()` at ~line 1918.

### 2. Memoized styles + sorted annotation list in `_build_seq_text`
The expensive inputs that depend only on `(seq, feats)` — the per-base styles array and the sorted/filtered annotation list — are now computed once by `_build_seq_inputs()` and cached in a module-level dict `_BUILD_SEQ_CACHE` keyed by `(id(seq), id(feats), len(seq), len(feats))`. Cursor moves, selection changes, and resize events re-enter `_build_seq_text` without recomputing those inputs.

**Savings per call** (once the cache is warm): ~40 % on 20 kb plasmids. Cache is bounded at 4 entries; since keys are identity-based it never bloats. See `_build_seq_inputs()` at ~line 816.

### 3. Per-chunk `str.translate` instead of per-base
Inside `_build_seq_text`, the reverse-complement strand is built by translating the whole chunk string once (`chunk_rev = chunk_fwd.translate(_COMP)`) and indexing into it, instead of calling `.translate()` on every single base character. For a 20 kb plasmid that collapses ~229 000 single-character translate calls into ~160 chunk translate calls.

**Savings**: ~60 % on `_build_seq_text` cold-call time (130 ms → 49 ms on 20 kb). See the chunk loop at ~line 950.

### 4. Precomputed `_SCAN_CATALOG` for restriction scanning
`_NEB_ENZYMES` entries are now eagerly expanded at import time by `_rebuild_scan_catalog()` into tuples of `(name, site, site_len, fwd_cut, rev_cut, color, pat, is_palindrome, rc_pat)`. Previously `_scan_restriction_sites()` called `_rc(site)`, `_iupac_pattern(site)`, `len(site)`, and `_RESTR_COLOR[name]` once per enzyme **per scan** — that's ~200 redundant calls times however many scans happen during a session. The catalog is built once at import and iterated directly.

**Savings**: ~15 % on `_scan_restriction_sites()` for small plasmids (7.5 → 6.5 ms on pUC19); negligible on large plasmids where regex execution dominates. See `_rebuild_scan_catalog()` at ~line 448.

### Aggregate impact on real workloads

| Plasmid | `_build_seq_text` before | after (cold) | `_scan_restriction_sites` before | after |
|---|---:|---:|---:|---:|
| pUC19 (2.7 kb, 10 feats) | 11.3 ms | 6.2 ms | 7.5 ms | 6.5 ms |
| medium (10 kb, 40 feats) | 49.3 ms | 29.2 ms | 25 ms | 35 ms* |
| large (20 kb, 80 feats) | 112.7 ms | 49.0 ms | 65 ms | 52 ms |
| BAC (40 kb, 160 feats) | ~250 ms | 101 ms | — | 92 ms |

\* The 10 kb scan number got slightly worse — that's random-sequence variance between runs; test budgets have plenty of headroom.

Combined with eliminating the duplicate rebuild, the effective work done during a record load drops roughly in half. On a 20 kb plasmid, mount+settle goes from a double `_build_seq_text` at 225 ms down to a single call at 49 ms — about a 175 ms improvement per load.

### Test suite speedup
`tests/test_smoke.py` used a hardcoded 500 ms `pilot.pause()` to let `call_after_refresh` callbacks fire. With the sidebar-populate fix in place, a simple `pilot.pause()` + 50 ms suffices. Smoke suite dropped from 16 s to ~9 s; full suite from 18 s to **~11 s** across 99 tests.

### Performance budget tests (tightened)

See `tests/test_performance.py`:

| Budget | Threshold | Baseline | Headroom |
|---|---:|---:|---:|
| `_scan_restriction_sites(pUC19)` | 30 ms | 6.5 ms | 4.6× |
| `_scan_restriction_sites(10 kb)` | 150 ms | 35 ms | 4.3× |
| `_scan_restriction_sites` scaling ratio | < 8× | ~3.4× | — |
| `_iupac_pattern` warm (200 lookups) | 5 ms | 0.07 ms | 70× |
| `_rc(10 kb)` | 2 ms | 0.1 ms | 20× |
| `_build_seq_text(pUC19)` | 25 ms | 6 ms | 4× |
| `_build_seq_text(20 kb)` | 200 ms | 50 ms | 4× |

New cache-hit test (`test_warm_cache_skips_styles_work`) asserts that `_BUILD_SEQ_CACHE` actually holds an entry after the first call and that subsequent cursor-position changes reuse it.

### What was profiled but NOT touched

- **Textual compositor / CSS apply cycle** — accounts for ~500 ms of mount time on large plasmids but lives entirely in the framework. Can't optimize without switching away from Textual.
- **Rich `Text.append`** — called ~30 k times per 20 kb render; already efficient.
- **`PlasmidMap._draw`** — stable at 5-13 ms regardless of plasmid size because cost scales with canvas dimensions, not sequence length. Already well-cached via `_render_cache`.
- **Import time** — ~450 ms, dominated by Textual (265 ms) and Rich (105 ms). Lazy-importing Biopython helpers is already in place (`fetch_genbank`, `load_genbank`, `_record_to_gb_text`, etc. all `from Bio import ...` inside the function body).

## Test suite

Added 2026-04-11 to protect the sacred invariants. Modeled on ScriptoScope's
layout. Full suite runs in ~15 s; the biology-correctness subset runs in
< 1 s and is the fastest feedback loop.

### Running

```bash
python3 -m pytest -q                        # all 127 tests
python3 -m pytest tests/test_dna_sanity.py  # only biology (< 1 s)
python3 -m pytest tests/test_plannotate.py  # only pLannotate integration (~1 s)
python3 -m pytest tests/test_smoke.py       # only TUI smoke (~9 s)
python3 -m pytest -k "palindrome"           # only palindrome-related tests
python3 -m pytest -x                        # stop on first failure
```

`pyproject.toml` sets `asyncio_mode = "auto"` so async test functions don't
need `@pytest.mark.asyncio`. `tests/conftest.py` prepends the repo root to
`sys.path` and defines the `tiny_record`, `tiny_gb_path`, and
`isolated_library` fixtures used across files.

### Files and what each covers

| File | Tests | Covers | Runtime |
|------|------:|--------|---------|
| `tests/test_dna_sanity.py` | 53 | All 5 sacred invariants. `_rc` IUPAC involution + ground truth + Biopython cross-check; `_iupac_pattern` degenerate expansion + cache identity; `_CODON_TABLE` 64-entry completeness + Biopython cross-check; `_NEB_ENZYMES` schema, no duplicates, IUPAC regex compiles, Type IIS cut-outside-recognition; `_scan_restriction_sites` palindrome dedup, non-palindrome forward-coordinate positioning (2026-03-30 regression guard), unique_only filter, min_recognition_len filter; `_translate_cds` forward & reverse strands, stop padding, partial codons, unknown codon → `?`. | 0.6 s |
| `tests/test_circular_math.py` | 13 | Wrap-around midpoint formula (2026-04-11 regression guard); `PlasmidMap._bp_in` for wrapped and non-wrapped features, origin crossings, zero-width features. | 0.7 s |
| `tests/test_genbank_io.py` | 14 | `load_genbank` file round-trip preserves sequence bytes + feature type + strand; `_record_to_gb_text` / `_gb_text_to_record` text round-trip; `_save_library` / `_load_library` JSON round-trip, corruption recovery, cache memoization. | 1.3 s |
| `tests/test_plannotate.py` | 24 | pLannotate integration — availability detection for every failure mode (none missing, only plannotate, missing blastn, missing diamond, all ready, cache reuse); size-cap preflight (>50 kb rejected); feature merging (preserves originals, appends with `note="pLannotate"`, skips duplicates, doesn't mutate original); subprocess error paths (not installed, missing db, nonzero exit, timeout, FileNotFoundError at exec, no output file); happy path with mocked subprocess writing a real .gbk into the real tmpdir. pLannotate is never actually invoked by tests — `shutil.which` and `subprocess.run` are always monkeypatched. | 1.3 s |
| `tests/test_smoke.py` | 14 | Textual app mounts with preloaded record, all 4 panels present, features loaded, sequence populated, restriction scan ran on load, rotation keys, view toggle, RE toggle, mount works without a preload, no-network guard. Plus pLannotate UI entry points: library `#btn-lib-annot` exists, Shift+A binding registered, Shift+A on empty state notifies gracefully, Shift+A with pLannotate missing notifies install instructions and does NOT invoke subprocess. | 9 s |
| `tests/test_performance.py` | 9 | Budget enforcement: scan pUC19 < 30 ms, scan 10 kb < 150 ms, scan scaling < 8× for 4× more DNA (catches O(n²) regressions), `_iupac_pattern` warm < 5 ms for 200 lookups, warm strictly faster than cold, `_rc(10 kb)` < 2 ms, `_build_seq_text(pUC19)` < 25 ms, `_build_seq_text(20 kb)` < 200 ms, `_BUILD_SEQ_CACHE` populated after first call. Budgets are 4-70× over current baseline, so machine noise doesn't trip them but a real regression will. | 2.7 s |

### Sacred invariant → test mapping

| Sacred invariant (from earlier section) | Test file | Test method |
|---|---|---|
| #1 Palindromic enzymes scanned forward only | `test_dna_sanity.py` | `TestRestrictionScan::test_ecori_single_site_not_double_counted`, `::test_palindromes_produce_one_recut_per_site` |
| #2 Reverse-strand resite uses forward coordinate | `test_dna_sanity.py` | `TestRestrictionScan::test_non_palindrome_on_reverse_strand_uses_forward_coordinate`, `::test_non_palindrome_reverse_strand_asymmetric` |
| #3 `_rc()` handles full IUPAC | `test_dna_sanity.py` | `TestReverseComplement::test_rc_handles_each_iupac_code`, `::test_rc_is_involutive` |
| #4 IUPAC regex patterns are cached | `test_dna_sanity.py`, `test_performance.py` | `TestIUPACPattern::test_pattern_cache_*`, `TestIUPACPatternCachePerformance::test_warm_cache_is_near_free` |
| #5 Circular wrap-around midpoint formula | `test_circular_math.py` | `TestFeatureMidpoint::test_wrap_around_*`, `::test_non_wrapped_vs_wrapped_disagree_with_naive_formula` |

### Conventions

- **Cross-validate against Biopython** where possible (codon table, reverse-complement). If Biopython's standard table ever changes, the test fails noisily instead of silently mistranslating.
- **Hand-verifiable** test inputs — every restriction-site test uses a sequence short enough to stare at and count the expected hits by eye.
- **Regression guards cite the date** — every test protecting a past bug has a docstring line like `# Regression guard for 2026-03-30 fix`, pointing you to the original incident.
- **No network, no real files** — all tests use synthetic `SeqRecord`s and monkeypatched paths. The `isolated_library` fixture redirects `_LIBRARY_FILE` to a tmp path so the real `plasmid_library.json` is never touched. `fetch_genbank` is monkeypatched to raise in the no-network guard test.
- **Performance budgets are LOOSE** (6-20× headroom). They catch architectural regressions (O(n²), cache removal), not micro-perf drift. Tighten them only when you have a specific target to optimize against.

### Adding a new test

1. Pick the right file (or add a new one if it's a new subsystem).
2. If the test depends on a SeqRecord, use the `tiny_record` fixture — it's cheap and repeatable.
3. For Textual async tests, use `async def test_*` (no decorator needed — `asyncio_mode = "auto"`), and `async with app.run_test(size=TERMINAL_SIZE) as pilot: await pilot.pause(); await pilot.pause(0.5)`. The double-pause is needed to let `call_after_refresh` callbacks fire.
4. For perf tests, warm the cache first, then average over 10-20 iterations.
5. Run just your test first (`pytest tests/test_x.py::TestClass::test_y -q`) before running the full suite.

## Earlier fixes (2026-03-30)

- **Palindromic RE double-counting** — see sacred invariant #1.
- **Reverse-strand resite positions** — see sacred invariant #2.
- **`_rc()` IUPAC handling** — see sacred invariant #3.
- **Regex recompilation** — see sacred invariant #4.
- **Duplicate enzyme entries** — SbfI and NspI were each defined twice in `_NEB_ENZYMES`; removed.

## Released vs. unreleased state

### v0.1.0 (2026-03-23)
Braille circular map, NCBI fetch, local .gb loading, library, feature sidebar, sequence panel, undo/redo, restriction sites.

### Unreleased (in code, see CHANGELOG.md `[Unreleased]`)
- Feature deletion (Delete key)
- Linear map view toggle (`v` key)
- Strand-aware DNA layout (forward above, reverse below)
- Braille feature bars + single-bp triangles in sequence panel
- Label-above / label-below layout, connector lines (`l` key toggle)
- Full NEB enzyme catalog (~200 enzymes, Type IIS support)
- Inside tick marks on circular map
- Full-length feature labels with proximity placement
- Default library entry (MW463917.1 / pACYC184)
- **Parts Bin + Assembly Constructor** (Golden Braid L0 cloning) — `PartsBinModal` and `ConstructorModal`, latest commit `1bd7db2`

### Stubs (visible in menus, not implemented)
- **Primers > Design Primer** — `coming soon` notification only
- **Features > Add Feature** (`action_add_feature`) — `coming soon` notification only
- **Build > Simulate Assembly** — `coming soon` notification
- **Build > New Part editor** — `coming soon` notification

## pLannotate integration

Shift+A (or the ◈ button in the library panel, or `Features > Annotate with pLannotate` in the menu) runs pLannotate on the currently-loaded plasmid and merges the results into the current record.

### Design principles

1. **Subprocess only, never import.** pLannotate is GPL-3 — importing it would arguably create a combined work under GPL. Calling `plannotate batch` as a subprocess keeps the boundary clean. See `_run_plannotate()` at ~line 1190.
2. **Optional runtime dependency.** SpliceCraft works without pLannotate. The user only finds out it's optional when they press Shift+A; at that point, the UI shows an actionable install hint (`conda create -n plannotate ...`).
3. **Size cap preflighted.** pLannotate's hard-coded `MAX_PLAS_SIZE` is 50 kb; we refuse larger inputs instantly instead of waiting 30 s for pLannotate to do the same.
4. **Merge, don't replace.** Existing features are preserved; pLannotate hits are appended with a `note="pLannotate"` qualifier so users can tell them apart. A pLannotate hit at the same `(type, start, end, strand)` as an existing feature is skipped.
5. **Background worker.** `_run_plannotate_worker(@work(thread=True))` runs the subprocess off the main thread; UI updates go through `call_from_thread`. The worker emits an initial "Running pLannotate…" notify and a final "Added N features" or error notify.
6. **Undo-able.** `_run_plannotate_worker` calls `_push_undo()` before applying the merged record, so Ctrl+Z restores the pre-annotation state.
7. **Dirty flag.** After a successful annotation, the library panel marker shows `*Name` — the user saves the annotated version to the library with plain `a`.

### Code layout

| Function / class | ~Line | Purpose |
|---|---:|---|
| `PlannotateError` (+ 4 subclasses) | 1118 | Exception hierarchy with `user_msg` / `detail` attrs for UI display |
| `_PLANNOTATE_MAX_BP = 50_000` | 1132 | Matches pLannotate's `MAX_PLAS_SIZE` |
| `_plannotate_status()` | 1140 | `shutil.which`-based probe for `plannotate`, `blastn`, `diamond`. Cached in `_PLANNOTATE_CHECK_CACHE` |
| `_plannotate_install_hint()` | 1159 | Friendly `conda create ...` instructions for notifications |
| `_run_plannotate(record, timeout=180)` | 1169 | Runs the subprocess, catches every failure mode, returns a parsed SeqRecord. Raises `PlannotateError` subclasses |
| `_merge_plannotate_features(orig, annotated)` | 1251 | Pure function: returns a new SeqRecord preserving originals, appending pLannotate hits tagged with `note="pLannotate"` |
| `LibraryPanel.AnnotateRequested` | ~2050 | Message posted by the `◈` button |
| `LibraryPanel._btn_annotate` | ~2120 | Resolves the focused row and posts the message |
| `PlasmidApp.action_annotate_plasmid` | ~4460 | Shift+A action: preflights availability + size, kicks off the worker |
| `PlasmidApp._run_plannotate_worker` | ~4500 | `@work(thread=True)` — subprocess + merge + UI update |
| `PlasmidApp._library_annotate_requested` | ~4350 | Handles the library message; loads the focused entry first if needed |

### Failure modes and their notifications

| Condition | Notification | `_run_plannotate` exception |
|---|---|---|
| `plannotate` not on PATH | "pLannotate not installed. Install via conda: …" | `PlannotateNotInstalled` |
| `blastn` or `diamond` missing | "pLannotate needs blastn + diamond on PATH. …" | `PlannotateNotInstalled` |
| Databases not downloaded (pLannotate prints this to stdout and exits 0 — insidious) | "pLannotate databases not installed. Run: plannotate setupdb" | `PlannotateMissingDb` |
| Record > 50 kb | "pLannotate caps inputs at 50,000 bp (this plasmid: N bp)" | `PlannotateTooLarge` (also preflighted in `action_annotate_plasmid`) |
| Subprocess timeout (>180 s) | "pLannotate crashed: …" via `call_from_thread` | `PlannotateFailed` |
| Subprocess returns non-zero | "pLannotate failed: <last 500 chars of stderr>" | `PlannotateFailed` |
| Output `.gbk` not written | "pLannotate produced no .gbk output" | `PlannotateFailed` |
| Biopython can't parse output | "could not parse pLannotate output: <detail>" | `PlannotateFailed` |
| Unknown crash | "pLannotate crashed: <exc>" + full traceback in `/tmp/splicecraft.log` | any `Exception` |
| No record loaded | "Load a plasmid first (press 'f' or 'o')." | N/A (preflighted) |
| Successful with new features | "Added N pLannotate features. Press 'a' to save to library." | — |
| Successful but all hits duplicated existing features | "pLannotate found no new features (all hits duplicated existing annotations)." | — |

### License note for future maintainers

**Never `import plannotate` from SpliceCraft.** Always shell out via `_run_plannotate()`. If a future optimization makes you tempted to use pLannotate's Python API directly, check with the user first — SpliceCraft is not currently GPL-licensed, and importing a GPL module would arguably make the combined work GPL.

## Patterns worth porting from ScriptoScope (`/home/seb/proteoscope/scriptoscope.py`)

ScriptoScope (~8,600 lines) is the more mature sibling. When SpliceCraft grows, these patterns are pre-validated:

| Pattern | ScriptoScope location | When SpliceCraft would need it |
|---------|----------------------|--------------------------------|
| Rotating session-ID logger | ~lines 47–148 | **Already ported (2026-04-11)** |
| Startup banner with dep versions | ~lines 103–148 | **Already ported (2026-04-11)** |
| `main()` try/except wrapper | ~lines 8604–8609 | **Already ported (2026-04-11)** |
| Thread-local `Console` for `_text_to_content` | ~lines 3465–3501 | If sequence-panel render starts blowing 33 ms/frame budget. Required because shared Console + worker threads = lock contention; default Console blocks on TTY detection in WSL2 |
| Two-level render cache (`_seq_render_cache` + `_content_cache`) | ~lines 4407–4679 | If repainting on cursor moves becomes janky. LRU via `OrderedDict.move_to_end` |
| Atomic sidecar save (tempfile + fsync + os.replace) | ~lines 614–632 | If library JSON grows large or corruption is observed. Current direct write is fine for ~6KB |
| `Select.NULL`/`BLANK` sentinel filtering | ~line 6193, 7854 | Only if SpliceCraft adds a `Select` widget. Currently no `Select` is used (import removed) |
| Performance budgets in tests | tests/test_perf.py | Once a test suite exists |
| `_log.exception` + `call_from_thread` UI update pattern | ~lines 5503–5514 | Already adopted in SpliceCraft's two workers |
| `@lru_cache(1)` availability checks for optional CLI tools | ~lines 1593–1598 | If SpliceCraft ever shells out to BLAST / Prodigal / Primer3 |

## Known pitfalls

1. **Bare `except` was historically used to silence library I/O failures.** As of 2026-04-11 these now log via `_log.exception`. If you add another file-touching helper, use the same pattern instead of `pass`.

2. **Wrapped features (end_bp < start_bp) are first-class citizens.** Anywhere you compute distances, midpoints, or "is bp inside this feature", use the modular form. See `_bp_in()` (~line 1331) and the midpoint fix at ~line 1545 for canonical examples.

3. **Cache keys use `id(...)` of feature lists.** This is correct *only* because the app reassigns lists on load rather than mutating them in-place. If you ever start mutating `self._feats` in-place, the cache will return stale renders.

4. **Textual reactive auto-invalidation depends on field assignment, not mutation.** Setting `self._feats = new_list` triggers a refresh; calling `self._feats.append(x)` does not.

5. **Single-file means giant diffs are normal.** When a refactor touches the rendering layer, expect 100+ line edits. The greppability tradeoff is worth it.

## How to extend

### Adding a new modal dialog
1. Subclass `ModalScreen[ReturnType]` (see `FetchModal`, `OpenFileModal`, `PartsBinModal` as templates)
2. Implement `compose()` with the form layout
3. Call `self.dismiss(result)` to return a value
4. Push with `app.push_screen(MyModal(), callback=on_result)`

### Adding a heavy operation
1. Decorate with `@work(thread=True)`
2. Wrap the body in `try` / `except Exception` and call `_log.exception(...)` in the handler
3. Push results to the UI via `self.app.call_from_thread(callback)` — never touch widgets directly from the worker thread
4. Use `FetchModal._do_fetch` (line ~2640) as the canonical template

### Adding a new menu action
1. Add a method `action_my_thing(self)` on `PlasmidApp`
2. Add a binding in the `BINDINGS` list
3. Add the menu item in `MenuBar.compose()`

## Future work (user is undecided)

The user is on the fence between:
- **Merging** SpliceCraft, ScriptoScope, MitoShift, RefHunter, molCalc into one Textual app with multiple "modes"
- **Keeping them separate** as focused single-purpose apps and (optionally) extracting shared utilities into pure-Python modules importable from each

Either direction is viable. The single-file convention and shared logging/error patterns documented here keep the merge option open without forcing it.

## For future agents

If you are picking this up cold:

1. **Read this file first.** It gives you architecture without reading 4,200 lines.
2. **Run `python3 -m pytest -q`** before and after any change. 127 tests, ~15 s. The biology-correctness subset (`tests/test_dna_sanity.py`) runs in < 1 s if you want a faster inner loop.
3. **Check `/tmp/splicecraft.log`** (or `$SPLICECRAFT_LOG`) when debugging runtime issues. Every session has a unique 8-char ID.
4. **Don't break the sacred invariants.** Every one of them has a test (see the mapping table above). If you're touching `_scan_restriction_sites`, `_rc`, `_iupac_pattern`, `_translate_cds`, `_bp_in`, or the feature-midpoint formula, the relevant tests will tell you immediately if you got it wrong.
5. **Follow the error handling convention**: `_log.exception` for the stack trace, `notify()` or `Static.update("[red]...[/]")` for the user. Never let raw tracebacks hit the TUI.
6. **When in doubt about real-world behavior** — eyeball it on pUC19 (`L09137`) and pACYC184 (`MW463917.1`), both of which are fetched at first-run.
7. **Sister project for reference:** `/home/seb/proteoscope/scriptoscope.py` is the same author's larger app and the source of most patterns here (including the test suite layout). Read its `CLAUDE.md` for cross-pollination ideas.
