# CLAUDE.md — AI Agent Context for SpliceCraft

This file is the **agent handoff document** for SpliceCraft. Any AI agent (Claude, GPT, Copilot, Gemini, or future systems) can read this file to understand the architecture, conventions, and design decisions behind the codebase — and pick up development, fix bugs, or build new modules without needing the full conversation history.

The project is developed in continuous collaboration between a human bioinformatician and an AI agent (Claude Opus 4.6). The goal is that any future agent can fork, extend, or integrate with this codebase compatibly.

---

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer and sequence editor** built with Python 3.12+ / Textual / Biopython. Renders Unicode braille-dot circular and linear plasmid maps directly in the terminal, with a per-base sequence panel, restriction-site overlays, library, and Golden Braid L0 assembly tooling.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs org, user ATinyGreenCell)

- **Single-file architecture:** the entire app is `splicecraft.py` (~7,100 lines). Intentional — avoids import complexity and keeps the codebase greppable. (Sibling project ScriptoScope follows the same convention at ~8,600 lines.)
- **Test suite:** 388 tests across 10 files in `tests/` (last refresh 2026-04-13). Full run ~75 s, biology subset (`test_dna_sanity.py`) < 1 s. See the **Test suite** section below.
- **Dependencies:** `textual>=8.2.3`, `biopython>=1.87`, `primer3-py>=2.3.0`, `platformdirs>=4.2`, plus `pytest>=9.0` / `pytest-asyncio>=1.3` for tests. `Bio.Seq` and `Bio.SeqRecord` are the only Biopython surfaces touched in hot paths. Users install via `pipx install splicecraft` (recommended on PEP 668 systems) or `pip install splicecraft` inside a venv. Developers working on the repo run `python3 splicecraft.py` directly from the clone. **Optional runtime:** `pLannotate` (conda, GPL-3) for the Shift+A annotation feature — SpliceCraft works fine without it and notifies the user how to install if they press Shift+A.
- **Published on PyPI** as `splicecraft`. Releases cut via `./release.sh X.Y.Z` (bumps version in both files, runs tests, builds, commits+tags+pushes; GitHub Actions `publish.yml` then publishes to PyPI via Trusted Publishing / OIDC). As of 2026-04-12 the latest published version is **0.2.2**.

## How to run

```bash
cd ~/SpliceCraft
python3 splicecraft.py              # empty canvas (dev repo, fastest feedback loop)
python3 splicecraft.py L09137       # fetch pUC19 from NCBI
python3 splicecraft.py myplasmid.gb # open local GenBank file
python3 -m pytest -q                # run the 250-test sanity suite

# End users (not devs) install from PyPI:
pipx install splicecraft            # or: pip install splicecraft (inside a venv)
splicecraft                         # same commands, no cd needed
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

### Top-level structure (line numbers ±30, accurate as of 2026-04-12)

| Lines | Section |
|-------|---------|
| 1–45 | Module docstring, stdlib imports |
| 46–100 | User data dir (`platformdirs`), legacy migration (`_user_data_dir`, `_migrate_legacy_data`) |
| 101–200 | Dependency check, rotating-file logger, session ID, startup banner |
| 200–330 | **Atomic JSON persistence** (`_safe_save_json` / `_safe_load_json` — tempfile + `.bak` + shrink guard), library loader/saver |
| 331–600 | NEB restriction enzyme catalog (~204 enzymes), IUPAC tables, cached regex, IUPAC-aware `_rc` |
| 600–760 | `_scan_restriction_sites()` — both strands, palindrome-aware, **circular wrap-around** (emits wrap sites as tail + head pieces; 2026-04-12) |
| 760–1055 | Sequence panel rendering helpers — `_assign_chunk_features`, `_render_feature_row_pair`, `_chunk_lane_groups`, `_build_seq_inputs` (memoized), `_build_seq_text` |
| 1055–1250 | OSC52 clipboard, CDS translation, GenBank I/O (`fetch_genbank`, `load_genbank`, record ↔ text round-trip) |
| 1350–1550 | **pLannotate integration** — subprocess-only, `_run_plannotate`, `_merge_plannotate_features`, error hierarchy |
| 1550–1660 | `_Canvas` + `_BrailleCanvas` (sub-character braille resolution) |
| 1660–2350 | `PlasmidMap` widget — circular + linear rendering, arcs, RE overlays, proximity label placement |
| 2350–2450 | `FeatureSidebar` — DataTable + detail panel |
| 2450–2620 | `LibraryPanel` — persistent plasmid collection with add/remove/rename/annotate |
| 2620–3065 | `SequencePanel` — DNA viewer, click-to-cursor, drag selection, double-stranded display |
| 3065–3375 | Modal dialogs — `EditSeqDialog`, `FetchModal`, `OpenFileModal`, `DropdownScreen` |
| 3375–3480 | `MenuBar` — File / Edit / Enzymes / Features / Primers / Parts / Constructor |
| 3480–3830 | **Primer design functions** — `_pick_binding_region`, `_design_gb_primers`, `_design_detection_primers`, `_design_cloning_primers_raw`, `_design_generic_primers` |
| 3830–4535 | **Golden Braid L0 UI** — `PartsBinModal`, `DomesticatorModal`, `ConstructorModal` |
| 4535–5440 | `PrimerDesignScreen` — Detection / Cloning / Golden Braid / Generic modes, primer library with mark/status/save |
| 5440–5660 | Small modals — `UnsavedQuitModal`, `PlasmidPickerModal`, `RenamePlasmidModal`, `LibraryDeleteConfirmModal` |
| 5660–7040 | `PlasmidApp` — main controller, keybindings, undo/redo, record loading, feature selection coord, menu actions, `@work` threads (NCBI fetch, seed, pLannotate), `_check_data_files` startup validation |
| 7040–end | `main()` entry point — `--version` / `--help` short-circuit, CLI arg parsing, try/except wrap for clean logging |

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

Every invariant below has at least one test protecting it. Breaking any will silently corrupt biology or user data. See the **Sacred invariant → test mapping** section below for specific test locations.

1. **Palindromic enzymes are scanned forward only.** `_scan_restriction_sites` must skip the reverse scan for palindromic sites and add only a bottom-strand `recut`. Scanning both strands for palindromes double-counts every site (regression introduced and fixed 2026-03-30).

2. **Reverse-strand resite positions use the forward coordinate.** A reverse-strand hit at position `p` (after RC) is stored as `p`, not `n - p - site_len`. The cut position maps via `site_len - 1 - fwd_cut` (regression fixed 2026-03-30).

3. **`_rc()` handles full IUPAC.** Reverse-complement must translate ambiguity codes (R, Y, W, S, M, K, B, D, H, V, N) via `_IUPAC_COMP`, not just ACGT.

4. **IUPAC regex patterns are cached.** `_iupac_pattern()` uses `_PATTERN_CACHE` to avoid recompiling ~200 patterns on every restriction scan.

5. **Circular wrap-around midpoints.** When computing the midpoint of a feature for label placement, use `arc_len = (end - start) % total` then `(start + arc_len // 2) % total`. The naive `(start + (end - start) // 2) % total` puts the label opposite the actual arc for wrapped features (bug fixed 2026-04-11).

6. **Circular wrap-around restriction scanning.** `_scan_restriction_sites(circular=True)` (default) scans `seq + seq[:max_site_len-1]` so recognition sequences that span the origin are found. Each wrap-around hit is emitted as **two resite pieces** — a labeled tail `[p, n)` plus an unlabeled head `[0, (p+site_len) - n)` — and **one recut** at the true absolute cut position `(p + fwd_cut) % n`. Downstream code that counts resites for filtering must count only labeled pieces (added 2026-04-12).

7. **Data-file saves always back up.** `_safe_save_json` writes a `.bak` of the existing file before replacing it, via `tempfile.mkstemp` + `os.fsync` + `os.replace`. If the about-to-be-written list is smaller than the existing one, a SHRINK GUARD warning is logged. Never bypass `_safe_save_json` by writing JSON directly — it is the user's only recovery path from a buggy save.

## Recent fixes (2026-04-13 session — unreleased)

Continued wrap-feature audit after the 2026-04-12 pass. 10 more verified bugs fixed across sequence editing, primer design, and circular-map rendering. Tests grew 361 → 388.

### Sequence editing
- **Wrap-feature mangling in `_rebuild_record_with_edit`.** `int(CompoundLocation.start)` returns `min(parts.start)` and `int(.end)` returns `max(parts.end)`, so reading wrap features as `(fs, fe) = (int(loc.start), int(loc.end))` silently flattened them into whole-plasmid FeatureLocations on every insert/replace. Now per-part shift with collapse-to-FeatureLocation when only 1 part survives. See `_rebuild_record_with_edit` at ~line 6675.
- **1-bp ghost stubs after delete.** `max(new_fs + 1, min(new_fe, new_len))` forced every post-edit feature to be at least 1 bp wide. Features fully consumed by a replace/delete survived as 1-bp stubs instead of being dropped. Clamp removed; zero-width features are now dropped.

### Undo/redo coherence
- **`_source_path` cleared on every in-place record change.** `_apply_record` unconditionally cleared `_source_path`, so after a pLannotate merge or primer-add, Ctrl+S no longer targeted the original file. The clear is now tied to `clear_undo=True` (i.e. fresh loads only).
- **Primer add-to-map aliased the undo snapshot.** `_add_selected_to_map` mutated `_current_record.features` directly; undo snapshots shared the object, so Ctrl+Z gave back a "pristine" state that wasn't. Now builds a fresh record via deepcopy and calls `_apply_record(clear_undo=False)`.

### Concurrent workers
- **pLannotate re-entry guard.** Second Shift+A while a first run is in flight could spawn parallel workers wasting 5-30s CPU. Added `_plannotate_running` flag with `finally` cleanup so only one run can be in flight at a time.
- **pLannotate merge missed `_mark_dirty`.** Worker called `lib.set_dirty(True)` (library panel marker) but not `self._unsaved=True`, so the quit prompt was silent on unsaved pLannotate edits. Now calls `_mark_dirty()` which flips both.

### Primer design on wrap regions (the big one)
- **`_design_gb/cloning/generic_primers` used `seq[start:end]` unconditionally.** Returned "" for wrap regions (end < start). Now use `_slice_circular(seq, start, end)` helper which returns `seq[start:] + seq[:end]` for wrap. Primer positions for wrap regions are computed via modular arithmetic: `(start + len(fwd_bind)) % total`, `(end - len(rev_bind)) % total`.
- **`_design_detection_primers` is Primer3-backed.** Primer3 is linear-only, so wrap regions now rotate the template to `seq[start:] + seq[:start]` before calling Primer3, then unrotate the returned positions back to original-template coordinates via `(coord + rotation) % total`.
- **Validators rejected wrap regions outright.** `_read_region_from` and `DomesticatorModal._design` both had `end <= start` as a "bad region" guard. Now accept `end < start` as a wrap indicator (but still reject `end == start` as empty).
- **`_add_selected_to_map` silently dropped wrap primers.** `if p_end <= p_start: continue` skipped every wrap primer with no notification. Now builds a `CompoundLocation([FL(p_start, total), FL(0, p_end)])` for wrap primers so they land on the map with both tail and head arcs visible.
- **`_feat_selected` computed feat_len as `end - start`.** Selecting a wrap feature in the primer-design dropdown gave a negative length, which the detection range auto-adjuster then set to bizarre values. Now uses `_feat_len(start, end, total)`.

### Circular-map rendering
- **Feature sort key used `end - start` for z-order.** Wrap features sorted to the *largest* bucket (due to negative key being most negative), rendering them on top of everything — often hiding smaller features they overlapped. Now all three sort sites (`PlasmidMap._draw`, `SequencePanel._click_to_bp`, `_annot_feats_sorted`, `_build_seq_inputs`) use `_feat_len(start, end, total)` which returns the correct biological length.
- **`FeatureSidebar.show_detail` displayed negative bp counts for wrap features.** Replaced `span = end - start` with `span = _feat_len(start, end, total)`.

### New helper
- **`_feat_len(start, end, total)`** — single source of truth for circular-aware feature length. Returns `(total - start) + end` when `end < start`, else `end - start`. All sort keys and length displays route through it.

## Recent fixes (2026-04-12 session — shipped as v0.2.2)

Six verified bugs fixed after a parallel review pass. Tests grew 246 → 250.

- **Golden Braid primer length validation.** `_design_gb_primers` silently returned garbage (binding < 18 bp, Tm = 0.0) when the selected region was too short. Now returns a clear `error` dict; `_run_goldenbraid` surfaces it in red in the results pane. The other three primer modes (detection / cloning / generic) already had this pattern.
- **pLannotate stale-record race.** If the user loaded a different plasmid while pLannotate was running, the worker's `_apply()` callback silently replaced the new plasmid with the merged OLD one. Now the callback checks `self._current_record is record` and discards the stale result with a warning.
- **Undo stack leaked across plasmid loads.** Ctrl+Z after switching plasmids yanked the user back to an unrelated edit. `_apply_record(clear_undo=True)` now clears stacks by default; in-place record changes (pLannotate merge, sequence edits) opt out via `clear_undo=False` so their `_push_undo()` entries remain intact.
- **Circular wrap-around restriction scanning.** Enzyme sites whose recognition sequence spanned the origin were invisible. See sacred invariant #6 for the exact contract — emits wrap sites as labeled tail + unlabeled head + one recut at `(p + fwd_cut) % n`.
- **Linear-view click matched zero-width features everywhere.** A malformed feature with `start == end` fell into the wrap branch `(bp >= s or bp <= e)`, always true. Linear click handler now calls `_bp_in(bp, f)` (half-open `[start, end)`), matching the circular-view behavior.
- **Shrink guard widened.** Previously only fired when writing zero entries to a populated file. Now fires on any shrink (`len(entries) < existing_count`), still logs-not-blocks — a deletion is legitimate but needs to be visible in `/tmp/splicecraft.log`.

## Earlier fixes (2026-04-11 session)

- **Removed unused imports** — `Coordinate`, `TabbedContent`, `TabPane`, `Select`.
- **Wrapped-feature label midpoint** (sacred invariant #5).
- **Logging infrastructure** — rotating file logger, session ID, startup banner.
- **Worker error reporting** — `FetchModal._do_fetch` and `PlasmidApp._seed_default_library` now log exceptions.
- **`main()` try/except wrap** — unhandled exceptions logged before re-raising; session end logged in `finally`.
- **Test suite first cut** (99 tests); performance pass (below).

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

Originally added 2026-04-11 to protect the sacred invariants; expanded in
every subsequent session. Full suite runs in ~55 s; the biology-correctness
subset runs in < 1 s and is the fastest feedback loop.

### Running

```bash
python3 -m pytest -q                        # all 388 tests
python3 -m pytest tests/test_dna_sanity.py  # only biology (< 1 s)
python3 -m pytest tests/test_plannotate.py  # only pLannotate integration (~1 s)
python3 -m pytest tests/test_primers.py     # only primer design (~2 s)
python3 -m pytest tests/test_domesticator.py # only Golden Braid (~2 s)
python3 -m pytest tests/test_smoke.py       # only TUI smoke (~12 s)
python3 -m pytest -k "palindrome"           # only palindrome-related tests
python3 -m pytest -x                        # stop on first failure
```

`pyproject.toml` sets `asyncio_mode = "auto"` so async test functions don't
need `@pytest.mark.asyncio`. `tests/conftest.py` prepends the repo root to
`sys.path`, defines `tiny_record` / `tiny_gb_path` / `isolated_library`
fixtures, and installs the **autouse** `_protect_user_data` fixture that
monkeypatches `_LIBRARY_FILE`, `_PARTS_BIN_FILE`, `_PRIMERS_FILE`, and their
caches to tmp paths for every test — no test can write to the real user
library files even by accident.

### Files and what each covers

| File | Tests | Covers |
|------|------:|--------|
| `tests/test_dna_sanity.py` | 74 | Sacred invariants 1-6 (RE scan, `_rc`, `_iupac_pattern`, codon table, **circular wrap scan**). `_NEB_ENZYMES` schema, no duplicates, IUPAC regex compiles, Type IIS cut-outside-recognition. `_translate_cds` forward & reverse strands, stop padding, partial codons, unknown codon → `?`. |
| `tests/test_circular_math.py` | 38 | Sacred invariant #5 (wrap midpoint). `PlasmidMap._bp_in` for wrapped / non-wrapped / zero-width features. `_feat_len` circular-aware length helper. |
| `tests/test_edit_record.py` | 14 | Regression guard for `_rebuild_record_with_edit`: wrap features survive insert/replace as CompoundLocation, collapse to FeatureLocation when only 1 part remains, features fully consumed by replace are dropped (not left as 1-bp stubs). Reverse-strand wraps keep strand. |
| `tests/test_genbank_io.py` | 59 | `load_genbank` round-trip preserves sequence + features + strand. `_save_library` / `_load_library` JSON round-trip, corruption recovery, cache memoization. |
| `tests/test_data_safety.py` | 28 | Sacred invariant #7 (atomic saves). `_safe_save_json` preserves `.bak`, `_safe_load_json` restores from `.bak` on corrupt main, `_protect_user_data` fixture confirmed working. |
| `tests/test_primers.py` | 60 | Detection primers (SEQUENCE_INCLUDED_REGION inside region), cloning primers (RE tails + GCGC pad), Golden Braid primers (BsaI domestication, overhang correctness), generic primers. **Wrap-region primer design** (`_slice_circular`, Primer3 template rotation, modular position mapping). |
| `tests/test_domesticator.py` | 41 | Golden Braid L0 positions and overhangs. Part validation, assembly lane construction, overhang compatibility. |
| `tests/test_plannotate.py` | 24 | Integration — availability detection, size-cap preflight, feature merging, subprocess error paths. pLannotate never actually invoked; `shutil.which` and `subprocess.run` monkeypatched. |
| `tests/test_smoke.py` | 41 | Textual app mounts with preloaded record, all 4 panels present, features loaded, restriction scan ran, rotation / view toggle / RE toggle work, mount without preload, no-network guard. pLannotate UI entry points + re-entry guard. Primer design screen mounts. `_apply_record` in-place semantics (clear_undo=False preserves source_path and undo). |
| `tests/test_performance.py` | 9 | Budget enforcement: scan pUC19 < 30 ms, scan 10 kb < 150 ms, scan scaling < 8× for 4× DNA, `_iupac_pattern` warm < 5 ms for 200 lookups, `_rc(10 kb)` < 2 ms, `_build_seq_text(pUC19)` < 25 ms, `_build_seq_text(20 kb)` < 200 ms, `_BUILD_SEQ_CACHE` populated after first call. |

### Sacred invariant → test mapping

| Sacred invariant | Test file | Test method |
|---|---|---|
| #1 Palindromic enzymes scanned forward only | `test_dna_sanity.py` | `TestRestrictionScan::test_ecori_single_site_not_double_counted`, `::test_palindromes_produce_one_recut_per_site` |
| #2 Reverse-strand resite uses forward coordinate | `test_dna_sanity.py` | `TestRestrictionScan::test_non_palindrome_on_reverse_strand_uses_forward_coordinate`, `::test_non_palindrome_reverse_strand_asymmetric` |
| #3 `_rc()` handles full IUPAC | `test_dna_sanity.py` | `TestReverseComplement::test_rc_handles_each_iupac_code`, `::test_rc_is_involutive` |
| #4 IUPAC regex patterns are cached | `test_dna_sanity.py`, `test_performance.py` | `TestIUPACPattern::test_pattern_cache_*`, `TestIUPACPatternCachePerformance::test_warm_cache_is_near_free` |
| #5 Circular wrap-around midpoint formula | `test_circular_math.py` | `TestFeatureMidpoint::test_wrap_around_*`, `::test_non_wrapped_vs_wrapped_disagree_with_naive_formula` |
| #6 Circular wrap-around RE scanning | `test_dna_sanity.py` | `TestRestrictionScan::test_circular_wraparound_ecori_found`, `::test_circular_wraparound_not_found_when_linear`, `::test_circular_wraparound_recut_position`, `::test_circular_wraparound_unique_filter` |
| #7 Data-file saves always back up | `test_data_safety.py` | `TestSafeSaveJson::test_*`, `TestSafeLoadJson::test_*`, `TestRealFilesNeverTouched` |

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

Every release after v0.1.0 is cut via `./release.sh X.Y.Z`. Version strings live
in two places — `pyproject.toml` `[project] version` and `splicecraft.py`
`__version__` — and release.sh updates both via sed.

### v0.2.2 (2026-04-12) — latest on PyPI
Bug-fix release shipping six fixes from the 2026-04-12 review pass:
circular RE scan wrap-around; Golden Braid primer length validation;
pLannotate stale-record guard; undo-stack clears on plasmid load;
widened shrink guard; linear-click consistency with `_bp_in`. See
"Recent fixes" section above for details.

### v0.2.1 (2026-04-12)
Golden Braid assembly UI (`PartsBinModal`, `DomesticatorModal`,
`ConstructorModal`); full Primer3-backed primer design workbench
(`PrimerDesignScreen` with Detection / Cloning / Golden Braid / Generic
modes); per-platform user data dir via `platformdirs`; atomic saves
with `.bak`; Textual 8.2.3 upgrade (fixes Footer reactive bug);
rename + focus-aware delete for library entries.

### v0.2.0 (2026-04-10)
First PyPI release. Feature deletion (Delete key); linear map view
toggle (`v`); strand-aware DNA layout; braille feature bars + single-bp
triangles; label-above / label-below layout; connector lines (`l`);
full NEB enzyme catalog (~204 enzymes, Type IIS); inside tick marks;
full-length feature labels with proximity placement; default library
entry (pACYC184).

### v0.1.0 (2026-03-23)
Initial prototype. Braille circular map, NCBI fetch, local .gb loading,
library, feature sidebar, sequence panel, undo/redo, basic restriction
site overlay.

### Stubs still in menus (not implemented)
- **Features > Add Feature** (`action_add_feature`) — `coming soon` only.
- **Build > Simulate Assembly** — `coming soon`.
- **Build > New Part editor** — `coming soon`.

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
4. **If the worker captures state that can change under it** (e.g. `self._current_record`), guard the callback with `if self._current_record is captured_record` and drop the result gracefully. See `_run_plannotate_worker` as the template — stale results are discarded with a warning.
5. Use `FetchModal._do_fetch` (line ~3190) or `_run_plannotate_worker` (line ~6955) as canonical templates

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

1. **Read this file first.** It gives you architecture without reading 7,100 lines.
2. **Run `python3 -m pytest -q`** before and after any change. 388 tests, ~75 s. The biology-correctness subset (`tests/test_dna_sanity.py`) runs in < 1 s if you want a faster inner loop.
3. **Check `/tmp/splicecraft.log`** (or `$SPLICECRAFT_LOG`) when debugging runtime issues. Every session has a unique 8-char ID.
4. **Don't break the sacred invariants.** Every one of them has a test (see the mapping table above). If you're touching `_scan_restriction_sites`, `_rc`, `_iupac_pattern`, `_translate_cds`, `_bp_in`, or the feature-midpoint formula, the relevant tests will tell you immediately if you got it wrong.
5. **Follow the error handling convention**: `_log.exception` for the stack trace, `notify()` or `Static.update("[red]...[/]")` for the user. Never let raw tracebacks hit the TUI.
6. **When in doubt about real-world behavior** — eyeball it on pUC19 (`L09137`) and pACYC184 (`MW463917.1`), both of which are fetched at first-run.
7. **Sister project for reference:** `/home/seb/proteoscope/scriptoscope.py` is the same author's larger app and the source of most patterns here (including the test suite layout). Read its `CLAUDE.md` for cross-pollination ideas.
