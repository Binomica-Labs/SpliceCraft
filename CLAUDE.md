# CLAUDE.md — AI Agent Context for SpliceCraft

This file is the **agent handoff document** for SpliceCraft. Any AI agent can read this file to understand the architecture, conventions, and design decisions behind the codebase — and pick up development without needing the full conversation history.

The project is developed in continuous collaboration between a human bioinformatician and an AI agent (Claude Opus 4.6).

---

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer and sequence editor** built with Python 3.12+ / Textual / Biopython. Renders Unicode braille-dot circular and linear plasmid maps directly in the terminal, with a per-base sequence panel, restriction-site overlays, library, and Golden Braid L0 assembly tooling.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs org, user ATinyGreenCell)

- **Single-file architecture:** the entire app is `splicecraft.py` (~7,100 lines). Intentional — avoids import complexity and keeps the codebase greppable. (Sibling project ScriptoScope follows the same convention at ~8,600 lines.)
- **Test suite:** 392 tests across 10 files in `tests/` (last refresh 2026-04-13). Full run ~75 s, biology subset (`test_dna_sanity.py`) < 1 s.
- **Dependencies:** `textual>=8.2.3`, `biopython>=1.87`, `primer3-py>=2.3.0`, `platformdirs>=4.2`, plus `pytest>=9.0` / `pytest-asyncio>=1.3` for tests. Users install via `pipx install splicecraft` or `pip install splicecraft` inside a venv. **Optional runtime:** `pLannotate` (conda, GPL-3) for the Shift+A annotation feature.
- **Published on PyPI** as `splicecraft`. Releases cut via `./release.sh X.Y.Z` (bumps version in both `pyproject.toml` and `splicecraft.py`, runs tests, builds, commits+tags+pushes; GitHub Actions `publish.yml` then publishes via Trusted Publishing / OIDC). Latest published: **v0.3.0** (per recent git log).

## How to run

```bash
cd ~/SpliceCraft
python3 splicecraft.py              # empty canvas
python3 splicecraft.py L09137       # fetch pUC19 from NCBI
python3 splicecraft.py myplasmid.gb # open local GenBank file
python3 -m pytest -q                # full test suite

# End users:
pipx install splicecraft
splicecraft
```

Logs: `/tmp/splicecraft.log` (override with `$SPLICECRAFT_LOG`). Each line is prefixed with an 8-char session ID for multi-run grepping.

### Optional: pLannotate for automatic annotation

Press **Shift+A** (or click ◈ in the library panel) to run pLannotate on the current plasmid. SpliceCraft only calls it as a subprocess — it is never imported (pLannotate is GPL-3; subprocess boundary avoids license entanglement).

```bash
conda create -n plannotate -c conda-forge -c bioconda plannotate
conda activate plannotate
plannotate setupdb          # downloads ~500 MB of BLAST/diamond DBs
```

If pLannotate is not on `PATH`, Shift+A just notifies the user and returns — nothing crashes.

## Architecture (single file: `splicecraft.py`)

### Top-level structure (line numbers ±30)

| Lines | Section |
|-------|---------|
| 1–100 | Docstring, imports, user data dir (`platformdirs`), legacy migration |
| 100–330 | Dependency check, rotating logger, atomic JSON persistence (`_safe_save_json` / `_safe_load_json` — tempfile + `.bak` + shrink guard) |
| 331–760 | NEB enzyme catalog (~204), IUPAC tables, cached regex, IUPAC-aware `_rc`, `_scan_restriction_sites()` (both strands, palindrome-aware, **circular wrap-around**) |
| 760–1055 | Sequence panel rendering helpers (`_assign_chunk_features`, `_render_feature_row_pair`, `_build_seq_inputs` (memoized), `_build_seq_text`) |
| 1055–1550 | OSC52 clipboard, CDS translation, GenBank I/O, **pLannotate integration** (subprocess-only) |
| 1550–2350 | `_Canvas` + `_BrailleCanvas` (sub-character braille resolution); `PlasmidMap` widget |
| 2350–2620 | `FeatureSidebar`, `LibraryPanel` |
| 2620–3375 | `SequencePanel` (DNA viewer, click-to-cursor, drag selection); modal dialogs |
| 3375–3830 | `MenuBar`; primer design functions (`_design_gb/cloning/detection/generic_primers`, `_pick_binding_region`) |
| 3830–4535 | Golden Braid L0 UI (`PartsBinModal`, `DomesticatorModal`, `ConstructorModal`) |
| 4535–5660 | `PrimerDesignScreen`; small modals (Quit, Picker, Rename, DeleteConfirm) |
| 5660–7040 | `PlasmidApp` — main controller, keybindings, undo/redo, `@work` threads (NCBI, seed, pLannotate) |
| 7040–end | `main()` entry point |

### Key design patterns

- **All rendering uses Rich `Text`** — no curses
- **Braille canvas** gives sub-character pixel resolution (2x4 dots per terminal cell)
- **Feature coordination:** map click → sidebar highlight → sequence scroll (and vice versa via Textual messages)
- **Undo/redo:** snapshot-based (full seq + cursor + deepcopy of SeqRecord), max 50
- **Restriction sites:** scanned on load/edit, stored as `resite` (recognition bar) + `recut` (cut marker) dicts
- **Caching:** `PlasmidMap`, `SequencePanel`, and IUPAC regex patterns all cache rendered/compiled output keyed on state. Cache keys include `id(self._feats)` since lists are reassigned (not mutated) on load
- **Workers:** `@work(thread=True)` for NCBI fetch and first-run library seed. Both use `call_from_thread` to push results back to the UI

## Logging convention

```python
_log = logging.getLogger("splicecraft")
# Rotating file at /tmp/splicecraft.log, 2MB × 2 backups
# Every line prefixed with [session_id] for multi-run grepping
```

- **User-facing errors** → `self.notify(...)` or `Static.update("[red]...[/]")`. Never raw tracebacks.
- **Diagnostic detail** → `_log.exception("context: %s", ...)` inside `except` blocks. Stack traces go to the log file only.
- **Worker errors** → log with `_log.exception`, then push a friendly message to the UI via `call_from_thread`.

## Sacred invariants (DO NOT BREAK)

Every invariant below has at least one test protecting it. See the **Sacred invariant → test mapping** section below.

1. **Palindromic enzymes are scanned forward only.** `_scan_restriction_sites` must skip the reverse scan for palindromic sites and add only a bottom-strand `recut`. Scanning both strands for palindromes double-counts every site.

2. **Reverse-strand resite positions use the forward coordinate.** A reverse-strand hit at position `p` (after RC) is stored as `p`, not `n - p - site_len`. The cut maps via `site_len - 1 - fwd_cut`.

3. **`_rc()` handles full IUPAC.** Reverse-complement must translate ambiguity codes (R, Y, W, S, M, K, B, D, H, V, N) via `_IUPAC_COMP`, not just ACGT.

4. **IUPAC regex patterns are cached.** `_iupac_pattern()` uses `_PATTERN_CACHE` to avoid recompiling ~200 patterns on every restriction scan.

5. **Circular wrap-around midpoints.** When computing the midpoint of a feature for label placement, use `arc_len = (end - start) % total` then `(start + arc_len // 2) % total`. The naive `(start + (end - start) // 2) % total` puts the label opposite the actual arc for wrapped features.

6. **Circular wrap-around restriction scanning.** `_scan_restriction_sites(circular=True)` (default) scans `seq + seq[:max_site_len-1]` so recognition sequences spanning the origin are found. Each wrap-around hit is emitted as **two resite pieces** (labeled tail `[p, n)` + unlabeled head `[0, (p+site_len) - n)`) and **one recut** at `(p + fwd_cut) % n`. Downstream code that counts resites for filtering must count only labeled pieces.

7. **Data-file saves always back up.** `_safe_save_json` writes a `.bak` of the existing file before replacing it, via `tempfile.mkstemp` + `os.fsync` + `os.replace`. Shrink guard logs a warning if writing fewer entries than exist. Never bypass `_safe_save_json` — it is the user's only recovery path.

8. **Wrap-aware feature length everywhere.** Use `_feat_len(start, end, total)` — returns `(total - start) + end` when `end < start`, else `end - start`. All sort keys, length displays, and biological-length checks must route through it. Naive `end - start` gives negative values for wrap features and breaks z-order, primer design, and sidebar displays.

9. **Wrap-feature integrity in record edits.** `int(CompoundLocation.start)` returns `min(parts.start)` and `int(.end)` returns `max(parts.end)`, silently flattening wrap features into whole-plasmid FeatureLocations. `_rebuild_record_with_edit` must per-part shift wrap features and only collapse to FeatureLocation when 1 part survives. Zero-width post-edit features must be dropped (no 1-bp ghost stubs).

10. **Undo snapshots must be deepcopied.** `_push_undo`, `_action_undo`, `_action_redo` all `deepcopy(self._current_record)` so future in-place mutations can't poison the stack.

## pLannotate integration

Shift+A (or ◈ in the library panel, or `Features > Annotate with pLannotate`) runs pLannotate as a subprocess and merges results into the current record.

### Design principles

1. **Subprocess only, never import.** pLannotate is GPL-3 — importing would arguably create a combined work under GPL. **Never `import plannotate`.**
2. **Optional runtime dependency.** SpliceCraft works without it. UI shows install hint when missing.
3. **Size cap preflighted** at 50 kb (matches pLannotate's `MAX_PLAS_SIZE`).
4. **Merge, don't replace.** Existing features preserved; pLannotate hits appended with `note="pLannotate"` qualifier. Hits matching `(type, start, end, strand)` of an existing feature are skipped.
5. **Background worker** with stale-record guard: callback checks `self._current_record is captured_record` and discards stale results.
6. **Re-entry guard** via `_plannotate_running` flag (with `finally` cleanup).
7. **Undo-able.** Worker calls `_push_undo()` before applying merged record.
8. **Dirty flag.** Marks both `lib.set_dirty(True)` and `self._unsaved=True` via `_mark_dirty()`.

### Code layout

| Function / class | Purpose |
|---|---|
| `PlannotateError` (+ 4 subclasses) | Exception hierarchy with `user_msg` / `detail` attrs |
| `_PLANNOTATE_MAX_BP = 50_000` | Matches pLannotate's `MAX_PLAS_SIZE` |
| `_plannotate_status()` | `shutil.which`-based probe, cached |
| `_run_plannotate(record, timeout=180)` | Subprocess runner, raises `PlannotateError` subclasses |
| `_merge_plannotate_features(orig, annotated)` | Pure function: returns new SeqRecord |
| `PlasmidApp.action_annotate_plasmid` | Shift+A action: preflights + kicks off worker |
| `PlasmidApp._run_plannotate_worker` | `@work(thread=True)` — subprocess + merge + UI update |

Failure modes (`PlannotateNotInstalled`, `PlannotateMissingDb`, `PlannotateTooLarge`, `PlannotateFailed`) map to actionable user notifications. Full traceback always written to `/tmp/splicecraft.log`.

## Test suite

Originally added 2026-04-11 to protect the sacred invariants; expanded each session. Full suite runs in ~75 s; biology-correctness subset runs in < 1 s.

### Running

```bash
python3 -m pytest -q                        # all 392 tests
python3 -m pytest tests/test_dna_sanity.py  # only biology (< 1 s)
python3 -m pytest -k "palindrome"           # filter by name
python3 -m pytest -x                        # stop on first failure
```

`pyproject.toml` sets `asyncio_mode = "auto"` so async tests don't need `@pytest.mark.asyncio`. `tests/conftest.py` defines `tiny_record` / `tiny_gb_path` / `isolated_library` fixtures, and installs the **autouse** `_protect_user_data` fixture that monkeypatches `_LIBRARY_FILE`, `_PARTS_BIN_FILE`, `_PRIMERS_FILE` and their caches to tmp paths. **No test can write to real user files.**

### Files

| File | Tests | Covers |
|------|------:|--------|
| `test_dna_sanity.py` | 74 | Sacred invariants 1-6 (RE scan, `_rc`, `_iupac_pattern`, codon table, circular wrap scan); Type IIS cut-outside-recognition; `_translate_cds` forward & reverse |
| `test_circular_math.py` | 38 | Sacred invariant #5 (wrap midpoint); `_bp_in` for wrapped/non-wrapped/zero-width; `_feat_len` |
| `test_edit_record.py` | 14 | Sacred invariant #9: wrap features survive insert/replace as CompoundLocation, collapse to FeatureLocation when 1 part remains, fully-consumed features dropped (no 1-bp stubs) |
| `test_genbank_io.py` | 59 | `load_genbank` round-trip; `_save_library` / `_load_library` JSON round-trip + corruption recovery |
| `test_data_safety.py` | 28 | Sacred invariant #7 (atomic saves, `.bak` recovery); `_protect_user_data` fixture confirmation |
| `test_primers.py` | 60 | Detection/cloning/Golden Braid/generic; **wrap-region primer design** (`_slice_circular`, Primer3 template rotation, modular position mapping) |
| `test_domesticator.py` | 41 | Golden Braid L0 positions/overhangs, part validation, assembly lanes |
| `test_plannotate.py` | 24 | Availability detection, size-cap preflight, feature merging, subprocess error paths (subprocess never actually invoked) |
| `test_smoke.py` | 45 | Textual app mounts; panels present; rotation/view-toggle/RE-toggle; pLannotate UI + re-entry guard; `_apply_record` semantics; sidebar wrap-coord display; undo snapshot independence |
| `test_performance.py` | 9 | Budget enforcement (loose, 4-20× headroom): scan pUC19 < 30 ms, scan 10 kb < 150 ms, scaling < 8×, `_iupac_pattern` warm < 5 ms, `_rc(10 kb)` < 2 ms, `_build_seq_text(20 kb)` < 200 ms, `_BUILD_SEQ_CACHE` populated after first call |

### Sacred invariant → test mapping

| Invariant | Test file | Test method |
|---|---|---|
| #1 Palindrome forward only | `test_dna_sanity.py` | `TestRestrictionScan::test_ecori_single_site_not_double_counted`, `::test_palindromes_produce_one_recut_per_site` |
| #2 Reverse-strand forward coord | `test_dna_sanity.py` | `TestRestrictionScan::test_non_palindrome_on_reverse_strand_uses_forward_coordinate` |
| #3 `_rc()` IUPAC | `test_dna_sanity.py` | `TestReverseComplement::test_rc_handles_each_iupac_code`, `::test_rc_is_involutive` |
| #4 Regex cache | `test_dna_sanity.py`, `test_performance.py` | `TestIUPACPattern::test_pattern_cache_*`, `TestIUPACPatternCachePerformance::test_warm_cache_is_near_free` |
| #5 Wrap midpoint | `test_circular_math.py` | `TestFeatureMidpoint::test_wrap_around_*` |
| #6 Circular wrap RE scan | `test_dna_sanity.py` | `TestRestrictionScan::test_circular_wraparound_*` |
| #7 Atomic saves | `test_data_safety.py` | `TestSafeSaveJson::*`, `TestSafeLoadJson::*`, `TestRealFilesNeverTouched` |
| #8 `_feat_len` | `test_circular_math.py` | `TestFeatLen::*` |
| #9 Wrap edit integrity | `test_edit_record.py` | (whole file) |
| #10 Undo deepcopy | `test_smoke.py` | `test_undo_snapshot_independence_under_in_place_mutation` |

### Conventions

- **Cross-validate against Biopython** where possible (codon table, reverse-complement). If Biopython's standard table changes, the test fails noisily.
- **Hand-verifiable** test inputs — every restriction-site test uses a sequence short enough to count expected hits by eye.
- **Regression guards cite the date** — every test protecting a past bug has a docstring like `# Regression guard for 2026-03-30 fix`.
- **No network, no real files** — all tests use synthetic `SeqRecord`s and monkeypatched paths.
- **Performance budgets are LOOSE** (6-20× headroom). They catch architectural regressions, not micro-perf drift.

### Adding a new test

1. Pick the right file (or add a new one).
2. For SeqRecord-based tests, use `tiny_record` fixture.
3. For Textual async tests: `async def test_*` (no decorator), `async with app.run_test(size=TERMINAL_SIZE) as pilot: await pilot.pause(); await pilot.pause(0.5)`. Double-pause is needed for `call_after_refresh` callbacks.
4. For perf tests, warm the cache then average 10-20 iterations.

## Performance notes

After the 2026-04-11 optimization pass, key wins:

1. **Sidebar populate cascade suppressed** via `_populating` flag + `call_after_refresh` deferred reset — eliminates duplicate `_build_seq_text` per record load. Saves ~50-180 ms depending on plasmid size.
2. **Memoized styles + sorted annotations** in `_build_seq_inputs()` cached in module-level `_BUILD_SEQ_CACHE` (4-entry, identity-keyed). Cursor moves don't recompute. ~40% savings warm.
3. **Per-chunk `str.translate`** for reverse strand instead of per-base. ~60% savings on `_build_seq_text` cold-call.
4. **`_SCAN_CATALOG`** precomputed at import time eliminates per-scan `_rc`/`_iupac_pattern`/`len` calls. ~15% savings on small plasmids.

What was profiled but **not touched**: Textual compositor (framework), Rich `Text.append` (already efficient), `PlasmidMap._draw` (canvas-bound, well-cached), import time (Textual + Rich dominate).

## Released vs. unreleased state

Versions live in `pyproject.toml` and `splicecraft.py` `__version__`; `release.sh` updates both via sed. See `git log --oneline` for full release history. Recent: v0.3.0 (Mutagenize modal with codon registry/harmonization), v0.2.8 (deep-copy record in undo/redo snapshots).

### Stubs still in menus (not implemented)
- **Features > Add Feature** — `coming soon`
- **Build > Simulate Assembly** — `coming soon`
- **Build > New Part editor** — `coming soon`

## Patterns worth porting from ScriptoScope (`/home/seb/proteoscope/scriptoscope.py`)

ScriptoScope (~8,600 lines) is the more mature sibling. Pre-validated patterns to lift if SpliceCraft grows:

| Pattern | When SpliceCraft needs it |
|---------|---------------------------|
| Thread-local `Console` for `_text_to_content` | If sequence-panel render starts blowing 33 ms/frame budget. Required because shared Console + worker threads = lock contention |
| Two-level render cache (`_seq_render_cache` + `_content_cache`) | If repainting on cursor moves becomes janky. LRU via `OrderedDict.move_to_end` |
| `Select.NULL`/`BLANK` sentinel filtering | Only if SpliceCraft adds a `Select` widget (currently none) |
| `@lru_cache(1)` availability checks for optional CLI tools | If SpliceCraft ever shells out to BLAST / Prodigal beyond pLannotate |

## Known pitfalls

1. **Bare `except` is forbidden.** All file-touching helpers must `_log.exception(...)` instead of `pass`.
2. **Wrapped features (`end < start`) are first-class citizens.** Anywhere you compute distances, midpoints, or "is bp inside this feature", use the modular form via `_bp_in()` or `_feat_len()`. See sacred invariants #5, #6, #8, #9.
3. **Cache keys use `id(...)` of feature lists.** Correct *only* because the app reassigns lists on load rather than mutating them in-place. If you start mutating `self._feats` in-place, caches return stale renders.
4. **Textual reactive auto-invalidation depends on field assignment, not mutation.** `self._feats = new_list` triggers refresh; `self._feats.append(x)` does not.
5. **Single-file means giant diffs are normal.** When a refactor touches the rendering layer, expect 100+ line edits. The greppability tradeoff is worth it.
6. **Primer3 is linear-only.** For wrap regions, rotate template to `seq[start:] + seq[:start]` before calling, then unrotate positions via `(coord + rotation) % total`. See `_design_detection_primers`.
7. **`_source_path` is preserved through in-place edits.** Only cleared when `clear_undo=True` (fresh loads). Otherwise Ctrl+S after pLannotate or primer-add would forget the original file.

## How to extend

### Adding a modal dialog
1. Subclass `ModalScreen[ReturnType]` (see `FetchModal`, `OpenFileModal`, `PartsBinModal`)
2. Implement `compose()` with the form layout
3. Call `self.dismiss(result)` to return a value
4. Push with `app.push_screen(MyModal(), callback=on_result)`

### Adding a heavy operation
1. Decorate with `@work(thread=True)`
2. Wrap body in `try` / `except Exception` and call `_log.exception(...)` in handler
3. Push results to UI via `self.app.call_from_thread(callback)` — never touch widgets directly from worker thread
4. **If the worker captures state that can change under it** (e.g. `self._current_record`), guard the callback with `if self._current_record is captured_record` and drop stale results gracefully. Template: `_run_plannotate_worker`.

### Adding a menu action
1. Add `action_my_thing(self)` on `PlasmidApp`
2. Add binding in `BINDINGS`
3. Add menu item in `MenuBar.compose()`

## Future work (user is undecided)

The user is on the fence between:
- **Merging** SpliceCraft, ScriptoScope, MitoShift, RefHunter, molCalc into one Textual app with multiple "modes"
- **Keeping them separate** as focused single-purpose apps and (optionally) extracting shared utilities into pure-Python modules

Either direction is viable. The single-file convention and shared logging/error patterns documented here keep the merge option open without forcing it.

## For future agents

1. **Read this file first.** It gives you architecture without reading 7,100 lines.
2. **Run `python3 -m pytest -q`** before and after any change. 392 tests, ~75 s. Biology subset (`tests/test_dna_sanity.py`) runs in < 1 s for a faster inner loop.
3. **Check `/tmp/splicecraft.log`** (or `$SPLICECRAFT_LOG`) when debugging. Every session has a unique 8-char ID.
4. **Don't break the sacred invariants.** Each has a test (see mapping table). If you touch `_scan_restriction_sites`, `_rc`, `_iupac_pattern`, `_translate_cds`, `_bp_in`, `_feat_len`, the midpoint formula, or `_rebuild_record_with_edit`, the relevant tests will tell you immediately.
5. **Follow the error handling convention**: `_log.exception` for stack traces, `notify()` or `Static.update("[red]...[/]")` for the user. Never let raw tracebacks hit the TUI.
6. **When in doubt about real-world behavior** — eyeball it on pUC19 (`L09137`) and pACYC184 (`MW463917.1`), both fetched at first-run.
7. **Sister project for reference:** `/home/seb/proteoscope/scriptoscope.py` is the same author's larger app and source of most patterns here.
8. **Past fix history lives in git.** Use `git log --oneline` and `git show <hash>` rather than restoring fix-log sections to this file.
