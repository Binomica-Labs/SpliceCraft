# CLAUDE.md — AI Agent Context for SpliceCraft

This file is the **agent handoff document** for SpliceCraft. Any AI agent can read this file to pick up development without needing the full conversation history.

The project is developed in continuous collaboration between a human bioinformatician and an AI agent (Claude Opus 4.6+).

---

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer, sequence editor, and cloning/mutagenesis workbench** built with Python 3.10+ / Textual / Biopython. Renders Unicode braille-dot plasmid maps directly in the terminal, with a per-base sequence panel, restriction-site overlays, a plasmid library, Golden Braid L0 assembly tooling, Primer3-backed primer design, and SOE-PCR site-directed mutagenesis.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs org, user ATinyGreenCell)

- **Single-file architecture:** the entire app is `splicecraft.py` (~10,620 lines). Intentional — avoids import complexity and keeps the codebase greppable. Sibling project ScriptoScope follows the same convention at ~8,600 lines.
- **Test suite:** 593 tests across 15 files in `tests/` (last refresh 2026-04-20). Full run ~205 s, biology subset (`test_dna_sanity.py`) < 1 s. `test_invariants_hypothesis.py` adds property-based fuzzing on top of hand-written regression tests.
- **Dependencies:** `textual>=8.2.3`, `biopython>=1.87`, `primer3-py>=2.3.0`, `platformdirs>=4.2`, plus `pytest>=9.0` / `pytest-asyncio>=1.3` / `hypothesis>=6.100` for tests. Users install via `pipx install splicecraft`. **Optional runtime:** `pLannotate` (conda, GPL-3) for the Shift+A annotation feature.
- **Published on PyPI** as `splicecraft`. Releases cut via `./release.sh X.Y.Z` (bumps version in both `pyproject.toml` and `splicecraft.py`, runs tests, builds, commits+tags+pushes; GitHub Actions `publish.yml` then publishes via Trusted Publishing / OIDC). Latest published: **v0.3.1**.

## How to run

```bash
cd ~/SpliceCraft
python3 splicecraft.py              # empty canvas
python3 splicecraft.py L09137       # fetch pUC19 from NCBI
python3 splicecraft.py myplasmid.gb # open local GenBank file (.gb/.gbk/.dna)
python3 -m pytest -q                # full test suite

# End users:
pipx install splicecraft
splicecraft
```

Logs: `~/.local/share/splicecraft/logs/splicecraft.log` (override with `$SPLICECRAFT_LOG`). Each line is prefixed with an 8-char session ID for multi-run grepping.

### Optional: pLannotate for automatic annotation

Press **Shift+A** (or click ◈ in the library panel) to run pLannotate on the current plasmid. SpliceCraft only calls it as a subprocess — it is never imported (pLannotate is GPL-3; subprocess boundary avoids license entanglement).

```bash
conda create -n plannotate -c conda-forge -c bioconda plannotate
conda activate plannotate
plannotate setupdb          # downloads ~500 MB of BLAST/diamond DBs
```

If pLannotate is not on `PATH`, Shift+A notifies the user and returns — nothing crashes.

## Architecture (single file: `splicecraft.py`)

### Top-level structure (line numbers ±30, current 2026-04-20)

| Lines | Section |
|-------|---------|
| 1–200 | Docstring, imports, user data dir (`platformdirs`), legacy migration, dependency check, rotating session-tagged logger (log file in `_DATA_DIR/logs`), feature-colour palette |
| 201–385 | Atomic JSON persistence (`_safe_save_json` / `_safe_load_json` + `_extract_entries` — schema-envelope format `{"_schema_version": 1, "entries": [...]}` with legacy bare-list back-compat; tempfile + `os.replace` + `.bak` + shrink guard) |
| 386–408 | Library cache loaders (`_load_library` / `_save_library`) |
| 409–1448 | NEB enzyme catalog (~204), IUPAC tables + cached regex, `_rc`, `_scan_restriction_sites` (palindrome-aware, wrap-around), `_assign_chunk_features`, `_render_feature_row_pair`, memoized `_build_seq_inputs` and `_build_seq_text`, OSC-52 clipboard, `_translate_cds` |
| 1449–1521 | Char-aspect detection + label helpers |
| 1522–1659 | GenBank I/O (`fetch_genbank`, `load_genbank` auto-detecting `.gb`/`.dna`, `_record_to_gb_text`, `_gb_text_to_record`) |
| 1660–1875 | **pLannotate** subprocess integration (`PlannotateError` hierarchy, `_run_plannotate`, `_merge_plannotate_features`) |
| 1876–1985 | `_Canvas` + `_BrailleCanvas` (sub-cell braille resolution) |
| 1986–2753 | `PlasmidMap` widget — circular/linear draw, label placement, `_draw_cache` |
| 2754–2868 | `FeatureSidebar` — scrollable feature table with click-to-select |
| 2869–3036 | `LibraryPanel` — plasmid library list, rename/delete buttons |
| 3037–3485 | `SequencePanel` — DNA viewer, click-to-cursor, drag selection |
| 3486–3825 | Core modals (`EditSeqDialog`, `FetchModal` with in-flight staleness guard, `OpenFileModal`, `DropdownScreen`) |
| 3826–3867 | `MenuBar` widget |
| 3868–4076 | Golden Braid L0 position catalog (BsaI overhangs, position constraints) |
| 4077–4130 | Parts-bin + primer-library persistence |
| 4131–4925 | Codon-usage registry (`_codon_*`), Kazusa parser, NCBI taxid search (`_safe_xml_parse`), harmonization, CAI/GC. Crash-recovery config (`_CRASH_RECOVERY_DIR`) sits at the top of this slab |
| 4926–5437 | SOE-PCR site-directed mutagenesis primer design (`_mut_*`) |
| 5438–5598 | `PartsBinModal` |
| 5599–6147 | `DomesticatorModal` + `ConstructorModal` (Golden Braid L0 UI) |
| 6148–6485 | `NcbiTaxonPickerModal` + `SpeciesPickerModal` (codon-table picker) |
| 6486–6716 | Mutagenize helpers (`_MutPreview`, `AminoAcidPickerModal`) |
| 6717–7361 | `MutagenizeModal` — full mutagenesis workflow |
| 7362–8460 | `PrimerDesignScreen` — full-screen primer workbench |
| 8461–8681 | Small modals (`UnsavedQuitModal`, `PlasmidPickerModal`, `RenamePlasmidModal`, `LibraryDeleteConfirmModal`) |
| 8682–10551 | `PlasmidApp` — main controller, keybindings, per-plasmid undo/redo stashes, crash-recovery autosave, `@work` threads |
| 10552–end | `main()` entry point |

### Key design patterns

- **Rich `Text` for all rendering** — no curses.
- **Braille canvas** gives sub-character pixel resolution (2×4 dots per terminal cell).
- **Feature coordination:** map click → sidebar highlight → sequence scroll (and back via Textual messages).
- **Undo/redo:** snapshot-based (full seq + cursor + `deepcopy` of SeqRecord), max 50. **Per-plasmid stashes**: switching plasmids stashes the outgoing history under the old `record.id` and restores the incoming plasmid's history (LRU-capped at 10 plasmids). Ctrl+Z never yanks you to an unrelated edit.
- **Crash-recovery autosave:** every dirty edit debounces (3 s) a write of the current record to `_DATA_DIR/crash_recovery/{safe_id}.gb`. Cleared on successful save or explicit abandon. On startup a non-empty dir notifies the user so they can recover.
- **Restriction sites:** scanned on load/edit, stored as `resite` (recognition bar) + `recut` (cut marker) dicts.
- **Caching:** `PlasmidMap._draw_cache`, `_BUILD_SEQ_CACHE`, `_PATTERN_CACHE`, `_SCAN_CATALOG` — all keyed on inputs (including `id(self._feats)` since lists are reassigned, not mutated, on load).
- **Workers:** `@work(thread=True)` for NCBI fetch, library seed, pLannotate, Kazusa codon fetch. Results pushed back via `call_from_thread`, with stale-record guards where the worker captures `self._current_record`.

## Logging convention

```python
_log = logging.getLogger("splicecraft")
# Rotating file at _DATA_DIR/logs/splicecraft.log (platform-specific), 2MB × 2 backups
# Every line prefixed with [session_id] for multi-run grepping
```

- **User-facing errors** → `self.notify(...)` or `Static.update("[red]...[/]")`. Never raw tracebacks.
- **Diagnostic detail** → `_log.exception("context: %s", ...)` inside `except` blocks. Stack traces go to the log file only.
- **Worker errors** → log with `_log.exception`, then push a friendly message to the UI via `call_from_thread`.
- **Narrow exception types.** Use `except NoMatches:` around `query_one` lookups, `except ET.ParseError:` around XML, `except (OSError, json.JSONDecodeError):` around file I/O. Reserve bare `except Exception` for worker bodies where anything can happen — and always log there.

## Sacred invariants (DO NOT BREAK)

Every invariant below has at least one test protecting it. See the **Sacred invariant → test mapping** section below.

1. **Palindromic enzymes are scanned forward only.** `_scan_restriction_sites` must skip the reverse scan for palindromic sites and add only a bottom-strand `recut`. Scanning both strands for palindromes double-counts every site.

2. **Reverse-strand resite positions use the forward coordinate.** A reverse-strand hit at position `p` (after RC) is stored as `p`, not `n - p - site_len`. The cut maps via `site_len - 1 - fwd_cut`.

3. **`_rc()` handles full IUPAC.** Reverse-complement must translate ambiguity codes (R, Y, W, S, M, K, B, D, H, V, N) via `_IUPAC_COMP`, not just ACGT.

4. **IUPAC regex patterns are cached.** `_iupac_pattern()` uses `_PATTERN_CACHE` to avoid recompiling ~200 patterns on every restriction scan.

5. **Circular wrap-around midpoints.** When computing the midpoint of a feature for label placement, use `arc_len = (end - start) % total` then `(start + arc_len // 2) % total`. The naive `(start + (end - start) // 2) % total` puts the label opposite the actual arc for wrapped features.

6. **Circular wrap-around restriction scanning.** `_scan_restriction_sites(circular=True)` (default) scans `seq + seq[:max_site_len-1]` so recognition sequences spanning the origin are found. Each wrap-around hit is emitted as **two resite pieces** (labeled tail `[p, n)` + unlabeled head `[0, (p+site_len) - n)`) and **one recut** at `(p + fwd_cut) % n`. Downstream code that counts resites for filtering must count only labeled pieces.

7. **Data-file saves always back up.** `_safe_save_json` writes a `.bak` of the existing file before replacing it, via `tempfile.mkstemp` + `os.fsync` + `os.replace`. Shrink guard logs a warning if writing fewer entries than exist. Writes envelope format `{"_schema_version": 1, "entries": [...]}` — loaders accept both envelope and legacy bare-list (pre-0.3.1) via `_extract_entries`, so upgrades never lose data. Future-version writes warn but still load. Never bypass `_safe_save_json` — it is the user's only recovery path.

8. **Wrap-aware feature length everywhere.** Use `_feat_len(start, end, total)` — returns `(total - start) + end` when `end < start`, else `end - start`. All sort keys, length displays, and biological-length checks must route through it. Naive `end - start` gives negative values for wrap features and breaks z-order, primer design, and sidebar displays.

9. **Wrap-feature integrity in record edits.** `int(CompoundLocation.start)` returns `min(parts.start)` and `int(.end)` returns `max(parts.end)`, silently flattening wrap features into whole-plasmid FeatureLocations. `_rebuild_record_with_edit` must per-part shift wrap features and only collapse to FeatureLocation when 1 part survives. Zero-width post-edit features must be dropped (no 1-bp ghost stubs).

10. **Undo snapshots must be deepcopied.** `_push_undo`, `_action_undo`, `_action_redo` all `deepcopy(self._current_record)` so future in-place mutations can't poison the stack.

## Core helper catalog

These are the load-bearing pure functions other code depends on. Most are at module level, a few are methods. Read these first before touching rendering, primer design, or the record pipeline.

| Helper | Line | Purpose |
|---|---:|---|
| `_safe_save_json` / `_safe_load_json` / `_extract_entries` | 251 / 331 / 228 | Atomic JSON I/O with `.bak` recovery and schema-envelope format. All four libraries go through these. |
| `_iupac_pattern` | 680 | IUPAC→regex compiler, cached in `_PATTERN_CACHE`. |
| `_IUPAC_COMP`, `_DNA_COMP_PRESERVE_CASE` | ~690 | Module-level `str.maketrans` tables (hot-path complement). |
| `_rc` | 697 | IUPAC-aware reverse complement. |
| `_feat_len`, `_slice_circular`, `_bp_in` | 701 / 707 / — | Wrap-aware geometry. Any "is bp X in feature?" or "how long is this feature" uses these. |
| `_scan_restriction_sites` | 749 | Palindrome-aware, wrap-aware restriction scan. Returns `(resites, recuts)` lists. |
| `_build_seq_inputs` / `_build_seq_text` | 1220 / 1253 | Sequence-panel renderer, memoized via `_BUILD_SEQ_CACHE`. |
| `_translate_cds` | 1423 | Forward and reverse CDS → protein. Cross-validated against Biopython. |
| `fetch_genbank` / `load_genbank` | 1545 / 1587 | NCBI Entrez fetch + local `.gb`/`.dna` load. |
| `_record_to_gb_text` / `_gb_text_to_record` | 1634 / 1654 | Serialise/deserialise SeqRecords as GenBank text. Caller's record is never mutated. |
| `_run_plannotate`, `_merge_plannotate_features` | 1719 / 1819 | pLannotate subprocess + merge. |
| `_pick_binding_region` | 3936 | Primer3-compatible region selection. |
| `_design_*_primers` | 3971+ | Detection, cloning, Golden Braid, generic primer design. |
| `_codon_*` | 4212+ | Codon-usage registry, harmonization, NCBI taxid search with `_safe_xml_parse` guard. |
| `_mut_*` | 4965+ | SOE-PCR mutagenesis primers, AA picker helpers. |
| `_rebuild_record_with_edit` | in `PlasmidApp` | Edit pipeline that preserves wrap features. Sacred invariant #9. |
| `_autosave_*` / `_stash_current_undo_and_load` | in `PlasmidApp` | Crash-recovery autosave + per-plasmid undo/redo stack stashing. |

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

Failure modes (`PlannotateNotInstalled`, `PlannotateMissingDb`, `PlannotateTooLarge`, `PlannotateFailed`) map to actionable user notifications. Full traceback always written to `~/.local/share/splicecraft/logs/splicecraft.log`.

## On-disk JSON format (schema v1)

All four persisted libraries (`library.json`, `parts_bin.json`, `primers.json`, `codon_tables.json`) use the envelope shape:

```json
{"_schema_version": 1, "entries": [...]}
```

**Legacy compatibility.** SpliceCraft < 0.3.1 wrote a bare JSON list. `_extract_entries` accepts both; a legacy file is silently rewritten as an envelope on the next save. When bumping `_CURRENT_SCHEMA_VERSION`, teach `_extract_entries` how to migrate entries forward *in the loader* so old files keep working. Files written by a newer SpliceCraft (higher version) still load but emit a warning so users know fields may drop on save.

## Crash-recovery autosave

Dirty edits trigger a 3-second debounced write of the current record to `_CRASH_RECOVERY_DIR/{safe_id}.gb` (default `~/.local/share/splicecraft/crash_recovery/`). The file is deleted on successful save (`_mark_clean`) or explicit abandon. On startup `_check_crash_recovery()` scans the dir and notifies the user if any `.gb` files survive — that means the prior session crashed before saving. The user recovers via File > Open on the named file.

Design notes:
- **`_autosave_path(record)`** sanitises `record.id` with `re.sub(r'[^A-Za-z0-9._-]', '_', ...)` and caps at 80 chars.
- **Atomic write** — `tempfile.mkstemp` in the target dir + `os.replace`, matching `_safe_save_json`'s guarantees.
- **Best-effort only** — `except Exception: _log.exception(...)` so a write failure never interrupts the user. Autosave is a safety net, not a source of truth.
- **Debounced via `self.set_timer`** — rapid edits coalesce into one write. `_mark_dirty` restarts the countdown; `_mark_clean` cancels it implicitly by deleting the target.

## Per-plasmid undo/redo stashes

`_apply_record(clear_undo=True)` (the "switch plasmid" path) stashes the outgoing plasmid's undo/redo stacks under its `record.id` in `_stashed_undo_stacks` / `_stashed_redo_stacks`, and restores the incoming plasmid's own history if it was edited before. LRU-capped at `_MAX_PLASMIDS_WITH_UNDO = 10` so opening dozens of plasmids can't balloon memory. The `_current_undo_key` tracks which plasmid's stack is live. `clear_undo=False` (in-place edits — pLannotate merge, primer-add) leaves the stacks intact.

## Test suite

Originally added 2026-04-11 to protect the sacred invariants; expanded each session.

### Running

```bash
python3 -m pytest -q                                # all 508 tests
python3 -m pytest tests/test_dna_sanity.py          # only biology (< 1 s)
python3 -m pytest tests/test_invariants_hypothesis.py  # property-based fuzzing
python3 -m pytest -k "palindrome"                   # filter by name
python3 -m pytest -x                                # stop on first failure
```

`pyproject.toml` sets `asyncio_mode = "auto"` so async tests don't need `@pytest.mark.asyncio`. `tests/conftest.py` defines `tiny_record` / `tiny_gb_path` / `isolated_library` fixtures, and installs the **autouse** `_protect_user_data` fixture that monkeypatches `_LIBRARY_FILE`, `_PARTS_BIN_FILE`, `_PRIMERS_FILE`, `_CODON_TABLES_FILE`, `_CRASH_RECOVERY_DIR`, and their caches to tmp paths. **No test can write to real user files.**

### Files

| File | Tests | Covers |
|------|------:|--------|
| `test_dna_sanity.py` | 74 | Sacred invariants 1–6; Type IIS cut-outside-recognition; `_translate_cds` forward & reverse |
| `test_primers.py` | 60 | Detection / cloning / Golden Braid / generic; **wrap-region primer design** (template rotation, modular position mapping) |
| `test_genbank_io.py` | 59 | `load_genbank` round-trip (GenBank + SnapGene `.dna`); `_save_library` / `_load_library` JSON round-trip + corruption recovery |
| `test_smoke.py` | 52 | Textual app mounts; panels present; rotation / view-toggle / RE-toggle; pLannotate UI + re-entry guard; `_apply_record` semantics; sidebar wrap-coord display; undo snapshot independence; **per-plasmid undo stashes + LRU eviction**; **crash-recovery autosave** |
| `test_mutagenize.py` | 49 | SOE-PCR primer design, codon substitution, `_mut_revcomp` / translate / CAI round-trips |
| `test_codon.py` | 42 | Codon registry persistence, harmonization, Kazusa parser, NCBI taxid XML safety, CAI/GC math |
| `test_domesticator.py` | 41 | Golden Braid L0 positions / overhangs, part validation, assembly lanes |
| `test_circular_math.py` | 38 | Sacred invariant #5 (wrap midpoint); `_bp_in` / `_feat_len` for wrapped / non-wrapped / zero-width |
| `test_data_safety.py` | 37 | Sacred invariant #7 (atomic saves, `.bak` recovery); **schema-envelope round-trip + legacy bare-list back-compat + future-version warning + shrink-guard counting both formats**; `features.json` redirected by `_protect_user_data`; `_protect_user_data` fixture confirmation |
| `test_add_feature.py` | 24 | **AddFeatureModal + insert pipeline**: qualifier parsing round-trip, `_extract_feature_entries_from_record` strand/wrap handling, modal form validation (empty name / invalid bases / IUPAC), save-to-library dedup, insert-at-cursor (fwd / rev / coord shift / dirty flag) |
| `test_plannotate.py` | 24 | Availability detection, size-cap preflight, feature merging, subprocess error paths (subprocess never actually invoked) |
| `test_modal_boundaries.py` | 23 | **Every modal stays inside the terminal**: root-container bounds + non-scrollable descendants fit at the baseline 160×48; AddFeatureModal-specific regression guards at 160×48 / 120×40 / 100×30 (regression for 2026-04-20 textbox-offscreen bug) |
| `test_features_library.py` | 15 | Persistent feature-library JSON round-trip, schema envelope, corruption recovery, cache invalidation, `_GENBANK_FEATURE_TYPES` curation (CDS / gene / promoter present, `source` excluded) |
| `test_edit_record.py` | 14 | Sacred invariant #9: wrap features survive insert/replace as CompoundLocation; fully-consumed features dropped (no 1-bp stubs) |
| `test_invariants_hypothesis.py` | 11 | Property-based fuzzing of sacred invariants #3, #5, #8: `_rc` involution + IUPAC closure + Biopython cross-check; `_feat_len` bounds + linear/wrap formulas; `_bp_in` count matches `_feat_len`; wrap midpoint lies on arc |
| `test_performance.py` | 9 | Budget enforcement (loose, 4–20× headroom): scan pUC19 < 30 ms, scan 10 kb < 150 ms, `_iupac_pattern` warm < 5 ms, `_rc(10 kb)` < 2 ms, `_build_seq_text(20 kb)` < 200 ms, `_BUILD_SEQ_CACHE` populated after first call |

### Sacred invariant → test mapping

| Invariant | Test file | Test method |
|---|---|---|
| #1 Palindrome forward only | `test_dna_sanity.py` | `TestRestrictionScan::test_ecori_single_site_not_double_counted`, `::test_palindromes_produce_one_recut_per_site` |
| #2 Reverse-strand forward coord | `test_dna_sanity.py` | `TestRestrictionScan::test_non_palindrome_on_reverse_strand_uses_forward_coordinate` |
| #3 `_rc()` IUPAC | `test_dna_sanity.py`, `test_invariants_hypothesis.py` | `TestReverseComplement::test_rc_handles_each_iupac_code`, `::test_rc_is_involutive`; `TestReverseComplementProperties::*` (fuzzed) |
| #4 Regex cache | `test_dna_sanity.py`, `test_performance.py` | `TestIUPACPattern::test_pattern_cache_*`, `TestIUPACPatternCachePerformance::test_warm_cache_is_near_free` |
| #5 Wrap midpoint | `test_circular_math.py`, `test_invariants_hypothesis.py` | `TestFeatureMidpoint::test_wrap_around_*`; `TestWrapMidpointProperties::*` (fuzzed) |
| #6 Circular wrap RE scan | `test_dna_sanity.py` | `TestRestrictionScan::test_circular_wraparound_*` |
| #7 Atomic saves | `test_data_safety.py` | `TestSafeSaveJson::*`, `TestSafeLoadJson::*`, `TestSchemaVersioning::*`, `TestRealFilesNeverTouched` |
| #8 `_feat_len` | `test_circular_math.py`, `test_invariants_hypothesis.py` | `TestFeatLen::*`; `TestFeatLenProperties::*`, `TestBpInProperties::*` (fuzzed) |
| #9 Wrap edit integrity | `test_edit_record.py` | (whole file) |
| #10 Undo deepcopy | `test_smoke.py` | `TestUndoSnapshotIndependence::*` |

### Test conventions

- **Cross-validate against Biopython** where possible (codon table, reverse-complement). If Biopython's standard table changes, the test fails noisily.
- **Hand-verifiable** test inputs — every restriction-site test uses a sequence short enough to count expected hits by eye.
- **Regression guards cite the date** — every test protecting a past bug has a docstring like `# Regression guard for 2026-03-30 fix`.
- **No network, no real files** — all tests use synthetic `SeqRecord`s and monkeypatched paths.
- **Performance budgets are LOOSE** (6–20× headroom). They catch architectural regressions, not micro-perf drift.
- **Property-based fuzzing** (`test_invariants_hypothesis.py`) complements hand-written regression tests. Use `@given` + `@settings(max_examples=..., deadline=None)` and `assume(...)` for filtering. Anchor every property to a sacred invariant so a Hypothesis failure maps to a concrete design contract.

### Adding a new test

1. Pick the right file (or add a new one).
2. For SeqRecord-based tests, use `tiny_record` fixture.
3. For Textual async tests: `async def test_*` (no decorator), `async with app.run_test(size=TERMINAL_SIZE) as pilot: await pilot.pause(); await pilot.pause(0.5)`. Double-pause is needed for `call_after_refresh` callbacks.
4. For perf tests, warm the cache then average 10–20 iterations.

## Performance notes

Key optimizations in place:

1. **Sidebar populate cascade suppressed** via `_populating` flag + `call_after_refresh` deferred reset — eliminates duplicate `_build_seq_text` per record load.
2. **Memoized `_build_seq_inputs()`** cached in module-level `_BUILD_SEQ_CACHE` (4-entry, identity-keyed). Cursor moves don't recompute.
3. **Per-chunk `str.translate`** for reverse strand instead of per-base. Module-level `_DNA_COMP_PRESERVE_CASE` avoids rebuilding the table each render.
4. **`_SCAN_CATALOG`** precomputed at import time eliminates per-scan `_rc` / `_iupac_pattern` / `len` calls.
5. **`_draw_cache`** on `PlasmidMap` — map render is only recomputed on size / mode / feature / RE-state change.

What was profiled but deliberately **not touched**: Textual compositor (framework), Rich `Text.append` (already efficient), import time (Textual + Rich dominate).

## Release + versioning

Versions live in `pyproject.toml` and `splicecraft.py::__version__`; `release.sh` updates both via sed. See `git log --oneline` for full release history. Recent: v0.3.1 (schema-versioned JSON envelope + crash-recovery autosave + per-plasmid undo stashes + Hypothesis property tests), v0.3.0 (Mutagenize modal with codon registry/harmonization), v0.2.8 (deep-copy record in undo/redo snapshots).

### Stubs still in menus (not implemented)
- **Build > Simulate Assembly** — `coming soon`
- **Build > New Part editor** — `coming soon`

## Known pitfalls

1. **Bare `except` is forbidden.** Use `except NoMatches` around `query_one`, `except ET.ParseError` around XML, `except (OSError, json.JSONDecodeError)` around file I/O. If you must catch `Exception`, `_log.exception(...)` it.
2. **Wrapped features (`end < start`) are first-class citizens.** Anywhere you compute distances, midpoints, or "is bp inside this feature", use the modular form via `_bp_in()` or `_feat_len()`. See sacred invariants #5, #6, #8, #9.
3. **Cache keys use `id(...)` of feature lists.** Correct *only* because the app reassigns lists on load rather than mutating them in-place. If you start mutating `self._feats` in-place, caches return stale renders.
4. **Textual reactive auto-invalidation depends on field assignment, not mutation.** `self._feats = new_list` triggers refresh; `self._feats.append(x)` does not.
5. **Single-file means giant diffs are normal.** When a refactor touches the rendering layer, expect 100+ line edits. The greppability tradeoff is worth it.
6. **Primer3 is linear-only.** For wrap regions, rotate template to `seq[start:] + seq[:start]` before calling, then unrotate positions via `(coord + rotation) % total`. See `_design_detection_primers`.
7. **`_source_path` is preserved through in-place edits.** Only cleared when `clear_undo=True` (fresh loads). Otherwise Ctrl+S after pLannotate or primer-add would forget the original file.
8. **NCBI responses go through `_safe_xml_parse`.** It rejects DOCTYPE/ENTITY before `ET.fromstring`. Don't add a new NCBI endpoint call without routing through it.

## How to extend — modular recipes

SpliceCraft is a single file on purpose, but new capabilities should still be **self-contained slabs** so the file stays navigable. Follow one of the recipes below.

### A. New pure helper function

Use case: new sequence transform, new analysis, new format. Pick this whenever the new code has no UI.

1. Place module-level helpers in the logically nearest section (use the Top-level structure table).
2. Name it `_snake_case` — leading underscore signals "internal, no public API guarantee".
3. Keep it **pure**: no globals, no logging, no UI.
4. Add a test in the matching `test_*.py` file. For bio logic, cross-validate against Biopython where possible.
5. If it's hot-path, add a `_performance.py` budget test.

### B. New persisted JSON library

Use case: a new user-facing collection (like parts bin, primers, codon tables).

1. Define `_MYTHING_FILE = _USER_DATA_DIR / "mything.json"` near the other four.
2. Write `_load_mything()` and `_save_mything(entries)` that route through `_safe_load_json` / `_safe_save_json` — **never** bypass these (sacred invariant #7). Envelope format + legacy back-compat come for free.
3. Filter `isinstance(entry, dict)` after load so hand-edited files can't crash `.get()` callers.
4. Add `_MYTHING_FILE` to the `_protect_user_data` autouse fixture in `tests/conftest.py`, plus a `_mything_cache` reset.
5. Cover corruption recovery in `test_data_safety.py` or a new `test_mything_io.py`.

### C. New modal dialog

Use case: a self-contained form that returns a result (file open, confirmation, parameter picker).

1. Subclass `ModalScreen[ReturnType]` (templates: `FetchModal`, `OpenFileModal`, `AminoAcidPickerModal`).
2. Implement `compose()` with the form layout (Horizontal / Vertical containers).
3. Use `query_one("#widget-id", WidgetType)` to read inputs. Wrap these in `except NoMatches` if mount order is unclear.
4. Call `self.dismiss(result)` to return. Escape should dismiss with `None`.
5. Push from the app: `self.push_screen(MyModal(args), callback=on_result)`.
6. Cover the modal in `test_smoke.py` — mount under `app.run_test`, assert widgets exist, drive `pilot.click` / `pilot.press` for a happy path.

### D. New heavy / background operation

Use case: anything that shouldn't block the UI loop — network fetch, subprocess, long compute.

1. Decorate with `@work(thread=True)` on a method of `PlasmidApp` (or the modal that owns it).
2. Wrap the body in `try / except Exception as exc`, log via `_log.exception(...)`, and push a user-friendly message with `self.app.call_from_thread(self._notify_err, exc)`.
3. Never touch widgets directly from the worker — always `call_from_thread`.
4. **If the worker captures mutable state** (e.g. `self._current_record`), capture the identity at entry and guard the callback with `if self._current_record is captured_record: ...`. Otherwise a fast user can apply your stale result on top of their newer record. Template: `PlasmidApp._run_plannotate_worker`.
5. **Re-entry guard** any worker the user can spam (like an "Annotate" button): set `self._myop_running = True` at entry, reset in a `finally` block.

### E. New menu action / keybinding

Use case: exposing a feature to the top menu or a global shortcut.

1. Add `action_my_thing(self)` on `PlasmidApp`.
2. Add a `Binding("key", "my_thing", "description")` to `PlasmidApp.BINDINGS`.
3. Add a menu item to the relevant entry in `MenuBar.compose()` — keep the letter-shortcut consistent with the existing style.
4. If the action opens a modal, delegate to recipe C. If it starts a worker, recipe D.

### F. New full-screen workbench (rare)

Use case: a standalone, modal-free workspace (like `PrimerDesignScreen`, `MutagenizeModal`).

1. Subclass `Screen` (not `ModalScreen`) for a permanent space; subclass `ModalScreen` for something dismissable.
2. Push with `self.push_screen(MyScreen(seq, feats, name))` from a menu action. Pop with `self.app.pop_screen()` or Escape.
3. Compose panels inside `Horizontal`/`Vertical` containers. Reuse widgets from the main app rather than cloning them.
4. Register keybindings on the screen itself via `BINDINGS` — they're scoped to the screen.

## Sister project reference (ScriptoScope)

`/home/seb/proteoscope/scriptoscope.py` (~8,600 lines) is the more mature sibling by the same author and source of most patterns here. When SpliceCraft hits scaling problems, check there first for pre-validated solutions:

| Pattern | When SpliceCraft would need it |
|---------|---|
| Thread-local `Console` for `_text_to_content` | If sequence-panel render starts blowing the 33 ms/frame budget |
| Two-level render cache (`_seq_render_cache` + `_content_cache`, LRU via `OrderedDict.move_to_end`) | If repainting on cursor moves becomes janky |
| `@lru_cache(1)` availability probes for optional CLI tools | If SpliceCraft shells out beyond pLannotate (e.g. BLAST, Prodigal) |

## Future work (user is undecided)

The user is weighing:
- **Merging** SpliceCraft, ScriptoScope, MitoShift, RefHunter, molCalc into one Textual app with multiple "modes"
- **Keeping them separate** as focused single-purpose apps and (optionally) extracting shared utilities into pure-Python modules

Either direction is viable. The single-file convention and shared logging/error patterns documented here keep the merge option open without forcing it.

## For future agents

1. **Read this file first.** It gives you architecture without reading 10k lines.
2. **Run `python3 -m pytest -q`** before and after any change. 508 tests, ~95 s. Biology subset (`tests/test_dna_sanity.py`) runs in < 1 s for a fast inner loop.
3. **Check `~/.local/share/splicecraft/logs/splicecraft.log`** (or `$SPLICECRAFT_LOG`) when debugging. Every session has a unique 8-char ID.
4. **Don't break the sacred invariants.** Each has a test (see mapping table). If you touch `_scan_restriction_sites`, `_rc`, `_iupac_pattern`, `_translate_cds`, `_bp_in`, `_feat_len`, the midpoint formula, or `_rebuild_record_with_edit`, the relevant tests will tell you immediately.
5. **Follow the error-handling convention**: `_log.exception` for stack traces, `notify()` or `Static.update("[red]...[/]")` for the user. Narrow `except` types. Never let raw tracebacks hit the TUI.
6. **When in doubt about real-world behavior** — eyeball it on pUC19 (`L09137`) and pACYC184 (`MW463917.1`), both fetched at first-run.
7. **Past fix history lives in git.** Use `git log --oneline` and `git show <hash>` rather than restoring fix-log sections to this file.
