# CLAUDE.md — AI Agent Context for SpliceCraft

Agent handoff document. Read before touching the codebase.

The project is developed by a human bioinformatician with an AI agent (Claude Opus 4.6+).

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer, sequence editor, and cloning/mutagenesis workbench** built with Python 3.10+ / Textual / Biopython. Renders Unicode braille-dot plasmid maps in the terminal, with a per-base sequence panel, restriction-site overlays, a **collection-driven plasmid library**, Golden Braid L0 assembly tooling, Primer3-backed primer design, and SOE-PCR site-directed mutagenesis.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs, user ATinyGreenCell). **PyPI:** `splicecraft`. Latest: **v0.4.0**.

- **Single-file architecture:** entire app is `splicecraft.py` (~16,200 lines). Intentional — keeps the codebase greppable. Sibling project ScriptoScope (~8,600 lines) follows the same convention.
- **Test suite:** 977 tests across 17 files in `tests/`. `pytest -n auto` ~3 min on 8 cores; sequential ~7 min. Biology subset (`test_dna_sanity.py`) < 1 s. `test_invariants_hypothesis.py` adds property-based fuzzing.
- **Dependencies:** `textual>=8.2.3`, `biopython>=1.87`, `primer3-py>=2.3.0`, `platformdirs>=4.2`. Tests: `pytest`, `pytest-asyncio`, `pytest-xdist`, `hypothesis`. No optional runtime deps — pLannotate integration was removed in 0.4.0.
- Releases via `./release.py X.Y.Z` (bumps version, runs tests, builds, tags, pushes; `publish.yml` uploads to PyPI via OIDC). Pure-Python — no bash/sed/grep dependencies.

## How to run

```bash
cd ~/SpliceCraft
python3 splicecraft.py              # empty canvas
python3 splicecraft.py L09137       # fetch pUC19 from NCBI
python3 splicecraft.py myplasmid.gb # open local GenBank file (.gb/.gbk/.dna)
python3 -m pytest -n auto -q        # full test suite (~3 min on 8 cores)
```

End users: `pipx install splicecraft && splicecraft`.

Logs: `~/.local/share/splicecraft/logs/splicecraft.log` (override with `$SPLICECRAFT_LOG`). Every line prefixed with an 8-char session ID for multi-run grepping.

## Architecture (single file: `splicecraft.py`)

### Top-level structure (line numbers ±50)

| Lines | Section |
|-------|---------|
| 1–200 | Docstring, imports (incl. module-level `from datetime import date as _date`), user data dir, dep check, rotating session-tagged logger, feature-colour palette |
| 201–460 | Atomic JSON persistence (`_safe_save_json` / `_safe_load_json` / `_extract_entries`; envelope schema `{"_schema_version":1,"entries":[...]}` + legacy bare-list back-compat) + library cache loaders |
| 461–575 | **Collections persistence** — `_load_collections` / `_save_collections` (deepcopy-on-read), `_get/_set_active_collection_name`, `_find_collection`, `_collection_name_taken`, `_ensure_default_collection` (Main Collection migration), `_sync_active_collection_plasmids`, `_restore_library_from_active_collection` |
| 576–1820 | NEB enzyme catalog (~204), IUPAC tables + cached regex, `_rc`, `_scan_restriction_sites`, `_assign_chunk_features`, `_render_feature_row_pair`, memoized `_build_seq_inputs`/`_build_seq_text`, OSC-52 clipboard, `_translate_cds` |
| 1821–1870 | Char-aspect detection + label helpers + `_cursor_row_key` DataTable utility |
| 1871–2110 | GenBank I/O (`fetch_genbank`, `load_genbank` auto-detect `.gb`/`.dna`, `_record_to_gb_text`, `_gb_text_to_record`, `_normalize_for_genbank`, `_export_genbank_to_path`, `_export_fasta_to_path`) |
| 2111–2200 | _(pLannotate slot — removed in 0.4.0; placeholder comment only)_ |
| 2201–2350 | `_Canvas` + `_BrailleCanvas` (sub-cell braille resolution) + `PlasmidMap` start |
| 2200–3160 | `PlasmidMap` widget — circular/linear draw, label placement, `_draw_cache`; `FeatureSidebar` |
| 3160–3730 | **`LibraryPanel`** — two-mode panel (collections list ↔ active-collection plasmids). Click a collection → enter plasmids view; back button (`←` in plasmids-button row) → return to collections. Per-mode `+ − ✎` CRUD; `set_dirty` does single-cell `update_cell_at` instead of full repopulate. |
| 3730–3990 | `SequencePanel` |
| 3990–4400 | Core modals (`EditSeqDialog`, `FetchModal`, `OpenFileModal`, `ExportGenBankModal`, `FastaExportModal`, `DropdownScreen`, `MenuBar`) |
| 4400–4900 | Golden Braid L0 position catalog (Esp3I/BsmBI overhangs); cloning grammar registry (`_BUILTIN_GRAMMARS`, `_load_custom_grammars`, `_get_active_grammar`); settings.json key/value persistence |
| 4900–5800 | Codon-usage registry, Kazusa parser, NCBI taxid search (`_safe_xml_parse`), CAI/GC. SOE-PCR mutagenesis primer design (`_mut_*`) |
| 5800–7100 | Feature-library workbench, `PlasmidFeaturePickerModal`, `AddFeatureModal`, `ColorPickerModal`, `_FeatureSnippetPanel`, `FeatureLibraryScreen` |
| 7100–7600 | `PartsBinModal`, FASTA file picker (`_FastaAwareDirectoryTree`, `FastaFilePickerModal`) |
| 7600–8400 | `_feats_for_domesticator` + `DomesticatorModal` (4-source part picker), `ConstructorModal` (Golden Braid L0 assembly UI) |
| 8400–8900 | `NcbiTaxonPickerModal`, `SpeciesPickerModal`, Mutagenize helpers (`_MutPreview`, `AminoAcidPickerModal`) |
| 8900–11800 | **`MutagenizeModal`** — 4-source CDS picker (map / library / parts bin / protein-harmonize); excludes "map" when no plasmid loaded |
| 11800–12750 | `PrimerDesignScreen` |
| 12750–13700 | Small modals (`UnsavedQuitModal`, `UnsavedNavigateModal`, `PlasmidPickerModal`, `RenamePlasmidModal`, `LibraryDeleteConfirmModal`, **`CollectionsModal`** for snapshot save, **`CollectionNameModal`** for create/rename, **`CollectionDeleteConfirmModal`**) |
| 13700–end | `PlasmidApp` — main controller, keybindings, undo/redo stashes (per-plasmid LRU), autosave, `@work` threads (NCBI fetch, library seed, Kazusa); `_discard_changes` (preserves `_source_path`); migration call site in `compose()`; `main()` |

### Key design patterns

- **Rich `Text` for all rendering** — no curses.
- **Braille canvas** gives sub-character pixel resolution (2×4 dots per terminal cell).
- **Feature coordination:** map click → sidebar highlight → sequence scroll (and back via Textual messages).
- **Undo/redo:** snapshot-based (full seq + cursor + `deepcopy` of SeqRecord), max 50. **Per-plasmid stashes** — switching plasmids stashes outgoing history under the old `record.id`, restores incoming. LRU-capped at 10 plasmids.
- **Crash-recovery autosave:** dirty edits debounce a 3 s write to `_DATA_DIR/crash_recovery/{safe_id}.gb`. Cleared on save/abandon. Startup notifies on survivors.
- **Caching:** `PlasmidMap._draw_cache`, `_BUILD_SEQ_CACHE`, `_PATTERN_CACHE`, `_SCAN_CATALOG` — keyed on inputs (using `id(self._feats)` since lists are reassigned, not mutated, on load). `_collections_cache` deepcopies on read so callers can mutate entries safely.
- **Workers:** `@work(thread=True)` for NCBI fetch, library seed, Kazusa codon fetch. Results pushed back via `call_from_thread` with stale-record guards.

## Logging & error handling

```python
_log = logging.getLogger("splicecraft")  # rotating 2MB×2 at _DATA_DIR/logs/splicecraft.log
```

- **User-facing errors** → `self.notify(...)` or `Static.update("[red]...[/]")`. Never raw tracebacks.
- **Diagnostic detail** → `_log.exception("context: %s", ...)` inside `except` blocks.
- **Worker errors** → `_log.exception` then push a friendly message via `call_from_thread`.
- **Narrow `except`.** Use `except NoMatches:` around `query_one`, `except ET.ParseError:` around XML, `except (OSError, json.JSONDecodeError):` around file I/O. Reserve bare `except Exception` for worker bodies — and always log there.

## Sacred invariants (DO NOT BREAK)

Each has at least one test. Mapping at end of file.

1. **Palindromic enzymes are scanned forward only.** `_scan_restriction_sites` skips reverse scan for palindromes; adds bottom-strand `recut`. Scanning both strands double-counts every site.

2. **Reverse-strand resite positions use the forward coordinate.** A reverse-strand hit at position `p` (after RC) is stored as `p`, not `n - p - site_len`. Cut maps via `site_len - 1 - fwd_cut`.

3. **`_rc()` handles full IUPAC.** Reverse-complement translates ambiguity codes (R, Y, W, S, M, K, B, D, H, V, N) via `_IUPAC_COMP`, not just ACGT.

4. **IUPAC regex patterns are cached.** `_iupac_pattern()` uses `_PATTERN_CACHE` to avoid recompiling ~200 patterns on every restriction scan.

5. **Circular wrap-around midpoints.** `arc_len = (end - start) % total` then `(start + arc_len // 2) % total`. Naive `(start + (end - start) // 2) % total` puts the label opposite the actual arc for wrapped features.

6. **Circular wrap-around restriction scanning.** `_scan_restriction_sites(circular=True)` scans `seq + seq[:max_site_len-1]`. Each wrap-around hit emits **two resite pieces** (labeled tail `[p, n)` + unlabeled head `[0, (p+site_len) - n)`) and **one recut** at `(p + fwd_cut) % n`. Code that counts resites for filtering must count only labeled pieces.

7. **Data-file saves always back up.** `_safe_save_json` writes `.bak` of existing file before replacing, via `tempfile.mkstemp` + `os.fsync` + `os.replace`. Shrink guard warns on writing fewer entries. Envelope format `{"_schema_version": 1, "entries": [...]}`; loaders accept both envelope and legacy bare-list (pre-0.3.1) via `_extract_entries`. Future-version writes warn but still load. **Never bypass `_safe_save_json`.**

8. **Wrap-aware feature length everywhere.** Use `_feat_len(start, end, total)` — returns `(total - start) + end` when `end < start`, else `end - start`. All sort keys, length displays, biology checks must route through it.

9. **Wrap-feature integrity in record edits.** `int(CompoundLocation.start)` returns `min(parts.start)` and `int(.end)` returns `max(parts.end)`, silently flattening wrap features. `_rebuild_record_with_edit` per-part shifts wrap features and only collapses to FeatureLocation when 1 part survives. Zero-width post-edit features dropped.

10. **Undo snapshots must be deepcopied.** `_push_undo`, `_action_undo`, `_action_redo` all `deepcopy(self._current_record)`.

## Plasmid collections (the new top-level organizational unit, 0.4.0+)

A **collection** is a named bucket of plasmids. Users create one on demand from the LibraryPanel; the active collection's plasmids are what shows in the panel's plasmids view.

**On-disk model:**
- `collections.json` is the single source of truth for collection identity + their plasmid contents.
- `plasmid_library.json` is a **live mirror** of the active collection — every `_save_library` call also calls `_sync_active_collection_plasmids`, which writes the same entries back into the active collection in `collections.json`.
- The active-collection name lives in `settings.json` as `{"key": "active_collection", "value": "Main Collection"}`.

**Migration / first-run** (`_ensure_default_collection`): runs in `App.compose()` (BEFORE children mount, since Textual fires mount events leaves→root). On first launch with no `collections.json`, wraps whatever's in `plasmid_library.json` in a "Main Collection" and points active at it. On subsequent launches, sets active to the first collection if it's missing (orphaned from a delete-active-collection flow).

**Restore** (`_restore_library_from_active_collection`): also runs in `compose()`. Writes the active collection's plasmids into `plasmid_library.json` so the panel renders the user's stored data even if they edited collections.json externally. Bypasses `_save_library`'s mirror — collection is the source.

**LibraryPanel two-mode UI:**
- *Collections view*: DataTable with `Name | Plasmids`. `+` (new collection prompt), `−` (delete with confirm), `✎` (rename). Click a row → enter that collection's plasmids view + set active.
- *Plasmids view*: existing per-plasmid table. `+` (save loaded record), `−` (remove), `←` (back to collections; prompts via `UnsavedNavigateModal` if `app._unsaved` is true), `✎` (rename plasmid).
- Mode persists in `_view_mode`. Returning user starts in plasmids view if active is set; first-runs / no-active start in collections view.
- The dirty marker `*` shows in both view headers — even from collections view, the user knows there are unsaved edits.

**Cache contract:** `_load_collections` and `_save_collections` `deepcopy` so a caller mutating a returned dict (e.g. `_btn_coll_rename` editing `c["name"]` in place before save) cannot poison `_collections_cache`. Same contract as `_load_features` per the existing convention.

**`UnsavedNavigateModal`** is the sibling of `UnsavedQuitModal` — same Save/Discard/Cancel structure but the verbs say "go back to collections" instead of "quit". `App._discard_changes` reloads the record from the library copy, clears undo/redo, and **preserves `_source_path`** (calls `_apply_record(record, clear_undo=False)` then restores source_path) so post-discard Ctrl+S still targets the user's original .gb file.

## Mutagenize source picker (0.3.10+)

`MutagenizeModal` accepts CDS DNA from any of four sources via a top-of-modal `Select`:

1. **Current map features** — only shown when `template_seq` is non-empty (modal launchable from blank canvas).
2. **Plasmid library** — pick a library entry, then pick a CDS feature on it.
3. **Parts bin** — pick a domesticated GB part directly. Filtered to `len % 3 == 0 && >= 30 bp` (same gate as map/library sources). The part's stored `sequence` is treated as a single-CDS pseudo-plasmid spanning `[0, len)`.
4. **Protein sequence (harmonize)** — paste 1-letter AA, harmonize via the active codon table (`_codon_harmonize` + `_codon_fix_sites` for BsaI scrubbing).

`_initial_source` and `_src_options` are computed in `__init__` (was previously set in compose, read in on_mount via `getattr` fallback — refactored 0.4.0). `_reset_cds_state(info_msg="")` consolidates the "clear CDS + primer state + refresh preview" block used by source-switch and parts-deselect.

## Cloning grammars (Golden Braid, MoClo, custom)

Every Type IIS-aware tool reads its overhangs / enzyme / forbidden-sites / coding-types / type→INSDC map from a **grammar dict**. Two ship as built-ins (`_BUILTIN_GRAMMARS["gb_l0"]`, `_BUILTIN_GRAMMARS["moclo_plant"]`); user-defined grammars persist to `cloning_grammars.json` and become editable in `GrammarEditorModal`. The active grammar id lives in `settings.json` (`{"key": "active_grammar", "value": "gb_l0"}`); `_get_active_grammar()` resolves the id to a dict and falls back to gb_l0 if the persisted id no longer exists (e.g., a custom grammar was deleted while still selected) — and writes the recovery back to settings so we don't keep falling back forever.

**Grammar schema:** `id`, `name`, `enzyme`, `site`, `spacer`, `pad`, `forbidden_sites: {enzyme: site}`, `positions: [{name, type, oh5, oh3, color}]`, `coding_types: [str]` (eligible for codon-fix repair), `type_to_insdc: {gb_type: insdc_type}`, `catalog: [(name, type, position, oh5, oh3, backbone, marker)]`, `editable: bool`. Built-ins set `editable=False` so the editor refuses to save over them — fork via "Duplicate as Custom" in the Parts Bin first.

**PartsBinModal integration.** Top-of-modal Select widget (`#parts-grammar-select`) flips the active grammar; the overhang table (`#parts-overhangs`) shows the active grammar's positions; `_all_rows` filters catalog + user parts to the active grammar's id. New parts saved via DomesticatorModal pick up `grammar=active_id` automatically; **legacy parts without a `grammar` field default to `gb_l0`** so existing v0.3.x data migrates intact. The Edit button opens `GrammarEditorModal` (read-only on built-ins); Duplicate forks the active grammar with a user-supplied name.

**DomesticatorModal integration.** The Position dropdown is built from `active_grammar["positions"]`; `_design_gb_primers` accepts a `grammar=` kwarg threaded through from the modal so the primer tail uses the active grammar's `pad + site + spacer + oh5` and the codon-fix repair scans `forbidden_sites` (which differs per grammar — GB L0 forbids Esp3I + BsaI; MoClo Plant forbids BsaI + BpiI). `_simulate_primed_amplicon` likewise accepts a `grammar=` kwarg; the Parts Bin "Copy Primed Sequence" button looks up the part's stored grammar so an old MoClo part still gets BsaI tails even after the user has flipped the active grammar to GB.

**Cache + deepcopy.** `_load_custom_grammars` and `_load_collections` both deepcopy on read so caller-side mutations of returned dicts don't poison the cache.

## On-disk JSON format (schema v1)

All seven libraries (`plasmid_library.json`, `parts_bin.json`, `primers.json`, `codon_tables.json`, `features.json`, `feature_colors.json`, `collections.json`) plus `settings.json` and `cloning_grammars.json` use:

```json
{"_schema_version": 1, "entries": [...]}
```

**Legacy compatibility.** Pre-0.3.1 wrote bare lists; `_extract_entries` accepts both, silently rewrites as envelope on next save. When bumping `_CURRENT_SCHEMA_VERSION`, teach `_extract_entries` to migrate forward in the loader. Newer-version files load with a warning so users know fields may drop on save.

## Crash-recovery autosave

Dirty edits trigger a 3 s debounced write to `_CRASH_RECOVERY_DIR/{safe_id}.gb` (default `~/.local/share/splicecraft/crash_recovery/`). Deleted on `_mark_clean` or abandon. `_check_crash_recovery()` at startup notifies on survivors; user recovers via File > Open.

- `_autosave_path(record)` sanitises `record.id` with `re.sub(r'[^A-Za-z0-9._-]', '_', ...)`, caps at 80 chars, and appends a 6-char sha256 to disambiguate IDs that collide after sanitisation.
- Atomic write — `tempfile.mkstemp` + `os.replace`.
- Best-effort — `except Exception: _log.exception(...)`. Autosave is a safety net, not source of truth.
- Debounced via `self.set_timer`. `_mark_dirty` restarts countdown; `_mark_clean` cancels by deleting target.

## Per-plasmid undo/redo stashes

`_apply_record(clear_undo=True)` (switch-plasmid path) stashes outgoing stacks under `record.id` in `_stashed_undo_stacks` / `_stashed_redo_stacks`, restores incoming history if previously edited. LRU-capped at `_MAX_PLASMIDS_WITH_UNDO = 10`. `_current_undo_key` tracks the live stack. `clear_undo=False` (in-place edits — primer-add, feature merge, **discard-from-library**) leaves stacks intact.

## Test suite

```bash
python3 -m pytest -n auto -q                          # full, parallel (~3 min)
python3 -m pytest -q                                  # serial (~7 min) — debugging
python3 -m pytest tests/test_dna_sanity.py            # biology only (< 1 s)
python3 -m pytest tests/test_invariants_hypothesis.py # property-based fuzzing
python3 -m pytest -k "palindrome"                     # filter
python3 -m pytest -x                                  # stop on first failure
```

Parallel runs rely on `pytest-xdist` + the autouse `_protect_user_data` fixture (per-test `tmp_path` isolation; monkeypatches `_LIBRARY_FILE`, `_PARTS_BIN_FILE`, `_PRIMERS_FILE`, `_CODON_TABLES_FILE`, `_FEATURES_FILE`, `_FEATURE_COLORS_FILE`, `_GRAMMARS_FILE`, `_SETTINGS_FILE`, `_COLLECTIONS_FILE`, `_CRASH_RECOVERY_DIR` and caches). **No test can write to real user files.** Module-level read-only caches (`_BUILD_SEQ_CACHE`, `_PATTERN_CACHE`, `_SCAN_CATALOG`) are safe — nothing writes them at test time.

`pyproject.toml` sets `asyncio_mode = "auto"` so async tests don't need `@pytest.mark.asyncio`. `tests/conftest.py` provides `tiny_record` / `tiny_gb_path` / `isolated_library` fixtures.

| File | Tests | Covers |
|------|------:|--------|
| `test_dna_sanity.py` | 74 | Sacred invariants 1–6; Type IIS cut-outside-recognition; `_translate_cds` |
| `test_primers.py` | 60 | Detection / cloning / Golden Braid / generic; wrap-region template rotation |
| `test_genbank_io.py` | 68 | `load_genbank` round-trip (GenBank + SnapGene `.dna`); JSON corruption recovery; `_export_fasta_to_path` |
| `test_smoke.py` | ~50 | Textual mounts; rotation / view / RE toggles; per-plasmid undo stashes; crash-recovery autosave |
| `test_mutagenize.py` | 49 | SOE-PCR primers, codon substitution, CAI round-trips |
| `test_codon.py` | 50+ | Codon registry, harmonization, Kazusa parser, NCBI XML safety, CAI/GC math; Mutagenize 4-source flow (incl. Parts Bin source); Mutagenize-without-plasmid |
| `test_domesticator.py` | 258 | Golden Braid L0 positions; 4-source picker; `_feats_for_domesticator`; FASTA picker; cloning simulator; codon-fix repair (multi-site, cascade-prevention, binding-region advisory); Save Primers (`pairs` list, DOM suffix); Save As Feature button (GB→INSDC type map, builtin reject, persist round-trip); `_feature_library_match` helper + `_features_generation` counter + `_build_feature_library_index`; "Feat Lib" column (exact-green / name-yellow / empty) with index-cached lookups (rebuilt only on gen advance); Save-As-Feature warnings on collision; **grammar abstraction** (`_BUILTIN_GRAMMARS`, `_all_grammars`, `_get_active_grammar` fallback, `_grammar_position_by_type`); settings.json + cloning_grammars.json round-trips with deepcopy isolation; Parts Bin grammar dropdown filters by active grammar (legacy parts → gb_l0); DomesticatorModal honours active grammar (overhangs / enzyme / forbidden-site scrub) |
| `test_collections.py` | 47 | Collections persistence + envelope schema + legacy bare-list back-compat; Main Collection migration (first-run, idempotent, fallback active); `_save_library` mirrors to active; `_sync_active_collection_plasmids` no-op when no active; LibraryPanel two-mode toggle; collection CRUD via panel (add/remove/rename); Save-loaded-plasmid-to-active-collection; Back-button unsaved-prompt (Save/Discard/Cancel); CollectionsModal Save/Load/Delete |
| `test_circular_math.py` | 38 | Sacred invariant #5; `_bp_in` / `_feat_len` |
| `test_data_safety.py` | 47 | Sacred invariant #7; envelope round-trip + legacy back-compat + future-version warning; `_atomic_write_text`; `_do_save` atomicity; collections.json isolation |
| `test_add_feature.py` | 24 | AddFeatureModal: qualifier round-trip, validation, save-to-library dedup, insert-at-cursor |
| `test_modal_boundaries.py` | 29 | Every modal fits in 160×48 (and AddFeatureModal at 100×30); CollectionsModal, CollectionNameModal, CollectionDeleteConfirmModal, UnsavedNavigateModal included |
| `test_feature_library_screen.py` | 95 | Workbench CRUD + 4-step strand cycle; deferred-save / dirty-tracking / UnsavedQuitModal-on-close; Edit-button prefill round-trip; AddFeatureModal Orientation + Color; Ctrl+Shift+F capture (drag-matches-feature enrichment); ColorPickerModal xterm grid + drag-to-preview; Export-FASTA |
| `test_features_library.py` | 29 | JSON round-trip; `_GENBANK_FEATURE_TYPES`; per-entry `color` + `strand=0`; `_resolve_feature_color` precedence |
| `test_edit_record.py` | 14 | Sacred invariant #9: wrap features survive insert/replace as CompoundLocation |
| `test_invariants_hypothesis.py` | 11 | Property-based fuzzing of invariants #3, #5, #8 |
| `test_performance.py` | 9 | Loose budgets (4–20× headroom): pUC19 scan < 30 ms, 10 kb scan < 150 ms, etc. |

### Sacred invariant → test mapping

| Inv | File | Method |
|---|---|---|
| #1 Palindrome forward | `test_dna_sanity.py` | `TestRestrictionScan::test_ecori_single_site_not_double_counted`, `::test_palindromes_produce_one_recut_per_site` |
| #2 Reverse-strand fwd coord | `test_dna_sanity.py` | `::test_non_palindrome_on_reverse_strand_uses_forward_coordinate` |
| #3 `_rc()` IUPAC | `test_dna_sanity.py`, `test_invariants_hypothesis.py` | `TestReverseComplement::*`; `TestReverseComplementProperties::*` |
| #4 Regex cache | `test_dna_sanity.py`, `test_performance.py` | `TestIUPACPattern::*`, `TestIUPACPatternCachePerformance::*` |
| #5 Wrap midpoint | `test_circular_math.py`, `test_invariants_hypothesis.py` | `TestFeatureMidpoint::*`; `TestWrapMidpointProperties::*` |
| #6 Circular wrap RE scan | `test_dna_sanity.py` | `TestRestrictionScan::test_circular_wraparound_*` |
| #7 Atomic saves | `test_data_safety.py` | `TestSafeSaveJson::*`, `TestSchemaVersioning::*`, `TestRealFilesNeverTouched` |
| #8 `_feat_len` | `test_circular_math.py`, `test_invariants_hypothesis.py` | `TestFeatLen::*`; `TestFeatLenProperties::*`, `TestBpInProperties::*` |
| #9 Wrap edit | `test_edit_record.py` | (whole file) |
| #10 Undo deepcopy | `test_smoke.py` | `TestUndoSnapshotIndependence::*` |

### Conventions

- Cross-validate against Biopython where possible. Hand-verifiable inputs (short enough to count hits by eye).
- Regression guards cite the date in their docstring (`# Regression guard for 2026-03-30 fix`).
- No network, no real files. Synthetic `SeqRecord` + monkeypatched paths.
- Performance budgets are loose (4–20× headroom). They catch architectural regressions, not micro-perf drift.
- Property-based fuzzing in `test_invariants_hypothesis.py` — anchor every property to a sacred invariant.
- Async Textual tests: `async def test_*`, `async with app.run_test(size=...) as pilot: await pilot.pause(); await pilot.pause(0.5)` (double pause for `call_after_refresh`).

## Performance notes

1. Sidebar populate cascade suppressed via `_populating` flag + `call_after_refresh` deferred reset.
2. `_build_seq_inputs()` cached in 4-entry identity-keyed `_BUILD_SEQ_CACHE`.
3. Per-chunk `str.translate` for reverse strand (module-level `_DNA_COMP_PRESERVE_CASE`).
4. `_SCAN_CATALOG` precomputed at import — eliminates per-scan `_rc` / `_iupac_pattern`.
5. `PlasmidMap._draw_cache` — only recomputed on size / mode / feature / RE-state change.
6. **`LibraryPanel.set_dirty`** — early-returns when dirty state didn't change; updates only the active row's Name cell via `update_cell_at(Coordinate(row, 0))` instead of rebuilding the whole DataTable on every keystroke. Falls back to `_repopulate_plasmids` if the incremental API isn't available.

Profiled but **not touched**: Textual compositor, Rich `Text.append`, import time.

## Release + versioning

Versions in `pyproject.toml` and `splicecraft.py::__version__`; `release.py` rewrites both via in-file regex (one match each — refuses to bump if the file's formatting drifted enough that zero or >1 lines match). See `git log --oneline` for full release history.

**Stubs in menus (not implemented):**
- Build > Simulate Assembly — `coming soon`
- Build > New Part editor — `coming soon`

## Known pitfalls

1. **Bare `except` is forbidden.** Use narrow types. `_log.exception` if catching `Exception`.
2. **Wrapped features (`end < start`) are first-class.** Use `_bp_in()` / `_feat_len()` for any distance, midpoint, or "is bp inside" check. See invariants #5, #6, #8, #9.
3. **Cache keys use `id(...)` of feature lists.** Correct only because lists are reassigned on load, not mutated. Don't start mutating `self._feats` in-place.
4. **Textual reactive auto-invalidation requires assignment, not mutation.** `self._feats = new_list` triggers refresh; `self._feats.append(x)` does not.
5. **Single-file means giant diffs are normal.** Rendering-layer refactors touch 100+ lines.
6. **Primer3 is linear-only.** For wrap regions, rotate template to `seq[start:] + seq[:start]` then unrotate via `(coord + rotation) % total`. See `_design_detection_primers`.
7. **`_source_path` survives in-place edits.** Cleared only when `clear_undo=True` (fresh loads). Otherwise Ctrl+S after primer-add or **Discard-from-library** still targets the original .gb file. `_discard_changes` explicitly stashes/restores `_source_path` around its `_apply_record(clear_undo=False)` call to honor this.
8. **NCBI responses go through `_safe_xml_parse`.** It rejects DOCTYPE/ENTITY before `ET.fromstring`. Don't add a new NCBI endpoint without it.
9. **Migration runs in `App.compose()`, not `on_mount`.** Mount events fire leaves→root, so anything in `App.on_mount` runs AFTER `LibraryPanel.on_mount`. Collections + active-collection setup must be done before children mount or the panel reads stale state.
10. **`_save_library` mirrors to the active collection.** Every panel CRUD writes BOTH `plasmid_library.json` and `collections.json` (each with its own `.bak`). This is intentional — the two files are kept in sync. Routing a write around `_save_library` (e.g. `_restore_library_from_active_collection`) bypasses the mirror; do that only when the collection IS the source.

## How to extend — modular recipes

### A. New pure helper

Place in nearest section per Top-level structure. Snake-case with leading underscore. Pure: no globals, no logging, no UI. Add a test (cross-validate against Biopython where biological). Hot-path? Add a `test_performance.py` budget.

### B. New persisted JSON library

Define `_MYTHING_FILE = _USER_DATA_DIR / "mything.json"`. Route load/save through `_safe_load_json` / `_safe_save_json` — never bypass (invariant #7). Filter `isinstance(entry, dict)` after load. Add the file + its cache to `_protect_user_data` in `tests/conftest.py` and to `_check_data_files` in `PlasmidApp`. Cover corruption recovery in `test_data_safety.py`. If callers will mutate returned entries (rename in place, etc.), `deepcopy` on read like `_load_collections` / `_load_features`.

### C. New modal

Subclass `ModalScreen[ReturnType]` (templates: `FetchModal`, `OpenFileModal`, `CollectionNameModal`). `query_one("#widget-id", WidgetType)` reads, `self.dismiss(result)` returns. Push via `self.push_screen(MyModal(args), callback=on_result)`. Use `_cursor_row_key(table)` for DataTable cursor reads. Cover happy path in `test_smoke.py` or a topic-specific test file. Add a row to `test_modal_boundaries.py::_MODAL_CASES`.

### D. New background worker

`@work(thread=True)` on `PlasmidApp` or owning modal. Wrap body in `try / except Exception as exc`, `_log.exception`, push friendly message via `call_from_thread`. **Stale-record guard:** capture `self._current_record` identity at entry; `if self._current_record is captured_record:` in callback. **Re-entry guard:** `self._myop_running` flag with `finally`. Template: `PlasmidApp._seed_default_library` (background NCBI fetch with stale-record guard).

### E. New menu action / keybinding

Add `action_my_thing(self)` on `PlasmidApp`. Add `Binding("key", "my_thing", "desc")` to `BINDINGS`. Add menu item to `MenuBar.compose()` (or to the dropdown dict in `open_menu`). Modal → recipe C; worker → recipe D.

### F. New full-screen workbench

Subclass `Screen` (or `ModalScreen` if dismissable). Push from menu action with state passed in. Reuse main-app widgets. Register screen-scoped `BINDINGS`.

## Sister project (ScriptoScope)

`/home/seb/proteoscope/scriptoscope.py` (~8,600 lines) — same author, same single-file convention. Patterns to crib if SpliceCraft scales:

- Thread-local `Console` for `_text_to_content` (if seq-panel render blows the 33 ms/frame budget).
- Two-level render cache (`_seq_render_cache` + `_content_cache`, LRU via `OrderedDict.move_to_end`).
- `@lru_cache(1)` availability probes for optional CLI tools (BLAST, Prodigal).

User is undecided whether to merge SpliceCraft / ScriptoScope / MitoShift / RefHunter / molCalc into one Textual app with modes, or keep them separate. Either is viable — the single-file convention keeps the option open.

## For future agents

1. Read this file first.
2. `python3 -m pytest -n auto -q` before and after any change. Biology subset (`tests/test_dna_sanity.py`) gives a < 1 s inner loop.
3. Logs: `~/.local/share/splicecraft/logs/splicecraft.log` (or `$SPLICECRAFT_LOG`). 8-char session ID per run.
4. Don't break sacred invariants — touching `_scan_restriction_sites`, `_rc`, `_iupac_pattern`, `_translate_cds`, `_bp_in`, `_feat_len`, the midpoint formula, or `_rebuild_record_with_edit` will trip the relevant tests immediately.
5. Follow the error convention: `_log.exception` for stack traces, `notify()` / `Static.update("[red]...[/]")` for the user. Never raw tracebacks.
6. Eyeball real-world behaviour on pUC19 (`L09137`) and pACYC184 (`MW463917.1`).
7. Past fix history lives in git — `git log --oneline` and `git show <hash>`.
