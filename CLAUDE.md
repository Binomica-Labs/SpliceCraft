# CLAUDE.md — AI Agent Context for SpliceCraft

Agent handoff document. Read before touching the codebase.

The project is developed by a human bioinformatician with an AI agent (Claude Opus 4.6+).

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer, sequence editor, and cloning/mutagenesis workbench** built with Python 3.10+ / Textual / Biopython. Renders Unicode braille-dot plasmid maps in the terminal, with a per-base sequence panel, restriction-site overlays, a plasmid library, Golden Braid L0 assembly tooling, Primer3-backed primer design, and SOE-PCR site-directed mutagenesis.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs, user ATinyGreenCell). **PyPI:** `splicecraft`. Latest: **v0.3.3**.

- **Single-file architecture:** entire app is `splicecraft.py` (~14,300 lines). Intentional — keeps the codebase greppable. Sibling project ScriptoScope (~8,600 lines) follows the same convention.
- **Test suite:** 874 tests across 16 files in `tests/`. `pytest -n auto` ~125 s on 8 cores; sequential ~400 s. Biology subset (`test_dna_sanity.py`) < 1 s. `test_invariants_hypothesis.py` adds property-based fuzzing.
- **Dependencies:** `textual>=8.2.3`, `biopython>=1.87`, `primer3-py>=2.3.0`, `platformdirs>=4.2`. Tests: `pytest`, `pytest-asyncio`, `pytest-xdist`, `hypothesis`. **Optional runtime:** `pLannotate` (conda, GPL-3) for Shift+A annotation.
- Releases via `./release.sh X.Y.Z` (bumps version, runs tests, builds, tags, pushes; `publish.yml` uploads to PyPI via OIDC).

## How to run

```bash
cd ~/SpliceCraft
python3 splicecraft.py              # empty canvas
python3 splicecraft.py L09137       # fetch pUC19 from NCBI
python3 splicecraft.py myplasmid.gb # open local GenBank file (.gb/.gbk/.dna)
python3 -m pytest -n auto -q        # full test suite (~2 min on 8 cores)
```

End users: `pipx install splicecraft && splicecraft`.

Logs: `~/.local/share/splicecraft/logs/splicecraft.log` (override with `$SPLICECRAFT_LOG`). Every line prefixed with an 8-char session ID for multi-run grepping.

### pLannotate (optional)

Press **Shift+A** (or click ◈ in the library panel) to annotate via pLannotate. SpliceCraft calls it as a subprocess only — never imports it (GPL-3 license isolation).

```bash
conda create -n plannotate -c conda-forge -c bioconda plannotate
conda activate plannotate
plannotate setupdb          # downloads ~500 MB of BLAST/diamond DBs
```

If pLannotate is not on PATH, Shift+A notifies the user — nothing crashes.

## Architecture (single file: `splicecraft.py`)

### Top-level structure (line numbers ±30)

| Lines | Section |
|-------|---------|
| 1–200 | Docstring, imports, user data dir, dep check, rotating session-tagged logger, feature-colour palette |
| 201–385 | Atomic JSON persistence (`_safe_save_json` / `_safe_load_json` / `_extract_entries`; envelope schema `{"_schema_version":1,"entries":[...]}` + legacy bare-list back-compat) |
| 386–408 | Library cache loaders |
| 409–1448 | NEB enzyme catalog (~204), IUPAC tables + cached regex, `_rc`, `_scan_restriction_sites`, `_assign_chunk_features`, `_render_feature_row_pair`, memoized `_build_seq_inputs`/`_build_seq_text`, OSC-52 clipboard, `_translate_cds` |
| 1449–1521 | Char-aspect detection + label helpers |
| 1522–1659 | GenBank I/O (`fetch_genbank`, `load_genbank` auto-detect `.gb`/`.dna`, `_record_to_gb_text`, `_gb_text_to_record`) |
| 1660–1875 | pLannotate subprocess integration |
| 1876–1985 | `_Canvas` + `_BrailleCanvas` (sub-cell braille resolution) |
| 1986–2753 | `PlasmidMap` widget — circular/linear draw, label placement, `_draw_cache` |
| 2754–2868 | `FeatureSidebar` |
| 2869–3036 | `LibraryPanel` |
| 3037–3485 | `SequencePanel` |
| 3486–3825 | Core modals (`EditSeqDialog`, `FetchModal`, `OpenFileModal`, `DropdownScreen`) |
| 3826–3867 | `MenuBar` |
| 3868–4076 | Golden Braid L0 position catalog (Esp3I/BsmBI overhangs) |
| 4077–4130 | Parts-bin + primer-library persistence |
| 4131–4925 | Codon-usage registry, Kazusa parser, NCBI taxid search (`_safe_xml_parse`), CAI/GC. Crash-recovery config sits at top of slab |
| 4926–5437 | SOE-PCR mutagenesis primer design (`_mut_*`) |
| 5953–6394 | `PlasmidFeaturePickerModal`, `AddFeatureModal` |
| 6395–7162 | Feature library workbench (`ColorPickerModal`, `_FeatureSnippetPanel`, `FeatureLibraryScreen`) |
| 7163–7345 | `PartsBinModal` |
| 7346–7550 | FASTA file picker (`_FastaAwareDirectoryTree`, `FastaFilePickerModal`) |
| 7494–8003 | `_feats_for_domesticator` + `DomesticatorModal` (4-source part picker) |
| 8004–8322 | `ConstructorModal` (Golden Braid L0 assembly UI) |
| 8323–8660 | `NcbiTaxonPickerModal`, `SpeciesPickerModal` |
| 8661–8891 | Mutagenize helpers (`_MutPreview`, `AminoAcidPickerModal`) |
| 8892–9536 | `MutagenizeModal` |
| 9537–10635 | `PrimerDesignScreen` |
| 10636–10856 | Small modals (`UnsavedQuitModal`, `PlasmidPickerModal`, `RenamePlasmidModal`, `LibraryDeleteConfirmModal`) |
| 10857–end | `PlasmidApp` — main controller, keybindings, undo/redo stashes, autosave, `@work` threads; `main()` |

### Key design patterns

- **Rich `Text` for all rendering** — no curses.
- **Braille canvas** gives sub-character pixel resolution (2×4 dots per terminal cell).
- **Feature coordination:** map click → sidebar highlight → sequence scroll (and back via Textual messages).
- **Undo/redo:** snapshot-based (full seq + cursor + `deepcopy` of SeqRecord), max 50. **Per-plasmid stashes** — switching plasmids stashes outgoing history under the old `record.id`, restores incoming. LRU-capped at 10 plasmids.
- **Crash-recovery autosave:** dirty edits debounce a 3 s write to `_DATA_DIR/crash_recovery/{safe_id}.gb`. Cleared on save/abandon. Startup notifies on survivors.
- **Caching:** `PlasmidMap._draw_cache`, `_BUILD_SEQ_CACHE`, `_PATTERN_CACHE`, `_SCAN_CATALOG` — keyed on inputs (using `id(self._feats)` since lists are reassigned, not mutated, on load).
- **Workers:** `@work(thread=True)` for NCBI fetch, library seed, pLannotate, Kazusa codon fetch. Results pushed back via `call_from_thread` with stale-record guards.

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

## Core helper catalog

Load-bearing pure functions. Read these first before touching rendering, primer design, or the record pipeline.

| Helper | Line | Purpose |
|---|---:|---|
| `_safe_save_json` / `_safe_load_json` / `_extract_entries` | 282 / 362 / 228 | Atomic JSON I/O with `.bak` recovery + schema envelope. All six libraries route through these. |
| `_atomic_write_text` | 251 | Module-level atomic write (tempfile + fsync + os.replace + cleanup). `_do_save`, `_do_autosave`, `_export_genbank_to_path`, `_export_fasta_to_path` use it. |
| `_iupac_pattern` | 680 | IUPAC→regex compiler, cached in `_PATTERN_CACHE`. |
| `_IUPAC_COMP`, `_DNA_COMP_PRESERVE_CASE` | ~690 | `str.maketrans` tables (hot-path complement). |
| `_rc` | 697 | IUPAC-aware reverse complement. |
| `_feat_len`, `_slice_circular`, `_bp_in` | 701 / 707 / — | Wrap-aware geometry. |
| `_scan_restriction_sites` | 749 | Palindrome- and wrap-aware. Returns `(resites, recuts)`. |
| `_build_seq_inputs` / `_build_seq_text` | 1220 / 1253 | Sequence-panel renderer, memoized via `_BUILD_SEQ_CACHE`. |
| `_translate_cds` | 1423 | Forward + reverse CDS → protein. Cross-validated against Biopython. |
| `fetch_genbank` / `load_genbank` | 1545 / 1587 | NCBI Entrez fetch + local `.gb`/`.dna` load. |
| `_record_to_gb_text` / `_gb_text_to_record` | 1634 / 1654 | Round-trip SeqRecords as GenBank text. Caller's record never mutated. |
| `_run_plannotate`, `_merge_plannotate_features` | 1719 / 1819 | pLannotate subprocess + merge. |
| `_pick_binding_region` | 3936 | Primer3-compatible region selection. |
| `_design_*_primers` | 3971+ | Detection / cloning / Golden Braid / generic primer design. |
| `_codon_*` | 4212+ | Codon registry, harmonization, NCBI taxid search via `_safe_xml_parse`. |
| `_mut_*` | 4965+ | SOE-PCR mutagenesis primers. |
| `_rebuild_record_with_edit` | in `PlasmidApp` | Edit pipeline preserving wrap features. Sacred invariant #9. |
| `_autosave_*` / `_stash_current_undo_and_load` | in `PlasmidApp` | Crash-recovery autosave + per-plasmid undo stack stashing. |

## pLannotate integration

Shift+A (or ◈ in library panel, or `Edit > Annotate with pLannotate`) runs pLannotate as a subprocess and merges results.

- **Subprocess only, never import.** pLannotate is GPL-3. Never `import plannotate`.
- **Optional runtime dependency.** UI shows install hint when missing.
- **Size cap preflighted** at 50 kb (matches pLannotate's `MAX_PLAS_SIZE`).
- **Merge, don't replace.** Existing features preserved; pLannotate hits get `note="pLannotate"`. Hits matching `(type, start, end, strand)` are skipped.
- **Background worker** with stale-record guard (`self._current_record is captured_record`).
- **Re-entry guard** via `_plannotate_running` flag with `finally` cleanup.
- **Undo-able + dirty-flagged.** `_push_undo()` before merge; `_mark_dirty()` after.

Failures (`PlannotateNotInstalled`, `PlannotateMissingDb`, `PlannotateTooLarge`, `PlannotateFailed`) map to user notifications. Tracebacks → log file only.

## Feature library workbench

Clicking **Features** in the menu bar pushes `FeatureLibraryScreen` (no dropdown). It's the sole place to browse, rename, recolor, or delete persistent feature entries. Per-plasmid feature *enumeration* stays on `FeatureSidebar`.

Entries carry optional `color` (`#RRGGBB`; `None` falls through to type default) and `strand` (`1`/`-1`/`0`/`2` = forward/reverse/arrowless/double-headed; Cycle-Strand walks `1 → -1 → 0 → 2 → 1`). Arrowless suits `rep_origin`, `misc_feature`, stem-loops; double-headed suits inverted repeats.

Snippet preview (`_FeatureSnippetPanel`) feeds a synthesized full-span feature dict through `_build_seq_text` — the same renderer the main `SequencePanel` uses, so previews match post-insertion display. `_render_feature_row_pair` branches on strand: `0`→solid `▒`, `2`→`◀▒…▒▶`, `≥1`→`▒…▒▶`, else→`◀▒…▒`.

Color precedence (`_resolve_feature_color`): entry's `color` → user default in `feature_colors.json` → `_DEFAULT_TYPE_COLORS[type]` → `_FEATURE_PALETTE[0]`. Always returns non-empty so Rich never barfs.

`Add Feature…` and `Annotate with pLannotate` live under the **Edit** menu. Keybindings: `Shift+A` pLannotate, `Ctrl+F` Add Feature, `Ctrl+Shift+F` capture flow.

**Ctrl+Shift+F capture.** Grabs Shift+drag selection (`sp._user_sel`, priority 1) or highlighted feature (`pm.selected_idx`, priority 2), opens `AddFeatureModal` prefilled with slice/name/type/strand/color/qualifiers. **If a drag selection's `(start, end)` matches a feature exactly, capture inherits that feature's full metadata** via `_prefill_from_feature`. Palette colors (`color(N)`) normalised to hex. Insert-at-cursor disabled (bases already in record). Restriction-site overlays (`type == "resite"`) rejected. Save → `_persist_feature_entry` then push `FeatureLibraryScreen`.

`AddFeatureModal` Orientation row: four radios (`#addfeat-strand-fwd/rev/none/both`) → strand `1 / -1 / 0 / 2`. Color row (`#addfeat-color-swatch` + Pick Color… / Auto buttons). `_gather` and `_apply_prefill` round-trip together.

`ColorPickerModal` carries the full xterm 256-color grid (16 ANSI + 216-cube + 24 grayscale) + free-form custom input (accepts `#RGB`, `#RRGGBB`, `0..255`, `color(N)`). `_normalise_color_input` canonicalises to uppercase `#RRGGBB`; `_xterm_index_to_hex` uses canonical `(0, 95, 135, 175, 215, 255)` cube ramp + `8 + 10*k` grayscale. Capability warning surfaces `console.color_system`. `_markup_safe_color` converts stray `color(N)` to hex before render. **Drag-to-preview:** `on_mouse_down` arms `_drag_active` on a `colorpick-x-*` cell (hit-test via `get_widget_at`), `on_mouse_move` repaints the big `#colorpick-preview-swatch`, `on_mouse_up` disarms. Non-left buttons + non-grid mouse-downs ignored.

## Parts Bin source picker

`PartsBinModal` "New Part" opens `DomesticatorModal` with four sources via top-of-modal `RadioSet` (`#dom-src`, `layout: horizontal` + `width: 1fr` + `overflow: hidden` so all four fit on one row). Modal sized `width: 110; max-width: 95%; min-width: 80`.

1. **Direct input** — free-form `TextArea`. `_resolve_source` strips non-IUPAC.
2. **Feature library** — dropdown from `_load_features()`; uses stored `sequence`.
3. **Feature from plasmid** — defaults to current plasmid; `Change…` pushes `PlasmidPickerModal`. On select, `_gb_text_to_record` + `_feats_for_domesticator(rec)` repopulates feature `Select` via `set_options(...)`.
4. **Open FASTA** — `Browse…` pushes `FastaFilePickerModal` (`_FastaAwareDirectoryTree` paints `.fa/.fasta/.fna/.ffn/.frn/.fas/.mpfa/.faa` lime green `#BFFF00`, others white). Parses via `_parse_fasta_single(path)` — validates IUPAC + **rejects multi-record FASTAs**. Errors notify with severity="error".

Panel visibility via `widget.display` toggle. `_feats_for_domesticator(record)` flattens compound/wrap features to outer bounds, drops `source`/`resite`/`recut`/zero-width. Keep in sync with `_feats_in_chunk` / `_extract_feature_entries_from_record`.

### Silent-mutation repair of internal BsaI / Esp3I sites

`DomesticatorModal` carries a codon-table picker (`#dom-codon-row`) seeded with `_codon_tables_get("83333")` (E. coli K12; shared with `MutagenizeModal`). `_design_gb_primers(..., codon_raw=None)` accepts a `{codon: (aa, count)}` dict.

- `_gb_find_forbidden_hits(seq)` returns `(enzyme, site, position)` triples on **both strands**, every occurrence (multi-site contamination must surface fully). `_GB_DOMESTICATION_FORBIDDEN = {"BsaI": "GGTCTC", "Esp3I": "CGTCTC"}`.
- For coding part types (`_GB_CODING_PART_TYPES = CDS / CDS-NS / C-tag`) with truthy `codon_raw` and `len(insert) % 3 == 0`, `_codon_fix_sites` swaps synonymous codons. Reuses the MutagenizeModal harmonizer.
- `_codon_fix_sites` cross-checks before/after hit-sets via `_forbidden_hit_set` — rejects any swap that introduces a new forbidden pattern (no cascade where fixing BsaI spawns Esp3I).
- Partial repair returns an error dict with the partial `mutations` list.
- Non-coding parts, out-of-frame inserts, missing codon table → reject with explanatory reason.

Why both BsaI and Esp3I are forbidden at L0: Esp3I self-cuts during L0 domestication; surviving BsaI re-cuts at L1 assembly. Both must be clean.

The mutated `insert_seq` is what to order as a **gBlock** — primers only change amplicon ends. **Binding-region advisory:** when a mutation lands in the first 18–25 bp (forward binding) or last 18–25 bp (reverse binding), `binding_region_mutations` flags it so the user knows the original plasmid CANNOT be PCR template — they must order the mutated insert and PCR from that.

### Primer naming + pairs list

`_design_gb_primers` returns a **`pairs`** list (1 entry currently; extensibility hook for future SOE-PCR splitting on un-repairable internal sites). Top-level keys mirror `pairs[0]` for back-compat.

`DomesticatorModal`'s **Save Primers** persists each pair to `primers.json` via `_save_primers`:

| Role | Suffix | Example |
|---|---|---|
| Detection (diagnostic PCR) | DET | `myGene-DET-F` / `myGene-DET-R` |
| Cloning (RE tails + GCGC pad) | CLO | `myGene-CLO-F` / `myGene-CLO-R` |
| Golden Braid L0 Domestication | DOM | `myPart-DOM-1-F` / `myPart-DOM-1-R` |

Only domestication primers carry `#` pair number. Dup-sequence guard: existing primer skipped (user notified); other entries in the batch still save. `PrimerDesignScreen` uses the same suffix table.

## On-disk JSON format (schema v1)

All six libraries (`library.json`, `parts_bin.json`, `primers.json`, `codon_tables.json`, `features.json`, `feature_colors.json`) use:

```json
{"_schema_version": 1, "entries": [...]}
```

**Legacy compatibility.** Pre-0.3.1 wrote bare lists; `_extract_entries` accepts both, silently rewrites as envelope on next save. When bumping `_CURRENT_SCHEMA_VERSION`, teach `_extract_entries` to migrate forward in the loader. Newer-version files load with a warning so users know fields may drop on save.

## Crash-recovery autosave

Dirty edits trigger a 3 s debounced write to `_CRASH_RECOVERY_DIR/{safe_id}.gb` (default `~/.local/share/splicecraft/crash_recovery/`). Deleted on `_mark_clean` or abandon. `_check_crash_recovery()` at startup notifies on survivors; user recovers via File > Open.

- `_autosave_path(record)` sanitises `record.id` with `re.sub(r'[^A-Za-z0-9._-]', '_', ...)`, caps at 80 chars.
- Atomic write — `tempfile.mkstemp` + `os.replace`.
- Best-effort — `except Exception: _log.exception(...)`. Autosave is a safety net, not source of truth.
- Debounced via `self.set_timer`. `_mark_dirty` restarts countdown; `_mark_clean` cancels by deleting target.

## Per-plasmid undo/redo stashes

`_apply_record(clear_undo=True)` (switch-plasmid path) stashes outgoing stacks under `record.id` in `_stashed_undo_stacks` / `_stashed_redo_stacks`, restores incoming history if previously edited. LRU-capped at `_MAX_PLASMIDS_WITH_UNDO = 10`. `_current_undo_key` tracks the live stack. `clear_undo=False` (in-place edits — pLannotate, primer-add) leaves stacks intact.

## FASTA export (Parts Bin + Feature Library)

Both screens carry **Export FASTA…** alongside CRUD. Routes through `_export_fasta_to_path(name, sequence, path) -> dict` (atomic; parent dirs created). User sees `FastaExportModal` (mirrors `ExportGenBankModal`). Returns `{"path", "bp", "name"}`. Empty-sequence entries warn instead of opening the modal.

## Parts Bin sequence view + cloning simulator

`PartsBinModal` carries a read-only `TextArea` (`#parts-seq-view`) for the highlighted insert. Click selects all (Ctrl+C-ready). Built-in catalog rows show a placeholder.

Three Copy buttons (all via `_copy_to_clipboard_osc52`):

| Button | Sequence |
|---|---|
| Copy Raw Sequence | `sequence` — insert only, no tails |
| Copy Primed Sequence | `_simulate_primed_amplicon(insert, oh5, oh3)` — `pad + Esp3I + spacer + oh5 + insert + oh3 + rc(spacer+Esp3I+pad)` |
| Copy Cloned Sequence | `_simulate_cloned_plasmid(insert, oh5, oh3)` — `oh5 + insert + oh3 + _PUPD2_BACKBONE_STUB` |

Cloning simulator math sits next to `_GB_L0_ENZYME_SITE` / `_GB_SPACER` / `_GB_PAD`. Golden Braid uses Esp3I/BsmBI at L0, BsaI at L1+ — same N(1)/N(5) geometry, same simulator math. `_PUPD2_BACKBONE_STUB` is a deterministic 420-bp ACGT placeholder scrubbed of every Type IIS site (`GGTCTC`, `GAGACC`, `CGTCTC`, `GAGACG`) on both strands. Replace with licensed pUPD2 and no callers change.

`DomesticatorModal._save` persists `primed_seq` and `cloned_seq` on the part dict; Parts Bin buttons prefer stored values, fall back to simulator at read time for legacy parts.

## Test suite

```bash
python3 -m pytest -n auto -q                          # full, parallel (~2 min)
python3 -m pytest -q                                  # serial (~7 min) — debugging
python3 -m pytest tests/test_dna_sanity.py            # biology only (< 1 s)
python3 -m pytest tests/test_invariants_hypothesis.py # property-based fuzzing
python3 -m pytest -k "palindrome"                     # filter
python3 -m pytest -x                                  # stop on first failure
```

Parallel runs rely on `pytest-xdist` + the autouse `_protect_user_data` fixture (per-test `tmp_path` isolation; monkeypatches `_LIBRARY_FILE`, `_PARTS_BIN_FILE`, `_PRIMERS_FILE`, `_CODON_TABLES_FILE`, `_FEATURES_FILE`, `_FEATURE_COLORS_FILE`, `_CRASH_RECOVERY_DIR` and caches). **No test can write to real user files.** Module-level read-only caches (`_BUILD_SEQ_CACHE`, `_PATTERN_CACHE`, `_SCAN_CATALOG`) are safe — nothing writes them at test time.

`pyproject.toml` sets `asyncio_mode = "auto"` so async tests don't need `@pytest.mark.asyncio`. `tests/conftest.py` provides `tiny_record` / `tiny_gb_path` / `isolated_library` fixtures.

| File | Tests | Covers |
|------|------:|--------|
| `test_dna_sanity.py` | 74 | Sacred invariants 1–6; Type IIS cut-outside-recognition; `_translate_cds` |
| `test_primers.py` | 60 | Detection / cloning / Golden Braid / generic; wrap-region template rotation |
| `test_genbank_io.py` | 68 | `load_genbank` round-trip (GenBank + SnapGene `.dna`); JSON corruption recovery; `_export_fasta_to_path` |
| `test_smoke.py` | 52 | Textual mounts; rotation / view / RE toggles; pLannotate UI + re-entry guard; per-plasmid undo stashes; crash-recovery autosave |
| `test_mutagenize.py` | 49 | SOE-PCR primers, codon substitution, CAI round-trips |
| `test_codon.py` | 42 | Codon registry, harmonization, Kazusa parser, NCBI XML safety, CAI/GC math |
| `test_domesticator.py` | 193 | Golden Braid L0 positions; 4-source picker; `_feats_for_domesticator`; FASTA picker; cloning simulator; codon-fix repair (multi-site, cascade-prevention, binding-region advisory); Save Primers (`pairs` list, DOM suffix) |
| `test_circular_math.py` | 38 | Sacred invariant #5; `_bp_in` / `_feat_len` |
| `test_data_safety.py` | 45 | Sacred invariant #7; envelope round-trip + legacy back-compat + future-version warning; `_atomic_write_text`; `_do_save` atomicity |
| `test_add_feature.py` | 24 | AddFeatureModal: qualifier round-trip, validation, save-to-library dedup, insert-at-cursor |
| `test_plannotate.py` | 24 | Availability, size cap, feature merging, error paths (no real subprocess) |
| `test_modal_boundaries.py` | 26 | Every modal fits in 160×48 (and AddFeatureModal at 100×30) |
| `test_feature_library_screen.py` | 86 | Workbench CRUD + 4-step strand cycle; AddFeatureModal Orientation + Color; Ctrl+Shift+F capture (drag-matches-feature enrichment); ColorPickerModal xterm grid + drag-to-preview; Export-FASTA |
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

Profiled but **not touched**: Textual compositor, Rich `Text.append`, import time.

## Release + versioning

Versions in `pyproject.toml` and `splicecraft.py::__version__`; `release.sh` updates both. See `git log --oneline` for full release history.

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
7. **`_source_path` survives in-place edits.** Cleared only when `clear_undo=True` (fresh loads). Otherwise Ctrl+S after pLannotate or primer-add forgets the original file.
8. **NCBI responses go through `_safe_xml_parse`.** It rejects DOCTYPE/ENTITY before `ET.fromstring`. Don't add a new NCBI endpoint without it.

## How to extend — modular recipes

### A. New pure helper

Place in nearest section per Top-level structure. Snake-case with leading underscore. Pure: no globals, no logging, no UI. Add a test (cross-validate against Biopython where biological). Hot-path? Add a `test_performance.py` budget.

### B. New persisted JSON library

Define `_MYTHING_FILE = _USER_DATA_DIR / "mything.json"`. Route load/save through `_safe_load_json` / `_safe_save_json` — never bypass (invariant #7). Filter `isinstance(entry, dict)` after load. Add the file + its cache to `_protect_user_data` in `tests/conftest.py`. Cover corruption recovery in `test_data_safety.py`.

### C. New modal

Subclass `ModalScreen[ReturnType]` (templates: `FetchModal`, `OpenFileModal`). `query_one("#widget-id", WidgetType)` reads, `self.dismiss(result)` returns. Push via `self.push_screen(MyModal(args), callback=on_result)`. Cover happy path in `test_smoke.py`.

### D. New background worker

`@work(thread=True)` on `PlasmidApp` or owning modal. Wrap body in `try / except Exception as exc`, `_log.exception`, push friendly message via `call_from_thread`. **Stale-record guard:** capture `self._current_record` identity at entry; `if self._current_record is captured_record:` in callback. **Re-entry guard:** `self._myop_running` flag with `finally`. Template: `PlasmidApp._run_plannotate_worker`.

### E. New menu action / keybinding

Add `action_my_thing(self)` on `PlasmidApp`. Add `Binding("key", "my_thing", "desc")` to `BINDINGS`. Add menu item to `MenuBar.compose()`. Modal → recipe C; worker → recipe D.

### F. New full-screen workbench

Subclass `Screen` (or `ModalScreen` if dismissable). Push from menu action with state passed in. Reuse main-app widgets. Register screen-scoped `BINDINGS`.

## Sister project (ScriptoScope)

`/home/seb/proteoscope/scriptoscope.py` (~8,600 lines) — same author, same single-file convention. Patterns to crib if SpliceCraft scales:

- Thread-local `Console` for `_text_to_content` (if seq-panel render blows the 33 ms/frame budget).
- Two-level render cache (`_seq_render_cache` + `_content_cache`, LRU via `OrderedDict.move_to_end`).
- `@lru_cache(1)` availability probes for optional CLI tools (BLAST, Prodigal beyond pLannotate).

User is undecided whether to merge SpliceCraft / ScriptoScope / MitoShift / RefHunter / molCalc into one Textual app with modes, or keep them separate. Either is viable — the single-file convention keeps the option open.

## For future agents

1. Read this file first.
2. `python3 -m pytest -n auto -q` before and after any change. Biology subset (`tests/test_dna_sanity.py`) gives a < 1 s inner loop.
3. Logs: `~/.local/share/splicecraft/logs/splicecraft.log` (or `$SPLICECRAFT_LOG`). 8-char session ID per run.
4. Don't break sacred invariants — touching `_scan_restriction_sites`, `_rc`, `_iupac_pattern`, `_translate_cds`, `_bp_in`, `_feat_len`, the midpoint formula, or `_rebuild_record_with_edit` will trip the relevant tests immediately.
5. Follow the error convention: `_log.exception` for stack traces, `notify()` / `Static.update("[red]...[/]")` for the user. Never raw tracebacks.
6. Eyeball real-world behaviour on pUC19 (`L09137`) and pACYC184 (`MW463917.1`).
7. Past fix history lives in git — `git log --oneline` and `git show <hash>`.
