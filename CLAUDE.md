# CLAUDE.md — AI Agent Context for SpliceCraft

This file is the **agent handoff document** for SpliceCraft. Any AI agent can read this file to pick up development without needing the full conversation history.

The project is developed in continuous collaboration between a human bioinformatician and an AI agent (Claude Opus 4.6+).

---

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer, sequence editor, and cloning/mutagenesis workbench** built with Python 3.10+ / Textual / Biopython. Renders Unicode braille-dot plasmid maps directly in the terminal, with a per-base sequence panel, restriction-site overlays, a plasmid library, Golden Braid L0 assembly tooling, Primer3-backed primer design, and SOE-PCR site-directed mutagenesis.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs org, user ATinyGreenCell)

- **Single-file architecture:** the entire app is `splicecraft.py` (~13,200 lines). Intentional — avoids import complexity and keeps the codebase greppable. Sibling project ScriptoScope follows the same convention at ~8,600 lines.
- **Test suite:** 857 tests across 16 files in `tests/` (last refresh 2026-04-21). Sequential run ~400 s; parallel run (`pytest -n auto`) ~145 s on 8 cores (~2.8× speedup). Biology subset (`test_dna_sanity.py`) < 1 s. `test_invariants_hypothesis.py` adds property-based fuzzing on top of hand-written regression tests.
- **Dependencies:** `textual>=8.2.3`, `biopython>=1.87`, `primer3-py>=2.3.0`, `platformdirs>=4.2`, plus `pytest>=9.0` / `pytest-asyncio>=1.3` / `pytest-xdist>=3.6` / `hypothesis>=6.100` for tests. Users install via `pipx install splicecraft`. **Optional runtime:** `pLannotate` (conda, GPL-3) for the Shift+A annotation feature.
- **Published on PyPI** as `splicecraft`. Releases cut via `./release.sh X.Y.Z` (bumps version in both `pyproject.toml` and `splicecraft.py`, runs tests, builds, commits+tags+pushes; GitHub Actions `publish.yml` then publishes via Trusted Publishing / OIDC). Latest published: **v0.3.1**.

## How to run

```bash
cd ~/SpliceCraft
python3 splicecraft.py              # empty canvas
python3 splicecraft.py L09137       # fetch pUC19 from NCBI
python3 splicecraft.py myplasmid.gb # open local GenBank file (.gb/.gbk/.dna)
python3 -m pytest -n auto -q        # full test suite (parallel, ~2 min on 8 cores)

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
| 3868–4076 | Golden Braid L0 position catalog (Esp3I/BsmBI overhangs, position constraints) |
| 4077–4130 | Parts-bin + primer-library persistence |
| 4131–4925 | Codon-usage registry (`_codon_*`), Kazusa parser, NCBI taxid search (`_safe_xml_parse`), harmonization, CAI/GC. Crash-recovery config (`_CRASH_RECOVERY_DIR`) sits at the top of this slab |
| 4926–5437 | SOE-PCR site-directed mutagenesis primer design (`_mut_*`) |
| 5953–6394 | `PlasmidFeaturePickerModal`, `AddFeatureModal` |
| 6395–7162 | **Feature library workbench**: `_normalise_color_input`, `_xterm_index_to_hex`, `ColorPickerModal`, `_FeatureSnippetPanel`, `FeatureLibraryScreen` (full-screen; entered by clicking `Features` in the menu bar) |
| 7163–7345 | `PartsBinModal` |
| 7346–7550 | **FASTA file picker** (`_FASTA_EXTS`, `_is_fasta_path`, `_parse_fasta_single` — rejects multi-record, `_FastaAwareDirectoryTree`, `FastaFilePickerModal` — lime-green FASTA highlight, white otherwise) |
| 7494–8003 | `_feats_for_domesticator` helper + `DomesticatorModal` (4-source part picker: Direct / Feature Library / Feature-from-Plasmid / Open FASTA) |
| 8004–8322 | `ConstructorModal` (Golden Braid L0 assembly UI) |
| 8323–8660 | `NcbiTaxonPickerModal` + `SpeciesPickerModal` (codon-table picker) |
| 8661–8891 | Mutagenize helpers (`_MutPreview`, `AminoAcidPickerModal`) |
| 8892–9536 | `MutagenizeModal` — full mutagenesis workflow |
| 9537–10635 | `PrimerDesignScreen` — full-screen primer workbench |
| 10636–10856 | Small modals (`UnsavedQuitModal`, `PlasmidPickerModal`, `RenamePlasmidModal`, `LibraryDeleteConfirmModal`) |
| 10857–end | `PlasmidApp` — main controller, keybindings, per-plasmid undo/redo stashes, crash-recovery autosave, `@work` threads; `main()` entry point |

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

Shift+A (or ◈ in the library panel, or `Edit > Annotate with pLannotate`) runs pLannotate as a subprocess and merges results into the current record.

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

## Feature library workbench

Clicking **Features** in the top menu bar pushes `FeatureLibraryScreen` directly — no dropdown. The screen is the sole place to browse, rename, recolor, or delete persistent feature entries. Per-plasmid feature *enumeration* remains on the right-hand `FeatureSidebar` (unchanged).

Entries gained two optional fields in v0.3.2:
- `color`: per-entry `"#RRGGBB"` override. `None` → fall through to the type default.
- `strand`: `1` (forward, `▶`) / `-1` (reverse, `◀`) / `0` (arrowless, `▒`) / `2` (double-headed, `◀…▶`). The Cycle-Strand button walks the 4-step loop `1 → -1 → 0 → 2 → 1`. Arrowless is meaningful for `rep_origin`, `misc_feature`, stem-loops; double-headed for inverted repeats / palindromic regulatory elements.

The snippet preview (`_FeatureSnippetPanel`) synthesizes a single full-span feature dict from the selected entry and feeds it through the shared **`_build_seq_text`** pipeline — the exact same renderer the main `SequencePanel` uses. Arrow direction comes from `_render_feature_row_pair`'s strand handling, so the preview always matches what the user sees in the main app after insertion. `_render_feature_row_pair` branches on strand: `0` → solid `▒` bar, `2` → `◀▒…▒▶`, `≥1` → `▒…▒▶`, else → `◀▒…▒`.

Color precedence (implemented in `_resolve_feature_color`): entry's `color` field → user default in `feature_colors.json` → built-in `_DEFAULT_TYPE_COLORS[type]` → `_FEATURE_PALETTE[0]`. Always returns a non-empty string so Rich never barfs.

`feature_colors.json` stores `{feature_type → hex}` as a list of `{"feature_type": ..., "color": ...}` dicts under the standard schema envelope. Missing / empty / corrupt → `{}`, and callers degrade to the built-in defaults.

`Add Feature…` and `Annotate with pLannotate` moved to the **Edit** menu when the Features dropdown was eliminated. Keybindings (`Shift+A` for pLannotate, `Ctrl+F` for Add Feature, `Ctrl+Shift+F` for the capture flow) live there.

**Ctrl+Shift+F capture shortcut.** From the main view, `Ctrl+Shift+F` invokes `PlasmidApp.action_capture_to_features`, which grabs either the Shift+drag DNA selection (`sp._user_sel`, priority 1) or the highlighted feature (`pm.selected_idx`, priority 2) and opens `AddFeatureModal` prefilled with the slice/name/type/strand/color/qualifiers. **If the drag selection's `(start, end)` exactly matches a feature's range, the capture inherits that feature's full metadata** (type, strand, color, qualifiers) via `_prefill_from_feature` — sidebar-click (which sets both `_user_sel` and `selected_idx`) and drag-selecting a feature both produce the same rich prefill. Palette-style colors (`color(N)`) are normalised to hex at capture so stored library entries and Rich markup previews never choke. Insert-at-cursor is disabled in this flow (the bases already live in the record). On Save, the app persists via `_persist_feature_entry` and pushes `FeatureLibraryScreen` so the user lands in the workbench with the new entry visible. Restriction-site overlays (`type == "resite"`) are rejected with a notification — they aren't real features. The shared helper `_persist_feature_entry` is also used by the regular Add-Feature path.

`AddFeatureModal`'s direction row is labelled **Orientation** (not "Strand") and holds four radios (`#addfeat-strand-fwd/rev/none/both`) backed by strand values `1 / -1 / 0 / 2`. The modal also carries a **Color** row (`#addfeat-color-swatch` + `Pick Color…` / `Auto` buttons) so captured or manually-set colors survive through Save. `_gather` and `_apply_prefill` round-trip strand, color, and the standard fields together.

**Expanded color picker.** `ColorPickerModal` grew a full xterm 256-color grid (16 ANSI + 216-cell cube + 24 grayscale) and a free-form custom input that accepts any of `#RGB`, `#RRGGBB`, `0..255` (xterm index), or `color(N)`. The helper `_normalise_color_input` canonicalises these to uppercase `#RRGGBB`, and `_xterm_index_to_hex` converts indices via the canonical `(0, 95, 135, 175, 215, 255)` cube ramp + `8 + 10*k` grayscale. A capability warning (`#colorpick-capability`) surfaces the terminal's `console.color_system` so users on 8/16-color terminals know truecolor picks will be approximated. `_markup_safe_color` (paired with `_resolve_feature_color`) converts any stray `color(N)` palette value to hex before rendering so Rich's markup lexer never trips on the parens. A large preview swatch (`#colorpick-preview-swatch`, 24w × 5h with tall border) dominates the header so the picked color is obvious at a glance, and **drag-to-preview** works across the xterm grid: `on_mouse_down` arms `self._drag_active` if the click lands on a `colorpick-x-*` cell (hit-tested via `get_widget_at`), `on_mouse_move` repaints the swatch every time the cursor enters a new cell, and `on_mouse_up` disarms. Non-left buttons and non-grid mouse-downs are ignored so the regular Save/Cancel/Apply buttons still work normally.

## Parts Bin source picker

The **New Part** flow in `PartsBinModal` opens `DomesticatorModal`, which now offers four mutually-exclusive sources selected via a top-of-modal `RadioSet` (`#dom-src`). The RadioSet uses `layout: horizontal` + `width: 1fr` per button + `overflow: hidden` so all four radios fit on **one row** with no scrollbar; the whole modal uses `width: 110; max-width: 95%; min-width: 80` so it flexes on narrow terminals.

1. **Direct input** (`#dom-panel-direct`) — free-form `TextArea` paste. `_resolve_source` strips anything outside IUPAC before handing a `(cleaned, 0, len(cleaned))` tuple to the Golden Braid validator.
2. **Feature library** (`#dom-panel-featlib`) — dropdown populated from `_load_features()`. Selection pulls the stored `sequence` and fires the validator as the whole span (`0 → len(seq)`).
3. **Feature from plasmid** (`#dom-panel-plasmid`) — defaults to the currently-open plasmid (name threaded in as `current_plasmid_name`) but the `Change…` button pushes `PlasmidPickerModal(current_id=...)`. On selection the modal calls `_gb_text_to_record` + `_feats_for_domesticator(rec)` and repopulates the feature `Select` via `sel.set_options(...)`. The picked feature's `start`/`end` are reused verbatim so the Golden Braid coordinate checks run on the real plasmid slice.
4. **Open FASTA** (`#dom-panel-fasta`) — `Browse…` button pushes `FastaFilePickerModal`, which renders a `DirectoryTree` subclass (`_FastaAwareDirectoryTree`) that paints FASTA files (`.fa / .fasta / .fna / .ffn / .frn / .fas / .mpfa / .faa`, case-insensitive) in **bold lime green** (`#BFFF00`) and every other file in **white** (`#FFFFFF`). On Open the modal dismisses with an absolute path; the domesticator parses it via `_parse_fasta_single(path)` (validates IUPAC + **rejects multi-record FASTAs** — we refuse to silently guess which record the user wanted) and stashes the sequence on `self._fasta_seq`. `_resolve_source` then returns `(self._fasta_seq, 0, len(self._fasta_seq))`. Parse errors surface as `app.notify(..., severity="error")` so the user sees "Multi-sequence FASTA not supported (N records found)", "Non-IUPAC characters…", or "No FASTA records found" rather than a traceback.

Panel visibility is handled by toggling `widget.display` on the four `Vertical` panels — there is still exactly one Save button and one Design button. The helper `_feats_for_domesticator(record)` (module-level, just above `DomesticatorModal`) flattens compound/wrap features to their outer bounds, drops `source`/`resite`/`recut`/zero-width entries, and returns the `{label, type, start, end, strand}` dicts the Golden Braid UI already speaks. Keep that helper in sync with `_feats_in_chunk` / `_extract_feature_entries_from_record` — the three are the canonical ways to translate `SeqFeature`s into dict-shape for UI consumers.

### Silent-mutation repair of internal BsaI / Esp3I sites

`DomesticatorModal` carries a codon-table picker (`#dom-codon-row` → `#dom-codon-label` + `#btn-dom-codon`) mirroring the one in `MutagenizeModal`. `on_mount` seeds `self._codon_entry` via `_codon_tables_get("83333")` (E. coli K12, the shared registry default) and `Change…` pushes `SpeciesPickerModal` with `_codon_picked` as the callback — the picker is shared across both modals so the user can add / reuse codon tables once.

`_design_gb_primers(..., codon_raw=None)` accepts an optional `{codon: (aa, count)}` dict. When the insert contains an internal **BsaI (GGTCTC / GAGACC)** or **Esp3I (CGTCTC / GAGACG)** site:
- The module-level helper `_gb_find_forbidden_hits(seq)` returns `(enzyme_name, site_found, position)` triples on both strands (module-level constant `_GB_DOMESTICATION_FORBIDDEN = {"BsaI": "GGTCTC", "Esp3I": "CGTCTC"}`).
- If `part_type ∈ _GB_CODING_PART_TYPES` (CDS / CDS-NS / C-tag) **and** `codon_raw` is truthy **and** `len(insert) % 3 == 0`, the function translates the insert and calls `_codon_fix_sites(insert, protein, codon_raw, sites=_GB_DOMESTICATION_FORBIDDEN)` — reusing the exact helper the MutagenizeModal's harmonizer already uses. On full repair, `insert_seq` in the result is the *mutated* sequence and `mutations` is a list of `"BsaI at nt N: GGT→GGC (codon C aa A, freq=F)"` strings.
- If the fix is partial (leftover sites overlap codons with no synonymous swap), the function returns an error dict **with the partial `mutations`** so the user sees what was repaired before giving up.
- Non-coding parts, out-of-frame inserts, and calls without a codon table still reject with an explanatory reason — these can't be fixed synonymously, so the user must pick a different template region or redesign manually.

Why both BsaI and Esp3I are forbidden at L0: Esp3I self-cuts during L0 domestication, but a surviving BsaI site would re-cut during the downstream L1 assembly. Both must be clean for the part to round-trip through Golden Braid cleanly.

The mutated `insert_seq` is what the user should order as a **gBlock / synthetic fragment** — primers only change the amplicon ends, not the middle, so silent mutations inside the insert only apply if the PCR template already carries them.

**Hardening (2026-04-21).** Three failure modes were tightened because "silent pass with a site still in the final insert" = wasted synthesis budget:

1. **All occurrences reported.** `_gb_find_forbidden_hits` walks every match (not just the first per enzyme), so multi-site contamination surfaces completely in error messages and in the pre-fix scan.
2. **No swap cascades.** `_codon_fix_sites` computes `_forbidden_hit_set(seq, all_forbidden)` before each candidate swap and rejects any swap whose after-set contains an entry not present in the before-set. A BsaI fix can never silently spawn an Esp3I (or the RC of either) anywhere in the sequence.
3. **Binding-region advisory.** When a mutation lands inside the first 18–25 bp (forward binding) or last 18–25 bp (reverse binding) of the insert, `_design_gb_primers` populates `binding_region_mutations` — a list of `{text, region: "fwd"|"rev", codon_start}` dicts. The Domesticator result panel surfaces this in red so the user knows the original plasmid CANNOT be used as PCR template; they must order the mutated insert as a gBlock and PCR from that.

### Primer naming + pairs list

`_design_gb_primers` returns a dict with a **`pairs`** list — currently exactly one entry, shaped like the top-level result (legacy callers that read `result["fwd_full"]` directly still work since the top-level keys mirror `pairs[0]`). The list is the extensibility hook for a future SOE-PCR splitting path: when an internal Type IIS site can't be silently repaired, the insert will be split at the bad site and each sub-amplicon will contribute its own pair.

`DomesticatorModal`'s **Save Primers** button writes every designed pair to `primers.json` via `_save_primers`, using the project-wide naming convention:

| Primer role | Suffix | Example |
|---|---|---|
| Detection (diagnostic PCR) | **DET** | `myGene-DET-F` / `myGene-DET-R` |
| Cloning (RE tails + GCGC pad) | **CLO** | `myGene-CLO-F` / `myGene-CLO-R` |
| Golden Braid L0 Domestication | **DOM** | `myPart-DOM-1-F` / `myPart-DOM-1-R` (pair 1), `-DOM-2-F/R` (pair 2), … |

Only domestication primers carry the `#` pair number, since Detection and Cloning always ship as a single pair. Dup-sequence guard: if a primer's sequence already exists in the library, that one entry is skipped (and the user is notified) — the other entries in the batch still save. `PrimerDesignScreen` uses the same suffix table for its auto-fill of the Save-Primer name inputs.

## On-disk JSON format (schema v1)

All six persisted libraries (`library.json`, `parts_bin.json`, `primers.json`, `codon_tables.json`, `features.json`, `feature_colors.json`) use the envelope shape:

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
python3 -m pytest -n auto -q                        # full suite, parallel (~2 min on 8 cores)
python3 -m pytest -q                                # full suite, serial (~7 min) — use when debugging
python3 -m pytest tests/test_dna_sanity.py          # only biology (< 1 s)
python3 -m pytest tests/test_invariants_hypothesis.py  # property-based fuzzing
python3 -m pytest -k "palindrome"                   # filter by name
python3 -m pytest -x                                # stop on first failure (implies serial)
```

**Parallel runs** (`-n auto`) rely on `pytest-xdist` and the autouse
`_protect_user_data` fixture's per-test `tmp_path` isolation. Workers share
the module-level read-only caches (`_BUILD_SEQ_CACHE`, `_PATTERN_CACHE`,
`_SCAN_CATALOG`) — nothing writes to them at test time. Use serial mode
(`-x`, `--pdb`, `-s`) when you need ordered output or debugger attach.

`pyproject.toml` sets `asyncio_mode = "auto"` so async tests don't need `@pytest.mark.asyncio`. `tests/conftest.py` defines `tiny_record` / `tiny_gb_path` / `isolated_library` fixtures, and installs the **autouse** `_protect_user_data` fixture that monkeypatches `_LIBRARY_FILE`, `_PARTS_BIN_FILE`, `_PRIMERS_FILE`, `_CODON_TABLES_FILE`, `_FEATURES_FILE`, `_FEATURE_COLORS_FILE`, `_CRASH_RECOVERY_DIR`, and their caches to tmp paths. **No test can write to real user files.**

### Files

| File | Tests | Covers |
|------|------:|--------|
| `test_dna_sanity.py` | 74 | Sacred invariants 1–6; Type IIS cut-outside-recognition; `_translate_cds` forward & reverse |
| `test_primers.py` | 60 | Detection / cloning / Golden Braid / generic; **wrap-region primer design** (template rotation, modular position mapping) |
| `test_genbank_io.py` | 68 | `load_genbank` round-trip (GenBank + SnapGene `.dna`); `_save_library` / `_load_library` JSON round-trip + corruption recovery; **`_export_fasta_to_path` atomic-write round-trip + empty-name / empty-seq rejection + overwrite + `.tmp` cleanup** |
| `test_smoke.py` | 52 | Textual app mounts; panels present; rotation / view-toggle / RE-toggle; pLannotate UI + re-entry guard; `_apply_record` semantics; sidebar wrap-coord display; undo snapshot independence; **per-plasmid undo stashes + LRU eviction**; **crash-recovery autosave** |
| `test_mutagenize.py` | 49 | SOE-PCR primer design, codon substitution, `_mut_revcomp` / translate / CAI round-trips |
| `test_codon.py` | 42 | Codon registry persistence, harmonization, Kazusa parser, NCBI taxid XML safety, CAI/GC math |
| `test_domesticator.py` | 193 | Golden Braid L0 positions / overhangs, part validation, assembly lanes; **Parts Bin 4-source picker (2026-04-20)**: RadioSet layout + `display`-based panel swap, Direct Input cleaning / Feature Library lookup / Feature-from-Plasmid lookup / **Open-FASTA** paths in `_resolve_source`, `_feats_for_domesticator` flattens compound wraps and drops `source`/resite/recut/zero-width, `PlasmidPickerModal` swap refreshes the feature Select via `set_options`; **horizontal-radio layout regression guards** (same y-coord, no vertical scrollbar, modal fits 90-col terminal); **FASTA picker** (`_is_fasta_path` extension matrix, `_parse_fasta_single` happy / **multi-record reject** / error paths, `_FastaAwareDirectoryTree` paints `.fa/.fasta/.fna/...` lime green and others white); **Parts Bin Export-FASTA button** (present, pushes `FastaExportModal` for user part, warns for built-in catalog rows); **Cloning simulator** (`_PUPD2_BACKBONE_STUB` deterministic + free of BsaI/BsmBI sites on both strands, `_simulate_primed_amplicon` digest-carves back to `oh5+insert+oh3`, `_simulate_cloned_plasmid` yields `oh5+insert+oh3+backbone`); **Parts Bin sequence TextArea + 3 Copy buttons** (Raw / Primed / Cloned via OSC 52, fallback when `primed_seq`/`cloned_seq` missing on legacy parts, warn on built-in catalog rows); **Silent-mutation repair of internal BsaI / Esp3I sites** (coding CDS/CDS-NS/C-tag route through `_codon_fix_sites` via shared codon registry; non-coding parts + out-of-frame inserts + missing codon table still reject; `DomesticatorModal` codon-picker UI defaults to E. coli K12, threads `self._codon_entry['raw']` into `_design_gb_primers`); **Save Primers to Library (2026-04-21)**: `_design_gb_primers` now returns a `pairs` list (1 entry currently, extensible for future SOE-PCR splitting) with top-level keys mirroring `pairs[0]` for back-compat; `DomesticatorModal` "Save Primers" button persists each pair as `{partName}-DOM-{n}-F/R` via `_save_primers` with dup-sequence guard; `PrimerDesignScreen` goldenbraid mode now uses DOM suffix (vs CLO for cloning / DET for detection); **Multi-site / cascade hardening (2026-04-21)**: `_gb_find_forbidden_hits` reports EVERY occurrence (not just first per enzyme) so multi-site contamination can't slip past the error path; `_codon_fix_sites` swap loop rejects any candidate that would introduce a new forbidden pattern anywhere (before/after hit-set cross-check via `_forbidden_hit_set`), preventing cascade failures where fixing BsaI accidentally spawns Esp3I; `binding_region_mutations` flags silent mutations that land inside the 5′ or 3′ primer binding windows so the user knows they must order the mutated insert as a gBlock and cannot PCR from the original template; edge-case coverage includes multi-BsaI, multi-Esp3I, mixed enzymes, reverse-strand (GAGACC / GAGACG), sites at 5′ / 3′ / interior, and cascade-prevention probes |
| `test_circular_math.py` | 38 | Sacred invariant #5 (wrap midpoint); `_bp_in` / `_feat_len` for wrapped / non-wrapped / zero-width |
| `test_data_safety.py` | 37 | Sacred invariant #7 (atomic saves, `.bak` recovery); **schema-envelope round-trip + legacy bare-list back-compat + future-version warning + shrink-guard counting both formats**; `features.json` redirected by `_protect_user_data`; `_protect_user_data` fixture confirmation |
| `test_add_feature.py` | 24 | **AddFeatureModal + insert pipeline**: qualifier parsing round-trip, `_extract_feature_entries_from_record` strand/wrap handling, modal form validation (empty name / invalid bases / IUPAC), save-to-library dedup, insert-at-cursor (fwd / rev / coord shift / dirty flag) |
| `test_plannotate.py` | 24 | Availability detection, size-cap preflight, feature merging, subprocess error paths (subprocess never actually invoked) |
| `test_modal_boundaries.py` | 26 | **Every modal stays inside the terminal**: root-container bounds + non-scrollable descendants fit at the baseline 160×48 (covers ColorPickerModal + FastaFilePickerModal + **FastaExportModal**); AddFeatureModal-specific regression guards at 160×48 / 120×40 / 100×30 (regression for 2026-04-20 textbox-offscreen bug) |
| `test_feature_library_screen.py` | 86 | **Features-tab workbench rework (2026-04-20)**: Menu click routes to FeatureLibraryScreen; CRUD actions (add / duplicate / remove / cycle-strand) persist via `_save_features`; **four-step strand cycle (1 → -1 → 0 → 2 → 1)**; ColorPickerModal returns expected dict shape; snippet DNA panel routes through `_build_seq_text` so `▶` / `◀` / `▒` / `◀…▶` arrows reflect strand 1 / -1 / 0 / 2; **AddFeatureModal 4-way Orientation radios** (Forward/Reverse/Arrowless/Double) with prefill+save round-trip; **Ctrl+Shift+F capture** (drag selection or highlighted feature → prefilled modal → Save → FeatureLibraryScreen); **drag-matches-feature enrichment** (exact-range drag inherits type/strand/color/qualifiers); **AddFeatureModal Color field** (prefill round-trip, Auto clears, capture threads color through); `_normalise_color_input` / `_xterm_index_to_hex` parametrized validation; expanded `ColorPickerModal` xterm-cell click + custom hex/index apply + capability warning; **drag-to-preview** (MouseDown arms `_drag_active`, MouseMove across cells repaints the big `#colorpick-preview-swatch` in real time, MouseUp disarms; non-left buttons ignored, non-grid MouseDown never arms); **Export-FASTA button** (present, pushes `FastaExportModal` threaded with selected entry, warns on empty library + empty-sequence entry) |
| `test_features_library.py` | 29 | Persistent feature-library JSON round-trip, schema envelope, corruption recovery, cache invalidation, `_GENBANK_FEATURE_TYPES` curation (CDS / gene / promoter present, `source` excluded); **per-entry `color` field + `strand=0` round-trip**; `_load_feature_colors` / `_save_feature_colors` persistence; `_resolve_feature_color` precedence (entry → user default → `_DEFAULT_TYPE_COLORS` → palette) |
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
2. **Run `python3 -m pytest -n auto -q`** before and after any change. 857 tests, ~145 s on 8 cores (or ~400 s serial). Biology subset (`tests/test_dna_sanity.py`) runs in < 1 s for a fast inner loop.

## FASTA export (Parts Bin + Feature Library)

Both library screens carry an **Export FASTA…** button alongside the usual CRUD buttons. The button routes through `_export_fasta_to_path(name, sequence, path) -> dict` (atomic `tempfile.mkstemp` + `os.replace`; fsync best-effort; parent dirs created). The user sees `FastaExportModal`, which mirrors `ExportGenBankModal` — Input field for path, Export/Cancel buttons, inline error status. On dismiss the caller gets `{"path", "bp", "name"}` and notifies. Entries without a sequence (built-in Golden Braid catalog parts in `_GB_L0_PARTS`, or library entries with empty `sequence`) warn via `app.notify(..., severity="warning")` instead of pushing an empty modal.

## Parts Bin sequence view + cloning simulator

`PartsBinModal` carries a scrollable **read-only TextArea** (`#parts-seq-view`) that holds the full insert of the highlighted row. Clicking anywhere on the TextArea selects every character (`TextArea.select_all`) so single-click → Ctrl+C is enough. Built-in catalog rows (no sequence) show a placeholder message instead of looking empty.

Three **Copy buttons** sit below the TextArea, all routed through `_copy_to_clipboard_osc52(text)`:

| Button | Sequence copied |
|---|---|
| Copy Raw Sequence    | `sequence` — just the insert, no primer tails |
| Copy Primed Sequence | `_simulate_primed_amplicon(insert, oh5, oh3)` — `pad + Esp3I + spacer + oh5 + insert + oh3 + rc(spacer+Esp3I+pad)` |
| Copy Cloned Sequence | `_simulate_cloned_plasmid(insert, oh5, oh3)` — `oh5 + insert + oh3 + _PUPD2_BACKBONE_STUB` (circular, linearised at the 5′ overhang) |

**Cloning simulator math** lives next to the `_GB_L0_ENZYME_SITE` / `_GB_SPACER` / `_GB_PAD` constants. Golden Braid splits enzymes across assembly levels: **L0 parts are domesticated with Esp3I (CGTCTC) / BsmBI**, while L1+ transcriptional units are assembled with BsaI (GGTCTC). The two have identical N(1)/N(5) geometry (→ 4-nt 5′ overhangs), so the same simulator math works for both; the constant just picks the recognition sequence. `_PUPD2_BACKBONE_STUB` is a **deterministic 420-bp ACGT placeholder** (seeded via `_build_pupd2_backbone_stub`) scrubbed of every Type IIS site on both strands — `GGTCTC`, `GAGACC` (BsaI), `CGTCTC`, `GAGACG` (Esp3I / BsmBI) — so the simulated cloned plasmid is guaranteed not to re-cut in either L0 or L1 assembly. Replace `_PUPD2_BACKBONE_STUB` with a licensed real pUPD2 sequence and no callers change.

`DomesticatorModal._save` persists `primed_seq` and `cloned_seq` on the part dict alongside the raw insert and primers; the Parts Bin buttons prefer those stored values but fall back to the simulator at read time for parts saved before the simulator existed.
3. **Check `~/.local/share/splicecraft/logs/splicecraft.log`** (or `$SPLICECRAFT_LOG`) when debugging. Every session has a unique 8-char ID.
4. **Don't break the sacred invariants.** Each has a test (see mapping table). If you touch `_scan_restriction_sites`, `_rc`, `_iupac_pattern`, `_translate_cds`, `_bp_in`, `_feat_len`, the midpoint formula, or `_rebuild_record_with_edit`, the relevant tests will tell you immediately.
5. **Follow the error-handling convention**: `_log.exception` for stack traces, `notify()` or `Static.update("[red]...[/]")` for the user. Narrow `except` types. Never let raw tracebacks hit the TUI.
6. **When in doubt about real-world behavior** — eyeball it on pUC19 (`L09137`) and pACYC184 (`MW463917.1`), both fetched at first-run.
7. **Past fix history lives in git.** Use `git log --oneline` and `git show <hash>` rather than restoring fix-log sections to this file.
