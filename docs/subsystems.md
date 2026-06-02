# Subsystems — feature-area deep dives

Companion to `CLAUDE.md`. Holds the long-form documentation for self-contained subsystems whose details are only relevant when touching that area.

## How agents should use this file

**Before editing one of these subsystems, grep this file for the tag.** Each section carries a `[SUB-xxx]` anchor + topic keywords.

```bash
grep -n '\[SUB-plasmidsaurus\]' docs/subsystems.md
grep -ni 'experiment\|project\|gel' docs/subsystems.md
```

## Tag index (grep table)

| Tag | Topic keywords |
|---|---|
| `[SUB-plasmidsaurus]` | Plasmidsaurus zip ingestion; pairwise alignment; SequencingScreen; sub-tabs |
| `[SUB-experiments]` | Experiments lab-notebook; projects; cross-refs `@`/`!`/`&`; spellcheck; ImageAttachModal |
| `[SUB-gels]` | Gels; gels.json; GelLibraryModal; SimulatorScreen |
| `[SUB-hmm-db]` | HMM database registry; Pfam-A / NCBIfam / custom URL download; `HmmDbCatalogModal`; `hmm_db_catalog.json`; pyhmmer hmmpress; cross-internet hardening |

---

## [SUB-hmm-db] HMM database registry + download manager (sweep #28, 2026-05-25)

Replaces the "paste an absolute path to Pfam-A.hmm" UX with a built-in
catalog + downloader. Lives under the BLAST modal's HMM-database
picker.

**Catalog file** `<DATA_DIR>/hmm_db_catalog.json` (envelope-v1) holds
the registered databases. Two builtins (`pfam-a`, `ncbifam`) are
re-injected on every load if missing, so a hand-edited catalog can't
permanently nuke the well-known defaults. User can add custom URLs
via `HmmDbAddEditModal`; user entries are renameable / removable.

**Per-DB on disk** `<DATA_DIR>/hmm_databases/<id>/`:
* `db.hmm.gz` — the streamed download (kept for re-press if pyhmmer
  changes its on-disk format)
* `db.hmm` — decompressed HMMER3 file
* `db.hmm.h3i` / `.h3m` / `.h3p` / `.h3f` — pyhmmer's `hmmpress` index
  files (created post-download)
* `meta.json` — schema-v1 envelope with `{id, version, downloaded_at,
  sha256, source_url, n_profiles, pressed, last_remote_version,
  last_remote_version_check_ts}`

**Modal stack**
* `BlastModal` HMM-database row: `Select` of registered DBs + Manage
  button + status banner (`✓ ready` / `update available` /
  `not downloaded`). On modal mount + on picker change, schedules a
  24h-cached background version check via `_auto_version_check_worker`.
* `HmmDbCatalogModal`: DataTable of every registered DB; Download,
  Check for updates, Add custom URL, Edit, Remove entry, Delete
  download buttons. Workers `_download_worker` (`@work(thread=True,
  exclusive=True, group="hmm_db_download")`) and `_version_check_
  worker`. `_blocks_undo=True`; `_dismiss_once`.
* `HmmDbAddEditModal`: Inputs for name, source URL, optional version
  URL, description. Builtin-edit mode locks name + builtin flag,
  only URL editable. `_blocks_undo=True`; `_dismiss_once`.

**Sacred network discipline** (see [INV-84] for the full list of layers):
HTTPS-only-by-default; bounded redirects (5 max) via custom
`OpenerDirector`; explicit `ssl.create_default_context()`;
Content-Type guard refuses HTML/JSON; gzip + HMMER3 magic
verification on first chunk; disk-space pre-check; one retry with
250ms backoff; cancel-aware (`cancel_check_cb` polled between
chunks; modal-close ⇒ cancellation); URL credential redaction in
every log line; cross-modal-instance download slot
(`_HMM_DB_DOWNLOAD_INFLIGHT`); atomic hmmpress cleanup; n_profiles
== 0 treated as download failure.

**Caps:**
* `_HMM_DB_DOWNLOAD_MAX_BYTES = 4 GB` — refuses a typo'd URL that
  points at a 100 GB tarball (or a hostile mirror serving /dev/zero).
* `_HMM_DB_VERSION_MAX_BYTES = 64 KB` — version file is tiny (Pfam's
  is ~150 bytes); a malicious mirror serving a huge version file
  can't OOM us.
* `_HMM_DB_DISK_RESERVE_MULTIPLIER = 2.5` — free space must be ≥
  Content-Length × 2.5 (covers decompression + hmmpress headroom).
* `_HMM_DB_NETWORK_TIMEOUT_S = 60` — connect + read timeout per op.
* `_HMM_DB_VERSION_CHECK_TTL_S = 86400` — 24h cache so reopening the
  modal doesn't re-poll on every click.

**Tests:** `tests/test_hmm_db_catalog.py` (57 cases) cover catalog
persistence (builtins always present; user entries round-trip;
corrupted entries dropped), sanitisation (ASCII-only id, no path
traversal/NUL/control chars; URL whitespace rejected pre-strip;
scheme allowlist), URL credential redaction, content-type guard,
disk-space pre-check (stubs `shutil.disk_usage`), magic-byte
verification (rejects non-gzip + non-HMMER3 + zip-bomb + corrupt
gzip), hmmpress cleanup, network retry (one retry + persistent
failure), version-file parsing, per-DB local state with 24h cache,
cross-modal download slot (concurrent acquire only one wins), and
modal smoke (fits 160×48).

---

## [SUB-plasmidsaurus] Plasmidsaurus ingestion + pairwise alignment (0.5.3+, sub-tabs 0.9.5+)

Two-stage pipeline:

1. **Zip ingestion** — `_list_gbk_members_in_zip(path)` lists `.gbk`/`.gb`/`.genbank`; `_extract_gbk_member(path, name)` reads one out. Size-capped (`_PLASMIDSAURUS_ZIP_MAX_BYTES = 500 MB`, `_PLASMIDSAURUS_MEMBER_MAX_BYTES = 50 MB`, `_PLASMIDSAURUS_MAX_MEMBERS = 2000`).

2. **Structured parse** — `_parse_plasmidsaurus_zip(path)` groups per sample (`{gbk, fasta, summary, perbase, histogram, coverage_plot, interactive_map, ab1_files, summary_text, perbase_coverage}`). Run-level extras (`<run>_gel.png`, README) → `run_files`. Category folders matched on `_<suffix>` anchor (`_genbank-files`, `_fasta-files`, `_summary-files`, `_per-base-data`, `_histograms`, `_coverage-plots`, `_interactive-map`, `_ab1-files`); run-ID by majority vote. Standalone `.gbk` outside category folders still surfaced. Summary bodies (≤`_PLASMIDSAURUS_SUMMARY_MAX_BYTES = 4 KB`) streamed inline; per-base TSVs stream-summarised into `{mean, min, max, n_pos, above_20x}`.

3. **Alignment** — `_pairwise_align(query, target, mode='global'|'local')`. **Global mode uses edlib** (bit-parallel Myers, ~90–600× faster — `[INV-91]`) when `_EDLIB_AVAILABLE`, else falls back to `Bio.Align.PairwiseAligner`; local mode always uses Biopython. Returns `{mode, score, identity_pct, aligned_q, aligned_t, n_matches, n_mismatches, n_gaps, n_gap_opens_q, n_gap_opens_t, q_len, t_len}`. Capped at `_PAIRWISE_MAX_LEN = 200_000`. Counting is IUPAC-aware (`_iupac_compatible`); edlib gets the same IUPAC pairs via `additionalEqualities` so engine + counting agree. Biopython aligned strings come from `Alignment[0]`/`Alignment[1]`, NOT `format()`-parsing.

Entry: `Sequencing → Plasmidsaurus` on `SequencingScreen`. Nested `TabbedContent` with 4 sub-sub-tabs: **General** (always enabled, `_ZipAwareDirectoryTree` + overview), **Samples** (DataTable `#align-members`, row keys carry gbk member name), **Quality** (`#plasmidsaurus-quality-table` + `#plasmidsaurus-runfiles-table`), **Align** (query + target Select `#align-target` + `#btn-align-go`). Gating via `_apply_subtab_gating(enabled)`. Same-path re-pick short-circuits at top of `_on_zip_picked`. `PlasmidsaurusAlignModal` module-level alias of `SequencingScreen` for back-compat.

**Align flow (target-as-canvas, post-2026-05-24):** on click `#btn-align-go`, the worker (a) routes through `_pick_best_rotation(query_seq, target_seq, is_circular=…, canvas_axis="target")` which tries plain alignment first and — if circular and plain identity < `_ROTATION_TRY_THRESHOLD_PCT` (default 80) — also tries query-rotation via `_find_circular_alignment_offset(target, query)` and target-rotation via `_find_circular_alignment_offset(query, target)`. Picks whichever candidate has the highest overall `identity_pct` (more aligned bp = more informative overlay bar). When the picked rotation is on the canvas-axis side, aq/at are shifted back to the canvas's original frame via `_rotate_aligned_to_original_target_frame` / `_rotate_aligned_to_original_query_frame` so downstream segments land at the right bp on the linear map. (b) Loads the picked library target into the canvas via `_apply_record(target_record)` if it isn't already loaded — mirrors the Alt+A convention where the canvas plasmid is the reference and overlays are the queries. (c) Registers the alignment with `axis="target"` so the blue bar lands at the right bp on the now-active target's linear view. (d) `_flush_active_alignments` persists onto the target's library entry (canvas == target ⇒ correct entry). Display label is `f"{order_num} {gbk_basename}"` (1-based row index in `_parsed_run["samples"]` + `.gbk` leaf with extension stripped, e.g. `1 RUN42_1_MAV34`). Stored as `source="sequencing"` for the manager modal's batch-delete.

**Toast** surfaces aligned-bp / target-bp + coverage% + overall identity (gap-inclusive — predicts how much blue/red the band shows) + matched-region identity in parens + rotation note when one was applied. Pre-2026-05-24 it showed only the matched-region identity ("99.9% (aligned region)"), which misled the user for length-mismatched plasmids where the bulk of the target was gray (no read coverage) but the matched portion was near-perfect. See `[INV-68]`.

**Pre-2026-05-24:** worker rotated the *target* instead and registered against the rotated record — bars landed in rotated-frame coords on whatever the canvas happened to show (often unrelated to the target), and `_flush_active_alignments` persisted onto the canvas plasmid's library entry, not the target's. Re-loading the target wouldn't restore the alignment. See `[INV-67]`.

**Same rotation picker now powers Alt+\\ (`_diff_align_worker`)** — pre-2026-05-24 the diff-plasmid worker ran plain global alignment with no rotation, paying edge gaps for circular targets whose origin didn't match the canvas. Now routes through `_pick_best_rotation(canvas_axis="query")` (canvas = query for the diff flow). Matches the agent endpoint `_h_diff_plasmid` behaviour for circular targets.

**Honest identity + zoomed-out mismatch visibility + manager counts (2026-05-31, `[INV-87]`):** identity percentages in the AlignmentManagerModal + VerificationReportModal route through `_format_identity_pct` — a sub-100% value never rounds up to "100%" (escalates decimals until it reads < 100, e.g. "99.99%"), so the number can't contradict `_identity_pct_color`'s strict `>= 100.0` → light-blue tier (a one-bp mismatch read "100.0%" in green pre-fix). The bar-mode (< 1 col/bp) linear overlay collapses segments to per-column WORST state via `_alignment_bar_columns` (`_ALIGN_STATE_PRIORITY`: mismatch > gap > match) so a lone 1-bp mismatch paints RED even when one terminal column spans hundreds of bp — pre-fix its column range was degenerate AND the following match run repainted it blue. The AlignmentManagerModal gained **"Mism"** (mismatched bases, `n_mismatches`) + **"Gaps"** (indel *events* — gap runs — via `_alignment_indel_events`, **not** gapped bp `n_gap_cols`, so the count matches the Verification Report's "Indels" for the same alignment) columns for a fast "how bad is it" read.

**Bulk-confirm real quality columns (2026-05-31, `[INV-88]`):** the BulkAlignConfirmModal's "k-mer" column is a *matching* score (canonical-k-mer Jaccard), not alignment identity — a circular plasmid's origin-junction k-mers alone hold a perfect read to ~99.8% Jaccard, so it can't distinguish a 1-bp mismatch and `:.0%` rounds everything to "100%". The modal now also carries **Ident / Mism / Gaps** columns. These are aligned **once**, up front: `_compute_bulk_quality` runs on the `_match_samples_worker` button-press thread (one `_pick_best_rotation` per target-bearing sample, `canvas_axis="target"`, progress on the bar) and caches both `_aln` (display `{ident, mism, gaps}`) and `_aln_result` (full result) on each match dict. The modal just renders `_aln` (no in-modal alignment); `_bulk_align_worker` reuses `_aln_result` on commit so the read is never aligned twice. Rendering is the pure `_bulk_quality_cells(aln, *, action, has_target)` (em-dash / "…" / "?" / honest `_format_identity_pct`). `_aln` is stripped from the confirm payload; `_aln_result` rides through for the commit reuse. (First cut aligned the read 3× — modal-open preview worker + commit; reworked to align-once-up-front before shipping.)

**Bulk auto-align (2026-05-24, `[INV-69]`):** the **2. Pick sample** sub-tab carries a `Bulk auto-align all samples` button. `_bulk_align` calls `_match_samples_to_library(samples, library, sequence_fallback=True)` — for each sample, fuzzy-match by name (`_normalize_for_match` strips Plasmidsaurus run+order prefix + extension + punctuation, lowercases) with a k-mer Jaccard fallback (k=20) when name matching is ambiguous. `BulkAlignConfirmModal` shows one row per sample with the matcher's proposed `action` (align / add / skip); Space rotates the cursor row's action. On confirm, `_bulk_align_worker` batches: aligned rows route through `_pick_best_rotation` + `_register_alignment` + `_flush_active_alignments` (same path as the per-sample flow); add-as-new rows append to the library with `source: plasmidsaurus:<run_id>:<sample_name>` provenance stamp.

**Verification report (2026-05-24, `[INV-69]`):** `VerificationReportModal` collects every stored alignment across the active library, sorted by status priority (worst first: divergent → partial → near → verified). Per-row: plasmid name, read label, status badge (`_alignment_quality_status`), identity %, coverage %, # SNPs, # indels, source. Variants computed on-demand from `aligned_q`/`aligned_t` via `_extract_variants_from_alignment` (gap runs merged into single indel records). Click-row dismisses with `("open", entry_id, first_variant_pos)`; caller's `_jump_to_library_entry_at_pos` loads the plasmid + scrolls the seq-panel cursor to the first variant.

**LibraryPanel "Seq" column (2026-05-24, `[INV-69]`):** new 3rd column shows per-entry sequencing badge: ✓ verified (green) / ⚠ near-match (yellow) / ~ partial (yellow) / ✗ divergent (red) / — none (dim). Driven by `_library_entry_alignment_summary` walking the entry's stored `alignments`. Updates incrementally via `LibraryPanel.refresh_seq_cell(entry_id)` from `_flush_active_alignments` and the AlignmentManagerModal save path — no full repopulate per alignment flush.

Hardening (0.9.4):
* **Per-base TSV zip-bomb defence** — `_PLASMIDSAURUS_PERBASE_MAX_BYTES = 100 MB` two-layer cap (central-directory check + chunked `codecs.getincrementaldecoder` 64 KB).
* **Single-pass zip-open** — `_batch_extract_gbk_meta` reads every sample in one `ZipFile` open.
* **NUL-anchored sentinels** — `_NO_GBK_KEY_PREFIX = "\x00no-gbk\x00"`, `_EMPTY_LIBRARY_SENTINEL = "\x00no-library\x00"` (collision-proof).

## [SUB-experiments] Experiments lab-notebook (0.9.6+, projects refactor 0.9.7+)

Top-level (Menu → Experiments) → `ExperimentProjectsPickerModal` first (mirrors `PartsBinPickerModal` → `PartsBinModal`). Picking pushes full-screen `ExperimentsScreen`. Split-pane:

* **Top row** — active project label + `Projects… [^P]` button (Open/New/Rename/Duplicate/Delete/Close).
* **Left pane** — entries DataTable `Updated` + `Title`, natural-sort `updated_at` desc. Width 1fr (~20%) `min-width: 24`. New/Open/Rename/Delete.
* **Right pane** — `TabbedContent[Compose | Attachments]`. Disabled until entry loaded.
  * **Compose** — `TextArea` markdown source. Live preview dropped 2026-05-18; `_render_plasmid_refs` preserved for export.
  * **Attachments** — per-entry image grid via `ImageAttachModal`; Win/Mac clipboard paste via `Pillow.ImageGrab.grabclipboard()` (Linux/WSL disabled — no pure-Python API).

**Projects layer (projects:experiments :: collections:plasmids):** `experiment_projects.json` holds named projects, each with `experiments: list[dict]`.
* `_load_experiment_projects`/`_save_experiment_projects` — cache + `_cache_lock` + deepcopy.
* `_get_active_project_name`/`_set_active_project_name` — `settings["active_project"]`.
* `_ensure_default_project` — first-run wraps existing `experiments.json` into `_DEFAULT_PROJECT_NAME = "Main Project"`. Called from `compose()` per pitfall #9.
* `_sync_active_project_experiments` — called by `_save_experiments` after every save.

**Sacred — Experiments mirror:** every save MUST go through `_save_experiments` → `_sync_active_project_experiments(entries)`. Routing around bypasses mirror (same threat as pitfall #10).

`experiments.json` envelope-v1 schema:

```
{
  "id":                   "exp-<8 hex>",
  "title":                str,               # <= 200 chars
  "body_md":              str,               # <= 1 MB
  "created_at":           ISO-8601 w/ tz,
  "updated_at":           ISO-8601 w/ tz,
  "tags":                 list[str],         # max 20, <= 60 chars
  "attached_plasmid_ids": list[str],
  "image_paths":          list[str],         # relative to attach dir
}
```

**Plasmid cross-refs.** `@<id>` inline. Lookbehind `(?<![\w@])` rejects emails / double-`@`; id must start `[A-Za-z]`. `_render_plasmid_refs` → markdown links with `splicecraft://plasmid/<id>` href.

**Action cross-refs.** `!<id>` (distinct sigil). `(?<![\w!])` blocks word-adjacent / double-`!`; next char must be letter so `![alt](url)` doesn't false-match. `_EXPERIMENT_ACTIONS` curated catalog (19 entries) via `ActionsPickerModal`; free-form ids accepted.

**Gel cross-refs (2026-05-19).** `&<id>` references saved gels in `gels.json`. Orange chip `_GEL_CHIP_COLOR = "#FFB347"`. `_extract_gel_refs` denormalises into `attached_gel_ids`. Pick via `Gel ref` button → `GelLibraryModal`.

**Click-to-open / Ctrl+G.** Ctrl+G / double-click scans line for `@`/`!`/`&` tag spanning cursor column. Plasmid hit → auto-save dirty compose → search every collection → switch + load → dismiss. Gel hit → `GelLibraryModal(initial_gel_id=<id>)`. Action hit → `ActionsPickerModal(initial_action=<id>)`.

**Legacy tag migration.** Pre-2026-05-18 `@plasmid:<id>`/`@actions:<id>` rewritten on every `_load_experiments` via `_migrate_legacy_tag_format`. One-way; once saved through `_save_experiments`, old format gone.

**In-editor token coloring.** `_ExperimentMarkdownTextArea` overrides `_build_highlight_map` for `@<id>` (lime `_PLASMID_CHIP_COLOR = "#9AFF80"`), `!<id>` (purple `_ACTIONS_CHIP_COLOR = "#C77FFF"`), `&<id>` (orange `_GEL_CHIP_COLOR = "#FFB347"`). Highlight names `splicecraft.plasmid_ref`/`.action_ref`/`.gel_ref` injected into theme's `syntax_styles` via `setdefault`. ASCII fast-path vs non-ASCII path (codepoint→byte table). Backspace at end of tag deletes whole tag. `on_click` intercepts only double-clicks (event.chain ≥ 2) for `TagOpenRequested`.

**Sacred sizing caps:**
* `_EXPERIMENT_BODY_MAX_BYTES = 1_000_000` per entry (deterministic truncate).
* `_EXPERIMENT_IMAGE_MAX_BYTES = 10_000_000` per image.
* `_EXPERIMENT_DIR_MAX_BYTES = 100_000_000` per entry.

**Filesystem invariants** (mirror .dna sidecar):
* `_sanitize_experiment_id` rejects empty/NUL/`..`/`/`/`\`/shell metas/>64 chars.
* `_experiment_attach_dir` walks FULL ancestor chain via `is_symlink()`. No `resolve()` divergence (would trip on macOS `/tmp` → `/private/tmp`).
* `_save_experiment_image` via `_atomic_write_bytes`. Filename `img-<ts>-<rand>.<ext>`. Clipboard tmpfiles (prefix `_EXPERIMENT_CLIP_TMP_PREFIX = "exp-clip-"`) unlinked after bytes copied.
* `_save_experiments` takes `_cache_lock` for save+cache-reassign, then `_sync_active_project_experiments`.
* `_persist_current` detects body-over-cap BEFORE save (notifies). Save path dedup-by-id replaces ALL matches.

**Spellcheck.** pyspellchecker-backed (pure-Python English wordlist). F7 / "Spellcheck" → `_spellcheck_body(body_md)` masks non-prose markdown regions and tokenises via `_SPELLCHECK_WORD_RE` (alphabetic + apostrophe + hyphen, ≥2 chars). `SpellcheckModal`: Replace/Add-to-dict/Skip per row. Custom dict via `experiments_custom_dict` settings key; `_clear_spellcheck_engine` invalidates cached engine after add.

**Hard deps:** `Pillow>=10.0`, `pyspellchecker>=0.8.0`, `rich-pixels>=3.0.0` (pure-Python wheels).

**Modal `_blocks_undo=True`** on `ExperimentsScreen` + `SpellcheckModal`.

**Unsaved-changes guard.** `ExperimentsScreen.action_cancel` pushes `ExperimentUnsavedChangesModal` (Save/Abandon/Close; default Close, Esc → cancel) when buffer dirty. Screen callback stays on top if Save fails so user can retry. Delete paths use `ExperimentDeleteConfirmModal` default No (sacred — stray Enter cannot delete data). `ExperimentProjectsPickerModal._do_delete` re-checks last-project guard inside confirm callback.

**Logging events:** `experiments.*` (new/save/delete/attach.image/remove.image/insert.{plasmid,action,gel}_ref/spellcheck.applied/tag.migrated), `project.*` (switched/created/renamed/duplicated/deleted), `gel.*` (created/renamed/deleted/loaded/ref.opened), `plasmid.ref.opened`, `action.ref.opened`. `_log_event` sanitises body (200-char truncation).

## [SUB-gels] Gels (saved agarose-gel snapshots, 2026-05-19+)

`gels.json` (`_GELS_FILE`) holds saved Simulator gel configurations. Schema envelope-v1:

```
{
  "id":          "gel-<8 hex>",
  "name":        str,               # <= 200 chars
  "lanes":       list[dict],        # [{name, source, detail}, ...] cap 20
  "agarose_pct": float,             # clamped 0.3–5.0; NaN/inf rejected
  "notes":       str,               # <= 2000 chars
  "created_at":  ISO-8601 w/ tz,
  "updated_at":  ISO-8601 w/ tz,
}
```

Helpers mirror experiments/projects: `_load_gels`/`_save_gels` with `_cache_lock` + deepcopy-on-read+save, `_safe_save_json` for atomic write. `_sanitize_gel_id` rejects empty/NUL/`..`/`/`/`\`/>64. `_normalise_gel_entry` caps all string fields, drops non-dict lanes, clamps agarose, replaces invalid ids with fresh `gel-<hex>`. `_find_gel(id)` returns None for unsanitisable ids (defensive against `&../etc`). `_gel_name_taken(name)` strip-compare case-sensitive.

`GelLibraryModal` dual-context picker:
* From `SimulatorScreen.Gel` (Library button) — opens with `current_lanes` + `current_agarose_pct`; Save current enabled; dismiss-with-id restores lanes/agarose and re-renders.
* From `ExperimentsScreen` (Gel ref button) — no current snapshot; Save current disabled; dismiss-with-id inserts `&<id>`.
* From click-to-open — opens with `initial_gel_id` scrolled to row.

`SimulatorScreen` has no persistent state; live `self._lanes`/`self._agarose_pct` are in-memory, only written via `_save_gels` on explicit save. Delete-last-gel allowed (no "active gel" concept).
