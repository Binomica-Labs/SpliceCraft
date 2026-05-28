# Features

What you can do without leaving the terminal.

## View

- **Braille dot-matrix circular maps** — plasmids rendered as crisp
  Unicode braille rings with per-strand feature arcs, directional
  arrowheads, and proximity-placed labels. `v` toggles linear view.
- **Per-base sequence panel** with two-strand display, wrap-aware
  feature lanes, restriction-site overlays, and inline AA translation
  (one letter per codon midpoint, in the CDS's colour, with wrap-CDS
  support across the origin). Click an AA letter to highlight the
  codon's three bases on the strand.
- **Per-strand restriction-cut visualisation** — clicking a sticky
  cutter (EcoRI, HindIII, BsaI, BsmBI, BbsI, …) tints upstream bases
  blue and downstream red, with the staggered overhang showing as
  different colours on the two strands.
- **200+ NEB enzymes** including Type IIS scanners; toggle restriction
  overlays with `r`, filter to unique cutters / 6+ bp / connectors.

## Edit

- **In-place sequence edits** with full undo / redo (50-deep snapshot
  stack, deepcopied SeqRecord). Per-plasmid undo stashes — switch
  records, edit, switch back, undo history is restored.
- **Feature CRUD**: add / merge / split / delete / rename / recolor
  features; clipboard copies (top strand or reverse-complement bottom
  strand). Mouse-drag selects ranges; `Enter` highlights the smallest
  feature enclosing the cursor.
- **Crash-recovery autosave** writes a 3-second-debounced `.gb`
  snapshot to the data dir; survivors surface on next launch.

## Synthesis

- **DNA synthesis composer** (`Synthesis → DNA`) — paste / type a
  sequence target, set topology, then iteratively annotate using the
  built-in feature library or motifs found via inline pattern match.
  Synthesised records save to the library with feature annotations,
  topology, and a `synthesised in SpliceCraft` provenance note.
- **Protein synthesis composer** (`Synthesis → Protein`) — design at
  the AA level. Paste a protein sequence, pick a codon-usage table
  (E. coli K12 / S. cerevisiae / H. sapiens / +100 Kazusa tables), and
  the composer back-translates to optimised DNA while scrubbing
  forbidden Type IIS sites against the active cloning grammar.
  In-editor protein-motif highlighting (signal peptide, NLS, mito
  targeting, tags) flags regions of biological interest as you type.
- **Codon-table picker** — the **Codon table** dropdown in the Synthesis
  (and Mutato) protein tools opens a picker to browse the Kazusa
  usage-table catalog by organism / taxid and persist a preferred table
  (also via `set-active-codon-table`). **Import TSV** pastes a custom
  usage table (codon, optional amino acid, count) straight into your
  library.

## Cloning

- **Cloning grammars** — GB L0 (Esp3I) and MoClo Plant (BsaI) ship as
  built-ins; user-defined grammars persist to `cloning_grammars.json`
  and are editable in `GrammarEditorModal`. The active grammar
  parameterises the Domesticator, Parts Bin, and Constructor — change
  enzyme / overhang / forbidden-site set without code edits.
- **Domesticator** — 4-source part picker (current map, library,
  Parts Bin, FASTA file). Auto-scrubs forbidden Type IIS sites in the
  CDS body via codon swap with cascade-prevention; primer tails follow
  the active grammar's pad / site / spacer / overhang.
- **Parts Bin** — domesticated parts catalog with per-grammar filtering;
  legacy parts default to GB L0; "Copy primed sequence" preserves the
  part's stored grammar. **Load Part** auto-classifies the currently-
  open plasmid by digesting it with each grammar's Type IIS enzyme and
  matching the released fragment's overhangs against the grammar's
  position table — register an externally-domesticated part without
  manually picking grammar / position.
- **Constructor** — multi-tab assembly UI: Traditional restriction
  cloning, Golden Braid / MoClo Type IIS assembly, and **Gibson
  assembly**. The Gibson tab stages N linear fragments, detects the
  longest exact-match overlap at every junction (incl. the wrap
  junction for circular topology), validates against a configurable
  minimum, and produces a single assembled product with each overlap
  appearing once. Reverse-orientation fragments surface a "did you
  mean to flip" hint instead of silently failing.
- **Traditional cloning** — restriction-digest + ligation simulator
  with three insert sources (current plasmid, library entry, free-form
  PCR product). 2-enzyme directional cuts produce both forward and
  reverse-orientation products; non-ligatable orientations are flagged
  rather than silently dropped. Save the simulated product back to the
  library with full **construction-history XML** (`<HistoryTree>`
  matching the popular commercial editor's format) so the lineage of
  multi-step builds is preserved across import/export.

## Primer design

- **Primer design** — detection / cloning / Golden Braid / generic via
  Primer3; primers can be added to the map as `primer_bind` features
  or saved to the persistent primer library (Designed → Ordered →
  Validated lifecycle).

## Mutagenesis

- **SOE-PCR site-directed mutagenesis** — design 4-primer SOE sets for
  any W140F-style point mutation. CDS source can be the loaded plasmid,
  a library entry, a Parts Bin part, or a free-form protein sequence
  (auto-optimised via the active codon table). Edge cases (mutation
  within 60 nt of a CDS end) auto-fall back to a 2-primer modified-outer
  PCR.

## Simulate

- **In-silico PCR + agarose gel** (Simulator menu) — design a primer
  pair against the loaded plasmid and the simulator enumerates every
  legal amplicon (exact-match binding model, wrap-aware on circular
  templates, capped at 50 results to flag mispriming runaway).
  Amplicons round-trip to the library as linear DNA entries with
  `primer_bind` features at both ends.
- **Agarose gel renderer** — paint up to 8 lanes (ladder / uncut
  plasmid / restriction digest / PCR amplicon) on a virtual gel at
  user-selectable agarose % (0.5 → 4.0). Mobility uses the Helling-
  Goodman-Boyer empirical curve (distance ∝ −log₁₀ bp within each
  agarose's resolution window) plus the standard form corrections —
  supercoiled migrates faster than linear, nicked / open-circle
  slower. Lane sources share the screen's template, so the amplicon
  designed in the PCR tab is immediately runnable in the Gel tab.

## Search

- **In-process BLAST** (`Ctrl+B`):
  - **BLASTN** (DNA → DNA) and **BLASTP** (protein → protein) via
    `pyhmmer.hmmer.nhmmer` / `phmmer` (HMMER 3 in-process at C speed);
    pure-Python ungapped fallback for queries below the HMMER profile-
    builder minimum (20 bp / 6 aa).
  - **HMMscan** reads any HMMER 3 `.hmm` / `.h3m` / `.h3p` file
    directly — point it at Pfam-A or any custom profile DB. Lazy file
    read so Pfam-scale (~1 GB) DBs don't pre-fetch into RAM.
  - DB build + search run in a `@work(thread=True)` worker; UI stays
    responsive on a 50-plasmid index. 4-entry LRU DB cache, auto-
    invalidated on `_save_collections`.
- **Six-frame ORF indexing** (opt-in checkbox) for BLASTP against
  unannotated regions of plasmid backbones.
- **Cross-collection plasmid search** — Edit → Find plasmid… opens a
  fuzzy / substring search over every plasmid in every collection,
  natural-sorted by `(collection, plasmid)` so `pBin2` lands before
  `pBin10`. One click opens the entry without manually switching
  collections.
- **Pairwise alignment of sequencing runs** — File → Align sequencing
  run loads a Plasmidsaurus `.zip` (or any `.gbk` / `.gb`), pairwise-
  aligns it against the loaded plasmid, and renders a full-screen
  alignment viewer with target-feature lane, parallel target/query
  rows, match track, and mismatch-red highlighting. Length-capped at
  200 kb per side; cancellable via the standard worker pattern.
- **Sanger trace viewer (`.ab1`)** — Sequencing → Sanger tab.
  Browse a directory for AB1 traces (highlighted sky-blue), pick
  one to see the base-called length, mean Phred quality, and a
  preview of the first 200 bases. Load straight onto the canvas or
  add to the library for downstream alignment. BioPython base-calls
  the trace channel; quality scores survive on
  `letter_annotations["phred_quality"]` for any caller that wants
  them.
- **New Plasmid modal** (`Ctrl+N`) — paste a sequence, optionally name
  + set topology, then either Create / Annotate-from-library
  (substring match) / Annotate-via-BLAST (≥90% identity → `misc_feature`).

## Library

- **Plasmid collections** — named buckets (e.g. "yeast project",
  "E. coli toolkit"); the panel toggles between a collection list and
  the active collection's plasmids. Atomic writes, `.bak` per change.
  Save the loaded record with `Ctrl+Shift+A`.
- **Bulk import a folder** — from the collections-list view, click `+`,
  type a name, and pick a folder via the embedded directory tree.
  Every `.dna` / `.gb` / `.gbk` / `.genbank` file inside is loaded
  independently into a new collection; failures are isolated per file
  and surfaced in a notify summary. Designed for migrating a
  popular-commercial-plasmid-editor archive in one shot.
- **`.dna` round-trip.** SpliceCraft reads the popular commercial
  plasmid editor's binary format (sequence + features + notes +
  primers + construction history) and writes it back — including the
  default `Primers` and `AdditionalSequenceProperties` packets the
  editor itself emits — so files round-trip through SpliceCraft
  cleanly into the editor's Viewer / Inspector panels. Imported
  primers feed into the persistent primer library (de-duplicated by
  sequence), and per-feature colours are recovered alongside.
  Construction history XML is preserved on import and synthesised on
  save for any product built via the Traditional cloning simulator.
- **Construction history viewer.** `File → View construction history`
  renders any record's `<HistoryTree>` lineage — fragments, enzymes,
  parent products — as a navigable tree so the provenance of a
  multi-step build is auditable at a glance.
- **Library fuzzy search** — subsequence match (case-insensitive,
  non-contiguous) against the visible table; natural-sorted so
  `pBin2` lands before `pBin10`.
- **Feature library** — reusable feature snippets (per-entry colour
  and strand) with a centralised browse / edit / rename / recolor /
  delete workbench. Display rows natural-sort independently of the
  on-disk order so `pPart-2` sits next to `pPart-10` rather than
  scattered alphabetically; entry indices remain stable across the
  re-sort so dirty-edit markers don't desync.
- **Bulk export collection** (`File → Export collection (bulk)…`) —
  pick a collection + format (GenBank / EMBL / FASTA / `.dna`) + a
  target folder. Each plasmid is written as `<name>.<ext>` with
  filesystem-safe sanitisation (path-traversal characters scrubbed,
  Windows reserved device names prefixed, case-insensitive collision
  defence on APFS/NTFS). Per-entry failures don't abort the run; the
  summary toast reports written / failed counts.

## File formats

| Format | Extensions | Import | Export | Preserves |
|---|---|---|---|---|
| GenBank | `.gb` / `.gbk` / `.genbank` | yes | yes | features, qualifiers, wrap topology |
| EMBL | `.embl` | yes | yes | features, qualifiers, wrap topology |
| CommercialSaaS `.dna` | `.dna` | yes | yes | features + colours + primers + construction history (round-trip) |
| FASTA | `.fa` / `.fasta` / `.fna` / `.ffn` / `.frn` / `.fas` / `.mpfa` / `.faa` | yes | yes (sequence only) | sequence only |
| Sanger trace | `.ab1` / `.abi` | yes | — | base-called sequence + Phred quality |
| FASTQ multi-read | `.fastq` / `.fq` | yes (≤ 1000 reads per file) | — | one library entry per read |
| GFF3 | `.gff` / `.gff3` | yes (standalone via `##FASTA`; or apply features to loaded canvas via `apply-gff3`) | yes | wrap features as same-`ID=` split rows; `Is_circular=true` on region row |
| Plasmidsaurus zip | `.zip` | yes (Sequencing → Plasmidsaurus tab) | — | consensus + run-level QC + AB1 traces |

All file reads route through size-cap + symlink-refusal checks; all
file writes route through `_atomic_write_text` / `_atomic_write_bytes`
(tempfile + fsync + replace + symlink refusal). Bulk-import and
single-file Open share one dispatch table so the agent CLI, GUI Open,
and folder-import all accept the same set.

## Experiments lab notebook

- **Projects layer** (`Menu → Experiments`) — named projects (e.g.
  "Yeast Y32 strain", "E. coli toolkit") each hold a list of
  experiment entries. Active project persists across sessions; first
  launch wraps the user's existing experiments into a default
  project.
- **Compose + Attachments** — per-entry markdown body (1 MB cap),
  plus an image attachment grid. Win/Mac clipboard paste via
  `Pillow.ImageGrab.grabclipboard()` (Linux/WSL disabled by Pillow).
- **Cross-refs** — `@<plasmid-id>` inlines a coloured chip linking
  to a library plasmid; `!<action-id>` references the curated
  `_EXPERIMENT_ACTIONS` catalog; `&<gel-id>` references a saved
  gel snapshot. Ctrl+G / double-click on any tag opens the
  referenced entity.
- **Spellcheck** (F7) — pyspellchecker-backed (pure-Python English
  wordlist) with markdown-aware masking. Custom dict per-user.

## Gels

- **Save gel snapshots** — gel images created in `Simulator → Gel`
  can be saved to `gels.json` for later reference, side-by-side
  comparison across timepoints, or for citation from experiment
  entries via `&<gel-id>` references.
- **Gel library** — `Simulator → Gels…` browses every saved gel
  with thumbnail rendering, lane composition, and gel-percentage
  metadata.

## Protein motifs

- **Curated motif catalog** — 30+ patterns (NLS, signal peptide,
  mitochondrial targeting sequence, common tags, phosphorylation
  sites). Surfaced in the Synthesis protein composer as you type,
  and in the AA lane on the main sequence panel under loaded CDSs.
- **User-overrides** — add (**New…**), edit, and remove (**Delete**)
  motifs persistently from the Synthesis protein tab's motif pane;
  overrides round-trip through `protein_motifs.json` with the same
  atomic-save + backup discipline as every other persisted file.
  Built-in motifs are protected — deleting your edit of one restores
  the original.

## Recovery + data safety

- **Four-layer JSON safety net** — every persisted file gets atomic
  write + single-gen `.bak` + rotating timestamped `.bak.<ts>` (10
  retained) + daily snapshot in `<DATA_DIR>/snapshots/` (30 retained).
  Suspicious-shrink guard (≥50% data loss) spills affected entries
  to `<DATA_DIR>/lost_entries/` BEFORE the overwrite lands.
- **Restore from backup** (`Settings → Restore … from backup…`) —
  one modal lists every recoverable copy of any user-data file
  across the four storage tiers. Damaged rows are surfaced tagged
  `[damaged]` instead of silently dropped — you can see what was
  there even if it can't be restored.
- **Master Delete** (`Settings → Master Delete…`) — a typed-`YES`
  guarded recovery affordance that resets the entire SpliceCraft
  data dir to a fresh-install state. Two-stage confirmation +
  enumerated cache reset so a partial wipe doesn't leave the in-
  memory caches stale relative to disk.
- **Pre-update snapshots** before any pip / pipx / uv subprocess;
  stored in a sibling directory so a hypothetical recursive-wipe
  bug in a new version cannot reach the snapshots.

## Drive it from outside the GUI

See [Agent API](agent-api.md) and [CLI sidecar](cli.md) for the
details. In short:

- **Agent API** (`splicecraft --agent`) exposes a localhost JSON API
  with bearer-token auth, covering every GUI action external AI
  agents need. 90+ endpoints; symlink-guarded write paths;
  length/range/shape validation at the boundary.
- **`splicecraft-cli`** — stdlib-only sidecar (~50 ms cold start)
  that reads connection details from the running session's token
  file. Intended for Claude Code, Cursor, aider, hand-rolled
  scripts, or any external automation.
