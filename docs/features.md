# Features

What you can do without leaving the terminal.

## View

- **Braille dot-matrix circular maps** — plasmids rendered as crisp
  Unicode braille rings with per-strand feature arcs, directional
  arrowheads, and proximity-placed labels. `v` toggles linear view.
- **Export the map as a publication-quality image** — File → *Export map
  image* renders the current plasmid as a circular
  **PNG or SVG**: feature arrowheads with leader-line labels, a
  de-collided restriction-site tick ring, and a centre name + bp block.
  Options for size (300–6000 px), a transparent background (for figure
  overlays), and toggling labels / sites. SVG is true vector; PNG is
  super-sampled + LANCZOS-downscaled for clean anti-aliasing.
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
- **Orientation**: `Alt+Shift+R` flips the whole record (reverse
  complement) — every feature is re-framed so it still covers the same
  bases read from the other strand, and ▶/◀ arrows swap; arrowless (▒)
  and double-stranded (◀▶) features keep their marker, since neither has
  a direction to reverse. `Alt+Shift+O` re-cuts a **circular** record so
  the cursor becomes base 1, from either map view — a real edit that
  persists into saves and exports, unlike the map's `Alt+O` / `[` `]`,
  which only spin the display. It refuses on a linear record: rotating
  one would join two ends that aren't joined. Both are undoable, and
  both clear on-screen alignments, whose coordinates no longer apply.

## Synthesis

- **DNA synthesis composer** (`Synthesis → DNA`) — paste / type a
  sequence target, set topology, then iteratively annotate using your
  own feature library, the built-in preset catalogue (rows marked
  `pre`), or motifs found via inline pattern match. **Insert** splices a
  preset's bases at the cursor and annotates them in one step, carrying
  the GenBank accession it came from into the feature's note, so a
  construct built here records where each element originated.
  **Annotate** overlays a preset onto the current selection without
  changing any DNA. **Edit** on a preset saves an editable copy into
  your library and leaves the shipped catalogue alone.
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
  library. **Build** generates a fresh table from an NCBI genome
  assembly or a local CDS FASTA (`cds_from_genomic.fna` / `.gz`, fully
  offline) in highly-expressed-genes or whole-genome mode (also via
  `add-codon-table` with `source: "genome" | "file"`). **Chart** views
  any table as the classic genetic-code grid —
  each codon annotated with its within-family usage (relative synonymous
  usage), each amino-acid family's single most-used codon highlighted in
  bold green (family-wide; stops counted as one family).

## Cloning

- **Cloning grammars** — GB L0 (Esp3I) and MoClo Plant (BsaI) ship as
  built-ins; user-defined grammars persist to `cloning_grammars.json`
  and are editable in `GrammarEditorModal`. The active grammar
  parameterises the Domesticator, Parts Bin, and Constructor — change
  enzyme / overhang / forbidden-site set without code edits.
- **Domesticator** — 4-source part picker (current map, library,
  Parts Bin, FASTA file). Auto-scrubs forbidden Type IIS sites in the
  CDS body via codon swap with cascade-prevention; primer tails follow
  the active grammar's pad / site / spacer / overhang. On save, one
  dialog names the cloned plasmid + picks its collection + the parts
  bin (independently of the linear synthesis fragment); the cloned
  plasmid is mirrored to the library carrying its domestication primers
  as `primer_bind` features (bound anneal+overhang vs unbound enzyme
  tail), routed through the shared `_commit_library_entry_to_collection`.
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
  mean to flip" hint instead of silently failing. A circular plasmid
  source can be **linearised at any base** ("Linearize at"), and
  **Design overlaps** auto-appends the homology arms (each downstream
  fragment gets a 5' tail from its upstream neighbour's 3' end) so
  every junction reaches the minimum — idempotent, re-run after a
  reorder.
- **Traditional cloning** — restriction-digest + ligation simulator
  with three insert sources (current plasmid, library entry, free-form
  PCR product). 2-enzyme directional cuts produce both forward and
  reverse-orientation products; non-ligatable orientations are flagged
  rather than silently dropped. **Self-ligation risk** is reported
  alongside the product: a backbone whose two ends match each other —
  one enzyme, two blunt cutters, or a compatible-cohesive pair such as
  SalI and XhoI, which both leave TCGA — re-closes without any insert
  and is the usual source of empty colonies. The report names the fix
  (dephosphorylate, or pick different overhangs) and says whether a
  diagnostic digest could even spot an empty-vector colony. It also
  lists the other circles the tube can make: two inserts in tandem or
  head-to-head, a shorter chain of fragments closing without the rest,
  an insert flipping in place, neighbours swapping order. A reaction
  with no off-target route says so explicitly, naming which ends
  disagree, so silence never has to be read as "not checked". Each junction's sticky overhang is
  annotated at the bases that actually anneal there, derived from the
  enzyme's own 5'/3' geometry, and a cloning site the ligation
  regenerates comes back as one intact feature instead of two halves
  marked "(disrupted)" — a feature the cloning genuinely broke still is.
  Save the simulated product back to the
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
- **Scrub (clone-free restriction-site removal)** — the Mutato modal's
  second tab cures the plasmid on the canvas of chosen recognition sites
  (default BsaI / Esp3I / BbsI) using the minimal point substitutions
  that destroy each site: synonymous swaps inside any overlapping CDS
  (verified against every reading frame, incl. opposite-strand genes, and
  preferring the chosen codon table's frequent synonyms), minimal base
  changes elsewhere, never introducing a new forbidden site.
  Substitution-only, so feature coordinates never move. Produces the
  cured circular plasmid (Apply to canvas, undoable) plus primers sliced
  from the cured template (binding == display; save to the primer library
  or add straight to the map). A **Re-circularize** selector picks the
  lab route:
  - **QuikChange** — one improved partial-overlap pair per locus, PCR →
    DpnI → transform (no ligase, no assembly).
  - **Golden Braid** — split the plasmid at each cure into BsaI-tailed PCR
    fragments that a one-pot Golden Gate (BsaI) reaction reassembles. The
    4 nt junction overhangs are native sequence, so the product is seamless
    except the cures. BsaI is the assembly enzyme, so **every** BsaI site is
    force-cured (an internal one would be cut mid-reaction); a BsaI site
    that can't be silently removed makes Golden Braid refuse (use
    QuikChange). The design self-checks by a digest+ligate simulation that
    must reassemble to the cured plasmid before commit.
  Un-curable sites are reported, not forced. Headless via the
  `scrub-plasmid` agent endpoint (`method: "quikchange" | "golden_braid"`).

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
- **Online BLAST / HMM** (second tab) — submit DNA / protein, a whole
  library plasmid, or a single feature to NCBI BLAST (blastn / blastp /
  blastx / tblastn / tblastx) or EMBL-EBI Pfam (hmmscan), with a live
  poll counter and a real Cancel. **Add to collection** fetches a
  highlighted nucleotide hit's full GenBank record (efetch) and saves it
  as a plasmid into a collection you pick — name + destination prompt,
  collision-safe via the universal commit helper. Enabled only for
  nucleotide-subject hits (blastn / tblastn / tblastx); a protein or
  Pfam hit has no DNA record to land in a plasmid collection.
- **Cross-collection plasmid search** — Edit → Find plasmid… opens a
  fuzzy / substring search over every plasmid in every collection,
  natural-sorted by `(collection, plasmid)` so `pBin2` lands before
  `pBin10`. One click opens the entry without manually switching
  collections.
- **Inverted-region detection** — a construct whose insert went in
  backwards used to align as plain divergence, because comparing a
  flipped region forwards scores like random sequence. Poorly-matching
  stretches are re-checked against the reverse complement; a stretch
  that matches is marked as inverted on the linear map (magenta), given
  its own grade in the Verification Report, and reported with its exact
  span and identity. A read that aligns wholly as the reverse
  complement is called out too.
- **Pairwise alignment of sequencing runs** — File → Align sequencing
  run loads a Plasmidsaurus `.zip` (or any `.gbk` / `.gb`), pairwise-
  aligns it against the loaded plasmid, and renders a full-screen
  alignment viewer with target-feature lane, parallel target/query
  rows, match track, and mismatch-red highlighting. Length-capped at
  200 kb per side; cancellable via the standard worker pattern.
- **Fetch sequencing runs from the Plasmidsaurus API** — Sequencing →
  "Fetch by item code from Plasmidsaurus API". Pulls a finished run from
  your Plasmidsaurus account over the official OAuth2 REST API and
  imports its `.gbk` assemblies into the library as new entries
  (tagged `plasmidsaurus:<code>:<sample>`, never overwriting). Opening it
  ENUMERATES your orders (newest first, with the order name you typed)
  into a table — select a row and press Enter to download it, or type a
  code for an order someone shared with you. Orders with no results —
  shipping labels and canceled orders, roughly 40% of a real listing —
  are hidden behind the "Show orders with no results" checkbox rather
  than left to 404 at download. The listing and its access token are
  cached for two minutes, because the account is capped at 10 API
  requests per minute; Refresh forces a re-fetch. Credentials
  resolve env-first (`PLASMIDSAURUS_CLIENT_ID` / `_SECRET`) then
  Settings ▸ Plasmidsaurus API. Scriptable via the `plasmidsaurus-items`
  / `download-plasmidsaurus` agent endpoints.
- **New L0 part from a synthetic fragment** — Constructor → any modular
  tab → **New Part from Syn Frag**. Pick a saved FRAG (one built by
  Synthesis → L0 Fragment) and the part type it will become; SpliceCraft
  loads the entry vector assigned to that tab's grammar and *simulates
  the cloning* — digests both the fragment and the vector with the
  grammar's Type IIS enzyme, ligates the released insert into the
  backbone and closes the circle — then files the result as a Level-0
  part, ready in the palette, **and saves the circular L0 plasmid to
  your library** as *"<part> (L0)"* with its construction lineage
  attached, so History shows the insert and the entry vector it went
  into. The part's fusion overhangs come from the
  grammar's own position table, so a two-tier fragment (entry overhangs
  wrapping the category pair) files the *category* pair, which is what
  L0→L1 actually uses. Only that grammar's positions are offered, and a
  fragment that doesn't carry the chosen type's overhangs is refused
  rather than filed as a part that can't assemble. Unwrapped fragments,
  a missing entry vector and duplicate names each fail with the specific
  reason. Scriptable — and available to Babs — via the
  `make-l0-part-from-fragment` agent endpoint.
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
- **ORF finder** (File ▸ Find ORFs, or the command palette) — six-frame
  scan with a minimum length in **amino acids**; ATG-only by default,
  with GTG / TTG opt-in for bacterial starts. Wrap-aware on circular
  plasmids: an ORF crossing the origin is listed with its end before its
  start and marked *(wrap)*. Picking a row highlights that ORF in the
  sequence panel.

  Read the **aa** column (or the API's `length_aa` / `nt_len`) for
  length — never the start/end pair. On a circular plasmid a frame can
  run the whole way round without hitting a stop, and no start/end pair
  on an n-bp circle can express an ORF of a full lap or more. Those are
  listed as **full lap** with their span pinned to the near-full circle;
  the reported length stays exact.

## Library

- **Plasmid collections** — named buckets (e.g. "yeast project",
  "E. coli toolkit"); the panel toggles between a collection list and
  the active collection's plasmids. Atomic writes, `.bak` per change.
  Add the loaded record to your library with `Alt+K` ("keep").
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
- **Built-in feature presets** (`Presets` in the Feature Library, or
  `p`) — a shipped catalogue of 117 common plasmid elements across 20
  categories: resistance markers (AmpR, KanR, NeoR, CmR, TetA, SpecR,
  HygR, PuroR, BlastR, ZeoR), origins (pUC/pMB1, p15A, f1, SV40, 2μ,
  CEN/ARS, R6Kγ, oriT), promoters for bacteria, mammalian cells, yeast
  and plants (T7, T3, SP6, lac, tac, trc, pBAD, J23119, CMV, EF-1α,
  PGK, U6, GAL1, 35S, NOS), terminators and polyA signals, reporters
  (EGFP, sfGFP, mCherry, mRFP1, firefly luciferase, lacZ-α), affinity
  tags and protease sites, 2A peptides, and recombination sites (loxP,
  FRT, attR, attB, I-SceI). Browse with search and a category filter,
  tick as many as you want, and copy them into your own library —
  where they become ordinary entries you can rename, recolour or edit.
  Every sequence was extracted from a public GenBank record rather
  than typed from memory, and each entry cites the accession plus how
  many independent submissions carry the identical sequence. Coding
  entries were checked for a real reading frame. Presets live in the
  program, never in your data directory, so they cost you nothing and
  cannot be corrupted; an entry of your own with the same name and
  type always replaces the preset in every list.
- **Annotate the loaded plasmid** (`File → Annotate from library +
  presets`) — matches the plasmid on screen against your own snippets and
  the preset catalogue, on both strands and across the origin, then shows
  every proposed feature before any of them land. One `Ctrl+Z` removes the
  whole batch. An element the record already carries under the same name
  and coordinates is not offered again, so re-running changes nothing. The
  scan runs on a worker thread, so a 200 kb record doesn't stall the UI.
- **Annotate a pasted sequence** — `Library → New plasmid` takes a
  pasted sequence and its **Annotate from library** button
  substring-matches it against your own snippets *and* the preset
  catalogue on both strands, wrap-aware on a circular template, so a
  freshly pasted vector comes back with its markers, origin, promoters
  and polylinker already drawn. Matching is exact: a vector carrying an
  allelic variant of an element (pGEX-4T-1's *bla* differs from the
  preset at two bases) will not match, which is what **Annotate via
  BLAST** is for.
- **Bulk export collection** (`File → Export collection (bulk)…`) —
  pick a collection + format (GenBank / EMBL / FASTA / `.dna` /
  circular-map **PNG** / **SVG**) + a target folder. Each plasmid is
  written as `<name>.<ext>` with filesystem-safe sanitisation
  (path-traversal characters scrubbed, Windows reserved device names
  prefixed, case-insensitive collision defence on APFS/NTFS). Per-entry
  failures don't abort the run; the summary toast reports written /
  failed counts.
- **Export marked plasmids as images** — mark rows in the library, press
  `p`, and every marked plasmid is rendered to a circular-map PNG / SVG
  in a folder you pick. Runs on a worker thread with per-entry failure
  isolation; ids resolve across every collection. Takes rows carrying
  *either* mark — export doesn't distinguish move from copy.
- **Staging plasmids (Ⓜ / Ⓒ / Ⓧ)** — **Space** cycles the highlighted
  library row through **none → Ⓜ (move) → Ⓒ (copy) → Ⓧ (delete) →
  none**, so the mark itself says what will happen to that plasmid —
  the same idea as the primer library's ★ / $ / M cycle. Press **`m`**
  to send the Ⓜ rows to a collection you pick, **`y`** to send the
  Ⓒ rows, and **Delete** to remove the Ⓧ rows. A mixed batch is staged
  in one pass and committed in as many, and committing one kind leaves
  the others marked. With nothing marked, `m` / `y` / Delete still act
  on the row under the cursor. Copies are deep-copied into the target
  and renamed on collision (`COPY` for the same sequence, `ALT` for a
  different one); `c` clears every mark, and marks reset on a
  collection switch.
- **Deleting several plasmids at once** — with rows marked Ⓧ, Delete
  opens one confirmation that leads with **how many** plasmids are
  about to go and names them, and — the number that actually matters —
  how many **sequences** would be held by no other entry afterwards.
  That count is per sequence, not per row, so two marked copies of the
  same construct report one sequence at risk rather than two. The whole
  delete is a single undo step: press **`u`** once and every plasmid
  returns to its original slot. Undo is all-or-nothing — if you switched
  collections in between, it restores nothing and says so, rather than
  bringing back part of the batch.

## File formats

| Format | Extensions | Import | Export | Preserves |
|---|---|---|---|---|
| GenBank | `.gb` / `.gbk` / `.genbank` | yes | yes | features, qualifiers, wrap topology, arrowless strand, display name (COMMENT marker) |
| EMBL | `.embl` | yes | yes | features, qualifiers, wrap topology, arrowless strand |
| CommercialSaaS `.dna` | `.dna` | yes | yes | features + colours + primers + construction history (round-trip) |
| FASTA | `.fa` / `.fasta` / `.fna` / `.ffn` / `.frn` / `.fas` / `.mpfa` / `.faa` | yes | yes (sequence only, wrapped at 60 columns) | sequence only |
| Sanger trace | `.ab1` / `.abi` | yes | — | base-called sequence + Phred quality |
| FASTQ multi-read | `.fastq` / `.fq` | yes (≤ 1000 reads per file) | — | one library entry per read |
| GFF3 | `.gff` / `.gff3` | yes (standalone via `##FASTA`; or apply features to loaded canvas via `apply-gff3`) | yes | wrap features as same-`ID=` split rows; `Is_circular=true` on region row |
| Plasmidsaurus zip | `.zip` | yes (Sequencing → Plasmidsaurus tab) | — | consensus + run-level QC + AB1 traces |
| Circular map image | `.png` / `.svg` | — | yes (File menu, bulk, or agent) | rendered circular map (picture, not re-importable) |

All file reads route through size-cap + symlink-refusal checks; all
file writes route through `_atomic_write_text` / `_atomic_write_bytes`
(tempfile + fsync + replace + symlink refusal). Bulk-import and
single-file Open share one dispatch table so the agent CLI, GUI Open,
and folder-import all accept the same set.

**Conformance.** Every text format SpliceCraft writes is checked against
its published specification — the NCBI GenBank flat-file layout and the
INSDC Feature Table definition, the ENA EMBL ID line, GFF3 1.26, RFC 4180
for CSV — by validators written from those specs rather than from the
serialiser (`tests/test_format_conformance.py`). In practice that means:
every line fits 80 columns and is 7-bit ASCII (non-ASCII text is
transliterated on the way out, while the library keeps exactly what you
typed); the LOCUS line keeps its classic column positions whatever your
plasmid is called; CDS translations wrap the way NCBI wraps them; GFF3
attribute tags are unique per line and an origin-spanning feature keeps
its arcs in 5'→3' order; FASTA wraps at 60 columns; and files are written
with LF line endings on every platform so a `.gb` from a Windows machine
is byte-identical to one from a Mac. Names that the GenBank LOCUS field
cannot hold (spaces, punctuation, anything long) are preserved in a
`SpliceCraft-name:` COMMENT marker, so an exported plasmid re-imports
under the name you gave it.

Passing our own validators is necessary and not sufficient, so every
export is also re-read by software written by someone else: an independent
GFF3 parser, an independent SnapGene `.dna` reader, and a different
Biopython release than the bundled one. That is how the `.dna` writer is
known to place every feature at the coordinates a third-party reader sees,
origin-spanning features included, and it is what a corpus of 27 files
written by other tools — NCBI nucleotide and protein records, European
Nucleotide Archive files, a 4.6-million-base genome — is checked against
on the way in. Protein records export too, with `aa` units and a wrapped
DBSOURCE line the way NCBI publishes them.

On the way in, GenBank / EMBL / FASTA / GFF3 files are read tolerantly: a
UTF-8 byte-order mark (what Windows editors add when they "save as
UTF-8"), a CP1252-encoded author name, and CRLF or CR-only line endings
all load rather than failing. Anything the parser had to guess at is
written to the log.

### Is that mismatch real? (basecall quality)

A Sanger read is unreliable at both ends — that is simply how the
chemistry runs out — so a disagreement there is the instrument guessing,
not a mutation in your plasmid. **Sequencing → Sanger (.ab1) → Check
against canvas** compares the trace against whatever plasmid is loaded
and answers on the read's own per-base quality scores: how many of the
differences the read actually stands behind.

The toast says it directly — *"1 real change, 1 in unreliable
basecalls"* — and the Verification Report gains a **Real** column
showing the same split (`2 (+35 noise)`). It also reports the usable
window after end-trimming, so a read that failed outright says so
instead of claiming a clean result over nothing.

A read with no recorded quality (an older instrument, a consensus FASTA)
shows an em-dash rather than a verdict: unknown and zero are different
answers. The threshold lives in **Settings → Sequencing → Min basecall
quality**, defaulting to Phred 20 (one wrong base in a hundred, the
Sanger convention); set it to 0 to trust every base equally and see
every mismatch as real.

## CRISPR guides

**CRISPR** on the menu (or `Alt+G`) scans whatever plasmid you have open for
guide sites and ranks them. Pick the nuclease — SpCas9 (NGG), SaCas9
(NNGRRT) or LbCas12a (TTTV, which sits on the *other* side of the spacer and
cuts with a stagger) — and optionally a bp window to aim at one exon rather
than the whole plasmid. Every guide comes back with its strand, the predicted
cut position, and the concerns that matter, named rather than scored: a
`TTTT` that will stop Pol III mid-transcript, a GC window outside 40–70%, a
homopolymer run, a spacer that folds back on itself.

**Two things it deliberately will not tell you.** There is no genome-wide
off-target prediction — that needs a genome index this program does not have,
and it will not pretend otherwise. **Check off-targets** searches the plasmid
in front of you and says how many base pairs it looked at, so "none found"
can never be mistaken for "specific". And there is no efficiency score:
published on-target models are regressions over large screens, and a number
that *looked* like one would be worse than the named flags you can check
yourself.

Pick a guide and **Cloning oligos** gives you the annealed pair for a Type IIS
guide vector, copied to the clipboard — with a 5' G added for U6 only when the
spacer does not already start with one, and said so, because that off-by-one is
how the pair is usually got wrong. **Annotate on map** drops the protospacer
onto the plasmid as a feature.

## Editing a residue, not a base

`Alt+Shift+G` changes an amino acid and lets the DNA follow. Pick a CDS, a
residue number and the new residue, and the resolved codon is shown before
anything happens — `CAT → TGG at 177, 178, 179` — because the mapping from
"residue 143" to three base pairs runs through the strand, the `/codon_start`
offset, any introns and the origin, and that is exactly where a hand-done edit
goes wrong. The codon is chosen from your active codon-usage table, with the
synonymous alternatives listed so you can pick one that avoids a site you care
about. A silent change is allowed and labelled as silent; introducing or
removing a stop is called out. `Ctrl+Z` undoes it.

## Ordering and running it

The primer export sheet now has a **plate** layout alongside the generic and
IDT-bulk ones: Plate / Well Position / Name / Sequence, filled row-major the way
vendor templates expect, rolling onto a second plate past 96 oligos instead of
quietly dropping the rest.

After **Run PCR** in the Simulator, **Thermocycler program** turns the result
into the block you type into the machine — annealing temperature from the primer
Tms the simulation already measured, extension time from the product length at
your polymerase's own rate, and every number labelled with the rule that
produced it. Pick Q5, Phusion or Taq. If no Tm is available it anneals at a
conventional 55 °C and *says so*, rather than presenting a guess as a
calculation.

## Do the reads agree?

One read agreeing with itself is not evidence. **Combine reads** in the
Verification Report rolls up every read on a plasmid and asks, for each
difference, how many reads that cover that base actually show it — and how many
looked and disagreed. A change two reads agree on is a mutation; one read
showing it while another covers the same base and disagrees is almost always a
miscall. It also reports what nobody read: a plasmid can score "100% identity"
over the half that was sequenced, and the largest unread stretch is named
outright.

## Experiments lab notebook

- **Projects layer** (`Menu → Experiments`) — named projects (e.g.
  "Yeast Y32 strain", "E. coli toolkit") each hold a list of
  experiment entries. Active project persists across sessions; first
  launch wraps the user's existing experiments into a default
  project.
- **Compose + Attachments** — per-entry markdown body (1 MB cap),
  plus an attachment grid. Win/Mac clipboard paste via
  `Pillow.ImageGrab.grabclipboard()` (Linux/WSL disabled by Pillow).
  Images preview inline; traces (`.ab1`), tabular exports
  (`.csv`/`.tsv`/`.txt`), datasheets (`.pdf`), sequence records
  (`.gb`/`.fasta`/`.fastq`/`.dna`) and `.xlsx`/`.zip` attach as files
  and list without a preview. Caps unchanged at 10 MB per file and
  100 MB per entry. Inserting a non-image into the body writes a plain
  link rather than an image tag.
- **Cross-refs** — `@<plasmid-id>` inlines a coloured chip linking
  to a library plasmid; `!<action-id>` references the curated
  `_EXPERIMENT_ACTIONS` catalog; `&<gel-id>` references a saved
  gel snapshot. Ctrl+G / double-click on any tag opens the
  referenced entity.
- **Spellcheck** (F7) — pyspellchecker-backed (pure-Python English
  wordlist) with markdown-aware masking. Custom dict per-user.
- **Filter + cross-project search** — the filter box above the entries
  list narrows the OPEN project as you type (`#tag` filters by tag and
  composes with the words), showing "3 of 41" so an over-narrow filter
  doesn't read as an empty notebook. `Ctrl+F` searches EVERY project:
  terms must all appear, in the title, tags or body, and each hit shows
  its project plus the matching line. Same two-tier split as the plasmid
  side (library panel filters the active collection, `LibrarySearchModal`
  goes cross-collection).
- **Backlinks** — `n` on a library row lists the notebook entries that
  mention that plasmid, and opens the one you pick (switching project
  first if it lives elsewhere). Queries both the entry id and the display
  name, since the `Plasmid ref` button writes the id while a hand-typed
  `@ref` is usually the name.
- **Protocol** — pulls a plasmid's recorded build steps into the entry as
  numbered markdown, from the same data the History viewer's Protocol
  pane shows. Only unambiguous operations become `!action` tags (Golden
  Gate, Gibson, PCR, mutagenesis) and only names that resolve in the
  library become `@refs`, so the inserted protocol never claims a bench
  step that wasn't recorded or plants a ref that goes nowhere.
  Transformation / miniprep / sequencing are absent by construction — the
  simulation never saw them.
- **Steps** — the forward half: pick bench actions (space to toggle) and
  get a heading per step with its `!action` tag, in catalog order.
  Headings and blank space only — no imposed "Cells: / Plates:" form.
- **Copy** — duplicates the selected entry as a starting point.
  Attachments are deliberately not carried over (they belong to the
  original entry's directory) and the app says so.
- **Export** (`Ctrl+E`) — one entry or the whole project, to Markdown or a
  standalone HTML page. Markdown keeps the body verbatim so it pastes
  back into a new entry with working refs, and adds a metadata header, a
  resolved References block and the attachment list around it. HTML
  inlines the stylesheet and embeds images up to 25 MB, falling back to
  local paths *and saying so* rather than producing a file that looks
  complete and breaks when sent.

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
- **User-overrides** — add (**New**), edit, and remove (**Delete**)
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

## BABS assistant

An in-app chat with a **local Ollama model** — no cloud, no API key,
nothing leaves the machine. Streaming markdown, `<think>` reasoning
hidden by default, a **context lifebar**, a copy-pasteable transcript
(`Ctrl+E` exports to markdown), and slash commands (`/help`, `/model`,
`/system`, `/temp`, `/reset`, `/retry`, `/agent`, `/autonomy`,
`/agentmodel`, `/agentprotocol`,
`/recall`, `/ingest`, `/learn`, `/forget`, `/remember`, `/memory`).

- **Agent mode** — Babs calls the same endpoints the `--agent` API
  exposes, and every change shows up live in the app. She is always
  told which plasmid is open. Default `ask` autonomy prompts before
  every write (a multi-step workflow lists all its changes and is
  approved once); `/autonomy auto` runs unattended, `readonly`
  forbids writes. Whole-library wipes are never reachable, and
  **physical robot motion always asks first, even in `auto`**. Once a
  turn has read anything from the internet (a web page, search results,
  a remote database record), **changes to your data ask first even in
  `auto`** — including saving a memory — because text from the web can
  try to steer her; the prompt says why it is asking. Agent
  turns run on a fast tool-capable model (qwen2.5:7b default) while
  ordinary chat stays on your chosen model — override with
  `/agentmodel <name>|chat|auto`. A model **without native tool
  calling** still works: Babs switches it to one JSON action per step,
  a shape Ollama enforces, so it can't invent a tool, and its answer
  still streams as it's written (`/agentprotocol auto|native|json` pins
  either mode). Without the Babs engine installed, the tools that need it
  (corpus recall, reference tables, memory) are left out rather than
  offered and left to fail. Web pages,
  search results and papers she reads are handed to the model as data,
  never as instructions. A long tool result is shortened with a note
  saying what was left out, and when a long task outgrows the model's
  context window the oldest tool results are set aside first, so the
  original request is never lost. If she runs out of tool steps she
  tells you what she did and what is left rather than stopping cold, and
  if she ends a reply by announcing a tool she never called ("let's go
  ahead and check it"), she is asked once to make the call. If
  a change she tried failed, the turn ends with a warning naming it; if
  you asked for a change and none ran at all, it ends with "Nothing was
  changed this turn" — both from what actually ran, however confident
  her answer sounds.
- **One-call plasmid tools** — *plasmid overview* (every feature and
  every enzyme that cuts once, with where each cut lands), *find a
  feature by name* (coordinates, sequence, protein, and the enzymes that
  cut once inside it; forgiving about spelling and Greek letters),
  *amplify a feature* (PCR primers checked to bind exactly where
  designed, across the origin and on either strand — and, asked to save
  them, saved to the primer library under your names), *digest* (exact
  fragment sizes for any set of enzymes on the open plasmid) and *add a
  feature* (positions exactly as you say them — "bases 100 to 200" —
  converted to SpliceCraft's coordinates in code). Also agent endpoints:
  `plasmid-overview`, `find-feature`, `amplify-feature`, `digest`.
- **Fast on a CPU** — the open plasmid's description rides in the
  conversation instead of Babs's instructions, so after an edit the
  model reuses its cached prompt instead of re-reading all of it
  (measured on a CPU-only laptop: 255 s → 14 s before she starts
  answering), and the model stays loaded for 30 minutes between
  messages rather than Ollama's default 5.
- **Test bench** — from a source checkout, `python scripts/babs_eval.py`
  runs ten everyday requests through the real agent against your local
  models and grades the answers, to compare models, protocols or prompt
  changes. Runs in a throwaway sandbox; slow on a CPU (minutes per task).
- **Corpus grounding** — with **Corpus** on, answers are grounded in
  the Babs research corpus with cited sources, folded into the turn
  inside a token budget scaled to the model's real context window.
  `/recall <query>` searches it directly; `/ingest <url>` stores an
  open-licence page permanently; `/learn <topic>` starts a focused
  crawl.
- **Learn** — a drift-resistant crawl of open scientific databases:
  seeds from open-access search, scores every candidate, and only
  follows citations out of strongly on-topic ones. Each topic gets an
  isolated corpus that can be folded into the main one afterwards.
  Also drivable over the agent API (`learn-start` / `-status` /
  `-results` / `-list` / `-merge`).
- **Index my library** — indexes your plasmids, primers, parts and
  notebook entries into her corpus. Rebuilt from scratch each run (no
  ghost records) and **private by construction**: every record flagged
  private, a `.no-egress` marker on the directory, skipped by the
  corpus-sync script, and never used to derive a web query.
- **Online lookups** — FPbase, UniProt, Europe PMC, NCBI/GenBank,
  Wikipedia, general web search and patents, including opening a
  result and reading the page text. Off until *Settings → "Allow Babs
  online database lookups"*; only the query string (or the URL being
  read) is ever sent, never a sequence. Web/patent search use Brave
  Search / PatentsView keys if configured, else keyless fallbacks.
- **Model manager** — models organise into collections like plasmids:
  mark rows to move or bulk-uninstall, search HuggingFace for GGUF
  models, and pull with a live bytes/speed progress bar.
- **Paper scraper** — runs the real [Babs](https://github.com/ATinyGreenCell/babs)
  background search (Europe PMC, OpenAlex, CORE, CGSpace, DOAJ) and
  re-index, with a jobs view and log tail. Shown only when the Babs
  repo is present (`~/babs`, or `$SPLICECRAFT_BABS_HOME`). Changing
  the *embedding* model warns that the corpus needs re-ingesting.
- **Persistent memory** — `/remember <fact>` survives across sessions
  as hand-editable markdown; `/memory` lists it, `/forget <slug>`
  removes one.

The tab is persistent — conversation, model choice and agent mode
survive switching away, and a model download keeps running while you
work elsewhere. Requires `ollama serve` locally (or
`$SPLICECRAFT_OLLAMA_HOST`); `splicecraft babs-setup` clones and
bootstraps the engine in one command.

## AUTOLAB (Opentrons OT-2)

A five-panel protocol designer for an OT-2 liquid handler.

- **Find Robots** scans the network (and USB) and lists what it finds,
  or auto-connects to the first.
- **Deck** — a top-down diagram of the robot drawn at the real slot
  aspect ratio, scaled and centred to the terminal, colour-filled by
  what each bay holds. Click a bay to place or clear labware.
- **Designer** — an ordered step sequence: transfer, distribute,
  consolidate, mix, delay, pause, comment.
- **Labware** — a library of custom labware defined from a grid form.
- **Library** — binds a deck plate to a plasmid collection so wells map
  to your plasmids; cherry-pick or replate by identity, or normalise
  DNA concentration to a target ng or ng/µL.
- **Calibrate** — deck / per-pipette offset / tip-length status, home
  the gantry, set per-slot offsets, and a **position check** that moves
  the gantry over each labware without ever actuating the plunger.

Designs save as named protocols in collections, compile to a real
Opentrons protocol, analyze on the robot's built-in simulate, and run
behind Arm + a clean analysis + health, calibration, attached-pipette
and door checks — with live Pause / Resume / Abort, a telemetry panel
(connection, light / door / motor / pipette state, progress bar with
ETA, fault banner), and Lights / Disengage Motors buttons. A build can
be logged straight to the Experiments notebook, cross-linked to the
plasmids it touched.

Scriptable via `ot2-compile` (a `transfers` list or multi-step `steps`
sequence), `ot2-analyze`, `ot2-status`, `ot2-calibration`, the gated
`ot2-run` / `ot2-position-check` / `ot2-home`, `ot2-run-control`,
`ot2-lights`, `ot2-disengage`, `ot2-normalize`, `ot2-plate-map`, plus
protocol and custom-labware library CRUD.

## Drive it from outside the GUI

See [Agent API](agent-api.md) and [CLI sidecar](cli.md) for the
details. In short:

- **Agent API** (`splicecraft --agent`) exposes a localhost JSON API
  with bearer-token auth, covering every GUI action external AI
  agents need. ~267 endpoints; symlink-guarded write paths;
  length/range/shape validation at the boundary.
- **`splicecraft-cli`** — stdlib-only sidecar (~50 ms cold start)
  that reads connection details from the running session's token
  file. Intended for Claude Code, Cursor, aider, hand-rolled
  scripts, or any external automation.
