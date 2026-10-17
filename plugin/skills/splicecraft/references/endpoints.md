# The endpoint surface, organised by job

291 endpoints as of SpliceCraft 1.3.11. **`splicecraft-cli tools --json` is
always authoritative** — this file helps you find the right name and warns you
about the sharp edges; the endpoint's own `doc_full` is the contract.

Naming convention, worth learning once: `list-X`/`get-X` read · `create-X`/
`update-X`/`delete-X` write · `simulate-X` is a read-only dry run whose write
twin takes the *same payload*.

---

## The core set

| # | Endpoint | Why it's core |
|---|---|---|
| 1 | `tools` | Discovery root — exact body shapes, one round-trip instead of N 400s |
| 2 | `status` | What's loaded, is it dirty (every write is dirty-guarded) |
| 3 | `plasmid-overview` | Best first look: name, length, topology, features, **unique cutters grouped by the feature they land in** |
| 4 | `search-library` | Cross-collection fuzzy find; turns a human name into a loadable id (`puc` matches `pUC19`) |
| 5 | `load-entry` | Open a library plasmid; resolves across collections |
| 6 | `features` | The `idx` handles every annotation write takes |
| 7 | `find-feature` | Name → coordinates + sequence + translation + unique cutters inside; forgiving of case/spacing/hyphens/Greek (`lacZ alpha` finds `LacZ-α`) |
| 8 | `get-sequence` | Pull bases from any span, wrap-aware (`bottom` for reverse complement) |
| 9 | `list-restriction-sites` | "What cuts, where, once?" with real cut offsets |
| 10 | `add-features` | Batch annotation under one lock — **not** N × `add-feature` |
| 11 | `add-current-to-library` | The explicit persist step. Without it the work evaporates |
| 12 | `simulate-golden-gate` → `golden-gate-assemble` | Dominant modern cloning path |
| 13 | `simulate-traditional-cloning` → `traditional-clone` | Restriction cloning, with the self-ligation off-target report |
| 14 | `design-primers` / `amplify-feature` | Primer design with binding verified |
| 15 | `check-primer-pair` / `check-primer` | Is the pair sound — Tm gap, hairpins, dimers — where does each bind, and what does it produce |
| 16 | `optimize-protein` | Codon optimization with host-hazard scrubbing |
| 17 | `verify-against-reads` | The other end: reads came back, does the clone match |
| 18 | `digest` | Overhang-aware digest of a **raw sequence you supply** — QC any junction you just designed |
| 19 | `create-experiment` | Write the bench record |
| 20 | `reference-lookup` | Curated local registries for named reagents — prefer over recall |

---

## Session, state, introspection

`tools` `status` `plasmid-overview` `get-settings` `set-setting` `undo` `redo`
`discard-changes` `restart` `shutdown` `capture-snapshot`

- `discard-changes` is the recovery path for a stuck-dirty canvas. After it,
  `status.dirty` is false and dirty-guarded mutations proceed without `force`.
- `undo` returning `undone: false` means the stack was empty — not an error.
- `set-setting` does not refresh live in-memory mirrors; most toggles take effect
  next session.
- `shutdown` / `restart` are **headless-only** (both 409 in an interactive
  `--agent` session). `restart` drops the canvas and undo stack, and is how a
  daemon picks up a `splicecraft update`.

## Getting DNA in and out

`fetch` `genbank-search` `load-entry` `load-file` `new-plasmid` `save`
`add-current-to-library` `bulk-import-folder` `bulk-export-collection`
`export-genbank` `export-fasta` `export-embl` `export-gff`
`export-commercialsaas` `export-map-image` `export-migrate-archive`
`apply-gff3` `import-feature-preset` `ingest-url` `read-url`

- **The fetch → save trap.** `fetch`, `load-file` and `new-plasmid` load for
  inspection/editing only; `saved` is always `false` to make that explicit. A
  bare `save` on such a record returns `409` — creating a library entry needs
  `{"create": true}`, or use `add-current-to-library`.
- `load-entry` resolves across collections (active first). Ambiguity → `409`
  listing holders; pass `collection` to pin.
- `load-file` / `add-current-to-library` take `{id, name}` overrides — essential
  in batch loops (see SKILL.md §3). `load-file` is capped at 50 MB, `force`
  overrides.
- `apply-gff3` requires a GFF3 with **no** `##FASTA` section (that case is
  `load-file`).
- `export-migrate-archive` has no inverse over the API by design — import
  replaces all data, the same risk class as a master wipe.

## Looking at a plasmid

`plasmid-overview` `features`(=`list-features`) `get-feature` `find-feature`
`find-sequence` `get-sequence` `span-contains` `list-restriction-sites` `digest`
`list-enzymes` `find-orfs` `get-history` `get-plasmid-protocol` `diff-plasmid`
`simulate-gel` `replace-sequence` `get-feature-colors` `set-feature-color`
`set-origin` `reverse-complement`

- `set-origin` re-cuts a CIRCULAR plasmid so a chosen base (or a named
  feature's start) becomes base 1; `reverse-complement` flips the whole record.
  Both re-frame every feature, keep spliced exons and origin wraps intact, are
  undoable with `undo`, and CLEAR stored alignments — every coordinate moved, so
  a band drawn against the old frame now points at the wrong bases. `set-origin`
  **422s on a linear record**: rotating it would join two ends that aren't
  joined. Flipping a linear record is fine.
- `list-restriction-sites` respects the active enzyme collection; the full
  catalogue returns hundreds of hits (309 on pUC19). Returns `cut_bp` /
  `bottom_cut_bp` so fragment sizes need no private offset table.
- **"What cuts it once?" is `{"unique_only": true}`** (38 on pUC19), not a scan you
  filter yourself. Add `respect_active_collection: false` to ask the whole
  catalogue regardless of the user's bench collection. An enzyme with one site
  *inside your region of interest* may still cut elsewhere on the plasmid and so
  won't linearise it — `unique_only` is what separates those.
- **An unknown enzyme name is a `400`, never an empty result** — the failure it
  replaced was wrong-answer-shaped, a typo returning `count: 0` that reads as "no
  sites here" rather than "I never looked".
- Isoschizomers and commercial synonyms resolve by the name you asked for
  (`Esp3I`=`BsaI`-family `BsmBI`, `Eco31I`=`BsaI`, `LguI`=`SapI`, `AarI`=`PaqCI`).
  **Neoschizomers stay apart**: XmaI (`C^CCGGG`) and SmaI (`CCC^GGG`) read the
  same six bases and leave different ends.
- `digest` works on a raw sequence in the body, not the canvas — use it to verify
  a junction or dropout overhang you just designed.
- `span-contains` exists because `a <= x < b` reports every base of an
  origin-spanning marker or operon as *outside* it.
- `find-orfs` takes `min_aa` in **amino acids**; `min_length`/`min_bp` are
  rejected 400 so a unit mix-up can't silently return a default-length result.
- `replace-sequence` drops features consumed entirely by the replace.
- `get-plasmid-protocol` renders construction history as paste-ready bench
  markdown; `recover-history-from-dna` rebuilds history from a `.dna` file
  (dry-run by default).

## Annotation

`add-feature` `add-features` `update-feature` `delete-feature`
`transfer-annotations` `annotate-from-presets` `apply-gff3`
`list-feature-presets` `get-feature-preset` `import-feature-preset`
`list-feature-library` `get-feature-library` `create-feature-library`
`update-feature-library` `delete-feature-library`

- `add-features` batches under one UI-thread lock and one stale-check. Prefer it.
- `add-feature`'s `strand` is the **arrow type**: `1` forward, `-1` reverse, `0`
  arrowless, `2` double-stranded. A `CDS` must span whole codons or it's rejected.
- `update-feature` is a partial update by `idx`; only supplied fields change.
- Both `transfer-annotations` and `annotate-from-presets` are **dry-run by
  default** (`{"apply": true}` to commit). Preset matching is **exact**, so an
  allelic variant won't be called — use `blast` for near-identical hits.
- `apply-gff3` is **not** dry-run: it writes on the first call and returns
  `{ok, features_added}` with no `applied` field. See the commit-key table in
  SKILL.md §3 — the dry-run endpoints do not all take the same key.
- `transfer-annotations`: **read `skipped`** (see SKILL.md §4).

## Library management

`list-library` `search-library` `index-library` `list-collections`
`get-active-collection` `set-active-collection` `create-collection`
`rename-collection` `delete-collection` `copy-plasmid` `copy-plasmids`
`move-plasmid` `move-plasmids` `rename-plasmid` `delete-from-library`
`restore-deleted-plasmid` `set-plasmid-status` `list-plasmid-statuses`

- `rename-plasmid`, `set-plasmid-status`, `delete-from-library` and
  `diff-plasmid` act on the **active collection only**.
- `move-plasmid`/`move-plasmids` are atomic under one lock — use them instead of
  copy + delete, which leaves a window where the plasmid exists nowhere.
- `copy-plasmid`'s `to` must **not** be the active collection (copying into the
  active set goes through `add-current-to-library`).
- `delete-from-library` matches by **name**, so one call can remove several
  entries. `restore-deleted-plasmid` undoes the most recent deletion —
  including your own — putting each entry back in the SLOT it came from. A
  bulk delete restores as ONE unit: all of it or none, because restoring part
  of a batch and reporting success reads as complete recovery. It 409s if the
  active collection changed since the delete (the records stay recoverable
  until you switch back).
- **There is no agent endpoint for the master wipe, deliberately.** Deleting
  and undoing individual plasmids is yours; erasing the whole data directory
  is the user's, at the keyboard. Don't offer to do it.

## Primers

`design-primers` `amplify-feature` `design-mutagenesis` `check-primer`
`check-primer-pair` `list-tm-profiles`
`check-primer-duplicates` `simulate-pcr` `save-amplicon` `pcr-program` `list-polymerases`
`create-primer` `get-primer` `update-primer` `delete-primer` `delete-primers`
`move-primer` `list-primers` `export-primers`
`create-primer-collection` `rename-primer-collection`
`delete-primer-collection` `list-primer-collections`
`get-active-primer-collection` `set-active-primer-collection`
`update-stored-primer` `import-primers`

- `design-primers` has **four** modes: `detection` (diagnostic amplicon inside
  the region), `cloning` (needs `site_5`+`site_3`, appends RE tails), `generic`
  (binding-only), and `allele_specific` (SNP / genotyping — it pins the primer's
  3' TERMINUS on the discriminating base, where the other modes slide the window;
  an A/T discriminating base is reported as precluding a GC clamp rather than
  silently scored down). For Type IIS grammar overhangs use `design-gb-part`.
  Each returns `qc` — the pair check below — so read `qc.warnings` too.
- `check-primer-pair` is the QC for a pair you were GIVEN (a colleague's, a
  paper's, a supplier's): pass the oligos as ordered and, for a tailed primer,
  its 3' annealing arm as `forward_binding` / `reverse_binding` — otherwise the
  tail is scored as mismatches and the Tm is the whole oligo's. Report
  `warnings` verbatim, and `pair.product.product_bp` (tails included — the band
  on the gel) rather than `template_bp`. `single_product: null` means the
  search was capped, not that there is one product.
- A Tm is only comparable under the same conditions. If the user quotes a Tm
  from Benchling, pass `tm_profile: "benchling_compatible"` (it reads ~4-5 °C
  lower than SpliceCraft's default); `list-tm-profiles` gives the conditions.
  Never convert one tool's Tm to another's by hand.
- Always surface `binds_as_designed` and `other_sites` — a primer that displays
  somewhere other than where it anneals is a wrong experiment.
- `save-amplicon` files one of `simulate-pcr`'s products as a linear library
  entry (it re-runs the same payload, so preview and save cannot disagree).
  `amplicon_index` 0 is the LONGEST product — pass the index you saw, don't
  assume there is only one. `name` is required.
- `update-stored-primer` renames a STORED primer, moves it through
  Designed→Ordered→Validated, or sets the `$` cart flag; `export-primers
  {cart: true}` then writes the cross-collection cart. Note `update-primer` is
  a DIFFERENT endpoint — it edits a `primer_bind` feature on the canvas.
- `import-primers` is the inverse of `export-primers {format: "generic"}`.
  Invalid rows are skipped and reported in `skipped`, never filed as broken
  primers — read it, a partial import is normal for a hand-edited sheet.
- `simulate-pcr`'s binding model is **exact-match**: no mismatch tolerance, no
  Tm-aware annealing. Don't present it as a thermodynamic prediction.
- `pcr-program` with no Tms anneals at a conventional 55 °C and **says so** in
  `notes` rather than presenting an invented number as derived.
- `export-primers` emits an order sheet (`generic` / `idt` / `plate`).
- `update-primer` requires the feature at `idx` to be a `primer_bind`; coordinate
  edits are excluded — relocate via delete + re-add.
- `delete-primer(s)` refuses (409) to delete from a named collection while it's
  **active** — switch away first. (Same rule for `delete-part`/bins and
  `delete-experiment`/projects; last-one guards apply too.)

## Cloning, assembly, and the Golden Braid / MoClo system

Preview → commit: `simulate-golden-gate`→`golden-gate-assemble` ·
`simulate-gibson`→`gibson-assemble` ·
`simulate-traditional-cloning`→`traditional-clone`

Pipeline: `design-gb-part` `design-synthesis-fragment` `domesticate-part`
`domesticate-parts` `make-l0-part-from-fragment` `assemble-into-entry-vector`
`classify-part` `auto-detect-entry-vectors` `scrub-plasmid` `lint-synthesis`
`assemble-operon` `simulate-pcr` `pcr-program`

Entry vectors / grammars / parts bins: `list-entry-vectors` `get-entry-vector`
`set-entry-vector` `clear-entry-vectors-for-grammar` `list-grammars`
`get-grammar` `create-grammar` `update-grammar` `delete-grammar` `list-parts`
`get-part` `create-part` `update-part` `delete-part` `move-part`
`list-parts-bins` `create-parts-bin` `rename-parts-bin` `delete-parts-bin`
`set-active-parts-bin` `get-active-parts-bin`

- **Nothing is ever picked by fragment size.** `simulate-traditional-cloning`
  tries every (vector fragment × insert fragment × orientation) and returns all
  compatible products; a stacked-TU insert can outgrow its carrier vector, which
  is exactly why size is banned as a heuristic. Choose by what the fragment
  carries (origin, resistance marker).
- Read `self_ligation` on your chosen product: `vector_self_closes` (expect
  background colonies), `double_insert` (tandem/inverted), `orientation_ambiguous`,
  and `clean` with a `clean_note` when nothing off-target closes.
- A circular vector or insert must cut **exactly twice**; an ambiguous ≥3-cut
  digest is refused with the cut list. Type IIS enzymes are refused for the
  insert — use the Golden Braid path.
- `simulate-golden-gate`: parts may be in **any order**, overhangs set the order;
  each part must release exactly one fragment. The vector **may** be cut more
  than twice and all pieces go back in. When nothing closes, the error lists
  every fragment's two overhangs so the dangling one is visible.
- Pass `carry_annotations: true` or the assembled product saves **feature-bare**;
  then check `carried` / `carry_warnings`.
- `domesticate-part` and `assemble-into-entry-vector` **fail loud (422)** rather
  than emit a stub backbone — a stub L0 is a wrong construct an automated caller
  can't distinguish from a real clone. Treat the 422 as information.
- `lint-synthesis` is the pre-order pre-flight (internal Type IIS sites, GC
  extremes, homopolymers, repeats, ORF check, 0–100 score).
- `scrub-plasmid` plans a clone-free, protein-preserving site scrub; an
  unrecognised enzyme name is a 400, never a scrub that quietly removes nothing.

## Mutagenesis and residue editing

`design-mutagenesis` `plan-residue-edit` `edit-residue` `replace-sequence`
`undo` `redo` `discard-changes`

`plan-residue-edit` resolves the three plasmid bp a codon occupies **through
strand, `/codon_start`, introns and the origin** — the mapping you'd otherwise do
by hand, which is where frame errors come from. It lists synonymous
`alternatives`. `edit-residue` applies it as a 3-for-3 substitution, so nothing
shifts, every feature keeps its coordinates, and it's undoable.
`design-mutagenesis` gives SOE-PCR outer/inner pairs for e.g. `"W140F"`.

## Codon optimization and protein

`optimize-protein` `analyse-cds` `list-codon-tables` `add-codon-table`
`delete-codon-table` `set-active-codon-table` `get-active-codon-table`
`list-protein-motifs` `set-protein-motif` `delete-protein-motif`
`taxon-search`
`list-protein-collections` `get-protein-collection`
`create-protein-collection` `rename-protein-collection`
`delete-protein-collection` `save-protein` `delete-protein`

- `optimize-protein`: `mode` = `frequency` (default) or `max_cai`, plus
  `hazard_hosts` (plant/mammalian/bacterial), `forbidden_motifs`, GC banding and
  `avoid_repeats_with` to diversify a panel. GC-window and diversify passes are
  capped at 1500 codons. Check `remaining_motifs`/`remaining_repeats`/`fixes`.
- Only *E. coli* K12 ships as a table; `add-codon-table` fetches from Kazusa by
  taxid, builds from an NCBI genome, or imports a local CDS FASTA offline.
- `analyse-cds` is the read half: CAI, GC3, `rare_pct`, `tandem_rare_runs`,
  `worst_window`, 5′ `ramp`.

## Large chromosomes

`plan-chromosome-segments` `list-assembly-strategies`

- A 1 Gbp target can't be ordered as one piece. `plan-chromosome-segments`
  cuts it into orderable segments and **proves the junctions are
  unambiguous** — read `n_ambiguous_junctions` and `ok` before the segment
  table. Two junctions sharing an overlap (or one being another's reverse
  complement) give an assembly that is chemically happy and structurally
  wrong: segments swap or insert inverted, and it surfaces after sequencing.
- `strategy` is REQUIRED — the chemistry sets the overlap size (Gibson ~30 bp,
  yeast HR ~120 bp, Golden Gate 4 nt), so there is no sensible default. See
  `list-assembly-strategies`.
- **Boundaries are nudged, not just taken.** Each cut moves within `nudge_bp`
  to the nearest clean junction; Golden Gate searches widest because 4 nt
  leaves only 256 overhangs and arbitrary cuts collide constantly. Segment
  lengths therefore vary. `nudge_bp: 0` takes the arithmetic cuts and shows
  you the collisions a fixed layout would have had.
- Forbidden sites come from the STRATEGY: only Golden Gate consults an
  internal Type IIS site. Pass `forbidden_enzymes` to override.
- `offtarget_checked: false` means that search was over budget and SKIPPED —
  not that nothing was found. Junction uniqueness is always checked.
- Pass `name`/`id` to read a library entry; a chromosome won't fit in the
  1 MiB body. `include_sequences` is off by default (a 2,000-segment plan
  carrying every base is tens of MB).
- The MAP aggregates above ~2,000 features (`map_lod_feature_threshold`):
  rotation went from 40 s to 8 ms at 200,000 features. A 76-column braille map
  has ~16,000 dots, so individual features can't be drawn at that scale —
  zoom in to resolve them.

## Sequence analysis

`blast` `blast-online` `hmmscan` `hmmer-web` `multi-align` `diff-plasmid`
`find-orfs` `map-transcription` `predict-transcript` `scan-promoters`
`scan-terminators` `rbs-strength` `design-rbs` `assemble-operon` `fold-rna`
`cofold-rna`
`add-hmm-database` `get-hmm-database` `list-hmm-databases`
`download-hmm-database` `delete-hmm-database`
`get-active-hmm-database` `set-active-hmm-database`

- `blast` is in-process against the user's own collections and **never leaves the
  machine**. A whole large library can exceed the DB build cap and 400 — pass
  `collections` to scope it. `blast-online`/`hmmer-web` ship data out and are
  403-gated (SKILL.md §6).
- `map-transcription` walks the circle from every promoter to every CDS: "what
  *else* transcribes this gene?" Checking terminators downstream of a cassette —
  the natural move — cannot see a read-through, because the read-through arrives
  from behind. `silent` vs `silent_in_host` are different claims and the gap is
  the point: a T7-only cassette is silent in a strain with no T7 RNAP *unless*
  something reads through into it.
- `predict-transcript` gives the mature mRNA of one unit plus Kozak and uORFs,
  separating real uORFs from ATGs `removed_by_splicing`. Introns come from
  **annotation only**; an unavailable model reports `skipped`, never an empty
  list that would read as clean.
- `scan-promoters`: `score` **ranks, it does not predict** — no trained model,
  not a transcription rate. `rbs-strength`'s `rel_strength` is relative-only.
- `scan-terminators`: one record per hairpin (a scan reporting every register
  separately turns one hairpin into twenty "terminators").
- `design-rbs` runs ~80 fold/cofold evaluations and blocks for a few seconds.
- `diff-plasmid` is rotation- and RC-aware and reports `inverted_segments`;
  `multi-align` is the batch form (≤20 targets).

## CRISPR

`list-cas-variants` `design-guides` `score-guide` `guide-offtargets`
`guide-cloning-oligos`

Variants differ in PAM geometry (SpCas9 NGG 3′, SaCas9 NNGRRT 3′, LbCas12a TTTV
5′ with a staggered cut). `score-guide` returns **named mechanistic flags**, not
a fabricated efficiency number. `offtarget_searched_bp` is the honest scope of
any off-target claim; omit the search and it's `0`, so silence can't be read as
"clean". No genome-wide prediction exists here and none is implied.
`guide-cloning-oligos` builds the annealed pair — writing the bottom oligo as
`AAAC + rc(guide)` is the classic error.

## Sequencing reads

`verify-against-reads` `read-consensus` `analyse-read-heterogeneity`
`read-trace` `store-alignment` `list-alignments` `update-alignment`
`delete-alignment`
`plasmidsaurus-items`(=`list-plasmidsaurus-items`) `list-plasmidsaurus-members`
`download-plasmidsaurus` `align-plasmidsaurus-zip`

- See SKILL.md §4 for the `identity_pct` / `uncovered_spans` traps.
- `read-trace` base-calls a Sanger `.ab1`/`.abi` and returns `bases` +
  `phred_quality` — the input `verify-against-reads`' `read_quality` wants.
  Feed them straight across to get the trace-quality split: a mismatch under a
  low-quality basecall is an instrument artefact, not a mutation.
- `store-alignment` aligns a read against a library plasmid and PERSISTS it,
  which is what `read-consensus` reads. **You name the read and the target;
  the engine computes the alignment** — there is deliberately no way to hand
  in gapped strings, because a caller-built alignment is a caller-built
  verdict. Re-storing the same read is idempotent. `list-/update-/delete-
  alignment` manage the stored rows (a `visible: false` row stays on disk and
  still counts in `read-consensus`).
- **Judge a read on `covered_identity_pct`, never raw `identity_pct`.** The
  raw figure counts the plasmid the read didn't reach as gaps, so a PERFECT
  800 bp Sanger read of a 5 kb plasmid scores ~16 %. Read it with
  `covered_pct`: "perfect over the third it covered" is the honest sentence.
- `analyse-read-heterogeneity` answers "is this culture one thing?" —
  `clonal` / `minor_variants` / `mixed`. A Plasmidsaurus consensus is depth 1, so
  `read-consensus` cannot see a heterogeneous population; it consenses to wild
  type and reads as clean.
- Plasmidsaurus: **pick by `has_results`, not `status`** — ~40% of a real listing
  is shipping-label orders marked `complete` with no results zip, and they sort
  to the top. `download-plasmidsaurus` runs synchronously (download → parse →
  save), so size your timeout. A `429` means the 10 req/min budget is spent.

## Lab notebook, gels, protocols

`create-experiment` `get-experiment` `update-experiment` `delete-experiment`
`duplicate-experiment` `move-experiment` `list-experiments` `search-experiments`
`experiment-backlinks` `export-experiment` `attach-experiment-image`
`experiment-step-template` `list-experiment-actions`
`create-experiment-project` `rename-experiment-project`
`delete-experiment-project` `list-experiment-projects`
`get-active-experiment-project` `set-active-experiment-project`
`list-gels` `get-gel` `create-gel` `update-gel` `delete-gel` `simulate-gel`
`save-protocol` `get-protocol` `list-protocols` `delete-protocol`
`create-protocol-collection` `list-protocol-collections`
`delete-protocol-collection`

Entry bodies carry cross-reference sigils `@plasmid`, `!action`, `&gel`,
extracted automatically. `search-experiments` requires **all** query terms and
**all** named tags. `experiment-backlinks` should be passed both the id and the
display name (a button writes the id, a human types the name).
`experiment-step-template` is the steps you're *about* to do;
`get-plasmid-protocol` is what the program *did*. `simulate-gel` renders an
in-silico agarose gel from ladder/plasmid/digest/pcr lanes.

## Robotics — Opentrons OT-2

`ot2-status` `ot2-calibration` `ot2-compile` `ot2-analyze` `ot2-run`
`ot2-run-control` `ot2-position-check` `ot2-home` `ot2-disengage` `ot2-lights`
`ot2-normalize` `ot2-plate-map`
`save-custom-labware` `get-custom-labware` `list-custom-labware`
`delete-custom-labware` `create-labware-collection`
`list-labware-collections` `delete-labware-collection`

Canonical order: `ot2-plate-map`/`ot2-normalize` (plan) → `ot2-compile` (offline
validate) → `ot2-analyze` (the robot's own simulate, no motion) →
`ot2-position-check` (motion, no liquid) → `ot2-run` → `ot2-status`.

`ot2-run` and `ot2-position-check` are **gated physical actuation**: both need
`{"confirm": true}` *and* a passing analysis, and refuse on an already-faulted
robot. Without `confirm` you get an analysis-only result and no motion.
`ot2-normalize` computes per-sample volumes to equalise DNA mass — pure compute,
no robot needed. "A plate is a collection": `ot2-plate-map` maps a plasmid
collection onto labware wells.

## Enzymes

`list-enzymes` `list-custom-enzymes` `get-custom-enzyme` `create-custom-enzyme`
`update-custom-enzyme` `delete-custom-enzyme`
`create-enzyme-collection` `update-enzyme-collection`
`delete-enzyme-collection` `get-enzyme-collection`
`list-enzyme-collections` `get-active-enzyme-collection`
`set-active-enzyme-collection`

Enzyme collections scope a restriction scan to what's actually on the bench.

## Knowledge, literature, reference lookup

`reference-lookup` `recall-knowledge` `index-library` `ingest-url`
`remember-fact` `memory-list` `memory-forget` `learn-start` `learn-status`
`learn-results` `learn-list` `learn-merge` `literature-search` `patent-search`
`web-search` `wikipedia-search` `read-url` `uniprot-search` `fpbase-search`

- `reference-lookup` is **local and instant**: curated registries of antibiotics
  with working concentrations, promoters, terminators, origins, Type IIS enzymes,
  strains and genotypes, binary/CRISPR vectors, base and prime editors, media,
  PGRs, recipes. Prefer it over your own recall for any named reagent.
- `recall-knowledge` returns ranked **cited passages** from the local corpus, not
  a generated answer. `index-library` indexes the user's own plasmids, primers,
  parts and notebook into a private local pack, so recall can answer "which of my
  plasmids carries a KanR and a T7 promoter?"
- **Every search in this cluster is egress-gated (403)** behind
  `allow_online_lookups` ("Allow Babs online database lookups"):
  `fpbase-search`, `uniprot-search`, `literature-search`, `genbank-search`,
  `wikipedia-search`, `web-search`, `patent-search`, `read-url`, plus
  `learn-start`, `ingest-url` and `recall-knowledge {web:true}`. That is a
  DIFFERENT setting from `allow_online_search`, which gates only `blast-online`
  and `hmmer-web` — see SKILL.md §7. `reference-lookup` and a plain
  `recall-knowledge` are local and ungated.

## Backup and maintenance

`list-backups` `restore-backup` `list-pre-update-snapshots`
`restore-pre-update-snapshot` `export-migrate-archive`

`list-backups` with no `label` summarises every user-data file across the
four-layer safety net (`legacy_bak`, `rotating_bak`, `snapshot`,
`lost_entries`). `restore-backup`'s `source_path` must come from a prior
`list-backups` response, and takes a fresh backup of current state first so the
pre-restore data isn't lost. Only ever on explicit user request.
