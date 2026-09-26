# Agent API

`splicecraft --agent` (alias `--agent-api`) exposes a localhost JSON
HTTP API with bearer-token auth, covering every GUI action external
AI agents need. The server binds `127.0.0.1` only and rejects any
request whose `Host` header isn't loopback (`127.0.0.1` / `localhost`
/ `::1`) — a DNS-rebinding defense, so a browser page can't reach the
API even to enumerate endpoints.

## Why it exists

Local AI coding agents (Claude Code, Cursor, aider, hand-rolled
scripts) work best when they can *do* things in the user's existing
environment, not just generate text. SpliceCraft's agent API lets an
agent drive the running GUI session through the side-door without
leaving its terminal.

## Quick start

```bash
splicecraft --agent                  # default port 6701
splicecraft --agent --agent-port 6800  # alternative port
splicecraft --headless               # agent API only, NO terminal UI (CI / no-pty)
splicecraft --agent --read-only      # attach ALONGSIDE a running GUI (reads only)
```

### Attaching while the GUI is open (`--read-only`)

The data dir is single-instance-locked for cache coherence, so an ordinary
`--agent` launch **refuses to start** while the GUI holds it. That used to mean
an agent could not look at the library at all during a session — and the way
round it was hand-parsing `collections.json` and `plasmid_blobs/*.gb`, which is
exactly what this API exists to prevent.

`--read-only` attaches without taking the lock. It serves every endpoint
`tools` reports as `write: false` and refuses every `write: true` one with
**409**, naming the PID that holds the lock so you know what to quit:

```json
{"error": "endpoint 'add-current-to-library' writes, and this SpliceCraft is
  attached --read-only while PID 4242 holds the data-dir lock. …",
 "read_only": true, "endpoint": "add-current-to-library", "holder_pid": 4242}
```

It implies `--headless` (a second TUI over one data dir is what the lock
prevents, and it would render caches that go stale the moment the GUI saves),
and it publishes its port + token to **`agent_token.readonly`**, not
`agent_token` — a guest must never repoint, or on quit delete, the host's
side-door. `splicecraft-cli` prefers the host token and falls back to the
guest, so `splicecraft-cli status` finds a read-only session with no flags;
set `SPLICECRAFT_READ_ONLY=1` to pin the guest when both exist. If the default
port is busy it walks the next few (a pinned `--agent-port` is never moved).

Refusing writes is the *legible* half of the guarantee. The half that actually
protects the data is that a read-only process never arms the save-authorisation
chokepoint at all, so any `_save_*` reaching it raises — and for the same
reason the launch-time collection / primer / parts-bin / project migrations are
skipped rather than attempted. The process holding the lock owns those mirrors;
the guest reads what the host wrote. Attaching leaves every data file
byte-identical.

`--headless` (alias `--agent-headless`, or env `SPLICECRAFT_HEADLESS=1`)
runs the JSON API under Textual's headless driver — no pty required, so
it backgrounds cleanly in CI / agent contexts (no `script -qfc` wrapper).
It implies `--agent`, prints a `SpliceCraft agent API ready on …` line to
stdout once the socket is accepting connections, and serves the same
endpoints as the UI launch. For readiness, poll the unauthenticated
`GET /healthz` (returns `{ok, status:"ready", version, headless}`)
instead of racing the token-file write + first real call.

A long-running daemon keeps serving the code it launched with, so after a
`splicecraft update` in another terminal it silently runs stale. `status`
surfaces this: it reports `installed_version` (what's on disk now) next to
the running `version`, plus a `stale` boolean (`installed_version` is `null`
when it can't be read, e.g. running from a source checkout). Every
authenticated response also carries a `_stale` warning when the running
version is behind, so you notice without polling `status` (the pre-auth
`/healthz` and `/tools` probes are served ahead of that stamp, so they stay
bare). A **headless** daemon can `POST /restart`
itself — it re-execs and the API returns on the same port (poll `/healthz`);
a GUI `--agent` session refuses (it would lose its live view — restart it
manually).

To stop a headless daemon, `POST /shutdown` rather than killing it. The app
quits through its normal teardown, so pending collection-sync and settings
writes are flushed, in-flight workers are drained, and the data-dir lock is
released; a `SIGKILL` skips all of that. The reply returns immediately and
the process exits about half a second later, so poll the returned `pid` (or
the port) to confirm — don't expect the socket to outlive the response. It
refuses with `409` when the canvas has unsaved edits (pass `{"force": true}`
to quit anyway), and, like `restart`, refuses on a GUI `--agent` session,
which has a keyboard and its own quit flow. Note that `--agent` already
auto-headlesses when it has no controlling TTY, so a daemon launched from
an agent context is eligible.

The server writes a token file at
`<DATA_DIR>/agent_token` containing the port + bearer token on the
first two lines. Hand the token to any client that needs to call the
API; the [CLI sidecar](cli.md) reads this file automatically.

```bash
# manual cURL
TOKEN=$(tail -1 ~/.local/share/splicecraft/agent_token)
PORT=$(head -1 ~/.local/share/splicecraft/agent_token)
curl -s -H "Authorization: Bearer $TOKEN" \
     -X POST -H "Content-Type: application/json" \
     -d '{"accession":"L09137"}' \
     http://127.0.0.1:$PORT/fetch
```

## Endpoint inventory

~267 endpoints across:

- **Records** — `new-plasmid` (create from a raw sequence, the Ctrl+N
  flow), get / set sequence, add / update / delete features (with the full
  arrow type — forward ▶ / reverse ◀ / arrowless / double-stranded ◀▶ —
  plus colour and arbitrary GenBank qualifiers, matching the Insert/Edit
  Feature dialog. **Colour takes whatever `list-features` reports** — hex
  (`#1f77b4`) or the terminal-palette form a GUI-drawn feature stores
  (`color(160)`), normalised to hex on the way in, so a read-features /
  write-features loop round-trips instead of refusing its own reader's
  output; `add-features` inserts MANY at once under one lock + one
  dirty-check; `name` is accepted as a synonym of `label` on add-feature and
  update-feature), list features (`features`, alias `list-features`), find ORFs
  (length cutoff in AMINO ACIDS —
  `min_aa`; `min_length`/`min_bp` are rejected so a bp-vs-aa mix-up can't
  silently return a default-length result. Each ORF carries `start`, `end`,
  `strand`, `length_aa`, `aa_seq`, plus `nt_len` — the exact coding length in
  bases INCLUDING the stop codon when there is one — `has_stop`, and
  `exceeds_one_lap`. **Take length from
  `nt_len` / `length_aa`, never from `end - start`:** on a circular plasmid a
  frame can run a full lap or more without an in-frame stop, and no
  start/end pair on an n-bp circle can express that. Such an ORF sets
  `exceeds_one_lap: true` and has its span pinned to the near-full circle.
  A frame with NO in-frame stop anywhere on the molecule sets
  `has_stop: false`; for it every codon is a residue, so `nt_len` is
  `length_aa * 3` rather than `(length_aa + 1) * 3`),
  `undo` / `redo` the last edit,
  `discard-changes` (revert the canvas to its library-stored copy / clear a
  stuck-dirty flag so the next load / new-plasmid proceeds without `force`),
  transfer annotations, apply GFF3 features to the loaded record
  (`apply-gff3`).
- **Workflows (one call each)** — `plasmid-overview` (the loaded plasmid's
  name, length, topology, every feature, and every enzyme that cuts it exactly
  once with the feature(s) that cut lands in); `find-feature` `{name, type?,
  max_seq?}` (a feature by NAME — forgiving of case, spacing, hyphens and Greek
  letters, best match first — with its coordinates both ways: 0-based
  end-exclusive `start`/`end` as every endpoint takes them, plus 1-based
  `start_1based`/`end_1based`; its own sequence 5'→3' on its strand, spliced
  parts joined; for a CDS the translation of the bases as they stand, with
  `translation_matches_annotation` when the file carries a /translation, and
  `protein` — the protein as the CDS's own annotation reads it. The two differ
  only at residue 1 of a CDS whose GTG / TTG start the file declares as Met
  (its `/translation` begins M and agrees with the bases over its N-terminus,
  or a `/transl_except` puts Met there): then `translation` keeps the codon
  table's letter, `protein` has M, and `initiator` says so —
  `{codon, codon_reads, declared_by: "translation" | "transl_except"}`; and
  `unique_cutters_inside` — every enzyme that cuts the whole plasmid once with
  that cut inside the feature, as `[{enzyme, cut_bp}]`);
  `amplify-feature` `{name | idx, type?, target_tm?}` (binding-only PCR primers
  that amplify one feature end to end — across the origin and on either
  strand — then locate each primer on the whole plasmid: `binds_as_designed`,
  `other_sites` ≥80% that could misprime, and `primer_at_feature_start`, which
  for a − strand feature is the REVERSE primer). `plasmid-overview`'s
  `unique_cutters` is grouped by the feature each cut falls in —
  `{label: [{enzyme, cut_bp}]}`, plus "(outside every feature)", with
  `unique_cutter_count` alongside — so "which enzymes cut once inside X" is
  a lookup. An ambiguous name is a 409
  listing candidates to pick by `idx`. All read-only; they exist so a client
  gets exact numbers without chaining endpoints and copying coordinates
  between them.
- **Files** — load (chromosome-scale safe via the path-based loader;
  supports `.gb` / `.gbk` / `.genbank` / `.dna` / `.embl` /
  FASTA / `.ab1` / single-record `.fastq` / `.gff3`),
  export GenBank / GFF3 / FASTA / EMBL / CommercialSaaS `.dna`
  (symlink-guarded), bulk import a folder (`bulk-import-folder` answers
  `422` / `ok: false` when every file failed, and then creates no collection,
  so a retry is not met with a 409; when only some failed it is `ok: true`
  with `partial: true`, the `failures` and a `warnings` line), bulk export a
  collection — to any of those formats or a circular-map image
  (`png` / `svg`) — via `bulk-export-collection`, or export the
  LOADED plasmid's circular map alone (`png` / `svg`) via
  `export-map-image`.
- **Library + collections** — list, search across collections,
  load an entry by name or id (`load-entry` resolves cross-collection:
  active collection first, then the others; pass `collection` to
  disambiguate), rename a plasmid's display name (`rename-plasmid` —
  fixes a name that got slugged to underscores), delete entries,
  copy-plasmid (copy an entry into another collection, name + features
  intact, in one call), move-plasmid (RELOCATE an entry between
  collections atomically — no copy-then-delete data-loss round-trip;
  handles the active collection on either side, re-staging the live
  library mirror), copy-plasmids / move-plasmids (the BULK form — a
  whole `{names:[…]}` list into one collection in a single locked save,
  the machine-friendly path that sidesteps the rate limiter; per-item
  results report copied/moved, not_found, ambiguous, and conflict, never
  silently renaming or dropping). A move is committed to the collections
  file first; if refreshing the live library mirror then fails, the move
  answers `ok: true` with a `warnings` entry saying so (the same holds for
  `move-part`) instead of reporting a failure for a move that happened.
  Create / rename / delete collections (deleting the ACTIVE one switches to
  the next and returns it as `active`), get / set the active collection,
  list / set plasmid statuses.
- **Parts** — list-parts, get-part, create-part / update-part (full
  parity with the Part editor — name / grammar / type / level / position /
  overhangs / sequence plus `backbone`, selection `marker`, and
  `fwd_primer` / `rev_primer` domestication primers, Tms derived; a part
  whose clone would re-create a Type IIS site at a junction comes back with
  `junction_sites` and a warning, since it would be cut in its own assembly),
  delete-part, move-part (reassign a part to another bin in one atomic
  call), classify-part (overhang-pair lookup against every grammar),
  make-l0-part-from-fragment (turn a saved synthetic FRAG into a
  Level-0 part *and the L0 plasmid it lives on*, by SIMULATING the
  clone into the grammar's entry vector — body `{fragment, part_type,
  [grammar], [name], [bin], [collection], [product_name],
  [save_plasmid]}`; the part goes to the bin and the circular plasmid
  is saved to the library as `<part> (L0)` with its lineage, via
  `domesticate-part` so there is one definition of an L0 clone
  (`save_plasmid: false` files only the part); the part's
  overhangs come from that grammar's position table, so a two-tier
  nested fragment files its CATEGORY pair rather than the entry pair a
  digest would report; 400 names the cause for an unwrapped fragment, a
  part type the grammar doesn't define — the error lists the ones it
  does — a part type whose overhangs aren't on the fragment, or a
  grammar with no entry vector assigned; 409 on a duplicate
  `(name, grammar)`. The filed `sequence` is the body BETWEEN the
  overhangs, exactly as `create-part` expects, so the part re-clones
  into the plasmid the call simulated; the reply carries `part_len`,
  `insert_len`, `plasmid_len`, `entry_vector`, `nested`, and
  `saved_name` / `saved_id` / `collection` for the plasmid — plus
  `plasmid_error` in the rare case the part filed but the plasmid
  didn't, which is reported rather than swallowed).
  `create-part` / `list-parts` / `delete-part` accept a `{bin}` to file
  into / read / prune just that named bin's own parts (a bin is a real
  partition, not only for writes). Parts-bin containers:
  list-parts-bins, create-parts-bin (the parts-side parallel to
  create-collection), set-active-parts-bin (switch the active named bin;
  mirrors the bin into the live parts file), rename-parts-bin,
  delete-parts-bin (the deleted bin's parts go with it; the last bin
  can't be removed; the active bin auto-promotes + re-mirrors).
- **Design** — gibson-assemble, simulate-gibson, traditional-clone /
  simulate-traditional-cloning (restriction digest + ligation: excise the
  insert, digest the vector, try every vector-fragment × insert-fragment ×
  orientation, and save the product. Each product carries a `self_ligation`
  report — whether the EMPTY vector re-closes (a single-enzyme cut, a blunt
  cut, or a compatible-cohesive pair such as SalI/XhoI, which both leave
  TCGA), how big it would be, and whether a diagnostic digest can even see
  it, plus double-insert / partial-assembly / orientation risks each with a
  severity and the fix; `clean: true` carries a note naming which ends
  disagree, so an empty risk list never reads as "not checked" — refuses to guess when more than one
  ligation is possible, pick with `vector_frag_idx` / `insert_frag_idx`. Pass
  `insert_circular:true` to cut a cassette OUT of a *plasmid* insert — e.g. an
  Ω multigene into a binary vector — so its two digest fragments are both
  offered instead of treating the insert as a linear PCR product. Pass
  `carry_annotations:true` to lift the `vector_name` / `insert_name` entries'
  own features onto the ligated product — split at the cuts and shifted to
  product coordinates — so the clone isn't feature-bare; a side whose passed
  sequence doesn't exactly match its named entry is reported in
  `carry_warnings` rather than mis-placed. Both names resolve **across every
  collection** (pin either with `vector_collection` / `insert_collection`;
  `409` on a genuine ambiguity), so a backbone and an insert filed apart
  need no staging copy. **Check the per-side `carried`
  `{vector, insert, any}`** — a skipped side otherwise reads as `ok: true`
  with half the annotations missing; `carry_strict:true` makes that skip a
  `422` instead. A restriction site the ligation REGENERATES (the normal
  outcome at a directional clone's junctions) comes back as one intact
  feature, not two `(disrupted)` halves, and each junction's sticky overhang
  is annotated at the bases that actually anneal),
  golden-gate-assemble / simulate-golden-gate (Type IIS — BsaI /
  BsmBI / BbsI / SapI / Esp3I — overhang-directed N-part assembly: parts in
  any order, chained by their 4-nt overhangs into a circle, with a
  unique-overhang + no-residual-site fidelity check. Pass
  `carry_annotations:true` with each part / the vector given as
  `{sequence, name, collection?}` to get an ANNOTATED product: each named
  entry's own features are slotted onto its released fragment by the digest
  and shifted into product coordinates by the chaining ligation, so the
  feature map falls out of the same simulation that built the sequence.
  `simulate-golden-gate` accepts the same payload, so a preview shows the
  same map the save will store. Per-input `carried`
  `{vector, parts: [...], any}` plus `carry_warnings` say what landed;
  `carry_strict:true` turns a sequence mismatch into a `422`. Omit the flag
  and the product saves feature-bare as before. **The vector may be cut
  more than twice** — a destination plasmid with background enzyme sites is
  released as several pieces and all of them go back into the circle, which
  `n_vector_fragments` / `n_vector_fragments_used` / `vector_cut_bp` report.
  When nothing closes, the error names the vector's cut positions and every
  fragment's two overhangs (`CGCT→TTAC, …`) so you can see which end has no
  partner, instead of being pointed at parts that were fine; more than one
  possible circle is warned about rather than silently picked),
  lint-synthesis (an "is it safe to order + assemble?" pre-flight over a bare
  `sequence` or a library `id`/`name`: internal Type IIS sites, extreme
  overall + windowed GC, long homopolymer runs, tandem repeats, degenerate
  bases, and — with `expect_cds` — a full-length ORF check, rolled into a
  0-100 `score` + line-item `warnings`; read-only), design-mutagenesis,
  scrub-plasmid (clone-free restriction-site removal: silent / synonymous
  cures inside CDSes + minimal swaps elsewhere; scrubs the loaded record or
  an explicit `seq`+`features`, optional `codon_taxid` biases coding cures to
  a host's frequent codons, never mutates the canvas. `method` picks the
  route: `"quikchange"` (default — cured sequence + an improved-QuikChange
  primer pair per locus) or `"golden_braid"` (split into BsaI-tailed
  fragments that Golden-Gate back together; force-cures every BsaI site,
  returns per-fragment primers + native junction overhangs + a digest+ligate
  `verified` flag). Cures prefer the silent option that does NOT lengthen a
  homopolymer run — ranked above the codon-usage tie-break, because a run
  long enough to stall synthesis is a hard failure and a rarer synonymous
  codon is not; `max_homopolymer_run` comes back on every plan, with a
  warning when no alternative avoided lengthening one. Each QuikChange round
  reports `hairpin_dg` /
  `homodimer_dg` in **kcal/mol** (`dg_units` says so), taking the worse of
  the two primers so the number works as a gate; `null` means it could not be
  measured — never "no structure" — and the round carries a warning saying
  the pair is not secondary-structure checked), design-gb-part (Golden Braid / MoClo domestication
  primers; pass `check_entry_vector:true` — plus optional `entry_vector_role`
  for a Golden-Braid role — to validate at DESIGN time that the part's
  overhangs actually match the configured entry vector's acceptor, so a
  mismatch surfaces here instead of at the clone step; verdict rides under
  `result.entry_vector_check`. **Two-tier aware:** if the configured L0 entry
  vector exposes EXTERNAL acceptor overhangs that differ from the position's
  category overhangs — a UPD/pUPD vector, e.g. external `CTCG`/`TGAG` vs category
  `AATG`/`GCTT` — the category pair is auto-NESTED inside the external pair so
  the fragment enters on Esp3I (external) and later releases on BsaI (category);
  actual cut overhangs return as `result.entry_oh5`/`entry_oh3`),
  design-synthesis-fragment (the direct-synthesis counterpart: wrap a `sequence`
  in the SAME nested overhangs so the SYNTHESISED fragment — ordered as a gBlock,
  no PCR — clones into the UPD vector as an L0 part. Pass `{sequence, part_type?,
  oh5?, oh3?, grammar?}` — a grammar POSITION type OR a custom 4-nt overhang
  pair; auto-nests from the configured entry vector like design-gb-part. Returns
  `result.fragment` (+ `length`, category `oh5`/`oh3`, external
  `entry_oh5`/`entry_oh3`, `nested`, `entry_vector`, and `internal_sites`
  flagging any grammar-forbidden site in the insert that would self-cut)),
  domesticate-part
  (the real Synthesis-tab L0 clone: digest the grammar's CONFIGURED entry
  vector + the part's primed amplicon at the Type IIS enzyme, ligate the
  insert into the backbone, and save the circular L0 plasmid with its
  insert + entry-vector lineage — pass `{sequence, oh5, oh3, name,
  grammar?, type?, fwd_primer?, rev_primer?, collection?}` (saves to the
  active collection, or the non-active `collection` if given). **Fails
  loud** (422) when no
  compatible entry vector is configured rather than silently emitting a
  pUPD2-stub backbone, so an agent never files a wrong construct;
  `domesticate-parts` runs the same engine over a `{parts:[…]}` batch in
  ONE call — dirty-guard checked once, per-item results so one part that
  can't clone doesn't abort the rest. Each item gets the routing-key check a
  single call gets (an item naming a key the engine ignores is that item's
  `code: 400` + `unsupported`); a batch that built nothing answers `422`,
  one that built some `ok: true` with `partial: true` + `warnings`),
  assemble-into-entry-vector (the multi-source level-up clone: chain L0
  parts into an α/L1 TU, or TUs into an Ω/L2 module, by their fusion
  overhangs and ligate into the configured `role` acceptor at the level-up
  Type IIS enzyme — pass `{sources, grammar?, source_level?, role, name,
  collection?}` (active collection by default, or the named one).
  Also **fails loud** (422) on a missing acceptor or a chain that can't
  close; the L2→binary hop is a known engine gap), design-primers
  (Primer3 primer-pair design over a region: `detection` picks a
  diagnostic amplicon, `cloning` appends RE-site tails, `generic`
  returns binding-only primers — no tails / overhangs; and
  `allele_specific`, which is anchored on a single `variant_pos` rather than a
  region and pins the primer's **3'-terminal base** on the variant, because
  that pinning IS allele-specific PCR — the other modes slide to whatever
  window scores best, and the Golden-Gate SDM designers anchor the 5' end, so
  neither can express it. Sweeps 18-27 bp and both orientations around the
  fixed terminus; returns `specificity`, the 3'-seeded binding sites on the
  template and on any named `off_targets`, because "this primer is
  allele-specific" is a claim about where it does NOT bind. When the
  discriminating base is A or T it SAYS there can be no GC clamp instead of
  only scoring it down, and `destabilize` adds the conventional second
  mismatch three bases from the 3' end — opt-in, since it changes the oligo
  you order), check-primer
  (a SINGLE oligo vs a template: Tm, GC%, and every 3'-anchored binding
  site — both strands, wrap-aware on a circular template, the tail
  scored as mismatches — to confirm a designed primer binds exactly once
  before saving; `template` defaults to the LOADED plasmid with its own
  topology, and with no template and nothing loaded it still returns the
  Tm and GC% with `binds: null` and a `note`; `sequence` is the template
  when `primer` is given and the primer when it isn't — flagged
  `primer_from: "sequence"`, refused past 1000 characters. Both are read
  strictly: a FASTA header line, spacing, numbering and 5'/3' notation are
  fine, but text mixed into the bases — a label (`kanR-R: ACGT…`), a second
  FASTA record, a letter that is not a base — is a 400 rather than being
  read as bases. `U` IS a base here: a deoxyuridine primer is read as T for
  binding and Tm, and the reply carries `read_u_as_t` so the reinterpretation
  is never silent. A degenerate oligo gets a real nearest-neighbour `tm` —
  its WEAKEST variant, which is what sets the anneal — plus `tm_range`
  spanning the mix. Sites come back ranked, and `truncated` says whether
  `max_sites` cut the list short), optimize-protein
  (codon-optimise an AA sequence to a chosen table — `table` is a taxid
  (integer or string) OR a table NAME, and one it cannot resolve is a `404`,
  never a quiet fall-back to the E. coli default; optional `stops`
  0–3 appends that many stop codons, and a trailing `*` run in the
  protein is honored as-is and overrides it; optional `transl_table`
  is an NCBI genetic-code id — default 1 = standard — so a host with a
  reassigned codon, e.g. 4 *Mycoplasma* or 6 ciliate, optimises and
  terminates against the code it actually reads; `mode` is `frequency`
  (default, matches the host's codon distribution) or `max_cai`
  (every position gets its best synonym); optional scrub passes —
  `hazard_hosts` (built-in plant/mammalian/bacterial expression
  hazards), `forbidden_motifs` (your own IUPAC patterns), `min_gc` /
  `max_gc` over a `gc_window`, and `avoid_repeats_with` (break shared
  runs against constructs you've already built) — each opt-in, with the
  response reporting `fixes`, `remaining_motifs`, `remaining_repeats`,
  and the achieved `gc_window_min` / `gc_window_max` so an unreachable
  target is visible, not implied met. The response also carries `gc3`
  (third-position "wobble" GC) and a `warnings` list: `max_cai` on an
  AT- or GC-biased host stacks that host's single top codon everywhere,
  collapsing GC3 — a silencing / mRNA-instability risk — while overall
  `gc` still looks healthy, so it warns; on such a host prefer `frequency`
  mode (matches the native GC3) or set `min_gc`). Only E. coli K12 ships
  as a table — `add-codon-table` fetches your host from Kazusa by taxid,
  builds one from an NCBI genome, or imports a local CDS FASTA.
- **Simulate** — simulate-pcr (exact-match in-silico amplification,
  wrap-aware on circular templates) and simulate-gel (per-lane band
  positions + optional rendered ASCII gel image; ladder / plasmid /
  digest / PCR-amplicon sources).
- **Restriction sites** — `list-restriction-sites` scans the LOADED record
  (singular `enzyme` accepted; omit it to scan the catalog, honouring your
  active enzyme collection). Each hit carries `start` / `end` / `strand`,
  both strand cuts (`cut_bp` / `bottom_cut_bp`), the recognition `site`, and
  `wraps`. Two rules keep a zero-site answer honest: **an unrecognised
  enzyme name is a `400`**, never `count: 0` (a typo used to read as "no
  sites here" instead of "I never looked"), and **isoschizomers resolve by
  the name you asked for** — the scan collapses enzymes that share a
  recognition site and cut (BsmBI / Esp3I / BsmBI-v2, BsaI / BspTNI, EcoRI /
  EcoRI-HF) onto one label so the map stays readable, so the endpoint
  re-expands it and reports `equivalent_enzymes` telling you which other
  spellings the answer already covers. **Commercial synonyms resolve too** —
  the same enzyme is sold under different names (Thermo's `Eco31I` is NEB's
  `BsaI`, `LguI` is `SapI`, `AarI` is `PaqCI`, `Asp718I` is `Acc65I`), and
  hits come back labelled with the name YOU asked for. **NEOschizomers stay
  apart**: XmaI (`C^CCGGG`) and SmaI (`CCC^GGG`) read the same six bases and
  leave different ends, so this endpoint reports both even though the plasmid
  map draws one bar for the pair — same for Acc65I/KpnI, ApaI/PspOMI,
  NheI/BmtI, KasI/NarI, AatII/ZraI, SacI/Eco53kI. Two further refusals exist for the same
  reason — both would otherwise return a confident `count: 0`: a `min_length`
  longer than the longest recognition site in the catalog (a units mix-up —
  `min_length` is in BASE PAIRS), and an `enzymes` list longer than 500 names
  (which also bounds the near-miss search behind the error message).
  Each hit also carries **`cuts`**: `false` marks a site sitting so close to
  an end of a LINEAR molecule that the enzyme's cut falls outside it. The
  recognition sequence is really there — and cuts as soon as the fragment is
  in a vector — it simply severs nothing in this molecule; such sites used to
  be dropped entirely, so a part with a BsaI site three bases from its end
  answered `count: 0` and a scrub called it clean. A `warnings` entry names
  them. When you rely on your **active enzyme collection** (the default),
  `warnings` also names any enzyme in it the catalog does not know; a
  collection whose names ALL fail now scans nothing and says so, rather than
  silently widening to the whole catalog.
- **Digest** — digest (cut a RAW sequence with named enzymes and report
  the cuts + resulting fragments with their **overhangs** — overhang-aware
  QC for a Golden-Braid / restriction junction without loading the
  sequence onto the canvas; with no `sequence` it digests the LOADED plasmid
  with its own topology and says so in `sequence_from`; `circular` defaults
  true for a raw sequence, a singular `enzyme`
  is accepted, names are matched case-insensitively so `bsai` cuts,
  commercial synonyms resolve (Thermo's `Eco31I` is NEB's `BsaI`) with
  `resolved_enzymes` reporting what each one became, and names the catalog
  genuinely doesn't know are reported under `unknown_enzymes` rather than
  silently dropped — this endpoint reports rather than refuses, so a long
  list still returns the cuts it can make. Only a request in which NO name
  resolves is refused, `400` with `unknown_enzymes`: answering it would
  present the uncut molecule as the digest).
- **Transcripts** — `predict-transcript` reconstructs the MATURE mRNA of an
  annotated transcription unit on the loaded record and scores translation
  initiation on it. Pick the unit with `promoter` / `cds` / `terminator` (a
  feature index or a label) or let it take the longest CDS and the nearest
  promoter/terminator on the same strand; `tx_start` + `tx_end` bound it
  directly when nothing is annotated. Returns `pre_mrna` / `mature_mrna`, the
  `exons` and `introns` between them (with donor/acceptor and a `canonical`
  GT–AG flag), `five_utr` / `cds` / `three_utr`, the `kozak` context, and
  `uorfs` classified `upstream` / `in_frame_extension` /
  `out_of_frame_overlap`. **The field to read is `removed_by_splicing`:**
  every upstream ATG that is in the DNA and is NOT in the message, with the
  intron that deletes it — the false positives a uATG screen run against the
  plasmid sequence reports as hazards. Introns come from ANNOTATION only
  (`intron` features, or the gaps between a spliced location's parts); nothing
  here predicts where a spliceosome would cut. `splice` carries a separate
  ADVISORY scan of the pre-mRNA with the calibrated plant PWM, marking sites
  that coincide with an annotated boundary so only genuinely cryptic ones are
  listed — and reporting `skipped` with a reason if the model is unavailable,
  never an empty list that would read as clean. `include_sequences: false`
  drops `pre_mrna` / `mature_mrna` / the UTR sequences and keeps the
  coordinates and verdicts, for a screening loop over a whole collection.
- **Regulatory elements** — `scan-promoters` finds sigma-70 (-35 / spacer /
  -10) promoters and `scan-terminators` finds intrinsic (Rho-independent)
  terminators, on `{sequence | id | name}`, wrap-aware, with FORWARD-strand
  coordinates for hits on either strand. A promoter hit carries both hexamers,
  the spacer, the predicted `tss` (`null` on a linear molecule when the start
  site would fall past its end — it never wraps to bp 0), a `score` and its
  `components`; **the score
  RANKS, it does not predict** — it is a composite of named parts with no
  trained model behind it, good for "does this insert carry a cryptic promoter
  and how does it compare with the others", not for a transcription rate. The
  default `min_score` is calibrated rather than picked: 0.75 yields about one
  chance hit per kb per strand, while real promoters from the built-in preset
  catalogue score 0.69-1.00. A terminator hit carries the stem/loop/U-tract
  geometry, a real Turner-2004 fold `dg`, and a coarse `efficiency` tier
  (`strong` / `moderate` / `weak`) — never a percentage. **One hairpin is ONE
  record:** overlapping stem/loop registers of the same structure are merged,
  because a scan that reports each register separately turns one terminator
  into twenty-odd. Rho-DEPENDENT termination is not modelled and is not
  guessed at. `min_stem` and `min_u_tract` are bounded by the scanner's own
  search limits: a floor above them could only ever return zero hits, and
  "0 terminators" reads as "this construct has none" rather than "that filter
  excluded every candidate", so it is refused with a 400.
- **Transcription mapping** — `map-transcription` walks the circle from every
  promoter to every CDS and answers **what ELSE transcribes this gene**. Body
  `{id? | name?, include_predicted?=true, min_promoter_score?, max_distance?}`;
  needs an annotated record, so a bare sequence is refused rather than silently
  answered with no genes. Per gene: `reachable_from` — each promoter that
  reaches it, the `distance_bp` around the circle, the `terminators_between`,
  and a `readthrough` of `clear` / `attenuated` / `blocked`.
  `predict-transcript` cannot answer this: it reconstructs ONE unit and needs
  the promoter named. A plasmid is a circle, every promoter on it transcribes
  until something stops it, and checking terminators DOWNSTREAM of a cassette —
  the natural thing to do — cannot see a read-through, because the read-through
  arrives from behind. **`silent` and `silent_in_host` are different claims and
  the gap between them is the point:** `silent` means nothing is aimed at the
  gene at all, `silent_in_host` means nothing the HOST polymerase can read
  reaches it. A T7-only cassette is `silent_in_host` in a strain with no T7
  RNAP — unless something else on the circle reads through into it, which is
  exactly the case that looks like a mystery at the bench. Each promoter
  carries `host_recognised` plus the `host_basis` that decided it. A promoter
  or terminator annotated with no strand is read FORWARD, flagged
  `strand_unknown: true`, and named in `warnings` — it is never silently
  dropped, which would make an arrowless terminator block nothing.
- **Read heterogeneity** — `analyse-read-heterogeneity` reports per-base
  ALLELE FRACTIONS from raw reads: **is what I sequenced one thing?** Body
  `{reference | reference_id | reference_name, reads_path (a FASTQ) | reads,
  platform?="unknown", min_fraction?=0.01, min_phred?=20, noise_floor?,
  homopolymer_filter?, max_reads?}`.

  **Pass `platform`** — `ont` / `illumina` / `sanger`. It sets the noise floor
  (the fraction below which a call gets no vote in the verdict, distinct from
  `min_fraction`, which controls what is *reported*) and whether homopolymer
  indels are muted. Judged as `unknown`, a nanopore run is measured against a
  floor built for short reads; the response warns when the data carries the
  long-read signature (most calls being indels inside homopolymers) and names
  the share that led it to say so. `read-consensus` cannot
  answer this — it reads the alignments STORED on a plasmid, and a
  Plasmidsaurus consensus is one read, so depth is 1 everywhere and anything
  below the consensus is invisible. A population of escapers each carrying a
  different inactivating mutation, none dominant, consenses back to wild type
  and reads as a clean clonal culture. Returns `positions` (each `{pos, ref,
  alt, alt_reads, depth, fraction, mean_phred, contradicted}`) and a `verdict`
  of `clonal` / `minor_variants` / `mixed`. **The denominator is every read
  that COVERED the base**, including the ones that disagreed; observations
  whose backing basecall is below `min_phred` are excluded and counted in
  `n_filtered_low_quality`, because at 1% every instrument's error floor looks
  like a sub-population. **A no-call is not an allele:** an `N` in a read means
  the basecaller could not read that base, so those observations are excluded
  and counted in `n_no_call_observations` rather than reported as a
  sub-population — long reads carry N runs routinely, and counting them made a
  single ordinary read turn a clonal culture into `mixed`.
  The thresholds behind the verdict are echoed in
  `thresholds` — they are conventions, not a calibrated model. Reads are
  bounded by a `reads x reference length` work budget and by `max_reads`;
  whenever either trims the set — the DEFAULT `max_reads` included —
  `capped_by` names it and a warning says so, rather than quietly changing
  every fraction (`capped_by` is null when every read was used). From a FASTQ
  only the first `max_reads` reads are aligned, but `reads_available` counts
  every read in the file. A mode is judged on the calls AT OR ABOVE the noise
  floor, never on a bin that merely straddles it, and a reported
  `mode_fraction` is never below the floor.

  **The verdict keys off the SHAPE of the distribution, not a count.** A count
  rule ("5+ positions above 1%") calls every long-read dataset `mixed`:
  measured against a simulated clonal nanopore run it gave 1,347 positions,
  78% inside homopolymers, top fraction 0.30 — clearing both the count rule
  and the dominant-fraction rule on a single clone. A verdict that always says
  `mixed` is worth what one that always says `clean` is worth. So `mixed` now
  needs either one substantial sub-population (a position at/above
  `mixed_fraction`) or a MODE — a bin standing above the one below it, which
  is what a linked haplotype produces, because every read of an escaper
  carries all of that clone's differences at one shared fraction. Independent
  per-position error decays instead. `distribution` returns the bins,
  `monotonic_decay` and `mode_fraction` so the verdict can be checked rather
  than trusted, and each position carries `counted_in_verdict` plus `muted_by`
  (`homopolymer` / `noise_floor` / null) so "why did this not count" is
  answerable from the response alone.

  **Homopolymer indels are reported but do not vote** (on by default; every
  platform miscalls run length, long reads most of all). They stay in
  `positions` tagged `context: "homopolymer"` with `homopolymer_len`, are
  counted in `n_muted_homopolymer`, and a warning fires whenever muting
  removed anything — because a real frameshift escaper can live in a polyA
  tract, and that is the one case wanting `homopolymer_filter: false`. A
  SUBSTITUTION inside a homopolymer is not a run-length miscall and keeps its
  vote. Whenever the noise floor drops a call and the answer is not
  `mixed`, a warning names the highest fraction it excluded — the floor fails
  in the dangerous direction (a false `mixed` costs a re-screen; a false
  `clonal` lets an escaper through), so it is never silent. Boolean fields are
  validated rather than coerced: `bool("false")` is True, and a flag whose
  whole purpose is turning a filter OFF must not invert on a string.
- **CDS codon analysis** — `analyse-cds` is the READ half of
  `optimize-protein`, which until now only wrote. Body
  `{sequence | id + feature, taxid?, rare_w?, ramp_codons?}`; returns
  `n_codons`, `gc`, `gc3`, `cai`, `rare_pct`, `tandem_rare_runs`,
  `worst_window` and `ramp`. `cai` comes from the same function
  `optimize-protein` reports, against the same table, so the two are directly
  comparable. `tandem_rare_runs` matters separately from `rare_pct` because
  consecutive rare codons stall a ribosome far more than the same codons
  scattered, and `ramp` reports the 5' end on its own because a deliberately
  slow translational ramp there is normal — averaging it in flags healthy
  designs and hides unhealthy bodies. **A spliced CDS is joined from its
  exons** and a `/codon_start` of 2 or 3 is honoured, both reported in
  `warnings`: reading a spliced gene as one contiguous block runs through its
  introns, so every codon after the first is out of frame and every number is
  then confidently wrong. An origin-WRAPPING CDS is a two-part location too
  (sacred invariant #9) and is deliberately NOT read as spliced.
- **Enzymes** — `list-enzymes` reads the COMBINED catalog (built-in NEB ∪
  your custom enzymes; `list-custom-enzymes` returns only the latter) with
  no sequence involved: recognition site, `fwd_cut` / `rev_cut` offsets,
  REBASE-style `notation` (`G^TCGAC`, `GGTCTC(1/5)`), overhang length +
  kind, a `type_iis` flag, and `aliases`. Filter with `names` (exact,
  case-insensitive), `search` (substring over the name AND the recognition
  site, so `GGTCTC` finds BsaI and BspTNI), `type_iis_only` or
  `custom_only`. Exists so a script can check a digest's fragment
  arithmetic against the enzyme's real offsets instead of hardcoding a
  private table of them.
- **Coordinates** — `span-contains` answers "is this base / feature inside
  that region?" on a molecule where **either** span may cross the origin.
  Body `{outer, inner, length?}`; spans are half-open `[start, end)` — the
  convention `list-features` and `list-restriction-sites` report — and one
  with `end <= start` wraps. `inner` takes one span or a list, and the
  response shape is the same either way (`results[]` + `all_contained`), so
  callers never branch on it. **Use it instead of `a <= x < b`:** the linear
  form reports every base of an origin-spanning T-DNA, marker or operon as
  OUTSIDE it, failing constructs that are actually correct. `list-features`
  reports the same thing per feature via `wraps` and a wrap-aware `length`.
- **Alignment** — diff-plasmid (one target, circular rotation
  auto-detected), multi-align (batch: the loaded plasmid or a given
  sequence vs many targets at once — the Alt+A overlay; rotation-aware
  per target, circular auto-detected from each target's topology or
  forced with a `circular` flag, with `picked_rotation` reported per
  row), list-plasmidsaurus-members, align-plasmidsaurus-zip,
  verify-against-reads (the "I built X, I got reads back — do they match?"
  check: a bare or library-resolved `reference` vs a list of raw `reads`
  (Nanopore / Sanger / a consensus), each aligned rotation + RC-aware; returns
  per-read identity% + a `match`/`inverted`/`mismatch` `verdict` against
  `min_identity`. `inverted` is its own verdict because a construct built
  BACKWARDS is a different next step from one that simply failed: each read
  reports the spans matching the reference's reverse complement in
  `inversions` + `inverted_bp`. `multi-align` rows and `diff-plasmid` carry
  the same `inverted_segments`).
- **History** — get-history returns the parsed `<HistoryTree>`
  lineage as nested JSON — the FULL per-node record, matching what the
  History tab shows: name / operation / length / topology / date /
  node_id / resurrectable, the input summaries and regenerated sites,
  the PCR `oligos`, the `primers` block (each with its binding sites:
  location, strand, annealed bases, Tm, and the unannealed 5′ flap),
  `hybridization_params`, run `parameters`, `end_modifications`,
  `sticky_ends`, `custom_map_label`, `feature_snapshot_count`, and
  `extra_attributes` / `other_children` for anything the parser doesn't
  model yet (so a field from a newer writer still reaches you). A
  `warnings` array lists claims the record makes that the plasmid's own
  sequence does not support — an absent regenerated enzyme site, a primer
  that moved or no longer binds, a length-preserving step that changed
  the length, deletion arithmetic that doesn't balance. Checked against
  the ROOT node only (an ancestor describes an earlier molecule, so
  checking it against today's sequence reports errors that aren't real).
  Read-only: nothing is ever rewritten to clear a warning.
  **recover-history-from-dna** restores lineage to plasmids whose saved
  `.dna` original documents more of their build than the library holds:
  it matches on EXACT sequence identity (so it works across a rename),
  adopts the source history only when it has strictly more nodes, and
  writes `history_xml` and nothing else. `dry_run` defaults to `true` —
  the first call reports what would change; pass
  `{"dry_run": false}` to apply. Optional `collection` / `name` narrow
  the scan itself: only `.dna` files whose sequence matches a plasmid in
  scope are indexed, so a run that stopped at the memory budget
  (`truncated: true`) is finished a collection at a time.
  Agent assemblies (`gibson-assemble`,
  `traditional-clone`, `golden-gate-assemble`) attach real parent
  lineage — each input fragment / part / vector becomes a parent node, so
  the product reads as a genuine `insertFragment` assembly instead of a
  flat `createDocument` leaf. Name your inputs (a `{name}` on a gibson
  fragment / golden-gate part, or `vector_name` / `insert_name` on a
  traditional clone) and a parent that matches a saved library entry
  nests that entry's own sub-lineage.
- **Codon tables** — list, add (Kazusa fetch or raw dict), delete.
- **Search** — blast (in-process BLASTN / BLASTP against your
  collections; pass `collections` to scope the DB to a subset and stay
  under the build cap on a large library), hmmscan. **Online** (off by
  default): blast-online (remote NCBI BLAST vs `nt`/`nr`/…) and
  hmmer-web (remote hmmscan vs Pfam at EBI) SHIP the query sequence to
  an external server, so they are refused `403` unless the SpliceCraft
  user has armed Settings → "Allow agent online BLAST/HMMER". That
  setting is the HUMAN half of a two-key gate — it is deliberately NOT
  in the agent settings allowlist, so an autonomous agent can never flip
  it on itself; the in-process blast / hmmscan stay local and never
  consult it. Both block for minutes on the remote job (502 on a
  network / server error) and are concurrency-capped.
- **Online reference-database lookups** (off by default) — fpbase-search
  (FPbase fluorescent proteins: spectra, oligomerization, xrefs, sequence),
  uniprot-search (UniProtKB proteins: function, organism, keywords),
  literature-search (Europe PMC papers: title/authors/journal/DOI/abstract),
  genbank-search (NCBI `nucleotide`/`protein` term search → accessions you can
  `fetch`; returns `total_matches`), wikipedia-search (lead-section summaries),
  web-search (general web), patent-search, and read-url (open ONE page by URL
  and return its readable text — HTML stripped to text, no JavaScript run; the
  companion to web-search, so Babs can open a result and actually read it, and
  it refuses binary / PDF pages with a `502`). All read-only. The `*-search`
  endpoints send only a QUERY STRING to the public database — read-url sends
  only the URL you give it — never your sequence (distinct from the BLAST/HMMER
  egress above), but egress still requires the SpliceCraft user to
  arm Settings → "Allow Babs online database lookups" (`403` otherwise). That
  setting is the HUMAN half of the gate: deliberately NOT in the agent settings
  allowlist, so an agent can't self-arm. web-search uses the Brave Search API
  when a `brave_search_api_key` is set (else keyless DuckDuckGo, best-effort);
  patent-search uses PatentsView when a `patentsview_api_key` is set (else
  keyless Google Patents, best-effort). Both keys are secret — redacted in logs
  and unreadable via the agent settings surface. A keyless provider that is
  rate-limited returns `502` with a hint to add a key. Upstream/network failures
  → `502`. Every fetch routes through the SSRF-hardened opener with the
  fail-closed DEMO_MODE egress gate.

- **Background Learning** (`learn-start`, `learn-status`, `learn-results`,
  `learn-list`) — grow Babs' corpus about a TOPIC with a focused, drift-resistant
  crawl of open databases + open-licensed web pages. `learn-start {topic, max_docs?,
  max_minutes?}` is the only **write** (it launches a crawl): refused `403` unless
  the user has armed the same "Allow Babs online database lookups" setting, and
  write-flagged so the in-app Babs prompts before an autonomous agent starts one.
  It returns immediately with a `slug`; poll `learn-status {slug}` for live counters
  (kept / fetched / dropped / frontier / status / drift), read `learn-results {slug,
  limit?}` for the scored kept papers, and `learn-list` for all sessions on the box.
  Only topic query strings leave the machine — never your sequences. The crawl
  engine lives in the separate Babs repo (`bb_learn.py`); each topic is isolated in
  its own corpus + vector collection.
- **Corpus, memory + reference** (`recall-knowledge`, `index-library`,
  `ingest-url`, `remember-fact`, `memory-list`, `memory-forget`,
  `reference-lookup`, `learn-merge`) — what Babs KNOWS, as opposed to what she can
  look up. `recall-knowledge {query, max_passages?}` retrieves scored passages
  (with sources) from the local corpus. `reference-lookup {query, registry?,
  max_hits?}` (or `{list_registries: true}`) reads Babs' CURATED registries —
  antibiotics/selection agents, promoters, terminators, UTRs, origins, Type IIS
  enzymes, assembly standards, markers, *E. coli*/*Agrobacterium* strains, binary
  and CRISPR vectors, base/prime editors, localization signals, chassis, media,
  growth regulators, bench recipes, troubleshooting — returning WHOLE structured
  records (working concentrations, cut notation, host ranges, genotypes) rather
  than retrieved prose. Prefer it over corpus recall for a named reagent, enzyme,
  strain or part: deterministic, instant, no embedding step. Both read-only and
  fully LOCAL — they send nothing anywhere.

  The **writes** all grow or prune what Babs believes, so each is write-flagged
  and the in-app Babs asks first. `index-library {}` indexes the user's OWN work
  (plasmids, primers, parts, lab notebook) into a **private** pack so recall can
  answer "which of my plasmids carries a KanR and a T7 promoter?" — it READS the
  SpliceCraft data dir and writes only into the Babs repo, rebuilding the pack
  each run so a renamed or deleted plasmid can't linger as a ghost. That pack is
  private by construction: every record flagged private, a `.no-egress` marker in
  its directory, excluded from the corpus-sync script, and never used to derive a
  web query. `ingest-url {url, title?}` is the only one with EGRESS — it persists
  ONE page into the corpus (SSRF-gated, redirects checked, open-licence pages
  only so the corpus stays redistributable) and is refused `403` unless the user
  armed Settings → "Allow Babs online database lookups". `remember-fact {fact,
  title?}` saves a durable fact; `memory-list` returns `[{slug, title}]` (check it
  before saving a duplicate, and to get slugs); `memory-forget {slug}` deletes one
  — forgetting changes what Babs believes in every future session, so it goes
  through the same approval as any other write. `learn-merge {slug}` folds a
  finished Background Learning session into the MAIN corpus.

- **RNA structure + RBS** — fold-rna (minimum-free-energy secondary
  structure: dot-bracket fold + ΔG in kcal/mol, pure-Python Turner-2004,
  no external dependency); cofold-rna (bound-state heterodimer ΔG of two
  strands, e.g. a 16S anti-SD : mRNA hybrid); rbs-strength (relative
  E. coli translation-initiation strength of a ribosome binding site —
  weighs SD:anti-SD match, 5'UTR occlusion, SD-to-start spacing, and the
  start codon; returns a RELATIVE ranking score, not an absolute rate);
  design-rbs (reverse-design a 5'UTR / Shine-Dalgarno + spacer to a
  target relative strength, with the achievable range + an on-target
  flag); assemble-operon (assemble a contiguous operon from a list of
  CDSs each with a target strength — context-aware RBS design + an
  element layout + per-gene achieved-vs-target report). DNA `T` read as
  `U`; ambiguous / over-length / bad input → 400.
- **HMM databases** — list / get / set-active / delete / add /
  download-hmm-database (the registry that backs hmmscan). `add`
  registers a custom `.hmm.gz` URL (mirrors the GUI "Add" form);
  `download` streams → decompresses → hmmpresses a catalog entry
  (builtin like `pfam-a` / `ncbifam`, or a custom one) into
  `<DATA_DIR>/hmm_databases/<id>/` exactly as a GUI download would, so
  `set-active` + `hmmscan` can then use it. `download` runs
  synchronously, so a large builtin (Pfam-A ~300 MB, NCBIfam ~600 MB)
  blocks the request for minutes; a 409 means a download for that id is
  already in flight. `delete` un-downloads the files but keeps the
  catalog entry so `download` can re-fetch it.
- **Plasmidsaurus** — list-plasmidsaurus-items (list your sequencing
  orders, most-recent first; `plasmidsaurus-items` is a back-compat
  alias) / download-plasmidsaurus (fetch a run's
  results zip by its 6-character item code over Plasmidsaurus's
  official OAuth2 REST API and import the run's `.gbk` assemblies into
  the library as new entries, tagged
  `source: plasmidsaurus:<code>:<sample>`, never overwriting existing
  entries). Credentials resolve env-first
  (`PLASMIDSAURUS_CLIENT_ID` / `_SECRET`) then the Settings values, and
  are deliberately NOT exposed through get-settings / set-setting.
  `download` runs synchronously and imports only `kind=results` (the
  archive that carries assemblies); no credentials → 400, a download /
  parse failure → 502, an archive with no samples → 422, and a 429 means
  the account's 10-requests-per-minute budget is spent.
  Each listed item is `{code, status, order_name, product_name,
  quantity, done_date, gross, has_results}`. **Choose a run by
  `has_results`, not by `status`** — roughly 40% of a real listing is
  shipping-label orders that are marked `status: "complete"` with a null
  `done_date` and sort to the top, but have no results zip (feeding one
  to `download-plasmidsaurus` 404s); `status` also takes `"canceled"`.
  `downloadable` counts the items that do have results. The server
  truncates the listing at 1000 items and offers no pagination, so a
  `truncated: true` response means older orders exist that the API
  cannot reach.
- **Custom enzymes + enzyme collections** — list / get / create /
  update / delete-custom-enzyme; list / get / create / update /
  delete-enzyme-collection; get / set-active-enzyme-collection.
- **Feature library** — list / get / create / update /
  delete-feature-library (reusable annotation snippets).
- **Feature presets** — `list-feature-presets` (the shipped, read-only
  catalogue of common plasmid elements; `{category}` / `{search}` filter
  it, `{include_sequence: true}` adds the bases, and the response also
  returns the `categories` list), `get-feature-preset` `{name}` for one
  full entry including its sequence and GenBank provenance, and
  `import-feature-preset` `{name}` or `{names:[…]}` to copy presets into
  the user's own feature library. An import replaces any entry with the
  same `(name, feature_type)` and reports how many it `replaced`; an
  unknown name makes the whole call a 404 and writes nothing, so a batch
  is never half-applied. Presets are NOT part of `list-feature-library`,
  which keeps returning only the user's own entries.
  `annotate-from-presets` matches the LOADED record against the feature
  library plus the catalogue, on both strands and across the origin, and is
  a dry run by default like `transfer-annotations` — pass `{"apply": true}`
  to append the features, `{"include_presets": false}` to scan only your own
  entries. It reports `already_annotated` for matches the record already
  carries under the same label over the same span, which are never offered
  again, so re-running is a no-op rather than a second copy of everything.
  Its `transfers` list is capped at `max_transfers` (default 500) so a
  tandem-repeat sequence cannot return thousands of rows; `total_found` and
  `truncated` report exactly what you are not being shown. `import-feature-preset`
  takes `name` OR `names`, never both, and collapses repeats so `count` always
  equals the number of entries the library gained.
- **Primer collections** — list-primer-collections, create-primer-collection
  (the primer-side parallel to create-collection), rename-primer-collection,
  delete-primer-collection (the last collection can't be removed; the active
  one auto-promotes + re-mirrors), set-active-primer-collection,
  move-primer (reassign a primer to another collection in one atomic call),
  plus per-primer CRUD. `list-primers` / `get-primer` / `delete-primer` accept
  a `{collection}` to read or prune just that collection's own primers (a
  collection is a real partition, not only for writes). `delete-primers`
  prunes a whole `{names:[…]}` batch in ONE locked save (the machine-friendly
  path for a large cleanup that would otherwise trip the rate limiter),
  reporting `removed` / `not_found` / `ambiguous`; a batch that removed
  nothing answers `404` / `ok: false` with an `error` that counts both (`409`
  when every name was ambiguous), and one that removed some answers `ok: true`
  with `partial: true` and a `warnings` line for what it left. `copy-plasmids`
  / `move-plasmids` answer the same way.
- **Data safety** — list-backups, restore-backup,
  list-pre-update-snapshots, restore-pre-update-snapshot. Restoring a
  settings backup whose active collection, primer library, parts bin or
  project no longer exists answers `409` with nothing changed; unsaved work
  is saved into its own collection before anything is restored, and the
  whole restore holds the data lock, so no save can land in the middle.
- **Settings** — get-settings, set-setting (allowlisted toggles only;
  a handful of settings — the Plasmidsaurus secret, the online-search
  egress gate `allow_online_search` — are deliberately NOT in the
  allowlist, so an agent can neither read nor flip them; those are
  human-armed via the GUI Settings dialog);
  get-feature-colors / set-feature-color (the per-feature-type render
  colour overrides).
- **Experiments lab notebook** — list / get / create / update /
  delete experiment entries; move-experiment (relocate an entry to
  another project in one atomic call); attach-experiment-image (attach a
  server-side file + embed it in the entry body — images AND data
  attachments such as `.ab1` / `.csv` / `.pdf` / GenBank / FASTA; the
  embedded reference is `![name](name)` for an image and `[name](name)`
  for anything else, and the reply says which via `is_image`. The endpoint
  keeps its `-image` name for compatibility); list / create /
  rename / delete projects; set active project (full notebook-layer CRUD).
  `create-experiment` / `list-experiments` / `delete-experiment` accept a
  `{project}` to file into / read / prune just that named project's own
  entries (a project is a real partition, not only for writes; deleting
  from the active project by name is refused — switch away first so the
  live mirror stays consistent).
- **Notebook read side** (2026-09-17) — the counterpart to the
  create/update/delete endpoints above, since nothing could previously
  find an entry by its contents:
  - `search-experiments` — `{query?, tags?, project?, all_projects?, limit?}`.
    Query terms must ALL appear, in the title, tags or body (AND across
    terms, OR across fields), so `"gibson failed"` finds the Gibson entry
    that mentions failure rather than every entry that says "gibson".
    `tags` narrows to entries carrying ALL the named tags and works on its
    own (the "every Gibson entry" browse). Each hit carries `project`,
    `matched_in` and a `snippet` of body context but **not** `body_md` —
    fetch that with `get-experiment`. Passing neither query nor tags is a
    400, not a dump of the notebook.
  - `experiment-backlinks` — `{ref, kind?, project?, all_projects?}`, the
    reverse of the `@plasmid` / `!action` / `&gel` body tags. `ref` takes a
    **list** so a plasmid can be asked for by both its library entry id
    (what the `Plasmid ref` button writes) and its display name (what a
    hand-typed ref usually is). Case-insensitive, and refs are
    re-extracted from each body rather than read from the stored
    `attached_*` xref — a hand-edited `experiments.json` can carry a stale
    index, and "nothing references that" when something does is the
    failure this endpoint exists to prevent.
  - `get-plasmid-protocol` — `{id, collection?}` → the plasmid's recorded
    construction steps, structured **and** as `markdown` ready to paste
    into an entry. Same `_history_build_steps` the History viewer renders,
    so it is what was recorded, not a reconstruction; only names that
    resolve in the library become `@refs`, so the markdown never plants a
    dangling one. Steps the program never saw (transformation, miniprep,
    sequencing) are absent by construction. 404 with no history, 422 on
    malformed history.
  - `experiment-step-template` — `{actions, heading?}` → a heading-per-step
    skeleton carrying each `!action` tag, in catalog order. The forward
    half of the pair: what you are about to do, rather than what was
    recorded. Unknown action ids are accepted (the catalog is curated, not
    enforced) and reported in `unknown_actions`.
  - `list-experiment-actions` — the curated bench-action catalog behind the
    `!<id>` tags, as `[{group, id, description}]`.
  - `export-experiment` — `{id?, project?, format?}` → the document TEXT
    (`markdown` or `html`), not a file, so the caller picks the
    destination. Markdown keeps each body verbatim so an export pastes
    back into a new entry with working refs; HTML is a standalone page
    with the stylesheet inlined.
  - `duplicate-experiment` — `{id, title?}`, a write. Image attachments are
    **not** copied (they live in a per-entry directory keyed by entry id,
    so carrying the paths over would point the copy into the original's
    directory); `attachments_copied` is always false and says so. Works on
    the **active** project only and refuses a `project` key rather than
    duplicating whatever shares that id there — unlike
    `create-experiment` it both reads and writes, and doing that against
    a non-active project while the live mirror holds the active one is
    the desync `set-active-experiment-project` exists to avoid.
- **Gels** — list / get / create / update / delete saved gel
  snapshots (in addition to simulate-gel for one-shot runs).
- **Protein motifs** — list (built-ins + user overrides),
  set (copy-on-write override), delete user overrides.
- **Entry vectors** — list (carries each vector's `role`; a grammar can
  hold several — Alpha1 / Alpha2 / Omega1 / Omega2 + the L0 vector), get,
  set, plus auto-detect across the full library and clear-for-grammar.
- **Active pointers** — every `set-active-*` (collection, codon-table,
  primer-collection, parts-bin, experiment-project, enzyme-collection,
  hmm-database) has a matching `get-active-*` so a client can read the
  current selection before changing it.
- **Utility** — check-primer-duplicates, capture-snapshot.
- **OT-2 / Opentrons** (liquid-handler control) — `ot2-compile` turns a
  plate-transfer plan (a pipette, labware on deck slots, and `from → to` well
  transfers with µL volumes) into an Opentrons Protocol API v2 `.py` text,
  validating it against a built-in deck catalog first (friendly labware aliases
  like `tiprack_300` / `eppi_24` / `plate_24` / `plate_96` / `reservoir_12`
  expand to canonical load names; volumes are range-checked against the
  pipette). For a **multi-step protocol designer**, pass an ordered `steps` list
  instead of `transfers` (it takes precedence): each step is a typed operation —
  `transfer` (with optional `mix_before` / `mix_after` / `blow_out` / `touch_tip`),
  `distribute` (one source → many wells), `consolidate` (many wells → one),
  `mix`, `delay`, `pause`, `comment` — compiled into one protocol; a control-only
  protocol (delay/pause/comment) needs no tips or labware. A plan that fails
  validation answers `422` / `ok: false` (with `valid: false` and its
  `errors`). Validation follows what the robot's own planner does: a tip rack
  from another pipette's family is an error, and so is a `distribute` /
  `consolidate` whose single trip will not fit the tip (the planner splits
  against the PIPETTE's maximum but loads against the TIP's, and silently
  ends the step at the first piece that does not fit); a multichannel step
  may start only in the rows its first channel reaches (row A on a 96-well
  plate, A–B on a 384); a custom definition with a well below the deck is
  refused. Wells
  are compiled the way the robot names them — `a1` and `A01` are written `A1`,
  a custom definition's own key is used verbatim — so a plan the validator
  passes is one the robot can load. `ot2-analyze` uploads a plan/protocol to the robot for its **built-in
  simulate** — server-side validation with NO motion, the pre-flight for a run.
  `ot2-status` returns a full state snapshot for **crash monitoring**:
  reachability + versions, pipette OK flags + volume specs, motor engagement,
  deck / instrument calibration health, module sensor readings, robot settings
  (variables), the status light, and — for the active run or a passed
  `run_id` — live run / command state, with a detected `faults` list + an `ok`
  verdict. `ot2-run` (a **write** endpoint) is GATED physical actuation: it moves
  the gantry only when the robot's own analysis passes AND the body carries
  `{"confirm": true}`, refuses on an already-faulted robot (pre-flight check), and
  monitors state throughout — halting and reporting the instant a fault is
  detected. Pass `{"wait": false}` to start a run and poll `ot2-status` with the
  returned `run_id` to watch it live. `ok` says whether the run did its job: a
  dry call (no `confirm`) whose analysis passed is `ok`, as is a run started
  with `wait: false`; a pre-flight refusal, a crash, or a run that ENDED
  `failed` or `stopped` is `422` / `ok: false` with an `error` naming why — so
  check `ok`, not just `ran`. (`422`, never a 5xx or 409: those are statuses a
  client may retry on its own, and a retried run moves the robot again. A run
  that started but could not be watched to its end — the play request
  failed, the monitor lost it, or it outran its timeout — is STOPPED and
  answers `422` with `crashed: true` the same way. `ot2-position-check`
  answers the same way.) The pre-run pipette check
  compares GENERATIONS: a
  protocol for a GEN2 pipette is refused with a GEN1 attached, while a GEN1
  protocol runs on the GEN2 pipette Opentrons back-maps it to. Host comes from
  the body's `host` or the persisted `ot2_host` setting, and must be on your
  own network: a private, link-local (USB) or CGNAT/Tailscale address, an
  address in one of this computer's own /24 networks (campus networks hand
  lab devices public addresses), or the address saved in AUTOLAB (`ot2_host`
  — only AUTOLAB writes it). Every name, `.local` included, is resolved,
  every address it returns must pass, and the request is sent to the address
  that was checked — never looked up a second time. A cloud-metadata address
  is refused however it is written (an IPv4-mapped `::ffff:` form included),
  as is a URL carrying credentials or a port outside 1–65535. (Compiler + client live in the app-free
  `splicecraft_opentrons` sibling.)
- **OT-2 run control** — `ot2-run-control` (a **write** endpoint) pauses, resumes,
  or stops a live run: `{"host": ..., "action": "pause"|"resume"|"stop",
  "run_id"?: ...}`. The active run is resolved automatically when `run_id` is
  omitted (so it also controls a run started from the Opentrons App). This is the
  manual counterpart to `ot2-run`'s automatic stop-on-fault.
- **OT-2 concentration normalisation** — `ot2-normalize` computes per-sample
  volumes to equalise DNA mass or hit a target concentration, no robot needed:
  `{"items": [{"name", "well", "concentration" (ng/µL)}, ...], "target_ng" | ("target_conc" + "final_volume"), "pipette"?}`.
  It respects the pipette's volume floor/ceiling (from `pipette`, or `min_vol` /
  `max_vol`), flags samples too dilute/concentrated to reach target (never dropping
  them), and — when a `src` + `dst` (+ `dst_wells` or `dst_labware`, optional
  `diluent_ref`) is given — also returns compile-ready `steps`. One labware may
  be both `src` and `dst` only when `dst_wells` is given and none of those
  wells is a sample or the diluent well (`400` otherwise — the defaults fill
  from A1, over samples still to be read).
- **OT-2 plate map** — `ot2-plate-map` maps a plasmid collection onto a labware's
  wells row-major (the "a plate IS a collection" link): `{"collection", "labware"}`
  → `{map: {well: {id, name}}, wells, n, overflow}`. Feed the result into
  `ot2-normalize` / `ot2-compile` to cherry-pick or replate a collection by
  identity.
- **OT-2 calibration + safe motion** — `ot2-calibration` reports deck +
  per-pipette-offset + tip-length calibration with a `ready` verdict and a
  `needs_calibration` list (calibrate those in the Opentrons App). `ot2-home` (a
  **write** endpoint) homes the gantry — safe motion, no plunger. `ot2-position-check`
  (a **write** endpoint) runs a MOTION-ONLY tour: the gantry moves to the *top* of
  each loaded labware's wells so you can verify alignment / labware offsets, emitting
  no aspirate/dispense/tip pick-up (the plunger never actuates, the pipette never
  descends into a well). It is gated like a run (`{"confirm": true}` + a clean
  analysis) but does NOT require a calibrated pipette — that is what you're checking;
  deck calibration + reachability are still required. Body:
  `{host, plan|<plan fields>, wells?, confirm?, wait?}`. **Labware offsets:** add an
  `offset: {x, y, z}` (mm) to any plan labware entry and it is applied to the run's
  motion (`ot2-run` and `ot2-position-check` both honour it); a liquid `ot2-run` also
  now refuses to move an *uncalibrated* pipette, refuses when the **attached pipette
  does not match** what the protocol loads (`reason: "pipette-mismatch"`), and refuses
  when the robot's **door-safety switch is enabled and the door is open** (`reason:
  "door-open"`). During a monitored run the robot's rail lights come on as a
  "gantry moving" indicator and are restored afterwards, and a run that overruns the
  poll timeout is **stopped on the robot** rather than left moving.
- **OT-2 lights + motor disengage** (zero-/low-motion controls) — `ot2-lights` (a
  **write** endpoint) toggles the robot's rail (status) lights: `{host, "on":
  true|false}` (defaults on) — handy to blink/identify the robot or turn the lights
  off for the night. `ot2-disengage` (a **write** endpoint) de-energises the gantry
  motors so the carriage can be moved by hand (`{host, "axes"?: [...]}`, all six axes
  by default); it homes nothing and is refused while a run is active. `ot2-status`
  additionally surfaces the **door** state (`{"status", "required_closed"}`) alongside
  lights, motors, and modules.
- **OT-2 protocol + custom-labware libraries** — manage saved designs and custom
  labware headlessly, mirroring the AUTOLAB Deck/Labware tabs. Protocols:
  `list-protocols`, `get-protocol` (returns a compile-able `plan`), `save-protocol`
  (**write**; validates the plan first), `delete-protocol`, and the collection set
  `list-protocol-collections` / `create-protocol-collection` (**write**) /
  `delete-protocol-collection` (**write**). Custom labware: the identical
  `list/get/save/delete-custom-labware` + `list/create/delete-labware-collection`
  set, where `save-custom-labware` takes an Opentrons `definition` (a `wells` map).
  Both stores are single-file collections that embed their items — backed up,
  restorable, and swept by Master Delete like every other user-data store.

Call `/tools` for the live discovery endpoint. Each entry is
`{name, method, write, doc, doc_full}` — `doc_full` is the endpoint's
COMPLETE docstring, which documents the request body (required / optional
keys, aliases, enums, size caps), so a client forms a correct call in one
round-trip instead of N trial-and-error 400s. `/tools` is authoritative:
it lists every registered endpoint, including the ~49 app-coupled
handlers that live in `splicecraft.py` rather than `splicecraft_agent.py`
— don't rely on grepping a single file for `@_agent_endpoint`.

Every **authenticated** success response also carries a uniform **`ok: true`**
flag and a predictable **`data`** field. `ok: true` is one consistent success
signal across every endpoint (an explicit `ok` a handler already sets is never
overwritten). `data` saves you from knowing each endpoint's ad-hoc key (`seq` /
`library` / `sites` / `matches` / …): it's the result with the envelope/metadata
stripped, unwrapped to the bare value when there's a single content key (a
scalar or list lands directly under `data`). The original keys stay too, so
it's a superset — read whichever you prefer.

One caveat worth knowing: on the endpoints that MODEL a reaction —
`simulate-golden-gate`, `simulate-gibson` — `ok` reports whether the
**assembly** worked, not whether the request was served. A dry run that
can't close a circle answers HTTP 200 with `ok: false` and the reason in
`result.errors`. (These two used to stamp `ok: true` regardless, so the
natural `if r["ok"]:` accepted a design that cannot be built and callers
had to write `r.get("ok") or r.get("result", {}).get("ok")`.)

Two **unauthenticated** endpoints are served *ahead* of the envelope wrapper
and stay bare — neither carries `data` (nor the `_stale` warning): the
`/healthz` readiness probe (`{ok, status, version, headless}`) and the
`/tools` self-describe (`{endpoints: [...]}`, read it by name). Both answer
before a client has a token, so they short-circuit the authenticated response
pipeline; for `/tools`, mirroring its `doc_full` corpus under `data` would
also needlessly double the API's largest response. Don't blindly read
`resp["data"]` on these two — check by endpoint.

- **CRISPR** — `list-cas-variants` (SpCas9 / SaCas9 / LbCas12a — three
  distinct PAM geometries), `design-guides` (over the loaded plasmid or a
  supplied `sequence`; `region` aims at an exon, `variant` picks the nuclease),
  `score-guide`, `guide-offtargets`, `guide-cloning-oligos`. Guides come back in
  FORWARD-strand coordinates with an explicit strand and cut site, wrap-aware on
  a circular plasmid. **Two deliberate non-claims:** the triage reports NAMED
  properties (GC window, Pol III terminator, homopolymers,
  self-complementarity) and not a fabricated efficiency score; and off-target
  search covers only what you hand it, with `offtarget_searched_bp` in every
  response so silence can never read as "specific". There is no genome-wide
  prediction here.
- **Protein-level editing** — `plan-residue-edit` (what changing residue *n*
  would do, without doing it) and `edit-residue` (do it). The residue is
  resolved to its three plasmid base pairs through strand, `/codon_start`,
  introns and the origin — the mapping you would otherwise do by hand — and the
  codon is chosen from the active codon-usage table, with `alternatives` ranked
  so you can pick one that avoids a site you care about. `silent`,
  `creates_stop` and `removes_stop` are reported rather than refused.
  `initiator` is set when residue 1 is a start codon the CDS's annotation
  reads as Met (see `find-feature`): `wt_aa` is then `M`, and changing the
  codon to ATG is a silent edit.
- **Reads + verification** — `read-consensus` combines every read STORED on a
  plasmid: for each difference, how many reads that COVER that base show it and
  how many looked and disagreed (`confirmed` / `single_read` / `unconfirmed` /
  `clean`). `uncovered_spans` is part of the answer, because a plasmid can read
  "100% identity" over the half that was sequenced. `verify-against-reads` takes
  read SEQUENCES instead, and its optional `read_quality` (per-read Phred
  arrays) splits each read's differences into the ones it supports and the ones
  sitting in unreliable basecalls. Each read is judged on its COVERED identity —
  over the stretch of reference the read spans, where a deletion or insertion
  inside that stretch counts against it and an `N` is a no-call, neither match
  nor mismatch. A read that covers nothing (all `N`) is never a pass, and a
  read that does not span the whole plasmid is marked `partial` (with its
  `covered_bp` / `coverage_pct`): a pass then confirms that window only. A
  read with fewer than 20 aligned bases is `too_short` and never passes. The
  reference's topology is honoured throughout — a linear reference is never
  judged as a circle — and a local-mode read's clipped ends are not counted
  as deletions.
- **Ordering + the bench** — `export-primers` writes an order sheet
  (`generic` / `idt` / `plate` — row-major A1..A12 wells, rolling onto plate 2
  past 96 oligos), closing the gap where only the GUI could order oligos.
  `pcr-program` + `list-polymerases` give the thermocycler block for a product
  length and primer Tms; every number carries the rule that produced it, and
  with no Tm supplied the program says it is using a conventional default
  instead of presenting an invented number as derived.

## Security posture

- **Bearer-token auth on EVERY endpoint**, reads included — the token
  check fires before the handler lookup, so an unauthenticated probe for
  any path gets a uniform `401` rather than a `404` that would confirm
  which endpoints exist. (The two pre-envelope endpoints below,
  `health` and `tools`, are the documented exceptions: a client has to
  be able to self-describe before it holds a token.) This section used
  to say reads were unauthenticated; they have not been for some time,
  and a reader planning around that would have mis-set their own
  expectations about what a co-resident local process can read.
- **Localhost only** (`127.0.0.1`) — single-tenant by design. Do not
  expose on a LAN.
- **Inputs are length-, range-, and shape-validated at the boundary.**
  A request body must arrive with a `Content-Length`: one sent with a
  transfer coding (chunked, gzip…) answers `411`, as does an HTTP/1.0
  request that carries a body without one — it used to run as an empty
  request. With an idempotency key, a body too deeply nested to fingerprint
  is refused (`400`) rather than risking a replay of another request's
  answer.
- **Symlink refusal**: write paths go through
  `_check_agent_write_path` which walks the full ancestor chain via
  `resolve()` divergence + per-segment `is_symlink()`. Pre-fix this
  only checked the immediate parent — see the symlink ancestor-chain
  regression in [`docs/invariants.md`](invariants.md).
- **Read-dir traversal** uses `lstat` + `S_ISDIR` to refuse
  directory-symlink escapes.
- **Per-handler size caps** — `_h_load_file` 50 MB
  (`force=true` override), agent paths capped via
  `_safe_file_size_check`, manifest reads capped at
  `_PRE_UPDATE_MANIFEST_MAX_BYTES`.
- **Secrets stay out of the agent surface.** The Plasmidsaurus Client
  Secret is excluded from the settings allowlist (so get-settings can't
  read it back and set-setting can't change it) and is redacted to
  `<redacted>` in the change log + event stream. The Plasmidsaurus
  download endpoint streams over the same hardened, HTTPS-only,
  size-capped fetch path as the rest of the network layer and verifies
  a zip-archive magic header before writing to disk.

## Cross-collection lookups

`load-entry` resolves **across collections**: it checks the active
collection first, then the others. A unique cross-collection hit loads
(and the response names its `collection`); an ambiguous name returns
`409` listing the holders; pass `collection` to pin the search.
`search-library` is cross-collection too, and the `id` it returns is a
valid `load-entry` key.

The mutation endpoints `rename-plasmid`, `set-plasmid-status`,
`delete-from-library`, and the lookup in `diff-plasmid` operate on the
**active collection only** (via `_load_library()`). To target a plasmid
elsewhere, `set-active-collection` to its home first — `search-library`
shows which collection holds it.

`delete-from-library` is **undoable by the human**: every entry it removes
is pushed onto the running app's library-undo stack, so the user can press
`u` (or Ctrl+Z) in the library panel and take back what an agent deleted —
restored to its original slot, newest first, for the rest of that session.
It matches by NAME, so one call can remove several entries; they come back
in their original order. Nothing is recorded if the save fails. Managing
the library IS an agent's job; what stays unreachable over HTTP is the
whole-data **master wipe**, which has no endpoint by design (enforced by
`tests/test_master_delete.py::test_no_agent_endpoint_exposes_wipe`).

`transfer-annotations` resolves its `source_id` **across collections**
(like `load-entry`): the active collection first, then the others, by
display **name or id**. Pass `source_collection` to pin the search, and an
ambiguous key across collections returns `409`. So a backbone in a
non-active collection transfers onto the loaded record without a temp copy.
To apply (not preview) pass `apply: true` — or the legacy `dry_run: false`;
the default is a dry run.

**Read the `skipped` block.** The response carries
`skipped: {min_len, considered, below_min_len, shortest_skipped_len,
no_sequence_match}`, because the default `min_len` of 30 drops every
promoter, operator, RBS, fusion overhang and restriction site *before* the
search runs — so a small `count` usually means the threshold, not a
mismatch. Lowering `min_len` picks those up at the cost of short features
matching in many places. For "same backbone, one block swapped" prefer
`carry_annotations` on the assembly endpoint, which maps coordinates
instead of searching sequence. A palindromic feature (any restriction
site) is searched forward only — scanning both strands used to emit the
identical span twice with opposite arrows.

`list-library` lists the active-collection library by default; pass
`{collection}` to scope it to any one collection's plasmids without
switching the active collection (`404` if there's no such collection).

### Name integrity

`load-entry` stamps the entry's stored display name onto the loaded
record, so a later `save` / `add-current-to-library` round-trips the real
name (spaces preserved) instead of the underscored GenBank LOCUS. A
record pulled in with `fetch` (inspect-only, `saved:false`) is **not**
auto-filed on `save` — creating a new library entry requires an explicit
`{create:true}` so an inspection can't silently pollute the active
collection. Unknown body keys are echoed back under `ignored` rather than
silently dropped. A stricter guard applies to a closed set of **routing /
selection** params (`collection`, `source_collection`, `bin`, `parts_bin`,
`enzyme`, `enzymes`, `orientation`, `rename`, `carry_annotations`): passing
one to an endpoint that
doesn't accept it is a hard **400** (it would change *where/how* the op applies,
so a silent drop is the SC-D footgun) — while any *other* unknown key stays soft
for forward-compat. The check is central, on every WRITE endpoint: the
dispatcher derives the keys each handler actually reads from its own code
(following `payload` into the helpers it is passed to) and refuses a routing
key the handler would ignore with `400` and `unsupported: [...]` —
`delete-part {"parts_bin": "Archive"}` used to answer `200` and delete from
the ACTIVE bin. The analysis follows a key read inside a constant-tuple loop
and reads an f-string key as its shape — so `traditional-clone` is known to
take `vector_collection` / `insert_collection` and a bare `collection` is
refused, rather than the product quietly landing in the active collection. A
routing key sent EMPTY (`null`, `""`, `false`, `[]`, `{}`) is inert and never
refused — typed clients send optional fields that way. A handler whose reads
cannot be fully traced is left alone rather than guessed at.

`load-file` and `add-current-to-library` accept optional `{id, name}`
overrides to stamp the record's identity directly. This is the
deterministic fix for a **batch build loop**: a from-scratch `.dna`
carries no name packet, so without an override every such file loads as
id `Cloned` and a loop of `add-current` calls silently REPLACES each
prior entry (the library dedups by id) down to a single survivor. Pass a
unique `id` (scrubbed to `[A-Za-z0-9_-]`, 32-char cap) per build — and a
spaced `name` for the display name — and all N persist. The override is
applied to a copy, so the live canvas record is untouched.

## Concurrency

- Heavy ops (BLAST build, BLAST search, HMMscan, alignment) run in
  `@work(thread=True)` workers; the API returns immediately with a
  status the client can poll, OR blocks the request until the worker
  completes — endpoint-specific.
- A fixed set of CPU- or network-heavy endpoints (`blast`, `hmmscan`,
  `multi-align`, `optimize-protein`, `verify-against-reads`,
  `analyse-read-heterogeneity`, `fetch`, `blast-online`, `design-guides`,
  `scan-terminators`, `map-transcription`, … — `_AGENT_HEAVY_ENDPOINTS`) share
  a concurrency cap of half the CPU cores (at least 2). A call past the cap is
  answered at once with `503` and `retry_after: 2` rather than queued; nothing
  ran, so retrying it is safe.
- The agent server uses `_agent_save_or_500(save_fn, label)` for
  every `_save_*` call so an OSError / RuntimeError becomes a 500 +
  in-app notify, not a silent in-memory / disk desync.

## Discovery + introspection

```bash
splicecraft-cli tools             # list every endpoint + one-line doc
splicecraft-cli status            # current record snapshot
splicecraft-cli features          # features on the loaded record
splicecraft-cli call <endpoint> --json '{...}'   # call ANY endpoint
```

`call` is a generic passthrough that reuses the same token / host / port
plumbing as the named subcommands (method defaults to POST when `--json`
is given, else GET), so a client never has to hand-roll HTTP against the
private `_request`. An HTTP error surfaces as structured JSON (with
`http_code`) plus a non-zero exit, instead of a hard exit that would kill
a batch mid-run.

See the [CLI sidecar](cli.md) for the full convenience wrapper.
