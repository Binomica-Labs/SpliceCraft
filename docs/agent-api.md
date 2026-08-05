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
```

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

~220 endpoints across:

- **Records** — `new-plasmid` (create from a raw sequence, the Ctrl+N
  flow), get / set sequence, add / update / delete features (with the full
  arrow type — forward ▶ / reverse ◀ / arrowless / double-stranded ◀▶ —
  plus colour and arbitrary GenBank qualifiers, matching the Insert/Edit
  Feature dialog; `add-features` inserts MANY at once under one lock + one
  dirty-check), list features (`features`, alias `list-features`), find ORFs
  (length cutoff in AMINO ACIDS —
  `min_aa`; `min_length`/`min_bp` are rejected so a bp-vs-aa mix-up can't
  silently return a default-length result. Each ORF carries `start`, `end`,
  `strand`, `length_aa`, `aa_seq`, plus `nt_len` — the exact coding length in
  bases INCLUDING the stop codon — and `exceeds_one_lap`. **Take length from
  `nt_len` / `length_aa`, never from `end - start`:** on a circular plasmid a
  frame can run a full lap or more without an in-frame stop, and no
  start/end pair on an n-bp circle can express that. Such an ORF sets
  `exceeds_one_lap: true` and has its span pinned to the near-full circle),
  `undo` / `redo` the last edit,
  `discard-changes` (revert the canvas to its library-stored copy / clear a
  stuck-dirty flag so the next load / new-plasmid proceeds without `force`),
  transfer annotations, apply GFF3 features to the loaded record
  (`apply-gff3`).
- **Files** — load (chromosome-scale safe via the path-based loader;
  supports `.gb` / `.gbk` / `.genbank` / `.dna` / `.embl` /
  FASTA / `.ab1` / single-record `.fastq` / `.gff3`),
  export GenBank / GFF3 / FASTA / EMBL / CommercialSaaS `.dna`
  (symlink-guarded), bulk import a folder, bulk export a
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
  silently renaming or dropping), create / rename / delete collections,
  get / set the active collection, list / set plasmid statuses.
- **Parts** — list-parts, get-part, create-part / update-part (full
  parity with the Part editor — name / grammar / type / level / position /
  overhangs / sequence plus `backbone`, selection `marker`, and
  `fwd_primer` / `rev_primer` domestication primers, Tms derived),
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
  orientation, and save the product — refuses to guess when more than one
  ligation is possible, pick with `vector_frag_idx` / `insert_frag_idx`. Pass
  `insert_circular:true` to cut a cassette OUT of a *plasmid* insert — e.g. an
  Ω multigene into a binary vector — so its two digest fragments are both
  offered instead of treating the insert as a linear PCR product. Pass
  `carry_annotations:true` to lift the `vector_name` / `insert_name` entries'
  own features onto the ligated product — split at the cuts and shifted to
  product coordinates — so the clone isn't feature-bare; a side whose passed
  sequence doesn't exactly match its named entry is reported in
  `carry_warnings` rather than mis-placed),
  golden-gate-assemble / simulate-golden-gate (Type IIS — BsaI /
  BsmBI / BbsI / SapI / Esp3I — overhang-directed N-part assembly: parts in
  any order, chained by their 4-nt overhangs into a circle, with a
  unique-overhang + no-residual-site fidelity check),
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
  `verified` flag)), design-gb-part (Golden Braid / MoClo domestication
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
  can't clone doesn't abort the rest),
  assemble-into-entry-vector (the multi-source level-up clone: chain L0
  parts into an α/L1 TU, or TUs into an Ω/L2 module, by their fusion
  overhangs and ligate into the configured `role` acceptor at the level-up
  Type IIS enzyme — pass `{sources, grammar?, source_level?, role, name,
  collection?}` (active collection by default, or the named one).
  Also **fails loud** (422) on a missing acceptor or a chain that can't
  close; the L2→binary hop is a known engine gap), design-primers
  (Primer3 primer-pair design over a region: `detection` picks a
  diagnostic amplicon, `cloning` appends RE-site tails, `generic`
  returns binding-only primers — no tails / overhangs), check-primer
  (a SINGLE oligo vs a template: Tm, GC%, and every 3'-anchored binding
  site — both strands, wrap-aware on a circular template, the tail
  scored as mismatches — to confirm a designed primer binds exactly once
  before saving), optimize-protein
  (codon-optimise an AA sequence to a chosen table; optional `stops`
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
- **Digest** — digest (cut a RAW sequence with named enzymes and report
  the cuts + resulting fragments with their **overhangs** — overhang-aware
  QC for a Golden-Braid / restriction junction without loading the
  sequence onto the canvas; `circular` defaults true, a singular `enzyme`
  is accepted, and names the catalog doesn't know are reported under
  `unknown_enzymes` rather than silently dropped).
- **Alignment** — diff-plasmid (one target, circular rotation
  auto-detected), multi-align (batch: the loaded plasmid or a given
  sequence vs many targets at once — the Alt+A overlay; rotation-aware
  per target, circular auto-detected from each target's topology or
  forced with a `circular` flag, with `picked_rotation` reported per
  row), list-plasmidsaurus-members, align-plasmidsaurus-zip,
  verify-against-reads (the "I built X, I got reads back — do they match?"
  check: a bare or library-resolved `reference` vs a list of raw `reads`
  (Nanopore / Sanger / a consensus), each aligned rotation + RC-aware; returns
  per-read identity% + a `match`/`mismatch` `verdict` against `min_identity`).
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
  the scan.
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
  reporting `removed` / `not_found` / `ambiguous`.
- **Data safety** — list-backups, restore-backup,
  list-pre-update-snapshots, restore-pre-update-snapshot.
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
  server-side image file + embed it in the entry body); list / create /
  rename / delete projects; set active project (full notebook-layer CRUD).
  `create-experiment` / `list-experiments` / `delete-experiment` accept a
  `{project}` to file into / read / prune just that named project's own
  entries (a project is a real partition, not only for writes; deleting
  from the active project by name is refused — switch away first so the
  live mirror stays consistent).
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
  protocol (delay/pause/comment) needs no tips or labware. `ot2-analyze` uploads a plan/protocol to the robot for its **built-in
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
  returned `run_id` to watch it live. Host comes from the body's `host` or the
  persisted `ot2_host` setting. (Compiler + client live in the app-free
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
  `diluent_ref`) is given — also returns compile-ready `steps`.
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

Two **unauthenticated** endpoints are served *ahead* of the envelope wrapper
and stay bare — neither carries `data` (nor the `_stale` warning): the
`/healthz` readiness probe (`{ok, status, version, headless}`) and the
`/tools` self-describe (`{endpoints: [...]}`, read it by name). Both answer
before a client has a token, so they short-circuit the authenticated response
pipeline; for `/tools`, mirroring its `doc_full` corpus under `data` would
also needlessly double the API's largest response. Don't blindly read
`resp["data"]` on these two — check by endpoint.

## Security posture

- **Bearer-token auth** on every write endpoint; reads are
  unauthenticated to keep scripted introspection ergonomic.
- **Localhost only** (`127.0.0.1`) — single-tenant by design. Do not
  expose on a LAN.
- **Inputs are length-, range-, and shape-validated at the boundary.**
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
`enzyme`, `enzymes`, `orientation`, `rename`): passing one to an endpoint that
doesn't accept it is a hard **400** (it would change *where/how* the op applies,
so a silent drop is the SC-D footgun) — while any *other* unknown key stays soft
for forward-compat. Currently enforced on the `move-primer` / `move-part` /
`move-experiment` routers.

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
