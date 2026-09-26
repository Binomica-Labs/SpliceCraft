# Agent API — complete endpoint index

Every registered endpoint, generated from the live handler registry
(`_state._AGENT_HANDLERS`) so it cannot drift from the code. `R` reads,
`W` writes (a write also takes the canvas dirty-guard). The prose guide,
with request/response shapes and worked examples, is `agent-api.md`;
this file exists so no endpoint is invisible.

Regenerate with `python3 scripts/gen_agent_endpoint_index.py`.

| Endpoint | R/W | What it does |
|---|---|---|
| `add-codon-table` | W | Add or replace a codon-usage table. Body either fetches from Kazusa or stamps a raw dict directly: |
| `add-current-to-library` | W | Persist the currently-loaded record into the active collection (panel equivalent of `Alt+K`). |
| `add-feature` | W | Add a feature — full parity with the AddFeatureModal. Body: `{start, end, label? (or name), type?, strand?, color?, qualifiers?}`. |
| `add-features` | W | Batch feature insert — add MANY features under ONE UI-thread lock and ONE stale-load-counter check (agent-API feedback P1-10), so annotating a… |
| `add-hmm-database` | W | Register a CUSTOM HMM database in the catalog — the headless equivalent of the GUI catalog modal's "Add" form (`HmmDbAddEditModal`). Body: `{name,… |
| `align-plasmidsaurus-zip` | R | Align a Plasmidsaurus zip member against a library plasmid. Body: `{path, member, target_id?, target_name?, mode?, circular?}`. |
| `amplify-feature` | R | Design PCR primers that amplify ONE feature of the LOADED plasmid end to end, then check where each primer actually binds. Body: `{name \| idx,… |
| `analyse-cds` | R | Read the codon usage of an EXISTING CDS — the read half of `optimize-protein`, which until now only wrote. |
| `analyse-read-heterogeneity` | R | Per-base ALLELE FRACTIONS from raw reads — "is this culture one thing?" |
| `annotate-from-presets` | W | Annotate the LOADED record from the feature library and the built-in preset catalogue by exact sequence match, on both strands and across the… |
| `apply-gff3` | W | Apply features from a sequence-less GFF3 file onto the currently-loaded record. Body: `{path}`. The GFF3 file MUST NOT carry a ##FASTA section —… |
| `assemble-into-entry-vector` | W | Assemble multiple sources into a configured Golden-Braid / MoClo acceptor — the multi-source level-up clone (L0 parts → α/L1 TU, or TUs → Ω/L2… |
| `assemble-operon` | R | Assemble a contiguous bacterial operon from CDSs, each RBS reverse-designed (context-aware) to its target relative strength. Body: `{genes: [{cds,… |
| `attach-experiment-image` | W | Attach a file to a lab-notebook experiment and embed a markdown reference in its body — the agent counterpart to the experiment screen's "Attach… |
| `auto-detect-entry-vectors` | W | Run the Type-IIS-digest auto-detection pass across every library entry and bind newly-detected acceptors to their grammar roles. Body:… |
| `blast` | R | In-process BLASTN / BLASTP against the user's plasmid collections. Body: `{query, program?, collections?, max_hits?, six_frame?, backend?}`.… |
| `blast-online` | R | Remote NCBI BLAST via the public URL API. Body: `{query, program?, database?, max_hits?}`. |
| `bulk-export-collection` | W | Write every plasmid in `collection` to `dir` as files of `format`. Body: `{collection, format, dir}`. `format` ∈ {"genbank", "embl", "fasta",… |
| `bulk-import-folder` | W | Import every `.dna` / `.gb` / `.gbk` / `.genbank` file in a folder into a target collection. Body: `{folder, collection}`. The collection is… |
| `capture-snapshot` | W | Write a Markdown UI snapshot (same content as Alt+D) and return the file path. Useful for an agent assembling a bug report: snapshot first, then… |
| `check-primer` | R | Check ONE primer against a template: melting temp, GC%, and every 3'-anchored binding site (both strands, wrap-aware on a circular template).… |
| `check-primer-duplicates` | R | Scan the primer library for duplicate sequences. Returns a list of groups where multiple entries share the same canonical primer sequence — useful… |
| `classify-part` | R | Classify a candidate part by digest-overhang matching, without persisting. Body: `{sequence, circular?: bool = true, features?: list}` — same… |
| `clear-entry-vectors-for-grammar` | W | Drop every entry-vector binding for one grammar id. Body: `{grammar_id: str}`. Used by grammar-delete flows; agents rarely need this directly but… |
| `cofold-rna` | R | Bound-state heterodimer ΔG of two strands. Body: `{seq_a, seq_b}`. Returns `{ok, dg}` — the free energy (kcal/mol) of strand B bound to strand A… |
| `copy-plasmid` | W | Copy a library plasmid into a collection, name + features intact. Body: `{name\|id, to, from?}`. |
| `copy-plasmids` | W | Batch-copy library plasmids into a collection in ONE locked save — the machine-friendly bulk path (agent-API feedback P1-6), so a large… |
| `create-collection` | W | Create a new (empty) collection. Body: `{name, folder?}`. If `folder` is supplied, every importable file in that directory is bulk-imported into… |
| `create-custom-enzyme` | W | Add a new custom enzyme. Body: `{name, site, fwd_cut, rev_cut, type?, supplier?}`. Returns 409 if `name` collides with an existing built-in OR… |
| `create-enzyme-collection` | W | Create a new enzyme collection. Body: `{name, enzymes?: [str, ...] = []}`. Returns 409 if a collection with the same name already exists. Unknown… |
| `create-experiment` | W | Create a new notebook entry in the active project. Body: `{title?: str = "Untitled entry", body_md?: str = "", tags?: list[str] = []}`. |
| `create-experiment-project` | W | Create a new (empty) experiment project. Body: `{name: str, description?: str = ""}`. Active-project pointer is NOT modified by this endpoint —… |
| `create-feature-library` | W | Add a feature snippet to the persistent feature library. Body: see `_agent_feature_dict` for the schema. Returns 409 if a feature with the same… |
| `create-gel` | W | Save a new gel snapshot. Body: `{name: str, lanes: list[dict], agarose_pct?: float = 1.0, notes?: str = ""}`. Routes through… |
| `create-grammar` | W | Create a new custom grammar. Body: every field in the schema listed above. Returns 400 on validation failure, 409 if a grammar with the same id… |
| `create-labware-collection` | W | Create an empty named custom-labware collection. Body: `{name}`. |
| `create-part` | W | Add a part to a parts bin. Body: see `_agent_part_dict` for the part schema, plus optional `{bin}`. |
| `create-parts-bin` | W | Create a new (empty) named parts bin (SC-M) — the parts-side parallel to `create-collection` / `create-primer-collection`. Body: `{name,… |
| `create-primer` | W | Add a primer to the library. Body: `{name, sequence, source?, status?, type?, tm?, notes?, collection?}`. |
| `create-primer-collection` | W | Create a new (empty) named primer collection. Body: `{name, description?}`. |
| `create-protocol-collection` | W | Create an empty named protocol collection. Body: `{name}`. |
| `delete-codon-table` | W | Remove a codon-usage table by `taxid` or `name`. Built-in tables (source='builtin') cannot be removed. |
| `delete-collection` | W | Delete a collection by exact name. Body: `{name}`. If the deleted collection is the active one, the next collection becomes active (reported as… |
| `delete-custom-enzyme` | W | Delete a custom enzyme by `name`. Built-in NEB enzymes are refused. Any enzyme collection that referenced the deleted name keeps the stale row —… |
| `delete-custom-labware` | W | Delete a custom labware — the FIRST match for `name` (case-insensitive), optionally scoped to `collection`. Body: `{name, collection?}`. |
| `delete-enzyme-collection` | W | Delete an enzyme collection by `name`. If the deleted collection was the active one, clears the active pointer. |
| `delete-experiment` | W | Delete a notebook entry by id, optionally scoped to a `{project}`. Without `project` it deletes from the active project; with `project` it removes… |
| `delete-experiment-project` | W | Delete a project + its entries. Body: `{name: str}`. Refuses to delete the last remaining project (matches the GUI last-project guard). If the… |
| `delete-feature` | W | Delete feature at `idx` (= row index in `features`). Body: `{idx}`. Returns the deleted feature's label + bp range so the agent can undo via its… |
| `delete-feature-library` | W | Remove a feature-library entry by (name, feature_type). Body: `{name, feature_type?}`. Named `delete-feature-library` so it doesn't collide with… |
| `delete-from-library` | W | Remove a plasmid library entry by `name`. Body: `{name}`. Mirrors `_save_library`'s sync-to-active-collection so the collection bookkeeping stays… |
| `delete-gel` | W | Delete a saved gel snapshot by id. |
| `delete-grammar` | W | Delete a custom grammar by id. Built-ins are refused. Also clears any entry vector bound to the grammar so a future create-grammar with the same… |
| `delete-hmm-database` | W | Delete the DOWNLOADED files for an HMM database by `id` (un-downloads it; the catalog entry stays so it can be re-fetched). Files removal is… |
| `delete-labware-collection` | W | Delete a custom-labware collection AND every labware in it. Body: `{name}`. |
| `delete-part` | W | Remove a part by `name` (optionally `grammar`), optionally scoped to a `{bin}`. Without `bin` it deletes from the active bin; with `bin` it… |
| `delete-parts-bin` | W | Delete a named parts bin + its parts. Body: `{name}`. Refuses to delete the last remaining bin (one must always exist). Matching is… |
| `delete-primer` | W | Remove a primer by `{name}` (or `{id}`) or exact `{sequence}`, optionally scoped to a `{collection}`. Body: one of `{name}`/`{id}`/`{sequence}` (+… |
| `delete-primer-collection` | W | Delete a named primer collection + its primers. Body: `{name}`. Refuses to delete the last remaining named collection (matches the GUI… |
| `delete-primers` | W | Batch-delete primers by name in ONE locked save — the machine-friendly path for a large cleanup (agent-API feedback SC-L), so a 600-item prune… |
| `delete-protein-motif` | W | Delete a user-stored protein motif override. Built-in motifs cannot be deleted — only user overrides. If `name` matches a built-in but the user… |
| `delete-protocol` | W | Delete a saved protocol — the FIRST match for `name` (case-insensitive), optionally scoped to `collection`. Body: `{name, collection?}`. |
| `delete-protocol-collection` | W | Delete a protocol collection AND every protocol in it. Body: `{name}`. |
| `design-gb-part` | R | Design Golden Braid / MoClo domestication primers. Body: `{template, start, end, part_type, grammar?, target_tm?, codon_taxid?,… |
| `design-guides` | R | CRISPR guide design over the loaded plasmid or a supplied sequence. |
| `design-mutagenesis` | R | Design SOE-PCR primers for a single-site mutation. Body: `{cds_dna, mutation, codon_taxid?}`. `mutation` is a string like `"W140F"` (WT-aa,… |
| `design-primers` | R | Primer-pair design over a target region. Body: `{template, start, end, mode?: "detection"\|"cloning"\|"generic", target_tm?: float, site_5?: str,… |
| `design-rbs` | R | Reverse-design a 5'UTR (Shine-Dalgarno + spacer) for a target relative RBS strength. Body: `{cds, target}` — `cds` begins with the start codon,… |
| `design-synthesis-fragment` | R | Wrap a sequence in the correct nested overhangs so the directly- SYNTHESISED fragment (order it as a gBlock — no PCR) becomes a Level-0 part once… |
| `diff-plasmid` | R | Pairwise alignment of the loaded record against another plasmid in the library. Body: `{target_id, mode?, circular?}`. |
| `digest` | R | Digest a RAW sequence with restriction enzymes and report the cuts + resulting fragments with their overhangs — the overhang-aware QC a… |
| `discard-changes` | W | Discard unsaved canvas edits and revert to the library-stored copy — the agent recovery path for a stuck-dirty canvas (agent-API feedback P1-10).… |
| `domesticate-part` | W | Domesticate a part into its configured Golden-Braid / MoClo entry vector — the real Synthesis-tab L0 clone, now agent-reachable (agent-API… |
| `domesticate-parts` | W | Batch domestication — run `domesticate-part` for each item in ONE call (agent-API feedback P1-6), so building N L0 parts doesn't round-trip the… |
| `download-hmm-database` | W | Download + index an HMM database by `id` so hmmscan can use it — the headless equivalent of the GUI catalog modal's Download button. Body: `{id}`… |
| `download-plasmidsaurus` | W | Download a Plasmidsaurus item's `results` zip by item code and import every .gbk sample into the library. Body: `{item_code, [kind]}`. |
| `duplicate-experiment` | W | Copy a notebook entry under a fresh id — "start from the write-up I already have". Body: `{id: str, title?: str}`. |
| `edit-residue` | W | Change one residue of an annotated CDS ON THE PLASMID. |
| `experiment-backlinks` | R | Which notebook entries reference an object — the reverse of the `@plasmid` / `!action` / `&gel` body tags. Body: `{ref: str \| list[str], kind?:… |
| `experiment-step-template` | R | A heading-per-step skeleton for the named bench actions. Body: `{actions: list[str], heading?: str}`. |
| `export-commercialsaas` | W | Write the current record to `path` as CommercialSaaS `.dna`. Body: `{path}`. Requires the loaded record to be in the active library (`.dna` writer… |
| `export-embl` | W | Write the current record to `path` as EMBL flatfile. Body: `{path}`. EMBL is a near-equivalent of GenBank — same feature model, different text… |
| `export-experiment` | R | Render a notebook entry — or a whole project — as a document. Body: `{id?: str, project?: str, format?: "markdown"\|"html"}`. |
| `export-fasta` | W | Write the current sequence to `path` as FASTA. Body: `{path}`. Single-record FASTA with the plasmid name as the header. |
| `export-genbank` | W | Write the current record to `path` as GenBank. Body: `{path}`. Path may be relative (resolved against $HOME). Doesn't change the record's… |
| `export-gff` | W | Write the current record to `path` as GFF3. Body: `{path}`. Wrap features become two GFF3 rows joined by a shared `ID=...` attribute. Returns… |
| `export-map-image` | W | Render the LOADED plasmid's circular map to `path` as PNG or SVG. Body: `{path, format?, size?, transparent?, labels?, sites?}`. `format` ∈… |
| `export-migrate-archive` | W | Package ALL user data into one portable, compressed `.zip` (backup, or moving to another machine). Body: `{path, include_hmm?}`. `path` must end… |
| `export-primers` | W | Write primers to a CSV order sheet. Body: `{path, format?: "generic"\|"idt"\|"plate", collection?, ids?, names?, scale?, purification?}`. |
| `features` | R | List features on the loaded record (idx, label, type, start, end, strand). |
| `fetch` | W | Fetch a GenBank record from NCBI by accession and load it into the canvas. Body: `{accession, force?}`. |
| `find-feature` | R | Find features on the LOADED plasmid by NAME — forgiving of case, spacing, hyphens and Greek letters (`lacZ alpha` finds `LacZ-α`), best match… |
| `find-orfs` | R | Six-frame ORF scan over the loaded record. Body: `{min_aa?: int, include_alt_starts?: bool}`. Defaults: 30 aa, ATG only. The length cutoff is in… |
| `find-sequence` | R | Locate a DNA subsequence / IUPAC motif on the loaded record. Body: `{query, both_strands?, max_hits?}`. |
| `fold-rna` | R | Fold a sequence to its minimum-free-energy RNA secondary structure. Body: `{sequence}` (alias `{seq}`). DNA `T` is read as RNA `U`. Returns `{ok,… |
| `fpbase-search` | R | Search FPbase (fpbase.org) for FLUORESCENT PROTEINS by name. Body: `{query, max_hits?}` (`query` is a substring, e.g. "mCherry", "GFP"). |
| `genbank-search` | R | Search NCBI for SEQUENCE RECORDS by term. Body: `{query, db?, max_hits?}` — `db` is "nucleotide" (default) or "protein". |
| `get-active-codon-table` | R | Return the active codon-table taxid preference (empty string if unset → Synthesis falls back to the first registry entry, K12 by convention).… |
| `get-active-collection` | R | Return the active plasmid collection name (or `None` if none is active). Read-only counterpart to `set-active-collection`. |
| `get-active-enzyme-collection` | R | Return the currently active enzyme collection name (or `None` if no collection is active — the scanner then uses the full catalog). |
| `get-active-experiment-project` | R | Return the active experiment-project name (or `None`). Counterpart to `set-active-experiment-project`. |
| `get-active-hmm-database` | R | Return the active HMM database id (empty string if unset). Counterpart to `set-active-hmm-database`. |
| `get-active-parts-bin` | R | Return the active parts-bin name (or `None`). Counterpart to `set-active-parts-bin`. |
| `get-active-primer-collection` | R | Return the active primer collection name (or `None`). Counterpart to `set-active-primer-collection`. |
| `get-custom-enzyme` | R | Return a single custom enzyme by `name`. Body: `{name}`. |
| `get-custom-labware` | R | The Opentrons `definition` for a custom labware. Body: `{name}`. Reference it from a plan's labware entry as `{"labware": name, "slot": N,… |
| `get-entry-vector` | R | Return the entry vector for a grammar. Body: `{grammar_id}`. Includes the full `gb_text` in the response so the agent can parse + render it… |
| `get-enzyme-collection` | R | Return one enzyme collection by `name`. Body: `{name}`. |
| `get-experiment` | R | Fetch one notebook entry (full body + metadata) by id. Body: `{id: str}`. |
| `get-feature` | R | Detail for feature at `idx`. Returns label, type, start, end, strand, color, qualifiers (= the full SeqFeature qualifiers dict). |
| `get-feature-colors` | R | Return the user's custom feature-type → colour overrides. Body: `{}`. |
| `get-feature-library` | R | Fetch a feature-library entry by (name, feature_type). Body: `{name, feature_type?}`. If feature_type is omitted, matches the first entry by name.… |
| `get-feature-preset` | R | Fetch one built-in preset by name, including its bases. Body: `{name, feature_type?}`. Name match is exact but case-insensitive. 404 when there is… |
| `get-gel` | R | Fetch one saved gel snapshot (full lane payload + notes). |
| `get-grammar` | R | Return the full grammar dict for `grammar_id`. Body: `{grammar_id}`. |
| `get-history` | R | Return the construction-history tree for a library entry as nested JSON. Body: `{name}` (or `{id}`). Returns `history=null` when the entry has no… |
| `get-hmm-database` | R | Detail for one HMM database by `id`. |
| `get-part` | R | Fetch a single part by `name` (and optional `grammar` to disambiguate when two grammars carry a part with the same name). Returns the full entry… |
| `get-plasmid-protocol` | R | A plasmid's recorded construction history as a bench protocol. Body: `{id: str, collection?: str}`. |
| `get-primer` | R | Return a primer-library entry by `{name\|id}` or exact (case-insensitive) `{sequence}`. Optional `{collection}` scopes the lookup to that… |
| `get-protocol` | R | The saved `plan` for a protocol. Body: `{name}`. The plan is compile-able by `ot2-compile` / `ot2-run`. |
| `get-sequence` | R | Return DNA from `[start, end)`. Body: `{start, end, bottom?}`. `bottom: true` returns the reverse-complement (5'→3' on the bottom strand).… |
| `get-settings` | R | Return every allowlisted user-toggle setting with its current value, default, and validator hint. Useful for an agent that wants to inspect which… |
| `gibson-assemble` | W | Run a Gibson simulation AND save the product to the library. Body extends `simulate-gibson` with: |
| `golden-gate-assemble` | W | Run a Golden Gate / MoClo (Type IIS) assembly AND save the product to the library — the Type IIS counterpart to `gibson-assemble`. Body extends… |
| `guide-cloning-oligos` | R | The annealed oligo pair that clones a spacer into a Type IIS guide vector. Body: `{guide: str, overhang_top?, overhang_bottom?, add_g?: bool}`. |
| `guide-offtargets` | R | Mismatch-tolerant search for a spacer in a sequence. Body: `{guide: str, sequence?, circular?, variant?, max_mismatch?: int, require_pam?: bool}`;… |
| `hmmer-web` | R | Remote hmmscan vs Pfam via the EBI HMMER web API. Body: `{query, max_hits?}` — `query` is a PROTEIN sequence. |
| `hmmscan` | R | Scan a query AA sequence against a HMMER 3 profile DB. Body: `{query, hmm_path, max_hits?}`. `hmm_path` must point at a readable .hmm / .h3m /… |
| `import-feature-preset` | W | Copy one or more built-in presets into the user's feature library. |
| `index-library` | W | Index the user's OWN work - plasmids, primers, parts and lab-notebook entries - into a PRIVATE Babs pack, so Babs can answer questions about their… |
| `ingest-url` | W | Persist ONE web page into Babs' knowledge corpus so she knows it in every future session. Body: `{url, title?}`. |
| `learn-list` | R | List Background Learning sessions on this machine: `{slug, topic, kept, status}` each. Read-only. |
| `learn-merge` | W | Fold a finished Background Learning session's corpus into Babs' MAIN corpus, making what it learned permanent. Body: `{slug}`. |
| `learn-results` | R | What a Background Learning session has KEPT so far. Body: `{slug, limit?}`. Per-paper records (title, url, relevance_score, hop) sorted by score,… |
| `learn-start` | W | Launch a Background Learning session: grow Babs' corpus about a TOPIC via a focused, drift-resistant crawl of open-access databases. Body:… |
| `learn-status` | R | Progress of a Background Learning session. Body: `{slug}`. Returns the live counters (kept / fetched / dropped / frontier / status / drift),… |
| `lint-synthesis` | R | Synthesis / assembly-readiness pre-flight for a construct — the "is it safe to order + assemble?" check. Body: `{sequence \| id \| name,… |
| `list-backups` | R | List recoverable backups for user-data files. Body: `{label?}`. |
| `list-cas-variants` | R | The Cas nucleases `design-guides` supports. Returns `{variants: [{key, label, pam, pam_side, guide_len, cut_offsets, notes}, ...]}`. |
| `list-codon-tables` | R | Available codon usage tables. Returns name, taxid, source (builtin/kazusa/user) and date added for each entry, plus the `active_taxid` field… |
| `list-collections` | R | Plasmid collection buckets. Returns name + plasmid count for each collection, and which one is currently active. |
| `list-custom-enzymes` | R | List every user-added custom enzyme. Built-in NEB enzymes are NOT included — fetch the combined view via list-restriction-sites or… |
| `list-custom-labware` | R | Every custom labware definition as `[{name, collection}]`. |
| `list-entry-vectors` | R | Return every grammar/role that currently has an assigned entry vector. Each item carries `{grammar_id, role, name, size, source}` — `gb_text` is… |
| `list-enzyme-collections` | R | List every enzyme collection (named subset of the master catalog). Each item carries `{name, enzymes: [str, ...]}`. |
| `list-enzymes` | R | Look up the restriction-enzyme catalog: recognition site, cut offsets, overhang and isoschizomers, with no sequence involved. |
| `list-experiment-actions` | R | The curated bench-action catalog behind the `!<id>` body tags — `[{group, id, description}, …]`. Curated, not enforced: a body may carry any… |
| `list-experiment-projects` | R | List experiment projects + the currently active one. Returns project name, description, save-stamp, and entry count per project (no body content —… |
| `list-experiments` | R | List entries in an experiment project's notebook. Pass `{project: "X"}` to list only project X's own entries (SC-M); omit it for the active… |
| `list-feature-library` | R | Return the feature library. Each item: `{name, feature_type, strand, sequence_length, color}`. Sequences themselves can be long; agents needing… |
| `list-feature-presets` | R | List the built-in feature presets — the shipped, GenBank-verified catalogue of common plasmid elements (markers, origins, promoters, terminators,… |
| `list-features` | R | List features on the loaded record (idx, label, type, start, end, strand). |
| `list-gels` | R | List every saved gel snapshot. Returns id, name, lane count, agarose %, and timestamps per gel (omits per-lane detail — fetch via `get-gel`). |
| `list-grammars` | R | Return every grammar (built-in + custom). Each item carries `{id, name, enzyme, level_up_enzyme, editable, n_positions, catalog_size}`. Fetch the… |
| `list-hmm-databases` | R | Every registered HMM database. Each item: `{id, name, ready, builtin}` (`ready` = pressed/usable on disk). `active` is the picked id used by hmmscan. |
| `list-labware-collections` | R | Named custom-labware collections as `[{name, count}]`. |
| `list-library` | R | Plasmid library entries. Returns name, id, length (bp), n_features, topology, source for each. Subset of the on-disk JSON tuned for an agent's… |
| `list-parts` | R | List parts in a parts bin. Optional body: `{bin?, grammar?: str, level?: int, position?: str}`. Pass `{bin: "X"}` to list only bin X's OWN parts… |
| `list-parts-bins` | R | Every named parts bin. Each item: `{name, n_parts, description}`; `active` is the currently-selected bin name (empty if none). |
| `list-plasmid-statuses` | R | Return the canonical workflow-status vocabulary. Discoverability helper so agents don't have to hard-code `DESIGNING`/`CLONING`/… |
| `list-plasmidsaurus-items` | R | List your Plasmidsaurus items (sequencing orders), most-recent first, over the official OAuth2 REST API. Read-only; no body. |
| `list-plasmidsaurus-members` | R | List the GenBank-format members of a Plasmidsaurus result zip. Body: `{path}`. Returns `{members: [{name, size}, ...]}`. |
| `list-polymerases` | R | The polymerases `pcr-program` knows. Returns `{polymerases: [{key, label, denature_c, extend_c, s_per_kb, anneal_rule}, ...]}` — the… |
| `list-pre-update-snapshots` | R | List pre-update snapshots created by `splicecraft update`. Returns each snapshot's id, from_version, mtime, and counts — same data the… |
| `list-primer-collections` | R | Every named primer collection. Each item: `{name, n_primers}`. The unnamed "default" collection (= top-level primers.json) is surfaced as `{name:… |
| `list-primers` | R | Return the primer library. Each item: `{name, sequence, source, tm, status, type, date, notes}`. Sequence + notes can be long; clients on small… |
| `list-protein-motifs` | R | List the merged protein-motif library (built-ins + user overrides). Each entry carries name, feature_type, sequence, color, description. Sourced… |
| `list-protocol-collections` | R | Named protocol collections as `[{name, count}]`. |
| `list-protocols` | R | Every saved OT-2 protocol as `[{name, collection}]`. Load a plan with `get-protocol`; save one with `save-protocol`. |
| `list-restriction-sites` | R | Scan the loaded record for restriction sites. Body: `{enzymes?: [str, ...], min_length?: int, unique_only?: bool, respect_active_collection?: bool… |
| `literature-search` | R | Search the scientific literature via Europe PMC (PubMed + PMC + preprints). Body: `{query, max_hits?}`. |
| `load-entry` | W | Load a plasmid library entry by name or id. |
| `load-file` | W | Load a `.gb` / `.gbk` / `.dna` file from a server-side path. Body: `{path, name?, id?, force?}`. Bypasses the 1 MiB JSON-body cap by reading the… |
| `make-l0-part-from-fragment` | W | Turn a saved synthetic fragment into a Level-0 PART **and the L0 plasmid it lives on**, by cloning it into the grammar's entry vector. Body:… |
| `map-transcription` | R | Walk the circle from every promoter to every CDS — "what ELSE transcribes this gene?" |
| `memory-forget` | W | Delete one memory by slug. Body: `{slug}`. Get slugs from memory-list. |
| `memory-list` | R | Everything in Babs' persistent memory: `[{slug, title}]`. Use it to check what she already knows before saving a duplicate, and to get the slug… |
| `move-experiment` | W | Move a notebook entry between projects without a create/delete round-trip (SC-M, the experiment analog of `move-primer`/`move-part`). Body: `{id,… |
| `move-part` | W | Move a part between bins without a create/delete round-trip (SC-M, the parts analog of `move-primer`). Body: `{name\|id, to, from?, grammar?}`. |
| `move-plasmid` | W | Move a library plasmid between collections without a copy/delete round-trip (the plasmid analog of `move-primer`/`move-part`/ `move-experiment`;… |
| `move-plasmids` | W | Batch-move library plasmids into a collection in ONE locked save — the bulk analog of `move-plasmid` (agent-API feedback P1-6). Body: `{names:… |
| `move-primer` | W | Move a primer between collections WITHOUT a create/delete round-trip (SC-K — that round-trip lost primers: `create-primer` rejected the duplicate,… |
| `multi-align` | R | Pairwise-align a query against MULTIPLE targets at once — the batch form of `diff-plasmid` (the Alt+A multi-align overlay). Body: `{targets:… |
| `new-plasmid` | W | Create a new plasmid from a raw sequence and load it onto the canvas — the agent counterpart to Ctrl+N (New Plasmid). Body: `{name, sequence,… |
| `optimize-protein` | R | Codon-optimize a 1-letter AA sequence to DNA using a codon table from the registry. Body: `{protein, table?, stops?, transl_table?, mode?,… |
| `ot2-analyze` | R | Compile (if given a plan) and analyse a protocol on the robot itself — the robot's built-in simulate. Catches deck / labware / tip errors a run… |
| `ot2-calibration` | R | OT-2 calibration status — deck calibration + per-pipette offset calibration + tip-length calibrations. Body: `{host}`. Returns the `calibration`… |
| `ot2-compile` | R | Compile a plate-transfer plan into an Opentrons OT-2 protocol (offline — no robot needed). Validates the plan against the built-in deck catalog first. |
| `ot2-disengage` | W | De-energise the OT-2 gantry motors so the carriage can be moved by hand and no motor holds torque — a zero-descent power-down (homes nothing; no… |
| `ot2-home` | W | Home the OT-2 gantry (all axes to their reference). Safe motion — no plunger, no labware descent. Body: `{host}`. A `write` endpoint (bearer… |
| `ot2-lights` | W | Toggle the OT-2 rail (status) lights — a zero-motion control, handy to blink / identify the robot, confirm you're talking to the right one, or… |
| `ot2-normalize` | R | Compute per-sample volumes to NORMALISE concentration across a set of wells — the everyday plate-prep task (equalise DNA mass, or hit a target… |
| `ot2-plate-map` | R | Map a plasmid collection onto a labware's wells (row-major) — the "a plate IS a collection" link. Body: `{"collection": "MyPlasmids", "labware":… |
| `ot2-position-check` | W | Run a MOTION-ONLY position check — the gantry tours the deck (moves to the top of each loaded labware's wells) so you can verify alignment +… |
| `ot2-run` | W | Run a plate-transfer protocol on the OT-2 — GATED PHYSICAL ACTUATION. |
| `ot2-run-control` | W | Pause / resume / stop the OT-2's current run — MANUAL run control (the counterpart to `ot2-run`'s automatic stop-on-fault). Lets you hold a run to… |
| `ot2-status` | R | Full OT-2 state snapshot — the monitoring surface. Reports reachability + versions, pipette OK flags + volume specs, motor engagement, deck /… |
| `patent-search` | R | Search patents. Body: `{query, max_hits?}`. |
| `pcr-program` | R | A thermocycler program for one PCR. Body: `{product_bp, primer_tms?: [float, ...], primers?: [name\|id, ...], polymerase?: "q5"\|"phusion"\|"taq",… |
| `plan-residue-edit` | R | What changing one residue of an annotated CDS would do — WITHOUT doing it. |
| `plasmid-overview` | R | One call for a first look at the LOADED plasmid: name, length, topology, unsaved state, its features, and every enzyme that cuts it exactly ONCE… |
| `plasmidsaurus-items` | R | List your Plasmidsaurus items (sequencing orders), most-recent first, over the official OAuth2 REST API. Read-only; no body. |
| `predict-transcript` | R | Reconstruct the MATURE mRNA of an annotated transcription unit on the loaded record, and score translation initiation on it. |
| `rbs-strength` | R | Relative E. coli translation-initiation strength of a ribosome binding site. Body: `{mrna, start}` — `start` is the 0-based index of the start… |
| `read-consensus` | R | Combine every sequencing read STORED on a plasmid into one answer. |
| `read-url` | R | Open ONE web page and return its readable text (HTML → text; NO JavaScript). Body: `{url, max_chars?}`. |
| `recall-knowledge` | R | Recall passages from Babs' knowledge corpus (papers, protocols and curated reference data) for a question. Body: `{query, max_passages?, packs?,… |
| `recover-history-from-dna` | W | Restore construction history to plasmids whose saved `.dna` original documents more of their build than the library currently holds. Body:… |
| `redo` | W | Redo the last undone canvas mutation — the agent counterpart to Ctrl+Y. Body: `{}`. Returns `{ok, redone, redo_remaining, undo_available,… |
| `reference-lookup` | R | Look a term up in Babs' CURATED reference registries - antibiotics and selection agents, promoters, terminators, UTRs/enhancers, plasmid origins,… |
| `remember-fact` | W | Save a durable fact to Babs' persistent memory. Body: `{fact, title?}`. |
| `rename-collection` | W | Rename a collection. Body: `{old, new}`. |
| `rename-experiment-project` | W | Rename a project. Body: `{name: str, new_name: str}`. If the renamed project is the active one, the active-pointer follows so subsequent… |
| `rename-parts-bin` | W | Rename a parts bin. Body: `{name, new_name}`. Matching is case-insensitive; 404 if `name` is unknown, 409 if `new_name` collides with another bin.… |
| `rename-plasmid` | W | Rename a library plasmid's display name. Body: `{old\|id, new}`. |
| `rename-primer-collection` | W | Rename a primer collection. Body: `{name, new_name}` (both named collections; the reserved default is the empty string and can't be renamed).… |
| `replace-sequence` | W | Replace `[start, end)` with `bases` (mutagenesis / edit). Body: `{start, end, bases}`. Pass `{"force": true}` to bypass the dirty-state guard.… |
| `restart` | W | Restart the `--agent` daemon so it picks up an on-disk update (SC-A). Body: `{}`. |
| `restore-backup` | W | Restore one user-data file from a chosen backup. Body: `{label, source_path}` where `source_path` must come from a prior `list-backups` response… |
| `restore-pre-update-snapshot` | W | Restore a pre-update snapshot by id (or "latest"). Body: `{id?: str}`. When `id` is missing or "latest", restores the newest snapshot. |
| `save` | W | Save the current record to its source file (if any) and library. |
| `save-custom-labware` | W | Save (or overwrite) a custom labware. Body: `{name, definition, collection?}` — `definition` is an Opentrons labware-definition object with a… |
| `save-protocol` | W | Save (or overwrite) a protocol. Body: `{name, plan, collection?}`. `plan` is a compile-able plan (see `ot2-compile`) — it is validated before… |
| `scan-promoters` | R | Scan for sigma-70 (-35 / spacer / -10) promoters. |
| `scan-terminators` | R | Scan for intrinsic (Rho-independent) terminators: a stem-loop followed by a U-tract. |
| `score-guide` | R | Triage a spacer you already have. Body: `{guide: str, u6_driven?: bool}`. |
| `scrub-plasmid` | R | Plan a clone-free restriction-site scrub. Body: `{seq?, features?, enzymes?, overlap?, method?, circular?, codon_taxid?}`. |
| `search-experiments` | R | Full-text search across notebook entries. Body: `{query?: str, tags?: list[str], project?: str, all_projects?: bool, limit?: int}`. |
| `search-library` | R | Cross-collection plasmid name search. Body: `{query?: str, limit?: int}`. Empty `query` returns the first `limit` plasmids across all collections… |
| `set-active-codon-table` | W | Persist the active codon-table preference. Body: `{taxid: str}` — pass empty string to clear (Synthesis then falls back to the first registry… |
| `set-active-collection` | W | Switch the active collection. Body: `{name}`. Mirrors the panel's Click-to-enter flow (writes the active pointer to settings.json + refreshes the… |
| `set-active-enzyme-collection` | W | Set or clear the active enzyme collection pointer. Body: `{name: str \| null}` — `null` clears the pointer. Refuses names that don't resolve to an… |
| `set-active-experiment-project` | W | Switch the active experiment project. Body: `{name: str}`. Mirrors the GUI flow (invariant #43 + sweep #9): updates the active-pointer, forces a… |
| `set-active-hmm-database` | W | Pick the active HMM database. Body: `{id}` (one of `list-hmm-databases`). Persists the `hmm_db_active_id` setting. |
| `set-active-parts-bin` | W | Switch the active parts bin. Body: `{name}` (one of `list-parts-bins`). The named bin in `parts_bin_collections.json` is the source of truth; its… |
| `set-active-primer-collection` | W | Switch the active primer collection AND re-point the live `primers.json` mirror at it. Body: `{name}` where `name` is one of the collections from… |
| `set-entry-vector` | W | Set or clear a grammar's entry vector. Body to set: `{grammar_id, name, gb_text, source?}`; body to clear: `{grammar_id, clear: true}`. `gb_text`… |
| `set-feature-color` | W | Set or clear the render colour for a feature type. Body: `{feature_type, color}` — `color` is a hex string (`#1f77b4` or the short `#abc` form).… |
| `set-plasmid-status` | W | Set / clear a library entry's workflow status. Body: `{name\|id, status}`. `status` must be one of the values from `list-plasmid-statuses`, or… |
| `set-protein-motif` | W | Create or override a protein motif. Body: `{name: str, sequence: str, feature_type?: str = "Motif", color?: str = "", description?: str = ""}`. |
| `set-setting` | W | Persist a single user-toggle setting. Body: `{key, value}`. `key` must be in the allowlist (see `get-settings`); `value` is type-checked +… |
| `shutdown` | W | Stop the `--agent` daemon gracefully. Body: `{}`, or `{"force": true}` to quit with unsaved canvas edits. |
| `simulate-gel` | R | Render an in-silico agarose gel. Body: `{lanes: [{name?, source, detail?}, ...], agarose_pct?: float = 1.0, template_seq?: str = "",… |
| `simulate-gibson` | R | Dry-run a Gibson assembly without saving. Body: `{fragments: [{name, sequence, features?}, ...], min_overlap?: int = 15, circular?: bool = true}`. |
| `simulate-golden-gate` | R | Dry-run a Golden Gate / MoClo (Type IIS) one-pot assembly. Body: `{parts: [seq, …], vector, enzyme?="BsaI"}`. |
| `simulate-pcr` | R | Run an in-silico PCR. Body: `{template_seq, fwd_primer, rev_primer, circular?: bool = true, max_amplicon?: int = 20000}`. Natural-guess aliases… |
| `simulate-traditional-cloning` | R | Dry-run a traditional restriction-cloning reaction: digest a vector and an insert with the chosen enzymes, then ligate. Body: `{vector_seq,… |
| `span-contains` | R | Circular-aware containment: does a span contain a base or another span, on a molecule where either may cross the origin? |
| `status` | R | Current session state: loaded record, dirty flag, source path. |
| `tools` | R | Self-describe: every endpoint with its read/write mode and full request-body documentation. |
| `traditional-clone` | W | Run a traditional restriction-cloning reaction AND save the product to the library — the restriction-cloning analog of `gibson-assemble`. Body… |
| `transfer-annotations` | W | Transfer features from a source plasmid (`source_id` in the library) onto the loaded record by sequence-identity matching. Body: `{source_id: str,… |
| `undo` | W | Undo the last canvas mutation (sequence edit, feature add/update/ delete, replace-sequence, new-plasmid load, …) — the agent counterpart to… |
| `uniprot-search` | R | Search UniProtKB for PROTEINS. Body: `{query, max_hits?}` (`query` is a UniProt query, e.g. "Cas9", "NPTII", "gene:BAR AND organism_id:3702"). |
| `update-custom-enzyme` | W | Replace an existing custom enzyme by `name`. Built-in NEB enzymes are not editable via this endpoint (refuse with 400) — add a custom enzyme with… |
| `update-enzyme-collection` | W | Replace an existing enzyme collection by `name`. Body: `{name, enzymes?: [str, ...], new_name?: str}`. `new_name` renames the collection (and… |
| `update-experiment` | W | Update an existing notebook entry. Body: `{id: str, title?: str, body_md?: str, tags?: list[str]}`. Fields not supplied keep their prior value.… |
| `update-feature` | W | Update feature at `idx` — full parity with the edit dialog. Body: `{idx, label? (or name), type?, strand?, color?, qualifiers?}`. Only the… |
| `update-feature-library` | W | Update a feature-library entry, looked up by (name, feature_type). Body: every field in `_agent_feature_dict`'s schema; missing fields preserve… |
| `update-gel` | W | Replace a saved gel's fields by id. Body: `{id: str, name?: str, lanes?: list[dict], agarose_pct?: float, notes?: str}`. Fields not supplied keep… |
| `update-grammar` | W | Replace an existing custom grammar by id. Built-in ids are refused (mirrors the UI's `editable=False` gate). |
| `update-part` | W | Update an existing part. Body: `{name, grammar?, ...new fields}`. Lookup is by `(name, grammar)`; if grammar is omitted matches the first part by… |
| `update-primer` | W | Edit a `primer_bind` feature's label / sequence / strand / notes. Body: `{idx, label?, primer_seq?, strand?, notes?}`. `idx` is the 0-based index… |
| `verify-against-reads` | R | Verify a designed construct against sequencing reads — the "I built X, I got reads back, do they match?" check. Body: `{reference \| reference_id \|… |
| `web-search` | R | General web search. Body: `{query, max_hits?}`. |
| `wikipedia-search` | R | Search Wikipedia for encyclopedic background. Body: `{query, max_hits?}`. |
