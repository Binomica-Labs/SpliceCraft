# Worked workflows

Every command here was run against a live SpliceCraft 1.3.0 session; the outputs
are real. Read the recipe matching the user's job — you don't need the others.

Set up once per shell if you're working against a sandbox rather than the user's
real library:

```bash
export SPLICECRAFT_DATA_DIR=/path/to/sandbox/splicecraft
```

---

## Annotate a bare NCBI record

Many GenBank records — including pUC19 as deposited — carry **zero** features.
The record is correct; nobody annotated it. `annotate-from-presets` matches the
built-in catalogue of ~117 curated vector elements plus the user's own feature
library by exact sequence.

```bash
splicecraft-cli fetch L09137                              # pUC19 → canvas (saved:false)
splicecraft-cli call annotate-from-presets --json '{}'    # DRY RUN — proposes only
```

```
7 proposals: AmpR promoter (2486–2591, −1), lac promoter (507–591, −1),
AmpR (bla) CDS (1625–2486, −1), lacZ-alpha CDS (145–469, −1),
M13 forward primer site (351–375, +1), M13 reverse (464–486, −1),
pUC19 MCS EcoRI-HindIII (395–452, +1)
```

Show the user the proposal, then apply and persist:

```bash
splicecraft-cli call annotate-from-presets --json '{"apply":true}'
splicecraft-cli call save --json '{"create":true}'     # create:true required — see below
```

Notes that matter:

- The response's `applied` field is the truth. `{}` returns `applied: false`.
- `already_annotated` counts matches the record already has, so re-running is a
  no-op rather than a way to stack duplicates.
- Matching is **exact**. An allelic variant won't be called — that's a refusal to
  guess, not a miss. Use `blast` for near-identical hits.
- `save` on a freshly-fetched record returns `409` unless you pass
  `{"create": true}`, because a `fetch` is inspect-only and silently filing it
  would pollute the active collection.
- The display name will be the GenBank LOCUS (`SYNPUC19CV`). Rename it to
  something a human recognises:
  `splicecraft-cli call rename-plasmid --json '{"old":"SYNPUC19CV","new":"pUC19"}'`

---

## Design primers and check they actually bind

Designing a primer is easy; proving it anneals where you think it does is the
part that matters. `amplify-feature` does both.

```bash
splicecraft-cli call amplify-feature --json '{"name":"lacZ-alpha"}'
```

```json
{"forward": {"seq":"CTATGCGGCATCAGAGCAGA","tm":60.0,"start":145,"end":165,
             "binds_as_designed":true,"other_sites":[]},
 "reverse": {"seq":"ATGACCATGATTACGCCAAGCT","tm":60.4,"start":447,"end":469,
             "binds_as_designed":true,"other_sites":[]},
 "product_size": 324,
 "primer_at_feature_start": "reverse",
 "strand_note": "this feature is on the − strand: its start (5' end) is at the
                 REVERSE primer, which begins with the feature's first bases"}
```

Always read `binds_as_designed` and `other_sites` before reporting a primer.
A non-empty `other_sites` means it anneals somewhere else too — a real
mispriming risk, not a formality.

Related: `check-primer` (one primer vs a template, both strands, wrap-aware),
`design-primers`, `pcr-program` (thermocycler conditions), `simulate-pcr`.
Persist with `create-primer` so the primer joins the user's ordering library.

---

## Traditional restriction cloning

**The simulation returns every compatible product and refuses to pick one.**
This is the single most important thing to understand about this endpoint.

```bash
splicecraft-cli call simulate-traditional-cloning --json '{
  "vector_seq": "<pUC19>", "vector_enzymes": ["EcoRI","HindIII"],
  "insert_seq": "<PCR product with flanking sites>", "insert_enzymes": ["EcoRI","HindIII"],
  "vector_circular": true, "insert_circular": false }'
```

Cutting pUC19 with EcoRI + HindIII yields two backbone fragments, so two
products come back:

| `vector_frag_idx` | orientation | length | what it is |
|---|---|---|---|
| 0 | reverse | 351 bp | the 51 bp MCS stuffer + insert — **not** what you want |
| 1 | forward | 2935 bp | the 2635 bp backbone + insert — the real construct |

Choose by **what the fragment carries** (origin, resistance marker), never by
size. Check `n_features`, or look for the selection marker. The same discipline
applies to `insert_circular: true`, where cutting a cassette out of a plasmid
leaves two candidate inserts, each with an `insert_frag_idx`, and neither is
picked for you.

Then read `self_ligation` on your chosen product — it's the off-target report:

```json
{"vector_self_closes": false, "empty_vector_bp": 2635,
 "vector_ends": {"left":"HindIII (5' AGCT)","right":"EcoRI (5' AATT)"},
 "double_insert": {"tandem": [], "inverted": []},
 "clean": true,
 "clean_note": "the vector's ends don't match each other (EcoRI (5' AATT) vs
                HindIII (5' AGCT)) and neither do the insert's"}
```

`vector_self_closes: true` means empty vector re-circularises — expect
background colonies, and tell the user (single-enzyme, blunt, or a
compatible-cohesive pair like SalI/XhoI all do this). `double_insert` flags
tandem/inverted multi-copy products. `orientation_ambiguous` means you cannot
tell the two orientations apart by digest.

To save a chosen product, pass the same payload to `traditional-clone` (the
write twin) and identify which product you want.

---

## Golden Gate / MoClo / Golden Braid

```bash
splicecraft-cli call simulate-golden-gate --json '{
  "parts": ["<part1>","<part2>"], "vector": "<destination>", "enzyme": "BsaI" }'
```

Parts may be in **any order** — the 4 nt overhangs set the order. Each part must
release exactly one fragment, so flank it with two inward-facing enzyme sites.

**Check `ok`, not the HTTP status.** A design that cannot close a circle is
served correctly (HTTP `200`) and reports the chemistry failing:

```
HTTP: 200   top-level ok: False
result.errors: ["part 1 released 0 fragment(s) when cut with BsaI
                (expected exactly 1 between two BsaI sites — check the part's
                 flanking sites / orientation)."]
```

Read `warnings` too: it flags a non-unique junction overhang, more than one circle
that could form from the same fragments, and residual enzyme sites (which would
be re-cut in the one-pot reaction). When nothing closes, the error lists every
fragment's two overhangs so the dangling one is visible — that's a design bug
you can point at, not a mystery.

A destination plasmid with background enzyme sites is released as several
pieces and the real reaction puts them all back; `n_vector_fragments` /
`n_vector_fragments_used` tell you what happened.

Pass `carry_annotations: true` to preview the product's feature map — the
preview and `golden-gate-assemble` (the write twin) take the same payload on
purpose, so what you previewed is what gets saved.

For the structured grammars, see `assemble-into-entry-vector` (multi-source
level-up into a configured acceptor), `design-gb-part`, `classify-part`, and
`scrub-plasmid` (domestication — removing internal enzyme sites synonymously).

---

## Codon-optimize a protein

```bash
splicecraft-cli call optimize-protein --json '{"protein":"MSKGEELFT..."}'
```

```json
{"dna":"ATGAGCAAAGGCGAAGAG…","cai":0.6912,"gc":50.26,"gc3":…,
 "table":"83333","n_codons":64,"warnings":[],"remaining_motifs":[],
 "remaining_repeats":[],"max_repeat":…}
```

Check `warnings`, `remaining_motifs` and `remaining_repeats` — an empty
`remaining_motifs` is what tells you the forbidden-site scrub actually
succeeded. Useful options (read `doc_full` for the full set): `taxid`/`table` to
pick the host, `mode`, `forbidden_motifs`, `min_gc`/`max_gc`, `avoid_repeats_with`.

A very high CAI target in an AT-rich host can push GC3 to an extreme that
synthesises badly — the response's GC3 warning is there for that reason; relay
it rather than reporting CAI alone.

The read counterpart is `analyse-cds`, which reports the codon usage of an
existing CDS.

---

## Look at what's on the canvas

```bash
splicecraft-cli status                       # name, length, topology, n_features, dirty
splicecraft-cli features                     # every feature with coordinates
splicecraft-cli get-sequence 100 200         # --bottom for the reverse complement
splicecraft-cli call find-orfs --json '{"min_aa":80}'
splicecraft-cli call list-restriction-sites --json '{}'
```

`list-restriction-sites` on pUC19 returns **309** hits across the full enzyme
catalogue. That's rarely what a user wants — filter to unique cutters, or scope
to an enzyme collection. Each hit carries `cut_bp`, `bottom_cut_bp`, `wraps`,
`aliases` (isoschizomers — `Esp3I` and `BsmBI` are the same enzyme) and
sometimes `methylation_context` (e.g. `Dam`), which matters when the plasmid
comes out of a standard *E. coli* strain.

---

## Verify a build against sequencing reads

```bash
splicecraft-cli call verify-against-reads --json '{...}'
splicecraft-cli call analyse-read-heterogeneity --json '{...}'
```

`analyse-read-heterogeneity` answers "is this culture one thing?" — per-base
allele fractions from the raw reads. A clone that matches the reference at the
consensus level can still be a mixed population; this is what catches it. For
Plasmidsaurus runs there's a dedicated cluster (`list-plasmidsaurus-*`,
`align-plasmidsaurus-zip`, `import-plasmidsaurus-*`).

---

## Record it in the lab notebook

SpliceCraft has a notebook keyed to the library, so the construct and the record
of building it live together.

```bash
splicecraft-cli call create-experiment --json '{
  "title":"pUC19-GFP build", "body_md":"Cut @pUC19 with EcoRI/HindIII…",
  "tags":["cloning"]}'
```

`@plasmid`, `!action` and `&gel` tokens in the body become cross-references.
`protocol-from-history` turns an actual construction history into protocol
steps — the provenance is already recorded, so the protocol doesn't have to be
written from memory.
