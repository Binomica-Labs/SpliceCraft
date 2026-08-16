# SpliceCraft

![SpliceCraft Logo](https://raw.githubusercontent.com/Binomica-Labs/SpliceCraft/master/splicecraftLogo.png)

[![PyPI](https://img.shields.io/pypi/v/splicecraft.svg)](https://pypi.org/project/splicecraft/)
[![DOI](https://zenodo.org/badge/1190666059.svg)](https://zenodo.org/badge/latestdoi/1190666059)
[![100% Python](https://img.shields.io/badge/100%25-Python-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![TUI: Textual](https://img.shields.io/badge/TUI-Textual-5A45FF?logo=python&logoColor=white)](https://textual.textualize.io/)
[![Tests](https://github.com/Binomica-Labs/SpliceCraft/actions/workflows/test.yml/badge.svg)](https://github.com/Binomica-Labs/SpliceCraft/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: stable](https://img.shields.io/badge/status-stable-4EBF71.svg)](https://github.com/Binomica-Labs/SpliceCraft)

## Your whole cloning workflow, in the terminal.

SpliceCraft is a plasmid workbench that runs where you already work. Open a
map, edit the sequence, design primers, plan a Golden Braid or MoClo assembly,
BLAST a hit, check your Sanger reads, and keep a lab notebook — all from the
keyboard, in one place, no browser tab and no cloud account. Circular and
linear maps render as crisp Unicode braille graphics in any modern terminal,
and export to publication-quality PNG or SVG.

It's built by a practicing bioengineer for daily bench work: the bug reports
come from real cloning, and so do the fixes.

![SpliceCraft screenshot](https://raw.githubusercontent.com/Binomica-Labs/SpliceCraft/master/splicecraftScreenshot.png)

**Why give it a try:**

- **Fast and local.** No Electron, no web app, no login. `pipx install splicecraft` and you're designing in seconds.
- **It does the whole job.** View → edit → design → clone → simulate → verify → document — one tool that understands how those steps connect.
- **It guards your data like it's irreplaceable** (because it is — see below).
- **It's scriptable.** A 230+ endpoint local API and a stdlib CLI let an agent or a shell script drive every workflow.

## Quick start

```bash
pipx install splicecraft
splicecraft                      # empty canvas
splicecraft L09137               # fetch pUC19 from NCBI on launch
splicecraft myplasmid.gb         # local GenBank or .dna
```

Press `?` once running for the keyboard reference, or `Ctrl+K` for a fuzzy
command palette that jumps straight to any tool by name.

x86-64 Linux, Intel macOS, and Windows install entirely from prebuilt wheels.
On **ARM64 Linux** and **Apple Silicon**, one dependency (`primer3-py`) has no
ARM wheel and compiles at install, so install a C toolchain first
(`sudo apt install build-essential python3-dev`, or `xcode-select --install`).
On **Windows**, use Windows Terminal — the console is auto-configured for
UTF-8 + ANSI at startup, though this path is CI-tested via mocks and not yet
confirmed on real Windows hardware. If braille shows as boxes, toggle
**'ASCII plasmid map'** in Settings → Display.

Full matrix: [`docs/PLATFORMS.md`](docs/PLATFORMS.md) · other install methods:
[`docs/install.md`](docs/install.md).

## A workhorse that just works

Your plasmid library is months — sometimes years — of work, so SpliceCraft is
built to be a daily driver you never have to worry about:

- **Your data is sacred.** Every save is atomic, backed up (`.bak` + rotating
  timestamps + daily snapshots), and guarded by a "suspicious shrink" refusal
  that won't replace a 156 MB library with an empty file. Name collisions
  always ask — skip / copy / overwrite — and self-updates snapshot everything
  first. Bulk actions preview what they'll do, and one press of **`u`** puts
  the whole batch back.
- **The biology is correct, and proven.** Palindromes, Type IIS,
  origin-spanning cuts, wrap-around features, non-standard genetic codes, and
  IUPAC are pinned to the base — behind **4,000+ tests** plus property-based
  fuzzing on the biology, crash-injection on the save path, and concurrency
  fuzzing on the data layer.
- **We go looking for trouble.** A long list of sacred invariants
  ([`CLAUDE.md`](CLAUDE.md)) and multi-pass pre-release audits hunt edge cases,
  data-loss windows, races, and security gaps before they reach you.

Details: [`docs/data-safety.md`](docs/data-safety.md) ·
[`SECURITY.md`](SECURITY.md).

## A guided tour

Everything hangs off a menu bar across the top, read left to right. Full
reference: [`docs/features.md`](docs/features.md).

### BLAST

Search without leaving the app (`Ctrl+B`). **Local** runs BLASTN / BLASTP /
HMMscan against your own library in-process via `pyhmmer` — no external
`blast+` to install — with a one-click Pfam-A / NCBIfam downloader. **Online**
sends DNA, protein, a feature, or a whole plasmid to NCBI or EMBL-EBI and
tables the hits, with a Cancel that really stops. Agents can only search online
once you tick the setting for it, so a script can never silently ship your
sequences off-box.

### Enzymes

Drive the restriction overlay — all sites, unique cutters, 6+/4+ bp, or just
the Golden Braid connectors. Multi-cutters wear a live **superscript
cut-count** (EcoRI², BsaI³) that ticks down as you edit a site out. Build named
**enzyme collections** from the 200+ NEB catalog plus your own customs; the
active collection scopes every scan.

### Features

A library for your reusable annotations — promoters, RBSs, tags, CDSs. Capture
a region off any plasmid, then drop it onto another to annotate or splice it
in. **Ctrl+F** finds a subsequence fuzzily on both strands; **Ctrl+C** on a CDS
copies the protein rather than the DNA. **Alt+Shift+R** flips the whole record
end-for-end with every feature re-framed; **Alt+Shift+O** re-cuts a circular
plasmid so the cursor becomes base 1 — a real, undoable edit, refused on a
linear record. **File ▸ Find ORFs** runs a wrap-aware six-frame scan.

### Primers

A full-screen Primer3 designer for detection, cloning, Golden Braid, and
generic primers, each with a **Designed → Ordered → Validated** lifecycle. The
**Primer Check** tab runs in-silico PCR across your library: one primer lists
every plasmid it anneals to with % identity, strand, and position; two add the
amplicon length and the feature amplified, ranked ✓ / ⚠ / ~ / ✗. Binding is
judged on the 3′ end, so a 5′ cloning tail lowers identity rather than
vanishing. The library organises into collections with fuzzy search, and
exports order-ready **CSV** — generic or the **IDT bulk-upload** template.

### Mutato

Site-directed mutagenesis. Point at a CDS, name the change (`L54A`), and
SpliceCraft designs the SOE-PCR primers — falling back to a 2-primer
modified-outer strategy near the ends, and only offering the shortcut when the
primer genuinely carries the change, so you never amplify wild-type by
accident. It also turns a pasted protein into an orderable CDS with
frequency-matched codon optimization.

Its **Scrub** tab cures a whole plasmid of restriction sites with no cloning:
pick the enzymes and SpliceCraft finds the minimal point changes that kill each
site — silent across every overlapping reading frame, never spawning a new one,
and reported when a site can't be cured silently. **Apply cure** re-circularizes
by QuikChange or Golden Braid, really digesting and ligating each amplicon, so
**History** reads as a genuine assembly.

### Synthesis

A gene-synthesis composer in three tabs. **DNA** is a scrolling linear editor
with strand markers, feature stripes, live translation, feature-aware paste,
click-to-highlight restriction sites, and zoom-to-overview. **Protein** fills
codons in from your chosen table as you type, with a built-in motif library
(His6, FLAG, HA, TEV, P2A, NLS, +30) and a tabbed codon-table manager that
builds tables from an NCBI genome, a local CDS file, Kazusa, or TSV.
**Operon Design** reverse-designs every RBS in its real assembled context, and
domesticates a *natural* operon by curing forbidden Type IIS sites with
primer-encoded synonymous edits.

Compose a part, hit **Clone Fragment**, and pick a path: a modular grammar
hands it to the Domesticator; Gibson or Traditional open the Constructor.
Ordering synthetic instead? **L0 Fragment** wraps the sequence in the correct
nested overhangs for direct synthesis — two-tier aware, so a category pair
nests inside the entry vector's external pair automatically.

### Parts

Your **Parts Bin** — the Level-0 blocks for grammar-based assembly, in
per-grammar bins. Multiple bins live side by side as collections, so a yeast
toolkit and a plant toolkit never get mixed up.

### Constructor

The assembly bench: Traditional, Gibson, Golden Braid, MoClo, or your own
grammar, driven by a 4-source part picker. Every assembly lands as one library
entry carrying every parent feature forward, so you can trace a finished L3
construct back to its L0 parts. Deleted something by mistake? **`u`** brings it
back — the last 100 deletions of the session are undoable.

Ordered a synthetic fragment and it arrived? **New Part from Syn Frag** *runs
the cloning* — cutting fragment and entry vector, ligating, closing the circle
— then files the part and saves the finished L0 plasmid with its History. A
fragment whose ends don't match the type you picked is turned away rather than
filed as a part that can't assemble.

### Simulator

In-silico PCR and agarose gels. Pick a template, run the PCR, then save the
amplicon or send it to a gel lane. Gels render at 0.5–4% on a real
Helling–Goodman–Boyer mobility curve; stack lanes, save a gel to reload later,
or cite it as `&<gel>` in your notebook.

### Sequencing

Verify constructs against real reads. Drop in a Plasmidsaurus `.zip` or fetch a
run straight from their API, then walk run → sample → target and **Align**: the
read lands as a colored bar (blue match / red mismatch / gray gap) on the
linear map — click it to jump the sequence panel to that exact spot. **Bulk
auto-align** matches a whole folder in one pass. The **Verification Report**
grades every construct (✓ verified / ⚠ near / ~ partial / ✗ divergent), and a
true sub-100% identity never rounds up to "100%".

### Experiments

A genuine lab notebook in markdown: a split-pane editor, entries grouped into
**projects**, and live colored cross-references — type `@plasmid`, `!action`,
or `&gel` and `Ctrl+G` jumps to the source. Attach images, previewed inline,
and spellcheck with `F7`.

### History

Every plasmid remembers how it was made. **History** opens with a **Protocol** —
a numbered recipe that reads left → right like the bench — above a lineage tree
you can drill into as deep as you like. Pick any step and you get everything the
record holds: what it did and where, the conditions, the primers with each
binding site's position, strand and Tm, the fragment's end chemistry, and the
enzymes regenerated.

If a record claims something the plasmid doesn't bear out — a site that isn't
there, a primer that no longer binds — the view says so instead of presenting
it as fact. It only ever reports: your sequences and your history are never
edited to make a warning disappear. Lineage rides along through `.dna` import /
export, and `recover-history-from-dna` finds the original by sequence — even
under a different name — to put a lost build record back.

### BABS

A chat assistant that lives next to your plasmids: a direct conversation with a
**local [Ollama](https://ollama.com) model** — no cloud, no API key, nothing
leaves your machine. Streaming markdown answers, a **❤ context lifebar**, slash
commands, and a copy-pasteable transcript (**Ctrl+E** exports it).

Flip **Agent** on and Babs can *drive SpliceCraft herself*, calling the same
endpoints the `--agent` API exposes — read the plasmid, design primers, run a
digest, clone, manage the library, even drive the OT-2 — with everything
showing up live in the app. She **asks before every write** by default;
destructive whole-library wipes are never reachable, and **physical robot
motion always asks first**. Turn **Corpus** on and answers come grounded in a
research corpus with cited sources, grown by the built-in **paper scraper** and
topic-focused **Learn** crawls, or by indexing **your own library** into a pack
that is private by construction. Online database lookups (FPbase, UniProt,
Europe PMC, NCBI, patents) stay off until you tick the setting, and even then
only your query string is sent — never your sequence.

Needs Ollama running locally; `splicecraft babs-setup` bootstraps the engine in
one command. Details: [`docs/features.md`](docs/features.md).

### AUTOLAB

SpliceCraft drives an **Opentrons OT-2** liquid handler. **Find Robots** scans
your network for OT-2s; the **Deck** draws the robot top-down as a clickable
grid you fill with labware; the **Designer** builds an ordered step sequence
(transfer, distribute, consolidate, mix, delay, pause). A **Library** tab binds
a deck plate to a plasmid collection so wells map to your plasmids — then
cherry-pick, replate by identity, or normalise DNA to a target ng/µL. Compile
to a real Opentrons protocol, analyze on the robot's simulate, and run it
behind Arm + a clean analysis + health, calibration, pipette and door checks,
with live Pause / Resume / Abort and a telemetry panel (state, progress + ETA,
fault banner). Every bit of it is scriptable.

### File & Settings

**File** opens / fetches (NCBI) / saves / exports (GenBank · FASTA · GFF3 ·
map image as PNG/SVG, one plasmid or a whole collection), bulk-imports a
folder, and restores from backup; every GenBank it writes stamps a traceable
`Created by SpliceCraft v…` COMMENT. It also hosts the **selection → cloning
hub** (**Alt+Shift+P**): highlight any DNA and pick Traditional, Golden Braid /
MoClo, or Gibson — each opens pre-loaded with the selection *and its features*.
**Migrate Data** packages your entire setup into one checksum-verified `.zip`
to move between machines; **Master Delete** is a triple-gated full wipe.
**Settings** collects every toggle plus the grammar, entry-vector,
enzyme-collection, and codon-table editors.

### Scripting

A 230+ endpoint localhost JSON API (`splicecraft --agent`, or `--headless` for
a no-UI server with a `/healthz` probe) and a stdlib-only CLI
(`splicecraft-cli`, including a `call` passthrough to every endpoint) drive
every workflow. `/tools` self-describes each endpoint's request schema. See
[`docs/agent-api.md`](docs/agent-api.md) and [`docs/cli.md`](docs/cli.md).

## Documentation

| Topic                         | Where                                                                |
|-------------------------------|----------------------------------------------------------------------|
| Install methods               | [`docs/install.md`](docs/install.md)                                |
| First five seconds with pUC19 | [`docs/getting-started.md`](docs/getting-started.md)                |
| Full feature list             | [`docs/features.md`](docs/features.md)                              |
| Keybindings + menus           | [`docs/keybindings.md`](docs/keybindings.md)                        |
| Data safety + backups         | [`docs/data-safety.md`](docs/data-safety.md)                        |
| Agent API (HTTP)              | [`docs/agent-api.md`](docs/agent-api.md)                            |
| CLI sidecar                   | [`docs/cli.md`](docs/cli.md)                                        |
| Architecture                  | [`docs/architecture.md`](docs/architecture.md)                      |
| Sacred invariants             | [`CLAUDE.md`](CLAUDE.md)                                            |
| Contributing                  | [`CONTRIBUTING.md`](CONTRIBUTING.md)                                |
| Security policy               | [`SECURITY.md`](SECURITY.md)                                        |
| Changelog                     | [`CHANGELOG.md`](CHANGELOG.md)                                      |

## Tests

```bash
python3 -m pytest -n auto -q                  # full suite (~5–6 min on 8 cores)
python3 -m pytest tests/test_dna_sanity.py    # biology correctness only (< 2 s)
```

All tests run offline against synthetic `SeqRecord`s and monkeypatched data
paths; the autouse `_protect_user_data` fixture guarantees no test can write to
real user files.

## Maintenance

Actively maintained by a practicing bioengineer running real cloning workflows
in it daily; releases typically go out the same week a problem surfaces at the
bench. Issues and PRs welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md) before
opening a non-trivial one.

## Citing SpliceCraft

If SpliceCraft did part of the work in something you publish, please cite it.
Every release is archived on [Zenodo](https://zenodo.org/) with its own DOI, so
a methods section can point at the exact version that produced the design. The
concept DOI [`10.5281/zenodo.21960400`](https://doi.org/10.5281/zenodo.21960400)
always resolves to the newest release.

```bash
splicecraft --citation     # APA reference + BibTeX, pinned to your version
```

The repository also ships a [`CITATION.cff`](CITATION.cff), so GitHub's **Cite
this repository** button exports APA or BibTeX directly.

## License

MIT
