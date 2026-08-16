# SpliceCraft

![SpliceCraft Logo](https://raw.githubusercontent.com/Binomica-Labs/SpliceCraft/master/splicecraftLogo.png)

**A plasmid workbench you live in.** SpliceCraft is a terminal-native
viewer, sequence editor, primer + mutagenesis designer, Golden Braid /
MoClo cloning workbench, and in-process BLAST / HMMscan engine — all
rendered as crisp Unicode braille graphics in any modern terminal.

Fetch from NCBI, load `.gb` / `.gbk` files or `.dna` files from the
popular commercial plasmid editor file format (single or in bulk),
organize plasmids into named collections, design diagnostic / cloning /
Golden Braid primers via Primer3, run SOE-PCR site-directed mutagenesis
on any CDS, and search your own plasmid library by sequence similarity
— without ever leaving the shell.

**Built for daily lab work.** SpliceCraft is actively maintained by a
practicing bioengineer who uses it as their primary day-to-day tool
for plasmid design, cloning planning, and sequence triage. Bug reports
come from the bench; releases ship from the bench. Every feature has a
real-world job.

![SpliceCraft screenshot](https://raw.githubusercontent.com/Binomica-Labs/SpliceCraft/master/splicecraftScreenshot.png)

## Where to go next

- [**Install**](install.md) — pipx, pip, uv, conda, dev checkout.
- [**Getting started**](getting-started.md) — first commands, key
  bindings primer, opening pUC19 in five seconds.
- [**Features overview**](features.md) — what the workbench does
  without leaving the terminal: viewing, editing, cloning, primer
  design, mutagenesis, simulation, search, library, agent API.
- [**Keybindings and menus**](keybindings.md) — the full reference.
- [**Data safety and backups**](data-safety.md) — four-layer
  per-file backup, daily snapshots, pre-update snapshots, lock-file
  hardening.
- [**Agent API**](agent-api.md) — ~135 HTTP endpoints for external AI
  agents to drive the running session.
- [**CLI sidecar**](cli.md) — `splicecraft-cli`, the stdlib-only
  client (~50 ms cold start).
- [**Architecture**](architecture.md) — the hub + layered-siblings
  layout and how to navigate it.
- [**Citing SpliceCraft**](citation.md) — the Zenodo DOI, and
  `splicecraft --citation` for an APA reference plus BibTeX pinned to
  your version.

## At a glance

- **Hub + layered siblings.** A ~99k-line hub (`splicecraft.py`) plus
  ~25 flat `splicecraft_*.py` modules in clean import layers (L0→L7) —
  still greppable, with a no-cycle guard enforcing the layering.
- **5,000+ tests** anchored on 10 sacred invariants (see [CLAUDE.md](
  https://github.com/Binomica-Labs/SpliceCraft/blob/master/CLAUDE.md)
  for the registry). Property-based fuzzing on biology primitives.
- **No network in the hot path.** NCBI fetch + Kazusa codon-table
  fetch are user-initiated; everything else is local.
- **Atomic writes everywhere.** Every save goes through the
  four-layer backup system; pre-update snapshots are taken before
  any pip/pipx/uv subprocess.
- **No external BLAST / HMMER install.** `pyhmmer` ships HMMER 3
  compiled into the wheel; BLASTN, BLASTP, and HMMscan all run
  in-process.

See [Architecture](architecture.md) for the longer story.
