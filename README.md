# SpliceCraft

```
╔═══════════════════════════════════════════════════════════════════════════════╗
║ ⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶ ║
║                                                                               ║
║  ________       __________           _________             ____________       ║
║  __  ___/__________  /__(_)____________  ____/____________ ___  __/_  /_      ║
║  _____ \___  __ \_  /__  /_  ___/  _ \  /    __  ___/  __ `/_  /_ _  __/      ║
║  ____/ /__  /_/ /  / _  / / /__ /  __/ /___  _  /   / /_/ /_  __/ / /_        ║
║  /____/ _  .___//_/  /_/  \___/ \___/\____/  /_/    \__,_/ /_/    \__/        ║
║         /_/                                                                   ║
║                                                                               ║
║        ·  I n - T e r m i n a l   P l a s m i d   W o r k b e n c h  ·        ║
║                                                                               ║
║ ⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶⠶ ║
╚═══════════════════════════════════════════════════════════════════════════════╝
```

[![PyPI](https://img.shields.io/pypi/v/splicecraft.svg)](https://pypi.org/project/splicecraft/)
[![Python](https://img.shields.io/pypi/pyversions/splicecraft.svg)](https://pypi.org/project/splicecraft/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: Beta](https://img.shields.io/badge/status-beta-orange.svg)](https://github.com/Binomica-Labs/SpliceCraft)

![SpliceCraft screenshot](https://raw.githubusercontent.com/Binomica-Labs/SpliceCraft/master/screenshot.jpg)

> ⚠️ **Beta software.** SpliceCraft is under active development and the API,
> UI, and on-disk file formats may change between releases. Your data files
> (`plasmid_library.json`, `parts_bin.json`, `primers.json`) are auto-backed
> up to `*.bak` on every save, but please keep your own off-disk copies of
> anything critical. **Do not use SpliceCraft as your sole system of record
> for cloning work until it exits beta.**

A terminal-based circular plasmid map viewer, sequence editor, **primer design
workbench**, and **Golden Braid parts domesticator** — rendered entirely in your
shell. Fetch any GenBank record by accession, load local files, annotate features
with pLannotate, design diagnostic / cloning / Golden Braid primers with Primer3,
and edit sequences — without ever leaving the terminal.

## Quick start

```bash
pipx install splicecraft
splicecraft              # empty canvas
splicecraft L09137       # fetch pUC19 from NCBI
splicecraft myplasmid.gb # open a local GenBank (.gb/.gbk) or SnapGene (.dna) file
```

**That's it.** `pipx` creates an isolated virtual environment for SpliceCraft
and its dependencies, so it won't clash with system Python packages. If you
don't have pipx yet: `sudo apt install pipx` on Debian/Ubuntu/WSL2,
`brew install pipx` on macOS, or `python -m pip install --user pipx` elsewhere.
User data lives in the platform-appropriate data directory (see the
Installation section below).

---

## Features

### Core visualization & editing
- **Braille dot-matrix circular map** — plasmids rendered as crisp Unicode braille
  rings with per-strand feature arcs, directional arrowheads, and proximity-placed
  labels
- **Linear map view** — toggle with `v` for a horizontal strip layout
- **Dithered sequence panel** — per-base DNA viewer with feature bars, restriction
  site overlays, and double-stranded display
- **Live NCBI fetch** — pull any GenBank record by accession number on demand
- **Local file support** — open `.gb` / `.gbk` (GenBank) or `.dna` (SnapGene native) files directly from disk
- **Free rotation** — spin the origin left or right with `[` / `]`
- **Restriction enzyme overlay** — 200+ NEB enzymes including Type IIS
  (BsaI, BsmBI, BbsI, SapI) with visible recognition arcs + cut markers

### Libraries (all persist to JSON)
- **Plasmid library** — SnapGene-style collection, auto-saves on import, survives
  restarts, supports rename and handslip-protected delete
- **Parts Bin** — Golden Braid L0 parts catalog with user-domesticated parts
  including sequences and primer pairs
- **Primer library** — all designed primers with Tm, length, date, status
  (Designed / Ordered / Validated), multi-select for batch operations

### Primer design (Primer3)
- **Detection primers** — diagnostic PCR; Primer3 picks the ideal pair within
  a selected region, 450-550 bp product by default (configurable)
- **Cloning primers** — RE-site tails + GCGC padding; 30+ common enzymes or
  type a custom recognition sequence
- **Golden Braid primers** — Esp3I / BsmBI domestication for all L0
  positions (Promoter, 5' UTR, CDS, CDS-NS, C-tag, Terminator). Splitting
  L0 (Esp3I) from L1 (BsaI) lets domesticated parts survive the L1
  Golden Gate reaction without re-cutting.
- **Generic primers** — simple binding primers, no tails
- Primers can be added to the plasmid map as `primer_bind` features
- Scrollable `TextArea` for custom sequence input; highlighted text = target

### Annotation
- **pLannotate integration** — press `Shift+A` (or use the `◈` library button)
  to auto-annotate a plasmid against pLannotate's curated feature database.
  Optional — see install notes below

### Feature operations
- **Feature sidebar** — click a row to highlight on map; click the map to
  select the feature under the cursor
- **Undo/redo** — 50-deep snapshot stack for all sequence edits
- **Delete protection** — focus-aware Delete key; confirmation modal (default
  focus = No) for library entries
- **Clipboard** — OSC-52 copy works in Windows Terminal, iTerm2, modern WSL

### Data safety
- **Atomic saves** — all JSON files written via tempfile + `os.replace`
- **Automatic backups** — every save writes `*.json.bak` before overwriting
- **Corrupt-file recovery** — missing files don't crash; corrupt files auto-
  restore from `.bak` with a warning notification on startup

---

## Installation

Requires **Python 3.10+**.

### With pipx (recommended)

```bash
pipx install splicecraft
```

`pipx` installs SpliceCraft (and its Textual / Biopython / primer3-py /
platformdirs deps) into an isolated virtual environment and places the
`splicecraft` command on your `PATH`. This is the right approach on
modern Debian, Ubuntu, Fedora, and WSL2, where `pip install` into the
system Python is blocked by [PEP 668](https://peps.python.org/pep-0668/).

If you don't already have pipx:

```bash
sudo apt install pipx           # Debian / Ubuntu / WSL2
brew install pipx               # macOS
python -m pip install --user pipx  # everywhere else
pipx ensurepath                 # one-time; adds ~/.local/bin to PATH
```

### With pip inside a venv

If you prefer a plain pip workflow, use a virtual environment:

```bash
python3 -m venv ~/.venvs/splicecraft
~/.venvs/splicecraft/bin/pip install splicecraft
~/.venvs/splicecraft/bin/splicecraft
```

(Plain `pip install splicecraft` into system Python works on older
distros and inside conda envs, but will be rejected by PEP 668 on any
recent Debian-family system — use `pipx` or a venv instead.)

### User data location

User data (library, parts bin, primers) lives in the platform-appropriate
data directory:

  | Platform | Path |
  |---|---|
  | Linux   | `~/.local/share/splicecraft/` |
  | macOS   | `~/Library/Application Support/splicecraft/` |
  | Windows | `%APPDATA%\splicecraft\` |

Override with `SPLICECRAFT_DATA_DIR=/path/to/dir splicecraft`.

### From source

```bash
git clone https://github.com/Binomica-Labs/SpliceCraft.git
cd SpliceCraft
pip install -e .        # inside a venv, or pass --break-system-packages
```

### Optional dependencies

**pLannotate** (for automatic plasmid annotation via `Shift+A`) — requires conda:

```bash
conda create -n plannotate -c conda-forge -c bioconda plannotate
conda activate plannotate
plannotate setupdb          # one-time ~500 MB BLAST database download
# then run SpliceCraft from the same conda env
```

SpliceCraft runs fine without pLannotate — the annotation feature just
notifies the user how to install it if pressed.

---

## Usage

```bash
splicecraft              # empty canvas
splicecraft L09137       # fetch pUC19 from NCBI on launch
splicecraft myplasmid.gb # open a local GenBank (.gb/.gbk) or SnapGene (.dna) file
splicecraft --version    # print version
splicecraft --help       # quick usage hint
```

If running from a git clone (`pip install -e .`), the same commands work;
you can also still run `python3 splicecraft.py` directly.

---

## Key Bindings

### Main screen

| Key            | Description                            |
|----------------|----------------------------------------|
| `[` / `]`      | Rotate map origin left / right         |
| `Shift+[/]`    | Rotate coarse (10× step)               |
| `Home`         | Reset origin to 0                      |
| `,` / `.`      | Circular map aspect wider / taller     |
| `v`            | Toggle circular ↔ linear map           |
| `l`            | Toggle feature label connector lines   |
| `r`            | Toggle restriction-site overlay        |
| `f`            | Fetch a record from NCBI by accession  |
| `Ctrl+O`       | Open a `.gb` file from disk            |
| `Ctrl+Shift+A` | Add current plasmid to the library     |
| `Shift+A`      | Annotate plasmid with pLannotate       |
| `Ctrl+E`       | Enter sequence editor mode             |
| `Ctrl+S`       | Save edits to file                     |
| `Ctrl+F`       | Add a new feature (from cursor or blank) |
| `Ctrl+Shift+F` | Capture selection / feature → Feature library |
| `Delete`       | Context-aware delete (feature or library entry) |
| `Ctrl+Z`       | Undo                                   |
| `Ctrl+Shift+Z` | Redo                                   |
| `Ctrl+C`       | Copy selection to clipboard            |
| `q`            | Quit                                   |

### Primer Design screen

| Key   | Description                                    |
|-------|------------------------------------------------|
| `esc` | Close primer screen                            |
| `m`   | Mark / unmark primer under cursor (`★`)        |
| `M`   | Mark / unmark all primers                      |
| `S`   | Cycle status: Designed → Ordered → Validated   |
| `Tab` | Cycle focus between fields                     |

### Mouse

| Action         | Description                                 |
|----------------|---------------------------------------------|
| Click          | Place cursor / select feature under pointer |
| Double-click   | Select full feature span                    |
| Drag           | Select a sequence range                     |
| Scroll wheel   | Rotate map (when over map panel)            |

---

## Menus

| Menu       | Items                                                              |
|------------|--------------------------------------------------------------------|
| File       | Open .gb file · Fetch from NCBI · Add to Library · Save · Quit     |
| Edit       | Edit Sequence · Undo · Redo · Delete Feature                       |
| Enzymes    | Show RE sites · Unique cutters · 6+/4+ bp sites · Connectors       |
| Features   | Add Feature · Delete Feature · Annotate with pLannotate            |
| Primers    | Opens the full-screen Primer Design workbench                      |
| Parts      | Opens the Parts Bin (Golden Braid L0 parts catalog)                |
| Constructor| Opens the Assembly Constructor for TU building                     |

---

## Requirements

| Package           | Version           | Purpose                                      |
|-------------------|-------------------|----------------------------------------------|
| Python            | ≥ 3.10            | Runtime                                      |
| Textual           | ≥ 8.2.3           | TUI framework and rendering engine           |
| Rich              | ≥ 14.0            | Terminal rendering (Textual dependency)      |
| Biopython         | ≥ 1.87            | GenBank parsing and NCBI Entrez fetch        |
| primer3-py        | ≥ 2.3.0           | Primer design (Tm, thermodynamic screening)  |
| pytest            | ≥ 9.0             | Test suite (dev only)                        |
| pytest-asyncio    | ≥ 1.3             | Async test support (dev only)                |
| **pLannotate**    | optional, conda   | Automatic plasmid annotation (`Shift+A`)     |
| **BLAST+**        | optional, conda   | Required by pLannotate                       |
| **Primer3 CLI**   | optional, `apt`   | Not used directly — primer3-py bundles it    |

---

## Data files

All user data persists as human-readable JSON in the repo directory:

| File                  | Purpose                                     |
|-----------------------|---------------------------------------------|
| `plasmid_library.json`| Saved plasmid collection (GenBank + metadata) |
| `parts_bin.json`      | User-domesticated Golden Braid parts        |
| `primers.json`        | Designed primer library                     |
| `*.json.bak`          | Automatic backup — written before each save |

These files are in `.gitignore` — they're user-local data, not repo content.
A manual backup rotation happens on every save so accidental data loss is
always recoverable via the `.bak` file.

---

## Tests

```bash
python3 -m pytest -q          # full suite (483 tests, ~75 s)
python3 -m pytest tests/test_dna_sanity.py    # biology correctness only (< 1s)
```

---

## Codebase tour

SpliceCraft is a single-file Python app (`splicecraft.py`, ~10,360 lines) built on
Textual + Biopython. The single-file layout is intentional — no import puzzles,
and everything is greppable from one place.

### Layout at a glance

| Region | Contents |
|--------|----------|
| Top of file | imports, user data dir (`platformdirs`), rotating session-tagged logger, `_safe_save_json` / `_safe_load_json` (atomic writes + `.bak` recovery) |
| Biology core | NEB enzyme catalog, IUPAC-aware `_rc`, `_scan_restriction_sites` (palindrome- and wrap-aware), sequence panel renderer, `_translate_cds` |
| I/O | `fetch_genbank` (NCBI Entrez), `load_genbank` (`.gb` / `.gbk` / `.dna`), pLannotate subprocess integration |
| Rendering | `_Canvas` + `_BrailleCanvas` (sub-cell dot matrix), `PlasmidMap`, `SequencePanel`, `FeatureSidebar`, `LibraryPanel`, `MenuBar` |
| Workbenches | `PrimerDesignScreen` (detection / cloning / Golden Braid / generic), `MutagenizeModal` (SOE-PCR mutagenesis), `DomesticatorModal` + `ConstructorModal` (Golden Braid L0) |
| Registries | codon-usage tables, parts bin, primer library — all persisted as JSON with automatic `.bak` backups |
| Controller | `PlasmidApp` — keybindings, undo/redo (50-deep snapshot stack), `@work` threads for NCBI / pLannotate / Kazusa |

### Contributor docs

`CLAUDE.md` at the repo root is the **agent + contributor handover document**.
It covers the 10 sacred invariants (biology correctness rules with test coverage),
the error-handling convention, performance notes, and modular recipes for adding
modals, workers, menus, or persisted libraries without tripping the invariants.
Read it before touching the rendering layer, record pipeline, or primer design.

### Running the suite

```bash
python3 -m pytest -q                   # full suite (~75 s)
python3 -m pytest tests/test_dna_sanity.py    # biology only (< 1 s)
```

All tests run offline against synthetic `SeqRecord`s and monkeypatched data paths;
an autouse fixture in `tests/conftest.py` guarantees no test can write to real
user files.

---

## License

MIT
