# SpliceCraft — Claude Development Memo

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer and sequence editor** built with Python/Textual/Biopython. Single-file app (`splicecraft.py`, ~3337 lines). Renders Unicode braille-dot circular and linear plasmid maps directly in the terminal.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs org, user ATinyGreenCell)

## Architecture (single file: `splicecraft.py`)

### Top-level structure (line numbers approximate):
- **Lines 1–95:** Imports, dependency check, library persistence (`plasmid_library.json`)
- **Lines 96–370:** NEB restriction enzyme catalog (~200 enzymes), IUPAC regex, reverse-complement helper
- **Lines 371–470:** `_scan_restriction_sites()` — scans both strands, returns resite + recut dicts
- **Lines 473–750:** Sequence panel rendering — `_assign_chunk_features()`, `_render_feature_row_pair()`, `_build_seq_text()` — forward-strand features above DNA, reverse below, braille bars with arrowheads
- **Lines 753–900:** Codon table, clipboard, CDS translation, GenBank I/O (fetch from NCBI, load local .gb)
- **Lines 901–1014:** `_Canvas` (2D char grid) and `_BrailleCanvas` (sub-character resolution via Unicode braille U+2800–U+28FF)
- **Lines 1024–1690:** `PlasmidMap` widget — circular + linear map rendering, feature arcs, restriction site overlays, label placement algorithm, tick marks
- **Lines 1692–1770:** `FeatureSidebar` — DataTable of features with detail panel
- **Lines 1770–1870:** `LibraryPanel` — persistent plasmid collection (JSON), add/remove entries
- **Lines 1872–2207:** `SequencePanel` — DNA viewer with click-to-cursor, drag selection, double-stranded display, feature annotation bars
- **Lines 2209–2500:** Modal dialogs — `EditSeqDialog`, `FetchModal`, `OpenFileModal`, `DropdownScreen`
- **Lines 2503–2542:** `MenuBar` widget — File, Edit, Enzymes, Features, Primers, Genes
- **Lines 2544–2572:** `UnsavedQuitModal`
- **Lines 2575–3337:** `PlasmidApp` (main app) — keybindings, undo/redo stack, record loading, feature selection coordination between map/sidebar/sequence panel, menu actions, entry point

### Key design patterns:
- **All rendering uses Rich `Text` objects** — no curses
- **Braille canvas** gives sub-character pixel resolution (2x4 dots per terminal cell)
- **Feature coordination:** map click -> sidebar highlight -> sequence scroll (and vice versa via messages)
- **Undo/redo:** snapshot-based (stores full seq + cursor + SeqRecord), max 50
- **Restriction sites:** scanned on load/edit, stored as `resite` (recognition bar) + `recut` (cut marker) dicts
- **Caching:** both PlasmidMap and SequencePanel cache rendered output keyed on state

## Current state (as of latest commit)

### Released features (v0.1.0, 2026-03-23):
- Braille circular map, NCBI fetch, local .gb loading, library, feature sidebar, sequence panel, undo/redo, restriction sites

### Unreleased features (in code, listed in CHANGELOG.md [Unreleased]):
- Feature deletion (Delete key)
- Linear map view toggle (v key)
- Strand-aware DNA layout (fwd above, rev below)
- Braille feature bars in sequence panel
- Single-bp feature triangles
- Label-above/label-below layout
- Feature connector lines (l key toggle)
- Full NEB enzyme catalog (~200 enzymes, Type IIS support)
- Inside tick marks on circular map
- Full-length feature labels (no 16-char truncation)
- Proximity label placement algorithm
- Default library entry (MW463917.1 / pACYC184)

### Menu items marked "coming soon":
- **Primers > Design Primer** — not implemented
- **Genes > Annotate from NCBI** — not implemented
- **Features > Add Feature** — stub only (`action_add_feature` just shows notification)

## How to run

```bash
cd ~/SpliceCraft
python3 splicecraft.py              # empty canvas
python3 splicecraft.py L09137       # fetch pUC19
python3 splicecraft.py myplasmid.gb # open local file
```

## Development notes

- **Single-file app** — all code in `splicecraft.py`, no package structure
- **No tests** — no test suite exists
- **Dependencies:** textual, biopython (installed system-wide via `--break-system-packages`)
- **WSL environment** — Ubuntu on WSL2, Python 3.12
- **Git auth:** gh CLI authenticated as ATinyGreenCell, push access to Binomica-Labs org via browser OAuth
