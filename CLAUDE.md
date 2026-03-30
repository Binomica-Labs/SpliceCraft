# SpliceCraft — Claude Development Memo

## What is SpliceCraft?

A **terminal-based circular plasmid map viewer and sequence editor** built with Python/Textual/Biopython. Single-file app (`splicecraft.py`, ~3440 lines). Renders Unicode braille-dot circular and linear plasmid maps directly in the terminal.

**Repo:** `github.com/Binomica-Labs/SpliceCraft` (Binomica Labs org, user ATinyGreenCell)

## Architecture (single file: `splicecraft.py`)

### Top-level structure (line numbers approximate):
- **Lines 1–95:** Imports, dependency check, library persistence (`plasmid_library.json`)
- **Lines 96–375:** NEB restriction enzyme catalog (~200 enzymes), IUPAC regex, cached pattern compilation, IUPAC-aware reverse-complement helper
- **Lines 386–480:** `_scan_restriction_sites()` — scans both strands, returns resite + recut dicts; palindrome-aware (skips reverse scan for palindromic sites)
- **Lines 483–800:** Sequence panel rendering — `_assign_chunk_features()`, `_render_feature_row_pair()`, `_build_seq_text()` — forward-strand features above DNA, reverse below, braille bars with arrowheads
- **Lines 800–950:** Codon table, clipboard, CDS translation, GenBank I/O (fetch from NCBI, load local .gb)
- **Lines 950–1070:** `_Canvas` (2D char grid) and `_BrailleCanvas` (sub-character resolution via Unicode braille U+2800–U+28FF)
- **Lines 1070–1740:** `PlasmidMap` widget — circular + linear map rendering, feature arcs, restriction site overlays, label placement algorithm, tick marks
- **Lines 1740–1820:** `FeatureSidebar` — DataTable of features with detail panel
- **Lines 1820–1920:** `LibraryPanel` — persistent plasmid collection (JSON), add/remove entries
- **Lines 1920–2260:** `SequencePanel` — DNA viewer with click-to-cursor, drag selection, double-stranded display, feature annotation bars
- **Lines 2260–2550:** Modal dialogs — `EditSeqDialog`, `FetchModal`, `OpenFileModal`, `DropdownScreen`
- **Lines 2550–2600:** `MenuBar` widget — File, Edit, Enzymes, Features, Primers, Genes
- **Lines 2600–2630:** `UnsavedQuitModal`
- **Lines 2630–3440:** `PlasmidApp` (main app) — keybindings, undo/redo stack, record loading, feature selection coordination between map/sidebar/sequence panel, menu actions, entry point

### Key design patterns:
- **All rendering uses Rich `Text` objects** — no curses
- **Braille canvas** gives sub-character pixel resolution (2x4 dots per terminal cell)
- **Feature coordination:** map click -> sidebar highlight -> sequence scroll (and vice versa via messages)
- **Undo/redo:** snapshot-based (stores full seq + cursor + SeqRecord), max 50
- **Restriction sites:** scanned on load/edit, stored as `resite` (recognition bar) + `recut` (cut marker) dicts
- **Caching:** PlasmidMap, SequencePanel, and regex patterns all cache rendered/compiled output keyed on state

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

### Bugs fixed (2026-03-30):
- **Palindromic RE double-counting** — `_scan_restriction_sites()` was scanning both strands for palindromic enzymes, creating 2 resites per physical site. The `unique_only` filter excluded all common palindromic enzymes (EcoRI, BamHI, HindIII, etc.). **Fix:** skip reverse scan for palindromes; add bottom-strand recut only.
- **Reverse-strand resite positions** — non-palindromic reverse-strand hits were placed at `n - p - site_len` instead of `p` (forward-strand coordinate). **Fix:** use `p` directly; map cut positions via `site_len - 1 - fwd_cut`.
- **`_rc()` IUPAC handling** — reverse-complement only translated ACGT, leaving ambiguity codes (R/Y/W/S/M/K etc.) unchanged. **Fix:** added full IUPAC complement table (`_IUPAC_COMP`).
- **Regex recompilation** — `_iupac_pattern()` recompiled ~200 patterns on every restriction scan. **Fix:** added `_PATTERN_CACHE` dict.
- **Duplicate enzyme entries** — SbfI and NspI each defined twice in `_NEB_ENZYMES`. **Fix:** removed duplicates.

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
