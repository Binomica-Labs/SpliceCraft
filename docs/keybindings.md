# Keybindings and menus

Press `?` in-app for the live reference (rendered as Markdown so you
can drag-select a combo to copy).

## Main screen

| Key            | Description                                  |
|----------------|----------------------------------------------|
| `[` / `]`      | Rotate map origin left / right (when map focused) |
| `← / →`        | Same as `[` / `]` (when map focused)         |
| `↑`            | Reset origin to 0 (when map focused)         |
| `Shift+[/]`    | Rotate coarse (10× step)                     |
| `,` / `.`      | Circular map aspect wider / taller           |
| `v`            | Toggle circular ↔ linear map                 |
| `l`            | Toggle feature label connector lines         |
| `r`            | Toggle restriction-site overlay              |
| `f`            | Fetch a record from NCBI by accession        |
| `Ctrl+O`       | Open a `.gb` / `.gbk` / `.dna` file from disk |
| `Ctrl+N`       | New Plasmid (paste sequence + optional annotate) |
| `Ctrl+B`       | BLAST modal (BLASTN / BLASTP / HMMscan)      |
| `Ctrl+Shift+A` | Add current plasmid to the library           |
| `Ctrl+A`       | Select-all sequence                          |
| `Ctrl+E`       | Enter sequence editor mode                   |
| `Ctrl+S`       | Save edits to file                           |
| `Ctrl+F`       | Add a new feature (from cursor or blank)     |
| `Ctrl+Shift+F` | Capture selection / feature → Feature library |
| `Ctrl+P`       | Primer Design workbench                      |
| `Enter`        | Highlight the feature enclosing the seq cursor |
| `Delete`       | Context-aware delete (feature or library entry) |
| `Ctrl+Z`       | Undo                                         |
| `Ctrl+Shift+Z` / `Ctrl+Y` | Redo                              |
| `Ctrl+C`       | Copy selection (top strand 5'→3', or AA when CDS highlighted) |
| `Alt+C`        | Copy selection (bottom strand, reverse-complement) |
| `F1` – `F4`    | Focus mode: library / map / features / sequence |
| `F5`           | Restore all panels (split-window layout)     |
| `F6` / `Ctrl+H` | Construction-history viewer (full-screen)   |
| `Alt+D`        | Capture UI snapshot to `<DATA_DIR>/ui_snapshots/` (bug-report attach) |
| `Alt+Shift+D`  | Toggle hover-status diagnostic row           |
| `?`            | Help modal                                   |
| `Ctrl+Q`       | Quit                                         |

## Mouse

| Action               | Description                                        |
|----------------------|----------------------------------------------------|
| Click DNA row        | Place cursor at that base                          |
| Click feature bar    | Highlight the feature, set cursor at its 5' end    |
| Click AA letter      | Highlight that codon's three bases on the strand   |
| Click restriction site | Highlight recognition span; tint upstream blue / downstream red per strand |
| Double-click         | Select full feature span                           |
| Drag                 | Select a sequence range                            |
| Scroll wheel         | Rotate map (when over map panel)                   |
| Click backbone       | Clear all panel highlights                         |

## Terminal-specific notes

The Ctrl/Shift/Alt key namespace is heavily intercepted by terminal
emulators. The current keymap was chosen to avoid common collisions
(see the 0.5.5.x churn in `CHANGELOG.md` for the history). If a
binding doesn't reach the app:

- `Alt+M` toggles **click-debug mode**: every keystroke + click is
  reported as a toast with the modifier set that actually arrived at
  the app. Use this to identify what the terminal swallowed.
- For Shift+click terminals where Shift is consumed by selection,
  Ctrl+click is registered as a Shift+click synonym.
- For Alt+combo terminals that send `Esc + key`, the app accepts
  both `Alt+X` and `Esc X`.

See `RELEASE_CHECKLIST.md` for the per-terminal smoke matrix the
maintainer runs before each release.

## Menus

| Menu        | Items                                                                            |
|-------------|----------------------------------------------------------------------------------|
| File        | Open · Fetch from NCBI · New Plasmid · Add to Library · Save · Export GenBank / GFF3 / FASTA · Align sequencing run (Plasmidsaurus) · Bulk import folder · Restore from backup · Quit |
| Settings    | Persisted toggles (RE overlay, primer binding length, custom enzyme list, …)     |
| Edit        | Edit Sequence · Undo · Redo · Add Feature · Capture → feat-lib · Delete Feature · Find plasmid… |
| Enzymes     | Show RE sites · Unique cutters · 6+/4+ bp sites · Connectors · Edit custom enzyme list… |
| Features    | Feature Library workbench                                                        |
| Primers     | Full-screen Primer Design workbench                                              |
| Mutato      | SOE-PCR site-directed mutagenesis designer (4-source CDS picker)                 |
| Parts       | Parts Bin (per-grammar; multi-bin via Parts Bin collections)                     |
| Constructor | Traditional cloning · Gibson assembly · Golden Braid / MoClo / custom grammar assembly |
| Simulator   | In-silico PCR (exact-match binding) + agarose gel rendering (0.5–4.0%, ladder / uncut / digest / amplicon lanes) |
| History     | Construction-history viewer (`<HistoryTree>` for the loaded plasmid)             |
| BLAST       | BLAST / HMMscan modal (Ctrl+B)                                                   |
