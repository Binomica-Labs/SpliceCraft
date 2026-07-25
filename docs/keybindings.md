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
| `Ctrl+O` / `o` | Open a `.gb` / `.gbk` / `.dna` file from disk |
| `Ctrl+N`       | New Plasmid (paste sequence + optional annotate) |
| `Ctrl+B`       | BLAST modal (BLASTN / BLASTP / HMMscan)      |
| `Ctrl+K`       | Command palette — fuzzy-jump to any tool     |
| `Alt+K`        | Add ("keep") the current plasmid in the library |
| `Ctrl+A`       | Select-all sequence                          |
| `Ctrl+E`       | Enter sequence editor mode                   |
| `Ctrl+S`       | Save edits to file                           |
| `Ctrl+F`       | Find a DNA subsequence (fuzzy, both strands) |
| `Alt+Shift+F`  | Add a new feature (from the current selection) |
| `Alt+Shift+C` | Capture selection / feature → Feature library |
| `Ctrl+P`       | Primer Design workbench                      |
| `Enter`        | Highlight the feature enclosing the seq cursor |
| `Delete`       | Context-aware delete (feature or library entry) |
| `Ctrl+Z`       | Undo — the sequence editor, or (when it has nothing left and the library has focus) the last plasmid delete |
| `Ctrl+Y`       | Redo (`Ctrl+Shift+Z` also works on terminals that report it) |
| `Ctrl+C`       | Copy selection (top strand 5'→3', or AA when CDS highlighted) |
| `Alt+C`        | Copy selection (bottom strand, reverse-complement) |
| `F1` – `F4`    | Focus mode: library / map / features / sequence |
| `F5`           | Restore all panels (split-window layout)     |
| `F6`           | Construction-history viewer (full-screen)    |
| `Alt+D`        | Capture UI snapshot to `<DATA_DIR>/ui_snapshots/` (bug-report attach) |
| `Alt+Shift+D`  | Toggle hover-status diagnostic row           |
| `?`            | Help modal                                   |
| `Ctrl+Q`       | Quit                                         |

## Library panel (when the plasmid list has focus)

| Key       | Description                                            |
|-----------|-------------------------------------------------------|
| `Space`   | Mark / unmark the highlighted plasmid                 |
| `c`       | Clear every mark                                      |
| `m` / `y` | Move / copy the marked plasmids to another collection |
| `p`       | Export the marked plasmids as circular-map images (PNG / SVG) |
| `s` / `h` | Set status / view history for the highlighted plasmid |
| `u`       | Undo the last plasmid delete of this session (restores it to its original slot) |

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

- **`Ctrl+Shift`+letter (and `Ctrl+I` / `Ctrl+H`) don't work on most
  terminals.** Without the Kitty keyboard protocol (VTE / Ptyxis, macOS
  Terminal, basic xterm), `Ctrl+Shift`+letter collapses to plain
  `Ctrl`+letter, and `Ctrl+I` / `Ctrl+H` arrive as Tab / Backspace. So
  every action uses a terminal-safe primary — `Alt+K` add-to-library,
  plain `c` clear-marks, `Alt+I` attach-image (Experiments), `Ctrl+Y`
  redo, `F6` history, `Alt+C` bottom-strand copy — with the old combos
  kept as aliases where the terminal reports them.
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

The top bar has 16 menus, left to right. All but **File** open their tool
directly (no dropdown). Press `Ctrl+K` for a fuzzy command palette that reaches
every one by name. Most menus also have an `Alt`+letter (shown in `?` Help).

| Menu        | Opens                                                                            |
|-------------|----------------------------------------------------------------------------------|
| File        | Open · Fetch from NCBI · New Plasmid · Add to Library (`Alt+K`) · Find plasmid · Diff · Find ORFs · Transfer annotations · Send selection to cloning · Save · Export (GenBank / FASTA / GFF3 / `.dna` / map image PNG · SVG) · Export collection · Collections · Migrate Data · What's New · Master Delete · Quit |
| Settings    | Persisted toggles (RE overlay, primer binding length, online lookups, ASCII map, …) |
| BLAST       | BLAST / HMMscan modal (`Ctrl+B`)                                                 |
| Enzymes     | Enzyme collections + settings (catalog, custom enzymes, unique / 6+ / 4+ cutters) |
| Features    | Feature Library workbench                                                        |
| Primers     | Full-screen Primer Design workbench (`Ctrl+P`)                                   |
| Mutato      | SOE-PCR site-directed mutagenesis + restriction-site Scrub                       |
| Synthesis   | Gene-synthesis composer (DNA · Protein · Operon design)                          |
| Parts       | Parts Bin (per-grammar; multi-bin via Parts Bin collections)                     |
| Constructor | Traditional · Gibson · Golden Braid · MoClo assembly                             |
| Simulator   | In-silico PCR + agarose gel rendering (0.5–4.0 %, ladder / uncut / digest / amplicon lanes) |
| Sequencing  | Plasmidsaurus run alignment / verification overlay                              |
| Experiments | Lab notebook (projects, entries, inline image attachments)                      |
| History     | Construction-history viewer (`F6` / `Alt+H`)                                     |
| AUTOLAB     | Opentrons OT-2 robot control (`Alt+U`)                                           |
| BABS        | Local-AI assistant (chat · model browser · learn · memory)                      |
