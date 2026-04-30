# SpliceCraft Changelog

---

## [0.4.5] — 2026-04-30

### Added

- **Inline amino-acid translation in the sequence panel.** Each CDS
  feature now has an extra row of one-letter AA codes drawn at codon
  midpoints, directly above (forward) or below (reverse) its bar. No
  more popping the translation strip in/out — the protein is always
  visible alongside the bases. Wrap-around CDS features (those that
  span the origin) translate correctly across the join.
- **Click an AA letter → highlight that codon's three bases on the
  DNA strand.** Cursor parks at the codon centre; Ctrl+C copies the
  3 bp. Empty cells between AA letters are no-ops by design.
- **Per-strand restriction-cut visualization.** Clicking a sticky
  cutter (EcoRI, HindIII, …) in the lane art now tints the upstream
  bases on each strand blue and the downstream bases red, showing
  the staggered overhang correctly — top and bottom strands carry
  different bg colours over the offset bps.
- **Library search input.** Pre-fills "Search"; clears on focus;
  Enter applies a fuzzy subsequence filter to the visible table
  (collections or plasmids); empty Enter clears the filter and
  restores the prefill.
- **Bottom-strand copy** — Alt+C (and Ctrl+Shift+C as an alias for
  terminals that distinguish it from Ctrl+C) reverse-complements
  the current selection before copying.
- **Enter on the seq cursor** highlights the smallest feature
  enclosing that bp — keyboard equivalent of clicking a feature.
- **Up arrow on the focused map** resets the origin to bp 1
  (keyboard partner to Home).
- **Pure-black UI theme** (`splicecraft-black`) — pinned at startup
  so panels and modals match the logo's true-black backdrop instead
  of textual-dark's near-black greys.
- **Toast notifications carry semantic colour.** Saves/loads/copies
  flash green ("success"); information stays neutral; warnings amber;
  errors red. Notifications fired while the splash is up are queued
  and replayed on dismiss so startup messages aren't lost.

### Changed

- **2D feature-lane packer** replaces the three-tier RE / 1bp /
  multi-bp lane stack. Every feature now sits in lane 0 (adjacent
  to the DNA strand) by default; only bp-range collisions push a
  feature up. Restriction sites participate in the same lanes as
  ordinary features — the parens row prints far from DNA, the cut
  arrow close. Lane depth is uncapped — features pile up as deep
  as the data demands.
- **Layout rework.** The library, plasmid map, and feature sidebar
  share one horizontal top row; the sequence panel sits beneath
  them and spans the full window width. The old per-feature detail
  box in the sidebar was removed (info still surfaces via the row
  + map highlight); the redundant `Sequence` header strip in the
  seq panel is gone too.
- **Map / sidebar feature picks now park the cursor at the feature's
  5' end** rather than its midpoint. Long CDS rows used to land the
  cursor mid-feature; the new behaviour anchors at the feature's
  start so users read top-down.
- **Lane clicks no longer scroll the seq panel.** The user clicked
  something they were already looking at, so jumping the viewport
  away from their cursor would be jarring.
- **Map rotation keys are focus-gated.** `[` / `]` and arrow keys
  rotate only when the plasmid map has focus — they no longer
  fire from modal screens or the seq panel.
- **Arrow keys clear the active RE highlight** and park the cursor
  immediately upstream (Left) or downstream (Right / Up / Down) of
  the top-strand cut.
- **Arrow keys exit a feature highlight** at the matching end and
  step one base in the arrow's direction, instead of being absorbed
  by the highlight.
- **Backbone clicks on the map** (or anywhere outside the four main
  panels) now clear every panel's highlight in one go.
- **Loading a library entry that's already loaded is now a no-op,**
  instead of clobbering undo/redo and any unsaved edits.
- **Performance budgets bumped** for the inline AA row + inter-chunk
  gap (`50 KB cursor ≤ 50 ms`, `150 KB cursor ≤ 120 ms`).

### Fixed

- **Wrap-CDS inline AA painting.** The new AA row was placing letters
  at the wrong bps with the wrong reading frame for any CDS that
  crosses the origin (head halves were translated as if they were
  fresh 0-indexed CDS fragments). `_feats_in_chunk` now stamps the
  original `(start, end)` on each split half as `_orig_start` /
  `_orig_end`, and `_paint_cds_aa` / `_cds_aa_list` / the AA-letter
  click handler all use those for codon math. Regression test in
  `TestWrapCDSInlineTranslation`.
- **`SequencePanel.on_mouse_down` AttributeError on first click.**
  The new lane-click skip-scroll logic read `self._last_lane_click`
  before any prior `on_click` had a chance to initialise it. Now
  set in `__init__` and reset before every `_click_to_bp` call so
  the flag reflects only the current click.

### Removed

- Three vestigial helpers (`_build_chunk_translation`, `_emit_aa_row`,
  `_chunk_has_cds`, ~125 lines) from an earlier AA-row prototype that
  was superseded by `_paint_cds_aa`. No callers, no tests.
- Sequence-panel header strip and translation footer (`#seq-hdr`,
  `#seq-trans`).
- Feature-sidebar detail box (`#detail-box`); `show_detail` is now
  a no-op kept for caller compatibility.

---

## [Unreleased]

### Added

- **SnapGene (.dna) file import** — `File → Open` and the `o` hotkey
  now accept SnapGene's native binary `.dna` format via Biopython's
  built-in parser. No manual GenBank export step required. Files are
  dispatched by extension (`.gb`, `.gbk`, `.genbank` → GenBank;
  `.dna` → SnapGene), case-insensitively. Malformed `.dna` files
  produce a user-friendly error pointing to the likely cause.

### Fixed

- **Golden Braid primer validation** — `_design_gb_primers` now returns a
  clear error when the selected region is shorter than 18 bp, instead of
  silently producing a too-short primer with `Tm=0.0`. `_run_goldenbraid`
  surfaces that error in red in the results pane.
- **pLannotate race condition** — if the user loaded a different plasmid
  while pLannotate was still running, the worker would silently replace
  the newly-loaded plasmid with the merged old one. The worker now checks
  `self._current_record is record` before applying and drops the stale
  result with a warning.
- **Undo stack leaked across plasmid loads** — pressing `Ctrl+Z` after
  switching plasmids could yank the user back to an unrelated edit on
  the previous plasmid. `_apply_record` now clears undo/redo on a fresh
  load (fetch / file open / library pick). In-place record changes
  (pLannotate merge, sequence edits) keep their undo entries intact.
- **Wrap-around restriction sites** — enzymes whose recognition sequence
  spans the origin of a circular plasmid are now found and rendered as
  two linked pieces (labeled tail + unlabeled head). Previously those
  sites were silently invisible.
- **Zero-width feature click detection** — a malformed feature with
  `start == end` used to match every click on the backbone in linear
  view. The linear click handler now shares `_bp_in`'s half-open
  `[start, end)` semantics, making zero-width features unclickable.
- **Shrink-guard widened** — the data-safety guard now logs any library
  shrink (not just nukes to zero entries), making accidental entry
  deletion easier to audit in `/tmp/splicecraft.log`.

### Added

- **Feature deletion** — press `Delete` to remove the selected feature (annotation only,
  sequence is untouched); fully undo/redo-able with `Ctrl+Z` / `Ctrl+Shift+Z`.

- **Toggleable linear map view** — press `v` to switch the circular map panel between
  circular and horizontal linear views.  Linear view uses the same braille-pixel rendering
  with per-strand feature bars, arrowheads, lane stacking, and feature labels.

- **Strand-aware DNA sequence panel layout** — forward-strand features always appear
  *above* the DNA sequence line; reverse-strand features always appear *below*, making
  strand identity immediately apparent.  Overlapping features on the same strand stack
  into additional lanes on their respective side.

- **Braille feature bars in sequence panel** — annotation bars now use solid braille
  block characters (`⣿`) matching the aesthetic of the map viewer, with `▶`/`◀`
  arrowheads at the true start/end of each feature.

- **Single-bp feature triangles** — features that are one base-pair wide render as `▼`
  (above DNA) or `▲` (below DNA), pointing inward toward the sequence line.

- **Label-above / label-below layout** — feature names appear outside the bar (above the
  bar for forward features, below for reverse), keeping the braille bar itself clean.
  Multiple non-overlapping features share a single horizontal row pair.

- **Feature connector lines** (`l` key toggle) — draws a `┊` connector between each
  feature label and its braille bar in the sequence panel, and a dotted radial leader
  line from the arc to the label in the circular map.  Both panels respond to the same
  toggle.

- **Full NEB restriction enzyme catalog** — ~200 enzymes from New England Biolabs,
  including Type IIS (BsaI, BsmBI, BbsI, …) with non-palindromic cut sites.  Each hit
  is visualized as two distinct overlays:
  - **Recognition sequence bar** (`resite`) — thin braille arc outside the backbone for
    forward-strand hits, inside for reverse-strand hits; same strand-above/below layout
    in the sequence panel.
  - **Cut site marker** (`recut`) — `↓` (forward) or `↑` (reverse) arrow in the
    sequence panel; radial `┼` tick on the circular and linear map at the exact cut
    position.  Type IIS cut sites appear displaced from the recognition sequence as
    expected.
  - Recognition sequence IUPAC codes (R, Y, W, S, M, K, B, D, H, V, N) are handled
    via regex; both strands are scanned.  Enzyme labels appear in the circular map
    alongside regular feature labels using the same proximity placement algorithm.

- **Circular map: inside tick marks** — bp graduation marks and labels now sit *inside*
  the backbone ring rather than outside, keeping the outer ring clean for feature labels.
  Two constants (`TICK_DR_MARK`, `TICK_DR_LABEL`) control the inset depth and scale
  automatically with the `,` / `.` aspect-ratio keys.

- **Circular map: full-length feature labels** — removed the 16-character truncation;
  labels now display their full name.

- **Circular map: proximity label placement** — labels are placed as close to the arc as
  possible, greedy-stepping radially outward only when a label would overlap an
  already-placed one.  `LABEL_DR_MIN` (default `9`) sets the minimum clearance.

- **Default library entry** — MW463917.1 (pACYC184) is fetched and added to the library
  automatically on first launch.  The NCBI fetch dialog pre-fills with this accession.

---

## [0.1.0] — 2026-03-23

### Added

- Initial release: braille-canvas circular plasmid map, NCBI live fetch, local `.gb`
  file loading, persistent plasmid library, feature sidebar with CDS translation,
  sequence panel with click-to-cursor, drag selection, undo/redo, and restriction-site
  overlay.
- ASCII logo and README.
