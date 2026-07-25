"""splicecraft_gels — agarose-gel simulation + render (Phase D, layer L1).

The [SUB-gels] simulation core, extracted from the hub: gel-entry normalisation
+ id minting, the agarose-mobility model, per-lane band computation, and the
text/Rich gel-image render. Pure (app-free) — operates on plain gel dicts; the
GelLibraryModal / SimulatorScreen pass them in. Depends only on lower layers:
`_digest_with_enzymes` (biology L0, for restriction-band sizing) + the util
sanitisers / `_now_iso`. Re-exported by the hub so `sc.<name>` + every existing
call site (modals, agent gel endpoints) resolves unchanged.
"""
from __future__ import annotations

import math
import uuid as _uuid

from rich.text import Text

from splicecraft_biology import _digest_with_enzymes
from splicecraft_util import (
    _now_iso, _sanitize_gel_id, _sanitize_label, _sanitize_note,
)


_GEL_NAME_MAX_LEN = 200

_GEL_NOTES_MAX_LEN = 2000

_GEL_LANE_NAME_MAX_LEN = 60

_GEL_LANE_DETAIL_MAX_LEN = 200

_GEL_LANE_SOURCE_MAX_LEN = 64

_GEL_LANES_MAX = 20            # mirrors SimulatorScreen._MAX_LANES

_GEL_AGAROSE_MIN = 0.3

_GEL_AGAROSE_MAX = 5.0

# Empirical agarose-gel resolution windows. Within each window, distance
# migrated is approximately linear in -log10(bp) (the Helling-Goodman-
# Boyer 1974 observation; see Sambrook & Russell, "Molecular Cloning",
# 3e Table 5-1). Outside the window, bands either run with the dye front
# (small fragments) or stick near the well (large fragments).
#
# Keys are agarose percentages (w/v). Values are (bp_lower_resolution,
# bp_upper_resolution). The lower bound is the smallest fragment that
# doesn't run with the dye front; the upper bound is the largest fragment
# that still resolves from the well.
_AGAROSE_RANGES: dict[float, tuple[int, int]] = {
    0.5:  (1000, 30_000),
    0.7:  ( 800, 12_000),
    0.8:  ( 800, 12_000),
    1.0:  ( 500, 10_000),
    1.2:  ( 400,  7_000),
    1.5:  ( 200,  4_000),
    2.0:  ( 100,  2_000),
    2.5:  ( 100,  1_500),
    3.0:  (  50,  1_000),
    4.0:  (  25,    500),
}

_AGAROSE_CHOICES: tuple[float, ...] = tuple(sorted(_AGAROSE_RANGES.keys()))

# Fraction of the lane reserved at EACH end for fragments outside the gel's
# resolution window — sub-resolution bands running past `bp_min` toward the
# dye front, and above-resolution bands sitting behind `bp_max` toward the
# well. `_agarose_mobility` compresses the in-window Helling-Goodman-Boyer
# map into the interior `[_GEL_EDGE_BAND, 1 - _GEL_EDGE_BAND]` so the two
# regions can never overlap (an overlap inverts the size ordering at each
# window boundary).
_GEL_EDGE_BAND = 0.03

# Effective MW multiplier per DNA form. Supercoiled runs faster than
# linear of equal size → effectively migrates as a smaller linear.
# Nicked / open-circle (relaxed) migrates slower → effectively a larger
# linear. The exact factors vary by gel %, voltage, and ionic strength;
# 0.7×/1.4× are the textbook midline values (Lewis & Slater 1986).
_GEL_FORM_FACTOR: dict[str, float] = {
    "linear":      1.0,
    "supercoiled": 0.7,
    "nicked":      1.4,
    "relaxed":     1.4,   # synonym for nicked / open-circle
}

# Standard agarose-gel size ladders. Each entry maps a ladder name to a
# list of band sizes in bp (top-to-bottom on the gel — largest first).
# Curated to span the common bench workflow: NEB-style 1 kb Plus and
# 1 kb for routine cloning; NEB 100 bp for small-fragment work; Lambda
# digests for legacy / large-fragment work.
_GEL_LADDERS: dict[str, list[int]] = {
    "1 kb Plus":  [15000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000,
                   850, 650, 500, 400, 300, 200, 100],
    "1 kb":       [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500,
                   1000, 750, 500, 250],
    "100 bp":     [1517, 1200, 1000, 900, 800, 700, 600, 500, 400,
                   300, 200, 100],
    "Lambda/HindIII": [23130, 9416, 6557, 4361, 2322, 2027, 564, 125],
}

_LADDER_NAMES: tuple[str, ...] = tuple(_GEL_LADDERS.keys())

# Display-parameter caps for `simulate-gel`. Match `_render_gel_image`
# defaults but bound the agent-facing knobs so a hostile body can't
# request a million-row gel.
_GEL_HEIGHT_MIN     = 4

_GEL_HEIGHT_MAX     = 200

_GEL_LANE_WIDTH_MIN = 1

_GEL_LANE_WIDTH_MAX = 32

_GEL_MAX_LANES      = 16   # 2× the in-UI cap of 8; agent flows may batch


def _new_gel_id(existing: "set[str] | None" = None) -> str:
    seen = existing or set()
    for _ in range(64):
        gid = f"gel-{_uuid.uuid4().hex[:8]}"
        if gid not in seen:
            return gid
    return f"gel-{_uuid.uuid4().hex}"


def _normalise_gel_entry(entry: dict, *, fresh: bool = False) -> dict:
    """Normalise: cap name + notes + lanes + lane-fields, stamp
    timestamps, sanitise id, clamp agarose % to a sane range
    (0.3 – 5.0 %). Mirrors `_normalise_experiment_entry`."""
    out = dict(entry) if isinstance(entry, dict) else {}
    gid = _sanitize_gel_id(out.get("id"))
    if gid is None:
        gid = _new_gel_id()
    out["id"] = gid
    raw_name = out.get("name")
    name = raw_name if isinstance(raw_name, str) else ""
    # Strip control bytes (terminal-escape defence) on top of the length
    # cap — gel + lane names render in the gel-picker DataTable.
    out["name"] = (_sanitize_label(name, max_len=_GEL_NAME_MAX_LEN)
                   or "Untitled gel")
    raw_notes = out.get("notes")
    notes = raw_notes if isinstance(raw_notes, str) else ""
    # Sweep #30 (2026-05-28): strip control bytes (preserving \t/\n for
    # multi-line notes) — gel name + lane fields already go through
    # _sanitize_label; notes was a raw slice, so an agent could persist a
    # terminal-escape that renders when the gel is opened. [INV-85]
    out["notes"] = _sanitize_note(notes, max_len=_GEL_NOTES_MAX_LEN)
    try:
        agar = float(out.get("agarose_pct", 1.0))
    except (TypeError, ValueError):
        agar = 1.0
    # NaN / inf rejection — `float("nan") < anything` is False so
    # `max/min` would let them pass through.
    if agar != agar or agar in (float("inf"), float("-inf")):
        agar = 1.0
    out["agarose_pct"] = max(_GEL_AGAROSE_MIN, min(_GEL_AGAROSE_MAX, agar))
    raw_lanes = out.get("lanes") or []
    if not isinstance(raw_lanes, list):
        raw_lanes = []
    lanes: "list[dict]" = []
    for ln in raw_lanes[:_GEL_LANES_MAX]:
        if not isinstance(ln, dict):
            continue
        raw_nm = ln.get("name")
        raw_src = ln.get("source")
        raw_det = ln.get("detail")
        nm  = raw_nm  if isinstance(raw_nm,  str) else ""
        src = raw_src if isinstance(raw_src, str) else "empty"
        det = raw_det if isinstance(raw_det, str) else ""
        lanes.append({
            "name":   _sanitize_label(nm,  max_len=_GEL_LANE_NAME_MAX_LEN),
            "source": (_sanitize_label(src, max_len=_GEL_LANE_SOURCE_MAX_LEN)
                       or "empty"),
            "detail": _sanitize_label(det, max_len=_GEL_LANE_DETAIL_MAX_LEN),
        })
    out["lanes"] = lanes
    now = _now_iso()
    if fresh or not isinstance(out.get("created_at"), str):
        out["created_at"] = now
    out["updated_at"] = now
    return out


def _agarose_mobility(bp: int, gel_pct: float,
                       dna_form: str = "linear") -> float:
    """Return relative mobility in [0, 1]: 0 = at the well (origin), 1
    = at the dye front. Distance migrated on a rendered gel of height
    H rows is `round(mobility * (H - 1))`.

    Within each gel's resolution window, mobility is linear in
    `-log10(eff_bp)` — the Helling-Goodman-Boyer empirical observation
    (Sambrook & Russell, "Molecular Cloning" 3e, Table 5-1).

    Outside the resolution window (refactor 2026-05-19): mobility
    continues to extrapolate linearly in log10 with the SAME slope,
    but with a damped soft-asymptote so very small fragments
    (`bp << bp_min`) stack near the dye front WITHOUT collapsing to
    the same row, and very large fragments (`bp >> bp_max`) stack
    near the well likewise without collapse. Two below-window
    fragments will still order by size (smaller faster); the same
    holds above-window (larger slower). Pre-fix the boundary hard-
    clamped to 0.97 / 0.03 so multiple sub-resolution bands piled on
    the same row regardless of relative size, which lost real
    ordering information.

    `dna_form` ∈ {"linear", "supercoiled", "nicked", "relaxed"}. Unknown
    forms are treated as linear.
    """
    if bp is None or bp <= 0:
        return 1.0
    factor = _GEL_FORM_FACTOR.get(dna_form, 1.0)
    eff_bp = max(1, int(round(bp * factor)))
    # Snap to nearest configured gel %.
    gel_pct = min(_AGAROSE_CHOICES, key=lambda g: abs(g - gel_pct))
    bp_min, bp_max = _AGAROSE_RANGES[gel_pct]
    log_lo = math.log10(bp_min)
    log_hi = math.log10(bp_max)
    log_x  = math.log10(eff_bp)
    # In-window raw mobility, 0..1 across the resolution window.
    raw = (log_hi - log_x) / (log_hi - log_lo)
    if 0.0 <= raw <= 1.0:
        # Compress into the INTERIOR band [_GEL_EDGE_BAND, 1-_GEL_EDGE_BAND].
        # `bp_min` is the smallest RESOLVABLE fragment, not the dye front, and
        # `bp_max` the largest resolvable, not the well — so lane space has to
        # be reserved beyond both for the out-of-window branches below.
        # Pre-fix the in-window map spanned the FULL 0..1 lane while the
        # out-of-window branches were squeezed into [0.97, 0.995] / [0.005,
        # 0.03], so the two regions OVERLAPPED and the order inverted at each
        # boundary: on a 1% gel a 499 bp fragment scored 0.983 against the
        # 500 bp fragment's 1.000, rendering the SMALLER band ABOVE the larger
        # one (and a 10 001 bp band below a 10 000 bp one). Visible on any
        # ladder that straddles a window edge — the 100 bp ladder on a 1% gel
        # put its 500 bp band below 400/300/200/100.
        return _GEL_EDGE_BAND + (1.0 - 2.0 * _GEL_EDGE_BAND) * raw
    # Out-of-window: damped extrapolation into the reserved edge band so
    # multiple below-window (or above-window) fragments retain size ordering
    # instead of hard-clamping onto one row. Each additional log10 unit past
    # the boundary halves the remaining gap to the asymptote (dye front /
    # well), which it approaches but never reaches. `0.5 ** excess` (not
    # `0.5 ** (1 + excess)`) makes the curve CONTINUOUS with the in-window
    # edge — the old form jumped half the band the instant a fragment fell
    # one bp outside the window.
    if raw > 1.0:
        excess = raw - 1.0       # > 0 — below bp_min, runs past the window
        damped = 1.0 - 0.5 ** excess          # 0 → 1
        return (1.0 - _GEL_EDGE_BAND) + _GEL_EDGE_BAND * damped
    # raw < 0.0 — above bp_max, sits behind the window toward the well
    deficit = -raw            # > 0
    damped = 1.0 - 0.5 ** deficit             # 0 → 1
    return _GEL_EDGE_BAND - _GEL_EDGE_BAND * damped


def _gel_bands_for_lane(
    lane:          dict,
    *,
    template_seq:  str,
    template_circular: bool,
    pcr_amplicon: "dict | None",
) -> list[tuple[int, str]]:
    """Resolve a lane's source descriptor into a list of `(bp, form)`
    bands. Pure function for testability — UI rendering separately maps
    each (bp, form) to a row index.

    Returns an empty list for empty / unrecognised sources.
    """
    src    = (lane.get("source") or "empty").lower()
    detail = (lane.get("detail") or "").strip()
    bands: list[tuple[int, str]] = []
    if src == "ladder":
        name = detail if detail in _GEL_LADDERS else _LADDER_NAMES[0]
        for bp in _GEL_LADDERS[name]:
            bands.append((bp, "linear"))
    elif src == "plasmid":
        seq_len = len(template_seq or "")
        if seq_len <= 0:
            return []
        if template_circular:
            # Uncut circular plasmid presents three bands: supercoiled
            # (fastest), linear (rare — from nicking during prep), nicked
            # / open-circle (slowest). Bench reality: a fresh prep is
            # mostly supercoiled with a faint nicked band; show both so
            # the user reads the rendering as a real gel image.
            bands.append((seq_len, "supercoiled"))
            bands.append((seq_len, "nicked"))
        else:
            bands.append((seq_len, "linear"))
    elif src == "digest":
        if not template_seq:
            return []
        enz_list = [e.strip() for e in detail.split(",") if e.strip()]
        if not enz_list:
            return []
        try:
            frags = _digest_with_enzymes(template_seq, enz_list,
                                          circular=template_circular)
        except (ValueError, KeyError, RuntimeError):
            return []
        for f in frags:
            bp = len(f.get("top_seq", "") or "")
            if bp > 0:
                bands.append((bp, "linear"))
    elif src == "pcr":
        # A `pcr` lane can freeze its own amplicon size in `_pcr_bp`
        # (stamped by "Send to Gel lane" so multiple amplicons coexist
        # on one gel). Prefer it; fall back to the screen's currently-
        # selected amplicon for a `pcr` lane the user added manually via
        # the source dropdown (no frozen size).
        bp = 0
        frozen = lane.get("_pcr_bp")
        if isinstance(frozen, int) and not isinstance(frozen, bool) \
                and frozen > 0:
            bp = frozen
        elif isinstance(pcr_amplicon, dict):
            # Defensive: agent endpoint accepts an arbitrary dict for
            # `pcr_amplicon`; a hostile / malformed payload could carry
            # a non-numeric `length`. `int()` on the bad value would
            # surface as a 500 — better to render an empty lane than
            # crash the gel.
            try:
                bp = int(pcr_amplicon.get("length", 0))
            except (TypeError, ValueError):
                bp = 0
        if bp > 0:
            bands.append((bp, "linear"))
    return bands


def _render_gel_image(
    lane_specs:    list[dict],
    *,
    template_seq:  str,
    template_circular: bool,
    pcr_amplicon: "dict | None",
    agarose_pct:   float,
    height:        int = 22,
    lane_width:    int = 7,
    label_col:     int = 7,
) -> Text:
    """Render the gel as Rich `Text`. One column per lane, well-at-top
    to dye-front-at-bottom. Returns a ready-to-`Static.update()` Text.

    `lane_specs` is the live list of lane config dicts from the gel
    tab. `pcr_amplicon` is the currently-selected PCR result from the
    PCR tab (or None if no PCR has been run). All migration math
    routes through `_agarose_mobility`.
    """
    rt = Text()
    n_lanes = len(lane_specs)
    if n_lanes == 0:
        rt.append("(no lanes — add at least one to render a gel)\n",
                   style="dim italic")
        return rt
    # Resolve lane bands.
    lane_bands: list[list[tuple[int, str]]] = []
    ladder_lane_idx = -1
    for li, lane in enumerate(lane_specs):
        bands = _gel_bands_for_lane(
            lane,
            template_seq=template_seq,
            template_circular=template_circular,
            pcr_amplicon=pcr_amplicon,
        )
        lane_bands.append(bands)
        if (lane.get("source") or "").lower() == "ladder" and ladder_lane_idx == -1:
            ladder_lane_idx = li

    # Pre-compute row indices for each band. The position is a
    # FLOATING-POINT row index — `mob * (height - 1)` — and we keep
    # the fractional part so adjacent bands whose log10(bp) differ
    # by less than 1 row can still resolve via a faint
    # "anti-aliased" tail on the adjacent row.
    #
    # Two grids:
    #   * `band_grid` — primary cells. `(row, lane) → count`. Heavy
    #                   `━` (1) / `▆` (2) / `█` (3+) glyph chosen
    #                   from count.
    #   * `band_faint` — secondary cells where a single band's
    #                   fractional position leans toward this row.
    #                   `(row, lane) → True` (presence only — no
    #                   pile-up logic). Renders as a light `─` glyph
    #                   so the eye reads the band's true position
    #                   as "between rows".
    band_grid: dict[tuple[int, int], int] = {}   # (row, lane) → band count
    band_faint: set[tuple[int, int]] = set()
    ladder_rows: dict[int, int] = {}             # row → bp size
    # Sub-row fractional rendering kicks in once the band's offset
    # from row-center exceeds this threshold — below it the band
    # sits squarely on its row and a faint tail would only add
    # visual noise (refactor 2026-05-19).
    _FAINT_FRAC_THRESHOLD = 0.25
    for li, bands in enumerate(lane_bands):
        for bp, form in bands:
            mob = _agarose_mobility(bp, agarose_pct, dna_form=form)
            row_float = mob * (height - 1)
            row = max(0, min(height - 1, int(round(row_float))))
            band_grid[(row, li)] = band_grid.get((row, li), 0) + 1
            if li == ladder_lane_idx:
                # The largest bp at this row wins the tick label.
                ladder_rows[row] = max(ladder_rows.get(row, 0), bp)
            # Fractional offset from the row's center.
            frac = row_float - row
            if abs(frac) > _FAINT_FRAC_THRESHOLD:
                # Pull a faint tail toward the adjacent row.
                row_sec = row + (1 if frac > 0 else -1)
                if 0 <= row_sec < height:
                    band_faint.add((row_sec, li))

    # Header row: lane numbers.
    head = " " * label_col
    for li in range(n_lanes):
        head += f"{li+1:^{lane_width}} "
    rt.append(head.rstrip() + "\n", style="bold white")
    # Lane names (truncated/padded).
    names = " " * label_col
    for li in range(n_lanes):
        label = (lane_specs[li].get("name") or "")[:lane_width]
        names += f"{label:^{lane_width}} "
    rt.append(names.rstrip() + "\n", style="cyan")
    # Wells row.
    wells = " " * label_col
    for li in range(n_lanes):
        wells += "█" * lane_width + " "
    rt.append(wells.rstrip() + "\n", style="grey50")
    # Body rows.
    for row in range(height):
        if row in ladder_rows:
            bp = ladder_rows[row]
            if bp >= 1000:
                label = f"{bp/1000:>4.1f}k"
            else:
                label = f"{bp:>5}"
            # Padded to exactly `label_col` chars so labelled rows
            # don't shift the lane columns relative to unlabelled
            # rows. Pre-fix `f"{label} "` was 6 chars and `label_col`
            # was 7 — bands jumped one column left whenever a ladder
            # label landed on the same row, breaking lane-to-well
            # alignment (refactor 2026-05-19).
            line_left = f"{label} ".ljust(label_col)
        else:
            line_left = " " * label_col
        line = line_left
        for li in range(n_lanes):
            count = band_grid.get((row, li), 0)
            if count == 0:
                # Empty cell — but if a band's fractional position
                # leans toward this row, render a LIGHT glyph for
                # sub-row resolution (refactor 2026-05-19). Adjacent
                # bands separated by less than a full row would
                # otherwise collapse on rounding; the faint tail
                # surfaces the true position as "between rows".
                if (row, li) in band_faint:
                    line += "─" * lane_width + " "
                else:
                    line += " " * lane_width + " "
            elif count == 1:
                # Thin horizontal-line glyph for a single band keeps
                # the gel image readable — solid blocks read as a
                # wall rather than a band (user UX call 2026-05-19).
                line += "━" * lane_width + " "
            elif count == 2:
                line += "▆" * lane_width + " "
            else:
                line += "█" * lane_width + " "
        rt.append(line.rstrip() + "\n",
                   style="bright_white" if line_left.strip() else "white")
    # Dye-front row.
    front = " " * label_col
    for li in range(n_lanes):
        front += "░" * lane_width + " "
    rt.append(front.rstrip(), style="dim cyan")
    # Surface a hint about PCR lanes that have no amplicon — the empty
    # lane is otherwise indistinguishable from a digest that failed or
    # a plasmid lane the user forgot to configure. Caller hint, not a
    # `_gel_bands_for_lane` change — the pure function still returns
    # an empty list per its documented contract.
    pcr_empty = [
        li + 1 for li, lane in enumerate(lane_specs)
        if (lane.get("source") or "").lower() == "pcr"
        and not lane_bands[li]
    ]
    if pcr_empty:
        nums  = ", ".join(str(n) for n in pcr_empty)
        label = "lanes" if len(pcr_empty) > 1 else "lane"
        rt.append(
            f"\n  PCR {label} {nums}: no amplicon — run a PCR first to populate.",
            style="dim italic yellow",
        )
    return rt
