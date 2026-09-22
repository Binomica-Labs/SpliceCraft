"""Circular plasmid-map image export (layer 1) — SVG + PNG.

A pure, app-free renderer that turns a plasmid (a normalised feature list
+ length, or a ``SeqRecord``) into a publication-quality circular map:

* the backbone circle,
* colour-coded feature arcs (arrowheads for stranded features), stacked
  onto concentric lanes when they overlap,
* feature names drawn ON the arc (rotated to the tangent) when they fit,
  else on an outer ellipse hugging the circle with a short radial leader,
* restriction-site ticks just outside the backbone, their enzyme names
  joining the same outer-ellipse label system,
* the plasmid name + bp count centred inside the circle when it fits the
  clear inner diameter, otherwise dropped just below the map.

Every external label auto-shrinks its font to clear the canvas edge, and
a dense plasmid's labels are de-collided per side and adaptively capped to
what physically fits the canvas height (the overflow is dropped + logged,
never silently truncated).

Two rendering back-ends consume ONE shared geometry pass
(:func:`_build_primitives`):

* :func:`render_plasmid_map_svg` — pure-stdlib string emission. True
  vector, infinitely scalable, and a transparent background is just the
  absence of the backing rect — ideal for figure editing.
* :func:`render_plasmid_map_png` — rasterised with Pillow (a hard
  dependency of the app). Drawn at a super-sampled resolution and
  down-scaled with LANCZOS so the aliased arcs/lines come out smooth.

The module is deliberately decoupled from ``PlasmidApp``: the live map
hands us its already-parsed ``pm._feats`` / ``pm._restr_feats`` so the
export matches the screen exactly, while :func:`_map_feats_from_record`
re-derives the same shape from a bare record for bulk/off-screen export.

Layer 1: imports biology / util / logging (L0) only — NOT the hub.
"""
from __future__ import annotations

import math
from pathlib import Path
from xml.sax.saxutils import escape as _xml_escape

import splicecraft_biology as _bio
from splicecraft_logging import _log
from splicecraft_util import _feat_label_full, _feature_traversal, _sanitize_label

# ── Tunables ─────────────────────────────────────────────────────────────────

_MAP_IMAGE_DEFAULT_SIZE = 1400   # px, square canvas (SVG viewBox units too)
_MAP_IMAGE_MIN_SIZE     = 300
_MAP_IMAGE_MAX_SIZE     = 6000   # guard against a runaway allocation
_PNG_SUPERSAMPLE        = 2      # draw at N× then LANCZOS-downscale for AA

# Site-tick cap keeps dense plasmids legible instead of a hairball of ticks.
# Feature/site NAME labels are capped adaptively per side by how many physically
# fit the canvas height (see `_place_radial_labels`), not by a fixed number.
_MAP_MAX_SITE_TICKS     = 60
_MAP_MAX_LANES          = 14   # bounds lane assignment on pathological plasmids
_LABEL_MAX_CHARS        = 28   # adaptive font shrinks to fit, so allow long names

# Radii + thicknesses as a fraction of the canvas side. Features stack INWARD
# from the backbone; every external label sits on an OUTER ELLIPSE hugging the
# circle, each joined to its feature/site by a leader. The circle is kept
# COMPACT so the outer zone holds radial labels + long text WITHOUT clipping
# the canvas edge, and leaders stay short + radial so they don't cross.
_R_BACKBONE   = 0.230
_BAND         = 0.046   # feature band thickness
_LANE_STEP    = 0.056   # inward radius drop per stacked lane
_R_MIN_LANE   = 0.088   # innermost a feature band may reach (protect centre)
_ARROW_HEAD   = 0.020   # arrowhead angular reach, radians (clamped to arc)
_SITE_TICK    = 0.018   # tick length just outside the backbone
# Outer label ellipse: taller than wide so the full canvas height is usable for
# stacked labels while long text still clears the left/right edges. Text sits at
# the ellipse and extends OUTWARD; the font auto-shrinks so it never clips.
_LABEL_ELLIPSE_A = 0.290   # horizontal semi-axis (label inner edge at mid-height)
_LABEL_ELLIPSE_B = 0.455   # vertical semi-axis
_LABEL_BASE_FS   = 0.0165  # label font size (fraction of s); shrinks to fit
_CHAR_W          = 0.56    # mean glyph width / font-size (label-fit estimate)
_MAP_MARGIN      = 0.022   # min clearance to the canvas edge (fraction of s)

_BACKBONE_COLOR = "#5A6472"
_SITE_COLOR     = "#FF69B4"   # mirrors hub _RESTR_SITE_FEATURE_COLOR
_TITLE_COLOR    = "#1A1D24"
_SUBTITLE_COLOR = "#5A6472"
_BG_COLOR       = "#FFFFFF"

# Mirror of splicecraft_widgets._FEATURE_PALETTE (xterm-256 indices). Kept
# inline so the pure renderer needn't import the Textual widget layer; a
# test (test_mapimage.py::test_palette_matches_widgets) guards against drift.
_FALLBACK_PALETTE: list[str] = [
    "color(39)",   "color(118)",  "color(208)",  "color(213)",  "color(51)",
    "color(220)",  "color(196)",  "color(46)",   "color(201)",  "color(129)",
    "color(166)",  "color(33)",   "color(226)",  "color(160)",  "color(87)",
    "color(105)",  "color(154)",  "color(203)",  "color(81)",   "color(185)",
]

# Export formats this module handles: fmt -> (extension, human label).
_MAP_EXPORT_FORMATS: dict = {
    "png": (".png", "PNG image"),
    "svg": (".svg", "SVG vector"),
}


# ── Colour handling ──────────────────────────────────────────────────────────

_NAMED_RGB = {
    "black": (0, 0, 0), "white": (255, 255, 255), "red": (255, 0, 0),
    "green": (0, 128, 0), "blue": (0, 0, 255), "yellow": (255, 255, 0),
    "cyan": (0, 255, 255), "magenta": (255, 0, 255), "gray": (128, 128, 128),
    "grey": (128, 128, 128), "orange": (255, 165, 0), "purple": (128, 0, 128),
    "pink": (255, 105, 180), "brown": (165, 42, 42), "navy": (0, 0, 128),
}


def _xterm_to_rgb(n: int) -> "tuple[int, int, int]":
    """Convert an xterm-256 palette index to a 24-bit RGB tuple (standard
    16 system colours + 6×6×6 cube + 24-step grayscale ramp)."""
    n = max(0, min(255, int(n)))
    if n < 16:
        base = [
            (0, 0, 0), (128, 0, 0), (0, 128, 0), (128, 128, 0),
            (0, 0, 128), (128, 0, 128), (0, 128, 128), (192, 192, 192),
            (128, 128, 128), (255, 0, 0), (0, 255, 0), (255, 255, 0),
            (0, 0, 255), (255, 0, 255), (0, 255, 255), (255, 255, 255),
        ]
        return base[n]
    if n < 232:
        n -= 16
        r, g, b = n // 36, (n % 36) // 6, n % 6
        conv = lambda v: 0 if v == 0 else 40 * v + 55
        return conv(r), conv(g), conv(b)
    v = 8 + (n - 232) * 10
    return v, v, v


def _to_rgb(color: "str | None") -> "tuple[int, int, int]":
    """Best-effort conversion of any SpliceCraft colour token to RGB.

    Accepts ``#RGB`` / ``#RRGGBB`` / ``color(N)`` / bare xterm index /
    ``rgb(r,g,b)`` / a small set of CSS names. Falls back to a neutral
    grey so a malformed colour never aborts a render."""
    if not isinstance(color, str):
        return (150, 150, 150)
    c = color.strip().lower()
    if not c:
        return (150, 150, 150)
    if c in _NAMED_RGB:
        return _NAMED_RGB[c]
    if c.startswith("#"):
        h = c[1:]
        if len(h) == 3:
            h = "".join(ch * 2 for ch in h)
        if len(h) == 6:
            try:
                return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
            except ValueError:
                return (150, 150, 150)
        return (150, 150, 150)
    if c.startswith("rgb(") and c.endswith(")"):
        try:
            parts = [int(p) for p in c[4:-1].split(",")]
            if len(parts) == 3:
                return tuple(max(0, min(255, p)) for p in parts)  # type: ignore[return-value]
        except ValueError:
            return (150, 150, 150)
    if c.startswith("color(") and c.endswith(")"):
        try:
            return _xterm_to_rgb(int(c[6:-1]))
        except ValueError:
            return (150, 150, 150)
    if c.isdigit():
        return _xterm_to_rgb(int(c))
    return (150, 150, 150)


def _hex(rgb: "tuple[int, int, int]") -> str:
    return "#%02X%02X%02X" % rgb


def _luminance(rgb: "tuple[int, int, int]") -> float:
    r, g, b = (v / 255.0 for v in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _ensure_readable(rgb: "tuple[int, int, int]") -> "tuple[int, int, int]":
    """Darken a very-light colour so label text stays legible on white.
    Leaves already-readable colours untouched."""
    if _luminance(rgb) <= 0.62:
        return rgb
    return tuple(int(v * 0.55) for v in rgb)  # type: ignore[return-value]


# ── Geometry ─────────────────────────────────────────────────────────────────

def _angle(bp: float, total: int, origin_bp: int) -> float:
    """bp → radians, matching the on-screen map: origin at the top of the
    circle (−π/2), increasing bp sweeps clockwise."""
    if total <= 0:
        return -math.pi / 2
    return 2 * math.pi * ((bp - origin_bp) % total) / total - math.pi / 2


def _in_arc(x: float, start: float, length: float, total: int) -> bool:
    """Does bp ``x`` fall inside the circular arc ``[start, start+length)``
    on a plasmid of ``total`` bp? Wrap-safe."""
    if total <= 0:
        return False
    return ((x - start) % total) < length


def _assign_lanes(feats: list, total: int) -> "dict[int, int]":
    """Greedy circular-arc colouring: map each feature index to a lane so
    overlapping arcs stack onto separate concentric rings. Returns
    ``{feat_index: lane}`` (lane 0 = on the backbone, higher = inward).

    Longest arcs are placed first so the dominant features own the outer
    lanes — matches how the eye reads a plasmid map."""
    order = sorted(
        range(len(feats)),
        key=lambda i: _bio._feat_len(feats[i]["start"], feats[i]["end"], total),
        reverse=True,
    )
    lanes: list[list[tuple[int, int]]] = []   # per lane: list of (start, arclen)
    result: dict[int, int] = {}
    for i in order:
        f = feats[i]
        s = int(f["start"])
        ln = _bio._feat_len(s, f["end"], total)
        placed = False
        for lane_idx, occ in enumerate(lanes):
            if all(
                not (_in_arc(s, os_, ol, total)
                     or _in_arc(os_, s, ln, total))
                for os_, ol in occ
            ):
                occ.append((s, ln))
                result[i] = lane_idx
                placed = True
                break
        if not placed:
            if len(lanes) >= _MAP_MAX_LANES:
                # Cap the ring count: a plasmid with a huge overlapping-feature
                # pile-up would otherwise spawn unbounded (invisibly thin) lanes
                # and make lane assignment O(F²). Stack the overflow onto the
                # innermost lane (accepting some visual overlap) instead.
                lanes[-1].append((s, ln))
                result[i] = len(lanes) - 1
            else:
                lanes.append([(s, ln)])
                result[i] = len(lanes) - 1
    return result


def _lane_radius(lane: int, size: float) -> float:
    """Centre radius (px) of feature ``lane``, clamped so deep stacks don't
    collide with the centre title block."""
    r = (_R_BACKBONE - lane * _LANE_STEP) * size
    return max(r, _R_MIN_LANE * size)


def _arc_points(cx: float, cy: float, r: float,
                a0: float, a1: float, steps: int) -> list:
    """Sample a circular arc from angle ``a0`` to ``a1`` (a1 ≥ a0) into
    ``steps``+1 points at radius ``r``."""
    if steps < 1:
        steps = 1
    out = []
    for k in range(steps + 1):
        a = a0 + (a1 - a0) * (k / steps)
        out.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    return out


def _band_polygon(cx: float, cy: float, r_mid: float, band: float,
                  a0: float, a1: float, strand: int) -> list:
    """Closed polygon (list of (x,y)) for a feature arc band from angle
    ``a0`` to ``a1`` at centre radius ``r_mid``, thickness ``band``. A
    stranded feature grows an arrowhead at its leading end."""
    r_out = r_mid + band / 2.0
    r_in = r_mid - band / 2.0
    span = a1 - a0
    steps = max(2, int(span / (2 * math.pi) * 240))
    head = min(_ARROW_HEAD, span * 0.5) if strand else 0.0

    if strand < 0:
        # Arrowhead at the low-bp (a0) end, pointing counter-clockwise.
        tip = (cx + r_mid * math.cos(a0), cy + r_mid * math.sin(a0))
        body_out = _arc_points(cx, cy, r_out, a0 + head, a1, steps)
        body_in = _arc_points(cx, cy, r_in, a0 + head, a1, steps)
        pts = [tip] + body_out + list(reversed(body_in))
        return pts
    if strand > 0:
        # Arrowhead at the high-bp (a1) end, pointing clockwise.
        tip = (cx + r_mid * math.cos(a1), cy + r_mid * math.sin(a1))
        body_out = _arc_points(cx, cy, r_out, a0, a1 - head, steps)
        body_in = _arc_points(cx, cy, r_in, a0, a1 - head, steps)
        pts = body_out + [tip] + list(reversed(body_in))
        return pts
    # Unstranded — a blunt band.
    body_out = _arc_points(cx, cy, r_out, a0, a1, steps)
    body_in = _arc_points(cx, cy, r_in, a0, a1, steps)
    return body_out + list(reversed(body_in))


def _spread_labels(items: list, min_gap: float, lo: float, hi: float) -> None:
    """In-place vertical de-collision. ``items`` is a list of dicts each
    carrying a ``y`` (desired) — nudged apart so neighbours are ≥
    ``min_gap`` apart, kept within [lo, hi]. Sorted by y first."""
    if not items:
        return
    items.sort(key=lambda d: d["y"])
    # Forward pass: push down.
    for i in range(1, len(items)):
        if items[i]["y"] - items[i - 1]["y"] < min_gap:
            items[i]["y"] = items[i - 1]["y"] + min_gap
    # If we overflowed the bottom, shift the whole stack up and clamp.
    overflow = items[-1]["y"] - hi
    if overflow > 0:
        for d in items:
            d["y"] -= overflow
    if items[0]["y"] < lo:
        shift = lo - items[0]["y"]
        for d in items:
            d["y"] += shift


# ── Primitive model (shared by both back-ends) ───────────────────────────────

def _fmt_bp(n: int) -> str:
    return f"{int(n):,} bp"


def _build_primitives(feats: list, total: int, *,
                      title: str = "", sites: "list | None" = None,
                      size: int = _MAP_IMAGE_DEFAULT_SIZE,
                      origin_bp: int = 0,
                      show_labels: bool = True,
                      show_sites: bool = True) -> dict:
    """Turn the plasmid data into back-end-agnostic drawing primitives.

    Returns a dict with keys ``size``, ``circles``, ``polys``, ``lines``,
    ``texts`` — each a list of tuples the SVG/PNG emitters iterate. Kept
    pure so it is trivially unit-testable without any raster/vector lib.
    """
    s = float(size)
    cx = cy = s / 2.0
    circles: list = []
    polys: list = []
    lines: list = []
    texts: list = []

    # Backbone.
    circles.append((cx, cy, _R_BACKBONE * s, _BACKBONE_COLOR,
                    max(1.5, s * 0.0022), None))

    external: list = []   # labels that don't fit on their band → outer ellipse
    inner_clear_r = _R_BACKBONE * s   # clear centre radius for the title
    if total > 0 and feats:
        lanes = _assign_lanes(feats, total)
        deepest = max(lanes.values()) if lanes else 0
        inner_clear_r = _lane_radius(deepest, s) - (_BAND / 2.0) * s
        for i, f in enumerate(feats):
            start = int(f["start"])
            arclen = _bio._feat_len(start, f["end"], total)
            if arclen <= 0:
                continue
            a0 = _angle(start, total, origin_bp)
            a1 = a0 + 2 * math.pi * (arclen / total)
            lane = lanes.get(i, 0)
            r_mid = _lane_radius(lane, s)
            strand = int(f.get("strand") or 0)
            rgb = _to_rgb(f.get("color"))
            poly = _band_polygon(cx, cy, r_mid, _BAND * s, a0, a1, strand)
            polys.append((poly, _hex(rgb), "#00000022", max(0.5, s * 0.0009)))

            if not show_labels:
                continue
            label = _sanitize_label(str(f.get("label") or f.get("type") or ""),
                                    max_len=_LABEL_MAX_CHARS)
            if not label:
                continue
            a_mid = a0 + (a1 - a0) / 2.0
            span = a1 - a0
            # Name drawn ON the arc when it fits the band's length — rotated to
            # the tangent, coloured for contrast against the band.
            band_px = _BAND * s
            onband_fs = min(_LABEL_BASE_FS * s, band_px * 0.66)
            if span > 0 and (r_mid * span) >= len(label) * onband_fs * _CHAR_W * 1.12:
                rot = math.degrees(a_mid) + 90.0
                while rot > 90.0:
                    rot -= 180.0
                while rot < -90.0:
                    rot += 180.0
                txt_c = "#101010" if _luminance(rgb) > 0.55 else "#FFFFFF"
                texts.append((cx + r_mid * math.cos(a_mid),
                              cy + r_mid * math.sin(a_mid),
                              label, onband_fs, txt_c, "middle", "bold", rot))
            else:
                # Anchor the leader at the ring perimeter (not the band's own
                # radius) so stacked inner-lane features in a cluster fan out
                # cleanly from the circle instead of converging on the centre.
                anchor_r = max(r_mid + (_BAND / 2.0) * s,
                               _R_BACKBONE * s + _SITE_TICK * s * 0.4)
                external.append({
                    "a": a_mid,
                    "ax": cx + anchor_r * math.cos(a_mid),
                    "ay": cy + anchor_r * math.sin(a_mid),
                    "label": label, "color": _hex(_ensure_readable(rgb)),
                })

    # Restriction ticks. Each unique cut's NAME joins the outer radial label
    # system so ticks + feature labels de-collide together (no separate pile).
    if show_sites and sites and total > 0:
        seen: set = set()
        nsites = 0
        for rf in sites:
            # Keep CUT positions only: a live map hands us its full
            # `_restr_feats` (recognition-span "resite" + single-bp "recut"
            # dicts); a pre-filtered list (no `type` key) is taken as-is.
            if rf.get("type") not in (None, "recut"):
                continue
            enzyme = _sanitize_label(str(rf.get("label") or ""),
                                     max_len=_LABEL_MAX_CHARS)
            pos = int(rf.get("start", 0))
            key = (enzyme, pos)
            if not enzyme or key in seen:
                continue
            seen.add(key)
            a = _angle(pos, total, origin_bp)
            r_out = _R_BACKBONE * s + _SITE_TICK * s
            lines.append((cx + _R_BACKBONE * s * math.cos(a),
                          cy + _R_BACKBONE * s * math.sin(a),
                          cx + r_out * math.cos(a), cy + r_out * math.sin(a),
                          _SITE_COLOR, max(1.0, s * 0.0015)))
            if show_labels:
                external.append({
                    "a": a, "ax": cx + r_out * math.cos(a),
                    "ay": cy + r_out * math.sin(a),
                    "label": enzyme, "color": _SITE_COLOR,
                })
            nsites += 1
            if nsites >= _MAP_MAX_SITE_TICKS:
                break

    if show_labels and external:
        dropped = _place_radial_labels(external, cx, cy, s, lines, texts)
        if dropped:
            _log.info("map image: %d labels dropped (canvas full)", dropped)

    _place_title(title, total, cx, cy, s, inner_clear_r, texts)

    return {"size": size, "circles": circles, "polys": polys,
            "lines": lines, "texts": texts}


def _place_radial_labels(items: list, cx: float, cy: float, s: float,
                         lines: list, texts: list) -> int:
    """Lay external labels around the OUTSIDE on a tall ellipse, each on the
    LEFT or RIGHT half (by its anchor), stacked in angular order and
    de-collided VERTICALLY — text extends outward into the free side-space so
    labels never overlap. Each label sits at its anchor's own angular height,
    so its leader is short + near-radial; and because a side's leaders keep the
    shared angular→vertical order they stay nested and don't cross. The font
    auto-shrinks so the longest clears the canvas edge and the side fits the
    height with slack (no jam-packed overflow that would drag a cluster away
    from its feature); excess is dropped EVENLY by angle so the rest keep their
    spread. Returns the count dropped."""
    A = _LABEL_ELLIPSE_A * s
    B = _LABEL_ELLIPSE_B * s
    margin = _MAP_MARGIN * s
    usable = s - 2.0 * margin
    dropped = 0
    for direction in (+1, -1):
        side = [d for d in items
                if (math.cos(d["a"]) >= 0) == (direction > 0)]
        if not side:
            continue
        n_side = len(side)
        maxchars = max(len(d["label"]) for d in side)
        # Font that (1) lets the longest label reach from its inner edge
        # (cx ± A, mid-height) to the margin and (2) fits the whole side's
        # stack into the height with a little slack.
        room = (s - margin) - (cx + A) if direction > 0 else (cx - A) - margin
        fs = _LABEL_BASE_FS * s
        if maxchars > 0 and room > 0:
            fs = min(fs, room / (maxchars * _CHAR_W))
        fs = min(fs, usable / (n_side * 1.62))
        fs = max(fs, s * 0.0072)
        lh = fs * 1.62
        # Order the side top→bottom by the anchor's vertical height. Both the
        # anchors and the labels then share this order, so their leaders stay
        # nested (never cross within a side).
        side.sort(key=lambda d: math.sin(d["a"]))
        # Cap with slack: leave ~6% of the height free so de-collision only
        # nudges neighbours apart instead of packing edge-to-edge and shifting
        # the whole stack (which is what drags a dense cluster's labels far
        # from their features → long crossing leaders).
        cap = max(1, int(usable / lh * 0.94))
        if n_side > cap:
            side = [side[(k * n_side) // cap] for k in range(cap)]
            dropped += n_side - cap
        for d in side:
            d["y"] = cy + B * math.sin(d["a"])
        _spread_labels(side, lh, margin + lh / 2.0, s - margin - lh / 2.0)
        anchor = "start" if direction > 0 else "end"
        for d in side:
            y = d["y"]
            t = max(-1.0, min(1.0, (y - cy) / B))
            x = cx + direction * A * math.sqrt(max(0.0, 1.0 - t * t))
            lines.append((d["ax"], d["ay"], x, y, d["color"],
                          max(0.5, s * 0.0009)))
            texts.append((x + direction * fs * 0.3, y, d["label"], fs,
                          d["color"], anchor, "normal", 0.0))
    return dropped


def _place_title(title: str, total: int, cx: float, cy: float, s: float,
                 inner_clear_r: float, texts: list) -> None:
    """Plasmid name + bp. Centred INSIDE the circle when the name fits the
    clear inner diameter; otherwise moved BELOW the circle (name, then bp),
    shrunk to fit the canvas width. bp always sits directly under the name."""
    name = _sanitize_label(title or "", max_len=60)
    bp = _fmt_bp(total) if total > 0 else ""
    margin = _MAP_MARGIN * s
    inner_d = max(0.0, 2.0 * inner_clear_r) * 0.9
    name_fs = s * 0.030
    if name and len(name) * name_fs * _CHAR_W <= inner_d:
        texts.append((cx, cy - s * 0.016, name, name_fs, _TITLE_COLOR,
                      "middle", "bold", 0.0))
        if bp:
            texts.append((cx, cy + s * 0.026, bp, s * 0.021, _SUBTITLE_COLOR,
                          "middle", "normal", 0.0))
        return
    # Below the circle: name shrunk to span at most the canvas width, bp under.
    y = cy + _R_BACKBONE * s + s * 0.06
    if name:
        nfs = min(s * 0.030, (s - 2.0 * margin) / max(1, len(name) * _CHAR_W))
        nfs = max(nfs, s * 0.012)
        texts.append((cx, y, name, nfs, _TITLE_COLOR, "middle", "bold", 0.0))
        y += nfs * 1.6
    if bp:
        texts.append((cx, y, bp, s * 0.020, _SUBTITLE_COLOR,
                      "middle", "normal", 0.0))


# ── SVG back-end ─────────────────────────────────────────────────────────────

_SVG_FONT = ("-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, "
             "'Helvetica Neue', Arial, sans-serif")


def render_plasmid_map_svg(feats: list, total: int, *, title: str = "",
                           sites: "list | None" = None,
                           size: int = _MAP_IMAGE_DEFAULT_SIZE,
                           origin_bp: int = 0, transparent: bool = False,
                           show_labels: bool = True,
                           show_sites: bool = True) -> str:
    """Render the plasmid map as an SVG document string (pure stdlib)."""
    size = _clamp_size(size)
    prims = _build_primitives(
        feats, total, title=title, sites=sites, size=size,
        origin_bp=origin_bp, show_labels=show_labels, show_sites=show_sites,
    )
    parts: list = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" '
        f'height="{size}" viewBox="0 0 {size} {size}" '
        f'font-family="{_xml_escape(_SVG_FONT)}">'
    ]
    if not transparent:
        parts.append(f'<rect width="{size}" height="{size}" fill="{_BG_COLOR}"/>')
    for cx, cy, r, stroke, w, fill in prims["circles"]:
        parts.append(
            f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r:.2f}" '
            f'fill="{fill or "none"}" stroke="{stroke}" '
            f'stroke-width="{w:.2f}"/>'
        )
    for poly, fill, stroke, w in prims["polys"]:
        pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in poly)
        parts.append(
            f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" '
            f'stroke-width="{w:.2f}" stroke-linejoin="round"/>'
        )
    for x0, y0, x1, y1, color, w in prims["lines"]:
        parts.append(
            f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" '
            f'stroke="{color}" stroke-width="{w:.2f}"/>'
        )
    for x, y, text, fs, color, anchor, weight, rot in prims["texts"]:
        tr = f' transform="rotate({rot:.2f} {x:.2f} {y:.2f})"' if rot else ""
        parts.append(
            f'<text x="{x:.2f}" y="{y:.2f}" font-size="{fs:.2f}" '
            f'fill="{color}" text-anchor="{anchor}" '
            f'dominant-baseline="central" '
            f'font-weight="{weight}"{tr}>{_xml_escape(text)}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


# ── PNG back-end (Pillow) ────────────────────────────────────────────────────

def _load_font(px: int):
    """A TrueType font at ``px`` pixels. Prefers DejaVu (ships with Pillow),
    then Pillow's scalable default; never raises."""
    from PIL import ImageFont
    px = max(6, int(px))
    for name in ("DejaVuSans.ttf", "Arial.ttf", "Helvetica.ttf"):
        try:
            return ImageFont.truetype(name, px)
        except OSError:
            continue
    try:
        return ImageFont.load_default(size=px)
    except TypeError:      # Pillow < 10 — non-scalable bitmap default
        return ImageFont.load_default()


_PIL_ANCHOR = {"start": "lm", "middle": "mm", "end": "rm"}


def _draw_rotated_text(base_img, x: float, y: float, text: str, font,
                       fill, rot_deg: float) -> None:
    """Draw ``text`` centred at ``(x, y)`` rotated ``rot_deg``° onto the RGBA
    ``base_img`` via a temp text tile → rotate → alpha-composite. ``rot_deg``
    is SVG-positive (clockwise), so the PIL rotate (CCW-positive) is negated,
    keeping the two back-ends visually identical for on-arc feature names."""
    from PIL import Image, ImageDraw
    probe = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
    left, top, right, bottom = probe.textbbox((0, 0), text, font=font)
    tw = max(1, int(math.ceil(right - left)))
    th = max(1, int(math.ceil(bottom - top)))
    pad = 2
    tile = Image.new("RGBA", (tw + 2 * pad, th + 2 * pad), (0, 0, 0, 0))
    ImageDraw.Draw(tile).text((pad - left, pad - top), text, font=font, fill=fill)
    tile = tile.rotate(-rot_deg, expand=True, resample=Image.Resampling.BICUBIC)
    base_img.alpha_composite(
        tile, (int(round(x - tile.width / 2.0)),
               int(round(y - tile.height / 2.0))))


def render_plasmid_map_png(feats: list, total: int, *, title: str = "",
                           sites: "list | None" = None,
                           size: int = _MAP_IMAGE_DEFAULT_SIZE,
                           origin_bp: int = 0, transparent: bool = False,
                           show_labels: bool = True,
                           show_sites: bool = True,
                           supersample: int = _PNG_SUPERSAMPLE) -> bytes:
    """Render the plasmid map as PNG bytes via Pillow. Drawn at
    ``supersample``× then LANCZOS-downscaled so arcs/lines are anti-aliased."""
    from io import BytesIO
    from PIL import Image, ImageDraw
    size = _clamp_size(size)
    ss = max(1, min(4, int(supersample)))
    big = size * ss
    prims = _build_primitives(
        feats, total, title=title, sites=sites, size=big,
        origin_bp=origin_bp, show_labels=show_labels, show_sites=show_sites,
    )
    bg = (255, 255, 255, 0) if transparent else (255, 255, 255, 255)
    img = Image.new("RGBA", (big, big), bg)
    draw = ImageDraw.Draw(img)

    for cx, cy, r, stroke, w, _fill in prims["circles"]:
        draw.ellipse([cx - r, cy - r, cx + r, cy + r],
                     outline=_to_rgb(stroke), width=max(1, int(round(w))))
    for poly, fill, stroke, w in prims["polys"]:
        draw.polygon(poly, fill=_to_rgb(fill),
                     outline=(0, 0, 0), width=max(1, int(round(w))))
    for x0, y0, x1, y1, color, w in prims["lines"]:
        draw.line([x0, y0, x1, y1], fill=_to_rgb(color),
                  width=max(1, int(round(w))))
    _font_cache: dict = {}
    for x, y, text, fs, color, anchor, weight, rot in prims["texts"]:
        key = int(fs)
        font = _font_cache.get(key)
        if font is None:
            font = _load_font(key)
            _font_cache[key] = font
        rgb = _to_rgb(color)
        pa = _PIL_ANCHOR.get(anchor, "lm")
        if rot:
            _draw_rotated_text(img, x, y, text, font, rgb, rot)
            continue
        try:
            draw.text((x, y), text, font=font, fill=rgb, anchor=pa)
            if weight == "bold":   # fake-bold: 1px over-stamp at supersample
                draw.text((x + ss, y), text, font=font, fill=rgb, anchor=pa)
        except (ValueError, OSError):
            # Anchored draw can reject certain default-font/anchor combos on
            # old Pillow — fall back to a top-left placement.
            draw.text((x, y), text, font=font, fill=rgb)

    if ss > 1:
        img = img.resize((size, size), Image.Resampling.LANCZOS)
    if not transparent:
        img = img.convert("RGB")
    buf = BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


# ── Record → data extractors (for off-screen / bulk export) ──────────────────

def _feat_bounds(feat, total: int, *,
                 circular: bool = True) -> "tuple[int, int] | None":
    """Wrap-aware (start, end) for a Biopython SeqFeature, or None if it
    can't be resolved. ``end`` may be < ``start`` (origin-crossing).

    The origin-wrap test matches `splicecraft_util._feat_bounds`: a TWO-part
    location covering ``[0, x)`` and ``[y, total)`` is the wrap, whatever
    order the parts are stored in. Reading `parts[0].start` / `parts[-1].end`
    positionally got the plus strand right by luck (INSDC writes it
    tail-first) and the MINUS strand wrong — Biopython reverses the parts of
    a ``complement(join(...))``, so a minus-strand wrap came back as
    ``(0, total)`` and drew as a full circle instead of its real arc.

    ``circular=False`` skips the re-encoding: a linear molecule has no origin,
    so that shape is a spliced feature there and flattens to outer bounds."""
    try:
        # One reading geometry for the whole app: `_feature_traversal` owns
        # the wrap-vs-splice decision (declared part order = reading order),
        # so the exported image can't disagree with the on-screen map.
        trav = _feature_traversal(getattr(feat, "location", None), total,
                                  circular=circular)
        if trav is None:
            return None
        start, end = int(trav[0]), int(trav[1])
        if end < start:
            return start, end
    except (TypeError, ValueError, AttributeError, IndexError):
        return None
    if total > 0:
        start %= total
        # Don't fold a full-length feature's end (== total) down to 0 — that
        # would zero its arc length and drop it. Only wrap ends past the origin.
        if end > total:
            end %= total
    return start, end


def _map_feats_from_record(record) -> "tuple[list[dict], int]":
    """Re-derive the map feature list (same shape as ``PlasmidMap._feats``)
    from a bare ``SeqRecord`` — for exporting entries that aren't the
    currently-loaded plasmid. Colours honour a feature's stored
    ``ApEinfo_fwdcolor`` / ``color`` qualifier, else cycle the palette."""
    if record is None:
        return [], 0
    total = len(getattr(record, "seq", "") or "")
    circular = str(
        (getattr(record, "annotations", {}) or {}).get("topology", "")
    ).strip().lower() != "linear"
    feats: list[dict] = []
    skip = {"source"}
    idx = 0
    for feat in getattr(record, "features", []) or []:
        if getattr(feat, "type", "") in skip:
            continue
        bounds = _feat_bounds(feat, total, circular=circular)
        if bounds is None:
            continue
        start, end = bounds
        strand = getattr(getattr(feat, "location", None), "strand", None)
        color = ""
        for q in ("ApEinfo_fwdcolor", "ApEinfo_revcolor", "color"):
            vals = feat.qualifiers.get(q, []) if hasattr(feat, "qualifiers") else []
            if vals:
                cand = str(vals[0]).strip()
                if cand.startswith("#") and len(cand) in (4, 7):
                    color = cand
                    break
        feats.append({
            "type": getattr(feat, "type", "misc"),
            "start": start, "end": end,
            # 0, not 1: Biopython's `strand is None` means "no direction",
            # and the on-screen path preserves it as 0 (arrowless). Coercing
            # to 1 grew a direction arrowhead on exported PNG/SVG maps that
            # the screen does not show.
            "strand": 0 if strand is None else int(strand),
            "color": color or _FALLBACK_PALETTE[idx % len(_FALLBACK_PALETTE)],
            "label": _feat_label_full(feat),
        })
        idx += 1
    return feats, total


def _map_sites_from_record(record, *, allowed_enzymes=None,
                           unique_only: bool = True) -> "list[dict]":
    """Restriction cut sites (recut dicts: ``{start, label}``) for a record,
    for the outer tick ring. Defensive: an unbuilt scan catalog or a bad
    sequence yields ``[]`` rather than raising."""
    seq = str(getattr(record, "seq", "") or "")
    if not seq:
        return []
    try:
        hits = _bio._scan_restriction_sites(
            seq, unique_only=unique_only, circular=True,
            allowed_enzymes=allowed_enzymes,
        )
    except Exception:      # scan catalog may be unbuilt off-hub; never fail export
        _log.debug("map image: restriction scan unavailable", exc_info=True)
        return []
    return [h for h in hits if h.get("type") == "recut"]


# ── High-level file export ───────────────────────────────────────────────────

def _clamp_size(size: int) -> int:
    try:
        size = int(size)
    except (TypeError, ValueError):
        size = _MAP_IMAGE_DEFAULT_SIZE
    return max(_MAP_IMAGE_MIN_SIZE, min(_MAP_IMAGE_MAX_SIZE, size))


def render_map_bytes(feats: list, total: int, *, fmt: str, **kw) -> bytes:
    """Render to encoded bytes in ``fmt`` ('png' or 'svg')."""
    fmt = (fmt or "png").lower()
    if fmt == "svg":
        return render_plasmid_map_svg(feats, total, **kw).encode("utf-8")
    if fmt == "png":
        return render_plasmid_map_png(feats, total, **kw)
    raise ValueError(f"unknown map image format {fmt!r}; "
                     f"choose one of {sorted(_MAP_EXPORT_FORMATS)}")


def export_plasmid_map(out_path, *, feats: "list | None" = None,
                       total: "int | None" = None, record=None,
                       fmt: str = "png", title: str = "",
                       sites: "list | None" = None,
                       size: int = _MAP_IMAGE_DEFAULT_SIZE,
                       origin_bp: int = 0, transparent: bool = False,
                       show_labels: bool = True,
                       show_sites: bool = True) -> dict:
    """Render a plasmid map and atomically write it to ``out_path``.

    Supply EITHER an already-parsed (``feats``, ``total``) pair (the live
    map hands us ``pm._feats`` / ``pm._total`` for on-screen fidelity) OR a
    ``record`` (re-derived via :func:`_map_feats_from_record`). Returns a
    summary ``{path, fmt, bp, features, bytes}``.

    Writes go to a user-chosen path OUTSIDE the sacred data dir (exactly
    like the GenBank/FASTA exporters), so no chokepoint applies — but the
    write is still atomic (temp file + fsync + replace) to avoid a
    half-written image on a crash/full-disk."""
    fmt = (fmt or "png").lower()
    if fmt not in _MAP_EXPORT_FORMATS:
        raise ValueError(f"unknown map image format {fmt!r}; "
                         f"choose one of {sorted(_MAP_EXPORT_FORMATS)}")
    if feats is None or total is None:
        if record is None:
            raise ValueError("export_plasmid_map: pass feats+total or record")
        feats, total = _map_feats_from_record(record)
        if show_sites and sites is None:
            sites = _map_sites_from_record(record)
    feats = feats or []
    total = int(total or 0)

    data = render_map_bytes(
        feats, total, fmt=fmt, title=title, sites=sites, size=size,
        origin_bp=origin_bp, transparent=transparent,
        show_labels=show_labels, show_sites=show_sites,
    )
    out = Path(out_path).expanduser()
    _atomic_write_bytes(out, data)
    _log.info("Exported plasmid map → %s (%s, %d bp, %d feats, %d bytes)",
              out, fmt, total, len(feats), len(data))
    return {"path": str(out), "fmt": fmt, "bp": total,
            "features": len(feats), "bytes": len(data)}


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    """Temp-file + fsync + os.replace write. Mirrors the persistence
    engine's atomicity for a user-chosen export path (which is outside the
    data dir, so it needs no chokepoint authorisation)."""
    import os
    import tempfile
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".sc-map-", suffix=path.suffix,
                               dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
