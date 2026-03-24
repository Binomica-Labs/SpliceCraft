#!/usr/bin/env python3
"""
splicecraft.py
==============
SpliceCraft — terminal circular plasmid map viewer.

Features:
  - Fetch any GenBank record by accession (pUC19 = L09137)
  - Load local .gb / .gbk files
  - Circular map with per-strand feature rings and arrowheads
  - Rotate origin freely with ← → keys or mouse scroll
  - Click map to select feature; click sidebar row to highlight on map
  - Feature detail panel
  - Plasmid library panel (left, SnapGene-style collection, persistent JSON)
  - DNA sequence viewer / editor (bottom, press e to edit, Ctrl+S to save)
  - CDS amino acid translation shown on feature selection

Run standalone:
    python3 splicecraft.py
    python3 splicecraft.py L09137          # fetch pUC19 on launch
    python3 splicecraft.py myplasmid.gb    # open local file
"""

import json
import math
import sys
from io import StringIO
from pathlib import Path

# ── Dependency check ───────────────────────────────────────────────────────────

_REQUIRED = {
    "textual":  "textual",
    "Bio":      "biopython",
}

def _check_deps():
    missing = []
    for module, package in _REQUIRED.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(package)
    if missing:
        print("Missing dependencies — install with:")
        print(f"  pip install {' '.join(missing)}")
        sys.exit(1)

_check_deps()

from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, ScrollableContainer
from textual.coordinate import Coordinate
from textual.events import Click, MouseDown, MouseMove, MouseUp, MouseScrollDown, MouseScrollUp
from textual.message import Message
from textual.reactive import reactive
from textual.screen import ModalScreen
from textual.widget import Widget
from textual.widgets import (
    Button, DataTable, Footer, Header, Input, Label,
    Select, Static, TabbedContent, TabPane,
)
from rich.text import Text

# ── Feature appearance ─────────────────────────────────────────────────────────

# Visually distinct per-feature palette (assigned by feature index, not type)
_FEATURE_PALETTE: list[str] = [
    "color(39)",   "color(118)",  "color(208)",  "color(213)",  "color(51)",
    "color(220)",  "color(196)",  "color(46)",   "color(201)",  "color(129)",
    "color(166)",  "color(33)",   "color(226)",  "color(160)",  "color(87)",
    "color(105)",  "color(154)",  "color(203)",  "color(81)",   "color(185)",
]

# ── Library persistence ────────────────────────────────────────────────────────

_LIBRARY_FILE = Path(__file__).parent / "plasmid_library.json"

def _load_library() -> list[dict]:
    if _LIBRARY_FILE.exists():
        try:
            return json.loads(_LIBRARY_FILE.read_text())
        except Exception:
            pass
    return []

def _save_library(entries: list[dict]) -> None:
    try:
        _LIBRARY_FILE.write_text(json.dumps(entries, indent=2))
    except Exception:
        pass

# ── Restriction sites ──────────────────────────────────────────────────────────

_RESTRICTION_SITES: dict[str, str] = {
    "EcoRI":   "GAATTC",
    "BamHI":   "GGATCC",
    "HindIII": "AAGCTT",
    "NcoI":    "CCATGG",
    "NdeI":    "CATATG",
    "XhoI":    "CTCGAG",
    "SalI":    "GTCGAC",
    "KpnI":    "GGTACC",
    "SacI":    "GAGCTC",
    "SpeI":    "ACTAGT",
    "XbaI":    "TCTAGA",
    "BsaI":    "GGTCTC",
    "BsmBI":   "CGTCTC",
    "BbsI":    "GAAGAC",
    "NotI":    "GCGGCCGC",
    "PstI":    "CTGCAG",
    "SphI":    "GCATGC",
    "ClaI":    "ATCGAT",
}

_RESTR_PALETTE: list[str] = [
    "color(220)", "color(208)", "color(196)", "color(160)",
    "color(105)", "color(129)", "color(57)",  "color(21)",
    "color(33)",  "color(39)",  "color(51)",  "color(87)",
    "color(118)", "color(154)", "color(190)", "color(226)",
    "color(185)", "color(180)",
]
_RESTR_COLOR: dict[str, str] = {
    name: _RESTR_PALETTE[i % len(_RESTR_PALETTE)]
    for i, name in enumerate(_RESTRICTION_SITES)
}


def _scan_restriction_sites(seq: str) -> list[dict]:
    """Scan both strands; return feature-like dicts for every hit."""
    seq_u = seq.upper()
    comp  = str.maketrans("ACGT", "TGCA")
    seq_rc = seq_u.translate(comp)[::-1]
    n = len(seq_u)
    feats: list[dict] = []
    for name, site in _RESTRICTION_SITES.items():
        color = _RESTR_COLOR[name]
        pos = 0
        while (p := seq_u.find(site, pos)) != -1:
            feats.append({"type": "site", "start": p, "end": p + len(site),
                          "strand": 1,  "color": color, "label": name})
            pos = p + 1
        pos = 0
        while (p := seq_rc.find(site, pos)) != -1:
            start = n - p - len(site)
            feats.append({"type": "site", "start": start, "end": start + len(site),
                          "strand": -1, "color": color, "label": name})
            pos = p + 1
    return feats


def _build_seq_text(seq: str, feats: list[dict], line_width: int = 60,
                    sel_range: "tuple[int,int] | None" = None,
                    user_sel:  "tuple[int,int] | None" = None,
                    cursor_pos: int = -1) -> Text:
    """Rich Text of the sequence with per-position feature coloring.

    sel_range  — feature highlight: bold + underline on feature bases
    user_sel   — shift-click selection: subtle background, used by edit dialog
    cursor_pos — click cursor: │ inserted before cursor_pos
    Annotation bars appear below each DNA line (one bar per overlapping
    non-site feature, capped at 4, largest first).
    """
    n = len(seq)
    styles = ["color(252)"] * n
    for f in reversed(feats):          # reversed so first feature wins
        col = f["color"]
        for i in range(f["start"], min(f["end"], n)):
            styles[i] = col

    sel_s  = sel_range[0] if sel_range else -1
    sel_e  = sel_range[1] if sel_range else -1
    usr_s  = user_sel[0]  if user_sel  else -1
    usr_e  = user_sel[1]  if user_sel  else -1

    seq_upper = seq.upper()
    result    = Text(no_wrap=False)

    # Annotation-bar features: exclude restriction sites, largest first
    annot_feats = sorted(
        [f for f in feats if f.get("type") != "site"],
        key=lambda f: -(f["end"] - f["start"]),
    )

    for chunk_start in range(0, n, line_width):
        chunk_end = min(chunk_start + line_width, n)
        result.append(f"{chunk_start + 1:>8}  ", style="color(245)")

        # RLE sequence line — priority: user_sel > feat_sel > plain
        # Cursor is a │ inserted *before* cursor_pos
        run_chars: list[str] = []
        run_style = ""
        for i in range(chunk_start, chunk_end):
            # Insert │ cursor before this position
            if cursor_pos == i:
                if run_chars:
                    result.append("".join(run_chars), style=run_style)
                    run_chars = []
                    run_style = ""
                result.append("│", style="bold white")
            base   = styles[i]
            in_usr = (usr_s <= i < usr_e)
            in_sel = (sel_s <= i < sel_e)
            if in_usr:
                sty = base + " on color(237)"
            elif in_sel:
                sty = "bold underline " + base
            else:
                sty = base
            if sty == run_style:
                run_chars.append(seq_upper[i])
            else:
                if run_chars:
                    result.append("".join(run_chars), style=run_style)
                run_style = sty
                run_chars = [seq_upper[i]]
        if run_chars:
            result.append("".join(run_chars), style=run_style)
        # Cursor at end of this line (after last char of chunk)
        if chunk_start <= cursor_pos == chunk_end:
            result.append("│", style="bold white")
        result.append("\n")

        # Annotation bars (up to 4 features overlapping this chunk)
        chunk_feats = [
            f for f in annot_feats
            if f["start"] < chunk_end and f["end"] > chunk_start
        ]
        for f in chunk_feats[:4]:
            bar_s = max(f["start"], chunk_start) - chunk_start
            bar_e = min(f["end"],   chunk_end)   - chunk_start
            bar_len = bar_e - bar_s
            if bar_len <= 0:
                continue
            starts_here = f["start"] >= chunk_start
            ends_here   = f["end"]   <= chunk_end
            strand      = f["strand"]
            if strand >= 0:
                lc = "━" if starts_here else "╌"
                rc = "▶" if ends_here   else "╌"
            else:
                lc = "◀" if starts_here else "╌"
                rc = "━" if ends_here   else "╌"
            if bar_len == 1:
                bar = ("▶" if ends_here   else "━") if strand >= 0 else \
                      ("◀" if starts_here else "━")
            elif bar_len == 2:
                bar = lc + rc
            else:
                inner = bar_len - 2
                lbl   = f["label"][:inner]
                pad   = inner - len(lbl)
                pl    = pad // 2
                bar   = lc + "━" * pl + lbl + "━" * (pad - pl) + rc
            result.append(" " * (10 + bar_s))
            result.append(bar, style=f["color"])
            result.append("\n")

    return result


# Standard genetic code for CDS translation (no biopython dependency)
_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

def _copy_to_clipboard(text: str) -> bool:
    """Copy text to system clipboard via xclip / xsel / wl-copy. Returns True on success."""
    import subprocess
    for cmd in (
        ["xclip", "-selection", "clipboard"],
        ["xsel", "--clipboard", "--input"],
        ["wl-copy"],
    ):
        try:
            subprocess.run(cmd, input=text.encode(), check=True, timeout=3,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return True
        except Exception:
            continue
    try:
        import pyperclip
        pyperclip.copy(text)
        return True
    except Exception:
        return False


def _translate_cds(full_seq: str, start: int, end: int, strand: int) -> str:
    """Translate a CDS region to single-letter AA string (stop codon → *)."""
    sub = full_seq[start:end].upper()
    if strand == -1:
        sub = sub.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    aa = [_CODON_TABLE.get(sub[i:i+3], "?") for i in range(0, len(sub) - 2, 3)]
    result = "".join(aa)
    if result and not result.endswith("*"):
        result += "*"
    return result


# ── Char aspect detection ──────────────────────────────────────────────────────

def _detect_char_aspect() -> float:
    """
    Return char_height / char_width for the current terminal by reading
    pixel dimensions via TIOCGWINSZ.  Falls back to 2.0 if unavailable.
    """
    try:
        import fcntl, os, struct, termios
        fds: list[tuple[int, bool]] = []
        try:
            fds.append((os.open("/dev/tty", os.O_RDWR | os.O_NOCTTY), True))
        except OSError:
            pass
        fds += [(0, False), (1, False), (2, False)]
        try:
            for fd, _ in fds:
                try:
                    buf = struct.pack("HHHH", 0, 0, 0, 0)
                    res = fcntl.ioctl(fd, termios.TIOCGWINSZ, buf)
                    rows, cols, px_w, px_h = struct.unpack("HHHH", res)
                    if px_w > 0 and px_h > 0 and rows > 0 and cols > 0:
                        ratio = (px_h / rows) / (px_w / cols)
                        if 0.8 <= ratio <= 5.0:
                            return round(ratio, 3)
                except OSError:
                    continue
        finally:
            for fd, should_close in fds:
                if should_close:
                    try:
                        os.close(fd)
                    except OSError:
                        pass
    except Exception:
        pass
    return 2.0


def _format_bp(bp: int) -> str:
    """Human-readable bp label with up to 2 decimal places in 'k' units."""
    if bp < 1000:
        return str(bp)
    if bp % 1000 == 0:
        return f"{bp // 1000}k"
    if bp % 100 == 0:
        return f"{bp / 1000:.1f}k"
    return f"{bp / 1000:.2f}k"


def _feat_label(feat) -> str:
    for q in ("label", "gene", "product", "standard_name", "note", "bound_moiety"):
        if q in feat.qualifiers:
            v = feat.qualifiers[q]
            s = v[0] if isinstance(v, list) else v
            return s[:28]
    return feat.type

def _nice_tick(total: int) -> int:
    """A tick interval that gives ~6-10 ticks for this plasmid size."""
    for t in [50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 25000, 50000]:
        if 4 <= total // t <= 14:
            return t
    return max(1, total // 8)

# ── GenBank I/O ────────────────────────────────────────────────────────────────

def fetch_genbank(accession: str, email: str = "splicecraft@local"):
    """Fetch a GenBank record by accession from NCBI Entrez. Returns SeqRecord."""
    from Bio import Entrez, SeqIO
    Entrez.email = email
    with Entrez.efetch(
        db="nucleotide", id=accession, rettype="gb", retmode="text"
    ) as handle:
        record = SeqIO.read(handle, "genbank")
    return record

def load_genbank(path: str):
    """Load a GenBank (.gb/.gbk) file. Returns SeqRecord."""
    from Bio import SeqIO
    return SeqIO.read(path, "genbank")

def _record_to_gb_text(record) -> str:
    """Serialize a SeqRecord to GenBank format text."""
    from Bio import SeqIO
    buf = StringIO()
    SeqIO.write(record, buf, "genbank")
    return buf.getvalue()

def _gb_text_to_record(text: str):
    """Parse GenBank format text back to a SeqRecord."""
    from Bio import SeqIO
    return SeqIO.read(StringIO(text), "genbank")

# ── Core drawing ───────────────────────────────────────────────────────────────

class _Canvas:
    """A mutable 2-D character grid that renders to a Rich Text object."""

    def __init__(self, width: int, height: int):
        self.w = width
        self.h = height
        self._chars:  list[list[str]] = [[" "] * width for _ in range(height)]
        self._styles: list[list[str]] = [[""]  * width for _ in range(height)]

    def put(self, col: int, row: int, ch: str, style: str = ""):
        if 0 <= col < self.w and 0 <= row < self.h:
            self._chars[row][col]  = ch
            self._styles[row][col] = style

    def put_text(self, col: int, row: int, text: str, style: str = ""):
        for j, ch in enumerate(text):
            self.put(col + j, row, ch, style)

    def to_rich_text(self) -> Text:
        result = Text(no_wrap=True, overflow="crop")
        for r in range(self.h):
            for c in range(self.w):
                ch = self._chars[r][c]
                st = self._styles[r][c]
                result.append(ch, style=st) if st else result.append(ch)
            if r < self.h - 1:
                result.append("\n")
        return result


class _BrailleCanvas:
    """
    Sub-character resolution canvas using Unicode braille (U+2800–U+28FF).

    Each terminal cell (col, row) encodes a 2-wide × 4-tall dot grid —
    8 pixels per character cell.  Braille dot layout:

        px%2=0  px%2=1
        dot1    dot4    ← py%4=0   (bits 0, 3)
        dot2    dot5    ← py%4=1   (bits 1, 4)
        dot3    dot6    ← py%4=2   (bits 2, 5)
        dot7    dot8    ← py%4=3   (bits 6, 7)

    Codepoint = 0x2800 + bitmask of active dots.
    Colors: higher-priority write wins per cell.
    """

    _DOT_BITS: list[list[int]] = [
        [0, 3],
        [1, 4],
        [2, 5],
        [6, 7],
    ]

    def __init__(self, cols: int, rows: int):
        self.cols = cols
        self.rows = rows
        self._bits:   list[list[int]] = [[0]  * cols for _ in range(rows)]
        self._colors: list[list[str]] = [[" "] * cols for _ in range(rows)]
        self._prio:   list[list[int]] = [[0]  * cols for _ in range(rows)]

    def set_pixel(self, px: int, py: int,
                  color: str = "", priority: int = 1) -> None:
        col, row = px // 2, py // 4
        if not (0 <= col < self.cols and 0 <= row < self.rows):
            return
        self._bits[row][col] |= 1 << self._DOT_BITS[py % 4][px % 2]
        if color and priority >= self._prio[row][col]:
            self._colors[row][col] = color
            self._prio[row][col]   = priority

    def combine(self, text_canvas: "_Canvas") -> Text:
        """
        Return a Rich Text object.
        Non-space cells from *text_canvas* are drawn on top;
        braille pixels fill the rest.
        Consecutive blank cells are batched into a single append call.
        """
        result = Text(no_wrap=True, overflow="crop")
        rows   = min(self.rows, text_canvas.h)
        cols   = min(self.cols, text_canvas.w)
        tc_chars  = text_canvas._chars
        tc_styles = text_canvas._styles
        bc_bits   = self._bits
        bc_colors = self._colors
        for row in range(rows):
            blank_run = 0
            tc_row  = tc_chars[row]
            tcs_row = tc_styles[row]
            bc_bits_row   = bc_bits[row]
            bc_colors_row = bc_colors[row]
            for col in range(cols):
                tc_ch = tc_row[col]
                if tc_ch == " " and not bc_bits_row[col]:
                    blank_run += 1
                else:
                    if blank_run:
                        result.append(" " * blank_run)
                        blank_run = 0
                    if tc_ch != " ":
                        st = tcs_row[col]
                        result.append(tc_ch, style=st) if st else result.append(tc_ch)
                    else:
                        ch = chr(0x2800 + bc_bits_row[col])
                        c  = bc_colors_row[col]
                        result.append(ch, style=c) if c != " " else result.append(ch)
            if blank_run:
                result.append(" " * blank_run)
            if row < rows - 1:
                result.append("\n")
        return result


def _arrow_char(tangent_angle: float) -> str:
    """Pick a directional arrow char for an arrowhead."""
    t = tangent_angle % (2 * math.pi)
    sector = int((t + math.pi / 4) / (math.pi / 2)) % 4
    return ["▶", "▼", "◀", "▲"][sector]


# ── PlasmidMap widget ──────────────────────────────────────────────────────────

class PlasmidMap(Widget):
    """
    Circular plasmid map.

    Keyboard: ← → to rotate when focused (click map first); [ ] always rotate.
    Mouse:    scroll to rotate; click to select a feature.
    """

    DEFAULT_CSS = """
    PlasmidMap {
        width: 1fr;
        height: 1fr;
        background: $background;
    }
    PlasmidMap:focus { border: solid $accent; }
    """

    can_focus = True

    BINDINGS = [
        Binding("left",        "rotate_cw",      "Rotate ←",     show=True),
        Binding("right",       "rotate_ccw",      "Rotate →",     show=True),
        Binding("shift+left",  "rotate_cw_lg",   "Rotate ←←",    show=False),
        Binding("shift+right", "rotate_ccw_lg",  "Rotate →→",    show=False),
        Binding("home",        "reset_origin",   "Reset",        show=False),
        Binding("comma",       "aspect_dec",     "Circle wider",  show=False),
        Binding("full_stop",   "aspect_inc",     "Circle taller", show=False),
    ]

    origin_bp:    reactive[int]   = reactive(0)
    selected_idx: reactive[int]   = reactive(-1)
    _aspect:      reactive[float] = reactive(2.0)

    # ── Messages ───────────────────────────────────────────────────────────────

    class FeatureSelected(Message):
        def __init__(self, idx: int, feat_dict: dict | None):
            self.idx       = idx
            self.feat_dict = feat_dict
            super().__init__()

    # ── Lifecycle ──────────────────────────────────────────────────────────────

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.record  = None
        self._feats:       list[dict] = []
        self._restr_feats: list[dict] = []   # restriction site overlay
        self._total: int = 0

    def on_mount(self) -> None:
        detected = _detect_char_aspect()
        if detected != self._aspect:
            self._aspect = detected

    def load_record(self, record) -> None:
        self.record       = record
        self._total       = len(record.seq)
        self.origin_bp    = 0
        self.selected_idx = -1
        self._feats       = self._parse(record)
        self._restr_feats = []
        self.refresh()

    def _parse(self, record) -> list[dict]:
        feats = []
        for feat in record.features:
            if feat.type in ("source",):
                continue
            start  = int(feat.location.start)
            end    = int(feat.location.end)
            strand = getattr(feat.location, "strand", 1) or 1
            idx    = len(feats)
            feats.append({
                "type":   feat.type,
                "start":  start,
                "end":    end,
                "strand": strand,
                "color":  _FEATURE_PALETTE[idx % len(_FEATURE_PALETTE)],
                "label":  _feat_label(feat),
            })
        return feats

    # ── Actions ────────────────────────────────────────────────────────────────

    def _step(self, coarse: bool = False) -> int:
        if not self._total:
            return 1
        return max(1, self._total // (10 if coarse else 200))

    def action_rotate_cw(self):
        if self._total:
            self.origin_bp = (self.origin_bp - self._step()) % self._total

    def action_rotate_ccw(self):
        if self._total:
            self.origin_bp = (self.origin_bp + self._step()) % self._total

    def action_rotate_cw_lg(self):
        if self._total:
            self.origin_bp = (self.origin_bp - self._step(True)) % self._total

    def action_rotate_ccw_lg(self):
        if self._total:
            self.origin_bp = (self.origin_bp + self._step(True)) % self._total

    def action_reset_origin(self):
        self.origin_bp = 0

    def action_aspect_inc(self):
        self._aspect = round(min(5.0, self._aspect + 0.05), 3)
        self.notify(f"Aspect {self._aspect:.2f}  (press , to widen)", timeout=1.5)

    def action_aspect_dec(self):
        self._aspect = round(max(0.5, self._aspect - 0.05), 3)
        self.notify(f"Aspect {self._aspect:.2f}  (press . to heighten)", timeout=1.5)

    def select_feature(self, idx: int) -> None:
        self.selected_idx = idx
        self.refresh()

    # ── Mouse ──────────────────────────────────────────────────────────────────

    def on_mouse_scroll_up(self, _: MouseScrollUp):
        self.action_rotate_cw()

    def on_mouse_scroll_down(self, _: MouseScrollDown):
        self.action_rotate_ccw()

    def _feat_at(self, x: int, y: int) -> int:
        if not self.record or not self._total:
            return -1
        w, h = self.size.width, self.size.height
        cx, cy, rx, ry = self._geometry(w, h)
        if rx == 0 or ry == 0:
            return -1
        dc_n = (x - cx) / rx
        dr_n = (y - cy) / ry
        r_norm = math.sqrt(dc_n ** 2 + dr_n ** 2)
        if r_norm < 0.75 or r_norm > 1.35:
            return -1
        angle = math.atan2(dr_n, dc_n)
        angle_norm = (angle + math.pi / 2) % (2 * math.pi)
        bp = int(self.origin_bp + self._total * angle_norm / (2 * math.pi)) % self._total
        for i, f in enumerate(self._feats):
            if self._bp_in(bp, f):
                return i
        return -1

    def on_click(self, event: Click):
        if not self.record:
            return
        idx = self._feat_at(event.x, event.y)
        self.selected_idx = idx
        f = self._feats[idx] if idx >= 0 else None
        self.post_message(self.FeatureSelected(idx, f))

    def _bp_in(self, bp: int, f: dict) -> bool:
        s, e = f["start"], f["end"]
        return (s <= bp < e) if e >= s else (bp >= s or bp < e)

    # ── Geometry ───────────────────────────────────────────────────────────────

    def _geometry(self, w: int, h: int) -> tuple[int, int, int, int]:
        cx     = w // 2
        cy     = h // 2
        aspect = self._aspect
        rx_from_w = cx - 16
        rx_from_h = int((cy - 3) * aspect)
        rx = max(8, min(rx_from_w, rx_from_h))
        ry = max(4, round(rx / aspect))
        return cx, cy, rx, ry

    def _angle_to_xy(self, angle: float, cx: int, cy: int,
                     rx: int, ry: int, dr: int = 0) -> tuple[int, int]:
        scale = ry / rx if rx else 0.5
        x = round(cx + (rx + dr) * math.cos(angle))
        y = round(cy + (ry + dr * scale) * math.sin(angle))
        return x, y

    def _bp_to_angle(self, bp: int) -> float:
        return 2 * math.pi * ((bp - self.origin_bp) % self._total) / self._total - math.pi / 2

    # ── Render ─────────────────────────────────────────────────────────────────

    # cache: (w, h, origin_bp, selected_idx, _aspect, n_feats, n_restr) → Text
    _render_cache: "tuple | None" = None

    def render(self) -> Text:
        if not self.record:
            return Text(
                "\n\n   No record loaded.\n"
                "   Press  f  to fetch from GenBank (e.g. L09137 for pUC19)\n"
                "   Press  o  to open a local .gb file\n",
                style="dim italic",
            )
        w, h = self.size.width, self.size.height
        if w < 30 or h < 14:
            return Text(f"  Window too small ({w}×{h})", style="dim red")
        key = (w, h, self.origin_bp, self.selected_idx, self._aspect,
               len(self._feats), len(self._restr_feats))
        if self._render_cache and self._render_cache[0] == key:
            return self._render_cache[1]
        result = self._draw(w, h)
        self._render_cache = (key, result)
        return result

    def _draw(self, w: int, h: int) -> Text:
        canvas = _Canvas(w, h)
        bc     = _BrailleCanvas(w, h)

        total  = self._total
        cx, cy, rx, ry = self._geometry(w, h)

        _ry_rx   = ry / rx if rx else 0.5
        TWO_PI   = 2 * math.pi
        HALF_PI  = math.pi / 2
        # Braille-space offsets for the circle centre
        cx2, cy4 = cx * 2, cy * 4
        rx2, ry4 = rx * 2, ry * 4

        # Enough steps to fill every braille column (2 dots/cell → circumference
        # in braille-px ≈ 2π*rx*2; ×2 safety margin → 4π*rx ≈ 13*rx)
        n_steps = max(360, int(13 * rx))

        # Local refs avoid repeated attribute lookups in hot loops
        bc_bits   = bc._bits
        bc_colors = bc._colors
        bc_prio   = bc._prio
        bc_cols   = bc.cols
        bc_rows   = bc.rows
        # Flat braille dot-bit lookup: index = (py & 3) << 1 | (px & 1)
        _BIT = (0, 3, 1, 4, 2, 5, 6, 7)

        def bp2a(bp: int) -> float:
            return self._bp_to_angle(bp)

        def a2xy(angle: float, dr: int = 0) -> tuple[int, int]:
            return self._angle_to_xy(angle, cx, cy, rx, ry, dr)

        # ── Backbone (inlined — no function call overhead) ─────────────────────
        _step_angle = TWO_PI / n_steps
        for step in range(n_steps):
            a  = step * _step_angle - HALF_PI
            px = round(cx2 + rx2 * math.cos(a))
            py = round(cy4 + ry4 * math.sin(a))
            pcol, prow = px >> 1, py >> 2
            if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                bc_bits[prow][pcol]   |= 1 << _BIT[(py & 3) << 1 | (px & 1)]
                if bc_prio[prow][pcol] == 0:
                    bc_colors[prow][pcol] = "color(238)"
                    bc_prio[prow][pcol]   = 1

        # ── Position ticks ────────────────────────────────────────────────────
        tick_int = _nice_tick(total)
        bp = 0
        while bp < total:
            angle = bp2a(bp)
            tx, ty = a2xy(angle, dr=2)
            canvas.put(tx, ty, "┼", "color(250)")
            label = _format_bp(bp)
            lx, ly = a2xy(angle, dr=4)
            if math.cos(angle) >= 0:
                canvas.put_text(lx, ly, label, "color(245)")
            else:
                canvas.put_text(lx - len(label) + 1, ly, label, "color(245)")
            bp += tick_int

        # ── Restriction site marks ─────────────────────────────────────────────
        for rf in self._restr_feats:
            mid_bp = (rf["start"] + rf["end"]) // 2
            tx, ty = a2xy(bp2a(mid_bp), dr=5)
            canvas.put(tx, ty, "▪", rf["color"])

        # ── Features (large → small) ──────────────────────────────────────────
        feats_sorted = sorted(
            enumerate(self._feats),
            key=lambda iv: -(iv[1]["end"] - iv[1]["start"]),
        )
        label_slots: list[tuple[float, str, str]] = []
        N_DR = 4   # radial samples per arc (was 6; 4 is visually identical)

        for orig_idx, f in feats_sorted:
            is_sel   = (orig_idx == self.selected_idx)
            color    = f["color"]
            style    = ("reverse " + color) if is_sel else color
            f_prio   = 3 if is_sel else 2

            start_bp, end_bp = f["start"], f["end"]
            strand           = f["strand"]
            dr_lo = 1 if strand >= 0 else -3
            dr_hi = dr_lo + 2

            start_a = bp2a(start_bp)
            end_a   = bp2a(end_bp)
            span_a  = (end_a - start_a) % TWO_PI or TWO_PI

            arc_steps = max(8, int(span_a / TWO_PI * n_steps))
            inv_arc   = 1.0 / arc_steps if arc_steps else 1.0

            # Precompute cos/sin once; reuse across all N_DR radial passes
            arc_cos = [math.cos(start_a + s * inv_arc * span_a) for s in range(arc_steps + 1)]
            arc_sin = [math.sin(start_a + s * inv_arc * span_a) for s in range(arc_steps + 1)]

            inv_ndr1 = 1.0 / max(1, N_DR - 1)
            dr_span  = dr_hi - dr_lo
            for dr_i in range(N_DR):
                dr     = dr_lo + dr_span * dr_i * inv_ndr1
                _rx2dr = (rx + dr) * 2
                _ry4dr = (ry + dr * _ry_rx) * 4
                for ca, sa in zip(arc_cos, arc_sin):
                    px   = round(cx2 + _rx2dr * ca)
                    py   = round(cy4 + _ry4dr * sa)
                    pcol = px >> 1
                    prow = py >> 2
                    if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                        bc_bits[prow][pcol]   |= 1 << _BIT[(py & 3) << 1 | (px & 1)]
                        if f_prio >= bc_prio[prow][pcol]:
                            bc_colors[prow][pcol] = style
                            bc_prio[prow][pcol]   = f_prio

            if strand >= 0:
                tip_a, tip_tan = end_a, end_a + HALF_PI
            else:
                tip_a, tip_tan = start_a, start_a - HALF_PI
            ax, ay = a2xy(tip_a, dr=dr_lo + 1)
            canvas.put(ax, ay, _arrow_char(tip_tan), "bold " + style)

            mid_bp = (start_bp + (end_bp - start_bp) // 2) % total
            label_slots.append((bp2a(mid_bp), f["label"][:16], color))

        # ── Labels ────────────────────────────────────────────────────────────
        label_slots.sort(key=lambda t: t[0])
        last_ly = -99
        for angle, lbl, color in label_slots:
            lx, ly = a2xy(angle, dr=rx // 3 + 5)
            if abs(ly - last_ly) < 2:
                continue
            last_ly = ly
            dot_x, dot_y = a2xy(angle, dr=3)
            canvas.put(dot_x, dot_y, "·", color)
            if math.cos(angle) >= 0:
                canvas.put_text(lx, ly, lbl, color)
            else:
                canvas.put_text(lx - len(lbl) + 1, ly, lbl, color)

        # ── Center info ───────────────────────────────────────────────────────
        name     = (self.record.name or self.record.id or "?")[:w // 3]
        size_txt = f"{total:,} bp"
        orig_txt = f"▲ {self.origin_bp:,}"
        for i, (txt, sty) in enumerate([
            (name,     "bold white"),
            (size_txt, "color(245)"),
            (orig_txt, "dim cyan"),
        ]):
            canvas.put_text(cx - len(txt) // 2, cy - 1 + i, txt, sty)

        return bc.combine(canvas)


# ── Sidebar feature table ──────────────────────────────────────────────────────

class FeatureSidebar(Widget):
    DEFAULT_CSS = """
    FeatureSidebar {
        width: 32;
        border-left: solid $primary;
    }
    #feat-table  { height: 1fr; }
    #detail-box  { height: 8; border-top: solid $accent; padding: 0 1; }
    #sidebar-hdr { background: $primary; padding: 0 1; }
    """

    class RowActivated(Message):
        def __init__(self, idx: int):
            self.idx = idx
            super().__init__()

    def compose(self) -> ComposeResult:
        yield Static(" Features", id="sidebar-hdr")
        yield DataTable(id="feat-table", cursor_type="row", zebra_stripes=True)
        yield Static("", id="detail-box")

    def on_mount(self):
        t = self.query_one("#feat-table", DataTable)
        t.add_columns("Type", "Label", "bp", "±")

    def populate(self, feats: list[dict]) -> None:
        t = self.query_one("#feat-table", DataTable)
        t.clear()
        for f in feats:
            strand_sym = "+" if f["strand"] == 1 else ("−" if f["strand"] == -1 else "·")
            bp_str     = f"{f['start']+1}‥{f['end']}"
            t.add_row(
                Text(f["type"][:12],  style=f["color"]),
                Text(f["label"][:14], style=f["color"]),
                bp_str,
                strand_sym,
            )

    def show_detail(self, f: dict | None) -> None:
        box = self.query_one("#detail-box", Static)
        if f is None:
            box.update(Text(""))
            return
        strand_sym = "+" if f["strand"] == 1 else ("−" if f["strand"] == -1 else "·")
        span = f["end"] - f["start"]
        t = Text()
        t.append(f["type"],  style=f"bold {f['color']}")
        t.append("\n")
        t.append(f["label"], style="white")
        t.append("\n")
        t.append(f"{f['start']+1}‥{f['end']} ({span:,} bp)", style="dim")
        t.append("\n")
        t.append(f"Strand: {strand_sym}", style="dim")
        box.update(t)

    _prog_row: int = -1   # cursor moves driven by highlight_row, not the user

    def highlight_row(self, idx: int) -> None:
        """Move cursor to row; suppresses the resulting RowActivated echo."""
        if idx < 0:
            return
        t = self.query_one("#feat-table", DataTable)
        self._prog_row = idx
        try:
            t.move_cursor(row=idx)
        except Exception:
            self._prog_row = -1

    @on(DataTable.RowHighlighted, "#feat-table")
    def _row_highlighted(self, event: DataTable.RowHighlighted):
        if event.cursor_row == self._prog_row:
            self._prog_row = -1
            return
        self.post_message(self.RowActivated(event.cursor_row))


# ── Library panel ──────────────────────────────────────────────────────────────

class LibraryPanel(Widget):
    """Left-hand plasmid library — persistent SnapGene-style collection."""

    DEFAULT_CSS = """
    LibraryPanel {
        width: 24;
        border-right: solid $primary;
    }
    #lib-hdr   { background: $primary; padding: 0 1; }
    #lib-table { height: 1fr; }
    #lib-btns  { height: 3; }
    #lib-btns Button { min-width: 5; margin: 0 0 0 1; }
    """

    class PlasmidLoad(Message):
        """User selected a library entry to load."""
        def __init__(self, entry: dict):
            self.entry = entry
            super().__init__()

    class AddCurrentRequested(Message):
        """User pressed '+' to add the currently-loaded record."""
        pass

    def compose(self) -> ComposeResult:
        yield Static(" Library", id="lib-hdr")
        yield DataTable(id="lib-table", cursor_type="row", zebra_stripes=True)
        with Horizontal(id="lib-btns"):
            yield Button("+", id="btn-lib-add", variant="success",
                         tooltip="Add current plasmid")
            yield Button("−", id="btn-lib-del", variant="error",
                         tooltip="Remove selected")

    def on_mount(self):
        t = self.query_one("#lib-table", DataTable)
        t.add_columns("Name", "bp")
        self._repopulate()

    def _repopulate(self):
        t = self.query_one("#lib-table", DataTable)
        t.clear()
        for entry in _load_library():
            t.add_row(
                entry["name"][:14],
                f"{entry['size']:,}",
                key=entry["id"],
            )

    def add_entry(self, record) -> None:
        """Serialize record and persist to library JSON."""
        gb_text = _record_to_gb_text(record)
        entries = _load_library()
        entries = [e for e in entries if e.get("id") != record.id]
        import datetime
        entries.insert(0, {
            "name":    record.name or record.id,
            "id":      record.id,
            "size":    len(record.seq),
            "n_feats": len([f for f in record.features if f.type != "source"]),
            "source":  getattr(record, "_tui_source", f"id:{record.id}"),
            "added":   datetime.date.today().isoformat(),
            "gb_text": gb_text,
        })
        _save_library(entries)
        self._repopulate()

    @on(DataTable.RowSelected, "#lib-table")
    def _row_selected(self, event: DataTable.RowSelected):
        key = event.row_key.value if event.row_key else None
        if key is None:
            return
        for entry in _load_library():
            if entry.get("id") == key:
                self.post_message(self.PlasmidLoad(entry))
                return

    @on(Button.Pressed, "#btn-lib-add")
    def _btn_add(self):
        self.post_message(self.AddCurrentRequested())

    def set_dirty(self, dirty: bool) -> None:
        """Show/hide unsaved-changes marker in the panel header."""
        self.query_one("#lib-hdr", Static).update(
            " * Library" if dirty else " Library"
        )

    @on(Button.Pressed, "#btn-lib-del")
    def _btn_del(self):
        t = self.query_one("#lib-table", DataTable)
        if t.row_count == 0:
            return
        row_keys = list(t.rows.keys())
        if not (0 <= t.cursor_row < len(row_keys)):
            return
        rk = row_keys[t.cursor_row]
        entries = [e for e in _load_library() if e.get("id") != rk.value]
        _save_library(entries)
        self._repopulate()


# ── Sequence panel ─────────────────────────────────────────────────────────────

class SequencePanel(Widget):
    """
    Bottom DNA sequence viewer.

    Click  on the sequence → place cursor + select feature at that position.
    Shift+click           → extend selection for editing.
    Shift+E               → open insert/replace dialog at cursor / selection.
    """

    DEFAULT_CSS = """
    SequencePanel {
        height: 14;
        border-top: solid $primary;
    }
    #seq-hdr    { background: $primary; padding: 0 1; height: 1; }
    #seq-scroll { height: 1fr; }
    #seq-trans  {
        height: 3; border-top: solid $primary-darken-2;
        padding: 0 1; display: none;
    }
    SequencePanel.has-trans #seq-trans { display: block; }
    """

    # ── Messages ───────────────────────────────────────────────────────────────

    class SequenceChanged(Message):
        """Emitted when the sequence is modified (commit=True = full rebuild)."""
        def __init__(self, seq: str, commit: bool = False):
            self.seq    = seq
            self.commit = commit
            super().__init__()

    class SequenceClick(Message):
        """User clicked on a base; app should select feature there."""
        def __init__(self, bp: int, double: bool = False):
            self.bp     = bp
            self.double = double
            super().__init__()

    # ── Lifecycle ──────────────────────────────────────────────────────────────

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._seq:          str                     = ""
        self._feats:        list[dict]              = []
        self._sel_range:    "tuple[int,int] | None" = None  # feature highlight
        self._user_sel:     "tuple[int,int] | None" = None  # drag/shift selection
        self._cursor_pos:   int                     = -1    # -1 = no cursor
        self._view_cache_key: "tuple | None"        = None
        self._view_cache_txt: "Text | None"         = None
        # Drag-to-select state
        self._drag_start_bp:    int  = -1
        self._has_dragged:      bool = False
        self._mouse_button_held: bool = False
        self._drag_was_shift:   bool = False
        self._last_was_drag:    bool = False

    def compose(self) -> ComposeResult:
        yield Static(
            " Sequence  (click: select · Shift+click: select region · Shift+E: edit)",
            id="seq-hdr",
        )
        with ScrollableContainer(id="seq-scroll"):
            yield Static("", id="seq-view")
        yield Static("", id="seq-trans")

    # ── Public API ─────────────────────────────────────────────────────────────

    def update_seq(self, seq: str, feats: list[dict]) -> None:
        """Called after loading a record or committing an edit."""
        self._seq        = seq
        self._feats      = feats
        self._sel_range  = None
        self._user_sel   = None
        self._cursor_pos = -1
        self.remove_class("has-trans")
        self._refresh_view()

    def highlight_feature(self, feat: "dict | None") -> None:
        """Highlight a feature's region in the sequence; show CDS translation."""
        if feat is None or not self._seq:
            self._sel_range = None
            self.remove_class("has-trans")
            self._refresh_view()
            return

        start, end = feat["start"], min(feat["end"], len(self._seq))
        self._sel_range = (start, end)
        self._user_sel  = None          # clear shift-selection on programmatic highlight
        self._refresh_view()

        # Scroll to feature start (accounting for annotation bar rows)
        self._scroll_to_row(self._bp_to_content_row(start))

        trans_box = self.query_one("#seq-trans", Static)
        if feat.get("type") == "CDS":
            aa   = _translate_cds(self._seq, start, end, feat.get("strand", 1))
            n_aa = len(aa.rstrip("*"))
            t    = Text()
            t.append(f" {feat['label']}  ", style=f"bold {feat['color']}")
            t.append(f"({n_aa} aa)  ", style="dim")
            t.append(aa[:self.size.width - 4], style=feat["color"])
            trans_box.update(t)
            self.add_class("has-trans")
        else:
            self.remove_class("has-trans")

    def select_feature_range(self, feat: dict) -> None:
        """Double-click highlight: set user_sel to the entire feature span."""
        if not self._seq or feat is None:
            return
        start = feat["start"]
        end   = min(feat["end"], len(self._seq))
        self._user_sel   = (start, end)
        self._sel_range  = None
        self._cursor_pos = start
        self._refresh_view()
        self._scroll_to_row(self._bp_to_content_row(start))

    # ── Mouse / click ──────────────────────────────────────────────────────────

    def on_mouse_down(self, event: MouseDown) -> None:
        if event.button != 1:
            return
        bp = self._click_to_bp(event.screen_x, event.screen_y)
        if bp < 0:
            return
        self._mouse_button_held = True
        self._drag_start_bp     = bp
        self._has_dragged       = False
        self._drag_was_shift    = event.shift
        if event.shift and self._cursor_pos >= 0:
            # Shift+click: extend selection from cursor to here
            s = min(self._cursor_pos, bp)
            e = max(self._cursor_pos, bp) + 1
            self._user_sel   = (s, e)
            self._cursor_pos = bp
            self._sel_range  = None
        else:
            # Plain click: place cursor, clear selection
            self._cursor_pos = bp
            self._user_sel   = None
        self._refresh_view()
        self._scroll_to_row(self._bp_to_content_row(bp))

    def on_mouse_move(self, event: MouseMove) -> None:
        if not self._mouse_button_held or self._drag_start_bp < 0:
            return
        bp = self._click_to_bp(event.screen_x, event.screen_y)
        if bp < 0 or bp == self._drag_start_bp and not self._has_dragged:
            return
        self._has_dragged = True
        s = min(self._drag_start_bp, bp)
        e = max(self._drag_start_bp, bp) + 1
        self._user_sel   = (s, e)
        self._cursor_pos = bp
        self._sel_range  = None
        self._refresh_view()

    def on_mouse_up(self, event: MouseUp) -> None:
        if event.button != 1:
            return
        self._last_was_drag     = self._has_dragged
        self._mouse_button_held = False
        self._drag_start_bp     = -1
        self._has_dragged       = False
        self._drag_was_shift    = False

    def on_click(self, event: Click) -> None:
        """Handle single/double click — fires after mouse_up at same position."""
        if self._last_was_drag:
            self._last_was_drag = False
            return
        self._last_was_drag = False
        bp = self._click_to_bp(event.screen_x, event.screen_y)
        if bp < 0:
            return
        double = event.chain >= 2
        self.post_message(self.SequenceClick(bp, double=double))

    def _click_to_bp(self, screen_x: int, screen_y: int) -> int:
        """Map absolute screen coords to a bp index, or -1 if not on a base."""
        if not self._seq:
            return -1
        try:
            scroll      = self.query_one("#seq-scroll", ScrollableContainer)
            reg         = scroll.region          # absolute screen position of viewport
            vp_x        = screen_x - reg.x       # column within the viewport
            vp_y        = screen_y - reg.y       # visible row within the viewport
            if vp_x < 0 or vp_y < 0 or vp_x >= reg.width or vp_y >= reg.height:
                return -1
            content_row = vp_y + int(scroll.scroll_y)
        except Exception:
            return -1

        line_width  = max(20, self.size.width - 14)
        annot_feats = sorted(
            [f for f in self._feats if f.get("type") != "site"],
            key=lambda f: -(f["end"] - f["start"]),
        )
        n   = len(self._seq)
        row = 0
        for chunk_start in range(0, n, line_width):
            chunk_end = min(chunk_start + line_width, n)
            if row == content_row:
                # Clicked on a DNA row
                seq_col = vp_x - 10   # 10 = 8-char line num + 2 spaces
                if 0 <= seq_col <= (chunk_end - chunk_start):
                    return chunk_start + seq_col
                return -1
            row += 1
            # Annotation bar rows
            chunk_annot = [
                f for f in annot_feats
                if f["start"] < chunk_end and f["end"] > chunk_start
            ]
            for f in chunk_annot[:4]:
                if row == content_row:
                    return (f["start"] + f["end"]) // 2
                row += 1
            if row > content_row:
                break
        return -1

    # ── Helpers ────────────────────────────────────────────────────────────────

    def _annot_feats_sorted(self) -> list:
        return sorted(
            [f for f in self._feats if f.get("type") != "site"],
            key=lambda f: -(f["end"] - f["start"]),
        )

    def _bp_to_content_row(self, bp: int) -> int:
        """Return the content row index (0-based) that contains bp."""
        line_width  = max(20, self.size.width - 14)
        annot_feats = self._annot_feats_sorted()
        n   = len(self._seq)
        row = 0
        for chunk_start in range(0, n, line_width):
            chunk_end = min(chunk_start + line_width, n)
            if bp < chunk_end:
                return row
            row += 1
            n_bars = sum(
                1 for f in annot_feats
                if f["start"] < chunk_end and f["end"] > chunk_start
            )
            row += min(n_bars, 4)
        return row

    def _scroll_to_row(self, row: int) -> None:
        try:
            self.query_one("#seq-scroll", ScrollableContainer).scroll_to(
                0, row, animate=False
            )
        except Exception:
            pass

    def _refresh_view(self) -> None:
        view = self.query_one("#seq-view", Static)
        if not self._seq:
            view.update(Text("  No sequence loaded.", style="dim italic"))
            return
        line_width = max(20, self.size.width - 14)
        key = (id(self._seq), id(self._feats), line_width,
               self._sel_range, self._user_sel, self._cursor_pos)
        if key != self._view_cache_key:
            self._view_cache_txt = _build_seq_text(
                self._seq, self._feats,
                line_width = line_width,
                sel_range  = self._sel_range,
                user_sel   = self._user_sel,
                cursor_pos = self._cursor_pos,
            )
            self._view_cache_key = key
        view.update(self._view_cache_txt)

    def on_resize(self, _) -> None:
        self._refresh_view()


# ── Sequence edit dialog ───────────────────────────────────────────────────────

class EditSeqDialog(ModalScreen):
    """Insert or replace a sequence region via a small dialog.

    mode="insert"  — inserts new_seq at position start
    mode="replace" — replaces seq[start:end] with new_seq
    Dismisses with (new_seq, mode, start, end) on OK, or None on cancel.
    """

    _VALID = frozenset("ATCGNRYSWKMBDHV")   # IUPAC DNA codes

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, mode: str, existing: str = "",
                 start: int = 0, end: int = 0):
        super().__init__()
        self._mode     = mode
        self._existing = existing
        self._start    = start
        self._end      = end

    def compose(self) -> ComposeResult:
        with Vertical(id="edit-dlg"):
            if self._mode == "insert":
                yield Static(
                    f" Insert sequence at position {self._start + 1} ",
                    id="edit-title",
                )
            else:
                excerpt = self._existing[:35] + ("…" if len(self._existing) > 35 else "")
                yield Static(
                    f" Replace {self._end - self._start} bp  "
                    f"({self._start + 1}‥{self._end}) ",
                    id="edit-title",
                )
                yield Static(f"  Current: {excerpt}", id="edit-current")
            yield Label("New sequence  (A T C G N and IUPAC codes only):")
            yield Input(id="edit-input", placeholder="ATCG…")
            yield Static("", id="edit-err")
            with Horizontal(id="edit-btns"):
                yield Button("OK", id="btn-ok", variant="primary")
                yield Button("Cancel", id="btn-cancel")

    def on_mount(self) -> None:
        self.query_one("#edit-input", Input).focus()

    @on(Input.Changed, "#edit-input")
    def _live_validate(self, event: Input.Changed) -> None:
        val = event.value.strip().upper()
        err = self.query_one("#edit-err", Static)
        if not val:
            err.update("")
            return
        bad = sorted(set(c for c in val if c not in self._VALID))
        if bad:
            chars = "  ".join(repr(c) for c in bad)
            err.update(Text(
                f"Invalid: {chars} — only A T C G N (IUPAC) allowed",
                style="bold red",
            ))
        else:
            err.update(Text(f"{len(val)} bp", style="dim green"))

    def _try_submit(self) -> None:
        val = self.query_one("#edit-input", Input).value.strip().upper()
        err = self.query_one("#edit-err", Static)
        if not val:
            err.update(Text("Please enter a sequence.", style="bold red"))
            return
        bad = [c for c in val if c not in self._VALID]
        if bad:
            return   # live validation already shows the error
        self.dismiss((val, self._mode, self._start, self._end))

    @on(Button.Pressed, "#btn-ok")
    def _ok(self, _) -> None:
        self._try_submit()

    @on(Button.Pressed, "#btn-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    @on(Input.Submitted, "#edit-input")
    def _submitted(self, _) -> None:
        self._try_submit()

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Fetch modal ────────────────────────────────────────────────────────────────

class FetchModal(ModalScreen):
    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def compose(self) -> ComposeResult:
        with Vertical(id="fetch-box"):
            yield Static(" Fetch GenBank Record ", id="fetch-title")
            yield Label("NCBI Accession  (e.g.  M77789  for pUC19,  Y14837  for pUC57):")
            yield Input(placeholder="MW463917.1", id="fetch-acc", value="MW463917.1")
            yield Label("Email  (required by NCBI Entrez):")
            yield Input(placeholder="you@example.com", id="fetch-email")
            with Horizontal(id="fetch-btns"):
                yield Button("Fetch", id="btn-fetch", variant="primary")
                yield Button("Cancel", id="btn-cancel-fetch")
            yield Static("", id="fetch-status", markup=True)

    @on(Button.Pressed, "#btn-fetch")
    def _fetch(self):
        acc   = self.query_one("#fetch-acc",   Input).value.strip()
        email = self.query_one("#fetch-email", Input).value.strip() or "splicecraft@local"
        if not acc:
            self.query_one("#fetch-status", Static).update("[red]Enter an accession.[/red]")
            return
        self.query_one("#fetch-status", Static).update(
            f"[dim]Fetching {acc!r} from NCBI…[/dim]"
        )
        self._do_fetch(acc, email)

    @work(thread=True)
    def _do_fetch(self, acc: str, email: str):
        try:
            record = fetch_genbank(acc, email)
            self.app.call_from_thread(self.dismiss, record)
        except Exception as exc:
            def _err():
                self.query_one("#fetch-status", Static).update(
                    f"[red]Error: {exc}[/red]"
                )
            self.app.call_from_thread(_err)

    @on(Input.Submitted)
    def _submitted(self):
        self.query_one("#btn-fetch", Button).press()

    @on(Button.Pressed, "#btn-cancel-fetch")
    def _cancel_btn(self):
        self.dismiss(None)

    def action_cancel(self):
        self.dismiss(None)


# ── Open-file modal ────────────────────────────────────────────────────────────

class OpenFileModal(ModalScreen):
    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def compose(self) -> ComposeResult:
        with Vertical(id="open-box"):
            yield Static(" Open GenBank File ", id="open-title")
            yield Label("File path  (.gb / .gbk):")
            yield Input(placeholder="/path/to/plasmid.gb", id="open-path")
            with Horizontal(id="open-btns"):
                yield Button("Open", id="btn-open", variant="primary")
                yield Button("Cancel", id="btn-cancel-open")
            yield Static("", id="open-status", markup=True)

    @on(Button.Pressed, "#btn-open")
    def _open(self):
        path = self.query_one("#open-path", Input).value.strip()
        if not path:
            self.query_one("#open-status", Static).update("[red]Enter a file path.[/red]")
            return
        try:
            record = load_genbank(path)
            record._tui_source = path   # remember where it came from
            self.dismiss(record)
        except Exception as exc:
            self.query_one("#open-status", Static).update(f"[red]{exc}[/red]")

    @on(Input.Submitted)
    def _submitted(self):
        self.query_one("#btn-open", Button).press()

    @on(Button.Pressed, "#btn-cancel-open")
    def _cancel_btn(self):
        self.dismiss(None)

    def action_cancel(self):
        self.dismiss(None)


# ── Unsaved-changes quit dialog ────────────────────────────────────────────────

class UnsavedQuitModal(ModalScreen):
    """Shown when the user tries to quit with unsaved edits."""

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def compose(self) -> ComposeResult:
        with Vertical(id="quit-dlg"):
            yield Static(" Unsaved Changes ", id="quit-title")
            yield Static(
                "  You have unsaved edits. What would you like to do?",
                id="quit-msg",
            )
            with Horizontal(id="quit-btns"):
                yield Button("Save & Quit",      id="btn-save-quit", variant="primary")
                yield Button("Abandon Changes",  id="btn-abandon",   variant="error")
                yield Button("Cancel",           id="btn-cancel-quit")

    @on(Button.Pressed, "#btn-save-quit")
    def _save_quit(self, _): self.dismiss("save")

    @on(Button.Pressed, "#btn-abandon")
    def _abandon(self, _):   self.dismiss("abandon")

    @on(Button.Pressed, "#btn-cancel-quit")
    def _cancel_btn(self, _): self.dismiss(None)

    def action_cancel(self): self.dismiss(None)


# ── Main app ───────────────────────────────────────────────────────────────────

class PlasmidApp(App):
    TITLE       = "SpliceCraft"
    TRANSITIONS = {}          # instant screen open/close — no slide animations
    _preload_record = None
    _current_record = None   # last-loaded SeqRecord
    _source_path:   "str | None" = None   # file the current record was loaded from
    _unsaved:        bool         = False  # True when there are unsaved edits
    _MAX_UNDO = 50

    CSS = """
Screen { background: $background; }

/* ── Layout ─────────────────────────────────────────────── */
#main-row   { height: 1fr; }
#center-col { width: 1fr; height: 1fr; }
#map-row    { height: 1fr; }

#status-bar {
    height: 1;
    background: $primary-darken-2;
    color: $text-muted;
    padding: 0 1;
}

/* ── Fetch modal ─────────────────────────────────────────── */
FetchModal { align: center middle; }
#fetch-box {
    width: 70; height: auto; max-height: 90%;
    background: $surface; border: solid $accent; padding: 1 2;
}
#fetch-title { background: $accent; padding: 0 1; margin-bottom: 1; }
#fetch-box Label { color: $text-muted; margin-top: 1; }
#fetch-btns { height: 3; margin-top: 1; }
#fetch-btns Button { margin-right: 1; }
#fetch-status { height: 1; margin-top: 1; }

/* ── Open-file modal ─────────────────────────────────────── */
OpenFileModal { align: center middle; }
#open-box {
    width: 70; height: auto;
    background: $surface; border: solid $primary; padding: 1 2;
}
#open-title { background: $primary; padding: 0 1; margin-bottom: 1; }
#open-box Label { color: $text-muted; margin-top: 1; }
#open-btns { height: 3; margin-top: 1; }
#open-btns Button { margin-right: 1; }
#open-status { height: 1; margin-top: 1; }

/* ── Edit-sequence dialog ────────────────────────────────── */
EditSeqDialog { align: center middle; }
#edit-dlg {
    width: 72; height: auto;
    background: $surface; border: solid $warning; padding: 1 2;
}
#edit-title   { background: $warning-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#edit-current { color: $text-muted; margin-bottom: 1; }
#edit-dlg Label { color: $text-muted; margin-top: 1; }
#edit-err  { height: 2; margin-top: 1; }
#edit-btns { height: 3; margin-top: 1; }
#edit-btns Button { margin-right: 1; }

/* ── Unsaved-quit dialog ─────────────────────────────────── */
UnsavedQuitModal { align: center middle; }
#quit-dlg {
    width: 60; height: auto;
    background: $surface; border: solid $error; padding: 1 2;
}
#quit-title { background: $error-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#quit-msg   { color: $text-muted; margin-bottom: 1; }
#quit-btns  { height: 3; margin-top: 1; }
#quit-btns Button { margin-right: 1; }
"""

    BINDINGS = [
        Binding("f",           "fetch",            "Fetch GenBank", show=True),
        Binding("o",           "open_file",        "Open .gb file", show=True),
        Binding("a",           "add_to_library",   "Add to lib",    show=True),
        Binding("E",           "edit_seq",         "Edit seq",      show=True,  key_display="E"),
        Binding("S",           "save",             "Save",          show=True,  key_display="S"),
        Binding("[",           "rotate_cw",        "Rotate ←",      show=True,  priority=True),
        Binding("]",           "rotate_ccw",       "Rotate →",      show=True,  priority=True),
        Binding("shift+[",     "rotate_cw_lg",     "Rotate ←←",     show=False, priority=True),
        Binding("shift+]",     "rotate_ccw_lg",    "Rotate →→",     show=False, priority=True),
        Binding("home",        "reset_origin",     "Reset origin",  show=True,  priority=True),
        Binding("q",           "quit",             "Quit",          show=True),
    ]

    def compose(self) -> ComposeResult:
        yield Header()
        with Horizontal(id="main-row"):
            yield LibraryPanel(id="library")
            with Vertical(id="center-col"):
                with Horizontal(id="map-row"):
                    yield PlasmidMap(id="plasmid-map")
                    yield FeatureSidebar(id="sidebar")
                yield SequencePanel(id="seq-panel")
        yield Static(
            Text(
                "  [ ] rotate   ← → cursor/map   Shift coarse   Home reset"
                "   f fetch   o open   a add-to-lib   E edit seq   S save   , / . circle",
                style="color(245)",
                no_wrap=True,
            ),
            id="status-bar",
        )
        yield Footer()

    # ── Delegate rotation keys to PlasmidMap ───────────────────────────────────

    def action_rotate_cw(self):
        self.query_one("#plasmid-map", PlasmidMap).action_rotate_cw()

    def action_rotate_ccw(self):
        self.query_one("#plasmid-map", PlasmidMap).action_rotate_ccw()

    def action_rotate_cw_lg(self):
        self.query_one("#plasmid-map", PlasmidMap).action_rotate_cw_lg()

    def action_rotate_ccw_lg(self):
        self.query_one("#plasmid-map", PlasmidMap).action_rotate_ccw_lg()

    def action_reset_origin(self):
        self.query_one("#plasmid-map", PlasmidMap).action_reset_origin()

    def action_edit_seq(self) -> None:
        sp = self.query_one("#seq-panel", SequencePanel)
        if not sp._seq:
            self.notify("No sequence loaded.", severity="warning")
            return
        if sp._user_sel is not None:
            # Replace the shift-selected region
            s, e     = sp._user_sel
            existing = sp._seq[s:e]
            self.push_screen(
                EditSeqDialog("replace", existing, s, e),
                callback=self._edit_dialog_result,
            )
        elif sp._cursor_pos >= 0:
            # Insert at cursor position
            pos = sp._cursor_pos
            self.push_screen(
                EditSeqDialog("insert", start=pos, end=pos),
                callback=self._edit_dialog_result,
            )
        else:
            self.notify(
                "Click on the sequence to place a cursor, "
                "or Shift+click to select a region.",
                severity="information",
            )

    def _edit_dialog_result(self, result) -> None:
        if result is None:
            return
        self._push_undo()
        new_bases, mode, s, e = result
        sp      = self.query_one("#seq-panel", SequencePanel)
        pm      = self.query_one("#plasmid-map", PlasmidMap)
        old_seq = sp._seq

        if mode == "insert":
            new_seq = old_seq[:s] + new_bases + old_seq[s:]
        else:
            new_seq = old_seq[:s] + new_bases + old_seq[e:]

        new_cursor = s + len(new_bases)

        if self._current_record is not None:
            new_record = self._rebuild_record_with_edit(new_seq, mode, s, e, new_bases)
            self._current_record = new_record
            pm.load_record(new_record)
            self.query_one("#sidebar", FeatureSidebar).populate(pm._feats)
            restr = _scan_restriction_sites(new_seq)
            pm._restr_feats = restr
            pm.refresh()
            sp.update_seq(new_seq, pm._feats + restr)
            self.notify(f"Sequence updated  ({len(new_seq):,} bp)")
        else:
            restr = _scan_restriction_sites(new_seq)
            pm._restr_feats = restr
            pm.refresh()
            sp.update_seq(new_seq, pm._feats + restr)

        self._mark_dirty()

        # Restore cursor after update_seq resets it
        sp._cursor_pos = new_cursor
        sp._user_sel   = None
        sp._refresh_view()

    def _rebuild_record_with_edit(self, new_seq: str, mode: str,
                                   s: int, e: int, new_bases: str):
        """Rebuild SeqRecord after an insert/replace, shifting feature coords precisely."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        ins_len  = len(new_bases)
        del_len  = 0 if mode == "insert" else (e - s)
        delta    = ins_len - del_len
        new_len  = len(new_seq)
        # For insert: the edit region is [s, s); e == s
        # For replace: the edit region is [s, e)
        edit_end = s if mode == "insert" else e

        new_record = SeqRecord(
            Seq(new_seq),
            id=self._current_record.id,
            name=self._current_record.name,
            description=self._current_record.description,
            annotations=dict(self._current_record.annotations),
        )

        for feat in self._current_record.features:
            fs = int(feat.location.start)
            fe = int(feat.location.end)

            if mode == "insert":
                if fe <= s:
                    new_fs, new_fe = fs, fe           # entirely before insert
                elif fs >= s:
                    new_fs, new_fe = fs + ins_len, fe + ins_len  # entirely after
                else:
                    new_fs, new_fe = fs, fe + ins_len  # spans insert → expand
            else:  # replace [s, e)
                if fe <= s:
                    new_fs, new_fe = fs, fe             # entirely before
                elif fs >= e:
                    new_fs, new_fe = fs + delta, fe + delta  # entirely after
                elif fs <= s and fe >= e:
                    new_fs, new_fe = fs, fe + delta     # spans replaced region
                elif fs < s:
                    new_fs, new_fe = fs, s + ins_len    # overlaps start of region
                else:
                    new_fs, new_fe = s, fe + delta      # overlaps end of region

            new_fe = max(new_fs + 1, min(new_fe, new_len))
            new_record.features.append(SeqFeature(
                FeatureLocation(new_fs, new_fe,
                                strand=getattr(feat.location, "strand", 1)),
                type=feat.type,
                qualifiers=dict(feat.qualifiers),
            ))

        return new_record

    def action_add_to_library(self):
        if self._current_record is None:
            self.notify("No record loaded to add.", severity="warning")
            return
        lib = self.query_one("#library", LibraryPanel)
        lib.add_entry(self._current_record)
        self.notify(f"Added {self._current_record.name} to library.")

    # ── Mount: auto-load preloaded record ──────────────────────────────────────

    def on_mount(self) -> None:
        self._undo_stack: list = []
        self._redo_stack: list = []
        if self._preload_record is not None:
            def _load_preload():
                record = self._preload_record
                path   = getattr(record, "_tui_source", None)
                self._apply_record(record)
                self._source_path = path
            self.call_after_refresh(_load_preload)
        elif not _load_library():
            self._seed_default_library()

    @work(thread=True)
    def _seed_default_library(self) -> None:
        """Fetch MW463917.1 and pre-populate the library on first run."""
        try:
            record = fetch_genbank("MW463917.1")
            def _add():
                lib = self.query_one("#library", LibraryPanel)
                lib.add_entry(record)
                self._apply_record(record)
            self.call_from_thread(_add)
        except Exception:
            pass

    # ── Keyboard: cursor movement, copy, undo/redo ─────────────────────────────

    def on_key(self, event) -> None:
        sp = self.query_one("#seq-panel", SequencePanel)

        # ── Ctrl+C: copy selection ────────────────────────────────────────────
        if event.key == "ctrl+c":
            seq = sp._seq
            if seq:
                sel = sp._user_sel or sp._sel_range
                if sel:
                    text = seq[sel[0]:sel[1]].upper()
                    if _copy_to_clipboard(text):
                        self.notify(f"Copied {len(text)} bp to clipboard")
                    else:
                        self.notify("Clipboard unavailable (install xclip / xsel / wl-copy)",
                                    severity="warning")
                else:
                    self.notify("No selection — Shift+click or double-click a feature first",
                                severity="information")
            event.stop()
            return

        # ── Ctrl+Z: undo ──────────────────────────────────────────────────────
        if event.key == "ctrl+z":
            self._action_undo()
            event.stop()
            return

        # ── Ctrl+Shift+Z / Ctrl+Y: redo ───────────────────────────────────────
        if event.key in ("ctrl+shift+z", "ctrl+Z", "ctrl+y"):
            self._action_redo()
            event.stop()
            return

        # ── Arrow keys: move sequence cursor ──────────────────────────────────
        if sp._cursor_pos < 0 or not sp._seq:
            return
        if event.key == "left":
            new_pos = max(0, sp._cursor_pos - 1)
        elif event.key == "right":
            new_pos = min(len(sp._seq), sp._cursor_pos + 1)
        else:
            return
        event.stop()
        sp._cursor_pos = new_pos
        sp._user_sel   = None
        sp._refresh_view()
        sp._scroll_to_row(sp._bp_to_content_row(new_pos))

    # ── Undo / redo ────────────────────────────────────────────────────────────

    def _push_undo(self) -> None:
        sp = self.query_one("#seq-panel", SequencePanel)
        if not sp._seq:
            return
        snapshot = (sp._seq, sp._cursor_pos, self._current_record)
        self._undo_stack.append(snapshot)
        if len(self._undo_stack) > self._MAX_UNDO:
            self._undo_stack.pop(0)
        self._redo_stack.clear()

    def _apply_snapshot(self, seq: str, cursor_pos: int, record) -> None:
        pm      = self.query_one("#plasmid-map",  PlasmidMap)
        sidebar = self.query_one("#sidebar",      FeatureSidebar)
        sp      = self.query_one("#seq-panel",    SequencePanel)
        self._current_record = record
        if record is not None:
            pm.load_record(record)
            sidebar.populate(pm._feats)
            restr = _scan_restriction_sites(seq)
            pm._restr_feats = restr
            pm.refresh()
            sp.update_seq(seq, pm._feats + restr)
        else:
            sp.update_seq(seq, [])
        sp._cursor_pos = cursor_pos
        sp._refresh_view()
        self._mark_dirty()

    def _action_undo(self) -> None:
        if not self._undo_stack:
            self.notify("Nothing to undo", severity="information")
            return
        sp = self.query_one("#seq-panel", SequencePanel)
        # Save current state to redo stack
        self._redo_stack.append((sp._seq, sp._cursor_pos, self._current_record))
        seq, cursor_pos, record = self._undo_stack.pop()
        self._apply_snapshot(seq, cursor_pos, record)
        remaining = len(self._undo_stack)
        self.notify(f"Undo  ({remaining} left)")

    def _action_redo(self) -> None:
        if not self._redo_stack:
            self.notify("Nothing to redo", severity="information")
            return
        sp = self.query_one("#seq-panel", SequencePanel)
        self._undo_stack.append((sp._seq, sp._cursor_pos, self._current_record))
        seq, cursor_pos, record = self._redo_stack.pop()
        self._apply_snapshot(seq, cursor_pos, record)
        remaining = len(self._redo_stack)
        self.notify(f"Redo  ({remaining} left)")

    # ── Dirty-state helpers ────────────────────────────────────────────────────

    def _mark_dirty(self) -> None:
        self._unsaved = True
        if self._current_record:
            n = len(self._current_record.seq)
            self.title = f"SpliceCraft — *{self._current_record.name}  ({n:,} bp)"
        try:
            self.query_one("#library", LibraryPanel).set_dirty(True)
        except Exception:
            pass

    def _mark_clean(self) -> None:
        self._unsaved = False
        if self._current_record:
            n = len(self._current_record.seq)
            self.title = f"SpliceCraft — {self._current_record.name}  ({n:,} bp)"
        try:
            self.query_one("#library", LibraryPanel).set_dirty(False)
        except Exception:
            pass

    def _do_save(self) -> bool:
        """Write current record to its source file. Returns True on success."""
        if self._current_record is None:
            self.notify("Nothing to save.", severity="warning")
            return False
        if not self._source_path:
            self.notify(
                "No source file — open a .gb file first, or use 'a' to add to library.",
                severity="warning",
            )
            return False
        try:
            Path(self._source_path).write_text(_record_to_gb_text(self._current_record))
            self._mark_clean()
            self.notify(f"Saved → {self._source_path}")
            return True
        except Exception as exc:
            self.notify(f"Save failed: {exc}", severity="error")
            return False

    def action_save(self) -> None:
        self._do_save()

    def action_quit(self) -> None:
        if self._unsaved:
            self.push_screen(UnsavedQuitModal(), callback=self._on_quit_response)
        else:
            self.exit()

    def _on_quit_response(self, result) -> None:
        if result == "save":
            if self._do_save():
                self.exit()
        elif result == "abandon":
            self.exit()
        # None = cancel → stay

    # ── Fetch / open ───────────────────────────────────────────────────────────

    def action_fetch(self):
        self.push_screen(FetchModal(), callback=self._apply_record)

    def action_open_file(self):
        def _cb(record):
            if record is None:
                return
            path = getattr(record, "_tui_source", None)
            self._apply_record(record)
            self._source_path = path   # restore after _apply_record clears it
        self.push_screen(OpenFileModal(), callback=_cb)

    # ── Central record loader ──────────────────────────────────────────────────

    def _apply_record(self, record) -> None:
        """Load a SeqRecord into all panels."""
        if record is None:
            return
        self._current_record = record
        self._source_path    = None   # caller sets this if it came from a file

        pm      = self.query_one("#plasmid-map", PlasmidMap)
        sidebar = self.query_one("#sidebar",     FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel",   SequencePanel)

        pm.load_record(record)
        sidebar.populate(pm._feats)

        seq_str = str(record.seq)
        restr   = _scan_restriction_sites(seq_str)

        # Store restriction sites on the map for visual overlay
        pm._restr_feats = restr
        pm.refresh()

        # Sequence panel: feature coloring = record feats + restriction sites
        seq_pnl.update_seq(seq_str, pm._feats + restr)

        self._mark_clean()
        self.notify(
            f"Loaded {record.name}  ({len(record.seq):,} bp, "
            f"{len(pm._feats)} features, {len(restr)} restriction sites)"
        )

    # ── Feature selection: map ↔ sidebar ↔ sequence panel ─────────────────────

    @on(SequencePanel.SequenceClick)
    def _seq_click(self, event: SequencePanel.SequenceClick) -> None:
        """Single or double click on the sequence — select the smallest feature at bp."""
        pm      = self.query_one("#plasmid-map", PlasmidMap)
        sidebar = self.query_one("#sidebar",     FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel",   SequencePanel)
        bp      = event.bp
        best_idx  = -1
        best_span = float("inf")
        for i, f in enumerate(pm._feats):
            s, e = f["start"], f["end"]
            if s <= bp < e and (e - s) < best_span:
                best_span = e - s
                best_idx  = i
        if best_idx >= 0:
            f = pm._feats[best_idx]
            pm.select_feature(best_idx)
            sidebar.show_detail(f)
            sidebar.highlight_row(best_idx)
            if event.double:
                # Double-click: highlight entire feature span as editable selection
                seq_pnl.select_feature_range(f)
            else:
                seq_pnl.highlight_feature(f)

    @on(PlasmidMap.FeatureSelected)
    def _map_feat_selected(self, event: PlasmidMap.FeatureSelected):
        sidebar = self.query_one("#sidebar",   FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel", SequencePanel)
        sidebar.show_detail(event.feat_dict)
        if event.idx >= 0:
            sidebar.highlight_row(event.idx)
        seq_pnl.highlight_feature(event.feat_dict)

    @on(FeatureSidebar.RowActivated)
    def _sidebar_row_activated(self, event: FeatureSidebar.RowActivated):
        pm      = self.query_one("#plasmid-map", PlasmidMap)
        sidebar = self.query_one("#sidebar",     FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel",   SequencePanel)
        pm.select_feature(event.idx)
        f = pm._feats[event.idx] if 0 <= event.idx < len(pm._feats) else None
        sidebar.show_detail(f)
        seq_pnl.highlight_feature(f)

    # ── Library events ─────────────────────────────────────────────────────────

    @on(LibraryPanel.PlasmidLoad)
    def _library_load(self, event: LibraryPanel.PlasmidLoad):
        gb_text = event.entry.get("gb_text", "")
        if not gb_text:
            self.notify(f"Library entry has no stored sequence.", severity="warning")
            return
        try:
            record = _gb_text_to_record(gb_text)
            self._apply_record(record)
        except Exception as exc:
            self.notify(f"Failed to load from library: {exc}", severity="error")

    @on(LibraryPanel.AddCurrentRequested)
    def _library_add_current(self, _):
        self.action_add_to_library()

    # ── Sequence edits ─────────────────────────────────────────────────────────

    @on(SequencePanel.SequenceChanged)
    def _seq_changed(self, event: SequencePanel.SequenceChanged):
        # Update restriction site overlay whenever sequence changes
        pm    = self.query_one("#plasmid-map", PlasmidMap)
        restr = _scan_restriction_sites(event.seq)
        pm._restr_feats = restr
        pm.refresh()


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    app = PlasmidApp()

    if arg:
        if Path(arg).exists():
            try:
                record = load_genbank(arg)
                record._tui_source = str(Path(arg).resolve())
            except Exception as exc:
                print(f"Could not load {arg!r}: {exc}", file=sys.stderr)
                sys.exit(1)
        else:
            print(f"Fetching {arg!r} from NCBI…", flush=True)
            try:
                record = fetch_genbank(arg)
                print(f"  Got: {record.name}  ({len(record.seq)} bp)")
            except Exception as exc:
                print(f"Fetch failed: {exc}", file=sys.stderr)
                sys.exit(1)
        app._preload_record = record

    app.run()


if __name__ == "__main__":
    main()
