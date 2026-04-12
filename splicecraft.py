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
import logging
import math
import os
import platform
import re
import sys
import uuid as _uuid
from io import StringIO
from logging.handlers import RotatingFileHandler
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

# ── Logging ────────────────────────────────────────────────────────────────────
# Borrowed from ScriptoScope: rotating file log with an 8-char session ID prefix
# on every line so multi-run logs are greppable. Default path /tmp/splicecraft.log,
# overridable via $SPLICECRAFT_LOG. UI never sees raw tracebacks — they go here.

_LOG_PATH   = os.environ.get("SPLICECRAFT_LOG") or "/tmp/splicecraft.log"
_SESSION_ID = _uuid.uuid4().hex[:8]

class _SessionFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        record.session = _SESSION_ID
        return True

_log = logging.getLogger("splicecraft")
_log.setLevel(logging.INFO)
_log.propagate = False
if not _log.handlers:
    try:
        _handler = RotatingFileHandler(
            _LOG_PATH, maxBytes=2 * 1024 * 1024, backupCount=2, encoding="utf-8"
        )
        _handler.setFormatter(logging.Formatter(
            fmt="%(asctime)s [%(session)s] %(levelname)-5s "
                "%(name)s.%(funcName)s:%(lineno)d %(message)s"
        ))
        _handler.addFilter(_SessionFilter())
        _log.addHandler(_handler)
    except OSError:
        # Fall back to a no-op handler if /tmp is read-only (rare; shouldn't crash UI)
        _log.addHandler(logging.NullHandler())

def _log_startup_banner() -> None:
    def _ver(import_name: str) -> str:
        try:
            mod = __import__(import_name)
            return getattr(mod, "__version__", "unknown")
        except ImportError:
            return "NOT INSTALLED"
    _log.info("=" * 60)
    _log.info("SpliceCraft session %s starting", _SESSION_ID)
    _log.info("python    : %s", sys.version.split()[0])
    _log.info("platform  : %s", platform.platform())
    _log.info("textual   : %s", _ver("textual"))
    _log.info("biopython : %s", _ver("Bio"))
    _log.info("log path  : %s", _LOG_PATH)
    _log.info("=" * 60)

from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, ScrollableContainer
from textual.events import Click, MouseDown, MouseMove, MouseUp, MouseScrollDown, MouseScrollUp
from textual.message import Message
from textual.reactive import reactive
from textual.screen import ModalScreen, Screen
from textual.widget import Widget
from textual.widgets import (
    Button, DataTable, Footer, Header, Input, Label, Select, Static,
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

# ── Safe JSON persistence ──────────────────────────────────────────────────────
#
# All user data (plasmid library, parts bin, primer library) goes through
# _safe_save_json which:
#   1. Backs up the existing file to *.bak BEFORE overwriting
#   2. Writes via tempfile + os.replace (atomic on POSIX — the file is either
#      fully written or not at all; no partial-write corruption)
#   3. Logs every write with entry count for post-mortem debugging
#
# _safe_load_json handles the read side:
#   - Missing file → [] (first run, not an error)
#   - Corrupt file → attempt restore from .bak; if .bak also corrupt → []
#   - Returns (entries, warning_message_or_None)

def _safe_save_json(path: Path, entries: list, label: str) -> None:
    """Atomically write `entries` as JSON to `path`, backing up first.

    The .bak file is the user's safety net — if a write goes wrong or the
    app crashes mid-save, the previous version survives as path.bak.
    """
    import os
    import tempfile

    # 1. Back up the existing file (if it has content)
    if path.exists():
        try:
            existing = path.read_bytes()
            if existing.strip():
                bak = path.with_suffix(path.suffix + ".bak")
                bak.write_bytes(existing)
        except OSError:
            _log.warning("Could not create backup for %s", path)

    # 2. Atomic write: tempfile in same dir → os.replace
    try:
        fd, tmp_name = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                json.dump(entries, fh, indent=2)
                fh.flush()
                try:
                    os.fsync(fh.fileno())
                except OSError:
                    pass
            os.replace(tmp_name, str(path))
            _log.info("Saved %s: %d entries to %s", label, len(entries), path)
        except Exception:
            # Clean up the temp file on failure
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
            raise
    except Exception:
        _log.exception("Failed to save %s to %s", label, path)


def _safe_load_json(path: Path, label: str) -> "tuple[list, str | None]":
    """Load a JSON array from `path`. Returns (entries, warning_or_None).

    - Missing file → ([], None) — normal first run, no warning.
    - Valid file   → (entries, None).
    - Corrupt file → attempt .bak restore; if .bak is valid →
      (bak_entries, warning). If .bak also corrupt → ([], warning).
    """
    if not path.exists():
        return [], None

    # Try the main file
    try:
        entries = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(entries, list):
            return entries, None
        _log.warning("%s: expected list, got %s", path, type(entries).__name__)
    except Exception:
        _log.exception("Corrupt %s file: %s", label, path)

    # Main file is corrupt — try the .bak
    bak = path.with_suffix(path.suffix + ".bak")
    if bak.exists():
        try:
            entries = json.loads(bak.read_text(encoding="utf-8"))
            if isinstance(entries, list):
                _log.info("Restored %s from backup %s (%d entries)",
                          label, bak, len(entries))
                # Overwrite the corrupt main file with the good backup
                try:
                    import shutil
                    shutil.copy2(str(bak), str(path))
                except OSError:
                    pass
                return entries, (
                    f"{label} was corrupt — restored {len(entries)} entries "
                    f"from backup."
                )
        except Exception:
            _log.exception("Backup %s also corrupt: %s", label, bak)

    return [], f"{label} is corrupt and no valid backup was found. Starting empty."


# ── Library persistence ────────────────────────────────────────────────────────

_LIBRARY_FILE = Path(__file__).parent / "plasmid_library.json"
_library_cache: "list | None" = None

def _load_library() -> list[dict]:
    global _library_cache
    if _library_cache is not None:
        return list(_library_cache)
    entries, warning = _safe_load_json(_LIBRARY_FILE, "Plasmid library")
    if warning:
        _log.warning(warning)
    _library_cache = entries
    return list(_library_cache)

def _save_library(entries: list[dict]) -> None:
    global _library_cache
    _safe_save_json(_LIBRARY_FILE, entries, "Plasmid library")
    _library_cache = list(entries)

# ── Restriction sites ──────────────────────────────────────────────────────────

# Comprehensive NEB enzyme catalog.
# Each entry maps enzyme name → (recognition_sequence, fwd_cut, rev_cut).
#
# Cut position convention (0-based, relative to start of recognition sequence):
#   0          = before first base of recognition seq
#   len(seq)   = after last base of recognition seq
#   >len(seq)  = downstream of recognition seq (Type IIS)
#   negative   = upstream of recognition seq
#
# IUPAC ambiguity codes used in recognition sequences:
#   R=A/G  Y=C/T  W=A/T  S=G/C  M=A/C  K=G/T
#   B=not-A  D=not-C  H=not-G  V=not-T  N=any
#
# Examples:
#   EcoRI  G^AATTC / CTTAA^G  → fwd=1, rev=5
#   BamHI  G^GATCC / CCTAG^G  → fwd=1, rev=5
#   SmaI   CCC^GGG / GGG^CCC  → fwd=3, rev=3  (blunt)
#   BsaI   GGTCTC(1/5)        → recognition=6bp, fwd=6+1=7, rev=6+5=11
#
_NEB_ENZYMES: dict[str, tuple[str, int, int]] = {

    # ── Common Type IIP — 6-bp palindromic cutters ─────────────────────────────
    "EcoRI":     ("GAATTC",       1,  5),  # G^AATTC   / CTTAA^G     5' overhang
    "EcoRV":     ("GATATC",       3,  3),  # GAT^ATC   / GAT^ATC     blunt
    "BamHI":     ("GGATCC",       1,  5),  # G^GATCC   / CCTAG^G     5' overhang
    "HindIII":   ("AAGCTT",       1,  5),  # A^AGCTT   / TTCGA^A     5' overhang
    "NcoI":      ("CCATGG",       1,  5),  # C^CATGG   / GGTAC^C     5' overhang
    "NdeI":      ("CATATG",       2,  4),  # CA^TATG   / GTAT^AC     5' overhang
    "XhoI":      ("CTCGAG",       1,  5),  # C^TCGAG   / GAGCT^C     5' overhang
    "SalI":      ("GTCGAC",       1,  5),  # G^TCGAC   / CAGCT^G     5' overhang
    "KpnI":      ("GGTACC",       5,  1),  # GGTAC^C   / G^GTACC     3' overhang
    "SacI":      ("GAGCTC",       5,  1),  # GAGCT^C   / G^AGCTC     3' overhang
    "SacII":     ("CCGCGG",       4,  2),  # CCGC^GG   / CC^GCGG     3' overhang
    "SpeI":      ("ACTAGT",       1,  5),  # A^CTAGT   / TGATC^A     5' overhang
    "XbaI":      ("TCTAGA",       1,  5),  # T^CTAGA   / AGATC^T     5' overhang
    "NotI":      ("GCGGCCGC",     2,  6),  # GC^GGCCGC / CGCCGG^CG   5' overhang (8-cutter)
    "PstI":      ("CTGCAG",       5,  1),  # CTGCA^G   / G^CTGCA     3' overhang
    "SphI":      ("GCATGC",       5,  1),  # GCATG^C   / C^GCATG     3' overhang
    "ClaI":      ("ATCGAT",       2,  4),  # AT^CGAT   / CGAT^AT     5' overhang
    "NheI":      ("GCTAGC",       1,  5),  # G^CTAGC   / CGATC^G     5' overhang
    "AvaI":      ("CYCGRG",       1,  5),  # C^YCGRG                 5' overhang (degenerate)
    "AvaII":     ("GGWCC",        1,  4),  # G^GWCC    / CCWG^G      5' overhang
    "AvrII":     ("CCTAGG",       1,  5),  # C^CTAGG   / GGATC^C     5' overhang
    "BclI":      ("TGATCA",       1,  5),  # T^GATCA   / AGTCA^T     5' overhang (dam-sensitive)
    "BglII":     ("AGATCT",       1,  5),  # A^GATCT   / TCTAG^A     5' overhang
    "BsiWI":     ("CGTACG",       1,  5),  # C^GTACG   / GCATG^C     5' overhang
    "BspEI":     ("TCCGGA",       1,  5),  # T^CCGGA   / AGGCC^T     5' overhang
    "BsrGI":     ("TGTACA",       1,  5),  # T^GTACA   / ACATG^T     5' overhang
    "BssHII":    ("GCGCGC",       1,  5),  # G^CGCGC   / CGCGC^G     5' overhang
    "BstBI":     ("TTCGAA",       2,  4),  # TT^CGAA   / AAGC^TT     5' overhang
    "BstEII":    ("GGTNACC",      1,  5),  # G^GTNACC  / CCANT^GG    5' overhang
    "BstXI":     ("CCANNNNNTGG",  8,  4),  # CCANN4^TGG/ CCA^N5TGG   3' overhang
    "BstYI":     ("RGATCY",       1,  5),  # R^GATCY   / YCTAG^R     5' overhang
    "CpoI":      ("CGGWCCG",      2,  5),  # CG^GWCCG  / CGCC^WGG    5' overhang
    "DraI":      ("TTTAAA",       3,  3),  # TTT^AAA   / TTT^AAA     blunt
    "DraIII":    ("CACNNNGTG",    6,  3),  # CACNNN^GTG/ GTG^NNNGTG  3' overhang
    "EagI":      ("CGGCCG",       1,  5),  # C^GGCCG   / GCCGG^C     5' overhang (NotI subset)
    "Eco47III":  ("AGCGCT",       3,  3),  # AGC^GCT                 blunt
    "Eco53kI":   ("GAGCTC",       5,  1),  # GAGCT^C                 3' overhang (SacI isoschizomer)
    "EcoNI":     ("CCTNNNNNAGG",  5,  6),  # CCTNN^NNNAGG            5' overhang
    "FseI":      ("GGCCGGCC",     6,  2),  # GGCCGG^CC / CC^GGCCGG   3' overhang (8-cutter)
    "FspI":      ("TGCGCA",       3,  3),  # TGC^GCA                 blunt
    "HaeII":     ("RGCGCY",       5,  1),  # RGCGC^Y   / R^GCGCY     3' overhang
    "HaeIII":    ("GGCC",         2,  2),  # GG^CC                   blunt (4-cutter)
    "HincII":    ("GTYRAC",       3,  3),  # GTY^RAC                 blunt
    "HindII":    ("GTYRAC",       3,  3),  # GTY^RAC                 blunt (HincII isoschizomer)
    "HpaI":      ("GTTAAC",       3,  3),  # GTT^AAC                 blunt
    "HpaII":     ("CCGG",         1,  3),  # C^CGG     / CGG^C        5' overhang (4-cutter)
    "MfeI":      ("CAATTG",       1,  5),  # C^AATTG   / GTTAA^C     EcoRI-compatible ends
    "MluI":      ("ACGCGT",       1,  5),  # A^CGCGT   / TGCGC^A     5' overhang
    "MscI":      ("TGGCCA",       3,  3),  # TGG^CCA                 blunt
    "MspI":      ("CCGG",         1,  3),  # C^CGG     / CGG^C        5' overhang (HpaII isoschizomer)
    "MunI":      ("CAATTG",       1,  5),  # C^AATTG                 MfeI isoschizomer
    "NarI":      ("GGCGCC",       2,  4),  # GG^CGCC   / CGCG^CC     3' overhang
    "NruI":      ("TCGCGA",       3,  3),  # TCG^CGA                 blunt
    "NsiI":      ("ATGCAT",       5,  1),  # ATGCA^T   / T^ATGCA     PstI-compatible ends
    "NspI":      ("RCATGY",       5,  1),  # RCATG^Y   / R^CATGY     3' overhang
    "PacI":      ("TTAATTAA",     5,  3),  # TTAAT^TAA / TTA^ATTAA   3' overhang (8-cutter)
    "PaeR7I":    ("CTCGAG",       1,  5),  # C^TCGAG                 XhoI isoschizomer
    "PciI":      ("ACATGT",       1,  5),  # A^CATGT   / TGTAC^A     5' overhang
    "PmeI":      ("GTTTAAAC",     4,  4),  # GTTT^AAAC               blunt (8-cutter)
    "PmlI":      ("CACGTG",       3,  3),  # CAC^GTG                 blunt
    "PscI":      ("ACATGT",       1,  5),  # A^CATGT                 PciI isoschizomer
    "PvuI":      ("CGATCG",       4,  2),  # CGATC^G   / G^CGATC     3' overhang
    "PvuII":     ("CAGCTG",       3,  3),  # CAG^CTG                 blunt
    "RsrII":     ("CGGWCCG",      2,  5),  # CG^GWCCG                CpoI isoschizomer
    "SbfI":      ("CCTGCAGG",     6,  2),  # CCTGCA^GG / CC^TGCAGG   PstI-compatible (8-cutter)
    "ScaI":      ("AGTACT",       3,  3),  # AGT^ACT                 blunt
    "SfiI":      ("GGCCNNNNNGGCC",8,  4),  # GGCCN4^NGGCC            3' overhang (13-bp)
    "SgrAI":     ("CRCCGGYG",     2,  6),  # CR^CCGGYG / GCCGGR^C    5' overhang
    "SmaI":      ("CCCGGG",       3,  3),  # CCC^GGG                 blunt
    "SnaBI":     ("TACGTA",       3,  3),  # TAC^GTA                 blunt
    "SrfI":      ("GCCCGGGC",     4,  4),  # GCCC^GGGC               blunt (8-cutter)
    "StuI":      ("AGGCCT",       3,  3),  # AGG^CCT                 blunt
    "SwaI":      ("ATTTAAAT",     4,  4),  # ATTT^AAAT               blunt (8-cutter)
    "Tth111I":   ("GACNNNGTC",    4,  5),  # GACN^NNGTC              1-base 3' overhang
    "XmaI":      ("CCCGGG",       1,  5),  # C^CCGGG   / GGGCC^C     5' overhang (SmaI isoschizomer)
    "XmnI":      ("GAANNNNTTC",   5,  5),  # GAANN^NNTTC             blunt

    # ── Rare 8-cutters ─────────────────────────────────────────────────────────
    "AscI":      ("GGCGCGCC",     2,  6),  # GG^CGCGCC / CGCGCC^GG   5' overhang
    "AsiSI":     ("GCGATCGC",     5,  3),  # GCGAT^CGC / GCG^ATCGC   3' overhang

    # ── Degenerate / IUPAC recognition sequences ───────────────────────────────
    "AccI":      ("GTYRAC",       2,  4),  # GT^YRAC                 3' overhang
    "AclI":      ("AACGTT",       2,  4),  # AA^CGTT                 3' overhang
    "AfeI":      ("AGCGCT",       3,  3),  # AGC^GCT                 blunt (Eco47III isoschizomer)
    "AflII":     ("CTTAAG",       1,  5),  # C^TTAAG                 MfeI-compatible ends
    "AflIII":    ("ACRYGT",       1,  5),  # A^CRYGT                 MluI-compatible ends
    "AgeI":      ("ACCGGT",       1,  5),  # A^CCGGT   / TGGCC^A     5' overhang
    "AhdI":      ("GACNNNNNGTC",  6,  5),  # GACNNNN^NGTC            1-base 3' overhang
    "AluI":      ("AGCT",         2,  2),  # AG^CT                   blunt (4-cutter)
    "ApaI":      ("GGGCCC",       5,  1),  # GGGCC^C   / G^GGCCC     3' overhang
    "ApaLI":     ("GTGCAC",       1,  5),  # G^TGCAC                 SphI-compatible ends
    "ApoI":      ("RAATTY",       1,  5),  # R^AATTY                 EcoRI isoschizomer (degenerate)
    "AatII":     ("GACGTC",       5,  1),  # GACGT^C   / G^ACGTC     3' overhang
    "BaeGI":     ("GKGCMC",       5,  1),  # GKGCM^C   / G^KGCMC     3' overhang
    "BglI":      ("GCCNNNNNGGC",  7,  4),  # GCCNNNN^NGGC            3' overhang
    "BmgBI":     ("CACGTC",       3,  3),  # CAC^GTC                 blunt
    "BsaAI":     ("YACGTR",       3,  3),  # YAC^GTR                 blunt
    "BsaBI":     ("GATNNNNATC",   5,  5),  # GATN4^ATC               blunt
    "BsaHI":     ("GRCGYC",       2,  4),  # GR^CGYC                 3' overhang
    "BsaWI":     ("WCCGGW",       1,  5),  # W^CCGGW                 5' overhang
    "BseYI":     ("CCCAGC",       5,  1),  # CCCAG^C   / C^CCAGC     3' overhang
    "BsiEI":     ("CGRYCG",       4,  2),  # CGRY^CG                 3' overhang
    "BsiHKAI":   ("GWGCWC",       5,  1),  # GWGCW^C                 3' overhang
    "BsrFI":     ("RCCGGY",       1,  5),  # R^CCGGY                 5' overhang
    "Bsp1286I":  ("GDGCHC",       5,  1),  # GDGCH^C                 3' overhang
    "BspHI":     ("TCATGA",       1,  5),  # T^CATGA                 NcoI-compatible ends
    "BsrI":      ("ACTGG",        1,  1),  # ACT^GG                  degenerate Type IIS-like
    "BstAPI":    ("GCANNNNNTGC",  7,  5),  # GCANNNN^NTGC            3' overhang
    "BstNI":     ("CCWGG",        2,  3),  # CC^WGG    / WGG^CC      3' overhang
    "BstUI":     ("CGCG",         2,  2),  # CG^CG                   blunt (4-cutter; methylation-sensitive)
    "BstZ17I":   ("GTATAC",       3,  3),  # GTA^TAC                 blunt
    "BtgI":      ("CCRYGG",       1,  5),  # C^CRYGG                 5' overhang
    "Cac8I":     ("GCNNGC",       3,  3),  # GCN^NGC                 blunt
    "CviAII":    ("CATG",         1,  3),  # C^ATG                   NcoI-compatible ends (4-cutter)
    "CviQI":     ("GTAC",         1,  3),  # G^TAC                   KpnI subset (4-cutter)
    "DpnI":      ("GATC",         2,  2),  # GA^TC                   blunt; cuts only methylated
    "DpnII":     ("GATC",         0,  4),  # ^GATC     / GATC^       5' overhang (4-cutter)
    "DrdI":      ("GACNNNNNNGTC", 7,  5),  # GACNNNNN^NGTC           3' overhang
    "EcoO109I":  ("RGGNCCY",      2,  5),  # RG^GNCCY                5' overhang
    "HphI":      ("GGTGA",        8, 12),  # GGTGA(8/7) downstream   Type IIS
    "KasI":      ("GGCGCC",       1,  5),  # G^GCGCC                 5' overhang (NarI isoschizomer)
    "MboI":      ("GATC",         0,  4),  # ^GATC                   DpnII isoschizomer
    "MboII":     ("GAAGA",        8, 12),  # GAAGA(8/7) downstream   Type IIS
    "MlyI":      ("GAGTC",       10, 10),  # GAGTC(5/5) downstream   blunt, Type IIS
    "MmeI":      ("TCCRAC",      20, 18),  # TCCRAC(20/18)           Type IIS far-cutter
    "MspA1I":    ("CMGCKG",       3,  3),  # CMG^CKG                 blunt
    "NgoMIV":    ("GCCGGC",       1,  5),  # G^CCGGC                 5' overhang (EagI-compatible)
    "NmeAIII":   ("GCCGAG",      21, 19),  # GCCGAG(15/13)           Type IIS far-cutter
    "PflMI":     ("CCANNNNNTGG",  7,  4),  # CCANN4^NTGG             3' overhang (BstXI isoschizomer)
    "PspOMI":    ("GGGCCC",       1,  5),  # G^GGCCC                 ApaI isoschizomer (5' overhang)
    "Sau3AI":    ("GATC",         0,  4),  # ^GATC                   BamHI-compatible ends (4-cutter)
    "SfcI":      ("CTRYAG",       1,  5),  # C^TRYAG                 5' overhang
    "SspI":      ("AATATT",       3,  3),  # AAT^ATT                 blunt
    "TaqI":      ("TCGA",         1,  3),  # T^CGA     / CGT^A       5' overhang (heat-stable)
    "Van91I":    ("CCANNNNNTGG",  7,  4),  # PflMI isoschizomer
    "ZraI":      ("GACGTC",       3,  3),  # GAC^GTC                 blunt (AatII-related)

    # ── Type IIS — cut outside recognition sequence ────────────────────────────
    # fwd/rev positions are still offsets from start of recognition seq.
    # For an n-bp recognition sequence cutting d1/d2 downstream:
    #   fwd = n + d1,  rev = n + d2
    "BaeI":      ("ACNNNNGTAYYC",-10, -7), # cuts 10/15 upstream (negative = upstream)
    "BbsI":      ("GAAGAC",       8, 12),  # GAAGAC(2/6)  BpiI isoschizomer
    "BcoDI":     ("GTCTC",        6, 10),  # GTCTC(1/5)   BsaI 5-bp variant
    "BceAI":     ("ACGGC",        9, 11),  # ACGGC(4/6)
    "BciVI":     ("GTATCC",      12,  6),  # GTATCC(6/5)  asymmetric
    "BfuAI":     ("ACCTGC",      10, 14),  # ACCTGC(4/8)  BspMI isoschizomer
    "BmrI":      ("ACTGGN",       6,  5),  # cuts 1 past end of recog
    "BpiI":      ("GAAGAC",       8, 12),  # BbsI isoschizomer
    "BsaI":      ("GGTCTC",       7, 11),  # GGTCTC(1/5)  Golden Gate workhorse
    "BsaXI":     ("ACNNNNNCTCC",  3, 12),  # cuts upstream and downstream (unusual)
    "BsbI":      ("CAACAC",      17, 15),  # far-cutter
    "BseJI":     ("GAAGAC",       8, 12),  # BbsI isoschizomer
    "BseLI":     ("CCNNNNNNNGG",  7,  4),  # 3' overhang
    "BseMII":    ("CTCAG",       10,  8),  # CTCAG(5/3)
    "BseRI":     ("GAGGAG",      10, 14),  # GAGGAG(10/8)
    "BsgI":      ("GTGCAG",      22, 20),  # GTGCAG(16/14) far-cutter
    "BslI":      ("CCNNNNNNNGG",  7,  4),  # 3' overhang (BseLI variant)
    "BsmAI":     ("GTCTC",        6, 10),  # GTCTC(1/5)   BsaI isoschizomer (5-bp)
    "BsmBI":     ("CGTCTC",       7, 11),  # CGTCTC(1/5)  Esp3I isoschizomer
    "BsmFI":     ("GGGAC",       15, 19),  # GGGAC(10/14)
    "BsmI":      ("GAATGC",       7,  1),  # GAATGC(1/-1) asymmetric Type IIS
    "BspLU11III":("ACATGT",       1,  5),  # PciI isoschizomer (not strictly IIS)
    "BspMI":     ("ACCTGC",      10, 14),  # ACCTGC(4/8)
    "BspQI":     ("GCTCTTC",      8, 11),  # SapI isoschizomer
    "BspTNI":    ("GGTCTC",       7, 11),  # BsaI isoschizomer
    "BsrBI":     ("CCGCTC",       3,  3),  # cuts within recog (special case)
    "BsrDI":     ("GCAATG",       8,  6),  # GCAATG(2/0)
    "BssSI":     ("CACGAG",      -5, -1),  # cuts upstream of recognition
    "BtgZI":     ("GCGATG",      16, 20),  # GCGATG(10/14)
    "BtsCI":     ("GGATG",        5,  3),  # GGATG(2/0)
    "BtsI":      ("GCAGTG",       8,  6),  # GCAGTG(2/0)
    "BtsImutI":  ("CAGTG",        5,  3),  # BtsCI variant
    "EarI":      ("CTCTTC",      10,  6),  # CTCTTC(4/1)  SapI-related
    "Esp3I":     ("CGTCTC",       7, 11),  # BsmBI isoschizomer
    "PaqCI":     ("CACCTGC",     11, 15),  # CACCTGC(4/8)
    "SapI":      ("GCTCTTC",      8, 11),  # GCTCTTC(1/4)
    "BsmBI-v2":  ("CGTCTC",       7, 11),  # v2/HF variant

    # ── High-Fidelity (HF) and v2 variants — same recognition/cut as canonical ─
    "AgeI-HF":   ("ACCGGT",       1,  5),
    "BamHI-HF":  ("GGATCC",       1,  5),
    "BclI-HF":   ("TGATCA",       1,  5),
    "BmtI":      ("GCTAGC",       1,  5),  # NheI isoschizomer / HF variant
    "BsiWI-HF":  ("CGTACG",       1,  5),
    "BsrFI-v2":  ("RCCGGY",       1,  5),
    "BsrGI-HF":  ("TGTACA",       1,  5),
    "BssSI-v2":  ("CACGAG",      -5, -1),
    "BstEII-HF": ("GGTNACC",      1,  5),
    "BstZ17I-HF":("GTATAC",       3,  3),
    "DraIII-HF": ("CACNNNGTG",    6,  3),
    "EcoRI-HF":  ("GAATTC",       1,  5),
    "EcoRV-HF":  ("GATATC",       3,  3),
    "HindIII-HF":("AAGCTT",       1,  5),
    "KpnI-HF":   ("GGTACC",       5,  1),
    "MfeI-HF":   ("CAATTG",       1,  5),
    "MluI-HF":   ("ACGCGT",       1,  5),
    "MunI-HF":   ("CAATTG",       1,  5),
    "NcoI-HF":   ("CCATGG",       1,  5),
    "NheI-HF":   ("GCTAGC",       1,  5),
    "NotI-HF":   ("GCGGCCGC",     2,  6),
    "NruI-HF":   ("TCGCGA",       3,  3),
    "NsiI-HF":   ("ATGCAT",       5,  1),
    "PstI-HF":   ("CTGCAG",       5,  1),
    "PvuI-HF":   ("CGATCG",       4,  2),
    "PvuII-HF":  ("CAGCTG",       3,  3),
    "SacI-HF":   ("GAGCTC",       5,  1),
    "SalI-HF":   ("GTCGAC",       1,  5),
    "SbfI-HF":   ("CCTGCAGG",     6,  2),
    "ScaI-HF":   ("AGTACT",       3,  3),
    "SpeI-HF":   ("ACTAGT",       1,  5),
    "SphI-HF":   ("GCATGC",       5,  1),
    "TaqI-v2":   ("TCGA",         1,  3),
    "XhoI-HF":   ("CTCGAG",       1,  5),
}

# Flat recognition-sequence-only dict derived from _NEB_ENZYMES (used by
# _scan_restriction_sites and _RESTR_COLOR below).
_RESTRICTION_SITES: dict[str, str] = {
    name: seq for name, (seq, _, _) in _NEB_ENZYMES.items()
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


_IUPAC_RE: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "W": "[AT]", "S": "[CG]",
    "M": "[AC]", "K": "[GT]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


_PATTERN_CACHE: dict[str, "re.Pattern[str]"] = {}

def _iupac_pattern(site: str) -> "re.Pattern[str]":
    if site not in _PATTERN_CACHE:
        _PATTERN_CACHE[site] = re.compile(
            "".join(_IUPAC_RE.get(c, c) for c in site.upper())
        )
    return _PATTERN_CACHE[site]


_IUPAC_COMP = str.maketrans(
    "ACGTRYWSMKBDHVN",
    "TGCAYRWSKMVHDBN",
)

def _rc(seq: str) -> str:
    return seq.upper().translate(_IUPAC_COMP)[::-1]


# Pre-built per-enzyme scan records: immutable derived values (compiled pattern,
# palindrome flag, RC pattern, color) computed once at import time rather than
# on every _scan_restriction_sites call. Iterating the pre-built list avoids
# ~200 per-call _iupac_pattern + _rc + dict-lookup + len calls per scan.
#
# Each entry is a tuple:
#   (name, site, site_len, fwd_cut, rev_cut, color, pat, is_palindrome, rc_pat)
# rc_pat is None for palindromic enzymes (no reverse scan needed).
_SCAN_CATALOG: "list[tuple]" = []


def _rebuild_scan_catalog() -> None:
    """(Re)populate `_SCAN_CATALOG` from `_NEB_ENZYMES`. Called at import
    time; also exposed so tests / future catalog edits can refresh it."""
    _SCAN_CATALOG.clear()
    for name, (site, fwd_cut, rev_cut) in _NEB_ENZYMES.items():
        site_u  = site.upper()
        pat     = _iupac_pattern(site_u)
        rc_site = _rc(site_u)
        is_pal  = (rc_site == site_u)
        rc_pat  = None if is_pal else _iupac_pattern(rc_site)
        _SCAN_CATALOG.append((
            name, site_u, len(site_u), fwd_cut, rev_cut,
            _RESTR_COLOR[name], pat, is_pal, rc_pat,
        ))

_rebuild_scan_catalog()


def _scan_restriction_sites(
    seq: str,
    min_recognition_len: int = 6,
    unique_only: bool = True,
) -> list[dict]:
    """Scan both strands; return resite + recut dicts for every hit.

    resite — the recognition sequence span (colored bar)
    recut  — the cut position (single-bp marker: ↓ above or ↑ below DNA)

    min_recognition_len — skip enzymes whose recognition site is shorter than this
                          (default 6 to reduce noise from 4-cutters)
    unique_only         — if True, only include enzymes that cut exactly once
                          (forward + reverse strand combined; default True)
    """
    seq_u = seq.upper()
    n = len(seq_u)
    # Per-enzyme results collected first so we can filter to unique cutters
    by_enzyme: dict[str, list[dict]] = {}
    seen: set[tuple[str, int, int]] = set()   # deduplicate palindromes

    for entry in _SCAN_CATALOG:
        name, site, site_len, fwd_cut, rev_cut, color, pat, is_palindrome, rc_pat = entry
        if site_len < min_recognition_len:
            continue
        hits: list[dict] = []

        # Forward strand scan
        for m in pat.finditer(seq_u):
            p = m.start()
            key = (name, p, 1)
            if key in seen:
                continue
            seen.add(key)
            # ext_cut_bp: absolute cut position when cut falls outside recognition
            _ext = (p + fwd_cut) if (fwd_cut <= 0 or fwd_cut >= site_len) else None
            hits.append({
                "type":       "resite",
                "start":      p,
                "end":        p + site_len,
                "strand":     1,
                "color":      color,
                "label":      name,
                "cut_col":    fwd_cut if 0 < fwd_cut < site_len else None,
                "ext_cut_bp": _ext,
            })
            cut_bp = min(p + fwd_cut, n - 1)
            hits.append({
                "type":   "recut",
                "start":  cut_bp,
                "end":    cut_bp + 1,
                "strand": 1,
                "color":  color,
                "label":  name,
            })

        # Reverse strand handling — uses precomputed is_palindrome and rc_pat
        # from _SCAN_CATALOG (no per-call _rc / _iupac_pattern work).
        if not is_palindrome:
            # Non-palindromic: scan for RC on forward strand to find
            # reverse-strand binding sites at their correct positions.
            for m in rc_pat.finditer(seq_u):
                p = m.start()
                key = (name, p, -1)
                if key in seen:
                    continue
                seen.add(key)
                # Cut column within the bar: enzyme's fwd_cut mapped to
                # the reversed orientation displayed on the forward strand
                rev_cut_col = site_len - 1 - fwd_cut
                _top_cut_bp = p + site_len - 1 - rev_cut   # top-strand cut in fwd coords
                _top_cut_outside = (_top_cut_bp < p or _top_cut_bp >= p + site_len)
                hits.append({
                    "type":       "resite",
                    "start":      p,
                    "end":        p + site_len,
                    "strand":     -1,
                    "color":      color,
                    "label":      name,
                    "cut_col":    rev_cut_col if 0 <= rev_cut_col < site_len else None,
                    "ext_cut_bp": _top_cut_bp if _top_cut_outside else None,
                })
                # Bottom-strand cut (enzyme's fwd_cut mapped to fwd coords)
                cut_bp = p + site_len - 1 - fwd_cut
                cut_bp = max(0, min(cut_bp, n - 1))
                hits.append({
                    "type":   "recut",
                    "start":  cut_bp,
                    "end":    cut_bp + 1,
                    "strand": -1,
                    "color":  color,
                    "label":  name,
                })

        if hits:
            by_enzyme[name] = hits

    feats: list[dict] = []
    placed: set[tuple[int, int]] = set()   # (start, end) of resites already shown
    for name, hits in by_enzyme.items():
        # Count recognition-sequence hits (resite only) across both strands
        if unique_only:
            n_sites = sum(1 for h in hits if h["type"] == "resite")
            if n_sites != 1:
                continue
        # Skip isoschizomers / HF-variants that land on an already-placed site
        positions = {(h["start"], h["end"]) for h in hits if h["type"] == "resite"}
        if positions & placed:
            continue
        placed |= positions
        feats.extend(hits)
    return feats


def _assign_chunk_features(
    chunk_feats: list[dict], chunk_start: int, chunk_end: int
) -> tuple[list[list[dict]], list[list[dict]]]:
    """
    Forward-strand features always go above DNA; reverse-strand always below.
    Within each group, overlapping features are stacked into greedy lanes
    (each lane is a list of non-overlapping features).  Capped at 3 lanes/side.
    """
    def _greedy_lanes(feats: list[dict]) -> list[list[dict]]:
        sorted_f = sorted(feats, key=lambda f: max(f["start"], chunk_start))
        lanes: list[list[dict]] = []
        lane_ends: list[int]    = []
        for f in sorted_f:
            bar_s = max(f["start"], chunk_start)
            bar_e = min(f["end"],   chunk_end)
            if bar_e - bar_s <= 0:
                continue
            placed = False
            for i, end in enumerate(lane_ends):
                if bar_s >= end:
                    lanes[i].append(f)
                    lane_ends[i] = bar_e
                    placed = True
                    break
            if not placed:
                lanes.append([f])
                lane_ends.append(bar_e)
        return lanes[:3]

    fwd = [f for f in chunk_feats if f["strand"] >= 0]
    rev = [f for f in chunk_feats if f["strand"] <  0]
    return _greedy_lanes(fwd), _greedy_lanes(rev)


def _render_feature_row_pair(
    result: "Text",
    feats: list[dict],
    chunk_start: int,
    chunk_end: int,
    prefix_w: int,
    is_below_dna: bool,
    show_connectors: bool,
    flip_label_bar: bool = False,
    single_row: bool = False,
) -> None:
    """
    Append one label row + optional connector row + one feature-bar row to result.
    For above-DNA: label / [connector] / bar.
    For below-DNA: bar / [connector] / label.

    Feature bars use ▒ (medium-shade / dither) as the fill glyph so the bar
    reads as a single flat coloured tone rather than a pattern of distinct dots.
    This is DELIBERATELY different from PlasmidMap, which still uses braille
    sub-character rendering (U+2800–U+28FF) for the circular-map feature arcs —
    only the sequence panel switched off braille.

    When single_row=True (RE lanes), collapse to one content row: the cut arrow is
    placed just outside the bracket (left of '(' above DNA, right of ')' below) only
    if that cell is a space, so it never overwrites another character.
    Multiple non-overlapping features share the same pair of rows horizontally.
    """
    content_w = chunk_end - chunk_start
    label_arr: list[tuple[str, str]] = [(" ", "")] * content_w
    bar_arr:   list[tuple[str, str]] = [(" ", "")] * content_w
    conn_arr:  list[tuple[str, str]] = [(" ", "")] * content_w

    for f in feats:
        bar_s = max(f["start"], chunk_start) - chunk_start
        bar_e = min(f["end"],   chunk_end)   - chunk_start
        bar_len = bar_e - bar_s
        if bar_len <= 0:
            continue
        starts_here = f["start"] >= chunk_start
        ends_here   = f["end"]   <= chunk_end
        strand      = f["strand"]
        color       = f["color"]
        label       = f.get("label", f.get("type", ""))

        feat_type = f.get("type", "")

        if feat_type == "resite":
            # ── Parenthesis-style RE site: ( EnzymeName ) ─────────────────
            # Bar row:   ( EnzymeName )  — bold white name, colored parens
            # Label row: ↓ or ↑ at the intra-site cut column (if cut is
            #            within the recognition site; Type IIS cuts elsewhere)
            cut_col = f.get("cut_col")   # 0-based offset from f["start"], or None

            # Opening / closing parens
            if starts_here and bar_len >= 1:
                bar_arr[bar_s] = ("(", color)
            if ends_here and bar_len >= 1:
                bar_arr[bar_s + bar_len - 1] = (")", color)

            # Enzyme name: bold white, centered in the interior
            interior_start = (1 if starts_here else 0)
            interior_end   = (bar_len - 1 if ends_here else bar_len)
            interior_len   = interior_end - interior_start
            if interior_len > 0:
                name_str  = label[:interior_len]
                name_pad  = interior_len - len(name_str)
                name_lpad = name_pad // 2
                for j, ch in enumerate(name_str):
                    pos = bar_s + interior_start + name_lpad + j
                    if 0 <= pos < content_w:
                        bar_arr[pos] = (ch, "bold white")

            # Cut marker in the label row so it doesn't obscure the name
            cut_ch = "↑" if is_below_dna else "↓"
            if cut_col is not None:
                visible_offset = cut_col - max(0, chunk_start - f["start"])
                cut_pos = bar_s + visible_offset
                if 0 <= cut_pos < content_w:
                    label_arr[cut_pos] = (cut_ch, "bold " + color)

            # Type IIS: dashed bridge from recognition bar to cut site
            ext_cut_bp = f.get("ext_cut_bp")
            if ext_cut_bp is not None and chunk_start <= ext_cut_bp < chunk_end:
                cut_abs = ext_cut_bp - chunk_start
                if 0 <= cut_abs < content_w:
                    label_arr[cut_abs] = (cut_ch, "bold " + color)
                # Bridge in bar_arr: from recognition end rightward, or from
                # cut leftward to recognition start (upstream cutters)
                if cut_abs >= bar_s + bar_len:       # downstream cut
                    for j in range(bar_s + bar_len, cut_abs):
                        if 0 <= j < content_w and bar_arr[j][0] == " ":
                            bar_arr[j] = ("╌", color)
                elif cut_abs < bar_s:                # upstream cut
                    for j in range(cut_abs + 1, bar_s):
                        if 0 <= j < content_w and bar_arr[j][0] == " ":
                            bar_arr[j] = ("╌", color)

            # Connector tick at midpoint
            mid = bar_s + bar_len // 2
            if 0 <= mid < content_w:
                conn_arr[mid] = ("┊", color)
            continue   # skip the regular label/bar/conn logic below

        elif feat_type == "recut":
            continue   # cut position is rendered inside the resite bar; skip here

        # ── Standard feature (non-RE) ──────────────────────────────────────

        # Label (centered in feature span)
        lbl = label[:bar_len]
        pad = bar_len - len(lbl)
        pl  = pad // 2
        lbl_str = " " * pl + lbl + " " * (pad - pl)
        for i, ch in enumerate(lbl_str):
            if 0 <= bar_s + i < content_w:
                label_arr[bar_s + i] = (ch, color)

        # Dithered feature bar — medium-shade fill + directional arrowhead at
        # the feature end (strand direction preserved). 1 bp features get a
        # triangle pointing toward the DNA row.
        if bar_len == 1:
            bar_str = "▲" if is_below_dna else "▼"
        elif strand >= 0:
            bar_str = "▒" * (bar_len - (1 if ends_here   else 0)) + ("▶" if ends_here   else "")
        else:
            bar_str = ("◀" if starts_here else "") + "▒" * (bar_len - (1 if starts_here else 0))
        for i, ch in enumerate(bar_str):
            if 0 <= bar_s + i < content_w:
                bar_arr[bar_s + i] = (ch, color)

        # Connector tick at midpoint of feature span
        mid = bar_s + bar_len // 2
        if 0 <= mid < content_w:
            conn_arr[mid] = ("┊", color)

    def _write_arr(arr: list[tuple[str, str]]) -> None:
        result.append(" " * prefix_w)
        run: list[str] = []
        sty = ""
        for ch, s in arr:
            if s == sty:
                run.append(ch)
            else:
                if run:
                    result.append("".join(run), style=sty)
                run = [ch]
                sty = s
        if run:
            result.append("".join(run), style=sty)
        result.append("\n")

    if single_row:
        # Place cut arrow adjacent to the bracket — never overlapping a name char.
        for f in feats:
            if f.get("type") != "resite":
                continue
            bar_s  = max(f["start"], chunk_start) - chunk_start
            bar_e  = min(f["end"],   chunk_end)   - chunk_start
            cut_ch = "↑" if is_below_dna else "↓"
            color  = f["color"]
            if not is_below_dna:
                # Above DNA: try left of opening paren, then right of closing paren
                if bar_s > 0 and bar_arr[bar_s - 1][0] == " ":
                    bar_arr[bar_s - 1] = (cut_ch, "bold " + color)
                elif bar_e < content_w and bar_arr[bar_e][0] == " ":
                    bar_arr[bar_e] = (cut_ch, "bold " + color)
            else:
                # Below DNA: try right of closing paren, then left of opening paren
                if bar_e < content_w and bar_arr[bar_e][0] == " ":
                    bar_arr[bar_e] = (cut_ch, "bold " + color)
                elif bar_s > 0 and bar_arr[bar_s - 1][0] == " ":
                    bar_arr[bar_s - 1] = (cut_ch, "bold " + color)
            # Type IIS: ext_cut_bp arrow goes into bar_arr at the actual cut position
            ext_cut_bp = f.get("ext_cut_bp")
            if ext_cut_bp is not None and chunk_start <= ext_cut_bp < chunk_end:
                cut_abs = ext_cut_bp - chunk_start
                if 0 <= cut_abs < content_w and bar_arr[cut_abs][0] == " ":
                    bar_arr[cut_abs] = (cut_ch, "bold " + color)
        if not is_below_dna:
            _write_arr(bar_arr)
            if show_connectors:
                _write_arr(conn_arr)
        else:
            if show_connectors:
                _write_arr(conn_arr)
            _write_arr(bar_arr)
        return

    first_arr  = bar_arr   if flip_label_bar else label_arr
    second_arr = label_arr if flip_label_bar else bar_arr
    if not is_below_dna:
        _write_arr(first_arr)
        if show_connectors:
            _write_arr(conn_arr)
        _write_arr(second_arr)
    else:
        _write_arr(second_arr)
        if show_connectors:
            _write_arr(conn_arr)
        _write_arr(first_arr)


def _chunk_lane_groups(
    chunk_feats: list[dict], chunk_start: int, chunk_end: int,
) -> "tuple[list, list, list, list, list, list]":
    """Split chunk features into separate lane groups for rendering order.

    Returns (re_above, onebp_above, reg_above, reg_below, onebp_below, re_below).
    Rendering order top→bottom:
      re_above → onebp_above → reg_above → DNA → reg_below → onebp_below → re_below
    Multi-bp regular features are closest to DNA; 1bp features next; RE sites farthest.
    """
    resites, onebp, multibp = [], [], []
    for f in chunk_feats:
        if f.get("type") == "resite":
            resites.append(f)
        elif f["end"] - f["start"] == 1:
            onebp.append(f)
        else:
            multibp.append(f)
    reg_above,   reg_below   = _assign_chunk_features(multibp, chunk_start, chunk_end)
    onebp_above, onebp_below = _assign_chunk_features(onebp,   chunk_start, chunk_end)
    re_above,    re_below    = _assign_chunk_features(resites, chunk_start, chunk_end)
    return re_above, onebp_above, reg_above, reg_below, onebp_below, re_below


# Per-(seq_id, feats_id) cache for expensive inputs of _build_seq_text that
# only depend on sequence and features, not on cursor/selection/line_width.
# Cache holds (styles_list, annot_feats_sorted). Invalidated by id — lists
# are reassigned on load, never mutated in place (see CLAUDE.md).
_BUILD_SEQ_CACHE: dict = {}

def _build_seq_inputs(seq: str, feats: list[tuple]) -> tuple:
    """Return (styles, annot_feats) for a given sequence/feature pair,
    memoised by identity. Both outputs are expensive on large plasmids and
    independent of cursor/selection state."""
    key = (id(seq), id(feats), len(seq), len(feats))
    hit = _BUILD_SEQ_CACHE.get(key)
    if hit is not None:
        return hit
    n = len(seq)
    styles = ["color(252)"] * n
    for f in reversed(feats):          # reversed so first feature wins
        col = f["color"]
        for i in range(f["start"], min(f["end"], n)):
            styles[i] = col
    annot_feats = sorted(
        [f for f in feats if f.get("type") not in ("site", "recut")],
        key=lambda f: -(f["end"] - f["start"]),
    )
    # Cap the cache at 4 entries (one active + a few stale) — we're keying
    # on id() so size stays tiny; this is just belt-and-braces.
    if len(_BUILD_SEQ_CACHE) >= 4:
        _BUILD_SEQ_CACHE.clear()
    _BUILD_SEQ_CACHE[key] = (styles, annot_feats)
    return styles, annot_feats


def _build_seq_text(seq: str, feats: list[dict], line_width: int = 60,
                    sel_range: "tuple[int,int] | None" = None,
                    user_sel:  "tuple[int,int] | None" = None,
                    cursor_pos: int = -1,
                    show_connectors: bool = False,
                    re_highlight: "dict | None" = None) -> Text:
    """Rich Text of the sequence with per-position feature coloring.

    sel_range    — feature highlight: bold + underline on feature bases
    user_sel     — shift-click selection: subtle background, used by edit dialog
    cursor_pos   — click cursor: reverse-video highlight on base at cursor_pos
    re_highlight — dict with keys: start, end, fwd_cut_bp, rev_cut_bp, color, name
                   When set, highlights the recognition bases on both strands
                   and marks cut positions with reverse-video.

    Rendering order (closest to DNA first):
      RE sites (far) → regular feature bars (close) → DNA → regular (close) → RE (far)
    """
    n     = len(seq)
    num_w = len(str(n)) if n else 1    # minimum digits needed for line numbers
    styles, annot_feats = _build_seq_inputs(seq, feats)

    sel_s  = sel_range[0] if sel_range else -1
    sel_e  = sel_range[1] if sel_range else -1
    usr_s  = user_sel[0]  if user_sel  else -1
    usr_e  = user_sel[1]  if user_sel  else -1

    # RE highlight ranges
    reh_s       = re_highlight["start"]      if re_highlight else -1
    reh_e       = re_highlight["end"]        if re_highlight else -1
    reh_fwd_cut = re_highlight["fwd_cut_bp"] if re_highlight else -1
    reh_rev_cut = re_highlight["rev_cut_bp"] if re_highlight else -1
    reh_color   = re_highlight["color"]      if re_highlight else ""

    seq_upper = seq.upper()
    _COMP     = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")   # base complement
    result    = Text(no_wrap=True, overflow="crop")

    for chunk_start in range(0, n, line_width):
        chunk_end = min(chunk_start + line_width, n)

        # ── Assign features to lane groups ──
        chunk_feats = [
            f for f in annot_feats
            if f["start"] < chunk_end and f["end"] > chunk_start
        ]
        re_above, onebp_above, reg_above, reg_below, onebp_below, re_below = (
            _chunk_lane_groups(chunk_feats, chunk_start, chunk_end)
        )

        # ── Rows ABOVE DNA (far → close): RE → 1bp → multi-bp ──
        for lane in re_above:
            _render_feature_row_pair(result, lane, chunk_start, chunk_end,
                                     num_w + 2, False, show_connectors,
                                     flip_label_bar=True)
        for lane in onebp_above:
            _render_feature_row_pair(result, lane, chunk_start, chunk_end,
                                     num_w + 2, False, show_connectors)
        for lane in reg_above:
            _render_feature_row_pair(result, lane, chunk_start, chunk_end,
                                     num_w + 2, False, show_connectors)

        # ── Double-stranded DNA block ─────────────────────────────────────
        # The cursor is shown as a reverse-video highlight on the base IN
        # PLACE rather than as an inserted glyph, so column counts match.

        def _strand_chars(bases: "list[str]") -> None:
            """Append base chars with RLE styling into result."""
            run: list[str] = []
            sty = ""
            for ch, s in bases:
                if s == sty:
                    run.append(ch)
                else:
                    if run:
                        result.append("".join(run), style=sty)
                    run, sty = [ch], s
            if run:
                result.append("".join(run), style=sty)

        # Perf: translate the whole chunk once instead of per-base.
        chunk_fwd = seq_upper[chunk_start:chunk_end]
        chunk_rev = chunk_fwd.translate(_COMP)
        chunk_len = chunk_end - chunk_start
        fwd_bases: list[tuple[str, str]] = []
        rc_bases:  list[tuple[str, str]] = []
        for j in range(chunk_len):
            i      = chunk_start + j
            base   = styles[i]
            in_usr = (usr_s <= i < usr_e)
            in_sel = (sel_s <= i < sel_e)
            in_re  = (reh_s <= i < reh_e)
            is_cur = (cursor_pos == i)

            if is_cur:
                fwd_sty = "reverse bold white"
                rev_sty = fwd_sty
            elif in_re:
                # Entire recognition region: white background, black text
                fwd_sty = f"reverse bold {reh_color}"
                rev_sty = f"reverse bold {reh_color}"
            elif in_usr:
                fwd_sty = base + " on color(237)"
                rev_sty = fwd_sty
            elif in_sel:
                fwd_sty = "bold underline " + base
                rev_sty = fwd_sty
            else:
                fwd_sty = base
                rev_sty = base
            fwd_bases.append((chunk_fwd[j], fwd_sty))
            rc_bases.append( (chunk_rev[j], rev_sty))

        # Forward strand
        result.append(f"{chunk_start + 1:>{num_w}}  ", style="color(245)")
        _strand_chars(fwd_bases)
        result.append("\n")

        # Reverse-complement strand (aligned column-for-column)
        result.append(" " * (num_w + 2), style="color(245)")
        _strand_chars(rc_bases)
        result.append("\n")

        # ── Rows BELOW DNA (close → far): multi-bp → 1bp → RE ──
        for lane in reg_below:
            _render_feature_row_pair(result, lane, chunk_start, chunk_end,
                                     num_w + 2, True, show_connectors)
        for lane in onebp_below:
            _render_feature_row_pair(result, lane, chunk_start, chunk_end,
                                     num_w + 2, True, show_connectors)
        for lane in re_below:
            _render_feature_row_pair(result, lane, chunk_start, chunk_end,
                                     num_w + 2, True, show_connectors,
                                     flip_label_bar=True)

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

def _copy_to_clipboard_osc52(text: str) -> bool:
    """Copy text via OSC 52 escape sequence — works in Windows Terminal, iTerm2, most modern terminals."""
    import base64
    encoded = base64.b64encode(text.encode()).decode()
    seq = f"\033]52;c;{encoded}\007"
    try:
        with open("/dev/tty", "w") as tty:
            tty.write(seq)
            tty.flush()
        return True
    except Exception:
        return False


def _translate_cds(full_seq: str, start: int, end: int, strand: int) -> str:
    """Translate a CDS region to single-letter AA string (stop codon → *).

    Uses _IUPAC_COMP for the reverse-complement step so IUPAC ambiguity codes
    (N, R, Y, etc.) are handled correctly. An earlier version used a bare
    ACGT-only maketrans which would silently pass degenerate bases through
    unchanged — producing wrong codons and silent mistranslation.
    """
    sub = full_seq[start:end].upper()
    if strand == -1:
        sub = sub.translate(_IUPAC_COMP)[::-1]
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


# ── External annotation (pLannotate) ──────────────────────────────────────────
#
# Integration with pLannotate (https://github.com/mmcguffi/pLannotate) as an
# OPTIONAL runtime dependency. pLannotate is GPL-3 so we only call its CLI
# as a subprocess — we never `import plannotate` (which would arguably create
# a combined work under GPL).
#
# If pLannotate is not installed (shutil.which returns None), the UI entry
# points notify the user with install instructions and no error propagates.
#
# Install path (from the pLannotate README):
#     conda create -n plannotate -c conda-forge -c bioconda plannotate
#     conda activate plannotate
#     plannotate setupdb        # ~500 MB database, one-time
#
# pLannotate refuses inputs larger than 50 kb (its MAX_PLAS_SIZE constant),
# so we preflight that too and give a specific error.

class PlannotateError(Exception):
    """Base class for pLannotate errors carrying a user-facing message."""
    def __init__(self, user_msg: str, detail: str = ""):
        super().__init__(user_msg if not detail else f"{user_msg}: {detail}")
        self.user_msg = user_msg
        self.detail   = detail

class PlannotateNotInstalled(PlannotateError): pass
class PlannotateMissingDb(PlannotateError):   pass
class PlannotateTooLarge(PlannotateError):    pass
class PlannotateFailed(PlannotateError):      pass

# pLannotate's hard-coded maximum plasmid size (MAX_PLAS_SIZE in its resources).
_PLANNOTATE_MAX_BP = 50_000

# Cached availability probe — cleared by setting to None.
_PLANNOTATE_CHECK_CACHE: "dict | None" = None

def _plannotate_status() -> dict:
    """Check whether pLannotate + BLAST+ + diamond are on PATH. Cached."""
    global _PLANNOTATE_CHECK_CACHE
    if _PLANNOTATE_CHECK_CACHE is not None:
        return _PLANNOTATE_CHECK_CACHE
    import shutil
    status = {
        "installed": shutil.which("plannotate") is not None,
        "blast":     shutil.which("blastn")     is not None,
        "diamond":   shutil.which("diamond")    is not None,
    }
    status["ready"] = all((status["installed"], status["blast"], status["diamond"]))
    _PLANNOTATE_CHECK_CACHE = status
    return status

def _plannotate_install_hint() -> str:
    """User-friendly install instructions for notifications."""
    return (
        "Install via conda:  conda create -n plannotate "
        "-c conda-forge -c bioconda plannotate && "
        "conda activate plannotate && plannotate setupdb"
    )

def _run_plannotate(record, timeout: int = 180):
    """Run `plannotate batch` on a temporary copy of `record` and return the
    parsed output as a Biopython SeqRecord. Raises a PlannotateError subclass
    on any failure; callers are expected to catch and surface `err.user_msg`.

    Never imports `plannotate` — invokes the CLI as a subprocess so the GPL
    boundary stays clean.
    """
    import os
    import subprocess
    import tempfile

    status = _plannotate_status()
    if not status["installed"]:
        raise PlannotateNotInstalled(
            "pLannotate is not installed", _plannotate_install_hint()
        )
    if not status["blast"] or not status["diamond"]:
        missing = []
        if not status["blast"]:   missing.append("blastn")
        if not status["diamond"]: missing.append("diamond")
        raise PlannotateNotInstalled(
            f"pLannotate requires {' + '.join(missing)} on PATH",
            _plannotate_install_hint(),
        )

    n = len(record.seq)
    if n > _PLANNOTATE_MAX_BP:
        raise PlannotateTooLarge(
            f"pLannotate max input is {_PLANNOTATE_MAX_BP:,} bp "
            f"(this record: {n:,} bp)"
        )

    with tempfile.TemporaryDirectory(prefix="splicecraft_plan_") as tmp:
        in_path = os.path.join(tmp, "input.gb")
        with open(in_path, "w", encoding="utf-8") as fh:
            fh.write(_record_to_gb_text(record))

        try:
            result = subprocess.run(
                [
                    "plannotate", "batch",
                    "-i", in_path,
                    "-o", tmp,
                    "-f", "annotated",
                    "-s", "",           # no "_pLann" suffix
                ],
                capture_output=True, text=True, timeout=timeout,
            )
        except FileNotFoundError:
            # PATH said plannotate existed, but the binary disappeared between
            # the `which` check and the subprocess. Rare; report as not-installed.
            raise PlannotateNotInstalled(
                "pLannotate not found on PATH", _plannotate_install_hint()
            )
        except subprocess.TimeoutExpired:
            raise PlannotateFailed(
                f"pLannotate timed out after {timeout}s",
                f"input was {n:,} bp",
            )

        combined_out = (result.stdout or "") + (result.stderr or "")
        if "Databases not downloaded" in combined_out:
            raise PlannotateMissingDb(
                "pLannotate databases not installed",
                "Run: plannotate setupdb",
            )

        if result.returncode != 0:
            err_tail = combined_out[-500:].strip() or "(no output)"
            raise PlannotateFailed("pLannotate failed", err_tail)

        # Locate output GenBank. `plannotate batch` writes <file_name>.gbk; we
        # passed `-f annotated` so look for annotated.gbk first, then any .gbk
        # that isn't our input.
        out_gb = os.path.join(tmp, "annotated.gbk")
        if not os.path.exists(out_gb):
            candidates = [
                f for f in os.listdir(tmp)
                if f.endswith(".gbk") and f != "input.gb"
            ]
            if not candidates:
                raise PlannotateFailed("pLannotate produced no .gbk output")
            out_gb = os.path.join(tmp, candidates[0])

        try:
            annotated = load_genbank(out_gb)
        except Exception as exc:
            raise PlannotateFailed("could not parse pLannotate output", str(exc))

    return annotated


def _merge_plannotate_features(original, annotated):
    """Return a NEW SeqRecord with `original`'s sequence and features, plus
    any non-duplicate features from `annotated` tagged with a "pLannotate"
    note qualifier.

    - Preserves the original sequence (annotated may round-trip differently).
    - Skips feature type "source" (GenBank boilerplate).
    - Skips a pLannotate feature if (type, start, end, strand) matches an
      existing feature — avoids duplicating features the user already has.
    """
    from copy import deepcopy
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    merged = SeqRecord(
        Seq(str(original.seq)),
        id=original.id,
        name=original.name,
        description=original.description,
        annotations=dict(original.annotations),
    )
    for feat in original.features:
        merged.features.append(deepcopy(feat))

    existing_keys = {
        (f.type, int(f.location.start), int(f.location.end), f.location.strand)
        for f in original.features if f.type != "source"
    }
    n_added = 0
    for feat in annotated.features:
        if feat.type == "source":
            continue
        key = (
            feat.type,
            int(feat.location.start),
            int(feat.location.end),
            feat.location.strand,
        )
        if key in existing_keys:
            continue
        new_feat = deepcopy(feat)
        # Tag with a "pLannotate" note so users (and future loads) can tell
        # where each feature came from.
        notes = new_feat.qualifiers.get("note", [])
        if isinstance(notes, str):
            notes = [notes]
        if not any("pLannotate" in n for n in notes):
            notes = ["pLannotate"] + list(notes)
        new_feat.qualifiers["note"] = notes
        merged.features.append(new_feat)
        n_added += 1

    # Stash the count so the UI caller can surface it in a notification.
    merged._plannotate_added = n_added
    return merged


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
        Binding("left",        "rotate_cw",        "Rotate ←",      show=True),
        Binding("right",       "rotate_ccw",        "Rotate →",      show=True),
        Binding("shift+left",  "rotate_cw_lg",     "Rotate ←←",     show=False),
        Binding("shift+right", "rotate_ccw_lg",    "Rotate →→",     show=False),
        Binding("home",        "reset_origin",     "Reset",         show=False),
        Binding("comma",       "aspect_dec",       "Circle wider",   show=False),
        Binding("full_stop",   "aspect_inc",       "Circle taller",  show=False),
        Binding("v",           "toggle_map_view",  "Toggle view",    show=False),
    ]

    origin_bp:    reactive[int]   = reactive(0)
    selected_idx: reactive[int]   = reactive(-1)
    _aspect:      reactive[float] = reactive(2.0)
    _map_mode:    reactive[str]   = reactive("circular")

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
        self._feats:          list[dict] = []
        self._restr_feats:    list[dict] = []   # restriction site overlay
        self._total:          int  = 0
        self._show_connectors: bool = False

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
            # Compound / joined locations (e.g. join(100..200,300..400)) are
            # flattened to their outer bounds. Plasmid features are virtually
            # never compound (no introns), but if an imported GenBank file has
            # one, we render the full span rather than silently dropping it.
            # The log records which features were flattened for debugging.
            try:
                from Bio.SeqFeature import CompoundLocation
                if isinstance(feat.location, CompoundLocation):
                    _log.info(
                        "Flattened compound feature %s (%d..%d) to outer bounds",
                        _feat_label(feat), start, end,
                    )
            except ImportError:
                pass
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

    def action_toggle_map_view(self):
        self._map_mode = "linear" if self._map_mode == "circular" else "circular"
        self.refresh()

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
        if self._map_mode == "linear":
            idx = self._feat_at_linear(event.x, event.y)
        else:
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

    # (w, h, origin_bp, selected_idx, _aspect, n_feats, n_restr, ...) → Text
    # NOTE: do NOT rename this to `_render_cache` — Textual.Widget has its
    # own internal `_render_cache` attribute (a _RenderCache dataclass for
    # strip caching) and the two would collide, producing confusing bugs
    # when either side's cache is inspected.
    _draw_cache: "tuple | None" = None

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
               id(self._feats), id(self._restr_feats), self._map_mode,
               self._show_connectors, self.record.name)
        if self._draw_cache and self._draw_cache[0] == key:
            return self._draw_cache[1]
        result = self._draw_linear(w, h) if self._map_mode == "linear" else self._draw(w, h)
        self._draw_cache = (key, result)
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

        # ── Position ticks (inside the circle) ───────────────────────────────
        # TICK_DR_MARK  — radial inset for the ┼ graduation mark
        # TICK_DR_LABEL — radial inset for the bp number label
        # Both are negative so they sit inside the backbone ring.
        # These values scale automatically with the aspect ratio (, / . keys).
        TICK_DR_MARK  = -2
        TICK_DR_LABEL = -5

        tick_int = _nice_tick(total)
        bp = 0
        while bp < total:
            angle  = bp2a(bp)
            tx, ty = a2xy(angle, dr=TICK_DR_MARK)
            canvas.put(tx, ty, "┼", "color(250)")
            label  = _format_bp(bp)
            lx, ly = a2xy(angle, dr=TICK_DR_LABEL)
            # Inside the circle: text points inward, so alignment is flipped
            # vs outside placement (right-side labels hang left, left-side hang right).
            if math.cos(angle) >= 0:
                canvas.put_text(lx - len(label) + 1, ly, label, "color(245)")
            else:
                canvas.put_text(lx, ly, label, "color(245)")
            bp += tick_int

        # ── Restriction site marks ─────────────────────────────────────────────
        # resite: draw a thin arc just outside (+) or inside (-) the backbone
        # recut:  draw a radial tick crossing the backbone at the cut position
        restr_labels: list[tuple[float, str, str]] = []
        for rf in self._restr_feats:
            color = rf["color"]
            if rf["type"] == "resite":
                dr_lo = 4 if rf["strand"] >= 0 else -5
                dr_hi = dr_lo + 1
                start_a = bp2a(rf["start"])
                end_a   = bp2a(rf["end"])
                span_a  = (end_a - start_a) % TWO_PI or TWO_PI
                arc_steps = max(4, int(span_a / TWO_PI * n_steps))
                inv_arc   = 1.0 / arc_steps if arc_steps else 1.0
                for dr_i in range(2):
                    dr     = dr_lo + (dr_hi - dr_lo) * dr_i
                    _rx2dr = (rx + dr) * 2
                    _ry4dr = (ry + dr * _ry_rx) * 4
                    for s in range(arc_steps + 1):
                        a  = start_a + s * inv_arc * span_a
                        px = round(cx2 + _rx2dr * math.cos(a))
                        py = round(cy4 + _ry4dr * math.sin(a))
                        pcol, prow = px >> 1, py >> 2
                        if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                            bc_bits[prow][pcol]   |= 1 << _BIT[(py & 3) << 1 | (px & 1)]
                            if 1 >= bc_prio[prow][pcol]:
                                bc_colors[prow][pcol] = color
                                bc_prio[prow][pcol]   = 1
                # Collect label at midpoint
                mid_bp = (rf["start"] + rf["end"]) // 2
                restr_labels.append((bp2a(mid_bp), rf["label"], color))
            elif rf["type"] == "recut":
                # Radial tick crossing the backbone at cut position
                cut_a = bp2a(rf["start"])
                for dr in range(-1, 3):
                    tx, ty = a2xy(cut_a, dr=dr)
                    if 0 <= tx < w and 0 <= ty < h:
                        canvas.put(tx, ty, "┼" if dr == 0 else "·", color)

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

            # Arc length on a circular plasmid: ((end - start) mod total) handles
            # wrap-around (end < start) without putting the label opposite the arc.
            arc_len = (end_bp - start_bp) % total
            mid_bp  = (start_bp + arc_len // 2) % total
            label_slots.append((bp2a(mid_bp), f["label"], color))

        # Add restriction site labels (from resite entries only)
        for angle, lbl, color in restr_labels:
            label_slots.append((angle, lbl, color))

        # ── Labels: place each as close to the arc as possible ───────────────
        # Greedily try increasing dr until the label's bounding box doesn't
        # overlap any already-placed label.
        # ↓ Tune this to control how far labels sit from the arc.
        LABEL_DR_MIN = 9          # minimum radial clearance from arc edge
        dr_min = LABEL_DR_MIN
        dr_max = max(rx // 2 + 6, LABEL_DR_MIN + 10)

        # placed: dict keyed by row → list of (x0, x1) bounding boxes
        placed_by_row: dict[int, list] = {}
        final_labels: list[tuple[float, str, str, int, int, int]] = []
        # angle, lbl, color, chosen_dr, lx, ly

        for angle, lbl, color in label_slots:
            on_right = math.cos(angle) >= 0
            chosen = None
            for dr in range(dr_min, dr_max + 1):
                lx, ly = a2xy(angle, dr=dr)
                if not (0 <= ly < h):
                    continue
                lbl_x0 = lx if on_right else max(0, lx - len(lbl) + 1)
                lbl_x1 = lx + len(lbl) - 1 if on_right else lx
                # Check against already-placed boxes on the same row only
                ok = True
                for bx0, bx1 in placed_by_row.get(ly, []):
                    if not (lbl_x1 < bx0 or lbl_x0 > bx1):
                        ok = False
                        break
                if ok:
                    chosen = (dr, lx, ly, lbl_x0, lbl_x1)
                    break
            if chosen is None:
                # Couldn't fit without overlap — place at max dr anyway
                lx, ly = a2xy(angle, dr=dr_max)
                lbl_x0 = lx if on_right else max(0, lx - len(lbl) + 1)
                lbl_x1 = lx + len(lbl) - 1 if on_right else lx
                chosen = (dr_max, lx, ly, lbl_x0, lbl_x1)
            dr_c, lx, ly, lbl_x0, lbl_x1 = chosen
            placed_by_row.setdefault(ly, []).append((lbl_x0, lbl_x1))
            final_labels.append((angle, lbl, color, dr_c, lx, ly))

        # Render
        for angle, lbl, color, dr_c, lx, ly in final_labels:
            on_right = math.cos(angle) >= 0

            # Dot just outside the arc
            dot_x, dot_y = a2xy(angle, dr=3)
            canvas.put(dot_x, dot_y, "·", color)

            # Optional connector line from dot to label
            if self._show_connectors:
                lbl_mid_x = lx + len(lbl) // 2 if on_right else lx - len(lbl) // 2
                steps = max(1, abs(lbl_mid_x - dot_x) + abs(ly - dot_y))
                for t in range(1, steps):
                    px = dot_x + (lbl_mid_x - dot_x) * t // steps
                    py = dot_y + (ly        - dot_y) * t // steps
                    if 0 <= px < w and 0 <= py < h and canvas._chars[py][px] == " ":
                        canvas.put(px, py, "·", color)

            if on_right:
                canvas.put_text(lx, ly, lbl, color)
            else:
                canvas.put_text(max(0, lx - len(lbl) + 1), ly, lbl, color)

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

    # ── Linear map ─────────────────────────────────────────────────────────────

    def _feat_at_linear(self, x: int, y: int) -> int:
        """Return feature index at terminal cell (x, y) in linear view, or -1."""
        if not self._total:
            return -1
        w, h      = self.size.width, self.size.height
        margin_l  = 5
        margin_r  = 2
        usable_w  = w - margin_l - margin_r
        backbone_row = max(4, h // 2)
        if x < margin_l or x >= w - margin_r or usable_w <= 0:
            return -1
        bp = int((x - margin_l) / usable_w * self._total)
        above = y < backbone_row
        below = y > backbone_row
        for i, f in enumerate(self._feats):
            s, e = f["start"], f["end"]
            in_range = (s <= bp <= e) if e > s else (bp >= s or bp <= e)
            if not in_range:
                continue
            if above and f["strand"] >= 0:
                return i
            if below and f["strand"] < 0:
                return i
        return -1

    def _draw_linear(self, w: int, h: int) -> Text:
        """Render a horizontal linear plasmid map."""
        canvas = _Canvas(w, h)
        bc     = _BrailleCanvas(w, h)
        total  = self._total

        if not total:
            canvas.put_text(w // 2 - 9, h // 2, "No record loaded", "dim")
            return bc.combine(canvas)

        # ── Layout ──
        margin_l     = 5
        margin_r     = 2
        usable_w     = w - margin_l - margin_r
        px_w         = usable_w * 2
        px_start     = margin_l * 2
        backbone_row = max(4, h // 2)
        backbone_py  = backbone_row * 4 + 1   # braille pixel row

        _BIT      = (0, 3, 1, 4, 2, 5, 6, 7)
        bc_bits   = bc._bits
        bc_colors = bc._colors
        bc_prio   = bc._prio
        bc_cols   = bc.cols
        bc_rows   = bc.rows

        def bp_to_px(bp: int) -> int:
            return px_start + int(bp / total * px_w)

        def _set(px: int, py: int, color: str, prio: int) -> None:
            pcol, prow = px >> 1, py >> 2
            if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                bc_bits[prow][pcol]   |= 1 << _BIT[(py & 3) << 1 | (px & 1)]
                if prio >= bc_prio[prow][pcol]:
                    bc_colors[prow][pcol] = color
                    bc_prio[prow][pcol]   = prio

        # ── Backbone ──
        for px in range(px_start, px_start + px_w + 1):
            _set(px, backbone_py, "color(238)", 1)

        # ── Ticks + bp labels ──
        tick_int  = _nice_tick(total)
        label_row = backbone_row + 1
        bp = 0
        while bp <= total:
            tx = margin_l + bp * usable_w // total
            canvas.put(tx, backbone_row, "┼", "color(250)")
            lbl = _format_bp(bp)
            if label_row < h:
                canvas.put_text(tx - len(lbl) // 2, label_row, lbl, "color(245)")
            bp += tick_int
        # Right cap
        end_tx = min(margin_l + usable_w, w - 1)
        canvas.put(end_tx, backbone_row, "┤", "color(250)")

        # ── Restriction site marks ──
        # resite: thin braille arc above (fwd) or below (rev) backbone
        # recut:  radial tick at cut position
        for rf in self._restr_feats:
            color = rf["color"]
            if rf["type"] == "resite":
                x0 = margin_l + rf["start"] * usable_w // total
                x1 = margin_l + rf["end"]   * usable_w // total
                x0 = max(margin_l, min(x0, w - margin_r - 1))
                x1 = max(margin_l, min(x1, w - margin_r - 1))
                dy = -2 if rf["strand"] >= 0 else 2
                bar_py = backbone_py + dy * 4
                for bx in range(x0 * 2, x1 * 2 + 1):
                    pcol, prow = bx >> 1, bar_py >> 2
                    if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                        bc_bits[prow][pcol]   |= 1 << _BIT[(bar_py & 3) << 1 | (bx & 1)]
                        if 1 >= bc_prio[prow][pcol]:
                            bc_colors[prow][pcol] = color
                            bc_prio[prow][pcol]   = 1
            elif rf["type"] == "recut":
                cut_x = margin_l + rf["start"] * usable_w // total
                if margin_l <= cut_x < w - margin_r:
                    row_above = backbone_row - 1
                    row_below = backbone_row + 1
                    canvas.put(cut_x, backbone_row, "┼", color)
                    if 0 <= row_above < h:
                        canvas.put(cut_x, row_above, "↓" if rf["strand"] >= 0 else " ", color)
                    if 0 <= row_below < h:
                        canvas.put(cut_x, row_below, "↑" if rf["strand"] < 0 else " ", color)

        # ── Lane assignment (greedy interval scheduling) ──
        fwd_ends: list[int] = []   # rightmost braille-x used per fwd lane
        rev_ends: list[int] = []
        feat_meta: list[tuple[bool, int]] = []   # (is_fwd, lane_idx)

        for feat in self._feats:
            is_fwd = feat["strand"] >= 0
            x0, x1 = bp_to_px(feat["start"]), bp_to_px(feat["end"])
            ends   = fwd_ends if is_fwd else rev_ends
            lane   = len(ends)
            for li, ex in enumerate(ends):
                if x0 > ex + 2:
                    lane = li
                    ends[li] = x1
                    break
            else:
                ends.append(x1)
            feat_meta.append((is_fwd, lane))

        # ── Draw features ──
        lane_gap    = 3    # braille rows between backbone and first lane
        bar_half    = 2    # half-height of feature bar in braille rows
        lane_stride = 7    # braille rows between lane centres

        for i, (feat, (is_fwd, lane)) in enumerate(zip(self._feats, feat_meta)):
            start_bp = feat["start"]
            end_bp   = feat["end"]
            strand   = feat["strand"]
            color    = feat["color"]
            label    = feat.get("label", feat.get("type", ""))

            # Handle wrap-around features
            if end_bp > start_bp:
                segments = [(bp_to_px(start_bp), bp_to_px(end_bp))]
            else:
                segments = [(bp_to_px(start_bp), px_start + px_w),
                            (px_start, bp_to_px(end_bp))]

            center_py = (backbone_py - lane_gap - bar_half - lane * lane_stride
                         if is_fwd else
                         backbone_py + lane_gap + bar_half + lane * lane_stride)
            bar_top    = center_py - bar_half
            bar_bottom = center_py + bar_half

            if bar_top < 0 or bar_bottom >= h * 4:
                continue

            is_sel = (i == self.selected_idx)
            prio   = 3 if is_sel else 2
            style  = ("reverse " + color) if is_sel else color
            feat_ty = center_py >> 2

            # Label row: above bar for forward, below for reverse
            label_ty = feat_ty - 1 if is_fwd else feat_ty + 1

            for sx0, sx1 in segments:
                sx0 = max(px_start, min(sx0, px_start + px_w))
                sx1 = max(px_start, min(sx1, px_start + px_w))
                if sx1 <= sx0:
                    continue

                # Fill bar in braille
                for py in range(bar_top, bar_bottom + 1):
                    for px in range(sx0, sx1):
                        _set(px, py, style, prio)

                # Arrowhead character
                if strand >= 0:
                    canvas.put(min(sx1 >> 1, w - margin_r - 1), feat_ty, "▶", style)
                else:
                    canvas.put(max(sx0 >> 1, margin_l), feat_ty, "◀", style)

                # Label above (fwd) or below (rev) the bar
                x0c, x1c = sx0 >> 1, sx1 >> 1
                max_lbl_w = x1c - x0c
                if max_lbl_w >= 1 and 0 <= label_ty < h:
                    lbl = label[:max_lbl_w]
                    lx  = x0c + (max_lbl_w - len(lbl)) // 2
                    canvas.put_text(lx, label_ty, lbl, style)

        # ── Header ──
        name = (self.record.name or self.record.id or "?")[:w // 3]
        canvas.put_text(margin_l, 0, f"{name}  {total:,} bp", "bold white")
        hint = "[ linear  ·  v = circular ]"
        canvas.put_text(w - len(hint) - 1, 0, hint, "dim")

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
        # Suppress the RowHighlighted cascade that fires when DataTable auto-
        # moves the cursor to row 0 after clear()+add_row. Without this guard,
        # every record load triggers a redundant SequencePanel rebuild (the
        # cascade routes through PlasmidApp._sidebar_row_activated → seq_pnl.
        # highlight_feature → _refresh_view). Measured savings: ~50 ms per
        # load on 10 kb plasmids, ~110 ms on 20 kb. See CLAUDE.md perf notes.
        #
        # Textual posts RowHighlighted asynchronously, so we can't use a
        # simple try/finally: the handler would see _populating == False
        # because populate() has already returned by the time the message is
        # dispatched. Instead we set the flag True here and schedule its
        # reset via call_after_refresh, which runs AFTER all pending messages.
        self._populating = True
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
        def _clear_populating():
            self._populating = False
        self.call_after_refresh(_clear_populating)

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
    _populating: bool = False   # suppress RowActivated cascade during populate()

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
        # Ignore every auto-highlight that fires while populate() is running.
        if self._populating:
            return
        if event.cursor_row == self._prog_row:
            self._prog_row = -1
            return
        self.post_message(self.RowActivated(event.cursor_row))


# ── Library panel ──────────────────────────────────────────────────────────────

class LibraryPanel(Widget):
    """Left-hand plasmid library — persistent SnapGene-style collection."""

    DEFAULT_CSS = """
    LibraryPanel {
        width: 26;
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

    class AnnotateRequested(Message):
        """User pressed the annotate button on a selected library row.
        entry_id is the library key of the row with the cursor."""
        def __init__(self, entry_id: "str | None"):
            self.entry_id = entry_id
            super().__init__()

    class RenameRequested(Message):
        """User pressed the rename (✎) button with a library row focused.
        entry_id is the library key of the row with the cursor."""
        def __init__(self, entry_id: "str | None"):
            self.entry_id = entry_id
            super().__init__()

    class GainedFocus(Message):
        """Posted when the library panel (or any descendant) gains focus.
        The app handles this by clearing the currently-selected feature so
        pressing Delete in the library cannot accidentally hit a feature
        that's no longer visually in focus."""
        pass

    def on_descendant_focus(self, _event):
        # DataTable inside the library panel is what actually gets focus;
        # Textual reports it to us via a DescendantFocus event. Propagate
        # to the app as a panel-level signal.
        self.post_message(self.GainedFocus())

    def compose(self) -> ComposeResult:
        yield Static(" Library", id="lib-hdr")
        yield DataTable(id="lib-table", cursor_type="row", zebra_stripes=True)
        with Horizontal(id="lib-btns"):
            yield Button("+", id="btn-lib-add", variant="success",
                         tooltip="Add current plasmid")
            yield Button("−", id="btn-lib-del", variant="error",
                         tooltip="Remove selected")
            yield Button("◈", id="btn-lib-annot", variant="primary",
                         tooltip="Annotate selected (pLannotate)  —  shortcut: Shift+A")
            yield Button("✎", id="btn-lib-rename", variant="default",
                         tooltip="Rename selected plasmid")

    def on_mount(self):
        self._active_id:    "str | None" = None
        self._active_dirty: bool         = False
        t = self.query_one("#lib-table", DataTable)
        t.add_columns("Name", "bp")
        self._repopulate()

    def _repopulate(self):
        t = self.query_one("#lib-table", DataTable)
        t.clear()
        for entry in _load_library():
            is_dirty = (entry["id"] == self._active_id and self._active_dirty)
            name_disp = ("*" + entry["name"])[:14] if is_dirty else entry["name"][:14]
            t.add_row(
                name_disp,
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

    @on(Button.Pressed, "#btn-lib-annot")
    def _btn_annotate(self):
        # Annotate the currently-focused row if any, else the active record.
        t = self.query_one("#lib-table", DataTable)
        entry_id: "str | None" = None
        if t.row_count > 0 and 0 <= t.cursor_row < t.row_count:
            row_keys = list(t.rows.keys())
            if 0 <= t.cursor_row < len(row_keys):
                entry_id = row_keys[t.cursor_row].value
        self.post_message(self.AnnotateRequested(entry_id))

    @on(Button.Pressed, "#btn-lib-rename")
    def _btn_rename(self):
        # Rename only works on the row with the DataTable cursor — if the
        # library is empty or no row is focused, we send None and the app
        # will notify the user.
        t = self.query_one("#lib-table", DataTable)
        entry_id: "str | None" = None
        if t.row_count > 0 and 0 <= t.cursor_row < t.row_count:
            row_keys = list(t.rows.keys())
            if 0 <= t.cursor_row < len(row_keys):
                entry_id = row_keys[t.cursor_row].value
        self.post_message(self.RenameRequested(entry_id))

    def set_active(self, entry_id: "str | None") -> None:
        """Mark which library entry is currently loaded (clears dirty flag)."""
        self._active_id    = entry_id
        self._active_dirty = False

    def set_dirty(self, dirty: bool) -> None:
        """Show unsaved-changes marker on the active row and in the panel header."""
        self._active_dirty = dirty
        self._repopulate()
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
        self._show_connectors:  bool = False
        self._re_highlight: "dict | None" = None  # RE cut visualization
        self._sel_anchor:   int         = -1    # anchor for Shift+arrow extension
        # Drag-to-select state
        self._drag_start_bp:    int  = -1
        self._has_dragged:      bool = False
        self._mouse_button_held: bool = False
        self._drag_was_shift:   bool = False
        self._last_was_drag:    bool = False
        # Set by _click_to_bp when the click lands on a resite bar row
        self._last_resite_click: "dict | None" = None
        self._sorted_feats_cache: "list | None" = None

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
        self._seq          = seq
        self._feats        = feats
        self._sorted_feats_cache = None
        self._sel_range    = None
        self._user_sel     = None
        self._cursor_pos   = -1
        self._re_highlight = None
        self._sel_anchor   = -1
        self.remove_class("has-trans")
        self._refresh_view()

    def highlight_feature(self, feat: "dict | None", cursor_bp: int = -1) -> None:
        """Highlight a feature's region in the sequence; show CDS translation.

        cursor_bp: if >= 0, anchor the cursor (and scroll) at this bp position.
                   Use this for sequence-panel clicks so scroll stays at the
                   clicked position rather than jumping to the feature start.
        """
        self._re_highlight = None
        if feat is None or not self._seq:
            self._sel_range = None
            self.remove_class("has-trans")
            self._refresh_view()
            return

        start, end = feat["start"], min(feat["end"], len(self._seq))
        self._sel_range  = (start, end)
        self._user_sel   = None          # clear shift-selection on programmatic highlight
        self._sel_anchor = -1
        if cursor_bp >= 0:
            self._cursor_pos = cursor_bp
        self._refresh_view()
        if cursor_bp >= 0:
            self._ensure_cursor_visible()

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

    def select_feature_range(self, feat: dict, cursor_bp: int = -1) -> None:
        """Highlight the entire feature span as a copyable selection.

        cursor_bp: if >= 0 (sequence-panel click), keep cursor at that bp and
                   anchor the Shift+arrow selection there. Otherwise cursor goes
                   to the feature end so Shift+arrow naturally extends outward.
        """
        if not self._seq or feat is None:
            return
        start = feat["start"]
        end   = min(feat["end"], len(self._seq))
        self._user_sel  = (start, end)
        self._sel_range = None
        self._re_highlight = None
        if cursor_bp >= 0:
            # Sequence-panel click: cursor stays at click; anchor at feature start
            self._cursor_pos = cursor_bp
            self._sel_anchor = start
        else:
            # Map/sidebar/double-click: cursor at feature end, anchor at start
            self._cursor_pos = max(start, end - 1)
            self._sel_anchor = start
        self._refresh_view()
        self._ensure_cursor_visible()

        # Show CDS translation when applicable
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
            # Shift+click: extend selection from anchor (or cursor) to here
            anchor = self._sel_anchor if self._sel_anchor >= 0 else self._cursor_pos
            self._sel_anchor = anchor
            s = min(anchor, bp)
            e = max(anchor, bp) + 1
            self._user_sel   = (s, e)
            self._cursor_pos = bp
            self._sel_range  = None
        else:
            # Plain click: place cursor, clear selection and anchor
            self._cursor_pos = bp
            self._user_sel   = None
            self._sel_anchor = -1
        self._refresh_view()
        self._ensure_cursor_visible()

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
        self._last_resite_click = None
        bp = self._click_to_bp(event.screen_x, event.screen_y)
        if bp < 0:
            return

        # If the click landed on a restriction site bar, highlight the recognition span
        resite = self._last_resite_click
        self._last_resite_click = None
        if resite is not None:
            hi_start = resite["start"]
            hi_end   = min(resite["end"], len(self._seq))
            ext_cut  = resite.get("ext_cut_bp")
            if ext_cut is not None:
                if ext_cut >= hi_end:
                    hi_end   = min(ext_cut + 1, len(self._seq))
                elif ext_cut < hi_start:
                    hi_start = ext_cut
            self._re_highlight = {
                "start":      hi_start,
                "end":        hi_end,
                "fwd_cut_bp": -1,
                "rev_cut_bp": -1,
                "color":      resite["color"],
                "name":       resite["label"],
            }
            self._sel_range  = None
            self._user_sel   = None
            self._cursor_pos = -1
            self._refresh_view()
            return

        # Regular click: clear any RE highlight
        self._re_highlight = None
        double = event.chain >= 2
        self.post_message(self.SequenceClick(bp, double=double))

    def _seq_render_width(self) -> int:
        """Character width of the render area, minus the 2-col vertical scrollbar."""
        return max(20, self.size.width - 2)

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

        n           = len(self._seq)
        num_w       = len(str(n)) if n else 1
        line_width  = max(20, self._seq_render_width() - (num_w + 2))
        annot_feats = sorted(
            [f for f in self._feats if f.get("type") not in ("site", "recut")],
            key=lambda f: -(f["end"] - f["start"]),
        )
        rpg     = 2 + (1 if self._show_connectors else 0)  # rows per feature group
        row     = 0
        seq_col = vp_x - (num_w + 2)   # offset past the num+2-space prefix

        def _check_lane(lane):
            """Check if click hit a feature in this lane; return bp or None."""
            for f in lane:
                bar_s = max(f["start"], chunk_start) - chunk_start
                bar_e = min(f["end"],   chunk_end)   - chunk_start
                if bar_s <= seq_col < bar_e:
                    if f.get("type") == "resite":
                        self._last_resite_click = f
                    return (f["start"] + f["end"]) // 2
            return -1

        for chunk_start in range(0, n, line_width):
            chunk_end   = min(chunk_start + line_width, n)
            chunk_feats = [f for f in annot_feats
                           if f["start"] < chunk_end and f["end"] > chunk_start]
            re_above, onebp_above, reg_above, reg_below, onebp_below, re_below = (
                _chunk_lane_groups(chunk_feats, chunk_start, chunk_end)
            )

            # Above: RE (far) → 1bp → multi-bp (close to DNA)
            for lane in (*re_above, *onebp_above, *reg_above):
                for _ in range(rpg):
                    if row == content_row:
                        return _check_lane(lane)
                    row += 1

            # DNA rows: fwd strand + RC strand (2 rows)
            for _ in range(2):
                if row == content_row:
                    if 0 <= seq_col < (chunk_end - chunk_start):
                        return chunk_start + seq_col
                    return -1
                row += 1

            # Below: multi-bp (close) → 1bp → RE (far)
            for lane in (*reg_below, *onebp_below, *re_below):
                for _ in range(rpg):
                    if row == content_row:
                        return _check_lane(lane)
                    row += 1

            if row > content_row:
                break
        return -1

    # ── Helpers ────────────────────────────────────────────────────────────────

    def _annot_feats_sorted(self) -> list:
        if self._sorted_feats_cache is None:
            self._sorted_feats_cache = sorted(
                [f for f in self._feats if f.get("type") not in ("site", "recut")],
                key=lambda f: -(f["end"] - f["start"]),
            )
        return self._sorted_feats_cache

    def _bp_to_content_row(self, bp: int) -> int:
        """Return the content row index (0-based) of the DNA line containing bp."""
        n           = len(self._seq)
        num_w       = len(str(n)) if n else 1
        line_width  = max(20, self._seq_render_width() - (num_w + 2))
        annot_feats = self._annot_feats_sorted()
        rpg = 2 + (1 if self._show_connectors else 0)
        row = 0
        for chunk_start in range(0, n, line_width):
            chunk_end   = min(chunk_start + line_width, n)
            chunk_feats = [f for f in annot_feats
                           if f["start"] < chunk_end and f["end"] > chunk_start]
            re_above, onebp_above, reg_above, reg_below, onebp_below, re_below = (
                _chunk_lane_groups(chunk_feats, chunk_start, chunk_end)
            )
            above_rows = (len(re_above) + len(onebp_above) + len(reg_above)) * rpg
            if bp < chunk_end:
                return row + above_rows   # forward-strand DNA row within this chunk
            below_rows = (len(reg_below) + len(onebp_below) + len(re_below)) * rpg
            row += above_rows + 2 + below_rows
        return row

    def _line_width(self) -> int:
        """Number of bp per displayed line."""
        n = len(self._seq)
        num_w = len(str(n)) if n else 1
        return max(20, self._seq_render_width() - (num_w + 2))

    def _scroll_to_row(self, row: int) -> None:
        try:
            self.query_one("#seq-scroll", ScrollableContainer).scroll_to(
                0, row, animate=False
            )
        except Exception:
            pass

    def _ensure_cursor_visible(self) -> None:
        """Scroll just enough so the cursor's DNA line is fully visible.

        Deferred to run after the next refresh so that the scroll container
        has already processed any content update (which can reset scroll_y).
        """
        if self._cursor_pos < 0 or not self._seq:
            return
        row = self._bp_to_content_row(self._cursor_pos)
        row_bottom = row + 1   # fwd + rc strand = 2 rows

        def _do_scroll():
            try:
                scroll = self.query_one("#seq-scroll", ScrollableContainer)
                vp_top = int(scroll.scroll_y)
                vp_h   = scroll.size.height
                vp_bottom = vp_top + vp_h - 1
                if row < vp_top:
                    scroll.scroll_to(0, row, animate=False)
                elif row_bottom > vp_bottom:
                    scroll.scroll_to(0, row_bottom - vp_h + 1, animate=False)
            except Exception:
                pass

        self.call_after_refresh(_do_scroll)

    def _refresh_view(self) -> None:
        view = self.query_one("#seq-view", Static)
        try:
            scroll = self.query_one("#seq-scroll", ScrollableContainer)
        except Exception:
            scroll = None
        if not self._seq:
            view.update(Text("  No sequence loaded.", style="dim italic"))
            return
        # num_w-char line number + "  " (2) + seq = num_w + 2 overhead
        # Use actual scroll-container content width (excludes 2-col vertical scrollbar)
        num_w      = len(str(len(self._seq))) if self._seq else 1
        line_width = max(20, self._seq_render_width() - (num_w + 2))
        reh_key = (
            self._re_highlight["start"], self._re_highlight["end"]
        ) if self._re_highlight else None
        key = (id(self._seq), id(self._feats), line_width,
               self._sel_range, self._user_sel, self._cursor_pos,
               self._show_connectors, reh_key)
        if key != self._view_cache_key:
            self._view_cache_txt = _build_seq_text(
                self._seq, self._feats,
                line_width      = line_width,
                sel_range       = self._sel_range,
                user_sel        = self._user_sel,
                cursor_pos      = self._cursor_pos,
                show_connectors = self._show_connectors,
                re_highlight    = self._re_highlight,
            )
            self._view_cache_key = key

        # Preserve scroll position across content update
        saved_y = scroll.scroll_y if scroll is not None else None
        view.update(self._view_cache_txt)
        if saved_y is not None and saved_y > 0:
            def _restore():
                if scroll is not None:
                    scroll.scroll_to(0, saved_y, animate=False)
            self.call_after_refresh(_restore)

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
            _log.exception("NCBI fetch failed for %s", acc)
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


# ── Dropdown menu modal ────────────────────────────────────────────────────────

class DropdownScreen(ModalScreen):
    """Lightweight overlay showing a positioned dropdown menu.

    Uses a near-transparent backdrop so the main app stays visible — the
    dropdown looks like a real popup anchored to the menu bar, not a
    separate "screen". Click outside the box dismisses.
    """

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    DEFAULT_CSS = """
    DropdownScreen {
        background: rgba(0, 0, 0, 0.15);
    }
    """

    def __init__(self, items: list, x: int, y: int) -> None:
        super().__init__()
        self._items = items   # (label, action_str | None)
        self._x = x
        self._y = y
        self._cursor = next(
            (i for i, (_, a) in enumerate(items) if a is not None), 0
        )

    def compose(self) -> ComposeResult:
        yield Static(
            self._render_content(),
            id="dropdown-box",
        )

    def on_mount(self) -> None:
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        box_h   = len(self._items) + 2
        box = self.query_one("#dropdown-box", Static)
        box.styles.offset = (self._x, self._y)
        box.styles.width  = inner_w
        box.styles.height = box_h
        box.styles.border = ("solid", "#555555")
        box.styles.background = "#1e1e1e"

    def _render_content(self) -> Text:
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        sep_line = "\u2500" * (inner_w - 2)
        result = Text()
        for i, (label, action) in enumerate(self._items):
            is_sep      = (label == "---")
            is_selected = (i == self._cursor and not is_sep and action is not None)
            is_disabled = (action is None and not is_sep)

            if is_sep:
                line = Text(sep_line + "\n", style="white")
            else:
                padded = f" {label:<{inner_w - 3}}"
                if is_selected:
                    line = Text(padded + "\n", style="reverse white")
                elif is_disabled:
                    line = Text(padded + "\n", style="dim white")
                else:
                    line = Text(padded + "\n", style="white")
            result.append_text(line)
        return result

    def _refresh_box(self) -> None:
        box = self.query_one("#dropdown-box", Static)
        box.update(self._render_content())

    def on_key(self, event) -> None:
        items = self._items
        if event.key == "up":
            pos = self._cursor - 1
            while pos >= 0 and (items[pos][0] == "---" or items[pos][1] is None):
                pos -= 1
            if pos >= 0:
                self._cursor = pos
                self._refresh_box()
            event.stop()
        elif event.key == "down":
            pos = self._cursor + 1
            while pos < len(items) and (items[pos][0] == "---" or items[pos][1] is None):
                pos += 1
            if pos < len(items):
                self._cursor = pos
                self._refresh_box()
            event.stop()
        elif event.key == "enter":
            label, action = items[self._cursor]
            if action is not None:
                self.dismiss(action)
            event.stop()

    def on_click(self, event: Click) -> None:
        bx = self._x
        by = self._y
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        bh = len(self._items) + 2
        cx, cy = event.screen_x, event.screen_y
        if bx <= cx < bx + inner_w and by <= cy < by + bh:
            row_in_box = cy - by - 1  # -1 for top border
            if 0 <= row_in_box < len(self._items):
                label, action = self._items[row_in_box]
                if label == "---" or action is None:
                    event.stop()
                    return
                self.dismiss(action)
        else:
            self.dismiss(None)
        event.stop()

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Menu bar widget ────────────────────────────────────────────────────────────

class MenuBar(Widget):
    DEFAULT_CSS = """
    MenuBar {
        height: 1;
        background: $primary-darken-3;
        layout: horizontal;
    }
    MenuBar Static {
        width: auto;
        padding: 0 2;
        color: $text;
    }
    MenuBar Static:hover {
        background: $primary;
    }
    """

    MENUS = ["File", "Edit", "Enzymes", "Features", "Primers", "Parts", "Constructor"]

    def compose(self) -> ComposeResult:
        for name in self.MENUS:
            yield Static(name, classes="menu-item", id=f"menu-{name.lower()}")

    def on_click(self, event: Click) -> None:
        for name in self.MENUS:
            widget_id = f"menu-{name.lower()}"
            try:
                item = self.query_one(f"#{widget_id}", Static)
                region = item.region
                if (region.x <= event.screen_x < region.x + region.width and
                        region.y <= event.screen_y < region.y + region.height):
                    x = region.x
                    y = region.y + 1
                    self.app.open_menu(name, x, y)
                    break
            except Exception:
                pass


# ── Golden Braid L0 part catalog (shared by Parts Bin and Constructor) ─────────

# (Name, Type, Position, 5' OH, 3' OH, Backbone, Selection Marker)
_GB_L0_PARTS: list[tuple] = [
    # ── Promoters (Position 1: GGAG → TGAC) ────────────────────────────
    ("CaMV 35S",          "Promoter",   "Pos 1",   "GGAG", "TGAC", "pUPD2", "Spectinomycin"),
    ("Nos",               "Promoter",   "Pos 1",   "GGAG", "TGAC", "pUPD2", "Spectinomycin"),
    ("AtUBQ10",           "Promoter",   "Pos 1",   "GGAG", "TGAC", "pUPD2", "Spectinomycin"),
    ("ZmUBI1",            "Promoter",   "Pos 1",   "GGAG", "TGAC", "pUPD2", "Spectinomycin"),
    ("AtRPS5a",           "Promoter",   "Pos 1",   "GGAG", "TGAC", "pUPD2", "Spectinomycin"),
    # ── 5′ UTRs (Position 2: TGAC → AATG) ─────────────────────────────
    ("TMV Omega",         "5' UTR",     "Pos 2",   "TGAC", "AATG", "pUPD2", "Spectinomycin"),
    ("Nos 5'UTR",         "5' UTR",     "Pos 2",   "TGAC", "AATG", "pUPD2", "Spectinomycin"),
    ("ADH1 5'UTR",        "5' UTR",     "Pos 2",   "TGAC", "AATG", "pUPD2", "Spectinomycin"),
    # ── CDS with stop (Positions 3-4: AATG → GCTT) ─────────────────────
    ("eGFP",              "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("mCherry",           "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("mVenus",            "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("mTurquoise2",       "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("GUS (uidA)",        "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("Luciferase (LUC+)", "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("NptII (KanR)",      "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("hptII (HygR)",      "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("Bar (BastaR)",      "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("Cas9 (SpCas9)",     "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    # ── CDS without stop (Position 3: AATG → TTCG) ─────────────────────
    ("eGFP (no stop)",    "CDS-NS",     "Pos 3",   "AATG", "TTCG", "pUPD2", "Spectinomycin"),
    ("mCherry (no stop)", "CDS-NS",     "Pos 3",   "AATG", "TTCG", "pUPD2", "Spectinomycin"),
    # ── C-terminal tags (Position 4: TTCG → GCTT) ──────────────────────
    ("GFP C-tag",         "C-tag",      "Pos 4",   "TTCG", "GCTT", "pUPD2", "Spectinomycin"),
    ("HA tag",            "C-tag",      "Pos 4",   "TTCG", "GCTT", "pUPD2", "Spectinomycin"),
    ("6xHis tag",         "C-tag",      "Pos 4",   "TTCG", "GCTT", "pUPD2", "Spectinomycin"),
    # ── Terminators (Position 5: GCTT → CGCT) ──────────────────────────
    ("Nos terminator",    "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
    ("CaMV 35S term",     "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
    ("OCS terminator",    "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
    ("rbcS terminator",   "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
    ("HSP18.2 term",      "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
]

_GB_TYPE_COLORS: dict[str, str] = {
    "Promoter":   "green",
    "5' UTR":     "cyan",
    "CDS":        "yellow",
    "CDS-NS":     "dark_orange",
    "C-tag":      "magenta",
    "Terminator": "blue",
}

# Canonical Golden Braid L0 part positions.
# Each maps a part-type name → (position label, 5' overhang, 3' overhang).
# Overhangs follow the published GB2.0 standard (Sarrion-Perdigones et al. 2013).
_GB_POSITIONS: dict[str, tuple[str, str, str]] = {
    "Promoter":    ("Pos 1",   "GGAG", "TGAC"),
    "5' UTR":      ("Pos 2",   "TGAC", "AATG"),
    "CDS":         ("Pos 3-4", "AATG", "GCTT"),
    "CDS-NS":      ("Pos 3",   "AATG", "TTCG"),
    "C-tag":       ("Pos 4",   "TTCG", "GCTT"),
    "Terminator":  ("Pos 5",   "GCTT", "CGCT"),
}

# BsaI recognition + tail used for all Golden Braid domestication primers.
# Padding bases improve BsaI digestion efficiency near DNA ends.
_GB_BSAI_SITE = "GGTCTC"
_GB_SPACER    = "A"           # 1 nt between recognition and the overhang
_GB_PAD       = "GCGC"        # 4 nt of extra bases for efficient end-cutting


def _pick_binding_region(seq: str, target_tm: float = 57.0,
                         min_len: int = 18, max_len: int = 25) -> tuple[str, float]:
    """Return the prefix of `seq` (length min_len..max_len) whose Tm is
    closest to `target_tm`. Uses primer3-py's SantaLucia Tm calculation.

    Returns (binding_sequence, tm). If primer3-py is not installed, falls
    back to a crude 2+4 rule estimate.
    """
    try:
        import primer3
        _tm = primer3.calc_tm
    except ImportError:
        # Crude fallback so the code still runs without primer3-py; the UI
        # will warn the user that Tm values are approximate.
        def _tm(s):
            gc = sum(1 for c in s.upper() if c in "GC")
            at = sum(1 for c in s.upper() if c in "AT")
            return 2 * at + 4 * gc

    best_seq, best_tm, best_diff = seq[:min_len], 0.0, float("inf")
    for n in range(min_len, min(max_len + 1, len(seq) + 1)):
        candidate = seq[:n]
        tm = _tm(candidate)
        diff = abs(tm - target_tm)
        if diff < best_diff:
            best_seq, best_tm, best_diff = candidate, tm, diff
    return best_seq, best_tm


def _design_gb_primers(
    template_seq: str,
    start: int,
    end: int,
    part_type: str,
    target_tm: float = 57.0,
) -> dict:
    """Design Golden Braid L0 domestication primers for a template region.

    The amplified product, after BsaI digestion, will carry the correct 4-nt
    overhangs for the chosen `part_type` and slot directly into a GB L0
    assembly.

    Primer structure (5'→3'):

        Forward: [pad] [BsaI] [spacer] [5' overhang]    [binding →]
        Reverse: [pad] [BsaI] [spacer] [RC 3' overhang] [← binding RC]

    Returns a dict with keys: part_type, position, oh5, oh3, insert_seq,
    fwd_binding, rev_binding, fwd_full, rev_full, fwd_tm, rev_tm,
    amplicon_len.
    """
    pos_label, oh5, oh3 = _GB_POSITIONS[part_type]
    insert = template_seq[start:end].upper()

    # Forward binding: first 18-25 bp of the insert
    fwd_bind, fwd_tm = _pick_binding_region(insert, target_tm)

    # Reverse binding: first 18-25 bp of the reverse-complement of the insert
    # (i.e. the last 18-25 bp of the insert, reverse-complemented)
    rev_bind, rev_tm = _pick_binding_region(_rc(insert), target_tm)

    # Assemble full primers
    fwd_tail = _GB_PAD + _GB_BSAI_SITE + _GB_SPACER + oh5
    rev_tail = _GB_PAD + _GB_BSAI_SITE + _GB_SPACER + _rc(oh3)

    fwd_full = fwd_tail + fwd_bind
    rev_full = rev_tail + rev_bind

    # Amplicon = full fwd + insert body + full rev (minus double-counted bindings)
    amplicon_len = len(fwd_full) + (len(insert) - len(fwd_bind)) + len(rev_full) - len(rev_bind) + len(rev_bind)
    # Simpler: amplicon = pad+bsai+spacer+oh + insert + oh_rc+spacer+bsai_rc+pad
    amplicon_len = len(fwd_tail) + len(insert) + len(rev_tail)

    return {
        "part_type":   part_type,
        "position":    pos_label,
        "oh5":         oh5,
        "oh3":         oh3,
        "insert_seq":  insert,
        "fwd_binding": fwd_bind,
        "rev_binding": rev_bind,
        "fwd_full":    fwd_full,
        "rev_full":    rev_full,
        "fwd_tm":      round(fwd_tm, 1),
        "rev_tm":      round(rev_tm, 1),
        "amplicon_len": amplicon_len,
    }


# ── Parts bin persistence ─────────────────────────────────────────────────────
# User-created parts (from the domesticator) are stored in parts_bin.json next
# to the main script. Each entry is a dict with at least the 7 canonical fields
# plus sequence, primers, and Tm values.

_PARTS_BIN_FILE = Path(__file__).parent / "parts_bin.json"
_parts_bin_cache: "list | None" = None

def _load_parts_bin() -> list[dict]:
    global _parts_bin_cache
    if _parts_bin_cache is not None:
        return list(_parts_bin_cache)
    entries, warning = _safe_load_json(_PARTS_BIN_FILE, "Parts bin")
    if warning:
        _log.warning(warning)
    _parts_bin_cache = entries
    return list(_parts_bin_cache)

def _save_parts_bin(entries: list[dict]) -> None:
    global _parts_bin_cache
    _safe_save_json(_PARTS_BIN_FILE, entries, "Parts bin")
    _parts_bin_cache = list(entries)


# ── Primer design (Primer3-backed) ─────────────────────────────────────────────
#
# Two primer types:
#   Detection — diagnostic PCR primers flanking a target region.
#               Primer3 picks the thermodynamically ideal pair.
#   Cloning   — primers with restriction-enzyme tails + GCGC padding for
#               cloning a region into a new vector.
#
# Primer library persists to primers.json (same dir as plasmid_library.json).

_PRIMERS_FILE = Path(__file__).parent / "primers.json"
_primers_cache: "list | None" = None

def _load_primers() -> list[dict]:
    global _primers_cache
    if _primers_cache is not None:
        return list(_primers_cache)
    entries, warning = _safe_load_json(_PRIMERS_FILE, "Primer library")
    if warning:
        _log.warning(warning)
    _primers_cache = entries
    return list(_primers_cache)

def _save_primers(entries: list[dict]) -> None:
    global _primers_cache
    _safe_save_json(_PRIMERS_FILE, entries, "Primer library")
    _primers_cache = list(entries)


# Common cloning enzymes for the RE-site dropdown in the cloning primer panel.
# Sorted alphabetically; each tuple is (display_label, enzyme_name).
_CLONING_RE_OPTIONS: list[tuple[str, str]] = sorted([
    (f"{name}  ({site})", name)
    for name, (site, fc, rc) in _NEB_ENZYMES.items()
    if name in {
        "EcoRI", "BamHI", "XhoI", "NdeI", "NcoI", "XbaI", "SpeI", "PstI",
        "HindIII", "SalI", "NotI", "BglII", "KpnI", "SacI", "NheI", "BsaI",
        "BsmBI", "BbsI", "SapI", "AgeI", "EcoRV", "ClaI", "MfeI", "MluI",
        "NruI", "SphI", "SfiI", "AvrII", "BsiWI", "BsrGI", "BstBI",
    }
], key=lambda t: t[0])


def _design_detection_primers(
    template_seq: str,
    target_start: int,
    target_end: int,
    product_min: int = 450,
    product_max: int = 550,
    target_tm: float = 60.0,
    primer_len: int = 25,
) -> dict:
    """Design diagnostic PCR primers flanking a target region using Primer3.

    Primer3 picks the most thermodynamically ideal pair that amplifies a
    product of `product_min`..`product_max` bp spanning the target region.

    Returns a dict with keys: fwd_seq, rev_seq, fwd_tm, rev_tm, fwd_pos,
    rev_pos, product_size, or an 'error' key on failure.
    """
    import primer3
    seq = template_seq.upper()

    target_len = target_end - target_start
    if target_len < 1:
        return {"error": "Target region is empty."}

    try:
        result = primer3.design_primers(
            seq_args={
                "SEQUENCE_TEMPLATE": seq,
                "SEQUENCE_TARGET": [target_start, target_len],
            },
            global_args={
                "PRIMER_TASK": "generic",
                "PRIMER_PICK_LEFT_PRIMER": 1,
                "PRIMER_PICK_RIGHT_PRIMER": 1,
                "PRIMER_OPT_SIZE": primer_len,
                "PRIMER_MIN_SIZE": max(15, primer_len - 5),
                "PRIMER_MAX_SIZE": min(36, primer_len + 5),
                "PRIMER_OPT_TM": target_tm,
                "PRIMER_MIN_TM": target_tm - 3,
                "PRIMER_MAX_TM": target_tm + 3,
                "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
                "PRIMER_NUM_RETURN": 1,
            },
        )
    except (OSError, Exception) as exc:
        # Primer3 raises OSError on internally-inconsistent parameters
        # (e.g. PRIMER_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE).
        return {"error": f"Primer3 rejected parameters: {exc}"}

    n_found = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    if n_found == 0:
        explain = result.get("PRIMER_LEFT_EXPLAIN", "")
        return {"error": f"Primer3 found no valid pair. {explain}"}

    fwd_pos = result["PRIMER_LEFT_0"]     # (start, length)
    rev_pos = result["PRIMER_RIGHT_0"]    # (start, length) — start is 3' end

    return {
        "fwd_seq":      result["PRIMER_LEFT_0_SEQUENCE"],
        "rev_seq":      result["PRIMER_RIGHT_0_SEQUENCE"],
        "fwd_tm":       round(result["PRIMER_LEFT_0_TM"], 1),
        "rev_tm":       round(result["PRIMER_RIGHT_0_TM"], 1),
        "fwd_pos":      (fwd_pos[0], fwd_pos[0] + fwd_pos[1]),   # (start, end)
        "rev_pos":      (rev_pos[0] - rev_pos[1] + 1, rev_pos[0] + 1),
        "product_size": result["PRIMER_PAIR_0_PRODUCT_SIZE"],
    }


def _design_cloning_primers(
    template_seq: str,
    start: int,
    end: int,
    re_5prime: str,
    re_3prime: str,
    target_tm: float = 57.0,
    padding: str = "GCGC",
) -> dict:
    """Design cloning primers with restriction-enzyme tails + GCGC padding.

    Structure (5'→3'):
        Forward: [padding] [5' RE site]    [binding region →]
        Reverse: [padding] [3' RE site RC] [← binding region RC]

    Returns dict with keys: fwd_full, rev_full, fwd_binding, rev_binding,
    fwd_tm, rev_tm, re_5prime, re_3prime, insert_seq, or 'error'.
    """
    if re_5prime not in _NEB_ENZYMES:
        return {"error": f"Unknown enzyme: {re_5prime}"}
    if re_3prime not in _NEB_ENZYMES:
        return {"error": f"Unknown enzyme: {re_3prime}"}

    site_5, _, _ = _NEB_ENZYMES[re_5prime]
    site_3, _, _ = _NEB_ENZYMES[re_3prime]

    insert = template_seq[start:end].upper()
    if len(insert) < 18:
        return {"error": "Region too short (< 18 bp)."}

    fwd_bind, fwd_tm = _pick_binding_region(insert, target_tm)
    rev_bind, rev_tm = _pick_binding_region(_rc(insert), target_tm)

    fwd_full = padding + site_5 + fwd_bind
    rev_full = padding + _rc(site_3) + rev_bind

    return {
        "fwd_full":    fwd_full,
        "rev_full":    rev_full,
        "fwd_binding": fwd_bind,
        "rev_binding": rev_bind,
        "fwd_tm":      round(fwd_tm, 1),
        "rev_tm":      round(rev_tm, 1),
        "re_5prime":   re_5prime,
        "re_3prime":   re_3prime,
        "site_5":      site_5,
        "site_3":      site_3,
        "insert_seq":  insert,
        "fwd_pos":     (start, start + len(fwd_bind)),
        "rev_pos":     (end - len(rev_bind), end),
    }


# ── Parts bin modal ────────────────────────────────────────────────────────────

class PartsBinModal(Screen):
    """Golden Braid-compatible L0 parts library — full-screen view.

    Uses Screen (not ModalScreen) so it fills the terminal cleanly instead
    of floating a fixed-size box on a dark overlay. Escape or the Close
    button pops back to the main app.

    Shows both the built-in reference catalog (_GB_L0_PARTS) and user-created
    parts from parts_bin.json. User parts appear first and include sequence
    + primer data; built-in parts have no sequence (shown as "—").
    """

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(id="parts-box"):
            yield Static(" Parts Bin  —  Golden Braid L0 Parts ", id="parts-title")
            yield DataTable(id="parts-table", cursor_type="row", zebra_stripes=True)
            yield Static("", id="parts-detail")
            with Horizontal(id="parts-btns"):
                yield Button("New Part", id="btn-new-part", variant="success")
                yield Button("Close",    id="btn-parts-close")
        yield Footer()

    def _all_rows(self) -> list[dict]:
        """Combine user-created parts (first) + built-in catalog into a
        uniform list of dicts for the table and detail panel."""
        rows: list[dict] = []
        for p in _load_parts_bin():
            rows.append({
                "name":     p.get("name", "?"),
                "type":     p.get("type", "?"),
                "position": p.get("position", "?"),
                "oh5":      p.get("oh5", ""),
                "oh3":      p.get("oh3", ""),
                "backbone": p.get("backbone", "pUPD2"),
                "marker":   p.get("marker", "Spectinomycin"),
                "sequence": p.get("sequence", ""),
                "fwd_primer": p.get("fwd_primer", ""),
                "rev_primer": p.get("rev_primer", ""),
                "fwd_tm":   p.get("fwd_tm", 0.0),
                "rev_tm":   p.get("rev_tm", 0.0),
                "user":     True,
            })
        for row in _GB_L0_PARTS:
            name, ptype, pos, oh5, oh3, backbone, marker = row
            rows.append({
                "name": name, "type": ptype, "position": pos,
                "oh5": oh5, "oh3": oh3, "backbone": backbone,
                "marker": marker, "sequence": "", "fwd_primer": "",
                "rev_primer": "", "fwd_tm": 0.0, "rev_tm": 0.0,
                "user": False,
            })
        return rows

    def on_mount(self) -> None:
        t = self.query_one("#parts-table", DataTable)
        t.add_columns(
            "Name", "Type", "Pos", "5'OH", "3'OH", "Sequence",
        )
        self._populate()

    def _populate(self) -> None:
        t = self.query_one("#parts-table", DataTable)
        t.clear()
        self._rows = self._all_rows()
        for r in self._rows:
            color = _GB_TYPE_COLORS.get(r["type"], "white")
            seq_preview = r["sequence"][:28] + "…" if len(r["sequence"]) > 28 else r["sequence"]
            if not seq_preview:
                seq_preview = "—"
            usr_mark = "★ " if r["user"] else ""
            t.add_row(
                Text(usr_mark + r["name"], style=color),
                Text(r["type"], style=f"dim {color}"),
                r["position"],
                Text(r["oh5"], style="bold cyan"),
                Text(r["oh3"], style="bold cyan"),
                Text(seq_preview, style="dim color(252)"),
            )

    @on(DataTable.RowHighlighted, "#parts-table")
    def _row_highlighted(self, event: DataTable.RowHighlighted) -> None:
        idx = event.cursor_row
        if idx < 0 or not hasattr(self, "_rows") or idx >= len(self._rows):
            return
        r = self._rows[idx]
        color = _GB_TYPE_COLORS.get(r["type"], "white")
        detail = Text()
        detail.append(r["name"], style=f"bold {color}")
        detail.append(f"  [{r['type']}]", style=f"dim {color}")
        if r["user"]:
            detail.append("  ★ user part", style="dim green")
        detail.append("\n")
        detail.append(f"Position: {r['position']}   ", style="white")
        detail.append("5′ OH: ", style="dim")
        detail.append(r["oh5"], style="bold cyan")
        detail.append("   3′ OH: ", style="dim")
        detail.append(r["oh3"], style="bold cyan")
        detail.append(f"   Backbone: {r['backbone']}   Sel: {r['marker']}", style="dim")
        if r["sequence"]:
            detail.append("\n")
            detail.append(f"Sequence ({len(r['sequence'])} bp): ", style="dim")
            detail.append(r["sequence"][:60], style="color(252)")
            if len(r["sequence"]) > 60:
                detail.append("…", style="dim")
        if r["fwd_primer"]:
            detail.append("\n")
            detail.append("Fwd: ", style="dim green")
            detail.append(r["fwd_primer"], style="green")
            detail.append(f"  Tm {r['fwd_tm']:.1f}°C", style="dim")
        if r["rev_primer"]:
            detail.append("   Rev: ", style="dim red")
            detail.append(r["rev_primer"], style="red")
            detail.append(f"  Tm {r['rev_tm']:.1f}°C", style="dim")
        self.query_one("#parts-detail", Static).update(detail)

    @on(Button.Pressed, "#btn-new-part")
    def _new_part(self, _) -> None:
        # Opens the domesticator modal. The current record's sequence is
        # passed so the domesticator can use it as template.
        rec = getattr(self.app, "_current_record", None)
        seq = str(rec.seq) if rec else ""
        feats = []
        try:
            pm = self.app.query_one("#plasmid-map", PlasmidMap)
            feats = pm._feats
        except Exception:
            pass

        def _on_result(part_dict):
            if part_dict is None:
                return
            entries = _load_parts_bin()
            entries.insert(0, part_dict)
            _save_parts_bin(entries)
            self._populate()
            self.app.notify(
                f"Saved '{part_dict['name']}' to Parts Bin "
                f"({len(part_dict.get('sequence', ''))} bp).",
            )

        self.app.push_screen(
            DomesticatorModal(seq, feats),
            callback=_on_result,
        )

    @on(Button.Pressed, "#btn-parts-close")
    def _close(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Constructor modal ──────────────────────────────────────────────────────────

class DomesticatorModal(ModalScreen):
    """Golden Braid L0 Parts Domesticator.

    Takes a template sequence + region, designs domestication primers with
    the correct BsaI sites + positional overhangs, and returns a part dict
    ready for saving to the Parts Bin.

    Primer structure (5'→3'):
        Forward: GCGC GGTCTC A [5' overhang] [binding region →]
        Reverse: GCGC GGTCTC A [RC 3' OH]    [← binding region RC]

    After BsaI digestion the amplicon carries the correct 4-nt sticky ends
    for Golden Braid L0 assembly.
    """

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, template_seq: str, feats: list[dict]):
        super().__init__()
        self._template = template_seq.upper()
        self._feats    = feats   # from PlasmidMap._feats, for the feature picker
        self._design:  "dict | None" = None   # result of _design_gb_primers

    def compose(self) -> ComposeResult:
        # Build the feature dropdown: "(start-end) label" for each non-RE feature
        feat_options: list[tuple[str, str]] = []
        for f in self._feats:
            if f.get("type") in ("resite", "recut"):
                continue
            label = f.get("label", f.get("type", "?"))
            val   = f"{f['start']}-{f['end']}"
            feat_options.append((f"{label}  ({f['start']+1}‥{f['end']})", val))

        # Part-type dropdown options
        type_options = [
            (f"{k}  ({v[0]}: {v[1]}→{v[2]})", k) for k, v in _GB_POSITIONS.items()
        ]

        with Vertical(id="dom-box"):
            yield Static(
                " Domesticate Part  —  Golden Braid L0 ",
                id="dom-title",
            )
            # ── Row 1: Part name + type ──
            with Horizontal(id="dom-row1"):
                with Vertical(id="dom-name-col"):
                    yield Label("Part name")
                    yield Input(placeholder="e.g. my-promoter", id="dom-name")
                with Vertical(id="dom-type-col"):
                    yield Label("Part type")
                    yield Select(type_options, id="dom-type", value="CDS")
            # ── Row 2: overhang info (auto-updated from type) ──
            yield Static("", id="dom-oh-info", markup=True)
            # ── Row 3: template region ──
            with Horizontal(id="dom-region-row"):
                with Vertical(id="dom-feat-col"):
                    yield Label("From feature")
                    yield Select(
                        feat_options,
                        id="dom-feat",
                        prompt="(select feature or enter manually)",
                    )
                with Vertical(id="dom-start-col"):
                    yield Label("Start (bp)")
                    yield Input(placeholder="1", id="dom-start", type="integer")
                with Vertical(id="dom-end-col"):
                    yield Label("End (bp)")
                    yield Input(
                        placeholder=str(len(self._template)) if self._template else "100",
                        id="dom-end",
                        type="integer",
                    )
            # ── Primer results ──
            yield Static("", id="dom-primer-results", markup=True)
            # ── Buttons ──
            with Horizontal(id="dom-btns"):
                yield Button(
                    "Design Primers", id="btn-dom-design", variant="primary",
                )
                yield Button(
                    "Save to Parts Bin", id="btn-dom-save", variant="success",
                    disabled=True,
                )
                yield Button("Cancel", id="btn-dom-cancel")

    def on_mount(self) -> None:
        self._update_oh_display()
        # Focus the name input
        self.query_one("#dom-name", Input).focus()

    # ── Feature selection fills start/end ──────────────────────────────────

    @on(Select.Changed, "#dom-feat")
    def _feat_selected(self, event: Select.Changed) -> None:
        val = event.value
        if not isinstance(val, str):
            return
        if "-" in val:
            parts = val.split("-", 1)
            try:
                self.query_one("#dom-start", Input).value = str(int(parts[0]) + 1)
                self.query_one("#dom-end",   Input).value = parts[1]
            except ValueError:
                pass

    # ── Part type changes update the overhang info ─────────────────────────

    @on(Select.Changed, "#dom-type")
    def _type_changed(self, _event) -> None:
        self._update_oh_display()

    def _update_oh_display(self) -> None:
        sel = self.query_one("#dom-type", Select)
        val = sel.value
        if not isinstance(val, str) or val not in _GB_POSITIONS:
            self.query_one("#dom-oh-info", Static).update("")
            return
        pos, oh5, oh3 = _GB_POSITIONS[val]
        self.query_one("#dom-oh-info", Static).update(
            f"  [dim]{pos}[/dim]   "
            f"5′ overhang: [bold cyan]{oh5}[/bold cyan]   →   "
            f"3′ overhang: [bold cyan]{oh3}[/bold cyan]   "
            f"[dim](BsaI domestication)[/dim]"
        )

    # ── Design primers ─────────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-dom-design")
    def _design(self, _) -> None:
        status = self.query_one("#dom-primer-results", Static)
        # Validate inputs
        part_type = self.query_one("#dom-type", Select).value
        if not isinstance(part_type, str) or part_type not in _GB_POSITIONS:
            status.update("[red]Select a part type.[/red]")
            return
        if not self._template:
            status.update("[red]No template sequence loaded.[/red]")
            return
        try:
            start = int(self.query_one("#dom-start", Input).value) - 1  # 1-based → 0-based
            end   = int(self.query_one("#dom-end",   Input).value)
        except ValueError:
            status.update("[red]Enter valid start and end positions.[/red]")
            return
        if start < 0 or end <= start or end > len(self._template):
            status.update(
                f"[red]Invalid region: {start+1}–{end} "
                f"(plasmid is {len(self._template)} bp)[/red]"
            )
            return
        if end - start < 20:
            status.update("[red]Region too short (< 20 bp).[/red]")
            return

        try:
            self._design = _design_gb_primers(self._template, start, end, part_type)
        except Exception as exc:
            _log.exception("Primer design failed")
            status.update(f"[red]Primer design failed: {exc}[/red]")
            return

        d = self._design
        t = Text()
        t.append("── Primers designed ─────────────────────────────────\n",
                 style="dim")
        t.append("\nForward (5'→3'):\n", style="bold green")
        # Show with structure annotation
        tail_len = len(_GB_PAD + _GB_BSAI_SITE + _GB_SPACER) + 4  # pad+bsai+spacer+oh
        t.append(f"  {d['fwd_full'][:tail_len]}", style="dim green")
        t.append(d["fwd_full"][tail_len:], style="bold green")
        t.append(f"   Tm {d['fwd_tm']:.1f}°C\n", style="dim")
        t.append(f"  {'─'*4}{'BsaI──':>7}{'─OH':>3}{'─── binding region':>20}\n",
                 style="dim")
        t.append("\nReverse (5'→3'):\n", style="bold red")
        t.append(f"  {d['rev_full'][:tail_len]}", style="dim red")
        t.append(d["rev_full"][tail_len:], style="bold red")
        t.append(f"   Tm {d['rev_tm']:.1f}°C\n", style="dim")
        t.append(f"  {'─'*4}{'BsaI──':>7}{'─OH':>3}{'─── binding region':>20}\n",
                 style="dim")
        t.append(f"\nInsert: {len(d['insert_seq'])} bp   "
                 f"Amplicon: {d['amplicon_len']} bp\n",
                 style="white")
        status.update(t)
        self.query_one("#btn-dom-save", Button).disabled = False

    # ── Save to parts bin ──────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-dom-save")
    def _save(self, _) -> None:
        if self._design is None:
            return
        name = self.query_one("#dom-name", Input).value.strip()
        if not name:
            self.query_one("#dom-primer-results", Static).update(
                "[red]Enter a part name before saving.[/red]"
            )
            return
        d = self._design
        part = {
            "name":        name,
            "type":        d["part_type"],
            "position":    d["position"],
            "oh5":         d["oh5"],
            "oh3":         d["oh3"],
            "backbone":    "pUPD2",
            "marker":      "Spectinomycin",
            "sequence":    d["insert_seq"],
            "fwd_primer":  d["fwd_full"],
            "rev_primer":  d["rev_full"],
            "fwd_tm":      d["fwd_tm"],
            "rev_tm":      d["rev_tm"],
        }
        self.dismiss(part)

    @on(Button.Pressed, "#btn-dom-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class ConstructorModal(ModalScreen):
    """Golden Braid TU Constructor — assemble L0 parts into a transcription unit."""

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    # L1 destination backbone info
    _BACKBONES: dict = {
        "Alpha1": {"id": "pDGB1_alpha1", "selection": "Spectinomycin", "note": "L1 alpha orientation"},
        "Alpha2": {"id": "pDGB1_alpha2", "selection": "Spectinomycin", "note": "L1 alpha orientation"},
        "Omega1": {"id": "pDGB1_omega1", "selection": "Kanamycin",     "note": "L1 omega orientation"},
        "Omega2": {"id": "pDGB1_omega2", "selection": "Kanamycin",     "note": "L1 omega orientation"},
    }

    # Golden Braid L1 boundary overhangs
    _TU_START = "GGAG"
    _TU_END   = "CGCT"

    # Part types that occupy each positional slot (for duplicate detection)
    _POS_SLOT: dict = {
        "Promoter":   1,
        "5' UTR":     2,
        "CDS":        3,
        "CDS-NS":     3,
        "C-tag":      4,
        "Terminator": 5,
    }

    def __init__(self) -> None:
        super().__init__()
        self._lane:     list[tuple] = []
        self._backbone: str         = "Alpha1"

    # ── Layout ────────────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:
        with Vertical(id="ctor-box"):
            yield Static(
                " Constructor  —  Golden Braid TU Assembly ", id="ctor-title"
            )
            with Horizontal(id="ctor-main"):
                # Left: parts palette
                with Vertical(id="ctor-palette-col"):
                    yield Static(" Parts Palette ", id="ctor-palette-hdr")
                    yield DataTable(
                        id="ctor-palette", cursor_type="row", zebra_stripes=True
                    )
                    yield Button(
                        "→  Add to Lane", id="btn-ctor-add", variant="primary"
                    )
                # Right: assembly lane
                with Vertical(id="ctor-lane-col"):
                    yield Static(" Assembly Lane ", id="ctor-lane-hdr")
                    yield DataTable(
                        id="ctor-lane", cursor_type="row", zebra_stripes=True
                    )
                    with Horizontal(id="ctor-lane-btns"):
                        yield Button("↑",        id="btn-lane-up")
                        yield Button("↓",        id="btn-lane-down")
                        yield Button("✕ Remove", id="btn-lane-remove", variant="error")
            # Backbone selector
            with Horizontal(id="ctor-backbone-row"):
                yield Static("Backbone:", id="ctor-backbone-label")
                for bb in self._BACKBONES:
                    classes = "bb-btn bb-active" if bb == self._backbone else "bb-btn"
                    yield Button(bb, id=f"btn-bb-{bb}", classes=classes)
            # Overhang chain + validation messages
            yield Static("", id="ctor-validation")
            # Bottom actions
            with Horizontal(id="ctor-btns"):
                yield Button(
                    "Simulate Assembly", id="btn-ctor-simulate",
                    variant="success", disabled=True
                )
                yield Button("Clear Lane", id="btn-ctor-clear", variant="warning")
                yield Button("Close",      id="btn-ctor-close")

    def on_mount(self) -> None:
        # Populate palette
        pt = self.query_one("#ctor-palette", DataTable)
        pt.add_columns("Name", "Type", "Pos", "5' OH", "3' OH")
        for row in _GB_L0_PARTS:
            name, ptype, pos, oh5, oh3, *_ = row
            color = _GB_TYPE_COLORS.get(ptype, "white")
            pt.add_row(
                Text(name,  style=color),
                Text(ptype, style=f"dim {color}"),
                pos,
                Text(oh5,   style="bold cyan"),
                Text(oh3,   style="bold cyan"),
            )
        # Set up lane columns
        lt = self.query_one("#ctor-lane", DataTable)
        lt.add_columns("#", "Name", "Type", "5' OH", "3' OH")
        self._refresh_validation()

    # ── Lane management ───────────────────────────────────────────────────────

    def _refresh_lane(self, restore_cursor: int = -1) -> None:
        lt = self.query_one("#ctor-lane", DataTable)
        lt.clear()
        for i, row in enumerate(self._lane):
            name, ptype, pos, oh5, oh3, *_ = row
            color = _GB_TYPE_COLORS.get(ptype, "white")
            lt.add_row(
                str(i + 1),
                Text(name,  style=color),
                Text(ptype, style=f"dim {color}"),
                Text(oh5,   style="bold cyan"),
                Text(oh3,   style="bold cyan"),
            )
        if restore_cursor >= 0 and restore_cursor < len(self._lane):
            try:
                lt.move_cursor(row=restore_cursor)
            except Exception:
                pass

    def _add_selected_part(self) -> None:
        pt  = self.query_one("#ctor-palette", DataTable)
        idx = pt.cursor_row
        if 0 <= idx < len(_GB_L0_PARTS):
            self._lane.append(_GB_L0_PARTS[idx])
            self._refresh_lane(restore_cursor=len(self._lane) - 1)
            self._refresh_validation()

    # ── Grammar validation ────────────────────────────────────────────────────

    def _validate(self) -> tuple[bool, list[str]]:
        """Return (is_valid, error_list). Valid = complete, correctly-chained TU."""
        if not self._lane:
            return False, ["Lane is empty — add L0 parts to build a TU."]

        errors: list[str] = []

        # 1. Boundary overhangs
        if self._lane[0][3] != self._TU_START:
            errors.append(
                f"First part must carry the {self._TU_START} overhang "
                f"(Promoter, Pos 1). Got {self._lane[0][3]!r}."
            )
        if self._lane[-1][4] != self._TU_END:
            errors.append(
                f"Last part must carry the {self._TU_END} overhang "
                f"(Terminator, Pos 5). Got {self._lane[-1][4]!r}."
            )

        # 2. Overhang continuity
        for i in range(len(self._lane) - 1):
            oh3 = self._lane[i][4]
            oh5 = self._lane[i + 1][3]
            if oh3 != oh5:
                errors.append(
                    f"Overhang mismatch at junction {i+1}→{i+2}: "
                    f"{self._lane[i][0]!r} ends {oh3!r} but "
                    f"{self._lane[i+1][0]!r} starts {oh5!r}."
                )

        # 3. Duplicate positional slots
        seen: dict[int, str] = {}
        for row in self._lane:
            name, ptype = row[0], row[1]
            slot = self._POS_SLOT.get(ptype)
            if slot is not None:
                if slot in seen:
                    errors.append(
                        f"Slot {slot} occupied twice: {seen[slot]!r} and {name!r}."
                    )
                else:
                    seen[slot] = name

        # 4. CDS-NS ↔ C-tag pairing
        for i, row in enumerate(self._lane):
            if row[1] == "CDS-NS":
                nxt = self._lane[i + 1][1] if i + 1 < len(self._lane) else None
                if nxt != "C-tag":
                    errors.append(
                        f"{row[0]!r} has no stop codon — must be immediately "
                        f"followed by a C-terminal tag (Pos 4)."
                    )
            elif row[1] == "C-tag":
                prv = self._lane[i - 1][1] if i > 0 else None
                if prv != "CDS-NS":
                    errors.append(
                        f"{row[0]!r} (C-tag, Pos 4) must follow a no-stop CDS "
                        f"(Pos 3). Found: {prv!r}."
                    )

        # 5. Mandatory parts
        types = {r[1] for r in self._lane}
        if "Promoter" not in types:
            errors.append("Missing Promoter (Pos 1, 5' OH: GGAG).")
        if "CDS" not in types and not ("CDS-NS" in types and "C-tag" in types):
            errors.append(
                "Missing CDS (Pos 3-4). Add a CDS, or a CDS-NS + C-tag pair."
            )
        if "Terminator" not in types:
            errors.append("Missing Terminator (Pos 5, 3' OH: CGCT).")

        return len(errors) == 0, errors

    def _build_chain(self) -> Text:
        """Render the overhang chain with colour-coded junctions."""
        t = Text()
        if not self._lane:
            t.append("(empty)", style="dim")
            return t

        # Opening backbone overhang
        start_ok = (self._lane[0][3] == self._TU_START)
        t.append("5'-", style="dim")
        t.append(self._TU_START, style="bold green" if start_ok else "bold red")

        for i, row in enumerate(self._lane):
            name, ptype, pos, oh5, oh3, *_ = row
            color   = _GB_TYPE_COLORS.get(ptype, "white")
            # incoming junction colour
            exp_in  = self._TU_START if i == 0 else self._lane[i - 1][4]
            junc_ok = (oh5 == exp_in)
            dash    = "—" if junc_ok else "≠"
            t.append(dash, style="white" if junc_ok else "bold red")
            t.append(f"[{name}]", style=color)
            t.append("—", style="white")
            # outgoing OH colour
            exp_out  = self._lane[i + 1][3] if i + 1 < len(self._lane) else self._TU_END
            oh3_ok   = (oh3 == exp_out)
            t.append(oh3, style="bold cyan" if oh3_ok else "bold red")

        t.append("-3'", style="dim")
        return t

    def _refresh_validation(self) -> None:
        is_valid, errors = self._validate()
        bb     = self._BACKBONES[self._backbone]
        vbox   = self.query_one("#ctor-validation", Static)
        sim    = self.query_one("#btn-ctor-simulate", Button)
        sim.disabled = not is_valid

        t = Text()
        t.append_text(self._build_chain())
        t.append("\n")
        if is_valid:
            t.append(
                f"✓  Valid TU — assembles into {self._backbone} "
                f"({bb['id']}, {bb['selection']} selection, {bb['note']})",
                style="bold green",
            )
        else:
            for err in errors:
                t.append(f"✗  {err}\n", style="bold red")
        vbox.update(t)

    # ── Button handlers ───────────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-ctor-add")
    def _on_add(self, _) -> None:
        self._add_selected_part()

    @on(Button.Pressed, "#btn-lane-up")
    def _on_up(self, _) -> None:
        lt  = self.query_one("#ctor-lane", DataTable)
        idx = lt.cursor_row
        if idx <= 0 or idx >= len(self._lane):
            return
        self._lane[idx - 1], self._lane[idx] = self._lane[idx], self._lane[idx - 1]
        self._refresh_lane(restore_cursor=idx - 1)
        self._refresh_validation()

    @on(Button.Pressed, "#btn-lane-down")
    def _on_down(self, _) -> None:
        lt  = self.query_one("#ctor-lane", DataTable)
        idx = lt.cursor_row
        if idx < 0 or idx >= len(self._lane) - 1:
            return
        self._lane[idx], self._lane[idx + 1] = self._lane[idx + 1], self._lane[idx]
        self._refresh_lane(restore_cursor=idx + 1)
        self._refresh_validation()

    @on(Button.Pressed, "#btn-lane-remove")
    def _on_remove(self, _) -> None:
        lt  = self.query_one("#ctor-lane", DataTable)
        idx = lt.cursor_row
        if 0 <= idx < len(self._lane):
            self._lane.pop(idx)
            self._refresh_lane(restore_cursor=min(idx, len(self._lane) - 1))
            self._refresh_validation()

    @on(Button.Pressed, ".bb-btn")
    def _on_backbone(self, event: Button.Pressed) -> None:
        bb = (event.button.id or "").replace("btn-bb-", "")
        if bb not in self._BACKBONES:
            return
        self._backbone = bb
        for name in self._BACKBONES:
            btn = self.query_one(f"#btn-bb-{name}", Button)
            btn.set_class(name == bb, "bb-active")
        self._refresh_validation()

    @on(Button.Pressed, "#btn-ctor-simulate")
    def _on_simulate(self, _) -> None:
        self.app.notify("Simulate Assembly: coming soon.", severity="information")

    @on(Button.Pressed, "#btn-ctor-clear")
    def _on_clear(self, _) -> None:
        self._lane.clear()
        self._refresh_lane()
        self._refresh_validation()

    @on(Button.Pressed, "#btn-ctor-close")
    def _on_close(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Primer design screen (full-screen) ─────────────────────────────────────────

class PrimerDesignScreen(Screen):
    """Full-screen Primer3-backed primer design workbench.

    Two workflows:
      Detection — diagnostic PCR primers (Primer3 picks ideal pair).
      Cloning   — primers with restriction-enzyme tails + GCGC padding.

    Designed primers are saved to primers.json and optionally added as
    primer_bind features to the currently-loaded plasmid.
    """

    BINDINGS = [Binding("escape", "cancel", "Close")]

    def __init__(self, template_seq: str, feats: list[dict],
                 plasmid_name: str = ""):
        super().__init__()
        self._template     = template_seq.upper()
        self._feats        = feats
        self._plasmid_name = plasmid_name
        self._det_result:  "dict | None" = None
        self._clo_result:  "dict | None" = None

    def compose(self) -> ComposeResult:
        # Feature dropdown
        feat_opts: list[tuple[str, str]] = []
        for f in self._feats:
            if f.get("type") in ("resite", "recut"):
                continue
            label = f.get("label", f.get("type", "?"))
            feat_opts.append(
                (f"{label}  ({f['start']+1}‥{f['end']})", f"{f['start']}-{f['end']}")
            )
        # RE site dropdown
        re_opts = _CLONING_RE_OPTIONS

        yield Header()
        with Vertical(id="pd-box"):
            yield Static(" Primer Design  —  Primer3 ", id="pd-title")

            # ── Source / region row ────────────────────────────────────────
            with Horizontal(id="pd-source-row"):
                with Vertical(id="pd-feat-col"):
                    yield Label("Feature")
                    yield Select(feat_opts, id="pd-feat",
                                 prompt="(select or enter manually)")
                with Vertical(id="pd-start-col"):
                    yield Label("Start")
                    yield Input(placeholder="1", id="pd-start", type="integer")
                with Vertical(id="pd-end-col"):
                    yield Label("End")
                    yield Input(
                        placeholder=str(len(self._template)) if self._template else "",
                        id="pd-end", type="integer")
                with Vertical(id="pd-name-col"):
                    yield Label("Part name")
                    yield Input(value=self._plasmid_name, id="pd-part-name")

            # ── Detection primers ──────────────────────────────────────────
            yield Static(
                " [bold]Detection Primers[/bold]  [dim](diagnostic PCR)[/dim]",
                id="pd-det-hdr", markup=True,
            )
            with Horizontal(id="pd-det-row"):
                yield Label("Product ")
                yield Input(value="450", id="pd-det-min", type="integer")
                yield Label("–")
                yield Input(value="550", id="pd-det-max", type="integer")
                yield Label(" bp   Tm ")
                yield Input(value="60", id="pd-det-tm", type="integer")
                yield Label("°C   Len ")
                yield Input(value="25", id="pd-det-len", type="integer")
                yield Label(" bp")
                yield Button("Design Detection", id="btn-det-design",
                             variant="primary")

            # ── Cloning primers ────────────────────────────────────────────
            yield Static(
                " [bold]Cloning Primers[/bold]  [dim](RE tails + GCGC padding)[/dim]",
                id="pd-clo-hdr", markup=True,
            )
            with Horizontal(id="pd-clo-row"):
                with Vertical(id="pd-clo-5col"):
                    yield Label("5' RE site")
                    yield Select(re_opts, id="pd-re5", value="EcoRI")
                with Vertical(id="pd-clo-3col"):
                    yield Label("3' RE site")
                    yield Select(re_opts, id="pd-re3", value="BamHI")
                with Vertical(id="pd-clo-tmcol"):
                    yield Label("Binding Tm")
                    yield Input(value="57", id="pd-clo-tm", type="integer")
                yield Button("Design Cloning", id="btn-clo-design",
                             variant="primary")

            # ── Results panel ──────────────────────────────────────────────
            yield Static("", id="pd-results", markup=True)
            with Horizontal(id="pd-result-names"):
                with Vertical(id="pd-fn-col"):
                    yield Label("Fwd name")
                    yield Input(id="pd-fwd-name", placeholder="fwd primer name")
                with Vertical(id="pd-rn-col"):
                    yield Label("Rev name")
                    yield Input(id="pd-rev-name", placeholder="rev primer name")

            # ── Action buttons ─────────────────────────────────────────────
            with Horizontal(id="pd-btns"):
                yield Button("Save to Primer Library", id="btn-pd-save",
                             variant="success", disabled=True)
                yield Button("Add as Features to Map", id="btn-pd-feat",
                             variant="warning", disabled=True)
                yield Button("Close", id="btn-pd-close")

            # ── Primer library table ───────────────────────────────────────
            yield Static(" Primer Library ", id="pd-lib-hdr")
            yield DataTable(id="pd-lib-table", cursor_type="row",
                            zebra_stripes=True)
        yield Footer()

    def on_mount(self) -> None:
        t = self.query_one("#pd-lib-table", DataTable)
        t.add_columns("Name", "Sequence", "Tm", "Type", "Source")
        self._refresh_library_table()

    def _refresh_library_table(self) -> None:
        t = self.query_one("#pd-lib-table", DataTable)
        t.clear()
        for p in _load_primers():
            t.add_row(
                Text(p.get("name", "?"), style="bold"),
                Text(p.get("sequence", "")[:40], style="dim color(252)"),
                f"{p.get('tm', 0):.1f}°C",
                p.get("primer_type", "?"),
                p.get("source", ""),
            )

    # ── Feature dropdown → fill start/end ──────────────────────────────────

    @on(Select.Changed, "#pd-feat")
    def _feat_selected(self, event: Select.Changed) -> None:
        val = event.value
        if not isinstance(val, str) or "-" not in val:
            return
        parts = val.split("-", 1)
        try:
            self.query_one("#pd-start", Input).value = str(int(parts[0]) + 1)
            self.query_one("#pd-end", Input).value = parts[1]
            # Auto-set part name from the feature label
            for f in self._feats:
                if f"{f['start']}-{f['end']}" == val:
                    name_inp = self.query_one("#pd-part-name", Input)
                    if not name_inp.value.strip():
                        name_inp.value = f.get("label", "")
                    break
        except ValueError:
            pass

    # ── Helpers ────────────────────────────────────────────────────────────

    def _read_region(self) -> "tuple[int, int, str] | None":
        """Read and validate start/end/part-name from the inputs.
        Returns (start_0based, end, part_name) or None after notifying."""
        try:
            start = int(self.query_one("#pd-start", Input).value) - 1
            end   = int(self.query_one("#pd-end", Input).value)
        except ValueError:
            self.app.notify("Enter valid start and end positions.", severity="error")
            return None
        if start < 0 or end <= start or end > len(self._template):
            self.app.notify(
                f"Invalid region: {start+1}–{end} "
                f"(sequence is {len(self._template)} bp).", severity="error")
            return None
        name = self.query_one("#pd-part-name", Input).value.strip() or "primer"
        return start, end, name

    def _show_result(self, design: dict, primer_type: str,
                     fwd_key: str, rev_key: str) -> None:
        """Display a primer pair in the results panel and fill the name
        inputs with the default naming scheme."""
        status = self.query_one("#pd-results", Static)
        t = Text()
        t.append("Forward (5'→3'):\n", style="bold green")
        t.append(f"  {design[fwd_key]}\n", style="green")
        t.append(f"  Tm {design['fwd_tm']:.1f}°C   "
                 f"{len(design[fwd_key])} nt\n", style="dim")
        t.append("Reverse (5'→3'):\n", style="bold red")
        t.append(f"  {design[rev_key]}\n", style="red")
        t.append(f"  Tm {design['rev_tm']:.1f}°C   "
                 f"{len(design[rev_key])} nt\n", style="dim")
        if "product_size" in design:
            t.append(f"Product: {design['product_size']} bp\n", style="white")
        if "re_5prime" in design:
            t.append(
                f"5' RE: {design['re_5prime']} ({design['site_5']})   "
                f"3' RE: {design['re_3prime']} ({design['site_3']})\n",
                style="cyan",
            )
        status.update(t)

        name = self.query_one("#pd-part-name", Input).value.strip() or "primer"
        suffix = "DET" if primer_type == "detection" else "CLO"
        self.query_one("#pd-fwd-name", Input).value = f"{name}-{suffix}-F"
        self.query_one("#pd-rev-name", Input).value = f"{name}-{suffix}-R"
        self.query_one("#btn-pd-save", Button).disabled = False
        self.query_one("#btn-pd-feat", Button).disabled = False

    # ── Detection design ───────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-det-design")
    def _design_detection(self, _) -> None:
        region = self._read_region()
        if region is None:
            return
        start, end, name = region
        self._clo_result = None
        try:
            p_min = int(self.query_one("#pd-det-min", Input).value)
            p_max = int(self.query_one("#pd-det-max", Input).value)
            tm    = float(self.query_one("#pd-det-tm", Input).value)
            plen  = int(self.query_one("#pd-det-len", Input).value)
        except ValueError:
            self.app.notify("Invalid detection primer parameters.", severity="error")
            return
        result = _design_detection_primers(
            self._template, start, end,
            product_min=p_min, product_max=p_max,
            target_tm=tm, primer_len=plen,
        )
        if "error" in result:
            self.query_one("#pd-results", Static).update(
                f"[red]{result['error']}[/red]")
            self._det_result = None
            return
        self._det_result = result
        self._det_result["_type"] = "detection"
        self._show_result(result, "detection", "fwd_seq", "rev_seq")

    # ── Cloning design ─────────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-clo-design")
    def _design_cloning(self, _) -> None:
        region = self._read_region()
        if region is None:
            return
        start, end, name = region
        self._det_result = None
        re5 = self.query_one("#pd-re5", Select).value
        re3 = self.query_one("#pd-re3", Select).value
        if not isinstance(re5, str) or not isinstance(re3, str):
            self.app.notify("Select both restriction enzymes.", severity="error")
            return
        try:
            tm = float(self.query_one("#pd-clo-tm", Input).value)
        except ValueError:
            tm = 57.0
        result = _design_cloning_primers(
            self._template, start, end, re5, re3, target_tm=tm,
        )
        if "error" in result:
            self.query_one("#pd-results", Static).update(
                f"[red]{result['error']}[/red]")
            self._clo_result = None
            return
        self._clo_result = result
        self._clo_result["_type"] = "cloning"
        self._show_result(result, "cloning", "fwd_full", "rev_full")

    # ── Save to primer library ─────────────────────────────────────────────

    @on(Button.Pressed, "#btn-pd-save")
    def _save_primers(self, _) -> None:
        result = self._det_result or self._clo_result
        if result is None:
            return
        fwd_name = self.query_one("#pd-fwd-name", Input).value.strip()
        rev_name = self.query_one("#pd-rev-name", Input).value.strip()
        if not fwd_name or not rev_name:
            self.app.notify("Enter primer names before saving.", severity="error")
            return
        ptype = result.get("_type", "?")
        fwd_key = "fwd_seq" if ptype == "detection" else "fwd_full"
        rev_key = "rev_seq" if ptype == "detection" else "rev_full"

        entries = _load_primers()
        for pname, seq, tm, pos in [
            (fwd_name, result[fwd_key], result["fwd_tm"], result["fwd_pos"]),
            (rev_name, result[rev_key], result["rev_tm"], result["rev_pos"]),
        ]:
            # Remove any existing primer with the same name (update in place)
            entries = [e for e in entries if e.get("name") != pname]
            entries.insert(0, {
                "name":        pname,
                "sequence":    seq,
                "tm":          tm,
                "primer_type": ptype,
                "source":      self._plasmid_name,
                "pos_start":   pos[0],
                "pos_end":     pos[1],
                "strand":      1 if pname.endswith("-F") else -1,
            })
        _save_primers(entries)
        self._refresh_library_table()
        self.app.notify(f"Saved {fwd_name} + {rev_name} to primer library.")

    # ── Add as features to the map ─────────────────────────────────────────

    @on(Button.Pressed, "#btn-pd-feat")
    def _add_features(self, _) -> None:
        result = self._det_result or self._clo_result
        if result is None:
            return
        rec = getattr(self.app, "_current_record", None)
        if rec is None:
            self.app.notify("No plasmid loaded.", severity="warning")
            return
        fwd_name = self.query_one("#pd-fwd-name", Input).value.strip() or "fwd"
        rev_name = self.query_one("#pd-rev-name", Input).value.strip() or "rev"
        fwd_pos  = result["fwd_pos"]
        rev_pos  = result["rev_pos"]

        from Bio.SeqFeature import SeqFeature, FeatureLocation

        rec.features.append(SeqFeature(
            FeatureLocation(fwd_pos[0], fwd_pos[1], strand=1),
            type="primer_bind",
            qualifiers={"label": [fwd_name]},
        ))
        rec.features.append(SeqFeature(
            FeatureLocation(rev_pos[0], rev_pos[1], strand=-1),
            type="primer_bind",
            qualifiers={"label": [rev_name]},
        ))
        # Re-apply the record to refresh map + sidebar + sequence panel
        try:
            self.app._apply_record(rec)
            # Also save the updated record back to the library
            lib = self.app.query_one("#library")
            lib.add_entry(rec)
        except Exception:
            _log.exception("Failed to add primer features to map")
        self.app.notify(
            f"Added {fwd_name} + {rev_name} as primer_bind features.")

    # ── Close ──────────────────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-pd-close")
    def _close(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
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


class RenamePlasmidModal(ModalScreen):
    """Prompt for a new name for a library entry.

    Dismisses with the new name (a non-empty string) or None on cancel.
    Input validation (non-empty, trimmed, collision check) lives in the
    app-side handler — the modal just collects a value.
    """

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, current_name: str, entry_id: str):
        super().__init__()
        self.current_name = current_name
        self.entry_id     = entry_id

    def compose(self) -> ComposeResult:
        with Vertical(id="rename-dlg"):
            yield Static(" Rename plasmid ", id="rename-title")
            yield Label(f"Current name:  {self.current_name}")
            yield Label("New name:")
            yield Input(
                value=self.current_name,
                placeholder="enter a new name",
                id="rename-input",
            )
            yield Static("", id="rename-status", markup=True)
            with Horizontal(id="rename-btns"):
                yield Button("Save",   id="btn-rename-save",   variant="primary")
                yield Button("Cancel", id="btn-rename-cancel")

    def on_mount(self) -> None:
        # Default focus on the Input, text pre-selected via select_on_focus
        # (Textual Input defaults to selecting all when focused, which is
        # what you want for a rename — typing replaces the old name).
        inp = self.query_one("#rename-input", Input)
        inp.focus()

    @on(Button.Pressed, "#btn-rename-save")
    def _save(self, _):
        self._try_submit()

    @on(Input.Submitted, "#rename-input")
    def _submitted(self, _):
        self._try_submit()

    def _try_submit(self) -> None:
        new_name = self.query_one("#rename-input", Input).value.strip()
        status   = self.query_one("#rename-status", Static)
        if not new_name:
            status.update("[red]Name cannot be empty.[/red]")
            return
        if new_name == self.current_name:
            # No-op rename — treat as cancel so the app doesn't bother writing.
            self.dismiss(None)
            return
        self.dismiss(new_name)

    @on(Button.Pressed, "#btn-rename-cancel")
    def _cancel_btn(self, _):
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class LibraryDeleteConfirmModal(ModalScreen):
    """Shown when the user presses Delete while focused on the library panel.

    Default focus is on the [No] button — the Delete key is easy to press
    by accident after tabbing into the library, so we require a deliberate
    left-arrow + Enter (or Tab + Enter) to confirm the removal. Escape or
    No dismisses with False; Yes dismisses with True.
    """

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, name: str, size: int, entry_id: str):
        super().__init__()
        self.entry_name = name
        self.entry_size = size
        self.entry_id   = entry_id

    def compose(self) -> ComposeResult:
        with Vertical(id="libdel-dlg"):
            yield Static(" Remove from library ", id="libdel-title")
            yield Static(
                f"  Remove [bold]{self.entry_name}[/bold] "
                f"({self.entry_size:,} bp) from the library?\n\n"
                f"  [dim]This removes the stored GenBank copy from "
                f"plasmid_library.json.\n"
                f"  The current on-screen record is not affected.[/dim]",
                id="libdel-msg",
                markup=True,
            )
            with Horizontal(id="libdel-btns"):
                yield Button("No",           id="btn-libdel-no",  variant="primary")
                yield Button("Yes, remove",  id="btn-libdel-yes", variant="error")

    def on_mount(self) -> None:
        # Default focus on "No" — the whole point of this dialog is to catch
        # a handslip, so Enter should do nothing destructive.
        self.query_one("#btn-libdel-no", Button).focus()

    @on(Button.Pressed, "#btn-libdel-no")
    def _no(self, _):
        self.dismiss(False)

    @on(Button.Pressed, "#btn-libdel-yes")
    def _yes(self, _):
        self.dismiss(True)

    def action_cancel(self) -> None:
        self.dismiss(False)


# ── Main app ───────────────────────────────────────────────────────────────────

class PlasmidApp(App):
    TITLE       = "SpliceCraft"
    TRANSITIONS = {}          # instant screen open/close — no slide animations
    _preload_record = None
    _current_record = None   # last-loaded SeqRecord
    _source_path:   "str | None" = None   # file the current record was loaded from
    _unsaved:        bool         = False  # True when there are unsaved edits
    _MAX_UNDO = 50
    _restr_unique_only: bool = True
    _restr_min_len: int = 6
    _show_restr: bool = True
    _restr_cache: "list" = []

    CSS = """
Screen { background: $background; }

/* ── Layout ─────────────────────────────────────────────── */
MenuBar { height: 1; dock: top; }
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

/* ── Library-delete confirmation ─────────────────────────── */
LibraryDeleteConfirmModal { align: center middle; }
#libdel-dlg {
    width: 64; height: auto;
    background: $surface; border: solid $error; padding: 1 2;
}
#libdel-title { background: $error-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#libdel-msg   { color: $text-muted; margin-bottom: 1; }
#libdel-btns  { height: 3; margin-top: 1; }
#libdel-btns Button { margin-right: 1; min-width: 14; }

/* ── Rename plasmid dialog ───────────────────────────────── */
RenamePlasmidModal { align: center middle; }
#rename-dlg {
    width: 60; height: auto;
    background: $surface; border: solid $primary; padding: 1 2;
}
#rename-title  { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#rename-input  { margin-top: 0; margin-bottom: 1; }
#rename-status { height: 1; color: $text-muted; }
#rename-btns   { height: 3; margin-top: 1; }
#rename-btns Button { margin-right: 1; }

/* ── Constructor modal ───────────────────────────────────── */
ConstructorModal { align: center middle; }
#ctor-box {
    width: 116; height: 42;
    background: $surface; border: solid $accent; padding: 1 2;
}
#ctor-title       { background: $accent-darken-2; padding: 0 1; margin-bottom: 1; }
#ctor-main        { height: 20; }
#ctor-palette-col { width: 1fr; border-right: solid $primary-darken-2; padding-right: 1; }
#ctor-palette-hdr { background: $primary-darken-2; padding: 0 1; }
#ctor-palette     { height: 1fr; }
#ctor-lane-col    { width: 1fr; padding-left: 1; }
#ctor-lane-hdr    { background: $primary-darken-2; padding: 0 1; }
#ctor-lane        { height: 1fr; }
#ctor-lane-btns   { height: 3; margin-top: 0; }
#ctor-lane-btns Button { min-width: 5; margin-right: 1; }
#ctor-backbone-row { height: 3; margin-top: 1; align: left middle; }
#ctor-backbone-label { width: auto; padding: 0 1; color: $text-muted; }
.bb-btn           { min-width: 9; margin-right: 1; }
.bb-active        { background: $accent; color: $text; }
#ctor-validation  {
    height: 5; border: solid $primary-darken-2;
    padding: 0 1; margin-top: 1; overflow-x: auto;
}
#ctor-btns        { height: 3; margin-top: 1; }
#ctor-btns Button { margin-right: 1; }

/* ── Parts bin (full-screen) ─────────────────────────────── */
#parts-box {
    width: 100%; height: 1fr;
    background: $surface; padding: 0 2;
}
#parts-title  { background: $success-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#parts-table  { height: 1fr; }
#parts-detail { height: 7; border-top: solid $accent; padding: 0 1; color: $text-muted; }
#parts-btns   { height: 3; margin-top: 1; }
#parts-btns Button { margin-right: 1; }

/* ── Domesticator modal ─────────────────────────────────── */
DomesticatorModal { align: center middle; }
#dom-box {
    width: 100; height: auto; max-height: 42;
    background: $surface; border: solid $accent; padding: 1 2;
}
#dom-title  { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#dom-row1   { height: 5; }
#dom-name-col { width: 1fr; padding-right: 1; }
#dom-type-col { width: 1fr; }
#dom-oh-info  { height: 1; margin-bottom: 1; }
#dom-region-row { height: 5; }
#dom-feat-col  { width: 2fr; padding-right: 1; }
#dom-start-col { width: 1fr; padding-right: 1; }
#dom-end-col   { width: 1fr; }
#dom-primer-results {
    height: auto; max-height: 14;
    border: solid $primary-darken-2; padding: 0 1; margin-top: 1;
    overflow-y: auto;
}
#dom-btns   { height: 3; margin-top: 1; }
#dom-btns Button { margin-right: 1; }

/* ── Primer design screen (full-screen) ─────────────────── */
#pd-box {
    width: 100%; height: 1fr;
    background: $surface; padding: 0 2;
}
#pd-title     { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#pd-source-row { height: 5; }
#pd-feat-col  { width: 2fr; padding-right: 1; }
#pd-start-col { width: 1fr; padding-right: 1; }
#pd-end-col   { width: 1fr; padding-right: 1; }
#pd-name-col  { width: 2fr; }
#pd-det-hdr   { height: 1; margin-top: 1; }
#pd-det-row   { height: 3; align: left middle; }
#pd-det-row Label { width: auto; padding: 0 0; }
#pd-det-row Input { width: 6; }
#pd-det-row Button { margin-left: 1; min-width: 20; }
#pd-clo-hdr   { height: 1; margin-top: 1; }
#pd-clo-row   { height: 5; }
#pd-clo-5col  { width: 1fr; padding-right: 1; }
#pd-clo-3col  { width: 1fr; padding-right: 1; }
#pd-clo-tmcol { width: auto; padding-right: 1; }
#pd-clo-row Button { margin-top: 2; min-width: 20; }
#pd-results   {
    height: auto; max-height: 10;
    border: solid $primary-darken-2; padding: 0 1; margin-top: 1;
}
#pd-result-names { height: 4; margin-top: 1; }
#pd-fn-col    { width: 1fr; padding-right: 1; }
#pd-rn-col    { width: 1fr; }
#pd-btns      { height: 3; margin-top: 1; }
#pd-btns Button { margin-right: 1; }
#pd-lib-hdr   { background: $accent-darken-2; color: $text; padding: 0 1; margin-top: 1; }
#pd-lib-table { height: 1fr; min-height: 4; }
"""

    BINDINGS = [
        Binding("f",           "fetch",            "Fetch GenBank", show=True),
        Binding("o",           "open_file",        "Open .gb file", show=True),
        Binding("a",           "add_to_library",   "Add to lib",    show=True),
        Binding("A",           "annotate_plasmid", "Annotate",      show=True,  key_display="A"),
        Binding("E",           "edit_seq",         "Edit seq",      show=True,  key_display="E"),
        Binding("S",           "save",             "Save",          show=True,  key_display="S"),
        Binding("[",           "rotate_cw",        "Rotate ←",      show=True,  priority=True),
        Binding("]",           "rotate_ccw",       "Rotate →",      show=True,  priority=True),
        Binding("shift+[",     "rotate_cw_lg",     "Rotate ←←",     show=False, priority=True),
        Binding("shift+]",     "rotate_ccw_lg",    "Rotate →→",     show=False, priority=True),
        Binding("home",        "reset_origin",     "Reset origin",  show=True,  priority=True),
        Binding("v",           "toggle_map_view",  "⊙/─ View",      show=True,  priority=True),
        Binding("l",           "toggle_connectors","Connectors",    show=True,  priority=True),
        Binding("r",           "toggle_restr",     "RE sites",      show=True,  priority=True),
        Binding("delete",      "delete_feature",   "Del feature",   show=True,  priority=True),
        Binding("q",           "quit",             "Quit",          show=True),
        Binding("ctrl+c",      "copy_selection",   "",              show=False, priority=True),
    ]

    # Actions that remain available even when a screen is pushed on top.
    # Everything else is suppressed to prevent confusing cross-screen
    # side-effects (e.g. pressing 'f' in the Primer Design screen should
    # not open a GenBank fetch modal underneath it).
    _ALWAYS_ALLOWED_ACTIONS: set[str] = {"quit"}

    def check_action(self, action: str, parameters: tuple) -> bool | None:
        """Textual calls this before executing any App-level binding action.

        When a Screen or ModalScreen is pushed on top of the default screen,
        we suppress all app-level actions EXCEPT quit. Screen-level bindings
        (defined on the pushed Screen itself) are unaffected — they go
        through the Screen's own check_action, not this one.
        """
        if action in self._ALWAYS_ALLOWED_ACTIONS:
            return True
        # If we're NOT on the default (main) screen, block the action.
        if self.screen.id != "_default":
            return False
        return True

    def compose(self) -> ComposeResult:
        yield Header()
        yield MenuBar()
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
                "   f fetch   o open   a add-to-lib   A annotate   E edit seq   S save   , / . circle",
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

    def action_toggle_map_view(self):
        self.query_one("#plasmid-map", PlasmidMap).action_toggle_map_view()

    def action_toggle_connectors(self):
        sp = self.query_one("#seq-panel", SequencePanel)
        pm = self.query_one("#plasmid-map", PlasmidMap)
        sp._show_connectors  = not sp._show_connectors
        pm._show_connectors  = sp._show_connectors
        sp._view_cache_key   = None   # invalidate seq panel cache
        pm._draw_cache       = None   # invalidate map draw cache
        sp._refresh_view()
        pm.refresh()
        state = "on" if sp._show_connectors else "off"
        self.notify(f"Label connectors {state}")

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

        self._restr_cache = _scan_restriction_sites(
            new_seq, min_recognition_len=self._restr_min_len, unique_only=self._restr_unique_only
        )
        displayed = self._restr_cache if self._show_restr else []
        if self._current_record is not None:
            new_record = self._rebuild_record_with_edit(new_seq, mode, s, e, new_bases)
            self._current_record = new_record
            pm.load_record(new_record)
            self.query_one("#sidebar", FeatureSidebar).populate(pm._feats)
            pm._restr_feats = displayed
            pm.refresh()
            sp.update_seq(new_seq, pm._feats + displayed)
            self.notify(f"Sequence updated  ({len(new_seq):,} bp)")
        else:
            pm._restr_feats = displayed
            pm.refresh()
            sp.update_seq(new_seq, pm._feats + displayed)

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

    def _rebuild_record_without_feature(self, feat_idx: int):
        """Create a new SeqRecord with the feat_idx-th non-source feature removed."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        new_record = SeqRecord(
            Seq(str(self._current_record.seq)),
            id=self._current_record.id,
            name=self._current_record.name,
            description=self._current_record.description,
            annotations=dict(self._current_record.annotations),
        )
        non_source_idx = 0
        for feat in self._current_record.features:
            if feat.type == "source":
                new_record.features.append(SeqFeature(
                    feat.location, type=feat.type,
                    qualifiers=dict(feat.qualifiers),
                ))
                continue
            if non_source_idx != feat_idx:
                new_record.features.append(SeqFeature(
                    feat.location, type=feat.type,
                    qualifiers=dict(feat.qualifiers),
                ))
            non_source_idx += 1
        return new_record

    def _focus_is_in_library(self) -> bool:
        """True when keyboard focus sits on the LibraryPanel or a descendant
        (e.g. the library DataTable). Used by the Delete key handler to decide
        whether the user meant to delete a feature or a library entry."""
        node = self.focused
        while node is not None:
            if isinstance(node, LibraryPanel):
                return True
            node = getattr(node, "parent", None)
        return False

    def action_delete_feature(self):
        # Focus-aware routing: when the keyboard focus is inside the library,
        # Delete targets the highlighted library row (with a confirmation
        # dialog defaulting to No), NOT the currently-selected feature. This
        # prevents a handslip after tabbing into the library from silently
        # deleting a feature the user forgot they had selected.
        if self._focus_is_in_library():
            self._request_library_delete()
            return
        pm = self.query_one("#plasmid-map", PlasmidMap)
        if pm.selected_idx < 0 or not pm._feats:
            self.notify("No feature selected.", severity="warning")
            return
        feat = pm._feats[pm.selected_idx]
        label = feat.get("label") or feat.get("type", "feature")
        self._push_undo()
        new_record = self._rebuild_record_without_feature(pm.selected_idx)
        sp = self.query_one("#seq-panel", SequencePanel)
        self._apply_snapshot(str(new_record.seq), sp._cursor_pos, new_record)
        self.notify(f"Deleted '{label}'  (Ctrl+Z to undo)")

    def _request_library_delete(self) -> None:
        """Called from action_delete_feature when the library has focus.
        Pops LibraryDeleteConfirmModal on the currently-highlighted row;
        on confirmation, removes the entry from plasmid_library.json."""
        lib = self.query_one("#library", LibraryPanel)
        t = lib.query_one("#lib-table", DataTable)
        if t.row_count == 0:
            self.notify("Library is empty.", severity="warning")
            return
        if not (0 <= t.cursor_row < t.row_count):
            self.notify("No library row highlighted.", severity="warning")
            return
        row_keys = list(t.rows.keys())
        if t.cursor_row >= len(row_keys):
            return
        entry_id = row_keys[t.cursor_row].value

        # Look up name and size for a user-friendly dialog message
        name, size = "?", 0
        for e in _load_library():
            if e.get("id") == entry_id:
                name = e.get("name") or e.get("id") or "?"
                size = e.get("size", 0)
                break

        def _on_confirm(result: "bool | None") -> None:
            if result is not True:
                return
            entries = [e for e in _load_library() if e.get("id") != entry_id]
            _save_library(entries)
            lib._repopulate()
            # If we just deleted the currently-loaded record's entry, mark
            # the header so the user knows they no longer have an active
            # library binding for it.
            if (self._current_record is not None
                    and self._current_record.id == entry_id):
                lib.set_active(None)
            self.notify(f"Removed '{name}' from library.")

        self.push_screen(
            LibraryDeleteConfirmModal(name, size, entry_id),
            callback=_on_confirm,
        )

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
        # Validate all user-data files before anything else. Corrupt files
        # are auto-restored from .bak if possible; the user is notified
        # either way so they know the state of their data.
        self._check_data_files()
        if self._preload_record is not None:
            def _load_preload():
                self._import_and_persist(self._preload_record)
            self.call_after_refresh(_load_preload)
        elif not _load_library():
            self._seed_default_library()

    def _check_data_files(self) -> None:
        """Validate plasmid library, parts bin, and primer library on startup.

        For each file:
          - Missing → first run, no warning.
          - Valid JSON array → all good.
          - Corrupt → _safe_load_json attempts .bak restore and returns a
            warning message that we surface to the user via notify().
          - Manually deleted mid-session → next load returns [], no crash.

        This runs BEFORE any load call so the caches are cold and
        _safe_load_json actually reads the files.
        """
        global _library_cache, _parts_bin_cache, _primers_cache
        for path, label, cache_attr in [
            (_LIBRARY_FILE,    "Plasmid library", "_library_cache"),
            (_PARTS_BIN_FILE,  "Parts bin",       "_parts_bin_cache"),
            (_PRIMERS_FILE,    "Primer library",  "_primers_cache"),
        ]:
            # Force a cold read (bypass cache) so we actually check the file
            globals()[cache_attr] = None
            entries, warning = _safe_load_json(path, label)
            globals()[cache_attr] = entries
            if warning:
                self.notify(warning, severity="warning", timeout=12)

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
            # First-run seed is best-effort: silent in UI, logged for debugging.
            _log.exception("Default library seed (MW463917.1) failed")

    # ── Keyboard: cursor movement, copy, undo/redo ─────────────────────────────

    def action_copy_selection(self) -> None:
        sp  = self.query_one("#seq-panel", SequencePanel)
        seq = sp._seq
        if not seq:
            return
        sel = sp._user_sel or sp._sel_range
        if sel:
            text = seq[sel[0]:sel[1]].upper()
            try:
                self.copy_to_clipboard(text)
                self.notify(f"Copied {len(text)} bp to clipboard")
            except Exception:
                if _copy_to_clipboard_osc52(text):
                    self.notify(f"Copied {len(text)} bp to clipboard")
                else:
                    self.notify("Clipboard unavailable", severity="warning")
        else:
            self.notify("No selection — click a feature or drag to select",
                        severity="information")

    def on_key(self, event) -> None:
        sp = self.query_one("#seq-panel", SequencePanel)

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

        # ── Arrow keys: move cursor; Shift+arrow extends selection ───────────
        if sp._cursor_pos < 0 or not sp._seq:
            return
        n  = len(sp._seq)
        k  = event.key
        lw = sp._line_width()
        if k in ("left", "shift+left"):
            new_pos = max(0, sp._cursor_pos - 1)
        elif k in ("right", "shift+right"):
            new_pos = min(n - 1, sp._cursor_pos + 1)
        elif k in ("up", "shift+up"):
            new_pos = max(0, sp._cursor_pos - lw)
        elif k in ("down", "shift+down"):
            new_pos = min(n - 1, sp._cursor_pos + lw)
        else:
            return
        event.stop()
        if k.startswith("shift+"):
            # Extend or shrink selection; anchor is set on first Shift+arrow press
            if sp._sel_anchor < 0:
                sp._sel_anchor = sp._cursor_pos
            sp._cursor_pos = new_pos
            s = min(sp._sel_anchor, new_pos)
            e = max(sp._sel_anchor, new_pos) + 1
            sp._user_sel  = (s, e)
            sp._sel_range = None
        else:
            sp._cursor_pos = new_pos
            sp._user_sel   = None
            sp._sel_anchor = -1
        sp._refresh_view()
        sp._ensure_cursor_visible()

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
            self._restr_cache = _scan_restriction_sites(
                seq, min_recognition_len=self._restr_min_len, unique_only=self._restr_unique_only
            )
            displayed = self._restr_cache if self._show_restr else []
            pm._restr_feats = displayed
            pm.refresh()
            sp.update_seq(seq, pm._feats + displayed)
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
        """Save current record to its source file and/or library. Returns True on success."""
        if self._current_record is None:
            self.notify("Nothing to save.", severity="warning")
            return False

        # Write to source file if one is known
        if self._source_path:
            try:
                Path(self._source_path).write_text(_record_to_gb_text(self._current_record))
            except Exception as exc:
                self.notify(f"Save failed: {exc}", severity="error")
                return False

        # Always update the library entry (add or overwrite)
        try:
            lib = self.query_one("#library", LibraryPanel)
            lib.add_entry(self._current_record)
        except Exception as exc:
            self.notify(f"Library update failed: {exc}", severity="error")
            return False

        self._mark_clean()
        if self._source_path:
            self.notify(f"Saved → {self._source_path}")
        else:
            self.notify(f"Saved {self._current_record.name} to library")
        return True

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
        # Auto-persist: fetched records go straight into the library so the
        # user never has to remember to press `a`. See _import_and_persist.
        self.push_screen(FetchModal(), callback=self._import_and_persist)

    def action_open_file(self):
        # Same auto-persist policy as fetch; _import_and_persist preserves
        # the record's _tui_source path for later "Save" operations.
        self.push_screen(OpenFileModal(), callback=self._import_and_persist)

    # ── Central record loader ──────────────────────────────────────────────────

    def _import_and_persist(self, record) -> None:
        """Apply a freshly-imported record to the UI AND save it to the library.

        Used by the three "user imported a plasmid" entry points — NCBI fetch,
        open local file, and CLI preload. Library loads, pLannotate merges,
        and undo/redo go through `_apply_record` directly so they don't
        re-save the same record.

        `add_entry` dedupes by `record.id`, so re-importing an existing entry
        updates it in place rather than creating a duplicate.
        """
        if record is None:
            return
        # _apply_record clears self._source_path; preserve the file path if
        # the record came from a local .gb file (it's stashed on the record
        # by load_genbank callers).
        source_path = getattr(record, "_tui_source", None)
        self._apply_record(record)
        if source_path is not None:
            self._source_path = source_path
        try:
            lib = self.query_one("#library", LibraryPanel)
            lib.add_entry(record)
            self.notify(f"Saved {record.name} to library.", timeout=4)
        except Exception:
            # UI already loaded the record; log and warn but don't hide it.
            _log.exception("auto-persist on import failed")
            self.notify(
                "Loaded record but could not save to library (see log).",
                severity="warning",
            )

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
        self._restr_cache = _scan_restriction_sites(
            seq_str, min_recognition_len=self._restr_min_len, unique_only=self._restr_unique_only
        )
        displayed = self._restr_cache if self._show_restr else []

        # Store restriction sites on the map for visual overlay
        pm._restr_feats = displayed
        pm.refresh()

        # Sequence panel: feature coloring = record feats + restriction sites
        seq_pnl.update_seq(seq_str, pm._feats + displayed)

        try:
            self.query_one("#library", LibraryPanel).set_active(record.id)
        except Exception:
            pass
        self._mark_clean()
        self.notify(
            f"Loaded {record.name}  ({len(record.seq):,} bp, "
            f"{len(pm._feats)} features, {len(self._restr_cache)} restriction sites)"
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
            # Single or double click: select full feature range (copyable highlight)
            seq_pnl.select_feature_range(f, cursor_bp=bp)

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

    @on(LibraryPanel.GainedFocus)
    def _library_gained_focus(self, _event) -> None:
        """Library panel (or its DataTable) just gained focus — clear any
        currently-selected feature so the Delete key can't accidentally hit
        one. The feature stays in the record; the user just re-clicks it in
        the map or sidebar to re-select."""
        try:
            pm = self.query_one("#plasmid-map", PlasmidMap)
        except Exception:
            return
        if pm.selected_idx < 0:
            return   # nothing to clear
        pm.selected_idx = -1
        pm.refresh()
        try:
            self.query_one("#sidebar", FeatureSidebar).show_detail(None)
        except Exception:
            pass

    @on(LibraryPanel.RenameRequested)
    def _library_rename_requested(self, event: LibraryPanel.RenameRequested):
        """Library's ✎ button was clicked. Opens RenamePlasmidModal; on Save,
        updates the library JSON (both the `name` field and the stored
        gb_text's LOCUS line) and — if the renamed entry is currently loaded
        — mutates `self._current_record.name` in place so the circular map,
        title bar, and next save all show the new name immediately."""
        if event.entry_id is None:
            self.notify("Highlight a library row first.", severity="warning")
            return
        current_name: "str | None" = None
        for e in _load_library():
            if e.get("id") == event.entry_id:
                current_name = e.get("name") or e.get("id")
                break
        if current_name is None:
            self.notify("Library entry not found.", severity="warning")
            return

        def _on_result(new_name: "str | None") -> None:
            if new_name is None:
                return   # user cancelled or no-op rename
            # Collision check: refuse if another entry already has that name
            # (different id). Same-entry is already filtered by the modal.
            for e in _load_library():
                if (e.get("id") != event.entry_id
                        and e.get("name") == new_name):
                    self.notify(
                        f"A plasmid named {new_name!r} already exists.",
                        severity="error",
                    )
                    return
            self._rename_library_entry(event.entry_id, new_name)

        self.push_screen(
            RenamePlasmidModal(current_name, event.entry_id),
            callback=_on_result,
        )

    def _rename_library_entry(self, entry_id: str, new_name: str) -> None:
        """Persist a rename:
          1. Update the `name` field in the library JSON for this entry
          2. Re-serialize the stored gb_text so the LOCUS line carries the
             new name — otherwise re-loading the library would show the old
             name from the GenBank record's internal header
          3. If this entry is the currently-loaded record, mutate
             `self._current_record.name` in place, invalidate the
             PlasmidMap render cache, refresh the map, and update the
             window title bar
        """
        entries = _load_library()
        for e in entries:
            if e.get("id") == entry_id:
                e["name"] = new_name
                # Re-serialize the stored gb_text with the new LOCUS name.
                # If the gb_text can't be parsed for any reason, fall back
                # to just updating the JSON `name` field — the library row
                # will show the new name and the next explicit save will
                # fix the gb_text.
                try:
                    rec = _gb_text_to_record(e.get("gb_text", ""))
                    rec.name = new_name
                    rec.id   = entry_id   # don't let SeqIO rewrite the id
                    e["gb_text"] = _record_to_gb_text(rec)
                except Exception:
                    _log.exception(
                        "rename: could not re-serialize gb_text for %s",
                        entry_id,
                    )
                break
        else:
            self.notify("Library entry vanished.", severity="warning")
            return
        _save_library(entries)

        # Refresh the library table so the new name shows.
        lib = self.query_one("#library", LibraryPanel)
        lib._repopulate()

        # If this is the currently-loaded record, update the in-memory object
        # so the map, title bar, and every other bit of UI tracking `record.
        # name` pick up the new name without a full reload.
        if (self._current_record is not None
                and self._current_record.id == entry_id):
            self._current_record.name = new_name
            try:
                pm = self.query_one("#plasmid-map", PlasmidMap)
                # record.name is in the cache key, but nuke it explicitly
                # for belt-and-braces (in case future refactors drop the
                # name from the key).
                pm._draw_cache = None
                pm.refresh()
            except Exception:
                pass
            # Refresh the window title via _mark_clean (which rebuilds it
            # from self._current_record.name).
            self._mark_clean()

        self.notify(f"Renamed to {new_name}.")

    @on(LibraryPanel.AnnotateRequested)
    def _library_annotate_requested(self, event: LibraryPanel.AnnotateRequested):
        """Library's ◈ button was clicked. If the focused row isn't the
        currently-loaded record, load it first, then run annotation."""
        if event.entry_id is None:
            # No row focused — fall back to annotating the current record.
            self.action_annotate_plasmid()
            return
        need_load = (
            self._current_record is None
            or self._current_record.id != event.entry_id
        )
        if need_load:
            for entry in _load_library():
                if entry.get("id") == event.entry_id:
                    try:
                        record = _gb_text_to_record(entry.get("gb_text", ""))
                        self._apply_record(record)
                    except Exception as exc:
                        _log.exception("library load for annotation failed")
                        self.notify(
                            f"Failed to load library entry: {exc}",
                            severity="error",
                        )
                        return
                    break
            else:
                self.notify("Library entry not found.", severity="warning")
                return
        self.action_annotate_plasmid()

    # ── Menu bar ───────────────────────────────────────────────────────────────

    def _rescan_restrictions(self) -> None:
        """Re-scan restriction sites with current settings and update UI."""
        sp = self.query_one("#seq-panel", SequencePanel)
        pm = self.query_one("#plasmid-map", PlasmidMap)
        if not sp._seq:
            return
        self._restr_cache = _scan_restriction_sites(
            sp._seq,
            min_recognition_len=self._restr_min_len,
            unique_only=self._restr_unique_only,
        )
        displayed = self._restr_cache if self._show_restr else []
        pm._restr_feats = displayed
        pm.refresh()
        sp.update_seq(sp._seq, pm._feats + displayed)

    def _apply_restr_visibility(self) -> None:
        """Push current cache to map/sequence panel respecting _show_restr flag."""
        sp = self.query_one("#seq-panel", SequencePanel)
        pm = self.query_one("#plasmid-map", PlasmidMap)
        if not sp._seq:
            return
        # Rescan if cache is stale or empty (e.g. loaded before cache was wired up)
        if self._show_restr and not self._restr_cache:
            self._restr_cache = _scan_restriction_sites(
                sp._seq,
                min_recognition_len=self._restr_min_len,
                unique_only=self._restr_unique_only,
            )
        displayed = self._restr_cache if self._show_restr else []
        pm._restr_feats = displayed
        pm.refresh()
        sp.update_seq(sp._seq, pm._feats + displayed)

    def open_menu(self, name: str, x: int, y: int) -> None:
        """Open a menu. Single-action menus skip the dropdown and fire their
        action directly; multi-action menus show the dropdown list."""

        # ── Direct-action menus (no dropdown) ──────────────────────────────
        # These open their target screen / notification instantly, avoiding
        # the jarring dark-overlay intermediate.
        if name == "Parts":
            self.action_open_parts_bin()
            return
        if name == "Constructor":
            self.action_open_constructor()
            return
        if name == "Primers":
            self.action_open_primer_design()
            return

        # ── Multi-action menus (dropdown) ──────────────────────────────────
        ck = "\u2713"  # checkmark
        nc = " "
        u  = ck if self._restr_unique_only else nc
        m6 = ck if self._restr_min_len == 6  else nc
        m4 = ck if self._restr_min_len == 4  else nc
        rs = ck if self._show_restr        else nc

        menus = {
            "File": [
                ("Open .gb file",   "open_file"),
                ("Fetch from NCBI", "fetch"),
                ("---",             None),
                ("Add to Library",  "add_to_library"),
                ("Save",            "save"),
                ("---",             None),
                ("Quit",            "quit"),
            ],
            "Edit": [
                ("Edit Sequence",   "edit_seq"),
                ("---",             None),
                ("Undo",            "undo"),
                ("Redo",            "redo"),
                ("---",             None),
                ("Delete Feature",  "delete_feature"),
            ],
            "Enzymes": [
                (f"[{rs}] Show RE sites  [r]",   "toggle_restr"),
                ("---",                            None),
                (f"[{u}] Unique cutters",         "toggle_restr_unique"),
                (f"[{m6}] 6+ bp sites",           "toggle_restr_min6"),
                (f"[{m4}] 4+ bp sites",           "toggle_restr_min4"),
                ("---",                            None),
                ("Toggle connectors",              "toggle_connectors"),
            ],
            "Features": [
                ("Add Feature...",              "add_feature"),
                ("Delete Feature",              "delete_feature"),
                ("---",                         None),
                ("Annotate with pLannotate",    "annotate_plasmid"),
                ("---",                         None),
                ("Toggle connectors",           "toggle_connectors"),
            ],
        }
        items = menus.get(name, [])
        if not items:
            return
        self.push_screen(
            DropdownScreen(items, x, y),
            callback=self._menu_action,
        )

    def _menu_action(self, action: "str | None") -> None:
        if action is None:
            return
        # Handle toggle actions directly since they need state updates
        if action in ("toggle_restr", "toggle_restr_unique", "toggle_restr_min6", "toggle_restr_min4"):
            getattr(self, f"action_{action}")()
        else:
            getattr(self, f"action_{action}")()

    def action_toggle_restr(self) -> None:
        self._show_restr = not self._show_restr
        self._apply_restr_visibility()
        state = "shown" if self._show_restr else "hidden"
        self.notify(f"Restriction enzymes {state}")

    def action_toggle_restr_unique(self) -> None:
        self._restr_unique_only = not self._restr_unique_only
        self._rescan_restrictions()
        state = "on" if self._restr_unique_only else "off"
        self.notify(f"Unique cutters {state}")

    def action_toggle_restr_min6(self) -> None:
        self._restr_min_len = 6
        self._rescan_restrictions()
        self.notify("Showing 6+ bp recognition sites")

    def action_toggle_restr_min4(self) -> None:
        self._restr_min_len = 4
        self._rescan_restrictions()
        self.notify("Showing 4+ bp recognition sites")

    def action_add_feature(self) -> None:
        self.notify("Add feature: coming soon", severity="information")

    # ── pLannotate annotation ──────────────────────────────────────────────────

    def action_annotate_plasmid(self) -> None:
        """Run pLannotate on the currently-loaded record (shortcut: Shift+A)."""
        if self._current_record is None:
            self.notify(
                "Load a plasmid first (press 'f' to fetch or 'o' to open).",
                severity="warning",
            )
            return
        status = _plannotate_status()
        if not status["ready"]:
            # Specific, actionable error. Detect which piece is missing so the
            # user knows what to install.
            if not status["installed"]:
                self.notify(
                    "pLannotate not installed. " + _plannotate_install_hint(),
                    severity="error", timeout=10,
                )
            else:
                missing = [k for k in ("blast", "diamond") if not status[k]]
                self.notify(
                    f"pLannotate needs {' + '.join(missing)} on PATH. "
                    + _plannotate_install_hint(),
                    severity="error", timeout=10,
                )
            return
        # Preflight the size cap so the user gets the error instantly instead
        # of waiting for pLannotate to reject it.
        n = len(self._current_record.seq)
        if n > _PLANNOTATE_MAX_BP:
            self.notify(
                f"pLannotate caps inputs at {_PLANNOTATE_MAX_BP:,} bp "
                f"(this plasmid: {n:,} bp).",
                severity="warning",
            )
            return
        self.notify(
            f"Running pLannotate on {self._current_record.name} "
            f"({n:,} bp)… this takes 5-30 s.",
            timeout=15,
        )
        self._run_plannotate_worker(self._current_record)

    @work(thread=True)
    def _run_plannotate_worker(self, record) -> None:
        """Background worker: runs pLannotate subprocess, merges, applies.
        Errors are logged and surfaced to the UI via notify(); nothing raw
        reaches the user."""
        try:
            annotated = _run_plannotate(record)
            merged    = _merge_plannotate_features(record, annotated)
        except PlannotateError as exc:
            _log.info("pLannotate: %s", exc)
            def _notify_err():
                self.notify(exc.user_msg, severity="error", timeout=10)
            self.call_from_thread(_notify_err)
            return
        except Exception as exc:
            _log.exception("pLannotate worker crashed")
            def _notify_crash():
                self.notify(
                    f"pLannotate crashed: {exc}", severity="error", timeout=10,
                )
            self.call_from_thread(_notify_crash)
            return

        n_added = getattr(merged, "_plannotate_added", 0)
        def _apply():
            if n_added == 0:
                self.notify(
                    "pLannotate found no new features (all hits duplicated "
                    "existing annotations).",
                    severity="information",
                )
                return
            self._push_undo()          # annotation is undo-able
            self._apply_record(merged)
            self.query_one("#library", LibraryPanel).set_dirty(True)
            self.notify(
                f"Added {n_added} pLannotate feature"
                f"{'s' if n_added != 1 else ''}. "
                "Press 'a' to save to library.",
                timeout=6,
            )
        self.call_from_thread(_apply)

    def action_open_parts_bin(self) -> None:
        self.push_screen(PartsBinModal())

    def action_open_constructor(self) -> None:
        self.push_screen(ConstructorModal())

    def action_open_primer_design(self) -> None:
        """Open the full-screen Primer Design workbench. Passes the current
        plasmid's sequence and features so the user can select regions."""
        rec = self._current_record
        seq   = str(rec.seq) if rec else ""
        name  = rec.name if rec else ""
        feats = []
        try:
            feats = self.query_one("#plasmid-map", PlasmidMap)._feats
        except Exception:
            pass
        self.push_screen(PrimerDesignScreen(seq, feats, name))

    def action_undo(self) -> None:
        self._action_undo()

    def action_redo(self) -> None:
        self._action_redo()

    # ── Sequence edits ─────────────────────────────────────────────────────────

    @on(SequencePanel.SequenceChanged)
    def _seq_changed(self, event: SequencePanel.SequenceChanged):
        # Update restriction site overlay whenever sequence changes
        pm = self.query_one("#plasmid-map", PlasmidMap)
        self._restr_cache = _scan_restriction_sites(
            event.seq,
            min_recognition_len=self._restr_min_len,
            unique_only=self._restr_unique_only,
        )
        displayed = self._restr_cache if self._show_restr else []
        pm._restr_feats = displayed
        pm.refresh()


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    _log_startup_banner()
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    app = PlasmidApp()

    if arg:
        if Path(arg).exists():
            try:
                record = load_genbank(arg)
                record._tui_source = str(Path(arg).resolve())
            except Exception as exc:
                _log.exception("Failed to load %s", arg)
                print(f"Could not load {arg!r}: {exc}", file=sys.stderr)
                sys.exit(1)
        else:
            print(f"Fetching {arg!r} from NCBI…", flush=True)
            try:
                record = fetch_genbank(arg)
                print(f"  Got: {record.name}  ({len(record.seq)} bp)")
            except Exception as exc:
                _log.exception("NCBI fetch failed for %s", arg)
                print(f"Fetch failed: {exc}", file=sys.stderr)
                sys.exit(1)
        app._preload_record = record

    try:
        app.run()
    except Exception:
        _log.exception("App terminated with unhandled exception")
        raise
    finally:
        _log.info("SpliceCraft session %s ending", _SESSION_ID)


if __name__ == "__main__":
    main()
