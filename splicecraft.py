#!/usr/bin/env python3
"""
splicecraft.py
==============
SpliceCraft — terminal circular plasmid map viewer.

Features:
  - Fetch any GenBank record by accession (pUC19 = L09137)
  - Load local .gb / .gbk (GenBank) or .dna (SnapGene) files
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

__version__ = "0.3.1"

# ── User data directory ────────────────────────────────────────────────────────
# All user-writable state (library, parts bin, primers, .bak files) lives in
# the platform-appropriate data dir:
#   Linux:   ~/.local/share/splicecraft/
#   macOS:   ~/Library/Application Support/splicecraft/
#   Windows: %APPDATA%\splicecraft\
# Override with $SPLICECRAFT_DATA_DIR (useful for tests and portable installs).

def _user_data_dir() -> Path:
    override = os.environ.get("SPLICECRAFT_DATA_DIR")
    if override:
        p = Path(override).expanduser()
    else:
        try:
            from platformdirs import user_data_dir
            p = Path(user_data_dir("splicecraft", appauthor=False, roaming=False))
        except ImportError:
            p = Path.home() / ".local" / "share" / "splicecraft"
    p.mkdir(parents=True, exist_ok=True)
    return p

_DATA_DIR = _user_data_dir()


def _migrate_legacy_data() -> None:
    """One-shot migration from Path(__file__).parent → _DATA_DIR.
    Idempotent: only copies files whose destination doesn't already exist.
    Preserves the source files (copy, not move) so a dev running from the
    repo checkout can still use them via $SPLICECRAFT_DATA_DIR."""
    import shutil
    legacy_root = Path(__file__).parent
    if legacy_root.resolve() == _DATA_DIR.resolve():
        return   # running from the data dir itself (no migration needed)
    names = [
        "plasmid_library.json", "parts_bin.json", "primers.json",
        "plasmid_library.json.bak", "parts_bin.json.bak", "primers.json.bak",
    ]
    migrated = []
    for name in names:
        src = legacy_root / name
        dst = _DATA_DIR / name
        if src.exists() and not dst.exists():
            try:
                shutil.copy2(src, dst)
                migrated.append(name)
            except OSError:
                pass
    if migrated:
        try:
            (_DATA_DIR / ".migrated").write_text("\n".join(migrated))
        except OSError:
            pass

_migrate_legacy_data()


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
# Rotating file log with an 8-char session ID prefix on every line so multi-run
# logs are greppable. Default path is `_DATA_DIR/logs/splicecraft.log` so logs
# survive reboots on systemd-tmpfiles distros (which wipe /tmp on boot).
# Overridable via $SPLICECRAFT_LOG. UI never sees raw tracebacks — they go here.

def _default_log_path() -> str:
    override = os.environ.get("SPLICECRAFT_LOG")
    if override:
        return override
    try:
        log_dir = _DATA_DIR / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        return str(log_dir / "splicecraft.log")
    except OSError:
        # Fall back to /tmp if the data dir is somehow unwritable (read-only
        # home, quota, etc.) — better than crashing at import time.
        return "/tmp/splicecraft.log"

_LOG_PATH   = _default_log_path()
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
        # Last-resort no-op handler if the log dir is read-only.
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
from textual.css.query import NoMatches
from textual.events import Click, MouseDown, MouseMove, MouseUp, MouseScrollDown, MouseScrollUp
from textual.message import Message
from textual.reactive import reactive
from textual.screen import ModalScreen, Screen
from textual.widget import Widget
from textual.widgets import (
    Button, DataTable, Footer, Header, Input, Label, ListItem, ListView,
    RadioButton, RadioSet, Select, Static, TextArea,
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
# All user data (plasmid library, parts bin, primer library, codon tables)
# goes through _safe_save_json which:
#   1. Backs up the existing file to *.bak BEFORE overwriting
#   2. Writes via tempfile + os.replace (atomic on POSIX — the file is either
#      fully written or not at all; no partial-write corruption)
#   3. Logs every write with entry count for post-mortem debugging
#
# _safe_load_json handles the read side:
#   - Missing file → [] (first run, not an error)
#   - Corrupt file → attempt restore from .bak; if .bak also corrupt → []
#   - Returns (entries, warning_message_or_None)
#
# ── On-disk format ────────────────────────────────────────────────────────────
# Current format (schema version 1):
#
#     {"_schema_version": 1, "entries": [...]}
#
# Legacy format (pre-0.3.1) was a bare JSON list. Loaders accept both so users
# upgrading in-place don't need to regenerate their libraries; the legacy file
# gets silently rewritten as an envelope on the next save.

_CURRENT_SCHEMA_VERSION = 1


def _extract_entries(raw, label: str) -> "tuple[list | None, str | None]":
    """Return (entries, warning) from a parsed-JSON payload.

    Accepts both the envelope format `{"_schema_version": N, "entries": [...]}`
    and the legacy bare-list format. Returns (None, warning) on unknown shape
    so the caller can fall through to the .bak.
    """
    if isinstance(raw, list):
        return raw, None
    if isinstance(raw, dict) and isinstance(raw.get("entries"), list):
        version = raw.get("_schema_version")
        if version is not None and version > _CURRENT_SCHEMA_VERSION:
            # Written by a newer SpliceCraft. Load the entries but warn so
            # the user knows fields may be silently dropped on re-save.
            return list(raw["entries"]), (
                f"{label} was written by a newer SpliceCraft "
                f"(schema v{version} > v{_CURRENT_SCHEMA_VERSION}) — some "
                f"fields may be lost on save."
            )
        return list(raw["entries"]), None
    return None, f"{label}: unexpected JSON shape ({type(raw).__name__})"


def _safe_save_json(path: Path, entries: list, label: str,
                    schema_version: int = _CURRENT_SCHEMA_VERSION) -> None:
    """Atomically write `entries` as JSON to `path`, backing up first.

    Writes the envelope format `{"_schema_version": N, "entries": [...]}`.

    The .bak file is the user's safety net — if a write goes wrong or the
    app crashes mid-save, the previous version survives as path.bak.

    **Shrink guard**: if the file currently has N entries and we're about to
    write M < N, we still write (the user may have legitimately deleted
    entries) BUT we log a loud warning so accidental nukes are visible in
    the session log for post-mortem debugging. The .bak always preserves
    the pre-write state regardless.
    """
    import os
    import tempfile

    # 1. Back up the existing file (if it has content)
    existing_count = 0
    if path.exists():
        try:
            existing = path.read_bytes()
            if existing.strip():
                bak = path.with_suffix(path.suffix + ".bak")
                bak.write_bytes(existing)
                # Count existing entries for the shrink guard. Accept both
                # envelope and legacy formats so upgrades don't trip the
                # shrink guard on their first save.
                try:
                    prev = json.loads(existing)
                    prev_entries, _ = _extract_entries(prev, label)
                    if prev_entries is not None:
                        existing_count = len(prev_entries)
                except Exception:
                    pass
        except OSError:
            _log.warning("Could not create backup for %s", path)

    # Shrink guard: log a loud warning if we're about to lose entries.
    # Any shrink is logged so accidental nukes are auditable; going to zero
    # from a populated file is almost always a bug (user never legitimately
    # empties the whole library at once) so the .bak is explicitly preserved
    # and the warning cites the path for recovery.
    if existing_count > 0 and len(entries) < existing_count:
        _log.warning(
            "SHRINK GUARD: %s is being overwritten with %d entries "
            "(was %d). If this is unexpected, restore from %s.bak",
            label, len(entries), existing_count, path,
        )

    payload = {"_schema_version": schema_version, "entries": entries}

    # 2. Atomic write: tempfile in same dir → os.replace
    try:
        fd, tmp_name = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
                fh.flush()
                try:
                    os.fsync(fh.fileno())
                except OSError:
                    pass
            os.replace(tmp_name, str(path))
            _log.info("Saved %s: %d entries to %s (schema v%d)",
                      label, len(entries), path, schema_version)
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
    """Load a JSON payload from `path`. Returns (entries, warning_or_None).

    Accepts both the current envelope format (`{"_schema_version": N,
    "entries": [...]}`) and the legacy flat-list format written by
    SpliceCraft < 0.3.1. The legacy file gets silently rewritten as an
    envelope on the next save.

    - Missing file → ([], None) — normal first run, no warning.
    - Valid file   → (entries, None).
    - Corrupt file → attempt .bak restore; if .bak is valid →
      (bak_entries, warning). If .bak also corrupt → ([], warning).
    """
    if not path.exists():
        return [], None

    # Try the main file
    main_warning: "str | None" = None
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
        entries, shape_warn = _extract_entries(raw, label)
        if entries is not None:
            return entries, shape_warn
        _log.warning("%s: %s", path, shape_warn)
        main_warning = shape_warn
    except Exception:
        _log.exception("Corrupt %s file: %s", label, path)

    # Main file is corrupt — try the .bak
    bak = path.with_suffix(path.suffix + ".bak")
    if bak.exists():
        try:
            raw = json.loads(bak.read_text(encoding="utf-8"))
            entries, _ = _extract_entries(raw, label)
            if entries is not None:
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

    return [], (main_warning
                or f"{label} is corrupt and no valid backup was found. "
                   "Starting empty.")


# ── Library persistence ────────────────────────────────────────────────────────

_LIBRARY_FILE = _DATA_DIR / "plasmid_library.json"
_library_cache: "list | None" = None

def _load_library() -> list[dict]:
    global _library_cache
    if _library_cache is not None:
        return list(_library_cache)
    entries, warning = _safe_load_json(_LIBRARY_FILE, "Plasmid library")
    if warning:
        _log.warning(warning)
    # Drop non-dict entries (hand-edited JSON, schema drift) so that
    # .get() calls downstream don't raise AttributeError.
    entries = [e for e in entries if isinstance(e, dict)]
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

# Case-preserving ACGT complement used by the sequence-panel renderer.
_DNA_COMP_PRESERVE_CASE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _rc(seq: str) -> str:
    return seq.upper().translate(_IUPAC_COMP)[::-1]


def _feat_len(start: int, end: int, total: int) -> int:
    """Circular-aware feature length. A wrap feature (end < start) is
    (total - start) + end bp long; a linear feature is end - start."""
    return (total - start) + end if end < start else end - start


def _slice_circular(seq: str, start: int, end: int) -> str:
    """Circular-aware slice. If end > start this is a normal slice; if
    end < start the slice wraps the origin and returns seq[start:] + seq[:end].
    end == start is treated as empty (not "wrap whole plasmid") — callers
    that want the latter should pass explicit boundaries. Used by the
    primer-design helpers so a region straddling the origin can be
    primer-designed without special casing at every call site.
    """
    if end >= start:
        return seq[start:end]
    return seq[start:] + seq[:end]


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
    circular: bool = True,
) -> list[dict]:
    """Scan both strands; return resite + recut dicts for every hit.

    resite — the recognition sequence span (colored bar)
    recut  — the cut position (single-bp marker: ↓ above or ↑ below DNA)

    min_recognition_len — skip enzymes whose recognition site is shorter than this
                          (default 6 to reduce noise from 4-cutters)
    unique_only         — if True, only include enzymes that cut exactly once
                          (forward + reverse strand combined; default True)
    circular            — if True (default), recognition sequences that span
                          the origin (bp n-1 → bp 0) are also found. SpliceCraft
                          is a plasmid viewer, so circularity is on by default.

    Wrap-around resites are emitted as TWO pieces so the existing linear-span
    rendering in the map / sequence panel stays correct: one piece on the
    "tail" (start..n) with the enzyme label, one unlabeled piece on the
    "head" (0..tail_len). The single-bp recut marker is placed at its real
    absolute position modulo n.
    """
    seq_u = seq.upper()
    n = len(seq_u)
    # Per-enzyme results collected first so we can filter to unique cutters
    by_enzyme: dict[str, list[dict]] = {}
    seen: set[tuple[str, int, int]] = set()   # deduplicate palindromes

    # For circular sequences, scan an augmented copy that includes up to
    # site_len-1 bp re-attached from the beginning so matches starting near
    # the end (that would otherwise be truncated) are found too.
    max_site_len = max((e[2] for e in _SCAN_CATALOG), default=0)
    scan_seq = (seq_u + seq_u[: max_site_len - 1]) if (circular and n > 0) else seq_u

    def _emit_resite(hits, p, site_len, strand, color, name,
                     cut_col, ext_cut_bp):
        """Emit one or two resite dicts depending on wrap. Labels only on the
        first piece so the map doesn't double-print. For wrapped sites, the
        cut_col / ext_cut_bp fields are only meaningful on the piece that
        actually contains the cut; we attach them to the tail piece by default
        and clear them on the head piece."""
        if p + site_len <= n:
            hits.append({
                "type":       "resite",
                "start":      p,
                "end":        p + site_len,
                "strand":     strand,
                "color":      color,
                "label":      name,
                "cut_col":    cut_col,
                "ext_cut_bp": ext_cut_bp,
            })
            return
        # Wraps origin: tail [p, n) + head [0, (p + site_len) - n).
        tail_len = n - p
        head_len = (p + site_len) - n
        # cut_col (bar-relative) maps to whichever piece actually contains the
        # cut. ext_cut_bp (absolute) is unrelated to the tail/head split — it's
        # only meaningful when cut_col is None (Type IIS cuts outside the
        # recognition sequence). Attach it to both pieces so the cut arrow is
        # drawn regardless of which chunk contains the external cut position;
        # the chunk-level `chunk_start <= ext_cut_bp < chunk_end` test makes
        # the render idempotent. Regression guard added 2026-04-13.
        tail_cut_col = cut_col if (cut_col is not None and cut_col < tail_len) else None
        head_cut_col = ((cut_col - tail_len) if (cut_col is not None and cut_col >= tail_len)
                        else None)
        hits.append({
            "type":       "resite",
            "start":      p,
            "end":        n,
            "strand":     strand,
            "color":      color,
            "label":      name,
            "cut_col":    tail_cut_col,
            "ext_cut_bp": ext_cut_bp,
        })
        hits.append({
            "type":       "resite",
            "start":      0,
            "end":        head_len,
            "strand":     strand,
            "color":      color,
            "label":      "",     # unlabeled continuation
            "cut_col":    head_cut_col,
            "ext_cut_bp": ext_cut_bp,
        })

    for entry in _SCAN_CATALOG:
        name, site, site_len, fwd_cut, rev_cut, color, pat, is_palindrome, rc_pat = entry
        if site_len < min_recognition_len:
            continue
        hits: list[dict] = []

        # Forward strand scan (over augmented sequence if circular)
        for m in pat.finditer(scan_seq):
            p = m.start()
            if p >= n:
                continue   # duplicate of match already found at p - n
            key = (name, p, 1)
            if key in seen:
                continue
            seen.add(key)
            # ext_cut_bp: absolute cut position when cut falls outside recognition
            _ext = ((p + fwd_cut) % n) if (fwd_cut <= 0 or fwd_cut >= site_len) else None
            _cc  = fwd_cut if 0 < fwd_cut < site_len else None
            _emit_resite(hits, p, site_len, 1, color, name, _cc, _ext)
            cut_bp = (p + fwd_cut) % n if n > 0 else 0
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
            for m in rc_pat.finditer(scan_seq):
                p = m.start()
                if p >= n:
                    continue   # duplicate of match already found at p - n
                key = (name, p, -1)
                if key in seen:
                    continue
                seen.add(key)
                # Cut column within the bar: enzyme's fwd_cut mapped to
                # the reversed orientation displayed on the forward strand
                rev_cut_col = site_len - 1 - fwd_cut
                _top_cut_bp = (p + site_len - 1 - rev_cut) % n   # top-strand cut in fwd coords
                _top_cut_outside = ((_top_cut_bp - p) % n) >= site_len
                _cc  = rev_cut_col if 0 <= rev_cut_col < site_len else None
                _ext = _top_cut_bp if _top_cut_outside else None
                _emit_resite(hits, p, site_len, -1, color, name, _cc, _ext)
                # Bottom-strand cut (enzyme's fwd_cut mapped to fwd coords)
                cut_bp = (p + site_len - 1 - fwd_cut) % n if n > 0 else 0
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
        # Count LABELED resites only — a wrap-around hit is emitted as one
        # labeled piece + one unlabeled continuation, but counts as 1 site.
        if unique_only:
            n_sites = sum(
                1 for h in hits if h["type"] == "resite" and h.get("label")
            )
            if n_sites != 1:
                continue
        # Skip isoschizomers / HF-variants that land on an already-placed site
        positions = {
            (h["start"], h["end"]) for h in hits
            if h["type"] == "resite" and h.get("label")
        }
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


def _feats_in_chunk(
    feats: list[dict], chunk_start: int, chunk_end: int, total: int
) -> list[dict]:
    """Return features overlapping [chunk_start, chunk_end), with wrap-around
    features (end < start) split into tail [start, total) + head [0, end)
    virtual pieces. The tail keeps the label; the head is unlabeled so the
    name only appears once per chunk row. Non-wrapped features pass through
    unchanged (same dict identity).

    Needed because the naive overlap test `start < chunk_end and end > chunk_start`
    drops wrap features from every chunk (both halves fail the conjunction).
    """
    out: list[dict] = []
    for f in feats:
        s, e = f["start"], f["end"]
        if e >= s:
            if s < chunk_end and e > chunk_start:
                out.append(f)
            continue
        # Wrap feature: split into tail + head.
        if s < chunk_end and total > chunk_start:
            out.append({**f, "end": total})
        if 0 < chunk_end and e > chunk_start:
            out.append({**f, "start": 0, "label": ""})
    return out


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
        if f["end"] >= f["start"]:
            for i in range(f["start"], min(f["end"], n)):
                styles[i] = col
        else:
            # Wrap feature: colour tail [start, n) + head [0, end).
            for i in range(f["start"], n):
                styles[i] = col
            for i in range(0, min(f["end"], n)):
                styles[i] = col
    annot_feats = sorted(
        [f for f in feats if f.get("type") not in ("site", "recut")],
        key=lambda f: -_feat_len(f["start"], f["end"], max(n, 1)),
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
    result    = Text(no_wrap=True, overflow="crop")

    for chunk_start in range(0, n, line_width):
        chunk_end = min(chunk_start + line_width, n)

        # ── Assign features to lane groups ──
        # _feats_in_chunk handles wrap features (end < start) by splitting
        # into tail + head virtual pieces before the overlap test.
        chunk_feats = _feats_in_chunk(annot_feats, chunk_start, chunk_end, n)
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
        chunk_rev = chunk_fwd.translate(_DNA_COMP_PRESERVE_CASE)
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

    Wrapped CDSes (end < start) are represented that way by the GenBank loader
    for `join(tail..end, 0..head)` origin-spanning features on circular
    plasmids. We concatenate the tail + head before translating so the
    protein comes out correctly. Regression guard added 2026-04-13.
    """
    if end < start:
        sub = (full_seq[start:] + full_seq[:end]).upper()
    else:
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
            if not isinstance(s, str):
                continue
            # Collapse whitespace characters (newline, tab, vertical tab)
            # into single spaces so a multi-line /note="…" qualifier
            # doesn't break the sidebar row or clobber the map label.
            # Then strip and fall through if the result is empty.
            s = " ".join(s.split())
            if s:
                return s[:28]
    return feat.type

def _nice_tick(total: int) -> int:
    """A tick interval that gives ~6-10 ticks for this plasmid size."""
    for t in [50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 25000, 50000]:
        if 4 <= total // t <= 14:
            return t
    return max(1, total // 8)

# ── GenBank I/O ────────────────────────────────────────────────────────────────

def _pick_single_record(records: list, source: str):
    """Given a list of SeqRecords, return the single one if there's exactly
    one, else raise ValueError with a user-friendly message. Used by both
    NCBI fetch and file load so the error text is consistent.
    """
    if not records:
        raise ValueError(
            f"{source} contained no GenBank records. Is it a valid .gb/.gbk file?"
        )
    if len(records) > 1:
        ids = ", ".join(r.id for r in records[:3])
        more = f" (and {len(records) - 3} more)" if len(records) > 3 else ""
        raise ValueError(
            f"{source} contains {len(records)} records — SpliceCraft loads "
            f"one plasmid at a time. Split the file or extract a single "
            f"record first (found: {ids}{more})."
        )
    return records[0]

_NCBI_TIMEOUT_S = 30   # cap long NCBI hangs; the UI worker can't otherwise cancel

def fetch_genbank(accession: str, email: str = "splicecraft@local"):
    """Fetch a GenBank record by accession from NCBI Entrez. Returns SeqRecord.

    Raises ValueError with a user-friendly message if NCBI returns no
    records (obsolete accession) or multiple records. A 30 s socket timeout
    is applied so a silent network stall surfaces as an error instead of
    pinning the worker thread forever.
    """
    import socket
    from Bio import Entrez, SeqIO
    Entrez.email = email
    prev_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(_NCBI_TIMEOUT_S)
    try:
        with Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="text"
        ) as handle:
            records = list(SeqIO.parse(handle, "genbank"))
    finally:
        socket.setdefaulttimeout(prev_timeout)
    return _pick_single_record(records, f"NCBI accession {accession!r}")

def _detect_plasmid_format(path: str) -> str:
    """Pick a Biopython SeqIO format key from a file path's extension.

    Supported:
      - GenBank        (.gb, .gbk, .genbank)       → "genbank"
      - SnapGene       (.dna)                       → "snapgene"

    Extensions are matched case-insensitively. Unknown extensions
    default to "genbank" since that's the most common plasmid format;
    the parser will then raise a clear error if the contents don't
    match.
    """
    from pathlib import Path
    suffix = Path(path).suffix.lower()
    if suffix == ".dna":
        return "snapgene"
    # .gb, .gbk, .genbank, or anything else — try GenBank.
    return "genbank"


def load_genbank(path: str):
    """Load a plasmid file (GenBank .gb/.gbk or SnapGene .dna). Returns
    SeqRecord.

    Despite the name (kept for backward compatibility), this also
    handles SnapGene native .dna files via Biopython's `snapgene`
    parser. Dispatch is based on file extension.

    For SnapGene files, Biopython leaves `record.id` / `record.name`
    as `<unknown id>` / `<unknown name>` sentinels (SnapGene's own
    name/title metadata is not exposed through SeqIO). Backfill both
    from the file stem so the library and map title show something
    human-readable instead of `<unknown name>`.

    Raises ValueError with a user-friendly message if the file has no
    records or multiple records.
    """
    from Bio import SeqIO
    from pathlib import Path as _P
    fmt = _detect_plasmid_format(path)
    try:
        records = list(SeqIO.parse(path, fmt))
    except ValueError as exc:
        # SnapGene parser raises ValueError on malformed files; rewrap
        # with a more useful message.
        if fmt == "snapgene":
            raise ValueError(
                f"Could not parse SnapGene file {path}: {exc}. "
                f"If this file was exported from an old SnapGene version, "
                f"try re-exporting as .dna from a current SnapGene release."
            ) from exc
        raise
    rec = _pick_single_record(records, path)

    # SnapGene (and occasionally minimally-annotated GenBank) records
    # leave id/name as Biopython sentinels. Fall back to the filename
    # so the UI has something meaningful to display.
    stem = _P(path).stem or "plasmid"
    # Sanitize: GenBank LOCUS names can't contain spaces; replace them
    # so round-tripping through _record_to_gb_text doesn't explode.
    safe_stem = stem.replace(" ", "_")[:16] or "plasmid"
    if not rec.id or rec.id.startswith("<unknown"):
        rec.id = safe_stem
    if not rec.name or rec.name.startswith("<unknown"):
        rec.name = safe_stem
    return rec

def _record_to_gb_text(record) -> str:
    """Serialize a SeqRecord to GenBank format text.

    Biopython's genbank writer requires `molecule_type` in annotations
    — if the record came from elsewhere and doesn't have it, default to
    "DNA" rather than crashing. The fill-in happens on a shallow
    SeqRecord copy so the caller's record is never mutated (avoids
    subtle races with concurrent readers and surprise side effects for
    callers that compare records by annotation contents).
    """
    from Bio import SeqIO
    from copy import copy as _shallow_copy
    anns = dict(getattr(record, "annotations", None) or {})
    anns.setdefault("molecule_type", "DNA")
    rec = _shallow_copy(record)
    rec.annotations = anns
    buf = StringIO()
    SeqIO.write(rec, buf, "genbank")
    return buf.getvalue()

def _gb_text_to_record(text: str):
    """Parse GenBank format text back to a SeqRecord."""
    from Bio import SeqIO
    return SeqIO.read(StringIO(text), "genbank")


# ── GenBank export (NCBI-compliant normalization + atomic write) ──────────────
#
# Biopython's SeqIO.write produces compliant GenBank output when annotations
# include the INSDC-mandated fields. `_normalize_for_genbank` fills in the
# minimum set a plasmid record needs to round-trip cleanly through any NCBI-
# compliant parser. The normalization is non-destructive: the caller's
# record is never mutated (shallow copy with fresh annotations dict).
#
# Defaults chosen for a synthetic plasmid of unknown provenance:
#   topology        = circular   (SpliceCraft only works on plasmids)
#   molecule_type   = DNA
#   data_file_division = SYN     (synthetic; matches what NCBI assigns to
#                                 user-submitted plasmids)
#   date            = today, formatted DD-MMM-YYYY
#   organism        = "synthetic construct"
#   taxonomy        = ["other sequences", "artificial sequences"]
#
# These match what NCBI emits for synthetic constructs submitted via BankIt.

_GB_LOCUS_NAME_MAX = 28  # NCBI relaxed LOCUS name length (spec is 16)


def _normalize_for_genbank(record):
    """Return a shallow copy of `record` with NCBI-required fields filled in.

    Idempotent — existing values are preserved. Only fills gaps. Caller's
    record is never mutated.
    """
    from copy import copy as _shallow_copy
    from datetime import datetime as _dt

    rec = _shallow_copy(record)
    anns = dict(getattr(record, "annotations", None) or {})

    anns.setdefault("molecule_type", "DNA")
    anns.setdefault("topology", "circular")
    anns.setdefault("data_file_division", "SYN")

    if not anns.get("date"):
        anns["date"] = _dt.now().strftime("%d-%b-%Y").upper()

    if not anns.get("accessions"):
        acc = rec.id if rec.id and rec.id != "<unknown id>" else ""
        anns["accessions"] = [acc]

    if not anns.get("organism"):
        anns["organism"] = "synthetic construct"
    if not anns.get("source"):
        anns["source"] = anns["organism"]
    if not anns.get("taxonomy"):
        anns["taxonomy"] = ["other sequences", "artificial sequences"]

    rec.annotations = anns

    # LOCUS name: spec is 16 chars; NCBI accepts up to 28 in practice.
    # Biopython itself warns if >16 but does not fail.
    if rec.name and len(rec.name) > _GB_LOCUS_NAME_MAX:
        rec.name = rec.name[:_GB_LOCUS_NAME_MAX]
    if not rec.name or rec.name == "<unknown name>":
        rec.name = (rec.id or "PLASMID")[:_GB_LOCUS_NAME_MAX] or "PLASMID"

    if not rec.description or rec.description == "<unknown description>":
        rec.description = rec.name

    if not rec.id or rec.id == "<unknown id>":
        rec.id = rec.name

    return rec


def _export_genbank_to_path(record, path) -> dict:
    """Write `record` to `path` as a GenBank file. Atomic + round-trip verified.

    Returns a small summary dict `{"path", "bp", "features"}` for UI reporting.

    Raises:
      OSError on filesystem failures (write, replace, fsync).
      ValueError if the round-trip parse fails or the parsed record
        disagrees with the source on sequence length, sequence content,
        or feature count — meaning the export is not byte-safe.

    The round-trip happens BEFORE the target file is touched, so a failed
    export never leaves a half-written / corrupt .gb at `path`.
    """
    import os
    import tempfile
    from pathlib import Path as _Path

    p = _Path(path).expanduser()
    normalized = _normalize_for_genbank(record)
    text = _record_to_gb_text(normalized)

    # Round-trip verify before touching the filesystem
    try:
        parsed = _gb_text_to_record(text)
    except Exception as exc:
        raise ValueError(f"export round-trip parse failed: {exc}") from exc
    if len(parsed.seq) != len(normalized.seq):
        raise ValueError(
            f"export round-trip sequence length mismatch "
            f"({len(parsed.seq)} vs {len(normalized.seq)})"
        )
    if str(parsed.seq).upper() != str(normalized.seq).upper():
        raise ValueError("export round-trip sequence content mismatch")
    if len(parsed.features) != len(normalized.features):
        raise ValueError(
            f"export round-trip feature count mismatch "
            f"({len(parsed.features)} vs {len(normalized.features)})"
        )

    # Atomic write — tempfile in the target's directory, then os.replace.
    p.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(
        prefix=f".{p.name}.", suffix=".tmp", dir=str(p.parent)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            try:
                os.fsync(fh.fileno())
            except OSError:
                pass
        os.replace(tmp, str(p))
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise

    _log.info(
        "Exported GenBank to %s (%d bp, %d features)",
        p, len(normalized.seq), len(normalized.features),
    )
    return {"path": str(p), "bp": len(normalized.seq),
            "features": len(normalized.features)}


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
            # Strip control chars so ANSI/escape sequences from pLannotate's
            # stderr can't corrupt the Textual rendering when surfaced via
            # notify(). Keep tabs/newlines for readability.
            raw_tail = combined_out[-500:]
            err_tail = "".join(
                ch for ch in raw_tail
                if ch in "\t\n" or (ch.isprintable() and ord(ch) >= 0x20)
            ).strip() or "(no output)"
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
        def __init__(self, idx: int, feat_dict: dict | None, bp: int = -1):
            self.idx       = idx
            self.feat_dict = feat_dict
            self.bp        = bp   # bp at click point, or -1 if unknown
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
        # Import once outside the loop to avoid repeated import work.
        try:
            from Bio.SeqFeature import CompoundLocation
        except ImportError:
            CompoundLocation = None
        # Counters the caller can inspect after load_record() to decide
        # whether to notify the user.
        self._n_flattened   = 0
        self._n_skipped     = 0
        self._n_clamped     = 0
        total = len(record.seq) if getattr(record, "seq", None) is not None else 0
        for feat in record.features:
            if feat.type in ("source",):
                continue
            # Biopython's UnknownPosition / BetweenPosition can't be cast to int.
            # Skip such features with a log entry rather than crashing the whole
            # import. Compound locations with fuzzy endpoints (e.g. `<100..200`)
            # ARE castable, so this only catches genuinely unknown coords.
            try:
                start = int(feat.location.start)
                end   = int(feat.location.end)
            except (TypeError, ValueError):
                self._n_skipped += 1
                _log.warning(
                    "Skipped feature %s with non-integer coords (type=%s)",
                    _feat_label(feat), feat.type,
                )
                continue
            strand = getattr(feat.location, "strand", 1) or 1

            # Compound / joined locations need special handling:
            #   * parts that are CONTIGUOUS (each part.end == next.start):
            #     not really split — SnapGene emits these for some
            #     annotation internals; treat as a single plain feature
            #     from [first.start, last.end) with no "flattened" warning
            #     because no information is lost.
            #   * join(450..500, 1..50) on a 500 bp seq → WRAPPING feature;
            #     rebuild as start=450, end=50 (end < start) so the existing
            #     _bp_in / wrap-midpoint machinery renders the correct arc.
            #   * join(100..200, 300..400) on a 500 bp seq (truly split,
            #     e.g. exons of an mRNA) → flatten to 100..400 + warn.
            is_compound = (CompoundLocation is not None
                           and isinstance(feat.location, CompoundLocation))
            if is_compound:
                parts = sorted(
                    feat.location.parts,
                    key=lambda p: int(p.start),
                )
                is_contiguous = all(
                    int(parts[i].end) == int(parts[i + 1].start)
                    for i in range(len(parts) - 1)
                )
                is_wrap = (
                    total > 0 and len(parts) == 2
                    and int(parts[0].start) == 0
                    and int(parts[-1].end) == total
                    and int(parts[0].end) < int(parts[-1].start)
                )
                if is_wrap:
                    # Origin-spanning wrap: head at [0, parts[0].end),
                    # tail at [parts[-1].start, total). Represent as
                    # start=tail, end=head so end<start signals wrap.
                    start = int(parts[-1].start)
                    end   = int(parts[0].end)
                    _log.info(
                        "Detected wrap feature %s (%d..%d → wraps origin)",
                        _feat_label(feat), start, end,
                    )
                elif is_contiguous:
                    # Adjacent sub-parts — outer bounds ARE the real span,
                    # no info lost, no need to warn the user.
                    pass
                else:
                    self._n_flattened += 1
                    _log.info(
                        "Flattened compound feature %s (%d..%d) to outer bounds",
                        _feat_label(feat), start, end,
                    )

            # Clamp coords that exceed the sequence length — rendering
            # math assumes start, end ∈ [0, total].
            if total > 0:
                clamped_start = max(0, min(start, total))
                clamped_end   = max(0, min(end,   total))
                if (clamped_start, clamped_end) != (start, end):
                    self._n_clamped += 1
                    _log.warning(
                        "Clamped feature %s coords (%d..%d → %d..%d) to "
                        "sequence length %d",
                        _feat_label(feat), start, end,
                        clamped_start, clamped_end, total,
                    )
                    start, end = clamped_start, clamped_end

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

    def _feat_at(self, x: int, y: int) -> tuple[int, int]:
        """Return (feature_idx, click_bp) at terminal cell (x, y), or (-1, -1)."""
        if not self.record or not self._total:
            return -1, -1
        w, h = self.size.width, self.size.height
        cx, cy, rx, ry = self._geometry(w, h)
        if rx == 0 or ry == 0:
            return -1, -1
        dc_n = (x - cx) / rx
        dr_n = (y - cy) / ry
        r_norm = math.sqrt(dc_n ** 2 + dr_n ** 2)
        if r_norm < 0.75 or r_norm > 1.35:
            return -1, -1
        angle = math.atan2(dr_n, dc_n)
        angle_norm = (angle + math.pi / 2) % (2 * math.pi)
        bp = int(self.origin_bp + self._total * angle_norm / (2 * math.pi)) % self._total
        for i, f in enumerate(self._feats):
            if self._bp_in(bp, f):
                return i, bp
        return -1, -1

    def on_click(self, event: Click):
        if not self.record:
            return
        if self._map_mode == "linear":
            idx, bp = self._feat_at_linear(event.x, event.y)
        else:
            idx, bp = self._feat_at(event.x, event.y)
        self.selected_idx = idx
        f = self._feats[idx] if idx >= 0 else None
        self.post_message(self.FeatureSelected(idx, f, bp))

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
            key=lambda iv: -_feat_len(iv[1]["start"], iv[1]["end"], total),
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

    def _feat_at_linear(self, x: int, y: int) -> tuple[int, int]:
        """Return (feature_idx, click_bp) at terminal cell (x, y) in linear view,
        or (-1, -1) if outside the map region / no feature matched."""
        if not self._total:
            return -1, -1
        w, h      = self.size.width, self.size.height
        margin_l  = 5
        margin_r  = 2
        usable_w  = w - margin_l - margin_r
        backbone_row = max(4, h // 2)
        if x < margin_l or x >= w - margin_r or usable_w <= 0:
            return -1, -1
        bp = int((x - margin_l) / usable_w * self._total)
        above = y < backbone_row
        below = y > backbone_row
        for i, f in enumerate(self._feats):
            # Use half-open [start, end) to match _bp_in elsewhere. This
            # also makes zero-width features (s == e) unclickable instead
            # of matching every column on the backbone.
            if not self._bp_in(bp, f):
                continue
            if above and f["strand"] >= 0:
                return i, bp
            if below and f["strand"] < 0:
                return i, bp
        return -1, -1

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
        try:
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
        except Exception:
            # Malformed feature dict (missing color/type/label) would otherwise
            # leave _populating=True forever, blocking every RowHighlighted
            # event and breaking the sidebar until app restart. Reset the flag
            # synchronously before re-raising so the cascade can recover.
            self._populating = False
            _log.exception("FeatureSidebar.populate failed mid-loop")
            raise
        def _clear_populating():
            self._populating = False
        self.call_after_refresh(_clear_populating)

    def show_detail(self, f: dict | None) -> None:
        box = self.query_one("#detail-box", Static)
        if f is None:
            box.update(Text(""))
            return
        strand_sym = "+" if f["strand"] == 1 else ("−" if f["strand"] == -1 else "·")
        # Wrap features need modular length (end < start means the feature
        # crosses the origin); naive end-start is negative and misleading.
        rec = getattr(self.app, "_current_record", None)
        total = len(rec.seq) if rec else max(f["end"], f["start"]) + 1
        span = _feat_len(f["start"], f["end"], total)
        if f["end"] < f["start"]:
            coord_str = f"{f['start']+1}‥{total},1‥{f['end']}"
        else:
            coord_str = f"{f['start']+1}‥{f['end']}"
        t = Text()
        t.append(f["type"],  style=f"bold {f['color']}")
        t.append("\n")
        t.append(f["label"], style="white")
        t.append("\n")
        t.append(f"{coord_str} ({span:,} bp)", style="dim")
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
            yield Button("+", id="btn-lib-add", variant="primary",
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
            key=lambda f: -_feat_len(f["start"], f["end"], n),
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
            chunk_feats = _feats_in_chunk(annot_feats, chunk_start, chunk_end, n)
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
            n = len(self._seq)
            self._sorted_feats_cache = sorted(
                [f for f in self._feats if f.get("type") not in ("site", "recut")],
                key=lambda f: -_feat_len(f["start"], f["end"], n),
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
            chunk_feats = _feats_in_chunk(annot_feats, chunk_start, chunk_end, n)
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
        except NoMatches:
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
            except NoMatches:
                return
            vp_top = int(scroll.scroll_y)
            vp_h   = scroll.size.height
            vp_bottom = vp_top + vp_h - 1
            if row < vp_top:
                scroll.scroll_to(0, row, animate=False)
            elif row_bottom > vp_bottom:
                scroll.scroll_to(0, row_bottom - vp_h + 1, animate=False)

        self.call_after_refresh(_do_scroll)

    def _refresh_view(self) -> None:
        view = self.query_one("#seq-view", Static)
        try:
            scroll = self.query_one("#seq-scroll", ScrollableContainer)
        except NoMatches:
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

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

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
    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

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
        except Exception as exc:
            _log.exception("NCBI fetch failed for %s", acc)
            def _err():
                # Modal may have been dismissed while the fetch was in flight;
                # query_one would then raise NoMatches. Fall back to a toast.
                if not self.is_mounted:
                    try:
                        self.app.notify(f"NCBI fetch failed: {exc}",
                                        severity="error", timeout=8)
                    except Exception:
                        _log.exception("notify fallback for fetch error failed")
                    return
                try:
                    self.query_one("#fetch-status", Static).update(
                        f"[red]Error: {exc}[/red]"
                    )
                except NoMatches:
                    try:
                        self.app.notify(f"NCBI fetch failed: {exc}",
                                        severity="error", timeout=8)
                    except Exception:
                        _log.exception("notify fallback for fetch error failed")
            self.app.call_from_thread(_err)
            return

        # Staleness guard: user may have hit Escape (which runs `dismiss(None)`
        # via action_cancel) while the HTTP round-trip was in flight. Calling
        # `self.dismiss(record)` on an already-dismissed modal raises, and even
        # if it didn't, we'd be stomping whatever the user did next. Silently
        # drop the fetched record if the modal is no longer mounted.
        def _apply():
            if not self.is_mounted:
                _log.info(
                    "Fetch %s completed after modal was dismissed; discarding.",
                    acc,
                )
                return
            self.dismiss(record)
        self.app.call_from_thread(_apply)

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
    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def compose(self) -> ComposeResult:
        with Vertical(id="open-box"):
            yield Static(" Open Plasmid File ", id="open-title")
            yield Label("File path  (.gb / .gbk / .dna):")
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


# ── Export-GenBank modal ───────────────────────────────────────────────────────

class ExportGenBankModal(ModalScreen):
    """Prompt for a target path, write the current record as GenBank.

    On submit, dismisses with a summary dict from `_export_genbank_to_path`
    (keys: path, bp, features) or None if cancelled. The caller is
    expected to show a notification on success.
    """

    BINDINGS = [
        Binding("escape", "cancel",      "Cancel"),
        Binding("tab",    "focus_next",  "Next",   show=False),
    ]

    def __init__(self, record, default_path: str = "") -> None:
        super().__init__()
        self._record = record
        self._default_path = default_path

    def compose(self) -> ComposeResult:
        with Vertical(id="export-box"):
            yield Static(" Export as GenBank ", id="export-title")
            name = getattr(self._record, "name", "") or "plasmid"
            bp = len(getattr(self._record, "seq", "") or "")
            feats = len(getattr(self._record, "features", []) or [])
            yield Label(f"[{name}]  {bp} bp, {feats} features")
            yield Label("Output path (.gb):")
            yield Input(
                value=self._default_path,
                placeholder="/path/to/plasmid.gb",
                id="export-path",
            )
            with Horizontal(id="export-btns"):
                yield Button("Export",  id="btn-export", variant="primary")
                yield Button("Cancel",  id="btn-cancel-export")
            yield Static("", id="export-status", markup=True)

    def on_mount(self) -> None:
        try:
            self.query_one("#export-path", Input).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-export")
    def _do_export(self) -> None:
        try:
            inp = self.query_one("#export-path", Input)
            status = self.query_one("#export-status", Static)
        except NoMatches:
            return
        path = inp.value.strip()
        if not path:
            status.update("[red]Enter an output path.[/red]")
            return
        try:
            summary = _export_genbank_to_path(self._record, path)
        except (OSError, ValueError) as exc:
            _log.exception("GenBank export to %s failed", path)
            status.update(f"[red]Export failed: {exc}[/red]")
            return
        self.dismiss(summary)

    @on(Input.Submitted)
    def _submitted(self) -> None:
        try:
            self.query_one("#btn-export", Button).press()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-cancel-export")
    def _cancel_btn(self) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
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

    MENUS = ["File", "Edit", "Enzymes", "Features", "Primers", "Mutagenize",
             "Parts", "Constructor"]

    def compose(self) -> ComposeResult:
        for name in self.MENUS:
            yield Static(name, classes="menu-item", id=f"menu-{name.lower()}")

    def on_click(self, event: Click) -> None:
        for name in self.MENUS:
            widget_id = f"menu-{name.lower()}"
            try:
                item = self.query_one(f"#{widget_id}", Static)
            except NoMatches:
                continue
            region = item.region
            if (region.x <= event.screen_x < region.x + region.width and
                    region.y <= event.screen_y < region.y + region.height):
                x = region.x
                y = region.y + 1
                self.app.open_menu(name, x, y)
                break


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


def _pick_binding_region(seq: str, target_tm: float = 60.0,
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

    # Defensive init: if the caller forgot the len(seq) >= min_len guard,
    # the loop below won't execute and we'd otherwise return Tm=0 with a
    # too-short binding. Compute Tm for whatever is there so downstream
    # validation (low Tm, short primer) still trips honestly.
    best_seq = seq[:max(min_len, 1)]
    best_tm  = _tm(best_seq) if best_seq else 0.0
    best_diff = float("inf")
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
    target_tm: float = 60.0,
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
    total  = len(template_seq)
    insert = _slice_circular(template_seq.upper(), start, end)
    wraps  = end < start

    # Need at least 18 bp to pick a proper binding region — otherwise
    # _pick_binding_region returns the whole (too-short) insert with Tm=0.
    if len(insert) < 18:
        return {
            "error": f"Golden Braid region is too short ({len(insert)} bp). "
                     f"Select at least 18 bp (recommended 25+ bp for a "
                     f"robust binding region).",
        }

    # Internal BsaI site check — a GGTCTC (or GAGACC on the bottom strand)
    # inside the insert will be cut during domestication and fragment the
    # part. GB L0 parts must be domesticated to remove these; catching it
    # here lets the user pick another region or redesign their template.
    bsai_rc  = _rc(_GB_BSAI_SITE)
    bad_fwd  = insert.find(_GB_BSAI_SITE)
    bad_rev  = insert.find(bsai_rc)
    if bad_fwd != -1 or bad_rev != -1:
        hits = []
        if bad_fwd != -1:
            hits.append(f"{_GB_BSAI_SITE} at +{bad_fwd + 1}")
        if bad_rev != -1:
            hits.append(f"{bsai_rc} at +{bad_rev + 1}")
        return {
            "error": f"Internal BsaI site found ({'; '.join(hits)}). "
                     f"The part would self-digest during domestication — "
                     f"pick a different region or silently mutate the site first.",
        }

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

    # Positions of the primer binding regions on the TEMPLATE (not the
    # full amplicon). The forward primer binds the top strand at the
    # start of the insert; the reverse primer binds the bottom strand at
    # the end of the insert (positions are reported in forward-strand
    # coordinates). Save-to-library needs these to add primer_bind
    # features to the map. For wrap regions, compute positions with
    # modular arithmetic so they land on the real plasmid coordinates.
    if wraps:
        fwd_pos = (start, (start + len(fwd_bind)) % total)
        rev_pos = ((end - len(rev_bind)) % total, end)
    else:
        fwd_pos = (start, start + len(fwd_bind))
        rev_pos = (end - len(rev_bind), end)

    return {
        "part_type":    part_type,
        "position":     pos_label,
        "oh5":          oh5,
        "oh3":          oh3,
        "insert_seq":   insert,
        "fwd_binding":  fwd_bind,
        "rev_binding":  rev_bind,
        "fwd_full":     fwd_full,
        "rev_full":     rev_full,
        "fwd_tm":       round(fwd_tm, 1),
        "rev_tm":       round(rev_tm, 1),
        "fwd_pos":      fwd_pos,
        "rev_pos":      rev_pos,
        "amplicon_len": amplicon_len,
    }


# ── Parts bin persistence ─────────────────────────────────────────────────────
# User-created parts (from the domesticator) are stored in parts_bin.json next
# to the main script. Each entry is a dict with at least the 7 canonical fields
# plus sequence, primers, and Tm values.

_PARTS_BIN_FILE = _DATA_DIR / "parts_bin.json"
_parts_bin_cache: "list | None" = None

def _load_parts_bin() -> list[dict]:
    global _parts_bin_cache
    if _parts_bin_cache is not None:
        return list(_parts_bin_cache)
    entries, warning = _safe_load_json(_PARTS_BIN_FILE, "Parts bin")
    if warning:
        _log.warning(warning)
    entries = [e for e in entries if isinstance(e, dict)]
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

_PRIMERS_FILE = _DATA_DIR / "primers.json"
_primers_cache: "list | None" = None

def _load_primers() -> list[dict]:
    global _primers_cache
    if _primers_cache is not None:
        return list(_primers_cache)
    entries, warning = _safe_load_json(_PRIMERS_FILE, "Primer library")
    if warning:
        _log.warning(warning)
    entries = [e for e in entries if isinstance(e, dict)]
    _primers_cache = entries
    return list(_primers_cache)

def _save_primers(entries: list[dict]) -> None:
    global _primers_cache
    _safe_save_json(_PRIMERS_FILE, entries, "Primer library")
    _primers_cache = list(entries)


# ── Feature library persistence ───────────────────────────────────────────────
#
# User-created features (from the Add Feature modal or extracted from existing
# plasmids) persist to features.json. An entry is a dict:
#
#     {
#       "name":         "lacZ-alpha",           # label for the sidebar
#       "feature_type": "CDS",                  # INSDC feature-table type
#       "sequence":     "ATG...TAA",            # 5'→3' forward-strand DNA
#       "strand":       1 | -1,                 # direction when inserted
#       "qualifiers":   {"gene": ["lacZ"], ...}, # GenBank-compatible qualifiers
#       "description":  "(free text)"           # optional
#     }
#
# The `feature_type` is validated against a curated list of INSDC-standard
# types (see `_GENBANK_FEATURE_TYPES`). Unknown types are accepted with a
# warning — some in-house projects use non-standard types and we shouldn't
# block them, but the warning flags risk of downstream tool incompatibility.

_FEATURES_FILE = _DATA_DIR / "features.json"
_features_cache: "list | None" = None

# Curated INSDC / GenBank feature-table types relevant to plasmid work.
# Full spec: https://www.insdc.org/submitting-standards/feature-table/
# Ordered by frequency so the modal dropdown puts CDS / gene / promoter at the
# top. The `source` type is excluded because each record already has exactly
# one `source` feature spanning the whole molecule — adding more would be
# invalid GenBank.
_GENBANK_FEATURE_TYPES: tuple = (
    "CDS", "gene", "mRNA", "tRNA", "rRNA", "ncRNA", "misc_RNA",
    "promoter", "terminator", "RBS", "polyA_signal", "regulatory",
    "5'UTR", "3'UTR", "intron", "exon", "operon",
    "primer_bind", "protein_bind", "misc_binding",
    "repeat_region", "LTR", "mobile_element",
    "rep_origin", "oriT",
    "sig_peptide", "mat_peptide", "transit_peptide", "propeptide",
    "misc_feature", "misc_recomb", "stem_loop", "variation",
)


def _load_features() -> list[dict]:
    global _features_cache
    if _features_cache is not None:
        return list(_features_cache)
    entries, warning = _safe_load_json(_FEATURES_FILE, "Feature library")
    if warning:
        _log.warning(warning)
    entries = [e for e in entries if isinstance(e, dict)]
    _features_cache = entries
    return list(_features_cache)


def _save_features(entries: list[dict]) -> None:
    global _features_cache
    _safe_save_json(_FEATURES_FILE, entries, "Feature library")
    _features_cache = list(entries)


def _parse_qualifier_string(raw: str) -> dict[str, list[str]]:
    """Parse a user-typed qualifier line like ``gene=lacZ; product=LacZ alpha``
    into the GenBank-style ``{key: [value]}`` dict.

    - Separator between pairs is ``;`` (newlines also accepted).
    - Separator between key and value is the first ``=``.
    - Keys and values are stripped; empty keys are skipped silently.
    - Duplicate keys are merged into a single list (GenBank allows repeated
      qualifiers — e.g. multiple ``/note=`` lines).
    """
    out: dict[str, list[str]] = {}
    if not raw:
        return out
    chunks = re.split(r"[;\n]", raw)
    for chunk in chunks:
        if "=" not in chunk:
            continue
        key, _, val = chunk.partition("=")
        key = key.strip()
        val = val.strip()
        if not key:
            continue
        # Permissive: accept an empty value. GenBank flags (e.g. /pseudo)
        # are technically value-less, but splicecraft always writes a string
        # form, so keep them as empty-string qualifiers.
        out.setdefault(key, []).append(val)
    return out


def _qualifiers_to_string(quals: dict) -> str:
    """Inverse of `_parse_qualifier_string`. Used to prefill the modal when
    the user imports an existing feature. Multi-value qualifiers get one
    `key=value` pair per occurrence."""
    if not quals:
        return ""
    pieces: list[str] = []
    for key, vals in quals.items():
        if isinstance(vals, (list, tuple)):
            for v in vals:
                pieces.append(f"{key}={v}")
        else:
            pieces.append(f"{key}={vals}")
    return "; ".join(pieces)


def _extract_feature_entries_from_record(record) -> list[dict]:
    """Return one feature-library entry dict per non-source feature.

    Wrap features (origin-spanning CompoundLocations) are flattened into the
    forward-strand genomic sequence before export, so the entry's ``sequence``
    is always the 5'→3' DNA that would be re-inserted. Reverse-strand
    features store the revcomp (i.e. the 5'→3' of the feature as read), which
    matches how the Add Feature modal expects input.
    """
    try:
        from Bio.SeqFeature import CompoundLocation
    except ImportError:
        CompoundLocation = tuple()  # type: ignore[assignment]
    seq = str(getattr(record, "seq", "") or "").upper()
    total = len(seq)
    entries: list[dict] = []
    for feat in getattr(record, "features", []) or []:
        if feat.type == "source":
            continue
        loc = feat.location
        strand = getattr(loc, "strand", 1) or 1
        # Assemble the forward-strand genomic sequence under the feature,
        # respecting wrap/compound locations.
        if isinstance(loc, CompoundLocation):
            parts_seq = []
            for part in loc.parts:
                s = int(part.start) % total if total else 0
                e = int(part.end)   % (total or 1) if total else 0
                if total and e <= s:
                    parts_seq.append(seq[s:] + seq[:e])
                else:
                    parts_seq.append(seq[s:e])
            fwd = "".join(parts_seq)
        else:
            s = int(loc.start)
            e = int(loc.end)
            fwd = seq[s:e]
        # Store 5'→3' of the feature as read. For reverse-strand CDS that is
        # the revcomp of the genomic slice.
        if strand == -1 and fwd:
            feat_seq = _rc(fwd)
        else:
            feat_seq = fwd
        entries.append({
            "name":         _feat_label(feat),
            "feature_type": feat.type,
            "sequence":     feat_seq,
            "strand":       1 if strand != -1 else -1,
            "qualifiers":   {k: list(v) if isinstance(v, (list, tuple)) else [v]
                             for k, v in (feat.qualifiers or {}).items()},
            "description":  "",
        })
    return entries


# ── Codon usage registry + harmonization (shared across modals) ───────────────
#
# Persistent JSON library of codon usage tables, plus pure-function harmonizer
# and restriction-site fixer. Vendored from superfolder_aeblue/codon_tables.py
# and codon_harmonize.py (2026-04-13) to keep the single-file convention.
#
# Storage schema (one entry per species in codon_tables.json):
#   {"name": "E. coli K12", "taxid": "83333", "source": "builtin"|"kazusa"|"user",
#    "added": "YYYY-MM-DD", "raw": {"GCT": ["A", 55], ...}}
# In-memory form uses tuples (aa, count); the "raw" JSON list is converted on
# load. Dedup key is the taxid when present, else the display name.

_CODON_GENETIC_CODE: dict[str, str] = {
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

_CODON_BUILTIN_K12: dict[str, tuple[str, int]] = {
    "GGG": ("G",  44), "GGA": ("G",  47), "GGT": ("G", 109), "GGC": ("G", 171),
    "GAG": ("E",  94), "GAA": ("E", 224), "GAT": ("D", 194), "GAC": ("D", 105),
    "GTG": ("V", 135), "GTA": ("V",  59), "GTT": ("V",  86), "GTC": ("V",  60),
    "GCG": ("A", 197), "GCA": ("A", 108), "GCT": ("A",  55), "GCC": ("A", 162),
    "AGG": ("R",   8), "AGA": ("R",   7), "AGT": ("S",  37), "AGC": ("S",  85),
    "AAG": ("K",  62), "AAA": ("K", 170), "AAT": ("N", 112), "AAC": ("N", 125),
    "ATG": ("M", 127), "ATA": ("I",  19), "ATT": ("I", 156), "ATC": ("I",  93),
    "ACG": ("T",  59), "ACA": ("T",  33), "ACT": ("T",  41), "ACC": ("T", 117),
    "TGG": ("W",  55), "TGT": ("C",  30), "TGC": ("C",  41),
    "TAT": ("Y",  86), "TAC": ("Y",  75),
    "TTG": ("L",  61), "TTA": ("L",  78), "TTT": ("F", 101), "TTC": ("F",  77),
    "TCG": ("S",  41), "TCA": ("S",  40), "TCT": ("S",  29), "TCC": ("S",  28),
    "CGG": ("R",  21), "CGA": ("R",  22), "CGT": ("R", 108), "CGC": ("R", 133),
    "CAG": ("Q", 142), "CAA": ("Q",  62), "CAT": ("H",  81), "CAC": ("H",  67),
    "CTG": ("L", 240), "CTA": ("L",  27), "CTT": ("L",  61), "CTC": ("L",  54),
    "CCG": ("P", 137), "CCA": ("P",  34), "CCT": ("P",  43), "CCC": ("P",  33),
    "TAA": ("*",   9), "TAG": ("*",   0), "TGA": ("*",   5),
}

# Forbidden sites for the harmonizer's restriction-site fixer. Keys are the
# forward site only; the fixer adds the reverse complement automatically if
# the site is non-palindromic.
_CODON_DEFAULT_FORBIDDEN: dict[str, str] = {
    "BsaI":    "GGTCTC",
    "BsmBI":   "CGTCTC",
    "BbsI":    "GAAGAC",
    "EcoRI":   "GAATTC",
    "NdeI":    "CATATG",
    "XhoI":    "CTCGAG",
    "BamHI":   "GGATCC",
    "HindIII": "AAGCTT",
    "NcoI":    "CCATGG",
    "SalI":    "GTCGAC",
    "KpnI":    "GGTACC",
    "SacI":    "GAGCTC",
}

_CODON_TABLES_FILE = _DATA_DIR / "codon_tables.json"
_codon_tables_cache: "list | None" = None

# ── Crash-recovery autosave ────────────────────────────────────────────────────
# Every edit marks the record dirty, which (debounced) writes the current record
# to _CRASH_RECOVERY_DIR/{safe_id}.gb. A successful save or clean-state flip
# deletes the file. On startup we scan the dir and notify the user if any
# leftover autosaves exist (i.e. the last session crashed before save).
_CRASH_RECOVERY_DIR = _DATA_DIR / "crash_recovery"


def _codon_raw_to_json(raw: dict) -> dict:
    """Convert in-memory {codon: (aa, count)} → JSON-safe {codon: [aa, count]}."""
    return {c: [aa, int(n)] for c, (aa, n) in raw.items()}


def _codon_raw_from_json(blob: dict) -> dict:
    """Inverse of _codon_raw_to_json. Accepts either tuples or 2-item lists."""
    out: dict = {}
    for c, v in blob.items():
        if isinstance(v, (list, tuple)) and len(v) == 2:
            out[c.upper()] = (str(v[0]), int(v[1]))
    return out


def _codon_tables_load() -> list[dict]:
    """Load codon-table registry. Returns a list of dicts; each entry has
    an extra 'raw' in-memory form with tuples. Seeds built-in E. coli K12
    on first run so the library is never empty."""
    global _codon_tables_cache
    if _codon_tables_cache is not None:
        return list(_codon_tables_cache)
    entries, warning = _safe_load_json(_CODON_TABLES_FILE, "Codon table library")
    if warning:
        _log.warning("Codon table library: %s", warning)
    fixed: list = []
    for e in entries:
        if not isinstance(e, dict):
            continue
        raw = _codon_raw_from_json(e.get("raw", {}))
        if not raw:
            continue
        fixed.append({
            "name":   e.get("name", "?"),
            "taxid":  str(e.get("taxid", "")),
            "source": e.get("source", "user"),
            "added":  e.get("added", ""),
            "raw":    raw,
        })
    # Seed built-in K12 if not present
    if not any(e.get("taxid") == "83333" for e in fixed):
        import datetime
        fixed.insert(0, {
            "name":   "E. coli K12",
            "taxid":  "83333",
            "source": "builtin",
            "added":  datetime.date.today().isoformat(),
            "raw":    dict(_CODON_BUILTIN_K12),
        })
        _codon_tables_cache = fixed
        _codon_tables_save(fixed)
    else:
        _codon_tables_cache = fixed
    return list(_codon_tables_cache)


def _codon_tables_save(entries: list[dict]) -> None:
    """Persist registry to disk via _safe_save_json (atomic, .bak)."""
    global _codon_tables_cache
    serializable = [{
        "name":   e.get("name", "?"),
        "taxid":  str(e.get("taxid", "")),
        "source": e.get("source", "user"),
        "added":  e.get("added", ""),
        "raw":    _codon_raw_to_json(e.get("raw", {})),
    } for e in entries]
    _safe_save_json(_CODON_TABLES_FILE, serializable, "Codon table library")
    _codon_tables_cache = list(entries)


def _codon_tables_add(name: str, taxid: str, raw: dict,
                      source: str = "user") -> dict:
    """Insert or replace a table in the registry. Dedup key is taxid when
    non-empty, else name. Returns the stored entry."""
    import datetime
    entries = _codon_tables_load()
    taxid = str(taxid or "").strip()
    name  = (name or "?").strip() or "?"
    entry = {
        "name":   name,
        "taxid":  taxid,
        "source": source,
        "added":  datetime.date.today().isoformat(),
        "raw":    dict(raw),
    }
    def _same(e):
        if taxid and e.get("taxid") == taxid:
            return True
        if not taxid and e.get("name") == name:
            return True
        return False
    kept = [e for e in entries if not _same(e)]
    kept.append(entry)
    _codon_tables_save(kept)
    return entry


def _codon_tables_get(key: str) -> "dict | None":
    """Look up a table by taxid or name (case-insensitive). Returns the
    in-memory entry (with 'raw' as tuples) or None."""
    key = (key or "").strip().lower()
    if not key:
        return None
    for e in _codon_tables_load():
        if str(e.get("taxid", "")).lower() == key:
            return e
        if str(e.get("name", "")).lower() == key:
            return e
    return None


def _codon_name_parts(name: str) -> tuple[str, str]:
    """Return (genus, species) as lowercased tokens from an entry name.

    Genus = first whitespace-delimited token; species = second token (or "").
    Names like "E. coli K12" yield ("e.", "coli"); "Escherichia coli" yields
    ("escherichia", "coli"). No normalization between abbreviated and
    unabbreviated genera — users search what they see.
    """
    parts = (name or "").strip().split()
    genus   = parts[0].lower() if parts else ""
    species = parts[1].lower() if len(parts) > 1 else ""
    return genus, species


def _codon_search(query: str, entries: "list | None" = None) -> list[dict]:
    """Ranked search over taxid, genus, species, and full name.

    Rank 0: taxid exact match
    Rank 1: taxid prefix match
    Rank 2: genus prefix (first whitespace token of name)
    Rank 3: species prefix (second whitespace token)
    Rank 4: substring anywhere in the full name

    Results are sorted by (rank, name) so same-genus entries cluster and
    the strongest match wins. An empty/whitespace query returns the
    registry in its persisted order, unchanged, to preserve caller
    expectations (the filter field shows "all tables" when cleared).
    """
    if entries is None:
        entries = _codon_tables_load()
    q = (query or "").strip().lower()
    if not q:
        return list(entries)
    ranked: list[tuple[int, str, dict]] = []
    for e in entries:
        name_lc  = str(e.get("name", "")).lower()
        taxid_lc = str(e.get("taxid", "")).lower()
        genus, species = _codon_name_parts(e.get("name", ""))
        if taxid_lc and taxid_lc == q:
            rank = 0
        elif taxid_lc and taxid_lc.startswith(q):
            rank = 1
        elif genus and genus.startswith(q):
            rank = 2
        elif species and species.startswith(q):
            rank = 3
        elif q in name_lc:
            rank = 4
        else:
            continue
        ranked.append((rank, name_lc, e))
    ranked.sort(key=lambda t: (t[0], t[1]))
    return [t[2] for t in ranked]


def _codon_parse_kazusa_html(html: str) -> "dict | None":
    """Parse Kazusa showcodon.cgi GCG-format HTML. Returns {codon: (aa, count)}
    or None on failure."""
    pre = re.search(r"<[Pp][Rr][Ee]>(.*?)</[Pp][Rr][Ee]>", html, re.DOTALL)
    text = pre.group(1) if pre else html
    pat = re.compile(r"\b([ACGTU]{3})\b\s+(\d+(?:\.\d+)?)", re.IGNORECASE)
    raw: dict = {}
    for m in pat.finditer(text):
        rna = m.group(1).upper()
        dna = rna.replace("U", "T")
        if dna not in _CODON_GENETIC_CODE or dna in raw:
            continue
        try:
            count = round(float(m.group(2)))
        except ValueError:
            continue
        raw[dna] = (_CODON_GENETIC_CODE[dna], count)
    if len([c for c in raw if raw[c][0] != "?"]) < 60:
        return None
    return raw


def _safe_xml_parse(xml_data: str):
    """Parse XML with basic defense against billion-laughs / XXE tricks.

    Python's stdlib ET (expat) already refuses to fetch external entities
    since 3.7.1, so the remaining attack surface is DTD-declared entity
    expansion. Rejecting any DOCTYPE/ENTITY up front is a dep-free guard
    against a compromised NCBI mirror or MITM. NCBI's real responses
    have no DTD, so this never false-positives.
    """
    import xml.etree.ElementTree as ET
    head = xml_data[:4096].lower()
    if "<!doctype" in head or "<!entity" in head:
        raise ET.ParseError("XML contains DTD/ENTITY — refusing to parse")
    return ET.fromstring(xml_data)


def _ncbi_prep_term(query: str) -> str:
    """Turn a user query into an NCBI Entrez term.

    * Single token (typed-as-genus): combine an exact-taxon subtree search
      restricted to species rank with a wildcard prefix search via OR, so
      typing 'Escherichia' returns every Escherichia species (subtree hit)
      AND typing a partial like 'Escher' still matches via the wildcard.
    * Multi-word query: append '*' to the trailing token so 'Homo sapien'
      matches 'Homo sapiens' etc.
    * User-supplied wildcards or field tags pass through untouched.
    """
    q = (query or "").strip()
    if not q or "*" in q or "[" in q:
        return q
    tokens = q.split()
    if len(tokens) == 1:
        t = tokens[0]
        return f"({t}[Subtree] AND species[Rank]) OR {t}*"
    tokens[-1] = tokens[-1] + "*"
    return " ".join(tokens)


def _ncbi_taxid_search(query: str, retmax: int = 200,
                      timeout: float = 15.0) -> tuple:
    """Search NCBI taxonomy for candidates matching `query`. Returns
    (hits, total_count, status_message) where each hit is
    {'taxid': str, 'name': str}. Names come from a batched esummary call
    (one round-trip for up to `retmax` ids). Partial queries are auto-
    wildcarded via `_ncbi_prep_term`. Pure network — run from a worker."""
    import urllib.request, urllib.parse
    import xml.etree.ElementTree as ET
    q = (query or "").strip()
    if not q:
        return [], 0, "Empty query"
    term = _ncbi_prep_term(q)
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    params = urllib.parse.urlencode({
        "db": "taxonomy", "term": term,
        "retmax": str(retmax), "retmode": "xml",
    })
    try:
        req = urllib.request.Request(f"{base}/esearch.fcgi?{params}",
                                     headers={"User-Agent": "SpliceCraft/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            xml_data = r.read().decode("utf-8", errors="replace")
    except Exception as exc:
        _log.exception("NCBI esearch failed for %r", q)
        return [], 0, f"Network error: {exc}"
    try:
        root = _safe_xml_parse(xml_data)
    except ET.ParseError as exc:
        return [], 0, f"Could not parse NCBI response: {exc}"
    ids = [e.text for e in root.findall(".//Id") if e.text]
    count_elem = root.find(".//Count")
    try:
        total = int(count_elem.text) if count_elem is not None and count_elem.text else len(ids)
    except ValueError:
        total = len(ids)
    if not ids:
        return [], 0, f"No NCBI taxonomy entry for '{q}'"
    # Batched esummary: one round-trip for all retrieved ids
    names_by_id: dict[str, str] = {}
    try:
        sparams = urllib.parse.urlencode({
            "db": "taxonomy", "id": ",".join(ids), "retmode": "xml",
        })
        req = urllib.request.Request(f"{base}/esummary.fcgi?{sparams}",
                                     headers={"User-Agent": "SpliceCraft/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            sxml = r.read().decode("utf-8", errors="replace")
        sroot = _safe_xml_parse(sxml)
        for doc in sroot.findall(".//DocSum"):
            did_el = doc.find("Id")
            if did_el is None or not did_el.text:
                continue
            did = did_el.text
            for item in doc.findall("Item"):
                if item.get("Name") == "ScientificName" and item.text:
                    names_by_id[did] = item.text
                    break
    except Exception:
        _log.exception("NCBI esummary failed for ids %s", ids[:3])
    hits = [{"taxid": tid,
             "name":  names_by_id.get(tid, f"(taxid {tid})")}
            for tid in ids]
    msg = f"{total} hit(s) for '{q}'"
    if total > len(hits):
        msg = f"Showing {len(hits)} of {total} hits for '{q}' (refine to narrow)"
    return hits, total, msg


def _codon_fetch_kazusa(taxid: str, timeout: float = 15.0) -> tuple:
    """Fetch codon usage from Kazusa for an NCBI taxid. Returns
    (raw_dict_or_None, status_message). Pure network call — callers should
    invoke from a worker thread."""
    import urllib.request
    taxid = str(taxid).strip()
    if not taxid.isdigit():
        return None, f"Invalid taxid '{taxid}' (must be numeric)"
    url = (f"https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi"
           f"?species={taxid}&aa=1&style=GCG")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            html = r.read().decode("utf-8", errors="replace")
    except Exception as exc:
        _log.exception("Kazusa fetch failed for taxid %s", taxid)
        return None, f"Network error: {exc}"
    low = html.lower()
    if "not found" in low or "no data" in low:
        return None, f"Taxid {taxid} not found in Kazusa database"
    raw = _codon_parse_kazusa_html(html)
    if raw is None:
        return None, f"Could not parse Kazusa table for taxid {taxid}"
    return raw, f"Fetched from Kazusa: {len(raw)} codons (taxid {taxid})"


def _codon_build_aa_map(raw: dict) -> tuple[dict, dict]:
    """Given {codon: (aa, count)}, return (aa_codons, codon_frac) where
    aa_codons[aa] = [(codon, frac), ...] sorted by fraction descending, and
    codon_frac[codon] = fractional usage for its amino acid."""
    from collections import defaultdict
    aa_total: dict = defaultdict(int)
    for codon, (aa, count) in raw.items():
        aa_total[aa] += int(count)
    codon_frac: dict = {}
    for codon, (aa, count) in raw.items():
        total = aa_total.get(aa, 0) or 1
        codon_frac[codon] = count / total
    aa_codons: dict = defaultdict(list)
    for codon, (aa, count) in raw.items():
        if aa == "*":
            continue
        aa_codons[aa].append((codon, codon_frac[codon]))
    for aa in aa_codons:
        aa_codons[aa].sort(key=lambda x: -x[1])
    return dict(aa_codons), codon_frac


def _codon_harmonize(protein: str, raw: dict) -> str:
    """Distribute synonymous codons across the protein so each amino acid's
    codon distribution matches the target organism (Angov 2011, adapted).
    Appends a TAA stop. Raises ValueError on unknown amino acids."""
    aa_codons, _ = _codon_build_aa_map(raw)
    aa_positions: dict = {}
    for i, aa in enumerate(protein):
        aa = aa.upper()
        aa_positions.setdefault(aa, []).append(i)
    codon_at = [""] * len(protein)
    for aa, positions in aa_positions.items():
        n = len(positions)
        codons_for_aa = aa_codons.get(aa, [])
        if not codons_for_aa:
            raise ValueError(f"No codons for amino acid '{aa}' in this table")
        if len(codons_for_aa) == 1:
            for pos in positions:
                codon_at[pos] = codons_for_aa[0][0]
            continue
        targets: list = []
        remainders: list = []
        allocated = 0
        for codon, frac in codons_for_aa:
            exact   = n * frac
            floored = int(exact)
            targets.append(floored)
            remainders.append((exact - floored, len(targets) - 1))
            allocated += floored
        shortage = n - allocated
        remainders.sort(key=lambda x: -x[0])
        for k in range(shortage):
            targets[remainders[k][1]] += 1
        queues = [[codon] * cnt
                  for (codon, _), cnt in zip(codons_for_aa, targets) if cnt > 0]
        interleaved: list = []
        i = 0
        while any(queues):
            q = queues[i % len(queues)]
            if q:
                interleaved.append(q.pop(0))
            i += 1
        for pos, codon in zip(positions, interleaved):
            codon_at[pos] = codon
    return "".join(codon_at) + "TAA"


def _codon_fix_sites(dna: str, protein: str, raw: dict,
                     sites: "dict | None" = None) -> tuple:
    """Substitute synonymous codons to remove internal restriction sites.
    `sites` is a forward-strand {name: site} dict; reverse complements are
    added automatically for non-palindromic sites. Returns (new_dna, fixes)."""
    if sites is None:
        sites = _CODON_DEFAULT_FORBIDDEN
    expanded: dict = {}
    for name, site in sites.items():
        site = site.upper()
        expanded[name] = site
        rc = _mut_revcomp(site)
        if rc != site:
            expanded[f"{name}_rc"] = rc
    aa_codons, _ = _codon_build_aa_map(raw)
    dna_list = list(dna)
    fixes: list[str] = []
    for enzyme, site in expanded.items():
        pos = 0
        while True:
            seq = "".join(dna_list)
            idx = seq.find(site, pos)
            if idx == -1:
                break
            fixed = False
            lo_codon = max(0, (idx // 3) - 1)
            hi_codon = (idx + len(site)) // 3 + 2
            for codon_idx in range(lo_codon, hi_codon):
                codon_start = codon_idx * 3
                if codon_start + 3 > len(dna_list) - 3:  # skip stop codon
                    break
                if codon_idx >= len(protein):
                    break
                aa = protein[codon_idx].upper()
                current = "".join(dna_list[codon_start:codon_start + 3])
                for alt, frac in aa_codons.get(aa, []):
                    if alt == current:
                        continue
                    test = dna_list[:]
                    test[codon_start:codon_start + 3] = list(alt)
                    if site not in "".join(test):
                        dna_list = test
                        strand = " (rc)" if enzyme.endswith("_rc") else ""
                        fixes.append(
                            f"{enzyme.replace('_rc', '')}{strand} at nt {idx+1}: "
                            f"{current}→{alt} (codon {codon_idx+1} {aa}, "
                            f"freq={frac:.3f})"
                        )
                        fixed = True
                        break
                if fixed:
                    break
            if not fixed:
                pos = idx + 1
    return "".join(dna_list), fixes


def _codon_cai(dna: str, raw: dict) -> float:
    """Codon Adaptation Index (geometric mean of per-codon freq ÷ peak freq
    of its amino-acid synonymy group). Skips stops and unknown codons."""
    import math
    aa_codons, codon_frac = _codon_build_aa_map(raw)
    w: list[float] = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        entry = raw.get(codon)
        if not entry or entry[0] == "*":
            continue
        peak = aa_codons[entry[0]][0][1] if entry[0] in aa_codons else 0.0
        if peak > 0:
            w.append(codon_frac.get(codon, 0.0) / peak)
    if not w:
        return 0.0
    return math.exp(sum(math.log(max(v, 1e-10)) for v in w) / len(w))


def _codon_gc(dna: str) -> float:
    """GC%. Empty string → 0."""
    if not dna:
        return 0.0
    gc = sum(1 for c in dna.upper() if c in "GC")
    return gc / len(dna) * 100.0


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
    """Design diagnostic PCR primers WITHIN a selected region using Primer3.

    Both primers bind INSIDE the region (target_start..target_end) and the
    amplicon is product_min..product_max bp. This is the standard approach
    for detection/screening primers: you pick a gene or feature and want a
    ~500 bp diagnostic band from within it.

    Uses SEQUENCE_INCLUDED_REGION (not SEQUENCE_TARGET) so Primer3 places
    both primers inside the selected region rather than trying to flank it.

    Returns a dict with keys: fwd_seq, rev_seq, fwd_tm, rev_tm, fwd_pos,
    rev_pos, product_size, or an 'error' key on failure.
    """
    import primer3
    seq   = template_seq.upper()
    total = len(seq)
    wraps = target_end < target_start

    # Primer3 is linear-only. For a wrap region we rotate the template
    # so the region becomes contiguous at [0, region_len), run Primer3,
    # then unrotate the returned positions via (coord + rotation) % total.
    if wraps:
        rotation    = target_start
        p3_seq      = seq[target_start:] + seq[:target_start]
        region_len  = (total - target_start) + target_end
        p3_start    = 0
    else:
        rotation    = 0
        p3_seq      = seq
        region_len  = target_end - target_start
        p3_start    = target_start

    if region_len < 1:
        return {"error": "Target region is empty."}
    if region_len < product_min:
        return {
            "error": f"Region ({region_len} bp) is shorter than minimum "
                     f"product size ({product_min} bp). Select a larger "
                     f"region or reduce the product size."
        }

    try:
        result = primer3.design_primers(
            seq_args={
                "SEQUENCE_TEMPLATE": p3_seq,
                # INCLUDED_REGION: primers must bind WITHIN this region.
                # This is the key difference from SEQUENCE_TARGET (which
                # would require primers to sit OUTSIDE the target).
                "SEQUENCE_INCLUDED_REGION": [p3_start, region_len],
            },
            global_args={
                "PRIMER_TASK": "generic",
                "PRIMER_PICK_LEFT_PRIMER": 1,
                "PRIMER_PICK_RIGHT_PRIMER": 1,
                # primer_len is the OPTIMAL length — Primer3 will expand
                # or contract within the min/max range to find the best Tm.
                "PRIMER_OPT_SIZE": primer_len,
                "PRIMER_MIN_SIZE": max(15, primer_len - 8),
                "PRIMER_MAX_SIZE": min(36, primer_len + 8),
                "PRIMER_OPT_TM": target_tm,
                "PRIMER_MIN_TM": target_tm - 3,
                "PRIMER_MAX_TM": target_tm + 3,
                "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
                "PRIMER_NUM_RETURN": 1,
            },
        )
    except (OSError, Exception) as exc:
        return {"error": f"Primer3 rejected parameters: {exc}"}

    n_found = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    if n_found == 0:
        explain = result.get("PRIMER_LEFT_EXPLAIN", "")
        return {"error": f"Primer3 found no valid pair. {explain}"}

    fwd_pos = result["PRIMER_LEFT_0"]     # (start, length) on p3_seq
    rev_pos = result["PRIMER_RIGHT_0"]    # (start, length) — start is 3' end on p3_seq

    # Unrotate positions back to original-template coordinates.
    fwd_start = (fwd_pos[0] + rotation) % total
    fwd_end   = (fwd_pos[0] + fwd_pos[1] + rotation) % total
    rev_start = (rev_pos[0] - rev_pos[1] + 1 + rotation) % total
    rev_end   = (rev_pos[0] + 1 + rotation) % total

    return {
        "fwd_seq":      result["PRIMER_LEFT_0_SEQUENCE"],
        "rev_seq":      result["PRIMER_RIGHT_0_SEQUENCE"],
        "fwd_tm":       round(result["PRIMER_LEFT_0_TM"], 1),
        "rev_tm":       round(result["PRIMER_RIGHT_0_TM"], 1),
        "fwd_pos":      (fwd_start, fwd_end),
        "rev_pos":      (rev_start, rev_end),
        "product_size": result["PRIMER_PAIR_0_PRODUCT_SIZE"],
    }


def _design_cloning_primers_raw(
    template_seq: str,
    start: int,
    end: int,
    site_5: str,
    site_3: str,
    name_5: str = "5'site",
    name_3: str = "3'site",
    target_tm: float = 60.0,
    padding: str = "GCGC",
) -> dict:
    """Design cloning primers with arbitrary recognition-site tails + padding.

    Accepts raw site sequences (not just NEB enzyme names) so users can
    enter custom cutter sequences.

    Structure (5'→3'):
        Forward: [padding] [5' site]    [binding region →]
        Reverse: [padding] [RC 3' site] [← binding region RC]

    Returns dict with keys: fwd_full, rev_full, fwd_binding, rev_binding,
    fwd_tm, rev_tm, re_5prime, re_3prime, site_5, site_3, insert_seq,
    fwd_pos, rev_pos, or 'error'.
    """
    site_5 = site_5.upper()
    site_3 = site_3.upper()
    if not site_5 or not set(site_5) <= set("ACGTRYWSMKBDHVN"):
        return {"error": f"Invalid 5' site sequence: {site_5!r}"}
    if not site_3 or not set(site_3) <= set("ACGTRYWSMKBDHVN"):
        return {"error": f"Invalid 3' site sequence: {site_3!r}"}

    total  = len(template_seq)
    insert = _slice_circular(template_seq.upper(), start, end)
    wraps  = end < start
    if len(insert) < 18:
        return {"error": "Region too short (< 18 bp)."}

    fwd_bind, fwd_tm = _pick_binding_region(insert, target_tm)
    rev_bind, rev_tm = _pick_binding_region(_rc(insert), target_tm)

    fwd_full = padding + site_5 + fwd_bind
    rev_full = padding + _rc(site_3) + rev_bind

    if wraps:
        fwd_pos = (start, (start + len(fwd_bind)) % total)
        rev_pos = ((end - len(rev_bind)) % total, end)
    else:
        fwd_pos = (start, start + len(fwd_bind))
        rev_pos = (end - len(rev_bind), end)

    return {
        "fwd_full":    fwd_full,
        "rev_full":    rev_full,
        "fwd_binding": fwd_bind,
        "rev_binding": rev_bind,
        "fwd_tm":      round(fwd_tm, 1),
        "rev_tm":      round(rev_tm, 1),
        "re_5prime":   name_5,
        "re_3prime":   name_3,
        "site_5":      site_5,
        "site_3":      site_3,
        "insert_seq":  insert,
        "fwd_pos":     fwd_pos,
        "rev_pos":     rev_pos,
    }


def _design_cloning_primers(
    template_seq: str,
    start: int,
    end: int,
    re_5prime: str,
    re_3prime: str,
    target_tm: float = 60.0,
    padding: str = "GCGC",
) -> dict:
    """Design cloning primers using NEB enzyme names. Delegates to
    _design_cloning_primers_raw after looking up recognition sites."""
    if re_5prime not in _NEB_ENZYMES:
        return {"error": f"Unknown enzyme: {re_5prime}"}
    if re_3prime not in _NEB_ENZYMES:
        return {"error": f"Unknown enzyme: {re_3prime}"}
    site_5, _, _ = _NEB_ENZYMES[re_5prime]
    site_3, _, _ = _NEB_ENZYMES[re_3prime]
    return _design_cloning_primers_raw(
        template_seq, start, end, site_5, site_3,
        name_5=re_5prime, name_3=re_3prime,
        target_tm=target_tm, padding=padding,
    )


def _design_generic_primers(
    template_seq: str,
    start: int,
    end: int,
    target_tm: float = 60.0,
) -> dict:
    """Design simple binding primers (no tails, no RE sites, no overhangs).

    Forward primer: optimal binding region at the start of the region.
    Reverse primer: optimal binding region at the end (reverse-complement).
    """
    total  = len(template_seq)
    insert = _slice_circular(template_seq.upper(), start, end)
    wraps  = end < start
    if len(insert) < 18:
        return {"error": "Region too short (< 18 bp)."}
    fwd_bind, fwd_tm = _pick_binding_region(insert, target_tm)
    rev_bind, rev_tm = _pick_binding_region(_rc(insert), target_tm)
    if wraps:
        fwd_pos = (start, (start + len(fwd_bind)) % total)
        rev_pos = ((end - len(rev_bind)) % total, end)
    else:
        fwd_pos = (start, start + len(fwd_bind))
        rev_pos = (end - len(rev_bind), end)
    return {
        "fwd_seq":  fwd_bind,
        "rev_seq":  rev_bind,
        "fwd_tm":   round(fwd_tm, 1),
        "rev_tm":   round(rev_tm, 1),
        "fwd_pos":  fwd_pos,
        "rev_pos":  rev_pos,
    }


# ── SOE-PCR mutagenesis primer design ──────────────────────────────────────────
#
# Ported from superfolder_aeblue/mutagenesis_primers.py. Designs a 4-primer /
# 3-reaction SOE (Splicing by Overlap Extension) workflow to introduce a
# single-residue point mutation into a CDS and produce a Golden Braid
# B3-B5-compatible amplicon in one pot.
#
# OUTER primers (constant per CDS):
#   FWD:  CCCC-GGTCTC-A-AATG-[CDS body from codon 2] → B3 overhang AATG
#   REV:  CCCC-GGTCTC-A-AACG-[revcomp of CDS end]    → B5 overhang CGTT
# INNER primers (one pair per mutation): mutant codon centered in anneal region,
# REV = revcomp(FWD). Edge cases (mutation < 60 nt from either CDS end) swap the
# inner pair for a single "modified outer" primer and a 2-primer direct PCR.
#
# Codon table is E. coli K12 (Kazusa taxid 83333). This is a single-organism
# convenience; callers designing for other hosts should review mut_codon choice.

_MUT_CODON_USAGE = {
    "GGG": ("G", 44),  "GGA": ("G", 47),  "GGT": ("G", 109), "GGC": ("G", 171),
    "GAG": ("E", 94),  "GAA": ("E", 224), "GAT": ("D", 194), "GAC": ("D", 105),
    "GTG": ("V", 135), "GTA": ("V", 59),  "GTT": ("V", 86),  "GTC": ("V", 60),
    "GCG": ("A", 197), "GCA": ("A", 108), "GCT": ("A", 55),  "GCC": ("A", 162),
    "AGG": ("R", 8),   "AGA": ("R", 7),   "AGT": ("S", 37),  "AGC": ("S", 85),
    "AAG": ("K", 62),  "AAA": ("K", 170), "AAT": ("N", 112), "AAC": ("N", 125),
    "ATG": ("M", 127), "ATA": ("I", 19),  "ATT": ("I", 156), "ATC": ("I", 93),
    "ACG": ("T", 59),  "ACA": ("T", 33),  "ACT": ("T", 41),  "ACC": ("T", 117),
    "TGG": ("W", 55),  "TGT": ("C", 30),  "TGC": ("C", 41),
    "TAT": ("Y", 86),  "TAC": ("Y", 75),
    "TTG": ("L", 61),  "TTA": ("L", 78),  "TTT": ("F", 101), "TTC": ("F", 77),
    "TCG": ("S", 41),  "TCA": ("S", 40),  "TCT": ("S", 29),  "TCC": ("S", 28),
    "CGG": ("R", 21),  "CGA": ("R", 22),  "CGT": ("R", 108), "CGC": ("R", 133),
    "CAG": ("Q", 142), "CAA": ("Q", 62),  "CAT": ("H", 81),  "CAC": ("H", 67),
    "CTG": ("L", 240), "CTA": ("L", 27),  "CTT": ("L", 61),  "CTC": ("L", 54),
    "CCG": ("P", 137), "CCA": ("P", 34),  "CCT": ("P", 43),  "CCC": ("P", 33),
    "TAA": ("*", 9),   "TAG": ("*", 0),   "TGA": ("*", 5),
}
_MUT_CODON_TO_AA = {c: aa for c, (aa, _) in _MUT_CODON_USAGE.items()}
_MUT_STOPS       = {"TAA", "TAG", "TGA"}

def _mut_aa_to_codons() -> dict:
    from collections import defaultdict
    totals = defaultdict(int)
    for c, (aa, n) in _MUT_CODON_USAGE.items():
        totals[aa] += n
    result: dict = defaultdict(list)
    for c, (aa, n) in _MUT_CODON_USAGE.items():
        if aa == "*":
            continue
        result[aa].append((c, n / totals[aa] if totals[aa] else 0.0))
    for aa in result:
        result[aa].sort(key=lambda x: -x[1])
    return dict(result)

_MUT_AA_TO_CODONS = _mut_aa_to_codons()

_MUT_BSAI_FWD_TAIL = "CCCC" + "GGTCTCA" + "AATG"   # 15 nt; AATG = A(extra)+ATG ovhg
_MUT_BSAI_REV_TAIL = "CCCC" + "GGTCTCA" + "AACG"   # 15 nt; AACG = revcomp(CGTT)
_MUT_MIN_SOE_FRAG  = 60                             # nt; below this → edge case

_MUT_P3 = dict(mv_conc=50.0, dv_conc=1.5, dntp_conc=0.2, dna_conc=250.0)


def _mut_parse(s: str) -> tuple:
    """Parse a mutation string like 'W140F'. Returns (wt_aa, pos_1based, mut_aa)."""
    m = re.fullmatch(r"([A-Za-z\*])(\d+)([A-Za-z\*])", s.strip())
    if not m:
        raise ValueError(f"Cannot parse '{s}'. Use format: [WT][pos][MUT], e.g. W140F")
    return m.group(1).upper(), int(m.group(2)), m.group(3).upper()


def _mut_revcomp(seq: str) -> str:
    return seq.upper().translate(_IUPAC_COMP)[::-1]


def _mut_translate(dna: str) -> str:
    aa: list = []
    for i in range(0, len(dna) - 2, 3):
        c = dna[i:i+3].upper()
        if c in _MUT_STOPS:
            break
        aa.append(_MUT_CODON_TO_AA.get(c, "?"))
    return "".join(aa)


_MUT_PREVIEW_DNA_COLOR = "color(118)"   # CDS-green, main-app-like
_MUT_PREVIEW_AA_COLOR  = "color(141)"   # purple
_MUT_PREVIEW_MUT_COLOR = "color(208)"   # orange (mutated codon + AA)


def _mut_build_preview_text(cds_dna: str,
                            protein_override: str = "",
                            mutation: "dict | None" = None,
                            line_width: int = 90,
                            cursor_aa: int = -1,
                            dna_color: str = _MUT_PREVIEW_DNA_COLOR,
                            aa_color:  str = _MUT_PREVIEW_AA_COLOR,
                            mut_color: str = _MUT_PREVIEW_MUT_COLOR) -> Text:
    """Build a Rich Text preview of the (optionally mutagenized) CDS.

    Two modes:

    * **DNA + AA**: `cds_dna` is non-empty. DNA is rendered main-app-style
      (green by default) with line numbers. Beneath each DNA row sits an
      AA row — one letter per codon, centered under the middle base of
      that codon, in purple. Stop codons show as `*` (not truncated,
      unlike `_mut_translate`, so the user can see the end of the CDS).
    * **AA only**: `cds_dna` is empty but `protein_override` is set. Used
      while the user is typing a protein in the custom-protein source
      before harmonization produces any DNA.

    `mutation` (dict with `wt_codon` / `mut_codon` / `nt_position` 1-based,
    as produced by `_mut_design_inner`) substitutes the mutant codon into
    the displayed DNA and highlights the three mutated bases plus the AA
    letter below in solid orange.

    `cursor_aa` (>= 0) marks the AA the user is focused on — it gets a
    reverse-video highlight across both DNA (the whole codon) and AA
    rows. Reverse style stacks with mutation: cursor-on-mutation ends up
    as reversed orange, which still reads as "cursor here, mutant here".
    """
    t = Text(no_wrap=True, overflow="crop")
    if not cds_dna:
        aa = (protein_override or "").upper()
        if not aa:
            return t
        for i in range(0, len(aa), line_width):
            chunk = aa[i:i + line_width]
            for j, ch in enumerate(chunk):
                idx = i + j
                if idx == cursor_aa:
                    t.append(ch, style=f"bold reverse {aa_color}")
                else:
                    t.append(ch, style=f"bold {aa_color}")
            t.append("\n")
        return t

    # Line width must be a multiple of 3 so codons don't straddle lines.
    lw = max(3, (line_width // 3) * 3)

    dna = cds_dna.upper()
    mut_lo = mut_hi = -1
    if mutation:
        wt_c  = (mutation.get("wt_codon")  or "").upper()
        mut_c = (mutation.get("mut_codon") or "").upper()
        try:
            nt_pos = int(mutation.get("nt_position") or 0)
        except (TypeError, ValueError):
            nt_pos = 0
        if wt_c and mut_c and 1 <= nt_pos <= len(dna) - 2:
            lo = nt_pos - 1
            dna = dna[:lo] + mut_c + dna[lo + 3:]
            mut_lo, mut_hi = lo, lo + 3

    n     = len(dna)
    num_w = len(str(n)) if n else 1
    cur_dna_lo = cursor_aa * 3 if cursor_aa >= 0 else -1
    cur_dna_hi = cur_dna_lo + 3 if cursor_aa >= 0 else -1

    def _base_style(i: int) -> str:
        is_mut = mut_lo <= i < mut_hi
        is_cur = cur_dna_lo <= i < cur_dna_hi
        if is_mut and is_cur:
            return f"bold reverse {mut_color}"
        if is_mut:
            return f"bold {mut_color}"
        if is_cur:
            return f"reverse {dna_color}"
        return dna_color

    def _aa_style(aa_idx: int) -> str:
        mid_i  = aa_idx * 3 + 1
        is_mut = mut_lo <= mid_i < mut_hi
        is_cur = (aa_idx == cursor_aa)
        if is_mut and is_cur:
            return f"bold reverse {mut_color}"
        if is_mut:
            return f"bold {mut_color}"
        if is_cur:
            return f"bold reverse {aa_color}"
        return f"bold {aa_color}"

    for chunk_start in range(0, n, lw):
        chunk_end = min(chunk_start + lw, n)
        t.append(f"{chunk_start + 1:>{num_w}d}  ", style="dim")
        # DNA row
        for i in range(chunk_start, chunk_end):
            t.append(dna[i], style=_base_style(i))
        t.append("\n")
        # AA row — one letter centered under each codon's middle base
        t.append(" " * (num_w + 2))
        for i in range(chunk_start, chunk_end):
            if i % 3 == 1:
                aa_idx = i // 3
                codon  = dna[aa_idx * 3:aa_idx * 3 + 3]
                aa_ch  = _MUT_CODON_TO_AA.get(codon, "?")
                t.append(aa_ch, style=_aa_style(aa_idx))
            else:
                t.append(" ")
        t.append("\n")

    return t


def _mut_next_cursor(current: int, protein_len: int, line_width: int,
                     dna_mode: bool, direction: str) -> int:
    """Compute the next cursor position given an arrow-key direction.

    * `current` — current cursor AA index; -1 means no cursor yet.
    * `direction` ∈ {"left", "right", "up", "down"}.
    * Up/Down step by one row's worth of amino acids (`line_width // 3`
      in DNA mode, `line_width` in AA-only mode).
    * Result is clamped to `[0, protein_len)`. -1 is returned for empty
      proteins. First press after no-cursor snaps to index 0 regardless
      of direction so arrow keys "wake up" the cursor intuitively.
    """
    if protein_len <= 0:
        return -1
    if current < 0:
        return 0
    step = max(1, (line_width // 3) if dna_mode else line_width)
    if   direction == "left":  new_idx = current - 1
    elif direction == "right": new_idx = current + 1
    elif direction == "up":    new_idx = current - step
    elif direction == "down":  new_idx = current + step
    else:                      return current
    return max(0, min(protein_len - 1, new_idx))


def _mut_click_to_aa_index(dna_mode: bool, dna_len: int, protein_len: int,
                           line_width: int, pad: int,
                           vp_x: int, content_row: int) -> int:
    """Translate a click at (vp_x, content_row) inside the preview widget
    to an amino-acid index, or -1 if the click missed an AA.

    Pure arithmetic — factored out so it can be unit-tested without
    standing up a Textual event loop. `content_row` is the click's
    row relative to the *content* (viewport row + scroll offset), not
    the raw event.y. `vp_x` is the column relative to the widget.
    """
    if protein_len <= 0 or line_width <= 0:
        return -1
    if dna_mode:
        # Two rendered rows per logical codon-line (DNA + AA)
        logical_line = content_row // 2
        bp_start = logical_line * line_width
        if bp_start < 0 or bp_start >= dna_len:
            return -1
        c_data = vp_x - pad
        if c_data < 0 or c_data >= line_width:
            return -1
        codon_idx_in_line = c_data // 3
        aa_idx = bp_start // 3 + codon_idx_in_line
    else:
        if content_row < 0:
            return -1
        aa_idx = content_row * line_width + (vp_x - pad)
    if aa_idx < 0 or aa_idx >= protein_len:
        return -1
    return aa_idx


def _mut_tm(seq: str) -> float:
    try:
        import primer3
        return primer3.calc_tm(seq, **_MUT_P3)
    except Exception:
        gc = sum(1 for c in seq.upper() if c in "GC")
        at = sum(1 for c in seq.upper() if c in "AT")
        return 2 * at + 4 * gc


def _mut_hairpin_dg(seq: str) -> float:
    try:
        import primer3
        return primer3.calc_hairpin(seq, **_MUT_P3).dg
    except Exception:
        return 0.0


def _mut_homodimer_dg(seq: str) -> float:
    try:
        import primer3
        return primer3.calc_homodimer(seq, **_MUT_P3).dg
    except Exception:
        return 0.0


def _mut_gc_pct(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) * 100 if seq else 0.0


def _mut_ends_gc(seq: str) -> bool:
    return bool(seq) and seq[-1].upper() in "GC"


def _mut_score_outer(anneal: str, target_tm: float = 60.0) -> float:
    t  = _mut_tm(anneal)
    gc = _mut_gc_pct(anneal)
    hp = _mut_hairpin_dg(anneal)
    return (
        abs(t - target_tm) * 2.0
        + (0 if _mut_ends_gc(anneal) else 4.0)
        + max(0, -hp - 1000) / 400.0
        + abs(gc - 50) * 0.1
    )


def _mut_design_fwd_anneal(dna: str) -> "dict | None":
    body = dna[3:]
    best = None
    for length in range(18, 28):
        anneal = body[:length]
        if len(anneal) < 18:
            continue
        s = _mut_score_outer(anneal)
        if best is None or s < best["score"]:
            best = {
                "anneal": anneal,
                "full":   _MUT_BSAI_FWD_TAIL + anneal,
                "tm_anneal": _mut_tm(anneal),
                "gc":     _mut_gc_pct(anneal),
                "score":  s,
            }
    return best


def _mut_design_rev_anneal(dna: str) -> "dict | None":
    end_rc = _mut_revcomp(dna)
    best = None
    for length in range(18, 28):
        anneal = end_rc[:length]
        if len(anneal) < 18:
            continue
        s = _mut_score_outer(anneal)
        if best is None or s < best["score"]:
            best = {
                "anneal": anneal,
                "full":   _MUT_BSAI_REV_TAIL + anneal,
                "tm_anneal": _mut_tm(anneal),
                "gc":     _mut_gc_pct(anneal),
                "score":  s,
            }
    return best


def _mut_design_outer(dna: str) -> dict:
    """Constant FWD/REV outer primers with BsaI-AATG / BsaI-AACG tails."""
    fwd = _mut_design_fwd_anneal(dna)
    rev = _mut_design_rev_anneal(dna)
    if fwd is None or rev is None:
        raise RuntimeError("CDS is too short to design outer primers (need ≥ 21 nt).")
    return {
        "fwd": fwd, "rev": rev,
        "b3_overhang": "AATG",
        "b5_overhang": "CGTT",
        "fwd_anneal_start": 3,
    }


def _mut_design_modified_outer(dna_mut: str, near_start: bool) -> dict:
    """Edge-case: mutation < 60 nt from a CDS end → fold mutant codon into a
    single outer primer. PCR becomes a 2-primer direct reaction, no SOE."""
    if near_start:
        p = _mut_design_fwd_anneal(dna_mut)
        if p is None:
            raise RuntimeError("Modified FWD outer design failed.")
        p["label"]    = "modified_FWD_outer"
        p["partner"]  = "REV_outer (unchanged)"
        p["replaces"] = "FWD_outer"
    else:
        p = _mut_design_rev_anneal(dna_mut)
        if p is None:
            raise RuntimeError("Modified REV outer design failed.")
        p["label"]    = "modified_REV_outer"
        p["partner"]  = "FWD_outer (unchanged)"
        p["replaces"] = "REV_outer"
    return p


def _mut_design_inner(dna: str, mut_pos_1: int, mut_aa: str, wt_aa: str,
                      codon_table: "dict | None" = None) -> dict:
    """Inner mutagenic pair (FWD carries mutant codon; REV = revcomp(FWD)).

    `codon_table` is an optional {codon: (aa, count)} map used to pick the
    mutant codon. Defaults to E. coli K12 (_MUT_AA_TO_CODONS)."""
    idx      = mut_pos_1 - 1
    nt_start = idx * 3

    wt_codon  = dna[nt_start:nt_start + 3]
    if len(wt_codon) < 3:
        raise ValueError(
            f"Position {mut_pos_1} is past the end of the CDS."
        )
    wt_actual = _MUT_CODON_TO_AA.get(wt_codon, "?")
    if wt_actual != wt_aa:
        raise ValueError(
            f"Position {mut_pos_1}: mutation says WT='{wt_aa}' but DNA codon "
            f"'{wt_codon}' encodes '{wt_actual}'."
        )

    if codon_table:
        aa_map, _ = _codon_build_aa_map(codon_table)
    else:
        aa_map = _MUT_AA_TO_CODONS

    if mut_aa == "*":
        mut_codon = "TAA"
    else:
        mut_codon = next(
            (c for c, _f in aa_map.get(mut_aa, []) if c != wt_codon),
            None,
        )
        if mut_codon is None:
            raise ValueError(f"No alternative codon available for '{mut_aa}' "
                             "in the selected codon table")

    mut_dna = dna[:nt_start] + mut_codon + dna[nt_start + 3:]

    TM_TARGET      = 60.0
    TM_MIN, TM_MAX = 55.0, 75.0
    GC_MIN, GC_MAX = 35.0, 68.0
    seq_len = len(mut_dna)

    candidates: list = []
    for left_ext in range(5, 28):
        for right_ext in range(5, 28):
            lo  = max(0, nt_start - left_ext)
            hi  = min(seq_len, nt_start + 3 + right_ext)
            fwd = mut_dna[lo:hi]
            if len(fwd) < 15 or len(fwd) > 58:
                continue
            t  = _mut_tm(fwd)
            gc = _mut_gc_pct(fwd)
            if not (TM_MIN <= t <= TM_MAX):
                continue
            if not (GC_MIN <= gc <= GC_MAX):
                continue
            hp = _mut_hairpin_dg(fwd)
            hd = _mut_homodimer_dg(fwd)
            score = (
                abs(t - TM_TARGET) * 2.0
                + (0 if _mut_ends_gc(fwd) else 4.0)
                + max(0, -hp - 1000) / 400.0
                + max(0, -hd - 2000) / 400.0
                + abs(gc - 50) * 0.1
                - (len(fwd) * 0.15 if abs(t - TM_TARGET) <= 1.0 else 0)
            )
            candidates.append({
                "fwd": fwd, "rev": _mut_revcomp(fwd),
                "tm": t, "gc": gc, "length": len(fwd),
                "hairpin_dg": hp, "homodimer_dg": hd, "score": score, "lo": lo,
            })

    if not candidates:
        raise RuntimeError(
            f"No valid inner primers found for {wt_aa}{mut_pos_1}{mut_aa}. "
            "Mutation may be too close to sequence ends."
        )

    seen: dict = {}
    for c in sorted(candidates, key=lambda x: x["score"]):
        if c["fwd"] not in seen:
            seen[c["fwd"]] = c
    ranked = sorted(seen.values(), key=lambda x: x["score"])[:5]
    for i, c in enumerate(ranked):
        c["rank"] = i + 1

    best_lo = ranked[0]["lo"]
    best_hi = best_lo + ranked[0]["length"]
    fwd_anneal_start = 3
    frag_a = best_hi - fwd_anneal_start
    frag_b = seq_len - best_lo

    near_start = frag_a < _MUT_MIN_SOE_FRAG
    near_end   = frag_b < _MUT_MIN_SOE_FRAG

    edge_case = None
    if near_start or near_end:
        modified_outer = _mut_design_modified_outer(mut_dna, near_start=near_start)
        edge_case = {
            "near_start":     near_start,
            "near_end":       near_end,
            "frag_a":         frag_a,
            "frag_b":         frag_b,
            "modified_outer": modified_outer,
        }

    return {
        "mutation":    f"{wt_aa}{mut_pos_1}{mut_aa}",
        "nt_position": nt_start + 1,
        "wt_codon":    wt_codon,
        "mut_codon":   mut_codon,
        "nt_changes":  sum(a != b for a, b in zip(wt_codon, mut_codon)),
        "candidates":  ranked,
        "edge_case":   edge_case,
    }


def _mut_extract_cds(full_seq: str, start: int, end: int, strand: int) -> str:
    """Return the CDS DNA in its biological 5'→3' orientation, handling
    origin-wrap (end < start) and reverse-strand features."""
    if end < start:
        sub = full_seq[start:] + full_seq[:end]
    else:
        sub = full_seq[start:end]
    sub = sub.upper()
    if strand == -1:
        sub = _mut_revcomp(sub)
    return sub


# ── Feature picker: pick one feature from a chosen plasmid ────────────────────

class PlasmidFeaturePickerModal(ModalScreen):
    """Scrollable list of non-source features from a specific library entry.

    Dismisses with a feature-library-style entry dict
    ``{name, feature_type, sequence, strand, qualifiers, description}``,
    or None on cancel. No persistence side effects — the caller decides
    whether to save the picked entry or just use it to prefill a form.
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, entries: list[dict], plasmid_name: str = ""):
        super().__init__()
        self._entries = list(entries)
        self._plasmid_name = plasmid_name or "plasmid"

    def compose(self) -> ComposeResult:
        with Vertical(id="featpick-dlg"):
            yield Static(f" Feature from [{self._plasmid_name}] ",
                         id="featpick-title")
            yield DataTable(id="featpick-table", cursor_type="row",
                            zebra_stripes=True)
            with Horizontal(id="featpick-btns"):
                yield Button("Select", id="btn-featpick-ok", variant="primary")
                yield Button("Cancel", id="btn-featpick-cancel")

    def on_mount(self) -> None:
        t = self.query_one("#featpick-table", DataTable)
        t.add_columns("Name", "Type", "Strand", "Length")
        for i, e in enumerate(self._entries):
            strand_str = "+" if e.get("strand", 1) == 1 else "−"
            t.add_row(
                e.get("name", "?"),
                e.get("feature_type", "?"),
                strand_str,
                f"{len(e.get('sequence', ''))} bp",
                key=str(i),
            )
        if self._entries:
            t.move_cursor(row=0)
            t.focus()

    @on(Button.Pressed, "#btn-featpick-ok")
    def _select(self, _):
        self._dismiss_cursor()

    @on(DataTable.RowSelected, "#featpick-table")
    def _row_selected(self, event):
        if event.row_key and event.row_key.value is not None:
            try:
                idx = int(event.row_key.value)
            except (TypeError, ValueError):
                return
            if 0 <= idx < len(self._entries):
                self.dismiss(dict(self._entries[idx]))

    def _dismiss_cursor(self) -> None:
        t = self.query_one("#featpick-table", DataTable)
        if t.row_count == 0:
            self.dismiss(None)
            return
        row_keys = list(t.rows.keys())
        if 0 <= t.cursor_row < len(row_keys):
            key = row_keys[t.cursor_row].value
            try:
                idx = int(key)
            except (TypeError, ValueError):
                self.dismiss(None)
                return
            if 0 <= idx < len(self._entries):
                self.dismiss(dict(self._entries[idx]))
                return
        self.dismiss(None)

    @on(Button.Pressed, "#btn-featpick-cancel")
    def _cancel_btn(self, _):
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Add-feature modal ──────────────────────────────────────────────────────────

class AddFeatureModal(ModalScreen):
    """Create or edit a feature-library entry.

    The modal collects: name, feature_type (from `_GENBANK_FEATURE_TYPES`),
    strand, 5'→3' sequence, and a qualifier line (``key=value; key=value``).
    Three terminal actions:

      * **Save to Library** — dismisses with ``{"action": "save", "entry": {...}}``.
      * **Insert at cursor** — dismisses with ``{"action": "insert", "entry": {...}}``.
      * **Cancel** — dismisses with ``None``.

    An **Import from plasmid** button opens the PlasmidPickerModal
    → PlasmidFeaturePickerModal chain and prefills the form.

    Validation (non-empty name, ACGT/IUPAC-only sequence, known feature type)
    lives on the app side; the modal just surfaces the error strings.
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, prefill: "dict | None" = None,
                 have_cursor: bool = False) -> None:
        super().__init__()
        self._prefill = dict(prefill) if prefill else {}
        self._have_cursor = have_cursor

    def compose(self) -> ComposeResult:
        p = self._prefill
        name        = p.get("name", "")
        feat_type   = p.get("feature_type", "CDS")
        sequence    = p.get("sequence", "")
        strand      = p.get("strand", 1)
        quals_str   = _qualifiers_to_string(p.get("qualifiers") or {})
        description = p.get("description", "")

        type_options = [(t, t) for t in _GENBANK_FEATURE_TYPES]
        if feat_type not in _GENBANK_FEATURE_TYPES:
            type_options.insert(0, (f"{feat_type} (non-standard)", feat_type))

        with Vertical(id="addfeat-dlg"):
            yield Static(" Add Feature ", id="addfeat-title")
            with ScrollableContainer(id="addfeat-body"):
                yield Label("Name:")
                yield Input(value=name, placeholder="e.g. lacZ-alpha",
                            id="addfeat-name")

                with Horizontal(id="addfeat-row1"):
                    with Vertical(id="addfeat-type-col"):
                        yield Label("Feature type:")
                        yield Select(type_options, value=feat_type,
                                     id="addfeat-type", allow_blank=False)
                    with Vertical(id="addfeat-strand-col"):
                        yield Label("Strand:")
                        with RadioSet(id="addfeat-strand"):
                            yield RadioButton("Forward (+)",
                                              value=(strand != -1),
                                              id="addfeat-strand-fwd")
                            yield RadioButton("Reverse (−)",
                                              value=(strand == -1),
                                              id="addfeat-strand-rev")

                yield Label("Sequence  (5'→3', ACGT/IUPAC; whitespace ignored):")
                yield TextArea(sequence, id="addfeat-seq")
                yield Label("Qualifiers  (e.g.  gene=lacZ; product=LacZ alpha):")
                yield Input(value=quals_str,
                            placeholder="key=value; key=value",
                            id="addfeat-quals")
                yield Label("Description  (optional):")
                yield Input(value=description, placeholder="free text",
                            id="addfeat-desc")
            yield Static("", id="addfeat-status", markup=True)
            with Horizontal(id="addfeat-btns"):
                yield Button("Import from plasmid…",
                             id="btn-addfeat-import")
                yield Button("Save to Library",
                             id="btn-addfeat-save",
                             variant="primary")
                yield Button("Insert at cursor",
                             id="btn-addfeat-insert",
                             variant="success",
                             disabled=not self._have_cursor)
                yield Button("Cancel", id="btn-addfeat-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#addfeat-name", Input).focus()
        except NoMatches:
            pass

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _gather(self) -> "dict | None":
        """Read form → entry dict, or write a red error and return None."""
        try:
            name   = self.query_one("#addfeat-name",  Input).value.strip()
            ftype  = self.query_one("#addfeat-type",  Select).value
            seqraw = self.query_one("#addfeat-seq",   TextArea).text
            quals  = self.query_one("#addfeat-quals", Input).value
            desc   = self.query_one("#addfeat-desc",  Input).value.strip()
            fwd_rb = self.query_one("#addfeat-strand-fwd", RadioButton)
        except NoMatches:
            return None
        status = self.query_one("#addfeat-status", Static)

        if not name:
            status.update("[red]Name cannot be empty.[/red]")
            return None
        # Normalise + validate sequence
        seq_clean = "".join(ch for ch in seqraw.upper() if not ch.isspace())
        if not seq_clean:
            status.update("[red]Sequence cannot be empty.[/red]")
            return None
        allowed = set("ACGTRYWSMKBDHVN")
        bad = [ch for ch in seq_clean if ch not in allowed]
        if bad:
            status.update(
                f"[red]Sequence has invalid bases: {''.join(sorted(set(bad)))[:10]}[/red]"
            )
            return None
        if ftype is None or ftype == Select.BLANK:
            status.update("[red]Choose a feature type.[/red]")
            return None
        strand = 1 if fwd_rb.value else -1
        return {
            "name":         name,
            "feature_type": str(ftype),
            "sequence":     seq_clean,
            "strand":       strand,
            "qualifiers":   _parse_qualifier_string(quals),
            "description":  desc,
        }

    # ── Buttons ──────────────────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-addfeat-save")
    def _save(self, _) -> None:
        entry = self._gather()
        if entry is None:
            return
        self.dismiss({"action": "save", "entry": entry})

    @on(Button.Pressed, "#btn-addfeat-insert")
    def _insert(self, _) -> None:
        entry = self._gather()
        if entry is None:
            return
        self.dismiss({"action": "insert", "entry": entry})

    @on(Button.Pressed, "#btn-addfeat-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    @on(Button.Pressed, "#btn-addfeat-import")
    def _import(self, _) -> None:
        entries = _load_library()
        if not entries:
            self.query_one("#addfeat-status", Static).update(
                "[yellow]Library is empty — save a plasmid first.[/yellow]"
            )
            return

        def _on_plasmid(entry_id):
            if not entry_id:
                return
            match = next((e for e in _load_library()
                          if e.get("id") == entry_id), None)
            if not match:
                self.query_one("#addfeat-status", Static).update(
                    "[red]Entry not found.[/red]"
                )
                return
            gb_text = match.get("gb_text", "")
            if not gb_text:
                self.query_one("#addfeat-status", Static).update(
                    "[red]Library entry has no sequence.[/red]"
                )
                return
            try:
                rec = _gb_text_to_record(gb_text)
            except Exception as exc:    # noqa: BLE001
                _log.exception("Feature import: failed to parse library entry")
                self.query_one("#addfeat-status", Static).update(
                    f"[red]Failed to load plasmid: {exc}[/red]"
                )
                return
            feat_entries = _extract_feature_entries_from_record(rec)
            if not feat_entries:
                self.query_one("#addfeat-status", Static).update(
                    "[yellow]Plasmid has no non-source features.[/yellow]"
                )
                return

            def _on_feat(picked):
                if not picked:
                    return
                self._apply_prefill(picked)

            self.app.push_screen(
                PlasmidFeaturePickerModal(feat_entries,
                                          plasmid_name=match.get("name", "")),
                callback=_on_feat,
            )

        self.app.push_screen(PlasmidPickerModal(None), callback=_on_plasmid)

    def _apply_prefill(self, entry: dict) -> None:
        """Fill the form fields from a picked feature. Leaves user's current
        name alone if it's non-empty so accidental imports don't clobber
        typed data."""
        try:
            name_inp = self.query_one("#addfeat-name", Input)
            type_sel = self.query_one("#addfeat-type", Select)
            seq_ta   = self.query_one("#addfeat-seq", TextArea)
            quals_in = self.query_one("#addfeat-quals", Input)
            desc_in  = self.query_one("#addfeat-desc", Input)
            fwd_rb   = self.query_one("#addfeat-strand-fwd", RadioButton)
            rev_rb   = self.query_one("#addfeat-strand-rev", RadioButton)
            status   = self.query_one("#addfeat-status", Static)
        except NoMatches:
            return
        if not name_inp.value.strip():
            name_inp.value = entry.get("name", "")
        ftype = entry.get("feature_type", "misc_feature")
        # Select raises if the value isn't a known option; prepend it if so.
        current_options = [v for _, v in getattr(type_sel, "_options", [])]
        if ftype not in current_options:
            new_opts = [(f"{ftype} (non-standard)", ftype)] + \
                       [(t, t) for t in _GENBANK_FEATURE_TYPES]
            type_sel.set_options(new_opts)
        type_sel.value = ftype
        seq_ta.text = entry.get("sequence", "")
        quals_in.value = _qualifiers_to_string(entry.get("qualifiers") or {})
        desc_in.value = entry.get("description", "") or ""
        is_rev = entry.get("strand", 1) == -1
        fwd_rb.value = not is_rev
        rev_rb.value = is_rev
        status.update(
            f"[green]Imported '{entry.get('name', '?')}' — "
            f"review and Save or Insert.[/green]"
        )

    def action_cancel(self) -> None:
        self.dismiss(None)


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

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(id="parts-box"):
            yield Static(" Parts Bin  —  Golden Braid L0 Parts ", id="parts-title")
            yield DataTable(id="parts-table", cursor_type="row", zebra_stripes=True)
            yield Static("", id="parts-detail")
            with Horizontal(id="parts-btns"):
                yield Button("New Part", id="btn-new-part", variant="primary")
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
        except NoMatches:
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

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

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
                    "Save to Parts Bin", id="btn-dom-save", variant="primary",
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
        total = len(self._template)
        if start < 0 or end < 0 or start >= total or end > total:
            status.update(
                f"[red]Invalid region: {start+1}–{end} "
                f"(plasmid is {total} bp)[/red]"
            )
            return
        if start == end:
            status.update("[red]Region is empty.[/red]")
            return
        region_len = _feat_len(start, end, total)
        if region_len < 20:
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

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

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
                    variant="primary", disabled=True
                )
                yield Button("Clear Lane", id="btn-ctor-clear", variant="default")
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


# ── NCBI taxon picker (sub-modal for species-name → taxid lookup) ─────────────

class NcbiTaxonPickerModal(ModalScreen):
    """Shown when the user types a non-numeric query in the Kazusa fetch
    field. Searches NCBI taxonomy (with an auto-wildcard so partial names
    like 'Escher' match) and lists candidates with their scientific names.

    Dismiss: {'taxid': str, 'name': str} when the user picks an entry, or
    None on cancel. The parent (SpeciesPickerModal) then drives the actual
    Kazusa fetch.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    def __init__(self, initial_query: str) -> None:
        super().__init__()
        self._initial_query = initial_query
        self._searching     = False
        self._hits: list[dict] = []

    def compose(self) -> ComposeResult:
        with Vertical(id="ncbi-box"):
            yield Static(" NCBI Taxonomy  —  Pick a Species ", id="ncbi-title")
            yield Label("Refine search  (partial names OK — 'Escher' → Escherichia*)")
            yield Input(value=self._initial_query,
                        placeholder="genus or species (e.g. Escherichia, Homo sapiens)",
                        id="ncbi-query")
            yield Static("", id="ncbi-info", markup=True)
            yield ListView(id="ncbi-list")
            with Horizontal(id="ncbi-btns"):
                yield Button("Fetch Selected", id="btn-ncbi-use",
                             variant="primary", disabled=True)
                yield Button("Cancel  [Esc]", id="btn-ncbi-cancel")

    def on_mount(self) -> None:
        if self._initial_query:
            self._start_search(self._initial_query)
        self.query_one("#ncbi-query", Input).focus()

    def _start_search(self, query: str) -> None:
        if self._searching or not query.strip():
            return
        self._searching = True
        self.query_one("#ncbi-info", Static).update(
            f"[yellow]Searching NCBI for '{query}'…[/yellow]"
        )
        self.query_one("#btn-ncbi-use", Button).disabled = True
        self._do_search(query.strip())

    @work(thread=True)
    def _do_search(self, query: str) -> None:
        try:
            hits, total, msg = _ncbi_taxid_search(query)
        except Exception as exc:
            _log.exception("NCBI taxonomy search worker failed for %r", query)
            hits, total, msg = [], 0, f"Search failed: {exc}"
        self.app.call_from_thread(self._search_done, hits, total, msg)

    def _search_done(self, hits: list, total: int, msg: str) -> None:
        self._searching = False
        if not self.is_mounted:
            return
        self._hits = hits
        lv = self.query_one("#ncbi-list", ListView)
        lv.clear()
        for h in hits:
            lv.append(ListItem(Label(
                f"{h['name']}  [dim](taxid {h['taxid']})[/dim]",
                markup=True,
            )))
        info = self.query_one("#ncbi-info", Static)
        if hits:
            info.update(f"[dim]{msg}[/dim]")
        else:
            info.update(f"[red]{msg}[/red]")

    @on(Input.Submitted, "#ncbi-query")
    def _on_submit(self, _event: Input.Submitted) -> None:
        self._start_search(self.query_one("#ncbi-query", Input).value)

    @on(ListView.Highlighted, "#ncbi-list")
    def _list_highlighted(self, _) -> None:
        lv = self.query_one("#ncbi-list", ListView)
        if lv.index is None:
            return
        self.query_one("#btn-ncbi-use", Button).disabled = False

    @on(ListView.Selected, "#ncbi-list")
    def _list_selected(self, _) -> None:
        lv = self.query_one("#ncbi-list", ListView)
        if lv.index is None:
            return
        self.dismiss(self._hits[lv.index])

    @on(Button.Pressed, "#btn-ncbi-use")
    def _use(self, _) -> None:
        lv = self.query_one("#ncbi-list", ListView)
        if lv.index is None:
            return
        self.dismiss(self._hits[lv.index])

    @on(Button.Pressed, "#btn-ncbi-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Species picker (shared modal for codon-table selection) ────────────────────

class SpeciesPickerModal(ModalScreen):
    """Reusable codon-table picker — any modal that cares about codon usage
    (Mutagenize, future codon-optimize, future gene-synthesis) can
    ``push_screen(SpeciesPickerModal(), callback=...)``. The callback
    receives the selected entry dict (with 'raw' as tuples) or None.

    Shows the persistent registry with a substring filter. Users can fetch
    any NCBI taxid from Kazusa via a worker thread; the fetched table is
    added to the registry automatically so it's available in future
    sessions.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    def __init__(self) -> None:
        super().__init__()
        self._fetching = False

    def compose(self) -> ComposeResult:
        with Vertical(id="sp-box"):
            yield Static(" Codon Usage Table  —  Pick or Fetch ", id="sp-title")
            yield Label("Filter species")
            yield Input(placeholder="search by genus, species, or taxid "
                                    "(e.g. Escherichia, coli, 9606)",
                        id="sp-filter")
            yield ListView(id="sp-list")
            yield Static("", id="sp-info", markup=True)
            with Horizontal(id="sp-fetch-row"):
                yield Input(placeholder="taxid or name (e.g. 9606, Homo sapiens, "
                                        "Escherichia coli)",
                            id="sp-taxid")
                yield Input(placeholder="Display name (optional)",
                            id="sp-name")
                yield Button("Fetch from Kazusa", id="btn-sp-fetch",
                             variant="primary")
            with Horizontal(id="sp-btns"):
                yield Button("Use Selected", id="btn-sp-use", variant="primary",
                             disabled=True)
                yield Button("Delete", id="btn-sp-delete", disabled=True)
                yield Button("Cancel  [Esc]", id="btn-sp-cancel")

    def on_mount(self) -> None:
        self._refresh_list("")
        self.query_one("#sp-filter", Input).focus()

    def _refresh_list(self, query: str) -> None:
        lv = self.query_one("#sp-list", ListView)
        lv.clear()
        self._entries: list = _codon_search(query)
        for e in self._entries:
            tax = f" (taxid {e['taxid']})" if e.get("taxid") else ""
            src = e.get("source", "user")
            tag = f"[{src}]"
            lv.append(ListItem(Label(f"{e['name']}{tax}   [dim]{tag}[/dim]",
                                     markup=True)))
        info = self.query_one("#sp-info", Static)
        if not self._entries:
            info.update("[dim]No matching entries. Use the fetch row below to "
                        "import a new table from Kazusa.[/dim]")
        else:
            summary = self._genus_summary(query, self._entries)
            if summary:
                info.update(f"[dim]{summary}[/dim]")
            else:
                info.update(f"[dim]{len(self._entries)} table(s) in library.[/dim]")
        self.query_one("#btn-sp-use", Button).disabled = True
        self.query_one("#btn-sp-delete", Button).disabled = True

    @staticmethod
    def _genus_summary(query: str, entries: list) -> str:
        """When the filter matches a genus and all shown entries share it,
        return a 'N species of <Genus>' string. Otherwise ''."""
        q = (query or "").strip()
        if not q or not entries:
            return ""
        genera = {str(e.get("name", "")).split()[0]
                  for e in entries
                  if str(e.get("name", "")).split()}
        if len(genera) != 1:
            return ""
        genus = next(iter(genera))
        if not genus.lower().startswith(q.lower()):
            return ""
        n = len(entries)
        noun = "entry" if n == 1 else "entries"
        return f"{n} {noun} in genus {genus}."

    @on(Input.Changed, "#sp-filter")
    def _filter_changed(self, event: Input.Changed) -> None:
        self._refresh_list(event.value)

    @on(ListView.Selected, "#sp-list")
    def _list_selected(self, _) -> None:
        lv = self.query_one("#sp-list", ListView)
        if lv.index is None:
            return
        self.query_one("#btn-sp-use", Button).disabled = False
        entry = self._entries[lv.index]
        self.query_one("#btn-sp-delete", Button).disabled = (
            entry.get("source") == "builtin"
        )

    @on(ListView.Highlighted, "#sp-list")
    def _list_highlighted(self, _) -> None:
        lv = self.query_one("#sp-list", ListView)
        if lv.index is None:
            return
        entry = self._entries[lv.index]
        self.query_one("#btn-sp-use", Button).disabled = False
        self.query_one("#btn-sp-delete", Button).disabled = (
            entry.get("source") == "builtin"
        )

    @on(Button.Pressed, "#btn-sp-use")
    def _use(self, _) -> None:
        lv = self.query_one("#sp-list", ListView)
        if lv.index is None:
            return
        self.dismiss(self._entries[lv.index])

    @on(Button.Pressed, "#btn-sp-delete")
    def _delete(self, _) -> None:
        lv = self.query_one("#sp-list", ListView)
        if lv.index is None:
            return
        entry = self._entries[lv.index]
        if entry.get("source") == "builtin":
            return
        all_entries = _codon_tables_load()
        kept = [e for e in all_entries
                if (e.get("taxid") or e.get("name")) !=
                   (entry.get("taxid") or entry.get("name"))]
        _codon_tables_save(kept)
        self._refresh_list(self.query_one("#sp-filter", Input).value)

    @on(Button.Pressed, "#btn-sp-fetch")
    def _fetch(self, _) -> None:
        if self._fetching:
            return
        query = self.query_one("#sp-taxid", Input).value.strip()
        name  = self.query_one("#sp-name", Input).value.strip()
        info  = self.query_one("#sp-info", Static)
        if not query:
            info.update("[red]Enter an NCBI taxid or species/genus name.[/red]")
            return
        if query.isdigit():
            # Numeric taxid: go straight to Kazusa
            self._fetching = True
            self.query_one("#btn-sp-fetch", Button).disabled = True
            info.update(f"[yellow]Fetching taxid {query} from Kazusa…[/yellow]")
            self._do_fetch(query, name)
            return
        # Non-numeric: push the NCBI picker sub-modal. Button stays enabled
        # during the sub-modal so Esc-cancel can return to a clean state.
        def _picked(hit: "dict | None") -> None:
            if hit is None:
                try:
                    self._refresh_list(self.query_one("#sp-filter", Input).value)
                except NoMatches:
                    pass
                return
            taxid = hit["taxid"]
            display = name or hit.get("name") or f"Species (taxid {taxid})"
            self._fetching = True
            try:
                self.query_one("#btn-sp-fetch", Button).disabled = True
                self.query_one("#sp-info", Static).update(
                    f"[yellow]Fetching taxid {taxid} ({display}) from Kazusa…[/yellow]"
                )
            except NoMatches:
                pass
            self._do_fetch(taxid, display)

        self.app.push_screen(NcbiTaxonPickerModal(query), callback=_picked)

    @work(thread=True)
    def _do_fetch(self, taxid: str, name: str) -> None:
        try:
            raw, msg = _codon_fetch_kazusa(taxid)
        except Exception as exc:
            _log.exception("Kazusa fetch worker failed for taxid %s", taxid)
            raw, msg = None, f"Fetch failed: {exc}"
        self.app.call_from_thread(self._fetch_done, taxid, name, raw, msg)

    def _fetch_done(self, taxid: str, name: str,
                    raw: "dict | None", msg: str) -> None:
        self._fetching = False
        # If the user dismissed the modal mid-fetch, persist the result (so
        # they don't lose a successful HTTP round-trip) but skip the UI calls.
        if not self.is_mounted:
            if raw is not None:
                display = name or f"Species (taxid {taxid})"
                try:
                    _codon_tables_add(display, taxid, raw, source="kazusa")
                except Exception:
                    _log.exception("Codon-table add failed for taxid %s", taxid)
            return
        try:
            info = self.query_one("#sp-info", Static)
            btn  = self.query_one("#btn-sp-fetch", Button)
            btn.disabled = False
            if raw is None:
                info.update(f"[red]{msg}[/red]")
                return
            display = name or f"Species (taxid {taxid})"
            _codon_tables_add(display, taxid, raw, source="kazusa")
            info.update(f"[green]{msg} — added as '{display}'.[/green]")
            self._refresh_list(self.query_one("#sp-filter", Input).value)
            for i, e in enumerate(self._entries):
                if str(e.get("taxid")) == str(taxid):
                    self.query_one("#sp-list", ListView).index = i
                    break
        except Exception:
            _log.exception("SpeciesPickerModal fetch-callback failed")

    @on(Button.Pressed, "#btn-sp-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Mutagenize helpers (preview widget + AA picker sub-modal) ─────────────────

class _MutPreview(Static):
    """Focus-and-click-aware Static for the Mutagenize CDS preview.

    Single click places the cursor on the clicked AA (and takes focus so
    subsequent keys are routed here). Double click OR Enter posts
    `AARequested(aa_index, aa_letter)`, which the parent handles by
    opening the AA picker. Arrow keys move the cursor — Left/Right step
    by one AA, Up/Down step by one displayed row's worth of AAs.

    Owns its own render state (DNA, mutation, cursor) so it can redraw
    itself after any cursor change without the parent having to re-run
    its full `_update_preview` pipeline.
    """

    can_focus = True

    BINDINGS = [
        Binding("left",  "cursor_left",    "Prev AA",       show=False),
        Binding("right", "cursor_right",   "Next AA",       show=False),
        Binding("up",    "cursor_up",      "Prev row",      show=False),
        Binding("down",  "cursor_down",    "Next row",      show=False),
        Binding("enter", "cursor_request", "Mutate",        show=False),
    ]

    class AARequested(Message):
        """Emitted when the user commits to mutating the focused AA
        (via double-click or Enter)."""
        def __init__(self, aa_index: int, aa_letter: str) -> None:
            self.aa_index  = aa_index
            self.aa_letter = aa_letter
            super().__init__()

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        # Source content (what parent passes in)
        self._cds_dna_src:      str  = ""
        self._mutation_src:     "dict | None" = None
        self._protein_override: str  = ""
        self._line_width:       int  = 90
        # Derived, post-mutation display state
        self._eff_dna:     str  = ""
        self._protein:     str  = ""
        self._pad:         int  = 0
        self._dna_mode:    bool = False
        # Cursor (-1 = not placed yet)
        self._cursor_aa:   int  = -1

    def bind_content(self, *, dna: str = "", mutation: "dict | None" = None,
                     protein_override: str = "", line_width: int = 90) -> None:
        """Replace the content being previewed. Resets cursor to -1 so it
        doesn't point into a stale AA when the CDS changes."""
        self._cds_dna_src      = dna or ""
        self._mutation_src     = mutation
        self._protein_override = protein_override or ""
        self._line_width       = max(1, line_width)
        self._cursor_aa        = -1
        self._recompute_display()
        self._render_and_update()

    def _recompute_display(self) -> None:
        if self._cds_dna_src:
            lw  = max(3, (self._line_width // 3) * 3)
            dna = self._cds_dna_src.upper()
            if self._mutation_src:
                mut_c  = (self._mutation_src.get("mut_codon") or "").upper()
                try:
                    nt_pos = int(self._mutation_src.get("nt_position") or 0)
                except (TypeError, ValueError):
                    nt_pos = 0
                if mut_c and 1 <= nt_pos <= len(dna) - 2:
                    lo  = nt_pos - 1
                    dna = dna[:lo] + mut_c + dna[lo + 3:]
            self._eff_dna  = dna
            self._protein  = "".join(
                _MUT_CODON_TO_AA.get(dna[i:i + 3].upper(), "?")
                for i in range(0, len(dna) - 2, 3)
            )
            self._pad      = len(str(len(dna))) + 2
            self._dna_mode = True
        else:
            self._eff_dna  = ""
            self._protein  = (self._protein_override or "").upper()
            self._pad      = 0
            self._dna_mode = False
        if self._cursor_aa >= len(self._protein):
            self._cursor_aa = -1

    def _render_and_update(self) -> None:
        t = _mut_build_preview_text(
            self._cds_dna_src,
            protein_override=self._protein_override,
            mutation=self._mutation_src,
            line_width=self._line_width,
            cursor_aa=self._cursor_aa,
        )
        self.update(t)

    # ── Mouse ──────────────────────────────────────────────────────────

    def on_click(self, event: Click) -> None:
        if not self._protein:
            return
        self.focus()
        try:
            reg  = self.region
            vp_x = event.screen_x - reg.x
            vp_y = event.screen_y - reg.y
            if vp_x < 0 or vp_y < 0 or vp_x >= reg.width or vp_y >= reg.height:
                return
            content_row = vp_y + int(self.scroll_y)
        except Exception:
            return
        lw = (max(3, (self._line_width // 3) * 3) if self._dna_mode
              else max(1, self._line_width))
        aa_idx = _mut_click_to_aa_index(
            self._dna_mode, len(self._eff_dna), len(self._protein),
            lw, self._pad, vp_x, content_row,
        )
        if aa_idx < 0:
            return
        aa_letter = self._protein[aa_idx]
        if aa_letter == "?":
            return
        # Always place cursor on click
        self._cursor_aa = aa_idx
        self._render_and_update()
        # Double-click commits to mutation
        if event.chain >= 2:
            self.post_message(self.AARequested(aa_idx, aa_letter))

    # ── Keyboard navigation ────────────────────────────────────────────

    def _move_cursor(self, direction: str) -> None:
        if not self._protein:
            return
        new_idx = _mut_next_cursor(
            self._cursor_aa, len(self._protein),
            self._line_width, self._dna_mode, direction,
        )
        if new_idx != self._cursor_aa:
            self._cursor_aa = new_idx
            self._render_and_update()

    def action_cursor_left(self)  -> None: self._move_cursor("left")
    def action_cursor_right(self) -> None: self._move_cursor("right")
    def action_cursor_up(self)    -> None: self._move_cursor("up")
    def action_cursor_down(self)  -> None: self._move_cursor("down")

    def action_cursor_request(self) -> None:
        if not self._protein or self._cursor_aa < 0:
            return
        aa_letter = self._protein[self._cursor_aa]
        if aa_letter == "?":
            return
        self.post_message(self.AARequested(self._cursor_aa, aa_letter))


class AminoAcidPickerModal(ModalScreen):
    """Tiny picker shown when the user clicks an AA in the Mutagenize
    preview. Returns the selected one-letter AA on dismiss, or None on
    cancel. The WT amino at the clicked position is filtered out so
    the user can't pick a no-op mutation."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    # 20 standard amino acids + stop. Ordered alphabetically by one-letter.
    _AA_CATALOG: list[tuple[str, str, str]] = [
        ("A", "Ala", "Alanine"),       ("C", "Cys", "Cysteine"),
        ("D", "Asp", "Aspartate"),     ("E", "Glu", "Glutamate"),
        ("F", "Phe", "Phenylalanine"), ("G", "Gly", "Glycine"),
        ("H", "His", "Histidine"),     ("I", "Ile", "Isoleucine"),
        ("K", "Lys", "Lysine"),        ("L", "Leu", "Leucine"),
        ("M", "Met", "Methionine"),    ("N", "Asn", "Asparagine"),
        ("P", "Pro", "Proline"),       ("Q", "Gln", "Glutamine"),
        ("R", "Arg", "Arginine"),      ("S", "Ser", "Serine"),
        ("T", "Thr", "Threonine"),     ("V", "Val", "Valine"),
        ("W", "Trp", "Tryptophan"),    ("Y", "Tyr", "Tyrosine"),
        ("*", "***", "Stop codon"),
    ]

    def __init__(self, position: int, wt_aa: str) -> None:
        super().__init__()
        self._position = position
        self._wt_aa    = (wt_aa or "").upper()
        self._choices: list[str] = [
            a for (a, _, _) in self._AA_CATALOG if a != self._wt_aa
        ]

    def compose(self) -> ComposeResult:
        with Vertical(id="aa-pick-box"):
            yield Static(f" Mutate {self._wt_aa}{self._position}  →  ? ",
                         id="aa-pick-title")
            yield Label("[dim]Pick the replacement amino acid. "
                        "Esc to cancel.[/dim]", markup=True)
            items: list = []
            for (a, tl, fn) in self._AA_CATALOG:
                if a == self._wt_aa:
                    continue
                items.append(ListItem(Label(
                    f"[bold]{a}[/bold]   {tl}   [dim]{fn}[/dim]",
                    markup=True,
                )))
            yield ListView(*items, id="aa-pick-list")
            with Horizontal(id="aa-pick-btns"):
                yield Button("Cancel  [Esc]", id="btn-aa-pick-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#aa-pick-list", ListView).focus()
        except NoMatches:
            pass

    @on(ListView.Selected, "#aa-pick-list")
    def _selected(self, _event) -> None:
        lv = self.query_one("#aa-pick-list", ListView)
        if lv.index is None or lv.index >= len(self._choices):
            return
        self.dismiss(self._choices[lv.index])

    @on(Button.Pressed, "#btn-aa-pick-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Mutagenize modal ───────────────────────────────────────────────────────────

class MutagenizeModal(ModalScreen):
    """SOE-PCR site-directed mutagenesis primer designer.

    Pick a CDS feature from the loaded plasmid, enter a mutation string
    (e.g. W140F), and design the 4-primer SOE set — 2 constant outer primers
    (BsaI-AATG / BsaI-AACG tails → GB B3/B5 overhangs) and 1 inner pair for
    this mutation. Edge cases (mutation within 60 nt of either CDS end) are
    resolved by swapping the inner pair for a single modified outer primer
    and a 2-primer direct PCR.

    Save to primer library persists via `_save_primers` (atomic JSON with
    .bak, sacred invariant #7).
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, template_seq: str, feats: list, plasmid_name: str = ""):
        super().__init__()
        self._template     = (template_seq or "").upper()
        self._feats        = feats or []
        self._plasmid_name = plasmid_name
        self._outer:  "dict | None" = None
        self._inner:  "dict | None" = None
        self._cds_dna:   str        = ""   # CDS (5'→3', post-strand, post-harmonize)
        self._cds_meta:  "dict | None" = None
        # Codon-table for mut_codon picking + harmonization (library entry dict)
        self._codon_entry: "dict | None" = None
        # For library source — current plasmid's features + template
        self._lib_template: str = ""
        self._lib_feats:    list = []

    def compose(self) -> ComposeResult:
        with Vertical(id="mut-box"):
            yield Static(
                " Mutagenize  —  Golden Braid SOE-PCR Site-Directed Mutagenesis ",
                id="mut-title",
            )
            # ── Source selector ──
            yield Label("CDS source")
            yield Select(
                [
                    ("Current map features", "map"),
                    ("Plasmid library",       "lib"),
                    ("Protein sequence (harmonize)", "prot"),
                ],
                id="mut-source", value="map", allow_blank=False,
            )

            # ── Source-specific bodies (toggled via .display) ──
            with Vertical(id="mut-src-map"):
                yield Label("CDS feature  (from the loaded plasmid)")
                yield Select(self._build_cds_options(self._template, self._feats),
                             id="mut-cds", prompt="(select a CDS feature)")

            with Vertical(id="mut-src-lib"):
                yield Label("Plasmid  (from your library)")
                yield Select(self._build_library_options(),
                             id="mut-lib", prompt="(select a plasmid)")
                yield Label("CDS feature")
                yield Select([("(load a plasmid first)", "_none")],
                             id="mut-lib-cds", prompt="(select a CDS feature)")

            with Vertical(id="mut-src-prot"):
                yield Label("Protein sequence  (AA, 1-letter; stops optional)")
                yield TextArea("", id="mut-prot-aa")
                with Horizontal(id="mut-prot-row"):
                    yield Input(placeholder="Name (e.g. aeBlue)", id="mut-prot-name")
                    yield Button("Harmonize → CDS", id="btn-mut-harmonize",
                                 variant="primary")

            yield Static("", id="mut-cds-info", markup=True)

            # ── Codon table picker ──
            with Horizontal(id="mut-codon-row"):
                yield Static("Codon table: [bold]E. coli K12[/bold] (taxid 83333)",
                             id="mut-codon-label", markup=True)
                yield Button("Change…", id="btn-mut-codon", variant="default")

            # ── CDS preview  (DNA + AA, or AA only for the protein source) ──
            # Click-aware: clicking an AA opens AminoAcidPickerModal and
            # seeds #mut-input with the resulting W140F-style shorthand.
            yield _MutPreview("", id="mut-preview", markup=False)

            # ── Mutation + Design ──
            with Horizontal(id="mut-row2"):
                with Vertical(id="mut-mut-col"):
                    yield Label("Mutation  (e.g. W140F)")
                    yield Input(placeholder="W140F", id="mut-input")
                with Vertical(id="mut-btn-col"):
                    yield Label(" ")
                    yield Button("Design SOE Primers", id="btn-mut-design",
                                 variant="primary")

            yield Static("", id="mut-results", markup=True)
            with Horizontal(id="mut-btns"):
                yield Button("Save to Primer Library", id="btn-mut-save",
                             variant="primary", disabled=True)
                yield Button("Cancel  [Esc]", id="btn-mut-cancel")

    # ── Option builders ────────────────────────────────────────────────────

    def _build_cds_options(self, template: str, feats: list) -> list:
        opts: list = []
        total = len(template)
        for f in feats:
            if f.get("type") not in ("CDS", "gene"):
                continue
            start = f.get("start")
            end   = f.get("end")
            if start is None or end is None:
                continue
            strand = f.get("strand", 1)
            span_bp = _feat_len(start, end, total) if total else (end - start)
            if span_bp < 30 or span_bp % 3 != 0:
                continue
            label = f.get("label", f.get("type", "CDS"))
            strand_s = "+" if strand == 1 else "−"
            lo_disp = start + 1
            hi_disp = end if end > start else f"{total},1..{end}"
            opt_label = f"{label}  ({strand_s}{lo_disp}‥{hi_disp}, {span_bp} bp)"
            opts.append((opt_label, f"{start}:{end}:{strand}"))
        if not opts:
            opts = [("(no CDS features on this plasmid)", "_none")]
        return opts

    def _build_library_options(self) -> list:
        try:
            entries = _load_library()
        except Exception:
            _log.exception("Mutagenize: failed to load plasmid library")
            entries = []
        opts: list = []
        for i, e in enumerate(entries):
            nm = e.get("name", f"entry_{i}")
            if not e.get("gb_text"):
                continue
            opts.append((nm, str(i)))
        if not opts:
            opts = [("(plasmid library is empty)", "_none")]
        return opts

    def on_mount(self) -> None:
        # Seed built-in K12 via registry load
        try:
            _codon_tables_load()
            self._codon_entry = _codon_tables_get("83333")
        except Exception:
            _log.exception("Mutagenize: codon registry load failed")
            self._codon_entry = None
        self._apply_source("map")
        info = self.query_one("#mut-cds-info", Static)
        info.update("[dim]Pick a source and CDS, then enter a mutation like "
                    "[bold]W140F[/bold] and press Design.[/dim]")
        self._update_preview()

    # ── Source switching ──────────────────────────────────────────────────

    @on(Select.Changed, "#mut-source")
    def _source_changed(self, event: Select.Changed) -> None:
        if isinstance(event.value, str):
            self._apply_source(event.value)

    def _apply_source(self, src: str) -> None:
        self.query_one("#mut-src-map", Vertical).display  = (src == "map")
        self.query_one("#mut-src-lib", Vertical).display  = (src == "lib")
        self.query_one("#mut-src-prot", Vertical).display = (src == "prot")
        # Clear any previously-loaded CDS so the user knows they need to re-pick.
        # Also drop the last designed primers so the preview doesn't keep
        # highlighting a stale mutation from the previous source.
        self._cds_dna  = ""
        self._cds_meta = None
        self._outer    = None
        self._inner    = None
        self.query_one("#mut-cds-info", Static).update("")
        self.query_one("#btn-mut-save", Button).disabled = True
        self._update_preview()

    # ── Map source ────────────────────────────────────────────────────────

    @on(Select.Changed, "#mut-cds")
    def _map_cds_changed(self, event: Select.Changed) -> None:
        self._load_cds_from_feature(event.value, self._template, self._feats,
                                    origin="map")

    # ── Library source ────────────────────────────────────────────────────

    @on(Select.Changed, "#mut-lib")
    def _lib_changed(self, event: Select.Changed) -> None:
        val = event.value
        cds_select = self.query_one("#mut-lib-cds", Select)
        info = self.query_one("#mut-cds-info", Static)
        if not isinstance(val, str) or val == "_none" or not val.isdigit():
            cds_select.set_options([("(load a plasmid first)", "_none")])
            self._lib_template = ""
            self._lib_feats    = []
            return
        entries = _load_library()
        try:
            entry = entries[int(val)]
        except (IndexError, ValueError):
            info.update("[red]Library entry not found.[/red]")
            return
        gb = entry.get("gb_text", "")
        if not gb:
            info.update("[red]This library entry has no GenBank text.[/red]")
            return
        try:
            rec = _gb_text_to_record(gb)
        except Exception as exc:
            _log.exception("Mutagenize: library entry parse failed")
            info.update(f"[red]Could not parse library GenBank: {exc}[/red]")
            return
        self._lib_template = str(rec.seq).upper()
        total = len(self._lib_template)
        self._lib_feats    = []
        for f in rec.features:
            if f.type not in ("CDS", "gene"):
                continue
            loc = f.location
            # CompoundLocation (origin-wrap): int(loc.start)/int(loc.end) would
            # flatten parts into min/max, losing the wrap. Use parts explicitly
            # so wrap CDS features (end < start) round-trip correctly through
            # _mut_extract_cds. Single-part FeatureLocation.parts == [self].
            parts = list(loc.parts) if hasattr(loc, "parts") and loc.parts else [loc]
            start = int(parts[0].start) % total if total else int(parts[0].start)
            end   = int(parts[-1].end)   % total if total else int(parts[-1].end)
            strand = 1 if (loc.strand in (None, 1)) else -1
            label  = (f.qualifiers.get("label") or
                      f.qualifiers.get("gene")  or
                      f.qualifiers.get("product") or [f.type])[0]
            self._lib_feats.append({
                "type":   f.type, "label": label, "strand": strand,
                "start":  start,  "end":   end,
            })
        self._plasmid_name = entry.get("name", "")
        opts = self._build_cds_options(self._lib_template, self._lib_feats)
        cds_select.set_options(opts)

    @on(Select.Changed, "#mut-lib-cds")
    def _lib_cds_changed(self, event: Select.Changed) -> None:
        self._load_cds_from_feature(event.value, self._lib_template,
                                    self._lib_feats, origin="lib")

    def _load_cds_from_feature(self, val, template: str, feats: list,
                               origin: str) -> None:
        info = self.query_one("#mut-cds-info", Static)
        if not isinstance(val, str) or val == "_none" or ":" not in val:
            self._cds_dna  = ""
            self._cds_meta = None
            info.update("")
            return
        try:
            s, e, st = val.split(":")
            start, end, strand = int(s), int(e), int(st)
        except ValueError:
            info.update("[red]Malformed CDS selection.[/red]")
            return
        cds = _mut_extract_cds(template, start, end, strand)
        if len(cds) < 30 or len(cds) % 3 != 0:
            info.update("[red]CDS too short or not a multiple of 3.[/red]")
            self._cds_dna  = ""
            self._cds_meta = None
            return
        protein = _mut_translate(cds)
        label = "CDS"
        for f in feats:
            if (f.get("start") == start and f.get("end") == end
                    and f.get("strand", 1) == strand):
                label = f.get("label", f.get("type", "CDS"))
                break
        self._cds_dna  = cds
        self._cds_meta = {"start": start, "end": end, "strand": strand,
                          "name": label, "origin": origin}
        # A new CDS invalidates any mutation designed against the previous one
        self._outer = None
        self._inner = None
        self._update_cds_info(cds, protein, strand_label=(
            "+" if strand == 1 else "−"
        ))
        self._update_preview()

    def _update_cds_info(self, cds: str, protein: str,
                         strand_label: str = "·") -> None:
        atg = ("[green]ATG[/green]" if cds.startswith("ATG")
               else "[yellow]no ATG at 5'[/yellow]")
        stop_tag = cds[-3:] if len(cds) >= 3 else ""
        stop_s = (f"[green]{stop_tag}[/green]" if stop_tag in _MUT_STOPS
                  else f"[yellow]{stop_tag} (no stop)[/yellow]")
        info = self.query_one("#mut-cds-info", Static)
        info.update(
            f"  [dim]{len(cds)} nt · {len(protein)} aa · strand "
            f"{strand_label}[/dim]   start {atg}   stop {stop_s}"
        )

    # ── Protein-input source ──────────────────────────────────────────────

    @on(Button.Pressed, "#btn-mut-harmonize")
    def _harmonize(self, _) -> None:
        info = self.query_one("#mut-cds-info", Static)
        if self._codon_entry is None:
            info.update("[red]No codon table selected.[/red]")
            return
        aa_raw = self.query_one("#mut-prot-aa", TextArea).text
        # Strip whitespace, digits, FASTA header markers, and separators —
        # but NOT '*' (stop): that's a meaningful, invalid-in-middle char
        # we want to flag explicitly rather than silently drop.
        aa = re.sub(r"[\s\d>_\-]", "", aa_raw).upper()
        if not aa:
            info.update("[red]Enter a protein sequence.[/red]")
            return
        # Allow a single trailing stop but reject mid-sequence stops.
        if "*" in aa[:-1]:
            info.update("[red]Stop codon '*' not allowed in the middle of the "
                        "protein sequence.[/red]")
            return
        if aa.endswith("*"):
            aa = aa[:-1]
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        bad = sorted({c for c in aa if c not in valid_aa})
        if bad:
            info.update(f"[red]Invalid amino-acid letters: {''.join(bad)}[/red]")
            return
        try:
            cds = _codon_harmonize(aa, self._codon_entry["raw"])
            cds, fixes = _codon_fix_sites(
                cds, aa, self._codon_entry["raw"],
                sites={"BsaI": "GGTCTC"},  # only guard the tail enzyme
            )
        except Exception as exc:
            _log.exception("Mutagenize: harmonize failed")
            info.update(f"[red]Harmonization failed: {exc}[/red]")
            return
        name = self.query_one("#mut-prot-name", Input).value.strip() or "protein"
        self._cds_dna  = cds
        self._cds_meta = {"start": 0, "end": len(cds), "strand": 1,
                          "name": name, "origin": "prot"}
        self._plasmid_name = name
        # Fresh harmonized CDS invalidates any previously designed primers
        self._outer = None
        self._inner = None
        protein = _mut_translate(cds)
        atg = "[green]ATG[/green]" if cds.startswith("ATG") else "[yellow]no ATG[/yellow]"
        stop_tag = cds[-3:] if len(cds) >= 3 else ""
        stop_s = (f"[green]{stop_tag}[/green]" if stop_tag in _MUT_STOPS
                  else f"[yellow]{stop_tag}[/yellow]")
        fix_note = f" · {len(fixes)} BsaI fix(es)" if fixes else ""
        info.update(
            f"  [dim]{len(cds)} nt · {len(protein)} aa · harmonized · "
            f"CAI {_codon_cai(cds, self._codon_entry['raw']):.2f} · "
            f"GC {_codon_gc(cds):.1f}%{fix_note}[/dim]   "
            f"start {atg}   stop {stop_s}"
        )
        self._update_preview()

    # ── CDS preview  (DNA + AA; AA-only while typing protein) ────────────

    def _update_preview(self) -> None:
        """Refresh the `#mut-preview` pane by handing current state to the
        preview widget, which owns its own rendering + cursor management.
        Called on every state change that could alter what the user sees:
        source switch, CDS load, harmonize, design, or live typing in
        the protein textarea.
        """
        try:
            preview = self.query_one("#mut-preview", _MutPreview)
        except Exception:
            return
        lw = 90
        if self._cds_dna:
            preview.bind_content(
                dna=self._cds_dna,
                mutation=self._inner or None,
                line_width=lw,
            )
        else:
            aa_raw = ""
            try:
                if self.query_one("#mut-source", Select).value == "prot":
                    aa_raw = self.query_one("#mut-prot-aa", TextArea).text
            except Exception:
                aa_raw = ""
            aa_clean = re.sub(r"[\s\d>_\-*]", "", aa_raw or "").upper()
            preview.bind_content(protein_override=aa_clean, line_width=lw)

    @on(TextArea.Changed, "#mut-prot-aa")
    def _on_prot_aa_changed(self, _event) -> None:
        # Only live-update while we have no harmonized CDS yet; once the
        # user harmonizes, self._cds_dna takes over and further edits to
        # the textarea shouldn't clobber the DNA preview.
        if not self._cds_dna:
            self._update_preview()

    @on(_MutPreview.AARequested)
    def _on_preview_aa_requested(self, event: "_MutPreview.AARequested") -> None:
        """Preview fired AARequested (double-click or Enter on a focused
        AA) → open the picker → seed the mutation textbox with
        `{WT}{pos}{NEW}` (e.g. 'W140F')."""
        position = event.aa_index + 1   # biology convention: 1-based
        wt_aa    = event.aa_letter

        def _picked(new_aa: "str | None") -> None:
            if not new_aa:
                return
            try:
                self.query_one("#mut-input", Input).value = (
                    f"{wt_aa}{position}{new_aa}"
                )
            except Exception:
                _log.exception("Mutagenize: failed to set mutation input")

        self.app.push_screen(
            AminoAcidPickerModal(position, wt_aa), callback=_picked,
        )

    # ── Codon-table picker ───────────────────────────────────────────────

    @on(Button.Pressed, "#btn-mut-codon")
    def _pick_codon_table(self, _) -> None:
        self.app.push_screen(SpeciesPickerModal(), callback=self._codon_picked)

    def _codon_picked(self, entry: "dict | None") -> None:
        if not entry:
            return
        self._codon_entry = entry
        lbl = self.query_one("#mut-codon-label", Static)
        tax = f" (taxid {entry['taxid']})" if entry.get("taxid") else ""
        lbl.update(f"Codon table: [bold]{entry['name']}[/bold]{tax}")

    # ── Design ───────────────────────────────────────────────────────────

    @on(Input.Submitted, "#mut-input")
    def _mut_enter(self, _) -> None:
        self.query_one("#btn-mut-design", Button).press()

    @on(Button.Pressed, "#btn-mut-design")
    def _design(self, _) -> None:
        status = self.query_one("#mut-results", Static)
        if not self._cds_dna:
            status.update("[red]Load a CDS first (pick a source above).[/red]")
            return
        mut_val = self.query_one("#mut-input", Input).value.strip()
        if not mut_val:
            status.update("[red]Enter a mutation (e.g. W140F).[/red]")
            return
        try:
            wt_aa, pos, mut_aa = _mut_parse(mut_val)
        except ValueError as exc:
            status.update(f"[red]{exc}[/red]")
            return
        if mut_aa == wt_aa:
            status.update(f"[red]Position {pos} is already '{wt_aa}' — "
                          "WT and mutant are identical.[/red]")
            return
        protein = _mut_translate(self._cds_dna)
        if pos < 1 or pos > len(protein):
            status.update(f"[red]Position {pos} out of range "
                          f"(protein is {len(protein)} aa).[/red]")
            return
        actual = protein[pos - 1]
        if actual != wt_aa:
            status.update(f"[red]Position {pos} is '{actual}', not '{wt_aa}'.[/red]")
            return
        codon_raw = (self._codon_entry or {}).get("raw")
        try:
            outer = _mut_design_outer(self._cds_dna)
            inner = _mut_design_inner(self._cds_dna, pos, mut_aa, wt_aa,
                                      codon_table=codon_raw)
        except Exception as exc:
            _log.exception("Mutagenesis primer design failed")
            status.update(f"[red]Primer design failed: {exc}[/red]")
            return

        self._outer = outer
        self._inner = inner
        status.update(self._render_results())
        self.query_one("#btn-mut-save", Button).disabled = False
        self._update_preview()

    def _render_results(self) -> Text:
        t = Text()
        outer = self._outer or {}
        inner = self._inner or {}
        mutation = inner.get("mutation", "?")
        wt_codon = inner.get("wt_codon", "???")
        mut_codon = inner.get("mut_codon", "???")
        nt_pos   = inner.get("nt_position", 0)
        nt_chg   = inner.get("nt_changes", 0)
        t.append(f"── {mutation}  ·  codon {wt_codon}→{mut_codon} "
                 f"({nt_chg} nt change{'s' if nt_chg != 1 else ''}) "
                 f"·  CDS nt {nt_pos}–{nt_pos + 2} ──\n", style="bold")

        fwd_o = outer.get("fwd", {})
        rev_o = outer.get("rev", {})
        t.append("\nOuter primers  ", style="bold")
        t.append("(constant for this CDS — order once)\n", style="dim")
        t.append("  FWD  ", style="green bold")
        full = fwd_o.get("full", "")
        tl = len(_MUT_BSAI_FWD_TAIL)
        t.append(full[:tl], style="dim green")
        t.append(full[tl:], style="green")
        t.append(f"   Tm {fwd_o.get('tm_anneal', 0):.1f}°C  "
                 f"{len(full)} bp\n", style="dim")
        t.append("  REV  ", style="red bold")
        full = rev_o.get("full", "")
        t.append(full[:tl], style="dim red")
        t.append(full[tl:], style="red")
        t.append(f"   Tm {rev_o.get('tm_anneal', 0):.1f}°C  "
                 f"{len(full)} bp\n", style="dim")

        ec = inner.get("edge_case")
        if ec:
            side = "CDS start" if ec["near_start"] else "CDS end"
            mod  = ec["modified_outer"]
            which = "A" if ec["near_start"] else "B"
            frag_bp = ec["frag_a"] if ec["near_start"] else ec["frag_b"]
            t.append(f"\n⚠  Edge case — mutation too close to {side} "
                     f"(SOE fragment {which} = {frag_bp} nt < {_MUT_MIN_SOE_FRAG}).\n",
                     style="yellow bold")
            t.append("   Use the modified outer primer in a 2-primer direct PCR "
                     "(no SOE needed):\n", style="dim")
            t.append(f"  {mod['label']}  ", style="magenta bold")
            t.append(mod["full"], style="magenta")
            t.append(f"   Tm {mod['tm_anneal']:.1f}°C  {len(mod['full'])} bp\n",
                     style="dim")
            t.append(f"   Partner: {mod['partner']}\n", style="dim")
        else:
            best = inner["candidates"][0]
            t.append("\nInner pair  ", style="bold")
            t.append(f"(one per mutation — Tm {best['tm']:.1f}°C, "
                     f"{best['length']} bp)\n", style="dim")
            t.append("  FWD  ", style="green bold")
            t.append(best["fwd"], style="green")
            t.append(f"   Tm {best['tm']:.1f}°C  GC {best['gc']:.1f}%\n",
                     style="dim")
            t.append("  REV  ", style="red bold")
            t.append(best["rev"], style="red")
            t.append("   (revcomp of FWD)\n", style="dim")

        t.append("\nProtocol: ", style="bold")
        if ec:
            t.append("PCR with the modified outer + unchanged partner → "
                     "BsaI Golden Gate.\n", style="dim")
        else:
            t.append("PCR1 FWD_outer+REV_inner, PCR2 FWD_inner+REV_outer, "
                     "PCR3 joins A+B → BsaI Golden Gate.\n", style="dim")
        return t

    # ── Save to primer library ────────────────────────────────────────────

    @on(Button.Pressed, "#btn-mut-save")
    def _save(self, _) -> None:
        if not (self._outer and self._inner and self._cds_meta):
            return
        import datetime
        today = datetime.date.today().isoformat()
        construct = self._cds_meta.get("name") or self._plasmid_name or "CDS"
        # Cap length so an unusually long feature label can't produce
        # primer names that blow past filesystem/name limits.
        safe_construct = (re.sub(r"[^\w\-]+", "_", construct)
                          .strip("_")[:32]) or "CDS"
        mutation = self._inner["mutation"]

        existing  = _load_primers()
        seen_seqs = {e.get("sequence", "").upper() for e in existing}
        entries   = list(existing)

        def _upsert(name: str, seq: str, tm: float, ptype: str, strand: int) -> bool:
            """Insert a primer record, replacing any existing row with the
            same name. Returns False if the sequence is already stored under
            a different name (skips that primer so we don't silently clobber
            user-edited variants)."""
            if seq.upper() in seen_seqs and not any(
                e.get("name") == name and e.get("sequence", "").upper() == seq.upper()
                for e in entries
            ):
                return False
            nonlocal_entries = [e for e in entries if e.get("name") != name]
            nonlocal_entries.insert(0, {
                "name":        name,
                "sequence":    seq,
                "tm":          round(tm, 1),
                "primer_type": ptype,
                "source":      f"mutagenize:{safe_construct}:{mutation}",
                "pos_start":   -1,
                "pos_end":     -1,
                "strand":      strand,
                "date":        today,
                "status":      "Designed",
            })
            entries.clear()
            entries.extend(nonlocal_entries)
            seen_seqs.add(seq.upper())
            return True

        saved: list = []
        skipped: list = []
        outer_fwd = self._outer["fwd"]
        outer_rev = self._outer["rev"]
        # Outer primers — named per-construct (shared across mutations)
        for (name, p, ptype, strand) in [
            (f"OUTER_FWD_{safe_construct}", outer_fwd, "mutagenesis_outer_fwd", 1),
            (f"OUTER_REV_{safe_construct}", outer_rev, "mutagenesis_outer_rev", -1),
        ]:
            ok = _upsert(name, p["full"], p["tm_anneal"], ptype, strand)
            (saved if ok else skipped).append(name)

        ec = self._inner.get("edge_case")
        if ec:
            mod = ec["modified_outer"]
            strand = 1 if ec["near_start"] else -1
            name = f"{mod['label']}_{safe_construct}_{mutation}"
            ok = _upsert(name, mod["full"], mod["tm_anneal"],
                         "mutagenesis_modified_outer", strand)
            (saved if ok else skipped).append(name)
        else:
            best = self._inner["candidates"][0]
            for (name, seq, ptype, strand) in [
                (f"INNER_FWD_{safe_construct}_{mutation}", best["fwd"],
                 "mutagenesis_inner_fwd", 1),
                (f"INNER_REV_{safe_construct}_{mutation}", best["rev"],
                 "mutagenesis_inner_rev", -1),
            ]:
                ok = _upsert(name, seq, best["tm"], ptype, strand)
                (saved if ok else skipped).append(name)

        _save_primers(entries)
        parts = [f"Saved {len(saved)} primer{'s' if len(saved) != 1 else ''} to library"]
        if skipped:
            parts.append(f"({len(skipped)} duplicate sequence(s) skipped)")
        self.app.notify(" ".join(parts))
        self.dismiss(True)

    @on(Button.Pressed, "#btn-mut-cancel")
    def _cancel_btn(self, _) -> None:
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

    BINDINGS = [
        Binding("escape",  "cancel",     "Close",         show=True),
        Binding("m",       "noop",       "Mark (★)",      show=True, key_display="m"),
        Binding("shift+m", "noop",       "Mark All",      show=True, key_display="M"),
        Binding("shift+s", "noop",       "Change Status", show=True, key_display="S"),
        Binding("tab",     "focus_next", "",               show=False),
    ]

    def action_noop(self) -> None:
        """Placeholder — m/M are handled in on_key, but declared in
        BINDINGS so they appear in the Footer."""
        pass

    def check_action(self, action: str, parameters: tuple) -> bool | None:
        """Always allow Screen-level actions. Without this, the App's
        check_action (which blocks non-default-screen actions) was
        suppressing the Footer display of our m/M/S bindings."""
        return True

    def __init__(self, template_seq: str, feats: list[dict],
                 plasmid_name: str = ""):
        super().__init__()
        self._template     = template_seq.upper()
        self._feats        = feats
        self._plasmid_name = plasmid_name
        # Default part name = first non-RE feature label, NOT the plasmid name.
        # Users expect to design primers for a specific feature, not the whole
        # plasmid. Falls back to plasmid name if no features exist.
        self._default_part_name = plasmid_name
        for f in feats:
            if f.get("type") not in ("resite", "recut", "source"):
                self._default_part_name = f.get("label", plasmid_name)
                break
        self._det_result:  "dict | None" = None
        self._clo_result:  "dict | None" = None
        self._lib_selected: set[int] = set()   # multi-selected library rows

    def compose(self) -> ComposeResult:
        # Feature dropdown options
        feat_opts: list[tuple[str, str]] = []
        for f in self._feats:
            if f.get("type") in ("resite", "recut"):
                continue
            label = f.get("label", f.get("type", "?"))
            feat_opts.append(
                (f"{label}  ({f['start']+1}‥{f['end']})", f"{f['start']}-{f['end']}")
            )
        re_opts = _CLONING_RE_OPTIONS
        gb_opts = [
            (f"{k}  ({v[0]}: {v[1]}→{v[2]})", k) for k, v in _GB_POSITIONS.items()
        ]
        source_opts = [
            ("Feature from map",    "feature"),
            ("Custom sequence",     "custom"),
        ]

        yield Header()
        with Vertical(id="pd-box"):
            yield Static(" Primer Design  —  Primer3 ", id="pd-title")

            # ═══ Open-book split: workflow on the left, library on the right
            with Horizontal(id="pd-book"):

                # ╔════════════════════════ LEFT PAGE ════════════════════╗
                with Vertical(id="pd-left-page"):

                    # ─── 1. TEMPLATE ────────────────────────────────────
                    with Vertical(classes="pd-section",
                                  id="pd-template-section"):
                        yield Static("1. TEMPLATE", classes="pd-section-hdr")

                        # Source toggle + plasmid chooser
                        with Horizontal(id="pd-src-row"):
                            yield Label("Source:", classes="pd-fld-lbl")
                            yield Select(source_opts, id="pd-source",
                                         value="feature", allow_blank=False)
                            yield Label("Plasmid:", classes="pd-fld-lbl",
                                        id="pd-plasmid-lbl")
                            yield Static(
                                self._plasmid_name or "(no plasmid loaded)",
                                id="pd-plasmid-name",
                            )
                            yield Button("Change",
                                         id="btn-pd-pickplasmid",
                                         variant="default")

                        # Feature-from-map row
                        with Horizontal(id="pd-src-feature",
                                        classes="pd-source-panel"):
                            with Vertical(id="pd-feat-col"):
                                yield Label("Feature")
                                yield Select(feat_opts, id="pd-feat",
                                             prompt="(select feature)")
                            with Vertical(id="pd-start-col"):
                                yield Label("Start")
                                yield Input(placeholder="1", id="pd-start",
                                            type="integer")
                            with Vertical(id="pd-end-col"):
                                yield Label("End")
                                yield Input(
                                    placeholder=str(len(self._template))
                                    if self._template else "",
                                    id="pd-end", type="integer")
                            with Vertical(id="pd-name-col"):
                                yield Label("Part name")
                                yield Input(
                                    value="", id="pd-part-name",
                                    placeholder=self._default_part_name)

                        # Custom sequence row
                        with Vertical(id="pd-src-custom",
                                      classes="pd-source-panel"):
                            yield Label("Paste sequence — highlight to "
                                        "target a region, or the full "
                                        "entry will be used:")
                            yield TextArea("", id="pd-custom-seq")

                        yield Static(
                            "[dim]Tip: enter Start > End to design primers "
                            "across the origin (e.g. 2900..200 on a 3 kb "
                            "plasmid).[/dim]",
                            id="pd-wrap-hint", markup=True,
                        )
                        yield Static("", id="pd-feat-info", markup=True)

                    # ─── 2. MODE ────────────────────────────────────────
                    with Vertical(classes="pd-section",
                                  id="pd-mode-section"):
                        yield Static("2. MODE", classes="pd-section-hdr")
                        with RadioSet(id="pd-mode-radio"):
                            yield RadioButton(
                                "Detection  [dim](diagnostic PCR)[/dim]",
                                id="rb-detection", value=True)
                            yield RadioButton(
                                "Cloning  [dim](RE tails + GCGC)[/dim]",
                                id="rb-cloning")
                            yield RadioButton(
                                "Golden Braid  [dim](L0 domestication)[/dim]",
                                id="rb-goldenbraid")
                            yield RadioButton(
                                "Generic  [dim](binding only)[/dim]",
                                id="rb-generic")

                    # ─── 3. PARAMETERS (+ Design button) ────────────────
                    # Cloning / GB are restacked to 2-3 internal rows so
                    # the whole panel fits inside the narrower left page.
                    with Vertical(classes="pd-section",
                                  id="pd-params-section"):
                        yield Static("3. PARAMETERS",
                                     classes="pd-section-hdr")

                        # Detection — still one row
                        with Horizontal(id="pd-panel-det",
                                        classes="pd-mode-panel"):
                            yield Label("Product")
                            yield Input(value="450", id="pd-det-min",
                                        type="integer")
                            yield Label("–")
                            yield Input(value="550", id="pd-det-max",
                                        type="integer")
                            yield Label("bp")
                            yield Label(" Tm")
                            yield Input(value="60", id="pd-det-tm",
                                        type="integer")
                            yield Label("°C")
                            yield Label(" Len")
                            yield Input(value="25", id="pd-det-len",
                                        type="integer")

                        # Cloning — stacked in 3 rows (fits half-width)
                        with Vertical(id="pd-panel-clo",
                                      classes="pd-mode-panel"):
                            with Horizontal(classes="pd-mode-row"):
                                yield Label("5' RE  ")
                                yield Select(re_opts, id="pd-re5",
                                             value="EcoRI")
                                yield Label("  or  ")
                                yield Input(placeholder="GAATTC",
                                            id="pd-cust5")
                            with Horizontal(classes="pd-mode-row"):
                                yield Label("3' RE  ")
                                yield Select(re_opts, id="pd-re3",
                                             value="BamHI")
                                yield Label("  or  ")
                                yield Input(placeholder="GGATCC",
                                            id="pd-cust3")
                            with Horizontal(classes="pd-mode-row"):
                                yield Label("Tm     ")
                                yield Input(value="60", id="pd-clo-tm",
                                            type="integer")
                                yield Label(" °C")

                        # Golden Braid — stacked in 2 rows
                        with Vertical(id="pd-panel-gb",
                                      classes="pd-mode-panel"):
                            with Horizontal(classes="pd-mode-row"):
                                yield Label("Part type  ")
                                yield Select(gb_opts, id="pd-gb-type",
                                             value="CDS")
                            yield Static("", id="pd-gb-oh-info",
                                         markup=True)

                        # Generic — still one row
                        with Horizontal(id="pd-panel-gen",
                                        classes="pd-mode-panel"):
                            yield Label("Tm")
                            yield Input(value="60", id="pd-gen-tm",
                                        type="integer")
                            yield Label("°C")
                            yield Label(" Source ID")
                            yield Input(placeholder="e.g. synthetic frag",
                                        id="pd-gen-source")

                        # Design button — docked at the bottom of section 3
                        with Horizontal(id="pd-design-row"):
                            yield Button("Design primers",
                                         id="btn-pd-design",
                                         variant="primary")

                # ╔═══════════════════════ RIGHT PAGE ════════════════════╗
                with Vertical(id="pd-right-page"):

                    # ─── RESULTS — sits above the library so newly-
                    # designed primers appear next to where you'll save them
                    with Vertical(classes="pd-section",
                                  id="pd-results-section"):
                        yield Static("RESULTS", classes="pd-section-hdr")
                        yield Static(
                            self._RESULTS_EMPTY_HINT,
                            id="pd-results", markup=True,
                        )
                        # Names on their own row so inputs fill the full
                        # right-page width (~2× longer than when they
                        # shared a row with the buttons).
                        with Horizontal(id="pd-result-names"):
                            yield Input(id="pd-fwd-name",
                                        placeholder="fwd primer name")
                            yield Input(id="pd-rev-name",
                                        placeholder="rev primer name")
                        with Horizontal(id="pd-result-actions"):
                            yield Button("Save to Library",
                                         id="btn-pd-save",
                                         variant="primary", disabled=True)
                            yield Button("Add to Map",
                                         id="btn-pdlib-addmap",
                                         variant="default")

                    # ─── PRIMER LIBRARY ─────────────────────────────────
                    with Horizontal(id="pd-lib-hdr-row"):
                        yield Static("PRIMER LIBRARY", id="pd-lib-hdr")
                        yield Button("Rename", id="btn-pdlib-rename",
                                     variant="default")
                        yield Button("Delete", id="btn-pdlib-del",
                                     variant="error")
                        yield Button("Close", id="btn-pd-close",
                                     variant="default")
                    yield DataTable(id="pd-lib-table", cursor_type="row",
                                    zebra_stripes=True)
        yield Footer()

    # RadioButton id → internal mode name. Each mode has a matching
    # parameter panel in the 3. PARAMETERS section; only one is shown.
    _RB_TO_MODE = {
        "rb-detection":   "detection",
        "rb-cloning":     "cloning",
        "rb-goldenbraid": "goldenbraid",
        "rb-generic":     "generic",
    }
    _MODE_PANELS = {
        "detection":   "#pd-panel-det",
        "cloning":     "#pd-panel-clo",
        "goldenbraid": "#pd-panel-gb",
        "generic":     "#pd-panel-gen",
    }

    _SOURCE_PANELS = {
        "feature":  "#pd-src-feature",
        "custom":   "#pd-src-custom",
    }

    def _current_mode(self) -> str:
        """Resolve the active RadioButton → mode name."""
        try:
            rs = self.query_one("#pd-mode-radio", RadioSet)
            # RadioSet's pressed_button is the currently-selected button
            pressed = rs.pressed_button
            if pressed is not None:
                return self._RB_TO_MODE.get(pressed.id or "", "detection")
        except Exception:
            pass
        return "detection"

    def _switch_mode(self, mode: str) -> None:
        """Show the parameter panel for `mode`, hide the others."""
        for m, sel in self._MODE_PANELS.items():
            try:
                self.query_one(sel).display = (m == mode)
            except NoMatches:
                pass
        if mode == "goldenbraid":
            try:
                self._update_gb_oh()
            except Exception:
                _log.exception("goldenbraid overhang refresh failed")

    def _switch_source(self, src: str) -> None:
        for s, sel in self._SOURCE_PANELS.items():
            try:
                self.query_one(sel).display = (s == src)
            except NoMatches:
                pass
        # Hide the plasmid chooser when source=custom (it's irrelevant)
        try:
            for wid in ("pd-plasmid-lbl", "pd-plasmid-name",
                        "btn-pd-pickplasmid"):
                self.query_one(f"#{wid}").display = (src == "feature")
        except NoMatches:
            pass

    @on(Select.Changed, "#pd-source")
    def _source_changed(self, event: Select.Changed) -> None:
        val = event.value
        if isinstance(val, str) and val in self._SOURCE_PANELS:
            self._switch_source(val)

    @on(RadioSet.Changed, "#pd-mode-radio")
    def _mode_changed(self, event: RadioSet.Changed) -> None:
        """User clicked a different mode in the wizard."""
        rb_id = (event.pressed.id or "") if event.pressed else ""
        mode = self._RB_TO_MODE.get(rb_id, "detection")
        self._switch_mode(mode)

    @on(Button.Pressed, "#btn-pd-pickplasmid")
    def _pick_plasmid(self, _) -> None:
        """Open the plasmid picker modal; on selection, swap the template
        and feature list to the chosen plasmid."""
        current_id = None
        rec = getattr(self.app, "_current_record", None)
        if rec is not None:
            current_id = rec.id

        def _on_result(entry_id):
            if entry_id is None:
                return
            for entry in _load_library():
                if entry.get("id") == entry_id:
                    gb = entry.get("gb_text", "")
                    if not gb:
                        self.app.notify("Library entry has no sequence.",
                                        severity="warning")
                        return
                    try:
                        new_rec = _gb_text_to_record(gb)
                    except Exception as exc:
                        _log.exception("Library load for primer-design failed")
                        self.app.notify(f"Failed to load: {exc}",
                                        severity="error")
                        return
                    self._template = str(new_rec.seq).upper()
                    self._plasmid_name = new_rec.name
                    self._feats = self._parse_features_from_record(new_rec)
                    self.query_one("#pd-plasmid-name", Static).update(
                        new_rec.name)
                    self._update_feature_dropdown()
                    self.app.notify(f"Loaded {new_rec.name} as primer template.")
                    return
            self.app.notify("Entry not found.", severity="warning")

        self.app.push_screen(PlasmidPickerModal(current_id),
                             callback=_on_result)

    def _parse_features_from_record(self, record) -> list[dict]:
        """Minimal feature parse from a SeqRecord — matches the keys used
        by the compose-time feat_opts list."""
        feats = []
        for feat in record.features:
            if feat.type == "source":
                continue
            s = int(feat.location.start)
            e = int(feat.location.end)
            strand = getattr(feat.location, "strand", 1) or 1
            feats.append({
                "type":   feat.type,
                "start":  s,
                "end":    e,
                "strand": strand,
                "label":  _feat_label(feat),
                "color":  "white",
            })
        return feats

    def _update_feature_dropdown(self) -> None:
        """Rebuild the feature dropdown options from self._feats."""
        sel = self.query_one("#pd-feat", Select)
        feat_opts: list[tuple[str, str]] = []
        for f in self._feats:
            if f.get("type") in ("resite", "recut"):
                continue
            label = f.get("label", f.get("type", "?"))
            feat_opts.append(
                (f"{label}  ({f['start']+1}‥{f['end']})",
                 f"{f['start']}-{f['end']}")
            )
        sel.set_options(feat_opts)

    def on_mount(self) -> None:
        t = self.query_one("#pd-lib-table", DataTable)
        t.add_columns("Name", "Sequence", "Len", "Tm", "Type", "Source", "Date", "Status")
        self._refresh_library_table()
        # Source defaults to Feature from map; mode defaults to Detection
        # (the Detection RadioButton has value=True).
        self._switch_source("feature")
        self._switch_mode("detection")

    _STATUS_COLORS = {
        "Designed":  "cyan",
        "Ordered":   "yellow",
        "Validated": "green",
    }

    _RESULTS_EMPTY_HINT = (
        "[dim]Set a template and parameters, then "
        "press  [bold]Design primers[/bold]  to "
        "generate a primer pair.[/dim]"
    )

    def _reset_for_new_design(self) -> None:
        """Clear the state specific to a single primer-pair design so
        the user can start fresh for the next pair. Called after a
        successful save-to-library.

        Reset:
          - primer-pair output (results cache, names, Save button,
            Results pane hint)
          - feature selection + start/end + part name + feature info
            (so the user picks a new region for the next design)

        Preserved:
          - plasmid (Source, Change-plasmid button)
          - mode + mode-specific parameters (Tm, product size, RE
            sites, etc.) — a common workflow is to tweak params and
            redesign with the same settings
          - custom-sequence TextArea content (if the user pasted one)
        """
        self._det_result = None
        self._clo_result = None
        try:
            # Per-pair output
            self.query_one("#pd-fwd-name", Input).value = ""
            self.query_one("#pd-rev-name", Input).value = ""
            self.query_one("#btn-pd-save", Button).disabled = True
            self.query_one("#pd-results", Static).update(
                self._RESULTS_EMPTY_HINT
            )

            # Feature + region selection. Select.clear() sets the value
            # to Select.BLANK (= False), resetting to the prompt state.
            self.query_one("#pd-feat", Select).clear()
            self.query_one("#pd-start", Input).value = ""
            self.query_one("#pd-end",   Input).value = ""
            self.query_one("#pd-part-name", Input).value = ""
            self.query_one("#pd-feat-info", Static).update("")
        except Exception:
            # Widgets may not exist if the screen is being dismissed
            # concurrently; ignore silently.
            pass

    def _refresh_library_table(self) -> None:
        t = self.query_one("#pd-lib-table", DataTable)
        saved_cursor = t.cursor_row if t.row_count > 0 else 0
        t.clear()
        primers = _load_primers()
        self._lib_selected &= set(range(len(primers)))
        for i, p in enumerate(primers):
            seq    = p.get("sequence", "")
            marked = i in self._lib_selected
            mark   = "★ " if marked else "  "
            status = p.get("status", "Designed")
            s_color = self._STATUS_COLORS.get(status, "white")
            t.add_row(
                Text(mark + p.get("name", "?"), style="bold"),
                Text(seq[:30], style="dim color(252)"),
                f"{len(seq)} nt",
                f"{p.get('tm', 0):.1f}°C",
                p.get("primer_type", "?"),
                p.get("source", ""),
                p.get("date", ""),
                Text(status, style=s_color),
            )
        if primers and 0 <= saved_cursor < len(primers):
            t.move_cursor(row=saved_cursor)

    # ── GB type selector ─────────────────────────────────────────────────

    @on(Select.Changed, "#pd-gb-type")
    def _gb_type_changed(self, _event) -> None:
        self._update_gb_oh()

    def _update_gb_oh(self) -> None:
        sel = self.query_one("#pd-gb-type", Select)
        val = sel.value
        if not isinstance(val, str) or val not in _GB_POSITIONS:
            return
        pos, oh5, oh3 = _GB_POSITIONS[val]
        self.query_one("#pd-gb-oh-info", Static).update(
            f"  [dim]{pos}[/dim]   "
            f"5′: [bold cyan]{oh5}[/bold cyan]  →  "
            f"3′: [bold cyan]{oh3}[/bold cyan]"
        )

    # ── Custom sequence → override template ────────────────────────────

    @on(TextArea.Changed, "#pd-custom-seq")
    def _custom_seq_changed(self, event: TextArea.Changed) -> None:
        """When the user pastes a custom sequence, auto-fill start=1 and
        end=len. Selected text within the TextArea will be used as the
        target region when designing (handled in _do_design)."""
        ta = self.query_one("#pd-custom-seq", TextArea)
        seq = ta.text.strip().upper().replace("\n", "").replace(" ", "")
        if seq and set(seq) <= set("ACGTRYWSMKBDHVN"):
            self.query_one("#pd-start", Input).value = "1"
            self.query_one("#pd-end", Input).value = str(len(seq))

    # ── Feature dropdown → fill start/end ──────────────────────────────────

    @on(Select.Changed, "#pd-feat")
    def _feat_selected(self, event: Select.Changed) -> None:
        val = event.value
        if not isinstance(val, str) or "-" not in val:
            return
        parts = val.split("-", 1)
        try:
            start = int(parts[0])
            end   = int(parts[1])
            self.query_one("#pd-start", Input).value = str(start + 1)
            self.query_one("#pd-end", Input).value = parts[1]
            # Wrap features have end < start; _feat_len handles that.
            total = len(self._template)
            feat_len = _feat_len(start, end, total)

            # Always set part name to the feature label
            feat_label = ""
            for f in self._feats:
                if f"{f['start']}-{f['end']}" == val:
                    feat_label = f.get("label", "")
                    self.query_one("#pd-part-name", Input).value = feat_label
                    break

            # Show feature length info + auto-adjust detection product range
            info = self.query_one("#pd-feat-info", Static)
            min_primeable = 50   # 2 × 25 bp primers — absolute minimum

            if feat_len < min_primeable:
                info.update(
                    f"  [bold]{feat_label}[/bold]  "
                    f"[red]{feat_len} bp — too short to prime "
                    f"(need ≥{min_primeable} bp for two ~25 bp primers)[/red]"
                )
            elif feat_len < 450:
                # Feature is smaller than the default 450-550 range —
                # auto-adjust to fit inside the feature.
                rec_min = max(min_primeable, feat_len - 100)
                rec_max = feat_len
                self.query_one("#pd-det-min", Input).value = str(rec_min)
                self.query_one("#pd-det-max", Input).value = str(rec_max)
                info.update(
                    f"  [bold]{feat_label}[/bold]  "
                    f"[yellow]{feat_len} bp[/yellow]  "
                    f"[dim]— adjusted detection range to "
                    f"{rec_min}–{rec_max} bp to fit inside feature[/dim]"
                )
            else:
                # Feature is large enough for the default range — reset to defaults
                self.query_one("#pd-det-min", Input).value = "450"
                self.query_one("#pd-det-max", Input).value = "550"
                info.update(
                    f"  [bold]{feat_label}[/bold]  "
                    f"[green]{feat_len} bp[/green]  "
                    f"[dim]— detection range 450–550 bp[/dim]"
                )
        except ValueError:
            pass

    # ── Primer library mark/unmark (m / ctrl+m) ──────────────────────────

    def on_key(self, event) -> None:
        """Handle m (mark/unmark cursor row) and ctrl+m (mark/unmark all)
        when the library table is focused."""
        try:
            t = self.query_one("#pd-lib-table", DataTable)
        except Exception:
            return
        if self.app.focused is not t:
            return

        primers = _load_primers()
        if event.key == "m":
            row = t.cursor_row
            if 0 <= row < len(primers):
                if row in self._lib_selected:
                    self._lib_selected.discard(row)
                else:
                    self._lib_selected.add(row)
                self._refresh_library_table()
            event.stop()
        elif event.key in ("M", "shift+m"):
            # Shift+M = mark all / unmark all (ctrl+m doesn't work in
            # terminals — Ctrl+M sends CR which is indistinguishable
            # from Enter)
            if len(self._lib_selected) == len(primers) and len(primers) > 0:
                self._lib_selected.clear()
            else:
                self._lib_selected = set(range(len(primers)))
            self._refresh_library_table()
            event.stop()
        elif event.key in ("S", "shift+s"):
            # Shift+S = cycle status: Designed → Ordered → Validated → Designed
            row = t.cursor_row
            if 0 <= row < len(primers):
                _CYCLE = ["Designed", "Ordered", "Validated"]
                cur = primers[row].get("status", "Designed")
                nxt = _CYCLE[(_CYCLE.index(cur) + 1) % 3] if cur in _CYCLE else _CYCLE[0]
                entries = _load_primers()
                if row < len(entries):
                    entries[row]["status"] = nxt
                    _save_primers(entries)
                    self._refresh_library_table()
                    name = primers[row].get("name", "?")
                    self.app.notify(f"{name}: {nxt}")
            event.stop()
        elif event.key == "delete":
            self._delete_marked_or_cursor()
            event.stop()

    def _delete_marked_or_cursor(self) -> None:
        """If primers are marked, confirm deletion of all marked. Otherwise
        confirm deletion of the cursor row."""
        primers = _load_primers()
        if self._lib_selected:
            count = len(self._lib_selected)
            names = [primers[i].get("name", "?") for i in sorted(self._lib_selected)
                     if 0 <= i < len(primers)]
            label = f"{count} marked primer{'s' if count != 1 else ''}"
        else:
            name = self._selected_primer_name()
            if name is None:
                self.app.notify("No primer selected.", severity="warning")
                return
            names = [name]
            label = f"primer {name!r}"

        def _on_confirm(result):
            if result is not True:
                return
            entries = _load_primers()
            name_set = set(names)
            entries = [e for e in entries if e.get("name") not in name_set]
            _save_primers(entries)
            self._lib_selected.clear()
            self._refresh_library_table()
            self.app.notify(f"Deleted {len(names)} primer{'s' if len(names) != 1 else ''}.")

        self.app.push_screen(
            LibraryDeleteConfirmModal(label, 0, ""),
            callback=_on_confirm,
        )

    # ── Helpers ────────────────────────────────────────────────────────────

    def _read_region(self) -> "tuple[int, int, str] | None":
        """Read and validate start/end/part-name from the inputs.
        Returns (start_0based, end, part_name) or None after notifying.

        `end <= start` is accepted as a wrap region on a circular plasmid;
        see _read_region_from for details.
        """
        return self._read_region_from(self._template)

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

    # ── Unified design button ─────────────────────────────────────────────

    @on(Button.Pressed, "#btn-pd-design")
    def _do_design(self, _) -> None:
        """Single Design button dispatches to the active primer mode."""
        mode_val = self._current_mode()

        # Resolve template: custom sequence overrides loaded plasmid.
        # If text is selected inside the TextArea, use ONLY the selection
        # as the template (and reset start/end to cover it). Otherwise use
        # the full TextArea text, or fall back to the loaded plasmid.
        ta = self.query_one("#pd-custom-seq", TextArea)
        raw_text = ta.text.strip().upper().replace("\n", "").replace(" ", "")
        selected = ta.selected_text.strip().upper().replace("\n", "").replace(" ", "")
        if selected and set(selected) <= set("ACGTRYWSMKBDHVN"):
            template = selected
            # Override start/end to the full selection
            self.query_one("#pd-start", Input).value = "1"
            self.query_one("#pd-end", Input).value = str(len(selected))
        elif raw_text and set(raw_text) <= set("ACGTRYWSMKBDHVN"):
            template = raw_text
        else:
            template = self._template
        if not template:
            self.app.notify("No sequence available. Load a plasmid or paste a custom sequence.", severity="error")
            return

        region = self._read_region_from(template)
        if region is None:
            return
        start, end, name = region
        self._det_result = None
        self._clo_result = None

        if mode_val == "detection":
            self._run_detection(template, start, end)
        elif mode_val == "cloning":
            self._run_cloning(template, start, end)
        elif mode_val == "goldenbraid":
            self._run_goldenbraid(template, start, end)
        elif mode_val == "generic":
            self._run_generic(template, start, end)

    def _read_region_from(self, template: str) -> "tuple[int, int, str] | None":
        """Like _read_region but uses the given template length for validation.

        `end <= start` is accepted as a wrap-around region on a circular
        plasmid — the primer-design helpers slice template[start:] + template[:end]
        and map primer positions back with modular arithmetic. We still
        reject negative start, end past template, and same-position empties.
        """
        try:
            start = int(self.query_one("#pd-start", Input).value) - 1
            end   = int(self.query_one("#pd-end", Input).value)
        except ValueError:
            self.app.notify("Enter valid start and end positions.", severity="error")
            return None
        total = len(template)
        if start < 0 or end < 0 or start >= total or end > total:
            self.app.notify(
                f"Invalid region: {start+1}–{end} (sequence is {total} bp).",
                severity="error")
            return None
        if end == start:
            self.app.notify(
                f"Region is empty (start and end both at {start+1}).",
                severity="error")
            return None
        name = self.query_one("#pd-part-name", Input).value.strip() or "primer"
        return start, end, name

    def _run_detection(self, template: str, start: int, end: int) -> None:
        try:
            p_min = int(self.query_one("#pd-det-min", Input).value)
            p_max = int(self.query_one("#pd-det-max", Input).value)
            tm    = float(self.query_one("#pd-det-tm", Input).value)
            plen  = int(self.query_one("#pd-det-len", Input).value)
        except ValueError:
            self.app.notify("Invalid detection parameters.", severity="error")
            return
        result = _design_detection_primers(
            template, start, end, product_min=p_min, product_max=p_max,
            target_tm=tm, primer_len=plen,
        )
        if "error" in result:
            self.query_one("#pd-results", Static).update(f"[red]{result['error']}[/red]")
            return
        self._det_result = result
        self._det_result["_type"] = "detection"
        self._show_result(result, "detection", "fwd_seq", "rev_seq")

    def _run_cloning(self, template: str, start: int, end: int) -> None:
        cust5 = self.query_one("#pd-cust5", Input).value.strip().upper()
        cust3 = self.query_one("#pd-cust3", Input).value.strip().upper()
        re5 = self.query_one("#pd-re5", Select).value
        re3 = self.query_one("#pd-re3", Select).value
        if cust5 and set(cust5) <= set("ACGTRYWSMKBDHVN"):
            site_5, name_5 = cust5, f"custom({cust5})"
        elif isinstance(re5, str) and re5 in _NEB_ENZYMES:
            site_5, name_5 = _NEB_ENZYMES[re5][0], re5
        else:
            self.app.notify("Select a 5' RE or enter a custom sequence.", severity="error")
            return
        if cust3 and set(cust3) <= set("ACGTRYWSMKBDHVN"):
            site_3, name_3 = cust3, f"custom({cust3})"
        elif isinstance(re3, str) and re3 in _NEB_ENZYMES:
            site_3, name_3 = _NEB_ENZYMES[re3][0], re3
        else:
            self.app.notify("Select a 3' RE or enter a custom sequence.", severity="error")
            return
        try:
            tm = float(self.query_one("#pd-clo-tm", Input).value)
        except ValueError:
            tm = 60.0
        result = _design_cloning_primers_raw(
            template, start, end, site_5, site_3, name_5, name_3, target_tm=tm,
        )
        if "error" in result:
            self.query_one("#pd-results", Static).update(f"[red]{result['error']}[/red]")
            return
        self._clo_result = result
        self._clo_result["_type"] = "cloning"
        self._show_result(result, "cloning", "fwd_full", "rev_full")

    def _run_goldenbraid(self, template: str, start: int, end: int) -> None:
        pt = self.query_one("#pd-gb-type", Select).value
        if not isinstance(pt, str) or pt not in _GB_POSITIONS:
            self.app.notify("Select a GB part type.", severity="error")
            return
        result = _design_gb_primers(template, start, end, pt)
        if "error" in result:
            self.query_one("#pd-results", Static).update(f"[red]{result['error']}[/red]")
            return
        self._clo_result = result
        self._clo_result["_type"] = "goldenbraid"
        self._show_result(result, "goldenbraid", "fwd_full", "rev_full")

    def _run_generic(self, template: str, start: int, end: int) -> None:
        try:
            tm = float(self.query_one("#pd-gen-tm", Input).value)
        except ValueError:
            tm = 60.0
        result = _design_generic_primers(template, start, end, target_tm=tm)
        if "error" in result:
            self.query_one("#pd-results", Static).update(f"[red]{result['error']}[/red]")
            return
        self._det_result = result
        self._det_result["_type"] = "generic"
        self._show_result(result, "generic", "fwd_seq", "rev_seq")

    # ── Save to primer library ─────────────────────────────────────────────

    @on(Button.Pressed, "#btn-pd-save")
    def _save_primers_btn(self, _) -> None:
        result = self._det_result or self._clo_result
        if result is None:
            self.app.notify("Design primers first.", severity="warning")
            return
        fwd_name = self.query_one("#pd-fwd-name", Input).value.strip()
        rev_name = self.query_one("#pd-rev-name", Input).value.strip()
        if not fwd_name or not rev_name:
            self.app.notify("Enter primer names before saving.", severity="error")
            return

        ptype = result.get("_type", "?")
        fwd_key = "fwd_seq" if ptype in ("detection", "generic") else "fwd_full"
        rev_key = "rev_seq" if ptype in ("detection", "generic") else "rev_full"

        # Source = plasmid name (not feature name). For generic mode,
        # use the source-ID input if filled.
        if ptype == "generic":
            source = self.query_one("#pd-gen-source", Input).value.strip()
            if not source:
                source = self._plasmid_name or "custom"
        else:
            source = self._plasmid_name or "custom"

        # Check for duplicate SEQUENCES (not names) already in the library
        entries = _load_primers()
        fwd_seq = result[fwd_key]
        rev_seq = result[rev_key]
        existing_seqs = {e.get("sequence", "").upper() for e in entries}
        dupes = []
        if fwd_seq.upper() in existing_seqs:
            dupes.append(fwd_name)
        if rev_seq.upper() in existing_seqs:
            dupes.append(rev_name)
        if dupes:
            self.app.notify(
                f"Duplicate sequence already in library: {', '.join(dupes)}. "
                f"Primer not saved — rename or modify the design.",
                severity="warning", timeout=8,
            )
            return

        import datetime
        today = datetime.date.today().isoformat()

        for pname, seq, tm, pos in [
            (fwd_name, fwd_seq, result["fwd_tm"], result["fwd_pos"]),
            (rev_name, rev_seq, result["rev_tm"], result["rev_pos"]),
        ]:
            entries = [e for e in entries if e.get("name") != pname]
            entries.insert(0, {
                "name":        pname,
                "sequence":    seq,
                "tm":          tm,
                "primer_type": ptype,
                "source":      source,
                "pos_start":   pos[0],
                "pos_end":     pos[1],
                "strand":      1 if pname.endswith("-F") else -1,
                "date":        today,
                "status":      "Designed",
            })
        _save_primers(entries)
        self._refresh_library_table()
        self.app.notify(f"Saved {fwd_name} + {rev_name} to primer library.")
        self._reset_for_new_design()

    # ── Add selected library primers as features ──────────────────────────

    @on(Button.Pressed, "#btn-pdlib-addmap")
    def _add_selected_to_map(self, _) -> None:
        """Add ALL multi-selected primers from the library as primer_bind
        features on the currently-loaded plasmid."""
        if not self._lib_selected:
            self.app.notify(
                "No primers marked. Press m to mark primers first.",
                severity="warning",
            )
            return
        rec = getattr(self.app, "_current_record", None)
        if rec is None:
            self.app.notify("No plasmid loaded.", severity="warning")
            return

        primers = _load_primers()
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from copy import deepcopy

        # Build a fresh record rather than mutating self.app._current_record:
        # the undo stack aliases the live record, and appending features in
        # place would silently corrupt prior undo snapshots.
        new_rec = SeqRecord(
            Seq(str(rec.seq)),
            id=rec.id, name=rec.name, description=rec.description,
            annotations=dict(rec.annotations),
        )
        for f in rec.features:
            new_rec.features.append(deepcopy(f))

        total = len(new_rec.seq)
        added = []
        for idx in sorted(self._lib_selected):
            if idx < 0 or idx >= len(primers):
                continue
            p = primers[idx]
            p_start = p.get("pos_start", 0)
            p_end   = p.get("pos_end", 0)
            strand  = p.get("strand", 1)
            name    = p.get("name", "primer")
            if p_end == p_start:
                continue
            # Wrap primer: pos_end < pos_start means the binding region
            # crosses the origin. Represent as a CompoundLocation so
            # downstream parsers (and GenBank exports) keep the two pieces.
            if p_end < p_start and 0 <= p_end and p_start < total:
                loc = CompoundLocation([
                    FeatureLocation(p_start, total, strand=strand),
                    FeatureLocation(0,       p_end, strand=strand),
                ])
            elif 0 <= p_start < p_end <= total:
                loc = FeatureLocation(p_start, p_end, strand=strand)
            else:
                # Positions out of range for the current plasmid — skip.
                continue

            # Don't duplicate: skip if a primer_bind with the same label
            # and position already exists on the record
            def _loc_matches(f_loc) -> bool:
                if isinstance(f_loc, CompoundLocation) and isinstance(loc, CompoundLocation):
                    return [(int(p.start), int(p.end)) for p in f_loc.parts] == \
                           [(int(p.start), int(p.end)) for p in loc.parts]
                if isinstance(f_loc, FeatureLocation) and isinstance(loc, FeatureLocation):
                    return int(f_loc.start) == p_start and int(f_loc.end) == p_end
                return False

            already = any(
                f.type == "primer_bind" and _loc_matches(f.location)
                for f in new_rec.features
            )
            if already:
                continue
            new_rec.features.append(SeqFeature(
                loc, type="primer_bind",
                qualifiers={"label": [name]},
            ))
            added.append(name)

        if not added:
            self.app.notify("Selected primers are already on the map.",
                            severity="information")
            return

        try:
            # clear_undo=False keeps the primer-add in the undo stack AND
            # preserves _source_path so Ctrl+S still targets the right file.
            self.app._push_undo()
            self.app._apply_record(new_rec, clear_undo=False)
            self.app._mark_dirty()
            lib = self.app.query_one("#library")
            lib.add_entry(new_rec)
        except Exception:
            _log.exception("Failed to add primer features to map")

        self.app.notify(
            f"Added {len(added)} primer{'s' if len(added) != 1 else ''} "
            f"as features: {', '.join(added)}"
        )

    # ── Primer library management ─────────────────────────────────────────

    def _selected_primer_name(self) -> "str | None":
        """Return the name of the currently-highlighted primer in the library
        table, or None if nothing is selected."""
        t = self.query_one("#pd-lib-table", DataTable)
        primers = _load_primers()
        if t.row_count == 0 or not (0 <= t.cursor_row < len(primers)):
            return None
        return primers[t.cursor_row].get("name")

    @on(Button.Pressed, "#btn-pdlib-rename")
    def _rename_primer(self, _) -> None:
        old_name = self._selected_primer_name()
        if old_name is None:
            self.app.notify("Highlight a primer to rename.", severity="warning")
            return

        def _on_result(new_name: "str | None") -> None:
            if new_name is None or new_name == old_name:
                return
            entries = _load_primers()
            # Check collision
            if any(e.get("name") == new_name for e in entries):
                self.app.notify(
                    f"A primer named {new_name!r} already exists.",
                    severity="error")
                return
            for e in entries:
                if e.get("name") == old_name:
                    e["name"] = new_name
                    break
            _save_primers(entries)
            self._refresh_library_table()
            self.app.notify(f"Renamed {old_name!r} → {new_name!r}")

        self.app.push_screen(
            RenamePlasmidModal(old_name, old_name),
            callback=_on_result,
        )

    @on(Button.Pressed, "#btn-pdlib-del")
    def _delete_primer(self, _) -> None:
        self._delete_marked_or_cursor()

    # ── Close ──────────────────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-pd-close")
    def _close(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Unsaved-changes quit dialog ────────────────────────────────────────────────

class UnsavedQuitModal(ModalScreen):
    """Shown when the user tries to quit with unsaved edits."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next button", show=False),
    ]

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


class PlasmidPickerModal(ModalScreen):
    """Scrollable plasmid-picker modal. Shows all entries from the library.
    Dismisses with the selected entry's id, or None on cancel.
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, current_id: "str | None" = None):
        super().__init__()
        self._current_id = current_id

    def compose(self) -> ComposeResult:
        with Vertical(id="pick-dlg"):
            yield Static(" Select plasmid from library ", id="pick-title")
            yield DataTable(id="pick-table", cursor_type="row",
                            zebra_stripes=True)
            with Horizontal(id="pick-btns"):
                yield Button("Select",  id="btn-pick-ok",     variant="primary")
                yield Button("Cancel",  id="btn-pick-cancel")

    def on_mount(self) -> None:
        t = self.query_one("#pick-table", DataTable)
        t.add_columns("Name", "ID", "Size", "Features")
        cursor = 0
        entries = _load_library()
        for i, e in enumerate(entries):
            t.add_row(
                Text(e.get("name", "?"), style="bold"),
                e.get("id", "?"),
                f"{e.get('size', 0):,} bp",
                f"{e.get('n_feats', 0)}",
                key=e.get("id"),
            )
            if self._current_id and e.get("id") == self._current_id:
                cursor = i
        if entries:
            t.move_cursor(row=cursor)
            t.focus()

    @on(Button.Pressed, "#btn-pick-ok")
    def _select(self, _):
        t = self.query_one("#pick-table", DataTable)
        if t.row_count == 0:
            self.dismiss(None)
            return
        row_keys = list(t.rows.keys())
        if 0 <= t.cursor_row < len(row_keys):
            self.dismiss(row_keys[t.cursor_row].value)
        else:
            self.dismiss(None)

    @on(DataTable.RowSelected, "#pick-table")
    def _row_selected(self, event):
        # Enter-key selection = same as clicking Select
        if event.row_key and event.row_key.value:
            self.dismiss(event.row_key.value)

    @on(Button.Pressed, "#btn-pick-cancel")
    def _cancel_btn(self, _):
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class RenamePlasmidModal(ModalScreen):
    """Prompt for a new name for a library entry.

    Tab cycles between the Input and Save/Cancel buttons.
    Dismisses with the new name (a non-empty string) or None on cancel.
    Input validation (non-empty, trimmed, collision check) lives in the
    app-side handler — the modal just collects a value.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

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
    """Generic delete-confirmation modal. Used by the plasmid library,
    primer library, and any future list that needs handslip protection.

    Default focus is on [No]. Tab cycles between [No] and [Yes, remove].
    Escape dismisses as False (cancel).
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next button", show=False),
    ]

    def __init__(self, name: str, size: int, entry_id: str):
        super().__init__()
        self.entry_name = name
        self.entry_size = size
        self.entry_id   = entry_id

    def compose(self) -> ComposeResult:
        size_str = f" ({self.entry_size:,} bp)" if self.entry_size > 0 else ""
        with Vertical(id="libdel-dlg"):
            yield Static(" Remove from library ", id="libdel-title")
            yield Static(
                f"  Remove [bold]{self.entry_name}[/bold]"
                f"{size_str} from the library?\n\n"
                f"  [dim]This cannot be undone from within the app.\n"
                f"  A backup (.bak) of the library file is kept.[/dim]",
                id="libdel-msg",
                markup=True,
            )
            with Horizontal(id="libdel-btns"):
                yield Button("No",           id="btn-libdel-no",  variant="default")
                yield Button("Yes, remove",  id="btn-libdel-yes", variant="error")

    def on_mount(self) -> None:
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
    _show_restr: bool = False
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

/* ── Export-GenBank modal ────────────────────────────────── */
ExportGenBankModal { align: center middle; }
#export-box {
    width: 72; height: auto;
    background: $surface; border: solid $primary; padding: 1 2;
}
#export-title { background: $primary; padding: 0 1; margin-bottom: 1; }
#export-box Label { color: $text-muted; margin-top: 1; }
#export-btns { height: 3; margin-top: 1; }
#export-btns Button { margin-right: 1; }
#export-status { height: 2; margin-top: 1; }

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

/* ── Plasmid picker modal ────────────────────────────────── */
PlasmidPickerModal { align: center middle; }
#pick-dlg {
    width: 80; height: 26;
    background: $surface; border: solid $primary; padding: 1 2;
}
#pick-title  { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#pick-table  { height: 1fr; }
#pick-btns   { height: 3; margin-top: 1; }
#pick-btns Button { margin-right: 1; min-width: 14; }

/* ── Feature picker modal ────────────────────────────────── */
PlasmidFeaturePickerModal { align: center middle; }
#featpick-dlg {
    width: 80; height: 26;
    background: $surface; border: solid $primary; padding: 1 2;
}
#featpick-title { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#featpick-table { height: 1fr; }
#featpick-btns  { height: 3; margin-top: 1; }
#featpick-btns Button { margin-right: 1; min-width: 14; }

/* ── Add-feature modal ───────────────────────────────────── */
AddFeatureModal { align: center middle; }
#addfeat-dlg {
    width: 82; height: 90%; max-height: 40;
    background: $surface; border: solid $primary; padding: 1 2;
}
#addfeat-title { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#addfeat-dlg Label { color: $text-muted; margin-top: 1; }
#addfeat-body       { height: 1fr; overflow-y: auto; }
#addfeat-row1       { height: 6; }
#addfeat-type-col   { width: 1fr; height: 6; margin-right: 2; }
#addfeat-strand-col { width: 30; height: 6; }
#addfeat-type-col Label, #addfeat-strand-col Label { margin-top: 0; height: 1; }
#addfeat-seq    { height: 6; min-height: 4; margin-top: 0; }
#addfeat-status { height: 2; margin-top: 1; }
#addfeat-btns   { height: 3; margin-top: 1; }
#addfeat-btns Button { margin-right: 1; }

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

/* ── Mutagenize modal ───────────────────────────────────── */
MutagenizeModal { align: center middle; }
#mut-box {
    width: 115; height: auto; max-height: 46;
    background: $surface; border: solid $accent; padding: 1 2;
}
#mut-title    { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#mut-box Label { color: $text-muted; margin-top: 1; }
#mut-src-map, #mut-src-lib, #mut-src-prot { height: auto; }
#mut-prot-aa  { height: 6; border: solid $primary-darken-2; }
#mut-prot-row { height: 3; margin-top: 1; }
#mut-prot-row Input { width: 2fr; margin-right: 1; }
#mut-prot-row Button { width: 1fr; }
#mut-cds-info { height: 1; margin: 0 0 1 0; }
#mut-codon-row { height: 3; margin: 1 0; }
#mut-codon-label { width: 4fr; padding: 1 1 0 1; }
#mut-codon-row Button { width: 1fr; }
#mut-row2     { height: 5; }
#mut-mut-col  { width: 3fr; padding-right: 1; }
#mut-btn-col  { width: 2fr; }
#btn-mut-design { width: 100%; }
#mut-preview {
    height: auto; max-height: 10;
    border: solid $primary-darken-2; padding: 0 1; margin-top: 1;
    overflow-y: auto; overflow-x: auto;
}
#mut-results {
    height: auto; max-height: 12;
    border: solid $primary-darken-2; padding: 0 1; margin-top: 1;
    overflow-y: auto;
}
#mut-btns     { height: 3; margin-top: 1; }
#mut-btns Button { margin-right: 1; }

/* ── AA picker sub-modal (from clicking an AA in the preview) ────────── */
AminoAcidPickerModal { align: center middle; }
#aa-pick-box {
    width: 60; height: auto; max-height: 30;
    background: $surface; border: solid $accent; padding: 1 2;
}
#aa-pick-title { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#aa-pick-box Label { color: $text-muted; margin: 0 0 1 0; }
#aa-pick-list  { height: auto; max-height: 20; border: solid $primary-darken-2; }
#aa-pick-btns  { height: 3; margin-top: 1; }
#aa-pick-btns Button { margin-right: 1; }

/* ── Species picker modal ───────────────────────────────── */
SpeciesPickerModal { align: center middle; }
#sp-box {
    width: 90; height: auto; max-height: 34;
    background: $surface; border: solid $accent; padding: 1 2;
}
#sp-title    { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#sp-box Label { color: $text-muted; margin-top: 1; }
#sp-list     { height: 12; border: solid $primary-darken-2; }
#sp-info     { height: 1; margin: 1 0; }
#sp-fetch-row { height: 3; margin-top: 1; }
#sp-fetch-row Input { width: 2fr; margin-right: 1; }
#sp-fetch-row Button { width: 1fr; }
#sp-btns     { height: 3; margin-top: 1; }
#sp-btns Button { margin-right: 1; }

/* ── Primer design screen (full-screen, Option B tabbed layout) ────── */
#pd-box {
    width: 100%; height: 1fr;
    background: $surface; padding: 0 1;
    overflow-y: auto;
}
#pd-title { background: $primary-darken-2; color: $text; padding: 0 1; }

/* Section wrapper: titled box around each logical group. Height hugs
   content (no flex-grow) so sections stack tightly without leftover rows. */
.pd-section {
    width: 100%;
    height: auto;
    border: round $primary-darken-2;
    padding: 0 1;
    margin: 1 0 0 0;
}
.pd-section-hdr {
    width: 100%; height: 1;
    color: $text;
    background: $primary-darken-2;
    text-style: bold;
    padding: 0 1;
    margin: 0 0 0 0;
}

/* ── Open-book split: workflow (left page) + library (right page).
   Left gets 3fr because it has more content (template + mode + params +
   results); right gets 2fr (library datatable). Mins keep each page
   usable on narrow terminals. ── */
#pd-book { width: 100%; height: 1fr; }
#pd-left-page {
    width: 3fr;
    min-width: 60;
    height: 1fr;
    padding-right: 1;
}
#pd-right-page {
    width: 2fr;
    min-width: 44;
    height: 1fr;
}

/* Sections stack vertically inside the left page, all full-width.
   Template and Mode get margin-top: 0 so the page starts compactly:
   Panel 1 moves up 1 row (closes the title→panel-1 gap), Panel 2 moves
   up 2 rows (same gap closed plus cascade from Panel 1 moving up),
   Panel 3 stays on its margin-top:1 and cascades up 2 rows too. */
#pd-template-section,
#pd-mode-section,
#pd-params-section { width: 100%; height: auto; }
#pd-template-section { margin-top: 0; }
#pd-mode-section     { margin-top: 0; }

/* Results section (right page, above library). Gets more height for
   the roomy 2-row name/button layout. */
#pd-results-section { width: 100%; height: auto; margin: 0 0 1 0; }

/* TEMPLATE */
#pd-src-row { height: 3; align: left middle; }
#pd-src-row Label { width: auto; padding: 0 1 0 0; content-align: center middle; }
#pd-src-row #pd-source { width: 22; }
#pd-plasmid-lbl { margin-left: 2; }
#pd-plasmid-name {
    width: auto; max-width: 40;
    padding: 0 1;
    content-align: left middle;
    color: $accent;
}
#pd-src-row #btn-pd-pickplasmid { min-width: 10; margin-left: 1; }

.pd-source-panel { height: auto; }
#pd-src-feature { height: 4; margin-top: 1; }
#pd-feat-col  { width: 3fr; padding-right: 1; }
#pd-start-col { width: 9;  padding-right: 1; }
#pd-end-col   { width: 9;  padding-right: 1; }
#pd-name-col  { width: 2fr; min-width: 18; }
#pd-feat      { width: 100%; }
#pd-start, #pd-end { width: 100%; }
#pd-part-name { width: 100%; }

#pd-src-custom { height: auto; margin-top: 1; }
#pd-custom-seq { height: 6; min-height: 4; }

#pd-feat-info { height: 1; margin-top: 1; }
#pd-wrap-hint { height: 1; margin-top: 1; }

/* MODE — stacked radio set (all 4 options visible) */
#pd-mode-radio { height: auto; padding: 0 1; background: transparent; border: none; }
#pd-mode-radio RadioButton { padding: 0 1; margin: 0; background: transparent; }

/* PARAMETERS — mode panels + docked Design button, one active at a time.
   height: auto so single-row panels (Detection, Generic) stay 3 rows and
   multi-row panels (Cloning = 3 inner rows, GB = 2) expand accordingly. */
.pd-mode-panel { height: auto; padding: 0; }
/* Each inner row inside a multi-row panel. */
.pd-mode-row { height: 3; align: left middle; padding: 0; }
/* Labels & inputs uniform across all rows. */
.pd-mode-panel Label, .pd-mode-row Label {
    width: auto; padding: 0 0 0 1; content-align: center middle;
}
.pd-mode-panel Input { width: 10; margin: 0 0; }
/* Single-row (Horizontal) panels also need an align rule. */
#pd-panel-det, #pd-panel-gen { height: 3; align: left middle; }

/* Detection: product min/max (3-4 digit bp); Tm, Len (2-3 digit). */
#pd-det-min, #pd-det-max { width: 10; }
#pd-det-tm, #pd-det-len  { width: 10; }

/* Cloning: RE Select fits longest label (~21 chars) + chrome;
   custom-RE Input roomy for 6-12 bp recognition sites + padding. */
#pd-re5, #pd-re3     { width: 28; }
#pd-cust5, #pd-cust3 { width: 18; }
#pd-clo-tm           { width: 10; }

/* Golden Braid: part-type Select (labels up to 30 chars); oh-info fills. */
#pd-gb-type    { width: 40; }
#pd-gb-oh-info { width: 1fr; content-align: left middle; padding-left: 2; }

/* Generic: source ID is free text. */
#pd-gen-tm     { width: 10; }
#pd-gen-source { width: 32; }

/* DESIGN button — lives inside section 3 now; just one row tall, centered.
   No external margins — the section's own margin-top handles separation. */
#pd-design-row { height: 3; align: center middle; margin: 1 0 0 0; }
#pd-design-row Button { min-width: 26; }

/* ── RESULTS section ── */
#pd-results { height: auto; min-height: 4; max-height: 14; padding: 0 1; }
/* Name inputs: full-width row so each fwd/rev box is ~2× longer than
   when it shared a row with the Save/Add-to-Map buttons. */
#pd-result-names { height: 3; align: left middle; margin-top: 1; }
#pd-result-names Input { width: 1fr; margin-right: 1; }
/* Action buttons now sit on their own row below the name inputs. */
#pd-result-actions { height: 3; align: right middle; margin-top: 0; }
#pd-result-actions Button { min-width: 18; margin-left: 1; }

/* ── PRIMER LIBRARY (right page) — header bar + DataTable ── */
/* Sits at the top of pd-right-page, no outer margin (the book's own gap
   between left and right pages is pd-left-page's padding-right). */
#pd-lib-hdr-row {
    height: 3;
    align: left middle;
    background: $accent-darken-2;
    padding: 0 1;
    margin-top: 0;
}
#pd-lib-hdr {
    width: 1fr;
    color: $text;
    text-style: bold;
    content-align: left middle;
}
#pd-lib-hdr-row Button { min-width: 10; margin-left: 1; }
#pd-lib-table { width: 100%; height: 1fr; min-height: 6; }

.pd-fld-lbl { width: auto; padding: 0 1 0 0; }
"""

    BINDINGS = [
        Binding("f",           "fetch",            "Fetch GenBank", show=True),
        Binding("o",           "open_file",        "Open file",     show=True),
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
    # Actions that are always allowed — even on non-default screens.
    # These include Screen-level actions that every screen needs (cancel,
    # focus_next, noop placeholders). The critical guard below checks
    # self.screen.id; only app-level (main-screen) shortcuts like f, o, a,
    # A, E, r, v, etc. get blocked.
    _ALWAYS_ALLOWED_ACTIONS: set[str] = {
        "cancel", "focus_next", "noop",
    }

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
        """Rebuild SeqRecord after an insert/replace, shifting feature coords precisely.

        Wrap features (origin-spanning CompoundLocations of the canonical
        `join(tail..total, 1..head)` form) and other compound locations
        are shifted per-part so the wrap structure survives the edit.
        Features consumed entirely by a replace are dropped — we never
        leave 1-bp ghost stubs behind.
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

        ins_len  = len(new_bases)
        del_len  = 0 if mode == "insert" else (e - s)
        delta    = ins_len - del_len
        new_len  = len(new_seq)

        def _shift_range(fs: int, fe: int):
            """Return post-edit (new_fs, new_fe) for a simple [fs, fe) range,
            or None if the range was consumed entirely by a replace."""
            if mode == "insert":
                if fe <= s:
                    new_fs, new_fe = fs, fe
                elif fs >= s:
                    new_fs, new_fe = fs + ins_len, fe + ins_len
                else:
                    new_fs, new_fe = fs, fe + ins_len
            else:  # replace [s, e)
                if fe <= s:
                    new_fs, new_fe = fs, fe
                elif fs >= e:
                    new_fs, new_fe = fs + delta, fe + delta
                elif fs <= s and fe >= e:
                    new_fs, new_fe = fs, fe + delta
                elif fs < s:
                    new_fs, new_fe = fs, s + ins_len
                else:
                    new_fs, new_fe = s, fe + delta
            new_fs = max(0, min(new_fs, new_len))
            new_fe = max(0, min(new_fe, new_len))
            if new_fe <= new_fs:
                return None
            return (new_fs, new_fe)

        new_record = SeqRecord(
            Seq(new_seq),
            id=self._current_record.id,
            name=self._current_record.name,
            description=self._current_record.description,
            annotations=dict(self._current_record.annotations),
        )

        for feat in self._current_record.features:
            loc = feat.location
            if isinstance(loc, CompoundLocation):
                new_parts = []
                for part in loc.parts:
                    shifted = _shift_range(int(part.start), int(part.end))
                    if shifted is None:
                        continue
                    n_fs, n_fe = shifted
                    new_parts.append(FeatureLocation(
                        n_fs, n_fe,
                        strand=getattr(part, "strand", None),
                    ))
                if not new_parts:
                    continue  # feature entirely consumed
                if len(new_parts) == 1:
                    new_loc = new_parts[0]
                else:
                    new_loc = CompoundLocation(
                        new_parts,
                        operator=getattr(loc, "operator", "join"),
                    )
            else:
                shifted = _shift_range(int(loc.start), int(loc.end))
                if shifted is None:
                    continue  # feature entirely consumed
                n_fs, n_fe = shifted
                new_loc = FeatureLocation(
                    n_fs, n_fe,
                    strand=getattr(loc, "strand", 1),
                )

            new_record.features.append(SeqFeature(
                new_loc,
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
        # Per-plasmid undo: switching plasmids stashes the old stacks under
        # the old record.id and restores this plasmid's own history if it
        # was edited before. See _stash_current_undo_and_load.
        self._undo_stack: list = []
        self._redo_stack: list = []
        self._stashed_undo_stacks: "dict[str, list]" = {}
        self._stashed_redo_stacks: "dict[str, list]" = {}
        self._stash_order: list[str] = []  # LRU
        self._MAX_PLASMIDS_WITH_UNDO = 10
        self._current_undo_key: "str | None" = None
        # Re-entrancy guard for pLannotate: spawning a second subprocess
        # while the first is still running wastes 5-30 s of CPU and risks
        # the stale-check discarding the newer result. See action_annotate_plasmid.
        self._plannotate_running: bool = False
        # Crash-recovery autosave: debounced so rapid edits coalesce into one
        # write. Cleared whenever the record is saved / marked clean.
        self._autosave_timer = None
        self._AUTOSAVE_DEBOUNCE_S = 3.0
        # Validate all user-data files before anything else. Corrupt files
        # are auto-restored from .bak if possible; the user is notified
        # either way so they know the state of their data.
        self._check_data_files()
        self._check_crash_recovery()
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

    # ── Crash-recovery autosave ────────────────────────────────────────────────

    def _autosave_path(self, record) -> "Path | None":
        """Return the autosave file path for `record`, or None if no id."""
        if record is None or not getattr(record, "id", ""):
            return None
        import re
        safe = re.sub(r'[^A-Za-z0-9._-]', '_', record.id)[:80]
        if not safe:
            return None
        return _CRASH_RECOVERY_DIR / f"{safe}.gb"

    def _schedule_autosave(self) -> None:
        """Debounce: restart the countdown to the next autosave write."""
        try:
            if self._autosave_timer is not None:
                self._autosave_timer.stop()
        except Exception:
            pass
        try:
            self._autosave_timer = self.set_timer(
                self._AUTOSAVE_DEBOUNCE_S, self._do_autosave
            )
        except Exception:
            # set_timer can fail during test teardown; autosave is best-effort.
            self._autosave_timer = None

    def _do_autosave(self) -> None:
        """Write the current record to its autosave file (atomic)."""
        if self._current_record is None or not self._unsaved:
            return
        path = self._autosave_path(self._current_record)
        if path is None:
            return
        try:
            import os
            import tempfile
            _CRASH_RECOVERY_DIR.mkdir(parents=True, exist_ok=True)
            text = _record_to_gb_text(self._current_record)
            fd, tmp = tempfile.mkstemp(
                prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
            )
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as fh:
                    fh.write(text)
                    fh.flush()
                    try:
                        os.fsync(fh.fileno())
                    except OSError:
                        pass
                os.replace(tmp, str(path))
                _log.info("Autosaved %s to %s (%d bp)",
                          self._current_record.name, path,
                          len(self._current_record.seq))
            except Exception:
                try:
                    os.unlink(tmp)
                except OSError:
                    pass
                raise
        except Exception:
            # Autosave is a safety net — never interrupt the user if it fails.
            _log.exception("Autosave to %s failed", path)

    def _clear_autosave(self, record) -> None:
        """Delete the autosave file for `record` (called after save / abandon)."""
        path = self._autosave_path(record)
        if path is None or not path.exists():
            return
        try:
            path.unlink()
            _log.info("Cleared autosave %s", path)
        except OSError:
            _log.exception("Failed to clear autosave %s", path)

    def _check_crash_recovery(self) -> None:
        """On startup, warn the user if leftover autosaves exist.

        A .gb file in `_CRASH_RECOVERY_DIR` means the previous session
        made unsaved edits and didn't cleanly save or abandon. The user
        can recover via File > Open or by inspecting the directory.
        """
        try:
            if not _CRASH_RECOVERY_DIR.exists():
                return
            leftovers = sorted(_CRASH_RECOVERY_DIR.glob("*.gb"))
        except OSError:
            _log.exception("Could not scan crash-recovery dir")
            return
        if not leftovers:
            return
        names = ", ".join(p.stem for p in leftovers[:3])
        if len(leftovers) > 3:
            names += f" (+{len(leftovers) - 3} more)"
        self.notify(
            f"Unsaved recovery files from a prior session: {names}. "
            f"Open via File > Open from {_CRASH_RECOVERY_DIR}",
            severity="warning", timeout=15,
        )

    @work(thread=True)
    def _seed_default_library(self) -> None:
        """Fetch MW463917.1 and pre-populate the library on first run."""
        try:
            record = fetch_genbank("MW463917.1")
            def _add():
                # If the user loaded or fetched a different plasmid while the
                # seed fetch was in flight, _apply_record would silently stomp
                # their record with the seed. Add the entry to the library so
                # they can pick it later, but skip the apply in that case.
                lib = self.query_one("#library", LibraryPanel)
                lib.add_entry(record)
                if self._current_record is None:
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

    def _stash_current_undo_and_load(self, new_key: "str | None") -> None:
        """Move the live undo/redo stacks into `_stashed_*_stacks[old_key]` and
        load `new_key`'s stashed stacks (empty if never seen).

        Called only from `_apply_record(clear_undo=True)` — the "switch
        plasmid" path. Lets the user flip between open plasmids without
        losing either one's edit history. LRU-capped at
        `_MAX_PLASMIDS_WITH_UNDO` so opening dozens of plasmids can't
        balloon memory.
        """
        if self._current_undo_key == new_key:
            return
        # Stash the outgoing plasmid's stacks if non-empty
        old_key = self._current_undo_key
        if old_key is not None and (self._undo_stack or self._redo_stack):
            self._stashed_undo_stacks[old_key] = list(self._undo_stack)
            self._stashed_redo_stacks[old_key] = list(self._redo_stack)
            if old_key in self._stash_order:
                self._stash_order.remove(old_key)
            self._stash_order.append(old_key)
            while len(self._stash_order) > self._MAX_PLASMIDS_WITH_UNDO:
                evict = self._stash_order.pop(0)
                self._stashed_undo_stacks.pop(evict, None)
                self._stashed_redo_stacks.pop(evict, None)
        # Load the incoming plasmid's stacks (empty list if never seen)
        if new_key is not None:
            self._undo_stack = self._stashed_undo_stacks.pop(new_key, [])
            self._redo_stack = self._stashed_redo_stacks.pop(new_key, [])
            if new_key in self._stash_order:
                self._stash_order.remove(new_key)
        else:
            self._undo_stack = []
            self._redo_stack = []
        self._current_undo_key = new_key

    def _push_undo(self) -> None:
        sp = self.query_one("#seq-panel", SequencePanel)
        if not sp._seq:
            return
        # Deep-copy the record so future in-place mutations of the live
        # _current_record cannot retroactively poison this snapshot. All
        # current edit paths build a fresh SeqRecord (via _rebuild_*), so
        # the copy is defensive — but cheap (~5 ms on a 50 kb plasmid)
        # and worth the safety margin against future contributors.
        from copy import deepcopy
        snapshot = (sp._seq, sp._cursor_pos, deepcopy(self._current_record))
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
        from copy import deepcopy
        sp = self.query_one("#seq-panel", SequencePanel)
        # Deep-copy current state to redo stack — same independence guarantee
        # as _push_undo.
        self._redo_stack.append(
            (sp._seq, sp._cursor_pos, deepcopy(self._current_record))
        )
        seq, cursor_pos, record = self._undo_stack.pop()
        self._apply_snapshot(seq, cursor_pos, record)
        remaining = len(self._undo_stack)
        self.notify(f"Undo  ({remaining} left)")

    def _action_redo(self) -> None:
        if not self._redo_stack:
            self.notify("Nothing to redo", severity="information")
            return
        from copy import deepcopy
        sp = self.query_one("#seq-panel", SequencePanel)
        self._undo_stack.append(
            (sp._seq, sp._cursor_pos, deepcopy(self._current_record))
        )
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
        except NoMatches:
            pass
        self._schedule_autosave()

    def _mark_clean(self) -> None:
        self._unsaved = False
        if self._current_record:
            n = len(self._current_record.seq)
            self.title = f"SpliceCraft — {self._current_record.name}  ({n:,} bp)"
        try:
            self.query_one("#library", LibraryPanel).set_dirty(False)
        except NoMatches:
            pass
        # Successful save / explicit clean → delete the recovery file
        self._clear_autosave(self._current_record)

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
                _log.exception("Save to %s failed", self._source_path)
                self.notify(f"Save failed: {exc}", severity="error")
                return False

        # Always update the library entry (add or overwrite)
        try:
            lib = self.query_one("#library", LibraryPanel)
            lib.add_entry(self._current_record)
        except Exception as exc:
            _log.exception("Library update failed during save")
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
            # User chose to discard unsaved work — clear the autosave so
            # next startup doesn't flag recovery for a file they abandoned.
            self._clear_autosave(self._current_record)
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

    def action_export_genbank(self) -> None:
        """Prompt for a path and write the current record as GenBank.

        Uses the record's `_tui_source` path as the default if present, or
        `{record.name}.gb` in the cwd otherwise. Export is round-trip
        verified before any file is touched (see `_export_genbank_to_path`).
        """
        if self._current_record is None:
            self.notify("No plasmid loaded.", severity="warning")
            return
        # Default: prefer the record's source file if it was opened from
        # disk and was already .gb; fall back to a file beside the cwd.
        src = getattr(self._current_record, "_tui_source", None) \
              or getattr(self, "_source_path", None)
        if src and str(src).lower().endswith((".gb", ".gbk")):
            default = str(src)
        else:
            name = self._current_record.name or "plasmid"
            default = f"{name}.gb"

        def _on_done(result):
            if result is None:
                return
            path = result.get("path", "?")
            bp = result.get("bp", 0)
            feats = result.get("features", 0)
            self.notify(
                f"Exported {feats} features / {bp} bp → {path}",
                timeout=8,
            )

        self.push_screen(
            ExportGenBankModal(self._current_record, default_path=default),
            callback=_on_done,
        )

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

    def _apply_record(self, record, *, clear_undo: bool = True) -> None:
        """Load a SeqRecord into all panels.

        `clear_undo=True` (default) is for fresh loads (fetch, file open,
        library pick) — it **stashes** the outgoing plasmid's undo/redo
        stacks under its record.id and restores the incoming plasmid's
        stacks if it was edited before (see `_stash_current_undo_and_load`).
        Users flipping between open plasmids keep per-plasmid history so
        Ctrl+Z never yanks you to an unrelated edit, but switching back to
        a previously-edited plasmid resurrects its history. Also clears
        `_source_path` so Ctrl+S can't accidentally overwrite the old file.

        `clear_undo=False` is for in-place record changes (pLannotate merge,
        primer-add) — the stacks stay intact and the edit remains undo-able.
        """
        if record is None:
            return
        if clear_undo:
            new_key = record.id if record.id else None
            self._stash_current_undo_and_load(new_key)
            self._source_path = None   # caller sets this if it came from a file
        self._current_record = record

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
        except NoMatches:
            pass
        self._mark_clean()
        self.notify(
            f"Loaded {record.name}  ({len(record.seq):,} bp, "
            f"{len(pm._feats)} features, {len(self._restr_cache)} restriction sites)"
        )

        # Heads-up for non-fatal import oddities that the user should know
        # about: skipped features (UnknownPosition), compound locations
        # flattened to outer bounds, and linear topology auto-detection.
        n_skip = getattr(pm, "_n_skipped", 0)
        if n_skip:
            self.notify(
                f"⚠ Skipped {n_skip} feature(s) with unknown coordinates — "
                f"see {_LOG_PATH} for details.",
                severity="warning", timeout=8,
            )
        n_flat = getattr(pm, "_n_flattened", 0)
        if n_flat:
            self.notify(
                f"⚠ {n_flat} feature(s) have joined coordinates (e.g. exons); "
                f"rendered as outer-bounds span.",
                severity="warning", timeout=8,
            )
        n_clamp = getattr(pm, "_n_clamped", 0)
        if n_clamp:
            self.notify(
                f"⚠ {n_clamp} feature(s) had coordinates outside the "
                f"sequence length — clamped to fit.",
                severity="warning", timeout=8,
            )
        topology = (record.annotations or {}).get("topology", "").lower()
        if topology == "linear" and pm._map_mode != "linear":
            pm._map_mode = "linear"
            self.notify(
                "File declares linear topology — switched to linear view "
                "(press 'v' to toggle).",
                severity="information", timeout=8,
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
        if event.feat_dict is not None and event.bp >= 0:
            # Cursor lands on the clicked base (already inside the feature via
            # _bp_in); feature span becomes the copyable selection.
            seq_pnl.select_feature_range(event.feat_dict, cursor_bp=event.bp)
        else:
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
            _log.exception("Library load failed for entry %r",
                           event.entry.get("name", "?"))
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
        except NoMatches:
            return
        if pm.selected_idx < 0:
            return   # nothing to clear
        pm.selected_idx = -1
        pm.refresh()
        try:
            self.query_one("#sidebar", FeatureSidebar).show_detail(None)
        except NoMatches:
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
            except NoMatches:
                pm = None
            if pm is not None:
                # record.name is in the cache key, but nuke it explicitly
                # for belt-and-braces (in case future refactors drop the
                # name from the key).
                pm._draw_cache = None
                pm.refresh()
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
        if name == "Mutagenize":
            self.action_open_mutagenize()
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
                ("Open file (.gb / .dna)", "open_file"),
                ("Fetch from NCBI", "fetch"),
                ("---",             None),
                ("Add to Library",  "add_to_library"),
                ("Save",            "save"),
                ("Export as GenBank (.gb)...", "export_genbank"),
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
        """Open the AddFeatureModal. If the sequence panel has an active
        cursor, the Insert button is enabled; otherwise only Save is."""
        sp = None
        cursor_pos = -1
        try:
            sp = self.query_one("#seq-panel", SequencePanel)
            cursor_pos = getattr(sp, "_cursor_pos", -1)
        except NoMatches:
            cursor_pos = -1
        have_cursor = (self._current_record is not None and cursor_pos >= 0)
        self.push_screen(
            AddFeatureModal(prefill=None, have_cursor=have_cursor),
            callback=self._add_feature_result,
        )

    def _add_feature_result(self, result) -> None:
        """Callback for AddFeatureModal. `result` is either None (cancel) or
        ``{"action": "save"|"insert", "entry": {...}}``."""
        if not result:
            return
        action = result.get("action")
        entry  = result.get("entry") or {}
        if action == "save":
            entries = _load_features()
            # De-dup on (name, feature_type); latest write wins.
            key = (entry.get("name"), entry.get("feature_type"))
            entries = [e for e in entries
                       if (e.get("name"), e.get("feature_type")) != key]
            entries.append(entry)
            try:
                _save_features(entries)
            except (OSError, ValueError) as exc:
                _log.exception("Failed to save feature to library")
                self.notify(f"Save failed: {exc}", severity="error")
                return
            self.notify(f"Saved '{entry.get('name')}' to feature library.")
            return
        if action == "insert":
            try:
                self._insert_feature_at_cursor(entry)
            except (ValueError, RuntimeError) as exc:
                _log.exception("Failed to insert feature")
                self.notify(f"Insert failed: {exc}", severity="error")

    def _insert_feature_at_cursor(self, entry: dict) -> None:
        """Insert the feature's DNA at the current sequence-panel cursor,
        shift existing feature coords via `_rebuild_record_with_edit`, and
        append a new SeqFeature spanning the inserted region.

        Reverse-strand entries have their `sequence` interpreted as the 5'→3'
        of the feature as read, so the genomic bases inserted are the RC.
        """
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        if self._current_record is None:
            raise RuntimeError("Load a plasmid first.")
        sp = self.query_one("#seq-panel", SequencePanel)
        pm = self.query_one("#plasmid-map", PlasmidMap)
        pos = getattr(sp, "_cursor_pos", -1)
        if pos < 0:
            raise RuntimeError(
                "Click on the sequence to place a cursor before inserting.")
        strand = -1 if entry.get("strand") == -1 else 1
        feat_seq = (entry.get("sequence") or "").upper()
        if not feat_seq:
            raise ValueError("Feature has no sequence.")
        # Genomic bases inserted are the forward-strand; reverse-strand
        # entries store the revcomp, so flip them back.
        if strand == -1:
            genomic = _rc(feat_seq)
        else:
            genomic = feat_seq

        old_seq = str(self._current_record.seq)
        new_seq = old_seq[:pos] + genomic + old_seq[pos:]
        self._push_undo()
        new_record = self._rebuild_record_with_edit(
            new_seq, "insert", pos, pos, genomic,
        )
        new_feat = SeqFeature(
            FeatureLocation(pos, pos + len(genomic), strand=strand),
            type=entry.get("feature_type") or "misc_feature",
            qualifiers={k: list(v) if isinstance(v, (list, tuple)) else [v]
                        for k, v in (entry.get("qualifiers") or {}).items()},
        )
        label = entry.get("name") or ""
        if label and "label" not in new_feat.qualifiers:
            new_feat.qualifiers["label"] = [label]
        new_record.features.append(new_feat)

        self._current_record = new_record
        pm.load_record(new_record)
        self._restr_cache = _scan_restriction_sites(
            new_seq,
            min_recognition_len=self._restr_min_len,
            unique_only=self._restr_unique_only,
        )
        displayed = self._restr_cache if self._show_restr else []
        pm._restr_feats = displayed
        pm.refresh()
        self.query_one("#sidebar", FeatureSidebar).populate(pm._feats)
        sp.update_seq(new_seq, pm._feats + displayed)
        sp._cursor_pos = pos + len(genomic)
        sp._user_sel   = None
        sp._refresh_view()
        self._mark_dirty()
        self.notify(
            f"Inserted '{label or entry.get('feature_type')}' "
            f"({len(genomic)} bp) at {pos + 1}."
        )

    # ── pLannotate annotation ──────────────────────────────────────────────────

    def action_annotate_plasmid(self) -> None:
        """Run pLannotate on the currently-loaded record (shortcut: Shift+A)."""
        if getattr(self, "_plannotate_running", False):
            self.notify(
                "pLannotate is already running — wait for the current run "
                "to finish.",
                severity="information",
            )
            return
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
        self._plannotate_running = True
        self._run_plannotate_worker(self._current_record)

    @work(thread=True)
    def _run_plannotate_worker(self, record) -> None:
        """Background worker: runs pLannotate subprocess, merges, applies.
        Errors are logged and surfaced to the UI via notify(); nothing raw
        reaches the user. The `_plannotate_running` re-entry flag is
        cleared unconditionally in `_finally` so a crashed run does not
        lock the user out of future annotation attempts."""
        merged = None
        err = None
        try:
            try:
                annotated = _run_plannotate(record)
                merged    = _merge_plannotate_features(record, annotated)
            except PlannotateError as exc:
                _log.info("pLannotate: %s", exc)
                err = ("error", exc.user_msg)
            except Exception as exc:
                _log.exception("pLannotate worker crashed")
                err = ("crash", str(exc))
        finally:
            def _finally():
                self._plannotate_running = False
                if err is not None:
                    kind, msg = err
                    if kind == "error":
                        self.notify(msg, severity="error", timeout=10)
                    else:
                        self.notify(f"pLannotate crashed: {msg}",
                                    severity="error", timeout=10)
                    return
                # Guard against races: if the user loaded a different plasmid
                # while pLannotate was running, silently applying the merged
                # OLD record would clobber their newer work.
                if self._current_record is not record:
                    self.notify(
                        "pLannotate finished, but you've loaded a different "
                        "plasmid in the meantime — discarding annotation result.",
                        severity="warning", timeout=8,
                    )
                    return
                n_added = getattr(merged, "_plannotate_added", 0)
                if n_added == 0:
                    self.notify(
                        "pLannotate found no new features (all hits duplicated "
                        "existing annotations).",
                        severity="information",
                    )
                    return
                self._push_undo()          # annotation is undo-able
                self._apply_record(merged, clear_undo=False)
                # Mark dirty AFTER _apply_record (which calls _mark_clean
                # internally) so the user gets prompted on quit and sees the
                # * in the title.
                self._mark_dirty()
                self.notify(
                    f"Added {n_added} pLannotate feature"
                    f"{'s' if n_added != 1 else ''}. "
                    "Press 'a' to save to library.",
                    timeout=6,
                )
            self.call_from_thread(_finally)

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
        except NoMatches:
            pass
        self.push_screen(PrimerDesignScreen(seq, feats, name))

    def action_open_mutagenize(self) -> None:
        """Open the SOE-PCR site-directed mutagenesis primer designer.
        Requires a loaded plasmid with at least one CDS feature."""
        rec = self._current_record
        if rec is None:
            self.notify("Load a plasmid first (press 'f' or 'o').",
                        severity="warning")
            return
        seq = str(rec.seq)
        feats: list = []
        try:
            feats = self.query_one("#plasmid-map", PlasmidMap)._feats
        except NoMatches:
            pass
        if not any(f.get("type") in ("CDS", "gene") for f in feats):
            self.notify(
                "No CDS features on this plasmid — nothing to mutagenize.",
                severity="warning",
            )
            return
        self.push_screen(MutagenizeModal(seq, feats, rec.name or ""))

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
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    # Handle --version / -V without loading the TUI
    if arg in ("--version", "-V"):
        print(f"splicecraft {__version__}")
        return
    if arg in ("--help", "-h"):
        print(
            f"splicecraft {__version__}\n"
            "Usage: splicecraft [ACCESSION | FILE.gb]\n\n"
            "  splicecraft             # empty canvas\n"
            "  splicecraft L09137      # fetch pUC19 from NCBI\n"
            "  splicecraft my.gb       # open a local GenBank file\n\n"
            "Data files (library, parts, primers) live in:\n"
            f"  {_DATA_DIR}\n"
            "Override with $SPLICECRAFT_DATA_DIR."
        )
        return
    if len(sys.argv) > 2:
        print(
            f"splicecraft takes at most one argument (got {len(sys.argv) - 1}: "
            f"{' '.join(sys.argv[1:])}). Pass a single accession or file.",
            file=sys.stderr,
        )
        sys.exit(2)
    _log_startup_banner()
    app = PlasmidApp()

    if arg:
        looks_like_file = arg.lower().endswith((".gb", ".gbk", ".genbank", ".dna"))
        if Path(arg).exists():
            try:
                record = load_genbank(arg)
                record._tui_source = str(Path(arg).resolve())
            except Exception as exc:
                _log.exception("Failed to load %s", arg)
                print(f"Could not load {arg!r}: {exc}", file=sys.stderr)
                sys.exit(1)
        elif looks_like_file:
            # Clear "file not found" message instead of confusingly trying
            # to fetch a file-looking string from NCBI.
            print(f"File not found: {arg}", file=sys.stderr)
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
