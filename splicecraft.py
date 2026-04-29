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
from datetime import date as _date
from io import StringIO
from logging.handlers import RotatingFileHandler
from pathlib import Path

__version__ = "0.4.2"

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
    Button, DataTable, DirectoryTree, Footer, Header, Input, Label, ListItem,
    ListView, RadioButton, RadioSet, Select, Static, TextArea,
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


def _atomic_write_text(path: Path, text: str, encoding: str = "utf-8") -> None:
    """Atomically write *text* to *path* via ``tempfile`` + ``os.replace``.

    Guarantees: a concurrent crash leaves either the previous file intact
    or the new file in place — never a partial write. Callers that need
    a ``.bak`` should use :func:`_safe_save_json` instead (it layers the
    envelope, shrink-guard, and schema handling on top of this).
    """
    import os
    import tempfile
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
    )
    try:
        with os.fdopen(fd, "w", encoding=encoding) as fh:
            fh.write(text)
            fh.flush()
            try:
                os.fsync(fh.fileno())
            except OSError:
                pass
        os.replace(tmp, str(path))
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


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
                except (json.JSONDecodeError, ValueError):
                    # Corrupt / unreadable — shrink guard will treat as 0.
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
    # Keep the active collection's plasmids in sync with the library — every
    # add / remove / rename on the panel feeds through here, so a single
    # mirror call covers all CRUD without changing call sites.
    _sync_active_collection_plasmids(entries)


# ── Plasmid collections ────────────────────────────────────────────────────────
#
# A "collection" is a named snapshot of a plasmid library — the user can keep
# several themed sets (e.g. "yeast project", "E. coli toolkit", "MoClo plant")
# and switch between them. Switching loads the snapshot into `plasmid_library.json`
# wholesale, so the rest of the app keeps working off the single live library.
#
# On-disk shape (envelope schema v1):
#   {"_schema_version": 1, "entries": [
#       {"name": "...", "description": "...", "plasmids": [<library entry>, ...]},
#       ...
#   ]}
# Each `plasmids` entry mirrors a library row exactly (id, name, gb_text, ...),
# so saving = `_load_library()`; loading = `_save_library(plasmids)`.

_COLLECTIONS_FILE = _DATA_DIR / "collections.json"
_collections_cache: "list | None" = None

def _load_collections() -> list[dict]:
    """Return a deepcopy of the collections list so callers can mutate
    entries (rename, edit plasmids list) without poisoning the in-memory
    cache. Matches the `_load_features` contract documented in CLAUDE.md.
    """
    global _collections_cache
    if _collections_cache is None:
        entries, warning = _safe_load_json(_COLLECTIONS_FILE, "Plasmid collections")
        if warning:
            _log.warning(warning)
        _collections_cache = [e for e in entries if isinstance(e, dict)]
    from copy import deepcopy
    return deepcopy(_collections_cache)

def _save_collections(entries: list[dict]) -> None:
    global _collections_cache
    _safe_save_json(_COLLECTIONS_FILE, entries, "Plasmid collections")
    from copy import deepcopy
    _collections_cache = deepcopy(entries)


# Active collection — which named collection is the panel currently showing.
# Stored in settings.json (key: "active_collection") so it persists across
# launches. Collection identity is the user-facing name; names are unique
# per the modal's dup-name guard.

_DEFAULT_COLLECTION_NAME = "Main Collection"


def _get_active_collection_name() -> "str | None":
    val = _get_setting("active_collection", None)
    return val if isinstance(val, str) and val else None


def _set_active_collection_name(name: "str | None") -> None:
    """Persist (or clear) the active-collection pointer."""
    _set_setting("active_collection", name or "")


def _find_collection(name: str) -> "dict | None":
    for c in _load_collections():
        if c.get("name") == name:
            return c
    return None


def _collection_name_taken(name: str) -> bool:
    """Dup-name guard for create / rename. Pure check, no side effects."""
    return _find_collection(name) is not None


def _ensure_default_collection() -> None:
    """Idempotent: guarantee at least one collection exists.

    First-run users have a non-empty `plasmid_library.json` but no
    `collections.json` — migrate by wrapping their existing plasmids in a
    "Main Collection" and marking it active. Empty-library first-runs just
    get an empty Main Collection so the panel always has something to show.
    """
    colls = _load_collections()
    if colls:
        if not _get_active_collection_name():
            first = colls[0].get("name")
            if first:
                _set_active_collection_name(first)
        return
    plasmids = _load_library()
    _save_collections([{
        "name":        _DEFAULT_COLLECTION_NAME,
        "description": "Default collection",
        "plasmids":    plasmids,
        "saved":       _date.today().isoformat(),
    }])
    _set_active_collection_name(_DEFAULT_COLLECTION_NAME)


def _sync_active_collection_plasmids(entries: list[dict]) -> None:
    """Mirror the live library's contents into the active collection so the
    on-disk collection record never drifts from the panel's view.

    Silent no-op if no collection is active, or if the active name has
    been deleted (e.g. user removed it via the manager) — the next
    explicit Load/Save will re-establish a target.
    """
    name = _get_active_collection_name()
    if not name:
        return
    colls = _load_collections()
    for c in colls:
        if c.get("name") == name:
            c["plasmids"] = [dict(e) for e in entries if isinstance(e, dict)]
            _save_collections(colls)
            return


def _restore_library_from_active_collection() -> None:
    """Refresh `plasmid_library.json` with the active collection's plasmids
    so the panel renders what the user expects after a restart or after
    edits made in another session.

    Called once during app startup (after `_ensure_default_collection`).
    Silent no-op if no collection is active or the active one was deleted.
    Bypasses `_save_library`'s mirror — the collection is the source.
    """
    global _library_cache
    name = _get_active_collection_name()
    if not name:
        return
    coll = _find_collection(name)
    if coll is None:
        return
    plasmids = [dict(p) for p in (coll.get("plasmids") or [])
                if isinstance(p, dict)]
    _safe_save_json(_LIBRARY_FILE, plasmids, "Plasmid library")
    _library_cache = list(plasmids)


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
                # the reversed orientation displayed on the forward strand.
                # Symmetry: a forward-strand cut at offset `c` from the
                # recognition's 5' end appears on a reverse-bound site at
                # offset `site_len - c` from the bar's left edge.
                rev_cut_col = site_len - fwd_cut
                _top_cut_bp = (p + site_len - rev_cut) % n   # top-strand cut in fwd coords
                _top_cut_outside = ((_top_cut_bp - p) % n) >= site_len
                _cc  = rev_cut_col if 0 <= rev_cut_col < site_len else None
                _ext = _top_cut_bp if _top_cut_outside else None
                _emit_resite(hits, p, site_len, -1, color, name, _cc, _ext)
                # Bottom-strand cut (enzyme's fwd_cut mapped to fwd coords)
                cut_bp = (p + site_len - fwd_cut) % n if n > 0 else 0
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
        # triangle pointing toward the DNA row. strand==0 renders as an
        # arrowless bar; strand==2 renders as a double-headed bar.
        if bar_len == 1:
            bar_str = "▲" if is_below_dna else "▼"
        elif strand == 0:
            bar_str = "▒" * bar_len
        elif strand == 2:
            head = "◀" if starts_here else "▒"
            tail = "▶" if ends_here   else "▒"
            middle = "▒" * max(0, bar_len - 2)
            bar_str = head + middle + tail
        elif strand >= 1:
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


# Per-(seq_id, feats_id, line_width) cache for chunk decomposition. Keyed
# identically to _BUILD_SEQ_CACHE plus line_width (which only changes on
# terminal resize). Independent of cursor/selection/show_connectors —
# show_connectors only changes the row-multiplier rpg, applied at lookup
# time via prefix-sum arithmetic. Holds:
#   chunks       — list of (chunk_start, chunk_end, lane_groups, above, below)
#                  where above/below are lane-pair counts (NOT row counts).
#   prefix_dna2  — list[int]: total rows up to chunk i with rpg=2 (so cursor
#                  moves on a 200 kb plasmid don't re-run _feats_in_chunk
#                  for 1500 chunks every keystroke).
#   prefix_lanes — list[int]: total lane pairs (above+below) up to chunk i.
# Total rows before chunk i for arbitrary rpg:
#   prefix_dna2[i] + (rpg - 2) * prefix_lanes[i]
_CHUNK_LAYOUT_CACHE: dict = {}


# Per-(seq_id, feats_id, line_width, show_connectors) cache of pre-rendered
# Rich Text per chunk, ASSUMING NO OVERLAYS (no cursor / selection / RE
# highlight). Cursor moves only re-render the chunk under the cursor; the
# other ~1500 chunks of a 200 kb BAC are reused from this cache. Profile-
# driven: `_render_feature_row_pair` was 78% of cursor-move time on a 150 kb
# plasmid; lane art is fully deterministic per chunk so caching it
# eliminates that hot path. Stored as `dict[key, list[Text|None]]`.
_CHUNK_STATIC_CACHE: dict = {}


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


def _build_seq_inputs(seq: str, feats: list[dict]) -> tuple[list[str], list[dict]]:
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


def _chunk_layout(seq: str, feats: list[dict], line_width: int):
    """Return per-chunk lane decomposition + prefix sums, cached by id.

    Used by both `_bp_to_content_row` (cursor → row index) and
    `_build_seq_text` (chunk-by-chunk render). Computing this once per
    (seq, feats, line_width) instead of per cursor move turns BAC/cosmid
    cursor scrolling from ~120 ms/keystroke into low single digits — the
    `_feats_in_chunk` + `_chunk_lane_groups` calls per chunk dominate
    otherwise (200 kb / line_width 127 = ~1.5k chunks × N features).

    Returns ``(chunks, prefix_dna2, prefix_lanes)`` where:
      - ``chunks[i] = (chunk_start, chunk_end, lane_groups, above_pairs, below_pairs)``
      - ``prefix_dna2[i]`` = total rows from row 0 up to (not including) chunk i
        assuming rpg = 2 (no connector row).
      - ``prefix_lanes[i]`` = total lane pairs (above+below) up to chunk i.
    Total rows before chunk i for any rpg:
      ``prefix_dna2[i] + (rpg - 2) * prefix_lanes[i]``
    """
    n = len(seq)
    key = (id(seq), id(feats), line_width, n, len(feats))
    hit = _CHUNK_LAYOUT_CACHE.get(key)
    if hit is not None:
        return hit

    # Reuse _build_seq_inputs' filtered+sorted annot_feats so cache hits
    # there propagate here (and so the two paths see identical lane order).
    _, annot_feats = _build_seq_inputs(seq, feats)

    chunks: list = []
    prefix_dna2  = [0]   # rows with rpg=2; index 0 = before chunk 0
    prefix_lanes = [0]
    if n == 0 or line_width <= 0:
        layout = (chunks, prefix_dna2, prefix_lanes)
        _CHUNK_LAYOUT_CACHE[key] = layout
        return layout

    for chunk_start in range(0, n, line_width):
        chunk_end   = min(chunk_start + line_width, n)
        chunk_feats = _feats_in_chunk(annot_feats, chunk_start, chunk_end, n)
        groups      = _chunk_lane_groups(chunk_feats, chunk_start, chunk_end)
        re_above, onebp_above, reg_above, reg_below, onebp_below, re_below = groups
        above_pairs = len(re_above) + len(onebp_above) + len(reg_above)
        below_pairs = len(re_below) + len(onebp_below) + len(reg_below)
        chunks.append((chunk_start, chunk_end, groups, above_pairs, below_pairs))
        # Per-chunk rows with rpg=2: above*2 + 2 (DNA pair) + below*2
        prefix_dna2.append(prefix_dna2[-1] + above_pairs * 2 + 2 + below_pairs * 2)
        prefix_lanes.append(prefix_lanes[-1] + above_pairs + below_pairs)

    layout = (chunks, prefix_dna2, prefix_lanes)
    if len(_CHUNK_LAYOUT_CACHE) >= 4:
        _CHUNK_LAYOUT_CACHE.clear()
    _CHUNK_LAYOUT_CACHE[key] = layout
    return layout


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

    # Cached chunk decomposition — eliminates per-chunk _feats_in_chunk +
    # _chunk_lane_groups recomputation on every cursor-move re-render. On a
    # 200 kb plasmid that's ~1500 chunks × N features → milliseconds saved
    # per keystroke when scrolling through cosmid/BAC-scale records.
    chunks_layout, _pf_dna2, _pf_lanes = _chunk_layout(seq, feats, line_width)

    # Per-chunk static-render cache. The first render fills it; subsequent
    # cursor moves only re-render the chunk under the cursor and reuse the
    # other ~N pre-rendered Text objects. Without this, the BAC-scale
    # `_render_feature_row_pair` cost (~78 % of render time) recurs every
    # keystroke. Cache stays valid for cursor/selection changes; it
    # invalidates only on (seq, feats, line_width, show_connectors) change.
    static_key   = (id(seq), id(feats), line_width, show_connectors)
    static_cache = _CHUNK_STATIC_CACHE.get(static_key)
    if static_cache is None or len(static_cache) != len(chunks_layout):
        static_cache = [None] * len(chunks_layout)
        _CHUNK_STATIC_CACHE[static_key] = static_cache
        if len(_CHUNK_STATIC_CACHE) > 4:
            _CHUNK_STATIC_CACHE.pop(next(iter(_CHUNK_STATIC_CACHE)))

    for i, (chunk_start, chunk_end, groups, _ab_pairs, _be_pairs) in enumerate(chunks_layout):
        # If any overlay touches this chunk, render fresh. Otherwise reuse
        # (or first-time-populate) the static cache.
        chunk_has_overlay = (
            (usr_s < chunk_end and usr_e > chunk_start)
            or (sel_s < chunk_end and sel_e > chunk_start)
            or (reh_s < chunk_end and reh_e > chunk_start)
            or (chunk_start <= cursor_pos < chunk_end)
        )

        if chunk_has_overlay:
            _render_chunk(result, chunk_start, chunk_end, groups, styles,
                          num_w, seq_upper, show_connectors,
                          sel_s, sel_e, usr_s, usr_e,
                          reh_s, reh_e, reh_color, cursor_pos)
        else:
            cached = static_cache[i]
            if cached is None:
                cached = Text(no_wrap=True, overflow="crop")
                _render_chunk(cached, chunk_start, chunk_end, groups, styles,
                              num_w, seq_upper, show_connectors,
                              -1, -1, -1, -1, -1, -1, "", -1)
                static_cache[i] = cached
            result.append(cached)

    return result


def _render_chunk(result: "Text", chunk_start: int, chunk_end: int,
                   groups: tuple, styles: list[str], num_w: int,
                   seq_upper: str, show_connectors: bool,
                   sel_s: int, sel_e: int, usr_s: int, usr_e: int,
                   reh_s: int, reh_e: int, reh_color: str,
                   cursor_pos: int) -> None:
    """Render one chunk into `result`. The DNA pair takes a fast RLE path
    when no overlay (cursor / selection / RE highlight) intersects the
    chunk; otherwise the per-base path applies overlay styles. Lane rows
    above/below the DNA pair are independent of overlay so they always
    render the same way — which is what makes the static-cache reuse safe.
    """
    re_above, onebp_above, reg_above, reg_below, onebp_below, re_below = groups

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
    chunk_fwd = seq_upper[chunk_start:chunk_end]
    chunk_rev = chunk_fwd.translate(_DNA_COMP_PRESERVE_CASE)
    chunk_len = chunk_end - chunk_start

    chunk_has_overlay = (
        (usr_s < chunk_end and usr_e > chunk_start)
        or (sel_s < chunk_end and sel_e > chunk_start)
        or (reh_s < chunk_end and reh_e > chunk_start)
        or (chunk_start <= cursor_pos < chunk_end)
    )

    if not chunk_has_overlay:
        # RLE the styles slice. Long runs are common because most bases
        # outside features carry the default `color(252)`.
        runs: list[tuple[int, int, str]] = []
        if chunk_len > 0:
            cur_sty = styles[chunk_start]
            rs = 0
            for j in range(1, chunk_len):
                s_j = styles[chunk_start + j]
                if s_j != cur_sty:
                    runs.append((rs, j, cur_sty))
                    rs = j
                    cur_sty = s_j
            runs.append((rs, chunk_len, cur_sty))

        result.append(f"{chunk_start + 1:>{num_w}}  ", style="color(245)")
        for rs, re_end, sty in runs:
            result.append(chunk_fwd[rs:re_end], style=sty)
        result.append("\n")

        result.append(" " * (num_w + 2), style="color(245)")
        for rs, re_end, sty in runs:
            result.append(chunk_rev[rs:re_end], style=sty)
        result.append("\n")
    else:
        def _strand_chars(bases: "list[tuple[str, str]]") -> None:
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
                # All overlays use the same high-contrast palette
                # (white background, black foreground) so the user's
                # selection / cursor never blends into a feature tint.
                # Cursor is bold to distinguish it from a passive selection.
                fwd_sty = "black on white bold"
                rev_sty = fwd_sty
            elif in_re:
                # RE highlight stays color-tinted — the colour carries
                # enzyme identity, which is information we'd lose by
                # collapsing every overlay to white-on-black.
                fwd_sty = f"reverse bold {reh_color}"
                rev_sty = fwd_sty
            elif in_usr:
                fwd_sty = "black on white"
                rev_sty = fwd_sty
            elif in_sel:
                fwd_sty = "black on white bold underline"
                rev_sty = fwd_sty
            else:
                fwd_sty = base
                rev_sty = base
            fwd_bases.append((chunk_fwd[j], fwd_sty))
            rc_bases.append( (chunk_rev[j], rev_sty))

        result.append(f"{chunk_start + 1:>{num_w}}  ", style="color(245)")
        _strand_chars(fwd_bases)
        result.append("\n")

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
    except (OSError, UnicodeEncodeError):
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
    except (ImportError, OSError):
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


def _cursor_row_key(table) -> "str | None":
    """Return the value of a DataTable's cursor row key, or None when
    the table is empty or the cursor is out of bounds.

    Centralises the boilerplate `list(t.rows.keys())` + bounds-check
    that was open-coded at ~10 sites (library panel buttons, picker
    modals, primer table). Always pair with the empty-table branch in
    the caller — this helper is read-only.
    """
    if table.row_count == 0:
        return None
    row_keys = list(table.rows.keys())
    if not (0 <= table.cursor_row < len(row_keys)):
        return None
    rk = row_keys[table.cursor_row]
    return rk.value if rk else None


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

    _atomic_write_text(p, text)

    _log.info(
        "Exported GenBank to %s (%d bp, %d features)",
        p, len(normalized.seq), len(normalized.features),
    )
    return {"path": str(p), "bp": len(normalized.seq),
            "features": len(normalized.features)}


def _export_fasta_to_path(name: str, sequence: str, path) -> dict:
    """Write `sequence` to `path` as a single-record FASTA. Atomic write.

    Returns `{"path", "bp", "name"}` on success. Raises:
      ValueError  — empty name or empty sequence.
      OSError     — filesystem failures (write, replace, fsync).

    The sequence is written on a single line (no hard-wrap at 80 chars);
    that matches what Biopython's default SeqIO writer emits for us
    elsewhere and keeps downstream `grep`/`awk` one-liners simple.
    """
    from pathlib import Path as _Path

    header = (name or "").strip()
    seq = (sequence or "").strip().upper()
    if not header:
        raise ValueError("FASTA export needs a non-empty record name.")
    if not seq:
        raise ValueError("FASTA export needs a non-empty sequence.")

    p = _Path(path).expanduser()
    _atomic_write_text(p, f">{header}\n{seq}\n")

    _log.info("Exported FASTA to %s (%s, %d bp)", p, header, len(seq))
    return {"path": str(p), "bp": len(seq), "name": header}


# pLannotate integration removed — the panel's annotate button became a
# back-to-collections button, and the Shift+A keybinding + Edit menu entry
# are gone. Users who want auto-annotation can still run pLannotate
# externally and import the resulting GenBank file via File > Open.


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

    @on(DataTable.RowSelected, "#feat-table")
    def _row_selected(self, event: DataTable.RowSelected):
        """Click on a row whose cursor is already on it doesn't fire
        RowHighlighted (the cursor didn't move), so a user clicking the same
        sidebar row twice — or clicking a row that was already
        highlight-driven by a map/seq-panel click — would otherwise be
        silent. RowSelected fires on every click; combine with
        RowHighlighted's `_prog_row` gate to avoid double-firing on
        first-time clicks (which fire both events).
        """
        if self._populating:
            return
        # `_row_highlighted` will already have fired when the cursor moved,
        # so skip the duplicate. Only emit when the click landed on the
        # already-current cursor row (no preceding RowHighlighted).
        try:
            t = self.query_one("#feat-table", DataTable)
        except NoMatches:
            return
        if event.cursor_row != t.cursor_row:
            return
        # No `_prog_row` check here — programmatic moves don't trigger
        # RowSelected, only real clicks/Enter do, so we always want to react.
        self.post_message(self.RowActivated(event.cursor_row))


# ── Library panel ──────────────────────────────────────────────────────────────

class LibraryPanel(Widget):
    """Left-hand plasmid panel — toggles between a *collections* list and
    the active collection's *plasmids* view.

    Two modes:
      * ``"collections"`` — a list of named collections. ``+ / − / ✎``
        buttons add / remove / rename collections in-place; clicking a
        row enters that collection's plasmid view.
      * ``"plasmids"`` — the existing per-plasmid view (the SnapGene-
        style library list). ``←`` returns to the collections view;
        ``+`` saves the currently-loaded record into the active
        collection (via ``_save_library`` → ``_sync_active_collection_plasmids``).

    The two on-disk files (``plasmid_library.json`` + ``collections.json``)
    stay in sync because every mutation routes through ``_save_library``,
    which mirrors entries into the active collection.
    """

    DEFAULT_CSS = """
    LibraryPanel {
        width: 26;
        border-right: solid $primary;
    }
    #lib-hdr        { background: $primary; padding: 0 1; }
    #lib-table      { height: 1fr; }
    #lib-coll-table { height: 1fr; }
    #lib-btns       { height: 3; }
    #lib-btns Button       { min-width: 5; margin: 0 0 0 1; }
    #lib-coll-btns         { height: 3; }
    #lib-coll-btns Button  { min-width: 5; margin: 0 0 0 1; }
    """

    class PlasmidLoad(Message):
        """User selected a library entry to load."""
        def __init__(self, entry: dict):
            self.entry = entry
            super().__init__()

    class AddCurrentRequested(Message):
        """User pressed '+' to add the currently-loaded record."""
        pass

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

    class CollectionSwitched(Message):
        """User picked a different collection — library has been swapped.
        The app refreshes its caches in response (no record changes)."""
        def __init__(self, name: str):
            self.name = name
            super().__init__()

    def on_descendant_focus(self, _event):
        # DataTable inside the library panel is what actually gets focus;
        # Textual reports it to us via a DescendantFocus event. Propagate
        # to the app as a panel-level signal.
        self.post_message(self.GainedFocus())

    def compose(self) -> ComposeResult:
        yield Static(" Collections", id="lib-hdr")
        # Collections-view widgets ────────────────────────────────────
        yield DataTable(id="lib-coll-table", cursor_type="row",
                        zebra_stripes=True)
        with Horizontal(id="lib-coll-btns"):
            yield Button("+", id="btn-coll-add", variant="primary",
                         tooltip="New collection")
            yield Button("−", id="btn-coll-del", variant="error",
                         tooltip="Remove selected collection")
            yield Button("✎", id="btn-coll-rename", variant="default",
                         tooltip="Rename selected collection")
        # Plasmids-view widgets ────────────────────────────────────
        yield DataTable(id="lib-table", cursor_type="row",
                        zebra_stripes=True)
        with Horizontal(id="lib-btns"):
            yield Button("+", id="btn-lib-add", variant="primary",
                         tooltip="Save loaded plasmid to this collection")
            yield Button("−", id="btn-lib-del", variant="error",
                         tooltip="Remove selected plasmid")
            yield Button("←", id="btn-lib-back", variant="primary",
                         tooltip="Back to collections")
            yield Button("✎", id="btn-lib-rename", variant="default",
                         tooltip="Rename selected plasmid")

    def on_mount(self):
        self._active_id:    "str | None" = None
        self._active_dirty: bool         = False
        # Start in plasmids view if the user already has an active
        # collection (returning user picks up where they left off);
        # else show the collections list so they can pick.
        self._view_mode: str = (
            "plasmids" if _get_active_collection_name() else "collections"
        )
        # Two-click activation guard for collection rows — armed by the
        # first RowSelected, disarmed by a 1.5 s timer or by switching
        # rows. See `_coll_row_selected`.
        self._coll_armed_name: "str | None" = None
        self._coll_arm_timer = None
        coll = self.query_one("#lib-coll-table", DataTable)
        coll.add_columns("Name", "Plasmids")
        plas = self.query_one("#lib-table", DataTable)
        plas.add_columns("Name", "bp")
        self._apply_view_mode()
        self._repopulate()

    # ── View-mode toggle ────────────────────────────────────────────────────

    def _apply_view_mode(self) -> None:
        is_coll = (self._view_mode == "collections")
        try:
            self.query_one("#lib-coll-table").display = is_coll
            self.query_one("#lib-coll-btns").display  = is_coll
            self.query_one("#lib-table").display      = not is_coll
            self.query_one("#lib-btns").display       = not is_coll
        except NoMatches:
            return
        self._update_header()

    def _update_header(self) -> None:
        try:
            hdr = self.query_one("#lib-hdr", Static)
        except NoMatches:
            return
        # Dirty marker shows in BOTH views — even in collections view
        # the user should see at a glance that there are unsaved edits
        # somewhere in the active collection.
        prefix = "* " if self._active_dirty else " "
        if self._view_mode == "collections":
            hdr.update(f"{prefix}Collections")
            return
        active = _get_active_collection_name() or "Library"
        # Cap to keep the 26-cell-wide panel from overflowing; reserve
        # 2 cells for the dirty marker.
        hdr.update(f"{prefix}{active[:22]}")

    # ── Repopulate dispatch ────────────────────────────────────────────────

    def _repopulate(self) -> None:
        if self._view_mode == "collections":
            self._repopulate_collections()
        else:
            self._repopulate_plasmids()

    def _repopulate_collections(self) -> None:
        t = self.query_one("#lib-coll-table", DataTable)
        t.clear()
        for c in _load_collections():
            name = c.get("name") or "?"
            n_plas = len(c.get("plasmids", []) or [])
            t.add_row(name[:14], str(n_plas), key=name)

    def _repopulate_plasmids(self) -> None:
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

    # ── Plasmid view: existing flow + back button ──────────────────────────

    def add_entry(self, record) -> None:
        """Serialize record and persist into the active collection."""
        gb_text = _record_to_gb_text(record)
        entries = _load_library()
        entries = [e for e in entries if e.get("id") != record.id]
        entries.insert(0, {
            "name":    record.name or record.id,
            "id":      record.id,
            "size":    len(record.seq),
            "n_feats": len([f for f in record.features if f.type != "source"]),
            "source":  getattr(record, "_tui_source", f"id:{record.id}"),
            "added":   _date.today().isoformat(),
            "gb_text": gb_text,
        })
        _save_library(entries)
        if self._view_mode == "plasmids":
            self._repopulate_plasmids()

    @on(DataTable.RowSelected, "#lib-table")
    def _row_selected(self, event: DataTable.RowSelected):
        if self._view_mode != "plasmids":
            return
        key = event.row_key.value if event.row_key else None
        if key is None:
            return
        for entry in _load_library():
            if entry.get("id") == key:
                self.post_message(self.PlasmidLoad(entry))
                return

    @on(Button.Pressed, "#btn-lib-back")
    def _btn_back(self):
        # If the loaded plasmid has unsaved edits, prompt before
        # navigating away. Without this, the asterisks in the table /
        # header are a heads-up but it's still easy to forget; the
        # modal forces the user to make a deliberate choice.
        app = self.app
        if not getattr(app, "_unsaved", False):
            self._do_back()
            return

        def _on_response(result: "str | None") -> None:
            if result == "save":
                # _do_save marks clean on success and notifies. If the
                # save fails (e.g. write error), stay in plasmids view
                # so the user can retry.
                if hasattr(app, "_do_save") and app._do_save():
                    self._do_back()
            elif result == "discard":
                if hasattr(app, "_discard_changes"):
                    app._discard_changes()
                self._do_back()
            # None → cancel; user stays in plasmids view.

        app.push_screen(
            UnsavedNavigateModal("go back to collections"),
            callback=_on_response,
        )

    def _do_back(self) -> None:
        """Switch to collections view. Called either directly (no
        unsaved edits) or after the user resolves the unsaved-prompt."""
        self._view_mode = "collections"
        self._apply_view_mode()
        self._repopulate_collections()

    @on(Button.Pressed, "#btn-lib-add")
    def _btn_add(self):
        self.post_message(self.AddCurrentRequested())

    @on(Button.Pressed, "#btn-lib-rename")
    def _btn_rename(self):
        # Rename only works on the row with the DataTable cursor — if the
        # library is empty or no row is focused, we send None and the app
        # will notify the user.
        entry_id = _cursor_row_key(self.query_one("#lib-table", DataTable))
        self.post_message(self.RenameRequested(entry_id))

    def set_active(self, entry_id: "str | None") -> None:
        """Mark which library entry is currently loaded (clears dirty flag)."""
        self._active_id    = entry_id
        self._active_dirty = False

    def set_dirty(self, dirty: bool) -> None:
        """Show unsaved-changes marker on the active row and in the panel header.

        Called from `_mark_dirty`/`_mark_clean` on EVERY keystroke that
        flips `_unsaved`, so we early-return when the state is unchanged
        and skip the table refresh in collections-view mode (where the
        plasmid table isn't visible anyway).
        """
        if self._active_dirty == dirty:
            return
        self._active_dirty = dirty
        if self._view_mode == "plasmids":
            self._refresh_active_row()
        self._update_header()

    def _refresh_active_row(self) -> None:
        """Update just the active plasmid's Name cell, not the whole table.

        Falls back to a full repopulate if the DataTable's incremental
        API isn't available (older Textual) or the active row can't be
        located. The fallback is correctness-preserving but defeats the
        keystroke-hot-path optimisation; it should never fire under
        Textual ≥ 8.2.3 (our pinned floor).
        """
        if not self._active_id:
            return
        try:
            t = self.query_one("#lib-table", DataTable)
        except NoMatches:
            return
        active_entry = next(
            (e for e in _load_library() if e.get("id") == self._active_id),
            None,
        )
        if active_entry is None:
            return
        nm = active_entry.get("name", "?")
        display = ("*" + nm)[:14] if self._active_dirty else nm[:14]
        try:
            from textual.coordinate import Coordinate
            for i, row_key in enumerate(t.rows.keys()):
                if row_key.value == self._active_id:
                    t.update_cell_at(Coordinate(i, 0), display)
                    return
        except Exception:
            self._repopulate_plasmids()

    @on(Button.Pressed, "#btn-lib-del")
    def _btn_del(self):
        self.request_delete_under_cursor()

    def request_delete_under_cursor(self) -> None:
        """Single entry point for both the `−` button and the keyboard
        Delete key — dispatches to the right confirm flow based on
        view mode. Plasmid view → one confirm; collections view → two
        confirms (the second one is the loud red warning) so a single
        keypress can never wipe out a whole library.
        """
        if self._view_mode == "plasmids":
            self._request_plasmid_delete()
        else:
            self._request_collection_delete()

    def _request_plasmid_delete(self) -> None:
        entry_id = _cursor_row_key(self.query_one("#lib-table", DataTable))
        if entry_id is None:
            return
        entry = next(
            (e for e in _load_library() if e.get("id") == entry_id),
            None,
        )
        if entry is None:
            return
        name = entry.get("name") or entry_id
        size = entry.get("size", 0) or 0

        def _on_confirm(yes: "bool | None") -> None:
            if not yes:
                return
            entries = [e for e in _load_library() if e.get("id") != entry_id]
            _save_library(entries)
            self._repopulate_plasmids()
            # If we just deleted the loaded record's library entry, drop
            # the panel's active-row binding so the dirty asterisk doesn't
            # point at a row that no longer exists.
            app = self.app
            cur = getattr(app, "_current_record", None)
            if cur is not None and cur.id == entry_id:
                self.set_active(None)

        self.app.push_screen(
            LibraryDeleteConfirmModal(name, size, entry_id),
            callback=_on_confirm,
        )

    # ── Collections view: enter / + / − / ✎ ────────────────────────────────

    @on(DataTable.RowSelected, "#lib-coll-table")
    def _coll_row_selected(self, event: DataTable.RowSelected):
        """Loading a collection swaps the entire library + active pointer,
        so we require an explicit double-activation: the first
        RowSelected on a row arms a confirmation hint, and only a second
        activation on the same row inside the timeout actually loads.
        Mirrors the "click twice to commit" pattern users get on the
        plasmid table (where the first click moves the cursor and the
        second activates).
        """
        if self._view_mode != "collections":
            return
        rk = event.row_key
        name = rk.value if rk else None
        if not name:
            return

        # Cancel any prior arm-timer; only the most recent click matters.
        prior_timer = getattr(self, "_coll_arm_timer", None)
        if prior_timer is not None:
            try:
                prior_timer.stop()
            except Exception:
                pass
            self._coll_arm_timer = None

        if getattr(self, "_coll_armed_name", None) != name:
            # First activation on this row — arm a 1.5 s window during
            # which a second activation will actually load.
            self._coll_armed_name = name
            self.app.notify(
                f"Click '{name}' again to load it.",
                timeout=2,
            )

            def _disarm() -> None:
                self._coll_armed_name = None
                self._coll_arm_timer = None

            try:
                self._coll_arm_timer = self.set_timer(1.5, _disarm)
            except Exception:
                # Timer setup can fail during teardown; fall back to a
                # single-shot arm with no auto-disarm.
                self._coll_arm_timer = None
            return

        # Second activation within the window — load it.
        self._coll_armed_name = None
        coll = _find_collection(name)
        if coll is None:
            return
        # Set active BEFORE writing the library so _save_library's mirror
        # writes back to the correct collection.
        _set_active_collection_name(name)
        plasmids = [dict(p) for p in (coll.get("plasmids") or [])
                    if isinstance(p, dict)]
        _save_library(plasmids)
        self._view_mode = "plasmids"
        self._apply_view_mode()
        self._repopulate_plasmids()
        self.post_message(self.CollectionSwitched(name))

    @on(Button.Pressed, "#btn-coll-add")
    def _btn_coll_add(self):
        def _picked(name: "str | None") -> None:
            if not name:
                return
            if _collection_name_taken(name):
                self.app.notify(
                    f"Collection '{name}' already exists.",
                    severity="warning",
                )
                return
            existing = _load_collections()
            existing.append({
                "name":        name,
                "description": "",
                "plasmids":    [],
                "saved":       _date.today().isoformat(),
            })
            _save_collections(existing)
            self._repopulate_collections()
        self.app.push_screen(
            CollectionNameModal("New collection", "", "Collection name"),
            callback=_picked,
        )

    @on(Button.Pressed, "#btn-coll-del")
    def _btn_coll_del(self):
        self._request_collection_delete()

    def _request_collection_delete(self) -> None:
        """Two-stage confirm. First modal asks "delete this collection?"
        with default-No focus; only on Yes does the second loud-red
        modal fire ("ARE YOU ABSOLUTELY SURE?", also default-No).
        Either No (or Esc) keeps the collection. Bound to both the
        `−` button and the keyboard Delete key when the panel is in
        collections view.
        """
        name = _cursor_row_key(self.query_one("#lib-coll-table", DataTable))
        if not name:
            return
        coll = _find_collection(name)
        n_plas = len((coll or {}).get("plasmids", []) or [])

        def _on_first(yes: "bool | None") -> None:
            if not yes:
                return

            def _on_second(yes2: "bool | None") -> None:
                if not yes2:
                    return
                remaining = [c for c in _load_collections()
                             if c.get("name") != name]
                _save_collections(remaining)
                if _get_active_collection_name() == name:
                    _set_active_collection_name(None)
                self._repopulate_collections()
                self.app.notify(f"Deleted collection '{name}'.")

            self.app.push_screen(
                ScaryDeleteConfirmModal(name, n_plas),
                callback=_on_second,
            )

        self.app.push_screen(
            CollectionDeleteConfirmModal(name, n_plas),
            callback=_on_first,
        )

    @on(Button.Pressed, "#btn-coll-rename")
    def _btn_coll_rename(self):
        old = _cursor_row_key(self.query_one("#lib-coll-table", DataTable))
        if not old:
            return

        def _picked(new_name: "str | None") -> None:
            if not new_name or new_name == old:
                return
            if _collection_name_taken(new_name):
                self.app.notify(
                    f"Collection '{new_name}' already exists.",
                    severity="warning",
                )
                return
            existing = _load_collections()
            for c in existing:
                if c.get("name") == old:
                    c["name"] = new_name
                    break
            _save_collections(existing)
            if _get_active_collection_name() == old:
                _set_active_collection_name(new_name)
            self._repopulate_collections()
            self._update_header()

        self.app.push_screen(
            CollectionNameModal("Rename collection", old, "Collection name"),
            callback=_picked,
        )


# ── Sequence panel ─────────────────────────────────────────────────────────────

class SequencePanel(Widget):
    """
    Bottom DNA sequence viewer.

    Click  on the sequence → place cursor + select feature at that position.
    Shift+click           → extend selection for editing.
    Ctrl+E                → open insert/replace dialog at cursor / selection.
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
        """User clicked on a base or a feature lane; app routes the
        outcome based on `from_lane`. Lane clicks → highlight the
        whole feature. DNA-row clicks → just place the cursor on `bp`,
        no feature-wide highlight even if `bp` happens to fall inside
        a feature's range.
        """
        def __init__(self, bp: int, double: bool = False,
                     from_lane: bool = False):
            self.bp        = bp
            self.double    = double
            self.from_lane = from_lane
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
            " Sequence  (click: select · Shift+click: select region · Ctrl+E: edit)",
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

    def highlight_feature(self, feat: "dict | None", cursor_bp: int = -1,
                            scroll: bool = True) -> None:
        """Highlight a feature's region in the sequence; show CDS translation.

        cursor_bp: if >= 0, anchor the cursor (and scroll) at this bp position.
                   Use this for sequence-panel clicks so scroll stays at the
                   clicked position rather than jumping to the feature start.
        scroll:    when True (default), call `_ensure_cursor_visible` so the
                   cursor is on screen. Pass False when the caller will
                   issue its own `center_on_bp` afterwards — otherwise two
                   sequential scrolls become a visible jitter.
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
            if scroll:
                self._ensure_cursor_visible()   # scroll BEFORE refresh
        self._refresh_view()

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

    def select_feature_range(self, feat: dict, cursor_bp: int = -1,
                              scroll: bool = True) -> None:
        """Highlight the entire feature span as a copyable selection.

        cursor_bp: if >= 0 (sequence-panel click), keep cursor at that bp and
                   anchor the Shift+arrow selection there. Otherwise cursor goes
                   to the feature end so Shift+arrow naturally extends outward.
        scroll:    when True, also adjust the viewport so the cursor is on
                   screen. Pass False when the caller will issue its own
                   `center_on_bp` immediately after — otherwise the user
                   sees a partial scroll (just-visible) followed by a
                   centring scroll, which reads as a perceptible jitter on
                   sidebar arrow-key navigation.
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
        # Scroll BEFORE refresh — see _ensure_cursor_visible docstring.
        if scroll:
            self._ensure_cursor_visible()
        self._refresh_view()

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
        # Scroll BEFORE refresh — see _ensure_cursor_visible docstring.
        self._ensure_cursor_visible()
        self._refresh_view()

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
        # `_last_lane_click` is set as a side-effect of `_click_to_bp` when
        # the click lands on a feature lane (the bar/arrow art) rather
        # than the DNA strand. Sequence-row clicks intentionally do NOT
        # trigger a whole-feature highlight — only the lane art does.
        self._last_lane_click = False
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
        self.post_message(self.SequenceClick(
            bp, double=double, from_lane=self._last_lane_click,
        ))

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
            """Check if click hit a feature in this lane; return bp or None.
            Sets `_last_lane_click=True` so on_click can route a lane hit
            to a whole-feature highlight while a DNA-row hit just places
            the cursor."""
            for f in lane:
                bar_s = max(f["start"], chunk_start) - chunk_start
                bar_e = min(f["end"],   chunk_end)   - chunk_start
                if bar_s <= seq_col < bar_e:
                    if f.get("type") == "resite":
                        self._last_resite_click = f
                    else:
                        self._last_lane_click = True
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
        """Return the content row index (0-based) of the DNA line containing bp.

        O(1) via `_chunk_layout` prefix sums — direct index by `bp //
        line_width` plus arithmetic on cached above/below pair counts. The
        previous chunk-by-chunk re-scan was the bottleneck for cursor
        scrolling on cosmid/BAC-scale records (~50 ms/keystroke at 50 kb).
        """
        n = len(self._seq)
        if n == 0:
            return 0
        line_width = self._line_width()
        if line_width <= 0:
            return 0
        chunks_layout, prefix_dna2, prefix_lanes = _chunk_layout(
            self._seq, self._feats, line_width
        )
        if not chunks_layout:
            return 0
        rpg = 2 + (1 if self._show_connectors else 0)
        chunk_idx = min(bp // line_width, len(chunks_layout) - 1)
        # Total rows up to (not including) chunk_idx, scaled for current rpg.
        rows_before = prefix_dna2[chunk_idx] + (rpg - 2) * prefix_lanes[chunk_idx]
        above_rows  = chunks_layout[chunk_idx][3] * rpg
        return rows_before + above_rows

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

    def center_on_bp(self, bp: int) -> None:
        """Scroll the sequence panel so `bp` lands at the vertical centre of
        the viewport. Used by sidebar/map click handlers — the user has
        deliberately clicked, so don't make them hunt for the result.

        Defers via `call_after_refresh` so the scroll runs after the queued
        `view.update()` tick from the preceding highlight call. Uses
        `scroll_to(force=True)` rather than `scroll_y = ...`: a direct
        attribute set on a Reactive can be silently reverted by Textual's
        scroll-target watcher when it fires next, but `scroll_to` goes
        through the proper code path and `force=True` skips clamping.
        Symptom of the broken assignment: scroll moves to centre, then
        reverts on the next event — visible as a "snap back" jitter.

        If the viewport hasn't been laid out yet on the first refresh
        (`vp_h == 0`), retry once on the next refresh. Without the retry
        the centring silently no-ops under suite load — release.py's
        serial test path was hitting this ~30 % of the time.
        """
        if not self._seq or bp < 0:
            return
        row = self._bp_to_content_row(bp)

        def _do_scroll(remaining_retries: int = 1) -> None:
            try:
                scroll = self.query_one("#seq-scroll", ScrollableContainer)
            except NoMatches:
                return
            vp_h = scroll.size.height
            if vp_h <= 0:
                if remaining_retries > 0:
                    self.call_after_refresh(
                        lambda: _do_scroll(remaining_retries - 1)
                    )
                return
            target_top = max(0, row - vp_h // 2)
            scroll.scroll_to(0, target_top, animate=False, force=True)

        self.call_after_refresh(_do_scroll)

    def _ensure_cursor_visible(self) -> None:
        """Follow the cursor: scroll the minimum necessary, no recentering.

        "Chunk" = the cursor's DNA pair plus the feature lanes (with
        labels) above and below it. Behaviour:

          - Chunk already fully in view → do nothing. The scroll bar
            stays put for every arrow press whose cursor lands in a
            row already on screen — no "pull back to recenter".
          - Chunk extends above the viewport → scroll up by exactly
            `vp_top - chunk_top` so `chunk_top` (the topmost label row
            of the above-lanes) lands at the viewport top.
          - Chunk extends below → scroll down by exactly
            `chunk_bottom - vp_bottom` so `chunk_bottom` (the bottommost
            below-lane row) lands at the viewport bottom.
          - Chunk taller than the viewport (very dense lane stacking) →
            track the cursor's DNA row only; whichever side the chunk
            spills off-screen on gets clipped. Direction-aware: if the
            DNA is above the viewport we pin it to vp_top (Up arrow
            into a tall chunk), if below we pin it to vp_bottom.

        No snap-to-extreme on first/last chunks. The user arrow nav
        wants the scroll bar to track row-by-row, not jump to either
        end of the document — and `chunk_top == 0` for chunk 0 (so
        target_top naturally lands at 0 when scrolling up into it),
        `chunk_bottom - vp_h + 1 == max_scroll_y` for the final chunk.
        The earlier explicit-snap branches always fired regardless of
        whether scrolling was needed, which read as a sudden jump.

        Uses `scroll_to(animate=False, immediate=True)` rather than the
        reactive-bypass `set_scroll`. `set_scroll` sets `scroll_y` via
        `set_reactive`, which skips the `watch_scroll_y` watcher — and
        that watcher is what propagates the new position to
        `vertical_scrollbar.position`. The scrollbar widget therefore
        renders one tick stale: content scrolls, the bar doesn't, the
        next event jolts the bar to where it should have been. That
        was the visible "borked tracking". `scroll_to` runs the full
        reactive pipeline (validators + watchers) so the bar repaints
        in lock-step with the content. `immediate=True` skips the
        default `call_after_refresh` defer; `animate=False` skips the
        smooth-scroll path entirely.

        Must be called BEFORE `_refresh_view()` so the new content
        paints with the new scroll position in the same tick.
        """
        if self._cursor_pos < 0 or not self._seq:
            return
        line_width = self._line_width()
        if line_width <= 0:
            return
        chunks_layout, prefix_dna2, prefix_lanes = _chunk_layout(
            self._seq, self._feats, line_width
        )
        if not chunks_layout:
            return

        rpg = 2 + (1 if self._show_connectors else 0)
        chunk_idx = min(self._cursor_pos // line_width, len(chunks_layout) - 1)
        above_pairs = chunks_layout[chunk_idx][3]
        below_pairs = chunks_layout[chunk_idx][4]

        chunk_top    = prefix_dna2[chunk_idx] + (rpg - 2) * prefix_lanes[chunk_idx]
        chunk_bottom = chunk_top + (above_pairs + below_pairs) * rpg + 1

        try:
            scroll = self.query_one("#seq-scroll", ScrollableContainer)
        except NoMatches:
            return
        vp_top = int(scroll.scroll_y)
        vp_h   = scroll.size.height
        if vp_h <= 0:
            return
        vp_bottom = vp_top + vp_h - 1
        max_y    = scroll.max_scroll_y

        chunk_height = chunk_bottom - chunk_top + 1
        dna_row      = chunk_top + above_pairs * rpg

        if chunk_height > vp_h:
            # Chunk doesn't fit — track the DNA row, lanes get clipped.
            # Direction-aware so an Up arrow into a tall chunk pins the
            # DNA at vp_top (above-lanes fill the viewport, below
            # clipped) and a Down arrow pins it at vp_bottom (below
            # clipped, above-lanes still as visible as they fit).
            if dna_row < vp_top:
                target_top = dna_row
            elif dna_row + 1 > vp_bottom:
                target_top = dna_row - vp_h + 2  # DNA pair at viewport bottom
            else:
                return
        elif chunk_top >= vp_top and chunk_bottom <= vp_bottom:
            # Already fully in view — do not scroll.
            return
        elif chunk_top < vp_top:
            # Scroll up just enough to expose the above-lanes + DNA.
            target_top = chunk_top
        else:
            # chunk_bottom > vp_bottom: scroll down just enough to fit
            # below-lanes.
            target_top = chunk_bottom - vp_h + 1

        target_top = max(0, min(target_top, max_y))
        if int(vp_top) == int(target_top):
            return
        scroll.scroll_to(0, target_top, animate=False, immediate=True)

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

        # Don't try to "preserve scroll across content update" here. An
        # earlier version captured `scroll.scroll_y` before `view.update`
        # and re-applied it via `call_after_refresh`. That fought with
        # `_ensure_cursor_visible`, which sets a NEW scroll target right
        # before `_refresh_view`: Textual hadn't propagated the new
        # `scroll_y` yet, so we'd capture the OLD value, then the
        # deferred would restore the OLD position — visible as a quick
        # scroll up followed by a snap back. Callers that change content
        # should set scroll before calling `_refresh_view`; the sync set
        # survives the refresh tick on its own.
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


class FastaExportModal(ModalScreen):
    """Prompt for a target path and write `(name, sequence)` as FASTA.

    Dismisses with the summary dict from `_export_fasta_to_path`
    (`path`, `bp`, `name`) or `None` if cancelled. Caller handles
    the success notification.
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "focus_next", "Next",   show=False),
    ]

    def __init__(self, name: str, sequence: str,
                 default_path: str = "", subtitle: str = "") -> None:
        super().__init__()
        self._name = name
        self._sequence = sequence
        self._default_path = default_path
        self._subtitle = subtitle

    def compose(self) -> ComposeResult:
        with Vertical(id="fasta-export-box"):
            yield Static(" Export as FASTA ", id="fasta-export-title")
            bp = len(self._sequence or "")
            sub = self._subtitle or f"[{self._name}]  {bp} bp"
            yield Label(sub)
            yield Label("Output path (.fa / .fasta):")
            yield Input(
                value=self._default_path,
                placeholder="/path/to/sequence.fa",
                id="fasta-export-path",
            )
            with Horizontal(id="fasta-export-btns"):
                yield Button("Export", id="btn-fasta-export-ok", variant="primary")
                yield Button("Cancel", id="btn-fasta-export-cancel")
            yield Static("", id="fasta-export-status", markup=True)

    def on_mount(self) -> None:
        try:
            self.query_one("#fasta-export-path", Input).focus()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-fasta-export-ok")
    def _do_export(self) -> None:
        try:
            inp = self.query_one("#fasta-export-path", Input)
            status = self.query_one("#fasta-export-status", Static)
        except NoMatches:
            return
        path = inp.value.strip()
        if not path:
            status.update("[red]Enter an output path.[/red]")
            return
        try:
            summary = _export_fasta_to_path(self._name, self._sequence, path)
        except (OSError, ValueError) as exc:
            _log.exception("FASTA export to %s failed", path)
            status.update(f"[red]Export failed: {exc}[/red]")
            return
        self.dismiss(summary)

    @on(Input.Submitted)
    def _submitted(self) -> None:
        try:
            self.query_one("#btn-fasta-export-ok", Button).press()
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-fasta-export-cancel")
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
                # "Features" is a direct-open workbench (no dropdown).
                # Every other menu surfaces items via DropdownScreen.
                if name == "Features":
                    self.app.push_screen(FeatureLibraryScreen())
                    break
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

# Coding-DNA part types: these are the only types where silent (synonymous)
# codon substitution can be used to remove an internal Type IIS site during
# domestication. Non-coding parts (promoters, UTRs, terminators) have no
# reading frame, so internal sites must be fixed manually or by picking a
# different template region.
_GB_CODING_PART_TYPES: frozenset[str] = frozenset({"CDS", "CDS-NS", "C-tag"})

# Parts-bin TitleCase type → INSDC feature_type. Used by
# `PartsBinModal._save_as_feature` to round-trip a Golden Braid part into
# the feature library, where types follow the INSDC vocabulary.
# CDS-NS / C-tag have no INSDC equivalent — they're GB-specific
# coding-DNA shapes — so they collapse to plain CDS; the GB position
# survives in the feature's description string instead.
_GB_PART_TYPE_TO_INSDC: dict[str, str] = {
    "Promoter":   "promoter",
    "5' UTR":     "5'UTR",
    "CDS":        "CDS",
    "CDS-NS":     "CDS",
    "C-tag":      "CDS",
    "Terminator": "terminator",
}

# Type IIS recognition + tail used for all Golden Braid L0 domestication
# primers. Golden Braid splits assembly across two enzymes: **L0 parts use
# Esp3I / BsmBI (CGTCTC(1/5))**, while downstream L1+ transcriptional units
# use BsaI (GGTCTC(1/5)). The two have identical N(1)/N(5) geometry (→ 4-nt
# 5' overhangs) but different recognition sequences, so a domesticated L0
# part survives a second round of Golden Gate in L1 without re-cutting.
# Padding bases improve Type IIS digestion efficiency near DNA ends.
_GB_L0_ENZYME_NAME = "Esp3I"       # Esp3I is the isoschizomer of BsmBI
_GB_L0_ENZYME_SITE = "CGTCTC"      # recognition; rc = "GAGACG"
_GB_SPACER         = "A"           # 1 nt between recognition and the overhang
_GB_PAD            = "GCGC"        # 4 nt of extra bases for efficient end-cutting


# ── Cloning simulation (Parts Bin preview) ─────────────────────────────────────
#
# The three "Copy …" buttons in the Parts Bin let the user grab the insert,
# the full PCR amplicon (insert + primer tails), or the simulated cloned
# plasmid (insert slotted into a pUPD2-shaped backbone) as plain text.
#
# The backbone here is a **placeholder** — not a real pUPD2 sequence. It is
# a deterministic pseudo-random ACGT string, scrubbed of every Type IIS site
# SpliceCraft touches (BsaI GGTCTC/GAGACC, Esp3I/BsmBI CGTCTC/GAGACG) so the
# simulated plasmid never contains a stray recognition site that would
# re-cut during a real assembly. Once we gain a licensed / verified pUPD2
# sequence we can drop it in here without changing any callers.

def _build_pupd2_backbone_stub(seed: int = 0xBACDBAC0, length: int = 420) -> str:
    """Return a deterministic ACGT string, free of BsaI/Esp3I/BsmBI sites on
    both strands, for use as a pUPD2-shaped placeholder backbone.

    Deterministic because the same insert must produce the same cloned
    sequence across sessions (otherwise the "Copy Cloned Sequence" output
    would silently drift). Seeded with a fixed constant.
    """
    import random as _random_mod
    rng = _random_mod.Random(seed)
    bases = [rng.choice("ACGT") for _ in range(length)]
    # Scrub both strands — the linear backbone becomes part of a circular
    # product, so a top-strand CGTCTC (Esp3I/BsmBI) and a bottom-strand
    # GAGACG are biologically equivalent and both must be absent.
    forbidden = ("GGTCTC", "GAGACC", "CGTCTC", "GAGACG")
    i = 0
    while i <= length - 6:
        window = "".join(bases[i:i + 6])
        if window in forbidden:
            # Flip the middle base to something that can't re-hit any
            # forbidden site; ACGT minus the current base leaves 3 choices.
            middle = i + 3
            current = bases[middle]
            for replacement in "ACGT":
                if replacement != current:
                    bases[middle] = replacement
                    break
            # Rewind a bit to catch any new site created at the boundary.
            i = max(0, i - 5)
            continue
        i += 1
    return "".join(bases)


_PUPD2_BACKBONE_STUB: str = _build_pupd2_backbone_stub()


def _simulate_primed_amplicon(
    insert: str, oh5: str, oh3: str,
    grammar: "dict | None" = None,
) -> str:
    """PCR amplicon top strand (5'→3'), as it would run on a pre-digest gel.

    Structure:  [pad] [enzyme site] [spacer] [oh5] [insert] [oh3]
                [rc(spacer)] [rc(enzyme site)] [rc(pad)]

    Matches the primer geometry in :func:`_design_gb_primers`. Defaults to
    Golden Braid L0 (Esp3I); pass ``grammar`` to use a different cloning
    grammar's enzyme/pad/spacer (e.g., MoClo Plant uses BsaI). Used by
    both DomesticatorModal (active grammar at design time) and
    PartsBinModal "Copy Primed Sequence" (the part's stored grammar).
    """
    g = grammar if isinstance(grammar, dict) else _BUILTIN_GRAMMARS["gb_l0"]
    pad    = g.get("pad",    _GB_PAD)
    site   = g.get("site",   _GB_L0_ENZYME_SITE)
    spacer = g.get("spacer", _GB_SPACER)
    left_tail  = pad + site + spacer
    right_tail = _rc(spacer) + _rc(site) + _rc(pad)
    return left_tail + oh5 + insert + oh3 + right_tail


def _simulate_cloned_plasmid(insert: str, oh5: str, oh3: str) -> str:
    """Simulated cloned circular plasmid, linearised at the 5' overhang.

    After the cloning grammar's enzyme cuts both the amplicon and the
    backbone, the insert fragment carries `oh5…oh3` on its 4-nt sticky
    ends and ligates into the backbone in a single orientation. The
    circular product, read starting at `oh5`, is:

        [oh5] [insert] [oh3] [backbone_body]

    The backbone here is `_PUPD2_BACKBONE_STUB` — a scrubbed placeholder
    that contains no BsaI/Esp3I sites on either strand, so the simulated
    plasmid is guaranteed not to re-cut in either L0 or L1 assembly.
    """
    return oh5 + insert + oh3 + _PUPD2_BACKBONE_STUB


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


_GB_DOMESTICATION_FORBIDDEN: dict[str, str] = {
    # Esp3I self-cuts during L0 domestication; BsaI would re-cut during any
    # downstream L1 assembly — both must be absent from the final part.
    "BsaI":  "GGTCTC",
    "Esp3I": "CGTCTC",
}


# ── Modular cloning grammars (Golden Braid, MoClo, custom) ─────────────────────
#
# A "grammar" parameterises every Type IIS-aware tool in the app: the
# Parts Bin (which catalog + which user parts are visible), the
# Domesticator (which positions/overhangs/enzyme/forbidden-sites apply
# to primer design), and downstream constructors. Each grammar carries:
#
#   - **enzyme + site + spacer + pad** — the Type IIS recognition + the
#     primer-tail bookkeeping bytes (`pad + site + spacer + oh5 + …`).
#   - **forbidden_sites** — every Type IIS recognition that must be
#     scrubbed from a domesticated part. For Golden Braid L0 that's
#     Esp3I (current step) + BsaI (next step); for MoClo Plant it's
#     BsaI (current) + BpiI (next). The codon-fix repair pipeline
#     iterates over this dict.
#   - **positions** — ordered list of `{name, type, oh5, oh3, color}`
#     dicts. The `type` is the user-facing part-shape label (e.g.
#     "Promoter", "CDS-NS"), shared across catalog rows so a single
#     position lookup powers the Domesticator's position dropdown,
#     the parts table coloring, and the overhang-table display.
#   - **coding_types** — subset of position types that have a reading
#     frame (so the codon-fix repair pipeline can swap synonymous
#     codons to remove an internal forbidden site).
#   - **type_to_insdc** — Title-cased part type → INSDC feature_type
#     for the Save-As-Feature flow. Lossy mappings (CDS-NS → CDS) are
#     covered in the description string.
#   - **catalog** — built-in reference parts as
#     ``(name, type, position_label, oh5, oh3, backbone, marker)``
#     tuples. Empty for grammars that ship without a reference set.
#   - **editable** — built-ins are read-only; user-defined grammars
#     are loaded from `cloning_grammars.json` and can be edited in
#     `GrammarEditorModal`.
#
# The ID space is flat: built-in IDs start with their grammar family
# (``gb_l0``, ``moclo_plant``); custom-grammar IDs are slugged from
# user-supplied names (``custom_my_assembly``, etc.) at creation.

_BUILTIN_GRAMMARS: dict[str, dict] = {
    "gb_l0": {
        "id":              "gb_l0",
        "name":            "Golden Braid L0",
        "enzyme":          _GB_L0_ENZYME_NAME,
        "site":            _GB_L0_ENZYME_SITE,
        "spacer":          _GB_SPACER,
        "pad":             _GB_PAD,
        "forbidden_sites": dict(_GB_DOMESTICATION_FORBIDDEN),
        "positions": [
            {"name": pos, "type": ptype, "oh5": oh5, "oh3": oh3,
             "color": _GB_TYPE_COLORS.get(ptype, "white")}
            for ptype, (pos, oh5, oh3) in _GB_POSITIONS.items()
        ],
        "coding_types":    sorted(_GB_CODING_PART_TYPES),
        "type_to_insdc":   dict(_GB_PART_TYPE_TO_INSDC),
        "catalog":         list(_GB_L0_PARTS),
        "editable":        False,
    },
    # Plant MoClo (Weber et al. 2011, Engler et al. 2014). BsaI at L0,
    # BpiI/BbsI at L1 — both scrubbed during domestication. Ships
    # without a built-in catalog because Plant MoClo's reference parts
    # depend heavily on the user's host system; users seed via "New
    # Part" or by duplicating into a custom grammar.
    "moclo_plant": {
        "id":              "moclo_plant",
        "name":            "MoClo Plant (Weber 2011)",
        "enzyme":          "BsaI",
        "site":            "GGTCTC",
        "spacer":          "A",
        "pad":             "GCGC",
        # BsaI for the current L0 cut; BpiI (= BbsI) for the next-level
        # MoClo assembly, which uses a different Type IIS site so the
        # L0 part survives the L1 reaction without re-cutting.
        "forbidden_sites": {"BsaI": "GGTCTC", "BpiI": "GAAGAC"},
        "positions": [
            {"name": "Pos 1", "type": "Promoter",   "oh5": "GGAG", "oh3": "AATG", "color": "green"},
            {"name": "Pos 2", "type": "5' UTR",     "oh5": "AATG", "oh3": "AGGT", "color": "cyan"},
            {"name": "Pos 3", "type": "CDS",        "oh5": "AGGT", "oh3": "GCTT", "color": "yellow"},
            {"name": "Pos 4", "type": "C-tag",      "oh5": "GCTT", "oh3": "GGTA", "color": "magenta"},
            {"name": "Pos 5", "type": "Terminator", "oh5": "GGTA", "oh3": "CGCT", "color": "blue"},
        ],
        "coding_types":    ["CDS", "C-tag"],
        "type_to_insdc": {
            "Promoter":   "promoter",
            "5' UTR":     "5'UTR",
            "CDS":        "CDS",
            "C-tag":      "CDS",
            "Terminator": "terminator",
        },
        "catalog":         [],
        "editable":        False,
    },
}


# Custom grammars persist to `cloning_grammars.json` (envelope schema,
# sacred invariant #7). Schema per entry mirrors the built-in dict
# above; ``editable`` is implicitly True for everything in this file.
_GRAMMARS_FILE = _DATA_DIR / "cloning_grammars.json"
_grammars_cache: "list | None" = None


def _load_custom_grammars() -> list[dict]:
    global _grammars_cache
    from copy import deepcopy
    if _grammars_cache is not None:
        return deepcopy(_grammars_cache)
    entries, warning = _safe_load_json(_GRAMMARS_FILE, "Cloning grammars")
    if warning:
        _log.warning(warning)
    entries = [e for e in entries if isinstance(e, dict) and isinstance(e.get("id"), str)]
    _grammars_cache = entries
    return deepcopy(_grammars_cache)


def _save_custom_grammars(entries: list[dict]) -> None:
    global _grammars_cache
    from copy import deepcopy
    _safe_save_json(_GRAMMARS_FILE, entries, "Cloning grammars")
    _grammars_cache = deepcopy(entries)


def _all_grammars() -> dict[str, dict]:
    """Return all grammars (built-in + user-defined) keyed by id.

    Built-ins come first; user-defined grammars override builtin IDs
    if they ever collide (defensive — UI prevents this on save). The
    returned dicts are independent copies, so callers may mutate them
    without poisoning the cache.
    """
    from copy import deepcopy
    out: dict[str, dict] = {gid: deepcopy(g) for gid, g in _BUILTIN_GRAMMARS.items()}
    for g in _load_custom_grammars():
        gid = g.get("id")
        if isinstance(gid, str):
            # Custom grammars are always editable regardless of what
            # the JSON file says — stops a mis-flagged file from
            # locking the user out of their own definitions.
            g = dict(g)
            g["editable"] = True
            out[gid] = g
    return out


def _get_active_grammar() -> dict:
    """Return the currently-active grammar dict. Falls back to GB L0
    if the persisted ``active_grammar`` id no longer resolves (e.g.,
    a custom grammar was deleted while still selected)."""
    grammars = _all_grammars()
    gid = _get_setting("active_grammar", "gb_l0")
    if gid in grammars:
        return grammars[gid]
    # Recover gracefully — flip the setting back to gb_l0 so we don't
    # keep falling back forever.
    _set_setting("active_grammar", "gb_l0")
    return grammars["gb_l0"]


def _grammar_position_by_type(grammar: dict, ptype: str) -> "dict | None":
    """Helper: find the position spec for a given part type within a
    grammar. ``None`` if the grammar doesn't define that type — which
    e.g. means CDS-NS isn't a valid pick under MoClo Plant."""
    for pos in grammar.get("positions", []):
        if pos.get("type") == ptype:
            return pos
    return None


def _grammar_dropdown_options() -> list[tuple[str, str]]:
    """Return ``[(display_name, id), …]`` for every grammar, in the
    canonical order used wherever a Select dropdown lists grammars
    (DomesticatorModal "Grammar" picker today; future menus likely):

      1. **Golden Braid L0 first** — the default reference grammar.
         Pinned at position 1 regardless of any other ordering
         shenanigans (e.g., a custom grammar id-sorted before
         ``gb_l0``).
      2. Other built-in grammars (MoClo Plant, etc.) in
         ``_BUILTIN_GRAMMARS`` insertion order.
      3. Custom grammars from ``cloning_grammars.json`` last, tagged
         ``(custom)`` for visual disambiguation.
    """
    grammars = _all_grammars()
    out: list[tuple[str, str]] = []
    if "gb_l0" in grammars:
        g = grammars["gb_l0"]
        out.append((f"{g.get('name', 'Golden Braid L0')}", "gb_l0"))
    for gid in _BUILTIN_GRAMMARS:
        if gid == "gb_l0" or gid not in grammars:
            continue
        g = grammars[gid]
        out.append((f"{g.get('name', gid)}", gid))
    for gid, g in grammars.items():
        if gid in _BUILTIN_GRAMMARS:
            continue
        out.append((f"{g.get('name', gid)}  (custom)", gid))
    return out


# ── App-wide settings (active grammar, future preferences) ─────────────────────

_SETTINGS_FILE = _DATA_DIR / "settings.json"
_settings_cache: "dict | None" = None


def _load_settings() -> dict:
    """Return the persistent settings dict. Stored on disk as a list of
    ``{"key": ..., "value": ...}`` envelope entries so it shares the
    schema layout (sacred invariant #7) with every other JSON file."""
    global _settings_cache
    if _settings_cache is not None:
        return dict(_settings_cache)
    entries, warning = _safe_load_json(_SETTINGS_FILE, "Settings")
    if warning:
        _log.warning(warning)
    settings: dict = {}
    for e in entries:
        if not isinstance(e, dict):
            continue
        k, v = e.get("key"), e.get("value")
        if isinstance(k, str):
            settings[k] = v
    _settings_cache = settings
    return dict(_settings_cache)


def _save_settings(settings: dict) -> None:
    global _settings_cache
    entries = [{"key": k, "value": v} for k, v in settings.items()]
    _safe_save_json(_SETTINGS_FILE, entries, "Settings")
    _settings_cache = dict(settings)


def _get_setting(key: str, default=None):
    return _load_settings().get(key, default)


def _set_setting(key: str, value) -> None:
    settings = _load_settings()
    settings[key] = value
    _save_settings(settings)


def _gb_find_forbidden_hits(
    seq: str,
    sites: "dict[str, str] | None" = None,
) -> list[tuple[str, str, int]]:
    """Return ``(enzyme_name, site_found, position)`` for **every** internal
    Type IIS recognition in *seq*, on both forward and reverse strands.

    The ``sites`` map (``{enzyme_name: recognition}``) defaults to the
    Golden Braid L0 forbidden set (Esp3I + BsaI). Pass a different
    grammar's ``forbidden_sites`` to scan against MoClo (BsaI + BpiI),
    a custom grammar, or any other Type IIS combination. Returns every
    occurrence — not just the first per enzyme. Critical for accurate
    reporting when an insert contains multiple sites: the user must
    know about all of them before paying for a gBlock synthesis.
    Results are sorted by position to aid downstream reporting.
    """
    if sites is None:
        sites = _GB_DOMESTICATION_FORBIDDEN
    out: list[tuple[str, str, int]] = []
    for name, site in sites.items():
        if not isinstance(site, str) or not site:
            continue
        rc = _rc(site)
        needles = [site] if rc == site else [site, rc]
        for needle in needles:
            start = 0
            while True:
                pos = seq.find(needle, start)
                if pos == -1:
                    break
                out.append((name, needle, pos))
                start = pos + 1
    out.sort(key=lambda t: (t[2], t[0], t[1]))
    return out


_CODON_FIX_POS_RE = re.compile(r"codon (\d+)")


def _codon_fix_mutation_positions(mutations: list[str]) -> list[int]:
    """Given the string list returned by :func:`_codon_fix_sites`, return
    each mutation's 0-based codon-start nucleotide position in the insert.

    The mutation format is fixed by ``_codon_fix_sites`` — ``(codon N …)``
    where ``N`` is 1-based. A missing / malformed entry gets ``-1`` so
    callers can filter without raising.
    """
    out: list[int] = []
    for m in mutations:
        match = _CODON_FIX_POS_RE.search(m) if isinstance(m, str) else None
        if match:
            out.append((int(match.group(1)) - 1) * 3)
        else:
            out.append(-1)
    return out


def _gb_binding_region_advisory(
    mutations: list[str],
    insert_len: int,
    fwd_bind_len: int,
    rev_bind_len: int,
) -> list[dict]:
    """Return one entry per mutation that lands inside a primer binding
    window. Each entry is ``{"text", "region", "codon_start"}`` where
    ``region`` is ``"fwd"`` or ``"rev"``. Empty list when every mutation
    is safely inside the amplicon's interior.

    Why this matters: if a mutation falls inside the 5′ or 3′ binding
    window, the PCR primer won't bind perfectly to the user's original
    plasmid template — they must order the *mutated* insert as a gBlock
    and use that as the PCR template (or redesign around the site).
    """
    if not mutations or insert_len <= 0:
        return []
    positions = _codon_fix_mutation_positions(mutations)
    fwd_hi = max(0, fwd_bind_len)             # [0, fwd_hi) covers fwd binding
    rev_lo = insert_len - max(0, rev_bind_len)  # [rev_lo, insert_len) covers rev binding
    out: list[dict] = []
    for text, codon_start in zip(mutations, positions):
        if codon_start < 0:
            continue
        codon_end = codon_start + 3
        # A 3-nt codon overlaps the fwd window if any nt is in [0, fwd_hi).
        in_fwd = codon_start < fwd_hi
        # Overlaps the rev window if any nt is in [rev_lo, insert_len).
        in_rev = codon_end > rev_lo
        if in_fwd:
            out.append({"text": text, "region": "fwd",
                        "codon_start": codon_start})
        if in_rev:
            out.append({"text": text, "region": "rev",
                        "codon_start": codon_start})
    return out


def _design_gb_primers(
    template_seq: str,
    start: int,
    end: int,
    part_type: str,
    target_tm: float = 60.0,
    codon_raw: "dict | None" = None,
    grammar: "dict | None" = None,
) -> dict:
    """Design modular cloning domestication primers for a template region.

    Defaults to Golden Braid L0 (Esp3I, GGAG/TGAC/AATG/GCTT/CGCT
    overhangs); pass ``grammar`` to use a different cloning grammar
    (MoClo Plant, custom user-defined). The amplified product, after
    digestion with the grammar's enzyme, carries the 4-nt overhangs
    associated with ``part_type`` in that grammar's position table.

    Primer structure (5'→3'):

        Forward: [pad] [enzyme site] [spacer] [5' overhang]    [binding →]
        Reverse: [pad] [enzyme site] [spacer] [RC 3' overhang] [← binding RC]

    When *codon_raw* (a ``{codon: (aa, count)}`` dict from the codon-table
    registry) is supplied and *part_type* is a coding type (CDS / CDS-NS /
    C-tag), internal BsaI or Esp3I sites are silently repaired by
    substituting synonymous codons before the primers are designed — the
    returned ``insert_seq`` is then the *mutated* sequence (which is what
    the user should order as a gBlock / synthetic fragment for PCR).
    The list of substitutions made is returned under ``mutations``.

    Returns a dict with keys: part_type, position, oh5, oh3, insert_seq,
    fwd_binding, rev_binding, fwd_full, rev_full, fwd_tm, rev_tm,
    amplicon_len, mutations, and a ``pairs`` list. ``pairs`` holds one
    dict per amplicon — each with ``fwd_full``, ``rev_full``, ``fwd_tm``,
    ``rev_tm``, ``fwd_pos``, ``rev_pos``, ``fwd_binding``, ``rev_binding``,
    ``amplicon_len``. Callers that only need the first pair can continue
    to read the legacy top-level keys; multi-pair savers should iterate
    ``pairs`` (future SOE-PCR splitting for non-repairable internal sites
    will extend this list beyond one pair).
    """
    g = grammar if isinstance(grammar, dict) else _BUILTIN_GRAMMARS["gb_l0"]
    pos_spec = _grammar_position_by_type(g, part_type)
    if pos_spec is None:
        return {
            "error": f"Part type {part_type!r} is not defined in grammar "
                     f"{g.get('name', '?')}. Available types: "
                     f"{', '.join(p.get('type', '?') for p in g.get('positions', []))}.",
            "mutations": [],
        }
    pos_label, oh5, oh3 = pos_spec.get("name", "?"), pos_spec.get("oh5", ""), pos_spec.get("oh3", "")
    forbidden_sites = g.get("forbidden_sites", _GB_DOMESTICATION_FORBIDDEN) or _GB_DOMESTICATION_FORBIDDEN
    coding_types = set(g.get("coding_types", []) or _GB_CODING_PART_TYPES)
    enzyme_pad = g.get("pad", _GB_PAD)
    enzyme_site = g.get("site", _GB_L0_ENZYME_SITE)
    enzyme_spacer = g.get("spacer", _GB_SPACER)

    total  = len(template_seq)
    insert = _slice_circular(template_seq.upper(), start, end)
    wraps  = end < start

    # Need at least 18 bp to pick a proper binding region — otherwise
    # _pick_binding_region returns the whole (too-short) insert with Tm=0.
    if len(insert) < 18:
        return {
            "error": f"Cloning region is too short ({len(insert)} bp). "
                     f"Select at least 18 bp (recommended 25+ bp for a "
                     f"robust binding region).",
            "mutations": [],
        }

    # Internal Type IIS check. The grammar's `forbidden_sites` lists
    # every recognition that must be absent from the final part — for
    # GB L0 that's Esp3I (current cut) + BsaI (next-level reuse); for
    # MoClo Plant, BsaI + BpiI. Coding parts can be repaired via
    # synonymous codon substitution; non-coding parts have no reading
    # frame so internal sites must be fixed manually.
    mutations: list[str] = []
    initial_hits = _gb_find_forbidden_hits(insert, sites=forbidden_sites)
    if initial_hits:
        hit_str = ", ".join(
            f"{name} {site} at +{pos + 1}"
            for name, site, pos in initial_hits
        )
        can_attempt_fix = (
            part_type in coding_types
            and bool(codon_raw)
            and len(insert) % 3 == 0
        )
        if can_attempt_fix:
            protein = _mut_translate(insert)
            if protein:
                fixed_insert, mutations = _codon_fix_sites(
                    insert, protein, codon_raw,
                    sites=forbidden_sites,
                )
                remaining = _gb_find_forbidden_hits(
                    fixed_insert, sites=forbidden_sites,
                )
                if remaining:
                    remain_str = ", ".join(
                        f"{name} {site} at +{pos + 1}"
                        for name, site, pos in remaining
                    )
                    return {
                        "error": f"Internal Type IIS site(s) remain after "
                                 f"silent-mutation attempt ({remain_str}). "
                                 f"The sites overlap codons with no "
                                 f"synonymous alternative in this codon "
                                 f"table — pick a different region or "
                                 f"redesign.",
                        "mutations": mutations,
                    }
                insert = fixed_insert
            else:
                return {
                    "error": f"Internal Type IIS site(s) found ({hit_str}) "
                             f"but the insert could not be translated for "
                             f"silent mutation — pick a different region.",
                    "mutations": [],
                }
        else:
            reasons: list[str] = []
            if part_type not in coding_types:
                reasons.append(f"{part_type} is non-coding")
            else:
                if not codon_raw:
                    reasons.append("no codon table selected")
                if len(insert) % 3 != 0:
                    reasons.append(f"insert length {len(insert)} bp is "
                                   f"not a multiple of 3")
            extra = f" ({'; '.join(reasons)})" if reasons else ""
            return {
                "error": f"Internal Type IIS site(s) found: {hit_str}. "
                         f"Silent-mutation repair unavailable{extra}. "
                         f"Pick a different region or redesign.",
                "mutations": [],
            }

    # Forward binding: first 18-25 bp of the insert
    fwd_bind, fwd_tm = _pick_binding_region(insert, target_tm)

    # Reverse binding: first 18-25 bp of the reverse-complement of the insert
    # (i.e. the last 18-25 bp of the insert, reverse-complemented)
    rev_bind, rev_tm = _pick_binding_region(_rc(insert), target_tm)

    # Assemble full primers using the grammar's enzyme/pad/spacer.
    fwd_tail = enzyme_pad + enzyme_site + enzyme_spacer + oh5
    rev_tail = enzyme_pad + enzyme_site + enzyme_spacer + _rc(oh3)

    fwd_full = fwd_tail + fwd_bind
    rev_full = rev_tail + rev_bind

    # Amplicon = pad + Esp3I + spacer + oh + insert + oh_rc + spacer_rc
    #          + Esp3I_rc + pad_rc
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

    pair = {
        "fwd_full":     fwd_full,
        "rev_full":     rev_full,
        "fwd_binding":  fwd_bind,
        "rev_binding":  rev_bind,
        "fwd_tm":       round(fwd_tm, 1),
        "rev_tm":       round(rev_tm, 1),
        "fwd_pos":      fwd_pos,
        "rev_pos":      rev_pos,
        "amplicon_len": amplicon_len,
    }
    # Binding-region advisory: flag any silent mutation that lands inside
    # the forward or reverse primer binding window. When non-empty, the
    # user must order the mutated insert as a gBlock and use that — not
    # the original template — as the PCR template.
    binding_region_mutations = _gb_binding_region_advisory(
        mutations, len(insert), len(fwd_bind), len(rev_bind),
    )
    return {
        "part_type":    part_type,
        "position":     pos_label,
        "oh5":          oh5,
        "oh3":          oh3,
        "insert_seq":   insert,
        "mutations":    mutations,
        "binding_region_mutations": binding_region_mutations,
        "pairs":        [pair],
        # Legacy top-level mirror of pairs[0] for callers (cloning simulator,
        # PrimerDesignScreen) that don't iterate the list yet.
        **pair,
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
#       "strand":       1 | -1 | 0,             # +1 forward, -1 reverse, 0 arrowless
#       "color":        "#RRGGBB" | None,       # optional per-entry override
#       "qualifiers":   {"gene": ["lacZ"], ...}, # GenBank-compatible qualifiers
#       "description":  "(free text)"           # optional
#     }
#
# The `feature_type` is validated against a curated list of INSDC-standard
# types (see `_GENBANK_FEATURE_TYPES`). Unknown types are accepted with a
# warning — some in-house projects use non-standard types and we shouldn't
# block them, but the warning flags risk of downstream tool incompatibility.
#
# Visual styling precedence when rendering a library entry:
#   1. entry["color"] if set (per-entry user override)
#   2. user-edited type default from feature_colors.json
#   3. built-in _DEFAULT_TYPE_COLORS for the feature_type
#   4. _FEATURE_PALETTE[0] as a last-resort fallback
# Use `_resolve_feature_color(entry)` rather than reading the field directly.

_FEATURES_FILE = _DATA_DIR / "features.json"
_features_cache: "list | None" = None

# Monotonic counter bumped on every change to the feature library cache
# (writes via `_save_features`, plus first-load / cache-miss reads via
# `_load_features`). Lets long-lived screens (PartsBinModal, etc.) cache
# derived indices and detect stale state without scanning the whole
# library every populate. Strict bump-on-change — never decremented.
_features_generation: int = 0

# User-edited type → default color map persisted as entries of
# {"feature_type": "<type>", "color": "#RRGGBB"}. Kept separate from the
# entry file so the defaults survive even when the library is empty.
_FEATURE_COLORS_FILE = _DATA_DIR / "feature_colors.json"
_feature_colors_cache: "dict[str, str] | None" = None

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

# Built-in default color per feature type. Used when the user has not
# customised a default in feature_colors.json. Hex strings so they render
# identically in Rich (which accepts "#RRGGBB").
_DEFAULT_TYPE_COLORS: dict[str, str] = {
    "CDS":             "#FFA500",
    "gene":            "#FFD700",
    "mRNA":            "#FFA07A",
    "tRNA":            "#FF69B4",
    "rRNA":            "#FF1493",
    "ncRNA":           "#DA70D6",
    "misc_RNA":        "#BA55D3",
    "promoter":        "#00CED1",
    "terminator":      "#DC143C",
    "RBS":             "#00FF7F",
    "polyA_signal":    "#FF6347",
    "regulatory":      "#7FFFD4",
    "5'UTR":           "#87CEEB",
    "3'UTR":           "#4682B4",
    "intron":          "#A9A9A9",
    "exon":            "#90EE90",
    "operon":          "#DDA0DD",
    "primer_bind":     "#00BFFF",
    "protein_bind":    "#F08080",
    "misc_binding":    "#FF8C00",
    "repeat_region":   "#CD853F",
    "LTR":             "#8B4513",
    "mobile_element":  "#8B008B",
    "rep_origin":      "#9370DB",
    "oriT":            "#BA55D3",
    "sig_peptide":     "#ADFF2F",
    "mat_peptide":     "#9ACD32",
    "transit_peptide": "#7CFC00",
    "propeptide":      "#6B8E23",
    "misc_feature":    "#20B2AA",
    "misc_recomb":     "#48D1CC",
    "stem_loop":       "#FF4500",
    "variation":       "#800080",
}


def _load_features() -> list[dict]:
    """Return an independent (deep-copied) list of feature library
    entries. Callers can mutate the returned dicts freely without
    poisoning the cache — important for ``FeatureLibraryScreen``
    which buffers in-place edits (rename / color / strand / etc.) and
    then either persists or abandons. A shallow ``list(_features_cache)``
    used to share dict refs with the cache, so an abandoned mutation
    would survive in the cache and leak into the next ``_load_features``
    consumer (a freshly opened FeatureLibraryScreen, the
    DomesticatorModal feature picker, etc.) as if it had been saved.
    """
    global _features_cache, _features_generation
    from copy import deepcopy
    if _features_cache is not None:
        return deepcopy(_features_cache)
    entries, warning = _safe_load_json(_FEATURES_FILE, "Feature library")
    if warning:
        _log.warning(warning)
    entries = [e for e in entries if isinstance(e, dict)]
    _features_cache = entries
    # A fresh disk read is the result of either first-load or an
    # external invalidation (test harness setting `_features_cache =
    # None`, or a hand-edit of features.json). Either way the contents
    # may have changed since the last write, so bump the generation so
    # consumers know to rebuild any derived indices.
    _features_generation += 1
    return deepcopy(_features_cache)


def _save_features(entries: list[dict]) -> None:
    """Persist `entries` and seed the in-memory cache with a deepcopy
    so subsequent caller-side mutations of `entries` (or any dict
    inside it) cannot leak into the cache after the save returns.
    Without the deepcopy, the dicts in `_features_cache` would alias
    the dicts in the caller's list — so e.g. a FeatureLibraryScreen
    that saved, then made another change, then abandoned, would leave
    the post-save mutations stuck in the cache.
    """
    global _features_cache, _features_generation
    from copy import deepcopy
    _safe_save_json(_FEATURES_FILE, entries, "Feature library")
    _features_cache = deepcopy(entries)
    _features_generation += 1


def _build_feature_library_index() -> dict[tuple[str, str], str]:
    """Return ``{(name, feature_type): sequence_upper}`` for the entire
    feature library — single sweep, used for O(1) part-vs-feature
    lookups in the Parts Bin "Feat Lib" column. Doing this per row
    inside ``_populate`` would be O(parts × features); building once
    and looking up is O(parts + features) and reusable across the
    entire populate.
    """
    index: dict[tuple[str, str], str] = {}
    for e in _load_features():
        name = e.get("name", "")
        ftype = e.get("feature_type", "")
        if not isinstance(name, str) or not isinstance(ftype, str):
            continue
        index[(name, ftype)] = (e.get("sequence", "") or "").upper()
    return index


def _classify_feature_library_match(
    index: dict[tuple[str, str], str],
    name: str, feature_type: str, sequence: str,
) -> str:
    """Look up ``(name, feature_type)`` in a pre-built feature-library
    ``index`` and return ``"exact"`` / ``"name"`` / ``""``.

    Same semantics as :func:`_feature_library_match`, but takes a
    pre-built index so callers iterating over many parts pay the
    library-scan cost once instead of once per part. Sequences are
    compared case-insensitively (case-folded on both sides).
    """
    key = (name, feature_type)
    if key not in index:
        return ""
    existing = index[key]
    return "exact" if existing == (sequence or "").upper() else "name"


def _feature_library_match(name: str, feature_type: str,
                           sequence: str) -> str:
    """Classify whether ``(name, feature_type, sequence)`` is already in the
    feature library. Returns one of:

      - ``"exact"`` — entry with same name + type + sequence exists.
        Saving would be a no-op modulo qualifier/description tweaks.
      - ``"name"``  — entry with same name + type exists but sequence
        differs. Saving will replace the stored sequence.
      - ``""``       — no entry with this (name, type) pair.

    Convenience wrapper around :func:`_build_feature_library_index` +
    :func:`_classify_feature_library_match` for one-off lookups (e.g.
    the warning notify in ``PartsBinModal._save_as_feature``). Loops
    that scan many parts should build the index once and call
    :func:`_classify_feature_library_match` directly to avoid
    rebuilding the index on every part.
    """
    return _classify_feature_library_match(
        _build_feature_library_index(), name, feature_type, sequence,
    )


def _load_feature_colors() -> dict[str, str]:
    """Return the user's customised type → color map. Missing file / empty
    entries → empty dict. Callers should combine this with
    ``_DEFAULT_TYPE_COLORS`` — that precedence is handled by
    ``_resolve_feature_color``."""
    global _feature_colors_cache
    if _feature_colors_cache is not None:
        return dict(_feature_colors_cache)
    entries, warning = _safe_load_json(_FEATURE_COLORS_FILE, "Feature colors")
    if warning:
        _log.warning(warning)
    result: dict[str, str] = {}
    for e in entries:
        if not isinstance(e, dict):
            continue
        ft  = e.get("feature_type")
        col = e.get("color")
        if isinstance(ft, str) and ft and isinstance(col, str) and col:
            result[ft] = col
    _feature_colors_cache = dict(result)
    return dict(_feature_colors_cache)


def _save_feature_colors(mapping: dict[str, str]) -> None:
    """Persist the type → color map. Written as a list of ``{"feature_type":
    ..., "color": ...}`` dicts so it shares the schema-envelope shape with
    the other libraries (sacred invariant #7)."""
    global _feature_colors_cache
    entries = [{"feature_type": ft, "color": col}
               for ft, col in mapping.items()]
    _safe_save_json(_FEATURE_COLORS_FILE, entries, "Feature colors")
    _feature_colors_cache = dict(mapping)


def _resolve_feature_color(entry: dict) -> str:
    """Effective render color for a library entry. Never returns an empty
    string — falls back through entry override → user default → built-in
    default → palette[0]. Always produces something Rich will parse.

    Palette values that use ``color(N)`` syntax are safe inside Rich
    Style.parse but blow up Rich's markup lexer when rewrapped as
    ``[color(N)]...[/]``. Normalise them to their hex equivalent before
    returning so every downstream caller can use plain markup templating.
    """
    col = entry.get("color") if isinstance(entry, dict) else None
    if isinstance(col, str) and col:
        return _markup_safe_color(col)
    ft = entry.get("feature_type", "") if isinstance(entry, dict) else ""
    user_defaults = _load_feature_colors()
    if ft in user_defaults:
        return _markup_safe_color(user_defaults[ft])
    if ft in _DEFAULT_TYPE_COLORS:
        return _markup_safe_color(_DEFAULT_TYPE_COLORS[ft])
    return _markup_safe_color(_FEATURE_PALETTE[0])


def _markup_safe_color(col: str) -> str:
    """Return ``col`` unchanged if it's already safe for Rich markup (hex,
    named color, or ``rgb(r,g,b)``). Palette-style ``color(N)`` values
    are converted to hex via :func:`_normalise_color_input` because the
    parens collide with Rich's ``[...]`` lexer. Falls back to returning
    the raw string if normalisation fails — better to render something
    than to crash the preview."""
    if not isinstance(col, str) or not col:
        return col
    if col.startswith("color("):
        canon = _normalise_color_input(col)
        if canon:
            return canon
    return col


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
        fixed.insert(0, {
            "name":   "E. coli K12",
            "taxid":  "83333",
            "source": "builtin",
            "added":  _date.today().isoformat(),
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
    entries = _codon_tables_load()
    taxid = str(taxid or "").strip()
    name  = (name or "?").strip() or "?"
    entry = {
        "name":   name,
        "taxid":  taxid,
        "source": source,
        "added":  _date.today().isoformat(),
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


def _forbidden_hit_set(seq: str, patterns) -> set[tuple[str, int]]:
    """Return ``{(pattern, position)}`` for every occurrence of every
    pattern in *seq*. Codon swaps are 3→3 bases so positions are stable
    under the swap, which lets us compare before/after hit sets directly."""
    out: set[tuple[str, int]] = set()
    for p in patterns:
        start = 0
        while True:
            i = seq.find(p, start)
            if i == -1:
                break
            out.add((p, i))
            start = i + 1
    return out


def _codon_fix_sites(dna: str, protein: str, raw: dict,
                     sites: "dict | None" = None) -> tuple:
    """Substitute synonymous codons to remove internal restriction sites.

    ``sites`` is a forward-strand ``{name: site}`` dict; reverse complements
    are added automatically for non-palindromic sites. Returns
    ``(new_dna, fixes)``.

    Hardening (2026-04-21) — a candidate swap is accepted only if it:
      1. Actually removes the target site at the current position, AND
      2. Introduces **no new** forbidden site (forward or RC) anywhere
         in the full sequence — counted against the full input site set,
         not just the enzyme currently being iterated. This guards
         against the classic failure mode of fixing BsaI by accidentally
         spawning an Esp3I (or the RC of either) a few bases away.

    Multiple occurrences of the same site are processed left-to-right;
    each swap only needs to remove its own position, so repeated sites
    of the same enzyme are handled correctly (pre-2026-04-21 the check
    was ``site not in test`` which failed when two copies were present).
    """
    if sites is None:
        sites = _CODON_DEFAULT_FORBIDDEN
    expanded: dict = {}
    for name, site in sites.items():
        site = site.upper()
        expanded[name] = site
        rc = _mut_revcomp(site)
        if rc != site:
            expanded[f"{name}_rc"] = rc
    # Flat tuple of every forbidden pattern (forward + RC). Used by the
    # per-swap cross-check to veto swaps that would introduce a NEW
    # pattern anywhere (different enzyme, different strand, different
    # position — the check is global).
    all_forbidden = tuple(expanded.values())
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
            before_hits = _forbidden_hit_set(seq, all_forbidden)
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
                    test_seq = "".join(test)
                    after_hits = _forbidden_hit_set(test_seq, all_forbidden)
                    # (1) Target site at idx must be gone.
                    if (site, idx) in after_hits:
                        continue
                    # (2) No new forbidden hit appears anywhere.
                    #     (Existing hits elsewhere are fine — later
                    #     iterations of this loop will process them.)
                    if after_hits - before_hits:
                        continue
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
        # Current color override for this entry. None = Auto (type default).
        self._color: "str | None" = self._prefill.get("color") or None

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
                        yield Label("Orientation:")
                        with RadioSet(id="addfeat-strand"):
                            yield RadioButton("Forward (▶)",
                                              value=(strand == 1),
                                              id="addfeat-strand-fwd")
                            yield RadioButton("Reverse (◀)",
                                              value=(strand == -1),
                                              id="addfeat-strand-rev")
                            yield RadioButton("Arrowless (▒)",
                                              value=(strand == 0),
                                              id="addfeat-strand-none")
                            yield RadioButton("Double (◀▶)",
                                              value=(strand == 2),
                                              id="addfeat-strand-both")

                with Horizontal(id="addfeat-color-row"):
                    yield Label("Color:", id="addfeat-color-label")
                    yield Static("", id="addfeat-color-swatch", markup=True)
                    yield Button("Pick Color", id="btn-addfeat-color")
                    yield Button("Auto",        id="btn-addfeat-color-clear")

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
                yield Button("Import from plasmid",
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
        self._refresh_color_swatch()

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _current_feature_type(self) -> str:
        """Read the feature-type Select, falling back to the prefill's type
        (or 'misc_feature') if the Select isn't mounted yet. Used so the
        ColorPicker knows which type's default color to show in the preview."""
        try:
            val = self.query_one("#addfeat-type", Select).value
            if isinstance(val, str) and val and val != Select.BLANK:
                return val
        except NoMatches:
            pass
        return self._prefill.get("feature_type", "misc_feature") or "misc_feature"

    def _refresh_color_swatch(self) -> None:
        """Repaint the color preview after _color changes. Always resolves
        through the same precedence chain used at render time, so Auto shows
        the effective default (not an empty box).

        Uses ``Text.append`` rather than markup-string interpolation so
        palette colors like ``color(39)`` render safely — Rich's Style
        parser accepts them, but its markup lexer chokes on the parens
        inside ``[color(39)]...[/]``.
        """
        try:
            swatch = self.query_one("#addfeat-color-swatch", Static)
        except NoMatches:
            return
        ftype = self._current_feature_type()
        synth = {"feature_type": ftype, "color": self._color}
        resolved = _resolve_feature_color(synth)
        t = Text()
        t.append("███ ", style=resolved)
        if self._color:
            t.append(resolved)
        else:
            t.append(f"Auto → {resolved}", style="dim")
        swatch.update(t)

    def _gather(self) -> "dict | None":
        """Read form → entry dict, or write a red error and return None."""
        try:
            name   = self.query_one("#addfeat-name",  Input).value.strip()
            ftype  = self.query_one("#addfeat-type",  Select).value
            seqraw = self.query_one("#addfeat-seq",   TextArea).text
            quals  = self.query_one("#addfeat-quals", Input).value
            desc   = self.query_one("#addfeat-desc",  Input).value.strip()
            fwd_rb  = self.query_one("#addfeat-strand-fwd",  RadioButton)
            rev_rb  = self.query_one("#addfeat-strand-rev",  RadioButton)
            none_rb = self.query_one("#addfeat-strand-none", RadioButton)
            both_rb = self.query_one("#addfeat-strand-both", RadioButton)
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
        if   rev_rb.value:  strand = -1
        elif none_rb.value: strand = 0
        elif both_rb.value: strand = 2
        else:               strand = 1   # forward (default)
        return {
            "name":         name,
            "feature_type": str(ftype),
            "sequence":     seq_clean,
            "strand":       strand,
            "color":        self._color,
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

    @on(Button.Pressed, "#btn-addfeat-color")
    def _open_color_picker(self, _) -> None:
        """Push ColorPickerModal; the dismiss payload updates self._color."""
        ftype = self._current_feature_type()

        def _on_color(result) -> None:
            if not isinstance(result, dict):
                return
            new_col = result.get("color")
            self._color = new_col if isinstance(new_col, str) and new_col else None
            if result.get("set_default") and isinstance(new_col, str) and new_col:
                defaults = _load_feature_colors()
                defaults[ftype] = new_col
                try:
                    _save_feature_colors(defaults)
                except (OSError, ValueError) as exc:
                    _log.exception("Saving type-default color failed")
                    self.notify(f"Saving default failed: {exc}",
                                severity="error")
            self._refresh_color_swatch()

        self.app.push_screen(ColorPickerModal(ftype, self._color),
                             callback=_on_color)

    @on(Button.Pressed, "#btn-addfeat-color-clear")
    def _clear_color(self, _) -> None:
        self._color = None
        self._refresh_color_swatch()

    @on(Select.Changed, "#addfeat-type")
    def _type_changed(self, _) -> None:
        # Auto-mode swatch tracks the effective default for the chosen type.
        self._refresh_color_swatch()

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
            fwd_rb   = self.query_one("#addfeat-strand-fwd",  RadioButton)
            rev_rb   = self.query_one("#addfeat-strand-rev",  RadioButton)
            none_rb  = self.query_one("#addfeat-strand-none", RadioButton)
            both_rb  = self.query_one("#addfeat-strand-both", RadioButton)
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
        strand = entry.get("strand", 1)
        fwd_rb.value  = (strand == 1)
        rev_rb.value  = (strand == -1)
        none_rb.value = (strand == 0)
        both_rb.value = (strand == 2)
        col = entry.get("color")
        self._color = col if isinstance(col, str) and col else None
        self._refresh_color_swatch()
        status.update(
            f"[green]Imported '{entry.get('name', '?')}' — "
            f"review and Save or Insert.[/green]"
        )

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Color picker modal ────────────────────────────────────────────────────────

_HEX3_RE = re.compile(r"^#[0-9A-Fa-f]{3}$")
_HEX6_RE = re.compile(r"^#[0-9A-Fa-f]{6}$")
_XTERM_RE = re.compile(r"^(?:color\()?(\d{1,3})\)?$")


def _normalise_color_input(raw: str) -> "str | None":
    """Accept any of ``#RGB`` / ``#RRGGBB`` / ``0..255`` / ``color(N)`` and
    return a canonical hex string. Returns ``None`` if the input isn't
    parseable — callers surface a validation error.

    Hex is preferred over ``color(N)`` on save because hex survives round-
    trip through JSON unambiguously and renders identically on every
    terminal that supports truecolor. xterm indices are converted to their
    24-bit RGB equivalent via the standard 6×6×6 cube / grayscale ramp.
    """
    if not isinstance(raw, str):
        return None
    raw = raw.strip()
    if not raw:
        return None
    if _HEX6_RE.match(raw):
        return "#" + raw[1:].upper()
    if _HEX3_RE.match(raw):
        r, g, b = raw[1], raw[2], raw[3]
        return f"#{r}{r}{g}{g}{b}{b}".upper()
    m = _XTERM_RE.match(raw)
    if m:
        try:
            idx = int(m.group(1))
        except ValueError:
            return None
        if 0 <= idx <= 255:
            return _xterm_index_to_hex(idx)
    return None


# Standard 16 ANSI colors as rendered by most terminals. These are the only
# colors guaranteed on an 8/16-color terminal; truecolor hex will be
# approximated to the nearest ANSI on such terminals by Rich.
_ANSI16_HEX: list[str] = [
    "#000000", "#800000", "#008000", "#808000",
    "#000080", "#800080", "#008080", "#C0C0C0",
    "#808080", "#FF0000", "#00FF00", "#FFFF00",
    "#0000FF", "#FF00FF", "#00FFFF", "#FFFFFF",
]


def _xterm_index_to_hex(idx: int) -> str:
    """Convert an xterm-256 color index (0..255) to the closest 24-bit RGB
    hex. Matches the xterm default palette — terminals may remap these but
    the vast majority follow the spec. Cube levels use the canonical
    ``[0, 95, 135, 175, 215, 255]`` ramp; grayscale uses
    ``8 + 10 * k`` for k in 0..23."""
    idx = max(0, min(255, int(idx)))
    if idx < 16:
        return _ANSI16_HEX[idx]
    if idx < 232:
        n = idx - 16
        levels = (0, 95, 135, 175, 215, 255)
        r = levels[(n // 36) % 6]
        g = levels[(n // 6)  % 6]
        b = levels[ n        % 6]
        return f"#{r:02X}{g:02X}{b:02X}"
    v = 8 + 10 * (idx - 232)
    return f"#{v:02X}{v:02X}{v:02X}"


class ColorPickerModal(ModalScreen):
    """Pick a display color for a feature-library entry.

    Returns via ``dismiss``:
      * ``None`` — user cancelled.
      * ``{"color": "#RRGGBB", "set_default": False}`` — set entry color.
      * ``{"color": "#RRGGBB", "set_default": True}`` — also save as the
        default for this feature type.
      * ``{"color": None, "set_default": False}`` — clear the override and
        fall back to the type default.

    Three ways to pick a color:

      1. **Curated quick-picks** — the 20-color ``_SWATCHES`` that reuse
         the main map palette.
      2. **xterm 256-color grid** — full 256-cell grid (16 ANSI + 216 cube
         + 24 grayscale). Renders as tiny colored buttons; click one to
         load into the preview.
      3. **Custom hex / index input** — free-form ``#RGB`` / ``#RRGGBB`` /
         ``0..255`` / ``color(N)``. Validated via ``_normalise_color_input``
         which also converts xterm indices to their RGB equivalent so the
         stored value is always a canonical uppercase hex string.

    If the terminal only supports 8/16 colors (``console.color_system`` is
    ``"standard"`` or ``None``), a yellow warning explains that truecolor
    choices will be approximated. The picker still works — you just can't
    visually distinguish similar hex colors on that terminal.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]

    # Curated quick-picks — hex-encoded so they round-trip through JSON.
    _SWATCHES: list[str] = [
        "#FF6347", "#FFA500", "#FFD700", "#FFFF00", "#ADFF2F",
        "#7CFC00", "#00FF7F", "#00CED1", "#1E90FF", "#4169E1",
        "#9370DB", "#BA55D3", "#FF69B4", "#FF1493", "#DC143C",
        "#A0522D", "#CD853F", "#20B2AA", "#708090", "#2F4F4F",
    ]

    def __init__(self, feature_type: str, current_color: "str | None") -> None:
        super().__init__()
        self._feature_type = feature_type
        self._current      = current_color or ""
        self._pending:  "str | None" = current_color or None
        self._drag_active: bool = False

    def compose(self) -> ComposeResult:
        type_default = _DEFAULT_TYPE_COLORS.get(self._feature_type, "")
        user_default = _load_feature_colors().get(self._feature_type, "")
        effective_default = _markup_safe_color(
            user_default or type_default or "#808080")

        with Vertical(id="colorpick-dlg"):
            yield Static(f" Pick color for {self._feature_type} ",
                         id="colorpick-title")
            yield Label(
                f"Current: {self._current or '[dim](auto — using type default)[/]'}   "
                f"Type default: [{effective_default}]███[/]",
                markup=True, id="colorpick-current",
            )
            yield Static("", id="colorpick-capability", markup=True)

            with Horizontal(id="colorpick-preview-row"):
                yield Static("", id="colorpick-preview-swatch")
                yield Static("", id="colorpick-preview-label", markup=True)

            with ScrollableContainer(id="colorpick-scroll"):
                yield Label("Curated", classes="colorpick-section-hdr")
                with Horizontal(id="colorpick-row"):
                    for i, hex_col in enumerate(self._SWATCHES):
                        yield Button("  ", id=f"colorpick-swatch-{i}",
                                     classes="colorpick-swatch")

                yield Label("xterm 256  (click any cell)",
                            classes="colorpick-section-hdr")
                with Vertical(id="colorpick-xterm-grid"):
                    # 16 ANSI colors — one row
                    with Horizontal(classes="colorpick-xterm-row"):
                        for i in range(16):
                            yield Button(" ", id=f"colorpick-x-{i}",
                                         classes="colorpick-xterm-cell")
                    # 216-color cube — 6 rows × 36 cols
                    for row in range(6):
                        with Horizontal(classes="colorpick-xterm-row"):
                            for col in range(36):
                                idx = 16 + row * 36 + col
                                yield Button(" ", id=f"colorpick-x-{idx}",
                                             classes="colorpick-xterm-cell")
                    # 24-color grayscale — one row
                    with Horizontal(classes="colorpick-xterm-row"):
                        for i in range(232, 256):
                            yield Button(" ", id=f"colorpick-x-{i}",
                                         classes="colorpick-xterm-cell")

                yield Label("Custom", classes="colorpick-section-hdr")
                with Horizontal(id="colorpick-custom-row"):
                    yield Label("Hex / xterm idx:",
                                id="colorpick-custom-label")
                    yield Input(
                        value=(self._current if _HEX6_RE.match(self._current)
                               else ""),
                        placeholder="#FF6347, F63, 208, or color(208)",
                        id="colorpick-hex-input",
                    )
                    yield Button("Apply", id="btn-colorpick-apply",
                                 variant="primary")

            yield Static("", id="colorpick-status", markup=True)
            with Horizontal(id="colorpick-btns"):
                yield Button("Auto (clear override)",
                             id="btn-colorpick-auto")
                yield Button("Save",
                             id="btn-colorpick-save",
                             variant="primary")
                yield Button("Save + set as type default",
                             id="btn-colorpick-default",
                             variant="success")
                yield Button("Cancel", id="btn-colorpick-cancel")

    def on_mount(self) -> None:
        # Paint curated swatches with their own background so the palette
        # is visible at a glance.
        for i, hex_col in enumerate(self._SWATCHES):
            try:
                btn = self.query_one(f"#colorpick-swatch-{i}", Button)
            except NoMatches:
                continue
            btn.styles.background = hex_col
        # Paint each xterm cell with its index's RGB value. If the terminal
        # is 8/16-color, Rich will downsample — that's expected, and the
        # capability warning already explains it.
        for idx in range(256):
            try:
                btn = self.query_one(f"#colorpick-x-{idx}", Button)
            except NoMatches:
                continue
            btn.styles.background = _xterm_index_to_hex(idx)
        self._refresh_capability_warning()
        self._refresh_status()

    def _refresh_capability_warning(self) -> None:
        """Warn the user if the terminal is 8/16-color — they can still
        pick truecolor hex, it'll just be approximated to the nearest
        ANSI when rendered."""
        try:
            cap = self.query_one("#colorpick-capability", Static)
        except NoMatches:
            return
        try:
            sys_name = (self.app.console.color_system or "").lower()
        except Exception:       # noqa: BLE001
            sys_name = ""
        if sys_name in ("truecolor", "256"):
            cap.update(
                f"[dim]Terminal palette: {sys_name} — full range available.[/]"
            )
        else:
            label = sys_name or "unknown / 8-color"
            cap.update(
                f"[yellow]Terminal palette: {label}. Truecolor choices "
                f"will be approximated to the nearest ANSI color.[/]"
            )

    def _refresh_status(self) -> None:
        """Repaint the big preview swatch + hex label. Called whenever
        ``self._pending`` changes — including during a live drag across
        xterm cells."""
        try:
            swatch = self.query_one("#colorpick-preview-swatch", Static)
            label  = self.query_one("#colorpick-preview-label",  Static)
        except NoMatches:
            return
        if self._pending:
            swatch.styles.background = self._pending
            t = Text()
            t.append("Selected: ", style="bold")
            t.append(self._pending, style=self._pending)
            label.update(t)
        else:
            swatch.styles.background = "transparent"
            label.update("[dim]Selected: Auto (use type default)[/]")

    def _set_pending(self, value: "str | None") -> None:
        """Central entry point for any source that changes the pending
        color — keeps the preview swatch in lock-step and clears any
        stale error message in #colorpick-status."""
        if value == self._pending:
            return
        self._pending = value
        self._refresh_status()
        try:
            self.query_one("#colorpick-status", Static).update("")
        except NoMatches:
            pass

    def _cell_index_at(self, sx: int, sy: int) -> "int | None":
        """Hit-test screen coords against the xterm grid. Returns the
        xterm index (0..255) under the cursor, or ``None`` if the point
        is outside any cell. Used by the live-drag preview so the user
        can sweep across the grid with a held click."""
        try:
            widget, _ = self.get_widget_at(sx, sy)
        except Exception:  # noqa: BLE001 — NoWidget + anything defensive
            return None
        btn_id = getattr(widget, "id", None) or ""
        if not btn_id.startswith("colorpick-x-"):
            return None
        try:
            idx = int(btn_id.rsplit("-", 1)[1])
        except ValueError:
            return None
        if 0 <= idx <= 255:
            return idx
        return None

    def on_mouse_down(self, event: MouseDown) -> None:
        """Entering a drag: if the mouse-down lands on an xterm cell,
        arm drag-mode and load that cell's color into the preview
        immediately. Other targets are left alone so regular button
        clicks (Save, Cancel, Apply) still work."""
        if event.button != 1:
            return
        idx = self._cell_index_at(event.screen_x, event.screen_y)
        if idx is not None:
            self._drag_active = True
            self._set_pending(_xterm_index_to_hex(idx))

    def on_mouse_move(self, event: MouseMove) -> None:
        """During a drag, follow the cursor across xterm cells and keep
        the preview in sync. No effect outside of drag-mode."""
        if not self._drag_active:
            return
        idx = self._cell_index_at(event.screen_x, event.screen_y)
        if idx is not None:
            self._set_pending(_xterm_index_to_hex(idx))

    def on_mouse_up(self, event: MouseUp) -> None:
        self._drag_active = False

    @on(Button.Pressed, ".colorpick-swatch")
    def _swatch(self, event: Button.Pressed) -> None:
        btn_id = event.button.id or ""
        if not btn_id.startswith("colorpick-swatch-"):
            return
        try:
            idx = int(btn_id.rsplit("-", 1)[1])
        except ValueError:
            return
        if 0 <= idx < len(self._SWATCHES):
            self._set_pending(self._SWATCHES[idx])

    @on(Button.Pressed, ".colorpick-xterm-cell")
    def _xterm_cell(self, event: Button.Pressed) -> None:
        btn_id = event.button.id or ""
        if not btn_id.startswith("colorpick-x-"):
            return
        try:
            idx = int(btn_id.rsplit("-", 1)[1])
        except ValueError:
            return
        if 0 <= idx <= 255:
            self._set_pending(_xterm_index_to_hex(idx))

    @on(Button.Pressed, "#btn-colorpick-apply")
    def _apply_custom(self, _) -> None:
        try:
            inp = self.query_one("#colorpick-hex-input", Input)
            status = self.query_one("#colorpick-status", Static)
        except NoMatches:
            return
        raw = inp.value.strip()
        canonical = _normalise_color_input(raw)
        if canonical is None:
            status.update(
                f"[red]Invalid color '{raw}' — use #RGB, #RRGGBB, or 0..255.[/]"
            )
            return
        self._set_pending(canonical)

    @on(Input.Submitted, "#colorpick-hex-input")
    def _hex_submitted(self, _) -> None:
        self._apply_custom(None)

    @on(Button.Pressed, "#btn-colorpick-auto")
    def _auto(self, _) -> None:
        self._set_pending(None)

    @on(Button.Pressed, "#btn-colorpick-save")
    def _save(self, _) -> None:
        self.dismiss({"color": self._pending, "set_default": False})

    @on(Button.Pressed, "#btn-colorpick-default")
    def _save_default(self, _) -> None:
        if not self._pending:
            self.query_one("#colorpick-status", Static).update(
                "[yellow]Pick a specific color before setting as default.[/]"
            )
            return
        self.dismiss({"color": self._pending, "set_default": True})

    @on(Button.Pressed, "#btn-colorpick-cancel")
    def _cancel(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Feature library workbench (full-screen) ───────────────────────────────────

class _FeatureSnippetPanel(Widget):
    """Visualization of a single library entry: header line (name/type/
    strand/length/color swatch), a double-stranded DNA block rendered via
    the shared ``_build_seq_text`` pipeline (so the preview exactly matches
    the main SequencePanel — dithered ▒ bar + directional arrowhead), and
    a qualifiers list.

    Not interactive — re-render by calling ``show(entry)``.
    """

    DEFAULT_CSS = """
    _FeatureSnippetPanel {
        height: 1fr;
        overflow-y: auto;
        padding: 1 2;
    }
    _FeatureSnippetPanel > #snippet-header,
    _FeatureSnippetPanel > #snippet-quals {
        width: 100%;
        height: auto;
    }
    _FeatureSnippetPanel > #snippet-dna {
        width: 100%;
        height: auto;
        margin: 1 0;
    }
    """

    def __init__(self) -> None:
        super().__init__()
        self._entry: "dict | None" = None

    def compose(self) -> ComposeResult:
        yield Static("", id="snippet-header", markup=True)
        yield Static("", id="snippet-dna", markup=False)
        yield Static("", id="snippet-quals", markup=True)

    def show(self, entry: "dict | None") -> None:
        self._entry = entry
        try:
            hdr  = self.query_one("#snippet-header", Static)
            dna  = self.query_one("#snippet-dna",    Static)
            qual = self.query_one("#snippet-quals",  Static)
        except NoMatches:
            return
        if entry is None:
            hdr.update("[dim]Select a feature on the left to view it.[/]")
            dna.update("")
            qual.update("")
            return
        hdr.update(self._format_header(entry))
        dna.update(self._render_dna(entry))
        qual.update(self._format_qualifiers(entry))

    @staticmethod
    def _format_header(entry: dict) -> str:
        name   = entry.get("name", "?") or "?"
        ftype  = entry.get("feature_type", "?") or "?"
        seq    = (entry.get("sequence", "") or "").upper()
        strand = entry.get("strand", 1)
        desc   = entry.get("description", "") or ""
        color  = _resolve_feature_color(entry)
        length = len(seq)
        strand_tag = {1: "→ forward", -1: "← reverse",
                      0: "· arrowless", 2: "↔ double"}.get(strand, "→ forward")
        lines: list[str] = [
            f"[bold {color}]{name}[/]   "
            f"[dim]type:[/] [{color}]{ftype}[/]   "
            f"[dim]strand:[/] {strand_tag}   "
            f"[dim]length:[/] {length} bp   "
            f"[dim]color:[/] [{color}]███[/] {color}",
        ]
        if desc:
            lines.append(f"[dim]{desc}[/]")
        return "\n".join(lines)

    @staticmethod
    def _render_dna(entry: "dict") -> "Text":
        """Build a full double-stranded DNA visualization by synthesizing a
        single full-span feature and running it through ``_build_seq_text`` —
        the exact same pipeline the main SequencePanel uses. The feature's
        ``strand`` decides which arrowhead (if any) is drawn."""
        seq = (entry.get("sequence", "") or "").upper()
        if not seq:
            return Text("(no sequence)", style="dim")
        name   = entry.get("name", "") or ""
        ftype  = entry.get("feature_type", "misc_feature") or "misc_feature"
        strand = entry.get("strand", 1)
        color  = _resolve_feature_color(entry)
        synth = [{
            "type":   ftype,
            "start":  0,
            "end":    len(seq),
            "strand": strand,
            "color":  color,
            "label":  name,
        }]
        return _build_seq_text(seq, synth, line_width=60)

    @staticmethod
    def _format_qualifiers(entry: dict) -> str:
        quals = entry.get("qualifiers") or {}
        if not quals:
            return "[dim]No qualifiers.[/]"
        lines = ["[bold]Qualifiers[/]"]
        for k in sorted(quals.keys()):
            vals = quals.get(k) or []
            if isinstance(vals, list):
                val_s = "; ".join(str(v) for v in vals)
            else:
                val_s = str(vals)
            lines.append(f"  [cyan]{k}[/] = {val_s}")
        return "\n".join(lines)


class FeatureLibraryScreen(Screen):
    """Full-screen library browser for persistent feature entries.

    Left column: DataTable listing every entry (Name / Type / Strand / bp
    / Color swatch). Right column: ``_FeatureSnippetPanel`` visualizing the
    selected entry. Bottom buttons handle CRUD + styling.

    Edits are buffered in ``self._entries`` and only written to disk by
    ``action_save`` (Save button or Ctrl+S). Dirty entries get an
    ``*`` prefix in the table; the title bar shows ``*`` when there's
    anything pending (covers deletions that leave no row to flag).
    Closing with pending changes triggers ``UnsavedQuitModal`` so the
    user can choose Save/Abandon/Cancel rather than silently losing
    work. Routes through ``_load_features`` / ``_save_features`` which
    enforce the schema envelope (sacred invariant #7).
    """

    BINDINGS = [
        Binding("escape", "close",     "Close"),
        Binding("a",      "add",       "Add"),
        Binding("e",      "edit",      "Edit"),
        Binding("r",      "rename",    "Rename"),
        Binding("d",      "duplicate", "Duplicate"),
        Binding("delete", "remove",    "Remove"),
        Binding("c",      "color",     "Color"),
        Binding("s",      "strand",    "Cycle Strand"),
        Binding("ctrl+s", "save",      "Save"),
    ]

    def check_action(self, action: str, parameters: tuple) -> bool | None:
        return True

    def __init__(self) -> None:
        super().__init__()
        self._entries: list[dict] = list(_load_features())
        self._selected_index: int = 0 if self._entries else -1
        # Dirty tracking. `_dirty_indices` flags entries with unsaved
        # field changes (asterisk prefix in the table). It can't cover
        # deletions because the row no longer exists, so
        # `_has_pending_changes` is the umbrella signal used by
        # `action_close` to decide whether to prompt.
        self._dirty_indices: set[int] = set()
        self._has_pending_changes: bool = False

    def compose(self) -> ComposeResult:
        yield Header()
        with Vertical(id="flib-box"):
            yield Static(" Feature Library ", id="flib-title")
            with Horizontal(id="flib-main"):
                with Vertical(id="flib-left"):
                    yield Static("Entries", classes="flib-section-hdr")
                    yield DataTable(id="flib-table",
                                    cursor_type="row", zebra_stripes=True)
                with Vertical(id="flib-right"):
                    yield Static("Preview", classes="flib-section-hdr")
                    yield _FeatureSnippetPanel()
            with Horizontal(id="flib-btns"):
                yield Button("Add",             id="btn-flib-add")
                yield Button("Edit",            id="btn-flib-edit")
                yield Button("Rename",          id="btn-flib-rename")
                yield Button("Duplicate",       id="btn-flib-dup")
                yield Button("Remove",          id="btn-flib-remove",
                             variant="error")
                yield Button("Color",           id="btn-flib-color")
                yield Button("Cycle Strand",    id="btn-flib-strand")
                yield Button("Export FASTA",    id="btn-flib-export-fasta")
                yield Button("Save",            id="btn-flib-save",
                             variant="primary")
                yield Button("Close  [Esc]",    id="btn-flib-close")
        yield Footer()

    def on_mount(self) -> None:
        tbl = self.query_one("#flib-table", DataTable)
        tbl.add_columns("Name", "Type", "±", "bp", "Color")
        self._repopulate_table()

    # ── rendering ────────────────────────────────────────────────────────────

    def _repopulate_table(self) -> None:
        try:
            tbl = self.query_one("#flib-table", DataTable)
        except NoMatches:
            return
        tbl.clear(columns=False)
        for i, entry in enumerate(self._entries):
            color = _resolve_feature_color(entry)
            strand = entry.get("strand", 1)
            strand_tag = {1: "+", -1: "−", 0: "·", 2: "↔"}.get(strand, "+")
            bp = len((entry.get("sequence") or ""))
            # Use Rich Text for the Color cell so the swatch actually tints
            swatch = Text("███ ", style=color)
            swatch.append(color, style="dim")
            base_name = entry.get("name", "?")
            name_cell = ("*" + base_name) if i in self._dirty_indices else base_name
            tbl.add_row(
                name_cell,
                entry.get("feature_type", "?"),
                strand_tag,
                str(bp),
                swatch,
            )
        # Keep the selection in range.
        if self._entries:
            if self._selected_index < 0 or self._selected_index >= len(self._entries):
                self._selected_index = 0
            tbl.move_cursor(row=self._selected_index)
        else:
            self._selected_index = -1
        self._refresh_title()
        self._refresh_preview()

    def _refresh_title(self) -> None:
        try:
            title = self.query_one("#flib-title", Static)
        except NoMatches:
            return
        title.update(
            " Feature Library *" if self._has_pending_changes
            else " Feature Library "
        )

    def _refresh_preview(self) -> None:
        try:
            snip = self.query_one(_FeatureSnippetPanel)
        except NoMatches:
            return
        if 0 <= self._selected_index < len(self._entries):
            snip.show(self._entries[self._selected_index])
        else:
            snip.show(None)

    # ── events ───────────────────────────────────────────────────────────────

    @on(DataTable.RowHighlighted, "#flib-table")
    def _row_highlighted(self, event: DataTable.RowHighlighted) -> None:
        try:
            tbl = self.query_one("#flib-table", DataTable)
        except NoMatches:
            return
        self._selected_index = tbl.cursor_row
        self._refresh_preview()

    @on(Button.Pressed, "#btn-flib-export-fasta")
    def _export_fasta_btn(self, _) -> None: self.action_export_fasta()

    def action_export_fasta(self) -> None:
        """Export the highlighted library entry as single-record FASTA."""
        entry = self._current()
        if entry is None:
            self.app.notify("Select a feature first.", severity="warning")
            return
        seq = entry.get("sequence") or ""
        if not seq:
            self.app.notify(
                "This entry has no sequence to export.",
                severity="warning",
            )
            return
        name = entry.get("name") or "feature"
        ftype = entry.get("feature_type") or "?"
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", name) or "feature"
        default_path = str(Path.home() / f"{safe}.fa")

        def _on_done(summary):
            if not summary:
                return
            self.app.notify(
                f"Exported '{summary['name']}' to {summary['path']} "
                f"({summary['bp']} bp).",
            )

        self.app.push_screen(
            FastaExportModal(
                name=name,
                sequence=seq,
                default_path=default_path,
                subtitle=f"[{name}]  [{ftype}]  {len(seq)} bp",
            ),
            callback=_on_done,
        )

    @on(Button.Pressed, "#btn-flib-close")
    def _close_btn(self, _) -> None:
        self.action_close()

    def action_close(self) -> None:
        """Pop back to the main app, prompting on unsaved changes.

        With pending edits, push ``UnsavedQuitModal`` (the same dialog
        the main quit path uses) and wait for the user's choice:
        ``"save"`` → persist + pop, ``"abandon"`` → pop without saving,
        ``None`` → stay open. Without changes, pop immediately.
        """
        if not self._has_pending_changes:
            self.app.pop_screen()
            return

        def _cb(result):
            if result == "save":
                if self._persist_all():
                    self.app.pop_screen()
                # On save failure stay open so the user can retry.
            elif result == "abandon":
                self.app.pop_screen()
            # None → cancel; do nothing.

        self.app.push_screen(UnsavedQuitModal(), callback=_cb)

    @on(Button.Pressed, "#btn-flib-save")
    def _save_btn(self, _) -> None: self.action_save()

    def action_save(self) -> None:
        """Persist self._entries → features.json. No-op if nothing pending."""
        if not self._has_pending_changes:
            self.app.notify("No changes to save.", severity="information")
            return
        if self._persist_all():
            self.app.notify(
                f"Saved {len(self._entries)} feature(s) to library."
            )

    # ── persistence helpers ──────────────────────────────────────────────────

    def _persist_all(self) -> bool:
        """Write self._entries → features.json and clear the dirty marks.
        Returns True on success."""
        try:
            _save_features(self._entries)
        except (OSError, ValueError) as exc:
            _log.exception("Feature library save failed")
            self.app.notify(f"Save failed: {exc}", severity="error")
            return False
        self._dirty_indices.clear()
        self._has_pending_changes = False
        self._repopulate_table()
        return True

    def _mark_dirty(self, idx: int) -> None:
        """Tag an entry as having unsaved field changes (asterisk in the
        table) and set the umbrella pending-changes flag.
        """
        if 0 <= idx < len(self._entries):
            self._dirty_indices.add(idx)
        self._has_pending_changes = True

    def _shift_dirty_after_remove(self, removed_idx: int) -> None:
        """Rewrite ``_dirty_indices`` after deleting entry ``removed_idx``.

        Drops ``removed_idx`` itself, shifts every higher index down by
        one. Without this, asterisks would stick to the wrong row after
        a delete (or, worse, point past the end of the list).
        """
        new_dirty: set[int] = set()
        for i in self._dirty_indices:
            if i < removed_idx:
                new_dirty.add(i)
            elif i > removed_idx:
                new_dirty.add(i - 1)
        self._dirty_indices = new_dirty

    def _current(self) -> "dict | None":
        if 0 <= self._selected_index < len(self._entries):
            return self._entries[self._selected_index]
        return None

    # ── actions ──────────────────────────────────────────────────────────────

    @on(Button.Pressed, "#btn-flib-add")
    def _add_btn(self, _) -> None: self.action_add()

    def action_add(self) -> None:
        def _cb(result):
            if not result:
                return
            entry = result.get("entry") if isinstance(result, dict) else None
            if not entry:
                return
            self._upsert_entry(entry, notice="Added")
        self.app.push_screen(AddFeatureModal(have_cursor=False), callback=_cb)

    @on(Button.Pressed, "#btn-flib-edit")
    def _edit_btn(self, _) -> None: self.action_edit()

    def action_edit(self) -> None:
        """Open AddFeatureModal pre-filled with the current entry. The
        modal already round-trips its prefill (`_apply_prefill` mirrors
        `_gather`) so editing is just "Add with the existing dict."
        Replaces the entry at the current index on save and marks it
        dirty; doesn't write to disk until Save.
        """
        entry = self._current()
        if entry is None:
            self.app.notify("Select a feature first.", severity="warning")
            return
        target_idx = self._selected_index

        def _cb(result):
            if not result:
                return
            new_entry = result.get("entry") if isinstance(result, dict) else None
            if not new_entry:
                return
            # If the user kept the same (name, feature_type), this is a
            # plain edit-in-place. If they changed the name/type to one
            # that already exists at a different index, dedup that one
            # too (Add-style "latest write wins").
            self._replace_entry(target_idx, new_entry)

        self.app.push_screen(
            AddFeatureModal(prefill=entry, have_cursor=False),
            callback=_cb,
        )

    def _upsert_entry(self, entry: dict, notice: str) -> None:
        """Append ``entry`` (or replace the entry with the same
        (name, feature_type) key — "latest write wins"). Marks the new
        index dirty and updates the selection.
        """
        key = (entry.get("name"), entry.get("feature_type"))
        existing_idx = next(
            (i for i, e in enumerate(self._entries)
             if (e.get("name"), e.get("feature_type")) == key),
            -1,
        )
        if existing_idx >= 0:
            self._entries[existing_idx] = entry
            self._selected_index = existing_idx
            self._mark_dirty(existing_idx)
        else:
            self._entries.append(entry)
            self._selected_index = len(self._entries) - 1
            self._mark_dirty(self._selected_index)
        self._repopulate_table()
        self.app.notify(f"{notice} '{entry.get('name')}' (unsaved).")

    def _replace_entry(self, target_idx: int, new_entry: dict) -> None:
        """Replace entry at ``target_idx``, deduping any other entry that
        ends up sharing the new (name, feature_type) key.
        """
        if not (0 <= target_idx < len(self._entries)):
            return
        new_key = (new_entry.get("name"), new_entry.get("feature_type"))
        # Drop any OTHER entry that now collides with the new key.
        for i in range(len(self._entries) - 1, -1, -1):
            if i == target_idx:
                continue
            e = self._entries[i]
            if (e.get("name"), e.get("feature_type")) == new_key:
                del self._entries[i]
                self._shift_dirty_after_remove(i)
                if i < target_idx:
                    target_idx -= 1
        self._entries[target_idx] = new_entry
        self._selected_index = target_idx
        self._mark_dirty(target_idx)
        self._repopulate_table()
        self.app.notify(f"Edited '{new_entry.get('name')}' (unsaved).")

    @on(Button.Pressed, "#btn-flib-rename")
    def _rename_btn(self, _) -> None: self.action_rename()

    def action_rename(self) -> None:
        entry = self._current()
        if entry is None:
            self.app.notify("Select a feature first.", severity="warning")
            return
        old = entry.get("name", "")
        idx = self._selected_index

        def _cb(new_name):
            if not new_name:
                return
            if new_name == old:
                return
            entry["name"] = str(new_name)
            self._mark_dirty(idx)
            self._repopulate_table()
            self.app.notify(f"Renamed '{old}' → '{new_name}' (unsaved).")

        self.app.push_screen(RenamePlasmidModal(old, ""), callback=_cb)

    @on(Button.Pressed, "#btn-flib-dup")
    def _dup_btn(self, _) -> None: self.action_duplicate()

    def action_duplicate(self) -> None:
        entry = self._current()
        if entry is None:
            self.app.notify("Select a feature first.", severity="warning")
            return
        import copy as _copy
        dup = _copy.deepcopy(entry)
        # Suffix a unique "(copy)" / "(copy N)" name so we don't dedup it out.
        base = dup.get("name", "feature")
        existing_names = {e.get("name") for e in self._entries}
        cand = f"{base} (copy)"
        n = 2
        while cand in existing_names:
            cand = f"{base} (copy {n})"
            n += 1
        dup["name"] = cand
        self._entries.append(dup)
        self._selected_index = len(self._entries) - 1
        self._mark_dirty(self._selected_index)
        self._repopulate_table()
        self.app.notify(f"Duplicated as '{cand}' (unsaved).")

    @on(Button.Pressed, "#btn-flib-remove")
    def _remove_btn(self, _) -> None: self.action_remove()

    def action_remove(self) -> None:
        entry = self._current()
        if entry is None:
            return
        name = entry.get("name", "?")
        removed_idx = self._selected_index
        del self._entries[removed_idx]
        self._shift_dirty_after_remove(removed_idx)
        self._has_pending_changes = True
        if self._selected_index >= len(self._entries):
            self._selected_index = len(self._entries) - 1
        self._repopulate_table()
        self.app.notify(f"Removed '{name}' (unsaved).")

    @on(Button.Pressed, "#btn-flib-color")
    def _color_btn(self, _) -> None: self.action_color()

    def action_color(self) -> None:
        entry = self._current()
        if entry is None:
            self.app.notify("Select a feature first.", severity="warning")
            return
        ftype = entry.get("feature_type", "")
        current = entry.get("color")
        idx = self._selected_index

        def _cb(result):
            if not result:
                return
            new_color = result.get("color")
            set_default = bool(result.get("set_default"))
            entry["color"] = new_color   # None → auto
            # User-defaults map (feature_colors.json) is a separate
            # file from the per-entry feature library — saving it
            # immediately is correct; the deferred-save model only
            # applies to features.json itself.
            if set_default and isinstance(new_color, str) and new_color:
                defaults = _load_feature_colors()
                defaults[ftype] = new_color
                try:
                    _save_feature_colors(defaults)
                except (OSError, ValueError) as exc:
                    _log.exception("Feature color default save failed")
                    self.app.notify(f"Save default failed: {exc}",
                                    severity="error")
            self._mark_dirty(idx)
            self._repopulate_table()
            shown = new_color if new_color else "auto"
            self.app.notify(f"Color set to {shown} (unsaved).")

        self.app.push_screen(ColorPickerModal(ftype, current), callback=_cb)

    @on(Button.Pressed, "#btn-flib-strand")
    def _strand_btn(self, _) -> None: self.action_strand()

    def action_strand(self) -> None:
        """Cycle strand direction: +1 (▶) → -1 (◀) → 0 (arrowless ▒) →
        2 (double ◀▒▶) → +1."""
        entry = self._current()
        if entry is None:
            self.app.notify("Select a feature first.", severity="warning")
            return
        cur = entry.get("strand", 1)
        nxt = {1: -1, -1: 0, 0: 2, 2: 1}.get(cur, 1)
        entry["strand"] = nxt
        self._mark_dirty(self._selected_index)
        self._repopulate_table()
        tag = {1:  "forward (→)",
               -1: "reverse (←)",
               0:  "arrowless (·)",
               2:  "double (↔)"}.get(nxt, "+")
        self.app.notify(f"Strand → {tag} (unsaved).")


# ── Cloning grammar editor ─────────────────────────────────────────────────────

class GrammarEditorModal(ModalScreen):
    """View or edit a cloning grammar (overhang table + enzyme + tail).

    Built-in grammars (``gb_l0``, ``moclo_plant``) open here in
    read-only mode — every input is disabled and Save/Delete are
    greyed out so the canonical references can't be corrupted. To
    modify a built-in, use Duplicate in the Parts Bin to fork it
    into an editable custom grammar first.

    Custom grammars open editable. Save validates IUPAC bases on
    every overhang field, requires at least one position, and writes
    the result to ``cloning_grammars.json``. Delete (only enabled
    for custom grammars) drops the entry permanently.

    Dismisses with:
      - ``"saved"`` — grammar persisted (caller should refresh).
      - ``"deleted"`` — grammar removed (caller should flip active
        back to ``gb_l0`` if this was the active one).
      - ``None`` — user cancelled, no changes.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, grammar_id: str) -> None:
        super().__init__()
        self._grammar_id = grammar_id
        self._is_builtin = grammar_id in _BUILTIN_GRAMMARS
        # Snapshot the grammar dict so cancel discards in-flight edits.
        from copy import deepcopy
        grammars = _all_grammars()
        self._grammar = deepcopy(
            grammars.get(grammar_id, _BUILTIN_GRAMMARS["gb_l0"])
        )

    def compose(self) -> ComposeResult:
        g = self._grammar
        editable = not self._is_builtin
        with Vertical(id="ged-dlg"):
            title = f" Grammar: {g.get('name', '?')}"
            if self._is_builtin:
                title += "   [built-in, read-only]"
            yield Static(title, id="ged-title")

            if self._is_builtin:
                yield Static(
                    "  Built-in grammars are read-only. "
                    "Duplicate from the Parts Bin to fork an editable copy.",
                    id="ged-builtin-banner",
                )

            with ScrollableContainer(id="ged-body"):
                yield Label("Name:")
                yield Input(value=g.get("name", ""),
                            id="ged-name", disabled=not editable)

                with Horizontal(id="ged-enzyme-row"):
                    yield Label("Enzyme:", classes="ged-inline-label")
                    yield Input(value=g.get("enzyme", ""),
                                id="ged-enzyme", disabled=not editable)
                    yield Label("Recognition:", classes="ged-inline-label")
                    yield Input(value=g.get("site", ""),
                                id="ged-site", disabled=not editable)

                with Horizontal(id="ged-tail-row"):
                    yield Label("Spacer:", classes="ged-inline-label")
                    yield Input(value=g.get("spacer", ""),
                                id="ged-spacer", disabled=not editable)
                    yield Label("Pad:", classes="ged-inline-label")
                    yield Input(value=g.get("pad", ""),
                                id="ged-pad", disabled=not editable)

                yield Label("Forbidden sites (one per line, NAME=SITE):")
                forbidden = g.get("forbidden_sites", {}) or {}
                if isinstance(forbidden, dict):
                    forb_text = "\n".join(
                        f"{k}={v}" for k, v in forbidden.items()
                    )
                else:
                    forb_text = ""
                yield TextArea(forb_text, id="ged-forbidden",
                               read_only=not editable, soft_wrap=False,
                               show_line_numbers=False)

                yield Label(
                    "Positions (one per line: name | type | 5'OH | 3'OH | color):"
                )
                pos_lines = []
                for p in g.get("positions", []) or []:
                    if not isinstance(p, dict):
                        continue
                    color_field = p.get("color") or ""
                    pos_lines.append(
                        f"{p.get('name','')} | {p.get('type','')} | "
                        f"{p.get('oh5','')} | {p.get('oh3','')} | {color_field}"
                    )
                yield TextArea("\n".join(pos_lines), id="ged-positions",
                               read_only=not editable, soft_wrap=False,
                               show_line_numbers=False)

                yield Label("Coding types (comma-separated, eligible for codon repair):")
                coding = g.get("coding_types", []) or []
                yield Input(value=", ".join(coding),
                            id="ged-coding", disabled=not editable)

            yield Static("", id="ged-status", markup=True)

            with Horizontal(id="ged-btns"):
                yield Button("Save",   id="btn-ged-save",
                             variant="primary", disabled=not editable)
                yield Button("Cancel", id="btn-ged-cancel")
                yield Button("Delete", id="btn-ged-delete",
                             variant="error", disabled=not editable)

    @on(Button.Pressed, "#btn-ged-cancel")
    def _cancel(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)

    @on(Button.Pressed, "#btn-ged-save")
    def _save(self, _) -> None:
        if self._is_builtin:
            return
        parsed = self._gather()
        if parsed is None:
            return
        entries = _load_custom_grammars()
        for i, e in enumerate(entries):
            if e.get("id") == self._grammar_id:
                entries[i] = parsed
                break
        else:
            entries.append(parsed)
        try:
            _save_custom_grammars(entries)
        except (OSError, ValueError) as exc:
            _log.exception("Grammar save failed")
            self.app.notify(f"Save failed: {exc}", severity="error")
            return
        self.app.notify(f"Saved grammar '{parsed.get('name')}'.")
        self.dismiss("saved")

    @on(Button.Pressed, "#btn-ged-delete")
    def _delete(self, _) -> None:
        if self._is_builtin:
            return
        entries = _load_custom_grammars()
        new_entries = [e for e in entries if e.get("id") != self._grammar_id]
        if len(new_entries) == len(entries):
            self.app.notify(
                "Grammar not found in custom grammars file.",
                severity="warning",
            )
            return
        try:
            _save_custom_grammars(new_entries)
        except (OSError, ValueError) as exc:
            _log.exception("Grammar delete failed")
            self.app.notify(f"Delete failed: {exc}", severity="error")
            return
        self.app.notify(
            f"Deleted grammar '{self._grammar.get('name', self._grammar_id)}'."
        )
        self.dismiss("deleted")

    # ── Parsing ──────────────────────────────────────────────────────────────

    def _gather(self) -> "dict | None":
        """Form → grammar dict. Returns None and writes a red status
        line on validation failure; never raises."""
        try:
            name      = self.query_one("#ged-name",      Input).value.strip()
            enzyme    = self.query_one("#ged-enzyme",    Input).value.strip()
            site      = self.query_one("#ged-site",      Input).value.strip().upper()
            spacer    = self.query_one("#ged-spacer",    Input).value.strip().upper()
            pad       = self.query_one("#ged-pad",       Input).value.strip().upper()
            forb_text = self.query_one("#ged-forbidden", TextArea).text
            pos_text  = self.query_one("#ged-positions", TextArea).text
            coding    = self.query_one("#ged-coding",    Input).value
        except NoMatches:
            return None
        status = self.query_one("#ged-status", Static)
        valid_iupac = set("ACGTRYWSMKBDHVN")

        if not name:
            status.update("[red]Name cannot be empty.[/red]")
            return None
        for label, val in (("Enzyme", enzyme),):
            if not val:
                status.update(f"[red]{label} cannot be empty.[/red]")
                return None
        for label, val in (
            ("Recognition site", site),
            ("Spacer", spacer),
            ("Pad",   pad),
        ):
            if not val:
                status.update(f"[red]{label} cannot be empty.[/red]")
                return None
            bad = [c for c in val if c not in valid_iupac]
            if bad:
                status.update(
                    f"[red]{label} contains invalid bases: "
                    f"{''.join(sorted(set(bad)))[:10]}[/red]"
                )
                return None

        # Forbidden sites: NAME=SITE per line.
        forbidden: dict[str, str] = {}
        for raw in forb_text.splitlines():
            line = raw.strip()
            if not line:
                continue
            if "=" not in line:
                status.update(
                    f"[red]Forbidden line must be NAME=SITE: {line!r}[/red]"
                )
                return None
            k, v = line.split("=", 1)
            k, v = k.strip(), v.strip().upper()
            if not k or not v:
                continue
            bad = [c for c in v if c not in valid_iupac]
            if bad:
                status.update(
                    f"[red]Forbidden site {k!r} has invalid bases: "
                    f"{''.join(sorted(set(bad)))[:10]}[/red]"
                )
                return None
            forbidden[k] = v

        # Positions: name | type | oh5 | oh3 [| color]
        positions: list[dict] = []
        for raw in pos_text.splitlines():
            line = raw.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                status.update(
                    f"[red]Position line needs ≥4 fields "
                    f"(name | type | 5'OH | 3'OH): {line!r}[/red]"
                )
                return None
            pname, ptype = parts[0], parts[1]
            poh5, poh3 = parts[2].upper(), parts[3].upper()
            if not pname or not ptype:
                status.update(
                    f"[red]Position name and type cannot be empty: {line!r}[/red]"
                )
                return None
            for label, val in (("5'OH", poh5), ("3'OH", poh3)):
                bad = [c for c in val if c not in valid_iupac]
                if bad:
                    status.update(
                        f"[red]Position {pname!r} {label} has invalid bases: "
                        f"{''.join(sorted(set(bad)))[:10]}[/red]"
                    )
                    return None
            color = parts[4] if len(parts) >= 5 and parts[4] else None
            entry: dict = {
                "name": pname, "type": ptype, "oh5": poh5, "oh3": poh3,
            }
            if color:
                entry["color"] = color
            positions.append(entry)

        if not positions:
            status.update("[red]Need at least one position.[/red]")
            return None

        coding_types = [t.strip() for t in coding.split(",") if t.strip()]

        # Auto-derive type_to_insdc — best-effort heuristic so Save As
        # Feature works without forcing the user to hand-edit JSON.
        # Coding types collapse to "CDS"; obvious labels (Promoter,
        # Terminator, 5'UTR/3'UTR) get their INSDC analogue; everything
        # else falls back to "misc_feature".
        type_to_insdc: dict[str, str] = {}
        for pos in positions:
            ptype = pos["type"]
            lower = ptype.lower()
            if ptype in coding_types:
                type_to_insdc[ptype] = "CDS"
            elif lower == "promoter":
                type_to_insdc[ptype] = "promoter"
            elif "terminator" in lower:
                type_to_insdc[ptype] = "terminator"
            elif lower.startswith("5") and "utr" in lower:
                type_to_insdc[ptype] = "5'UTR"
            elif lower.startswith("3") and "utr" in lower:
                type_to_insdc[ptype] = "3'UTR"
            else:
                type_to_insdc[ptype] = "misc_feature"

        return {
            "id":              self._grammar_id,
            "name":            name,
            "enzyme":          enzyme,
            "site":            site,
            "spacer":          spacer,
            "pad":             pad,
            "forbidden_sites": forbidden,
            "positions":       positions,
            "coding_types":    coding_types,
            "type_to_insdc":   type_to_insdc,
            "catalog":         self._grammar.get("catalog", []) or [],
            "editable":        True,
        }


# ── Parts bin modal ────────────────────────────────────────────────────────────

class PartsBinModal(Screen):
    """Modular cloning parts library — full-screen view.

    Uses Screen (not ModalScreen) so it fills the terminal cleanly instead
    of floating a fixed-size box on a dark overlay. Escape or the Close
    button pops back to the main app.

    The active **grammar** (Golden Braid L0, MoClo Plant, or any
    user-defined custom) gates what's visible: only catalog and
    user-saved parts tagged with the active grammar's id show up.
    The Grammar dropdown at the top of the modal switches grammars
    (persisted in `settings.json` as ``active_grammar``); the
    overhang table beneath it is the canonical Position / Type / 5'OH
    / 3'OH map for the active grammar so the user always knows which
    sticky ends apply. ``Edit`` and ``Duplicate`` open the grammar
    editor; built-in grammars are read-only there, but
    "Duplicate as Custom" forks them into a new editable grammar.

    The "Feat Lib" column flags parts already registered in the
    persistent feature library. We build a single
    ``{(name, feature_type): sequence_upper}`` index on first
    populate and reuse it across renders. The cache is gated on
    ``_features_generation`` — bumped by every ``_save_features``
    call — so the column stays in sync with edits made by Save As
    Feature, Ctrl+Shift+F capture, or the Feature Library workbench
    without paying the scan cost on every populate.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self) -> None:
        super().__init__()
        # Pre-built feature-library index for the "Feat Lib" column.
        # Re-derived only when `_features_generation` differs from the
        # snapshot we last saw, so opening the parts bin from a
        # session where the feature library hasn't changed avoids the
        # whole scan.
        self._feat_lib_index: dict[tuple[str, str], str] = {}
        self._feat_lib_gen_seen: int = -1

    def compose(self) -> ComposeResult:
        """Single-pane loadout: parts table dominates, detail + sequence
        peek out at the bottom, all buttons live on a single bottom row.
        Grammar selection happens inside the New Part modal — every part
        in the table carries its own ``grammar`` id, surfaced as a
        column rather than a top-of-modal filter.
        """
        yield Header()
        with Vertical(id="parts-box"):
            yield Static(" Parts Bin ", id="parts-title")
            yield DataTable(id="parts-table", cursor_type="row", zebra_stripes=True)
            yield Static("", id="parts-detail")
            # Read-only TextArea so the full sequence is visible with a
            # scrollbar and selectable/copyable with the standard terminal
            # gestures (click inside to focus → Ctrl+A → Ctrl+C). Click
            # anywhere in the area auto-selects the whole sequence so
            # single-click → Ctrl+C is enough.
            yield TextArea(
                "", id="parts-seq-view",
                read_only=True, soft_wrap=True, show_line_numbers=False,
            )
            with Horizontal(id="parts-btns"):
                yield Button("Copy Raw",        id="btn-parts-copy-raw")
                yield Button("Copy Primed",     id="btn-parts-copy-primed")
                yield Button("Copy Cloned",     id="btn-parts-copy-cloned")
                yield Button("New Part",        id="btn-new-part",    variant="primary")
                yield Button("Save As Feature", id="btn-parts-save-as-feature")
                yield Button("Export FASTA",    id="btn-parts-export-fasta")
                yield Button("Close",           id="btn-parts-close")
        yield Footer()

    def _active_grammar_id(self) -> str:
        """Currently-active grammar id, defaulting to gb_l0 when the
        setting hasn't been written yet. Used by ``_new_part`` to tag
        freshly-saved parts with whichever grammar the New Part modal
        was set to."""
        return _get_setting("active_grammar", "gb_l0")

    # ── Row data ─────────────────────────────────────────────────────────────

    def _all_rows(self) -> list[dict]:
        """Every part in the bin, regardless of grammar — the table now
        shows them all and tags each with a Grammar column. User parts
        first, then every built-in grammar's catalog appended in
        registry order. Legacy user parts (no ``grammar`` field, from
        pre-grammar versions of SpliceCraft) default to ``gb_l0``."""
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
                "grammar":  p.get("grammar", "gb_l0"),
            })
        # Built-in catalogs from every grammar — concatenated, with each
        # row tagged so the Grammar column makes it clear which assembly
        # standard the row belongs to.
        for gid, grammar in _all_grammars().items():
            for row in grammar.get("catalog", []):
                try:
                    name, ptype, pos, oh5, oh3, backbone, marker = row
                except (TypeError, ValueError):
                    continue
                rows.append({
                    "name": name, "type": ptype, "position": pos,
                    "oh5": oh5, "oh3": oh3, "backbone": backbone,
                    "marker": marker, "sequence": "", "fwd_primer": "",
                    "rev_primer": "", "fwd_tm": 0.0, "rev_tm": 0.0,
                    "user": False,
                    "grammar": gid,
                })
        return rows

    def _type_color_map(self) -> dict[str, str]:
        """Return ``{type_name: rich_color}`` merged across every
        grammar's positions. Used by ``_populate`` to colour the Type
        cell for any row regardless of which grammar it belongs to —
        legacy GB types fall through to ``_GB_TYPE_COLORS``."""
        out = dict(_GB_TYPE_COLORS)
        for grammar in _all_grammars().values():
            for pos in grammar.get("positions", []):
                ptype = pos.get("type")
                color = pos.get("color")
                if isinstance(ptype, str) and ptype and isinstance(color, str) and color:
                    out[ptype] = color
        return out

    def on_mount(self) -> None:
        t = self.query_one("#parts-table", DataTable)
        t.add_columns(
            "Name", "Type", "Pos", "5'OH", "3'OH", "Sequence",
            "Feat Lib", "Grammar",
        )
        self._populate()

    def _populate(self) -> None:
        t = self.query_one("#parts-table", DataTable)
        t.clear()
        self._rows = self._all_rows()
        # Refresh the feature-library index up-front so every row's
        # cell render is an O(1) dict lookup. Without this each row
        # would call `_feature_library_match`, which itself calls
        # `_load_features()` and walks the entire library list — O(N×M)
        # per populate.
        self._refresh_feat_lib_index()
        type_colors = self._type_color_map()
        # Pre-resolve grammar display names so the Grammar column
        # shows human-readable labels instead of raw ids.
        grammars = _all_grammars()
        grammar_label = {
            gid: g.get("name", gid) for gid, g in grammars.items()
        }
        for r in self._rows:
            color = type_colors.get(r["type"], "white")
            seq_preview = r["sequence"][:28] + "…" if len(r["sequence"]) > 28 else r["sequence"]
            if not seq_preview:
                seq_preview = "—"
            gid = r.get("grammar", "gb_l0")
            grammar_cell = Text(
                grammar_label.get(gid, gid),
                style="dim",
            )
            t.add_row(
                Text(r["name"], style=color),
                Text(r["type"], style=f"dim {color}"),
                r["position"],
                Text(r["oh5"], style="bold cyan"),
                Text(r["oh3"], style="bold cyan"),
                Text(seq_preview, style="dim color(252)"),
                self._lib_status_cell(r),
                grammar_cell,
            )

    def _refresh_feat_lib_index(self) -> None:
        """Rebuild the feature-library index iff the global
        ``_features_generation`` advanced since the last build. Reused
        across every ``_populate`` so a no-op repaint (e.g., the
        parts table re-rendering for an unrelated reason) doesn't
        re-scan the whole feature library. The first call after mount
        always rebuilds since `_feat_lib_gen_seen` defaults to -1.
        """
        current = _features_generation
        if current == self._feat_lib_gen_seen:
            return
        self._feat_lib_index = _build_feature_library_index()
        self._feat_lib_gen_seen = current

    def _lib_status_cell(self, row: dict) -> Text:
        """Render the "Feat Lib" column cell for a parts-bin row.

        - User part with exact match in feature library → green ✓.
        - User part with same (name, type) but a different sequence
          (so a Save would replace) → yellow ✓ as a "stale" marker.
        - Anything else (built-in catalog, unsaved part) → empty.

        Built-in rows have no `sequence`, so they always render empty;
        comparing their type/name to library entries would surface
        false positives if a user happened to give a feature library
        entry the same name as a catalog row.

        Lookup goes through the pre-built `_feat_lib_index` (refreshed
        in `_populate` only when `_features_generation` advances), so
        this method is O(1) per row.
        """
        if not row.get("user") or not row.get("sequence"):
            return Text("")
        insdc = _GB_PART_TYPE_TO_INSDC.get(row["type"], "misc_feature")
        match = _classify_feature_library_match(
            self._feat_lib_index,
            row.get("name", ""), insdc, row.get("sequence", ""),
        )
        if match == "exact":
            return Text("✓", style="bold green")
        if match == "name":
            return Text("✓", style="bold yellow")
        return Text("")

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
        detail.append("\n")
        detail.append(f"Position: {r['position']}   ", style="white")
        detail.append("5′ OH: ", style="dim")
        detail.append(r["oh5"], style="bold cyan")
        detail.append("   3′ OH: ", style="dim")
        detail.append(r["oh3"], style="bold cyan")
        detail.append(f"   Backbone: {r['backbone']}   Sel: {r['marker']}", style="dim")
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

        # Full sequence drops into the TextArea below. Built-in catalog
        # parts have no sequence; show a friendly placeholder so the
        # scroll area doesn't look broken.
        seq_view = self.query_one("#parts-seq-view", TextArea)
        if r["sequence"]:
            header = (
                f"> {r['name']} | {r['type']} | pos {r['position']} | "
                f"{r['oh5']}…{r['oh3']} | {len(r['sequence'])} bp\n"
            )
            seq_view.text = header + r["sequence"]
        else:
            seq_view.text = (
                "(Built-in catalog entry — no sequence attached. "
                "Create a new part to see the insert, primed amplicon, "
                "and cloned plasmid here.)"
            )

    # Clicking the sequence area primes it for copy-all: first focus it,
    # then select every character so the user can Ctrl+C immediately.
    @on(TextArea.SelectionChanged, "#parts-seq-view")
    def _seq_selection_changed(self, event) -> None:
        # No-op — Textual wants the handler registered but we drive
        # selection explicitly via _select_sequence below. Kept to suppress
        # "unhandled message" noise on terminals that emit it.
        pass

    def on_click(self, event) -> None:
        """Clicking inside the read-only sequence TextArea selects the
        entire sequence (ready for Ctrl+C). Any other click bubbles up to
        Textual's default handlers so buttons and the table still work."""
        try:
            seq_view = self.query_one("#parts-seq-view", TextArea)
        except NoMatches:
            return
        widget = getattr(event, "widget", None)
        if widget is seq_view or (widget is not None and seq_view in widget.ancestors):
            self._select_sequence_in_view(seq_view)

    def _select_sequence_in_view(self, seq_view: TextArea) -> None:
        """Select the entire TextArea content so Ctrl+C in the terminal
        copies everything. TextArea.select_all is the supported API."""
        try:
            seq_view.focus()
            seq_view.select_all()
        except Exception:
            _log.exception("parts-bin: failed to select sequence text")

    # ── Copy helpers ──────────────────────────────────────────────────────

    def _selected_user_row(self) -> dict | None:
        """Return the currently-highlighted row IF it's a user part with a
        sequence, else notify and return None. Built-in catalog rows
        don't carry sequence/primer data, so every Copy action needs the
        same guard."""
        try:
            t = self.query_one("#parts-table", DataTable)
        except NoMatches:
            return None
        idx = t.cursor_row
        if (idx is None or idx < 0
                or not hasattr(self, "_rows") or idx >= len(self._rows)):
            self.app.notify("Select a part first.", severity="warning")
            return None
        r = self._rows[idx]
        if not r.get("sequence"):
            self.app.notify(
                "Built-in catalog parts have no sequence to copy. "
                "Create a new part first.",
                severity="warning",
            )
            return None
        return r

    def _copy_and_notify(self, text: str, label: str, bp_note: str) -> None:
        """Push `text` to the terminal clipboard via OSC 52 and notify the
        user. Falls back to a notification-only path if the terminal
        refused the escape sequence."""
        ok = _copy_to_clipboard_osc52(text)
        if ok:
            self.app.notify(f"Copied {label} to clipboard ({bp_note}).")
        else:
            self.app.notify(
                f"Could not access clipboard — select the sequence "
                f"panel and press Ctrl+C instead ({label}, {bp_note}).",
                severity="warning",
            )

    @on(Button.Pressed, "#btn-parts-copy-raw")
    def _copy_raw(self, _) -> None:
        """Copy just the insert (raw part sequence, no primer tails)."""
        r = self._selected_user_row()
        if r is None:
            return
        seq = r["sequence"]
        self._copy_and_notify(seq, "raw sequence", f"{len(seq)} bp")

    @on(Button.Pressed, "#btn-parts-copy-primed")
    def _copy_primed(self, _) -> None:
        """Copy the full PCR amplicon (insert + primer tails = pad + Esp3I
        + spacer + oh5 + insert + oh3 + rc(spacer+Esp3I+pad)).

        Older saved parts may predate the simulator, so recompute on the
        fly if `primed_seq` is missing. Parts created before the Esp3I
        switch (v0.3.2) still carry their original BsaI-primed sequence
        in `primed_seq` — the fallback only fires when that field is
        absent, at which point the current L0 enzyme (Esp3I) is used."""
        r = self._selected_user_row()
        if r is None:
            return
        # Use the part's stored grammar so a copy of an old MoClo part
        # uses BsaI tails even after the user has flipped the active
        # grammar to GB; legacy parts (no grammar field) fall back to
        # gb_l0 which preserves v0.3.x behaviour.
        part_grammar = _all_grammars().get(
            r.get("grammar", "gb_l0"), _BUILTIN_GRAMMARS["gb_l0"],
        )
        seq = r.get("primed_seq") or _simulate_primed_amplicon(
            r["sequence"], r.get("oh5", ""), r.get("oh3", ""),
            grammar=part_grammar,
        )
        self._copy_and_notify(seq, "primed amplicon", f"{len(seq)} bp")

    @on(Button.Pressed, "#btn-parts-copy-cloned")
    def _copy_cloned(self, _) -> None:
        """Copy the simulated cloned plasmid (insert ligated into pUPD2
        backbone stub, linearised at the 5' overhang)."""
        r = self._selected_user_row()
        if r is None:
            return
        seq = r.get("cloned_seq") or _simulate_cloned_plasmid(
            r["sequence"], r.get("oh5", ""), r.get("oh3", ""),
        )
        self._copy_and_notify(
            seq, "cloned plasmid", f"{len(seq)} bp circular (linearised at 5′ OH)",
        )

    @on(Button.Pressed, "#btn-new-part")
    def _new_part(self, _) -> None:
        # Opens the domesticator modal. The current record's sequence + name
        # are passed so the "Feature from plasmid" source can default to the
        # plasmid the user already has open.
        rec = getattr(self.app, "_current_record", None)
        seq = str(rec.seq) if rec else ""
        feats = []
        try:
            pm = self.app.query_one("#plasmid-map", PlasmidMap)
            feats = pm._feats
        except NoMatches:
            pass
        current_name = (getattr(rec, "name", "") or getattr(rec, "id", "") or "") if rec else ""

        def _on_result(part_dict):
            if part_dict is None:
                return
            # Tag the part with the active grammar so the parts bin
            # filter can route it to the right tab next time. Legacy
            # parts (saved before grammars existed) default to gb_l0
            # in `_all_rows`.
            part_dict.setdefault("grammar", self._active_grammar_id())
            entries = _load_parts_bin()
            entries.insert(0, part_dict)
            _save_parts_bin(entries)
            self._populate()
            self.app.notify(
                f"Saved '{part_dict['name']}' to Parts Bin "
                f"({len(part_dict.get('sequence', ''))} bp).",
            )

        self.app.push_screen(
            DomesticatorModal(seq, feats, current_plasmid_name=current_name),
            callback=_on_result,
        )

    @on(Button.Pressed, "#btn-parts-export-fasta")
    def _export_fasta(self, _) -> None:
        """Export the highlighted row's sequence as single-record FASTA.

        Built-in catalog parts have no sequence attached, so we bail
        with a friendly notify rather than pushing an empty modal."""
        try:
            t = self.query_one("#parts-table", DataTable)
        except NoMatches:
            return
        idx = t.cursor_row
        if (idx is None or idx < 0
                or not hasattr(self, "_rows") or idx >= len(self._rows)):
            self.app.notify("Select a part first.", severity="warning")
            return
        r = self._rows[idx]
        seq = r.get("sequence", "") or ""
        if not seq:
            self.app.notify(
                "Built-in catalog parts have no sequence to export. "
                "Create a new part first.",
                severity="warning",
            )
            return
        name = r.get("name", "part")
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", name) or "part"
        default_path = str(Path.home() / f"{safe}.fa")

        def _on_done(summary):
            if not summary:
                return
            self.app.notify(
                f"Exported '{summary['name']}' to {summary['path']} "
                f"({summary['bp']} bp).",
            )

        self.app.push_screen(
            FastaExportModal(
                name=name,
                sequence=seq,
                default_path=default_path,
                subtitle=f"[{name}]  [{r.get('type', '?')}]  {len(seq)} bp",
            ),
            callback=_on_done,
        )

    @on(Button.Pressed, "#btn-parts-save-as-feature")
    def _save_as_feature(self, _) -> None:
        """Register the highlighted user part as a feature-library entry.

        Built-in catalog rows have no sequence — they're descriptive
        position/overhang bookkeeping with no insert — so they're
        rejected. For real user parts, build a prefill dict mapping
        the GB part shape onto INSDC feature vocabulary
        (`_GB_PART_TYPE_TO_INSDC`) and open `AddFeatureModal` so the
        user can adjust the name / qualifiers before committing.
        Saving routes through ``app._persist_feature_entry`` (the same
        helper Ctrl+Shift+F capture uses) so the latest write wins on
        (name, feature_type) collisions.

        Before opening the modal, ``_feature_library_match`` is
        consulted: if the (name, INSDC type) pair is already in the
        library the user gets a yellow warning notify so a silent
        replace doesn't surprise them. The modal opens in either case
        — the user retains the option to cancel.

        Insert-at-cursor is disabled here even when a record is open —
        the parts bin is a different mental model (registering the
        part's *idea* into the library), and inserting from a modal
        layered on top of another modal would leave the parts bin
        stranded over the freshly-edited record.
        """
        r = self._selected_user_row()
        if r is None:
            return

        ftype = _GB_PART_TYPE_TO_INSDC.get(r["type"], "misc_feature")
        sequence = r.get("sequence", "")
        oh5, oh3 = r.get("oh5", ""), r.get("oh3", "")
        pos = r.get("position", "")
        backbone = r.get("backbone", "")
        bits: list[str] = ["Golden Braid L0 part"]
        if pos:
            bits.append(f"Position {pos}")
        if oh5 or oh3:
            bits.append(f"5' OH {oh5} / 3' OH {oh3}")
        if backbone:
            bits.append(f"backbone {backbone}")
        # CDS-NS / C-tag both collapse to plain "CDS" in INSDC. The
        # NS / C-tag distinction is meaningful (NS = no stop codon,
        # C-tag = C-terminal fusion fragment), so preserve it in the
        # description rather than losing it on the round-trip.
        if r["type"] in {"CDS-NS", "C-tag"}:
            bits.append(f"GB type: {r['type']}")
        description = "; ".join(bits)

        # Warn before opening so the user can spot a collision they
        # didn't intend. `_feature_library_match` distinguishes exact
        # (no-op save) from name-only (Save will replace the stored
        # sequence) — the wording differs because the consequences do.
        name = r.get("name", "")
        match = _feature_library_match(name, ftype, sequence)
        if match == "exact":
            self.app.notify(
                f"'{name}' is already in the feature library "
                f"(exact match). Saving again is a no-op.",
                severity="warning",
            )
        elif match == "name":
            self.app.notify(
                f"A feature named '{name}' (type {ftype}) is already "
                f"in the feature library with a different sequence. "
                f"Saving will replace it.",
                severity="warning",
            )

        prefill = {
            "name":         name,
            "feature_type": ftype,
            "sequence":     sequence,
            "strand":       1,
            "color":        None,
            "qualifiers":   {},
            "description":  description,
        }

        def _on_done(result):
            if not result:
                return
            entry = result.get("entry") if isinstance(result, dict) else None
            if not entry or result.get("action") != "save":
                return
            persist = getattr(self.app, "_persist_feature_entry", None)
            if persist is None or not persist(entry):
                return
            self.app.notify(
                f"Saved '{entry.get('name')}' as a "
                f"{entry.get('feature_type')} feature."
            )
            # Refresh so the "Feat Lib" column flips to ✓ for the row
            # we just registered.
            self._populate()

        self.app.push_screen(
            AddFeatureModal(prefill=prefill, have_cursor=False),
            callback=_on_done,
        )

    @on(Button.Pressed, "#btn-parts-close")
    def _close(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── FASTA file picker ──────────────────────────────────────────────────────────
#
# Used by DomesticatorModal's "Open FASTA" source. The picker is a standard
# DirectoryTree wrapped in a modal, but with custom label rendering that
# paints FASTA files lime green and everything else white so the user can
# scan a mixed directory quickly. Returns an absolute path string on
# dismiss, or None on cancel.

# Lowercased extensions (including the leading dot) that count as FASTA
# for the picker's highlight rule. Kept broad so GenBank-style FASTA
# dumps (.fna for nucleotide, .faa for amino acid, .ffn for coding,
# etc.) all light up.
_FASTA_EXTS: frozenset[str] = frozenset({
    ".fa", ".fasta", ".fna", ".ffn", ".frn", ".fas", ".mpfa", ".faa",
})

# Lime green and plain white — deliberately high contrast so the eye
# finds FASTA files immediately in a mixed directory listing.
_FASTA_PICKER_FASTA_STYLE = "bold #BFFF00"
_FASTA_PICKER_OTHER_STYLE = "#FFFFFF"


def _is_fasta_path(path) -> bool:
    """True if ``path`` looks like a FASTA file by extension. Accepts
    anything with a ``suffix`` attribute (``pathlib.Path`` or ``DirEntry``)
    or a plain string."""
    try:
        suffix = getattr(path, "suffix", None)
        if suffix is None:
            suffix = Path(str(path)).suffix
    except Exception:
        return False
    return suffix.lower() in _FASTA_EXTS


def _parse_fasta_single(path: str) -> tuple[str, str]:
    """Parse a FASTA file that must contain **exactly one** record and
    return ``(record_id, sequence)``.

    Multi-record FASTA files are rejected: the domesticator / parts bin
    flow only makes sense for a single part, so we surface a helpful
    error rather than silently picking the first record.

    Raises ``ValueError`` with a user-friendly message on any failure
    (read errors, zero records, multiple records, empty or non-IUPAC
    sequence). The sequence is upper-cased on success and validated
    against the IUPAC alphabet plus ``-``/``*``/``X`` for gap / stop /
    unknown."""
    from Bio import SeqIO
    try:
        records = list(SeqIO.parse(path, "fasta"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Failed to read FASTA: {exc}") from exc
    if not records:
        raise ValueError("No FASTA records found in file.")
    if len(records) > 1:
        raise ValueError(
            f"Multi-sequence FASTA not supported ({len(records)} records "
            "found). Please provide a single-record FASTA."
        )
    rec = records[0]
    seq = str(rec.seq).upper()
    if not seq:
        raise ValueError("FASTA record has empty sequence.")
    valid = set("ACGTURYMKSWBDHVN-X*")
    bad = sorted(set(seq) - valid)
    if bad:
        raise ValueError(
            f"Non-IUPAC characters in sequence: {''.join(bad[:8])}"
        )
    return (rec.id or "fasta", seq)


class _FastaAwareDirectoryTree(DirectoryTree):
    """DirectoryTree variant that colours FASTA files lime green and
    every other file white. Directories are left alone so Textual's
    default folder styling still applies."""

    def render_label(self, node, base_style, style):
        label = super().render_label(node, base_style, style)
        data = node.data
        if data is None:
            return label
        p = getattr(data, "path", None)
        if p is None:
            return label
        try:
            if not p.is_file():
                return label
        except OSError:
            return label
        styled = label.copy()
        if _is_fasta_path(p):
            styled.stylize(_FASTA_PICKER_FASTA_STYLE)
        else:
            styled.stylize(_FASTA_PICKER_OTHER_STYLE)
        return styled


class FastaFilePickerModal(ModalScreen):
    """Modal file browser that returns the path to a selected FASTA file.

    Dismisses with ``str`` (absolute path) on Open, or ``None`` on Cancel /
    Escape. FASTA files are painted lime green in the tree; other files
    are white so the user can scan a mixed directory quickly. The tree
    starts in ``start_path`` when given (and readable), else ``$HOME``."""

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, start_path: "str | None" = None) -> None:
        super().__init__()
        start = Path(start_path).expanduser() if start_path else Path.home()
        try:
            if not start.is_dir():
                start = Path.home()
        except OSError:
            start = Path.home()
        self._start = str(start)
        self._selected: "str | None" = None

    def compose(self) -> ComposeResult:
        with Vertical(id="fasta-box"):
            yield Static(" Open FASTA File ", id="fasta-title")
            yield Static(
                f"[dim]{self._start}[/dim]", id="fasta-header", markup=True
            )
            yield _FastaAwareDirectoryTree(self._start, id="fasta-tree")
            yield Static(
                "[dim]FASTA files are highlighted in lime green. "
                "Click a file, then Open.[/dim]",
                id="fasta-hint", markup=True,
            )
            yield Static("", id="fasta-status", markup=True)
            with Horizontal(id="fasta-btns"):
                yield Button("Open", id="btn-fasta-open",
                             variant="primary", disabled=True)
                yield Button("Cancel", id="btn-fasta-cancel")

    def on_mount(self) -> None:
        try:
            self.query_one("#fasta-tree", _FastaAwareDirectoryTree).focus()
        except NoMatches:
            pass

    @on(DirectoryTree.FileSelected)
    def _on_file_selected(self, event) -> None:
        self._selected = str(event.path)
        try:
            self.query_one("#fasta-header", Static).update(
                f"[dim]{self._selected}[/dim]"
            )
            self.query_one("#btn-fasta-open", Button).disabled = False
            self.query_one("#fasta-status", Static).update("")
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-fasta-open")
    def _open(self) -> None:
        if self._selected:
            self.dismiss(self._selected)
            return
        try:
            self.query_one("#fasta-status", Static).update(
                "[red]Pick a file first.[/red]"
            )
        except NoMatches:
            pass

    @on(Button.Pressed, "#btn-fasta-cancel")
    def _cancel_btn(self) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Constructor modal ──────────────────────────────────────────────────────────

def _feats_for_domesticator(record) -> list[dict]:
    """Parse a SeqRecord into the shape the DomesticatorModal feature picker
    expects: ``{label, type, start, end, strand}`` per non-source feature.

    Kept deliberately simpler than ``PlasmidMap._parse`` — compound/wrap
    features are flattened to their outer bounds (the domesticator only
    needs slice coordinates, not rendering geometry), and restriction-site
    overlays are skipped. Features with non-integer coords or zero width
    are dropped rather than raising."""
    out: list[dict] = []
    total = len(getattr(record, "seq", "") or "")
    for feat in getattr(record, "features", []) or []:
        if feat.type in ("source", "resite", "recut"):
            continue
        try:
            start = int(feat.location.start)
            end   = int(feat.location.end)
        except (TypeError, ValueError):
            continue
        if total > 0:
            start = max(0, min(start, total))
            end   = max(0, min(end,   total))
        if start == end:
            continue
        out.append({
            "label":  _feat_label(feat),
            "type":   feat.type,
            "start":  start,
            "end":    end,
            "strand": getattr(feat.location, "strand", 1) or 1,
        })
    return out


class DomesticatorModal(ModalScreen):
    """Golden Braid L0 Parts Domesticator.

    Takes a template sequence + region, designs domestication primers with
    the correct Esp3I / BsmBI sites + positional overhangs, and returns a
    part dict ready for saving to the Parts Bin.

    Primer structure (5'→3'):
        Forward: GCGC CGTCTC A [5' overhang] [binding region →]
        Reverse: GCGC CGTCTC A [RC 3' OH]    [← binding region RC]

    After Esp3I digestion the amplicon carries the correct 4-nt sticky
    ends for Golden Braid L0 assembly. L0 uses Esp3I (CGTCTC) so the
    domesticated part survives the downstream L1+ BsaI (GGTCTC) assembly
    without re-cutting.
    """

    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, template_seq: str, feats: list[dict],
                 current_plasmid_name: str = ""):
        super().__init__()
        self._template = template_seq.upper()
        self._feats    = feats   # from PlasmidMap._feats (the *current* plasmid)
        self._design:  "dict | None" = None   # result of _design_gb_primers
        # ── Source-picker state ────────────────────────────────────────────
        # Four sources for the part's DNA:
        #   "direct"  : user types/pastes into a TextArea
        #   "featlib" : pick from persistent features.json (feature library)
        #   "plasmid" : pick a plasmid, then a feature from it
        #   "fasta"   : browse the filesystem for a FASTA file, parse first record
        # The "plasmid" source defaults to the plasmid the user has open —
        # swapped via the picker button if they want a library entry instead.
        self._source: str = "direct"
        # Plasmid-source state. `_plasmid_pick_id = None` means "use the
        # current in-app plasmid (seq+feats passed to __init__)".
        self._plasmid_pick_id:    "str | None"  = None
        self._plasmid_pick_name:  str           = (
            current_plasmid_name or ("— current plasmid —" if self._template else "")
        )
        self._plasmid_pick_seq:   str           = self._template
        self._plasmid_pick_feats: list[dict]    = [
            f for f in self._feats if f.get("type") not in ("resite", "recut")
        ]
        # FASTA-source state. Populated by `_on_fasta_picked` after the user
        # picks a file in `FastaFilePickerModal`.
        self._fasta_path: "str | None" = None
        self._fasta_name: str           = ""
        self._fasta_seq:  str           = ""
        # Codon table for silent-mutation repair of internal BsaI / Esp3I
        # sites in coding parts. Seeded to E. coli K12 in on_mount; the user
        # can swap via the "Change…" button.
        self._codon_entry: "dict | None" = None

    # ── Option builders (used by compose + on re-pick) ─────────────────────

    def _featlib_options(self) -> list[tuple[str, str]]:
        """Build Select options for the Feature Library source. Value is the
        integer index into ``_load_features()`` (as a str, since Select values
        must be strings) so the design step can look the entry back up."""
        opts: list[tuple[str, str]] = []
        for i, e in enumerate(_load_features()):
            name = e.get("name", "?") or "?"
            ft   = e.get("feature_type", "misc")
            blen = len(e.get("sequence", "") or "")
            opts.append((f"{name}  [{ft}, {blen} bp]", str(i)))
        return opts

    def _type_options_for(self, grammar: dict) -> tuple[list[tuple[str, str]], str]:
        """Build ``(type_options, default_value)`` for the Part Type
        Select based on a given grammar. Each option is
        ``"{type}  ({position}: 5'OH→3'OH)"`` so the user picks a
        type and immediately sees what overhangs they'll get. Default
        prefers the first coding type (codon-fix repair only works on
        coding parts), falling back to the first listed type."""
        type_options: list[tuple[str, str]] = []
        for pos in grammar.get("positions", []):
            ptype = pos.get("type")
            if not isinstance(ptype, str) or not ptype:
                continue
            label = (
                f"{ptype}  ({pos.get('name','?')}: "
                f"{pos.get('oh5','')}→{pos.get('oh3','')})"
            )
            type_options.append((label, ptype))
        if not type_options:
            # Pathological grammar with no positions — surface
            # something rather than crashing the Select widget.
            type_options = [("(no positions defined)", "")]
        coding_types = grammar.get("coding_types", []) or []
        default_type = next(
            (t for _label, t in type_options if t in coding_types),
            type_options[0][1] if type_options else "",
        )
        return type_options, default_type

    def _plasmid_feat_options(self) -> list[tuple[str, str]]:
        """Build Select options for the currently-picked plasmid's features.
        Value is the index (as str) into ``self._plasmid_pick_feats``."""
        opts: list[tuple[str, str]] = []
        for i, f in enumerate(self._plasmid_pick_feats):
            label = f.get("label") or f.get("type", "?")
            s, e  = f.get("start", 0), f.get("end", 0)
            opts.append((f"{label}  ({s+1}‥{e})", str(i)))
        return opts

    def compose(self) -> ComposeResult:
        # Grammar selection lives inside this modal now (formerly on
        # the Parts Bin). Picking a different grammar from the
        # dropdown re-persists the active-grammar setting and rebuilds
        # the Type Select with that grammar's positions.
        active_grammar = _get_active_grammar()
        active_gid = active_grammar.get("id", "gb_l0")
        type_options, default_type = self._type_options_for(active_grammar)

        with Vertical(id="dom-box"):
            yield Static(
                f" Domesticate Part  —  {active_grammar.get('name', '?')} ",
                id="dom-title",
            )
            # Scrollable body — everything between title and buttons. Primer
            # design results expand vertically, so the body needs to scroll
            # on narrow terminals rather than overflow off-screen.
            with ScrollableContainer(id="dom-body"):
                # ── Row 0: Cloning grammar picker ──
                with Horizontal(id="dom-grammar-row"):
                    yield Label("Cloning grammar")
                    yield Select(
                        _grammar_dropdown_options(),
                        value=active_gid,
                        id="dom-grammar-select",
                        allow_blank=False,
                    )
                # ── Row 1: Part name + type ──
                with Horizontal(id="dom-row1"):
                    with Vertical(id="dom-name-col"):
                        yield Label("Part name")
                        yield Input(placeholder="e.g. my-promoter", id="dom-name")
                    with Vertical(id="dom-type-col"):
                        yield Label("Part type")
                        yield Select(type_options, id="dom-type",
                                     value=default_type)
                # ── Row 2: overhang info (auto-updated from type) ──
                yield Static("", id="dom-oh-info", markup=True)
                # ── Codon table picker (for silent-mutation repair) ──
                with Horizontal(id="dom-codon-row"):
                    yield Static(
                        "Codon table: [bold]E. coli K12[/bold] (taxid 83333)",
                        id="dom-codon-label", markup=True,
                    )
                    yield Button("Change", id="btn-dom-codon",
                                 variant="default")
                # ── Row 3: Source picker ──
                yield Label("Source")
                with RadioSet(id="dom-src"):
                    yield RadioButton("Direct input",         id="dom-src-direct",  value=True)
                    yield RadioButton("Feature library",      id="dom-src-featlib")
                    yield RadioButton("Feature from plasmid", id="dom-src-plasmid")
                    yield RadioButton("Open FASTA",           id="dom-src-fasta")
                # ── Direct-input panel ──
                with Vertical(id="dom-panel-direct", classes="dom-src-panel"):
                    yield Label("Paste or type the part sequence (5'→3'):")
                    yield TextArea("", id="dom-direct-seq")
                # ── Feature-library panel ──
                with Vertical(id="dom-panel-featlib", classes="dom-src-panel"):
                    yield Label("Pick a feature from the feature library:")
                    yield Select(
                        self._featlib_options(),
                        id="dom-featlib-select",
                        prompt="(select feature)",
                    )
                    yield Static("", id="dom-featlib-preview", markup=True)
                # ── Feature-from-plasmid panel ──
                with Vertical(id="dom-panel-plasmid", classes="dom-src-panel"):
                    with Horizontal(id="dom-plasmid-hdr"):
                        yield Label("Plasmid:")
                        yield Static(
                            self._plasmid_pick_name or "(none loaded)",
                            id="dom-plasmid-name",
                        )
                        yield Button("Change", id="btn-dom-pick-plasmid")
                    yield Label("Pick feature:")
                    yield Select(
                        self._plasmid_feat_options(),
                        id="dom-plasmid-feat-select",
                        prompt="(select feature)",
                    )
                    yield Static("", id="dom-plasmid-feat-preview", markup=True)
                # ── Open-FASTA panel ──
                with Vertical(id="dom-panel-fasta", classes="dom-src-panel"):
                    with Horizontal(id="dom-fasta-hdr"):
                        yield Label("File:")
                        yield Static("(no file selected)", id="dom-fasta-name")
                        yield Button("Browse", id="btn-dom-pick-fasta")
                    yield Static("", id="dom-fasta-preview", markup=True)
                # ── Primer results ──
                yield Static("", id="dom-primer-results", markup=True)
            # ── Buttons (pinned below scroll body so they're always reachable) ──
            with Horizontal(id="dom-btns"):
                yield Button(
                    "Design Primers", id="btn-dom-design", variant="primary",
                )
                yield Button(
                    "Save to Parts Bin", id="btn-dom-save", variant="primary",
                    disabled=True,
                )
                yield Button(
                    "Save Primers", id="btn-dom-save-primers",
                    variant="primary", disabled=True,
                )
                yield Button("Cancel", id="btn-dom-cancel")

    def on_mount(self) -> None:
        self._update_oh_display()
        self._refresh_source_panels()
        # Seed codon table registry with the built-in E. coli K12 entry
        # (shared with Mutagenize — registry caches across modals).
        try:
            _codon_tables_load()
            self._codon_entry = _codon_tables_get("83333")
        except Exception:
            _log.exception("Domesticator: codon registry load failed")
            self._codon_entry = None
        self._update_codon_label()
        # Focus the name input
        self.query_one("#dom-name", Input).focus()

    def _update_codon_label(self) -> None:
        try:
            lbl = self.query_one("#dom-codon-label", Static)
        except NoMatches:
            return
        entry = self._codon_entry
        if not entry:
            lbl.update("[red]Codon table: none selected[/red]")
            return
        tax = f" (taxid {entry['taxid']})" if entry.get("taxid") else ""
        lbl.update(f"Codon table: [bold]{entry['name']}[/bold]{tax}")

    @on(Button.Pressed, "#btn-dom-codon")
    def _pick_codon_table(self, _) -> None:
        self.app.push_screen(SpeciesPickerModal(),
                             callback=self._codon_picked)

    def _codon_picked(self, entry: "dict | None") -> None:
        if not entry:
            return
        self._codon_entry = entry
        self._update_codon_label()

    # ── Source-panel visibility ────────────────────────────────────────────

    def _refresh_source_panels(self) -> None:
        """Show the panel matching ``self._source``; hide the others."""
        panels = {
            "direct":  "#dom-panel-direct",
            "featlib": "#dom-panel-featlib",
            "plasmid": "#dom-panel-plasmid",
            "fasta":   "#dom-panel-fasta",
        }
        for key, sel in panels.items():
            try:
                self.query_one(sel).display = (key == self._source)
            except NoMatches:
                pass

    @on(RadioSet.Changed, "#dom-src")
    def _source_changed(self, event: RadioSet.Changed) -> None:
        rb_id = getattr(event.pressed, "id", "") or ""
        mapping = {
            "dom-src-direct":  "direct",
            "dom-src-featlib": "featlib",
            "dom-src-plasmid": "plasmid",
            "dom-src-fasta":   "fasta",
        }
        self._source = mapping.get(rb_id, "direct")
        self._refresh_source_panels()
        # Wipe stale primer results when the user changes source — the
        # previously-designed primers no longer reflect what's on screen.
        try:
            self.query_one("#dom-primer-results", Static).update("")
            self.query_one("#btn-dom-save", Button).disabled = True
            self.query_one("#btn-dom-save-primers", Button).disabled = True
        except NoMatches:
            pass

    # ── Feature-library preview ────────────────────────────────────────────

    @on(Select.Changed, "#dom-featlib-select")
    def _featlib_preview(self, event: Select.Changed) -> None:
        val = event.value
        preview = self.query_one("#dom-featlib-preview", Static)
        if not isinstance(val, str) or not val.isdigit():
            preview.update("")
            return
        entries = _load_features()
        idx = int(val)
        if idx < 0 or idx >= len(entries):
            preview.update("")
            return
        e = entries[idx]
        seq = (e.get("sequence") or "").upper()
        strand = e.get("strand", 1)
        head = seq[:40] + ("…" if len(seq) > 40 else "")
        preview.update(
            f"  [dim]{len(seq)} bp · strand {strand}[/dim]   {head}"
        )

    # ── Plasmid-source handlers ────────────────────────────────────────────

    @on(Button.Pressed, "#btn-dom-pick-plasmid")
    def _pick_plasmid(self, _) -> None:
        """Open the plasmid picker; on selection, load that plasmid's record
        from the library, rebuild the feature dropdown, and refresh the label."""
        def _on_picked(pid):
            if not pid:
                return
            entries = _load_library()
            match = next((e for e in entries if e.get("id") == pid), None)
            if match is None:
                self.app.notify(f"Plasmid '{pid}' not in library.",
                                severity="warning")
                return
            try:
                rec = _gb_text_to_record(match.get("gb_text", "") or "")
            except Exception as exc:
                _log.exception("Failed to parse library entry %s", pid)
                self.app.notify(f"Failed to load '{pid}': {exc}",
                                severity="error")
                return
            self._plasmid_pick_id    = pid
            self._plasmid_pick_name  = match.get("name") or pid
            self._plasmid_pick_seq   = str(getattr(rec, "seq", "") or "").upper()
            self._plasmid_pick_feats = _feats_for_domesticator(rec)
            try:
                self.query_one("#dom-plasmid-name", Static).update(
                    self._plasmid_pick_name
                )
                sel = self.query_one("#dom-plasmid-feat-select", Select)
                sel.set_options(self._plasmid_feat_options())
                self.query_one("#dom-plasmid-feat-preview", Static).update("")
            except NoMatches:
                pass

        self.app.push_screen(
            PlasmidPickerModal(current_id=self._plasmid_pick_id),
            callback=_on_picked,
        )

    # ── FASTA-source handlers ──────────────────────────────────────────────

    @on(Button.Pressed, "#btn-dom-pick-fasta")
    def _pick_fasta(self, _) -> None:
        """Open the FASTA file picker; on selection, parse the file and
        display a preview. The picker starts in ``$HOME`` when first opened,
        then remembers the last picked directory via `self._fasta_path`."""
        def _on_picked(path: "str | None") -> None:
            if not path:
                return
            try:
                name, seq = _parse_fasta_single(path)
            except ValueError as exc:
                self.app.notify(str(exc), severity="error")
                return
            self._fasta_path = path
            self._fasta_name = name
            self._fasta_seq  = seq
            try:
                display = Path(path).name or path
                self.query_one("#dom-fasta-name", Static).update(display)
                head = seq[:40] + ("…" if len(seq) > 40 else "")
                self.query_one("#dom-fasta-preview", Static).update(
                    f"  [dim]{name} · {len(seq)} bp[/dim]   {head}"
                )
                self.query_one("#dom-primer-results", Static).update("")
                self.query_one("#btn-dom-save", Button).disabled = True
            except NoMatches:
                pass

        start_dir = (
            str(Path(self._fasta_path).parent)
            if self._fasta_path else None
        )
        self.app.push_screen(
            FastaFilePickerModal(start_path=start_dir),
            callback=_on_picked,
        )

    @on(Select.Changed, "#dom-plasmid-feat-select")
    def _plasmid_feat_preview(self, event: Select.Changed) -> None:
        val = event.value
        preview = self.query_one("#dom-plasmid-feat-preview", Static)
        if not isinstance(val, str) or not val.isdigit():
            preview.update("")
            return
        idx = int(val)
        if idx < 0 or idx >= len(self._plasmid_pick_feats):
            preview.update("")
            return
        f = self._plasmid_pick_feats[idx]
        s, e = f.get("start", 0), f.get("end", 0)
        total = len(self._plasmid_pick_seq)
        blen  = _feat_len(s, e, total) if total else 0
        preview.update(
            f"  [dim]{f.get('type','?')} · {s+1}‥{e} · {blen} bp · "
            f"strand {f.get('strand', 1)}[/dim]"
        )

    # ── Grammar / part-type changes update the overhang info ──────────────

    @on(Select.Changed, "#dom-grammar-select")
    def _grammar_changed(self, event: Select.Changed) -> None:
        """Switching the cloning grammar from this modal persists the
        choice (so the next New Part / Save As Feature flow defaults
        to it) and rebuilds the Type Select with that grammar's
        positions. The title bar and overhang info reflect the new
        grammar immediately."""
        if event.value == Select.BLANK:
            return
        new_gid = str(event.value)
        if new_gid == _get_setting("active_grammar", "gb_l0"):
            return
        _set_setting("active_grammar", new_gid)
        new_grammar = _get_active_grammar()
        type_options, default_type = self._type_options_for(new_grammar)
        try:
            type_sel = self.query_one("#dom-type", Select)
            type_sel.set_options(type_options)
            type_sel.value = default_type
        except NoMatches:
            pass
        try:
            title = self.query_one("#dom-title", Static)
            title.update(
                f" Domesticate Part  —  {new_grammar.get('name', '?')} "
            )
        except NoMatches:
            pass
        self._update_oh_display()

    @on(Select.Changed, "#dom-type")
    def _type_changed(self, _event) -> None:
        self._update_oh_display()

    def _update_oh_display(self) -> None:
        sel = self.query_one("#dom-type", Select)
        val = sel.value
        grammar = _get_active_grammar()
        pos_spec = (
            _grammar_position_by_type(grammar, val)
            if isinstance(val, str) else None
        )
        if pos_spec is None:
            self.query_one("#dom-oh-info", Static).update("")
            return
        self.query_one("#dom-oh-info", Static).update(
            f"  [dim]{pos_spec.get('name','?')}[/dim]   "
            f"5′ overhang: [bold cyan]{pos_spec.get('oh5','')}[/bold cyan]   →   "
            f"3′ overhang: [bold cyan]{pos_spec.get('oh3','')}[/bold cyan]   "
            f"[dim]({grammar.get('enzyme','?')} domestication)[/dim]"
        )

    # ── Design primers ─────────────────────────────────────────────────────

    def _resolve_source(self) -> "tuple[str, int, int] | str":
        """Return ``(template, start, end)`` for the active source, or a
        short error string describing why it can't be resolved."""
        if self._source == "direct":
            raw = self.query_one("#dom-direct-seq", TextArea).text
            # Strip whitespace and common paste artefacts (line breaks, FASTA
            # headers, numbers). Anything non-ACGTU-or-IUPAC is dropped with
            # a warning so a raw paste "works".
            cleaned = "".join(
                c for c in raw.upper()
                if c in "ACGTURYSWKMBDHVN"
            )
            if not cleaned:
                return "Paste a sequence first."
            return cleaned, 0, len(cleaned)
        if self._source == "featlib":
            val = self.query_one("#dom-featlib-select", Select).value
            if not isinstance(val, str) or not val.isdigit():
                return "Pick a feature from the library."
            entries = _load_features()
            idx = int(val)
            if idx < 0 or idx >= len(entries):
                return "Selected feature is no longer in the library."
            seq = (entries[idx].get("sequence") or "").upper()
            if not seq:
                return "Selected feature has no sequence."
            return seq, 0, len(seq)
        if self._source == "plasmid":
            if not self._plasmid_pick_seq:
                return "Pick a plasmid from the library first."
            val = self.query_one("#dom-plasmid-feat-select", Select).value
            if not isinstance(val, str) or not val.isdigit():
                return "Pick a feature from the plasmid."
            idx = int(val)
            if idx < 0 or idx >= len(self._plasmid_pick_feats):
                return "Selected feature is no longer available."
            f = self._plasmid_pick_feats[idx]
            return self._plasmid_pick_seq, f.get("start", 0), f.get("end", 0)
        if self._source == "fasta":
            if not self._fasta_seq:
                return "Pick a FASTA file first."
            return self._fasta_seq, 0, len(self._fasta_seq)
        return "Unknown source."

    @on(Button.Pressed, "#btn-dom-design")
    def _design(self, _) -> None:
        status = self.query_one("#dom-primer-results", Static)
        part_type = self.query_one("#dom-type", Select).value
        grammar = _get_active_grammar()
        if (not isinstance(part_type, str)
                or _grammar_position_by_type(grammar, part_type) is None):
            status.update("[red]Select a part type.[/red]")
            return
        resolved = self._resolve_source()
        if isinstance(resolved, str):
            status.update(f"[red]{resolved}[/red]")
            return
        template, start, end = resolved
        total = len(template)
        if start < 0 or end < 0 or start > total or end > total:
            status.update(
                f"[red]Invalid region: {start+1}–{end} "
                f"(sequence is {total} bp)[/red]"
            )
            return
        if start == end:
            status.update("[red]Region is empty.[/red]")
            return
        region_len = _feat_len(start, end, total)
        if region_len < 20:
            status.update(f"[red]Region too short ({region_len} bp, need ≥ 20).[/red]")
            return

        codon_raw = (self._codon_entry or {}).get("raw")
        try:
            self._design = _design_gb_primers(
                template, start, end, part_type, codon_raw=codon_raw,
                grammar=grammar,
            )
        except Exception as exc:
            _log.exception("Primer design failed")
            status.update(f"[red]Primer design failed: {exc}[/red]")
            return

        if "error" in self._design:
            msg = self._design["error"]
            muts = self._design.get("mutations") or []
            body = f"[red]{msg}[/red]"
            if muts:
                body += "\n\n[dim]Silent mutations that were applied before "
                body += "giving up:[/dim]\n"
                for m in muts:
                    body += f"  · {m}\n"
            status.update(body)
            return

        d = self._design
        pairs = d.get("pairs") or []
        t = Text()
        t.append("── Primers designed ─────────────────────────────────\n",
                 style="dim")
        # If silent mutations were applied to remove internal sites, call
        # them out so the user knows to order the (mutated) insert as a
        # gBlock rather than PCR from the raw template.
        muts = d.get("mutations") or []
        if muts:
            t.append(
                f"\n[{len(muts)}] silent mutation(s) applied to remove "
                f"internal BsaI / Esp3I sites:\n",
                style="bold yellow",
            )
            for m in muts:
                t.append(f"  · {m}\n", style="yellow")
            t.append(
                "  [dim]The insert shown below (and the saved part) reflects "
                "these changes.\n  Order as a gBlock — the primers will not "
                "introduce these mutations during PCR.[/dim]\n",
                style="dim",
            )
        # Flag mutations landing inside primer binding windows. The primers
        # are designed against the mutated insert, so they won't bind the
        # user's original template there — a gBlock is mandatory.
        br_muts = d.get("binding_region_mutations") or []
        if br_muts:
            t.append(
                f"\n⚠ [{len(br_muts)}] mutation(s) land inside primer "
                f"binding region(s):\n",
                style="bold red",
            )
            for entry in br_muts:
                region = "5′ (forward)" if entry["region"] == "fwd" else "3′ (reverse)"
                t.append(f"  · {region}: {entry['text']}\n", style="red")
            t.append(
                "  [dim]The original plasmid CANNOT be used as the PCR "
                "template — the primers would\n  mismatch at the mutated "
                "bases. Order the mutated insert as a gBlock and PCR\n  "
                "from that, or redesign to avoid silent mutations inside "
                "the binding windows.[/dim]\n",
                style="dim",
            )
        # Show every designed primer pair. For the current single-amplicon
        # design this is one; when SOE-PCR splitting is added later, pairs
        # will contain N+1 entries (one per sub-amplicon).
        gtail = (
            grammar.get("pad", _GB_PAD)
            + grammar.get("site", _GB_L0_ENZYME_SITE)
            + grammar.get("spacer", _GB_SPACER)
        )
        tail_len = len(gtail) + 4
        enzyme_label = (grammar.get("enzyme", "?") or "?")[:5]
        legend = (
            f"  {'─'*4}{enzyme_label + '─':>7}{'─OH':>3}"
            f"{'─── binding region':>20}\n"
        )
        n_pairs = len(pairs)
        for i, p in enumerate(pairs, start=1):
            if n_pairs > 1:
                t.append(f"\n── Pair {i} of {n_pairs} ──\n", style="bold cyan")
            t.append(f"\nPair {i} Forward (5'→3'):\n", style="bold green")
            t.append(f"  {p['fwd_full'][:tail_len]}", style="dim green")
            t.append(p["fwd_full"][tail_len:], style="bold green")
            t.append(f"   Tm {p['fwd_tm']:.1f}°C\n", style="dim")
            t.append(legend, style="dim")
            t.append(f"\nPair {i} Reverse (5'→3'):\n", style="bold red")
            t.append(f"  {p['rev_full'][:tail_len]}", style="dim red")
            t.append(p["rev_full"][tail_len:], style="bold red")
            t.append(f"   Tm {p['rev_tm']:.1f}°C\n", style="dim")
            t.append(legend, style="dim")
            t.append(f"\nAmplicon: {p['amplicon_len']} bp\n", style="white")
        t.append(f"\nInsert: {len(d['insert_seq'])} bp   "
                 f"{n_pairs} primer pair(s) total\n",
                 style="white")
        status.update(t)
        self.query_one("#btn-dom-save", Button).disabled = False
        self.query_one("#btn-dom-save-primers", Button).disabled = False

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
        insert = d["insert_seq"]
        oh5    = d["oh5"]
        oh3    = d["oh3"]
        grammar = _get_active_grammar()
        part = {
            "name":        name,
            "type":        d["part_type"],
            "position":    d["position"],
            "oh5":         oh5,
            "oh3":         oh3,
            "backbone":    "pUPD2",
            "marker":      "Spectinomycin",
            "sequence":    insert,
            "fwd_primer":  d["fwd_full"],
            "rev_primer":  d["rev_full"],
            "fwd_tm":      d["fwd_tm"],
            "rev_tm":      d["rev_tm"],
            "primed_seq":  _simulate_primed_amplicon(
                insert, oh5, oh3, grammar=grammar,
            ),
            "cloned_seq":  _simulate_cloned_plasmid(insert, oh5, oh3),
            "grammar":     grammar.get("id", "gb_l0"),
        }
        self.dismiss(part)

    # ── Save primers to library ────────────────────────────────────────────

    @on(Button.Pressed, "#btn-dom-save-primers")
    def _save_primers_to_library(self, _) -> None:
        """Persist every designed primer pair to primers.json.

        Naming follows the project-wide convention:
            {partName}-DOM-{n}-F / {partName}-DOM-{n}-R
        where DOM tags the primer type (domestication, vs CLO cloning /
        DET detection) and {n} is the 1-indexed pair number within this
        domestication run.
        """
        status = self.query_one("#dom-primer-results", Static)
        if self._design is None:
            status.update("[red]Design primers first.[/red]")
            return
        part_name = self.query_one("#dom-name", Input).value.strip()
        if not part_name:
            status.update("[red]Enter a part name before saving primers.[/red]")
            return
        pairs = self._design.get("pairs") or []
        if not pairs:
            status.update("[red]No primer pairs to save.[/red]")
            return

        source = part_name  # domesticator parts don't carry a plasmid context
        today = _date.today().isoformat()
        entries = _load_primers()
        existing_seqs = {e.get("sequence", "").upper() for e in entries}

        new_rows: list[dict] = []
        dupes: list[str] = []
        for idx, p in enumerate(pairs, start=1):
            fwd_name = f"{part_name}-DOM-{idx}-F"
            rev_name = f"{part_name}-DOM-{idx}-R"
            fwd_seq = p["fwd_full"]
            rev_seq = p["rev_full"]
            if fwd_seq.upper() in existing_seqs:
                dupes.append(fwd_name)
            else:
                existing_seqs.add(fwd_seq.upper())
                new_rows.append({
                    "name":        fwd_name,
                    "sequence":    fwd_seq,
                    "tm":          p["fwd_tm"],
                    "primer_type": "goldenbraid",
                    "source":      source,
                    "pos_start":   p["fwd_pos"][0],
                    "pos_end":     p["fwd_pos"][1],
                    "strand":      1,
                    "date":        today,
                    "status":      "Designed",
                })
            if rev_seq.upper() in existing_seqs:
                dupes.append(rev_name)
            else:
                existing_seqs.add(rev_seq.upper())
                new_rows.append({
                    "name":        rev_name,
                    "sequence":    rev_seq,
                    "tm":          p["rev_tm"],
                    "primer_type": "goldenbraid",
                    "source":      source,
                    "pos_start":   p["rev_pos"][0],
                    "pos_end":     p["rev_pos"][1],
                    "strand":      -1,
                    "date":        today,
                    "status":      "Designed",
                })

        if not new_rows:
            self.app.notify(
                f"All {len(pairs) * 2} primer sequences already exist in the "
                f"library — nothing saved.",
                severity="warning", timeout=8,
            )
            return

        # Prepend new rows so the library table surfaces them first.
        by_name = {r["name"] for r in new_rows}
        entries = [e for e in entries if e.get("name") not in by_name]
        entries = new_rows + entries
        _save_primers(entries)
        msg = f"Saved {len(new_rows)} primer(s) to library."
        if dupes:
            msg += f" Skipped {len(dupes)} duplicate(s): {', '.join(dupes)}"
        self.app.notify(msg, timeout=8)
        self.query_one("#btn-dom-save-primers", Button).disabled = True

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
        # Source-dropdown options computed once here so on_mount can read
        # the same default the Select was rendered with. The "Current map
        # features" option only appears when a plasmid is actually loaded;
        # otherwise the modal is launchable from a blank canvas (lib /
        # parts / prot still produce a CDS).
        self._src_options: list[tuple[str, str]] = []
        if self._template:
            self._src_options.append(("Current map features", "map"))
        self._src_options.extend([
            ("Plasmid library",              "lib"),
            ("Parts bin",                    "parts"),
            ("Protein sequence (harmonize)", "prot"),
        ])
        self._initial_source = self._src_options[0][1]

    def compose(self) -> ComposeResult:
        with Vertical(id="mut-box"):
            yield Static(
                " Mutagenize  —  Golden Braid SOE-PCR Site-Directed Mutagenesis ",
                id="mut-title",
            )
            yield Label("CDS source")
            yield Select(
                self._src_options,
                id="mut-source", value=self._initial_source, allow_blank=False,
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

            with Vertical(id="mut-src-parts"):
                yield Label("Part  (from your Parts Bin)")
                yield Select(self._build_parts_options(),
                             id="mut-parts", prompt="(select a part)")

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
                yield Button("Change", id="btn-mut-codon", variant="default")

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

    def _build_parts_options(self) -> list:
        """Build Select options for the Parts Bin source. Filters to parts
        whose stored insert sequence is a valid CDS shape (≥ 30 bp, length
        a multiple of 3) — same gate the map/library sources apply. The
        value is the integer index into ``_load_parts_bin()`` (as a str)."""
        try:
            entries = _load_parts_bin()
        except Exception:
            _log.exception("Mutagenize: failed to load parts bin")
            entries = []
        opts: list = []
        for i, e in enumerate(entries):
            seq = (e.get("sequence") or "").upper()
            if len(seq) < 30 or len(seq) % 3 != 0:
                continue
            nm = e.get("name", f"part_{i}")
            ptype = e.get("type", "?")
            opts.append((f"{nm}  [{ptype}, {len(seq)} bp]", str(i)))
        if not opts:
            opts = [("(no eligible parts in bin)", "_none")]
        return opts

    def on_mount(self) -> None:
        # Seed built-in K12 via registry load
        try:
            _codon_tables_load()
            self._codon_entry = _codon_tables_get("83333")
        except Exception:
            _log.exception("Mutagenize: codon registry load failed")
            self._codon_entry = None
        # Falls back to lib/parts/prot if no plasmid is loaded — the source
        # dropdown excludes "map" in that case (see compose).
        self._apply_source(self._initial_source)
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
        self.query_one("#mut-src-map", Vertical).display   = (src == "map")
        self.query_one("#mut-src-lib", Vertical).display   = (src == "lib")
        self.query_one("#mut-src-parts", Vertical).display = (src == "parts")
        self.query_one("#mut-src-prot", Vertical).display  = (src == "prot")
        # Switching source invalidates any previously-loaded CDS and the
        # last designed primers — otherwise the preview keeps showing a
        # stale mutation from the previous source.
        self._reset_cds_state()

    def _reset_cds_state(self, info_msg: str = "") -> None:
        """Clear CDS + primer state and refresh the preview. Used by
        source-switch and by a deselect on any of the CDS dropdowns."""
        self._cds_dna  = ""
        self._cds_meta = None
        self._outer    = None
        self._inner    = None
        try:
            self.query_one("#mut-cds-info", Static).update(info_msg)
            self.query_one("#btn-mut-save", Button).disabled = True
        except NoMatches:
            return
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

    # ── Parts-bin source ──────────────────────────────────────────────────

    @on(Select.Changed, "#mut-parts")
    def _parts_changed(self, event: Select.Changed) -> None:
        val = event.value
        info = self.query_one("#mut-cds-info", Static)
        if not isinstance(val, str) or val == "_none" or not val.isdigit():
            self._reset_cds_state()
            return
        try:
            entries = _load_parts_bin()
        except Exception:
            _log.exception("Mutagenize: failed to load parts bin")
            info.update("[red]Could not load Parts Bin.[/red]")
            return
        try:
            entry = entries[int(val)]
        except (IndexError, ValueError):
            info.update("[red]Parts bin entry not found.[/red]")
            return
        seq = (entry.get("sequence") or "").upper()
        if len(seq) < 30 or len(seq) % 3 != 0:
            info.update("[red]Part sequence is too short or not a multiple "
                        "of 3.[/red]")
            return
        name = entry.get("name") or f"part_{val}"
        self._plasmid_name = name
        # A part's insert is already in 5'→3' biological orientation, so
        # treat it as a single-CDS pseudo-plasmid spanning [0, len(seq)).
        synthetic_feat = {
            "type": "CDS", "label": name, "strand": 1,
            "start": 0, "end": len(seq),
        }
        self._load_cds_from_feature(f"0:{len(seq)}:1", seq, [synthetic_feat],
                                    origin="parts")

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
        today = _date.today().isoformat()
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
        # Suffix tags the primer role: DET = detection (diagnostic PCR),
        # CLO = cloning (with RE tails), DOM = Golden Braid L0 domestication.
        if primer_type == "detection":
            suffix = "DET"
        elif primer_type == "goldenbraid":
            suffix = "DOM"
        else:
            suffix = "CLO"
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

        today = _date.today().isoformat()

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


class UnsavedNavigateModal(ModalScreen):
    """Shown when the user tries to navigate (e.g. Back to Collections)
    with unsaved edits. Sibling of `UnsavedQuitModal` — kept separate
    because the button labels and verb differ ("go back" vs "quit"),
    and the wording matters for users to understand the consequence.

    Dismisses with ``"save"`` (caller saves then proceeds), ``"discard"``
    (caller reverts the in-memory record from the library copy then
    proceeds), or ``None`` (cancel — stay).
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "focus_next", "Next button", show=False),
    ]

    def __init__(self, action_phrase: str = "leave"):
        super().__init__()
        self._action_phrase = action_phrase

    def compose(self) -> ComposeResult:
        with Vertical(id="navunsv-dlg"):
            yield Static(" Unsaved Changes ", id="navunsv-title")
            yield Static(
                f"  The loaded plasmid has unsaved edits.\n"
                f"  Save before you {self._action_phrase}?",
                id="navunsv-msg",
            )
            with Horizontal(id="navunsv-btns"):
                yield Button("Save",            id="btn-navunsv-save",
                             variant="primary")
                yield Button("Discard Changes", id="btn-navunsv-discard",
                             variant="error")
                yield Button("Cancel",          id="btn-navunsv-cancel")

    @on(Button.Pressed, "#btn-navunsv-save")
    def _save(self, _):     self.dismiss("save")

    @on(Button.Pressed, "#btn-navunsv-discard")
    def _discard(self, _):  self.dismiss("discard")

    @on(Button.Pressed, "#btn-navunsv-cancel")
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
        self.dismiss(_cursor_row_key(self.query_one("#pick-table", DataTable)))

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


class CollectionsModal(ModalScreen):
    """Browse, save, and load named collections of plasmids.

    A collection is a saved snapshot of the plasmid library — the user
    can keep several themed sets (e.g. yeast project, E. coli toolkit,
    MoClo plant) and switch between them. Loading a collection replaces
    the current library wholesale; ``_save_library`` writes a `.bak`
    via the same atomic-save invariant the rest of the app uses, so the
    previous library can be recovered manually if the swap was a mistake.

    Dismisses with ``{"loaded": name, "n_plasmids": int}`` when a
    collection was loaded, else None.
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Close"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def compose(self) -> ComposeResult:
        with Vertical(id="coll-dlg"):
            yield Static(" Plasmid Collections ", id="coll-title")
            yield Label("Save the current library as a new collection:")
            with Horizontal(id="coll-save-row"):
                yield Input(placeholder="Collection name",
                            id="coll-save-name")
                yield Button("Save", id="btn-coll-save", variant="primary")
            yield Static("", id="coll-status", markup=True)
            yield Label("Existing collections:")
            yield DataTable(id="coll-table", cursor_type="row",
                            zebra_stripes=True)
            with Horizontal(id="coll-btns"):
                yield Button("Load Selected", id="btn-coll-load",
                             variant="warning")
                yield Button("Delete", id="btn-coll-del", variant="error")
                yield Button("Close", id="btn-coll-close")

    def on_mount(self) -> None:
        t = self.query_one("#coll-table", DataTable)
        t.add_columns("Name", "# Plasmids", "Description")
        self._repopulate()
        self.query_one("#coll-save-name", Input).focus()

    def _repopulate(self) -> None:
        t = self.query_one("#coll-table", DataTable)
        t.clear()
        for c in _load_collections():
            name = c.get("name") or "?"
            n_plas = len(c.get("plasmids", []) or [])
            desc = (c.get("description") or "")[:40]
            t.add_row(name, str(n_plas), desc, key=name)

    @on(Button.Pressed, "#btn-coll-save")
    def _save(self, _) -> None:
        name = self.query_one("#coll-save-name", Input).value.strip()
        status = self.query_one("#coll-status", Static)
        if not name:
            status.update("[red]Enter a collection name.[/red]")
            return
        if _collection_name_taken(name):
            status.update(
                f"[red]A collection named '{name}' already exists.[/red]"
            )
            return
        plasmids = _load_library()
        existing = _load_collections()
        existing.append({
            "name":        name,
            "description": f"Saved {len(plasmids)} plasmid(s)",
            "plasmids":    plasmids,
            "saved":       _date.today().isoformat(),
        })
        _save_collections(existing)
        self.query_one("#coll-save-name", Input).value = ""
        status.update(
            f"[green]Saved '{name}' ({len(plasmids)} plasmid(s)).[/green]"
        )
        self._repopulate()

    @on(Input.Submitted, "#coll-save-name")
    def _save_submitted(self, _) -> None:
        self._save(None)

    @on(Button.Pressed, "#btn-coll-load")
    def _load(self, _) -> None:
        status = self.query_one("#coll-status", Static)
        name = _cursor_row_key(self.query_one("#coll-table", DataTable))
        if not name:
            status.update("[red]No collection selected.[/red]")
            return
        coll = _find_collection(name)
        if coll is None:
            status.update(f"[red]Collection '{name}' not found.[/red]")
            return
        plasmids = [p for p in (coll.get("plasmids") or [])
                    if isinstance(p, dict)]
        # Set active BEFORE the library write so the upcoming sync mirror
        # targets the newly-selected collection. Since `plasmids` came
        # from `coll`, the mirror writes the same content back — a true
        # no-op only because we read `plasmids` before the mirror runs.
        _set_active_collection_name(name)
        _save_library(plasmids)
        self.dismiss({"loaded": name, "n_plasmids": len(plasmids)})

    @on(Button.Pressed, "#btn-coll-del")
    def _delete(self, _) -> None:
        status = self.query_one("#coll-status", Static)
        name = _cursor_row_key(self.query_one("#coll-table", DataTable))
        if not name:
            status.update("[red]No collection selected.[/red]")
            return
        existing = [c for c in _load_collections() if c.get("name") != name]
        _save_collections(existing)
        self._repopulate()
        status.update(f"[dim]Deleted collection '{name}'.[/dim]")

    @on(Button.Pressed, "#btn-coll-close")
    def _close_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class CollectionNameModal(ModalScreen):
    """Tiny prompt modal for creating or renaming a collection.

    Dismisses with the trimmed name string, or None on cancel.
    Caller is responsible for collision-checking before persisting.
    """

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, title: str, current: str = "",
                 placeholder: str = "Collection name") -> None:
        super().__init__()
        self.title_text = title
        self.current = current
        self.placeholder_text = placeholder

    def compose(self) -> ComposeResult:
        with Vertical(id="collname-dlg"):
            yield Static(f" {self.title_text} ", id="collname-title")
            yield Label("Name:")
            yield Input(value=self.current,
                        placeholder=self.placeholder_text,
                        id="collname-input")
            yield Static("", id="collname-status", markup=True)
            with Horizontal(id="collname-btns"):
                yield Button("OK",     id="btn-collname-ok",     variant="primary")
                yield Button("Cancel", id="btn-collname-cancel")

    def on_mount(self) -> None:
        self.query_one("#collname-input", Input).focus()

    @on(Button.Pressed, "#btn-collname-ok")
    def _ok(self, _) -> None:
        self._submit()

    @on(Input.Submitted, "#collname-input")
    def _submitted(self, _) -> None:
        self._submit()

    def _submit(self) -> None:
        name = self.query_one("#collname-input", Input).value.strip()
        if not name:
            self.query_one("#collname-status", Static).update(
                "[red]Name cannot be empty.[/red]"
            )
            return
        self.dismiss(name)

    @on(Button.Pressed, "#btn-collname-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    def action_cancel(self) -> None:
        self.dismiss(None)


class CollectionDeleteConfirmModal(ModalScreen):
    """Confirm-on-delete modal for collections — different copy from
    LibraryDeleteConfirmModal (which talks about library entries).

    Default focus on [No] to protect against handslip-deletes.
    Dismisses True (delete) or False (keep)."""

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, name: str, n_plasmids: int) -> None:
        super().__init__()
        self.coll_name = name
        self.n_plas = n_plasmids

    def compose(self) -> ComposeResult:
        plural = "" if self.n_plas == 1 else "s"
        with Vertical(id="colldel-dlg"):
            yield Static(" Delete collection ", id="colldel-title")
            yield Static(
                f"  Delete collection [bold]{self.coll_name}[/bold]?\n"
                f"  ({self.n_plas} plasmid{plural})\n\n"
                f"  [dim]A backup is written to\n"
                f"  collections.json.bak before the change.[/dim]",
                id="colldel-msg", markup=True,
            )
            with Horizontal(id="colldel-btns"):
                yield Button("No",          id="btn-colldel-no",  variant="default")
                yield Button("Yes, delete", id="btn-colldel-yes", variant="error")

    def on_mount(self) -> None:
        self.query_one("#btn-colldel-no", Button).focus()

    @on(Button.Pressed, "#btn-colldel-no")
    def _no(self, _) -> None:
        self.dismiss(False)

    @on(Button.Pressed, "#btn-colldel-yes")
    def _yes(self, _) -> None:
        self.dismiss(True)

    def action_cancel(self) -> None:
        self.dismiss(False)


class ScaryDeleteConfirmModal(ModalScreen):
    """Second-stage confirmation for collection delete — deliberately
    visually loud (red border + warning banner + emphatic copy) to make
    the user pause. Default focus on [No] like every confirm modal in
    the app. Dismisses True (delete) or False (keep)."""

    BINDINGS = [
        Binding("escape", "cancel",     "Cancel"),
        Binding("tab",    "focus_next", "Next", show=False),
    ]

    def __init__(self, name: str, n_plasmids: int) -> None:
        super().__init__()
        self.coll_name = name
        self.n_plas = n_plasmids

    def compose(self) -> ComposeResult:
        plural = "" if self.n_plas == 1 else "s"
        with Vertical(id="scarydel-dlg"):
            yield Static(
                "  ⚠   ARE YOU ABSOLUTELY SURE?   ⚠  ",
                id="scarydel-title", markup=False,
            )
            yield Static(
                f"\n  This will [bold red]permanently delete[/bold red] the "
                f"collection\n"
                f"  [bold]{self.coll_name}[/bold] and its "
                f"[bold red]{self.n_plas} plasmid{plural}[/bold red].\n\n"
                f"  [yellow]The plasmids inside will also be removed from\n"
                f"  the library mirror.[/yellow]\n\n"
                f"  A backup of [italic]collections.json[/italic] is written "
                f"to\n"
                f"  [italic]collections.json.bak[/italic] in your data "
                f"directory —\n"
                f"  recover from it manually if you change your mind.\n",
                id="scarydel-msg", markup=True,
            )
            with Horizontal(id="scarydel-btns"):
                yield Button("No, keep it", id="btn-scarydel-no",
                             variant="default")
                yield Button("Yes, delete forever", id="btn-scarydel-yes",
                             variant="error")

    def on_mount(self) -> None:
        self.query_one("#btn-scarydel-no", Button).focus()

    @on(Button.Pressed, "#btn-scarydel-no")
    def _no(self, _) -> None:
        self.dismiss(False)

    @on(Button.Pressed, "#btn-scarydel-yes")
    def _yes(self, _) -> None:
        self.dismiss(True)

    def action_cancel(self) -> None:
        self.dismiss(False)


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

/* ── FASTA file picker modal ─────────────────────────────── */
FastaFilePickerModal { align: center middle; }
#fasta-box {
    width: 90; max-width: 95%; min-width: 60;
    height: 32; max-height: 90%;
    background: $surface; border: solid $accent; padding: 1 2;
}
#fasta-title  { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#fasta-header { height: 1; margin-bottom: 1; color: $text-muted; }
#fasta-tree   { height: 1fr; border: solid $primary-darken-2; }
#fasta-hint   { height: 1; margin-top: 1; color: $text-muted; }
#fasta-status { height: 1; margin-top: 1; }
#fasta-btns   { height: 3; margin-top: 1; }
#fasta-btns Button { margin-right: 1; }

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

/* ── Export-FASTA modal ──────────────────────────────────── */
FastaExportModal { align: center middle; }
#fasta-export-box {
    width: 72; height: auto;
    background: $surface; border: solid $primary; padding: 1 2;
}
#fasta-export-title { background: $primary; padding: 0 1; margin-bottom: 1; }
#fasta-export-box Label { color: $text-muted; margin-top: 1; }
#fasta-export-btns { height: 3; margin-top: 1; }
#fasta-export-btns Button { margin-right: 1; }
#fasta-export-status { height: 2; margin-top: 1; }

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

/* ── Unsaved-navigate dialog (back to collections, switch tab, ...) ── */
UnsavedNavigateModal { align: center middle; }
#navunsv-dlg {
    width: 60; height: auto;
    background: $surface; border: solid $warning; padding: 1 2;
}
#navunsv-title { background: $warning-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#navunsv-msg   { color: $text-muted; margin-bottom: 1; }
#navunsv-btns  { height: 3; margin-top: 1; }
#navunsv-btns Button { margin-right: 1; }

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

/* ── Plasmid collections modal ───────────────────────────── */
CollectionsModal { align: center middle; }
#coll-dlg {
    width: 100; height: 32;
    background: $surface; border: solid $primary; padding: 1 2;
}
#coll-title    { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#coll-dlg Label { color: $text-muted; margin-top: 1; }
#coll-save-row { height: 3; margin-top: 1; }
#coll-save-row Input  { width: 1fr; margin-right: 1; }
#coll-save-row Button { width: 14; }
#coll-status   { height: 1; margin: 0 0 1 0; }
#coll-table    { height: 1fr; }
#coll-btns     { height: 3; margin-top: 1; }
#coll-btns Button { margin-right: 1; min-width: 12; }

/* ── Collection name prompt + delete confirm ───────────────── */
CollectionNameModal { align: center middle; }
#collname-dlg {
    width: 60; height: 14;
    background: $surface; border: solid $primary; padding: 1 2;
}
#collname-title  { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#collname-dlg Label { color: $text-muted; margin-top: 1; }
#collname-input  { margin-top: 1; }
#collname-status { height: 1; margin-top: 1; }
#collname-btns   { height: 3; margin-top: 1; }
#collname-btns Button { margin-right: 1; min-width: 10; }

CollectionDeleteConfirmModal { align: center middle; }
#colldel-dlg {
    width: 60; height: 16;
    background: $surface; border: solid $primary; padding: 1 2;
}
#colldel-title { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#colldel-msg   { height: 1fr; }
#colldel-btns  { height: 3; margin-top: 1; }
#colldel-btns Button { margin-right: 1; min-width: 12; }

/* ── Scary second-confirm: deliberately loud, can't miss it ──── */
ScaryDeleteConfirmModal { align: center middle; }
#scarydel-dlg {
    width: 70; height: 20;
    background: $surface;
    border: thick $error;
    padding: 1 2;
}
#scarydel-title {
    background: $error;
    color: $text;
    text-style: bold;
    text-align: center;
    padding: 0 1;
    margin-bottom: 1;
}
#scarydel-msg  { height: 1fr; }
#scarydel-btns { height: 3; margin-top: 1; }
#scarydel-btns Button { margin-right: 1; min-width: 18; }

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
#addfeat-color-row       { height: 3; margin-top: 1; }
#addfeat-color-label     { width: 8; margin-top: 1; }
#addfeat-color-swatch    { width: 1fr; margin-top: 1; }
#addfeat-color-row Button { margin-right: 1; min-width: 14; }
#addfeat-seq    { height: 6; min-height: 4; margin-top: 0; }
#addfeat-status { height: 2; margin-top: 1; }
#addfeat-btns   { height: 3; margin-top: 1; }
#addfeat-btns Button { margin-right: 1; }

/* ── Color picker modal ──────────────────────────────────── */
ColorPickerModal { align: center middle; }
#colorpick-dlg {
    width: 120; height: 90%; max-height: 40;
    background: $surface; border: solid $primary; padding: 1 2;
}
#colorpick-title   { background: $primary-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#colorpick-current { height: 1; margin-bottom: 0; }
#colorpick-capability { height: 1; margin-bottom: 1; }
#colorpick-preview-row    { height: 5; margin-bottom: 1; align: left middle; }
#colorpick-preview-swatch { width: 24; height: 5; border: tall $primary; margin-right: 2; }
#colorpick-preview-label  { width: 1fr; height: 5; content-align: left middle; }
#colorpick-scroll  { height: 1fr; overflow-y: auto; }
.colorpick-section-hdr { background: $primary-darken-2; padding: 0 1; margin-top: 1; height: 1; color: $text; }
#colorpick-row     { height: 3; margin-bottom: 0; }
.colorpick-swatch  { min-width: 4; width: 4; height: 3; margin: 0 0 0 0; border: none; }
#colorpick-xterm-grid { height: auto; margin-bottom: 1; }
.colorpick-xterm-row { height: 1; }
.colorpick-xterm-cell { min-width: 3; width: 3; height: 1; margin: 0 0 0 0; border: none; padding: 0; }
#colorpick-custom-row { height: 3; margin-bottom: 1; }
#colorpick-custom-label { width: 18; margin-top: 1; }
#colorpick-hex-input { width: 1fr; margin-right: 1; }
#colorpick-custom-row Button { margin-right: 1; min-width: 10; }
#colorpick-status  { height: 1; margin: 0 0 1 0; }
#colorpick-btns    { height: 3; }
#colorpick-btns Button { margin-right: 1; }

/* ── Feature library full-screen ─────────────────────────── */
#flib-box {
    width: 100%; height: 1fr;
    background: $surface; padding: 0 2;
}
#flib-title        { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#flib-main         { height: 1fr; }
#flib-left         { width: 1fr; padding-right: 1; border-right: solid $primary-darken-2; }
#flib-right        { width: 2fr; padding-left: 1; }
.flib-section-hdr  { background: $primary-darken-2; padding: 0 1; }
#flib-table        { height: 1fr; }
#flib-btns         { height: 3; margin-top: 1; }
#flib-btns Button  { margin-right: 1; }

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

/* ── Grammar editor modal ────────────────────────────────── */
GrammarEditorModal { align: center middle; }
#ged-dlg {
    width: 110; max-width: 95%; min-width: 80;
    height: 90%; max-height: 44;
    background: $surface; border: solid $accent; padding: 1 2;
}
#ged-title { background: $primary-darken-1; color: $text; padding: 0 1; }
#ged-builtin-banner {
    background: $warning-darken-2; color: $text;
    padding: 0 1; margin-top: 1;
}
#ged-body { height: 1fr; padding: 0 1; }
#ged-enzyme-row, #ged-tail-row { height: 3; margin-top: 1; }
#ged-enzyme-row Input, #ged-tail-row Input { width: 22; margin-right: 2; }
.ged-inline-label { padding: 1 1 0 0; width: 14; }
#ged-forbidden { height: 5; border: solid $primary-darken-2; }
#ged-positions { height: 9; border: solid $primary-darken-2; }
#ged-status { height: auto; min-height: 1; padding: 0 1; }
#ged-btns { height: 3; margin-top: 1; }
#ged-btns Button { margin-right: 1; }

/* ── Parts bin (full-screen) ─────────────────────────────── */
#parts-box {
    width: 100%; height: 1fr;
    background: $surface; padding: 0 2;
}
#parts-title  { background: $success-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#parts-table  { height: 1fr; }
#parts-detail { height: 3; border-top: solid $accent; padding: 0 1; color: $text-muted; }
#parts-seq-view {
    height: 10; border: solid $accent; padding: 0 1;
    background: $surface-darken-1;
}
/* Single bottom button row holds Copy + action buttons together. */
#parts-btns   { height: 3; margin-top: 1; }
#parts-btns Button { margin-right: 1; }

/* ── Domesticator modal ─────────────────────────────────── */
DomesticatorModal { align: center middle; }
#dom-box {
    width: 110; max-width: 95%; min-width: 80;
    height: 90%; max-height: 46;
    background: $surface; border: solid $accent; padding: 1 2;
}
#dom-title  { background: $accent-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#dom-body   { height: 1fr; overflow-y: auto; }
#dom-grammar-row { height: 3; margin-bottom: 1; }
#dom-grammar-row Label { padding: 1 1 0 0; width: 18; }
#dom-grammar-row Select { width: 1fr; }
#dom-row1   { height: 5; }
#dom-name-col { width: 1fr; padding-right: 1; }
#dom-type-col { width: 1fr; }
#dom-oh-info  { height: 1; margin-bottom: 1; }
#dom-codon-row { height: 3; margin-bottom: 1; }
#dom-codon-label { width: 4fr; padding: 1 1 0 1; }
#dom-codon-row Button { width: 1fr; }
/* Source picker: four radios on ONE row, no scrollbar, flex with modal width */
#dom-src    {
    layout: horizontal;
    height: 3; width: 100%;
    margin-bottom: 1;
    overflow: hidden;
}
#dom-src > RadioButton {
    width: 1fr; height: 1;
    margin: 0 1 0 0; padding: 0;
    background: transparent; border: none;
}
.dom-src-panel { width: 100%; height: auto; margin-bottom: 1; }
#dom-direct-seq { width: 100%; height: 6; }
#dom-plasmid-hdr, #dom-fasta-hdr { height: 3; width: 100%; }
#dom-plasmid-hdr Label, #dom-fasta-hdr Label {
    width: auto; padding-right: 1; content-align: left middle;
}
#dom-plasmid-name, #dom-fasta-name {
    width: 1fr; content-align: left middle; color: $text;
}
#dom-plasmid-hdr Button, #dom-fasta-hdr Button { margin-left: 1; }
#dom-featlib-select, #dom-plasmid-feat-select { width: 100%; }
#dom-featlib-preview, #dom-plasmid-feat-preview, #dom-fasta-preview {
    height: 1; width: 100%;
}
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
#mut-src-map, #mut-src-lib, #mut-src-parts, #mut-src-prot { height: auto; }
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
        Binding("ctrl+o",      "open_file",        "Open file",     show=True),
        Binding("ctrl+shift+a","add_to_library",   "Add to lib",    show=True),
        Binding("ctrl+e",      "edit_seq",         "Edit seq",      show=True),
        Binding("ctrl+s",      "save",             "Save",          show=True),
        Binding("ctrl+f",      "add_feature",      "Add feature",   show=True),
        Binding("ctrl+shift+f","capture_to_features", "→ Feat lib", show=True,  priority=True),
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
    # self.screen.id; only app-level (main-screen) shortcuts like f,
    # Ctrl+O, Ctrl+Shift+A, A, Ctrl+E, r, v, etc. get blocked.
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
        # Migration must run BEFORE children compose+mount — Textual fires
        # mount events leaves→root, so anything in App.on_mount runs AFTER
        # LibraryPanel.on_mount, by which point the panel has already read
        # `active_collection` (None) and locked into the wrong view mode.
        # Both helpers are idempotent.
        _ensure_default_collection()
        _restore_library_from_active_collection()
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
                "   f fetch   ^O open   ^S save   ^E edit   ^F add-feat   ^⇧F →feat-lib   ^⇧A add-to-lib",
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
        # Focus-aware routing: when the keyboard focus is inside the
        # library, Delete targets whatever's under the cursor in the
        # current view — a plasmid (one-stage confirm) or a collection
        # (two-stage confirm with a loud red second modal). The same
        # entry point is used by the panel's `−` button so button and
        # keyboard behave identically.
        if self._focus_is_in_library():
            self.query_one("#library", LibraryPanel).request_delete_under_cursor()
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
        # Crash-recovery autosave: debounced so rapid edits coalesce into one
        # write. Cleared whenever the record is saved / marked clean.
        self._autosave_timer = None
        self._AUTOSAVE_DEBOUNCE_S = 3.0
        # Validate all user-data files before anything else. Corrupt files
        # are auto-restored from .bak if possible; the user is notified
        # either way so they know the state of their data.
        self._check_data_files()
        self._check_crash_recovery()
        # Migration to the collection-driven model already ran in compose
        # (so child panels see the correct active collection on mount).
        if self._preload_record is not None:
            def _load_preload():
                self._import_and_persist(self._preload_record)
            self.call_after_refresh(_load_preload)
        else:
            lib = _load_library()
            if lib:
                # User has at least one library entry — auto-load the first
                # so the canvas isn't blank on startup. Falls through silently
                # if the entry's gb_text is missing or unparsable; the user
                # can still pick another row from the panel.
                first = lib[0]
                gb_text = first.get("gb_text", "")
                if gb_text:
                    try:
                        record = _gb_text_to_record(gb_text)
                        def _load_first(r=record):
                            self._apply_record(r)
                        self.call_after_refresh(_load_first)
                    except Exception:
                        _log.exception(
                            "Auto-load of first library entry %r failed",
                            first.get("name", "?"),
                        )
            else:
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
        global _collections_cache, _settings_cache
        # Force a cold read on every JSON registry so a corrupt file is
        # detected NOW (with .bak recovery + a user notify) rather than at
        # the first lazy-load when something breaks downstream.
        for path, label, cache_attr in [
            (_LIBRARY_FILE,     "Plasmid library",     "_library_cache"),
            (_PARTS_BIN_FILE,   "Parts bin",           "_parts_bin_cache"),
            (_PRIMERS_FILE,     "Primer library",      "_primers_cache"),
            (_COLLECTIONS_FILE, "Plasmid collections", "_collections_cache"),
            (_SETTINGS_FILE,    "Settings",            "_settings_cache"),
        ]:
            globals()[cache_attr] = None
            _, warning = _safe_load_json(path, label)
            if warning:
                self.notify(warning, severity="warning", timeout=12)
        # Caches stay None so the next typed loader (e.g. _load_settings)
        # rebuilds them through its own filter/shape logic.

    # ── Crash-recovery autosave ────────────────────────────────────────────────

    def _autosave_path(self, record) -> "Path | None":
        """Return the autosave file path for `record`, or None if no id.

        The filename is `{safe}-{hash6}.gb` where `safe` is the sanitised
        record.id and `hash6` is a 6-char hex of sha256(record.id). Two
        records with `id` like 'foo/bar' and 'foo_bar' both sanitise to
        'foo_bar' but get distinct hashes, so they no longer overwrite each
        other in the crash-recovery directory.
        """
        if record is None or not getattr(record, "id", ""):
            return None
        import hashlib
        safe = re.sub(r'[^A-Za-z0-9._-]', '_', record.id)[:80]
        if not safe:
            return None
        h = hashlib.sha256(record.id.encode("utf-8")).hexdigest()[:6]
        return _CRASH_RECOVERY_DIR / f"{safe}-{h}.gb"

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
            _atomic_write_text(
                path, _record_to_gb_text(self._current_record),
            )
            _log.info("Autosaved %s to %s (%d bp)",
                      self._current_record.name, path,
                      len(self._current_record.seq))
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
        # Skip if the user is arrowing through one of the navigable
        # DataTables (feature sidebar, plasmid library, etc.). Otherwise
        # this handler races with the table's own arrow handler — every
        # Down keystroke would BOTH advance the seq cursor here AND fire
        # the table's row-highlight cascade, scrolling the sequence panel
        # twice (visible jitter on sidebar-feature browsing). The
        # SequencePanel itself is not focusable, so we can't use
        # `sp.has_focus`; instead, check whether some OTHER focused
        # widget owns the keystroke.
        focused = self.focused
        if focused is not None:
            from textual.widgets import DataTable
            if isinstance(focused, DataTable):
                return
        if sp._cursor_pos < 0 or not sp._seq:
            return
        n  = len(sp._seq)
        k  = event.key
        lw = sp._line_width()
        if k in ("left", "shift+left"):
            new_pos = max(0, sp._cursor_pos - 1)
        elif k in ("right", "shift+right"):
            # Cap at n (one past last base) so the dialog's "insert at cursor"
            # path can append. cursor_pos in [0, n] inclusive matches Python's
            # half-open slicing convention used by `_edit_dialog_result`.
            new_pos = min(n, sp._cursor_pos + 1)
        elif k in ("up", "shift+up"):
            new_pos = max(0, sp._cursor_pos - lw)
        elif k in ("down", "shift+down"):
            # Cap at n - 1 (last base). Down is meant to navigate visually
            # row-to-row; on the last row it should land on the last
            # basepair, not on n (which has no base to highlight and looks
            # like the cursor disappeared off-screen). Insert-at-end is
            # still reachable via Right arrow.
            new_pos = min(max(0, n - 1), sp._cursor_pos + lw)
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
        # Order matters: scroll BEFORE refresh. Textual's `view.update()`
        # tick resets `scroll_y` when content changes, so a sync set
        # afterwards gets clobbered. Setting scroll_y before the refresh
        # means the new content paints at the already-correct viewport.
        sp._ensure_cursor_visible()
        sp._refresh_view()

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

    def _discard_changes(self) -> None:
        """Revert the in-memory record to whatever the library has stored.

        Used by `UnsavedNavigateModal` when the user picks "Discard
        Changes" — reloads the saved copy and clears the undo stack so
        the discarded edits cannot be revived via Ctrl+Z. If the record
        was never saved (no library entry), just marks clean.
        """
        if self._current_record is None:
            return
        record_id = getattr(self._current_record, "id", None)
        match = next(
            (e for e in _load_library() if e.get("id") == record_id),
            None,
        ) if record_id else None
        if match is None or not match.get("gb_text"):
            self._mark_clean()
            return
        try:
            record = _gb_text_to_record(match["gb_text"])
        except Exception:
            _log.exception("Discard reload failed for %r", record_id)
            self.notify("Failed to revert from library; marking clean.",
                        severity="warning")
            self._mark_clean()
            return
        # Wipe undo so the discarded edits can't be reapplied.
        self._undo_stack.clear()
        self._redo_stack.clear()
        # _apply_record(clear_undo=True) would null _source_path — preserve
        # the file path so a subsequent Ctrl+S still targets the original
        # .gb file. The undo-stack clear above already covers the discard
        # semantics; we just need an in-place record swap here.
        saved_source = self._source_path
        self._apply_record(record, clear_undo=False)
        self._source_path = saved_source

    def _do_save(self) -> bool:
        """Save current record to its source file and/or library. Returns True on success."""
        if self._current_record is None:
            self.notify("Nothing to save.", severity="warning")
            return False

        # Write to source file if one is known
        if self._source_path:
            try:
                _atomic_write_text(
                    Path(self._source_path),
                    _record_to_gb_text(self._current_record),
                )
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
        open local file, and CLI preload. Library loads and undo/redo go
        through `_apply_record` directly so they don't re-save the same record.

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

        `clear_undo=False` is for in-place record changes (primer-add,
        feature-merge) — the stacks stay intact and the edit remains undo-able.
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

    def _wrap_aware_midpoint(self, f: dict, n: int) -> int:
        """Return the bp at the visual midpoint of a feature, honouring the
        circular `end < start` wrap convention. A naive `(start + end) // 2`
        puts the midpoint on the OPPOSITE side of the plasmid for wrap features."""
        s, e = f["start"], f["end"]
        if not n:
            return 0
        arc_len = (e - s) % n
        return (s + arc_len // 2) % n

    def _focus_feature(self, f: "dict | None", bp: int) -> None:
        """Single-source UX for "user picked a feature": highlight it as a
        copyable user_sel range, place the cursor on `bp`, and centre the
        sequence panel on `bp`. Called from every feature-pick entry point
        (sequence-panel lane click, plasmid-map click, sidebar row click)
        so the three feel identical.

        Backbone clicks pass `f=None` — only the centring + clearing of
        prior highlights applies in that case.

        Both highlight calls go in with `scroll=False`: their own
        `_ensure_cursor_visible` would scroll the cursor "just into view"
        first, then `center_on_bp` would scroll again to centre. The two
        scrolls in sequence read as a perceptible jitter when arrowing
        through the sidebar feature list.
        """
        seq_pnl = self.query_one("#seq-panel", SequencePanel)
        if f is not None and bp >= 0:
            seq_pnl.select_feature_range(f, cursor_bp=bp, scroll=False)
        else:
            seq_pnl.highlight_feature(f, scroll=False)
        if bp >= 0:
            seq_pnl.center_on_bp(bp)

    @on(SequencePanel.SequenceClick)
    def _seq_click(self, event: SequencePanel.SequenceClick) -> None:
        """Click on the sequence panel.

        - **Lane art click** (`from_lane=True`): the user clicked a feature
          bar / arrow / label, so we highlight the whole feature's DNA
          and surface it in the sidebar / map.
        - **DNA-row click** (`from_lane=False`): the user clicked a base —
          just move the cursor there and clear any whole-feature
          highlight. The base may live inside a feature, but the user
          asked for a single-base operation, not a feature pick.
        """
        pm      = self.query_one("#plasmid-map", PlasmidMap)
        sidebar = self.query_one("#sidebar",     FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel",   SequencePanel)
        bp      = event.bp

        if not event.from_lane:
            # Plain DNA-row click. The cursor is already at `bp` (set
            # during mouse_down before this click event fired); we just
            # clear any lingering whole-feature highlight from a prior
            # lane click, leave map/sidebar selection alone, and let the
            # caller's mouse-down refresh stand.
            if seq_pnl._sel_range is not None:
                seq_pnl._sel_range = None
                seq_pnl._refresh_view()
            return

        # Route through PlasmidMap._bp_in + _feat_len so wrap features
        # (end < start) are honoured. A naive `s <= bp < e` misses every
        # origin-spanning feature; `e - s` is also negative for them.
        total = len(seq_pnl._seq)
        best_idx  = -1
        best_span = float("inf")
        for i, f in enumerate(pm._feats):
            if not pm._bp_in(bp, f):
                continue
            span = _feat_len(f["start"], f["end"], total) if total else 0
            if span < best_span:
                best_span = span
                best_idx  = i
        if best_idx >= 0:
            f = pm._feats[best_idx]
            pm.select_feature(best_idx)
            sidebar.show_detail(f)
            sidebar.highlight_row(best_idx)
            self._focus_feature(f, bp)

    @on(PlasmidMap.FeatureSelected)
    def _map_feat_selected(self, event: PlasmidMap.FeatureSelected):
        sidebar = self.query_one("#sidebar",   FeatureSidebar)
        sidebar.show_detail(event.feat_dict)
        if event.idx >= 0:
            sidebar.highlight_row(event.idx)
        self._focus_feature(event.feat_dict, event.bp)

    @on(FeatureSidebar.RowActivated)
    def _sidebar_row_activated(self, event: FeatureSidebar.RowActivated):
        pm      = self.query_one("#plasmid-map", PlasmidMap)
        sidebar = self.query_one("#sidebar",     FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel",   SequencePanel)
        pm.select_feature(event.idx)
        f = pm._feats[event.idx] if 0 <= event.idx < len(pm._feats) else None
        sidebar.show_detail(f)
        # Sidebar has no click-bp of its own — anchor at the feature's
        # wrap-aware midpoint so the cursor lands inside the feature and
        # the centred view shows it symmetrically.
        bp = self._wrap_aware_midpoint(f, len(seq_pnl._seq)) if f is not None else -1
        self._focus_feature(f, bp)

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
                ("Open file (.gb / .dna)  [^O]", "open_file"),
                ("Fetch from NCBI  [f]",         "fetch"),
                ("---",                          None),
                ("Add to Library  [^⇧A]",        "add_to_library"),
                ("Save  [^S]",                   "save"),
                ("Export as GenBank (.gb)...",   "export_genbank"),
                ("---",                          None),
                ("Collections...",               "open_collections"),
                ("---",                          None),
                ("Quit  [q]",                    "quit"),
            ],
            "Edit": [
                ("Edit Sequence  [^E]",            "edit_seq"),
                ("---",                             None),
                ("Undo",                            "undo"),
                ("Redo",                            "redo"),
                ("---",                             None),
                ("Add Feature...  [^F]",            "add_feature"),
                ("Capture selection → feat-lib  [^⇧F]", "capture_to_features"),
                ("Delete Feature",                  "delete_feature"),
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

    def action_capture_to_features(self) -> None:
        """Ctrl+Shift+F: grab the drag-selected DNA *or* the highlighted feature
        from the main view, open the AddFeatureModal prefilled, and after
        Save transport the user to the FeatureLibraryScreen so they see the
        new entry in context.

        Priority:
          1. ``sp._user_sel`` (Shift+drag / Shift+arrow region / sidebar-click
             selection). If the range exactly matches one of ``pm._feats``,
             the prefill carries that feature's metadata (type, strand,
             color, qualifiers). Otherwise it's a raw DNA capture with
             generic ``misc_feature`` / strand=1 defaults.
          2. ``pm.selected_idx`` (feature clicked on map / sidebar with no
             active selection) — always full metadata.

        Insert-at-cursor is intentionally disabled for the capture workflow
        because the bases are already in the record."""
        try:
            sp = self.query_one("#seq-panel",   SequencePanel)
            pm = self.query_one("#plasmid-map", PlasmidMap)
        except NoMatches:
            return
        seq = sp._seq or ""
        if not seq:
            self.notify("No sequence loaded.", severity="warning")
            return

        prefill: "dict | None" = None

        if sp._user_sel is not None:
            s, e = sp._user_sel
            sub = seq[s:e].upper()
            if not sub:
                self.notify("Selection is empty.", severity="warning")
                return
            matched = -1
            for i, f in enumerate(pm._feats):
                if f.get("type") in ("resite", "recut"):
                    continue
                if (f.get("start"), f.get("end")) == (s, e):
                    matched = i
                    break
            if matched >= 0:
                prefill = self._prefill_from_feature(pm._feats[matched], seq)
            else:
                prefill = {
                    "name":         "",
                    "feature_type": "misc_feature",
                    "sequence":     sub,
                    "strand":       1,
                    "color":        None,
                    "qualifiers":   {},
                    "description":  f"Captured from bases {s + 1}..{e}",
                }
        elif 0 <= pm.selected_idx < len(pm._feats):
            feat = pm._feats[pm.selected_idx]
            if feat.get("type") in ("resite", "recut"):
                self.notify(
                    "Select a real feature (not a restriction site).",
                    severity="warning",
                )
                return
            prefill = self._prefill_from_feature(feat, seq)
        else:
            self.notify(
                "Select a feature (click the map) or Shift+drag a DNA region first.",
                severity="information",
            )
            return

        self.push_screen(
            AddFeatureModal(prefill=prefill, have_cursor=False),
            callback=self._capture_feature_result,
        )

    def _prefill_from_feature(self, feat: dict, seq: str) -> dict:
        """Build an AddFeatureModal prefill from a ``pm._feats`` entry. Pulls
        qualifiers from the matching SeqFeature on the current record so
        gene/product/note values ride along with the capture."""
        s, e = feat["start"], feat["end"]
        fwd_slice = (seq[s:] + seq[:e]) if (e < s) else seq[s:e]
        strand = feat.get("strand", 1)
        feat_seq = _rc(fwd_slice) if strand == -1 else fwd_slice
        quals: dict = {}
        if self._current_record is not None:
            for bf in getattr(self._current_record, "features", []) or []:
                if bf.type != feat.get("type"):
                    continue
                try:
                    bs = int(bf.location.start)
                    be = int(bf.location.end)
                except (TypeError, ValueError):
                    continue
                if (bs, be) in ((s, e), (e, s)):
                    quals = {k: list(v) if isinstance(v, (list, tuple)) else [v]
                             for k, v in (bf.qualifiers or {}).items()}
                    break
        # The map palette uses Rich's ``color(N)`` syntax which round-trips
        # through Style.parse fine, but chokes Rich's markup lexer when
        # rewrapped as ``[color(N)]...[/]``. Convert to canonical hex on
        # capture so the stored library entry and every downstream preview
        # can use simple markup without escaping gymnastics.
        raw_color = feat.get("color")
        if isinstance(raw_color, str) and raw_color.startswith("color("):
            raw_color = _normalise_color_input(raw_color) or raw_color
        return {
            "name":         feat.get("label") or feat.get("type", ""),
            "feature_type": feat.get("type", "misc_feature"),
            "sequence":     feat_seq.upper(),
            "strand":       strand,
            "color":        raw_color,
            "qualifiers":   quals,
            "description":  "",
        }

    def _capture_feature_result(self, result) -> None:
        """Callback for the Ctrl+Shift+F capture flow. On Save, persist and push
        FeatureLibraryScreen. Insert is unreachable (have_cursor=False)."""
        if not result:
            return
        action = result.get("action")
        entry  = result.get("entry") or {}
        if action != "save":
            return
        if self._persist_feature_entry(entry):
            self.notify(f"Added '{entry.get('name')}' to feature library.")
            self.push_screen(FeatureLibraryScreen())

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

    def _persist_feature_entry(self, entry: dict) -> bool:
        """Append `entry` to features.json, de-duping on (name, feature_type)
        so the latest write wins. Returns True on success; on failure, logs
        and notifies the user so callers can bail cleanly."""
        entries = _load_features()
        key = (entry.get("name"), entry.get("feature_type"))
        entries = [e for e in entries
                   if (e.get("name"), e.get("feature_type")) != key]
        entries.append(entry)
        try:
            _save_features(entries)
        except (OSError, ValueError) as exc:
            _log.exception("Failed to save feature to library")
            self.notify(f"Save failed: {exc}", severity="error")
            return False
        return True

    def _add_feature_result(self, result) -> None:
        """Callback for AddFeatureModal. `result` is either None (cancel) or
        ``{"action": "save"|"insert", "entry": {...}}``."""
        if not result:
            return
        action = result.get("action")
        entry  = result.get("entry") or {}
        if action == "save":
            if self._persist_feature_entry(entry):
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

    def action_open_parts_bin(self) -> None:
        self.push_screen(PartsBinModal())

    def action_open_collections(self) -> None:
        """Open the plasmid-collections manager."""
        self.push_screen(CollectionsModal(),
                         callback=self._on_collections_dismissed)

    def _on_collections_dismissed(self, result) -> None:
        """If a collection was loaded, replace the current library wholesale
        and refresh the panel into the plasmids view of the new collection.
        Keep the active record (the user might be editing) — they can
        switch via the panel if needed."""
        if not isinstance(result, dict) or not result.get("loaded"):
            return
        try:
            panel = self.query_one("#library", LibraryPanel)
            panel._view_mode = "plasmids"
            panel._apply_view_mode()
            panel._repopulate()
        except NoMatches:
            pass
        self.notify(
            f"Loaded collection '{result['loaded']}' "
            f"({result['n_plasmids']} plasmid(s)). "
            f"Previous library backed up to plasmid_library.json.bak.",
            timeout=10,
        )

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

        Always opens — the modal supports plasmid library, parts bin, and
        protein-input sources, so a CDS source is reachable even with no
        plasmid loaded. The 'Current map features' option just shows
        '(no CDS features on this plasmid)' if the canvas is empty."""
        rec = self._current_record
        seq = str(rec.seq) if rec is not None else ""
        name = (rec.name or "") if rec is not None else ""
        feats: list = []
        try:
            feats = self.query_one("#plasmid-map", PlasmidMap)._feats
        except NoMatches:
            pass
        self.push_screen(MutagenizeModal(seq, feats, name))

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
