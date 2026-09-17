"""splicecraft_persistence — the domain-agnostic save/load engine (layer 1).

Extracted from the splicecraft.py hub (Phase B-main). This is the persistence
*service*: atomic write + multi-generation backup rotation + lost-entries
spillover + schema-envelope parse/migrate + the read-back-validated safe save,
plus the save-authorization chokepoint that refuses writes to the real data dir
unless a sanctioned caller opted in. Domain `_save_X` / `_load_X` wrappers stay
in the hub and call this engine; they migrate to their own modules later.

Dependencies are deliberately minimal so this layer never imports the hub:
stdlib + `splicecraft_state` (paths, caches, the `_SAVES_AUTHORIZED` flag, the
backup/lost-entries tunables, the blob-dehydration hooks) + the logger. The hub
re-exports every name below (``from splicecraft_persistence import _safe_save_json
as _safe_save_json`` ...) so ``sc._safe_save_json`` and every existing call site
resolve unchanged.

SACRED (see CLAUDE.md): `_safe_save_json` always backs up (.bak + mkstemp +
fsync + os.replace), re-raises on failure, and gates on
`_state._SAVES_AUTHORIZED` via `_refuse_unauthorized_write` at its first line.
Blob dehydration is invoked through `_state._dehydrate_*_hook` (registered by the
hub) so the engine stays decoupled from the blob subsystem.
"""
from __future__ import annotations

import json
import threading
from typing import Callable as _Callable
from contextlib import contextmanager
from datetime import datetime as _datetime
from pathlib import Path

import splicecraft_state as _state
from splicecraft_logging import _log, _log_event


_ENTRY_MIGRATIONS: "dict[str, dict[tuple[int, int], _Callable[[dict], dict]]]" = {
    # No migrations registered yet — schema is at v1 across every file.
    # Populate this when the schema bumps.
}

def _migrate_entries(entries: list,
                     from_version: int,
                     to_version: int,
                     label: str) -> "tuple[list, list[str]]":
    """Walk registered migrations from `from_version` up to
    `to_version` for the given file `label`. Returns
    `(migrated_entries, warnings)`. Each missing intermediate step
    is a no-op (the entries pass through unchanged).

    Pure (no I/O), idempotent for any range with no registered
    migrations. The function never returns the input list directly
    — it always returns a fresh list (even if no transformation
    happened) so callers can assume independence from the cache."""
    warnings: list[str] = []
    if from_version >= to_version:
        return list(entries), warnings
    if not isinstance(entries, list):
        warnings.append(f"{label}: refusing to migrate non-list entries")
        return [], warnings
    out: list = list(entries)
    migrations = _ENTRY_MIGRATIONS.get(label, {})
    current = from_version
    while current < to_version:
        next_step = current + 1
        migrator = migrations.get((current, next_step))
        if migrator is None:
            # No registered migration for this step — assume the
            # bump was additive and entries pass through. Defensive
            # default-handling at field-read time covers missing keys.
            current = next_step
            continue
        _log_event(
            "migration.step", label=label,
            from_v=current, to_v=next_step,
            n_entries=len(out),
        )
        migrated: list = []
        n_failed = 0
        for entry in out:
            if not isinstance(entry, dict):
                # Non-dict entries can't be migrated; surface a
                # warning and drop them rather than crash the loader.
                warnings.append(
                    f"{label}: dropped non-dict entry during "
                    f"v{current} → v{next_step} migration: {type(entry).__name__}"
                )
                n_failed += 1
                continue
            try:
                # Migrators MUST return a fresh dict; we don't deep-
                # copy here (the cache layer above does that on its
                # own read/save path).
                migrated.append(migrator(entry))
            except (KeyError, ValueError, TypeError) as exc:
                warnings.append(
                    f"{label}: v{current} → v{next_step} migration "
                    f"failed for entry: {exc}"
                )
                _log_event(
                    "migration.failed", label=label,
                    from_v=current, to_v=next_step,
                    exc_type=type(exc).__name__,
                )
                n_failed += 1
                # Keep the entry as-is rather than drop it; better to
                # surface a v1-shaped entry in a v2 list than to lose
                # the user's data outright.
                migrated.append(entry)
        out = migrated
        current = next_step
        if n_failed > 0:
            _log_event(
                "migration.step.done", label=label,
                to_v=next_step, n_entries=len(out),
                n_failed=n_failed,
            )
    return out, warnings

def _extract_entries(raw, label: str) -> "tuple[list | None, str | None]":
    """Return (entries, warning) from a parsed-JSON payload.

    Accepts both the envelope format `{"_schema_version": N, "entries": [...]}`
    and the legacy bare-list format. Returns (None, warning) on unknown shape
    so the caller can fall through to the .bak.

    Forward-compat: a higher schema_version loads with a warning (entries
    flow through unchanged; deepcopy preserves unknown fields across a
    save round-trip). The observed version is recorded so the next save
    for the same label preserves it on disk.
    Backward-compat: a lower schema_version is funnelled through the
    `_migrate_entries` registry so every consumer sees the current shape.
    """
    if isinstance(raw, list):
        # Legacy bare-list format predates `_schema_version` so the
        # entries are version-0 by definition. Run them through any
        # registered 0→1 migration before handing back.
        if 0 < _state._CURRENT_SCHEMA_VERSION:
            migrated, _ = _migrate_entries(raw, 0, _state._CURRENT_SCHEMA_VERSION, label)
            return migrated, None
        return list(raw), None
    if isinstance(raw, dict) and isinstance(raw.get("entries"), list):
        version = raw.get("_schema_version")
        entries_raw = list(raw["entries"])
        if version is not None and isinstance(version, int) and \
                version > _state._CURRENT_SCHEMA_VERSION:
            # Written by a newer SpliceCraft. Load the entries but warn
            # so the user knows fields may be silently dropped on
            # re-save. The version recording happens in
            # `_safe_load_json` where the path is available — keying
            # the registry by path is more accurate than by label.
            return entries_raw, (
                f"{label} was written by a newer SpliceCraft "
                f"(schema v{version} > v{_state._CURRENT_SCHEMA_VERSION}) — some "
                f"fields may be lost on save."
            )
        # Apply any registered migrations to bring older payloads up
        # to the current schema. v=None falls back to 1 (the first
        # explicit version) — older bare-list files never had it.
        from_v = int(version) if isinstance(version, int) else 1
        if from_v < _state._CURRENT_SCHEMA_VERSION:
            migrated, warns = _migrate_entries(
                entries_raw, from_v, _state._CURRENT_SCHEMA_VERSION, label
            )
            warning = warns[0] if warns else None
            return migrated, warning
        return entries_raw, None
    return None, f"{label}: unexpected JSON shape ({type(raw).__name__})"

def _safe_file_size_check(path: Path, max_bytes: int, label: str
                            ) -> "tuple[bool, str | None]":
    """Verify that `path` is a regular file (not symlink / FIFO /
    device) AND that its size is ≤ `max_bytes`. Returns
    `(ok, reason_if_not)`.

    `path.stat()` follows symlinks; a symlink → `/dev/zero` reports
    `st_size=0`, passes the byte cap, then a subsequent read consumes
    RAM until the OS kills the worker. `lstat` + `S_ISREG` rejects
    those up front. FIFOs / character devices / sockets also report
    `st_size=0` and are similarly hostile; the same `S_ISREG` check
    catches them.

    Returns `(False, reason)` on any of:
      * path doesn't exist
      * path is a symlink (regardless of target)
      * path is not a regular file (FIFO / device / socket / dir)
      * file size > max_bytes
    """
    import os, stat
    try:
        st = os.lstat(str(path))
    except OSError as exc:
        return False, f"{label} could not stat: {exc}"
    mode = st.st_mode
    if stat.S_ISLNK(mode):
        return False, (
            f"{label} is a symlink — refusing for safety "
            f"(symlinks to character devices report 0 bytes and OOM "
            f"on read)"
        )
    if not stat.S_ISREG(mode):
        return False, f"{label} is not a regular file (mode {oct(mode)})"
    if st.st_size > max_bytes:
        return False, (
            f"{label} file is {st.st_size:,} bytes (cap "
            f"{max_bytes:,}); refusing to load"
        )
    return True, None

def _fsync_parent_dir(path: Path) -> None:
    """Fsync `path.parent` so the rename's directory entry update is
    journalled. `os.replace` is atomic for the *inode* on POSIX, but
    the directory entry change is journalled separately — a power loss
    between the rename and the parent-dir flush can leave the inode
    containing the new data while the directory entry still points at
    the OLD inode after fsck. Linux-only (Windows doesn't expose
    directory fsync; opening a directory for fsync fails). Best-effort:
    a failure is logged at WARNING but NOT re-raised — the rename has
    already succeeded; only crash-durability of the dir entry is at
    stake, so a recurring fault should be visible, not silent."""
    import os
    try:
        dir_fd = os.open(str(path.parent), os.O_RDONLY)
    except OSError:
        return
    try:
        try:
            os.fsync(dir_fd)
        except OSError as exc:
            # Not re-raised (the rename already succeeded); logged so a
            # recurring directory-durability fault surfaces in the
            # diagnostic bundle instead of vanishing (2026-06-12 audit).
            _log.warning(
                "fsync of parent dir %s failed (rename succeeded; dir-entry "
                "durability is best-effort): %s", path.parent, exc,
            )
    finally:
        try:
            os.close(dir_fd)
        except OSError:
            pass

def _atomic_write_text(path: Path, text: str, encoding: str = "utf-8",
                       newline: "str | None" = None) -> None:
    """Atomically write *text* to *path* via ``tempfile`` + ``os.replace``.

    ``newline`` is handed straight to :func:`open`. The default (``None``)
    keeps Python's universal-newline translation, which turns every ``\n``
    into ``os.linesep`` — right for a platform-native text file, WRONG for an
    interchange format whose line endings are part of the contract. Two
    callers pass it explicitly:

      * the sequence-format exporters (``.gb`` / ``.embl`` / ``.gff3`` /
        ``.fa``) pass ``"\n"``, so a GenBank file written on Windows is
        byte-identical to the same file written on Linux;
      * the CSV writer passes ``""``, because :mod:`csv` already emits the
        RFC 4180 ``\r\n`` itself — translation on top of that produced
        ``\r\r\n``, which spreadsheet importers read as blank rows.

    Guarantees: a concurrent crash leaves either the previous file intact
    or the new file in place — never a partial write. The parent
    directory is also fsynced after the rename so the directory entry
    change is journalled (see `_fsync_parent_dir`). Callers that need
    a ``.bak`` should use :func:`_safe_save_json` instead (it layers the
    envelope, shrink-guard, and schema handling on top of this).

    Refuses to write when *path* is a symlink — ``os.replace`` would
    silently replace the symlink itself with a regular file, breaking
    the link chain and leaving the original target stale. Users who
    legitimately want to overwrite through a symlink should resolve
    the path themselves before calling.
    """
    import os
    import stat as _stat
    import tempfile
    try:
        st = path.lstat()
    except FileNotFoundError:
        st = None
    if st is not None and _stat.S_ISLNK(st.st_mode):
        raise OSError(
            f"refusing to write to {path}: target is a symlink "
            f"(would break the link chain). Resolve the path and "
            f"retry if this is intentional."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
    )
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline=newline) as fh:
            fh.write(text)
            fh.flush()
            # 2026-05-27 (audit-2 H2): re-raise fsync failures. EIO /
            # ENOSPC on fsync means the data is NOT on stable storage;
            # silently proceeding to `os.replace` produces a rename
            # that points at an inode whose data never reached the
            # platters — power loss after that = lost save with UI
            # already showing "saved". Loud failure is correct here.
            os.fsync(fh.fileno())
        os.replace(tmp, str(path))
        _fsync_parent_dir(path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise

def _atomic_write_bytes(path: Path, data: bytes) -> None:
    """Atomically write *data* to *path* via ``tempfile`` + ``os.replace``.

    Byte-mode counterpart to :func:`_atomic_write_text`. Used by the
    `_safe_save_json` backup rotation (legacy `.bak` + timestamped
    `.bak.<ts>`) and the daily-snapshot copy so a mid-write crash
    cannot truncate the recovery files that invariant #31's four-layer
    safety net depends on. Raises ``OSError`` on disk failure so
    callers can decide to surface or log.

    Refuses symlinked targets — see `_atomic_write_text` for rationale.
    """
    import os
    import stat as _stat
    import tempfile
    try:
        st = path.lstat()
    except FileNotFoundError:
        st = None
    if st is not None and _stat.S_ISLNK(st.st_mode):
        raise OSError(
            f"refusing to write to {path}: target is a symlink "
            f"(would break the link chain). Resolve the path and "
            f"retry if this is intentional."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
    )
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
            fh.flush()
            # 2026-05-27 (audit-2 H2): re-raise fsync — see the text-
            # writer twin above for rationale (silent EIO = lost save).
            os.fsync(fh.fileno())
        os.replace(tmp, str(path))
        _fsync_parent_dir(path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise

def _backup_filename_patterns(path: Path) -> "tuple[str, ...]":
    """Glob patterns matching every timestamped backup of `path`.

    Two patterns are returned and the caller is expected to union them:
      * ``<name>.bak.????????-??????`` — the base timestamp form
      * ``<name>.bak.????????-??????.*`` — collision-bumped variants
        (``.bak.20260515-123456.1``, ``.2``, …) which `_safe_save_json`
        emits when two saves land in the same wall-second

    Pre-0.8.9 returned only the base pattern, so collision-bumped
    backups leaked forever — a rapid-Ctrl+S burst on a slow disk would
    accumulate `.bak.<ts>.<N>` files that the pruner never matched.
    """
    base = f"{path.name}.bak.????????-??????"
    return (base, base + ".*")

def _backup_filename_pattern(path: Path) -> str:
    """Deprecated single-pattern accessor (returned only the base
    timestamp form). Kept for any external probe that may have
    referenced it; new code should use `_backup_filename_patterns`."""
    return _backup_filename_patterns(path)[0]

def _backup_sort_key(name: str) -> "tuple":
    """Oldest→newest ordering key for one backup filename.

    The collision bump MUST be compared numerically. As a bare integer in
    the name it sorts lexicographically — `.1 < .10 < .11 < .12 < .2` — so
    once a file was saved more than nine times inside one wall-second, the
    "newest" backup was bump 9: retention then pruned the genuinely newest
    generations and the recovery chain restored a stale one. Reachable from
    a bulk delete, a scripted agent loop, or repeated settings writes.
    """
    rest = name.split(".bak.", 1)[-1]
    if rest.endswith(".gz"):
        rest = rest[:-3]
    parts = rest.split(".")
    ts = parts[0]
    bump = 0
    if len(parts) > 1:
        try:
            bump = int(parts[1])
        except (TypeError, ValueError):
            bump = 0
    return (ts, bump, name)


def _iter_backups(path: Path) -> "list[Path]":
    """Return all timestamped backups for `path` across both glob
    patterns, de-duplicated, OLDEST FIRST — so `[-1]` is the newest and
    `reversed(...)` walks newest→oldest, which is what every consumer
    (retention, compression, the recovery chain, the save-time dedup)
    assumes. Ordering comes from `_backup_sort_key`, never from the raw
    filename."""
    seen: "dict[str, Path]" = {}
    for pat in _backup_filename_patterns(path):
        try:
            for p in path.parent.glob(pat):
                seen[p.name] = p
        except OSError:
            continue
    return [seen[k] for k in sorted(seen.keys(), key=_backup_sort_key)]

def _prune_backups(path: Path, keep: "int | None" = None) -> None:
    """Delete all but the most recent `keep` timestamped backups of
    `path`, then enforce an aggregate byte ceiling.

    Two-stage retention:
      1. COUNT — keep the newest `keep` (default `_state._BACKUP_RETENTION_COUNT`)
         timestamped backups; unlink the rest.
      2. SIZE — if the survivors still total more than
         `_state._BACKUP_TOTAL_SIZE_CAP_BYTES`, drop the oldest remaining ones
         until the total falls under the cap, but never below
         `_state._BACKUP_MIN_KEEP` newest generations.

    The plain `<file>.bak` (legacy single-generation) is never pruned by
    this function — it's overwritten on every save. `keep` and both caps
    are read at call time (not bound as default arguments) so a test can
    monkeypatch the module constants and observe the new values."""
    if keep is None:
        keep = _state._BACKUP_RETENTION_COUNT
    # A retention of 0 would delete EVERY rotated backup (and a negative one
    # all but the last), silently turning the rotation off. Always keep at
    # least the newest generation.
    keep = max(1, int(keep))
    candidates = _iter_backups(path)
    # Sort descending (newest first by lex-sortable timestamp), keep
    # the head, prune the tail.
    candidates.reverse()
    for old in candidates[keep:]:
        try:
            old.unlink()
        except OSError:
            _log.debug("Could not prune old backup %s", old)
    # Stage 2: aggregate byte cap. Re-scan survivors (the count-prune
    # above may have unlinked some) newest-first, summing sizes; once
    # over the cap, drop oldest until under it. Always retain at least
    # `_state._BACKUP_MIN_KEEP` newest so a single oversized file can't leave
    # the user with no rollback point.
    survivors = _iter_backups(path)
    survivors.reverse()   # newest first
    sized: "list[tuple[Path, int]]" = []
    total = 0
    for b in survivors:
        try:
            sz = b.stat().st_size
        except OSError:
            sz = 0
        sized.append((b, sz))
        total += sz
    cap = _state._BACKUP_TOTAL_SIZE_CAP_BYTES
    idx = len(sized) - 1
    while total > cap and idx >= _state._BACKUP_MIN_KEEP:
        victim, sz = sized[idx]
        try:
            victim.unlink()
            total -= sz
            _log.info("Byte-cap pruned old backup %s (%d bytes)", victim, sz)
        except OSError:
            _log.debug("Could not byte-prune old backup %s", victim)
        idx -= 1

def _read_backup_bytes(bak: Path) -> bytes:
    """Read a backup's raw JSON bytes, transparently decompressing a
    gzipped backup (`.bak.<ts>.gz`). `_compress_old_backups` gzips OLDER
    timestamped backups to reclaim disk; the newest timestamped backup
    and the legacy `.bak` stay PLAIN, so the primary recovery path never
    depends on gzip. Raises OSError on a read / decompress failure so
    callers can fall through to the next backup."""
    if bak.name.endswith(".gz"):
        import gzip
        with gzip.open(bak, "rb") as fh:
            return fh.read()
    return bak.read_bytes()

def _compress_old_backups(path: Path) -> None:
    """Gzip every timestamped backup of `path` EXCEPT the newest, to
    reclaim disk — JSON of plasmid `gb_text` compresses ~4-5× (the user's
    backups were 4.7 GB uncompressed across the two big files). The
    newest timestamped backup and the legacy `.bak` are deliberately left
    PLAIN so the primary recovery path (main → legacy `.bak` → newest
    timestamped) never has to decompress; compression only affects
    tertiary, older recovery points. Idempotent and best-effort:
    already-`.gz` backups are skipped, failures logged and skipped, never
    raised. Runs off the interactive path (startup housekeeping) so the
    CPU of compressing a 274 MB backup never adds latency to a save."""
    import gzip
    baks = _iter_backups(path)   # oldest → newest (lex sort)
    if len(baks) <= 1:
        return
    for bak in baks[:-1]:        # leave the newest plain
        if bak.name.endswith(".gz"):
            continue
        gz_path = bak.with_name(bak.name + ".gz")
        if gz_path.exists():
            # A prior run wrote the .gz but crashed before unlinking the
            # plain. Drop the redundant plain copy now.
            try:
                bak.unlink()
            except OSError:
                pass
            continue
        try:
            raw = bak.read_bytes()
            # Atomic .gz write (tempfile + replace) so a crash mid-compress
            # can't leave a truncated .gz masquerading as a good backup.
            buf = gzip.compress(raw, compresslevel=6)
            _atomic_write_bytes(gz_path, buf)
        except OSError:
            _log.debug("Could not compress backup %s", bak)
            continue
        # Only remove the plain backup once the .gz is safely on disk.
        try:
            bak.unlink()
            _log.info("Compressed backup %s → %s (%d → %d bytes)",
                      bak.name, gz_path.name, len(raw), len(buf))
        except OSError:
            _log.debug("Could not remove plain backup %s after compress", bak)

def _spill_lost_entries(path: Path, lost: list, label: str) -> "Path | None":
    """Persist `lost` (entries that the new save would discard) under
    `<DATA_DIR>/lost_entries/<file_stem>-<ts>.json`. Returns the
    written path or None on failure. Never raises — this is the
    last-resort safety net and must not interfere with the parent
    save call."""
    if not lost:
        return None
    try:
        lost_dir = path.parent / _state._LOST_ENTRIES_DIR_NAME
        lost_dir.mkdir(parents=True, exist_ok=True)
        ts = _datetime.now().strftime("%Y%m%d-%H%M%S")
        out = lost_dir / f"{path.stem}-{ts}.json"
        # Bump on collision — two saves in the same wall-second both
        # tripping the shrink guard would otherwise have the second
        # spill silently overwrite the first. Spillover is the
        # last-ditch recovery; losing it is a hard regression.
        bump = 0
        while out.exists():
            bump += 1
            out = lost_dir / f"{path.stem}-{ts}.{bump}.json"
        # Route through `_atomic_write_text` so a mid-write crash
        # (disk full, RO mount, power loss) leaves either nothing or
        # a complete recovery dump — never a half-written file that
        # masquerades as evidence. This is the safety-net for the
        # safety-net, but a corrupt spill is worse than no spill.
        _atomic_write_text(
            out,
            json.dumps({
                "_schema_version": _state._CURRENT_SCHEMA_VERSION,
                "_label":          label,
                "_recovered_from": str(path),
                "_recovered_at":   ts,
                "entries":         lost,
            }, indent=2),
        )
        # Bound the directory right after writing so a long-running
        # session that spills repeatedly can't grow it without limit.
        # The just-written `out` is newest, so it's never the prune
        # victim. Best-effort — `_prune_lost_entries` never raises.
        _prune_lost_entries(lost_dir)
        return out
    except Exception:
        _log.exception("Could not spill lost entries for %s", label)
        return None

def _spill_raw_bytes(path: Path, data: bytes, label: str,
                      *, reason: str) -> "Path | None":
    """Sweep #35 (2026-05-26): unstructured byte-level spill for cases
    where `_extract_entries` can't parse the prior file (corrupt JSON,
    unknown schema). Drops the raw bytes into
    `<DATA_DIR>/lost_entries/<file_stem>-raw-<ts>.<ext>` so a user
    doing forensic recovery can `cat` or `grep` the content even
    when the timestamped `.bak.<ts>` rotation has been pruned.
    Never raises — best-effort safety net for the safety net."""
    if not data:
        return None
    try:
        lost_dir = path.parent / _state._LOST_ENTRIES_DIR_NAME
        lost_dir.mkdir(parents=True, exist_ok=True)
        ts = _datetime.now().strftime("%Y%m%d-%H%M%S")
        out = lost_dir / f"{path.stem}-raw-{ts}{path.suffix}"
        bump = 0
        while out.exists():
            bump += 1
            out = lost_dir / f"{path.stem}-raw-{ts}.{bump}{path.suffix}"
        _atomic_write_bytes(out, data)
        _log.info(
            "%s: raw-bytes spill (%s) → %s (%d bytes)",
            label, reason, out, len(data),
        )
        _prune_lost_entries(lost_dir)
        return out
    except Exception:
        _log.exception("Could not raw-spill prior bytes for %s", label)
        return None

def _diff_lost_entries(prev_entries: list, new_entries: list) -> list:
    """Return entries that exist in `prev_entries` but not in
    `new_entries`. Match by `id` first (every plasmid library entry
    carries one); fall back to identity equality for entries that
    don't have an id."""
    if not prev_entries:
        return []
    new_ids = {e.get("id") for e in new_entries
               if isinstance(e, dict) and e.get("id")}
    new_no_id = [e for e in new_entries
                 if isinstance(e, dict) and not e.get("id")]
    lost = []
    for e in prev_entries:
        if not isinstance(e, dict):
            # Pass through non-dict legacy entries — better to over-
            # spill than under-spill here.
            lost.append(e)
            continue
        eid = e.get("id")
        if eid:
            if eid not in new_ids:
                lost.append(e)
        else:
            # Best-effort identity match for id-less entries.
            if e not in new_no_id:
                lost.append(e)
    return lost

def _prune_lost_entries(lost_dir: "Path | None" = None) -> None:
    """Bound the `lost_entries/` spillover directory.

    `_spill_lost_entries` / `_spill_raw_bytes` drop recovery copies here.
    Pre-1.0.22 nothing pruned them, so collection switching (each switch
    shrinks the active library, tripping the suspicious-shrink spill)
    accumulated a full ~150 MB dump per switch — 1.5 GB of residue whose
    contents were never at risk (they live in collections.json). Two
    stages mirror `_prune_backups`:

      1. COUNT — keep the newest `_state._LOST_ENTRIES_RETENTION_COUNT` files.
      2. SIZE — drop oldest until under `_state._LOST_ENTRIES_TOTAL_SIZE_CAP_BYTES`,
         but never below one file (the most-recent recovery point).

    Ranked by mtime: spill names mix `<stem>-<ts>.json` and
    `<stem>-raw-<ts>.<ext>`, so mtime is the simplest total order.
    Best-effort — never raises; a locked / read-only data dir just means
    residue lingers until the next launch. `lost_dir` defaults to
    `_state._DATA_DIR/lost_entries` (sandboxed in tests via conftest's
    `_state._DATA_DIR` monkeypatch); callers with a known file path pass
    `<path>.parent / _state._LOST_ENTRIES_DIR_NAME` explicitly."""
    if lost_dir is None:
        lost_dir = _state._DATA_DIR / _state._LOST_ENTRIES_DIR_NAME
    try:
        files = [p for p in lost_dir.iterdir() if p.is_file()]
    except OSError:
        return

    def _mtime(p: Path) -> float:
        try:
            return p.stat().st_mtime
        except OSError:
            return 0.0

    files.sort(key=_mtime, reverse=True)   # newest first
    keep = _state._LOST_ENTRIES_RETENTION_COUNT
    # Stage 1: count.
    for old in files[keep:]:
        try:
            old.unlink()
            _log.info("Pruned old lost-entries spill %s", old)
        except OSError:
            _log.debug("Could not prune lost-entries file %s", old)
    # Stage 2: aggregate byte cap over the survivors.
    survivors = files[:keep]
    sized: "list[tuple[Path, int]]" = []
    total = 0
    for p in survivors:
        try:
            sz = p.stat().st_size
        except OSError:
            sz = 0
        sized.append((p, sz))
        total += sz
    cap = _state._LOST_ENTRIES_TOTAL_SIZE_CAP_BYTES
    idx = len(sized) - 1
    while total > cap and idx >= 1:
        victim, sz = sized[idx]
        try:
            victim.unlink()
            total -= sz
            _log.info("Byte-cap pruned lost-entries spill %s (%d bytes)",
                      victim, sz)
        except OSError:
            _log.debug("Could not byte-prune lost-entries file %s", victim)
        idx -= 1

# L3 catastrophic-shrink token. The shrink guard refuses writes that
# would discard >90% of an existing-cached file's entries (base
# population ≥10) unless this token is positive. Legitimate
# catastrophic shrinks (`_restore_from_backup`, programmatic wipes)
# arm the token via `_allow_catastrophic_shrink()` for the duration
# of the save. Refcount semantics so nested re-entry is safe.
# Thread-confined (2026-06-12 data-loss fix): these bypass refcounts MUST
# be visible only to the thread that armed them. They were module-global
# ints, but `_safe_save_json` takes no lock of its own and the app issues
# concurrent saves from several threads (agent server + LibraryPanel disk-
# writer `@work` groups). A mirror-swap armed on one thread (a collection
# switch) then disarmed the catastrophic-shrink REFUSAL for an UNRELATED
# destructive save on another thread, so an accidental wipe coinciding with
# a collection switch could slip past the guard built for the 2026-05-22
# incident. `threading.local()` confines each arming to its arming thread;
# arming and the guarded `_safe_save_json` are always same-thread (the
# active-slot write runs synchronously on the caller's thread inside
# `_cache_lock`; only the non-shrinking collection mirror defers to a worker).
_shrink_guard_tokens = threading.local()

def _catastrophic_shrink_depth() -> int:
    """Per-thread arm count for the L3 catastrophic-shrink bypass."""
    return getattr(_shrink_guard_tokens, "catastrophic", 0)

def _mirror_swap_depth() -> int:
    """Per-thread arm count for the expected-mirror-swap bypass."""
    return getattr(_shrink_guard_tokens, "mirror_swap", 0)

@contextmanager
def _expected_mirror_swap():
    """Arm the thread-local mirror-swap refcount for the `with` block.

    Wrap any active-slot write whose data is mirrored in a sibling
    collections file — the four mirror pairs in `_safe_save_json_mirror`,
    plus the active-collection switch via `_switch_active_collection_library`.
    The shrink guard then skips BOTH the redundant `lost_entries/` spill
    and the >90% catastrophic refusal, because the "dropped" entries are
    provably safe in the collections file. Refcounted; the decrement runs
    in `finally` so an exception can't leave the gate stuck open."""
    _shrink_guard_tokens.mirror_swap = _mirror_swap_depth() + 1
    try:
        yield
    finally:
        _shrink_guard_tokens.mirror_swap = _mirror_swap_depth() - 1

@contextmanager
def _allow_catastrophic_shrink():
    """Context manager that arms the L3 shrink guard's bypass token
    for the duration of the `with` block. Legitimate large-shrink
    paths (Restore-from-backup writing a tiny backup over a large
    live file; programmatic data wipes that route through
    `_safe_save_json` rather than `unlink`) wrap their save call:

        with _allow_catastrophic_shrink():
            _save_library(small_list)

    Refcounted so nested calls (e.g. a Restore that triggers a
    cache-mirror save during its own write) compose. The token
    decrement runs in `finally` so an exception inside the block
    can't leave the gate stuck open."""
    _shrink_guard_tokens.catastrophic = _catastrophic_shrink_depth() + 1
    try:
        yield
    finally:
        _shrink_guard_tokens.catastrophic = _catastrophic_shrink_depth() - 1

def _safe_save_json_mirror(path: Path, entries: list, label: str,
                            *, schema_version: "int | None" = None
                            ) -> None:
    """[INV-83, sweep #27] Mirror-write helper. Use this **instead of**
    a bare `_safe_save_json` for any write whose source-of-truth lives
    in a sibling *collections* file.

    Caught-live failure (2026-05-25 incident): switching from
    "DemoColl Parts" (26 entries) to an empty "FFE Parts" bin called
    ``_safe_save_json(_state._PARTS_BIN_FILE, [], "Parts bin")`` directly.
    The L3 catastrophic shrink guard saw 26 → 0 entries and refused
    the write — saving the user from data loss but also breaking the
    legitimate bin switch. The "old" data was NEVER at risk: the 26
    DemoColl entries live in ``parts_bin_collections.json`` under the
    "DemoColl Parts" key. The shrink guard correctly treats parts_bin.json
    as the primary store; this helper tells it explicitly that this
    particular write is a *mirror swap* — the primary is the collections
    file, the named bin's data is intact, so a 100% shrink is a
    legitimate state transition.

    Symmetric for every mirror pair:

      * `parts_bin.json`   <- `parts_bin_collections.json`     (active bin)
      * `plasmid_library.json` <- `collections.json`           (active collection)
      * `primers.json`     <- `primer_collections.json`        (active primer coll)
      * `experiments.json` <- `experiment_projects.json`       (active project)

    Implementation: wraps `_safe_save_json` in `_expected_mirror_swap()`
    (2026-06-03). Previously it used `_allow_catastrophic_shrink()`,
    which bypassed the >90% REFUSAL but STILL spilled the "dropped"
    entries to `lost_entries/` — a redundant copy of data that already
    lives in the collections file, and for the plasmid_library mirror a
    ~150 MB write on every collection switch (the cause of the 1.5 GB
    `lost_entries/` residue). `_expected_mirror_swap()` skips that
    redundant spill too. The wrapper exists so:
    (a) callers don't have to remember the context manager pattern;
    (b) a `grep _safe_save_json_mirror` enumerates every cross-mirror
    write — a future audit can confirm none have drifted back to the
    bare `_safe_save_json`."""
    with _expected_mirror_swap():
        _safe_save_json(path, entries, label, schema_version=schema_version)

def _refuse_unauthorized_write(path: Path, label: str) -> None:
    """Raise `RuntimeError` if the current process has not opted in
    to data-dir writes via `_authorize_writes` /
    `_authorize_writes_for_sandbox`. The message names the exact path
    that was about to be written so the caller can see what would
    have been clobbered.

    Routed through `_safe_save_json`'s first line so every persisted-
    file write (library / collections / primers / parts / features /
    grammars / entry vectors / codon tables / experiments / gels /
    protein motifs / settings) is gated. Pre-fix, an ad-hoc
    `import splicecraft; sc._save_collections([])` was enough to nuke
    160 MB of user data."""
    if _state._SAVES_AUTHORIZED:
        return
    raise RuntimeError(
        f"refusing to write {label!r} → {path}: data-dir writes are "
        f"not authorised in this process. If you're running a verifier "
        f"or probe, sandbox `XDG_DATA_HOME=$(mktemp -d)` BEFORE "
        f"`import splicecraft` and call "
        f"`splicecraft._authorize_writes_for_sandbox(splicecraft._DATA_DIR)`. "
        f"See CLAUDE.md sacred block + "
        f"`.claude/skills/verifier-splicecraft.md`. "
        f"(Authorisation is set automatically by `main()`, the pytest "
        f"`_protect_user_data` fixture, and the agent HTTP server.)"
    )

def _recovery_write_allowed(path: Path, label: str) -> bool:
    """May `_safe_load_json`'s recovery path WRITE to the data dir?

    Loading is supposed to be read-only, but the recovery chain mutates: it
    renames a corrupt live file aside and rewrites it from a backup, and it
    re-creates a missing main file. Those are data-dir writes, and they were
    the one hole in the L2 chokepoint — `_safe_save_json` refuses an
    unauthorised process on its first line, while a bare
    `import splicecraft; _load_library()` could still rename and overwrite
    live files. An unauthorised caller still gets its recovered entries back;
    it just does not get to touch the disk."""
    try:
        _refuse_unauthorized_write(path, label)
        return True
    except RuntimeError:
        _log.warning(
            "%s: recovery for %s would write to the data dir, but writes are "
            "not authorised in this process — returning the recovered entries "
            "WITHOUT repairing the file on disk.", label, path.name,
        )
        return False


def _refuse_unauthorized_delete(path: Path, label: str) -> None:
    """[INV-75, sweep #27] Mirror of `_refuse_unauthorized_write` for
    delete paths. Sweep #26 hardened the write chokepoint to cover
    every JSON save helper + the `.dna` sidecar, experiment images
    and crash-recovery autosaves; deletes are the equally-dangerous
    half of the data-dir API and were not gated.

    Threat model: a maintenance helper added in a future sweep that
    walks `_state._DATA_DIR` to prune rotated backups, expired snapshots, or
    stale crash-recovery files could be imported by an unsandboxed
    probe and tricked into running `unlink()` / `shutil.rmtree()`
    against the user's real files. The write gate alone wouldn't
    catch that — it only fires on `_safe_save_json` / the named
    chokepoint helpers.

    Sanctioned callers route through this helper before calling
    `unlink` / `rmtree` on a path under `_state._DATA_DIR`. The authorisation
    state is shared with the write gate (`_state._SAVES_AUTHORIZED`); a
    process that's authorised to write is by definition authorised
    to delete (the four sanctioned callers all need both)."""
    if _state._SAVES_AUTHORIZED:
        return
    raise RuntimeError(
        f"refusing to delete {label!r} → {path}: data-dir deletes are "
        f"not authorised in this process. Same gate as "
        f"`_refuse_unauthorized_write` — sandbox `XDG_DATA_HOME` or "
        f"call `_authorize_writes_for_sandbox` before invoking any "
        f"helper that unlinks/rmtrees under `_state._DATA_DIR`."
    )

def _safe_save_json(path: Path, entries: list, label: str,
                    schema_version: "int | None" = None) -> None:
    """Atomically write `entries` as JSON to `path`, backing up first.

    Writes the envelope format `{"_schema_version": N, "entries": [...]}`.

    Triple-layer data safety:

      1. **Atomic write** (tempfile + fsync + os.replace).
      2. **Multi-generation backup**: the previous file content is
         copied to `<file>.bak` (single, latest — legacy compat) AND
         to `<file>.bak.YYYYMMDD-HHMMSS` (timestamped). Last
         `_state._BACKUP_RETENTION_COUNT` (10) timestamped generations are
         retained on disk; older ones pruned.
      3. **Lost-entries spillover**: if the new write would shrink
         the entry count by >50% (with at least 5 prior entries), the
         discarded entries are dumped to
         `<DATA_DIR>/lost_entries/<file_stem>-<ts>.json` BEFORE the
         save proceeds. The save still runs (the user may have
         legitimately deleted entries), but the dropped data is never
         silently destroyed.

    **L2 write authorisation** (2026-05-22): first line refuses with
    `RuntimeError` if the process has not flipped `_state._SAVES_AUTHORIZED`
    via one of the four sanctioned callers (`main()`, pytest fixture,
    agent server, sandboxed verifier). Caught failure: a probe that
    `import`'d splicecraft from `/tmp/sc_probe.py` and called
    `_save_collections([])` nuked the user's real 160 MB collections
    file. The gate makes that scenario raise instead of write.

    Errors (disk full, RO mount, permission denied) are logged AND
    re-raised so callers can `notify` the user. Silent swallow used
    to desync UI state from disk — sacred invariant #7.

    Refuses to write through a symlink. Pre-0.8.9 the symlink check
    only lived at agent endpoints (`_check_agent_write_path`); the
    library / collections / parts-bin etc. save path itself trusted
    `path` to be a regular file. A symlink at the target (planted
    accidentally or otherwise) would let the backup-read step at the
    top of this function — `path.read_bytes()` follows symlinks — copy
    arbitrary file content into the user-readable `.bak`, and the
    subsequent atomic-write would overwrite the link target. Refuse
    up front so neither leak nor overwrite can happen.
    """
    _refuse_unauthorized_write(path, label)
    import os
    import tempfile

    if path.is_symlink():
        msg = (
            f"Refusing to save {label} through symlink at {path}. "
            f"Move/remove the symlink and rerun, or set the data "
            f"dir to a path with no symlinks in it."
        )
        _log.error(msg)
        raise OSError(msg)

    # Sweep #22: ancestor-chain symlink walk, mirroring the defense
    # `_check_agent_write_path` already had (sweep #10 invariant #50).
    # Pre-sweep, `_safe_save_json` only checked the path itself —
    # a symlink at any DEEPER ancestor (e.g. `~/.local` → `/etc`)
    # could redirect every write under the data dir. Walking the
    # ancestor chain via per-segment `is_symlink()` closes the gap.
    # Errors mid-walk surface as OSError so the caller's notify path
    # tells the user something concrete went wrong.
    cur = path.parent
    seen: set = set()
    while True:
        try:
            if cur.is_symlink():
                msg = (
                    f"Refusing to save {label}: ancestor directory "
                    f"is a symlink at {cur}. Move/remove the symlink "
                    f"and rerun."
                )
                _log.error(msg)
                raise OSError(msg)
        except OSError as exc:
            # Re-raise our own refusal as-is; wrap other stat errors
            # (permission, ENOENT mid-walk) as a refusal too — we
            # can't tell whether an ancestor is safe.
            if "Refusing to save" in str(exc):
                raise
            msg = (
                f"Refusing to save {label}: could not stat ancestor "
                f"{cur}: {exc}"
            )
            _log.error(msg)
            raise OSError(msg) from exc
        if cur.parent == cur or str(cur) in seen:
            break
        seen.add(str(cur))
        cur = cur.parent

    # Schema-version stamp: preserve the highest observed version for
    # this file. A file written by a newer SpliceCraft (recorded via
    # `_state._OBSERVED_SCHEMA_VERSIONS` from `_safe_load_json`) must stamp
    # back with that higher version so a v2 file edited under a v1
    # binary doesn't demote to v1 on save (which would re-trigger the
    # v1→v2 migrator next time a v2 binary opens it, double-migrating).
    if schema_version is None:
        observed = _state._OBSERVED_SCHEMA_VERSIONS.get(str(path), 0)
        schema_version = max(_state._CURRENT_SCHEMA_VERSION, observed)

    # OVERSIZE GUARD (sacred invariant — DO NOT BYPASS).
    #
    # If the file already on disk is larger than `_safe_load_json` can
    # read, the in-memory cache may be EMPTY (the load returned [] +
    # a warning). Overwriting the file here would silently destroy
    # whatever real data was in it — that's exactly how the user's
    # FFE collection got nuked on 2026-05-10 when their `gb_text`
    # bloat pushed `collections.json` past the old 50 MB cap.
    #
    # Refuse: raise OSError so the caller's `try/except` surfaces a
    # real error to the user instead of silently writing the empty
    # in-memory state. The cap was bumped to 1 GB; tripping this
    # guard now indicates a genuinely runaway file (or a future cap
    # bump is needed) — either way, "abort + scream" is the right
    # answer over "silently nuke".
    if path.exists():
        try:
            existing_size = path.stat().st_size
        except OSError:
            existing_size = 0
        if existing_size > _state._SAFE_LOAD_JSON_MAX_BYTES:
            msg = (
                f"Refusing to overwrite oversized {label} "
                f"({existing_size:,} bytes > "
                f"{_state._SAFE_LOAD_JSON_MAX_BYTES:,} cap). The previous "
                f"load was refused (see `_safe_load_json`); writing "
                f"the empty in-memory state would silently destroy "
                f"the file. Move {path} aside before retrying."
            )
            _log.error(msg)
            raise OSError(msg)

    # Blob-store dehydration (v1.0.23): for the two gb_text-heavy files,
    # replace each entry's inline `gb_text` with a content-addressed
    # `gb_ref` (the blob is written + VERIFIED first). Done here — at the
    # single write chokepoint — so EVERY library/collections save path is
    # covered, present and future (the in-place LibraryPanel workers,
    # Restore-from-backup, the active-collection mirror, agent endpoints).
    # A blob-write failure raises BEFORE any backup/metadata write, so the
    # previous-good file + its blobs stay intact (abort-don't-corrupt).
    # Entry COUNT is preserved, so the shrink guard / spillover below still
    # see true counts. `_state._LIBRARY_FILE` / `_state._COLLECTIONS_FILE` resolve as
    # module globals at call time (late binding) so conftest's path
    # monkeypatch applies in tests; non-blob files (primers, parts, …)
    # are untouched.
    #
    # The transform is invoked through a _state-registered hook
    # (`_state._dehydrate_{entries,collections}_hook`, wired by the hub at
    # import) instead of a direct call, so the persistence engine stays
    # decoupled from the blob-store subsystem (which remains hub-side): the
    # engine sibling imports only stdlib + _state + the logger, never the hub.
    # The hooks are always registered before any save can run (guarded by
    # test_blob_store); the `is not None` check only covers the import window,
    # where no save occurs. Path dispatch is unchanged from the direct form.
    if path == _state._LIBRARY_FILE:
        if _state._dehydrate_entries_hook is not None:
            entries = _state._dehydrate_entries_hook(entries)
    elif path == _state._COLLECTIONS_FILE:
        if _state._dehydrate_collections_hook is not None:
            entries = _state._dehydrate_collections_hook(entries)

    # Step 1: read prior content for backup + shrink-guard analysis.
    existing_count = 0
    prev_entries: "list | None" = None
    if path.exists():
        try:
            existing = path.read_bytes()
            if existing.strip():
                # Backup dedup (2026-06-03): if the prior file is
                # byte-identical to the most recent existing backup it is
                # already preserved, so BOTH backup writes below are
                # redundant. Fixes the "N identical 274 MB backups in one
                # wall-second" waste from rapid re-saves of unchanged data
                # (observed: 3 collections.json backups within 18 s).
                # Compare size first (cheap stat); read bytes only on a
                # size match. Dedups only against a PLAIN newest backup —
                # the newest is always plain (compression touches only
                # OLDER backups), so a `.gz` newest just means no dedup
                # this round (conservative, never wrong). When skipped the
                # legacy `.bak` already equals `existing` (it held the
                # prior save's content, which == `existing` precisely
                # because nothing changed), so recovery stays correct.
                skip_redundant_backup = False
                _existing_baks = _iter_backups(path)
                if _existing_baks:
                    _newest_bak = _existing_baks[-1]
                    try:
                        if (not _newest_bak.name.endswith(".gz")
                                and _newest_bak.stat().st_size == len(existing)
                                and _newest_bak.read_bytes() == existing):
                            skip_redundant_backup = True
                            _log.debug(
                                "Backup dedup: %s unchanged since %s — "
                                "skipping redundant backup writes",
                                path.name, _newest_bak.name,
                            )
                    except OSError:
                        skip_redundant_backup = False
                if not skip_redundant_backup:
                    # 2026-05-27 (audit-2 H1): write the timestamped
                    # backup FIRST so a verified-good rotated copy exists
                    # before we touch the legacy `.bak`. Pre-fix the
                    # order was legacy → timestamped: on ENOSPC during
                    # the legacy write the original-good `.bak` was
                    # clobbered with a truncated file, and the
                    # timestamped write that would have caught the
                    # condition came AFTER. Now: timestamped first
                    # (raises on failure → save aborts before legacy
                    # is touched), legacy second (also raises now to
                    # avoid the silent-corrupt-legacy-bak case).
                    ts = _datetime.now().strftime("%Y%m%d-%H%M%S")
                    bak_ts = path.with_name(f"{path.name}.bak.{ts}")
                    bump = 0
                    while bak_ts.exists():
                        bump += 1
                        # Zero-padded so the NAME also sorts correctly for
                        # any external tool reading the directory; the
                        # in-process ordering goes through `_backup_sort_key`,
                        # which parses the bump either way.
                        bak_ts = path.with_name(
                            f"{path.name}.bak.{ts}.{bump:03d}")
                    try:
                        _atomic_write_bytes(bak_ts, existing)
                    except OSError as exc:
                        _log.exception(
                            "Backup rotation failed for %s — aborting "
                            "save to preserve existing data", path,
                        )
                        raise OSError(
                            f"backup rotation failed for {label} "
                            f"({path.name}): {exc}. Save aborted to "
                            f"keep the previous-good file intact."
                        ) from exc
                    # Legacy single-generation backup — `_safe_load_json`'s
                    # recovery path still reads this exact name, and many
                    # tests depend on it. Keep it overwriting per save.
                    # 2026-05-27 (audit-2 H1): re-raise on failure too. A
                    # half-written legacy `.bak` is worse than no legacy
                    # `.bak` (recovery path tries it first); the
                    # timestamped backup written above is the safety
                    # net if this raises.
                    bak_legacy = path.with_suffix(path.suffix + ".bak")
                    try:
                        _atomic_write_bytes(bak_legacy, existing)
                    except OSError as exc:
                        _log.exception(
                            "Legacy .bak write failed for %s — aborting "
                            "save to preserve existing data", path,
                        )
                        raise OSError(
                            f"legacy .bak write failed for {label} "
                            f"({path.name}): {exc}. Save aborted to "
                            f"keep the previous-good file intact."
                        ) from exc
                # Count + extract entries for the shrink guard. Accept
                # both envelope and legacy bare-list formats so the
                # first save after an upgrade doesn't false-positive.
                #
                # Sweep #35 (2026-05-26): if the prior file is valid
                # JSON but `_extract_entries` returns `None` (unknown
                # schema shape, mangled envelope, hand-edited file),
                # `existing_count` stays 0 → the shrink guard at the
                # next step won't fire even when the new save drops
                # data on the floor. The timestamped `.bak.<ts>` we
                # already wrote above is the user's primary recovery
                # path; we additionally drop the raw bytes into
                # `lost_entries/` so a grep-by-content recovery has
                # something to chew on without needing to find the
                # rotation by date. Best-effort only — failure here
                # never blocks the save.
                try:
                    prev = json.loads(existing)
                    extracted, _ = _extract_entries(prev, label)
                    if extracted is not None:
                        prev_entries = extracted
                        existing_count = len(extracted)
                    else:
                        _log.warning(
                            "%s: prior file parsed as JSON but yielded "
                            "no extractable entries (schema mismatch). "
                            "Shrink guard cannot compare counts; "
                            "spilling raw bytes for recovery.",
                            label,
                        )
                        _spill_raw_bytes(path, existing, label,
                                          reason="extract-returned-none")
                except (json.JSONDecodeError, ValueError):
                    _log.warning(
                        "%s: prior file is non-empty but not valid "
                        "JSON. Shrink guard cannot compare counts; "
                        "spilling raw bytes for recovery.",
                        label,
                    )
                    _spill_raw_bytes(path, existing, label,
                                      reason="invalid-json")
        except OSError as exc:
            # 2026-05-28 (sweep #30): the file EXISTS (we're inside
            # `if path.exists()`) but could not be read to back it up.
            # Falling through to Step 3 would overwrite real,
            # un-backed-up data with the in-memory state — and if that
            # state is empty/short because the SAME I/O fault made the
            # cache load `[]`, the good file is destroyed with no
            # recovery copy AND the shrink guard can't fire
            # (existing_count stayed 0). A save we cannot back up is not
            # safe: refuse so the caller surfaces it. (A legitimately
            # ABSENT file never reaches here — only the present-but-
            # unreadable case does.)
            msg = (
                f"Refusing to save {label}: prior file at {path} exists "
                f"but could not be read to back it up ({exc}). "
                f"Overwriting would destroy un-backed-up data. Resolve "
                f"the read error (permissions / disk) and retry."
            )
            _log.error(msg)
            raise OSError(msg) from exc

    # Step 2: shrink guard with spillover + L3 catastrophic-shrink
    # refusal (2026-05-22).
    #
    # **Three tiers** as the loss ratio escalates:
    #
    #   * ANY shrink → log a warning.
    #   * SUSPICIOUS shrink (>50% loss, base population >=5) → spill
    #     the discarded entries to `lost_entries/` so data is never
    #     silently destroyed. Save still proceeds (user may have
    #     legitimately bulk-deleted).
    #   * CATASTROPHIC shrink (>90% loss, base population >=10) →
    #     **REFUSE the save** with `RuntimeError` unless wrapped in
    #     `_allow_catastrophic_shrink()`. The 00:37 incident — running
    #     app at startup writing 43 bytes over a 156 MB library — is
    #     exactly this signature. Legitimate catastrophic shrinks
    #     (`_restore_from_backup` restoring a tiny backup over a large
    #     live file, hypothetical Master-Delete-via-save flows) opt in
    #     via the context manager.
    if existing_count > 0 and len(entries) < existing_count:
        mirror_swap = _mirror_swap_depth() > 0
        _log.warning(
            "SHRINK GUARD: %s is being overwritten with %d entries "
            "(was %d)%s. If this is unexpected, restore from %s.bak "
            "or %s.bak.<timestamp>.",
            label, len(entries), existing_count,
            " [expected mirror swap — the dropped entries remain in the "
            "sibling collections file, not lost]" if mirror_swap else "",
            path, path,
        )
        if mirror_swap:
            # Expected mirror swap (active-collection / parts-bin / primer
            # / experiment switch). The "dropped" entries are NOT lost —
            # they remain in the sibling collections file that is the real
            # source of truth — so we neither spill a redundant
            # `lost_entries/` copy (pure waste: a ~150 MB write per
            # collection switch, the cause of the 1.5 GB residue) nor
            # refuse the write on a >90% shrink (the latent bug where a
            # big→tiny collection switch raised RuntimeError). See
            # `_expected_mirror_swap` / `_safe_save_json_mirror` /
            # `_switch_active_collection_library`.
            pass
        else:
            suspicious = (existing_count >= 5
                          and len(entries) < existing_count // 2)
            catastrophic = (existing_count >= 10
                            and len(entries) * 10 < existing_count)
            if suspicious and prev_entries is not None:
                lost = _diff_lost_entries(prev_entries, entries)
                spilled = _spill_lost_entries(path, lost, label)
                if spilled is not None:
                    _log.warning(
                        "SHRINK GUARD: %d %s entries dumped to %s before "
                        "overwrite — recoverable on user request.",
                        len(lost), label, spilled,
                    )
            if catastrophic and _catastrophic_shrink_depth() <= 0:
                raise RuntimeError(
                    f"refusing to write {label!r}: catastrophic shrink "
                    f"({existing_count} → {len(entries)} entries, "
                    f"{100 * (existing_count - len(entries)) / existing_count:.1f}% loss). "
                    f"This signature matches the 2026-05-22 incident "
                    f"(running app at startup writing 43 bytes over a "
                    f"156 MB library). The discarded entries have been "
                    f"spilled to `lost_entries/` so nothing is lost. "
                    f"If this save is genuinely intentional (Restore-from-"
                    f"backup, programmatic data wipe), wrap the call in "
                    f"`with splicecraft._allow_catastrophic_shrink():`. "
                    f"For an active-collection/bin/primer switch use "
                    f"`_switch_active_collection_library` / "
                    f"`_safe_save_json_mirror` instead."
                )

    # Step 3: atomic write — tempfile in same dir → os.replace.
    payload = {"_schema_version": schema_version, "entries": entries}
    try:
        fd, tmp_name = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2)
                fh.flush()
                # 2026-05-28 (sweep #30): re-raise fsync failures here
                # too — mirrors the H2 hardening in `_atomic_write_text`
                # / `_atomic_write_bytes`. This inline writer was the
                # last save path that still swallowed `os.fsync` errors:
                # on EIO / ENOSPC at fsync the data never reached stable
                # storage, yet `os.replace` proceeded and the function
                # returned success — the UI showed "saved" over a write
                # a power-loss would lose. The surrounding handler unlinks
                # the temp file and the error propagates (sacred #7).
                os.fsync(fh.fileno())
            # Read-back validation BEFORE the atomic swap (2026-06-03):
            # confirm the temp file we just wrote + fsynced is actually
            # parseable before it replaces the live file. `json.dump`
            # raising is already handled (temp unlinked, main untouched),
            # but a silent truncation / FS fault between write and rename
            # would otherwise be promoted to the live file. Size-gated:
            # files <= `_state._SAVE_READBACK_FULL_PARSE_MAX_BYTES` get a full
            # `json.loads` + entry-count check; larger files get a cheap
            # tail check (non-empty + closes with `}`) that catches
            # truncation without a multi-second re-parse. On failure raise
            # so the surrounding handler unlinks the temp and leaves the
            # live file intact (sacred #7).
            try:
                tmp_size = os.path.getsize(tmp_name)
            except OSError:
                tmp_size = -1
            if tmp_size == 0:
                raise OSError(
                    f"refusing to commit {label}: temp file is empty "
                    f"after write ({len(entries)} entries expected)"
                )
            if 0 < tmp_size <= _state._SAVE_READBACK_FULL_PARSE_MAX_BYTES:
                with open(tmp_name, "r", encoding="utf-8") as _vfh:
                    _reparsed = json.load(_vfh)
                if not (isinstance(_reparsed, dict)
                        and isinstance(_reparsed.get("entries"), list)
                        and len(_reparsed["entries"]) == len(entries)):
                    raise OSError(
                        f"refusing to commit {label}: temp-file read-back "
                        f"did not match the payload "
                        f"({len(entries)} entries written)"
                    )
            elif tmp_size > _state._SAVE_READBACK_FULL_PARSE_MAX_BYTES:
                with open(tmp_name, "rb") as _vfh:
                    _vfh.seek(max(0, tmp_size - 64))
                    _tail = _vfh.read()
                if not _tail.rstrip().endswith(b"}"):
                    raise OSError(
                        f"refusing to commit {label}: temp file appears "
                        f"truncated (tail does not close the JSON object)"
                    )
            os.replace(tmp_name, str(path))
            # Fsync the parent directory so the rename's directory entry
            # update is journalled — see `_fsync_parent_dir`. Pre-fix
            # `_safe_save_json` skipped this, leaving the
            # rename-then-power-loss window where the inode held the new
            # data but the directory entry still pointed at the OLD inode
            # after fsck.
            _fsync_parent_dir(path)
            _log.info("Saved %s: %d entries to %s (schema v%d)",
                      label, len(entries), path, schema_version)
        except Exception:
            try:
                os.unlink(tmp_name)
            except OSError:
                pass
            raise
    except Exception:
        _log.exception("Failed to save %s to %s", label, path)
        raise

    # Step 4: prune timestamped backups outside the retention window.
    # Done AFTER a successful save so a write failure doesn't prune
    # the very backup we'd want to recover from.
    _prune_backups(path)

def _backup_info(path: Path) -> "dict | None":
    """Parse `path` as a SpliceCraft persistence file (envelope or
    legacy bare-list) and return ``{n_entries, mtime_str, error}``.
    A damaged backup returns ``n_entries=None`` and a non-empty
    ``error`` string so the Restore UI can surface it tagged
    ``[damaged]`` instead of silently dropping it — users trying to
    recover need to see what was there, even if parsing fails.

    Returns ``None`` only when the file is structurally absent or
    so wrong (e.g. cannot even stat) that no row should be shown.

    Size-capped at `_state._SAFE_LOAD_JSON_MAX_BYTES` (mirrors `_safe_load_json`'s
    1 GB cap) — a corrupted/oversized legacy `.bak` could otherwise OOM
    the Restore modal on open. Symlink-rejected via the same lstat path.
    """
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return None
    ts = _datetime.fromtimestamp(mtime).strftime("%Y-%m-%d %H:%M:%S")
    ok, reason = _safe_file_size_check(
        path, _state._SAFE_LOAD_JSON_MAX_BYTES, "backup",
    )
    if not ok:
        return {"n_entries": None, "mtime_str": ts,
                "error": reason or "size/symlink check failed"}
    try:
        raw = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"n_entries": None, "mtime_str": ts,
                "error": f"parse failed: {exc}"}
    if isinstance(raw, list):
        n = len(raw)
    elif isinstance(raw, dict):
        entries = raw.get("entries", [])
        n = len(entries) if isinstance(entries, list) else 0
    else:
        return {"n_entries": None, "mtime_str": ts,
                "error": "unexpected JSON shape (not list or envelope)"}
    return {"n_entries": n, "mtime_str": ts, "error": ""}

def _safe_load_json(path: Path, label: str) -> "tuple[list, str | None]":
    """Load a JSON payload from `path`. Returns (entries, warning_or_None).

    Accepts both the current envelope format (`{"_schema_version": N,
    "entries": [...]}`) and the legacy flat-list format written by
    SpliceCraft < 0.3.1. The legacy file gets silently rewritten as an
    envelope on the next save.

    Capped at ``_state._SAFE_LOAD_JSON_MAX_BYTES`` (1 GB) — same cap as bulk
    import. A corrupted / mis-restored / hostile shared library would
    otherwise OOM on read.

    - Missing file → ([], None) — normal first run, no warning.
    - Valid file   → (entries, None).
    - Oversized file → ([], warning). The matching guard in
      `_safe_save_json` REFUSES to overwrite an oversized file so
      the in-memory empty state can't silently nuke the user's data.
    - Corrupt file → attempt .bak restore; if .bak is valid →
      (bak_entries, warning). If .bak also corrupt → ([], warning).
    """
    if not path.exists():
        # Normal first run (no file AND no backups) → empty, no warning.
        # But a surviving `.bak` / rotated backup while the main file is GONE is
        # the signature of a failed recovery-rewrite (a prior load renamed the
        # corrupt main aside, then the bak→main copy hit ENOSPC / RO mount, so
        # the main is absent but the good data is intact in the backup). The old
        # early-return then reported an EMPTY library though the data survived.
        # Recover from the backup chain (legacy `.bak` first, then rotations
        # newest→oldest) and re-materialise the main file for the next save.
        for _cand in ([path.with_suffix(path.suffix + ".bak")]
                      + list(reversed(_iter_backups(path)))):
            if not _cand.exists():
                continue
            ok_c, _r = _safe_file_size_check(
                _cand, _state._SAFE_LOAD_JSON_MAX_BYTES, label)
            if not ok_c:
                continue
            try:
                _raw = json.loads(_read_backup_bytes(_cand))
            except Exception:
                continue
            _entries, _ = _extract_entries(_raw, label)
            if _entries is None:
                continue
            _log.warning(
                "%s: main file %s is missing but recovered %d entries from "
                "backup %s.", label, path, len(_entries), _cand.name)
            try:
                if _recovery_write_allowed(path, label):
                    _atomic_write_bytes(path, _read_backup_bytes(_cand))
            except OSError:
                _log.warning("Could not rewrite %s from backup %s",
                             path, _cand.name)
            return _entries, (
                f"{label} main file was missing — restored "
                f"{len(_entries)} entries from backup {_cand.name}.")
        return [], None

    # Size cap + symlink rejection — refuse to read multi-GB files
    # and refuse to follow symlinks (a symlink → /dev/zero reports
    # 0 bytes through `path.stat()` and a subsequent read would
    # consume RAM until the kernel OOM-killed us).
    ok, reason = _safe_file_size_check(
        path, _state._SAFE_LOAD_JSON_MAX_BYTES, label,
    )
    if not ok:
        _log.warning("%s: %s", path, reason)
        return [], reason

    # Try the main file
    main_warning: "str | None" = None
    try:
        raw = json.loads(path.read_text(encoding="utf-8-sig"))
        entries, shape_warn = _extract_entries(raw, label)
        if entries is not None:
            # If this file was written by a newer SpliceCraft, stash
            # the observed schema version keyed by absolute path so
            # the next save preserves it (no schema-stamp demotion on
            # downgrade-then-edit round trips). The version is on the
            # envelope; pull it out here so we don't need to re-parse.
            if isinstance(raw, dict):
                v = raw.get("_schema_version")
                if isinstance(v, int) and v > _state._CURRENT_SCHEMA_VERSION:
                    _state._OBSERVED_SCHEMA_VERSIONS[str(path)] = v
            return entries, shape_warn
        _log.warning("%s: %s", path, shape_warn)
        main_warning = shape_warn
    except Exception:
        _log.exception("Corrupt %s file: %s", label, path)

    # Main file is corrupt — try the .bak. Apply the same size cap
    # + symlink rejection here: if the main file was small/corrupt
    # but the `.bak` happens to be oversized, the recovery path
    # would otherwise OOM while the main load was safely refused.
    bak = path.with_suffix(path.suffix + ".bak")
    if bak.exists():
        ok_bak, _reason = _safe_file_size_check(
            bak, _state._SAFE_LOAD_JSON_MAX_BYTES, label,
        )
        if not ok_bak:
            _log.warning("Backup %s rejected: oversized/symlink", bak)
            return [], (
                main_warning
                or f"{label} is corrupt and the backup was rejected. "
                "Starting empty."
            )
        try:
            raw = json.loads(bak.read_text(encoding="utf-8"))
            entries, _ = _extract_entries(raw, label)
            if entries is not None:
                _log.info("Restored %s from backup %s (%d entries)",
                          label, bak, len(entries))
                # 2026-05-27 (audit-2 H3): also preserve the schema-
                # version stamp from the backup. Pre-fix
                # `_state._OBSERVED_SCHEMA_VERSIONS` was only populated when
                # the main file parsed cleanly, so a `.bak` from a
                # newer SpliceCraft that we recovered through here
                # would be re-saved at `_state._CURRENT_SCHEMA_VERSION`,
                # silently demoting the schema. Now we capture the
                # bak's envelope version too.
                if isinstance(raw, dict):
                    bak_v = raw.get("_schema_version")
                    if isinstance(bak_v, int) and bak_v > _state._CURRENT_SCHEMA_VERSION:
                        _state._OBSERVED_SCHEMA_VERSIONS[str(path)] = bak_v
                # 2026-05-27 (audit-2 M5): rename the corrupt main
                # aside to ``<file>.corrupt-<ts>`` BEFORE the bak-
                # to-main copy so a forensic look at what corrupted
                # the main file is still possible. Best-effort; on
                # failure proceed (the corrupt file gets overwritten
                # as it always did).
                try:
                    corrupt_ts = _datetime.now().strftime("%Y%m%d-%H%M%S")
                    corrupt_aside = path.with_name(
                        f"{path.name}.corrupt-{corrupt_ts}",
                    )
                    bump = 0
                    while corrupt_aside.exists():
                        bump += 1
                        corrupt_aside = path.with_name(
                            f"{path.name}.corrupt-{corrupt_ts}.{bump}",
                        )
                    if _recovery_write_allowed(path, label):
                        path.rename(corrupt_aside)
                except OSError:
                    _log.warning(
                        "Could not preserve corrupt main file aside "
                        "for forensic inspection — overwriting in place",
                    )
                # Atomically overwrite the corrupt main file with the
                # good backup. Pre-0.8.9 used `shutil.copy2(bak, path)`
                # which is a non-atomic stream copy: a power loss or
                # process kill mid-copy would leave the main file
                # truncated — paradoxically *less* recoverable than
                # the corrupt state we were rescuing from. Routing
                # through `_atomic_write_bytes` (tempfile + replace +
                # parent-dir fsync) preserves the recovery's guarantee.
                try:
                    if _recovery_write_allowed(path, label):
                        _atomic_write_bytes(path, bak.read_bytes())
                except OSError:
                    _log.warning(
                        "Could not rewrite main file %s from backup",
                        path,
                    )
                return entries, (
                    f"{label} was corrupt — restored {len(entries)} entries "
                    f"from backup."
                )
        except Exception:
            _log.exception("Backup %s also corrupt: %s", label, bak)

    # Recovery chain (2026-06-03): the legacy `.bak` is the fastest
    # recovery, but if BOTH the main file and `.bak` are corrupt (e.g.
    # two bad saves in a row — `.bak` is overwritten every save), the
    # timestamped rotation still holds older good copies. Walk them
    # newest → oldest, transparently decompressing gzipped backups
    # (`_compress_old_backups`), and restore from the first that parses.
    # The main file is still present here (the legacy block only renames
    # it aside on SUCCESS), so we preserve it for forensics then rewrite.
    for chain_bak in reversed(_iter_backups(path)):
        ok_chain, _r = _safe_file_size_check(
            chain_bak, _state._SAFE_LOAD_JSON_MAX_BYTES, label,
        )
        if not ok_chain:
            continue
        try:
            chain_raw = json.loads(_read_backup_bytes(chain_bak))
        except Exception:
            continue
        chain_entries, _shape = _extract_entries(chain_raw, label)
        if chain_entries is None:
            continue
        _log.info("Restored %s from rotated backup %s (%d entries)",
                  label, chain_bak.name, len(chain_entries))
        if isinstance(chain_raw, dict):
            cv = chain_raw.get("_schema_version")
            if isinstance(cv, int) and cv > _state._CURRENT_SCHEMA_VERSION:
                _state._OBSERVED_SCHEMA_VERSIONS[str(path)] = cv
        try:
            corrupt_ts = _datetime.now().strftime("%Y%m%d-%H%M%S")
            corrupt_aside = path.with_name(f"{path.name}.corrupt-{corrupt_ts}")
            bump = 0
            while corrupt_aside.exists():
                bump += 1
                corrupt_aside = path.with_name(
                    f"{path.name}.corrupt-{corrupt_ts}.{bump}")
            if path.exists() and _recovery_write_allowed(path, label):
                path.rename(corrupt_aside)
        except OSError:
            _log.warning(
                "Could not preserve corrupt main file aside for %s — "
                "overwriting in place", path,
            )
        try:
            if _recovery_write_allowed(path, label):
                _atomic_write_bytes(path, _read_backup_bytes(chain_bak))
        except OSError:
            _log.warning(
                "Could not rewrite main file %s from rotated backup", path,
            )
        return chain_entries, (
            f"{label} and its .bak were corrupt — restored "
            f"{len(chain_entries)} entries from rotated backup "
            f"{chain_bak.name}."
        )

    return [], (main_warning
                or f"{label} is corrupt and no valid backup was found. "
                   "Starting empty.")
