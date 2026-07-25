"""splicecraft_backup — user-data backup / restore / migrate engine (Phase D, layer L1).

The whole-data-directory lifecycle operations, extracted from the hub so the
data-safety machinery is one bounded, testable unit:

  * Daily launch-time snapshots (`_snapshot_data_files` / `_prune_old_snapshots`)
    + data-dir housekeeping (`_run_data_dir_housekeeping`).
  * The user-data file/dir REGISTRIES (`_USER_DATA_FILE_ATTRS` / `_USER_DATA_DIR_ATTRS`
    / `_OPERATIONAL_FILE_ATTRS`) — the canonical "what is user data" definition that
    backup, Master Delete, and migrate-archive all drive from.
  * Master Delete TARGET ENUMERATION (`_master_delete_*_targets` / `_residual_paths`
    / `_log_files`) — the pure "which paths are user data" lists; the actual
    app-coupled wipe (`_perform_master_delete`) stays hub-side.
  * Pre-update snapshots for `splicecraft update` rollback (`_create_pre_update_snapshot`
    / `_list_pre_update_snapshots` / `_validate_snapshot_member`
    / `_restore_pre_update_snapshot` / `_enforce_pre_update_retention`
    + `_resolve_pre_update_backup_dir`).
  * Migrate-archive export/import for moving data between machines
    (`_export_migrate_archive` / `_import_migrate_archive`).

SACRED: every write goes through persistence's atomic-write + chokepoint
(`_atomic_write_*` re-exported from L0). Two hub-pinned bits are reached through
`_state` hooks (no upward import, behaviour-identical): `_resolve_data_attr_hook`
(the hub's `_resolve_state_or_hub`, whose `globals()` fallback must run hub-side)
and `_gc_orphan_blobs_hook` (the blob store stays hub-side). Re-exported by the
hub so `sc.<name>` + every call site resolves unchanged.

Imports only L0 siblings (persistence / util / logging / state) → layer L1.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import sys
from collections import namedtuple as _namedtuple
from datetime import date as _date
from pathlib import Path

import splicecraft_state as _state
from splicecraft_logging import _log, _log_event, _timed
from splicecraft_persistence import (
    _allow_catastrophic_shrink,
    _atomic_write_bytes,
    _atomic_write_text,
    _backup_info,
    _refuse_unauthorized_delete,
    _refuse_unauthorized_write,
    _compress_old_backups,
    _extract_entries,
    _iter_backups,
    _prune_backups,
    _prune_lost_entries,
    _safe_file_size_check,
    _safe_save_json,
)
from splicecraft_util import _RUNTIME_PLATFORM, _is_windows_reserved_stem


# ── Daily launch-time snapshot ────────────────────────────────────────────────
#
# Belt-and-braces above `_safe_save_json`'s rotating backups: on every
# new calendar day the user starts the app, we copy each persistent
# JSON file to `<DATA_DIR>/snapshots/<file>-YYYY-MM-DD.json`. That way a
# bad upgrade or worker race that writes empty data over the weekend can
# always be reverted from a known-good snapshot up to a month old.
#
# Snapshots are pruned once they're older than `_SNAPSHOT_RETENTION_DAYS`
# (30) so the directory doesn't grow without bound. Per-file size is
# tiny (~100 KB typical library), so 30 days × 4 files ≈ 12 MB at most.
#
# Files snapshotted: plasmid library, collections, parts bin, primers.
# (.dna sidecars and feature library are deliberately excluded — the
# former is byte-identical to the user's source files, the latter is
# rebuildable from the active library.)
# (snapshot subdir name moved to splicecraft_state, Phase B-main.)
_SNAPSHOT_RETENTION_DAYS = 30


# Per-file cap for daily snapshots, expressed as a fraction of the
# aggregate ceiling below (kept in sync by construction, see [INV-167]).
#
# This used to be a flat 50 MB, which meant the BIGGEST and most
# expensive-to-lose library was the one file with NO daily snapshot —
# the worst possible place for the safety net to have a hole, and a real
# 2026-07-25 data-loss recovery turned on precisely that snapshot
# existing. Disk growth is already bounded properly by
# `_SNAPSHOT_TOTAL_SIZE_CAP` in `_prune_old_snapshots` (oldest-first),
# so the only job left for a per-file cap is to refuse a file so large
# that snapshotting it would evict every other snapshot for no net gain.
# Half the aggregate budget guarantees at least two generations can
# coexist. Anything skipped here still has `_safe_save_json`'s .bak
# rotation (10 generations) and the pre-update snapshots.
#
# NB: compression was evaluated as an alternative to skipping and
# REJECTED — `_snapshot_data_files` runs SYNCHRONOUSLY on the launch
# path, and gzip level 6 measures ~1.2 s per 43 MB (≈4.5 s for a 160 MB
# library) versus ~0.3 s for a plain copy. That is exactly the hang the
# sibling housekeeping pass spawns a daemon thread to avoid. Do not
# re-add compression here without moving the work off the launch path.
_SNAPSHOT_FILE_SIZE_CAP  = 250 * 1024 * 1024  # = _SNAPSHOT_TOTAL_SIZE_CAP // 2


# Aggregate cap across the entire snapshots directory. Per-file cap
# alone allows 4 files × 49 MB × 30 days = ~5.9 GB of accumulated
# snapshots. The aggregate cap drops the oldest snapshots until total
# usage falls below this ceiling; nothing structural is lost because
# the `.bak` rotation (per-file, 10 generations) covers recent rollback
# and the pre-update snapshots cover upgrade rollback. Date-based
# retention runs first; the aggregate cap is the last-line defense
# for genuinely-heavy users.
_SNAPSHOT_TOTAL_SIZE_CAP = 500 * 1024 * 1024  # 500 MB


# ── The canonical persistent-cache registry (single source of truth) ──────────
# One entry per persistent cache / user-data file. The parallel enumerations
# across the data-safety machinery are all VIEWS of this ONE table:
#   * `_USER_DATA_FILE_ATTRS` (backup / wipe / migrate)  = the `file_attr`s where `user_data`
#   * `_MASTER_DELETE_CACHE_ATTRS` (agent; cache-bust)   = the `cache_attr`s where `master_delete`
#   * `_AGENT_BACKUP_LABELS` (backup list/restore)        = `{label: file_attr}` where `label`
#   * the pytest sandbox pairs (conftest `_protect_user_data`) = `(file_attr, cache_attr)` where `user_data`
# The runtime lists below are kept EXPLICIT (not derived) so the data-safety
# code stays auditable, but the `test_registry_*` cases in
# `tests/test_data_safety.py` (class TestRealFilesNeverTouched) assert each one
# matches this registry EXACTLY — so the recurring cache-bust / backup-label
# enumeration drift (INV-43 / 50 / 64 / 66 / 142) can't return: add a cache HERE
# and the failing test names precisely which views still need it. The pytest
# sandbox (conftest `_SANDBOXED_DATA_FILE_ATTRS`) is its OWN hand-maintained
# list — NOT derived from this registry — cross-checked against it by
# `test_registry_matches_sandbox_pairs` (so a new cache still needs a conftest
# entry, and the guard fails loudly until it has one; no silent data-loss).
#
# Order follows `_USER_DATA_FILE_ATTRS`; `_MASTER_DELETE_CACHE_ATTRS` preserves it
# and interleaves the two index-only caches (`file_attr`/`label` = None — a
# derived lookup index with no on-disk file, busted on restore but never backed up).
_CacheSpec = _namedtuple(
    "_CacheSpec", ("label", "file_attr", "cache_attr", "user_data", "master_delete"))
_PERSISTENT_CACHE_REGISTRY: "tuple[_CacheSpec, ...]" = (
    _CacheSpec("plasmid_library",       "_LIBRARY_FILE",               "_library_cache",               True,  True),
    _CacheSpec("collections",           "_COLLECTIONS_FILE",           "_collections_cache",           True,  True),
    _CacheSpec("parts_bin",             "_PARTS_BIN_FILE",             "_parts_bin_cache",             True,  True),
    _CacheSpec("parts_bin_collections", "_PARTS_BIN_COLLECTIONS_FILE", "_parts_bin_collections_cache", True,  True),
    _CacheSpec("primers",               "_PRIMERS_FILE",               "_primers_cache",               True,  True),
    _CacheSpec("primer_collections",    "_PRIMER_COLLECTIONS_FILE",    "_primer_collections_cache",    True,  True),
    _CacheSpec(None,                    None,                          "_primer_usage_cache",          False, True),   # index (sweep #20)
    _CacheSpec("features",              "_FEATURES_FILE",              "_features_cache",              True,  True),
    _CacheSpec("feature_colors",        "_FEATURE_COLORS_FILE",        "_feature_colors_cache",        True,  True),
    _CacheSpec(None,                    None,                          "_feature_library_index_cache", False, True),   # index (sweep #20)
    _CacheSpec("grammars",              "_GRAMMARS_FILE",              "_grammars_cache",              True,  True),
    _CacheSpec("entry_vectors",         "_ENTRY_VECTORS_FILE",         "_entry_vectors_cache",         True,  True),
    _CacheSpec("codon_tables",          "_CODON_TABLES_FILE",          "_codon_tables_cache",          True,  True),
    _CacheSpec("settings",              "_SETTINGS_FILE",              "_settings_cache",              True,  True),
    _CacheSpec("experiments",           "_EXPERIMENTS_FILE",           "_experiments_cache",           True,  True),
    _CacheSpec("experiment_projects",   "_EXPERIMENT_PROJECTS_FILE",   "_experiment_projects_cache",   True,  True),
    _CacheSpec("gels",                  "_GELS_FILE",                  "_gels_cache",                  True,  True),
    _CacheSpec("protein_motifs",        "_PROTEIN_MOTIFS_FILE",        "_protein_motifs_cache",        True,  True),
    _CacheSpec("custom_enzymes",        "_CUSTOM_ENZYMES_FILE",        "_custom_enzymes_cache",        True,  True),
    _CacheSpec("enzyme_collections",    "_ENZYME_COLLECTIONS_FILE",    "_enzyme_collections_cache",    True,  True),
    _CacheSpec("hmm_db_catalog",        "_HMM_DB_CATALOG_FILE",        "_hmm_db_catalog_cache",        True,  True),
    _CacheSpec("protein_collections",   "_PROTEIN_COLLECTIONS_FILE",   "_protein_collections_cache",   True,  True),
    _CacheSpec("model_collections",     "_MODEL_COLLECTIONS_FILE",     "_model_collections_cache",     True,  True),
    _CacheSpec("protocol_collections",  "_PROTOCOL_COLLECTIONS_FILE",  "_protocol_collections_cache",  True,  True),
    _CacheSpec("custom_labware",        "_CUSTOM_LABWARE_FILE",        "_custom_labware_cache",        True,  True),
)


# Canonical user-data files. Order matters for the manifest so a future
# audit can replay the snapshot in deterministic sequence; do not reorder.
# (Must equal the `user_data` `file_attr`s of `_PERSISTENT_CACHE_REGISTRY` above —
# enforced by test_registry_matches_user_data_file_attrs.)
_USER_DATA_FILE_ATTRS: tuple = (
    "_LIBRARY_FILE",          # plasmid_library.json — the user's plasmids
    "_COLLECTIONS_FILE",      # collections.json — collection definitions
    "_PARTS_BIN_FILE",        # parts_bin.json — active bin's parts (mirror)
    "_PARTS_BIN_COLLECTIONS_FILE",   # parts_bin_collections.json — all bins
    "_PRIMERS_FILE",          # primers.json — primer designs
    "_PRIMER_COLLECTIONS_FILE",  # primer_collections.json — named primer collections
    "_FEATURES_FILE",         # features.json — feature library
    "_FEATURE_COLORS_FILE",   # feature_colors.json — feature colours
    "_GRAMMARS_FILE",         # cloning_grammars.json — custom grammars
    "_ENTRY_VECTORS_FILE",    # entry_vectors.json — entry vectors
    "_CODON_TABLES_FILE",     # codon_tables.json — codon usage tables
    "_SETTINGS_FILE",         # settings.json — persisted user toggles
    "_EXPERIMENTS_FILE",      # experiments.json — lab-notebook entries
    "_EXPERIMENT_PROJECTS_FILE",  # experiment_projects.json — all projects
    "_GELS_FILE",             # gels.json — saved agarose-gel snapshots
    "_PROTEIN_MOTIFS_FILE",   # protein_motifs.json — user motif overrides (sweep #15)
    "_CUSTOM_ENZYMES_FILE",   # custom_enzymes.json — user-added restriction enzymes
    "_ENZYME_COLLECTIONS_FILE",  # enzyme_collections.json — named subsets of master catalog
    "_HMM_DB_CATALOG_FILE",   # hmm_db_catalog.json — registry of HMM databases (sweep #28)
    "_PROTEIN_COLLECTIONS_FILE",  # protein_collections.json — named protein-sequence collections
    "_MODEL_COLLECTIONS_FILE",  # model_collections.json — BABS model-picker collections (INV-139)
    "_PROTOCOL_COLLECTIONS_FILE",  # protocol_collections.json — saved OT-2 protocol designs
    "_CUSTOM_LABWARE_FILE",  # custom_labware.json — custom Opentrons labware definitions
)


# User-data sub-directories — autosaved unsaved-edits + .dna sidecars.
# Both are user-generated content the user might be relying on for recovery.
_USER_DATA_DIR_ATTRS: tuple = (
    "_CRASH_RECOVERY_DIR",   # autosaved .gb files for unsaved records
    "_DNA_ORIGINALS_DIR",    # .dna sidecars (CommercialSaaS round-trip)
    "_EXPERIMENTS_DIR",      # lab-notebook image attachments (per-entry)
    "_PLUGINS_DIR",          # reserved for future plugin storage (empty for now)
    "_HMM_DATABASES_DIR",    # downloaded HMM databases (sweep #28) — Pfam-A etc.
    "_PLASMID_BLOB_DIR",     # content-addressed gb_text blobs (v1.0.23) — sequences
)


# Operational files: NOT user data. These are listed explicitly so that
# the `test_every_data_file_constant_is_classified` audit catches any
# new `_*_FILE` constant that a future contributor adds without
# deciding which list it belongs in. If you add a new `_*_FILE` and
# the audit fails, either:
#   * append it to `_USER_DATA_FILE_ATTRS` (it carries data the user
#     would lose if it were destroyed), OR
#   * append it here (it's regenerated/transient — agent tokens,
#     migration markers, ephemeral state).
_OPERATIONAL_FILE_ATTRS: tuple = (
    "_AGENT_TOKEN_FILE",      # ephemeral; regenerated each --agent-api launch
    "_DATA_VERSION_FILE",     # last-touched-version stamp; tiny, regenerated
)


# How many pre-update snapshots to keep before pruning oldest. 5 covers
# the realistic "I updated last week and now things look weird" window
# while bounding worst-case disk usage to ~5× the user's library size.
_PRE_UPDATE_SNAPSHOT_RETENTION = 5


# Filenames inside each snapshot directory.
_PRE_UPDATE_MANIFEST_NAME = "manifest.json"


_PRE_UPDATE_STAGING_PREFIX = ".tmp-"


# Manifests are pure JSON metadata (file list + checksums + small
# headers). Real manifests run ~2-16 KB; a few-MB cap leaves slack for
# unforeseen growth while refusing a planted multi-GB manifest that
# would OOM the launch-time listing or the restore handler. The
# directory is user-writable (sibling of _state._DATA_DIR by default), so a
# defensive cap matters even though normal flows never approach it.
_PRE_UPDATE_MANIFEST_MAX_BYTES = 4 * 1024 * 1024


# Snapshot manifest schema version. Bump this when the on-disk shape
# changes incompatibly (e.g., a new required field). Restore code
# refuses any manifest whose `schema_version` exceeds this — better
# to surface a clear error than to risk a silent partial restore. A
# manifest with a *lower* version still loads (additive fields are
# read with `.get(...)` defaults).
_PRE_UPDATE_SCHEMA_VERSION = 1


# ── Migrate Data: one portable file for ALL user data ─────────────────────────
#
# A migrate archive is a single ZIP wrapping a pre-update-style snapshot of the
# canonical user-data registry (`_USER_DATA_FILE_ATTRS` + `_USER_DATA_DIR_ATTRS`)
# — the SAME source of truth the pre-update snapshot uses, so any data file is
# migrated the moment it joins the registry (the classification test in
# tests/test_smoke.py FAILS the build if a new persisted file is left out). The
# export is atomic (temp file + os.replace); the import reuses
# `_restore_pre_update_snapshot` (automatic pre-import backup + sha256-verified
# atomic per-file replace). See `[INV-116]` in docs/invariants.md — keep this in
# lockstep with the registry whenever a new data structure is added (INV-116).
_MIGRATE_MARKER_NAME    = "splicecraft-migrate.json"  # top-level archive id file


_MIGRATE_FORMAT_VERSION = 1                            # bump on layout change


_MIGRATE_ARCHIVE_SUBDIR = "data"                       # snapshot dir inside zip


# The downloaded HMM databases are large, re-downloadable REFERENCE data (Pfam
# et al.), not the user's custom work — their catalog json IS migrated, so the
# new machine re-downloads them on demand. Excluded by default; the export modal
# can opt in. (Any other re-derivable heavy tree added later goes here.)
_MIGRATE_DEFAULT_EXCLUDE_ATTRS = frozenset({"_HMM_DATABASES_DIR"})


# Runaway-extract / zip-bomb guards, checked BEFORE extraction. A real migration
# of a heavy library is far under these; a small archive claiming to expand past
# them is refused before it can fill the disk. (sha256 verification in the
# restore is the second line of defence against tampered contents.)
_MIGRATE_MAX_TOTAL_UNCOMPRESSED = 64 * 1024 ** 3       # 64 GiB ceiling


_MIGRATE_MAX_MEMBERS            = 2_000_000            # absurd file-count guard


# Anchored regex for snapshot directory names. Retention pruning uses
# this so a maliciously symlinked or misconfigured `backup_dir` (e.g.
# pointing at `/`, `~`, or a system path with bin/etc/home subdirs)
# still cannot rmtree directories we did not create. The pattern
# matches every name `_create_pre_update_snapshot` ever produces —
# `<8-digit date>-<6-digit time>-<8-hex random>__from-<token>` plus
# an optional `.<int>` collision-bump suffix.
_PRE_UPDATE_NAME_RE = re.compile(
    r"^\d{8}-\d{6}-[0-9a-f]{8}__from-[A-Za-z0-9._-]+(?:\.\d+)?$"
)


def _snapshot_data_files(data_dir: Path,
                          paths: "list[Path] | None" = None) -> "list[Path]":
    """Copy each path in `paths` to `data_dir/snapshots/<stem>-YYYY-MM-DD.<ext>`
    if today's snapshot doesn't already exist. Skips missing / empty
    files (no point preserving an empty library). Returns the list of
    snapshot paths that were freshly written this call (empty if all
    were already up-to-date).

    Best-effort by design: any OSError is caught and logged so a
    locked / read-only data dir never crashes the launch path. The
    user's worst-case is a missed snapshot day, not an aborted app.

    `paths` defaults to every persisted user-data file (the resolved
    contents of `_USER_DATA_FILE_ATTRS` — kept in sync with the same
    table the pre-update snapshot consults, so the daily snapshot net
    can't drift behind the user's actual data footprint). Passing a
    custom list keeps the helper test-friendly without exposing fragile
    module globals to the test.
    """
    # L2 chokepoint: these snapshot copies land under the data dir, so gate
    # them like every other data-dir writer — an unsandboxed probe must not be
    # able to write into the user's real data dir through this helper. A no-op
    # in production (main() / the pytest fixture / the agent server authorise
    # writes before this ever runs).
    _refuse_unauthorized_write(data_dir, "daily snapshot")
    if paths is None:
        paths = []
        g = globals()
        for attr in g.get("_USER_DATA_FILE_ATTRS", ()):
            p = _state._resolve_data_attr_hook(attr)
            if isinstance(p, Path):
                paths.append(p)
    snap_dir = data_dir / _state._SNAPSHOT_DIR_NAME
    today = _date.today().isoformat()  # YYYY-MM-DD
    written: list[Path] = []
    try:
        snap_dir.mkdir(parents=True, exist_ok=True)
    except OSError:
        _log.warning("Could not create snapshot dir %s", snap_dir)
        return written
    for src in paths:
        if not src.exists():
            continue
        try:
            src_size = src.stat().st_size
        except OSError:
            continue
        if src_size == 0:
            continue
        # Refuse only a file so large that snapshotting it would evict
        # every other snapshot (see the cap's own comment). This is now
        # a genuinely rare case, so say so LOUDLY — the old flat 50 MB
        # cap skipped the user's most valuable file in silence
        # ([INV-167]). Checked BEFORE the read so an oversized file is
        # never pulled into memory.
        if src_size > _SNAPSHOT_FILE_SIZE_CAP:
            _log.warning(
                "NO daily snapshot for %s: %d bytes exceeds the "
                "per-file cap of %d (half the %d-byte snapshot-dir "
                "budget). Rollback for this file is limited to "
                "_safe_save_json's .bak rotation and the pre-update "
                "snapshots.",
                src.name, src_size, _SNAPSHOT_FILE_SIZE_CAP,
                _SNAPSHOT_TOTAL_SIZE_CAP,
            )
            continue
        dest = snap_dir / f"{src.stem}-{today}{src.suffix}"
        # A `.gz` sibling also counts as "today is already snapshotted".
        # Nothing here WRITES one — compression is deliberately not done
        # on the launch path — but a user who hand-compresses a snapshot
        # to reclaim disk must not have it silently re-copied (and must
        # not lose the pre-damage snapshot to a same-day relaunch).
        if dest.exists() or dest.with_name(dest.name + ".gz").exists():
            continue
        try:
            _atomic_write_bytes(dest, src.read_bytes())
            written.append(dest)
            _log.info(
                "Snapshotted %s → %s (%d bytes)",
                src.name, dest.name, dest.stat().st_size,
            )
        except OSError:
            _log.exception("Could not snapshot %s", src)
    _prune_old_snapshots(snap_dir)
    return written


def _prune_old_snapshots(snap_dir: Path,
                          retain_days: "int | None" = None) -> None:
    """Delete snapshots whose date stamp is older than
    `retain_days`. `retain_days` defaults to
    `_SNAPSHOT_RETENTION_DAYS`, but is read at call time so a test
    can monkeypatch the module constant. Best-effort — never raises."""
    import re as _re
    if retain_days is None:
        retain_days = _SNAPSHOT_RETENTION_DAYS
    today = _date.today()
    # Snapshot filenames look like `<stem>-YYYY-MM-DD.<ext>`. The
    # regex pulls the date out of the trailing portion of the stem so
    # arbitrary stems (with their own dashes) round-trip cleanly.
    date_re = _re.compile(r"-(\d{4})-(\d{2})-(\d{2})$")

    def _dated_stem(p: Path) -> str:
        """Stem with a trailing `.gz` peeled off first.

        `_snapshot_data_files` never writes `.gz` itself, but a user may
        hand-compress a snapshot to reclaim disk. `Path.stem` of
        `x-<date>.json.gz` still ends in `.json`, so the `$`-anchored
        regex would miss it — the file would look UNDATED, and undated
        files are deleted before any dated one by the aggregate cap
        below. That would quietly destroy the snapshots a user cared
        enough about to compress ([INV-167]).
        """
        name = p.name[:-3] if p.name.endswith(".gz") else p.name
        return Path(name).stem

    try:
        candidates = list(snap_dir.iterdir())
    except OSError:
        return
    for snap in candidates:
        if not snap.is_file():
            continue
        m = date_re.search(_dated_stem(snap))
        if m is None:
            continue
        try:
            snap_date = _date(int(m.group(1)), int(m.group(2)),
                               int(m.group(3)))
        except ValueError:
            continue
        if (today - snap_date).days > retain_days:
            try:
                snap.unlink()
            except OSError:
                _log.debug("Could not prune old snapshot %s", snap)

    # Aggregate-size cap: after date-based pruning, if the directory
    # still exceeds `_SNAPSHOT_TOTAL_SIZE_CAP`, drop the oldest
    # remaining snapshots until usage falls below the cap. Files are
    # ranked by their embedded date (oldest first); a file with an
    # unparseable date is treated as oldest and is removed before any
    # dated file. Pre-fix the snapshot dir had no aggregate ceiling —
    # 4 files × 49 MB × 30 days could grow to ~5.9 GB.
    try:
        survivors = [s for s in snap_dir.iterdir() if s.is_file()]
    except OSError:
        return
    total = 0
    sized: list[tuple[int, "_date | None", Path]] = []
    for s in survivors:
        try:
            sz = s.stat().st_size
        except OSError:
            continue
        total += sz
        m = date_re.search(_dated_stem(s))
        if m is None:
            sized.append((sz, None, s))
        else:
            try:
                sized.append(
                    (sz,
                     _date(int(m.group(1)), int(m.group(2)),
                            int(m.group(3))),
                     s),
                )
            except ValueError:
                sized.append((sz, None, s))
    if total <= _SNAPSHOT_TOTAL_SIZE_CAP:
        return
    # Oldest-first sort: undated files first, then ascending by date.
    sized.sort(key=lambda t: (t[1] is not None,
                                t[1] or _date.min))
    for sz, _d, snap in sized:
        if total <= _SNAPSHOT_TOTAL_SIZE_CAP:
            break
        try:
            snap.unlink()
            total -= sz
            _log.info(
                "Pruned snapshot %s (%d bytes) — aggregate cap",
                snap.name, sz,
            )
        except OSError:
            _log.debug("Could not prune snapshot %s", snap)


def _run_data_dir_housekeeping(data_dir: Path) -> None:
    """Reclaim disk left by the backup + spillover safety nets.

    For every user-data file: compress older timestamped backups
    (`_compress_old_backups`) then enforce the count + byte-cap retention
    (`_prune_backups`). Then bound the `lost_entries/` spillover
    directory (`_prune_lost_entries`).

    This is what reclaims the EXISTING residue — a user upgrading with
    4.7 GB of uncompressed backups + 1.5 GB of un-pruned spillover sees
    it shrink on the next launch; `_safe_save_json` keeps it bounded
    thereafter. Runs OFF the UI thread (the caller spawns it as a daemon)
    because compressing a 274 MB backup is CPU-heavy and would otherwise
    add a multi-second hang to launch. Best-effort throughout — every
    step swallows its own errors so a locked / read-only data dir never
    disrupts launch. `_USER_DATA_FILE_ATTRS` is read via globals() so
    this can be defined ahead of that table (same pattern as
    `_snapshot_data_files`)."""
    g = globals()
    for attr in g.get("_USER_DATA_FILE_ATTRS", ()):
        p = _state._resolve_data_attr_hook(attr)
        if not isinstance(p, Path):
            continue
        # Compress BEFORE pruning so the byte cap measures compressed
        # sizes and therefore retains more generations.
        try:
            _compress_old_backups(p)
        except Exception:
            _log.debug("housekeeping: compress failed for %s", p)
        try:
            _prune_backups(p)
        except Exception:
            _log.debug("housekeeping: backup prune failed for %s", p)
    try:
        _prune_lost_entries(data_dir / _state._LOST_ENTRIES_DIR_NAME)
    except Exception:
        _log.debug("housekeeping: lost-entries prune failed")
    # Reclaim disk from superseded/deleted plasmid-sequence blobs. Safe by
    # construction (quarantines, never deletes; aborts on unreadable
    # metadata; grace window protects in-flight saves) — see
    # `_gc_orphan_blobs`.
    try:
        _gc = getattr(_state, "_gc_orphan_blobs_hook", None)
        if _gc is not None:
            _gc()
    except Exception:
        _log.debug("housekeeping: blob GC failed")


def _fsync_path(path: Path) -> None:
    """Best-effort fsync of a file or directory. Silently no-ops on
    platforms / filesystems that can't fsync the target (Windows
    refuses fsync on directories; some networked filesystems return
    ENOSYS). Used by the pre-update snapshot writer so a power loss
    between copy + os.replace can't leave a half-written manifest.

    Mirrors the durability convention `_safe_save_json` follows for
    JSON registry writes."""
    try:
        flag = os.O_RDONLY
        # Some POSIX flavours need O_DIRECTORY for dir fsync; missing
        # on Windows so we guard with getattr.
        if path.is_dir():
            flag |= getattr(os, "O_DIRECTORY", 0)
        fd = os.open(str(path), flag)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
    except (OSError, AttributeError, ValueError):
        # fsync is durability insurance; never let it crash a backup.
        pass


# Well-known SHARED system directories that Master Delete's recursive rmtree of
# the pre-update backup dir must never target. The auto-derived path + the
# filesystem-root guard already refuse `/`, depth-1 dirs (whose parent IS the
# root), the home dir, and the data dir / its ancestors; this closes the
# depth-2+ gap (`/usr/local`, `/var/lib`, `/opt/x` would otherwise pass). A
# DEDICATED subdir under one of these (e.g. `/opt/splicecraft-backups`) is still
# allowed — only the shared dir itself is refused.
_UNSAFE_BACKUP_DIRS: "frozenset[str]" = frozenset(
    str(Path(p)) for p in (
        "/usr", "/usr/local", "/usr/bin", "/usr/sbin", "/usr/lib",
        "/usr/share", "/opt", "/etc", "/var", "/var/lib", "/var/log",
        "/var/tmp", "/bin", "/sbin", "/lib", "/lib64", "/boot", "/sys",
        "/proc", "/dev", "/run", "/root", "/home", "/tmp", "/mnt",
        "/media", "/srv",
        "/Applications", "/System", "/Library", "/Users",   # macOS
    )
)


def _resolve_pre_update_backup_dir(data_dir: "Path | None" = None) -> Path:
    """Where pre-update snapshots are stored. Sibling of `_state._DATA_DIR` by
    default so a recursive wipe of `_state._DATA_DIR` doesn't take the
    snapshots with it. Overridable via $SPLICECRAFT_UPDATE_BACKUP_DIR
    (also honoured by tests so they never write to the user's real
    home directory).

    Hardened: resolves both inputs to canonical absolute paths so
    relative `_state._DATA_DIR` values, trailing slashes, and intermediate
    symlinks don't subtly change where snapshots land. Refuses paths
    that already exist as a file (won't silently overwrite), or
    derive into a filesystem root (parent.parent == parent) which
    almost certainly indicates a misconfigured data dir.
    """
    override = os.environ.get("SPLICECRAFT_UPDATE_BACKUP_DIR", "").strip()
    if override:
        candidate = Path(override).expanduser()
        # Resolve, but tolerate the common "doesn't exist yet" case —
        # `Path.resolve()` is non-strict by default since 3.6.
        candidate = candidate.resolve()
        if candidate.exists() and not candidate.is_dir():
            raise OSError(
                f"$SPLICECRAFT_UPDATE_BACKUP_DIR points at a non-directory "
                f"({candidate}). Set it to a writable directory path."
            )
        # Master Delete recursively rmtrees this directory (it's a sibling-
        # backup target in `_master_delete_sibling_targets`). The auto-derived
        # path below is guarded against deriving into a filesystem root; the
        # override skipped that guard entirely, so a value like `$HOME` or the
        # data dir would be wiped on a Master Delete. Refuse the catastrophic
        # targets: a filesystem root, the home directory, and the data dir
        # itself or any ancestor of it.
        _parent = candidate.parent
        _data = _state._DATA_DIR.resolve()
        if (_parent == candidate or _parent.parent == _parent
                or candidate == Path.home().resolve()
                or candidate == _data or candidate in _data.parents
                or str(candidate) in _UNSAFE_BACKUP_DIRS):
            raise OSError(
                f"Refusing $SPLICECRAFT_UPDATE_BACKUP_DIR={str(candidate)!r}: "
                "it is a filesystem root, your home directory, or the data dir "
                "/ an ancestor of it — Master Delete would recursively wipe it. "
                "Point it at a dedicated backup directory."
            )
        return candidate
    base = (data_dir if data_dir is not None else _state._DATA_DIR).resolve()
    parent = base.parent
    # Defensive guard: refuse to write a snapshot into a filesystem
    # root or any path where the parent is the same as the path
    # itself. `parent.parent == parent` is the cross-platform "is
    # filesystem root" check (covers `/`, `C:\`, network shares).
    if parent == base or parent.parent == parent:
        raise OSError(
            f"Refusing to derive a backup location next to {base!r} — "
            "set $SPLICECRAFT_UPDATE_BACKUP_DIR to a writable directory."
        )
    candidate = parent / f"{base.name}-update-backups"
    if candidate.exists() and not candidate.is_dir():
        raise OSError(
            f"Default backup location {candidate} is not a directory. "
            "Set $SPLICECRAFT_UPDATE_BACKUP_DIR to override."
        )
    return candidate


def _iter_user_data_paths() -> "list[tuple[str, str, Path, str]]":
    """Yield `(attr_name, kind, path, name)` for every user-data file
    or directory that currently exists on disk.

    `kind` is `"file"` or `"dir"`. `name` is the basename used inside
    the snapshot. Reads the live module attribute each call so the
    autouse `_protect_user_data` test fixture's monkeypatched paths
    are picked up correctly.
    """
    out: list[tuple[str, str, Path, str]] = []
    for attr in _USER_DATA_FILE_ATTRS:
        p = _state._resolve_data_attr_hook(attr)
        if isinstance(p, Path) and p.is_file():
            out.append((attr, "file", p, p.name))
    for attr in _USER_DATA_DIR_ATTRS:
        p = _state._resolve_data_attr_hook(attr)
        if isinstance(p, Path) and p.is_dir():
            out.append((attr, "dir", p, p.name))
    return out


def _master_delete_file_targets() -> "list[Path]":
    """Every regular file under `_state._DATA_DIR` that Master Delete removes.

    Built from `_USER_DATA_FILE_ATTRS` + `_OPERATIONAL_FILE_ATTRS`
    plus each file's `.bak` / `.bak.<ts>[.N]` siblings (`_safe_save_json`
    writes both single-gen + timestamped backups; cf. invariant #31).
    Reads live module attributes so the autouse `_protect_user_data`
    fixture's monkeypatch flows through.
    """
    out: list[Path] = []
    for attr in _USER_DATA_FILE_ATTRS + _OPERATIONAL_FILE_ATTRS:
        p = _state._resolve_data_attr_hook(attr)
        if not isinstance(p, Path):
            continue
        out.append(p)
        # Single-generation backup.
        out.append(p.with_suffix(p.suffix + ".bak"))
        # Timestamped rotating backups + same-second collision bumps.
        parent = p.parent
        if parent.is_dir():
            try:
                for sib in parent.glob(p.name + ".bak.*"):
                    out.append(sib)
            except OSError:
                pass
    return out


def _master_delete_dir_targets() -> "list[Path]":
    """Every directory directly under `_state._DATA_DIR` that Master Delete
    removes.

    `_USER_DATA_DIR_ATTRS` covers crash_recovery / dna_originals /
    experiments / plugins. The four ad-hoc dirs created inline by
    other code paths (snapshots / lost_entries / clipboard /
    ui_snapshots) are appended here so the wipe is exhaustive.

    Excludes the logs directory — the active log file is held open
    by `RotatingFileHandler`; rotated backups are handled separately
    by `_master_delete_log_files`. Excludes `_state._DATA_DIR` itself so the
    process-held lockfile survives.
    """
    out: list[Path] = []
    seen: set[Path] = set()

    def _add(p: Path) -> None:
        try:
            key = p.resolve(strict=False)
        except OSError:
            key = p
        if key in seen:
            return
        seen.add(key)
        out.append(p)

    for attr in _USER_DATA_DIR_ATTRS:
        p = _state._resolve_data_attr_hook(attr)
        if isinstance(p, Path):
            _add(p)
    # Ad-hoc dirs created inline elsewhere (no module constant for
    # most). Resolve `_state._DATA_DIR` at call time so monkeypatched test
    # paths are honoured.
    data_root = _state._DATA_DIR
    if isinstance(data_root, Path):
        for name in ("snapshots", "lost_entries", "clipboard"):
            _add(data_root / name)
    ui = _state._UI_SNAPSHOTS_DIR
    if isinstance(ui, Path):
        _add(ui)
    return out


def _master_delete_sibling_targets() -> "list[Path]":
    """Pre-update backup directory (sibling of `_state._DATA_DIR`).

    Resolved via `_resolve_pre_update_backup_dir` so
    `$SPLICECRAFT_UPDATE_BACKUP_DIR` (set by tests) takes precedence.
    Returns an empty list if the resolver refuses (filesystem root
    edge case) — nothing to wipe there.
    """
    out: list[Path] = []
    try:
        cand = _resolve_pre_update_backup_dir()
    except OSError:
        return out
    if isinstance(cand, Path):
        out.append(cand)
    return out


def _master_delete_extra_root_files() -> "list[Path]":
    """Stragglers directly under `_state._DATA_DIR` that aren't covered by
    `_USER_DATA_FILE_ATTRS` / `_OPERATIONAL_FILE_ATTRS` but are still
    user-impacting state. As of 2026-05-20 this is just the
    `.migrated` marker from `_migrate_legacy_data`; future additions
    can be appended here so they're explicitly enumerated rather
    than relying on the residual-sweep fallback.
    """
    out: list[Path] = []
    data_root = _state._DATA_DIR
    if isinstance(data_root, Path):
        out.append(data_root / ".migrated")
    return out


def _master_delete_residual_paths() -> "tuple[list[Path], list[Path]]":
    """Return (files, dirs) under `_state._DATA_DIR` left over after the
    named-target wipe finished. Returns empty lists if the data dir
    no longer exists.

    `splicecraft.lock` and the `logs/` directory are excluded — the
    lockfile is held by the running process; the logs dir contains
    the active log handle that `RotatingFileHandler` writes into.
    Both are handled by other parts of `_perform_master_delete`.

    Designed to be re-runnable: a second sweep should find nothing.
    Used as the final defense-in-depth pass so any future code that
    writes a new file under `_state._DATA_DIR` (without registering it in
    `_USER_DATA_FILE_ATTRS` / `_OPERATIONAL_FILE_ATTRS`) still gets
    wiped on Master Delete.
    """
    files: list[Path] = []
    dirs: list[Path] = []
    data_root = _state._DATA_DIR
    if not isinstance(data_root, Path):
        return files, dirs
    if not data_root.is_dir():
        return files, dirs
    # Anchor to the active lockfile + log dir resolved values so a
    # symlink-mediated path mismatch can't bypass the skip.
    lock_path = data_root / "splicecraft.lock"
    log_dir = data_root / "logs"
    try:
        lock_resolved = lock_path.resolve(strict=False)
    except OSError:
        lock_resolved = lock_path
    try:
        log_resolved = log_dir.resolve(strict=False)
    except OSError:
        log_resolved = log_dir
    try:
        children = list(data_root.iterdir())
    except OSError:
        return files, dirs
    for child in children:
        try:
            child_resolved = child.resolve(strict=False)
        except OSError:
            child_resolved = child
        if child_resolved == lock_resolved:
            continue
        if child_resolved == log_resolved:
            continue
        try:
            if child.is_symlink():
                # Symlinks under DATA_DIR are unusual; remove them
                # via unlink (don't follow into shared filesystems).
                files.append(child)
            elif child.is_file():
                files.append(child)
            elif child.is_dir():
                dirs.append(child)
        except OSError:
            continue
    return files, dirs


def _master_delete_log_files() -> "list[Path]":
    """Rotated log backups (`splicecraft.log.1`, `.2`, …) ready for
    deletion. The currently-open `_LOG_PATH` is excluded — closing the
    file handle out from under the running `RotatingFileHandler`
    works on POSIX but raises `PermissionError` on Windows. Restart
    rolls the active log into a fresh file naturally.
    """
    out: list[Path] = []
    log_path = globals().get("_LOG_PATH")
    if not log_path:
        return out
    try:
        log_file = Path(log_path)
    except (TypeError, ValueError):
        return out
    log_dir = log_file.parent
    if not log_dir.is_dir():
        return out
    try:
        active_resolved = log_file.resolve(strict=False)
    except OSError:
        active_resolved = log_file
    try:
        children = list(log_dir.iterdir())
    except OSError:
        return out
    for child in children:
        try:
            if not child.is_file():
                continue
            if child.resolve(strict=False) == active_resolved:
                continue
            out.append(child)
        except OSError:
            continue
    return out


def _sha256_file(path: Path) -> str:
    """SHA-256 of a file, computed in 64 KB chunks. Used in the
    manifest so a partial / truncated copy in a snapshot can be
    detected at restore time."""
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(64 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _safe_snapshot_token(s: str, max_len: int = 32) -> str:
    """Sanitise a string for use inside a snapshot directory name.
    Strips path separators, leading dots, and anything outside a
    conservative ASCII set so the resulting filename is portable
    across Linux/Mac/Windows."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]", "_", s)
    cleaned = cleaned.lstrip(".") or "x"
    if _is_windows_reserved_stem(cleaned):
        cleaned = f"_{cleaned}"
    return cleaned[:max_len]


def _enforce_pre_update_retention(backup_dir: Path,
                                    keep: int = _PRE_UPDATE_SNAPSHOT_RETENTION
                                    ) -> None:
    """Prune oldest snapshots in `backup_dir` so at most `keep`
    remain. Best-effort: a failed delete is logged but never raised
    (the snapshot we just took is more important than reaping old
    ones).

    SAFETY: rmtree is restricted to subdirectories whose names match
    the strict regex `_PRE_UPDATE_NAME_RE`. Without that filter, a
    misconfigured (or symlinked) `backup_dir` pointing at `/` or `~`
    could match `bin`/`etc`/`home`/etc. against the previous loose
    `not name.startswith(staging-prefix)` check — and rmtree the
    user's entire system. The regex names cannot occur anywhere
    except in directories `_create_pre_update_snapshot` itself
    produced. This is enforced by
    `test_retention_only_rmtrees_snapshot_named_dirs`.

    A symlinked backup_dir is also refused outright: if the user (or
    an attacker) replaces the backup location with a symlink, we
    don't follow it for cleanup purposes.
    """
    if keep < 0:
        keep = 0
    try:
        if not backup_dir.is_dir():
            return
        # Refuse to walk into a symlinked backup_dir. `is_symlink`
        # works on POSIX + Windows (NTFS junctions read as symlinks).
        if backup_dir.is_symlink():
            _log.warning(
                "Refusing retention sweep: %s is a symlink", backup_dir
            )
            return
        # Resolve to canonical path and re-check it isn't a system
        # root. Belt + suspenders against `SPLICECRAFT_UPDATE_BACKUP_DIR=/`.
        resolved = backup_dir.resolve()
        if resolved.parent == resolved:
            _log.warning(
                "Refusing retention sweep: %s resolves to a filesystem root",
                backup_dir,
            )
            return
        snaps = []
        for p in backup_dir.iterdir():
            if not p.is_dir():
                continue
            if p.is_symlink():
                # Don't follow symlinks during cleanup either.
                continue
            if not _PRE_UPDATE_NAME_RE.match(p.name):
                # Foreign directory — not one of ours. Leave alone.
                continue
            try:
                snaps.append((p.stat().st_mtime, p))
            except OSError:
                continue
        snaps.sort(reverse=True)  # newest first
        for _, old in snaps[keep:]:
            try:
                shutil.rmtree(old, ignore_errors=False)
            except OSError as exc:
                _log.warning("Could not prune old snapshot %s: %s", old, exc)
    except OSError as exc:
        _log.warning("Pre-update retention sweep failed: %s", exc)


@_timed("op.create_pre_update_snapshot")
def _create_pre_update_snapshot(version_from: str,
                                 *,
                                 backup_dir: "Path | None" = None,
                                 retention: int = _PRE_UPDATE_SNAPSHOT_RETENTION,
                                 exclude_attrs: "frozenset[str] | None" = None,
                                 ) -> Path:
    """Atomically snapshot every user-data file + directory into a
    fresh subdirectory of `backup_dir`. Returns the path to the
    completed snapshot. Raises `OSError` (or a subclass) if anything
    goes wrong — callers MUST treat that as "abort the upgrade".

    Atomicity: copying happens under a hidden staging directory on
    the same filesystem. The staging directory is renamed to its
    final name (`os.replace`) only after every file copy and the
    manifest write succeed. A partially-written staging directory is
    cleaned up before re-raising. The temp-dir + atomic-rename pattern
    mirrors `_safe_save_json`.

    The snapshot lives OUTSIDE `_state._DATA_DIR` (default sibling) so a
    bug in a new version that recursively wipes `_state._DATA_DIR` cannot
    destroy the recovery copy.
    """
    import datetime as _dt
    import secrets as _secrets

    if backup_dir is None:
        backup_dir = _resolve_pre_update_backup_dir()
    backup_dir.mkdir(parents=True, exist_ok=True)
    # Hardened: refuse to write into a symlinked backup_dir. Combined
    # with the symlink-aware retention sweep, this means even a
    # malicious environment variable can't trick us into writing to —
    # or later deleting from — a system path.
    if backup_dir.is_symlink():
        raise OSError(
            f"Refusing to use {backup_dir} as backup location: "
            "the path is a symlink."
        )

    ts = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    rand = _secrets.token_hex(4)
    safe_ver = _safe_snapshot_token(version_from)
    staging = backup_dir / f"{_PRE_UPDATE_STAGING_PREFIX}{ts}-{rand}"
    final_name = f"{ts}-{rand}__from-{safe_ver}"
    final_path = backup_dir / final_name
    # Even with random hex, defensively bump the suffix if the final
    # name already exists (would only happen on a clock-rewind + hash
    # collision, but cheap insurance).
    bump = 0
    while final_path.exists():
        bump += 1
        final_path = backup_dir / f"{final_name}.{bump}"

    try:
        staging.mkdir(parents=True, exist_ok=False)
        sources = _iter_user_data_paths()
        if exclude_attrs:
            # Migrate-export omits large, re-derivable trees (the
            # downloaded HMM databases — their catalog json IS still
            # copied, so the new machine re-downloads them on demand).
            # Every existing caller passes nothing, so the pre-update
            # snapshot footprint is byte-for-byte unchanged.
            sources = [s for s in sources if s[0] not in exclude_attrs]
        # 2026-05-27 (audit-2 H4): hold `_state._cache_lock` across the
        # entire snapshot copy loop. Pre-fix `shutil.copy2` ran
        # without the lock — a concurrent `_save_library` /
        # `_save_collections` could `os.replace` the source mid-
        # read, producing a Frankenstein with new data after old.
        # The manifest then sha256'd the destination, locking in
        # the corruption. RLock allows the inner `_safe_save_json`
        # calls (none expected during snapshot, but defensive)
        # to re-enter freely.
        manifest: dict = {
            "schema_version": _PRE_UPDATE_SCHEMA_VERSION,
            "from_version": version_from,
            "timestamp": ts,
            "data_dir": str(_state._DATA_DIR),
            # Capture the runtime Python + platform that wrote the
            # snapshot. Restore-time consumers (and `--list-snapshots`
            # display) use this to flag cross-platform / cross-Python
            # restores that may have unexpected side-effects (e.g.
            # binary-pickle state we may add later won't survive a
            # 3.10 → 3.14 jump). Additive — older readers ignore via
            # `.get(...)` defaults so this stays back-compatible.
            #
            # `_RUNTIME_PLATFORM` is computed once at module import so
            # we don't re-call `platform.platform()` every snapshot —
            # on some OSes that helper shells out via `subprocess`,
            # which conflicts with tests that monkeypatch `subprocess.run`
            # to capture the upgrade command.
            "from_python_version": "{}.{}.{}".format(
                sys.version_info[0], sys.version_info[1], sys.version_info[2]
            ),
            "from_platform": _RUNTIME_PLATFORM,
            "files": [],
            "directories": [],
        }
        # H4: take the cache lock around the entire copy loop so the
        # source files are quiesced against concurrent `_save_library`
        # / `_save_collections`. Pre-fix `shutil.copy2` could capture
        # a mid-rename inode → manifest sha256 locked in a half-old/
        # half-new "Frankenstein" file. RLock allows the inner
        # `_atomic_write_text(manifest_path, ...)` to re-enter freely.
        with _state._cache_lock:
            for attr, kind, src, name in sources:
                dst = staging / name
                if kind == "file":
                    shutil.copy2(src, dst)
                    # Durability: fsync the copy before the atomic
                    # rename so a power loss between copy and replace
                    # can't leave a half-written file with a manifest
                    # that claims it's complete. Mirrors `_safe_save_json`.
                    _fsync_path(dst)
                    # Hash the destination (post-copy) so a silent
                    # corrupt would be visible in the manifest.
                    manifest["files"].append({
                        "attr": attr,
                        "name": name,
                        "size": dst.stat().st_size,
                        "sha256": _sha256_file(dst),
                    })
                else:  # kind == "dir"
                    shutil.copytree(src, dst)
                    # Best-effort fsync of every file inside the
                    # copied tree. We don't error out if a single
                    # file can't be fsynced (read-only mounts,
                    # network FS quirks).
                    try:
                        for f in dst.rglob("*"):
                            if f.is_file():
                                _fsync_path(f)
                    except OSError:
                        pass
                    # Counting + size for directories is informational
                    # only; restoring uses the directory tree as-is.
                    file_count = sum(
                        1 for _ in dst.rglob("*") if _.is_file()
                    )
                    manifest["directories"].append({
                        "attr": attr,
                        "name": name,
                        "file_count": file_count,
                    })
        # Manifest goes in last so its presence is itself a "snapshot
        # is complete" marker. Restore code refuses to operate on a
        # snapshot dir without a manifest. Sweep #35 (2026-05-26):
        # route through `_atomic_write_text` so a crash mid-write
        # leaves either the previous (non-existent) state OR a
        # complete manifest — never a truncated one that passes the
        # "manifest exists" gate but fails JSON parsing during
        # restore. `_atomic_write_text` does its own fsync of the
        # tmp file before rename, so we no longer need the
        # standalone `_fsync_path(manifest_path)` call.
        manifest_path = staging / _PRE_UPDATE_MANIFEST_NAME
        _atomic_write_text(
            manifest_path,
            json.dumps(manifest, indent=2),
        )
        # Fsync the staging directory so the new entries (manifest +
        # copies) are durable in the directory's metadata before we
        # rename it to its final visible name. POSIX-only; the helper
        # silently no-ops on Windows.
        _fsync_path(staging)
        # Atomic rename — single filesystem op on POSIX, transactional
        # MoveFileEx on Windows. After this point the snapshot is
        # visible at its final name; before it, only `.tmp-...` exists.
        os.replace(staging, final_path)
        # Fsync the parent so the rename itself is durable.
        _fsync_path(backup_dir)
    except (OSError, shutil.Error):
        # Best-effort cleanup of the partial staging dir.
        try:
            if staging.exists():
                shutil.rmtree(staging, ignore_errors=True)
        except OSError:
            pass
        raise

    # Retention: prune oldest snapshots beyond N. Best-effort — never
    # let a retention failure invalidate the snapshot we just took.
    _enforce_pre_update_retention(backup_dir, keep=retention)
    return final_path


def _list_pre_update_snapshots(backup_dir: "Path | None" = None
                                ) -> "list[dict]":
    """Return a list of snapshot summaries, newest first. Each entry:
        {id, path, mtime, from_version, n_files, n_dirs, total_size}
    `id` is the directory's basename (what the user passes to
    --restore-pre-update). Snapshots without a manifest are skipped
    (those are partially-written or hand-corrupted)."""
    if backup_dir is None:
        try:
            backup_dir = _resolve_pre_update_backup_dir()
        except OSError:
            return []
    if not backup_dir.is_dir():
        return []
    out: list[dict] = []
    for p in backup_dir.iterdir():
        if not p.is_dir():
            continue
        if p.name.startswith(_PRE_UPDATE_STAGING_PREFIX):
            continue
        manifest = p / _PRE_UPDATE_MANIFEST_NAME
        if not manifest.is_file():
            continue
        ok_size, _reason = _safe_file_size_check(
            manifest, _PRE_UPDATE_MANIFEST_MAX_BYTES, "pre-update manifest",
        )
        if not ok_size:
            _log.warning("pre-update snapshot %s: manifest rejected (%s)",
                         p.name, _reason)
            continue
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError, UnicodeDecodeError):
            continue
        try:
            mtime = p.stat().st_mtime
        except OSError:
            continue
        total = 0
        try:
            for f in p.rglob("*"):
                if f.is_file():
                    total += f.stat().st_size
        except OSError:
            pass
        out.append({
            "id": p.name,
            "path": p,
            "mtime": mtime,
            "from_version": data.get("from_version", "?"),
            # Newer fields — older snapshots return defaults so the
            # listing stays backward-compatible.
            "from_python_version": data.get("from_python_version", "?"),
            "from_platform": data.get("from_platform", "?"),
            "schema_version": data.get("schema_version", 1),
            "n_files": len(data.get("files", [])) if isinstance(data.get("files"), list) else 0,
            "n_dirs": len(data.get("directories", [])) if isinstance(data.get("directories"), list) else 0,
            "total_size": total,
        })
    out.sort(key=lambda e: e["mtime"], reverse=True)
    return out


def _validate_snapshot_member(snap_path: Path, name: str) -> "Path | None":
    """Return the absolute path of a manifest-named snapshot member,
    or None if the name would escape the snapshot directory (path
    traversal) or contain a separator. Pure validation — never
    touches the filesystem.

    A manifest is JSON, and JSON is plaintext: anyone with write
    access to the backup directory could craft `name: "../../foo"`
    to read or overwrite arbitrary files. This guard reduces every
    `name` to a single basename and verifies the resolved path is
    still under the snapshot directory.
    """
    if not isinstance(name, str) or not name:
        return None
    # Reject anything that smells like a path component vs. a basename.
    if "/" in name or "\\" in name or name in (".", ".."):
        return None
    candidate = (snap_path / name).resolve()
    snap_resolved = snap_path.resolve()
    try:
        # Python 3.9+ has Path.is_relative_to. We're on 3.10+ per
        # pyproject; safe to use directly.
        if not candidate.is_relative_to(snap_resolved):
            return None
    except (AttributeError, ValueError):
        return None
    return candidate


@_timed("op.restore_pre_update_snapshot")
def _restore_pre_update_snapshot(snap_id_or_path: "str | Path",
                                  *,
                                  backup_dir: "Path | None" = None,
                                  ) -> dict:
    """Restore a pre-update snapshot. Each canonical JSON is replaced
    via the temp-file + os.replace pattern (same as `_safe_save_json`).
    Sub-directories are restored by atomically renaming the current
    directory aside, copying the snapshot in, then removing the
    stash. A pre-restore snapshot is taken first so the user can
    undo a bad restore.

    Hardening (sacred — do not relax without a corresponding test):
      * `schema_version` must be ≤ `_PRE_UPDATE_SCHEMA_VERSION`. A
        higher version means the snapshot was written by a newer
        SpliceCraft we don't know how to read; refuse rather than
        silently restore an unfamiliar shape.
      * Each `attr` in the manifest must be in `_USER_DATA_FILE_ATTRS`
        or `_USER_DATA_DIR_ATTRS`. A tampered manifest cannot target
        arbitrary `_*_FILE` attributes (e.g. agent token, install
        marker) for overwrite.
      * Each `name` must be basename-only and resolve under the
        snapshot directory. Blocks path-traversal reads.
      * Each file's SHA-256 (per the manifest) is recomputed on the
        copy BEFORE the atomic os.replace overwrites the user's
        current file. A bit-rotted snapshot cannot silently corrupt
        good live data.
      * Directory restore rolls back fully on partial copytree:
        delete the partial target, then rename stash back.

    Returns a summary dict:
        {pre_restore_snapshot, restored_files, restored_dirs,
         failed: [(name, reason), ...]}
    Raises FileNotFoundError if the snapshot id doesn't exist.
    Raises ValueError on unsupported schema version.
    """
    # L2 chokepoint: this OVERWRITES and RMTREES live user data. Gate it like
    # every other data-dir writer/deleter so an unsandboxed probe can't be
    # tricked into clobbering the real data dir — the exact disaster the
    # chokepoint exists to prevent. No-op in production (main() / pytest fixture
    # / agent server authorise first); the pre-restore snapshot below keeps an
    # authorised restore undoable.
    _refuse_unauthorized_write(_state._DATA_DIR, "pre-update restore")
    _refuse_unauthorized_delete(_state._DATA_DIR, "pre-update restore")
    if backup_dir is None:
        backup_dir = _resolve_pre_update_backup_dir()
    if isinstance(snap_id_or_path, Path):
        snap_path = snap_id_or_path
    else:
        # Allow either a bare snapshot id (basename) or 'latest'.
        snap_id = str(snap_id_or_path).strip()
        if snap_id.lower() == "latest":
            snaps = _list_pre_update_snapshots(backup_dir)
            if not snaps:
                raise FileNotFoundError(
                    f"No pre-update snapshots in {backup_dir}"
                )
            snap_path = snaps[0]["path"]
        else:
            snap_path = backup_dir / snap_id
    if not snap_path.is_dir():
        raise FileNotFoundError(f"snapshot not found: {snap_path}")
    manifest_path = snap_path / _PRE_UPDATE_MANIFEST_NAME
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"snapshot is incomplete (no manifest): {snap_path}"
        )
    ok_size, reason = _safe_file_size_check(
        manifest_path, _PRE_UPDATE_MANIFEST_MAX_BYTES, "pre-update manifest",
    )
    if not ok_size:
        raise ValueError(f"manifest rejected: {reason}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError(f"manifest is not a JSON object: {manifest_path}")

    # Schema-version negotiation: forward-compat refuses anything
    # newer than we understand (would require migration logic we
    # don't yet have). Older schemas load with `.get(...)` defaults.
    snap_schema = manifest.get("schema_version", 1)
    try:
        snap_schema_int = int(snap_schema)
    except (TypeError, ValueError):
        raise ValueError(
            f"manifest schema_version is not an int: {snap_schema!r}"
        )
    if snap_schema_int > _PRE_UPDATE_SCHEMA_VERSION:
        raise ValueError(
            f"snapshot {snap_path.name} was created by a newer SpliceCraft "
            f"(manifest schema_version={snap_schema_int}, this build supports "
            f"≤{_PRE_UPDATE_SCHEMA_VERSION}). Upgrade SpliceCraft and retry, "
            "or use the matching version to restore."
        )

    # Whitelists: the manifest CANNOT name an attribute or path
    # outside our published user-data set. Tamper-resistance.
    file_attrs = set(_USER_DATA_FILE_ATTRS)
    dir_attrs = set(_USER_DATA_DIR_ATTRS)

    # Pre-restore snapshot — gives the user one click of undo if the
    # restored data turns out to be older than what they had.
    pre_restore = _create_pre_update_snapshot(
        f"{_state._sc_version}-pre-restore",
        backup_dir=backup_dir,
    )

    summary: dict = {
        "pre_restore_snapshot": str(pre_restore),
        "restored_files": [],
        "restored_dirs": [],
        "failed": [],
    }

    files_list = manifest.get("files")
    if not isinstance(files_list, list):
        files_list = []
    for entry in files_list:
        if not isinstance(entry, dict):
            continue
        attr = entry.get("attr", "")
        name = entry.get("name", "")
        # Hardening #1: attr must be in the published user-data set.
        if attr not in file_attrs:
            summary["failed"].append(
                (name or "(no name)",
                 f"manifest attr {attr!r} is not a recognised user-data file")
            )
            continue
        # Hardening #2: name must be a safe basename inside snap_path.
        src = _validate_snapshot_member(snap_path, name)
        if src is None or not src.is_file():
            summary["failed"].append(
                (name, "snapshot member name is invalid or src absent")
            )
            continue
        target = _state._resolve_data_attr_hook(attr)
        if not isinstance(target, Path):
            # Production code never has this happen (attrs are module
            # constants), but defending anyway: skip the entry rather
            # than write to whatever non-Path object happens to live
            # at that name in globals().
            summary["failed"].append(
                (name, f"runtime attr {attr!r} is not a Path")
            )
            continue
        # Initialised pre-try so the except cleanup can see it even
        # if the OSError fires on mkdir before the staging path is
        # assigned (otherwise pyright flags `tmp` as possibly unbound).
        # The staging name is randomised via `mkstemp` so concurrent
        # restore attempts (UI + agent API in the same process) can't
        # collide on a deterministic `<target>.restoring` path and
        # truncate each other mid-copy.
        tmp: "Path | None" = None
        tmp_fd: "int | None" = None
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            import tempfile as _tempfile
            tmp_fd, tmp_str = _tempfile.mkstemp(
                prefix=f".{target.name}.", suffix=".restoring",
                dir=str(target.parent),
            )
            os.close(tmp_fd)
            tmp_fd = None
            tmp = Path(tmp_str)
            shutil.copy2(src, tmp)
            _fsync_path(tmp)
            # Hardening #3: verify SHA-256 against the manifest BEFORE
            # the atomic overwrite of the live file. A bit-rotted /
            # tampered snapshot cannot silently corrupt good live data.
            # The sha256 entry is MANDATORY (one of invariant #39's
            # sacred-four checks); pre-0.8.9 silently skipped verify
            # when the field was missing/empty, which a tampered manifest
            # in the user-writable backup dir could exploit. Refuse the
            # restore rather than fall through to a blind rename.
            expected = entry.get("sha256")
            if not isinstance(expected, str) or not expected:
                try:
                    tmp.unlink()
                except OSError:
                    pass
                summary["failed"].append(
                    (name,
                     "manifest entry is missing the mandatory sha256 "
                     "field — refusing restore (snapshot may be tampered "
                     "or written by an unsupported tool)")
                )
                continue
            actual = _sha256_file(tmp)
            if actual != expected:
                # Discard the staged copy and skip the rename. The
                # user's live file remains intact.
                try:
                    tmp.unlink()
                except OSError:
                    pass
                summary["failed"].append(
                    (name,
                     f"sha256 mismatch (snapshot corrupted): "
                     f"expected {expected[:16]}…, got {actual[:16]}…")
                )
                continue
            os.replace(tmp, target)
            _fsync_path(target.parent)
            summary["restored_files"].append(name)
        except (OSError, shutil.Error) as exc:
            # Best-effort: clean up any half-written staging file.
            try:
                if tmp is not None and tmp.exists():
                    tmp.unlink()
            except OSError:
                pass
            summary["failed"].append((name, str(exc)))

    dirs_list = manifest.get("directories")
    if not isinstance(dirs_list, list):
        dirs_list = []
    for entry in dirs_list:
        if not isinstance(entry, dict):
            continue
        attr = entry.get("attr", "")
        name = entry.get("name", "")
        if attr not in dir_attrs:
            summary["failed"].append(
                (name or "(no name)",
                 f"manifest attr {attr!r} is not a recognised user-data dir")
            )
            continue
        src = _validate_snapshot_member(snap_path, name)
        if src is None or not src.is_dir():
            summary["failed"].append(
                (name, "snapshot member name is invalid or src absent")
            )
            continue
        target = _state._resolve_data_attr_hook(attr)
        if not isinstance(target, Path):
            summary["failed"].append(
                (name, f"runtime attr {attr!r} is not a Path")
            )
            continue
        # Sweep #25 (2026-05-23): cap total snapshot-dir restore size
        # at 5 GB. Pre-fix the dir branch had no size verification
        # (unlike the file branch's mandatory sha256 check), so a
        # tampered manifest could declare `n_files=0` for a multi-GB
        # crash_recovery dir and exhaust user disk on restore. The
        # snapshot dir lives under the sibling-update-backups path
        # which has looser perms than `_state._DATA_DIR`, so a local-
        # attacker-with-write to that parent has a plausible attack.
        _SNAPSHOT_DIR_RESTORE_MAX_BYTES = 5 * 1024 * 1024 * 1024
        try:
            total_size = 0
            for f in src.rglob("*"):
                if f.is_file():
                    total_size += f.stat().st_size
                    if total_size > _SNAPSHOT_DIR_RESTORE_MAX_BYTES:
                        break
            if total_size > _SNAPSHOT_DIR_RESTORE_MAX_BYTES:
                summary["failed"].append(
                    (name, f"snapshot dir exceeds 5 GB restore cap "
                     f"({total_size:,} bytes); refused"),
                )
                continue
        except OSError as exc:
            summary["failed"].append(
                (name, f"size-check failed: {exc}"),
            )
            continue
        # Randomised staging (mirrors the file branch's `mkstemp`): move the
        # live dir INTO a UNIQUE `mkdtemp` stash so two concurrent restores
        # (UI + agent API in the same process) can't collide on a deterministic
        # `.restoring-old` path and rmtree each other's only copy — and a
        # leftover stash from a killed prior restore is never blindly deleted.
        # The pre-restore snapshot taken above is the ultimate backstop.
        import tempfile as _tempfile
        stash_parent: "Path | None" = None
        stashed: "Path | None" = None
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                stash_parent = Path(_tempfile.mkdtemp(
                    prefix=f".{target.name}.restoring-",
                    dir=str(target.parent)))
                stashed = stash_parent / target.name
                os.rename(target, stashed)
            shutil.copytree(src, target)
            _fsync_path(target.parent)
            if stash_parent is not None:
                shutil.rmtree(stash_parent, ignore_errors=True)
            summary["restored_dirs"].append(name)
        except (OSError, shutil.Error) as exc:
            # Rollback (handles a partial copytree): remove the partial target,
            # move the stashed original back, then clean up the stash parent.
            # Only touch `target` when we actually stashed a copy — never
            # destroy the live dir when there was nothing to stash.
            try:
                if stashed is not None and stashed.exists():
                    if target.exists():
                        shutil.rmtree(target, ignore_errors=True)
                    os.rename(stashed, target)
                if stash_parent is not None and stash_parent.exists():
                    shutil.rmtree(stash_parent, ignore_errors=True)
            except OSError:
                pass
            summary["failed"].append((name, str(exc)))

    return summary


def _export_migrate_archive(dest_path: "str | Path", *,
                            include_hmm: bool = False) -> dict:
    """Package ALL user data into one portable, compressed ``.zip`` at
    ``dest_path``, atomically. Returns a summary dict
    ``{path, bytes, n_files, n_dirs, included_hmm}``.

    The archive wraps a pre-update-style snapshot (sha256 manifest of
    every registry file + the content-addressed blob store + the
    construction history embedded in the data) under
    ``_MIGRATE_ARCHIVE_SUBDIR/``, with a top-level ``_MIGRATE_MARKER_NAME``
    identifying the format. It is written to a temp file on the
    destination filesystem and ``os.replace``-d into place, so a crash
    mid-write never leaves a half-archive at the target.

    READ-ONLY with respect to the data dir — it only copies out (the
    throwaway snapshot lives in the OS tempdir and is removed here). The
    HMM databases are omitted unless ``include_hmm`` (large, re-derivable
    from the always-included catalog). Raises OSError on failure.
    """
    import zipfile as _zipfile
    import tempfile as _tempfile
    import datetime as _dt
    dest = Path(dest_path).expanduser()
    dest.parent.mkdir(parents=True, exist_ok=True)
    exclude = frozenset() if include_hmm else _MIGRATE_DEFAULT_EXCLUDE_ATTRS

    snap_parent = Path(_tempfile.mkdtemp(prefix="sc-migrate-export-"))
    tmp_zip: "Path | None" = None
    try:
        snap = _create_pre_update_snapshot(
            f"migrate-export-{_state._sc_version}",
            backup_dir=snap_parent,
            retention=10 ** 9,            # throwaway dir — never prune it
            exclude_attrs=exclude,
        )
        # Registry counts for the user-facing summary (from the manifest).
        try:
            _sm = json.loads(
                (snap / _PRE_UPDATE_MANIFEST_NAME).read_text(encoding="utf-8"))
            n_files = len(_sm.get("files", []) or [])
            n_dirs = len(_sm.get("directories", []) or [])
        except (OSError, ValueError):
            n_files = n_dirs = 0
        marker = {
            "format": "splicecraft-migrate",
            "format_version": _MIGRATE_FORMAT_VERSION,
            "app_version": _state._sc_version,
            "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "snapshot_subdir": _MIGRATE_ARCHIVE_SUBDIR,
            "included_hmm": bool(include_hmm),
        }
        fd, tmp_str = _tempfile.mkstemp(
            prefix=f".{dest.name}.", suffix=".scmig-tmp", dir=str(dest.parent))
        os.close(fd)
        tmp_zip = Path(tmp_str)
        with _zipfile.ZipFile(tmp_zip, "w",
                              compression=_zipfile.ZIP_DEFLATED,
                              compresslevel=6, allowZip64=True) as zf:
            zf.writestr(_MIGRATE_MARKER_NAME, json.dumps(marker, indent=2))
            for p in sorted(snap.rglob("*")):
                if p.is_file():
                    arc = (f"{_MIGRATE_ARCHIVE_SUBDIR}/"
                           + p.relative_to(snap).as_posix())
                    zf.write(p, arcname=arc)
        # Durability: fsync the temp zip, atomic rename, fsync the dir.
        with open(tmp_zip, "rb") as _f:
            os.fsync(_f.fileno())
        os.replace(tmp_zip, dest)
        tmp_zip = None
        _fsync_path(dest.parent)
    finally:
        shutil.rmtree(snap_parent, ignore_errors=True)
        if tmp_zip is not None:
            try:
                tmp_zip.unlink()
            except OSError:
                pass
    size = dest.stat().st_size
    _log_event("migrate.export.ok", bytes=size, n_files=n_files,
               n_dirs=n_dirs, included_hmm=bool(include_hmm))
    return {"path": str(dest), "bytes": size, "n_files": n_files,
            "n_dirs": n_dirs, "included_hmm": bool(include_hmm)}


def _import_migrate_archive(zip_path: "str | Path") -> dict:
    """Import a migrate archive produced by `_export_migrate_archive`,
    atomically replacing the current user data, and return the
    `_restore_pre_update_snapshot` summary (which carries
    ``pre_restore_snapshot`` — the automatic, reversible rollback point)
    plus the archive ``marker``.

    Bulletproof, in this order:
      1. Validate the top-level marker + ``format_version`` (refuse a
         newer/unknown format rather than guess at its shape).
      2. Refuse path-traversal / absolute / drive-letter member names and
         cap the total uncompressed size (zip-bomb guard) BEFORE writing a
         single byte to disk.
      3. Extract ONLY the wrapped snapshot to a temp staging dir, then
         hand it to `_restore_pre_update_snapshot` — which snapshots the
         CURRENT data first (so the import is undoable) and verifies every
         file's sha256 before the atomic per-file replace. A mismatched
         file is skipped, never half-written.

    Raises ValueError on a malformed / unsupported / unsafe archive,
    FileNotFoundError when the path is missing, OSError on IO failure.
    The caller is responsible for busting in-memory caches afterwards
    (the on-disk data has changed underneath them).
    """
    # L2 chokepoint: this replaces live user data (its own staging plus the
    # _restore_pre_update_snapshot call below). Gate at entry so an unsandboxed
    # probe can't import an archive over the real data dir. No-op once writes
    # are authorised (main() / pytest fixture / agent server).
    _refuse_unauthorized_write(_state._DATA_DIR, "migrate-archive import")
    import zipfile as _zipfile
    import tempfile as _tempfile
    zp = Path(zip_path).expanduser()
    if not zp.is_file():
        raise FileNotFoundError(f"migrate archive not found: {zp}")
    if not _zipfile.is_zipfile(zp):
        raise ValueError(f"{zp.name!r} is not a valid .zip archive.")
    staging = Path(_tempfile.mkdtemp(prefix="sc-migrate-import-"))
    staging_real = str(staging.resolve())
    try:
        with _zipfile.ZipFile(zp, "r") as zf:
            try:
                raw = zf.read(_MIGRATE_MARKER_NAME)
            except KeyError:
                raise ValueError(
                    "This .zip is not a SpliceCraft migrate archive "
                    f"(missing {_MIGRATE_MARKER_NAME}).") from None
            try:
                marker = json.loads(raw.decode("utf-8"))
            except (json.JSONDecodeError, UnicodeDecodeError):
                raise ValueError(
                    "The migrate archive's marker is corrupt.") from None
            if not (isinstance(marker, dict)
                    and marker.get("format") == "splicecraft-migrate"):
                raise ValueError(
                    "This .zip is not a SpliceCraft migrate archive.")
            fmt = marker.get("format_version")
            if not isinstance(fmt, int) or fmt > _MIGRATE_FORMAT_VERSION:
                raise ValueError(
                    f"This archive's format (v{fmt}) is newer than this "
                    f"SpliceCraft supports (≤ v{_MIGRATE_FORMAT_VERSION}). "
                    "Update SpliceCraft and try again.")
            subdir = marker.get("snapshot_subdir") or _MIGRATE_ARCHIVE_SUBDIR
            if (not isinstance(subdir, str) or "/" in subdir
                    or "\\" in subdir or subdir in ("", ".", "..")):
                raise ValueError("The migrate archive's marker is malformed.")

            infos = zf.infolist()
            if len(infos) > _MIGRATE_MAX_MEMBERS:
                raise ValueError("Migrate archive has too many members.")
            total = 0
            for info in infos:
                name = info.filename
                parts = Path(name).parts
                if (name.startswith("/") or name.startswith("\\")
                        or ".." in parts
                        or (len(name) >= 2 and name[1] == ":")):
                    raise ValueError(
                        f"Refusing unsafe member path in archive: {name!r}")
                total += int(getattr(info, "file_size", 0) or 0)
                if total > _MIGRATE_MAX_TOTAL_UNCOMPRESSED:
                    raise ValueError(
                        "Migrate archive expands beyond the safety cap "
                        "(possible zip-bomb) — refusing to extract.")
            prefix = subdir + "/"
            members = [i for i in infos
                       if i.filename.startswith(prefix) and not i.is_dir()]
            if not members:
                raise ValueError("Migrate archive contains no data payload.")
            # Bomb guard #2: the pre-loop tally above trusted the central-
            # directory `file_size`; a member that UNDER-claims its size would
            # slip past it and then stream unbounded data here. Cap the ACTUAL
            # bytes written across all members and abort the instant they
            # exceed the same ceiling. (`staging` is rmtree'd in `finally`, so
            # a tripped cap leaves no partial extraction behind.)
            written = 0
            for info in members:
                rel = info.filename[len(prefix):]
                out = staging / rel
                if not str(out.resolve()).startswith(staging_real):
                    raise ValueError(f"Refusing unsafe extract target: {rel!r}")
                out.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(info, "r") as src, open(out, "wb") as dst:
                    while True:
                        chunk = src.read(1 << 20)
                        if not chunk:
                            break
                        written += len(chunk)
                        if written > _MIGRATE_MAX_TOTAL_UNCOMPRESSED:
                            raise ValueError(
                                "Migrate archive expands beyond the safety cap "
                                "(possible zip-bomb) — refusing to extract.")
                        dst.write(chunk)
        # The extracted dir IS a pre-update snapshot (manifest.json at its
        # root). Hand it to the battle-tested restore: it snapshots the
        # CURRENT data first (reversible) and sha256-verifies every file
        # before the atomic replace.
        summary = _restore_pre_update_snapshot(staging)
        summary["marker"] = marker
        _log_event(
            "migrate.import.ok",
            restored=len(summary.get("restored_files", []) or []),
            dirs=len(summary.get("restored_dirs", []) or []),
            failed=len(summary.get("failed", []) or []),
        )
        return summary
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _format_pre_update_snapshot_table(snaps: "list[dict]") -> str:
    """Render `_list_pre_update_snapshots()` output as a human-readable
    table for printing to the terminal. Returns the empty string for
    no snapshots so the caller can compose its own "no snapshots
    found" message."""
    import datetime as _dt
    if not snaps:
        return ""
    lines = ["    Snapshot ID                                 From       Files  Dirs   Size",
             "    " + "─" * 80]
    for s in snaps:
        when = _dt.datetime.fromtimestamp(s["mtime"]).strftime("%Y-%m-%d %H:%M")
        size_kb = s["total_size"] / 1024.0
        if size_kb < 1024:
            size_str = f"{size_kb:7.1f} KB"
        else:
            size_str = f"{size_kb / 1024.0:7.1f} MB"
        # Truncate long ids only in display; they remain copy-pasteable
        # because we print the full id below the table.
        sid_disp = s["id"] if len(s["id"]) <= 44 else s["id"][:41] + "…"
        lines.append(
            f"    {sid_disp:<44} {s['from_version']:<10} "
            f"{s['n_files']:>4}  {s['n_dirs']:>4}  {size_str:>10}    "
            f"({when})"
        )
    return "\n".join(lines)


# ── Settings -> Restore-from-backup (the per-file recovery UI's engine; distinct
# from the pre-update snapshot restore above). Discovers every recoverable copy of
# a target JSON across the four storage tiers + applies one back to the live file
# (chokepoint-guarded via persistence). `_AGENT_BACKUP_LABELS` maps the agent /
# GUI friendly backup-target names to their _state attr; the hub GUI restore reads
# it via the re-export. `_resolve_backup_label` resolves attrs via the _state hook.
def _list_recoverable_backups(target_path: Path) -> "list[dict]":
    """Return every recoverable copy of `target_path` across the four
    storage tiers, sorted newest first.

    Each entry: ``{kind, source_path, n_entries, mtime_str}`` where
    ``kind`` is one of:
      * ``"legacy_bak"``    — single-generation `<file>.bak`
      * ``"rotating_bak"``  — timestamped `<file>.bak.YYYYMMDD-HHMMSS`
      * ``"snapshot"``      — daily snapshot in `<DATA_DIR>/snapshots/`
      * ``"lost_entries"``  — shrink-guard spillover in
                               `<DATA_DIR>/lost_entries/`
    """
    found: list[dict] = []
    data_dir = target_path.parent
    # Legacy single-gen backup.
    legacy = target_path.with_suffix(target_path.suffix + ".bak")
    if legacy.exists():
        info = _backup_info(legacy)
        if info:
            found.append({"kind": "legacy_bak",
                          "source_path": legacy, **info})
    # Rotating multi-gen backups (covers both `.bak.<ts>` and
    # collision-bumped `.bak.<ts>.<N>`).
    for bak in _iter_backups(target_path):
        info = _backup_info(bak)
        if info:
            found.append({"kind": "rotating_bak",
                          "source_path": bak, **info})
    # Daily snapshots.
    snap_dir = data_dir / _state._SNAPSHOT_DIR_NAME
    if snap_dir.exists():
        try:
            stem = target_path.stem
            ext = target_path.suffix
            for snap in snap_dir.glob(f"{stem}-*{ext}"):
                info = _backup_info(snap)
                if info:
                    found.append({"kind": "snapshot",
                                  "source_path": snap, **info})
        except OSError:
            pass
    # Lost-entries spillover.
    lost_dir = data_dir / _state._LOST_ENTRIES_DIR_NAME
    if lost_dir.exists():
        try:
            for lost in lost_dir.glob(f"{target_path.stem}-*.json"):
                info = _backup_info(lost)
                if info:
                    found.append({"kind": "lost_entries",
                                  "source_path": lost, **info})
        except OSError:
            pass
    # Newest first by mtime string (lex-sortable).
    found.sort(key=lambda x: x["mtime_str"], reverse=True)
    return found


def _restore_from_backup(target_path: Path, source_path: Path,
                          label: str) -> int:
    """Read entries from `source_path` and write them onto
    `target_path` via `_safe_save_json` so the current state is
    automatically backed up under the new rotating-backup discipline
    before being overwritten. Returns the number of entries
    restored. Raises ValueError on unparseable / oversized source,
    OSError on write failure.

    Size-capped at `_state._SAFE_LOAD_JSON_MAX_BYTES` to match
    `_safe_load_json` — a co-resident attacker who planted a 50 GB
    file under `<DATA_DIR>/` would otherwise OOM the restore worker.
    Symlink-rejected via lstat.
    """
    ok, reason = _safe_file_size_check(
        source_path, _state._SAFE_LOAD_JSON_MAX_BYTES, "backup",
    )
    if not ok:
        raise ValueError(reason or "backup file rejected")
    try:
        raw = json.loads(source_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"unreadable backup: {exc}") from exc
    # If the source backup carries a higher-than-current schema version,
    # stamp it under the *target* path so the subsequent `_safe_save_json`
    # write preserves the v stamp instead of demoting it to
    # `_state._CURRENT_SCHEMA_VERSION`. Mirrors `_safe_load_json`'s observation
    # logic — restoring a v2 backup onto a v1 SpliceCraft no longer
    # silently downgrades the user's data on the next save.
    if isinstance(raw, dict):
        try:
            observed_v = int(raw.get("_schema_version", 0))
        except (TypeError, ValueError):
            observed_v = 0
        if observed_v > _state._CURRENT_SCHEMA_VERSION:
            _state._OBSERVED_SCHEMA_VERSIONS[str(target_path)] = observed_v
    entries, _shape_warn = _extract_entries(raw, label)
    if entries is None:
        raise ValueError(
            f"backup {source_path.name} is not a recognisable "
            f"{label} payload"
        )
    # Restore can legitimately write a small backup over a large live
    # file (e.g. user picks `library.json.bak.<old>` after a recent
    # bulk add). Opt in to the L3 catastrophic-shrink bypass so the
    # restore isn't blocked by the same guard that catches accidental
    # wipes.
    # Hold `_state._cache_lock` across the write so this restore serialises
    # against the domain `_save_X` mirror saves (which take the same lock around
    # `_safe_save_json` + the sibling mirror) — otherwise a concurrent agent-
    # thread save could interleave and drift the library<->collections mirror
    # (INV-50). `_safe_save_json` takes no lock itself, so there is no re-entry.
    with _state._cache_lock, _allow_catastrophic_shrink():
        _safe_save_json(target_path, entries, label)
    return len(entries)


# Mapping of user-friendly data-file labels to the splicecraft attr
# that holds the live path. Agents pass the friendly label so we don't
# leak filesystem-paths into the wire surface, and so a relocated
# DATA_DIR doesn't require an agent-side change.
_AGENT_BACKUP_LABELS: dict = {
    "plasmid_library":         "_LIBRARY_FILE",
    "collections":             "_COLLECTIONS_FILE",
    "parts_bin":               "_PARTS_BIN_FILE",
    "parts_bin_collections":   "_PARTS_BIN_COLLECTIONS_FILE",
    "primers":                 "_PRIMERS_FILE",
    "primer_collections":      "_PRIMER_COLLECTIONS_FILE",
    "features":                "_FEATURES_FILE",
    "feature_colors":          "_FEATURE_COLORS_FILE",
    "grammars":                "_GRAMMARS_FILE",
    "entry_vectors":           "_ENTRY_VECTORS_FILE",
    "codon_tables":            "_CODON_TABLES_FILE",
    "settings":                "_SETTINGS_FILE",
    # Parity with `RestoreFromBackupModal._TARGETS` — agent should
    # be able to list/restore every persisted user-data file the GUI
    # can. Pre-fix the agent rejected Experiments/Gels/Protein-motifs
    # labels even though their `.bak` files exist on disk.
    "experiments":             "_EXPERIMENTS_FILE",
    "experiment_projects":     "_EXPERIMENT_PROJECTS_FILE",
    "gels":                    "_GELS_FILE",
    "protein_motifs":          "_PROTEIN_MOTIFS_FILE",
    # Enzyme catalog extensions (2026-05-22). Parity with
    # `RestoreFromBackupModal._TARGETS` and `_USER_DATA_FILE_ATTRS`.
    "custom_enzymes":          "_CUSTOM_ENZYMES_FILE",
    "enzyme_collections":      "_ENZYME_COLLECTIONS_FILE",
    # Sweep #28: HMM database registry (covers the catalog JSON, not
    # the per-DB downloads under hmm_databases/<id>/ — those are large
    # binaries handled by the file-level Master Delete sweep, not by
    # the .bak rotation backup chain).
    "hmm_db_catalog":          "_HMM_DB_CATALOG_FILE",
    "protocol_collections":    "_PROTOCOL_COLLECTIONS_FILE",
    "custom_labware":          "_CUSTOM_LABWARE_FILE",
    "protein_collections":     "_PROTEIN_COLLECTIONS_FILE",
    # BABS model-picker collections (INV-139, 2026-06-29).
    "model_collections":       "_MODEL_COLLECTIONS_FILE",
}


def _resolve_backup_label(label: str) -> "tuple[Path | None, str]":
    """Map a label → Path + human label. Returns (None, error_msg)
    if the label isn't registered. Used by both list/restore agent
    endpoints so the validation is uniform."""
    if not isinstance(label, str) or not label:
        return (None, "missing 'label'")
    attr = _AGENT_BACKUP_LABELS.get(label)
    if attr is None:
        return (None,
                f"unknown label {label!r} (valid: "
                f"{sorted(_AGENT_BACKUP_LABELS)})")
    path = _state._resolve_data_attr_hook(attr)
    if not isinstance(path, Path):
        return (None, f"label {label!r} not mapped to a Path")
    return (path, label)
