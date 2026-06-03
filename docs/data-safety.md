# Data safety and backups

A workbench you trust your day to has to behave like one. SpliceCraft
takes data safety and predictable behaviour as first-class design
constraints.

## Four-layer data-safety net for every JSON write

For every save through `_safe_save_json` (the only sanctioned path):

1. **Atomic write** via `tempfile.mkstemp` + `os.fsync` + `os.replace`
   (plus a parent-dir `fsync` so the rename is journalled), with the
   prior version copied to `*.json.bak` first. Before the swap the
   freshly-written temp file is **read back and validated** — a full
   `json.loads` + entry-count check for files ≤ 32 MB, a cheap
   "non-empty and closes with `}`" tail check for larger ones — so a
   silently-truncated write is caught and the live file is left intact
   instead of being overwritten with garbage.
2. **Rotating timestamped backups** (`*.json.bak.YYYYMMDD-HHMMSS`,
   last 10 retained per file) so an *old* good copy is recoverable,
   not just the most recent. Two refinements keep this from eating
   disk on big libraries:
   - **Compression** — at launch, every timestamped backup except the
     newest is gzipped (`*.bak.<ts>.gz`); plasmid JSON compresses
     ~4-5×. The newest timestamped backup and the legacy `.bak` stay
     **plain** so the fast recovery path never has to decompress.
   - **Aggregate byte cap** — after the count cap, the oldest backups
     are dropped until the per-file total falls under 1 GB (never
     below 2 generations). A 274 MB `collections.json` used to accrete
     10 × 274 MB ≈ 2.7 GB of backups per file.
   - **Dedup** — a save whose prior file is byte-identical to the most
     recent backup skips both backup writes (no more "N identical
     274 MB backups in one wall-second").
3. **Daily per-file snapshots** to `<data dir>/snapshots/` (30 days
   retained) — written once per calendar day at launch.
4. **Suspicious-shrink guard**: if a save would discard >50% of
   entries (with ≥5 prior), the discarded entries are spilled to
   `<data dir>/lost_entries/` **before** the overwrite proceeds; a
   >90% discard (≥10 prior) is **refused** outright unless an explicit
   bypass is armed. **Mirror swaps are exempt** — switching the active
   collection / parts-bin / primer-set / project rewrites the active
   mirror file (e.g. `plasmid_library.json`) from one collection's
   contents to another's, which legitimately shrinks it; the dropped
   entries are *not* lost (they live in the sibling `collections.json`
   / `*_collections.json`), so those writes go through
   `_switch_active_collection_library` / `_safe_save_json_mirror`,
   which neither spill a redundant copy nor refuse a big shrink. (This
   also fixes a latent bug where switching from a large to a tiny
   collection could be refused outright.) The `lost_entries/` directory
   is now itself bounded (last 5 files, ≤ 500 MB) — pre-1.0.22 it grew
   without limit, reaching 1.5 GB of switch-driven spills.

**Recovery.** On load, a corrupt main file is restored from the legacy
`.bak`; if **both** the main file and `.bak` are corrupt (e.g. two bad
saves in a row), `_safe_load_json` walks the timestamped rotation
newest→oldest — transparently decompressing `.gz` backups — and
restores from the first that parses. The corrupt main is renamed aside
to `*.corrupt-<ts>` for forensics first.

**Launch housekeeping.** A background thread at startup
(`_run_data_dir_housekeeping`) compresses older backups, enforces the
byte caps, and prunes `lost_entries/` — so the disk savings above also
reclaim *existing* residue on the next launch, not just future writes.
It runs off the UI thread because compressing a 274 MB backup is
CPU-heavy.

**Settings → Restore from backup…** surfaces every recoverable copy
across all four tiers; pick a row, get a one-click restore (the
pre-restore state goes through the same backup chain, so even an
accidental restore is reversible).

## Pre-update snapshots

Every `splicecraft update` snapshots your full library, collections,
parts bin, primers, feature library, grammars, codon tables, settings,
crash-recovery autosaves, and `.dna` sidecars **before** invoking
pip / pipx / uv / pixi. If the snapshot can't be taken (disk full,
permissions), the upgrade aborts.

Snapshots live in a **sibling directory** of the data dir
(`<DATA_DIR>/../<DATA_DIR.name>-update-backups/`, override
`$SPLICECRAFT_UPDATE_BACKUP_DIR`) so a hypothetical recursive-wipe bug
in a new version cannot touch them.

```bash
splicecraft update --list-snapshots
splicecraft update --restore-pre-update latest
splicecraft update --restore-pre-update 20260514-143022-abc123__from-0.8.5
```

The pre-update snapshot is itself reversible (a pre-restore snapshot
is taken before any restore), so even an accidental rollback can be
undone.

### Sacred restore checks

`_restore_pre_update_snapshot` enforces four checks on every restore
candidate before any `os.replace` runs:

1. **schema_version** ≤ `_PRE_UPDATE_SCHEMA_VERSION`
2. **attr** in the user-data whitelist
3. **name** rejects path separators / `..`
4. **SHA-256** re-verified against the manifest before `os.replace`

A manifest with `sha256` missing / empty is refused outright.

## Crash-recovery autosave

Dirty edits debounce a 3-second write to a per-record `.gb` snapshot
in `<DATA_DIR>/crash_recovery/`. Power-cut your laptop mid-edit; the
next launch surfaces the survivors via a toast.

## Lock + concurrency hardening

- **Per-data-dir lockfile** at `<DATA_DIR>/splicecraft.lock`
  (POSIX `fcntl.flock`, Win `msvcrt.locking`). PID is `fsync`-ed
  before acquire returns. `SPLICECRAFT_SKIP_LOCK=1` bypass for CI.
- **Stale-PID detection** (`os.kill(pid, 0)`) lets a SpliceCraft
  killed on a shared filesystem release its lock on the next launch.
- **`_cache_lock` (RLock)** wraps every `_save_*` + cache
  reassignment so concurrent saves can't land
  `os.replace` A→B while cache reassignments land B→A. RLock
  because chains nest (`_save_library` ⇒ `_sync_active_collection_plasmids`
  ⇒ `_save_collections`).
- **Modal cap** dispatches `callback(None)` on overflow so
  modal-push fanouts can't pile up indefinitely.

## Defence-in-depth size caps

Every external input has a documented ceiling. Sample:

| Surface                | Cap                          |
|------------------------|------------------------------|
| `_safe_load_json`      | 1 GB                         |
| `_h_load_file` (agent) | 50 MB (`force=true` override) |
| `_gb_text_to_record`   | 64 MB                        |
| Plasmidsaurus zip      | 500 MB, 50 MB / member, 2000 members |
| `.dna` history XML     | streaming LZMA with `max_length` |
| NCBI / PyPI / Kazusa   | per-fetch `_*_MAX_RESPONSE_BYTES` |
| Pre-update manifest    | 4 MB                         |
| CLI sidecar response   | 50 MB                        |
| CLI token file         | 1 KB                         |

See [SECURITY.md](
https://github.com/Binomica-Labs/SpliceCraft/blob/master/SECURITY.md)
for the full threat-model writeup.

## What lives where

All user data persists as human-readable JSON in the user data
directory.

| File                            | Purpose                                                  |
|---------------------------------|----------------------------------------------------------|
| `collections.json`              | Named collections of plasmids — source of truth          |
| `plasmid_library.json`          | Live mirror of the active collection's plasmids          |
| `parts_bin.json`                | Active parts-bin's user-domesticated cloning parts       |
| `parts_bin_collections.json`    | Named parts-bin snapshots (multi-bin storage)            |
| `primers.json`                  | Designed primer library                                  |
| `features.json`                 | Reusable feature snippets                                |
| `feature_colors.json`           | Per-type feature color overrides                         |
| `codon_tables.json`             | Cached codon-usage tables fetched from Kazusa            |
| `cloning_grammars.json`         | User-defined cloning grammars                            |
| `entry_vectors.json`            | Entry vectors bound to grammars                          |
| `settings.json`                 | App preferences                                          |
| `crash_recovery/*.gb`           | Per-record crash-recovery autosaves                      |
| `dna_originals/*.dna`           | Sidecars for round-tripping `.dna` files                 |
| `logs/splicecraft.log`          | Rotating per-session log (5 MB × 4)                      |
| `ui_snapshots/*.md`             | Alt+D bug-report dumps                                   |
| `snapshots/`, `*.bak.*`, `lost_entries/` | Four-layer JSON safety net                      |
| `../splicecraft-update-backups/` | Pre-update snapshots                                    |

The schema envelope (`{"_schema_version": 1, "entries": [...]}`)
silently accepts the legacy bare-list format (pre-0.3.1) and rewrites
it on the next save. Newer-version files load with a warning rather
than crashing.

## Diagnostic logging + UI snapshot + bundle

Three surfaces for bug-report archives:

- **Rotating log** at `<DATA_DIR>/logs/splicecraft.log` (override
  `$SPLICECRAFT_LOG`). 5 MB × 4 backups, 8-char session ID prefix.
  **Never logs sequence content** — `_repr_for_log` truncates /
  summarises any DNA-shaped payload.
- **`Alt+D` UI snapshot** → `<DATA_DIR>/ui_snapshots/ui-snapshot-<ts>.md`.
  Version, Python, platform, screen stack, focused widget, terminal
  size, settings, active collection / grammar, 200-line log tail
  with `/home/<user>` → `~`. Retention 20.
- **`splicecraft logs --bundle [--out PATH]`** atomically zips logs +
  last 5 UI snapshots + sanitized settings + system info + README
  into a single ZIP for emailing. Sequence content **never** leaks.
