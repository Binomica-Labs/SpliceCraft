#!/usr/bin/env python3
"""Bulk-import CommercialSaaS .dna (or GenBank .gb/.gbk) files into a SpliceCraft
collection.

Usage:
    python3 scripts/bulk_import.py --dir /path/to/folder
    python3 scripts/bulk_import.py --dir tests --collection "FFE Trial"
    SPLICECRAFT_DATA_DIR=/tmp/sc-trial python3 scripts/bulk_import.py --dir tests

Each file is imported independently — a single corrupt file doesn't abort the
batch. The script writes only to collections.json (which auto-backs-up via
_safe_save_json); the active collection and plasmid_library.json are NOT
touched, so a botched import is reversible by deleting the new collection
from the LibraryPanel.

By design: refuses to overwrite an existing collection of the same name.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make `splicecraft` importable when running from anywhere.
_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

import splicecraft as sc


def main() -> int:
    p = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    p.add_argument("--dir", "-d", required=True, type=Path,
                   help="Folder to scan for .dna / .gb / .gbk files")
    p.add_argument("--collection", "-c", default=None,
                   help="Target collection name (default: folder basename)")
    args = p.parse_args()

    folder: Path = args.dir.expanduser().resolve()
    if not folder.is_dir():
        print(f"error: {folder} is not a directory", file=sys.stderr)
        return 2

    coll_name = args.collection or folder.name
    if sc._collection_name_taken(coll_name):
        print(f"error: collection '{coll_name}' already exists. "
              f"Pick a different --collection name or delete the existing one.",
              file=sys.stderr)
        return 2

    print(f"Scanning {folder}")
    print(f"Target data dir: {sc._state._DATA_DIR}")
    print(f"Target collection: '{coll_name}'\n")

    entries, failures = sc._bulk_import_folder(folder)

    # Folder-level notices (a skipped subfolder, the entry cap) are not
    # failed files: they must not print as FAIL, inflate the file total,
    # or flip the exit code. See issue #31.
    bad_files, folder_notes = sc._bulk_import_split_failures(folder, failures)

    for entry in entries:
        print(f"  ok   {entry['source'][5:]}  "
              f"({entry['size']:,} bp, {entry['n_feats']} feats)")
    for path, reason in bad_files:
        print(f"  FAIL {path.name}: {reason}")
    for _, reason in folder_notes:
        print(f"  NOTE {reason}")

    if not entries:
        print("\nNothing imported.", file=sys.stderr)
        return 1

    colls = sc._load_collections()
    colls.append({
        "name":        coll_name,
        "description": f"Bulk imported from {folder}",
        "plasmids":    entries,
        "saved":       sc._date.today().isoformat(),
    })
    sc._save_collections(colls)

    total = len(entries) + len(bad_files)
    print(f"\n{'='*60}")
    print(f"Imported  {len(entries)} / {total} file(s)")
    if bad_files:
        print(f"Failed    {len(bad_files)} (see log for tracebacks)")
    for _, reason in folder_notes:
        print(f"Note      {reason}")
    print(f"Written to collection '{coll_name}' in {sc._COLLECTIONS_FILE}")
    return 0 if not bad_files else 1


if __name__ == "__main__":
    sys.exit(main())
