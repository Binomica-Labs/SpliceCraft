"""test_blob_store — content-addressed plasmid blob store (v1.0.23).

The blob store holds each unique ``gb_text`` ONCE, content-addressed by
sha256, immutable + write-once. These tests cover the primitives;
dehydrate/rehydrate + migration + orphan-GC live in test_data_safety.py.

The autouse ``_protect_user_data`` fixture (conftest) sandboxes
``_DATA_DIR`` to ``tmp_path``, so the store writes under
``tmp_path/plasmid_blobs`` — never the user's real data dir.
"""
from __future__ import annotations

import hashlib
import json
import os

import pytest

import splicecraft as sc


def _age(path, seconds: float):
    """Backdate a file's mtime by `seconds` (to step past GC grace /
    quarantine-retention windows deterministically)."""
    st = path.stat()
    old = st.st_mtime - seconds
    os.utime(path, (old, old))

GB = "LOCUS       x         10 bp\nORIGIN\n1 acgtacgtac\n//\n"


class TestBlobHash:
    def test_hash_matches_sha256(self):
        assert sc._blob_hash(GB) == hashlib.sha256(
            GB.encode("utf-8")).hexdigest()

    def test_hash_deterministic(self):
        assert sc._blob_hash(GB) == sc._blob_hash(GB)

    def test_hash_is_64_lowercase_hex(self):
        h = sc._blob_hash(GB)
        assert len(h) == 64
        assert all(c in "0123456789abcdef" for c in h)


class TestBlobPath:
    def test_valid_ref_path(self):
        ref = sc._blob_hash(GB)
        p = sc._blob_path(ref)
        assert p.name == f"{ref}.gb"
        assert p.parent == sc._plasmid_blob_dir()

    @pytest.mark.parametrize("bad", [
        "", "xyz", "../../etc/passwd", "A" * 64,      # uppercase rejected
        "g" * 64, "0" * 63, "0" * 65, "/abs", "ab/cd",
        "../" + "0" * 61,
    ])
    def test_invalid_ref_rejected(self, bad):
        with pytest.raises(ValueError):
            sc._blob_path(bad)

    def test_non_str_ref_rejected(self):
        with pytest.raises(ValueError):
            sc._blob_path(None)  # type: ignore[arg-type]


class TestBlobWriteRead:
    def test_write_returns_correct_ref(self):
        assert sc._blob_write(GB) == sc._blob_hash(GB)

    def test_round_trip(self):
        ref = sc._blob_write(GB)
        assert sc._blob_read(ref) == GB

    def test_write_creates_file(self):
        ref = sc._blob_write(GB)
        assert sc._blob_path(ref).is_file()

    def test_write_idempotent_no_rewrite(self):
        ref1 = sc._blob_write(GB)
        mtime1 = sc._blob_path(ref1).stat().st_mtime_ns
        ref2 = sc._blob_write(GB)                  # identical content
        assert ref1 == ref2
        # Immutable: the existing blob is left untouched (write skipped).
        assert sc._blob_path(ref2).stat().st_mtime_ns == mtime1

    def test_distinct_content_distinct_blobs(self):
        r1 = sc._blob_write(GB)
        r2 = sc._blob_write(GB + "X")
        assert r1 != r2
        assert sc._blob_read(r1) == GB
        assert sc._blob_read(r2) == GB + "X"

    def test_empty_string_blob(self):
        ref = sc._blob_write("")
        assert ref == hashlib.sha256(b"").hexdigest()
        assert sc._blob_read(ref) == ""

    def test_unicode_round_trip(self):
        s = GB + '  /note="éü中文"\n'
        ref = sc._blob_write(s)
        assert sc._blob_read(ref) == s

    def test_write_non_str_raises(self):
        with pytest.raises(ValueError):
            sc._blob_write(b"bytes")  # type: ignore[arg-type]

    def test_exists(self):
        ref = sc._blob_hash(GB)
        assert sc._blob_exists(ref) is False
        sc._blob_write(GB)
        assert sc._blob_exists(ref) is True

    def test_exists_invalid_ref_is_false(self):
        assert sc._blob_exists("nope") is False


class TestBlobReadFailures:
    """Reads must NEVER return a wrong/empty sequence in place of a real
    one — every failure mode returns None (loud log) so callers can tell
    'unavailable' from 'empty'."""

    def test_missing_returns_none(self):
        ref = sc._blob_hash("never written " + GB)   # valid shape, absent
        assert sc._blob_read(ref) is None

    def test_invalid_ref_returns_none(self):
        assert sc._blob_read("../escape") is None
        assert sc._blob_read("") is None

    def test_corrupt_blob_detected(self):
        ref = sc._blob_write(GB)
        sc._blob_path(ref).write_bytes(b"TAMPERED - different content")
        # Hash no longer matches the ref → refuse rather than return junk.
        assert sc._blob_read(ref) is None

    def test_truncated_blob_detected(self):
        ref = sc._blob_write(GB)
        sc._blob_path(ref).write_bytes(GB.encode("utf-8")[:-5])
        assert sc._blob_read(ref) is None

    def test_post_write_verification_catches_bad_write(self, monkeypatch):
        """If the atomic writer lands bytes that don't hash to the ref,
        _blob_write must detect it, unlink the bad file, and raise — a
        corrupt write never becomes a trusted blob."""
        def _bad_write(path, data):
            path.write_bytes(b"corrupted-does-not-match-hash")
        monkeypatch.setattr(sc, "_atomic_write_bytes", _bad_write)
        payload = "fresh content " + GB
        with pytest.raises(OSError):
            sc._blob_write(payload)
        # The bad file was cleaned up.
        ref = sc._blob_hash(payload)
        assert not sc._blob_path(ref).exists()


def _entry(eid="p1", **extra):
    e = {"id": eid, "name": "Plasmid 1", "size": 10, "gb_text": GB}
    e.update(extra)
    return e


class TestDehydrateEntry:
    def test_strips_gb_text_to_ref(self):
        out = sc._dehydrate_entry(_entry())
        assert "gb_text" not in out
        assert out["gb_ref"] == sc._blob_hash(GB)
        assert sc._blob_path(out["gb_ref"]).is_file()

    def test_preserves_other_fields(self):
        out = sc._dehydrate_entry(_entry(status="verified", n_feats=3))
        assert out["id"] == "p1" and out["name"] == "Plasmid 1"
        assert out["size"] == 10 and out["status"] == "verified"
        assert out["n_feats"] == 3

    def test_does_not_mutate_input(self):
        e = _entry()
        sc._dehydrate_entry(e)
        assert e["gb_text"] == GB and "gb_ref" not in e

    def test_idempotent(self):
        once = sc._dehydrate_entry(_entry())
        twice = sc._dehydrate_entry(once)
        assert twice == once                       # already a ref → unchanged

    def test_empty_gb_text_unchanged(self):
        e = _entry(gb_text="")
        assert sc._dehydrate_entry(e) == e

    def test_missing_gb_text_unchanged(self):
        e = {"id": "p1", "name": "x"}
        assert sc._dehydrate_entry(e) == e

    def test_non_dict_passthrough(self):
        assert sc._dehydrate_entry("nope") == "nope"  # type: ignore[arg-type]

    def test_aborts_if_blob_write_fails(self, monkeypatch):
        def _boom(_):
            raise OSError("disk full")
        monkeypatch.setattr(sc, "_blob_write", _boom)
        with pytest.raises(OSError):
            sc._dehydrate_entry(_entry())


class TestRehydrateEntry:
    def test_resolves_ref_to_gb_text(self):
        dh = sc._dehydrate_entry(_entry())
        rh = sc._rehydrate_entry(dh)
        assert rh["gb_text"] == GB
        assert "gb_ref" not in rh

    def test_round_trip_preserves_everything(self):
        e = _entry(status="cloned", source="import")
        rt = sc._rehydrate_entry(sc._dehydrate_entry(e))
        assert rt == e

    def test_inline_gb_text_wins_and_drops_stray_ref(self):
        e = {"id": "p1", "gb_text": GB, "gb_ref": "0" * 64}
        out = sc._rehydrate_entry(e)
        assert out["gb_text"] == GB
        assert "gb_ref" not in out

    def test_does_not_mutate_input(self):
        dh = sc._dehydrate_entry(_entry())
        ref = dh["gb_ref"]
        sc._rehydrate_entry(dh)
        assert dh["gb_ref"] == ref and "gb_text" not in dh

    def test_missing_blob_keeps_ref_and_empties_text(self):
        dh = sc._dehydrate_entry(_entry())
        ref = dh["gb_ref"]
        sc._blob_path(ref).unlink()                # blob vanishes
        out = sc._rehydrate_entry(dh)
        assert out["gb_text"] == ""
        assert out["gb_ref"] == ref                # retained for recovery

    def test_auto_recovers_when_blob_restored(self):
        dh = sc._dehydrate_entry(_entry())
        ref = dh["gb_ref"]
        blob = sc._blob_path(ref)
        saved = blob.read_bytes()
        blob.unlink()
        unresolved = sc._rehydrate_entry(dh)
        assert unresolved["gb_text"] == "" and unresolved["gb_ref"] == ref
        blob.write_bytes(saved)                    # restore from "backup"
        recovered = sc._rehydrate_entry(unresolved)
        assert recovered["gb_text"] == GB and "gb_ref" not in recovered

    def test_no_text_no_ref_unchanged(self):
        e = {"id": "p1", "name": "x"}
        assert sc._rehydrate_entry(e) == e


class TestDehydrateDedupAndLists:
    def test_same_content_dedups_to_one_blob(self):
        a = sc._dehydrate_entry(_entry("a"))
        b = sc._dehydrate_entry(_entry("b"))          # same gb_text
        assert a["gb_ref"] == b["gb_ref"]
        blobs = list(sc._plasmid_blob_dir().glob("*.gb"))
        assert len(blobs) == 1                         # stored once

    def test_entries_round_trip(self):
        entries = [_entry("a"), _entry("b", gb_text=GB + "B"),
                   {"id": "c", "name": "no-seq"}]
        rt = sc._rehydrate_entries(sc._dehydrate_entries(entries))
        assert rt == entries

    def test_dehydrate_entries_aborts_on_failure(self, monkeypatch):
        calls = {"n": 0}
        real = sc._blob_write

        def _flaky(text):
            calls["n"] += 1
            if calls["n"] == 2:
                raise OSError("disk full mid-list")
            return real(text)
        monkeypatch.setattr(sc, "_blob_write", _flaky)
        with pytest.raises(OSError):
            sc._dehydrate_entries([_entry("a"), _entry("b", gb_text=GB + "B")])

    def test_collections_round_trip(self):
        colls = [
            {"name": "C1", "plasmids": [_entry("a"), _entry("b", gb_text=GB + "B")]},
            {"name": "C2", "plasmids": [_entry("a")]},   # shares blob with C1
            {"name": "empty", "plasmids": []},
        ]
        dh = sc._dehydrate_collections(colls)
        # Cross-collection dedup: "a" appears in C1 + C2 → one blob.
        refs = {p["gb_ref"] for c in dh for p in c["plasmids"]}
        assert len(list(sc._plasmid_blob_dir().glob("*.gb"))) == len(refs)
        rt = sc._rehydrate_collections(dh)
        assert rt == colls

    def test_inline_detectors(self):
        assert sc._entries_have_inline_gb_text([_entry()]) is True
        assert sc._entries_have_inline_gb_text(
            sc._dehydrate_entries([_entry()])) is False
        assert sc._collections_have_inline_gb_text(
            [{"name": "C", "plasmids": [_entry()]}]) is True
        assert sc._collections_have_inline_gb_text(
            sc._dehydrate_collections(
                [{"name": "C", "plasmids": [_entry()]}])) is False


class TestLibraryBoundaryIntegration:
    """End-to-end through the LIVE save/load path: dehydrate-on-save +
    rehydrate-on-load via _save_library / _load_library. conftest
    sandboxes _DATA_DIR + authorizes writes + forces sync mirror."""

    def test_save_dehydrates_disk_load_rehydrates_cache(self):
        gb = "LOCUS       p1         5 bp\nORIGIN\n1 acgta\n//\n"
        sc._save_library([{"id": "P1", "name": "p one", "gb_text": gb}])
        # On disk: dehydrated (gb_ref, no inline gb_text).
        raw = json.loads(sc._LIBRARY_FILE.read_text())
        disk = raw["entries"][0]
        assert "gb_text" not in disk
        assert disk["gb_ref"] == sc._blob_hash(gb)
        assert sc._blob_path(disk["gb_ref"]).is_file()
        # Bust cache + reload: fully materialised again (gb_text back).
        sc._library_cache = None
        loaded = sc._load_library()
        assert loaded[0]["gb_text"] == gb
        assert "gb_ref" not in loaded[0]
        assert loaded[0]["id"] == "P1" and loaded[0]["name"] == "p one"

    def test_backward_compat_loads_old_inline_format(self):
        gb = "LOCUS       p2         5 bp\nORIGIN\n1 acgta\n//\n"
        # Old-format file: inline gb_text, no gb_ref, no blob.
        sc._LIBRARY_FILE.write_text(json.dumps(
            {"_schema_version": 1, "entries": [{"id": "P2", "gb_text": gb}]}))
        sc._library_cache = None
        loaded = sc._load_library()
        assert loaded[0]["gb_text"] == gb            # inline still resolves
        # A subsequent save migrates it to the blob format on disk.
        sc._save_library(loaded)
        raw = json.loads(sc._LIBRARY_FILE.read_text())
        assert "gb_text" not in raw["entries"][0]
        assert raw["entries"][0]["gb_ref"] == sc._blob_hash(gb)
        assert sc._blob_path(sc._blob_hash(gb)).is_file()

    def test_missing_blob_keeps_entry_does_not_drop(self):
        gb = "LOCUS       p3         5 bp\nORIGIN\n1 acgta\n//\n"
        sc._save_library([{"id": "P3", "name": "p3", "gb_text": gb}])
        ref = sc._blob_hash(gb)
        sc._blob_path(ref).unlink()                  # blob lost externally
        sc._library_cache = None
        loaded = sc._load_library()
        assert len(loaded) == 1                      # entry NOT dropped
        assert loaded[0]["id"] == "P3"
        assert loaded[0]["gb_text"] == ""            # unavailable, not fake
        assert loaded[0]["gb_ref"] == ref            # retained for recovery

    def test_dedup_same_sequence_one_blob(self):
        gb = "LOCUS       d          5 bp\nORIGIN\n1 acgta\n//\n"
        sc._save_library([
            {"id": "A", "name": "a", "gb_text": gb},
            {"id": "B", "name": "b", "gb_text": gb},   # identical sequence
        ])
        # Two entries, identical gb_text → exactly one blob on disk.
        assert len(list(sc._plasmid_blob_dir().glob("*.gb"))) == 1
        raw = json.loads(sc._LIBRARY_FILE.read_text())
        refs = {e["gb_ref"] for e in raw["entries"]}
        assert refs == {sc._blob_hash(gb)}

    def test_collections_save_dehydrates(self):
        gb = "LOCUS       c          5 bp\nORIGIN\n1 acgta\n//\n"
        sc._save_collections([
            {"name": "C1", "plasmids": [{"id": "P1", "gb_text": gb}]},
        ])
        raw = json.loads(sc._COLLECTIONS_FILE.read_text())
        plas = raw["entries"][0]["plasmids"][0]
        assert "gb_text" not in plas
        assert plas["gb_ref"] == sc._blob_hash(gb)
        sc._collections_cache = None
        loaded = sc._load_collections()
        assert loaded[0]["plasmids"][0]["gb_text"] == gb


def _write_lib_refs(refs):
    sc._LIBRARY_FILE.write_text(json.dumps({
        "_schema_version": 1,
        "entries": [{"id": f"e{i}", "gb_ref": r} for i, r in enumerate(refs)],
    }))


class TestOrphanBlobGC:
    """GC quarantines unreferenced blobs — NEVER deletes, NEVER touches a
    referenced (or recently-written, or backup-referenced) blob, and aborts
    entirely if it can't fully trust the ref-set."""

    def test_quarantines_orphan_not_deletes(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        assert sc._gc_orphan_blobs() == 1
        assert not sc._blob_path(ref).exists()            # moved out of store
        q = (sc._plasmid_blob_dir() / sc._BLOB_ORPHAN_DIR_NAME / f"{ref}.gb")
        assert q.is_file()                                 # quarantined
        assert q.read_bytes() == GB.encode("utf-8")        # content intact

    def test_keeps_blob_referenced_by_library(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        _write_lib_refs([ref])
        assert sc._gc_orphan_blobs() == 0
        assert sc._blob_path(ref).is_file()

    def test_keeps_blob_referenced_by_collection(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        sc._COLLECTIONS_FILE.write_text(json.dumps({
            "_schema_version": 1,
            "entries": [{"name": "C",
                         "plasmids": [{"id": "p", "gb_ref": ref}]}],
        }))
        assert sc._gc_orphan_blobs() == 0
        assert sc._blob_path(ref).is_file()

    def test_keeps_blob_referenced_only_by_backup(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        bak = sc._LIBRARY_FILE.with_name(
            sc._LIBRARY_FILE.name + ".bak.20260101-000000")
        bak.write_text(json.dumps(
            {"_schema_version": 1, "entries": [{"id": "e", "gb_ref": ref}]}))
        assert sc._gc_orphan_blobs() == 0                  # rollback-protected
        assert sc._blob_path(ref).is_file()

    def test_keeps_blob_referenced_only_by_daily_snapshot(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        snap = sc._DATA_DIR / sc._SNAPSHOT_DIR_NAME
        snap.mkdir(parents=True, exist_ok=True)
        (snap / "plasmid_library-2026-06-03.json").write_text(json.dumps(
            {"_schema_version": 1, "entries": [{"id": "e", "gb_ref": ref}]}))
        assert sc._gc_orphan_blobs() == 0                  # snapshot-restore safe
        assert sc._blob_path(ref).is_file()

    def test_keeps_blob_referenced_only_by_lost_entries_spill(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        ld = sc._DATA_DIR / sc._LOST_ENTRIES_DIR_NAME
        ld.mkdir(parents=True, exist_ok=True)
        (ld / "collections-20260603-000000.json").write_text(json.dumps(
            {"_schema_version": 1,
             "entries": [{"name": "C",
                          "plasmids": [{"id": "p", "gb_ref": ref}]}]}))
        assert sc._gc_orphan_blobs() == 0                  # spill-recovery safe
        assert sc._blob_path(ref).is_file()

    def test_respects_grace_period(self):
        ref = sc._blob_write(GB)                            # fresh mtime
        assert sc._gc_orphan_blobs() == 0                  # too recent
        assert sc._blob_path(ref).is_file()

    def test_aborts_when_current_metadata_corrupt(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        sc._LIBRARY_FILE.write_text("{ corrupt json not closeable")
        # Eligible orphan, but an unreadable current file → ABORT, touch nothing.
        assert sc._gc_orphan_blobs() == 0
        assert sc._blob_path(ref).is_file()

    def test_no_blob_dir_is_noop(self):
        assert sc._gc_orphan_blobs() == 0                  # no dir → 0, no raise

    def test_idempotent(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        assert sc._gc_orphan_blobs() == 1
        assert sc._gc_orphan_blobs() == 0                  # nothing left

    def test_quarantine_retention_prune(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        sc._gc_orphan_blobs()
        q = (sc._plasmid_blob_dir() / sc._BLOB_ORPHAN_DIR_NAME / f"{ref}.gb")
        assert q.is_file()
        _age(q, sc._BLOB_ORPHAN_RETENTION_DAYS * 86400 + 100)
        sc._prune_blob_quarantine(q.parent)
        assert not q.exists()                              # pruned past retention

    def test_quarantine_keeps_recent(self):
        ref = sc._blob_write(GB)
        _age(sc._blob_path(ref), sc._BLOB_GC_GRACE_SECONDS + 100)
        sc._gc_orphan_blobs()
        q = (sc._plasmid_blob_dir() / sc._BLOB_ORPHAN_DIR_NAME / f"{ref}.gb")
        sc._prune_blob_quarantine(q.parent)               # fresh → kept
        assert q.is_file()
