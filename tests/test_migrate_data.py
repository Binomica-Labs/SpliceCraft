"""
test_migrate_data — portable export / import of ALL user data (File ▸
Migrate Data, 2026-06-09).

The feature wraps a pre-update-style snapshot in one portable ``.zip``:
export = snapshot → atomic zip; import = unzip → `_restore_pre_update_snapshot`
(automatic pre-import backup + sha256-verified atomic replace). These tests
pin the data-safety contract: byte-for-byte round-trips (incl. the
content-addressed blob store + embedded construction history), and the
refusal of foreign / future-format / corrupt / path-traversal archives with
the live data left untouched. The migrate manifest derives from the same
`_USER_DATA_FILE_ATTRS` registry that test_smoke's classification test
enforces — so a new data file is migrated automatically.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import zipfile

import pytest

import splicecraft as sc


def _gb(name: str, seq: str) -> str:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    r = SeqRecord(Seq(seq), id=name, name=name, description=name,
                  annotations={"molecule_type": "DNA", "topology": "circular"})
    return sc._record_to_gb_text(r)


def _seed() -> None:
    sc._save_collections([{"name": "Default", "plasmids": []}])
    sc._commit_library_entry_to_collection({
        "id": "p1", "name": "Plasmid 1", "size": 1200, "n_feats": 1,
        "source": "test", "added": "2026-06-09",
        "gb_text": _gb("P1", "ATGC" * 300),
        "history_xml": "<HistoryTree><Node>built</Node></HistoryTree>",
    }, "Default")
    sc._save_primers([{"name": "PR1", "sequence": "ACGTACGTACGT"}])
    sc._save_parts_bin([{"name": "part1", "type": "CDS", "sequence": "ATGAAATAA"}])


def _bust() -> None:
    for c in sc._MASTER_DELETE_CACHE_ATTRS:
        owner = sc._state if hasattr(sc._state, c) else sc
        if hasattr(owner, c):
            setattr(owner, c, None)


def _fingerprint() -> dict:
    fp = {}
    for attr in sc._USER_DATA_FILE_ATTRS:
        p = sc._resolve_state_or_hub(attr)
        if isinstance(p, sc.Path) and p.is_file():
            fp[attr] = hashlib.sha256(p.read_bytes()).hexdigest()
    bd = sc._plasmid_blob_dir()
    if bd.is_dir():
        for b in sorted(bd.rglob("*")):
            if b.is_file():
                fp["blob:" + b.name] = hashlib.sha256(b.read_bytes()).hexdigest()
    return fp


def _wipe() -> None:
    for attr in sc._USER_DATA_FILE_ATTRS:
        p = sc._resolve_state_or_hub(attr)
        if isinstance(p, sc.Path) and p.is_file():
            p.unlink()
    bd = sc._plasmid_blob_dir()
    if bd.is_dir():
        shutil.rmtree(bd)
    _bust()


class TestMigrateRoundTrip:
    def test_export_import_byte_for_byte(self, tmp_path):
        _seed()
        before = _fingerprint()
        assert len(before) >= 4
        zp = str(tmp_path / "data.zip")
        summ = sc._export_migrate_archive(zp)
        assert summ["bytes"] > 0
        assert zipfile.is_zipfile(zp)
        _wipe()
        assert not sc._state._COLLECTIONS_FILE.is_file()
        isumm = sc._import_migrate_archive(zp)
        assert isumm.get("pre_restore_snapshot"), "import must take a backup"
        assert not isumm.get("failed"), isumm.get("failed")
        _bust()
        assert _fingerprint() == before    # byte-for-byte, incl. blobs

    def test_history_and_sequences_survive(self, tmp_path):
        _seed()
        zp = str(tmp_path / "d.zip")
        sc._export_migrate_archive(zp)
        _wipe()
        sc._import_migrate_archive(zp)
        _bust()
        colls = sc._extract_entries(
            json.loads(sc._state._COLLECTIONS_FILE.read_text()), "c")[0] or []
        plasmids = next((c["plasmids"] for c in colls
                         if c.get("name") == "Default"), [])
        assert plasmids
        reh = [sc._rehydrate_entry(p) for p in plasmids]
        assert all("LOCUS" in (r.get("gb_text") or "") for r in reh)  # blobs
        # history_xml is now blob-dehydrated too (history_ref on disk), so check
        # the REHYDRATED entries — same as the gb_text assertion above.
        assert any("HistoryTree" in (r.get("history_xml") or "")
                   for r in reh)

    def test_archive_layout(self, tmp_path):
        _seed()
        zp = str(tmp_path / "d.zip")
        sc._export_migrate_archive(zp)
        with zipfile.ZipFile(zp) as zf:
            names = zf.namelist()
        assert sc._MIGRATE_MARKER_NAME in names
        assert any(n.endswith("manifest.json") for n in names)
        assert any("plasmid_blobs" in n for n in names)   # sequences included

    def test_hmm_excluded_by_default_but_optional(self, tmp_path):
        _seed()
        hd = getattr(sc._state, "_HMM_DATABASES_DIR", None)
        if not isinstance(hd, sc.Path):
            pytest.skip("no HMM databases dir attr")
        hd.mkdir(parents=True, exist_ok=True)
        (hd / "Pfam.hmm").write_text("FAKE")
        z1 = str(tmp_path / "lean.zip")
        sc._export_migrate_archive(z1)
        with zipfile.ZipFile(z1) as zf:
            assert not any("Pfam.hmm" in n for n in zf.namelist())
        z2 = str(tmp_path / "full.zip")
        sc._export_migrate_archive(z2, include_hmm=True)
        with zipfile.ZipFile(z2) as zf:
            assert any("Pfam.hmm" in n for n in zf.namelist())


class TestMigrateRefusals:
    def test_foreign_zip_refused(self, tmp_path):
        z = str(tmp_path / "f.zip")
        with zipfile.ZipFile(z, "w") as zf:
            zf.writestr("x.txt", "no")
        with pytest.raises(ValueError, match="not a SpliceCraft migrate"):
            sc._import_migrate_archive(z)

    def test_future_format_refused(self, tmp_path):
        z = str(tmp_path / "fu.zip")
        with zipfile.ZipFile(z, "w") as zf:
            zf.writestr(sc._MIGRATE_MARKER_NAME, json.dumps(
                {"format": "splicecraft-migrate", "format_version": 999}))
        with pytest.raises(ValueError, match="newer than this"):
            sc._import_migrate_archive(z)

    def test_path_traversal_refused(self, tmp_path):
        z = str(tmp_path / "ev.zip")
        with zipfile.ZipFile(z, "w") as zf:
            zf.writestr(sc._MIGRATE_MARKER_NAME, json.dumps(
                {"format": "splicecraft-migrate", "format_version": 1,
                 "snapshot_subdir": "data"}))
            zf.writestr("data/../../etc/evil", "pwn")
        with pytest.raises(ValueError, match="(?i)unsafe"):
            sc._import_migrate_archive(z)

    def test_corrupt_zip_refused_data_intact(self, tmp_path):
        _seed()
        before = _fingerprint()
        good = str(tmp_path / "g.zip")
        sc._export_migrate_archive(good)
        bad = str(tmp_path / "bad.zip")
        data = open(good, "rb").read()
        with open(bad, "wb") as f:
            f.write(data[: len(data) // 2])
        with pytest.raises((ValueError, OSError, zipfile.BadZipFile)):
            sc._import_migrate_archive(bad)
        _bust()
        assert _fingerprint() == before    # refused import left data intact

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            sc._import_migrate_archive(str(tmp_path / "nope.zip"))


class TestMigrateRegistryInvariant:
    def test_export_covers_registry(self, tmp_path):
        """Every registry FILE present on disk (minus default exclusions)
        appears in the export manifest — so a new data file added to the
        registry is migrated automatically (the futureproof invariant)."""
        _seed()
        zp = str(tmp_path / "d.zip")
        sc._export_migrate_archive(zp)
        with zipfile.ZipFile(zp) as zf:
            man = json.loads(zf.read("data/manifest.json"))
        manifest_attrs = {f["attr"] for f in man.get("files", [])}
        for attr in sc._USER_DATA_FILE_ATTRS:
            if attr in sc._MIGRATE_DEFAULT_EXCLUDE_ATTRS:
                continue
            p = sc._resolve_state_or_hub(attr)
            if isinstance(p, sc.Path) and p.is_file():
                assert attr in manifest_attrs, (
                    f"{attr} is in the user-data registry but missing from the "
                    "migrate manifest — keep _export_migrate_archive driven by "
                    "the registry.")


# ── audit 2026-10-23: the import is all or nothing, and deletes no blob ─────

class TestMigrateImportIsAllOrNothing:
    """The import replaced the metadata files first and the folders after,
    each on its own. A full disk during the blob copy left the archive's
    collections pointing at blobs the rollback had just removed — every
    imported plasmid without a sequence — and the blob store was replaced
    wholesale although it is content-addressed, so this machine's blobs went
    with it. Now: blobs are MERGED in first (never deleted), the other
    folders replaced next, and the metadata only once every folder landed."""

    def _archive_from_machine_a(self, tmp_path):
        _seed()
        cr = sc._state._CRASH_RECOVERY_DIR
        cr.mkdir(parents=True, exist_ok=True)
        (cr / "a-autosave.gb").write_text(_gb("AS", "ACGT" * 50))
        zp = str(tmp_path / "a.zip")
        sc._export_migrate_archive(zp)
        return zp

    def _machine_b(self):
        _wipe()
        cr = sc._state._CRASH_RECOVERY_DIR
        if cr.is_dir():
            shutil.rmtree(cr)
        sc._save_collections([{"name": "Default", "plasmids": []}])
        sc._commit_library_entry_to_collection({
            "id": "pB", "name": "B plasmid", "size": 800, "n_feats": 0,
            "source": "test", "added": "2026-10-24",
            "gb_text": _gb("PB", "GATC" * 200)}, "Default")
        _bust()
        files = {k: v for k, v in _fingerprint().items()
                 if not k.startswith("blob:")}
        blobs = {p.name for p in sc._plasmid_blob_dir().glob("*.gb")}
        assert blobs, "machine B has no blob to lose"
        return files, blobs

    @staticmethod
    def _fail_copies_into(monkeypatch, live_dir):
        import splicecraft_backup as B
        real_copy2, real_tree = B.shutil.copy2, B.shutil.copytree

        def _inside(dst):
            try:
                return sc.Path(dst).resolve().is_relative_to(live_dir.resolve())
            except OSError:
                return False

        def copy2(src, dst, *a, **k):
            if _inside(dst):
                raise OSError(28, "No space left on device")
            return real_copy2(src, dst, *a, **k)

        def copytree(src, dst, *a, **k):
            if _inside(dst) or sc.Path(dst) == live_dir:
                sc.Path(dst).mkdir(parents=True, exist_ok=True)
                raise OSError(28, "No space left on device")
            return real_tree(src, dst, *a, **k)
        monkeypatch.setattr(B.shutil, "copy2", copy2)
        monkeypatch.setattr(B.shutil, "copytree", copytree)

    def test_a_full_disk_during_the_blob_merge_changes_nothing(
            self, tmp_path, monkeypatch):
        zp = self._archive_from_machine_a(tmp_path)
        files, blobs = self._machine_b()
        self._fail_copies_into(monkeypatch, sc._plasmid_blob_dir())
        summary = sc._import_migrate_archive(zp)
        assert summary["failed"] and not summary["restored_files"], summary
        _bust()
        assert {k: v for k, v in _fingerprint().items()
                if not k.startswith("blob:")} == files
        assert blobs <= {p.name for p in sc._plasmid_blob_dir().glob("*.gb")}

    def test_a_failed_folder_after_the_blobs_replaces_no_metadata(
            self, tmp_path, monkeypatch):
        zp = self._archive_from_machine_a(tmp_path)
        files, blobs = self._machine_b()
        self._fail_copies_into(monkeypatch, sc._state._CRASH_RECOVERY_DIR)
        summary = sc._import_migrate_archive(zp)
        assert summary["failed"] and not summary["restored_files"], summary
        _bust()
        assert {k: v for k, v in _fingerprint().items()
                if not k.startswith("blob:")} == files
        assert blobs <= {p.name for p in sc._plasmid_blob_dir().glob("*.gb")}
        assert not (sc._state._CRASH_RECOVERY_DIR / "a-autosave.gb").exists()

    def test_an_import_keeps_this_machines_blobs(self, tmp_path):
        zp = self._archive_from_machine_a(tmp_path)
        _files, blobs = self._machine_b()
        summary = sc._import_migrate_archive(zp)
        assert not summary["failed"], summary["failed"]
        _bust()
        assert blobs <= {p.name for p in sc._plasmid_blob_dir().glob("*.gb")}
        plasmids = next(c["plasmids"] for c in sc._load_collections()
                        if c.get("name") == "Default")
        assert [p["id"] for p in plasmids] == ["p1"]
        assert "LOCUS" in (sc._rehydrate_entry(plasmids[0]).get("gb_text") or "")

    def test_a_queued_mirror_cannot_land_on_the_import(self, tmp_path):
        """A deferred collection mirror queued before the import used to land
        AFTER it: with the caches busted, it wrote its queued snapshot — this
        machine's old library — into the imported collection of the same
        name. The import transaction lands it first."""
        zp = self._archive_from_machine_a(tmp_path)
        self._machine_b()
        stale = [dict(e) for e in sc._load_library()]
        seq = sc._mark_mirror_pending()
        sc._state._collection_sync_pending = ("Default", stale, seq)
        sc._migrate_import_transaction(zp)
        payload = sc._state._collection_sync_pending
        if payload:                       # what the background drain does next
            sc._land_deferred_mirror(*payload)
        _bust()
        plasmids = next(c["plasmids"] for c in sc._load_collections()
                        if c.get("name") == "Default")
        assert [p["id"] for p in plasmids] == ["p1"]


# ── all-or-nothing really is (audit 2026-11-09) ─────────────────────────────
# `_restore_snapshot_all_or_nothing` staged and sha256-verified every file,
# then replaced them ONE AT A TIME and continued past a failure, never setting
# `aborted`. `plasmid_library.json` and `collections.json` are a MIRROR PAIR,
# so a partial is not "most of the restore" — the next launch resolves the
# mismatch by discarding the half that landed. Measured: the archive's library
# landed, collections failed EACCES, and the imported plasmids were gone at the
# next start. EACCES is the realistic trigger on Windows, a supported
# platform, where `os.replace` fails if the target is held open by a sync
# client, antivirus, or a second SpliceCraft.

def _mig_entry(name: str, seq: str = "ACGT" * 40) -> dict:
    return {"id": f"id-{name}", "name": name, "size": len(seq),
            "n_feats": 0, "gb_text": "", "source": ""}


def _mig_write_generation(tag: str, n: int) -> list:
    entries = [_mig_entry(f"p{i}") for i in range(n)]
    sc._save_collections([{"name": "C", "description": tag,
                           "plasmids": entries}])
    sc._save_library(entries)
    sc._set_active_collection_name("C")
    return entries


def _mig_state() -> tuple:
    lib = sc._load_library()
    col = next((c for c in sc._load_collections()
                if c.get("name") == "C"), {})
    return len(lib), len(col.get("plasmids") or []), col.get("description")


class TestAllOrNothingRestoreIsAtomic:
    def _snapshot_of(self, tag: str, n: int, backup_dir):
        import splicecraft_backup as B
        _mig_write_generation(tag, n)
        return B._create_pre_update_snapshot("1.3.18", backup_dir=backup_dir)

    def test_a_failed_metadata_replace_restores_everything(self, tmp_path):
        import os as _os
        import splicecraft_backup as B
        bk = tmp_path / "bk"
        bk.mkdir()
        snap = self._snapshot_of("ARCHIVE", 8, bk)
        _mig_write_generation("THIS-MACHINE", 3)
        assert _mig_state() == (3, 3, "THIS-MACHINE"), "fixture premise"

        real_replace = _os.replace
        fired: list = []

        def one_shot(src, dst, *a, **k):
            # Fail the FIRST replace of collections.json only, so the
            # rollback is allowed to succeed as it would in reality.
            if str(dst).endswith("collections.json") and not fired:
                fired.append(1)
                raise OSError(13, "Permission denied")
            return real_replace(src, dst, *a, **k)

        # A SEPARATE MonkeyPatch, never the test's own `monkeypatch` +
        # `.undo()`: that would also revert conftest's sandboxing of the data
        # dir (they share it) and point the next save at the user's real
        # library. The L2 chokepoint refuses that — which is how this very
        # test caught itself doing it — but the test would then be proving
        # nothing. Same rule as `test_save_safety.py`.
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(_os, "replace", one_shot)
            summary = B._restore_pre_update_snapshot(snap, backup_dir=bk,
                                                      all_or_nothing=True)
        assert fired, "premise: the injected failure never fired"

        # The caller must be able to TELL. Pre-fix `aborted` was never set.
        assert summary.get("aborted") is True, summary
        assert not summary.get("restored_files"), summary
        # And the generation it found must be intact — a mixed mirror pair is
        # the state the next launch resolves by throwing data away.
        with sc._cache_lock:
            for attr in sc._MASTER_DELETE_CACHE_ATTRS:
                sc._reset_master_delete_cache(attr)
        n_lib, n_col, desc = _mig_state()
        assert (n_lib, n_col, desc) == (3, 3, "THIS-MACHINE"), (
            f"left a mixed generation: library={n_lib} collections={n_col} "
            f"({desc})")

    def test_no_stash_files_are_left_behind(self, tmp_path):
        """The rollback copies each original aside; none of those copies may
        survive a successful restore or an aborted one."""
        import os as _os
        import splicecraft_backup as B
        bk = tmp_path / "bk"
        bk.mkdir()
        snap = self._snapshot_of("ARCHIVE", 8, bk)
        _mig_write_generation("THIS-MACHINE", 3)
        real_replace = _os.replace
        fired: list = []

        def one_shot(src, dst, *a, **k):
            if str(dst).endswith("collections.json") and not fired:
                fired.append(1)
                raise OSError(13, "Permission denied")
            return real_replace(src, dst, *a, **k)

        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(_os, "replace", one_shot)
            B._restore_pre_update_snapshot(snap, backup_dir=bk,
                                            all_or_nothing=True)
        leftovers = sorted(p.name for p in sc._state._DATA_DIR.iterdir()
                            if ".restore-orig" in p.name)
        assert not leftovers, leftovers

    def test_a_clean_restore_still_lands_both_halves(self, tmp_path):
        """Control: the atomicity must not have broken the happy path."""
        import splicecraft_backup as B
        bk = tmp_path / "bk"
        bk.mkdir()
        snap = self._snapshot_of("ARCHIVE", 8, bk)
        _mig_write_generation("THIS-MACHINE", 3)
        summary = B._restore_pre_update_snapshot(snap, backup_dir=bk,
                                                  all_or_nothing=True)
        assert not summary.get("aborted"), summary
        with sc._cache_lock:
            for attr in sc._MASTER_DELETE_CACHE_ATTRS:
                sc._reset_master_delete_cache(attr)
        n_lib, n_col, desc = _mig_state()
        assert (n_lib, n_col, desc) == (8, 8, "ARCHIVE"), (
            f"library={n_lib} collections={n_col} ({desc})")
