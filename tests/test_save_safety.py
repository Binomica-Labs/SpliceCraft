"""
test_save_safety — Save must never destroy what the user already had.

Audit 2026-10-23 drove real workflows through the event loop and compared
memory, panels and disk after every step. These tests pin what it found:

  * Ctrl+S wrote GenBank text over the `.dna` / FASTA file a plasmid was
    opened from, keeping the old name — the user's original was gone.
"""
from __future__ import annotations

import asyncio

import pytest

import splicecraft as sc

TERMINAL_SIZE = (160, 48)


async def _settle(pilot, n=4):
    for _ in range(n):
        await pilot.pause()


async def _wait_workers(app, pilot, timeout=20.0):
    loop = asyncio.get_running_loop()
    end = loop.time() + timeout
    while any(w.is_running for w in app.workers) and loop.time() < end:
        await asyncio.sleep(0.02)
    await _settle(pilot, 2)


def _fasta(tmp_path):
    p = tmp_path / "user_seq.fasta"
    p.write_text(">user_seq some description\n" + "ACGT" * 300 + "\n")
    return p


class TestSaveNeverRewritesAForeignSourceFile:

    async def test_sync_save_leaves_a_fasta_source_untouched(self, tmp_path):
        src = _fasta(tmp_path)
        before = src.read_bytes()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            rec = sc._fasta_path_to_record(str(src))   # OpenFileModal's path
            rec._tui_source = str(src)
            app._import_and_persist(rec)
            await _wait_workers(app, pilot)
            assert app._source_path == str(src)
            assert app._do_save() is True
            await _wait_workers(app, pilot)
        assert src.read_bytes() == before, "Save rewrote the FASTA file"
        ids = [e.get("id") for e in sc._load_library()]
        assert rec.id in ids, "the record must still reach the library"

    async def test_ctrl_s_leaves_a_dna_source_untouched(self, tmp_path):
        # A real SnapGene `.dna`, made by the app's own exporter.
        demo = sc._make_demo_record()
        dna = tmp_path / "user_vector.dna"
        sc._export_commercialsaas_dna(
            {"id": demo.id, "name": demo.name,
             "gb_text": sc._record_to_gb_text(demo)}, str(dna))
        before = dna.read_bytes()
        assert not before.startswith(b"LOCUS")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            rec = sc.load_genbank(str(dna))
            rec._tui_source = str(dna)
            app._import_and_persist(rec)
            await _wait_workers(app, pilot)
            app.action_save()
            await _wait_workers(app, pilot)
        assert dna.read_bytes() == before, "Ctrl+S rewrote the .dna file"

    async def test_a_genbank_source_is_still_written_back(self, tmp_path):
        src = tmp_path / "mine.gb"
        src.write_text(sc._record_to_gb_text(sc._make_demo_record()))
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            rec = sc.load_genbank(str(src))
            rec._tui_source = str(src)
            app._import_and_persist(rec)
            await _wait_workers(app, pilot)
            src.write_text("stale")            # Save must replace this
            assert app._do_save() is True
        assert src.read_text().startswith("LOCUS")


# ═══════════════════════════════════════════════════════════════════════════
# D1 — a kill before the deferred collection mirror lands must not revert
# ═══════════════════════════════════════════════════════════════════════════
# The library file is written durably and its active-collection mirror is
# deferred to a background thread. The startup restore rewrites the library
# FROM the collection unless the mirror-dirty marker is up — and the marker
# used to be raised only when a mirror FAILED. A SIGKILL / OOM / power loss in
# between therefore reverted the save at the next launch, silently (measured
# lost at 0.05 s and 0.15 s after a save of a 2,500-plasmid library).

def _entry(i, seq="ACGT" * 30):
    return {"id": f"p{i}", "name": f"p {i}", "size": len(seq), "n_feats": 0,
            "source": "test", "added": "2026-01-01",
            "gb_text": sc._record_to_gb_text(
                __import__("Bio.SeqRecord", fromlist=["SeqRecord"]).SeqRecord(
                    __import__("Bio.Seq", fromlist=["Seq"]).Seq(seq),
                    id=f"p{i}", name=f"p{i}",
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"}))}


def _active_collection_with(entries):
    sc._save_collections([{"name": "Main", "description": "",
                           "plasmids": [dict(e) for e in entries],
                           "saved": "2026-01-01"}])
    sc._set_active_collection_name("Main")
    sc._save_library([dict(e) for e in entries])
    assert not sc._mirror_is_dirty()


def _production_async_mirror(monkeypatch):
    """conftest forces the mirror synchronous for determinism; these tests
    are ABOUT the deferred path production takes."""
    monkeypatch.setattr(sc, "_collection_sync_force_sync", False)


def _kill_before_the_mirror_lands(monkeypatch):
    """The process 'dies' before the background drain runs: the drain is
    never allowed to write."""
    _production_async_mirror(monkeypatch)
    monkeypatch.setattr(sc, "_drain_collection_sync_loop", lambda: None)
    monkeypatch.setattr(sc._state, "_collection_sync_thread", None)
    monkeypatch.setattr(sc._state, "_collection_sync_pending", None)


@pytest.fixture(autouse=True)
def _reset_collection_sync_shutdown():
    """`_flush_pending_collection_sync` sets a process-wide shutdown flag
    that nothing clears; don't leak it into later tests on this worker."""
    yield
    sc._collection_sync_shutdown.clear()


class TestDeferredMirrorSurvivesAKill:

    def test_save_library_async(self, monkeypatch):
        base = [_entry(i) for i in range(3)]
        _active_collection_with(base)
        _kill_before_the_mirror_lands(monkeypatch)
        sc._save_library(base + [_entry(99)], async_sync=True)
        assert sc._mirror_is_dirty(), "no marker while the mirror is pending"
        # Next launch: the restore must keep the newer library file.
        sc._state._library_cache = None
        sc._restore_library_from_active_collection()
        assert "p99" in {e.get("id") for e in sc._load_library()}

    def test_library_panel_add_worker_path(self, monkeypatch):
        base = [_entry(i) for i in range(3)]
        _active_collection_with(base)
        _kill_before_the_mirror_lands(monkeypatch)
        with sc._cache_lock:
            sc._state._library_cache = [dict(e) for e in base + [_entry(7)]]
        sc.LibraryPanel._write_library_then_defer_mirror(
            sc._state._library_cache)
        assert sc._mirror_is_dirty()
        sc._state._library_cache = None
        sc._restore_library_from_active_collection()
        assert "p7" in {e.get("id") for e in sc._load_library()}

    def test_the_marker_clears_once_the_mirror_lands(self, monkeypatch):
        base = [_entry(i) for i in range(3)]
        _active_collection_with(base)
        _production_async_mirror(monkeypatch)
        sc._save_library(base + [_entry(42)], async_sync=True)
        sc._flush_pending_collection_sync(timeout_s=10.0)
        assert not sc._mirror_is_dirty()
        coll = sc._find_collection("Main")
        assert "p42" in {e.get("id") for e in coll["plasmids"]}

    def test_no_marker_is_left_when_no_collection_receives_it(self):
        """A stray marker would make the next launch flush this library
        into whichever collection becomes active (D10)."""
        sc._save_collections([])
        sc._set_active_collection_name("")
        sc._save_library([_entry(1)], async_sync=True)
        assert not sc._mirror_is_dirty()

    def test_a_late_deferred_mirror_never_lands_older_content(
            self, monkeypatch):
        """D6: the drain wrote the snapshot captured at QUEUE time. A
        synchronous commit landing newer content in between (Constructor /
        Gibson / PCR saves) was erased by it — and the marker then cleared,
        so the next launch reverted the library too."""
        base = [_entry(i) for i in range(3)]
        _active_collection_with(base)
        _kill_before_the_mirror_lands(monkeypatch)          # queue, no drain
        sc._save_library(base + [_entry(50)], async_sync=True)
        monkeypatch.setattr(sc, "_collection_sync_force_sync", True)
        sc._save_library(base + [_entry(50), _entry(51)])    # newer, sync
        sc._flush_pending_collection_sync(timeout_s=5.0)     # stale one lands
        ids = {e.get("id") for e in sc._find_collection("Main")["plasmids"]}
        assert {"p50", "p51"} <= ids, ids


# ═══════════════════════════════════════════════════════════════════════════
# W2 — Save & Quit must persist before it exits
# ═══════════════════════════════════════════════════════════════════════════

def _disk_entry(entry_id):
    entries, _w = sc._safe_load_json(sc._state._LIBRARY_FILE,
                                     "Plasmid library")
    return next((e for e in entries if e.get("id") == entry_id), None)


class TestSaveAndQuitKeepsTheEdits:

    async def test_save_and_quit_writes_the_library_before_exiting(
            self, monkeypatch):
        """`_do_save` handed the library write to a worker and quit exited in
        the same tick, cancelling it — after `_mark_clean` had deleted the
        autosave (lost 8 runs in 10). The worker is made to never run here,
        which is what the exit's cancellation did."""
        rec = sc._make_demo_record()
        n0 = len(rec.seq)
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            app._import_and_persist(rec)
            await _wait_workers(app, pilot)
            monkeypatch.setattr(sc.LibraryPanel, "_add_save_to_disk",
                                lambda self, entries: None)
            app._edit_dialog_result(("GGGGGGG", "insert", 5, 5))
            await _settle(pilot)
            assert app._unsaved
            rid = app._current_record.id
            app._on_quit_response("save")
            await _settle(pilot)
        on_disk = _disk_entry(rid)
        assert on_disk is not None and on_disk.get("size") == n0 + 7, on_disk
        coll = sc._find_collection(sc._get_active_collection_name())
        mirrored = next(e for e in coll["plasmids"] if e.get("id") == rid)
        assert mirrored.get("size") == n0 + 7, "the collection mirror lags"

    def test_exit_net_writes_a_cache_change_that_never_reached_disk(self):
        rec = sc._make_demo_record()
        sc._save_library([_entry(1)])
        assert sc._flush_library_cache_if_unsaved() is False   # in sync
        with sc._cache_lock:
            sc._state._library_cache = [_entry(1), _entry(2)]  # write lost
        assert sc._flush_library_cache_if_unsaved() is True
        assert _disk_entry("p2") is not None
        assert sc._flush_library_cache_if_unsaved() is False


# ═══════════════════════════════════════════════════════════════════════════
# W3 / W12 — the dirty flag must tell the truth
# ═══════════════════════════════════════════════════════════════════════════

class TestDirtyFlagTellsTheTruth:

    async def test_renaming_a_dirty_plasmid_keeps_it_dirty(self):
        """The rename refreshed the title with `_mark_clean`, which dropped
        the unsaved flag: clicking another plasmid then asked nothing and the
        edit was gone (undo empty too)."""
        rec = sc._make_demo_record()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            app._import_and_persist(rec)
            await _wait_workers(app, pilot)
            app._edit_dialog_result(("GGGGGGG", "insert", 5, 5))
            await _settle(pilot)
            assert app._unsaved
            app._rename_library_entry(app._current_record.id,
                                      "Renamed Plasmid")
            await _wait_workers(app, pilot)
            assert app._unsaved, "rename declared unsaved edits saved"
            assert app.title.startswith("SpliceCraft — *Renamed Plasmid")

    async def test_clearing_the_canvas_clears_the_unsaved_marks(self):
        rec = sc._make_demo_record()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            app._import_and_persist(rec)
            await _wait_workers(app, pilot)
            app._edit_dialog_result(("GGG", "insert", 5, 5))
            await _settle(pilot)
            assert "*" in app.title
            app._clear_canvas()
            await _settle(pilot)
            assert "*" not in app.title


# ═══════════════════════════════════════════════════════════════════════════
# W4 / W5 — a plasmid belongs to the collection it was opened in
# ═══════════════════════════════════════════════════════════════════════════

def _same_id_entry(seq):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    r = SeqRecord(Seq(seq), id="pSame", name="pSame",
                  annotations={"molecule_type": "DNA",
                               "topology": "circular"})
    return {"id": "pSame", "name": "pSame", "size": len(seq), "n_feats": 0,
            "source": "test", "added": "2026-01-01",
            "gb_text": sc._record_to_gb_text(r)}


async def _switch(app, pilot, name):
    lib = app.query_one("#library", sc.LibraryPanel)
    lib._commit_collection_switch(name, sc._find_collection(name))
    await _wait_workers(app, pilot)
    await _settle(pilot, 4)


class TestCollectionIdentity:

    async def test_each_collection_shows_its_own_copy_of_a_shared_id(self):
        seq_a, seq_b = "ACGT" * 75, "ACGT" * 75 + "GATTACA"
        sc._save_collections([
            {"name": "A", "description": "", "saved": "2026-01-01",
             "plasmids": [_same_id_entry(seq_a)]},
            {"name": "B", "description": "", "saved": "2026-01-01",
             "plasmids": [_same_id_entry(seq_b)]}])
        sc._set_active_collection_name("A")
        sc._save_library([_same_id_entry(seq_a)])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _switch(app, pilot, "B")
            assert len(app._current_record.seq) == len(seq_b)
            await _switch(app, pilot, "A")
            assert len(app._current_record.seq) == len(seq_a), (
                "collection A is showing collection B's copy")

    async def test_ctrl_s_saves_into_the_collection_it_was_pressed_in(
            self, monkeypatch):
        """The worker resolved the collection (and the record) when it ran:
        a switch landing first filed the save into the new collection."""
        seq_a, seq_b = "ACGT" * 75, "TTGCA" * 60
        other = dict(_same_id_entry(seq_b), id="pOther", name="pOther")
        sc._save_collections([
            {"name": "A", "description": "", "saved": "2026-01-01",
             "plasmids": [_same_id_entry(seq_a)]},
            {"name": "B", "description": "", "saved": "2026-01-01",
             "plasmids": [other]}])
        sc._set_active_collection_name("A")
        sc._save_library([_same_id_entry(seq_a)])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _switch(app, pilot, "A")
            app._edit_dialog_result(("GGGGGGG", "insert", 5, 5))
            await _settle(pilot)
            # Hold the worker until the switch has landed — the race.
            real = sc.PlasmidApp._save_worker
            held = []
            monkeypatch.setattr(sc.PlasmidApp, "_save_worker",
                                lambda self, *a: held.append(a))
            app.action_save()
            await _switch(app, pilot, "B")
            monkeypatch.setattr(sc.PlasmidApp, "_save_worker", real)
            app._save_worker(*held[0])
            await _wait_workers(app, pilot)
        a_ids = {e["id"]: e for e in sc._find_collection("A")["plasmids"]}
        b_ids = {e["id"] for e in sc._find_collection("B")["plasmids"]}
        assert a_ids["pSame"]["size"] == len(seq_a) + 7, "edit not in A"
        assert "pSame" not in b_ids, "the save was filed into B"


# ═══════════════════════════════════════════════════════════════════════════
# AG2 — an export can never land on the app's own storage
# ═══════════════════════════════════════════════════════════════════════════

class TestExportsStayOutOfTheDataDir:

    def _targets(self):
        blobs = sc._state._DATA_DIR / "plasmid_blobs"
        blobs.mkdir(parents=True, exist_ok=True)
        victim = blobs / ("ab" * 32 + ".gb")
        victim.write_text("STORED PLASMID — MUST SURVIVE")
        return victim

    def test_every_export_engine_refuses(self):
        import pytest
        import splicecraft_mapimage as mi
        victim = self._targets()
        rec = sc._make_demo_record()
        entry = {"id": rec.id, "name": rec.name,
                 "gb_text": sc._record_to_gb_text(rec)}
        calls = [
            lambda: sc._export_genbank_to_path(rec, victim),
            lambda: sc._export_gff_to_path(rec, victim),
            lambda: sc._export_embl_to_path(rec, victim),
            lambda: sc._export_fasta_to_path("x", "ACGT", victim),
            lambda: sc._export_commercialsaas_dna(entry, victim),
            lambda: sc._export_primers_to_csv([], victim),
            lambda: mi.export_plasmid_map(victim, feats=[], total=100,
                                          title="x"),
        ]
        for call in calls:
            with pytest.raises(OSError, match="data folder"):
                call()
        assert victim.read_text() == "STORED PLASMID — MUST SURVIVE"

    def test_agent_write_check_refuses(self):
        import splicecraft_agent as sca
        from pathlib import Path
        victim = self._targets()
        assert "data folder" in (sca._check_agent_write_path(Path(victim))
                                 or "")

    def test_exports_elsewhere_still_work(self, tmp_path):
        # (tmp_path IS the data dir in tests — a user file at its top level,
        # or in a folder of the user's own, is not managed storage.)
        rec = sc._make_demo_record()
        for out in (tmp_path / "ok.gb", tmp_path / "exports" / "ok2.gb"):
            out.parent.mkdir(exist_ok=True)
            sc._export_genbank_to_path(rec, out)
            assert out.read_text().startswith("LOCUS")


# ═══════════════════════════════════════════════════════════════════════════
# U3 — a failing background task must not end the session
# ═══════════════════════════════════════════════════════════════════════════

class TestAWorkerErrorDoesNotEndTheSession:

    async def test_autosave_worker_error_is_reported_not_fatal(
            self, monkeypatch):
        """Textual's default `exit_on_error=True` ended the WHOLE app when an
        exception escaped a worker — here the autosave worker's deepcopy,
        which sits outside its try (the error a concurrent in-place mutation
        raises). The app's own log + toast net only ever ran after the
        shutdown had begun, and the unsaved edit went with the app."""
        notes = []
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            app._apply_record(sc._make_demo_record())
            await _wait_workers(app, pilot)
            real_notify = app.notify
            monkeypatch.setattr(
                app, "notify",
                lambda m, **k: notes.append(str(m)) or real_notify(m, **k))
            app._edit_dialog_result(("GGG", "insert", 5, 5))
            await _settle(pilot)

            def boom(*a, **k):
                raise RuntimeError("dictionary changed size during iteration")
            # A SEPARATE MonkeyPatch: `monkeypatch.undo()` would also revert
            # conftest's sandboxing of the data files (they share it).
            with pytest.MonkeyPatch.context() as mp:
                mp.setattr(sc, "deepcopy", boom)
                w = app._do_autosave()
                for _ in range(60):
                    if w.is_finished:
                        break
                    await asyncio.sleep(0.05)
            await _settle(pilot, 3)
            assert app.is_running, "a worker error ended the session"
            assert app._unsaved, "the unsaved edit must still be there"
        assert app.return_code in (0, None)
        assert any("background task failed" in n for n in notes)


# ═══════════════════════════════════════════════════════════════════════════
# D3 — an unseen crash autosave is never pruned
# ═══════════════════════════════════════════════════════════════════════════

class TestCrashRecoveryNeverDiscardsUnseenWork:

    class _App:
        def __init__(self):
            self.notes = []

        def notify(self, msg, **kw):
            self.notes.append(str(msg))

    def _autosave(self, name, age_days):
        import os
        import time
        d = sc._state._CRASH_RECOVERY_DIR
        d.mkdir(parents=True, exist_ok=True)
        p = d / f"{name}-0123456789ab.gb"
        p.write_text(sc._record_to_gb_text(sc._make_demo_record()))
        t = time.time() - age_days * 86400
        os.utime(p, (t, t))
        return p

    def test_an_old_unseen_autosave_is_shown_not_deleted(self):
        """Files older than 30 days were unlinked BEFORE the notice was
        built: a crash whose user didn't relaunch within a month lost the
        only copy of the unsaved edits, without ever being told."""
        old = self._autosave("pMyUnsavedWork", 31)
        app = self._App()
        sc.PlasmidApp._check_crash_recovery(app)
        assert old.exists()
        assert any("pMyUnsavedWork" in n for n in app.notes)

    def test_a_shown_old_autosave_is_still_pruned(self):
        old = self._autosave("pSeenBefore", 31)
        sc._set_setting("crash_recovery_seen",
                        [f"{old.name}|{int(old.stat().st_mtime)}"])
        sc.PlasmidApp._check_crash_recovery(self._App())
        assert not old.exists()


# ═══════════════════════════════════════════════════════════════════════════
# D4 — a file that cannot be READ is not a corrupt file
# ═══════════════════════════════════════════════════════════════════════════

class TestUnreadableIsNotCorrupt:

    def _three_saves(self):
        base = [_entry(i) for i in range(6)]
        sc._save_library(base)
        sc._save_library(base + [_entry(100)])
        sc._save_library(base + [_entry(100), _entry(200)])   # newest
        sc._state._library_cache = None
        return sc._state._LIBRARY_FILE

    def _flaky(self, mp, path, times):
        import errno
        import pathlib
        real = pathlib.Path.read_text
        left = {"n": times}

        def read_text(self, *a, **k):
            if self == path and left["n"] > 0:
                left["n"] -= 1
                raise OSError(errno.EIO, "Input/output error")
            return real(self, *a, **k)
        mp.setattr(pathlib.Path, "read_text", read_text)

    def test_a_transient_read_error_loads_the_newest_data(self, monkeypatch):
        """One EIO on a VALID library renamed it aside as `.corrupt-<ts>`,
        rewrote it from the one-save-older `.bak` and called it corrupt."""
        lib = self._three_saves()
        self._flaky(monkeypatch, lib, 1)
        entries, warn = sc._safe_load_json(lib, "Plasmid library")
        assert "p200" in {e.get("id") for e in entries}
        assert not list(lib.parent.glob(lib.name + ".corrupt-*"))

    def test_a_persistent_read_error_touches_nothing(self):
        lib = self._three_saves()
        before = lib.read_bytes()
        # A SEPARATE MonkeyPatch, never `monkeypatch.undo()`: that would also
        # revert conftest's sandboxing of the data files (they share it).
        with pytest.MonkeyPatch.context() as mp:
            self._flaky(mp, lib, 99)
            entries, warn = sc._safe_load_json(lib, "Plasmid library")
        assert entries == [] and "could not be read" in (warn or "")
        assert not list(lib.parent.glob(lib.name + ".corrupt-*"))
        assert lib.read_bytes() == before
        # The session never saw the contents, so it may not replace them.
        with pytest.raises(OSError, match="could not be read"):
            sc._safe_save_json(lib, [_entry(1)], "Plasmid library")
        sc._safe_load_json(lib, "Plasmid library")   # readable again
        sc._safe_save_json(lib, [_entry(i) for i in range(6)] +
                           [_entry(100), _entry(200)], "Plasmid library")


def test_no_test_undoes_the_shared_monkeypatch():
    """conftest's `_protect_user_data` sandboxes every data file through the
    test's OWN `monkeypatch`, so `monkeypatch.undo()` in a test reverts the
    sandbox mid-test: the data-file paths point at the user's REAL data dir
    again (only the write-authorisation flag reverting with them kept writes
    off it). Use `with pytest.MonkeyPatch.context() as mp:` for a patch that
    must end early (audit 2026-10-23)."""
    import pathlib
    import re
    hits = []
    for f in sorted(pathlib.Path(__file__).parent.glob("test_*.py")):
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if re.match(r"\s*monkeypatch\.undo\(\)\s*(#.*)?$", line):
                hits.append(f"{f.name}:{i}")
    assert not hits, hits


# ═══════════════════════════════════════════════════════════════════════════
# Audit 2026-10-23, round 2 — workflow findings W6–W13 + NEW-1, each driven
# through the real app path. Every `test_*` below failed on the pre-fix tree;
# the `test_guard_*` ones pin behaviour that was already right.
# ═══════════════════════════════════════════════════════════════════════════

def _plain_entry(rid, seq, *, name=None, **extra):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    r = SeqRecord(Seq(seq), id=rid, name=rid,
                  annotations={"molecule_type": "DNA", "topology": "circular"})
    e = {"id": rid, "name": name or rid, "size": len(seq), "n_feats": 0,
         "source": "test", "added": "2026-01-01",
         "gb_text": sc._record_to_gb_text(r)}
    e.update(extra)
    return e


def _seed(colls, active):
    sc._save_collections([
        {"name": n, "description": "", "saved": "2026-01-01",
         "plasmids": [dict(e) for e in es]} for n, es in colls.items()])
    sc._set_active_collection_name(active)
    sc._save_library([dict(e) for e in colls[active]])


def _disk_library():
    entries, _w = sc._safe_load_json(sc._state._LIBRARY_FILE, "Plasmid library")
    return [e for e in entries if isinstance(e, dict)]


def _disk_coll(name):
    colls, _w = sc._safe_load_json(sc._state._COLLECTIONS_FILE, "Collections")
    c = next((c for c in colls if isinstance(c, dict) and c.get("name") == name), {})
    return {e.get("id"): e for e in (c.get("plasmids") or []) if isinstance(e, dict)}


def _lib(app):
    return app.query_one("#library", sc.LibraryPanel)


async def _open_row(app, pilot, rid):
    _lib(app).post_message(sc.LibraryPanel.PlasmidLoad(
        sc._find_library_entry_by_id(rid)))
    await _settle(pilot, 4)


# ═══════════════════════════════════════════════════════════════════════════
# W6 — moving / copying the loaded plasmid
# ═══════════════════════════════════════════════════════════════════════════

def _seed_move():
    _seed({"Main": [_plain_entry("aMove", "ACGT" * 150),
                    _plain_entry("zOther", "TTGCA" * 100)],
           "B": [_plain_entry("pB", "GATC" * 120)]}, "Main")


class TestMovingTheLoadedPlasmid:

    async def test_a_dirty_loaded_plasmid_is_offered_for_saving_first(self):
        """No prompt: the copy that landed in B lacked the edit, and the next
        Ctrl+S re-created the plasmid in Main with it — two diverged copies."""
        _seed_move()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            assert app._current_record.id == "aMove"
            app._edit_dialog_result(("GGGGG", "insert", 5, 5))
            await _settle(pilot)
            _lib(app).action_move_marked()          # `m` on the cursor row
            await _settle(pilot, 6)
            assert isinstance(app.screen, sc.UnsavedNavigateModal), \
                type(app.screen).__name__
            app.screen.query_one("#btn-navunsv-save").press()
            await _settle(pilot, 6)
            await _wait_workers(app, pilot)
            assert isinstance(app.screen, sc.MoveCopyToCollectionModal)
            app.screen.query_one("#btn-movecopy-go").press()   # B: only target
            await _settle(pilot, 6)
            await _wait_workers(app, pilot)
        assert _disk_coll("B")["aMove"]["size"] == 605
        assert "aMove" not in _disk_coll("Main")

    async def test_ctrl_s_after_the_move_does_not_refile_it_in_the_source(self):
        """A CLEAN move, then the user keeps working on what the canvas shows:
        Ctrl+S filed it back into Main beside the moved copy."""
        _seed_move()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            _lib(app).action_move_marked()
            await _settle(pilot, 6)
            app.screen.query_one("#btn-movecopy-go").press()
            await _settle(pilot, 6)
            await _wait_workers(app, pilot)
            if app._current_record is not None:     # (the fix may clear it)
                app._edit_dialog_result(("GGGGG", "insert", 5, 5))
                await _settle(pilot)
                app.action_save()
                await _wait_workers(app, pilot)
        assert "aMove" not in _disk_coll("Main"), "re-created in the source"
        assert _disk_coll("B")["aMove"]["size"] == 600


# ═══════════════════════════════════════════════════════════════════════════
# W7 — .dna export of a canvas with unsaved edits
# ═══════════════════════════════════════════════════════════════════════════

async def _export_dna(app, pilot, out_dir, fname):
    app.action_export_commercialsaas()
    await _settle(pilot)
    scr = app.screen
    assert isinstance(scr, sc.ExportCommercialSaaSModal), type(scr).__name__
    scr._selected_dir = str(out_dir)
    scr.query_one("#export-filename").value = fname
    scr.query_one("#btn-export").press()
    await _settle(pilot)
    return str(sc.load_genbank(str(out_dir / fname)).seq).upper()


class TestDnaExportShowsTheCanvas:

    async def test_unsaved_edits_reach_the_exported_dna(self, tmp_path):
        out = tmp_path / "exports"
        out.mkdir()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            app._import_and_persist(sc._make_demo_record())
            await _wait_workers(app, pilot)
            app._edit_dialog_result(("GGATCC", "insert", 5, 5))
            await _settle(pilot)
            want = str(app._current_record.seq).upper()
            got = await _export_dna(app, pilot, out, "x.dna")
        assert got == want, f"exported {len(got)} bp, canvas {len(want)} bp"

    async def test_unsaved_edits_reach_an_export_of_a_dna_import(self, tmp_path):
        """Sidecar (splice-mode) entry: the stored sidecar no longer matches the
        canvas, so the export must not hand back the original bytes."""
        out = tmp_path / "exports"
        out.mkdir()
        demo = sc._make_demo_record()
        src = tmp_path / "in" / "vector.dna"
        src.parent.mkdir()
        sc._export_commercialsaas_dna({"id": demo.id, "name": demo.name,
                                       "gb_text": sc._record_to_gb_text(demo)},
                                      str(src))
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            r = sc.load_genbank(str(src))
            r._tui_source = str(src)
            app._import_and_persist(r)
            await _wait_workers(app, pilot)
            app._edit_dialog_result(("TTTTTTTT", "insert", 20, 20))
            await _settle(pilot)
            want = str(app._current_record.seq).upper()
            got = await _export_dna(app, pilot, out, "y.dna")
        assert got == want, f"exported {len(got)} bp, canvas {len(want)} bp"


# ═══════════════════════════════════════════════════════════════════════════
# W8 — an import declined at the name-collision prompt
# ═══════════════════════════════════════════════════════════════════════════

def _seed_collision(tmp_path):
    _seed({"Main": [_plain_entry("aBase", "ACGT" * 150),
                    _plain_entry("pColl", "TTGCA" * 100)]}, "Main")
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    p = tmp_path / "in" / "pColl.gb"
    p.parent.mkdir(exist_ok=True)
    p.write_text(sc._record_to_gb_text(SeqRecord(
        Seq("GATC" * 130), id="pColl", name="pColl",
        annotations={"molecule_type": "DNA", "topology": "circular"})))
    return p


async def _open_colliding_file(app, pilot, path, button):
    r = sc.load_genbank(str(path))
    r._tui_source = str(path)                  # what OpenFileModal hands back
    app._on_open_file_done(r)
    await _settle(pilot)
    assert isinstance(app.screen, sc.NameCollisionModal), type(app.screen).__name__
    app.screen.query_one(button).press()
    await _settle(pilot)
    await _wait_workers(app, pilot)


class TestADeclinedImportLeavesTheCanvasAlone:

    @pytest.mark.parametrize("button", ["#btn-namecoll-keep",
                                        "#btn-namecoll-cancel"])
    async def test_keep_or_cancel_restores_the_previous_canvas(
            self, tmp_path, button):
        path = _seed_collision(tmp_path)
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            assert app._current_record.id == "aBase"
            await _open_colliding_file(app, pilot, path, button)
            assert app._current_record.id == "aBase", (
                f"the declined import is still on the canvas "
                f"({app._current_record.id}, unsaved={app._unsaved})")
            app.action_save()                       # Ctrl+S on what is shown
            await _wait_workers(app, pilot)
        names = [e.get("name") for e in _disk_library()]
        assert names.count("pColl") == 1, names

    async def test_guard_overwrite_still_imports(self, tmp_path):
        path = _seed_collision(tmp_path)
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            await _open_colliding_file(app, pilot, path, "#btn-namecoll-overwrite")
            assert len(app._current_record.seq) == 520
        lib = [e for e in _disk_library() if e.get("name") == "pColl"]
        assert len(lib) == 1 and lib[0]["size"] == 520


# ═══════════════════════════════════════════════════════════════════════════
# W10 — the prompt's "Save" keeps what Ctrl+S keeps
# ═══════════════════════════════════════════════════════════════════════════

_META = {"source": "constructor:moclo_plant:backbone", "added": "2020-01-01",
         "kind": "amplicon", "x_future_field": {"keep": "me"}}


class TestBothSavePathsKeepEntryMetadata:

    @pytest.mark.parametrize("how", ["ctrl_s", "prompt_save"])
    async def test_a_save_keeps_what_it_does_not_own(self, how):
        e = _plain_entry("aGib", "ACGT" * 150, **_META)
        _seed({"Main": [e, _plain_entry("zOther", "TTGCA" * 100)]}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            assert app._current_record.id == "aGib"
            app._edit_dialog_result(("GGGGG", "insert", 5, 5))
            await _settle(pilot)
            if how == "ctrl_s":
                app.action_save()
            else:
                await _open_row(app, pilot, "zOther")
                assert isinstance(app.screen, sc.UnsavedNavigateModal)
                app.screen.query_one("#btn-navunsv-save").press()
                await _settle(pilot)
            await _wait_workers(app, pilot)
        saved = _disk_entry("aGib")
        assert saved["size"] == 605
        assert {k: saved.get(k) for k in _META} == _META

    async def test_prompt_save_keeps_the_delete_cascade_grammar_scoped(self):
        e = _plain_entry("aGib", "ACGT" * 150, **_META)
        _seed({"Main": [e, _plain_entry("zOther", "TTGCA" * 100)]}, "Main")
        sc._save_parts_bin([
            {"name": "aGib", "grammar": "moclo_plant", "type": "CDS",
             "sequence": "ATGAAATAA"},
            {"name": "aGib", "grammar": "gb_l0", "type": "CDS",
             "sequence": "ATGCCCTAA"}])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            app._edit_dialog_result(("GGGGG", "insert", 5, 5))
            await _settle(pilot)
            await _open_row(app, pilot, "zOther")
            app.screen.query_one("#btn-navunsv-save").press()
            await _settle(pilot)
            await _wait_workers(app, pilot)
            t = _lib(app).query_one("#lib-table")
            t.move_cursor(row=[k.value for k in t.rows.keys()].index("aGib"))
            await _settle(pilot, 2)
            _lib(app)._request_plasmid_delete()
            await _settle(pilot)
            app.screen.dismiss(True)                 # confirm the delete
            await _settle(pilot)
            await _wait_workers(app, pilot)
            left = [(p.get("name"), p.get("grammar")) for p in sc._load_parts_bin()]
        assert ("aGib", "gb_l0") in left, left


# ═══════════════════════════════════════════════════════════════════════════
# W11 — per-plasmid undo history survives save + switch away + back
# ═══════════════════════════════════════════════════════════════════════════

def _seed_two():
    _seed({"Main": [_plain_entry("aOne", "ACGT" * 150),
                    _plain_entry("zTwo", "TTGCA" * 100)]}, "Main")


class TestUndoHistoryFollowsTheSavedPlasmid:

    async def test_prompt_save_then_switch_back_keeps_undo(self):
        _seed_two()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            for i in range(3):
                app._edit_dialog_result(("G" * (i + 2), "insert", 5, 5))
                await _settle(pilot, 2)
            await _open_row(app, pilot, "zTwo")
            app.screen.query_one("#btn-navunsv-save").press()
            await _settle(pilot)
            await _wait_workers(app, pilot)
            await _open_row(app, pilot, "aOne")
            assert len(app._undo_stack) == 3, len(app._undo_stack)

    async def test_ctrl_s_then_switch_back_keeps_undo(self):
        _seed_two()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            app._edit_dialog_result(("GG", "insert", 5, 5))
            await _settle(pilot, 2)
            app.action_save()
            await _wait_workers(app, pilot)
            await _open_row(app, pilot, "zTwo")
            await _open_row(app, pilot, "aOne")
            assert len(app._undo_stack) == 1, len(app._undo_stack)

    async def test_rename_then_switch_back_keeps_undo(self):
        """The stash key carries the record id, and a rename changes the id:
        the history was filed under the old one and never found again."""
        _seed_two()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            app._edit_dialog_result(("GG", "insert", 5, 5))
            await _settle(pilot, 2)
            app.action_save()
            await _wait_workers(app, pilot)
            app._rename_library_entry("aOne", "aOne renamed")
            await _wait_workers(app, pilot)
            await _open_row(app, pilot, "zTwo")
            await _open_row(app, pilot, "aOne_renamed")   # id follows the name
            assert len(app._undo_stack) == 1, len(app._undo_stack)

    async def test_guard_unsaved_edits_still_come_back_as_redo(self):
        """The H3 net: an unprompted switch (agent load) parks the live state
        as a redo step under the STORED content's key. Must keep working."""
        _seed_two()
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            app._edit_dialog_result(("GG", "insert", 5, 5))
            await _settle(pilot, 2)
            app._apply_record(sc._gb_text_to_record(
                sc._find_library_entry_by_id("zTwo")["gb_text"]))
            await _settle(pilot, 2)
            await _open_row(app, pilot, "aOne")
            assert (len(app._undo_stack), len(app._redo_stack)) == (1, 1)


# ═══════════════════════════════════════════════════════════════════════════
# W13 — a rename past the 200-character cap
# ═══════════════════════════════════════════════════════════════════════════

async def _rename(app, pilot, rid, new_name):
    t = _lib(app).query_one("#lib-table")
    t.move_cursor(row=[k.value for k in t.rows.keys()].index(rid))
    await _settle(pilot, 2)
    _lib(app).query_one("#btn-lib-rename").press()
    await _settle(pilot)
    assert isinstance(app.screen, sc.RenamePlasmidModal)
    app.screen.query_one("#rename-input").value = new_name
    app.screen.query_one("#btn-rename-save").press()
    await _settle(pilot)
    await _wait_workers(app, pilot)


class TestRenameNeverSilentlyChangesTheName:

    async def test_an_over_long_name_is_refused_not_cut(self):
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            await _rename(app, pilot, "aOne", "pBig " + "x" * 336)
            still_open = isinstance(app.screen, sc.RenamePlasmidModal)
        stored = [e.get("name") for e in _disk_library()]
        assert still_open and stored == ["aOne"], (still_open, [len(n) for n in stored])

    async def test_the_collision_check_sees_the_name_that_would_be_stored(self):
        long_name = "A" * 200
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150),
                        _plain_entry("zLong", "TTGCA" * 100, name=long_name)]},
              "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            await _rename(app, pilot, "aOne", long_name + "B" * 141)
        names = [e.get("name") for e in _disk_library()]
        assert names.count(long_name) == 1, "two plasmids share one name"

    async def test_a_control_character_cannot_sneak_a_duplicate_name_past(self):
        """The check compared the TYPED text; `_rename_library_entry` then
        stripped the control character and stored a name another plasmid
        already has (pasted text can carry one)."""
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150),
                        _plain_entry("pDup", "TTGCA" * 100)]}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            await _rename(app, pilot, "aOne", "pDup\x07")
        names = [e.get("name") for e in _disk_library()]
        assert names.count("pDup") == 1, names


# ═══════════════════════════════════════════════════════════════════════════
# NEW-1 — Ctrl+S must not write the library's id into the user's own file
# ═══════════════════════════════════════════════════════════════════════════

def _snapgene_style_gb(seq):
    lines = [f"LOCUS       pSnap{len(seq):>23} bp    DNA     circular SYN 01-JAN-2026",
             "DEFINITION  pSnap.", "ACCESSION   .", "VERSION     .",
             "KEYWORDS    .", "SOURCE      synthetic DNA construct",
             "  ORGANISM  synthetic DNA construct", "            other sequences.",
             "FEATURES             Location/Qualifiers",
             f"     source          1..{len(seq)}",
             '                     /mol_type="other DNA"', "ORIGIN"]
    s = seq.lower()
    for i in range(0, len(s), 60):
        chunk = " ".join(s[i + j:i + j + 10] for j in range(0, 60, 10)
                         if s[i + j:i + j + 10])
        lines.append(f"{i + 1:>9} {chunk}")
    return "\n".join(lines + ["//"]) + "\n"


async def test_ctrl_s_keeps_the_files_own_accession(tmp_path):
    """`VERSION .` is a placeholder id, so the import re-derives a library id
    (`pSnap`); Ctrl+S then wrote it into the user's file as `ACCESSION pSnap`
    / `VERSION pSnap` — an accession that file never had."""
    _seed({"Main": [_plain_entry("aBase", "TTGCA" * 100)]}, "Main")
    path = tmp_path / "pSnap.gb"
    path.write_text(_snapgene_style_gb("ACGT" * 150))
    app = sc.PlasmidApp()
    async with app.run_test(size=TERMINAL_SIZE) as pilot:
        await _settle(pilot)
        await _wait_workers(app, pilot)
        r = sc.load_genbank(str(path))
        r._tui_source = str(path)
        app._on_open_file_done(r)
        await _settle(pilot)
        await _wait_workers(app, pilot)
        assert app._current_record.id not in (".", "")
        app._edit_dialog_result(("GGG", "insert", 5, 5))
        await _settle(pilot)
        app.action_save()
        await _wait_workers(app, pilot)
    text = path.read_text()
    head = [ln for ln in text.splitlines()
            if ln.startswith(("ACCESSION", "VERSION"))]
    assert head == ["ACCESSION   .", "VERSION     ."], head
    assert sc.load_genbank(str(path)).seq.upper().startswith("ACGTAGGG")


# ═══════════════════════════════════════════════════════════════════════════
# `v` (map view): no disk write per press, and Ctrl+S keeps the toggle
# (audit 2026-10-23: each press rewrote library + collections on the UI thread
# and burned a timestamped backup slot in both; Ctrl+S then REVERTED it)
# ═══════════════════════════════════════════════════════════════════════════

def _data_dir_state():
    d = sc._state._DATA_DIR
    files = sorted(str(p.relative_to(d)) for p in d.rglob("*") if p.is_file())
    lib = sc._state._LIBRARY_FILE
    return files, lib.stat().st_mtime_ns


class TestMapViewToggle:

    async def test_v_writes_nothing_and_ctrl_s_keeps_the_toggle(self):
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150, map_mode="circular"),
                        _plain_entry("zTwo", "TTGCA" * 100)]}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            await _open_row(app, pilot, "aOne")
            await _wait_workers(app, pilot)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            assert pm._map_mode == "circular"
            before = _data_dir_state()
            pm.action_toggle_map_view()
            await _settle(pilot)
            await _wait_workers(app, pilot)
            assert pm._map_mode == "linear"
            assert _data_dir_state() == before, "a view toggle wrote to disk"
            app._edit_dialog_result(("GG", "insert", 5, 5))
            await _settle(pilot)
            app.action_save()
            await _wait_workers(app, pilot)
        assert _disk_entry("aOne")["map_mode"] == "linear"

    async def test_the_exit_flush_persists_a_toggle_nothing_else_saved(self):
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150, map_mode="circular"),
                        _plain_entry("zTwo", "TTGCA" * 100)]}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            await _open_row(app, pilot, "aOne")
            app.query_one("#plasmid-map", sc.PlasmidMap).action_toggle_map_view()
            await _settle(pilot)
            await _wait_workers(app, pilot)
        assert sc._flush_library_cache_if_unsaved() is True
        assert _disk_entry("aOne")["map_mode"] == "linear"

    async def test_launch_opens_the_first_entry_in_its_stored_view(self):
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150, map_mode="linear"),
                        _plain_entry("zTwo", "TTGCA" * 100)]}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            assert app._current_record.id == "aOne"
            assert app.query_one("#plasmid-map", sc.PlasmidMap)._map_mode == "linear"


# ═══════════════════════════════════════════════════════════════════════════
# Ctrl+S redraws only the saved row (audit 2026-10-23: rebuilding all rows
# cost ~0.3 s per save at 2,000 entries)
# ═══════════════════════════════════════════════════════════════════════════

def _row_cells_text(t, key):
    row = t.get_row_index(key)
    return [str(c) for c in t.get_row_at(row)]


class TestSaveRedrawsOneRow:

    async def test_ctrl_s_updates_the_row_in_place(self, monkeypatch):
        entries = [_plain_entry(f"p{i:03d}", "ACGT" * (100 + i)) for i in range(40)]
        _seed({"Main": entries}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            await _open_row(app, pilot, "p020")
            await _wait_workers(app, pilot)
            t = _lib(app).query_one("#lib-table", sc.DataTable)
            cursor_before = t.cursor_row
            clears = [0]
            real_clear = t.clear

            def counting_clear(*a, **k):
                clears[0] += 1
                return real_clear(*a, **k)

            with pytest.MonkeyPatch.context() as mp:
                mp.setattr(t, "clear", counting_clear)
                app._edit_dialog_result(("GGGGGGG", "insert", 5, 5))
                await _settle(pilot)
                app.action_save()
                await _wait_workers(app, pilot)
                await _settle(pilot)
            assert clears[0] == 0, "Ctrl+S rebuilt the whole table"
            assert t.cursor_row == cursor_before
            in_place = _row_cells_text(t, "p020")
            assert "487" in in_place[4]        # 480 + 7 bp, shown in place
            _lib(app)._repopulate_plasmids()  # what a full rebuild draws
            await _settle(pilot)
            assert _row_cells_text(t, "p020") == in_place

    async def test_a_new_plasmid_still_repopulates(self):
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            await _wait_workers(app, pilot)
            rec = sc._make_demo_record()
            rec.id = rec.name = "zNew"
            app._apply_record(rec)
            await _settle(pilot)
            app._mark_dirty()
            app.action_save()
            await _wait_workers(app, pilot)
            await _settle(pilot)
            t = _lib(app).query_one("#lib-table", sc.DataTable)
            keys = [k.value for k in t.rows.keys()]
            assert "zNew" in keys and "aOne" in keys


# ═══════════════════════════════════════════════════════════════════════════
# Data-dir hygiene (audit 2026-10-23: D5 snapshot order, D7 symlink, D9 .tmp)
# ═══════════════════════════════════════════════════════════════════════════

def test_orphan_tmp_files_are_swept_from_blobs_and_crash_recovery():
    import os
    import time as _time
    old = _time.time() - (sc._ORPHAN_TMP_MIN_AGE_S + 10)
    made = []
    for d in (sc._plasmid_blob_dir(), sc._state._CRASH_RECOVERY_DIR):
        d.mkdir(parents=True, exist_ok=True)
        f = d / ".deadbeef.gb.abc123.tmp"
        f.write_text("half-written")
        os.utime(f, (old, old))
        made.append(f)
    sc._sweep_orphan_tmp_files(sc._state._DATA_DIR)
    assert not any(f.exists() for f in made), made


def test_a_symlinked_data_file_is_named_and_the_launch_refused(
        monkeypatch, tmp_path, capsys):
    import sys as _sys
    real = tmp_path / "elsewhere.json"
    real.write_text('{"_schema_version": 1, "entries": []}')
    coll = sc._state._COLLECTIONS_FILE
    if coll.exists():
        coll.unlink()
    coll.symlink_to(real)
    monkeypatch.setattr(_sys, "argv", ["splicecraft"])
    # main() refuses a terminal under 100x30 before it reaches the guard;
    # report a wide one so an 80x24 shell tests the symlink, not the size.
    import shutil as _shutil
    monkeypatch.setattr(_shutil, "get_terminal_size",
                        lambda *_a, **_k: _shutil.os.terminal_size((160, 48)))
    monkeypatch.setattr(sc.PlasmidApp, "run", lambda *a, **k: (_ for _ in ()).throw(
        AssertionError("must not launch over a symlinked data file")))
    with pytest.raises(SystemExit) as ei:
        sc.main()
    assert ei.value.code == 2
    err = capsys.readouterr().err
    assert "symlink" in err and str(coll) in err, err


async def test_the_daily_snapshot_predates_the_launch_restores():
    import json
    # The library holds `xOnlyHere`, the active collection does not, and no
    # dirty marker says the library is newer: the launch restore rewrites the
    # library from the collection. Today's snapshot must still hold it.
    _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
    # Written WITHOUT the mirror: `_save_library` would copy `xOnlyHere` into
    # Main too, and then there is nothing for the launch to restore over.
    sc._safe_save_json(sc._state._LIBRARY_FILE,
                       [_plain_entry("aOne", "ACGT" * 150),
                        _plain_entry("xOnlyHere", "TTGCA" * 100)],
                       "Plasmid library")
    sc._state._library_cache = None
    assert set(_disk_coll("Main")) == {"aOne"}
    assert not sc._mirror_is_dirty()
    app = sc.PlasmidApp()
    app._skip_snapshot = False
    async with app.run_test(size=TERMINAL_SIZE) as pilot:
        await _settle(pilot)
        await _wait_workers(app, pilot)
    # The scenario really happened: the launch restored over the library.
    assert {e.get("id") for e in _disk_library()} == {"aOne"}
    snaps = list((sc._state._DATA_DIR / sc._state._SNAPSHOT_DIR_NAME)
                 .glob("plasmid_library-*.json"))
    assert snaps, "no snapshot taken"
    held = {e.get("id") for p in snaps
            for e in json.loads(p.read_text()).get("entries", [])}
    assert "xOnlyHere" in held, held


def test_an_import_and_a_concurrent_save_both_land(monkeypatch):
    """The Plasmidsaurus fetch, the agent's `download-plasmidsaurus` and the
    Sanger "Add to library" each read the library, appended, and saved with
    no lock held — a save from another thread landing in between was
    overwritten by the importer's older copy (audit 2026-10-23). A writer
    that obeys the lock now waits for the import instead of being erased by
    it."""
    import threading
    import splicecraft_dataaccess as da
    _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
    real = da._load_library
    done = threading.Event()

    def other_writer():
        with sc._state._cache_lock:
            lib = real()
            lib.append(_plain_entry("fromUI", "TTGCA" * 100))
            sc._save_library(lib)
        done.set()
    t = threading.Thread(target=other_writer)

    def read_then_race():
        out = real()
        t.start()
        done.wait(0.5)          # lets an unlocked importer lose the race
        return out
    monkeypatch.setattr(da, "_load_library", read_then_race)
    sc._library_append_unique([_plain_entry("imported", "GATC" * 120),
                               _plain_entry("aOne", "CCGG" * 90)])
    t.join(10)
    assert not t.is_alive()
    sc._state._library_cache = None
    ids = {e.get("id") for e in _disk_library()}
    assert {"aOne", "aOne_2", "fromUI", "imported"} <= ids, ids
    assert set(_disk_coll("Main")) >= {"aOne", "aOne_2", "fromUI", "imported"}


class TestCollectionChangesAreOneTransaction:
    """Add / rename / delete of a collection cloned the collections, changed
    the clone and published it with no lock held — and a deferred collection
    mirror still queued for the old name then found no collection of that
    name and was dropped, its library change gone at the next launch (audit
    2026-10-23). The change now lands any queued mirror first, under one
    `_cache_lock` hold."""

    async def _queued_mirror_then(self, act):
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await _settle(pilot)
            lib = _lib(app)
            newer = [_plain_entry("aOne", "ACGT" * 150),
                     _plain_entry("queued", "TTGCA" * 100)]
            sc._safe_save_json(sc._state._LIBRARY_FILE, newer,
                               "Plasmid library")
            sc._state._library_cache = sc._typed_clone(newer)
            seq = sc._mark_mirror_pending()
            sc._state._collection_sync_pending = (
                "Main", [dict(e) for e in newer], seq)
            lib.query_one("#btn-lib-back").action_press()
            await _settle(pilot)
            await act(app, lib, pilot)
            await _wait_workers(app, pilot)
            payload = sc._state._collection_sync_pending
            if payload:                       # the background drain, next
                sc._land_deferred_mirror(*payload)

    async def test_a_rename_keeps_the_queued_mirrors_change(self):
        async def rename(app, lib, pilot):
            t = lib.query_one("#lib-coll-table", sc.DataTable)
            t.move_cursor(row=list(t.rows.keys()).index(
                next(k for k in t.rows.keys() if k.value == "Main")))
            lib._btn_coll_rename()
            await _settle(pilot)
            app.screen.dismiss("Main2")
            await _settle(pilot)
        await self._queued_mirror_then(rename)
        assert "queued" in _disk_coll("Main2"), _disk_coll("Main2")

    async def test_adding_a_collection_keeps_the_queued_mirrors_change(self):
        async def add(app, lib, pilot):
            lib._commit_collections_change(lambda colls: colls.append(
                {"name": "Other", "description": "", "plasmids": [],
                 "saved": "2026-10-24"}))
            await _settle(pilot)
        await self._queued_mirror_then(add)
        assert "queued" in _disk_coll("Main"), _disk_coll("Main")
        assert _disk_coll("Other") == {}

    def test_a_mirror_cannot_land_inside_a_collection_change(self,
                                                              monkeypatch):
        """The deferred mirror landing from its own thread between the
        change's clone and its save was reverted by the older clone — and
        the mirror had already cleared its marker, so the next launch lost
        it. Every collection change holds the lock across both now."""
        import threading
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        newer = [_plain_entry("aOne", "ACGT" * 150),
                 _plain_entry("queued", "TTGCA" * 100)]
        sc._safe_save_json(sc._state._LIBRARY_FILE, newer, "Plasmid library")
        sc._state._library_cache = sc._typed_clone(newer)
        payload = ("Main", [dict(e) for e in newer], sc._mark_mirror_pending())
        real = sc._load_collections
        done = threading.Event()
        t = threading.Thread(
            target=lambda: (sc._land_deferred_mirror(*payload), done.set()))
        raced = []

        def load_then_race():
            out = real()
            if not raced:
                raced.append(1)
                t.start()
                done.wait(0.5)
            return out
        monkeypatch.setattr(sc, "_load_collections", load_then_race)
        sc._collections_transaction(lambda colls: colls.append(
            {"name": "Other", "description": "", "plasmids": [],
             "saved": "2026-10-24"}))
        t.join(10)
        assert not t.is_alive()
        assert "queued" in _disk_coll("Main"), _disk_coll("Main")
        assert _disk_coll("Other") == {}


# ══════════════════════════════════════════════════════════════════════════
# The disk writer must persist the CACHE, not its dispatch-time payload
# ══════════════════════════════════════════════════════════════════════════

class TestAWorkerNeverOverwritesANewerEdit:
    """`_rename_save_to_disk` / `_delete_save_to_disk` run on a worker and
    were handed a parts-bin snapshot captured on the UI thread at dispatch.
    Every UI-thread parts-bin edit is an unlocked read-modify-write, so one
    landing between dispatch and the write was overwritten: gone from disk,
    still shown in the panel, and invisible until the next launch. The
    LIBRARY half of the same writer has persisted the authoritative cache
    since sweep #30; the parts bin was the half left behind (audit
    2026-10-24)."""

    def test_the_parts_bin_writer_keeps_an_edit_made_after_dispatch(self):
        dispatch = [{"id": "p1", "name": "partOne", "sequence": "ACGT"}]
        sc._save_parts_bin(dispatch)
        # The user adds a part while the worker is in flight.
        live = sc._load_parts_bin()
        live.insert(0, {"id": "p2", "name": "partTwo", "sequence": "TTTT"})
        sc._save_parts_bin(live)
        assert {e["id"] for e in sc._load_parts_bin()} == {"p1", "p2"}
        # Now the worker lands, writing what it captured at dispatch.
        with sc._cache_lock:
            sc._safe_save_json(
                sc._state._PARTS_BIN_FILE,
                sc._authoritative_parts_bin_snapshot(dispatch), "Parts bin")
        sc._state._parts_bin_cache = None
        assert {e["id"] for e in sc._load_parts_bin()} == {"p1", "p2"}

    def test_it_falls_back_to_the_payload_when_the_cache_is_unset(self):
        """Defensive: an unset cache must not write an empty bin."""
        fallback = [{"id": "p9", "name": "only", "sequence": "ACGT"}]
        sc._state._parts_bin_cache = None
        assert sc._authoritative_parts_bin_snapshot(fallback) == fallback


class TestABulkImportFlushCannotLoseAPrimer:
    """The `.dna` primer flush at the end of a bulk import is a read-modify-
    write on a WORKER thread. `_load_primers` and `_save_primers` each take
    `_cache_lock` on their own, so an unlocked flush released it in between —
    and a primer the user saved from the panel in that window was dropped
    when the flush wrote its list back. One `with _cache_lock` across both
    closes it; the lock is re-entrant, so the accessors still take it
    (audit 2026-10-24). Same shape `_library_append_unique` fixed for the
    library."""

    def test_a_concurrent_primer_save_survives_the_flush(self, tmp_path,
                                                         monkeypatch):
        import shutil
        import threading
        from pathlib import Path
        src = tmp_path / "dna_in"
        src.mkdir()
        here = Path(__file__).resolve().parent
        # DISCOVERED, not named: any `.dna` will do — the test is about the
        # import loop holding the cache lock, not about a particular file.
        # Naming one put a plasmid name from the maintainer's own library into
        # a repo that is pushed PUBLICLY (release.py's derived-name gate caught
        # it, [INV-228]), and a glob cannot leak one.
        present = sorted(here.glob("*.dna"))[:1] + \
            sorted(here.parent.glob("*.dna"))[:1]
        if not present:
            pytest.skip(".dna fixtures not present")
        for f in present:
            shutil.copy(f, src / f.name)
        sc._save_primers([])

        # A panel-style save fired from ANOTHER thread the moment the flush
        # has read the primer list. It is not joined here: the flush holds
        # the lock, so this thread blocks inside `_save_primers` until the
        # flush is done — waiting for it from inside the lock would deadlock.
        started = threading.Event()
        panel = {"id": "prPANEL", "name": "fromPanel",
                 "sequence": "AAAACCCCGGGGTTTTACGT"}

        def panel_save():
            started.set()
            cur = sc._load_primers()
            cur.insert(0, dict(panel))
            sc._save_primers(cur)

        real_load = sc._load_primers
        thread: "list[threading.Thread]" = []

        def load_then_race():
            out = real_load()
            if not thread:
                t = threading.Thread(target=panel_save)
                thread.append(t)
                t.start()
                started.wait(1.0)
            return out

        monkeypatch.setattr(sc, "_load_primers", load_then_race)
        entries, failures = sc._bulk_import_folder(src)
        monkeypatch.setattr(sc, "_load_primers", real_load)
        assert failures == [], failures
        assert entries
        if thread:
            thread[0].join(10)
            assert not thread[0].is_alive()
        sc._state._primers_cache = None
        on_disk = sc._load_primers()
        ids = {e.get("id") for e in on_disk}
        seqs = {(e.get("sequence") or "").upper() for e in on_disk}
        assert panel["sequence"] in seqs, sorted(seqs)
        assert "prPANEL" in ids, ids
        # And the import's own primers are still there.
        assert len(on_disk) > 1, on_disk


# ── The pop window the barrier could not see (audit 2026-11-09) ─────────────
# `_land_pending_collection_sync` — the barrier `_collections_transaction` and
# `_migrate_import_transaction` use to satisfy [INV-224] #3/#4 — only ever
# inspected `_state._collection_sync_pending`. The background worker POPS that
# slot and only THEN blocks on `_cache_lock`, so for essentially the whole life
# of every deferred mirror the barrier answered "nothing queued": the
# transaction went ahead, and the mirror was afterwards dropped as an orphan
# together with the `.mirror-dirty` marker that was the newer library file's
# only protection.
#
# Tested WITHOUT the live daemon, on purpose. conftest forces mirrors
# synchronous because a real deferred write bleeds across tests, and no test
# exercised the async path — which is why this survived. The worker's two
# steps are tiny and are reproduced here exactly, so the handoff is pinned
# deterministically.

def _simulate_worker_pop():
    """Exactly what `_drain_collection_sync_loop` does: take `_pending` and
    publish it as in-flight, in ONE hold, then (in the real worker) block on
    `_cache_lock`. Returns the popped payload."""
    with sc._collection_sync_lock:
        payload = sc._state._collection_sync_pending
        sc._state._collection_sync_pending = None
        sc._state._collection_sync_inflight = payload
    return payload


class TestTheBarrierSeesAPoppedMirror:
    def _queue(self):
        entries = [_plain_entry("aOne", "ACGT" * 150),
                   _plain_entry("queued", "TTGCA" * 100)]
        sc._safe_save_json(sc._state._LIBRARY_FILE, entries,
                           "Plasmid library")
        sc._state._library_cache = sc._typed_clone(entries)
        seq = sc._mark_mirror_pending()
        sc._state._collection_sync_pending = (
            "Main", [dict(e) for e in entries], seq)
        return entries, seq

    def test_the_barrier_lands_a_mirror_the_worker_already_popped(self):
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        entries, _seq = self._queue()
        payload = _simulate_worker_pop()
        assert payload is not None, "premise: a mirror was queued"
        assert sc._state._collection_sync_pending is None, \
            "premise: the worker emptied the queued slot"
        # Pre-fix this returned False — the barrier saw an empty queue and the
        # transaction proceeded as if nothing were in flight.
        assert sc._land_pending_collection_sync("test barrier") is True
        names = [p.get("name")
                 for c in sc._load_collections() if c.get("name") == "Main"
                 for p in (c.get("plasmids") or [])]
        assert "queued" in names, (
            f"the popped mirror never landed; Main holds {names}")

    def test_the_worker_declines_a_mirror_the_barrier_took_over(self):
        """The other half: once a barrier has claimed it, the worker must NOT
        write its own (now older) payload on top."""
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        _entries, seq = self._queue()
        payload = _simulate_worker_pop()
        assert sc._land_pending_collection_sync("test barrier") is True
        # The worker now reaches `_land_deferred_mirror` with `inflight=True`.
        assert sc._claim_inflight_mirror(seq) is False, (
            "the worker re-claimed a mirror the barrier had already taken")

    def test_a_claim_succeeds_when_nobody_stole_it(self):
        """Control — the claim must not refuse the normal path, or every
        deferred mirror would stop landing."""
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        _entries, seq = self._queue()
        _simulate_worker_pop()
        assert sc._claim_inflight_mirror(seq) is True
        assert sc._state._collection_sync_inflight is None

    def test_an_empty_queue_still_reports_nothing(self):
        """Control: with neither slot populated the barrier must say so, or
        every transaction would think it had work to do."""
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        with sc._collection_sync_lock:
            sc._state._collection_sync_pending = None
            sc._state._collection_sync_inflight = None
        assert sc._land_pending_collection_sync("test barrier") is False


class TestABustedCacheDropsTheMirror:
    """`_land_deferred_mirror`'s D6 rule is "write the AUTHORITATIVE library,
    never the queue-time snapshot". The `isinstance` guard simply fell through
    when `_library_cache` was None — precisely the post-
    `_bust_all_data_caches()` state of an import / restore / master-delete —
    so the stale queue-time snapshot was written and erased a newer save."""

    def test_a_none_cache_drops_rather_than_writing_the_snapshot(self):
        _seed({"Main": [_plain_entry("aOne", "ACGT" * 150)]}, "Main")
        stale = [_plain_entry("aOne", "ACGT" * 150),
                 _plain_entry("stale", "TTGCA" * 100)]
        newer = [_plain_entry("aOne", "ACGT" * 150),
                 _plain_entry("newer", "GGCCA" * 100)]
        # A newer synchronous commit lands `newer` in the collection.
        sc._save_library(newer)
        # Then the caches are busted, as an import / restore does.
        sc._state._library_cache = None
        seq = sc._mark_mirror_pending()
        sc._land_deferred_mirror("Main", [dict(e) for e in stale], seq)
        names = [p.get("name")
                 for c in sc._load_collections() if c.get("name") == "Main"
                 for p in (c.get("plasmids") or [])]
        assert "stale" not in names, (
            f"the queue-time snapshot was written over a newer save; "
            f"Main holds {names}")
        assert "newer" in names, f"the newer save was lost; Main holds {names}"
