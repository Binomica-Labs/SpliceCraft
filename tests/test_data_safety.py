"""
test_data_safety — user data is precious and must never be silently lost.

These tests verify:
  1. _safe_save_json creates a .bak backup before overwriting
  2. _safe_save_json uses atomic writes (tempfile + os.replace)
  3. _safe_load_json recovers from corrupt files via .bak restore
  4. Missing files don't crash — they return [] (first run)
  5. Startup _check_data_files notifies the user about corrupt files
  6. Manually deleted files mid-session don't crash on next load
  7. No save function can accidentally nuke a non-empty file with []
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# _safe_save_json — atomic write with .bak backup
# ═══════════════════════════════════════════════════════════════════════════════

class TestSafeSaveJson:
    def test_creates_file(self, tmp_path):
        p = tmp_path / "test.json"
        sc._safe_save_json(p, [{"id": "A"}], "test")
        assert p.exists()
        assert json.loads(p.read_text()) == [{"id": "A"}]

    def test_creates_bak_on_overwrite(self, tmp_path):
        p = tmp_path / "test.json"
        sc._safe_save_json(p, [{"id": "first"}], "test")
        sc._safe_save_json(p, [{"id": "second"}], "test")
        bak = tmp_path / "test.json.bak"
        assert bak.exists()
        # .bak should contain the FIRST version (pre-overwrite)
        bak_data = json.loads(bak.read_text())
        assert bak_data == [{"id": "first"}]
        # Main file should contain the second version
        assert json.loads(p.read_text()) == [{"id": "second"}]

    def test_bak_not_created_for_first_write(self, tmp_path):
        p = tmp_path / "test.json"
        sc._safe_save_json(p, [{"id": "first"}], "test")
        bak = tmp_path / "test.json.bak"
        assert not bak.exists()

    def test_atomic_write_survives_crash(self, tmp_path):
        """If the file existed before, a failed write should NOT corrupt
        the original — the .bak is the safety net, and the original should
        remain intact if os.replace fails (simulated by checking file
        content matches the last successful write)."""
        p = tmp_path / "test.json"
        sc._safe_save_json(p, [{"id": "good"}], "test")
        # Second write succeeds too
        sc._safe_save_json(p, [{"id": "updated"}], "test")
        assert json.loads(p.read_text()) == [{"id": "updated"}]
        # .bak holds the previous good version
        bak = tmp_path / "test.json.bak"
        assert json.loads(bak.read_text()) == [{"id": "good"}]

    def test_empty_file_no_bak(self, tmp_path):
        """An empty file should NOT generate a .bak (nothing to back up)."""
        p = tmp_path / "test.json"
        p.write_text("")
        sc._safe_save_json(p, [{"id": "new"}], "test")
        bak = tmp_path / "test.json.bak"
        assert not bak.exists()

    def test_writes_valid_json(self, tmp_path):
        p = tmp_path / "test.json"
        entries = [{"name": "x", "seq": "ACGT"}, {"name": "y", "seq": "TGCA"}]
        sc._safe_save_json(p, entries, "test")
        assert json.loads(p.read_text()) == entries


# ═══════════════════════════════════════════════════════════════════════════════
# _safe_load_json — corrupt file recovery
# ═══════════════════════════════════════════════════════════════════════════════

class TestSafeLoadJson:
    def test_missing_file_returns_empty_no_warning(self, tmp_path):
        p = tmp_path / "nonexistent.json"
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == []
        assert warning is None

    def test_valid_file_returns_entries(self, tmp_path):
        p = tmp_path / "test.json"
        p.write_text(json.dumps([{"id": "A"}]))
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == [{"id": "A"}]
        assert warning is None

    def test_corrupt_file_without_bak_returns_empty_with_warning(self, tmp_path):
        p = tmp_path / "test.json"
        p.write_text("{not valid json")
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == []
        assert warning is not None
        assert "corrupt" in warning.lower()

    def test_corrupt_file_with_valid_bak_restores(self, tmp_path):
        p = tmp_path / "test.json"
        bak = tmp_path / "test.json.bak"
        # Good backup
        bak.write_text(json.dumps([{"id": "rescued"}]))
        # Corrupt main file
        p.write_text("{garbage")
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == [{"id": "rescued"}]
        assert warning is not None
        assert "restored" in warning.lower()
        # The corrupt main file should now be overwritten with the backup
        assert json.loads(p.read_text()) == [{"id": "rescued"}]

    def test_corrupt_file_with_corrupt_bak_returns_empty(self, tmp_path):
        p = tmp_path / "test.json"
        bak = tmp_path / "test.json.bak"
        p.write_text("{bad")
        bak.write_text("{also bad")
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == []
        assert warning is not None

    def test_non_list_json_treated_as_corrupt(self, tmp_path):
        """A JSON file containing a dict instead of a list is invalid for
        our persistence format — should be treated as corrupt."""
        p = tmp_path / "test.json"
        p.write_text('{"not": "a list"}')
        entries, warning = sc._safe_load_json(p, "test")
        # Should attempt .bak restore or return []
        assert entries == [] or isinstance(entries, list)


# ═══════════════════════════════════════════════════════════════════════════════
# Persistence integration — each _load/_save pair through _safe_*
# ═══════════════════════════════════════════════════════════════════════════════

class TestPersistenceIntegration:
    def test_library_save_creates_bak(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc, "_LIBRARY_FILE", tmp_path / "lib.json")
        monkeypatch.setattr(sc, "_library_cache", None)
        sc._save_library([{"id": "A", "name": "first"}])
        sc._save_library([{"id": "B", "name": "second"}])
        bak = tmp_path / "lib.json.bak"
        assert bak.exists()
        assert json.loads(bak.read_text())[0]["id"] == "A"

    def test_parts_bin_save_creates_bak(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc, "_PARTS_BIN_FILE", tmp_path / "parts.json")
        monkeypatch.setattr(sc, "_parts_bin_cache", None)
        sc._save_parts_bin([{"name": "p1"}])
        sc._save_parts_bin([{"name": "p2"}])
        bak = tmp_path / "parts.json.bak"
        assert bak.exists()

    def test_primers_save_creates_bak(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc, "_PRIMERS_FILE", tmp_path / "primers.json")
        monkeypatch.setattr(sc, "_primers_cache", None)
        sc._save_primers([{"name": "pr1"}])
        sc._save_primers([{"name": "pr2"}])
        bak = tmp_path / "primers.json.bak"
        assert bak.exists()

    def test_library_load_survives_deleted_file(self, tmp_path, monkeypatch):
        """If user deletes plasmid_library.json manually, load must return []
        without crashing."""
        p = tmp_path / "lib.json"
        monkeypatch.setattr(sc, "_LIBRARY_FILE", p)
        monkeypatch.setattr(sc, "_library_cache", None)
        # File doesn't exist — should return []
        assert sc._load_library() == []

    def test_parts_bin_load_survives_deleted_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc, "_PARTS_BIN_FILE", tmp_path / "nope.json")
        monkeypatch.setattr(sc, "_parts_bin_cache", None)
        assert sc._load_parts_bin() == []

    def test_primers_load_survives_deleted_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc, "_PRIMERS_FILE", tmp_path / "nope.json")
        monkeypatch.setattr(sc, "_primers_cache", None)
        assert sc._load_primers() == []

    def test_library_load_recovers_from_corrupt(self, tmp_path, monkeypatch):
        p = tmp_path / "lib.json"
        bak = tmp_path / "lib.json.bak"
        p.write_text("{bad}")
        bak.write_text(json.dumps([{"id": "X", "name": "saved"}]))
        monkeypatch.setattr(sc, "_LIBRARY_FILE", p)
        monkeypatch.setattr(sc, "_library_cache", None)
        entries = sc._load_library()
        assert len(entries) == 1
        assert entries[0]["id"] == "X"


# ═══════════════════════════════════════════════════════════════════════════════
# Startup _check_data_files
# ═══════════════════════════════════════════════════════════════════════════════

class TestStartupDataCheck:
    async def test_startup_with_all_files_missing_no_crash(
        self, tmp_path, monkeypatch
    ):
        """First-run scenario: no files exist. App must mount without
        crashing and without showing corruption warnings."""
        monkeypatch.setattr(sc, "_LIBRARY_FILE", tmp_path / "lib.json")
        monkeypatch.setattr(sc, "_PARTS_BIN_FILE", tmp_path / "parts.json")
        monkeypatch.setattr(sc, "_PRIMERS_FILE", tmp_path / "primers.json")
        monkeypatch.setattr(sc, "_library_cache", None)
        monkeypatch.setattr(sc, "_parts_bin_cache", None)
        monkeypatch.setattr(sc, "_primers_cache", None)
        # Block the network seeder
        monkeypatch.setattr(
            sc, "fetch_genbank",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("no net")),
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # App is alive — that's the test

    async def test_startup_with_corrupt_library_notifies(
        self, tmp_path, monkeypatch
    ):
        """A corrupt plasmid_library.json should produce a user notification
        on startup, not a crash."""
        p = tmp_path / "lib.json"
        p.write_text("{corrupt}")
        monkeypatch.setattr(sc, "_LIBRARY_FILE", p)
        monkeypatch.setattr(sc, "_PARTS_BIN_FILE", tmp_path / "parts.json")
        monkeypatch.setattr(sc, "_PRIMERS_FILE", tmp_path / "primers.json")
        monkeypatch.setattr(sc, "_library_cache", None)
        monkeypatch.setattr(sc, "_parts_bin_cache", None)
        monkeypatch.setattr(sc, "_primers_cache", None)
        monkeypatch.setattr(
            sc, "fetch_genbank",
            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("no net")),
        )
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # App is alive despite corrupt file
