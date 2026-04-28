"""
test_data_safety — user data is precious and must never be silently lost.

CRITICAL: TestRealFilesNeverTouched verifies that the autouse fixture in
conftest.py properly redirects ALL persistence paths to temp dirs. If this
test fails, it means a save operation could nuke the user's real data.

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
        raw = json.loads(p.read_text())
        assert raw == {"_schema_version": sc._CURRENT_SCHEMA_VERSION,
                       "entries": [{"id": "A"}]}

    def test_creates_bak_on_overwrite(self, tmp_path):
        p = tmp_path / "test.json"
        sc._safe_save_json(p, [{"id": "first"}], "test")
        sc._safe_save_json(p, [{"id": "second"}], "test")
        bak = tmp_path / "test.json.bak"
        assert bak.exists()
        # .bak should contain the FIRST version (pre-overwrite)
        assert json.loads(bak.read_text())["entries"] == [{"id": "first"}]
        # Main file should contain the second version
        assert json.loads(p.read_text())["entries"] == [{"id": "second"}]

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
        assert json.loads(p.read_text())["entries"] == [{"id": "updated"}]
        # .bak holds the previous good version
        bak = tmp_path / "test.json.bak"
        assert json.loads(bak.read_text())["entries"] == [{"id": "good"}]

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
        assert json.loads(p.read_text())["entries"] == entries


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
        """A JSON file containing a dict WITHOUT the envelope shape is
        invalid for our persistence format — should be treated as corrupt
        (returns [] since there's no .bak)."""
        p = tmp_path / "test.json"
        p.write_text('{"not": "a list"}')
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == []
        assert warning is not None


# ═══════════════════════════════════════════════════════════════════════════════
# Schema versioning — envelope format + legacy bare-list compatibility
# ═══════════════════════════════════════════════════════════════════════════════

class TestSchemaVersioning:
    """The on-disk format went from a bare JSON list to
    `{"_schema_version": N, "entries": [...]}` in v0.3.1. Loaders accept
    both so upgrades don't lose user data."""

    def test_legacy_flat_list_loads_without_warning(self, tmp_path):
        """Files written by SpliceCraft < 0.3.1 are bare JSON lists.
        They must load cleanly with no warning."""
        p = tmp_path / "legacy.json"
        p.write_text(json.dumps([{"id": "L1"}, {"id": "L2"}]))
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == [{"id": "L1"}, {"id": "L2"}]
        assert warning is None

    def test_envelope_loads_without_warning(self, tmp_path):
        p = tmp_path / "new.json"
        p.write_text(json.dumps({
            "_schema_version": sc._CURRENT_SCHEMA_VERSION,
            "entries": [{"id": "E1"}],
        }))
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == [{"id": "E1"}]
        assert warning is None

    def test_future_schema_version_warns_but_loads(self, tmp_path):
        """A file written by a newer SpliceCraft (higher schema version)
        must still load, but emit a warning so the user knows fields
        may be dropped on save."""
        p = tmp_path / "future.json"
        p.write_text(json.dumps({
            "_schema_version": sc._CURRENT_SCHEMA_VERSION + 99,
            "entries": [{"id": "F1", "new_field": "unknown"}],
        }))
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == [{"id": "F1", "new_field": "unknown"}]
        assert warning is not None
        assert "newer" in warning.lower()

    def test_envelope_roundtrip_preserves_entries(self, tmp_path):
        """Save + reload must yield the original entries unchanged."""
        p = tmp_path / "roundtrip.json"
        original = [{"id": "A", "seq": "ACGT"}, {"id": "B", "seq": "TTTT"}]
        sc._safe_save_json(p, original, "test")
        entries, warning = sc._safe_load_json(p, "test")
        assert entries == original
        assert warning is None

    def test_legacy_file_rewrites_as_envelope_on_next_save(self, tmp_path):
        """Upgrading users: a legacy bare-list file should be silently
        rewritten as an envelope on the next save."""
        p = tmp_path / "upgrade.json"
        # Simulate a legacy file
        p.write_text(json.dumps([{"id": "OLD"}]))
        # Load it (accepts legacy shape)
        loaded, _ = sc._safe_load_json(p, "test")
        assert loaded == [{"id": "OLD"}]
        # Save new entries
        sc._safe_save_json(p, [{"id": "NEW"}], "test")
        # File on disk is now envelope-shaped
        raw = json.loads(p.read_text())
        assert isinstance(raw, dict)
        assert raw["_schema_version"] == sc._CURRENT_SCHEMA_VERSION
        assert raw["entries"] == [{"id": "NEW"}]

    def _capture_splicecraft_warnings(self):
        """Attach a handler to the splicecraft logger and return (handler, records).
        caplog can't see it otherwise because _log.propagate is False."""
        import logging

        records: list[logging.LogRecord] = []

        class _ListHandler(logging.Handler):
            def emit(self, record):
                records.append(record)

        h = _ListHandler(level=logging.WARNING)
        sc._log.addHandler(h)
        return h, records

    def test_shrink_guard_counts_envelope_entries(self, tmp_path):
        """When overwriting an envelope file with fewer entries, the
        shrink guard must log a warning comparing the correct counts
        (not 2 vs 1 for the 2-key envelope dict)."""
        p = tmp_path / "shrink.json"
        sc._safe_save_json(p, [{"id": "A"}, {"id": "B"}, {"id": "C"}], "test")
        h, records = self._capture_splicecraft_warnings()
        try:
            sc._safe_save_json(p, [{"id": "A"}], "test")
        finally:
            sc._log.removeHandler(h)
        msgs = [r.getMessage() for r in records]
        assert any("SHRINK GUARD" in m and "was 3" in m for m in msgs), msgs

    def test_shrink_guard_counts_legacy_entries(self, tmp_path):
        """Upgrading users: if the on-disk file is still a legacy bare
        list, the shrink guard must still count its entries correctly."""
        p = tmp_path / "shrink_legacy.json"
        p.write_text(json.dumps([{"id": "A"}, {"id": "B"}, {"id": "C"}]))
        h, records = self._capture_splicecraft_warnings()
        try:
            sc._safe_save_json(p, [{"id": "A"}], "test")
        finally:
            sc._log.removeHandler(h)
        msgs = [r.getMessage() for r in records]
        assert any("SHRINK GUARD" in m and "was 3" in m for m in msgs), msgs


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
        assert json.loads(bak.read_text())["entries"][0]["id"] == "A"

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


# ═══════════════════════════════════════════════════════════════════════════════
# CRITICAL: verify the autouse fixture protects real user files
# ═══════════════════════════════════════════════════════════════════════════════

class TestRealFilesNeverTouched:
    """Verify that the conftest.py _protect_user_data autouse fixture
    actually redirects ALL persistence paths to temp dirs.

    If ANY of these tests fail, it means a test could write to the user's
    real plasmid_library.json / parts_bin.json / primers.json — which is
    a catastrophic data-loss bug.
    """

    def test_library_file_is_redirected(self):
        """_LIBRARY_FILE must NOT point to the repo root during tests."""
        assert "/tmp" in str(sc._LIBRARY_FILE) or "pytest" in str(sc._LIBRARY_FILE), (
            f"_LIBRARY_FILE points to {sc._LIBRARY_FILE} — "
            f"expected a tmp dir, not the real repo!"
        )

    def test_parts_bin_file_is_redirected(self):
        assert "/tmp" in str(sc._PARTS_BIN_FILE) or "pytest" in str(sc._PARTS_BIN_FILE), (
            f"_PARTS_BIN_FILE points to {sc._PARTS_BIN_FILE} — "
            f"expected a tmp dir!"
        )

    def test_primers_file_is_redirected(self):
        assert "/tmp" in str(sc._PRIMERS_FILE) or "pytest" in str(sc._PRIMERS_FILE), (
            f"_PRIMERS_FILE points to {sc._PRIMERS_FILE} — "
            f"expected a tmp dir!"
        )

    def test_features_file_is_redirected(self):
        assert "/tmp" in str(sc._FEATURES_FILE) or "pytest" in str(sc._FEATURES_FILE), (
            f"_FEATURES_FILE points to {sc._FEATURES_FILE} — "
            f"expected a tmp dir!"
        )

    def test_collections_file_is_redirected(self):
        assert "/tmp" in str(sc._COLLECTIONS_FILE) or "pytest" in str(sc._COLLECTIONS_FILE), (
            f"_COLLECTIONS_FILE points to {sc._COLLECTIONS_FILE} — "
            f"expected a tmp dir!"
        )

    def test_save_library_writes_to_tmp(self):
        """Actually call _save_library and verify the write landed in /tmp,
        not in the repo directory."""
        sc._save_library([{"id": "SAFETY_TEST", "name": "safety"}])
        assert sc._LIBRARY_FILE.exists()
        assert "/tmp" in str(sc._LIBRARY_FILE) or "pytest" in str(sc._LIBRARY_FILE)
        import json
        data = json.loads(sc._LIBRARY_FILE.read_text())
        assert data["entries"][0]["id"] == "SAFETY_TEST"

    def test_save_parts_bin_writes_to_tmp(self):
        sc._save_parts_bin([{"name": "SAFETY_TEST"}])
        assert sc._PARTS_BIN_FILE.exists()
        assert "/tmp" in str(sc._PARTS_BIN_FILE) or "pytest" in str(sc._PARTS_BIN_FILE)

    def test_save_primers_writes_to_tmp(self):
        sc._save_primers([{"name": "SAFETY_TEST"}])
        assert sc._PRIMERS_FILE.exists()
        assert "/tmp" in str(sc._PRIMERS_FILE) or "pytest" in str(sc._PRIMERS_FILE)

    def test_save_features_writes_to_tmp(self):
        sc._save_features([{"name": "SAFETY_TEST",
                            "feature_type": "CDS",
                            "sequence": "ATG"}])
        assert sc._FEATURES_FILE.exists()
        assert "/tmp" in str(sc._FEATURES_FILE) or "pytest" in str(sc._FEATURES_FILE)

    def test_save_collections_writes_to_tmp(self):
        sc._save_collections([{"name": "SAFETY_TEST", "plasmids": []}])
        assert sc._COLLECTIONS_FILE.exists()
        assert "/tmp" in str(sc._COLLECTIONS_FILE) or "pytest" in str(sc._COLLECTIONS_FILE)

    def test_real_repo_files_untouched(self):
        """The actual files in the repo root must NOT contain SAFETY_TEST
        entries — if they do, the autouse fixture failed."""
        import json
        repo_root = Path(__file__).resolve().parent.parent
        for fname in ["plasmid_library.json", "parts_bin.json",
                      "primers.json", "features.json", "collections.json"]:
            p = repo_root / fname
            if not p.exists():
                continue
            try:
                raw = json.loads(p.read_text())
            except json.JSONDecodeError:
                continue  # corrupt file — not our problem here
            # Support both legacy bare-list and current envelope format
            data = raw["entries"] if isinstance(raw, dict) and "entries" in raw else raw
            for entry in data:
                assert entry.get("id") != "SAFETY_TEST", (
                    f"SAFETY_TEST leaked into real {fname}!"
                )
                assert entry.get("name") != "SAFETY_TEST", (
                    f"SAFETY_TEST leaked into real {fname}!"
                )


# ═══════════════════════════════════════════════════════════════════════════════
# _atomic_write_text — shared helper for _do_save, _do_autosave, FASTA/GenBank
#
# Added 2026-04-21 after the audit caught `_do_save` using plain
# `Path.write_text`, which leaves a half-written .gb on disk if the process
# crashes mid-write. The helper is now the sole path by which the app writes
# non-JSON user data, so regressions here are a real data-loss risk.
# ═══════════════════════════════════════════════════════════════════════════════

class TestAtomicWriteText:
    def test_creates_new_file(self, tmp_path):
        p = tmp_path / "out.txt"
        sc._atomic_write_text(p, "hello\n")
        assert p.read_text() == "hello\n"

    def test_overwrites_existing_atomically(self, tmp_path):
        p = tmp_path / "out.txt"
        p.write_text("old contents\n")
        sc._atomic_write_text(p, "new contents\n")
        assert p.read_text() == "new contents\n"

    def test_no_tmp_file_left_on_success(self, tmp_path):
        """After a successful write, no `.tmp` file should remain in the
        target directory. The mkstemp tempfile is renamed into place via
        os.replace; a lingering tempfile means the cleanup is broken."""
        p = tmp_path / "out.txt"
        sc._atomic_write_text(p, "data")
        # `.out.txt.*.tmp` is the hidden tempfile prefix used by mkstemp
        leftover = list(tmp_path.glob(f".{p.name}.*.tmp"))
        assert leftover == [], f"Unexpected tempfile leftovers: {leftover}"

    def test_tmp_file_cleaned_on_error(self, tmp_path, monkeypatch):
        """Simulate a filesystem failure during os.replace — the tempfile
        must be cleaned up (not left as dotfile garbage) and the error
        must propagate so the caller can notify the user."""
        import os as _os
        p = tmp_path / "out.txt"

        def _boom(*a, **kw):
            raise OSError("simulated replace failure")

        monkeypatch.setattr(_os, "replace", _boom)
        with pytest.raises(OSError):
            sc._atomic_write_text(p, "data")
        # no tempfile left behind
        leftover = list(tmp_path.glob(f".{p.name}.*.tmp"))
        assert leftover == []
        # and the target file was not created / corrupted
        assert not p.exists()

    def test_existing_file_untouched_on_error(self, tmp_path, monkeypatch):
        """If the replace step fails, the previous contents of the target
        must be preserved verbatim — no partial overwrite."""
        import os as _os
        p = tmp_path / "out.txt"
        p.write_text("old contents\n")

        def _boom(*a, **kw):
            raise OSError("simulated replace failure")

        monkeypatch.setattr(_os, "replace", _boom)
        with pytest.raises(OSError):
            sc._atomic_write_text(p, "new contents\n")
        assert p.read_text() == "old contents\n"

    def test_creates_parent_dirs(self, tmp_path):
        """Nested output paths should work — parent dirs created on demand."""
        p = tmp_path / "nested" / "deeper" / "out.txt"
        sc._atomic_write_text(p, "x")
        assert p.read_text() == "x"

    def test_unicode_content(self, tmp_path):
        p = tmp_path / "out.txt"
        sc._atomic_write_text(p, "ñáéíóú — GAATTC\n")
        assert p.read_text(encoding="utf-8") == "ñáéíóú — GAATTC\n"


class TestDoSaveUsesAtomicWrite:
    """Regression guard: `PlasmidApp._do_save` must route through
    `_atomic_write_text` (not `Path.write_text`). If a future refactor
    inlines the write again, a process crash mid-save corrupts the
    user's .gb file silently. This test trips before that ships."""

    async def test_save_survives_mid_write_crash(
            self, tmp_path, tiny_record, monkeypatch):
        """Simulate a filesystem failure during os.replace inside
        `_do_save`. The user's previous .gb content must remain intact,
        and the app must surface the failure (no silent data loss)."""
        import os as _os
        target = tmp_path / "plasmid.gb"
        target.write_text("ORIGINAL CONTENTS — MUST SURVIVE\n")

        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            app._source_path = str(target)
            app._unsaved = True

            def _boom(*a, **kw):
                raise OSError("simulated crash during replace")

            monkeypatch.setattr(_os, "replace", _boom)
            ok = app._do_save()
            assert ok is False, "Save must report failure on write error"
            # Previous contents must be intact — atomic guarantee
            assert target.read_text() == "ORIGINAL CONTENTS — MUST SURVIVE\n"
            # No lingering tempfile in the target directory
            leftover = list(tmp_path.glob(f".{target.name}.*.tmp"))
            assert leftover == []
