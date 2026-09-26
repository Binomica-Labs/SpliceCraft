"""test_sweep9 — regression coverage for adversarial audit sweep #9 (2026-05-19).

Six parallel audit agents (concurrency / data-integrity / security /
biology / UI / robustness) examined the 0.9.4–0.9.6 surface and the
sacred-invariant catalog. This file locks in the HIGH + MEDIUM fixes
shipped under the sweep so a future refactor cannot quietly regress
them.

The `_protect_user_data` autouse fixture in `conftest.py` redirects every
`_*_FILE` constant to a temp dir, so module-global writes stay isolated
from the developer's real data.
"""
from __future__ import annotations

import io
import zipfile
from pathlib import Path

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# H1 — `_h_restore_pre_update_snapshot` busts every 0.9.6+ cache
# ═══════════════════════════════════════════════════════════════════════════════

class TestRestoreEndpointBustsNewCaches:
    """Pre-fix the cache-bust list under `_h_restore_pre_update_snapshot`
    enumerated 11 pre-0.9.6 caches but missed `_experiments_cache`,
    `_experiment_projects_cache`, and `_gels_cache`. A restore would
    write the old files back to disk, but the next UI mutation would
    silently overwrite them from the stale in-memory cache."""

    def test_handler_source_lists_all_three_new_caches(self):
        # Sweep #25 (2026-05-23): the handler no longer hand-lists
        # cache names — it iterates `_MASTER_DELETE_CACHE_ATTRS` (the
        # canonical source of truth). Verify both that the handler
        # iterates the master tuple AND that the master tuple
        # actually contains the three sweep-#9 caches.
        import inspect
        src = inspect.getsource(sc._h_restore_pre_update_snapshot)
        assert "_MASTER_DELETE_CACHE_ATTRS" in src, (
            "Restore handler must iterate the canonical "
            "_MASTER_DELETE_CACHE_ATTRS tuple, not a hand-list"
        )
        for name in (
            "_experiments_cache",
            "_experiment_projects_cache",
            "_gels_cache",
        ):
            assert name in sc._MASTER_DELETE_CACHE_ATTRS, (
                f"_MASTER_DELETE_CACHE_ATTRS missing {name}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# H2 — `RestoreFromBackupModal._TARGETS` includes new persisted files
# ═══════════════════════════════════════════════════════════════════════════════

class TestRestoreModalCoversNewFiles:
    """Without the new files in `_TARGETS`, a user with a backup on
    disk had no in-app path to recover from corruption. The agent-API
    cache-bust (H1) is necessary but not sufficient — the UI also
    needs to be able to PICK the target."""

    def test_targets_include_experiments_projects_gels(self):
        labels = [label for label, _attr in sc.RestoreFromBackupModal._TARGETS]
        assert "Experiments" in labels
        assert "Experiment projects" in labels
        assert "Gels" in labels

    def test_target_attrs_resolve_to_actual_constants(self):
        """The second tuple element is the attribute name on the
        module — confirm each resolves to a real Path constant."""
        for _label, attr in sc.RestoreFromBackupModal._TARGETS:
            val = sc._resolve_state_or_hub(attr)
            assert val is not None, f"{attr} missing from hub or _state"
            assert isinstance(val, Path), f"{attr} not a Path"


# ═══════════════════════════════════════════════════════════════════════════════
# H3 — Keyboard Ctrl+Z routes through `action_undo` (respects modal block)
# ═══════════════════════════════════════════════════════════════════════════════

class TestCtrlZRespectsBlocksUndo:
    """Pre-fix `App.on_key` called `self._action_undo()` (private)
    directly for Ctrl+Z, bypassing `_undo_blocked_by_modal()` — so
    modal `_blocks_undo` was dead code via the keyboard. Menu-driven
    Undo did respect the guard. Now keyboard also routes through the
    public `action_undo` path."""

    def test_on_key_source_uses_public_action_undo(self):
        import inspect
        # We don't have an easy unit-mockable surface for the App-
        # level `on_key` since it depends on a mounted Textual app.
        # Inspect the source to confirm the routing fix is in place.
        src = inspect.getsource(sc.PlasmidApp.on_key)
        # The Ctrl+Z branch must call the PUBLIC `action_undo`
        # (NOT the private `_action_undo`) so the modal-block guard
        # fires.
        # Cheap pattern: count both, public should appear (at
        # least) once in the ctrl+z and ctrl+shift+z branches.
        assert "self.action_undo()" in src
        assert "self.action_redo()" in src
        # The private-method bypass (`self._action_undo()`) MUST
        # NOT appear in `on_key`.
        assert "self._action_undo()" not in src
        assert "self._action_redo()" not in src


# ═══════════════════════════════════════════════════════════════════════════════
# H8 — Windows zip separator resolver
# ═══════════════════════════════════════════════════════════════════════════════

class TestWindowsZipSeparatorResolver:
    """A Plasmidsaurus zip built on Windows can ship member names like
    `cat\\file.gbk` in the central directory. Pre-fix our code stored
    the normalized (forward-slash) form and called
    `zf.getinfo(stored_name)` later, which KeyError'd because zipfile
    matches on exact stored name. `_zf_get_member_info` retries with
    the normalized form folded back to either side."""

    def _make_backslash_zip(self, tmp_path: Path) -> Path:
        """Build a zip whose central directory uses backslash
        separators (the Windows ecosystem case)."""
        p = tmp_path / "win_zip.zip"
        with zipfile.ZipFile(p, "w") as zf:
            # `zipfile.writestr` accepts a name with arbitrary chars —
            # passing `cat\\file.gbk` writes that exact string.
            zf.writestr("category\\sample.gbk", b"LOCUS test 100 bp\n//\n")
        return p

    def test_resolver_finds_backslash_member_via_normalized_lookup(
        self, tmp_path,
    ):
        zip_path = self._make_backslash_zip(tmp_path)
        with zipfile.ZipFile(zip_path, "r") as zf:
            # Direct getinfo with forward slashes would fail; our
            # resolver folds and finds it.
            info = sc._zf_get_member_info(zf, "category/sample.gbk")
            assert info is not None
            assert info.file_size > 0

    def test_resolver_passes_through_already_correct_name(
        self, tmp_path,
    ):
        # The common case: name matches exactly. No fold needed.
        p = tmp_path / "linux_zip.zip"
        with zipfile.ZipFile(p, "w") as zf:
            zf.writestr("category/sample.gbk", b"x")
        with zipfile.ZipFile(p, "r") as zf:
            info = sc._zf_get_member_info(zf, "category/sample.gbk")
            assert info.file_size == 1

    def test_resolver_raises_keyerror_when_missing(self, tmp_path):
        p = tmp_path / "z.zip"
        with zipfile.ZipFile(p, "w") as zf:
            zf.writestr("real.txt", b"hi")
        with zipfile.ZipFile(p, "r") as zf:
            with pytest.raises(KeyError):
                sc._zf_get_member_info(zf, "missing.txt")


# ═══════════════════════════════════════════════════════════════════════════════
# H9 — Tag regexes reject HTML entities + URL params
# ═══════════════════════════════════════════════════════════════════════════════

class TestTagRegexRejectsHTMLEntities:
    """Pre-fix `_GEL_REF_RE` matched the entity name inside any pasted
    HTML / markdown export (`&amp;`, `&nbsp;`, `&copy;`...), polluting
    `attached_gel_ids` on save, false-highlighting in the editor, and
    surfacing a misleading notify on Ctrl+G click-through. The atomic-
    group + trailing-reject pattern (sweep #9) rejects entire entity
    sequences instead of matching shorter prefixes."""

    @pytest.mark.parametrize("haystack", [
        "&amp;", "&nbsp;", "&copy;", "&lt;X&gt;", "foo &amp; bar",
        "&author=jane",  # URL param
    ])
    def test_gel_regex_does_not_match_entity_or_param(self, haystack):
        matches = list(sc._GEL_REF_RE.finditer(haystack))
        assert matches == [], f"unexpected match in {haystack!r}: {matches}"

    @pytest.mark.parametrize("haystack,expected", [
        ("&gel-abc123", ["&gel-abc123"]),
        ("&gel-abc &gel-def", ["&gel-abc", "&gel-def"]),
        ("see &gel-1, please", ["&gel-1"]),
        ("end with &gel", ["&gel"]),
    ])
    def test_gel_regex_still_matches_legitimate_tags(
        self, haystack, expected,
    ):
        matches = [m.group(0) for m in sc._GEL_REF_RE.finditer(haystack)]
        assert matches == expected

    def test_plasmid_regex_rejects_url_param_form(self):
        assert list(sc._PLASMID_REF_RE.finditer("@author=jane")) == []

    def test_plasmid_regex_still_blocks_emails(self):
        # Existing pre-sweep behaviour — confirm we didn't regress
        # it while reshaping the pattern.
        assert list(sc._PLASMID_REF_RE.finditer("user@example.com")) == []

    def test_action_regex_reads_a_trailing_semicolon_as_punctuation(self):
        # Only `&` forms an HTML entity; after `!` a semicolon is prose
        # ("!digest; then !ligate"), so the audit 2026-09-22 reversed this.
        assert [m.group(1) for m in sc._ACTIONS_REF_RE.finditer(
            "!digest; then !ligate")] == ["digest", "ligate"]

    def test_action_regex_still_matches_normal_tag(self):
        matches = [m.group(0) for m in sc._ACTIONS_REF_RE.finditer("!digest")]
        assert matches == ["!digest"]

    def test_extract_gel_refs_does_not_pollute_with_html(self):
        body = "Used &amp; and &nbsp; entities, then &gel-real-1 in prose."
        assert sc._extract_gel_refs(body) == ["gel-real-1"]


# ═══════════════════════════════════════════════════════════════════════════════
# H10 — Lockfile stale detection via argv (PID-recycle resistance)
# ═══════════════════════════════════════════════════════════════════════════════

class TestPidIsSplicecraft:
    """Pre-fix the staleness check was `_pid_alive` only — a long-uptime
    system whose PID counter wraps and re-assigns the recorded PID to
    an unrelated process (sshd, bash, etc.) locked the user out of
    their own data dir indefinitely. `_pid_is_splicecraft` inspects
    `/proc/<pid>/cmdline` so PID-recycled-to-unrelated is now detected
    as stale."""

    def test_returns_false_for_pytest_process(self):
        # Our own argv contains "pytest", NOT "splicecraft" — so the
        # function should report False (it's a NEGATIVE control:
        # pytest is not a splicecraft instance even though our test
        # imports the splicecraft module). On platforms where the
        # check isn't implementable, None is returned.
        import os as _os
        result = sc._pid_is_splicecraft(_os.getpid())
        if result is None:
            pytest.skip("/proc/<pid>/cmdline not available on this platform")
        assert result is False

    def test_returns_false_for_invalid_pid(self):
        # PID 0 / negative should short-circuit to False without
        # touching the FS.
        assert sc._pid_is_splicecraft(0) is False
        assert sc._pid_is_splicecraft(-1) is False

    def test_returns_false_for_nonexistent_pid(self):
        # PID very-high-number unlikely to be live.
        result = sc._pid_is_splicecraft(2_000_000)
        # On platforms without /proc, we return None (preserves the
        # pessimistic "assume live" stance).
        assert result in (False, None)

    def test_detects_splicecraft_via_mocked_cmdline(
        self, tmp_path, monkeypatch,
    ):
        # Stub `open()` for `/proc/<pid>/cmdline` so we can exercise
        # the substring-detection logic without depending on a real
        # live splicecraft process.
        import builtins
        real_open = builtins.open
        fake_pid = 999_999

        def _stub_open(path, *args, **kwargs):
            if str(path) == f"/proc/{fake_pid}/cmdline":
                return io.BytesIO(b"/usr/bin/python3\x00splicecraft\x00")
            return real_open(path, *args, **kwargs)

        monkeypatch.setattr(builtins, "open", _stub_open)
        result = sc._pid_is_splicecraft(fake_pid)
        if result is None:
            pytest.skip("/proc/<pid>/cmdline path not used on this platform")
        assert result is True

    def test_path_component_is_not_a_running_instance(self, monkeypatch):
        """Regression (2026-06-08): "splicecraft" appearing only as a
        PATH COMPONENT — the interpreter path, a venv dir, the CWD, or
        an editor's target dir — must NOT register as a live instance;
        only the actual command/script token counts. Otherwise the pipx
        venv path (.../pipx/venvs/splicecraft/...) or running the suite
        from a dev tree at ~/SpliceCraft/.venv/ reads as "running" and
        causes a spurious "another instance is already running" lock-out
        when the recorded PID is recycled."""
        import builtins, io
        real_open = builtins.open
        fake_pid = 999_001
        # The pytest-from-an-in-tree-venv shape: splicecraft is only in
        # the interpreter PATH; there is no splicecraft *command* token.
        cmdline = (b"/home/u/SpliceCraft/.venv/bin/python\x00-m\x00pytest"
                   b"\x00tests/test_x.py\x00")

        def _stub_open(path, *a, **k):
            if str(path) == f"/proc/{fake_pid}/cmdline":
                return io.BytesIO(cmdline)
            return real_open(path, *a, **k)

        monkeypatch.setattr(builtins, "open", _stub_open)
        result = sc._pid_is_splicecraft(fake_pid)
        if result is None:
            pytest.skip("/proc/<pid>/cmdline path not used on this platform")
        assert result is False

    def test_detects_every_real_launch_mode(self, monkeypatch):
        """The tightened matcher still detects EVERY real launch mode, so
        a genuinely live instance is never mistaken for dead (which would
        free the data-dir lock and let two instances clobber the shared
        library): console script, `python splicecraft.py`,
        `python -m splicecraft`, and the pipx entry-point."""
        import builtins, io
        real_open = builtins.open
        cases = {
            999_011: b"/home/u/.local/bin/splicecraft\x00",
            999_012: b"/usr/bin/python3\x00splicecraft.py\x00L09137\x00",
            999_013: b"/usr/bin/python3\x00-m\x00splicecraft\x00",
            999_014: (b"/home/u/.local/share/pipx/venvs/splicecraft/bin/"
                      b"python\x00/home/u/.local/share/pipx/venvs/"
                      b"splicecraft/bin/splicecraft\x00"),
        }

        def _make(pid, data):
            def _stub_open(path, *a, **k):
                if str(path) == f"/proc/{pid}/cmdline":
                    return io.BytesIO(data)
                return real_open(path, *a, **k)
            return _stub_open

        for pid, data in cases.items():
            monkeypatch.setattr(builtins, "open", _make(pid, data))
            result = sc._pid_is_splicecraft(pid)
            if result is None:
                pytest.skip("/proc/<pid>/cmdline path not used on this platform")
            assert result is True, f"missed real launch mode (pid {pid})"


# ═══════════════════════════════════════════════════════════════════════════════
# H11 — Clipboard image megapixel cap
# ═══════════════════════════════════════════════════════════════════════════════

class TestClipMegapixelCap:
    """A hostile / accidental multi-monitor screenshot would let
    `ImageGrab.grabclipboard()` decode a multi-GB bitmap into memory
    AND then write it to `/tmp` as PNG before the byte cap caught it.
    The new `_EXPERIMENT_CLIP_MAX_PIXELS` constant short-circuits
    before save."""

    def test_constant_exists_and_is_finite(self):
        assert hasattr(sc, "_EXPERIMENT_CLIP_MAX_PIXELS")
        cap = sc._EXPERIMENT_CLIP_MAX_PIXELS
        assert isinstance(cap, int)
        # Sanity bounds: must accept typical screenshots (~8 MP for
        # 4K) but reject a contrived multi-GB bitmap.
        assert 10_000_000 < cap < 1_000_000_000


# ═══════════════════════════════════════════════════════════════════════════════
# M1 — `_save_experiments` mirror runs inside `_cache_lock`
# ═══════════════════════════════════════════════════════════════════════════════

class TestSaveExperimentsMirrorInsideLock:
    """Pre-fix `_sync_active_project_experiments` ran AFTER the lock
    release in `_save_experiments`. Future concurrent writers could
    desync `experiments.json` vs `experiment_projects.json` in a
    persistent way. The fix moves the mirror INSIDE the `with`
    block (RLock supports re-entry)."""

    def test_save_experiments_source_holds_lock_through_mirror(self):
        import inspect
        src = inspect.getsource(sc._save_experiments)
        # The active-project mirror MUST fire before the cache-lock block
        # exits. Phase D: `_save_experiments` moved to splicecraft_dataaccess
        # and now (a) locks via `with _state._cache_lock:` and (b) fires the
        # mirror through `_state._sync_active_project_experiments_hook` (the
        # hub-side `_sync_active_project_experiments`, sweep #9 unchanged). The
        # invariant is the same: the mirror-fire appears at deeper indentation
        # than the `with` line, i.e. inside the locked block.
        lines = src.splitlines()
        with_indent = -1
        sync_indent = -1
        for ln in lines:
            if (with_indent < 0
                    and ln.lstrip().startswith("with ")
                    and "_cache_lock" in ln):
                with_indent = len(ln) - len(ln.lstrip())
            if (sync_indent < 0
                    and "_sync_active_project_experiments_hook" in ln):
                sync_indent = len(ln) - len(ln.lstrip())
        assert with_indent >= 0, "no `with ... _cache_lock` in source"
        assert sync_indent >= 0, "no active-project mirror hook fired"
        # Sync must be inside the `with` block: greater indent.
        assert sync_indent > with_indent, (
            "_sync_active_project_experiments must be INSIDE the "
            "_cache_lock with-block; sweep #9 fix regressed"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# M2 — `_check_data_files` covers new 0.9.6 files
# ═══════════════════════════════════════════════════════════════════════════════

class TestCheckDataFilesCoversNewFiles:
    """Without launch-time validation of the new files, corruption only
    surfaces on lazy first-load and the warning may not reach the user
    via `notify`."""

    def test_check_data_files_covers_new_files(self):
        # Sweep #26 (2026-05-23): `_check_data_files` no longer hand-
        # lists file attrs — it iterates `_USER_DATA_FILE_ATTRS`
        # (the canonical registry). Verify both that the function
        # source iterates the registry AND that the registry actually
        # contains the four sweep-#9 + sibling persisted files.
        import inspect
        src = inspect.getsource(sc.PlasmidApp._check_data_files)
        assert "_USER_DATA_FILE_ATTRS" in src, (
            "_check_data_files must drive its file list from the "
            "canonical _USER_DATA_FILE_ATTRS registry, not a hand-list"
        )
        for token in (
            "_EXPERIMENTS_FILE",
            "_EXPERIMENT_PROJECTS_FILE",
            "_GELS_FILE",
            "_PARTS_BIN_COLLECTIONS_FILE",
        ):
            assert token in sc._USER_DATA_FILE_ATTRS, (
                f"_USER_DATA_FILE_ATTRS missing {token}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# M3 — `_sweep_orphan_tmp_files` walks the experiments subdirs
# ═══════════════════════════════════════════════════════════════════════════════

class TestOrphanTmpSweepCoversExperiments:
    """Per-entry `_atomic_write_bytes` crashes leave `.tmp` files
    under `<DATA_DIR>/experiments/<entry_id>/`. The sweep must reach
    those subdirs (one level deep) so leaked tmpfiles don't slowly
    consume the per-entry 100 MB cap."""

    def test_sweep_finds_old_tmp_in_experiment_subdir(
        self, tmp_path, monkeypatch,
    ):
        import time as _time
        # Redirect both _DATA_DIR and _EXPERIMENTS_DIR so the sweep
        # walks our tmp_path instead of the user's actual data dir.
        monkeypatch.setattr(sc._state, "_DATA_DIR", tmp_path)
        exp_root = tmp_path / "experiments"
        exp_root.mkdir()
        monkeypatch.setattr(sc._state, "_EXPERIMENTS_DIR", exp_root)
        # Make a per-entry subdir with an orphan tmp file aged
        # beyond the cleanup threshold.
        entry_dir = exp_root / "exp-deadbeef"
        entry_dir.mkdir()
        tmp_file = entry_dir / "img.tmp"
        tmp_file.write_text("orphan")
        # Backdate so the >1 h age check passes.
        old = _time.time() - (sc._ORPHAN_TMP_MIN_AGE_S + 10)
        import os as _os
        _os.utime(tmp_file, (old, old))
        # Sweep should find and remove it.
        removed = sc._sweep_orphan_tmp_files(tmp_path)
        assert removed >= 1
        assert not tmp_file.exists()


# ═══════════════════════════════════════════════════════════════════════════════
# M6 — `_blocks_undo` on new modals (required for H3 to be effective)
# ═══════════════════════════════════════════════════════════════════════════════

class TestNewModalsBlockUndo:
    """After H3 (Ctrl+Z routed through public action_undo), modals that
    own text input MUST opt in via `_blocks_undo = True` or the
    keyboard shortcut will still unwind plasmid edits behind them."""

    @pytest.mark.parametrize("cls_name", [
        "ExperimentProjectsPickerModal",
        "ExperimentDeleteConfirmModal",
        "ExperimentUnsavedChangesModal",
        "ExperimentRenameModal",
        "GelLibraryModal",
        "ActionsPickerModal",
        "ImageAttachModal",
        "SpellcheckModal",
    ])
    def test_modal_has_blocks_undo(self, cls_name):
        cls = getattr(sc, cls_name)
        assert getattr(cls, "_blocks_undo", False) is True, (
            f"{cls_name} missing `_blocks_undo = True` (sweep #9)"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# M7 — `_summarize_perbase_tsv` accepts float `reads_all`
# ═══════════════════════════════════════════════════════════════════════════════

class TestPerbaseAcceptsFloat:
    """Pre-fix the first-line header probe used `int()` which raised
    ValueError on a float (`100\\tA\\t0.5`), discarding the first data
    row. Sweep #9 uses `float()` for the probe and `int(float(...))`
    for the value."""

    def test_first_line_with_float_reads_all_is_counted(self):
        # First data row has float reads_all — pre-fix this would be
        # treated as a header and skipped, dropping the first row.
        tsv = "100\tA\t10.0\n101\tC\t20\n102\tG\t30\n"
        fh = io.BytesIO(tsv.encode("utf-8"))
        out = sc._summarize_perbase_tsv(fh)
        # All 3 rows should be counted (n=3), not 2.
        assert out.get("n_pos") == 3

    def test_integer_reads_all_still_works(self):
        tsv = "100\tA\t10\n101\tC\t20\n"
        fh = io.BytesIO(tsv.encode("utf-8"))
        out = sc._summarize_perbase_tsv(fh)
        assert out.get("n_pos") == 2


# ═══════════════════════════════════════════════════════════════════════════════
# M8 — Legacy tag migration on save (not only on load)
# ═══════════════════════════════════════════════════════════════════════════════

class TestLegacyTagMigrationOnSave:
    """Pre-fix `_migrate_legacy_tag_format` ran only on load. A body
    that arrived in-memory via paste / import AFTER initial load would
    persist with the old format on save. Now `_normalise_experiment_entry`
    runs the migration in-place so saves always emit the single-sigil
    form."""

    def test_legacy_plasmid_token_rewritten_on_normalise(self):
        entry = {
            "id": "exp-12345678",
            "body_md": "See @plasmid:pBR322 for context.",
        }
        out = sc._normalise_experiment_entry(entry)
        assert "@plasmid:" not in out["body_md"]
        assert "@pBR322" in out["body_md"]

    def test_legacy_action_token_rewritten_on_normalise(self):
        entry = {
            "id": "exp-12345678",
            "body_md": "Used @actions:digest on the fragment.",
        }
        out = sc._normalise_experiment_entry(entry)
        assert "@actions:" not in out["body_md"]
        assert "!digest" in out["body_md"]

    def test_extracted_xref_uses_new_format_id(self):
        # After legacy migration on save, the denormalised
        # `attached_*` lists should contain the new-form ids.
        entry = {
            "id": "exp-12345678",
            "body_md": "Ref @plasmid:pBR322 here.",
        }
        out = sc._normalise_experiment_entry(entry)
        assert out["attached_plasmid_ids"] == ["pBR322"]


# ═══════════════════════════════════════════════════════════════════════════════
# M9 — Settings schema covers new persisted keys
# ═══════════════════════════════════════════════════════════════════════════════

class TestSettingsSchemaNewKeys:
    """The 0.8.x / 0.9.x added persisted keys (`active_parts_bin`,
    `active_project`, `experiments_custom_dict`) that weren't in
    `_SETTINGS_SCHEMA`. The unknown-key forward-compat passthrough
    caught them today, but the schema is the documented contract."""

    @pytest.mark.parametrize("key,expected_type,default", [
        ("active_parts_bin", str, ""),
        ("active_project", str, ""),
        ("experiments_custom_dict", list, []),
    ])
    def test_new_key_in_schema(self, key, expected_type, default):
        assert key in sc._SETTINGS_SCHEMA, (
            f"{key} missing from _SETTINGS_SCHEMA (sweep #9)"
        )
        allowed, default_val = sc._SETTINGS_SCHEMA[key]
        assert expected_type in allowed
        assert default_val == default

    def test_validate_settings_accepts_list_for_custom_dict(self):
        raw = {"experiments_custom_dict": ["foo", "bar"]}
        cleaned, warnings = sc._validate_settings(raw)
        assert cleaned["experiments_custom_dict"] == ["foo", "bar"]
        assert warnings == []

    def test_validate_settings_rejects_wrong_type_for_active_project(self):
        # Strict type check — an int value for a str-typed key should
        # land on the default, with a warning.
        raw = {"active_project": 42}
        cleaned, warnings = sc._validate_settings(raw)
        assert cleaned["active_project"] == ""
        assert any("active_project" in w for w in warnings)


# ═══════════════════════════════════════════════════════════════════════════════
# H7 — SequencingScreen same-path cache invalidates on mtime change
# ═══════════════════════════════════════════════════════════════════════════════

class TestSequencingZipCacheRespectsMtime:
    """The same-path short-circuit in `_on_zip_picked` previously
    skipped re-parse when the path matched, regardless of whether
    the file's content had changed. Real workflow: user re-runs
    Plasmidsaurus, overwrites the local zip — the Samples tab
    silently shows stale samples."""

    def test_screen_carries_signature_attribute(self):
        # Without an `_zip_signature` attribute on the SequencingScreen
        # instance, the new content-drift detection can't work.
        # Instantiate and confirm the attribute exists with the
        # documented default (`None`).
        ss = sc.SequencingScreen()
        assert hasattr(ss, "_zip_signature")
        assert ss._zip_signature is None

    def test_reset_clears_signature(self):
        ss = sc.SequencingScreen()
        ss._zip_signature = (123456789, 1024)
        ss._reset_zip_state()
        assert ss._zip_signature is None


# ═══════════════════════════════════════════════════════════════════════════════
# Wave 4 — Deferred Tier 1 fixes
# ═══════════════════════════════════════════════════════════════════════════════

class TestSigtermHandler:
    """Sweep #9 installs SIGTERM + SIGHUP handlers in `main()` that
    translate the signal to `KeyboardInterrupt` so the teardown
    `finally` block runs. Without this, systemd / container shutdowns
    skip lockfile release, settings flush, and log shutdown."""

    def test_main_source_installs_signal_handlers(self):
        import inspect
        src = inspect.getsource(sc.main)
        # Source must mention installing handlers for both signals;
        # the registration pattern is `signal.signal(SIGNAME, ...)`.
        assert "SIGTERM" in src
        assert "SIGHUP" in src
        assert "signal.signal" in src or "_signal.signal" in src


class TestNonAsciiTruncationSingleEncode:
    """Pre-fix `_normalise_experiment_entry`'s body-truncation loop
    re-encoded the entire body on every 1024-char shrink — quadratic
    on multi-MB non-ASCII content. Sweep #9 collapses to a single
    encode + slice + decode (`errors="ignore"` handles mid-multibyte
    truncation)."""

    def test_large_non_ascii_body_truncated_cleanly(self):
        # Build a body slightly over the byte cap using Chinese
        # characters (3 bytes each in UTF-8).
        big_char = "中"  # "中" — 3 bytes UTF-8
        body = big_char * (sc._EXPERIMENT_BODY_MAX_BYTES // 2)
        # That's ~1.5 MB of bytes — over the 1 MB cap.
        encoded = body.encode("utf-8")
        assert len(encoded) > sc._EXPERIMENT_BODY_MAX_BYTES
        entry = {"id": "exp-12345678", "body_md": body}
        out = sc._normalise_experiment_entry(entry)
        out_bytes = out["body_md"].encode("utf-8")
        # Must fit within cap and not include a torn multibyte
        # sequence (decode would have errored if any partial byte
        # leaked through, but `errors="ignore"` is the safety net).
        assert len(out_bytes) <= sc._EXPERIMENT_BODY_MAX_BYTES
        # Body should not be empty (we just truncated, not nuked).
        assert len(out["body_md"]) > 0


class TestGelInspectOnlyDisablesLoad:
    """Pre-fix the click-through `_open_gel_ref` path opened
    `GelLibraryModal` with callback=None, so the Load button silently
    no-op'd. The `inspect_only` flag now disables Load + Save in that
    mode so the user gets honest UI feedback."""

    def test_inspect_only_constructor_arg_accepted(self):
        m = sc.GelLibraryModal(inspect_only=True)
        assert m._inspect_only is True

    def test_default_constructor_is_not_inspect_only(self):
        m = sc.GelLibraryModal()
        assert m._inspect_only is False

    def test_actions_picker_also_supports_inspect_only(self):
        m = sc.ActionsPickerModal(initial_action="digest",
                                  inspect_only=True)
        assert m._inspect_only is True


class TestPersistCurrentRefreshesTextArea:
    """`ExperimentsScreen._persist_current` source must call
    `load_text` to re-sync the TextArea after a save-side normalise
    that may have truncated the body. Without this the on-screen
    text drifts from the on-disk version after the cap fires."""

    def test_source_calls_load_text_on_normalisation_diff(self):
        import inspect
        src = inspect.getsource(sc.ExperimentsScreen._persist_current)
        assert "load_text" in src
        assert "normalised_body" in src or "normalised body" in src


class TestDataDirWhitespaceEnv:
    """`splicecraft_cli._data_dir` must reject whitespace-only
    `$SPLICECRAFT_DATA_DIR` so a shell scripting bug doesn't land
    us on a literal whitespace-named relative directory."""

    def test_whitespace_only_env_falls_back_to_platform_default(
        self, monkeypatch,
    ):
        from splicecraft_cli import _data_dir
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", "   ")
        result = _data_dir()
        # Must not be the whitespace literal. The actual default
        # depends on the platform (`platformdirs` or `~/.local/...`),
        # so just check we got "something not whitespace".
        assert str(result).strip()
        assert str(result) != "   "
        assert str(result) != " "

    def test_legitimate_env_is_honoured(self, monkeypatch, tmp_path):
        from splicecraft_cli import _data_dir
        target = tmp_path / "custom-data"
        monkeypatch.setenv("SPLICECRAFT_DATA_DIR", str(target))
        result = _data_dir()
        assert result == target


class TestLibrarySearchNoMatchNotify:
    """`LibrarySearchModal._on_query_submitted` source must contain
    a notify call for the no-match path so the user gets feedback
    instead of silent no-op."""

    def test_source_notifies_on_empty_matches(self):
        import inspect
        src = inspect.getsource(
            sc.LibrarySearchModal._on_query_submitted,
        )
        assert "notify" in src
        assert "self._matches" in src


class TestExperimentsDeleteStickyCursor:
    """`_confirm_delete` source must capture cursor row BEFORE the
    confirm modal and restore via `move_cursor` after delete +
    repopulate."""

    def test_source_captures_cursor_row_before_callback(self):
        import inspect
        src = inspect.getsource(
            sc.ExperimentsScreen._confirm_delete,
        )
        assert "cur_row_before" in src
        assert "move_cursor" in src


class TestBlastCacheLock:
    """`_BLAST_CACHE_LOCK` exists and is used by both
    `_blast_get_db` and `_blast_clear_cache`."""

    def test_lock_exists(self):
        assert hasattr(sc, "_BLAST_CACHE_LOCK")
        # Sweep #27 (LOW-11): converted to RLock for convention parity
        # with `_cache_lock`. Either lock type is acceptable for the
        # critical sections (none re-enter today), but RLock is the
        # codebase-wide standard so a future caller that does end up
        # nesting can compose without deadlock risk. Accept either
        # type so a hypothetical revert wouldn't trip the test, but
        # do confirm it IS a lock-like object with `acquire`/`release`.
        lock = sc._BLAST_CACHE_LOCK
        assert hasattr(lock, "acquire") and hasattr(lock, "release")
        # Probe the context-manager interface — both Lock and RLock
        # support it; this catches a regression to a non-lock object.
        with lock:
            pass

    def test_get_db_source_uses_lock(self):
        import inspect
        src = inspect.getsource(sc._blast_get_db)
        assert "_BLAST_CACHE_LOCK" in src

    def test_clear_source_uses_lock(self):
        import inspect
        src = inspect.getsource(sc._blast_clear_cache)
        assert "_BLAST_CACHE_LOCK" in src


class TestSettingsCacheDeepcopy:
    """`_set_setting` must use `_typed_clone` (deepcopy) when
    reseating `_settings_cache`, not a shallow `dict(...)`. Without
    this, a nested-list value (`experiments_custom_dict`) would
    share the caller's reference."""

    def test_set_setting_source_uses_typed_clone(self):
        import inspect
        src = inspect.getsource(sc._set_setting)
        assert "_typed_clone(settings)" in src
        # Pre-fix used `dict(settings)` — must NOT appear now.
        # (Defensive — if a future change re-introduces the shallow
        # copy, this test catches it.)
        assert "_settings_cache = dict(settings)" not in src


class TestAgentSetSettingStrictBool:
    """`_settings_validator_int_range` must reject bool values so
    `set-setting min_primer_binding=true` is refused instead of
    silently coerced to `1`."""

    def test_validator_rejects_bool(self):
        validator = sc._settings_validator_int_range(1, 60)
        result, err = validator(True)
        assert result is None
        assert err is not None
        assert "bool" in err.lower()

    def test_validator_accepts_real_int(self):
        validator = sc._settings_validator_int_range(1, 60)
        result, err = validator(15)
        assert result == 15
        assert err is None


class TestRotateSeqRecordDeadCodeGone:
    """`_rotate_seq_record` had an unreachable `if new_e == 0 and
    new_s + flen == n` branch removed in sweep #9. Removing the
    dead branch must not change rotation correctness."""

    def test_rotation_handles_boundary_offset(self):
        # Functional check (not source-grep): with `flen == n` and
        # `offset == 0` the ternary correctly resolves new_e to n.
        # Pre-fix the unreachable if-branch was a no-op so output
        # is identical either way; this just locks in that the
        # rotation logic still works at the boundary.
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        from Bio.Seq import Seq
        seq = Seq("A" * 100)
        rec = SeqRecord(
            seq=seq, id="t", name="t",
            features=[
                SeqFeature(
                    FeatureLocation(0, 100, strand=1),
                    type="misc_feature",
                ),
            ],
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        rotated = sc._rotate_seq_record(rec, 0)
        assert rotated is not None
        # Feature must still span the full 100 bp post-rotation.
        f = rotated.features[0]
        s = int(f.location.start)
        e = int(f.location.end)
        assert (s, e) == (0, 100)


# ═══════════════════════════════════════════════════════════════════════════════
# Sacred-invariant labelled regression tests (V1_GATE hard gate #3)
# ═══════════════════════════════════════════════════════════════════════════════

class TestInvariant4PatternCacheBounded:
    """Sacred invariant #4: IUPAC regex patterns are cached in
    `_PATTERN_CACHE`, bounded by `_PATTERN_CACHE_MAX` (256 default).
    A long-lived process with thousands of distinct custom enzymes
    must NOT grow the cache without bound."""

    def test_pattern_cache_max_is_finite(self):
        from splicecraft_biology import _PATTERN_CACHE_MAX
        assert isinstance(_PATTERN_CACHE_MAX, int)
        assert 0 < _PATTERN_CACHE_MAX < 100_000

    def test_pattern_cache_evicts_above_cap(self):
        from splicecraft_biology import (
            _PATTERN_CACHE, _PATTERN_CACHE_MAX, _iupac_pattern,
        )
        # Snapshot + clear so we can measure deterministically.
        snapshot = dict(_PATTERN_CACHE)
        _PATTERN_CACHE.clear()
        try:
            # Compile cap+50 distinct unique recognition strings.
            # Use a deterministic generator: stride through IUPAC
            # 4-letter combinations to guarantee uniqueness without
            # collisions.
            n = _PATTERN_CACHE_MAX + 50
            iupac = "ACGTRYSWKMBDHVN"
            seen: set[str] = set()
            i = 0
            for a in iupac:
                for b in iupac:
                    for c in iupac:
                        for d in iupac:
                            pat = a + b + c + d
                            if pat in seen:
                                continue
                            seen.add(pat)
                            _iupac_pattern(pat)
                            i += 1
                            if i >= n:
                                break
                        if i >= n:
                            break
                    if i >= n:
                        break
                if i >= n:
                    break
            # Cache MUST NOT exceed the documented max.
            assert len(_PATTERN_CACHE) <= _PATTERN_CACHE_MAX
        finally:
            _PATTERN_CACHE.clear()
            _PATTERN_CACHE.update(snapshot)


class TestInvariant12ReHighlightSchemaCurrent:
    """Sacred invariant #12 (0.4.5+): `_re_highlight` schema is
    `{start, end, top_cut_bp, bottom_cut_bp, color, name}`. Legacy
    `fwd_cut_bp` / `rev_cut_bp` keys must NOT survive."""

    def test_no_callsite_uses_legacy_fwd_cut_bp_key(self):
        # Grep source for the legacy keys — they must not appear
        # anywhere in production code.
        src_path = sc.__file__
        with open(src_path, encoding="utf-8") as fh:
            text = fh.read()
        # Allow appearance in COMMENTS (sweep history may
        # reference them) but reject any actual dict-key usage.
        # Heuristic: legacy keys followed by `:` indicate dict
        # literal usage.
        assert '"fwd_cut_bp":' not in text
        assert '"rev_cut_bp":' not in text


class TestInvariant29IdLessEntriesSkipped:
    """Sacred invariant #29: cross-collection search skips id-less
    entries. Without this filter, a dismiss payload `(collection,
    "")` would alias every untagged entry to the first one in the
    active library on load."""

    def test_search_skips_entries_without_id(
        self, tmp_path, monkeypatch,
    ):
        # Seed two collections; one has an entry with empty id, the
        # other with a valid id matching the query.
        monkeypatch.setattr(sc._state, "_collections_cache", None)
        monkeypatch.setattr(sc._state, "_library_cache", None)
        valid_collection = {
            "name": "Valid",
            "plasmids": [{
                "id": "pBR322", "name": "pBR322",
                "gb_text": "LOCUS pBR322 4361 bp\n//\n",
            }],
        }
        empty_id_collection = {
            "name": "EmptyId",
            "plasmids": [{
                "id": "",  # empty — must be skipped
                "name": "pBR322",
                "gb_text": "LOCUS something 4361 bp\n//\n",
            }],
        }
        monkeypatch.setattr(
            sc, "_load_collections",
            lambda: [valid_collection, empty_id_collection],
        )
        results = sc._search_collections_library("pBR322", limit=10)
        # Only the valid-id entry should appear.
        for r in results:
            assert r.get("id"), f"id-less entry surfaced: {r!r}"


class TestInvariant34ClassifyMultiGrammarConsistent:
    """Sacred invariant #34/#40: `_classify_part_from_plasmid` loops
    over `_all_grammars()` in REGISTRY order and picks the first
    grammar whose Type IIS digest produces an overhang pair matching
    one of its positions. Two grammars sharing an overhang pair → the
    first-registered grammar's classification wins (deterministic)."""

    def test_classify_is_deterministic_across_calls(self):
        # Random non-grammar sequence — no digest will match any
        # grammar's overhang table, so the classifier returns None.
        # Determinism check: two calls return the same result.
        seq = "ATGCATGC" * 50
        a = sc._classify_part_from_plasmid(seq, circular=True)
        b = sc._classify_part_from_plasmid(seq, circular=True)
        assert a == b

    def test_classify_linear_skipped(self):
        # Linear records skip the digest path entirely per the
        # function's contract.
        seq = "ATGCATGC" * 50
        result = sc._classify_part_from_plasmid(seq, circular=False)
        assert result is None
