"""
test_experiment_projects — multi-project lab-notebook storage.

Covers the data layer that backs the experiment-projects modal:

  * `_load_experiment_projects` / `_save_experiment_projects`
    round-trip schema preservation + cache hygiene.
  * `_ensure_default_project` migration (existing `experiments.json`
    contents wrap into "Main Project" on first launch; idempotent
    on subsequent launches).
  * Active-project pointer setter / getter / `_find_project` /
    `_project_name_taken`.
  * `_save_experiments` mirrors into the active project (sacred
    contract: the multi-project record never drifts from
    `experiments.json`).
  * `_sync_active_project_experiments` is a silent no-op when there
    is no active project (first-launch race).
"""
from __future__ import annotations

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Round-trip + cache hygiene
# ═══════════════════════════════════════════════════════════════════════════════

class TestExperimentProjectsRoundTrip:
    """`_save_experiment_projects` / `_load_experiment_projects`
    preserve the full schema (including unknown forward-compat keys)."""

    def test_empty_load(self):
        # Fresh `_EXPERIMENT_PROJECTS_FILE` (autouse conftest redirects
        # to tmp) returns an empty list, no exception.
        assert sc._load_experiment_projects() == []

    def test_round_trip_preserves_fields(self):
        entries = [{
            "name": "Yeast project",
            "description": "Saccharomyces work",
            "experiments": [
                {"id": "exp-aaaaaaaa", "title": "Day 1", "body_md": "Hello"},
            ],
            "saved": "2026-05-12",
            "_plugin_data": {"some_plugin": {"x": 1}},  # reserved field
        }]
        sc._save_experiment_projects(entries)
        loaded = sc._load_experiment_projects()
        assert loaded == entries

    def test_load_deepcopies(self):
        """Per invariant #17 — caller mutations after load must not
        poison the in-memory cache."""
        sc._save_experiment_projects([
            {"name": "A", "description": "", "experiments": [], "saved": ""},
        ])
        first = sc._load_experiment_projects()
        first[0]["name"] = "MUTATED"
        second = sc._load_experiment_projects()
        assert second[0]["name"] == "A"

    def test_save_deepcopies(self):
        """Per invariant #17 — caller mutations after save must not
        leak into the next load via shared dict refs."""
        entries = [{"name": "A", "description": "",
                    "experiments": [], "saved": ""}]
        sc._save_experiment_projects(entries)
        entries[0]["name"] = "MUTATED"
        loaded = sc._load_experiment_projects()
        assert loaded[0]["name"] == "A"


# ═══════════════════════════════════════════════════════════════════════════════
# Active-project pointer
# ═══════════════════════════════════════════════════════════════════════════════

class TestActiveProjectPointer:
    def test_initial_value_is_none(self):
        assert sc._get_active_project_name() is None

    def test_set_and_get(self):
        sc._set_active_project_name("Yeast project")
        assert sc._get_active_project_name() == "Yeast project"

    def test_clear(self):
        sc._set_active_project_name("A")
        sc._set_active_project_name(None)
        assert sc._get_active_project_name() is None

    def test_find_project(self):
        sc._save_experiment_projects([
            {"name": "A", "description": "", "experiments": [], "saved": ""},
            {"name": "B", "description": "", "experiments": [], "saved": ""},
        ])
        assert sc._find_project("A")["name"] == "A"
        assert sc._find_project("B")["name"] == "B"
        assert sc._find_project("C") is None

    def test_name_taken(self):
        sc._save_experiment_projects([
            {"name": "A", "description": "", "experiments": [], "saved": ""},
        ])
        assert sc._project_name_taken("A") is True
        assert sc._project_name_taken("B") is False


class TestExperimentMirrorHardening:
    """Startup restore + set-active revert for experiment projects — the
    experiments twins of the primer-collection mirror hardening
    ([INV-161]/[INV-83]). Experiments previously had NO launch restore, so a
    switch/move-experiment whose live re-sync failed left a permanent desync."""

    _A = {"name": "A", "description": "", "saved": "",
          "experiments": [{"id": "e1", "body": "from-A", "tags": []}]}
    _B = {"name": "B", "description": "", "saved": "",
          "experiments": [{"id": "e2", "body": "from-B", "tags": []}]}

    def test_restore_rebuilds_live_file_from_active_project(self):
        import json
        sc._save_experiment_projects([dict(self._A)])
        sc._set_active_project_name("A")
        # Stash a stale live experiments.json directly, then assert the restore
        # rebuilds it from the active project (source of truth).
        sc._safe_save_json_mirror(
            sc._state._EXPERIMENTS_FILE,
            [{"id": "stale", "body": "stale", "tags": []}], "Experiments")
        sc._restore_experiments_from_active_project()
        with open(sc._state._EXPERIMENTS_FILE) as fh:
            data = json.load(fh)
        entries = data.get("entries", data) if isinstance(data, dict) else data
        ids = {e.get("id") for e in entries}
        assert ids == {"e1"} and "stale" not in ids

    def test_agent_set_active_reverts_on_mirror_failure(self, monkeypatch):
        import splicecraft_agent as sca
        sc._save_experiment_projects([dict(self._A), dict(self._B)])
        sc._h_set_active_experiment_project(None, {"name": "A"})
        def _boom(*a, **k):
            raise OSError("simulated disk failure")
        monkeypatch.setattr(sca, "_safe_save_json_mirror", _boom)
        r = sc._h_set_active_experiment_project(None, {"name": "B"})
        assert isinstance(r, tuple) and r[1] == 500
        assert (sc._get_active_project_name() or "") == "A"   # reverted


# ═══════════════════════════════════════════════════════════════════════════════
# Migration
# ═══════════════════════════════════════════════════════════════════════════════

class TestEnsureDefaultProject:
    """First-run migration: wrap existing `experiments.json` into a
    "Main Project" wrapper. Idempotent on subsequent calls."""

    def test_empty_first_run_creates_empty_main_project(self):
        # No experiments.json contents, no experiment_projects.json yet.
        assert sc._load_experiments() == []
        assert sc._load_experiment_projects() == []
        sc._ensure_default_project()
        projs = sc._load_experiment_projects()
        assert len(projs) == 1
        assert projs[0]["name"] == sc._DEFAULT_PROJECT_NAME
        assert projs[0]["experiments"] == []
        assert sc._get_active_project_name() == sc._DEFAULT_PROJECT_NAME

    def test_existing_entries_wrap_into_main_project(self):
        # Seed `experiments.json` with two entries BEFORE migration runs.
        sc._save_experiments([
            {"id": "exp-aaaaaaaa", "title": "Day 1",
             "body_md": "First entry"},
            {"id": "exp-bbbbbbbb", "title": "Day 2",
             "body_md": "Second entry"},
        ])
        # Migration not yet run — projects file is empty.
        # (`_save_experiments` mirrors silently when no active project.)
        assert sc._load_experiment_projects() == []
        sc._ensure_default_project()
        projs = sc._load_experiment_projects()
        assert len(projs) == 1
        assert projs[0]["name"] == sc._DEFAULT_PROJECT_NAME
        assert len(projs[0]["experiments"]) == 2
        assert projs[0]["experiments"][0]["id"] == "exp-aaaaaaaa"

    def test_idempotent_on_subsequent_calls(self):
        sc._ensure_default_project()
        before = sc._load_experiment_projects()
        sc._ensure_default_project()
        sc._ensure_default_project()
        after = sc._load_experiment_projects()
        assert before == after

    def test_restores_active_pointer_if_lost(self):
        # Projects exist on disk but the active-project setting was
        # cleared (could happen if a user hand-edited settings.json).
        sc._save_experiment_projects([
            {"name": "Custom", "description": "",
             "experiments": [], "saved": ""},
        ])
        sc._set_active_project_name(None)
        sc._ensure_default_project()
        # Should adopt the first existing project's name, not blow
        # away the user's custom project.
        assert sc._get_active_project_name() == "Custom"
        assert len(sc._load_experiment_projects()) == 1

    def test_repairs_orphaned_active_pointer(self):
        """Edge sweep 2026-05-18 — `settings.json::active_project`
        names a project that doesn't exist (renamed/deleted by hand).
        The mirror in `_sync_active_project_experiments` would
        silently skip on every save until the pointer is fixed.
        `_ensure_default_project` MUST reset to a real project."""
        sc._save_experiment_projects([
            {"name": "Real", "description": "",
             "experiments": [], "saved": ""},
        ])
        sc._set_active_project_name("Ghost")   # doesn't exist
        sc._ensure_default_project()
        # Repaired to the first existing project.
        assert sc._get_active_project_name() == "Real"


# ═══════════════════════════════════════════════════════════════════════════════
# Picker modal — delete edges
# ═══════════════════════════════════════════════════════════════════════════════

class TestPickerDeleteEdges:
    """Defensive checks on `ExperimentProjectsPickerModal._do_delete`
    (edge sweep 2026-05-18)."""

    async def test_do_delete_refuses_when_list_shrinks_mid_confirm(self):
        """Between the user clicking Delete (which opens the confirm
        modal) and clicking Yes, another path might have already
        deleted projects. `_do_delete` must re-check the
        last-project guard rather than blindly saving an empty
        list."""
        sc._save_experiment_projects([
            {"name": "A", "description": "", "experiments": [], "saved": ""},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause()
            modal = sc.ExperimentProjectsPickerModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause()
            # Simulate confirm fired but list now has only 1 entry.
            modal._do_delete("A")
            await pilot.pause()
            # Project still exists — refused.
            projs = sc._load_experiment_projects()
            assert {p["name"] for p in projs} == {"A"}

    async def test_do_delete_refuses_when_target_vanished(self):
        """Target project deleted by another path before confirm —
        `_do_delete` should refuse rather than save a no-op."""
        sc._save_experiment_projects([
            {"name": "A", "description": "", "experiments": [], "saved": ""},
            {"name": "B", "description": "", "experiments": [], "saved": ""},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause()
            modal = sc.ExperimentProjectsPickerModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause()
            modal._do_delete("Ghost")  # never existed
            await pilot.pause()
            projs = sc._load_experiment_projects()
            # Nothing got deleted.
            assert {p["name"] for p in projs} == {"A", "B"}

    async def test_do_delete_holds_cache_lock_across_rmw(
            self, monkeypatch,
    ):
        """Sweep #35 (2026-05-26): the full read→save→promote→mirror
        sequence in `_do_delete` runs under `_cache_lock`. Pre-fix,
        each helper acquired the lock individually and released
        between steps — a concurrent agent endpoint adding a new
        project between our load and save would have its addition
        silently overwritten by the `remaining` list computed pre-add.

        Verification strategy: monkeypatch `_save_experiment_projects`
        so that DURING the save call, a sibling thread tries
        non-blocking acquire of `_cache_lock`. If the outer
        `with _cache_lock:` is in place, the sibling thread can't
        get the lock (it's held by `_do_delete`'s thread). If the
        outer lock isn't held, the sibling acquires immediately —
        which would be the bug signature.
        """
        import threading as _threading

        sc._save_experiment_projects([
            {"name": "A", "description": "", "experiments": [], "saved": ""},
            {"name": "B", "description": "", "experiments": [], "saved": ""},
        ])
        sibling_could_acquire: list[bool] = []
        orig_save = sc._save_experiment_projects

        def _detect(entries):
            # Probe from a separate thread — non-blocking acquire
            # tells us whether `_cache_lock` is currently held by
            # the caller (the `_do_delete` thread).
            def _probe():
                got = sc._cache_lock.acquire(blocking=False)
                sibling_could_acquire.append(got)
                if got:
                    sc._cache_lock.release()
            t = _threading.Thread(target=_probe)
            t.start()
            t.join(timeout=2.0)
            return orig_save(entries)

        monkeypatch.setattr(sc, "_save_experiment_projects", _detect)
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause()
            modal = sc.ExperimentProjectsPickerModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause()
            modal._do_delete("A")
            await pilot.pause()
        assert sibling_could_acquire, (
            "Detector never fired — `_save_experiment_projects` was "
            "not called from `_do_delete`."
        )
        # Sibling MUST NOT have been able to acquire — the outer
        # `with _cache_lock:` in `_do_delete` should keep it locked.
        assert not any(sibling_could_acquire), (
            "Sibling thread acquired `_cache_lock` during the "
            "_do_delete RMW — the outer lock is missing or "
            "released too early."
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Mirror sync (sacred contract: every `_save_experiments` updates active)
# ═══════════════════════════════════════════════════════════════════════════════

class TestMirrorSync:
    def test_save_experiments_mirrors_into_active(self):
        sc._ensure_default_project()   # creates Main Project + sets active
        assert sc._get_active_project_name() == sc._DEFAULT_PROJECT_NAME
        sc._save_experiments([
            {"id": "exp-aaaaaaaa", "title": "Day 1", "body_md": "Hello"},
        ])
        projs = sc._load_experiment_projects()
        main = next(p for p in projs
                    if p["name"] == sc._DEFAULT_PROJECT_NAME)
        assert len(main["experiments"]) == 1
        assert main["experiments"][0]["id"] == "exp-aaaaaaaa"

    def test_save_with_no_active_project_is_silent_noop(self):
        """First-launch race: `_save_experiments` may fire before
        `_ensure_default_project` runs (test fixtures, agent API).
        The mirror call should silently do nothing rather than crash
        or create a phantom project."""
        sc._set_active_project_name(None)
        assert sc._load_experiment_projects() == []
        # This must not raise.
        sc._save_experiments([
            {"id": "exp-cccccccc", "title": "Free-floating",
             "body_md": ""},
        ])
        # And it must not have created a project out of nowhere.
        assert sc._load_experiment_projects() == []
        # The experiments file itself was still written.
        loaded = sc._load_experiments()
        assert len(loaded) == 1
        assert loaded[0]["id"] == "exp-cccccccc"

    def test_save_mirrors_into_only_the_active_project(self):
        # Two projects; only the active one should pick up new entries.
        sc._save_experiment_projects([
            {"name": "A", "description": "", "experiments": [], "saved": ""},
            {"name": "B", "description": "",
             "experiments": [{"id": "exp-kept", "title": "kept",
                              "body_md": ""}],
             "saved": ""},
        ])
        sc._set_active_project_name("A")
        sc._save_experiments([
            {"id": "exp-new-in-A", "title": "new", "body_md": ""},
        ])
        projs = sc._load_experiment_projects()
        a = next(p for p in projs if p["name"] == "A")
        b = next(p for p in projs if p["name"] == "B")
        assert len(a["experiments"]) == 1
        assert a["experiments"][0]["id"] == "exp-new-in-A"
        # B is untouched.
        assert len(b["experiments"]) == 1
        assert b["experiments"][0]["id"] == "exp-kept"

    def test_save_mirror_skips_when_active_name_was_deleted(self):
        """User deletes the active project out from under us (e.g.
        via the modal concurrently with a save). The mirror should
        silently skip rather than recreate the deleted project."""
        sc._save_experiment_projects([
            {"name": "A", "description": "", "experiments": [], "saved": ""},
        ])
        sc._set_active_project_name("Ghost")   # points at a deleted project
        sc._save_experiments([
            {"id": "exp-aaaaaaaa", "title": "x", "body_md": ""},
        ])
        projs = sc._load_experiment_projects()
        # No "Ghost" entry materialised.
        assert {p["name"] for p in projs} == {"A"}
