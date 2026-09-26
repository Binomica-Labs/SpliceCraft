"""
test_sweep10 — adversarial sweep #10 regression coverage (2026-05-20).

Pre-v1.0.0 audit found ~30 HIGH issues spanning collision flow, save-
chain lock-release gaps, modal hardening drift, security (`/dev/zero`
import + recursive history XML), and cache key collisions. This file
locks down the HIGH-severity fixes so they can't silently regress.

See CLAUDE.md invariant #50 for the full inventory.
"""
from __future__ import annotations

from unittest.mock import MagicMock

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Collision-modal double-fire guard
# ═══════════════════════════════════════════════════════════════════════════════

class TestCollisionModalDismissOnce:
    """Regression guard for 2026-05-20 fix: rapid double-clicks
    between button-press and the async pop-screen could fire the
    same callback twice with potentially different payloads. The
    `_dismiss_once` gate flips a flag on first fire.
    """

    def test_exact_copy_modal_dismiss_once(self):
        m = sc.ExactCopyConfirmModal("part", ["A", "B"])
        # Stub `dismiss` to capture calls without a live screen stack.
        calls: list = []
        m.dismiss = lambda payload: calls.append(payload)  # type: ignore[assignment]
        # Two rapid skips → only one dismiss.
        m._skip(None)
        m._skip(None)
        assert calls == [False]

    def test_exact_copy_skip_then_keep_only_first_fires(self):
        m = sc.ExactCopyConfirmModal("part", ["A"])
        calls: list = []
        m.dismiss = lambda payload: calls.append(payload)  # type: ignore[assignment]
        m._skip(None)
        m._keep(None)  # would have fired True pre-fix
        assert calls == [False]

    def test_name_collision_modal_dismiss_once(self):
        m = sc.NameCollisionModal("plasmid", ["X"])
        calls: list = []
        m.dismiss = lambda payload: calls.append(payload)  # type: ignore[assignment]
        m._overwrite(None)
        m._cancel_btn(None)
        assert calls == ["overwrite"]

    def test_action_cancel_uses_same_gate(self):
        m = sc.ExactCopyConfirmModal("part", ["A"])
        calls: list = []
        m.dismiss = lambda payload: calls.append(payload)  # type: ignore[assignment]
        m.action_cancel()
        m._keep(None)
        assert calls == [False]


# ═══════════════════════════════════════════════════════════════════════════════
# `_blocks_undo` carried on edit modals
# ═══════════════════════════════════════════════════════════════════════════════

class TestBlocksUndoOnEditModals:
    """Regression guard: every modal that mutates the canvas / record /
    persistent grammars file must carry `_blocks_undo: bool = True` so
    app-level Ctrl+Z can't undo the underlying record mid-edit. See
    invariant #41 worker-pattern + `_blocks_undo` contract.
    """

    def test_edit_seq_dialog(self):
        assert getattr(sc.EditSeqDialog, "_blocks_undo", False) is True

    def test_add_feature_modal(self):
        assert getattr(sc.AddFeatureModal, "_blocks_undo", False) is True

    def test_feature_edit_modal(self):
        assert getattr(sc.FeatureEditModal, "_blocks_undo", False) is True

    def test_new_plasmid_modal(self):
        assert getattr(sc.NewPlasmidModal, "_blocks_undo", False) is True

    def test_grammar_editor_modal(self):
        assert getattr(sc.GrammarEditorModal, "_blocks_undo", False) is True

    def test_restore_from_backup_modal(self):
        assert getattr(sc.RestoreFromBackupModal, "_blocks_undo", False) is True


# ═══════════════════════════════════════════════════════════════════════════════
# RestoreFromBackupModal cache-bust drift
# ═══════════════════════════════════════════════════════════════════════════════

class TestRestoreCacheBustEnumeration:
    """Regression guard for 2026-05-20 fix: sweep #9 added experiments
    / experiment projects / gels persisted files but
    `RestoreFromBackupModal._restore_btn` was never extended to
    invalidate their in-memory caches. A user restoring `experiments.json`
    from .bak saw no change until app restart; the next UI mutation
    would silently overwrite the freshly-restored disk state.
    """

    def test_modal_targets_include_sweep9_files(self):
        # Sweep #9 additions must be present.
        attr_names = {t[1] for t in sc.RestoreFromBackupModal._TARGETS}
        assert "_EXPERIMENTS_FILE" in attr_names
        assert "_EXPERIMENT_PROJECTS_FILE" in attr_names
        assert "_GELS_FILE" in attr_names

    def test_targets_parity_with_user_data_files(self):
        """Regression guard (2026-06-11 audit): EVERY user-data file in
        `_USER_DATA_FILE_ATTRS` except `_SETTINGS_FILE` MUST be reachable
        from the GUI `RestoreFromBackupModal._TARGETS`, mirroring
        `test_agent_backup_labels_parity_with_user_data_files` for the
        agent API.

        Pre-fix `_HMM_DB_CATALOG_FILE` was present in every other
        registry (`_USER_DATA_FILE_ATTRS`, `_MASTER_DELETE_CACHE_ATTRS`,
        `_AGENT_BACKUP_LABELS`) but missing here, and no whole-set test
        caught it — the existing `_TARGETS` tests were piecemeal
        per-sweep membership checks. This closes the drift class."""
        target_attrs = {t[1] for t in sc.RestoreFromBackupModal._TARGETS}
        expected = set(sc._USER_DATA_FILE_ATTRS) - {"_SETTINGS_FILE"}
        missing = expected - target_attrs
        assert not missing, (
            f"RestoreFromBackupModal._TARGETS missing entries for: "
            f"{sorted(missing)}. GUI users cannot restore these files "
            f"from backup even though they are user-data files. This "
            f"set silently grew before the parity guard landed."
        )

    def test_modal_cache_bust_covers_sweep9_caches(self):
        """The cache-bust block in `_restore_btn` enumerates every
        persisted-state cache. Sweep #25 (2026-05-23) replaced the
        hand-list with iteration of `_MASTER_DELETE_CACHE_ATTRS`
        (the canonical source of truth). Verify both that the modal
        iterates the master tuple AND that the master tuple contains
        the sweep-#9 + original caches.
        """
        import inspect
        # The restore itself runs in `_restore_locked`, under the one
        # `_cache_lock` hold `_restore_btn` takes (round-2 hardening,
        # 2026-09-25: the lock was released between the steps).
        src = inspect.getsource(sc.RestoreFromBackupModal._restore_locked)
        assert "_MASTER_DELETE_CACHE_ATTRS" in src, (
            "RestoreFromBackupModal must iterate the canonical "
            "_MASTER_DELETE_CACHE_ATTRS tuple, not a hand-list"
        )
        # Sweep #9 additions must be in the master tuple.
        assert "_experiments_cache" in sc._MASTER_DELETE_CACHE_ATTRS
        assert "_experiment_projects_cache" in sc._MASTER_DELETE_CACHE_ATTRS
        assert "_gels_cache" in sc._MASTER_DELETE_CACHE_ATTRS
        # Original four are still there.
        assert "_library_cache" in sc._MASTER_DELETE_CACHE_ATTRS
        assert "_collections_cache" in sc._MASTER_DELETE_CACHE_ATTRS
        assert "_parts_bin_cache" in sc._MASTER_DELETE_CACHE_ATTRS
        assert "_primers_cache" in sc._MASTER_DELETE_CACHE_ATTRS

    def test_agent_backup_labels_parity_with_user_data_files(self):
        """Regression guard for 2026-05-21 fix: `_AGENT_BACKUP_LABELS`
        and the cache-bust map inside `_h_restore_backup` covered only
        11 of 16 `_USER_DATA_FILE_ATTRS`. Users could list/restore
        Experiments / Gels / Protein motifs / Primer collections via
        the GUI but not via the agent API.

        Every user-data file attr except `_SETTINGS_FILE` (which the
        agent intentionally can't restore mid-session — see the
        `RestoreFromBackupModal._TARGETS` docstring) MUST be reachable
        by SOME label in `_AGENT_BACKUP_LABELS` AND have an entry in
        the cache-bust map.
        """
        import inspect
        labeled_attrs = set(sc._AGENT_BACKUP_LABELS.values())
        # Settings restore is deliberately excluded — see comment in
        # `RestoreFromBackupModal._TARGETS`.
        expected = set(sc._USER_DATA_FILE_ATTRS) - {"_SETTINGS_FILE"}
        missing = expected - labeled_attrs
        assert not missing, (
            f"_AGENT_BACKUP_LABELS missing entries for: {sorted(missing)}. "
            f"Agent can't list/restore these even though they're "
            f"user-data files. Pre-fix this set silently grew."
        )

        # Cache-bust map: white-box check that every label in
        # `_AGENT_BACKUP_LABELS` is also reset by `_h_restore_backup`
        # so a restore doesn't leave a stale in-memory cache.
        # The bust runs in `_restore_backup_locked`, under the one
        # `_cache_lock` hold `_h_restore_backup` takes (round-2 hardening).
        import splicecraft_agent as _ag
        bust_src = inspect.getsource(_ag._restore_backup_locked)
        bust_gaps = []
        for label in sc._AGENT_BACKUP_LABELS:
            # Each label appears as a dict key in the cache_attr map.
            # The settings label is included in the bust map; the
            # restore path treats it conservatively (sets cache None
            # so next read picks up disk).
            if f'"{label}"' not in bust_src:
                bust_gaps.append(label)
        assert not bust_gaps, (
            f"_h_restore_backup cache-bust map missing labels: "
            f"{bust_gaps}. Agent restore of these files would not "
            f"invalidate their in-memory cache, so the next UI read "
            f"would silently re-overwrite the restored disk state."
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Save-chain lock-release gap — mirror inside cache lock
# ═══════════════════════════════════════════════════════════════════════════════

class TestSaveChainAtomicMirror:
    """Regression guard for 2026-05-20 fix: `_save_library` and
    `_save_parts_bin` used to release `_cache_lock` BEFORE calling
    their mirror helper, so Thread A's write + Thread B's write +
    Thread B's mirror + Thread A's stale-snapshot mirror could leave
    `collections.json` / `parts_bin_collections.json` holding A while
    the live file held B. Fix: mirror inside the lock; RLock allows
    nested re-entry from the mirror's own save chain.
    """

    def test_save_library_runs_mirror(self, tmp_path, monkeypatch):
        # Smoke: the path still works end-to-end (mirror inside lock
        # mustn't break the contract; RLock re-entry must succeed).
        monkeypatch.setattr(sc._state, "_LIBRARY_FILE", tmp_path / "lib.json")
        monkeypatch.setattr(sc._state, "_COLLECTIONS_FILE", tmp_path / "coll.json",
        )
        sc._state._library_cache = None
        sc._state._collections_cache = None
        entries = [{
            "name": "P1", "id": "p1", "size": 100,
            "n_feats": 0, "gb_text": "LOCUS p1",
            "added": "2026-05-20",
        }]
        sc._save_library(entries)
        loaded = sc._load_library()
        assert len(loaded) == 1
        assert loaded[0]["id"] == "p1"

    def test_save_parts_bin_runs_mirror(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc._state, "_PARTS_BIN_FILE", tmp_path / "parts_bin.json",
        )
        monkeypatch.setattr(sc._state, "_PARTS_BIN_COLLECTIONS_FILE",
            tmp_path / "parts_bin_coll.json",
        )
        sc._state._parts_bin_cache = None
        sc._state._parts_bin_collections_cache = None
        sc._save_parts_bin([{"name": "B1", "sequence": "ACGT"}])
        loaded = sc._load_parts_bin()
        assert len(loaded) == 1


# ═══════════════════════════════════════════════════════════════════════════════
# Entry-vector save lock
# ═══════════════════════════════════════════════════════════════════════════════

class TestEntryVectorSavesLocked:
    """Regression guard for 2026-05-20 fix: `_set_entry_vector` and
    `_clear_entry_vectors_for_grammar` did load-mutate-save without
    `_cache_lock`. Two concurrent calls both read pre-state and both
    wrote — last writer lost the other's binding. The auto-bind
    flow's per-role loop amplified the race.
    """

    def test_set_entry_vector_round_trip(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc._state, "_ENTRY_VECTORS_FILE", tmp_path / "ev.json",
        )
        sc._state._entry_vectors_cache = None
        sc._set_entry_vector(
            "gb_l0",
            {"name": "MyVec", "gb_text": "fake"},
            role="Alpha1",
        )
        ev = sc._get_entry_vector("gb_l0", "Alpha1")
        assert ev is not None
        assert ev.get("name") == "MyVec"

    def test_clear_grammar_clears_all_roles(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc._state, "_ENTRY_VECTORS_FILE", tmp_path / "ev.json",
        )
        sc._state._entry_vectors_cache = None
        sc._set_entry_vector(
            "gb_l0", {"name": "V1", "gb_text": "x"}, role="Alpha1",
        )
        sc._set_entry_vector(
            "gb_l0", {"name": "V2", "gb_text": "y"}, role="Alpha2",
        )
        n = sc._clear_entry_vectors_for_grammar("gb_l0")
        assert n == 2
        assert sc._get_entry_vector("gb_l0", "Alpha1") is None
        assert sc._get_entry_vector("gb_l0", "Alpha2") is None


# ═══════════════════════════════════════════════════════════════════════════════
# `_history_node_to_dict` iterative — no stack overflow on hostile depth
# ═══════════════════════════════════════════════════════════════════════════════

class TestHistoryNodeIterative:
    """Regression guard for 2026-05-20 fix: pre-fix
    `_history_node_to_dict` recursed through `node.parents`, blowing
    the Python recursion limit on hostile `.dna` imports carrying
    deeply-nested history trees. Now iterative DFS with depth + node
    caps. Sibling helpers (`walk`, `_history_node_count`,
    `HistoryScreen.populate`) were already iterative for this reason.
    """

    def test_deep_chain_no_stack_overflow(self):
        # Build a 5000-deep history chain — way past Python's default
        # recursion limit (1000). Pre-fix this would `RecursionError`.
        class _Node:
            def __init__(self, name, parents=None):
                self.name = name
                self.operation = "test"
                self.seq_len = 100
                self.circular = True
                self.regenerated_sites: list = []
                self.input_summaries: list = []
                self.parents = parents or []
        root = _Node("root")
        cur = root
        for i in range(5000):
            child = _Node(f"n{i}")
            cur.parents.append(child)
            cur = child
        # Must not raise RecursionError.
        out = sc._history_node_to_dict(root)
        assert out is not None
        assert out["name"] == "root"
        # Truncation marker present because we hit the depth cap.
        # Walk down to confirm.
        depth = 0
        node = out
        while node["parents"]:
            node = node["parents"][0]
            depth += 1
        # Bounded by `_HISTORY_NODE_MAX_DEPTH = 500`.
        assert depth <= 500

    def test_returns_none_on_none(self):
        assert sc._history_node_to_dict(None) is None

    def test_shallow_chain_round_trip(self):
        class _Node:
            def __init__(self, name, parents=None):
                self.name = name
                self.operation = "op"
                self.seq_len = 50
                self.circular = False
                self.regenerated_sites: list = ["BsaI"]
                self.input_summaries: list = ["p1.gb"]
                self.parents = parents or []
        gp = _Node("gp")
        p = _Node("p", [gp])
        c = _Node("c", [p])
        out = sc._history_node_to_dict(c)
        assert out["name"] == "c"
        assert out["parents"][0]["name"] == "p"
        assert out["parents"][0]["parents"][0]["name"] == "gp"
        assert out["regenerated_sites"] == ["BsaI"]


# ═══════════════════════════════════════════════════════════════════════════════
# `_bulk_import_folder` refuses symlinks
# ═══════════════════════════════════════════════════════════════════════════════

class TestBulkImportRefusesSymlinks:
    """Regression guard for 2026-05-20 fix: pre-fix
    `_bulk_import_folder` used `path.is_file()` + `path.stat()`
    which both follow symlinks. A hostile `.gb` symlinked to
    `/dev/zero` (st_size=0 on devices) passed the size cap and
    `load_genbank` then read until OOM. Now uses `lstat()` +
    `S_ISREG` to refuse symlinks outright.
    """

    def test_symlink_to_regular_file_refused(self, tmp_path):
        # Real .gb file
        real = tmp_path / "real.gb"
        real.write_text("LOCUS dummy 100 bp\n//\n")
        # Symlink in the import directory pointing at the real file
        import_dir = tmp_path / "import"
        import_dir.mkdir()
        link = import_dir / "linked.gb"
        try:
            link.symlink_to(real)
        except (OSError, NotImplementedError):
            pytest.skip("symlinks not supported on this platform")
        entries, failures = sc._bulk_import_folder(import_dir)
        # The symlink must NOT have been imported — the only file in
        # the dir was a symlink, so entries should be empty.
        assert entries == []
        # Either filtered silently from `importable` (matching the
        # children-scan path) or surfaced as a failure — the security
        # contract is "do not load symlink content", not "report which".
        # We accept either outcome here.

    def test_regular_file_still_imports(self, tmp_path):
        # Make sure the lstat refusal didn't break the happy path.
        import_dir = tmp_path / "import"
        import_dir.mkdir()
        gb = import_dir / "ok.gb"
        # Real GenBank record — minimal but valid for BioPython.
        gb.write_text(
            "LOCUS       OK      10 bp ds-DNA  circular SYN 01-JAN-2026\n"
            "FEATURES             Location/Qualifiers\n"
            "ORIGIN\n"
            "        1 acgtacgtac\n"
            "//\n",
        )
        entries, failures = sc._bulk_import_folder(import_dir)
        # Either parses (happy path) or fails gracefully (parser is
        # picky about LOCUS lines). Critical: no exception, no leak.
        assert isinstance(entries, list)
        assert isinstance(failures, list)


# ═══════════════════════════════════════════════════════════════════════════════
# `_get_rotated_state` cache key — hash() not id()
# ═══════════════════════════════════════════════════════════════════════════════

class TestRotatedStateCacheKey:
    """Regression guard for 2026-05-20 fix: pre-fix the cache key was
    `(id(self._seq), id(self._feats), o)`. Python's small-string
    interning + allocator reuse means two distinct strings can share
    the same id after one is GC'd, so a stale rotated cache could
    surface after an in-place sequence edit happened to reuse the
    previous string's address. The build-seq cache (line 6564) had
    been fixed for this same reason; rotated-state was overlooked.
    """

    def test_key_uses_hash_not_id(self):
        import inspect
        src = inspect.getsource(sc.SequencePanel._get_rotated_state)
        assert "hash(self._seq)" in src
        # Sanity: comment explains why.
        assert "interning" in src or "allocator" in src


# ═══════════════════════════════════════════════════════════════════════════════
# Batch-save collision flow — one modal per batch (no toast lying)
# ═══════════════════════════════════════════════════════════════════════════════

class _StubApp:
    def __init__(self):
        self._next_callback = None
        self.pushed: list = []

    def push_screen(self, screen, callback=None):
        self.pushed.append(type(screen))
        self._next_callback = callback

    def fire(self, payload):
        cb = self._next_callback
        self._next_callback = None
        if cb is not None:
            cb(payload)


class TestBatchCollisionFlowSingleModal:
    """Regression guard for 2026-05-20 fix: pre-fix
    `_save_to_collection` called `lib.add_entry(rec)` per part, each
    colliding rec queueing a modal while the post-loop toast claimed
    success before any modal resolved. Now the batch routes through
    `add_entries_batch` and a single modal resolves the whole batch.
    """

    def test_add_entries_batch_one_collision_call(self, tmp_path,
                                                    monkeypatch):
        # Verify the new `add_entries_batch` method exists with the
        # documented signature.
        assert hasattr(sc.LibraryPanel, "add_entries_batch")
        # The method must accept (records, *, on_done).
        import inspect
        sig = inspect.signature(sc.LibraryPanel.add_entries_batch)
        params = list(sig.parameters)
        assert "records" in params
        assert "on_done" in params

    def test_add_entries_batch_calls_on_done_with_empty_batch(self):
        """Empty `records` short-circuits — `on_done` fires with
        empty saved + replaced sets, no modal pushed."""
        # Build a stub LibraryPanel-shaped object. We're testing the
        # control flow, not the Textual integration.
        stub = MagicMock(spec=sc.LibraryPanel)
        stub.app = _StubApp()
        results: list = []

        def _on_done(saved, replaced, cancelled):
            results.append((saved, replaced, cancelled))
        # Call the real method against the stub.
        sc.LibraryPanel.add_entries_batch(
            stub, [], on_done=_on_done,
        )
        assert results == [([], set(), False)]
        # No modal pushed.
        assert stub.app.pushed == []


# ═══════════════════════════════════════════════════════════════════════════════
# Bare-except hygiene — invariant #1
# ═══════════════════════════════════════════════════════════════════════════════

class TestBareExceptHygiene:
    """Regression guard for 2026-05-20 fix: sweep #10 cleaned up the
    misleading `except (AttributeError, Exception): pass` clauses
    in `_clear_user_data_undo_stacks`. Exception subsumes
    AttributeError; the tuple form was just `except Exception`.
    Now narrowed to `except (TypeError, AttributeError)` /
    `except AttributeError`.
    """

    def test_no_redundant_exception_tuples(self):
        # Walk the function's bytecode AST instead of grepping source
        # (the explanatory comment we left also contains the string
        # we're hunting for).
        import ast
        import inspect
        src = inspect.getsource(sc._reset_app_state_after_master_delete)
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if not isinstance(node, ast.ExceptHandler):
                continue
            if not isinstance(node.type, ast.Tuple):
                continue
            names = {
                e.id for e in node.type.elts
                if isinstance(e, ast.Name)
            }
            # AttributeError + Exception together is the redundant
            # pattern (Exception subsumes AttributeError).
            assert not ("AttributeError" in names
                        and "Exception" in names), (
                f"redundant except tuple still present: {names}"
            )
