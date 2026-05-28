"""
test_sweep11 — adversarial sweep #11 regression coverage (2026-05-20).

Sweep #11 is the deferred-items pass following sweep #10. Targets:

  * Active-pointer flip race at collection / parts-bin switch +
    delete-promote (4 sites).
  * Cross-group worker library writes — RMW under `_cache_lock`
    so Constructor / Traditional / Gibson / Domesticator-mirror
    saves don't silently overwrite each other.
  * Settings flush worker drain on Master Delete — daemon worker
    can no longer resurrect wiped settings.
  * UI freeze on bulk-import auto-bind — moved to @work + batched
    entry-vector writes.
  * `_find_library_entry_by_id` hot-path helper.
  * Agent API hardening: `_h_create_collection` atomic check,
    `_h_hmmscan` filesystem-state oracle collapse,
    `.dna` writer XML control-char strip.

See CLAUDE.md invariant #51 for the full inventory.
"""
from __future__ import annotations

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Active-pointer flip race — settings flush sync before file write
# ═══════════════════════════════════════════════════════════════════════════════

class TestActivePointerFlushSync:
    """Regression guard: every active-pointer flip that's followed by
    a sync file write must call `_settings_flush_sync()` between
    them. Pre-fix a power loss could leave settings.json saying OLD
    while the dependent file held NEW.
    """

    def test_collection_switch_calls_flush(self):
        import inspect
        # LibraryPanel.action_open_collection (or _on_select) →
        # _set_active_collection_name + _settings_flush_sync.
        # Look at the collection-switch click handler.
        src = inspect.getsource(sc.LibraryPanel)
        # The pattern "_set_active_collection_name(name)" followed
        # within a few lines by "_settings_flush_sync()" must exist.
        # Verify the literal pairing is present somewhere.
        assert "_set_active_collection_name(name)" in src
        assert "_settings_flush_sync()" in src
        # And the collection delete-promote path uses the same pair.
        # Pre-sweep #11 the delete-promote path used `_set_active_*`
        # alone with no flush.
        # White-box: count flush sites under the panel — should be
        # at least 2 (switch + promote).
        n_flush = src.count("_settings_flush_sync()")
        assert n_flush >= 2, (
            f"expected ≥2 flush_sync calls in LibraryPanel, got {n_flush}"
        )

    def test_parts_bin_switch_calls_flush(self):
        import inspect
        # PartsBinPickerModal._open + _delete pair must both flush.
        src = inspect.getsource(sc.PartsBinPickerModal)
        assert "_set_active_parts_bin_name(name)" in src
        assert "_set_active_parts_bin_name(promoted)" in src
        # Both sites must follow with flush_sync.
        assert "_settings_flush_sync()" in src
        # At least 2 flush sites (open + delete-promote).
        n_flush = src.count("_settings_flush_sync()")
        assert n_flush >= 2, (
            f"expected ≥2 flush_sync calls in PartsBinPickerModal, got {n_flush}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Cross-group worker library writes — RMW under _cache_lock
# ═══════════════════════════════════════════════════════════════════════════════

class TestCrossGroupWriteRMW:
    """Regression guard: every worker that does a library
    load-modify-save must wrap the sequence under `_cache_lock`.
    Pre-fix workers in different `@work` groups bypassed each
    other's exclusivity and the second writer silently dropped
    the first writer's new entry.
    """

    def test_constructor_persist_holds_lock(self):
        import inspect
        # Look at _persist_assembly. The save block runs under
        # `with _cache_lock:`.
        for name, obj in vars(sc).items():
            if not isinstance(obj, type):
                continue
            if not hasattr(obj, "_persist_assembly"):
                continue
            src = inspect.getsource(obj._persist_assembly)
            assert "with _cache_lock" in src, (
                f"{name}._persist_assembly missing _cache_lock"
            )
            assert "_save_library" in src

    def test_gibson_save_holds_lock(self):
        import inspect
        # Find _gibson_save_worker and verify _cache_lock is present.
        for name, obj in vars(sc).items():
            if not isinstance(obj, type):
                continue
            if not hasattr(obj, "_gibson_save_worker"):
                continue
            src = inspect.getsource(obj._gibson_save_worker)
            # The cross-group RMW race protection now lives in the shared
            # `_commit_library_entry_to_collection` helper (which holds
            # `_cache_lock`); accept either the inline lock or delegation
            # to that locked helper.
            assert ("with _cache_lock" in src
                    or "_commit_library_entry_to_collection" in src), (
                f"{name}._gibson_save_worker missing lock-protected commit"
            )

    def test_domesticator_mirror_holds_lock(self):
        import inspect
        for name, obj in vars(sc).items():
            if not isinstance(obj, type):
                continue
            if not hasattr(obj, "_domesticator_library_mirror_worker"):
                continue
            src = inspect.getsource(
                obj._domesticator_library_mirror_worker
            )
            assert "with _cache_lock" in src, (
                f"{name}._domesticator_library_mirror_worker "
                f"missing _cache_lock"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# Settings flush worker drain on Master Delete
# ═══════════════════════════════════════════════════════════════════════════════

class TestMasterDeleteFlushDrain:
    """Regression guard: `_perform_master_delete` must wait for an
    in-flight settings-flush worker to drain before unlinking files.
    Pre-fix the worker could resurrect wiped settings if it was
    mid-`_safe_save_json` when the wipe started.
    """

    def test_master_delete_drains_flush(self):
        import inspect
        src = inspect.getsource(sc._perform_master_delete)
        # The drain handshake polls `_settings_flush_running` under
        # the flush lock with a bounded deadline.
        assert "_settings_flush_running" in src
        assert "_settings_flush_lock" in src
        # The drain uses a deadline (not infinite wait).
        assert "deadline" in src or "timeout" in src.lower()


# ═══════════════════════════════════════════════════════════════════════════════
# UI freeze: auto-bind moved to @work
# ═══════════════════════════════════════════════════════════════════════════════

class TestAutoBindMovedToWorker:
    """Regression guard: bulk-import auto-bind must dispatch to a
    `@work` thread, not run inline on the UI thread.
    """

    def test_library_panel_has_auto_bind_worker(self):
        assert hasattr(sc.LibraryPanel, "_auto_bind_worker"), (
            "LibraryPanel._auto_bind_worker missing — pre-sweep #11 "
            "the auto-bind ran inline on UI thread"
        )

    def test_set_entry_vectors_batch_exists(self):
        # Batch-write helper exists.
        assert callable(getattr(sc, "_set_entry_vectors_batch", None))


class TestSetEntryVectorsBatchSemantics:
    """The new batch helper must take `_cache_lock` once + write
    once, AND honour the same per-tuple contract as
    `_set_entry_vector`.
    """

    def test_batch_round_trip(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            sc, "_ENTRY_VECTORS_FILE", tmp_path / "ev.json",
        )
        sc._entry_vectors_cache = None
        # 3 updates, one save call.
        updates = [
            ("gb_l0", "Alpha1",
             {"name": "V1", "gb_text": "x"}),
            ("gb_l0", "Alpha2",
             {"name": "V2", "gb_text": "y"}),
            ("gb_l0", "Omega1",
             {"name": "V3", "gb_text": "z"}),
        ]
        n_changed = sc._set_entry_vectors_batch(updates)
        assert n_changed == 3
        assert (sc._get_entry_vector("gb_l0", "Alpha1") or {}).get("name") == "V1"
        assert (sc._get_entry_vector("gb_l0", "Alpha2") or {}).get("name") == "V2"
        assert (sc._get_entry_vector("gb_l0", "Omega1") or {}).get("name") == "V3"

    def test_batch_empty_is_noop(self):
        # No-op for empty input.
        assert sc._set_entry_vectors_batch([]) == 0

    def test_batch_skip_invalid_grammar_id(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            sc, "_ENTRY_VECTORS_FILE", tmp_path / "ev.json",
        )
        sc._entry_vectors_cache = None
        updates = [
            ("", "Alpha1", {"name": "Bad", "gb_text": "x"}),
            ("gb_l0", "Alpha1", {"name": "Good", "gb_text": "y"}),
        ]
        n_changed = sc._set_entry_vectors_batch(updates)
        assert n_changed == 1
        ev = sc._get_entry_vector("gb_l0", "Alpha1")
        assert ev is not None
        assert ev.get("name") == "Good"


# ═══════════════════════════════════════════════════════════════════════════════
# `_find_library_entry_by_id` helper
# ═══════════════════════════════════════════════════════════════════════════════

class TestFindLibraryEntryByIdHelper:
    """Regression guard: the helper exists, returns deep-copies, and
    correctly handles edge cases (empty id, non-string id, no match).
    """

    def test_helper_exists(self):
        assert callable(getattr(sc, "_find_library_entry_by_id", None))

    def test_empty_id_returns_none(self):
        assert sc._find_library_entry_by_id("") is None
        # Non-string defensively returns None too.
        assert sc._find_library_entry_by_id(None) is None  # type: ignore[arg-type]
        assert sc._find_library_entry_by_id(42) is None  # type: ignore[arg-type]

    def test_no_match_returns_none(self, tmp_path, monkeypatch):
        monkeypatch.setattr(sc, "_LIBRARY_FILE", tmp_path / "lib.json")
        sc._library_cache = None
        # Empty library
        sc._save_library([])
        assert sc._find_library_entry_by_id("nonexistent") is None

    def test_match_returns_deep_clone(self, tmp_path, monkeypatch):
        # Test data uses `id == sanitize(name)` so the post-2026-05-24
        # id-name backfill (PIT-36) is a no-op and lookup by id stays
        # stable across the load.
        monkeypatch.setattr(sc, "_LIBRARY_FILE", tmp_path / "lib.json")
        sc._library_cache = None
        sc._save_library([
            {"id": "P1", "name": "P1", "gb_text": "LOCUS p1"},
            {"id": "P2", "name": "P2", "gb_text": "LOCUS p2"},
        ])
        e = sc._find_library_entry_by_id("P1")
        assert e is not None
        assert e["name"] == "P1"
        # Mutating the returned dict must NOT poison the cache.
        e["name"] = "MUTATED"
        e2 = sc._find_library_entry_by_id("P1")
        assert e2 is not None
        assert e2["name"] == "P1"


# ═══════════════════════════════════════════════════════════════════════════════
# Agent API hardening
# ═══════════════════════════════════════════════════════════════════════════════

class TestAgentApiHardening:
    """Regression guards for the agent-API MEDIUM fixes.
    """

    def test_h_create_collection_atomic_check(self):
        import inspect
        src = inspect.getsource(sc._h_create_collection)
        # Two name-taken checks now: the early-out + the atomic
        # re-check inside the cache lock.
        assert src.count("_collection_name_taken(name)") >= 2
        assert "_cache_lock" in src

    def test_h_hmmscan_collapses_error_oracle(self):
        import inspect
        src = inspect.getsource(sc._h_hmmscan)
        # The 404 / 400 differential is collapsed to a single
        # generic 400 rejection.
        assert "not acceptable" in src
        # No more differentiated "not found" vs "rejected" string
        # in the response body — the detail goes to the log only.
        # Verify the log path is present.
        assert "_log.info" in src

    def test_dna_writer_strips_control_chars(self):
        import inspect
        src = inspect.getsource(
            sc._build_commercialsaas_features_packet_from_record
        )
        # Sanitisation pass on qualifier `<V text=>` values.
        assert "sanitised" in src or "sanitized" in src
        # The strip rule keeps printable + tab + newline.
        assert "\\t" in src or "\\n" in src
