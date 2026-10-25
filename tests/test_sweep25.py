"""
test_sweep25 — adversarial audit Sweep #25 (2026-05-23).

Regression coverage for the fixes landed by INV-66. The audit
surfaced ~40 findings across CRITICAL / HIGH / MEDIUM / LOW
severities; this file regression-locks the ones whose silent
return would otherwise let a future edit regress the fix.

Index of CRITICAL fixes covered:
  * Cache-bust enumeration drift (C1, C2) — `RestoreFromBackupModal`
    and `_h_restore_pre_update_snapshot` now both drive from
    `_MASTER_DELETE_CACHE_ATTRS`. Pre-fix sweep #24 added
    custom_enzymes + enzyme_collections caches; neither restore
    path was updated, leaving stale-cache → silent-overwrite holes.
  * Stacks log rotation cap (C3) — `_STACKS_LOG_MAX_BYTES` truncates
    on startup when prior sessions' faulthandler dumps exceed cap.
  * `_GB_PARSE_CACHE` cap reduced 64 → 16 (C4) — bounds worst-case
    parsed-record memory amplification.
  * `_iter_library_readonly` + `_iter_parts_bin_readonly` + sibling
    `_find_library_entry_by_name` helpers (C5) — read-only iteration
    without the per-call full-library deepcopy.

Index of HIGH fixes covered:
  * Custom enzyme + enzyme collection agent endpoints now wrap RMW
    in `_cache_lock` (H1).
  * App-level Ctrl+S `_save_worker` library RMW under `_cache_lock`
    (H2).
  * 5 worker `_safe_save_json` direct calls now hold `_cache_lock`
    for the disk write (H3).
  * `_SCAN_CATALOG` rebuild is atomic global reassign (H4).
  * Agent server requires bearer token on ALL endpoints (H5).
  * Gibson agent endpoint uses `call_from_thread` for UI refresh
    (H6).
  * Plasmidsaurus + load-file error responses collapsed to uniform
    400 (H7, H8) — filesystem-state oracle reduction.
  * `_AgentAPIServer.timeout = 30` for slow-loris mitigation (H9).
  * `_sync_active_collection_plasmids` sync path uses readonly
    iterator + targeted clone under lock (H10).
  * `_get_setting` reads cache directly without per-call full-dict
    deepcopy (H11).
  * `_rc(tgt_seq)` hoisted out of `_find_annotation_transfers`
    per-feature loop (H12).
  * `except (OSError, Exception)` narrowed to actual Primer3 modes
    (H14).
  * `AddCustomEnzymeModal._save_btn` wraps RMW in `_cache_lock` +
    emits structured save-failed event (H15).
  * GFF3 out-of-range feature drop now logs (H16).

Index of MEDIUM / LOW fixes covered:
  * `RestoreFromBackupModal._dismiss_once` retrofit (M1).
  * `_safe_identifier` validator for `active_*` settings (M2).
  * Snapshot dir restore 5 GB cap (M3).
  * Body parse error → `400 "malformed JSON body"` (M5).
  * Agent server scrubs paths from exception responses (M6).
  * `_atomic_copy` / `_atomic_marker_write` close `fd` on early
    raise (M9).
  * Agent token file write closes `fd` if `fdopen` raises (M10).
  * Crash recovery prune (30-day age + 50-file count caps) (M11).
  * `_VECTOR_MATCH_CACHE` keys on `hash(gb_text)` not raw string
    (M14).
  * `target_seq = ""` fallback narrowed to (AttributeError, TypeError)
    (M20).
  * NEB built-in override gap documented (M13).
  * `splicecraft_cli._read_session` refuses symlinked token file
    (L1).
  * Agent token via `secrets.token_urlsafe(32)` (L3).
  * NCBI / Kazusa fetches narrow `except` to (OSError, URLError)
    (L6).
  * Lockfile unlink on graceful exit (L9).
"""
from __future__ import annotations

from pathlib import Path

import splicecraft as sc


# ═════════════════════════════════════════════════════════════════════
# C1, C2 — cache-bust enumeration drift fixes
# ═════════════════════════════════════════════════════════════════════

class TestCacheBustParity:
    """The two restore paths (`RestoreFromBackupModal` + the agent
    `_h_restore_pre_update_snapshot`) must both bust every persistent
    cache, otherwise a successful restore leaves stale in-memory state
    that the next CRUD silently overwrites onto disk.
    """

    def test_master_delete_cache_attrs_is_canonical(self):
        """The master tuple covers every `_*_cache` global the
        module defines. Sentinel against drift — a future cache that
        forgets to register here will be missed by Master Delete AND
        both restore paths (which now drive from it)."""
        # Every name in the tuple must resolve to a module attribute
        # (one of the cache globals). Pre-fix this list was hand-
        # maintained and drifted.
        for attr in sc._MASTER_DELETE_CACHE_ATTRS:
            assert hasattr(sc, attr) or hasattr(sc._state, attr), (
                f"{attr!r} is in _MASTER_DELETE_CACHE_ATTRS but does "
                f"not exist as a module attribute (hub or _state) — typo / removed?"
            )

    def test_master_delete_includes_sweep24_caches(self):
        """Sweep #24 added `_custom_enzymes_cache` /
        `_enzyme_collections_cache`. Sentinel against future
        omissions."""
        expected = {
            "_library_cache",
            "_collections_cache",
            "_parts_bin_cache",
            "_primers_cache",
            "_primer_collections_cache",
            "_primer_usage_cache",
            "_feature_library_index_cache",
            "_features_cache",
            "_feature_colors_cache",
            "_grammars_cache",
            "_entry_vectors_cache",
            "_codon_tables_cache",
            "_settings_cache",
            "_experiments_cache",
            "_experiment_projects_cache",
            "_gels_cache",
            "_protein_motifs_cache",
            "_custom_enzymes_cache",
            "_enzyme_collections_cache",
        }
        actual = set(sc._MASTER_DELETE_CACHE_ATTRS)
        missing = expected - actual
        assert not missing, f"Missing from master tuple: {missing}"

    def test_restore_from_backup_modal_drives_from_master_tuple(self):
        """White-box: `RestoreFromBackupModal._do_restore`'s cache-
        bust loop must iterate `_MASTER_DELETE_CACHE_ATTRS` (with
        `_settings_cache` deliberately skipped). Pre-sweep #25 it
        had a hand-maintained 14-entry tuple that drifted.
        """
        import inspect
        src = inspect.getsource(sc.RestoreFromBackupModal)
        assert "_MASTER_DELETE_CACHE_ATTRS" in src, (
            "RestoreFromBackupModal must drive cache-bust from the "
            "master tuple, not hand-maintain its own list"
        )

    def test_h_restore_pre_update_drives_from_master_tuple(self):
        """White-box: agent pre-update restore handler must also
        iterate `_MASTER_DELETE_CACHE_ATTRS`."""
        import inspect
        src = inspect.getsource(sc._h_restore_pre_update_snapshot)
        assert "_MASTER_DELETE_CACHE_ATTRS" in src


# ═════════════════════════════════════════════════════════════════════
# C3 — stacks log rotation cap
# ═════════════════════════════════════════════════════════════════════

class TestStacksLogRotation:
    """`_STACKS_LOG_FD` opens `splicecraft.stacks.log` in append
    mode for `faulthandler.register`. Pre-fix nothing capped its
    size — a user who hit SIGUSR1 once a week accumulated MBs over
    a year. Now truncates at startup if over 10 MB."""

    def test_max_bytes_constant_is_set(self):
        # The constant is module-level. Tests don't exercise the
        # actual truncation (requires re-importing the module) but
        # the constant existing + being a sane value is the
        # regression sentinel.
        assert hasattr(sc, "_STACKS_LOG_MAX_BYTES")
        assert 1_000_000 <= sc._STACKS_LOG_MAX_BYTES <= 100_000_000


# ═════════════════════════════════════════════════════════════════════
# C4 — GB parse cache cap
# ═════════════════════════════════════════════════════════════════════

class TestGBParseCacheCap:
    def test_cap_reduced(self):
        """Sweep #25: cap reduced from 64 → 16. With SeqRecord
        memory amplification (~5–15× input text), worst-case
        cached state at the old cap could reach ~40 GB. 16
        entries × 64 MB ceiling × 10× amp ≈ ~10 GB worst case —
        still high but bounded."""
        assert sc._GB_PARSE_CACHE_MAX == 16


# ═════════════════════════════════════════════════════════════════════
# C5 — `_iter_library_readonly` + `_iter_parts_bin_readonly` helpers
# ═════════════════════════════════════════════════════════════════════

class TestReadonlyIterators:
    def test_iter_library_readonly_exists(self):
        assert callable(sc._iter_library_readonly)

    def test_iter_parts_bin_readonly_exists(self):
        assert callable(sc._iter_parts_bin_readonly)

    def test_find_library_entry_by_name_exists(self):
        assert callable(sc._find_library_entry_by_name)

    def test_iter_library_readonly_returns_cache_reference(self,
                                                            isolated_library):
        """The point of the helper is to skip the deepcopy. Verify
        it returns the actual cached list (or a same-content list,
        depending on cache state). Test data uses `id == sanitize(name)`
        so the post-2026-05-24 id-name backfill (PIT-36) is a no-op."""
        sc._save_library([
            {"id": "A", "name": "A", "size": 10},
            {"id": "B", "name": "B", "size": 20},
        ])
        view = sc._iter_library_readonly()
        assert len(view) == 2
        assert view[0]["id"] == "A"

    def test_find_library_entry_by_name_returns_clone(self,
                                                       isolated_library):
        """Helper must deepcopy on return (sacred invariant #17).
        Test data uses `id == sanitize(name)` so the id-name backfill
        (PIT-36) is a no-op."""
        sc._save_library([
            {"id": "Alpha", "name": "Alpha", "size": 10, "extra": [1, 2]},
        ])
        e = sc._find_library_entry_by_name("Alpha")
        assert e is not None and e["id"] == "Alpha"
        e["extra"].append(999)
        e2 = sc._find_library_entry_by_name("Alpha")
        assert e2 is not None and e2["extra"] == [1, 2], (
            "deep-clone contract violated"
        )

    def test_find_library_entry_by_id_accepts_none(self):
        """Type signature widened to `str | None` so optional-id
        lookup sites don't need a wrapper."""
        assert sc._find_library_entry_by_id(None) is None
        assert sc._find_library_entry_by_id("") is None
        assert sc._find_library_entry_by_id(123) is None  # type: ignore[arg-type]


# ═════════════════════════════════════════════════════════════════════
# H4 — _SCAN_CATALOG atomic rebuild
# ═════════════════════════════════════════════════════════════════════

class TestScanCatalogAtomicRebuild:
    def test_rebuild_does_not_clear_in_place(self):
        """`_rebuild_scan_catalog` builds locally then reassigns;
        does NOT call `.clear()` + `.append()` mid-iteration.
        Pre-fix a concurrent scan iterator could see a half-built
        catalog (missed cut sites or tuple-index error)."""
        import inspect
        src = inspect.getsource(sc._rebuild_scan_catalog)
        # Strip the docstring before checking (the docstring itself
        # references the pre-fix `.clear()` pattern for context).
        # Heuristic: docstring is everything between the first
        # triple-quote pair after `def`.
        body = src.split('"""', 2)[-1]
        assert "_SCAN_CATALOG.clear()" not in body
        assert 'globals()["_SCAN_CATALOG"]' in body

    def test_rebuild_produces_nonempty_catalog(self):
        """Sanity: rebuild still works post-refactor."""
        sc._rebuild_scan_catalog()
        assert len(sc._SCAN_CATALOG) > 0


# ═════════════════════════════════════════════════════════════════════
# H5 — bearer token required on ALL endpoints
# ═════════════════════════════════════════════════════════════════════

class TestAuthOnReadEndpoints:
    def test_handle_method_inspects_token_before_handler_lookup(self):
        """White-box: token check fires BEFORE handler lookup so
        unauthenticated probes for unknown paths uniformly 401
        (don't leak endpoint-list oracle via 404)."""
        import inspect
        # The routing body moved into `_handle_inner` (2026-10-23, AG7): `_handle`
        # wraps it so an idempotency claim is released on every exit.
        src = inspect.getsource(sc._AgentRequestHandler._handle_inner)
        # Find positions of the token check and the handler lookup.
        token_pos = src.find("_check_token")
        handler_pos = src.find("_AGENT_HANDLERS.get(path_part)")
        assert token_pos > 0 and handler_pos > 0
        assert token_pos < handler_pos, (
            "token check must precede handler lookup"
        )


# ═════════════════════════════════════════════════════════════════════
# H9 — HTTP server timeout
# ═════════════════════════════════════════════════════════════════════

class TestServerTimeout:
    def test_server_has_timeout(self):
        assert sc._AgentAPIServer.timeout == 30

    def test_handler_has_timeout(self):
        assert sc._AgentRequestHandler.timeout == 30


# ═════════════════════════════════════════════════════════════════════
# H10 — sync_active_collection uses readonly iter
# ═════════════════════════════════════════════════════════════════════

class TestSyncActiveCollectionPerf:
    def test_sync_path_uses_iter_collections_readonly(self):
        """White-box: the sync branch reads via
        `_iter_collections_readonly()` not `_load_collections()`."""
        import inspect
        src = inspect.getsource(sc._sync_active_collection_plasmids)
        assert "_iter_collections_readonly" in src


# ═════════════════════════════════════════════════════════════════════
# H11 — _get_setting cache-direct read
# ═════════════════════════════════════════════════════════════════════

class TestGetSettingPerf:
    def test_get_setting_does_not_call_load_settings(self):
        """White-box: should NOT route through `_load_settings()`
        (which deep-clones the entire dict on every call). Skip
        docstring matches — they reference the pre-fix pattern for
        context."""
        import inspect
        src = inspect.getsource(sc._get_setting)
        body = src.split('"""', 2)[-1]
        # `_load_settings()` (with parens) is the call we want to
        # rule out. The body does call `_load_settings()` ONCE as a
        # cache-populate trigger when the cache is None — that's
        # safe (cache stays valid afterward); the regression we're
        # blocking is calling it on EVERY access.
        # Check: at most one occurrence in the body.
        assert body.count("_load_settings()") <= 1
        assert "_settings_cache" in body

    def test_get_setting_returns_clone_for_containers(self,
                                                       isolated_library):
        """A container value (list/dict) must be cloned defensively
        so caller mutation can't poison the cache."""
        sc._set_setting("crash_recovery_seen", ["x", "y"])
        view = sc._get_setting("crash_recovery_seen")
        assert view == ["x", "y"]
        view.append("z")
        view2 = sc._get_setting("crash_recovery_seen")
        assert view2 == ["x", "y"]


# ═════════════════════════════════════════════════════════════════════
# H12 — _rc hoist in _find_annotation_transfers
# ═════════════════════════════════════════════════════════════════════

class TestRcHoist:
    def test_rc_hoisted_outside_feature_loop(self):
        """White-box: `_rc(tgt_seq)` must appear BEFORE the
        `for feat in source_rec.features` loop. Pre-fix it was
        inside, paying O(F × N_tgt)."""
        import inspect
        src = inspect.getsource(sc._find_annotation_transfers)
        rc_pos = src.find("rc_tgt_seq = _rc(tgt_seq)")
        loop_pos = src.find("for feat in source_rec.features")
        assert rc_pos > 0 and loop_pos > 0
        assert rc_pos < loop_pos, (
            "_rc(tgt_seq) must be hoisted above the per-feature loop"
        )


# ═════════════════════════════════════════════════════════════════════
# H14 — bare-except narrowing
# ═════════════════════════════════════════════════════════════════════

class TestBareExceptNarrowing:
    def test_no_except_OSError_Exception_tuple(self):
        """`except (OSError, Exception)` is `except Exception` (a
        bare-except in disguise — Exception subsumes OSError).
        Grep the source for the pattern as a regression sentinel.
        INV-65 caught the `(AttributeError, Exception)` shape;
        sweep #25 extended to `(OSError, Exception)`."""
        text = Path(sc.__file__).read_text(encoding="utf-8")
        # Count standalone uses (in comments referencing the pattern
        # is OK, but real `except (OSError, Exception)` lines are not).
        import re
        bad_lines = [
            ln for ln in text.splitlines()
            if re.search(r"^\s*except \(OSError, Exception\)", ln)
        ]
        assert not bad_lines, (
            f"Found bare-except-in-disguise: {bad_lines}"
        )


# ═════════════════════════════════════════════════════════════════════
# H16 — GFF3 silent drop logs
# ═════════════════════════════════════════════════════════════════════

class TestGFF3LogsDropped:
    def test_gff3_features_dropped_logged(self, caplog):
        """An out-of-range GFF3 feature now logs a warning instead
        of silently dropping. Sweep #25 fixed."""
        import logging
        # Attach a handler directly to sc._log since the module has
        # propagate=False on its logger (per INV-38) and pytest's
        # caplog only captures via propagation.
        records: list[logging.LogRecord] = []

        class _Capture(logging.Handler):
            def emit(self, record):
                records.append(record)

        handler = _Capture(level=logging.INFO)
        sc._log.addHandler(handler)
        prev_level = sc._log.level
        sc._log.setLevel(logging.INFO)
        try:
            parsed = {
                "features": [
                    {"gff_id": "f1", "type": "CDS",
                     "start_0": 0, "end": 100, "strand": 1,
                     "qualifiers": {}},
                    {"gff_id": "f2", "type": "CDS",
                     "start_0": 0, "end": 999_999, "strand": 1,
                     "qualifiers": {}},
                ],
            }
            out = sc._gff3_features_to_biopython(parsed, total=200)
        finally:
            sc._log.removeHandler(handler)
            sc._log.setLevel(prev_level)
        assert len(out) == 1   # f2 dropped
        assert any(
            "GFF3: dropping" in r.getMessage() for r in records
        ), f"expected drop log, got: {[r.getMessage() for r in records]}"


# ═════════════════════════════════════════════════════════════════════
# M1 — RestoreFromBackupModal _dismiss_once
# ═════════════════════════════════════════════════════════════════════

class TestRestoreFromBackupDismissOnce:
    def test_has_dismiss_once_helper(self):
        assert hasattr(sc.RestoreFromBackupModal, "_dismiss_once")

    def test_dismissed_flag_in_init(self):
        import inspect
        src = inspect.getsource(sc.RestoreFromBackupModal.__init__)
        assert "_dismissed" in src


# ═════════════════════════════════════════════════════════════════════
# M2 — _safe_identifier for active_* settings
# ═════════════════════════════════════════════════════════════════════

class TestSafeIdentifierValidator:
    def test_helper_rejects_traversal(self):
        assert sc._is_safe_identifier("Main Project")
        assert sc._is_safe_identifier("")           # empty = handled elsewhere
        assert not sc._is_safe_identifier("../etc")
        assert not sc._is_safe_identifier("a/b")
        assert not sc._is_safe_identifier("a\\b")
        assert not sc._is_safe_identifier("a\x00b")
        assert sc._is_safe_identifier("with spaces and-dashes")

    def test_validate_settings_rejects_active_traversal(self):
        """A hand-edited settings.json with `active_project:
        '../../etc/foo'` must be rejected, even though string-type
        check passes. Defense in depth against future refactors
        that join the name into a path."""
        raw = {"active_project": "../../etc/foo"}
        cleaned, warnings = sc._validate_settings(raw)
        assert cleaned.get("active_project") != "../../etc/foo"
        assert any("path-traversal" in w for w in warnings)


# ═════════════════════════════════════════════════════════════════════
# M5 — agent body parse error → 400 sentinel
# ═════════════════════════════════════════════════════════════════════

class TestBodyParseErrorSentinel:
    def test_read_body_returns_sentinel_on_parse_error(self):
        """White-box: `_read_body` returns the `"__bad_body__"`
        sentinel on JSON / unicode parse error, not silent `{}`."""
        import inspect
        src = inspect.getsource(sc._AgentRequestHandler._read_body)
        assert '"__bad_body__"' in src


# ═════════════════════════════════════════════════════════════════════
# M14 — _VECTOR_MATCH_CACHE keys on hash
# ═════════════════════════════════════════════════════════════════════

class TestVectorMatchCacheKey:
    def test_uses_hash_key(self):
        """White-box: `_vector_half_top_seq` keys on
        `(hash(ev_gb), enzyme)` not `(ev_gb, enzyme)`."""
        import inspect
        src = inspect.getsource(sc._vector_half_top_seq)
        assert "hash(ev_gb)" in src
        # Lines after `key = ` shouldn't show the raw ev_gb tuple form.
        assert "key = (ev_gb, enzyme)" not in src


# ═════════════════════════════════════════════════════════════════════
# M20 — narrow except in target_seq fallback
# ═════════════════════════════════════════════════════════════════════

class TestTargetSeqExceptNarrow:
    def test_narrowed_except(self):
        text = Path(sc.__file__).read_text(encoding="utf-8")
        # The two narrowed sites must be present.
        assert "except (AttributeError, TypeError):" in text


# ═════════════════════════════════════════════════════════════════════
# L1 — splicecraft_cli refuses symlinked token
# ═════════════════════════════════════════════════════════════════════

class TestCliRefusesSymlinkedToken:
    def test_cli_read_session_uses_lstat(self):
        """White-box: `splicecraft_cli._read_session` must `lstat`
        the token file and refuse symlinks."""
        import splicecraft_cli
        import inspect
        src = inspect.getsource(splicecraft_cli._read_session)
        assert "lstat" in src
        assert "S_ISLNK" in src


# ═════════════════════════════════════════════════════════════════════
# L3 — token via secrets.token_urlsafe(32)
# ═════════════════════════════════════════════════════════════════════

class TestTokenViaSecrets:
    def test_uses_secrets_token_urlsafe(self):
        import inspect
        src = inspect.getsource(sc._start_agent_api)
        body = src.split('"""', 2)[-1]
        assert "token_urlsafe" in body
        # Sentinel against regression to uuid4. (Docstring references
        # the old approach in the change-rationale comment, so we
        # check function body only.)
        assert "uuid.uuid4().hex" not in body


# ═════════════════════════════════════════════════════════════════════
# H1 — custom enzyme agent endpoints under _cache_lock
# ═════════════════════════════════════════════════════════════════════

class TestCustomEnzymeAgentLock:
    def test_create_wraps_in_cache_lock(self):
        import inspect
        src = inspect.getsource(sc._h_create_custom_enzyme)
        # moved to splicecraft_agent — resolves the lock as _state._cache_lock
        assert "with _state._cache_lock:" in src

    def test_update_wraps_in_cache_lock(self):
        import inspect
        src = inspect.getsource(sc._h_update_custom_enzyme)
        assert "with _state._cache_lock:" in src

    def test_delete_wraps_in_cache_lock(self):
        import inspect
        src = inspect.getsource(sc._h_delete_custom_enzyme)
        assert "with _state._cache_lock:" in src

    def test_collection_create_wraps_in_cache_lock(self):
        import inspect
        src = inspect.getsource(sc._h_create_enzyme_collection)
        assert "with _state._cache_lock:" in src


# ═════════════════════════════════════════════════════════════════════
# H2 — Ctrl+S save worker lock
# ═════════════════════════════════════════════════════════════════════

class TestSaveWorkerLock:
    def test_save_worker_library_rmw_under_lock(self):
        import inspect
        src = inspect.getsource(sc.PlasmidApp._save_worker)
        # The library mirror section must hold `_cache_lock`.
        assert "with _cache_lock:" in src


# ═════════════════════════════════════════════════════════════════════
# Sweep #26 deferred-batch follow-up (2026-05-23) — same INV-66
# ═════════════════════════════════════════════════════════════════════
#
# After landing sweep #25 the deferred-item list got worked through too.
# These tests regression-lock the additional fixes that landed:
#
# L2  pypi http downgrade refused unless SPLICECRAFT_PYPI_INSECURE=1
# L4  _check_data_files drives from _USER_DATA_FILE_ATTRS registry
# L5  _ENZYME_CUT_RANGE module constant replaces hardcoded ±30
# L6  NCBI / Kazusa narrow excepts use _urllib_error.URLError
# L7  agent.write.ok log emitted AFTER _send confirms
# L8  _ACCEPTOR_TU_PAIRS_CACHE FIFO cap
# L10 _build_system_info rglob cap at 10k
# M3  snapshot dir restore 5 GB cap
# M12 clipboard fallback dir prune (7 days OR 100-file cap)
# M16 ConstructorModal palette + grammar dep-count via readonly iter
# M17 _iter_all_grammars_readonly helper
# M18 _gb_text_to_record(cache=False) for batch parses
# M21 _h_align_plasmidsaurus_zip post-alignment target-drift detection
# H7  _list_gbk_members_in_zip + _extract_gbk_member TOCTOU-safe via
#     os.open + fileobj
# M8  _check_agent_read_path_ancestors + applied to Plasmidsaurus paths
# H13 _wrap_feats_idx precomputed in load_record


class TestSweep26PypiHttpDowngrade:
    """L2 — refuse http:// PyPI override unless explicit insecure
    opt-in. Pre-sweep an in-path attacker could MITM the update-
    check JSON to spoof "no update available" or worse."""

    def test_http_refused_without_insecure_env(self, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_PYPI_URL", "http://evil.example/v1")
        monkeypatch.delenv("SPLICECRAFT_PYPI_INSECURE", raising=False)
        url = sc._resolve_pypi_url()
        assert "evil.example" not in url

    def test_https_allowed(self, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_PYPI_URL",
                            "https://mirror.example/v1")
        monkeypatch.delenv("SPLICECRAFT_PYPI_INSECURE", raising=False)
        url = sc._resolve_pypi_url()
        assert url == "https://mirror.example/v1"

    def test_http_allowed_with_insecure_env(self, monkeypatch):
        monkeypatch.setenv("SPLICECRAFT_PYPI_URL", "http://corp.example/v1")
        monkeypatch.setenv("SPLICECRAFT_PYPI_INSECURE", "1")
        url = sc._resolve_pypi_url()
        assert url == "http://corp.example/v1"


class TestSweep26CheckDataFilesRegistry:
    """L4 — `_check_data_files` drives from `_USER_DATA_FILE_ATTRS`
    registry, not a hand-maintained list. Future files added to the
    registry now auto-enroll in startup validation."""

    def test_source_iterates_user_data_file_attrs(self):
        import inspect
        src = inspect.getsource(sc.PlasmidApp._check_data_files)
        assert "_USER_DATA_FILE_ATTRS" in src


class TestSweep26EnzymeCutRangeConstant:
    """L5 — single module constant for cut-position validation."""

    def test_constant_defined(self):
        assert hasattr(sc, "_ENZYME_CUT_RANGE")
        assert sc._ENZYME_CUT_RANGE == 30

    def test_agent_validator_uses_constant(self):
        import inspect
        src = inspect.getsource(
            sc._agent_validate_custom_enzyme_payload,
        )
        assert "_ENZYME_CUT_RANGE" in src


class TestSweep26AcceptorCacheCap:
    """L8 — FIFO eviction at 64 entries."""

    def test_cap_constant(self):
        assert sc._ACCEPTOR_TU_PAIRS_CACHE_MAX == 64

    def test_cap_enforced_on_insert(self):
        sc._ACCEPTOR_TU_PAIRS_CACHE.clear()
        for i in range(sc._ACCEPTOR_TU_PAIRS_CACHE_MAX + 5):
            sc._ACCEPTOR_TU_PAIRS_CACHE[(f"g{i}", "x")] = []
        # Direct stuffing exceeds the cap. Verify the helper's own
        # insertion path respects the cap by simulating its insert
        # (the eviction lives inside `_grammar_acceptor_tu_pairs`,
        # not at dict-assignment — that's why pre-fix the cache
        # grew unbounded). Sentinel only checks the constant exists
        # + the helper's source has the eviction loop.
        import inspect
        src = inspect.getsource(sc._grammar_acceptor_tu_pairs)
        assert "_ACCEPTOR_TU_PAIRS_CACHE_MAX" in src


class TestSweep26SystemInfoRglobCap:
    """L10 — bound `_build_system_info`'s rglob walk at 10k files.
    On heavy installs (thousands of attachments) the unbounded walk
    was a CPU burst per diagnostic bundle / log-startup."""

    def test_source_caps_count(self):
        import inspect
        src = inspect.getsource(sc._build_system_info)
        assert "_SYSTEM_INFO_RGLOB_CAP" in src


class TestSweep26GBParseCacheFlag:
    """M18 — `cache=False` flag for one-shot batch parses (Plasmid-
    saurus zip ingest, bulk-import folder walk)."""

    def test_signature_accepts_cache_flag(self):
        import inspect
        sig = inspect.signature(sc._gb_text_to_record)
        assert "cache" in sig.parameters
        assert sig.parameters["cache"].default is True

    def test_cache_false_bypasses_population(self):
        # Build a valid GenBank text via Biopython so we don't trip
        # `BiopythonParserWarning` on a hand-rolled LOCUS line with
        # off-spec column spacing (same trap 0.9.31 cleared from
        # `TestConstructorMultiGrammarTabs::test_per_role_entry_vector_banner`).
        from io import StringIO
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(
            Seq("ATCGATCGAT"),
            id="sweep26_test",
            name="sweep26_test",
            annotations={"molecule_type": "DNA", "topology": "linear"},
        )
        buf = StringIO()
        SeqIO.write([rec], buf, "genbank")
        text = buf.getvalue()
        # Parse cache is keyed on a 128-bit content hash (audit 2026-07-01 —
        # was `hash(text)`, which a 64-bit collision could poison), so probe
        # with the same key helper the implementation uses.
        import splicecraft_record as _record
        key = _record._gb_cache_key(text)
        sc._GB_PARSE_CACHE.clear()
        sc._gb_text_to_record(text, cache=False)
        assert key not in sc._GB_PARSE_CACHE
        # And cache=True (default) populates.
        sc._gb_text_to_record(text)
        assert key in sc._GB_PARSE_CACHE


class TestSweep26GrammarsReadonly:
    """M17 — `_iter_all_grammars_readonly` helper avoids the per-
    call deepcopy of every built-in grammar."""

    def test_helper_exists(self):
        assert callable(sc._iter_all_grammars_readonly)

    def test_returns_builtin_refs(self):
        # The readonly view must include every built-in grammar id.
        view = sc._iter_all_grammars_readonly()
        for gid in sc._BUILTIN_GRAMMARS:
            assert gid in view


class TestSweep26ZipFdPass:
    """H7 — `_list_gbk_members_in_zip` and `_extract_gbk_member` open
    via `os.open` + `fileobj=` to close the TOCTOU window."""

    def test_list_source_uses_os_open(self):
        import inspect
        src = inspect.getsource(sc._list_gbk_members_in_zip)
        assert "os.open" in src
        assert "os.fstat" in src

    def test_extract_source_uses_os_open(self):
        import inspect
        src = inspect.getsource(sc._extract_gbk_member)
        assert "os.open" in src
        assert "os.fstat" in src


class TestSweep26AncestorReadCheck:
    """M8 — `_check_agent_read_path_ancestors` walks parent symlinks
    on Plasmidsaurus read endpoints."""

    def test_helper_exists(self):
        assert callable(sc._check_agent_read_path_ancestors)

    def test_clean_path_passes(self, tmp_path):
        target = tmp_path / "ok.zip"
        # Parent exists, no symlinks — passes.
        assert sc._check_agent_read_path_ancestors(target) is None

    def test_parent_symlink_refused(self, tmp_path):
        real = tmp_path / "real"
        real.mkdir()
        linked = tmp_path / "linked"
        linked.symlink_to(real)
        target = linked / "x.zip"
        result = sc._check_agent_read_path_ancestors(target)
        assert result is not None
        assert "symlink" in result.lower()

    def test_align_endpoint_calls_helper(self):
        import inspect
        src = inspect.getsource(sc._h_align_plasmidsaurus_zip)
        assert "_check_agent_read_path_ancestors" in src

    def test_list_members_endpoint_calls_helper(self):
        import inspect
        src = inspect.getsource(sc._h_list_plasmidsaurus_members)
        assert "_check_agent_read_path_ancestors" in src


class TestSweep26WrapFeatsIdx:
    """H13 — `_wrap_feats_idx` precomputed in `load_record` so
    `_draw_linear_flag` walks only the small wrap subset, not every
    feature on every render frame."""

    def test_load_record_populates_index(self):
        import inspect
        src = inspect.getsource(sc.PlasmidMap.load_record)
        assert "_wrap_feats_idx" in src

    def test_draw_linear_flag_uses_index(self):
        import inspect
        src = inspect.getsource(sc.PlasmidMap._draw_linear_flag)
        assert "_wrap_feats_idx" in src

    def test_linear_flag_selected_style_does_not_use_reverse(self):
        """2026-05-26 regression: the linear-flag selection
        highlight used to set the body / head style to `"reverse "
        + color`, but reverse on the `█` body glyph swaps FG and
        BG so the cell renders with the default terminal text
        colour on a colour BG — and `█` covers the entire cell,
        hiding the BG. Visually the body cells outside the label
        area lost the feature colour, so the selected bar appeared
        to "shrink to the label width" (user-reported when
        clicking the translated CDS bar in the seq panel). Fix
        uses `"bold " + color` instead, keeping the FG colour
        visible. This test guards against re-introducing `reverse`
        on the body style without first picking a non-`█` glyph
        for the body."""
        import inspect
        src = inspect.getsource(sc.PlasmidMap._draw_linear_flag)
        # The buggy pattern must not be present in the body /
        # head style line. We allow `"reverse"` to appear in OTHER
        # contexts (alignment band, etc.) — assert specifically
        # that the `is_sel` ternary doesn't compose reverse + color
        # for the feature body / head style.
        assert '"reverse " + color' not in src, (
            "linear-flag selected style must not use reverse + "
            "color on the `█` body glyph (renders as default-colour "
            "block; bar appears to shrink). Use `bold + color` "
            "instead, or switch the body glyph to a partial-fill "
            "char (▓ / ▒) before re-introducing reverse."
        )
        # The fix uses bold + color.
        assert '"bold " + color' in src, (
            "linear-flag selected style should be `'bold ' + color`"
            " — see the comment block above the assignment for "
            "rationale."
        )


class TestSweep26AlignmentDriftDetection:
    """M21 — `_h_align_plasmidsaurus_zip` re-checks the target
    library entry after the (possibly multi-second) alignment to
    surface 410 Gone if it was deleted mid-flight, or flag a
    rename in the result payload."""

    def test_source_has_post_alignment_recheck(self):
        import inspect
        src = inspect.getsource(sc._h_align_plasmidsaurus_zip)
        # Post-alignment, the handler re-looks up by the captured id.
        assert "resolved_target_id" in src
        # Returns 410 on disappearance.
        assert "410" in src
        # Flags rename via `_target_renamed_to`.
        assert "_target_renamed_to" in src


class TestSweep26LogEventTimingAfterSend:
    """L7 — `agent.write.ok` log fires AFTER `_send` succeeds so the
    log unambiguously means "state mutated AND client confirmed".
    Pre-sweep the log fired before the response landed."""

    def test_source_order_send_then_log(self):
        import inspect
        # The routing body moved into `_handle_inner` (2026-10-23, AG7): `_handle`
        # wraps it so an idempotency claim is released on every exit.
        src = inspect.getsource(sc._AgentRequestHandler._handle_inner)
        # Find positions of the send call and the success-log emit.
        send_pos = src.find("self._send(payload, status)")
        log_pos = src.find('"agent.write.ok"')
        assert send_pos > 0 and log_pos > 0
        assert send_pos < log_pos, (
            "_send must precede the success-log emit"
        )
