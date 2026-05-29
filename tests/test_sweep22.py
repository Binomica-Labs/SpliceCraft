"""
test_sweep22 — load-bearing function audit + file migration framework
(2026-05-21).

User asked for a deep audit of the most load-bearing functions
followed by a futureproofing pass for users bringing files from
previous SpliceCraft versions. This file regression-locks:

  * `_safe_save_json` ancestor-symlink walk — pre-sweep, only the
    target path was checked, NOT its ancestor chain. Sweep #10 added
    the walk to `_check_agent_write_path` but didn't propagate to
    `_safe_save_json`. A symlink at any deeper ancestor could
    redirect every save. Now mirrors the agent-side defense.

  * `_iupac_pattern` case-fold cache key — pre-sweep, `"gaattc"` and
    `"GAATTC"` occupied two cache slots producing identical regex
    objects. Wastes 1 slot per mixed-case variant. Cache key now
    uppercased on lookup AND store.

  * File migration framework (`_extract_entries` + `_migrate_entries`
    + `_CURRENT_SCHEMA_VERSION`) — infrastructure existed since at
    least 0.5.x (invariant #36) but has NEVER been exercised in
    production because no migrations are registered. When the schema
    eventually bumps, the framework needs to work correctly. These
    tests exercise every path so a future bump can ship with
    confidence.
"""
from __future__ import annotations

import json

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# _safe_save_json — ancestor-symlink walk
# ═══════════════════════════════════════════════════════════════════════════════

class TestSafeSaveJsonAncestorSymlink:
    """Pre-sweep regression: a symlink at a parent / grandparent
    directory would redirect writes. Now refused outright with a
    clear OSError that mentions the symlinked ancestor.
    """

    def test_direct_path_symlink_refused(self, tmp_path):
        # Direct symlink at the target is rejected (pre-existing
        # behaviour — invariant from 0.8.9 fix).
        target_dir = tmp_path / "data"
        target_dir.mkdir()
        real_file = tmp_path / "real.json"
        real_file.write_text("{}", encoding="utf-8")
        symlinked = target_dir / "library.json"
        symlinked.symlink_to(real_file)
        with pytest.raises(OSError, match=r"symlink"):
            sc._safe_save_json(symlinked, [{"id": "x"}], "Test")

    def test_parent_symlink_refused(self, tmp_path):
        # Sweep #22: symlink at the IMMEDIATE parent should refuse.
        real_dir = tmp_path / "real_data"
        real_dir.mkdir()
        linked_dir = tmp_path / "linked_data"
        linked_dir.symlink_to(real_dir)
        target = linked_dir / "library.json"
        with pytest.raises(OSError, match=r"symlink"):
            sc._safe_save_json(target, [{"id": "x"}], "Test")

    def test_grandparent_symlink_refused(self, tmp_path):
        # Sweep #22: symlink at any DEEPER ancestor must also refuse.
        real_grandparent = tmp_path / "real_gp"
        real_grandparent.mkdir()
        linked_gp = tmp_path / "linked_gp"
        linked_gp.symlink_to(real_grandparent)
        # Create the intermediate dir under the symlink so the path
        # is fully formed.
        (linked_gp / "subdir").mkdir()
        target = linked_gp / "subdir" / "library.json"
        with pytest.raises(OSError, match=r"symlink"):
            sc._safe_save_json(target, [{"id": "x"}], "Test")

    def test_no_symlink_succeeds(self, tmp_path):
        # Baseline: a path with NO symlinks anywhere in its ancestry
        # should save normally. Catches a too-aggressive narrowing.
        (tmp_path / "subdir").mkdir()
        target = tmp_path / "subdir" / "library.json"
        sc._safe_save_json(target, [{"id": "x"}], "Test")
        assert target.exists()


# ═══════════════════════════════════════════════════════════════════════════════
# _iupac_pattern — case-fold cache key
# ═══════════════════════════════════════════════════════════════════════════════

class TestIupacPatternCaseFold:
    """Pre-sweep regression: mixed-case calls to `_iupac_pattern`
    occupied separate cache slots despite producing identical regex.
    """

    def setup_method(self):
        # Test in isolation — clear the cache so prior calls from
        # other tests don't pollute the count.
        from splicecraft_biology import _PATTERN_CACHE
        _PATTERN_CACHE.clear()

    def test_lowercase_and_uppercase_hit_same_slot(self):
        from splicecraft_biology import _iupac_pattern, _PATTERN_CACHE
        p1 = _iupac_pattern("GAATTC")
        p2 = _iupac_pattern("gaattc")
        # Same compiled regex returned for either case.
        assert p1 is p2
        # And only ONE cache entry exists.
        assert len(_PATTERN_CACHE) == 1
        # Stored under the uppercased key.
        assert "GAATTC" in _PATTERN_CACHE

    def test_mixed_case_normalized(self):
        from splicecraft_biology import _iupac_pattern, _PATTERN_CACHE
        p1 = _iupac_pattern("GaAtTc")
        p2 = _iupac_pattern("gAaTtC")
        assert p1 is p2
        assert len(_PATTERN_CACHE) == 1

    def test_pattern_match_correct_after_normalization(self):
        from splicecraft_biology import _iupac_pattern
        pat = _iupac_pattern("gaattc")
        # Pattern matches the uppercase form (sites are uppercased
        # in the underlying regex regardless of cache-key case).
        assert pat.search("AAGAATTCAA") is not None


# ═══════════════════════════════════════════════════════════════════════════════
# File-migration framework — every load path runs migrations
# ═══════════════════════════════════════════════════════════════════════════════

class TestMigrationFramework:
    """The framework lives in `_extract_entries` + `_migrate_entries`.
    It's been in place since at least 0.5.x but has never been
    exercised in production because no migrations are registered.
    These tests verify each load path through the framework so a
    future schema bump can ship with confidence.
    """

    def test_bare_list_legacy_accepted(self):
        # Pre-0.3.1 format: top-level JSON was a bare list, no
        # envelope. `_extract_entries` recognizes this as v0 and
        # runs the v0→current migrator (which is currently a no-op).
        raw = [{"id": "a"}, {"id": "b"}]
        entries, warning = sc._extract_entries(raw, "Test")
        assert entries == raw
        assert warning is None

    def test_envelope_v1_loaded(self):
        raw = {"_schema_version": 1, "entries": [{"id": "a"}]}
        entries, warning = sc._extract_entries(raw, "Test")
        assert entries == [{"id": "a"}]
        assert warning is None

    def test_envelope_future_version_warns_but_loads(self):
        raw = {"_schema_version": 99, "entries": [{"id": "a"}]}
        entries, warning = sc._extract_entries(raw, "Test")
        assert entries == [{"id": "a"}]
        assert warning is not None
        assert "newer SpliceCraft" in warning

    def test_envelope_unknown_shape_returns_none(self):
        # Garbage shape — `_safe_load_json` falls through to .bak.
        raw = {"not_entries": []}
        entries, warning = sc._extract_entries(raw, "Test")
        assert entries is None
        assert warning is not None

    def test_envelope_with_extra_fields_preserved(self):
        # A future SpliceCraft might add `_origin_app` or similar
        # envelope fields. Loading them on older binaries should NOT
        # strip them — entries pass through deepcopy-aware.
        raw = {
            "_schema_version": 1,
            "_origin_app": "future-build",
            "entries": [{"id": "a"}],
        }
        entries, _ = sc._extract_entries(raw, "Test")
        # Entries themselves unchanged.
        assert entries == [{"id": "a"}]

    def test_migrate_entries_chains_versions(self, monkeypatch):
        # Inject hypothetical migrators 0→1 and 1→2 to verify the
        # walk applies them in order. Tests the chain mechanic
        # without committing to a real schema bump.
        v0_to_v1 = lambda e: {**e, "v1_field": True}  # noqa: E731
        v1_to_v2 = lambda e: {**e, "v2_field": True}  # noqa: E731
        monkeypatch.setitem(
            sc._ENTRY_MIGRATIONS, "Test",
            {(0, 1): v0_to_v1, (1, 2): v1_to_v2},
        )
        out, warnings = sc._migrate_entries(
            [{"id": "a"}], 0, 2, "Test",
        )
        assert out == [{"id": "a", "v1_field": True, "v2_field": True}]
        assert warnings == []

    def test_migrate_entries_missing_step_is_noop(self, monkeypatch):
        # Schema bumps that are purely additive (new optional field)
        # don't need an explicit migrator — entries pass through.
        v0_to_v1 = lambda e: {**e, "v1_field": True}  # noqa: E731
        # NO 1→2 registered.
        monkeypatch.setitem(
            sc._ENTRY_MIGRATIONS, "Test",
            {(0, 1): v0_to_v1},
        )
        out, _ = sc._migrate_entries([{"id": "a"}], 0, 2, "Test")
        # 0→1 applied; 1→2 no-op.
        assert out == [{"id": "a", "v1_field": True}]

    def test_migrate_entries_failed_migrator_keeps_entry(self, monkeypatch):
        # A migrator that raises must NOT drop the entry — keep it
        # in its pre-migration shape and surface a warning. Better
        # to load a stale-shaped entry than lose the user's data.
        def bad_migrator(entry):
            raise ValueError("intentional test failure")

        monkeypatch.setitem(
            sc._ENTRY_MIGRATIONS, "Test",
            {(0, 1): bad_migrator},
        )
        out, warnings = sc._migrate_entries(
            [{"id": "a"}], 0, 1, "Test",
        )
        # Entry preserved verbatim.
        assert out == [{"id": "a"}]
        # Warning surfaced.
        assert len(warnings) == 1
        assert "migration failed" in warnings[0].lower()

    def test_migrate_entries_non_dict_skipped(self, monkeypatch):
        # A bare-list legacy file might contain non-dict garbage
        # (string entries, e.g.). Skip them with a warning rather
        # than crashing the loader.
        v0_to_v1 = lambda e: {**e, "v1_field": True}  # noqa: E731
        monkeypatch.setitem(
            sc._ENTRY_MIGRATIONS, "Test",
            {(0, 1): v0_to_v1},
        )
        out, warnings = sc._migrate_entries(
            [{"id": "a"}, "not a dict", {"id": "b"}],
            0, 1, "Test",
        )
        # Bad entry skipped; good entries migrated.
        assert {"id": "a", "v1_field": True} in out
        assert {"id": "b", "v1_field": True} in out
        assert len(warnings) >= 1

    def test_migrate_entries_no_op_when_at_current(self):
        # `from_version >= to_version` — nothing to do, return a
        # fresh list (not the input itself, for caller independence).
        original = [{"id": "a"}]
        out, warnings = sc._migrate_entries(original, 1, 1, "Test")
        assert out == original
        assert out is not original   # Fresh list, no alias.

    def test_migrate_entries_non_list_input(self):
        # Defensive: non-list payload returns empty list + warning.
        out, warnings = sc._migrate_entries(
            {"not": "a list"}, 0, 1, "Test",  # type: ignore[arg-type]
        )
        assert out == []
        assert len(warnings) == 1


# ═══════════════════════════════════════════════════════════════════════════════
# Round-trip: legacy bare-list file → save → loaded as envelope
# ═══════════════════════════════════════════════════════════════════════════════

class TestLegacyFileRoundTrip:
    """End-to-end test for users importing pre-envelope SpliceCraft
    JSON files: the loader should normalize them through migration
    AND the next save should re-shape them into the envelope format.
    Closes the user's "files from previous versions" futureproofing
    concern.
    """

    def test_bare_list_load_then_save_emits_envelope(self, tmp_path):
        # Write a bare-list file (pre-0.3.1 SpliceCraft format).
        bare_path = tmp_path / "library.json"
        bare_path.write_text(
            json.dumps([{"id": "legacy_a"}, {"id": "legacy_b"}]),
            encoding="utf-8",
        )
        # Load via `_safe_load_json` — should accept the bare list.
        entries, _ = sc._safe_load_json(bare_path, "Test")
        assert len(entries) == 2
        # Now save back through `_safe_save_json` (which always
        # emits the envelope format).
        sc._safe_save_json(bare_path, entries, "Test")
        # Verify the on-disk format is now envelope.
        reloaded = json.loads(bare_path.read_text(encoding="utf-8"))
        assert isinstance(reloaded, dict)
        assert "_schema_version" in reloaded
        assert "entries" in reloaded
        assert len(reloaded["entries"]) == 2


# ═══════════════════════════════════════════════════════════════════════════════
# Coverage check: every persisted file routes through _safe_load_json
# ═══════════════════════════════════════════════════════════════════════════════

class TestMigrationCoverageEveryLoadPath:
    """White-box check: every `_load_*` helper should route through
    `_safe_load_json` so the migration framework applies uniformly.
    A future load helper that bypasses this would silently skip
    migrations when a schema bump lands.
    """

    def test_every_load_helper_uses_safe_load_json(self):
        import inspect
        # Map of (loader_name, expected source-contains substring).
        loaders = [
            "_load_library", "_load_collections", "_load_parts_bin",
            "_load_parts_bin_collections", "_load_primers",
            "_load_features", "_load_protein_motifs",
            "_load_feature_colors", "_load_custom_grammars",
            "_load_entry_vectors", "_load_settings",
            "_load_experiments", "_load_experiment_projects",
            "_load_gels",
            # Added 2026-05-29: persisted types that post-date the
            # original list — custom enzymes + enzyme collections
            # (INV-85) and codon tables — must route through the
            # framework too, else a future schema bump silently skips
            # them. NB `_codon_tables_load` does NOT follow the
            # `_load_*` prefix; loaders are resolved by name, not by
            # pattern, so the asymmetry is fine.
            "_load_custom_enzymes", "_load_enzyme_collections",
            "_codon_tables_load",
        ]
        # INV-73 (2026-05-25): walk the call graph one level to follow
        # cache-refactor indirections. `_load_library` no longer
        # inlines `_safe_load_json` — it routes through
        # `_ensure_library_cache_populated_and_migrated` (INV-72 cache
        # rewrite). The migration guarantee still holds; the test
        # needs to follow the indirection so it doesn't fail on a
        # legitimate refactor.
        for name in loaders:
            fn = getattr(sc, name, None)
            assert fn is not None, f"missing load helper: {name}"
            src = inspect.getsource(fn)
            if "_safe_load_json" in src:
                continue
            # Follow one level of indirection: find any helper called
            # in the body that itself routes through _safe_load_json.
            found = False
            for line in src.splitlines():
                # crude but effective: any `_some_helper_name(` call
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                for tok in stripped.replace("(", " ").split():
                    if not tok.startswith("_") or tok in {
                        name, "_typed_clone", "_log", "_cache_lock",
                    }:
                        continue
                    target = getattr(sc, tok, None)
                    if target is None or not callable(target):
                        continue
                    try:
                        sub_src = inspect.getsource(target)
                    except (TypeError, OSError):
                        continue
                    if "_safe_load_json" in sub_src:
                        found = True
                        break
                if found:
                    break
            assert found, (
                f"{name} doesn't route through _safe_load_json — "
                f"will silently skip migrations on schema bump. "
                f"Add the call or document the exception."
            )
