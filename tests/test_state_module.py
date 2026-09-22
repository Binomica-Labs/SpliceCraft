"""Guard the shared-state pattern (the Option-3 state decoupling).

Mutable process state that the monolith kept as bare module globals and mutated
at runtime (render-tier flags now; data paths / caches later) is migrating into
``splicecraft_state``, accessed by attribute (``state.X``) so the hub, every
sibling, and the tests all read and write ONE copy.

The failure this guards against is the *stale-shadow* trap: a migrated name also
exists as a hub-level global (a leftover ``_X = ...`` in splicecraft.py, or a
``from splicecraft_state import _X``), so a monkeypatch or a runtime write to one
copy silently misses readers using the other. That trap is exactly what makes
naive module extraction unsafe in this codebase, so it gets a dedicated test
that every future migration extends via ``_MIGRATED``.
"""
from __future__ import annotations

import splicecraft as sc
import splicecraft_state

# Names migrated out of the hub into splicecraft_state. APPEND to this as the
# decoupling proceeds; each entry is then held to the single-source-of-truth
# invariants below.
_MIGRATED = [
    "_ASCII_MODE",
    "_ASCII_FORCED",
    "_WIN_UTF8_CONSOLE",
    "_ESCAPE_ASPECT",
    # Phase A1: caches, generation counters, background-write coordination
    "_BLAST_CACHE_GENERATION",
    "_DANGLING_ACTIVE_COLLECTION_NAME",
    "_SPELLCHECK_ENGINE",
    "_WHATS_NEW_CACHE",
    "_collection_sync_pending",
    "_collection_sync_thread",
    "_feature_library_index_cache",
    "_features_generation",
    "_primer_usage_cache",
    "_primer_usage_cache_gen",
    "_settings_flush_pending",
    "_settings_flush_running",
    # Phase A2a: persisted-data caches
    "_codon_tables_cache", "_collections_cache", "_custom_enzymes_cache",
    "_entry_vectors_cache", "_enzyme_collections_cache", "_experiment_projects_cache",
    "_experiments_cache", "_feature_colors_cache", "_features_cache", "_gels_cache",
    "_grammars_cache", "_hmm_db_catalog_cache", "_library_cache",
    "_model_collections_cache", "_parts_bin_cache",
    "_parts_bin_collections_cache", "_primer_collections_cache", "_primers_cache",
    "_protein_collections_cache", "_protein_motifs_cache", "_settings_cache",
    # Phase A2b: one-shot backfill / name-trim flags
    "_collections_backfill_done", "_collections_origin_history_backfill_done",
    "_entry_vectors_name_trim_done", "_id_name_backfill_done",
    "_origin_history_backfill_done", "_parts_bin_collections_backfill_done",
    "_parts_bin_sequence_backfill_done", "_primer_collections_backfill_done",
    "_primers_name_trim_done",
    # Phase A2c: the save-authorization chokepoint flag
    "_SAVES_AUTHORIZED", "_SAVES_AUTHORIZED_REASON",
    # Phase B-prep: the data directory + the persisted-file path constants
    "_DATA_DIR", "_DATA_DIR_NOTICE",
    "_AGENT_TOKEN_FILE", "_CODON_TABLES_FILE", "_COLLECTIONS_FILE",
    "_CUSTOM_ENZYMES_FILE", "_DATA_VERSION_FILE", "_ENTRY_VECTORS_FILE",
    "_ENZYME_COLLECTIONS_FILE", "_EXPERIMENTS_FILE", "_EXPERIMENT_PROJECTS_FILE",
    "_FEATURES_FILE", "_FEATURE_COLORS_FILE", "_GELS_FILE", "_GRAMMARS_FILE",
    "_HMM_DB_CATALOG_FILE", "_LIBRARY_FILE", "_MODEL_COLLECTIONS_FILE",
    "_PARTS_BIN_COLLECTIONS_FILE",
    "_PARTS_BIN_FILE", "_PRIMERS_FILE", "_PRIMER_COLLECTIONS_FILE",
    "_PROTEIN_COLLECTIONS_FILE", "_PROTEIN_MOTIFS_FILE", "_SETTINGS_FILE",
    # Phase B-main: persistence-engine config tunables (monkeypatched by the
    # data-safety suite, so a single shared copy in _state is mandatory)
    "_CURRENT_SCHEMA_VERSION", "_OBSERVED_SCHEMA_VERSIONS",
    "_BACKUP_RETENTION_COUNT", "_BACKUP_TOTAL_SIZE_CAP_BYTES", "_BACKUP_MIN_KEEP",
    "_LOST_ENTRIES_DIR_NAME", "_LOST_ENTRIES_RETENTION_COUNT",
    "_LOST_ENTRIES_TOTAL_SIZE_CAP_BYTES", "_SAFE_LOAD_JSON_MAX_BYTES",
    "_SAVE_READBACK_FULL_PARSE_MAX_BYTES", "_SNAPSHOT_DIR_NAME",
]


def test_hub_exposes_state_module():
    """`sc._state` is the canonical handle the hub + tests patch through."""
    assert sc._state is splicecraft_state


def test_migrated_names_live_in_state():
    missing = [n for n in _MIGRATED if not hasattr(splicecraft_state, n)]
    assert not missing, f"migrated name(s) absent from splicecraft_state: {missing}"


def test_no_stale_hub_shadow():
    """A migrated name must NOT also exist in the hub's own namespace, where it
    would desync from `state.X` writes / monkeypatches."""
    shadows = [n for n in _MIGRATED if n in vars(sc)]
    assert not shadows, (
        f"`splicecraft` keeps stale shadow(s) of migrated state: {shadows}. "
        "Remove the hub-level binding and read via `_state.<name>` -- a "
        "leftover global or a by-value `from splicecraft_state import` desyncs "
        "from the live value (the trap this whole refactor exists to avoid)."
    )


def test_runtime_writer_mutates_the_shared_copy(monkeypatch):
    """`_set_ascii_mode` (a runtime writer) must flip the SAME storage readers
    see -- i.e. the state module's attribute, not a hub-local copy."""
    monkeypatch.setattr(splicecraft_state, "_ASCII_MODE", False)
    sc._set_ascii_mode(True)
    assert splicecraft_state._ASCII_MODE is True, (
        "_set_ascii_mode did not write the shared splicecraft_state copy"
    )
    sc._set_ascii_mode(False)
    assert splicecraft_state._ASCII_MODE is False


def test_all_state_hooks_registered():
    """Every `_state._*_hook` must be REGISTERED once the hub has imported -- it
    points at a real hub (or sibling) function, not its `splicecraft_state`
    fail-loud / None placeholder. The hub+siblings architecture leans on these:
    a sibling reaches a hub-pinned dependency through `_state._X_hook`, never by
    importing the hub (a cycle). If a future refactor adds a hook slot but forgets
    the `_state._X_hook = ...` registration (or drops an existing one), the
    dependent endpoints / accessors break at runtime; this catches it at import."""
    hooks = [n for n in dir(splicecraft_state) if n.endswith("_hook")]
    assert hooks, "expected _*_hook attributes on the state module"
    unregistered = [
        n for n in hooks
        if (fn := getattr(splicecraft_state, n)) is None
        or getattr(fn, "__module__", "") == "splicecraft_state"
    ]
    assert not unregistered, (
        f"_state hook(s) left unregistered after `import splicecraft`: "
        f"{unregistered}. The hub must register each at import "
        "(`_state._X_hook = _X`); an unregistered hook means the sibling "
        "endpoints / accessors that call it fail (fail-loud default)."
    )
