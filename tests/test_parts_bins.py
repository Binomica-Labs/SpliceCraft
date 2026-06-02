"""
test_parts_bins — multi-bin parts storage.

Covers the data layer that backs the parts-bin picker UI:

  * `_load_parts_bin_collections` / `_save_parts_bin_collections`
    round-trip schema preservation.
  * `_ensure_default_parts_bin` migration (existing `parts_bin.json`
    contents wrap into "Main Parts Bin" on first launch; idempotent
    on subsequent launches).
  * Active-bin pointer setter / getter / `_find_parts_bin` /
    `_parts_bin_name_taken`.
  * `_save_parts_bin` mirrors into the active bin (sacred contract:
    the multi-bin record never drifts from `parts_bin.json`).
  * `_sync_active_parts_bin_parts` is a silent no-op when there is
    no active bin (first-launch race).
"""
from __future__ import annotations

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# Round-trip + cache hygiene
# ═══════════════════════════════════════════════════════════════════════════════

class TestPartsBinCollectionsRoundTrip:
    """`_save_parts_bin_collections` / `_load_parts_bin_collections`
    preserve the full schema (including unknown forward-compat keys)."""

    def test_empty_load(self):
        # Fresh `_PARTS_BIN_COLLECTIONS_FILE` (autouse conftest
        # redirects to tmp) returns an empty list, no exception.
        assert sc._load_parts_bin_collections() == []

    def test_round_trip_preserves_fields(self):
        entries = [{
            "name": "Yeast toolkit",
            "description": "Saccharomyces parts",
            "parts": [{"name": "pYeast-1", "level": 0}],
            "saved": "2026-05-12",
            "_plugin_data": {"some_plugin": {"x": 1}},  # reserved field
        }]
        sc._save_parts_bin_collections(entries)
        loaded = sc._load_parts_bin_collections()
        assert loaded == entries

    def test_load_deepcopies(self):
        """Per invariant #17 — caller mutations after load must not
        poison the in-memory cache."""
        sc._save_parts_bin_collections([
            {"name": "A", "description": "", "parts": [], "saved": ""},
        ])
        first = sc._load_parts_bin_collections()
        first[0]["name"] = "MUTATED"
        second = sc._load_parts_bin_collections()
        assert second[0]["name"] == "A"

    def test_save_deepcopies(self):
        """Per invariant #17 — caller mutations after save must not
        leak into the next load via shared dict refs."""
        entries = [{"name": "A", "description": "", "parts": [], "saved": ""}]
        sc._save_parts_bin_collections(entries)
        entries[0]["name"] = "MUTATED"
        loaded = sc._load_parts_bin_collections()
        assert loaded[0]["name"] == "A"


# ═══════════════════════════════════════════════════════════════════════════════
# Active-bin pointer
# ═══════════════════════════════════════════════════════════════════════════════

class TestActiveBinPointer:
    def test_initial_value_is_none(self):
        assert sc._get_active_parts_bin_name() is None

    def test_set_and_get(self):
        sc._set_active_parts_bin_name("Yeast toolkit")
        assert sc._get_active_parts_bin_name() == "Yeast toolkit"

    def test_clear(self):
        sc._set_active_parts_bin_name("A")
        sc._set_active_parts_bin_name(None)
        assert sc._get_active_parts_bin_name() is None

    def test_find_bin(self):
        sc._save_parts_bin_collections([
            {"name": "A", "description": "", "parts": [], "saved": ""},
            {"name": "B", "description": "", "parts": [], "saved": ""},
        ])
        assert sc._find_parts_bin("A")["name"] == "A"
        assert sc._find_parts_bin("B")["name"] == "B"
        assert sc._find_parts_bin("C") is None

    def test_name_taken(self):
        sc._save_parts_bin_collections([
            {"name": "A", "description": "", "parts": [], "saved": ""},
        ])
        assert sc._parts_bin_name_taken("A") is True
        assert sc._parts_bin_name_taken("B") is False

    def test_switch_active_parts_bin(self):
        # Shared bin-switch helper (used by the bins picker's Open AND
        # the Domesticator store dialog): flips the active pointer +
        # mirrors the chosen bin's parts; safe no-op for a missing bin;
        # the outgoing bin's parts are never lost.
        sc._save_parts_bin_collections([
            {"name": "A", "description": "", "saved": "",
             "parts": [{"name": "pa", "type": "CDS", "sequence": "ATG"}]},
            {"name": "B", "description": "", "saved": "",
             "parts": [{"name": "pb", "type": "Promoter", "sequence": "TTG"}]},
        ])
        sc._set_active_parts_bin_name("A")
        sc._switch_active_parts_bin("A")          # prime the active mirror
        assert sc._switch_active_parts_bin("B") is True
        assert sc._get_active_parts_bin_name() == "B"
        sc._parts_bin_cache = None
        assert {p.get("name") for p in sc._load_parts_bin()} == {"pb"}
        colls = {b["name"]: b for b in sc._load_parts_bin_collections()}
        assert {p.get("name") for p in colls["A"]["parts"]} == {"pa"}
        assert {p.get("name") for p in colls["B"]["parts"]} == {"pb"}
        # Missing bin → no-op, active unchanged.
        assert sc._switch_active_parts_bin("Ghost") is False
        assert sc._get_active_parts_bin_name() == "B"


# ═══════════════════════════════════════════════════════════════════════════════
# Migration
# ═══════════════════════════════════════════════════════════════════════════════

class TestEnsureDefaultPartsBin:
    """First-run migration: wrap existing `parts_bin.json` into a
    "Main Parts Bin" wrapper. Idempotent on subsequent calls."""

    def test_empty_first_run_creates_empty_main_bin(self):
        # No parts_bin.json contents, no parts_bin_collections.json yet.
        assert sc._load_parts_bin() == []
        assert sc._load_parts_bin_collections() == []
        sc._ensure_default_parts_bin()
        bins = sc._load_parts_bin_collections()
        assert len(bins) == 1
        assert bins[0]["name"] == sc._DEFAULT_PARTS_BIN_NAME
        assert bins[0]["parts"] == []
        assert sc._get_active_parts_bin_name() == sc._DEFAULT_PARTS_BIN_NAME

    def test_existing_parts_wrap_into_main_bin(self):
        # Seed `parts_bin.json` with two parts BEFORE migration runs.
        sc._save_parts_bin([
            {"name": "p1", "level": 0, "sequence": "ATGC"},
            {"name": "p2", "level": 1, "sequence": "GGGG"},
        ])
        assert sc._load_parts_bin_collections() == []   # not yet migrated
        sc._ensure_default_parts_bin()
        bins = sc._load_parts_bin_collections()
        assert len(bins) == 1
        assert bins[0]["name"] == sc._DEFAULT_PARTS_BIN_NAME
        assert len(bins[0]["parts"]) == 2
        assert bins[0]["parts"][0]["name"] == "p1"

    def test_idempotent_on_subsequent_calls(self):
        sc._ensure_default_parts_bin()
        before = sc._load_parts_bin_collections()
        sc._ensure_default_parts_bin()
        sc._ensure_default_parts_bin()
        after = sc._load_parts_bin_collections()
        assert before == after

    def test_restores_active_pointer_if_lost(self):
        # Bins exist on disk but the active-bin setting was cleared
        # (could happen if a user hand-edited settings.json).
        sc._save_parts_bin_collections([
            {"name": "Custom", "description": "", "parts": [], "saved": ""},
        ])
        sc._set_active_parts_bin_name(None)
        sc._ensure_default_parts_bin()
        # Should adopt the first existing bin's name, not blow away
        # the user's custom bin.
        assert sc._get_active_parts_bin_name() == "Custom"
        assert len(sc._load_parts_bin_collections()) == 1


# ═══════════════════════════════════════════════════════════════════════════════
# Mirror sync (sacred contract: every `_save_parts_bin` updates the bin)
# ═══════════════════════════════════════════════════════════════════════════════

class TestMirrorSync:
    def test_save_parts_bin_mirrors_into_active(self):
        sc._ensure_default_parts_bin()   # creates Main Parts Bin + sets active
        assert sc._get_active_parts_bin_name() == sc._DEFAULT_PARTS_BIN_NAME
        sc._save_parts_bin([
            {"name": "p1", "level": 0, "sequence": "ATGC"},
        ])
        bins = sc._load_parts_bin_collections()
        main = next(b for b in bins if b["name"] == sc._DEFAULT_PARTS_BIN_NAME)
        assert len(main["parts"]) == 1
        assert main["parts"][0]["name"] == "p1"

    def test_save_with_no_active_bin_is_silent_noop(self):
        """First-launch race: `_save_parts_bin` may fire before
        `_ensure_default_parts_bin` runs (test fixtures, agent API).
        The mirror call should silently do nothing rather than crash
        or create a phantom bin."""
        sc._set_active_parts_bin_name(None)
        assert sc._load_parts_bin_collections() == []
        # This must not raise.
        sc._save_parts_bin([{"name": "p1"}])
        # And it must not have created a bin out of nowhere.
        assert sc._load_parts_bin_collections() == []
        # The parts file itself was still written (legacy behaviour).
        assert sc._load_parts_bin() == [{"name": "p1"}]

    def test_save_mirrors_into_only_the_active_bin(self):
        # Two bins; only the active one should pick up the new parts.
        sc._save_parts_bin_collections([
            {"name": "A", "description": "", "parts": [], "saved": ""},
            {"name": "B", "description": "",
             "parts": [{"name": "kept"}], "saved": ""},
        ])
        sc._set_active_parts_bin_name("A")
        sc._save_parts_bin([{"name": "new-in-A"}])
        bins = sc._load_parts_bin_collections()
        a = next(b for b in bins if b["name"] == "A")
        b = next(b for b in bins if b["name"] == "B")
        assert a["parts"] == [{"name": "new-in-A"}]
        assert b["parts"] == [{"name": "kept"}]   # untouched

    def test_save_mirror_skips_when_active_name_was_deleted(self):
        """User deletes the active bin out from under us (e.g. via
        the picker concurrently with a save). The mirror should
        silently skip rather than recreate the deleted bin."""
        sc._save_parts_bin_collections([
            {"name": "A", "description": "", "parts": [], "saved": ""},
        ])
        sc._set_active_parts_bin_name("Ghost")   # points at a deleted bin
        sc._save_parts_bin([{"name": "p1"}])
        bins = sc._load_parts_bin_collections()
        # No "Ghost" entry materialised.
        assert {b["name"] for b in bins} == {"A"}
