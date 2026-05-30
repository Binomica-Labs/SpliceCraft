"""Enzyme collections + custom enzymes + Settings modal (2026-05-22).

Covers:
* Persistence layer (`_load/_save_enzyme_collections`,
  `_load/_save_custom_enzymes`, `_all_enzymes()` combined accessor).
* One-shot CSV → collection migration
  (`_migrate_legacy_custom_enzyme_csv`).
* Restriction-scan integration: `_active_enzyme_allowed_set` filters
  against `_all_enzymes()` so custom enzymes participate.
* `EnzymeCollectionsModal` smoke: two right-pane modes
  (list ↔ editor), `+ Add new enzyme…` push, search field default.
* `SettingsModal` smoke: every toggle wired, sub-modal launcher
  buttons present.
"""
import pytest

import splicecraft as sc

pytestmark = [pytest.mark.usefixtures("_protect_user_data")]


# ── Persistence layer ──────────────────────────────────────────────────────

def test_enzyme_collections_round_trip():
    sc._save_enzyme_collections([
        {"name": "Common cloners", "enzymes": ["EcoRI", "BamHI"]},
    ])
    out = sc._load_enzyme_collections()
    assert len(out) == 1
    assert out[0]["name"] == "Common cloners"
    assert out[0]["enzymes"] == ["EcoRI", "BamHI"]


def test_find_enzyme_collection_missing_returns_none():
    sc._save_enzyme_collections([{"name": "A", "enzymes": []}])
    assert sc._find_enzyme_collection("nope") is None
    assert sc._find_enzyme_collection("A") is not None


def test_active_enzyme_collection_pointer():
    assert sc._get_active_enzyme_collection_name() is None
    sc._set_active_enzyme_collection_name("Foo")
    assert sc._get_active_enzyme_collection_name() == "Foo"
    sc._set_active_enzyme_collection_name(None)
    assert sc._get_active_enzyme_collection_name() is None


def test_custom_enzymes_round_trip():
    payload = {
        "name": "MyEnzI",
        "site": "GGTACC",
        "fwd_cut": 1,
        "rev_cut": 5,
        "type": "II_5overhang",
        "supplier": "Lab fridge",
    }
    sc._save_custom_enzymes([payload])
    out = sc._load_custom_enzymes()
    assert out == [payload]
    meta = sc._custom_enzyme_meta("MyEnzI")
    assert meta is not None and meta["supplier"] == "Lab fridge"


def test_all_enzymes_includes_custom():
    sc._save_custom_enzymes([
        {"name": "TestUnique1", "site": "AAATTT",
         "fwd_cut": 3, "rev_cut": 3, "type": "II_blunt", "supplier": ""},
    ])
    combined = sc._all_enzymes()
    assert "TestUnique1" in combined
    assert combined["TestUnique1"] == ("AAATTT", 3, 3)
    # Built-in enzymes still present.
    assert "EcoRI" in combined


def test_custom_overrides_builtin_on_name_collision():
    """User-added enzyme with same name as built-in overrides — gives
    the user the last word on cut definitions."""
    # EcoRI normally cuts G^AATTC (fwd_cut=1).
    sc._save_custom_enzymes([
        {"name": "EcoRI", "site": "GAATTC",
         "fwd_cut": 4, "rev_cut": 2, "type": "II_5overhang", "supplier": ""},
    ])
    combined = sc._all_enzymes()
    assert combined["EcoRI"] == ("GAATTC", 4, 2)


# ── Restriction-scan integration ──────────────────────────────────────────

def test_active_enzyme_allowed_set_none_when_no_active():
    assert sc._active_enzyme_allowed_set() is None


def test_active_enzyme_allowed_set_resolves_to_known_enzymes():
    sc._save_enzyme_collections([
        {"name": "TwoCutter",
         "enzymes": ["EcoRI", "BamHI", "BogusName123"]},
    ])
    sc._set_active_enzyme_collection_name("TwoCutter")
    allowed = sc._active_enzyme_allowed_set()
    # Bogus name dropped (not in master); known names kept.
    assert allowed == frozenset({"EcoRI", "BamHI"})


def test_active_enzyme_allowed_set_includes_custom():
    sc._save_custom_enzymes([
        {"name": "TestUnique1", "site": "AAATTT",
         "fwd_cut": 3, "rev_cut": 3, "type": "II_blunt", "supplier": ""},
    ])
    sc._save_enzyme_collections([
        {"name": "WithCustom",
         "enzymes": ["EcoRI", "TestUnique1"]},
    ])
    sc._set_active_enzyme_collection_name("WithCustom")
    allowed = sc._active_enzyme_allowed_set()
    assert allowed == frozenset({"EcoRI", "TestUnique1"})


# ── CSV → collection migration ────────────────────────────────────────────

def test_migrate_legacy_csv_empty_is_noop():
    sc._set_setting("restr_custom_enzymes", "")
    sc._migrate_legacy_custom_enzyme_csv()
    assert sc._load_enzyme_collections() == []


def test_migrate_legacy_csv_creates_collection():
    sc._set_setting("restr_custom_enzymes", "EcoRI, BamHI, HindIII")
    sc._set_setting("restr_use_custom_list", True)
    sc._migrate_legacy_custom_enzyme_csv()
    colls = sc._load_enzyme_collections()
    assert len(colls) == 1
    assert colls[0]["name"] == "Custom (legacy)"
    assert sorted(colls[0]["enzymes"]) == ["BamHI", "EcoRI", "HindIII"]
    # was_active=True → migrated collection should be the active one.
    assert sc._get_active_enzyme_collection_name() == "Custom (legacy)"
    # Legacy settings cleared.
    assert sc._get_setting("restr_custom_enzymes", "") == ""
    assert sc._get_setting("restr_use_custom_list", False) is False


def test_migrate_legacy_csv_idempotent():
    """Re-running the migration must not create a duplicate."""
    sc._set_setting("restr_custom_enzymes", "EcoRI")
    sc._migrate_legacy_custom_enzyme_csv()
    sc._set_setting("restr_custom_enzymes", "EcoRI")  # simulate re-run
    sc._migrate_legacy_custom_enzyme_csv()
    colls = sc._load_enzyme_collections()
    assert len(colls) == 1


def test_migrate_legacy_csv_unknown_names_only_clears_settings():
    """Pre-fix safeguard: a CSV holding only unknown names shouldn't
    create an empty `Custom (legacy)` row."""
    sc._set_setting("restr_custom_enzymes", "BogusA, BogusB")
    sc._migrate_legacy_custom_enzyme_csv()
    assert sc._load_enzyme_collections() == []
    # But the legacy settings should still be cleared.
    assert sc._get_setting("restr_custom_enzymes", "") == ""


# ── Modal smoke tests ─────────────────────────────────────────────────────

@pytest.mark.asyncio
async def test_enzyme_collections_modal_mounts():
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        app.push_screen(sc.EnzymeCollectionsModal())
        await pilot.pause()
        await pilot.pause()
        # Master table populated with the full catalog.
        from textual.widgets import DataTable
        master = pilot.app.screen.query_one(
            "#ec-master-table", DataTable,
        )
        assert master.row_count > 0
        # Catalogs table in list mode is reachable too.
        catalogs = pilot.app.screen.query_one(
            "#ec-catalogs-table", DataTable,
        )
        assert catalogs is not None


@pytest.mark.asyncio
async def test_enzyme_collections_modal_has_both_tabs():
    """The modal exposes two tabs: Enzyme Sets (catalog manager) +
    Enzyme Settings (the restriction-overlay toggles that used to
    live in the Enzymes dropdown). Both TabPanes must mount."""
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        app.push_screen(sc.EnzymeCollectionsModal())
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import TabPane, Checkbox, RadioSet
        screen = pilot.app.screen
        # Both panes mount.
        sets = screen.query_one("#tab-sets", TabPane)
        settings = screen.query_one("#tab-settings", TabPane)
        assert sets is not None and settings is not None
        # Settings-tab widgets all present.
        assert screen.query_one("#ec-set-show-restr", Checkbox) is not None
        assert screen.query_one("#ec-set-unique", Checkbox) is not None
        assert screen.query_one("#ec-set-min-len", RadioSet) is not None
        assert screen.query_one("#ec-set-connectors", Checkbox) is not None


@pytest.mark.asyncio
async def test_enzyme_collections_modal_opens_in_list_mode_when_no_active():
    """Default landing is the catalog list view (Mode A)."""
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.EnzymeCollectionsModal()
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        assert modal._mode == sc.EnzymeCollectionsModal._MODE_LIST


@pytest.mark.asyncio
async def test_enzyme_collections_modal_opens_in_editor_when_active_set():
    """If an active catalog exists, the modal opens in editor mode."""
    sc._save_enzyme_collections([
        {"name": "PreActive", "enzymes": ["EcoRI"]},
    ])
    sc._set_active_enzyme_collection_name("PreActive")
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.EnzymeCollectionsModal()
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        assert modal._mode == sc.EnzymeCollectionsModal._MODE_EDITOR
        assert modal._editor_catalog == "PreActive"


@pytest.mark.asyncio
async def test_add_custom_enzyme_modal_persists():
    """AddCustomEnzymeModal saves to custom_enzymes.json and the
    new enzyme appears in `_all_enzymes()` immediately."""
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.AddCustomEnzymeModal()
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import Input, Button
        modal.query_one("#ace-name", Input).value = "AceTestEnz"
        modal.query_one("#ace-site", Input).value = "ACGTAC"
        modal.query_one("#ace-fwd",  Input).value = "2"
        modal.query_one("#ace-rev",  Input).value = "4"
        await pilot.pause()
        modal.query_one("#ace-save", Button).press()
        await pilot.pause()
        await pilot.pause()
    assert "AceTestEnz" in sc._all_enzymes()
    assert sc._all_enzymes()["AceTestEnz"] == ("ACGTAC", 2, 4)


@pytest.mark.asyncio
async def test_add_custom_enzyme_rejects_name_collision_with_builtin():
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.AddCustomEnzymeModal()
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import Input, Button
        modal.query_one("#ace-name", Input).value = "EcoRI"  # built-in!
        modal.query_one("#ace-site", Input).value = "GAATTC"
        modal.query_one("#ace-fwd",  Input).value = "1"
        modal.query_one("#ace-rev",  Input).value = "5"
        await pilot.pause()
        modal.query_one("#ace-save", Button).press()
        await pilot.pause()
        # Modal is still mounted (rejected) — built-in collision must
        # NOT have been saved into custom_enzymes.json.
        assert modal.is_mounted
        saved = sc._load_custom_enzymes()
        assert not any(
            e.get("name") == "EcoRI" for e in saved
        ), "name-collision check failed — EcoRI got persisted"


@pytest.mark.asyncio
async def test_settings_modal_renders_every_group():
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.SettingsModal()
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import Checkbox, Button, Input
        # Every toggle widget present.
        for cid in (
            "#set-tooltips",
            "#set-click-debug",
            "#set-check-updates",
            "#set-constr-filter",
        ):
            assert modal.query_one(cid, Checkbox) is not None
        # Min primer binding numeric input + apply button.
        assert modal.query_one("#set-min-primer", Input) is not None
        assert modal.query_one("#set-min-primer-apply", Button) is not None
        # Sub-modal launchers — all five buttons.
        for bid in (
            "#set-grammars",
            "#set-entry-vectors",
            "#set-enzyme-collections",
            "#set-codon-tables",
            "#set-restore",
        ):
            assert modal.query_one(bid, Button) is not None


def test_settings_modal_blocks_undo():
    """SettingsModal mutates persistent state — stray Ctrl+Z on the
    canvas underneath must not race the save."""
    assert sc.SettingsModal._blocks_undo is True


def test_enzyme_collections_modal_blocks_undo():
    """Same reasoning as SettingsModal."""
    assert sc.EnzymeCollectionsModal._blocks_undo is True


def test_add_custom_enzyme_modal_blocks_undo():
    assert sc.AddCustomEnzymeModal._blocks_undo is True


# ── Menu wiring ──────────────────────────────────────────────────────────

def test_settings_menubar_opens_modal_directly():
    """Clicking the Settings menubar item opens `SettingsModal`
    straight away — no dropdown. The dropdown collapsed to a single
    entry first (sweep #24 landing) then was promoted to direct-open
    when the Enzymes menu was also direct-opened. Greppable so a
    refactor that re-introduces a Settings dropdown trips this test."""
    with open(sc.__file__, encoding="utf-8") as f:
        src = f.read()
    # Action exists.
    assert "def action_open_settings" in src
    # MenuBar.on_click routes Settings directly to the action.
    assert "self.app.action_open_settings()" in src
    # No dropdown items registered for "Settings".
    assert '"Settings": [' not in src


def test_enzymes_menubar_opens_modal_directly():
    """Same pattern as Settings — the Enzymes menubar item opens
    `EnzymeCollectionsModal` directly. The radio toggles that used to
    live in the dropdown moved to the modal's `Enzyme Settings` tab."""
    with open(sc.__file__, encoding="utf-8") as f:
        src = f.read()
    assert "def action_open_enzyme_collections" in src
    assert "self.app.action_open_enzyme_collections()" in src
    # No dropdown items registered for "Enzymes".
    assert '"Enzymes": [' not in src
    # Legacy custom-enzyme-list entry stays unwired.
    assert '"edit_custom_enzyme_list"' not in src


# ── Master-delete wiring ──────────────────────────────────────────────────

def test_enzyme_collections_in_user_data_file_attrs():
    """New file constants must appear in `_USER_DATA_FILE_ATTRS` so
    Master Delete / pre-update snapshots / restore UI cover them."""
    assert "_ENZYME_COLLECTIONS_FILE" in sc._USER_DATA_FILE_ATTRS
    assert "_CUSTOM_ENZYMES_FILE" in sc._USER_DATA_FILE_ATTRS


def test_enzyme_caches_in_master_delete_attrs():
    assert "_enzyme_collections_cache" in sc._MASTER_DELETE_CACHE_ATTRS
    assert "_custom_enzymes_cache" in sc._MASTER_DELETE_CACHE_ATTRS


def test_restore_modal_targets_include_enzyme_files():
    targets = dict(sc.RestoreFromBackupModal._TARGETS)
    assert "_ENZYME_COLLECTIONS_FILE" in targets.values()
    assert "_CUSTOM_ENZYMES_FILE" in targets.values()


def test_agent_backup_labels_include_enzyme_files():
    """[INV-64] sacred contract — agent restore-backup parity with
    `_USER_DATA_FILE_ATTRS`. Without these the agent API can't list
    or restore the two new files even though the GUI can."""
    labels = sc._AGENT_BACKUP_LABELS
    assert labels.get("custom_enzymes") == "_CUSTOM_ENZYMES_FILE"
    assert labels.get("enzyme_collections") == "_ENZYME_COLLECTIONS_FILE"


# ── JSON envelope + corruption recovery ──────────────────────────────────

def test_enzyme_collections_writes_schema_envelope():
    """Sacred invariant #7 — every persistent JSON file lands as
    ``{"_schema_version": 1, "entries": [...]}``."""
    import json
    sc._save_enzyme_collections([{"name": "Foo", "enzymes": ["EcoRI"]}])
    raw = sc._ENZYME_COLLECTIONS_FILE.read_text(encoding="utf-8")
    obj = json.loads(raw)
    assert isinstance(obj, dict)
    assert obj.get("_schema_version") == 1
    assert isinstance(obj.get("entries"), list)


def test_custom_enzymes_writes_schema_envelope():
    import json
    sc._save_custom_enzymes([
        {"name": "EnvX", "site": "GAATTC",
         "fwd_cut": 1, "rev_cut": 5, "type": "II_5overhang", "supplier": ""},
    ])
    obj = json.loads(sc._CUSTOM_ENZYMES_FILE.read_text(encoding="utf-8"))
    assert obj.get("_schema_version") == 1
    assert obj["entries"][0]["name"] == "EnvX"


def test_enzyme_collections_accepts_legacy_bare_list():
    """Pre-0.3.1 back-compat — `_extract_entries` must accept a bare
    list payload (no envelope)."""
    import json
    sc._ENZYME_COLLECTIONS_FILE.write_text(
        json.dumps([{"name": "Legacy", "enzymes": ["BamHI"]}]),
        encoding="utf-8",
    )
    # Bust the cache so the next load re-reads from disk.
    sc._enzyme_collections_cache = None
    out = sc._load_enzyme_collections()
    assert out == [{"name": "Legacy", "enzymes": ["BamHI"]}]


def test_enzyme_collections_corrupted_json_returns_empty():
    """`_safe_load_json` recovers from corruption by returning ``None``
    so `_load_enzyme_collections` falls back to ``[]``."""
    sc._ENZYME_COLLECTIONS_FILE.write_text("{not-json", encoding="utf-8")
    sc._enzyme_collections_cache = None
    out = sc._load_enzyme_collections()
    assert out == []


def test_save_enzyme_collections_propagates_on_failure(monkeypatch):
    """Sacred invariant #7 — `_safe_save_json` re-raises so callers
    can notify. Without re-raise the EnzymeCollectionsModal would
    silently desync UI from disk."""
    def boom(*args, **kwargs):
        raise OSError("disk full (synthetic)")
    monkeypatch.setattr(sc, "_safe_save_json", boom)
    with pytest.raises(OSError, match="disk full"):
        sc._save_enzyme_collections([{"name": "X", "enzymes": []}])


# ── AddCustomEnzymeModal validation ──────────────────────────────────────

@pytest.mark.asyncio
async def test_add_custom_enzyme_rejects_non_iupac():
    """Recognition site with non-IUPAC characters must be rejected
    BEFORE persistence — otherwise `_iupac_pattern` would build a
    bogus regex at scan time."""
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.AddCustomEnzymeModal()
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import Input, Button
        modal.query_one("#ace-name", Input).value = "BadIUPAC"
        # 'Z' and '*' are not in the IUPAC alphabet.
        modal.query_one("#ace-site", Input).value = "GAAZTC"
        modal.query_one("#ace-fwd",  Input).value = "1"
        modal.query_one("#ace-rev",  Input).value = "5"
        await pilot.pause()
        modal.query_one("#ace-save", Button).press()
        await pilot.pause()
        assert modal.is_mounted, "modal should not dismiss on invalid IUPAC"
        assert not any(
            e.get("name") == "BadIUPAC"
            for e in sc._load_custom_enzymes()
        )


@pytest.mark.asyncio
async def test_add_custom_enzyme_rejects_cut_position_out_of_range():
    app = sc.PlasmidApp()
    async with app.run_test(size=(160, 48)) as pilot:
        await pilot.pause()
        modal = sc.AddCustomEnzymeModal()
        app.push_screen(modal)
        await pilot.pause()
        await pilot.pause()
        from textual.widgets import Input, Button
        modal.query_one("#ace-name", Input).value = "OutOfRange"
        modal.query_one("#ace-site", Input).value = "GAATTC"
        # site length 6 → allowed range -30..36. 9999 is far outside.
        modal.query_one("#ace-fwd",  Input).value = "9999"
        modal.query_one("#ace-rev",  Input).value = "5"
        await pilot.pause()
        modal.query_one("#ace-save", Button).press()
        await pilot.pause()
        assert modal.is_mounted
        assert not any(
            e.get("name") == "OutOfRange"
            for e in sc._load_custom_enzymes()
        )


# ── Active-collection edge cases ─────────────────────────────────────────

def test_active_collection_returns_none_when_collection_deleted():
    """Pre-fix: if the active-collection pointer references a deleted
    collection, `_active_enzyme_allowed_set` must return None rather
    than raising or returning an empty set (which would silently scan
    NOTHING)."""
    sc._save_enzyme_collections([
        {"name": "WillDelete", "enzymes": ["EcoRI"]},
    ])
    sc._set_active_enzyme_collection_name("WillDelete")
    # Simulate a hand-edited collections.json that removed the row
    # without clearing the active pointer.
    sc._save_enzyme_collections([])
    assert sc._get_active_enzyme_collection_name() == "WillDelete"
    allowed = sc._active_enzyme_allowed_set()
    assert allowed is None


# ── Restriction-scan integration with custom enzymes ─────────────────────

def test_scan_finds_custom_enzyme_hit():
    """End-to-end: after `_save_custom_enzymes`, the scanner must
    surface hits for the new enzyme. Pre-fix the scanner read
    `_NEB_ENZYMES` directly and skipped every custom enzyme."""
    sc._save_custom_enzymes([
        {"name": "TestUniqHit", "site": "GGGCCC",
         "fwd_cut": 3, "rev_cut": 3, "type": "II_blunt", "supplier": ""},
    ])
    # Build a sequence with one occurrence of GGGCCC.
    seq = "AAAAAA" + "GGGCCC" + "AAAAAA"
    cuts = sc._enzyme_cuts(seq, ["TestUniqHit"], circular=False)
    assert len(cuts) == 1
    assert cuts[0].get("enzyme") == "TestUniqHit"


def test_scan_filter_via_active_collection_includes_custom():
    """`_scan_restriction_sites` with an allowed-set containing a
    custom enzyme must return that enzyme's hits."""
    sc._save_custom_enzymes([
        {"name": "TestUniqAllow", "site": "TTAATT",
         "fwd_cut": 3, "rev_cut": 3, "type": "II_blunt", "supplier": ""},
    ])
    seq = "CCCC" + "TTAATT" + "CCCC"
    sites = sc._scan_restriction_sites(
        seq, circular=False,
        allowed_enzymes=frozenset({"TestUniqAllow"}),
    )
    assert any(s.get("label") == "TestUniqAllow" for s in sites
                if s.get("type") == "resite")
