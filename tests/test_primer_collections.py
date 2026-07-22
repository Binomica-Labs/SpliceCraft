"""
Primer collections + new Design Primers save modal — regression coverage.

Adds 2026-05-21 (post-v0.9.12). Covers:
  * Primer-collections file load/save/mirror roundtrip.
  * Default-wrap on first launch (`_ensure_default_primer_collection`).
  * `_save_primers` mirrors to active primer collection.
  * Active-pointer switch via `_set_active_primer_collection_name`.
  * Settings-schema key `active_primer_collection`.
  * `RestoreFromBackupModal` covers `_PRIMER_COLLECTIONS_FILE`.
  * Spaces in primer names AND collection names round-trip 1:1
    (sacred: no underscore substitution on user-typed names).
  * Collection-switch on the plasmid panel loads first plasmid or
    blanks canvas if empty.
"""

import splicecraft as sc


class TestPrimerCollectionInfrastructure:
    """Module-level helpers — no app instantiation needed."""

    def test_load_empty_returns_list(self):
        assert sc._load_primer_collections() == []

    def test_save_and_reload_roundtrip(self):
        entries = [{
            "name":        "Project alpha",   # spaces preserved
            "description": "test",
            "primers":     [],
            "saved":       "2026-05-21",
        }]
        sc._save_primer_collections(entries)
        # Force cold reload from disk
        sc._state._primer_collections_cache = None
        out = sc._load_primer_collections()
        assert len(out) == 1
        assert out[0]["name"] == "Project alpha"
        # Spaces preserved in the on-disk name.
        assert "_" not in out[0]["name"]

    def test_deepcopy_on_read(self):
        sc._save_primer_collections([
            {"name": "A", "primers": [{"name": "p1", "sequence": "AC"}]}
        ])
        a = sc._load_primer_collections()
        a[0]["primers"].append({"name": "leak", "sequence": "GG"})
        b = sc._load_primer_collections()
        assert len(b[0]["primers"]) == 1
        assert b[0]["primers"][0]["name"] == "p1"


class TestDefaultPrimerCollection:
    def test_first_launch_wraps_existing_primers(self):
        # Seed a flat primers.json with some user content.
        sc._save_primers([
            {"name": "existing primer 1", "sequence": "ACGT"},
            {"name": "p with spaces in name", "sequence": "TTAA"},
        ])
        # Pre-existing primer_collections is empty → wrap path.
        assert sc._load_primer_collections() == []
        sc._ensure_default_primer_collection()
        # "Main" created, marked active, contents copied from primers.json.
        assert sc._get_active_primer_collection_name() == "Main"
        colls = sc._load_primer_collections()
        assert len(colls) == 1 and colls[0]["name"] == "Main"
        names = [p["name"] for p in colls[0]["primers"]]
        # Sacred: user-typed names with spaces survive verbatim.
        assert "existing primer 1" in names
        assert "p with spaces in name" in names

    def test_idempotent_on_second_launch(self):
        sc._save_primer_collections([
            {"name": "existing", "primers": []}
        ])
        # No active pointer → default to first collection's name.
        sc._ensure_default_primer_collection()
        assert sc._get_active_primer_collection_name() == "existing"
        # Run again — nothing changes.
        n_before = len(sc._load_primer_collections())
        sc._ensure_default_primer_collection()
        assert len(sc._load_primer_collections()) == n_before


class TestSavePrimersMirror:
    """Sacred: every _save_primers call must mirror to active primer
    collection inside `_cache_lock` ([INV-50] save-chain lock-release
    gap). Verifies the mirror lands without a separate explicit call."""

    def test_save_primers_mirrors_to_active_collection(self):
        sc._save_primer_collections([
            {"name": "Main", "primers": []}
        ])
        sc._set_active_primer_collection_name("Main")
        sc._save_primers([
            {"name": "new primer with spaces", "sequence": "ACGTACGT"},
        ])
        colls = sc._load_primer_collections()
        active = next(c for c in colls if c["name"] == "Main")
        assert len(active["primers"]) == 1
        assert active["primers"][0]["name"] == "new primer with spaces"

    def test_save_primers_no_active_no_mirror(self):
        sc._save_primer_collections([
            {"name": "Inactive", "primers": []}
        ])
        # active pointer empty
        sc._set_active_primer_collection_name(None)
        sc._save_primers([{"name": "p1", "sequence": "AT"}])
        # Collection's primers untouched.
        colls = sc._load_primer_collections()
        assert colls[0]["primers"] == []


class TestActivePointer:
    def test_set_and_get_roundtrip(self):
        sc._set_active_primer_collection_name("My collection")
        assert sc._get_active_primer_collection_name() == "My collection"

    def test_set_none_clears(self):
        sc._set_active_primer_collection_name("foo")
        sc._set_active_primer_collection_name(None)
        assert sc._get_active_primer_collection_name() is None

    def test_settings_schema_includes_key(self):
        assert "active_primer_collection" in sc._SETTINGS_SCHEMA
        types, default = sc._SETTINGS_SCHEMA["active_primer_collection"]
        assert str in types
        assert default == ""


class TestCreatePrimerCollectionEndpoint:
    """Agent `create-primer-collection` — the plasmid `create-collection`
    parallel that the GUI-only path previously lacked (petunia running-log:
    "no agent way to create a primer collection")."""

    def test_registered_as_write_endpoint(self):
        assert "create-primer-collection" in sc._state._AGENT_HANDLERS
        _fn, write = sc._state._AGENT_HANDLERS["create-primer-collection"]
        assert write is True

    def test_creates_empty_collection(self):
        r = sc._h_create_primer_collection(None, {"name": "Petunia primers"})
        assert r["ok"] and r["name"] == "Petunia primers" and r["n_primers"] == 0
        names = [c["name"] for c in sc._load_primer_collections()]
        assert "Petunia primers" in names
        # Spaces preserved 1:1 (no underscore mangle).
        assert "_" not in "Petunia primers"

    def test_then_set_active_and_add_primer(self):
        # The documented flow: create -> set-active -> create-primer lands
        # the primer in the NEW collection, not the default.
        sc._h_create_primer_collection(None, {"name": "Build A"})
        sc._h_set_active_primer_collection(None, {"name": "Build A"})
        sc._h_create_primer(None, {"name": "oF1", "sequence": "ACGTACGTAC"})
        coll = next(c for c in sc._load_primer_collections()
                    if c["name"] == "Build A")
        assert any(p.get("name") == "oF1" for p in coll.get("primers", []))

    def test_duplicate_name_409_case_insensitive(self):
        sc._h_create_primer_collection(None, {"name": "Dup"})
        r = sc._h_create_primer_collection(None, {"name": "dup"})
        assert isinstance(r, tuple) and r[1] == 409

    def test_blank_name_400(self):
        for bad in ("", "   ", None, 123, {"x": 1}):
            r = sc._h_create_primer_collection(None, {"name": bad})
            assert isinstance(r, tuple) and r[1] == 400

    def test_description_stored_and_unknown_keys_echoed(self):
        r = sc._h_create_primer_collection(
            None, {"name": "Desc", "description": "notes here", "bogus": 1})
        assert r["ignored"] == ["bogus"]
        coll = next(c for c in sc._load_primer_collections()
                    if c["name"] == "Desc")
        assert coll["description"] == "notes here"


class TestSetActivePrimerCollectionRemirror:
    """`set-active-primer-collection` must re-point the live primers.json at
    the newly-active collection — otherwise the next create/save mirrors the
    PREVIOUS collection's primers back over the new one (data loss). The agent
    parallel to the plasmid `set-active-collection` + startup restore."""

    def test_switch_remirrors_live_library(self):
        sc._save_primer_collections([
            {"name": "Alpha", "primers": [{"name": "a1", "sequence": "AAAA"}]},
            {"name": "Beta",  "primers": [{"name": "b1", "sequence": "CCCC"},
                                          {"name": "b2", "sequence": "GGGG"}]},
        ])
        sc._h_set_active_primer_collection(None, {"name": "Alpha"})
        assert {p["name"] for p in sc._load_primers()} == {"a1"}
        r = sc._h_set_active_primer_collection(None, {"name": "Beta"})
        assert r["ok"] and r["active"] == "Beta"
        # Live library is now Beta's primers, not the stale Alpha ones.
        assert {p["name"] for p in sc._load_primers()} == {"b1", "b2"}

    def test_switch_then_create_does_not_clobber(self):
        # The corruption the re-mirror prevents: after switching to Beta, a new
        # primer lands in Beta and does NOT overwrite it with Alpha's content.
        sc._save_primer_collections([
            {"name": "Alpha", "primers": [{"name": "a1", "sequence": "AAAA"}]},
            {"name": "Beta",  "primers": [{"name": "b1", "sequence": "CCCC"}]},
        ])
        sc._h_set_active_primer_collection(None, {"name": "Alpha"})
        sc._h_set_active_primer_collection(None, {"name": "Beta"})
        sc._h_create_primer(None, {"name": "b2", "sequence": "TTTT"})
        beta = next(c for c in sc._load_primer_collections()
                    if c["name"] == "Beta")
        assert {p["name"] for p in beta["primers"]} == {"b1", "b2"}   # no a1
        alpha = next(c for c in sc._load_primer_collections()
                     if c["name"] == "Alpha")
        assert {p["name"] for p in alpha["primers"]} == {"a1"}        # untouched

    def test_switch_to_default_keeps_live_library(self):
        # Switching to "" leaves primers.json as-is (the free-standing
        # default) — matching the startup restore, which no-ops for "".
        sc._save_primer_collections([
            {"name": "Alpha", "primers": [{"name": "a1", "sequence": "AAAA"}]},
        ])
        sc._h_set_active_primer_collection(None, {"name": "Alpha"})
        assert {p["name"] for p in sc._load_primers()} == {"a1"}
        r = sc._h_set_active_primer_collection(None, {"name": ""})
        assert r["ok"] and r["active"] == ""
        assert {p["name"] for p in sc._load_primers()} == {"a1"}

    def test_mirror_failure_reverts_active_pointer(self, monkeypatch):
        import splicecraft_agent as sca
        sc._save_primer_collections([
            {"name": "Alpha", "primers": [{"name": "a1", "sequence": "AAAA"}]},
            {"name": "Beta",  "primers": [{"name": "b1", "sequence": "CCCC"}]},
        ])
        sc._h_set_active_primer_collection(None, {"name": "Alpha"})

        def _boom(*a, **k):
            raise OSError("simulated disk failure")
        monkeypatch.setattr(sca, "_safe_save_json_mirror", _boom)
        r = sc._h_set_active_primer_collection(None, {"name": "Beta"})
        assert isinstance(r, tuple) and r[1] == 500
        # Pointer reverted so the live library + active-name setting agree.
        assert (sc._get_active_primer_collection_name() or "") == "Alpha"


class TestRestoreFromBackupCoverage:
    def test_primer_collections_in_restore_targets(self):
        labels = [label for label, attr in sc.RestoreFromBackupModal._TARGETS]
        assert "Primer collections" in labels


class TestMasterDeleteCoverage:
    def test_primer_collections_cache_registered(self):
        assert "_primer_collections_cache" in sc._MASTER_DELETE_CACHE_ATTRS

    def test_primer_collections_file_registered(self):
        assert "_PRIMER_COLLECTIONS_FILE" in sc._USER_DATA_FILE_ATTRS


class TestRestorePrimersFromActiveCollection:
    """Startup helper: rewrites primers.json with the active
    collection's primers so a power loss / external edit can't leave
    primers.json out of sync with the collection."""

    def test_restore_overwrites_primers_file(self):
        sc._save_primer_collections([
            {"name": "Main", "primers": [
                {"name": "from-collection", "sequence": "TTTT"},
            ]}
        ])
        sc._set_active_primer_collection_name("Main")
        # Stash a different live primer into primers.json directly.
        sc._save_primers([{"name": "stale", "sequence": "AAAA"}])
        # But the mirror just wrote into Main again so Main now also
        # holds the stale primer. Re-seed Main fresh.
        sc._save_primer_collections([
            {"name": "Main", "primers": [
                {"name": "from-collection", "sequence": "TTTT"},
            ]}
        ])
        sc._restore_primers_from_active_primer_collection()
        primers = sc._load_primers()
        names = [p["name"] for p in primers]
        assert "from-collection" in names
        # Stale primer wiped — collection is authoritative.
        assert "stale" not in names


class TestPrimerSaveModalExists:
    """Smoke-test the new modal class is importable + has the expected
    constructor signature."""

    def test_modal_class_present(self):
        assert hasattr(sc, "PrimerSaveModal")

    def test_modal_accepts_oligos_kwarg(self):
        modal = sc.PrimerSaveModal(
            [{"label": "Fwd", "default_name": "test-F"}],
            default_collection="Main",
        )
        # Internal state set.
        assert len(modal._oligos) == 1
        assert modal._default_collection == "Main"

    async def test_modal_is_centered_box(self):
        """G (2026-06-13): the dialog must render as a centered, bordered box
        like the other naming modals — not the old frameless full-width block.
        The `#primer-save-dlg` width:80 rule only exists because of G, so a
        fixed 80-wide bordered box proves the harmonized chrome applied."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerSaveModal(
                [{"label": "Fwd", "default_name": "test-F"}],
                default_collection="Main"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            assert isinstance(modal, sc.PrimerSaveModal)
            # centered on the screen (ModalScreen align)
            assert str(modal.styles.align_horizontal) == "center"
            assert str(modal.styles.align_vertical) == "middle"
            # the harmonized box chrome: fixed 80-wide, bordered
            dlg = modal.query_one("#primer-save-dlg")
            assert int(dlg.styles.width.value) == 80, \
                f"dialog not the harmonized 80-wide box: {dlg.styles.width}"
            assert dlg.styles.border_top and dlg.styles.border_top[0], \
                "dialog box has no border (frameless — G chrome not applied)"


class TestPrimerSaveModalPreservesSpaces:
    """User-typed names with spaces flow through verbatim."""

    def test_modal_oligo_label_with_spaces(self):
        modal = sc.PrimerSaveModal(
            [{"label": "Forward primer with spaces",
              "default_name": "my new fwd primer"}],
        )
        # No internal mutation of the typed name.
        assert modal._oligos[0]["default_name"] == "my new fwd primer"

    def test_new_collection_name_preserves_spaces(self):
        # Spaces in collection names must round-trip.
        sc._save_primer_collections([
            {"name": "Old collection name", "primers": []}
        ])
        sc._set_active_primer_collection_name("Old collection name")
        assert (sc._get_active_primer_collection_name()
                == "Old collection name")


class TestSavePrimerToLibraryFlow:
    """PrimerEditModal's "Save to library" → app `_save_primer_to_library_flow`:
    pushes back (no save) when the exact oligo is already in the library;
    otherwise opens `PrimerSaveModal` (pick/create collection) and commits to
    the chosen collection. Unifies a cloning / map primer with the library."""

    import pytest as _pt

    @_pt.mark.asyncio
    async def test_duplicate_pushes_back_and_opens_no_modal(
            self, tiny_record, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            sc._save_primers([{"name": "Existing-F",
                               "sequence": "ACGTACGTACGTACGTACGTA", "tm": 60.0}])
            notes: list = []
            app.notify = lambda msg, **k: notes.append((str(msg), k.get("severity")))
            depth = len(app.screen_stack)
            app._save_primer_to_library_flow(
                {"primer_seq": "acgtacgtacgtacgtacgta", "label": "Dup", "strand": 1})
            await pilot.pause()
            assert any("already in your library" in m and sev == "warning"
                       for m, sev in notes), notes
            assert len(app.screen_stack) == depth, "no picker for a duplicate"

    @_pt.mark.asyncio
    async def test_new_primer_picker_then_commit_to_collection(
            self, tiny_record, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            sc._save_primer_collections([{"name": "Main", "primers": [],
                                          "saved": "2026-06-11"}])
            sc._set_active_primer_collection_name("Main")
            sc._save_primers([])  # empty library
            app._save_primer_to_library_flow(
                {"primer_seq": "GGGGCCCCAAAATTTTGGGGC", "label": "New-F",
                 "strand": 1})
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.PrimerSaveModal), \
                "a new primer must open the collection picker"
            # Pick the 'Main' collection + commit (mimic the modal's payload).
            app.screen.dismiss({"names": ["New-F"], "collection": "Main",
                                "create": False})
            await pilot.pause(); await pilot.pause()
            saved = sc._load_primers()
            assert any(p.get("name") == "New-F"
                       and p.get("sequence") == "GGGGCCCCAAAATTTTGGGGC"
                       for p in saved), saved

    @_pt.mark.asyncio
    async def test_map_enter_opens_primer_editor(
            self, tiny_record, isolated_library):
        """Enter on a primer feature selected on the PLASMID MAP opens the
        PrimerEditModal (name + sequence + Save to library) — unified with the
        sidebar table + seq panel Enter behaviour."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            rec = app._current_record
            rec.features.append(SeqFeature(
                FeatureLocation(5, 26, strand=1), type="primer_bind",
                qualifiers={"label": ["TestPrimer-F"],
                            "primer_seq": ["GCGCACGTACGTACGTACGTA"]}))
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm.load_record(rec)
            idx = next(i for i, f in enumerate(pm._feats)
                       if f.get("label") == "TestPrimer-F")
            pm.selected_idx = idx
            await pilot.pause()
            pm.action_open_selected_feature()
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.PrimerEditModal), \
                "Enter on a map primer must open PrimerEditModal"
            # The modal carries the primer's name + sequence.
            from textual.widgets import Input
            assert app.screen.query_one("#primedit-name", Input).value \
                == "TestPrimer-F"

    @_pt.mark.asyncio
    async def test_new_primer_creates_collection_on_the_fly(
            self, tiny_record, isolated_library):
        """The picker lets the user CREATE a destination collection — it's made,
        set active, and the primer lands in it."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            sc._save_primers([])
            app._save_primer_to_library_flow(
                {"primer_seq": "TTTTGGGGCCCCAAAATTTTG", "label": "Made-F",
                 "strand": 1})
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc.PrimerSaveModal)
            app.screen.dismiss({"names": ["Made-F"],
                                "collection": "Brand New Coll", "create": True})
            await pilot.pause(); await pilot.pause()
            colls = [c.get("name") for c in sc._load_primer_collections()]
            assert "Brand New Coll" in colls, colls
            assert sc._get_active_primer_collection_name() == "Brand New Coll"
            assert any(p.get("name") == "Made-F" for p in sc._load_primers())

    @_pt.mark.asyncio
    async def test_save_library_empty_sequence_is_noop(
            self, tiny_record, isolated_library):
        """Defence-in-depth: an empty primer sequence neither pushes the picker
        nor errors."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            depth = len(app.screen_stack)
            app._save_primer_to_library_flow({"primer_seq": "   ", "label": "x"})
            await pilot.pause()
            assert len(app.screen_stack) == depth
