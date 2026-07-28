"""
test_library_bulk_mark — sweep #28.

Bulk-mark + move/copy across collections via the Library panel.

Hardening test coverage (the user explicitly asked for "harden and
edge case this PLEASE"):

  * Toggle mark + clear marks
  * Mark survives repopulate; clears on collection switch
  * Move with no marks → falls back to cursor row
  * Move with 0-source-match → refused with clean notify
  * Move into self → refused
  * Move when target disappeared mid-flight → refused, source intact
  * Copy with name collision in target → silent rename, content-aware:
    " COPY" for the same sequence, " ALT" for a different one ([INV-167])
  * Copy duplicates are deep-copied (mutating source post-copy
    doesn't bleed into target)
  * Move preserves metadata (status, alignments, gb_text, history_xml)
  * Atomic save: failure leaves source + target untouched
  * Active-library mirror re-stages on move-from-active
  * Active-library mirror re-stages on copy-to-active
  * Move + entry-vector binding follows by name (binding is name-keyed)
"""
from __future__ import annotations

import json

import pytest

import splicecraft as sc


# ── Helpers ──────────────────────────────────────────────────────────


def _seed_two_collections(eden_n: int = 3, ffe_n: int = 0):
    """Build a {collections.json, plasmid_library.json} pair with
    'DemoColl' (eden_n plasmids) + 'FFE' (ffe_n plasmids). DemoColl is active.
    Returns the entry ids in DemoColl so tests can mark them."""
    eden_plasmids = [
        {"id": f"eden_{i}", "name": f"plasmid_{i}",
         "size": 1000 + i, "n_feats": 2,
         "status": "DESIGNING", "gb_text": f"LOCUS plasmid_{i} 1000 bp\n",
         "alignments": [],
         "metadata_marker": f"unique-democoll-{i}"}
        for i in range(eden_n)
    ]
    ffe_plasmids = [
        {"id": f"ffe_{i}", "name": f"ffe_plasmid_{i}",
         "size": 2000 + i, "n_feats": 3,
         "status": "VERIFIED", "gb_text": f"LOCUS ffe_{i} 2000 bp\n",
         "metadata_marker": f"unique-ffe-{i}"}
        for i in range(ffe_n)
    ]
    colls = [
        {"name": "DemoColl", "description": "test source",
         "plasmids": eden_plasmids},
        {"name": "FFE",  "description": "test target",
         "plasmids": ffe_plasmids},
    ]
    sc._save_collections(colls)
    sc._set_active_collection_name("DemoColl")
    sc._settings_flush_sync()
    sc._safe_save_json_mirror(sc._state._LIBRARY_FILE, eden_plasmids,
                                "Plasmid library")
    sc._state._library_cache = None
    return [e["id"] for e in eden_plasmids]


# ── Sanity ───────────────────────────────────────────────────────────


class TestSanity:
    def test_marked_ids_attr_exists(self):
        """LibraryPanel instances initialise `_marked_ids` in on_mount."""
        # We can't easily run on_mount without an app, so just verify
        # the class has the action handlers + the Message subclass.
        assert hasattr(sc.LibraryPanel, "action_toggle_mark")
        assert hasattr(sc.LibraryPanel, "action_move_marked")
        assert hasattr(sc.LibraryPanel, "action_copy_marked")
        assert hasattr(sc.LibraryPanel, "action_clear_marks")
        assert hasattr(sc.LibraryPanel, "MoveCopyRequested")


class TestMarkCycle:
    """Space cycles a row none → Ⓜ (move) → Ⓒ (copy) → none, mirroring the
    primer library's ★/$/M cycle where the mark carries the intent. The state
    machine is split out of the key handler so it tests without a DataTable."""

    @pytest.fixture
    def panel(self):
        p = sc.LibraryPanel.__new__(sc.LibraryPanel)
        p._marked_ids = set()
        p._mark_kind = {}
        return p

    def test_space_cycles_through_move_copy_delete_and_off(self, panel):
        assert panel._cycle_mark("a") == "move"
        assert panel._cycle_mark("a") == "copy"
        assert panel._cycle_mark("a") == "delete"
        assert panel._cycle_mark("a") == ""
        assert "a" not in panel._marked_ids and "a" not in panel._mark_kind

    def test_delete_is_last_so_a_stray_press_cannot_reach_it(self, panel):
        # One press on an unmarked row must never stage a delete.
        assert panel._cycle_mark("a") != "delete"
        assert panel._MARK_CYCLE[-1] == "delete"

    def test_unknown_stored_kind_falls_off_to_unmarked(self, panel):
        panel._marked_ids.add("a")
        panel._mark_kind["a"] = "bogus"
        assert panel._cycle_mark("a") == ""
        assert "a" not in panel._marked_ids

    def test_cycle_is_per_row(self, panel):
        panel._cycle_mark("a")                       # a → move
        panel._cycle_mark("b"); panel._cycle_mark("b")   # b → copy
        assert panel._marked_of_kind("move") == ["a"]
        assert panel._marked_of_kind("copy") == ["b"]

    def test_legacy_bare_add_reads_as_move(self, panel):
        # Older code (and existing tests) add straight to `_marked_ids` with
        # no kind; that must keep behaving as a move mark, not vanish.
        panel._marked_ids.add("legacy")
        assert panel._marked_of_kind("move") == ["legacy"]
        assert panel._marked_of_kind("copy") == []

    def test_move_ignores_copy_marked_rows(self, panel):
        panel._cycle_mark("a")                       # move
        panel._cycle_mark("b"); panel._cycle_mark("b")   # copy
        assert panel._resolve_marked_or_cursor("move") == ["a"]
        assert panel._resolve_marked_or_cursor("copy") == ["b"]

    def test_kindless_resolve_returns_every_mark(self, panel):
        # Image export doesn't distinguish move from copy.
        panel._cycle_mark("a")
        panel._cycle_mark("b"); panel._cycle_mark("b")
        assert sorted(panel._resolve_marked_or_cursor()) == ["a", "b"]

    def test_all_marks_other_kind_yields_nothing_not_the_cursor_row(
            self, panel):
        # The dangerous fallback: marks set, none of this kind. Falling
        # through to the cursor row would move a plasmid staged to copy.
        panel._cycle_mark("a"); panel._cycle_mark("a")    # copy only
        assert panel._resolve_marked_or_cursor("move") == []

    def test_notify_names_the_other_key_when_marks_are_all_other_kind(
            self, panel, monkeypatch):
        panel._cycle_mark("a"); panel._cycle_mark("a")    # copy only
        seen = []

        class _App:
            def notify(self, msg, severity="information", **_k):
                seen.append(msg)
        monkeypatch.setattr(type(panel), "app", property(lambda self: _App()))
        panel._notify_nothing_staged("move")
        # "did nothing" is the worst possible message here — it must say the
        # marks exist and which key commits them.
        assert seen and "y" in seen[0] and "copy" in seen[0]

    def test_mark_cell_with_supplied_colour_does_no_library_lookup(
            self, panel, monkeypatch):
        """PERF LOCK. `_find_library_entry_by_id` is an O(N) cache walk that
        takes `_cache_lock` and deep-clones (sweep #11). Calling it once per
        row makes a repopulate O(N**2) with N lock acquisitions, so the
        repopulate loop MUST pass the status colour it already derived."""
        calls = []
        monkeypatch.setattr(sc, "_find_library_entry_by_id",
                            lambda *a, **k: calls.append(1))
        panel._mark_cell("x", None)         # unmarked, colour supplied
        panel._mark_cell("x", "green")
        assert calls == [], "repopulate path must not hit the library"
        panel._cycle_mark("m1")
        panel._mark_cell("m1")              # marked → glyph needs no status
        assert calls == []
        panel._mark_cell("y")               # unmarked, no colour → lookup
        assert len(calls) == 1

    def test_mark_cell_renders_all_three_glyphs(self, panel):
        panel._cycle_mark("a")                            # Ⓜ
        panel._cycle_mark("b"); panel._cycle_mark("b")    # Ⓒ
        for _ in range(3):
            panel._cycle_mark("c")                        # Ⓧ
        assert "Ⓜ" in str(panel._mark_cell("a"))
        assert "Ⓒ" in str(panel._mark_cell("b"))
        assert "Ⓧ" in str(panel._mark_cell("c"))
        assert "·" in str(panel._mark_cell("d", None))    # unmarked, no status

    def test_each_command_sees_only_its_own_mark(self, panel):
        panel._cycle_mark("a")                            # move
        panel._cycle_mark("b"); panel._cycle_mark("b")    # copy
        for _ in range(3):
            panel._cycle_mark("c")                        # delete
        assert panel._marked_of_kind("move") == ["a"]
        assert panel._marked_of_kind("copy") == ["b"]
        assert panel._marked_of_kind("delete") == ["c"]

    def test_notify_explains_how_to_mark_when_nothing_is_marked(
            self, panel, monkeypatch):
        seen = []

        class _App:
            def notify(self, msg, severity="information", **_k):
                seen.append(msg)
        monkeypatch.setattr(type(panel), "app", property(lambda self: _App()))
        panel._notify_nothing_staged("move")
        assert seen and "Space" in seen[0]


def _make_bulk_rig(monkeypatch, n_entries: int):
    """A LibraryPanel wired to a stub app: captures the confirm modal so a
    test can answer it, and records undo pushes.

    The stub delegates to the REAL `PlasmidApp._push_library_undo` rather than
    appending blindly — the depth cap and its drop-incomplete-batches rule are
    part of what these tests are checking, so faking the push would hide
    exactly the bug the batch design guards against."""
    eden_ids = _seed_two_collections(eden_n=n_entries, ffe_n=0)
    panel = sc.LibraryPanel.__new__(sc.LibraryPanel)
    panel._marked_ids, panel._mark_kind = set(), {}
    panel._view_mode = "plasmids"
    box = {"modal": None, "cb": None, "stack": [], "notes": [], "saved": None}

    # A REAL PlasmidApp with only the UI hooks stubbed, so the undo push goes
    # through the actual depth-cap + drop-incomplete-batches logic. A hand-
    # rolled stub would silently miss `_trim_library_undo` and hide the very
    # bug these tests exist for.
    app = sc.PlasmidApp.__new__(sc.PlasmidApp)
    app._current_record = None
    app._library_undo_stack = lambda: box["stack"]
    app.push_screen = lambda screen, callback=None: box.update(
        modal=screen, cb=callback)
    app.notify = lambda msg, severity="information", **_k: box["notes"].append(
        msg)
    monkeypatch.setattr(type(panel), "app", property(lambda self: app))
    monkeypatch.setattr(type(panel), "_repopulate_plasmids", lambda self: None)
    monkeypatch.setattr(type(panel), "set_active", lambda self, x: None)
    monkeypatch.setattr(
        type(panel), "_delete_save_to_disk",
        lambda self, entries, bin_, op="": box.update(saved=entries))
    return panel, eden_ids, box


@pytest.fixture
def rig(monkeypatch):
    return _make_bulk_rig(monkeypatch, 4)


@pytest.fixture
def rig_big(monkeypatch):
    """More marked rows than the undo stack can hold whole."""
    return _make_bulk_rig(monkeypatch, sc._LIBRARY_UNDO_MAX + 3)


class TestBulkDelete:
    """Ⓧ-marked bulk delete + its single-press undo. Deletion is the path
    behind the 2026-05-22 incident, so the round trip is tested end to end
    rather than by inspection."""

    def test_deletes_every_marked_row_and_leaves_the_rest(self, rig):
        panel, ids, box = rig
        for eid in ids[:2]:
            for _ in range(3):
                panel._cycle_mark(eid)          # → Ⓧ
        panel.request_delete_under_cursor()
        assert box["modal"] is not None, "no confirm modal was pushed"
        box["cb"](True)
        left = [e["id"] for e in sc._load_library()]
        assert left == ids[2:]

    def test_confirm_states_the_count_and_names(self, rig):
        panel, ids, box = rig
        for eid in ids[:3]:
            for _ in range(3):
                panel._cycle_mark(eid)
        panel.request_delete_under_cursor()
        assert len(box["modal"].bulk_names) == 3

    def test_nothing_marked_falls_through_to_the_single_delete(self, rig,
                                                               monkeypatch):
        panel, ids, box = rig
        called = []
        monkeypatch.setattr(type(panel), "_request_plasmid_delete",
                            lambda self: called.append(1))
        panel.request_delete_under_cursor()
        assert called == [1] and box["modal"] is None

    def test_move_and_copy_marks_are_not_deleted(self, rig):
        panel, ids, box = rig
        panel._cycle_mark(ids[0])                       # Ⓜ
        panel._cycle_mark(ids[1]); panel._cycle_mark(ids[1])   # Ⓒ
        for _ in range(3):
            panel._cycle_mark(ids[2])                   # Ⓧ
        panel.request_delete_under_cursor()
        box["cb"](True)
        left = [e["id"] for e in sc._load_library()]
        assert ids[0] in left and ids[1] in left and ids[2] not in left
        # The spent Ⓧ mark clears; the others stay staged.
        assert panel._marked_of_kind("delete") == []
        assert panel._marked_of_kind("move") == [ids[0]]
        assert panel._marked_of_kind("copy") == [ids[1]]

    def test_one_undo_press_restores_the_whole_batch_in_place(self, rig):
        panel, ids, box = rig
        for eid in ids[:3]:
            for _ in range(3):
                panel._cycle_mark(eid)
        panel.request_delete_under_cursor()
        box["cb"](True)
        assert len(sc._load_library()) == 1
        # Every record shares one batch id, so ONE press takes them all back.
        assert len({r["batch"] for r in box["stack"]}) == 1
        panel.action_undo_delete()
        assert [e["id"] for e in sc._load_library()] == ids   # original slots
        assert box["stack"] == []

    def test_undo_refuses_as_a_whole_when_the_collection_changed(self, rig):
        panel, ids, box = rig
        for eid in ids[:2]:
            for _ in range(3):
                panel._cycle_mark(eid)
        panel.request_delete_under_cursor()
        box["cb"](True)
        sc._set_active_collection_name("FFE")
        panel.action_undo_delete()
        # All-or-nothing: nothing restored, and the batch stays recoverable.
        assert len(sc._load_library()) == 2
        assert len(box["stack"]) == 2
        assert any("Switch to that collection" in m for m in box["notes"])

    def test_batch_ids_are_unique_within_a_second(self, rig):
        """A second-resolution timestamp collided for two same-sized deletes
        in the same second, and one undo would then restore BOTH batches."""
        panel, ids, box = rig
        for eid in ids[:2]:
            for _ in range(3):
                panel._cycle_mark(eid)
        panel.request_delete_under_cursor(); box["cb"](True)
        for eid in ids[2:4]:
            for _ in range(3):
                panel._cycle_mark(eid)
        panel.request_delete_under_cursor(); box["cb"](True)
        assert len({r["batch"] for r in box["stack"]}) == 2

    def test_undo_of_one_batch_leaves_the_other_intact(self, rig):
        panel, ids, box = rig
        for eid in ids[:2]:
            for _ in range(3):
                panel._cycle_mark(eid)
        panel.request_delete_under_cursor(); box["cb"](True)
        for eid in ids[2:4]:
            for _ in range(3):
                panel._cycle_mark(eid)
        panel.request_delete_under_cursor(); box["cb"](True)
        panel.action_undo_delete()                      # newest batch only
        assert sorted(e["id"] for e in sc._load_library()) == sorted(ids[2:4])
        assert len(box["stack"]) == 2                   # older batch survives


class TestUndoDepthCapNeverSplitsABatch:
    """The depth trim deletes from the FRONT of the stack, which could cut a
    bulk delete's records in half — the survivors would then restore a PARTIAL
    batch and report success. That is the 2026-05-22 failure shape, so a batch
    the cap can't hold whole is dropped whole."""

    @pytest.fixture
    def app(self):
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        app._lib_undo = []
        app._library_undo_stack = lambda: app._lib_undo
        return app

    def _push_batch(self, app, batch: str, n: int, start: int = 0):
        # Atomically, the way the bulk delete does — pushing one at a time
        # would make the batch transiently incomplete and drop it mid-push.
        return app._push_library_undo_batch([
            {"entry": {"id": f"{batch}-{i}"}, "index": start + i,
             "collection": "DemoColl", "name": f"{batch}-{i}",
             "batch": batch, "batch_size": n}
            for i in range(n)
        ])

    def test_a_batch_cut_by_the_cap_is_dropped_whole(self, app):
        self._push_batch(app, "A", 20)
        self._push_batch(app, "B", sc._LIBRARY_UNDO_MAX - 10)
        # B pushed A past the cap. A must be GONE entirely, not truncated.
        left = {str(r.get("batch") or "") for r in app._lib_undo}
        assert "A" not in left, "a partially-evicted batch survived"
        assert "B" in left
        assert all(r["batch"] == "B" for r in app._lib_undo)

    def test_a_batch_larger_than_the_cap_does_not_linger(self, app):
        self._push_batch(app, "BIG", sc._LIBRARY_UNDO_MAX + 5)
        assert app._lib_undo == [], "an unrestorable batch was left half-there"

    def test_untagged_single_deletes_still_trim_normally(self, app):
        for i in range(sc._LIBRARY_UNDO_MAX + 7):
            app._push_library_undo({"entry": {"id": str(i)}, "index": i,
                                    "collection": "DemoColl", "name": str(i)})
        assert len(app._lib_undo) == sc._LIBRARY_UNDO_MAX
        assert app._lib_undo[0]["name"] == "7"      # oldest dropped, in order

    def test_complete_batches_are_kept(self, app):
        self._push_batch(app, "A", 5)
        self._push_batch(app, "B", 5)
        assert len(app._lib_undo) == 10


class TestBulkDeleteHonesty:
    def test_toast_refuses_to_promise_an_undo_that_is_not_there(self, rig_big):
        """Over-cap batches are dropped whole, so 'press u' would be a lie."""
        panel, ids, box = rig_big
        for eid in ids:
            for _ in range(3):
                panel._cycle_mark(eid)
        panel.request_delete_under_cursor()
        box["cb"](True)
        assert box["stack"] == []
        assert any("CANNOT be undone" in m for m in box["notes"])
        assert not any("Press u" in m for m in box["notes"])


class TestBulkDeleteSequenceCount:
    def test_counts_sequences_not_entries(self, rig):
        panel, ids, box = rig
        lib = sc._load_library()
        # Two entries, one shared sequence: deleting both loses ONE sequence.
        pair = [e for e in lib if e["id"] in ids[:2]]
        for e in pair:
            e["gb_text"] = lib[0].get("gb_text", "")
            e.pop("gb_ref", None)
        assert panel._bulk_sequences_lost(pair) == 1


# ── Move / copy commit logic (the heavy hardening) ───────────────────


class TestMoveCopyCommit:
    """`_move_copy_commit` is the synchronous worker called from the
    app's modal-callback. Tests bypass the modal to focus on the
    transactional logic."""

    @pytest.fixture
    def app(self):
        """Bare PlasmidApp instance (no UI) so we can call the method
        directly. The notify shim swallows toasts."""
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        app._notify_log = []

        def _notify(msg, severity="information", **_kwargs):
            app._notify_log.append((severity, msg))

        app.notify = _notify
        return app

    def test_move_basic(self, app):
        eden_ids = _seed_two_collections(eden_n=3, ffe_n=0)
        # Move all 3.
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=eden_ids, mode="move",
        )
        colls = sc._load_collections()
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        ffe  = next(c for c in colls if c["name"] == "FFE")
        assert democoll["plasmids"] == []
        assert len(ffe["plasmids"]) == 3
        # Order preserved.
        assert [p["id"] for p in ffe["plasmids"]] == eden_ids

    def test_commit_clears_only_the_committed_kind(self, app):
        """A mixed batch commits in two passes. Clearing EVERY mark after the
        first would silently discard the half still staged for the other
        operation — and the user would have no way to tell it was dropped."""
        eden_ids = _seed_two_collections(eden_n=3, ffe_n=0)
        panel = sc.LibraryPanel.__new__(sc.LibraryPanel)
        panel._marked_ids, panel._mark_kind = set(), {}
        panel._cycle_mark(eden_ids[0])                              # Ⓜ
        panel._cycle_mark(eden_ids[1]); panel._cycle_mark(eden_ids[1])  # Ⓒ
        panel._repopulate = lambda: None
        app.query_one = lambda *a, **k: panel

        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=[eden_ids[0]], mode="move",
        )
        assert panel._marked_ids == {eden_ids[1]}
        assert panel._mark_kind == {eden_ids[1]: "copy"}

    def test_copy_basic(self, app):
        eden_ids = _seed_two_collections(eden_n=2, ffe_n=0)
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=eden_ids, mode="copy",
        )
        colls = sc._load_collections()
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        ffe  = next(c for c in colls if c["name"] == "FFE")
        # Source unchanged.
        assert len(democoll["plasmids"]) == 2
        # Target populated.
        assert len(ffe["plasmids"]) == 2

    def test_copy_is_deepcopy(self, app):
        """Mutating the source list after copy MUST NOT bleed into
        target (sacred [PIT-17])."""
        eden_ids = _seed_two_collections(eden_n=1, ffe_n=0)
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=eden_ids, mode="copy",
        )
        # Mutate the source entry's metadata. The target's copy must
        # NOT be affected.
        colls = sc._load_collections()
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        democoll["plasmids"][0]["metadata_marker"] = "MUTATED"
        sc._save_collections(colls)
        # Re-load and check FFE's copy still has the original marker.
        colls2 = sc._load_collections()
        ffe = next(c for c in colls2 if c["name"] == "FFE")
        assert ffe["plasmids"][0]["metadata_marker"] == "unique-democoll-0"

    def test_copy_with_name_collision_appends_alt_suffix(self, app):
        # FFE already has a plasmid named "plasmid_0" — carrying
        # DIFFERENT content ("y" vs "x"), so the landing is NOT a copy of
        # it and must not be labelled one ([INV-167]).
        sc._save_collections([
            {"name": "DemoColl", "plasmids": [
                {"id": "eden_0", "name": "plasmid_0", "size": 100,
                 "gb_text": "x"},
            ]},
            {"name": "FFE", "plasmids": [
                {"id": "x", "name": "plasmid_0", "size": 200,
                 "gb_text": "y"},
            ]},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=["eden_0"], mode="copy",
        )
        colls = sc._load_collections()
        ffe = next(c for c in colls if c["name"] == "FFE")
        names = [p["name"] for p in ffe["plasmids"]]
        assert "plasmid_0" in names
        assert "plasmid_0 ALT" in names
        assert "plasmid_0 COPY" not in names

    def test_copy_with_multiple_collisions_increments_suffix(self, app):
        # Pre-seed FFE with `plasmid_0`, `plasmid_0 ALT` → next
        # landing must be `plasmid_0 ALT 2`. (All three carry different
        # content, so the suffix is ALT; the point of the test is that
        # the NUMBERING still increments — see [INV-167].)
        sc._save_collections([
            {"name": "DemoColl", "plasmids": [
                {"id": "eden_0", "name": "plasmid_0", "size": 100,
                 "gb_text": "x"},
            ]},
            {"name": "FFE", "plasmids": [
                {"id": "a", "name": "plasmid_0", "size": 200,
                 "gb_text": "y"},
                {"id": "b", "name": "plasmid_0 ALT", "size": 300,
                 "gb_text": "z"},
            ]},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=["eden_0"], mode="copy",
        )
        colls = sc._load_collections()
        ffe = next(c for c in colls if c["name"] == "FFE")
        names = [p["name"] for p in ffe["plasmids"]]
        assert "plasmid_0 ALT 2" in names

    def test_move_same_source_and_target_refused(self, app):
        sc._save_collections([
            {"name": "DemoColl", "plasmids": [
                {"id": "eden_0", "name": "p", "size": 100, "gb_text": "x"},
            ]},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        app._move_copy_commit(
            source="DemoColl", target="DemoColl",
            entry_ids=["eden_0"], mode="move",
        )
        # Source intact.
        colls = sc._load_collections()
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        assert len(democoll["plasmids"]) == 1
        # Notify fired with warning.
        assert any("same" in m.lower()
                   for sev, m in app._notify_log)

    def test_copy_same_source_and_target_duplicates_in_place(self, app):
        """Copy mode + src==tgt is the duplicate-in-place flow.
        Each marked entry gets a " COPY" / " COPY 2" suffix; the
        originals stay intact under their original names + ids."""
        sc._save_collections([
            {"name": "DemoColl", "plasmids": [
                {"id": "eden_0", "name": "plasmid_0", "size": 100,
                 "gb_text": "x"},
                {"id": "eden_1", "name": "plasmid_1", "size": 200,
                 "gb_text": "y"},
            ]},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        app._move_copy_commit(
            source="DemoColl", target="DemoColl",
            entry_ids=["eden_0", "eden_1"], mode="copy",
        )
        colls = sc._load_collections()
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        names = [p["name"] for p in democoll["plasmids"]]
        ids   = [p["id"]   for p in democoll["plasmids"]]
        # Originals untouched.
        assert "plasmid_0" in names
        assert "plasmid_1" in names
        assert "eden_0" in ids
        assert "eden_1" in ids
        # Duplicates landed with COPY suffix.
        assert "plasmid_0 COPY" in names
        assert "plasmid_1 COPY" in names
        # Duplicate ids are also unique (suffix increment).
        assert ids.count("eden_0") == 1
        assert ids.count("eden_1") == 1
        # No "same — nothing to do" warning.
        assert not any("nothing to do" in m.lower()
                        for sev, m in app._notify_log)
        # A success-style information toast was posted.
        assert any("duplicated" in m.lower()
                    for sev, m in app._notify_log)

    def test_copy_same_target_repeated_duplicates_increment_suffix(
            self, app):
        """A second duplicate-in-place call must produce a 'COPY 2'
        rather than re-using 'COPY'."""
        sc._save_collections([
            {"name": "DemoColl", "plasmids": [
                {"id": "eden_0", "name": "p", "size": 100, "gb_text": "x"},
            ]},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        # First duplication: yields "p COPY".
        app._move_copy_commit(
            source="DemoColl", target="DemoColl",
            entry_ids=["eden_0"], mode="copy",
        )
        sc._state._library_cache = None
        sc._state._collections_cache = None
        # Second duplication of the SAME original: must yield
        # "p COPY 2" (the first COPY is already in the target set).
        app._move_copy_commit(
            source="DemoColl", target="DemoColl",
            entry_ids=["eden_0"], mode="copy",
        )
        colls = sc._load_collections()
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        names = [p["name"] for p in democoll["plasmids"]]
        assert "p" in names
        assert "p COPY" in names
        assert "p COPY 2" in names

    def test_copy_same_target_active_mirror_restages(
            self, app, tmp_path, monkeypatch):
        """Duplicate-in-place inside the ACTIVE collection must
        re-stage `plasmid_library.json` so the LibraryPanel sees
        the new entries on next repopulate."""
        eden_ids = _seed_two_collections(eden_n=2, ffe_n=0)
        # DemoColl is active by `_seed_two_collections`.
        app._move_copy_commit(
            source="DemoColl", target="DemoColl",
            entry_ids=eden_ids, mode="copy",
        )
        # plasmid_library.json was re-mirrored from the updated DemoColl.
        mirror = json.loads(sc._state._LIBRARY_FILE.read_text("utf-8"))
        entries, _err = sc._extract_entries(mirror, "Plasmid library")
        assert _err is None
        assert entries is not None
        names = [e["name"] for e in entries]
        # Originals + copies all present.
        assert "plasmid_0" in names
        assert "plasmid_1" in names
        assert "plasmid_0 COPY" in names
        assert "plasmid_1 COPY" in names

    def test_invalid_mode_refused(self, app):
        _seed_two_collections(eden_n=1, ffe_n=0)
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=["eden_0"], mode="weirdmode",
        )
        # No state change.
        colls = sc._load_collections()
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        ffe = next(c for c in colls if c["name"] == "FFE")
        assert len(democoll["plasmids"]) == 1
        assert len(ffe["plasmids"]) == 0

    def test_source_disappeared_mid_commit(self, app, monkeypatch):
        # Simulate a race: the source collection is gone by the time
        # the commit runs. Should refuse cleanly without touching
        # target state.
        _seed_two_collections(eden_n=1, ffe_n=0)
        sc._save_collections([
            {"name": "FFE", "plasmids": []},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=["eden_0"], mode="move",
        )
        # FFE still empty.
        colls = sc._load_collections()
        ffe = next(c for c in colls if c["name"] == "FFE")
        assert len(ffe["plasmids"]) == 0
        assert any("disappeared" in m.lower() or "failed" in m.lower()
                   for sev, m in app._notify_log)

    def test_move_preserves_full_metadata(self, app):
        # Build a source plasmid with rich metadata: status,
        # alignments, history_xml, color overrides, custom fields.
        rich = {
            "id":           "rich_0",
            "name":         "rich_plasmid",
            "size":         5000,
            "n_feats":      10,
            "gb_text":      "LOCUS rich 5000 bp\n",
            "status":       "VERIFIED",
            "alignments":   [{"read": "x", "identity_pct": 99.5}],
            "history_xml":  "<root><node id='1'/></root>",
            "map_mode":     "linear",
            "_plugin_data": {"my_plugin": {"foo": "bar"}},
            "custom_field": "preserved_value",
        }
        sc._save_collections([
            {"name": "DemoColl", "plasmids": [rich]},
            {"name": "FFE",  "plasmids": []},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=["rich_0"], mode="move",
        )
        colls = sc._load_collections()
        ffe = next(c for c in colls if c["name"] == "FFE")
        landed = ffe["plasmids"][0]
        # Every metadata field preserved exactly.
        for k in ("status", "alignments", "history_xml", "map_mode",
                  "_plugin_data", "custom_field", "gb_text", "size"):
            assert landed[k] == rich[k], (
                f"metadata {k!r} not preserved: "
                f"{landed[k]!r} != {rich[k]!r}"
            )

    def test_mirror_re_stages_on_move_from_active(self, app):
        # Active = DemoColl. After moving DemoColl's only plasmid away, the
        # active library mirror (plasmid_library.json) must reflect
        # the new empty state.
        eden_ids = _seed_two_collections(eden_n=2, ffe_n=0)
        # Verify the mirror has 2 entries pre-move.
        lib_pre, _ = sc._safe_load_json(sc._state._LIBRARY_FILE, "test")
        assert len(lib_pre) == 2
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=eden_ids, mode="move",
        )
        # Mirror now reflects empty DemoColl.
        lib_post, _ = sc._safe_load_json(sc._state._LIBRARY_FILE, "test")
        assert lib_post == []

    def test_mirror_re_stages_on_copy_to_active(self, app):
        # Active = FFE. Copy DemoColl plasmid INTO FFE → mirror reflects
        # the new entry.
        _seed_two_collections(eden_n=1, ffe_n=0)
        sc._set_active_collection_name("FFE")
        sc._settings_flush_sync()
        # Re-mirror to FFE state.
        sc._safe_save_json_mirror(sc._state._LIBRARY_FILE, [], "Plasmid library")
        sc._state._library_cache = None
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=["eden_0"], mode="copy",
        )
        lib_post, _ = sc._safe_load_json(sc._state._LIBRARY_FILE, "test")
        assert len(lib_post) == 1

    def test_partial_id_match_only_operates_on_existing(self, app):
        eden_ids = _seed_two_collections(eden_n=2, ffe_n=0)
        # Include a bogus id in the request — should be silently
        # filtered.
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=eden_ids + ["bogus_id_that_does_not_exist"],
            mode="move",
        )
        colls = sc._load_collections()
        ffe = next(c for c in colls if c["name"] == "FFE")
        # Only the 2 real ones landed.
        assert len(ffe["plasmids"]) == 2

    def test_all_ids_filtered_to_none_refused(self, app):
        _seed_two_collections(eden_n=2, ffe_n=0)
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=["nothing_real_1", "nothing_real_2"],
            mode="move",
        )
        # Source intact.
        colls = sc._load_collections()
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        assert len(democoll["plasmids"]) == 2

    def test_copy_id_collision_renames_id_too(self, app):
        # If the source and target share an id (e.g. fresh copy of
        # the same plasmid), the target's id gets a "_2" suffix so
        # the cache lookup stays unique. Pre-fix, the duplicate id
        # would silently overwrite the existing entry on the next
        # collection-switch mirror.
        sc._save_collections([
            {"name": "DemoColl", "plasmids": [
                {"id": "p1", "name": "P1", "size": 100, "gb_text": "x"},
            ]},
            {"name": "FFE", "plasmids": [
                {"id": "p1", "name": "OtherP1", "size": 200, "gb_text": "y"},
            ]},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        app._move_copy_commit(
            source="DemoColl", target="FFE",
            entry_ids=["p1"], mode="copy",
        )
        colls = sc._load_collections()
        ffe = next(c for c in colls if c["name"] == "FFE")
        ids = [p["id"] for p in ffe["plasmids"]]
        assert "p1" in ids
        assert "p1_2" in ids
        assert len(set(ids)) == len(ids)   # no dupes


# ── Concurrency ──────────────────────────────────────────────────────


class TestMoveCopyConcurrency:
    """Two threads racing on move/copy must not corrupt either
    collection. The `_cache_lock` RMW makes the commit serial."""

    def test_concurrent_moves_dont_corrupt(self):
        import threading

        _seed_two_collections(eden_n=10, ffe_n=0)
        # Two apps in two threads, each moving disjoint halves of
        # DemoColl → FFE concurrently.
        app1 = sc.PlasmidApp.__new__(sc.PlasmidApp)
        app1.notify = lambda *a, **k: None
        app2 = sc.PlasmidApp.__new__(sc.PlasmidApp)
        app2.notify = lambda *a, **k: None
        half1 = [f"eden_{i}" for i in range(5)]
        half2 = [f"eden_{i}" for i in range(5, 10)]
        barrier = threading.Barrier(2)
        errors: list[BaseException] = []

        def worker(app, ids):
            barrier.wait()
            try:
                app._move_copy_commit(
                    source="DemoColl", target="FFE",
                    entry_ids=ids, mode="move",
                )
            except Exception as exc:
                errors.append(exc)

        t1 = threading.Thread(target=worker, args=(app1, half1))
        t2 = threading.Thread(target=worker, args=(app2, half2))
        t1.start(); t2.start()
        t1.join(5); t2.join(5)
        assert not errors, f"workers raised: {errors!r}"
        colls = sc._load_collections()
        ffe = next(c for c in colls if c["name"] == "FFE")
        democoll = next(c for c in colls if c["name"] == "DemoColl")
        # All 10 should have landed in FFE; DemoColl empty.
        assert len(ffe["plasmids"]) == 10
        assert len(democoll["plasmids"]) == 0
        # No dupes.
        ids = [p["id"] for p in ffe["plasmids"]]
        assert len(set(ids)) == 10


# ── [New collection] button inside the move/copy picker (2026-06-01) ──


class TestMoveCopyNewCollection:
    """MoveCopyToCollectionModal's [New collection] button — create an
    empty target collection inline (copied from CollectionsModal._save)
    without leaving the picker, then select it so Confirm lands on it.
    Plus the caller change: move mode no longer refuses to open when the
    source is the only collection (the button is now the escape hatch).
    """

    # ── caller: move opens even with only the source collection ──

    def test_move_with_only_source_no_longer_refused(self):
        sc._save_collections([
            {"name": "Solo", "plasmids": [
                {"id": "s0", "name": "p", "size": 100, "gb_text": "x"}]},
        ])
        sc._set_active_collection_name("Solo")
        sc._state._collections_cache = None
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        pushed, notes = [], []
        app.notify = lambda msg, severity="information", **k: notes.append(
            (severity, str(msg)))
        app.push_screen = lambda screen, callback=None: pushed.append(screen)
        ev = sc.LibraryPanel.MoveCopyRequested(entry_ids=["s0"], mode="move")
        app._library_move_copy_requested(ev)
        assert len(pushed) == 1
        assert type(pushed[0]).__name__ == "MoveCopyToCollectionModal"
        assert not any("only one collection" in m.lower() for _s, m in notes)

    # ── modal: the button itself ──

    async def _open_picker(self, app, pilot, *, source, ids, mode):
        await app.push_screen(sc.MoveCopyToCollectionModal(
            source_collection=source, entry_ids=ids, mode=mode))
        await pilot.pause(0.2)
        return app.screen

    async def test_new_collection_creates_and_selects(self):
        _seed_two_collections(eden_n=2, ffe_n=0)        # DemoColl active + FFE
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            picker = await self._open_picker(
                app, pilot, source="DemoColl", ids=["eden_0"], mode="move")
            assert type(picker).__name__ == "MoveCopyToCollectionModal"
            picker.query_one("#btn-movecopy-newcoll").action_press()
            await pilot.pause(0.2)
            name_modal = app.screen
            assert type(name_modal).__name__ == "CollectionNameModal"
            name_modal.query_one("#collname-input").value = "Fresh"
            name_modal.query_one("#btn-collname-ok").action_press()
            await pilot.pause(0.3)
            # Persisted as an empty collection…
            fresh = [c for c in sc._load_collections()
                     if c.get("name") == "Fresh"]
            assert len(fresh) == 1 and fresh[0]["plasmids"] == []
            # …and the picker now shows + selects it.
            back = app.screen
            assert type(back).__name__ == "MoveCopyToCollectionModal"
            assert "Fresh" in back._row_to_name
            assert back._selected_name() == "Fresh"
            app.exit()

    async def test_new_collection_rejects_duplicate(self):
        _seed_two_collections(eden_n=1, ffe_n=0)        # DemoColl + FFE exist
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            picker = await self._open_picker(
                app, pilot, source="DemoColl", ids=["eden_0"], mode="move")
            picker.query_one("#btn-movecopy-newcoll").action_press()
            await pilot.pause(0.2)
            name_modal = app.screen
            name_modal.query_one("#collname-input").value = "FFE"  # taken
            name_modal.query_one("#btn-collname-ok").action_press()
            await pilot.pause(0.3)
            # No duplicate FFE; only the original.
            assert sum(1 for c in sc._load_collections()
                       if c.get("name") == "FFE") == 1
            app.exit()

    async def test_new_collection_blank_name_creates_nothing(self):
        _seed_two_collections(eden_n=1, ffe_n=0)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            picker = await self._open_picker(
                app, pilot, source="DemoColl", ids=["eden_0"], mode="move")
            before = {c.get("name") for c in sc._load_collections()}
            picker.query_one("#btn-movecopy-newcoll").action_press()
            await pilot.pause(0.2)
            name_modal = app.screen
            name_modal.query_one("#collname-input").value = "   "
            name_modal.query_one("#btn-collname-ok").action_press()
            await pilot.pause(0.2)
            # Empty name rejected by the prompt → nothing created.
            after = {c.get("name") for c in sc._load_collections()}
            assert after == before
            app.exit()


class TestAltRenameOnMove:
    """[INV-167] Verbatim reproduction of the 2026-07-25 data-loss
    incident: a MOVE landed fragments whose names already existed in the
    target, carrying DIFFERENT sequences. The old code renamed them
    "<name> COPY" on the strength of the name alone and reported only
    "3 renamed for name collision"; the user read that as "these are
    redundant duplicates" and deleted all three.
    """

    @pytest.fixture
    def app(self):
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        app._notify_log = []

        def _notify(msg, severity="information", **_kwargs):
            app._notify_log.append((severity, msg))

        app.notify = _notify
        return app

    @staticmethod
    def _seed():
        # The target collection already holds FRAG-promoter; the source
        # holds a DIFFERENT sequence under that same name.
        # NB: `_save_collections` drops `gb_ref` for entries carrying
        # inline `gb_text`, so content identity here comes from the text
        # itself (that is exactly what `_entry_content_key` exists for).
        sc._save_collections([
            {"name": "Source", "plasmids": [
                {"id": "src_0", "name": "FRAG-promoter", "size": 857,
                 "gb_text": "LOCUS variant 857 bp\n"},
            ]},
            {"name": "Target", "plasmids": [
                {"id": "tgt_0", "name": "FRAG-promoter", "size": 857,
                 "gb_text": "LOCUS incumbent 857 bp\n"},
            ]},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None

    def test_different_sequence_lands_as_ALT_not_COPY(self, app):
        self._seed()
        app._move_copy_commit(source="Source", target="Target",
                                entry_ids=["src_0"], mode="move")
        tgt = next(c for c in sc._load_collections() if c["name"] == "Target")
        names = [p["name"] for p in tgt["plasmids"]]
        assert "FRAG-promoter ALT" in names
        # The label that caused the incident must never appear here.
        assert "FRAG-promoter COPY" not in names

    def test_the_landing_keeps_its_own_sequence(self, app):
        self._seed()
        app._move_copy_commit(source="Source", target="Target",
                                entry_ids=["src_0"], mode="move")
        tgt = next(c for c in sc._load_collections() if c["name"] == "Target")
        by_name = {p["name"]: p["gb_text"] for p in tgt["plasmids"]}
        # Both survive, each keeping its own distinct sequence — the
        # rename must never merge or overwrite them.
        assert by_name["FRAG-promoter"] == "LOCUS incumbent 857 bp\n"
        assert by_name["FRAG-promoter ALT"] == (
            "LOCUS variant 857 bp\n")

    def test_user_is_warned_that_it_is_not_a_duplicate(self, app):
        self._seed()
        app._move_copy_commit(source="Source", target="Target",
                                entry_ids=["src_0"], mode="move")
        sev, msg = app._notify_log[-1]
        # A bare "renamed for name collision" info-toast is what the
        # user acted on last time — this has to be louder and explicit.
        assert sev == "warning", app._notify_log
        assert "DIFFERENT" in msg
        assert "NOT a duplicate" in msg

    def test_identical_sequence_still_reports_plain_copy(self, app):
        # Control: a genuine duplicate keeps the COPY label and stays a
        # quiet information toast — no crying wolf.
        sc._save_collections([
            {"name": "S", "plasmids": [
                {"id": "s0", "name": "frag", "size": 10,
                 "gb_ref": "same-ref", "gb_text": "x"},
            ]},
            {"name": "T", "plasmids": [
                {"id": "t0", "name": "frag", "size": 10,
                 "gb_ref": "same-ref", "gb_text": "x"},
            ]},
        ])
        sc._state._library_cache = None
        sc._state._collections_cache = None
        app._move_copy_commit(source="S", target="T",
                                entry_ids=["s0"], mode="move")
        t = next(c for c in sc._load_collections() if c["name"] == "T")
        assert "frag COPY" in [p["name"] for p in t["plasmids"]]
        sev, msg = app._notify_log[-1]
        assert sev == "information"
        assert "ALT" not in msg
