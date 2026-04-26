"""
test_feature_library_screen — FeatureLibraryScreen + ColorPickerModal.

Covers the full-screen library workbench introduced in the v0.3.2 rework.
Scope:
  * Clicking the `Features` menu bar item pushes FeatureLibraryScreen
    (direct open, no dropdown).
  * Screen mounts cleanly with zero / one / many entries.
  * Add / Edit / Duplicate / Remove / Color / Cycle-Strand actions
    mutate the in-memory list and mark dirty; persistence happens only
    on Save (action_save) or via the unsaved-quit prompt.
  * Closing with pending edits pushes UnsavedQuitModal.
  * ColorPickerModal returns the expected dict shape for each button path.

All tests use the `_protect_user_data` autouse fixture, so the real
features.json is untouched. The fixture also redirects
_FEATURE_COLORS_FILE.
"""
from __future__ import annotations

import pytest

import splicecraft as sc
from textual.widgets import (
    Button, DataTable, Input, Label, RadioButton, Static, TextArea,
)


_BASELINE = (160, 48)


# ═══════════════════════════════════════════════════════════════════════════════
# Menu routing
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeaturesMenuOpensLibrary:
    """Clicking the `Features` entry in the top menu bar must push
    FeatureLibraryScreen directly — no dropdown. Regression guard for the
    v0.3.2 rework where `Features` changed from a dropdown to a workbench."""

    async def test_click_features_menu_opens_library_screen(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            mb = app.query_one(sc.MenuBar)
            item = mb.query_one("#menu-features", Static)
            r = item.region
            await pilot.click(offset=(r.x + 1, r.y))
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.FeatureLibraryScreen)


# ═══════════════════════════════════════════════════════════════════════════════
# Screen mount + population
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatureLibraryScreenMount:
    """The screen must mount without error for realistic library states."""

    async def test_mounts_with_empty_library(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            tbl = app.screen.query_one("#flib-table", DataTable)
            assert tbl.row_count == 0

    async def test_mounts_with_seeded_entries(self, tiny_record):
        sc._save_features([
            {"name": "lacZ", "feature_type": "CDS",
             "sequence": "ATG", "strand": 1,
             "qualifiers": {}, "description": ""},
            {"name": "tac",  "feature_type": "promoter",
             "sequence": "TTG", "strand": 1, "color": "#00FF00",
             "qualifiers": {}, "description": ""},
            {"name": "ori",  "feature_type": "rep_origin",
             "sequence": "GCA", "strand": 0,
             "qualifiers": {}, "description": ""},
        ])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            tbl = app.screen.query_one("#flib-table", DataTable)
            assert tbl.row_count == 3

    async def test_preview_renders_first_entry(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATGACC", "strand": 1,
            "qualifiers": {"gene": ["lacZ"]}, "description": "",
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            snip = app.screen.query_one(sc._FeatureSnippetPanel)
            assert snip._entry is not None
            assert snip._entry["name"] == "lacZ"


# ═══════════════════════════════════════════════════════════════════════════════
# CRUD actions
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatureLibraryCrud:
    """Remove / Duplicate / Cycle-Strand mutate the in-memory list and
    set the dirty flags. Persistence only happens via action_save.
    """

    async def test_remove_buffers_until_save(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_remove()
            await pilot.pause()
            # In-memory list shows the removal; disk does not.
            assert app.screen._entries == []
            assert app.screen._has_pending_changes is True
            sc._features_cache = None
            assert len(sc._load_features()) == 1
            # action_save persists.
            app.screen.action_save()
            await pilot.pause()
            sc._features_cache = None
            assert sc._load_features() == []
            assert app.screen._has_pending_changes is False

    async def test_duplicate_adds_copy_suffix(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_duplicate()
            await pilot.pause()
            # Buffered: in-memory has the dup, disk still has just the
            # original until action_save runs.
            assert len(app.screen._entries) == 2
            assert "(copy)" in app.screen._entries[1]["name"]
            sc._features_cache = None
            assert len(sc._load_features()) == 1
            app.screen.action_save()
            await pilot.pause()
            sc._features_cache = None
            loaded = sc._load_features()
            assert len(loaded) == 2
            assert "(copy)" in loaded[1]["name"]

    async def test_cycle_strand_buffers_until_save(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            # Four-step cycle: 1 → -1 → 0 → 2 → 1.
            for expected in (-1, 0, 2, 1):
                app.screen.action_strand()
                await pilot.pause()
                # Disk still untouched until save.
                sc._features_cache = None
                assert sc._load_features()[0]["strand"] == 1
                # In-memory entry reflects the cycle step.
                assert app.screen._entries[0]["strand"] == expected
                # Edited entry is dirty.
                assert 0 in app.screen._dirty_indices
            app.screen.action_save()
            await pilot.pause()
            sc._features_cache = None
            assert sc._load_features()[0]["strand"] == 1
            assert app.screen._dirty_indices == set()

    async def test_close_pops_back_to_main(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.FeatureLibraryScreen)
            app.screen.action_close()
            await pilot.pause()
            await pilot.pause(0.05)
            assert not isinstance(app.screen, sc.FeatureLibraryScreen)


class TestFeatureLibraryUnsavedFlow:
    """The deferred-save model: dirty flags, asterisk prefix in the
    table, and the UnsavedQuitModal-on-close gate.
    """

    async def test_dirty_entry_gets_asterisk_prefix(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_strand()
            await pilot.pause()
            tbl = app.screen.query_one("#flib-table", DataTable)
            row_keys = list(tbl.rows.keys())
            cell = tbl.get_cell(row_keys[0], list(tbl.columns.keys())[0])
            assert str(cell).startswith("*"), (
                f"Dirty entry name should be prefixed with '*', got {cell!r}"
            )

    async def test_title_marks_pending_changes(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            title = app.screen.query_one("#flib-title", Static)
            assert "*" not in str(title.render()), (
                f"Clean title should not contain '*': {title.render()!r}"
            )
            app.screen.action_strand()
            await pilot.pause()
            assert "*" in str(title.render()), (
                f"Dirty title should contain '*': {title.render()!r}"
            )

    async def test_close_with_dirty_pushes_unsaved_quit_modal(
        self, tiny_record,
    ):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_strand()
            await pilot.pause()
            app.screen.action_close()
            await pilot.pause()
            await pilot.pause(0.1)
            # Modal pushed on top — library screen is still in the stack.
            assert isinstance(app.screen, sc.UnsavedQuitModal)

    async def test_unsaved_quit_save_persists_then_pops(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_strand()  # cycle 1 → -1
            await pilot.pause()
            app.screen.action_close()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.UnsavedQuitModal)
            # Click "Save & Quit".
            app.screen.query_one("#btn-save-quit", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Both modals popped → back to main app.
            assert not isinstance(app.screen, sc.UnsavedQuitModal)
            assert not isinstance(app.screen, sc.FeatureLibraryScreen)
            sc._features_cache = None
            assert sc._load_features()[0]["strand"] == -1

    async def test_unsaved_quit_abandon_pops_without_persisting(
        self, tiny_record,
    ):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_strand()  # cycle 1 → -1 in memory
            await pilot.pause()
            app.screen.action_close()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.UnsavedQuitModal)
            app.screen.query_one("#btn-abandon", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert not isinstance(app.screen, sc.FeatureLibraryScreen)
            sc._features_cache = None
            # Disk still shows original strand=1.
            assert sc._load_features()[0]["strand"] == 1

    async def test_unsaved_quit_cancel_keeps_screen_open(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_strand()
            await pilot.pause()
            app.screen.action_close()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.UnsavedQuitModal)
            app.screen.query_one("#btn-cancel-quit", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Modal popped, library screen still open with dirty state.
            assert isinstance(app.screen, sc.FeatureLibraryScreen)
            assert app.screen._has_pending_changes is True

    async def test_remove_shifts_dirty_indices(self, tiny_record):
        sc._save_features([
            {"name": "a", "feature_type": "CDS", "sequence": "ATG", "strand": 1},
            {"name": "b", "feature_type": "CDS", "sequence": "ATG", "strand": 1},
            {"name": "c", "feature_type": "CDS", "sequence": "ATG", "strand": 1},
        ])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            # Mark indices 0 and 2 dirty.
            screen._mark_dirty(0)
            screen._mark_dirty(2)
            assert screen._dirty_indices == {0, 2}
            # Remove middle entry (index 1).
            screen._selected_index = 1
            screen.action_remove()
            await pilot.pause()
            # Index 0 still dirty, former index 2 is now index 1.
            assert screen._dirty_indices == {0, 1}, (
                f"After removing idx=1 from {{0,2}}, expected {{0,1}}; "
                f"got {screen._dirty_indices}"
            )

    async def test_save_without_changes_is_noop(self, tiny_record):
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_save()  # nothing to save
            await pilot.pause()
            assert app.screen._has_pending_changes is False

    async def test_abandon_does_not_poison_cache(self, tiny_record):
        """Regression guard: in-place mutations of feature dicts during
        a deferred-save session must not leak into ``_features_cache``.

        Pre-fix, ``_load_features()`` returned a shallow ``list(...)``
        of the cache, so dict refs were shared. FeatureLibraryScreen
        mutated entries in place (rename/strand/color), and even after
        the user picked Abandon those mutations stuck in the cache —
        next consumer of ``_load_features()`` (a fresh
        FeatureLibraryScreen, the DomesticatorModal feature picker,
        etc.) would see the abandoned edit as if it had been saved.
        Fix: deepcopy on read and on write into the cache.

        This test deliberately does NOT clear ``_features_cache``
        before reloading — that's the whole point. The cache must
        stay clean on its own.
        """
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_strand()  # cycle 1 → -1 in memory
            await pilot.pause()
            app.screen.action_close()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.UnsavedQuitModal)
            app.screen.query_one("#btn-abandon", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Note: NO sc._features_cache = None reset.
            cached = sc._load_features()
            assert cached[0]["strand"] == 1, (
                f"After Abandon, _features_cache still carried the "
                f"in-memory mutation; expected strand=1, got "
                f"strand={cached[0]['strand']}. Cache leaked via "
                f"shared dict refs."
            )

    async def test_post_save_mutation_does_not_poison_cache(self, tiny_record):
        """A trickier variant: user saves once, then makes more edits,
        then abandons. The post-save edits must NOT survive in the
        cache. Pre-fix, ``_save_features`` did
        ``_features_cache = list(entries)`` which shared dict refs
        with the caller's list — so any mutation of the caller's list
        after the save also mutated the cache. Fix: deepcopy on save.
        """
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            # First mutation, save it.
            app.screen.action_strand()  # 1 → -1
            await pilot.pause()
            app.screen.action_save()
            await pilot.pause()
            # Disk + cache should both be -1 now.
            assert sc._load_features()[0]["strand"] == -1
            # Second mutation, abandon it.
            app.screen.action_strand()  # -1 → 0 (in memory only)
            await pilot.pause()
            app.screen.action_close()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.UnsavedQuitModal)
            app.screen.query_one("#btn-abandon", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # NO cache reset — the cache must already be clean.
            cached = sc._load_features()
            assert cached[0]["strand"] == -1, (
                f"Post-save mutation leaked into cache through shared "
                f"dict refs from _save_features; expected the saved "
                f"value (-1), got {cached[0]['strand']}."
            )

    async def test_edit_replaces_entry_and_marks_dirty(self, tiny_record):
        """The Edit button pre-fills AddFeatureModal with the current
        entry; saving from the modal replaces the entry at its current
        index and tags it dirty (no asterisk on disk until save).
        """
        sc._save_features([{
            "name": "lacZ", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_edit()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.AddFeatureModal)
            # Modal should be pre-filled with lacZ.
            name_input = app.screen.query_one("#addfeat-name", Input)
            assert name_input.value == "lacZ"
            # Change the name and save.
            name_input.value = "lacZ-edited"
            app.screen.query_one("#btn-addfeat-save", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Back on the library screen.
            assert isinstance(app.screen, sc.FeatureLibraryScreen)
            assert app.screen._entries[0]["name"] == "lacZ-edited"
            assert 0 in app.screen._dirty_indices
            assert app.screen._has_pending_changes is True
            # Disk untouched.
            sc._features_cache = None
            assert sc._load_features()[0]["name"] == "lacZ"


# ═══════════════════════════════════════════════════════════════════════════════
# ColorPickerModal
# ═══════════════════════════════════════════════════════════════════════════════

class TestColorPickerModal:
    """The color picker must mount, dismiss with the correct dict shape,
    and guard against "set as default" without a chosen color."""

    async def test_mount_at_baseline(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.ColorPickerModal("CDS", "#123456"))
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.ColorPickerModal)

    async def test_cancel_returns_none(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        captured = []
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(
                sc.ColorPickerModal("CDS", None),
                callback=captured.append,
            )
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.action_cancel()
            await pilot.pause()
            await pilot.pause(0.05)
        assert captured == [None]

    async def test_save_returns_color_dict(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        captured = []
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(
                sc.ColorPickerModal("CDS", "#FF0000"),
                callback=captured.append,
            )
            await pilot.pause()
            await pilot.pause(0.1)
            btn = app.screen.query_one("#btn-colorpick-save", Button)
            btn.press()
            await pilot.pause()
            await pilot.pause(0.05)
        assert captured and captured[0] == {"color": "#FF0000",
                                            "set_default": False}

    async def test_auto_returns_none_color(self, tiny_record):
        """Clicking Auto clears the override (color: None)."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        captured = []
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(
                sc.ColorPickerModal("CDS", "#FF0000"),
                callback=captured.append,
            )
            await pilot.pause()
            await pilot.pause(0.1)
            # Press Auto first (sets pending -> None)
            app.screen.query_one("#btn-colorpick-auto", Button).press()
            await pilot.pause()
            # Then press Save
            app.screen.query_one("#btn-colorpick-save", Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
        assert captured and captured[0]["color"] is None

    async def test_set_as_default_requires_color(self, tiny_record):
        """Clicking 'Save + set as type default' without a color picked
        must NOT dismiss — user needs to pick a color first."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        captured = []
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(
                sc.ColorPickerModal("CDS", None),
                callback=captured.append,
            )
            await pilot.pause()
            await pilot.pause(0.1)
            # Explicitly set pending to None via Auto
            app.screen.query_one("#btn-colorpick-auto", Button).press()
            await pilot.pause()
            app.screen.query_one("#btn-colorpick-default", Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
            # Still open — dismissal blocked
            assert isinstance(app.screen, sc.ColorPickerModal)
        assert captured == []   # never fired


# ═══════════════════════════════════════════════════════════════════════════════
# Snippet panel rendering
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatureSnippetPanelFormat:
    """``_FeatureSnippetPanel`` splits its output into three pieces:
    ``_format_header`` (markup with name/type/strand/length/color),
    ``_render_dna`` (a ``Rich.Text`` double-stranded DNA block that goes
    through the shared ``_build_seq_text`` pipeline — the same renderer
    the main SequencePanel uses), and ``_format_qualifiers`` (markup)."""

    # ── Header: name / type / strand / length / color swatch ──────────────────

    def test_header_renders_length_in_bp(self):
        p = sc._FeatureSnippetPanel()
        out = p._format_header({
            "name": "x", "feature_type": "CDS",
            "sequence": "ATGACC", "strand": 1,
        })
        assert "6 bp" in out

    def test_header_shows_color(self):
        p = sc._FeatureSnippetPanel()
        out = p._format_header({
            "name": "x", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1, "color": "#ABCDEF",
        })
        assert "#ABCDEF" in out

    def test_header_shows_strand_label(self):
        """Each strand value picks a distinct human-readable tag so the
        user knows which direction the snippet's arrow will render in."""
        p = sc._FeatureSnippetPanel()
        base = {"name": "x", "feature_type": "CDS", "sequence": "ATG"}
        assert "forward"   in p._format_header({**base, "strand":  1})
        assert "reverse"   in p._format_header({**base, "strand": -1})
        assert "arrowless" in p._format_header({**base, "strand":  0})
        assert "double"    in p._format_header({**base, "strand":  2})

    # ── Qualifiers ────────────────────────────────────────────────────────────

    def test_qualifiers_rendered(self):
        p = sc._FeatureSnippetPanel()
        out = p._format_qualifiers({
            "qualifiers": {"gene": ["lacZ"], "note": ["hello"]},
        })
        assert "gene" in out and "lacZ" in out
        assert "note" in out and "hello" in out

    def test_no_qualifiers_shows_placeholder(self):
        p = sc._FeatureSnippetPanel()
        out = p._format_qualifiers({"qualifiers": {}})
        assert "No qualifiers" in out

    # ── DNA block: ``_build_seq_text`` pipeline (dithered ▒ bar + arrow) ──────
    # The DNA Rich Text contains the feature-bar row drawn by
    # ``_render_feature_row_pair``. Its glyph set depends on strand:
    #   strand=1  → ▒…▒▶
    #   strand=-1 → ◀▒…▒
    #   strand=0  → ▒…▒       (no arrows)
    #   strand=2  → ◀▒…▒▶
    # "ATG" is 3 bp, so the 1-bp ▲/▼ shortcut in the bar renderer never
    # fires here.

    def test_dna_forward_strand_has_right_arrow(self):
        p = sc._FeatureSnippetPanel()
        out = p._render_dna({
            "name": "fwd", "feature_type": "CDS",
            "sequence": "ATG", "strand": 1,
        }).plain
        assert "▶" in out
        assert "◀" not in out

    def test_dna_reverse_strand_has_left_arrow(self):
        p = sc._FeatureSnippetPanel()
        out = p._render_dna({
            "name": "rev", "feature_type": "CDS",
            "sequence": "ATG", "strand": -1,
        }).plain
        assert "◀" in out
        assert "▶" not in out

    def test_dna_arrowless_strand_has_no_arrow(self):
        p = sc._FeatureSnippetPanel()
        out = p._render_dna({
            "name": "flat", "feature_type": "rep_origin",
            "sequence": "GCA", "strand": 0,
        }).plain
        assert "▶" not in out
        assert "◀" not in out
        assert "▒" in out

    def test_dna_double_strand_has_both_arrows(self):
        """strand == 2 puts ◀ at the start AND ▶ at the end of the bar."""
        p = sc._FeatureSnippetPanel()
        out = p._render_dna({
            "name": "two", "feature_type": "CDS",
            "sequence": "ATGACC", "strand": 2,
        }).plain
        assert "◀" in out
        assert "▶" in out

    def test_dna_empty_sequence_handled(self):
        """Empty sequence must render a placeholder instead of calling
        through to ``_build_seq_text`` (which has undefined behavior on
        a zero-length seq)."""
        p = sc._FeatureSnippetPanel()
        out = p._render_dna({
            "name": "empty", "feature_type": "CDS",
            "sequence": "", "strand": 1,
        }).plain
        assert "no sequence" in out


# ═══════════════════════════════════════════════════════════════════════════════
# AddFeatureModal — 4-way Orientation RadioSet
# ═══════════════════════════════════════════════════════════════════════════════

class TestAddFeatureModalOrientation:
    """The modal's Strand row was renamed to Orientation and expanded from 2
    radios to 4 (Forward / Reverse / Arrowless / Double) so the library can
    represent strand values 1 / -1 / 0 / 2 everywhere — matching the Cycle-
    Strand button on the Features workbench."""

    async def test_label_reads_orientation_not_strand(self, tiny_record):
        """Regression guard: the label text must say 'Orientation:' — the
        user renamed it to match the 4-way arrow semantics, and keeping the
        word 'Strand' here would confuse callers who only know the biological
        meaning."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal(have_cursor=False))
            await pilot.pause()
            await pilot.pause(0.1)
            labels = [str(lbl.render()) for lbl
                      in app.screen.query(Label).results()]
            assert any("Orientation" in t for t in labels), (
                f"Expected an 'Orientation:' label; saw {labels}"
            )
            for t in labels:
                # Strip Content()/markup wrappers, then normalise
                stripped = t.strip().rstrip(":")
                assert "strand" not in stripped.lower(), (
                    f"'Strand' must no longer appear as a label; saw {t}"
                )

    async def test_all_four_radios_present(self, tiny_record):
        """Each orientation needs its own stable id so prefill and gather
        logic can target the right radio without index arithmetic."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal(have_cursor=False))
            await pilot.pause()
            await pilot.pause(0.1)
            for rid in ("#addfeat-strand-fwd", "#addfeat-strand-rev",
                        "#addfeat-strand-none", "#addfeat-strand-both"):
                app.screen.query_one(rid, RadioButton)

    @pytest.mark.parametrize("rid,expected_strand", [
        ("#addfeat-strand-fwd",  1),
        ("#addfeat-strand-rev", -1),
        ("#addfeat-strand-none", 0),
        ("#addfeat-strand-both", 2),
    ])
    async def test_save_returns_matching_strand(self, rid, expected_strand,
                                                tiny_record):
        """Pressing Save while a given orientation radio is lit must emit a
        `strand` matching that radio's semantic value — 1/-1/0/2."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        captured: list = []
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(
                sc.AddFeatureModal(have_cursor=False),
                callback=captured.append,
            )
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.query_one("#addfeat-name", Input).value = "x"
            app.screen.query_one("#addfeat-seq", TextArea).text = "ATGACC"
            # Light up the target radio; RadioSet clears siblings on value=True.
            app.screen.query_one(rid, RadioButton).value = True
            await pilot.pause()
            app.screen.query_one("#btn-addfeat-save", Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
        assert captured, "modal never dismissed"
        assert captured[0]["entry"]["strand"] == expected_strand

    @pytest.mark.parametrize("prefill_strand,expected_lit_radio", [
        (1,  "#addfeat-strand-fwd"),
        (-1, "#addfeat-strand-rev"),
        (0,  "#addfeat-strand-none"),
        (2,  "#addfeat-strand-both"),
    ])
    async def test_prefill_lights_correct_radio(self, prefill_strand,
                                                expected_lit_radio,
                                                tiny_record):
        """Opening the modal with `prefill={'strand': X}` must light exactly
        the radio that matches X. Guards the Ctrl+Shift+F flow where a captured
        feature's strand drives which orientation the user sees selected."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal(
                prefill={"name": "x", "feature_type": "CDS",
                         "sequence": "ATG", "strand": prefill_strand},
                have_cursor=False,
            ))
            await pilot.pause()
            await pilot.pause(0.1)
            for rid in ("#addfeat-strand-fwd", "#addfeat-strand-rev",
                        "#addfeat-strand-none", "#addfeat-strand-both"):
                rb = app.screen.query_one(rid, RadioButton)
                if rid == expected_lit_radio:
                    assert rb.value, f"{rid} should be lit for strand={prefill_strand}"
                else:
                    assert not rb.value, (
                        f"{rid} must be off when {expected_lit_radio} is lit"
                    )


# ═══════════════════════════════════════════════════════════════════════════════
# Ctrl+Shift+F: capture current selection / highlighted feature → AddFeatureModal
# ═══════════════════════════════════════════════════════════════════════════════

class TestCaptureToFeatures:
    """Ctrl+Shift+F (``action_capture_to_features``) should:
      1. Read either ``sp._user_sel`` (drag selection) or the highlighted
         feature from the plasmid map and open AddFeatureModal prefilled.
      2. On Save, persist the entry to features.json (via the same
         ``_persist_feature_entry`` helper used by the regular Add path).
      3. Push FeatureLibraryScreen so the user lands in the workbench with
         their new entry visible — the user's stated intent."""

    async def test_empty_state_does_not_open_modal(self, tiny_record):
        """No selection + no highlighted feature → notify, no modal push.
        Opening a blank modal would force the user to retype the name,
        type, AND sequence from scratch — defeating the shortcut's purpose."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            # Ensure nothing is selected and no feature highlighted.
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp._user_sel = None
            pm.selected_idx = -1
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.05)
            assert not isinstance(app.screen, sc.AddFeatureModal)

    async def test_drag_selection_prefills_raw_sequence(self, tiny_record):
        """With a Shift+drag region set, the modal should open prefilled
        with the sliced DNA, a blank name, and feature_type=misc_feature."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm.selected_idx = -1
            sp._user_sel = (3, 9)
            expected_sub = sp._seq[3:9].upper()
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.AddFeatureModal)
            seq_ta = app.screen.query_one("#addfeat-seq", TextArea)
            name   = app.screen.query_one("#addfeat-name", Input).value
            ftype  = app.screen.query_one("#addfeat-type", sc.Select).value
            assert seq_ta.text.strip().upper() == expected_sub
            assert name == "", "drag-capture leaves the name blank"
            assert ftype == "misc_feature"

    async def test_highlighted_feature_prefills_all_fields(self, tiny_record):
        """With a feature selected on the map, the modal opens prefilled
        with that feature's label / type / sequence / strand. ``tiny_record``
        feature #0 is a forward CDS at [0, 27); feature #1 is a reverse
        misc_feature at [50, 80)."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp._user_sel = None
            # Select the forward CDS
            pm.selected_idx = next(i for i, f in enumerate(pm._feats)
                                   if f.get("type") == "CDS")
            feat = pm._feats[pm.selected_idx]
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.AddFeatureModal)
            # Sequence should match the forward-strand slice (+1 strand)
            seq_ta = app.screen.query_one("#addfeat-seq", TextArea)
            assert seq_ta.text.upper() == sp._seq[feat["start"]:feat["end"]].upper()
            ftype = app.screen.query_one("#addfeat-type", sc.Select).value
            assert ftype == "CDS"
            # Forward radio lit
            assert app.screen.query_one("#addfeat-strand-fwd",
                                        RadioButton).value

    async def test_reverse_feature_stores_revcomp_as_sequence(self, tiny_record):
        """Reverse-strand features store the 5'→3' of the feature as read —
        i.e. the reverse-complement of the genomic slice — matching what
        AddFeatureModal's Insert path expects so a round-trip through the
        library preserves biology."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp._user_sel = None
            pm.selected_idx = next(i for i, f in enumerate(pm._feats)
                                   if f.get("strand") == -1)
            feat    = pm._feats[pm.selected_idx]
            fwd_gen = sp._seq[feat["start"]:feat["end"]].upper()
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            seq_ta = app.screen.query_one("#addfeat-seq", TextArea)
            assert seq_ta.text.upper() == sc._rc(fwd_gen)
            # Reverse radio lit
            assert app.screen.query_one("#addfeat-strand-rev",
                                        RadioButton).value

    async def test_save_persists_and_jumps_to_library_screen(self, tiny_record):
        """The capture flow's raison d'être: Save writes to features.json
        AND the user lands on FeatureLibraryScreen so they see their new
        entry. This is the user-visible contract of the shortcut."""
        sc._save_features([])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm.selected_idx = -1
            sp._user_sel = (0, 6)
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.query_one("#addfeat-name", Input).value = "captured-1"
            app.screen.query_one("#btn-addfeat-save", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            # Must assert INSIDE the context — app shuts down on exit.
            assert isinstance(app.screen, sc.FeatureLibraryScreen)
            sc._features_cache = None
            entries = sc._load_features()
            assert len(entries) == 1
            assert entries[0]["name"] == "captured-1"

    async def test_selection_takes_priority_over_highlighted_feature(
            self, tiny_record):
        """If BOTH a drag selection and a highlighted feature exist at the
        same time, the raw drag selection wins. Rationale: the drag is the
        user's most recent explicit gesture, and it may span feature
        boundaries (so re-deriving from the feature would drop bases)."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm.selected_idx = next(i for i, f in enumerate(pm._feats)
                                   if f.get("type") == "CDS")
            sp._user_sel = (40, 46)
            expected_sub = sp._seq[40:46].upper()
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            seq_ta = app.screen.query_one("#addfeat-seq", TextArea)
            assert seq_ta.text.upper() == expected_sub
            # Name is blank (drag path), not the feature's label
            assert app.screen.query_one("#addfeat-name", Input).value == ""

    async def test_resite_overlay_is_rejected(self, tiny_record):
        """A highlighted restriction-site overlay (``type == "resite"``) is
        not a real feature; capturing one would produce a nonsense library
        entry. The shortcut must refuse and notify."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp._user_sel = None
            # Inject a fake resite into pm._feats and point selected_idx at it.
            pm._feats = list(pm._feats) + [{
                "type": "resite", "start": 0, "end": 6, "strand": 1,
                "color": "red", "label": "EcoRI",
            }]
            pm.selected_idx = len(pm._feats) - 1
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.05)
            assert not isinstance(app.screen, sc.AddFeatureModal)


# ═══════════════════════════════════════════════════════════════════════════════
# Ctrl+Shift+F: drag selection that matches a feature's exact range
# ═══════════════════════════════════════════════════════════════════════════════

class TestCaptureDragMatchesFeature:
    """If the drag selection's (start, end) exactly matches a real feature's
    range, the prefill carries that feature's full metadata — not the
    generic misc_feature defaults. Regression guard for the 2026-04-20
    enhancement where sidebar-click (which sets both ``_user_sel`` and
    ``selected_idx``) must produce the same rich prefill as the
    highlighted-feature path."""

    async def test_drag_matching_cds_carries_type_strand_name(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            cds_idx = next(i for i, f in enumerate(pm._feats)
                           if f.get("type") == "CDS")
            feat = pm._feats[cds_idx]
            pm.selected_idx = -1
            sp._user_sel = (feat["start"], feat["end"])
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.AddFeatureModal)
            ftype = app.screen.query_one("#addfeat-type", sc.Select).value
            name  = app.screen.query_one("#addfeat-name", Input).value
            assert ftype == "CDS"
            # Feature's label (or type fallback) carries into the name field
            assert name == (feat.get("label") or "CDS")
            # Forward radio still lit
            assert app.screen.query_one("#addfeat-strand-fwd",
                                        RadioButton).value

    async def test_drag_matching_reverse_feature_stores_revcomp(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            rev_idx = next(i for i, f in enumerate(pm._feats)
                           if f.get("strand") == -1)
            feat = pm._feats[rev_idx]
            pm.selected_idx = -1
            sp._user_sel = (feat["start"], feat["end"])
            fwd_gen = sp._seq[feat["start"]:feat["end"]].upper()
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            seq_ta = app.screen.query_one("#addfeat-seq", TextArea)
            assert seq_ta.text.upper() == sc._rc(fwd_gen)
            assert app.screen.query_one("#addfeat-strand-rev",
                                        RadioButton).value

    async def test_drag_not_matching_any_feature_uses_generic_defaults(
            self, tiny_record):
        """If the drag doesn't coincide with a feature, the old behavior
        still applies — blank name, misc_feature, strand=1. Otherwise
        dragging a random region would surprise the user with a
        leaking feature type from the last feature they clicked."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pm.selected_idx = -1
            # Pick a range that doesn't match any (start, end) pair in _feats.
            sp._user_sel = (7, 13)
            assert not any((f["start"], f["end"]) == (7, 13)
                           for f in pm._feats)
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            ftype = app.screen.query_one("#addfeat-type", sc.Select).value
            name  = app.screen.query_one("#addfeat-name", Input).value
            assert ftype == "misc_feature"
            assert name == ""


# ═══════════════════════════════════════════════════════════════════════════════
# AddFeatureModal: Color field
# ═══════════════════════════════════════════════════════════════════════════════

class TestAddFeatureModalColorField:
    """The Add Feature modal gained a Color field so captured-feature
    colors (from Ctrl+Shift+F) survive through Save, and manual entries can set
    a per-entry color without a round-trip through FeatureLibraryScreen."""

    async def test_color_prefill_round_trips_through_save(self, tiny_record):
        """Prefilled color must come back out of _gather unchanged."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AddFeatureModal(prefill={
                "name": "widget", "feature_type": "CDS",
                "sequence": "ATG", "strand": 1,
                "color": "#AABBCC",
            }, have_cursor=False)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            assert modal._color == "#AABBCC"
            entry = modal._gather()
            assert entry is not None
            assert entry["color"] == "#AABBCC"

    async def test_no_color_prefill_defaults_to_none(self, tiny_record):
        """Without a color in the prefill, _gather returns color=None so the
        library entry falls through to the type default."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AddFeatureModal(prefill={
                "name": "widget", "feature_type": "CDS",
                "sequence": "ATG", "strand": 1,
            }, have_cursor=False)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            assert modal._color is None
            entry = modal._gather()
            assert entry is not None
            assert entry["color"] is None

    async def test_auto_button_clears_prefilled_color(self, tiny_record):
        """Clicking 'Auto' while a color is set clears the override back to
        None so the render falls through to type default."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.AddFeatureModal(prefill={
                "name": "widget", "feature_type": "CDS",
                "sequence": "ATG", "strand": 1,
                "color": "#112233",
            }, have_cursor=False)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            assert modal._color == "#112233"
            modal.query_one("#btn-addfeat-color-clear", Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
            assert modal._color is None

    async def test_capture_threads_color_into_addfeature_modal(self, tiny_record):
        """Full path: Ctrl+Shift+F on a highlighted feature should carry the
        feature's ``color`` field (from ``_feats``) into the modal's
        ``_color``. Palette-format ``color(N)`` values get normalised to
        the equivalent hex so the stored library entry and every
        downstream preview can use plain markup."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            sp = app.query_one("#seq-panel", sc.SequencePanel)
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            sp._user_sel = None
            pm.selected_idx = next(i for i, f in enumerate(pm._feats)
                                   if f.get("type") == "CDS")
            raw_color = pm._feats[pm.selected_idx].get("color")
            expected = sc._normalise_color_input(raw_color) or raw_color
            app.action_capture_to_features()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.AddFeatureModal)
            assert app.screen._color == expected


# ═══════════════════════════════════════════════════════════════════════════════
# Color input helpers (pure functions — no Textual needed)
# ═══════════════════════════════════════════════════════════════════════════════

class TestNormaliseColorInput:
    """``_normalise_color_input`` accepts user-typed colors in multiple
    formats and returns a canonical ``#RRGGBB`` string. Any form the user
    might plausibly type for a terminal color must round-trip cleanly."""

    @pytest.mark.parametrize("raw,expected", [
        ("#FF6347", "#FF6347"),
        ("#ff6347", "#FF6347"),
        ("#F63",    "#FF6633"),
        ("#f63",    "#FF6633"),
        ("208",     "#FF8700"),
        ("0",       "#000000"),
        ("255",     "#EEEEEE"),
        ("color(39)", "#00AFFF"),
    ])
    def test_valid_inputs_produce_canonical_hex(self, raw, expected):
        assert sc._normalise_color_input(raw) == expected

    @pytest.mark.parametrize("raw", [
        "", "   ", "not-a-color", "#ZZZ", "#1234567", "#12345",
        "256", "-1", "color()", "rgb(1,2,3)",
    ])
    def test_invalid_inputs_return_none(self, raw):
        assert sc._normalise_color_input(raw) is None


class TestXtermIndexToHex:
    """``_xterm_index_to_hex`` maps xterm-256 indices to their spec RGB
    value. The three regions (0-15 ANSI, 16-231 cube, 232-255 grayscale)
    each have their own formula; this test pins a representative cell in
    each so future palette tweaks fail loudly."""

    def test_ansi_region(self):
        assert sc._xterm_index_to_hex(0)  == "#000000"
        assert sc._xterm_index_to_hex(15) == "#FFFFFF"

    def test_cube_region_red_corner(self):
        assert sc._xterm_index_to_hex(196) == "#FF0000"

    def test_cube_region_origin(self):
        # Index 16 is the base of the 6×6×6 cube — all channels at 0.
        assert sc._xterm_index_to_hex(16) == "#000000"

    def test_grayscale_region(self):
        # First gray (idx 232) is 0x08, last (idx 255) is 0xEE per spec.
        assert sc._xterm_index_to_hex(232) == "#080808"
        assert sc._xterm_index_to_hex(255) == "#EEEEEE"

    def test_out_of_range_clamps(self):
        # Negative clamps to 0 (black); >255 clamps to 255 (final gray).
        assert sc._xterm_index_to_hex(-5)  == "#000000"
        assert sc._xterm_index_to_hex(999) == "#EEEEEE"


# ═══════════════════════════════════════════════════════════════════════════════
# ColorPickerModal: expanded full-palette picker
# ═══════════════════════════════════════════════════════════════════════════════

class TestColorPickerExpanded:
    """Regression guards for the 2026-04-20 ColorPickerModal rework:
    the curated swatch grid, the xterm 256 grid, and the custom hex
    input must all drive the same ``_pending`` state and dismiss with
    an uppercase canonical hex string."""

    async def test_xterm_cell_click_sets_pending(self, tiny_record):
        """Clicking an xterm cell button loads that index's RGB into the
        pending color — no round-trip through the input field."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            # Simulate pressing xterm cell 196 (a bright red) directly.
            modal.query_one("#colorpick-x-196", Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
            assert modal._pending == "#FF0000"

    async def test_custom_hex_apply_sets_pending(self, tiny_record):
        """Typing a hex into the Input + clicking Apply updates _pending
        to the canonicalised uppercase form."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            modal.query_one("#colorpick-hex-input", Input).value = "#f63"
            modal.query_one("#btn-colorpick-apply", Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
            assert modal._pending == "#FF6633"

    async def test_invalid_hex_leaves_pending_unchanged(self, tiny_record):
        """If the user types a garbage string and clicks Apply, _pending
        must NOT change and the status bar surfaces a red error. This is
        the user-visible safety net for typo input."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", "#FF6347")
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            assert modal._pending == "#FF6347"
            modal.query_one("#colorpick-hex-input", Input).value = "not-a-color"
            modal.query_one("#btn-colorpick-apply", Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
            assert modal._pending == "#FF6347"

    async def test_xterm_index_applied_via_custom_input(self, tiny_record):
        """xterm indices typed into the custom input are accepted and
        converted to their hex equivalent — letting power users type '196'
        instead of hunting for the red corner in the grid."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            modal.query_one("#colorpick-hex-input", Input).value = "196"
            modal.query_one("#btn-colorpick-apply", Button).press()
            await pilot.pause()
            await pilot.pause(0.05)
            assert modal._pending == "#FF0000"

    async def test_capability_warning_mounts(self, tiny_record):
        """The capability warning widget exists and gets populated on
        mount — even if a test terminal returns None for color_system, we
        render a human-readable label rather than crashing. The test just
        checks that the widget was rendered (non-zero size) — the exact
        text depends on the terminal's reported color_system, which
        varies between CI environments."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            cap = modal.query_one("#colorpick-capability", Static)
            # The widget must be mounted; _refresh_capability_warning must
            # not have raised. Its exact rendered text is terminal-dependent.
            assert cap.is_mounted
            assert cap.region.width > 0


# ═══════════════════════════════════════════════════════════════════════════════
# ColorPickerModal — drag-to-preview across xterm cells
# ═══════════════════════════════════════════════════════════════════════════════

class TestColorPickerDragPreview:
    """The xterm 256 grid supports click-and-drag: holding the mouse button
    down and sweeping across cells updates the preview live, committing on
    release. These tests drive the handlers directly with stub events —
    Textual's pilot doesn't expose a clean mouse-drag primitive, so we go
    one level below the event loop."""

    def _make_event(self, sx: int, sy: int, button: int = 1):
        """Stub mouse event — the handlers only read ``screen_x``,
        ``screen_y``, and ``button``, so a SimpleNamespace is enough."""
        from types import SimpleNamespace
        return SimpleNamespace(screen_x=sx, screen_y=sy, button=button)

    async def test_mouse_down_on_cell_arms_drag_and_sets_pending(
            self, tiny_record):
        """MouseDown on an xterm cell should (a) flip `_drag_active` on and
        (b) immediately load that cell's color into the preview."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)

            cell = modal.query_one("#colorpick-x-196", Button)
            r = cell.region
            modal.on_mouse_down(self._make_event(
                r.x + r.width // 2, r.y + r.height // 2))
            assert modal._drag_active is True
            assert modal._pending == sc._xterm_index_to_hex(196)

    async def test_mouse_move_during_drag_updates_pending(self, tiny_record):
        """Sweeping across cells while holding the button down should
        rewrite the pending color every time the cursor enters a new cell."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)

            start = modal.query_one("#colorpick-x-21",  Button)  # blue cube
            over  = modal.query_one("#colorpick-x-82",  Button)  # green
            end   = modal.query_one("#colorpick-x-196", Button)  # red

            modal.on_mouse_down(self._make_event(
                start.region.x + 1, start.region.y))
            assert modal._pending == sc._xterm_index_to_hex(21)

            modal.on_mouse_move(self._make_event(
                over.region.x + 1, over.region.y))
            assert modal._pending == sc._xterm_index_to_hex(82)

            modal.on_mouse_move(self._make_event(
                end.region.x + 1, end.region.y))
            assert modal._pending == sc._xterm_index_to_hex(196)

    async def test_mouse_move_without_drag_is_noop(self, tiny_record):
        """Plain mouse-move (no drag active) must not change the preview —
        otherwise hovering over the grid would silently rewrite the user's
        choice."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", "#FF6347")
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)

            cell = modal.query_one("#colorpick-x-21", Button)
            assert modal._drag_active is False
            modal.on_mouse_move(self._make_event(
                cell.region.x + 1, cell.region.y))
            assert modal._pending == "#FF6347"

    async def test_mouse_up_clears_drag_flag(self, tiny_record):
        """Releasing the button must disarm drag-mode so later hovers don't
        keep painting the preview."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)

            start = modal.query_one("#colorpick-x-82",  Button)
            after = modal.query_one("#colorpick-x-196", Button)

            modal.on_mouse_down(self._make_event(
                start.region.x + 1, start.region.y))
            modal.on_mouse_up(self._make_event(
                start.region.x + 1, start.region.y))
            assert modal._drag_active is False
            before = modal._pending

            modal.on_mouse_move(self._make_event(
                after.region.x + 1, after.region.y))
            assert modal._pending == before

    async def test_mouse_down_outside_grid_does_not_arm_drag(
            self, tiny_record):
        """MouseDown landing anywhere but an xterm cell (e.g. on a Save
        button) must leave drag-mode disarmed so scrolling/clicking the
        non-grid UI keeps working normally."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)

            save = modal.query_one("#btn-colorpick-save", Button)
            modal.on_mouse_down(self._make_event(
                save.region.x + 1, save.region.y))
            assert modal._drag_active is False
            assert modal._pending is None

    async def test_non_left_button_mouse_down_is_ignored(self, tiny_record):
        """Right-click / middle-click on a cell must NOT arm drag-mode —
        only the primary (left) button should drive the drag workflow."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)

            cell = modal.query_one("#colorpick-x-196", Button)
            modal.on_mouse_down(self._make_event(
                cell.region.x + 1, cell.region.y, button=3))  # right click
            assert modal._drag_active is False
            assert modal._pending is None

    async def test_preview_swatch_repaints_when_pending_changes(
            self, tiny_record):
        """The dedicated big preview swatch must update its background
        color every time `_set_pending` fires — that's the user's only
        visual feedback during a drag."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ColorPickerModal("CDS", None)
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)

            swatch = modal.query_one("#colorpick-preview-swatch", Static)
            modal._set_pending("#123456")
            await pilot.pause()
            assert swatch.styles.background.hex.upper() == "#123456"

            modal._set_pending(None)
            await pilot.pause()
            # Clearing selection returns the swatch to transparent —
            # Textual represents that as Color(0, 0, 0, a=0).
            bg = swatch.styles.background
            assert bg.a == 0


# ═══════════════════════════════════════════════════════════════════════════════
# Export FASTA button
# ═══════════════════════════════════════════════════════════════════════════════

class TestFeatureLibraryExportFasta:
    """`FeatureLibraryScreen` grew an "Export FASTA…" button that writes
    the selected entry's sequence to a single-record FASTA file via
    `FastaExportModal`. Entries without a sequence must be rejected with
    a notification rather than opening an empty export modal."""

    async def test_export_button_present(self, tiny_record):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            assert screen.query_one("#btn-flib-export-fasta",
                                    Button) is not None

    async def test_export_pushes_fasta_modal_with_selected_entry(self, tiny_record):
        sc._save_features([{
            "name":         "lacZ",
            "feature_type": "CDS",
            "strand":       1,
            "sequence":     "ATGCATGCATGCATGC",
            "color":        None,
            "qualifiers":   {},
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            tbl = screen.query_one("#flib-table", DataTable)
            tbl.move_cursor(row=0)
            await pilot.pause()
            screen.query_one("#btn-flib-export-fasta", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            top = app.screen
            assert isinstance(top, sc.FastaExportModal)
            assert top._name == "lacZ"
            assert top._sequence == "ATGCATGCATGCATGC"

    async def test_export_with_empty_library_warns(self, tiny_record):
        """No entries → button press should notify and not push the export modal."""
        sc._save_features([])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            screen.query_one("#btn-flib-export-fasta", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.FeatureLibraryScreen)

    async def test_export_entry_without_sequence_warns(self, tiny_record):
        """Entries with an empty sequence field are rejected — no export modal."""
        sc._save_features([{
            "name":         "empty-entry",
            "feature_type": "misc_feature",
            "strand":       1,
            "sequence":     "",
            "color":        None,
            "qualifiers":   {},
        }])
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.FeatureLibraryScreen())
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            tbl = screen.query_one("#flib-table", DataTable)
            tbl.move_cursor(row=0)
            await pilot.pause()
            screen.query_one("#btn-flib-export-fasta", Button).press()
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.FeatureLibraryScreen)
