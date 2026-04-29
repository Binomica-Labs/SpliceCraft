"""
test_modal_boundaries — every modal's visible widgets stay inside the terminal.

**Why this exists**: on 2026-04-20 the AddFeatureModal's type/strand row
stretched to ~39 rows (child Verticals inferred ``height: 1fr`` inside a
Horizontal), which pushed the sequence TextArea to ``y=51`` on a 48-row
terminal — off screen. The fix was to pin the row height and move the body
into a ScrollableContainer. This test is the regression guard for that
class of bug.

**Contract**:
  * The modal's **root container** (the box with a border) must fit
    entirely within the terminal viewport.
  * Every descendant widget that is **not** inside a scrolling container
    must fit inside the terminal viewport.
  * Widgets inside a ScrollableContainer / VerticalScroll / HorizontalScroll
    are allowed to extend past the viewport — they'll scroll into view on
    demand.
"""
from __future__ import annotations

import pytest

import splicecraft as sc
from textual.containers import (
    ScrollableContainer, VerticalScroll, HorizontalScroll,
)


# 160×48 is the project's baseline terminal size (matches the rest of the
# test suite's TERMINAL_SIZE). Every modal must fit here; smaller sizes
# are opt-in per modal via TestAddFeatureModalRegressionGuards below.
_BASELINE = (160, 48)

# Extra sizes for AddFeatureModal's explicit regression guards. The
# original bug made the TextArea invisible on 48 rows, so we check a
# handful of realistic small-terminal shapes to keep it boxed in.
_ADDFEAT_REGRESSION_SIZES = [
    (160, 48),
    (120, 40),
    (100, 30),
]


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

_SCROLLING_TYPES = (ScrollableContainer, VerticalScroll, HorizontalScroll)


def _inside_scrolling_ancestor(widget, modal) -> bool:
    """True if any ancestor between `widget` (exclusive) and `modal`
    (exclusive) is a scrolling container. Such widgets may legitimately
    extend past the viewport."""
    node = widget.parent
    while node is not None and node is not modal:
        if isinstance(node, _SCROLLING_TYPES):
            return True
        node = node.parent
    return False


def _assert_widget_in_bounds(widget, term_w: int, term_h: int,
                             label: str) -> None:
    r = widget.region
    # Zero-size widgets are fine; they just don't render.
    if r.width == 0 or r.height == 0:
        return
    assert r.x >= 0, f"{label}: x={r.x} < 0"
    assert r.y >= 0, f"{label}: y={r.y} < 0"
    assert r.x + r.width  <= term_w, (
        f"{label}: x+w={r.x + r.width} > terminal width {term_w}"
    )
    assert r.y + r.height <= term_h, (
        f"{label}: y+h={r.y + r.height} > terminal height {term_h}  "
        f"(widget would render off screen)"
    )


async def _check_modal(app, pilot, modal, term_w: int, term_h: int) -> None:
    """Push the already-pushed modal's descendants through the bounds check.

    `modal` is expected to be the currently-focused screen.
    """
    await pilot.pause()
    await pilot.pause(0.05)

    # 1. The screen itself — Textual Screens fill the terminal, so its
    #    region should equal (0, 0, term_w, term_h). We don't assert
    #    equality (some Textual versions add padding) but the region
    #    must not exceed the terminal.
    _assert_widget_in_bounds(modal, term_w, term_h,
                             f"{type(modal).__name__}(screen)")

    # 2. Every descendant widget: if it's outside any scrolling ancestor,
    #    it must fit. Widgets inside a ScrollableContainer are exempt —
    #    they render into a virtual space that scrolls.
    for w in modal.walk_children(with_self=False):
        if _inside_scrolling_ancestor(w, modal):
            continue
        _assert_widget_in_bounds(
            w, term_w, term_h,
            f"{type(modal).__name__} > {type(w).__name__}#{w.id or '<?>'}",
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Factories: each returns a fresh modal ready to push
# ═══════════════════════════════════════════════════════════════════════════════

def _make_feats() -> list[dict]:
    """A minimal feats list the Mutagenize/Domesticator/PrimerDesign modals
    accept without erroring out during compose."""
    return [{
        "type": "CDS", "start": 0, "end": 30, "strand": 1,
        "label": "testA", "color": "white",
    }]


def _make_record():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq("ATG" * 40),
                    id="TEST", name="TEST",
                    description="boundary-test record")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = "circular"
    return rec


# (label, factory_callable) — factory returns a newly-constructed modal.
# Every modal listed here must be constructible without I/O.
_MODAL_CASES = [
    ("FetchModal",                 lambda: sc.FetchModal()),
    ("OpenFileModal",              lambda: sc.OpenFileModal()),
    ("ExportGenBankModal",         lambda: sc.ExportGenBankModal(_make_record())),
    ("EditSeqDialog.insert",       lambda: sc.EditSeqDialog("insert",   "",  5, 5)),
    ("EditSeqDialog.replace",      lambda: sc.EditSeqDialog("replace", "ATCG", 5, 9)),
    ("UnsavedQuitModal",           lambda: sc.UnsavedQuitModal()),
    ("UnsavedNavigateModal",       lambda: sc.UnsavedNavigateModal(
                                               "go back to collections")),
    ("RenamePlasmidModal",         lambda: sc.RenamePlasmidModal("my_plasmid", "abc123")),
    ("LibraryDeleteConfirmModal",  lambda: sc.LibraryDeleteConfirmModal(
                                               "test", 3000, "abc123")),
    ("PlasmidPickerModal",         lambda: sc.PlasmidPickerModal(None)),
    ("PlasmidFeaturePickerModal",  lambda: sc.PlasmidFeaturePickerModal(
                                               [{"name": "lacZ", "feature_type": "CDS",
                                                 "sequence": "ATG", "strand": 1}],
                                               plasmid_name="demo")),
    ("AddFeatureModal",            lambda: sc.AddFeatureModal(have_cursor=True)),
    ("AminoAcidPickerModal",       lambda: sc.AminoAcidPickerModal(42, "W")),
    ("NcbiTaxonPickerModal",       lambda: sc.NcbiTaxonPickerModal("")),
    ("SpeciesPickerModal",         lambda: sc.SpeciesPickerModal()),
    ("ConstructorModal",           lambda: sc.ConstructorModal()),
    ("DomesticatorModal",          lambda: sc.DomesticatorModal("ATG" * 50,
                                                                 _make_feats())),
    ("MutagenizeModal",            lambda: sc.MutagenizeModal("ATG" * 50,
                                                               _make_feats(),
                                                               plasmid_name="test")),
    ("ColorPickerModal",           lambda: sc.ColorPickerModal("CDS", "#FF0000")),
    ("FastaFilePickerModal",       lambda: sc.FastaFilePickerModal()),
    ("FastaExportModal",           lambda: sc.FastaExportModal(
                                               name="demo",
                                               sequence="ATGCATGC",
                                               default_path="/tmp/demo.fa")),
    ("CollectionsModal",           lambda: sc.CollectionsModal()),
    ("CollectionNameModal",        lambda: sc.CollectionNameModal(
                                               "New collection", "")),
    ("CollectionDeleteConfirmModal", lambda: sc.CollectionDeleteConfirmModal(
                                                  "MyCollection", 5)),
    ("ScaryDeleteConfirmModal",      lambda: sc.ScaryDeleteConfirmModal(
                                                  "MyCollection", 12)),
]


# ═══════════════════════════════════════════════════════════════════════════════
# The boundary check itself
# ═══════════════════════════════════════════════════════════════════════════════

class TestModalBoundaries:
    """Every modal's root container fits in the terminal, and every
    descendant widget outside a scrolling container does too.

    Runs at the project's baseline terminal size (160×48). Smaller
    terminals aren't guaranteed — the Mutagenize and Constructor modals,
    for instance, are full-screen workbenches that can legitimately
    exceed 100×30. Add a dedicated `TestFooModalRegressionGuards` class
    below if a specific modal needs to support smaller terminals."""

    @pytest.mark.parametrize("label,factory", _MODAL_CASES,
                             ids=[c[0] for c in _MODAL_CASES])
    async def test_modal_fits_in_baseline_terminal(
            self, label, factory, tiny_record, isolated_library):
        term_w, term_h = _BASELINE
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = factory()
            app.push_screen(modal)
            await _check_modal(app, pilot, modal, term_w, term_h)


class TestAddFeatureModalRegressionGuards:
    """Explicit regression guards for the AddFeatureModal layout bug —
    the textbox MUST be visible inside the modal body, not past the
    terminal bottom. Exists alongside the parametrized test above so a
    failure points at the AddFeatureModal specifically."""

    @pytest.mark.parametrize("size", _ADDFEAT_REGRESSION_SIZES,
                             ids=[f"{w}x{h}" for (w, h) in
                                  _ADDFEAT_REGRESSION_SIZES])
    async def test_sequence_textarea_inside_scrollable_body(
            self, size, tiny_record, isolated_library):
        """The TextArea is inside a ScrollableContainer, so its absolute
        region may legitimately extend beyond the terminal on tiny screens.
        The important invariant: the **modal box itself** and the
        **button row** (non-scrollable) fit in the terminal."""
        term_w, term_h = size
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal(have_cursor=True))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen

            box = modal.query_one("#addfeat-dlg")
            btns = modal.query_one("#addfeat-btns")
            body = modal.query_one("#addfeat-body")

            _assert_widget_in_bounds(box,  term_w, term_h,
                                     "AddFeatureModal root box")
            _assert_widget_in_bounds(btns, term_w, term_h,
                                     "AddFeatureModal button row")
            _assert_widget_in_bounds(body, term_w, term_h,
                                     "AddFeatureModal scroll body")

    @pytest.mark.parametrize("size", _ADDFEAT_REGRESSION_SIZES,
                             ids=[f"{w}x{h}" for (w, h) in
                                  _ADDFEAT_REGRESSION_SIZES])
    async def test_save_cancel_buttons_reachable(
            self, size, tiny_record, isolated_library):
        """Save and Cancel buttons must always be clickable regardless of
        terminal size — if they're offscreen, the user can't dismiss the
        modal."""
        term_w, term_h = size
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=size) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            app.push_screen(sc.AddFeatureModal(have_cursor=True))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen

            for btn_id in ("#btn-addfeat-save",
                           "#btn-addfeat-insert",
                           "#btn-addfeat-cancel",
                           "#btn-addfeat-import"):
                btn = modal.query_one(btn_id)
                _assert_widget_in_bounds(btn, term_w, term_h,
                                         f"AddFeatureModal {btn_id}")
