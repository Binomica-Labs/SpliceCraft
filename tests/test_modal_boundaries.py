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
    ("BulkExportCollectionModal",  lambda: sc.BulkExportCollectionModal()),
    ("EditSeqDialog.insert",       lambda: sc.EditSeqDialog("insert",   "",  5, 5)),
    ("EditSeqDialog.replace",      lambda: sc.EditSeqDialog("replace", "ATCG", 5, 9)),
    ("UnsavedQuitModal",           lambda: sc.UnsavedQuitModal()),
    ("UnsavedNavigateModal",       lambda: sc.UnsavedNavigateModal(
                                               "go back to collections")),
    ("RenamePlasmidModal",         lambda: sc.RenamePlasmidModal("my_plasmid", "abc123")),
    ("MinPrimerBindingModal",      lambda: sc.MinPrimerBindingModal(15)),
    ("MultiRecordFastaModal",      lambda: sc.MultiRecordFastaModal(
                                               5, default_name="demo_set")),
    ("PartsBinDeleteConfirmModal.single",
                                    lambda: sc.PartsBinDeleteConfirmModal(
                                                ["myCDS"])),
    ("PartsBinDeleteConfirmModal.multi",
                                    lambda: sc.PartsBinDeleteConfirmModal(
                                                ["a", "b", "c", "d", "e"])),
    ("LibraryDeleteConfirmModal",  lambda: sc.LibraryDeleteConfirmModal(
                                               "test", 3000, "abc123")),
    ("UpdateAvailableModal",       lambda: sc.UpdateAvailableModal(
                                               "0.9.1", "0.9.0")),
    ("PlasmidPickerModal",         lambda: sc.PlasmidPickerModal(None)),
    ("PlasmidFeaturePickerModal",  lambda: sc.PlasmidFeaturePickerModal(
                                               [{"name": "lacZ", "feature_type": "CDS",
                                                 "sequence": "ATG", "strand": 1}],
                                               plasmid_name="demo")),
    ("AddFeatureModal",            lambda: sc.AddFeatureModal(selection_range=(0, 10))),
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
    ("GffExportModal",             lambda: sc.GffExportModal(
                                               _make_record(),
                                               default_path="/tmp/demo.gff3")),
    ("ORFFinderModal",             lambda: sc.ORFFinderModal(
                                               "ATGGCCGCCGCCGCCTAA" * 5,
                                               circular=True)),
    ("LibrarySearchModal",         lambda: sc.LibrarySearchModal()),
    ("FeatureSearchModal",         lambda: sc.FeatureSearchModal(
                                               _make_feats(), total=3000)),
    ("LoadPartSourceModal",        lambda: sc.LoadPartSourceModal()),
    ("NamePlasmidModal",           lambda: sc.NamePlasmidModal(
                                               "test_plasmid",
                                               target_label="TU")),
    ("AnnotationTransferModal",    lambda: sc.AnnotationTransferModal(
                                               source_label="src",
                                               target_label="tgt",
                                               transfers=[])),
    ("CollectionsModal",           lambda: sc.CollectionsModal()),
    ("RestoreFromBackupModal",     lambda: sc.RestoreFromBackupModal()),
    ("PrimerDuplicatesModal",      lambda: sc.PrimerDuplicatesModal(
                                               total=671, seq_duplicates=0,
                                               name_collisions=156,
                                               final_kept=515)),
    ("PrimerPlasmidsModal",        lambda: sc.PrimerPlasmidsModal(
                                               primer_entry={
                                                   "name": "P-test-F",
                                                   "sequence": "GATC" * 5,
                                                   "tm": 60.0,
                                               },
                                               usages=[
                                                   {"collection": "TestColl",
                                                    "plasmid_id": "test1",
                                                    "plasmid_name": "pTest1",
                                                    "start": 100, "end": 120,
                                                    "strand": 1},
                                                   {"collection": "TestColl",
                                                    "plasmid_id": "test2",
                                                    "plasmid_name": "pTest2",
                                                    "start": 200, "end": 220,
                                                    "strand": -1},
                                               ])),
    ("CollectionNameModal",        lambda: sc.CollectionNameModal(
                                               "New collection", "")),
    ("AmpliconSaveModal",          lambda: sc.AmpliconSaveModal(
                                               default_name="pUC19_PCR_500bp_1-500",
                                               collections=["Default", "Eden"],
                                               active_collection="Default")),
    ("NewCollectionModal",         lambda: sc.NewCollectionModal()),
    ("CollectionDeleteConfirmModal", lambda: sc.CollectionDeleteConfirmModal(
                                                  "MyCollection", 5)),
    ("ScaryDeleteConfirmModal",      lambda: sc.ScaryDeleteConfirmModal(
                                                  "MyCollection", 12)),
    # Master Delete trio — File → wipe all user data. Three modals
    # in sequence: typed-YES gate, cooldown-gated confirm, single-
    # button result summary. All three must fit at the baseline.
    ("MasterDeleteModal",            lambda: sc.MasterDeleteModal(
                                                  files_count=14,
                                                  dirs_count=4,
                                                  pre_update_present=True)),
    ("MasterDeleteConfirmModal",     lambda: sc.MasterDeleteConfirmModal(
                                                  files_count=14,
                                                  dirs_count=4,
                                                  pre_update_present=True,
                                                  cooldown_s=0.0)),
    ("MasterDeleteResultModal",      lambda: sc.MasterDeleteResultModal({
                                                  "files_removed":     14,
                                                  "dirs_removed":       4,
                                                  "log_files_removed":  2,
                                                  "pre_update_removed": True,
                                                  "residual_files":     0,
                                                  "residual_dirs":      0,
                                                  "errors":             0,
                                              })),
    ("QuitConfirmModal",             lambda: sc.QuitConfirmModal()),
    ("SplashScreen",                 lambda: sc.SplashScreen()),
    ("HelpModal",                    lambda: sc.HelpModal()),
    ("WhatsNewModal",                lambda: sc.WhatsNewModal("0.5.11.0")),
    ("PlasmidStatusPickerModal",     lambda: sc.PlasmidStatusPickerModal(
                                                "pUC19", "VERIFIED")),
    ("NewPlasmidModal",              lambda: sc.NewPlasmidModal()),
    ("BlastModal",                   lambda: sc.BlastModal()),
    ("PlasmidsaurusAlignModal",      lambda: sc.PlasmidsaurusAlignModal()),
    ("MultiAlignPickerModal",        lambda: sc.MultiAlignPickerModal()),
    ("AlignmentManagerModal",        lambda: sc.AlignmentManagerModal([])),
    ("PartsBinPickerModal",          lambda: sc.PartsBinPickerModal()),
    ("ExperimentProjectsPickerModal", lambda: sc.ExperimentProjectsPickerModal()),
    ("ActionsPickerModal",            lambda: sc.ActionsPickerModal()),
    ("ExperimentUnsavedChangesModal", lambda: sc.ExperimentUnsavedChangesModal()),
    ("GelLibraryModal",               lambda: sc.GelLibraryModal()),
    ("LargeFileConfirmModal",        lambda: sc.LargeFileConfirmModal(
                                                "/some/big.gb", "12.3 MB",
                                                threshold_text="cap = 5 MB")),
    ("PrimerEditModal",              lambda: sc.PrimerEditModal(
                                                idx=0,
                                                feat={"label": "P-fwd",
                                                      "type":  "primer_bind",
                                                      "start": 100,
                                                      "end":   108,
                                                      "strand": 1,
                                                      "color": "#00BFFF"},
                                                total=2686,
                                                primer_seq="GAATTCATGAAACG",
                                                notes="See Bolivar 1977",
                                                template="A" * 2686)),
    ("FeatureEditModal",             lambda: sc.FeatureEditModal(
                                                idx=0,
                                                feat={"label": "lacZ",
                                                      "type":  "CDS",
                                                      "start": 100,
                                                      "end":   400,
                                                      "strand": 1,
                                                      "color": "#80FF80"},
                                                total=2686,
                                                sequence="ATG" * 100,
                                                notes="See [Bolivar 1977](https://doi.org/10.1016/0378-1119(77)90000-2)")),
    ("PartEditModal",                lambda: sc.PartEditModal(
                                                idx=0,
                                                part={"name": "myCDS",
                                                      "type": "CDS",
                                                      "position": "Pos 3",
                                                      "oh5": "AATG",
                                                      "oh3": "GCTT",
                                                      "backbone": "pUPD2",
                                                      "marker": "Spectinomycin",
                                                      "sequence": "ATG" * 50,
                                                      "fwd_primer": "GCGCCGTCTCAAATG",
                                                      "rev_primer": "GCGCCGTCTCAAAGC",
                                                      "fwd_tm": 60.5,
                                                      "rev_tm": 61.2,
                                                      "grammar": "gb_l0"})),
    ("HistoryViewerModal",           lambda: sc.HistoryViewerModal(
                                                "pUC19",
                                                sc._CommercialSaaSHistoryNode.new(
                                                    name="pUC19.dna",
                                                    seq_len=2686,
                                                    circular=True,
                                                    operation="insertFragment",
                                                    node_id=0))),
    ("ExportCommercialSaaSModal",          lambda: sc.ExportCommercialSaaSModal(
                                                {"id": "ex", "name": "ex",
                                                 "size": 2686,
                                                 "history_xml":
                                                     "<HistoryTree/>"},
                                                default_path="/tmp/ex.dna")),
    # DropdownScreen is a positioned popup overlay rather than a
    # centered dialog — its `compose` doesn't anchor to the standard
    # 160×48 layout box but it still has to render inside the canvas
    # bounds. Use a 3-item menu (one separator + two actionable rows)
    # at offset (10, 5) — enough to exercise the layout without
    # special-casing the boundary check.
    ("DropdownScreen",             lambda: sc.DropdownScreen(
                                                items=[("Open…", "app.open_file"),
                                                       ("---",   None),
                                                       ("Quit",  "app.quit")],
                                                x=10, y=5)),
    # GrammarEditorModal opens a built-in grammar in read-only mode.
    # Built-ins resolve from `_BUILTIN_GRAMMARS` so the modal builds
    # without a `cloning_grammars.json` on disk (the autouse
    # `_protect_user_data` fixture monkeypatches that path away anyway).
    # `gb_l0` is the default-active grammar so it's the most
    # representative case to pin.
    ("GrammarEditorModal.builtin", lambda: sc.GrammarEditorModal("gb_l0")),
    # Create-mode editor: empty grammar_id flips the modal into the
    # blank-fields path. Important to pin separately because the
    # compose branch differs (no entry-vector pickers, "New Grammar"
    # title, mandatory level-up enzyme input).
    ("GrammarEditorModal.new",     lambda: sc.GrammarEditorModal("")),
    ("GrammarManagerModal",        lambda: sc.GrammarManagerModal()),
    ("EntryVectorsModal",          lambda: sc.EntryVectorsModal("gb_l0")),
    # Enzyme catalog management + custom-enzyme add modals (2026-05-22).
    # Two-pane layout fits comfortably at 140 width; the smaller add
    # modal is single-column at 80.
    ("EnzymeCollectionsModal",     lambda: sc.EnzymeCollectionsModal()),
    ("AddCustomEnzymeModal",       lambda: sc.AddCustomEnzymeModal()),
    ("_EnzymeNamePromptModal",     lambda: sc._EnzymeNamePromptModal(
                                                "New catalog — name?")),
    # SettingsModal consolidates the Settings dropdown into one dialog
    # (2026-05-22). Width 80, vertical grouping.
    ("SettingsModal",              lambda: sc.SettingsModal()),
    ("EditGrammarConfirmModal",    lambda: sc.EditGrammarConfirmModal(
                                                "MyCustomGrammar",
                                                n_dependents=3)),
    ("_ConfirmDeleteGrammarModal", lambda: sc._ConfirmDeleteGrammarModal(
                                                "MyCustomGrammar",
                                                "  Test body markup line.\n"
                                                "  3 dependents.")),
    # SimulatorScreen is a full-screen workbench (like PrimerDesignScreen /
    # FeatureLibraryScreen). Tested with a small linear template + minimal
    # feats so compose has no I/O.
    ("SimulatorScreen",            lambda: sc.SimulatorScreen(
                                                "ATGC" * 100,
                                                _make_feats(),
                                                "test",
                                                "circular")),
    # Experiments lab-notebook modals (added 2026-05-18 with the
    # Experiments toolbar feature). ExperimentsScreen itself is a
    # full-screen Screen (like SimulatorScreen / SequencingScreen);
    # it's covered by its own tests in `test_experiments.py`. These
    # are the modal sidekicks that DO need the 160×48 fit guarantee.
    ("ExperimentDeleteConfirmModal", lambda: sc.ExperimentDeleteConfirmModal(
                                                  title="Delete: foo?",
                                                  body="body text")),
    ("ExperimentRenameModal",        lambda: sc.ExperimentRenameModal(
                                                  "Old title")),
    ("ImageAttachModal",             lambda: sc.ImageAttachModal()),
    ("SpellcheckModal.empty",        lambda: sc.SpellcheckModal([])),
    ("SpellcheckModal.populated",    lambda: sc.SpellcheckModal([
                                                  ("teste",
                                                   ["test", "tester"]),
                                                  ("mispelt", ["misspelt"]),
                                              ])),
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


class TestButtonsFitInsideDialogContent:
    """Regression guard 2026-05-23: a modal whose button-bar container
    uses `padding-top: 1` on a fixed `height: 3` only leaves 2 rows for
    the actual buttons, which need 3 (top border + label + bottom
    border). The bottom border bled into the dialog's padding-bottom
    and rendered as cropped.

    The fix swapped to `margin-top: 1` so the gap is allocated OUTSIDE
    the fixed-height container, preserving the 3-row budget for the
    button widget. This test pins both fixed modals — register every
    Button on the modal and assert its bottom row sits within the
    dialog's content area (excludes the dialog's own padding+border).
    """

    @pytest.mark.parametrize("label,factory", [
        ("AlignmentManagerModal",
         lambda: sc.AlignmentManagerModal([])),
        ("MultiAlignPickerModal",
         lambda: sc.MultiAlignPickerModal()),
    ])
    async def test_button_bottom_fits_dialog_content(
            self, label, factory, tiny_record, isolated_library):
        from textual.widgets import Button
        from textual.containers import Vertical
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = factory()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            # The modal's dialog is the only top-level Vertical with an
            # id (both modals follow the `#…-dlg` convention).
            dialog = next(
                (v for v in modal.query(Vertical)
                 if v.id and v.id.endswith("-dlg")),
                None,
            )
            assert dialog is not None, (
                f"{label}: could not locate dialog Vertical container"
            )
            # Dialog content area = dialog region minus border (1 each
            # side) and padding (1 each side via `padding: 1 2`).
            dlg_r = dialog.region
            content_bottom = dlg_r.y + dlg_r.height - 1 - 2  # border-bottom + padding-bottom
            for btn in modal.query(Button):
                br = btn.region
                btn_bottom = br.y + br.height - 1
                assert btn_bottom <= content_bottom, (
                    f"{label}: Button#{btn.id!r} bottom y={btn_bottom} "
                    f"overflows dialog content area (ends at y="
                    f"{content_bottom}). The button-bar container "
                    f"likely uses padding-top inside fixed height; "
                    f"swap to margin-top so the buttons keep their "
                    f"full row budget."
                )


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
            app.push_screen(sc.AddFeatureModal(selection_range=(0, 10)))
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
            app.push_screen(sc.AddFeatureModal(selection_range=(0, 10)))
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



# ═══════════════════════════════════════════════════════════════════════════════
# Centering invariant — sweep #13 (2026-05-20 UX audit)
# ═══════════════════════════════════════════════════════════════════════════════

class TestModalCenteringInvariant:
    """Every ModalScreen subclass — except DropdownScreen (which is
    a positioned popup anchored to a menubar item via explicit
    `styles.offset`) — must sit at the centre of the terminal
    window. Backed by the global `ModalScreen { align: center
    middle; }` rule in PlasmidApp.CSS that the 2026-05-20 audit
    added so future modals get centring for free.
    """

    def _all_modal_subclasses(self) -> list[type]:
        """Walk sc.* for every class derived from textual.screen.ModalScreen
        (excluding ModalScreen itself, abstract bases, and DropdownScreen).
        """
        from textual.screen import ModalScreen
        out: list[type] = []
        for name in dir(sc):
            obj = getattr(sc, name, None)
            if (isinstance(obj, type)
                    and issubclass(obj, ModalScreen)
                    and obj is not ModalScreen
                    and name != "DropdownScreen"):
                out.append(obj)
        return out

    def test_global_modalscreen_centering_rule_present(self):
        """The PlasmidApp.CSS must carry the global centering rule
        so every ModalScreen subclass inherits centre alignment
        without needing to repeat the rule per-class."""
        assert "ModalScreen { align: center middle; }" in sc.PlasmidApp.CSS

    def test_dropdownscreen_override_present(self):
        """DropdownScreen's explicit `styles.offset` positioning
        needs the global ModalScreen centring undone — verify the
        override is in PlasmidApp.CSS so the menubar dropdown
        keeps landing under its triggering menu item."""
        assert "DropdownScreen { align: left top; }" in sc.PlasmidApp.CSS

    def test_audit_subclass_coverage_substantial(self):
        """Sanity check the test harness — must discover the bulk
        of the modal subclasses, not a handful. Catches an import
        regression that would silently empty the audit set."""
        klasses = self._all_modal_subclasses()
        assert len(klasses) >= 30, (
            f"only {len(klasses)} modal subclasses discovered — "
            "expected 30+, possible import regression"
        )

    async def test_dropdownscreen_still_top_left_aligned(self):
        """The DropdownScreen override must keep its natural top-
        left layout so `styles.offset = (x, y)` continues to anchor
        the popup correctly to its menubar item."""
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            await pilot.pause()
            modal = sc.DropdownScreen(
                items=[("Item A", "noop"), ("Item B", "noop")],
                x=5, y=2,
            )
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause()
            scr = app.screen
            assert str(scr.styles.align_horizontal) == "left"
            assert str(scr.styles.align_vertical) == "top"
            modal.dismiss(None)
            await pilot.pause()
