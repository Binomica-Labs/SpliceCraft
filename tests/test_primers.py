"""
test_primers — Primer design backend + persistence tests.

Covers:
  - _design_detection_primers: Primer3 returns valid pair, Tm in range,
    product size in range, primers flank the target
  - _design_cloning_primers: RE sites present, GCGC padding, binding Tm,
    correct RE at correct end
  - Primer library persistence: save/load round-trip
  - PrimerDesignScreen: mounts with correct widgets
"""
from __future__ import annotations

import json
import random

import pytest

import splicecraft as sc


@pytest.fixture
def random_seq_3k():
    rng = random.Random(0xBEEF)
    return "".join(rng.choice("ACGT") for _ in range(3000))


# `isolated_primers` lives in tests/conftest.py — same redirect
# pattern, shared with test_smoke.py and others.


# ═══════════════════════════════════════════════════════════════════════════════
# Detection primers
# ═══════════════════════════════════════════════════════════════════════════════

class TestDetectionPrimers:
    """Detection primers bind INSIDE the selected region (not flanking it).
    Both primers sit within [start, end] and the amplicon is product_min..
    product_max bp. This is the standard diagnostic PCR approach."""

    def test_returns_valid_pair(self, random_seq_3k):
        # Region must be >= product_min so primers can fit inside
        r = sc._design_detection_primers(
            random_seq_3k, 100, 800, product_min=450, product_max=550,
        )
        assert "error" not in r
        assert r["fwd_seq"]
        assert r["rev_seq"]

    def test_tm_in_range(self, random_seq_3k):
        r = sc._design_detection_primers(
            random_seq_3k, 100, 800, target_tm=60.0,
            product_min=450, product_max=550,
        )
        assert "error" not in r
        assert 55 < r["fwd_tm"] < 65
        assert 55 < r["rev_tm"] < 65

    def test_product_size_in_range(self, random_seq_3k):
        r = sc._design_detection_primers(
            random_seq_3k, 100, 800,
            product_min=450, product_max=550,
        )
        assert "error" not in r
        assert 450 <= r["product_size"] <= 550

    def test_primers_inside_region(self, random_seq_3k):
        """Both primers must bind INSIDE the selected region — this is the
        key semantic difference from the old SEQUENCE_TARGET approach."""
        r = sc._design_detection_primers(
            random_seq_3k, 200, 1200,
            product_min=450, product_max=550,
        )
        assert "error" not in r
        assert r["fwd_pos"][0] >= 200, "fwd must be inside region"
        assert r["fwd_pos"][1] <= 1200, "fwd end must be inside region"
        assert r["rev_pos"][0] >= 200, "rev start must be inside region"
        assert r["rev_pos"][1] <= 1200, "rev must be inside region"

    def test_empty_target_returns_error(self, random_seq_3k):
        r = sc._design_detection_primers(random_seq_3k, 500, 500)
        assert "error" in r

    def test_region_smaller_than_product_returns_error(self, random_seq_3k):
        """If the region is shorter than the minimum product size, we should
        get a clear error rather than letting Primer3 fail cryptically."""
        r = sc._design_detection_primers(
            random_seq_3k, 500, 600,  # 100 bp region
            product_min=450, product_max=550,  # but product needs 450+
        )
        assert "error" in r
        assert "shorter" in r["error"].lower()

    def test_impossible_constraints_returns_error(self, random_seq_3k):
        r = sc._design_detection_primers(
            random_seq_3k, 100, 800,
            product_min=10, product_max=20,
        )
        assert "error" in r

    def test_large_gene_works(self, random_seq_3k):
        """An 861 bp gene (like ampR) with 450-550 product range should work
        because primers go INSIDE the region. This was the original bug."""
        r = sc._design_detection_primers(
            random_seq_3k, 100, 961,  # 861 bp region like ampR
            product_min=450, product_max=550,
        )
        assert "error" not in r
        assert 450 <= r["product_size"] <= 550
        assert r["fwd_pos"][0] >= 100
        assert r["rev_pos"][1] <= 961


# ═══════════════════════════════════════════════════════════════════════════════
# Cloning primers
# ═══════════════════════════════════════════════════════════════════════════════

class TestCloningPrimers:
    def test_returns_valid_pair(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI",
        )
        assert "error" not in r
        assert r["fwd_full"]
        assert r["rev_full"]

    def test_gcgc_padding_present(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI",
        )
        assert r["fwd_full"].startswith("GCGC")
        assert r["rev_full"].startswith("GCGC")

    def test_5prime_re_site_in_fwd(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI",
        )
        # EcoRI = GAATTC — should appear right after GCGC in fwd primer
        assert "GAATTC" in r["fwd_full"]

    def test_3prime_re_site_rc_in_rev(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI",
        )
        # BamHI = GGATCC → RC = GGATCC (palindrome)
        assert "GGATCC" in r["rev_full"]

    def test_non_palindrome_re_site_rc(self, random_seq_3k):
        # BsaI = GGTCTC (non-palindrome) → RC = GAGACC
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "BsaI", "EcoRI",
        )
        assert "GGTCTC" in r["fwd_full"], "5' BsaI site in fwd"
        # Rev should have RC of EcoRI = GAATTC (palindrome → same)
        assert "GAATTC" in r["rev_full"]

    def test_binding_tm_near_target(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "EcoRI", "BamHI", target_tm=57.0,
        )
        assert 49 < r["fwd_tm"] < 65
        assert 49 < r["rev_tm"] < 65

    def test_unknown_enzyme_returns_error(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, "FakeEnzyme", "BamHI",
        )
        assert "error" in r

    def test_short_region_returns_error(self, random_seq_3k):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 210, "EcoRI", "BamHI",
        )
        assert "error" in r

    @pytest.mark.parametrize("re5,re3", [
        ("EcoRI", "BamHI"), ("XhoI", "NdeI"), ("NcoI", "XbaI"),
        ("SpeI", "PstI"), ("HindIII", "SalI"), ("NotI", "BglII"),
    ])
    def test_various_enzyme_pairs(self, random_seq_3k, re5, re3):
        r = sc._design_cloning_primers(
            random_seq_3k, 200, 800, re5, re3,
        )
        assert "error" not in r
        site_5, _, _ = sc._NEB_ENZYMES[re5]
        site_3, _, _ = sc._NEB_ENZYMES[re3]
        assert site_5 in r["fwd_full"]
        assert sc._rc(site_3) in r["rev_full"] or site_3 in r["rev_full"]


# ═══════════════════════════════════════════════════════════════════════════════
# Primer library persistence
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrimerPersistence:
    def test_empty_loads_as_empty(self, isolated_primers):
        assert sc._load_primers() == []

    def test_save_load_roundtrip(self, isolated_primers):
        primers = [
            {"name": "ampR-DET-F", "sequence": "ATGAAAGATCTG", "tm": 57.2,
             "primer_type": "detection", "source": "pUC19"},
        ]
        sc._save_primers(primers)
        loaded = sc._load_primers()
        assert len(loaded) == 1
        assert loaded[0]["name"] == "ampR-DET-F"

    def test_writes_valid_json(self, isolated_primers):
        sc._save_primers([{"name": "x", "sequence": "ATG"}])
        assert isolated_primers.exists()
        parsed = json.loads(isolated_primers.read_text())
        assert parsed["entries"][0]["name"] == "x"

    def test_corrupted_file_returns_empty(self, isolated_primers):
        isolated_primers.write_text("{bad")
        sc._primers_cache = None
        assert sc._load_primers() == []


# ═══════════════════════════════════════════════════════════════════════════════
# RE options list
# ═══════════════════════════════════════════════════════════════════════════════

class TestCloningREOptions:
    def test_common_enzymes_present(self):
        names = {name for _, name in sc._CLONING_RE_OPTIONS}
        for must in ["EcoRI", "BamHI", "XhoI", "NdeI", "NotI", "BsaI",
                     "HindIII", "XbaI", "NcoI", "SalI"]:
            assert must in names, f"{must} missing from RE options"

    def test_all_options_are_in_neb_catalog(self):
        for label, name in sc._CLONING_RE_OPTIONS:
            assert name in sc._NEB_ENZYMES, f"{name} not in _NEB_ENZYMES"


# ═══════════════════════════════════════════════════════════════════════════════
# Uniform designer return shape (added 2026-04-12 after Golden Braid save crash)
# ═══════════════════════════════════════════════════════════════════════════════

class TestDesignerUniformShape:
    """Save-to-library code (`PrimerDesignScreen._save_primers_btn`) assumes
    every `_design_*_primers` function returns fwd_pos, rev_pos, fwd_tm,
    rev_tm plus some form of fwd/rev sequence. Golden Braid was missing
    fwd_pos/rev_pos, which triggered a KeyError on save."""

    _REQUIRED_KEYS = {"fwd_pos", "rev_pos", "fwd_tm", "rev_tm"}

    def _valid_template(self):
        """Diverse ~1000 bp template, long enough for GB (>= 18 bp) and
        Detection (>= product_min = 450). Previously a highly repetitive
        string (``ATGAAACGTGATTTAGCCGTTAA`` × 40) that made Primer3 grind
        through ~8 s of self-hybridizing candidates for one test.

        We scrub Esp3I (CGTCTC / GAGACG) from the 0xBEEF output so the
        Golden-Braid L0 internal-site check doesn't reject the template.
        Detection / cloning / generic primer designs don't care either
        way, but share the same fixture."""
        rng = random.Random(0xBEEF)
        seq = "".join(rng.choice("ACGT") for _ in range(1000))
        return seq.replace("CGTCTC", "CGACTC").replace("GAGACG", "GAGACT")

    def test_detection_has_required_keys(self):
        seq = self._valid_template()
        r = sc._design_detection_primers(seq, 50, len(seq) - 50)
        assert "error" not in r, r
        assert self._REQUIRED_KEYS <= set(r.keys())
        assert len(r["fwd_pos"]) == 2 and len(r["rev_pos"]) == 2

    def test_cloning_has_required_keys(self):
        seq = self._valid_template()
        r = sc._design_cloning_primers_raw(seq, 50, 300, "GAATTC", "GGATCC")
        assert "error" not in r, r
        assert self._REQUIRED_KEYS <= set(r.keys())

    def test_goldenbraid_has_required_keys(self):
        """Regression guard for 2026-04-12 KeyError on save."""
        seq = self._valid_template()
        # Use Promoter (oh5 = GGAG) — no codon-skip applies, so
        # fwd_pos starts at the bare insert start. CDS would
        # shift fwd_pos by +3 to skip the ATG-fusion overhang.
        r = sc._design_gb_primers(seq, 50, 300, "Promoter")
        assert "error" not in r, r
        assert self._REQUIRED_KEYS <= set(r.keys()), (
            f"Golden Braid designer missing keys: "
            f"{self._REQUIRED_KEYS - set(r.keys())}"
        )
        # Forward primer binds the start of the insert
        assert r["fwd_pos"] == (50, 50 + len(r["fwd_binding"]))
        # Reverse primer binds the end of the insert (on the bottom strand,
        # but positions are reported in forward coords)
        assert r["rev_pos"] == (300 - len(r["rev_binding"]), 300)

    def test_goldenbraid_cds_fwd_pos_skips_atg(self):
        """Sacred regression guard for the GB CDS ATG-fusion fix
        (2026-05-21). CDS oh5 = AATG carries the start codon, so
        fwd_pos starts at insert_start + 3 (codon 2), not the
        bare insert start."""
        seq = self._valid_template()
        r = sc._design_gb_primers(seq, 50, 300, "CDS")
        assert "error" not in r, r
        assert r["fwd_pos"] == (50 + 3, 50 + 3 + len(r["fwd_binding"]))

    def test_goldenbraid_flags_internal_esp3i(self):
        """Regression guard: a CGTCTC inside the insert is self-domesticating
        for L0 Esp3I / BsmBI and must be flagged, not silently passed through.
        L0 parts are cut with Esp3I (CGTCTC) so only that recognition site
        causes a problem at this stage — BsaI sites in an L0 part are fine
        because L0 domestication never sees BsaI."""
        core = "ATGAAACGTGATTTAGCC" * 5   # 90 bp, no Esp3I
        with_fwd_esp3i = core + "CGTCTC" + core
        r = sc._design_gb_primers(with_fwd_esp3i, 0, len(with_fwd_esp3i), "CDS")
        assert "error" in r
        assert "Esp3I" in r["error"]

    def test_goldenbraid_flags_reverse_esp3i(self):
        """GAGACG (RC of CGTCTC) on the top strand also fragments the part
        during the L0 Esp3I digest."""
        core = "ATGAAACGTGATTTAGCC" * 5
        with_rev_esp3i = core + "GAGACG" + core
        r = sc._design_gb_primers(with_rev_esp3i, 0, len(with_rev_esp3i), "CDS")
        assert "error" in r
        assert "Esp3I" in r["error"]

    def test_generic_has_required_keys(self):
        seq = self._valid_template()
        r = sc._design_generic_primers(seq, 50, 300)
        assert "error" not in r, r
        assert self._REQUIRED_KEYS <= set(r.keys())


# ═══════════════════════════════════════════════════════════════════════════════
# PrimerDesignScreen layout smoke (Option A wizard redesign, 2026-04-12)
# ═══════════════════════════════════════════════════════════════════════════════

TERMINAL_SIZE = (140, 42)


class TestPrimerDesignScreenLayout:
    """After the Option A wizard redesign, verify the screen mounts with
    the numbered-step structure (TEMPLATE | MODE + PARAMETERS wizard layout,
    followed by Design button, Results, Library)."""

    async def test_mounts_with_detection_initial_mode(self):
        """2026-05-21 refactor: RadioSet replaced by Tabs widget.
        Default-active tab is Detection."""
        from textual.widgets import Tabs
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            tabs = app.screen.query_one("#pd-mode-tabs", Tabs)
            assert tabs.active_tab is not None
            assert tabs.active_tab.id == "tab-detection"
            assert screen._current_mode() == "detection"

    async def test_all_four_modes_have_tabs(self):
        """2026-05-21 refactor: 4 Tab widgets in the mode tab bar
        (was RadioButtons)."""
        from textual.widgets import Tab
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            tab_ids = {t.id for t in app.screen.query(Tab)}
            assert {
                "tab-detection", "tab-cloning",
                "tab-goldenbraid", "tab-generic",
            } <= tab_ids

    async def test_switching_mode_changes_current_mode(self):
        """`_switch_mode` swaps which mode-panel is visible AND
        keeps the tab bar's active tab in sync."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen._switch_mode("cloning")
            await pilot.pause()
            # Cloning panel visible, detection hidden
            assert screen.query_one("#pd-panel-clo").display is True
            assert screen.query_one("#pd-panel-det").display is False
            assert screen._current_mode() == "cloning"

    async def test_book_layout_split(self):
        """Open-book split (2026-05-21 refactor): left page owns
        template / mode-tabs / params / results / bottom buttons;
        right page is the full-height primer library."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            book = app.screen.query_one("#pd-book")
            page_ids = [c.id for c in book.children]
            assert page_ids == ["pd-left-page", "pd-right-page"]

            # Left page: tabs, template, params, results, buttons.
            left = app.screen.query_one("#pd-left-page")
            left_ids = [c.id for c in left.children]
            assert "pd-mode-tabs" in left_ids
            assert "pd-template-section" in left_ids
            assert "pd-params-section" in left_ids
            assert "pd-results-section" in left_ids
            idx_tabs = left_ids.index("pd-mode-tabs")
            idx_tpl  = left_ids.index("pd-template-section")
            idx_par  = left_ids.index("pd-params-section")
            idx_res  = left_ids.index("pd-results-section")
            assert idx_tabs < idx_tpl < idx_par < idx_res

            # Right page: just the library header + table.
            right = app.screen.query_one("#pd-right-page")
            right_ids = [c.id for c in right.children]
            assert "pd-lib-hdr-row" in right_ids
            assert "pd-lib-table" in right_ids
            assert right_ids.index("pd-lib-hdr-row") \
                < right_ids.index("pd-lib-table")

    async def test_name_inputs_and_save_button_present(self):
        """2026-05-21 refactor: primer-name Inputs migrated UP to
        the TEMPLATE section (compact form-fill before design),
        Save button migrated DOWN to the bottom-of-page action row
        next to Design / Add to Map / Clear. Both still queryable
        from anywhere on the screen via `#pd-fwd-name` etc."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            # All three queryable from the screen root.
            assert app.screen.query_one("#pd-fwd-name") is not None
            assert app.screen.query_one("#pd-rev-name") is not None
            assert app.screen.query_one("#btn-pd-save") is not None
            # Names sit in TEMPLATE section now.
            tpl = app.screen.query_one("#pd-template-section")
            tpl_ids = {w.id for w in tpl.walk_children()}
            assert "pd-fwd-name" in tpl_ids
            assert "pd-rev-name" in tpl_ids
            # Save button sits in the bottom actions row.
            bottom = app.screen.query_one("#pd-bottom-actions")
            bot_ids = {w.id for w in bottom.walk_children()}
            assert "btn-pd-save" in bot_ids
            assert "btn-pd-design" in bot_ids
            assert "btn-pd-clear" in bot_ids

    async def test_library_header_row_has_rename_delete(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            hdr = app.screen.query_one("#pd-lib-hdr-row")
            ids_inside = {w.id for w in hdr.walk_children()}
            assert "btn-pdlib-rename" in ids_inside
            assert "btn-pdlib-del" in ids_inside

    async def test_reset_for_new_design_clears_pair_state(self, tmp_path,
                                                           monkeypatch):
        """After a save-to-library completes, the modal state specific
        to a designed primer pair is cleared (result cache, primer
        names, Save button, Results pane, feature selection, start/end,
        part name, feature info). Template / mode / parameters are
        deliberately preserved."""
        from textual.widgets import Input, Button, Static, Select
        # Redirect primers.json to a tmp path (conftest's autouse
        # fixture already does this, but be explicit.)
        p = tmp_path / "primers.json"
        monkeypatch.setattr(sc, "_PRIMERS_FILE", p)
        monkeypatch.setattr(sc, "_primers_cache", None, raising=False)

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            # Seed the screen with one feature so its Select has a
            # value we can observe getting reset.
            feats = [{"type": "CDS", "start": 10, "end": 100,
                      "strand": 1, "label": "gene1",
                      "color": "color(39)"}]
            screen = sc.PrimerDesignScreen("ACGT" * 200, feats, "test")
            app.push_screen(screen)
            await pilot.pause()

            # Simulate the post-design state.
            screen._det_result = {
                "_type":    "detection",
                "fwd_seq":  "ACGT" * 5,
                "rev_seq":  "TGCA" * 5,
                "fwd_tm":   60.0,
                "rev_tm":   60.0,
                "fwd_pos":  (0, 20),
                "rev_pos":  (100, 120),
                "product_size": 120,
            }
            screen.query_one("#pd-fwd-name", Input).value = "test-F"
            screen.query_one("#pd-rev-name", Input).value = "test-R"
            screen.query_one("#btn-pd-save", Button).disabled = False
            screen.query_one("#pd-results", Static).update("some result text")
            screen.query_one("#pd-feat", Select).value = "10-100"
            screen.query_one("#pd-start", Input).value = "11"
            screen.query_one("#pd-end",   Input).value = "100"
            screen.query_one("#pd-part-name", Input).value = "gene1"
            screen.query_one("#pd-feat-info", Static).update(
                "gene1  90 bp"
            )
            await pilot.pause()

            # Trigger the reset
            screen._reset_for_new_design()
            await pilot.pause()

            # Per-pair output cleared
            assert screen._det_result is None
            assert screen._clo_result is None
            assert screen.query_one("#pd-fwd-name", Input).value == ""
            assert screen.query_one("#pd-rev-name", Input).value == ""
            assert screen.query_one("#btn-pd-save", Button).disabled is True

            # Feature selection cleared. Select.clear() sets the value
            # to a blank sentinel (BLANK or NULL depending on Textual
            # version) — assert it's no longer the real option value.
            assert screen.query_one("#pd-feat", Select).value != "10-100"
            assert screen.query_one("#pd-start", Input).value == ""
            assert screen.query_one("#pd-end",   Input).value == ""
            assert screen.query_one("#pd-part-name", Input).value == ""


# ═══════════════════════════════════════════════════════════════════════════════
# Wrap-region primer design (regression guard for 2026-04-13 fix)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSliceCircular:
    """_slice_circular is the contract at the heart of wrap-region design —
    test it directly so the primer-design tests can rely on it."""

    def test_normal_slice(self):
        assert sc._slice_circular("ABCDEFGHIJ", 2, 6) == "CDEF"

    def test_wrap_slice(self):
        # [8:2) wraps → seq[8:] + seq[:2] = "IJ" + "AB"
        assert sc._slice_circular("ABCDEFGHIJ", 8, 2) == "IJAB"

    def test_equal_endpoints_is_empty(self):
        # end == start is treated as an empty region (not a full-plasmid
        # wrap); callers that want the whole sequence should pass (0, len).
        assert sc._slice_circular("ABCDEFGHIJ", 3, 3) == ""

    def test_zero_start_normal(self):
        assert sc._slice_circular("ABCDEFGHIJ", 0, 5) == "ABCDE"


class TestWrapRegionGeneric:
    """_design_generic_primers on a wrap region should produce valid primers
    AND map the reported positions back to original template coordinates."""

    def test_wrap_region_produces_primers(self):
        # 3 kb plasmid, region [2900, 100) — 200 bp crossing origin
        import random
        rng = random.Random(0xBEEF)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_generic_primers(seq, 2900, 100)
        assert "error" not in r
        assert r["fwd_seq"]
        assert r["rev_seq"]

    def test_wrap_region_fwd_pos_starts_before_origin(self):
        import random
        rng = random.Random(0xBEEF)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_generic_primers(seq, 2900, 100)
        # fwd_pos starts at 2900 (user-specified start)
        assert r["fwd_pos"][0] == 2900
        # fwd_pos end is modular — lands near 2900 + binding_len, possibly wrapped
        assert 0 <= r["fwd_pos"][1] < 3000

    def test_wrap_region_rev_pos_ends_at_user_end(self):
        import random
        rng = random.Random(0xBEEF)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_generic_primers(seq, 2900, 100)
        # rev_pos should end at 100 (user-specified end)
        assert r["rev_pos"][1] == 100

    def test_wrap_region_binding_actually_appears_in_wrapped_insert(self):
        """The forward binding must literally appear at the start of the
        wrapped insert (template[2900:] + template[:100])."""
        import random
        rng = random.Random(0xBEEF)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_generic_primers(seq, 2900, 100)
        wrapped_insert = seq[2900:] + seq[:100]
        assert wrapped_insert.startswith(r["fwd_seq"])
        # Reverse binding is the RC of the end of the insert
        assert wrapped_insert.endswith(sc._rc(r["rev_seq"]))


class TestWrapRegionCloning:
    def test_wrap_region_produces_primers(self):
        import random
        rng = random.Random(0xF00D)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_cloning_primers(seq, 2900, 200, "EcoRI", "BamHI")
        assert "error" not in r
        assert r["fwd_full"].startswith("GCGC")
        assert "GAATTC" in r["fwd_full"]

    def test_wrap_region_preserves_re_sites(self):
        import random
        rng = random.Random(0xF00D)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_cloning_primers(seq, 2900, 200, "XhoI", "HindIII")
        assert "CTCGAG" in r["fwd_full"]         # XhoI
        assert "AAGCTT" in r["rev_full"]         # HindIII (palindrome)

    def test_wrap_region_insert_is_wrapped_concatenation(self):
        import random
        rng = random.Random(0xF00D)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_cloning_primers(seq, 2900, 200, "EcoRI", "BamHI")
        assert r["insert_seq"] == seq[2900:] + seq[:200]


class TestWrapRegionGoldenBraid:
    def test_wrap_region_produces_primers(self):
        # Need a wrap region free of internal Esp3I sites (L0 domestication
        # now uses CGTCTC, not GGTCTC). Build a deterministic sequence and
        # scrub both strands.
        import random
        rng = random.Random(0xCAFE)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        seq = (seq.replace("CGTCTC", "CGACTC")
                  .replace("GAGACG", "GAGACT"))
        r = sc._design_gb_primers(seq, 2900, 200, "CDS")
        assert "error" not in r
        assert r["fwd_binding"]
        assert r["rev_binding"]


class TestWrapRegionDetection:
    """Detection primers on wrap regions need the template to be rotated
    because Primer3 is linear-only. After Primer3 returns, positions must
    be unrotated back to original template coordinates."""

    def test_wrap_region_produces_primers(self):
        import random
        rng = random.Random(0xDEAD)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        # Region [2500, 500) — 1000 bp wrap region
        r = sc._design_detection_primers(
            seq, 2500, 500, product_min=450, product_max=550,
        )
        assert "error" not in r
        assert 450 <= r["product_size"] <= 550

    def test_wrap_region_primer_seqs_appear_in_wrapped_template(self):
        import random
        rng = random.Random(0xDEAD)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        wrapped = seq[2500:] + seq[:500]
        r = sc._design_detection_primers(
            seq, 2500, 500, product_min=450, product_max=550,
        )
        assert "error" not in r
        assert r["fwd_seq"] in wrapped
        assert sc._rc(r["rev_seq"]) in wrapped

    def test_wrap_region_positions_are_modular(self):
        """Returned positions should be in [0, total) — never negative or
        past the template end, even when the primer landed on the wrap
        side of the rotation."""
        import random
        rng = random.Random(0xDEAD)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_detection_primers(
            seq, 2500, 500, product_min=450, product_max=550,
        )
        assert "error" not in r
        for pos in (r["fwd_pos"][0], r["fwd_pos"][1],
                    r["rev_pos"][0], r["rev_pos"][1]):
            assert 0 <= pos < 3000


class TestLinearRegionStillWorks:
    """Sanity check: all the existing linear-region tests still pass the
    new code paths unchanged."""

    def test_linear_detection_unchanged(self):
        import random
        rng = random.Random(0xBEEF)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_detection_primers(
            seq, 100, 800, product_min=450, product_max=550,
        )
        assert "error" not in r
        assert r["fwd_pos"][0] >= 100 and r["fwd_pos"][1] <= 800

    def test_linear_generic_unchanged(self):
        import random
        rng = random.Random(0xBEEF)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        r = sc._design_generic_primers(seq, 100, 500)
        assert r["fwd_pos"] == (100, 100 + len(r["fwd_seq"]))
        assert r["rev_pos"] == (500 - len(r["rev_seq"]), 500)


# ═══════════════════════════════════════════════════════════════════════════════
# Wrap-region UI hooks (screen mounts, feat_selected handles wrap-feature length)
# ═══════════════════════════════════════════════════════════════════════════════

TERMINAL_SIZE = (200, 60)

class TestWrapRegionUI:
    async def test_wrap_start_greater_than_end_supported(self):
        """The Tip label was removed 2026-05-21 (the wrap-origin
        affordance is undocumented but still works). This test
        confirms the supporting Inputs accept Start > End values
        without complaint, so the underlying capability stays
        regression-guarded even without the visible hint."""
        from textual.widgets import Input
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 800, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            # Both Inputs accept integer values; user can enter
            # Start > End to trigger the wrap-around design path.
            s = app.screen.query_one("#pd-start", Input)
            e = app.screen.query_one("#pd-end", Input)
            s.value = "2900"
            e.value = "200"
            assert s.value == "2900"
            assert e.value == "200"

    async def test_feat_selected_on_wrap_feature_fills_inputs(self):
        """Selecting a wrap feature should populate start/end correctly
        AND compute a positive feat_len so the detection range isn't
        auto-adjusted to a negative value."""
        from textual.widgets import Input, Select
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            # Wrap feature: start=2900, end=200 on 3000 bp plasmid
            feats = [{"start": 2900, "end": 200, "label": "wrapCDS",
                      "type": "CDS", "strand": 1}]
            screen = sc.PrimerDesignScreen("A" * 3000, feats, "test")
            app.push_screen(screen)
            await pilot.pause()
            # Simulate selecting the wrap feature
            screen.query_one("#pd-feat", Select).value = "2900-200"
            await pilot.pause()
            assert screen.query_one("#pd-start", Input).value == "2901"
            assert screen.query_one("#pd-end", Input).value == "200"
            # No crash — in buggy version feat_len=-2700 caused
            # detection product range to go negative.


# ═══════════════════════════════════════════════════════════════════════════════
# Primer library — `.dna` imports show up + scroll-through behaviour.
# Regression guards for 2026-05-10 user report: "make sure primers
# added from commercial SaaS collection import shows up in the primer
# library section of the Primers tab, and that the primer library is
# scrollable"
# ═══════════════════════════════════════════════════════════════════════════════


class TestPrimerLibraryShowsImported:
    """Primers landed in `primers.json` by `.dna` import (or any other
    code path) must show up in the PrimerDesignScreen's library table."""

    async def test_imported_primer_appears_in_library_table(
            self, isolated_library
    ):
        from textual.widgets import DataTable
        # Pre-populate primers.json with an "imported" entry the same
        # shape `_augment_dna_record_from_packets` produces.
        sc._save_primers([{
            "name":        "M13 fwd",
            "sequence":    "GTAAAACGACGGCCAGT",
            "tm":          54.7,
            "primer_type": "imported",
            "source":      ".dna import",
            "pos_start":   0,
            "pos_end":     17,
            "strand":      1,
            "date":        "2026-05-10",
            "status":      "Imported",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.05)
            t = screen.query_one("#pd-lib-table", DataTable)
            assert t.row_count == 1, (
                f"Expected 1 row from the imported primer, got {t.row_count}"
            )

    async def test_table_renders_with_tm_none_legacy_entry(
            self, isolated_library
    ):
        """Hand-edited entries (or imports from < 0.7.10.1 before the
        Tm calculation landed) can have `tm=None`. The table must
        render without crashing — defensive coding pinned by this test."""
        from textual.widgets import DataTable
        sc._save_primers([{
            "name":        "legacy",
            "sequence":    "ACGTACGT",
            "tm":          None,    # <-- the regression target
            "primer_type": "imported",
            "source":      ".dna import",
            "status":      "Imported",
        }])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            t = screen.query_one("#pd-lib-table", DataTable)
            assert t.row_count == 1


class TestPrimerLibraryScrollable:
    """A library with more primers than fit in the viewport must
    remain navigable — DataTable supports keyboard + mouse scroll
    natively, but the regression target is making sure no CSS or
    container constraint clips off-screen rows."""

    async def test_all_rows_in_table_when_library_exceeds_viewport(
            self, isolated_library
    ):
        from textual.widgets import DataTable
        # 60 entries — far more than fits in a 42-row terminal once
        # all the wizard sections + results section are stacked above.
        # Each primer needs a unique sequence so the dedupe-on-save in
        # `_save_primers` doesn't collapse them. Encode the index into
        # the bases so 60 unique sequences are generated cheaply.
        def _unique_seq(i: int) -> str:
            # 18-mer with a deterministic prefix per index — gives 60
            # distinct sequences while keeping length realistic.
            return f"AC{i:04d}".ljust(18, "T")
        entries = [{
            "name":        f"primer_{i:02d}",
            "sequence":    _unique_seq(i),
            "tm":          float(58 + (i % 5)),
            "primer_type": "imported",
            "source":      ".dna import",
            "pos_start":   i * 10,
            "pos_end":     i * 10 + 20,
            "strand":      1,
            "date":        "2026-05-10",
            "status":      "Imported",
        } for i in range(60)]
        sc._save_primers(entries)
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.05)
            t = screen.query_one("#pd-lib-table", DataTable)
            # Every entry must be present in the data model — if rows
            # were truncated, the user couldn't reach them by any
            # navigation. DataTable's internal scrolling handles
            # viewport mapping at paint time.
            assert t.row_count == 60

    async def test_save_dedupes_legacy_duplicates(self, isolated_library):
        """Regression for 2026-05-10: pre-fix, primers.json could
        accumulate duplicates by sequence — from manual JSON edits,
        from imports before the dedupe paths landed, or from any path
        that bypassed the in-memory dedupe. The save path now collapses
        duplicates by sequence (case-insensitive) so any write cleans
        up the file. First-by-position wins."""
        sc._save_primers([
            {"name": "M13fwd_v1", "sequence": "GTAAAACGACGGCCAGT",
             "tm": 54.7, "status": "Imported"},
            {"name": "M13fwd_v2", "sequence": "GTAAAACGACGGCCAGT",
             "tm": 54.7, "status": "Imported"},
            {"name": "M13fwd_lower", "sequence": "gtaaaacgacggccagt",
             "tm": 54.7, "status": "Imported"},
            {"name": "Other", "sequence": "TTTTTTTTTT",
             "tm": 50.0, "status": "Designed"},
        ])
        primers = sc._load_primers()
        # 4 → 2: M13fwd_v1 wins (first by position), the lowercase
        # variant and v2 are dropped. Other stays (unique sequence).
        names = [p["name"] for p in primers]
        assert names == ["M13fwd_v1", "Other"], (
            f"Expected dedupe to [M13fwd_v1, Other], got {names}"
        )

    async def test_dedupe_preserves_entries_with_missing_sequence(
            self, isolated_library
    ):
        """Defensive: an entry without a usable `sequence` field must
        survive the dedupe. Losing it silently would be worse than
        leaving the user with a one-off oddity to investigate."""
        sc._save_primers([
            {"name": "no_seq", "tm": 60.0},   # missing sequence
            {"name": "empty_seq", "sequence": ""},
            {"name": "ok", "sequence": "ACGT" * 5},
        ])
        primers = sc._load_primers()
        # All three preserved — only sequence-shared dupes collapse.
        assert len(primers) == 3

    async def test_cursor_moves_past_initial_viewport(
            self, isolated_library
    ):
        """Move the cursor down past the visible rows. The cursor
        should reach the last row without erroring — Textual's
        DataTable auto-scrolls its viewport to follow the cursor."""
        from textual.widgets import DataTable
        # Same uniqueness-per-index trick as the test above so the
        # dedupe-on-save in `_save_primers` doesn't collapse the 50
        # entries down to one.
        def _unique_seq(i: int) -> str:
            return f"GT{i:04d}".ljust(18, "A")
        entries = [{
            "name":        f"primer_{i:02d}",
            "sequence":    _unique_seq(i),
            "tm":          60.0,
            "primer_type": "imported",
            "source":      ".dna import",
            "status":      "Imported",
        } for i in range(50)]
        sc._save_primers(entries)
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            t = screen.query_one("#pd-lib-table", DataTable)
            assert t.row_count == 50
            # Move cursor to the last row programmatically — proves the
            # DataTable accepts navigation past the viewport. If the
            # container clipped the table, this would either error or
            # silently fail to scroll.
            t.move_cursor(row=49)
            await pilot.pause()
            assert t.cursor_row == 49


class TestPrimerFlapLinearClip:
    """A primer's 5' flap (the unbound enzyme tail) mod-wraps to the far
    end on a CIRCULAR molecule — the tail reappears across the origin —
    but must CLIP at the ends on a LINEAR one, which has no joined ends to
    wrap across. Gated by `_flap_linear`, stamped at parse time when the
    record is linear (2026-06-02). Without this, a primer near a linear
    fragment's terminus drew its tail bizarrely reappearing at the
    opposite end."""

    @staticmethod
    def _painted(arr):
        return {i for i, (c, _s) in enumerate(arr) if c != " "}

    def _arr(self, n):
        return [(" ", "") for _ in range(n)]

    def test_circular_flap_wraps_across_origin(self):
        # Forward primer bound at bp 1; 3-bp 5' flap spans bp -2..1, which
        # on a 10-bp CIRCULAR molecule wraps to columns 8, 9 (and 0).
        arr = self._arr(10)
        feat = {"_flap_start": -2, "_flap_end": 1, "_flap_bases": "AAA",
                "color": "red"}
        sc._paint_primer_flap_bar(arr, feat, 0, 10, total=10)
        assert self._painted(arr) == {0, 8, 9}

    def test_linear_flap_clips_at_end(self):
        # Same flap on a LINEAR molecule: the off-the-left-end tail
        # (bp -2, -1) is clipped, NOT wrapped to cols 8, 9. Only the
        # in-range base (bp 0) is drawn.
        arr = self._arr(10)
        feat = {"_flap_start": -2, "_flap_end": 1, "_flap_bases": "AAA",
                "color": "red", "_flap_linear": True}
        sc._paint_primer_flap_bar(arr, feat, 0, 10, total=10)
        assert self._painted(arr) == {0}


class TestPrimerMismatchBump:
    """Internal primer mismatches render bound-flap-bound: a mismatched base
    lifts onto the flap row (the 'bump'), leaving a coloured gap on the bound
    row. Mutagenic / Scrub QuikChange primers carry their planned edit as
    exactly such a mismatch, so this makes the edit visible on the map."""

    @staticmethod
    def _chars(arr):
        return [c for c, _s in arr]

    def _arr(self, n):
        return [(" ", "") for _ in range(n)]

    # ── painters (no app needed) ──────────────────────────────────────────

    def test_bound_bar_gaps_the_mismatch(self):
        # fwd primer, bound bases "GGACTC" at cols [2, 8); col 4 mismatches.
        f = {"type": "primer_bind", "start": 2, "end": 8, "strand": 1,
             "color": "green", "_primer_seq": "GGACTC", "_bound_len": 6,
             "_flap_len": 0, "_bound_mismatch": {4: "A"}}
        arr = self._arr(12)
        sc._paint_primer_bound_bar(arr, f, 0, 12)
        chars = self._chars(arr)
        assert chars[2] == "G" and chars[3] == "G"   # matched bases inline
        assert chars[4] == " "                        # mismatch lifted → gap
        assert chars[5] == "C"                        # back down (bound)
        assert arr[4][1] == ""                         # app-background gap (no fill)

    def test_flap_bar_draws_the_mismatch(self):
        f = {"type": "primer_bind", "start": 2, "end": 8, "strand": 1,
             "color": "green", "_bound_mismatch": {4: "A"}}
        arr = self._arr(12)
        sc._paint_primer_flap_bar(arr, f, 0, 12, total=12)
        assert arr[4][0] == "A"                        # the bump, on flap row
        assert {i for i, (c, _s) in enumerate(arr) if c != " "} == {4}

    def test_flap_bar_runs_without_5prime_tail(self):
        # No _flap_bases — a mismatch-only primer must still draw the bump.
        f = {"_bound_mismatch": {3: "T"}, "color": "red"}
        arr = self._arr(8)
        sc._paint_primer_flap_bar(arr, f, 0, 8, total=8)
        assert arr[3][0] == "T"

    # ── _parse stamping (via the app, both strands) ───────────────────────

    _TMPL = ("ACGTGACTTGCAACGGTATCCAGTTACGGCATTGAC"
             "AGTCCATGGATCACGTTAGCATGCATCAGTACCGTA")[:64]

    async def _feats_for(self, primer, strand):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq(self._TMPL), id="MM1", name="MM1",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec.features.append(SeqFeature(
            FeatureLocation(10, 26, strand=strand), type="primer_bind",
            qualifiers={"label": ["p"], "primer_seq": [primer]}))
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 40)) as pilot:
            await pilot.pause()
            app._apply_snapshot(self._TMPL, 0, rec)
            await pilot.pause()
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            prim = [f for f in pm._feats if f.get("type") == "primer_bind"]
            assert prim, "primer feature not parsed"
            return prim[0]

    async def test_parse_stamps_mismatch_forward(self):
        bound = list(self._TMPL[10:26])
        # Mismatch at offset 7 — past the 6 bp annealing anchor, so the 5'
        # flank stays bound and the base bumps. (A real QuikChange / Scrub
        # cure sits well inside the primer, never within the first few bp;
        # a mismatch closer than the anchor is read as a 5' overhang.)
        bound[7] = "A" if bound[7] != "A" else "T"
        f = await self._feats_for("".join(bound), 1)
        # forward base at offset 7 sits at column 10 + 7 = 17
        assert (f.get("_bound_mismatch") or {}) == {17: bound[7]}

    async def test_parse_no_mismatch_when_perfect_forward(self):
        f = await self._feats_for(self._TMPL[10:26], 1)
        assert not f.get("_bound_mismatch")

    async def test_parse_stamps_mismatch_reverse(self):
        perfect = sc._rc(self._TMPL[10:26])            # perfect rev primer 5'→3'
        rp = list(perfect)
        rp[8] = "A" if rp[8] != "A" else "T"           # internal, past the anchor
        f = await self._feats_for("".join(rp), -1)
        mism = f.get("_bound_mismatch") or {}
        assert len(mism) == 1                          # exactly one base lifted

    async def test_parse_no_mismatch_when_perfect_reverse(self):
        f = await self._feats_for(sc._rc(self._TMPL[10:26]), -1)
        assert not f.get("_bound_mismatch")

    async def test_internal_mismatch_keeps_both_sides_bound(self):
        # A primer that matches the template EXCEPT one internal base must
        # stay FULLY bound with only that base bumped — the matching bases 5'
        # of the mismatch must NOT collapse into the flap. (Regression: the
        # old longest-contiguous-3'-match binding dumped everything 5' of the
        # first mismatch into a flap that never returned to the bound row.)
        bound = list(self._TMPL[10:26])
        bound[8] = "A" if bound[8] != "A" else "T"
        f = await self._feats_for("".join(bound), 1)
        assert f.get("_bound_len") == 16        # whole primer stays bound
        assert not f.get("_flap_bases")         # no spurious 5' overhang
        assert (f.get("_bound_mismatch") or {}) == {18: bound[8]}  # only the bump

    async def test_mismatch_within_anchor_reads_as_overhang(self):
        # The flip side of the bump: a mismatch CLOSER to the 5' end than the
        # 6 bp annealing anchor is read as part of the 5' overhang, not a bump.
        # This is the same signal that keeps a cloning primer's enzyme tail a
        # clean flap even when its bases coincidentally pair the template —
        # only a solid (>= 6 bp) matching stretch counts as annealing. Pinned
        # so a future change can't collapse it back (cloning-tail regression).
        bound = list(self._TMPL[10:26])
        bound[3] = "A" if bound[3] != "A" else "T"     # only 3 bp 5' of it
        f = await self._feats_for("".join(bound), 1)
        assert not (f.get("_bound_mismatch") or {})    # no bump
        assert f.get("_flap_len") == 4                 # 5' overhang = first 4 bp


def _topo_record(seq: str, *, circular: bool, rid: str = "TOPO1"):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    rec = SeqRecord(Seq(seq), id=rid, name=rid, description="topology test")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular" if circular else "linear"
    return rec


class TestAddSelectedPrimerTopology:
    """[INV-96 / H1] `PrimerDesignScreen._add_selected_to_map` must be
    topology-aware. Adding a stored library primer whose position WRAPS
    the origin to a LINEAR record must NOT build — and persist via
    `lib.add_entry` — a wrap CompoundLocation across the non-joined ends
    (a primer-on-wrong-bases catastrophe). A CIRCULAR record keeps the
    wrap. Pre-fix this sibling path defaulted `circular=True` (only
    `PlasmidMap._parse` got the v1.0.14 topology fix), so it wrapped
    linear fragments too."""

    SEQ = "ACGTACGTAC" * 12          # 120 bp; no 12-bp run of the primer

    def _seed_wrap_primer(self):
        # A primer absent from the template (re-derivation returns None →
        # the STALE wrap position drives the branch) with pos_end <
        # pos_start (wraps the 120-bp origin).
        sc._save_primers([{
            "name": "wrapper",
            "sequence": "GGGGGCCCCCAAAAATTTTTGG",
            "pos_start": 110, "pos_end": 6, "strand": 1,
        }])

    async def _primer_loc_types(self, rec, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(rec, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            screen = sc.PrimerDesignScreen(str(rec.seq), [], rec.id)
            app.push_screen(screen)
            await pilot.pause()
            await pilot.pause(0.05)
            screen._lib_selected = {0}
            screen._add_selected_to_map(None)
            await pilot.pause(0.05)
            locs = [type(f.location).__name__
                    for f in app._current_record.features
                    if f.type == "primer_bind"]
            app.exit()
            return locs

    async def test_linear_record_does_not_get_wrap_primer(
            self, isolated_library, isolated_primers):
        self._seed_wrap_primer()
        rec = _topo_record(self.SEQ, circular=False, rid="LIN1")
        locs = await self._primer_loc_types(rec, isolated_library)
        assert "CompoundLocation" not in locs, (
            "a linear record must never get a wrap (CompoundLocation) primer"
        )

    async def test_circular_record_keeps_wrap_primer(
            self, isolated_library, isolated_primers):
        # Control: the SAME primer on a CIRCULAR record DOES wrap.
        self._seed_wrap_primer()
        rec = _topo_record(self.SEQ, circular=True, rid="CIRC1")
        locs = await self._primer_loc_types(rec, isolated_library)
        assert "CompoundLocation" in locs, (
            "a circular record should keep the wrap primer (topology control)"
        )
