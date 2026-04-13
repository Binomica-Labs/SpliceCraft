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


@pytest.fixture
def isolated_primers(tmp_path, monkeypatch):
    tmp_p = tmp_path / "primers.json"
    monkeypatch.setattr(sc, "_PRIMERS_FILE", tmp_p)
    monkeypatch.setattr(sc, "_primers_cache", None)
    return tmp_p


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
        assert parsed[0]["name"] == "x"

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
        # Long enough for GB (>= 18 bp) and Detection (>= product_min = 450).
        return ("A" * 50 + "ATGAAACGTGATTTAGCCGTTAA" * 40 + "T" * 50)

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
        r = sc._design_gb_primers(seq, 50, 300, "CDS")
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

    def test_goldenbraid_flags_internal_bsai(self):
        """Regression guard: a GGTCTC inside the insert is self-domesticating
        and must be flagged, not silently passed through."""
        core = "ATGAAACGTGATTTAGCC" * 5   # 90 bp, no BsaI
        with_fwd_bsai = core + "GGTCTC" + core
        r = sc._design_gb_primers(with_fwd_bsai, 0, len(with_fwd_bsai), "CDS")
        assert "error" in r
        assert "BsaI" in r["error"]

    def test_goldenbraid_flags_reverse_bsai(self):
        """GAGACC (RC of GGTCTC) on the top strand also fragments the part."""
        core = "ATGAAACGTGATTTAGCC" * 5
        with_rev_bsai = core + "GAGACC" + core
        r = sc._design_gb_primers(with_rev_bsai, 0, len(with_rev_bsai), "CDS")
        assert "error" in r
        assert "BsaI" in r["error"]

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
        from textual.widgets import RadioSet
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            rs = app.screen.query_one("#pd-mode-radio", RadioSet)
            assert rs.pressed_button is not None
            assert rs.pressed_button.id == "rb-detection"
            assert screen._current_mode() == "detection"

    async def test_all_four_modes_have_radio_buttons(self):
        from textual.widgets import RadioButton
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            rb_ids = {rb.id for rb in app.screen.query(RadioButton)}
            assert {
                "rb-detection", "rb-cloning",
                "rb-goldenbraid", "rb-generic",
            } <= rb_ids

    async def test_switching_mode_changes_current_mode(self):
        """Selecting a different RadioButton changes _current_mode() and
        swaps the visible parameter panel."""
        from textual.widgets import RadioSet
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

    async def test_book_layout_split(self):
        """Open-book split: left page owns input (template / mode /
        parameters); right page owns output (results + primer library).
        Results sits above Library so freshly-designed primers appear
        next to where you'll save them."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            book = app.screen.query_one("#pd-book")
            page_ids = [c.id for c in book.children]
            assert page_ids == ["pd-left-page", "pd-right-page"]

            # Left page: three input sections in order (no results).
            left = app.screen.query_one("#pd-left-page")
            left_ids = [c.id for c in left.children]
            assert "pd-results-section" not in left_ids
            idx_t = left_ids.index("pd-template-section")
            idx_m = left_ids.index("pd-mode-section")
            idx_p = left_ids.index("pd-params-section")
            assert idx_t < idx_m < idx_p

            # Right page: results above library.
            right = app.screen.query_one("#pd-right-page")
            right_ids = [c.id for c in right.children]
            idx_r  = right_ids.index("pd-results-section")
            idx_lh = right_ids.index("pd-lib-hdr-row")
            idx_lt = right_ids.index("pd-lib-table")
            assert idx_r < idx_lh < idx_lt

    async def test_results_section_has_name_and_save(self):
        """Regression guard: primer-name Inputs + Save button live
        INSIDE the results section, not as a separate row above it."""
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            results_section = app.screen.query_one("#pd-results-section")
            ids_inside = {w.id for w in results_section.walk_children()}
            assert "pd-fwd-name" in ids_inside
            assert "pd-rev-name" in ids_inside
            assert "btn-pd-save" in ids_inside

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
        # Need a wrap region free of internal BsaI sites. Build a
        # deterministic sequence and avoid GGTCTC / GAGACC.
        import random
        rng = random.Random(0xCAFE)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        # Ensure no BsaI sites by scanning
        if "GGTCTC" in seq or "GAGACC" in seq:
            seq = seq.replace("GGTCTC", "GCTGTC").replace("GAGACC", "GACCGA")
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
    async def test_wrap_hint_present_on_mount(self):
        """The hint telling users they can enter start > end for wrap
        regions should always be visible on mount."""
        from textual.widgets import Static
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 800, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            hint = app.screen.query_one("#pd-wrap-hint", Static)
            rendered = str(hint.render())
            assert "Start > End" in rendered
            assert "origin" in rendered

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
