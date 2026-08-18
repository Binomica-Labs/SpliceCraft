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
        sc._state._primers_cache = None
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
        # The skip collapses a DUPLICATED start codon, so the region has to
        # open on ATG for it to fire at all — otherwise there is nothing to
        # collapse and skipping would eat three of the user's own bases
        # ([INV-183]). The shared fixture is plain random DNA; plant it.
        base = self._valid_template()
        seq = base[:50] + "ATG" + base[53:]
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

    def test_goldenbraid_honors_target_tm(self):
        """The GB Tm field feeds `target_tm` into the binding-region
        optimizer: a higher target selects an arm at least as hot/long,
        and the requested target is recorded for the results panel."""
        seq = self._valid_template()
        lo = sc._design_gb_primers(seq, 50, 300, "Promoter", target_tm=52.0)
        hi = sc._design_gb_primers(seq, 50, 300, "Promoter", target_tm=68.0)
        assert "error" not in lo and "error" not in hi
        assert lo["target_tm"] == 52.0 and hi["target_tm"] == 68.0
        assert hi["fwd_tm"] >= lo["fwd_tm"]
        assert len(hi["fwd_binding"]) >= len(lo["fwd_binding"])


class TestPickBindingRegionNextBest:
    """`_pick_binding_region` returns the binding length whose Tm is
    CLOSEST to the target. When the target is unreachable within the
    18–25 nt window it returns the next-best (closest achievable) Tm —
    never a failure or Tm=0. This is the contract the Golden Braid Tm
    field relies on for "can't reach it → pick the closest"."""

    # 52 bp, mixed GC so Tm rises monotonically across 18–25 nt prefixes.
    _SEQ = "ACGTTGCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC"

    def test_in_range_picks_closest(self):
        primer3 = pytest.importorskip("primer3")
        bind, tm = sc._pick_binding_region(self._SEQ, target_tm=60.0)
        assert 18 <= len(bind) <= 25
        diffs = [abs(primer3.calc_tm(self._SEQ[:L]) - 60.0)
                 for L in range(18, 26)]
        assert abs(tm - 60.0) == pytest.approx(min(diffs))

    def test_next_best_when_target_too_high(self):
        # 95 C is unreachable in 18–25 nt — the longest (hottest, closest)
        # arm is the next-best.
        bind, tm = sc._pick_binding_region(self._SEQ, target_tm=95.0)
        assert len(bind) == 25
        assert tm > 0

    def test_next_best_when_target_too_low(self):
        # 20 C is below the whole window — the shortest (coolest, closest)
        # arm is the next-best.
        bind, _tm = sc._pick_binding_region(self._SEQ, target_tm=20.0)
        assert len(bind) == 18

    def test_default_target_is_60(self):
        assert (sc._pick_binding_region(self._SEQ)
                == sc._pick_binding_region(self._SEQ, target_tm=60.0))

    def test_short_seq_never_returns_zero_tm(self):
        # Below min_len the scan loop can't run; the defensive init must
        # still return a real (non-zero) Tm for whatever bases exist, not
        # a silent Tm=0 that downstream low-Tm checks would misread.
        _bind, tm = sc._pick_binding_region("ACGTACGTAC", target_tm=60.0)
        assert tm > 0


class TestPickBindingRegionAmbiguousBases:
    """Regression: `_pick_binding_region` used to catch only `ImportError`
    around primer3, but `primer3.calc_tm` RAISES `ValueError` on any
    non-ACGT base. A Plasmidsaurus consensus carries N at low-coverage
    positions, NCBI records carry IUPAC ambiguity codes, and the Synthesis
    composer lets you type IUPAC ([INV-53]) — so designing primers over
    such a region raised straight out of the worker instead of degrading
    to the 2+4 approximation the way every sibling Tm helper does."""

    _CLEAN = "ACGTTGCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC"

    def test_n_in_the_binding_window_does_not_raise(self):
        seq = "ACGTNNCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC"
        bind, tm = sc._pick_binding_region(seq, target_tm=60.0)
        assert 18 <= len(bind) <= 25
        assert tm > 0

    def test_iupac_in_the_binding_window_does_not_raise(self):
        seq = "ACGTRYCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC"
        bind, tm = sc._pick_binding_region(seq, target_tm=60.0)
        assert 18 <= len(bind) <= 25
        assert tm > 0

    def test_ambiguous_base_scores_the_at_gc_midpoint(self):
        # The fallback must not UNDERestimate an ambiguous arm into looking
        # unusable — N counts 3, between AT's 2 and GC's 4 (mirrors
        # `_primer_tm` / `_mut_tm`).
        _b, tm_n = sc._pick_binding_region("N" * 18, target_tm=60.0)
        _b, tm_at = sc._pick_binding_region("AT" * 9, target_tm=60.0)
        _b, tm_gc = sc._pick_binding_region("GC" * 9, target_tm=60.0)
        assert tm_at < tm_n < tm_gc

    def test_clean_template_still_uses_primer3(self):
        # The fallback must not silently take over the ACGT path.
        primer3 = pytest.importorskip("primer3")
        _bind, tm = sc._pick_binding_region(self._CLEAN, target_tm=60.0)
        assert any(tm == pytest.approx(primer3.calc_tm(self._CLEAN[:L]))
                   for L in range(18, 26))

    def test_generic_primer_design_survives_an_n(self):
        # End-to-end: the crash surfaced through the region designers.
        import random
        rng = random.Random(4)
        seq = "".join(rng.choice("ACGT") for _ in range(400))
        seq = seq[:50] + "NN" + seq[52:]
        res = sc._design_generic_primers(seq, 50, 250)
        assert res.get("fwd_seq") and res.get("rev_seq")


class TestPrimerOligoLengthCap:
    """The 50 bp total-oligo cap (2026-06-09 user report): designs grow the
    binding region to reach ~60 °C for low-GC (AT-rich) templates — the old
    fixed 25 bp binding cap stranded them at ~50 °C — while keeping the WHOLE
    oligo (5' tail + binding) within `_PRIMER_MAX_OLIGO_LEN` so synthesis
    stays cheap."""

    # ~30% GC (low-GC-host / codon-optimised style): a 25 bp binding tops out
    # near 53 °C, so reaching 60 °C requires the binding to grow to ~35 bp.
    _LOWGC = ("ATGAAAACATTAGAAAAATTAGCAGAAGAATTAGGTGTACCAAAATGGGTTATTAACGAT"
              "TTAGCAGAACAATTAGGTATTAAAGAAGCATTAGCAGATTTAGGTGAAGCATTAGAAAAA")

    def _design(self):
        tmpl = "ACACGTACGT" * 3 + self._LOWGC + "ACGTACACGT" * 3
        return sc._design_gb_primers(tmpl, 30, 30 + len(self._LOWGC), "Promoter")

    def test_binding_max_len_from_tail(self):
        assert sc._binding_max_len(0)  == sc._PRIMER_MAX_OLIGO_LEN      # no tail
        assert sc._binding_max_len(15) == sc._PRIMER_MAX_OLIGO_LEN - 15
        # An over-long tail never shrinks the binding below the 18 bp floor.
        assert sc._binding_max_len(40) == 18

    def test_gb_lowgc_binding_grows_to_reach_target(self):
        pytest.importorskip("primer3")   # accurate Tm needed for the °C asserts
        r = self._design()
        assert "error" not in r, r.get("error")
        p = r["pairs"][0]
        # Grew past the OLD 25 bp cap specifically to reach ~60 °C...
        assert len(p["fwd_binding"]) > 25, "binding didn't grow past the old cap"
        assert p["fwd_tm"] >= 58.0, f"fwd Tm still low: {p['fwd_tm']}"
        # ...without exceeding the total-oligo budget.
        assert len(p["fwd_full"]) <= sc._PRIMER_MAX_OLIGO_LEN
        assert len(p["rev_full"]) <= sc._PRIMER_MAX_OLIGO_LEN

    def test_all_gb_arms_within_oligo_cap(self):
        # SACRED budget: neither full primer may exceed the cap, ever.
        r = self._design()
        assert "error" not in r, r.get("error")
        for p in r["pairs"]:
            assert len(p["fwd_full"]) <= sc._PRIMER_MAX_OLIGO_LEN
            assert len(p["rev_full"]) <= sc._PRIMER_MAX_OLIGO_LEN

    def test_scrub_gb_shares_the_budget(self):
        # Scrub's binding cap is derived from the same total budget:
        # tail (pad + site + spacer) + binding == cap, and it grew past 25.
        tail = (len(sc._SCRUB_GB_PAD) + len(sc._SCRUB_GB_SITE)
                + len(sc._SCRUB_GB_SPACER))
        assert tail + sc._SCRUB_GB_BIND_MAX == sc._PRIMER_MAX_OLIGO_LEN
        assert sc._SCRUB_GB_BIND_MAX > 25


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
        monkeypatch.setattr(sc._state, "_PRIMERS_FILE", p)
        monkeypatch.setattr(sc._state, "_primers_cache", None, raising=False)

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

    @pytest.mark.asyncio
    async def test_lib_header_does_not_clip_right_edge(self):
        """The primer-library header must not overflow the right edge.
        The view toggle only flips `display` on a view CHANGE, so without a
        hidden-by-default the collections-action group rendered ALONGSIDE the
        primers group on first open and pushed the Close button off the right
        side of the terminal. Default (primers) view: collections buttons stay
        hidden and Close fits within the terminal width."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(185, 50)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            for _ in range(4):
                await pilot.pause()
            coll = app.screen.query_one("#pd-lib-coll-btns")
            close = app.screen.query_one("#btn-pd-close")
            row = app.screen.query_one("#pd-lib-hdr-row")
            assert not coll.display, \
                "collection-action buttons must start hidden in primers view"
            assert close.region.right <= 185, \
                f"Close clips the right edge (right={close.region.right})"
            assert row.region.right <= 185

    @pytest.mark.asyncio
    async def test_lib_header_close_reachable_on_narrow_terminal(self):
        """On a narrow terminal the primer-list buttons need more room than
        the panel has; the header scrolls horizontally and Close is docked
        right, so Close (and the contained row) stay on-screen rather than
        clipping off the right edge."""
        W = 130
        app = sc.PlasmidApp()
        async with app.run_test(size=(W, 50)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            for _ in range(4):
                await pilot.pause()
            close = app.screen.query_one("#btn-pd-close")
            row = app.screen.query_one("#pd-lib-hdr-row")
            assert close.region.width > 0 and close.region.right <= W, \
                f"Close not reachable at width {W} (right={close.region.right})"
            assert row.region.right <= W, \
                "header row overflows the terminal instead of scrolling"

    @pytest.mark.asyncio
    async def test_lib_header_collections_view_shows_its_actions(self):
        """Toggling to collections view must reveal its action buttons
        (the hidden-by-default is overridden by the toggle's inline display)
        and still fit within the terminal."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(185, 50)) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            for _ in range(4):
                await pilot.pause()
            screen._toggle_lib_view_mode(None)
            for _ in range(4):
                await pilot.pause()
            coll = app.screen.query_one("#pd-lib-coll-btns")
            prim = app.screen.query_one("#pd-lib-primers-btns")
            close = app.screen.query_one("#btn-pd-close")
            assert coll.display and not prim.display
            assert close.region.right <= 185


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


class TestPrimerTemplateRepickInv98:
    """[INV-98] Re-picking a library plasmid as the primer-design template
    must take the clean display name from the library entry, not the
    underscored GenBank LOCUS — else default primer names carry underscores."""

    _TERM = (180, 50)

    @pytest.mark.asyncio
    async def test_repick_uses_entry_name_not_locus(self, monkeypatch):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        gb = ("LOCUS       MY_PLASMID  12 bp ds-DNA circular SYN 01-JAN-2026\n"
              "ORIGIN\n        1 atgcatgcatgc\n//\n")
        entry = {"id": "p1", "name": "MY PLASMID", "gb_text": gb, "size": 12}
        monkeypatch.setattr(sc, "_load_library", lambda: [entry])
        app = sc.PlasmidApp()
        async with app.run_test(size=self._TERM) as pilot:
            await pilot.pause()
            await pilot.pause()
            # A current record is required for the designer to open.
            app._current_record = SeqRecord(Seq("ATGC" * 20), id="cur", name="cur")
            app._open_primer_design_impl(show_picker=False)
            await pilot.pause()
            await pilot.pause()
            screen = app.screen
            assert isinstance(screen, sc.PrimerDesignScreen)
            # Drive the picker callback directly with our entry id.
            orig_push = app.push_screen
            app.push_screen = lambda modal, callback=None, *a, **k: (
                callback("p1") if callback else None)
            try:
                screen._pick_plasmid(None)
                await pilot.pause()
            finally:
                app.push_screen = orig_push
            assert screen._plasmid_name == "MY PLASMID"
            assert "_" not in screen._plasmid_name

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

    async def test_two_scattered_internal_mismatches_both_bump(self):
        bound = list(self._TMPL[10:26])
        bound[7] = "A" if bound[7] != "A" else "T"
        bound[11] = "A" if bound[11] != "A" else "T"
        f = await self._feats_for("".join(bound), 1)
        assert (f.get("_bound_mismatch") or {}) == {17: bound[7], 21: bound[11]}
        assert f.get("_bound_len") == 16 and not f.get("_flap_bases")

    async def test_adjacent_three_base_cure_run_all_bump(self):
        # a triple-base Scrub cure (3 adjacent mismatches) — all three bump.
        bound = list(self._TMPL[10:26])
        for off in (7, 8, 9):
            bound[off] = "A" if bound[off] != "A" else "T"
        f = await self._feats_for("".join(bound), 1)
        assert set(f.get("_bound_mismatch") or {}) == {17, 18, 19}
        assert f.get("_bound_len") == 16 and not f.get("_flap_bases")

    async def test_arbitrary_count_internal_mismatches_all_bump(self):
        # No cap: five scattered internal mismatches all bump.
        bound = list(self._TMPL[10:26])
        offs = (6, 8, 10, 12, 14)
        for off in offs:
            bound[off] = "A" if bound[off] != "A" else "T"
        f = await self._feats_for("".join(bound), 1)
        assert set(f.get("_bound_mismatch") or {}) == {10 + o for o in offs}

    async def test_reverse_strand_multiple_mismatches_bump(self):
        perfect = sc._rc(self._TMPL[10:26])
        rp = list(perfect)
        for off in (7, 8):
            rp[off] = "A" if rp[off] != "A" else "T"
        f = await self._feats_for("".join(rp), -1)
        assert len(f.get("_bound_mismatch") or {}) == 2

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


class TestCloningPrimerFlapAfterPurify:
    """The unbound 5' flap of a cloning primer renders correctly ONLY once the
    insert digest purifies the primer pad out of the product (INV-119). With
    the pad still in the clone (the pre-fix carryover bug) the WHOLE primer —
    pad + binding — matches the product contiguously, so `_rederive_primer_
    binding` returns a full-length footprint and the flap collapses to zero
    (the primer looks fully bound). Purifying the off-cuts removes the pad, so
    only the 3' annealing arm matches and the pad re-emerges as the flap. No
    flap-render code changed — the flap was a symptom of the pad carryover."""

    BIND = "ATGAAACGTACGGCCTTAACG"          # 21 nt that anneals
    PRIMER = "GCGC" + "GTCGAC" + BIND        # GCGC pad + SalI + binding (31 nt)

    def test_pad_laden_product_binds_whole_primer_no_flap(self):
        # Pre-fix product still carries [pad][site] right before the binding.
        buggy = "TTTTTTTTTT" + self.PRIMER + "CCCCCCCCCC"
        rb = sc._rederive_primer_binding(self.PRIMER, 1, buggy, len(buggy),
                                         hint_start=0, circular=False)
        assert rb is not None and (rb[1] - rb[0]) == len(self.PRIMER), (
            "pad carryover makes the whole primer match → flap would vanish")

    def test_purified_product_leaves_pad_as_flap(self):
        # Fixed product: SalI site reforms at the junction, but the GCGC pad is
        # gone. The primer binds the reformed GTCGAC + arm; only the 4 nt pad
        # flaps.
        fixed = "AAAAAAAA" + "GTCGAC" + self.BIND + "TTTTTTTT"
        rb = sc._rederive_primer_binding(self.PRIMER, 1, fixed, len(fixed),
                                         hint_start=0, circular=False)
        assert rb is not None
        flap = len(self.PRIMER) - (rb[1] - rb[0])
        assert flap == 4, f"the 4 nt GCGC pad must read as the flap, got {flap}"


class TestRederiveStrictNotOffset:
    """`_rederive_primer_binding` must NOT slide a primer to a best-offset when
    the strict contiguous-suffix search fails — `_attach_pcr_primers_to_record`
    calls with hint_start=0, so a best-offset window at the origin could anchor
    a non-binding primer there (a mis-placed cloning primer is catastrophic). It
    returns None instead; the length-short *stored-feature* repair lives in
    `PlasmidMap._parse`."""

    def test_indel_primer_returns_none_not_offset_anchor(self):
        # A primer with no clean >=12 bp contiguous suffix (1 bp indel near the
        # 3' end) must come back None, NOT anchored at 0.
        primer = "GCGCCGTCTCAAATGACTGCATGCAGTACGTAGCTAGCAT"     # 40 bp, unique
        frag = primer[:36] + primer[37:] + "TTTTTTTTTTGGGGGGGGGG"  # drop base@36
        assert sc._rederive_primer_binding(
            primer, 1, frag, len(frag), hint_start=0, circular=False) is None

    def test_clean_full_match_uses_strict_path(self):
        primer = "GCGCCGTCTCAAATGACGTGCATCGATGCATCGTAGCATG"     # 40 bp
        frag = primer + "TTTTTTTTTTGGGGGGGGGGCCCCCCCCCC"
        rb = sc._rederive_primer_binding(primer, 1, frag, len(frag),
                                         hint_start=0, circular=False)
        assert rb == (0, 40), f"clean primer must anchor exactly, got {rb!r}"

    def test_foreign_primer_keeps_none(self):
        frag = "GCGCCGTCTCAAATGACGTACGTACGTACGTACGTTTTTTTTTTT"
        foreign = "ATATATATATATATATATATATATATATATATATAT"       # no real binding
        assert sc._rederive_primer_binding(
            foreign, 1, frag, len(frag), hint_start=0, circular=False) is None

    @pytest.mark.asyncio
    async def test_parse_reanchors_length_short_stored_feature(self):
        """A primer_bind feature whose stored location is 1 bp SHORTER than its
        primer_seq (a fragment built short of its primer) re-anchors flush in
        `_parse` — bound bar lays on its true bases, only the real indel shows
        as a couple of bumps — NOT every base slid one column (46 phantom
        mismatches)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        region = "GCGCCGTCTCAAATGACTGCATGCAGTACGTAGCTAGCATCGATCGATCAG"  # 51 bp
        primer = region[:len(region) - 3] + "G" + region[len(region) - 3:]  # 52
        frag = region + "TTTTTTTTTTGGGGGGGGGGCCCCCCCCCC"
        rec = SeqRecord(Seq(frag), id="SHORT", name="SHORT",
                        annotations={"molecule_type": "DNA",
                                     "topology": "linear"})
        rec.features.append(SeqFeature(
            FeatureLocation(0, len(region), strand=1), type="primer_bind",
            qualifiers={"label": ["P"], "primer_seq": [primer]}))
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 50)) as pilot:
            for _ in range(6):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                for _ in range(2):
                    await pilot.pause()
            app._apply_record(rec)
            for _ in range(6):
                await pilot.pause()
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            pf = [f for f in pm._feats if f.get("type") == "primer_bind"][0]
            assert pf.get("_bound_len") == len(primer)         # re-anchored full
            assert len(pf.get("_bound_mismatch") or {}) <= 6   # only the indel


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


class TestPrimerLibraryMarking:
    """Space marks the highlighted primer in the main Primers library, and a
    mark keeps the cursor row + scroll position — no jump to the bottom of the
    viewport — so marking many primers in a row doesn't jolt the list."""

    @staticmethod
    def _uniq_seq(i: int) -> str:
        """A unique 20 nt sequence per i (base-4 ACGT tail) so `_save_primers`'
        sequence-dedup keeps all N rows — a repeated sequence would collapse
        the table to a handful of rows and the scroll test would be moot."""
        tail = "".join("ACGT"[(i >> (2 * k)) & 3] for k in range(8))
        return "ACGTACGTACGT" + tail

    def _seed(self, n: int) -> None:
        sc._save_primers([
            {"name": f"P{i:03d}", "sequence": self._uniq_seq(i),
             "tm": 60.0, "primer_type": "generic", "source": "t",
             "pos_start": -1, "pos_end": -1, "strand": 1,
             "date": "2026-06-03", "status": "Designed"}
            for i in range(n)])

    async def test_space_marks_highlighted_keeps_cursor_and_scroll(
            self, isolated_primers):
        self._seed(60)
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 28)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            t = screen.query_one("#pd-lib-table", sc.DataTable)
            t.focus()
            await pilot.pause()
            t.move_cursor(row=45)
            await pilot.pause()
            await pilot.pause()
            cursor_before = t.cursor_row
            scroll_before = t.scroll_offset.y
            await pilot.press("space")
            await pilot.pause()
            await pilot.pause()
            assert len(screen._lib_selected) == 1            # exactly one marked
            assert t.cursor_row == cursor_before             # cursor didn't jump
            assert t.scroll_offset.y == scroll_before        # viewport stayed put

    async def test_space_marks_multiple_without_losing_prior(
            self, isolated_primers):
        # Invoke the handler directly with a synthetic Space (a real terminal
        # can double-fire a keypress, which would toggle a row off again — not
        # what we're testing here). We're checking that marks ACCUMULATE.
        class _SpaceKey:
            key = "space"
            def stop(self):
                pass
        self._seed(40)
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 28)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            t = screen.query_one("#pd-lib-table", sc.DataTable)
            assert t.row_count == 40                         # all rows present
            t.focus()
            await pilot.pause()
            counts = []
            for r in (8, 20, 33):
                t.move_cursor(row=r)
                await pilot.pause()
                screen.on_key(_SpaceKey())
                await pilot.pause()
                counts.append(len(screen._lib_selected))
            assert counts == [1, 2, 3]                       # each mark accumulates


class TestPrimerCsvExport:
    """`_export_primers_to_csv` writes a spreadsheet-pasteable order CSV
    (Name, Sequence, Length, Tm) for synthesis-company orders."""

    def test_writes_expected_rows(self, tmp_path):
        import csv as _csv
        primers = [
            {"name": "P1", "sequence": "acgtACGT", "tm": 58.3},
            {"name": "P2, weird", "sequence": "TTTT", "tm": None},
            {"name": "", "sequence": "GGGG"},          # empty name → fallback
            {"name": "noseq", "sequence": ""},         # no sequence → skipped
        ]
        out = tmp_path / "order.csv"
        res = sc._export_primers_to_csv(primers, out)
        assert res["count"] == 3 and res["path"] == str(out)
        rows = list(_csv.reader(out.open()))
        assert rows[0] == ["Name", "Sequence", "Length", "Tm"]
        assert rows[1] == ["P1", "ACGTACGT", "8", "58.3"]   # upper-cased, Tm 1dp
        assert rows[2] == ["P2, weird", "TTTT", "4", ""]    # comma stays intact
        assert rows[3] == ["primer", "GGGG", "4", ""]       # empty name fallback
        assert len(rows) == 4                                # no-seq row skipped

    def test_empty_raises(self, tmp_path):
        with pytest.raises(ValueError):
            sc._export_primers_to_csv([], tmp_path / "x.csv")
        with pytest.raises(ValueError):
            sc._export_primers_to_csv(
                [{"name": "a", "sequence": ""}], tmp_path / "y.csv")

    def test_refuses_malformed_oligo(self, tmp_path):
        """Catastrophic-class: refuse the whole export (naming the offender) if
        any oligo has non-DNA characters — never ship a wrong order."""
        primers = [{"name": "good", "sequence": "ACGT"},
                   {"name": "BadOne", "sequence": "ACGTXZ12"}]
        with pytest.raises(ValueError, match="BadOne"):
            sc._export_primers_to_csv(primers, tmp_path / "o.csv")

    def test_idt_format_columns(self, tmp_path):
        import csv as _csv
        primers = [
            {"name": "P1", "sequence": "ACGTACGT", "tm": 58.3},
            {"name": "P2", "sequence": "TTTT", "scale": "100nm",
             "purification": "HPLC"},                       # per-oligo override
        ]
        out = tmp_path / "idt.csv"
        res = sc._export_primers_to_csv(primers, out, order_format="idt")
        assert res["count"] == 2
        rows = list(_csv.reader(out.open()))
        assert rows[0] == ["Name", "Sequence", "Scale", "Purification"]
        assert rows[1] == ["P1", "ACGTACGT", "25nm", "STD"]     # defaults
        assert rows[2] == ["P2", "TTTT", "100nm", "HPLC"]       # overridden

    def test_unknown_format_raises(self, tmp_path):
        with pytest.raises(ValueError, match="order_format"):
            sc._export_primers_to_csv(
                [{"name": "a", "sequence": "ACGT"}], tmp_path / "x.csv",
                order_format="snapgene")

    def test_idt_still_refuses_malformed_oligo(self, tmp_path):
        # The catastrophic-class refusal holds regardless of order format.
        with pytest.raises(ValueError, match="Bad"):
            sc._export_primers_to_csv(
                [{"name": "Bad", "sequence": "ACGTXZ"}], tmp_path / "o.csv",
                order_format="idt")


class TestPrimerCsvImport:
    """`_import_primers_from_csv` round-trips the export and HARDENS against
    broken primers (skip + report, never apply)."""

    def test_round_trip(self, tmp_path):
        primers = [
            {"name": "Fwd-1", "sequence": "ACGTACGTACGT", "tm": 55.2},
            {"name": "Rev, A", "sequence": "TTTTGGGGCCCC", "tm": 61.0},
        ]
        out = tmp_path / "order.csv"
        sc._export_primers_to_csv(primers, out)
        res = sc._import_primers_from_csv(out)
        assert res["skipped"] == []
        got = {(p["name"], p["sequence"], p["tm"]) for p in res["primers"]}
        assert got == {("Fwd-1", "ACGTACGTACGT", 55.2),
                       ("Rev, A", "TTTTGGGGCCCC", 61.0)}
        assert all(p["primer_type"] == "imported" and p["status"] == "Imported"
                   for p in res["primers"])

    def test_skips_invalid_empty_and_duplicate(self, tmp_path):
        p = tmp_path / "mix.csv"
        p.write_text(
            "Name,Sequence,Length,Tm\n"
            "ok,ACGTACGT,8,58\n"
            "bad,ACGTXOXO,8,50\n"          # non-DNA → skipped
            "blank,,0,\n"                  # empty oligo → skipped
            "dup,ACGTACGT,8,58\n"          # duplicate oligo → skipped
        )
        res = sc._import_primers_from_csv(p)
        assert [x["name"] for x in res["primers"]] == ["ok"]
        assert len(res["skipped"]) == 3
        assert any("non-DNA" in s for s in res["skipped"])
        assert any("duplicate" in s for s in res["skipped"])

    def test_headerless_positional_tm_computed(self, tmp_path):
        # Elim-style 2-column (name, sequence), no header, no Tm → Tm computed.
        p = tmp_path / "elim.csv"
        p.write_text("MyFwd,ACGTACGTACGTACGTACGT\nMyRev,TTTTAAAACCCCGGGGTTTT\n")
        res = sc._import_primers_from_csv(p)
        assert [x["name"] for x in res["primers"]] == ["MyFwd", "MyRev"]
        assert all(isinstance(x["tm"], float) for x in res["primers"])

    def test_bom_tolerated(self, tmp_path):
        p = tmp_path / "bom.csv"
        p.write_bytes(b"\xef\xbb\xbfName,Sequence\nP,ACGTACGT\n")
        res = sc._import_primers_from_csv(p)
        assert [x["name"] for x in res["primers"]] == ["P"]

    def test_reordered_columns_and_tm_read(self, tmp_path):
        # Synthesis order forms vary the column order — detect by header name,
        # not position, and read the supplied Tm.
        p = tmp_path / "reorder.csv"
        p.write_text("Tm,Oligo,Name\n55.5,ACGTACGTACGT,Fwd\n")
        res = sc._import_primers_from_csv(p)
        assert len(res["primers"]) == 1
        prim = res["primers"][0]
        assert prim["name"] == "Fwd"
        assert prim["sequence"] == "ACGTACGTACGT"
        assert prim["tm"] == 55.5

    def test_empty_and_oversized_raise(self, tmp_path):
        empty = tmp_path / "e.csv"
        empty.write_text("   \n")
        with pytest.raises(ValueError):
            sc._import_primers_from_csv(empty)
        big = tmp_path / "big.csv"
        big.write_bytes(b"a" * (sc._PRIMER_CSV_MAX_BYTES + 100))
        with pytest.raises(ValueError, match="cap"):
            sc._import_primers_from_csv(big)


class TestPrimerCsvUiFlow:
    """The EXPORT / IMPORT buttons in the Primer Designer library panel open
    the CSV modals and round-trip primers into the active collection."""

    def _seed(self, n):
        sc._save_primers([
            {"name": f"P{i:03d}",
             "sequence": "ACGTACGTACGT" +
             "".join("ACGT"[(i >> (2 * k)) & 3] for k in range(8)),
             "tm": 60.0, "primer_type": "generic", "source": "t",
             "pos_start": -1, "pos_end": -1, "strand": 1,
             "date": "2026-06-03", "status": "Designed"}
            for i in range(n)])

    @pytest.mark.asyncio
    async def test_export_button_opens_modal_with_primers(
            self, isolated_primers):
        self._seed(3)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen._export_collection_csv(None)
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.PrimerCsvExportModal)
            assert len(app.screen._primers) == 3
            app.screen.dismiss(None)
            await pilot.pause()

    @pytest.mark.asyncio
    async def test_import_button_imports_valid_skips_invalid(
            self, isolated_primers, tmp_path):
        self._seed(2)
        csvp = tmp_path / "imp.csv"
        csvp.write_text("Name,Sequence\nNewOne,GGGGCCCCAAAATTTT\nBad,XOXOXO\n")
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            before = len(sc._load_primers())
            app.screen._import_collection_csv(None)
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.PrimerCsvImportModal)
            app.screen.dismiss(str(csvp))          # simulate picking the file
            for _ in range(5):
                await pilot.pause()
            after = sc._load_primers()
            assert len(after) == before + 1         # NewOne added, Bad skipped
            assert any(p["name"] == "NewOne" for p in after)

    @pytest.mark.asyncio
    async def test_collection_rename_and_delete(self, isolated_primers):
        self._seed(2)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            # A NON-empty 2nd collection (empty ones get pruned) so delete is
            # allowed — delete refuses the last remaining collection.
            colls = sc._load_primer_collections()
            colls.append({"name": "Order1", "description": "",
                          "primers": [{"name": "O1",
                                       "sequence": "TTTTGGGGCCCCAAAA",
                                       "tm": 60.0, "status": "Designed"}],
                          "saved": "2026-06-10"})
            sc._save_primer_collections(colls)
            assert len(sc._load_primer_collections()) == 2
            active = sc._get_active_primer_collection_name()
            screen._rename_primer_collection(None)
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.CollectionNameModal)
            app.screen.dismiss("Renamed-Coll")
            for _ in range(4):
                await pilot.pause()
            assert sc._get_active_primer_collection_name() == "Renamed-Coll"
            names = {c["name"] for c in sc._load_primer_collections()}
            assert "Renamed-Coll" in names and active not in names
            # Delete the (now renamed) active collection — NO-default confirm.
            screen._delete_primer_collection(None)
            for _ in range(4):
                await pilot.pause()
            assert isinstance(app.screen,
                              sc._PrimerCollectionDeleteConfirmModal)
            app.screen.dismiss(True)
            for _ in range(4):
                await pilot.pause()
            names2 = {c["name"] for c in sc._load_primer_collections()}
            assert "Renamed-Coll" not in names2 and "Order1" in names2

    @pytest.mark.asyncio
    async def test_bulk_delete_sequence_identity_count_and_full_removal(
            self, isolated_primers):
        """Bulk delete: the confirm enumerates the exact count, deletes by the
        unique SEQUENCE (so a shared NAME never over-deletes — the count can't
        go stale), and removes from BOTH the working library AND the active
        collection record."""
        sc._save_primers([
            {"name": "Dup", "sequence": "ACGTACGTACGTAAAA", "tm": 60.0,
             "status": "Designed"},
            {"name": "Dup", "sequence": "ACGTACGTACGTTTTT", "tm": 60.0,
             "status": "Designed"},          # SAME name, different oligo
            {"name": "Keep", "sequence": "GGGGCCCCGGGGCCCC", "tm": 60.0,
             "status": "Designed"},
        ])
        target = "ACGTACGTACGTAAAA"
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            primers = sc._load_primers()
            idx = next(i for i, p in enumerate(primers)
                       if p["sequence"] == target)
            screen._lib_selected = {idx}            # mark exactly ONE "Dup"
            screen._delete_marked_or_cursor()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.LibraryDeleteConfirmModal)
            assert "1 marked primer" in app.screen.entry_name   # exact count
            app.screen.dismiss(True)
            for _ in range(4):
                await pilot.pause()
            after = {p["sequence"] for p in sc._load_primers()}
            assert target not in after                          # marked gone
            assert "ACGTACGTACGTTTTT" in after                  # same-name kept
            assert "GGGGCCCCGGGGCCCC" in after
            assert len(after) == 2                              # no over-delete
            # …and gone from the active collection record too.
            coll = sc._find_primer_collection(
                sc._get_active_primer_collection_name())
            coll_seqs = {p["sequence"] for p in (coll.get("primers") or [])}
            assert target not in coll_seqs and len(coll_seqs) == 2

    @pytest.mark.asyncio
    async def test_cart_toggle_collect_and_export(self, isolated_primers):
        self._seed(3)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            # Mark 2 primers $ (cart). Space-cycling a primer onto $ sets the
            # persisted `in_cart` flag; set it directly here, then verify the
            # CART collect + export still work under the new mark model.
            entries = sc._load_primers()
            entries[0]["in_cart"] = True
            entries[1]["in_cart"] = True
            sc._save_primers(entries)
            screen._refresh_library_table()
            await pilot.pause()
            await pilot.pause()
            carted = sc._collect_cart_primers()
            assert len(carted) == 2 and all(p.get("in_cart") for p in carted)
            # persisted in the active collection record
            coll = sc._find_primer_collection(
                sc._get_active_primer_collection_name())
            assert sum(1 for p in coll["primers"] if p.get("in_cart")) == 2
            # Export Cart opens the CSV modal carrying the 2 carted primers.
            screen._export_cart_csv(None)
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc.PrimerCsvExportModal)
            assert len(app.screen._primers) == 2
            app.screen.dismiss(None)
            await pilot.pause()
            # un-cart all (Space-cycling off $ clears in_cart)
            entries = sc._load_primers()
            for e in entries:
                e.pop("in_cart", None)
            sc._save_primers(entries)
            await pilot.pause()
            assert sc._collect_cart_primers() == []

    @pytest.mark.asyncio
    async def test_move_marked_across_collections(self, isolated_primers):
        """The MOVE button moves the M-marked primers (Space-cycled to M) to a
        chosen collection. (Copy was dropped with the c/m/y keys.)"""
        self._seed(3)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            colls = sc._load_primer_collections()
            colls.append({"name": "Dest", "description": "",
                          "primers": [{"name": "X",
                                       "sequence": "GGGGGGGGGGGGCCCC",
                                       "tm": 60.0, "status": "Designed"}],
                          "saved": "2026-06-10"})
            sc._save_primer_collections(colls)
            before = sc._load_primers()
            moved = {before[0]["sequence"].upper(),
                     before[1]["sequence"].upper()}
            # M-mark 2 primers, press MOVE → Dest: source loses them.
            screen._move_marked = {0, 1}
            screen._move_marked_to_collection()
            await pilot.pause()
            await pilot.pause()
            assert isinstance(app.screen, sc._PrimerMoveCopyModal)
            app.screen.dismiss("Dest")
            for _ in range(4):
                await pilot.pause()
            after = {p["sequence"].upper() for p in sc._load_primers()}
            assert not (moved & after)                            # gone from src
            assert len(sc._load_primers()) == 1
            dest_seqs = {p["sequence"].upper()
                         for p in sc._find_primer_collection("Dest")["primers"]}
            assert moved <= dest_seqs

    @pytest.mark.asyncio
    async def test_space_cycles_primer_mark(self, isolated_primers):
        """Space cycles the cursor primer's SINGLE radio mark:
        none → ★ (select) → $ (cart) → M (move) → none. A primer holds at most
        one mark at a time (never ★ and M together)."""
        from textual.widgets import DataTable
        self._seed(2)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            t = screen.query_one("#pd-lib-table", DataTable)
            t.focus()
            t.move_cursor(row=0)
            await pilot.pause()
            pidx = screen._row_to_primer_idx[0]

            def mark():
                if pidx in screen._lib_selected:
                    return "star"
                ents = sc._load_primers()
                if 0 <= pidx < len(ents) and ents[pidx].get("in_cart"):
                    return "cart"
                if pidx in screen._move_marked:
                    return "move"
                return "none"

            assert mark() == "none"
            for expected in ("star", "cart", "move", "none"):
                await pilot.press("space")
                await pilot.pause()
                assert mark() == expected, \
                    f"after space expected {expected!r}, got {mark()!r}"
            # One mark at a time: cycle to ★ and confirm $/M are clear.
            await pilot.press("space")           # none → ★
            await pilot.pause()
            assert pidx in screen._lib_selected
            assert pidx not in screen._move_marked
            assert not sc._load_primers()[pidx].get("in_cart")


class TestPrimerCollectionWorkflow:
    """New-collection workflow + mark-state hygiene — the 'weird stuff making a
    new primer collection' report, plus the audit's `_move_marked` leaks."""

    @staticmethod
    def _seed(n):
        sc._save_primers([
            {"name": f"P{i:03d}",
             "sequence": "ACGTACGTACGT"
             + "".join("ACGT"[(i >> (2 * k)) & 3] for k in range(8)),
             "tm": 60.0, "primer_type": "generic", "source": "t",
             "pos_start": -1, "pos_end": -1, "strand": 1,
             "date": "2026-06-11", "status": "Designed"}
            for i in range(n)])

    @pytest.mark.asyncio
    async def test_new_collection_activates_and_lands_in_primers_view(
            self, isolated_primers):
        # Making a collection should ACTIVATE it and drop you into its primers
        # view — not leave it visible-but-inactive (the reported weirdness).
        self._seed(2)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause(); await pilot.pause(0.1)
            screen = app.screen
            screen._lib_view_mode = "collections"
            screen._new_primer_collection(None)
            await pilot.pause(); await pilot.pause()
            assert isinstance(app.screen, sc._PrimerCollectionNameModal)
            app.screen.dismiss("Batch2")
            for _ in range(5):
                await pilot.pause()
            assert sc._get_active_primer_collection_name() == "Batch2"
            assert screen._lib_view_mode == "primers"
            assert "Batch2" in {c.get("name")
                                for c in sc._load_primer_collections()}

    @pytest.mark.asyncio
    async def test_refresh_clips_stale_move_marks(self, isolated_primers):
        # A `_move_marked` index past the (shrunken) primer list must be culled
        # alongside `_lib_selected` — else MOVE acts on a phantom index.
        self._seed(2)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause(); await pilot.pause(0.1)
            screen = app.screen
            screen._move_marked = {0, 5}          # 5 is out of range (2 primers)
            screen._lib_selected = {9}            # also stale
            screen._refresh_library_table()
            await pilot.pause()
            assert screen._move_marked == {0}
            assert screen._lib_selected == set()

    @pytest.mark.asyncio
    async def test_shift_s_changes_status_of_marked_set(self, isolated_primers):
        # Shift+S acts on the ★-marked set (bulk), landing them all on ONE next
        # status (Designed → Ordered) derived from the cursor.
        from textual.widgets import DataTable
        self._seed(3)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause(); await pilot.pause(0.1)
            screen = app.screen
            t = screen.query_one("#pd-lib-table", DataTable)
            t.focus(); t.move_cursor(row=0)
            await pilot.pause()
            screen._lib_selected = {0, 1, 2}
            await pilot.press("S")
            await pilot.pause(); await pilot.pause()
            statuses = {p["status"] for p in sc._load_primers()}
            assert statuses == {"Ordered"}, \
                f"all ★-marked primers → one status, got {statuses}"

    @pytest.mark.asyncio
    async def test_space_cart_cycle_survives_dedupe_shrink(
            self, isolated_primers, monkeypatch):
        # A duplicate-sequence pair sitting in primers.json: cycling a LATER
        # primer to $ (cart) triggers _save_primers → dedupe, which drops the
        # dup and reindexes. The handler must re-clamp marks + rebuild the row
        # map (not leave a stale surgical paint), so a follow-up action hits the
        # right primer. Seed the cache directly — _save_primers would dedupe.
        from textual.widgets import DataTable
        dup = "ACGTACGTACGTACGTACGT"

        def _p(name, seq):
            return {"name": name, "sequence": seq, "tm": 60.0,
                    "primer_type": "generic", "source": "t", "pos_start": -1,
                    "pos_end": -1, "strand": 1, "date": "2026-06-11",
                    "status": "Designed"}
        monkeypatch.setattr(sc._state, "_primers_name_trim_done", True)
        monkeypatch.setattr(sc._state, "_primers_cache", [
            _p("P-A", "TTTTTTTTTTTTTTTTTTTT"),
            _p("P-DUP", dup),
            _p("P-B", dup),                    # duplicate sequence of P-DUP
            _p("P-C", "GGGGGGGGGGGGGGGGGGGG"),
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 40)) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause(); await pilot.pause(0.1)
            screen = app.screen
            t = screen.query_one("#pd-lib-table", DataTable)
            t.focus()
            assert len(sc._load_primers()) == 4        # dup present pre-save
            pc_row = next(
                r for r, pidx in enumerate(screen._row_to_primer_idx)
                if sc._load_primers()[pidx]["name"] == "P-C")
            t.move_cursor(row=pc_row)
            await pilot.pause()
            await pilot.press("space")                 # P-C: none → ★
            await pilot.pause()
            await pilot.press("space")                 # P-C: ★ → $ (dedupe shrinks)
            await pilot.pause(); await pilot.pause()
            now = sc._load_primers()
            assert len(now) == 3                       # the dup row was dropped
            # Row map rebuilt to match the deduped list (not left stale at 4).
            assert len(screen._row_to_primer_idx) == 3
            # The cart flag landed on the RIGHT primer, and only it.
            assert [pp["name"] for pp in now if pp.get("in_cart")] == ["P-C"]


class TestPrimerFlapBumped:
    """The unbound 5' flap renders on its OWN row, exactly ONE row farther from
    the DNA than the bound bar — bound stays bound, the overhang lifts off it —
    with ZERO gap: bound / flap / name are three CONSECUTIVE rows. Covers every
    orientation + flap length + mismatch-bump + origin-wrap / linear-clip case."""

    SEQ = "ACGTACGTAC" * 20          # 200 bp

    @staticmethod
    def _clear_seq_caches():
        # `_build_seq_text` keys its caches on id(feats); transient one-shot
        # feat lists in a big test run can reuse ids and collide. Clear first so
        # each render is fresh (the real app keeps one stable feats list).
        for c in ("_BUILD_SEQ_CACHE", "_CHUNK_STATIC_CACHE",
                  "_CHUNK_OVERLAY_CACHE", "_CHUNK_LAYOUT_CACHE"):
            cache = getattr(sc, c, None)
            if cache is not None:
                cache.clear()

    def _rows(self, feat):
        self._clear_seq_caches()
        txt = str(sc._build_seq_text(self.SEQ, [feat], line_width=len(self.SEQ)))
        return txt.splitlines()

    @staticmethod
    def _idx(rows, needle):
        return next((i for i, r in enumerate(rows) if needle in r), -1)

    def _assert_bumped(self, feat, arrow, flap_bases, label):
        """Bound stays on its row; the flap sits EXACTLY one row farther from the
        DNA than the bound (arrow) row; the name sits one row past the flap —
        three consecutive rows, ZERO gap — and the flap is NOT on the bound row."""
        rows = self._rows(feat)
        a  = self._idx(rows, arrow)             # bound bar (arrow) row
        f  = self._idx(rows, flap_bases)        # flap row
        nm = self._idx(rows, label)             # name row
        assert a >= 0 and f >= 0 and nm >= 0, \
            f"missing row (arrow={a} flap={f} name={nm})"
        # Bound is bound: the flap is NOT flush on the bound row.
        assert flap_bases not in rows[a], f"flap still on the bound row: {rows[a]!r}"
        step = -1 if arrow == "▶" else +1       # fwd lifts up, rev lifts down
        # Flap exactly ONE row out from the bound bar — zero gap.
        assert f == a + step, f"flap not 1 row out (arrow@{a} flap@{f})"
        # Name immediately past the flap — three consecutive rows, zero gap.
        assert nm == f + step, f"gap between flap and name (flap@{f} name@{nm})"

    @staticmethod
    def _fwd(start, end, flap_len, mism=None):
        fl = "GCGCCGTCTCT"[:flap_len]
        f = {"type": "primer_bind", "start": start, "end": end, "strand": 1,
             "label": "F", "color": "#ffaa00",
             "_primer_seq": fl + "A" * (end - start), "_bound_len": end - start,
             "_flap_len": flap_len, "_flap_bases": fl,
             "_flap_start": start - flap_len, "_flap_end": start}
        if mism:
            f["_bound_mismatch"] = mism
        return f

    @staticmethod
    def _rev(start, end, flap_len, mism=None):
        fl = "GCGCCGTCTCT"[:flap_len]
        f = {"type": "primer_bind", "start": start, "end": end, "strand": -1,
             "label": "R", "color": "#00aaff",
             "_primer_seq": fl + "T" * (end - start), "_bound_len": end - start,
             "_flap_len": flap_len, "_flap_bases": fl[::-1],
             "_flap_start": end, "_flap_end": end + flap_len}
        if mism:
            f["_bound_mismatch"] = mism
        return f

    def test_forward_flap_bumped(self):
        self._assert_bumped(self._fwd(40, 65, 11), "▶", "GCGCCGTCTCT", "F")

    def test_reverse_flap_bumped(self):
        self._assert_bumped(self._rev(40, 65, 11), "◀", "GCGCCGTCTCT"[::-1], "R")

    def test_forward_short_flap_bumped(self):
        self._assert_bumped(self._fwd(40, 65, 4), "▶", "GCGC", "F")  # pad-only

    def test_reverse_short_flap_bumped(self):
        self._assert_bumped(self._rev(40, 65, 4), "◀", "GCGC"[::-1], "R")

    def test_forward_long_flap_bumped(self):
        self._assert_bumped(self._fwd(60, 90, 11), "▶", "GCGCCGTCTCT", "F")

    def test_reverse_long_flap_bumped(self):
        self._assert_bumped(self._rev(60, 90, 11), "◀", "GCGCCGTCTCT"[::-1], "R")

    def test_reverse_flap_with_mismatch_bumped(self):
        # domestication-style: flap tail + an internal bump both ride the flap
        # row, exactly one row off the bound bar.
        self._assert_bumped(self._rev(40, 65, 11, mism={50: "G"}),
                            "◀", "GCGCCGTCTCT"[::-1], "R")

    def test_forward_flap_with_mismatch_bumped(self):
        self._assert_bumped(self._fwd(40, 65, 11, mism={50: "G"}),
                            "▶", "GCGCCGTCTCT", "F")

    def test_bump_only_primer_lifts_to_its_own_row(self):
        # A mutagenic primer with NO 5' tail but an internal cure: the bump
        # still lifts onto its own row, one off the bound bar (height 3 → name
        # two rows out), zero gap.
        f = self._fwd(40, 65, 0, mism={50: "G"})
        for k in ("_flap_bases", "_flap_start", "_flap_end", "_flap_len"):
            f.pop(k, None)
        rows = self._rows(f)
        a, nm = self._idx(rows, "▶"), self._idx(rows, "F")
        assert a >= 0 and nm >= 0, (a, nm)
        assert nm == a - 2, f"bump-only primer not height 3 (arrow@{a} name@{nm})"

    def test_circular_origin_wrap_flap_bumped(self):
        # Forward primer near bp 0: the 5' flap wraps the origin (negative
        # start). It still lifts to its own row off the bound bar — no crash.
        f = self._fwd(3, 28, 11)
        f["_flap_start"] = -8           # [-8, 3) wraps the origin
        rows = self._rows(f)
        a = self._idx(rows, "▶")
        assert a >= 0                   # rendered, no crash
        # The wrapped flap shows on the row ABOVE the bound bar, never on it.
        assert "GCG" in rows[a - 1] and "GCG" not in rows[a], \
            (rows[a - 1], rows[a])

    def test_linear_end_flap_clips_no_crash(self):
        # Reverse primer at the very END of a LINEAR fragment: the 5' flap
        # dangles past the end and CLIPS (not wraps) without a crash.
        n = len(self.SEQ)
        f = self._rev(n - 20, n, 11)    # flap [n, n+11) past the end
        f["_flap_linear"] = True
        rows = self._rows(f)
        assert self._idx(rows, "◀") >= 0   # bound bar rendered, flap clipped


class TestPrimerSeqGenBankWrapRoundTrip:
    """A long `/primer_seq` qualifier must survive a GenBank round-trip
    WITHOUT a space jammed into the sequence.

    BioPython's GenBank writer wraps a long qualifier value across lines
    and its parser rejoins the continuation with a space — fine for a
    free-text `/note`, corrupting for a DNA sequence. Any primer longer
    than the ~58-char value-wrap width round-tripped through `.gb` with a
    literal space mid-sequence, which the seq panel drew as a phantom gap
    that shifted the primer off the template. `_gb_text_to_record` now
    strips it on load (`_repair_wrapped_primer_seqs`)."""

    @staticmethod
    def _record_with_primer(primer_seq: str):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        # A template the primer binds, padded so it's a real plasmid.
        seq = ("AGCT" * 20) + primer_seq + ("TTAA" * 20)
        rec = SeqRecord(Seq(seq), id="rt", name="rt")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        start = 80
        rec.features.append(SeqFeature(
            FeatureLocation(start, start + len(primer_seq), strand=1),
            type="primer_bind",
            qualifiers={"label": ["p"], "primer_seq": [primer_seq]},
        ))
        return rec

    def test_long_primer_seq_has_no_injected_space(self):
        # 64 nt — comfortably past the GenBank value-wrap width, so the
        # writer splits it across two lines.
        primer = "GCGCCGTCTCAAATGAATAAATGTATTCCAATGATAATTAATGGAATGATCAGATTTTGATAA"
        assert len(primer) > 58
        gb = sc._record_to_gb_text(self._record_with_primer(primer))
        # The on-disk GenBank legitimately wraps the value across lines…
        assert "\n" in gb
        # …but the parsed-back primer must be whitespace-free and intact.
        rec2 = sc._gb_text_to_record(gb, cache=False)
        pf = next(f for f in rec2.features if f.type == "primer_bind")
        got = pf.qualifiers["primer_seq"][0]
        assert " " not in got and "\n" not in got, repr(got)
        assert got == primer

    def test_short_primer_seq_is_unchanged(self):
        primer = "GTAAAACGACGGCCAGT"          # 17 nt, never wraps
        gb = sc._record_to_gb_text(self._record_with_primer(primer))
        rec2 = sc._gb_text_to_record(gb, cache=False)
        pf = next(f for f in rec2.features if f.type == "primer_bind")
        assert pf.qualifiers["primer_seq"][0] == primer

    def test_repair_helper_strips_internal_whitespace(self):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        rec = SeqRecord(Seq("ACGT" * 30), id="x", name="x")
        rec.features.append(SeqFeature(
            FeatureLocation(0, 10, strand=1), type="primer_bind",
            qualifiers={"primer_seq": ["GCGCCGTCTCAAATG GGAA TGAT"]}))
        sc._repair_wrapped_primer_seqs(rec)
        assert rec.features[0].qualifiers["primer_seq"][0] == "GCGCCGTCTCAAATGGGAATGAT"


class TestPrimerSeqRenderChokepoint:
    """The render path (`PlasmidMap._parse`) is the universal chokepoint:
    ANY primer that reaches the canvas — gb-loaded, in-memory, agent-built,
    .dna-imported — must have its `/primer_seq` whitespace collapsed before
    it drives the bound/flap split, so no source can paint a phantom gap."""

    @staticmethod
    def _record(primer_seq_with_space: str):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        clean = "".join(primer_seq_with_space.split())
        seq = ("AGCT" * 25) + clean + ("TTAA" * 25)
        rec = SeqRecord(Seq(seq), id="r", name="r")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        start = 100
        # DIRTY qualifier carried directly — NOT via a GenBank round-trip — so
        # the render guard is exercised in isolation from the load-time repair.
        rec.features.append(SeqFeature(
            FeatureLocation(start, start + len(clean), strand=1),
            type="primer_bind",
            qualifiers={"label": ["p"], "primer_seq": [primer_seq_with_space]},
        ))
        return rec, clean

    async def _parsed_primer(self, dirty):
        rec, clean = self._record(dirty)
        sc.PlasmidApp._skip_seed = True
        sc.PlasmidApp._skip_snapshot = True
        sc.PlasmidApp._skip_update_check = True
        sc.PlasmidApp._preload_demo_record = None
        app = sc.PlasmidApp()
        async with app.run_test(size=(220, 50)) as pilot:
            for _ in range(5):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                await pilot.pause()
            app._apply_record(rec, clear_undo=True)
            for _ in range(3):
                await pilot.pause()
            pm = app.query_one("#plasmid-map", sc.PlasmidMap)
            f = next(x for x in pm._feats if x.get("type") == "primer_bind")
            return dict(f), clean

    @pytest.mark.asyncio
    async def test_render_strips_internal_space(self):
        dirty = "GCGCCGTCTCAAATGAATAAATGTATTCCAATGATAATTAATGGAA TGAT"
        f, clean = await self._parsed_primer(dirty)
        ps = f.get("_primer_seq", "")
        assert " " not in ps and "\n" not in ps, repr(ps)
        assert ps == clean
        # A perfectly-matching primer must produce NO bound-bar gap.
        assert not f.get("_bound_mismatch"), f.get("_bound_mismatch")

    @pytest.mark.asyncio
    async def test_render_strips_newline_and_tab(self):
        dirty = "GCGCCGTCTCAAATG\tAATAAATGTATTCCAATGAT\nAATTAATGGAATGAT"
        f, clean = await self._parsed_primer(dirty)
        ps = f.get("_primer_seq", "")
        assert ps == clean and not any(c.isspace() for c in ps)


class TestCartReset:
    """A successful cart-CSV export resets the cart: `_clear_all_cart_marks`
    strips every $ (in_cart) flag — across all collections (the cart source)
    and the primer library — so the next order starts empty. Idempotent."""

    def test_clears_cart_marks_and_is_idempotent(self, isolated_primers):
        sc._save_primer_collections([
            {"name": "Main", "primers": [
                {"name": "P1", "sequence": "AAAACCCCGGGGTTTT", "in_cart": True},
                {"name": "P3", "sequence": "GGGGTTTTAAAACCCC", "in_cart": True},
                {"name": "P2", "sequence": "TTTTGGGGCCCCAAAA"},   # not in cart
            ]},
        ])
        assert len(sc._collect_cart_primers()) == 2     # the pending order
        n = sc._clear_all_cart_marks()
        assert n == 2                                    # distinct cleared
        assert sc._collect_cart_primers() == []          # cart now empty
        for c in sc._load_primer_collections():
            assert all(not p.get("in_cart") for p in (c.get("primers") or []))
        # Idempotent — clearing an already-empty cart is a no-op.
        assert sc._clear_all_cart_marks() == 0

    def test_empty_cart_is_noop(self, isolated_primers):
        sc._save_primer_collections([
            {"name": "Main", "primers": [
                {"name": "P1", "sequence": "AAAACCCCGGGGTTTT"}]}])
        assert sc._clear_all_cart_marks() == 0
        assert sc._collect_cart_primers() == []

    def test_collect_includes_library_marks_not_just_collections(
            self, isolated_primers):
        """Regression (2026-06-13): the $ mark is set on the primer LIBRARY, so
        the export must collect from the library too — not only collections —
        else a marked primer that isn't in a scanned collection silently drops
        from the CSV (yet `_clear_all_cart_marks` would still clear it, so the
        mark just vanished). Now symmetric with the clear."""
        sc._save_primer_collections([])
        sc._set_active_primer_collection_name("")          # no mirror
        sc._save_primers([
            {"name": "L1", "sequence": "ACGTACGTACGTACGTACGT", "in_cart": True},
            {"name": "L2", "sequence": "TTTTGGGGCCCCAAAATTTT"},
        ])
        assert {p["name"] for p in sc._collect_cart_primers()} == {"L1"}, \
            "library $ mark not collected (collections-only scan regressed)"

    def test_set_primer_cart_flag_atomic_by_sequence(self, isolated_primers):
        """`_set_primer_cart_flag` sets/clears the $ flag BY SEQUENCE (case-
        insensitive), persists it, and surfaces it in the cart; unknown
        sequence → False. This is the race-safe persist behind the Space-cycle
        $ toggle (a stale whole-list save could clobber a concurrent mark)."""
        sc._save_primer_collections([])
        sc._set_active_primer_collection_name("")
        sc._save_primers([{"name": "A", "sequence": "ACGTACGTACGTACGTACGT"}])
        assert sc._set_primer_cart_flag("acgtacgtacgtacgtacgt", True) is True
        assert {p["name"] for p in sc._collect_cart_primers()} == {"A"}
        assert any(p.get("in_cart") for p in sc._load_primers())
        assert sc._set_primer_cart_flag("ACGTACGTACGTACGTACGT", False) is True
        assert sc._collect_cart_primers() == []
        assert sc._set_primer_cart_flag("NOSUCHSEQUENCE", True) is False


class TestPrimerSearchAndMarks:
    """Primer-library search bar (#pd-lib-search): fuzzy name filter, marks
    (★ select / M move — full-list-index keyed) survive filtering with no
    staleness, and CLEAR resets every mark."""

    async def test_search_filters_marks_survive_and_clear(self, isolated_primers):
        sc._save_primers([
            {"name": "alpha-DOM-1-F", "sequence": "AAAACCCCGGGGTTTTAAAA"},
            {"name": "beta-DOM-1-F",  "sequence": "CCCCGGGGTTTTAAAACCCC"},
            {"name": "gamma-DOM-1-F", "sequence": "GGGGTTTTAAAACCCCGGGG"},
        ])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 50, [], "t"))
            await pilot.pause()
            await pilot.pause(0.1)
            screen = app.screen
            screen._refresh_library_table()
            await pilot.pause()
            t = screen.query_one("#pd-lib-table", sc.DataTable)
            assert t.row_count == 3
            # Mark "alpha" (★) by its FULL-list index, then filter it OUT.
            primers = sc._load_primers()
            alpha_idx = next(i for i, p in enumerate(primers)
                             if p["name"] == "alpha-DOM-1-F")
            screen._lib_selected = {alpha_idx}
            screen._filter_text = "beta"
            screen._refresh_library_table()
            await pilot.pause()
            assert t.row_count == 1                 # only beta shows
            assert alpha_idx in screen._lib_selected  # mark survived the filter
            # Clearing the filter brings alpha back, STILL marked (no staleness).
            screen._filter_text = ""
            screen._refresh_library_table()
            await pilot.pause()
            assert t.row_count == 3
            assert alpha_idx in screen._lib_selected
            # A 0-match filter is safe (no out-of-range cursor move).
            screen._filter_text = "zzzznomatch"
            screen._refresh_library_table()
            await pilot.pause()
            assert t.row_count == 0
            # CLEAR resets the transient marks.
            screen._filter_text = ""
            screen._move_marked = {alpha_idx}
            screen._clear_primer_marks(None)
            await pilot.pause()
            assert screen._lib_selected == set() and screen._move_marked == set()
