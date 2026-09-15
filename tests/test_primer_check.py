"""
Primer Check — in-silico PCR across the library (PrimerDesignScreen 5th tab).

Covers the pure biology engine and the UI wiring:

  * `_circ_slice` / `_arc_intervals`            — wrap-aware slicing helpers
  * `_primer_binding_sites`                      — 3'-anchored fuzzy binding
  * `_insilico_pcr_amplicons`                    — amplicon pairing geometry,
                                                   incl. EXACT-primer parity with
                                                   `_simulate_pcr` (bp-for-bp)
  * `_amplicon_feature_summary`                  — feature an amplicon covers
  * `_primer_check_confidence`                   — identity → badge tiers
  * PrimerDesignScreen "Primer Check" tab        — present, mode-switch, scan e2e

The engine is deliberately stricter at the 3' end than `_search_subsequence`
and more tolerant in the 5' region than `_simulate_pcr` — so it answers
"would this primer prime here, and how well?".
"""
import asyncio
import random

import pytest

import splicecraft as sc

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


def _rand_dna(n: int, seed: int) -> str:
    return "".join(random.Random(seed).choices("ACGT", k=n))


TERMINAL_SIZE = (140, 42)


# ═══════════════════════════════════════════════════════════════════════════════
# Pure slicing helpers
# ═══════════════════════════════════════════════════════════════════════════════

class TestCircSlice:
    def test_no_wrap(self):
        assert sc._circ_slice("ABCDEFGH", 2, 3, 8) == "CDE"

    def test_wrap(self):
        assert sc._circ_slice("ABCDEFGH", 6, 4, 8) == "GHAB"

    def test_negative_start(self):
        # -2 mod 8 == 6 → same as the wrap case
        assert sc._circ_slice("ABCDEFGH", -2, 4, 8) == "GHAB"

    def test_degenerate(self):
        assert sc._circ_slice("", 0, 4, 8) == ""
        assert sc._circ_slice("ABC", 0, 0, 3) == ""


class TestArcIntervals:
    def test_no_wrap(self):
        assert sc._arc_intervals(100, 50, 600) == [(100, 150)]

    def test_wrap(self):
        assert sc._arc_intervals(580, 40, 600) == [(580, 600), (0, 20)]

    def test_full_circle(self):
        assert sc._arc_intervals(0, 600, 600) == [(0, 600)]

    def test_degenerate(self):
        assert sc._arc_intervals(0, 0, 600) == []
        assert sc._arc_intervals(0, 10, 0) == []


# ═══════════════════════════════════════════════════════════════════════════════
# Binding sites
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrimerBindingSites:
    T = _rand_dna(600, seed=1)

    def test_forward_exact(self):
        fwd = self.T[100:122]
        sites = sc._primer_binding_sites(fwd, self.T, len(self.T),
                                         circular=False)
        assert len(sites) == 1
        s = sites[0]
        assert s["strand"] == 1
        assert s["foot_start"] == 100
        assert s["ident_pct"] == 100.0
        assert s["mismatches"] == 0

    def test_reverse_exact(self):
        region = self.T[400:422]
        rev = sc._rc(region)             # reverse primer as ordered, 5'→3'
        sites = sc._primer_binding_sites(rev, self.T, len(self.T),
                                         circular=False)
        assert len(sites) == 1
        s = sites[0]
        assert s["strand"] == -1
        assert s["foot_start"] == 400    # 3'/left edge of the footprint
        assert s["ident_pct"] == 100.0

    def test_five_prime_tail_lowers_identity_but_still_found(self):
        # Cloning-style 5' flap that does NOT anneal — the 3' seed still
        # binds, so the site is found but identity drops below 100%.
        fwd = "GGGGGCCCCC" + self.T[100:122]      # 10 nt tail + 22 nt binding
        sites = sc._primer_binding_sites(fwd, self.T, len(self.T),
                                         circular=False)
        assert sites, "tailed primer should still bind via its 3' seed"
        s = sites[0]
        assert s["strand"] == 1
        assert 0.0 < s["ident_pct"] < 100.0
        assert s["mismatches"] > 0

    def test_three_prime_mismatch_not_found(self):
        # A mismatch inside the 3' seed kills the binding call — that primer
        # would not extend, so it is correctly NOT reported.
        fwd = list(self.T[100:122])
        fwd[-1] = {"A": "C", "C": "A", "G": "T", "T": "G"}[fwd[-1]]
        sites = sc._primer_binding_sites("".join(fwd), self.T, len(self.T),
                                         circular=False)
        assert sites == []

    def test_circular_wrap_site(self):
        # A primer that straddles the origin is found on a circular template.
        primer = self.T[590:] + self.T[:12]      # 10 + 12 = 22 nt, wraps origin
        sites = sc._primer_binding_sites(primer, self.T, len(self.T),
                                         circular=True)
        assert any(s["foot_start"] == 590 and s["strand"] == 1
                   for s in sites)
        # Linear template: no wrap, so the straddling site is absent.
        lin = sc._primer_binding_sites(primer, self.T, len(self.T),
                                       circular=False)
        assert all(s["foot_start"] != 590 for s in lin)

    def test_primer_longer_than_template(self):
        assert sc._primer_binding_sites("ACGT" * 50, "ACGTACGT", 8,
                                        circular=False) == []

    def test_foreign_char_raises(self):
        with pytest.raises(ValueError):
            sc._primer_binding_sites("ACGTZZZACGTACGT", self.T, len(self.T))

    def test_multiple_sites_repeat(self):
        # A primer whose binding region occurs twice → two forward sites.
        motif = "GACTAGCATGGATCCGTTACG"
        templ = ("AAAA" + motif + "TTTT" + motif + "CCCC")
        sites = sc._primer_binding_sites(motif, templ, len(templ),
                                         circular=False)
        starts = sorted(s["foot_start"] for s in sites if s["strand"] == 1)
        assert starts == [4, 4 + len(motif) + 4]


# ═══════════════════════════════════════════════════════════════════════════════
# Amplicon pairing — parity with _simulate_pcr on exact primers
# ═══════════════════════════════════════════════════════════════════════════════

class TestInsilicoPcrAmplicons:
    T = _rand_dna(600, seed=1)

    def _mine(self, p1, p2, *, circular, max_amplicon=20000):
        n = len(self.T)
        s1 = sc._primer_binding_sites(p1, self.T, n, circular=circular)
        s2 = sc._primer_binding_sites(p2, self.T, n, circular=circular)
        return sc._insilico_pcr_amplicons(s1, s2, n, circular=circular,
                                          max_amplicon=max_amplicon)

    def test_linear_parity_exact(self):
        fwd = self.T[100:122]
        rev = sc._rc(self.T[400:422])
        sim = sc._simulate_pcr(self.T, fwd, rev, circular=False)
        mine = self._mine(fwd, rev, circular=False)
        assert sim and mine
        assert mine[0]["start"] == sim[0]["start"] == 100
        assert mine[0]["length"] == sim[0]["length"] == 322

    def test_circular_parity_exact_nonwrap(self):
        fwd = self.T[100:122]
        rev = sc._rc(self.T[400:422])
        sim = sc._simulate_pcr(self.T, fwd, rev, circular=True)
        mine = self._mine(fwd, rev, circular=True)
        assert sim and mine
        assert mine[0]["start"] == sim[0]["start"] == 100
        assert mine[0]["length"] == sim[0]["length"] == 322
        assert mine[0]["wraps"] is False

    def test_circular_parity_exact_wrap(self):
        fwd = self.T[500:522]
        rev = sc._rc(self.T[50:72])
        sim = sc._simulate_pcr(self.T, fwd, rev, circular=True)
        mine = self._mine(fwd, rev, circular=True)
        assert sim and mine
        # Both place the amplicon at 500 spanning the origin.
        assert mine[0]["start"] == sim[0]["start"] == 500
        assert mine[0]["length"] == sim[0]["length"]
        assert mine[0]["wraps"] is True

    def test_both_tailed_length_parity(self):
        # When BOTH primers are tailed `_simulate_pcr`'s partial-binding
        # fallback fires; the product length (tails included) matches ours.
        fwd = "GGGGGCCCCC" + self.T[100:122]
        rev = "TTTTTAAAAA" + sc._rc(self.T[400:422])
        sim = sc._simulate_pcr(self.T, fwd, rev, circular=False)
        mine = self._mine(fwd, rev, circular=False)
        assert sim and mine
        assert mine[0]["length"] == sim[0]["length"]
        # The weaker primer caps certainty below a perfect match.
        assert mine[0]["certainty"] < 100.0

    def test_one_tailed_primer_still_found(self):
        # A pair where only ONE primer carries a 5' tail — a tailed cloning
        # primer against a stock sequencing primer, the commonest geometry
        # there is. Until the 2026-09-15 audit `_simulate_pcr` returned
        # NOTHING here (its partial-binding fallback required BOTH primers to
        # miss an exact match), so the two simulators disagreed on the same
        # pair; this test asserted that disagreement. They now agree.
        fwd = "GGGGGCCCCC" + self.T[100:122]      # tailed
        rev = sc._rc(self.T[400:422])             # exact
        mine = self._mine(fwd, rev, circular=False)
        assert mine
        expected_len = (400 + 22) - (122 - len(fwd))
        assert mine[0]["length"] == expected_len
        sim = sc._simulate_pcr(self.T, fwd, rev, circular=False)
        assert len(sim) == 1
        assert sim[0]["length"] == expected_len
        # The product carries the tail; the template arc under it is real.
        assert sim[0]["amplicon_seq"].startswith(fwd)
        assert sim[0]["amplicon_seq"].endswith(sc._rc(rev))

    def test_max_amplicon_filter(self):
        fwd = self.T[100:122]
        rev = sc._rc(self.T[400:422])
        # 322 bp product is filtered out by a 200 bp ceiling.
        assert self._mine(fwd, rev, circular=False, max_amplicon=200) == []
        assert self._mine(fwd, rev, circular=False, max_amplicon=400)

    def test_no_reverse_no_amplicon(self):
        fwd = self.T[100:122]
        # Two forward primers (no reverse-strand binder) → no product.
        fwd2 = self.T[300:322]
        s1 = sc._primer_binding_sites(fwd, self.T, len(self.T), circular=False)
        s2 = sc._primer_binding_sites(fwd2, self.T, len(self.T), circular=False)
        assert sc._insilico_pcr_amplicons(s1, s2, len(self.T),
                                          circular=False) == []

    def test_overlapping_primers_filtered_as_dimer(self):
        # Forward 5' at 100, reverse 3' at 105 → a 27 bp "product" shorter
        # than the two 22-mers laid end to end (44 bp) → filtered, matching
        # `_simulate_pcr`'s min_amp.
        fwd = self.T[100:122]
        rev = sc._rc(self.T[105:127])
        assert sc._simulate_pcr(self.T, fwd, rev, circular=False) == []
        assert self._mine(fwd, rev, circular=False) == []


# ═══════════════════════════════════════════════════════════════════════════════
# Feature mapping + confidence
# ═══════════════════════════════════════════════════════════════════════════════

class TestAmpliconFeatureSummary:
    T = _rand_dna(600, seed=1)

    def _rec(self, feats):
        rec = SeqRecord(Seq(self.T), id="x", name="x")
        for loc, label in feats:
            rec.features.append(SeqFeature(
                loc, type="CDS", qualifiers={"label": [label]}))
        return rec

    def test_picks_overlapping_feature(self):
        rec = self._rec([
            (FeatureLocation(150, 350, strand=1), "GeneX"),
            (FeatureLocation(500, 560, strand=1), "Elsewhere"),
        ])
        label, n = sc._amplicon_feature_summary(rec, 100, 322, len(self.T))
        assert label == "GeneX"
        assert n == 1

    def test_picks_largest_overlap(self):
        rec = self._rec([
            (FeatureLocation(150, 350, strand=1), "GeneX"),   # 200 bp overlap
            (FeatureLocation(105, 125, strand=1), "Tiny"),    # 20 bp overlap
        ])
        label, n = sc._amplicon_feature_summary(rec, 100, 322, len(self.T))
        assert label == "GeneX"
        assert n == 2

    def test_none_overlapping(self):
        rec = self._rec([(FeatureLocation(500, 560, strand=1), "Far")])
        label, n = sc._amplicon_feature_summary(rec, 100, 100, len(self.T))
        assert label == "—"
        assert n == 0

    def test_ignores_source_feature(self):
        rec = SeqRecord(Seq(self.T), id="x", name="x")
        rec.features.append(SeqFeature(
            FeatureLocation(0, len(self.T), strand=1), type="source"))
        rec.features.append(SeqFeature(
            FeatureLocation(150, 350, strand=1), type="CDS",
            qualifiers={"label": ["GeneX"]}))
        label, n = sc._amplicon_feature_summary(rec, 100, 322, len(self.T))
        assert label == "GeneX"
        assert n == 1

    def test_wrap_amplicon_overlaps_origin_feature(self):
        # Amplicon [560, 60) crosses the origin; a feature at [580, 600)
        # is covered.
        rec = self._rec([(FeatureLocation(580, 600, strand=1), "OriGene")])
        label, n = sc._amplicon_feature_summary(rec, 560, 100, len(self.T))
        assert label == "OriGene"
        assert n == 1


class TestPrimerCheckConfidence:
    @pytest.mark.parametrize("pct,glyph,color", [
        (100.0, "✓", "bright_cyan"),
        (95.0, "✓", "green"),
        (80.0, "⚠", "yellow"),
        (65.0, "~", "dark_orange"),
        (50.0, "✗", "red"),
    ])
    def test_tiers(self, pct, glyph, color):
        assert sc._primer_check_confidence(pct) == (glyph, color)

    def test_non_numeric(self):
        assert sc._primer_check_confidence(None) == ("?", "white")


# ═══════════════════════════════════════════════════════════════════════════════
# UI — the "Primer Check" tab
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrimerCheckTab:
    async def test_tab_and_widgets_present(self):
        from textual.widgets import Tab
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.PrimerDesignScreen("ACGT" * 200, [], "test"))
            await pilot.pause()
            tab_ids = {t.id for t in app.screen.query(Tab)}
            assert "tab-primercheck" in tab_ids
            for wid in ("#pd-pc-section", "#pd-pc-p1", "#pd-pc-p2",
                        "#pd-pc-scope", "#pd-pc-maxamp", "#btn-pd-pc-go",
                        "#pd-pc-results"):
                assert app.screen.query_one(wid) is not None

    async def test_switch_hides_design_chrome(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen._switch_mode("primercheck")
            await pilot.pause()
            assert screen.query_one("#pd-pc-section").display is True
            assert screen.query_one("#pd-template-section").display is False
            assert screen.query_one("#pd-results-section").display is False
            assert screen.query_one("#pd-bottom-actions").display is False
            assert screen._current_mode() == "primercheck"
            # Back to a design mode restores the chrome.
            screen._switch_mode("detection")
            await pilot.pause()
            assert screen.query_one("#pd-pc-section").display is False
            assert screen.query_one("#pd-template-section").display is True

    async def test_invalid_primer_shows_error(self):
        from textual.widgets import Input, Static
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.query_one("#pd-pc-p1", Input).value = "ACGTZZZACGTACGT"
            screen._pc_go()
            await pilot.pause()
            txt = str(screen.query_one("#pd-pc-status", Static).render())
            assert "Primer 1" in txt

    async def test_too_short_primer_shows_error(self):
        from textual.widgets import Input, Static
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.query_one("#pd-pc-p1", Input).value = "ACGT"
            screen._pc_go()
            await pilot.pause()
            txt = str(screen.query_one("#pd-pc-status", Static).render())
            assert "too short" in txt

    async def test_empty_primer1_hints_when_only_primer2_filled(self):
        from textual.widgets import Input, Static
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.query_one("#pd-pc-p2", Input).value = "ACGTACGTACGTACGTACGT"
            screen._pc_go()
            await pilot.pause()
            txt = str(screen.query_one("#pd-pc-status", Static).render())
            assert "Enter Primer 1" in txt
            assert "Primer 1 box" in txt

    async def test_scan_empty_library_reports_no_hits(self, monkeypatch):
        from textual.widgets import Input
        monkeypatch.setattr(sc, "_iter_library_readonly", lambda: [])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.query_one("#pd-pc-p1", Input).value = "ACGTACGTACGTACGTACGT"
            screen._pc_go()
            for _ in range(60):
                await pilot.pause()
                await asyncio.sleep(0.02)
                if getattr(screen, "_pc_rows", None) is not None:
                    break
            assert screen._pc_rows == []

    async def test_scan_skips_protein_and_malformed_entries(self, monkeypatch):
        from textual.widgets import Input
        T = _rand_dna(600, seed=1)
        fwd = T[100:122]

        def _dna_entry(eid, name, seq, kind="plasmid"):
            rec = SeqRecord(Seq(seq), id=eid, name=eid,
                            annotations={"molecule_type": "DNA",
                                         "topology": "circular"})
            return {"id": eid, "name": name, "kind": kind,
                    "topology": "circular", "size": len(seq),
                    "gb_text": sc._record_to_gb_text(rec)}

        good = _dna_entry("GOOD", "Good plasmid", T)
        # valid GenBank but flagged protein → skipped before the scan
        protein = _dna_entry("PROT", "Protein entry", T, kind="protein")
        # unparseable gb_text → per-entry guard skips it, scan continues
        malformed = {"id": "BAD", "name": "Bad entry", "kind": "plasmid",
                     "topology": "circular", "size": 50,
                     "gb_text": "NOT A GENBANK FILE"}
        monkeypatch.setattr(sc, "_iter_library_readonly",
                            lambda: [good, protein, malformed])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.query_one("#pd-pc-p1", Input).value = fwd
            screen._pc_go()
            for _ in range(80):
                await pilot.pause()
                await asyncio.sleep(0.03)
                if getattr(screen, "_pc_rows", None) is not None:
                    break
            rows = getattr(screen, "_pc_rows", None)
            assert rows is not None
            assert [r["name"] for r in rows] == ["Good plasmid"]

    async def test_scan_lists_hits_only(self, monkeypatch):
        from textual.widgets import Input
        T = _rand_dna(600, seed=1)
        T_miss = _rand_dna(600, seed=99)
        fwd = T[100:122]

        def _entry(eid, name, seq):
            rec = SeqRecord(Seq(seq), id=eid, name=eid,
                            annotations={"molecule_type": "DNA",
                                         "topology": "circular"})
            return {"id": eid, "name": name, "kind": "plasmid",
                    "topology": "circular", "size": len(seq),
                    "gb_text": sc._record_to_gb_text(rec)}

        entries = [_entry("HIT", "Hit plasmid", T),
                   _entry("MISS", "Miss plasmid", T_miss)]
        monkeypatch.setattr(sc, "_iter_library_readonly", lambda: entries)

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.query_one("#pd-pc-p1", Input).value = fwd
            screen._pc_go()
            for _ in range(80):
                await pilot.pause()
                await asyncio.sleep(0.03)
                if getattr(screen, "_pc_rows", None) is not None:
                    break
            rows = getattr(screen, "_pc_rows", None)
            assert rows is not None, "scan worker never produced results"
            names = [r["name"] for r in rows]
            assert "Hit plasmid" in names
            assert "Miss plasmid" not in names
            hit = next(r for r in rows if r["name"] == "Hit plasmid")
            assert hit["two"] is False
            assert hit["ident"] == 100.0
            assert hit["foot_start"] == 100

    async def test_ctrl_c_copies_selected_primer_seq(self, monkeypatch):
        from textual.widgets import DataTable
        captured = {}
        outcomes = []

        def _fake_copy(app, text, label="copy"):
            captured["text"] = text
            captured["label"] = label
            return ("clipboard", None)

        monkeypatch.setattr(sc, "_copy_to_clipboard_with_fallback", _fake_copy)
        monkeypatch.setattr(sc, "_load_primers", lambda: [
            {"name": "MyPrimer", "sequence": " acgtACGTacgtACGT ",
             "tm": 55.0, "primer_type": "detection", "status": "Designed"}])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            t = screen.query_one("#pd-lib-table", DataTable)
            screen.set_focus(t)
            t.move_cursor(row=0)
            await pilot.pause()
            monkeypatch.setattr(
                app, "_notify_copy_outcome",
                lambda n, what, mode, detail: outcomes.append(
                    (n, what, mode, detail)))
            screen.action_copy_primer_seq()
            await pilot.pause()
        # Sequence copied verbatim (normalised: trimmed + upper-cased).
        assert captured.get("text") == "ACGTACGTACGTACGT"
        assert "MyPrimer" in captured.get("label", "")
        # Toast reports it was copied + the base count.
        assert outcomes == [(16, "bases", "clipboard", None)]

    async def test_ctrl_c_skips_while_typing_in_an_input(self, monkeypatch):
        from textual.widgets import Input
        captured = {}
        monkeypatch.setattr(
            sc, "_copy_to_clipboard_with_fallback",
            lambda app, text, label="copy": captured.setdefault("text", text)
            or ("clipboard", None))
        monkeypatch.setattr(sc, "_load_primers", lambda: [
            {"name": "P", "sequence": "ACGTACGTACGTACGT"}])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.set_focus(screen.query_one("#pd-pc-p1", Input))
            await pilot.pause()
            screen.action_copy_primer_seq()
            await pilot.pause()
        assert "text" not in captured        # editing not hijacked

    async def test_ctrl_c_strips_control_chars_from_sequence(self, monkeypatch):
        from textual.widgets import DataTable
        captured = {}

        def _fake_copy(app, text, label="copy"):
            captured["text"] = text
            return ("clipboard", None)

        monkeypatch.setattr(sc, "_copy_to_clipboard_with_fallback", _fake_copy)
        # A crafted entry whose 'sequence' carries escape / control bytes,
        # whitespace and digits around the real bases (e.g. from a malicious
        # /primer_seq qualifier). Only IUPAC bases must reach the clipboard.
        monkeypatch.setattr(sc, "_load_primers", lambda: [
            {"name": "Crafted", "sequence": "acgt\x1b[2J\x07ACGT 123"}])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            t = screen.query_one("#pd-lib-table", DataTable)
            screen.set_focus(t)
            t.move_cursor(row=0)
            await pilot.pause()
            screen.action_copy_primer_seq()
            await pilot.pause()
        got = captured.get("text", "")
        assert got == "ACGTACGT"
        assert "\x1b" not in got and "\x07" not in got and " " not in got

    async def test_ctrl_c_no_primer_selected_is_safe(self, monkeypatch):
        monkeypatch.setattr(sc, "_load_primers", lambda: [])
        copied = {}
        monkeypatch.setattr(
            sc, "_copy_to_clipboard_with_fallback",
            lambda app, text, label="copy": copied.setdefault("text", text)
            or ("clipboard", None))
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.action_copy_primer_seq()    # empty library — must not crash
            await pilot.pause()
        assert "text" not in copied

    async def test_scan_two_primers_reports_amplicon(self, monkeypatch):
        from textual.widgets import Input
        T = _rand_dna(600, seed=1)
        fwd = T[100:122]
        rev = sc._rc(T[400:422])

        rec = SeqRecord(Seq(T), id="P", name="P",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec.features.append(SeqFeature(
            FeatureLocation(150, 350, strand=1), type="CDS",
            qualifiers={"label": ["GeneX"]}))
        entry = {"id": "P", "name": "Plasmid P", "kind": "plasmid",
                 "topology": "circular", "size": len(T),
                 "gb_text": sc._record_to_gb_text(rec)}
        monkeypatch.setattr(sc, "_iter_library_readonly", lambda: [entry])

        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen.query_one("#pd-pc-p1", Input).value = fwd
            screen.query_one("#pd-pc-p2", Input).value = rev
            screen._pc_go()
            for _ in range(80):
                await pilot.pause()
                await asyncio.sleep(0.03)
                if getattr(screen, "_pc_rows", None) is not None:
                    break
            rows = getattr(screen, "_pc_rows", None)
            assert rows, "two-primer scan produced no results"
            r = rows[0]
            assert r["two"] is True
            assert r["amp_len"] == 322
            assert r["feature"] == "GeneX"
            assert r["certainty"] == 100.0
