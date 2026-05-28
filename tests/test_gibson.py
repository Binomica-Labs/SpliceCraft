"""
test_gibson — Gibson assembly simulator + UI smoke.

Covers the module-level helpers that power ConstructorModal's "Gibson"
tab: `_gibson_overlap_len`, `_simulate_gibson_assembly`,
`_gibson_record_from_result`. The UI is exercised in a small async
smoke test against `GibsonAssemblyPane`.

Pure unit tests where possible — no Textual app required for the
simulator. The simulator is a deterministic transform on
fragment-dict lists; tests assert sequence equality, feature
shifting, and the failure-mode signalling that the UI surfaces.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


# ──────────────────────────────────────────────────────────────────────────────
# _gibson_overlap_len
# ──────────────────────────────────────────────────────────────────────────────


class TestGibsonOverlapLen:
    def test_exact_overlap_of_20_bp(self):
        a = "TTT" + "ACGTACGTACGTACGTACGT"   # 20 bp tail
        b = "ACGTACGTACGTACGTACGT" + "GGG"   # same 20 bp head
        assert sc._gibson_overlap_len(a, b, min_overlap=15) == 20

    def test_no_overlap_returns_zero(self):
        assert sc._gibson_overlap_len("AAAA" * 10, "CCCC" * 10,
                                        min_overlap=15) == 0

    def test_overlap_below_minimum_returns_zero(self):
        # 10 bp homology — below default min_overlap of 15.
        a = "TTT" + "ACGTACGTAC"
        b = "ACGTACGTAC" + "GGG"
        assert sc._gibson_overlap_len(a, b, min_overlap=15) == 0
        assert sc._gibson_overlap_len(a, b, min_overlap=10) == 10

    def test_prefers_longest_match(self):
        # 30 bp homology embedded in 50 bp tail — must pick the longest,
        # not stop at the first 15 bp.
        homology = "ACGT" * 10  # 40 bp
        a = "TTT" + homology
        b = homology + "GGG"
        assert sc._gibson_overlap_len(a, b, min_overlap=15) == 40

    def test_case_insensitive(self):
        assert sc._gibson_overlap_len(
            "ttt" + "acgtacgtacgtacgtacgt",
            "ACGTACGTACGTACGTACGT" + "GGG",
            min_overlap=15,
        ) == 20

    def test_empty_inputs(self):
        assert sc._gibson_overlap_len("", "ACGT", min_overlap=15) == 0
        assert sc._gibson_overlap_len("ACGT", "", min_overlap=15) == 0


# ──────────────────────────────────────────────────────────────────────────────
# _simulate_gibson_assembly — happy paths
# ──────────────────────────────────────────────────────────────────────────────


def _frag(name: str, sequence: str, features=None) -> dict:
    return {"name": name, "sequence": sequence,
            "features": list(features or [])}


class TestSimulateGibsonAssemblyHappyPath:
    def test_two_fragments_linear(self):
        oh = "ACGTACGTACGTACGTACGT"  # 20 bp
        f1 = _frag("F1", "AAAAA" + oh)
        f2 = _frag("F2", oh + "TTTTT")
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        # Linear: F1 + F2[overlap:]
        assert r["product_seq"] == "AAAAA" + oh + "TTTTT"
        assert len(r["overlaps"]) == 1
        assert r["overlaps"][0]["length"] == 20
        assert r["overlaps"][0]["ok"] is True

    def test_three_fragments_circular(self):
        oh_a = "AAAAAAAAAAAAAAAAAAAA"
        oh_b = "CCCCCCCCCCCCCCCCCCCC"
        oh_c = "GGGGGGGGGGGGGGGGGGGG"
        # F1's 3' overlaps with F2's 5' (oh_a).
        # F2's 3' overlaps with F3's 5' (oh_b).
        # F3's 3' overlaps with F1's 5' (oh_c) — wraps the circle.
        f1 = _frag("F1", oh_c + "TTTTTTTTTT" + oh_a)
        f2 = _frag("F2", oh_a + "TTTTTTTTTT" + oh_b)
        f3 = _frag("F3", oh_b + "TTTTTTTTTT" + oh_c)
        r = sc._simulate_gibson_assembly([f1, f2, f3], min_overlap=15,
                                            circular=True)
        assert r["success"] is True
        # Each fragment is 50 bp; pairwise overlap of 20 bp; circular
        # dedup drops the wrap overlap. Expected total: 90 bp.
        expected = (oh_c + "TTTTTTTTTT" + oh_a + "TTTTTTTTTT"
                     + oh_b + "TTTTTTTTTT")
        assert r["product_seq"] == expected
        assert len(r["product_seq"]) == 90
        assert len(r["overlaps"]) == 3
        assert all(o["ok"] for o in r["overlaps"])
        # Last junction is the wrap.
        assert r["overlaps"][-1]["is_wrap"] is True
        assert r["overlaps"][0]["is_wrap"] is False

    def test_single_fragment_linear_passthrough(self):
        f = _frag("solo", "AAAAACCCCCGGGGGTTTTT")
        r = sc._simulate_gibson_assembly([f], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        assert r["product_seq"] == "AAAAACCCCCGGGGGTTTTT"
        assert r["overlaps"] == []
        assert any("Single linear fragment" in w
                    for w in r["warnings"])


# ──────────────────────────────────────────────────────────────────────────────
# _simulate_gibson_assembly — feature shifting
# ──────────────────────────────────────────────────────────────────────────────


class TestSimulateGibsonAssemblyFeatures:
    def test_feature_in_first_fragment_keeps_position(self):
        oh = "ACGT" * 5  # 20 bp
        f1 = _frag("F1", "GGGGGGGGGG" + oh,  # 30 bp
                    features=[{"start": 0, "end": 10, "strand": 1,
                                "type": "CDS", "label": "gene1"}])
        f2 = _frag("F2", oh + "TTTTTTTTTT")
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        # Feature should be at the same position in the product.
        feats = [f for f in r["features"] if f.get("label") == "gene1"]
        assert len(feats) == 1
        assert feats[0]["start"] == 0
        assert feats[0]["end"] == 10

    def test_feature_in_second_fragment_shifts(self):
        oh = "ACGT" * 5  # 20 bp
        f1 = _frag("F1", "AAAAA" + oh)  # 25 bp
        f2 = _frag("F2", oh + "TTTTTGGGGG",  # 30 bp
                    features=[{"start": 25, "end": 30, "strand": 1,
                                "type": "CDS", "label": "gene2"}])
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        # F2's position 25 shifts to product position
        #   offset = len(F1) - oh = 25 - 20 = 5
        # so gene2 lives at [5+25, 5+30) = [30, 35)
        # Product length = 5 + 20 + 10 = 35
        assert len(r["product_seq"]) == 35
        feats = [f for f in r["features"] if f.get("label") == "gene2"]
        assert len(feats) == 1
        assert feats[0]["start"] == 30
        assert feats[0]["end"] == 35

    def test_feature_in_leading_overlap_of_later_fragment_is_skipped(self):
        # A feature entirely inside F2's leading-overlap region would be
        # a duplicate of one F1 already supplies. The simulator drops
        # it to avoid double-annotation.
        oh = "ACGT" * 5  # 20 bp
        f1 = _frag("F1", "AAAAA" + oh)
        f2 = _frag("F2", oh + "TTTTTTTTTT",
                    features=[{"start": 0, "end": 20, "strand": 1,
                                "type": "misc_feature",
                                "label": "homology"}])
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        labels = [f.get("label") for f in r["features"]]
        assert "homology" not in labels

    def test_feature_straddling_overlap_is_kept(self):
        # A feature that straddles the overlap boundary (start in
        # overlap, end past it) is kept — the user designed something
        # meaningful that crosses the junction.
        oh = "ACGT" * 5  # 20 bp
        f1 = _frag("F1", "AAAAA" + oh)
        f2 = _frag("F2", oh + "TTTTTTTTTT",
                    features=[{"start": 10, "end": 25, "strand": 1,
                                "type": "CDS", "label": "spanning"}])
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        feats = [f for f in r["features"] if f.get("label") == "spanning"]
        assert len(feats) == 1
        # offset = 25 - 20 = 5, so span shifts to [15, 30)
        assert feats[0]["start"] == 15
        assert feats[0]["end"] == 30


# ──────────────────────────────────────────────────────────────────────────────
# _simulate_gibson_assembly — failure modes
# ──────────────────────────────────────────────────────────────────────────────


class TestSimulateGibsonAssemblyFailures:
    def test_no_overlap_fails(self):
        f1 = _frag("F1", "A" * 50)
        f2 = _frag("F2", "C" * 50)
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is False
        assert r["errors"]
        assert len(r["overlaps"]) == 1
        assert r["overlaps"][0]["ok"] is False
        assert r["overlaps"][0]["length"] == 0

    def test_circular_wrap_failure(self):
        # Linear chain works (F1→F2, F2→F3 share overlap) but F3's tail
        # doesn't match F1's head — circular closure should fail and
        # the wrap-junction error names the wrap.
        oh = "ACGT" * 5  # 20 bp
        f1 = _frag("F1", "AAA" + oh)
        f2 = _frag("F2", oh + "BBB" + oh)
        f3 = _frag("F3", oh + "CCC")
        r = sc._simulate_gibson_assembly([f1, f2, f3], min_overlap=15,
                                            circular=True)
        assert r["success"] is False
        # Junction 3 is the wrap. The other two pass.
        oks = [o["ok"] for o in r["overlaps"]]
        assert oks == [True, True, False]
        assert r["overlaps"][-1]["is_wrap"] is True

    def test_empty_fragments_list(self):
        r = sc._simulate_gibson_assembly([], min_overlap=15,
                                            circular=True)
        assert r["success"] is False
        assert r["errors"]

    def test_fragment_fully_consumed_by_overlaps(self):
        # Edge case: F3 is exactly the wrap overlap of F1 + the leading
        # overlap with F2 — its body would be empty. Should refuse.
        oh = "ACGTACGTACGTACGTACGT"  # 20 bp
        f1 = _frag("F1", oh + "AAA" + oh)
        f2 = _frag("F2", oh + "BBB" + oh)
        f3 = _frag("F3", oh + oh)   # 40 bp = lead + trail overlap
        r = sc._simulate_gibson_assembly([f1, f2, f3], min_overlap=15,
                                            circular=True)
        assert r["success"] is False
        assert any("fully consumed" in e for e in r["errors"])


# ──────────────────────────────────────────────────────────────────────────────
# _gibson_record_from_result — SeqRecord construction
# ──────────────────────────────────────────────────────────────────────────────


class TestGibsonRecordFromResult:
    def test_circular_product_to_record(self):
        oh = "ACGT" * 5
        f1 = _frag("F1", oh + "AAAAAAAAAA" + oh)
        f2 = _frag("F2", oh + "TTTTTTTTTT" + oh)
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=True)
        assert r["success"]
        rec = sc._gibson_record_from_result(r, name="myGibson")
        assert rec is not None
        assert str(rec.seq) == r["product_seq"]
        assert rec.annotations["topology"] == "circular"
        assert rec.annotations["molecule_type"] == "DNA"

    def test_linear_product_to_record(self):
        oh = "ACGT" * 5
        f1 = _frag("F1", "AAAAA" + oh)
        f2 = _frag("F2", oh + "TTTTT")
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        rec = sc._gibson_record_from_result(r, name="lin")
        assert rec is not None
        assert rec.annotations["topology"] == "linear"

    def test_failed_result_returns_none(self):
        r = {"success": False, "product_seq": "", "circular": True,
             "features": [], "overlaps": [], "errors": ["nope"],
             "warnings": []}
        assert sc._gibson_record_from_result(r, name="x") is None

    def test_features_round_trip_through_genbank(self):
        oh = "ACGT" * 5
        f1 = _frag("F1", oh + "AAAAAAAAAA" + oh,
                    features=[{"start": 20, "end": 30, "strand": 1,
                                "type": "CDS", "label": "g1"}])
        f2 = _frag("F2", oh + "TTTTTTTTTT" + oh)
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=True)
        rec = sc._gibson_record_from_result(r, name="rt")
        assert rec is not None
        # Round-trip through GenBank text + back.
        gb_text = sc._record_to_gb_text(rec)
        rec2 = sc._gb_text_to_record(gb_text)
        # The g1 label must survive.
        labels = [(f.qualifiers.get("label") or [""])[0]
                   for f in rec2.features
                   if (f.qualifiers.get("label") or [""])[0]]
        assert "g1" in labels


# ──────────────────────────────────────────────────────────────────────────────
# GibsonAssemblyPane — UI smoke (full Textual app)
# ──────────────────────────────────────────────────────────────────────────────


class TestGibsonAssemblyPaneSmoke:
    @pytest.mark.asyncio
    async def test_constructor_modal_has_gibson_tab(self, tiny_record,
                                                       isolated_library):
        """Smoke test: open ConstructorModal, switch to the Gibson tab,
        confirm the pane mounts without errors and the lane table is
        empty by default."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            # Pane exists.
            pane = modal.query_one("#ctor-gib-pane",
                                      sc.GibsonAssemblyPane)
            assert pane is not None
            # Lane starts empty.
            assert pane._lane == []
            # Topology defaults to circular.
            assert pane._is_circular() is True
            # Min-overlap defaults to the global.
            assert pane._min_overlap() == sc._GIBSON_MIN_OVERLAP_BP

    @pytest.mark.asyncio
    async def test_pane_add_paste_fragment_and_simulate(
            self, tiny_record, isolated_library):
        """Drive the pane via direct method calls (faster + more
        deterministic than clicking through the UI). Adds two paste
        fragments with a 20 bp overlap and triggers Simulate."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            pane = modal.query_one("#ctor-gib-pane",
                                      sc.GibsonAssemblyPane)
            # Push two fragments into the lane directly.
            oh = "ACGTACGTACGTACGTACGT"  # 20 bp
            pane._lane.append({
                "name": "f1", "sequence": "AAAAA" + oh,
                "features": [], "source": "test",
            })
            pane._lane.append({
                "name": "f2", "sequence": oh + "TTTTT",
                "features": [], "source": "test",
            })
            # Run the validation pipeline.
            pane._refresh_lane_table()
            pane._refresh_overlap_view()
            # Topology defaults to circular but circular here requires
            # F2's tail to overlap F1's head — these fragments only have
            # one valid overlap (F1's tail → F2's head). Switch to linear
            # via the radio: easier to assert success.
            pane.query_one("#gib-topo-linear", sc.RadioButton).value = True
            await pilot.pause()
            # Trigger simulate by calling the handler with a dummy event.
            pane._on_simulate(None)
            await pilot.pause()
            assert pane._product is not None
            assert pane._product["success"] is True
            # Save button should be enabled.
            save_btn = pane.query_one("#btn-gib-save", sc.Button)
            assert save_btn.disabled is False


# ──────────────────────────────────────────────────────────────────────────────
# Edge cases — sanitisation, self-circularisation, wrap features
# ──────────────────────────────────────────────────────────────────────────────


class TestGibsonEdgeCases:
    # Non-repetitive 20 bp homology arm — distinct enough that the only
    # self-overlap within the fragment is the intended arm match.
    ARM = "AAACCCGGGTTTAAACCCGGG"[:20]   # 20 bp, no internal repeat

    def test_self_circularisation_n1_circular(self):
        """A single fragment whose 3' end overlaps its own 5' end
        forms a circle. Body = fragment − wrap overlap."""
        arm = self.ARM
        body_only = "TTTTTGGGGGAAAAA"   # 15 bp, no match to arm
        f = _frag("self", arm + body_only + arm)
        r = sc._simulate_gibson_assembly([f], min_overlap=15,
                                            circular=True)
        assert r["success"] is True
        # Product = arm + body_only (35 bp; wrap overlap dropped once)
        assert r["product_seq"] == arm + body_only
        assert len(r["overlaps"]) == 1
        assert r["overlaps"][0]["is_wrap"] is True
        assert r["overlaps"][0]["length"] == 20

    def test_self_circularisation_no_overlap_fails(self):
        # No self-overlap of 15+ bp: a 30 bp run of A's followed by 30
        # bp of C's. 3' end is C-only; 5' end is A-only.
        f = _frag("noself", "A" * 30 + "C" * 30)
        r = sc._simulate_gibson_assembly([f], min_overlap=15,
                                            circular=True)
        assert r["success"] is False
        assert r["errors"]

    def test_self_circularisation_overlap_eq_full_length_capped(self):
        """A fragment whose entire body matches its prefix would
        trivially produce a self-match equal to fragment length. The
        algorithm caps at len-1 so the user gets a real overlap (or
        a junction failure), not a degenerate "fully consumed" pass
        followed by an empty product."""
        arm = self.ARM    # 20 bp
        # Fragment is the arm by itself — len == 20, min_overlap == 15.
        # Cap forces max_check = 19, so we scan 19..15 for an exact
        # suffix==prefix match inside the arm. The arm has no such
        # match (non-repetitive), so the junction fails.
        f = _frag("ouroboros", arm)
        r = sc._simulate_gibson_assembly([f], min_overlap=15,
                                            circular=True)
        assert r["success"] is False
        assert r["errors"]

    def test_middle_fragment_fully_consumed(self):
        """A middle fragment whose entire sequence equals the leading
        overlap leaves no body to contribute. Caught by per-fragment
        body validation; surfaces a clear 'consumed' error."""
        arm = self.ARM   # 20 bp, non-repetitive
        f1 = _frag("F1", "GAAGAA" + arm)   # 26 bp
        f2 = _frag("F2", arm)              # 20 bp = arm exactly
        f3 = _frag("F3", arm + "TGTGTG")   # 26 bp
        r = sc._simulate_gibson_assembly([f1, f2, f3], min_overlap=15,
                                            circular=False)
        assert r["success"] is False
        # F2's leading overlap == its full length → "consumed" error.
        assert any("consumed" in e for e in r["errors"])

    def test_whitespace_and_lowercase_normalised(self):
        oh_lower = "acgtacgtacgtacgtacgt"
        oh_upper = "ACGTACGTACGTACGTACGT"
        f1 = _frag("F1", "aaa\n" + oh_lower)
        f2 = _frag("F2", oh_upper + "\tTTT")
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        # Output is always uppercase, whitespace stripped.
        assert r["product_seq"] == "AAA" + oh_upper + "TTT"
        assert "\n" not in r["product_seq"]
        assert "\t" not in r["product_seq"]

    def test_empty_sequence_rejected(self):
        f1 = _frag("F1", "ACGT" * 10)
        f_empty = _frag("F2", "")
        r = sc._simulate_gibson_assembly([f1, f_empty], min_overlap=15,
                                            circular=False)
        assert r["success"] is False
        assert any("no sequence" in e for e in r["errors"])

    def test_non_dict_fragment_rejected(self):
        r = sc._simulate_gibson_assembly(
            [{"name": "F1", "sequence": "ACGT" * 20}, "not a dict"],
            min_overlap=15, circular=False,
        )
        assert r["success"] is False
        assert any("dict" in e for e in r["errors"])

    def test_short_fragment_warning(self):
        """A fragment shorter than 3 × min_overlap surfaces a warning
        but doesn't block the assembly."""
        oh = "ACGTACGTACGTACGTACGT"  # 20 bp
        f1 = _frag("F1", "A" * 50 + oh)
        f2 = _frag("F2", oh + "TTTTT")   # 25 bp = short (< 3*15)
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        assert any("short" in w for w in r["warnings"])

    def test_wrap_feature_round_trip_in_product(self):
        """A feature in the last fragment that straddles the trailing
        overlap lands as a wrap feature (end < start) in the circular
        product — the bases at F_last's tail equal F1's head, so the
        feature crosses the wrap junction."""
        oh = "ACGTACGTACGTACGTACGT"        # 20 bp
        body = "TTTTTTTTTTGGGGGGGGGG"      # 20 bp
        # F1: oh + 10 A's + oh = 50 bp
        # F2: oh + body + oh   = 60 bp (trailing oh wraps to F1's head)
        f1 = _frag("F1", oh + "AAAAAAAAAA" + oh)
        # Feature spans F2 positions [30, 50): half is body's last 10 bp,
        # half is the trailing overlap (which dedups onto F1's head).
        f2 = _frag("F2", oh + body + oh,
                    features=[{"start": 30, "end": 50, "strand": 1,
                                "type": "CDS", "label": "wraps"}])
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=True)
        assert r["success"] is True
        prod_len = len(r["product_seq"])
        # Product = F1 (50) + body (20) = 70 bp.
        assert prod_len == 70
        feats = [f for f in r["features"] if f.get("label") == "wraps"]
        assert len(feats) == 1
        s, e = feats[0]["start"], feats[0]["end"]
        # Wrap encoding: end < start signals an origin-spanning feature.
        # Half lands at product[60:70] (body[10:20] = 'G'*10), the other
        # half at product[0:10] (F1[:10] = first half of `oh`).
        assert s == 60 and e == 10

    def test_overlap_capped_at_max(self):
        """Overlap detection caps at `max_overlap` even when a longer
        exact match exists — bounds the worst-case probe so a
        pathological repetitive sequence can't blow up runtime."""
        # 300 bp arm followed by distinct bodies. Without the cap the
        # whole arm could match; with max_overlap=200, capped to 200.
        arm = "ACGTACGTAC" * 30   # 300 bp; non-self-identical caller
        a = "AAA" + arm           # 303 bp, arm at the 3' end
        b = arm + "TTT"           # 303 bp, arm at the 5' end
        k = sc._gibson_overlap_len(a, b, min_overlap=15,
                                     max_overlap=200)
        assert k == 200
        # With max raised, the full arm (300 bp) is detected.
        k = sc._gibson_overlap_len(a, b, min_overlap=15,
                                     max_overlap=400)
        assert k == 300

    def test_self_overlap_capped_below_full_length(self):
        """The self-match cap (`len(a) - 1` when `a == b`) skips the
        trivial whole-string match and recovers the real homology arm
        embedded at fragment ends."""
        arm = self.ARM   # 20 bp non-repetitive
        body = "TTTTTTTTTTTTTTTTTTTT"  # 20 bp, no overlap with arm
        seq = arm + body + arm   # 60 bp
        # Self-probe: cap at 59 prevents the trivial 60-bp full-string
        # match; the real arm (20 bp) is found instead.
        k = sc._gibson_overlap_len(seq, seq, min_overlap=15)
        assert k == 20

    def test_min_overlap_clamped_to_one(self):
        a = "AAAACGT"
        b = "CGT" + "TTTT"   # 3 bp overlap
        # min_overlap negative — clamps to 1 internally.
        assert sc._gibson_overlap_len(a, b, min_overlap=-5) == 3
        assert sc._gibson_overlap_len(a, b, min_overlap=0) == 3

    @pytest.mark.asyncio
    async def test_pane_n1_self_circ_flow(self, tiny_record,
                                              isolated_library):
        """User adds a single self-circ fragment + circular topology
        → simulator detects self-overlap, save button enables."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            pane = modal.query_one("#ctor-gib-pane",
                                      sc.GibsonAssemblyPane)
            arm = "AAACCCGGGTTTAAACCCGG"   # 20 bp non-repetitive
            body = "GAGTAGTCATCATCAGTGT"   # 19 bp, no arm match
            # Single self-circ fragment.
            pane._lane.append({
                "name":     "selfcirc",
                "sequence": arm + body + arm,
                "features": [],
                "source":   "test",
            })
            pane._refresh_lane_table()
            # Default topology = circular, default min_overlap = 15.
            pane._on_simulate(None)
            await pilot.pause()
            assert pane._product is not None
            assert pane._product["success"] is True
            assert len(pane._product["product_seq"]) == 20 + 19   # arm + body
            assert pane._product["overlaps"][0]["is_wrap"] is True
            save_btn = pane.query_one("#btn-gib-save", sc.Button)
            assert save_btn.disabled is False

    @pytest.mark.asyncio
    async def test_pane_paste_size_cap_rejects(self, tiny_record,
                                                   isolated_library):
        """Pasting a sequence longer than _GIB_MAX_PASTE_BP is rejected
        with a user-facing error rather than silently pegging the
        textarea."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            pane = modal.query_one("#ctor-gib-pane",
                                      sc.GibsonAssemblyPane)
            pane._mode = "pcr"
            pane._apply_mode_visibility()
            # 1 MB + 1 bp of 'A' — exceeds the cap by 1.
            too_long = "A" * (pane._GIB_MAX_PASTE_BP + 1)
            pane.query_one("#gib-pcr-name", sc.Input).value = "huge"
            pane.query_one("#gib-pcr-seq", sc.TextArea).text = too_long
            await pilot.pause()
            with pytest.raises(ValueError) as exc:
                pane._build_fragment_from_inputs()
            assert "exceeds" in str(exc.value)

    def test_overlaps_list_populated_on_failure(self):
        """Failure path still exposes the overlaps list — UI uses it to
        render which junctions passed vs failed."""
        f1 = _frag("F1", "AAAA" + "ACGTACGTACGTACGTACGT")
        f2 = _frag("F2", "ACGTACGTACGTACGTACGT" + "TTTT")  # 20 bp ok
        f3 = _frag("F3", "X" * 50)                          # F2→F3 fails
        r = sc._simulate_gibson_assembly([f1, f2, f3], min_overlap=15,
                                            circular=False)
        assert r["success"] is False
        # Even on failure, every junction is reported.
        assert len(r["overlaps"]) == 2
        assert r["overlaps"][0]["ok"] is True
        assert r["overlaps"][1]["ok"] is False


# ──────────────────────────────────────────────────────────────────────────────
# Hardening regressions (2026-05-14)
# ──────────────────────────────────────────────────────────────────────────────


class TestGibsonHardening:
    """Regressions for the 2026-05-14 hardening pass:

    * H4 RC probe — flipped fragments surface a "did you mean to RC"
      hint instead of silently failing with "no homology".
    * H1 wrap math — refactored wrap-coord logic doesn't drop features
      whose span equals product length.
    * H3 wrap-pair merge — wrap-source halves re-merge into one wrap
      feature in the product when their halves still span a junction.
    * H2 negative-offset skip — features that fall before the product
      start (pathological middle-fragment exhaustion) are skipped, not
      silently clamped to 0.
    * Cleanup — internal ``_wrap_*`` sentinel fields never leak into
      product feature dicts.
    """

    ARM = "AAACCCGGGTTTAAACCCGGG"[:20]

    def test_rc_hint_when_second_fragment_flipped(self):
        """Junction failure with an RC'd second fragment surfaces a
        targeted hint naming the flipped fragment."""
        arm = self.ARM
        f1 = _frag("F1", "AAAAA" + arm)
        f2_correct_seq = arm + "TTTTT"
        # Flipped: user provided the RC of what they intended.
        f2_flipped = _frag("F2_flipped", sc._rc(f2_correct_seq))
        r = sc._simulate_gibson_assembly([f1, f2_flipped], min_overlap=15,
                                            circular=False)
        assert r["success"] is False
        assert r["overlaps"][0]["ok"] is False
        # rc_hint must mention F2_flipped specifically.
        hint = r["overlaps"][0]["rc_hint"]
        assert "F2_flipped" in hint
        # User-facing error inherits the hint.
        assert any("F2_flipped" in e for e in r["errors"])
        assert any("flip" in e.lower() for e in r["errors"])

    def test_rc_hint_when_first_fragment_flipped(self):
        """When the upstream fragment is flipped, the hint targets it
        (not the downstream)."""
        arm = self.ARM
        f1_correct_seq = "AAAAA" + arm
        f1_flipped = _frag("F1_flipped", sc._rc(f1_correct_seq))
        f2 = _frag("F2", arm + "TTTTT")
        r = sc._simulate_gibson_assembly([f1_flipped, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is False
        hint = r["overlaps"][0]["rc_hint"]
        assert "F1_flipped" in hint

    def test_no_rc_hint_when_genuinely_no_overlap(self):
        """When neither orientation matches, no misleading hint is
        appended."""
        f1 = _frag("F1", "A" * 50)
        f2 = _frag("F2", "C" * 50)
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is False
        assert r["overlaps"][0]["rc_hint"] == ""

    def test_wrap_sentinels_stripped_from_output(self):
        """The internal `_wrap_pair`, `_wrap_role`, `_wrap_total`
        marker fields used by the simulator to detect wrap-source
        halves never leak into output feature dicts."""
        oh = "ACGT" * 5
        # Wrap-source feature: end < start, source is whole-plasmid
        # so it's split into two halves by `_record_features`. We
        # construct the halved input directly here to stay independent
        # of GibsonAssemblyPane.
        f1 = _frag("F1", oh + "AAAAAAAAAA" + oh, features=[
            {"start": 25, "end": 30, "strand": 1, "type": "CDS",
             "label": "wrap_tail", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "tail",
             "_wrap_total": 30},
            {"start": 0, "end": 5, "strand": 1, "type": "CDS",
             "label": "wrap_tail", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "head",
             "_wrap_total": 30},
        ])
        f2 = _frag("F2", oh + "TTTTTTTTTT" + oh)
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=True)
        assert r["success"] is True
        for f in r["features"]:
            for key in f.keys():
                assert not key.startswith("_wrap_"), (
                    f"internal sentinel field {key!r} leaked into "
                    f"output: {f!r}"
                )

    def test_wrap_pair_halves_both_survive_when_split_in_product(self):
        """When a wrap-source feature's halves come from a fragment
        that's interior to the product (not spanning the product's
        own wrap junction), both halves survive as separate features
        in the product. The simulator must NOT silently drop one of
        the halves — both annotate biologically distinct bases of
        the source.

        Real-Gibson merge of the halves into one wrap feature would
        require the source plasmid to span the entire product
        without homology-arm trimming, which the chemistry doesn't
        produce; the merge code path in `_simulate_gibson_assembly`
        is defensive for future scenarios.
        """
        oh = "ACGT" * 5    # 20 bp
        f1_seq = oh + "AAAAAAAAAA" + oh   # 50 bp
        f1 = _frag("F1", f1_seq, features=[
            {"start": 40, "end": 50, "strand": 1, "type": "CDS",
             "label": "circular_feature", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "tail",
             "_wrap_total": 50},
            {"start": 0, "end": 5, "strand": 1, "type": "CDS",
             "label": "circular_feature", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "head",
             "_wrap_total": 50},
        ])
        f2 = _frag("F2", oh + "TTTTTTTTTT" + oh)
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=True)
        assert r["success"] is True
        feats = [f for f in r["features"]
                  if f.get("label") == "circular_feature"]
        # Both halves preserved (no merge — halves don't span product wrap).
        assert len(feats) == 2
        # Sentinel keys stripped from output.
        for f in feats:
            assert not any(k.startswith("_wrap_") for k in f)
        # Tail at [40, 50), head at [0, 5) — coords preserved through
        # the shift (offset 0 for F1).
        starts = sorted(f["start"] for f in feats)
        ends   = sorted(f["end"]   for f in feats)
        assert starts == [0, 40]
        assert ends   == [5, 50]

    def test_wrap_pair_remains_split_when_halves_separated(self):
        """When the two halves of a wrap-source feature end up at
        non-junction positions in the product (e.g. the source plasmid
        was used as a non-first fragment, so the halves no longer
        share the product wrap), both pieces remain as separate
        features."""
        oh = "ACGT" * 5
        # F1 is plain. F2 carries the wrap halves at its own internal
        # positions — they'll shift by F2's offset and no longer span
        # the product wrap.
        f1 = _frag("F1", "AAAAA" + oh + "AAAAA")   # 30 bp
        f2 = _frag("F2", oh + "GGGGGGGGGG" + oh, features=[
            # Halves at F2-local positions [25, 30) and [0, 5).
            {"start": 25, "end": 30, "strand": 1, "type": "CDS",
             "label": "circular_feature", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "tail",
             "_wrap_total": 30},
            {"start": 0, "end": 5, "strand": 1, "type": "CDS",
             "label": "circular_feature", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "head",
             "_wrap_total": 30},
        ])
        # Linear so we can reason about positions cleanly. F2 lands
        # at offset = len(F1) - oh = 30 - 20 = 10 (after lead-strip).
        # Wait — F1 doesn't end with `oh`; need to give it a tail
        # overlap with F2.
        f1 = _frag("F1", "AAAAA" + oh)              # 25 bp, tail oh
        # F2 leads with oh, then body, then trailing arm not relevant.
        f2 = _frag("F2", oh + "GGGGGGGGGG" + "TTT", features=[
            # Halves at F2-local positions [30, 33) and [0, 3).
            {"start": 30, "end": 33, "strand": 1, "type": "CDS",
             "label": "split_feature", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "tail",
             "_wrap_total": 33},
            {"start": 0, "end": 3, "strand": 1, "type": "CDS",
             "label": "split_feature", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "head",
             "_wrap_total": 33},
        ])
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        feats = [f for f in r["features"]
                  if f.get("label") == "split_feature"]
        # The head half at F2-local [0, 3) falls inside F2's leading
        # overlap (20 bp); skip filter drops it. Only the tail half
        # survives — 1 feature, not merged.
        assert len(feats) == 1
        # No sentinel keys leaked.
        assert not any(k.startswith("_wrap_") for k in feats[0])

    def test_wrap_pair_head_inside_leading_overlap_filtered(self):
        """When a wrap-source's head-half is entirely inside the
        leading-overlap region of a downstream fragment, the skip
        filter drops it (the preceding fragment supplies the same
        bases). The tail half passes through normally."""
        oh = "ACGT" * 5   # 20 bp
        f1 = _frag("F1", "AAAAA" + oh)
        f2 = _frag("F2", oh + "TTTTT" + oh, features=[
            # Head half is [0, 10) — entirely inside the 20 bp leading
            # overlap. Tail half is [20, 30) — past the overlap.
            {"start": 20, "end": 30, "strand": 1, "type": "CDS",
             "label": "split", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "tail",
             "_wrap_total": 30},
            {"start": 0, "end": 10, "strand": 1, "type": "CDS",
             "label": "split", "color": "",
             "_wrap_pair": "wrap:0", "_wrap_role": "head",
             "_wrap_total": 30},
        ])
        r = sc._simulate_gibson_assembly([f1, f2], min_overlap=15,
                                            circular=False)
        assert r["success"] is True
        feats = [f for f in r["features"] if f.get("label") == "split"]
        # Head dropped (inside leading overlap); tail kept.
        assert len(feats) == 1

    def test_record_features_marks_wrap_pair(self):
        """`GibsonAssemblyPane._record_features` tags both halves of a
        wrap-source feature with the same ``_wrap_pair`` id so the
        simulator can re-merge them later."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, CompoundLocation, FeatureLocation
        rec = SeqRecord(Seq("ACGT" * 10), id="wrap", name="wrap")
        rec.annotations["topology"] = "circular"
        rec.annotations["molecule_type"] = "DNA"
        # CompoundLocation wrap: [35, 40) + [0, 5)
        loc = CompoundLocation([
            FeatureLocation(35, 40, strand=1),
            FeatureLocation(0, 5, strand=1),
        ])
        rec.features.append(SeqFeature(loc, type="CDS",
                                       qualifiers={"label": ["w"]}))
        out = sc.GibsonAssemblyPane._record_features(rec)
        # Two halves, same _wrap_pair id, different roles.
        wrap_halves = [f for f in out if f.get("label") == "w"]
        assert len(wrap_halves) == 2
        pair_ids = {f["_wrap_pair"] for f in wrap_halves}
        assert len(pair_ids) == 1
        roles = {f["_wrap_role"] for f in wrap_halves}
        assert roles == {"head", "tail"}

    @pytest.mark.asyncio
    async def test_save_dispatches_worker(self, tiny_record,
                                              isolated_library):
        """Clicking Save opens the name+collection prompt; confirming it
        dispatches the background worker with the captured snapshot +
        the chosen name/collection (instead of writing synchronously on
        the UI thread). We confirm by patching ``_gibson_save_worker``."""
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            pane = modal.query_one("#ctor-gib-pane",
                                      sc.GibsonAssemblyPane)
            # Stage a single self-circ fragment so Simulate produces a
            # committed product the Save button can act on.
            arm = "AAACCCGGGTTTAAACCCGG"
            body = "GAGTAGTCATCATCAGTGT"
            pane._lane.append({
                "name":     "selfcirc",
                "sequence": arm + body + arm,
                "features": [],
                "source":   "test",
            })
            pane._refresh_lane_table()
            pane._on_simulate(None)
            await pilot.pause()
            assert pane._product is not None

            # Intercept the worker — assert it's the dispatch point.
            calls: list[dict] = []

            def fake_worker(*, product, lane, circular, entry_counter,
                            name=None, collection=None):
                calls.append({
                    "product_success": product.get("success"),
                    "lane_len": len(lane),
                    "circular": circular,
                    "entry_counter": entry_counter,
                    "name": name,
                    "collection": collection,
                })

            pane._gibson_save_worker = fake_worker  # type: ignore[assignment]
            pane._on_save(None)
            await pilot.pause()
            # Save now opens the name+collection prompt; the worker is NOT
            # dispatched until the user confirms.
            name_modal = app.screen
            assert isinstance(name_modal, sc.NamePlasmidModal)
            name_modal.query_one("#nameplasmid-input", sc.Input).value = \
                "My Gibson Build"
            await pilot.pause()
            name_modal._try_submit()
            await pilot.pause()
            assert len(calls) == 1
            assert calls[0]["product_success"] is True
            assert calls[0]["lane_len"] == 1
            assert calls[0]["name"] == "My Gibson Build"
            assert isinstance(calls[0]["collection"], str)
            # Save button disabled on confirm to prevent double-click
            # duplicate inserts (the worker hasn't finished).
            save_btn = pane.query_one("#btn-gib-save", sc.Button)
            assert save_btn.disabled is True
