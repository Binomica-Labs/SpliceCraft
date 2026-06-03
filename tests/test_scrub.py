"""
test_scrub — clone-free restriction-site removal ("Scrub" tab).

Scrub cures a plasmid of chosen recognition sites by introducing the
minimal point substitutions that destroy each site while keeping every
overlapping CDS's protein intact and never spawning a new forbidden site.
The lab method is improved-QuikChange whole-plasmid PCR (DpnI + transform,
no ligase, no cloning) — see splicecraft `_scrub_design` / `_scrub_qc_primers`
and docs/invariants [INV-97].

These guard the biology contract. If they fail, DO NOT SHIP — a wrong
scrub silently changes a protein or leaves a "cured" site behind.

Asserted invariants (the real guarantees, deterministic regardless of
which synonymous codon the scorer happens to pick):

  * length unchanged + ACGT-only + substitution-only (only the reported
    edit positions differ from the original)
  * EVERY CDS feature's translation is byte-identical before/after
  * any site still present in the cured sequence is one we REPORTED as
    skipped (never silently left behind / forced)
  * every site we reported removed is actually gone
"""
from __future__ import annotations

import splicecraft as sc


# ── helpers ────────────────────────────────────────────────────────────────


def _cds(start, end, strand=1, label="orf", **extra):
    f = {"type": "CDS", "start": start, "end": end, "strand": strand,
         "label": label}
    f.update(extra)
    return f


def _site_idents(seq, enzymes, circular=True):
    """{(enzyme, rec_start, strand)} for every instance of `enzymes`."""
    forward = sc._scrub_resolve_sites(enzymes)
    allowed = frozenset(forward.keys())
    return {(t["enzyme"], t["rec_start"], t["strand"])
            for t in sc._scrub_scan_targets(seq, allowed, circular)}


def _assert_invariants(seq, feats, enzymes, *, circular=True):
    """Run `_scrub_design` and assert every biology guarantee. Returns the
    plan so callers can make extra case-specific assertions."""
    seq = seq.upper()
    plan = sc._scrub_design(seq, feats, enzymes, circular=circular)
    assert plan["ok"] is True
    cured = plan["cured_seq"]

    # length + alphabet + substitution-only
    assert len(cured) == len(seq), "scrub must not change sequence length"
    assert set(cured) <= set("ACGT"), "cured sequence must be ACGT-only"
    diffs = {i for i in range(len(seq)) if seq[i] != cured[i]}
    edit_positions = {e["pos"] for e in plan["edits"]}
    assert diffs == edit_positions, (
        "every base difference must correspond to a reported edit "
        f"(diffs={sorted(diffs)} edits={sorted(edit_positions)})")
    for e in plan["edits"]:
        assert e["frm"] == seq[e["pos"]]
        assert e["to"] == cured[e["pos"]]
        assert e["frm"] != e["to"]

    # every CDS protein preserved (across wrap / strand / frame)
    for f in feats or []:
        if f.get("type") != "CDS":
            continue
        before = sc._scrub_cds_protein(seq, f)
        after = sc._scrub_cds_protein(cured, f)
        assert before == after, (
            f"CDS {f.get('label')} protein changed: {before} -> {after}")

    # residual sites must all be reported as skipped; removed sites gone
    residual = _site_idents(cured, enzymes, circular)
    skipped = {(s["enzyme"], s["pos"], s["strand"])
               for s in plan["sites_skipped"]}
    assert residual <= skipped, (
        f"site left behind without being reported skipped: "
        f"{residual - skipped}")
    for r in plan["sites_removed"]:
        assert (r["enzyme"], r["pos"], r["strand"]) not in residual, (
            f"site reported removed but still present: {r}")
    return plan


# ── CDS site (the marquee case) ──────────────────────────────────────────────


class TestScrubCds:
    def test_bsai_in_cds_cured_silently(self):
        # ATG GGT CTC AAA GGG CCC TTT GAC TAA  → BsaI (GGTCTC) at nt 3,
        # spanning the Gly|Leu codon boundary. A synonymous wobble kills it.
        seq = "ATGGGTCTCAAAGGGCCCTTTGACTAA"
        feats = [_cds(0, 27, 1)]
        assert sc._translate_cds(seq, 0, 27, 1) == "MGLKGPFD*"
        plan = _assert_invariants(seq, feats, ["BsaI"])
        assert len(plan["edits"]) == 1
        assert plan["sites_removed"] and not plan["sites_skipped"]
        # the single edit lands inside the recognition window [3, 9)
        assert 3 <= plan["edits"][0]["pos"] < 9
        assert plan["edits"][0]["region"].startswith("CDS")
        assert plan["n_rounds"] == 1

    def test_reverse_strand_cds(self):
        # Same ORF read on the - strand. A reverse-strand BsaI lives here as
        # GAGACC on the top strand; the cure must respect the - frame.
        cds_fwd = "ATGGGTCTCAAAGGGCCCTTTGACTAA"
        seq = sc._rc(cds_fwd)            # now the ORF is on the - strand
        feats = [_cds(0, 27, -1)]
        _assert_invariants(seq, feats, ["BsaI"])

    def test_codon_frequency_tiebreak(self):
        # GGT|CTC (Gly|Leu) — the Gly wobble (nt 5) and a Leu swap (nt 8) both
        # kill BsaI silently. A table that strongly prefers GGA steers the
        # cure to GGT→GGA; without it the deterministic tie-break picks T→C.
        cds = "ATGGGTCTCAAAGGGCCCTTTGACTAA"
        feats = [_cds(0, 27, 1)]
        raw = {"GGA": ("G", 90), "GGC": ("G", 5), "GGG": ("G", 5),
               "GGT": ("G", 0),
               "CTA": ("L", 1), "CTG": ("L", 1), "CTC": ("L", 1),
               "CTT": ("L", 1), "TTA": ("L", 1), "TTG": ("L", 1)}
        plan = sc._scrub_design(cds, feats, ["BsaI"], codon_raw=raw)
        assert len(plan["edits"]) == 1
        assert (plan["edits"][0]["pos"], plan["edits"][0]["to"]) == (5, "A")
        assert plan["cured_seq"][3:6] == "GGA"
        assert (sc._translate_cds(plan["cured_seq"], 0, 27, 1)
                == sc._translate_cds(cds, 0, 27, 1))
        # no table → deterministic tie-break lands on a different synonym
        plan2 = sc._scrub_design(cds, feats, ["BsaI"])
        assert plan2["edits"][0]["to"] == "C"


# ── non-coding sites ─────────────────────────────────────────────────────────


class TestScrubNonCoding:
    def test_noncoding_forward(self):
        seq = "AAAACCCCGGTCTCAAAACCCCGGGGTT"   # one forward BsaI, no CDS
        plan = _assert_invariants(seq, [], ["BsaI"])
        assert len(plan["sites_removed"]) == 1
        assert plan["edits"][0]["region"] == "non-coding"

    def test_noncoding_reverse_strand(self):
        seq = "AAAAACCCCCGAGACCAAAAACCCCCG"    # GAGACC == reverse BsaI
        plan = _assert_invariants(seq, [], ["BsaI"])
        assert len(plan["sites_removed"]) == 1
        assert plan["sites_removed"][0]["strand"] == -1

    def test_noncoding_in_annotated_feature_still_cures(self):
        # A promoter (non-CDS) overlapping the site: no synonymy constraint,
        # but the region label notes it sits inside the annotation.
        seq = "AAAACCCCGGTCTCAAAACCCCGGGGTT"
        feats = [{"type": "promoter", "start": 0, "end": 20, "label": "Pxyz"}]
        plan = _assert_invariants(seq, feats, ["BsaI"])
        assert len(plan["sites_removed"]) == 1
        assert "Pxyz" in plan["edits"][0]["region"]


# ── origin-spanning site ─────────────────────────────────────────────────────


class TestScrubWrap:
    def test_site_across_origin(self):
        # GGT at the 3' end + CTC at the 5' end → GGTCTC wraps the origin.
        seq = "CTCAAAAACCCCCAAAAACCCCCGGT"     # len 26
        assert _site_idents(seq, ["BsaI"]), "test setup: expected a wrap hit"
        plan = _assert_invariants(seq, [], ["BsaI"])
        assert len(plan["sites_removed"]) == 1
        # the cured base is one of the wrap recognition positions
        assert plan["edits"][0]["pos"] in {23, 24, 25, 0, 1, 2}


# ── two overlapping CDSes on opposite strands ────────────────────────────────


class TestScrubDualFrame:
    def test_overlapping_cds_both_frames_preserved(self):
        # Same span annotated as a CDS on BOTH strands. Whatever the scrubber
        # does, NEITHER protein may change — and if no change satisfies both
        # frames the site must be reported skipped, never forced.
        seq = "ATGGGTCTCAAAGGGCCCTTTAAATAG"   # forward ORF w/ BsaI
        feats = [_cds(0, 27, 1, label="fwd"), _cds(0, 27, -1, label="rev")]
        plan = _assert_invariants(seq, feats, ["BsaI"])
        # accounting: the one BsaI site is either removed or explicitly skipped
        assert len(plan["sites_removed"]) + len(plan["sites_skipped"]) == 1


# ── multi-enzyme: no new forbidden site of ANY scrubbed enzyme ────────────────


class TestScrubMultiEnzyme:
    def test_adjacent_bsai_and_esp3i(self):
        # GGTCTC (BsaI) and CGTCTC (Esp3I) close together. Curing one must not
        # spawn the other (or a BbsI) anywhere.
        seq = "AAAAGGTCTCAAAACGTCTCAAAAGGGGTT"
        enzymes = ["BsaI", "Esp3I", "BbsI"]
        before = _site_idents(seq, enzymes)
        assert len(before) == 2
        plan = _assert_invariants(seq, [], enzymes)
        assert len(plan["sites_removed"]) == 2
        assert not _site_idents(plan["cured_seq"], enzymes)


# ── no-op / empty paths ──────────────────────────────────────────────────────


class TestScrubEdges:
    def test_no_sites_present(self):
        seq = "AAAACCCCGGGGTTTTAAAACCCCGGGGTT"
        plan = sc._scrub_design(seq, [], ["BsaI"])
        assert plan["ok"] and plan["cured_seq"] == seq
        assert plan["edits"] == [] and plan["sites_removed"] == []
        assert plan["n_rounds"] == 0

    def test_empty_enzyme_set_warns(self):
        seq = "AAAAGGTCTCAAAA"
        plan = sc._scrub_design(seq, [], [])
        assert plan["edits"] == []
        assert any("enzyme" in w.lower() for w in plan["warnings"])

    def test_empty_sequence(self):
        plan = sc._scrub_design("", [], ["BsaI"])
        assert plan["ok"] and plan["cured_seq"] == ""
        assert plan["edits"] == []

    def test_large_plasmid_warns(self):
        # >8 kb earns the linear-amplification warning. Build clean filler.
        seq = ("ACGTACGT" * 1100)          # 8800 bp, no GGTCTC
        seq = seq.replace("GGTCTC", "GGAACC")
        plan = sc._scrub_design(seq, [], ["BsaI"])
        assert any("8 kb" in w or "linear" in w for w in plan["warnings"])


# ── clustering into QuikChange rounds ────────────────────────────────────────


class TestScrubClustering:
    def test_near_edits_one_round(self):
        assert sc._scrub_cluster_edits([5, 7, 9], 1000) == [[5, 7, 9]]

    def test_distant_edits_separate_rounds(self):
        clusters = sc._scrub_cluster_edits([5, 70], 1000)
        assert clusters == [[5], [70]]

    def test_origin_adjacent_edits_merge(self):
        # 98 and 2 are 4 bp apart across the origin of a 100 bp plasmid.
        clusters = sc._scrub_cluster_edits([2, 98], 100)
        assert len(clusters) == 1
        assert set(clusters[0]) == {2, 98}

    def test_design_reports_round_count(self):
        # Two BsaI sites far apart on BOTH arcs of the circle (>30 bp each
        # way) → two separate QuikChange rounds. (Sites close across the
        # origin would correctly merge into one round — see the cluster
        # unit tests above.)
        seq = ("AAAAGGTCTCAAAA"            # BsaI @4
               + "C" * 40
               + "GGTCTCAAAA"             # BsaI @54
               + "C" * 40)                # arcs ≈ 50 and 54 bp
        plan = sc._scrub_design(seq, [], ["BsaI"])
        assert len(plan["sites_removed"]) == 2
        assert plan["n_rounds"] == 2


# ── improved-QuikChange primer design ────────────────────────────────────────

# Balanced 50%-GC filler (like a real plasmid) that provably contains no
# default-set site: "GATCATGC" has no "GG"/"AA"/"AG"/"CG→T", so no GGTCTC,
# GAGACC, CGTCTC or GAAGAC can form (checked by the assert in _design).
# Inserting exactly one GGTCTC gives a clean single-site template.
_FILLER = ("GATCATGC" * 33)                  # 264 bp, 50% GC


def _one_site_template(at: int) -> str:
    """A ~264 bp template with a single BsaI site starting at `at`."""
    return (_FILLER[:at] + "GGTCTC" + _FILLER[at + 6:])[:len(_FILLER)]


class TestScrubQCPrimers:
    def _design(self, at, overlap="improved"):
        seq = _one_site_template(at)
        assert len(_site_idents(seq, ["BsaI"])) == 1
        plan = sc._scrub_design(seq, [], ["BsaI"])
        cured = plan["cured_seq"]
        assert plan["n_rounds"] == 1
        positions = plan["clusters"][0]["positions"]
        qc = sc._scrub_qc_primers(cured, positions, overlap=overlap)
        return seq, cured, positions, qc

    def test_binding_equals_display(self):
        # The reported primer sequence MUST equal the circular slice of the
        # cured template at the reported coordinates (catastrophic-class).
        seq, cured, positions, qc = self._design(120)
        assert "error" not in qc, qc.get("error")
        n = len(cured)
        assert qc["fwd_seq"] == sc._circ_extract(
            cured, qc["fwd_start"], qc["fwd_len"], n)
        assert qc["rev_seq"] == sc._mut_revcomp(sc._circ_extract(
            cured, qc["rev_start"], qc["rev_len"], n))

    def test_primer_geometry(self):
        seq, cured, positions, qc = self._design(120)
        n = len(cured)
        # lengths in range
        assert 25 <= qc["fwd_len"] <= 48 and 25 <= qc["rev_len"] <= 48
        # the cure sits inside BOTH primer footprints
        for g in positions:
            assert any((qc["fwd_start"] + i) % n == g
                       for i in range(qc["fwd_len"]))
            assert any((qc["rev_start"] + i) % n == g
                       for i in range(qc["rev_len"]))
        # carries the cure as an internal mismatch vs the parent
        assert qc["n_mismatch"] >= 1
        # improved = partial overlap → the two footprints are NOT identical
        assert (qc["fwd_start"], qc["fwd_len"]) != (qc["rev_start"], qc["rev_len"])
        assert qc["overlap_style"] == "improved"
        assert qc["fwd_tm"] and qc["rev_tm"]
        assert qc["fwd_tm_qc"] and qc["rev_tm_qc"]

    def test_primer_actually_carries_cure(self):
        # The forward primer, aligned to the PARENT at its coords, differs
        # exactly at the cured base(s) — proof it encodes the cure.
        seq, cured, positions, qc = self._design(120)
        n = len(cured)
        parent_footprint = sc._circ_extract(seq, qc["fwd_start"], qc["fwd_len"], n)
        mismatches = sum(1 for a, b in zip(qc["fwd_seq"], parent_footprint)
                         if a != b)
        assert mismatches == qc["n_mismatch"] >= 1

    def test_origin_adjacent_primer_wraps(self):
        # Site near the very start → primers straddle the origin. Binding
        # must still equal display through the wrap.
        seq, cured, positions, qc = self._design(2)
        assert "error" not in qc, qc.get("error")
        n = len(cured)
        assert qc["fwd_seq"] == sc._circ_extract(
            cured, qc["fwd_start"], qc["fwd_len"], n)
        assert qc["rev_seq"] == sc._mut_revcomp(sc._circ_extract(
            cured, qc["rev_start"], qc["rev_len"], n))

    def test_classic_full_overlap_option(self):
        seq, cured, positions, qc = self._design(120, overlap="classic")
        assert "error" not in qc, qc.get("error")
        assert qc["overlap_style"] == "classic"
        # classic = full overlap → fwd footprint == rev footprint
        assert (qc["fwd_start"], qc["fwd_len"]) == (qc["rev_start"], qc["rev_len"])

    def test_template_too_small_errors(self):
        qc = sc._scrub_qc_primers("ACGT" * 5, [4], overlap="improved")
        assert "error" in qc

    def test_cluster_span_wrap(self):
        # circular arc helper: positions either side of the origin
        assert sc._scrub_cluster_span([10, 12], 100) == (10, 12)
        s, e = sc._scrub_cluster_span([2, 98], 100)
        assert (s, e) == (98, 2)        # wraps: end < start


# ── MutagenizeModal "Scrub" tab (UI wiring) ──────────────────────────────────

_BASELINE = (160, 48)


def _record_with_bsai(n=120, at=60):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    filler = ("GATCATGC" * ((n // 8) + 2))[:n]
    seq = (filler[:at] + "GGTCTC" + filler[at + 6:])[:n]
    rec = SeqRecord(Seq(seq), id="scrubtest", name="scrubtest",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    return seq, rec


class TestScrubTabStructure:
    async def test_modal_exposes_both_tabs_and_scrub_widgets(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app.push_screen(sc.MutagenizeModal("ATGAAAGGG" * 4, [], "p"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # both tabs present
            assert modal.query_one("#mut-tab-soe") is not None
            assert modal.query_one("#mut-tab-scrub") is not None
            # scrub widgets present
            for wid in ("#btn-scrub-run", "#scrub-overlap", "#scrub-results",
                        "#scrub-results-body", "#btn-scrub-apply",
                        "#btn-scrub-saveprimers", "#btn-scrub-tomap",
                        "#btn-scrub-enzymes", "#btn-scrub-close",
                        "#btn-scrub-codon", "#scrub-codon-label"):
                assert modal.query_one(wid) is not None
            # Copy protocol was removed — it must no longer be in the tree.
            assert len(modal.query("#btn-scrub-protocol")) == 0
            # commit buttons start disabled (no scrub computed yet)
            for wid in ("#btn-scrub-apply", "#btn-scrub-saveprimers",
                        "#btn-scrub-tomap"):
                assert modal.query_one(wid, sc.Button).disabled is True
            # SOE tab still intact
            assert modal.query_one("#btn-mut-design", sc.Button) is not None
            assert modal.query_one("#mut-source", sc.Select) is not None


async def _assert_scrub_layout_fits(size):
    app = sc.PlasmidApp()
    async with app.run_test(size=size) as pilot:
        await pilot.pause()
        app.push_screen(sc.MutagenizeModal("ATG" * 200, [], "p",
                                           show_tab="scrub"))
        await pilot.pause()
        await pilot.pause(0.1)
        modal = app.screen
        box = modal.query_one("#mut-box").region
        # fullscreen: the box fills the whole terminal
        assert (box.x, box.y) == (0, 0), f"modal not at origin @ {size}: {box}"
        assert box.width == size[0] and box.height == size[1], (
            f"modal not fullscreen @ {size}: {box}")
        assert modal.query_one("#mut-title").region.y >= 0
        for wid in ("#btn-scrub-codon", "#btn-scrub-enzymes",
                    "#btn-scrub-run", "#btn-scrub-apply",
                    "#btn-scrub-saveprimers", "#btn-scrub-tomap",
                    "#btn-scrub-close"):
            r = modal.query_one(wid).region
            assert r.width > 0 and r.height > 0, f"{wid} not laid out @ {size}"
            # fully inside the box — horizontally AND vertically (the action
            # row was being clipped under the box's bottom border).
            assert box.x <= r.x and r.right <= box.right, (
                f"{wid} overflows box horizontally @ {size}: {r} vs {box}")
            assert box.y <= r.y and r.bottom <= box.bottom, (
                f"{wid} overflows box vertically @ {size}: {r} vs {box}")
        # bottom action buttons share the row evenly (1fr each)
        widths = [modal.query_one(w).region.width for w in
                  ("#btn-scrub-apply", "#btn-scrub-saveprimers",
                   "#btn-scrub-tomap", "#btn-scrub-close")]
        assert max(widths) - min(widths) <= 2, f"uneven buttons: {widths}"


class TestScrubTabLayout:
    # Regression for the 2026-06-02 snapshots: the modal first clipped its
    # title off the top, then the action-button row got clipped under the
    # box's bottom border. Check the whole thing fits — title, every Scrub
    # button (horizontally AND vertically) — at a tall and a short terminal.
    async def test_fits_tall_terminal(self):
        await _assert_scrub_layout_fits((171, 43))

    async def test_fits_short_terminal(self):
        await _assert_scrub_layout_fits((120, 30))


class TestScrubTabRender:
    def test_render_without_mount(self):
        # _render_scrub only touches its args, so an un-mounted shell
        # exercises it deterministically.
        seq, _rec = _record_with_bsai()
        plan = sc._scrub_design(seq, [], ["BsaI"])
        rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                       round_no=i)
                  for i, c in enumerate(plan["clusters"], 1)]
        m = sc.MutagenizeModal.__new__(sc.MutagenizeModal)
        rendered = m._render_scrub(plan, rounds).plain
        assert "Removed 1 site" in rendered
        assert "BsaI" in rendered
        assert "FWD" in rendered and "REV" in rendered
        # the PCR → DpnI → transform protocol summary is still shown in-report
        assert "DpnI" in rendered


class TestScrubTabApply:
    async def test_apply_to_canvas_then_undo(self):
        seq, rec = _record_with_bsai()
        rec._tui_display_name = "Scrub Test 1"   # spaces — must survive
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app._apply_snapshot(seq, 0, rec)         # load test plasmid
            await pilot.pause()
            feats = app.query_one("#plasmid-map", sc.PlasmidMap)._feats
            app.push_screen(sc.MutagenizeModal(seq, feats, "scrubtest"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            # compute the plan + primers and feed the UI callback directly
            # (avoids racing the @work thread in the headless harness)
            plan = sc._scrub_design(seq, feats, ["BsaI"], circular=True)
            rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                           round_no=i)
                      for i, c in enumerate(plan["clusters"], 1)]
            modal._scrub_apply_result(plan, rounds)
            await pilot.pause()
            assert modal.query_one("#btn-scrub-apply", sc.Button).disabled is False

            modal._scrub_apply_canvas(None)
            await pilot.pause()
            cured = str(app._current_record.seq)
            assert cured == plan["cured_seq"]
            assert cured != seq
            assert not sc._scrub_scan_targets(cured, frozenset(["BsaI"]), True)
            # the user-typed display name (with spaces) must survive the cure
            assert getattr(app._current_record, "_tui_display_name",
                           None) == "Scrub Test 1"
            # apply button re-disabled after a successful apply
            assert modal.query_one("#btn-scrub-apply", sc.Button).disabled is True

            # the apply is undoable
            app._action_undo()
            await pilot.pause()
            assert str(app._current_record.seq) == seq


# ── agent endpoint: scrub-plasmid ────────────────────────────────────────────


class TestScrubEndpoint:
    def test_seq_path_cures(self):
        seq, _rec = _record_with_bsai()
        res = sc._h_scrub_plasmid(None, {"seq": seq, "enzymes": ["BsaI"]})
        assert res["ok"] is True
        assert res["cured_seq"] and res["cured_seq"] != seq
        assert len(res["sites_removed"]) == 1
        assert not sc._scrub_scan_targets(res["cured_seq"],
                                          frozenset(["BsaI"]), True)
        assert res["rounds"] and "fwd_seq" in res["rounds"][0]

    def test_missing_seq_no_record_400(self):
        res = sc._h_scrub_plasmid(None, {})
        assert isinstance(res, tuple) and res[1] == 400

    def test_bad_overlap_400(self):
        res = sc._h_scrub_plasmid(None, {"seq": "ACGT" * 30, "overlap": "x"})
        assert isinstance(res, tuple) and res[1] == 400

    def test_features_not_list_400(self):
        res = sc._h_scrub_plasmid(None, {"seq": "ACGT" * 30, "features": "no"})
        assert isinstance(res, tuple) and res[1] == 400

    def test_codon_taxid_bad_type_400(self):
        res = sc._h_scrub_plasmid(None, {"seq": "ACGT" * 30, "codon_taxid": 83333})
        assert isinstance(res, tuple) and res[1] == 400

    def test_codon_taxid_unknown_404(self):
        res = sc._h_scrub_plasmid(
            None, {"seq": "ACGT" * 30, "codon_taxid": "999999999"})
        assert isinstance(res, tuple) and res[1] == 404

    def test_canvas_path_preserves_protein(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        cds = "ATGGGTCTCAAAGGGCCCTTTGACTAA"     # BsaI at nt 3, in-frame
        rec = SeqRecord(Seq(cds), id="t", name="t",
                        annotations={"molecule_type": "DNA"})
        rec.features.append(SeqFeature(FeatureLocation(0, 27, strand=1),
                                       type="CDS",
                                       qualifiers={"label": ["orf"]}))
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        app._current_record = rec
        res = sc._h_scrub_plasmid(app, {"enzymes": ["BsaI"]})
        assert res["ok"] is True
        cured = res["cured_seq"]
        # protein preserved (the canvas path read the CDS frame)
        assert (sc._translate_cds(cured, 0, 27, 1)
                == sc._translate_cds(cds, 0, 27, 1))
        assert not sc._scrub_scan_targets(cured, frozenset(["BsaI"]), True)
        assert res["edits"][0]["region"].startswith("CDS")


# ── primer output scroll + "Primers → Map" ───────────────────────────────────


def _multi_site_seq():
    # 600 bp, 50% GC, four well-separated BsaI sites → four QuikChange rounds
    # → a primer report tall enough to need scrolling.
    s = list(("GATCATGC" * 80)[:600])
    for at in (50, 200, 350, 500):
        s[at:at + 6] = list("GGTCTC")
    return "".join(s)


class TestScrubScroll:
    async def test_long_primer_report_scrolls(self):
        seq = _multi_site_seq()
        assert len(_site_idents(seq, ["BsaI"])) == 4
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app.push_screen(sc.MutagenizeModal(seq, [], "multi",
                                               show_tab="scrub"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            plan = sc._scrub_design(seq, [], ["BsaI"], circular=True)
            assert plan["n_rounds"] == 4
            rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                           round_no=i)
                      for i, c in enumerate(plan["clusters"], 1)]
            modal._scrub_apply_result(plan, rounds)
            await pilot.pause()
            results = modal.query_one("#scrub-results")
            # a real scroll container (not a bare clipping Static)
            assert isinstance(results, sc.VerticalScroll)
            assert results.max_scroll_y > 0, "tall report should be scrollable"


class TestScrubToMap:
    async def test_primers_added_to_map(self):
        seq, rec = _record_with_bsai()
        rec._tui_display_name = "Scrub Test 1"   # spaces — must survive save
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app._apply_snapshot(seq, 0, rec)
            await pilot.pause()
            feats = app.query_one("#plasmid-map", sc.PlasmidMap)._feats
            app.push_screen(sc.MutagenizeModal(seq, feats, "scrubtest",
                                               show_tab="scrub"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            plan = sc._scrub_design(seq, feats, ["BsaI"], circular=True)
            rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                           round_no=i)
                      for i, c in enumerate(plan["clusters"], 1)]
            modal._scrub_apply_result(plan, rounds)
            await pilot.pause()
            assert modal.query_one("#btn-scrub-tomap", sc.Button).disabled is False
            before = sum(1 for f in app._current_record.features
                         if f.type == "primer_bind")
            modal._scrub_save_to_map(None)
            await pilot.pause()
            scrub_primers = [
                f for f in app._current_record.features
                if f.type == "primer_bind"
                and (f.qualifiers.get("label", [""]) or [""])[0].startswith(
                    "SCRUB_")]
            # one round → a fwd + a rev primer_bind feature, each with seq
            assert len(scrub_primers) >= 2
            assert all("primer_seq" in f.qualifiers for f in scrub_primers)
            after = sum(1 for f in app._current_record.features
                        if f.type == "primer_bind")
            assert after >= before + 2
            # the user-typed display name (with spaces) must survive the
            # primer-add — a freshly-rebuilt record would drop it and the
            # next save would underscore it ("FFE 6" → "FFE_6").
            assert getattr(app._current_record, "_tui_display_name",
                           None) == "Scrub Test 1"
            # undoable
            app._action_undo()
            await pilot.pause()
            assert sum(1 for f in app._current_record.features
                       if f.type == "primer_bind") == before
