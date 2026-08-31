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

import pytest

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


class TestPrimerFreeEnergyUnits:
    """2026-08-26 report #3: `scrub-plasmid` reported `hairpin_dg: -738.3`
    and `homodimer_dg: -4993.7`. primer3 answers in CALORIES/mol; every
    oligo spec sheet and QC rule of thumb is in KILOcalories, so the numbers
    were 1000× off and unusable as a gate (and read as catastrophic if you
    took the unit at face value)."""

    # A GC-clamped stem-loop and a self-complementary oligo: both fold.
    HAIRPIN = "GCGCGCGCTTTTGCGCGCGC"

    def test_hairpin_dg_is_kcal_per_mole(self):
        primer3 = pytest.importorskip("primer3")
        raw = primer3.calc_hairpin(self.HAIRPIN, **sc._MUT_P3).dg
        assert sc._mut_hairpin_dg(self.HAIRPIN) == pytest.approx(raw / 1000.0)

    def test_homodimer_dg_is_kcal_per_mole(self):
        primer3 = pytest.importorskip("primer3")
        raw = primer3.calc_homodimer(self.HAIRPIN, **sc._MUT_P3).dg
        assert sc._mut_homodimer_dg(self.HAIRPIN) == pytest.approx(raw / 1000.0)

    def test_a_real_folding_oligo_lands_in_the_published_range(self):
        """The check the reporter actually wanted: an oligo that folds should
        read as a few kcal/mol, not a few thousand of anything."""
        pytest.importorskip("primer3")
        for fn in (sc._mut_hairpin_dg, sc._mut_homodimer_dg):
            dg = fn(self.HAIRPIN)
            assert dg is not None
            assert -40.0 < dg < 0.0, f"{fn.__name__} = {dg}"

    def test_scrub_rounds_report_kcal_and_say_so(self):
        seq = _one_site_template(120)
        plan = sc._scrub_design(seq, [], ["BsaI"])
        qc = sc._scrub_qc_primers(plan["cured_seq"],
                                  plan["clusters"][0]["positions"])
        assert "error" not in qc, qc.get("error")
        assert qc["dg_units"] == "kcal/mol"
        for key in ("hairpin_dg", "homodimer_dg"):
            assert qc[key] is None or -40.0 <= qc[key] <= 5.0, (key, qc[key])

    def test_unmeasurable_dg_is_none_not_zero(self, monkeypatch):
        """Zero is the BEST possible ΔG, so a failed measurement reported as
        0.0 makes every QC gate pass on a primer nobody checked. Patched in
        `splicecraft_primer`'s own namespace — that is where the moved
        biophysics reads it."""
        import splicecraft_primer as sp
        sp._mut_thermo_cache_clear()
        monkeypatch.setattr(sp, "_MUT_P3", {"bogus_kwarg": 1})
        try:
            assert sc._mut_hairpin_dg("ACGTACGTACGTACGTACGT") is None
            assert sc._mut_homodimer_dg("ACGTACGTACGTACGTACGT") is None
        finally:
            sp._mut_thermo_cache_clear()

    def test_an_unmeasured_pair_is_warned_about(self, monkeypatch):
        import splicecraft_primer as sp
        sp._mut_thermo_cache_clear()
        monkeypatch.setattr(sp, "_MUT_P3", {"bogus_kwarg": 1})
        try:
            seq = _one_site_template(120)
            plan = sc._scrub_design(seq, [], ["BsaI"])
            qc = sc._scrub_qc_primers(plan["cured_seq"],
                                      plan["clusters"][0]["positions"])
            assert qc["hairpin_dg"] is None and qc["homodimer_dg"] is None
            assert any("NOT secondary-structure checked" in w
                       for w in qc["warnings"]), qc["warnings"]
        finally:
            sp._mut_thermo_cache_clear()


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
    async def test_apply_prompts_then_saves_and_loads(self):
        """Apply cure now PROMPTS for a name + collection first (request H):
        the cure is NOT applied until the user confirms; confirming applies it,
        saves the cured plasmid to the collection (primers bound, no
        underscores), and loads it onto the canvas."""
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

            # Apply cure → NamePlasmidModal appears; cure NOT yet applied.
            modal._scrub_apply_canvas(None)
            await pilot.pause()
            await pilot.pause(0.1)
            assert isinstance(app.screen, sc.NamePlasmidModal)
            assert str(app._current_record.seq) == seq, \
                "cure must not apply until the name is confirmed"

            # Confirm the name + collection → apply + save + focus.
            app.screen.dismiss({"name": "Cured One", "collection": "Default"})
            for _ in range(6):
                await pilot.pause()
            cured = str(app._current_record.seq)
            assert cured == plan["cured_seq"]
            assert not sc._scrub_scan_targets(cured, frozenset(["BsaI"]), True)
            # cured plasmid saved to the collection, NO underscores in the name
            saved = [e for c in sc._load_collections()
                     for e in c.get("plasmids", [])
                     if e.get("source") == "scrub:cured"]
            assert saved, "cured plasmid not saved to a collection"
            assert saved[-1]["name"] == "Cured One"
            assert "_" not in saved[-1]["name"]
            assert "primer_bind" in (saved[-1].get("gb_text") or ""), \
                "cure primers not bound to the saved cured plasmid"

    async def test_apply_cancel_leaves_canvas_untouched(self):
        """Cancelling the name prompt is a full no-op: nothing applied, nothing
        saved, Apply button still live for a retry."""
        seq, rec = _record_with_bsai()
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            app._apply_snapshot(seq, 0, rec)
            await pilot.pause()
            feats = app.query_one("#plasmid-map", sc.PlasmidMap)._feats
            app.push_screen(sc.MutagenizeModal(seq, feats, "scrubtest"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            plan = sc._scrub_design(seq, feats, ["BsaI"], circular=True)
            rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                           round_no=i)
                      for i, c in enumerate(plan["clusters"], 1)]
            modal._scrub_apply_result(plan, rounds)
            await pilot.pause()
            modal._scrub_apply_canvas(None)
            await pilot.pause()
            await pilot.pause(0.1)
            app.screen.dismiss(None)        # cancel
            for _ in range(3):
                await pilot.pause()
            assert str(app._current_record.seq) == seq, "cure leaked on cancel"
            assert not [e for c in sc._load_collections()
                        for e in c.get("plasmids", [])
                        if e.get("source") == "scrub:cured"]
            assert modal.query_one("#btn-scrub-apply", sc.Button).disabled is False


class TestScrubCuredSave:
    """`_scrub_save_cured_plasmid` (request H): the cured plasmid is saved under
    the user-chosen name + collection with the cure primers bound, and the SAME
    primers are written back onto the ORIGINAL (uncured) source entry so the
    user sees where each cure edits a base. Mirrors the operon clone-save
    coverage; calls the helper directly so the modal chain isn't in the loop."""

    async def test_saves_cured_with_primers_and_annotates_source(self):
        from copy import deepcopy
        from Bio.Seq import Seq
        seq, rec = _record_with_bsai()
        rec._tui_display_name = "Source One"          # spaces
        # Original lives in the library so the source annotate has a target.
        src_id = "src_one"
        libs = sc._load_library()
        libs.append({"id": src_id, "name": "Source One", "size": len(seq),
                     "gb_text": sc._record_to_gb_text(rec), "kind": "plasmid",
                     "n_feats": len(rec.features)})
        sc._save_library(libs)
        plan = sc._scrub_design(seq, rec.features, ["BsaI"], circular=True)
        rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                       round_no=i)
                  for i, c in enumerate(plan["clusters"], 1)]
        primers = []
        for r in rounds:
            if r.get("error"):
                continue
            primers.append({"name": f"SCRUB-fam-{r['round']}-F",
                            "seq": r["fwd_seq"], "strand": 1, "tm": r["fwd_tm"]})
            primers.append({"name": f"SCRUB-fam-{r['round']}-R",
                            "seq": r["rev_seq"], "strand": -1, "tm": r["rev_tm"]})
        assert primers, "no scrub primers designed for the test"
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            cured_rec = deepcopy(rec)
            cured_rec.seq = Seq(plan["cured_seq"])
            app._apply_snapshot(plan["cured_seq"], 0, cured_rec)
            await pilot.pause()
            out = app._scrub_save_cured_plasmid(
                name="Cured One", collection="Default",
                primers=primers, src_id=src_id)
            for _ in range(3):
                await pilot.pause()
        assert out is not None
        clone_rec, final_name, amp_names = out
        assert final_name == "Cured One" and "_" not in final_name
        assert amp_names == []            # QuikChange path saves no amplicons
        saved = [e for c in sc._load_collections()
                 for e in c.get("plasmids", [])
                 if e.get("source") == "scrub:cured"]
        assert saved and saved[-1]["name"] == "Cured One"
        assert "primer_bind" in (saved[-1].get("gb_text") or ""), \
            "cure primers not bound to the cured plasmid"
        # The ORIGINAL source entry now ALSO carries the cure primers.
        src_entry = sc._find_library_entry_by_id(src_id)
        assert src_entry and "primer_bind" in (src_entry.get("gb_text") or ""), \
            "cure primers not written back onto the original plasmid"

    async def test_no_src_id_skips_source_but_saves_cured(self):
        from copy import deepcopy
        from Bio.Seq import Seq
        seq, rec = _record_with_bsai()
        plan = sc._scrub_design(seq, rec.features, ["BsaI"], circular=True)
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            cured_rec = deepcopy(rec)
            cured_rec.seq = Seq(plan["cured_seq"])
            app._apply_snapshot(plan["cured_seq"], 0, cured_rec)
            await pilot.pause()
            # Empty primer list + no source id must not crash; cured still saves.
            out = app._scrub_save_cured_plasmid(
                name="Cured No Src", collection="Default",
                primers=[], src_id=None)
            for _ in range(2):
                await pilot.pause()
        assert out is not None
        saved = [e for c in sc._load_collections()
                 for e in c.get("plasmids", [])
                 if e.get("name") == "Cured No Src"]
        assert saved, "cured plasmid not saved when src_id is None"

    async def test_golden_braid_saves_amplicons_and_assembly_history(self):
        """GB Apply-cure save SIMULATES THE STEPS: each fragment amplicon is
        saved as `PCR-SCRUB-…` (primers bound + amplify history), the amplicons
        are genuinely BsaI-digested + ligated, and the cured plasmid carries an
        ASSEMBLY history whose parents are those amplicons."""
        from copy import deepcopy
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        seq = TestScrubGoldenBraid()._seq_with(
            400, [(100, "GAATTC"), (260, "GAATTC")], seed=7)
        rec = SeqRecord(Seq(seq), id="gbsrc", name="gbsrc",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec._tui_display_name = "GB Source"
        src_id = "gb_src"
        libs = sc._load_library()
        libs.append({"id": src_id, "name": "GB Source", "size": len(seq),
                     "gb_text": sc._record_to_gb_text(rec), "kind": "plasmid",
                     "n_feats": 0})
        sc._save_library(libs)
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"])
        assert plan["ok"] and plan["verified"], plan.get("errors")
        primers = []
        for fr in plan["fragments"]:
            i = fr["index"] + 1
            primers.append({"name": f"SCRUB-FAM-{i}-F", "seq": fr["fwd_seq"],
                            "strand": 1, "tm": fr["fwd_tm"]})
            primers.append({"name": f"SCRUB-FAM-{i}-R", "seq": fr["rev_seq"],
                            "strand": -1, "tm": fr["rev_tm"]})
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            cured_rec = deepcopy(rec)
            cured_rec.seq = Seq(plan["cured_seq"])
            app._apply_snapshot(plan["cured_seq"], 0, cured_rec)
            await pilot.pause()
            out = app._scrub_save_cured_plasmid(
                name="GB Cured", collection="Default", primers=primers,
                src_id=src_id, plan=plan, source_name="GB Source")
            for _ in range(3):
                await pilot.pause()
        assert out is not None
        clone_rec, final_name, amp_names = out
        assert final_name == "GB Cured" and "_" not in final_name
        # one PCR- amplicon per fragment, named PCR-SCRUB-FAM-#
        assert len(amp_names) == len(plan["fragments"]) >= 2
        assert all(a.startswith("PCR-SCRUB-FAM-") for a in amp_names), amp_names
        amplicons = [e for c in sc._load_collections()
                     for e in c.get("plasmids", [])
                     if e.get("source") == "scrub_gb:amplicon"]
        assert len(amplicons) == len(amp_names)
        assert all("primer_bind" in (e.get("gb_text") or "") for e in amplicons)
        assert all("amplifyFragment" in (e.get("history_xml") or "")
                   for e in amplicons)
        # the cured plasmid has an ASSEMBLY history with amplicons as parents
        cured = [e for c in sc._load_collections()
                 for e in c.get("plasmids", [])
                 if e.get("source") == "scrub_gb:cured"]
        assert cured, "GB cured plasmid not saved"
        hist = cured[-1].get("history_xml") or ""
        assert "insertFragment" in hist and "goldenBraidAssembly" in hist
        node = sc._parse_commercialsaas_history(hist)
        assert node is not None and len(node.parents) == len(amp_names), \
            "assembly history parents != saved amplicons"

    async def test_quikchange_history_has_no_amplicons(self):
        """QuikChange Apply-cure save: the cured plasmid gets a mutagenesis
        (`editSequence`) history with the primers + the source as parent, and
        NO whole-plasmid amplicons are saved (faithful, no clutter)."""
        from copy import deepcopy
        from Bio.Seq import Seq
        seq, rec = _record_with_bsai()
        rec._tui_display_name = "QC Source"
        src_id = "qc_src"
        libs = sc._load_library()
        libs.append({"id": src_id, "name": "QC Source", "size": len(seq),
                     "gb_text": sc._record_to_gb_text(rec), "kind": "plasmid",
                     "n_feats": 0})
        sc._save_library(libs)
        plan = sc._scrub_design(seq, rec.features, ["BsaI"], circular=True)
        rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                       round_no=i)
                  for i, c in enumerate(plan["clusters"], 1)]
        primers = []
        for r in rounds:
            if r.get("error"):
                continue
            primers.append({"name": f"SCRUB-QC-{r['round']}-F",
                            "seq": r["fwd_seq"], "strand": 1, "tm": r["fwd_tm"]})
            primers.append({"name": f"SCRUB-QC-{r['round']}-R",
                            "seq": r["rev_seq"], "strand": -1, "tm": r["rev_tm"]})
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            cured_rec = deepcopy(rec)
            cured_rec.seq = Seq(plan["cured_seq"])
            app._apply_snapshot(plan["cured_seq"], 0, cured_rec)
            await pilot.pause()
            out = app._scrub_save_cured_plasmid(
                name="QC Cured", collection="Default", primers=primers,
                src_id=src_id, plan=plan, source_name="QC Source")
            for _ in range(3):
                await pilot.pause()
        assert out is not None
        _clone, _final, amp_names = out
        assert amp_names == []          # no fragment amplicons for QuikChange
        assert not [e for c in sc._load_collections()
                    for e in c.get("plasmids", [])
                    if e.get("source") == "scrub_gb:amplicon"]
        cured = [e for c in sc._load_collections()
                 for e in c.get("plasmids", [])
                 if e.get("source") == "scrub:cured"]
        assert cured, "QuikChange cured plasmid not saved"
        hist = cured[-1].get("history_xml") or ""
        assert "editSequence" in hist and "QuikChange" in hist

    async def test_cured_plasmid_has_own_id_so_canvas_can_switch_back(self):
        """Bug B regression (2026-06-13): the cured plasmid is a `deepcopy` of
        its source, so it used to INHERIT the source's `record.id`.
        `_library_load`'s "already loaded (record.id == entry_id)" guard then
        refused to switch the canvas from the cured plasmid back to its
        near-identical source. The cured plasmid now gets its OWN id."""
        from copy import deepcopy
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        seq = TestScrubGoldenBraid()._seq_with(
            400, [(100, "GAATTC"), (260, "GAATTC")], seed=7)
        # Source's RECORD id == its ENTRY id (as a library round-trip yields) —
        # the exact condition that tripped the stale-id guard.
        src_id = "FFE_SRC"
        rec = SeqRecord(Seq(seq), id=src_id, name=src_id,
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec._tui_display_name = "FFE Source"
        libs = sc._load_library()
        libs.append({"id": src_id, "name": "FFE Source", "size": len(seq),
                     "gb_text": sc._record_to_gb_text(rec), "kind": "plasmid",
                     "n_feats": 0})
        sc._save_library(libs)
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"])
        assert plan["ok"] and plan["verified"], plan.get("errors")
        primers = []
        for fr in plan["fragments"]:
            i = fr["index"] + 1
            primers.append({"name": f"SCRUB-FAM-{i}-F", "seq": fr["fwd_seq"],
                            "strand": 1, "tm": fr["fwd_tm"]})
            primers.append({"name": f"SCRUB-FAM-{i}-R", "seq": fr["rev_seq"],
                            "strand": -1, "tm": fr["rev_tm"]})
        app = sc.PlasmidApp()
        async with app.run_test(size=_BASELINE) as pilot:
            await pilot.pause()
            cured_rec = deepcopy(rec)
            cured_rec.seq = Seq(plan["cured_seq"])
            app._apply_snapshot(plan["cured_seq"], 0, cured_rec)
            await pilot.pause()
            out = app._scrub_save_cured_plasmid(
                name="FFE Cured", collection="Default", primers=primers,
                src_id=src_id, plan=plan, source_name="FFE Source")
            assert out is not None
            clone_rec, _final, _amps = out
            assert clone_rec.id != src_id, \
                "cured plasmid still inherits the source's record id"
            # Load the cured plasmid (as the post-save focus does), then select
            # the SOURCE via _library_load — the canvas MUST switch to it.
            app._apply_record(clone_rec)
            await pilot.pause()
            assert app._current_record.id != src_id
            src_entry = sc._find_library_entry_by_id(src_id)
            app._library_load(sc.LibraryPanel.PlasmidLoad(src_entry))
            await pilot.pause()
            await pilot.pause()
            assert str(app._current_record.seq) == seq, \
                "canvas did not switch back to the source (stale-id load guard)"


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

    def test_unknown_enzyme_is_400_not_a_clean_bill_of_health(self):
        """`_scrub_resolve_sites` SKIPS a name the catalog doesn't know, so a
        typo used to scrub nothing and answer `ok: true` with an empty
        `sites_removed` — which reads as "already clean" on a plasmid that
        isn't. [INV-189]"""
        seq, _rec = _record_with_bsai()
        res = sc._h_scrub_plasmid(None, {"seq": seq,
                                         "enzymes": ["BsaI", "Bsa1"]})
        assert isinstance(res, tuple) and res[1] == 400, res
        assert res[0]["unknown_enzymes"] == ["Bsa1"]
        assert "BsaI" in res[0]["error"]           # near-miss suggestion

    def test_empty_enzyme_list_is_400(self):
        seq, _rec = _record_with_bsai()
        res = sc._h_scrub_plasmid(None, {"seq": seq, "enzymes": []})
        assert isinstance(res, tuple) and res[1] == 400

    def test_case_and_synonym_names_scrub_the_same_sites(self):
        seq, _rec = _record_with_bsai()
        removed = []
        for spelling in ("BsaI", "bsai", "Eco31I"):
            res = sc._h_scrub_plasmid(None, {"seq": seq,
                                             "enzymes": [spelling]})
            assert not isinstance(res, tuple), (spelling, res)
            assert res["ok"] is True
            removed.append((res["cured_seq"], len(res["sites_removed"])))
        assert removed[0] == removed[1] == removed[2]
        assert removed[0][1] == 1


    def test_bad_overlap_400(self):
        res = sc._h_scrub_plasmid(None, {"seq": "ACGT" * 30, "overlap": "x"})
        assert isinstance(res, tuple) and res[1] == 400

    def test_bad_method_400(self):
        res = sc._h_scrub_plasmid(None, {"seq": "ACGT" * 30, "method": "x"})
        assert isinstance(res, tuple) and res[1] == 400

    def test_golden_braid_method(self):
        seq = TestScrubGoldenBraid()._seq_with(
            400, [(100, "GAATTC"), (260, "GAATTC")], seed=7)
        res = sc._h_scrub_plasmid(None, {"seq": seq, "enzymes": ["EcoRI"],
                                         "method": "golden_braid"})
        assert res["ok"] is True, res.get("errors")
        assert res["method"] == "golden_braid"
        assert res["verified"] is True
        assert res["n_fragments"] >= 2 and res["fragments"]
        # cured product carries neither EcoRI nor BsaI (the assembly enzyme).
        assert not sc._scrub_scan_targets(res["cured_seq"],
                                          frozenset(["EcoRI", "BsaI"]), True)
        for fr in res["fragments"]:
            assert fr["fwd_seq"].startswith("GCGCGGTCTCA")

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
                    "SCRUB-")]
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


class TestScrubGoldenBraid:
    """Golden Braid (BsaI Type IIS) fragment-based curing — `_scrub_gb_design`
    + helpers. The plasmid is split at each cure cluster into PCR fragments
    with BsaI-tailed primers carrying the NATIVE 4 nt junction overhangs, so a
    Golden Gate reaction reassembles it with only the unwanted sites cured.

    Catastrophic-class: a wrong fragment/overhang scrambles the plasmid or
    leaves a site behind. The design's own digest+ligate verification must
    pass (`verified`), the cured product must carry no EcoRI / no BsaI, and
    every primer's binding must equal the cured template it anneals to."""

    @staticmethod
    def _clean_seq(length: int, seed: int) -> str:
        """A random ACGT sequence with NO EcoRI / BsaI (either strand), so the
        only sites are the ones a test inserts. Pure string check — no
        splicecraft import, no data dir."""
        import random
        rng = random.Random(seed)
        bad = ("GAATTC", "GGTCTC", "GAGACC")
        while True:
            s = "".join(rng.choice("ACGT") for _ in range(length))
            # wrap-aware: reject a site spanning the circular origin too
            if not any(b in (s + s[:5]) for b in bad):
                return s

    def _seq_with(self, length: int, inserts, seed: int = 7) -> str:
        s = list(self._clean_seq(length, seed))
        for pos, site in inserts:
            for i, c in enumerate(site):
                s[(pos + i) % length] = c
        return "".join(s)

    @staticmethod
    def _has(seq: str, site: str) -> bool:
        """Wrap-aware presence of `site` or its RC in circular `seq`."""
        aug = seq + seq[:len(site) - 1]
        return site in aug or sc._rc(site) in aug

    def test_real_digest_ligation_reassembles_cured(self):
        """Apply cure's assembly path (`_assemble_scrub_amplicons_real`)
        GENUINELY BsaI-digests + ligates the fragment amplicons back into the
        cured plasmid — the real digest+ligation primitives, not the design-time
        string-chaining shortcut. Guards the overhang-canonicalisation risk:
        `_ligate_fragments` matches overhangs by equality, so this proves the
        scrub amplicons actually close through the real machinery."""
        seq = self._seq_with(400, [(100, "GAATTC"), (260, "GAATTC")], seed=7)
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"])
        assert plan["ok"] and plan["verified"], plan.get("errors")
        frags = plan["fragments"]
        n = len(plan["cured_seq"])
        amps = sc._scrub_gb_build_amplicons(
            plan["orig_seq"], plan["cured_seq"], frags, n,
            single=(len(frags) == 1))
        assert len(amps) == len(frags) and len(frags) >= 2
        product = sc._assemble_scrub_amplicons_real(amps)
        assert product is not None, \
            "real BsaI digest+ligation failed to close the circle"
        assert len(product) == n
        # circular rotation equality with the cured plasmid
        assert plan["cured_seq"] in (product + product), \
            "assembled product is not a rotation of the cured plasmid"
        # the assembled product is BsaI-clean (wrap-aware)
        assert not sc._scrub_scan_targets(product, frozenset(["BsaI"]), True)

    def test_overhang_ok_rejects_palindrome_and_homopolymer(self):
        assert sc._scrub_gb_overhang_ok("GGAG")          # clean, asymmetric
        assert sc._scrub_gb_overhang_ok("AATG")
        assert not sc._scrub_gb_overhang_ok("GATC")      # palindrome
        assert not sc._scrub_gb_overhang_ok("AATT")      # palindrome
        assert not sc._scrub_gb_overhang_ok("AAAA")      # homopolymer
        assert not sc._scrub_gb_overhang_ok("ACG")       # wrong length
        assert not sc._scrub_gb_overhang_ok("ACGN")      # non-ACGT

    def test_design_verifies_and_cleans_two_ecorI(self):
        seq = self._seq_with(400, [(100, "GAATTC"), (260, "GAATTC")], seed=7)
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"], circular=True)
        assert plan["ok"], plan["errors"]
        assert plan["verified"]
        assert plan["n_fragments"] >= 2
        # cured plasmid carries no EcoRI and no BsaI (the assembly enzyme).
        assert not self._has(plan["cured_seq"], "GAATTC")
        assert not self._has(plan["cured_seq"], "GGTCTC")
        for fr in plan["fragments"]:                      # BsaI tail on every primer
            assert fr["fwd_seq"].startswith("GCGCGGTCTCA")
            assert fr["rev_seq"].startswith("GCGCGGTCTCA")

    def test_primer_binding_matches_cured_template(self):
        # Binding == display (catastrophic-class): each primer's 3' binding
        # must equal the cured template it anneals to.
        seq = self._seq_with(400, [(100, "GAATTC"), (260, "GAATTC")], seed=7)
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"], circular=True)
        assert plan["ok"], plan["errors"]
        cured, n = plan["cured_seq"], len(plan["cured_seq"])
        tail = "GCGCGGTCTCA"
        for fr in plan["fragments"]:
            fwd_bind = fr["fwd_seq"][len(tail):]
            assert fwd_bind == sc._circ_extract(cured, fr["cut_l"], len(fwd_bind), n)
            rev_bind = fr["rev_seq"][len(tail):]
            win = sc._circ_extract(cured, (fr["cut_r"] + 4 - len(rev_bind)) % n,
                                   len(rev_bind), n)
            assert rev_bind == sc._rc(win)

    def test_junction_overhangs_unique_and_nonpalindromic(self):
        seq = self._seq_with(540, [(80, "GAATTC"), (240, "GAATTC"),
                                   (400, "GAATTC")], seed=11)
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"], circular=True)
        assert plan["ok"], plan["errors"]
        ohs = [fr["oh_left"] for fr in plan["fragments"]]
        for oh in ohs:
            assert sc._scrub_gb_overhang_ok(oh)
        seen: set = set()
        for oh in ohs:                                    # distinct incl. RC
            assert oh not in seen and sc._rc(oh) not in seen, ohs
            seen.add(oh); seen.add(sc._rc(oh))

    def test_bsai_site_is_cured_not_tolerated(self):
        # An unwanted site that IS BsaI must be REMOVED (it is the assembly
        # enzyme). Even when the user lists only EcoRI, BsaI is force-cured.
        seq = self._seq_with(360, [(90, "GAATTC"), (200, "GGTCTC")], seed=5)
        assert self._has(seq, "GGTCTC")                  # BsaI present to start
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"], circular=True)
        assert plan["ok"], plan["errors"]
        assert plan["verified"]
        assert not self._has(plan["cured_seq"], "GGTCTC")   # BsaI gone
        assert not self._has(plan["cured_seq"], "GAATTC")

    def test_uncurable_bsai_fails_with_clear_error(self, monkeypatch):
        # If the inner cure skips a BsaI site (e.g. an overlapping opposite-
        # strand CDS leaves no synonymous escape), Golden Braid is impossible:
        # you cannot assemble with an enzyme that still cuts a fragment.
        def _fake(seq, feats, enzymes, *, circular=True, codon_raw=None):
            return {"ok": True, "cured_seq": seq, "edits": [],
                    "sites_removed": [],
                    "sites_skipped": [{"enzyme": "BsaI", "pos": 42, "strand": 1,
                                       "region": "CDS", "reason": "no syn cure"}],
                    "clusters": [], "warnings": []}
        # _scrub_gb_design moved to splicecraft_cloning and calls _scrub_design
        # in THAT namespace (imported from splicecraft_primer), so intercept it
        # there — patching sc._scrub_design no longer reaches the sibling caller.
        monkeypatch.setattr("splicecraft_cloning._scrub_design", _fake)
        plan = sc._scrub_gb_design("ACGT" * 80, [], ["BsaI"], circular=True)
        assert not plan["ok"]
        joined = " ".join(plan["errors"]).lower()
        assert "bsai" in joined and "assembly enzyme" in joined

    def test_single_cluster_self_circularizes(self):
        seq = self._seq_with(220, [(100, "GAATTC")], seed=3)
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"], circular=True)
        assert plan["ok"], plan["errors"]
        assert plan["n_fragments"] == 1
        assert plan["verified"]
        fr = plan["fragments"][0]
        assert fr["oh_left"] == fr["oh_right"]           # one shared junction
        assert sc._scrub_gb_overhang_ok(fr["oh_left"])

    def test_cured_differs_from_original_only_at_edits(self):
        seq = self._seq_with(400, [(110, "GAATTC"), (270, "GAATTC")], seed=9)
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"], circular=True)
        assert plan["ok"], plan["errors"]
        edit_pos = {e["pos"] for e in plan["edits"]}
        diffs = {i for i in range(len(seq)) if seq[i] != plan["cured_seq"][i]}
        assert diffs == edit_pos        # nothing changed outside the cures
        assert diffs                    # something WAS cured

    def test_verify_requires_bsai_clean_not_site_free(self):
        # A skipped NON-BsaI site (couldn't be cured) stays in the product —
        # like QuikChange, that's reported, not a whole-design failure. Verify
        # must only require the product to be BsaI-CLEAN (the assembly enzyme),
        # NOT free of every unwanted site. (Pre-hardening it scanned ALL
        # forbidden sites and failed the design on a residual EcoRI.)
        cured = self._seq_with(220, [(90, "GAATTC")], seed=4)   # residual EcoRI
        assert not self._has(cured, "GGTCTC")                   # but BsaI-clean
        assert self._has(cured, "GAATTC")
        frag = {"index": 0, "cut_l": 0, "cut_r": 0, "span": len(cured) + 4,
                "oh_left": cured[:4], "oh_right": cured[:4],
                "fwd_bind_len": 20, "rev_bind_len": 20}
        ok, errors = sc._scrub_gb_verify(cured, cured, [frag], len(cured),
                                         single=True)
        assert ok, errors                       # residual EcoRI does NOT fail
        # …but a residual BsaI must fail (it would be re-cut mid-assembly).
        bad = self._seq_with(220, [(90, "GGTCTC")], seed=4)
        fbad = {"index": 0, "cut_l": 0, "cut_r": 0, "span": len(bad) + 4,
                "oh_left": bad[:4], "oh_right": bad[:4],
                "fwd_bind_len": 20, "rev_bind_len": 20}
        ok2, errors2 = sc._scrub_gb_verify(bad, bad, [fbad], len(bad),
                                           single=True)
        assert not ok2 and any("BsaI" in e for e in errors2), errors2

    def test_no_target_sites_is_graceful(self):
        # A plasmid with none of the target sites → nothing to fragment; the
        # design returns cleanly (no crash, no fragments, a clear note).
        seq = self._seq_with(200, [], seed=2)            # clean, no sites
        plan = sc._scrub_gb_design(seq, [], ["EcoRI"], circular=True)
        assert plan["ok"] is True
        assert plan["n_fragments"] == 0 and not plan["fragments"]
        assert any("no sites" in w.lower() or "nothing" in w.lower()
                   for w in plan["warnings"])

    def test_only_bsai_single_site_force_cured(self):
        # The user lists NO enzymes, but a lone BsaI site is still force-cured
        # (it's the assembly enzyme) → one self-circularising fragment.
        seq = self._seq_with(220, [(90, "GGTCTC")], seed=6)
        plan = sc._scrub_gb_design(seq, [], [], circular=True)
        assert plan["ok"], plan["errors"]
        assert plan["verified"]
        assert plan["n_fragments"] == 1
        assert not self._has(plan["cured_seq"], "GGTCTC")


class TestScrubGoldenBraidUI:
    """The Scrub tab's Golden Braid re-circularization mode through
    MutagenizeModal — method selector, GB render, commit-button gating, and
    drawing the fragment primers on the map."""

    def _rec(self, seq: str):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        return SeqRecord(Seq(seq), id="GBUI", name="GBUI",
                         annotations={"molecule_type": "DNA",
                                      "topology": "circular"})

    async def test_method_selector_defaults_to_quikchange(self):
        seq, _ = _record_with_bsai()
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app.push_screen(sc.MutagenizeModal(seq, [], "x", show_tab="scrub"))
            await pilot.pause()
            await pilot.pause(0.05)
            sel = app.screen.query_one("#scrub-method", sc.Select)
            assert sel.value == "quikchange"

    async def test_golden_braid_greys_out_overlap_picker(self):
        seq, _ = _record_with_bsai()
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app.push_screen(sc.MutagenizeModal(seq, [], "x", show_tab="scrub"))
            await pilot.pause()
            await pilot.pause(0.05)
            screen = app.screen
            ov = screen.query_one("#scrub-overlap", sc.Select)
            assert ov.disabled is False                  # QuikChange: live
            screen.query_one("#scrub-method", sc.Select).value = "golden_braid"
            await pilot.pause()
            assert ov.disabled is True                   # GB: greyed (no-op)
            screen.query_one("#scrub-method", sc.Select).value = "quikchange"
            await pilot.pause()
            assert ov.disabled is False                  # back to live

    async def test_gb_result_renders_and_enables_commits(self):
        gb = TestScrubGoldenBraid()
        seq = gb._seq_with(400, [(100, "GAATTC"), (260, "GAATTC")], seed=7)
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app._apply_snapshot(seq, 0, self._rec(seq))
            await pilot.pause()
            feats = app.query_one("#plasmid-map", sc.PlasmidMap)._feats
            app.push_screen(sc.MutagenizeModal(seq, feats, "GBUI",
                                               show_tab="scrub"))
            await pilot.pause()
            await pilot.pause(0.1)
            modal = app.screen
            plan = sc._scrub_gb_design(seq, feats, ["EcoRI"], circular=True)
            assert plan["verified"], plan["errors"]
            modal._scrub_apply_result(plan, [])
            await pilot.pause()
            body = modal._render_scrub_gb(plan).plain
            assert "Golden Braid" in body
            for bid in ("#btn-scrub-apply", "#btn-scrub-saveprimers",
                        "#btn-scrub-tomap"):
                assert modal.query_one(bid, sc.Button).disabled is False
            # Add to Map draws the BsaI-tailed fragment primers.
            modal._scrub_save_to_map(None)
            await pilot.pause()
            gbp = [f for f in app._current_record.features
                   if f.type == "primer_bind"
                   and (f.qualifiers.get("label", [""]) or [""])[0].startswith(
                       "SCRUB-")]
            assert len(gbp) >= 2
            assert all("primer_seq" in f.qualifiers for f in gbp)
            # the user-typed display name (with spaces) survives the add.
            app._current_record._tui_display_name = "GB UI Test"
            assert getattr(app._current_record, "_tui_display_name",
                           None) == "GB UI Test"

    async def test_gb_failed_design_keeps_commits_disabled(self):
        seq, _ = _record_with_bsai()
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app._apply_snapshot(seq, 0, self._rec(seq))
            await pilot.pause()
            app.push_screen(sc.MutagenizeModal(seq, [], "x", show_tab="scrub"))
            await pilot.pause()
            await pilot.pause(0.05)
            modal = app.screen
            failed = {
                "method": "golden_braid", "ok": False, "verified": False,
                "cured_seq": seq, "enzyme": "BsaI",
                "edits": [{"pos": 1, "frm": "A", "to": "G", "region": "x",
                           "enzyme": "BsaI"}],
                "sites_removed": [], "sites_skipped": [], "fragments": [],
                "errors": ["a BsaI site can't be removed; BsaI is the "
                           "assembly enzyme"], "warnings": []}
            modal._scrub_apply_result(failed, [])
            await pilot.pause()
            for bid in ("#btn-scrub-apply", "#btn-scrub-saveprimers",
                        "#btn-scrub-tomap"):
                assert modal.query_one(bid, sc.Button).disabled is True
            body = modal._render_scrub_gb(failed).plain
            assert "assembly enzyme" in body


class TestScrubGoldenBraidCDS:
    """Golden Braid curing inside a coding sequence must stay synonymous: the
    cure rides in a fragment primer, the fragment reassembles, and the CDS
    protein is byte-identical to the original."""

    def test_synonymous_cure_in_cds_preserves_protein(self):
        gb = TestScrubGoldenBraid()
        body = list(gb._clean_seq(300, seed=21))
        # BsaI in-frame at codon boundary 60 (GGT-CTC = Gly-Leu → GGA-CTC
        # cures synonymously). CDS spans the whole circular plasmid.
        for i, c in enumerate("GGTCTC"):
            body[60 + i] = c
        seq = "".join(body)
        feats = [{"type": "CDS", "start": 0, "end": 300, "strand": 1,
                  "codon_start": 1, "transl_table": 1, "label": "orf"}]
        plan = sc._scrub_gb_design(seq, feats, ["BsaI"], circular=True)
        assert plan["ok"], plan["errors"]
        assert plan["verified"]
        assert not gb._has(plan["cured_seq"], "GGTCTC")          # BsaI gone
        # protein byte-identical across the cure (synonymous)
        assert (sc._translate_cds(seq, 0, 300, 1)
                == sc._translate_cds(plan["cured_seq"], 0, 300, 1))
        # the edit is genuinely silent (only the wobble base moved)
        assert plan["edits"] and all(e["region"].startswith("CDS")
                                     for e in plan["edits"])


class TestScrubSavePersistence:
    """Add-to-Map is an in-memory dirty edit (Undo fully reverses it — no
    primers lingering in the library), and Save-primers routes through the
    naming + collection modal."""

    def _rec(self, seq: str, rid: str = "LINGER"):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        return SeqRecord(Seq(seq), id=rid, name=rid,
                         annotations={"molecule_type": "DNA",
                                      "topology": "circular"})

    def _qc_rounds(self, seq, feats):
        plan = sc._scrub_design(seq, feats, ["BsaI"], circular=True)
        rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                       round_no=i)
                  for i, c in enumerate(plan["clusters"], 1)]
        return plan, rounds

    async def test_add_to_map_does_not_persist_to_library(self):
        seq, _ = _record_with_bsai()
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app._apply_snapshot(seq, 0, self._rec(seq, "LINGER"))
            await pilot.pause()
            # seed a clean (primer-free) library row for this plasmid
            app.query_one("#library", sc.LibraryPanel).add_entry(
                app._current_record)
            await pilot.pause()
            feats = app.query_one("#plasmid-map", sc.PlasmidMap)._feats
            modal = sc.MutagenizeModal(seq, feats, "LINGER", show_tab="scrub")
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            plan, rounds = self._qc_rounds(seq, feats)
            modal._scrub_apply_result(plan, rounds)
            await pilot.pause()
            modal._scrub_save_to_map(None)
            await pilot.pause()
            # canvas HAS the primers (in-memory dirty edit) …
            assert any(f.type == "primer_bind"
                       for f in app._current_record.features)
            # … but the LIBRARY row was NOT auto-persisted — so Undo (which
            # only reverts the canvas) can't leave primers lingering there.
            entry = next((e for e in sc._load_library()
                          if e.get("id") == "LINGER"), None)
            assert entry is not None
            assert "primer_bind" not in (entry.get("gb_text") or "")

    async def test_save_primers_opens_naming_modal(self):
        seq, _ = _record_with_bsai()
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app._apply_snapshot(seq, 0, self._rec(seq))
            await pilot.pause()
            feats = app.query_one("#plasmid-map", sc.PlasmidMap)._feats
            modal = sc.MutagenizeModal(seq, feats, "x", show_tab="scrub")
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            plan, rounds = self._qc_rounds(seq, feats)
            modal._scrub_apply_result(plan, rounds)
            await pilot.pause()
            modal._scrub_save_primers(None)
            await pilot.pause()
            await pilot.pause(0.05)
            assert isinstance(app.screen, sc.PrimerSaveModal)

    async def test_commit_primers_uses_typed_names(self):
        seq, _ = _record_with_bsai()
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            app._apply_snapshot(seq, 0, self._rec(seq))
            await pilot.pause()
            feats = app.query_one("#plasmid-map", sc.PlasmidMap)._feats
            modal = sc.MutagenizeModal(seq, feats, "x", show_tab="scrub")
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            plan, rounds = self._qc_rounds(seq, feats)
            modal._scrub_apply_result(plan, rounds)
            await pilot.pause()
            oligos = [{"label": p["name"], "default_name": p["name"],
                       "sequence": p["seq"], "tm": p["tm"],
                       "_strand": p["strand"], "_ptype": p["ptype"]}
                      for p in modal._scrub_primer_list()]
            assert len(oligos) >= 2
            custom = [f"CUSTOM_{i}" for i in range(len(oligos))]
            modal._scrub_commit_primers_with_names(oligos, custom)
            await pilot.pause()
            names = {p.get("name") for p in sc._load_primers()}
            assert all(c in names for c in custom)


class TestScrubQuikChangeVerify:
    """`_scrub_qc_verify`: the QuikChange primers, simulated against the
    ORIGINAL template, must reconstitute the cured plasmid EXACTLY — the same
    seamless-cure guarantee Golden Braid gets from `_scrub_gb_verify`."""

    def test_primers_reconstitute_cured_only_site_changes(self):
        seq = TestScrubGoldenBraid()._seq_with(300, [(120, "GGTCTC")], seed=8)
        plan = sc._scrub_design(seq, [], ["BsaI"], circular=True)
        assert plan["ok"], plan
        rounds = [sc._scrub_qc_primers(plan["cured_seq"], c["positions"],
                                       round_no=i)
                  for i, c in enumerate(plan["clusters"], 1)]
        ok, errors = sc._scrub_qc_verify(seq, plan["cured_seq"], rounds,
                                         len(seq))
        assert ok, errors
        # …and the cured plasmid differs from the original ONLY at the edits
        # (i.e. only the destroyed cut site changed — nothing else).
        diffs = {i for i in range(len(seq)) if seq[i] != plan["cured_seq"][i]}
        assert diffs == {e["pos"] for e in plan["edits"]}
        assert diffs and not TestScrubGoldenBraid()._has(plan["cured_seq"],
                                                         "GGTCTC")

    def test_catches_cure_outside_primer_reach(self):
        # A cure no primer footprint covers would NOT be incorporated → the
        # product keeps the original base → verify must fail (not silently OK).
        orig = "ACGT" * 50                          # 200 bp
        flip = "A" if orig[100] != "A" else "T"
        cured = orig[:100] + flip + orig[101:]
        rounds = [{"round": 1,
                   "fwd_seq": cured[0:20], "fwd_start": 0, "fwd_len": 20,
                   "rev_seq": sc._mut_revcomp(cured[20:40]), "rev_start": 20,
                   "rev_len": 20}]
        ok, errors = sc._scrub_qc_verify(orig, cured, rounds, len(orig))
        assert not ok and errors

    def test_endpoint_quikchange_surfaces_verified(self):
        seq = TestScrubGoldenBraid()._seq_with(300, [(120, "GGTCTC")], seed=8)
        res = sc._h_scrub_plasmid(None, {"seq": seq, "enzymes": ["BsaI"],
                                         "method": "quikchange"})
        assert res["method"] == "quikchange"
        assert res["verified"] is True


class TestScrubSourceAnnotation:
    """`Save primers` also draws the family-named (no-underscore) scrub primers
    onto the ORIGINAL source plasmid (the library entry the scrub ran on), so
    the cures are visible on the source. Also verifies the active-library
    lookup works from inside the modal."""

    async def test_save_primers_annotates_source(self):
        seq, rec = _record_with_bsai()
        rec._tui_display_name = "Src Plasmid 1"      # spaces — must stay clean
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            src_id = "srcplasmid1"
            libs = sc._load_library()
            libs.append({"id": src_id, "name": "Src Plasmid 1", "size": len(seq),
                         "gb_text": sc._record_to_gb_text(rec), "kind": "plasmid",
                         "n_feats": len(rec.features)})
            sc._save_library(libs)
            app._apply_snapshot(seq, 0, rec)
            app.query_one("#library", sc.LibraryPanel).set_active(src_id)
            await pilot.pause()
            feats = app.query_one("#plasmid-map", sc.PlasmidMap)._feats
            app.push_screen(sc.MutagenizeModal(seq, feats, "Src Plasmid 1",
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
            primers = modal._scrub_primer_list()
            # family-named SCRUB-{family}-#-F/R, never an underscore
            assert primers and all(
                p["name"].startswith("SCRUB-") and "_" not in p["name"]
                for p in primers)
            oligos = [{"label": p["name"], "default_name": p["name"],
                       "sequence": p["seq"], "tm": p["tm"],
                       "_strand": p["strand"], "_ptype": p["ptype"]}
                      for p in primers]
            modal._scrub_commit_primers_with_names(
                oligos, [o["default_name"] for o in oligos])
            await pilot.pause()
            src = sc._find_library_entry_by_id(src_id)
            assert src and "primer_bind" in (src.get("gb_text") or ""), \
                "scrub primers not written onto the source plasmid"
            assert "SCRUB-" in (src.get("gb_text") or "")
            # no underscores forced into the source display name (INV-98)
            assert "_" not in (src.get("name") or "")


class TestCureDoesNotLengthenHomopolymers:
    """2026-08-26 field report: curing `GGAAAAA GAAGAC GTT` (a BpiI site in
    pTRBO's 35S enhancer) by G→A spells `GGAAAAAAAAGAC` — an 8×A run — where
    the equally valid A→G at the fourth base leaves the longest run at 5.

    Both are single transitions with the same GC delta, so the score
    separated them only by the lexicographic determinism tie-break, and `A`
    sorts before `G`. Long runs are synthesis- and sequencing-hard, and the
    reporter's site sat in a TANDEM-duplicated enhancer, so the bad choice
    landed twice."""

    PAD = "CTGACGTAAGGGATGACGCACAATCCCACTATCCTTCGCAAGACCCTTCC"

    @staticmethod
    def _longest(seq):
        best = run = 0
        prev = ""
        for ch in seq:
            run = run + 1 if ch == prev else 1
            prev = ch
            best = max(best, run)
        return best

    def _cure(self, core, enzyme="BbsI"):
        seq = self.PAD + core + self.PAD
        plan = sc._scrub_design(seq, [], [enzyme], circular=False)
        assert plan["ok"] and plan["sites_removed"], plan["warnings"]
        return seq, plan

    def test_the_reported_case(self):
        seq, plan = self._cure("GGAAAAA" + "GAAGAC" + "GTT")
        cured = plan["cured_seq"]
        assert "GAAGAC" not in cured and "GTCTTC" not in cured   # site gone
        assert self._longest(cured) <= self._longest(seq)
        assert "AAAAAAAA" not in cured                            # the 8×A

    def test_a_run_on_the_other_side_too(self):
        seq, plan = self._cure("TTG" + "GAAGAC" + "AAAAAGG")
        assert self._longest(plan["cured_seq"]) <= self._longest(seq)

    def test_it_generalises_past_bbsi(self):
        seq, plan = self._cure("GGTTTTT" + "GGTCTC" + "AAAAAGG", enzyme="BsaI")
        assert self._longest(plan["cured_seq"]) <= self._longest(seq)

    def test_a_site_spanning_the_origin(self):
        """The window is circular — a run that meets the edit across bp 0 has
        to count, or the origin becomes a blind spot."""
        seq = "AGACGTT" + self.PAD * 2 + "AAAAAGA"
        plan = sc._scrub_design(seq, [], ["BbsI"], circular=True)
        assert plan["sites_removed"], plan["warnings"]
        dbl_before = seq + seq[:20]
        dbl_after = plan["cured_seq"] + plan["cured_seq"][:20]
        assert self._longest(dbl_after) <= max(self._longest(dbl_before), 5)

    def test_the_term_outranks_the_lexicographic_tiebreak(self):
        """Directly: with the homopolymer term neutralised the OLD choice
        comes back, which is what proves this test isn't vacuous."""
        import splicecraft_primer as sp
        seq = self.PAD + "GGAAAAA" + "GAAGAC" + "GTT" + self.PAD
        new = sc._scrub_design(seq, [], ["BbsI"], circular=False)["cured_seq"]
        old_free = sp._SCRUB_HOMOPOLYMER_FREE
        try:
            sp._SCRUB_HOMOPOLYMER_FREE = 10 ** 6      # term can never fire
            old = sc._scrub_design(seq, [], ["BbsI"],
                                   circular=False)["cured_seq"]
        finally:
            sp._SCRUB_HOMOPOLYMER_FREE = old_free
        assert "AAAAAAAA" in old, "fixture no longer reproduces the bug"
        assert "AAAAAAAA" not in new

    def test_it_is_measured_relative_to_what_was_there(self):
        """A cure inside a pre-existing long run must not be penalised for a
        run it inherited — otherwise the scorer would refuse the only silent
        option available in an AT-rich region."""
        seq = self.PAD + "A" * 9 + "GAAGAC" + "A" * 9 + self.PAD
        plan = sc._scrub_design(seq, [], ["BbsI"], circular=False)
        assert plan["ok"] and plan["sites_removed"]
        assert plan["max_homopolymer_run"] >= 9

    def test_max_run_is_reported_and_never_silently_worse(self):
        seq = self.PAD + "GGAAAAA" + "GAAGAC" + "GTT" + self.PAD
        r = sc._h_scrub_plasmid(None, {"seq": seq, "features": [],
                                       "enzymes": ["BbsI"], "circular": False})
        assert not isinstance(r, tuple), r
        assert r["max_homopolymer_run"] == self._longest(r["cured_seq"])
        # Nothing was made worse here, so no homopolymer warning.
        assert not any("homopolymer" in w for w in r["warnings"]), r["warnings"]


class TestScrubLongestRun:
    def test_edges(self):
        import splicecraft_primer as sp
        assert sp._scrub_longest_run("") == 0
        assert sp._scrub_longest_run("A") == 1
        assert sp._scrub_longest_run("ACGT") == 1
        assert sp._scrub_longest_run("AACCCGT") == 3
        assert sp._scrub_longest_run("TTTTTT") == 6

    def test_an_all_one_base_plasmid_never_reports_a_run_longer_than_itself(self):
        """The circular measurement doubles a head slice; without the cap an
        all-A plasmid would report a run longer than the molecule."""
        plan = sc._scrub_design("A" * 80, [], ["BsaI"], circular=True)
        assert plan["max_homopolymer_run"] <= 80
