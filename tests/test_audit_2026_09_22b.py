"""Regression tests for the RESIDUALS of the 2026-09-22 deep audit.

The first pass (`test_audit_2026_09_22.py`, shipped v1.2.71) took the findings
that could put wrong biology in front of a user. This file covers the tail that
was triaged and left: the same discipline — every assertion checks against
something outside the code under test (a brute-force count computed in the
guide's own frame, a coordinate derived from the variant table by hand, a
Biopython re-parse) rather than against the implementation's own idea of itself.
"""
from __future__ import annotations

import random

import pytest

import splicecraft as sc
import splicecraft_crispr as _crispr
import splicecraft_regulatory as _reg
import splicecraft_seqanalysis as _sa


# ══════════════════════════════════════════════════════════════════════════
# CRISPR
# ══════════════════════════════════════════════════════════════════════════

# ── N2: only the guide's OWN site is dropped from its off-target list ──────

class TestDuplicatedTargetIsAnOffTarget:
    """A target that occurs twice cuts twice. Dropping every zero-mismatch hit
    (rather than the on-target one by coordinate) reported `n_offtargets: 0`
    for exactly the guide that is most dangerous."""

    SPACER = "ACGTTGCAAGGCTTACGTCA"

    def _dup_plasmid(self):
        random.seed(11)
        pad = lambda n: "".join(random.choice("ACGT") for _ in range(n))
        return (pad(300) + self.SPACER + "AGG" + pad(300)
                + self.SPACER + "AGG" + pad(300))

    def test_a_second_perfect_site_is_reported(self):
        seq = self._dup_plasmid()
        r = _crispr._design_guides(seq, variant="spcas9", circular=False,
                                   offtarget_in=seq, limit=200)
        hits = [g for g in r["guides"] if g["guide"] == self.SPACER]
        assert len(hits) == 2, "probe failed to plant two copies"
        for g in hits:
            assert g["n_offtargets"] >= 1
            perfect = [h["start"] for h in g["offtargets"]
                       if h["mismatches"] == 0]
            # the OTHER copy, never this one
            assert perfect and g["start"] not in perfect

    def test_the_on_target_hit_is_still_excluded(self):
        seq = self._dup_plasmid()
        r = _crispr._design_guides(seq, variant="spcas9", circular=False,
                                   offtarget_in=seq, limit=200)
        for g in r["guides"]:
            assert not any(h["start"] == g["start"]
                           and h["strand"] == g["strand"]
                           and h["mismatches"] == 0
                           for h in g.get("offtargets", []))

    def test_a_perfect_hit_in_a_foreign_sequence_is_kept(self):
        # Searching a host genome: there is no on-target hit to drop at all.
        seq = self._dup_plasmid()
        r = _crispr._design_guides(seq, variant="spcas9", circular=False,
                                   offtarget_in=self.SPACER + "AGG",
                                   limit=200)
        g = next(x for x in r["guides"] if x["guide"] == self.SPACER)
        assert g["n_offtargets"] == 1

    def test_a_lowercase_copy_of_the_molecule_is_the_same_molecule(self):
        # The sequence is normalised to upper case; the off-target input was
        # compared raw, so a lowercase caller got every guide's own site back
        # as a perfect off-target.
        seq = self._dup_plasmid().lower()
        upper = _crispr._design_guides(seq.upper(), circular=True,
                                       offtarget_in=seq.upper(), limit=500)
        lower = _crispr._design_guides(seq, circular=True,
                                       offtarget_in=seq, limit=500)
        count = lambda r: {(g["start"], g["strand"]): g["n_offtargets"]
                           for g in r["guides"]}
        assert count(lower) == count(upper)
        assert lower["offtarget_searched_bp"] == len(seq)


# ── AA11: the dialog's off-target search does not freeze the app ──────────────

class TestOfftargetSearchLeavesTheAppResponsive:
    """The off-target pass is the scan's inner loop once per guide (about a
    minute and a half at 20 kb) and the dialog ran it on the UI thread."""

    SEQ = "ACGTACGTAGGT" * 40

    async def _open(self, app, pilot):
        await pilot.pause()
        modal = sc.CrisprGuideModal(self.SEQ, circular=True, label="pTest")
        app.push_screen(modal)
        await pilot.pause()
        return modal

    async def test_the_search_runs_off_the_ui_thread(
            self, monkeypatch, tiny_record, isolated_library):
        import threading
        import splicecraft_modals as md
        real, on_ui = md._design_guides, []

        def spy(*a, **k):
            if k.get("offtarget_in"):
                on_ui.append(threading.current_thread()
                             is threading.main_thread())
            return real(*a, **k)
        monkeypatch.setattr(md, "_design_guides", spy)
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            modal = await self._open(app, pilot)
            modal._on_offtarget(None)
            await app.workers.wait_for_complete()
            await pilot.pause()
            assert on_ui == [False]
            assert modal._guides and all(
                g.get("n_offtargets") is not None for g in modal._guides)

    async def test_a_late_search_does_not_overwrite_a_newer_design(
            self, monkeypatch, tiny_record, isolated_library):
        import threading
        import splicecraft_modals as md
        real, gate = md._design_guides, threading.Event()

        def slow(*a, **k):
            if k.get("offtarget_in"):
                gate.wait(10)
            return real(*a, **k)
        monkeypatch.setattr(md, "_design_guides", slow)
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=(160, 48)) as pilot:
            modal = await self._open(app, pilot)
            modal._on_offtarget(None)
            await pilot.pause()
            assert modal.query_one("#crispr-btn-offtarget").disabled
            modal._run_design(offtarget=False)     # the user re-designs
            gate.set()
            await app.workers.wait_for_complete()
            await pilot.pause()
            assert modal._guides and all(
                g.get("n_offtargets") is None for g in modal._guides)
            assert not modal.query_one("#crispr-btn-offtarget").disabled


# ── N3: seed mismatches are counted at the PAM-proximal end on both strands ─

def _ref_seed_mismatches(guide, window_fwd, strand, seed_len, three_prime):
    """Reference: count in the GUIDE's own frame, with no forward-strand
    slicing anywhere — the calculation the module is trying to approximate."""
    observed = window_fwd if strand == 1 else sc._rc(window_fwd)
    sl = (slice(len(guide) - seed_len, len(guide)) if three_prime
          else slice(0, seed_len))
    return sum(1 for a, b in zip(guide[sl], observed[sl]) if a != b)


class TestSeedMismatchesMeasuredOnThePamProximalEnd:
    def test_minus_strand_hits_match_a_guide_frame_reference(self):
        random.seed(7)
        checked = wrong = 0
        for _ in range(250):
            g = "".join(random.choice("ACGT") for _ in range(20))
            var = list(sc._rc(g))
            for i in random.sample(range(20), random.randint(1, 3)):
                var[i] = random.choice([c for c in "ACGT" if c != var[i]])
            pre = "".join(random.choice("ACGT") for _ in range(50))
            post = "".join(random.choice("ACGT") for _ in range(50))
            seq = pre + "CCT" + "".join(var) + post   # CCT = rc(AGG)
            for h in _crispr._guide_offtargets(
                    g, seq, variant="spcas9", circular=False,
                    max_mismatch=4)["hits"]:
                if h["strand"] != -1:
                    continue
                window = seq[h["start"]:h["start"] + 20]
                checked += 1
                if h["seed_mismatches"] != _ref_seed_mismatches(
                        g, window, -1, _crispr._SEED_LEN, True):
                    wrong += 1
        assert checked > 50, "probe planted too few minus-strand hits"
        assert wrong == 0

    def test_plus_strand_hits_are_unchanged(self):
        random.seed(8)
        g = "ACGTTGCAAGGCTTACGTCA"
        var = list(g)
        var[1] = "T" if g[1] != "T" else "A"       # one seed-distal mismatch
        seq = "TTTT" + "".join(var) + "AGG" + "TTTT"
        h = next(x for x in _crispr._guide_offtargets(
            g, seq, variant="spcas9", circular=False)["hits"]
            if x["strand"] == 1)
        assert h["mismatches"] == 1
        assert h["seed_mismatches"] == _ref_seed_mismatches(
            g, seq[h["start"]:h["start"] + 20], 1, _crispr._SEED_LEN, True)

    def test_a_five_prime_pam_variant_seeds_at_the_other_end(self):
        # lbcas12a's PAM is 5' of the spacer, so its seed is the spacer's START
        # on the plus strand and its END on the minus strand — the mirror of
        # the Cas9 case, which is why one slice for both could not be right.
        assert (_crispr._CAS_VARIANTS["lbcas12a"]["pam_side"] == "5prime")
        random.seed(9)
        g = "".join(random.choice("ACGT") for _ in range(23))
        var = list(sc._rc(g))
        var[0] = random.choice([c for c in "ACGT" if c != var[0]])
        seq = ("GGGG" + "".join(var) + sc._rc("TTTA") + "GGGG")
        hits = [h for h in _crispr._guide_offtargets(
            g, seq, variant="lbcas12a", circular=False,
            max_mismatch=3)["hits"] if h["strand"] == -1]
        assert hits
        for h in hits:
            window = seq[h["start"]:h["start"] + 23]
            assert h["seed_mismatches"] == _ref_seed_mismatches(
                g, window, -1, _crispr._SEED_LEN, False)


# ── N8: the off-target cap must not swallow a whole strand ─────────────────

class TestOffTargetCapIsStrandFair:
    def _many_plus_one_perfect_minus(self):
        g = "ACGTACGTACGTACGTACGT"
        random.seed(3)
        blocks = []
        for _ in range(40):
            v = list(g)
            for i in random.sample(range(20), 3):
                v[i] = random.choice([c for c in "ACGT" if c != v[i]])
            blocks.append("".join(v) + "AGG"
                          + "".join(random.choice("ACGT") for _ in range(10)))
        return g, "".join(blocks) + "CCT" + sc._rc(g)

    def test_the_minus_strand_is_still_scanned_when_the_cap_fills(self):
        g, seq = self._many_plus_one_perfect_minus()
        r = _crispr._guide_offtargets(g, seq, variant="spcas9",
                                      circular=False, max_mismatch=3,
                                      max_hits=10)
        assert r["truncated"] is True
        assert -1 in {h["strand"] for h in r["hits"]}

    def test_the_worst_hit_survives_the_cap(self):
        g, seq = self._many_plus_one_perfect_minus()
        r = _crispr._guide_offtargets(g, seq, variant="spcas9",
                                      circular=False, max_mismatch=3,
                                      max_hits=10)
        assert len(r["hits"]) <= 10
        assert any(h["strand"] == -1 and h["mismatches"] == 0
                   for h in r["hits"]), "a perfect hit lost to the cap"

    def test_truncation_is_reported_through_design_guides(self):
        g, seq = self._many_plus_one_perfect_minus()
        r = _crispr._design_guides(seq, variant="spcas9", circular=False,
                                   offtarget_in=seq, limit=5)
        assert any("offtarget_truncated" in x for x in r["guides"])


# ── N9: a staggered cutter reports BOTH nicks and the overhang ─────────────

class TestStaggeredCutIsReportedOnBothStrands:
    def _both_strand_staggered(self):
        random.seed(5)
        pad = lambda n: "".join(random.choice("ACGT") for _ in range(n))
        core = pad(23)
        seq = (pad(40) + "TTTA" + core + pad(40)
               + sc._rc("TTTA" + core) + pad(40))
        guides = _crispr._find_guides(seq, variant="lbcas12a", circular=False)
        plus = next(g for g in guides
                    if g["strand"] == 1 and g["guide"] == core)
        minus = next(g for g in guides
                     if g["strand"] == -1 and g["guide"] == core)
        return plus, minus

    def test_nicks_match_the_variant_table_by_hand(self):
        top_off, bot_off = _crispr._CAS_VARIANTS["lbcas12a"]["cut_offsets"]
        glen = _crispr._CAS_VARIANTS["lbcas12a"]["guide_len"]
        plus, minus = self._both_strand_staggered()
        # Plus: offsets count forward from the protospacer 5' end.
        assert plus["cut_site"] == plus["start"] + top_off
        assert plus["cut_site_bottom"] == plus["start"] + bot_off
        # Minus: the guide strand runs right-to-left, so the guide-strand nick
        # is the forward-frame BOTTOM one and the offsets swap roles.
        assert minus["cut_site"] == minus["start"] + glen - bot_off
        assert minus["cut_site_bottom"] == minus["start"] + glen - top_off

    def test_the_overhang_has_the_same_sign_on_both_strands(self):
        plus, minus = self._both_strand_staggered()
        assert plus["overhang"] == minus["overhang"] == 5

    def test_a_blunt_variant_reports_no_stagger(self):
        random.seed(6)
        seq = ("".join(random.choice("ACGT") for _ in range(200))
               + "TTTTTTTTTTTTTTTTTTTTAGG")
        guides = _crispr._find_guides(seq, variant="spcas9", circular=False)
        assert guides
        for g in guides:
            assert g["overhang"] == 0
            assert g["cut_site"] == g["cut_site_bottom"]


# ── N11: a wrapped region on a linear molecule is refused, not answered 0 ──

class TestWrappedRegionOnALinearMolecule:
    def test_it_refuses_instead_of_returning_no_guides(self):
        random.seed(12)
        seq = "".join(random.choice("ACGT") for _ in range(800))
        with pytest.raises(ValueError, match="linear"):
            _crispr._design_guides(seq, circular=False, region=(500, 100))

    def test_a_circular_molecule_still_wraps(self):
        random.seed(12)
        seq = "".join(random.choice("ACGT") for _ in range(800))
        r = _crispr._design_guides(seq, circular=True, region=(700, 100))
        n = len(seq)
        assert r["n_found"] >= 1
        for g in r["guides"]:
            assert g["cut_site"] >= 700 or g["cut_site"] < 100, \
                f"cut {g['cut_site']} outside the wrapped region"
        assert n == 800

    def test_the_endpoint_turns_it_into_a_400(self):
        import splicecraft_agent as _agent
        random.seed(12)
        seq = "".join(random.choice("ACGT") for _ in range(800))
        got = _agent._h_design_guides(
            None, {"sequence": seq, "circular": False, "region": [500, 100]})
        assert isinstance(got, tuple) and got[1] == 400
        assert "linear" in got[0]["error"]


# ── N11a: the endpoint's documented body keys are the ones it reads ────────

def test_guide_cloning_oligos_documents_the_keys_it_reads():
    import splicecraft_agent as _agent
    doc = _agent._h_guide_cloning_oligos.__doc__ or ""
    assert "overhang_top" in doc and "overhang_bottom" in doc
    assert "top_prefix" not in doc and "bottom_prefix" not in doc
    # and they really are the keys that take effect
    got = _agent._h_guide_cloning_oligos(
        None, {"guide": "GACGTTGCAAGGCTTACGTCA",
               "overhang_top": "AAAA", "overhang_bottom": "TTTT"})
    assert got["top"].startswith("AAAA")
    assert got["bottom"].startswith("TTTT")


# ══════════════════════════════════════════════════════════════════════════
# Regulatory
# ══════════════════════════════════════════════════════════════════════════

# ── N10: the documented reporting floor is inclusive ───────────────────────

class TestTerminatorReportingFloorIsInclusive:
    """`_TERM_MAX_DG`'s calibration comment documents `dg <= -8` as the floor.
    The guard read `>=`, so a hairpin folding at exactly the floor was dropped —
    the scan and its own published threshold disagreed by one value."""

    SEQ = "GGGGCCCCGCGCAAGCGCGGGGCCCCTTTTTTTT"

    def _scan_with_fold(self, monkeypatch, dg):
        monkeypatch.setattr(_reg, "_rna_fold", lambda s: ("", dg))
        return _reg._scan_terminators(self.SEQ, circular=False,
                                      both_strands=False)

    def test_a_hairpin_at_exactly_the_floor_is_reported(self, monkeypatch):
        hits = self._scan_with_fold(monkeypatch, _reg._TERM_MAX_DG)
        assert hits, f"dg == {_reg._TERM_MAX_DG} (the documented floor) dropped"
        assert all(h["dg"] == pytest.approx(_reg._TERM_MAX_DG) for h in hits)

    def test_a_hairpin_just_above_the_floor_is_still_dropped(self, monkeypatch):
        assert self._scan_with_fold(monkeypatch, _reg._TERM_MAX_DG + 0.1) == []

    def test_a_hairpin_below_the_floor_is_reported(self, monkeypatch):
        assert self._scan_with_fold(monkeypatch, _reg._TERM_MAX_DG - 0.1)


# ── N7: one origin-spanning hairpin is ONE record, whatever the rotation ───

class TestArcsOverlap:
    """The circular-overlap primitive the merge now rests on, checked against a
    brute-force base-set intersection."""

    @staticmethod
    def _bases(start, length, total):
        return {(start + i) % total for i in range(length)}

    def test_matches_a_brute_force_base_intersection(self):
        random.seed(41)
        total = 37
        for _ in range(3000):
            a0, al = random.randrange(total), random.randint(1, total - 1)
            b0, bl = random.randrange(total), random.randint(1, total - 1)
            want = bool(self._bases(a0, al, total)
                        & self._bases(b0, bl, total))
            got = _reg._arcs_overlap(a0, al, b0, bl, total)
            assert got is want, (a0, al, b0, bl, total, got, want)

    def test_a_zero_length_arc_touches_nothing(self):
        assert _reg._arcs_overlap(5, 0, 5, 10, 100) is False

    def test_a_full_lap_touches_everything(self):
        assert _reg._arcs_overlap(0, 100, 55, 1, 100) is True


class TestTerminatorScanIsRotationInvariant:
    """One hairpin is ONE record at every rotation. Overlap used to be judged in
    PADDED coordinates, and the two registers of an origin-spanning hairpin are
    discovered in different padded frames — so at exactly the rotations that put
    a terminator across the origin, it came back twice."""

    TERM = "GGGGCCCCGCGCAAGCGCGGGGCCCCTTTTTTTT"

    @staticmethod
    def _canon(hits, n, shift):
        """Map a rotation-`shift` hit set back into frame 0."""
        return {((h["start"] + shift) % n, (h["end"] + shift) % n,
                 h["strand"], h["stem_len"], h["u_tract"]) for h in hits}

    def _seq(self, seed, pad_len=60):
        random.seed(seed)
        pad = "".join(random.choice("ACGT") for _ in range(pad_len * 2))
        return pad[:pad_len] + self.TERM + pad[pad_len:]

    def test_every_rotation_gives_the_same_hit_set(self):
        seq = self._seq(31)
        n = len(seq)
        base = self._canon(_reg._scan_terminators(seq, circular=True), n, 0)
        assert base, "probe planted no terminator"
        for shift in range(n):
            rot = seq[shift:] + seq[:shift]
            got = self._canon(
                _reg._scan_terminators(rot, circular=True), n, shift)
            assert got == base, (
                f"rotation {shift}: extra={sorted(got - base)} "
                f"missing={sorted(base - got)}")

    def test_no_rotation_produces_duplicate_records(self):
        seq = self._seq(32)
        n = len(seq)
        for shift in range(n):
            rot = seq[shift:] + seq[:shift]
            keys = [(h["start"], h["end"], h["strand"])
                    for h in _reg._scan_terminators(rot, circular=True)]
            assert len(keys) == len(set(keys)), (
                f"duplicate terminator records at rotation {shift}: {keys}")

    def test_two_terminators_stay_two_at_every_rotation(self):
        random.seed(33)
        pad = lambda k: "".join(random.choice("ACGT") for _ in range(k))
        seq = pad(40) + self.TERM + pad(40) + self.TERM + pad(40)
        n = len(seq)
        base = self._canon(_reg._scan_terminators(seq, circular=True), n, 0)
        assert len(base) == 2
        for shift in range(0, n, 3):
            rot = seq[shift:] + seq[:shift]
            got = self._canon(
                _reg._scan_terminators(rot, circular=True), n, shift)
            assert got == base, f"rotation {shift} disagrees"

    def test_a_linear_molecule_still_merges_registers(self):
        seq = self._seq(34)
        hits = _reg._scan_terminators(seq, circular=False)
        keys = [(h["start"], h["end"], h["strand"]) for h in hits]
        assert len(keys) == len(set(keys))


# ══════════════════════════════════════════════════════════════════════════
# Transcription map
# ══════════════════════════════════════════════════════════════════════════

# ── N12: an arrowless CDS has no reading direction, so both are tried ──────

class TestArrowlessGeneIsNotAssumedForward:
    """`map-transcription` read `strand` with `or 1`, so an arrowless CDS was
    silently forward — and the reverse promoter actually pointed at it was
    filtered out before it could be considered. The gene came back `silent`."""

    @staticmethod
    def _seq(n=600):
        random.seed(51)
        return "".join(random.choice("ACGT") for _ in range(n))

    def _feats(self, gene_strand):
        # A reverse promoter sitting to the RIGHT of the gene, firing leftwards
        # into its 3'-on-the-forward-strand end.
        return [
            {"start": 300, "end": 340, "strand": -1, "type": "promoter",
             "label": "pRev"},
            {"start": 100, "end": 280, "strand": gene_strand, "type": "CDS",
             "label": "orfX"},
        ]

    def test_a_reverse_promoter_reaches_an_arrowless_gene(self):
        r = _sa._map_transcription(self._seq(), self._feats(0),
                                   circular=True, include_predicted=False)
        gene = next(g for g in r["genes"] if g["cds"]["label"] == "orfX")
        assert gene["cds"]["strand_unknown"] is True
        assert gene["n_sources"] >= 1, "arrowless gene still reads as silent"
        assert gene["silent"] is False

    def test_the_unknown_strand_is_reported_not_invented(self):
        r = _sa._map_transcription(self._seq(), self._feats(0),
                                   circular=True, include_predicted=False)
        gene = next(g for g in r["genes"] if g["cds"]["label"] == "orfX")
        assert gene["cds"]["strand"] == 0
        assert any("no strand" in w for w in r["warnings"])

    def test_a_declared_strand_still_filters_by_strand(self):
        # The same layout with the gene declared FORWARD: the reverse promoter
        # must NOT be offered, or the fix has broken the ordinary case.
        r = _sa._map_transcription(self._seq(), self._feats(1),
                                   circular=True, include_predicted=False)
        gene = next(g for g in r["genes"] if g["cds"]["label"] == "orfX")
        assert gene["cds"]["strand"] == 1
        assert gene["cds"]["strand_unknown"] is False
        assert gene["n_sources"] == 0
        assert not any("no strand" in w for w in r["warnings"])

    def test_the_strand_helper_keeps_its_promoter_default(self):
        # A promoter must point somewhere to be modelled, so its caller keeps
        # the historical +1 for "no arrow" — absent OR an explicit 0 (harden
        # 2026-09-24: the explicit 0 used to escape the default, and an
        # arrowless promoter then drove nothing); only the gene caller asks
        # for 0.
        assert _sa._tx_feature_strand({}) == 1
        assert _sa._tx_feature_strand({}, default=0) == 0
        assert _sa._tx_feature_strand({"strand": 0}) == 1
        assert _sa._tx_feature_strand({"strand": 0}, default=0) == 0
        assert _sa._tx_feature_strand({"strand": -1}) == -1
        assert _sa._tx_feature_strand({"strand": 2}) == 1
        assert _sa._tx_feature_strand({"strand": None}, default=0) == 0


# ── AA11: the quadratic endpoints are concurrency-capped ──────────────────

def test_offtarget_and_fold_endpoints_are_heavy_classified():
    # Measured at ~40 s for `design-guides` with `offtarget: true` on a 20 kb
    # record; an uncapped burst pins every worker thread.
    for name in ("design-guides", "guide-offtargets", "fold-rna"):
        assert name in sc._AGENT_HEAVY_ENDPOINTS, f"{name} is not capped"
    # every heavy name must be a real endpoint, or the cap silently applies to
    # nothing
    for name in sc._AGENT_HEAVY_ENDPOINTS:
        assert name in sc._state._AGENT_HANDLERS, f"{name} is not an endpoint"


# ══════════════════════════════════════════════════════════════════════════
# Codon / translation
# ══════════════════════════════════════════════════════════════════════════

def _cds_record(seq, start, end, quals, strand=1, circular=True):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    r = SeqRecord(Seq(seq), id="probe", name="probe")
    r.annotations["molecule_type"] = "DNA"
    if circular:
        r.annotations["topology"] = "circular"
    r.features = [SeqFeature(FeatureLocation(start, end, strand=strand),
                             type="CDS", qualifiers=quals)]
    return r


def _cds_detail(rec):
    """The row `find-feature` builds, without needing a live app."""
    loc = rec.features[0].location
    row = {"idx": 0, "label": "cds", "type": "CDS",
           "start": int(loc.start), "end": int(loc.end),
           "strand": loc.strand, "length": len(loc), "wraps": False}
    return sc._agent_feature_detail(rec, row, max_seq=10_000)


# ── C5: an alternative initiator is not a mismatch ─────────────────────────

class TestAlternativeInitiatorIsNotAMismatch:
    """A bacterial CDS starting GTG is loaded with fMet, so the file says `M`
    where the codon encodes `V`. `_translate_cds` reports the codon honestly;
    the CONVENTION belongs to the comparison, which used to call a perfectly
    correct annotation a mismatch."""

    CDS = "GTGAAAGTTCTGGCC" + "TAA"          # V K V L A * by the table

    def _rec(self, translation="MKVLA", table="11"):
        seq = "AAAA" + self.CDS + "TTTT"
        return _cds_record(seq, 4, 4 + len(self.CDS),
                           {"translation": [translation],
                            "transl_table": [table]})

    def test_the_translation_itself_stays_honest(self):
        assert _cds_detail(self._rec())["translation"].rstrip("*") == "VKVLA"

    def test_a_bacterial_gtg_start_agrees_with_the_file(self):
        d = _cds_detail(self._rec())
        assert d["translation_matches_annotation"] is True
        assert "initiator" in d["translation_note"]

    def test_table_one_does_not_accept_a_gtg_initiator(self):
        # GTG is an initiator in table 11, not in table 1 — so the same bases
        # with the same annotation must still disagree there.
        assert "GTG" not in sc._codon_start_codons_for(1)
        assert "GTG" in sc._codon_start_codons_for(11)
        d = _cds_detail(self._rec(table="1"))
        assert d["translation_matches_annotation"] is False

    def test_a_real_disagreement_past_position_one_is_still_caught(self):
        d = _cds_detail(self._rec(translation="MKVLW"))
        assert d["translation_matches_annotation"] is False

    def test_the_start_codon_sets_match_biopython(self):
        from Bio.Data import CodonTable
        for tid in (1, 2, 4, 5, 6, 11, 25):
            want = frozenset(
                CodonTable.unambiguous_dna_by_id[tid].start_codons)
            assert sc._codon_start_codons_for(tid) == want


# ── C6: /transl_except declares the residue that may differ ────────────────

class TestTranslExceptIsHonoured:
    """A selenoprotein's file translation shows `U` where the bases read a
    premature `*`. The record DECLARES that position; a comparison that ignores
    the declaration gave a correctly recoded stop the same verdict as a real
    frameshift."""

    CDS = "ATGAAATGACTGGCC" + "TAA"          # M K * L A *, TGA at residue 3

    def _rec(self, quals=None, strand=1):
        seq = "GGGG" + (self.CDS if strand == 1 else sc._rc(self.CDS)) + "CCCC"
        q = {"translation": ["MKULA"]}
        q.update(quals or {})
        return _cds_record(seq, 4, 4 + len(self.CDS), q, strand=strand)

    def test_a_declared_recoded_stop_agrees(self):
        # residue 3 = CDS bases 7..9 = genomic 1-based 11..13
        d = _cds_detail(self._rec({"transl_except": ["(pos:11..13,aa:Sec)"]}))
        assert d["translation"].rstrip("*") == "MK*LA"
        assert d["translation_matches_annotation"] is True
        assert "transl_except" in d["translation_note"]

    def test_without_the_qualifier_it_is_still_a_mismatch(self):
        assert _cds_detail(self._rec())[
            "translation_matches_annotation"] is False

    def test_a_qualifier_at_another_position_excuses_nothing(self):
        d = _cds_detail(self._rec({"transl_except": ["(pos:20..22,aa:Sec)"]}))
        assert d["translation_matches_annotation"] is False
        assert "/transl_except" in d["translation_note"]

    def test_a_minus_strand_cds_maps_through_reading_order(self):
        # Reading order runs genomic high->low, so residue 3 is NOT at
        # (base/3); dividing the raw base number by three would name the wrong
        # residue on this record.
        n = len(self.CDS)
        hi0 = 4 + n - 1 - 3 * (3 - 1)          # first base of residue 3
        d = _cds_detail(self._rec(
            {"transl_except": [f"(pos:{hi0 - 1}..{hi0 + 1},aa:Sec)"]},
            strand=-1))
        assert d["translation"].rstrip("*") == "MK*LA"
        assert d["translation_matches_annotation"] is True

    def test_no_coding_map_means_no_exemption(self):
        # The helper relaxes a comparison, so with no reading-order map it must
        # relax nothing rather than guess.
        assert sc._transl_except_positions(["(pos:11..13,aa:Sec)"], None) == set()
        assert sc._transl_except_positions(["(pos:11..13,aa:Sec)"],
                                          list(range(4, 22))) == {3}

    def test_a_partial_codon_span_is_ignored(self):
        # Two bases cannot be one residue; guessing which is the wrong answer.
        assert sc._transl_except_positions(["(pos:11..12,aa:Sec)"],
                                           list(range(4, 22))) == set()


# ── C7: an ambiguous terminator is warned about ────────────────────────────

class TestAmbiguousTerminatorWarns:
    PROTEIN = "M" + "AK" * 40

    @staticmethod
    def _ambiguous_per_biopython(tid):
        from Bio.Data import CodonTable
        ct = CodonTable.unambiguous_dna_by_id[tid]
        return tuple(sorted(set(ct.forward_table) & set(ct.stop_codons)))

    def test_the_ambiguous_sets_match_biopython(self):
        import splicecraft_codon as _cod
        for tid in range(1, 34):
            try:
                want = self._ambiguous_per_biopython(tid)
            except KeyError:
                continue
            assert _cod._codon_ambiguous_stops(tid) == want, tid

    def test_only_three_tables_have_them(self):
        import splicecraft_codon as _cod
        from Bio.Data import CodonTable
        known = set(CodonTable.unambiguous_dna_by_id)
        got = {t for t in range(1, 34)
               if t in known and _cod._codon_ambiguous_stops(t)}
        assert got == {27, 28, 31}

    def _run(self, table, protein=None, stops=1):
        import splicecraft_agent as _agent
        payload = {"protein": protein or self.PROTEIN, "stops": stops}
        if table:
            payload["transl_table"] = table
        return _agent._h_optimize_protein(None, payload)

    def test_an_ambiguous_stop_warns(self):
        for tid in (27, 28, 31):
            got = self._run(tid)
            assert [w for w in got["warnings"] if "BOTH a stop" in w], tid

    def test_a_normal_table_does_not_warn(self):
        for tid in (1, 11, 6):
            got = self._run(tid)
            assert not [w for w in got["warnings"] if "BOTH a stop" in w], tid

    def test_a_code_with_no_safe_stop_says_so(self):
        # Table 28's three stops are ALL ambiguous, so there is no alternative
        # to offer — the warning has to say that rather than implying a fix.
        msg = [w for w in self._run(28)["warnings"] if "BOTH a stop" in w][0]
        assert "no unconditional terminator" in msg
        msg27 = [w for w in self._run(27)["warnings"] if "BOTH a stop" in w][0]
        assert "no unconditional terminator" in msg27   # TGA is its only stop


# ── C8: the GC3 floor and GC3 itself cover the coding body only ────────────

class TestGc3ExcludesAppendedStops:
    def _run(self, protein, stops):
        import splicecraft_agent as _agent
        return _agent._h_optimize_protein(
            None, {"protein": protein, "stops": stops})

    def test_appended_stops_do_not_pad_a_short_body_past_the_floor(self):
        import splicecraft_codon as _cod
        floor = _cod._CODON_GC3_MIN_CODONS
        got = self._run("M" * (floor - 1), 3)
        assert got["n_codons"] == floor + 2      # the stops ARE in the output
        assert not [w for w in got["warnings"] if "GC3" in w]

    def test_a_long_enough_body_still_warns(self):
        import splicecraft_codon as _cod
        floor = _cod._CODON_GC3_MIN_CODONS
        got = self._run("M" * (floor + 5), 1)
        assert [w for w in got["warnings"] if "GC3" in w]

    def test_gc3_is_measured_over_the_body(self):
        import splicecraft_codon as _cod
        got = self._run("MAKAKAKAKAKAKAKAKAKAK", 1)
        assert got["gc3"] == round(_cod._codon_gc3(got["dna"][:-3]), 2)


# ══════════════════════════════════════════════════════════════════════════
# Restriction
# ══════════════════════════════════════════════════════════════════════════

import splicecraft_biology as _bio          # noqa: E402


def _circ_record(seq, circular=True):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    r = SeqRecord(Seq(seq), id="p", name="p")
    r.annotations["molecule_type"] = "DNA"
    if circular is not None:
        r.annotations["topology"] = "circular" if circular else "linear"
    return r


def _app_with(rec):
    return type("A", (), {"_current_record": rec})()


def _filler(n, step=7):
    return "".join("ACGT"[(i * step) % 4] for i in range(n))


# ── X7: a full-lap span contains everything ────────────────────────────────

class TestFullLapSpanContainment:
    """`span-contains` reduced its coordinates `% total` before asking, which
    turned a whole-molecule `[0, n]` outer into the empty span `(0, 0)` — the
    single case `_span_in_span` had been hardened for."""

    N = 1000

    def _run(self, outer, inner):
        import splicecraft_agent as _agent
        return _agent._h_span_contains(
            None, {"length": self.N, "outer": outer, "inner": inner})

    def test_a_whole_molecule_outer_contains_every_span(self):
        got = self._run([0, self.N], [[10, 20], [900, 100], [0, self.N]])
        assert got["outer"]["length"] == self.N
        assert got["all_contained"] is True

    def test_a_full_lap_inner_is_not_a_point(self):
        got = self._run([500, 600], [0, self.N])
        assert got["results"][0]["length"] == self.N
        assert got["results"][0]["contains"] is False

    def test_an_ordinary_outer_still_excludes(self):
        got = self._run([100, 200], [[110, 120], [150, 250]])
        assert [r["contains"] for r in got["results"]] == [True, False]

    def test_span_in_span_matches_a_brute_force_base_set(self):
        # The reference is the actual SET of bases each span covers — computed
        # without any of the module's arithmetic.
        n = 11

        def bases(s, e):
            if e != s and (e - s) % n == 0:
                return set(range(n))
            return {(s + i) % n for i in range((e - s) % n)}

        for os_ in range(-n, 2 * n, 3):
            for oe in range(-n, 2 * n, 3):
                ob = bases(os_, oe)
                for is_ in range(-n, 2 * n, 2):
                    for ie in range(-n, 2 * n, 2):
                        ib = bases(is_, ie)
                        want = (ib <= ob) if ib else ((is_ % n) in ob)
                        assert sc._span_in_span(is_, ie, os_, oe, n) is want, (
                            is_, ie, os_, oe)

    def test_full_lap_length_helper(self):
        assert sc._span_full_lap_len(0, 100, 100) == 100
        assert sc._span_full_lap_len(0, 0, 100) == 0
        assert sc._span_full_lap_len(90, 10, 100) == 20
        assert sc._span_full_lap_len(50, 150, 100) == 100


# ── X9: IUPAC codes match by base-set containment, both directions ──────────

class TestIupacSubsetRule:
    """A site code `S` is satisfied by a sequence code `Q` when
    `bases(Q) ⊆ bases(S)`. The classes were written as the site's OWN base set,
    so a sequence carrying `GGTNACC` did not match BstEII's `GGTNACC` — the
    scan reported a site absent from a sequence that certainly contains it."""

    CODES = "ACGTRYWSMKBDHVN"

    def test_every_code_pair_follows_the_subset_rule(self):
        for site in self.CODES:
            pat = sc._iupac_pattern(site)
            for q in self.CODES:
                want = _bio._IUPAC_BASES[q] <= _bio._IUPAC_BASES[site]
                assert bool(pat.fullmatch(q)) is want, (site, q)

    def test_an_ambiguous_site_position_accepts_an_ambiguous_base(self):
        # BstEII GGTNACC — the enzyme does not care what that base is.
        assert sc._iupac_pattern("GGTNACC").fullmatch("GGTNACC")
        assert sc._iupac_pattern("GGTNACC").fullmatch("GGTAACC")

    def test_an_unknown_base_where_the_enzyme_needs_one_is_not_a_match(self):
        # Conservative in the other direction: an N might be anything, so it
        # must not manufacture a site.
        assert not sc._iupac_pattern("GAATTC").fullmatch("GAATTN")
        assert not sc._iupac_pattern("GAATTC").fullmatch("NNNNNN")

    def test_a_run_of_n_creates_no_sites(self):
        hits = sc._scan_restriction_sites("N" * 500, unique_only=False,
                                          circular=True)
        assert [h for h in hits if h.get("type") == "resite"] == []

    def test_the_regex_classes_are_derived_from_one_base_table(self):
        # Two hand-written maps is how one code comes to mean two things.
        for code, cls in _bio._IUPAC_RE.items():
            assert cls == _bio._iupac_class_for(code)


# ── X16: an empty allow-set is not "no filter" ─────────────────────────────

class TestEmptyAllowSetIsNotNoFilter:
    """`frozenset()` means "this collection resolved to no usable enzyme, scan
    nothing"; `None` means "no filter". The cache key folded the first onto the
    second, so depending on call order an empty set could be answered with the
    FULL catalog's sites — undoing the fix that made it return the empty set."""

    SEQ = "GAATTCAAAAGGATCCTTTTAAGCTTCCCCC" * 20
    KW = dict(circular=True, unique_only=False)

    def test_no_filter_first_then_empty(self):
        sc._state._RESTR_SCAN_CACHE.clear()
        full = sc._scan_restriction_sites(self.SEQ, **self.KW)
        empty = sc._scan_restriction_sites(
            self.SEQ, allowed_enzymes=frozenset(), **self.KW)
        assert len(full) > 0
        assert empty == []

    def test_empty_first_then_no_filter(self):
        sc._state._RESTR_SCAN_CACHE.clear()
        empty = sc._scan_restriction_sites(
            self.SEQ, allowed_enzymes=frozenset(), **self.KW)
        full = sc._scan_restriction_sites(self.SEQ, **self.KW)
        assert empty == []
        assert len(full) > 0

    def test_a_one_enzyme_set_is_its_own_slot(self):
        sc._state._RESTR_SCAN_CACHE.clear()
        full = sc._scan_restriction_sites(self.SEQ, **self.KW)
        one = sc._scan_restriction_sites(
            self.SEQ, allowed_enzymes=frozenset({"EcoRI"}), **self.KW)
        assert 0 < len(one) < len(full)


# ── X12: `wraps` means the recognition site crosses the origin ─────────────

class TestWrapsFlagMeansWraps:
    @staticmethod
    def _ecori_row(seq, circular=True):
        import splicecraft_agent as _agent
        got = _agent._h_list_restriction_sites(
            _app_with(_circ_record(seq, circular)),
            {"enzymes": ["EcoRI"], "unique_only": False})
        return [r for r in got["sites"] if r["enzyme"] == "EcoRI"]

    def test_a_site_across_the_origin_wraps(self):
        rows = self._ecori_row("TTC" + _filler(200) + "GAA")
        assert rows and all(r["wraps"] for r in rows)

    def test_a_site_ending_on_the_last_base_does_not_wrap(self):
        rows = self._ecori_row(_filler(200) + "GAATTC")
        assert rows and not any(r["wraps"] for r in rows)

    def test_a_site_at_the_start_does_not_wrap(self):
        rows = self._ecori_row("GAATTC" + _filler(200))
        assert rows and not any(r["wraps"] for r in rows)


# ── X13: the exported map honours the record's topology ────────────────────

class TestMapImageHonoursTopology:
    SEQ = "TTC" + _filler(400) + "GAA"      # EcoRI across the origin

    def test_a_linear_record_grows_no_origin_spanning_ticks(self):
        import splicecraft_mapimage as _mi
        circ = _mi._map_sites_from_record(_circ_record(self.SEQ, True),
                                          unique_only=False)
        lin = _mi._map_sites_from_record(_circ_record(self.SEQ, False),
                                         unique_only=False)
        assert len(circ) > len(lin), \
            "a linear molecule still shows the wrap-only cut"

    def test_an_unannotated_record_is_a_circle(self):
        # The app-wide rule (`_record_is_circular`): no topology line = circle,
        # which is what the on-screen map draws.
        import splicecraft_mapimage as _mi
        none = _mi._map_sites_from_record(_circ_record(self.SEQ, None),
                                         unique_only=False)
        circ = _mi._map_sites_from_record(_circ_record(self.SEQ, True),
                                         unique_only=False)
        assert len(none) == len(circ)


# ── X10: the map tick and the sequence-panel arrow mark the same base ──────

def test_reverse_bound_type_iis_tick_matches_the_cut_arrow():
    seq = _filler(120) + "GAGACC" + _filler(120)   # rc(GGTCTC) = BsaI reverse
    hits = sc._scan_restriction_sites(
        seq, circular=True, unique_only=False,
        allowed_enzymes=frozenset({"BsaI"}))
    resite = next(h for h in hits if h["type"] == "resite")
    recut = next(h for h in hits if h["type"] == "recut")
    assert resite["strand"] == -1
    # The panel draws its arrow at `ext_cut_bp`; the map draws its tick at the
    # recut's start. Marking different cuts put them 4 bp apart.
    assert recut["start"] == resite["ext_cut_bp"] == resite["top_cut_bp"]
    assert resite["bottom_cut_bp"] != resite["top_cut_bp"]   # still staggered


def _panel_marks(resite, n):
    """The bases the sequence panel draws a cut arrow at, for one resite."""
    marks = set()
    if resite.get("cut_col") is not None:
        marks.add((resite["start"] + resite["cut_col"]) % n)
    if resite.get("ext_cut_bp") is not None:
        marks.add(resite["ext_cut_bp"])
    return marks


def test_every_reverse_bound_site_ticks_where_its_arrow_points():
    """X10, all of it. Pinning the tick to the top cut fixed BsaI and moved it
    off the arrow for nine others (BsmI, BsrI, BtsI, …) whose arrow marks the
    bottom cut. Every non-palindromic catalog enzyme, bound in reverse: the
    tick must be a base the panel's arrow marks — the top cut whenever the
    arrow marks it."""
    enzymes = sc._all_enzymes()
    bad = []
    for name, (site, _fwd, _rev) in sorted(enzymes.items()):
        rc = sc._rc(site)
        if sc._iupac_pattern(rc).pattern == sc._iupac_pattern(site).pattern:
            continue                                   # palindromic
        concrete = "".join("A" if b not in "ACGT" else b for b in rc)
        seq = "A" * 90 + concrete + "A" * 90
        hits = sc._scan_restriction_sites(
            seq, circular=True, unique_only=False,
            allowed_enzymes=frozenset({name}))
        rev_sites = [h for h in hits
                     if h["type"] == "resite" and h["strand"] == -1
                     and h.get("label")]
        ticks = {h["start"] for h in hits
                 if h["type"] == "recut" and h["strand"] == -1}
        for rs in rev_sites:
            marks = _panel_marks(rs, len(seq))
            want = ({rs["top_cut_bp"]} if rs["top_cut_bp"] in marks
                    else marks)
            if not (ticks & want):
                bad.append((name, sorted(ticks), sorted(marks)))
    assert not bad, bad[:10]


# ── X5: isoschizomer spellings ride along on an unnamed listing ────────────

class TestIsoschizomersReachableWithoutNamingThem:
    """The scan labels one name per recognition site — right for the map, wrong
    for a query. Naming an enzyme already re-expanded the collapse (INV-187);
    an UNFILTERED listing still showed one spelling, so a caller hunting
    `Esp3I` read the plasmid as not carrying it."""

    SEQ = _filler(150) + "CGTCTC" + _filler(150)      # BsmBI / Esp3I

    def _sites(self):
        import splicecraft_agent as _agent
        return _agent._h_list_restriction_sites(
            _app_with(_circ_record(self.SEQ)), {"unique_only": False})["sites"]

    def test_the_other_spellings_are_reported(self):
        rows = [r for r in self._sites()
                if "Esp3I" in (r.get("aliases") or []) or r["enzyme"] == "Esp3I"]
        assert rows, "no row mentions Esp3I at all"

    def test_aliases_are_omitted_when_there_is_only_one_name(self):
        for r in self._sites():
            if r.get("aliases") is not None:
                assert r["aliases"], "an empty alias list was emitted"

    def test_the_alias_set_is_symmetric(self):
        for name in ("BsmBI", "Esp3I", "BsaI", "Eco31I"):
            for other in sc._enzyme_aliases(name):
                assert name in sc._enzyme_aliases(other), (name, other)


# ── X8: Dam / Dcm methylation is reported, never adjudicated ───────────────

class TestMethylationContextIsWarnedNotJudged:
    def test_methylated_positions_match_a_hand_count(self):
        #        0123456789...
        # seq =  AAGATCAACCAGGAA
        # GATC starts at 2 -> Dam methylates the A at 3 (top) and the A of the
        # bottom strand's GATC at 4.  CCAGG starts at 8 -> Dcm methylates the
        # inner C at 9 (top) and the bottom strand's inner C at 11.
        seq = "AAGATCAACCAGGAA"
        assert seq[2:6] == "GATC" and seq[8:13] == "CCAGG"
        got = _bio._methylated_base_positions(seq, circular=False)
        assert got == {3: "Dam", 4: "Dam", 9: "Dcm", 11: "Dcm"}

    def test_it_is_wrap_aware(self):
        # GATC split across the origin: "TC" + ... + "GA"
        seq = "TC" + _filler(60) + "GA"
        got = _bio._methylated_base_positions(seq, circular=True)
        n = len(seq)
        assert ((n - 1) % n) in got and 0 in got
        assert 0 not in _bio._methylated_base_positions(seq, circular=False)

    def test_an_enzyme_whose_site_always_contains_the_target(self):
        assert _bio._enzyme_site_methylation_targets("TGATCA") == ["Dam"]
        assert _bio._enzyme_site_methylation_targets("GATC") == ["Dam"]
        assert _bio._enzyme_site_methylation_targets("CCWGG") == ["Dcm"]
        assert _bio._enzyme_site_methylation_targets("GAATTC") == []

    def test_an_ambiguous_site_counts_only_when_it_cannot_avoid_the_target(self):
        # RGATCY (BstYI) always contains GATC; GGWCC never does.
        assert _bio._enzyme_site_methylation_targets("RGATCY") == ["Dam"]
        assert _bio._enzyme_site_methylation_targets("GGWCC") == []
        # An N cannot GUARANTEE a target base, so NGATCN is Dam and GNTC is not.
        assert _bio._enzyme_site_methylation_targets("GNTC") == []

    def test_the_endpoint_warns_without_claiming_a_verdict(self):
        import splicecraft_agent as _agent
        seq = _filler(100) + "TGATCA" + _filler(100) + "GAATTC" + _filler(100)
        got = _agent._h_list_restriction_sites(
            _app_with(_circ_record(seq)),
            {"enzymes": ["BclI", "EcoRI"], "unique_only": False})
        by = {r["enzyme"]: r for r in got["sites"]}
        assert by["BclI"]["methylation_context"] == ["Dam"]
        assert "methylation_context" not in by["EcoRI"]
        warn = " ".join(got.get("warnings", []))
        assert "Dam" in warn and "BclI" in warn
        # It must NOT claim BclI is blocked without naming the counter-examples,
        # because "contains the target" is not "is blocked".
        assert "DpnI requires it" in warn and "BamHI" in warn

    def test_list_enzymes_reports_the_enzyme_level_property(self):
        import splicecraft_agent as _agent
        rows = {r["name"]: r for r in
                _agent._h_list_enzymes(None, {})["enzymes"]}
        assert rows["BclI"]["methylation_targets"] == ["Dam"]
        assert rows["EcoRI"]["methylation_targets"] == []
        # BamHI CONTAINS GATC and is NOT blocked — the field must report the
        # target, so a consumer cannot read it as a verdict.
        assert rows["BamHI"]["methylation_targets"] == ["Dam"]


# ══════════════════════════════════════════════════════════════════════════
# Cloning
# ══════════════════════════════════════════════════════════════════════════

import splicecraft_cloning as _clo          # noqa: E402


def _clean_of(names, n, seed):
    """A random n-mer carrying none of `names`' sites on either strand."""
    rng = random.Random(seed)
    while True:
        s = "".join(rng.choice("ACGT") for _ in range(n))
        if not _clo._enzyme_cuts(s, list(names), circular=False):
            return s


# ── L3 / L3b: an origin-spanning marker tags ONE fragment ──────────────────

class TestWrapMarkerTagsOneFragment:
    """`int(loc.start)` on a CompoundLocation returns `min(parts.start)`
    (sacred #9), so an AmpR annotated across bp 0 flattened to a near-whole-
    plasmid feature. The splitter then put a piece of it on BOTH halves of the
    digest, both reported a backbone marker, and the caller fell back to the
    size heuristic these functions exist to replace."""

    def _vector(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        filler = _filler(600)
        seq = filler[:300] + "GGTCTCA" + "T" * 20 + "TGAGACC" + filler[300:]
        n = len(seq)
        rec = SeqRecord(Seq(seq), id="pEV", name="pEV")
        rec.annotations.update({"molecule_type": "DNA",
                                "topology": "circular"})
        rec.features = [SeqFeature(
            CompoundLocation([FeatureLocation(550, n, strand=1),
                              FeatureLocation(0, 120, strand=1)]),
            type="CDS", qualifiers={"label": ["AmpR"]})]
        return rec, seq, n

    def test_the_marshalled_feature_keeps_its_wrap_encoding(self):
        rec, _seq, _n = self._vector()
        for fn in (_clo._acceptor_vector_features,
                   _sa._ev_frag_input_features):
            row = next(f for f in fn(rec) if f.get("label") == "AmpR")
            assert row["end"] < row["start"], (
                f"{fn.__name__} flattened the wrap to "
                f"{row['start']}..{row['end']}")

    def test_exactly_one_digest_half_carries_the_marker(self):
        rec, seq, _n = self._vector()
        frags, err = _clo._excise_fragment_pair(
            seq, ["BsaI"], features=_sa._ev_frag_input_features(rec),
            circular=True)
        assert err is None and len(frags) == 2
        tagged = [f for f in frags if _sa._fragment_has_backbone_marker(f)]
        assert len(tagged) == 1
        # and it is the BACKBONE half, not the 16 bp dropout
        assert len(tagged[0]["top_seq"]) > 100

    def test_the_flattened_form_would_have_tagged_both(self):
        # Pins the mechanism, so a future flatten fails here rather than
        # silently returning to the size heuristic.
        _rec, seq, n = self._vector()
        flat = [{"start": 0, "end": n, "type": "CDS", "label": "AmpR"}]
        frags, _err = _clo._excise_fragment_pair(
            seq, ["BsaI"], features=flat, circular=True)
        assert sum(1 for f in frags
                   if _sa._fragment_has_backbone_marker(f)) == 2


# ── CL2 / CL3: Golden Gate tries both orientations ─────────────────────────

class TestGoldenGateTriesBothOrientations:
    """A Golden Gate fragment is double-stranded; its OVERHANGS make the
    assembly directional, not which strand the file stores. A part supplied
    reverse-complemented presents the same overhangs at swapped ends, and
    trying only the stored orientation answered "the overhangs don't chain"."""

    @staticmethod
    def _frag(top, left_oh, right_oh, label):
        return {
            "top_seq": top,
            "left":  {"kind": "5'", "overhang_seq": left_oh,
                      "enzyme": "BsaI"},
            "right": {"kind": "5'", "overhang_seq": right_oh,
                      "enzyme": "BsaI"},
            "features": [],
            "source_label": label,
        }

    def _parts(self):
        # vector: [GGAG .. CGCT]   part: [CGCT .. GGAG]  -> closes a circle
        vec = self._frag("GGAG" + _filler(200) + "CGCT", "GGAG", "CGCT", "vec")
        part = self._frag("CGCT" + _filler(120, 3) + "GGAG", "CGCT", "GGAG",
                          "part")
        return vec, part

    def test_a_forward_part_still_assembles(self):
        vec, part = self._parts()
        got = _clo._gg_close_chain(vec, [part], [])
        assert got, "the ordinary forward assembly broke"

    def test_a_reverse_complemented_part_assembles_too(self):
        vec, part = self._parts()
        flipped = _clo._rc_fragment(part)
        got = _clo._gg_close_chain(vec, [flipped], [])
        assert got, "an RC'd part still reads as 'overhangs don't chain'"

    def test_the_flip_is_tagged_so_one_circle_counts_once(self):
        # The tags are what keep a circle found from two seeds from reading as
        # two products; without them each seed's freshly built flip has a new
        # id and the ambiguity warning fires spuriously.
        vec, part = self._parts()
        _clo._gg_close_chain(vec, [part], [])
        assert part["_gg_src_id"] == id(part)
        assert part["_gg_flipped"] is False

    def test_a_palindromic_fragment_is_not_offered_twice(self):
        # Offering a self-complementary fragment in both orientations would
        # double every solution count — ambiguity that isn't.
        vec, part = self._parts()
        sols = _clo._gg_close_chain(vec, [part], [])
        flipped_sols = _clo._gg_close_chain(vec, [_clo._rc_fragment(part)], [])
        assert len(sols) == len(flipped_sols) == 1


# ── CL5: a junction's two orientations are told apart ──────────────────────

def test_junction_warnings_name_their_orientation_when_they_differ():
    import inspect
    src = inspect.getsource(_clo._annotate_scars_on_product)
    # The classifications are collected per orientation and merged, rather than
    # each pass appending under the same label.
    assert "junction_cls" in src
    assert "(forward orientation)" in src and "(reverse orientation)" in src


class TestJunctionWarningsAreNotDuplicated:
    def _product(self):
        # A plain EcoRI/BamHI insert into an EcoRI/BamHI vector.
        vec_seq = _filler(400)
        ins_seq = _filler(200, 3)
        vec = {"top_seq": vec_seq,
               "left":  {"kind": "5'", "overhang_seq": "AATT",
                         "enzyme": "EcoRI"},
               "right": {"kind": "5'", "overhang_seq": "GATC",
                         "enzyme": "BamHI"},
               "features": [], "source_label": "vec"}
        ins = {"top_seq": ins_seq,
               "left":  {"kind": "5'", "overhang_seq": "GATC",
                         "enzyme": "BamHI"},
               "right": {"kind": "5'", "overhang_seq": "AATT",
                         "enzyme": "EcoRI"},
               "features": [], "source_label": "ins"}
        return vec, ins

    def test_no_junction_label_appears_twice_unqualified(self):
        vec, ins = self._product()
        result = _clo._simulate_traditional_cloning(ins, vec)
        _clo._annotate_scars_on_product(result, [ins], vec)
        lines = [w for w in result.get("warnings", [])
                 if w.startswith("Junction ")]
        bare = [ln for ln in lines if "orientation)" not in ln]
        assert len(bare) == len(set(bare)), (
            f"the same junction was reported twice unqualified: {bare}")


# ── CL6: an overlap longer than designed is reported ──────────────────────

def test_long_measured_overlap_is_reported_not_trimmed():
    import inspect
    src = inspect.getsource(
        sc.GibsonAssemblyPane._design_homology_arms)
    assert "_GIB_LONG_OVERLAP_SLACK_BP" in src, \
        "the designer no longer measures how far past min_overlap it matched"
    # And it must still NOT subtract the excess — that was tried and reverted
    # (`[project_audit_2026_07_18]`): designed homology and a coincidental
    # terminal identity are indistinguishable here, so trimming deletes bases.
    assert "up_seq[-min_oh:]" in src
    assert "found - min_oh" not in src and "- found" not in src


def test_the_long_overlap_warning_reaches_the_user():
    import inspect
    src = inspect.getsource(sc.GibsonAssemblyPane._on_design_arms)
    assert "long_overlaps" in src
    assert "collapses the WHOLE match" in src


def test_the_long_overlap_slack_is_a_named_constant():
    assert isinstance(_clo._GIB_LONG_OVERLAP_SLACK_BP, int)
    assert 0 < _clo._GIB_LONG_OVERLAP_SLACK_BP < _clo._GIBSON_MIN_OVERLAP_BP


# ── CL7: an ambiguous run is not a homology arm ───────────────────────────

class TestAmbiguityDoesNotAnneal:
    ARM = "ACGTTGCAAGGCTTACGTCAGG"          # 22 definite bases

    def test_a_real_arm_is_still_found(self):
        assert _clo._gibson_overlap_len("TTTT" + self.ARM,
                                        self.ARM + "AAAA") == len(self.ARM)

    def test_an_all_n_run_is_not_an_overlap(self):
        # Two runs of N compare EQUAL as strings and anneal not at all.
        assert _clo._gibson_overlap_len("TTTT" + "N" * 20,
                                        "N" * 20 + "AAAA") == 0

    def test_a_long_arm_tolerates_a_couple_of_ns(self):
        arm = self.ARM[:9] + "NN" + self.ARM[11:]
        assert _clo._gibson_overlap_len("TTTT" + arm, arm + "AAAA") == len(arm)

    def test_an_arm_with_too_few_definite_bases_is_refused(self):
        arm = self.ARM[:12] + "N" * 8       # 12 definite, floor is 15
        assert _clo._gibson_overlap_len("TTTT" + arm, arm + "AAAA") == 0

    def test_the_definite_base_count_is_what_it_says(self):
        assert _clo._gibson_overlap_definite_bases("ACGTN") == 4
        assert _clo._gibson_overlap_definite_bases("acgt") == 4
        assert _clo._gibson_overlap_definite_bases("NNNN") == 0
        assert _clo._gibson_overlap_definite_bases("ACRYWSN") == 2


# ── CL8: a junction-regenerated Type IIS site is reported ─────────────────

class TestClonedPlasmidJunctionSites:
    """The stub backbone is scrubbed and the insert is domesticated, so a site
    on the simulated product was formed at a junction — and the part would then
    be cut inside the very assembly it was built for. The designer refuses this;
    the stub-backbone simulation said nothing."""

    ENZ = ("BsaI", "Esp3I")

    def test_a_clean_part_reports_nothing(self):
        ins = _clean_of(self.ENZ, 300, 4)
        prod = _clo._simulate_cloned_plasmid(ins, "AATG", "GCTT", "CDS")
        assert _clo._cloned_plasmid_regenerated_sites(prod, ins, "AATG") == []

    def test_a_site_formed_at_the_overhang_boundary_is_reported(self):
        # oh5 "AATG" + an insert starting "GTCTC" spells GGTCTC (BsaI) across
        # the boundary; neither parent carries it.
        ins = "GTCTC" + _clean_of(self.ENZ, 295, 5)
        prod = _clo._simulate_cloned_plasmid(ins, "AATG", "GCTT", "")
        got = _clo._cloned_plasmid_regenerated_sites(prod, ins, "AATG")
        assert got and {d["enzyme"] for d in got} == {"BsaI"}

    def test_a_site_the_insert_already_carried_is_not_blamed_on_a_junction(self):
        # An undomesticated insert's OWN site must not be reported as
        # junction-formed — that would send the user looking at the overhangs.
        body = _clean_of(self.ENZ, 140, 6)
        ins = body + "GGTCTC" + _clean_of(self.ENZ, 140, 7)
        prod = _clo._simulate_cloned_plasmid(ins, "AATG", "GCTT", "")
        got = _clo._cloned_plasmid_regenerated_sites(prod, ins, "AATG")
        assert got == [], f"the insert's own site was blamed on a junction: {got}"


# ══════════════════════════════════════════════════════════════════════════
# Sequencing reads
# ══════════════════════════════════════════════════════════════════════════

def _rand_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# ── Partial reads: identity over what was COVERED ─────────────────────────

class TestPartialReadIdentity:
    """`identity_pct` is BLAST-style over the whole alignment, and in a
    semi-global alignment every base a partial read never reached is a gap
    column. So a PERFECT 400 bp read of a 1,200 bp reference scored 33% and the
    default `min_identity` of 99 called it a mismatch — the alignment was taught
    about partial reads (INV-206) and the metric was not. A Sanger read is the
    documented input, and it could never pass."""

    REF = _rand_seq(1200, 17)

    def _run(self, read, **kw):
        import splicecraft_agent as _agent
        return _agent._h_verify_against_reads(None, dict(
            {"reference": self.REF, "reads": [read], "circular": False}, **kw))

    def test_a_perfect_partial_read_matches(self):
        got = self._run(self.REF[400:800])
        r = got["reads"][0]
        assert r["covered_identity_pct"] == 100.0
        assert r["covered_bp"] == 400
        assert r["partial"] is True
        assert r["passes"] is True
        assert got["verdict"] == "match"

    def test_the_blast_style_figure_is_unchanged(self):
        # Existing callers read `identity_pct`; its meaning must not shift.
        r = self._run(self.REF[400:800])["reads"][0]
        assert r["identity_pct"] < 40.0
        assert r["coverage_pct"] == pytest.approx(33.33, abs=0.01)

    def test_a_passing_partial_read_says_it_verified_only_its_window(self):
        # "match" on 33% coverage is not a verified plasmid, and saying so is
        # the whole point — this is the most dangerous sentence the endpoint
        # can produce.
        got = self._run(self.REF[400:800])
        warn = " ".join(got.get("warnings", []))
        assert "partial" in warn and "COVERS" in warn
        assert "consensus" in warn

    def test_a_full_length_perfect_read_is_not_flagged_partial(self):
        got = self._run(self.REF)
        r = got["reads"][0]
        assert r["partial"] is False
        assert r["identity_pct"] == r["covered_identity_pct"] == 100.0
        assert not got.get("warnings")

    def test_a_divergent_full_length_read_still_fails(self):
        mut = list(self.REF)
        for i in range(0, len(mut), 40):
            mut[i] = {"A": "C", "C": "A", "G": "T", "T": "G"}[mut[i]]
        got = self._run("".join(mut))
        assert got["verdict"] == "mismatch"
        assert got["reads"][0]["passes"] is False

    def test_an_unrelated_partial_read_still_fails(self):
        got = self._run(_rand_seq(400, 99))
        assert got["verdict"] == "mismatch"
        assert got["reads"][0]["passes"] is False

    def test_covered_identity_is_iupac_aware_like_the_other_figure(self):
        aln = {"aligned_q": "ACGTN---", "aligned_t": "ACGTAGGG"}
        pct, covered = sc._covered_identity(aln)
        # The 3 unread columns don't count, and neither does the N: a base the
        # read could not call is not coverage (R10). It is not a difference
        # either, so the identity over what WAS read stays 100%.
        assert covered == 4
        assert pct == 100.0
        # A read base that rules the reference base OUT is still evidence.
        pct2, cov2 = sc._covered_identity({"aligned_q": "ACGTY",
                                           "aligned_t": "ACGTA"})
        assert (cov2, pct2) == (5, 80.0)

    def test_no_aligned_base_is_no_answer_not_agreement(self):
        assert sc._covered_identity({"aligned_q": "----",
                                     "aligned_t": "ACGT"}) == (0.0, 0)


# ── R7: indels group on their left-normalised position ────────────────────

class TestIndelsGroupAcrossReads:
    """An aligner may place the same 1 bp deletion anywhere inside a
    homopolymer. The rollup keyed on the raw alignment position, so one event
    became three, each with a single read behind it."""

    @staticmethod
    def _homopolymer_ref():
        return _rand_seq(120, 23) + "AAAAAAAA" + _rand_seq(120, 24)

    def test_one_deletion_three_placements_one_variant(self):
        import splicecraft_agent as _agent
        ref = self._homopolymer_ref()
        run_start = 120
        reads = [ref[:run_start + k * 3] + ref[run_start + k * 3 + 1:]
                 for k in range(3)]
        got = _agent._h_verify_against_reads(
            None, {"reference": ref, "reads": reads, "circular": False})
        vs = [v for v in (got.get("consensus") or {}).get("variants") or []
              if v.get("type") == "deletion"]
        assert len(vs) == 1, f"support split across {len(vs)} buckets"
        assert vs[0]["support"] == 3

    def test_normalisation_is_reported_separately_from_the_alignment_position(self):
        # `target_pos` must stay where the alignment put it — that is the column
        # the Phred lookup uses — while `norm_pos` is the canonical one.
        tgt = "GGGAAAAACCC"
        for gap_at in (3, 5, 7):
            aq = list(tgt)
            aq[gap_at] = "-"
            v = sc._extract_variants_from_alignment("".join(aq), tgt)[0]
            assert v["target_pos"] == gap_at
            assert v["norm_pos"] == 3
            assert v["norm_ref"] == "A"

    def test_insertions_normalise_and_rotate(self):
        tgt = "GGGAAAAACCC"
        for ins_at in (3, 5, 7):
            at = tgt[:ins_at] + "-" + tgt[ins_at:]
            aq = tgt[:ins_at] + "A" + tgt[ins_at:]
            v = sc._extract_variants_from_alignment(aq, at)[0]
            assert v["norm_pos"] == 3 and v["norm_alt"] == "A"

    def test_a_substitution_gets_no_normalised_position(self):
        v = sc._extract_variants_from_alignment("GGGAAATACCC",
                                               "GGGAAAAACCC")[0]
        assert v["type"] == "snp"
        assert "norm_pos" not in v

    def test_the_group_key_falls_back_when_there_is_no_norm_pos(self):
        k = sc._variant_group_key({"type": "snp", "target_pos": 7,
                                   "ref": "A", "alt": "T"})
        assert k == ("snp", 7, "A", "T")


# ── R10: an ambiguity code is a no-call, not a change ─────────────────────

class TestAmbiguityIsNotAChange:
    """The segment renderer counted `N` as a match (via `_iupac_compatible`)
    while the trace tally counted it as a real change — the same alignment read
    two ways on two screens."""

    def test_the_variant_carries_its_compatibility(self):
        v = sc._extract_variants_from_alignment("GGGAANAACCC",
                                               "GGGAAAAACCC")[0]
        assert v["type"] == "snp" and v["iupac_compatible"] is True
        v2 = sc._extract_variants_from_alignment("GGGAATAACCC",
                                                "GGGAAAAACCC")[0]
        assert v2["iupac_compatible"] is False

    def test_an_incompatible_ambiguity_pair_is_still_a_change(self):
        v = sc._extract_variants_from_alignment("GGGAARAACCC",
                                               "GGGAACAACCC")[0]
        assert v["iupac_compatible"] is False     # R = A/G, never C

    def test_the_trace_verdict_does_not_call_an_n_a_real_change(self):
        ref = _rand_seq(200, 31)
        read = list(ref)
        read[100] = "N"
        summary = sc._trace_verification_summary(
            "".join(read), ref, [50] * len(ref), min_phred=20, circular=True)
        assert summary["n_snp"] == 0, "an N counted as a SNP"
        assert summary["n_no_call"] == 1, "the N was silently dropped instead"
        assert summary["verdict"] == "clean"

    def test_a_real_change_still_reaches_the_verdict(self):
        ref = _rand_seq(200, 32)
        read = list(ref)
        read[100] = {"A": "C", "C": "A", "G": "T", "T": "G"}[read[100]]
        summary = sc._trace_verification_summary(
            "".join(read), ref, [50] * len(ref), min_phred=20, circular=True)
        assert summary["n_snp"] == 1
        assert summary["verdict"] == "real_changes"

    def test_the_two_halves_now_agree_about_the_same_alignment(self):
        # The map's segment renderer and the variant tally, on one alignment.
        ref = _rand_seq(160, 33)
        read = list(ref)
        read[80] = "N"
        segs = sc._alignment_to_target_segments("".join(read), ref)
        assert not [s for s in segs if s[2] == "mismatch"]
        summary = sc._trace_verification_summary(
            "".join(read), ref, [50] * len(ref), min_phred=20, circular=True)
        assert summary["n_snp"] == 0


# ── R13: inversion detection must not depend on the stored rotation ───────

def test_the_wrap_edge_allowance_scales_with_the_block():
    # A fixed column count made the same inversion detectable at one origin and
    # invisible at another. The allowance is now a floor PLUS a fraction of the
    # joined block's own length.
    assert sc._INVERSION_WRAP_EDGE_SLACK > 0
    assert 0.0 < sc._INVERSION_WRAP_EDGE_FRACTION < 0.5
    import inspect
    src = inspect.getsource(sc._alignment_inverted_segments)
    assert "_INVERSION_WRAP_EDGE_FRACTION" in src
    # Which candidates the wrap join considers is pinned BEHAVIOURALLY, at 30
    # rotations, by `TestOriginInversionsAreRotationInvariant`.


# ── R14: a per-base quality array must cover the read ────────────────────

class TestReadQualityMustBeParallelToItsRead:
    REF = _rand_seq(300, 41)

    def test_a_short_quality_array_is_refused(self):
        import splicecraft_agent as _agent
        read = self.REF[:200]
        got = _agent._h_verify_against_reads(None, {
            "reference": self.REF, "reads": [read], "circular": False,
            "read_quality": [[50] * 100]})       # half the read
        assert isinstance(got, tuple) and got[1] == 400
        assert "one per base" in got[0]["error"]

    def test_a_correct_array_is_accepted(self):
        import splicecraft_agent as _agent
        read = self.REF[:200]
        got = _agent._h_verify_against_reads(None, {
            "reference": self.REF, "reads": [read], "circular": False,
            "read_quality": [[50] * 200]})
        assert not isinstance(got, tuple)

    def test_null_is_still_allowed(self):
        import splicecraft_agent as _agent
        got = _agent._h_verify_against_reads(None, {
            "reference": self.REF, "reads": [self.REF[:200]],
            "circular": False, "read_quality": [None]})
        assert not isinstance(got, tuple)

    def test_heterogeneity_checks_it_too(self):
        import splicecraft_agent as _agent
        got = _agent._h_analyse_read_heterogeneity(None, {
            "reference": self.REF, "reads": [self.REF[:200]] * 3,
            "circular": False, "read_quality": [[50] * 200, [50] * 5,
                                                [50] * 200]})
        assert isinstance(got, tuple) and got[1] == 400
        assert "one per base" in got[0]["error"]


def test_capped_by_is_none_when_nothing_was_capped():
    import splicecraft_agent as _agent
    ref = _rand_seq(400, 43)
    got = _agent._h_analyse_read_heterogeneity(None, {
        "reference": ref, "reads": [ref] * 3, "circular": False,
        "max_reads": 50})
    # The caller PASSED max_reads, but all three reads were used — reporting a
    # cap told them their fractions came from a truncated set.
    assert got["capped_by"] is None
    assert got["reads_available"] == 3


# ── R14b/d: a read cap is named whichever limit set it, and a FASTQ counts ──

class TestReadCapsAreNamed:
    """The default `max_reads` capped 450 inline reads to 400 while reporting
    `capped_by: None`; and the FASTQ loader stops AT `max_reads`, so a 1,500
    read file and a 400 read file both reported 400 available, no cap and no
    warning — a silent truncation under every fraction."""

    @staticmethod
    def _fastq(tmp_path, ref, n):
        p = tmp_path / "reads.fastq"
        p.write_text("".join(f"@r{i}\n{ref}\n+\n{'I' * len(ref)}\n"
                             for i in range(n)))
        return str(p)

    def test_the_default_max_reads_is_named_when_it_caps(self, monkeypatch):
        import splicecraft_agent as _agent
        monkeypatch.setattr(_agent, "_HETERO_READS_MAX", 5)
        ref = _rand_seq(400, 44)
        got = _agent._h_analyse_read_heterogeneity(None, {
            "reference": ref, "reads": [ref] * 8, "circular": False})
        assert got["reads_available"] == 8 and got["max_reads"] == 5
        assert got["capped_by"] == "max_reads"

    def test_a_fastq_larger_than_the_sample_says_so(self, tmp_path):
        import splicecraft_agent as _agent
        ref = _rand_seq(400, 45)
        got = _agent._h_analyse_read_heterogeneity(None, {
            "reference": ref, "reads_path": self._fastq(tmp_path, ref, 12),
            "circular": False, "max_reads": 5})
        assert got["reads_available"] == 12 and got["capped_by"] == "max_reads"
        assert any("only 5 were used" in w for w in got["warnings"])

    def test_a_fastq_of_exactly_the_sample_is_not_capped(self, tmp_path):
        import splicecraft_agent as _agent
        ref = _rand_seq(400, 46)
        got = _agent._h_analyse_read_heterogeneity(None, {
            "reference": ref, "reads_path": self._fastq(tmp_path, ref, 5),
            "circular": False, "max_reads": 5})
        assert got["reads_available"] == 5 and got["capped_by"] is None


# ── R14c: a heterogeneity mode needs calls at or above the floor ───────────

class TestAModeNeedsEvidenceAboveTheFloor:
    """The shape test judged a bin by its CENTRE. With a floor inside a bin,
    calls all below the floor made a mode (and a `mixed` verdict), and a real
    population just above it was skipped."""

    def test_sub_floor_calls_sharing_a_bin_are_not_a_mode(self):
        import splicecraft_agent as _agent
        fr = [0.051, 0.055, 0.058, 0.06, 0.062, 0.065] + [0.01] * 2
        assert _agent._hetero_shape(fr, 0.07)["mode_fraction"] is None

    def test_a_population_just_above_the_floor_is_a_mode(self):
        import splicecraft_agent as _agent
        out = _agent._hetero_shape([0.14] * 6 + [0.01] * 5, 0.13)
        assert out["mode_fraction"] is not None
        assert out["mode_fraction"] >= 0.13


# ══════════════════════════════════════════════════════════════════════════
# Data safety
# ══════════════════════════════════════════════════════════════════════════

# ── D11: a recovery is reported whoever loaded the file first ─────────────

class TestRecoveryIsRememberedNotJustReturned:
    """The startup check reports recoveries by loading each file again — but a
    recovery REWRITES the main file, so a load that happened earlier in the
    launch swallowed the notice and the user was never told their library had
    been corrupt and restored."""

    def test_a_recovery_is_recorded_at_the_point_it_happens(self, tmp_path):
        import json
        import splicecraft_persistence as _p
        sc._state._DATA_RECOVERIES.clear()
        main = tmp_path / "plasmid_library.json"
        bak = main.with_suffix(main.suffix + ".bak")
        bak.write_text(json.dumps({"_schema_version": 1,
                                   "entries": [{"id": "x", "name": "Kept"}]}))
        main.write_text("{ not json")
        entries, warning = _p._safe_load_json(main, "Plasmid library")
        assert len(entries) == 1 and warning
        assert [lbl for lbl, _m in sc._state._DATA_RECOVERIES] \
            == ["Plasmid library"]
        sc._state._DATA_RECOVERIES.clear()

    def test_the_log_is_bounded_and_deduplicated(self):
        import splicecraft_persistence as _p
        sc._state._DATA_RECOVERIES.clear()
        for _ in range(200):
            _p._note_data_recovery("Plasmid library", "restored")
        assert len(sc._state._DATA_RECOVERIES) == 1
        for i in range(200):
            _p._note_data_recovery(f"label {i}", "restored")
        assert len(sc._state._DATA_RECOVERIES) <= 64
        sc._state._DATA_RECOVERIES.clear()

    def test_the_startup_check_reports_a_recovery_it_did_not_see(self):
        import inspect
        src = inspect.getsource(sc.PlasmidApp._check_data_files)
        assert "_DATA_RECOVERIES" in src, \
            "the startup check no longer reports recoveries from earlier loads"


# ── D12: a read-only launch does not write ────────────────────────────────

class TestReadOnlyLaunchDoesNotWrite:
    def test_the_one_shot_migrations_are_skipped(self):
        import inspect
        src = inspect.getsource(sc.PlasmidApp.on_mount)
        i = src.index("_migrate_parts_bin_markers_from_vector")
        # The read-only gate must come BEFORE the migration calls.
        assert "_AGENT_READ_ONLY" in src[:i], \
            "the writing migrations are no longer gated on read-only mode"

    def test_the_crash_recovery_prune_survives_a_refused_delete(self):
        # The chokepoint raises RuntimeError, which is CORRECT in read-only
        # mode. Catching only OSError let it escape into on_mount.
        import inspect
        src = inspect.getsource(sc.PlasmidApp._check_crash_recovery)
        assert "except OSError:\n                    pass" not in src
        assert src.count("except (OSError, RuntimeError)") >= 3


# ── D10: a library that belongs to no collection is kept, not overwritten ──
#
# The first fix marked the mirror dirty when the pointer was cleared. That made
# the next launch KEEP the orphaned library — and then flush it into the
# collection it had just made active, replacing that collection's plasmids.
# These pin the replacement: the orphans are rescued into their own
# collection, and no existing collection is touched.

def _orphan_launch(monkeypatch):
    import json
    import splicecraft_state as st
    monkeypatch.setattr(sc, "_collection_sync_force_sync", True)
    monkeypatch.setenv("SPLICECRAFT_SKIP_SETTINGS_FLUSH", "1")
    sc._save_collections([
        {"name": "Main", "plasmids": [_d_ent("M1", 31)]},
        {"name": "Other", "plasmids": [_d_ent("O1", 32)]}])
    sc._set_active_collection_name("")               # no active collection
    sc._safe_save_json_mirror(st._LIBRARY_FILE,
                              [_d_ent("ORPHAN", 33), _d_ent("M1", 31)],
                              "Plasmid library")
    for attr in ("_library_cache", "_collections_cache", "_settings_cache"):
        setattr(st, attr, None)
    sc._ensure_default_collection()                  # the launch sequence
    sc._restore_library_from_active_collection()
    return json


def test_a_launch_rescues_the_orphans_into_their_own_collection(monkeypatch):
    _orphan_launch(monkeypatch)
    colls = _colls_by_name()
    assert colls["Recovered plasmids"] == ["ORPHAN"]
    assert colls["Main"] == ["M1"] and colls["Other"] == ["O1"]


def test_a_switch_after_that_does_not_overwrite_any_collection(monkeypatch):
    _orphan_launch(monkeypatch)
    sc._activate_collection("Other")                 # the next switch
    assert _colls_by_name()["Main"] == ["M1"]
    assert _colls_by_name()["Other"] == ["O1"]


def test_deleting_the_active_collection_through_the_agent_switches(
        monkeypatch):
    import splicecraft_state as st
    monkeypatch.setattr(sc, "_collection_sync_force_sync", True)
    sc._save_collections([{"name": "Keep", "plasmids": [_d_ent("K1", 41)]},
                          {"name": "Scratch", "plasmids": [_d_ent("S1", 42)]}])
    sc._activate_collection("Scratch")
    r = st._AGENT_HANDLERS["delete-collection"][0](_MirrorRig(),
                                                   {"name": "Scratch"})
    assert r["active"] == "Keep"
    assert sc._get_active_collection_name() == "Keep"
    assert [e["name"] for e in sc._load_library()] == ["K1"]
    assert "Recovered plasmids" not in _colls_by_name()   # a delete is a delete


# ── D13: backups are ordered by generation, never by the clock ─────────────

class _Clock:
    """A stand-in for `datetime` in the persistence module whose `now()` the
    test sets."""
    def __init__(self):
        from datetime import datetime
        self._real = datetime
        self.t = datetime(2026, 9, 23, 1, 30, 0)

    def now(self, *a, **k):
        return self.t

    def __getattr__(self, name):
        return getattr(self._real, name)


class TestBackupsOrderByGeneration:
    """Retention keeps "the newest N" backups. Ordering them by the clock
    stamp in the name went wrong three ways: a busy second reused the slot
    retention had just freed (the new backup sorted OLDEST and was pruned at
    once), a clock that stepped back filed new backups under old ones, and a
    clock that ran fast and was corrected kept the future-stamped ones and
    pruned everything written after the correction. Each check reads what the
    surviving backup FILES contain."""

    @pytest.fixture
    def rig(self, monkeypatch):
        import json
        import splicecraft_persistence as pe
        import splicecraft_state as st
        clock = _Clock()
        monkeypatch.setattr(pe, "_datetime", clock)
        path = st._DATA_DIR / "gen_probe.json"

        def save(i):
            pe._safe_save_json(path, [{"id": f"p{i}", "v": i}], "Gen probe")

        def kept():
            out = []
            for bak in pe._iter_backups(path):
                raw = json.loads(pe._read_backup_bytes(bak).decode())
                out.append(raw["entries"][0]["v"])
            return out

        return clock, save, kept, monkeypatch

    def test_a_busy_second_keeps_the_newest(self, rig):
        import splicecraft_state as st
        clock, save, kept, mp = rig
        mp.setattr(st, "_BACKUP_RETENTION_COUNT", 10)
        for i in range(16):                       # all in one wall-second
            save(i)
        # 15 prior states were backed up (payloads 0..14); keep the last 10.
        assert kept() == list(range(5, 15))

    def test_a_clock_that_steps_back(self, rig):
        from datetime import timedelta
        clock, save, kept, _mp = rig
        for i in range(4):
            save(i)
            clock.t += timedelta(seconds=5)
        clock.t -= timedelta(hours=1)            # DST fall-back / NTP step
        for i in range(4, 7):
            save(i)
            clock.t += timedelta(seconds=5)
        assert kept()[-1] == 5                   # the newest IS the newest
        assert kept() == sorted(kept())

    def test_a_fast_clock_that_was_corrected(self, rig):
        import splicecraft_state as st
        from datetime import timedelta
        clock, save, kept, mp = rig
        mp.setattr(st, "_BACKUP_RETENTION_COUNT", 3)
        real = clock.t
        clock.t = real + timedelta(days=365)     # a year fast
        for i in range(3):
            save(i)
            clock.t += timedelta(seconds=5)
        clock.t = real                           # corrected
        for i in range(3, 7):
            save(i)
            clock.t += timedelta(seconds=5)
        assert kept() == [3, 4, 5]

    def test_older_names_sort_before_every_generation(self):
        import splicecraft_persistence as pe
        names = ["x.json.bak.20990101-000000.g000002",
                 "x.json.bak.20260101-000000.012",
                 "x.json.bak.20250101-000000.g000001.gz",
                 "x.json.bak.20270101-000000",
                 "x.json.bak.20260101-000000.2"]
        assert sorted(names, key=pe._backup_sort_key) == [
            "x.json.bak.20260101-000000.2",
            "x.json.bak.20260101-000000.012",
            "x.json.bak.20270101-000000",
            "x.json.bak.20250101-000000.g000001.gz",
            "x.json.bak.20990101-000000.g000002"]


# ── D15: a Master Delete leaves nothing recoverable in memory ─────────────

def test_master_delete_clears_the_library_undo_buffer():
    app = sc.PlasmidApp.__new__(sc.PlasmidApp)
    app._library_undo_stack().append({"name": "Deleted plasmid",
                                      "gb_text": "LOCUS x"})
    assert app._library_undo_stack()
    sc._reset_app_state_after_master_delete(app)
    assert app._library_undo_stack() == [], \
        "the entries the user asked to destroy were one keypress from returning"


def test_the_lock_does_not_unlink_a_file_another_process_holds():
    import inspect
    src = inspect.getsource(sc._acquire_data_dir_lock)
    i = src.index("lock.contended")
    before = src[:i]
    assert "lockfile.unlink" not in before, \
        "unlinking on contention races the winner's inode — two instances can " \
        "then both hold the data dir"


# ── D9: a committed move is not reported as failed ───────────────────────

# (The committed-move mirror failure is pinned BEHAVIOURALLY by
# `TestACommittedMoveSurvivesAMirrorFailure` below — a source-text check here
# asserted the backwards `_mark_mirror_dirty` recovery that test caught.)


# ── D3: restoring settings re-stages the library it re-points ────────────

def test_restoring_settings_restages_the_library():
    import inspect
    src = inspect.getsource(sc._reconcile_mirror_after_restore)
    assert '_is("_SETTINGS_FILE")' in src, \
        "a settings restore can move the active-collection pointer with the " \
        "live mirror still holding the old collection's plasmids"
    assert "_mark_mirror_dirty" in src


# ══════════════════════════════════════════════════════════════════════════
# Formats, agent API, docs
# ══════════════════════════════════════════════════════════════════════════

def _prot_record(loc, seq="MKVLAAGTRWQ" * 4):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature
    r = SeqRecord(Seq(seq), id="PROT", name="PROT", description="p")
    r.annotations["molecule_type"] = "protein"
    r.features = [SeqFeature(loc, type="Region",
                             qualifiers={"note": ["domain"]})]
    return r


# ── FM9: a protein record exports whatever strand its features carry ──────

class TestProteinExportIgnoresStrand:
    """A GenBank PROTEIN file has no strand column, so BioPython reads every
    protein feature back strandless. Comparing strand in the round-trip guard
    refused every protein record whose features carried one — and, once the
    strand was stripped for that reason, refused them for lacking one."""

    def _locs(self):
        from Bio.SeqFeature import FeatureLocation, CompoundLocation
        return [
            ("order", CompoundLocation([FeatureLocation(0, 5, strand=1),
                                        FeatureLocation(10, 15, strand=1)],
                                       operator="order")),
            ("join", CompoundLocation([FeatureLocation(0, 5, strand=1),
                                       FeatureLocation(10, 15, strand=1)])),
            ("plain forward", FeatureLocation(0, 20, strand=1)),
            ("strandless", FeatureLocation(0, 20)),
        ]

    def test_every_strand_shape_exports(self, tmp_path):
        for label, loc in self._locs():
            out = tmp_path / f"{label.replace(' ', '_')}.gb"
            res = sc._export_genbank_to_path(_prot_record(loc), out)
            assert res["path"] == str(out), label
            assert out.read_text().startswith("LOCUS"), label

    @pytest.mark.parametrize("op", ["order", "bond"])
    def test_a_compound_operator_survives_export(self, tmp_path, op):
        # Biopython writes and re-reads `order(...)` / `bond(...)` verbatim;
        # the strand-stripping rebuild used to drop the operator, so every one
        # came out as `join(...)`. Checked against Biopython's own parse of
        # the file, with no warning claiming a loss that no longer happens.
        from Bio import SeqIO
        from Bio.SeqFeature import FeatureLocation, CompoundLocation
        loc = CompoundLocation([FeatureLocation(0, 5, strand=1),
                                FeatureLocation(10, 15, strand=1)],
                               operator=op)
        out = tmp_path / "o.gb"
        res = sc._export_genbank_to_path(_prot_record(loc), out)
        assert not res.get("warnings")
        back = SeqIO.read(str(out), "genbank").features
        assert [f.location.operator for f in back if f.type == "Region"] \
            == [op]

    def test_an_arrowless_order_feature_keeps_its_operator(self, tmp_path):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        r = SeqRecord(Seq("ACGT" * 40), id="N", name="N")
        r.annotations.update({"molecule_type": "DNA", "topology": "circular"})
        r.features = [SeqFeature(
            CompoundLocation([FeatureLocation(0, 5, strand=0),
                              FeatureLocation(10, 15, strand=0)],
                             operator="order"),
            type="misc_feature", qualifiers={"label": ["o"]})]
        out = tmp_path / "n.gb"
        sc._export_genbank_to_path(r, out)
        assert "order(" in out.read_text()
        back = [f for f in SeqIO.read(str(out), "genbank").features
                if f.type == "misc_feature"]
        assert back[0].location.operator == "order"

    def test_a_nucleotide_record_still_keeps_its_strand(self, tmp_path):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        r = SeqRecord(Seq("ACGT" * 40), id="N", name="N")
        r.annotations.update({"molecule_type": "DNA", "topology": "circular"})
        r.features = [SeqFeature(FeatureLocation(0, 20, strand=-1),
                                 type="CDS", qualifiers={"note": ["rev"]})]
        out = tmp_path / "n.gb"
        sc._export_genbank_to_path(r, out)
        assert "complement(" in out.read_text()


# ── L5b: the exported topology is the one the app believes ────────────────

class TestExportedTopologyMatchesTheApp:
    """`_record_is_circular` says an unannotated record is a CIRCLE and the map
    draws one; the exporter defaulted the same record to `linear`, so every
    other tool read the user's plasmid as a fragment."""

    @staticmethod
    def _rec(topology=None):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        r = SeqRecord(Seq("ACGT" * 40), id="P", name="P")
        r.annotations["molecule_type"] = "DNA"
        if topology:
            r.annotations["topology"] = topology
        return r

    @pytest.mark.parametrize("topology", [None, "linear", "circular"])
    def test_the_two_halves_agree(self, topology):
        rec = self._rec(topology)
        exported = sc._normalize_for_genbank(rec).annotations["topology"]
        app_says = "circular" if sc._record_is_circular(rec) else "linear"
        assert exported == app_says, (
            f"the app says {app_says} and the file says {exported}")

    def test_an_explicit_topology_is_never_overridden(self):
        for topology in ("linear", "circular"):
            rec = self._rec(topology)
            assert sc._normalize_for_genbank(
                rec).annotations["topology"] == topology

    def test_a_fasta_import_is_still_linear(self, tmp_path):
        # The reason the export used to default to `linear` — handled at INGEST
        # now, which is why the export can follow the app.
        fa = tmp_path / "x.fa"
        fa.write_text(">frag\n" + "ACGT" * 20 + "\n")
        rec = sc._fasta_path_to_record(str(fa))
        assert rec.annotations.get("topology") == "linear"
        assert sc._normalize_for_genbank(
            rec).annotations["topology"] == "linear"


# ── FM7: a long lineage still parses ──────────────────────────────────────

class TestDeepHistoryStillParses:
    """A `.dna` construction history nests one level per save, and the XML
    guard capped nesting at 256 — so the provenance of the most-worked-on
    plasmids was exactly what stopped loading."""

    @staticmethod
    def _nested(n):
        return ("<Hist>" + "".join(f'<Node ID="{i}">' for i in range(n))
                + "".join("</Node>" for _ in range(n)) + "</Hist>")

    @pytest.mark.parametrize("depth", [200, 260, 400, 500])
    def test_a_realistic_lineage_parses(self, depth):
        assert sc._safe_xml_parse(self._nested(depth)) is not None

    def test_the_guard_still_bounds(self):
        import xml.etree.ElementTree as ET
        with pytest.raises(ET.ParseError):
            sc._safe_xml_parse(self._nested(5000))

    def test_the_ceiling_clears_the_history_model(self):
        # The two caps must not disagree: the model accepts 500 levels, so the
        # parser has to reach at least that far.
        assert sc._safe_xml_parse(self._nested(
            sc._HISTORY_NODE_MAX_DEPTH)) is not None


# ── FM11: a corrupt embedded XML is an error, not a traceback ─────────────

def test_load_genbank_turns_a_malformed_xml_into_a_value_error():
    import inspect
    src = inspect.getsource(sc.load_genbank)
    assert "_ExpatError" in src and "RecursionError" in src, \
        "a corrupt .dna history escapes load_genbank as a raw traceback"


# ── FM10 / L10: the exported file carries no control bytes, no locale ─────

class TestExportedTextIsCleanAndLocaleFree:
    @staticmethod
    def _rec_with(note):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        r = SeqRecord(Seq("ACGT" * 40), id="P", name="P")
        r.annotations.update({"molecule_type": "DNA", "topology": "circular"})
        r.features = [SeqFeature(FeatureLocation(0, 40, strand=1),
                                 type="misc_feature",
                                 qualifiers={"note": [note]})]
        return r

    def test_a_lone_cr_or_tab_never_reaches_the_file(self, tmp_path):
        out = tmp_path / "x.gb"
        sc._export_genbank_to_path(
            self._rec_with("before\rafter\ttab and \x0b vtab"), out)
        text = out.read_text()
        assert "\r" not in text and "\t" not in text and "\x0b" not in text

    def test_the_locus_date_comes_from_a_fixed_table(self):
        import re as _re
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        r = SeqRecord(Seq("ACGT" * 20), id="P", name="P")
        r.annotations["molecule_type"] = "DNA"
        date = sc._normalize_for_genbank(r).annotations["date"]
        assert _re.fullmatch(
            r"\d{2}-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-\d{4}",
            date), date

    def test_it_does_not_call_strftime_for_the_month(self):
        # `%b` is locale-dependent and the flat-file format is not.
        import inspect
        import splicecraft_fileio as _fio
        src = inspect.getsource(_fio._normalize_for_genbank)
        assert 'strftime("%d-%b-%Y")' not in src, \
            "the LOCUS date is built with a locale-dependent month again"
        assert "_INSDC_MONTHS" in src


# ══════════════════════════════════════════════════════════════════════════
# Agent API
# ══════════════════════════════════════════════════════════════════════════

class TestDigestInputHandling:
    SEQ = "GAATTCAAAAGGATCC" * 5

    def _digest(self, enzymes):
        import splicecraft_agent as _agent
        return _agent._h_digest(None, {"sequence": self.SEQ,
                                       "enzymes": enzymes, "circular": False})

    def test_a_comma_separated_string_is_a_list(self):
        # A GET query string has no lists, so `?enzymes=EcoRI,BamHI` is the
        # natural spelling — and taking it as ONE name digested with nothing.
        got = self._digest("EcoRI,BamHI")
        assert got["n_cuts"] == 10
        assert got["unknown_enzymes"] == []

    def test_no_catalog_name_contains_a_comma(self):
        # What makes the split unambiguous.
        assert not [n for n in sc._all_enzymes() if "," in n]

    def test_an_unknown_name_is_still_reported(self):
        # Beside a name that resolves it is reported; ALONE there is nothing
        # to digest with, which is a 400 (see `TestDigestNeedsAnEnzyme`).
        got = self._digest(["EcoRI", "NotAnEnzyme"])
        assert got["unknown_enzymes"] == ["NotAnEnzyme"]
        alone = self._digest(["NotAnEnzyme"])
        assert isinstance(alone, tuple) and alone[1] == 400
        assert alone[0]["unknown_enzymes"] == ["NotAnEnzyme"]


class TestIdempotencyKeyIsBoundToTheBody:
    """A key replayed with a DIFFERENT body was answered with the first call's
    cached response, so the second operation never ran while the caller was
    told it had — the one failure an idempotency key exists to prevent."""

    def test_the_same_body_replays(self):
        sc._agent_idempotency_put("e", "k1", {"ok": True}, 200, "hashA")
        assert sc._agent_idempotency_get("e", "k1", "hashA") == ({"ok": True},
                                                                 200)

    def test_a_different_body_conflicts(self):
        sc._agent_idempotency_put("e", "k2", {"ok": True}, 200, "hashA")
        assert sc._agent_idempotency_get("e", "k2", "hashB") == "__conflict__"

    def test_no_hash_still_replays_for_older_callers(self):
        sc._agent_idempotency_put("e", "k3", {"ok": True}, 200, "hashA")
        assert sc._agent_idempotency_get("e", "k3", None) == ({"ok": True}, 200)


class TestEnvelopesReflectTheOutcome:
    """`ok` is what HAPPENED, not "the request parsed". Each of these answered
    `ok: true` beside a result that said nothing was done, and a caller
    checking the envelope carried on. (The domesticate-parts batch is covered
    behaviourally in test_agent_api.TestDomesticatePartsBatch.)"""

    @staticmethod
    def _h(name):
        return sc._state._AGENT_HANDLERS[name][0]

    @staticmethod
    def _gb_file(path, name="pGood"):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq("ATGC" * 60), id=name, name=name,
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        path.write_text(sc._record_to_gb_text(rec))

    def _import(self, folder, coll):
        return self._h("bulk-import-folder")(
            _MirrorRig(), {"folder": str(folder), "collection": coll})

    def _names(self):
        return [c["name"] for c in sc._load_collections()]

    # ── bulk-import-folder ──
    def test_an_import_that_read_nothing_creates_nothing(self, tmp_path):
        bad = tmp_path / "bad"
        bad.mkdir()
        (bad / "broken.gb").write_text("this is not a GenBank file\n")
        r = self._import(bad, "Batch")
        assert isinstance(r, tuple) and r[1] == 422
        assert r[0]["ok"] is False and r[0]["n_imported"] == 0
        assert r[0]["failures"] and "broken.gb" in r[0]["failures"][0]["path"]
        assert "Batch" not in self._names()

    def test_the_same_call_works_once_the_folder_is_fixed(self, tmp_path):
        folder = tmp_path / "f"
        folder.mkdir()
        (folder / "broken.gb").write_text("garbage\n")
        assert self._import(folder, "Batch")[1] == 422
        (folder / "broken.gb").unlink()
        self._gb_file(folder / "good.gb")
        r = self._import(folder, "Batch")
        assert not isinstance(r, tuple), r            # was 409 "already exists"
        assert r["ok"] is True and r["n_imported"] == 1
        assert "Batch" in self._names()

    def test_an_empty_folder_creates_nothing(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        r = self._import(empty, "Nothing")
        assert isinstance(r, tuple) and r[1] == 422
        assert "Nothing" not in self._names()

    def test_a_partial_import_keeps_what_it_read_and_says_so(self, tmp_path):
        folder = tmp_path / "mix"
        folder.mkdir()
        self._gb_file(folder / "good.gb")
        (folder / "bad.gb").write_text("garbage\n")
        r = self._import(folder, "Mixed")
        # A success with failures (harden 2026-09-24): `ok: false` here dropped
        # `data` and made a caller retry into 409 "already exists".
        assert not isinstance(r, tuple), r
        assert r["ok"] is True and r["partial"] is True and "error" not in r
        assert r["warnings"] and r["failures"]
        assert r["n_imported"] == 1 and r["n_failed"] == 1
        assert "Mixed" in self._names()
        assert sc._agent_data_envelope(r, 200).get("data")

    # ── ot2-compile ──
    @staticmethod
    def _plan(src_well="A1"):
        return {"pipette": "p300_single", "mount": "left",
                "tips": {"labware": "tiprack_300", "slot": 8},
                "labware": {"src": {"labware": "eppi_24", "slot": 1},
                            "dst": {"labware": "plate_24", "slot": 2}},
                "transfers": [{"from": f"src:{src_well}", "to": "dst:A1",
                               "volume": 50}]}

    def test_an_invalid_plan_is_not_ok(self):
        r, code = self._h("ot2-compile")(None, self._plan("Z9"))
        assert code == 422          # the CLI exits 0 on any 2xx
        assert r["valid"] is False and "protocol" not in r
        assert r["ok"] is False and r.get("error")

    def test_a_valid_plan_is_ok(self):
        r = self._h("ot2-compile")(None, self._plan("A1"))
        assert r["valid"] is True and r["ok"] is True and "protocol" in r
        assert "error" not in r

    # ── delete-primers ──
    def test_deleting_nothing_is_not_ok(self):
        sc._save_primers([])
        self._h("create-primer")(None, {"name": "dup", "sequence": "AAAACCCCGG"})
        self._h("create-primer")(None, {"name": "dup", "sequence": "TTTTGGGGAA"})
        r, code = self._h("delete-primers")(None, {"names": ["dup", "nope"]})
        # 409 when a name was ambiguous — those primers exist (the single
        # `delete-primer` answers 409 too); 404 only when none was found
        # (hardening 2026-09-25).
        assert code == 409
        assert r["removed"] == 0
        assert r["ok"] is False
        assert "1 not found" in r["error"] and "1 ambiguous" in r["error"]
        assert len(sc._load_primers()) == 2

    def test_a_partial_prune_is_ok_and_lists_what_it_skipped(self):
        sc._save_primers([])
        self._h("create-primer")(None, {"name": "p1", "sequence": "ACGTACGTAA"})
        r = self._h("delete-primers")(None, {"names": ["p1", "nope"]})
        assert r["ok"] is True and r["removed"] == 1
        assert r["not_found"] == ["nope"] and "error" not in r


# ── AA20b: a flag sent as TEXT reads the right way round ──────────────────

class TestTextFlagsReadTheRightWayRound:
    """A GET query string, shell-built JSON and small models all send booleans
    as text, and `bool("false")` is True. These sites still used it, so
    `apply: "false"` WROTE the annotations it was asked only to preview and
    `circular: "false"` analysed a linear molecule as a circle."""

    class _PresetApp:
        _record_load_counter = 0
        _unsaved = False

        def __init__(self, rec):
            self._current_record = rec
            self.applied = []

        def call_from_thread(self, fn, *a, **k):
            return fn(*a, **k)

        def _apply_annotation_transfers(self, transfers, rec):
            self.applied.append(len(transfers))

    @staticmethod
    def _h(name):
        return sc._state._AGENT_HANDLERS[name][0]

    def _app(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        import splicecraft_presets as P
        t7 = next(f for f in P._preset_features()
                  if f["name"] == "T7 promoter")["sequence"]
        seq = _rand_seq(300, 5) + t7 + _rand_seq(300, 6)
        rec = SeqRecord(Seq(seq), id="t", name="t",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        return self._PresetApp(rec)

    def test_apply_false_as_text_only_previews(self):
        app = self._app()
        r = self._h("annotate-from-presets")(app, {"apply": "false"})
        assert r["count"] >= 1, "probe: the planted T7 promoter was not found"
        assert r["applied"] is False and app.applied == []

    def test_apply_true_as_text_applies(self):
        app = self._app()
        r = self._h("annotate-from-presets")(app, {"apply": "true"})
        assert r["applied"] is True and app.applied

    def test_an_unreadable_apply_is_refused_not_guessed(self):
        app = self._app()
        r = self._h("annotate-from-presets")(app, {"apply": "maybe"})
        assert isinstance(r, tuple) and r[1] == 400 and app.applied == []

    def test_transfer_annotations_reads_apply_the_same_way(self):
        # The source lookup 404s AFTER the flag is read, so a refusal here can
        # only come from the flag.
        r = self._h("transfer-annotations")(
            self._app(), {"source_id": "nothing-by-this-name", "apply": "maybe"})
        assert isinstance(r, tuple) and r[1] == 400 and "apply" in r[0]["error"]

    def test_a_scan_target_sent_as_text_linear_is_linear(self):
        seq = _rand_seq(200, 7)
        import splicecraft_agent as _agent
        got = _agent._agent_resolve_scan_target(None, {"sequence": seq,
                                                       "circular": "false"})
        assert got[1] is False
        got = _agent._agent_resolve_scan_target(None, {"sequence": seq})
        assert got[1] is True                    # the documented default

    def test_a_read_reference_sent_as_text_linear_is_linear(self):
        seq = _rand_seq(200, 8)
        import splicecraft_agent as _agent
        _ref, circ, err = _agent._agent_resolve_read_reference(
            {"reference": seq, "circular": "false"})
        assert err is None and circ is False

    def test_all_projects_false_as_text_stays_in_the_active_project(self):
        r = self._h("search-experiments")(None, {"query": "gibson",
                                                 "all_projects": "false"})
        assert not isinstance(r, tuple), r
        assert r["scope"] != "all"


# ── S5: the Insert button's own `@id` reference opens its plasmid ────────────

class TestNotebookPlasmidRefsResolve:
    """The notebook's Insert button writes `@<entry id>`. The resolver searched
    display names fuzzily, so the id `pUC19_v2` never matched the entry it
    names (displayed `pUC19 v2`) and the click reported "No plasmid"."""

    @staticmethod
    def _ent(eid, name, seed):
        e = _d_ent(eid, seed)
        e["name"] = name
        return e

    @pytest.fixture
    def _lib(self):
        sc._save_collections([
            {"name": "Active", "plasmids": [self._ent("act-1", "pShared", 1)]},
            {"name": "Backbones", "plasmids": [
                self._ent("pUC19_v2", "pUC19 v2", 2),
                self._ent("bb-7", "pDup", 3)]},
            {"name": "Archive", "plasmids": [
                self._ent("arc-1", "pShared", 4),
                self._ent("arc-2", "pDup", 5)]},
        ])
        sc._set_active_collection_name("Active")
        sc._save_library([self._ent("act-1", "pShared", 1)])

    def test_the_buttons_id_spelling_opens_its_plasmid(self, _lib):
        chosen, _ = sc._resolve_plasmid_ref("pUC19_v2")
        assert chosen == {"collection": "Backbones", "id": "pUC19_v2"}

    def test_a_hand_typed_display_name_resolves_whatever_its_case(self, _lib):
        # A hand-typed ref is the DISPLAY name, and nobody types its case
        # exactly; "puc19 V2" is the entry the map shows as "pUC19 v2".
        chosen, _ = sc._resolve_plasmid_ref("puc19 V2")
        assert chosen == {"collection": "Backbones", "id": "pUC19_v2"}

    def test_a_same_name_in_the_active_collection_wins(self, _lib):
        chosen, _ = sc._resolve_plasmid_ref("pShared")
        assert chosen == {"collection": "Active", "id": "act-1"}

    def test_a_name_shared_elsewhere_is_ambiguous_not_guessed(self, _lib):
        chosen, matches = sc._resolve_plasmid_ref("pDup")
        assert chosen is None
        assert {(m["collection"], m["id"]) for m in matches} == {
            ("Backbones", "bb-7"), ("Archive", "arc-2")}

    def test_the_screen_loads_what_the_resolver_chose(self, _lib,
                                                      monkeypatch):
        loaded, notes = [], []

        class _App:
            _unsaved = False

            def notify(self, msg, **k):
                notes.append(msg)
        app = _App()
        monkeypatch.setattr(sc.ExperimentsScreen, "app",
                            property(lambda self: app), raising=False)
        scr = sc.ExperimentsScreen.__new__(sc.ExperimentsScreen)
        scr.dismiss = lambda *a, **k: None
        scr._switch_and_load_plasmid = lambda c, i: loaded.append((c, i))
        scr._open_plasmid_ref("pUC19_v2")
        assert loaded == [("Backbones", "pUC19_v2")] and not notes


# ── S10f: a notebook tag is read whole, or not at all ───────────────────────

class TestNotebookTagsReadWhole:
    """The `;` guard meant for HTML entities (`&amp;`) also dropped `@pUC19;`
    in prose; a token past the 64-character cap resolved as its first 64
    characters; and the Insert button wrote `@35S-GFP`, which no reader ever
    recognised because a tag starts with a letter."""

    def test_a_semicolon_after_a_plasmid_or_action_is_punctuation(self):
        assert sc._extract_plasmid_refs(
            "Transformed @pUC19; plated on LB.") == ["pUC19"]
        assert sc._extract_action_refs("!digest; then !ligate") == [
            "digest", "ligate"]

    def test_an_html_entity_is_still_not_a_gel(self):
        assert sc._extract_gel_refs("A &amp; B, then &gel-1; done") == []
        assert sc._extract_gel_refs("then &gel-1, done") == ["gel-1"]

    @pytest.mark.parametrize("extract,sigil", [
        (sc._extract_plasmid_refs, "@"), (sc._extract_action_refs, "!"),
        (sc._extract_gel_refs, "&")])
    def test_a_token_past_the_cap_is_not_cut_to_a_prefix(self, extract,
                                                        sigil):
        ok, long = "p" + "a" * 63, "p" + "a" * 64
        assert extract(f"see {sigil}{ok} here") == [ok]
        assert extract(f"see {sigil}{long} here") == []
        assert extract(f"see {sigil}{ok}.x here") == []

    @pytest.fixture
    def _lib(self):
        ents = []
        for eid, name, seed in (("35S-GFP", "35S-GFP", 1), ("1234", "pABC", 2),
                                ("pUC19_v2", "pUC19 v2", 3)):
            e = _d_ent(eid, seed)
            e["name"] = name
            ents.append(e)
        sc._save_collections([{"name": "Lab", "plasmids": ents}])
        sc._set_active_collection_name("Lab")
        sc._save_library(ents)

    def test_the_insert_button_writes_a_tag_that_reads_back(self, _lib):
        assert sc._plasmid_ref_token("Lab", "pUC19_v2") == "pUC19_v2"
        # an id that cannot be a tag falls back to a name that can
        assert sc._plasmid_ref_token("Lab", "1234") == "pABC"
        for tok in ("pUC19_v2", "pABC"):
            chosen, _ = sc._resolve_plasmid_ref(
                sc._extract_plasmid_refs(f"see @{tok} here")[0])
            assert chosen["collection"] == "Lab"

    def test_no_spelling_that_can_be_a_tag_writes_nothing(self, _lib):
        assert sc._plasmid_ref_token("Lab", "35S-GFP") is None


# ── S10g: a bare-string `tags` is no tags, never its letters ────────────────

class TestBareStringTagsAreNotLetters:
    """Loading does not normalise `tags`, so a hand-edited `"tags": "gibson"`
    reached every reader that iterated it: exports listed `g, i, b, s, o, n`,
    a duplicate got six one-letter tags, the agent returned the letters, and
    the editor showed them — for the next save to write back. Every reader,
    and the save normaliser, now reads a bare string the way the tags field
    takes it: comma-separated tags (reading it as NO tags, the first fix,
    blanked the field and the next save deleted the tag)."""

    ENTRY = {"id": "exp-aaaa0001", "title": "Prep", "body_md": "done",
             "tags": "gibson", "created_at": "2026-09-01T10:00:00",
             "updated_at": "2026-09-01T10:00:00"}

    def test_exports_do_not_list_letters(self):
        md = sc._experiment_markdown_document(dict(self.ENTRY))
        html = sc._experiment_html_document([dict(self.ENTRY)])
        assert "g, i, b" not in md and "g, i, b" not in html
        assert "gibson" in md and "gibson" in html

    def test_a_duplicate_does_not_get_letter_tags(self):
        d = sc._experiment_duplicate(dict(self.ENTRY), new_id="exp-aaaa0002")
        assert d["tags"] == ["gibson"]

    def test_the_agent_does_not_return_letters(self):
        import json
        sc._state._EXPERIMENTS_FILE.write_text(json.dumps(
            {"_schema_version": 1, "entries": [dict(self.ENTRY)]}))
        sc._state._experiments_cache = None
        e = dict(self.ENTRY, image_paths="img.png")
        sc._state._EXPERIMENTS_FILE.write_text(json.dumps(
            {"_schema_version": 1, "entries": [e]}))
        sc._state._experiments_cache = None
        r = sc._state._AGENT_HANDLERS["list-experiments"][0](None, {})
        [got] = r["experiments"]
        assert got["tags"] == ["gibson"] and got["image_paths"] == []


# ── S10h: a date that never happened renders as no date ─────────────────────

@pytest.mark.parametrize("stamp,shown", [
    ("2026-02-30", ""), ("2026-04-31", ""), ("2026-02-29", ""),
    ("0000-01-01", ""), ("2024-02-29", "FEB 29 2024"),
    ("2026.12.31", "DEC 31 2026"), ("2026-06-09T14:30", "JUN 9 2026 14:30"),
])
def test_only_real_calendar_dates_render(stamp, shown):
    assert sc._history_human_dt(stamp) == shown


# ── P10: a mutagenesis dead end names the ends only when an end is in the way ─

class TestMutagenesisDeadEndsNameTheRealCause:
    """The window search tries three extensions shorter than a primer by
    construction; each counted as a "too close to a sequence end" rejection,
    so a Tm / GC dead end in the middle of a gene blamed the ends too."""

    DNA = "ATG" + "AAT" * 99            # AT-only: every window fails Tm or GC

    def test_a_dead_end_mid_sequence_does_not_blame_the_ends(self):
        with pytest.raises(RuntimeError) as ei:
            sc._mut_design_inner(self.DNA, 40, "K", "N")
        assert "sequence end" not in str(ei.value)

    def test_a_mutation_at_an_end_still_says_so(self):
        with pytest.raises(RuntimeError) as ei:
            sc._mut_design_inner(self.DNA, 2, "K", "N")
        assert "sequence end" in str(ei.value)


# ── S10i: notebook entries order by time, not by wall-clock text ────────────

class TestEntriesOrderByTime:
    """Stamps carry the local UTC offset, so across a DST fall-back the later
    edit (`01:10-05:00`) sorted as OLDER than `01:30-04:00`, and entries made
    while travelling interleaved by time zone."""

    EARLIER = "2026-11-01T01:30:00-04:00"     # EDT: 05:30 UTC
    LATER = "2026-11-01T01:10:00-05:00"       # EST: 06:10 UTC

    def _entries(self):
        return [{"id": "exp-00000001", "title": "a", "body_md": "gibson @pUC19",
                 "tags": [], "updated_at": self.EARLIER},
                {"id": "exp-00000002", "title": "b", "body_md": "gibson @pUC19",
                 "tags": [], "updated_at": self.LATER},
                {"id": "exp-00000003", "title": "c", "body_md": "gibson @pUC19",
                 "tags": [], "updated_at": "not a date"}]

    def test_the_later_instant_is_later(self):
        assert sc._iso_instant(self.LATER) > sc._iso_instant(self.EARLIER)
        assert sc._iso_instant("2026-11-01T06:10:00Z") == \
            sc._iso_instant(self.LATER)
        assert sc._iso_instant("garbage") == float("-inf")

    def test_search_breaks_a_tie_by_the_newer_edit(self):
        hits = sc._experiment_search(self._entries(), "gibson")
        assert [h["entry"]["id"] for h in hits] == [
            "exp-00000002", "exp-00000001", "exp-00000003"]

    def test_backlinks_list_the_newer_edit_first(self):
        got = sc._experiments_referencing(self._entries(), ["pUC19"],
                                          kind="plasmid")
        assert [e["id"] for e in got] == [
            "exp-00000002", "exp-00000001", "exp-00000003"]


# ── S10k: a run that failed or was stopped is not "complete" ────────────────

class TestAnEndedRunIsNotASucceededRun:
    """A run that ENDED `failed` or `stopped` without a fault reading came back
    `crashed: False`; AUTOLAB printed "run failed — complete" in green over a
    full progress bar, and the agent envelope said `ok: true`."""

    FAILED = {"ran": True, "run_id": "R1", "run_status": "failed",
              "crashed": False, "run_errors": [{"detail": "no tip attached"}]}

    @pytest.mark.parametrize("res,ok", [
        ({"ran": False, "reason": "confirm-required"}, True),
        ({"ran": False, "reason": "door-open", "detail": "close it"}, False),
        ({"ran": True, "run_status": "running"}, True),
        ({"ran": True, "run_status": "succeeded", "crashed": False}, True),
        ({"ran": True, "run_status": "stopped", "crashed": False}, False),
        ({"ran": True, "run_status": "failed", "crashed": True}, False),
    ])
    def test_the_outcome_says_whether_the_run_did_its_job(self, res, ok):
        import splicecraft_opentrons as ot
        got, why = ot._ot2_run_outcome(res)
        assert got is ok and bool(why) is (not ok)

    def test_the_agent_reports_a_failed_run_as_not_ok(self, monkeypatch):
        import splicecraft_opentrons as ot
        from tests.test_opentrons import _good_plan
        monkeypatch.setattr(ot, "_ot2_run_protocol",
                            lambda *a, **k: dict(self.FAILED))
        r = sc._state._AGENT_HANDLERS["ot2-run"][0](
            None, {"host": "10.1.2.3", "confirm": True, **_good_plan()})
        # A non-2xx — a status-checking caller (the CLI exits 0 on any 2xx)
        # read the failed run as done — and a 4xx, never a 5xx: a 5xx is not
        # cached for an idempotent retry, so the retry would run it again.
        assert isinstance(r, tuple) and r[1] == 422
        assert r[0]["ok"] is False and "no tip attached" in r[0]["error"]

    def test_a_failed_position_check_is_not_ok_either(self, monkeypatch):
        import splicecraft_opentrons as ot
        from tests.test_opentrons import _good_plan
        monkeypatch.setattr(ot, "_ot2_run_protocol",
                            lambda *a, **k: dict(self.FAILED))
        r = sc._state._AGENT_HANDLERS["ot2-position-check"][0](
            None, {"host": "10.1.2.3", "confirm": True, **_good_plan()})
        assert isinstance(r, tuple) and r[1] == 422, r
        assert r[0]["ok"] is False and "no tip attached" in r[0]["error"]

    def test_a_dry_run_is_still_ok(self, monkeypatch):
        import splicecraft_opentrons as ot
        from tests.test_opentrons import _good_plan
        monkeypatch.setattr(ot, "_ot2_run_protocol", lambda *a, **k: {
            "ran": False, "reason": "confirm-required"})
        r = sc._state._AGENT_HANDLERS["ot2-position-check"][0](
            None, {"host": "10.1.2.3", **_good_plan()})
        assert not isinstance(r, tuple) and r["ok"] is True

    def _screen(self, monkeypatch):
        lines, notes, filled = [], [], []

        class _App:
            def notify(self, msg, **k):
                notes.append(k.get("severity"))
        app = _App()
        monkeypatch.setattr(sc.AutolabScreen, "app",
                            property(lambda self: app), raising=False)
        scr = sc.AutolabScreen.__new__(sc.AutolabScreen)
        scr._mon = lines.append
        scr._set_crash_banner = lambda faults: None
        scr._fill_progress_bar = lambda: filled.append(1)
        scr._finish_run_progress = lambda label: None
        return scr, lines, notes, filled

    def test_the_tab_does_not_call_a_failed_run_complete(self, monkeypatch):
        scr, lines, notes, filled = self._screen(monkeypatch)
        scr._run_done(dict(self.FAILED))
        assert not filled and notes == ["error"]
        assert not any("complete" in ln and "NOT" not in ln for ln in lines)
        assert any("no tip attached" in ln for ln in lines)

    def test_a_succeeded_run_still_completes(self, monkeypatch):
        scr, lines, notes, filled = self._screen(monkeypatch)
        scr._run_done({"ran": True, "run_status": "succeeded",
                       "crashed": False})
        assert filled and not notes
        assert any("complete" in ln for ln in lines)


# ── S10l: a pipette's generation is part of what the protocol asked for ─────

class TestPipetteGenerationsAreChecked:
    """The pre-run gate compared pipette names with the generation stripped: a
    GEN2 protocol passed on a GEN1 pipette (and failed after the run started),
    while a GEN1 protocol was refused on the GEN2 pipette Opentrons runs it on."""

    @staticmethod
    def _check(want, model):
        import splicecraft_opentrons as ot
        return ot._ot2_pipette_mismatch(
            [{"mount": "left", "pipetteName": want}],
            [{"mount": "left", "model": model}])

    @pytest.mark.parametrize("want,model", [
        ("p300_single_gen2", "p300_single_v1.5"),
        ("p20_single_gen2", "p20_multi_v2.1"),
        ("p300_single", "p1000_single_v2.0"),
        # opentrons-shared-data 9.1.2: a P300 GEN2's `backCompatNames` is
        # `p300_single` only — it does not stand in for a P50 (round-2
        # hardening, 2026-09-25; this pair was wrongly listed as accepted).
        ("p50_single", "p300_single_v2.1"),
    ])
    def test_a_pipette_that_cannot_run_it_is_refused(self, want, model):
        assert self._check(want, model)

    @pytest.mark.parametrize("want,model", [
        ("p10_single", "p20_single_v2.2"),
        ("p300_multi", "p300_multi_v2.1"),
        ("p300_single_gen2", "p300_single_v2.1"),
        ("p300_single_gen2", "p300_single_gen2"),
        ("p300_single", "p300_single_v1.5"),
    ])
    def test_a_pipette_that_runs_it_is_accepted(self, want, model):
        assert self._check(want, model) == []


# ── S7: the GUI will not plate a source onto itself ──────────────────────────

class TestAutolabRefusesTheSourceAsDestination:
    """Destination wells fill in order, independently of the source wells, so
    with one plate for both a step overwrites a sample a later step still has
    to read. The agent's `ot2-normalize` refused it; the AUTOLAB tab built it.
    Only the widget reads are stubbed — the builders run for real."""

    def _screen(self, monkeypatch, dst_slot):
        notes, built = [], []

        class _App:
            def notify(self, msg, **k):
                notes.append((msg, k.get("severity")))
        app = _App()
        monkeypatch.setattr(sc.AutolabScreen, "app",
                            property(lambda self: app), raising=False)
        scr = sc.AutolabScreen.__new__(sc.AutolabScreen)
        scr._deck = {
            1: {"id": "src", "labware": "plate_24",
                "map": {"A1": {"id": "p1", "name": "p1"},
                        "A2": {"id": "p2", "name": "p2"}}},
            2: {"id": "dst", "labware": "plate_24"},
        }
        scr._bound_slot = 1
        vals = {"#autolab-lib-dst": str(dst_slot)}
        nums = {"#autolab-lib-vol": 20, "#autolab-lib-target": 100,
                "#autolab-lib-final": 20}
        scr._val = lambda sel, default="": vals.get(sel, default)
        scr._read_num = lambda sel: nums.get(sel)
        scr._parse_concentrations = lambda: {"A1": 50.0, "A2": 80.0}
        scr._mon = lambda markup: None
        scr._append_steps = lambda steps, what: built.append(steps)
        return scr, notes, built

    @pytest.mark.parametrize("builder", ["_cherry_pick", "_normalize_build"])
    def test_the_source_plate_is_refused_as_destination(self, monkeypatch,
                                                        builder):
        scr, notes, built = self._screen(monkeypatch, dst_slot=1)
        getattr(scr, builder)()
        assert built == []
        assert any("source plate" in m and sev == "error" for m, sev in notes)

    @pytest.mark.parametrize("builder", ["_cherry_pick", "_normalize_build"])
    def test_a_second_plate_still_builds(self, monkeypatch, builder):
        scr, notes, built = self._screen(monkeypatch, dst_slot=2)
        getattr(scr, builder)()
        assert built and built[0], notes
        moves = [s for s in built[0] if s.get("type") == "transfer"]
        assert moves and all(str(s["to"]).startswith("dst:") for s in moves)


# ── S4 / SEC10: an exported notebook shows each attachment once, working ─────

class TestNotebookExportShowsEachAttachmentOnce:
    """The Attach button inserts `![name](name)` into the body AND records the
    file; the exporters then listed every recorded file again under
    Attachments. In HTML the body copy was also broken — its bare stored name
    resolves nowhere once the page leaves the app — so the user got a broken
    image followed by a working one."""

    DATA = "data:image/png;base64,AAAA"

    @staticmethod
    def _e(body, paths):
        return {"id": "exp-1", "title": "Gel day", "body_md": body,
                "tags": [], "image_paths": paths}

    def _srcs(self):
        return {"gel.png": self.DATA, "trace.ab1": "file:///lab/trace.ab1"}

    def test_markdown_shows_a_body_image_once_and_resolved(self):
        import splicecraft_experiments as E
        e = self._e("Ran @pUC19.\n\n![gel.png](gel.png)\n", ["gel.png",
                                                            "trace.ab1"])
        doc = E._experiment_markdown_document(e, image_srcs=self._srcs())
        assert doc.count("gel.png") == 1              # the alt text, once
        assert f"![gel.png]({self.DATA})" in doc      # pointed at the file
        assert "Ran @pUC19." in doc                   # the rest verbatim
        # the file the body never mentions is still listed
        att = doc.split("## Attachments", 1)[1]
        # a local PATH, not a `file:` URI — CommonMark viewers refuse that
        # scheme (hardening 2026-09-25); the HTML export keeps the URI
        assert "[trace.ab1](/lab/trace.ab1)" in att

    def test_html_shows_a_body_image_once_and_resolved(self):
        import splicecraft_experiments as E
        e = self._e("![gel.png](gel.png)", ["gel.png"])
        html = E._experiment_html_document([e], image_srcs=self._srcs())
        assert html.count("<img") == 1
        assert f'<img src="{self.DATA}"' in html
        assert "<h2>Attachments</h2>" not in html

    def test_a_mention_inside_code_is_not_showing_it(self):
        import splicecraft_experiments as E
        e = self._e("```\n![gel.png](gel.png)\n```\nand `![gel.png](gel.png)`",
                    ["gel.png"])
        doc = E._experiment_markdown_document(e, image_srcs=self._srcs())
        body, att = doc.split("## Attachments", 1)
        assert "![gel.png](gel.png)" in body          # code left as written
        assert self.DATA not in body and self.DATA in att

    def test_a_local_file_uri_survives_the_html_page(self):
        # A Windows attachment left on disk: handed over raw, `C:\...` read
        # as a URL scheme and the attachment vanished from the page.
        import splicecraft_experiments as E
        uri = "file:///C:/Users/lab/notebook/gel.png"   # Path.as_uri() on Windows
        e = self._e("", ["gel.png"])
        html = E._experiment_html_document([e], image_srcs={"gel.png": uri})
        assert f'<img src="{uri}"' in html

    def test_a_file_link_typed_into_a_body_is_still_refused(self):
        import splicecraft_experiments as E
        html = E._markdown_subset_to_html("[x](file:///etc/passwd)")
        assert "file:" not in html

    def test_unclosed_link_runs_render_in_linear_time(self):
        import time
        import splicecraft_experiments as E
        body = "![a](" * 200_000                      # the 1 MB body cap
        t = time.perf_counter()
        E._experiment_html_document([self._e(body, [])])
        E._experiment_markdown_document(self._e(body, ["gel.png"]),
                                        image_srcs=self._srcs())
        assert time.perf_counter() - t < 2.0          # was ~45 s at 160 KB

    def test_a_numbered_list_keeps_its_start(self):
        import splicecraft_experiments as E
        html = E._markdown_subset_to_html("1. one\n\nsome prose\n\n2. two")
        assert '<ol start="2">' in html and "<li>two</li>" in html
        assert "<ol>" in html                          # the first starts at 1


# ── SEC4: a name shows as itself inside the app's markup ─────────────────────

_TRICKY_NAMES = ["Backups\\", "pUC19 [v2]", "a\\[b]", "x [ y ] z",
                 "pUC19 [3' end]", "[@click=app.quit]q", "lab [/b] stock",
                 "\\\\", "[red]", "tail [", "] lead", "p[1]\\"]


class TestNamesSurviveMarkup:
    """Textual 8 cannot escape a backslash: it reads any `[` after one as
    literal. A name ENDING in a backslash therefore escaped the template's
    own closing tag ("Delete collection Backups\\[/bold]?"), and a bracketed
    word vanished from a toast ("pUC19 [v2]" → "pUC19 ")."""

    @pytest.mark.parametrize("name", _TRICKY_NAMES)
    def test_an_escaped_name_renders_as_itself_in_both_parsers(self, name):
        from rich.text import Text
        from textual.content import Content
        t = f"Delete [bold]{sc._markup_escape(name)}[/bold]?"
        for parsed in (Content.from_markup(t), Text.from_markup(t)):
            assert parsed.plain.replace("​", "") == f"Delete {name}?"
            spans = list(parsed.spans)
            assert len(spans) == 1, spans              # the template's bold
            assert spans[0].end - spans[0].start == len(
                parsed.plain) - len("Delete ?")        # …covering the name

    def test_a_toast_keeps_the_bracketed_word_and_the_colour(self, monkeypatch):
        from textual.app import App
        from textual.content import Content
        from textual.notifications import Notification
        from textual.widgets._toast import Toast
        sent = []
        monkeypatch.setattr(App, "notify",
                            lambda self, message, **k: sent.append((message, k)))
        monkeypatch.setattr(sc.PlasmidApp, "screen",
                            property(lambda self: object()), raising=False)
        app = sc.PlasmidApp()
        app.notify("[green]Imported 'pUC19 [v2]'[/green]")
        msg, kw = sent[-1]
        shown = Toast(Notification(msg, "", "information", 3,
                                   markup=kw.get("markup", True))).render()
        assert isinstance(shown, Content)
        assert shown.plain == "Imported 'pUC19 [v2]'"
        assert shown.spans                             # still green

    def test_a_dialog_line_keeps_the_bracketed_word_and_the_bold(self):
        from textual.widgets import Static
        from textual.widgets import _static
        shown = getattr(_static, "visualize")(
            Static(), "Delete collection [bold]pUC19 [v2][/bold]?")
        assert shown.plain == "Delete collection pUC19 [v2]?"
        assert shown.spans


# ── S1: two custom labware never share a load name on the robot ──────────────

class TestCustomLabwareLoadNamesStayDistinct:
    """The robot files an embedded definition under its load name. "Rack α"
    and "Rack β" both became `rack`, so two different racks were one labware
    type to it — the second driven with the first one's well geometry."""

    import re as _re
    _LOAD_NAME = _re.compile(r"^[a-z0-9._]+$")

    def test_letters_without_an_ascii_form_do_not_collide(self):
        import splicecraft_opentrons as ot
        a = ot._ot2_build_labware_def("Rack α", 4, 6, z_dim=80.0)["parameters"]["loadName"]
        b = ot._ot2_build_labware_def("Rack β", 4, 6, z_dim=80.0)["parameters"]["loadName"]
        assert a != b
        assert self._LOAD_NAME.match(a) and self._LOAD_NAME.match(b)

    def test_names_are_stable_and_readable(self):
        import splicecraft_opentrons as ot
        n = ot._ot2_custom_load_name
        assert n("Rack α") == n("Rack α")                # deterministic
        assert n("Crème rack") == "creme_rack"           # accents fold
        assert n("Deep well 96") == "deep_well_96"       # ASCII unchanged
        assert n("株式") != n("会社")                     # no ASCII at all
        assert self._LOAD_NAME.match(n("株式"))

    @staticmethod
    def _plan(def_a, def_b):
        return {"pipette": "p300_single", "mount": "left",
                "tips": {"labware": "tiprack_300", "slot": 8},
                "labware": {"a": {"labware": "custom", "slot": 1,
                                  "definition": def_a},
                            "b": {"labware": "custom", "slot": 2,
                                  "definition": def_b}},
                "transfers": [{"from": "a:A1", "to": "b:A1", "volume": 50}]}

    def test_two_saved_definitions_sharing_a_name_are_kept_apart(self):
        # Definitions saved BEFORE the fix can still share a load name.
        import splicecraft_opentrons as ot
        a = ot._ot2_build_labware_def("rack", 4, 6, z_dim=80.0)
        b = ot._ot2_build_labware_def("rack", 3, 4, spacing=20.0, z_dim=80.0)
        assert a["parameters"]["loadName"] == b["parameters"]["loadName"]
        plan = self._plan(a, b)
        report = ot._ot2_validate_plan(plan)
        assert not report["errors"], report
        assert any("rack_2" in w for w in report["warnings"])
        proto = ot._ot2_compile_protocol(plan)
        assert "'loadName': 'rack'" in proto and "'loadName': 'rack_2'" in proto
        # the saved definition itself is untouched
        assert b["parameters"]["loadName"] == "rack"

    def test_the_same_definition_twice_is_one_type(self):
        import splicecraft_opentrons as ot
        a = ot._ot2_build_labware_def("rack", 4, 6, z_dim=80.0)
        plan = self._plan(a, a)
        report = ot._ot2_validate_plan(plan)
        assert not any("load name" in w for w in report["warnings"])
        assert "rack_2" not in ot._ot2_compile_protocol(plan)


# ── SEC10a / SEC10d: exported maps are well-formed and never overwrite ───────

class TestExportedMapsAreWellFormedAndDistinct:
    @staticmethod
    def _rec(label):
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord
        r = SeqRecord(Seq("ACGT" * 200), id="p", name="p",
                      annotations={"molecule_type": "DNA",
                                   "topology": "circular"})
        r.features.append(SeqFeature(FeatureLocation(10, 200, strand=1),
                                     type="CDS", qualifiers={"label": [label]}))
        return r

    @pytest.mark.parametrize("bad", ["￾", "￿", "\x01"])
    def test_an_xml_illegal_character_leaves_the_svg_readable(self, bad):
        # XML 1.0 forbids these outright; one in a label or the title made
        # the whole SVG unreadable to every viewer.
        import xml.etree.ElementTree as ET
        import splicecraft_mapimage as mi
        feats, total = mi._map_feats_from_record(self._rec(f"GFP{bad}x"))
        feats[0]["label"] = f"GFP{bad}x"         # as if unsanitised upstream
        svg = mi.render_plasmid_map_svg(feats, total, title=f"p{bad}1",
                                        size=600)
        ET.fromstring(svg)
        assert "GFP" in svg

    def test_names_never_carry_a_noncharacter(self):
        assert sc._sanitize_label("pUC￾19￿") == "pUC19"

    def test_bulk_map_images_do_not_overwrite_a_decomposed_twin(self, tmp_path):
        # macOS stores names decomposed: "Plasmidé" and "Plasmidé" are
        # ONE file there, so the second image silently replaced the first.
        import unicodedata
        gb = sc._record_to_gb_text(self._rec("GFP"))
        entries = [{"name": "Plasmidé", "gb_text": gb},
                   {"name": "Plasmidé", "gb_text": gb}]
        r = sc._bulk_export_map_images(entries, tmp_path, fmt="svg")
        assert not r["failures"], r["failures"]
        keys = {unicodedata.normalize("NFC", w["path"]).casefold()
                for w in r["written"]}
        assert len(keys) == 2


# ── S10c: a well is written the way the robot names it ────────────────────

class TestWellsAreWrittenAsTheRobotNamesThem:
    """The validator read `a1`, `A01` and a fullwidth digit as A1, but the
    compiler wrote the well as typed and the robot looks wells up by exact
    name: a plan SpliceCraft called valid was refused by the robot's own
    analysis."""

    @staticmethod
    def _plan(src_well, dst_well="A1", definition=None):
        dst = ({"definition": definition, "slot": 2} if definition
               else {"labware": "plate_24", "slot": 2})
        return {"pipette": "p300_single", "mount": "left",
                "tips": {"labware": "tiprack_300", "slot": 8},
                "labware": {"src": {"labware": "eppi_24", "slot": 1},
                            "dst": dst},
                "transfers": [{"from": f"src:{src_well}",
                               "to": f"dst:{dst_well}", "volume": 50}]}

    @staticmethod
    def _wells_in(proto):
        import re
        return re.findall(r'(?:src|dst)\["([^"]*)"\]', proto)

    @pytest.mark.parametrize("spelled", ["a1", "A01", "A\uff11", " a1 "])
    def test_a_loose_spelling_compiles_to_the_robots_name(self, spelled):
        import splicecraft_opentrons as ot
        plan = self._plan(spelled, dst_well=spelled)
        assert not ot._ot2_validate_plan(plan)["errors"]
        assert self._wells_in(ot._ot2_compile_protocol(plan)) == ["A1", "A1"]

    def test_a_custom_definition_keeps_its_own_key(self):
        import splicecraft_opentrons as ot
        d = ot._ot2_build_labware_def("Tube block", rows=2, cols=3,
                                      volume=1500.0, depth=40.0, z_dim=80.0)
        plan = self._plan("A1", dst_well="b02", definition=d)
        assert not ot._ot2_validate_plan(plan)["errors"]
        assert self._wells_in(ot._ot2_compile_protocol(plan))[-1] == "B2"

    def test_a_letter_that_is_not_a_row_is_still_refused(self):
        import splicecraft_opentrons as ot
        errs = ot._ot2_validate_plan(self._plan("\uff21\uff11"))["errors"]
        assert errs and "\uff21\uff11" in errs[0]

    def test_the_load_budget_counts_one_well_once(self):
        import splicecraft_opentrons as ot
        plan = self._plan("a1")
        plan["transfers"].append({"from": "src:A01", "to": "dst:A2",
                                  "volume": 30})
        assert ot._ot2_plan_summary(plan)["source_volumes"] == {"src:A1": 80.0}


# ── S10d: a distribute budgets the liquid the robot actually draws ────────

class TestDistributeBudgetsWhatTheRobotDraws:
    """`distribute` aspirates the pipette's minimum volume on top of every
    trip and blows it out to the trash; the load budget left it out, so the
    source ran dry before the last well. The summary's total also counted a
    distribute's per-well volume once."""

    @staticmethod
    def _plan(pipette, tips, vol, n, extra=()):
        from splicecraft_opentrons import _ot2_wells
        dsts = [f"dst:{w}" for w in _ot2_wells(8, 12)[:n]]
        return {"pipette": pipette, "mount": "left",
                "tips": {"labware": tips, "slot": 8},
                "labware": {"src": {"labware": "eppi_24", "slot": 1},
                            "dst": {"labware": "nest_96_wellplate_2ml_deep",
                                    "slot": 2}},
                "steps": [{"type": "distribute", "from": "src:A1",
                           "to": dsts, "volume": vol}, *extra]}

    def test_every_trip_draws_the_disposal_volume(self):
        import splicecraft_opentrons as ot
        # p1000: 100 uL disposal, 18 x 50 uL per trip -> 24 wells take 2 trips
        s = ot._ot2_plan_summary(self._plan("p1000_single_gen2",
                                            "tiprack_1000", 50, 24))
        assert s["source_volumes"] == {"src:A1": 1200.0 + 2 * 100.0}

    def test_smaller_tips_mean_more_trips(self):
        import splicecraft_opentrons as ot
        # p300 on 200 uL tips: 20 uL disposal leaves room for 9 x 20 uL
        s = ot._ot2_plan_summary(self._plan(
            "p300_single_gen2", "opentrons_96_filtertiprack_200ul", 20, 12))
        assert s["source_volumes"] == {"src:A1": 240.0 + 2 * 20.0}

    def test_a_dispense_bigger_than_a_trip_is_split(self):
        import splicecraft_opentrons as ot
        # 290 uL does not fit beside a 20 uL disposal in 300: two trips a well
        s = ot._ot2_plan_summary(self._plan("p300_single_gen2",
                                            "tiprack_300", 290, 3))
        assert s["source_volumes"] == {"src:A1": 870.0 + 6 * 20.0}

    def test_the_total_counts_every_well_served_and_no_mix(self):
        import splicecraft_opentrons as ot
        plan = self._plan("p300_single_gen2", "tiprack_300", 20, 12, extra=[
            {"type": "mix", "at": "dst:A1", "volume": 100, "repetitions": 3},
            {"type": "consolidate", "from": ["dst:A1", "dst:A2", "dst:A3"],
             "to": "src:B1", "volume": 10}])
        assert ot._ot2_plan_summary(plan)["total_volume_ul"] == 240 + 30


# ── FM11d: an oversized .dna XML block is refused before anything parses it ──

class TestOversizedDnaBlocksAreRefusedUnread:
    """Biopython's `.dna` reader builds a minidom tree of every XML packet
    before SpliceCraft's own packet readers see any of them, so their per-block
    cap guarded nothing: a 200 MB features block was expanded in memory first.
    The pre-walk reads only the 5-byte packet headers."""

    @staticmethod
    def _dna(tmp_path, features_xml: bytes):
        import splicecraft_fileio as fio
        cookie = sc._COMMERCIALSAAS_COOKIE_MAGIC + bytes([0, 1, 0, 15, 0, 19])
        p = tmp_path / "big.dna"
        p.write_bytes(
            fio._build_commercialsaas_packet(0x09, cookie)
            + fio._build_commercialsaas_packet(0x00, b"\x03" + b"ACGT" * 50)
            + fio._build_commercialsaas_packet(0x0A, features_xml))
        return p

    @staticmethod
    def _features(n):
        body = b"".join(
            b'<Feature name="f%d" type="misc_feature" directionality="1">'
            b'<Segment range="%d-%d" type="standard"/></Feature>'
            % (i, 1 + i % 100, 50 + i % 100) for i in range(n))
        return (b'<?xml version="1.0"?><Features nextValidID="%d">' % n
                + body + b"</Features>")

    def test_an_oversized_block_is_refused_before_it_is_parsed(
            self, tmp_path, monkeypatch):
        import importlib
        # Biopython's `.dna` reader, named through the hex-encoded format
        # string ([PIT-18]).
        sgio = importlib.import_module(
            "Bio.SeqIO." + sc._BIOPYTHON_DNA_FMT.capitalize()
            .replace("gene", "Gene") + "IO")
        import splicecraft_fileio as fio
        # Built BEFORE the cap is lowered: the writer refuses a block over the
        # cap too (harden 2026-09-24), so SpliceCraft never writes a `.dna`
        # it would refuse to read.
        p = self._dna(tmp_path, self._features(200))
        assert p.stat().st_size > 4096
        monkeypatch.setattr(fio, "_COMMERCIALSAAS_PACKET_MAX_XML", 4096)
        parsed = []

        def _minidom(*a, **k):
            parsed.append(1)
            raise AssertionError("minidom parsed the oversized block")
        monkeypatch.setattr(sgio, "parseString", _minidom)
        with pytest.raises(ValueError, match=r"features block is [0-9.]+ MB"):
            sc.load_genbank(str(p))
        assert parsed == []

    def test_a_block_under_the_cap_still_opens(self, tmp_path):
        rec = sc.load_genbank(str(self._dna(tmp_path, self._features(20))))
        assert len(rec.seq) == 200
        assert sum(f.type == "misc_feature" for f in rec.features) == 20

    def test_a_file_that_is_not_a_dna_keeps_the_parsers_own_error(
            self, tmp_path):
        # No cookie: whatever its bytes say, "block too large" is not the story.
        p = tmp_path / "not.dna"
        p.write_bytes(b"\x0a\xff\xff\xff\xff" + b"x" * 64)
        with pytest.raises(ValueError) as ei:
            sc.load_genbank(str(p))
        assert "block is" not in str(ei.value)


def test_a_non_ascii_bearer_token_is_a_401_not_a_dead_connection():
    import inspect
    src = inspect.getsource(sc._AgentRequestHandler._check_token)
    # `secrets.compare_digest` RAISES TypeError on a non-ASCII str, and that
    # escaped the handler leaving the request with no response at all.
    assert "isascii()" in src


def test_a_chunked_body_is_refused_not_treated_as_empty():
    import inspect
    src = inspect.getsource(sc._AgentRequestHandler._read_body)
    assert "chunked" in src and "__chunked_body__" in src


class TestAgentPortValidation:
    def test_the_range_is_checked_at_parse_time(self):
        import inspect
        src = inspect.getsource(sc.main)
        assert "_valid_port" in src, \
            "an out-of-range port still reaches bind() inside on_mount"

    def test_zero_is_excluded(self):
        # Port 0 makes the kernel pick a random free port and nothing tells the
        # user which — the API would run at an address they cannot find.
        import inspect
        src = inspect.getsource(sc.main)
        assert "1 <= value <= 65535" in src


def test_a_qualifier_name_must_be_legal_genbank():
    import splicecraft_agent as _agent
    ok, err = _agent._agent_sanitize_qualifiers({"gene": "ampR"})
    assert err is None and ok == {"gene": ["ampR"]}
    bad, err = _agent._agent_sanitize_qualifiers({"my qualifier": "x"})
    assert bad is None and "not a legal GenBank qualifier" in err


def test_an_empty_tag_is_not_a_search_filter():
    import splicecraft_agent as _agent
    got = _agent._h_search_experiments(None, {"tags": [""]})
    assert isinstance(got, tuple) and got[1] == 400
    assert "empty tag" in got[0]["error"]


def test_a_read_only_guest_does_not_delete_the_hosts_token(tmp_path):
    import inspect
    src = inspect.getsource(sc._stop_agent_api)
    assert "_AGENT_READ_ONLY" in src and "_AGENT_TOKEN_READONLY_NAME" in src, \
        "a guest's shutdown still unlinks the shared agent_token"


# ══════════════════════════════════════════════════════════════════════════
# OT-2
# ══════════════════════════════════════════════════════════════════════════

class TestOt2HostIsNotAFetchPrimitive:
    def test_a_normal_robot_address_works(self, monkeypatch):
        import socket
        import splicecraft_opentrons as _ot
        monkeypatch.setattr(socket, "getaddrinfo", lambda *a, **k: [
            (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("192.168.1.56", 0))])
        assert _ot._ot2_base_url("192.168.1.56") == "http://192.168.1.56:31950"
        assert _ot._ot2_base_url("opentrons.local") == "http://192.168.1.56:31950"

    def test_credentials_in_the_host_are_refused(self):
        import splicecraft_opentrons as _ot
        with pytest.raises(_ot.OT2Error, match="credentials"):
            _ot._ot2_base_url("user:pw@elsewhere.example")

    def test_the_metadata_service_is_refused_however_it_is_spelled(self):
        # 2852039166 / 0xa9fea9fe / 169.254.43518 / octal all resolve to
        # 169.254.169.254 but are not IP literals to `ipaddress`, so they used
        # to pass as "hostnames".
        import splicecraft_opentrons as _ot
        for host in ("169.254.169.254", "http://169.254.169.254:80",
                     "2852039166", "0xa9fea9fe", "169.254.43518",
                     "0251.0376.0251.0376", "100.100.100.200"):
            with pytest.raises(_ot.OT2Error):
                _ot._ot2_base_url(host)

    def test_a_usb_robot_is_reachable(self):
        # Refusing the whole link-local range (the first fix) refused every
        # USB-connected OT-2 — the addresses the module's own discovery reports.
        import splicecraft_opentrons as _ot
        assert _ot._ot2_base_url("169.254.33.7") == "http://169.254.33.7:31950"
        assert _ot._ot2_base_url("[fe80::1]") == "http://[fe80::1]:31950"
        assert _ot._ot2_base_url("100.101.102.103").endswith(":31950")  # Tailscale

    def test_a_host_is_not_a_url(self):
        import splicecraft_opentrons as _ot
        for host in ("evil.example/x?s=1", "a b.local", "host%2f.local",
                     "8.8.8.8"):
            with pytest.raises(_ot.OT2Error):
                _ot._ot2_base_url(host)

    def test_an_ipv6_literal_gets_the_default_port(self):
        import splicecraft_opentrons as _ot
        assert _ot._ot2_base_url("[2001:db8::5]") == \
            "http://[2001:db8::5]:31950"

    def test_the_netloc_split_is_ipv6_aware(self):
        import splicecraft_opentrons as _ot
        assert _ot._ot2_split_netloc("[fe80::1]") == ("fe80::1", False)
        assert _ot._ot2_split_netloc("[fe80::1]:31950") == ("fe80::1", True)
        assert _ot._ot2_split_netloc("host:8080") == ("host", True)
        assert _ot._ot2_split_netloc("host") == ("host", False)


class TestCustomLabwareGeometry:
    def test_a_height_is_never_guessed(self):
        # The first fix derived "depth + 3 mm": a real 24-tube rack (≈80 mm
        # tall, tube bottoms ≈42 mm up) came out with its wells at 3 mm — the
        # tip driven ~38 mm into the tubes — and a 14.78 mm PCR plate's well
        # bottoms rose ~2 mm, so a small aspirate drew air.
        import splicecraft_opentrons as _ot
        with pytest.raises(_ot.OT2Error, match="overall height"):
            _ot._ot2_build_labware_def("Tube rack", 4, 6, category="tubeRack",
                                       depth=40.0, volume=1500.0)

    def test_a_stated_height_puts_the_wells_where_it_says(self):
        # Opentrons' own 24-tube 1.5 mL rack: 79.85 mm tall, 37.9 mm tubes.
        import splicecraft_opentrons as _ot
        d = _ot._ot2_build_labware_def("Tube rack", 4, 6, category="tubeRack",
                                       depth=37.9, volume=1500.0, z_dim=79.85)
        assert d["dimensions"]["zDimension"] == 79.85
        assert {w["z"] for w in d["wells"].values()} == {41.95}

    @pytest.mark.parametrize("junk", ["abc", float("nan"), float("inf"), 0.5])
    def test_a_junk_height_is_named_as_junk(self, junk):
        # It used to be replaced by 15 mm and then reported as "labware height
        # 15 mm" — a number the user never typed.
        import splicecraft_opentrons as _ot
        with pytest.raises(_ot.OT2Error, match="millimetres"):
            _ot._ot2_build_labware_def("X", 2, 2, depth=10.0, z_dim=junk)

    def test_saved_labware_with_wells_below_the_deck_is_refused(self):
        # Every tube rack the old builder saved (40 mm wells in a "15 mm"
        # labware, z = -25) is still in users' libraries: stopped before the
        # robot is ever asked to run it.
        import splicecraft_opentrons as _ot
        good = _ot._ot2_build_labware_def("Rack", 4, 6, depth=37.9, z_dim=79.85)
        bad = {**good, "dimensions": {**good["dimensions"], "zDimension": 15.0},
               "wells": {k: {**w, "z": 15.0 - w["depth"]}
                         for k, w in good["wells"].items()}}
        def plan(definition):
            return {"pipette": "p300_single", "tips": {"labware": "tiprack_300",
                                                       "slot": 8},
                    "labware": {"r": {"slot": 1, "definition": definition},
                                "p": {"labware": "plate_96", "slot": 2}},
                    "transfers": [{"from": "r:A1", "to": "p:A1", "volume": 50}]}
        assert not _ot._ot2_validate_plan(plan(good))["errors"]
        errs = _ot._ot2_validate_plan(plan(bad))["errors"]
        assert any("BELOW the deck" in e for e in errs)

    def test_an_explicit_contradiction_is_refused(self):
        import splicecraft_opentrons as _ot
        with pytest.raises(_ot.OT2Error, match="below the deck"):
            _ot._ot2_build_labware_def("Bad", 4, 6, depth=40.0, z_dim=15.0)

    def test_the_load_name_is_ascii(self):
        import re as _re
        import splicecraft_opentrons as _ot
        d = _ot._ot2_build_labware_def("Plaque Épaisse 96", 8, 12, z_dim=15.7)
        assert _re.fullmatch(r"[a-z0-9._]+", d["parameters"]["loadName"])


class TestOt2PlanChecks:
    LW = {"src": {"labware": "biorad_96_wellplate_200ul_pcr", "slot": 1},
          "dst": {"labware": "biorad_96_wellplate_200ul_pcr", "slot": 2}}
    TX = [{"from": "src:A1", "to": "dst:A1", "volume": 50}]

    def _validate(self, pipette, rack="tiprack_300", transfers=None):
        import splicecraft_opentrons as _ot
        return _ot._ot2_validate_plan({
            "pipette": pipette, "labware": self.LW,
            "transfers": transfers or self.TX,
            "tips": [{"labware": rack, "slot": 8}]})

    def test_a_smaller_tip_is_a_supported_pairing_not_an_error(self):
        # A p300 on 200 µL filter tips is standard; the robot caps each
        # aspirate at the tip and splits bigger volumes itself. The first fix
        # refused it — even for this 50 µL plan.
        v = self._validate("p300_single",
                           rack="opentrons_96_filtertiprack_200ul")
        assert not v["errors"], v["errors"]
        assert any("limited to 200" in w for w in v["warnings"])
        big = self._validate("p300_single",
                             rack="opentrons_96_filtertiprack_200ul",
                             transfers=[{"from": "src:A1", "to": "dst:A1",
                                         "volume": 280}])
        assert not big["errors"], big["errors"]

    def test_a_matched_pair_passes(self):
        assert not self._validate("p300_single", "tiprack_300")["errors"]
        assert not self._validate("p20_single_gen2",
                                  "opentrons_96_tiprack_20ul")["errors"]

    def test_a_tip_the_nozzle_cannot_take_is_an_error(self):
        # Opentrons' own pipette definitions (`supportedTips`): a P20 takes 10 /
        # 20 µL tips, a P300 200 / 300, a P1000 1000.
        for pip, rack in (("p20_single_gen2", "tiprack_300"),
                          ("p300_single_gen2", "tiprack_20"),
                          ("p1000_single_gen2", "tiprack_300")):
            v = self._validate(pip, rack)
            assert any("cannot pick up" in e for e in v["errors"]), (pip, rack)

    def test_a_p50_on_its_standard_300_ul_tips_is_not_flagged(self):
        # The old "much larger than its range" guess flagged exactly this.
        v = self._validate("p50_single", "tiprack_300",
                           transfers=[{"from": "src:A1", "to": "dst:A1",
                                       "volume": 20}])
        assert not v["errors"]
        assert not [w for w in v["warnings"] if "seat" in w or "tips" in w]

    def test_a_multichannel_says_it_moves_a_whole_column(self):
        v = self._validate("p300_multi_gen2")
        assert any("channels" in w and "whole column" in w
                   for w in v["warnings"])

    def test_a_multichannel_below_row_a_is_an_error(self):
        v = self._validate("p300_multi_gen2",
                           transfers=[{"from": "src:B3", "to": "dst:B3",
                                       "volume": 50}])
        assert any("dst:B3 (row A only)" in e for e in v["errors"])

    def test_a_single_channel_gets_neither(self):
        v = self._validate("p300_single")
        assert not [w for w in v["warnings"] if "channels" in w]


def test_normalise_reports_which_wells_missed_the_target():
    import inspect
    import splicecraft_opentrons as _ot
    src = inspect.getsource(_ot._ot2_normalize_volumes)
    assert "on_target" in src, \
        "a clamped, under-target well still reports only `ok`"


class TestNormaliseOnOnePlate:
    """One plate for source and destination is refused only when a well would
    be both read and written — the first fix refused normalising column 1
    into column 7 of the same plate."""

    ITEMS = [{"name": f"s{r}", "well": f"{r}1", "concentration": 50.0}
             for r in "ABCDEFGH"]

    def _call(self, **extra):
        import splicecraft as sc
        body = {"items": self.ITEMS, "target_ng": 100, "final_volume": 20,
                "src": "plate", "dst": "plate", **extra}
        return sc._state._AGENT_HANDLERS["ot2-normalize"][0](None, body)

    def test_disjoint_destination_wells_build(self):
        r = self._call(dst_wells=[f"{r}7" for r in "ABCDEFGH"])
        assert not isinstance(r, tuple), r
        assert r["steps"] and all(s["to"].endswith("7") for s in r["steps"])

    def test_a_destination_that_is_a_source_is_refused(self):
        r = self._call(dst_wells=["a01"] + [f"{r}7" for r in "BCDEFGH"])
        assert isinstance(r, tuple) and r[1] == 400
        assert "A1" in r[0]["error"]

    def test_the_diluent_well_is_not_a_destination(self):
        r = self._call(dst_wells=[f"{r}7" for r in "ABCDEFGH"],
                       diluent_ref="plate:C7", target_conc=10,
                       target_ng=None)
        assert isinstance(r, tuple) and r[1] == 400
        assert "C7" in r[0]["error"]

    def test_defaulted_destinations_on_one_plate_are_refused(self):
        r = self._call(dst_labware="plate_96")
        assert isinstance(r, tuple) and r[1] == 400


# ══════════════════════════════════════════════════════════════════════════
# Notebook / gels
# ══════════════════════════════════════════════════════════════════════════

def test_a_saved_gel_keeps_a_pcr_lane_size():
    import splicecraft_gels as _gels
    out = _gels._normalise_gel_entry(
        {"name": "g", "lanes": [{"name": "PCR", "source": "pcr",
                                 "detail": "F+R", "_pcr_bp": 1234}]},
        fresh=True)
    assert out["lanes"][0]["_pcr_bp"] == 1234


def test_a_bogus_pcr_size_is_dropped():
    import splicecraft_gels as _gels
    out = _gels._normalise_gel_entry(
        {"name": "g", "lanes": [
            {"name": "a", "source": "pcr", "detail": "", "_pcr_bp": "junk"},
            {"name": "b", "source": "pcr", "detail": "", "_pcr_bp": -5}]},
        fresh=True)
    assert all("_pcr_bp" not in ln for ln in out["lanes"])


def test_an_attached_image_is_recorded_on_the_entry():
    import inspect
    src = inspect.getsource(sc._h_attach_experiment_image)
    assert 'entries[idx]["image_paths"]' in src, \
        "the file is written and linked but never recorded, so exports embed " \
        "nothing"


# ══════════════════════════════════════════════════════════════════════════
# Docs
# ══════════════════════════════════════════════════════════════════════════

class TestAgentEndpointIndexDoesNotDrift:
    """`docs/agent-api.md` is prose and documented about half the endpoints; the
    half it missed was invisible to anyone reading the docs. The generated index
    is the completeness guarantee, and this is what keeps it honest."""

    @staticmethod
    def _index_path():
        from pathlib import Path as _P
        return _P(__file__).parent.parent / "docs" / "agent-endpoints.md"

    def test_the_index_exists(self):
        assert self._index_path().is_file()

    def test_every_endpoint_appears(self):
        text = self._index_path().read_text()
        missing = sorted(n for n in sc._state._AGENT_HANDLERS
                         if f"| `{n}` |" not in text)
        assert not missing, (
            f"{len(missing)} endpoint(s) missing from docs/agent-endpoints.md "
            f"({', '.join(missing[:8])}…) — regenerate with "
            f"`python3 scripts/gen_agent_endpoint_index.py`")

    def test_it_lists_nothing_that_is_not_an_endpoint(self):
        import re as _re
        text = self._index_path().read_text()
        listed = set(_re.findall(r"^\| `([^`]+)` \|", text, _re.M))
        assert listed <= set(sc._state._AGENT_HANDLERS)


def test_the_docs_do_not_claim_reads_are_unauthenticated():
    from pathlib import Path as _P
    doc = (_P(__file__).parent.parent / "docs" / "agent-api.md").read_text()
    # The token check fires BEFORE the handler lookup, so reads need one too.
    assert "reads are\n  unauthenticated" not in doc
    assert "Bearer-token auth on EVERY endpoint" in doc


# ══════════════════════════════════════════════════════════════════════════
# Edits, re-origin, flip — what a location carries besides its numbers
# ══════════════════════════════════════════════════════════════════════════

_H10_GB_HEAD = """LOCUS       TESTPART                 300 bp    DNA     circular SYN 01-JAN-2026
DEFINITION  test.
ACCESSION   TESTPART
VERSION     TESTPART
DBLINK      BioProject: PRJNA000001
KEYWORDS    .
SOURCE      synthetic
  ORGANISM  synthetic
            .
FEATURES             Location/Qualifiers
     source          1..300
                     /organism="synthetic"
     CDS             <1..>90
                     /label="partialCDS"
                     /codon_start=1
     misc_feature    order(120..130,140..150)
                     /label="orderfeat"
     gene            200..>250
                     /label="partialgene"
ORIGIN
"""


def _h10_record():
    seq = ("ATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCG"
           "GGTTTAAACCCGGGTTTAAACCCGGGTTTAAA" + "C" * 210)
    assert len(seq) == 300
    lines = []
    for i in range(0, 300, 60):
        chunk = seq[i:i + 60].lower()
        lines.append(f"{i + 1:>9} "
                     + " ".join(chunk[j:j + 10] for j in range(0, 60, 10)))
    return sc._gb_text_to_record(_H10_GB_HEAD + "\n".join(lines) + "\n//\n")


def _locs(rec):
    return {f.qualifiers.get("label", [f.type])[0]: str(f.location)
            for f in rec.features}


class _EditRig:
    _current_record = None


def _edit(rec, mode, s, e, bases):
    seq = str(rec.seq)
    new = seq[:s] + bases + (seq[s:] if mode == "insert" else seq[e:])
    return sc.PlasmidApp._rebuild_record_with_edit(
        _EditRig(), new, mode, s, e, bases, source_record=rec)


_H10_LABELS = ("partialCDS", "orderfeat", "partialgene")


class TestEditsKeepWhatALocationCarries:
    """A base edit, a re-origin and a flip each rebuilt every location from
    plain integers: a 5'-partial CDS (`<1..>90`) came back as a complete ORF
    that starts and stops wherever the record happens to, a re-origin turned
    `order()` into `join()`, and all three dropped the record's DBLINK. The
    reference for each is Biopython's own slicing / reverse complement, which
    carries fuzzy positions correctly for any feature the cut does not cross."""

    def test_an_insert_matches_biopython_slicing(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = _h10_record()
        out = _edit(rec, "insert", 280, 280, "GGG")
        ref = rec[:280] + SeqRecord(Seq("GGG")) + rec[280:]
        got, want = _locs(out), _locs(ref)
        for label in _H10_LABELS:
            assert got[label] == want[label], label

    def test_an_upstream_insert_shifts_a_partial_end_and_keeps_it(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = _h10_record()
        out = _edit(rec, "insert", 100, 100, "GGGG")
        ref = rec[:100] + SeqRecord(Seq("GGGG")) + rec[100:]
        for label in ("orderfeat", "partialgene", "partialCDS"):
            assert _locs(out)[label] == _locs(ref)[label], label
        assert _locs(out)["partialgene"] == "[203:>254](+)"

    def test_a_replace_keeps_the_markers_on_both_sides(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = _h10_record()
        out = _edit(rec, "replace", 160, 170, "AA")
        ref = rec[:160] + SeqRecord(Seq("AA")) + rec[170:]
        for label in _H10_LABELS:
            assert _locs(out)[label] == _locs(ref)[label], label

    def test_dblink_survives_an_edit(self):
        out = _edit(_h10_record(), "insert", 280, 280, "GGG")
        assert out.dbxrefs == ["BioProject:PRJNA000001"]
        assert "DBLINK" in sc._record_to_gb_text(out)

    @pytest.mark.parametrize("where", [0, 300])
    def test_a_prepend_or_append_grows_the_source(self, where):
        out = _edit(_h10_record(), "insert", where, where, "GGGG")
        src = [f for f in out.features if f.type == "source"][0]
        assert (int(src.location.start), int(src.location.end)) == (0, 304)

    def test_a_reorigin_matches_biopython_and_keeps_dblink(self):
        rec = _h10_record()
        rot = sc._rotate_seq_record(rec, 100, keep_source=True)
        ref = rec[100:] + rec[:100]          # no feature below crosses bp 100
        got, want = _locs(rot), _locs(ref)
        for label in _H10_LABELS:
            assert got[label] == want[label], label
        assert got["orderfeat"].startswith("order{")
        assert rot.dbxrefs == ["BioProject:PRJNA000001"]

    def test_a_partial_feature_split_by_the_new_origin_keeps_both_ends(self):
        # `<1..>90` rotated so base 51 becomes base 1: the 5' partial end is now
        # at 251, the 3' partial end at 40, and the cut between them is exact.
        rot = sc._rotate_seq_record(_h10_record(), 50, keep_source=True)
        assert _locs(rot)["partialCDS"] == "join{[<250:300](+), [0:>40](+)}"

    def test_a_flip_matches_biopython_reverse_complement(self):
        rec = _h10_record()
        fl = sc._reverse_complement_record(rec)
        ref = rec.reverse_complement(id=True, name=True, description=True,
                                     features=True, annotations=True,
                                     dbxrefs=True)
        got, want = _locs(fl), _locs(ref)
        for label in _H10_LABELS:
            assert got[label] == want[label], label
        assert got["partialCDS"] == "[<210:>300](-)"
        assert fl.dbxrefs == ["BioProject:PRJNA000001"]

    def test_deleting_a_feature_keeps_dblink(self):
        rig = _EditRig()
        rig._current_record = _h10_record()
        out = sc.PlasmidApp._rebuild_record_without_feature(rig, 0)
        assert out.dbxrefs == ["BioProject:PRJNA000001"]


# ── H3 follow-up: the unsaved-changes guard must not re-run its caller ─────

class TestUnsavedGuardContract:
    """`_guard_unsaved_then` ran `action` on the clean path as well as after
    the modal. Every caller passes ITSELF as `action`, so a clean canvas
    recursed until RecursionError — opening a library row or importing a file
    crashed the app."""

    class _Rig:
        _current_record = object()

        def __init__(self, unsaved):
            self._unsaved = unsaved
            self.pushed = []

        def push_screen(self, screen, callback=None):
            self.pushed.append((screen, callback))

    def test_clean_means_proceed_without_running_the_action(self):
        rig = self._Rig(unsaved=False)
        ran = []
        assert sc.PlasmidApp._guard_unsaved_then(
            rig, "open x", lambda: ran.append(1)) is True
        assert ran == [] and rig.pushed == []

    def test_unsaved_defers_to_the_modal_which_runs_the_action(self):
        rig = self._Rig(unsaved=True)
        ran = []
        rig._discard_changes = lambda: setattr(rig, "_unsaved", False)
        assert sc.PlasmidApp._guard_unsaved_then(
            rig, "open x", lambda: ran.append(1)) is False
        (_screen, callback), = rig.pushed
        callback("discard")
        assert ran == [1] and rig._unsaved is False
        callback(None)                       # cancel: nothing more happens
        assert ran == [1]


# ══════════════════════════════════════════════════════════════════════════
# The other editor's .dna — our own export, read back
# ══════════════════════════════════════════════════════════════════════════

def _dna_packets(path):
    """(type, payload) for every packet of a .dna file — read straight off the
    bytes, so assertions about the FILE do not route through our reader."""
    import struct
    data = path.read_bytes()
    out, i = [], 0
    while i + 5 <= len(data):
        (length,) = struct.unpack(">I", data[i + 1:i + 5])
        out.append((data[i], data[i + 5:i + 5 + length]))
        i += 5 + length
    return out


def _dna_features_xml(path):
    import xml.etree.ElementTree as ET
    for t, payload in _dna_packets(path):
        if t == 0x0A:
            return ET.fromstring(payload.decode("utf-8").split("?>", 1)[-1])
    return None


# ── FM6: a .dna we wrote reads back as the plasmid we wrote ────────────────

class TestOwnDnaRoundTrip:
    """Exporting to `.dna` and opening the file again lost the plasmid's name
    (it came back named after the first WORD of its description — Biopython's
    reading of the Comments note), every custom feature colour (the writer
    stamped the type default), any text between angle brackets (read back as
    markup and deleted) and the arrowlessness of an arrowless feature (read
    back as forward). The file also declared the molecule SINGLE-stranded."""

    NOTE = "use < 5 ng and > 2 ng; Amp<sup>R</sup> & co"

    @classmethod
    def _entry(cls):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        r = SeqRecord(Seq("ACGTTGCAAC" * 60), id="pColor", name="pColor",
                      description="Cloning vector pColor, complete sequence")
        r.annotations.update(molecule_type="DNA", topology="circular")
        r.features = [
            SeqFeature(FeatureLocation(10, 200, 1), type="CDS", qualifiers={
                "label": ["GFP"], "ApEinfo_fwdcolor": ["#00ff00"],
                "ApEinfo_revcolor": ["#00ff00"], "note": [cls.NOTE]}),
            SeqFeature(FeatureLocation(250, 300, 0), type="misc_feature",
                       qualifiers={"label": ["arrowless"], "note": ["R&D"]}),
            SeqFeature(FeatureLocation(320, 380, -1), type="misc_feature",
                       qualifiers={"label": ["reverse"]}),
            SeqFeature(FeatureLocation(400, 460, 0), type="misc_feature",
                       qualifiers={"label": ["both"],
                                   "SpliceCraft_strand": ["double"]}),
        ]
        return {"id": "pcolor", "name": "pColor GFP v2",
                "gb_text": sc._record_to_gb_text(r)}

    def _round_trip(self, tmp_path, fname="export.dna"):
        out = tmp_path / fname
        sc._export_commercialsaas_dna(self._entry(), out)
        return out, sc.load_genbank(str(out))

    @staticmethod
    def _by_label(rec):
        return {f.qualifiers["label"][0]: f for f in rec.features
                if f.type != "source"}

    def test_the_name_survives_a_renamed_file(self, tmp_path):
        _, back = self._round_trip(tmp_path, "renamed.dna")
        assert sc.PlasmidApp._record_display_name(back) == "pColor GFP v2"

    def test_the_description_is_not_the_name(self, tmp_path):
        _, back = self._round_trip(tmp_path, "renamed.dna")
        assert "Cloning" not in (back.name, back.id)
        assert back.description == "Cloning vector pColor, complete sequence"

    def test_a_custom_colour_survives(self, tmp_path):
        out, back = self._round_trip(tmp_path)
        seg = _dna_features_xml(out).find(".//Feature[@name='GFP']/Segment")
        assert seg.get("color").lower() == "#00ff00"
        assert self._by_label(back)["GFP"].qualifiers[
            "ApEinfo_fwdcolor"] == ["#00ff00"]

    def test_text_that_looks_like_markup_survives(self, tmp_path):
        _, back = self._round_trip(tmp_path)
        by = self._by_label(back)
        assert by["GFP"].qualifiers["note"] == [self.NOTE]
        assert by["arrowless"].qualifiers["note"] == ["R&D"]

    def test_every_strand_survives(self, tmp_path):
        _, back = self._round_trip(tmp_path)
        by = self._by_label(back)
        assert by["GFP"].location.strand == 1
        assert by["reverse"].location.strand == -1
        assert by["arrowless"].location.strand in (0, None)
        assert by["both"].location.strand in (0, None)
        assert by["both"].qualifiers.get("SpliceCraft_strand") == ["double"]

    def test_directions_use_the_formats_own_vocabulary(self, tmp_path):
        # No attribute = no arrow, "3" = both — the reading an independent
        # third-party parser gives the same attribute.
        out, _ = self._round_trip(tmp_path)
        dirs = {f.get("name"): f.get("directionality")
                for f in _dna_features_xml(out).iter("Feature")}
        assert dirs == {"GFP": "1", "reverse": "2", "arrowless": None,
                        "both": "3"}

    def test_the_file_declares_double_stranded_dna(self, tmp_path):
        # Flag byte: 0x01 circular, 0x02 double-stranded, 0x04/0x08/0x10
        # Dam/Dcm/EcoKI. Every real file in tests/*.dna carries 0x1f.
        out, _ = self._round_trip(tmp_path)
        flags = next(p[0] for t, p in _dna_packets(out) if t == 0x00)
        assert flags == 0x1f
        real = next(p[0] for t, p in _dna_packets(
            __import__("pathlib").Path("tests/FFE 2 ENTRY A1.dna")) if t == 0x00)
        assert flags == real

    def test_control_bytes_in_the_description_do_not_break_the_file(
            self, tmp_path):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        r = SeqRecord(Seq("ACGT" * 60), id="P", name="P",
                      description="bad\x01byte￾ here")
        r.annotations.update(molecule_type="DNA", topology="circular")
        out = tmp_path / "ctl.dna"
        out.write_bytes(sc._write_commercialsaas_dna_bytes(r))
        back = sc.load_genbank(str(out))
        assert str(back.seq).upper() == "ACGT" * 60


class TestThirdPartyDnaConventions:
    """What the other editor writes, read the way it means it — checked
    against the fixture's own XML, not against our reader."""

    FIXTURE = "tests/FFE 2 ENTRY A1.dna"

    def test_a_feature_drawn_without_an_arrow_imports_without_one(self):
        from pathlib import Path
        root = _dna_features_xml(Path(self.FIXTURE))
        undirected = sorted(
            (f.get("type"), f.get("name")) for f in root.iter("Feature")
            if f.get("directionality") is None and f.get("type") != "source")
        assert undirected, "fixture changed: no undirected features left"
        rec = sc.load_genbank(self.FIXTURE)
        got = sorted((f.type, f.qualifiers.get("label", [""])[0])
                     for f in rec.features
                     if f.type != "source" and f.location.strand in (0, None))
        assert got == undirected

    def test_a_directed_feature_keeps_its_arrow(self):
        from pathlib import Path
        root = _dna_features_xml(Path(self.FIXTURE))
        want = sorted((f.get("type"), f.get("name"),
                       1 if f.get("directionality") == "1" else -1)
                      for f in root.iter("Feature")
                      if f.get("directionality") in ("1", "2"))
        rec = sc.load_genbank(self.FIXTURE)
        got = sorted((f.type, f.qualifiers.get("label", [""])[0],
                      f.location.strand)
                     for f in rec.features if f.location.strand in (1, -1))
        assert got == want

    def test_a_rich_note_is_decoded(self, tmp_path):
        # The editor keeps a rich note as an HTML document; an entity inside
        # it is the character it stands for.
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from xml.sax.saxutils import quoteattr
        r = SeqRecord(Seq("ACGT" * 50), id="P", name="P")
        r.annotations.update(molecule_type="DNA", topology="circular")
        r.features = [SeqFeature(FeatureLocation(0, 40, 1), type="CDS",
                                 qualifiers={"label": ["g"]})]
        data = sc._write_commercialsaas_dna_bytes(r)
        rich = ("<html><body>Kan<sup>R</sup> &amp; 5&#x27; &lt;UTR&gt;"
                "</body></html>")
        packets = []
        for t, _n, p in sc._iter_commercialsaas_packets(data):
            if t == 0x0A:
                p = p.decode().replace(
                    "</Feature>",
                    f'<Q name="note"><V text={quoteattr(rich)}/></Q>'
                    "</Feature>", 1).encode()
            packets.append(sc._build_commercialsaas_packet(t, p))
        f = tmp_path / "rich.dna"
        f.write_bytes(b"".join(packets))
        cds = [x for x in sc.load_genbank(str(f)).features
               if x.type == "CDS"][0]
        assert cds.qualifiers["note"] == ["KanR & 5' <UTR>"]


# ── S3: a history node is its whole subtree ────────────────────────────────

def _hnode(name, bp, op="insertFragment"):
    return sc._CommercialSaaSHistoryNode.new(
        name=name, seq_len=bp, circular=True, operation=op)


def _hasm(name, bp, parents):
    node = _hnode(name, bp)
    for p in parents:
        node.add_parent(p)
    return node


class TestHistoryIdentityIsTheSubtree:
    """The viewer and the protocol collapse a repeated sub-assembly. Keying
    that on name + length + operation hid every later link of a same-length
    edit chain; keying it on the node ID as well (the first fix) stopped the
    collapse altogether, because every history we build is RENUMBERED — so
    these trees are renumbered before they are asked anything."""

    def _chain(self):
        node = _hnode("pX.dna", 3000, "invalid")
        for _ in range(3):
            node = _hasm("pX.dna", 3000, [node])      # same name, same length
        sc._history_renumber_node_ids(node)
        return node

    def _reused(self):
        def tu():
            return _hasm("TU.dna", 4000, [_hnode("pENTR.dna", 2000),
                                          _hnode("prom.dna", 500)])
        mod = _hasm("MOD.dna", 9000, [_hnode("pENTR2.dna", 2500), tu(), tu()])
        sc._history_renumber_node_ids(mod)
        return mod

    def test_every_link_of_an_edit_chain_is_its_own_step(self):
        assert len(sc._history_build_steps(self._chain())) == 3

    def test_edit_chain_links_never_share_a_signature(self):
        node, sigs = self._chain(), []
        while node is not None:
            sigs.append(sc._history_node_signature(node))
            node = (node.parents or [None])[0]
        assert len(set(sigs)) == len(sigs) == 4

    def test_a_reused_subassembly_is_still_one_step(self):
        steps = sc._history_build_steps(self._reused())
        assert [s["product"] for s in steps] == ["TU", "MOD"]

    def test_copies_match_although_their_ids_differ(self):
        _bb, tu1, tu2 = self._reused().parents
        assert tu1.node_id != tu2.node_id
        assert (sc._history_node_signature(tu1)
                == sc._history_node_signature(tu2))


# ── FM11: the COMMENT name marker owns its own line, nothing after it ──────

class TestDisplayNameMarkerOwnsOnlyItsLine:
    """The display name travels in a `SpliceCraft-name:` COMMENT line. The
    restore took EVERYTHING after the marker as the name, so a line another
    tool appended after ours ("Sequence verified …") became part of the
    plasmid's name — and the next save cut it out of the COMMENT."""

    LONG = ("A really quite long plasmid display name that will certainly "
            "wrap around the column limit")

    @staticmethod
    def _rec(comment, disp=None):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        r = SeqRecord(Seq("ACGT" * 50), id="My_Plasmid", name="My_Plasmid",
                      description="d")
        r.annotations.update(molecule_type="DNA", topology="circular",
                             comment=comment)
        if disp is not None:
            r._tui_display_name = disp
        return r

    def _third_party_after(self, name):
        # What another tool produces: our file, one line appended after ours.
        text = sc._record_to_gb_text(self._rec("", disp=name))
        return self._append_comment(
            text, "Sequence verified 2026-03-02 by Sanger")

    @staticmethod
    def _append_comment(gb_text, line):
        out, done = [], False
        rows = gb_text.split("\n")
        for i, row in enumerate(rows):
            out.append(row)
            nxt = rows[i + 1] if i + 1 < len(rows) else ""
            in_comment = row.startswith(("COMMENT", " " * 12)) and any(
                r.startswith("COMMENT") for r in rows[:i + 1])
            if (not done and in_comment and not nxt.startswith(" " * 12)):
                out.append(" " * 12 + line)
                done = True
        return "\n".join(out)

    @pytest.mark.parametrize("name", ["My Plasmid", LONG])
    def test_a_following_line_is_not_part_of_the_name(self, name):
        back = sc._gb_text_to_record(self._third_party_after(name), cache=False)
        assert back._tui_display_name == name

    @pytest.mark.parametrize("name", ["My Plasmid", LONG])
    def test_a_following_line_survives_the_next_save(self, name):
        back = sc._gb_text_to_record(self._third_party_after(name), cache=False)
        again = sc._gb_text_to_record(sc._record_to_gb_text(back), cache=False)
        assert "Sequence verified 2026-03-02 by Sanger" in \
            again.annotations["comment"]
        assert again._tui_display_name == name

    @pytest.mark.parametrize("name", ["My Plasmid", LONG])
    def test_our_own_round_trip_is_stable(self, name):
        once = sc._gb_text_to_record(
            sc._record_to_gb_text(self._rec("", disp=name)), cache=False)
        assert once._tui_display_name == name
        text = sc._record_to_gb_text(once)
        twice = sc._gb_text_to_record(text, cache=False)
        assert twice._tui_display_name == name
        assert sc._record_to_gb_text(twice) == text
        assert twice.annotations["comment"] == once.annotations["comment"]


# ── FM11 / D15: an existing blob is verified, and its age refreshed ────────

class TestBlobReuseIsVerified:
    """`_blob_write` returned early whenever a blob file already existed, so a
    truncated / bit-rotted blob was kept (and the good text discarded) and
    every later read reported the sequence unavailable. Reusing a blob also
    left its mtime alone, which is the only thing the orphan GC's grace
    window looks at."""

    TEXT = "LOCUS       X  10 bp DNA linear\nORIGIN\n        1 acgtacgtac\n//\n"

    @staticmethod
    def _forget_verification():
        getattr(sc, "_BLOB_VERIFIED_PATHS", set()).clear()

    @classmethod
    def _sandboxed(cls, monkeypatch, tmp_path):
        import splicecraft_state as _st
        monkeypatch.setattr(_st, "_PLASMID_BLOB_DIR", tmp_path / "blobs")
        cls._forget_verification()
        return sc._blob_path(sc._blob_hash_bytes(cls.TEXT.encode()))

    def test_a_truncated_blob_is_rewritten(self, monkeypatch, tmp_path):
        path = self._sandboxed(monkeypatch, tmp_path)
        ref = sc._blob_write(self.TEXT)
        path.write_bytes(self.TEXT.encode()[:5])        # damage it
        self._forget_verification()
        assert sc._blob_write(self.TEXT) == ref
        assert sc._blob_read(ref) == self.TEXT

    def test_same_size_corruption_is_caught_once_per_session(
            self, monkeypatch, tmp_path):
        path = self._sandboxed(monkeypatch, tmp_path)
        ref = sc._blob_write(self.TEXT)
        self._forget_verification()
        bad = bytearray(self.TEXT.encode())
        bad[-3] ^= 0x01                                   # one flipped bit
        path.write_bytes(bytes(bad))
        assert sc._blob_write(self.TEXT) == ref
        assert sc._blob_read(ref) == self.TEXT

    def test_reuse_refreshes_the_gc_grace_window(self, monkeypatch, tmp_path):
        import os
        import time
        path = self._sandboxed(monkeypatch, tmp_path)
        sc._blob_write(self.TEXT)
        old = time.time() - 10 * sc._BLOB_GC_GRACE_SECONDS
        os.utime(path, (old, old))
        sc._blob_write(self.TEXT)
        assert time.time() - path.stat().st_mtime < sc._BLOB_GC_GRACE_SECONDS


# ── R5: a LOCAL alignment reports positions in the plasmid's frame ─────────

class TestLocalModePositions:
    """A local alignment's rows held only the aligned stretch, and every
    consumer counts columns from base 0 — so a SNP at bp 1400, in a read that
    aligned from bp 1000, came back at bp 400 in `mode: "local"`, with the
    uncovered span shifted to match. Both modes must agree, and agree with the
    base that was actually changed."""

    @staticmethod
    def _case():
        import splicecraft_state as st
        rng = random.Random(5)
        plasmid = "".join(rng.choice("ACGT") for _ in range(3000))
        sub = {"A": "C", "C": "G", "G": "T", "T": "A"}
        read = list(plasmid[1000:1800])
        read[400] = sub[read[400]]              # plasmid bp 1400
        return st._AGENT_HANDLERS["verify-against-reads"][0], plasmid, \
            "".join(read)

    @pytest.mark.parametrize("mode", ["global", "local"])
    def test_the_snp_is_where_it_was_made(self, mode):
        fn, plasmid, read = self._case()
        r = fn(None, {"reference": plasmid, "reads": [read], "circular": True,
                      "mode": mode})
        snps = [v for v in r["consensus"]["variants"] if v["type"] == "snp"]
        assert [v["target_pos"] for v in snps] == [1400], mode
        assert r["consensus"]["covered_bp"] == 800

    def test_local_and_global_agree_on_what_was_not_read(self):
        fn, plasmid, read = self._case()
        spans = {m: [tuple(s) for s in fn(None, {
                    "reference": plasmid, "reads": [read], "circular": True,
                    "mode": m})["consensus"]["uncovered_spans"]]
                 for m in ("global", "local")}
        assert spans["local"] == spans["global"]

    def test_local_rows_carry_the_target_and_soft_clip_the_read(self):
        # The whole TARGET (so a column maps to its plasmid bp), and only the
        # ALIGNED stretch of the read: its clipped ends are soft clips. Put in
        # as overhang against target gaps (the first fix), a clipped prefix
        # read as an insertion before bp 0 (harden 2026-09-24).
        q = "TTTTTTTT" + "ACGTACGGACGTAAACCC" + "GGGGGG"
        t = "CCCCCCCCCCCC" + "ACGTACGTACGTAAACCC" + "AAAAAAAAAA"
        r = sc._pairwise_align(q, t, mode="local")
        assert r["aligned_q"].replace("-", "") == q[8:26]
        assert r["aligned_t"].replace("-", "") == t
        assert r["local_t_span"] == [12, 30] and r["local_q_span"] == [8, 26]
        # The identity is still the LOCAL one: 17 of 18 aligned columns.
        assert r["n_matches"] == 17 and r["n_mismatches"] == 1
        variants = sc._extract_variants_from_alignment(
            r["aligned_q"], r["aligned_t"])
        assert [v["target_pos"] for v in variants if v["type"] == "snp"] \
            == [12 + 7]
        assert not [v for v in variants if v["type"] == "insertion"]


# ── D15a / D3: restoring must never delete, and must stick ─────────────────

def _d_ent(name, seed):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    s = "".join(random.Random(seed).choice("ACGT") for _ in range(120))
    gb = sc._record_to_gb_text(SeqRecord(
        Seq(s), id=name, name=name, description=name,
        annotations={"molecule_type": "DNA", "topology": "circular"}))
    return {"id": name, "name": name, "size": 120, "n_feats": 0, "gb_text": gb}


class TestRestoringASpillMergesIt:
    """A lost-entries spill holds only what one shrinking save dropped.
    Restored like a full backup, it REPLACED the file — the library shrank to
    the spilled entries and everything else went."""

    def test_the_spill_comes_back_and_nothing_else_leaves(self):
        import splicecraft_state as st
        from pathlib import Path
        lib = [_d_ent(f"P{i}", i) for i in range(10)]
        sc._save_library(lib)
        sc._save_library(lib[:3])            # a 7-of-10 shrink -> spilled
        spill = next(r for r in sc._list_recoverable_backups(st._LIBRARY_FILE)
                     if r["kind"] == "lost_entries")
        sc._save_library(sc._load_library() + [_d_ent("NEW", 99)])
        n = sc._restore_from_backup(st._LIBRARY_FILE,
                                    Path(spill["source_path"]),
                                    "Plasmid library")
        st._library_cache = None
        names = sorted(e["name"] for e in sc._load_library())
        assert names == sorted([f"P{i}" for i in range(10)] + ["NEW"])
        assert n == 7

    def test_restoring_the_same_spill_twice_adds_nothing(self):
        import splicecraft_state as st
        from pathlib import Path
        lib = [_d_ent(f"Q{i}", 50 + i) for i in range(6)]
        sc._save_library(lib)
        sc._save_library(lib[:2])
        spill = Path(next(r for r in sc._list_recoverable_backups(
            st._LIBRARY_FILE) if r["kind"] == "lost_entries")["source_path"])
        sc._restore_from_backup(st._LIBRARY_FILE, spill, "Plasmid library")
        st._library_cache = None
        assert sc._restore_from_backup(st._LIBRARY_FILE, spill,
                                       "Plasmid library") == 0
        st._library_cache = None
        assert len(sc._load_library()) == 6

    def test_a_spilled_collection_merges_into_its_namesake(self):
        cur = [{"name": "A", "plasmids": [{"id": "x"}, {"id": "y"}]}]
        spilled = [{"name": "A", "plasmids": [{"id": "x"}, {"id": "gone"}]},
                   {"name": "B", "plasmids": [{"id": "b1"}]}]
        merged, n = sc._merge_spilled_entries(cur, spilled)
        assert [c["name"] for c in merged] == ["A", "B"]
        assert [p["id"] for p in merged[0]["plasmids"]] == ["x", "y", "gone"]
        assert n == 2


class TestRestoredCollectionsStick:
    """Restoring collections.json put the recovered plasmid back — and the
    next ordinary library save, still holding the pre-restore mirror, wrote
    it straight back out of the collection."""

    class _App:
        _unsaved = False

        def call_from_thread(self, fn, *a, **k):
            pass

    def test_the_next_library_save_keeps_the_restored_entry(self):
        import json
        import splicecraft_state as st
        monkey_sync = getattr(sc, "_collection_sync_force_sync", None)
        sc._collection_sync_force_sync = True
        try:
            a = [_d_ent("A1", 1), _d_ent("KEEP-ME", 2)]
            sc._save_collections([{"name": "A", "plasmids": a}])
            sc._set_active_collection_name("A")
            sc._save_library(list(a))
            sc._save_library([a[0]])                 # deleted by mistake

            def colls():
                raw = json.loads(st._COLLECTIONS_FILE.read_text())["entries"]
                return {c["name"]: sorted(p["name"] for p in c["plasmids"])
                        for c in raw}

            rows = st._AGENT_HANDLERS["list-backups"][0](
                self._App(), {"label": "collections"})
            src = next(r["source_path"] for r in rows["backups"]
                       if "KEEP-ME" in open(r["source_path"]).read())
            st._AGENT_HANDLERS["restore-backup"][0](
                self._App(), {"label": "collections", "source_path": src})
            assert colls()["A"] == ["A1", "KEEP-ME"]
            lib = sc._load_library()
            lib[0]["status"] = "CLONING"
            sc._save_library(lib)                    # the next ordinary save
            assert colls()["A"] == ["A1", "KEEP-ME"]
        finally:
            sc._collection_sync_force_sync = monkey_sync


# ── D9 / D10: work committed to a collection survives the next save ────────

class _MirrorRig:
    _unsaved = False

    def call_from_thread(self, fn, *a, **k):
        pass


def _colls_by_name():
    import json
    import splicecraft_state as st
    raw = json.loads(st._COLLECTIONS_FILE.read_text())["entries"]
    return {c["name"]: sorted(p["name"] for p in c["plasmids"]) for c in raw}


class TestACommittedMoveSurvivesAMirrorFailure:
    """A move writes collections.json (the source of truth) and then re-stages
    the library mirror. When that second write failed, the moved plasmid sat
    in the collection while the cache still held the pre-move view — and the
    next ordinary library save wrote that view back over the collection. A
    plasmid moved INTO the active collection was then in no collection."""

    @pytest.fixture
    def _setup(self, monkeypatch):
        import splicecraft_persistence as pe
        import splicecraft_state as st
        monkeypatch.setattr(sc, "_collection_sync_force_sync", True)
        active = [_d_ent("K1", 11)]
        sc._save_collections([{"name": "Active", "plasmids": active},
                              {"name": "Other",
                               "plasmids": [_d_ent("MOVE-ME", 12)]}])
        sc._set_active_collection_name("Active")
        sc._save_library(list(active))
        real = pe._safe_save_json

        def fail_library(path, *a, **k):
            if path == st._LIBRARY_FILE:
                raise OSError(28, "No space left on device (injected)")
            return real(path, *a, **k)

        return lambda: monkeypatch.setattr(pe, "_safe_save_json", fail_library), \
            lambda: monkeypatch.setattr(pe, "_safe_save_json", real)

    @pytest.mark.parametrize("endpoint,payload", [
        ("move-plasmid", {"name": "MOVE-ME", "to": "Active"}),
        ("move-plasmids", {"names": ["MOVE-ME"], "to": "Active"}),
    ])
    def test_the_next_save_keeps_the_moved_plasmid(self, _setup, endpoint,
                                                   payload):
        import splicecraft_state as st
        break_it, heal_it = _setup
        break_it()
        r = st._AGENT_HANDLERS[endpoint][0](_MirrorRig(), payload)
        heal_it()
        assert isinstance(r, dict) and r.get("ok"), r
        assert r.get("warnings"), "a partial failure went unreported"
        assert "MOVE-ME" in _colls_by_name()["Active"]
        lib = sc._load_library()                  # the in-session view
        assert "MOVE-ME" in [e["name"] for e in lib]
        lib[0]["status"] = "CLONING"
        sc._save_library(lib)                     # the next ordinary save
        assert "MOVE-ME" in _colls_by_name()["Active"]


class TestNoActiveCollectionKeepsTheLibraryFile:
    """With the active collection deleted, a library save mirrors into
    nothing — and still cleared the dirty marker, so the next launch
    overwrote the library from whichever collection became active and a
    fresh import was gone."""

    def test_an_import_survives_the_relaunch(self, monkeypatch):
        import json
        import splicecraft_state as st
        monkeypatch.setattr(sc, "_collection_sync_force_sync", True)
        monkeypatch.setenv("SPLICECRAFT_SKIP_SETTINGS_FLUSH", "1")
        sc._save_collections([
            {"name": "Keep", "plasmids": [_d_ent("K1", 21)]},
            {"name": "Scratch", "plasmids": [_d_ent("S1", 22)]}])
        sc._set_active_collection_name("Scratch")
        sc._save_library([_d_ent("S1", 22)])
        st._AGENT_HANDLERS["delete-collection"][0](_MirrorRig(),
                                                   {"name": "Scratch"})
        lib = sc._load_library()
        lib.insert(0, _d_ent("NEW-IMPORT", 23))
        sc._save_library(lib)
        for attr in ("_library_cache", "_collections_cache",
                     "_settings_cache"):
            setattr(st, attr, None)
        sc._ensure_default_collection()          # the launch sequence
        sc._restore_library_from_active_collection()
        on_disk = [p["name"] for p in
                   json.loads(st._LIBRARY_FILE.read_text())["entries"]]
        assert "NEW-IMPORT" in on_disk


def test_primers_saved_with_no_active_primer_collection_survive_a_relaunch(
        monkeypatch):
    import json
    import splicecraft_state as st
    monkeypatch.setenv("SPLICECRAFT_SKIP_SETTINGS_FLUSH", "1")
    p1 = {"name": "P1", "sequence": "ACGTACGTACGTACGTACGT", "tm": 50.0}
    new = {"name": "AGENT-F", "sequence": "GGGCCCAAATTTGGGCCCAA", "tm": 55.0}
    sc._save_primer_collections([{"name": "Main", "primers": [p1]}])
    sc._set_active_primer_collection_name("Main")
    sc._save_primers([p1])
    sc._set_active_primer_collection_name("")
    sc._save_primers([p1, new])                  # mirrors into nothing
    for attr in ("_primers_cache", "_primer_collections_cache",
                 "_settings_cache"):
        setattr(st, attr, None)
    sc._ensure_default_primer_collection()       # the launch sequence
    sc._restore_primers_from_active_primer_collection()
    colls = {c["name"]: [p["name"] for p in c["primers"]]
             for c in json.loads(
                 st._PRIMER_COLLECTIONS_FILE.read_text())["entries"]}
    assert colls["Recovered primers"] == ["AGENT-F"]
    assert colls["Main"] == ["P1"]                # untouched


# ── primer goto: switch collections atomically, and ask first ──────────────

def _goto_fixture():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    rng = random.Random(7)

    def rec(rid, n):
        r = SeqRecord(Seq("".join(rng.choice("ACGT") for _ in range(n))),
                      id=rid, name=rid, description=rid,
                      annotations={"molecule_type": "DNA",
                                   "topology": "circular"})
        r.features = [SeqFeature(FeatureLocation(5, 50, strand=1),
                                 type="primer_bind",
                                 qualifiers={"label": ["pr"]})]
        return r

    def ent(r, name):
        return {"id": r.id, "name": name, "size": len(r.seq), "n_feats": 1,
                "source": "t", "added": "2025-01-01",
                "gb_text": sc._record_to_gb_text(r)}

    a, b, p2 = rec("PLASMIDA", 400), rec("PLASMIDB", 500), rec("PLASMIDP2", 600)
    ea, eb, ep2 = ent(a, "Plasmid A"), ent(b, "Plasmid B"), ent(p2, "Plasmid P2")
    sc._save_collections([{"name": "C1", "plasmids": [ea, eb]},
                          {"name": "C2", "plasmids": [ep2]}])
    sc._set_active_collection_name("C1")
    sc._safe_save_json_mirror(sc._state._LIBRARY_FILE, [ea, eb],
                              "Plasmid library")
    sc._state._library_cache = None
    return a, p2


class TestPrimerGotoSwitchesAtomically:
    """Jumping from the primer library to a plasmid in ANOTHER collection
    moved only the active pointer. The library mirror still held the old
    collection, so the next save wrote C1's plasmids into C2 — over C2's own.
    And it replaced an edited canvas without asking."""

    async def _run(self, dirty):
        from splicecraft_modals import UnsavedNavigateModal
        a, p2 = _goto_fixture()
        sc.PlasmidApp._skip_seed = True
        sc.PlasmidApp._skip_snapshot = True
        app = sc.PlasmidApp()
        out = {}
        async with app.run_test(size=(200, 50)) as pilot:
            for _ in range(6):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                await pilot.pause()
            app._AUTOSAVE_DEBOUNCE_S = 10_000.0
            app._apply_record(a)
            await pilot.pause()
            if dirty:
                app._edit_dialog_result(("GGG", "insert", 100, 100))
                await pilot.pause()
            app._goto_primer_in_plasmid_inner({
                "collection": "C2", "plasmid_id": p2.id, "start": 5, "end": 50})
            for _ in range(6):
                await pilot.pause()
            out["modal"] = isinstance(app.screen, UnsavedNavigateModal)
            out["canvas"] = getattr(app._current_record, "id", None)
            out["active"] = sc._get_active_collection_name()
            out["library"] = sorted(e["id"] for e in sc._load_library())
        return out

    async def test_a_clean_canvas_switches_the_whole_library(self):
        out = await self._run(dirty=False)
        assert out["canvas"] == "PLASMIDP2" and out["active"] == "C2"
        assert out["library"] == ["PLASMIDP2"]
        sc._save_library(sc._load_library())      # the next ordinary save
        assert _colls_by_name() == {"C1": ["Plasmid A", "Plasmid B"],
                                    "C2": ["Plasmid P2"]}

    async def test_an_edited_canvas_is_offered_first(self):
        out = await self._run(dirty=True)
        assert out["modal"], "the edited canvas was replaced without asking"
        assert out["canvas"] == "PLASMIDA" and out["active"] == "C1"


# ── H3: every path that replaces the canvas offers unsaved edits first ─────

def _h3_row_select(app):
    app.post_message(sc.LibraryPanel.PlasmidLoad(
        sc._find_library_entry_by_id("PLASMIDB")))


def _h3_find_plasmid(app):
    cap = {}
    orig = app.push_screen
    app.push_screen = lambda screen, callback=None, *a, **k: cap.setdefault(
        "cb", callback)
    try:
        app.action_find_plasmid()
    finally:
        app.push_screen = orig
    cap["cb"](("C1", "PLASMIDB"))


def _h3_focus_saved(app):
    b = sc._gb_text_to_record(sc._find_library_entry_by_id("PLASMIDB")["gb_text"])
    app._focus_saved_library_plasmid(b, "Plasmid B")


def _h3_primer_goto(app):
    app._goto_primer_in_plasmid({"collection": "C1", "plasmid_id": "PLASMIDB",
                                 "start": 5, "end": 50})


class TestEveryLoadAsksFirst:
    """Row-select and import asked; a find-plasmid pick, a post-save focus
    and a primer goto replaced an edited canvas outright."""

    async def _drive(self, trigger, *, dirty, answer=None):
        from splicecraft_modals import UnsavedNavigateModal
        a, _p2 = _goto_fixture()
        sc.PlasmidApp._skip_seed = True
        sc.PlasmidApp._skip_snapshot = True
        app = sc.PlasmidApp()
        out = {}
        async with app.run_test(size=(200, 50)) as pilot:
            for _ in range(6):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                await pilot.pause()
            app._AUTOSAVE_DEBOUNCE_S = 10_000.0
            app._apply_record(a)
            await pilot.pause()
            if dirty:
                app._edit_dialog_result(("GGG", "insert", 100, 100))
                await pilot.pause()
            trigger(app)
            for _ in range(8):
                await pilot.pause()
            out["modal"] = isinstance(app.screen, UnsavedNavigateModal)
            out["canvas"] = getattr(app._current_record, "id", None)
            out["len"] = len(app._current_record.seq)
            if out["modal"] and answer is not None:
                app.screen._dismiss_once(answer)
                for _ in range(8):
                    await pilot.pause()
                out["after"] = getattr(app._current_record, "id", None)
        return out

    TRIGGERS = [_h3_row_select, _h3_find_plasmid, _h3_focus_saved,
                _h3_primer_goto]

    @pytest.mark.parametrize("trigger", TRIGGERS,
                             ids=lambda f: f.__name__[4:])
    async def test_an_edited_canvas_is_offered_first(self, trigger):
        out = await self._drive(trigger, dirty=True, answer="discard")
        assert out["modal"], "the edited canvas was replaced without asking"
        assert (out["canvas"], out["len"]) == ("PLASMIDA", 403)
        assert out["after"] == "PLASMIDB"

    @pytest.mark.parametrize("trigger", TRIGGERS,
                             ids=lambda f: f.__name__[4:])
    async def test_a_clean_canvas_just_loads(self, trigger):
        out = await self._drive(trigger, dirty=False)
        assert not out["modal"] and out["canvas"] == "PLASMIDB"

    async def test_open_file_asks_before_the_picker(self):
        from splicecraft_modals import UnsavedNavigateModal
        a, _p2 = _goto_fixture()
        sc.PlasmidApp._skip_seed = True
        sc.PlasmidApp._skip_snapshot = True
        app = sc.PlasmidApp()
        async with app.run_test(size=(200, 50)) as pilot:
            for _ in range(6):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                await pilot.pause()
            app._apply_record(a)
            await pilot.pause()
            app._edit_dialog_result(("GGG", "insert", 100, 100))
            await pilot.pause()
            app.action_open_file()
            await pilot.pause()
            assert isinstance(app.screen, UnsavedNavigateModal)

    async def test_an_unguarded_switch_keeps_the_edits_reachable(self):
        # The net under the prompts: switch away WITHOUT asking (a direct
        # `_apply_record`, as a path that bypasses the guard would) and come
        # back — redo restores the unsaved edit.
        a, p2 = _goto_fixture()
        sc.PlasmidApp._skip_seed = True
        sc.PlasmidApp._skip_snapshot = True
        app = sc.PlasmidApp()
        async with app.run_test(size=(200, 50)) as pilot:
            for _ in range(6):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                await pilot.pause()
            app._AUTOSAVE_DEBOUNCE_S = 10_000.0
            app._apply_record(a)
            await pilot.pause()
            app._edit_dialog_result(("GGG", "insert", 100, 100))
            await pilot.pause()
            app._apply_record(p2)
            await pilot.pause()
            app._apply_record(a)
            await pilot.pause()
            assert len(app._current_record.seq) == 400
            app._action_redo()
            await pilot.pause()
            assert len(app._current_record.seq) == 403


# ── AA5: a routing key an endpoint would ignore is refused, centrally ──────

class TestRoutingKeysAreReadOrRefused:
    """Only a handful of write endpoints refused a routing key they don't read,
    so `delete-part {"parts_bin": "Archive"}` answered 200 and deleted from
    the ACTIVE bin. The dispatcher now derives each handler's readable keys
    from its own code (following `payload` into helpers) and refuses the rest
    — and leaves any handler it cannot fully account for alone."""

    def _parts(self):
        import splicecraft_state as st
        part = {"name": "P1", "sequence": "ACGT" * 10, "grammar": "gb_l0",
                "part_type": "CDS", "oh5": "AATG", "oh3": "GCTT"}
        sc._save_parts_bin_collections([
            {"name": "Active", "parts": [dict(part)]},
            {"name": "Archive", "parts": [dict(part)]}])
        sc._set_active_parts_bin_name("Active")
        sc._save_parts_bin([dict(part)])
        return st

    def test_an_unread_bin_key_is_refused_and_nothing_is_deleted(self):
        st = self._parts()
        payload, status = sc._agent_invoke(
            None, "delete-part", {"name": "P1", "parts_bin": "Archive"})
        assert status == 400 and payload["unsupported"] == ["parts_bin"]
        bins = {b["name"]: [p["name"] for p in b["parts"]]
                for b in sc._load_parts_bin_collections()}
        assert bins == {"Active": ["P1"], "Archive": ["P1"]}

    def test_a_key_the_handler_reads_is_not_refused(self):
        # delete-part reads `bin` — the spelling it documents.
        self._parts()

        class _App:
            _unsaved = False

            def call_from_thread(self, fn, *a, **k):
                return fn(*a, **k)

        payload, status = sc._agent_invoke(
            _App(), "delete-part", {"name": "P1", "bin": "Nope"})
        assert status == 404, payload             # reached the handler
        assert "Nope" in payload["error"]

    def test_the_analysis_agrees_with_every_handlers_own_gets(self):
        # Every `payload.get("k")` literal in a write handler's own source must
        # be in the key set derived for it — a dropped read would turn into a
        # refusal of a key the handler honours.
        import inspect
        import re
        import splicecraft_agent as ag
        missing = []
        for name, (fn, write) in sc._state._AGENT_HANDLERS.items():
            if not write:
                continue
            keys = ag._agent_payload_keys(fn)
            if keys is None:
                continue
            src = inspect.getsource(fn)
            for k in re.findall(r"""payload\.get\(\s*["']([A-Za-z_0-9]+)["']""",
                                src):
                if k not in keys:
                    missing.append((name, k))
        assert not missing, missing[:10]

    def test_most_write_endpoints_are_covered(self):
        import splicecraft_agent as ag
        writes = [fn for fn, w in sc._state._AGENT_HANDLERS.values() if w]
        known = [fn for fn in writes if ag._agent_payload_keys(fn) is not None]
        assert len(known) >= 0.9 * len(writes)


def test_a_read_only_guest_may_shut_itself_down(monkeypatch):
    import splicecraft_state as st
    monkeypatch.setattr(st, "_AGENT_READ_ONLY", True)
    monkeypatch.setattr(st, "_AGENT_READ_ONLY_HOLDER_PID", 4242)
    monkeypatch.setattr(sc, "_agent_schedule_shutdown", lambda app: None)

    class _Guest:
        _headless = True
        _unsaved = False

        def call_from_thread(self, fn, *a, **k):
            pass

    payload, status = sc._agent_invoke(_Guest(), "shutdown", {})
    assert status == 200 and payload.get("shutting_down") is True


# ── FM10: a line break inside a value never reaches the stored text ────────

class TestEveryLineBreakSplitsAValue:
    """Only "\\n" was split out of qualifier values. A lone CR (GFF3 `%0D`, a
    `.dna` `&#13;`) or a VT went into the stored GenBank text raw; a bulk
    export copied it verbatim, and nothing could read the file — while the
    export reported 0 failures."""

    def _rec(self, note):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        r = SeqRecord(Seq("ACGT" * 30), id="P", name="P")
        r.annotations.update(molecule_type="DNA", topology="circular")
        r.features = [SeqFeature(FeatureLocation(0, 40, 1),
                                 type="misc_feature",
                                 qualifiers={"note": [note]})]
        return r

    @pytest.mark.parametrize("brk", ["\r", "\x0b", "\x0c", " "])
    def test_the_stored_text_reads_back(self, brk):
        import io
        from Bio import SeqIO
        text = sc._record_to_gb_text(self._rec(f"first{brk}second"))
        assert brk not in text
        back = SeqIO.read(io.StringIO(text), "genbank")
        assert back.features[0].qualifiers["note"] == ["first", "second"]

    def test_the_bulk_fast_path_refuses_control_characters(self):
        ok = sc._record_to_gb_text(self._rec("plain note"))
        assert sc._gb_text_is_spec_clean(ok)
        for bad in ("\r", "\t", "\x0b", "\x00"):
            assert not sc._gb_text_is_spec_clean(
                ok.replace("plain note", f"plain{bad}note"))


# ── FM11: files that parse into something unusable are refused cleanly ─────

_CONTIG_ONLY = """LOCUS       CTG1                    5000 bp    DNA     linear   CON 01-JAN-2026
DEFINITION  scaffold.
ACCESSION   CTG1
VERSION     CTG1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     source          1..5000
CONTIG      join(ABC000001.1:1..5000)
//
"""


class TestUnusableRecordsAreRefused:
    def test_a_contig_only_record_says_so(self, tmp_path):
        with pytest.raises(ValueError, match="CONTIG"):
            sc._gb_text_to_record(_CONTIG_ONLY, cache=False)
        f = tmp_path / "ctg.gb"
        f.write_text(_CONTIG_ONLY)
        with pytest.raises(ValueError, match="CONTIG"):
            sc.load_genbank(str(f))

    def _dna_with_segment(self, tmp_path, segment_attrs):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        r = SeqRecord(Seq("ACGT" * 50), id="P", name="P")
        r.annotations.update(molecule_type="DNA", topology="circular")
        r.features = [SeqFeature(FeatureLocation(0, 10, 1), type="misc_feature",
                                 qualifiers={"label": ["m"]})]
        data = sc._write_commercialsaas_dna_bytes(r)
        out = []
        for t, _n, pl in sc._iter_commercialsaas_packets(data):
            if t == 0x0A:
                pl = pl.decode().replace('range="1-10"', segment_attrs).encode()
            out.append(sc._build_commercialsaas_packet(t, pl))
        f = tmp_path / "seg.dna"
        f.write_bytes(b"".join(out))
        return f

    def test_a_zero_based_segment_is_clamped_not_exported_broken(
            self, tmp_path):
        f = self._dna_with_segment(tmp_path, 'range="0-10"')
        rec = sc.load_genbank(str(f))
        m = [x for x in rec.features if x.type == "misc_feature"][0]
        assert (int(m.location.start), int(m.location.end)) == (0, 10)
        sc._export_genbank_to_path(rec, tmp_path / "ok.gb")   # no refusal

    def test_a_segment_without_a_range_is_a_clean_error(self, tmp_path):
        f = self._dna_with_segment(tmp_path, 'nope="1"')
        with pytest.raises(ValueError, match="malformed"):
            sc.load_genbank(str(f))


# ── N9: the LbCas12a stagger is the enzyme's, wherever the site sits ───────

def test_a_staggered_cut_straddling_the_origin_reports_the_same_overhang():
    import splicecraft_crispr as cr
    rng = random.Random(3)
    seen = set()
    for trial in range(40):
        seq = "".join(rng.choice("ACGT") for _ in range(120))
        r = cr._design_guides(seq, variant="lbcas12a", circular=True,
                              limit=500)
        for g in r["guides"]:
            seen.add(g.get("overhang"))
    assert seen == {5}, seen


# ── FM11h: the NAME on the LOCUS line is not the topology ──────────────────

class TestLocusTopologyIsTheField:
    L = "LOCUS       {name:<16} 2686 bp    DNA     {topo:<9} SYN 01-JAN-2026\n"

    @pytest.mark.parametrize("name,topo", [("linear", "circular"),
                                           ("Circular", "linear"),
                                           ("bp", "linear")])
    def test_both_readers_read_the_field(self, name, topo):
        text = self.L.format(name=name, topo=topo)
        assert sc._topology_from_gb_text(text) == topo
        assert sc._gb_text_is_circular(text) is (topo == "circular")


# ── AA6: a digest with nothing to digest with is refused ───────────────────

class TestDigestNeedsAnEnzyme:
    SEQ = "GAATTC" + "ACGT" * 60

    def _digest(self, enzymes):
        import splicecraft_state as st
        return st._AGENT_HANDLERS["digest"][0](
            None, {"sequence": self.SEQ, "enzymes": enzymes})

    def test_only_unknown_names_is_a_400(self):
        r = self._digest(["EcoR1"])
        assert isinstance(r, tuple) and r[1] == 400
        assert r[0]["unknown_enzymes"] == ["EcoR1"]

    def test_a_list_with_a_typo_still_digests_and_reports_it(self):
        r = self._digest(["EcoRI", "EcoR1"])
        assert isinstance(r, dict) and r["n_cuts"] >= 1
        assert r["unknown_enzymes"] == ["EcoR1"]


# ── C6: a declared selenocysteine is neither a warning nor a stop ──────────

class TestDeclaredResiduesAreNotStops:
    CDS = "ATGAAATGACTGGCC" + "TAA"          # M K (TGA=Sec) L A *

    def _rec(self, declared=True):
        seq = "GGGG" + self.CDS + "CCCC"
        q = {"label": ["selP"], "translation": ["MKULA"]}
        if declared:
            q["transl_except"] = ["(pos:11..13,aa:Sec)"]
        return _cds_record(seq, 4, 4 + len(self.CDS), q)

    def _map_feat(self, rec):
        return [f for f in sc.PlasmidMap()._parse(rec)
                if f.get("type") == "CDS"][0]

    def test_the_map_does_not_flag_a_declared_recoding(self):
        assert not self._map_feat(self._rec(True)).get("_premature_stops")
        assert self._map_feat(self._rec(False)).get("_premature_stops")

    def test_the_edit_planner_reads_the_declared_residue(self):
        rec = self._rec(True)
        feat = [f for f in rec.features if f.type == "CDS"][0]
        plan = sc._protein_edit_plan(rec, feat, 3, "C")
        assert plan["wt_aa"] == "U" and plan["removes_stop"] is False


# ── R10: an unreadable base is a no-call everywhere, not only in one place ─

class TestNoCallsAreNotEvidence:
    """The single-read trace verdict treated an `N` as a no-call; everything
    else read it as coverage, as agreement, or — compared character by
    character — as a change. Two reads with the same `N` "confirmed" a
    mutation, a read with 300 bp of `N` was ✓ verified at 100%, and an NNK
    library reference reported its designed positions as confirmed changes."""

    @staticmethod
    def _ref(n=2000, seed=21):
        rng = random.Random(seed)
        return "".join(rng.choice("ACGT") for _ in range(n))

    @staticmethod
    def _verify(ref, reads):
        import splicecraft_state as st
        return st._AGENT_HANDLERS["verify-against-reads"][0](
            None, {"reference": ref, "reads": reads, "circular": True})

    def test_the_same_n_in_two_reads_confirms_nothing(self):
        ref = self._ref()
        read = ref[:1000] + "N" + ref[1001:]
        r = self._verify(ref, [read, read])
        assert r["consensus"]["n_confirmed"] == 0
        assert not [v for v in r["consensus"]["variants"]
                    if v.get("target_pos") == 1000]

    def test_a_run_of_n_is_neither_coverage_nor_agreement(self):
        ref = self._ref()
        read = ref[:800] + "N" * 300 + ref[1100:]
        r = self._verify(ref, [read])
        assert r["reads"][0]["covered_bp"] == 1700
        assert r["consensus"]["covered_bp"] == 1700
        res = sc._pairwise_align(read, ref)
        assert res["n_no_calls"] == 300
        code, _glyph, _col = sc._alignment_quality_status(res, len(ref))
        assert code != "verified"

    def test_a_designed_nnk_site_is_not_a_change(self):
        ref = self._ref()
        lib = ref[:600] + "NNK" + ref[603:]
        read = ref[:600] + "GCG" + ref[603:]          # one clone of the library
        r = self._verify(lib, [read, read])
        assert r["consensus"]["n_confirmed"] == 0

    def test_a_real_change_is_still_confirmed(self):
        ref = self._ref()
        sub = {"A": "C", "C": "G", "G": "T", "T": "A"}
        read = ref[:700] + sub[ref[700]] + ref[701:]
        r = self._verify(ref, [read, read])
        assert r["consensus"]["n_confirmed"] == 1


# ── R13: an inversion across the origin is one event at every rotation ─────

def _flip_block(seq, start, length):
    """`seq` with the `length`-bp block starting at `start` (wrapping)
    reverse-complemented in place."""
    n = len(seq)
    rot = seq[start:] + seq[:start]
    flipped = sc._rc(rot[:length]) + rot[length:]
    return flipped[n - start:] + flipped[:n - start]


def _inverted_bases(target, query):
    res = sc._pairwise_align(query, target)
    segs = sc._alignment_inverted_segments(res["aligned_q"], res["aligned_t"])
    return {b for s in segs for b in range(s["t_start"], s["t_end"])}


class TestOriginInversionsAreRotationInvariant:
    """Halves of a flip that straddles bp 0 were joined only when BOTH failed
    on their own. With unequal halves the bigger one passed alone, was reported
    as its own inversion, and the smaller half was lost — an 800 bp flip at
    2900 read as (0, 700) — and a 1,200 bp one answered differently at
    different rotations."""

    N = 3000
    SLACK = 25          # bases the block edges may be off by (chance matches)

    @pytest.mark.parametrize("start,length", [(2900, 800), (2950, 400),
                                              (2950, 1200)])
    def test_the_whole_block_is_found_at_every_rotation(self, start, length):
        rng = random.Random(start + length)
        plasmid = "".join(rng.choice("ACGT") for _ in range(self.N))
        construct = _flip_block(plasmid, start, length)
        truth = {(start + k) % self.N for k in range(length)}
        for rot in range(0, self.N, 100):        # 30 rotations
            p = plasmid[rot:] + plasmid[:rot]
            c = construct[rot:] + construct[:rot]
            found = {(b + rot) % self.N for b in _inverted_bases(p, c)}
            missing, extra = truth - found, found - truth
            assert len(missing) <= self.SLACK * 2, (rot, len(missing))
            assert len(extra) <= self.SLACK * 2, (rot, len(extra))


# ── CL8: a junction site is one the insert did not bring ───────────────────

def _brute_sites(seq, pat):
    ext = seq + seq[:len(pat) - 1]
    rc = sc._rc(pat)
    return [i for i in range(len(seq)) if ext[i:i + len(pat)] in (pat, rc)]


class TestJunctionSitesAreAttributedBySpan:
    """Comparing cut lists against the insert scanned on its own reported
    every Type IIS site the insert carries near an end as "formed at a
    junction" — 539 false in a fuzz against 13 real. The reference here is a
    brute-force site finder and the rule stated in the docstring: a site is a
    junction's when its recognition span is not wholly inside the insert."""

    @pytest.mark.parametrize("oh5,oh3,ptype", [
        ("GGAG", "TACT", "Promoter"), ("AATG", "GCTT", "CDS"),
        ("GCTT", "CGCT", "Terminator")])
    def test_agrees_with_brute_force(self, oh5, oh3, ptype):
        rng = random.Random(len(oh5 + oh3 + ptype))
        sites = ["GAGACC", "GGTCTC", "CGTCTC", "GAGACG"]
        for _ in range(300):
            ins = "".join(rng.choice("ACGT")
                          for _ in range(rng.randint(40, 120)))
            if ptype == "CDS":
                ins = "ATG" + ins[3:]
            if rng.random() < 0.4:
                ins = ins[:3] + rng.choice(sites) + ins[9:]
            if rng.random() < 0.3:
                ins = ins[:-6] + rng.choice(sites)
            cloned = sc._simulate_cloned_plasmid(ins, oh5, oh3, ptype)
            got = {d["enzyme"] for d in
                   sc._cloned_plasmid_regenerated_sites(cloned, ins, oh5)}
            lo = cloned.find(ins) if ins in cloned else cloned.find(ins[3:]) - 3
            want = set()
            for enz, pat in (("BsaI", "GGTCTC"), ("Esp3I", "CGTCTC")):
                for s0 in _brute_sites(cloned, pat):
                    if not (lo <= s0 and s0 + 6 <= lo + len(ins)):
                        want.add(enz)
            assert got == want, (ins, got, want)

    def test_an_insert_that_opens_with_its_own_site_is_not_a_junction(self):
        ins = "GAGACC" + "ACGT" * 15
        cloned = sc._simulate_cloned_plasmid(ins, "GGAG", "TACT", "Promoter")
        assert not sc._cloned_plasmid_regenerated_sites(cloned, ins, "GGAG")


# ── H7: an edit made while Ctrl+S takes its snapshot is not marked saved ───

async def test_an_edit_during_the_save_snapshot_stays_unsaved(monkeypatch):
    import copy as _copy
    import threading
    import time
    a, _p2 = _goto_fixture()
    real_deepcopy = _copy.deepcopy
    edited = threading.Event()

    def slow_copy(obj, *args, **kwargs):
        from Bio.SeqRecord import SeqRecord
        if (isinstance(obj, SeqRecord)
                and threading.current_thread() is not threading.main_thread()):
            edited.wait(3.0)          # the edit lands while we copy
            time.sleep(0.05)
        return real_deepcopy(obj, *args, **kwargs)

    monkeypatch.setattr(sc, "deepcopy", slow_copy)
    sc.PlasmidApp._skip_seed = True
    sc.PlasmidApp._skip_snapshot = True
    app = sc.PlasmidApp()
    async with app.run_test(size=(200, 50)) as pilot:
        for _ in range(6):
            await pilot.pause()
        while len(app.screen_stack) > 1:
            app.pop_screen()
            await pilot.pause()
        app._AUTOSAVE_DEBOUNCE_S = 10_000.0
        app._apply_record(a)
        await pilot.pause()
        app._edit_dialog_result(("GGG", "insert", 100, 100))
        await pilot.pause()
        app.action_save()
        await pilot.pause(0.1)
        app._edit_dialog_result(("TTTT", "insert", 20, 20))   # during the copy
        edited.set()
        await app.workers.wait_for_complete()
        for _ in range(6):
            await pilot.pause()
        assert len(app._current_record.seq) == 407
        assert app._unsaved, "an edit the save never wrote was marked saved"
