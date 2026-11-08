"""Plant splice-site scoring (`splicecraft_splice`).

Two classes of assertion here, and the distinction matters:

* **Biology** — the trained matrices must recover splice-site consensus that
  was never given to them explicitly. If the donor matrix does not rediscover
  `CAG` on the exon side and `AAGT` at intron +3..+6, the training data or the
  window offsets are wrong and every score is meaningless.
* **Calibration** — held-out discrimination and the flagging burden. A
  predictor that fires on every GT is a regex with extra steps.

Regression guard for the 2026-07-26 pseudocount fix: with a flat Laplace +1 and
class sizes differing ~18x, an IMPOSSIBLE base at the GT/AG anchor scored
+4.2 bits — higher than the correct one.
"""
import pytest

import splicecraft_splice as sp


def _trained():
    return bool(sp._SPLICE_DONOR_PWM and sp._SPLICE_ACCEPTOR_PWM)


pytestmark = pytest.mark.skipif(not _trained(),
                                reason="no vendored splice model")

# A canonical strong donor and acceptor, written out in window coordinates.
STRONG_DONOR = "CAGGTAAGT"                  # exon CAG | intron GTAAGT
STRONG_ACCEPTOR = "T" * 18 + "AGG"          # pyrimidine tract | AG | exon G


class TestModelRecoversBiology:
    def test_donor_consensus_matches_known_site(self):
        cons = sp._splice_consensus(sp._SPLICE_DONOR_PWM)
        assert cons[:3] == "CAG"             # exon -3..-1
        assert cons[5:9] == "AAGT"           # intron +3..+6 of GTAAGT

    def test_acceptor_has_polypyrimidine_tract(self):
        cons = sp._splice_consensus(sp._SPLICE_ACCEPTOR_PWM)
        tract = cons[:16]
        pyrimidine = sum(1 for c in tract if c in "CT")
        assert pyrimidine / len(tract) >= 0.75, tract

    def test_acceptor_prefers_G_at_exon_start(self):
        # The base immediately 3' of the AG is G-enriched in real acceptors.
        assert sp._splice_consensus(sp._SPLICE_ACCEPTOR_PWM)[20] == "G"

    def test_anchor_positions_carry_no_information(self):
        # Regression guard for 2026-07-26: the decoy population is itself
        # GT/AG-anchored, so those positions are constant in BOTH classes and
        # must contribute exactly 0 bits. A flat pseudocount made them score
        # +4.2 bits for an impossible base, because the two classes differ ~18x
        # in size and were therefore smoothed by different amounts.
        d_bits = sp._splice_position_bits(sp._SPLICE_DONOR_PWM)
        a_bits = sp._splice_position_bits(sp._SPLICE_ACCEPTOR_PWM)
        assert d_bits[3] == pytest.approx(0.0, abs=1e-6)
        assert d_bits[4] == pytest.approx(0.0, abs=1e-6)
        assert a_bits[18] == pytest.approx(0.0, abs=1e-6)
        assert a_bits[19] == pytest.approx(0.0, abs=1e-6)

    def test_informative_positions_actually_carry_signal(self):
        # The anchors being 0 must not mean the whole matrix is flat.
        d_bits = sp._splice_position_bits(sp._SPLICE_DONOR_PWM)
        assert max(d_bits) > 1.5


class TestAnchorPrecondition:
    """The GT/AG is necessary; the model only ranks among GTs."""

    def test_non_GT_is_not_scored_as_a_donor(self):
        assert sp._splice_score_donor("CAGATAAGT", 3) is None

    def test_non_AG_is_not_scored_as_an_acceptor(self):
        assert sp._splice_score_acceptor("T" * 18 + "CGG", 18) is None

    def test_real_donor_scores(self):
        # Assert not-None separately: `None > 0` is a TypeError, not a
        # failure message that tells you the site went unscored.
        score = sp._splice_score_donor(STRONG_DONOR, 3)
        assert score is not None and score > 0

    def test_out_of_range_returns_none_not_zero(self):
        # A 0.0 would read as "perfectly neutral site", not "no answer".
        assert sp._splice_score_donor("GT", 0) is None
        assert sp._splice_score_window("ACGT", sp._SPLICE_DONOR_PWM) is None

    def test_ambiguity_code_returns_none(self):
        assert sp._splice_score_window("CAGGTAANT", sp._SPLICE_DONOR_PWM) is None


class TestCalibration:
    def test_strong_site_clears_threshold(self):
        thr = sp._splice_threshold("donor")
        assert sp._splice_score_donor(STRONG_DONOR, 3) >= thr

    def test_poor_site_does_not_clear_threshold(self):
        # GT present, everything else wrong for a donor.
        poor = "TTTGTCCCC"
        assert sp._splice_score_donor(poor, 3) < sp._splice_threshold("donor")

    def test_does_not_fire_on_every_GT(self):
        # The whole point: a regex flags ~62 GT/kb. At the calibrated
        # threshold this must flag a small fraction of them.
        import random
        rng = random.Random(11)
        seq = "".join(rng.choice("ACGT") for _ in range(20_000))
        gts = sum(1 for i in range(len(seq) - 1) if seq[i:i + 2] == "GT")
        flagged = [h for h in sp._splice_scan(seq, kind="donor")]
        assert gts > 800                      # plenty of raw dinucleotides
        assert len(flagged) < gts * 0.10      # but few survive scoring

    def test_reported_threshold_is_the_one_applied(self):
        hits = sp._splice_scan(STRONG_DONOR * 40, kind="donor")
        assert hits
        assert all(h["score"] >= h["threshold"] for h in hits)


class TestScan:
    def test_finds_planted_donor_and_acceptor(self):
        seq = "ACGT" * 20 + STRONG_DONOR + "ACGT" * 30 + STRONG_ACCEPTOR + "ACGT" * 20
        hits = sp._splice_scan(seq)
        kinds = {h["kind"] for h in hits}
        assert kinds == {"donor", "acceptor"}

    def test_positions_are_forward_coordinates(self):
        seq = "ACGT" * 20 + STRONG_DONOR + "ACGT" * 20
        donor = next(h for h in sp._splice_scan(seq, kind="donor")
                     if h["score"] > 5)
        assert seq[donor["position"]:donor["position"] + 2] == "GT"

    def test_sorted_by_descending_score(self):
        seq = "ACGT" * 20 + STRONG_DONOR + "ACGT" * 30 + STRONG_ACCEPTOR + "ACGT" * 20
        scores = [h["score"] for h in sp._splice_scan(seq)]
        assert scores == sorted(scores, reverse=True)

    def test_reverse_strand_off_by_default(self):
        # mRNA is single-stranded; only the sense strand is spliced.
        seq = "ACGT" * 20 + STRONG_DONOR + "ACGT" * 20
        assert all(h["strand"] == "+" for h in sp._splice_scan(seq))

    def test_bad_kind_raises(self):
        with pytest.raises(ValueError, match="kind must be"):
            sp._splice_scan("ACGT", kind="exon")

    def test_untrained_model_raises_rather_than_reporting_clean(self):
        # An empty result from an unloaded model would read as "no sites".
        saved = sp._SPLICE_DONOR_PWM, sp._SPLICE_ACCEPTOR_PWM, sp._SPLICE_MODEL_META
        try:
            sp._splice_load_model([], [], {})
            with pytest.raises(RuntimeError, match="not trained"):
                sp._splice_scan("ACGTGTACGT")
        finally:
            sp._splice_load_model(*saved)


class TestPairRisk:
    """A lone donor is far less dangerous than a donor + acceptor pair at
    plausible intron spacing — the pair is what removes sequence."""

    def _seq(self, gap):
        return ("ACGT" * 20 + STRONG_DONOR + "ACGT" * gap
                + STRONG_ACCEPTOR + "ACGT" * 20)

    def test_pair_detected_at_plausible_intron_spacing(self):
        pairs = sp._splice_pair_risk(self._seq(30))
        assert pairs
        assert 60 <= pairs[0]["intron_len"] <= 3000

    def test_no_pair_when_too_close_to_be_an_intron(self):
        # Under 60 nt is below the plant minimum.
        assert not sp._splice_pair_risk("ACGT" * 5 + STRONG_DONOR
                                        + STRONG_ACCEPTOR + "ACGT" * 5)

    def test_in_frame_loss_flagged(self):
        # An excision that is a multiple of 3 deletes whole codons: the
        # product stays in frame and can still be ELISA-positive while being
        # biologically wrong. That is the insidious case.
        for gap in range(20, 60):
            pairs = sp._splice_pair_risk(self._seq(gap))
            if pairs and pairs[0]["intron_len"] % 3 == 0:
                assert pairs[0]["in_frame_loss"] is True
                return
        pytest.skip("no in-frame excision in the swept range")

    def test_combined_score_is_the_sum(self):
        pairs = sp._splice_pair_risk(self._seq(30))
        p = pairs[0]
        assert p["combined_score"] == pytest.approx(
            p["donor"]["score"] + p["acceptor"]["score"], abs=0.01)


class TestModelSummary:
    def test_summary_reports_provenance(self):
        s = sp._splice_model_summary()
        assert s["trained"] is True
        assert s["donor_width"] == sp._SPLICE_DONOR_LEN
        assert s["acceptor_width"] == sp._SPLICE_ACCEPTOR_LEN
        # Provenance is recorded; the clade is deliberately left generic.
        assert "verified" in s.get("source", "")

    def test_held_out_discrimination_recorded(self):
        s = sp._splice_model_summary()
        assert s["donor"]["auc"] > 0.90
        assert s["acceptor"]["auc"] > 0.90
        assert s["donor"]["decoy_fpr"] <= 0.02


# ──────────────────────────────────────────────────────────────────────────
# Field report 2026-11-08: the scan's advertised sensitivity, the pairing
# stage's pre-filter, and the annotated-end flag.
# ──────────────────────────────────────────────────────────────────────────

class TestAdvertisedSensitivityMatchesCalibration:
    """The docstring is the sentence a reader uses to decide what `0 sites`
    means, so it must not promise a sensitivity the model does not have.

    Pre-fix it claimed the cutoff was the 5th percentile of REAL site scores
    ("as strong as 95% of genuine splice sites"). The shipped calibration is
    FPR-anchored at `fpr_target` 0.01, giving 0.4593 / 0.3325 sensitivity —
    roughly 2x and 3x overstated, which is what made a sequence carrying a
    usable donor read as clean."""

    def test_docstring_does_not_claim_a_sensitivity_anchor(self):
        # The old CLAIM, not the historical note that records it was wrong:
        # the correction paragraph quotes the phrase on purpose, so match the
        # affirmative wording instead of the bare words.
        doc = " ".join((sp.__doc__ or "").split())
        assert "Thresholds are sensitivity-anchored" not in doc
        assert "a candidate scoring above it is as strong as 95%" not in doc
        assert "The default is the 5th percentile of REAL site scores" not in doc
        # and it must state the real anchoring
        assert "FALSE-POSITIVE-anchored" in doc
        assert "a clean scan is not a clean sequence" in doc.lower()
        assert "fpr_target" in doc

    def test_docstring_sensitivity_figures_match_the_shipped_meta(self):
        doc = sp.__doc__ or ""
        for kind in ("donor", "acceptor"):
            sens = sp._SPLICE_MODEL_META[kind]["sensitivity"]
            assert f"{sens:.4f}" in doc, (
                f"{kind} sensitivity {sens} is not the figure the docstring "
                f"quotes — prose and calibration have drifted apart")

    def test_basis_ships_sensitivity_with_the_threshold(self):
        b = sp._splice_scan_basis()
        for kind in ("donor", "acceptor"):
            meta = sp._SPLICE_MODEL_META[kind]
            assert b[kind]["threshold"] == meta["threshold"]
            assert b[kind]["sensitivity"] == meta["sensitivity"]
        assert b["fpr_target"] == sp._SPLICE_MODEL_META["fpr_target"]
        assert "not 'no usable site'" in b["note"]

    def test_basis_withholds_sensitivity_when_the_cutoff_is_overridden(self):
        # The recorded figure was measured AT the calibrated threshold. Reusing
        # it for a different bar would be a measurement claim nobody made.
        b = sp._splice_scan_basis(donor_threshold=1.0)
        assert b["donor"]["threshold"] == 1.0
        assert b["donor"]["sensitivity"] is None
        assert b["donor"]["decoy_fpr"] is None
        assert b["donor"]["calibrated_threshold"] == \
            sp._SPLICE_MODEL_META["donor"]["threshold"]
        # the untouched kind keeps its figures
        assert b["acceptor"]["sensitivity"] == \
            sp._SPLICE_MODEL_META["acceptor"]["sensitivity"]


class TestAboveCalibratedIsStampedByTheScan:
    """A relaxed scan must not have to re-derive 'did this clear the real bar'
    from the reported score: `score` is rounded to 3 dp, so a site at 3.76599
    reads as exactly 3.766 and compares as >= the 3.766 bar it actually
    fails."""

    def test_every_hit_carries_the_stamp(self):
        hits = sp._splice_scan("ACGT" * 10 + STRONG_DONOR + "ACGT" * 10)
        assert hits
        assert all(h["above_calibrated"] is True for h in hits)

    def test_stamp_is_independent_of_the_scan_threshold(self):
        seq = "ACGT" * 10 + STRONG_DONOR + "ACGT" * 10
        relaxed = sp._splice_scan(seq, kind="donor", donor_threshold=-99.0)
        assert len(relaxed) > 1, "a relaxed scan should admit weak GTs too"
        assert any(h["above_calibrated"] is False for h in relaxed), \
            "weak GTs admitted by a relaxed bar must be marked sub-threshold"
        assert any(h["above_calibrated"] is True for h in relaxed)

    def test_stamp_uses_the_unrounded_score(self):
        cal = sp._splice_threshold("donor")
        # A window whose score rounds TO the threshold but sits below it.
        borderline = [h for h in sp._splice_scan(
            "C" * 20 + "CAAGTAAGA" + "C" * 20, kind="donor",
            donor_threshold=-99.0) if h["score"] == pytest.approx(cal, abs=5e-4)]
        if not borderline:
            pytest.skip("no borderline window in this fixture")
        for h in borderline:
            raw = sp._splice_score_donor("C" * 20 + "CAAGTAAGA" + "C" * 20,
                                         h["position"])
            assert h["above_calibrated"] is (raw >= cal), (
                "the stamp disagrees with the unrounded comparison — this is "
                "exactly the rounding trap it exists to close")


class TestPairStageScansBelowTheReportingStage:
    """A pair carries two independent pieces of evidence, so gating each end on
    a cutoff sized for a LONE site discards the shape that actually bites: a
    mediocre donor plus a strong acceptor. Measured case from the field
    report: donor 3.596 (below the 3.766 bar) + acceptor 5.135 (above 4.439),
    invisible to the default scan and only reachable by pushing
    `donor_threshold` through `**kw` by hand."""

    # donor whose score sits just BELOW the calibrated bar
    WEAK_DONOR = "CAAGTAAGA"
    # acceptor well above it
    STRONG_ACC = "TCCTCTTCTTGCTTTTGCAGGCT"

    def _seq(self, gap=120):
        filler = ("CTCTCTCTCT" * ((gap // 10) + 1))[:gap]
        return "C" * 20 + self.WEAK_DONOR + filler + self.STRONG_ACC + "C" * 20

    def test_the_weak_donor_really_is_below_its_own_cutoff(self):
        s = sp._splice_score_window(self.WEAK_DONOR, sp._SPLICE_DONOR_PWM)
        assert s < sp._splice_threshold("donor")
        a = sp._splice_score_window(self.STRONG_ACC, sp._SPLICE_ACCEPTOR_PWM)
        assert a > sp._splice_threshold("acceptor")
        assert s + a >= (sp._splice_threshold("donor")
                         + sp._splice_threshold("acceptor")), \
            "the fixture must clear the combined gate while failing one end"

    def test_mediocre_donor_plus_strong_acceptor_is_reported_by_default(self):
        pairs = sp._splice_pair_risk(self._seq(), kind="both")
        assert pairs, ("a sub-threshold donor reaching a strong acceptor at "
                       "intron spacing must be visible without the caller "
                       "hand-lowering donor_threshold")

    def test_the_relaxed_end_is_named_not_hidden(self):
        p = sp._splice_pair_risk(self._seq(), kind="both")[0]
        assert p["both_above_threshold"] is False
        assert p["sub_threshold_ends"] == ["donor"]

    def test_a_pair_below_the_combined_gate_is_not_reported(self):
        gate = (sp._splice_threshold("donor")
                + sp._splice_threshold("acceptor"))
        pairs = sp._splice_pair_risk(self._seq(), kind="both")
        assert all(p["combined_score"] >= gate for p in pairs)
        # raising the gate past the fixture's own score drops it
        high = sp._splice_pair_risk(self._seq(), kind="both",
                                    combined_threshold=99.0)
        assert high == []

    def test_strict_pairs_are_still_labelled_as_such(self):
        both = ("ACGT" * 20 + STRONG_DONOR + "ACGT" * 30
                + STRONG_ACCEPTOR + "ACGT" * 20)
        pairs = sp._splice_pair_risk(both, kind="both")
        assert pairs
        assert pairs[0]["both_above_threshold"] is True
        assert pairs[0]["sub_threshold_ends"] == []

    def test_rule_is_a_superset_of_the_old_both_ends_pass_rule(self):
        """Relaxing the floor may only ADD pairs. Every pair the old rule
        reported — both ends over their own cutoff — must still be there."""
        import random
        gate = (sp._splice_threshold("donor")
                + sp._splice_threshold("acceptor"))
        lost = 0
        for r in range(12):
            random.seed(2000 + r)
            s = "".join(random.choices("ACGT", k=6000))
            strict = sp._splice_pair_risk(
                s, kind="both",
                donor_threshold=sp._splice_threshold("donor"),
                acceptor_threshold=sp._splice_threshold("acceptor"),
                combined_threshold=gate)
            got = {(p["donor"]["position"], p["acceptor"]["position"],
                    p["donor"]["strand"]) for p in
                   sp._splice_pair_risk(s, kind="both")}
            for p in strict:
                if (p["donor"]["position"], p["acceptor"]["position"],
                        p["donor"]["strand"]) not in got:
                    lost += 1
        assert lost == 0

    def test_floor_below_break_even_adds_nothing(self):
        """MEASURED saturation: nothing under 0.0 bits can reach the combined
        gate, which is why `_SPLICE_PAIR_FLOOR` is a compute bound rather than
        a judgement call. If this fails, the floor was picked wrong."""
        import random
        for r in range(4):
            random.seed(700 + r)
            s = "".join(random.choices("ACGT", k=8000))
            at_floor = sp._splice_pair_risk(s, kind="both")
            lower = sp._splice_pair_risk(s, kind="both", pair_floor=-3.0)
            assert len(lower) == len(at_floor)

    def test_explicit_threshold_still_wins_over_the_floor(self):
        tight = sp._splice_pair_risk(self._seq(), kind="both",
                                     donor_threshold=99.0)
        assert tight == []

    def test_uncalibrated_model_refuses_rather_than_inventing_a_gate(self):
        saved = (sp._SPLICE_DONOR_PWM, sp._SPLICE_ACCEPTOR_PWM,
                 sp._SPLICE_MODEL_META)
        try:
            sp._splice_load_model(saved[0], saved[1],
                                  {"source": "no thresholds"})
            with pytest.raises(RuntimeError, match="calibrated"):
                sp._splice_pair_risk("ACGT" * 100)
        finally:
            sp._splice_load_model(*saved)


class TestPairAnnotatedEnds:
    """A pair with one end on a REAL annotated boundary is a different animal
    from a fully cryptic one — usually the real donor reaching past its own
    partner because the model missed the real acceptor, which a
    0.33-sensitivity test does two times in three. Pre-fix the key did not
    exist and nothing distinguished them."""

    def _seq(self):
        return ("ACGT" * 20 + STRONG_DONOR + "ACGT" * 30
                + STRONG_ACCEPTOR + "ACGT" * 20)

    def test_absent_annotation_means_no_annotated_ends(self):
        for p in sp._splice_pair_risk(self._seq(), kind="both"):
            assert p["annotated_ends"] == []

    def test_donor_on_an_annotated_edge_is_flagged(self):
        pairs = sp._splice_pair_risk(self._seq(), kind="both")
        assert pairs
        dpos = pairs[0]["donor"]["position"]
        flagged = sp._splice_pair_risk(self._seq(), kind="both",
                                       annotated_edges={dpos})
        assert flagged[0]["annotated_ends"] == ["donor"]

    def test_both_ends_annotated_is_the_real_intron(self):
        pairs = sp._splice_pair_risk(self._seq(), kind="both")
        dpos = pairs[0]["donor"]["position"]
        apos = pairs[0]["acceptor"]["position"]
        flagged = sp._splice_pair_risk(self._seq(), kind="both",
                                       annotated_edges={dpos, apos})
        assert flagged[0]["annotated_ends"] == ["donor", "acceptor"]

    def test_edge_matches_the_far_base_of_the_dinucleotide_too(self):
        pairs = sp._splice_pair_risk(self._seq(), kind="both")
        dpos = pairs[0]["donor"]["position"]
        flagged = sp._splice_pair_risk(self._seq(), kind="both",
                                       annotated_edges={dpos + 2})
        assert flagged[0]["annotated_ends"] == ["donor"]


class TestPairWindowingIsExhaustiveEquivalent:
    """The pairing loop walks a bisected SLICE of the sorted acceptors rather
    than all of them. That is a pure optimization on a correctness-critical
    path, so it is asserted against the exhaustive walk it replaced — not
    merely timed. Without it the loop is O(donors x acceptors) over the whole
    molecule: 1 Mbp of AT-rich sequence put 8,352 donors against 6,346
    acceptors (53M combinations, 6.6 s), and relaxing the pair floor is what
    multiplied both counts ~10x."""

    @staticmethod
    def _exhaustive(seq, **kw):
        kw.setdefault("donor_threshold", sp._SPLICE_PAIR_FLOOR)
        kw.setdefault("acceptor_threshold", sp._SPLICE_PAIR_FLOOR)
        gate = (sp._splice_threshold("donor")
                + sp._splice_threshold("acceptor"))
        hits = sp._splice_scan(seq, **kw)
        out = []
        for strand in ("+", "-"):
            ds = sorted((h for h in hits if h["kind"] == "donor"
                         and h.get("strand", "+") == strand),
                        key=lambda h: h["position"])
            acc = sorted((h for h in hits if h["kind"] == "acceptor"
                          and h.get("strand", "+") == strand),
                         key=lambda h: h["position"])
            for d in ds:
                for a in acc:
                    span = (a["position"] + 2 - d["position"] if strand == "+"
                            else d["position"] + 2 - a["position"])
                    if 60 <= span <= 3000 and d["score"] + a["score"] >= gate:
                        out.append((strand, d["position"], a["position"], span))
        return sorted(out)

    @pytest.mark.parametrize("both_strands", [False, True])
    def test_matches_the_exhaustive_walk(self, both_strands):
        import random
        for r in range(24):
            random.seed(3000 + r)
            gc = (0.3, 0.45, 0.6)[r % 3]
            n = (400, 2000, 7000)[r % 3]
            seq = "".join(random.choices(
                "ACGT", weights=[(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2],
                k=n))
            kw = {"kind": "both"}
            if both_strands:
                kw["both_strands"] = True
            got = sorted((p["donor"]["strand"], p["donor"]["position"],
                          p["acceptor"]["position"], p["intron_len"])
                         for p in sp._splice_pair_risk(seq, **kw))
            assert got == self._exhaustive(seq, **kw), f"seed {3000 + r}"

    def test_a_pair_at_each_end_of_the_intron_window_survives(self):
        """Boundary check on the bisect bounds: exactly 60 and exactly 3000
        must both still pair, and 59 / 3001 must not."""
        for gap, expect in ((60, True), (3000, True), (59, False),
                            (3001, False)):
            # donor GT at offset 23; acceptor AG placed so the span is `gap`
            head = "C" * 20 + STRONG_DONOR
            apos = 23 + gap - 2            # AG start for this span
            pad = apos - len(head) - (len(STRONG_ACCEPTOR) - 3)
            if pad < 0:
                pytest.skip("geometry does not fit")
            seq = head + "C" * pad + STRONG_ACCEPTOR + "C" * 20
            spans = [p["intron_len"] for p in sp._splice_pair_risk(seq,
                                                                   kind="both")]
            assert (gap in spans) is expect, (gap, spans)


class TestPairStageRefusesUnusableNumbers:
    """A NaN cutoff compares False against every score, so every pair would be
    dropped and the answer would read as 'no excisable intron'. Silence
    standing in for a refusal is the failure mode this module keeps being
    audited for."""

    def test_nan_combined_threshold_refuses(self):
        with pytest.raises(ValueError, match="NaN"):
            sp._splice_pair_risk("ACGT" * 100,
                                 combined_threshold=float("nan"))

    def test_nan_pair_floor_refuses(self):
        with pytest.raises(ValueError, match="NaN"):
            sp._splice_pair_risk("ACGT" * 100, pair_floor=float("nan"))
