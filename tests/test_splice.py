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
