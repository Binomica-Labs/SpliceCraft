"""Plant translation-initiation context scoring (`splicecraft_kozak`).

Four classes of assertion, and the distinction is the point:

* **Geometry** — the window offsets must be right, and the data proves it
  rather than a comment: both trained populations are ATG-anchored, so the
  three ATG columns MUST carry ~0 bits. A non-zero there means the window
  slipped or the smoothing is wrong, and every score below is meaningless.
* **Biology** — the matrix must recover a plant initiation consensus it was
  never told: A-rich upstream, A at -3, G at +4.
* **Calibration** — the held-out numbers recorded in the metadata are the
  claim this model makes. A silent retrain that moves them has to fail here.
* **Refusals** — the three things the model deliberately does not do. Each is
  a decision backed by a measurement (see `scripts/calibrate_kozak.py`), so
  each gets a test: no partial windows, no score without an ATG, and `None`
  rather than `0.0` for "no answer" — 0.0 is a real score meaning "as likely
  a start as an upstream ATG".

The rule this sits beside is NOT deprecated: `score` / `strength` are still
the -3 purine / +4 G consensus, kept because it is interpretable. The two can
disagree, and `TestRuleAndModelCoexist` pins that they each carry their own
basis so a reader is never left guessing which one answered.
"""
import pytest

import splicecraft_kozak as kz

pytestmark = pytest.mark.skipif(not kz._KOZAK_PWM,
                                reason="no vendored initiation model")


def _score(window: str) -> float:
    """Score a window that MUST be scoreable. `_kozak_score_window` returns
    `float | None` by design, so a test comparing its result has to say which
    it expects — pyright flagged every bare comparison here."""
    s = kz._kozak_score_window(window)
    assert s is not None, window
    return s


def _pct(window: str) -> float:
    p = kz._kozak_percentile(_score(window))
    assert p is not None, window
    return p

# The consensus context the matrix should like, and the AT-rich leader context
# it should not. Both are written in window coordinates: 12 nt, ATG, +4.
STRONG = "AAAAAAGCAACA" + "ATG" + "G"
WEAK = "TTTTTTTTTTTT" + "ATG" + "T"


class TestGeometry:
    def test_window_constants_agree_with_the_vendored_matrix(self):
        assert kz._KOZAK_WIDTH == kz._KOZAK_UP + kz._KOZAK_DOWN
        assert len(kz._KOZAK_PWM) == kz._KOZAK_WIDTH
        meta = kz._KOZAK_MODEL_META
        assert (meta["up"], meta["down"]) == (kz._KOZAK_UP, kz._KOZAK_DOWN)
        assert meta["atg_offset"] == kz._KOZAK_UP

    def test_the_three_atg_columns_carry_no_information(self):
        """THE coordinate proof. Both populations are ATG-anchored, so these
        columns cannot discriminate and must score ~0 for EVERY base — which
        only holds if the pseudocount is proportional (a flat Laplace +1 over
        unequal class sizes makes an impossible base outscore the right one,
        the `splicecraft_splice` 2026-07-26 fix) AND the window offsets are
        right. If the ATG does not land here, it lands somewhere else."""
        for off in range(3):
            col = kz._KOZAK_PWM[kz._KOZAK_UP + off]
            assert max(col.values()) - min(col.values()) < 0.01, off
            for base, lod in col.items():
                assert abs(lod) < 0.01, (off, base, lod)

    def test_plus4_is_the_most_informative_position(self):
        bits = [max(c.values()) - min(c.values()) for c in kz._KOZAK_PWM]
        assert bits.index(max(bits)) == kz._KOZAK_UP + 3


class TestModelRecoversBiology:
    def test_consensus_is_the_plant_initiation_context(self):
        """Most frequent base per position among real starts — recorded at
        training time, because the log-odds MATRIX reads GC-rich (the uATGs it
        discriminates against are even more AT-rich) and would mislead anyone
        reading it as biology."""
        cons = kz._KOZAK_MODEL_META["real_start_consensus"]
        up = kz._KOZAK_UP
        assert cons[up:up + 3] == "ATG"
        assert cons[up + 3] == "G"                       # +4
        assert cons[up - 3] == "A"                       # -3 purine
        assert cons[:6].count("A") >= 5                  # A-rich leader

    def test_matrix_reads_gc_rich_because_the_decoys_are_at_rich(self):
        at = kz._KOZAK_MODEL_META["upstream_at_fraction"]
        assert at["utr_uatgs"] > at["real_starts"]

    def test_consensus_context_outscores_an_at_rich_one(self):
        assert _score(STRONG) > _score(WEAK) + 5


class TestRefusals:
    def test_no_atg_no_score(self):
        seq = "A" * 12 + "CTG" + "G"
        assert kz._kozak_score_at(seq, 12) is None
        ctx = kz._kozak_context(seq, 12)
        assert ctx["pwm"]["score"] is None
        assert "no ATG" in ctx["pwm"]["note"]

    def test_partial_window_is_refused_not_padded(self):
        """A shorter window sums fewer log-odds terms, so it is not comparable
        to the calibrated distribution. 4.2% of real plant 5'UTR upstream ATGs
        land here, so this path is exercised, not theoretical."""
        seq = "AT" + "ATG" + "G"
        assert kz._kozak_score_at(seq, 2) is None
        ctx = kz._kozak_context(seq, 2)
        assert ctx["pwm"]["score"] is None
        assert ctx["pwm"]["window"] == ""
        assert str(kz._KOZAK_UP) in ctx["pwm"]["note"]
        # ...and the rule still answers, from the bases it does have.
        assert ctx["score"] == 1 and ctx["minus3"] == ""

    def test_runs_off_the_three_prime_end(self):
        assert kz._kozak_score_at("A" * 12 + "ATG", 12) is None

    def test_unscoreable_is_none_never_zero(self):
        """0.0 is a REAL score — 'as likely a start as an upstream ATG'. A
        caller cannot distinguish it from 'no answer', so no path may use it
        as one."""
        for seq, pos in (("AT" + "ATGG", 2), ("A" * 12 + "CTGG", 12),
                         ("A" * 12 + "ATGG", 99)):
            assert kz._kozak_score_at(seq, pos) is None
        assert kz._kozak_score_window("ACGT") is None          # wrong length
        assert kz._kozak_score_window("N" + STRONG[1:]) is None  # ambiguity

    def test_no_rate_is_offered_anywhere(self):
        refused = kz._KOZAK_MODEL_META["refused"]
        assert {"efficiency_rate", "partial_windows",
                "window_past_plus4"} <= set(refused)
        ctx = kz._kozak_context(STRONG, kz._KOZAK_UP)
        assert not {"rate", "efficiency", "yield"} & set(ctx["pwm"])


class TestCalibration:
    def test_recorded_numbers_are_the_measured_ones(self):
        """These ARE the claim. A retrain that moves them must update this
        test deliberately, the same discipline as the preset catalogue digest."""
        m = kz._KOZAK_MODEL_META
        assert m["window"] == "-12..+4"
        assert m["auc_vs_utr_uatg"] == 0.766
        assert m["rule_auc_vs_utr_uatg"] == 0.712
        assert m["delta_auc_over_rule"] == 0.054
        assert m["delta_auc_ci95"] == [0.035, 0.074]
        # ...and the gain must be real: the CI cannot straddle zero.
        assert m["delta_auc_ci95"][0] > 0
        # The advantage must survive on the decoy population that the
        # coding/UTR composition shortcut cannot explain.
        assert m["auc_vs_cds_atg"] > m["rule_auc_vs_cds_atg"]

    def test_generalises_across_clades(self):
        loso = kz._KOZAK_MODEL_META["auc_leave_one_species_out"]
        assert len(loso) >= 7
        assert min(loso.values()) > 0.70

    def test_not_learned_from_a_gene_predictor(self):
        """Gnomon-predicted starts carry the predictor's own Kozak-like prior.
        If the model only worked there it would be learning the predictor."""
        by_class = kz._KOZAK_MODEL_META["auc_by_record_class"]
        assert by_class["NM"] >= by_class["XM"]

    def test_matrix_digest_is_pinned(self):
        """One constant to update on a deliberate retrain; any accidental edit
        of GENERATED data fails here instead of silently shipping."""
        import hashlib
        blob = ";".join(f"{i}:{b}:{col[b]:+.6f}"
                        for i, col in enumerate(kz._KOZAK_PWM)
                        for b in "ACGT")
        digest = hashlib.sha256(blob.encode()).hexdigest()[:16]
        assert digest == "ff973575de7e0a5d", digest


class TestPercentileAndTiers:
    def test_percentile_is_bounded_and_monotone(self):
        q = kz._KOZAK_QUANTILES
        assert kz._kozak_percentile(q[0] - 10) == 0.0
        assert kz._kozak_percentile(q[-1] + 10) == 100.0
        assert kz._kozak_percentile(None) is None
        seen = [kz._kozak_percentile(s) for s in range(-12, 9)]
        assert seen == sorted(seen)

    def test_percentile_is_against_real_starts_not_the_decoys(self):
        """The consensus context must beat most real plant starts; an AT-rich
        one must beat almost none."""
        assert _pct(STRONG) > 80
        assert _pct(WEAK) < 5

    def test_tiers_sit_exactly_on_the_measured_cuts(self):
        tiers = kz._KOZAK_MODEL_META["tiers"]
        strong, adequate = tiers["strong"]["cut"], tiers["adequate"]["cut"]
        assert adequate < strong
        assert kz._kozak_tier(strong)[0] == "strong"
        assert kz._kozak_tier(strong - 1e-9)[0] == "adequate"
        assert kz._kozak_tier(adequate)[0] == "adequate"
        assert kz._kozak_tier(adequate - 1e-9)[0] == "weak"
        assert kz._kozak_tier(None) == (None, None)

    def test_strong_tier_reports_the_burden_it_admits(self):
        """A tier with no false-positive rate beside it is an unanchored
        label. `strong` is the median real start, and 15% of upstream ATGs
        reach it — that number ships with the verdict."""
        tier, rate = kz._kozak_tier(kz._KOZAK_MODEL_META["tiers"]["strong"]["cut"])
        assert tier == "strong" and rate is not None
        assert rate == kz._KOZAK_MODEL_META["tiers"]["strong"]["utr_uatgs_at_or_above"]
        assert 0.0 < rate < 0.5


class TestRuleAndModelCoexist:
    @pytest.mark.parametrize("seq,pos,expect", [
        ("TTTATT" + "ATG" + "G", 6, 2),      # A at -3, G at +4
        ("TTTGTT" + "ATG" + "G", 6, 2),      # G at -3 is a purine too
        ("TTTATT" + "ATG" + "T", 6, 1),
        ("TTTTTT" + "ATG" + "G", 6, 1),
        ("TTTTTT" + "ATG" + "T", 6, 0),
        ("TTTCTT" + "ATG" + "C", 6, 0),
    ])
    def test_the_rule_is_unchanged(self, seq, pos, expect):
        """Ported verbatim from `test_transcript.py`: the consensus rule's
        arithmetic must not have moved when the matrix landed beside it."""
        assert kz._kozak_context(seq, pos)["score"] == expect

    def test_context_window_is_still_the_literature_one(self):
        ctx = kz._kozak_context(STRONG, kz._KOZAK_UP)
        assert ctx["context"] == "GCAACA" + "ATG" + "G"      # -6..+4
        assert len(ctx["context"]) == 10
        assert ctx["pwm"]["window"] == STRONG                # -12..+4, scored

    def test_both_halves_name_their_own_basis(self):
        ctx = kz._kozak_context("TTTTTTTTTTTT" + "ATG" + "G", kz._KOZAK_UP)
        # The rule is satisfied once (+4 G) and calls this adequate...
        assert ctx["score"] == 1 and ctx["strength"] == "adequate"
        # ...while the matrix puts it below every real plant start.
        assert ctx["pwm"]["tier"] == "weak"
        assert ctx["pwm"]["percentile"] == 0.0
        # Neither reads as "the answer" without saying what produced it.
        assert "consensus" in ctx["rule_model"]
        assert "PWM" in ctx["model"]

    def test_truncated_still_means_the_printed_window(self):
        assert kz._kozak_context("AT" + "ATGG", 2)["truncated"] is True
        assert kz._kozak_context(STRONG, kz._KOZAK_UP)["truncated"] is False


class TestTrainAndLoad:
    def test_training_recovers_the_sign_of_an_enrichment(self):
        pwm = kz._kozak_train(["AATG"] * 100, ["TATG"] * 100, 4)
        assert pwm[0]["A"] > 0 > pwm[0]["T"]

    def test_a_position_constant_in_both_classes_scores_zero(self):
        """The proportional-pseudocount property, stated as a test: an
        invariant position cannot discriminate, so every base there — even one
        that never occurs — must score ~0."""
        pwm = kz._kozak_train(["ATG"] * 500, ["ATG"] * 9000, 3)
        for col in pwm:
            for lod in col.values():
                assert abs(lod) < 0.01

    def test_no_usable_windows_raises(self):
        with pytest.raises(ValueError):
            kz._kozak_train(["AT"], ["AT"], 4)

    def test_without_a_model_the_rule_still_answers(self):
        """A missing vendored model must not make the module useless or, worse,
        silently report 0.0. Mirrors `_splice_scan`'s untrained guard, except
        the rule here can still say something true."""
        saved = (kz._KOZAK_PWM, kz._KOZAK_QUANTILES, kz._KOZAK_MODEL_META)
        try:
            kz._kozak_load_model([], [], {})
            ctx = kz._kozak_context(STRONG, kz._KOZAK_UP)
            assert ctx["pwm"]["score"] is None
            assert ctx["pwm"]["note"] == "no initiation model loaded"
            assert ctx["pwm"]["tier"] is None
            assert ctx["score"] == 2 and ctx["strength"] == "strong"
            assert "no trained model" in ctx["model"]
            assert kz._kozak_model_summary()["trained"] is False
        finally:
            kz._kozak_load_model(*saved)
        assert kz._kozak_model_summary()["trained"] is True


class TestModelSummary:
    def test_summary_carries_provenance_and_refusals(self):
        s = kz._kozak_model_summary()
        assert s["trained"] is True
        assert s["width"] == kz._KOZAK_WIDTH
        assert "RefSeq" in s["source"] and "one transcript per gene" in s["source"]
        assert "5'UTR upstream ATGs" in s["decoy"]
        assert s["uatg_coverage"] > 0.9
        assert set(s["refused"]) >= {"efficiency_rate", "partial_windows",
                                     "window_past_plus4"}


class TestDegenerateInput:
    """What the scorer does with input it cannot have expected.

    Each case here was a real answer the module gave before it was hardened,
    and in three of the four the answer LOOKED fine: a `weak` rule verdict on
    a lower-case message, a confident `adequate` for a position the message
    does not have, and a note blaming the sequence for the model's own
    mismatched width. A crash is easier to notice than any of those; a
    plausible wrong number is the failure mode that reaches a bench."""

    def test_lowercase_message_scores_as_upper_case(self):
        """The matrix half always up-cased its window; the RULE half did not.

        So a GenBank record's own lower-case ORIGIN sequence came back `weak`
        from the rule and `strong` from the matrix — two halves of one dict
        disagreeing about CASE, which is the one thing they must never
        disagree about. (`_predict_transcript` up-cases first, so this was
        reachable only through a direct call — which is exactly what the
        agent API's `score-initiation-context` and every future caller do.)"""
        up = kz._kozak_context(STRONG, kz._KOZAK_UP)
        lo = kz._kozak_context(STRONG.lower(), kz._KOZAK_UP)
        assert up["score"] == 2 and up["strength"] == "strong"
        assert (lo["score"], lo["strength"]) == (up["score"], up["strength"])
        assert lo["pwm"]["score"] == up["pwm"]["score"]
        assert lo["context"] == up["context"] == "GCAACAATGG", (
            "the printed context is DNA in upper case, as every other "
            "sequence this app shows is")

    def test_rna_spelling_scores_as_dna(self):
        """The parameter is named `mrna`. A message spelled in RNA has to
        score the same as the same message spelled in DNA — before, the rule
        scored it 2 while the matrix refused it with 'no ATG at this
        position (AUG)'."""
        dna = kz._kozak_context(STRONG, kz._KOZAK_UP)
        rna = kz._kozak_context(STRONG.replace("T", "U"), kz._KOZAK_UP)
        assert rna["pwm"]["score"] == dna["pwm"]["score"] is not None
        assert rna["pwm"]["note"] == ""
        assert rna["score"] == dna["score"]

    @pytest.mark.parametrize("pos", [-1, -5, 40, 90])
    def test_position_off_the_message_raises(self, pos):
        """A negative index used to read `plus4` from the FAR END of the
        string via Python's own wrap and report `adequate` for a base the
        message does not have at that offset; a large one raised IndexError
        from a bare slice deep inside. One check, one exception."""
        with pytest.raises(IndexError):
            kz._kozak_context("ACGT" * 10, pos)

    def test_single_quantile_model_has_no_percentile_rather_than_crashing(self):
        """`_kozak_load_model` takes any quantile list — that escape hatch is
        for a clade-matched retrain and for tests' tiny models. One quantile
        is no distribution, and the interpolation divided by zero."""
        saved = kz._KOZAK_QUANTILES
        try:
            kz._KOZAK_QUANTILES = [1.0]
            assert kz._kozak_percentile(0.5) is None
            kz._KOZAK_QUANTILES = []
            assert kz._kozak_percentile(0.5) is None
        finally:
            kz._KOZAK_QUANTILES = saved
        assert kz._kozak_percentile(_score(STRONG)) is not None

    def test_wrong_width_model_blames_the_model_not_the_sequence(self):
        """The fall-through note said "window carries a base outside ACGT",
        which is undebuggable when the bases are plainly ACGT and it is the
        loaded matrix that is the wrong size."""
        saved = kz._KOZAK_PWM
        try:
            kz._KOZAK_PWM = [{b: 0.0 for b in "ACGT"} for _ in range(8)]
            ctx = kz._kozak_context(STRONG, kz._KOZAK_UP)
            assert ctx["pwm"]["score"] is None
            assert "8 nt wide" in ctx["pwm"]["note"]
            assert "outside ACGT" not in ctx["pwm"]["note"]
        finally:
            kz._KOZAK_PWM = saved
