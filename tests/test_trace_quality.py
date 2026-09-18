"""Basecall-quality-aware trace verification (2026-09-17).

The alignment pipeline scored every mismatch the same wherever it landed. A
Sanger read does not work that way: both ends are ragged, and a disagreement
there is the instrument guessing. These tests pin the Phred-aware layer that
separates "the read stands behind this difference" from "this is trace noise".

`[INV-88]`'s quality columns are ALIGNMENT quality (identity / coverage /
indels). Everything here is BASECALL quality, the other axis.

The pure helpers are driven with REAL `_pairwise_align` output rather than
hand-written aligned strings wherever the coordinate mapping is what's under
test — a hand-built pair would let a query/target frame mix-up pass.
"""

from __future__ import annotations

import pytest

import splicecraft as sc


# A trace shaped like a real one: junk 5' lead-in, a long trustworthy middle,
# ragged 3' fall-off. Matches the profile measured on Biopython's `3730.ab1`
# (12 bases trimmed from the 5' end, 72 from the 3').
def _trace(n=300, lead=15, tail=40, good=55, bad=6):
    return [bad] * lead + [good] * (n - lead - tail) + [bad] * tail


def _seq(n=300, seed=7):
    import random
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _mutate(seq, pos, to=None):
    base = seq[pos]
    alt = to or ("A" if base != "A" else "G")
    return seq[:pos] + alt + seq[pos + 1:]


class TestPhredIsRecorded:

    def test_empty_is_not_recorded(self):
        assert sc._phred_is_recorded([]) is False
        assert sc._phred_is_recorded(None) is False

    def test_all_zero_is_not_recorded(self):
        """Regression guard for 2026-09-17: Biopython's real `Abi/310.ab1` has
        868 bases all at Phred 0 — an instrument that recorded no quality, not
        a catastrophic read. Calling it 'unusable' would tell the user their
        DNA is bad when the file simply cannot answer."""
        assert sc._phred_is_recorded([0] * 868) is False

    def test_one_real_score_is_recorded(self):
        assert sc._phred_is_recorded([0, 0, 0, 12, 0]) is True


class TestPhredTrimBounds:

    def test_trims_both_ragged_ends(self):
        ph = _trace()
        start, end = sc._phred_trim_bounds(ph, min_phred=20)
        assert 0 < start <= 15
        assert 260 <= end < 300
        # everything inside the window really is good
        inside = ph[start:end]
        assert sum(inside) / len(inside) >= 20

    def test_all_high_quality_keeps_everything(self):
        start, end = sc._phred_trim_bounds([40] * 100, min_phred=20)
        assert (start, end) == (0, 100)

    def test_all_low_quality_returns_an_empty_window(self):
        """The honest answer for a failed trace. Reporting 100% identity over
        an empty window would be the worst output this feature could give."""
        assert sc._phred_trim_bounds([5] * 100, min_phred=20) == (0, 0)

    def test_empty_input(self):
        assert sc._phred_trim_bounds([], min_phred=20) == (0, 0)
        assert sc._phred_trim_bounds(None, min_phred=20) == (0, 0)

    def test_read_shorter_than_the_window_is_still_judged(self):
        """A 6-base read must not be discarded for being smaller than the
        10-base sliding window."""
        assert sc._phred_trim_bounds([40] * 6, min_phred=20) == (0, 6)
        assert sc._phred_trim_bounds([4] * 6, min_phred=20) == (0, 0)

    def test_a_single_bad_base_does_not_end_the_read(self):
        """Window MEAN, not per-base: one dropout inside a clean stretch is
        not the end of the usable read."""
        ph = [40] * 50 + [2] + [40] * 50
        start, end = sc._phred_trim_bounds(ph, min_phred=20)
        assert (start, end) == (0, 101)

    def test_threshold_zero_keeps_the_whole_read(self):
        assert sc._phred_trim_bounds([0, 1, 2, 3] * 10, min_phred=0) == (0, 40)

    def test_non_numeric_entries_are_ignored(self):
        assert sc._phred_trim_bounds([40] * 20 + [None, "x"],  # type: ignore[list-item]
                                     min_phred=20) == (0, 20)


class TestAlignmentPositionMaps:

    def test_gap_columns_carry_no_query_base(self):
        cols = sc._alignment_query_positions("AC-GT", "ACGGT")
        assert cols == [0, 1, None, 2, 3]

    def test_target_columns_skip_target_gaps(self):
        t2c = sc._alignment_target_columns("ACGGT", "AC-GT")
        assert t2c == {0: 0, 1: 1, 2: 3, 3: 4}

    @pytest.mark.parametrize("aq,at", [
        ("", ""), ("ACGT", ""), ("", "ACGT"), ("ACGT", "ACG"),
    ])
    def test_degenerate_pairs_return_empty(self, aq, at):
        """Same refusal as `_extract_variants_from_alignment`, so the two
        walks can never disagree about what counts as usable."""
        assert sc._alignment_query_positions(aq, at) == []
        assert sc._alignment_target_columns(aq, at) == {}


class TestVariantPhred:

    def test_snp_takes_the_read_base_score(self):
        seq = _seq()
        ph = _trace()
        pos = 150                                    # inside the good middle
        ref = _mutate(seq, pos)
        r = sc._pairwise_align(seq, ref, mode="global")
        vs = sc._extract_variants_from_alignment(r["aligned_q"], r["aligned_t"])
        snps = [v for v in vs if v["type"] == "snp"]
        assert len(snps) == 1
        got = sc._variant_phred(snps[0], r["aligned_q"], r["aligned_t"], ph)
        assert got == ph[pos]

    def test_deletion_uses_the_flanking_bases(self):
        """The read has no base inside a deletion, so its confidence has to
        come from the bases either side — a deletion called between two
        Phred-6 bases is noise."""
        seq = _seq(200)
        ref = seq[:100] + "TTTT" + seq[100:]        # read is missing 4 bp
        ph = [50] * len(seq)
        r = sc._pairwise_align(seq, ref, mode="global")
        vs = sc._extract_variants_from_alignment(r["aligned_q"], r["aligned_t"])
        dels = [v for v in vs if v["type"] == "deletion"]
        assert dels, "expected a deletion"
        assert sc._variant_phred(dels[0], r["aligned_q"],
                                 r["aligned_t"], ph) == 50

    def test_insertion_takes_the_worst_inserted_base(self):
        seq = _seq(200)
        ins_at = 100
        read = seq[:ins_at] + "GGGG" + seq[ins_at:]
        ph = [50] * len(read)
        ph[ins_at + 2] = 3                          # one shaky inserted base
        r = sc._pairwise_align(read, seq, mode="global")
        vs = sc._extract_variants_from_alignment(r["aligned_q"], r["aligned_t"])
        ins = [v for v in vs if v["type"] == "insertion"]
        assert ins, "expected an insertion"
        assert sc._variant_phred(ins[0], r["aligned_q"],
                                 r["aligned_t"], ph) == 3

    def test_no_quality_array_returns_none_not_zero(self):
        """`None` and 0 must stay distinguishable: a `.gbk` consensus with no
        quality must not masquerade as a failed trace."""
        seq = _seq(120)
        ref = _mutate(seq, 60)
        r = sc._pairwise_align(seq, ref, mode="global")
        v = [x for x in sc._extract_variants_from_alignment(
            r["aligned_q"], r["aligned_t"]) if x["type"] == "snp"][0]
        assert sc._variant_phred(v, r["aligned_q"], r["aligned_t"], []) is None


class TestAnnotateVariants:

    def test_confident_is_tri_state(self):
        seq = _seq(200)
        ref = _mutate(seq, 100)
        r = sc._pairwise_align(seq, ref, mode="global")
        vs = sc._extract_variants_from_alignment(r["aligned_q"], r["aligned_t"])
        hi = sc._annotate_variants_with_quality(
            vs, r["aligned_q"], r["aligned_t"], [50] * 200, min_phred=20)
        lo = sc._annotate_variants_with_quality(
            vs, r["aligned_q"], r["aligned_t"], [8] * 200, min_phred=20)
        unk = sc._annotate_variants_with_quality(
            vs, r["aligned_q"], r["aligned_t"], [], min_phred=20)
        assert [v["confident"] for v in hi] == [True]
        assert [v["confident"] for v in lo] == [False]
        assert [v["confident"] for v in unk] == [None]

    def test_does_not_mutate_the_caller_list(self):
        """Stored alignment variants are re-read from the library; annotating
        in place would write a threshold-dependent verdict into user data."""
        seq = _seq(120)
        ref = _mutate(seq, 60)
        r = sc._pairwise_align(seq, ref, mode="global")
        vs = sc._extract_variants_from_alignment(r["aligned_q"], r["aligned_t"])
        before = [dict(v) for v in vs]
        sc._annotate_variants_with_quality(
            vs, r["aligned_q"], r["aligned_t"], [50] * 120)
        assert vs == before
        assert all("phred" not in v for v in vs)

    def test_truncated_sentinel_passes_through_unjudged(self):
        out = sc._annotate_variants_with_quality(
            [{"type": "truncated", "target_pos": 9, "length": 0,
              "ref": "", "alt": "", "omitted_after_pos": 9}],
            "ACGT", "ACGT", [40, 40, 40, 40])
        assert out[0]["type"] == "truncated"
        assert "confident" not in out[0]


class TestTraceVerificationSummary:
    """The headline behaviour: the same single mismatch reads as a real change
    or as noise depending only on the basecall quality under it."""

    def _summary(self, mut_pos, ph=None, **kw):
        seq = _seq()
        ref = _mutate(seq, mut_pos)
        r = sc._pairwise_align(seq, ref, mode="global")
        return sc._trace_verification_summary(
            r["aligned_q"], r["aligned_t"], ph if ph is not None else _trace(),
            **kw)

    def test_mismatch_in_the_good_middle_is_a_real_change(self):
        s = self._summary(150)
        assert s["verdict"] == "real_changes"
        assert (s["n_snp_confident"], s["n_snp_lowq"]) == (1, 0)

    def test_the_same_mismatch_in_the_ragged_tail_is_noise(self):
        s = self._summary(295)
        assert s["verdict"] == "low_conf"
        assert (s["n_snp_confident"], s["n_snp_lowq"]) == (0, 1)

    def test_no_quality_array_says_so_rather_than_guessing(self):
        s = self._summary(150, ph=[])
        assert s["verdict"] == "no_quality"
        assert s["has_quality"] is False

    def test_all_zero_quality_is_no_quality_not_unusable(self):
        s = self._summary(150, ph=[0] * 300)
        assert s["verdict"] == "no_quality"

    def test_a_failed_trace_is_unusable(self):
        s = self._summary(150, ph=[4] * 300)
        assert s["verdict"] == "unusable"
        assert s["usable_bp"] == 0

    def test_a_perfect_read_is_clean(self):
        seq = _seq()
        r = sc._pairwise_align(seq, seq, mode="global")
        s = sc._trace_verification_summary(
            r["aligned_q"], r["aligned_t"], _trace())
        assert s["verdict"] == "clean"
        assert s["n_confident"] == 0

    def test_lowering_the_floor_promotes_tail_noise_to_a_real_change(self):
        """The floor is a user setting, so it has to actually change the call.
        At 0 every base is trusted, so the ragged-tail mismatch counts."""
        assert self._summary(295, min_phred=0)["verdict"] == "real_changes"

    def test_a_floor_above_the_whole_read_makes_it_unusable(self):
        """Not `low_conf`: when NO window clears the bar the read proves
        nothing, and that outranks describing its variant calls. Checked
        before the variant branches on purpose."""
        assert self._summary(150, min_phred=93)["verdict"] == "unusable"

    def test_a_mismatch_in_a_dip_is_low_confidence(self):
        """The discriminating case: the read is comfortably usable overall,
        but the mismatch itself sits in a local quality dip — so there IS a
        discrepancy and the read does not stand behind it."""
        ph = _trace()
        ph[150] = 12                       # local dropout at the mutated base
        s = self._summary(150, ph=ph, min_phred=20)
        assert s["verdict"] == "low_conf"
        assert (s["n_snp_confident"], s["n_snp_lowq"]) == (0, 1)
        assert s["usable_bp"] > 200        # ...and the read was fine otherwise

    def test_reports_the_usable_window_and_mean(self):
        s = self._summary(150)
        assert s["read_bp"] == 300
        assert 0 < s["usable_bp"] < 300
        assert s["trim_start"] > 0 and s["trim_end"] < 300
        assert s["mean_phred"] > 20
        assert 0 <= s["pct_at_or_above"] <= 100


class TestSettingIsWired:

    def test_threshold_is_in_the_settings_schema(self):
        """[PREFS] — an unschema'd key is dropped by `_validate_settings`, so
        the preference would silently fail to persist."""
        assert "sanger_min_phred" in sc._SETTINGS_SCHEMA
        types, default = sc._SETTINGS_SCHEMA["sanger_min_phred"]
        assert int in types
        assert default == sc._SANGER_MIN_PHRED_DEFAULT

    def test_survives_validation_round_trip(self):
        cleaned, _warnings = sc._validate_settings({"sanger_min_phred": 30})
        assert cleaned.get("sanger_min_phred") == 30

    def test_class_level_default_exists_before_compose(self):
        """The agent API and a headless run never reach `compose()`."""
        assert sc.PlasmidApp._sanger_min_phred == sc._SANGER_MIN_PHRED_DEFAULT


class TestVerdictPhrase:

    def _summary(self, **over):
        base = {"has_quality": True, "min_phred": 20, "verdict": "clean",
                "n_confident": 0, "n_lowq": 0, "usable_bp": 700}
        base.update(over)
        return base

    def test_no_quality_is_silent(self):
        assert sc._trace_verdict_phrase(self._summary(has_quality=False)) == ""
        assert sc._trace_verdict_phrase({}) == ""
        assert sc._trace_verdict_phrase(None) == ""  # type: ignore[arg-type]

    def test_leads_with_the_actionable_number(self):
        got = sc._trace_verdict_phrase(self._summary(
            verdict="real_changes", n_confident=2, n_lowq=35))
        assert got.startswith("2 real changes")
        assert "35" in got

    def test_singular_reads_correctly(self):
        assert sc._trace_verdict_phrase(self._summary(
            verdict="real_changes", n_confident=1, n_lowq=0)
        ) == "1 real change"

    def test_low_confidence_names_the_threshold(self):
        got = sc._trace_verdict_phrase(self._summary(
            verdict="low_conf", n_lowq=4))
        assert "none well-supported" in got and "20" in got

    def test_unusable_does_not_claim_cleanliness(self):
        got = sc._trace_verdict_phrase(self._summary(verdict="unusable"))
        assert "clean" not in got and "usable" in got

    def test_clean_reports_the_window(self):
        assert "700" in sc._trace_verdict_phrase(self._summary(verdict="clean"))


class TestQualitySurvivesStorage:
    """Regression guard for 2026-09-17: the verdict is computed while the trace
    record is in hand and read back from the LIBRARY long after, so both halves
    of the persistence pair have to agree about the field. The serialiser is a
    fixed-key whitelist — the first cut stored the verdict in memory and the
    Verification Report showed an em-dash for every checked read."""

    def _entry(self, quality=None):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        target = SeqRecord(Seq("ACGT" * 20), id="pT", name="pT",
                           description="t",
                           annotations={"molecule_type": "DNA",
                                        "topology": "linear"})
        e = {"name": "read", "query_label": "read", "target_label": "pT",
             "target_record": target,
             "result": {"aligned_q": "ACGT", "aligned_t": "ACGT",
                        "identity_pct": 100.0}}
        if quality is not None:
            e["quality"] = quality
        return e

    def test_serialiser_carries_the_verdict(self):
        q = {"has_quality": True, "verdict": "real_changes",
             "n_confident": 2, "n_lowq": 35, "min_phred": 20}
        out = sc._serialize_alignment_for_storage(self._entry(q))
        assert out.get("quality") == q

    def test_absent_verdict_adds_no_dead_key(self):
        """A `.gbk` consensus has no quality; every stored alignment on disk
        must not grow a null field for it."""
        out = sc._serialize_alignment_for_storage(self._entry())
        assert "quality" not in out

    def test_empty_verdict_adds_no_dead_key(self):
        assert "quality" not in sc._serialize_alignment_for_storage(self._entry({}))

    def test_serialiser_copies_rather_than_aliasing(self):
        """The stored dict must not alias the in-memory one, or a later
        threshold change would silently rewrite what is already on disk."""
        q = {"has_quality": True, "verdict": "clean", "n_confident": 0,
             "n_lowq": 0, "min_phred": 20}
        out = sc._serialize_alignment_for_storage(self._entry(q))
        q["verdict"] = "MUTATED"
        assert out["quality"]["verdict"] == "clean"


class TestTraceCheckEndToEnd:
    """Drives the real button through the real worker — the quality split is
    worthless if it never reaches a surface the user looks at."""

    @staticmethod
    def _read_and_reference(n=300):
        """A read with real ragged ends, and a reference carrying one mutation
        in its trustworthy middle and one in the junk tail."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        seq = _seq(n)
        ph = _trace(n)
        read = SeqRecord(Seq(seq), id="rd", name="rd")
        read.letter_annotations["phred_quality"] = ph
        ref = list(seq)
        mid, tail = n // 2, n - 2
        for pos in (mid, tail):
            ref[pos] = "A" if ref[pos] != "A" else "G"
        plasmid = SeqRecord(Seq("".join(ref)), id="pREF", name="pREF",
                            description="ref",
                            annotations={"molecule_type": "DNA",
                                         "topology": "linear"})
        return read, plasmid, ph, mid, tail

    async def test_check_against_canvas_splits_real_from_noise(self):
        read, plasmid, ph, mid, tail = self._read_and_reference()
        assert ph[mid] >= 20 and ph[tail] < 20, "fixture must straddle Q20"

        app = sc.PlasmidApp()
        app._preload_record = plasmid
        async with app.run_test(size=(180, 50)) as pilot:
            for _ in range(5):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                await pilot.pause()

            screen = sc.SequencingScreen(initial_tab="sanger")
            app.push_screen(screen)
            for _ in range(6):
                await pilot.pause()
            screen._sanger_record = read
            screen._sanger_path = "rd.ab1"
            screen.query_one("#btn-sanger-verify").disabled = False
            await pilot.pause()
            await pilot.click("#btn-sanger-verify")
            for _ in range(200):
                if getattr(app, "_alignments", None):
                    break
                await pilot.pause(0.05)
            for _ in range(4):
                await pilot.pause()

            alns = getattr(app, "_alignments", []) or []
            assert len(alns) == 1, "the read should have registered once"
            q = alns[0].get("quality")
            assert isinstance(q, dict), "no quality verdict on the alignment"
            # Two planted mismatches, one either side of the threshold.
            assert q["n_snp"] == 2
            assert q["n_snp_confident"] == 1
            assert q["n_snp_lowq"] == 1
            assert q["verdict"] == "real_changes"
            assert "1 real change" in sc._trace_verdict_phrase(q)

    async def test_refuses_politely_with_an_empty_canvas(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(180, 50)) as pilot:
            for _ in range(5):
                await pilot.pause()
            while len(app.screen_stack) > 1:
                app.pop_screen()
                await pilot.pause()
            app._current_record = None
            read, _plasmid, _ph, _m, _t = self._read_and_reference()
            screen = sc.SequencingScreen(initial_tab="sanger")
            app.push_screen(screen)
            for _ in range(6):
                await pilot.pause()
            screen._sanger_record = read
            screen.query_one("#btn-sanger-verify").disabled = False
            await pilot.pause()
            await pilot.click("#btn-sanger-verify")
            for _ in range(4):
                await pilot.pause()
            assert not getattr(app, "_alignments", [])
            msg = str(screen.query_one("#sanger-preview").content)
            assert "canvas" in msg.lower()


class TestAgentParity:
    """The trace-quality split must be reachable through the side-door, not
    only from the Sanger tab — an agent that pulled reads down cannot open a
    GUI to judge them."""

    class _App:
        _current_record = None

    def _call(self, name, payload):
        import splicecraft_state as st
        fn, _w = st._AGENT_HANDLERS[name]
        return fn(self._App(), payload)

    def _ref_and_read(self, n=300):
        seq = _seq(n)
        ph = _trace(n)
        mid, tail = n // 2, n - 2
        ref = list(seq)
        for p in (mid, tail):
            ref[p] = "A" if ref[p] != "A" else "G"
        return "".join(ref), seq, ph, mid, tail

    def test_quality_splits_real_from_noise_through_the_api(self):
        ref, read, ph, mid, tail = self._ref_and_read()
        assert ph[mid] >= 20 > ph[tail], "fixture must straddle Q20"
        out = self._call("verify-against-reads",
                         {"reference": ref, "reads": [read],
                          "read_quality": [ph]})
        assert isinstance(out, dict), out
        q = (out.get("read_quality") or [None])[0]
        assert isinstance(q, dict), "no per-read quality verdict"
        assert q["n_snp"] == 2
        assert q["n_snp_confident"] == 1
        assert q["n_snp_lowq"] == 1

    def test_omitting_quality_makes_no_claim(self):
        ref, read, _ph, _m, _t = self._ref_and_read()
        out = self._call("verify-against-reads",
                         {"reference": ref, "reads": [read]})
        assert out["read_quality"] == [None]

    def test_null_entry_is_allowed_per_read(self):
        ref, read, ph, _m, _t = self._ref_and_read()
        out = self._call("verify-against-reads",
                         {"reference": ref, "reads": [read, read],
                          "read_quality": [ph, None]})
        assert isinstance(out["read_quality"][0], dict)
        assert out["read_quality"][1] is None

    @pytest.mark.parametrize("rq,needle", [
        ("nope", "must be a list"),
        ([[30], [30]], "parallel"),
        ([["hot"]], "numbers or null"),
    ])
    def test_malformed_quality_is_400(self, rq, needle):
        ref, read, _ph, _m, _t = self._ref_and_read()
        body, code = self._call("verify-against-reads",
                                {"reference": ref, "reads": [read],
                                 "read_quality": rq})
        assert code == 400 and needle in body["error"]

    def test_min_phred_is_validated_and_effective(self):
        ref, read, ph, _m, tail = self._ref_and_read()
        body, code = self._call("verify-against-reads",
                                {"reference": ref, "reads": [read],
                                 "read_quality": [ph], "min_phred": 99})
        assert code == 400 and "0-93" in body["error"]
        out = self._call("verify-against-reads",
                         {"reference": ref, "reads": [read],
                          "read_quality": [ph], "min_phred": 0})
        # At 0 every base is trusted, so the tail mismatch counts as real.
        assert (out["read_quality"][0])["n_snp_confident"] == 2
