"""Multi-read consensus — combining several reads on one plasmid.

One read agreeing with itself is not evidence. These tests pin the question the
rollup answers: for each difference, how many reads that COVER that base show
it, and how many looked and disagreed.

Alignment inputs are built with the REAL `_pairwise_align` and the REAL segment
builder rather than hand-written dicts, so a coordinate or segment-state change
breaks these instead of sliding past a fixture that encodes the old shape.
"""

from __future__ import annotations

import random

import pytest

import splicecraft as sc


def _seq(n=600, seed=3):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def _mutate(seq, pos, to=None):
    base = seq[pos]
    return seq[:pos] + (to or ("A" if base != "A" else "G")) + seq[pos + 1:]


def _read(plasmid: str, read_seq: str, label: str, *,
          quality: "dict | None" = None) -> dict:
    """An alignment entry shaped exactly like `_register_alignment` builds one,
    using the real aligner and the real segment builder."""
    result = sc._pairwise_align(read_seq, plasmid, mode="global")
    aq, at = result["aligned_q"], result["aligned_t"]
    entry = {
        "name": label, "query_label": label, "target_label": "pTest",
        "result": result, "aligned_q": aq, "aligned_t": at,
        "axis": "target",
        "segments": sc._alignment_to_target_segments(aq, at),
    }
    if quality is not None:
        entry["quality"] = quality
    return entry


def _partial(seq, lo, hi):
    """A read covering only [lo, hi) — what a real Sanger read does."""
    return seq[lo:hi]


class TestCoverage:

    def test_no_reads(self):
        s = sc._multi_read_summary([], 600)
        assert s["verdict"] == "no_reads"
        assert s["uncovered_spans"] == [(0, 600)]

    def test_one_full_read_covers_everything(self):
        p = _seq()
        s = sc._multi_read_summary([_read(p, p, "r1")], len(p))
        assert s["covered_pct"] == 100.0
        assert s["depth2_pct"] == 0.0
        assert s["uncovered_spans"] == []
        assert s["verdict"] == "clean"

    def test_partial_reads_report_the_hole(self):
        """A plasmid can be '100% identity' over the half that was read.
        Coverage is the other half of the claim."""
        p = _seq()
        reads = [_read(p, _partial(p, 0, 200), "fwd"),
                 _read(p, _partial(p, 400, 600), "rev")]
        s = sc._multi_read_summary(reads, len(p))
        assert s["covered_pct"] < 100.0
        assert s["uncovered_spans"], "the unread middle must be reported"
        lo, hi = s["uncovered_spans"][0]
        assert 180 <= lo <= 220 and 380 <= hi <= 420

    def test_overlapping_reads_report_depth_two(self):
        p = _seq()
        reads = [_read(p, _partial(p, 0, 400), "fwd"),
                 _read(p, _partial(p, 200, 600), "rev")]
        s = sc._multi_read_summary(reads, len(p))
        assert s["depth2_bp"] > 0
        assert s["depth2_bp"] <= s["covered_bp"]

    def test_a_gap_segment_is_not_coverage(self):
        """Counting an unsequenced stretch as coverage would let a hole look
        like confirmed agreement — the most dangerous error available here."""
        p = _seq()
        r = _read(p, _partial(p, 0, 200), "fwd")
        spans = sc._read_covered_spans(r)
        assert all(hi <= 260 for _lo, hi in spans), spans


class TestVariantSupport:

    def test_two_reads_agreeing_is_confirmed(self):
        p = _seq()
        mutant = _mutate(p, 300)
        reads = [_read(p, mutant, "fwd"), _read(p, mutant, "rev")]
        s = sc._multi_read_summary(reads, len(p))
        assert s["verdict"] == "confirmed"
        assert s["n_confirmed"] == 1
        v = s["variants"][0]
        assert v["target_pos"] == 300
        assert v["support"] == 2 and v["depth"] == 2
        assert v["contradicted"] == 0
        assert sorted(v["reads"]) == ["fwd", "rev"]

    def test_one_read_contradicted_by_another_is_single_read(self):
        """The artefact case: one read shows it, another covers the same base
        and disagrees."""
        p = _seq()
        reads = [_read(p, _mutate(p, 300), "fwd"), _read(p, p, "rev")]
        s = sc._multi_read_summary(reads, len(p))
        assert s["verdict"] == "single_read"
        v = s["variants"][0]
        assert v["support"] == 1 and v["depth"] == 2
        assert v["contradicted"] == 1

    def test_uncorroborated_variant_says_sequence_again(self):
        """Neither confirmed nor refuted — the honest middle answer, not a
        verdict the data cannot support."""
        p = _seq()
        reads = [_read(p, _partial(_mutate(p, 100), 0, 200), "fwd"),
                 _read(p, _partial(p, 400, 600), "rev")]
        s = sc._multi_read_summary(reads, len(p))
        assert s["verdict"] == "unconfirmed"
        v = s["variants"][0]
        assert v["support"] == 1 and v["contradicted"] == 0

    def test_different_alts_at_one_base_stay_separate(self):
        """A C>T and a C>G at the same position are two calls, not one variant
        with support 2."""
        p = _seq()
        base = p[300]
        alt1 = "A" if base != "A" else "G"
        alt2 = "T" if base != "T" else "C"
        reads = [_read(p, _mutate(p, 300, alt1), "r1"),
                 _read(p, _mutate(p, 300, alt2), "r2")]
        s = sc._multi_read_summary(reads, len(p))
        at300 = [v for v in s["variants"] if v["target_pos"] == 300]
        assert len(at300) == 2
        assert all(v["support"] == 1 for v in at300)
        assert s["verdict"] != "confirmed"

    def test_min_support_is_honoured(self):
        p = _seq()
        m = _mutate(p, 300)
        reads = [_read(p, m, "r1"), _read(p, m, "r2")]
        assert sc._multi_read_summary(reads, len(p),
                                      min_support=3)["verdict"] != "confirmed"
        assert sc._multi_read_summary(reads, len(p),
                                      min_support=2)["verdict"] == "confirmed"

    def test_variants_sort_most_supported_first(self):
        p = _seq()
        m2 = _mutate(_mutate(p, 100), 300)
        reads = [_read(p, m2, "r1"), _read(p, m2, "r2"),
                 _read(p, _mutate(p, 500), "r3")]
        s = sc._multi_read_summary(reads, len(p))
        supports = [v["support"] for v in s["variants"]]
        assert supports == sorted(supports, reverse=True)

    def test_basecall_confidence_composes(self):
        """The trace check stores WHICH variants a read backs; the rollup reads
        those positions. The per-base Phred array cannot be persisted, so this
        is the only way the two features can compose."""
        p = _seq()
        m = _mutate(p, 300)
        q = {"has_quality": True, "verdict": "real_changes",
             "n_confident": 1, "n_lowq": 0, "min_phred": 20,
             "confident_positions": [300], "lowq_positions": []}
        reads = [_read(p, m, "fwd", quality=q), _read(p, m, "rev", quality=q)]
        s = sc._multi_read_summary(reads, len(p))
        assert s["variants"][0]["confident_reads"] == 2

    def test_read_without_quality_contributes_no_confidence(self):
        p = _seq()
        m = _mutate(p, 300)
        reads = [_read(p, m, "fwd"), _read(p, m, "rev")]
        assert s_v(sc._multi_read_summary(reads, len(p)))["confident_reads"] == 0


def s_v(summary):
    return summary["variants"][0]


class TestAxisHandling:

    def test_query_axis_variants_are_numbered_in_the_query_frame(self):
        """`_extract_variants_from_alignment` numbers in its SECOND argument's
        frame, so an axis='query' alignment must pass the query string second —
        the trap `_collect_rows` documents. Getting it backwards puts every
        variant of a diff alignment at the wrong bp."""
        p = _seq(400)
        m = _mutate(p, 200)
        a = _read(p, m, "diff")
        a["axis"] = "query"
        vs = sc._alignment_variants_in_axis(a)
        # In the query (read) frame the changed base is still at 200 here
        # because the two sequences are the same length with one substitution.
        assert any(v["target_pos"] == 200 for v in vs)

    def test_degenerate_alignment_yields_nothing(self):
        assert sc._alignment_variants_in_axis({}) == []
        assert sc._alignment_variants_in_axis(
            {"result": {"aligned_q": "", "aligned_t": ""}}) == []


class TestPhrase:

    def test_phrase_names_the_evidence(self):
        p = _seq()
        m = _mutate(p, 300)
        s = sc._multi_read_summary([_read(p, m, "a"), _read(p, m, "b")],
                                   len(p))
        txt = sc._multi_read_phrase(s)
        assert "confirmed" in txt and "2" in txt

    def test_phrase_for_no_reads(self):
        assert "no reads" in sc._multi_read_phrase(
            sc._multi_read_summary([], 100))

    def test_phrase_tolerates_junk(self):
        assert sc._multi_read_phrase(None) == ""      # type: ignore[arg-type]
        assert sc._multi_read_phrase({}) == ""


class TestCoverageIsNotTheHull:
    """Regression guard for 2026-09-17: an edit-distance aligner has
    astronomically many co-optimal placements for a short read in a long
    target, and one was measured scattering a 200-mer into 160 fragments spread
    over 580 bp. Coverage must count bases actually aligned, never the hull, or
    that read would claim 580 bp of evidence for 200 bp of data."""

    def _scattered(self):
        p = _seq(600, seed=3)
        read = _mutate(p, 100)[0:200]
        a = _read(p, read, "scattered")
        return p, a

    def test_a_scattered_alignment_does_not_inflate_coverage(self):
        p, a = self._scattered()
        spans = sc._read_covered_spans(a)
        covered = sum(hi - lo for lo, hi in spans)
        extent = sc._read_observed_extent(a)
        assert extent is not None
        hull = extent[1] - extent[0]
        # The read is 200 bases long; coverage must reflect that, not the hull.
        assert covered <= 220, f"claimed {covered} bp of coverage"
        if hull > covered:
            assert covered < hull        # the split is doing real work here

    def test_covered_spans_are_merged_and_sorted(self):
        """Depth counting must not double-count one read at a segment
        boundary."""
        _p, a = self._scattered()
        spans = sc._read_covered_spans(a)
        assert spans == sorted(spans)
        for (lo1, hi1), (lo2, _hi2) in zip(spans, spans[1:]):
            assert hi1 < lo2, "touching spans should have been merged"

    def test_a_clean_partial_read_reports_no_variants(self):
        """The case the extent model exists for: a perfect 200 bp read on a
        400 bp plasmid is 50% coverage and zero variants, not one 200 bp
        deletion."""
        p = _seq(400, seed=11)
        s = sc._multi_read_summary([_read(p, p[:200], "half")], len(p))
        assert s["variants"] == []
        assert s["verdict"] == "clean"
        assert 40.0 <= s["covered_pct"] <= 60.0

    def test_stored_alignment_without_segments_still_has_coverage(self):
        """A stored alignment carries only the aligned strings; segments are
        rebuilt. Without that the agent endpoint would report every read as
        zero-coverage — 'nothing was sequenced' rather than 'I did not look'."""
        p = _seq(400, seed=12)
        full = _read(p, p[:200], "half")
        stored = {"result": full["result"], "axis": "target",
                  "query_label": "half"}          # no 'segments' key
        assert sc._read_covered_spans(stored) == sc._read_covered_spans(full)


class TestAgentEndpoints:
    """Agent parity: the rollup is reachable through the side-door."""

    class _App:
        _current_record = None
        _alignments: list = []

    def _call(self, name, payload, *, rec=None, alignments=None):
        import splicecraft_state as st
        app = self._App()
        app._current_record = rec
        app._alignments = list(alignments or [])
        fn, _w = st._AGENT_HANDLERS[name]
        return fn(app, payload)

    def _rec(self, seq):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        return SeqRecord(Seq(seq), id="pT", name="pT",
                         annotations={"molecule_type": "DNA",
                                      "topology": "circular"})

    def test_endpoints_are_registered(self):
        import splicecraft_state as st
        assert "read-consensus" in st._AGENT_HANDLERS
        assert "verify-against-reads" in st._AGENT_HANDLERS

    def test_read_consensus_on_the_loaded_plasmid(self):
        p = _seq()
        m = _mutate(p, 300)
        out = self._call("read-consensus", {}, rec=self._rec(p),
                         alignments=[_read(p, m, "fwd"), _read(p, m, "rev")])
        assert out["verdict"] == "confirmed"
        assert out["n_confirmed"] == 1
        assert out["plasmid"] == "pT"
        assert "confirmed" in out["phrase"]

    def test_read_consensus_reports_uncovered(self):
        """Coverage is part of the answer, not a footnote."""
        p = _seq()
        out = self._call("read-consensus", {}, rec=self._rec(p),
                         alignments=[_read(p, _partial(p, 0, 200), "fwd")])
        assert out["uncovered_spans"]
        assert out["covered_pct"] < 100.0

    def test_read_consensus_with_no_reads(self):
        p = _seq()
        out = self._call("read-consensus", {}, rec=self._rec(p))
        assert out["verdict"] == "no_reads"

    def test_read_consensus_needs_a_target(self):
        body, code = self._call("read-consensus", {})
        assert code == 422

    def test_min_support_validation(self):
        body, code = self._call("read-consensus", {"min_support": 0},
                                rec=self._rec(_seq()))
        assert code == 400 and ">= 1" in body["error"]

    def test_unknown_entry_is_404(self):
        body, code = self._call("read-consensus", {"id": "nope"})
        assert code == 404

    def test_verify_against_reads_carries_the_consensus(self):
        """The rollup rides on the endpoint that already has the alignments —
        no second alignment pass."""
        p = _seq(400, seed=21)
        m = _mutate(p, 200)
        out = self._call("verify-against-reads",
                         {"reference": p, "reads": [m, m]})
        assert isinstance(out, dict), out
        cons = out.get("consensus")
        assert isinstance(cons, dict), "consensus block missing"
        assert cons["n_reads"] == 2
        assert cons["verdict"] == "confirmed"
        assert any(v["support"] == 2 for v in cons["variants"])
