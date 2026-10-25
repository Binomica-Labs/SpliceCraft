"""Regression tests for the 2026-10-24 read-verification / alignment audit.

A read-only auditor probed the whole read path against outside references —
brute-force edit distance, Biopython, the planting that built each read — and
every finding below was reproduced here before it was fixed. Each test fails
on the tree as it stood before this pass (v1.3.11).

The shape most of these share: TWO HALVES OF ONE ANSWER DISAGREEING. The same
response said `verdict: "match"` and `consensus.verdict: "confirmed"` about one
100 bp deletion; the batch report counted a perfect read as 2 indels while the
consensus called it clean; `read-trace` reported a trim the GUI does not do.
An endpoint that contradicts itself is worse than one that refuses, because
only one of the two numbers gets read.
"""
from __future__ import annotations

import math

import pytest

import splicecraft as sc


def _handler(name):
    return sc._state._AGENT_HANDLERS[name][0]


def _rnd(n, seed):
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _verify(ref, reads, **kw):
    r = _handler("verify-against-reads")(None, {"reference": ref,
                                                "reads": reads, **kw})
    assert not isinstance(r, tuple), r
    return r


def _hetero(**body):
    r = _handler("analyse-read-heterogeneity")(None, body)
    assert not isinstance(r, tuple), r
    return r


def _other(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


# ══════════════════════════════════════════════════════════════════════════
# Long partial reads: the cap was quadratic in read length
# ══════════════════════════════════════════════════════════════════════════

class TestALongPartialReadCanBeVerified:
    """`_semiglobal_align`'s window is `read + ~664` wide, so `read x window`
    passed the 60M-cell cap at about 7.4 kb — against a documented 200 kb
    limit. 16 of 16 perfect 7.5-9.5 kb reads were unverifiable: 8 refused, 7
    given the WRONG reason, and one returned `mismatch` with 1,837 phantom
    variants because the only candidate left was a chance seed in the other
    orientation. 8-10 kb nanopore reads are ordinary."""

    REF = _rnd(15_000, 11)

    @pytest.mark.parametrize("read_len", [7_600, 9_500, 12_000])
    def test_a_perfect_long_read_is_clean(self, read_len):
        read = self.REF[1_000:1_000 + read_len]
        r = _verify(self.REF, [read], circular=True)
        assert r["verdict"] == "match", r.get("warnings")
        assert r["reads"][0]["covered_identity_pct"] == 100.0
        assert r["reads"][0]["covered_bp"] == read_len
        assert r["consensus"]["variants"] == []

    def test_a_reverse_complemented_long_read_is_clean(self):
        read = sc._rc(self.REF[1_000:10_000])
        r = _verify(self.REF, [read], circular=True)
        assert r["verdict"] == "match" and r["reads"][0]["rc"] is True
        assert r["consensus"]["variants"] == []

    def test_a_long_read_across_the_origin_is_clean(self):
        read = self.REF[11_000:] + self.REF[:5_000]
        r = _verify(self.REF, [read], circular=True)
        assert r["verdict"] == "match"
        assert r["reads"][0]["covered_bp"] == len(read)
        assert r["consensus"]["variants"] == []

    def test_one_snp_in_a_long_read_is_the_only_variant(self):
        at = 5_000
        read = (self.REF[1_000:at] + _other(self.REF[at])
                + self.REF[at + 1:10_000])
        r = _verify(self.REF, [read], circular=True)
        got = [(v["type"], v["target_pos"]) for v in r["consensus"]["variants"]]
        assert got == [("snp", at)], got

    @pytest.mark.parametrize("size", [10, 500, 900])
    def test_an_indel_a_long_read_spans_is_sized_exactly(self, size):
        """The band follows the read's seed-chain FOOTPRINT, not its length: a
        read missing 900 bp covers 900 MORE reference than it has bases, and a
        band sized from the read reported the real deletion plus a spurious
        insertion where the band ran out."""
        at = 5_000
        read = self.REF[1_000:at] + self.REF[at + size:10_000 + size]
        r = _verify(self.REF, [read], circular=True)
        got = [(v["type"], v["target_pos"], v["length"])
               for v in r["consensus"]["variants"]]
        # ONE deletion, sized exactly. Its position may left-normalise a base
        # or two where the flanks repeat, which is a different event from a
        # band that clipped the deletion and invented an insertion.
        assert len(got) == 1, got
        kind, pos, length = got[0]
        assert (kind, length) == ("deletion", size), got
        assert abs(pos - at) <= 2, got

    def test_one_long_read_does_not_fail_a_whole_heterogeneity_run(self):
        """A single over-cap read made `analyse-read-heterogeneity` 400 for
        every read in the request."""
        reads = [self.REF] * 3 + [self.REF[1_000:10_000]]
        r = _hetero(reference=self.REF, reads=reads, circular=True)
        assert r["n_reads"] == 4 and r["verdict"] == "clonal"


# ══════════════════════════════════════════════════════════════════════════
# Local mode never rotated
# ══════════════════════════════════════════════════════════════════════════

class TestLocalModePlacesAReadAcrossTheOrigin:
    """Rotation was gated on `identity_pct`, which in LOCAL mode is the
    identity of the aligned BLOCK (~100% for a clean read) — and the terminal-
    tail trigger reads gap columns local mode does not produce. Both triggers
    were blind, so a partial read crossing bp 0 was soft-clipped to one side
    and a real SNP in the clipped half was missed 40 times out of 40, on
    forward reads as much as reverse-complemented ones."""

    REF = _rnd(3_000, 13)

    @pytest.mark.parametrize("rc", [False, True])
    def test_a_snp_in_the_clipped_half_is_found(self, rc):
        at = 2_900                      # before bp 0, in the clipped arc
        mutated = self.REF[:at] + _other(self.REF[at]) + self.REF[at + 1:]
        read = mutated[2_750:] + mutated[:550]
        r = _verify(self.REF, [sc._rc(read) if rc else read],
                    circular=True, mode="local")
        got = [(v["type"], v["target_pos"]) for v in r["consensus"]["variants"]]
        assert got == [("snp", at)], got

    def test_the_whole_read_is_placed_not_clipped(self):
        read = self.REF[2_750:] + self.REF[:550]
        r = _verify(self.REF, [read], circular=True, mode="local")
        assert r["verdict"] == "match"
        assert r["reads"][0]["covered_bp"] == len(read)


class TestAPartialReadIsNeverQueryRotated:
    """Rotating a PARTIAL read splices its halves together in an order no
    molecule has, and the chimera could out-score the honest placement: one
    read came back `covered_bp: 4000` at 48% identity while
    `uncovered_spans` said [(1050, 3114)] — one answer describing itself two
    ways. Target rotation is the meaningful move for a partial read."""

    def test_the_picker_does_not_rotate_a_partial_read(self):
        t = _rnd(4_000, 17)
        read = t[3_600:] + t[:400]          # 800 bp across the origin
        res = sc._pick_best_rotation(read, t, is_circular=True)
        assert res["picked_rotation"] != "query", res["picked_rotation"]
        assert res["query_rotation"] == 0
        cov, cov_bp = sc._covered_identity(res, len(t), circular=True)
        assert (cov, cov_bp) == (100.0, len(read))


# ══════════════════════════════════════════════════════════════════════════
# Heterogeneity: the false-`clonal` directions
# ══════════════════════════════════════════════════════════════════════════

class TestHeterogeneityNeverMutesASubPopulationSilently:
    """Sub-`min_phred` observations were dropped from `alt_reads` while those
    reads STAYED in the rollup's `depth`, so nothing was left to form a
    fraction from and the position vanished from the answer entirely: a true
    50/50 mixture basecalled at Q15 — ordinary for older ONT — came back
    `clonal`, `positions: []`, with no warning. That is the false-clonal
    direction this endpoint's own comments call dangerous."""

    REF = _rnd(2_000, 23)
    POS = 700

    def _reads(self, n_alt=10, n_ref=10):
        alt = _other(self.REF[self.POS])
        return ([self.REF[:self.POS] + alt + self.REF[self.POS + 1:]] * n_alt
                + [self.REF] * n_ref)

    def test_a_fifty_fifty_mixture_at_q15_is_not_clonal(self):
        reads = self._reads()
        r = _hetero(reference=self.REF, reads=reads, circular=False,
                    read_quality=[[15] * len(reads[0])] * len(reads))
        assert r["verdict"] == "inconclusive", r["verdict"]
        assert r["n_muted_low_quality"] == 1
        assert any("Q20" in w for w in r["warnings"]), r["warnings"]

    def test_the_muted_position_is_still_reported_with_its_reason(self):
        reads = self._reads()
        r = _hetero(reference=self.REF, reads=reads, circular=False,
                    read_quality=[[15] * len(reads[0])] * len(reads))
        [pos] = [p for p in r["positions"] if p["pos"] == self.POS]
        assert pos["muted_by"] == "low_quality"
        assert pos["counted_in_verdict"] is False
        # Scored on the basecalls that passed (none), beside what the reads
        # actually showed (half of them).
        assert pos["alt_reads"] == 0 and pos["low_quality_reads"] == 10
        assert pos["fraction"] == 0.0 and pos["observed_fraction"] == 0.5

    def test_a_genuinely_clonal_sample_at_q15_is_still_clonal(self):
        """The guard must not fire on low quality alone — only where the poor
        basecalls were the ONLY support for a difference."""
        reads = [self.REF] * 20
        r = _hetero(reference=self.REF, reads=reads, circular=False,
                    read_quality=[[15] * len(self.REF)] * 20)
        assert r["verdict"] == "clonal"
        assert r["n_muted_low_quality"] == 0

    def test_good_quality_still_reads_the_mixture(self):
        reads = self._reads()
        r = _hetero(reference=self.REF, reads=reads, circular=False,
                    read_quality=[[40] * len(reads[0])] * len(reads))
        assert r["verdict"] == "mixed"
        [pos] = [p for p in r["positions"] if p["pos"] == self.POS]
        assert pos["fraction"] == 0.5 and pos["muted_by"] is None


class TestHeterogeneityWithNothingCoveredIsNotClonal:
    """`clonal` means "the reads are one sequence". Three all-N reads aligned
    and then reported `clonal` at `covered_bp: 0` with no warning, which reads
    exactly like a clean result; reads covering 3.8% of the plasmid said the
    same thing as a fully covered clonal culture."""

    REF = _rnd(2_000, 29)

    def test_no_coverage_is_no_reads(self):
        r = _hetero(reference=self.REF, reads=["N" * 2_000] * 3,
                    circular=False)
        assert r["covered_bp"] == 0
        assert r["verdict"] == "no_reads"
        assert any("nothing to judge" in w for w in r["warnings"])

    def test_thin_coverage_says_so(self):
        r = _hetero(reference=self.REF, reads=[self.REF[:25]] * 3,
                    circular=False)
        assert r["verdict"] == "clonal"
        assert any("cover only" in w for w in r["warnings"]), r["warnings"]


# ══════════════════════════════════════════════════════════════════════════
# verify-against-reads: a verdict that contradicted its own response
# ══════════════════════════════════════════════════════════════════════════

class TestMatchNeverContradictsTheConsensus:
    """A 100 bp deletion in a 10 kb plasmid is 99.0% covered identity, which
    clears the default threshold — so two reads both carrying it verified as
    `verdict: "match"` while `consensus.verdict` in the SAME response read
    `confirmed` for that deletion and both trace verdicts said
    `real_changes`. An identity threshold cannot express "no indel"."""

    REF = _rnd(10_000, 31)

    def _read_with_deletion(self, at=4_000, size=100):
        return self.REF[:at] + self.REF[at + size:]

    @pytest.mark.parametrize("n_reads", [1, 2])
    def test_a_confident_deletion_blocks_match(self, n_reads):
        r = _verify(self.REF, [self._read_with_deletion()] * n_reads,
                    circular=True)
        assert r["reads"][0]["passes"] is True       # still clears identity
        assert r["verdict"] == "mismatch"
        assert any("cannot express" in w for w in r.get("warnings") or [])

    def test_a_clean_read_still_matches(self):
        r = _verify(self.REF, [self.REF], circular=True)
        assert r["verdict"] == "match" and not r["consensus"]["variants"]

    def test_a_contradicted_difference_does_not_block_match(self):
        """One read showing a deletion another read LOOKED AT and disagreed
        with is not evidence; `min_identity` is still the tolerance there."""
        r = _verify(self.REF, [self._read_with_deletion(), self.REF],
                    circular=True)
        assert r["verdict"] == "match", r.get("warnings")


class TestInvertedMeansOnlyTheOrientationIsWrong:
    """The verdict accepted any read that merely HAD an inversion, so a read
    with an inversion PLUS an unrelated 400 bp deletion was reported as
    `inverted` with covered identity down to 82.7%. The docs define the
    verdict as "the bases are right, the orientation isn't"."""

    REF = _rnd(4_000, 37)

    def _read(self, with_deletion=0):
        s, e = 900, 1_400
        read = self.REF[:s] + sc._rc(self.REF[s:e]) + self.REF[e:]
        if with_deletion:
            d = 2_600
            read = read[:d] + read[d + with_deletion:]
        return read

    def test_an_inversion_alone_is_inverted(self):
        r = _verify(self.REF, [self._read()], circular=True)
        assert r["verdict"] == "inverted"
        assert r["reads"][0]["inversions_only"] is True

    @pytest.mark.parametrize("size", [150, 400])
    def test_an_inversion_plus_an_unrelated_deletion_is_a_mismatch(self, size):
        r = _verify(self.REF, [self._read(with_deletion=size)], circular=True)
        assert r["verdict"] == "mismatch"
        read = r["reads"][0]
        assert read["inversions"], "the inversion is still reported"
        assert read["inversions_only"] is False
        # Identity outside the inverted spans is what the verdict rests on.
        assert read["identity_outside_inversions_pct"] < 99.0


class TestReadOverhangPastALinearEndIsNotAnInsertion:
    """On a LINEAR reference the read may hang off either end, and those
    columns are gaps in the target row exactly like an insertion — but there
    is no reference base there to insert into. A perfect read plus 30 bp of
    adapter scored 98.52% and came back `mismatch`, while the consensus in
    the same response was clean and the trace verdict agreed with it. The
    same overhang rule `_variant_gate_pos` already applies."""

    REF = _rnd(2_000, 41)
    FLANK = "GGATCCTTAATTAAGCGGCCGCTTAATTAAGC"[:30]

    @pytest.mark.parametrize("read,label", [
        (FLANK + REF, "lead"),
        (REF + FLANK, "trail"),
        (FLANK + REF + FLANK, "both"),
    ])
    def test_overhang_still_matches(self, read, label):
        r = _verify(self.REF, [read], circular=False)
        assert r["verdict"] == "match", (label, r["reads"][0])
        assert r["reads"][0]["covered_identity_pct"] == 100.0
        assert r["consensus"]["variants"] == []

    def test_an_insertion_past_a_circle_is_still_an_insertion(self):
        """There is no end to hang off a circle: an insertion after the last
        base is at bp 0 and counts."""
        r = _verify(self.REF, [self.REF + self.FLANK], circular=True)
        assert r["verdict"] == "mismatch"


# ══════════════════════════════════════════════════════════════════════════
# Counts that have to agree ([INV-87])
# ══════════════════════════════════════════════════════════════════════════

class TestTheTraceVerdictCountsAnOriginIndelOnce:
    """`_trace_verification_summary` extracted its own variants with no
    `circular=` and no origin-split merge, so ONE indel crossing bp 0 arrived
    as the two halves the rows split it into: the toast said "2 real changes"
    where the consensus listed one deletion. Both callers — the GUI's "Check
    against canvas" and the agent — let `variants` default."""

    REF = _rnd(6_000, 43)

    def test_a_deletion_across_bp_zero_is_one_event(self):
        size = 9
        # A full-length read missing `size` bases straddling the origin.
        read = self.REF[5:-4]
        res = sc._pick_best_rotation(read, self.REF, is_circular=True,
                                     mode="global", canvas_axis="target")
        q = sc._phred_in_alignment_frame([40] * len(read), res)
        summary = sc._trace_verification_summary(
            res["aligned_q"], res["aligned_t"], q, circular=True,
            frame_shift=sc._read_rotation_in_rows(res))
        assert summary["n_indel"] == 1, summary
        assert summary["n_indel_confident"] == 1
        assert len(summary["confident_positions"]) == 1
        cons = sc._multi_read_summary([{"result": res, "axis": "target"}],
                                      len(self.REF), circular=True)
        assert len(cons["variants"]) == summary["n_indel"]
        assert cons["variants"][0]["length"] == size


class TestTheBatchReportAgreesWithThePerReadVerdict:
    """`VerificationReportModal._collect_rows` extracted variants raw: no
    extent gate, so a perfect 800 bp Sanger read showed Indels=2 (its two
    unread arcs) and jumped the cursor to bp 0, and a 20-base N run showed
    SNPs=20 while the consensus was clean. The Alignment Manager's "Gaps"
    column promises to agree with that number and said 1."""

    REF = _rnd(5_000, 47)

    @pytest.mark.parametrize("read,label", [
        (REF[2_000:2_800], "partial"),
        (REF[4_600:] + REF[:400], "across the origin"),
        (REF[2_000:2_300] + "N" * 20 + REF[2_320:2_800], "with no-calls"),
    ])
    def test_a_clean_read_counts_zero_everywhere(self, read, label):
        res = sc._pick_best_rotation(read, self.REF, is_circular=True,
                                     mode="global", canvas_axis="target")
        align = {"result": res, "axis": "target"}
        variants = sc._alignment_reportable_variants(
            align, len(self.REF), circular=True)
        cons = sc._multi_read_summary([align], len(self.REF), circular=True)
        assert variants == [], (label, variants)
        assert sc._alignment_reportable_indels(
            align, len(self.REF), circular=True) == 0, label
        assert cons["variants"] == [], label

    def test_a_real_snp_is_still_counted(self):
        at = 2_400
        read = (self.REF[2_000:at] + _other(self.REF[at])
                + self.REF[at + 1:2_800])
        res = sc._pick_best_rotation(read, self.REF, is_circular=True,
                                     mode="global", canvas_axis="target")
        variants = sc._alignment_reportable_variants(
            {"result": res, "axis": "target"}, len(self.REF), circular=True)
        assert [(v["type"], v["target_pos"]) for v in variants] == [
            ("snp", at)]

    def test_the_stored_verdict_wins_when_there_is_one(self):
        """A row with a per-read verdict stored at align time reports THAT
        count: it saw the per-base Phred array the report cannot, and two
        numbers derived twice are two numbers that can disagree."""
        res = sc._pick_best_rotation(self.REF[2_000:2_800], self.REF,
                                     is_circular=True, mode="global",
                                     canvas_axis="target")
        align = {"result": res, "axis": "target",
                 "quality": {"n_snp": 3, "n_indel": 2}}
        assert sc._alignment_reportable_indels(
            align, len(self.REF), circular=True) == 2


class TestOneDeletionBesideItsNeighboursIsStillOneDeletion:
    """`_refine_indel_clusters` bounded each window at the midpoint to its
    neighbour's UNGROWN edge, so a neighbour that had grown 12 bp towards its
    own SNP left this window 6 of the 8 context columns it needs — with 21
    clean columns between the two events — and the cluster went unrefined: a
    3 bp deletion stayed reported as two 1 bp ones, which [INV-213] says must
    not happen."""

    def test_a_deletion_with_a_snp_and_a_deletion_nearby_is_not_split(self):
        ref = _rnd(400, 70017)
        # The spacing that isolates it: a 3 bp deletion at 106, an SNP at
        # 130, and the 2 bp deletion at 142 — 21 clean columns from the SNP.
        read = (ref[:106] + ref[109:130] + _other(ref[130])
                + ref[131:142] + ref[144:])
        aq, at = sc._myers_align_global(read, ref)
        raw = [(v["type"], v["target_pos"], v["length"])
               for v in sc._extract_variants_from_alignment(aq, at)]
        # The unit-cost engine splits it, which is what refinement is for.
        assert ("deletion", 140, 1) in raw and ("deletion", 143, 1) in raw, raw
        rq, rt = sc._refine_indel_clusters(aq, at)
        variants = [(v["type"], v["target_pos"], v["length"])
                    for v in sc._extract_variants_from_alignment(rq, rt)]
        assert ("deletion", 142, 2) in variants, variants
        assert [v for v in variants if v[0] == "deletion"] == [
            ("deletion", 106, 3), ("deletion", 142, 2)], variants


# ══════════════════════════════════════════════════════════════════════════
# read-trace / align-plasmidsaurus-zip
# ══════════════════════════════════════════════════════════════════════════

class TestTheTrimReportedIsTheTrimTheGuiDoes:
    """`trimmed_span` took the outermost base at or above the floor, so one
    lucky Q25 base held 75 Q5 tail bases inside the "usable" span — [15, 576]
    where `_phred_trim_bounds` trims to (9, 500). The docstring promises the
    GUI's trim, and two trims that disagree make it false."""

    def test_the_span_matches_phred_trim_bounds(self, tmp_path, monkeypatch):
        qual = [5] * 20 + [35] * 480 + [5] * 75
        qual[540] = 25                      # the one lucky base in the tail
        want = sc._phred_trim_bounds(qual, min_phred=20)
        # The naive rule the endpoint used kept 35 more tail bases on the
        # strength of that single base.
        naive = [i for i, q in enumerate(qual) if q >= 20]
        assert [naive[0], naive[-1] + 1] != list(want)
        assert want[1] < naive[-1] + 1

        class _Rec:
            id = "trace"
            seq = "A" * len(qual)
            letter_annotations = {"phred_quality": qual}

        p = tmp_path / "t.ab1"
        p.write_bytes(b"ABIF")
        monkeypatch.setattr(sc._state, "_DUMMY", None, raising=False)
        import splicecraft_agent as _ag
        monkeypatch.setattr(_ag, "_ab1_path_to_record", lambda _p: _Rec())
        monkeypatch.setattr(_ag, "_get_setting",
                            lambda k, d=None: 20 if k == "sanger_min_phred"
                            else d)
        r = _handler("read-trace")(None, {"path": str(p), "trim": True})
        assert not isinstance(r, tuple), r
        assert r["trimmed_span"] == list(want), (r["trimmed_span"], want)

    def test_a_failed_trace_reports_no_usable_span(self):
        """`_phred_trim_bounds` returns (0, 0) when nothing clears the bar,
        and that must read as "no usable window", not as the whole read."""
        assert sc._phred_trim_bounds([5] * 100, min_phred=20) == (0, 0)


class TestAnAlignmentNeverReportsARenameThatDidNotHappen:
    """The rename check compared the library's display name to
    `target_record.name`, the GenBank LOCUS form of it — spaces become
    underscores ([INV-98]) — so `_target_renamed_to` was reported for every
    entry whose name holds a space, and `target_name` came back in a form
    that does not resolve the entry it just named."""

    def test_the_locus_form_is_not_the_display_name(self):
        """The premise: a display name with a space does not survive into the
        LOCUS field unchanged, so comparing the two can only ever report a
        drift that did not happen. The endpoint-level assertion lives in
        `test_agent_api.py::test_align_resolves_target_by_name`."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        name = "Test plasmid v2"
        rec = SeqRecord(Seq("ACGT" * 50), id=name, name=name,
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        parsed = sc._gb_text_to_record(sc._record_to_gb_text(rec), cache=False)
        assert parsed.name != name, parsed.name
        assert " " not in parsed.name


# ══════════════════════════════════════════════════════════════════════════
# read_quality validation
# ══════════════════════════════════════════════════════════════════════════

class TestPhredArraysAreValidatedOnce:
    """Two copies of this validation disagreed. NaN and +-Infinity are floats
    and passed the `isinstance` test, then `int(v)` raised and the request
    became a scrubbed 500; a negative Phred was accepted outright; and the
    length was compared against the RAW read in one endpoint and the CLEANED
    read in the other, so a read with a trailing newline got a 400 on a
    correctly sized array."""

    REF = _rnd(500, 59)

    @pytest.mark.parametrize("endpoint,key", [
        ("verify-against-reads", "reads"),
        ("analyse-read-heterogeneity", "reads"),
    ])
    @pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf, -5, 94])
    def test_an_unusable_phred_is_a_clean_400(self, endpoint, key, bad):
        q = [30] * len(self.REF)
        q[10] = bad
        r = _handler(endpoint)(None, {"reference": self.REF,
                                      key: [self.REF],
                                      "read_quality": [q],
                                      "circular": False})
        assert isinstance(r, tuple) and r[1] == 400, r
        assert "read_quality[0]" in r[0]["error"], r[0]

    def test_a_read_with_whitespace_still_matches_its_array(self):
        r = _verify(self.REF, [self.REF + "\n"], circular=False,
                    read_quality=[[30] * len(self.REF)])
        assert r["verdict"] == "match"

    def test_a_short_array_is_still_refused(self):
        r = _handler("verify-against-reads")(None, {
            "reference": self.REF, "reads": [self.REF],
            "read_quality": [[30] * 10], "circular": False})
        assert isinstance(r, tuple) and r[1] == 400
        assert "one per base" in r[0]["error"]


# ── Hardening pass (edge cases, 2026-10-25) ──────────────────────────────────


class TestAStoredAlignmentWeCannotReadBreaksNothing:
    """`alignments.json` is a file on disk, and `align.get("result") or {}`
    let a TRUTHY non-dict through to `.get`. Four readers then raised
    AttributeError — the variant extraction behind the Verification Report
    and the Alignment Manager's "Gaps" column, and the two span helpers the
    overlay paints from. (`_alignment_quality_status`, which paints the
    LibraryPanel glyph, already coerced; `_alignment_result_dict` makes the
    guarantee uniform instead of one guard per reader.)"""

    JUNK = ["oops", [1, 2], 7, 3.5, True]

    @pytest.mark.parametrize("junk", JUNK)
    def test_the_variant_extraction_contributes_nothing(self, junk):
        assert sc._alignment_variants_in_axis({"result": junk},
                                              circular=True) == []

    @pytest.mark.parametrize("junk", JUNK)
    def test_the_overlay_span_helpers_contribute_nothing(self, junk):
        assert sc._read_segments({"result": junk}) == []
        assert sc._read_covered_spans({"result": junk}) == []

    @pytest.mark.parametrize("junk", JUNK)
    def test_the_report_and_the_gaps_column_contribute_nothing(self, junk):
        align = {"result": junk, "visible": True}
        assert sc._alignment_reportable_variants(align, 100,
                                                 circular=True) == []
        assert sc._alignment_reportable_indels(align, 100, circular=True) == 0

    def test_a_real_result_still_reads(self):
        res = {"aligned_q": "ACGT", "aligned_t": "ACGT"}
        assert sc._alignment_result_dict({"result": res}) is res

    def test_one_unreadable_alignment_does_not_blank_the_library_column(self):
        # The end-to-end property, driven the way LibraryPanel drives it.
        entry = {"id": "e1", "name": "x", "size": 100,
                 "alignments": [{"result": "oops", "visible": True}]}
        summary = sc._library_entry_alignment_summary(entry)
        assert summary is not None
        assert isinstance(summary.get("glyph"), str)


class TestACorruptIndelCountIsNotAnIndel:
    """`isinstance(True, int)` is True, so a corrupted `"n_indel": true` on
    disk reported one indel on a clean read. The same bool guard
    `_merge_origin_split_indels` already applies to a stored position."""

    ROWS = {"aligned_q": "ACGTACGTAC", "aligned_t": "ACGTACGTAC"}

    def test_a_boolean_indel_count_is_ignored(self):
        align = {"quality": {"n_indel": True}, "result": self.ROWS}
        assert sc._alignment_reportable_indels(align, 10, circular=True) == 0

    def test_a_real_count_is_still_preferred(self):
        align = {"quality": {"n_indel": 3}, "result": self.ROWS}
        assert sc._alignment_reportable_indels(align, 10, circular=True) == 3

    @pytest.mark.parametrize("bad", ["3", 3.0, None, [3]])
    def test_a_non_int_count_falls_back_to_the_extraction(self, bad):
        align = {"quality": {"n_indel": bad}, "result": self.ROWS}
        assert sc._alignment_reportable_indels(align, 10, circular=True) == 0


class TestAReadWithNothingToCompareIsNotAMismatch:
    """A trace that came back all `N` was reported `verdict: "mismatch"` at
    0% identity — "your construct is wrong" from a read that says nothing.
    `_covered_identity`'s own contract is that `(0.0, 0)` means no answer,
    never agreement; reading it as DISAGREEMENT costs a good clone."""

    REF = _rnd(1200, 4242)

    def test_an_all_no_call_read_is_inconclusive(self):
        r = _verify(self.REF, ["N" * 400], circular=True)
        assert r["verdict"] == "inconclusive"
        assert r["reads"][0]["no_data"] is True
        assert r["reads"][0]["covered_bp"] == 0

    def test_the_warning_says_the_reads_are_evidence_of_nothing(self):
        r = _verify(self.REF, ["N" * 400], circular=True)
        assert any("no comparable base" in w for w in r["warnings"])

    def test_one_good_read_beside_it_still_decides(self):
        r = _verify(self.REF, ["N" * 400, self.REF[100:600]], circular=True)
        assert r["verdict"] != "inconclusive"
        assert r["reads"][0].get("no_data") is True
        assert r["reads"][1].get("no_data") is None

    def test_a_real_mismatch_is_still_a_mismatch(self):
        r = _verify(self.REF, [_rnd(400, 99)], circular=True)
        assert r["verdict"] == "mismatch"
        assert r["reads"][0].get("no_data") is None

    def test_a_perfect_read_is_unaffected(self):
        r = _verify(self.REF, [self.REF[100:600]], circular=True)
        assert r["verdict"] == "match"
        assert r["reads"][0].get("no_data") is None
