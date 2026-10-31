"""Spliced-transcript prediction (`_predict_transcript` + `predict-transcript`).

The failure this subsystem exists to remove, in one sentence: when a construct
carries a genomic 5'UTR intron, every check written against the PLASMID sees
the unspliced leader while the biology depends on the SPLICED one, so an
upstream-ATG screen reports intronic ATGs that the spliceosome deletes — false
positives shaped exactly like real uORF hazards (reported 2026-08-24 by an
agent driving the headless API through a binary-vector build).

Two things are cross-validated against something OUTSIDE the module, because
the coordinate math here spans the origin AND both strands — the fault class
behind [INV-182]/[INV-183]/[INV-184]:

  * the translated CDS is checked against Biopython, not against the module's
    own translator;
  * the whole answer is checked for ROTATION and STRAND invariance — the same
    physical construct must answer identically however the plasmid happens to
    be drawn, which no amount of internally-consistent arithmetic can fake.
"""
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

import splicecraft as sc
import splicecraft_seqanalysis as sa


# ── A transcription unit shaped like the one that produced the bug ─────────
# promoter · leader-exon-1 · INTRON (two ATGs inside) · leader-exon-2 (one
# real uORF) · CDS · terminator. On the DNA a uATG screen sees three upstream
# ATGs; on the message there is one.
PROM = "TTGACA" + "C" * 88 + "TATAAT"
LEAD1 = "CACACACACACACACACACACACACACACACACACACAC"          # no ATG
INTRON = ("GT" + "AAGTTTCTTCTTATCTTCTT" + "ATGCCCTAA"      # intronic ATG #1
          + "TTTCTTTCTTTCTT" + "ATGAAA"                    # intronic ATG #2
          + "TTTTTTTTTTTTTTTTTTTTTTTTTGCAG")[:-2] + "AG"
LEAD2 = "ATGCCCTGATCAAACC"                                 # a REAL uORF: M P *
CDS = ("ATGGCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGGT"
       "GATGTTAATTAA")
TERM = "AATAAA" + "G" * 44
SEQ = PROM + LEAD1 + INTRON + LEAD2 + CDS + TERM
N = len(SEQ)
P1 = len(PROM)
P2 = P1 + len(LEAD1)
P3 = P2 + len(INTRON)
P4 = P3 + len(LEAD2)
P5 = P4 + len(CDS)
P6 = P5 + len(TERM)

FEATS = [
    {"start": 0,  "end": P1, "strand": 1, "type": "promoter",   "label": "35S"},
    {"start": P2, "end": P3, "strand": 1, "type": "intron",     "label": "leader intron"},
    {"start": P4, "end": P5, "strand": 1, "type": "CDS",        "label": "GFP"},
    {"start": P5, "end": P6, "strand": 1, "type": "terminator", "label": "NOS"},
]


def _predict(**kw):
    r = sa._predict_transcript(SEQ, FEATS, circular=True, **kw)
    assert r["ok"], r
    return r


def _rotate(seq, feats, k):
    n = len(seq)
    return (seq[k:] + seq[:k],
            [dict(f, start=(f["start"] - k) % n, end=(f["end"] - k) % n)
             for f in feats])


def _flip(seq, feats):
    """Reverse-complement the whole plasmid and flip every feature onto the
    other strand — the same construct, drawn the other way round."""
    n = len(seq)
    return (sc._rc(seq),
            [dict(f, start=(n - f["end"]) % n, end=(n - f["start"]) % n,
                  strand=-f["strand"]) for f in feats])


def _answer(r):
    """Everything that must not depend on how the plasmid is drawn."""
    return (r["mature_mrna"], r["cds"]["protein"], r["five_utr"]["seq"],
            tuple((u["kind"], u["protein"]) for u in r["uorfs"]),
            len(r["removed_by_splicing"]), r["kozak"]["strength"])


class TestTheReportedBug:
    """The uATG false positives, and the field that separates them out."""

    def test_the_unspliced_leader_carries_three_atgs(self):
        """The premise. If this ever stops being true the rest of the class
        is testing nothing."""
        assert SEQ[P1:P4].count("ATG") == 3

    def test_the_mature_message_carries_one(self):
        r = _predict()
        assert len(r["uorfs"]) == 1
        u = r["uorfs"][0]
        assert u["kind"] == "upstream" and u["protein"] == "MP"

    def test_the_other_two_are_named_as_removed_by_splicing(self):
        """Not merely absent — reported, with the intron that removes them, so
        a screen can show its work instead of silently disagreeing with one
        run against the DNA."""
        r = _predict()
        assert len(r["removed_by_splicing"]) == 2
        for x in r["removed_by_splicing"]:
            assert x["intron"] == "leader intron"
            assert SEQ[x["genomic"]:x["genomic"] + 3] == "ATG"
            assert P2 <= x["genomic"] < P3      # genuinely inside the intron

    def test_removed_atgs_are_absent_from_the_mature_leader(self):
        r = _predict()
        assert r["five_utr"]["seq"].count("ATG") == 1


class TestMatureMessage:
    def test_it_is_the_exons_concatenated(self):
        """Ground truth built by hand, not by re-running the module."""
        r = _predict()
        assert r["mature_mrna"] == LEAD1 + LEAD2 + CDS + TERM
        assert r["spliced"] is True

    def test_the_intron_is_found_with_its_real_boundaries(self):
        r = _predict()
        assert len(r["introns"]) == 1
        i = r["introns"][0]
        assert (i["donor"], i["acceptor"], i["canonical"]) == ("GT", "AG", True)
        assert i["length"] == len(INTRON)

    def test_protein_matches_biopython(self):
        r = _predict()
        assert r["cds"]["protein"] == str(Seq(CDS).translate(to_stop=True))

    def test_utrs_partition_the_message(self):
        r = _predict()
        assert (r["five_utr"]["seq"] + r["mature_mrna"][r["cds"]["start"]:
                                                        r["cds"]["end"]]
                + r["three_utr"]["seq"]) == r["mature_mrna"]

    def test_an_unspliced_unit_still_works(self):
        """No intron annotated → mature == pre, and nothing is claimed to be
        removed."""
        feats = [f for f in FEATS if f["type"] != "intron"]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] and r["spliced"] is False
        assert r["mature_mrna"] == r["pre_mrna"]
        assert r["removed_by_splicing"] == []
        assert len(r["uorfs"]) == 3         # what the DNA-only screen sees


class TestCoordinateInvariance:
    """The origin-and-strand fault class. Internally-consistent arithmetic can
    still be wrong; answering identically from every drawing cannot."""

    def test_rotation_invariant(self):
        ref = _answer(_predict())
        for k in range(0, N, 7):
            seq, feats = _rotate(SEQ, FEATS, k)
            r = sa._predict_transcript(seq, feats, circular=True)
            assert r["ok"], (k, r)
            assert _answer(r) == ref, f"rotation {k} changed the answer"

    def test_strand_invariant(self):
        ref = _answer(_predict())
        seq, feats = _flip(SEQ, FEATS)
        r = sa._predict_transcript(seq, feats, circular=True)
        assert r["ok"] and r["unit"]["strand"] == -1
        assert _answer(r) == ref

    def test_minus_strand_across_the_origin(self):
        """Both at once — the case that actually breaks implementations."""
        ref = _answer(_predict())
        for k in range(0, N, 11):
            seq, feats = _rotate(*_flip(SEQ, FEATS), k)
            r = sa._predict_transcript(seq, feats, circular=True)
            assert r["ok"], (k, r)
            assert _answer(r) == ref, f"minus-strand rotation {k} differs"

    def test_a_linear_record_refuses_a_wrapping_unit(self):
        seq, feats = _rotate(SEQ, FEATS, P5 + 20)
        r = sa._predict_transcript(seq, feats, circular=False)
        assert r["ok"] is False and "LINEAR" in r["error"]


class TestALinearUnitMayEndOnTheLastBase:
    """A linear record's last edge is `len(seq)` itself. Folding it `% len`
    moved a terminator that ends on the last base to bp 0, and the unit was
    refused as "crossing the origin" (audit 2026-10-23)."""

    def test_the_reference_unit_reads_the_same_linear(self):
        assert FEATS[-1]["end"] == N        # the premise: it ends on the end
        r = sa._predict_transcript(SEQ, FEATS, circular=False)
        assert r["ok"], r
        assert _answer(r) == _answer(_predict())

    def test_its_mirror_reads_the_same_linear(self):
        flipped = [dict(f, start=N - f["end"], end=N - f["start"],
                        strand=-f["strand"]) for f in FEATS]
        r = sa._predict_transcript(sc._rc(SEQ), flipped, circular=False)
        assert r["ok"] and r["unit"]["strand"] == -1, r
        assert _answer(r) == _answer(_predict())

    @pytest.mark.parametrize("strand", [1, -1])
    def test_an_explicit_tx_end_at_the_length(self, strand):
        feats = [f for f in FEATS if f["type"] in ("CDS", "intron")]
        seq = SEQ
        lo, hi = P1, N
        if strand < 0:
            seq = sc._rc(SEQ)
            feats = [dict(f, start=N - f["end"], end=N - f["start"],
                          strand=-1) for f in feats]
            lo, hi = 0, N - P1
        r = sa._predict_transcript(seq, feats, circular=False,
                                   tx_start=lo, tx_end=hi)
        assert r["ok"], r
        assert _answer(r) == _answer(_predict())

    def test_a_cds_filling_the_whole_molecule(self):
        feats = [{"start": 0, "end": len(CDS), "strand": 1, "type": "CDS",
                  "label": "whole"}]
        r = sa._predict_transcript(CDS, feats, circular=False,
                                   tx_start=0, tx_end=len(CDS))
        assert r["ok"], r
        assert r["mature_mrna"] == CDS
        assert r["cds"]["protein"] == _predict()["cds"]["protein"]


class TestIntronSources:
    def test_gap_between_location_parts_is_an_intron(self):
        """The standard GenBank way to spell a spliced feature — no `intron`
        feature involved."""
        feats = [f for f in FEATS if f["type"] != "intron"]
        feats = [dict(f, start=P1, end=P5, parts=[[P1, P2], [P3, P5]])
                 if f["type"] == "CDS" else f for f in feats]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] and len(r["introns"]) == 1
        assert r["mature_mrna"] == _predict()["mature_mrna"]

    def test_overlapping_intron_annotations_merge(self):
        feats = FEATS + [{"start": P2 + 10, "end": P3, "strand": 1,
                          "type": "intron", "label": "dup"}]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] and len(r["introns"]) == 1

    def test_a_cds_start_inside_an_intron_is_refused(self):
        """An inconsistent annotation must fail, not yield half a message."""
        feats = [dict(f, start=P2 + 6, end=P3 + 6) if f["type"] == "intron"
                 else f for f in FEATS]
        feats = [dict(f, start=P3, end=P5) if f["type"] == "CDS" else f
                 for f in feats]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] is False and "intron" in r["error"]


class TestKozak:
    def test_the_real_start_scores_strong(self):
        r = _predict()
        k = r["kozak"]
        assert k["context"].endswith("ATGG")
        assert (k["minus3"], k["plus4"], k["strength"]) == ("A", "G", "strong")

    def test_scoring_is_the_documented_two_position_rule(self):
        # The -3 base is index 3 of the six upstream bases; +4 is the one
        # after the ATG. Purine at -3 and G at +4 are worth one point each.
        assert sa._kozak_context("TTTATT" + "ATG" + "G", 6)["score"] == 2
        assert sa._kozak_context("TTTGTT" + "ATG" + "G", 6)["score"] == 2
        assert sa._kozak_context("TTTATT" + "ATG" + "T", 6)["score"] == 1
        assert sa._kozak_context("TTTTTT" + "ATG" + "G", 6)["score"] == 1
        assert sa._kozak_context("TTTTTT" + "ATG" + "T", 6)["score"] == 0
        assert sa._kozak_context("TTTCTT" + "ATG" + "C", 6)["score"] == 0

    def test_the_rule_and_the_matrix_each_name_their_own_basis(self):
        """Both answers ship in one dict and they CAN disagree, so neither may
        read as "the answer" on its own. `score` / `strength` stay the -3/+4
        consensus rule; `pwm` is the calibrated matrix. The contract and the
        matrix itself live in `splicecraft_kozak` (see `tests/test_kozak.py`);
        what this pins is that the transcript reader surfaces both."""
        k = sa._kozak_context("AAAAAAGCAACA" + "ATG" + "G", 12)
        assert "consensus" in k["rule_model"]
        assert "PWM" in k["model"] and "AUC" in k["model"]
        assert k["pwm"]["percentile"] > 80
        assert k["pwm"]["tier"] == "strong"
        # A context the rule calls adequate and the matrix calls weak.
        weak = sa._kozak_context("TTTTTTTTTTTT" + "ATG" + "G", 12)
        assert weak["strength"] == "adequate" and weak["pwm"]["tier"] == "weak"

    def test_the_matrix_abstains_rather_than_scoring_a_partial_window(self):
        k = sa._kozak_context("AT" + "ATGG", 2)
        assert k["pwm"]["score"] is None and k["pwm"]["note"]
        assert k["score"] == 1          # the rule still answers

    def test_a_short_leader_is_marked_truncated_not_padded(self):
        k = sa._kozak_context("AT" + "ATGG", 2)
        assert k["truncated"] is True and len(k["context"]) < 10


class TestUorfRanking:
    """`rank` 1 is the upstream ATG a scanning ribosome is most likely to stop
    at — the ordering the -3/+4 rule structurally cannot produce, since 38% of
    real starts and 13% of upstream ATGs share its top bucket."""

    def _leader(self, leader):
        seq = PROM + leader + CDS + TERM
        p1, p2 = len(PROM), len(PROM) + len(leader)
        feats = [
            {"start": 0, "end": p1, "strand": 1, "type": "promoter", "label": "p"},
            {"start": p2, "end": p2 + len(CDS), "strand": 1, "type": "CDS", "label": "c"},
            {"start": p2 + len(CDS), "end": len(seq), "strand": 1,
             "type": "terminator", "label": "t"},
        ]
        r = sa._predict_transcript(seq, feats, circular=True)
        assert r["ok"], r
        return r

    def test_rank_follows_the_calibrated_score_not_the_position(self):
        # Three uORFs: a weak context FIRST, a strong one SECOND. Ranking by
        # position would put them the other way round.
        r = self._leader("TTTTTTTTTTTT" + "ATGCCCTGA"
                         + "AAAAAAGCAACA" + "ATGCCCTGA" + "CACACACA")
        ranked = sorted(r["uorfs"], key=lambda u: u["rank"])
        assert [u["rank"] for u in ranked] == [1, 2]
        assert ranked[0]["position"] > ranked[1]["position"]
        assert (ranked[0]["kozak"]["pwm"]["score"]
                > ranked[1]["kozak"]["pwm"]["score"])

    def test_the_rule_cannot_order_what_the_matrix_can(self):
        """Both uORFs here sit in the SAME rule bucket, so the rule is silent
        and the matrix is not. This is the capability the model adds."""
        r = self._leader("TTTTTTTTTATT" + "ATGCCCTGA"
                         + "AAAAAAGCAACA" + "ATGCCCTGA" + "CACACACA")
        assert len(r["uorfs"]) == 2
        assert len({u["kozak"]["score"] for u in r["uorfs"]}) == 1
        assert len({u["rank"] for u in r["uorfs"]}) == 2

    def test_an_atg_too_close_to_the_cap_is_unranked_and_warned_about(self):
        """A partial window sums fewer log-odds terms, so it cannot be ordered
        against a full one. The count is warned about rather than left to be
        inferred from a gap in the numbering."""
        r = self._leader("ATGCCCTGA" + "CACACACA")
        u = r["uorfs"][0]
        assert u["position"] < 12
        assert u["rank"] is None
        assert u["kozak"]["pwm"]["score"] is None
        assert any("cannot be scored" in w for w in r["warnings"]), r["warnings"]
        # The warning counts them; the CAUSE is per-ATG, because an ambiguous
        # base in the window produces the same `None` as a short 5' end.
        assert any("kozak.pwm.note" in w for w in r["warnings"])

    def test_a_scoreable_uorf_is_not_warned_about(self):
        r = self._leader("AAAAAAGCAACA" + "ATGCCCTGA" + "CACACACA")
        assert r["uorfs"][0]["rank"] == 1
        assert not any("cannot score" in w for w in r["warnings"])


class TestWeakStartWarning:
    """A weak context on the CDS ITSELF is the finding worth surfacing: it is
    fixable by changing the leader. The trigger is anchored on the measured
    real-start distribution, not picked."""

    def _with_leader(self, leader):
        seq = PROM + leader + CDS + TERM
        p1, p2 = len(PROM), len(PROM) + len(leader)
        feats = [
            {"start": 0, "end": p1, "strand": 1, "type": "promoter", "label": "p"},
            {"start": p2, "end": p2 + len(CDS), "strand": 1, "type": "CDS", "label": "c"},
            {"start": p2 + len(CDS), "end": len(seq), "strand": 1,
             "type": "terminator", "label": "t"},
        ]
        return sa._predict_transcript(seq, feats, circular=True)

    def test_it_fires_on_a_context_poorer_than_real_plant_starts(self):
        r = self._with_leader("T" * 20)
        assert r["kozak"]["pwm"]["tier"] == "weak", r["kozak"]["pwm"]
        assert any("real plant start codons" in w for w in r["warnings"]), \
            r["warnings"]

    def test_it_stays_quiet_on_a_good_context(self):
        r = self._with_leader("CACACACAAAAAAAGCAACA")
        assert r["kozak"]["pwm"]["tier"] != "weak"
        assert not any("real plant start" in w for w in r["warnings"])


class TestUorfClassification:
    def _with_leader(self, leader):
        seq = PROM + leader + CDS + TERM
        p1, p2 = len(PROM), len(PROM) + len(leader)
        feats = [
            {"start": 0, "end": p1, "strand": 1, "type": "promoter", "label": "p"},
            {"start": p2, "end": p2 + len(CDS), "strand": 1, "type": "CDS", "label": "c"},
            {"start": p2 + len(CDS), "end": len(seq), "strand": 1,
             "type": "terminator", "label": "t"},
        ]
        r = sa._predict_transcript(seq, feats, circular=True)
        assert r["ok"], r
        return r

    def test_a_uorf_that_stops_before_the_cds(self):
        r = self._with_leader("ATGCCCTGA" + "CACACACA")
        assert [u["kind"] for u in r["uorfs"]] == ["upstream"]

    def test_an_out_of_frame_overlap_is_called_out(self):
        # ATG with no stop before the CDS start, and out of frame with it.
        r = self._with_leader("ATG" + "CACACACAC")     # 12 nt, 12 % 3 == 0…
        kinds = [u["kind"] for u in r["uorfs"]]
        assert kinds and kinds[0] in ("in_frame_extension",
                                      "out_of_frame_overlap")

    def test_an_in_frame_upstream_atg_is_an_n_terminal_extension(self):
        r = self._with_leader("ATG" + "CACACA")        # 9 nt → in frame
        assert r["uorfs"][0]["kind"] == "in_frame_extension"

    def test_min_uorf_aa_filters(self):
        r = self._with_leader("ATGCCCTGA" + "CACACACA")
        assert len(r["uorfs"]) == 1
        seq_feats = sa._predict_transcript(
            PROM + "ATGCCCTGA" + "CACACACA" + CDS + TERM,
            [{"start": 0, "end": len(PROM), "strand": 1, "type": "promoter", "label": "p"},
             {"start": len(PROM) + 17, "end": len(PROM) + 17 + len(CDS),
              "strand": 1, "type": "CDS", "label": "c"},
             {"start": len(PROM) + 17 + len(CDS),
              "end": len(PROM + "ATGCCCTGA" + "CACACACA" + CDS + TERM),
              "strand": 1, "type": "terminator", "label": "t"}],
            circular=True, min_uorf_aa=5)
        assert seq_feats["ok"] and seq_feats["uorfs"] == []


class TestSpliceAdvisory:
    def test_annotated_boundaries_are_not_reported_as_cryptic(self):
        """The annotated donor/acceptor obviously score well — reporting them
        as cryptic hazards would bury any real one."""
        r = _predict()
        assert r["splice"]["skipped"] is False
        assert r["splice"]["n_sites"] >= 2
        assert r["splice"]["n_cryptic"] == 0

    def test_it_can_be_turned_off(self):
        r = _predict(check_splice=False)
        assert r["splice"]["skipped"] is True

    def test_an_unavailable_model_reads_as_not_checked(self, monkeypatch):
        """Never an empty list, which would read as clean."""
        def _boom(*a, **k):
            raise RuntimeError("splice model is not trained")
        monkeypatch.setattr(sa, "_splice_scan", _boom)
        r = _predict()
        assert r["splice"]["skipped"] is True
        assert "not trained" in r["splice"]["reason"]
        assert any("SKIPPED" in w for w in r["warnings"])


class TestFeatureSelection:
    def test_by_label(self):
        r = _predict(cds="GFP", promoter="35S", terminator="NOS")
        assert r["unit"]["cds"] == "GFP"

    def test_an_explicit_label_that_misses_is_an_error(self):
        """Never a quiet fall-back to a type guess — that would report a
        transcript for a different feature than the one asked for."""
        r = sa._predict_transcript(SEQ, FEATS, circular=True, cds="nope")
        assert r["ok"] is False and "nope" in r["error"]

    def test_explicit_bounds_replace_the_anchors(self):
        feats = [f for f in FEATS if f["type"] in ("intron", "CDS")]
        r = sa._predict_transcript(SEQ, feats, circular=True,
                                    tx_start=P1, tx_end=P6)
        assert r["ok"] and r["mature_mrna"] == _predict()["mature_mrna"]

    def test_no_anchors_and_no_bounds_says_what_to_do(self):
        feats = [f for f in FEATS if f["type"] in ("intron", "CDS")]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] is False and "tx_start" in r["error"]

    def test_no_cds_at_all(self):
        feats = [f for f in FEATS if f["type"] != "CDS"]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] is False and "CDS" in r["error"]


# ── The agent endpoint ─────────────────────────────────────────────────────
class _App:
    def __init__(self, rec):
        self._current_record = rec


def _record(rot=0):
    s = SEQ[rot:] + SEQ[:rot]
    rec = SeqRecord(Seq(s), id="tx", name="tx",
                    annotations={"molecule_type": "DNA", "topology": "circular"})
    def loc(a, b):
        a, b = (a - rot) % N, (b - rot) % N
        if b > a:
            return FeatureLocation(a, b, strand=1)
        return CompoundLocation([FeatureLocation(a, N, strand=1),
                                 FeatureLocation(0, b, strand=1)])
    for a, b, t, lab in ((0, P1, "promoter", "35S"),
                         (P2, P3, "intron", "leader intron"),
                         (P4, P5, "CDS", "GFP"),
                         (P5, P6, "terminator", "NOS")):
        rec.features.append(
            SeqFeature(loc(a, b), type=t, qualifiers={"label": [lab]}))
    return rec


class TestPredictTranscriptEndpoint:
    def test_it_answers_on_the_loaded_record(self):
        r = sc._h_predict_transcript(_App(_record()), {})
        assert not isinstance(r, tuple), r
        assert r["mature_mrna"] == LEAD1 + LEAD2 + CDS + TERM
        assert len(r["uorfs"]) == 1 and len(r["removed_by_splicing"]) == 2

    def test_a_wrapping_feature_is_not_mistaken_for_a_spliced_one(self):
        """THE trap in the record→dict step: an origin-spanning feature is a
        two-part CompoundLocation exactly like a spliced one (sacred invariant
        #9), so a naive parts-are-exons reading invents an intron spanning the
        whole backbone."""
        rec = _record(rot=P5 + 20)          # drags the terminator across bp 0
        assert any(len(getattr(f.location, "parts", [1])) > 1
                   for f in rec.features), "fixture no longer has a wrap"
        dicts = sc._agent_transcript_feature_dicts(rec)
        wrapped = [d for d in dicts if d["end"] < d["start"]]
        assert wrapped, "expected an end<start feature"
        assert all("parts" not in d for d in wrapped)
        r = sc._h_predict_transcript(_App(rec), {})
        assert not isinstance(r, tuple), r
        assert len(r["introns"]) == 1
        assert r["mature_mrna"] == LEAD1 + LEAD2 + CDS + TERM

    def test_no_record_is_422(self):
        assert sc._h_predict_transcript(_App(None), {})[1] == 422

    def test_an_unreadable_unit_is_422_with_the_reason(self):
        out = sc._h_predict_transcript(_App(_record()), {"cds": "nope"})
        assert isinstance(out, tuple) and out[1] == 422
        assert "nope" in out[0]["error"]

    def test_half_a_bound_pair_is_400(self):
        out = sc._h_predict_transcript(_App(_record()), {"tx_start": 5})
        assert isinstance(out, tuple) and out[1] == 400

    @pytest.mark.parametrize("bad", [{"cds": True}, {"cds": 1.5},
                                     {"min_uorf_aa": "lots"}])
    def test_bad_parameter_types_are_400(self, bad):
        out = sc._h_predict_transcript(_App(_record()), bad)
        assert isinstance(out, tuple) and out[1] == 400, bad

    def test_registered_and_read_only(self):
        fn, write = sc._state._AGENT_HANDLERS["predict-transcript"]
        assert fn is sc._h_predict_transcript and write is False


class TestMalformedInput:
    """Edge-case sweep. A record can carry a fuzzy or hand-edited location,
    and a caller can build the feature list itself — neither may crash out of
    the middle of the coordinate walk."""

    def test_a_non_integer_coordinate_is_skipped_not_raised(self):
        feats = [dict(f, start="oops") if f["type"] == "CDS" else f
                 for f in FEATS]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] is False                      # no CDS left to use
        assert any("skipped" in w for w in r["warnings"])

    def test_one_bad_feature_does_not_sink_the_rest(self):
        feats = FEATS + [{"start": None, "end": 5, "strand": 1,
                          "type": "misc_feature", "label": "fuzzy"}]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] is True
        assert r["mature_mrna"] == LEAD1 + LEAD2 + CDS + TERM
        assert any("skipped" in w for w in r["warnings"])

    def test_non_dict_entries_are_skipped(self):
        r = sa._predict_transcript(SEQ, FEATS + ["nonsense", None, 42],
                                    circular=True)
        assert r["ok"] is True and any("skipped" in w for w in r["warnings"])

    @pytest.mark.parametrize("parts", [
        [["a", "b"]], [[1]], "notalist", [{"start": 1}], None, [],
    ])
    def test_garbage_parts_fall_back_to_the_single_span(self, parts):
        """A malformed `parts` must read as 'not spliced', never as an intron
        somewhere arbitrary."""
        feats = [dict(f, parts=parts) if f["type"] == "CDS" else f
                 for f in FEATS if f["type"] != "intron"]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] is True and r["introns"] == []

    def test_an_intron_swallowing_the_cds_is_refused(self):
        feats = FEATS + [{"start": P1, "end": P6, "strand": 1,
                          "type": "intron", "label": "everything"}]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] is False and "intron" in r["error"]

    def test_a_zero_length_cds_is_refused(self):
        feats = [dict(f, end=f["start"]) if f["type"] == "CDS" else f
                 for f in FEATS]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert r["ok"] is False

    def test_a_cds_with_no_utr_on_either_side(self):
        seq = PROM + CDS + TERM
        p1 = len(PROM)
        feats = [
            {"start": 0, "end": p1, "strand": 1, "type": "promoter", "label": "p"},
            {"start": p1, "end": p1 + len(CDS), "strand": 1, "type": "CDS", "label": "c"},
            {"start": p1 + len(CDS), "end": len(seq), "strand": 1,
             "type": "terminator", "label": "t"},
        ]
        r = sa._predict_transcript(seq, feats, circular=True)
        assert r["ok"] and r["five_utr"]["length"] == 0
        assert r["uorfs"] == [] and r["kozak"]["truncated"] is True

    def test_a_large_unit_stays_fast(self):
        """The advisory PWM scan runs over the whole pre-mRNA; a 60 kb unit
        must not turn a read endpoint into a stall."""
        import time
        lead = ("ACGTTGCA" * 7500)[:60_000]
        seq = PROM + lead + CDS + TERM
        p1 = len(PROM)
        feats = [
            {"start": 0, "end": p1, "strand": 1, "type": "promoter", "label": "p"},
            {"start": p1 + len(lead), "end": p1 + len(lead) + len(CDS),
             "strand": 1, "type": "CDS", "label": "c"},
            {"start": p1 + len(lead) + len(CDS), "end": len(seq), "strand": 1,
             "type": "terminator", "label": "t"},
        ]
        t = time.perf_counter()
        r = sa._predict_transcript(seq, feats, circular=True)
        assert r["ok"] and time.perf_counter() - t < 5.0


class TestIncludeSequences:
    def test_off_drops_the_bulk_but_keeps_every_verdict(self):
        rec = _record()
        full = sc._h_predict_transcript(_App(rec), {})
        lean = sc._h_predict_transcript(_App(rec), {"include_sequences": False})
        assert "mature_mrna" not in lean and "pre_mrna" not in lean
        assert "seq" not in lean["five_utr"] and "seq" not in lean["three_utr"]
        # The answers themselves are untouched.
        for key in ("uorfs", "removed_by_splicing", "introns", "exons",
                    "kozak", "cds", "unit"):
            assert lean[key] == full[key], key
        assert lean["five_utr"]["length"] == full["five_utr"]["length"]

    def test_on_by_default(self):
        r = sc._h_predict_transcript(_App(_record()), {})
        assert r["mature_mrna"] and r["pre_mrna"]


# ── A longer unit, for the premature-3'-end scan ───────────────────────────
# The default fixture's mature message is under the 160 nt scoring window in
# places and its hand-typed "CDS" is AT-rich, which trips the scanner for the
# wrong reason. This one is built from the GenBank-verified preset catalogue
# instead: EGFP (720 bp, 62% GC — scores no 3'-end regions on its own) and the
# NOS terminator (266 bp, 36% GC — scores one), which is what a real plant
# cassette actually looks like.
import splicecraft_presets as _presets


def _preset_seq(name: str) -> str:
    for f in _presets._preset_features():
        if f["name"] == name:
            return f["sequence"]
    raise AssertionError(f"preset {name!r} missing")


_L_PROM = _preset_seq("CaMV 35S promoter")
_L_LEAD1 = "CGCCGCCACCGCCGCCACCGCCGCCACCGCCGCCACC"        # GC-rich, no ATG
# A plant intron: GT..AG, long, and AT/U-rich — which is exactly what the
# 3'-end score keys on. A PRE-mRNA scan would flag this; the mature one
# must not.
_L_INTRON = ("GT" + ("AAGTTTCTTCTTATCTTCTTTTTCTTTCTTTCTT"
                     "ATTTATTTATTTTTATTCTTTTCTTTCTTATCT" * 4) + "AG")
_L_LEAD2 = "CCACCGCCACC"
_L_CDS = _preset_seq("EGFP")                             # ends in TAA
# In-frame, stop-free and 100% AT — 30 codons of Asn-Leu, a stretch of CDS
# that reads like the end of a message.
_L_CDS_ATRICH = "AATTTA" * 30
_L_TERM = _preset_seq("NOS terminator")


def _long_unit(*, at_in_cds: bool = False):
    """``(seq, feats)`` for the long cassette, optionally with an AT-rich
    block spliced into the middle of the CDS, in frame."""
    cds = _L_CDS
    if at_in_cds:
        cut = 300                                        # a codon boundary
        cds = cds[:cut] + _L_CDS_ATRICH + cds[cut:]
    seq = _L_PROM + _L_LEAD1 + _L_INTRON + _L_LEAD2 + cds + _L_TERM
    p1 = len(_L_PROM)
    p2 = p1 + len(_L_LEAD1)
    p3 = p2 + len(_L_INTRON)
    p4 = p3 + len(_L_LEAD2)
    p5 = p4 + len(cds)
    feats = [
        {"start": 0, "end": p1, "strand": 1, "type": "promoter", "label": "35S"},
        {"start": p2, "end": p3, "strand": 1, "type": "intron", "label": "leader intron"},
        {"start": p4, "end": p5, "strand": 1, "type": "CDS", "label": "EGFP"},
        {"start": p5, "end": p5 + len(_L_TERM), "strand": 1,
         "type": "terminator", "label": "NOS"},
    ]
    return seq, feats


def _predict_long(*, at_in_cds: bool = False, **kw):
    seq, feats = _long_unit(at_in_cds=at_in_cds)
    r = sa._predict_transcript(seq, feats, circular=True, **kw)
    assert r["ok"], r
    return r


def _short_unit():
    """A cassette whose mature message is under the 160 nt scoring window."""
    prom, cds, term = _L_PROM, "ATG" + "GCT" * 20 + "TAA", "GCGC" * 10
    seq = prom + cds + term
    feats = [
        {"start": 0, "end": len(prom), "strand": 1, "type": "promoter", "label": "p"},
        {"start": len(prom), "end": len(prom) + len(cds), "strand": 1,
         "type": "CDS", "label": "tiny"},
        {"start": len(prom) + len(cds), "end": len(seq), "strand": 1,
         "type": "terminator", "label": "t"},
    ]
    return seq, feats


class TestPrematureThreePrimeEnd:
    """A stretch of the message that a plant reads as the END of it.

    The protein is never finished and the DNA sequences perfectly — the same
    silent failure mode the cryptic-splice advisory exists for."""

    def test_the_three_utr_is_not_reported_as_premature(self):
        """The terminator is SUPPOSED to look like a message end.

        Flagging it would bury every real hit under one that is working as
        designed, on every correctly-built cassette."""
        r = _predict_long()
        tp = r["three_prime"]
        assert tp["skipped"] is False
        assert tp["n_regions"] >= 1, "the AT-rich NOS terminator should score"
        assert tp["n_premature"] == 0, (
            f"a 3'UTR hit is the terminator doing its job: {tp['premature']}")
        assert not any("3' end" in w for w in r["warnings"])

    def test_an_at_rich_block_inside_the_cds_is_flagged(self):
        r = _predict_long(at_in_cds=True)
        tp = r["three_prime"]
        assert tp["n_premature"] >= 1, tp
        cds_end = r["cds"]["end"]
        assert all(x["start"] < cds_end for x in tp["premature"])
        assert any("upstream of the stop codon" in w for w in r["warnings"])

    def test_an_at_rich_intron_is_not_flagged(self):
        """THE design decision, asserted rather than assumed.

        Plant introns are AT- and U-rich — exactly what the 3'-end score keys
        on — so scanning the PRE-mRNA would report every intron as a hazard.
        This proves the mature message is what gets scanned, by showing the
        pre-mRNA scan WOULD have fired on the intron."""
        from splicecraft_regulettuce import _scan_3prime_risk
        r = _predict_long()
        pre_hits = _scan_3prime_risk(r["pre_mrna"], circular=False,
                                     both_strands=False)
        intron = r["introns"][0]
        in_intron = [h for h in pre_hits
                     if h["start"] < intron["end"]
                     and h["end"] > intron["start"]]
        assert in_intron, (
            "fixture is not exercising the trap — the intron must be AT-rich "
            "enough that a pre-mRNA scan would flag it")
        assert r["three_prime"]["n_premature"] == 0, (
            "the intron was flagged, so the scan is reading the pre-mRNA")

    def test_a_message_shorter_than_the_window_says_so(self):
        """Not checked must never read as clean."""
        seq, feats = _short_unit()
        r = sa._predict_transcript(seq, feats, circular=True)
        assert r["ok"], r
        tp = r["three_prime"]
        assert tp["skipped"] is True
        assert "shorter than" in tp["reason"]
        assert "160" in tp["reason"]

    def test_it_can_be_turned_off(self):
        r = _predict_long(check_three_prime=False)
        assert r["three_prime"]["skipped"] is True
        assert r["three_prime"]["reason"] == "not requested"

    def test_a_failing_scan_reads_as_not_checked(self, monkeypatch):
        import splicecraft_regulettuce as rl

        def _boom(*a, **k):
            raise ValueError("scanner exploded")
        monkeypatch.setattr(rl, "_scan_3prime_risk", _boom)
        r = _predict_long()
        assert r["three_prime"]["skipped"] is True
        assert "exploded" in r["three_prime"]["reason"]
        assert any("SKIPPED" in w for w in r["warnings"])

    def test_the_scan_is_single_stranded(self):
        """An mRNA is read one way, so a 3'-end-like composition on the
        complement is not a hazard — the rule `_splice_scan` already
        applies."""
        from splicecraft_regulettuce import _scan_3prime_risk
        r = _predict_long(at_in_cds=True)
        one = _scan_3prime_risk(r["mature_mrna"], circular=False,
                                both_strands=False)
        got = {(h["start"], h["length"]) for h in r["three_prime"]["regions"]}
        assert got == {(h["start"], h["length"]) for h in one}, (
            "the transcript scan must use the single-strand result")


class TestThreePrimeEndpoint:
    def test_endpoint_passes_the_flag_through(self):
        app = _App(_record())
        on = sc._state._AGENT_HANDLERS["predict-transcript"][0](app, {})
        off = sc._state._AGENT_HANDLERS["predict-transcript"][0](
            app, {"check_three_prime": False})
        for res in (on, off):
            assert not isinstance(res, tuple), res
        assert off["three_prime"]["skipped"] is True
        assert off["three_prime"]["reason"] == "not requested"
        assert "three_prime" in on


# ══════════════════════════════════════════════════════════════════════════
#  The initiation-model endpoints
# ══════════════════════════════════════════════════════════════════════════

def _call(name, payload, app=None):
    return sc._state._AGENT_HANDLERS[name][0](app, payload)


class TestInitiationEndpoints:
    def test_both_are_registered_read_only(self):
        """A `write=True` endpoint is refused outright by a `--read-only`
        attach, which is how the GUI-is-open case is served."""
        for name in ("score-initiation-context", "initiation-model"):
            entry = sc._state._AGENT_HANDLERS.get(name)
            assert entry is not None, f"{name} not registered"
            assert entry[1] is False, f"{name} is registered as a write"

    def test_the_model_endpoint_publishes_its_refusals(self):
        m = _call("initiation-model", {})["model"]
        assert m["trained"] is True and m["window"] == "-12..+4"
        assert m["auc_vs_utr_uatg"] > m["rule_auc_vs_utr_uatg"]
        assert set(m["refused"]) >= {"efficiency_rate", "partial_windows",
                                     "window_past_plus4"}

    def test_scoring_defaults_to_the_first_atg_and_ranks_the_rest(self):
        seq = ("AAAAAAGCAACA" + "ATGGCT"          # strong context, pos 12
               + "TTTTTTTTTTTT" + "ATGCCC")       # weak context, pos 30
        r = _call("score-initiation-context", {"sequence": seq})
        assert r["position"] == 12 and r["codon"] == "ATG"
        assert r["n_atgs"] == 2
        strong = next(a for a in r["atgs"] if a["position"] == 12)
        weak = next(a for a in r["atgs"] if a["position"] == 30)
        assert (strong["rank"], weak["rank"]) == (1, 2)
        assert strong["percentile"] > weak["percentile"]
        assert r["context"]["pwm"]["tier"] == "strong"

    def test_a_position_without_an_atg_answers_rather_than_refusing(self):
        """The caller asked about an offset; the honest answer names the codon
        that is actually there instead of a 400 that hides it."""
        r = _call("score-initiation-context",
                  {"sequence": "A" * 12 + "CTGG" + "ATGG", "position": 12})
        assert r["codon"] == "CTG"
        assert r["context"]["pwm"]["score"] is None
        assert "no ATG" in r["context"]["pwm"]["note"]

    def test_an_atg_too_close_to_the_five_prime_end_is_unranked(self):
        r = _call("score-initiation-context", {"sequence": "AT" + "ATGGCCGGC"})
        row = r["atgs"][0]
        assert row["rank"] is None and row["pwm_score"] is None
        assert str(sa._KOZAK_UP) in row["note"]

    @pytest.mark.parametrize("payload,frag", [
        ({}, "sequence"),
        ({"sequence": "ATG"}, "at least 4"),
        ({"sequence": "ATGGZZ"}, "non-IUPAC"),
        ({"sequence": "CCCCCCCC"}, "no ATG"),
        ({"sequence": "ATGGCC", "position": 99}, "'position' must be"),
        ({"sequence": "ATGGCC", "position": -1}, "'position' must be"),
    ])
    def test_malformed_input_is_refused_with_a_reason(self, payload, frag):
        res = _call("score-initiation-context", payload)
        assert isinstance(res, tuple) and res[1] == 400, res
        assert frag in res[0]["error"], res[0]["error"]

    def test_the_atg_list_is_capped_and_says_so(self):
        r = _call("score-initiation-context",
                  {"sequence": "ATGGCC" * 50, "limit": 10})
        assert len(r["atgs"]) == 10 and r["n_atgs"] == 50
        assert "rank` orders only those" in r["truncated"]


class TestStartCodonIsAlwaysInTheMessage:
    """`_predict_transcript` reads the Kozak context, the translation and
    every uORF frame off `mature[cds_start_m:]`.

    If `cds_start_m` could ever equal `len(mature)` — which is what
    `_mature_pos` returns for a pre-mRNA offset at or past the transcript's
    end — all three would be describing the bases at the far END of the
    message as though they were the start codon, and the old `_kozak_context`
    quietly obliged by reading `minus3` from there. The guard that refuses
    that leans entirely on `_tx_map_span` never handing back a start at or
    past `pre_len`, so THAT is what gets fuzzed here: the property, not the
    branch it protects."""

    def test_tx_map_span_never_starts_at_or_past_the_unit_end(self):
        import itertools
        rng = range(-3, 9)
        for gs, ge, a, pre_len, strand, circ in itertools.product(
                rng, rng, (0, 2, 5), (0, 1, 4, 8), (1, -1), (True, False)):
            got = sa._tx_map_span(gs, ge, a, 8, strand, pre_len,
                                  circular=circ)
            if got is None:
                continue
            s, e = got
            assert 0 <= s < pre_len, (gs, ge, a, pre_len, strand, circ, got)
            assert s < e <= pre_len, (gs, ge, a, pre_len, strand, circ, got)

    def test_the_cds_start_always_indexes_into_the_mature_message(self):
        """Same claim, through the real function on every rotation: the start
        codon the result reports must be a base the message actually has."""
        for k in (0, 1, 37, len(SEQ) // 2, len(SEQ) - 1):
            seq, feats = _rotate(SEQ, FEATS, k)
            r = sa._predict_transcript(seq, feats, circular=True)
            assert r["ok"], (k, r)
            mature = r["mature_mrna"]
            start = r["cds"]["start"]
            assert 0 <= start < len(mature)
            assert mature[start:start + 3] == "ATG", k
