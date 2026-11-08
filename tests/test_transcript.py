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
        t = time.process_time()
        r = sa._predict_transcript(seq, feats, circular=True)
        spent = time.process_time() - t
        # CPU time: the claim is about the algorithm, and a wall-clock budget
        # under `-n auto` fails on starvation instead.
        assert r["ok"]
        assert spent < 5.0, f"predict took {spent:.3f} s CPU (budget 5.0)"


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


# ══════════════════════════════════════════════════════════════════════════
#  Field report 2026-11-08 — three things the answer used to leave unsaid
# ══════════════════════════════════════════════════════════════════════════
#
# Reported from a real plant construct driven through the headless API:
#
#   * a genomic promoter feature that SPANS the transcribed 5'UTR and its
#     intron (which is how plant polyubiquitin promoters are usually
#     annotated) made the function answer `introns: []`,
#     `removed_by_splicing: 0`, 5'UTR 1 nt — the one thing it exists to
#     report, silently invisible;
#   * a CDS ending in the `TGA TGA` tandem-stop idiom — used throughout the
#     reporter's own library — came back warned as "splicing may have shifted
#     the frame, or the annotation is stale";
#   * `tx_start` / `tx_end` bounded the unit while the response still named
#     the annotated promoter and terminator, reading as though those features
#     had defined it.


class TestPromoterSwallowingTheIntron:
    """A promoter feature drawn over the transcribed leader."""

    # Same construct, but the promoter annotation runs to the end of the
    # leader intron instead of stopping at the TSS.
    FEATS_SWALLOWED = [
        {"start": 0,  "end": P3, "strand": 1, "type": "promoter",
         "label": "Pdemo"},
        {"start": P2, "end": P3, "strand": 1, "type": "intron",
         "label": "leader intron"},
        {"start": P4, "end": P5, "strand": 1, "type": "CDS", "label": "GFP"},
        {"start": P5, "end": P6, "strand": 1, "type": "terminator",
         "label": "NOS"},
    ]

    def test_the_intron_really_does_vanish(self):
        """The premise: the TSS is the promoter's END, so an intron inside the
        promoter is upstream of the transcript and cannot be spliced. That is
        correct behaviour and it is why the warning has to exist."""
        r = sa._predict_transcript(SEQ, self.FEATS_SWALLOWED, circular=True)
        assert r["ok"], r
        assert r["introns"] == []
        assert r["removed_by_splicing"] == []

    def test_and_it_is_now_said_out_loud(self):
        r = sa._predict_transcript(SEQ, self.FEATS_SWALLOWED, circular=True)
        hits = [w for w in r["warnings"]
                if "INSIDE the chosen promoter" in w]
        assert len(hits) == 1, r["warnings"]
        w = hits[0]
        assert "leader intron" in w and "Pdemo" in w
        # and it names the way out, which is the only actionable part
        assert "tx_start" in w

    def test_an_exon_feature_inside_the_promoter_counts_too(self):
        feats = self.FEATS_SWALLOWED + [
            {"start": P1, "end": P2, "strand": 1, "type": "exon",
             "label": "leader exon 1"},
        ]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert any("exon" in w and "INSIDE the chosen promoter" in w
                   for w in r["warnings"]), r["warnings"]

    def test_an_intron_on_the_other_strand_is_not_warned_about(self):
        """It belongs to a gene transcribed the other way; warning about it
        would be noise, and the existing opposite-strand rule already says so
        for introns inside the unit."""
        feats = [dict(f) for f in self.FEATS_SWALLOWED]
        feats[1] = dict(feats[1], strand=-1)
        r = sa._predict_transcript(SEQ, feats, circular=True)
        assert not [w for w in r["warnings"]
                    if "INSIDE the chosen promoter" in w], r["warnings"]

    def test_the_normal_annotation_stays_silent(self):
        """The promoter that stops at the TSS must not pick up a warning."""
        r = _predict()
        assert not [w for w in r["warnings"]
                    if "INSIDE the chosen promoter" in w], r["warnings"]

    def test_an_intron_outside_the_unit_is_reported_too(self):
        """Any other reason an annotated intron misses the unit is the same
        class of silence: `introns: []` with no explanation."""
        feats = [f for f in FEATS if f["type"] != "intron"] + [
            {"start": P5 + 10, "end": P5 + 30, "strand": 1, "type": "intron",
             "label": "downstream intron"},
        ]
        # bound the unit short so the extra intron falls outside it
        r = sa._predict_transcript(SEQ, feats, circular=True,
                                   tx_start=P1, tx_end=P5)
        assert r["ok"], r
        assert any("lies OUTSIDE this transcription unit" in w
                   for w in r["warnings"]), r["warnings"]


class TestBoundedUnitDoesNotContradictItself:
    """Field report 2026-11-08. The promoter-collision warning fired on the
    OVERLAP alone, so a caller who did exactly what it told them — pass
    `tx_start` at the real TSS — got back `introns: 1, spliced: true,
    removed_by_splicing: N` together with warnings insisting the intron "lies
    upstream of the transcript and is NOT spliced out of this message", whose
    remedy was the thing they had already done. One response, two halves
    contradicting each other."""

    FEATS = TestPromoterSwallowingTheIntron.FEATS_SWALLOWED

    def _bounded(self):
        # tx_start at the REAL transcription start, inside the promoter
        # feature and upstream of the leader intron.
        return sa._predict_transcript(SEQ, self.FEATS, circular=True,
                                      tx_start=P1, tx_end=P6)

    def test_the_intron_is_actually_spliced_now(self):
        r = self._bounded()
        assert r["ok"], r
        assert len(r["introns"]) == 1, r["introns"]
        assert r["spliced"] is True

    def test_no_warning_claims_it_is_unspliced(self):
        r = self._bounded()
        bad = [w for w in r["warnings"]
               if "NOT spliced out of this message" in w]
        assert bad == [], (
            "the answer says the intron IS spliced; a warning saying it is "
            f"not contradicts it: {bad}")

    def test_the_promoter_collision_warning_is_gone(self):
        r = self._bounded()
        assert not [w for w in r["warnings"]
                    if "INSIDE the chosen promoter" in w], r["warnings"]

    def test_the_override_is_still_stated(self):
        """Silencing the contradiction must not silence the real notice that
        annotation did not define the unit."""
        r = self._bounded()
        assert any("did NOT define it" in w for w in r["warnings"]), \
            r["warnings"]

    def test_the_warning_survives_where_it_is_still_true(self):
        """Gated on the OUTCOME, not removed: with the unit bounded so the
        intron really does fall outside it, the diagnostic must still fire —
        and must not advise passing tx_start, which the caller just did."""
        r = sa._predict_transcript(SEQ, self.FEATS, circular=True,
                                   tx_start=P3, tx_end=P6)
        assert r["ok"], r
        assert r["introns"] == []
        hits = [w for w in r["warnings"] if "INSIDE the chosen promoter" in w]
        assert len(hits) == 1, r["warnings"]
        assert "outside this transcription unit" in hits[0]
        assert "pass tx_start at the real transcription start" not in hits[0]
        assert "widen tx_start/tx_end" in hits[0]

    def test_annotation_bounded_still_gets_the_tx_start_remedy(self):
        r = sa._predict_transcript(SEQ, self.FEATS, circular=True)
        hits = [w for w in r["warnings"] if "INSIDE the chosen promoter" in w]
        assert len(hits) == 1, r["warnings"]
        assert "tx_start" in hits[0]


class TestSpliceBlockCarriesItsDetectionPower:
    """A site count is unreadable without the sensitivity behind it, and a
    pair list is unreadable without knowing which ends are REAL boundaries."""

    def test_basis_is_reported_beside_the_counts(self):
        r = sa._predict_transcript(SEQ, FEATS, circular=True,
                                   check_splice=True)
        sb = r["splice"]
        assert sb["skipped"] is False, sb
        assert sb["basis"]["donor"]["sensitivity"] is not None
        assert sb["basis"]["fpr_target"] is not None

    def test_pairs_are_reported_and_carry_annotated_ends(self):
        r = sa._predict_transcript(SEQ, FEATS, circular=True,
                                   check_splice=True)
        sb = r["splice"]
        assert "n_pairs" in sb and "pairs" in sb
        assert "n_pairs_fully_cryptic" in sb
        for pair in sb["pairs"]:
            assert "annotated_ends" in pair
            assert "both_above_threshold" in pair

    def test_an_empty_scan_says_what_empty_means(self):
        """`n_cryptic: 0` on a false-positive-anchored model is a weaker claim
        than 'clean', and the response has to say so."""
        r = sa._predict_transcript(SEQ, FEATS, circular=True,
                                   check_splice=True)
        sb = r["splice"]
        if sb.get("n_cryptic") or sb.get("n_pairs_fully_cryptic"):
            pytest.skip("this fixture is not clean")
        assert any("weaker statement than a clean sequence" in w
                   for w in r["warnings"]), r["warnings"]


class TestTandemStopCodons:
    """`TGA TGA` is a terminator, not a stale annotation."""

    # CDS as above but with a second in-frame stop, and the annotation drawn
    # over both — exactly what the reporter's library does.
    CDS2 = CDS + "TGA"
    SEQ2 = PROM + LEAD1 + INTRON + LEAD2 + CDS2 + TERM
    Q4 = P4
    Q5 = Q4 + len(CDS2)
    Q6 = Q5 + len(TERM)
    FEATS2 = [
        {"start": 0,  "end": P1, "strand": 1, "type": "promoter", "label": "35S"},
        {"start": P2, "end": P3, "strand": 1, "type": "intron",
         "label": "leader intron"},
        {"start": Q4, "end": Q5, "strand": 1, "type": "CDS", "label": "GFP"},
        {"start": Q5, "end": Q6, "strand": 1, "type": "terminator",
         "label": "NOS"},
    ]

    def _run(self):
        r = sa._predict_transcript(self.SEQ2, self.FEATS2, circular=True)
        assert r["ok"], r
        return r

    def test_the_premise_two_in_frame_stops(self):
        assert self.CDS2.endswith("TAATGA")
        assert len(self.CDS2) % 3 == 0

    def test_no_stale_annotation_warning(self):
        r = self._run()
        assert not [w for w in r["warnings"] if "annotation is stale" in w], \
            r["warnings"]

    def test_the_count_is_reported_instead(self):
        r = self._run()
        assert r["cds"]["stop_codons"] == 2

    def test_a_single_stop_still_reports_one(self):
        r = _predict()
        assert r["cds"]["stop_codons"] == 1

    def test_a_genuinely_stale_annotation_is_still_warned_about(self):
        """The warning must keep working for what it was written for: trailing
        bases that are NOT stops."""
        feats = [dict(f) for f in self.FEATS2]
        feats[2] = dict(feats[2], end=self.Q5 + 9)   # 3 extra sense codons
        r = sa._predict_transcript(self.SEQ2, feats, circular=True)
        assert any("annotation is stale" in w for w in r["warnings"]), \
            r["warnings"]

    def test_three_tandem_stops_are_fine_too(self):
        cds3 = CDS + "TGATAG"
        seq3 = PROM + LEAD1 + INTRON + LEAD2 + cds3 + TERM
        q5 = P4 + len(cds3)
        feats = [
            {"start": 0,  "end": P1, "strand": 1, "type": "promoter", "label": "P"},
            {"start": P4, "end": q5, "strand": 1, "type": "CDS", "label": "GFP"},
            {"start": q5, "end": q5 + len(TERM), "strand": 1,
             "type": "terminator", "label": "T"},
        ]
        r = sa._predict_transcript(seq3, feats, circular=True)
        assert r["ok"], r
        assert r["cds"]["stop_codons"] == 3
        assert not [w for w in r["warnings"] if "annotation is stale" in w]


class TestExplicitBoundsAreAttributed:
    """`tx_start` / `tx_end` win, and the response says so."""

    def test_bounds_from_annotation_by_default(self):
        assert _predict()["unit"]["bounds_from"] == "annotation"

    def test_bounds_from_tx_when_given(self):
        r = sa._predict_transcript(SEQ, FEATS, circular=True,
                                   tx_start=P1, tx_end=P6)
        assert r["ok"], r
        assert r["unit"]["bounds_from"] == "tx_start/tx_end"

    def test_the_overridden_features_are_named(self):
        r = sa._predict_transcript(SEQ, FEATS, circular=True,
                                   tx_start=P1, tx_end=P6)
        assert any("did NOT define it" in w for w in r["warnings"]), \
            r["warnings"]

    def test_one_sided_bounds_still_refuse_with_a_reason(self):
        """Not a silent failure — the report suspected one, and the refusal is
        what proves it isn't."""
        r = sa._predict_transcript(SEQ, FEATS, circular=True, tx_start=P1)
        assert r["ok"] is False
        assert "tx_start and tx_end must BOTH be integers" in r["error"]

    def test_every_refusal_carries_an_error_string(self):
        """Audited because the report read `ok: False` as "no reason given"."""
        for kw in ({"tx_start": P1},
                   {"tx_end": P6},
                   {"tx_start": 10, "tx_end": 10},
                   {"tx_start": P6, "tx_end": P1},
                   {"tx_start": P5 + 5, "tx_end": P6}):
            r = sa._predict_transcript(SEQ, FEATS, circular=True, **kw)
            if r.get("ok"):
                continue
            assert r.get("error"), kw


class TestOutOfUnitWarningsAreBounded:
    """One line each is right for a plasmid and wrong for a chromosome
    record, where every gene's introns are "outside this unit"."""

    def _run(self, n_in, n_out):
        pad = "ACGT" * 200
        seq = SEQ + pad
        feats = [
            {"start": 0, "end": P3, "strand": 1, "type": "promoter",
             "label": "P"},
            {"start": P4, "end": P5, "strand": 1, "type": "CDS", "label": "C"},
            {"start": P5, "end": P6, "strand": 1, "type": "terminator",
             "label": "T"},
        ]
        for i in range(n_in):
            feats.append({"start": P2 + i, "end": P2 + i + 5, "strand": 1,
                          "type": "intron", "label": f"in{i}"})
        for i in range(n_out):
            feats.append({"start": P6 + 10 * i, "end": P6 + 10 * i + 5,
                          "strand": 1, "type": "intron", "label": f"out{i}"})
        r = sa._predict_transcript(seq, feats, circular=True)
        assert r["ok"], r
        return r["warnings"]

    def test_a_handful_each_get_their_own_line(self):
        w = self._run(3, 2)
        assert sum("INSIDE the chosen promoter" in x for x in w) == 3
        assert sum("lies OUTSIDE this" in x for x in w) == 2
        assert not [x for x in w if x.startswith("…and")]

    def test_past_the_cap_they_are_counted_not_dropped(self):
        w = self._run(30, 30)
        cap = sa._TX_MAX_OUT_OF_UNIT
        assert sum("INSIDE the chosen promoter" in x for x in w) == cap
        assert sum("lies OUTSIDE this" in x for x in w) == cap
        # the remainder is stated, so a truncated list never reads as the
        # whole list
        assert f"…and {30 - cap} more intron/exon feature(s)" in " ".join(w)
        assert f"…and {30 - cap} more annotated intron(s)" in " ".join(w)


class TestPairsThatEatTheCdsAreNeverSilent:
    """Field report 2026-11-08, round 3. The pair warning was gated on
    "neither end annotated", on the reasoning that one annotated end means
    the model missed the real partner. That was true of the case that
    prompted it and wrong as a rule: a REAL donor — one the spliceosome
    demonstrably uses — reaching a cryptic acceptor that OUTSCORES its own
    partner is exon skipping, and an annotated end makes it MORE credible,
    not less. On the reporting construct the strongest pair in the unit
    (combined 10.964, both ends over threshold, in-frame loss) would have
    spliced the entire CDS out while the response said nothing.

    Dedicated fixture: the annotated intron OPENS with a strong donor, and a
    strong cryptic acceptor sits inside the CDS 124 nt downstream — codon
    aligned as `TTT`x6 + `AGG` so it spells a real pyrimidine tract without
    a stop. `n_pairs_fully_cryptic` is 0 here, so the pre-fix gate made this
    silent."""

    PROM = "TTGACA" + "C" * 40 + "TATAAT"
    LEAD1 = "CAC" * 13                              # no ATG
    INTRON = "GTAAGT" + "CTTCTTCTTCTTCTTCTT" * 3 + "TGCAG"
    LEAD2 = "CCACCACCA"
    CDS = ("ATGGCTAAAGGTGAAGAATTATTCACTGGT"
           + "TTT" * 6 + "AGG"
           + "GTTGTCCCAATTTTGGTTGAATTAGATGGTGATGTTAAT" + "TAA")
    TERM = "AATAAA" + "G" * 30
    SEQ = PROM + LEAD1 + INTRON + LEAD2 + CDS + TERM
    Q1 = len(PROM)
    Q2 = Q1 + len(LEAD1)
    Q3 = Q2 + len(INTRON)
    Q4 = Q3 + len(LEAD2)
    Q5 = Q4 + len(CDS)
    Q6 = Q5 + len(TERM)
    FEATS = [
        {"start": 0,   "end": Q1, "strand": 1, "type": "promoter",   "label": "P"},
        {"start": Q2,  "end": Q3, "strand": 1, "type": "intron",     "label": "ivr"},
        {"start": Q4,  "end": Q5, "strand": 1, "type": "CDS",        "label": "gene"},
        {"start": Q5,  "end": Q6, "strand": 1, "type": "terminator", "label": "T"},
    ]

    def _res(self):
        r = sa._predict_transcript(self.SEQ, self.FEATS, circular=True,
                                   check_splice=True)
        assert r["ok"], r
        return r

    def test_the_fixture_poses_the_case_the_old_gate_missed(self):
        """The premise. If `fully_cryptic` ever becomes non-zero here the
        fixture has drifted and the rest of the class proves less."""
        sb = self._res()["splice"]
        assert sb["n_pairs_excising_cds"] >= 1
        assert sb["n_pairs_fully_cryptic"] == 0
        into = sb["pairs_excising_cds"]
        assert any(p["annotated_ends"] == ["donor"] for p in into), into

    def test_an_annotated_end_does_not_buy_silence(self):
        r = self._res()
        assert any("would excise part of the CDS" in w
                   for w in r["warnings"]), r["warnings"]
        assert any("MORE credible" in w for w in r["warnings"])

    def test_the_warning_names_the_worst_pair(self):
        r = self._res()
        worst = max(r["splice"]["pairs_excising_cds"],
                    key=lambda p: p["combined_score"])
        assert any(str(worst["combined_score"]) in w
                   for w in r["warnings"]), r["warnings"]

    def test_every_listed_pair_is_a_plausible_intron(self):
        sb = self._res()["splice"]
        assert sb["pairs_excising_cds"], "fixture must list one"
        for pair in sb["pairs_excising_cds"]:
            d = pair["donor"]["position"]
            a = pair["acceptor"]["position"] + 2
            assert max(d, a) > min(d, a)
            assert 60 <= pair["intron_len"] <= 3000

    def test_a_clean_unit_still_reports_the_all_clear_caveat(self):
        """The empty-scan caveat must account for the new class too, not
        fire while a CDS-excising pair is listed."""
        r = self._res()
        sb = r["splice"]
        if sb["n_pairs_excising_cds"] or sb["n_cryptic"]:
            assert not any("weaker statement than a clean sequence" in w
                           for w in r["warnings"])


class TestUnitCoordinatesStateTheirBasis:
    """Third coordinate-basis misread in three rounds. `{start: 2070, end:
    4384}` reads as 1-based inclusive; it is 0-based half-open, and a caller
    passing 2070 where 2069 was meant silently clipped a base off exon 1."""

    def test_genomic_span_is_labelled_and_given_both_ways(self):
        r = _predict()
        g = r["unit"]["genomic"]
        assert g["basis"] == "0-based half-open"
        assert g["start_1based"] == g["start"] + 1
        assert g["end_1based"] == g["end"]

    def test_the_0_based_values_are_what_tx_start_takes(self):
        """The echoed `start` must be exactly what you passed, so the two
        cannot drift apart."""
        r = sa._predict_transcript(SEQ, FEATS, circular=True,
                                   tx_start=P1, tx_end=P5)
        assert r["unit"]["genomic"]["start"] == P1
        assert r["unit"]["genomic"]["end"] == P5


class TestExonOutsideTheUnitIsMissingNotRetained:
    """One template served both feature types and the consequence was
    backwards for an exon: an exon left outside the unit is not "retained
    instead of spliced", it is absent from the message."""

    def test_an_intron_says_not_spliced_out(self):
        r = sa._predict_transcript(
            SEQ, TestPromoterSwallowingTheIntron.FEATS_SWALLOWED,
            circular=True)
        hits = [w for w in r["warnings"] if "INSIDE the chosen promoter" in w]
        assert hits and "NOT spliced out of this message" in hits[0]

    def test_an_exon_says_missing_from_the_message(self):
        feats = TestPromoterSwallowingTheIntron.FEATS_SWALLOWED + [
            {"start": P1, "end": P2, "strand": 1, "type": "exon",
             "label": "leader exon 1"},
        ]
        r = sa._predict_transcript(SEQ, feats, circular=True)
        exon_w = [w for w in r["warnings"]
                  if "INSIDE the chosen promoter" in w and "exon" in w]
        assert exon_w, r["warnings"]
        assert "MISSING from this message" in exon_w[0]
        assert "NOT spliced out" not in exon_w[0]
