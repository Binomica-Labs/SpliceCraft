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

    def test_it_never_claims_to_be_a_trained_matrix(self):
        """The splice model IS calibrated; this is not, and a caller must be
        able to tell them apart."""
        assert "not a trained PWM" in sa._kozak_context("A" * 10, 6)["model"]

    def test_a_short_leader_is_marked_truncated_not_padded(self):
        k = sa._kozak_context("AT" + "ATGG", 2)
        assert k["truncated"] is True and len(k["context"]) < 10


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
