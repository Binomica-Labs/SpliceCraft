"""Eukaryotic transcription-unit assembler (`splicecraft_cassette`).

The plant-nuclear counterpart to `test_rna_fold.py::TestAssembleOperon`. The
biological invariant that matters most here is PROTEIN PRESERVATION: every
motif-scrub pass is a synonymous codon swap, so a scrub that changes the
translation is a silent-biology failure of exactly the class [INV-109] and the
codon-optimizer property tests exist to catch. Cross-validated against
Biopython's translator rather than the module's own.
"""
import pytest
from Bio.Seq import Seq

import splicecraft_cassette as cas
from splicecraft_dataaccess import _CODON_BUILTIN_K12


# A promoter/terminator pair with no coding meaning — non-ACGT-free, and short
# enough to keep the scrub passes fast.
PROM = "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGCACCTGCA" * 2
TERM = "GGCGCGCCTTAATTAAGGCCGGCCTTCGAAGGATCCGTCGAC" * 2
# Signal peptide + mature CDS, both codon-aligned.
SP = "ATGGCTAGCTTGCTGCTTCTCGCCCTGCTGGTGCTGAGCCTG"          # 42 nt = 14 codons
CDS = ("GCTACATTTGAAATCGTTAACCGTTGCTCTTATACAGTTTGGGCTGCTGCT"
       "AGCAAAGGTGATGCTGCTCTGGATGCTGGTGGTCGTCAGCTGAACAGC")      # 98 -> pad below
CDS = CDS[:96]                                                  # 96 nt = 32 codons
KDEL = "AAGGATGAACTG"                                           # 12 nt = 4 codons


def _translate(dna: str) -> str:
    """Independent reference translation (Biopython, not the module's hook)."""
    return str(Seq(dna).translate())


class TestLayoutContract:
    """The tiling guarantee `_assemble_operon` also makes: layout covers the
    sequence exactly, no gap, no overlap."""

    def test_layout_tiles_sequence_exactly(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, utr5="AGGAGACAACAACAACAAAA", signal_peptide=SP,
            cds=CDS, ctag=KDEL, terminator=TERM, scrub_motifs=False)
        seq = r["sequence"]
        cursor = 0
        for el in r["layout"]:
            assert el["start"] == cursor            # contiguous
            assert seq[el["start"]:el["end"]]       # non-empty
            cursor = el["end"]
        assert cursor == len(seq)                   # covers the whole sequence
        assert [el["kind"] for el in r["layout"]] == [
            "promoter", "utr5", "signal_peptide", "cds", "ctag", "terminator"]

    def test_each_part_lands_verbatim_at_its_slot(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, cds=CDS, terminator=TERM, scrub_motifs=False)
        seq = r["sequence"]
        by_kind = {el["kind"]: el for el in r["layout"]}
        assert seq[by_kind["promoter"]["start"]:by_kind["promoter"]["end"]] == PROM
        assert seq[by_kind["cds"]["start"]:by_kind["cds"]["end"]] == CDS
        assert seq[by_kind["terminator"]["start"]:by_kind["terminator"]["end"]] == TERM

    def test_element_order_is_fixed_not_call_order(self):
        # Parts supplied out of order still assemble 5'->3' by grammar.
        r = cas._assemble_expression_cassette(
            terminator=TERM, cds=CDS, promoter=PROM, scrub_motifs=False)
        assert [el["kind"] for el in r["layout"]] == ["promoter", "cds", "terminator"]

    def test_empty_cassette_raises(self):
        with pytest.raises(ValueError, match="empty"):
            cas._assemble_expression_cassette(promoter=None, cds=None)

    def test_ambiguity_codes_rejected(self):
        # An N in a part would translate to X and silently corrupt the scrub.
        with pytest.raises(ValueError, match="unambiguous"):
            cas._assemble_expression_cassette(promoter="ACGTNNNN", cds=CDS)


class TestFrame:
    """Frame across the signal-peptide -> CDS -> C-tag ORF. A part whose
    length is not a multiple of 3 frameshifts everything 3' of it."""

    def test_aligned_parts_are_in_frame(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, signal_peptide=SP, cds=CDS, ctag=KDEL,
            terminator=TERM, scrub_motifs=False)
        assert all(f["in_frame"] for f in r["frame"])
        assert [f["element"] for f in r["frame"]] == \
            ["signal_peptide", "cds", "ctag"]
        assert r["frame"][-1]["offset_out"] == 0

    def test_offset_is_cumulative_across_junctions(self):
        # A 20-nt signal peptide (not %3) must shift BOTH downstream elements.
        r = cas._assemble_expression_cassette(
            promoter=PROM, signal_peptide=SP[:20], cds=CDS, ctag=KDEL,
            terminator=TERM, scrub_motifs=False)
        frame = {f["element"]: f for f in r["frame"]}
        assert frame["signal_peptide"]["offset_in"] == 0
        assert frame["signal_peptide"]["codon_aligned"] is False
        assert frame["cds"]["offset_in"] == 20 % 3        # == 2, shifted
        assert frame["cds"]["in_frame"] is False          # aligned but wrong offset
        assert frame["cds"]["codon_aligned"] is True

    def test_scrubbing_out_of_frame_block_raises(self):
        # Refuse to hand a frameshifted block to the codon layer: it would
        # translate in the wrong frame and "preserve" a protein that was
        # never real.
        with pytest.raises(ValueError, match="out-of-frame"):
            cas._assemble_expression_cassette(
                promoter=PROM, signal_peptide=SP[:20], cds=CDS,
                terminator=TERM, codon_table=_CODON_BUILTIN_K12,
                scrub_motifs=True)


class TestMotifScan:
    """Regression guards for two silent-all-clear bugs found 2026-07-25."""

    def test_dict_of_named_sites_is_scanned_not_its_keys(self):
        # Regression guard for 2026-07-25 fix: `_forbidden_hit_set` takes an
        # iterable of SITE strings. Passing a {name: site} dict iterates the
        # KEYS, which aren't valid IUPAC, so every motif is skipped and a
        # cassette full of BsaI sites scans clean.
        hits = cas._cassette_scan_motifs("AAAAGGTCTCAAAA", {"BsaI": "GGTCTC"})
        assert [h["motif"] for h in hits] == ["BsaI"]
        assert hits[0]["position"] == 4

    def test_reverse_strand_sites_are_found(self):
        # Regression guard for 2026-07-25 fix: a Type IIS site is
        # non-palindromic and cuts from the bottom strand too. GAGACC is the
        # reverse complement of BsaI's GGTCTC.
        hits = cas._cassette_scan_motifs("AAAAGAGACCAAAA", {"BsaI": "GGTCTC"})
        assert [h["motif"] for h in hits] == ["BsaI"]
        assert hits[0]["strand"] == "-"

    def test_palindromic_site_not_double_counted(self):
        # A palindrome is its own reverse complement; report it once.
        hits = cas._cassette_scan_motifs("AAAGGATCCAAA", {"BamHI": "GGATCC"})
        assert len(hits) == 1

    def test_motifs_outside_coding_are_flagged_unfixable(self):
        # AATAAA inside a terminator cannot be synonymously repaired.
        r = cas._assemble_expression_cassette(
            promoter=PROM, cds=CDS, terminator="GGATCCAATAAAGGTACCTT",
            scrub_motifs=False)
        polya = [h for h in r["motif_hits"] if h["motif"] == "polyA_AATAAA"]
        assert polya, "AATAAA in the terminator was not detected"
        assert all(h["fixable_synonymously"] is False for h in polya)
        assert any("outside the coding block" in w for w in r["warnings"])


class TestScrubPreservesProtein:
    """The mission-critical invariant: scrubbing is synonymous, always."""

    def test_protein_identical_after_scrub(self):
        before = _translate(SP + CDS + KDEL)
        r = cas._assemble_expression_cassette(
            promoter=PROM, signal_peptide=SP, cds=CDS, ctag=KDEL,
            terminator=TERM, codon_table=_CODON_BUILTIN_K12,
            scrub_motifs=True)
        coding = "".join(
            r["sequence"][el["start"]:el["end"]] for el in r["layout"]
            if el["kind"] in ("signal_peptide", "cds", "ctag"))
        assert _translate(coding) == before        # cross-validated, Biopython
        assert r["scrub"]["skipped"] is False

    def test_scrub_preserves_element_lengths_so_layout_stays_valid(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, signal_peptide=SP, cds=CDS, ctag=KDEL,
            terminator=TERM, codon_table=_CODON_BUILTIN_K12,
            scrub_motifs=True)
        by_kind = {el["kind"]: el for el in r["layout"]}
        assert by_kind["signal_peptide"]["end"] - by_kind["signal_peptide"]["start"] == len(SP)
        assert by_kind["cds"]["end"] - by_kind["cds"]["start"] == len(CDS)
        assert by_kind["ctag"]["end"] - by_kind["ctag"]["start"] == len(KDEL)
        cursor = 0                                  # still tiles exactly
        for el in r["layout"]:
            assert el["start"] == cursor
            cursor = el["end"]
        assert cursor == len(r["sequence"])

    def test_enzyme_site_planted_in_coding_is_removed(self):
        # GGTCTC placed in-frame and synonymously escapable.
        planted = "ATGGGTCTCGCTACATTTGAAATCGTTAACCGTTGCTCTTATACAGTT"
        assert len(planted) % 3 == 0
        r = cas._assemble_expression_cassette(
            promoter=PROM, cds=planted, terminator=TERM,
            codon_table=_CODON_BUILTIN_K12, scrub_motifs=True, enzyme="BsaI")
        coding = next(r["sequence"][el["start"]:el["end"]]
                      for el in r["layout"] if el["kind"] == "cds")
        assert _translate(coding) == _translate(planted)     # protein held
        assert "GGTCTC" not in coding and "GAGACC" not in coding

    def test_scrub_skipped_without_codon_table_and_warns(self):
        # Never silently rewrite against a guessed frame.
        r = cas._assemble_expression_cassette(
            promoter=PROM, cds=CDS, terminator=TERM, scrub_motifs=True)
        assert r["scrub"]["skipped"] is True
        assert "codon_table" in r["scrub"]["reason"]
        assert any("codon_table" in w for w in r["warnings"])

    def test_unknown_enzyme_raises(self):
        with pytest.raises(ValueError, match="unknown enzyme"):
            cas._assemble_expression_cassette(
                promoter=PROM, cds=CDS, terminator=TERM,
                codon_table=_CODON_BUILTIN_K12, enzyme="NotAnEnzyme")

    def test_internal_stop_codon_refused(self):
        stopped = "ATGGCTTAAGCTACATTTGAA"          # TAA at codon 3
        with pytest.raises(ValueError, match="internal stop"):
            cas._assemble_expression_cassette(
                promoter=PROM, cds=stopped, terminator=TERM,
                codon_table=_CODON_BUILTIN_K12, scrub_motifs=True)


class TestInternalRepeats:
    """A clonal host cannot segregate a silenced insert away, so an internal
    repeat is permanent rather than an F2 inconvenience."""

    UNIT = ("ACGTTGCAGGCCTTAAGGCATCGATCGGATCCTTAAGGCCTAGCTAGCTAGGCATCGGAT"
            "CCAAGGTTCCAAGGTTCCAAGCCTTAAGGCCTTAAGGCCTA")

    def test_direct_repeat_detected(self):
        dup = self.UNIT + "TTTTGGGGCCCCAAAA" + self.UNIT
        hits = cas._cassette_internal_repeats(dup, 100)
        direct = [h for h in hits if h["kind"] == "direct"]
        assert direct and direct[0]["length"] >= 100

    def test_reports_maximal_run_not_every_window(self):
        # A 101-bp repeat must be ONE finding, not 2 overlapping windows.
        dup = self.UNIT + "TTTTGGGGCCCCAAAA" + self.UNIT
        direct = [h for h in cas._cassette_internal_repeats(dup, 100)
                  if h["kind"] == "direct"]
        assert len(direct) == 1

    def test_inverted_repeat_detected(self):
        from splicecraft_biology import _rc
        inv = self.UNIT + "TTTTGGGGCCCCAAAA" + _rc(self.UNIT)
        hits = cas._cassette_internal_repeats(inv, 100)
        assert any(h["kind"] == "inverted" for h in hits)

    def test_excluded_mar_pair_suppressed(self):
        # The intentional flanking MAR pair must not be reported...
        dup = self.UNIT + "TTTTGGGGCCCCAAAA" + self.UNIT
        spans = [(0, len(self.UNIT)), (len(dup) - len(self.UNIT), len(dup))]
        assert not [h for h in cas._cassette_internal_repeats(dup, 100, exclude=spans)
                    if h["kind"] == "direct"]

    def test_exclusion_does_not_mask_a_third_copy(self):
        # ...but a MAR that also repeats against the payload still is.
        dup = self.UNIT + "TTTTGGGG" + self.UNIT + "CCCCAAAA" + self.UNIT
        spans = [(0, len(self.UNIT)), (len(dup) - len(self.UNIT), len(dup))]
        assert [h for h in cas._cassette_internal_repeats(dup, 100, exclude=spans)
                if h["kind"] == "direct"]

    def test_no_false_positive_on_unique_sequence(self):
        # Seeded PRNG, not an arithmetic pattern: `"ACGT"[(i * 7 + i // 4) % 4]`
        # looks scrambled but has period 16, and the detector correctly found a
        # 384-bp repeat in it.
        import random
        rng = random.Random(20260725)
        uniq = "".join(rng.choice("ACGT") for _ in range(4000))
        assert not [h for h in cas._cassette_internal_repeats(uniq, 100)
                    if h["kind"] == "direct"]


class TestOverhangs:
    """Type IIS fusion sites: a skipped position is a build instruction, a
    duplicated fusion site is a defect."""

    def test_full_stack_chains_on_grammar_defaults(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, utr5="AGGAGACAACAACAACAAAA", signal_peptide=SP,
            cds=CDS, ctag=KDEL, terminator=TERM, scrub_motifs=False)
        assert all(c["nominal_match"] for c in r["overhangs"]["chain"])
        assert not r["overhangs"]["duplicates"]

    def test_shared_slot_elements_are_one_fragment(self):
        # signal_peptide + cds both occupy the CDS position -> one junction,
        # not two, and no bogus duplicate on their shared AGGT.
        r = cas._assemble_expression_cassette(
            promoter=PROM, signal_peptide=SP, cds=CDS, terminator=TERM,
            scrub_motifs=False)
        froms = [c["from"] for c in r["overhangs"]["chain"]]
        assert "signal_peptide+cds" in froms
        assert not r["overhangs"]["duplicates"]

    def test_skipped_ctag_reports_required_fusion_site(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, cds=CDS, terminator=TERM, scrub_motifs=False)
        last = r["overhangs"]["chain"][-1]
        assert last["required"] == "GGTA"        # terminator's oh5
        assert last["nominal"] == "GCTT"         # CDS position's default oh3
        assert last["nominal_match"] is False
        assert any("domesticate" in w for w in r["warnings"])


class TestGenBankEmission:
    """Features must be correctly TYPED, not misc_feature with a label -- that
    is what makes the file render in SnapGene and survive a round-trip."""

    def _record(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, utr5="AGGAGACAACAACAACAAAA", signal_peptide=SP,
            cds=CDS, ctag=KDEL, terminator=TERM, scrub_motifs=False)
        return r, cas._cassette_to_record(r, "THM-TEST", description="unit test")

    def test_feature_types_are_insdc(self):
        _, rec = self._record()
        types = {f.type for f in rec.features}
        assert {"promoter", "5'UTR", "sig_peptide", "CDS", "terminator"} <= types
        assert "misc_feature" not in types

    def test_coordinates_match_layout(self):
        r, rec = self._record()
        by_label = {f.qualifiers["label"][0]: f for f in rec.features}
        for el in r["layout"]:
            feat = by_label[el["name"]]
            assert int(feat.location.start) == el["start"]
            assert int(feat.location.end) == el["end"]

    def test_whole_orf_feature_spans_sp_through_ctag(self):
        r, rec = self._record()
        orf = [f for f in rec.features if f.qualifiers["label"][0].endswith("ORF")]
        assert len(orf) == 1
        starts = [el["start"] for el in r["layout"]
                  if el["kind"] == "signal_peptide"]
        ends = [el["end"] for el in r["layout"] if el["kind"] == "ctag"]
        assert int(orf[0].location.start) == starts[0]
        assert int(orf[0].location.end) == ends[0]

    def test_orf_translates_without_internal_stop(self):
        r, rec = self._record()
        orf = next(f for f in rec.features
                   if f.qualifiers["label"][0].endswith("ORF"))
        protein = _translate(str(rec.seq[orf.location.start:orf.location.end]))
        assert "*" not in protein.rstrip("*")

    def test_notes_ride_through_to_the_feature(self):
        r, _ = self._record()
        rec = cas._cassette_to_record(
            r, "THM-03", notes={"ctag": "KDEL is NOT cleaved -- not FEMA GRAS 3814"})
        ctag = next(f for f in rec.features
                    if f.qualifiers["label"][0] == "ctag")
        assert "FEMA GRAS" in ctag.qualifiers["note"][0]

    def test_round_trips_through_genbank_text(self):
        from io import StringIO
        from Bio import SeqIO
        _, rec = self._record()
        handle = StringIO()
        SeqIO.write(rec, handle, "genbank")
        handle.seek(0)
        back = SeqIO.read(handle, "genbank")
        assert str(back.seq) == str(rec.seq)
        assert {f.type for f in back.features} == {f.type for f in rec.features}


class TestSpliceIntegration:
    """§1.4.4 — cryptic splice sites are scored over the assembled cassette."""

    STRONG_DONOR = "CAGGTAAGT"
    STRONG_ACCEPTOR = "T" * 18 + "AGG"

    def test_clean_cassette_reports_a_ran_check(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, cds=CDS, terminator=TERM, scrub_motifs=False)
        assert r["splice"]["skipped"] is False
        assert "sites" in r["splice"]

    def test_planted_donor_acceptor_pair_is_flagged(self):
        # A pair at plausible intron spacing can excise sequence from the
        # transcript -- the failure mode that leaves the DNA looking perfect.
        payload = (self.STRONG_DONOR + "ACGT" * 30 + self.STRONG_ACCEPTOR)
        r = cas._assemble_expression_cassette(
            promoter=PROM, cds=payload, terminator=TERM, scrub_motifs=False)
        assert r["splice"]["pairs"], "donor/acceptor pair not detected"
        assert any("PAIR" in w for w in r["warnings"])

    def test_check_can_be_disabled(self):
        r = cas._assemble_expression_cassette(
            promoter=PROM, cds=CDS, terminator=TERM, scrub_motifs=False,
            check_splice=False)
        assert r["splice"]["skipped"] is True

    def test_absent_model_warns_rather_than_reporting_clean(self):
        import splicecraft_splice as sp
        saved = (sp._SPLICE_DONOR_PWM, sp._SPLICE_ACCEPTOR_PWM,
                 sp._SPLICE_MODEL_META)
        try:
            sp._splice_load_model([], [], {})
            r = cas._assemble_expression_cassette(
                promoter=PROM, cds=CDS, terminator=TERM, scrub_motifs=False)
            assert r["splice"]["skipped"] is True
            assert any("SKIPPED" in w for w in r["warnings"])
        finally:
            sp._splice_load_model(*saved)


class TestDryRun:
    def test_dry_run_flag_is_reported_and_changes_nothing(self):
        kw = dict(promoter=PROM, cds=CDS, terminator=TERM, scrub_motifs=False)
        wet = cas._assemble_expression_cassette(**kw)
        dry = cas._assemble_expression_cassette(**kw, dry_run=True)
        assert dry["dry_run"] is True and wet["dry_run"] is False
        assert dry["sequence"] == wet["sequence"]
        assert dry["layout"] == wet["layout"]
