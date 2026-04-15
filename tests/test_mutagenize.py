"""
test_mutagenize — SOE-PCR mutagenesis primer design.

Covers the module-level helpers ported from mutagenesis_primers.py:
_mut_parse, _mut_translate, _mut_revcomp, _mut_design_outer,
_mut_design_inner, _mut_extract_cds, and the edge-case branch for
mutations within 60 nt of either CDS end.
"""
import pytest

import splicecraft as sc


# ── Hand-verifiable fixtures ──────────────────────────────────────────────────

# 246-nt realistic CDS (~50 % GC, varied codons) — 81 aa protein + stop.
# Hand-verified aa positions used in tests: 2=A, 3=E, 40=V, 78=A.
_CDS_LONG = (
    "ATG"
    "GCTGAAGTTCAGGATAACCTGGCGAAAGTTCAGGAAGCGGTTGATACCCTGAAACGTGGT"
    "CTGGAAGCGGCGAAAGCGACCCTGGAAAAAGCGGGTGAAGATATCGCGAAAGCGGTTGAT"
    "GGTAAACGTAAAGGCGATCTGGAAAAACTGGCGGAAGCGCTGCAGAAAGTTGAAGCGGAT"
    "ATCGCGAAAGCGGTTGATGGTAAACGTAAAGGCGATCTGGAAAAACTGGCGGAAGCGCTG"
    "TAA"
)


# ── _mut_parse ────────────────────────────────────────────────────────────────

class TestParseMutation:
    def test_basic(self):
        assert sc._mut_parse("W140F") == ("W", 140, "F")

    def test_lowercase_accepted(self):
        assert sc._mut_parse("w140f") == ("W", 140, "F")

    def test_stop_codon(self):
        assert sc._mut_parse("W140*") == ("W", 140, "*")

    def test_rejects_bad_format(self):
        with pytest.raises(ValueError):
            sc._mut_parse("W140")
        with pytest.raises(ValueError):
            sc._mut_parse("140F")
        with pytest.raises(ValueError):
            sc._mut_parse("not a mutation")


# ── _mut_revcomp ──────────────────────────────────────────────────────────────

class TestRevcomp:
    def test_simple(self):
        assert sc._mut_revcomp("ATGC") == "GCAT"

    def test_involutive(self):
        seq = "ATGCATGCGGTTAA"
        assert sc._mut_revcomp(sc._mut_revcomp(seq)) == seq


# ── _mut_translate ────────────────────────────────────────────────────────────

class TestTranslate:
    def test_stops_at_stop_codon(self):
        # "ATG AAA TAA GGG" — stop at codon 3, GGG never translated.
        assert sc._mut_translate("ATGAAATAAGGG") == "MK"

    def test_long_cds(self):
        protein = sc._mut_translate(_CDS_LONG)
        assert len(protein) == 81
        assert protein[0]  == "M"
        assert protein[1]  == "A"   # GCT
        assert protein[2]  == "E"   # GAA


# ── Outer primer design — BsaI tails (regression guard) ───────────────────────

class TestDesignOuter:
    """The outer primers are constant per CDS and must carry the BsaI-AATG
    (FWD) / BsaI-AACG (REV) tails that produce GB B3 / B5 overhangs after
    digestion. Changes here will break every Golden Braid L0 assembly."""

    def test_fwd_carries_bsai_aatg_tail(self):
        outer = sc._mut_design_outer(_CDS_LONG)
        assert outer["fwd"]["full"].startswith("CCCCGGTCTCAAATG")
        assert outer["b3_overhang"] == "AATG"

    def test_rev_carries_bsai_aacg_tail(self):
        outer = sc._mut_design_outer(_CDS_LONG)
        assert outer["rev"]["full"].startswith("CCCCGGTCTCAAACG")
        # The vector-side overhang name is CGTT (= revcomp of AACG on insert)
        assert outer["b5_overhang"] == "CGTT"

    def test_fwd_anneal_starts_after_atg(self):
        """FWD_outer anneal must begin at nt 4 (index 3) of the CDS so the
        AATG overhang reconstitutes the start codon in the assembled part."""
        outer = sc._mut_design_outer(_CDS_LONG)
        assert outer["fwd_anneal_start"] == 3
        # Anneal region is taken from _CDS_LONG[3:3+length]
        anneal = outer["fwd"]["anneal"]
        assert _CDS_LONG[3:3 + len(anneal)] == anneal

    def test_rev_anneal_is_revcomp_of_cds_end(self):
        outer = sc._mut_design_outer(_CDS_LONG)
        anneal = outer["rev"]["anneal"]
        end_rc = sc._mut_revcomp(_CDS_LONG)
        assert end_rc.startswith(anneal)


# ── Inner pair — revcomp invariant, WT codon check ────────────────────────────

class TestDesignInner:
    def test_rev_is_revcomp_of_fwd(self):
        """Inner REV must be the exact revcomp of inner FWD — this is the
        whole point of the SOE joint primer pair.
        Signature: _mut_design_inner(dna, mut_pos_1, mut_aa, wt_aa)."""
        inner = sc._mut_design_inner(_CDS_LONG, 40, "F", "V")  # V40F
        best = inner["candidates"][0]
        assert sc._mut_revcomp(best["fwd"]) == best["rev"]

    def test_wt_codon_mismatch_raises(self):
        """If the caller says WT='W' but the DNA codon at that position
        doesn't actually encode W, we must error rather than produce a
        nonsense mutation primer."""
        # Position 2 is A (GCT); caller claims WT='W' → error.
        with pytest.raises(ValueError, match="mutation says WT='W'"):
            sc._mut_design_inner(_CDS_LONG, 2, "F", "W")

    def test_mut_codon_differs_from_wt(self):
        """mut_codon must encode the requested mutant aa and differ from the
        wt codon so the DNA actually changes."""
        inner = sc._mut_design_inner(_CDS_LONG, 40, "F", "V")  # V40F
        assert inner["wt_codon"] == "GTT"                 # codon 40 of CDS
        assert sc._MUT_CODON_TO_AA[inner["mut_codon"]] == "F"
        assert inner["mut_codon"] != inner["wt_codon"]

    def test_mutation_string_format(self):
        """Mutation string format is WT_AA + pos + MUT_AA."""
        inner = sc._mut_design_inner(_CDS_LONG, 40, "F", "V")  # mut=F, wt=V
        assert inner["mutation"] == "V40F"

    def test_no_alt_codon_for_single_codon_aa(self):
        """Met has only one codon (ATG). Asking to mutate an interior Met
        back to Met must error — there is no alternative codon."""
        cds = "ATG" + ("GCT" * 30) + "ATG" + ("GCT" * 30) + "TAA"
        assert sc._mut_translate(cds)[31] == "M"
        with pytest.raises(ValueError):
            sc._mut_design_inner(cds, 32, "M", "M")


# ── Edge-case branch ──────────────────────────────────────────────────────────

class TestEdgeCase:
    """Mutations within _MUT_MIN_SOE_FRAG (60 nt) of either CDS end must
    trigger the modified-outer branch and skip the inner pair."""

    def test_near_start_triggers_modified_fwd(self):
        # Position 3 is E (codon 3 = GAA). Fragment A ≈ 9 nt → far below 60.
        inner = sc._mut_design_inner(_CDS_LONG, 3, "F", "E")   # E3F
        ec = inner["edge_case"]
        assert ec is not None
        assert ec["near_start"] is True
        assert ec["near_end"] is False
        assert ec["modified_outer"]["label"] == "modified_FWD_outer"
        # The modified FWD carries the BsaI-AATG tail like the normal FWD.
        assert ec["modified_outer"]["full"].startswith("CCCCGGTCTCAAATG")

    def test_near_end_triggers_modified_rev(self):
        # Position 78 is A (codon 78 = GCG). Fragment B ≈ 12 nt → below 60.
        inner = sc._mut_design_inner(_CDS_LONG, 78, "F", "A")  # A78F
        ec = inner["edge_case"]
        assert ec is not None
        assert ec["near_end"] is True
        assert ec["modified_outer"]["label"] == "modified_REV_outer"
        assert ec["modified_outer"]["full"].startswith("CCCCGGTCTCAAACG")

    def test_middle_mutation_no_edge_case(self):
        # Position 40 is V, well away from both ends (~120 nt from either).
        inner = sc._mut_design_inner(_CDS_LONG, 40, "F", "V")  # V40F
        assert inner["edge_case"] is None


# ── CDS extraction — strand and wrap handling ─────────────────────────────────

class TestExtractCds:
    """_mut_extract_cds must return the CDS in its biological 5'→3'
    orientation regardless of strand or origin-wrap."""

    def test_forward_strand_simple(self):
        seq = "AAAA" + _CDS_LONG + "TTTT"
        cds = sc._mut_extract_cds(seq, 4, 4 + len(_CDS_LONG), 1)
        assert cds == _CDS_LONG

    def test_reverse_strand_is_revcomp(self):
        """A CDS on the reverse strand at plasmid[a:b] should come back as
        revcomp(plasmid[a:b]) so the first codon is ATG."""
        rc = sc._mut_revcomp(_CDS_LONG)
        seq = "AAAA" + rc + "TTTT"
        cds = sc._mut_extract_cds(seq, 4, 4 + len(rc), -1)
        assert cds == _CDS_LONG
        assert cds.startswith("ATG")

    def test_wrap_around_origin(self):
        """A feature with end < start spans the origin. The extracted CDS
        must be tail + head, in order."""
        # Build a "plasmid" where the CDS starts near the end and wraps:
        # CDS = ATG + GCT*5 + TAA = 21 nt. Place first 15 nt at the end of
        # the plasmid and the last 6 nt at the start.
        cds = "ATG" + ("GCT" * 5) + "TAA"
        assert len(cds) == 21
        padding = "N" * 30
        # plasmid layout: [last 6 nt of cds][padding][first 15 nt of cds]
        plasmid = cds[15:] + padding + cds[:15]
        start = len(plasmid) - 15  # where the CDS head lives
        end   = 6                  # where the CDS tail ends (wrapped)
        extracted = sc._mut_extract_cds(plasmid, start, end, 1)
        assert extracted == cds


# ── CDS preview renderer  (_mut_build_preview_text) ──────────────────────────

class TestPreviewText:
    """The Mutagenize modal's preview widget — DNA row with per-codon AA
    letters centered beneath. The rendered Text should line up so that
    every AA sits directly under the middle base of its codon (column
    offset `codon_idx * 3 + 1` relative to the DNA chunk).
    """

    def _plain_lines(self, text) -> list[str]:
        """Return the rendered Text as stripped-of-styles lines."""
        return text.plain.split("\n")

    def test_empty_inputs_render_nothing(self):
        t = sc._mut_build_preview_text("", protein_override="")
        assert t.plain == ""

    def test_aa_only_when_no_dna(self):
        """Protein-input source before harmonization: AA wraps, no DNA row."""
        aa = "MALAK" * 4      # 20 aa
        t = sc._mut_build_preview_text("", protein_override=aa, line_width=12)
        lines = [l for l in self._plain_lines(t) if l]
        # Wraps at 12 → 20 aa → 2 lines of 12 + 8
        assert lines == [aa[:12], aa[12:]]

    def test_dna_and_aa_alignment(self):
        """For a 9-bp CDS 'ATGGCCAGC' with translation 'MAS', the AA row
        must place M/A/S at columns 1, 4, 7 of the DNA chunk — directly
        under the middle base of each codon."""
        cds = "ATGGCCAGC"
        t = sc._mut_build_preview_text(cds, line_width=9)
        lines = self._plain_lines(t)
        # Two rendered lines (DNA, AA) + terminating blank.
        # Drop the line-number prefix (" 1  " = 4 chars for a 9-digit width
        # of 1 + 2 spaces) — num_w = len("9") = 1, so prefix width = 3.
        num_w = 1
        pad = num_w + 2
        dna_line = lines[0][pad:]
        aa_line  = lines[1][pad:]
        assert dna_line == cds
        # AA positioned at codon middle: columns 1, 4, 7
        assert len(aa_line) == 9
        assert aa_line[1] == "M"
        assert aa_line[4] == "A"
        assert aa_line[7] == "S"
        # All other columns are whitespace
        for i in (0, 2, 3, 5, 6, 8):
            assert aa_line[i] == " "

    def test_aa_alignment_across_line_wrap(self):
        """A 12-bp CDS with line_width=6 (2 codons per line) must still
        place AAs centered under their codons on each wrapped line."""
        cds = "ATGGCCAGCAAA"  # M A S K
        t = sc._mut_build_preview_text(cds, line_width=6)
        lines = self._plain_lines(t)
        num_w = len(str(len(cds)))   # "12" → 2
        pad = num_w + 2
        # Line 0 = DNA row 0; line 1 = AA row 0; line 2 = DNA row 1; line 3 = AA row 1
        assert lines[0][pad:] == "ATGGCC"
        assert lines[1][pad:] == " M  A "
        assert lines[2][pad:] == "AGCAAA"
        assert lines[3][pad:] == " S  K "

    def test_stop_codon_shown_as_asterisk(self):
        """A trailing TAA stop should render as '*' in the AA row (not
        truncated like _mut_translate does)."""
        cds = "ATGTAA"  # M *
        t = sc._mut_build_preview_text(cds, line_width=6)
        lines = self._plain_lines(t)
        num_w = 1
        pad = num_w + 2
        assert lines[0][pad:] == "ATGTAA"
        assert lines[1][pad:] == " M  * "

    def test_mutation_substitutes_mutant_codon(self):
        """When mutation dict is passed, the DNA row must show the mutant
        codon (not the WT) and the AA row must show the mutant AA."""
        cds = "ATGTGGGCCTAA"   # M W A *
        mutation = {
            "wt_codon":    "TGG",
            "mut_codon":   "TTT",
            "nt_position": 4,   # 1-based: codon 2 occupies nt 4..6
        }
        t = sc._mut_build_preview_text(cds, mutation=mutation, line_width=12)
        lines = self._plain_lines(t)
        num_w = len(str(len(cds)))   # 2
        pad = num_w + 2
        assert lines[0][pad:] == "ATGTTTGCCTAA"   # TGG → TTT
        # AA row: positions 1=M, 4=F (mutated), 7=A, 10=*
        assert lines[1][pad:][1]  == "M"
        assert lines[1][pad:][4]  == "F"
        assert lines[1][pad:][7]  == "A"
        assert lines[1][pad:][10] == "*"

    def test_line_width_rounded_to_multiple_of_three(self):
        """Line width must be a multiple of 3 so codons don't straddle
        wrap boundaries — passing 10 should behave like 9."""
        cds = "ATGGCCAGCAAATTT"   # M A S K F — 15 bp
        t9  = sc._mut_build_preview_text(cds, line_width=9)
        t10 = sc._mut_build_preview_text(cds, line_width=10)
        # Both should render identically (10 rounds down to 9)
        assert t9.plain == t10.plain

    def test_aa_colored_purple_by_default(self):
        """AAs in the preview should carry the purple style so they visually
        separate from the green DNA backbone."""
        cds = "ATGGCCAGC"   # M A S
        t = sc._mut_build_preview_text(cds, line_width=9)
        # Walk the Text's style spans and find the AA chars
        rendered = t.render(console=None) if False else None   # unused
        purple = sc._MUT_PREVIEW_AA_COLOR
        # Each AA char should have a span whose style contains the purple color
        aa_chars_found: set[str] = set()
        for span in t.spans:
            if span.style and purple in str(span.style):
                # Extract the character(s) under this span
                aa_chars_found.update(t.plain[span.start:span.end])
        assert {"M", "A", "S"}.issubset(aa_chars_found)

    def test_mutation_highlight_uses_orange(self):
        cds = "ATGTGGGCCTAA"
        mutation = {"wt_codon": "TGG", "mut_codon": "TTT", "nt_position": 4}
        t = sc._mut_build_preview_text(cds, mutation=mutation, line_width=12)
        orange = sc._MUT_PREVIEW_MUT_COLOR
        orange_chars: list[str] = []
        for span in t.spans:
            if span.style and orange in str(span.style):
                orange_chars.extend(t.plain[span.start:span.end])
        # The three mutated DNA bases (TTT) and the mutant AA (F) should
        # all appear under an orange-colored span at least once.
        joined = "".join(orange_chars)
        assert joined.count("T") >= 3
        assert "F" in joined

    def test_cursor_adds_reverse_style(self):
        cds = "ATGGCCAGCAAA"   # M A S K
        t = sc._mut_build_preview_text(cds, cursor_aa=2, line_width=12)
        # AA index 2 is 'S' (codon AGC). The AA char should be styled with
        # reverse-video (on top of the purple color).
        reverse_aa_chars: list[str] = []
        for span in t.spans:
            if span.style and "reverse" in str(span.style):
                reverse_aa_chars.extend(t.plain[span.start:span.end])
        # The 3 DNA bases of the cursor codon AND the AA letter should
        # all be highlighted with reverse-video.
        joined = "".join(reverse_aa_chars)
        assert "S" in joined
        # Codon AGC occupies positions 6-8 in the DNA → A, G, C all reversed
        assert "A" in joined
        assert "G" in joined
        assert "C" in joined

    def test_cursor_minus_one_renders_no_reverse(self):
        """Default cursor_aa=-1 means no cursor — no reverse styling
        anywhere (except if there's a mutation, which uses its own
        style)."""
        cds = "ATGGCCAGCAAA"
        t = sc._mut_build_preview_text(cds, line_width=12)
        for span in t.spans:
            if span.style:
                assert "reverse" not in str(span.style)


# ── Click-to-AA index math  (_mut_click_to_aa_index) ──────────────────────────

class TestClickToAA:
    """Pure arithmetic that backs `_MutPreview.on_click`. No Textual
    required — we poke the helper directly with (vp_x, content_row)."""

    # A 12-bp CDS with 4 codons: ATG GCC AGC AAA → M A S K
    # line_width = 12 → everything fits on one logical line
    # pad = num_w + 2 = len("12") + 2 = 4

    def test_click_dna_row_hits_correct_codon(self):
        # content_row 0 == DNA row; columns 4..6 = codon 0 (ATG → M)
        for x in (4, 5, 6):
            assert sc._mut_click_to_aa_index(
                True, 12, 4, 12, 4, x, 0,
            ) == 0
        # columns 7..9 = codon 1 (GCC → A)
        for x in (7, 8, 9):
            assert sc._mut_click_to_aa_index(
                True, 12, 4, 12, 4, x, 0,
            ) == 1

    def test_click_aa_row_hits_same_codon_as_dna_row(self):
        # content_row 1 == AA row — same column math
        assert sc._mut_click_to_aa_index(True, 12, 4, 12, 4, 4, 1) == 0
        assert sc._mut_click_to_aa_index(True, 12, 4, 12, 4, 7, 1) == 1

    def test_click_on_prefix_returns_minus_one(self):
        # Columns 0..3 are the line-number prefix
        for x in range(0, 4):
            assert sc._mut_click_to_aa_index(True, 12, 4, 12, 4, x, 0) == -1

    def test_click_past_end_returns_minus_one(self):
        # Column 16 is past the DNA (12 cols after prefix ends at 16)
        assert sc._mut_click_to_aa_index(True, 12, 4, 12, 4, 16, 0) == -1

    def test_click_on_second_logical_line(self):
        # 24-bp CDS, line_width 12 → codons 0-3 on logical line 0,
        # codons 4-7 on logical line 1 (content_row 2 = line 1 DNA row)
        # pad = len("24") + 2 = 4
        # content_row 2 → logical_line 1 → bp_start = 12 → aa offset = 4
        assert sc._mut_click_to_aa_index(True, 24, 8, 12, 4, 4, 2) == 4
        assert sc._mut_click_to_aa_index(True, 24, 8, 12, 4, 7, 3) == 5  # AA row

    def test_aa_only_mode(self):
        # 20 aa protein wrapped at 10 per line → row 0 = aa 0..9, row 1 = aa 10..19
        assert sc._mut_click_to_aa_index(False, 0, 20, 10, 0, 0, 0)  == 0
        assert sc._mut_click_to_aa_index(False, 0, 20, 10, 0, 9, 0)  == 9
        assert sc._mut_click_to_aa_index(False, 0, 20, 10, 0, 3, 1)  == 13
        # Out of range
        assert sc._mut_click_to_aa_index(False, 0, 20, 10, 0, 0, 2)  == -1
        assert sc._mut_click_to_aa_index(False, 0, 20, 10, 0, 10, 1) == -1


# ── Cursor keyboard navigation  (_mut_next_cursor) ────────────────────────────

class TestCursorNav:
    def test_first_keypress_snaps_to_zero(self):
        # cursor=-1 (no cursor yet) → any direction places it at 0
        for d in ("left", "right", "up", "down"):
            assert sc._mut_next_cursor(-1, 50, 30, True,  d) == 0
            assert sc._mut_next_cursor(-1, 50, 10, False, d) == 0

    def test_left_right_by_one(self):
        assert sc._mut_next_cursor(5, 50, 30, True, "left")  == 4
        assert sc._mut_next_cursor(5, 50, 30, True, "right") == 6

    def test_left_right_clamp(self):
        assert sc._mut_next_cursor(0,  50, 30, True, "left")  == 0
        assert sc._mut_next_cursor(49, 50, 30, True, "right") == 49

    def test_up_down_step_dna_mode(self):
        # line_width=30 bp → 10 AAs per row
        assert sc._mut_next_cursor(15, 50, 30, True, "up")   == 5
        assert sc._mut_next_cursor(5,  50, 30, True, "down") == 15

    def test_up_down_clamp(self):
        # Can't go past row 0 or past protein length
        assert sc._mut_next_cursor(3,  50, 30, True, "up")   == 0
        assert sc._mut_next_cursor(45, 50, 30, True, "down") == 49

    def test_up_down_step_aa_only_mode(self):
        # line_width=10 AAs per row
        assert sc._mut_next_cursor(15, 50, 10, False, "up")   == 5
        assert sc._mut_next_cursor(5,  50, 10, False, "down") == 15

    def test_empty_protein_returns_minus_one(self):
        assert sc._mut_next_cursor(0, 0, 30, True, "right") == -1


# ── AA picker sub-modal ──────────────────────────────────────────────────────

class TestAAPicker:
    def test_catalog_contains_20_plus_stop(self):
        # 20 proteinogenic amino acids + stop = 21 entries total
        assert len(sc.AminoAcidPickerModal._AA_CATALOG) == 21
        codes = {a for (a, _, _) in sc.AminoAcidPickerModal._AA_CATALOG}
        assert codes == set("ACDEFGHIKLMNPQRSTVWY") | {"*"}

    def test_wt_aa_excluded_from_choices(self):
        modal = sc.AminoAcidPickerModal(position=140, wt_aa="W")
        assert "W" not in modal._choices
        # All other AAs should still be pickable
        assert set(modal._choices) == (
            set("ACDEFGHIKLMNPQRSTVWY") | {"*"}
        ) - {"W"}