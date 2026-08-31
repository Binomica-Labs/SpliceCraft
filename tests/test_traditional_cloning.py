"""
test_traditional_cloning — Restriction-digest + ligation engine.

Covers the module-level helpers that power ConstructorModal's "Traditional"
tab: `_enzyme_cuts`, `_digest_with_enzymes`, `_make_synthetic_fragment`,
`_ligate_fragments`, `_close_circular`, `_rc_fragment`,
`_simulate_traditional_cloning`. The UI is exercised separately in
test_smoke.py.

These are pure unit tests: no Textual app, no async setup. Each function
is a deterministic transform on a small synthetic SeqRecord-equivalent
seq string + feature dict list.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


# A synthetic 28 bp circular plasmid carrying one EcoRI site (GAATTC at 4..10)
# and one BamHI site (GGATCC at 18..24). Cuts:
#   EcoRI:  top=5,  bot=9   → 5' overhang AATT
#   BamHI:  top=19, bot=23  → 5' overhang GATC
TINY_PLASMID = "AAAAGAATTCAAAAAAAAGGATCCAAAA"


# ──────────────────────────────────────────────────────────────────────────────
# _enzyme_cuts
# ──────────────────────────────────────────────────────────────────────────────


class TestEnzymeCuts:
    def test_finds_palindromic_cut(self):
        cuts = sc._enzyme_cuts("AAAGAATTCAAA", ["EcoRI"], circular=False)
        assert len(cuts) == 1
        c = cuts[0]
        assert c["enzyme"] == "EcoRI"
        assert c["top"] == 4   # G^AATTC, recognition starts at 3
        assert c["bot"] == 8
        assert c["kind"] == "5'"
        assert c["overhang_seq"] == "AATT"

    def test_finds_three_prime_overhang(self):
        # KpnI: GGTAC^C (fwd=5, rev=1) → 3' overhang GTAC
        cuts = sc._enzyme_cuts("AAAGGTACCAAA", ["KpnI"], circular=False)
        assert len(cuts) == 1
        assert cuts[0]["kind"] == "3'"
        assert cuts[0]["overhang_seq"] == "GTAC"

    def test_finds_blunt_cut(self):
        # EcoRV: GAT^ATC (fwd=3, rev=3) → blunt
        cuts = sc._enzyme_cuts("AAAGATATCAAA", ["EcoRV"], circular=False)
        assert len(cuts) == 1
        assert cuts[0]["kind"] == "blunt"
        assert cuts[0]["overhang_seq"] == ""

    def test_returns_empty_for_unknown_enzyme(self):
        assert sc._enzyme_cuts("AAAA", ["NotARealEnzyme"]) == []

    def test_returns_empty_for_empty_seq(self):
        assert sc._enzyme_cuts("", ["EcoRI"]) == []

    def test_multiple_enzymes_sorted_by_position(self):
        cuts = sc._enzyme_cuts(TINY_PLASMID, ["BamHI", "EcoRI"], circular=True)
        assert [c["enzyme"] for c in cuts] == ["EcoRI", "BamHI"]
        assert cuts[0]["top"] < cuts[1]["top"]

    def test_circular_finds_origin_spanning_site(self):
        # Place GAATTC straddling the origin: TTCxxxxxxGAA. The 6 bp
        # recognition would be seq[-3:] + seq[:3] = "GAA" + "TTC" = "GAATTC".
        # _enzyme_cuts must find it via the augmented scan.
        seq = "TTCAAAAAAAAGAA"   # 14 bp; "GAATTC" at indices 11,12,13,0,1,2
        cuts = sc._enzyme_cuts(seq, ["EcoRI"], circular=True)
        assert len(cuts) == 1
        # Top cut at 11 + 1 = 12 (mod 14)
        assert cuts[0]["top"] == 12

    def test_linear_skips_origin_spanning_site(self):
        seq = "TTCAAAAAAAAGAA"
        cuts = sc._enzyme_cuts(seq, ["EcoRI"], circular=False)
        assert cuts == []

    def test_origin_wrapping_cut_overhang_is_correct(self):
        """Regression guard for 2026-05-05 fix: a cut whose top + bot
        positions straddle the origin (post-modulo, top > bot even
        though the cut is a 5' overhang) used to compute
        `seq[bot:top]` — a slice that can be the entire plasmid minus
        4 bp instead of the 4-bp overhang. Now the function uses
        pre-modulo values to detect wrap and emits `seq[lo:] + seq[:hi]`
        for the canonical 4-bp overhang."""
        # Place GAATTC straddling the origin: site spans [11..14, 0..1].
        # Top cut at (11+1)%14 = 12, bot cut at (11+5)%14 = 2 (mod 14).
        # Naive impl would have lo=2, hi=12 → "AAAAAAAAA" (10 bp).
        # Correct impl: overhang_len = 4, anchor at top=12, end=2 (wrap)
        # → seq[12:] + seq[:2] = "AA" + "TT" = "AATT".
        seq = "TTCAAAAAAAAAGAA"   # 15 bp; hmm, let's recompute exactly.
        # Use 14 bp: "TTCAAAAAAAAGAA" — site GAATTC at pos 11, 12, 13, 0, 1, 2.
        seq = "TTCAAAAAAAAGAA"   # 14 bp
        cuts = sc._enzyme_cuts(seq, ["EcoRI"], circular=True)
        assert len(cuts) == 1
        assert cuts[0]["overhang_seq"] == "AATT", (
            f"expected AATT for wrap cut, got {cuts[0]['overhang_seq']!r}"
        )
        # Kind detection: pre-modulo top < pre-modulo bot for fwd<rev
        # (5' overhang). The post-modulo values can flip apparent
        # ordering but kind must still report "5'".
        assert cuts[0]["kind"] == "5'"

    def test_dedupes_palindromic_double_emission(self):
        """Palindromic enzymes match on both strands at the same site;
        the function should report each physical cut exactly once."""
        cuts = sc._enzyme_cuts("AAAGAATTCAAA", ["EcoRI"], circular=False)
        # Even though EcoRI matches on top AND bottom strands at this
        # site, the result is one cut.
        assert len(cuts) == 1

    def test_non_palindromic_finds_reverse_strand_match(self):
        # BsaI (GGTCTC, type IIS, fwd=7, rev=11) is non-palindromic; its RC
        # is GAGACC, so GAGACC on the top strand = BsaI bound on the reverse
        # strand. With the site placed far enough from the 5' end that the
        # downstream-of-recognition cut still lands INSIDE the molecule, the
        # scanner must report that real cut at its forward-strand coords.
        #   top cut = p + site_len - rev_cut = 8 + 6 - 11 = 3
        #   bot cut = p + site_len - fwd_cut = 8 + 6 -  7 = 7
        seq = "A" * 8 + "GAGACC" + "A" * 8   # n=22, GAGACC at p=8
        cuts = sc._enzyme_cuts(seq, ["BsaI"], circular=False)
        assert len(cuts) == 1
        c = cuts[0]
        assert c["enzyme"] == "BsaI"
        assert c["top"] == 3
        assert c["bot"] == 7
        assert c["kind"] == "5'"
        assert c["overhang_seq"] == "AAAA"

    def test_reverse_strand_cut_off_5prime_end_dropped_on_linear(self):
        """A reverse-strand Type IIS enzyme whose cut would fall PAST the
        5' end produces no real cut on a LINEAR molecule — the enzyme binds
        but its scissile bond is off the end. Pre-fix the cut wrapped via
        `% n` into a phantom boundary near the 3' end (here top=-2 → 13)."""
        # GAGACC at p=3: top cut = 3 + 6 - 11 = -2 (off the 5' end).
        seq = "AAAGAGACCAAAAAA"   # n=15
        cuts = sc._enzyme_cuts(seq, ["BsaI"], circular=False)
        assert cuts == [], (
            f"expected no in-bounds cut on the linear molecule, got {cuts}"
        )

    def test_reverse_strand_cut_off_origin_wraps_on_circular(self):
        """The SAME off-the-end site IS a real cut on a CIRCULAR molecule:
        the scissile bond wraps around the origin. The guard that drops the
        linear phantom must NOT touch the circular wrap."""
        seq = "AAAGAGACCAAAAAA"   # n=15, GAGACC at p=3
        cuts = sc._enzyme_cuts(seq, ["BsaI"], circular=True)
        assert len(cuts) == 1
        c = cuts[0]
        assert c["enzyme"] == "BsaI"
        # top = (3 + 6 - 11) % 15 = 13 ; bot = (3 + 6 - 7) % 15 = 2
        assert c["top"] == 13
        assert c["bot"] == 2


# ──────────────────────────────────────────────────────────────────────────────
# _digest_with_enzymes — fragment slicing
# ──────────────────────────────────────────────────────────────────────────────


class TestDigestWithEnzymes:
    def test_circular_two_cuts_yields_two_fragments(self):
        frags = sc._digest_with_enzymes(TINY_PLASMID, ["EcoRI", "BamHI"],
                                          circular=True)
        assert len(frags) == 2
        # Total bp should equal the original.
        assert sum(len(f["top_seq"]) for f in frags) == len(TINY_PLASMID)
        # Both fragments have 5' overhangs at both ends.
        for f in frags:
            assert f["left"]["kind"] == "5'"
            assert f["right"]["kind"] == "5'"

    def test_circular_zero_cuts_yields_one_uncut_fragment(self):
        # No EcoRI in the input
        seq = "AAAACCCCGGGG"
        frags = sc._digest_with_enzymes(seq, ["EcoRI"], circular=True)
        assert len(frags) == 1
        assert frags[0]["top_seq"] == seq
        assert frags[0]["left"]["kind"] == "linear"

    def test_linear_two_cuts_yields_three_fragments(self):
        seq = "AAAAGAATTCAAAAAAAAGGATCCAAAA"
        frags = sc._digest_with_enzymes(seq, ["EcoRI", "BamHI"],
                                          circular=False)
        assert len(frags) == 3
        # First fragment: linear left edge, EcoRI right edge.
        assert frags[0]["left"]["kind"]  == "linear"
        assert frags[0]["right"]["enzyme"] == "EcoRI"
        # Last fragment: BamHI left, linear right.
        assert frags[-1]["left"]["enzyme"] == "BamHI"
        assert frags[-1]["right"]["kind"]   == "linear"

    def test_round_trip_through_ligate_and_close(self):
        """Cut the plasmid → re-ligate two fragments → close circle.
        The closed top_seq must contain the original (rotated) as a
        substring of itself doubled."""
        frags = sc._digest_with_enzymes(TINY_PLASMID, ["EcoRI", "BamHI"],
                                          circular=True)
        joined  = sc._ligate_fragments(frags[0], frags[1])
        closed  = sc._close_circular(joined)
        assert closed is not None
        assert closed["circular"] is True
        # Doubling the original lets us match any rotation of the closed seq.
        doubled = TINY_PLASMID + TINY_PLASMID
        assert closed["top_seq"] in doubled

    def test_wrap_feature_splits_into_head_and_tail(self):
        """Regression guard for 2026-05-05 fix: a feature with end < start
        (origin-spanning, e.g., oriV) used to route through the slotting
        algorithm assuming start ≤ end and end up in the wrong fragment.
        Now wrap features pre-split into a tail [start, n) + head [0, end)
        so each half lands in the correct fragment."""
        # 28 bp circular plasmid. Cuts at top=5 (EcoRI) and top=19 (BamHI)
        # produce two fragments: [5..19) (14 bp) and [19..5) (14 bp wrap).
        # Wrap feature [22..3) covers positions 22..27 + 0..2.
        feats = [{"start": 22, "end": 3, "label": "wrap-feat", "strand": 1}]
        frags = sc._digest_with_enzymes(TINY_PLASMID, ["EcoRI", "BamHI"],
                                          circular=True, features=feats)
        assert len(frags) == 2
        # The wrap feature should appear in the wrap fragment (the one
        # whose left edge is BamHI). Its halves: tail=[22, 28) → local
        # [3, 9); head=[0, 3) → local [9, 12).
        wrap_frag = next(f for f in frags
                          if f["left"]["enzyme"] == "BamHI")
        labeled = [f for f in wrap_frag["features"]
                    if f.get("label") == "wrap-feat"]
        # Two halves, both tagged with `_wrap_origin_split`.
        assert len(labeled) == 2
        splits = {f["_wrap_origin_split"] for f in labeled}
        assert splits == {"tail", "head"}

    def test_features_get_slotted_into_correct_fragment(self):
        # Plasmid with one feature in the EcoRI..BamHI window (insert region)
        feats = [{"start": 8, "end": 15, "label": "ins-feat", "strand": 1}]
        frags = sc._digest_with_enzymes(TINY_PLASMID, ["EcoRI", "BamHI"],
                                          circular=True, features=feats)
        # The feature at 8..15 falls into the fragment whose top range
        # starts at the EcoRI cut (top=5) and runs to the BamHI cut (top=19).
        labeled = [f for f in frags if any(
            ft.get("label") == "ins-feat" for ft in f["features"])]
        assert len(labeled) == 1
        # And the local coords are shifted: 8-5=3 ... 15-5=10.
        ft = next(ft for ft in labeled[0]["features"]
                    if ft.get("label") == "ins-feat")
        assert ft["start"] == 3
        assert ft["end"]   == 10

    def test_feature_crossing_the_wrap_fragment_start_is_not_inverted(self):
        """Regression: a feature spanning the cut at the ORIGIN-CROSSING
        fragment's 5' end used to come back INVERTED (`start > end`).

        The wrap fragment maps absolute→local modularly, so a piece that
        begins BEFORE the fragment lands ABOVE `frag_len`, and the old
        `min(new_s, frag_len)` clamp pinned that START to the fragment's
        END instead of to 0. The inverted feature rode straight through
        `_simulate_traditional_cloning` (which shifts start/end by the same
        offset) into the shipped product, where `_feat_len` reads it as a
        wrap — a 49 bp AmpR remnant rendered as a 489 bp "AmpR"."""
        # 28 bp circular. Cuts: EcoRI top=5, BamHI top=19.
        #   frag A = [5, 19)  len 14
        #   frag B = [19, 5)  len 14  (wraps the origin)
        # A feature [15, 24) crosses the BamHI cut at 19, so its TAIL half
        # belongs to the wrap fragment at local [0, 5)  (abs 19..24).
        feats = [{"start": 15, "end": 24, "label": "span", "strand": 1}]
        frags = sc._digest_with_enzymes(TINY_PLASMID, ["EcoRI", "BamHI"],
                                        circular=True, features=feats)
        wrap_frag = next(f for f in frags if f["left"]["enzyme"] == "BamHI")
        tail = next(ft for ft in wrap_frag["features"]
                    if ft.get("label") == "span")
        assert tail["start"] == 0
        assert tail["end"] == 5
        assert tail["start"] <= tail["end"]

    def test_single_cut_linearisation_splits_a_spanning_feature(self):
        """One cut on a circular plasmid yields ONE linear fragment. A
        feature spanning that cut is genuinely split across the new 5' and
        3' ends — it must come back as two in-bounds pieces, never as a
        single `start > end` "wrap" (meaningless on a linear fragment)."""
        # TINY_PLASMID is 28 bp with a single EcoRI cut at top=5.
        feats = [{"start": 2, "end": 9, "label": "span", "strand": 1}]
        frags = sc._digest_with_enzymes(TINY_PLASMID, ["EcoRI"],
                                        circular=True, features=feats)
        assert len(frags) == 1
        n = len(frags[0]["top_seq"])
        pieces = [ft for ft in frags[0]["features"]
                  if ft.get("label") == "span"]
        assert len(pieces) == 2
        for p in pieces:
            assert 0 <= p["start"] <= p["end"] <= n
        # abs [2, 5) → local [25, 28); abs [5, 9) → local [0, 4).
        assert sorted((p["start"], p["end"]) for p in pieces) == [(0, 4), (25, 28)]
        # No annotated base is lost by the split.
        assert sum(p["end"] - p["start"] for p in pieces) == 7

    def test_fragment_wholly_inside_a_feature_still_gets_it(self):
        """Regression: features were slotted by ENDPOINT — a piece went to the
        fragment holding `start` and the one holding `end-1`, and nowhere else.
        A feature longer than one fragment therefore lost every fragment
        BETWEEN its ends, even though those are 100% inside the gene."""
        seq = "".join("ACGT"[i % 4] for i in range(300))
        seq = seq[:40] + "GAATTC" + seq[46:]      # EcoRI   → cut 41
        seq = seq[:140] + "GGATCC" + seq[146:]    # BamHI   → cut 141
        seq = seq[:240] + "AAGCTT" + seq[246:]    # HindIII → cut 241
        enz = ["EcoRI", "BamHI", "HindIII"]
        feats = [{"start": 30, "end": 260, "label": "longGene", "strand": 1}]
        frags = sc._digest_with_enzymes(seq, enz, circular=True, features=feats)
        assert len(frags) == 3
        for i, f in enumerate(frags):
            got = [ft for ft in f["features"] if ft.get("label") == "longGene"]
            assert got, f"fragment {i} is inside the gene but lost the annotation"
        # The two fully-enclosed fragments are covered end to end.
        mids = [f for f in frags
                if any(ft.get("_split") == "mid" for ft in f["features"])]
        assert len(mids) == 2
        for f in mids:
            ft = next(x for x in f["features"] if x.get("label") == "longGene")
            assert (ft["start"], ft["end"]) == (0, len(f["top_seq"]))

    def test_feature_wrapping_past_every_cut_is_not_read_as_contained(self):
        """A circular feature long enough to run past ALL the cuts lands both
        endpoints back in the SAME fragment. Endpoint-slotting read that as
        'contained in one fragment' and dropped it from every other fragment
        it actually covers."""
        # 28 bp, cuts at 5 (EcoRI) and 19 (BamHI). A feature [22, 18) wraps
        # the origin, so it is pre-split into [22, 28) + [0, 18) — the latter
        # runs from before cut 5, past it, to inside fragment [5, 19).
        feats = [{"start": 22, "end": 18, "label": "big", "strand": 1}]
        frags = sc._digest_with_enzymes(TINY_PLASMID, ["EcoRI", "BamHI"],
                                        circular=True, features=feats)
        for i, f in enumerate(frags):
            assert any(ft.get("label") == "big" for ft in f["features"]), (
                f"fragment {i} lost the wrap feature")

    def test_no_fragment_feature_is_ever_inverted(self):
        """Property guard across both topologies and cut counts: a
        fragment is a LINEAR molecule, so every feature slotted onto it
        must satisfy `0 <= start <= end <= len(top_seq)`."""
        import random
        rng = random.Random(20260724)
        enz = ["EcoRI", "BamHI", "HindIII", "PstI", "XhoI"]
        for _ in range(300):
            seq = "".join(rng.choice("ACGT") for _ in range(120))
            for name in rng.sample(enz, rng.randint(1, 3)):
                site = sc._NEB_ENZYMES[name][0]
                p = rng.randrange(len(seq) - len(site))
                seq = seq[:p] + site + seq[p + len(site):]
            fs, fe = rng.randrange(len(seq)), rng.randrange(len(seq))
            for circ in (True, False):
                # `end < start` means "wraps the origin" — well-formed only on
                # a circular molecule; on a linear one it is malformed input.
                feats = [{"start": fs, "end": fe if circ else max(fs, fe),
                          "label": "f", "strand": 1}]
                for frag in sc._digest_with_enzymes(seq, enz, circular=circ,
                                                    features=feats):
                    L = len(frag["top_seq"])
                    for ft in frag["features"]:
                        assert 0 <= ft["start"] <= ft["end"] <= L, (
                            f"circular={circ} feat={feats} -> {ft} on a "
                            f"{L} bp fragment")


# ──────────────────────────────────────────────────────────────────────────────
# _make_synthetic_fragment (modes b + c — PCR product / feature with tails)
# ──────────────────────────────────────────────────────────────────────────────


class TestMakeSyntheticFragment:
    def test_stamps_canonical_overhangs(self):
        frag = sc._make_synthetic_fragment(
            "GAGCATGAAACGGCCAAGTAA",
            enz_left="EcoRI", enz_right="BamHI",
            source_label="myPCR",
        )
        assert frag["left"]["overhang_seq"]  == "AATT"
        assert frag["left"]["enzyme"]        == "EcoRI"
        assert frag["right"]["overhang_seq"] == "GATC"
        assert frag["right"]["enzyme"]       == "BamHI"

    def test_three_prime_overhang_kind(self):
        frag = sc._make_synthetic_fragment(
            "AAA", enz_left="KpnI", enz_right="EcoRI",
        )
        assert frag["left"]["kind"]  == "3'"
        assert frag["right"]["kind"] == "5'"

    def test_blunt_kind(self):
        frag = sc._make_synthetic_fragment(
            "AAA", enz_left="EcoRV", enz_right="EcoRV",
        )
        assert frag["left"]["kind"]  == "blunt"
        assert frag["right"]["kind"] == "blunt"
        assert frag["left"]["overhang_seq"]  == ""
        assert frag["right"]["overhang_seq"] == ""

    def test_unknown_enzyme_raises(self):
        with pytest.raises(ValueError):
            sc._make_synthetic_fragment("AAA", enz_left="NotReal",
                                          enz_right="EcoRI")

    def test_type_iis_rejected_with_clear_message(self):
        """Regression guard for 2026-05-05 fix: BsaI / Esp3I / BsmBI
        cut OUTSIDE their recognition (`fwd_cut > site_len`), so the
        synthetic-fragment model can't know what the overhang bases
        will be — those depend on the surrounding context the user
        supplies. Used to silently produce an empty-overhang fragment
        that wouldn't ligate; now raises ValueError pointing the
        user at "From plasmid" mode."""
        with pytest.raises(ValueError) as exc:
            sc._make_synthetic_fragment("AAACCCGGG",
                                          enz_left="BsaI",
                                          enz_right="EcoRI")
        assert "Type IIS" in str(exc.value)
        # Both ends checked — also rejects when the right enzyme is IIS.
        with pytest.raises(ValueError) as exc:
            sc._make_synthetic_fragment("AAACCCGGG",
                                          enz_left="EcoRI",
                                          enz_right="Esp3I")
        assert "Type IIS" in str(exc.value)


# ──────────────────────────────────────────────────────────────────────────────
# _excise_pcr_insert — digest + gel-purify the insert (NO primer-pad carryover)
# ──────────────────────────────────────────────────────────────────────────────


class TestExcisePcrInsert:
    """The bench truth: a tailed PCR product cut with the cloning enzymes is
    purified down to its MIDDLE fragment — the GCGC primer pad and the
    recognition-site bases OUTSIDE the cut are off-cuts you throw away, so they
    NEVER reach the clone. Regression for the user-reported "extra bases from
    the GCGC SalI padding got integrated into the payload"."""

    PAYLOAD = "ATGAAACGTACGGCCTTAACGTACGT"        # 26 bp, no internal sites

    def test_tailed_amplicon_drops_pad_and_keeps_overhangs(self):
        """SalI both ends (the user's enzyme — palindromic 5' overhang)."""
        amp = "GCGC" + "GTCGAC" + self.PAYLOAD + "GTCGAC" + sc._rc("GCGC")
        frag, err = sc._excise_pcr_insert(amp, "SalI", "SalI")
        assert err is None, err
        assert "GCGC" not in frag["top_seq"], "GCGC pad leaked into the insert"
        assert self.PAYLOAD in frag["top_seq"], "payload mangled"
        # G^TCGAC: insert retains TCGAC at the 5' end, is recessed to G at 3'.
        assert frag["top_seq"].startswith("TCGAC")
        assert frag["top_seq"].endswith("G")
        assert frag["left"]["overhang_seq"] == "TCGA"
        assert frag["right"]["overhang_seq"] == "TCGA"
        assert frag["left"]["kind"] == "5'" and frag["right"]["kind"] == "5'"

    def test_bare_payload_is_wrapped_then_purified(self):
        """The manual 'Paste DNA' donor is the BARE region — wrap it in the
        tails the primers would add, then digest to the same purified insert.
        No pad, real EcoRI/BamHI overhangs."""
        frag, err = sc._excise_pcr_insert(self.PAYLOAD, "EcoRI", "BamHI")
        assert err is None, err
        assert "GCGC" not in frag["top_seq"]
        assert self.PAYLOAD in frag["top_seq"]
        assert frag["top_seq"].startswith("AATTC")     # EcoRI G^AATTC
        assert frag["top_seq"].endswith("G")           # BamHI G^GATCC
        assert frag["left"]["overhang_seq"] == "AATT"
        assert frag["right"]["overhang_seq"] == "GATC"

    def test_three_prime_overhang_reforms_site(self):
        """KpnI (GGTAC^C, 3' overhang) — the retained site bases live on the
        OPPOSITE fragment, but the recognition site still reforms on ligation
        and no pad survives."""
        amp = "GCGC" + "GGTACC" + self.PAYLOAD + "GGTACC" + sc._rc("GCGC")
        frag, err = sc._excise_pcr_insert(amp, "KpnI", "KpnI")
        assert err is None, err
        assert "GCGC" not in frag["top_seq"]
        assert self.PAYLOAD in frag["top_seq"]
        assert frag["left"]["kind"] == "3'" and frag["left"]["overhang_seq"] == "GTAC"

    def test_features_land_on_payload_after_purify(self):
        """An amplicon-local feature rides onto the purified insert at the
        offset the retained 5' site bases imply (EcoRI keeps AATTC = 5 nt)."""
        lead = len("GCGC") + len("GAATTC")             # 10
        amp = "GCGC" + "GAATTC" + self.PAYLOAD + "GGATCC" + sc._rc("GCGC")
        feats = [{"start": lead + 12, "end": lead + 21, "type": "CDS",
                  "label": "MyGene", "strand": 1}]
        frag, err = sc._excise_pcr_insert(amp, "EcoRI", "BamHI", features=feats)
        assert err is None, err
        g = [f for f in frag["features"] if f.get("label") == "MyGene"]
        assert g, "feature dropped from purified insert"
        assert g[0]["end"] - g[0]["start"] == 9        # length preserved
        assert g[0]["start"] == 5 + 12                 # AATTC lead + payload off

    def test_internal_site_refused(self):
        """A copy of the enzyme site INSIDE the region would cut the amplicon
        internally — refuse rather than ship a wrong/short clone."""
        amp = ("GCGC" + "GAATTC" + "ATG" + "GAATTC" + "AAACGT"
               + "GGATCC" + sc._rc("GCGC"))
        frag, err = sc._excise_pcr_insert(amp, "EcoRI", "BamHI")
        assert frag is None and err and "inside" in err.lower()

    def test_type_iis_refused(self):
        frag, err = sc._excise_pcr_insert(self.PAYLOAD, "BsaI", "EcoRI")
        assert frag is None and err and "Type IIS" in err

    def test_no_pad_in_final_clone_end_to_end(self):
        """The whole point: simulate the ligation and prove the FINISHED clone
        carries no GCGC pad and reforms the SalI site at the junction."""
        amp = "GCGC" + "GTCGAC" + self.PAYLOAD + "GTCGAC" + sc._rc("GCGC")
        ins, err = sc._excise_pcr_insert(amp, "SalI", "SalI")
        assert err is None, err
        vec = ("TTTTTTTTGTCGACAAAACCCCGGGGTTTTGTCGAC"
               + "AAAAAAAAAAAA")                        # 2 SalI sites
        vfrags, verr = sc._excise_fragment_pair(vec, ["SalI"], circular=True)
        assert verr is None, verr
        backbone = sc._pick_backbone_fragment(vfrags) or vfrags[0]
        res = sc._simulate_traditional_cloning(ins, backbone)
        prod = res["forward"]["top_seq"]
        assert "GCGC" not in prod, "pad leaked into the final clone"
        assert self.PAYLOAD in prod
        # One SalI reforms at the internal junction; the other wraps the origin.
        assert (prod + prod[:5]).count("GTCGAC") >= 2


# ──────────────────────────────────────────────────────────────────────────────
# _ligate_fragments + _close_circular (compatibility predicate)
# ──────────────────────────────────────────────────────────────────────────────


class TestLigateAndClose:
    def test_compatible_palindromic_ligates(self):
        a = sc._make_synthetic_fragment("AAA", enz_left="EcoRV",
                                          enz_right="EcoRI")
        b = sc._make_synthetic_fragment("CCC", enz_left="EcoRI",
                                          enz_right="EcoRV")
        merged = sc._ligate_fragments(a, b)
        assert merged is not None
        assert merged["top_seq"] == "AAACCC"
        # Outer ends preserved.
        assert merged["left"]  == a["left"]
        assert merged["right"] == b["right"]

    def test_incompatible_kinds_dont_ligate(self):
        # 5' overhang can't ligate to 3' overhang even if the bases match.
        a = sc._make_synthetic_fragment("AAA", enz_left="EcoRV",
                                          enz_right="EcoRI")
        # Hand-craft a 3' end with the same sequence as EcoRI's 5' overhang.
        b = {
            "top_seq": "CCC",
            "left":  {"overhang_seq": "AATT", "kind": "3'", "enzyme": "fake"},
            "right": {"overhang_seq": "", "kind": "linear", "enzyme": ""},
            "features": [], "source_label": "",
        }
        assert sc._ligate_fragments(a, b) is None

    def test_linear_edge_never_ligates(self):
        a = sc._make_synthetic_fragment("AAA", enz_left="EcoRI",
                                          enz_right="EcoRI")
        b = {
            "top_seq": "CCC",
            "left":  {"overhang_seq": "", "kind": "linear", "enzyme": ""},
            "right": {"overhang_seq": "", "kind": "linear", "enzyme": ""},
            "features": [], "source_label": "",
        }
        assert sc._ligate_fragments(a, b) is None

    def test_close_circular_palindromic_succeeds(self):
        # Synthetic fragment with EcoRI on both ends → ligates to a circle.
        frag = sc._make_synthetic_fragment("AAACCC", enz_left="EcoRI",
                                              enz_right="EcoRI")
        closed = sc._close_circular(frag)
        assert closed is not None
        assert closed["circular"] is True

    def test_close_circular_mismatched_ends_fails(self):
        frag = sc._make_synthetic_fragment("AAA", enz_left="EcoRI",
                                              enz_right="BamHI")
        assert sc._close_circular(frag) is None


# ──────────────────────────────────────────────────────────────────────────────
# _rc_fragment
# ──────────────────────────────────────────────────────────────────────────────


class TestRcFragment:
    def test_top_seq_reverse_complemented(self):
        frag = sc._make_synthetic_fragment(
            "ATGAAACG", enz_left="EcoRI", enz_right="BamHI",
        )
        rc = sc._rc_fragment(frag)
        assert rc["top_seq"] == "CGTTTCAT"

    def test_left_right_swapped(self):
        frag = sc._make_synthetic_fragment(
            "AAA", enz_left="EcoRI", enz_right="BamHI",
        )
        rc = sc._rc_fragment(frag)
        # After RC: the EcoRI end (was left) is now on the right; BamHI is on left.
        assert rc["left"]["enzyme"]  == "BamHI"
        assert rc["right"]["enzyme"] == "EcoRI"
        # Overhangs are RC'd; both palindromic so values unchanged.
        assert rc["left"]["overhang_seq"]  == "GATC"
        assert rc["right"]["overhang_seq"] == "AATT"

    def test_features_flipped_correctly(self):
        frag = sc._make_synthetic_fragment(
            "AAACCCGGG",  # 9 bp
            enz_left="EcoRI", enz_right="BamHI",
            features=[{"start": 1, "end": 4, "label": "x", "strand": 1}],
        )
        rc = sc._rc_fragment(frag)
        # 9 - end (4) → 5; 9 - start (1) → 8
        ft = rc["features"][0]
        assert ft["start"]  == 5
        assert ft["end"]    == 8
        assert ft["strand"] == -1


# ──────────────────────────────────────────────────────────────────────────────
# _simulate_traditional_cloning (full vector + insert pipeline)
# ──────────────────────────────────────────────────────────────────────────────


class TestSimulateTraditionalCloning:
    def test_directional_cloning_only_one_orientation_works(self):
        # EcoRI + BamHI: different non-cross-compatible overhangs → directional.
        insert = sc._make_synthetic_fragment(
            "GAGCATGAAACGGCCAAGTAA",
            enz_left="EcoRI", enz_right="BamHI",
            source_label="insert",
        )
        # Vector with EcoRI on right, BamHI on left so the insert ligates
        # forward.
        vector = {
            "top_seq": "TGGCCCC" * 10,   # 70 bp dummy backbone
            "left":  {"overhang_seq": "GATC", "kind": "5'", "enzyme": "BamHI"},
            "right": {"overhang_seq": "AATT", "kind": "5'", "enzyme": "EcoRI"},
            "features":     [],
            "source_label": "vector",
        }
        result = sc._simulate_traditional_cloning(insert, vector)
        # Only forward should be compatible: vector.right(AATT) ↔
        # insert.left(AATT) ✓; vector.left(GATC) ↔ insert.right(GATC) ✓.
        assert result["forward"]["compatible"] is True
        assert result["reverse"]["compatible"] is False
        # Warnings include a directional-cloning message.
        assert any("Directional" in w for w in result["warnings"])
        assert result["errors"] == []

    def test_ambiguous_orientation_emits_warning(self):
        # Single-enzyme palindromic cut: insert can ligate either way.
        insert = sc._make_synthetic_fragment(
            "GAGCATGAAACGGCCAAGTAA",
            enz_left="EcoRI", enz_right="EcoRI",
            source_label="insert",
        )
        vector = {
            "top_seq": "TGGCCCC" * 10,
            "left":  {"overhang_seq": "AATT", "kind": "5'", "enzyme": "EcoRI"},
            "right": {"overhang_seq": "AATT", "kind": "5'", "enzyme": "EcoRI"},
            "features":     [],
            "source_label": "vector",
        }
        result = sc._simulate_traditional_cloning(insert, vector)
        assert result["forward"]["compatible"] is True
        assert result["reverse"]["compatible"] is True
        assert any("Ambiguous" in w for w in result["warnings"])

    def test_no_compatibility_raises_error(self):
        insert = sc._make_synthetic_fragment(
            "AAA", enz_left="EcoRI", enz_right="EcoRI",
        )
        vector = {
            "top_seq": "CCC",
            "left":  {"overhang_seq": "GATC", "kind": "5'", "enzyme": "BamHI"},
            "right": {"overhang_seq": "GATC", "kind": "5'", "enzyme": "BamHI"},
            "features":     [],
            "source_label": "vector",
        }
        result = sc._simulate_traditional_cloning(insert, vector)
        assert result["forward"]["compatible"] is False
        assert result["reverse"]["compatible"] is False
        assert any("Neither orientation" in e for e in result["errors"])

    def test_compat_compatible_overhang_enzymes(self):
        """BamHI (GATC) and BglII (GATC) produce the same canonical
        overhang — fragments cut with either should ligate."""
        insert = sc._make_synthetic_fragment(
            "AAA", enz_left="BamHI", enz_right="BamHI",
        )
        vector = {
            "top_seq": "CCC",
            "left":  {"overhang_seq": "GATC", "kind": "5'", "enzyme": "BglII"},
            "right": {"overhang_seq": "GATC", "kind": "5'", "enzyme": "BglII"},
            "features":     [],
            "source_label": "vector",
        }
        result = sc._simulate_traditional_cloning(insert, vector)
        assert result["forward"]["compatible"] is True

    def test_excise_fragment_pair_circular_two_cuts(self):
        frags, err = sc._excise_fragment_pair(
            TINY_PLASMID, ["EcoRI", "BamHI"], circular=True,
        )
        assert err is None
        assert len(frags) == 2

    def test_excise_fragment_pair_no_cuts_errors_with_message(self):
        frags, err = sc._excise_fragment_pair(
            "AAACCCGGG", ["EcoRI"], circular=True,
        )
        assert err is not None
        assert "no cut sites" in err["error"]

    def test_excise_fragment_pair_one_cut_circular_errors(self):
        # Plasmid with exactly one EcoRI: can't excise an insert.
        seq = "AAAGAATTCAAAACCCCGGGGTTTTAAAA"
        frags, err = sc._excise_fragment_pair(
            seq, ["EcoRI"], circular=True,
        )
        assert err is not None
        assert "≥2 cuts" in err["error"]

    def test_excise_fragment_pair_three_cuts_circular_errors(self):
        """Regression guard for 2026-05-06 fix: ≥3 cuts on a circular
        plasmid is ambiguous — the helper used to silently return the
        full fragment list with err=None, and a future caller blindly
        picking ``fragments[0:2]`` would ship a wrong product. Helper
        now surfaces a clear error so the user can pick a different
        enzyme pair."""
        # Three EcoRI sites on a circular plasmid.
        seq = ("GAATTC" + "A" * 10 +
               "GAATTC" + "A" * 10 +
               "GAATTC" + "A" * 10)
        frags, err = sc._excise_fragment_pair(
            seq, ["EcoRI"], circular=True,
        )
        assert err is not None
        assert "exactly 2" in err["error"]
        assert "3" in err["error"]   # tells the user how many were found
        # Caller must surface this; the fragment list is still returned
        # for diagnostic display (don't break the existing return shape).
        assert isinstance(frags, list)

    @pytest.mark.parametrize("n_sites", [3, 4, 5, 7])
    def test_excise_fragment_pair_many_cuts_circular_errors(self, n_sites):
        """Sacred invariant #25 generalised: ANY count >2 must error.
        Parameterised so a regression that special-cases 3 (e.g.
        `if n_cuts == 3:`) can't slip through against 4+. Each test
        builds a circular sequence with `n_sites` EcoRI sites + 10 bp
        spacers."""
        seq = "".join("GAATTC" + "A" * 10 for _ in range(n_sites))
        frags, err = sc._excise_fragment_pair(
            seq, ["EcoRI"], circular=True,
        )
        assert err is not None, (
            f"{n_sites} cuts on circular should error; got err=None"
        )
        assert "exactly 2" in err["error"]
        assert str(n_sites) in err["error"]

    def test_excise_fragment_pair_two_enzymes_three_total_cuts_errors(self):
        """Invariant #25: total cut count is what matters, not per-
        enzyme count. EcoRI×2 + BamHI×1 = 3 total → must error.
        Catches a hypothetical regression that only checks per-enzyme
        counts."""
        # Two EcoRI + one BamHI, all in a circular sequence.
        seq = (
            "GAATTC" + "A" * 10 +
            "GAATTC" + "A" * 10 +
            "GGATCC" + "A" * 10
        )
        frags, err = sc._excise_fragment_pair(
            seq, ["EcoRI", "BamHI"], circular=True,
        )
        assert err is not None
        assert "exactly 2" in err["error"]
        # Per-enzyme breakdown is in the message.
        assert "EcoRI" in err["error"]
        assert "BamHI" in err["error"]

    def test_excise_fragment_pair_linear_three_cuts_does_not_error(self):
        """Invariant #25's strict 2-cut requirement applies only to
        CIRCULAR plasmids. Linear sequences with N cuts produce N+1
        fragments naturally, no ambiguity in selecting "the insert"
        because the caller knows linear has ends. The check must NOT
        false-positive on the linear path."""
        seq = (
            "GAATTC" + "A" * 10 +
            "GAATTC" + "A" * 10 +
            "GAATTC" + "A" * 10
        )
        frags, err = sc._excise_fragment_pair(
            seq, ["EcoRI"], circular=False,
        )
        # Linear is allowed any cut count without the >2 hard-stop.
        # (`err` may still be set for OTHER reasons — e.g. zero cuts —
        # but should NOT be the "exactly 2" message.)
        if err is not None:
            assert "exactly 2" not in err["error"], (
                f"linear path should NOT trigger exactly-2 hard-stop; "
                f"got err={err['error']!r}"
            )

    # ──────────────────────────────────────────────────────────────────
    # End-to-end UI tests — open the modal, switch tabs, simulate.
    # ──────────────────────────────────────────────────────────────────

    @staticmethod
    async def _setup_lane(pilot, modal, app, *,
                          vec_name="pVec", donor_pcr_seq=None,
                          donor_pcr_name="myInsert",
                          donor_plasmid_name=None,
                          e1="EcoRI", e2="BamHI"):
        """Wire the Constructor's Traditional pane through the new
        lane-baked-vector flow (2026-05-23):
          1. Plasmid mode → cursor pVec → Add → becomes backbone.
          2. Set backbone enzymes via the edit panel.
          3. PCR / plasmid mode → add donor → set its enzymes.
        Returns the trad pane for follow-up assertions."""
        from textual.widgets import (TabbedContent, RadioButton,
                                      Input, TextArea, Select,
                                      DataTable, Button)
        modal.query_one("#ctor-tabs", TabbedContent).active = \
            "ctor-tab-traditional"
        await pilot.pause(); await pilot.pause(0.05)
        entries = sorted(
            (e for e in sc._load_library() if isinstance(e, dict)),
            key=lambda e: sc._natural_sort_key(
                e.get("name") or e.get("id") or ""
            ),
        )
        # ── Backbone (plasmid mode) ──
        modal.query_one(
            "#trad-mode-plasmid", RadioButton,
        ).value = True
        await pilot.pause()
        vec_idx = next(i for i, e in enumerate(entries)
                        if e.get("name") == vec_name)
        modal.query_one(
            "#trad-source-table", DataTable,
        ).move_cursor(row=vec_idx)
        await pilot.pause()
        modal.query_one(
            "#btn-trad-add-frag", Button,
        ).press()
        await pilot.pause()
        modal.query_one("#trad-edit-enz-1", Select).value = e1
        await pilot.pause()
        modal.query_one("#trad-edit-enz-2", Select).value = e2
        await pilot.pause()
        # ── Donor ──
        if donor_plasmid_name is not None:
            modal.query_one(
                "#trad-mode-plasmid", RadioButton,
            ).value = True
            await pilot.pause()
            donor_idx = next(i for i, e in enumerate(entries)
                              if e.get("name") == donor_plasmid_name)
            modal.query_one(
                "#trad-source-table", DataTable,
            ).move_cursor(row=donor_idx)
            await pilot.pause()
        else:
            modal.query_one(
                "#trad-mode-pcr", RadioButton,
            ).value = True
            await pilot.pause()
            modal.query_one(
                "#trad-pcr-name", Input,
            ).value = donor_pcr_name
            modal.query_one(
                "#trad-pcr-seq", TextArea,
            ).text = (donor_pcr_seq
                       or "GAGCATGAAACGGCCAAGTAA")
        modal.query_one(
            "#btn-trad-add-frag", Button,
        ).press()
        await pilot.pause()
        modal.query_one("#trad-edit-enz-1", Select).value = e1
        await pilot.pause()
        modal.query_one("#trad-edit-enz-2", Select).value = e2
        await pilot.pause()
        return modal.query_one("#ctor-trad-pane",
                                  sc.TraditionalCloningPane)

    async def test_constructor_opens_with_tabbed_content(
            self, tiny_record, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            from textual.widgets import TabbedContent
            tabs = modal.query_one("#ctor-tabs", TabbedContent)
            # Three tabs registered (2026-05-07): the single modular
            # tab was split into per-grammar tabs (Golden Braid,
            # MoClo Plant) so each gets its own parts palette filtered
            # by the active grammar.
            tab_ids = [p.id for p in modal.query("TabPane")]
            assert "ctor-tab-traditional" in tab_ids
            assert "ctor-tab-gb_l0"       in tab_ids
            assert "ctor-tab-moclo_plant" in tab_ids
            # Traditional opens by default.
            assert tabs.active == "ctor-tab-traditional"

    async def test_source_filter_maps_row_to_correct_plasmid(
            self, tiny_record, isolated_library):
        """[E2 safety] The source-filter box must keep the row→entry mapping
        exact: a filtered row resolves to the plasmid SHOWN on that row, never
        the entry at the same index in the UNfiltered list (which would
        silently clone the wrong fragment)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from tests.test_smoke import _build_app, TERMINAL_SIZE

        def _entry(name, ln):
            r = SeqRecord(Seq("ACGT" * (ln // 4)), id=name, name=name)
            r.annotations["molecule_type"] = "DNA"
            r.annotations["topology"] = "circular"
            return {"id": name, "name": name, "size": ln,
                    "gb_text": sc._record_to_gb_text(r)}
        # sorted: Alpha(100), Beta(200), Delta(300), Gamma(400)
        sc._save_library([_entry("Alpha", 100), _entry("Beta", 200),
                          _entry("Delta", 300), _entry("Gamma", 400)])
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane", sc.TraditionalCloningPane)
            mine = {"Alpha", "Beta", "Delta", "Gamma"}
            # The mapping INVARIANT: under any filter, the record resolved for
            # display row i is EXACTLY the entry shown on row i (its bp size is
            # unique per entry, so a mismatch means the resolver indexed the
            # unfiltered list — the wrong-plasmid bug). 'e' matches Beta+Delta
            # (a strict subset whose row indices differ from the full list).
            for flt in ("", "e", "a", "lph"):
                pane._trad_source_filter = flt
                pane._populate_library_tables()
                await pilot.pause()
                filtered = pane._filtered_source_entries()
                for i, ent in enumerate(filtered):
                    if ent["name"] not in mine:
                        continue
                    r = pane._record_for_table_row(i, "#trad-source-table")
                    assert r is not None and len(r.seq) == ent["size"], (
                        flt, i, ent["name"], len(r.seq) if r else None)

    async def test_traditional_pane_pcr_mode_simulate(
            self, tiny_record, isolated_library):
        """End-to-end: backbone-as-first-lane-row + PCR-mode donor →
        Simulate produces forward + reverse products."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        vec = SeqRecord(
            Seq("AAAAAAAAGAATTCSTUFFERSTUFFERGGATCCTTTTTTTT"
                .replace("STUFFER", "AAAAA")),
            id="pVec", name="pVec",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        buf = io.StringIO()
        SeqIO.write(vec, buf, "genbank")
        sc._save_library([{
            "id": "pVec", "name": "pVec",
            "gb_text": buf.getvalue(),
            "size": len(vec.seq),
        }])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from textual.widgets import Button
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            pane = await self._setup_lane(pilot, modal, app)
            # Sanity: lane has backbone + donor configured.
            assert len(pane._lane_inserts) == 2
            assert pane._lane_inserts[0]["role"] == "backbone"
            assert pane._lane_inserts[0]["enz_left"] == "EcoRI"
            assert pane._lane_inserts[1]["role"] == "donor"
            assert pane._lane_inserts[1]["enz_left"] == "EcoRI", (
                f"donor enzymes: {pane._lane_inserts[1]}"
            )
            modal.query_one("#btn-trad-simulate", Button).press()
            await pilot.pause(); await pilot.pause(0.2)
            from textual.widgets import Static
            results = str(
                modal.query_one("#trad-results-text", Static).content
            )
            assert isinstance(pane._fwd_product, dict), (
                f"results: {results!r}"
            )
            assert isinstance(pane._rev_product, dict)
            assert "top_seq" in pane._fwd_product
            assert pane._fwd_product["top_seq"]
            assert "compatible" in pane._fwd_product

    async def test_traditional_pane_save_forward_to_library(
            self, tiny_record, isolated_library):
        """After Simulate the enabled Save button opens
        ``NamePlasmidModal``; accepting the user-chosen name writes
        the entry through the same library-save chain. The persisted
        GenBank text round-trips and the entry carries
        ``source: 'traditional:fwd'`` so future cascade hooks can
        identify trad-origin entries."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        vec = SeqRecord(
            Seq("AAAAAAAAGAATTCAAAAAAAAAAGGATCCTTTTTTTT"),
            id="pVec", name="pVec",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        buf = io.StringIO()
        SeqIO.write(vec, buf, "genbank")
        sc._save_library([{"id": "pVec", "name": "pVec",
                            "gb_text": buf.getvalue(),
                            "size": len(vec.seq)}])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from textual.widgets import Button
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            pane = await self._setup_lane(pilot, modal, app)
            modal.query_one("#btn-trad-simulate", Button).press()
            await pilot.pause(); await pilot.pause(0.2)
            save_fwd = modal.query_one("#btn-trad-save-fwd", Button)
            save_rev = modal.query_one("#btn-trad-save-rev", Button)
            # Exactly one of them is enabled (directional cloning).
            assert save_fwd.disabled is not save_rev.disabled
            (save_fwd if not save_fwd.disabled else save_rev).press()
            await pilot.pause()
            name_modal = app.screen
            assert isinstance(name_modal, sc.NamePlasmidModal), (
                f"expected NamePlasmidModal; "
                f"got {type(name_modal).__name__}"
            )
            chosen_name = "my-trad-product"
            name_modal.dismiss(chosen_name)
            await pilot.pause(); await pilot.pause(0.05)
            after = sc._load_library()
            saved_entries = [
                e for e in after
                if isinstance(e, dict) and e.get("name") == chosen_name
            ]
            assert len(saved_entries) == 1
            saved = saved_entries[0]
            assert saved.get("size", 0) > 0
            assert "gb_text" in saved
            assert "LOCUS" in saved["gb_text"]
            assert str(saved.get("source", "")).startswith("traditional:")
            assert save_fwd.disabled and save_rev.disabled
            assert pane._fwd_product is None
            assert pane._rev_product is None

    async def test_traditional_pane_save_records_history_xml(
            self, tiny_record, isolated_library):
        """Regression guard for 2026-05-06 Phase 4b wiring:
        a successful Traditional cloning save attaches a CommercialSaaS-
        compatible `<HistoryTree>` to the new library entry. The
        XML must parse via `_parse_commercialsaas_history`, name the new
        plasmid as the top node, and link both source plasmids as
        parent fragments."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        # Set up a vector with EcoRI + BamHI flanking a stuffer.
        vec = SeqRecord(
            Seq("AAAAAAAAGAATTCAAAAAAAAAAGGATCCTTTTTTTT"),
            id="pVec", name="pVec",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        buf = io.StringIO()
        SeqIO.write(vec, buf, "genbank")
        sc._save_library([{"id": "pVec", "name": "pVec",
                            "gb_text": buf.getvalue(),
                            "size": len(vec.seq)}])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from textual.widgets import Button
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            await self._setup_lane(pilot, modal, app)
            modal.query_one("#btn-trad-simulate", Button).press()
            await pilot.pause(); await pilot.pause(0.2)
            save_fwd = modal.query_one("#btn-trad-save-fwd", Button)
            save_rev = modal.query_one("#btn-trad-save-rev", Button)
            (save_fwd if not save_fwd.disabled else save_rev).press()
            await pilot.pause()
            name_modal = app.screen
            assert isinstance(name_modal, sc.NamePlasmidModal)
            name_modal.dismiss("my-trad-history-test")
            await pilot.pause(); await pilot.pause(0.05)
            saved = next(
                (e for e in sc._load_library()
                  if isinstance(e, dict)
                  and e.get("name") == "my-trad-history-test"),
                None,
            )
            assert saved is not None
            assert "history_xml" in saved, (
                f"history not attached to entry {saved.get('name')!r}; "
                f"keys: {sorted(saved)}")
            # Parse it — root should name the new product, with two
            # parent nodes (the synthesised PCR insert + the pVec
            # vector).
            root = sc._parse_commercialsaas_history(saved["history_xml"])
            assert root is not None
            assert saved["name"] in root.name
            assert root.operation == "insertFragment"
            parent_names = [p.name for p in root.parents]
            assert any("PCR-product" in pn or "myInsert" in pn
                        for pn in parent_names), parent_names
            assert any("pVec" in pn for pn in parent_names), parent_names
            # Regenerated sites carry both enzymes used.
            site_names = [s["name"] for s in root.regenerated_sites]
            assert {"EcoRI", "BamHI"} <= set(site_names)
            # Input summary marks which orientation was saved.
            ops = [s["manipulation"] for s in root.input_summaries]
            assert any(op in ("ligateFwd", "ligateRev") for op in ops)

    async def test_save_buttons_redisable_on_input_change(
            self, tiny_record, isolated_library):
        """Regression guard for 2026-05-05 fix: after a successful
        Simulate the Save buttons enable, but if the user then
        changes any input (enzyme, source row, mode, PCR name) the
        cached product is stale — Save must re-disable until the
        user re-Simulates."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        vec = SeqRecord(
            Seq("AAAAAAAAGAATTCAAAAAAAAAAGGATCCTTTTTTTT"),
            id="pVec", name="pVec",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        buf = io.StringIO()
        SeqIO.write(vec, buf, "genbank")
        sc._save_library([{"id": "pVec", "name": "pVec",
                            "gb_text": buf.getvalue(),
                            "size": len(vec.seq)}])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from textual.widgets import Button, Select, DataTable
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            await self._setup_lane(pilot, modal, app)
            modal.query_one("#btn-trad-simulate", Button).press()
            await pilot.pause(); await pilot.pause(0.2)
            save_fwd = modal.query_one("#btn-trad-save-fwd", Button)
            save_rev = modal.query_one("#btn-trad-save-rev", Button)
            assert not (save_fwd.disabled and save_rev.disabled)
            # Park cursor on the backbone row (idx 0) so the edit
            # panel mirrors the backbone's enzymes, then change one —
            # both Save buttons must re-disable as the cached product
            # is now stale.
            lt = modal.query_one("#trad-lane", DataTable)
            lt.move_cursor(row=0)
            await pilot.pause()
            modal.query_one(
                "#trad-edit-enz-2", Select,
            ).value = "HindIII"
            await pilot.pause()
            assert save_fwd.disabled and save_rev.disabled

    async def test_traditional_pane_simulate_without_source_errors(
            self, tiny_record, isolated_library):
        """Hitting Simulate in plasmid mode without picking an insert
        plasmid surfaces a clear error rather than crashing. Library
        starts empty (the autouse fixture wipes it) so the source
        DataTable has no rows."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            from textual.widgets import TabbedContent, Button
            modal.query_one("#ctor-tabs", TabbedContent).active = \
                "ctor-tab-traditional"
            await pilot.pause()
            # Click Simulate in default (plasmid) mode without any
            # plasmid in the library — should set a red error message
            # and leave the cached products empty.
            modal.query_one("#btn-trad-simulate", Button).press()
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane",
                                     sc.TraditionalCloningPane)
            # No products should be cached on a failed simulate.
            assert pane._fwd_product is None
            assert pane._rev_product is None

    async def test_lane_add_remove_reorder(
            self, tiny_record, isolated_library):
        """Lane-management buttons: Add to Lane queues a fragment; ↑/↓
        reorder via cursor; ✕ Remove drops the row under the cursor.
        Empty lane → Simulate fails with a clear message rather than
        crashing.
        Regression for the 2026-05-23 lane rebuild."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        vec = SeqRecord(
            Seq("AAAAAAAAGAATTCAAAAAAAAAAGGATCCTTTTTTTT"),
            id="pVec", name="pVec",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        buf = io.StringIO()
        SeqIO.write(vec, buf, "genbank")
        sc._save_library([{"id": "pVec", "name": "pVec",
                            "gb_text": buf.getvalue(),
                            "size": len(vec.seq)}])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            from textual.widgets import (TabbedContent, RadioButton,
                                          Input, TextArea, Select,
                                          DataTable, Button, Static)
            modal.query_one("#ctor-tabs", TabbedContent).active = \
                "ctor-tab-traditional"
            await pilot.pause()
            await pilot.pause(0.05)
            modal.query_one("#trad-mode-pcr", RadioButton).value = True
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane",
                                     sc.TraditionalCloningPane)
            # Empty lane + Simulate → red error message.
            modal.query_one("#btn-trad-simulate", Button).press()
            await pilot.pause()
            results = str(
                modal.query_one("#trad-results-text", Static).content
            )
            assert "Lane is empty" in results, results
            # Queue three PCR fragments with different enzyme pairs so
            # reorder + remove can be distinguished by lane row state.
            # `await pilot.pause()` between adds is required: Button
            # .press() posts a message handled asynchronously, so
            # without a pause the next iteration overwrites the
            # Input value before the previous handler captures it
            # and all three lane rows end up identical to the last.
            async def _add(name, e1, e2, seq):
                modal.query_one(
                    "#trad-pcr-name", Input,
                ).value = name
                modal.query_one(
                    "#trad-pcr-seq",  TextArea,
                ).text = seq
                modal.query_one(
                    "#btn-trad-add-frag", Button,
                ).press()
                await pilot.pause()
                # Enzymes set in the master/detail editor after Add
                # to Lane (the new flow defers enzyme choice until the
                # donor row exists in the lane).
                modal.query_one(
                    "#trad-edit-enz-1", Select,
                ).value = e1
                modal.query_one(
                    "#trad-edit-enz-2", Select,
                ).value = e2
                await pilot.pause()
            await _add("alpha", "EcoRI",   "BamHI",   "AAACCC")
            await _add("beta",  "BamHI",   "SalI",    "GGGTTT")
            await _add("gamma", "SalI",    "HindIII", "CCCAAA")
            assert [s["name"] for s in pane._lane_inserts] == \
                ["alpha", "beta", "gamma"]
            # Reorder: move middle row (beta) up.
            lt = modal.query_one("#trad-lane", DataTable)
            lt.move_cursor(row=1)
            await pilot.pause()
            modal.query_one("#btn-trad-lane-up", Button).press()
            await pilot.pause()
            assert [s["name"] for s in pane._lane_inserts] == \
                ["beta", "alpha", "gamma"]
            # Remove the cursor row (now beta at row 0).
            lt.move_cursor(row=0)
            await pilot.pause()
            modal.query_one("#btn-trad-lane-remove", Button).press()
            await pilot.pause()
            assert [s["name"] for s in pane._lane_inserts] == \
                ["alpha", "gamma"]
            # Clear Lane drops everything.
            modal.query_one("#btn-trad-lane-clear", Button).press()
            await pilot.pause()
            assert pane._lane_inserts == []

    async def test_save_pushes_name_modal_with_default(
            self, tiny_record, isolated_library):
        """Clicking Save Forward after a successful Simulate must
        push ``NamePlasmidModal`` with an auto-default name derived
        from ``{vector} · {fragment names} (suffix)`` — same shape
        Golden Braid uses. Cancelling the modal is a no-op (no
        library entry added)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        vec = SeqRecord(
            Seq("AAAAAAAAGAATTCAAAAAAAAAAGGATCCTTTTTTTT"),
            id="pVec", name="pVec",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        buf = io.StringIO()
        SeqIO.write(vec, buf, "genbank")
        sc._save_library([{"id": "pVec", "name": "pVec",
                            "gb_text": buf.getvalue(),
                            "size": len(vec.seq)}])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from textual.widgets import Button
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause(); await pilot.pause(0.05)
            await self._setup_lane(pilot, modal, app)
            modal.query_one("#btn-trad-simulate", Button).press()
            await pilot.pause(); await pilot.pause(0.2)
            save_fwd = modal.query_one("#btn-trad-save-fwd", Button)
            save_rev = modal.query_one("#btn-trad-save-rev", Button)
            (save_fwd if not save_fwd.disabled else save_rev).press()
            await pilot.pause()
            name_modal = app.screen
            assert isinstance(name_modal, sc.NamePlasmidModal), (
                f"expected NamePlasmidModal; "
                f"got {type(name_modal).__name__}"
            )
            default = name_modal._default_name
            assert "pVec" in default, default
            assert "myInsert" in default, default
            assert "fwd" in default or "rev" in default, default
            lib_before = len(sc._load_library())
            name_modal.dismiss(None)
            await pilot.pause(); await pilot.pause(0.05)
            assert len(sc._load_library()) == lib_before

    async def test_master_detail_editor_populates_on_lane_select(
            self, tiny_record, isolated_library):
        """After Add to Lane, the master/detail editor auto-parks on
        the new row. Setting enzymes in the editor flows into the
        lane row's spec and the lane DataTable's E1/E2 columns update
        accordingly. Regression for the 2026-05-23 lane rebuild's
        master/detail wiring."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        vec = SeqRecord(
            Seq("AAAAAAAAGAATTCAAAAAAAAAAGGATCCTTTTTTTT"),
            id="pVec", name="pVec",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        buf = io.StringIO()
        SeqIO.write(vec, buf, "genbank")
        sc._save_library([{"id": "pVec", "name": "pVec",
                            "gb_text": buf.getvalue(),
                            "size": len(vec.seq)}])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            from textual.widgets import (TabbedContent, RadioButton,
                                          Input, TextArea, Select,
                                          Button)
            modal.query_one("#ctor-tabs", TabbedContent).active = \
                "ctor-tab-traditional"
            await pilot.pause(); await pilot.pause(0.05)
            modal.query_one("#trad-mode-pcr", RadioButton).value = True
            await pilot.pause()
            modal.query_one("#trad-pcr-name", Input).value = "myFrag"
            modal.query_one("#trad-pcr-seq",  TextArea).text = (
                "GAGCATGAAACGGCCAAGTAA")
            modal.query_one("#btn-trad-add-frag", Button).press()
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane",
                                     sc.TraditionalCloningPane)
            assert pane._edit_row_idx == 0, (
                "edit panel should auto-park on the new row"
            )
            assert pane._lane_inserts[0]["enz_left"] == ""
            modal.query_one("#trad-edit-enz-1", Select).value = "EcoRI"
            await pilot.pause()
            modal.query_one("#trad-edit-enz-2", Select).value = "BamHI"
            await pilot.pause()
            assert pane._lane_inserts[0]["enz_left"]  == "EcoRI"
            assert pane._lane_inserts[0]["enz_right"] == "BamHI"

    async def test_donor_frag_radio_overrides_auto_pick(
            self, tiny_record, isolated_library):
        """For a plasmid-mode donor with valid 2-fragment digest, the
        donor-fragment radio toggles ``donor_frag_idx`` between 0 and
        1, overriding the feature-aware auto-pick used when the user
        hasn't touched the radio. Verifies that the Simulate path
        consumes the override (different ligated sequence depending
        on which fragment is picked)."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        # Donor plasmid with EcoRI + BamHI sites producing two clearly
        # different-sized fragments.
        donor = SeqRecord(
            Seq("AAAAAAAAGAATTCTTTTTTTTTTTTTTTTGGATCCCCCCCCCC"),
            id="pDonor", name="pDonor",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        vec = SeqRecord(
            Seq("ACGTACGTGAATTCAAAAAAAAAAAAAAGGATCCACGTACGT"),
            id="pVec", name="pVec",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        gb = {}
        for r in (donor, vec):
            b = io.StringIO()
            SeqIO.write(r, b, "genbank")
            gb[r.id] = b.getvalue()
        sc._save_library([
            {"id": "pDonor", "name": "pDonor",
              "gb_text": gb["pDonor"], "size": len(donor.seq)},
            {"id": "pVec",   "name": "pVec",
              "gb_text": gb["pVec"],   "size": len(vec.seq)},
        ])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            from textual.widgets import (TabbedContent, Select,
                                          DataTable, Button, RadioButton)
            modal.query_one("#ctor-tabs", TabbedContent).active = \
                "ctor-tab-traditional"
            await pilot.pause(); await pilot.pause(0.05)
            entries = sorted(
                (e for e in sc._load_library() if isinstance(e, dict)),
                key=lambda e: sc._natural_sort_key(
                    e.get("name") or e.get("id") or ""
                ),
            )
            donor_idx = next(i for i, e in enumerate(entries)
                              if e.get("name") == "pDonor")
            # Pick donor + Add.
            modal.query_one(
                "#trad-source-table", DataTable,
            ).move_cursor(row=donor_idx)
            await pilot.pause()
            modal.query_one("#btn-trad-add-frag", Button).press()
            await pilot.pause()
            # Set enzymes via edit panel.
            modal.query_one(
                "#trad-edit-enz-1", Select,
            ).value = "EcoRI"
            await pilot.pause()
            modal.query_one(
                "#trad-edit-enz-2", Select,
            ).value = "BamHI"
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane",
                                     sc.TraditionalCloningPane)
            assert pane._edit_frags is not None, (
                "digest should populate _edit_frags for plasmid mode"
            )
            assert len(pane._edit_frags) == 2
            # Toggle to fragment B (idx=1) — overrides auto-pick.
            modal.query_one(
                "#trad-edit-frag-1", RadioButton,
            ).value = True
            await pilot.pause()
            assert pane._lane_inserts[0]["donor_frag_idx"] == 1, (
                f"radio toggle should set donor_frag_idx=1; "
                f"got {pane._lane_inserts[0]['donor_frag_idx']}"
            )
            # Toggle to fragment A (idx=0).
            modal.query_one(
                "#trad-edit-frag-0", RadioButton,
            ).value = True
            await pilot.pause()
            assert pane._lane_inserts[0]["donor_frag_idx"] == 0

    async def test_digest_cache_avoids_redundant_work(
            self, tiny_record, isolated_library):
        """`_cached_digest` memoises `(entry_id, enzymes)` digests so
        repeat lookups (the gel preview's per-event re-render) don't
        re-scan the source plasmid. Verifies cache hit on second call
        with the same enzymes."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        donor = SeqRecord(
            Seq("AAAAAAAAGAATTCTTTTTTTTTTTTTTTTGGATCCCCCCCCCC"),
            id="pDonor", name="pDonor",
            annotations={"molecule_type": "DNA",
                          "topology": "circular"},
        )
        b = io.StringIO()
        SeqIO.write(donor, b, "genbank")
        sc._save_library([{"id": "pDonor", "name": "pDonor",
                            "gb_text": b.getvalue(),
                            "size": len(donor.seq)}])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            from textual.widgets import TabbedContent
            modal.query_one("#ctor-tabs", TabbedContent).active = \
                "ctor-tab-traditional"
            await pilot.pause(); await pilot.pause(0.05)
            pane = modal.query_one("#ctor-trad-pane",
                                     sc.TraditionalCloningPane)
            frags_1, err_1 = pane._cached_digest(
                "pDonor", ["EcoRI", "BamHI"],
            )
            assert err_1 == "", err_1
            assert frags_1 is not None and len(frags_1) == 2
            cache_key = ("pDonor", tuple(sorted(["EcoRI", "BamHI"])))
            assert cache_key in pane._digest_cache
            # Second call with same enzymes → cache hit (same list
            # identity, no re-digest).
            frags_2, err_2 = pane._cached_digest(
                "pDonor", ["EcoRI", "BamHI"],
            )
            assert err_2 == ""
            assert frags_2 is frags_1, "second call should hit cache"
            # Different enzyme combo → cache miss → new entry.
            _, err_3 = pane._cached_digest(
                "pDonor", ["EcoRI", "BamHI", "HindIII"],
            )
            # 3+ cuts on a circular plasmid surfaces as an error
            # (ambiguous excise) — still a valid memo lookup, just
            # err returned without a frags list.
            assert err_3 != "" or len(pane._digest_cache) >= 2

    async def test_gel_preview_renders_for_donor_lane(
            self, tiny_record, isolated_library):
        """Gel preview shows non-empty text once at least one donor
        is queued. Pristine state shows the empty-state hint."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        import io
        donor = SeqRecord(
            Seq("AAAAAAAAGAATTCTTTTTTTTTTTTTTTTGGATCCCCCCCCCC"),
            id="pDonor", name="pDonor",
            annotations={"molecule_type": "DNA", "topology": "circular"},
        )
        b = io.StringIO()
        SeqIO.write(donor, b, "genbank")
        sc._save_library([{"id": "pDonor", "name": "pDonor",
                            "gb_text": b.getvalue(),
                            "size": len(donor.seq)}])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause(); await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            from textual.widgets import (TabbedContent, Static,
                                          DataTable, Button, Select)
            modal.query_one("#ctor-tabs", TabbedContent).active = \
                "ctor-tab-traditional"
            await pilot.pause(); await pilot.pause(0.05)
            gel = modal.query_one("#trad-gel-text", Static)
            # Empty state: no donors queued — gel shows hint.
            assert "empty lane" in str(gel.content)
            # Add a donor + set enzymes.
            entries = sorted(
                (e for e in sc._load_library() if isinstance(e, dict)),
                key=lambda e: sc._natural_sort_key(
                    e.get("name") or e.get("id") or ""
                ),
            )
            donor_idx = next(i for i, e in enumerate(entries)
                              if e.get("name") == "pDonor")
            modal.query_one(
                "#trad-source-table", DataTable,
            ).move_cursor(row=donor_idx)
            await pilot.pause()
            modal.query_one("#btn-trad-add-frag", Button).press()
            await pilot.pause()
            modal.query_one("#trad-edit-enz-1", Select).value = "EcoRI"
            await pilot.pause()
            modal.query_one("#trad-edit-enz-2", Select).value = "BamHI"
            await pilot.pause()
            # Gel now non-empty + carries the lane label (truncated
            # to lane_width=7 chars + the backbone "[B]" prefix that
            # `_refresh_gel_preview` adds to backbone rows).
            content = str(modal.query_one(
                "#trad-gel-text", Static,
            ).content)
            assert "empty lane" not in content
            # Lane label is "[B]pDonor" truncated to 7 chars → "[B]pDon".
            assert "pDon" in content, (
                f"donor lane label missing from gel render; "
                f"content head: {content[:200]!r}"
            )

    def test_features_carry_through_with_correct_offsets(self):
        # Insert with one feature at positions 2..5 (relative to insert).
        insert = sc._make_synthetic_fragment(
            "AAACCCGGG", enz_left="EcoRI", enz_right="BamHI",
            features=[{"start": 2, "end": 5, "label": "in-feat", "strand": 1}],
        )
        vector = {
            "top_seq": "T" * 20,
            "left":  {"overhang_seq": "GATC", "kind": "5'", "enzyme": "BamHI"},
            "right": {"overhang_seq": "AATT", "kind": "5'", "enzyme": "EcoRI"},
            "features": [{"start": 5, "end": 10, "label": "vec-feat",
                            "strand": 1}],
            "source_label": "vector",
        }
        result = sc._simulate_traditional_cloning(insert, vector)
        fwd_feats = result["forward"]["features"]
        # Vector feat at 5..10 stays put; insert feat at 2..5 shifts by 20.
        labels = {f["label"]: (f["start"], f["end"]) for f in fwd_feats}
        assert labels["vec-feat"] == (5, 10)
        assert labels["in-feat"]  == (22, 25)


class TestSimulateTraditionalCloningMulti:
    """N-way ligation wrapper. ``_simulate_traditional_cloning_multi``
    pre-chains 1..N insert fragments in lane order, then delegates to
    the 2-fragment engine for the final vector + chained-insert
    ligation."""

    def _vec(self, oh_left: str, oh_right: str,
              enz_left: str = "EcoRI", enz_right: str = "EcoRI"):
        return {
            "top_seq": "TGGCCCC" * 10,
            "left":  {"overhang_seq": oh_left,  "kind": "5'",
                       "enzyme": enz_left},
            "right": {"overhang_seq": oh_right, "kind": "5'",
                       "enzyme": enz_right},
            "features": [], "source_label": "vec",
        }

    def test_n1_matches_two_fragment_engine(self):
        """N=1 chemistry (top_seq + compatibility) must be identical
        between the 2-fragment engine and the multi-wrapper. The
        multi-wrapper adds scar-detection warnings + misc_feature
        annotations on top, so the test compares the chemistry-
        relevant fields rather than full-dict equality."""
        insert = sc._make_synthetic_fragment(
            "GAGCATGAAACGGCCAAGTAA",
            enz_left="EcoRI", enz_right="BamHI",
            source_label="insert",
        )
        vector = self._vec("GATC", "AATT",
                            enz_left="BamHI", enz_right="EcoRI")
        single = sc._simulate_traditional_cloning(insert, vector)
        multi  = sc._simulate_traditional_cloning_multi([insert], vector)
        for orient in ("forward", "reverse"):
            assert single[orient]["top_seq"] == multi[orient]["top_seq"]
            assert single[orient]["compatible"] == \
                multi[orient]["compatible"]
        assert single["errors"] == multi["errors"]

    def test_empty_lane_errors(self):
        """An empty insert list short-circuits with a clear error
        rather than nuking the vector — defensive guard for the UI
        path where the user clicks Simulate before adding fragments."""
        vector = self._vec("AATT", "AATT")
        result = sc._simulate_traditional_cloning_multi([], vector)
        assert result["forward"]["compatible"] is False
        assert any("No insert" in e for e in result["errors"])

    def test_three_way_chain_ligates(self):
        """Three inserts chained by matching sticky ends + a vector
        whose outer ends match the chain's flanks → forward
        orientation ligates; the chained insert lands at the right
        offset on the product."""
        # Chain: insert1 (E..B) → insert2 (B..S) → insert3 (S..H)
        # Sticky ends: BamHI=GATC, SalI=TCGA, HindIII=AGCT, EcoRI=AATT
        i1 = sc._make_synthetic_fragment(
            "AAAA", enz_left="EcoRI",   enz_right="BamHI",
            source_label="i1",
        )
        i2 = sc._make_synthetic_fragment(
            "CCCC", enz_left="BamHI",   enz_right="SalI",
            source_label="i2",
        )
        i3 = sc._make_synthetic_fragment(
            "GGGG", enz_left="SalI",    enz_right="HindIII",
            source_label="i3",
        )
        # Vector outer ends: left=AGCT (HindIII), right=AATT (EcoRI)
        # so the chain's outer flanks (EcoRI on i1.left, HindIII on
        # i3.right) ligate to the vector's matching ends.
        vector = self._vec(
            "AGCT", "AATT",
            enz_left="HindIII", enz_right="EcoRI",
        )
        result = sc._simulate_traditional_cloning_multi(
            [i1, i2, i3], vector,
        )
        assert result["errors"] == [], (
            f"unexpected errors: {result['errors']}"
        )
        assert result["forward"]["compatible"] is True
        # Sequence carries vector + chained inserts (no gaps); top
        # strand should at minimum include each insert body verbatim.
        fwd_seq = result["forward"]["top_seq"]
        for body in ("AAAA", "CCCC", "GGGG"):
            assert body in fwd_seq, f"missing {body!r} in {fwd_seq!r}"

    def test_biobrick_spei_xbai_junction_is_scar(self):
        """SpeI (A^CTAGT) + XbaI (T^CTAGA) produce the same CTAG 5'
        overhang so they ligate, but the resulting junction sequence
        is ACTAGA — neither SpeI (ACTAGT) nor XbaI (TCTAGA). That's
        the iGEM BioBrick idempotent property: ligated joints are
        uncuttable by either parent enzyme. `_classify_junction`
        must detect this and surface a warning so the user knows
        the joint is irreversible. The scar must also appear as a
        misc_feature on the saved product."""
        # SpeI cut leaves CTAG overhang; XbaI cut leaves CTAG
        # overhang. Use synthetic fragments stamped with these
        # enzyme overhangs so the engine sees ACTAGA at the joint.
        insert = sc._make_synthetic_fragment(
            "GAGCATG", enz_left="SpeI", enz_right="XbaI",
            source_label="biobrick-part",
        )
        vector = {
            "top_seq": "TTTTTTTT" * 5,
            "left":  {"overhang_seq": "CTAG", "kind": "5'",
                       "enzyme": "XbaI"},
            "right": {"overhang_seq": "CTAG", "kind": "5'",
                       "enzyme": "SpeI"},
            "features": [], "source_label": "biobrick-vec",
        }
        result = sc._simulate_traditional_cloning_multi(
            [insert], vector,
        )
        assert result["forward"]["compatible"] is True
        # Warnings include at least one "scar" notice — the SpeI/
        # XbaI joint is uncuttable.
        joined = " ".join(result["warnings"])
        assert "scar" in joined.lower(), (
            f"expected SpeI/XbaI scar warning; got: {result['warnings']!r}"
        )
        # The product is NO LONGER annotated with a "LIGATION SCAR" feature
        # (user request 2026-06-09: leave the scars as-is in the sequence,
        # not annotated). The scar CLASSIFICATION still rides the warnings
        # above; the junction itself is tagged as a light-blue, arrowless
        # 4 bp overhang `misc_feature` instead.
        feats = result["forward"]["features"]
        assert not any("SCAR" in str(f.get("label", "")).upper()
                       for f in feats), (
            f"LIGATION SCAR feature should be gone; features: {feats!r}"
        )
        overhang_feats = [
            f for f in feats
            if f.get("type") == "misc_feature"
            and f.get("strand") == 0
            and f.get("color") == "#ADD8E6"
        ]
        assert overhang_feats, (
            f"expected a light-blue arrowless overhang tag; "
            f"features: {feats!r}"
        )

    def test_origin_junction_overhang_is_4bp_wrap(self):
        """The closing junction at the origin tags a FULL 4 bp overhang as a
        wrap feature (end < start → CompoundLocation on save), not the 2 bp
        head a flat [0,2) clamp gives (adversarial review F6)."""
        vec = {"top_seq": "AAAATTTT" * 6, "left": {"enzyme": "EcoRI"},
               "right": {"enzyme": "BamHI"}}
        ins = [{"top_seq": "C" * 24, "left": {"enzyme": "BamHI"},
                "right": {"enzyme": "EcoRI"}}]
        prod = vec["top_seq"] + ins[0]["top_seq"]
        res = {"forward": {"top_seq": prod, "features": [], "compatible": True},
               "reverse": {"top_seq": prod, "features": [], "compatible": True},
               "warnings": []}
        sc._annotate_scars_on_product(res, ins, vec)
        wrap = [f for f in res["forward"]["features"] if f["start"] > f["end"]]
        assert wrap, f"no origin wrap overhang: {res['forward']['features']!r}"
        assert len(wrap[0]["label"]) == 4
        assert wrap[0]["color"] == "#ADD8E6" and wrap[0]["strand"] == 0

    def test_pcr_insert_carries_features_to_product(self):
        """A PCR / Clone-region insert's own features must reach the cloned
        product. `_build_insert_from_pcr` used to drop them, so the insert
        ligated in as a featureless black box (a cloned TU lost its features)."""
        pad, s5, s3 = "GCGC", "GAATTC", "GGATCC"        # EcoRI / BamHI
        insert = "ATGAAACGT" + "ACTGCATGCAGTACGTAGCT" * 4
        amplicon = pad + s5 + insert + s3 + sc._rc(pad)
        lead = len(pad) + len(s5)
        ifeat = {"start": lead + 12, "end": lead + 60, "type": "CDS",
                 "label": "MyGene", "strand": 1}
        ins_frag = sc._make_synthetic_fragment(
            amplicon, enz_left="EcoRI", enz_right="BamHI", features=[ifeat])
        assert any(f["label"] == "MyGene" for f in ins_frag["features"])
        vbody = "TTACGGATCAGCTAGGCATTAGC" * 6
        vfrag = sc._make_synthetic_fragment(
            s3 + vbody + s5, enz_left="BamHI", enz_right="EcoRI",
            source_label="vec")
        res = sc._simulate_traditional_cloning(ins_frag, vfrag)
        gene = [f for f in res["forward"]["features"]
                if f.get("label") == "MyGene"]
        assert gene, ("insert feature dropped from product: "
                      f"{[f.get('label') for f in res['forward']['features']]}")
        assert gene[0]["end"] - gene[0]["start"] == 48   # length preserved

    def test_disrupted_split_feature_labeling(self):
        """`_label_disrupted_split_features` suffixes "(disrupted)" + a note on
        cut-split (`_split`) halves only, names the enzyme(s), and is
        idempotent — leaving un-split (carried insert) features alone."""
        feats = [
            {"label": "lacZ", "type": "CDS", "_split": "head", "start": 0, "end": 8},
            {"label": "lacZ", "type": "CDS", "_split": "tail", "start": 50, "end": 58},
            {"label": "AmpR", "type": "CDS", "start": 20, "end": 40},   # not split
        ]
        sc._label_disrupted_split_features(feats, ["BamHI", "EcoRI"])
        assert feats[0]["label"] == "lacZ (disrupted)"
        assert feats[1]["label"] == "lacZ (disrupted)"
        assert "BamHI/EcoRI" in feats[0]["note"]        # sorted enzyme set
        assert feats[2]["label"] == "AmpR"              # un-split untouched
        assert "note" not in feats[2]
        # Idempotent — a second pass doesn't double-suffix.
        sc._label_disrupted_split_features(feats, ["EcoRI"])
        assert feats[0]["label"] == "lacZ (disrupted)"

    def test_simulate_disrupts_vector_feature_at_cut(self):
        """A vector feature the cut sites sit inside (e.g. lacZα across the MCS)
        is split by the digest (`_split` head/tail) and must surface in the
        cloned product as "(disrupted)" with a note naming the cut enzyme —
        the blue/white-screening knock-out you'd see cloning into a lacZ MCS."""
        pad, s5, s3 = "GCGC", "GAATTC", "GGATCC"        # EcoRI / BamHI
        insert = "ATGAAACGT" + "ACTGCATGCAGTACGTAGCT" * 3
        ins = sc._make_synthetic_fragment(
            pad + s5 + insert + s3 + sc._rc(pad),
            enz_left="EcoRI", enz_right="BamHI")
        vbody = "TTACGGATCAGCTAGGCATTAGC" * 6
        # Vector digested at the MCS inside lacZα → the backbone fragment
        # carries a split half of lacZα (what `_split_features_at_cuts` emits).
        lacz_head = {"start": 0, "end": 8, "type": "CDS", "label": "lacZα",
                     "strand": 1, "_split": "head"}
        vfrag = sc._make_synthetic_fragment(
            s3 + vbody + s5, enz_left="BamHI", enz_right="EcoRI",
            source_label="vec", features=[lacz_head])
        res = sc._simulate_traditional_cloning(ins, vfrag)
        lac = [f for f in res["forward"]["features"]
               if "lacZ" in str(f.get("label"))]
        assert lac, "vector lacZα feature lost from product"
        assert all("(disrupted)" in f["label"] for f in lac)
        assert any("EcoRI" in (f.get("note") or "")
                   or "BamHI" in (f.get("note") or "") for f in lac)

    def test_disrupts_feature_with_excised_middle(self):
        """The classic lacZ-MCS case: a feature whose MIDDLE is excised by two
        cuts that both sit inside it — so both ends survive in the same
        (backbone) fragment — is still flagged disrupted, not just features
        split across two fragments. Regression for the gap found verifying on
        real pUC19 (cloning into lacZα's MCS used to leave lacZα looking
        intact)."""
        feats = [{"start": 10, "end": 60, "strand": 1,
                  "type": "CDS", "label": "lacZ"}]
        # Circular vector, cuts at 25 + 40 (both strictly inside [10, 60)).
        slots = sc._split_features_at_cuts(feats, 100, [25, 40], circular=True)
        allf = [f for fs in slots.values() for f in fs]
        lac = [f for f in allf if f.get("label") == "lacZ"]
        assert lac and any(f.get("_split") == "whole" for f in lac), (
            "lacZ spanning an excised MCS not flagged: "
            f"{[f.get('_split') for f in lac]}")
        sc._label_disrupted_split_features(allf, ["EcoRI"])
        assert all("(disrupted)" in f["label"] for f in lac)

    def test_excise_inside_feature_does_not_flag_neighbours(self):
        """An excise INSIDE one feature must not flag a neighbouring feature
        that merely shares the backbone fragment but has no cut inside it."""
        feats = [
            {"start": 10, "end": 60,  "strand": 1, "type": "CDS", "label": "lacZ"},
            {"start": 70, "end": 95,  "strand": 1, "type": "CDS", "label": "keep"},
        ]
        slots = sc._split_features_at_cuts(feats, 100, [25, 40], circular=True)
        allf = [f for fs in slots.values() for f in fs]
        keep = [f for f in allf if f.get("label") == "keep"]
        assert keep and all(f.get("_split") is None for f in keep)

    def test_internal_junction_mismatch_surfaces_pair(self):
        """If two adjacent inserts have incompatible sticky ends, the
        error message names the failing pair by source_label so the
        user can fix that specific junction."""
        i1 = sc._make_synthetic_fragment(
            "AAA", enz_left="EcoRI", enz_right="BamHI",
            source_label="i1",
        )
        # i2's LEFT is SalI (TCGA), but i1's RIGHT is BamHI (GATC)
        # — mismatched. Junction 1 → 2 should fail.
        i2 = sc._make_synthetic_fragment(
            "CCC", enz_left="SalI", enz_right="EcoRI",
            source_label="i2",
        )
        vector = self._vec("AATT", "AATT")
        result = sc._simulate_traditional_cloning_multi([i1, i2], vector)
        assert result["forward"]["compatible"] is False
        joined_err = " ".join(result["errors"])
        assert "Junction 1 → 2" in joined_err
        assert "'i1'" in joined_err and "'i2'" in joined_err


# ═══════════════════════════════════════════════════════════════════════════════
# Constructor modal — multi-grammar tabs (2026-05-07)
# ═══════════════════════════════════════════════════════════════════════════════
#
# The Constructor was refactored from a single Modular tab (gb_l0
# only) into per-grammar tabs (Golden Braid + MoClo Plant) plus the
# pre-existing Traditional tab. Each modular tab pulls its parts
# palette from `parts_bin.json` filtered by the tab's grammar id
# when the `constructor_filter_by_grammar` setting is on (default).
# These tests cover the wiring: tabs exist, palette filters, palette
# reflects user parts, per-grammar entry vector is independent.

import pytest as _ctor_pytest


@_ctor_pytest.fixture
def isolated_parts_bin(tmp_path, monkeypatch):
    """Redirect `_PARTS_BIN_FILE` to a tmp path so the Constructor
    palette tests don't touch the real parts_bin.json. Mirrors the
    fixture in test_domesticator.py."""
    tmp_bin = tmp_path / "parts_bin.json"
    monkeypatch.setattr(sc._state, "_PARTS_BIN_FILE", tmp_bin)
    monkeypatch.setattr(sc._state, "_parts_bin_cache", None)
    return tmp_bin


class TestConstructorMultiGrammarTabs:
    """The Constructor exposes Traditional + Golden Braid + MoClo
    Plant tabs. Each modular tab uses its own grammar context for
    palette filtering, validation, and entry-vector banner."""

    def test_palette_helper_excludes_builtin_catalog(self, isolated_parts_bin):
        """Regression guard for 2026-05-07: the palette must NOT
        include the built-in `_GB_L0_PARTS` catalog rows. Those are
        placeholder entries with no real sequence — they can't
        actually assemble, so showing them as palette options is
        misleading. An empty parts bin → an empty palette is the
        honest state."""
        # Empty parts bin → empty palette under either filter setting.
        sc._save_parts_bin([])
        for filter_enabled in (True, False):
            rows = sc._palette_rows_for_grammar(
                "gb_l0", filter_enabled=filter_enabled,
            )
            assert rows == [], (
                f"Empty parts bin must produce empty palette "
                f"(filter={filter_enabled}); got {len(rows)} rows"
            )
        # Even with the gb_l0 grammar's catalog populated, none of
        # the catalog names should leak into the palette.
        names_with_no_user_parts = {
            r[0] for r in sc._palette_rows_for_grammar(
                "gb_l0", filter_enabled=True,
            )
        }
        catalog_names = {row[0] for row in sc._GB_L0_PARTS}
        assert not (names_with_no_user_parts & catalog_names), (
            "Built-in catalog names leaked into the palette: "
            f"{names_with_no_user_parts & catalog_names}"
        )

    def test_palette_helper_filters_by_grammar(self, isolated_parts_bin):
        # Save two parts with different grammar ids; the helper
        # filtered by grammar should only return the matching one.
        sc._save_parts_bin([
            {
                "name": "myProm", "type": "Promoter", "position": "Pos 1",
                "oh5": "GGAG", "oh3": "TGAC", "grammar": "gb_l0",
                "sequence": "ACGT" * 10,
            },
            {
                "name": "myMocloProm", "type": "Promoter",
                "position": "Pos 1", "oh5": "GGAG", "oh3": "AATG",
                "grammar": "moclo_plant", "sequence": "ACGT" * 10,
            },
        ])
        gb_rows = sc._palette_rows_for_grammar(
            "gb_l0", filter_enabled=True,
        )
        names = {r[0] for r in gb_rows}
        assert "myProm" in names
        assert "myMocloProm" not in names
        moclo_rows = sc._palette_rows_for_grammar(
            "moclo_plant", filter_enabled=True,
        )
        names = {r[0] for r in moclo_rows}
        assert "myMocloProm" in names
        assert "myProm" not in names

    def test_palette_helper_unfiltered_shows_all_user_parts(
            self, isolated_parts_bin):
        sc._save_parts_bin([
            {
                "name": "myProm", "type": "Promoter", "position": "Pos 1",
                "oh5": "GGAG", "oh3": "TGAC", "grammar": "gb_l0",
                "sequence": "ACGT" * 10,
            },
            {
                "name": "myMocloProm", "type": "Promoter",
                "position": "Pos 1", "oh5": "GGAG", "oh3": "AATG",
                "grammar": "moclo_plant", "sequence": "ACGT" * 10,
            },
        ])
        rows = sc._palette_rows_for_grammar(
            "gb_l0", filter_enabled=False,
        )
        names = {r[0] for r in rows}
        assert "myProm" in names
        assert "myMocloProm" in names

    def test_grammar_tu_overhangs_derives_from_positions(self):
        # gb_l0: first position oh5=GGAG (Promoter), last oh3=CGCT (Term).
        gb_l0 = sc._BUILTIN_GRAMMARS["gb_l0"]
        assert sc._grammar_tu_overhangs(gb_l0) == ("GGAG", "CGCT")
        # MoClo Plant: first oh5=GGAG, last oh3=CGCT (same TU
        # boundaries by Engler 2014 convention, different junctions).
        moclo = sc._BUILTIN_GRAMMARS["moclo_plant"]
        assert sc._grammar_tu_overhangs(moclo) == ("GGAG", "CGCT")

    def test_grammar_pos_slots_includes_cds_ns_alias(self):
        gb_l0 = sc._BUILTIN_GRAMMARS["gb_l0"]
        slots = sc._grammar_pos_slots(gb_l0)
        # CDS-NS shares the CDS slot via alias-fallback so the
        # duplicate-slot detection in `ConstructorModal._validate`
        # treats a no-stop CDS as occupying the same logical position
        # as a full CDS (you can't ligate both into the same TU lane).
        # Slot indices shift when the grammar is expanded — assert
        # the alias relationship, not a hardcoded slot number.
        assert "CDS" in slots
        assert slots["CDS-NS"] == slots["CDS"]
        # Promoter-only shares the Promoter slot for the same reason.
        assert "Promoter" in slots
        assert slots["Promoter-only"] == slots["Promoter"]

    async def test_constructor_modal_has_three_tabs(
            self, tiny_record, isolated_library):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            tab_ids = [p.id for p in modal.query("TabPane")]
            assert "ctor-tab-traditional" in tab_ids
            assert "ctor-tab-gb_l0"       in tab_ids
            assert "ctor-tab-moclo_plant" in tab_ids
            # Each grammar tab has its own palette + lane.
            from textual.widgets import DataTable
            assert modal.query_one("#ctor-palette-gb_l0", DataTable)
            assert modal.query_one("#ctor-lane-gb_l0",    DataTable)
            assert modal.query_one("#ctor-palette-moclo_plant", DataTable)
            assert modal.query_one("#ctor-lane-moclo_plant",    DataTable)

    async def test_modular_tab_palette_reflects_parts_bin(
            self, tiny_record, isolated_library, isolated_parts_bin,
    ):
        # User saves a Golden Braid CDS + a MoClo CDS. The GB
        # tab's palette should show only the GB part by default
        # (filter ON); the MoClo tab shows only the MoClo part.
        sc._save_parts_bin([
            {
                "name": "MyGBgene", "type": "CDS", "position": "Pos 3-4",
                "oh5": "AATG", "oh3": "GCTT", "grammar": "gb_l0",
                "sequence": "ATG" + "AAA" * 30 + "TAA",
            },
            {
                "name": "MyMocloGene", "type": "CDS",
                "position": "Pos 3", "oh5": "AGGT", "oh3": "GCTT",
                "grammar": "moclo_plant",
                "sequence": "ATG" + "GGG" * 30 + "TAA",
            },
        ])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            gb_names = {r[0] for r in modal._palette_rows.get("gb_l0", [])}
            mc_names = {r[0] for r in modal._palette_rows.get("moclo_plant", [])}
            assert "MyGBgene"     in gb_names
            assert "MyMocloGene"  not in gb_names
            assert "MyMocloGene"  in mc_names
            assert "MyGBgene"     not in mc_names

    async def test_modular_tab_per_grammar_lane_state(
            self, tiny_record, isolated_library, isolated_parts_bin,
    ):
        # Adding a part on the GB tab doesn't pollute the MoClo
        # tab's lane — each tab owns its own state.
        sc._save_parts_bin([
            {
                "name": "GBProm", "type": "Promoter", "position": "Pos 1",
                "oh5": "GGAG", "oh3": "TGAC", "grammar": "gb_l0",
                "sequence": "ACGT" * 10,
            },
        ])
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            modal._add_selected_part("gb_l0")
            assert len(modal._lanes["gb_l0"])       == 1
            assert len(modal._lanes["moclo_plant"]) == 0

    async def test_filter_setting_default_on(
            self, tiny_record, isolated_library, isolated_parts_bin,
    ):
        # The persisted setting defaults to True; the filter checkbox
        # in each modular tab inherits that default value.
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            from textual.widgets import Checkbox
            chk_gb = modal.query_one(
                "#chk-ctor-filter-gb_l0", Checkbox,
            )
            chk_mc = modal.query_one(
                "#chk-ctor-filter-moclo_plant", Checkbox,
            )
            assert chk_gb.value is True
            assert chk_mc.value is True

    async def test_filter_setting_persists_to_disk(
            self, tiny_record, isolated_library, isolated_parts_bin,
    ):
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            assert sc._get_setting(
                "constructor_filter_by_grammar", True,
            ) is True
            app.action_toggle_constructor_filter()
            await pilot.pause()
            sc._settings_flush_sync()
            assert sc._get_setting(
                "constructor_filter_by_grammar", True,
            ) is False

    async def test_per_role_entry_vector_banner(
            self, tiny_record, isolated_library, isolated_parts_bin,
            tmp_path, monkeypatch,
    ):
        # Configure different L1 acceptors for the four GB roles
        # (Alpha1 / Alpha2 / Omega1 / Omega2). After the user picks
        # a backbone via `_select_backbone`, the banner-style
        # summary reflects the active role's bound vector.
        ev_file = tmp_path / "entry_vectors.json"
        monkeypatch.setattr(sc._state, "_ENTRY_VECTORS_FILE", ev_file)
        monkeypatch.setattr(sc._state, "_entry_vectors_cache", None)
        for role, name in (
            ("Alpha1", "FFE2_test"),
            ("Alpha2", "FFE3_test"),
            ("Omega1", "FFE4_test"),
            ("Omega2", "FFE5_test"),
        ):
            sc._set_entry_vector("gb_l0", {
                "name":   name, "size": 60, "source": "test",
                "gb_text": (
                    f"LOCUS       {name}_locus           60 bp    DNA     circular SYN 01-JAN-2026\n"
                    "FEATURES             Location/Qualifiers\n"
                    "ORIGIN\n        1 " + "a" * 60 + "\n//\n"
                ),
            }, role=role)
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            # No backbone active by default — the user has to press
            # one to set it (matches the 2026-05-08 click-to-set UX).
            modal._select_backbone("gb_l0", "Alpha1")
            banner = modal._entry_vector_summary_for_grammar("gb_l0")
            assert "FFE2_test" in banner
            assert "Alpha1"    in banner
            # Switch to Omega1 → banner swaps.
            modal._select_backbone("gb_l0", "Omega1")
            banner = modal._entry_vector_summary_for_grammar("gb_l0")
            assert "FFE4_test" in banner
            assert "Omega1"    in banner

    async def test_clicking_backbone_button_routes_to_select_backbone(
            self, tiny_record, isolated_library, isolated_parts_bin,
    ):
        """Regression guard for 2026-05-08: the dispatcher's
        `_gid_from_button` only matched buttons whose id ENDED with
        `-{gid}`. Backbone buttons have an extra `-{role}` suffix
        (`btn-bb-{gid}-{role}`), so clicks were silently dropped —
        the buttons looked clickable but did nothing. Clicking each
        bb button must now invoke `_select_backbone(gid, role)`."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            calls: list[tuple] = []

            def _spy(gid, role):
                calls.append((gid, role))
                # Don't call real_select — it would try to push the
                # picker modal for unbound roles, which the test
                # doesn't need.

            modal._select_backbone = _spy
            from textual.widgets import TabbedContent
            modal.query_one(
                "#ctor-tabs", TabbedContent,
            ).active = "ctor-tab-gb_l0"
            await pilot.pause()
            await pilot.pause(0.1)
            await pilot.click("#btn-bb-gb_l0-Alpha1")
            await pilot.pause()
            await pilot.click("#btn-bb-gb_l0-Omega2")
            await pilot.pause()
            assert ("gb_l0", "Alpha1") in calls
            assert ("gb_l0", "Omega2") in calls

    async def test_role_buttons_inside_modal_horizontal_extent(
            self, tiny_record, isolated_library, isolated_parts_bin,
    ):
        """Regression guard for 2026-05-08: a leading
        ``Static("Backbone:")`` inside the row's Horizontal grabbed
        the full width and pushed the role columns past the modal's
        right edge (modal x=20..140 but Alpha2 at x=152). The fix
        drops the inline label in favour of a section header above
        the row. This test pins every role button's right edge
        inside the modal's container so it can't regress."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.1)
            from textual.widgets import TabbedContent
            modal.query_one(
                "#ctor-tabs", TabbedContent,
            ).active = "ctor-tab-gb_l0"
            await pilot.pause()
            await pilot.pause(0.1)
            box_region = modal.query_one("#ctor-box").region
            modal_right = box_region.x + box_region.width
            from textual.widgets import Button
            for role in ("Alpha1", "Alpha2", "Omega1", "Omega2"):
                btn = modal.query_one(f"#btn-bb-gb_l0-{role}", Button)
                btn_right = btn.region.x + btn.region.width
                assert btn_right <= modal_right, (
                    f"{role} button extends past modal right edge "
                    f"({btn.region} vs box {box_region})"
                )

    async def test_role_button_label_shows_bound_vector_name(
            self, tiny_record, isolated_library, isolated_parts_bin,
            tmp_path, monkeypatch,
    ):
        # Per-role static label ABOVE each backbone button shows the
        # bound vector's name (truncated for the column width). With
        # no binding the label reads "(none)".
        ev_file = tmp_path / "entry_vectors.json"
        monkeypatch.setattr(sc._state, "_ENTRY_VECTORS_FILE", ev_file)
        monkeypatch.setattr(sc._state, "_entry_vectors_cache", None)
        sc._set_entry_vector("gb_l0", {
            "name":   "MyAcceptor", "size": 3000, "source": "test",
            "gb_text": "",
        }, role="Alpha1")
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            # Alpha1 label shows the bound vector name.
            assert "MyAcceptor" in modal._backbone_name_label("gb_l0", "Alpha1")
            # Omega1 has nothing bound — label says "(none)".
            assert "(none)" in modal._backbone_name_label("gb_l0", "Omega1")
            from textual.widgets import Static
            lbl_a1 = modal.query_one(
                "#lbl-bb-gb_l0-Alpha1", Static,
            )
            lbl_o1 = modal.query_one(
                "#lbl-bb-gb_l0-Omega1", Static,
            )
            assert "MyAcceptor" in str(lbl_a1.render())
            assert "(none)"     in str(lbl_o1.render())

    async def test_unbound_button_press_opens_picker(
            self, tiny_record, isolated_library, isolated_parts_bin,
            tmp_path, monkeypatch,
    ):
        # When the user clicks a role button that has no plasmid
        # bound yet, `_select_backbone` should hand off to
        # `_pick_acceptor_for_role(then_activate=True)` rather than
        # silently activating an empty role. We verify by stubbing
        # `_pick_acceptor_for_role` and checking it gets called with
        # the right role + activate flag.
        ev_file = tmp_path / "entry_vectors.json"
        monkeypatch.setattr(sc._state, "_ENTRY_VECTORS_FILE", ev_file)
        monkeypatch.setattr(sc._state, "_entry_vectors_cache", None)
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            await pilot.pause(0.05)
            calls: list[tuple] = []

            def _stub(gid, role="", *, then_activate=False):
                calls.append((gid, role, then_activate))

            modal._pick_acceptor_for_role = _stub
            modal._select_backbone("gb_l0", "Alpha2")
            assert calls == [("gb_l0", "Alpha2", True)], (
                f"Unbound click should open picker; got {calls}"
            )

    def test_no_hardcoded_acceptor_ids_in_constructor_backbones(self):
        """Regression guard for 2026-05-08 (do-not-hardcode):
        `_CONSTRUCTOR_BACKBONES` must not contain a fixed ``id``
        field for any role. Roles are slots — the user binds a
        plasmid to each via "Change…" — never pre-bound to a
        specific vector by the source code."""
        for gid, roles in sc._CONSTRUCTOR_BACKBONES.items():
            for role, meta in roles.items():
                assert "id" not in meta, (
                    f"{gid}.{role} hardcodes id={meta.get('id')!r} "
                    f"— roles must not be pre-bound."
                )

    def test_set_get_entry_vector_role_isolation(
            self, tmp_path, monkeypatch,
    ):
        # Setting a vector for Alpha1 doesn't disturb Alpha2 / Omega1
        # / Omega2 / the legacy singleton (role="").
        ev_file = tmp_path / "entry_vectors.json"
        monkeypatch.setattr(sc._state, "_ENTRY_VECTORS_FILE", ev_file)
        monkeypatch.setattr(sc._state, "_entry_vectors_cache", None)
        sc._set_entry_vector("gb_l0", {
            "name": "L0_singleton", "size": 1, "source": "t",
            "gb_text": "",
        })
        sc._set_entry_vector("gb_l0", {
            "name": "alpha1_vec", "size": 2, "source": "t",
            "gb_text": "",
        }, role="Alpha1")
        sc._set_entry_vector("gb_l0", {
            "name": "omega1_vec", "size": 3, "source": "t",
            "gb_text": "",
        }, role="Omega1")
        # All three coexist.
        assert sc._get_entry_vector("gb_l0")["name"] == "L0_singleton"
        assert sc._get_entry_vector("gb_l0", "Alpha1")["name"] == "alpha1_vec"
        assert sc._get_entry_vector("gb_l0", "Omega1")["name"] == "omega1_vec"
        assert sc._get_entry_vector("gb_l0", "Alpha2") is None
        # Clearing one role doesn't disturb the others.
        sc._set_entry_vector("gb_l0", None, role="Alpha1")
        assert sc._get_entry_vector("gb_l0", "Alpha1") is None
        assert sc._get_entry_vector("gb_l0")["name"] == "L0_singleton"
        assert sc._get_entry_vector("gb_l0", "Omega1")["name"] == "omega1_vec"


class TestClonePartWrapperFeature:
    """INV-121: a part carrying lifted sub-features (a native operon's
    individual genes) clones WITHOUT the redundant whole-part wrapper that
    spans the entire insert — only the genes show. A part with NO
    sub-features keeps its single span so the cloned insert isn't an
    unannotated black box (user-reported: a cloned operon carried both the
    genes AND a feature named after the plasmid spanning the whole operon)."""

    def test_wrapper_dropped_when_genes_attach(self):
        # 5' "TTT" leader so the body doesn't start with ATG (avoids the
        # AATG-overhang start-codon collapse) — the insert anchors verbatim.
        body = ("TTT" + "ATG" + "GCT" * 60 + "TAA"   # geneA
                + "AAACCCGGGTTT"                      # spacer
                + "ATG" + "GCA" * 60 + "TGA")         # geneB
        g1_start, g1_end = 3, 3 + 3 + 180 + 3
        g2_start = g1_end + 12
        part = {
            "name": "myOperon", "type": "OPERON", "position": "Pos 3-4",
            "oh5": "AATG", "oh3": "GCTT", "sequence": body,
            "grammar": "gb_l0", "level": 0,
            "insert_feats": [
                {"start": g1_start, "end": g1_end, "strand": 1,
                 "type": "CDS", "label": "geneA"},
                {"start": g2_start, "end": len(body), "strand": 1,
                 "type": "CDS", "label": "geneB"},
            ],
        }
        rec = sc._part_to_cloned_seqrecord(part)
        labels = [str(f.qualifiers.get("label", [""])[0]) for f in rec.features]
        notes = [str(n) for f in rec.features
                 for n in (f.qualifiers.get("note") or [])]
        assert "geneA" in labels and "geneB" in labels       # genes carried
        assert not any(n.startswith("GB part type:") for n in notes), \
            f"whole-operon wrapper not dropped: {notes!r}"

    def test_wrapper_kept_without_subfeatures(self):
        part = {
            "name": "plainCDS", "type": "CDS", "position": "Pos 2-3",
            "oh5": "AATG", "oh3": "GCTT",
            "sequence": "TTTGCT" + "GCA" * 40 + "TAA",
            "grammar": "gb_l0", "level": 0,
        }
        rec = sc._part_to_cloned_seqrecord(part)
        notes = [str(n) for f in rec.features
                 for n in (f.qualifiers.get("note") or [])]
        assert any(n.startswith("GB part type:") for n in notes), \
            "featureless part lost its single-span wrapper"


# ══════════════════════════════════════════════════════════════════════════════
# Fragment-pick ambiguity + "that's the donor's backbone" guard
#
# Field report 2026-08-05: cloning a payload out of a donor and into a
# designated destination vector produced a construct carrying the DONOR's
# resistance cassette instead of the expression circuit.
#
# Mechanism (all three parts had to line up):
#   1. Cut donor and vector with the SAME enzyme pair and the donor's two
#      halves come off the circle with MIRRORED overhangs — the payload
#      reads BamHI→…→XbaI while the backbone half reads XbaI→…→BamHI. Only
#      one of them matches the vector in the FORWARD orientation, and which
#      one is a coin flip. Here the payload ligated only in REVERSE.
#   2. The results pane always led with "Forward orientation / ✗ no
#      ligation", which reads as a REJECTED FRAGMENT at a glance.
#   3. Switching to the half that DID ligate forward produced a perfectly
#      valid — and biologically wrong — plasmid, with no warning at all.
#
# Plus a latent hole found while tracing it: both pickers filter on the
# backbone-marker predicate and fall back to size when the filter doesn't
# leave exactly one candidate, but only warned when NO fragment carried a
# marker. An ALL-marked digest (very common: a vector's tiny inter-MCS
# stuffer still carries a slice of the annotated origin) emptied the
# candidate list just as thoroughly and fell through in total silence.
# ══════════════════════════════════════════════════════════════════════════════


def _frag(nbp: int, labels: list[str], ftype: str = "misc_feature") -> dict:
    """A minimal fragment dict — only the fields the pickers read."""
    return {"top_seq": "A" * nbp,
            "features": [{"start": 0, "end": 3, "type": ftype,
                          "label": lb} for lb in labels]}


class TestFragmentPickAmbiguity:
    def test_all_fragments_marked_is_flagged_ambiguous(self):
        """The silent hole. Both halves carry a marker => the marker filter
        leaves ZERO candidates => the pick is made by size and MUST say so."""
        small, big = _frag(500, ["AmpR"]), _frag(4000, ["p15A ori"])
        picked = sc._pick_insert_fragment([small, big])
        assert picked is small                      # size fallback: smallest
        assert picked.get("_size_fallback_no_marker") is True, (
            "an all-marked digest fell back to size with no ambiguity flag — "
            "the Constructor's warning never fires and the wrong half clones "
            "silently"
        )

    def test_no_fragments_marked_still_flagged(self):
        """The case that already worked — must keep working."""
        small, big = _frag(500, ["geneA"]), _frag(4000, ["geneB"])
        picked = sc._pick_insert_fragment([small, big])
        assert picked is small
        assert picked.get("_size_fallback_no_marker") is True

    def test_exactly_one_marked_is_unambiguous(self):
        """Marker evidence resolved it — no fallback, so no warning noise."""
        payload, backbone = _frag(4000, ["geneA"]), _frag(500, ["AmpR"])
        picked = sc._pick_insert_fragment([payload, backbone])
        assert picked is payload, "marker evidence must beat the size heuristic"
        assert "_size_fallback_no_marker" not in picked

    def test_backbone_pick_flags_all_marked(self):
        """`_pick_backbone_fragment`'s mirror of the same hole: several
        candidates is just as inconclusive as none."""
        stuffer = _frag(40, ["p15A ori"])
        backbone = _frag(6000, ["ermB"], ftype="gene")
        picked = sc._pick_backbone_fragment([stuffer, backbone])
        assert picked is backbone                   # size fallback: largest
        assert picked.get("_size_fallback_no_marker") is True

    def test_backbone_pick_unambiguous_when_one_marked(self):
        payload, backbone = _frag(4000, ["geneA"]), _frag(500, ["KanR"])
        picked = sc._pick_backbone_fragment([payload, backbone])
        assert picked is backbone
        assert "_size_fallback_no_marker" not in picked

    def test_single_fragment_is_never_ambiguous(self):
        """Nothing to choose between — don't cry wolf."""
        only = _frag(500, ["AmpR"])
        assert sc._pick_insert_fragment([only]) is only
        assert "_size_fallback_no_marker" not in only


class TestBackboneMarkerLabels:
    def test_names_the_markers_it_matched(self):
        frag = _frag(500, ["Ori*", "AmpR", "AmpR promoter", "geneA"])
        labels = sc._fragment_backbone_marker_labels(frag)
        assert labels == ["Ori*", "AmpR", "AmpR promoter"]
        assert "geneA" not in labels

    def test_empty_when_predicate_is_false(self):
        frag = _frag(500, ["geneA", "terminator"])
        assert sc._fragment_backbone_marker_labels(frag) == []
        assert sc._fragment_has_backbone_marker(frag) is False

    def test_agrees_with_the_predicate(self):
        """The list is the predicate's evidence — they must never disagree."""
        for labels in (["AmpR"], ["geneA"], ["Ori*", "x"], [], ["rep_origin"]):
            frag = _frag(100, labels)
            assert bool(sc._fragment_backbone_marker_labels(frag)) is \
                sc._fragment_has_backbone_marker(frag), labels

    def test_matches_on_feature_type_with_no_label(self):
        frag = {"top_seq": "A" * 100,
                "features": [{"start": 0, "end": 3, "type": "rep_origin",
                              "label": ""}]}
        assert sc._fragment_has_backbone_marker(frag) is True
        assert sc._fragment_backbone_marker_labels(frag) == ["rep_origin"]


class TestMirroredOverhangsAreDirectional:
    """The chemistry behind the field report, on synthetic sequence: the two
    halves of a doubly-cut circle carry their overhangs in OPPOSITE order,
    so exactly one of them ligates 'forward' into a like-cut vector. Both
    are valid products — only which DNA rides along differs."""

    def test_both_halves_ligate_but_in_opposite_orientations(self):
        # 40 bp circle, EcoRI at 4 and BamHI at 18 => two halves.
        donor = "AAAAGAATTCAAAAAAAAGGATCCAAAAAAAAAAAAAAAA"
        frags, err = sc._excise_fragment_pair(
            donor, ["EcoRI", "BamHI"], circular=True,
            features=[{"start": 11, "end": 17, "type": "CDS",
                       "label": "payload"},
                      {"start": 30, "end": 38, "type": "CDS",
                       "label": "AmpR"}],
            source_label="donor")
        assert err is None and len(frags) == 2
        # Mirrored ends is the whole point — assert it explicitly.
        a, b = frags
        assert (a["left"]["overhang_seq"], a["right"]["overhang_seq"]) == \
               (b["right"]["overhang_seq"], b["left"]["overhang_seq"]), \
            "halves of a doubly-cut circle should carry mirrored overhangs"
        vector = sc._pick_backbone_fragment(list(frags))
        # Each half ligates into that vector in exactly ONE orientation, and
        # they disagree about which — the trap the user fell into.
        seen = []
        for f in frags:
            out = sc._simulate_traditional_cloning(f, vector)
            fwd = out["forward"]["compatible"]
            rev = out["reverse"]["compatible"]
            assert fwd != rev, "expected a directional (one-way) result"
            assert not out["errors"], out["errors"]
            seen.append(fwd)
        assert seen[0] != seen[1], (
            "the two halves should prefer OPPOSITE orientations — that is "
            "why 'Forward ✗' on the right fragment must not read as 'this "
            "fragment cannot be cloned'"
        )


class TestBackboneMarkerVocabulary:
    """Short resistance-gene names must match as WHOLE TOKENS, never as
    substrings — and the vocabulary has to actually cover the markers real
    shuttle vectors use.

    Found 2026-08-05: `ermB` (the standard Gram-positive MLS marker, and the
    selection on common Gram-positive shuttle vectors) was absent from both
    this predicate's keywords and `_SELECTION_MARKER_KEYWORDS`. A backbone
    carrying it read as marker-free, so insert/backbone selection fell back
    to fragment size with nothing logged."""

    @pytest.mark.parametrize("label", ["ermB", "ermC", "erm", "ErmB",
                                       "ermB resistance", "MLS (ermAM)"])
    def test_erm_family_is_a_backbone_marker(self, label):
        frag = _frag(500, [label], ftype="gene")
        assert sc._fragment_has_backbone_marker(frag) is True, \
            f"{label!r} not recognised as a selection marker"

    @pytest.mark.parametrize("label", [
        "Terminator", "Term 908", "Lambda T0 Terminator",
        "T7 terminator", "rrnB T1 terminator",
    ])
    def test_terminators_are_not_backbone_markers(self, label):
        """The reason "erm" can't be a substring: every construct has a
        terminator, and marking them all would disable the heuristic."""
        frag = _frag(500, [label], ftype="gene")
        assert sc._fragment_has_backbone_marker(frag) is False, \
            f"{label!r} false-positived as a backbone marker"

    @pytest.mark.parametrize("label", ["category", "blast hit", "catalase"])
    def test_short_name_substrings_do_not_false_positive(self, label):
        frag = _frag(500, [label], ftype="gene")
        assert sc._fragment_has_backbone_marker(frag) is False

    @pytest.mark.parametrize("label", ["cat", "bla", "aadA", "nptII"])
    def test_short_names_match_as_whole_tokens(self, label):
        frag = _frag(500, [label], ftype="gene")
        assert sc._fragment_has_backbone_marker(frag) is True

    @pytest.mark.parametrize("ftype,label", [
        # Real annotations from the field report. Both mention a marker but
        # neither IS one — a promoter annotated after another plasmid, and
        # oligos that anneal near the gene. Marking these tagged the PAYLOAD
        # half as backbone, which made both halves look like backbone and
        # collapsed the pick back to raw fragment size.
        ("regulatory",  "as_annotated_pXX-ermB"),
        ("primer_bind", "ermB-CLO-F"),
        ("primer_bind", "PermB-CLO-R"),
        ("primer_bind", "PROM ermB"),
        ("misc_feature", "ermB homology arm"),
    ])
    def test_gene_name_outside_a_gene_feature_is_not_a_marker(
            self, ftype, label):
        frag = _frag(500, [label], ftype=ftype)
        assert sc._fragment_has_backbone_marker(frag) is False, (
            f"{ftype} /label={label!r} read as a backbone marker — it names "
            f"something NEAR the gene, not the gene"
        )

    @pytest.mark.parametrize("ftype", ["gene", "CDS"])
    def test_real_marker_gene_is_matched(self, ftype):
        frag = _frag(500, ["ermB"], ftype=ftype)
        assert sc._fragment_has_backbone_marker(frag) is True

    def test_long_form_keywords_ignore_feature_type(self):
        """The pre-existing substring keywords are specific enough that they
        stay unrestricted — don't narrow them and lose coverage."""
        for ftype in ("misc_feature", "regulatory", "CDS"):
            assert sc._fragment_has_backbone_marker(
                _frag(500, ["AmpR"], ftype=ftype)) is True
            assert sc._fragment_has_backbone_marker(
                _frag(500, ["p15A ori"], ftype=ftype)) is True

    def test_cloning_mirror_agrees(self):
        """`splicecraft_cloning._frag_carries_backbone_marker` is a hand-kept
        copy (L3 can't import seqanalysis — import cycle). Drift between the
        two means the Constructor and the modular resolver disagree about
        which half is the backbone. (The full corpus runs through both in
        `TestBackboneVocabularyCorpus`; these are the erm-family cases.)"""
        import splicecraft_cloning as scc
        cases = [("gene", ["ermB"]), ("gene", ["Terminator"]),
                 ("misc_feature", ["AmpR"]), ("misc_feature", ["p15A ori"]),
                 ("gene", ["cat"]), ("gene", ["category"]),
                 ("CDS", ["geneA"]), ("rep_origin", [""]),
                 ("regulatory", ["as_annotated_pXX-ermB"]),
                 ("primer_bind", ["ermB-CLO-F"]),
                 ("gene", ["Term 908"]), ("CDS", ["nptII"]),
                 ("misc_feature", ["blast"]),
                 ("misc_feature", ["EcoRI"]), ("primer_bind", ["EcoRI-F"]),
                 ("CDS", ["GFP (Aequorea victoria)"]),
                 ("misc_feature", ["oriT"]), ("misc_feature", ["pMB1ori"]),
                 ("misc_feature", ["origin of replication"])]
        for ftype, labels in cases:
            frag = _frag(100, labels, ftype=ftype)
            assert scc._frag_carries_backbone_marker(frag) is \
                sc._fragment_has_backbone_marker(frag), \
                f"marker mirrors disagree on {ftype} {labels!r}"

    def test_detect_selection_marker_knows_erm(self):
        """The user-visible marker name must resolve too — not just the
        internal backbone heuristic."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        rec = SeqRecord(Seq("ACGT" * 100), id="pShuttle", name="pShuttle")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        rec.features.append(SeqFeature(
            FeatureLocation(10, 100), type="CDS",
            qualifiers={"label": ["ermB"]}))
        assert sc._detect_selection_marker(
            sc._record_to_gb_text(rec)) == "Erythromycin"


class TestOriNeedsALeftBoundary:
    """`ori` is three letters, and it is the TAIL of names that get annotated
    on payload fragments constantly — `EcoRI` above all.

    Field report 2026-08-17 (agent-driven Golden Gate): an alpha/omega digest
    whose INSERT half carried nothing but an `EcoRI` site annotation read as
    backbone, so BOTH halves did, and the pick fell through to the size
    heuristic the marker test exists to replace. The reporter worked around
    it by selecting on the `rep_origin` feature TYPE.
    """

    # The exact digest from the report: the 886 nt insert half is annotated
    # only with an EcoRI site; the 1822 nt backbone half carries the real
    # markers (and an EcoRI site of its own).
    INSERT_LABELS = ["EcoRI"]
    BACKBONE_LABELS = ["AmpR", "AmpR promoter", "Ori*", "EcoRI"]

    def test_the_field_report_digest_picks_the_right_half(self):
        insert = _frag(886, self.INSERT_LABELS)
        backbone = _frag(1822, self.BACKBONE_LABELS)
        assert sc._fragment_has_backbone_marker(insert) is False, (
            "an EcoRI site annotation made the insert half read as backbone"
        )
        assert sc._fragment_has_backbone_marker(backbone) is True
        picked = sc._pick_insert_fragment([insert, backbone])
        assert picked is insert
        assert "_size_fallback_no_marker" not in picked, (
            "marker evidence singled out one fragment — the pick must NOT be "
            "a size guess"
        )

    def test_the_field_report_digest_picks_the_right_backbone(self):
        insert = _frag(886, self.INSERT_LABELS)
        backbone = _frag(1822, self.BACKBONE_LABELS)
        picked = sc._pick_backbone_fragment([insert, backbone])
        assert picked is backbone
        assert "_size_fallback_no_marker" not in picked

    def test_report_does_not_name_a_restriction_site_as_a_marker(self):
        """The reporter's `markers=['EcoRI']` diagnosis line came from here."""
        assert sc._fragment_backbone_marker_labels(
            _frag(886, self.INSERT_LABELS)) == []
        assert sc._fragment_backbone_marker_labels(
            _frag(1822, self.BACKBONE_LABELS)) == ["AmpR", "AmpR promoter",
                                                   "Ori*"]

    @pytest.mark.parametrize("ftype,label", [
        ("misc_feature", "EcoRI"),
        ("misc_feature", "EcoRI site"),
        ("misc_feature", "EcoRI/BamHI MCS"),
        ("primer_bind", "EcoRI-F"),
        ("CDS", "GFP (Aequorea victoria)"),
        ("CDS", "mCherry (Discosoma sp.)"),
        ("misc_feature", "a priori design"),
        ("misc_feature", "historic clone"),
        # The two English words that OPEN with it — a left boundary alone
        # lets both through.
        ("misc_feature", "insert (original)"),
        ("misc_feature", "original construct"),
        ("misc_feature", "GFP, reverse orientation"),
        ("misc_feature", "orientation marker"),
        ("CDS", "originally from pBR322"),
    ])
    def test_letter_joined_ori_is_not_an_origin(self, ftype, label):
        frag = _frag(500, [label], ftype=ftype)
        assert sc._fragment_has_backbone_marker(frag) is False, (
            f"{label!r} false-positived as an origin of replication"
        )
        assert sc._fragment_backbone_marker_labels(frag) == []

    @pytest.mark.parametrize("label", [
        "ori", "Ori*", "ORI", "ori (pMB1)", "ColE1 ori", "pUC ori", "f1 ori",
        "p15A ori", "pSC101 ori", "R6K gamma ori", "oriT", "oriV", "oriC",
        "ori2", "oriP", "oriR6K", "oriColE1", "rep_origin",
        "origin of replication", "origins of replication", "pBR322 origin",
        "pMB1ori",
    ])
    def test_every_real_origin_spelling_still_matches(self, label):
        """The boundary must not cost coverage — a false NEGATIVE here drops
        the pick back to the same size heuristic, just quietly."""
        frag = _frag(500, [label])
        assert sc._fragment_has_backbone_marker(frag) is True, (
            f"{label!r} stopped being recognised as an origin"
        )

    def test_a_marker_past_the_display_clamp_still_reports(self):
        """`_fragment_backbone_marker_labels` clamps names to 80 chars for the
        one-line warning. Clamping BEFORE matching would make the predicate
        say True while the evidence list came back empty."""
        label = "x" * 100 + " ColE1 ori"
        frag = _frag(500, [label])
        assert sc._fragment_has_backbone_marker(frag) is True
        names = sc._fragment_backbone_marker_labels(frag)
        assert names, "predicate said backbone, evidence list said nothing"
        assert len(names[0]) <= 80

    def test_cloning_mirror_agrees_on_ori(self):
        import splicecraft_cloning as scc
        for ftype, label in [("misc_feature", "EcoRI"),
                             ("primer_bind", "EcoRI-F"),
                             ("CDS", "GFP (Aequorea victoria)"),
                             ("misc_feature", "ColE1 ori"),
                             ("misc_feature", "oriT"),
                             ("misc_feature", "pMB1ori"),
                             ("misc_feature", "insert (original)"),
                             ("misc_feature", "reverse orientation"),
                             ("misc_feature", "origin of replication")]:
            frag = _frag(100, [label], ftype=ftype)
            assert scc._frag_carries_backbone_marker(frag) is \
                sc._fragment_has_backbone_marker(frag), \
                f"marker mirrors disagree on {ftype} /label={label!r}"


# ── The adversarial corpus behind TestBackboneVocabularyCorpus ───────────────
# Real annotation text, not invented strings: SnapGene / pLannotate / Benchling
# feature names as they actually appear on plasmid maps. Every entry here was
# run through the predicate on 2026-08-18 while chasing the `ori`-in-`EcoRI`
# field report, and the ones marked below FAILED before that sweep.
_PAYLOAD_ANNOTATIONS = [
    # restriction sites — `misc_feature` is how every tool writes them
    ("misc_feature", "EcoRI"),                 # ← the field report
    ("misc_feature", "EcoRI site"),
    ("misc_feature", "EcoRI/BamHI MCS"),
    ("primer_bind",  "EcoRI-F"),
    ("misc_feature", "BamHI"), ("misc_feature", "XbaI"),
    ("misc_feature", "EcoRV"), ("misc_feature", "Esp3I"),
    ("misc_feature", "Tth111I"), ("misc_feature", "PaqCI"),
    # species / domain names with unlucky substrings
    ("CDS", "GFP (Aequorea victoria)"),        # ← vict-ORI-a
    ("CDS", "mCherry (Discosoma sp.)"),
    ("CDS", "tetratricopeptide repeat"),       # ← TETR-atricopeptide
    ("CDS", "TPR (tetratricopeptide) domain"),
    ("CDS", "tetramerization domain"),         # ← TETR-amerization
    ("CDS", "tetraspanin"),
    ("CDS", "Tetrahymena ribozyme"),
    ("misc_feature", "tetraloop"),
    ("CDS", "campR-like protein"),             # ← c-AMPR-like
    # tet regulatory parts — payload in every Tet-On / Tet-Off construct
    ("protein_bind", "tetracycline operator"),
    ("protein_bind", "tetO"),
    ("promoter", "tetracycline-responsive promoter"),
    ("promoter", "TRE (tet response element)"),
    ("misc_feature", "tetR binding site"),
    ("primer_bind", "TetR-F"),
    # ordinary English on hand-annotated features
    ("misc_feature", "insert (original)"),
    ("misc_feature", "original construct"),
    ("misc_feature", "GFP, reverse orientation"),
    ("CDS", "originally from pBR322"),
    ("misc_feature", "a priori design"),
    ("misc_feature", "historic clone"),
    # parts / tags / sites that must stay payload
    ("CDS", "lacZ-alpha"), ("CDS", "LacI"), ("CDS", "AraC"),
    ("promoter", "T7 promoter"), ("terminator", "T7 terminator"),
    ("terminator", "rrnB T1 terminator"), ("misc_feature", "MCS"),
    ("CDS", "6xHis"), ("CDS", "FLAG tag"), ("misc_feature", "loxP"),
    ("primer_bind", "M13 fwd"), ("CDS", "rop"), ("misc_feature", "bom"),
    ("CDS", "catalase"), ("gene", "category B part"), ("gene", "blast hit"),
    ("CDS", "spectrin repeat"), ("CDS", "kanamycin-insensitive kinase"),
]

_MARKER_ANNOTATIONS = [
    # origins, every spelling that turns up in the wild
    ("rep_origin", "ori"), ("rep_origin", "Ori*"), ("misc_feature", "ColE1 ori"),
    ("rep_origin", "pUC ori"), ("rep_origin", "f1 ori"),
    ("rep_origin", "p15A ori"), ("rep_origin", "pVS1 oriV"),
    ("misc_feature", "oriT"), ("misc_feature", "oriP"),
    ("misc_feature", "oriR6K"), ("misc_feature", "ori2"),
    ("misc_feature", "origin of replication"),
    ("misc_feature", "origins of replication"),
    ("misc_feature", "pBR322 origin"), ("misc_feature", "pMB1ori"),
    ("misc_feature", "rep_origin"),
    # resistance
    ("CDS", "AmpR"), ("CDS", "KanR"), ("CDS", "NeoR/KanR"), ("CDS", "SpecR"),
    ("CDS", "CmR"), ("CDS", "ermB"), ("CDS", "cat"), ("CDS", "bla"),
    ("CDS", "aadA"), ("CDS", "nptII"), ("CDS", "hygR"),
    ("promoter", "AmpR promoter"), ("promoter", "KanR Promoter"),
    ("CDS", "chloramphenicol acetyltransferase"),
    ("CDS", "carbenicillin resistance"), ("CDS", "spectinomycin resistance"),
    # tet — on a gene-bearing feature type, which is what makes it the gene
    ("CDS", "TetR"), ("gene", "tetR"),
    ("CDS", "tetracycline resistance protein"), ("gene", "tetA"),
]


class TestBackboneVocabularyCorpus:
    """The whole vocabulary against real annotation text, in both mirrors.

    A false POSITIVE here is the expensive direction: it makes a payload half
    look like backbone, and once both halves look like backbone every consumer
    falls back to guessing by size — which is the heuristic the marker test
    exists to replace. A false NEGATIVE only drops to that same fallback, with
    a warning. So the corpus is scored asymmetrically on purpose.
    """

    @pytest.mark.parametrize("ftype,label", _PAYLOAD_ANNOTATIONS)
    def test_payload_annotation_is_not_a_marker(self, ftype, label):
        frag = _frag(500, [label], ftype=ftype)
        assert sc._fragment_has_backbone_marker(frag) is False, (
            f"{ftype} /label={label!r} read as backbone — a payload fragment "
            f"carrying it would be mistaken for the vector half"
        )
        assert sc._fragment_backbone_marker_labels(frag) == []

    @pytest.mark.parametrize("ftype,label", _MARKER_ANNOTATIONS)
    def test_real_marker_is_recognised(self, ftype, label):
        frag = _frag(500, [label], ftype=ftype)
        assert sc._fragment_has_backbone_marker(frag) is True, (
            f"{ftype} /label={label!r} stopped being recognised as a backbone "
            f"marker — the pick silently degrades to fragment size"
        )

    @pytest.mark.parametrize(
        "ftype,label", _PAYLOAD_ANNOTATIONS + _MARKER_ANNOTATIONS)
    def test_cloning_mirror_agrees(self, ftype, label):
        """`splicecraft_cloning._frag_carries_backbone_marker` is a hand-kept
        copy (L3 can't import seqanalysis — import cycle). Drift means the
        Constructor and the domestication resolver disagree about which half
        of a digest is the backbone."""
        import splicecraft_cloning as scc
        frag = _frag(100, [label], ftype=ftype)
        assert scc._frag_carries_backbone_marker(frag) is \
            sc._fragment_has_backbone_marker(frag), \
            f"marker mirrors disagree on {ftype} /label={label!r}"

    @pytest.mark.parametrize(
        "ftype,label", _PAYLOAD_ANNOTATIONS + _MARKER_ANNOTATIONS)
    def test_evidence_list_agrees_with_the_predicate(self, ftype, label):
        """The list IS the predicate's evidence — one shared matcher, so they
        cannot disagree by construction. This pins that they still share it."""
        frag = _frag(100, [label], ftype=ftype)
        assert bool(sc._fragment_backbone_marker_labels(frag)) is \
            sc._fragment_has_backbone_marker(frag)


class TestKeywordMatchingRule:
    """The word-start rule itself, independent of any particular keyword."""

    def test_keyword_matches_at_a_word_start_only(self):
        import splicecraft_seqanalysis as sa
        rx = sa._label_keyword_re(("zzmarker",), {})
        assert rx.search("zzmarker")            # bare
        assert rx.search("pUC zzmarker")        # after a space
        assert rx.search("rep-zzmarker")        # after a separator
        assert rx.search("p15zzmarker")         # after a DIGIT — pMB1ori
        assert rx.search("ZZMarker")            # case-insensitive
        assert not rx.search("xzzmarker")       # mid-word — EcoRI's disease
        assert not rx.search("Ecozzmarker")

    def test_traps_block_the_words_a_keyword_opens(self):
        import splicecraft_seqanalysis as sa
        rx = sa._label_keyword_re(("zz",), {"zz": ("top", "bar")})
        assert rx.search("zz") and rx.search("zzebra")
        assert not rx.search("zztop") and not rx.search("zzbarrel")

    def test_regex_is_built_from_the_keyword_tuple(self):
        """Adding a keyword to the tuple must take effect — the regex is not
        a hand-written second copy of the vocabulary."""
        import splicecraft_seqanalysis as sa
        for kw in sa._BACKBONE_LABEL_KEYWORDS:
            assert sa._BACKBONE_LABEL_RE.search(kw), \
                f"{kw!r} is in the vocabulary but its own regex misses it"

    def test_an_empty_vocabulary_matches_nothing(self):
        """Fail CLOSED. An empty alternation compiles to a pattern that
        matches the empty string at every position — i.e. every fragment
        would read as backbone, which is the catastrophic direction."""
        import splicecraft_seqanalysis as sa
        for kws in ((), ("",), ("", "")):
            rx = sa._label_keyword_re(kws, {})
            assert not rx.search("")
            assert not rx.search("AmpR")
            assert not rx.search("anything at all")

    def test_word_start_rule_implies_substring(self):
        """Soundness of the hot-path gate: `_feature_is_backbone_marker` only
        runs the regex when a bare keyword substring is present, which is only
        valid because the regex can never match where the substring doesn't.
        Fuzzed against strings built from keyword fragments + noise."""
        import random
        import splicecraft_seqanalysis as sa
        rnd = random.Random(20260818)
        alphabet = ["ori", "tet", "amp", "r", "kan", "-", " ", "1", "x",
                    "gina", "ent", "ecoRI", "_", "rep", "cm", "spec"]
        for _ in range(4000):
            s = "".join(rnd.choice(alphabet)
                        for _ in range(rnd.randint(1, 8)))
            if sa._BACKBONE_LABEL_RE.search(s):
                low = s.lower()
                assert any(kw in low for kw in sa._BACKBONE_LABEL_KEYWORDS), (
                    f"{s!r} matches the rule but no bare keyword occurs in it "
                    f"— the substring pre-filter would skip it"
                )


class TestTetIsARegulatoryPartOutsideAGene:
    """tet is the one marker family whose repressor and operator are standard
    PAYLOAD parts (tetO / TRE / tet-responsive promoter). Same rule the erm
    family already uses: the name counts on a gene-bearing feature type."""

    @pytest.mark.parametrize("ftype", ["CDS", "gene"])
    @pytest.mark.parametrize("label", [
        "TetR", "tetR", "tetracycline resistance protein", "tetA", "tetM"])
    def test_tet_on_a_gene_is_a_marker(self, ftype, label):
        assert sc._fragment_has_backbone_marker(
            _frag(500, [label], ftype=ftype)) is True

    @pytest.mark.parametrize("ftype", [
        "protein_bind", "promoter", "primer_bind", "misc_feature",
        "regulatory"])
    @pytest.mark.parametrize("label", [
        "tetracycline operator", "tetR binding site", "TetR-F"])
    def test_tet_outside_a_gene_is_not(self, ftype, label):
        assert sc._fragment_has_backbone_marker(
            _frag(500, [label], ftype=ftype)) is False


class TestMarkerHelpersSurviveMalformedFragments:
    """Both helpers run inside Constructor workers, on digest output AND on
    caller-supplied JSON from the agent API. Malformed input must read as
    "no marker", never raise out of the worker."""

    @pytest.mark.parametrize("frag", [
        None, "not a fragment", 42, [],
        {}, {"features": None}, {"features": "AmpR"}, {"features": 7},
        {"features": [None, 3, "AmpR"]},
        {"features": [{"type": None, "label": None}]},
        {"features": [{}]},
        {"features": ({"type": "rep_origin", "label": ""},)},   # tuple
    ])
    def test_no_raise_on_junk(self, frag):
        import splicecraft_cloning as scc
        got = sc._fragment_has_backbone_marker(frag)
        assert isinstance(got, bool)
        assert isinstance(sc._fragment_backbone_marker_labels(frag), list)
        assert scc._frag_carries_backbone_marker(frag) is got, \
            "the mirrors disagree on malformed input"

    def test_a_tuple_of_features_still_matches(self):
        """The tuple case above is the only junk input that IS a marker."""
        assert sc._fragment_has_backbone_marker(
            {"features": ({"type": "rep_origin", "label": ""},)}) is True

    def test_non_string_label_does_not_crash_or_match(self):
        frag = {"features": [{"type": "CDS", "label": ["AmpR"]}]}
        # A list label stringifies to "['AmpR']" — the keyword is there, at a
        # word start (after `['`), so it matches. What matters is no raise.
        assert isinstance(sc._fragment_has_backbone_marker(frag), bool)


class TestResultsPaneHardening:
    """The results pane interpolates user-controlled text — lane-row names and
    the feature LABELS naming a backbone marker — into a `Static`, which
    interprets Rich markup. An ordinary label like "TU [draft]" parses as an
    unclosed tag: it swallows the rest of the message, or raises MarkupError
    and blanks the pane outright."""

    @staticmethod
    def _outcome(warning: str) -> dict:
        return {
            "forward": {"top_seq": "A" * 100, "features": [],
                        "compatible": True},
            "reverse": {"top_seq": "A" * 100, "features": [],
                        "compatible": False},
            "warnings": [warning],
            "errors": [],
        }

    @pytest.mark.parametrize("hostile", [
        "insert TU [draft] carries markers",
        "AmpR [old] / Ori* [v2]",
        "unclosed [bold and [/] stray close",
        "[link=http://x]not a link[/link]",
    ])
    async def test_markup_in_a_warning_is_escaped(
            self, hostile, tiny_record, isolated_library):
        from textual.widgets import Static

        from tests.test_smoke import TERMINAL_SIZE, _build_app

        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane", sc.TraditionalCloningPane)
            target = pane.query_one("#trad-results-text", Static)
            pane._render_results(self._outcome(hostile), target)
            await pilot.pause()
            # Must survive markup parsing AND still say what it meant to say.
            # The contract is "markup is NEUTRALISED, message intact" — not
            # byte-identical: on a pathological input (a stray `[/]` after an
            # escaped opener) Textual's Content parser leaves the escaping
            # backslash visible. That's cosmetic. What must never happen is
            # the old behaviour: the tag swallowing the rest of the sentence,
            # or MarkupError blanking the pane.
            plain = target.visual.plain.replace("\\", "")
            assert hostile.replace("\\", "") in plain, (
                f"warning text mangled by markup parsing: {plain!r}"
            )

    async def test_error_text_is_escaped_too(
            self, tiny_record, isolated_library):
        from textual.widgets import Static

        from tests.test_smoke import TERMINAL_SIZE, _build_app

        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane", sc.TraditionalCloningPane)
            target = pane.query_one("#trad-results-text", Static)
            out = self._outcome("ok")
            out["errors"] = ["no overhangs matched on [pFoo v2]"]
            pane._render_results(out, target)
            await pilot.pause()
            plain = target.visual.plain
            assert "[pFoo v2]" in plain

    @pytest.mark.parametrize("fwd_ok,rev_ok,expect", [
        (True,  False, "FORWARD"),
        (False, True,  "REVERSE"),
    ])
    async def test_directional_result_names_the_button(
            self, fwd_ok, rev_ok, expect, tiny_record, isolated_library):
        """The misread that produced a wrong construct: a one-way ligation
        must not open with the orientation that DOESN'T work."""
        from textual.widgets import Static

        from tests.test_smoke import TERMINAL_SIZE, _build_app

        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane", sc.TraditionalCloningPane)
            target = pane.query_one("#trad-results-text", Static)
            out = self._outcome("x")
            out["forward"]["compatible"] = fwd_ok
            out["reverse"]["compatible"] = rev_ok
            pane._render_results(out, target)
            await pilot.pause()
            plain = target.visual.plain
            head = plain.strip().splitlines()[0]
            assert expect in head, (
                f"one-way ligation did not lead with the working orientation: "
                f"{head!r}"
            )

    async def test_no_headline_when_both_or_neither_ligate(
            self, tiny_record, isolated_library):
        """Don't invent a recommendation when there isn't one to make."""
        from textual.widgets import Static

        from tests.test_smoke import TERMINAL_SIZE, _build_app

        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            modal = sc.ConstructorModal()
            await app.push_screen(modal)
            await pilot.pause()
            pane = modal.query_one("#ctor-trad-pane", sc.TraditionalCloningPane)
            target = pane.query_one("#trad-results-text", Static)
            for fwd_ok, rev_ok in ((True, True), (False, False)):
                out = self._outcome("x")
                out["forward"]["compatible"] = fwd_ok
                out["reverse"]["compatible"] = rev_ok
                pane._render_results(out, target)
                await pilot.pause()
                plain = target.visual.plain
                assert "Ligates as" not in plain, (fwd_ok, rev_ok, plain[:80])


class TestMarkerLabelsHardening:
    @pytest.mark.parametrize("bad", [None, "", 0, [], "notadict"])
    def test_non_dict_fragment_returns_empty(self, bad):
        assert sc._fragment_backbone_marker_labels(bad) == []

    def test_missing_and_malformed_features_are_skipped(self):
        assert sc._fragment_backbone_marker_labels({"top_seq": "A"}) == []
        assert sc._fragment_backbone_marker_labels(
            {"features": None}) == []
        assert sc._fragment_backbone_marker_labels(
            {"features": ["AmpR", None, 7, {"type": "gene",
                                            "label": "AmpR"}]}) == ["AmpR"]

    def test_pathological_label_is_clamped(self):
        """An imported GenBank /label can be enormous; it must not swamp the
        one-line warning it feeds."""
        frag = _frag(100, ["AmpR " + "z" * 5000])
        labels = sc._fragment_backbone_marker_labels(frag)
        assert len(labels) == 1
        assert len(labels[0]) <= 80

    def test_non_string_label_does_not_raise(self):
        frag = {"features": [{"type": "gene", "label": 12345},
                             {"type": "gene", "label": ["AmpR"]}]}
        # Must not raise; coerces via str().
        sc._fragment_backbone_marker_labels(frag)
        sc._fragment_has_backbone_marker(frag)

    def test_duplicates_collapse(self):
        frag = _frag(100, ["AmpR", "AmpR", "AmpR"])
        assert sc._fragment_backbone_marker_labels(frag) == ["AmpR"]


class TestOriginSpanningFeatureSurvivesTheClone:
    """2026-08-26 report #4 (the reporter's prediction that
    `carry_annotations` had the same origin-wrap hole as `list-features`,
    which it did).

    A feature that wraps bp 0 is stored `end < start` and `_split_features_
    at_cuts` HAS to halve it at the origin before slotting. Nothing put the
    halves back, so a 25 bp T-DNA right border nowhere near the cloning site
    came out of the clone as two adjacent 12 bp features with the same name —
    bases right, annotation saying the border was in pieces, which is exactly
    what an "is my border intact" check looks at."""

    PRE  = "TTGACCATGCATGGCCTAAG" * 5     # 100, holds bp 0
    MID  = "ACGTTAAGCCTTAAGGCCAA" * 3     # 60, between the two cuts
    POST = "AGGCCTTACGTACGTTAAGC" * 5     # 100

    def _vector(self):
        return self.PRE + "GAATTC" + self.MID + "GGATCC" + self.POST

    def _record(self, seq):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        n = len(seq)
        rec = SeqRecord(Seq(seq), id="vecW", name="vecW",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        border = SeqFeature(
            CompoundLocation([FeatureLocation(n - 20, n, 1),
                              FeatureLocation(0, 20, 1)]), type="misc_feature")
        border.qualifiers["label"] = ["RB border"]
        kan = SeqFeature(FeatureLocation(40, 90, 1), type="CDS")
        kan.qualifiers["label"] = ["KanR"]
        rec.features = [border, kan]
        return rec

    def _clone(self):
        vec = self._vector()
        sc._save_library([{"name": "vecW", "id": "vecW",
                           "gb_text": sc._record_to_gb_text(
                               self._record(vec))}])
        ins = "GAATTC" + "ACGTACGTAAGGCCTTAACC" * 4 + "GGATCC"
        r = sc._h_traditional_clone(None, {
            "vector_seq": vec, "insert_seq": ins,
            "vector_enzymes": ["EcoRI", "BamHI"],
            "insert_enzymes": ["EcoRI", "BamHI"],
            "vector_circular": True, "carry_annotations": True,
            "vector_name": "vecW", "product_name": "clonedW"})
        assert not isinstance(r, tuple), r
        ent = next(e for e in sc._load_library() if e["name"] == r["saved_name"])
        return vec, sc._gb_text_to_record(ent["gb_text"])

    def test_the_border_comes_through_as_one_feature(self):
        _vec, rec = self._clone()
        borders = [f for f in rec.features
                   if (f.qualifiers.get("label") or [""])[0] == "RB border"]
        assert len(borders) == 1, [
            (f.qualifiers.get("label"), sc._feat_bounds(f, len(rec.seq)))
            for f in rec.features]

    def test_and_it_still_spells_the_original_bases(self):
        """The coordinates being tidy is not the point — the annotated span
        has to be the same 40 bases the vector had."""
        vec, rec = self._clone()
        n_v, total = len(vec), len(rec.seq)
        want = vec[n_v - 20:] + vec[:20]
        border = next(f for f in rec.features
                      if (f.qualifiers.get("label") or [""])[0] == "RB border")
        s, e, _strand = sc._feat_bounds(border, total)
        got = (str(rec.seq[s:] + rec.seq[:e]) if e < s
               else str(rec.seq[s:e]))
        assert got == want
        assert sc._feat_len(s, e, total) == 40

    def test_a_border_the_cloning_really_cut_stays_two_pieces(self):
        """The rejoin must not paper over a genuine disruption: when a cut
        lands INSIDE the wrapping feature, two flanking pieces IS the truth,
        and the surviving half is labelled disrupted."""
        vec = "GAATTC" + self.PRE + "GGATCC" + self.POST   # EcoRI cuts at bp 1
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
        n = len(vec)
        rec = SeqRecord(Seq(vec), id="vecX", name="vecX",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        border = SeqFeature(
            CompoundLocation([FeatureLocation(n - 12, n, 1),
                              FeatureLocation(0, 4, 1)]), type="misc_feature")
        border.qualifiers["label"] = ["RB border"]
        rec.features = [border]
        sc._save_library([{"name": "vecX", "id": "vecX",
                           "gb_text": sc._record_to_gb_text(rec)}])
        ins = "GAATTC" + "ACGTACGTAAGGCCTTAACC" * 4 + "GGATCC"
        r = sc._h_traditional_clone(None, {
            "vector_seq": vec, "insert_seq": ins,
            "vector_enzymes": ["EcoRI", "BamHI"],
            "insert_enzymes": ["EcoRI", "BamHI"],
            "vector_circular": True, "carry_annotations": True,
            "vector_name": "vecX", "product_name": "clonedX",
            "vector_frag_idx": 0})
        assert not isinstance(r, tuple), r
        ent = next(e for e in sc._load_library() if e["name"] == r["saved_name"])
        prod = sc._gb_text_to_record(ent["gb_text"])
        labels = [(f.qualifiers.get("label") or [""])[0] for f in prod.features]
        assert any("RB border" in lb and "(disrupted)" in lb for lb in labels), \
            labels


class TestRejoinOriginSplitFeatures:
    """Unit-level cover for `_rejoin_origin_split_features` — the branches the
    end-to-end clone test can't reach on demand."""

    @staticmethod
    def _half(gid, which, start, end, label="RB border"):
        return {"start": start, "end": end, "strand": 1,
                "type": "misc_feature", "label": label,
                "_wrap_origin_split": which, "_wrap_origin_group": gid}

    def test_adjacent_halves_merge(self):
        feats = [self._half("w0", "tail", 85, 105),
                 self._half("w0", "head", 105, 125)]
        out = sc._rejoin_origin_split_features(feats, 292)
        assert len(out) == 1
        assert (out[0]["start"], out[0]["end"]) == (85, 125)
        # Tags are consumed — a rejoined feature is no longer "split".
        assert "_wrap_origin_split" not in out[0]
        assert "_wrap_origin_group" not in out[0]

    def test_halves_meeting_at_the_products_own_origin_stay_a_wrap(self):
        """Ligation can land the pieces either side of bp 0 of the NEW
        molecule; then the rejoined feature wraps, `end < start` — the same
        shape it had on the vector (sacred invariant #8)."""
        feats = [self._half("w0", "tail", 280, 292),
                 self._half("w0", "head", 0, 8)]
        out = sc._rejoin_origin_split_features(feats, 292)
        assert len(out) == 1
        assert (out[0]["start"], out[0]["end"]) == (280, 8)
        assert sc._feat_len(280, 8, 292) == 20

    def test_separated_halves_are_left_alone(self):
        """Cloning that genuinely pulled the two apart must keep reporting
        two pieces — merging them would invent a feature spanning an insert
        that was dropped in between."""
        feats = [self._half("w0", "tail", 10, 30),
                 self._half("w0", "head", 200, 220)]
        out = sc._rejoin_origin_split_features(feats, 292)
        assert len(out) == 2

    def test_two_wrapping_features_do_not_cross_pair(self):
        """The group id is what stops one feature's head merging onto
        another's tail — the failure a bare head/tail tag couldn't prevent."""
        feats = [self._half("w0", "tail", 0, 10, "RB"),
                 self._half("w1", "head", 10, 20, "LB"),
                 self._half("w0", "head", 40, 50, "RB"),
                 self._half("w1", "tail", 50, 60, "LB")]
        out = sc._rejoin_origin_split_features(feats, 100)
        spans = sorted((f["start"], f["end"]) for f in out)
        assert spans == [(0, 10), (10, 20), (40, 50), (50, 60)]

    def test_untagged_features_pass_through_untouched(self):
        feats = [{"start": 1, "end": 2, "type": "CDS", "label": "x"}]
        assert sc._rejoin_origin_split_features(feats, 100) == feats

    def test_a_three_way_split_is_left_alone(self):
        """A half that was cut AGAIN gives the group three members. Merging
        two of them would invent a span across the piece that was excised —
        so an ambiguous group is not merged at all."""
        feats = [self._half("w0", "tail", 0, 10),
                 self._half("w0", "head", 10, 20),
                 self._half("w0", "head", 20, 30)]
        assert len(sc._rejoin_origin_split_features(feats, 100)) == 3

    def test_a_lone_half_passes_through(self):
        """The other half went onto the fragment that was discarded — one
        piece is the truth."""
        feats = [self._half("w0", "tail", 0, 10)]
        out = sc._rejoin_origin_split_features(feats, 100)
        assert len(out) == 1 and out[0]["_wrap_origin_split"] == "tail"

    def test_opposite_strands_never_merge(self):
        a = self._half("w0", "tail", 0, 10)
        b = self._half("w0", "head", 10, 20)
        b["strand"] = -1
        assert len(sc._rejoin_origin_split_features([a, b], 100)) == 2

    def test_a_non_integer_coordinate_does_not_raise(self):
        """A hand-edited GenBank or a caller-built feature list can carry a
        fuzzy coordinate; it must drop out of the merge, not out of the
        cloning run."""
        a = self._half("w0", "tail", 0, 10)
        b = self._half("w0", "head", 10, 20)
        a["start"] = "<1"
        assert len(sc._rejoin_origin_split_features([a, b], 100)) == 2
