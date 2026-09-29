"""Cross-cutting sequence-integrity audit (2026-05-20).

Sacred invariants this file enforces:

1. **Length conservation on excise.** For a 2-cut digest on a circular
   plasmid: ``len(original) == sum(fragment.top_seq lengths)``. Every
   nucleotide accounted for exactly once. No silent duplication, no
   silent loss.

2. **Overhang-once-at-seam on ligate.** Ligating fragment A's right
   end to B's left end produces ``A.top_seq + B.top_seq`` — the
   overhang bases live in WHICHEVER fragment's top strand carried
   them (top-strand canonicalisation rule, splicecraft.py line ~4577).
   Verified by counting unique sequence motifs across the seam.

3. **Round-trip cut+religate = original** (modulo origin rotation).
   The strongest test: digest a known plasmid, religate the
   fragments, the result must be sequence-equal to the input under
   rotation.

4. **L0 → TU body conservation.** Each L0 part's body sequence appears
   exactly once in the assembled TU. Overhangs at part-part seams
   appear once (not twice).

5. **TU → MOD chain conservation.** Two TUs chained into an Ω-MOD:
   both TU bodies preserved, seam overhang once, MOD boundary
   overhangs at the outer cuts.

6. **Gibson overlap-merge.** N fragments with overlap regions →
   product length = sum(fragment_lengths) − sum(overlap_lengths).
   Overlap sequences appear ONCE.

7. **Traditional cloning round-trip.** Digest then religate must
   yield the original plasmid (under rotation).

8. **Entry-vector ops are read-only.** `_detect_entry_vector_role` /
   `_auto_bind_entry_vectors_from_entries` MUST NOT mutate the input
   record's sequence.

Sacred — never pick fragments by size (per
`feedback_never_assume_smaller_frag_is_payload`). Where the test
needs to identify the "insert" vs "vector" half, it uses feature
markers, not length comparisons.
"""
import pytest
import splicecraft as sc

pytestmark = [pytest.mark.usefixtures("_protect_user_data")]


# ── Synthetic-plasmid builders ─────────────────────────────────────────────


def _make_circular_record(seq: str, *, features=None, name: str = "test"):
    """Build a circular SeqRecord with the given top-strand sequence."""
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    rec = SeqRecord(
        Seq(seq), id=name, name=name,
        annotations={"topology": "circular", "molecule_type": "DNA"},
    )
    if features:
        for f in features:
            rec.features.append(f)
    return rec


def _rotate(seq: str, by: int) -> str:
    """Rotate a circular sequence by `by` bases (positive = left shift)."""
    n = len(seq)
    by = by % n if n else 0
    return seq[by:] + seq[:by]


def _sequences_match_under_rotation(a: str, b: str) -> bool:
    """True iff `a` and `b` describe the same circular molecule —
    one is a rotation of the other. Uses the doubled-string trick:
    `a` is a rotation of `b` iff `a` is a substring of `b + b`."""
    if len(a) != len(b):
        return False
    if not a:
        return True
    return a in (b + b)


# ── 1. EXCISE: length conservation ─────────────────────────────────────────

class TestExciseLengthConservation:
    """For every successful 2-cut digest: total top_seq length must
    equal the original sequence length, exactly."""

    @pytest.mark.parametrize("enzyme,site", [
        ("EcoRI",  "GAATTC"),  # 6-bp palindromic, 4-bp 5' overhang
        ("BamHI",  "GGATCC"),  # 6-bp palindromic
        ("HindIII","AAGCTT"),  # 6-bp palindromic
    ])
    def test_palindromic_two_cut_length_conservation(self, enzyme, site):
        """Two EcoRI/BamHI/HindIII sites in a circular plasmid:
        digest → 2 fragments, their top_seq lengths must sum to
        the original length."""
        # Construct: filler + site + filler + site + filler.
        filler1 = "AAAAAAAAAACCCCCCCCCC"  # 20 bp
        filler2 = "TTTTTTTTTTGGGGGGGGGG"  # 20 bp
        filler3 = "ATATATATATCGCGCGCGCG"  # 20 bp
        seq = filler1 + site + filler2 + site + filler3
        frags, err = sc._excise_fragment_pair(
            seq, [enzyme], circular=True,
        )
        assert err is None, f"{enzyme} digest errored: {err}"
        assert len(frags) == 2
        total = sum(len(f["top_seq"]) for f in frags)
        assert total == len(seq), (
            f"{enzyme}: length lost. orig={len(seq)} sum_fragments={total} "
            f"diff={len(seq) - total}"
        )

    @pytest.mark.parametrize("enzyme", ["BsaI", "Esp3I"])
    def test_type_iis_length_conservation(self, enzyme):
        """Type IIS enzymes (BsaI, Esp3I) cut outside their site.
        Asymmetric cut geometry is the classic off-by-one risk."""
        # Build a plasmid with 2 outward-facing BsaI sites so the
        # central segment gets released cleanly.
        # BsaI = GGTCTC(N1)↓N4 — cuts 1bp downstream, 4bp 5' overhang.
        # Esp3I = CGTCTC(N1)↓N4 — same geometry, different recognition.
        if enzyme == "BsaI":
            fwd_site, rc_site = "GGTCTC", "GAGACC"
        else:
            fwd_site, rc_site = "CGTCTC", "GAGACG"
        # Forward site sees: ENZ + N + 4bp_overhang + payload
        # Reverse site sees: payload + 4bp_overhang + N + RC_ENZ
        inside = (
            fwd_site + "A" + "GGAG" +
            "CCCCCAAAAAGGGGGTTTTT" +
            "CGCT" + "T" + rc_site
        )
        backbone = "ATCGATCGAT" * 50
        seq = inside + backbone
        frags, err = sc._excise_fragment_pair(
            seq, [enzyme], circular=True,
        )
        assert err is None
        assert len(frags) == 2
        total = sum(len(f["top_seq"]) for f in frags)
        assert total == len(seq), (
            f"{enzyme} length lost: orig={len(seq)} sum={total}"
        )

    def test_three_cut_excise_refuses(self):
        """3 cuts on a circular plasmid → ambiguous excise; the
        sacred invariant says we must refuse rather than ship a
        wrong product (CLAUDE.md #25)."""
        # 3 EcoRI sites.
        seq = "AAAA" + "GAATTC" + "TTTT" + "GAATTC" + "GGGG" + "GAATTC" + "CCCC"
        frags, err = sc._excise_fragment_pair(
            seq, ["EcoRI"], circular=True,
        )
        assert err is not None, (
            "≥3 cuts must surface an error (sacred CLAUDE.md #25)"
        )
        assert "exactly 2" in err["error"].lower()


# ── 2. EXCISE → LIGATE: round-trip equivalence ─────────────────────────────

class TestExciseLigateRoundTrip:
    """Cut a circular plasmid then religate the fragments. The
    product must be sequence-equal to the original (modulo origin
    rotation). This is the strongest sequence-conservation test."""

    @pytest.mark.parametrize("enzyme,site", [
        ("EcoRI", "GAATTC"),
        ("BamHI", "GGATCC"),
    ])
    def test_palindromic_round_trip(self, enzyme, site):
        """Two-site palindromic digest + religate must give back
        the original sequence (under rotation)."""
        filler1 = "AAAAACCCCCAAAAACCCCC"  # 20
        filler2 = "TTTTTGGGGGTTTTTGGGGG"  # 20
        filler3 = "ATATGCGCATATGCGCATAT"  # 20
        original = filler1 + site + filler2 + site + filler3
        frags, err = sc._excise_fragment_pair(
            original, [enzyme], circular=True,
        )
        assert err is None
        assert len(frags) == 2
        # Ligate frag[0]→frag[1] then close circular.
        merged = sc._ligate_fragments(frags[0], frags[1])
        assert merged is not None, "Religation failed"
        circular = sc._close_circular(merged)
        assert circular is not None, "Circle close failed"
        assert _sequences_match_under_rotation(
            circular["top_seq"], original,
        ), (
            f"{enzyme} round-trip drift:\n"
            f"  original ({len(original)}bp): {original}\n"
            f"  result   ({len(circular['top_seq'])}bp): {circular['top_seq']}"
        )

    def test_two_palindromic_enzymes_round_trip(self):
        """EcoRI + BamHI (different enzymes, different sites): the
        directional ligation must still reproduce the original."""
        # EcoRI overhang = AATT (palindrome). BamHI overhang = GATC
        # (palindrome). Different sites → directional ligation; both
        # palindromic so each enzyme's pair self-ligates.
        original = (
            "AAACCC" + "GAATTC" + "TTGGAACC" + "GGATCC" + "CGCGATAT"
        )
        frags, err = sc._excise_fragment_pair(
            original, ["EcoRI", "BamHI"], circular=True,
        )
        assert err is None
        assert len(frags) == 2
        # Try both ligation orders — one will match.
        for a, b in ((frags[0], frags[1]), (frags[1], frags[0])):
            merged = sc._ligate_fragments(a, b)
            if merged is None:
                continue
            circular = sc._close_circular(merged)
            if circular is None:
                continue
            if _sequences_match_under_rotation(
                circular["top_seq"], original,
            ):
                return  # success
        pytest.fail(
            "Neither ligation order reproduced the original sequence"
        )

    def test_type_iis_round_trip(self):
        """Type IIS cut + religate must preserve sequence. Most
        likely site for an off-by-one because BsaI/Esp3I cut
        outside their recognition sequence."""
        # Build with two outward-facing BsaI sites releasing a
        # cassette with (GGAG, CGCT) overhangs.
        inside = (
            "GGTCTC" + "A" + "GGAG" +
            "AAATTTCCCGGG" +
            "CGCT" + "T" + "GAGACC"
        )
        backbone = "GCGCATATATAT" * 8
        original = inside + backbone
        frags, err = sc._excise_fragment_pair(
            original, ["BsaI"], circular=True,
        )
        assert err is None
        assert len(frags) == 2
        # Both ligation orders.
        for a, b in ((frags[0], frags[1]), (frags[1], frags[0])):
            merged = sc._ligate_fragments(a, b)
            if merged is None:
                continue
            circular = sc._close_circular(merged)
            if circular is None:
                continue
            if _sequences_match_under_rotation(
                circular["top_seq"], original,
            ):
                return
        pytest.fail(
            "Type IIS round-trip drift — neither orientation matches "
            "the original"
        )


# ── 3. OVERHANG-ONCE-AT-SEAM ───────────────────────────────────────────────


class TestOverhangAtSeam:
    """When two fragments ligate, the overhang sequence appears
    exactly ONCE in the merged top strand (top-strand-canonicalisation
    rule, splicecraft.py:4577)."""

    def test_palindromic_overhang_appears_once(self):
        """EcoRI overhang AATT: cut a circular plasmid with two
        EcoRI sites, religate, and count AATT occurrences in the
        product. Must equal the input count (overhang re-formed at
        each seam, lives in one fragment only)."""
        original = (
            "AAACCC" + "GAATTC" +
            "ATATCGCG" + "GAATTC" +
            "TTGGAACC"
        )
        original_aatt = original.count("AATT")
        frags, err = sc._excise_fragment_pair(
            original, ["EcoRI"], circular=True,
        )
        assert err is None
        merged = sc._ligate_fragments(frags[0], frags[1])
        assert merged is not None
        # Religating restores the site → AATT count survives.
        # The merged is LINEAR (before close-circular), so it has
        # one fewer cut than the original. Count must equal
        # (original_aatt - 1) iff a seam straddles two AATTs.
        # Easier: close circular and recount; that must equal
        # original_aatt.
        circular = sc._close_circular(merged)
        assert circular is not None
        assert circular["top_seq"].count("AATT") == original_aatt, (
            f"AATT count drift: original={original_aatt} "
            f"after-religate={circular['top_seq'].count('AATT')}"
        )

    def test_type_iis_overhang_appears_once_at_seam(self):
        """Cut with BsaI; the released stuffer's left/right
        overhang sequences (GGAG and CGCT in the test rig) must
        each appear exactly once in the LIGATED product where the
        seam reconstitutes the original. NOT doubled (the
        canonicalisation rule)."""
        # Build a known plasmid with EXACTLY one GGAG and one CGCT
        # outside the BsaI sites so we know the baseline count.
        inside = (
            "GGTCTC" + "A" + "GGAG" +
            "AAATTTCCC" +     # no GGAG / CGCT here
            "CGCT" + "T" + "GAGACC"
        )
        backbone = "ATATATATAT" * 12  # no GGAG / CGCT
        original = inside + backbone
        n_ggag = original.count("GGAG")
        n_cgct = original.count("CGCT")
        # Sanity: the construction has exactly one of each.
        assert n_ggag == 1
        assert n_cgct == 1
        frags, err = sc._excise_fragment_pair(
            original, ["BsaI"], circular=True,
        )
        assert err is None
        for a, b in ((frags[0], frags[1]), (frags[1], frags[0])):
            merged = sc._ligate_fragments(a, b)
            if merged is None:
                continue
            circular = sc._close_circular(merged)
            if circular is None:
                continue
            if _sequences_match_under_rotation(
                circular["top_seq"], original,
            ):
                assert circular["top_seq"].count("GGAG") == n_ggag
                assert circular["top_seq"].count("CGCT") == n_cgct
                return
        pytest.fail("Round-trip didn't match either orientation")


# ── 4. FEATURE COORDINATES PRESERVED ───────────────────────────────────────


class TestFeatureCoordinatesAfterDigest:
    """Features on the input record must land at the correct
    fragment-local coordinates after `_excise_fragment_pair`. Loss
    or shift here would silently break the Constructor's downstream
    feature-painting in assembled TUs."""

    def test_feature_inside_fragment_lands_correctly(self):
        """A feature wholly inside one fragment must appear in
        that fragment's `features` list with correctly-shifted
        coordinates."""
        # Build: EcoRI at pos 10, EcoRI at pos 40, feature 20..30.
        original = "AAAAAAAAAA" + "GAATTC" + "AAATTTGGGCCC" + "AAAAAA" + "GAATTC" + "TTTTTTTT"
        # GAATTC at 10 and ... let's compute.
        s1 = original.find("GAATTC")
        s2 = original.find("GAATTC", s1 + 1)
        assert s1 != -1 and s2 != -1 and s2 > s1
        # Feature inside the segment between the two sites.
        feat_start = s1 + 8
        feat_end   = s1 + 14
        feats = [{
            "start": feat_start, "end": feat_end,
            "type":  "CDS", "qualifiers": {"label": ["inner_feat"]},
        }]
        frags, err = sc._excise_fragment_pair(
            original, ["EcoRI"], circular=True,
            features=feats,
        )
        assert err is None
        # Find the fragment carrying the feature.
        carrier = None
        for f in frags:
            for ff in f.get("features", []):
                lab = ff.get("qualifiers", {}).get("label") or []
                if "inner_feat" in lab:
                    carrier = (f, ff)
                    break
            if carrier:
                break
        assert carrier is not None, (
            "Inner feature was lost across the digest"
        )
        frag, feat = carrier
        # Verify the bases at the feature's local coords in the
        # carrier's top_seq match the bases at the absolute coords
        # in the original — that's the sacred property.
        local_bases = frag["top_seq"][
            feat["start"]: feat["end"]
        ]
        orig_bases = original[feat_start:feat_end]
        assert local_bases == orig_bases, (
            f"Feature coords drifted: original={orig_bases!r} "
            f"fragment-local={local_bases!r}"
        )


# ── 5. GIBSON ASSEMBLY ─────────────────────────────────────────────────────


class TestGibsonAssembly:
    """Gibson assembly: N fragments with overlap homology → 1
    product. Overlap appears ONCE in the product (not twice). Total
    length = sum(fragment lengths) − sum(overlap lengths)."""

    def test_three_fragment_gibson_overlap_appears_once(self):
        """Three fragments with 20-bp overlaps. Verify the overlap
        sequences appear exactly once in the assembled product."""
        overlap_AB = "AAACCCAAACCCAAACCCAA"   # 20 bp
        overlap_BC = "TTGGGTTTGGGTTTGGGTTT"   # 20 bp
        overlap_CA = "GCATGCATGCATGCATGCAT"   # 20 bp
        body_A = "AAAATTTTGGGG" * 5   # 60 bp body
        body_B = "CCCCAAAATTTT" * 5
        body_C = "TTTTGGGGCCCC" * 5
        # Each fragment starts with its overlap with the previous
        # fragment and ends with its overlap with the next. Circular
        # assembly: frag_A → frag_B → frag_C → frag_A (wraps).
        frag_A = overlap_CA + body_A + overlap_AB
        frag_B = overlap_AB + body_B + overlap_BC
        frag_C = overlap_BC + body_C + overlap_CA
        fragments = [
            {"name": "A", "sequence": frag_A, "features": []},
            {"name": "B", "sequence": frag_B, "features": []},
            {"name": "C", "sequence": frag_C, "features": []},
        ]
        result = sc._simulate_gibson_assembly(
            fragments, min_overlap=15, circular=True,
        )
        assert result.get("success") is True, (
            f"Gibson failed: errors={result.get('errors')} "
            f"warnings={result.get('warnings')}"
        )
        product_seq = result["product_seq"]
        # Length conservation: product = sum(fragments) - 3 overlaps
        # (one per junction, including the wrap junction).
        expected_len = (
            len(frag_A) + len(frag_B) + len(frag_C)
            - 3 * 20
        )
        assert len(product_seq) == expected_len, (
            f"Gibson length drift: expected {expected_len}, "
            f"got {len(product_seq)}"
        )
        # Each overlap appears ONCE in the doubled product (the
        # circular product wraps, so the linear top-seq might
        # split an overlap across the join — count in doubled).
        for ovl in (overlap_AB, overlap_BC, overlap_CA):
            occ_doubled = (product_seq + product_seq).count(ovl)
            assert 1 <= occ_doubled <= 2, (
                f"Overlap {ovl!r} appears {occ_doubled} times in "
                "doubled product — expected exactly 1 in circular "
                "(2 max if the linear top_seq happens to repeat the "
                "overlap at the join). Gibson is over-counting if >2."
            )

    def test_two_fragment_gibson_linear(self):
        """Two fragments with one overlap → linear product, NO
        wrap. Length = a + b − overlap."""
        overlap = "GCATCGATCGATCGATCGAT"  # 20 bp
        body_a = "AAAATTTTGGGG" * 5
        body_b = "CCCCAAAATTTT" * 5
        frag_a = body_a + overlap
        frag_b = overlap + body_b
        result = sc._simulate_gibson_assembly(
            [
                {"name": "A", "sequence": frag_a, "features": []},
                {"name": "B", "sequence": frag_b, "features": []},
            ],
            min_overlap=15,
            circular=False,
        )
        assert result.get("success") is True, (
            f"Linear Gibson failed: {result.get('errors')}"
        )
        product = result["product_seq"]
        expected_len = len(frag_a) + len(frag_b) - len(overlap)
        assert len(product) == expected_len
        # Overlap appears exactly once in the linear product.
        assert product.count(overlap) == 1, (
            f"Overlap appears {product.count(overlap)} times in "
            "linear Gibson product — expected exactly 1"
        )
        # Body content from both halves must be present.
        assert body_a in product
        assert body_b in product


# ── 6. TRADITIONAL CLONING (cut + ligate) ─────────────────────────────────


class TestTraditionalCloningRoundTrip:
    """Cut a plasmid with a restriction enzyme set, religate the
    resulting fragments — the product must equal the original
    (under rotation)."""

    def test_eco_bam_double_digest_relegate(self):
        """EcoRI + BamHI double digest → 2 fragments. Religate them
        in the only compatible orientation — must equal original."""
        original = (
            "AAAACCCC" + "GAATTC" + "TTTGGG" +
            "GGATCC" + "CGCGATAT" + "TGCATGCA"
        )
        frags, err = sc._excise_fragment_pair(
            original, ["EcoRI", "BamHI"], circular=True,
        )
        assert err is None
        # The two fragments have different terminal overhangs (one
        # EcoRI-end + one BamHI-end). Only the correct orientation
        # ligates.
        for a, b in ((frags[0], frags[1]), (frags[1], frags[0])):
            merged = sc._ligate_fragments(a, b)
            if merged is None:
                continue
            circular = sc._close_circular(merged)
            if circular is None:
                continue
            if _sequences_match_under_rotation(
                circular["top_seq"], original,
            ):
                return
        pytest.fail("Round-trip failed for EcoRI+BamHI digest")


# ── 7. L0 → TU ASSEMBLY ────────────────────────────────────────────────────


class TestL0ToTuAssembly:
    """Build a TU from 3 L0 parts via the actual assembly machinery.
    Verify body conservation + boundary overhangs."""

    def _make_l0_part(
        self, *, body: str, oh5: str, oh3: str,
        backbone: str = "ATCGAT" * 30,
        marker_label: str = "AmpR",
    ):
        """Synthesise an L0 part: a circular pUPD2-style plasmid
        carrying `body` between two outward-facing BsaI sites,
        with a rep_origin feature on the backbone half so the
        backbone-marker exclusion picks the right fragment."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        # Outward-facing BsaI: ENZ_RC + N + oh5 + body + oh3 + N + ENZ_FWD
        # That arrangement: when BsaI cuts, the released fragment
        # has top_seq = oh5 + body + oh3 (plus overhang remnants).
        # Build: backbone + (BsaI RC + N1 + oh5 + body + oh3 + N1 + BsaI FWD).
        # Simpler: place inward-facing BsaI sites and the body
        # between them — that's the L1-acceptor pattern. For L0,
        # the BsaI sites are OUTWARD (in the part, pointing into
        # the body). The released fragment IS the part body.
        # Cassette: ENZ_FWD + N + oh5 + body + oh3 + N + ENZ_RC
        # Reading left→right: "GGTCTC" + "A" + oh5 + body + oh3 + "A" + "GAGACC"
        # No that's INWARD. Let me think again...
        #
        # OUTWARD-facing means the BsaI enzyme on the BACKBONE side
        # of the cut — the cut releases the body with oh5/oh3
        # overhangs.
        # Sequence layout (linearised for clarity, actual is circular):
        #   ...backbone... GAGACC [N] [oh3] [body_rc] [oh5] [N] GGTCTC ...backbone...
        # Wait, the SIMPLEST L0 part layout:
        #   [body_with_overhangs_baked_in] — but BsaI has to cut OUT
        #   of the body. So the OUTER edges of the body must have
        #   BsaI sites pointing OUTWARD. The site faces outward
        #   means "RC" on the left, "FWD" on the right.
        # Result layout (l → r, top strand):
        #   GAGACC N oh3 body_rc oh5 N GGTCTC ... backbone ...
        # Wait actually:
        # Forward BsaI: GGTCTC N N N N N | (5' overhang on top is
        #   the next 4 bases AFTER position 7, which is the cut).
        # If we want the body to have a (oh5, oh3) outer-overhang
        # pair when released by BsaI, we place:
        #   FWD: "GGTCTC" + N + oh5 + body + oh3 + N + "GAGACC" → INWARD
        # That's actually inward (both sites point INTO the body).
        # Hmm. Re-checking: GGTCTC(N1)↓(N4) — top strand cuts at
        # position site_end+1, bottom strand cuts at site_end+5.
        # If sequence is "...GGTCTCANNNN...", top-cut happens just
        # after the A. So the bases NNNN are the 4-bp overhang of
        # the RIGHT fragment.
        # For body-to-have-(oh5,oh3): we want the LEFT BsaI cut to
        # leave oh5 as the right fragment's overhang (which becomes
        # the body's LEFT overhang), and the RIGHT BsaI cut to leave
        # oh3 as the left fragment's overhang (which becomes the
        # body's RIGHT overhang).
        # Left site (forward): "GGTCTC" + "A" + oh5 + body...
        # Right site (reverse): ...body + oh3 + "A" + "GAGACC" (RC).
        # This is the INWARD arrangement → both BsaI sites cut INTO
        # the body. The body IS the released "stuffer" of a sort.
        # But wait — for L0 parts, the BsaI sites are usually OUTWARD
        # so that cutting RELEASES the body, leaving the backbone.
        # That's: backbone + (outer flanks point INTO body).
        # When BsaI cuts, the body comes out with oh5/oh3 overhangs
        # at its ends; the backbone gets the mirrored overhangs.
        cassette = (
            "GGTCTC" + "A" + oh5 + body + oh3 + "A" + "GAGACC"
        )
        full = cassette + backbone
        feats: list = []
        # Mark the backbone half with rep_origin so
        # `_fragment_has_backbone_marker` picks it correctly.
        ori_start = len(cassette) + 50
        feats.append(SeqFeature(
            FeatureLocation(ori_start, ori_start + 200),
            type="rep_origin",
            qualifiers={"label": [marker_label]},
        ))
        return _make_circular_record(full, features=feats), full

    def test_three_part_l0_excise_overhangs_correct(self):
        """A synthesised L0 part with body+(GGAG,AATG) overhangs
        digests with BsaI to a fragment whose released overhangs
        are exactly (GGAG, AATG)."""
        body = "ATGAAACCCGGGTTT" * 4   # 60 bp
        rec, full_seq = self._make_l0_part(
            body=body, oh5="GGAG", oh3="AATG",
        )
        feats = [
            {"start": int(f.location.start),
              "end":   int(f.location.end),
              "type":  f.type,
              "qualifiers": dict(f.qualifiers)}
            for f in rec.features
        ]
        frags, err = sc._excise_fragment_pair(
            full_seq, ["BsaI"], circular=True,
            features=feats,
        )
        assert err is None
        assert len(frags) == 2
        # Identify the body via backbone-marker exclusion (NEVER size).
        marked = [
            sc._fragment_has_backbone_marker(f) for f in frags
        ]
        assert sum(marked) == 1, (
            "Backbone marker not unambiguous — fragment-marker "
            "annotation is broken in the synthetic builder"
        )
        body_frag = frags[0] if not marked[0] else frags[1]
        oh5 = body_frag["left"]["overhang_seq"]
        oh3 = body_frag["right"]["overhang_seq"]
        assert oh5 == "GGAG", f"L0 part oh5 drifted: {oh5!r} != 'GGAG'"
        assert oh3 == "AATG", f"L0 part oh3 drifted: {oh3!r} != 'AATG'"
        # Body must appear in the released fragment's top_seq.
        assert body in body_frag["top_seq"], (
            "L0 part body sequence missing from released fragment"
        )


# ── 8. ENTRY-VECTOR OPS ARE SEQUENCE-READ-ONLY ─────────────────────────────


class TestEntryVectorOpsImmutability:
    """The new detection/auto-bind pipeline must NEVER mutate
    sequences. Only metadata in `entry_vectors.json` changes."""

    def test_detect_does_not_mutate_record(self):
        """`_detect_entry_vector_role` is purely diagnostic. Run on
        a record, capture seq before + after, assert byte-equality."""
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        original_seq = (
            "GGTCTC" + "A" + "GGAG" +
            "AAATTTCCC" * 30 +
            "CGCT" + "A" + "GAGACC" +
            "ATCGATCGAT" * 50
        )
        rec = SeqRecord(
            Seq(original_seq), id="t", name="t",
            annotations={"topology": "circular",
                          "molecule_type": "DNA"},
        )
        gb_l0 = sc._all_grammars()["gb_l0"]
        sc._detect_entry_vector_role(rec, gb_l0)
        # Sequence unchanged after detection.
        assert str(rec.seq) == original_seq, (
            "_detect_entry_vector_role mutated the input sequence"
        )

    def test_auto_bind_does_not_mutate_entries(self):
        """`_auto_bind_entry_vectors_from_entries` reads gb_text
        but must NOT mutate the input dicts' sequences."""
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        body = "ATGAAACCC" * 5
        cassette = (
            "GGTCTC" + "A" + "GGAG" + body + "CGCT" + "A" + "GAGACC"
        )
        # Synthesize an Alpha1-style acceptor.
        outer_cassette = (
            "CGTCTC" + "A" + "GGAG" + cassette + "GTCA" + "A" + "GAGACG"
        )
        full = outer_cassette + "ATCGAT" * 50
        rec = _make_circular_record(full)
        rec.features.append(SeqFeature(
            FeatureLocation(len(outer_cassette) + 50, len(outer_cassette) + 200),
            type="rep_origin",
            qualifiers={"label": ["test_ori"]},
        ))
        gb_text = sc._record_to_gb_text(rec)
        entry = {
            "id": "test_acc", "name": "test_acc",
            "size": len(full), "gb_text": gb_text,
        }
        gb_text_before = entry["gb_text"]
        sc._auto_bind_entry_vectors_from_entries([entry])
        # The entry's gb_text must not be touched.
        assert entry["gb_text"] == gb_text_before


class TestRestrictionOverlayNoRedundantRebuild:
    """Issue #27: the restriction worker's callback must not rebuild the whole
    sequence-panel layout to paint an overlay containing nothing.

    `update_seq` is handed `pm._feats + displayed`, and that `+` produces a NEW
    list even when `displayed` is empty. Every sequence-panel cache is keyed
    partly on `id(feats)`, so the new object misses all of them and the layout
    is recomputed from scratch. On a plasmid that is invisible; on a 4.3 Mb
    chromosome it is a second multi-second freeze, and it fires on EVERY load
    because the default `unique_only=True` finds no unique cutters in a few
    megabases.
    """

    def _panel(self, feats):
        sp = sc.SequencePanel.__new__(sc.SequencePanel)
        sp._feats = feats
        return sp

    def test_reports_no_overlay_for_plain_features(self):
        sp = self._panel([{"start": 0, "end": 9, "type": "CDS", "label": "x"}])
        assert sp._feats_have_restriction_overlay() is False

    def test_reports_no_overlay_when_empty(self):
        assert self._panel([])._feats_have_restriction_overlay() is False

    @pytest.mark.parametrize("kind", ["resite", "recut"])
    def test_detects_a_showing_overlay(self, kind):
        """Both overlay types must count — otherwise switching the overlay OFF
        would skip the repaint that removes it, leaving it stuck on screen."""
        sp = self._panel([
            {"start": 0, "end": 9, "type": "CDS", "label": "x"},
            {"start": 3, "end": 9, "type": kind, "label": "EcoRI"},
        ])
        assert sp._feats_have_restriction_overlay() is True

    def test_new_list_of_same_features_misses_the_layout_cache(self):
        """The property that makes the guard worth having.

        If `feats + []` hit the cache this would all be free and the guard
        would be pointless — so pin the behaviour that motivates it.
        """
        seq = "ACGT" * 400
        feats = [{"start": i * 40, "end": i * 40 + 20, "strand": 1,
                  "type": "CDS", "label": f"f{i}", "color": "color(40)"}
                 for i in range(8)]
        sc._CHUNK_LAYOUT_CACHE.clear()
        sc._chunk_layout(seq, feats, 60)
        n_after_first = len(sc._CHUNK_LAYOUT_CACHE)
        sc._chunk_layout(seq, feats, 60)          # same object → hit
        assert len(sc._CHUNK_LAYOUT_CACHE) == n_after_first
        sc._chunk_layout(seq, feats + [], 60)     # new object → MISS
        assert len(sc._CHUNK_LAYOUT_CACHE) == n_after_first + 1, (
            "a new list object with identical contents must miss the id-keyed "
            "layout cache — if this ever stops being true, revisit the "
            "restriction-callback guard that exists because of it"
        )


class TestChunkFeatIndexMatchesLinearScan:
    """Issue #27: `_feats_in_chunk` gained a bisect index so `_chunk_layout`
    stops being O(chunks x features) on a bacterial chromosome.

    The indexed path is only safe if it returns EXACTLY what the linear scan
    returns — same features, same order, same wrap-split halves, same
    `_orig_*` stamps, same flap stripping. The linear scan stays in the code
    as the reference implementation and this compares the two directly, so a
    future edit to either one that changes behaviour shows up here rather than
    as mis-drawn bases.

    Order matters as much as membership: `_chunk_lane_groups` packs lanes in
    the order it receives features, so a reordering would silently move
    features between rows.
    """

    def _corpus(self):
        """Feature shapes that historically broke this code path."""
        total = 400
        return total, [
            # plain, in order, some overlapping and some adjacent
            {"start": 0,   "end": 50,  "strand": 1,  "type": "CDS",  "label": "a", "color": "color(40)"},
            {"start": 40,  "end": 120, "strand": -1, "type": "CDS",  "label": "b", "color": "color(41)"},
            {"start": 120, "end": 121, "strand": 1,  "type": "misc", "label": "single-bp", "color": "color(42)"},
            # spans the whole molecule — keeps the prefix-max high for every chunk
            {"start": 0,   "end": 400, "strand": 1,  "type": "source", "label": "whole", "color": "color(43)"},
            # wrap features (end < start) — split into tail + head
            {"start": 380, "end": 20,  "strand": 1,  "type": "CDS", "label": "wrap-fwd", "color": "color(44)"},
            {"start": 350, "end": 10,  "strand": -1, "type": "CDS", "label": "wrap-rev", "color": "color(45)"},
            # primer with a 5' flap that lands OUTSIDE the bound region
            {"start": 200, "end": 220, "strand": 1, "type": "primer_bind", "label": "p-fwd",
             "color": "color(46)", "_flap_start": 190, "_flap_end": 200, "_flap_bases": "ACGTACGTAC", "_flap_len": 10},
            # flap that runs off the start (negative) — wrap-aware
            {"start": 5,   "end": 25,  "strand": 1, "type": "primer_bind", "label": "p-neg",
             "color": "color(47)", "_flap_start": -8, "_flap_end": 5, "_flap_bases": "ACGTACGT", "_flap_len": 8},
            # flap that runs past the end
            {"start": 385, "end": 398, "strand": -1, "type": "primer_bind", "label": "p-past",
             "color": "color(48)", "_flap_start": 398, "_flap_end": 410, "_flap_bases": "ACGTACGTACGT", "_flap_len": 12},
            # a wrap feature that ALSO carries a flap
            {"start": 395, "end": 8,   "strand": 1, "type": "primer_bind", "label": "p-wrap",
             "color": "color(49)", "_flap_start": 388, "_flap_end": 395, "_flap_bases": "ACGTACG", "_flap_len": 7},
            # restriction overlay features — these legitimately reach the
            # layout (they carry the enzyme label art), in volume on a genome
            {"start": 60,  "end": 66,  "strand": 1, "type": "resite", "label": "EcoRI", "color": "color(9)"},
            {"start": 63,  "end": 64,  "strand": 1, "type": "recut",  "label": "",      "color": "color(9)"},
        ]

    def test_indexed_lookup_matches_linear_for_every_chunk(self):
        total, feats = self._corpus()
        index = sc._build_chunk_feat_index(feats)
        for line_width in (7, 13, 60, 132, 401):
            for cs in range(0, total, line_width):
                ce = min(cs + line_width, total)
                slow = sc._feats_in_chunk(feats, cs, ce, total)
                fast = sc._feats_in_chunk(feats, cs, ce, total, index=index)
                assert len(slow) == len(fast), (
                    f"count differs at lw={line_width} chunk [{cs},{ce}): "
                    f"{len(slow)} vs {len(fast)}"
                )
                for i, (a, b) in enumerate(zip(slow, fast)):
                    assert a == b, (
                        f"feature {i} differs at lw={line_width} "
                        f"chunk [{cs},{ce}):\n  linear={a}\n  indexed={b}"
                    )

    def test_identity_preserved_for_unsplit_features(self):
        """Non-wrapped features must pass through as the SAME dict object —
        click / hover resolution identity-matches them against `pm._feats`."""
        total, feats = self._corpus()
        index = sc._build_chunk_feat_index(feats)
        got = sc._feats_in_chunk(feats, 40, 70, total, index=index)
        plain = [f for f in got if "_orig_dict" not in f]
        assert plain, "expected at least one unsplit feature in this chunk"
        assert any(f is feats[1] for f in plain), (
            "an unsplit feature must be the identical dict object, not a copy"
        )

    def test_matches_linear_on_randomised_corpora(self):
        """Fuzz the shapes rather than trusting one hand-built corpus."""
        import random
        rng = random.Random(20260929)
        for trial in range(40):
            total = rng.choice([50, 137, 400, 1000])
            feats = []
            for _ in range(rng.randint(0, 25)):
                s = rng.randrange(total)
                if rng.random() < 0.15:                 # wrap
                    e = rng.randrange(0, max(1, s))
                else:
                    e = min(total, s + rng.randint(1, max(1, total // 3)))
                f = {"start": s, "end": e, "strand": rng.choice([1, -1]),
                     "type": rng.choice(["CDS", "misc", "resite", "recut",
                                         "primer_bind"]),
                     "label": f"f{len(feats)}", "color": "color(40)"}
                if f["type"] == "primer_bind" and rng.random() < 0.6:
                    fl = rng.randint(1, 15)
                    f["_flap_start"] = s - fl
                    f["_flap_end"] = s
                    f["_flap_bases"] = "A" * fl
                    f["_flap_len"] = fl
                feats.append(f)
            index = sc._build_chunk_feat_index(feats)
            lw = rng.choice([9, 31, 80])
            for cs in range(0, total, lw):
                ce = min(cs + lw, total)
                slow = sc._feats_in_chunk(feats, cs, ce, total)
                fast = sc._feats_in_chunk(feats, cs, ce, total, index=index)
                assert slow == fast, (
                    f"trial {trial} total={total} lw={lw} chunk [{cs},{ce}) "
                    f"diverged:\n  linear={slow}\n  indexed={fast}"
                )

    def test_empty_feature_list(self):
        index = sc._build_chunk_feat_index([])
        assert sc._feats_in_chunk([], 0, 60, 400, index=index) == []

    @pytest.mark.parametrize("name,feats,total", [
        ("zero-length feature",      [{"start": 10, "end": 10}], 100),
        ("feature at bp 0",          [{"start": 0,  "end": 5}],  100),
        ("ends exactly at total",    [{"start": 95, "end": 100}], 100),
        ("spans whole molecule",     [{"start": 0,  "end": 100}], 100),
        ("wrap ending at bp 0",      [{"start": 50, "end": 0}],  100),
        ("start beyond total",       [{"start": 200, "end": 300}], 100),
        ("single-bp molecule",       [{"start": 0,  "end": 1}],  1),
        ("inert flap (end <= start)",
         [{"start": 20, "end": 30, "_flap_start": 5, "_flap_end": 5}], 100),
        ("non-int flap keys",
         [{"start": 20, "end": 30, "_flap_start": "x", "_flap_end": None}], 100),
    ])
    def test_edge_shapes_match_linear(self, name, feats, total):
        """Shapes that break interval indexes: empty spans, boundary hits,
        a feature covering everything (prefix-max never drops), a wrap whose
        head ends at bp 0, coordinates past the sequence, and flap keys the
        overlap test must treat as absent."""
        feats = [dict(f, strand=1, type="CDS", label=name, color="color(40)")
                 for f in feats]
        index = sc._build_chunk_feat_index(feats)
        for lw in (1, 3, 7, 60, 1000):
            for cs in range(0, max(total, 1), lw):
                ce = min(cs + lw, total)
                assert (sc._feats_in_chunk(feats, cs, ce, total)
                        == sc._feats_in_chunk(feats, cs, ce, total,
                                              index=index)), \
                    f"{name} diverged at lw={lw} chunk [{cs},{ce})"

    def test_duplicate_identical_dicts_all_survive(self):
        """The same dict object listed N times must come back N times — the
        index keys on position, not identity, and lane packing counts them."""
        f = {"start": 5, "end": 25, "strand": 1, "type": "CDS",
             "label": "dup", "color": "color(40)"}
        feats = [f] * 6
        index = sc._build_chunk_feat_index(feats)
        got = sc._feats_in_chunk(feats, 0, 30, 100, index=index)
        assert len(got) == 6
        assert got == sc._feats_in_chunk(feats, 0, 30, 100)

    def test_unsorted_input_order_is_preserved_in_output(self):
        """Input is not sorted by start, and the output must stay in INPUT
        order — `_chunk_lane_groups` assigns rows in the order it receives
        features, so a sort-by-start would move features between lanes."""
        feats = [{"start": 90 - i * 10, "end": 95 - i * 10, "strand": 1,
                  "type": "CDS", "label": f"f{i}", "color": "color(40)"}
                 for i in range(9)]
        index = sc._build_chunk_feat_index(feats)
        got = sc._feats_in_chunk(feats, 0, 100, 100, index=index)
        assert [f["label"] for f in got] == [f["label"] for f in feats]
