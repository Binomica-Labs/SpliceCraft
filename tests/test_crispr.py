"""CRISPR guide design — `splicecraft_crispr` + its agent endpoints.

The guide finder is checked against an INDEPENDENT brute-force reference
written from the biology (spacer + PAM geometry), not from the implementation.
That is the only way a coordinate convention gets caught: reading the scanner
back to yourself agrees with whatever it already does.

The oligo pair is checked for ANNEALING (equal lengths, complementary duplex
core) and against the published pX330 `CACCG` / `AAAC…C` form, because a pair
that does not anneal is a wasted cloning week and the missing trailing C is the
classic way it is written down wrong.
"""

from __future__ import annotations

import random

import pytest

import splicecraft as sc
import splicecraft_crispr as cr


_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _brute_spcas9(seq: str, circular: bool = False) -> dict:
    """Every SpCas9 site from first principles.

    Plus strand reads  5'-[20 nt spacer][NGG]-3'.
    Minus strand is the same thing on the bottom strand, which on the FORWARD
    strand reads 5'-[CCN][revcomp of the spacer]-3'.
    """
    n = len(seq)
    out: dict = {}
    span = n if circular else n - 22

    def get(a: int, ln: int) -> str:
        if circular:
            return "".join(seq[(a + k) % n] for k in range(ln))
        return seq[a:a + ln] if 0 <= a and a + ln <= n else ""

    for i in range(max(0, span)):
        sp, pam = get(i, 20), get(i + 20, 3)
        if len(sp) == 20 and len(pam) == 3 and pam[1:] == "GG" \
                and all(c in "ACGT" for c in sp):
            out[(i % n, 1)] = (sp, pam)
    for i in range(max(0, span)):
        pam_f, sp_f = get(i, 3), get(i + 3, 20)
        if len(pam_f) == 3 and len(sp_f) == 20 and pam_f[:2] == "CC" \
                and all(c in "ACGT" for c in sp_f):
            out[((i + 3) % n, -1)] = (_rc(sp_f), _rc(pam_f))
    return out


def _seen(seq, **kw):
    return {(g["start"], g["strand"]): (g["guide"], g["pam"])
            for g in cr._find_guides(seq, **kw)}


class TestGuideFinderAgainstBruteForce:

    @pytest.mark.parametrize("n", [60, 120, 300, 900])
    def test_linear_matches_first_principles(self, n):
        rnd = random.Random(n)
        for _ in range(6):
            seq = "".join(rnd.choice("ACGT") for _ in range(n))
            assert _seen(seq, circular=False) == _brute_spcas9(seq, False)

    @pytest.mark.parametrize("n", [40, 80, 200])
    def test_circular_matches_first_principles(self, n):
        """Origin wraps are where a scanner usually diverges (sacred #6)."""
        rnd = random.Random(n + 1)
        for _ in range(6):
            seq = "".join(rnd.choice("ACGT") for _ in range(n))
            assert _seen(seq, circular=True) == _brute_spcas9(seq, True)

    def test_plus_strand_hand_checked(self):
        seq = "A" * 20 + "AGG" + "T" * 20
        g = next(x for x in cr._find_guides(seq) if x["strand"] == 1)
        assert (g["start"], g["end"], g["pam"]) == (0, 20, "AGG")
        assert g["cut_site"] == 17            # 3 bp 5' of the PAM

    def test_minus_strand_uses_forward_coordinates(self):
        """Sacred #2: a reverse-strand hit stores the FORWARD coordinate."""
        seq = "G" * 20 + "CCT" + "T" * 20
        g = next(x for x in cr._find_guides(seq) if x["strand"] == -1)
        assert g["start"] == 23               # forward frame, not 0
        assert g["guide"] == "A" * 20         # revcomp of the T-run
        assert g["pam"] == "AGG"             # oriented to the guide's strand
        assert g["cut_site"] == 26

    def test_overlapping_pams_are_all_found(self):
        """`NGG` sites overlap constantly; a non-overlapping finditer would
        silently drop every second one."""
        seq = "A" * 25 + "GGGGG" + "T" * 25
        starts = {g["start"] for g in cr._find_guides(seq) if g["strand"] == 1}
        assert len(starts) >= 3

    def test_ambiguous_protospacer_is_skipped(self):
        """A guide ordered against an N is wrong; skip rather than emit it."""
        seq = "A" * 10 + "N" + "A" * 9 + "AGG" + "T" * 20
        assert not [g for g in cr._find_guides(seq) if g["strand"] == 1]

    def test_linear_does_not_invent_a_wrapped_site(self):
        seq = "GG" + "A" * 40
        assert all(not g["wraps"] for g in cr._find_guides(seq, circular=False))

    def test_sequence_shorter_than_one_site(self):
        assert cr._find_guides("ACGT") == []


class TestCasVariants:

    def test_every_variant_scans(self):
        rnd = random.Random(3)
        seq = "".join(rnd.choice("ACGT") for _ in range(2000))
        for key in cr._cas_variant_names():
            hits = cr._find_guides(seq, variant=key)
            v = cr._CAS_VARIANTS[key]
            assert all(len(h["guide"]) == v["guide_len"] for h in hits)
            assert all(len(h["pam"]) == len(v["pam"]) for h in hits)

    def test_cpf1_pam_is_five_prime(self):
        """TTTV sits 5' of the spacer — the opposite geometry to Cas9, and the
        easiest thing to get backwards."""
        seq = "TTTA" + "ACGT" * 6 + "G" * 40
        g = next(x for x in cr._find_guides(seq, variant="lbcas12a")
                 if x["strand"] == 1)
        assert g["pam"] == "TTTA"
        assert g["pam_end"] == g["start"]       # PAM immediately 5' of spacer

    def test_unknown_variant_names_the_alternatives(self):
        with pytest.raises(ValueError, match="spcas9"):
            cr._cas_variant("cas13")

    def test_variant_lookup_is_forgiving_about_formatting(self):
        assert cr._cas_variant("SpCas9")["pam"] == "NGG"
        assert cr._cas_variant(" sp-cas9 ")["pam"] == "NGG"


class TestGuideTriage:

    def test_poliii_terminator_is_fatal(self):
        """U6/H1 stops on a T-run, so the guide is truncated before it is ever
        loaded. A mechanism, not a correlation."""
        s = cr._score_guide("GGGGTTTTGGGGCCCCAAAA")
        assert "Pol III terminator (TTTT)" in s["flags"]
        assert s["tier"] == "poor"

    def test_clean_guide_is_good(self):
        s = cr._score_guide("GACCTGAAGTCTAGGTCCCT")
        assert s["flags"] == [] and s["tier"] == "good"

    def test_missing_five_prime_g_alone_does_not_demote(self):
        """It is fixed by adding one base, so it must not outrank a real
        concern."""
        s = cr._score_guide("ACCTGAAGTCTAGGTCCCTA")
        assert s["flags"] == ["does not start with G"]
        assert s["tier"] == "good"

    def test_u6_driven_false_drops_the_g_advice(self):
        s = cr._score_guide("ACCTGAAGTCTAGGTCCCTA", u6_driven=False)
        assert s["flags"] == []

    @pytest.mark.parametrize("guide,needle", [
        ("AAAAAAAAAAAAAAAAAAAA", "homopolymer"),
        ("GCGCGCGCGCGCGCGCGCGC", "GC"),
        ("ATATATATATATATATATAT", "GC"),
    ])
    def test_named_concerns(self, guide, needle):
        assert any(needle in f for f in cr._score_guide(guide)["flags"])

    def test_empty_guide_does_not_crash(self):
        assert cr._score_guide("")["tier"] == "poor"

    def test_no_fabricated_efficiency_number(self):
        """Guard the design decision: a number that LOOKS like a trained
        on-target prediction and is not would be worse than no number."""
        s = cr._score_guide("GACCTGAAGTCTAGGTCCCT")
        for banned in ("efficiency", "score", "doench", "rule_set",
                       "on_target"):
            assert banned not in s, f"{banned!r} implies a model we do not have"


class TestCloningOligos:

    @pytest.mark.parametrize("guide", [
        "ACCTGAAGTCTAGGTCCCTA",     # needs a G added
        "GCCTGAAGTCTAGGTCCCTA",     # already starts with G
    ])
    def test_the_two_oligos_actually_anneal(self, guide):
        o = cr._guide_cloning_oligos(guide)
        assert len(o["top"]) == len(o["bottom"]), "a gap at one end"
        core_top = o["top"][len(o["overhang_top"]):]
        core_bot = o["bottom"][len(o["overhang_bottom"]):]
        assert core_top == _rc(core_bot)

    def test_matches_the_published_px330_form(self):
        """Regression guard for 2026-09-17: the first cut emitted
        `AAAC + rc(guide)` and dropped the trailing C that pairs with the G in
        CACCG, so the two oligos came out 25 and 24 nt."""
        o = cr._guide_cloning_oligos("ACCTGAAGTCTAGGTCCCTA")
        assert o["top"] == "CACCG" + "ACCTGAAGTCTAGGTCCCTA"
        assert o["bottom"].startswith("AAAC") and o["bottom"].endswith("C")
        assert len(o["top"]) == len(o["bottom"]) == 25
        assert o["added_g"] is True

    def test_a_g_initiated_guide_gets_no_second_g(self):
        o = cr._guide_cloning_oligos("GCCTGAAGTCTAGGTCCCTA")
        assert o["added_g"] is False
        assert o["duplex"] == "GCCTGAAGTCTAGGTCCCTA"
        assert not o["top"].startswith("CACCGG")

    def test_add_g_off_is_honoured(self):
        o = cr._guide_cloning_oligos("ACCTGAAGTCTAGGTCCCTA", add_g=False)
        assert o["added_g"] is False
        assert o["duplex"] == "ACCTGAAGTCTAGGTCCCTA"

    def test_custom_overhangs(self):
        o = cr._guide_cloning_oligos("GACCTGAAGTCTAGGTCCCT",
                                     overhang_top="AGAT",
                                     overhang_bottom="TTTG")
        assert o["top"].startswith("AGAT") and o["bottom"].startswith("TTTG")
        assert len(o["top"]) == len(o["bottom"])

    @pytest.mark.parametrize("bad", ["", "ACGTN", "not dna", None])
    def test_rejects_an_unusable_spacer(self, bad):
        with pytest.raises(ValueError):
            cr._guide_cloning_oligos(bad)  # type: ignore[arg-type]

    def test_rejects_ambiguous_overhang(self):
        with pytest.raises(ValueError, match="overhang_top"):
            cr._guide_cloning_oligos("GACCTGAAGTCTAGGTCCCT",
                                     overhang_top="CANN")


class TestOffTargets:

    def test_finds_its_own_site_with_zero_mismatches(self):
        rnd = random.Random(5)
        seq = "".join(rnd.choice("ACGT") for _ in range(400))
        g = cr._find_guides(seq)[0]
        res = cr._guide_offtargets(g["guide"], seq)
        assert any(h["mismatches"] == 0 for h in res["hits"])

    def test_reports_the_scope_of_the_search(self):
        """`searched_bp` is what stops "no off-targets" being read as a
        genome-wide claim."""
        res = cr._guide_offtargets("ACCTGAAGTCTAGGTCCCTA", "ACGT" * 100)
        assert res["searched_bp"] == 400

    def test_seed_mismatches_are_counted_separately(self):
        """Cas9 tolerates distal mismatches far more than PAM-proximal ones, so
        3 distal and 3 seed mismatches are not the same risk."""
        spacer = "ACCTGAAGTCTAGGTCCCTA"
        distal = "TTT" + spacer[3:]            # 3 mismatches, PAM-distal
        target = distal + "AGG"
        res = cr._guide_offtargets(spacer, target, max_mismatch=4)
        hit = next(h for h in res["hits"] if h["mismatches"] == 3)
        assert hit["seed_mismatches"] == 0

    def test_hits_are_ordered_most_dangerous_first(self):
        rnd = random.Random(9)
        seq = "".join(rnd.choice("ACGT") for _ in range(600))
        g = cr._find_guides(seq)[0]
        res = cr._guide_offtargets(g["guide"], seq, max_mismatch=4)
        keys = [(h["seed_mismatches"], h["mismatches"]) for h in res["hits"]]
        assert keys == sorted(keys)

    def test_require_pam_filters(self):
        spacer = "ACCTGAAGTCTAGGTCCCTA"
        no_pam = spacer + "TTT"                # no NGG
        assert cr._guide_offtargets(spacer, no_pam,
                                    require_pam=True)["hits"] == []
        assert cr._guide_offtargets(spacer, no_pam,
                                    require_pam=False)["hits"]

    def test_empty_inputs(self):
        assert cr._guide_offtargets("", "ACGT")["hits"] == []
        assert cr._guide_offtargets("ACGTACGTACGTACGTACGT", "")["hits"] == []


class TestDesignGuides:

    def _seq(self, n=1200, seed=4):
        rnd = random.Random(seed)
        return "".join(rnd.choice("ACGT") for _ in range(n))

    def test_region_restricts_by_cut_site(self):
        seq = self._seq()
        res = cr._design_guides(seq, region=(400, 600), limit=500)
        assert res["guides"]
        assert all(400 <= g["cut_site"] < 600 for g in res["guides"])

    def test_no_offtarget_request_makes_no_claim(self):
        """`offtarget_searched_bp == 0` so silence cannot read as 'clean'."""
        res = cr._design_guides(self._seq())
        assert res["offtarget_searched_bp"] == 0
        assert all("offtargets" not in g for g in res["guides"])

    def test_offtarget_drops_the_on_target_hit(self):
        seq = self._seq()
        res = cr._design_guides(seq, offtarget_in=seq, limit=5)
        assert res["offtarget_searched_bp"] == len(seq)
        for g in res["guides"]:
            assert all(h["mismatches"] > 0 for h in g["offtargets"])

    def test_good_guides_sort_first(self):
        res = cr._design_guides(self._seq(), limit=500)
        rank = {"good": 0, "ok": 1, "poor": 2}
        tiers = [rank[g["score"]["tier"]] for g in res["guides"]]
        assert tiers == sorted(tiers)

    def test_truncated_is_reported(self):
        res = cr._design_guides(self._seq(4000), limit=5)
        assert res["truncated"] is True and len(res["guides"]) == 5

    def test_real_plasmid_shape(self):
        """pUC19-like: circular, and every guide's coordinates stay in range."""
        seq = self._seq(2686)
        res = cr._design_guides(seq, circular=True, limit=500)
        assert res["n_found"] > 100
        for g in res["guides"]:
            assert 0 <= g["start"] < len(seq)
            assert 0 <= g["cut_site"] < len(seq)


class TestCrisprAgentEndpoints:
    """Every user-callable operation has an agent side-door."""

    class _App:
        _current_record = None

    def _call(self, name, payload, rec=None):
        import splicecraft_state as st
        app = self._App()
        app._current_record = rec
        fn, _write = st._AGENT_HANDLERS[name]
        return fn(app, payload)

    def _rec(self, n=1200):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rnd = random.Random(2)
        return SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(n))),
                         id="p", name="p",
                         annotations={"molecule_type": "DNA",
                                      "topology": "circular"})

    @pytest.mark.parametrize("name", [
        "list-cas-variants", "design-guides", "score-guide",
        "guide-offtargets", "guide-cloning-oligos",
    ])
    def test_endpoint_is_registered(self, name):
        import splicecraft_state as st
        assert name in st._AGENT_HANDLERS

    def test_list_cas_variants(self):
        out = self._call("list-cas-variants", {})
        keys = {v["key"] for v in out["variants"]}
        assert {"spcas9", "sacas9", "lbcas12a"} <= keys

    def test_design_guides_on_the_loaded_plasmid(self):
        out = self._call("design-guides", {"limit": 5}, rec=self._rec())
        assert out["n_found"] > 0 and len(out["guides"]) == 5
        assert out["circular"] is True

    def test_design_guides_accepts_a_raw_sequence(self):
        rnd = random.Random(8)
        seq = "".join(rnd.choice("ACGT") for _ in range(400))
        out = self._call("design-guides", {"sequence": seq})
        assert out["n_found"] > 0

    def test_raw_sequence_defaults_to_linear(self):
        """A pasted exon is the common case; wrapping it would invent guides
        across a junction that does not exist."""
        out = self._call("design-guides", {"sequence": "ACGT" * 100})
        assert out["circular"] is False

    def test_no_plasmid_and_no_sequence_is_422(self):
        body, code = self._call("design-guides", {})
        assert code == 422 and "no plasmid" in body["error"]

    @pytest.mark.parametrize("payload,needle", [
        ({"variant": "cas13"}, "unknown Cas variant"),
        ({"region": [5, 5]}, "empty"),
        ({"region": [0]}, "[start, end]"),
        ({"sequence": "ATGXYZ"}, "non-IUPAC"),
        ({"max_mismatch": 99}, "0-10"),
        ({"limit": 0}, ">= 1"),
    ])
    def test_bad_payloads_are_400(self, payload, needle):
        body, code = self._call("design-guides", payload, rec=self._rec())
        assert code == 400 and needle in body["error"]

    def test_offtarget_true_searches_the_target(self):
        rec = self._rec(600)
        out = self._call("design-guides",
                         {"offtarget": True, "limit": 3}, rec=rec)
        assert out["offtarget_searched_bp"] == 600

    def test_offtarget_accepts_another_sequence(self):
        out = self._call("design-guides",
                         {"offtarget": "ACGT" * 50, "limit": 2},
                         rec=self._rec(400))
        assert out["offtarget_searched_bp"] == 200

    def test_score_guide(self):
        out = self._call("score-guide", {"guide": "GGGGTTTTGGGGCCCCAAAA"})
        assert out["tier"] == "poor" and out["guide"]

    def test_score_guide_requires_a_guide(self):
        body, code = self._call("score-guide", {})
        assert code == 400

    def test_guide_cloning_oligos_endpoint(self):
        out = self._call("guide-cloning-oligos",
                         {"guide": "ACCTGAAGTCTAGGTCCCTA"})
        assert len(out["top"]) == len(out["bottom"]) == 25

    def test_guide_cloning_oligos_custom_overhangs(self):
        out = self._call("guide-cloning-oligos",
                         {"guide": "GACCTGAAGTCTAGGTCCCT",
                          "overhang_top": "AGAT", "overhang_bottom": "TTTG"})
        assert out["overhang_top"] == "AGAT"

    def test_guide_offtargets_endpoint(self):
        rec = self._rec(500)
        g = cr._find_guides(str(rec.seq), circular=True)[0]
        out = self._call("guide-offtargets",
                         {"guide": g["guide"]}, rec=rec)
        assert out["searched_bp"] == 500
        assert any(h["mismatches"] == 0 for h in out["hits"])


class TestHubSurface:

    @pytest.mark.parametrize("name", [
        "_find_guides", "_design_guides", "_score_guide",
        "_guide_offtargets", "_guide_cloning_oligos", "_cas_variant",
        "_CAS_VARIANTS", "_cas_variant_names",
    ])
    def test_reexported_on_the_hub(self, name):
        assert hasattr(sc, name), f"sc.{name} must resolve"

    def test_take_circular_does_not_shadow_biology_slice(self):
        """Two different contracts (length vs end boundary) deliberately keep
        two different names — sharing one made every call site ambiguous."""
        assert sc._slice_circular("ABCDEFGH", 6, 2) == "GHAB"      # end
        assert cr._take_circular("ABCDEFGH", 6, 4) == "GHAB"       # length
        assert not hasattr(sc, "_take_circular")
