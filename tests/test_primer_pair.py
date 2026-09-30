"""
Primer-pair QC, Tm profiles, and thread-safe primer3 access (PR #28, [INV-211]).

Covers:

  * `_p3_tm` / `_p3_thermo`         — the ONE way SpliceCraft calls primer3:
                                       parity with the module functions, a
                                       per-call condition reset, unknown
                                       keywords refused, and — the reason they
                                       exist — correct answers under threads
  * Tm profiles                      — registry, aliases, primer3's defaults
                                       pinned, Benchling reference Tms
  * `_primer_oligo_qc` / `_primer_pair_qc` — structure Tm vs Primer3's 47 °C
                                       limit, the 60-nt thal limit, degenerate
                                       oligos, the 3'-anchored dedupe
  * `_insilico_pcr_amplicons`        — ranks every pair BEFORE it caps,
                                       checked against a brute-force reference
  * `check-primer-pair`, `list-tm-profiles`, the `tm_profile` + `qc` parts of
    `design-primers` / `check-primer`
  * the Primer Design screen         — the pair check under a design result and
                                       in the Primer Check tab

Everything numeric is checked against primer3's own module functions or
against Primer3's DESIGN output, never against SpliceCraft's own numbers.
"""
from __future__ import annotations

import random
import threading
import types

import pytest

import splicecraft as sc

primer3 = pytest.importorskip("primer3")

TERMINAL_SIZE = (160, 50)
_P3 = dict(mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0)
_BENCH = dict(mv_conc=50.0, dv_conc=0.0, dntp_conc=0.0, dna_conc=250.0)


def _rand_dna(n: int, seed: int) -> str:
    return "".join(random.Random(seed).choices("ACGT", k=n))


def _rc(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


# Both fold into the same ~61 °C hairpin (a harmless-looking −2.6 to −3.2
# kcal/mol), and each is the other's reverse complement, so together they also
# make a strong primer-dimer. The difference is WHERE: the primer's hairpin is
# at its 5' end and leaves the 3' end free; the partner's hairpin holds its 3'
# end — the one that can't prime.
_HAIRPIN_PRIMER = "ACGTAGGCCTACGTTTTTTTTTTTT"
_HAIRPIN_PARTNER = "AAAAAAAAAAAACGTAGGCCTACGT"
_HAIRPIN_3P = "TTTTTTTTTTTACGTAGGCCTACGT"
# An ordinary, well-behaved pair (M13/pUC-style sequencing primers).
_CLEAN_FWD = "AGCGGATAACAATTTCACACAGGA"
_CLEAN_REV = "GTAAAACGACGGCCAGTGAATTGT"


# ═══════════════════════════════════════════════════════════════════════════════
# primer3 access: parity, reset, keywords, threads
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrimer3Access:
    SEQS = ["TGGCCTTTTTGCGTTTCTACAAACT", _HAIRPIN_PRIMER, _CLEAN_REV,
            "GCGCGAATTCATGACCATGATTACGCCAAGCTTGCATGCCTGCAGG"]

    @pytest.mark.parametrize("cond", [{}, _P3, _BENCH])
    def test_tm_matches_the_module_function(self, cond):
        for s in self.SEQS:
            assert sc._p3_tm(s, **cond) == primer3.calc_tm(s, **cond)

    def test_every_alignment_kind_matches_the_module_function(self):
        a, b = _HAIRPIN_PRIMER, _HAIRPIN_PARTNER
        pairs = (("hairpin", primer3.calc_hairpin, (a,)),
                 ("homodimer", primer3.calc_homodimer, (a,)),
                 ("heterodimer", primer3.calc_heterodimer, (a, b)),
                 ("end_stability", primer3.calc_end_stability, (a, b)))
        for kind, fn, args in pairs:
            ours = sc._p3_thermo(kind, *args, **_P3)
            ref = fn(*args, **_P3)
            assert (ours.dg, ours.tm) == (ref.dg, ref.tm), kind

    def test_conditions_do_not_leak_into_the_next_call(self):
        # The per-thread instance is RE-set every call: asking for Benchling's
        # conditions once must not leave them behind for a default call.
        s = self.SEQS[0]
        bench = sc._p3_tm(s, **_BENCH)
        assert sc._p3_tm(s) == primer3.calc_tm(s) != bench
        sc._p3_thermo("hairpin", _HAIRPIN_PRIMER, dv_conc=0.0)
        assert (sc._p3_thermo("hairpin", _HAIRPIN_PRIMER).dg
                == primer3.calc_hairpin(_HAIRPIN_PRIMER).dg)

    def test_unknown_keywords_are_refused_not_swallowed(self):
        # ThermoAnalysis.set_thermo_args takes **kwargs and ignores a typo —
        # which would quietly compute under the default for that condition.
        with pytest.raises(TypeError):
            sc._p3_tm("ACGTACGTACGTACGTACGT", mg_conc=0.0)
        with pytest.raises(TypeError):
            sc._p3_thermo("hairpin", "ACGTACGTACGTACGTACGT", tm_method="x")
        with pytest.raises(ValueError):
            sc._p3_thermo("pseudoknot", "ACGTACGTACGTACGTACGT")

    def test_long_pair_is_refused_like_primer3(self):
        with pytest.raises(RuntimeError):
            sc._p3_thermo("heterodimer", "A" * 61, "T" * 61)

    def test_concurrent_structures_are_the_serial_answers(self):
        """The reason `_p3_thermo` exists. primer3's alignment keeps its DP
        tables in C statics with the GIL released, so concurrent hairpin and
        homodimer calls corrupted each other (~20 % wrong on 60-mers across
        six threads). Distinct sequences per call, so no cache can hide it."""
        rnd = random.Random(11)
        per_thread, n_threads = 30, 6
        work = [[("".join(rnd.choice("ACGT") for _ in range(58)),
                  rnd.random() < 0.5) for _ in range(per_thread)]
                for _ in range(n_threads)]
        ref = {(s, hp): (primer3.calc_hairpin(s, **_P3).dg if hp
                         else primer3.calc_homodimer(s, **_P3).dg) / 1000.0
               for lot in work for s, hp in lot}
        import splicecraft_primer as sp
        sp._mut_thermo_cache_clear()
        wrong: "list[tuple]" = []
        start = threading.Barrier(n_threads)

        def run(lot):
            start.wait()
            for s, hp in lot:
                got = sc._mut_hairpin_dg(s) if hp else sc._mut_homodimer_dg(s)
                if got is None or abs(got - ref[(s, hp)]) > 1e-9:
                    wrong.append((s, hp, got, ref[(s, hp)]))

        threads = [threading.Thread(target=run, args=(lot,)) for lot in work]
        try:
            for t in threads:
                t.start()
            for t in threads:
                t.join(timeout=120)
        finally:
            sp._mut_thermo_cache_clear()
        assert not wrong, f"{len(wrong)} of {per_thread * n_threads} wrong"

    def test_concurrent_tms_under_different_profiles(self):
        # Two conditions at once on the SAME sequences: each thread must get
        # its own profile's number, never the other thread's.
        seqs = [_rand_dna(25, seed) for seed in range(40)]
        want = {name: {s: sc._primer_tm(s, name) for s in seqs}
                for name in sc._TM_PROFILES}
        wrong: list = []
        start = threading.Barrier(4)

        def run(name):
            start.wait()
            for _ in range(15):
                for s in seqs:
                    if sc._primer_tm(s, name) != want[name][s]:
                        wrong.append((name, s))

        names = list(sc._TM_PROFILES) * 2
        threads = [threading.Thread(target=run, args=(n,)) for n in names]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=120)
        assert not wrong


# ═══════════════════════════════════════════════════════════════════════════════
# Tm profiles
# ═══════════════════════════════════════════════════════════════════════════════

class TestTmProfiles:
    def test_default_profile_is_the_app_wide_conditions(self):
        assert sc._TM_PROFILE_DEFAULT == "primer3_default"
        assert sc._tm_profile_kwargs("primer3_default") == sc._MUT_P3

    def test_primer3_defaults_are_pinned(self):
        """The default profile, `_MUT_P3`, and every bare `calc_tm(s)` in the
        app must be ONE set of conditions. A primer3 release that moved its
        defaults would silently shift every default-profile Tm — this fails
        first."""
        for seed in range(30):
            s = _rand_dna(18 + seed, seed)
            assert primer3.calc_tm(s) == primer3.calc_tm(s, **sc._MUT_P3)
            assert sc._primer_tm(s) == round(primer3.calc_tm(s), 1)
            assert (primer3.calc_hairpin(s).dg
                    == primer3.calc_hairpin(s, **sc._MUT_P3).dg)

    @pytest.mark.parametrize(("sequence", "expected"), [
        # Read off Benchling (default SantaLucia settings) by the PR author.
        ("TACTGTTTCTCCATACCCGTTTTTTTGGG", 59.4),
        ("TGGCCTTTTTGCGTTTCTACAAACT", 58.2),
        ("CTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCC"
         "CGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTAT"
         "CGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGG", 80.3),
    ])
    def test_benchling_reference_tms(self, sequence, expected):
        assert sc._primer_tm(sequence, "benchling_compatible") == expected
        # …and the default profile is a different, higher number.
        assert sc._primer_tm(sequence) > expected + 3

    @pytest.mark.parametrize(("raw", "key"), [
        (None, "primer3_default"), ("", "primer3_default"),
        ("  ", "primer3_default"), ("primer3_default", "primer3_default"),
        ("PRIMER3_DEFAULT", "primer3_default"), ("default", "primer3_default"),
        ("primer3", "primer3_default"), ("benchling", "benchling_compatible"),
        ("Benchling-Compatible", "benchling_compatible"),
        (" benchling compatible ", "benchling_compatible"),
    ])
    def test_spellings_resolve(self, raw, key):
        assert sc._tm_profile_resolve(raw) == key

    @pytest.mark.parametrize("raw", ["snappy", "benchling_default", 1, True,
                                     ["benchling"], {"k": "v"}])
    def test_unknown_is_none_never_the_default(self, raw):
        assert sc._tm_profile_resolve(raw) is None
        with pytest.raises(ValueError):
            sc._tm_profile_kwargs(raw)

    def test_kwargs_are_a_copy(self):
        kw = sc._tm_profile_kwargs("benchling_compatible")
        kw["dv_conc"] = 99.0
        assert sc._tm_profile_kwargs("benchling_compatible")["dv_conc"] == 0.0

    def test_catalog_spells_out_the_conditions(self):
        cat = {p["key"]: p for p in sc._tm_profile_catalog()}
        assert set(cat) == set(sc._TM_PROFILES)
        assert [k for k, p in cat.items() if p["default"]] == ["primer3_default"]
        b = cat["benchling_compatible"]
        assert (b["monovalent_mm"], b["mg_mm"], b["dntp_mm"], b["oligo_nm"]) \
            == (50.0, 0.0, 0.0, 250.0)
        d = cat["primer3_default"]
        assert (d["monovalent_mm"], d["mg_mm"], d["dntp_mm"], d["oligo_nm"]) \
            == (50.0, 1.5, 0.6, 50.0)

    def test_list_endpoint(self):
        r = sc._h_list_tm_profiles(None, {})
        assert r["ok"] and r["default"] == "primer3_default"
        assert r["profiles"] == sc._tm_profile_catalog()
        assert "47" in r["note"]

    def test_binding_region_pick_follows_the_profile(self):
        # Benchling's scale reads ~4 °C lower, so the same 60 °C target needs
        # a LONGER arm — and the returned Tm is on the profile's scale.
        seq = _rand_dna(60, 5)
        d_seq, d_tm = sc._pick_binding_region(seq, 60.0, max_len=40)
        b_seq, b_tm = sc._pick_binding_region(seq, 60.0, max_len=40,
                                              tm_profile="benchling_compatible")
        assert len(b_seq) > len(d_seq)
        assert b_tm == pytest.approx(primer3.calc_tm(b_seq, **_BENCH))
        assert d_tm == pytest.approx(primer3.calc_tm(d_seq))


# ═══════════════════════════════════════════════════════════════════════════════
# The QC engine
# ═══════════════════════════════════════════════════════════════════════════════

class TestPrimerStructure:
    def test_matches_primer3_design_output(self):
        """The structure Tm IS Primer3's own design metric — the number its
        PRIMER_MAX_*_TH limits test — floored at 0 the way it reports it."""
        tmpl = _rand_dna(900, 7)
        res = primer3.design_primers(
            {"SEQUENCE_TEMPLATE": tmpl,
             "SEQUENCE_INCLUDED_REGION": [0, len(tmpl)]},
            {"PRIMER_TASK": "generic", "PRIMER_PICK_LEFT_PRIMER": 1,
             "PRIMER_PICK_RIGHT_PRIMER": 1, "PRIMER_OPT_SIZE": 22,
             "PRIMER_PRODUCT_SIZE_RANGE": [[400, 600]],
             "PRIMER_NUM_RETURN": 1})
        left = res["PRIMER_LEFT_0_SEQUENCE"]
        right = res["PRIMER_RIGHT_0_SEQUENCE"]
        want = {
            ("hairpin", left, ""): res["PRIMER_LEFT_0_HAIRPIN_TH"],
            ("homodimer", left, ""): res["PRIMER_LEFT_0_SELF_ANY_TH"],
            ("end_stability", left, left): res["PRIMER_LEFT_0_SELF_END_TH"],
            ("heterodimer", left, right): res["PRIMER_PAIR_0_COMPL_ANY_TH"],
        }
        for (kind, a, b), tm in want.items():
            got = sc._primer_structure(kind, a, b)
            assert got is not None
            assert got["tm"] == round(max(0.0, tm), 1), kind

    def test_no_structure_reads_zero(self):
        st = sc._primer_structure("hairpin", "AAAAAAAAAAAAAAAAAAAA")
        assert st == {"dg": 0.0, "tm": 0.0, "end_paired": False}

    def test_which_end_the_hairpin_holds(self):
        # The stem drawn by primer3 decides it, checked against the drawing.
        for seq, want in ((_HAIRPIN_PRIMER, False), (_HAIRPIN_3P, True)):
            res = primer3.calc_hairpin(seq, output_structure=True, **_P3)
            drawn = res.ascii_structure_lines[0].split("\t", 1)[1]
            assert len(drawn) == len(seq)
            assert any(c in "/\\" for c in drawn[-5:]) is want
            assert sc._primer_structure("hairpin", seq)["end_paired"] is want

    def test_unmeasured_is_none(self, monkeypatch):
        import splicecraft_primer as sp
        monkeypatch.setattr(sp, "_MUT_P3", {"bogus_kwarg": 1})
        assert sc._primer_structure("hairpin", _CLEAN_FWD) is None


class TestPrimerOligoQc:
    def test_clean_primer(self):
        q = sc._primer_oligo_qc(_CLEAN_FWD)
        assert q["warnings"] == []
        assert q["tm"] == sc._primer_tm(_CLEAN_FWD)
        assert q["binding"] == q["primer"] == _CLEAN_FWD and q["tail"] == ""
        for f in ("hairpin", "homodimer", "homodimer_3p"):
            assert isinstance(q[f"{f}_dg"], float)
            assert isinstance(q[f"{f}_tm"], float)

    def test_a_hairpin_holding_the_3p_end_is_flagged_by_its_tm(self):
        """ΔG −3.2 kcal/mol is nowhere near the −9 rule of thumb, but the
        hairpin melts at ~60 °C — still folded at the anneal — and its stem
        holds the 3' end, so the primer cannot prime."""
        q = sc._primer_oligo_qc(_HAIRPIN_3P, role="forward")
        assert q["hairpin_dg"] > -9.0
        assert q["hairpin_tm"] >= q["anneal_ref_c"]
        assert q["hairpin_3p"] is True
        assert any(w.startswith("the forward primer's 3' end is held in a "
                                "hairpin") for w in q["warnings"])

    def test_a_hairpin_leaving_the_3p_end_free_is_reported_not_flagged(self):
        """Measured on 200 designs: a hairpin anywhere else — above all a
        cloning tail folding onto the arm — flagged 82 % of cloning pairs under
        a flat 47 °C limit, in primers labs use every day. It is shown, with
        `hairpin_3p` false, and not warned about."""
        q = sc._primer_oligo_qc(_HAIRPIN_PRIMER, role="forward")
        assert q["hairpin_tm"] > 47.0 and q["hairpin_3p"] is False
        assert not any("hairpin" in w for w in q["warnings"])

    def test_tail_sets_structures_arm_sets_tm(self):
        arm = "TTCTTGTTGAATTAGATGGTGA"
        full = "CGGTCACCTGCCCTGGGAG" + arm
        q = sc._primer_oligo_qc(full, binding=arm)
        assert q["tail"] == "CGGTCACCTGCCCTGGGAG"
        assert q["tm"] == sc._primer_tm(arm)
        assert q["gc_pct"] == round(100 * sum(c in "GC" for c in arm)
                                    / len(arm), 1)
        want = primer3.calc_hairpin(full, **_P3)
        assert q["hairpin_dg"] == round(want.dg / 1000.0, 2)

    def test_binding_must_be_the_three_prime_end(self):
        with pytest.raises(ValueError):
            sc._primer_oligo_qc("GGGGACGTACGTACGTACGT", binding="GGGGACGTAC")
        with pytest.raises(ValueError):
            sc._primer_oligo_qc("")

    def test_over_sixty_nt_is_not_folded_and_says_why(self):
        q = sc._primer_oligo_qc("GCGC" * 16 + _CLEAN_FWD, role="forward")
        assert q["hairpin_dg"] is None and q["homodimer_tm"] is None
        assert q["homodimer_3p_dg"] is None
        assert any("primer3 folds oligos only up to 60 nt" in w
                   for w in q["warnings"])
        # The Tm still comes back — a long oligo is not an unmeasured one.
        assert isinstance(q["tm"], float)

    def test_degenerate_oligo(self):
        q = sc._primer_oligo_qc("ACGTRYACGTACGTACGTNN", role="forward")
        lo, hi = q["tm_range"]
        assert q["tm"] == lo < hi
        assert any("ambiguity codes" in w for w in q["warnings"])

    def test_profile_moves_tm_not_structure(self):
        d = sc._primer_oligo_qc(_HAIRPIN_PRIMER)
        b = sc._primer_oligo_qc(_HAIRPIN_PRIMER,
                                tm_profile="benchling_compatible")
        assert b["tm"] == sc._primer_tm(_HAIRPIN_PRIMER,
                                        "benchling_compatible") != d["tm"]
        for f in ("hairpin_dg", "hairpin_tm", "homodimer_dg", "homodimer_tm",
                  "homodimer_3p_dg", "homodimer_3p_tm"):
            assert b[f] == d[f], f

    def test_one_structure_one_warning(self):
        # The reverse primer's strongest self-dimer IS 3'-anchored: flagged
        # once, as the 3' (extendable) one.
        q = sc._primer_oligo_qc(_HAIRPIN_PARTNER, role="reverse")
        assert q["homodimer_dg"] == q["homodimer_3p_dg"]
        dimers = [w for w in q["warnings"] if "self-dimer" in w]
        assert len(dimers) == 1 and "3' end" in dimers[0]

    def test_unmeasured_is_warned_not_passed(self, monkeypatch):
        import splicecraft_primer as sp
        monkeypatch.setattr(sp, "_MUT_P3", {"bogus_kwarg": 1})
        q = sc._primer_oligo_qc(_CLEAN_FWD, role="forward")
        assert q["hairpin_dg"] is None
        assert any("NOT secondary-structure checked" in w
                   for w in q["warnings"])


class TestPrimerPairQc:
    def test_clean_pair(self):
        q = sc._primer_pair_qc(_CLEAN_FWD, _CLEAN_REV)
        assert q["warnings"] == []
        assert q["tm_profile"] == "primer3_default"
        assert q["dg_units"] == "kcal/mol"
        assert q["dimer_tm_limit"] == 47.0
        assert q["anneal_ref_c"] == round(
            min(sc._primer_tm(_CLEAN_FWD), sc._primer_tm(_CLEAN_REV)) - 5, 1)
        assert q["pair"]["tm_delta"] == round(
            abs(q["forward"]["tm"] - q["reverse"]["tm"]), 1)
        het = primer3.calc_heterodimer(_CLEAN_FWD, _CLEAN_REV, **_P3)
        assert q["pair"]["heterodimer_dg"] == round(het.dg / 1000.0, 2)
        ends = [primer3.calc_end_stability(a, b, **_P3) for a, b in
                ((_CLEAN_FWD, _CLEAN_REV), (_CLEAN_REV, _CLEAN_FWD))]
        worst = max(ends, key=lambda e: e.tm)
        assert q["pair"]["heterodimer_3p_tm"] == round(max(0.0, worst.tm), 1)
        assert "warnings" not in q["forward"] and "warnings" not in q["reverse"]

    def test_primer_dimer_flagged(self):
        q = sc._primer_pair_qc(_HAIRPIN_PRIMER, _HAIRPIN_PARTNER)
        assert q["pair"]["heterodimer_tm"] > 47.0
        assert any(w.startswith("the two primers pair") for w in q["warnings"])

    def test_tm_mismatch_uses_the_amplify_feature_words(self):
        q = sc._primer_pair_qc("ATATATTAATTAAATATAAT", _CLEAN_REV)
        assert q["pair"]["tm_delta"] > 5.0
        assert f"primer Tms differ by {q['pair']['tm_delta']:.1f} °C" \
            in q["warnings"]

    def test_both_long_skips_the_dimer_one_long_does_not(self):
        long_f = "GCGC" * 16 + _CLEAN_FWD
        long_r = "ATAT" * 16 + _CLEAN_REV
        both = sc._primer_pair_qc(long_f, long_r)
        assert both["pair"]["heterodimer_dg"] is None
        assert any("both primers are longer than 60 nt" in w
                   for w in both["warnings"])
        one = sc._primer_pair_qc(long_f, _CLEAN_REV)
        assert isinstance(one["pair"]["heterodimer_dg"], float)
        assert isinstance(one["pair"]["heterodimer_3p_dg"], float)

    def test_unknown_profile(self):
        with pytest.raises(ValueError):
            sc._primer_pair_qc(_CLEAN_FWD, _CLEAN_REV, tm_profile="nope")


# ═══════════════════════════════════════════════════════════════════════════════
# In-silico PCR: rank, THEN cap
# ═══════════════════════════════════════════════════════════════════════════════

def _site(strand, foot, ident, length=20):
    return {"strand": strand, "foot_start": foot, "length": length,
            "ident_pct": ident, "mismatches": 0}


def _brute_amplicons(sites_a, sites_b, total, circular, max_amp):
    """Independent reference: every forward×reverse pair, the same geometry
    written out longhand, sorted, no cap."""
    out = {}
    pooled = [(s, i) for i, lst in enumerate((sites_a, sites_b)) for s in lst]
    for fs, fi in pooled:
        if fs["strand"] != 1:
            continue
        for rs, ri in pooled:
            if rs["strand"] != -1:
                continue
            if circular:
                span = (rs["foot_start"] - fs["foot_start"]) % total \
                    + rs["length"]
            elif rs["foot_start"] >= fs["foot_start"] \
                    and rs["foot_start"] + rs["length"] <= total:
                span = rs["foot_start"] + rs["length"] - fs["foot_start"]
            else:
                continue
            if span < fs["length"] + rs["length"] or span > max_amp:
                continue
            key = (fs["foot_start"] % total,
                   (fs["foot_start"] + span) % total, span)
            out.setdefault(key, (min(fs["ident_pct"], rs["ident_pct"]),
                                 span, fs["foot_start"] % total))
    return sorted(out.values(), key=lambda c: (-c[0], c[1], c[2]))


class TestInsilicoRankThenCap:
    def test_best_product_survives_many_worse_ones(self):
        """The best forward site (100 %) comes SECOND in scan order behind a
        forward site with 60 mediocre partners. Capping in scan order (the old
        behaviour) filled all 50 slots with 80 %-certainty products and never
        reached the 99 % one."""
        total = 10_000
        f_many = _site(1, 6_000, 100.0)
        f_good = _site(1, 5_000, 99.0)          # sorts after f_many
        revs = [_site(-1, 6_300 + 40 * i, 80.0) for i in range(60)]
        r_good = _site(-1, 5_400, 100.0)        # upstream of f_many (linear)
        sites_a = [f_many, f_good]
        sites_b = revs + [r_good]
        amps = sc._insilico_pcr_amplicons(sites_a, sites_b, total,
                                          circular=False)
        assert len(amps) == sc._PCR_MAX_AMPLICONS
        assert amps[0]["certainty"] == 99.0
        assert (amps[0]["start"], amps[0]["length"]) == (5_000, 420)

    @pytest.mark.parametrize("seed", range(8))
    @pytest.mark.parametrize("circular", [False, True])
    def test_matches_brute_force(self, seed, circular):
        rnd = random.Random(seed)
        total = rnd.randint(500, 3_000)

        def sites(k):
            return [_site(rnd.choice((1, -1)), rnd.randrange(total),
                          rnd.choice((100.0, 95.0, 90.0, 85.0, 70.0)),
                          rnd.randint(15, 30)) for _ in range(k)]
        sa, sb = sites(rnd.randint(0, 40)), sites(rnd.randint(0, 40))
        want = _brute_amplicons(sa, sb, total, circular, 2_000)
        for cap in (5, 50, 10_000):
            got = sc._insilico_pcr_amplicons(sa, sb, total, circular=circular,
                                             max_amplicon=2_000,
                                             max_amplicons=cap)
            assert [(a["certainty"], a["length"], a["start"]) for a in got] \
                == want[:cap]


# ═══════════════════════════════════════════════════════════════════════════════
# check-primer-pair
# ═══════════════════════════════════════════════════════════════════════════════

_T = _rand_dna(1200, 101)
_F_ARM = _T[100:122]
_R_ARM = _rc(_T[700:722])
_F_TAIL = "GCGCGAATTC"
_R_TAIL = "GCGCGGATCCTT"


def _pair(**kw):
    body = {"forward_primer": _F_TAIL + _F_ARM, "reverse_primer": _R_TAIL + _R_ARM,
            "forward_binding": _F_ARM, "reverse_binding": _R_ARM,
            "template": _T, "circular": False}
    body.update(kw)
    return sc._h_check_primer_pair(None, {k: v for k, v in body.items()
                                          if v is not ...})


class TestCheckPrimerPair:
    def test_tailed_pair(self):
        r = _pair()
        assert r["ok"] and r["tm_profile"] == "primer3_default"
        f, rv, p = r["forward"], r["reverse"], r["pair"]
        assert f["tail"] == _F_TAIL and f["binding"] == _F_ARM
        assert f["tm"] == sc._primer_tm(_F_ARM)
        assert rv["tm"] == sc._primer_tm(_R_ARM)
        assert f["n_sites"] == rv["n_sites"] == 1 and f["best_identity"] == 100
        assert p["n_amplicons"] == 1 and p["single_product"] is True
        prod = p["product"]
        assert (prod["start"], prod["end"], prod["template_bp"]) == (100, 722, 622)
        assert prod["product_bp"] == 622 + len(_F_TAIL) + len(_R_TAIL)
        assert prod["primers"] == "forward+reverse"
        assert r["template_bp"] == len(_T) and r["circular"] is False
        assert r["warnings"] == []

    def test_product_bp_is_what_the_simulator_builds(self):
        # The PCR simulator splices the tails onto the product; ours must
        # agree on its length.
        sim = sc._simulate_pcr(_T, _F_TAIL + _F_ARM, _R_TAIL + _R_ARM)
        assert len(sim) == 1
        assert _pair()["pair"]["product"]["product_bp"] == sim[0]["length"]

    def test_no_arm_given_scores_the_tail_and_says_so(self):
        r = _pair(forward_binding=..., reverse_binding=...)
        f = r["forward"]
        assert f["binding"] == f["primer"] and f["best_identity"] < 100
        assert f["tm"] == sc._primer_tm(_F_TAIL + _F_ARM)
        assert any("give forward_binding" in w for w in r["warnings"])
        # Footprint scored over the whole oligo already spans the tail.
        assert r["pair"]["product"]["product_bp"] \
            == r["pair"]["product"]["template_bp"]

    def test_aliases(self):
        r = sc._h_check_primer_pair(None, {
            "fwd": _F_ARM, "rev_primer": _R_ARM, "template_seq": _T,
            "circular": False, "tm_profile": "Benchling"})
        assert r["ok"] and r["tm_profile"] == "benchling_compatible"
        assert r["forward"]["tm"] == sc._primer_tm(_F_ARM,
                                                   "benchling_compatible")
        r2 = sc._h_check_primer_pair(None, {
            "forward": _F_ARM, "reverse": _R_ARM, "sequence": _T,
            "fwd_binding": _F_ARM[-15:]})
        assert r2["forward"]["binding"] == _F_ARM[-15:]

    @pytest.mark.parametrize(("patch", "code", "needle"), [
        ({"forward_primer": ...}, 400, "forward_primer"),
        ({"forward_primer": 42}, 400, "forward_primer"),
        ({"reverse_primer": "ACGTZZACGT"}, 400, "not DNA bases"),
        ({"reverse_primer": "fwd-1: ACGTACGTACGTACGT"}, 400, "label"),
        ({"forward_binding": _F_ARM[:-1]}, 400, "3' suffix"),
        ({"forward_binding": ""}, 400, "non-empty"),
        ({"forward_binding": "   "}, 400, "non-empty"),
        ({"forward_binding": _F_ARM[-9:]}, 400, "at least 10"),
        ({"forward_primer": "ACGTACGTA", "forward_binding": ...}, 400,
         "at least 10"),
        ({"forward_primer": "A" * 1001, "forward_binding": ...}, 413,
         "too long"),
        ({"template": "ACGT" * 50_001}, 413, "cap"),
        ({"template": 7}, 400, "string"),
        ({"template": "NOT:A:TEMPLATE"}, 400, "label"),
        ({"max_amplicon": 0}, 400, "max_amplicon"),
        ({"max_amplicon": "lots"}, 400, "max_amplicon"),
        ({"tm_profile": "snappy"}, 400, "tm_profile"),
    ])
    def test_refusals(self, patch, code, needle):
        res = _pair(**patch)
        assert isinstance(res, tuple), res
        err, got = res
        assert got == code and needle in err["error"], err

    def test_u_is_read_as_t_but_a_header_is_not_bases(self):
        r = _pair(forward_primer=_F_TAIL + _F_ARM.replace("T", "U"),
                  forward_binding=...)
        assert r["forward"]["read_u_as_t"] is True
        assert "read_u_as_t" not in r["reverse"]
        r2 = _pair(forward_primer=">UNIQUE_PRIMER\n" + _F_TAIL + _F_ARM)
        assert "read_u_as_t" not in r2["forward"]

    def test_defaults_to_the_loaded_plasmid_and_its_topology(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        # Put the forward arm across the origin: only a CIRCULAR read finds a
        # product, so this also proves the record's topology was used.
        rot = _T[110:] + _T[:110]
        rec = SeqRecord(Seq(rot), id="x",
                        annotations={"topology": "circular"})
        app = types.SimpleNamespace(_current_record=rec)
        r = sc._h_check_primer_pair(app, {
            "forward_primer": _F_ARM, "reverse_primer": _R_ARM})
        assert r["template"] == "the loaded plasmid" and r["circular"] is True
        prod = r["pair"]["product"]
        assert prod["wraps"] is True and prod["template_bp"] == 622
        assert prod["start"] == len(rot) - 10
        rec.annotations["topology"] = "linear"
        lin = sc._h_check_primer_pair(app, {
            "forward_primer": _F_ARM, "reverse_primer": _R_ARM})
        assert lin["circular"] is False and lin["pair"]["product"] is None

    def test_no_template_anywhere_is_thermodynamics_only(self):
        app = types.SimpleNamespace(_current_record=None)
        r = sc._h_check_primer_pair(app, {"forward_primer": _CLEAN_FWD,
                                          "reverse_primer": _CLEAN_REV})
        assert r["ok"] and r["template_bp"] is None and r["circular"] is None
        assert r["forward"]["binds"] is None and r["forward"]["sites"] == []
        assert r["pair"]["n_amplicons"] is None
        assert r["pair"]["single_product"] is None
        assert "binding and products were not checked" in r["note"]
        assert r["pair"]["heterodimer_dg"] is not None

    def test_product_ending_on_the_origin_reads_n_not_zero(self):
        f_arm, r_arm = _T[400:422], _rc(_T[len(_T) - 22:])
        r = sc._h_check_primer_pair(None, {
            "forward_primer": f_arm, "reverse_primer": r_arm,
            "template": _T, "circular": True})
        assert r["pair"]["product"]["end"] == len(_T)

    def test_swapped_names_are_still_the_pairs_product(self):
        r = _pair(forward_primer=_R_ARM, reverse_primer=_F_ARM,
                  forward_binding=..., reverse_binding=...)
        assert r["pair"]["product"]["primers"] == "reverse+forward"
        assert r["pair"]["single_product"] is True

    def test_a_primer_priming_both_strands_is_named(self):
        # The forward arm's reverse complement sits downstream as well, so
        # the forward primer alone makes a product.
        tmpl = _T[:900] + _rc(_F_ARM) + _T[922:]
        r = _pair(template=tmpl)
        kinds = {a["primers"] for a in r["pair"]["amplicons"]}
        assert "forward+forward" in kinds
        assert r["pair"]["single_product"] is False
        assert any("the forward primer on its own" in w
                   for w in r["warnings"])

    def test_no_product_is_explained(self):
        # Both arms on the SAME strand: they bind, but face the same way.
        r = _pair(reverse_primer=_T[700:722], reverse_binding=...)
        assert r["forward"]["binds"] and r["reverse"]["binds"]
        assert r["pair"]["product"] is None
        assert r["pair"]["single_product"] is False
        assert any("make no product together" in w for w in r["warnings"])

    def test_non_binding_primer(self):
        r = _pair(reverse_primer=_rand_dna(22, 9), reverse_binding=...)
        assert r["reverse"]["binds"] is False
        assert any("the reverse primer does not anneal" in w
                   for w in r["warnings"])

    def test_misprime_site_warned(self):
        # A second, one-mismatch copy of the forward arm further in.
        near = _F_ARM[:3] + ("A" if _F_ARM[3] != "A" else "C") + _F_ARM[4:]
        tmpl = _T[:300] + near + _T[322:]
        r = _pair(template=tmpl)
        assert r["forward"]["n_sites"] == 2
        assert any("also matches 1 other site(s) ≥80%" in w
                   for w in r["warnings"])

    def test_structure_flags_reach_the_endpoint(self):
        r = sc._h_check_primer_pair(None, {
            "forward_primer": _HAIRPIN_PRIMER,
            "reverse_primer": _HAIRPIN_PARTNER})
        assert any("hairpin that melts at" in w for w in r["warnings"])
        assert r["dimer_tm_limit"] == 47.0 and r["anneal_ref_c"] > 40


# ═══════════════════════════════════════════════════════════════════════════════
# design-primers / check-primer: the profile everywhere, the QC attached
# ═══════════════════════════════════════════════════════════════════════════════

_DT = _rand_dna(900, 202)


class TestDesignPrimersProfile:
    @pytest.mark.parametrize("mode", ["detection", "cloning", "generic",
                                      "allele_specific"])
    def test_unknown_profile_is_refused_in_every_mode(self, mode):
        # allele_specific is dispatched early — it used to skip the check and
        # quietly design under the default conditions.
        r = sc._h_design_primers(None, {
            "template": _DT, "start": 0, "end": 600, "mode": mode,
            "site_5": "GAATTC", "site_3": "GGATCC", "variant_pos": 400,
            "tm_profile": "snappy"})
        assert isinstance(r, tuple) and r[1] == 400
        assert "tm_profile" in r[0]["error"]

    def test_allele_specific_on_the_profiles_scale(self):
        d = sc._h_design_primers(None, {"template": _DT, "variant_pos": 400,
                                        "mode": "allele_specific"})
        b = sc._h_design_primers(None, {"template": _DT, "variant_pos": 400,
                                        "mode": "allele_specific",
                                        "tm_profile": "benchling"})
        assert d["tm_profile"] == "primer3_default"
        assert b["tm_profile"] == "benchling_compatible"
        assert "tm_profile" not in b["ignored"]
        seq = b["result"]["seq"]
        assert b["result"]["tm"] == pytest.approx(
            primer3.calc_tm(seq, **_BENCH))
        # The default design is exactly the designer's own.
        direct = sc._design_allele_specific_primer(_DT, 400)
        assert d["result"]["seq"] == direct["seq"]
        assert d["result"]["tm"] == direct["tm"]

    @pytest.mark.parametrize("mode", ["detection", "cloning", "generic"])
    def test_design_carries_the_pair_qc(self, mode):
        r = sc._h_design_primers(None, {
            "template": _DT, "start": 0, "end": 800, "mode": mode,
            "site_5": "GAATTC", "site_3": "GGATCC"})
        res, qc = r["result"], r["qc"]
        assert qc and isinstance(qc["warnings"], list)
        if mode == "cloning":
            assert qc["forward"]["primer"] == res["fwd_full"]
            assert qc["forward"]["binding"] == res["fwd_binding"]
            assert qc["forward"]["tail"] == "GCGC" + "GAATTC"
        else:
            assert qc["forward"]["primer"] == res["fwd_seq"]
        assert qc["forward"]["tm"] == res["fwd_tm"]

    def test_default_designs_are_unchanged(self):
        r = sc._h_design_primers(None, {"template": _DT, "start": 0,
                                        "end": 800})
        direct = sc._design_detection_primers(_DT, 0, 800)
        assert r["result"] == direct and r["tm_profile"] == "primer3_default"


class TestCheckPrimerProfile:
    def test_alias_and_echo(self):
        r = sc._h_check_primer(None, {"primer": _CLEAN_FWD,
                                      "tm_profile": "Benchling"})
        assert r["tm_profile"] == "benchling_compatible"
        assert r["tm"] == sc._primer_tm(_CLEAN_FWD, "benchling_compatible")

    def test_null_profile_is_the_default(self):
        r = sc._h_check_primer(None, {"primer": _CLEAN_FWD,
                                      "tm_profile": None})
        assert r["tm_profile"] == "primer3_default"

    def test_header_u_is_not_read_as_a_base(self):
        r = sc._h_check_primer(None, {"primer": ">UP_1\n" + _CLEAN_FWD,
                                      "template": _CLEAN_FWD * 3})
        assert "read_u_as_t" not in r

    def test_template_less_answer_reports_u_too(self):
        app = types.SimpleNamespace(_current_record=None)
        r = sc._h_check_primer(app, {"primer": _CLEAN_FWD.replace("T", "U")})
        assert r["binds"] is None and r["read_u_as_t"] is True


# ═══════════════════════════════════════════════════════════════════════════════
# The Primer Design screen
# ═══════════════════════════════════════════════════════════════════════════════

class TestDesignResultQcHelper:
    def test_each_mode_names_its_oligos(self):
        tail_f, tail_r = "GCGCGGTCTCAGGAG", "GCGCGGTCTCAAAGC"
        gb = {"fwd_full": tail_f + _F_ARM, "rev_full": tail_r + _R_ARM,
              "fwd_binding": _F_ARM, "rev_binding": _R_ARM}
        q = sc._design_result_qc("goldenbraid", gb)
        assert q["forward"]["tail"] == tail_f and q["reverse"]["tail"] == tail_r
        q2 = sc._design_result_qc("generic", {"fwd_seq": _F_ARM,
                                              "rev_seq": _R_ARM})
        assert q2["forward"]["tail"] == ""
        assert sc._design_result_qc("cloning", {"fwd_full": "ACGT"}) is None

    def test_text_block(self):
        t = sc._primer_qc_text(sc._primer_pair_qc(_CLEAN_FWD, _CLEAN_REV))
        assert "Pair check" in t.plain and "✓" in t.plain
        assert "ΔTm" in t.plain and "primer-dimer" in t.plain
        bad = sc._primer_qc_text(
            sc._primer_pair_qc(_HAIRPIN_PRIMER, _HAIRPIN_PARTNER))
        assert "⚠" in bad.plain and "hairpin that melts at" in bad.plain
        one = sc._primer_qc_text(sc._primer_oligo_qc(_CLEAN_FWD),
                                 whole_oligo=True)
        assert "Primer check" in one.plain and "(whole oligo)" in one.plain
        # Warning text is data, never markup.
        weird = {"pair": {"tm_delta": None}, "forward": {}, "reverse": {},
                 "warnings": ["[bold red]not markup[/]"]}
        assert "[bold red]not markup[/]" in sc._primer_qc_text(weird).plain


class TestPrimerDesignScreenQc:
    async def test_design_result_shows_the_pair_check(self):
        from textual.widgets import Static, Tabs
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen(_DT, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            # Activate the tab the way a click does — its handler clears any
            # previous design, so it has to settle BEFORE the design runs.
            screen.query_one("#pd-mode-tabs", Tabs).active = "tab-generic"
            await pilot.pause()
            assert screen._current_mode() == "generic"
            screen._run_generic(_DT, 0, 600)
            await app.workers.wait_for_complete()
            await pilot.pause()
            assert screen._det_result is not None
            qc = screen._det_result.get("qc")
            assert qc and qc["forward"]["primer"] == screen._det_result["fwd_seq"]
            txt = str(screen.query_one("#pd-results", Static).render())
            assert "Pair check" in txt and "ΔG kcal/mol" in txt

    async def test_primer_check_tab_shows_qc_before_the_scan(self, monkeypatch):
        from textual.widgets import Input, Static
        monkeypatch.setattr(sc, "_iter_library_readonly", lambda: [])
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen._switch_mode("primercheck")
            qc = screen.query_one("#pd-pc-qc", Static)
            screen.query_one("#pd-pc-p1", Input).value = _HAIRPIN_PRIMER
            screen.query_one("#pd-pc-p2", Input).value = _HAIRPIN_PARTNER
            screen._pc_go()
            await app.workers.wait_for_complete()
            await pilot.pause()
            txt = str(qc.render())
            assert "Pair check" in txt and "hairpin that melts at" in txt
            # One primer → the single-oligo check.
            screen.query_one("#pd-pc-p2", Input).value = ""
            screen._pc_go()
            await app.workers.wait_for_complete()
            await pilot.pause()
            assert "Primer check" in str(qc.render())
            # A refused input clears it — no stale QC beside a red error.
            screen.query_one("#pd-pc-p1", Input).value = "ACGT"
            screen._pc_go()
            await pilot.pause()
            assert str(qc.render()).strip() == ""
            await app.workers.wait_for_complete()

    async def test_the_ui_thread_never_computes_the_qc(self, monkeypatch):
        # primer3's alignment shares a lock with a running Primer3 design; the
        # screen must never wait on it. Pressing Scan must not run the QC on
        # the UI thread — the scan worker does.
        import threading as _t
        from textual.widgets import Input
        monkeypatch.setattr(sc, "_iter_library_readonly", lambda: [])
        seen: list = []
        real = sc._primer_pair_qc

        def spy(*a, **k):
            seen.append(_t.current_thread() is _t.main_thread())
            return real(*a, **k)
        monkeypatch.setattr(sc, "_primer_pair_qc", spy)
        app = sc.PlasmidApp()
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            screen = sc.PrimerDesignScreen("ACGT" * 200, [], "test")
            app.push_screen(screen)
            await pilot.pause()
            screen._switch_mode("primercheck")
            screen.query_one("#pd-pc-p1", Input).value = _CLEAN_FWD
            screen.query_one("#pd-pc-p2", Input).value = _CLEAN_REV
            screen._pc_go()
            await app.workers.wait_for_complete()
            await pilot.pause()
        assert seen == [False]
