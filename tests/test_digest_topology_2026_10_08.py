"""Regression tests for `[INV-215]` — a digest of a LINEAR molecule run with
the `circular=True` default used to be indistinguishable from a real circular
digest in everything but the fragment COUNT.

From a round-3 agent field report (2026-10-08). The reporter's own repro is
the point: digest a 76 bp synthetic fragment with and without
`circular=False` and the insert-bearing internal fragment comes back
BYTE-IDENTICAL either way — only the fragment count differs (3 vs 2, because
the two end pieces fuse into one origin-spanning fragment). So the check
anyone would actually write, "the digest releases my intended insert", passes
whether or not you remembered the flag. A design script could carry a genuine
end-piece error, verify clean, and get ordered.

The fix does NOT change the default: an undeclared topology reads as circular
app-wide (the 2026-09-22 reversal — `_record_is_circular`, the rule the map
already draws). Instead the output stops hiding which molecule it answered
about: every fragment and every cut carries `wraps`, fragments carry absolute
`start`/`end`, and the agent endpoint warns when the topology was ASSUMED and
the answer turns on it.

Every test here fails on the tree as it stood before this change.
"""
from __future__ import annotations

import random

import pytest

import splicecraft as sc
import splicecraft_biology as _bio


# A linear synthetic fragment: 5' tail + EcoRI + payload + EcoRI + 3' tail.
_PAYLOAD = "ATGGCTAGCTAGGCATCAGCTACGATCGATCGTAGCTAGC"
_LINEAR = "ACGTACGTACGTA" + "GAATTC" + _PAYLOAD + "GAATTC" + "TTACGTACGTACGTAC"

# Two ends that spell GAATTC only when JOINED — EcoRI cuts this piece nowhere
# as a linear molecule, and once as a circle.
_JOIN_ONLY = "TTC" + "ACGTACGTACGTACGTAAGGCCTTAACCGGTT" + "GAA"


def _cut(top: int, *, kind: str = "blunt", enzyme: str = "FakeI") -> dict:
    """A hand-built cut for `_fragments_from_cuts` — the coordinate edges
    below need a cut at an exact bp, which no real enzyme obliges."""
    return {"top": top, "bot": top, "kind": kind,
            "overhang_seq": "", "enzyme": enzyme}


# ══════════════════════════════════════════════════════════════════════════
# The reporter's repro: the positive check passes either way
# ══════════════════════════════════════════════════════════════════════════

class TestTheSilentHalfIsTheFragmentCount:
    """Pin the behaviour that made this silent, so nobody "fixes" it by
    changing the default and calls the issue closed."""

    def test_insert_fragment_is_byte_identical_under_both_topologies(self):
        lin = sc._digest_with_enzymes(_LINEAR, ["EcoRI"], circular=False)
        circ = sc._digest_with_enzymes(_LINEAR, ["EcoRI"], circular=True)
        # The count differs...
        assert (len(lin), len(circ)) == (3, 2)
        # ...but the payload-bearing fragment does not, which is why
        # "my insert came out" is not a topology check.
        lin_hit = [f for f in lin if _PAYLOAD in f["top_seq"]]
        circ_hit = [f for f in circ if _PAYLOAD in f["top_seq"]]
        assert len(lin_hit) == len(circ_hit) == 1
        assert lin_hit[0]["top_seq"] == circ_hit[0]["top_seq"]

    def test_the_wrap_flag_is_what_tells_them_apart(self):
        lin = sc._digest_with_enzymes(_LINEAR, ["EcoRI"], circular=False)
        circ = sc._digest_with_enzymes(_LINEAR, ["EcoRI"], circular=True)
        assert [f["wraps"] for f in lin] == [False, False, False]
        assert sum(f["wraps"] for f in circ) == 1
        # The fused fragment is exactly the 3' tail + the 5' tail.
        fused = next(f for f in circ if f["wraps"])
        tail3 = next(f for f in lin if f is lin[-1])
        tail5 = lin[0]
        assert fused["top_seq"] == tail3["top_seq"] + tail5["top_seq"]


# ══════════════════════════════════════════════════════════════════════════
# Fragment provenance
# ══════════════════════════════════════════════════════════════════════════

class TestFragmentsSayWhereTheyCameFrom:
    def test_doubling_identity_locates_every_fragment(self):
        """`(seq+seq)[start : start+len]` is the fragment — including the
        wrapping one. Same contract `amplify-feature` primers already use."""
        for circ in (False, True):
            frags = sc._digest_with_enzymes(_LINEAR, ["EcoRI", "BamHI"],
                                            circular=circ)
            twice = _LINEAR + _LINEAR
            for f in frags:
                s, ts = f["start"], f["top_seq"]
                assert twice[s:s + len(ts)] == ts, (circ, f)

    def test_wrap_is_not_derivable_from_end_lt_start(self):
        """Two fragments both reduce to `end == start` or `end < start`, and
        only one of them is a join — so a caller cannot re-derive the flag,
        which is why it is stored. Driven through `_fragments_from_cuts` with
        hand-built cuts because these are coordinate edges, not enzymes."""
        seq = "ACGT" * 10                     # n = 40
        n = len(seq)

        # (a) ONE cut away from the origin → a full-lap fragment that IS the
        #     join of the 3' and 5' ends. `start == end`, and it wraps.
        lap = _bio._fragments_from_cuts(seq, [_cut(5)], circular=True)
        assert len(lap) == 1
        assert (lap[0]["start"], lap[0]["end"]) == (5, 5)
        assert lap[0]["wraps"] is True
        assert lap[0]["top_seq"] == seq[5:] + seq[:5]

        # (b) ONE cut AT the origin → the same `start == end` shape reduced,
        #     but the fragment is the plain contiguous molecule.
        whole = _bio._fragments_from_cuts(seq, [_cut(0)], circular=True)
        assert len(whole) == 1
        assert (whole[0]["start"], whole[0]["end"]) == (0, n)
        assert whole[0]["wraps"] is False
        assert whole[0]["top_seq"] == seq

    def test_a_fragment_ending_at_the_origin_does_not_wrap(self):
        """A cut at bp 0 makes the preceding fragment a contiguous 3'-end
        slice: it reduces to `end < start` yet joins nothing. Marking it a
        wrap would be a false alarm on a genuinely linear piece."""
        seq = "ACGT" * 10
        frags = _bio._fragments_from_cuts(seq, [_cut(0), _cut(20)],
                                          circular=True)
        tail = next(f for f in frags if f["start"] == 20)
        assert (tail["start"], tail["end"]) == (20, len(seq))
        assert tail["wraps"] is False
        assert tail["top_seq"] == seq[20:]
        assert not any(f["wraps"] for f in frags)

    def test_at_most_one_fragment_can_span_the_origin(self):
        frags = sc._digest_with_enzymes(
            _LINEAR, ["EcoRI", "BamHI", "HindIII"], circular=True)
        assert sum(f["wraps"] for f in frags) <= 1

    def test_linear_digest_never_claims_a_wrap(self):
        for enz in (["EcoRI"], ["BsaI"], ["EcoRI", "BamHI"]):
            for f in sc._digest_with_enzymes(_LINEAR, enz, circular=False):
                assert f["wraps"] is False
                assert f["end"] >= f["start"]

    def test_uncut_molecule_is_the_whole_input(self):
        for circ in (False, True):
            frags = sc._digest_with_enzymes(_LINEAR, ["NotI"], circular=circ)
            assert len(frags) == 1
            f = frags[0]
            assert (f["start"], f["end"], f["wraps"]) == (0, len(_LINEAR),
                                                          False)


# ══════════════════════════════════════════════════════════════════════════
# The phantom cut — the sharper half
# ══════════════════════════════════════════════════════════════════════════

class TestCutsThatOnlyExistOnACircle:
    def test_a_site_spanning_the_join_is_marked(self):
        """Ends that spell GAATTC only when joined. On a linear piece EcoRI
        cuts nowhere; as a circle it cuts once — and that cut used to be
        reported indistinguishably from a real one, so "is there an internal
        EcoRI site?" answered YES for a piece the enzyme never touches."""
        assert sc._enzyme_cuts(_JOIN_ONLY, ["EcoRI"], circular=False) == []
        circ = sc._enzyme_cuts(_JOIN_ONLY, ["EcoRI"], circular=True)
        assert len(circ) == 1
        assert circ[0]["wraps"] is True

    def test_an_interior_cut_is_not_marked(self):
        for c in sc._enzyme_cuts(_LINEAR, ["EcoRI"], circular=True):
            assert c["wraps"] is False

    def test_linear_scan_never_marks_a_cut(self):
        for enz in (["EcoRI"], ["BsaI"], ["BaeI"], ["DpnII"]):
            for c in sc._enzyme_cuts(_LINEAR, enz, circular=False):
                assert c["wraps"] is False

    def test_a_bond_two_enzymes_make_is_real_if_either_site_is_interior(self):
        """The coincident-bond merge must AND, not OR: one enzyme reaching a
        bond from a site wholly inside the molecule makes that cut real on a
        linear molecule too, whatever the other site's geometry."""
        # DpnII / Sau3AI / MboI all sever GATC — interior here.
        seq = "ACGTACGT" + "GATC" + "ACGTACGTACGTACGTACGT"
        cuts = sc._enzyme_cuts(seq, ["DpnII", "Sau3AI", "MboI"], circular=True)
        assert cuts, "GATC should be cut"
        for c in cuts:
            assert c["wraps"] is False, c


# ══════════════════════════════════════════════════════════════════════════
# Fuzz: the contract holds on every rotation / topology / enzyme geometry
# ══════════════════════════════════════════════════════════════════════════

class TestContractUnderFuzz:
    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_partition_and_wrap_contract(self, seed):
        rng = random.Random(seed)
        enz_sets = [["EcoRI"], ["BsaI"], ["DpnII"], ["BaeI"], ["AlwNI"],
                    ["EcoRI", "BamHI", "HindIII", "BsaI"]]
        saw_wrap = 0
        for _ in range(400):
            n = rng.randint(12, 160)
            s = "".join(rng.choice("ACGT") for _ in range(n))
            enz = rng.choice(enz_sets)
            circ = rng.random() < 0.5
            frags = sc._digest_with_enzymes(s, enz, circular=circ)
            twice = s + s
            total = 0
            for f in frags:
                st, ts, w = f["start"], f["top_seq"], f["wraps"]
                assert twice[st:st + len(ts)] == ts
                assert 0 <= st < n and 0 <= f["end"] <= n
                assert w == (st + len(ts) > n)
                assert f["end"] == (st + len(ts)) - (n if w else 0)
                if not circ:
                    assert not w
                total += len(ts)
                saw_wrap += bool(w)
            # A digest partitions the molecule exactly once.
            assert total == n, (s, enz, circ)
            assert sum(1 for f in frags if f["wraps"]) <= 1
            if circ:
                lin_bonds = {(c["top"], c["bot"])
                             for c in sc._enzyme_cuts(s, enz, circular=False)}
                for c in sc._enzyme_cuts(s, enz, circular=True):
                    # Every cut NOT marked `wraps` must survive linearisation;
                    # that is the flag's definition, checked against the real
                    # linear scan rather than against itself.
                    if not c["wraps"]:
                        assert (c["top"], c["bot"]) in lin_bonds, (s, enz, c)
        assert saw_wrap, "fuzz never produced a wrapping fragment"


# ══════════════════════════════════════════════════════════════════════════
# The agent endpoint
# ══════════════════════════════════════════════════════════════════════════

class TestDigestEndpointSurfacesTopology:
    def test_assumed_topology_that_changes_the_answer_warns(self):
        out = sc._h_digest(None, {"sequence": _LINEAR, "enzymes": ["EcoRI"]})
        assert out["circular"] is True
        warn = " ".join(out["warnings"]).lower()
        assert "circular" in warn and "linear" in warn
        assert any(f["wraps"] for f in out["fragments"])

    def test_declared_linear_is_not_warned_and_gives_three_fragments(self):
        out = sc._h_digest(None, {"sequence": _LINEAR, "enzymes": ["EcoRI"],
                                  "circular": False})
        assert out["n_fragments"] == 3
        assert "warnings" not in out
        assert not any(f["wraps"] for f in out["fragments"])

    def test_declared_circular_is_not_scolded(self):
        """A user who said `circular: true` knows what they have. Warning on
        every plasmid digest would just train the reader to skip `warnings`."""
        out = sc._h_digest(None, {"sequence": _LINEAR, "enzymes": ["EcoRI"],
                                  "circular": True})
        assert any(f["wraps"] for f in out["fragments"])
        assert "warnings" not in out

    def test_phantom_cut_is_named_in_the_warning(self):
        out = sc._h_digest(None, {"sequence": _JOIN_ONLY,
                                  "enzymes": ["EcoRI"]})
        assert out["cuts"][0]["wraps"] is True
        assert any("EcoRI" in w and "circular" in w
                   for w in out["warnings"]), out["warnings"]

    def test_fragments_carry_locatable_coordinates(self):
        out = sc._h_digest(None, {"sequence": _LINEAR, "enzymes": ["EcoRI"],
                                  "include_fragment_seq": True})
        twice = _LINEAR + _LINEAR
        for f in out["fragments"]:
            assert twice[f["start"]:f["start"] + f["length"]] == f["seq"]
            assert f["length"] == len(f["seq"])

    def test_loaded_record_topology_is_declared_not_assumed(self, monkeypatch):
        """Digesting the canvas must stay warning-free: the record's own
        topology (or the settled unannotated⇒circular rule) is a declaration,
        not a guess. Re-litigating that is the 2026-09-22 reversal."""
        class _Rec:
            seq = _LINEAR
            annotations: dict = {}

        out = sc._h_digest(type("A", (), {"_current_record": _Rec()})(),
                           {"enzymes": ["EcoRI"]})
        assert out["sequence_from"] == "the loaded plasmid"
        assert any(f["wraps"] for f in out["fragments"])
        assert "warnings" not in out


class TestSimulatePcrEchoesTopology:
    """`template_seq` is always body-supplied, so the circular default is an
    assumption every time — and a wrapping amplicon on a linear template is a
    product no tube contains."""

    _MID = "".join("ACGTTGCAAGGCTTAACCGGATCCGTAATTCAGGCCTTAAGGCCAATTCCGGTTAACC"[i % 58]
                   for i in range(300))
    _FWD = "GGCATTACGCATTGGCCATA"
    _REV_SITE = "TTGGCCAATTCGCAATTGCC"
    _TMPL = _REV_SITE + _MID + _FWD

    @property
    def _rev(self):
        return sc._rc(self._REV_SITE)

    def test_response_echoes_the_topology_it_simulated(self):
        for circ in (True, False):
            out = sc._h_simulate_pcr(None, {"template_seq": self._TMPL,
                                            "fwd_primer": self._FWD,
                                            "rev_primer": self._rev,
                                            "circular": circ})
            assert out["circular"] is circ

    def test_assumed_topology_with_a_wrapping_amplicon_warns(self):
        out = sc._h_simulate_pcr(None, {"template_seq": self._TMPL,
                                        "fwd_primer": self._FWD,
                                        "rev_primer": self._rev})
        assert out["n"] == 1 and out["amplicons"][0]["wraps"] is True
        assert any("circular" in w for w in out["warnings"])

    def test_declared_linear_yields_no_product_and_no_warning(self):
        out = sc._h_simulate_pcr(None, {"template_seq": self._TMPL,
                                        "fwd_primer": self._FWD,
                                        "rev_primer": self._rev,
                                        "circular": False})
        assert out["n"] == 0
        assert "warnings" not in out
