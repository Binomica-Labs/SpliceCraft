"""Allele-specific (ARMS) PCR primer design — `design-primers` mode
``allele_specific``.

Allele-specific PCR IS the act of pinning the primer's 3'-terminal base on the
variant: a matched 3' end extends, a mismatched one does not. None of the
existing designers could express that — `_design_detection_primers` and
`_design_generic_primers` slide freely to whatever window scores best, and
`_mut_design_fwd_anneal` / `_mut_design_rev_anneal` anchor the FIVE-prime end
(they build Golden-Gate SDM primers whose 5' tail carries the site), so sliding
one of those never fixes which base lands last.

The defining property is therefore what gets asserted: the last base IS the
allele, and the primer differs from the reference at exactly that position.
"""

from __future__ import annotations

import random

import pytest

import splicecraft as sc


_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
_TRANSVERSION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def _rand(n: int, seed: int) -> str:
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _rc(s: str) -> str:
    return "".join(_COMP[c] for c in reversed(s))


TEMPLATE = _rand(1200, 21)
POS = 600
REF_BASE = TEMPLATE[POS]
ALT_BASE = _TRANSVERSION[REF_BASE]


def _design(body):
    payload = {"template": TEMPLATE, "mode": "allele_specific"}
    payload.update(body)
    return sc._state._AGENT_HANDLERS["design-primers"][0](None, payload)


class TestThreePrimePinning:
    """The one property that makes it allele-specific at all."""

    def test_the_last_base_is_the_requested_allele(self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE})
        b = res["result"]
        want = (ALT_BASE if b["orientation"] == "fwd"
                else _COMP[ALT_BASE])
        assert b["seq"][-1] == want

    def test_the_primer_differs_from_the_reference_only_at_the_variant(self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "orientation": "fwd"})
        b = res["result"]
        native = TEMPLATE[POS - b["length"] + 1: POS + 1]
        diffs = [i for i, (x, y) in enumerate(zip(b["seq"], native)) if x != y]
        assert diffs == [b["length"] - 1]

    def test_the_reverse_orientation_anneals_its_three_prime_end_to_the_variant(
            self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "orientation": "rev"})
        b = res["result"]
        assert b["orientation"] == "rev"
        assert b["seq"][-1] == _COMP[ALT_BASE]
        # Reverse-complementing it back must reproduce the template window
        # with only the variant changed.
        win = _rc(b["seq"])
        native = TEMPLATE[POS: POS + b["length"]]
        diffs = [i for i, (x, y) in enumerate(zip(win, native)) if x != y]
        assert diffs == [0]

    def test_length_is_swept_within_the_documented_range(self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE})
        assert 18 <= res["result"]["length"] <= 27

    def test_auto_orientation_picks_the_better_scoring_side(self):
        auto = _design({"variant_pos": POS, "alt_base": ALT_BASE})
        fwd = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "orientation": "fwd"})
        rev = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "orientation": "rev"})
        assert auto["result"]["score"] == pytest.approx(
            min(fwd["result"]["score"], rev["result"]["score"]))

    def test_omitting_alt_base_targets_the_reference_allele(self):
        res = _design({"variant_pos": POS})
        assert res["result"]["matches_template"] is True
        assert res["result"]["allele"] == REF_BASE


class TestGcClamp:
    """An A/T discriminating base is a CONSTRAINT OF THE ALLELE, not a
    property of a badly chosen primer — so it is stated, not merely scored
    down. Quietly ranking it lower leaves the caller wondering why every
    candidate looked bad."""

    def _pos_with_base(self, base):
        return next(i for i in range(400, 800) if TEMPLATE[i] == base)

    def test_an_at_allele_is_reported_in_words(self):
        res = _design({"variant_pos": self._pos_with_base("G"),
                       "alt_base": "A"})
        b = res["result"]
        assert b["gc_clamp"] is False
        assert b["clamp_note"] and "GC clamp" in b["clamp_note"]
        assert any("GC clamp" in w for w in res["warnings"])

    def test_a_gc_allele_has_a_clamp_and_no_note(self):
        res = _design({"variant_pos": self._pos_with_base("A"),
                       "alt_base": "G"})
        b = res["result"]
        assert b["gc_clamp"] is True
        assert b["clamp_note"] is None
        assert not any("GC clamp" in w for w in res["warnings"])

    def test_the_note_names_the_actual_three_prime_base(self):
        res = _design({"variant_pos": self._pos_with_base("G"),
                       "alt_base": "A"})
        b = res["result"]
        assert f"is {b['three_prime_base']}" in b["clamp_note"]


class TestSpecificity:

    def test_an_alt_allele_primer_does_not_seed_on_the_reference(self):
        """That failure to bind IS the discrimination, so it is reported as a
        property rather than left as an unexplained empty site list."""
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE})
        assert res["discriminates"] is True
        assert res["n_intended_sites"] == 0

    def test_a_reference_allele_primer_does_seed_on_the_template(self):
        res = _design({"variant_pos": POS})
        assert res["n_intended_sites"] >= 1
        assert res["discriminates"] is False

    def test_a_three_prime_seed_on_a_named_off_target_is_found(self):
        first = _design({"variant_pos": POS, "alt_base": ALT_BASE})
        seed = first["result"]["seq"][-12:]
        decoy = _rand(200, 60) + seed + _rand(200, 61)
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "off_targets": {"pDECOY": decoy}})
        by_name = {s["name"]: s for s in res["specificity"]}
        assert by_name["pDECOY"]["n_sites"] >= 1
        assert res["n_offtarget_sites"] >= 1
        assert any("3'-seeded" in w for w in res["warnings"])

    def test_a_clean_off_target_raises_no_warning(self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "off_targets": {"pCLEAN": _rand(400, 62)}})
        by_name = {s["name"]: s for s in res["specificity"]}
        assert by_name["pCLEAN"]["n_sites"] == 0
        assert not any("3'-seeded" in w for w in res["warnings"])

    def test_off_targets_accept_a_list_of_objects(self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "off_targets": [{"name": "pX",
                                        "sequence": _rand(300, 63)}]})
        assert any(s["name"] == "pX" for s in res["specificity"])

    def test_sites_are_tagged_intended_or_not(self):
        res = _design({"variant_pos": POS})
        rows = [r for s in res["specificity"] for r in s["sites"]]
        assert rows and all("intended" in r for r in rows)


class TestDestabilize:
    """Opt-in, because it changes the primer you order."""

    def test_off_by_default(self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE})
        assert res["result"]["destabilized"] is False
        assert res["result"]["destabilize_pos"] is None

    def test_on_request_it_mismatches_the_third_base_from_the_end(self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "destabilize": True})
        b = res["result"]
        assert b["destabilized"] is True
        assert b["destabilize_pos"] == b["length"] - 3

    def test_it_never_touches_the_discriminating_base(self):
        plain = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                         "orientation": "fwd"})["result"]
        destab = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                          "orientation": "fwd",
                          "destabilize": True})["result"]
        assert destab["seq"][-1] == plain["seq"][-1]


class TestPartnerPrimer:

    def test_a_partner_primer_gets_a_product_size(self):
        disc = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                        "orientation": "fwd"})["result"]
        partner = _rc(TEMPLATE[POS + 200: POS + 222])
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE,
                       "orientation": "fwd", "partner_primer": partner})
        assert res["partner"]["n_sites"] >= 1
        assert res["partner"]["strand"] == -1
        assert res["partner"]["product_size"] == pytest.approx(
            222 - (POS - disc["length"] + 1 - POS), abs=260)
        assert res["partner"]["product_size"] > 0

    def test_no_partner_reports_none_rather_than_a_guess(self):
        res = _design({"variant_pos": POS, "alt_base": ALT_BASE})
        assert res["partner"] is None


class TestAlleleSpecificInput:

    @pytest.mark.parametrize("body,code", [
        ({}, 400),                                        # no variant_pos
        ({"variant_pos": 99999}, 400),
        ({"variant_pos": -1}, 400),
        ({"variant_pos": POS, "alt_base": "X"}, 400),
        ({"variant_pos": POS, "alt_base": "AG"}, 400),
        ({"variant_pos": POS, "orientation": "sideways"}, 400),
        ({"variant_pos": POS, "target_tm": 5}, 400),
        ({"variant_pos": POS, "off_targets": "nope"}, 400),
        ({"variant_pos": POS, "off_targets": {"x": "not bases!!"}}, 400),
        ({"variant_pos": POS, "partner_primer": "!!!"}, 400),
    ])
    def test_bad_input_is_refused_with_a_reason(self, body, code):
        res = _design(body)
        assert isinstance(res, tuple) and res[1] == code, res

    def test_a_variant_too_close_to_a_linear_end_is_refused_not_truncated(
            self):
        res = _design({"variant_pos": 2, "orientation": "fwd"})
        assert isinstance(res, tuple) and res[1] == 422
        assert "circular" in res[0]["error"]

    def test_the_same_variant_works_when_the_template_is_circular(self):
        res = _design({"variant_pos": 2, "orientation": "fwd",
                       "circular": True})
        assert not isinstance(res, tuple)
        assert res["result"]["length"] >= 18

    def test_the_other_modes_still_work(self):
        """The new mode is dispatched before start/end are parsed; the
        existing modes must be untouched by that."""
        out = sc._state._AGENT_HANDLERS["design-primers"][0](
            None, {"template": TEMPLATE, "start": 100, "end": 500,
                   "mode": "generic"})
        assert not isinstance(out, tuple), out
        assert out["result"]["fwd_seq"]

    def test_an_unknown_mode_is_still_refused(self):
        out = sc._state._AGENT_HANDLERS["design-primers"][0](
            None, {"template": TEMPLATE, "start": 0, "end": 100,
                   "mode": "nope"})
        assert isinstance(out, tuple) and out[1] == 400
        assert "allele_specific" in out[0]["error"]

    def test_it_is_registered_read_only(self):
        assert sc._state._AGENT_HANDLERS["design-primers"][1] is False
