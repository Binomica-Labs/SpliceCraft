# pyright: reportArgumentType=false, reportCallIssue=false, reportAttributeAccessIssue=false
#
# Tests pass deliberately wrong-typed inputs to `_gel_bands_for_lane`
# (malformed pcr_amplicon, etc.) to verify the new defensive paths, and
# touch BioPython SeqRecord fields (`record.seq`) where pyright's stubs
# under-narrow `Seq | MutableSeq | None`. `pyproject.toml` excludes
# `tests/**` from pyright for the same reason; the file-scope pragma
# keeps editor / harness diagnostics aligned with that policy.
"""
test_simulator — PCR sim + agarose gel physics + SimulatorScreen.

Regression guard for 2026-05-15 feature add. Three pure-function layers
(`_simulate_pcr`, `_agarose_mobility`, `_gel_bands_for_lane` /
`_render_gel_image`) plus the UI screen's smoke construction. Modal
boundary fit is covered by `test_modal_boundaries.py::_MODAL_CASES`.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


# ═══════════════════════════════════════════════════════════════════════════════
# _simulate_pcr — primer-binding model + amplicon enumeration
# ═══════════════════════════════════════════════════════════════════════════════

class TestPcrPrimerSourceToggle:
    """Sweep #37 (2026-05-27): the PCR pane's primer-source toggle
    lets the user pick primers from the saved library or type
    them free-form. Verifies the option builder + the default
    source mode."""

    def test_build_primer_library_options_empty_library(
            self, monkeypatch,
    ):
        import splicecraft as sc
        monkeypatch.setattr(sc, "_load_primers", lambda: [])
        ss = sc.SimulatorScreen("ATCG" * 200, [], "p", "circular")
        opts = ss._build_primer_library_options()
        assert len(opts) == 1
        label, value = opts[0]
        assert "empty" in label.lower()
        assert value == ""

    def test_build_primer_library_options_real_entries(
            self, monkeypatch,
    ):
        import splicecraft as sc
        def _stub():
            return [
                {"name": "p1", "sequence": "ATCGATCGATCG", "type": "Fwd"},
                {"name": "p2", "sequence": "GCTAGCTAGCTA", "type": ""},
                # Empty + non-ACGT entries get dropped.
                {"name": "broken", "sequence": "XXXX"},
                {"name": "noseq",  "sequence": ""},
            ]
        monkeypatch.setattr(sc, "_load_primers", _stub)
        ss = sc.SimulatorScreen("ATCG" * 200, [], "p", "circular")
        opts = ss._build_primer_library_options()
        assert len(opts) == 2
        # Sequence is the value (so the Select handler can write
        # it straight into the Input without an extra lookup).
        assert {v for _label, v in opts} == {
            "ATCGATCGATCG", "GCTAGCTAGCTA",
        }
        # Name + type tag both surface in the label.
        first_label = opts[0][0]
        assert "p1" in first_label
        assert "Fwd" in first_label

    def test_build_primer_library_options_dedupes_by_sequence(
            self, monkeypatch,
    ):
        """Two entries with the same sequence collapse to one
        option so the dropdown doesn't carry visual duplicates."""
        import splicecraft as sc
        def _stub():
            return [
                {"name": "first",  "sequence": "AAACCCGGGTTT"},
                {"name": "second", "sequence": "AAACCCGGGTTT"},
            ]
        monkeypatch.setattr(sc, "_load_primers", _stub)
        ss = sc.SimulatorScreen("ATCG" * 200, [], "p", "circular")
        opts = ss._build_primer_library_options()
        assert len(opts) == 1
        # First-wins (matches `_dedupe_primers_by_sequence`).
        assert "first" in opts[0][0]

    def test_default_source_mode_is_custom(self):
        """Existing free-text UX is preserved by default — the
        library Select only appears when the user explicitly
        switches to Library mode."""
        import splicecraft as sc
        ss = sc.SimulatorScreen("ATCG" * 200, [], "p", "circular")
        assert ss._pcr_source_mode == "custom"


class TestSimulatePcrBasics:
    """Forward/reverse exact-match binding on linear and circular templates."""

    def test_linear_one_amplicon(self):
        # 100 bp template; primers at positions [0, 20) and [80, 100)
        seq = ("ATGCGATCGATCGATCGCGT"   # fwd binding site 0..20
                + "A" * 60
                + "GCATCGTAGCTAGCTGATCG") # rev_rc binding site 80..100
        fwd = "ATGCGATCGATCGATCGCGT"
        rev = sc._rc("GCATCGTAGCTAGCTGATCG")
        amps = sc._simulate_pcr(seq, fwd, rev, circular=False)
        assert len(amps) == 1
        a = amps[0]
        assert a["start"] == 0
        assert a["end"] == 100
        assert a["length"] == 100
        assert a["wraps"] is False
        assert a["amplicon_seq"] == seq

    def test_circular_wrap_around_origin(self):
        # Place fwd near end, rev_rc near start — amplicon must cross origin
        seq = ("A" * 20 + "ATGCGATCGATCGATCGCGT"
                + "A" * 10 + "GCATCGTAGCTAGCTGATCG" + "A" * 30)
        fwd = "GCATCGTAGCTAGCTGATCG"
        rev = sc._rc("ATGCGATCGATCGATCGCGT")
        amps = sc._simulate_pcr(seq, fwd, rev, circular=True,
                                  max_amplicon=200)
        wrap_amps = [a for a in amps if a["wraps"]]
        assert wrap_amps, "expected at least one wrapping amplicon"
        a = wrap_amps[0]
        assert a["start"] == 50
        # rev_rc lands at position 20 in canonical seq; amplicon end on
        # template is 20 + 20 = 40 (right-exclusive)
        assert a["end"] == 40
        assert a["length"] == 90

    def test_mispriming_multiple_amplicons(self):
        # Forward primer appears twice; expect two amplicons
        repeat = "GATCGATCGATCGATCGATC"   # 20 bp
        rev_site = sc._rc("GTACGTACGTACGTACGTAC")  # rev binding seq on top
        seq = (repeat + "A" * 100 + repeat + "A" * 50 + rev_site)
        amps = sc._simulate_pcr(seq, repeat,
                                  sc._rc("GTACGTACGTACGTACGTAC"),
                                  circular=False, max_amplicon=500)
        assert len(amps) == 2
        # Sorted longest first
        assert amps[0]["length"] > amps[1]["length"]

    def test_no_match_returns_empty(self):
        seq = "ATGC" * 100
        amps = sc._simulate_pcr(seq, "AAAAAAAAAAAAAAAAAAAA",
                                  "TTTTTTTTTTTTTTTTTTTT")
        assert amps == []


class TestSimulatePcrInputValidation:
    """Input sanitation — bad inputs return [] rather than crash."""

    def test_empty_primers(self):
        assert sc._simulate_pcr("ATGC" * 100, "", "GCATGCATGCATGCAT") == []
        assert sc._simulate_pcr("ATGC" * 100, "GCATGCATGCATGCAT", "") == []

    def test_empty_template(self):
        assert sc._simulate_pcr("", "GCATGCATGCATGCAT",
                                  "GCATGCATGCATGCAT") == []

    def test_primer_too_short(self):
        # Floor is _PCR_MIN_PRIMER_LEN
        short = "A" * (sc._PCR_MIN_PRIMER_LEN - 1)
        assert sc._simulate_pcr("ATGC" * 100, short,
                                  "GCATGCATGCATGCAT") == []

    def test_primer_too_long(self):
        too_long = "A" * (sc._PCR_MAX_PRIMER_LEN + 1)
        assert sc._simulate_pcr("ATGC" * 100, too_long,
                                  "GCATGCATGCATGCAT") == []

    def test_non_acgt_primer(self):
        # 2026-05-27 (audit-5 primer H1): IUPAC chars used to silently
        # return []; the GUI couldn't distinguish "primer was
        # filtered" from "primer doesn't bind". Now raises ValueError
        # so the caller can surface a proper error to the user.
        import pytest
        with pytest.raises(ValueError, match="IUPAC"):
            sc._simulate_pcr("ATGC" * 100, "NNNNNNNNNNNNNNNNNNNN",
                              "GCATGCATGCATGCAT")
        with pytest.raises(ValueError, match="IUPAC"):
            sc._simulate_pcr("ATGC" * 100, "GCATGCATGCATGCAT",
                              "RYWSMKBDHVRYWSMKBDHV")

    def test_none_inputs(self):
        # Defensive: type guard
        assert sc._simulate_pcr(None, "GCATGCATGCATGCAT",         # type: ignore
                                  "GCATGCATGCATGCAT") == []
        assert sc._simulate_pcr("ATGC" * 100, None,                # type: ignore
                                  "GCATGCATGCATGCAT") == []

    def test_zero_or_negative_max_amplicon(self):
        seq = "ATGCGATCGATCGATCGCGT" + "A" * 50 + sc._rc("ATGCGATCGATCGATCGCGT")
        # max_amplicon clamped to min 1; below min_amp returns []
        assert sc._simulate_pcr(seq, "ATGCGATCGATCGATCGCGT",
                                  "ATGCGATCGATCGATCGCGT",
                                  max_amplicon=0) == []

    def test_template_size_cap(self):
        # Templates larger than _PCR_MAX_TEMPLATE_BP refused to avoid
        # freezing the UI on chromosome-scale find()s.
        big = "A" * (sc._PCR_MAX_TEMPLATE_BP + 1)
        assert sc._simulate_pcr(big, "ATGCGATCGATCGATCGCGT",
                                  "ATGCGATCGATCGATCGCGT") == []

    def test_amplicon_count_cap(self):
        # Construct a template with > _PCR_MAX_AMPLICONS legit pairings —
        # confirm the result list is capped.
        # 60 copies of a fwd binding site, all paired with one rev site
        # downstream → up to 60 amplicons. Cap is 50.
        fwd_site = "GATCGATCGATCGATCGATC"
        rev_site_on_top = "GTACGTACGTACGTACGTAC"
        spacer = "A" * 10
        seq = (fwd_site + spacer) * 60 + rev_site_on_top + "A" * 50
        amps = sc._simulate_pcr(seq, fwd_site, sc._rc(rev_site_on_top),
                                  max_amplicon=5000)
        assert len(amps) <= sc._PCR_MAX_AMPLICONS

    def test_amplicon_length_below_min_excluded(self):
        # min_amp = len(fwd) + len(rev). Anything shorter excluded.
        seq = "ATGCGATCGATCGATCGCGT" + "GCATGCATGCATGCATGCAT"
        rev = sc._rc("GCATGCATGCATGCATGCAT")
        # Place rev at position 20 directly after fwd → length = 40
        # length should equal len(fwd) + len(rev) = 40
        amps = sc._simulate_pcr(seq, "ATGCGATCGATCGATCGCGT", rev,
                                  max_amplicon=5000)
        assert len(amps) == 1
        assert amps[0]["length"] == 40

    def test_primer_hit_explosion_refused(self):
        # Regression guard for 2026-05-17 hardening: a primer that
        # matches the template thousands of times (the pathological
        # "all-A primer on all-A tract" case) used to push the inner
        # O(N²) loop into multi-second territory. Now we cap each
        # side's hit list at `_PCR_MAX_PRIMER_HITS` and refuse with
        # an empty result if exceeded.
        fwd = "A" * 12   # well below _PCR_MAX_PRIMER_LEN
        seq = "A" * (sc._PCR_MAX_PRIMER_HITS + 100)  # > cap hits
        rev = "T" * 12   # _rc = "A" * 12 — same blowup
        amps = sc._simulate_pcr(seq, fwd, rev, max_amplicon=5000)
        assert amps == []

    def test_primer_hit_at_cap_still_runs(self):
        # Just BELOW the cap should not refuse — verifies the cap is a
        # ceiling, not a floor.
        # A 12-bp ACGT primer has 4^12 ≈ 16.7M possible sequences;
        # a random ATGC template has expected hit count ≈ N/4^12.
        # Build a template with exactly 10 fwd-hits + 10 rev-hits.
        fwd = "ATGCATGCATGC"   # 12 bp
        rev = "GCATGCATGCAT"   # 12 bp; _rc = "ATGCATGCATGC" again (rev-comp of palindromic-feeling string)
        # Spaced 100 bp apart → 10 fwd, then 1 rev far downstream.
        spacer = "G" * 100
        rev_rc_target = sc._rc(rev)
        seq = (fwd + spacer) * 10 + rev_rc_target + "G" * 50
        amps = sc._simulate_pcr(seq, fwd, rev, max_amplicon=5000)
        # 10 fwd hits + 1 rev_rc hit = up to 10 amplicons; capped at 50
        # which is well above. Result must NOT be the empty refusal list.
        assert len(amps) >= 1


# ═══════════════════════════════════════════════════════════════════════════════
# _agarose_mobility — empirical migration model
# ═══════════════════════════════════════════════════════════════════════════════

class TestAgaroseMobility:
    """Helling-Goodman-Boyer log-linear migration with form corrections."""

    def test_within_window_monotone(self):
        # On a 1% gel, mobility decreases as bp grows (small fragments
        # run further toward the dye front).
        mobs = [sc._agarose_mobility(bp, 1.0) for bp in
                 (500, 1000, 2000, 5000, 10000)]
        for a, b in zip(mobs, mobs[1:]):
            assert a >= b - 1e-9   # monotone non-increasing

    def test_below_lower_resolution(self):
        # Very small fragments run with the dye front
        mob = sc._agarose_mobility(50, 1.0)
        assert mob > 0.9

    def test_above_upper_resolution(self):
        # Very large fragments stick near the well
        mob = sc._agarose_mobility(100_000, 1.0)
        assert mob < 0.1

    def test_zero_or_negative_bp(self):
        assert sc._agarose_mobility(0, 1.0) == 1.0
        assert sc._agarose_mobility(-5, 1.0) == 1.0
        assert sc._agarose_mobility(None, 1.0) == 1.0  # type: ignore

    def test_supercoiled_faster_than_linear(self):
        # Same bp → supercoiled mobility is HIGHER (closer to dye front)
        bp = 3000
        m_linear = sc._agarose_mobility(bp, 1.0, "linear")
        m_sc     = sc._agarose_mobility(bp, 1.0, "supercoiled")
        assert m_sc > m_linear

    def test_nicked_slower_than_linear(self):
        bp = 3000
        m_linear = sc._agarose_mobility(bp, 1.0, "linear")
        m_nick   = sc._agarose_mobility(bp, 1.0, "nicked")
        assert m_nick < m_linear

    def test_unknown_form_defaults_to_linear(self):
        bp = 3000
        m_linear = sc._agarose_mobility(bp, 1.0, "linear")
        m_unknown = sc._agarose_mobility(bp, 1.0, "no-such-form")
        assert m_linear == pytest.approx(m_unknown)

    def test_gel_pct_snaps_to_nearest(self):
        # 0.95% should snap to 1.0% (nearest configured)
        bp = 3000
        m_exact = sc._agarose_mobility(bp, 1.0)
        m_snap  = sc._agarose_mobility(bp, 0.95)
        assert m_exact == pytest.approx(m_snap)

    def test_below_window_still_orders_by_size(self):
        """Refactor 2026-05-19 — pre-fix the boundary hard-clamped
        to 0.97 so two below-window fragments stacked on the same
        row regardless of size. The soft-asymptote now keeps them
        ordered (smaller faster) while still piling near the dye
        front."""
        # 1% gel: window is 500..10000 bp. 50/100/200 are all below.
        m_50  = sc._agarose_mobility(50,  1.0)
        m_100 = sc._agarose_mobility(100, 1.0)
        m_200 = sc._agarose_mobility(200, 1.0)
        # All near the dye front
        assert all(m > 0.9 for m in (m_50, m_100, m_200))
        # Strict monotone — smaller bp runs further
        assert m_50 > m_100 > m_200

    def test_above_window_still_orders_by_size(self):
        """Two above-window fragments must not collapse to the same
        row either — larger sticks closer to the well."""
        # 1% gel: 20000, 50000, 100000 all above window
        m_20k  = sc._agarose_mobility(20_000,  1.0)
        m_50k  = sc._agarose_mobility(50_000,  1.0)
        m_100k = sc._agarose_mobility(100_000, 1.0)
        assert all(m < 0.1 for m in (m_20k, m_50k, m_100k))
        # Strict monotone — larger bp stays closer to the well
        assert m_20k > m_50k > m_100k

    def test_in_window_bounds_unchanged(self):
        """Sanity: the in-window branch is unchanged by the
        extrapolation refactor. A 1 kb band on a 1% gel still
        returns the same Helling-Goodman-Boyer linear result it
        always did."""
        import math
        m = sc._agarose_mobility(1000, 1.0)
        # raw = (log10(10000) - log10(1000)) / (log10(10000) - log10(500))
        log_lo = math.log10(500)
        log_hi = math.log10(10000)
        expected = (log_hi - math.log10(1000)) / (log_hi - log_lo)
        assert m == pytest.approx(expected, abs=1e-6)


# ═══════════════════════════════════════════════════════════════════════════════
# _gel_bands_for_lane — source → bands resolution
# ═══════════════════════════════════════════════════════════════════════════════

class TestGelBandsForLane:
    """Each lane source kind produces the expected list of (bp, form) pairs."""

    def test_empty_returns_no_bands(self):
        bands = sc._gel_bands_for_lane(
            {"source": "empty", "detail": ""},
            template_seq="", template_circular=False, pcr_amplicon=None,
        )
        assert bands == []

    def test_ladder_picks_named_ladder(self):
        bands = sc._gel_bands_for_lane(
            {"source": "ladder", "detail": "1 kb"},
            template_seq="", template_circular=False, pcr_amplicon=None,
        )
        bps = [bp for bp, _ in bands]
        assert 1000 in bps
        assert 250 in bps  # smallest band in NEB 1 kb

    def test_ladder_unknown_name_falls_back(self):
        bands = sc._gel_bands_for_lane(
            {"source": "ladder", "detail": "not-a-ladder"},
            template_seq="", template_circular=False, pcr_amplicon=None,
        )
        assert len(bands) > 0   # uses first ladder, not empty

    def test_plasmid_circular_yields_two_bands(self):
        # Circular uncut → supercoiled + nicked (linear-from-prep-nicking
        # is not modeled; user runs a digest lane to get the linear band)
        bands = sc._gel_bands_for_lane(
            {"source": "plasmid", "detail": ""},
            template_seq="ATGC" * 100, template_circular=True,
            pcr_amplicon=None,
        )
        forms = sorted(form for _, form in bands)
        assert forms == ["nicked", "supercoiled"]
        assert all(bp == 400 for bp, _ in bands)

    def test_plasmid_linear_yields_one_band(self):
        bands = sc._gel_bands_for_lane(
            {"source": "plasmid", "detail": ""},
            template_seq="ATGC" * 100, template_circular=False,
            pcr_amplicon=None,
        )
        assert bands == [(400, "linear")]

    def test_plasmid_empty_template_yields_no_bands(self):
        bands = sc._gel_bands_for_lane(
            {"source": "plasmid", "detail": ""},
            template_seq="", template_circular=True, pcr_amplicon=None,
        )
        assert bands == []

    def test_digest_with_known_enzyme(self):
        # A circular plasmid with two EcoRI (GAATTC) sites → 2 fragments
        seq = "G" * 100 + "GAATTC" + "A" * 200 + "GAATTC" + "C" * 100
        bands = sc._gel_bands_for_lane(
            {"source": "digest", "detail": "EcoRI"},
            template_seq=seq, template_circular=True, pcr_amplicon=None,
        )
        assert len(bands) == 2
        # Both linear forms
        assert all(f == "linear" for _, f in bands)
        total = sum(bp for bp, _ in bands)
        assert total == len(seq)

    def test_digest_empty_enzyme_list(self):
        bands = sc._gel_bands_for_lane(
            {"source": "digest", "detail": ""},
            template_seq="ATGC" * 100, template_circular=True,
            pcr_amplicon=None,
        )
        assert bands == []

    def test_digest_unknown_enzyme(self):
        # _digest_with_enzymes returns single uncut frag on unknown enzymes —
        # but the resulting lane is still well-defined (one band = template).
        # Defensive: should not crash.
        bands = sc._gel_bands_for_lane(
            {"source": "digest", "detail": "FakeEnzyme,AlsoFake"},
            template_seq="ATGC" * 100, template_circular=True,
            pcr_amplicon=None,
        )
        # No matches → single uncut fragment
        assert len(bands) == 1

    def test_pcr_amplicon_source(self):
        amp = {"length": 1234, "amplicon_seq": "ATGC" * 308 + "AT"}
        bands = sc._gel_bands_for_lane(
            {"source": "pcr", "detail": ""},
            template_seq="", template_circular=False, pcr_amplicon=amp,
        )
        assert bands == [(1234, "linear")]

    def test_pcr_no_amplicon_yields_no_bands(self):
        bands = sc._gel_bands_for_lane(
            {"source": "pcr", "detail": ""},
            template_seq="", template_circular=False, pcr_amplicon=None,
        )
        assert bands == []

    def test_pcr_amplicon_with_malformed_length_yields_no_bands(self):
        # Regression guard for 2026-05-17 hardening: the agent endpoint
        # accepts an arbitrary `pcr_amplicon` dict (not just one built by
        # `_simulate_pcr`). A hostile / mis-typed `length` field used to
        # crash `int()` mid-render. Now defensively coerces to 0 → empty
        # band list.
        for bad_length in ("not-a-number", None, [1, 2], {"x": 1}):
            bands = sc._gel_bands_for_lane(
                {"source": "pcr", "detail": ""},
                template_seq="", template_circular=False,
                pcr_amplicon={"length": bad_length},
            )
            assert bands == [], f"failed for length={bad_length!r}"

    def test_pcr_amplicon_non_dict_yields_no_bands(self):
        # Defensive: `pcr_amplicon` must be a dict. A bare list / int /
        # string used to AttributeError on `.get("length")`. Now skips.
        for bad in ([1, 2], 42, "string", True):
            bands = sc._gel_bands_for_lane(
                {"source": "pcr", "detail": ""},
                template_seq="", template_circular=False,
                pcr_amplicon=bad,
            )
            assert bands == [], f"failed for amplicon={bad!r}"


# ═══════════════════════════════════════════════════════════════════════════════
# _append_pcr_gel_lane — "Send to Gel lane" lane-append logic
# ═══════════════════════════════════════════════════════════════════════════════

class TestAppendPcrGelLane:
    """The pure core of the "Send to Gel lane" button: append a new
    `pcr` lane that freezes its OWN amplicon size (`_pcr_bp`) so multiple
    amplicons coexist on one gel. Always appends to the right — never
    reuses or overwrites an existing lane."""

    def test_appends_frozen_lane(self):
        lanes = [{"name": "Ladder", "source": "ladder", "detail": "1 kb"}]
        idx, at_cap = sc._append_pcr_gel_lane(lanes, "PCR 1,234 bp", 1234, 8)
        assert at_cap is False
        assert idx == 1
        assert lanes[idx] == {"name": "PCR 1,234 bp", "source": "pcr",
                               "detail": "", "_pcr_bp": 1234}

    def test_each_call_appends_distinct_lane(self):
        # Two sends → two pcr lanes, each pinned to its own size.
        lanes = [{"name": "Ladder", "source": "ladder", "detail": "1 kb"}]
        sc._append_pcr_gel_lane(lanes, "PCR 500 bp", 500, 8)
        sc._append_pcr_gel_lane(lanes, "PCR 900 bp", 900, 8)
        pcr = [ln for ln in lanes if ln["source"] == "pcr"]
        assert len(pcr) == 2
        assert [ln["_pcr_bp"] for ln in pcr] == [500, 900]

    def test_never_overwrites_existing_lane(self):
        # A user-added lane between sends is preserved; the new pcr lane
        # lands to its right.
        lanes = [
            {"name": "Ladder", "source": "ladder", "detail": "1 kb"},
            {"name": "PCR 500 bp", "source": "pcr", "detail": "",
             "_pcr_bp": 500},
            {"name": "my digest", "source": "digest", "detail": "EcoRI"},
        ]
        idx, at_cap = sc._append_pcr_gel_lane(lanes, "PCR 900 bp", 900, 8)
        assert idx == 3
        assert lanes[2]["name"] == "my digest"      # untouched
        assert lanes[3]["_pcr_bp"] == 900

    def test_at_cap_refuses(self):
        lanes = [{"name": f"L{i}", "source": "ladder", "detail": ""}
                 for i in range(8)]
        idx, at_cap = sc._append_pcr_gel_lane(lanes, "PCR 100 bp", 100, 8)
        assert (idx, at_cap) == (-1, True)
        assert len(lanes) == 8                       # nothing appended

    def test_frozen_size_renders_independent_of_selected_amplicon(self):
        # A lane's frozen `_pcr_bp` wins over whatever amplicon is
        # currently selected — that's what lets two sends show two
        # different bands on one gel.
        lanes = [{"name": "Ladder", "source": "ladder", "detail": "1 kb"}]
        idx, _ = sc._append_pcr_gel_lane(lanes, "PCR 777 bp", 777, 8)
        bands = sc._gel_bands_for_lane(
            lanes[idx],
            template_seq="", template_circular=False,
            pcr_amplicon={"length": 1234},   # different selected amplicon
        )
        assert bands == [(777, "linear")]

    def test_manual_pcr_lane_falls_back_to_selected(self):
        # A pcr lane WITHOUT a frozen size (added via the source
        # dropdown, not the button) still tracks the selected amplicon.
        lane = {"name": "PCR", "source": "pcr", "detail": ""}
        bands = sc._gel_bands_for_lane(
            lane, template_seq="", template_circular=False,
            pcr_amplicon={"length": 333},
        )
        assert bands == [(333, "linear")]


# ═══════════════════════════════════════════════════════════════════════════════
# _render_gel_image — visual rendering smoke
# ═══════════════════════════════════════════════════════════════════════════════

class TestRenderGelImage:
    """Renders without crashing; returns a Rich Text with expected anchors."""

    def test_empty_lanes_returns_notice(self):
        from rich.text import Text
        rt = sc._render_gel_image(
            [], template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0,
        )
        assert isinstance(rt, Text)
        assert "no lanes" in str(rt)

    def test_ladder_lane_emits_bp_labels(self):
        lanes = [{"name": "L", "source": "ladder", "detail": "1 kb"}]
        rt = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0,
        )
        text = str(rt)
        # NEB 1 kb has a 10000 band → "10.0k" label at the top of the gel
        assert "10.0k" in text or "10000" in text

    def test_dye_front_present(self):
        # Bottom-of-gel separator (░) renders unconditionally
        lanes = [{"name": "L", "source": "ladder", "detail": "1 kb"}]
        rt = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0,
        )
        assert "░" in str(rt)

    def test_multi_source_lanes(self):
        seq = "G" * 100 + "GAATTC" + "A" * 200 + "GAATTC" + "C" * 100
        lanes = [
            {"name": "L",      "source": "ladder",  "detail": "1 kb"},
            {"name": "Uncut",  "source": "plasmid", "detail": ""},
            {"name": "EcoRI",  "source": "digest",  "detail": "EcoRI"},
        ]
        rt = sc._render_gel_image(
            lanes, template_seq=seq, template_circular=True,
            pcr_amplicon=None, agarose_pct=1.0,
        )
        text = str(rt)
        # Three lane labels appear in the header
        assert "Uncut" in text
        assert "EcoRI" in text

    def test_height_param_respected(self):
        lanes = [{"name": "L", "source": "ladder", "detail": "1 kb"}]
        rt_short = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0, height=10,
        )
        rt_tall = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0, height=30,
        )
        # Taller gel has more lines
        assert str(rt_tall).count("\n") > str(rt_short).count("\n")

    def test_sub_row_faint_band_appears(self):
        """Refactor 2026-05-19 — bands whose fractional row offset
        exceeds 0.25 render a LIGHT glyph (`─`) on the adjacent row
        so the eye reads them as "between rows". Lets two bands
        separated by less than a full row of mobility still resolve
        visually instead of collapsing on rounding.

        Construction: pick a band whose mobility lands at a row
        boundary mid-cell. A 1 kb ladder on a 1% gel at height=22
        always produces at least one sub-row offset because the
        log10-spacing of ladder rungs doesn't align perfectly with
        whole rows."""
        lanes = [{"name": "L", "source": "ladder", "detail": "1 kb"}]
        rt = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0, height=22,
        )
        text = str(rt)
        # The faint glyph `─` (U+2500) must appear at least once in
        # the body rows — most ladder bands will produce a faint
        # tail on a 22-row render.
        assert "─" in text, (
            "expected at least one faint `─` glyph from a sub-row "
            f"band tail; got:\n{text}"
        )

    def test_band_glyph_is_thin_line(self):
        """User UX call 2026-05-19 — bands use the heavy horizontal
        line glyph `━` (U+2501) rather than a full block. Solid
        blocks read as a wall; thin lines read as proper gel bands.
        Wells stay as `█` so the visual distinction between well
        (where the DNA was loaded) and band (where it migrated to)
        remains clear."""
        lanes = [{"name": "L", "source": "ladder", "detail": "1 kb"}]
        rt = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0,
        )
        text = str(rt)
        # Single-band cells use `━`; wells use `█`.
        assert "━" in text
        assert "█" in text

    def test_bands_align_to_well_columns(self):
        """Each band must occupy the exact same column-range as its
        lane's well row. Pre-fix (alignment bug 2026-05-19) the body
        rows with a ladder bp label drifted one column left of
        unlabelled rows because `f"{label} "` was 6 chars while
        `label_col` was 7 — the labelled row's bands didn't align
        with the wells. The `ljust(label_col)` fix made every row's
        lane columns identical regardless of whether the row has a
        ladder bp tick."""
        lanes = [{"name": "L", "source": "ladder", "detail": "1 kb"}]
        rt = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0,
            height=20, lane_width=7, label_col=7,
        )
        lines = str(rt).splitlines()
        # Wells row uses `█`. Find it (must start with whitespace
        # then `█`, i.e. the lane-only well row, not a labelled
        # band row).
        well_row = next(
            ln for ln in lines
            if "█" in ln and ln.lstrip().startswith("█")
        )
        well_start = well_row.index("█")
        # Body band rows use `━` (single band) or `▆` / `█` (multi-
        # band) — find ones that contain `━`.
        body_rows = [ln for ln in lines if "━" in ln]
        assert body_rows, (
            "expected at least one body band row with `━` from the "
            "ladder lane"
        )
        for br in body_rows:
            assert br.index("━") == well_start, (
                f"band column {br.index('━')} != well column "
                f"{well_start}: row={br!r}"
            )

    def test_pcr_lane_with_no_amplicon_shows_hint(self):
        """Regression guard for 2026-05-17 audit fix: when a lane uses
        ``source='pcr'`` but ``pcr_amplicon`` is None, the render appends
        a dim-italic caption naming the affected lane(s). Pre-fix the
        empty PCR lane was indistinguishable from a digest that failed
        or a plasmid lane the user forgot to configure."""
        lanes = [
            {"name": "L",   "source": "ladder", "detail": "1 kb"},
            {"name": "Amp", "source": "pcr",    "detail": ""},
        ]
        rt = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0,
        )
        text = str(rt)
        assert "no amplicon" in text.lower()
        # Lane index (1-based) included so the user knows WHICH lane.
        assert "lane 2" in text.lower()

    def test_pcr_lane_with_amplicon_omits_hint(self):
        """Negative case for the 2026-05-17 hint: a populated PCR lane
        must NOT emit the 'no amplicon' caption."""
        lanes = [{"name": "Amp", "source": "pcr", "detail": ""}]
        rt = sc._render_gel_image(
            lanes, template_seq="ATGC" * 100, template_circular=False,
            pcr_amplicon={"length": 400, "amplicon_seq": "ATGC" * 100},
            agarose_pct=1.0,
        )
        assert "no amplicon" not in str(rt).lower()

    def test_multiple_empty_pcr_lanes_listed(self):
        """Two PCR-empty lanes → both numbers appear, with plural
        'lanes'."""
        lanes = [
            {"name": "A", "source": "pcr", "detail": ""},
            {"name": "L", "source": "ladder", "detail": "1 kb"},
            {"name": "B", "source": "pcr", "detail": ""},
        ]
        rt = sc._render_gel_image(
            lanes, template_seq="", template_circular=False,
            pcr_amplicon=None, agarose_pct=1.0,
        )
        text = str(rt).lower()
        assert "lanes 1, 3" in text
        assert "no amplicon" in text


# ═══════════════════════════════════════════════════════════════════════════════
# SimulatorScreen — smoke construction
# ═══════════════════════════════════════════════════════════════════════════════

class TestSimulatorScreenConstruction:
    """Screen builds with various input shapes without exceptions."""

    def test_with_full_template(self):
        s = sc.SimulatorScreen("ATGC" * 100, [
            {"type": "CDS", "start": 0, "end": 30, "strand": 1,
             "label": "lacZ", "color": "white"},
        ], "pUC19", "circular")
        assert s._template == "ATGC" * 100
        assert s._template_circular is True
        assert s._plasmid_name == "pUC19"
        # Default lane config seeded with 4 lanes
        assert len(s._lanes) == 4

    def test_with_empty_template(self):
        s = sc.SimulatorScreen()
        assert s._template == ""
        assert s._template_circular is True   # default
        assert s._plasmid_name == "(no plasmid)"

    def test_linear_topology(self):
        s = sc.SimulatorScreen("ATGC" * 100, [], "test", "linear")
        assert s._template_circular is False

    def test_none_args_handled_defensively(self):
        s = sc.SimulatorScreen(None, None, None, None)   # type: ignore
        assert s._template == ""
        assert s._plasmid_name == "(no plasmid)"
        assert s._template_circular is True


# ═══════════════════════════════════════════════════════════════════════════════
# Send to Gel lane — full handler wiring (pilot integration)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSendToGelLaneWiring:
    """End-to-end: "Send to Gel lane" clears the demo lanes on the first
    send (ladder stays in lane 1), appends a fresh `pcr` lane per send,
    switches to the Gel tab, and renders the gel. Pre-fix the handler
    only switched tabs + toasted; the amplicon never reached a lane."""

    def _host_app(self):
        from textual.app import App

        class _Host(App):
            def on_mount(self) -> None:
                self.push_screen(
                    sc.SimulatorScreen("ATGC" * 100, [], "pUC19", "circular")
                )
        return _Host()

    def _seed_amplicon(self, screen, length):
        screen._pcr_amplicons = [{
            "length": length, "amplicon_seq": "A" * length,
            "start": 0, "end": length, "wraps": False,
            "gc_pct": 50.0, "fwd_tm": 60.0, "rev_tm": 60.0,
        }]
        screen._selected_pcr_idx = 0

    async def test_first_send_clears_demo_lanes_ladder_first(self):
        from textual.widgets import TabbedContent, Static
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            assert isinstance(screen, sc.SimulatorScreen)
            # Starts with the 4 demo lanes.
            assert len(screen._lanes) == 4
            self._seed_amplicon(screen, 1234)
            screen._on_pcr_send_to_gel(None)
            await pilot.pause()
            await pilot.pause()
            # Demo lanes cleared → ladder + one pcr lane.
            assert len(screen._lanes) == 2
            assert screen._lanes[0]["source"] == "ladder"
            assert screen._lanes[1]["source"] == "pcr"
            assert screen._lanes[1]["_pcr_bp"] == 1234
            assert "1,234" in screen._lanes[1]["name"]
            # Switched to the Gel tab + auto-rendered (the Static now
            # holds the gel image — non-empty).
            tabs = screen.query_one("#sim-tabs", TabbedContent)
            assert tabs.active == "sim-tab-gel"
            img = screen.query_one("#sim-gel-image", Static)
            assert str(img.render()).strip() != ""

    async def test_subsequent_sends_append_and_preserve_user_lanes(self):
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            assert isinstance(screen, sc.SimulatorScreen)
            # First send: 500 bp → [Ladder, PCR 500].
            self._seed_amplicon(screen, 500)
            screen._on_pcr_send_to_gel(None)
            await pilot.pause()
            assert [ln["source"] for ln in screen._lanes] == ["ladder", "pcr"]
            # User adds a lane in between.
            screen._lanes.append(
                {"name": "my check", "source": "digest", "detail": "EcoRI"}
            )
            screen._refresh_lane_rows()
            await pilot.pause()
            # Second send: 900 bp appends to the RIGHT, user lane intact.
            self._seed_amplicon(screen, 900)
            screen._on_pcr_send_to_gel(None)
            await pilot.pause()
            srcs = [ln["source"] for ln in screen._lanes]
            assert srcs == ["ladder", "pcr", "digest", "pcr"]
            assert screen._lanes[2]["name"] == "my check"   # preserved
            # The two pcr lanes carry their own frozen sizes.
            pcr_bps = [ln["_pcr_bp"] for ln in screen._lanes
                       if ln["source"] == "pcr"]
            assert pcr_bps == [500, 900]

    async def test_button_warns_when_no_amplicon(self):
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            assert isinstance(screen, sc.SimulatorScreen)
            screen._pcr_amplicons = []
            screen._selected_pcr_idx = -1
            n_lanes_before = len(screen._lanes)
            screen._on_pcr_send_to_gel(None)
            await pilot.pause()
            # No amplicon → no lane mutation, demo lanes untouched.
            assert len(screen._lanes) == n_lanes_before
            assert screen._gel_is_demo is True


# ═══════════════════════════════════════════════════════════════════════════════
# Gel name field + demo-vs-user-gel state
# ═══════════════════════════════════════════════════════════════════════════════

class TestGelNameAndDemoState:
    """A freshly-opened gel is the auto-named "Demo Gel" until the user
    adds / removes / edits a lane (or sends an amplicon), after which it
    becomes their own gel and the auto-name is dropped."""

    def _host_app(self):
        from textual.app import App

        class _Host(App):
            def on_mount(self) -> None:
                self.push_screen(
                    sc.SimulatorScreen("ATGC" * 100, [], "pUC19", "circular")
                )
        return _Host()

    def test_fresh_gel_is_demo(self):
        s = sc.SimulatorScreen("ATGC" * 100, [], "p", "circular")
        assert s._gel_is_demo is True
        assert s._gel_name == "Demo Gel"

    async def test_name_field_shows_demo_gel_on_open(self):
        from textual.widgets import Input
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            name_w = screen.query_one("#sim-gel-name", Input)
            assert name_w.value == "Demo Gel"

    async def test_add_lane_exits_demo_and_renames(self):
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            screen._on_gel_add_lane(None)
            await pilot.pause()
            assert screen._gel_is_demo is False
            assert screen._gel_name == "Untitled Gel"

    async def test_send_to_gel_exits_demo_and_renames(self):
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            screen._pcr_amplicons = [{
                "length": 800, "amplicon_seq": "A" * 800,
                "start": 0, "end": 800, "wraps": False,
                "gc_pct": 50.0, "fwd_tm": 60.0, "rev_tm": 60.0,
            }]
            screen._selected_pcr_idx = 0
            screen._on_pcr_send_to_gel(None)
            await pilot.pause()
            assert screen._gel_is_demo is False
            assert screen._gel_name == "Untitled Gel"

    async def test_typing_gel_name_exits_demo(self):
        from textual.widgets import Input
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            name_w = screen.query_one("#sim-gel-name", Input)
            name_w.value = "My favourite gel"
            await pilot.pause()
            assert screen._gel_name == "My favourite gel"
            assert screen._gel_is_demo is False

    async def test_mount_does_not_exit_demo(self):
        # The mount-time Changed echoes from the demo lane widgets must
        # NOT flip the gel out of demo state (readiness-flag guard).
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            await pilot.pause()
            screen = app.screen
            assert screen._gel_is_demo is True
            assert screen._gel_name == "Demo Gel"
            # Demo ladder detail not auto-reset on open.
            assert screen._lanes[0]["detail"] == "1 kb"


# ═══════════════════════════════════════════════════════════════════════════════
# Library entry build (save-to-library path)
# ═══════════════════════════════════════════════════════════════════════════════

class TestBuildAmpliconLibraryEntry:
    """Save-to-library round-trip: amplicon dict → library-shaped entry."""

    def _make_screen(self):
        # __new__ instead of __init__ so we don't construct widgets
        # (avoids needing an active App context)
        s = sc.SimulatorScreen.__new__(sc.SimulatorScreen)
        s._plasmid_name = "pUC19"
        s._template = "ATGC" * 100
        s._template_circular = True
        return s

    def test_schema_fields(self):
        s = self._make_screen()
        amp = {
            "start": 0, "end": 50, "length": 50,
            "wraps": False, "fwd_seq": "ATGCGATCGATCGATCGCGT",
            "rev_seq": "ATGCGATCGATCGATCGCGT",
            "amplicon_seq": "ATGCGATCGATCGATCGCGT" + "A" * 30,
            "gc_pct": 50.0, "fwd_tm": 60.0, "rev_tm": 60.0,
        }
        entry = s._build_amplicon_library_entry(amp)
        assert set(entry.keys()) >= {"id", "name", "size", "n_feats",
                                       "source", "added", "gb_text"}
        assert entry["size"] == 50
        assert entry["source"] == "simulator:pcr"
        # Has primer features
        assert entry["n_feats"] >= 1
        # GenBank text declares linear topology
        assert "linear" in entry["gb_text"]

    def test_empty_amplicon_rejected(self):
        s = self._make_screen()
        amp = {
            "start": 0, "end": 0, "length": 0, "wraps": False,
            "fwd_seq": "AAAA", "rev_seq": "AAAA",
            "amplicon_seq": "",
            "gc_pct": 0.0, "fwd_tm": None, "rev_tm": None,
        }
        with pytest.raises(ValueError):
            s._build_amplicon_library_entry(amp)

    def test_name_sanitized(self):
        s = self._make_screen()
        s._plasmid_name = "weird name!@#$%^&*()"
        amp = {
            "start": 0, "end": 50, "length": 50, "wraps": False,
            "fwd_seq": "ATGCGATCGATCGATCGCGT",
            "rev_seq": "ATGCGATCGATCGATCGCGT",
            "amplicon_seq": "A" * 50,
            "gc_pct": 0.0, "fwd_tm": None, "rev_tm": None,
        }
        entry = s._build_amplicon_library_entry(amp)
        # ID has only A-Z0-9_- chars
        import re
        assert re.match(r"^[A-Za-z0-9_-]+$", entry["id"])

    def test_explicit_name_used_verbatim(self):
        # The collection-targeted save passes the user's name; it must
        # be used as-is (the commit path re-uniquifies against the
        # target collection, NOT the active library).
        s = self._make_screen()
        amp = {
            "start": 0, "end": 50, "length": 50, "wraps": False,
            "fwd_seq": "ATGCGATCGATCGATCGCGT",
            "rev_seq": "ATGCGATCGATCGATCGCGT",
            "amplicon_seq": "A" * 50,
            "gc_pct": 0.0, "fwd_tm": None, "rev_tm": None,
        }
        entry = s._build_amplicon_library_entry(amp, name="My Insert v2")
        assert entry["name"] == "My Insert v2"
        # id derived from the explicit name.
        assert entry["id"] == "My_Insert_v2"


# ═══════════════════════════════════════════════════════════════════════════════
# AmpliconSaveModal — name + collection picker (pure init logic)
# ═══════════════════════════════════════════════════════════════════════════════

class TestAmpliconSaveModal:
    """The modal's constructor normalises the collection list so the
    Select always mounts with a valid value."""

    def test_active_collection_first_and_deduped(self):
        m = sc.AmpliconSaveModal(
            default_name="amp", collections=["A", "B", "A", "Eden"],
            active_collection="B",
        )
        # Active first, duplicates removed, order otherwise preserved.
        assert m._collections == ["B", "A", "Eden"]
        assert m._active == "B"

    def test_empty_collections_falls_back_to_default(self):
        m = sc.AmpliconSaveModal(
            default_name="amp", collections=[], active_collection=None,
        )
        assert m._collections == ["Default"]
        assert m._active == "Default"

    def test_blank_default_name_falls_back(self):
        m = sc.AmpliconSaveModal(
            default_name="   ", collections=["X"], active_collection="X",
        )
        assert m._default_name == "PCR amplicon"


# ═══════════════════════════════════════════════════════════════════════════════
# _commit_amplicon_to_collection — collection-targeted save
# ═══════════════════════════════════════════════════════════════════════════════

def _seed_collections(active="Default", others=None, seed_entries=None):
    """Seed collections.json + the active-library mirror. `seed_entries`
    pre-loads the active collection (for collision tests)."""
    colls = [{"name": active, "description": "",
              "plasmids": list(seed_entries or [])}]
    for name in (others or []):
        colls.append({"name": name, "description": "", "plasmids": []})
    sc._save_collections(colls)
    sc._set_active_collection_name(active)
    sc._settings_flush_sync()
    sc._safe_save_json_mirror(
        sc._LIBRARY_FILE, list(seed_entries or []), "Plasmid library")
    sc._library_cache = None


def _amp_entry(name="amp1", entry_id="amp1", size=500):
    return {"id": entry_id, "name": name, "size": size, "n_feats": 1,
            "source": "simulator:pcr", "added": "2026-05-27",
            "gb_text": f"LOCUS {entry_id} {size} bp\n"}


def _lib_entry_with_seq(name, eid, seq="ATGC" * 50, circular=True):
    """A library entry carrying a real, parseable GenBank record (so the
    template picker can load its sequence)."""
    from io import StringIO
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(
        Seq(seq), id=eid[:16], name=eid[:16],
        annotations={"molecule_type": "DNA",
                      "topology": "circular" if circular else "linear"},
    )
    buf = StringIO()
    SeqIO.write([rec], buf, "genbank")
    return {"id": eid, "name": name, "size": len(seq), "n_feats": 0,
            "source": "test", "added": "2026-05-28",
            "gb_text": buf.getvalue()}


class TestCommitAmpliconToCollection:
    """`_commit_amplicon_to_collection` appends a linear amplicon into a
    chosen collection with name/id collision-rename + active-mirror
    sync. Tests bypass the modal to focus on the transactional logic."""

    def _screen(self):
        return sc.SimulatorScreen.__new__(sc.SimulatorScreen)

    def test_save_into_active_syncs_mirror(self):
        _seed_collections(active="Default")
        s = self._screen()
        final = s._commit_amplicon_to_collection(_amp_entry(), "Default")
        assert final == "amp1"
        colls = sc._load_collections()
        default = next(c for c in colls if c["name"] == "Default")
        assert [e["name"] for e in default["plasmids"]] == ["amp1"]
        # Active-mirror updated too.
        lib = sc._load_library()
        assert any(e.get("id") == "amp1" for e in lib)

    def test_save_into_other_collection_leaves_active_mirror(self):
        _seed_collections(active="Default", others=["Eden"])
        s = self._screen()
        s._commit_amplicon_to_collection(_amp_entry(), "Eden")
        colls = sc._load_collections()
        eden = next(c for c in colls if c["name"] == "Eden")
        default = next(c for c in colls if c["name"] == "Default")
        assert [e["name"] for e in eden["plasmids"]] == ["amp1"]
        assert default["plasmids"] == []
        # Active (Default) mirror untouched.
        assert sc._load_library() == []

    def test_name_collision_gets_copy_suffix(self):
        _seed_collections(active="Default",
                          seed_entries=[_amp_entry("amp1", "existing_id")])
        s = self._screen()
        final = s._commit_amplicon_to_collection(
            _amp_entry("amp1", "amp1"), "Default")
        assert final == "amp1 COPY"
        names = [e["name"] for e in
                 sc._load_library()]
        assert "amp1" in names and "amp1 COPY" in names

    def test_id_collision_gets_numeric_suffix(self):
        _seed_collections(active="Default",
                          seed_entries=[_amp_entry("other", "dup_id")])
        s = self._screen()
        s._commit_amplicon_to_collection(
            _amp_entry("fresh", "dup_id"), "Default")
        ids = [e["id"] for e in sc._load_library()]
        assert "dup_id" in ids and "dup_id_2" in ids

    def test_missing_collection_is_created(self):
        _seed_collections(active="Default")
        s = self._screen()
        s._commit_amplicon_to_collection(_amp_entry(), "Brand New")
        colls = sc._load_collections()
        assert any(c["name"] == "Brand New" for c in colls)
        new = next(c for c in colls if c["name"] == "Brand New")
        assert [e["name"] for e in new["plasmids"]] == ["amp1"]

    def test_entry_deepcopied_no_caller_leak(self):
        _seed_collections(active="Default")
        s = self._screen()
        entry = _amp_entry()
        s._commit_amplicon_to_collection(entry, "Default")
        # Mutating the caller's dict after the commit must not change
        # what landed on disk.
        entry["name"] = "mutated"
        default = next(c for c in sc._load_collections()
                       if c["name"] == "Default")
        assert default["plasmids"][0]["name"] == "amp1"


# ═══════════════════════════════════════════════════════════════════════════════
# Save-amplicon flow — modal → callback → commit (pilot integration)
# ═══════════════════════════════════════════════════════════════════════════════

class TestSaveAmpliconFlow:
    """End-to-end: clicking Save opens AmpliconSaveModal; confirming it
    commits the amplicon into the chosen collection."""

    def _host_app(self):
        from textual.app import App

        class _Host(App):
            def on_mount(self) -> None:
                self.push_screen(
                    sc.SimulatorScreen("ATGC" * 100, [], "pUC19", "circular")
                )
        return _Host()

    async def test_save_button_opens_modal_then_commits(self):
        from textual.widgets import Input, Select
        _seed_collections(active="Default", others=["Eden"])
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            assert isinstance(screen, sc.SimulatorScreen)
            screen._entry_counter = getattr(app, "_record_load_counter", 0)
            screen._pcr_amplicons = [{
                "start": 0, "end": 600, "length": 600, "wraps": False,
                "fwd_seq": "ATGCGATCGATCGATCGCGT",
                "rev_seq": "ATGCGATCGATCGATCGCGT",
                "amplicon_seq": "A" * 600,
                "gc_pct": 50.0, "fwd_tm": 60.0, "rev_tm": 60.0,
            }]
            screen._selected_pcr_idx = 0
            screen._on_pcr_save(None)
            await pilot.pause()
            # The naming modal is now on top.
            modal = app.screen
            assert isinstance(modal, sc.AmpliconSaveModal)
            modal.query_one("#ampsave-name", Input).value = "My amplicon"
            modal.query_one("#ampsave-collection", Select).value = "Eden"
            modal._submit()
            await pilot.pause()
            await pilot.pause()
            # Landed in Eden, not the active Default.
            colls = sc._load_collections()
            eden = next(c for c in colls if c["name"] == "Eden")
            assert [e["name"] for e in eden["plasmids"]] == ["My amplicon"]

    async def test_save_refuses_with_no_amplicon(self):
        _seed_collections(active="Default")
        app = self._host_app()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            screen._pcr_amplicons = []
            screen._selected_pcr_idx = -1
            screen._on_pcr_save(None)
            await pilot.pause()
            # No modal pushed — still on the Simulator screen.
            assert isinstance(app.screen, sc.SimulatorScreen)


# ═══════════════════════════════════════════════════════════════════════════════
# PCR template picker — dropdown defaulting to the active plasmid
# ═══════════════════════════════════════════════════════════════════════════════

class TestPcrTemplatePicker:
    """The Template is a dropdown plasmid picker. It pre-selects the
    plasmid active in the main app when the Simulator opened; with none
    active it adopts the first library plasmid."""

    def test_active_plasmid_preselected(self):
        s = sc.SimulatorScreen("ATGC" * 100, [], "pUC19", "circular")
        assert s._current_template_id == "__loaded__"
        assert s._loaded_name == "pUC19"

    def test_no_active_picks_first_library(self):
        _seed_collections(active="Default", seed_entries=[
            _lib_entry_with_seq("plasmidA", "pa"),
            _lib_entry_with_seq("plasmidB", "pb"),
        ])
        s = sc.SimulatorScreen("", [], "", "circular")
        assert s._current_template_id == "pa"   # first library plasmid

    def test_no_active_empty_library_sentinel(self):
        _seed_collections(active="Default")     # empty library
        s = sc.SimulatorScreen("", [], "", "circular")
        assert s._current_template_id == "__none__"

    def test_build_options_current_first_then_library(self):
        _seed_collections(active="Default", seed_entries=[
            _lib_entry_with_seq("plasmidA", "pa"),
        ])
        s = sc.SimulatorScreen("ATGC" * 100, [], "pUC19", "circular")
        opts = s._build_template_options()
        assert opts[0][1] == "__loaded__"
        assert opts[0][0].startswith("pUC19")
        assert any(v == "pa" for _, v in opts)

    def test_build_options_dedups_loaded_from_library(self):
        # When the active plasmid IS a library row (same name), it's
        # offered only once — as "(current)".
        _seed_collections(active="Default", seed_entries=[
            _lib_entry_with_seq("pUC19", "puc"),
        ])
        s = sc.SimulatorScreen("ATGC" * 100, [], "pUC19", "circular")
        opts = s._build_template_options()
        assert sum(1 for label, _ in opts if "pUC19" in label) == 1

    def test_meta_text(self):
        s = sc.SimulatorScreen("ATGC" * 25, [], "p", "linear")
        assert "100 bp" in s._template_meta_text()
        assert "linear" in s._template_meta_text()


class TestPcrTemplatePickerWiring:
    """Pilot-level: the dropdown loads the picked plasmid as the PCR
    template, clearing stale results."""

    def _host(self, seq="ATGC" * 100, name="pUC19", topo="circular"):
        from textual.app import App

        class _Host(App):
            def on_mount(self) -> None:
                self.push_screen(sc.SimulatorScreen(seq, [], name, topo))
        return _Host()

    async def test_dropdown_value_is_loaded_on_open(self):
        from textual.widgets import Select
        _seed_collections(active="Default", seed_entries=[
            _lib_entry_with_seq("plasmidA", "pa"),
        ])
        app = self._host()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            sel = app.screen.query_one("#sim-pcr-template-select", Select)
            assert sel.value == "__loaded__"

    async def test_no_active_adopts_first_library_plasmid(self):
        from textual.app import App
        _seed_collections(active="Default", seed_entries=[
            _lib_entry_with_seq("plasmidA", "pa", seq="TTTTAAAA" * 20),
        ])

        class _Host(App):
            def on_mount(self) -> None:
                self.push_screen(sc.SimulatorScreen("", [], "", "circular"))
        app = _Host()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            await pilot.pause()
            screen = app.screen
            assert screen._template == "TTTTAAAA" * 20
            assert screen._plasmid_name == "plasmidA"

    async def test_pick_changes_template_and_clears_results(self):
        from textual.widgets import Select
        _seed_collections(active="Default", seed_entries=[
            _lib_entry_with_seq("plasmidA", "pa",
                                 seq="GGGGCCCC" * 15, circular=False),
        ])
        app = self._host(seq="ATGC" * 100, name="pUC19")
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            assert screen._plasmid_name == "pUC19"
            # Pretend a prior PCR run left results.
            screen._pcr_amplicons = [{"length": 50}]
            screen._selected_pcr_idx = 0
            screen.query_one("#sim-pcr-template-select", Select).value = "pa"
            await pilot.pause()
            assert screen._template == "GGGGCCCC" * 15
            assert screen._template_circular is False
            assert screen._plasmid_name == "plasmidA"
            assert screen._current_template_id == "pa"
            # Stale results dropped.
            assert screen._pcr_amplicons == []
            assert screen._selected_pcr_idx == -1

    async def test_template_and_source_share_one_row(self):
        from textual.widgets import Select, Static
        _seed_collections(active="Default", seed_entries=[
            _lib_entry_with_seq("plasmidA", "pa"),
        ])
        app = self._host(seq="ATGC" * 100, name="pUC19")
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            screen = app.screen
            tpl = screen.query_one("#sim-pcr-template-select", Select)
            src = screen.query_one("#sim-pcr-source", Select)
            meta = screen.query_one("#sim-pcr-template-meta", Static)
            # All three on the same row.
            assert tpl.region.y == src.region.y == meta.region.y
            # Source picker sits to the RIGHT of the template picker.
            assert src.region.x >= tpl.region.x + tpl.region.width
            # Template picker is widened to fit full plasmid names.
            assert tpl.region.width >= 40

    async def test_amplicon_input_packs_after_label_and_labels_centered(self):
        from textual.widgets import Input, Label
        app = self._host()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            screen = app.screen
            row = screen.query_one("#sim-pcr-params-row")
            label = row.query_one(Label)
            inp = screen.query_one("#sim-pcr-maxamp", Input)
            # Input packs immediately after the label — not shoved to
            # the far right edge as it was when the label expanded.
            assert inp.region.x <= label.region.x + label.region.width + 1
            assert inp.region.x < 40
            # Labels span the full 3-tall row so their text vertically
            # centres on the input (was top-aligned before).
            assert label.region.height == 3
            tpl_label = screen.query_one(
                "#sim-pcr-template-row").query_one(Label)
            assert tpl_label.region.height == 3
            # Default is the common 500 bp, and the box is wide enough
            # for a 7-figure entry (Textual Input chrome is 6 cols, so a
            # 7-digit value + cursor needs total width >= 14).
            assert inp.value == "500"
            assert inp.region.width >= 14


# ═══════════════════════════════════════════════════════════════════════════════
# Simulator button double-fire hardening
# ═══════════════════════════════════════════════════════════════════════════════

class TestSimulatorButtonHardening:
    """Real terminals can deliver two `Pressed` events for one physical
    click. Close must not pop the screen twice (ScreenStackError), and
    the modal-opening buttons must not stack two dialogs (double-save)."""

    def _host(self):
        from textual.app import App

        class _Host(App):
            def on_mount(self) -> None:
                self.push_screen(
                    sc.SimulatorScreen("ATGC" * 100, [], "pUC19", "circular")
                )
        return _Host()

    def _seed_amp(self, screen):
        screen._entry_counter = 0
        screen._pcr_amplicons = [{
            "start": 0, "end": 500, "length": 500, "wraps": False,
            "fwd_seq": "ATGCGATCGATCGATCGCGT",
            "rev_seq": "ATGCGATCGATCGATCGCGT",
            "amplicon_seq": "A" * 500,
            "gc_pct": 50.0, "fwd_tm": 60.0, "rev_tm": 60.0,
        }]
        screen._selected_pcr_idx = 0

    async def test_close_double_fire_no_crash(self):
        app = self._host()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            assert isinstance(screen, sc.SimulatorScreen)
            screen._close(None)
            screen._close(None)   # second must be a no-op, not a crash
            await pilot.pause()
            assert screen._dismissed is True

    async def test_save_double_fire_opens_one_modal(self):
        _seed_collections(active="Default")
        app = self._host()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            self._seed_amp(screen)
            screen._on_pcr_save(None)
            screen._on_pcr_save(None)   # double-fire
            await pilot.pause()
            modals = [s for s in app.screen_stack
                      if isinstance(s, sc.AmpliconSaveModal)]
            assert len(modals) == 1

    async def test_child_modal_flag_clears_on_dismiss(self):
        from textual.widgets import Input
        _seed_collections(active="Default")
        app = self._host()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            self._seed_amp(screen)
            screen._on_pcr_save(None)
            await pilot.pause()
            assert screen._child_modal_open is True
            modal = app.screen
            assert isinstance(modal, sc.AmpliconSaveModal)
            modal.query_one("#ampsave-name", Input).value = "amp"
            modal._submit()
            await pilot.pause()
            # Flag cleared so a later Save / Library opens normally.
            assert screen._child_modal_open is False

    async def test_library_double_fire_opens_one_modal(self):
        _seed_collections(active="Default")
        app = self._host()
        async with app.run_test(size=(170, 48)) as pilot:
            await pilot.pause()
            screen = app.screen
            # Switch to the Gel tab so the lane widgets are queryable.
            screen.query_one("#sim-tabs").active = "sim-tab-gel"
            await pilot.pause()
            screen._on_gel_library(None)
            screen._on_gel_library(None)   # double-fire
            await pilot.pause()
            libs = [s for s in app.screen_stack
                    if isinstance(s, sc.GelLibraryModal)]
            assert len(libs) == 1


# ═══════════════════════════════════════════════════════════════════════════════
# _exact_match_positions — helper
# ═══════════════════════════════════════════════════════════════════════════════

class TestExactMatchPositions:
    """Sliding-window exact match including overlapping hits."""

    def test_no_matches(self):
        assert sc._exact_match_positions("ATGC", "TTTT") == []

    def test_single_match(self):
        assert sc._exact_match_positions("ATGCATGC", "ATGC") == [0, 4]

    def test_overlapping_match(self):
        # Overlapping AAA in AAAAA at positions 0,1,2
        assert sc._exact_match_positions("AAAAA", "AAA") == [0, 1, 2]

    def test_pattern_longer_than_text(self):
        assert sc._exact_match_positions("AT", "ATGC") == []

    def test_empty_inputs(self):
        assert sc._exact_match_positions("", "ATGC") == []
        assert sc._exact_match_positions("ATGC", "") == []


# ═══════════════════════════════════════════════════════════════════════════════
# Demo plasmid — default no-arg launch content
# ═══════════════════════════════════════════════════════════════════════════════

class TestDemoPlasmid:
    """`_make_demo_record` builds a deterministic 1 kb circular record."""

    def test_length_1000(self):
        rec = sc._make_demo_record()
        assert len(rec.seq) == 1000
        # Literal constant must match (avoids accidental drift from
        # editing the literal but forgetting the comment)
        assert len(sc._DEMO_PLASMID_SEQ) == 1000

    def test_topology_circular(self):
        rec = sc._make_demo_record()
        assert rec.annotations["topology"] == "circular"
        assert rec.annotations["molecule_type"] == "DNA"

    def test_feature_count(self):
        rec = sc._make_demo_record()
        # ori + MCS + CDS + terminator
        assert len(rec.features) == 4

    def test_deterministic(self):
        rec1 = sc._make_demo_record()
        rec2 = sc._make_demo_record()
        assert str(rec1.seq) == str(rec2.seq)
        assert rec1.id == rec2.id

    def test_embedded_restriction_sites(self):
        # Gel simulator demo expects EcoRI / HindIII / BamHI / XbaI sites
        # to all be present. If someone edits _DEMO_PLASMID_SEQ they
        # need to know the demo loses its "click Run gel for cuts"
        # affordance.
        seq = str(sc._make_demo_record().seq)
        assert "GAATTC" in seq    # EcoRI
        assert "AAGCTT" in seq    # HindIII
        assert "GGATCC" in seq    # BamHI
        assert "TCTAGA" in seq    # XbaI

    async def test_demo_loads_without_persisting_to_library(
            self, tiny_record, isolated_library):
        """No-arg launch path: `_preload_demo_record` is applied to the
        canvas via `_apply_record` (not `_import_and_persist`), so the
        user's plasmid library is untouched."""
        from textual.events import MouseMove  # noqa: F401  (imported via app)
        app = sc.PlasmidApp()
        # Simulate no-arg launch: no _preload_record, demo set instead
        app._preload_record = None
        demo = sc._make_demo_record()
        app._preload_demo_record = demo
        # Disable seed-from-NCBI fallback (tests never hit network)
        app._skip_seed = True
        async with app.run_test(size=(160, 48)) as pilot:
            await pilot.pause()
            await pilot.pause(0.1)
            # The canvas should be showing the demo plasmid
            assert app._current_record is not None
            assert app._current_record.id == sc._DEMO_PLASMID_NAME
            # Library MUST NOT contain a demo entry
            lib = sc._load_library()
            names = {e.get("name", "") for e in lib}
            assert sc._DEMO_PLASMID_NAME not in names
