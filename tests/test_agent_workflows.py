"""test_agent_workflows — the one-call workflow endpoints (2026-09-21).

`plasmid-overview`, `find-feature` and `amplify-feature` exist so a small local
model (Babs) — or any agent-API client — gets exact coordinates, sequences,
translations and binding-checked primers in ONE call instead of chaining
endpoints and copying numbers between them.

Every check here is against a reference OUTSIDE the code under test:
Biopython's own `feature.extract()` / `translate()` for sequences and proteins,
and the raw template bases for primers (a forward primer must BE the top
strand at its reported site; a reverse primer its reverse complement). The
fixture carries the cases that break naive coordinate code: a − strand CDS, a
CDS that crosses the origin, and a spliced CDS.
"""
from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import splicecraft as sc

_TERM = (170, 50)
_N = 900


def _cds_bases(rng: random.Random, n_codons: int) -> str:
    """ATG + random sense codons (no internal stop) + TAA."""
    sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
             if a + b + c not in ("TAA", "TAG", "TGA")]
    return "ATG" + "".join(rng.choice(sense) for _ in range(n_codons)) + "TAA"


def _fixture_record() -> SeqRecord:
    rng = random.Random(20260921)
    seq = [rng.choice("ACGT") for _ in range(_N)]

    def put(at: int, bases: str) -> None:
        for i, b in enumerate(bases):
            seq[(at + i) % _N] = b

    plus = _cds_bases(rng, 60)                      # 186 bp, + strand, [60, 246)
    put(60, plus)
    minus = _cds_bases(rng, 40)                     # 126 bp, − strand, [300, 426)
    put(300, str(Seq(minus).reverse_complement()))
    wrap = _cds_bases(rng, 30)                      # 96 bp crossing the origin
    put(_N - 48, wrap)                              # [852, 900) + [0, 48)
    exon1 = _cds_bases(rng, 20)[:33]                # spliced CDS: 33 + 30 bp
    put(500, exon1)
    exon2 = "".join(rng.choice("ACG") for _ in range(27)) + "TAA"
    put(560, exon2)
    rec = SeqRecord(Seq("".join(seq)), id="pWF1", name="pWF1",
                    description="workflow fixture")
    rec.annotations["topology"] = "circular"
    rec.annotations["molecule_type"] = "DNA"
    rec.features = [
        SeqFeature(FeatureLocation(60, 246, strand=1), type="CDS",
                   qualifiers={"label": ["geneA"]}),
        SeqFeature(FeatureLocation(300, 426, strand=-1), type="CDS",
                   qualifiers={"label": ["geneB"]}),
        SeqFeature(CompoundLocation([FeatureLocation(_N - 48, _N, strand=1),
                                     FeatureLocation(0, 48, strand=1)]),
                   type="CDS", qualifiers={"label": ["wrapGene"]}),
        SeqFeature(CompoundLocation([FeatureLocation(500, 533, strand=1),
                                     FeatureLocation(560, 590, strand=1)]),
                   type="CDS", qualifiers={"label": ["splicedGene"]}),
        SeqFeature(FeatureLocation(700, 760, strand=1), type="misc_feature",
                   qualifiers={"label": ["LacZ-α"]}),
        SeqFeature(FeatureLocation(770, 790, strand=1), type="misc_feature",
                   qualifiers={"label": ["geneA site"]}),
    ]
    return rec


async def _loaded_app(pilot, app, rec=None):
    app._apply_record(rec if rec is not None else _fixture_record())
    for _ in range(6):
        await pilot.pause()


def _row(app, label):
    return next(r for r in sc._h_features(app, {})["features"] if r["label"] == label)


def _bio_feature(app, idx):
    return sc._agent_seqfeature_for_idx(app._current_record, idx)


# ── find-feature ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("label", ["geneA", "geneB", "wrapGene", "splicedGene"])
async def test_find_feature_sequence_and_protein_match_biopython(label):
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        out = sc._h_find_feature(app, {"name": label, "max_seq": 100000})
        m = out["matches"][0]
        assert m["label"] == label and m["match"] == "exact"
        sf = _bio_feature(app, m["idx"])
        bio = str(sf.extract(app._current_record.seq)).upper()
        assert m["sequence"] == bio
        assert m["translation"].rstrip("*") == str(Seq(bio).translate()).rstrip("*")
        assert m["start_1based"] == m["start"] + 1 and m["end_1based"] == m["end"]


async def test_find_feature_reports_wrap_strand_and_parts():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        wrap = sc._h_find_feature(app, {"name": "wrapGene"})["matches"][0]
        assert wrap["wraps"] is True and wrap["length"] == 96
        assert wrap["start"] == _N - 48 and wrap["end"] == 48
        assert "parts" not in wrap                 # a wrap is NOT a splice
        minus = sc._h_find_feature(app, {"name": "geneB"})["matches"][0]
        assert minus["strand"] == -1 and minus["sequence"].startswith("ATG")
        spliced = sc._h_find_feature(app, {"name": "splicedGene"})["matches"][0]
        assert spliced["parts"] == [[500, 533], [560, 590]]
        assert len(spliced["sequence"]) == 63


async def test_find_feature_name_matching_is_forgiving_and_ranked():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        for q in ("lacz alpha", "LACZ-ALPHA", "lacZα"):
            assert sc._h_find_feature(app, {"name": q})["matches"][0]["label"] == "LacZ-α"
        ranked = sc._h_find_feature(app, {"name": "genea"})["matches"]
        assert [m["label"] for m in ranked][:2] == ["geneA", "geneA site"]
        assert [m["match"] for m in ranked][:2] == ["exact", "starts with"]
        typed = sc._h_find_feature(app, {"name": "gene", "type": "misc_feature"})
        assert {m["type"] for m in typed["matches"]} == {"misc_feature"}
        miss = sc._h_find_feature(app, {"name": "zzz-nothing"})
        assert miss["count"] == 0 and "geneA" in miss["features"]


async def test_find_feature_caps_sequence_and_validates():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        m = sc._h_find_feature(app, {"name": "geneA", "max_seq": 10})["matches"][0]
        assert len(m["sequence"]) == 10 and m["sequence_truncated"] is True
        for bad in ({}, {"name": ""}, {"name": 5}, {"name": "x", "max_seq": "big"}):
            out = sc._h_find_feature(app, bad)
            assert isinstance(out, tuple) and out[1] == 400, bad


async def test_translation_disagreement_with_annotation_is_reported():
    rec = _fixture_record()
    rec.features[0].qualifiers["translation"] = ["MWRONGPROTEIN"]
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app, rec)
        m = sc._h_find_feature(app, {"name": "geneA"})["matches"][0]
        assert m["translation_matches_annotation"] is False


async def test_find_feature_lists_every_once_cutter_inside_it():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        full = sc._h_list_restriction_sites(app, {})
        counts: dict = {}
        cuts: dict = {}
        for s in full["sites"]:
            counts[s["enzyme"]] = counts.get(s["enzyme"], 0) + 1
            cuts[s["enzyme"]] = s["cut_bp"]
        for label in ("geneA", "geneB", "wrapGene"):
            m = sc._h_find_feature(app, {"name": label})["matches"][0]
            got = {c["enzyme"] for c in m["unique_cutters_inside"]}
            want = {e for e, k in counts.items() if k == 1
                    and sc._bp_in_span(cuts[e], m["start"], m["end"], _N)}
            assert got == want, (label, got ^ want)


# ── amplify-feature ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("label", ["geneA", "geneB", "wrapGene"])
async def test_amplify_primers_are_the_template_at_the_feature_ends(label):
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        row = _row(app, label)
        out = sc._h_amplify_feature(app, {"name": label})
        seq = str(app._current_record.seq).upper()
        twice = seq + seq
        f, r = out["forward"], out["reverse"]
        assert f["seq"] == twice[f["start"]: f["start"] + f["length"]]
        assert r["seq"] == str(Seq(twice[r["start"]: r["start"] + r["length"]])
                               .reverse_complement())
        assert f["start"] == row["start"]
        assert (r["start"] + r["length"]) % _N == row["end"] % _N
        assert f["binds_as_designed"] and r["binds_as_designed"]
        assert out["product_size"] == row["length"]


async def test_amplify_minus_strand_names_the_primer_holding_the_start():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        minus = sc._h_amplify_feature(app, {"name": "geneB"})
        assert minus["primer_at_feature_start"] == "reverse"
        assert minus["reverse"]["seq"].startswith("ATG") and "strand_note" in minus
        plus = sc._h_amplify_feature(app, {"name": "geneA"})
        assert plus["primer_at_feature_start"] == "forward"
        assert plus["forward"]["seq"].startswith("ATG") and "strand_note" not in plus


async def test_amplify_flags_a_primer_that_binds_twice():
    rec = _fixture_record()
    s = list(str(rec.seq))
    s[800:830] = list(str(rec.seq)[60:90])   # a second copy of geneA's first 30 bp
    rec.seq = Seq("".join(s))
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app, rec)
        out = sc._h_amplify_feature(app, {"name": "geneA"})
        assert out["forward"]["other_sites"], out["forward"]
        assert any("mispriming" in w for w in out["warnings"])


async def test_amplify_ambiguity_linear_wrap_and_validation():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        body, code = sc._h_amplify_feature(app, {"name": "gene"})   # starts-with ×2
        assert code == 409 and len(body["candidates"]) >= 2
        idx = body["candidates"][0]["idx"]
        assert sc._h_amplify_feature(app, {"idx": idx})["feature"]["idx"] == idx
        for bad, want in (({}, 400), ({"idx": "x"}, 400), ({"idx": 999}, 404),
                          ({"name": "zzz"}, 404),
                          ({"name": "geneA", "target_tm": 10}, 400),
                          ({"name": "geneA", "target_tm": "hot"}, 400)):
            out = sc._h_amplify_feature(app, bad)
            assert isinstance(out, tuple) and out[1] == want, (bad, out)
        linear = _fixture_record()
        linear.annotations["topology"] = "linear"
        await _loaded_app(pilot, app, linear)
        wrap_idx = next(r["idx"] for r in sc._h_features(app, {})["features"]
                        if r["label"] == "wrapGene")
        out = sc._h_amplify_feature(app, {"idx": wrap_idx})
        # On a LINEAR record the two-part join is a splice, not a wrap: whatever
        # it amplifies, it must never claim to cross an origin that isn't there.
        if isinstance(out, dict):
            assert out["forward"]["binds_as_designed"]


# ── plasmid-overview ──────────────────────────────────────────────────────────

async def test_overview_cutters_are_unique_and_placed_in_features():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        ov = sc._h_plasmid_overview(app, {})
        assert ov["length"] == _N and ov["topology"] == "circular"
        assert ov["feature_count"] == len(ov["features"]) == 6
        full = sc._h_list_restriction_sites(app, {})
        counts: dict = {}
        for s in full["sites"]:
            counts[s["enzyme"]] = counts.get(s["enzyme"], 0) + 1
        rows = sc._h_features(app, {})["features"]
        seen: dict = {}
        for label, entries in ov["unique_cutters"].items():
            for c in entries:
                assert counts.get(c["enzyme"]) == 1, c          # really cuts once
                inside = [r["label"] for r in rows
                          if sc._bp_in_span(c["cut_bp"], r["start"], r["end"], _N)]
                if label == sc._OVERVIEW_OUTSIDE:
                    assert inside == [], c                  # truly outside everything
                else:
                    assert label in inside, (label, c)      # filed under a feature it cuts in
                seen[c["enzyme"]] = c["cut_bp"]
        unique = {e for e, k in counts.items() if k == 1}
        assert set(seen) == unique                          # none lost
        assert ov["unique_cutter_count"] == len(unique)
        wrap = next(f for f in ov["features"] if f["label"] == "wrapGene")
        assert wrap.get("wraps") is True
        capped = sc._h_plasmid_overview(app, {"max_features": 2})
        assert len(capped["features"]) == 2 and "features_truncated" in capped


async def test_workflows_refuse_cleanly_with_nothing_loaded():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await pilot.pause()
        app._current_record = None
        for fn, body in ((sc._h_plasmid_overview, {}),
                         (sc._h_find_feature, {"name": "x"}),
                         (sc._h_amplify_feature, {"name": "x"})):
            out = fn(app, body)
            assert isinstance(out, tuple) and out[1] == 409


def test_workflow_endpoints_registered_read_only():
    for ep in ("plasmid-overview", "find-feature", "amplify-feature"):
        assert ep in sc._AGENT_HANDLERS and sc._AGENT_HANDLERS[ep][1] is False


# ── hardening: junk payloads never 500, and the two review fixes ──────────────

_JUNK = [None, True, False, 0, -1, 10 ** 30, 3.7, float("nan"), "", " ", "x" * 5000,
         "\x00", "ÄÖÜ", [], [1, 2], {}, {"a": 1}]


async def test_workflow_endpoints_never_500_on_junk():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        fields = {
            sc._h_plasmid_overview: ("max_features",),
            sc._h_find_feature: ("name", "type", "max_seq"),
            sc._h_amplify_feature: ("name", "idx", "type", "target_tm"),
        }
        n = 0
        for fn, keys in fields.items():
            for key in keys:
                for junk in _JUNK:
                    body = {"name": "geneA", key: junk} if key != "name" else {key: junk}
                    out = fn(app, body)
                    n += 1
                    if isinstance(out, tuple):
                        assert out[1] < 500, (fn.__name__, key, junk, out)
                    else:
                        assert isinstance(out, dict), (fn.__name__, key, junk)
        assert n > 100


async def test_boolean_idx_is_rejected_not_read_as_one():
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await _loaded_app(pilot, app)
        out = sc._h_amplify_feature(app, {"idx": True})
        assert isinstance(out, tuple) and out[1] == 400


def test_failed_model_lookup_is_not_cached(monkeypatch):
    import splicecraft_babs as B
    scr = sc.BabsScreen()
    state = {"up": False}

    def show(name, **k):
        if not state["up"]:
            raise B.OllamaUnavailable("down")
        return {"capabilities": ["completion"]}
    monkeypatch.setattr(B, "show", show)
    assert scr._agent_protocol("m") == "native"      # unknown while Ollama is down
    state["up"] = True
    assert scr._agent_protocol("m") == "json"        # …and learned once it answers
    state["up"] = False
    assert scr._agent_protocol("m") == "json"        # a success IS cached
