"""test_entry_reads — the ``entry`` parameter on the read endpoints
(2026-10-31, workbench spec W3).

A read endpoint answers about the plasmid on the CANVAS. That made a library
unbrowsable from a ``--read-only`` attachment: ``load-entry`` is a write (it
reseats a canvas a running GUI also owns), so a guest could list every entry
and read none of them.

``{"entry": "<name or id>"}`` says "answer about this library entry, do not
touch the canvas". The contract these tests pin:

* **Parity.** The entry answer is the answer you get after loading it — the
  same rows, the same coordinates, the same cutters. Measured by loading the
  entry through `load-entry` and comparing, feature by feature, including the
  cases that break naive coordinate code (a − strand CDS, one that crosses
  the origin, a spliced one).
* **Independence.** The answer does not depend on what is on the canvas, and
  does not change it. Including with nothing loaded at all.
* **It is a read.** It works under ``--read-only`` with the save
  chokepoint disarmed, because nothing in the path writes.
* **Refusal, not a wrong answer.** An unknown entry is a 404 naming it; an
  endpoint that does not implement ``entry`` refuses the key rather than
  answering about the canvas; ``collection`` on its own selects nothing.
"""
from __future__ import annotations

import asyncio
import random

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import splicecraft as sc

_TERM = (170, 50)
_N = 900

# A display name with a space and mixed case, because that is the shape
# [INV-98] exists for: the GenBank LOCUS cannot carry it, so a read keyed on
# it proves the entry's stored name is what resolves.
ENTRY = "shelf plasmid one"
ENTRY_ID = "ent-1"
COLL = "Shelf"
OTHER = "Other shelf"

# The reads that take `entry`, as the response-shaping differs per endpoint.
_SIMPLE_READS = ("get-sequence", "find-orfs", "list-restriction-sites",
                 "features", "list-features", "get-feature", "find-feature",
                 "plasmid-overview", "simulate-gel", "get-history")


def _cds(rng: random.Random, n_codons: int) -> str:
    sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
             if a + b + c not in ("TAA", "TAG", "TGA")]
    return "ATG" + "".join(rng.choice(sense) for _ in range(n_codons)) + "TAA"


def _fixture_record(seed: int = 20261031, name: str = "pSHELF1") -> SeqRecord:
    """A plasmid with the coordinate-nasty features: + strand, − strand, one
    crossing the origin, one spliced, and a unique cutter inside a gene."""
    rng = random.Random(seed)
    seq = [rng.choice("ACGT") for _ in range(_N)]

    def put(at: int, bases: str) -> None:
        for i, b in enumerate(bases):
            seq[(at + i) % _N] = b

    put(60, _cds(rng, 60))                                   # [60, 246) +
    put(300, str(Seq(_cds(rng, 40)).reverse_complement()))   # [300, 426) −
    put(_N - 48, _cds(rng, 30))                              # crosses origin
    put(500, _cds(rng, 20)[:33])
    put(560, "".join(rng.choice("ACG") for _ in range(27)) + "TAA")
    put(150, "GAATTC")                                       # EcoRI, in geneA
    rec = SeqRecord(Seq("".join(seq)), id=name, name=name,
                    description="entry-read fixture")
    rec.annotations["topology"] = "circular"
    rec.annotations["molecule_type"] = "DNA"
    rec.features = [
        # A `source` feature first, deliberately: it occupies no row in the
        # map's table, so an index walk that forgot to skip it would shift
        # every qualifier lookup by one — the bug `_feat_is_map_indexed`
        # exists to prevent, and the parity assertions below would catch it
        # on the entry path.
        SeqFeature(FeatureLocation(0, _N, strand=1), type="source",
                   qualifiers={"organism": ["synthetic construct"]}),
        SeqFeature(FeatureLocation(60, 246, strand=1), type="CDS",
                   qualifiers={"label": ["geneA"], "note": ["first"]}),
        SeqFeature(FeatureLocation(300, 426, strand=-1), type="CDS",
                   qualifiers={"label": ["geneB"]}),
        SeqFeature(CompoundLocation([FeatureLocation(_N - 48, _N, strand=1),
                                     FeatureLocation(0, 48, strand=1)]),
                   type="CDS", qualifiers={"label": ["wrapGene"]}),
        SeqFeature(CompoundLocation([FeatureLocation(500, 533, strand=1),
                                     FeatureLocation(560, 590, strand=1)]),
                   type="CDS", qualifiers={"label": ["splicedGene"]}),
        SeqFeature(FeatureLocation(700, 760, strand=1), type="misc_feature",
                   qualifiers={"label": ["marker"]}),
    ]
    return rec


def _entry_for(rec, *, name: str = ENTRY, eid: str = ENTRY_ID) -> dict:
    return {"id": eid, "name": name, "size": len(rec.seq),
            "gb_text": sc._record_to_gb_text(rec)}


def _seed(*, extra_collection: bool = False, rec=None):
    """Put the fixture in the library as the only entry of the active
    collection. Returns the record."""
    rec = rec if rec is not None else _fixture_record()
    colls = [{"name": COLL, "plasmids": [_entry_for(rec)]}]
    if extra_collection:
        colls.append({"name": OTHER, "plasmids": []})
    sc._save_collections(colls)
    sc._activate_collection(COLL)
    return rec


class _Bare:
    """The surface the entry path can legitimately touch: a canvas that holds
    nothing. Deliberately NOT a mock of anything under test — an entry read
    must not reach the app at all, and an attribute these handlers invent a
    use for would raise here."""

    _current_record = None
    _unsaved = False

    def call_from_thread(self, fn, *a, **k):   # pragma: no cover - unused
        raise AssertionError("an entry read must not marshal to the UI thread")


def _basis(out: dict) -> tuple:
    """Pop the basis stamp and return it, so the rest can be compared to the
    canvas answer field for field."""
    return out.pop("entry", None), out.pop("collection", None)


# ── parity with the loaded plasmid ───────────────────────────────────────────

async def test_entry_reads_match_the_plasmid_once_it_is_loaded():
    """The whole promise: reading an entry in place answers what loading it
    would. Compared field for field, on the record a `load-entry` produced,
    so a divergence in either path fails here."""
    _seed()
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        loaded = await asyncio.to_thread(sc._h_load_entry, app,
                                         {"name": ENTRY})
        assert loaded.get("ok") is True, loaded
        for _ in range(6):
            await pilot.pause()

        for payload, fn in (
            ({}, sc._h_features),
            ({"idx": 0}, sc._h_get_feature),
            ({"name": "geneA", "max_seq": 100000}, sc._h_find_feature),
            ({"name": "wrapGene"}, sc._h_find_feature),
            ({"name": "splicedGene"}, sc._h_find_feature),
            ({}, sc._h_plasmid_overview),
            ({"start": 60, "end": 70}, sc._h_get_sequence),
            ({"start": _N - 20, "end": 10}, sc._h_get_sequence),
            ({"min_aa": 25}, sc._h_find_orfs),
            ({"enzymes": ["EcoRI"]}, sc._h_list_restriction_sites),
            ({"unique_only": True}, sc._h_list_restriction_sites),
        ):
            canvas = fn(app, dict(payload))
            entry = fn(app, dict(payload, entry=ENTRY))
            assert isinstance(entry, dict), (fn.__name__, entry)
            assert _basis(entry) == (ENTRY, COLL), fn.__name__
            assert entry == canvas, (fn.__name__, payload)


async def test_every_feature_row_matches_including_wrap_and_splice():
    """Row by row, because the rows are the part built by a different code
    path: the canvas reads the map widget's table, the entry re-runs the same
    parse off-thread. `idx`, wrap encoding, exons and colour all have to
    agree or the two views of one plasmid disagree."""
    _seed()
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        await asyncio.to_thread(sc._h_load_entry, app, {"name": ENTRY})
        for _ in range(6):
            await pilot.pause()
        canvas = sc._h_features(app, {})["features"]
        entry_res = sc._h_features(app, {"entry": ENTRY})
        # The basis stamp first: without it this comparison would pass just
        # as well if `entry` were IGNORED, since the entry is also what is
        # loaded. It is the one assertion that makes the rest mean anything.
        assert _basis(dict(entry_res)) == (ENTRY, COLL)
        entry = entry_res["features"]
        assert canvas == entry
        labels = [r["label"] for r in entry]
        assert labels == ["geneA", "geneB", "wrapGene", "splicedGene",
                          "marker"]
        wrap = next(r for r in entry if r["label"] == "wrapGene")
        assert wrap["wraps"] is True and wrap["end"] < wrap["start"]
        # Each row's detail matches too, via the idx the row carries.
        for row in entry:
            got = sc._h_get_feature(app, {"entry": ENTRY, "idx": row["idx"]})
            assert got["label"] == row["label"]
            assert (got["start"], got["end"]) == (row["start"], row["end"])
        spliced = sc._h_find_feature(app, {"entry": ENTRY,
                                           "name": "splicedGene"})
        assert spliced["matches"][0]["parts"] == [[500, 533], [560, 590]]


# ── independence from the canvas ─────────────────────────────────────────────

async def test_an_entry_read_ignores_and_does_not_disturb_the_canvas():
    """A different plasmid is loaded. The entry answer is about the entry,
    the canvas is left byte-identical, and its dirty flag does not move."""
    other = _fixture_record(seed=777, name="pOTHER")
    _seed()
    app = sc.PlasmidApp()
    async with app.run_test(size=_TERM) as pilot:
        # The mount-time auto-load of the first library entry is scheduled
        # with `call_after_refresh`; let it land before seating our own
        # record, or it clobbers the canvas afterwards and this test compares
        # the entry against itself.
        for _ in range(6):
            await pilot.pause()
        app._apply_record(other)
        for _ in range(6):
            await pilot.pause()
        before = app._current_record
        before_seq = str(before.seq)
        app._unsaved = True
        # The canvas must really be a DIFFERENT molecule, or "it ignored the
        # canvas" would be unfalsifiable.
        assert before_seq != str(_fixture_record().seq)

        out = sc._h_plasmid_overview(app, {"entry": ENTRY})
        assert out["entry"] == ENTRY and out["length"] == _N
        # `unsaved` describes what was reported, not the canvas beside it.
        assert out["unsaved"] is False
        assert sc._h_get_sequence(app, {"entry": ENTRY, "start": 150,
                                        "end": 156})["seq"] == "GAATTC"
        canvas_same = sc._h_get_sequence(app, {"start": 150, "end": 156})
        assert canvas_same["seq"] == before_seq[150:156]

        assert app._current_record is before
        assert str(app._current_record.seq) == before_seq
        assert app._unsaved is True


def test_entry_reads_need_no_canvas_at_all():
    """The read-only session's actual case: attached, nothing loaded, and the
    library still readable. Every entry-capable read answers; the canvas-
    targeted call still refuses, so the two are not being confused."""
    rec = _seed()
    app = _Bare()
    assert sc._h_get_sequence(app, {"entry": ENTRY, "start": 0,
                                    "end": 6})["seq"] == str(rec.seq)[:6]
    assert sc._h_features(app, {"entry": ENTRY})["features"][0]["label"] == "geneA"
    assert sc._h_plasmid_overview(app, {"entry": ENTRY})["length"] == _N
    assert sc._h_find_orfs(app, {"entry": ENTRY})["count"] > 0
    # Ranked, forgiving matching is unchanged — `geneA` rides along as a
    # similar spelling, exactly as it does for the canvas.
    assert sc._h_find_feature(
        app, {"entry": ENTRY, "name": "geneB"})["matches"][0]["label"] == "geneB"
    assert sc._h_get_feature(app, {"entry": ENTRY, "idx": 1})["label"] == "geneB"
    assert sc._h_list_restriction_sites(
        app, {"entry": ENTRY, "enzymes": ["EcoRI"]})["count"] == 1
    assert sc._h_get_history(app, {"entry": ENTRY}).get("ok") is True
    body, status = sc._h_get_sequence(app, {"start": 0, "end": 6})
    assert status == 422 and "no plasmid" in body["error"]


def test_the_id_resolves_as_well_as_the_name():
    _seed()
    app = _Bare()
    by_id = sc._h_get_sequence(app, {"entry": ENTRY_ID, "start": 0, "end": 9})
    by_name = sc._h_get_sequence(app, {"entry": ENTRY, "start": 0, "end": 9})
    # Both name the stored DISPLAY name back, whichever key was used.
    assert by_id == by_name and by_id["entry"] == ENTRY


# ── it is a read ─────────────────────────────────────────────────────────────

def test_an_entry_read_writes_nothing_under_read_only(monkeypatch):
    """`--read-only` never arms the save chokepoint, so a `_save_*` anywhere
    in this path would raise. Proven by making every save raise and reading
    anyway — the fixture authorises writes for the sandbox, which would have
    hidden a stray one."""
    _seed()

    def _boom(*a, **k):
        raise AssertionError("an entry read must not save")

    monkeypatch.setattr(sc, "_safe_save_json", _boom)
    monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY", True)
    monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY_HOLDER_PID", 4242)
    app = _Bare()
    assert sc._h_plasmid_overview(app, {"entry": ENTRY})["entry"] == ENTRY
    # Through the dispatcher, which is what a read-only guest actually hits.
    out, status = (sc._agent_invoke(app, "features", {"entry": ENTRY}), 200)
    body = out[0] if isinstance(out, tuple) else out
    assert isinstance(body, dict) and body.get("entry") == ENTRY
    refused, code = sc._agent_invoke(app, "load-entry", {"name": ENTRY})
    assert code == 409 and refused["holder_pid"] == 4242


def test_status_reports_read_only_and_the_lock_holder(monkeypatch):
    """W1: a served UI could only discover the mode by attempting a write and
    catching the 409. `status` now says so up front, naming the same PID the
    refusal names."""
    app = _Bare()
    out = sc._h_status(app, {})
    assert out["read_only"] is False and out["lock_holder_pid"] is None
    monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY", True)
    monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY_HOLDER_PID", 31337)
    out = sc._h_status(app, {})
    assert out["read_only"] is True and out["lock_holder_pid"] == 31337
    refused, code = sc._agent_invoke(app, "save", {})
    assert code == 409
    assert refused["holder_pid"] == out["lock_holder_pid"]


def test_a_recorded_pid_of_zero_reads_as_null(monkeypatch):
    """Nothing recorded is `null`, not `0` — a UI printing "quit PID 0" would
    be pointing at the kernel."""
    monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY", True)
    monkeypatch.setattr(sc._state, "_AGENT_READ_ONLY_HOLDER_PID", 0)
    assert sc._h_status(_Bare(), {})["lock_holder_pid"] is None


# ── refusals ─────────────────────────────────────────────────────────────────

def test_an_unknown_entry_is_a_404_that_names_it():
    """"No such entry" and "no sites on that entry" must not read alike —
    the distinction `digest` already makes for an unknown enzyme."""
    _seed()
    app = _Bare()
    body, status = sc._h_list_restriction_sites(app, {"entry": "nope"})
    assert status == 404 and "nope" in body["error"]
    body, status = sc._h_features(app, {"entry": "nope"})
    assert status == 404 and body.get("features") is None


@pytest.mark.parametrize("inert", [None, "", [], {}, False])
def test_an_inert_entry_names_no_entry_and_is_not_an_error(inert):
    """`null` / `""` / `[]` / `{}` / `false` name nothing. A typed client or
    an LLM tool layer fills every field it knows about, and refusing on a
    key's presence alone is how ~120 write endpoints broke for those callers
    in 2026-09. So an inert `entry` leaves the canvas as the target — here
    there is none, so the endpoint's own "no plasmid loaded" answers."""
    _seed()
    app = _Bare()
    body, status = sc._h_get_sequence(app, {"entry": inert, "start": 0,
                                            "end": 4})
    assert status == 422 and "no plasmid" in body["error"]
    assert sc._h_features(app, {"entry": inert}) == {"features": []}
    # And the dispatcher's refusal for endpoints that don't take `entry` must
    # not fire on one either.
    out = sc._agent_invoke(app, "digest", {"entry": inert,
                                            "sequence": "GAATTC",
                                            "enzymes": ["EcoRI"]})
    body = out[0] if isinstance(out, tuple) else out
    assert "does not take 'entry'" not in str(body.get("error", ""))


@pytest.mark.parametrize("bad", ["   ", "\t", 7, 3.5, True, ["pUC19"],
                                 {"name": "pUC19"}])
def test_a_selector_that_looks_real_and_selects_nothing_is_refused(bad):
    """The other half of the rule: a whitespace string or a number is a
    mistake, not an unset field, and the canvas is not a sane fallback for
    it — the caller plainly meant an entry."""
    _seed()
    body, status = sc._h_get_sequence(_Bare(), {"entry": bad, "start": 0,
                                                "end": 4})
    assert status == 400, (bad, body)
    assert "entry" in body["error"]


@pytest.mark.parametrize("inert", [None, "", []])
def test_an_inert_collection_is_ignored_too(inert):
    _seed()
    app = _Bare()
    assert sc._h_features(app, {"collection": inert}) == {"features": []}
    out = sc._h_get_sequence(app, {"entry": ENTRY, "collection": inert,
                                   "start": 0, "end": 4})
    assert out["collection"] == COLL


def test_collection_alone_names_nothing():
    _seed()
    body, status = sc._h_features(_Bare(), {"collection": COLL})
    assert status == 400 and "only qualifies 'entry'" in body["error"]


def test_collection_pins_the_search_and_a_wrong_one_is_a_404():
    rec = _fixture_record()
    sc._save_collections([
        {"name": COLL, "plasmids": [_entry_for(rec)]},
        {"name": OTHER, "plasmids": []},
    ])
    sc._activate_collection(OTHER)
    app = _Bare()
    # Not in the ACTIVE collection, found anyway, and named for what it is.
    out = sc._h_get_sequence(app, {"entry": ENTRY, "start": 0, "end": 4})
    assert out["collection"] == COLL
    out = sc._h_get_sequence(app, {"entry": ENTRY, "collection": COLL,
                                   "start": 0, "end": 4})
    assert out["collection"] == COLL
    body, status = sc._h_get_sequence(app, {"entry": ENTRY,
                                            "collection": OTHER,
                                            "start": 0, "end": 4})
    assert status == 404
    body, status = sc._h_get_sequence(app, {"entry": ENTRY,
                                            "collection": "no such",
                                            "start": 0, "end": 4})
    assert status == 404 and "no collection named" in body["error"]


def test_a_key_in_two_collections_is_a_409_listing_them():
    """Identical to `load-entry`'s answer, because it IS `load-entry`'s
    resolver — a key must not mean one molecule to a load and another to a
    read."""
    rec = _fixture_record()
    sc._save_collections([
        {"name": "Left", "plasmids": [_entry_for(rec, eid="l1")]},
        {"name": "Right", "plasmids": [_entry_for(rec, eid="r1")]},
    ])
    sc._activate_collection("Left")
    app = _Bare()
    # The active collection wins while it holds a copy (historical, silent).
    assert sc._h_get_sequence(app, {"entry": ENTRY, "start": 0,
                                    "end": 4})["collection"] == "Left"
    sc._save_collections([
        {"name": "Elsewhere", "plasmids": []},
        {"name": "Left", "plasmids": [_entry_for(rec, eid="l1")]},
        {"name": "Right", "plasmids": [_entry_for(rec, eid="r1")]},
    ])
    sc._activate_collection("Elsewhere")
    body, status = sc._h_get_sequence(app, {"entry": ENTRY, "start": 0,
                                            "end": 4})
    assert status == 409 and body["collections"] == ["Left", "Right"]
    load_body, load_status = sc._h_load_entry(app, {"name": ENTRY})
    assert (load_status, load_body["collections"]) == (409, ["Left", "Right"])


def test_an_entry_with_no_stored_sequence_is_a_422():
    sc._save_collections([{"name": COLL, "plasmids": [
        {"id": "empty-1", "name": "no bases here", "size": 0}]}])
    sc._activate_collection(COLL)
    body, status = sc._h_features(_Bare(), {"entry": "no bases here"})
    assert status == 422 and "no stored sequence" in body["error"]


def test_endpoints_that_do_not_take_entry_refuse_it():
    """The silent-wrong this closes: a read endpoint echoes no `ignored`
    list, so a dropped `entry` would have answered about the canvas and read
    as an answer about the entry."""
    _seed()
    app = _Bare()
    for name in ("digest", "analyse-cds", "list-library", "load-entry",
                 "delete-feature", "add-current-to-library"):
        body, status = sc._agent_invoke(app, name, {"entry": ENTRY})
        assert status == 400, (name, status, body)
        assert body["unsupported"] == ["entry"]
        assert "get-sequence" in body["error"]
        assert set(body["entry_accepted_by"]) == set(
            sc._AGENT_ENTRY_PARAM_ENDPOINTS)


def test_the_entry_capable_list_matches_the_handlers_that_read_it():
    """Drift guard. The dispatcher refuses `entry` on anything outside this
    frozenset, so an endpoint that grows `entry` support without joining the
    set becomes unreachable — and one that leaves the set while still reading
    the key would read it from a payload the dispatcher let through."""
    for name, (fn, _write) in sc._state._AGENT_HANDLERS.items():
        keys = sc._agent_payload_keys(fn)
        if keys is None:                 # analysis could not establish them
            continue
        reads_entry = "entry" in keys
        in_set = name in sc._AGENT_ENTRY_PARAM_ENDPOINTS
        assert reads_entry == in_set, (
            f"{name}: reads 'entry'={reads_entry} but in the entry-capable "
            f"set={in_set}")
    # Every member resolves to a registered READ endpoint.
    for name in sc._AGENT_ENTRY_PARAM_ENDPOINTS:
        fn, write = sc._state._AGENT_HANDLERS[name]
        assert write is False, f"{name} writes; `entry` is a read parameter"


# ── the two endpoints with a selector of their own ──────────────────────────

def test_simulate_gel_takes_its_template_from_the_entry():
    rec = _seed()
    app = _Bare()
    from_entry = sc._h_simulate_gel(app, {
        "entry": ENTRY, "lanes": [{"source": "digest", "detail": "EcoRI"}]})
    explicit = sc._h_simulate_gel(app, {
        "template_seq": str(rec.seq), "template_circular": True,
        "lanes": [{"source": "digest", "detail": "EcoRI"}]})
    assert _basis(from_entry) == (ENTRY, COLL)
    assert from_entry == explicit
    # One EcoRI site on a circle: one full-length linear band.
    assert [b["bp"] for b in from_entry["lanes"][0]["bands"]] == [_N]


def test_simulate_gel_refuses_two_templates():
    _seed()
    app = _Bare()
    body, status = sc._h_simulate_gel(
        app, {"entry": ENTRY, "template_seq": "ACGTACGT",
              "lanes": [{"source": "plasmid"}]})
    assert status == 422 and "EITHER" in body["error"]
    # An EMPTY `template_seq` names no template — inert, so not a conflict.
    ok = sc._h_simulate_gel(app, {"entry": ENTRY, "template_seq": "",
                                  "lanes": [{"source": "plasmid"}]})
    assert ok["entry"] == ENTRY


def test_the_gel_echoes_the_topology_it_ran():
    """`template_circular` is NOT refused alongside `entry` — `false` is the
    unset value of a boolean for most typed clients, so refusing it would
    reject requests that meant nothing by it. The entry's own topology wins
    and the answer says which it used, because a supercoiled band on what the
    caller believes is a linear fragment is the whole lane misread."""
    _seed()
    app = _Bare()
    out = sc._h_simulate_gel(app, {"entry": ENTRY, "template_circular": False,
                                   "lanes": [{"source": "plasmid"}]})
    assert out["template_circular"] is True        # the entry is circular
    forms = {b["form"] for b in out["lanes"][0]["bands"]}
    assert "supercoiled" in forms or "circular" in " ".join(forms)
    # Without `entry` the caller's flag is still the caller's flag.
    plain = sc._h_simulate_gel(app, {"template_seq": "ACGT" * 40,
                                     "template_circular": False,
                                     "lanes": [{"source": "plasmid"}]})
    assert plain["template_circular"] is False


def test_a_linear_entry_runs_as_linear_on_the_gel():
    """The entry's own topology is used — a supercoiled band for a linear
    fragment (or the reverse) would be the gel lying about the molecule."""
    lin = _fixture_record()
    lin.annotations["topology"] = "linear"
    _seed(rec=lin)
    out = sc._h_simulate_gel(_Bare(), {"entry": ENTRY,
                                       "lanes": [{"source": "plasmid"}]})
    forms = {b["form"] for b in out["lanes"][0]["bands"]}
    assert forms == {"linear"}


def test_get_history_entry_is_the_same_selector_not_a_second_one():
    rec = _fixture_record()
    sc._save_collections([
        {"name": "Elsewhere", "plasmids": []},
        {"name": COLL, "plasmids": [_entry_for(rec)]},
    ])
    sc._activate_collection("Elsewhere")
    app = _Bare()
    # `name` searches the active collection's mirror (unchanged); `entry`
    # resolves across collections like every other `entry`.
    body, status = sc._h_get_history(app, {"name": ENTRY})
    assert status == 404
    out = sc._h_get_history(app, {"entry": ENTRY})
    assert out.get("ok") is True and out["collection"] == COLL
    body, status = sc._h_get_history(app, {"entry": ENTRY, "name": ENTRY})
    assert status == 422 and "EITHER" in body["error"]


# ── the off-canvas feature parse ────────────────────────────────────────────

def test_the_off_canvas_parse_is_the_map_s_own_parse():
    """`_agent_map_rows_for_record` runs `PlasmidMap._parse` against a
    scratch `self` with `__slots__`: if a future `_parse` reaches for another
    attribute of the widget, this raises instead of returning a feature table
    built from half a map."""
    rec = _fixture_record()
    rows = sc._agent_map_rows_for_record(rec)
    assert [r["label"] for r in rows] == ["geneA", "geneB", "wrapGene",
                                          "splicedGene", "marker"]
    scratch = sc._FeatRowScratch()
    with pytest.raises(AttributeError):
        scratch.record = rec                     # type: ignore[attr-defined]
    assert (scratch._n_skipped, scratch._n_flattened) == (0, 0)


# ── hardening, 2026-11-01 ───────────────────────────────────────────────────

def test_topology_comes_from_the_record_that_was_read():
    """The deleted `_agent_record_is_circular(app)` read the CANVAS's
    topology. With a linear entry under a circular canvas, that would have
    answered an origin-wrap request about a molecule with no origin to wrap —
    the exact shape `entry` exists to stop."""
    linear = _fixture_record()
    linear.annotations["topology"] = "linear"
    _seed(rec=linear)
    app = sc.PlasmidApp()

    async def _run():
        async with app.run_test(size=_TERM) as pilot:
            # Let the mount-time auto-load of the first library entry land
            # FIRST — it is scheduled with `call_after_refresh`, so applying
            # before it fires has the autoload clobber the canvas afterwards.
            for _ in range(6):
                await pilot.pause()
            circular = _fixture_record(seed=99, name="pCIRC")
            circular.annotations["topology"] = "circular"
            app._apply_record(circular)
            for _ in range(6):
                await pilot.pause()
            assert app._current_record.annotations["topology"] == "circular"
            # The canvas is circular: a wrap answers.
            assert sc._h_get_sequence(app, {"start": _N - 5,
                                            "end": 5})["length"] == 10
            # The entry is linear: the same wrap is refused, naming why.
            body, status = sc._h_get_sequence(app, {"entry": ENTRY,
                                                    "start": _N - 5,
                                                    "end": 5})
            assert status == 400 and "linear" in body["error"]

    import asyncio
    asyncio.run(_run())


def test_a_stored_sequence_that_is_not_text_is_a_422(monkeypatch):
    """A hand-edited `collections.json` can put anything in `gb_text`.
    Handing it to the parser raised deep inside and came back as a scrubbed
    500, which reads as "SpliceCraft broke" rather than "that entry is not
    storing text". The save path normalises, so only the data source is
    stubbed here — never the resolution under test."""
    _seed()
    monkeypatch.setattr(
        "splicecraft_agent._agent_scan_library_for_key",
        lambda key, coll=None: [(COLL, {"id": "x", "name": ENTRY,
                                        "gb_text": {"not": "text"}})])
    body, status = sc._h_features(_Bare(), {"entry": ENTRY})
    assert status == 422 and "not text" in body["error"]


def test_the_off_canvas_parse_logs_one_line_not_one_per_feature():
    """A LOAD narrates every wrap and every flattened compound feature: one
    plasmid, once, with someone watching. An entry read can run on every
    click, and 800 spliced genes wrote 800 INFO lines per call into a
    ROTATING log — diagnostic history spent on the same sentence 800 times.
    The aggregate line carries the same counts."""
    import logging
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

    rec = _fixture_record()
    for i in range(40):
        s = 100 + i * 10
        rec.features.append(SeqFeature(
            CompoundLocation([FeatureLocation(s, s + 3, strand=1),
                              FeatureLocation(s + 5, s + 8, strand=1)]),
            type="CDS", qualifiers={"label": [f"sp{i}"]}))

    records: "list[logging.LogRecord]" = []

    class _Grab(logging.Handler):
        def emit(self, record):
            records.append(record)

    logger = logging.getLogger(sc._log.name)
    handler = _Grab()
    prev = logger.level
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    try:
        rows = sc._agent_map_rows_for_record(rec)
    finally:
        logger.removeHandler(handler)
        logger.setLevel(prev)
    assert len(rows) == len(rec.features) - 1          # `source` holds no row
    said = [r.getMessage() for r in records]
    assert not [m for m in said if "Compound feature" in m], said[:3]
    aggregate = [m for m in said if "off-canvas parse" in m]
    assert len(aggregate) == 1, said
    assert "41 flattened" in aggregate[0], aggregate


def test_a_load_still_narrates_every_feature():
    """The quiet flag is the off-canvas caller's, not the map's: a load says
    exactly what it said before."""
    import logging
    rec = _fixture_record()
    records: "list[logging.LogRecord]" = []

    class _Grab(logging.Handler):
        def emit(self, record):
            records.append(record)

    app = sc.PlasmidApp()

    async def _run():
        logger = logging.getLogger(sc._log.name)
        handler = _Grab()
        prev = logger.level
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
        try:
            async with app.run_test(size=_TERM) as pilot:
                app._apply_record(rec)
                for _ in range(6):
                    await pilot.pause()
        finally:
            logger.removeHandler(handler)
            logger.setLevel(prev)

    import asyncio
    asyncio.run(_run())
    said = [r.getMessage() for r in records]
    assert [m for m in said if "Compound feature" in m], said[:5]
    assert [m for m in said if "Detected wrap feature" in m], said[:5]
    assert not [m for m in said if "off-canvas parse" in m]


# ── check-acceptor-vector joins the `entry` family (2026-11-08) ────────────


@pytest.mark.usefixtures("_protect_user_data")
def test_check_acceptor_vector_answers_about_an_entry():
    """Screening a whole shelf of entry vectors is the motivating use, and it
    must not reseat the canvas to do it."""
    _seed()
    app = type("A", (), {"_current_record": None})()
    res = sc._h_check_acceptor_vector(app, {"entry": ENTRY})
    assert not isinstance(res, tuple), res
    assert res["basis"]["entry"] == ENTRY
    assert "enzymes" in res and "role" in res
    assert app._current_record is None


@pytest.mark.usefixtures("_protect_user_data")
def test_check_acceptor_vector_404s_an_unknown_entry():
    _seed()
    app = type("A", (), {"_current_record": None})()
    body, code = sc._h_check_acceptor_vector(app, {"entry": "no such plasmid"})
    assert code == 404 and "no such plasmid" in body["error"]


def test_check_acceptor_vector_is_in_the_entry_allow_set():
    import splicecraft_agent as _ag
    assert "check-acceptor-vector" in _ag._AGENT_ENTRY_PARAM_ENDPOINTS
