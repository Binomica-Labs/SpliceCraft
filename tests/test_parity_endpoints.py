"""
test_parity_endpoints — the UI→side-door parity endpoints added 2026-10-16.

A 2026-10-16 audit mapped all 338 `action_*` methods, 429 bindings and every
mutating button against the endpoint registry and found twelve capabilities
the GUI had and the agent API did not. These cover the nine that were built:

  set-origin / reverse-complement   whole-record reframes (Alt+Shift+O / +R)
  save-amplicon                     file a simulated PCR product
  import-primers                    the inverse `export-primers` advertised
  update-stored-primer              rename / order status / order-cart flag
  taxon-search                      organism name → NCBI taxid
  read-trace                        Sanger .ab1 bases + Phred
  protein collections (7)           an entire store with no HTTP reach
  codon_forbidden_enzymes           the setting the two surfaces disagreed on

The oracle for the reframes is the BASES, not the coordinates — the engines
move every coordinate at once, so only extracted sequence proves conservation.
"""
from __future__ import annotations


import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

import splicecraft as sc
import splicecraft_agent as sca


def _call(_ep, _app=None, **payload):
    """Invoke a registered handler directly, as the dispatcher would.

    Params are underscore-prefixed so a payload key called `name` (several
    endpoints take one) cannot collide with the helper's own signature."""
    entry = sc._state._AGENT_HANDLERS.get(_ep)
    assert entry is not None, f"endpoint {_ep!r} is not registered"
    return entry[0](_app, payload)


def _status(res):
    return res[1] if isinstance(res, tuple) else 200


def _body(res):
    return res[0] if isinstance(res, tuple) else res


# ═══════════════════════════════════════════════════════════════════════════
# set-origin / reverse-complement
# ═══════════════════════════════════════════════════════════════════════════

class _ReframeApp:
    """Only the UI plumbing is stubbed — `_apply_snapshot`, the alignment
    invalidation and the seq-panel lookup. The transforms under test
    (`_rotate_seq_record` / `_reverse_complement_record`) are module-level
    functions and run for real, so a conserved feature here is the engine's
    doing, not the harness's."""

    class _SP:
        _cursor_pos = 0

    def __init__(self, record):
        self._current_record = record
        self._undo: list = []
        self._cleared = 0

    def call_from_thread(self, fn, *a, **kw):
        return fn(*a, **kw)

    def _push_undo(self):
        self._undo.append(self._current_record)

    def _pop_undo_silently(self):
        if self._undo:
            self._current_record = self._undo.pop()

    def _invalidate_alignments_for_reframe(self, rec):
        return self._cleared

    def _apply_snapshot(self, seq, cursor, record):
        self._current_record = record

    def query_one(self, _sel, _cls=None):
        return self._SP()


def _plasmid(total=300, circular=True):
    rng = __import__("random").Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(total))
    rec = SeqRecord(Seq(seq), id="P", name="P",
                    annotations={"molecule_type": "DNA",
                                 "topology": ("circular" if circular
                                              else "linear")})
    rec.features = [
        SeqFeature(FeatureLocation(40, 90, 1), type="CDS",
                   qualifiers={"label": ["plain"]}),
        SeqFeature(CompoundLocation([FeatureLocation(120, 150, 1),
                                     FeatureLocation(200, 230, 1)]),
                   type="CDS", qualifiers={"label": ["spliced"]}),
        SeqFeature(CompoundLocation([FeatureLocation(total - 20, total, 1),
                                     FeatureLocation(0, 10, 1)]),
                   type="misc_feature", qualifiers={"label": ["wrap"]}),
    ]
    return rec


def _bases(rec):
    """{label: extracted bases} via Biopython — the outside oracle."""
    return {(f.qualifiers.get("label") or ["?"])[0]:
            str(f.extract(rec.seq)).upper()
            for f in rec.features if f.type != "source"}


@pytest.mark.parametrize("bp", [1, 45, 125, 205, 299])
def test_set_origin_conserves_every_feature(bp):
    rec = _plasmid()
    before = _bases(rec)
    app = _ReframeApp(rec)
    res = _call("set-origin", app, bp=bp)
    assert _status(res) == 200, res
    assert _body(res)["rotated"] is True
    after = _bases(app._current_record)
    assert after == before, f"rotation at {bp} changed feature bases"


def test_set_origin_by_feature_name():
    rec = _plasmid()
    before = _bases(rec)
    app = _ReframeApp(rec)
    res = _call("set-origin", app, feature="spliced")
    assert _status(res) == 200, res
    assert _body(res)["offset"] == 120
    assert _bases(app._current_record) == before


def test_set_origin_by_feature_idx():
    app = _ReframeApp(_plasmid())
    res = _call("set-origin", app, feature_idx=0)
    assert _status(res) == 200, res
    assert _body(res)["offset"] == 40


def test_set_origin_refuses_a_linear_record():
    """Rotating a linear molecule fabricates a junction between two ends that
    are not joined — the GUI refuses it and so must the side door."""
    app = _ReframeApp(_plasmid(circular=False))
    res = _call("set-origin", app, bp=5)
    assert _status(res) == 422
    assert "circular" in _body(res)["error"].lower()


@pytest.mark.parametrize("kw,code", [
    ({}, 400),                                  # neither
    ({"bp": 5, "feature": "plain"}, 400),       # both
    ({"bp": 99999}, 400),                       # out of range
    ({"bp": -1}, 400),
    ({"feature": "nope"}, 404),                 # unknown
    ({"feature_idx": 99}, 400),                 # out of range
])
def test_set_origin_refusals(kw, code):
    app = _ReframeApp(_plasmid())
    assert _status(_call("set-origin", app, **kw)) == code


def test_set_origin_at_the_existing_origin_is_a_no_op():
    app = _ReframeApp(_plasmid())
    res = _call("set-origin", app, bp=0)
    assert _status(res) == 200
    assert _body(res)["rotated"] is False


def test_set_origin_ambiguous_feature_name_is_409():
    rec = _plasmid()
    rec.features.append(SeqFeature(FeatureLocation(250, 260, 1), type="CDS",
                                    qualifiers={"label": ["plain"]}))
    app = _ReframeApp(rec)
    res = _call("set-origin", app, feature="plain")
    assert _status(res) == 409
    assert "feature_idx" in _body(res)["error"]


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


@pytest.mark.parametrize("circular", [True, False])
def test_reverse_complement_conserves_every_feature(circular):
    """`extract()` strand-corrects, so after flipping the whole molecule each
    feature must still extract the SAME bases."""
    rec = _plasmid(circular=circular) if circular else _plasmid(
        circular=False)
    if not circular:
        # a linear record can't carry the origin-wrapping fixture feature
        rec.features = [f for f in rec.features
                        if (f.qualifiers.get("label") or [""])[0] != "wrap"]
    before = _bases(rec)
    app = _ReframeApp(rec)
    res = _call("reverse-complement", app)
    assert _status(res) == 200, res
    assert _bases(app._current_record) == before


def test_reverse_complement_works_on_a_linear_record():
    """Unlike set-origin: flipping is meaningful for both topologies."""
    app = _ReframeApp(_plasmid(circular=False))
    assert _status(_call("reverse-complement", app)) == 200


def test_reframes_refuse_with_nothing_loaded():
    for name in ("set-origin", "reverse-complement"):
        assert _status(_call(name, _ReframeApp(None), bp=1)) == 422


# ═══════════════════════════════════════════════════════════════════════════
# taxon-search / read-trace — validation + the egress gate
# ═══════════════════════════════════════════════════════════════════════════

def test_taxon_search_is_egress_gated(monkeypatch):
    """Gated like every other `*-search`: a query string an autonomous agent
    composed is exactly what the gate holds."""
    monkeypatch.setattr(sca, "_online_lookups_armed", lambda: False)
    res = _call("taxon-search", None, query="Saccharomyces cerevisiae")
    assert _status(res) == 403


def test_taxon_search_refuses_a_bare_taxid(monkeypatch):
    monkeypatch.setattr(sca, "_online_lookups_armed", lambda: True)
    res = _call("taxon-search", None, query="4102")
    assert _status(res) == 400
    assert "already a taxid" in _body(res)["error"]


def test_taxon_search_shape(monkeypatch):
    monkeypatch.setattr(sca, "_online_lookups_armed", lambda: True)
    monkeypatch.setattr(
        sca, "_ncbi_taxid_search",
        lambda q, retmax=200: ([{"taxid": "4932",
                                 "name": "Saccharomyces cerevisiae"}],
                               1, "saccharomyces cerevisiae[All Names]"))
    res = _call("taxon-search", None, query="saccharomyces cerevisiae")
    assert _status(res) == 200
    b = _body(res)
    assert b["hits"][0]["taxid"] == "4932"
    assert b["total_matches"] == 1
    assert b["term_used"]


@pytest.mark.parametrize("kw,code", [
    ({}, 400),
    ({"path": "/tmp/x.fastq"}, 400),       # wrong extension
])
def test_read_trace_validation(kw, code):
    assert _status(_call("read-trace", None, **kw)) == code


def test_read_trace_missing_file_is_404(tmp_path):
    res = _call("read-trace", None, path=str(tmp_path / "nope.ab1"))
    assert _status(res) == 404


# ═══════════════════════════════════════════════════════════════════════════
# Protein collections — the whole store
# ═══════════════════════════════════════════════════════════════════════════

PROT = "MAEAKGFERLDLKAWQ"


def test_protein_collection_round_trip():
    assert _status(_call("create-protein-collection", None, name="Panel")) == 200
    res = _call("save-protein", None, collection="Panel", name="ProtX",
                sequence=PROT, description="test")
    assert _status(res) == 200, res
    assert _body(res)["overwrote"] is False

    got = _call("get-protein-collection", None, name="Panel")
    assert _status(got) == 200
    prots = _body(got)["proteins"]
    assert len(prots) == 1
    assert prots[0]["sequence"] == PROT
    assert prots[0]["length"] == len(PROT)

    lst = _body(_call("list-protein-collections", None))
    assert any(c["name"] == "Panel" and c["n_proteins"] == 1
               for c in lst["collections"])


def test_save_protein_upserts_by_name():
    _call("create-protein-collection", None, name="P2")
    _call("save-protein", None, collection="P2", name="x", sequence=PROT)
    res = _call("save-protein", None, collection="P2", name="x",
                sequence=PROT + "GG")
    assert _body(res)["overwrote"] is True
    prots = _body(_call("get-protein-collection", None, name="P2"))["proteins"]
    assert len(prots) == 1, "upsert must not duplicate"
    assert prots[0]["sequence"] == PROT + "GG"


def test_save_protein_creates_the_collection():
    res = _call("save-protein", None, collection="Fresh", name="p",
                sequence=PROT)
    assert _status(res) == 200
    assert _status(_call("get-protein-collection", None, name="Fresh")) == 200


def test_save_protein_refuses_non_canonical_residues():
    res = _call("save-protein", None, collection="P3", name="oops",
                sequence="MAEJJJOUZZ@@")
    assert _status(res) == 400, res


def test_save_protein_warns_but_accepts_pure_acgt():
    """A/C/G/T are all canonical residues (Ala/Cys/Gly/Thr), so a DNA string
    is indistinguishable from a real peptide of those four — refusing it
    would reject legitimate proteins. It is accepted and FLAGGED instead."""
    res = _call("save-protein", None, collection="P4", name="dna?",
                sequence="ACGTACGTACGTACGTACGT")
    assert _status(res) == 200, res
    assert _body(res)["warnings"], "a pure-ACGT 'protein' must be flagged"
    assert "amino acids" in _body(res)["warnings"][0]


def test_save_protein_does_not_warn_on_a_real_protein():
    res = _call("save-protein", None, collection="P5", name="real",
                sequence=PROT)
    assert "warnings" not in _body(res)


def test_delete_protein_collection_needs_confirm():
    _call("create-protein-collection", None, name="Doomed")
    _call("save-protein", None, collection="Doomed", name="a", sequence=PROT)
    res = _call("delete-protein-collection", None, name="Doomed")
    assert _status(res) == 400
    assert "confirm" in _body(res)["error"]
    # still there
    assert _status(_call("get-protein-collection", None, name="Doomed")) == 200
    ok = _call("delete-protein-collection", None, name="Doomed", confirm=True)
    assert _status(ok) == 200
    assert _status(_call("get-protein-collection", None, name="Doomed")) == 404


def test_protein_collection_rename_and_clash():
    _call("create-protein-collection", None, name="A")
    _call("create-protein-collection", None, name="B")
    assert _status(_call("rename-protein-collection", None,
                         name="A", new_name="C")) == 200
    assert _status(_call("rename-protein-collection", None,
                         name="C", new_name="B")) == 409
    assert _status(_call("rename-protein-collection", None,
                         name="ghost", new_name="D")) == 404


def test_delete_protein_distinguishes_the_two_404s():
    _call("create-protein-collection", None, name="K")
    _call("save-protein", None, collection="K", name="keep", sequence=PROT)
    miss_coll = _call("delete-protein", None, collection="ghost", name="keep")
    assert _status(miss_coll) == 404
    assert "collection" in _body(miss_coll)["error"]
    miss_prot = _call("delete-protein", None, collection="K", name="ghost")
    assert _status(miss_prot) == 404
    assert _status(_call("delete-protein", None,
                         collection="K", name="keep")) == 200


def test_missing_protein_collection_is_404_not_empty():
    """An empty list would read as "it exists and is empty"."""
    assert _status(_call("get-protein-collection", None, name="nope")) == 404


# ═══════════════════════════════════════════════════════════════════════════
# codon_forbidden_enzymes — the setting the two surfaces disagreed on
# ═══════════════════════════════════════════════════════════════════════════

def test_codon_forbidden_enzymes_is_settable():
    res = _call("set-setting", None, key="codon_forbidden_enzymes",
                value=["BsaI", "BsmBI"])
    assert _status(res) == 200, res
    got = _body(_call("get-settings", None))["settings"]
    assert "codon_forbidden_enzymes" in got
    assert got["codon_forbidden_enzymes"]["value"] == ["BsaI", "BsmBI"]


def test_codon_forbidden_enzymes_refuses_an_unknown_name():
    """REFUSED, not dropped: a forbidden site silently not scrubbed ends up in
    DNA the user orders. The opposite choice from the display-filter CSV."""
    res = _call("set-setting", None, key="codon_forbidden_enzymes",
                value=["BsaI", "NotAnEnzyme"])
    assert _status(res) == 400
    assert "NotAnEnzyme" in _body(res)["error"]


def test_codon_forbidden_enzymes_accepts_csv_and_dedupes():
    res = _call("set-setting", None, key="codon_forbidden_enzymes",
                value="BsaI, BsaI , BsmBI")
    assert _status(res) == 200, res
    assert _body(_call("get-settings", None))["settings"][
        "codon_forbidden_enzymes"]["value"] == ["BsaI", "BsmBI"]


# ═══════════════════════════════════════════════════════════════════════════
# tools --filter — budgeted discovery
# ═══════════════════════════════════════════════════════════════════════════

def test_tools_unfiltered_keeps_its_bare_shape():
    """Unfiltered `/tools` is served ahead of the dispatcher's envelope and
    has always been exactly `{endpoints: [...]}`. The filter's extra keys must
    NOT leak onto it — a client may rely on that shape."""
    b = _body(_call("tools", None))
    assert set(b) == {"endpoints"}
    assert len(b["endpoints"]) == len(sc._state._AGENT_HANDLERS)


def test_tools_filter_ranks_and_reports_an_honest_total():
    b = _body(_call("tools", None, filter="primer"))
    assert 0 < b["count"] <= 25
    assert b["total_matches"] >= b["count"]
    if b["total_matches"] > b["count"]:
        assert b["capped"] is True, "a capped page must say so"
    # a name hit outranks a doc-only hit
    assert any("primer" in e["name"] for e in b["endpoints"][:3])


def test_tools_filter_brief_drops_the_doc_corpus():
    b = _body(_call("tools", None, filter="primer", brief=True))
    assert b["endpoints"]
    assert all("doc_full" not in e for e in b["endpoints"])


def test_tools_filter_no_match_is_zero_not_everything():
    b = _body(_call("tools", None, filter="zzzznotathing"))
    assert b["count"] == 0 and b["total_matches"] == 0


def test_tools_filter_finds_the_new_endpoints():
    """The point of the filter: searching a VERB finds the capability."""
    for term, want in (("origin", "set-origin"),
                       ("taxid", "taxon-search"),
                       ("trace", "read-trace"),
                       ("amplicon", "save-amplicon")):
        names = [e["name"] for e in
                 _body(_call("tools", None, filter=term, brief=True))[
                     "endpoints"]]
        assert want in names, f"filter {term!r} did not surface {want}"


# ═══════════════════════════════════════════════════════════════════════════
# import-primers / update-stored-primer
# ═══════════════════════════════════════════════════════════════════════════

_P1 = "ACGTACGTACGTACGTACGT"
_P2 = "TTGACCGGATCCAATTCCGG"


def _csv(tmp_path, rows) -> str:
    p = tmp_path / "primers.csv"
    p.write_text("Name,Sequence\n" + "".join(f"{n},{s}\n" for n, s in rows))
    return str(p)


def test_import_primers_round_trip(tmp_path):
    path = _csv(tmp_path, [("fwd", _P1), ("rev", _P2)])
    res = _call("import-primers", None, path=path)
    assert _status(res) == 200, res
    assert _body(res)["imported"] == 2
    names = {p["name"] for p in _body(_call("list-primers", None))["primers"]}
    assert {"fwd", "rev"} <= names


def test_import_primers_reports_skipped_rows(tmp_path):
    path = _csv(tmp_path, [("good", _P1), ("bad", "NOTDNA!!"), ("empty", "")])
    res = _call("import-primers", None, path=path)
    assert _status(res) == 200, res
    assert _body(res)["imported"] == 1
    assert _body(res)["skipped"], "an invalid row must be REPORTED, not silent"


def test_import_primers_duplicate_modes(tmp_path):
    path = _csv(tmp_path, [("fwd", _P1)])
    assert _body(_call("import-primers", None, path=path))["imported"] == 1
    again = _body(_call("import-primers", None, path=path))
    assert again["imported"] == 0 and again["duplicates"] == 1
    path2 = _csv(tmp_path, [("renamed", _P1)])
    rep = _body(_call("import-primers", None, path=path2,
                      on_duplicate="replace"))
    assert rep["imported"] == 1
    names = {p["name"] for p in _body(_call("list-primers", None))["primers"]}
    assert "renamed" in names and "fwd" not in names


@pytest.mark.parametrize("kw,code", [
    ({}, 400),
    ({"path": "/tmp/x.txt"}, 400),
    ({"path": "/tmp/definitely-not-here-12345.csv"}, 404),
])
def test_import_primers_validation(kw, code):
    assert _status(_call("import-primers", None, **kw)) == code


def test_import_primers_unknown_collection_is_404(tmp_path):
    path = _csv(tmp_path, [("fwd", _P1)])
    res = _call("import-primers", None, path=path, collection="ghost")
    assert _status(res) == 404


def _seed_primer():
    r = _call("create-primer", None, name="seed", sequence=_P1)
    assert _status(r) == 200, r


def test_update_stored_primer_sets_status_and_name():
    _seed_primer()
    res = _call("update-stored-primer", None, sequence=_P1,
                new_name="renamed", status="Ordered")
    assert _status(res) == 200, res
    assert set(_body(res)["changed"]) == {"name", "status"}
    p = next(p for p in _body(_call("list-primers", None))["primers"]
             if p["sequence"] == _P1)
    assert p["name"] == "renamed" and p["status"] == "Ordered"


def test_update_stored_primer_rejects_a_bad_status():
    _seed_primer()
    assert _status(_call("update-stored-primer", None,
                         sequence=_P1, status="Shipped")) == 400


def test_update_stored_primer_noop_reports_empty_changed():
    """A no-op must be visible as `changed: []`, not read as a successful
    edit."""
    _seed_primer()
    res = _call("update-stored-primer", None, sequence=_P1, new_name="seed")
    assert _status(res) == 200
    assert _body(res)["changed"] == []


def test_update_stored_primer_404_and_identity():
    assert _status(_call("update-stored-primer", None,
                         sequence="GGGGGGGGGGGGGGGGGGGG",
                         status="Ordered")) == 404
    assert _status(_call("update-stored-primer", None,
                         status="Ordered")) == 400     # no identifier


def test_cart_flag_and_cart_export(tmp_path):
    _seed_primer()
    res = _call("update-stored-primer", None, sequence=_P1, in_cart=True)
    assert _status(res) == 200, res
    assert "in_cart" in _body(res)["changed"]
    out = tmp_path / "cart.csv"
    exp = _call("export-primers", None, path=str(out), cart=True,
                clear_cart_after=True)
    assert _status(exp) == 200, exp
    assert out.exists()
    assert _body(exp)["cart_cleared"] >= 1
    # cart is now empty, so a second cart export has nothing to write
    out2 = tmp_path / "cart2.csv"
    assert _status(_call("export-primers", None, path=str(out2),
                         cart=True)) == 404


def test_cart_export_with_an_empty_cart_is_404(tmp_path):
    _seed_primer()
    res = _call("export-primers", None, path=str(tmp_path / "c.csv"),
                cart=True)
    assert _status(res) == 404
    assert "cart is empty" in _body(res)["error"]


# ═══════════════════════════════════════════════════════════════════════════
# Stored alignments — the writer, and the asymmetry it closes
# ═══════════════════════════════════════════════════════════════════════════

_PLASMID = ("ATGGCTGAAGCTAAGGGTTTCGAACGTCTGGATCTGAAAGCTTGG"
            "ACGTACGTACGTACGTAAGGCCTTAACCGGTTAAGGCCTTAACCG"
            "GTTAAGGCCTTAACCGGTTAACCGGTTAAGGCCTTGACCGGATCC")


def _seed_entry(name="Target", coll=None):
    """One real library entry, present in BOTH places the alignment
    endpoints read: the active-collection library mirror (which
    `_iter_library_readonly` walks to resolve the target) and
    `collections.json` (which the persistence helper walks to find the
    collection that holds it). Seeded directly rather than through
    `create-collection`, which is app-coupled."""
    rec = SeqRecord(Seq(_PLASMID), id=name[:8], name=name[:8],
                    annotations={"molecule_type": "DNA",
                                 "topology": "circular"})
    entry = {"id": name, "name": name,
             "gb_text": sc._record_to_gb_text(rec)}
    coll = coll or (sc._get_active_collection_name() or "Default")
    colls = [c for c in sc._load_collections() if isinstance(c, dict)]
    tgt = next((c for c in colls if str(c.get("name") or "") == coll), None)
    if tgt is None:
        tgt = {"name": coll, "plasmids": []}
        colls.append(tgt)
    plasmids = [p for p in (tgt.get("plasmids") or []) if isinstance(p, dict)]
    plasmids = [p for p in plasmids if str(p.get("id") or "") != entry["id"]]
    plasmids.append(dict(entry))
    tgt["plasmids"] = plasmids
    sc._save_collections(colls)
    sc._save_library(plasmids)
    return entry


def test_store_alignment_persists_and_read_consensus_can_use_it():
    """THE point of the writer: `read-consensus` reads stored alignments, and
    nothing but the GUI could write them — so an agent could never set up
    that endpoint's own input. One call must now close that loop."""
    _seed_entry()
    before = _body(_call("read-consensus", None, name="Target"))
    assert before.get("n_reads", 0) == 0, "fixture should start with no reads"

    read = _PLASMID[10:70]                      # a clean 60 bp "read"
    res = _call("store-alignment", None, name="Target", read_seq=read,
                read_name="sanger-1")
    assert _status(res) == 200, res
    # A PERFECT short read scores low on raw identity by construction — the
    # plasmid it doesn't reach counts as gaps. `covered_identity_pct` is the
    # number to judge on ([INV-207]).
    b = _body(res)
    assert b["covered_identity_pct"] == 100.0, b
    assert b["identity_pct"] < 100.0, "raw identity should include the gaps"
    assert 0 < b["covered_pct"] <= 100

    lst = _body(_call("list-alignments", None, name="Target"))
    assert lst["count"] == 1
    assert lst["alignments"][0]["read_name"] == "sanger-1"

    after = _body(_call("read-consensus", None, name="Target"))
    assert after["n_reads"] == 1, "read-consensus still can't see it"


def test_store_alignment_is_idempotent():
    """Re-storing a byte-identical alignment updates the row instead of
    appending a duplicate (`_merge_stored_alignments` semantics)."""
    _seed_entry(name="Idem")
    read = _PLASMID[5:65]
    for _ in range(3):
        assert _status(_call("store-alignment", None, name="Idem",
                             read_seq=read, read_name="r")) == 200
    assert _body(_call("list-alignments", None, name="Idem"))["count"] == 1


def test_store_alignment_computes_it_rather_than_trusting_the_caller():
    """There is deliberately no way to hand in gapped strings or
    coordinates — a caller-built alignment is a caller-built VERDICT, and
    `read-consensus` reads it as ground truth."""
    import inspect
    src = inspect.getsource(sc._h_store_alignment)
    assert "_pairwise_align(" in src, "the engine must do the alignment"
    for forbidden in ('payload.get("aligned_q")', 'payload.get("aligned_t")',
                      'payload.get("result")', 'payload.get("segments")'):
        assert forbidden not in src, f"store-alignment must not accept {forbidden}"


@pytest.mark.parametrize("kw,code", [
    ({}, 400),                                       # no target
    ({"name": "ghost", "read_seq": "ACGTACGTAC"}, 404),
])
def test_store_alignment_refusals(kw, code):
    _seed_entry()
    assert _status(_call("store-alignment", None, **kw)) == code


def test_store_alignment_refuses_an_empty_read():
    _seed_entry()
    assert _status(_call("store-alignment", None, name="Target",
                         read_seq="")) == 400


def test_store_alignment_caps_an_oversized_read():
    _seed_entry()
    res = _call("store-alignment", None, name="Target",
                read_seq="ACGT" * 100_000)          # 400 kb
    assert _status(res) == 413


def test_update_alignment_toggles_visibility():
    _seed_entry(name="Vis")
    _call("store-alignment", None, name="Vis", read_seq=_PLASMID[:60],
          read_name="r")
    sid = _body(_call("list-alignments", None,
                      name="Vis"))["alignments"][0]["stored_id"]
    res = _call("update-alignment", None, name="Vis", stored_id=sid,
                visible=False)
    assert _status(res) == 200, res
    assert "visible" in _body(res)["changed"]
    lst = _body(_call("list-alignments", None, name="Vis"))
    assert lst["alignments"][0]["visible"] is False
    assert lst["count"] == 1, "a hidden alignment is still ON DISK and listed"
    assert lst["n_visible"] == 0
    # a no-op repeat reports no change rather than claiming an edit
    again = _call("update-alignment", None, name="Vis", stored_id=sid,
                  visible=False)
    assert _body(again)["changed"] == []


def test_update_alignment_404_on_unknown_stored_id():
    _seed_entry(name="U404")
    _call("store-alignment", None, name="U404", read_seq=_PLASMID[:60])
    res = _call("update-alignment", None, name="U404",
                stored_id="nope", visible=False)
    assert _status(res) == 404
    assert _body(res)["stored_ids"], "the error should name the real ids"


def test_delete_alignment_by_id_and_all():
    _seed_entry(name="Del")
    for i in range(3):
        _call("store-alignment", None, name="Del",
              read_seq=_PLASMID[i:60 + i], read_name=f"r{i}")
    lst = _body(_call("list-alignments", None, name="Del"))
    assert lst["count"] == 3
    sid = lst["alignments"][0]["stored_id"]
    one = _call("delete-alignment", None, name="Del", stored_id=sid)
    assert _status(one) == 200 and _body(one)["deleted"] == 1
    assert _body(_call("list-alignments", None, name="Del"))["count"] == 2
    allr = _call("delete-alignment", None, name="Del", all=True)
    assert _status(allr) == 200 and _body(allr)["deleted"] == 2
    assert _body(_call("list-alignments", None, name="Del"))["count"] == 0


def test_delete_alignment_requires_an_explicit_selector():
    """`all` is spelled out rather than reachable by omitting stored_id —
    dropping every alignment resets a verification."""
    _seed_entry(name="Sel")
    _call("store-alignment", None, name="Sel", read_seq=_PLASMID[:60])
    assert _status(_call("delete-alignment", None, name="Sel")) == 400
    assert _status(_call("delete-alignment", None, name="Sel",
                         stored_id="x", all=True)) == 400


# ═══════════════════════════════════════════════════════════════════════════
# restore-deleted-plasmid — the agent can now consume the undo it creates
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture
def undo_rig(monkeypatch):
    """A real `PlasmidApp` + real `LibraryPanel` with ONLY the UI hooks
    stubbed, so `_push_library_undo` (depth cap, drop-incomplete-batches) and
    `action_undo_delete` (slot restore, collection guard, batch all-or-
    nothing) all run for real. Stubbing either would hide the exact
    behaviour these tests exist to check."""
    lib = [{"id": f"P{i}", "name": f"Plasmid {i}", "gb_text": ""}
           for i in range(4)]
    coll = sc._get_active_collection_name() or "Default"
    sc._save_collections([{"name": coll, "plasmids": [dict(e) for e in lib]}])
    # A fresh sandbox has no ACTIVE collection, and the restore's collection
    # guard (the UI's own, `want and want != active`) then refuses every
    # restore. Set it so the rig matches a real session.
    sc._set_active_collection_name(coll)
    sc._save_library([dict(e) for e in lib])

    panel = sc.LibraryPanel.__new__(sc.LibraryPanel)
    panel._marked_ids, panel._mark_kind = set(), {}
    panel._view_mode = "plasmids"
    box = {"notes": [], "saved": None}

    app = sc.PlasmidApp.__new__(sc.PlasmidApp)
    app._current_record = None
    app.notify = lambda m, severity="information", **k: box["notes"].append(m)
    app.call_from_thread = lambda fn, *a, **k: fn(*a, **k)
    app.query_one = lambda *a, **k: panel
    monkeypatch.setattr(type(panel), "app", property(lambda self: app))
    monkeypatch.setattr(type(panel), "_repopulate_plasmids", lambda self: None)
    monkeypatch.setattr(type(panel), "set_active", lambda self, x: None)
    monkeypatch.setattr(
        type(panel), "_delete_save_to_disk",
        lambda self, entries, bin_, op="": sc._save_library(entries))
    return app, panel, [e["id"] for e in lib], box, coll


def _delete_one(app, idx, ids, coll):
    """Record one deletion the way the delete paths do, then remove it."""
    lib = sc._load_library()
    entry = next(e for e in lib if e["id"] == ids[idx])
    app._push_library_undo({"name": entry["name"], "entry": dict(entry),
                            "index": idx, "collection": coll})
    sc._save_library([e for e in lib if e["id"] != ids[idx]])


def test_restore_deleted_plasmid_puts_it_back_in_its_slot(undo_rig):
    app, panel, ids, box, coll = undo_rig
    _delete_one(app, 1, ids, coll)
    assert [e["id"] for e in sc._load_library()] == [ids[0], ids[2], ids[3]]
    res = _call("restore-deleted-plasmid", app)
    assert _status(res) == 200, res
    assert _body(res)["count"] == 1
    assert _body(res)["restored"] == ["Plasmid 1"]
    # Back in the SLOT it came from, not appended at the end.
    assert [e["id"] for e in sc._load_library()] == ids


def test_restore_with_an_empty_stack_is_404(undo_rig):
    app, panel, ids, box, coll = undo_rig
    res = _call("restore-deleted-plasmid", app)
    assert _status(res) == 404
    assert _body(res)["restored"] == []


def test_restore_refuses_when_the_collection_changed(undo_rig):
    """The library file IS the active collection's contents, so restoring
    into a different one would file the plasmid somewhere it never lived.
    The records must SURVIVE the refusal."""
    app, panel, ids, box, coll = undo_rig
    _delete_one(app, 0, ids, coll)
    stack = app._library_undo_stack()
    stack[-1]["collection"] = "SomeOtherCollection"
    res = _call("restore-deleted-plasmid", app)
    assert _status(res) == 409
    assert _body(res)["deleted_from"] == "SomeOtherCollection"
    assert len(app._library_undo_stack()) == 1, \
        "a refused restore must leave the record recoverable"


def test_restore_brings_a_whole_bulk_batch_back_in_one_call(undo_rig):
    """Restoring 2 of 3 and reporting success is the partial-recovery-reads-
    as-complete failure the batch design exists to prevent."""
    app, panel, ids, box, coll = undo_rig
    lib = sc._load_library()
    # Records are pushed highest-index-first, all sharing one batch id.
    for idx in (2, 1, 0):
        entry = next(e for e in lib if e["id"] == ids[idx])
        app._push_library_undo({"name": entry["name"], "entry": dict(entry),
                                "index": idx, "collection": coll,
                                "batch": "bulk-1"})
    sc._save_library([e for e in lib if e["id"] not in set(ids[:3])])
    assert [e["id"] for e in sc._load_library()] == [ids[3]]
    res = _call("restore-deleted-plasmid", app)
    assert _status(res) == 200, res
    assert _body(res)["count"] == 3, "one call must restore the WHOLE batch"
    assert [e["id"] for e in sc._load_library()] == ids
    assert app._library_undo_stack() == []


def test_master_delete_has_no_endpoint():
    """Authorised 2026-10-16: agents may delete and undo. The WIPE stays
    human-only — a deliberately unreachable capability, not an oversight."""
    names = set(sc._state._AGENT_HANDLERS)
    for forbidden in ("master-delete", "wipe", "wipe-all",
                      "delete-everything", "master-wipe"):
        assert forbidden not in names
    assert not [n for n in names if "master" in n and "delete" in n]


# ═══════════════════════════════════════════════════════════════════════════
# store-alignment hardening: what counts as a read at all
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("n_bp", [1, 5, 19])
def test_a_read_too_short_to_be_evidence_is_refused(n_bp):
    """`read-consensus` counts every stored row as a read, so a 1 bp "read"
    would raise depth at some base and could flip a verdict while supporting
    nothing. It used to be stored at 0% identity."""
    _seed_entry(name="Short")
    res = _call("store-alignment", None, name="Short",
                read_seq=_PLASMID[:n_bp])
    assert _status(res) == 400
    assert "cannot support or contradict" in _body(res)["error"]


def test_a_mostly_no_call_read_is_refused():
    """An ambiguity-dominated read aligns at ~0% and would be filed as a
    FAILED VERIFICATION rather than a bad read — the no-call-as-evidence
    shape [INV-200] removed from the heterogeneity caller."""
    _seed_entry(name="Ns")
    res = _call("store-alignment", None, name="Ns", read_seq="N" * 60)
    assert _status(res) == 400
    assert "no-calls" in _body(res)["error"]


def test_a_read_longer_than_its_template_is_flagged():
    """`covered_pct` saturates at 100, so without a warning a 3x concatemer
    reads as "perfect and fully covered"."""
    _seed_entry(name="Cat")
    res = _call("store-alignment", None, name="Cat",
                read_seq=_PLASMID * 3)
    assert _status(res) == 200, res
    b = _body(res)
    assert b["warnings"], "a read longer than the plasmid must say so"
    assert "concatemer" in b["warnings"][0]


def test_a_normal_read_carries_no_warnings():
    """The guards must not fire on ordinary input, or they become noise."""
    _seed_entry(name="Fine")
    b = _body(_call("store-alignment", None, name="Fine",
                    read_seq=_PLASMID[10:80]))
    assert b["warnings"] == []
    assert b["covered_identity_pct"] == 100.0
