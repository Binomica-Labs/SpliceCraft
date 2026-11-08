"""Minimum-viable coverage for the agent endpoints added this
session: grammar / primer / parts / feature-library CRUD. Each
test exercises the happy path + at least one error class so a
regression in payload validation surfaces in CI."""

from __future__ import annotations

import pytest

import splicecraft as sc


class _MockApp:
    """Minimal app stub for testing agent endpoints that route writes
    through `app.call_from_thread`. The real PlasmidApp uses the call
    to hop to the Textual UI thread; in tests we run inline."""

    def call_from_thread(self, fn):
        return fn()


@pytest.fixture
def mock_app():
    return _MockApp()


# ── Cloning grammar CRUD ──────────────────────────────────────────────────────


def _make_grammar(gid: str = "test_grammar") -> dict:
    return {
        "id":              gid,
        "name":            "Test Grammar",
        "enzyme":          "BsaI",
        "level_up_enzyme": "BpiI",
        "site":            "GGTCTC",
        "spacer":          "N",
        "pad":             "GCGC",
        "forbidden_sites": {"BsaI": "GGTCTC"},
        "positions": [
            {"name": "Pos 1", "type": "Promoter",
             "oh5": "GGAG", "oh3": "AATG", "color": "green"},
            {"name": "Pos 2", "type": "CDS",
             "oh5": "AATG", "oh3": "GCTT", "color": "yellow"},
        ],
        "coding_types": ["CDS"],
        "type_to_insdc": {"CDS": "CDS", "Promoter": "promoter"},
    }


def test_grammar_create_get_delete_roundtrip():
    g = _make_grammar()
    res = sc._h_create_grammar(None, g)
    assert res == {"ok": True, "grammar_id": "test_grammar"}

    got = sc._h_get_grammar(None, {"grammar_id": "test_grammar"})
    assert got["ok"]
    assert got["grammar"]["name"] == "Test Grammar"
    assert got["grammar"]["editable"]

    listed = sc._h_list_grammars(None, {})
    ids = {row["id"] for row in listed["grammars"]}
    assert "test_grammar" in ids
    assert "gb_l0" in ids   # built-ins present

    deleted = sc._h_delete_grammar(None, {"grammar_id": "test_grammar"})
    assert deleted == {"ok": True, "grammar_id": "test_grammar"}


def test_grammar_create_refuses_builtin_id():
    g = _make_grammar(gid="gb_l0")
    res = sc._h_create_grammar(None, g)
    assert isinstance(res, tuple) and res[1] == 400
    assert "built-in" in res[0]["error"]


def test_grammar_update_unknown_returns_404():
    g = _make_grammar(gid="nonexistent_grammar_xyz")
    res = sc._h_update_grammar(None, g)
    assert isinstance(res, tuple) and res[1] == 404


def test_grammar_delete_builtin_refused():
    res = sc._h_delete_grammar(None, {"grammar_id": "gb_l0"})
    assert isinstance(res, tuple) and res[1] == 400


# ── Primer CRUD ───────────────────────────────────────────────────────────────


def test_primer_create_list_get_delete(mock_app):
    res = sc._h_create_primer(mock_app, {
        "name":     "test_primer_fwd",
        "sequence": "ATGCATGCATGCATGCATGC",
        "tm":       62.5,
        "status":   "Designed",
    })
    assert res["ok"]
    assert res["name"] == "test_primer_fwd"

    listed = sc._h_list_primers(mock_app, {})
    assert listed["count"] >= 1

    got = sc._h_get_primer(mock_app, {"sequence": "ATGCATGCATGCATGCATGC"})
    assert got["ok"]
    assert got["primer"]["name"] == "test_primer_fwd"
    assert got["primer"]["tm"] == 62.5

    res = sc._h_delete_primer(mock_app, {
        "sequence": "ATGCATGCATGCATGCATGC",
    })
    assert res["ok"]


def test_primer_create_duplicate_returns_409(mock_app):
    sc._h_create_primer(mock_app, {
        "name": "p1", "sequence": "ACGTACGTACGTACGTACGT",
    })
    res = sc._h_create_primer(mock_app, {
        "name": "p2", "sequence": "ACGTACGTACGTACGTACGT",
    })
    assert isinstance(res, tuple) and res[1] == 409
    assert "existing_name" in res[0]


def test_primer_get_missing_returns_404():
    res = sc._h_get_primer(None, {"sequence": "NONEXISTENT_SEQ_ATGC"})
    assert isinstance(res, tuple) and res[1] == 404


# ── Parts bin CRUD ─────────────────────────────────────────────────────────────


def test_part_create_update_delete_roundtrip(mock_app):
    part = {
        "name":     "test_part_cds",
        "grammar":  "gb_l0",
        "type":     "CDS",
        "position": "PC",
        "level":    0,
        "oh5":      "AATG",
        "oh3":      "GCTT",
        "sequence": "ATGGCAACG" * 5,
    }
    res = sc._h_create_part(mock_app, part)
    assert res["ok"]
    assert res["name"] == "test_part_cds"

    got = sc._h_get_part(mock_app, {"name": "test_part_cds"})
    assert got["ok"]
    assert got["part"]["grammar"] == "gb_l0"

    # Update
    part["oh5"] = "GGAG"
    res = sc._h_update_part(mock_app, part)
    assert res["ok"]

    res = sc._h_delete_part(mock_app, {
        "name": "test_part_cds", "grammar": "gb_l0",
    })
    assert res["ok"]
    assert res["deleted"] == 1


def test_part_create_with_vector_and_primer_fields(mock_app):
    """Parameter parity with PartEditModal's Edit mode: an agent can set
    backbone, selection marker, and forward/reverse domestication primers
    (the Tms are derived on save, mirroring the modal)."""
    part = {
        "name": "parity_part", "grammar": "gb_l0", "type": "CDS",
        "level": 0, "sequence": "ATGGCAACG" * 5,
        "backbone": "pUPD2", "marker": "Cm",
        "fwd_primer": "ACGTACGTACGTACGT", "rev_primer": "TGCATGCATGCATGCA",
    }
    assert sc._h_create_part(mock_app, part)["ok"]
    p = sc._h_get_part(mock_app, {"name": "parity_part"})["part"]
    assert p["backbone"] == "pUPD2" and p["marker"] == "Cm"
    assert p["fwd_primer"] == "ACGTACGTACGTACGT"
    assert p["rev_primer"] == "TGCATGCATGCATGCA"
    assert isinstance(p["fwd_tm"], float) and p["fwd_tm"] > 0
    assert isinstance(p["rev_tm"], float) and p["rev_tm"] > 0


def test_part_rejects_non_iupac_primer(mock_app):
    res = sc._h_create_part(mock_app, {
        "name": "badprimer", "grammar": "gb_l0", "type": "CDS", "level": 0,
        "sequence": "ATGGCA" * 5, "fwd_primer": "ACGTXYZ"})
    assert isinstance(res, tuple) and res[1] == 400


def test_part_create_duplicate_returns_409(mock_app):
    part = {
        "name":     "dup_part",
        "grammar":  "gb_l0",
        "type":     "CDS",
        "level":    0,
        "sequence": "ATGGCA" * 5,
    }
    sc._h_create_part(mock_app, part)
    res = sc._h_create_part(mock_app, part)
    assert isinstance(res, tuple) and res[1] == 409


def test_part_update_unknown_returns_404(mock_app):
    part = {
        "name":     "ghost_part",
        "grammar":  "gb_l0",
        "type":     "CDS",
        "level":    0,
        "sequence": "ATG" * 10,
    }
    res = sc._h_update_part(mock_app, part)
    assert isinstance(res, tuple) and res[1] == 404


# ── Feature library CRUD ──────────────────────────────────────────────────────


def test_feature_library_create_get_delete():
    res = sc._h_create_feature_library(None, {
        "name":         "test_promoter",
        "feature_type": "promoter",
        "strand":       1,
        "color":        "green",
        "sequence":     "AAGCTTCCAGTACAGGCTTGCAGTAGCT",
    })
    assert res["ok"]

    got = sc._h_get_feature_library(None, {"name": "test_promoter"})
    assert got["ok"]
    assert got["feature"]["feature_type"] == "promoter"

    listed = sc._h_list_feature_library(None, {})
    assert any(r["name"] == "test_promoter" for r in listed["features"])

    res = sc._h_delete_feature_library(None, {
        "name": "test_promoter", "feature_type": "promoter",
    })
    assert res["ok"]


def test_feature_library_create_duplicate_returns_409():
    f = {
        "name":         "dup_feat",
        "feature_type": "CDS",
        "strand":       1,
        "sequence":     "ATGGCA" * 6,
    }
    sc._h_create_feature_library(None, f)
    res = sc._h_create_feature_library(None, f)
    assert isinstance(res, tuple) and res[1] == 409


def test_feature_library_update_unknown_returns_404():
    res = sc._h_update_feature_library(None, {
        "name":     "ghost_feat",
        "sequence": "AAATTT",
    })
    assert isinstance(res, tuple) and res[1] == 404


# ── Custom enzymes + enzyme collections CRUD (sweep #24) ─────────────────────

pytestmark_protect = pytest.mark.usefixtures("_protect_user_data")


@pytestmark_protect
def test_custom_enzyme_create_get_update_delete_roundtrip(mock_app):
    payload = {
        "name":     "TestAgentEnz",
        "site":     "GGTACC",
        "fwd_cut":  1, "rev_cut": 5,
        "type":     "II_5overhang",
        "supplier": "agent-test",
    }
    res = sc._h_create_custom_enzyme(mock_app, payload)
    assert res == {"ok": True, "name": "TestAgentEnz"}

    got = sc._h_get_custom_enzyme(None, {"name": "TestAgentEnz"})
    assert got["ok"] and got["enzyme"]["site"] == "GGTACC"

    listed = sc._h_list_custom_enzymes(None, {})
    names = {e["name"] for e in listed["enzymes"]}
    assert "TestAgentEnz" in names

    upd = sc._h_update_custom_enzyme(mock_app, {
        **payload, "supplier": "edited",
    })
    assert upd["ok"]
    assert sc._h_get_custom_enzyme(None, {"name": "TestAgentEnz"}) \
        ["enzyme"]["supplier"] == "edited"

    deleted = sc._h_delete_custom_enzyme(mock_app, {"name": "TestAgentEnz"})
    assert deleted == {"ok": True, "name": "TestAgentEnz"}


@pytestmark_protect
def test_custom_enzyme_create_rejects_builtin_collision(mock_app):
    res = sc._h_create_custom_enzyme(mock_app, {
        "name":     "EcoRI",   # built-in
        "site":     "GAATTC",
        "fwd_cut":  1, "rev_cut": 5,
    })
    assert isinstance(res, tuple) and res[1] == 409


@pytestmark_protect
def test_custom_enzyme_create_rejects_bad_iupac(mock_app):
    res = sc._h_create_custom_enzyme(mock_app, {
        "name":     "BadSite",
        "site":     "GAAZTC",   # Z is not IUPAC
        "fwd_cut":  1, "rev_cut": 5,
    })
    assert isinstance(res, tuple) and res[1] == 400


@pytestmark_protect
def test_custom_enzyme_update_unknown_returns_404(mock_app):
    res = sc._h_update_custom_enzyme(mock_app, {
        "name":     "GhostEnz",
        "site":     "GAATTC",
        "fwd_cut":  1, "rev_cut": 5,
    })
    assert isinstance(res, tuple) and res[1] == 404


@pytestmark_protect
def test_enzyme_collection_create_get_update_delete_roundtrip(mock_app):
    res = sc._h_create_enzyme_collection(mock_app, {
        "name":    "AgentCol",
        "enzymes": ["EcoRI", "BamHI"],
    })
    assert res == {"ok": True, "name": "AgentCol"}

    got = sc._h_get_enzyme_collection(None, {"name": "AgentCol"})
    assert got["ok"]
    assert sorted(got["collection"]["enzymes"]) == ["BamHI", "EcoRI"]

    upd = sc._h_update_enzyme_collection(mock_app, {
        "name":    "AgentCol",
        "enzymes": ["EcoRI", "BamHI", "HindIII"],
    })
    assert upd["ok"]

    renamed = sc._h_update_enzyme_collection(mock_app, {
        "name":     "AgentCol",
        "new_name": "AgentColRenamed",
    })
    assert renamed == {"ok": True, "name": "AgentColRenamed"}

    deleted = sc._h_delete_enzyme_collection(mock_app, {
        "name": "AgentColRenamed",
    })
    assert deleted == {"ok": True, "name": "AgentColRenamed"}


@pytestmark_protect
def test_enzyme_collection_create_duplicate_returns_409(mock_app):
    sc._h_create_enzyme_collection(mock_app, {
        "name":    "DupCol",
        "enzymes": ["EcoRI"],
    })
    res = sc._h_create_enzyme_collection(mock_app, {
        "name":    "DupCol",
        "enzymes": ["EcoRI"],
    })
    assert isinstance(res, tuple) and res[1] == 409


@pytestmark_protect
def test_active_enzyme_collection_get_set_clear(mock_app):
    sc._h_create_enzyme_collection(mock_app, {
        "name":    "ActiveTarget",
        "enzymes": ["EcoRI"],
    })
    assert sc._h_get_active_enzyme_collection(None, {}) \
        == {"ok": True, "name": None}

    set_res = sc._h_set_active_enzyme_collection(mock_app, {
        "name": "ActiveTarget",
    })
    assert set_res == {"ok": True, "name": "ActiveTarget"}
    assert sc._h_get_active_enzyme_collection(None, {}) \
        == {"ok": True, "name": "ActiveTarget"}

    clear_res = sc._h_set_active_enzyme_collection(mock_app, {
        "name": None,
    })
    assert clear_res == {"ok": True, "name": None}


@pytestmark_protect
def test_set_active_enzyme_collection_unknown_returns_404(mock_app):
    res = sc._h_set_active_enzyme_collection(mock_app, {
        "name": "DoesNotExist",
    })
    assert isinstance(res, tuple) and res[1] == 404


@pytestmark_protect
def test_delete_enzyme_collection_clears_active_pointer(mock_app):
    sc._h_create_enzyme_collection(mock_app, {
        "name":    "WillDelete",
        "enzymes": ["EcoRI"],
    })
    sc._h_set_active_enzyme_collection(mock_app, {"name": "WillDelete"})
    assert sc._get_active_enzyme_collection_name() == "WillDelete"
    sc._h_delete_enzyme_collection(mock_app, {"name": "WillDelete"})
    assert sc._get_active_enzyme_collection_name() is None


# ── Parts-bin collection endpoints ────────────────────────────────────────────


def test_list_parts_bins():
    sc._save_parts_bin_collections([
        {"name": "Bin A", "parts": [], "description": "d"},
    ])
    res = sc._h_list_parts_bins(None, {})
    assert res["ok"]
    names = {b["name"] for b in res["parts_bins"]}
    assert "Bin A" in names


def test_set_active_parts_bin_roundtrip():
    sc._save_parts_bin_collections([
        {"name": "Bin A", "parts": [], "description": ""},
        {"name": "Bin B",
         "parts": [{"name": "p1", "sequence": "ATGC"}],
         "description": ""},
    ])
    res = sc._h_set_active_parts_bin(None, {"name": "Bin B"})
    assert res["ok"] and res["active"] == "Bin B"
    assert res["n_parts"] == 1
    assert sc._get_active_parts_bin_name() == "Bin B"


def test_set_active_parts_bin_unknown_404():
    res = sc._h_set_active_parts_bin(None, {"name": "no-such-bin-xyz"})
    assert isinstance(res, tuple) and res[1] == 404


def test_set_active_parts_bin_missing_name_400():
    res = sc._h_set_active_parts_bin(None, {})
    assert isinstance(res, tuple) and res[1] == 400


# ── HMM-database endpoints ─────────────────────────────────────────────────────


def test_list_hmm_databases_has_builtins():
    res = sc._h_list_hmm_databases(None, {})
    assert res["ok"]
    ids = {d["id"] for d in res["hmm_databases"]}
    assert ids   # built-ins always present


def test_get_hmm_database_roundtrip():
    listed = sc._h_list_hmm_databases(None, {})["hmm_databases"]
    assert listed
    some_id = listed[0]["id"]
    res = sc._h_get_hmm_database(None, {"id": some_id})
    assert res["ok"] and res["id"] == some_id


def test_get_hmm_database_unknown_404():
    res = sc._h_get_hmm_database(None, {"id": "no-such-db"})
    assert isinstance(res, tuple) and res[1] == 404


def test_set_active_hmm_database_roundtrip():
    listed = sc._h_list_hmm_databases(None, {})["hmm_databases"]
    assert listed
    some_id = listed[0]["id"]
    res = sc._h_set_active_hmm_database(None, {"id": some_id})
    assert res["ok"] and res["active"] == some_id
    assert (sc._get_setting("hmm_db_active_id", "") or "") == some_id


def test_set_active_hmm_database_unknown_404():
    res = sc._h_set_active_hmm_database(None, {"id": "no-such-db"})
    assert isinstance(res, tuple) and res[1] == 404


def test_delete_hmm_database_unknown_404():
    res = sc._h_delete_hmm_database(None, {"id": "no-such-db"})
    assert isinstance(res, tuple) and res[1] == 404


def test_delete_hmm_database_not_downloaded_ok():
    listed = sc._h_list_hmm_databases(None, {})["hmm_databases"]
    assert listed
    some_id = listed[0]["id"]
    res = sc._h_delete_hmm_database(None, {"id": some_id})
    assert res["ok"] and res["files_removed"] == 0


# ── HMM database add + download (this session) ────────────────────────────────


def test_add_hmm_database_roundtrip():
    res = sc._h_add_hmm_database(
        None, {"name": "My Test HMMs",
               "url": "https://example.org/db.hmm.gz"})
    assert res["ok"] and res["id"] == "my_test_hmms"
    listed = sc._h_list_hmm_databases(None, {})["hmm_databases"]
    by_id = {d["id"]: d for d in listed}
    assert "my_test_hmms" in by_id
    assert by_id["my_test_hmms"]["builtin"] is False
    assert by_id["my_test_hmms"]["ready"] is False  # registered, not fetched


def test_add_hmm_database_missing_name_400():
    res = sc._h_add_hmm_database(
        None, {"url": "https://example.org/db.hmm.gz"})
    assert isinstance(res, tuple) and res[1] == 400


def test_add_hmm_database_bad_url_400():
    res = sc._h_add_hmm_database(
        None, {"name": "Bad URL DB", "url": "ftp://example.org/db.hmm.gz"})
    assert isinstance(res, tuple) and res[1] == 400


def test_add_hmm_database_duplicate_id_409():
    payload = {"name": "Dup DB", "url": "https://example.org/db.hmm.gz"}
    assert sc._h_add_hmm_database(None, payload)["ok"]
    res = sc._h_add_hmm_database(None, dict(payload))
    assert isinstance(res, tuple) and res[1] == 409


def test_add_hmm_database_builtin_collision_409():
    # "Pfam-A" → id "pfam-a", a shipped builtin: must refuse, not shadow.
    res = sc._h_add_hmm_database(
        None, {"name": "Pfam-A", "url": "https://example.org/db.hmm.gz"})
    assert isinstance(res, tuple) and res[1] == 409


def test_download_hmm_database_unknown_404():
    res = sc._h_download_hmm_database(None, {"id": "no-such-db"})
    assert isinstance(res, tuple) and res[1] == 404


def test_download_hmm_database_success(monkeypatch):
    sc._h_add_hmm_database(
        None, {"name": "Stub DB", "url": "https://example.org/db.hmm.gz"})
    calls = {}

    def _fake_perform(entry, **kw):
        calls["entry"] = entry
        return {"id": entry["id"], "n_profiles": 7, "pressed": True,
                "sha256": "abc", "version": "downloaded", "bytes": 1234}

    monkeypatch.setattr(
        "splicecraft_agent._hmm_db_perform_download", _fake_perform)
    res = sc._h_download_hmm_database(None, {"id": "stub_db"})
    assert res["ok"] and res["id"] == "stub_db"
    assert res["n_profiles"] == 7 and res["pressed"] is True
    assert calls["entry"]["id"] == "stub_db"
    # Slot freed after success → a retry isn't spuriously 409.
    assert sc._hmm_db_acquire_download_slot("stub_db")
    sc._hmm_db_release_download_slot("stub_db")


def test_download_hmm_database_inflight_409():
    # A held slot (e.g. a GUI download in flight) blocks the agent path.
    assert sc._hmm_db_acquire_download_slot("pfam-a")
    try:
        res = sc._h_download_hmm_database(None, {"id": "pfam-a"})
        assert isinstance(res, tuple) and res[1] == 409
    finally:
        sc._hmm_db_release_download_slot("pfam-a")


def test_download_hmm_database_failure_500_releases_slot(monkeypatch):
    def _boom(entry, **kw):
        raise RuntimeError("network exploded")

    monkeypatch.setattr("splicecraft_agent._hmm_db_perform_download", _boom)
    res = sc._h_download_hmm_database(None, {"id": "pfam-a"})
    assert isinstance(res, tuple) and res[1] == 500
    # Slot must be released even on failure.
    assert sc._hmm_db_acquire_download_slot("pfam-a")
    sc._hmm_db_release_download_slot("pfam-a")


# ── export-map-image (F6) + list-plasmidsaurus-items alias (F5) ──────────────

def test_export_map_image_endpoint(tmp_path):
    """`export-map-image` renders the loaded record to png / svg, rejects a
    bad format (400) and a missing record (422)."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq("ACGT" * 200), id="pMap", name="pMap")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    app = type("A", (), {"_current_record": rec})()
    fn = sc._state._AGENT_HANDLERS["export-map-image"][0]   # via the registry
    for fmt in ("png", "svg"):
        out = tmp_path / f"map.{fmt}"
        res = fn(app, {"path": str(out)})
        assert res.get("ok") is True, res
        assert res["fmt"] == fmt
        assert out.exists() and out.stat().st_size > 0
    # extension must match the format
    bad = fn(app, {"path": str(tmp_path / "x.png"), "format": "gif"})
    assert isinstance(bad, tuple) and bad[1] == 400
    # no plasmid loaded → 422
    empty = type("A", (), {"_current_record": None})()
    nores = fn(empty, {"path": str(tmp_path / "y.png")})
    assert isinstance(nores, tuple) and nores[1] == 422


def test_list_plasmidsaurus_items_alias_registered():
    """F5: the `list-`-prefixed name is an alias for the same handler."""
    h = sc._state._AGENT_HANDLERS
    assert "list-plasmidsaurus-items" in h
    assert h["list-plasmidsaurus-items"] == h["plasmidsaurus-items"]


# ── check-acceptor-vector + the optimizer's splice block (2026-11-08) ────────


def _acceptor_rec(extra_pair: "tuple[str, str] | None" = None,
                  extra_site: str = "GGTCTC"):
    """An α1-shaped acceptor; with `extra_pair`, a SECOND pair of `extra_site`
    cuts in the backbone releasing a fragment with those overhangs (the
    reported four-site vector)."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    bsa_f, bsa_r = "GGTCTC" + "A", "T" + sc._rc("GGTCTC")
    esp_f, esp_r = "CGTCTC" + "A", "T" + sc._rc("CGTCTC")
    inner = bsa_f + "GGAG" + "TTGGAACCAA" * 20 + "CGCT" + bsa_r
    outer = esp_f + "GGAG" + inner + "GTCA" + esp_r
    backbone = "ATCGATCGAT" * 60
    if extra_pair:
        x_f, x_r = extra_site + "A", "T" + sc._rc(extra_site)
        backbone += (x_f + extra_pair[0] + "CCTTGGAACC" * 30
                     + extra_pair[1] + x_r + "ATCGATCGAT" * 60)
    seq = outer + backbone
    rec = SeqRecord(Seq(seq), id="acc", name="acc",
                    annotations={"topology": "circular",
                                 "molecule_type": "DNA"})
    rec.features.append(SeqFeature(
        FeatureLocation(len(outer) + 100, len(outer) + 500),
        type="rep_origin", qualifiers={"label": ["ori"]}))
    return rec


def test_check_acceptor_vector_passes_a_clean_acceptor():
    app = type("A", (), {"_current_record": _acceptor_rec()})()
    res = sc._h_check_acceptor_vector(app, {})
    assert res["ok"] is True
    assert sorted(res["usable_enzymes"]) == ["BsaI", "Esp3I"]
    assert res["role"] == {"detected": True, "name": "Alpha1",
                           "label": "Alpha1", "confidence": "strict"}


def test_check_acceptor_vector_explains_a_failed_detection():
    """The silence this endpoint exists to end: auto-detect returned None and
    said nothing, for years."""
    app = type("A", (), {"_current_record": _acceptor_rec(("AATG", "GCTT"))})()
    res = sc._h_check_acceptor_vector(app, {})
    assert res["role"]["detected"] is False
    assert "4 BsaI sites" in res["role"]["reason"]
    row = next(r for r in res["enzymes"] if r["enzyme"] == "BsaI")
    assert row["extra"][0]["matches_slot"] == "CDS"
    assert any("COMPETES" in w for w in res["warnings"])


def test_check_acceptor_vector_rejects_unknown_enzyme_loudly():
    app = type("A", (), {"_current_record": _acceptor_rec()})()
    body, code = sc._h_check_acceptor_vector(app, {"enzymes": ["NotAnEnzyme"]})
    assert code == 400 and "unknown enzyme" in body["error"]


def test_check_acceptor_vector_needs_a_record():
    app = type("A", (), {"_current_record": None})()
    body, code = sc._h_check_acceptor_vector(app, {})
    assert code == 422


def test_check_acceptor_vector_rejects_unknown_grammar():
    app = type("A", (), {"_current_record": _acceptor_rec()})()
    body, code = sc._h_check_acceptor_vector(app, {"grammar": "nope"})
    assert code == 404


def test_optimize_protein_reports_cryptic_splice_sites_for_a_plant_host():
    """`hazard_hosts: ["plant"]` buys a 2-motif floor; the response has to say
    what the trained model sees, whether or not the repair was asked for."""
    res = sc._h_optimize_protein(None, {
        "protein": "MKPGSAKVTGAKPGSADETWFCTWFCWFCTL",
        "table": "4097", "hazard_hosts": ["plant"]})
    assert not isinstance(res, tuple), res
    assert res["splice"]["skipped"] is False
    assert res["splice"]["scrubbed"] is False
    assert "n_sites" in res["splice"] and "n_pairs" in res["splice"]
    assert res["splice"]["model"]


def test_optimize_protein_can_scrub_the_sites_it_reports():
    """A protein whose frequency-mode optimisation spells a real donor: the
    opt-in repair must clear it and the response must prove it."""
    body = {"protein": "MKVTGGVNGTAGGVNRTGGVNATGGLKDVNGTSGGVNGTYGGVNHTGG",
            "table": "4097", "hazard_hosts": ["plant"]}
    plain = sc._h_optimize_protein(None, dict(body))
    assert not isinstance(plain, tuple), plain
    scrubbed = sc._h_optimize_protein(None, dict(body, scrub_splice_sites=True))
    assert not isinstance(scrubbed, tuple), scrubbed
    assert scrubbed["splice"]["scrubbed"] is True
    assert scrubbed["splice"]["n_sites"] <= plain["splice"]["n_sites"]
    # whatever it changed, the protein is the invariant
    assert sc._mut_translate(scrubbed["dna"]).rstrip("*") == body["protein"]


def test_optimize_protein_skips_the_scan_when_told_to():
    res = sc._h_optimize_protein(None, {
        "protein": "MKPGSAKVTG", "table": "4097",
        "hazard_hosts": ["plant"], "check_splice": False})
    assert res["splice"]["skipped"] is True
    assert res["splice"]["reason"]


def test_optimize_protein_does_not_scan_for_a_non_plant_host():
    """The model is plant-trained — it is the wrong ruler for an animal cell,
    and a silent wrong answer is worse than no answer."""
    res = sc._h_optimize_protein(None, {
        "protein": "MKPGSAKVTG", "table": "83333",
        "hazard_hosts": ["bacterial"]})
    assert res["splice"]["skipped"] is True


def test_lint_synthesis_acceptor_enzymes_round_trip():
    rec = _acceptor_rec(("AATG", "GCTT"))
    res = sc._h_lint_synthesis(None, {"sequence": str(rec.seq),
                                      "acceptor_enzymes": ["BsaI", "Esp3I"]})
    assert res["ok"] is True
    assert res["stats"]["acceptor_sites"] == {"BsaI": 4, "Esp3I": 2}
    assert any(w["kind"] == "acceptor_site" and w["level"] == "error"
               for w in res["warnings"])


def test_lint_synthesis_rejects_an_unknown_acceptor_enzyme():
    body, code = sc._h_lint_synthesis(None, {"sequence": "ACGT" * 50,
                                             "acceptor_enzymes": ["Nope"]})
    assert code == 400 and "unknown enzyme" in body["error"]


@pytestmark_protect
def test_design_gb_part_says_when_its_entry_vector_is_not_a_clean_acceptor():
    """`_entry_vector_acceptor_overhangs` answers None for BOTH "nothing
    configured" and "the configured vector doesn't cut cleanly", and the
    designer treats both as "use the one-tier layout" — so a design built
    against a four-site vector came back looking normal."""
    # four ESP3I sites: that is the enzyme the domestication design digests
    # the entry vector with, so this is the vector shape that silently
    # degraded the design to one tier.
    rec = _acceptor_rec(("AATG", "GCTT"), extra_site="CGTCTC")
    from io import StringIO
    from Bio import SeqIO
    buf = StringIO()
    SeqIO.write(rec, buf, "genbank")
    gb_text = buf.getvalue()
    sc._set_entry_vector("gb_l0", {"name": "four-site UPD",
                                   "size": len(str(rec.seq)),
                                   "source": "test", "id": "four",
                                   "gb_text": gb_text}, "")
    try:
        res = sc._h_design_gb_part(None, {
            "template": "ATG" + "GCTAAACCTGGTTCTGCT" * 6 + "TAA",
            "start": 0, "end": 3 + 18 * 6 + 3, "part_type": "CDS"})
        assert not isinstance(res, tuple), res
        warns = res["result"].get("warnings") or []
        assert any("not a clean Esp3I acceptor" in w for w in warns), warns
        assert any("check-acceptor-vector" in w for w in warns), warns
    finally:
        sc._set_entry_vector("gb_l0", None, "")


@pytestmark_protect
def test_design_gb_part_is_quiet_when_no_entry_vector_is_configured():
    """No vector is the ordinary one-tier case, not a problem to report."""
    sc._set_entry_vector("gb_l0", None, "")
    res = sc._h_design_gb_part(None, {
        "template": "ATG" + "GCTAAACCTGGTTCTGCT" * 6 + "TAA",
        "start": 0, "end": 3 + 18 * 6 + 3, "part_type": "CDS"})
    assert not isinstance(res, tuple), res
    assert not [w for w in (res["result"].get("warnings") or [])
                if "acceptor" in w]


def test_optimize_protein_refuses_to_repair_without_reporting():
    """A repair whose result is not measured is a claim, not a result."""
    body, code = sc._h_optimize_protein(None, {
        "protein": "MKVTGG", "table": "4097",
        "scrub_splice_sites": True, "check_splice": False})
    assert code == 400 and "drop one of the two" in body["error"]


def test_check_acceptor_vector_dedups_a_repeated_enzyme():
    app = type("A", (), {"_current_record": _acceptor_rec()})()
    res = sc._h_check_acceptor_vector(app, {"enzymes": ["BsaI", "BsaI"]})
    assert [r["enzyme"] for r in res["enzymes"]] == ["BsaI"]
