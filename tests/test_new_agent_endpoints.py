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
