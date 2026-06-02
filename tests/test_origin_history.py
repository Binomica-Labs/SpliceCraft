"""
test_origin_history — every newly-generated plasmid lands with at least
a minimal construction-history root, so the History tab is never empty
for a plasmid the user made (the reported pei311 symptom: a Golden Braid
clone with no history).

Covers:
  * the origin-history helpers (`_gb_text_is_circular`,
    `_origin_history_descriptor`, `_build_origin_history_xml`,
    `_ensure_entry_origin_history`);
  * the Constructor enzyme-gate fix — `_build_history_for_assembly` no
    longer drops ALL history when the grammar's enzyme can't resolve;
  * the integration backfill: a record auto-persisted through the normal
    save path lands with `history_xml`.
"""
from __future__ import annotations

import pytest

import splicecraft as sc


# ── topology read ───────────────────────────────────────────────────────


class TestGbTextIsCircular:
    def test_linear_locus(self):
        assert sc._gb_text_is_circular(
            "LOCUS  FRAG  825 bp  DNA  linear  SYN") is False

    def test_circular_locus(self):
        assert sc._gb_text_is_circular(
            "LOCUS  Pei311  2559 bp  DNA  circular  SYN") is True

    def test_empty_or_none_defaults_circular(self):
        assert sc._gb_text_is_circular("") is True
        assert sc._gb_text_is_circular(None) is True

    def test_unmarked_defaults_circular(self):
        assert sc._gb_text_is_circular("LOCUS x 100 bp") is True

    def test_only_locus_line_consulted(self):
        # "linear" elsewhere in the body must not flip a circular LOCUS.
        gb = "LOCUS x 100 bp DNA circular\nCOMMENT linear plasmid map\n"
        assert sc._gb_text_is_circular(gb) is True


# ── source → descriptor mapping ─────────────────────────────────────────


class TestOriginHistoryDescriptor:
    @pytest.mark.parametrize("source,manip", [
        ("synthesis:FRAG-Pei311", "synthesizeDNA"),
        ("paste:foo", "pasteSequence"),
        ("simulator:pcr", "amplifyPCR"),
        ("file:thing.gb", "importFile"),
        ("auto-detect (import)", "importFile"),
        ("sequencing", "sequencingRead"),
        ("plasmidsaurus:run:s", "sequencingRead"),
        ("constructor:gb_l0:backbone", "assembly"),
        ("id:Pei311", "createDocument"),
        ("", "createDocument"),
        ("something-unexpected", "createDocument"),
    ])
    def test_mapping(self, source, manip):
        op, m = sc._origin_history_descriptor(source)
        assert m == manip
        assert op in ("createDocument", "importFile", "insertFragment")


# ── XML builder ─────────────────────────────────────────────────────────


class TestBuildOriginHistoryXml:
    def test_parses_round_trip(self):
        xml = sc._build_origin_history_xml(
            name="Pei311", seq_len=2559, circular=True, source="id:Pei311")
        assert xml
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None
        assert root.seq_len == 2559
        assert root.circular is True
        assert root.name.startswith("Pei311")
        assert root.input_summaries[0]["manipulation"] == "createDocument"

    def test_linear_synthesis_fragment(self):
        xml = sc._build_origin_history_xml(
            name="FRAG-Pei311", seq_len=825, circular=False,
            source="synthesis:FRAG-Pei311")
        root = sc._parse_commercialsaas_history(xml)
        assert root.circular is False
        assert root.input_summaries[0]["manipulation"] == "synthesizeDNA"

    def test_origin_root_is_a_leaf(self):
        # An origin root documents creation, not an assembly — no parents.
        xml = sc._build_origin_history_xml(
            name="x", seq_len=10, circular=True, source="paste:x")
        root = sc._parse_commercialsaas_history(xml)
        assert root.parents == []


# ── in-place ensure helper ──────────────────────────────────────────────


class TestEnsureEntryOriginHistory:
    def test_stamps_when_absent(self):
        entry = {"name": "Pei311", "id": "Pei311", "size": 2559,
                 "source": "id:Pei311",
                 "gb_text": "LOCUS Pei311 2559 bp DNA circular"}
        sc._ensure_entry_origin_history(entry)
        assert "history_xml" in entry
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root is not None and root.circular is True

    def test_noop_when_history_present(self):
        existing = "<HistoryTree><Node ID='0' name='x.dna'/></HistoryTree>"
        entry = {"name": "X", "size": 10, "source": "id:X",
                 "history_xml": existing}
        sc._ensure_entry_origin_history(entry)
        assert entry["history_xml"] == existing   # untouched

    def test_linear_topology_from_gb_text(self):
        entry = {"name": "amp", "size": 500, "source": "simulator:pcr",
                 "gb_text": "LOCUS amp 500 bp DNA linear"}
        sc._ensure_entry_origin_history(entry)
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root.circular is False

    def test_non_dict_is_safe(self):
        sc._ensure_entry_origin_history(None)     # must not raise
        sc._ensure_entry_origin_history("nope")


# ── Constructor enzyme-gate fix (the pei311 symptom) ────────────────────


class TestConstructorEnzymeGate:
    """`ConstructorModal._build_history_for_assembly` uses no instance
    state, so a bare object() stands in for `self`."""

    def _call(self, grammar):
        return sc.ConstructorModal._build_history_for_assembly(
            object(),
            name="MyAssembly", product_seq_len=3000,
            gid="gb_l0", grammar=grammar,
            entry_vector={"name": "pUPD"},
            parts=[{"name": "FRAG-Pei311", "size": 825}],
            source_level=0)

    def test_history_built_when_enzyme_missing(self):
        xml = self._call({})              # grammar with no "enzyme" → was None
        assert xml
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None
        # Lineage still recorded (vector + part as parents).
        assert len(root.parents) >= 1
        # No regenerated-site marker without an enzyme.
        assert root.regenerated_sites == []

    def test_history_built_when_enzyme_present(self):
        xml = self._call({"enzyme": "BsaI"})
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None
        assert any(s["name"] == "BsaI" for s in root.regenerated_sites)


# ── integration: normal save path backfills history ─────────────────────


class TestOriginHistoryIntegration:
    async def test_auto_persisted_record_gets_history(
            self, tiny_record, isolated_library):
        """A record persisted through the ordinary canvas-save path (the
        pei311 case — source `id:`) lands in the library WITH a
        construction-history root."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = sc._load_library()
            assert lib, "preloaded record should auto-persist"
            assert any(e.get("history_xml") for e in lib), (
                "a generated plasmid must land with construction history"
            )
            # And that history parses cleanly.
            hist = next(e["history_xml"] for e in lib if e.get("history_xml"))
            assert sc._parse_commercialsaas_history(hist) is not None
            app.exit()

    async def test_inplace_edit_appends_step_via_add_entry(
            self, tiny_record, isolated_library):
        """Driving the real chokepoint: saving v1 then re-saving an
        edited (longer) v2 under the SAME id appends an editSequence step
        chaining the pre-edit version, while origin-only saves don't."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            lib = app.query_one(sc.LibraryPanel)
            lib.add_entry(_rec("ATGCATGCATGCATGCATGC"))          # v1, 20 bp
            await pilot.pause(0.1)
            lib.add_entry(_rec("ATGCATGCATGCATGCATGCGGGGGG"))    # v2, 26 bp
            await pilot.pause(0.1)
            entry = sc._find_library_entry_by_id("P")
            assert entry and entry.get("history_xml")
            root = sc._parse_commercialsaas_history(entry["history_xml"])
            assert root.input_summaries[0]["manipulation"] == "editSequence"
            assert root.seq_len == 26
            assert root.parents, "edit must chain the pre-edit version"
            app.exit()


# ── in-place edit appends an editSequence step (user choice 2026-06-01) ──


def _rec(seq, *, circular=True):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    r = SeqRecord(Seq(seq), id="P", name="P")
    r.annotations["molecule_type"] = "DNA"
    r.annotations["topology"] = "circular" if circular else "linear"
    return r


class TestEditHistoryAppend:
    def _entry_for(self, record):
        gb = sc._record_to_gb_text(record)
        return {"name": "P", "id": "P", "size": len(record.seq),
                "gb_text": gb}

    def test_length_change_appends_edit_node(self):
        prev = _rec("ATGCATGCAT")            # 10 bp
        rec = _rec("ATGCATGCATGGGG")         # 14 bp — edited (insert)
        entry = self._entry_for(rec)
        sc._maybe_append_edit_history(
            entry, rec, sc._record_to_gb_text(prev),
            prev_size=10, prev_history="")
        assert "history_xml" in entry
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root.operation == "replace"
        assert root.input_summaries[0]["manipulation"] == "editSequence"
        assert root.seq_len == 14
        assert len(root.parents) >= 1        # pre-edit version is the parent

    def test_no_change_is_noop(self):
        rec = _rec("ATGCATGCAT")
        entry = self._entry_for(rec)
        sc._maybe_append_edit_history(
            entry, rec, sc._record_to_gb_text(rec),
            prev_size=10, prev_history="")
        assert "history_xml" not in entry    # identical sequence → nothing

    def test_same_length_substitution_appends(self):
        prev = _rec("AAAAAAAAAA")
        rec = _rec("AAAAAAAAAC")             # same length, one substitution
        entry = self._entry_for(rec)
        sc._maybe_append_edit_history(
            entry, rec, sc._record_to_gb_text(prev),
            prev_size=10, prev_history="")
        assert "history_xml" in entry
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root.input_summaries[0]["manipulation"] == "editSequence"

    def test_prior_history_becomes_parent_lineage(self):
        # The pre-edit version carried its own history → it must survive
        # as the parent subtree (lineage chains through edits).
        prev_hist = sc._build_origin_history_xml(
            name="P", seq_len=10, circular=True, source="synthesis:P")
        prev = _rec("ATGCATGCAT")
        rec = _rec("ATGCATGCATTT")
        entry = self._entry_for(rec)
        sc._maybe_append_edit_history(
            entry, rec, sc._record_to_gb_text(prev),
            prev_size=10, prev_history=prev_hist)
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root.parents, "edit node must have a parent"
        parent = root.parents[0]
        # parent is the synthesized pre-edit version
        assert parent.input_summaries[0]["manipulation"] == "synthesizeDNA"

    def test_non_dict_is_safe(self):
        sc._maybe_append_edit_history(None, _rec("AT"), "", 2, "")

    def test_empty_new_seq_is_noop(self):
        entry = {"name": "P", "size": 0, "gb_text": ""}
        sc._maybe_append_edit_history(entry, _rec(""), "x", 5, "")
        assert "history_xml" not in entry


# ── bloat guard: generated history can't grow without bound ──────────────


class TestEditHistoryBloatCap:
    def test_chain_stays_bounded_over_many_edits(self):
        """400 sequence-changing edits must NOT yield a 400-deep history
        — the recent-window cap keeps it bounded (stored in BOTH the
        library and collections files, so unbounded growth = real bloat)."""
        prev_hist = ""
        prev_len = 10
        for i in range(400):
            new_len = 11 + i
            rec = _rec("A" * new_len)
            entry = {"name": "P", "id": "P", "size": new_len,
                     "gb_text": f"LOCUS P {new_len} bp DNA circular"}
            sc._maybe_append_edit_history(
                entry, rec, "", prev_size=prev_len, prev_history=prev_hist)
            prev_hist = entry.get("history_xml", "")
            prev_len = new_len
        assert prev_hist
        root = sc._parse_commercialsaas_history(prev_hist)
        n = sc._history_node_count(root)
        assert n <= sc._HISTORY_TREE_MAX_NODES + 8, f"unbounded: {n} nodes"
        assert len(prev_hist) <= sc._HISTORY_XML_MAX_BYTES

    def test_truncate_marks_truncation(self):
        root = sc._CommercialSaaSHistoryNode.new(
            name="top.dna", seq_len=100, circular=True,
            operation="replace", node_id=0)
        cur = root
        for i in range(20):
            child = sc._CommercialSaaSHistoryNode.new(
                name=f"v{i}.dna", seq_len=99 - i, circular=True,
                operation="insertFragment", node_id=0)
            cur.add_parent(child)
            cur = child
        sc._truncate_history_to_recent(root, max_nodes=5)
        assert sc._history_node_count(root) <= 8
        assert any("truncated" in nd.name.lower() for nd in root.walk())

    def test_no_truncation_under_budget(self):
        root = sc._CommercialSaaSHistoryNode.new(
            name="x.dna", seq_len=10, circular=True,
            operation="replace", node_id=0)
        child = sc._CommercialSaaSHistoryNode.new(
            name="y.dna", seq_len=9, circular=True,
            operation="insertFragment", node_id=0)
        root.add_parent(child)
        sc._truncate_history_to_recent(root, max_nodes=256)
        assert sc._history_node_count(root) == 2          # untouched
        assert not any("truncated" in nd.name.lower()
                       for nd in root.walk())

    def test_huge_prev_history_collapses_to_leaf(self):
        huge = "x" * (sc._HISTORY_XML_MAX_BYTES + 1000)
        rec = _rec("ATGCATGCATGCAT")                       # 14 bp
        entry = {"name": "P", "id": "P", "size": 14,
                 "gb_text": "LOCUS P 14 bp DNA circular"}
        sc._maybe_append_edit_history(
            entry, rec, "", prev_size=10, prev_history=huge)
        assert "history_xml" in entry
        assert len(entry["history_xml"]) <= sc._HISTORY_XML_MAX_BYTES
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert sc._history_node_count(root) <= 3          # edited + leaf


# ── one-shot legacy backfill — DATA SAFETY IS CRITICAL ──────────────────


class TestOriginHistoryBackfill:
    """The migration that stamps origin history on legacy history-less
    plasmids. SACRED: it must NEVER drop, reorder, or rewrite entries —
    only ADD `history_xml`. These tests lock that down."""

    def test_pure_backfill_is_additive_and_count_invariant(self):
        entries = [
            {"id": "a", "name": "a", "size": 1, "source": "id:a",
             "gb_text": "", "status": "DESIGNING"},
            "not-a-dict",                                  # tolerated, skipped
            {"id": "b", "name": "b", "size": 2, "source": "id:b",
             "gb_text": "", "history_xml": "<HistoryTree/>"},
            {"id": "c", "name": "c", "size": 3, "source": "id:c",
             "gb_text": ""},
        ]
        n_before = len(entries)
        out, n_changed = sc._backfill_origin_history(entries)
        assert out is entries                              # in place
        assert len(out) == n_before                        # COUNT invariant
        assert n_changed == 2                              # a + c
        assert out[0].get("history_xml") and out[3].get("history_xml")
        assert out[0]["status"] == "DESIGNING"             # other fields kept
        assert out[2]["history_xml"] == "<HistoryTree/>"   # existing untouched
        # Idempotent.
        _, again = sc._backfill_origin_history(entries)
        assert again == 0

    def test_library_load_backfills_legacy_and_preserves_data(
            self, tmp_path, monkeypatch):
        import json
        tmp_lib = tmp_path / "library.json"
        monkeypatch.setattr(sc, "_LIBRARY_FILE", tmp_lib)
        monkeypatch.setattr(sc, "_library_cache", None)
        monkeypatch.setattr(sc, "_origin_history_backfill_done", False)
        legacy = [
            {"id": "Pei311", "name": "Pei311", "size": 2559, "n_feats": 3,
             "source": "id:Pei311", "status": "VERIFIED",
             "gb_text": "LOCUS Pei311 2559 bp DNA circular"},
            {"id": "FRAG", "name": "FRAG", "size": 825, "n_feats": 1,
             "source": "synthesis:frag",
             "gb_text": "LOCUS FRAG 825 bp DNA linear"},
        ]
        tmp_lib.write_text(json.dumps(
            {"_schema_version": 1, "entries": legacy}))
        loaded = sc._load_library()
        assert len(loaded) == 2                            # SACRED: count kept
        assert all(e.get("history_xml") for e in loaded)
        pei = next(e for e in loaded if e["id"] == "Pei311")
        assert pei["status"] == "VERIFIED" and pei["size"] == 2559
        assert sc._parse_commercialsaas_history(pei["history_xml"]) is not None
        # Persisted to disk, count preserved there too.
        on_disk = json.loads(tmp_lib.read_text())["entries"]
        assert len(on_disk) == 2
        assert all(e.get("history_xml") for e in on_disk)

    def test_library_backfill_does_not_overwrite_existing_history(
            self, tmp_path, monkeypatch):
        import json
        tmp_lib = tmp_path / "library.json"
        monkeypatch.setattr(sc, "_LIBRARY_FILE", tmp_lib)
        monkeypatch.setattr(sc, "_library_cache", None)
        monkeypatch.setattr(sc, "_origin_history_backfill_done", False)
        keep = sc._build_origin_history_xml(
            name="Keep", seq_len=10, circular=True, source="paste:keep")
        tmp_lib.write_text(json.dumps({"_schema_version": 1, "entries": [
            {"id": "Keep", "name": "Keep", "size": 10, "source": "paste:keep",
             "gb_text": "LOCUS Keep 10 bp DNA circular", "history_xml": keep},
        ]}))
        loaded = sc._load_library()
        assert loaded[0]["history_xml"] == keep            # untouched

    def test_collections_load_backfills_embedded_plasmids(
            self, tmp_path, monkeypatch):
        import json
        tmp_coll = tmp_path / "collections.json"
        monkeypatch.setattr(sc, "_COLLECTIONS_FILE", tmp_coll)
        monkeypatch.setattr(sc, "_collections_cache", None)
        monkeypatch.setattr(
            sc, "_collections_origin_history_backfill_done", False)
        colls = [{"name": "Phase Genomics", "plasmids": [
            {"id": "Pei311", "name": "Pei311", "size": 2559,
             "source": "id:Pei311",
             "gb_text": "LOCUS Pei311 2559 bp DNA circular"},
        ]}]
        tmp_coll.write_text(json.dumps(
            {"_schema_version": 1, "entries": colls}))
        out = sc._load_collections()
        plasmids = out[0]["plasmids"]
        assert len(plasmids) == 1                          # count preserved
        assert plasmids[0].get("history_xml")


# ── Domesticator L0 clone records insert + entry-vector lineage ─────────


class TestDomesticatedCloneHistory:
    """A Domesticator clone is a single-part L0 assembly — the part's
    amplicon ligated into its entry vector. `_build_history_for_l0_clone`
    must record that lineage so a synthesis→parts→clone plasmid's History
    reads like a Constructor build, NOT a bare `createDocument` leaf (the
    reported pei311 symptom, again: 2026-06-02). Uses no instance state."""

    def _call(self, grammar):
        return sc._build_history_for_l0_clone(
            name="Pei311", seq_len=2553,
            grammar_id="gb_l0", grammar=grammar,
            entry_vector={"name": "FFE_1_ENTRY_UPD"},
            part={"name": "Pei311 insert", "size": 819,
                  "sequence": "ATGC" * 205},
        )

    def test_clone_history_records_insert_and_vector(self):
        xml = self._call({"enzyme": "Esp3I"})
        assert xml
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None
        # An assembly node, NOT the bare createDocument leaf.
        assert root.operation == "insertFragment"
        # Both inputs recorded as parent fragments. (_parent_node_for_entry
        # stamps a ".dna" suffix on node names, same as the root + the
        # Constructor's parents — strip it for the identity check.)
        names = {p.name.removesuffix(".dna") for p in root.parents}
        assert "FFE_1_ENTRY_UPD" in names, "entry vector must be a parent"
        assert "Pei311 insert" in names, "insert part must be a parent"
        # The L0 Type IIS enzyme lands as a regenerated-site marker.
        assert any(s["name"] == "Esp3I" for s in root.regenerated_sites)

    def test_clone_history_built_when_enzyme_missing(self):
        # Enzyme is a DETAIL — lineage still recorded, just no site marker
        # (same policy as _build_history_for_assembly after the pei311 fix).
        xml = self._call({})
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None
        assert root.operation == "insertFragment"
        assert len(root.parents) == 2
        assert root.regenerated_sites == []


# ── over-cap prior history collapses WITH a truncation marker ───────────


class TestEditHistoryOverCapMarker:
    """When an in-place edit's prior history exceeds the per-edit re-nest
    budget it's collapsed to a leaf — but the leaf must carry an
    "(earlier history truncated)" marker so the History viewer signals
    dropped lineage instead of a bare, parent-less leaf (2026-06-02).
    Mirrors the marker the node-budget cap already leaves."""

    @classmethod
    def _has_marker(cls, node) -> bool:
        if node.name == "(earlier history truncated)":
            return True
        return any(cls._has_marker(c) for c in node.parents)

    def test_over_cap_prior_history_leaves_marker(self):
        entry = {"name": "Big", "id": "Big",
                 "gb_text": "LOCUS Big 30 bp DNA circular"}
        rec = _rec("ATGCATGCATGCATGCATGCATGCATGCAT")        # 30 bp (≠ prev)
        oversized = "<HistoryTree>" + "x" * (sc._HISTORY_XML_MAX_BYTES + 10)
        sc._maybe_append_edit_history(
            entry, rec, prev_gb_text="LOCUS Big 20 bp DNA circular",
            prev_size=20, prev_history=oversized,
        )
        assert "history_xml" in entry
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root is not None
        assert root.parents, "edit must chain the collapsed pre-edit version"
        assert self._has_marker(root), (
            "over-cap prior history must leave a truncation marker"
        )

    def test_within_cap_prior_history_has_no_spurious_marker(self):
        # A small prior history re-nests normally — no truncation marker.
        entry = {"name": "Sm", "id": "Sm",
                 "gb_text": "LOCUS Sm 30 bp DNA circular"}
        rec = _rec("ATGCATGCATGCATGCATGCATGCATGCAT")        # 30 bp
        small_prev = sc._build_origin_history_xml(
            name="Sm", seq_len=20, circular=True, source="paste:x")
        sc._maybe_append_edit_history(
            entry, rec, prev_gb_text="LOCUS Sm 20 bp DNA circular",
            prev_size=20, prev_history=small_prev,
        )
        assert "history_xml" in entry
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert not self._has_marker(root)
