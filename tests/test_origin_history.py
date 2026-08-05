"""
test_origin_history — every newly-generated plasmid lands with at least
a minimal construction-history root, so the History tab is never empty
for a plasmid the user made (the reported demo311 symptom: a Golden Braid
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
            "LOCUS  Demo311  2559 bp  DNA  circular  SYN") is True

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
        ("synthesis:FRAG-Demo311", "synthesizeDNA"),
        ("paste:foo", "pasteSequence"),
        ("simulator:pcr", "amplifyPCR"),
        ("file:thing.gb", "importFile"),
        ("auto-detect (import)", "importFile"),
        ("sequencing", "sequencingRead"),
        ("plasmidsaurus:run:s", "sequencingRead"),
        ("constructor:gb_l0:backbone", "assembly"),
        # P0-2 (agent-API feedback): agent assembly endpoints save with an
        # `agent:<method>` source — they must map to an assembly node, not
        # the flat `createDocument` they fell through to pre-fix.
        ("agent:golden-gate", "assembly"),
        ("agent:traditional", "assembly"),
        ("agent:gibson", "assembly"),
        ("id:Demo311", "createDocument"),
        ("", "createDocument"),
        ("something-unexpected", "createDocument"),
    ])
    def test_mapping(self, source, manip):
        op, m = sc._origin_history_descriptor(source)
        assert m == manip
        assert op in ("createDocument", "importFile", "insertFragment")

    @pytest.mark.parametrize("source", [
        "agent:golden-gate", "agent:traditional", "agent:gibson",
    ])
    def test_agent_assembly_is_insert_fragment(self, source):
        """The operation (not just the manipulation) must be the assembly
        `insertFragment` so an agent-built construct reads as an assembly
        in the History tab."""
        op, _m = sc._origin_history_descriptor(source)
        assert op == "insertFragment"


# ── XML builder ─────────────────────────────────────────────────────────


class TestBuildOriginHistoryXml:
    def test_parses_round_trip(self):
        xml = sc._build_origin_history_xml(
            name="Demo311", seq_len=2559, circular=True, source="id:Demo311")
        assert xml
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None
        assert root.seq_len == 2559
        assert root.circular is True
        assert root.name.startswith("Demo311")
        assert root.input_summaries[0]["manipulation"] == "createDocument"

    def test_linear_synthesis_fragment(self):
        xml = sc._build_origin_history_xml(
            name="FRAG-Demo311", seq_len=825, circular=False,
            source="synthesis:FRAG-Demo311")
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
        entry = {"name": "Demo311", "id": "Demo311", "size": 2559,
                 "source": "id:Demo311",
                 "gb_text": "LOCUS Demo311 2559 bp DNA circular"}
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


# ── Constructor enzyme-gate fix (the demo311 symptom) ────────────────────


class TestConstructorEnzymeGate:
    """`ConstructorModal._build_history_for_assembly` uses no instance
    state, so a bare object() stands in for `self`."""

    def _call(self, grammar):
        return sc.ConstructorModal._build_history_for_assembly(
            object(),
            name="MyAssembly", product_seq_len=3000,
            gid="gb_l0", grammar=grammar,
            entry_vector={"name": "pUPD"},
            parts=[{"name": "FRAG-Demo311", "size": 825}],
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
        demo311 case — source `id:`) lands in the library WITH a
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

    async def test_canvas_ctrl_s_preserves_and_appends_history(
            self, tiny_record, isolated_library):
        """The async canvas Ctrl+S worker (`action_save` -> `_save_worker`)
        must follow the SAME construction-history contract as add_entry:
        preserve a built plasmid's prior history across a re-save AND append
        an editSequence step when the canvas sequence changed. Regression
        for 'no history on clonings and mods' — this path used to write the
        library entry with NO history field at all, silently dropping a
        clone / mod / .dna lineage on the first Ctrl+S."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        from Bio.Seq import Seq
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            rec = app._current_record
            rid = rec.id
            # Stamp a rich construction history on the stored entry, the way
            # a clone / mod save would, so the re-save must PRESERVE it.
            prior = sc._build_origin_history_xml(
                name=str(rec.name), seq_len=len(rec.seq), circular=True,
                source="constructor:gb_l0")
            assert prior
            entries = sc._load_library()
            assert any(e.get("id") == rid for e in entries)
            for e in entries:
                if e.get("id") == rid:
                    e["history_xml"] = prior
            sc._save_library(entries)
            prev_len = len(rec.seq)
            # Edit the loaded sequence and fire the async Ctrl+S worker.
            rec.seq = rec.seq + Seq("GGGGGG")
            app._unsaved = True
            app._source_path = None        # library save, not a file
            app.action_save()
            for _ in range(12):
                await pilot.pause(0.05)
            entry = sc._find_library_entry_by_id(rid)
            assert entry and entry.get("history_xml"), \
                "canvas Ctrl+S dropped the construction history"
            root = sc._parse_commercialsaas_history(entry["history_xml"])
            assert root is not None
            # An actual edit -> editSequence root chaining the pre-edit
            # lineage (whose subtree preserves the prior history).
            assert root.input_summaries[0]["manipulation"] == "editSequence"
            assert root.seq_len == prev_len + 6
            assert root.parents, "edit must chain the pre-edit lineage"
            app.exit()

    async def test_canvas_ctrl_s_no_edit_keeps_prior_history(
            self, tiny_record, isolated_library):
        """A canvas Ctrl+S with NO sequence change (e.g. a feature-only
        tweak, or a reflexive save right after loading a built plasmid)
        must keep the prior construction history verbatim — not drop it,
        not spuriously append an edit node."""
        from tests.test_smoke import _build_app, TERMINAL_SIZE
        app = _build_app(tiny_record, isolated_library)
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            await pilot.pause(0.05)
            rec = app._current_record
            rid = rec.id
            prior = sc._build_origin_history_xml(
                name=str(rec.name), seq_len=len(rec.seq), circular=True,
                source="constructor:gb_l0")
            entries = sc._load_library()
            for e in entries:
                if e.get("id") == rid:
                    e["history_xml"] = prior
            sc._save_library(entries)
            app._unsaved = True
            app._source_path = None
            app.action_save()              # no seq edit
            for _ in range(12):
                await pilot.pause(0.05)
            entry = sc._find_library_entry_by_id(rid)
            assert entry and entry.get("history_xml") == prior, \
                "re-save without an edit must keep prior history verbatim"
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
        monkeypatch.setattr(sc._state, "_LIBRARY_FILE", tmp_lib)
        monkeypatch.setattr(sc._state, "_library_cache", None)
        monkeypatch.setattr(sc._state, "_origin_history_backfill_done", False)
        legacy = [
            {"id": "Demo311", "name": "Demo311", "size": 2559, "n_feats": 3,
             "source": "id:Demo311", "status": "VERIFIED",
             "gb_text": "LOCUS Demo311 2559 bp DNA circular"},
            {"id": "FRAG", "name": "FRAG", "size": 825, "n_feats": 1,
             "source": "synthesis:frag",
             "gb_text": "LOCUS FRAG 825 bp DNA linear"},
        ]
        tmp_lib.write_text(json.dumps(
            {"_schema_version": 1, "entries": legacy}))
        loaded = sc._load_library()
        assert len(loaded) == 2                            # SACRED: count kept
        assert all(e.get("history_xml") for e in loaded)
        demo = next(e for e in loaded if e["id"] == "Demo311")
        assert demo["status"] == "VERIFIED" and demo["size"] == 2559
        assert sc._parse_commercialsaas_history(demo["history_xml"]) is not None
        # Persisted to disk, count preserved there too.
        on_disk = json.loads(tmp_lib.read_text())["entries"]
        assert len(on_disk) == 2
        # history_xml is blob-dehydrated on disk now (history_ref + blob);
        # rehydrate to confirm the backfilled history survives the round-trip.
        assert all(sc._rehydrate_entry(e).get("history_xml") for e in on_disk)

    def test_library_backfill_does_not_overwrite_existing_history(
            self, tmp_path, monkeypatch):
        import json
        tmp_lib = tmp_path / "library.json"
        monkeypatch.setattr(sc._state, "_LIBRARY_FILE", tmp_lib)
        monkeypatch.setattr(sc._state, "_library_cache", None)
        monkeypatch.setattr(sc._state, "_origin_history_backfill_done", False)
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
        monkeypatch.setattr(sc._state, "_COLLECTIONS_FILE", tmp_coll)
        monkeypatch.setattr(sc._state, "_collections_cache", None)
        monkeypatch.setattr(sc._state, "_collections_origin_history_backfill_done", False)
        colls = [{"name": "Acme Labs", "plasmids": [
            {"id": "Demo311", "name": "Demo311", "size": 2559,
             "source": "id:Demo311",
             "gb_text": "LOCUS Demo311 2559 bp DNA circular"},
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
    reported demo311 symptom, again: 2026-06-02). Uses no instance state."""

    def _call(self, grammar):
        return sc._build_history_for_l0_clone(
            name="Demo311", seq_len=2553,
            grammar_id="gb_l0", grammar=grammar,
            entry_vector={"name": "FFE_1_ENTRY_UPD"},
            part={"name": "Demo311 insert", "size": 819,
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
        assert "Demo311 insert" in names, "insert part must be a parent"
        # The L0 Type IIS enzyme lands as a regenerated-site marker.
        assert any(s["name"] == "Esp3I" for s in root.regenerated_sites)

    def test_clone_history_built_when_enzyme_missing(self):
        # Enzyme is a DETAIL — lineage still recorded, just no site marker
        # (same policy as _build_history_for_assembly after the demo311 fix).
        xml = self._call({})
        root = sc._parse_commercialsaas_history(xml)
        assert root is not None
        assert root.operation == "insertFragment"
        assert len(root.parents) == 2
        assert root.regenerated_sites == []


# ── over-cap prior history collapses WITH a truncation marker ───────────


class TestEditHistoryOverCapMarker:
    """An in-place edit re-nests the pre-edit version's whole lineage.

    When the chained tree overshoots the byte ceiling, the reduction is
    SHEDDING COSMETIC PAYLOAD (`<Features>` / `<HistoryColors>` /
    `<StrandColors>` — annotation + colour snapshots that record no
    manipulation), NOT dropping ancestors. Pre-2026-08-04 a prior history
    over 256 KB collapsed to a bare leaf, discarding every imported step;
    the largest real imported history in the reference library is 233 KB,
    i.e. one edit away from that cliff.

    The "(earlier history truncated)" marker survives for the genuinely
    unfittable case, so a viewer never renders a silently parent-less leaf.
    """

    @classmethod
    def _has_marker(cls, node) -> bool:
        if node.name == "(earlier history truncated)":
            return True
        return any(cls._has_marker(c) for c in node.parents)

    @staticmethod
    def _names(node) -> set:
        return {n.name for n in node.walk()}

    @staticmethod
    def _bulky_prior_history(n_ancestors: int, feats_per_node: int) -> str:
        """A VALID prior history: a chain of `n_ancestors` real steps,
        each weighed down with a cosmetic `<Features>` snapshot."""
        import xml.etree.ElementTree as ET
        root = sc._CommercialSaaSHistoryNode.new(
            name="ANC-0.dna", seq_len=1000, circular=True,
            operation="goldenGateAssembly")
        cur = root
        for i in range(1, n_ancestors):
            cur.add_input_summary(manipulation="insert", name1="BsaI",
                                   name2="BsaI", val1=i, val2=i + 10)
            block = ET.SubElement(cur.element, "Features")
            for f in range(feats_per_node):
                el = ET.SubElement(block, "Feature")
                el.set("name", f"feat-{i}-{f}")
                el.set("padding", "z" * 200)      # cosmetic bulk
            nxt = sc._CommercialSaaSHistoryNode.new(
                name=f"ANC-{i}.dna", seq_len=1000 - i, circular=True,
                operation="insertFragment")
            cur.add_parent(nxt)
            cur = nxt
        return sc._serialize_commercialsaas_history(root)

    def test_bulky_prior_history_keeps_every_step(self, monkeypatch):
        """The whole ancestor chain survives an edit whose re-nested tree
        blows the ceiling — the cosmetic snapshots are what gets shed."""
        # Squeeze the ceiling rather than generating megabytes: the
        # function reads the module global at call time.
        monkeypatch.setattr(sc, "_HISTORY_XML_MAX_BYTES", 8 * 1024)
        prev = self._bulky_prior_history(n_ancestors=6, feats_per_node=12)
        assert len(prev) > 8 * 1024, "fixture must actually exceed the ceiling"
        entry = {"name": "Big", "id": "Big",
                 "gb_text": "LOCUS Big 30 bp DNA circular"}
        rec = _rec("ATGCATGCATGCATGCATGCATGCATGCAT")        # 30 bp (≠ prev)
        sc._maybe_append_edit_history(
            entry, rec, prev_gb_text="LOCUS Big 20 bp DNA circular",
            prev_size=20, prev_history=prev,
        )
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root is not None
        names = self._names(root)
        for i in range(6):
            assert f"ANC-{i}.dna" in names, (
                f"ancestor ANC-{i} was dropped — a manipulation was lost "
                f"to a byte budget"
            )
        assert not self._has_marker(root), (
            "nothing was truncated; the marker must not appear"
        )
        # …and the record of WHAT each step did survives too.
        assert any(n.input_summaries for n in root.walk())
        # The cosmetic payload is what went.
        assert "<Features" not in entry["history_xml"]

    def test_strip_cosmetic_keeps_the_manipulation_record(self):
        """`_strip_history_cosmetic_payload` touches ONLY the cosmetic
        tags — every element that records a manipulation stays."""
        import xml.etree.ElementTree as ET
        node = sc._CommercialSaaSHistoryNode.new(
            name="P.dna", seq_len=500, circular=False,
            operation="amplifyFragment")
        node.add_input_summary(manipulation="amplify", val1=10, val2=400)
        node.add_regenerated_site("BsaI", 42, 1)
        node.add_oligo(name="F", sequence="ACGT")
        node.add_parameter(name="polymerase", value="blunt")
        ET.SubElement(node.element, "Primers")
        for tag in ("Features", "HistoryColors", "StrandColors"):
            ET.SubElement(node.element, tag)
        removed = sc._strip_history_cosmetic_payload(node)
        assert removed == 3
        xml = sc._serialize_commercialsaas_history(node)
        for gone in ("<Features", "<HistoryColors", "<StrandColors"):
            assert gone not in xml
        for kept in ("InputSummary", "RegeneratedSite", "Oligo",
                     "Parameter", "Primers", "amplifyFragment"):
            assert kept in xml, f"{kept} must survive a cosmetic shed"

    def test_unfittable_history_still_leaves_marker(self, monkeypatch):
        """A tree that can't fit even after shedding falls back to the
        bare edit node — and still flags that lineage existed."""
        monkeypatch.setattr(sc, "_HISTORY_XML_MAX_BYTES", 512)
        # Names alone blow the budget, so shedding cosmetics can't help.
        prev = sc._build_origin_history_xml(
            name="L" * 400, seq_len=20, circular=True, source="paste:x")
        entry = {"name": "Big", "id": "Big",
                 "gb_text": "LOCUS Big 30 bp DNA circular"}
        rec = _rec("ATGCATGCATGCATGCATGCATGCATGCAT")
        sc._maybe_append_edit_history(
            entry, rec, prev_gb_text="LOCUS Big 20 bp DNA circular",
            prev_size=20, prev_history=prev,
        )
        root = sc._parse_commercialsaas_history(entry["history_xml"])
        assert root is not None and root.parents
        assert self._has_marker(root), (
            "an unfittable prior history must still leave a truncation marker"
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


# ── provenance recovery from saved `.dna` originals ─────────────────────


def _gb(name: str, seq: str, circular: bool = True) -> str:
    """Minimal but real GenBank text — `_gb_text_sequence` reads ORIGIN."""
    topo = "circular" if circular else "linear"
    body = "".join(
        f"{i * 60 + 1:>9} {seq[i * 60:(i + 1) * 60]}\n"
        for i in range((len(seq) + 59) // 60)
    )
    return (f"LOCUS       {name} {len(seq)} bp DNA {topo} SYN\n"
            f"FEATURES             Location/Qualifiers\n"
            f"ORIGIN\n{body}//\n")


class TestRecoverHistoryFromDna:
    """`recover-history-from-dna` restores lineage from the `.dna` files
    already on disk, matched by EXACT sequence identity.

    Motivating measurement (2026-08-04, reference library of 1,184
    plasmids): 41 entries were sequence-identical to a saved `.dna` whose
    history carried up to 18 MORE nodes than what the library held — real
    build records, recoverable, invisible."""

    SEQ_A = ("ATGCATGCAT" * 12)          # 120 bp
    SEQ_B = ("GGCCTTAAGG" * 12)          # 120 bp

    def _sidecar(self, entry_id: str, seq: str, history_xml: str) -> None:
        """Write a REAL .dna sidecar through the shipping writer, so the
        test exercises the same bytes the app would read back."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        rec = SeqRecord(Seq(seq), id=entry_id, name=entry_id)
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "circular"
        data = sc._write_commercialsaas_dna_bytes(rec, history_xml=history_xml)
        p = sc._dna_sidecar_path(entry_id)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)

    def _rich_history(self, name: str, n_ancestors: int) -> str:
        root = sc._CommercialSaaSHistoryNode.new(
            name=f"{name}.dna", seq_len=120, circular=True,
            operation="goldenGateAssembly")
        cur = root
        for i in range(1, n_ancestors):
            nxt = sc._CommercialSaaSHistoryNode.new(
                name=f"{name}-anc{i}.dna", seq_len=100 + i, circular=True,
                operation="insertFragment")
            cur.add_parent(nxt)
            cur = nxt
        return sc._serialize_commercialsaas_history(root)

    def _seed(self, *, entry_name="pRich", entry_id="pRich",
               stored_nodes=1, seq=None):
        seq = seq or self.SEQ_A
        stub = sc._build_origin_history_xml(
            name=entry_name, seq_len=len(seq), circular=True,
            source="id:x")
        sc._save_collections([
            {"name": "Work",
             "plasmids": [{"id": entry_id, "name": entry_name,
                           "gb_text": _gb(entry_name, seq),
                           "size": len(seq), "history_xml": stub}],
             "saved": "2026-01-01"},
        ])
        sc._set_active_collection_name("Work")
        sc._restore_library_from_active_collection()
        return stub

    def _stored(self, name="pRich"):
        c = next(c for c in sc._load_collections() if c["name"] == "Work")
        e = next(p for p in c["plasmids"] if p["name"] == name)
        return sc._rehydrate_entry(e)

    def test_dry_run_reports_without_writing(self):
        stub = self._seed()
        self._sidecar("pRich", self.SEQ_A, self._rich_history("pRich", 6))
        r = sc._h_recover_history_from_dna(None, {})     # dry_run defaults True
        assert r["ok"] and r["dry_run"] is True
        assert r["updated_count"] == 1
        hit = r["updated"][0]
        assert hit["name"] == "pRich"
        assert hit["nodes_before"] == 1 and hit["nodes_after"] == 6
        # Nothing on disk changed.
        assert self._stored()["history_xml"] == stub

    def test_apply_restores_the_richer_lineage(self):
        self._seed()
        self._sidecar("pRich", self.SEQ_A, self._rich_history("pRich", 6))
        r = sc._h_recover_history_from_dna(None, {"dry_run": False})
        assert r["ok"] and r["updated_count"] == 1
        root = sc._parse_commercialsaas_history(self._stored()["history_xml"])
        names = {n.name for n in root.walk()}
        for i in range(1, 6):
            assert f"pRich-anc{i}.dna" in names

    def test_only_history_is_touched(self):
        """SACRED: sequence, features, name, id and size must be untouched
        — this is an annotation-only repair."""
        self._seed()
        before = self._stored()
        self._sidecar("pRich", self.SEQ_A, self._rich_history("pRich", 6))
        sc._h_recover_history_from_dna(None, {"dry_run": False})
        after = self._stored()
        for field in ("id", "name", "size", "gb_text"):
            assert after[field] == before[field], f"{field} was modified"
        assert after["history_xml"] != before["history_xml"]

    def test_never_thins_an_existing_history(self):
        """A source with FEWER nodes than what's stored is ignored — the
        guard that makes running this repeatedly safe."""
        self._seed()
        rich = self._rich_history("pRich", 9)
        sc._save_collections([
            {"name": "Work",
             "plasmids": [{"id": "pRich", "name": "pRich",
                           "gb_text": _gb("pRich", self.SEQ_A),
                           "size": len(self.SEQ_A), "history_xml": rich}],
             "saved": "2026-01-01"}])
        sc._restore_library_from_active_collection()
        self._sidecar("pRich", self.SEQ_A, self._rich_history("pRich", 3))
        r = sc._h_recover_history_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 0
        assert sc._history_node_count_of_xml(
            self._stored()["history_xml"]) == 9

    def test_matches_across_a_rename(self):
        """The `.dna` is filed under a different name — an exact sequence
        match still recovers it. This is the common real case: a plasmid
        renamed after it was built."""
        self._seed(entry_name="pNewName", entry_id="pNewName")
        self._sidecar("SomeOldName", self.SEQ_A,
                       self._rich_history("SomeOldName", 5))
        r = sc._h_recover_history_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 1
        assert r["updated"][0]["name"] == "pNewName"
        root = sc._parse_commercialsaas_history(
            self._stored("pNewName")["history_xml"])
        assert "SomeOldName.dna" in root.name

    def test_different_sequence_is_never_matched(self):
        """A `.dna` whose sequence differs — even by one base — must not
        donate its lineage. Wrong provenance is worse than none."""
        self._seed()
        off_by_one = self.SEQ_A[:-1] + ("G" if self.SEQ_A[-1] != "G" else "C")
        self._sidecar("Impostor", off_by_one,
                       self._rich_history("Impostor", 8))
        r = sc._h_recover_history_from_dna(None, {"dry_run": False})
        assert r["updated_count"] == 0

    def test_collection_and_name_filters(self):
        self._seed()
        self._sidecar("pRich", self.SEQ_A, self._rich_history("pRich", 6))
        assert sc._h_recover_history_from_dna(
            None, {"collection": "Nope"})["updated_count"] == 0
        assert sc._h_recover_history_from_dna(
            None, {"name": "ghost"})["updated_count"] == 0
        assert sc._h_recover_history_from_dna(
            None, {"collection": "Work", "name": "pRich"})["updated_count"] == 1

    def test_invalid_dry_run_rejected(self):
        assert sc._h_recover_history_from_dna(None, {"dry_run": "yes"})[1] == 400

    def test_scan_reports_when_it_stops_early(self, monkeypatch):
        """A truncated scan MUST say so. Without the flag, a run that hit
        the memory budget and recovered nothing looks exactly like a
        complete run that found nothing to recover — and the user would
        stop looking."""
        self._seed()
        self._sidecar("pRich", self.SEQ_A, self._rich_history("pRich", 6))
        self._sidecar("pOther", self.SEQ_B, self._rich_history("pOther", 4))
        # Squeeze the byte budget so the index fills after the first file.
        monkeypatch.setattr(sc, "_HISTORY_RECOVER_MAX_INDEX_BYTES", 1)
        index, note = sc._scan_dna_originals_for_history()
        assert len(index) == 1, "budget must stop the scan after one entry"
        assert note and "memory budget" in note
        r = sc._h_recover_history_from_dna(None, {})
        assert r["truncated"] is True
        assert "memory budget" in r["note"]

    def test_complete_scan_is_not_flagged_truncated(self):
        self._seed()
        self._sidecar("pRich", self.SEQ_A, self._rich_history("pRich", 6))
        r = sc._h_recover_history_from_dna(None, {})
        assert r["truncated"] is False
        assert r["updated_count"] == 1
        assert "note" not in r

    def test_file_count_cap_is_also_reported(self, monkeypatch):
        self._seed()
        self._sidecar("pRich", self.SEQ_A, self._rich_history("pRich", 6))
        self._sidecar("pOther", self.SEQ_B, self._rich_history("pOther", 4))
        monkeypatch.setattr(sc, "_HISTORY_RECOVER_MAX_SIDECARS", 1)
        _index, note = sc._scan_dna_originals_for_history()
        assert note and ".dna files were scanned" in note


class TestGeneratedVsImportedHistorySerialisation:
    """`_finalize_generated_history` renumbers node IDs so a merged
    lineage can't hold duplicates. It must run on trees WE build and
    NEVER on an imported one — rewriting a source file's own identifiers
    would break the byte-faithful `.dna` round-trip.

    This regression exists because the split was briefly violated: a
    bulk edit routed `_stamp_history_root_date` — which touches ONLY the
    root `date` attribute of imported XML — through the builder path.
    """

    IMPORTED = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<HistoryTree>'
        '<Node name="prod.dna" type="DNA" seqLen="900" strandedness="double"'
        ' ID="7" circular="1" operation="goldenGateAssembly">'
        '<Node name="vec.dna" type="DNA" seqLen="500" strandedness="double"'
        ' ID="3" circular="1" operation="invalid"/>'
        '<Node name="ins.dna" type="DNA" seqLen="400" strandedness="double"'
        ' ID="9" circular="0" operation="invalid"/>'
        '</Node></HistoryTree>'
    )

    @staticmethod
    def _ids(xml):
        import xml.etree.ElementTree as ET
        return [n.get("ID") for n in ET.fromstring(xml).iter("Node")]

    def test_stamping_a_date_preserves_imported_node_ids(self):
        """The date stamp is the ONLY permitted change."""
        assert self._ids(self.IMPORTED) == ["7", "3", "9"]
        out = sc._stamp_history_root_date(self.IMPORTED, "2026-06-11")
        assert self._ids(out) == ["7", "3", "9"], (
            "imported node IDs were rewritten — the .dna round-trip is no "
            "longer byte-faithful"
        )
        root = sc._parse_commercialsaas_history(out)
        assert root.date == "2026-06-11"

    def test_stamp_is_a_no_op_without_a_date(self):
        assert sc._stamp_history_root_date(self.IMPORTED, None) is self.IMPORTED
        assert sc._stamp_history_root_date(self.IMPORTED, "") is self.IMPORTED

    def test_builder_path_renumbers_merged_lineages(self):
        """The other half of the contract: a tree we BUILD gets unique
        IDs, so nesting two imported subtrees (each numbered from its own
        origin) can't leave two nodes claiming the same identity."""
        root = sc._CommercialSaaSHistoryNode.new(
            name="product.dna", seq_len=900, circular=True,
            operation="insertFragment", node_id=7)
        for sub in (self.IMPORTED, self.IMPORTED):
            root.add_parent(sc._parse_commercialsaas_history(sub))
        merged = sc._serialize_commercialsaas_history(root)
        ids = self._ids(merged)
        assert len(ids) != len(set(ids)), "fixture must actually collide"
        out = sc._finalize_generated_history(root)
        ids = self._ids(out)
        assert ids == [str(i) for i in range(len(ids))]
        assert len(ids) == len(set(ids))

    def test_renumber_preserves_every_node_and_attribute(self):
        """Renumbering must touch the ID and nothing else."""
        import xml.etree.ElementTree as ET
        root = sc._parse_commercialsaas_history(self.IMPORTED)
        before = [(n.get("name"), n.get("seqLen"), n.get("operation"))
                  for n in ET.fromstring(self.IMPORTED).iter("Node")]
        sc._history_renumber_node_ids(root)
        after = [(n.get("name"), n.get("seqLen"), n.get("operation"))
                 for n in ET.fromstring(
                     sc._serialize_commercialsaas_history(root)).iter("Node")]
        assert after == before
