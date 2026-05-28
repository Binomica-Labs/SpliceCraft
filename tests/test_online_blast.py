"""
test_online_blast — Online BLAST / HMM tab of the BlastModal.

The networked engines (`_ncbi_blast_online`, `_hmmer_web_hmmscan`) are NOT
exercised against the live servers here — that would make the suite slow,
flaky, and dependent on NCBI / EBI uptime. Instead we test:

  * the pure parsers (`_ncbi_blast_parse_xml`, `_hmmer_web_parse_json`) on
    canned payloads, including the XXE/DOCTYPE refusal,
  * query cleaning + the program→database + query-length-limit maps,
  * the modal UI wiring (tab present, program toggle, oversize/empty
    refusal, render → DataTable) WITHOUT dispatching a real search.

The 160×48 boundary regression lives in test_modal_boundaries.py.
"""
from __future__ import annotations

import pytest

import splicecraft as sc

TERMINAL_SIZE = (160, 48)


# ═══════════════════════════════════════════════════════════════════════
# Pure helpers — program/database, query-length limits, query cleaning
# ═══════════════════════════════════════════════════════════════════════

class TestOnlineHelpers:
    def test_db_for_program(self):
        assert sc._ncbi_blast_db_for("blastp") == "nr"
        assert sc._ncbi_blast_db_for("blastx") == "nr"
        assert sc._ncbi_blast_db_for("blastn") == "nt"
        assert sc._ncbi_blast_db_for("tblastn") == "nt"
        assert sc._ncbi_blast_db_for("tblastx") == "nt"

    def test_max_query_len_matches_ncbi_limits(self):
        # NCBI: 1,000,000 for nucleotide queries, 100,000 for protein.
        assert sc._online_max_query_len("blastn") == 1_000_000
        assert sc._online_max_query_len("blastx") == 1_000_000
        assert sc._online_max_query_len("tblastx") == 1_000_000
        assert sc._online_max_query_len("blastp") == 100_000
        assert sc._online_max_query_len("tblastn") == 100_000
        assert sc._online_max_query_len("hmmscan") == 100_000

    def test_clean_query_strips_headers_and_whitespace(self):
        raw = ">seq1 description\nACGT acgt\nAC GT\n"
        assert sc._online_clean_query(raw, "blastn") == "ACGTACGTACGT"

    def test_clean_query_rna_to_dna_for_nucleotide_programs(self):
        assert sc._online_clean_query("acguacgu", "blastn") == "ACGTACGT"
        assert sc._online_clean_query("acguacgu", "blastx") == "ACGTACGT"
        # Protein programs must NOT touch U (it's a valid residue context).
        assert sc._online_clean_query("MKUU", "blastp") == "MKUU"
        assert sc._online_clean_query("MKUU", "hmmscan") == "MKUU"

    def test_clean_query_does_not_truncate(self):
        # The caller enforces the limit + refuses; cleaning never trims.
        big = "A" * 1_500_000
        assert len(sc._online_clean_query(big, "blastn")) == 1_500_000

    def test_program_query_kind(self):
        for p in ("blastn", "blastx", "tblastx"):
            assert sc._program_query_kind(p) == "nt"
        for p in ("blastp", "tblastn", "hmmscan"):
            assert sc._program_query_kind(p) == "protein"


# ═══════════════════════════════════════════════════════════════════════
# NCBI BLAST XML parser
# ═══════════════════════════════════════════════════════════════════════

_BLAST_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_id>gi|12345</Hit_id>
          <Hit_def>pUC19 cloning vector complete sequence</Hit_def>
          <Hit_accession>L09137</Hit_accession>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>120.5</Hsp_bit-score>
              <Hsp_evalue>1.2e-30</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>60</Hsp_query-to>
              <Hsp_hit-from>100</Hsp_hit-from>
              <Hsp_hit-to>159</Hsp_hit-to>
              <Hsp_identity>58</Hsp_identity>
              <Hsp_align-len>60</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_id>gi|67890</Hit_id>
          <Hit_def>some other vector</Hit_def>
          <Hit_accession>MW463917</Hit_accession>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>40.0</Hsp_bit-score>
              <Hsp_evalue>3.0</Hsp_evalue>
              <Hsp_query-from>5</Hsp_query-from>
              <Hsp_query-to>25</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>21</Hsp_hit-to>
              <Hsp_identity>18</Hsp_identity>
              <Hsp_align-len>21</Hsp_align-len>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"""


class TestNcbiBlastParseXml:
    def test_basic(self):
        hits = sc._ncbi_blast_parse_xml(_BLAST_XML, 50)
        assert len(hits) == 2
        h = hits[0]
        assert h["accession"] == "L09137"
        assert h["description"] == "pUC19 cloning vector complete sequence"
        assert h["identity_pct"] == pytest.approx(96.7, abs=0.05)
        assert h["aln_len"] == 60
        assert h["evalue"] == pytest.approx(1.2e-30)
        assert h["bit_score"] == pytest.approx(120.5)
        assert (h["q_start"], h["q_end"]) == (1, 60)
        assert (h["s_start"], h["s_end"]) == (100, 159)

    def test_max_hits_cap(self):
        hits = sc._ncbi_blast_parse_xml(_BLAST_XML, 1)
        assert len(hits) == 1
        assert hits[0]["accession"] == "L09137"

    def test_empty_output(self):
        xml = ("<?xml version='1.0'?><BlastOutput><BlastOutput_iterations>"
               "<Iteration><Iteration_hits></Iteration_hits></Iteration>"
               "</BlastOutput_iterations></BlastOutput>")
        assert sc._ncbi_blast_parse_xml(xml, 50) == []

    def test_external_doctype_allowed(self):
        # Real NCBI BLAST XML opens with an external <!DOCTYPE … .dtd>.
        # The parser must accept it (regression: _safe_xml_parse used to
        # refuse every DOCTYPE, rejecting all real NCBI responses).
        with_dtd = _BLAST_XML.replace(
            '<?xml version="1.0"?>\n',
            '<?xml version="1.0"?>\n<!DOCTYPE BlastOutput PUBLIC '
            '"-//NCBI//NCBI BlastOutput/EN" '
            '"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">\n', 1)
        hits = sc._ncbi_blast_parse_xml(with_dtd, 50)
        assert len(hits) == 2
        assert hits[0]["accession"] == "L09137"

    def test_rejects_internal_entity_subset(self):
        # A DOCTYPE carrying an internal ENTITY subset (billion-laughs) is
        # still refused → parser raises RuntimeError.
        evil = ('<?xml version="1.0"?><!DOCTYPE x [<!ENTITY a "boom">]>'
                "<BlastOutput></BlastOutput>")
        with pytest.raises(RuntimeError):
            sc._ncbi_blast_parse_xml(evil, 50)


class TestSafeXmlDtd:
    """The `allow_dtd` opt-in on _safe_xml_parse keeps XXE/billion-laughs
    defenses while permitting NCBI's external-DTD reference."""

    def test_external_dtd_only_with_flag(self):
        import xml.etree.ElementTree as ET
        xml = ('<?xml version="1.0"?><!DOCTYPE BlastOutput PUBLIC "x" '
               '"http://example/x.dtd"><BlastOutput><n>1</n></BlastOutput>')
        with pytest.raises(ET.ParseError):
            sc._safe_xml_parse(xml)                 # default: refused
        root = sc._safe_xml_parse(xml, allow_dtd=True)
        assert root.tag == "BlastOutput"

    def test_internal_subset_refused_even_with_flag(self):
        import xml.etree.ElementTree as ET
        xml = ('<?xml version="1.0"?><!DOCTYPE x [<!ENTITY a "boom">]>'
               "<x>hi</x>")
        with pytest.raises(ET.ParseError):
            sc._safe_xml_parse(xml, allow_dtd=True)

    def test_standalone_entity_refused(self):
        import xml.etree.ElementTree as ET
        with pytest.raises(ET.ParseError):
            sc._safe_xml_parse('<!ENTITY a "x"><x/>', allow_dtd=True)


class TestFeatureProtein:
    def test_prefers_translation_qualifier(self):
        entry = {"feature_type": "CDS", "sequence": "ATGAAATAA",
                 "qualifiers": {"translation": ["MKLV"]}}
        assert sc._feature_protein(entry) == "MKLV"

    def test_strips_trailing_stop_in_qualifier(self):
        entry = {"sequence": "", "qualifiers": {"translation": ["MK*"]}}
        assert sc._feature_protein(entry) == "MK"

    def test_translates_sequence_when_no_qualifier(self):
        entry = {"feature_type": "CDS", "sequence": "ATGAAATTTTAA",
                 "qualifiers": {}}
        assert sc._feature_protein(entry) == "MKF"   # ATG AAA TTT (TAA=*)

    def test_honours_codon_start(self):
        entry = {"feature_type": "CDS", "sequence": "CATGAAATTTTAA",
                 "qualifiers": {"codon_start": ["2"]}}
        assert sc._feature_protein(entry) == "MKF"   # skip 1 base → ATG…

    def test_empty(self):
        assert sc._feature_protein({"sequence": "", "qualifiers": {}}) == ""


# ═══════════════════════════════════════════════════════════════════════
# EBI HMMER JSON parser
# ═══════════════════════════════════════════════════════════════════════

class TestHmmerWebParseJson:
    def test_result_hits_shape(self):
        obj = {"status": "SUCCESS", "result": {"hits": [
            {"acc": "PF00069.27", "name": "Pkinase",
             "desc": "Protein kinase domain", "evalue": "3.1e-40",
             "score": "135.2", "ndom": 2},
            {"acc": "PF07714.20", "name": "PK_Tyr_Ser-Thr",
             "desc": "Protein tyrosine kinase", "evalue": "1e-10",
             "score": "40.0", "nincluded": 1},
        ]}}
        rows = sc._hmmer_web_parse_json(obj, 50)
        assert len(rows) == 2
        assert rows[0]["acc"] == "PF00069.27"
        assert rows[0]["name"] == "Pkinase"
        assert rows[0]["description"] == "Protein kinase domain"
        assert rows[0]["evalue"] == pytest.approx(3.1e-40)
        assert rows[0]["bit_score"] == pytest.approx(135.2)
        assert rows[0]["n_dom"] == 2
        # Domain count falls back to nincluded when ndom is absent.
        assert rows[1]["n_dom"] == 1

    def test_metadata_drives_name_and_description(self):
        # Live EBI v1 shape: top-level `name` is an internal numeric id and
        # `desc` is null; the human-readable family name + description live
        # in `metadata`. Regression for the blank-description / numeric-name
        # bug.
        obj = {"result": {"hits": [{
            "acc": "PF07714.23", "name": "000017756", "desc": None,
            "evalue": 2.7e-95, "score": 318.9, "ndom": 1,
            "is_included": True,
            "metadata": {
                "accession": "PF07714", "identifier": "PK_Tyr_Ser-Thr",
                "description": "Protein tyrosine and serine/threonine kinase",
                "clan": "CL0016", "type": "Domain",
                "external_link":
                    "https://www.ebi.ac.uk/interpro/entry/pfam/PF07714"},
        }]}}
        row = sc._hmmer_web_parse_json(obj, 50)[0]
        assert row["acc"] == "PF07714.23"
        assert row["name"] == "PK_Tyr_Ser-Thr"        # not "000017756"
        assert row["description"].startswith("Protein tyrosine")
        assert row["clan"] == "CL0016"
        assert row["type"] == "Domain"
        assert "interpro" in row["link"]
        assert row["included"] is True

    def test_included_flag_false(self):
        obj = {"result": {"hits": [{
            "acc": "PF1", "is_included": False,
            "metadata": {"identifier": "x", "description": "d"}}]}}
        assert sc._hmmer_web_parse_json(obj, 50)[0]["included"] is False

    def test_alternate_container_shapes(self):
        hit = {"acc": "PF00001", "name": "x", "desc": "d",
               "evalue": "1e-5", "score": "10"}
        assert len(sc._hmmer_web_parse_json(
            {"results": {"hits": [hit]}}, 50)) == 1
        assert len(sc._hmmer_web_parse_json({"hits": [hit]}, 50)) == 1

    def test_domain_count_from_domains_list(self):
        obj = {"result": {"hits": [
            {"acc": "PF1", "domains": [{}, {}, {}]},
        ]}}
        assert sc._hmmer_web_parse_json(obj, 50)[0]["n_dom"] == 3

    def test_max_hits_cap(self):
        obj = {"result": {"hits": [{"acc": f"PF{i}"} for i in range(10)]}}
        assert len(sc._hmmer_web_parse_json(obj, 3)) == 3

    def test_garbage_inputs(self):
        assert sc._hmmer_web_parse_json(None, 50) == []
        assert sc._hmmer_web_parse_json("not a dict", 50) == []
        assert sc._hmmer_web_parse_json({"result": {}}, 50) == []
        assert sc._hmmer_web_parse_json({"result": {"hits": "nope"}}, 50) == []


# ═══════════════════════════════════════════════════════════════════════
# Modal UI — no live search is ever dispatched in these tests
# ═══════════════════════════════════════════════════════════════════════

class TestOnlineTabUI:
    async def test_online_widgets_present(self, tiny_record,
                                          isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            for wid in ("#tab-blast-online", "#online-query",
                        "#online-program", "#online-db", "#online-maxhits",
                        "#btn-online-search", "#online-table",
                        "#online-status", "#btn-online-from-plasmid",
                        "#btn-online-from-feature"):
                modal.query_one(wid)

    async def test_program_hmmscan_toggles_ui(self, tiny_record,
                                              isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#online-program", sc.Select).value = "hmmscan"
            await pilot.pause()
            assert modal.query_one("#online-db-col").display is False
            assert str(modal.query_one(
                "#btn-online-search", sc.Button).label) == "Search Pfam"
            # Back to a BLAST program restores the DB control + label.
            modal.query_one("#online-program", sc.Select).value = "blastp"
            await pilot.pause()
            assert modal.query_one("#online-db-col").display is True
            assert modal.query_one("#online-db", sc.Input).value == "nr"
            assert str(modal.query_one(
                "#btn-online-search", sc.Button).label) == "Search NCBI"

    async def test_from_feature_protein_mode_tracks_program(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            # hmmscan → translate CDS (protein mode).
            modal.query_one("#online-program", sc.Select).value = "hmmscan"
            await pilot.pause()
            modal.query_one("#btn-online-from-feature", sc.Button).press()
            await pilot.pause()
            assert modal._online_feature_protein_mode is True
            await pilot.press("escape")          # close the plasmid picker
            await pilot.pause()
            # blastn → raw nucleotides (not protein mode).
            modal.query_one("#online-program", sc.Select).value = "blastn"
            await pilot.pause()
            modal.query_one("#btn-online-from-feature", sc.Button).press()
            await pilot.pause()
            assert modal._online_feature_protein_mode is False
            await pilot.press("escape")
            await pilot.pause()

    async def test_program_switch_resets_query_on_type_flip(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            ta = modal.query_one("#online-query", sc.TextArea)
            prog = modal.query_one("#online-program", sc.Select)
            # blastn (nt) → blastp (protein): query alphabet flips → cleared
            ta.text = "ACGTACGTACGT"
            prog.value = "blastp"
            await pilot.pause()
            assert ta.text == ""
            # blastp → tblastn (both protein): preserved
            ta.text = "MKTAYIAKQR"
            prog.value = "tblastn"
            await pilot.pause()
            assert ta.text == "MKTAYIAKQR"
            # tblastn (protein) → blastn (nt): flips → cleared
            prog.value = "blastn"
            await pilot.pause()
            assert ta.text == ""
            # blastn → blastx (both nt): preserved
            ta.text = "ACGTACGT"
            prog.value = "blastx"
            await pilot.pause()
            assert ta.text == "ACGTACGT"

    async def test_blast_detail_is_enriched(self, tiny_record,
                                            isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal._online_set_searching(True)
            modal._online_blast_done(
                "blastn", sc._ncbi_blast_parse_xml(_BLAST_XML, 50), None)
            await pilot.pause()
            table = modal.query_one("#online-table", sc.DataTable)
            table.move_cursor(row=0)
            await pilot.pause()
            detail = str(modal.query_one("#online-detail", sc.Static).render())
            assert "L09137" in detail
            assert "ncbi.nlm.nih.gov" in detail        # accession link
            assert "pUC19" in detail                    # description
            assert "bits" in detail                     # stats

    async def test_hmm_detail_is_enriched(self, tiny_record,
                                          isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#online-program", sc.Select).value = "hmmscan"
            await pilot.pause()
            hits = sc._hmmer_web_parse_json({"result": {"hits": [{
                "acc": "PF07714.23", "evalue": 1e-95, "score": 318.9,
                "ndom": 1, "is_included": True,
                "metadata": {"identifier": "PK_Tyr_Ser-Thr",
                             "description": "Protein tyrosine kinase",
                             "clan": "CL0016", "type": "Domain",
                             "external_link":
                             "https://www.ebi.ac.uk/interpro/entry/pfam/"
                             "PF07714"}}]}}, 50)
            modal._online_set_searching(True)
            modal._online_hmm_done(hits, None)
            await pilot.pause()
            table = modal.query_one("#online-table", sc.DataTable)
            table.move_cursor(row=0)
            await pilot.pause()
            detail = str(modal.query_one("#online-detail", sc.Static).render())
            assert "PK_Tyr_Ser-Thr" in detail
            assert "Protein tyrosine kinase" in detail
            assert "CL0016" in detail
            assert "interpro" in detail

    async def test_empty_query_refused(self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#btn-online-search", sc.Button).press()
            await pilot.pause()
            assert modal._online_busy is False
            txt = str(modal.query_one("#online-status", sc.Static).render())
            assert "first" in txt.lower() or "query" in txt.lower()

    async def test_oversize_query_refused_without_dispatch(
            self, tiny_record, isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#online-program", sc.Select).value = "blastp"
            await pilot.pause()
            # 100,001 > the 100,000 protein limit → refuse, never dispatch.
            modal.query_one("#online-query", sc.TextArea).text = "A" * 100_001
            await pilot.pause()
            modal.query_one("#btn-online-search", sc.Button).press()
            await pilot.pause()
            assert modal._online_busy is False
            txt = str(modal.query_one("#online-status", sc.Static).render())
            assert "100,000" in txt and "most" in txt.lower()

    async def test_blast_done_renders_table(self, tiny_record,
                                            isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal._online_set_searching(True)  # simulate in-flight
            hits = sc._ncbi_blast_parse_xml(_BLAST_XML, 50)
            modal._online_blast_done("blastn", hits, None)
            await pilot.pause()
            table = modal.query_one("#online-table", sc.DataTable)
            assert table.row_count == 2
            assert len(table.columns) == len(modal._ONLINE_BLAST_COLS)
            assert modal._online_busy is False
            txt = str(modal.query_one("#online-status", sc.Static).render())
            assert "2 hits" in txt

    async def test_hmm_done_renders_table(self, tiny_record,
                                          isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#online-program", sc.Select).value = "hmmscan"
            await pilot.pause()
            modal._online_set_searching(True)
            rows = sc._hmmer_web_parse_json(
                {"result": {"hits": [
                    {"acc": "PF00069", "name": "Pkinase", "desc": "kinase",
                     "evalue": "1e-9", "score": "55", "ndom": 1}]}}, 50)
            modal._online_hmm_done(rows, None)
            await pilot.pause()
            table = modal.query_one("#online-table", sc.DataTable)
            assert table.row_count == 1
            assert len(table.columns) == len(modal._ONLINE_HMM_COLS)
            assert modal._online_busy is False

    async def test_search_dispatch_blast_renders(
            self, monkeypatch, tiny_record, isolated_library):
        # Stub the network engine so the full press→worker→render path is
        # exercised without touching NCBI.
        captured = {}

        def _fake(query, program, database, max_hits,
                  progress_cb=None, cancel_event=None):
            captured.update(program=program, database=database,
                            max_hits=max_hits, qlen=len(query))
            if progress_cb:
                progress_cb("Submitting test…")
            return sc._ncbi_blast_parse_xml(_BLAST_XML, max_hits)

        monkeypatch.setattr(sc, "_ncbi_blast_online", _fake)
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#online-program", sc.Select).value = "blastp"
            await pilot.pause()
            modal.query_one("#online-query", sc.TextArea).text = "MKTAYIAKQR"
            await pilot.pause()
            modal.query_one("#btn-online-search", sc.Button).press()
            for _ in range(50):
                await pilot.pause()
                if not modal._online_busy:
                    break
            assert captured["program"] == "blastp"
            assert captured["database"] == "nr"   # auto DB for blastp
            table = modal.query_one("#online-table", sc.DataTable)
            assert table.row_count == 2
            assert modal._online_busy is False

    async def test_search_dispatch_hmmscan_renders(
            self, monkeypatch, tiny_record, isolated_library):
        def _fake(protein, max_hits, progress_cb=None, cancel_event=None):
            if progress_cb:
                progress_cb("Submitting test…")
            return sc._hmmer_web_parse_json(
                {"result": {"hits": [
                    {"acc": "PF00069", "name": "Pkinase", "desc": "kinase",
                     "evalue": "1e-9", "score": "55", "ndom": 1}]}}, max_hits)

        monkeypatch.setattr(sc, "_hmmer_web_hmmscan", _fake)
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            modal.query_one("#online-program", sc.Select).value = "hmmscan"
            await pilot.pause()
            modal.query_one("#online-query", sc.TextArea).text = "MKTAYIAKQR"
            await pilot.pause()
            modal.query_one("#btn-online-search", sc.Button).press()
            for _ in range(50):
                await pilot.pause()
                if not modal._online_busy:
                    break
            table = modal.query_one("#online-table", sc.DataTable)
            assert table.row_count == 1
            assert len(table.columns) == len(modal._ONLINE_HMM_COLS)
            assert modal._online_busy is False

    async def test_footer_disabled_on_online_tab(self, tiny_record,
                                                 isolated_library):
        app = sc.PlasmidApp()
        app._preload_record = tiny_record
        async with app.run_test(size=TERMINAL_SIZE) as pilot:
            await pilot.pause()
            app.push_screen(sc.BlastModal())
            await pilot.pause()
            modal = app.screen
            tabs = modal.query_one("#blast-tabs", sc.TabbedContent)
            run_btn = modal.query_one("#btn-blast-run", sc.Button)
            # Switch by focusing a widget in the Online tab — the same path
            # a real user's tab click takes (TabPane.Focused → active). A
            # bare `tabs.active = …` races the mount-time focus of the
            # Local query box, which reasserts the Local tab.
            modal.query_one("#online-query", sc.TextArea).focus()
            for _ in range(10):
                await pilot.pause()
                if run_btn.disabled:
                    break
            assert tabs.active == "tab-blast-online"
            assert run_btn.disabled
            assert modal.query_one("#btn-blast-build", sc.Button).disabled
            # Switching back to the Local tab re-enables them.
            modal.query_one("#blast-query", sc.TextArea).focus()
            for _ in range(10):
                await pilot.pause()
                if not run_btn.disabled:
                    break
            assert not run_btn.disabled
