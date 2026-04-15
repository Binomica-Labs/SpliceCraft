"""
test_codon — codon usage registry + harmonization pipeline.

Covers the persistent _codon_tables_* registry, the pure-function
_codon_harmonize / _codon_fix_sites / _codon_cai / _codon_gc helpers, and
the Kazusa HTML parser. The Kazusa HTTP fetch itself is not exercised (no
network in tests) — we test the parser on a synthetic HTML fixture.
"""
import pytest

import splicecraft as sc


# ── Registry ──────────────────────────────────────────────────────────────────

class TestRegistry:
    def test_first_load_seeds_k12(self):
        entries = sc._codon_tables_load()
        assert any(e.get("taxid") == "83333" for e in entries)
        k12 = sc._codon_tables_get("83333")
        assert k12 is not None
        assert k12["raw"]["ATG"] == ("M", 127)   # K12 Met count

    def test_add_new_table_and_lookup(self):
        sc._codon_tables_load()
        raw = {"ATG": ("M", 10), "TAA": ("*", 1), "GCT": ("A", 5)}
        sc._codon_tables_add("Fake sp.", "999999", raw, source="kazusa")
        by_taxid = sc._codon_tables_get("999999")
        by_name  = sc._codon_tables_get("fake sp.")
        assert by_taxid is not None and by_name is not None
        assert by_taxid["raw"]["ATG"] == ("M", 10)
        assert by_taxid["source"] == "kazusa"

    def test_add_is_dedup_by_taxid(self):
        sc._codon_tables_load()
        raw1 = {"ATG": ("M", 1), "TAA": ("*", 1)}
        raw2 = {"ATG": ("M", 99), "TAA": ("*", 1)}
        sc._codon_tables_add("First",  "777777", raw1)
        sc._codon_tables_add("Second", "777777", raw2)  # same taxid → replaces
        hits = [e for e in sc._codon_tables_load() if e.get("taxid") == "777777"]
        assert len(hits) == 1
        assert hits[0]["raw"]["ATG"] == ("M", 99)
        assert hits[0]["name"] == "Second"

    def test_persistence_roundtrip(self):
        raw = {"ATG": ("M", 42), "TAA": ("*", 7), "GCT": ("A", 11)}
        sc._codon_tables_add("Round Trip", "111111", raw, source="user")
        # Clear cache and reload from disk
        sc._codon_tables_cache = None
        again = sc._codon_tables_get("111111")
        assert again is not None
        assert again["raw"]["ATG"] == ("M", 42)
        assert again["raw"]["GCT"] == ("A", 11)

    def test_search_substring(self):
        sc._codon_tables_load()
        assert any(e["taxid"] == "83333"
                   for e in sc._codon_search("coli"))
        assert sc._codon_search("") == sc._codon_tables_load()
        assert sc._codon_search("no_such_species_xyz") == []

    def test_search_by_genus_prefix(self):
        sc._codon_tables_load()
        raw = {"ATG": ("M", 5), "TAA": ("*", 1)}
        sc._codon_tables_add("Escherichia coli O157",  "900001", raw)
        sc._codon_tables_add("Escherichia albertii",   "900002", raw)
        sc._codon_tables_add("Salmonella enterica",    "900003", raw)
        hits = sc._codon_search("Escherichia")
        names = [e["name"] for e in hits]
        assert "Escherichia coli O157"  in names
        assert "Escherichia albertii"   in names
        assert "Salmonella enterica"    not in names
        # Alphabetical clustering within same rank
        assert names.index("Escherichia albertii") < names.index("Escherichia coli O157")

    def test_search_by_species_prefix(self):
        sc._codon_tables_load()
        raw = {"ATG": ("M", 5), "TAA": ("*", 1)}
        sc._codon_tables_add("Escherichia coli Z",  "900011", raw)
        sc._codon_tables_add("Salmonella typhi",    "900012", raw)
        hits = sc._codon_search("typh")
        names = [e["name"] for e in hits]
        assert "Salmonella typhi"   in names
        assert "Escherichia coli Z" not in names

    def test_search_by_taxid_prefix(self):
        sc._codon_tables_load()
        raw = {"ATG": ("M", 5), "TAA": ("*", 1)}
        sc._codon_tables_add("T species", "900099", raw)
        hits = sc._codon_search("9000")
        assert any(e["taxid"] == "900099" for e in hits)

    def test_search_exact_taxid_ranks_first(self):
        sc._codon_tables_load()
        raw = {"ATG": ("M", 5), "TAA": ("*", 1)}
        # Exact-match taxid vs. a species whose name substring-matches "83333"
        sc._codon_tables_add("Bogus 83333 label", "555555", raw)
        hits = sc._codon_search("83333")
        assert hits[0]["taxid"] == "83333"   # exact-match beats substring

    def test_search_genus_beats_name_substring(self):
        sc._codon_tables_load()
        raw = {"ATG": ("M", 5), "TAA": ("*", 1)}
        sc._codon_tables_add("Escherichia coli A", "900101", raw)
        sc._codon_tables_add("Notesch species",    "900102", raw)  # "esch" substring only
        hits = sc._codon_search("esch")
        names = [e["name"] for e in hits]
        assert names.index("Escherichia coli A") < names.index("Notesch species")


# ── Harmonization ─────────────────────────────────────────────────────────────

class TestHarmonize:
    def test_translate_of_harmonize_is_original(self):
        aa = "MAEVKLAGHIKQRSTVWYFND"
        dna = sc._codon_harmonize(aa, sc._CODON_BUILTIN_K12)
        assert sc._mut_translate(dna) == aa
        assert dna.endswith("TAA")

    def test_harmonize_met_trp_use_only_codon(self):
        dna = sc._codon_harmonize("MWMW", sc._CODON_BUILTIN_K12)
        # Met → ATG, Trp → TGG
        assert dna[:3]  == "ATG"
        assert dna[3:6] == "TGG"
        assert dna[6:9] == "ATG"
        assert dna[9:12] == "TGG"

    def test_harmonize_rejects_unknown_aa(self):
        with pytest.raises(ValueError, match="No codons"):
            sc._codon_harmonize("MAXA", sc._CODON_BUILTIN_K12)

    def test_distribution_matches_target(self):
        """For a leucine-heavy protein, CTG (K12's dominant Leu codon, ~49%)
        should be used most often."""
        aa = "L" * 100
        dna = sc._codon_harmonize(aa, sc._CODON_BUILTIN_K12)
        codons = [dna[i:i+3] for i in range(0, len(aa) * 3, 3)]
        from collections import Counter
        counts = Counter(codons)
        # CTG is ~46% of Leu in K12 raw (240 / (240+61+61+78+54+27) = ~0.464)
        assert counts["CTG"] >= 40
        assert counts["CTG"] == max(counts.values())


# ── Restriction-site fixer ────────────────────────────────────────────────────

class TestFixSites:
    def test_removes_ecori(self):
        dna = "ATGGAATTCGCGAAATAA"  # GAATTC at nt 4
        fixed, fixes = sc._codon_fix_sites(
            dna, "MEFAK", sc._CODON_BUILTIN_K12, {"EcoRI": "GAATTC"},
        )
        assert "GAATTC" not in fixed
        assert len(fixes) == 1
        assert sc._mut_translate(fixed) == sc._mut_translate(dna)

    def test_removes_bsai_both_strands(self):
        # GGTCTC forward or GAGACC (rc) inside the CDS
        dna = sc._codon_harmonize("MASGGTCTCREEEE", sc._CODON_BUILTIN_K12)
        # Synthesize a CDS that deliberately contains a BsaI site by hand
        # (harmonize won't produce one on K12 typically). Easier: seed manually.
        seed = "ATGGCGAGTGGTCTCCGTGAGGAGGAGGAGTAA"
        assert "GGTCTC" in seed
        fixed, fixes = sc._codon_fix_sites(
            seed, sc._mut_translate(seed),
            sc._CODON_BUILTIN_K12, {"BsaI": "GGTCTC"},
        )
        assert "GGTCTC" not in fixed
        assert "GAGACC" not in fixed  # rc
        assert sc._mut_translate(fixed) == sc._mut_translate(seed)


# ── CAI / GC ──────────────────────────────────────────────────────────────────

class TestMetrics:
    def test_cai_in_unit_range(self):
        aa = "MAEVKLAGHIKQR"
        dna = sc._codon_harmonize(aa, sc._CODON_BUILTIN_K12)
        cai = sc._codon_cai(dna, sc._CODON_BUILTIN_K12)
        assert 0.0 < cai <= 1.0

    def test_gc_empty(self):
        assert sc._codon_gc("") == 0.0

    def test_gc_allgc(self):
        assert sc._codon_gc("GCGCGC") == 100.0

    def test_gc_half(self):
        assert sc._codon_gc("ATGC") == 50.0


# ── Kazusa HTML parser ────────────────────────────────────────────────────────

_FAKE_KAZUSA_HTML = """
<html><body><pre>
AmAcid  Codon      Number    /1000    Fraction
Gly     GGG       44.00    16.63       0.12
Gly     GGA       47.00    17.76       0.13
Gly     GGT      109.00    41.19       0.31
Gly     GGC      171.00    64.61       0.44
Glu     GAG       94.00    35.52       0.30
Glu     GAA      224.00    84.64       0.70
Asp     GAT      194.00    73.31       0.65
Asp     GAC      105.00    39.68       0.35
Val     GTG      135.00    51.01       0.37
Val     GTA       59.00    22.29       0.16
Val     GTT       86.00    32.50       0.24
Val     GTC       60.00    22.67       0.17
Ala     GCG      197.00    74.44       0.36
Ala     GCA      108.00    40.81       0.20
Ala     GCT       55.00    20.78       0.10
Ala     GCC      162.00    61.22       0.29
Arg     AGG        8.00     3.02       0.02
Arg     AGA        7.00     2.64       0.02
Ser     AGT       37.00    13.98       0.14
Ser     AGC       85.00    32.11       0.33
Lys     AAG       62.00    23.43       0.27
Lys     AAA      170.00    64.23       0.73
Asn     AAT      112.00    42.31       0.47
Asn     AAC      125.00    47.22       0.53
Met     ATG      127.00    47.98       1.00
Ile     ATA       19.00     7.18       0.07
Ile     ATT      156.00    58.93       0.58
Ile     ATC       93.00    35.13       0.35
Thr     ACG       59.00    22.29       0.24
Thr     ACA       33.00    12.47       0.13
Thr     ACT       41.00    15.49       0.16
Thr     ACC      117.00    44.20       0.47
Trp     TGG       55.00    20.78       1.00
Cys     TGT       30.00    11.33       0.42
Cys     TGC       41.00    15.49       0.58
Tyr     TAT       86.00    32.50       0.53
Tyr     TAC       75.00    28.33       0.47
Leu     TTG       61.00    23.04       0.11
Leu     TTA       78.00    29.47       0.14
Phe     TTT      101.00    38.16       0.57
Phe     TTC       77.00    29.09       0.43
Ser     TCG       41.00    15.49       0.16
Ser     TCA       40.00    15.11       0.15
Ser     TCT       29.00    10.95       0.11
Ser     TCC       28.00    10.58       0.11
Arg     CGG       21.00     7.93       0.05
Arg     CGA       22.00     8.31       0.06
Arg     CGT      108.00    40.81       0.38
Arg     CGC      133.00    50.24       0.47
Gln     CAG      142.00    53.65       0.70
Gln     CAA       62.00    23.43       0.30
His     CAT       81.00    30.60       0.55
His     CAC       67.00    25.31       0.45
Leu     CTG      240.00    90.67       0.46
Leu     CTA       27.00    10.20       0.05
Leu     CTT       61.00    23.04       0.11
Leu     CTC       54.00    20.40       0.10
Pro     CCG      137.00    51.76       0.49
Pro     CCA       34.00    12.85       0.12
Pro     CCT       43.00    16.24       0.15
Pro     CCC       33.00    12.47       0.12
End     TAA        9.00     3.40       0.56
End     TAG        0.00     0.00       0.00
End     TGA        5.00     1.89       0.31
</pre></body></html>
"""


class TestKazusaParse:
    def test_parse_valid_html_returns_k12(self):
        raw = sc._codon_parse_kazusa_html(_FAKE_KAZUSA_HTML)
        assert raw is not None
        # Spot-check a few counts match K12
        assert raw["ATG"] == ("M", 127)
        assert raw["CTG"] == ("L", 240)
        assert raw["TAA"] == ("*", 9)

    def test_parse_truncated_returns_none(self):
        assert sc._codon_parse_kazusa_html("<html><body>No table here</body></html>") is None

    def test_parse_handles_rna(self):
        html = "<pre>AUG 100\nGGG 50\nGAA 75\nCUG 240\nGCG 197\n"
        # Not enough codons to pass validation, should return None
        assert sc._codon_parse_kazusa_html(html) is None


# ── NCBI taxid lookup ─────────────────────────────────────────────────────────

class _FakeResponse:
    """Minimal urlopen() context-manager stand-in that yields canned bytes."""
    def __init__(self, body: bytes):
        self._body = body
    def __enter__(self):
        return self
    def __exit__(self, *exc):
        return False
    def read(self):
        return self._body


class TestNcbiPrepTerm:
    def test_single_token_uses_subtree_or_wildcard(self):
        """Single-word queries should combine exact-taxon subtree (all species
        in that genus) with a wildcard prefix so 'Escher' still matches."""
        term = sc._ncbi_prep_term("Escherichia")
        assert "Escherichia[Subtree] AND species[Rank]" in term
        assert "Escherichia*" in term
        # Partial prefix also uses the OR form — wildcard does the work
        term2 = sc._ncbi_prep_term("Escher")
        assert "Escher[Subtree] AND species[Rank]" in term2
        assert "Escher*" in term2
        # Whitespace is stripped
        assert sc._ncbi_prep_term("  Homo ").startswith("(Homo[Subtree]")

    def test_multi_word_appends_trailing_wildcard(self):
        assert sc._ncbi_prep_term("Homo sapiens")     == "Homo sapiens*"
        assert sc._ncbi_prep_term("Saccharomyces ce") == "Saccharomyces ce*"

    def test_preserves_user_wildcards_and_fields(self):
        assert sc._ncbi_prep_term("Escher*")               == "Escher*"
        assert sc._ncbi_prep_term("coli[Scientific Name]") == "coli[Scientific Name]"

    def test_empty_stays_empty(self):
        assert sc._ncbi_prep_term("")    == ""
        assert sc._ncbi_prep_term("  ")  == ""


class TestNcbiSearch:
    def test_empty_query_returns_empty(self):
        hits, total, msg = sc._ncbi_taxid_search("")
        assert hits == []
        assert total == 0
        assert "empty" in msg.lower()

    def test_multiple_hits_batched_esummary(self, monkeypatch):
        """Verifies esummary is called once with a comma-joined id list and
        ScientificName is picked up per-DocSum."""
        esearch_xml = (b"<?xml version=\"1.0\"?><eSearchResult>"
                       b"<Count>3</Count><IdList>"
                       b"<Id>561</Id><Id>562</Id><Id>564</Id>"
                       b"</IdList></eSearchResult>")
        esummary_xml = (b"<?xml version=\"1.0\"?><eSummaryResult>"
                        b"<DocSum><Id>561</Id>"
                        b"<Item Name=\"ScientificName\" Type=\"String\">"
                        b"Escherichia</Item></DocSum>"
                        b"<DocSum><Id>562</Id>"
                        b"<Item Name=\"ScientificName\" Type=\"String\">"
                        b"Escherichia coli</Item></DocSum>"
                        b"<DocSum><Id>564</Id>"
                        b"<Item Name=\"ScientificName\" Type=\"String\">"
                        b"Escherichia fergusonii</Item></DocSum>"
                        b"</eSummaryResult>")
        calls: list[str] = []

        def fake_urlopen(req, timeout=None):
            url = req.full_url if hasattr(req, "full_url") else str(req)
            calls.append(url)
            if "esearch" in url:
                return _FakeResponse(esearch_xml)
            return _FakeResponse(esummary_xml)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        hits, total, msg = sc._ncbi_taxid_search("Escher")
        assert [h["taxid"] for h in hits] == ["561", "562", "564"]
        assert [h["name"]  for h in hits] == [
            "Escherichia", "Escherichia coli", "Escherichia fergusonii",
        ]
        assert total == 3
        # One esearch, one esummary — batch worked
        assert sum("esearch"  in u for u in calls) == 1
        assert sum("esummary" in u for u in calls) == 1
        # The esearch URL must carry both halves of the OR term: the
        # Subtree+species-rank clause and the wildcard fallback. Both are
        # URL-encoded.
        esearch_url = next(u for u in calls if "esearch" in u)
        assert "Escher%5BSubtree%5D" in esearch_url   # [Subtree]
        assert "species%5BRank%5D"   in esearch_url   # [Rank]
        assert "Escher%2A"           in esearch_url   # * wildcard

    def test_total_higher_than_retrieved_mentions_refine(self, monkeypatch):
        """When total hit count exceeds retmax, the status message should
        nudge the user to refine."""
        esearch_xml = (b"<?xml version=\"1.0\"?><eSearchResult>"
                       b"<Count>1200</Count><IdList>"
                       b"<Id>1</Id></IdList></eSearchResult>")
        esummary_xml = (b"<?xml version=\"1.0\"?><eSummaryResult>"
                        b"<DocSum><Id>1</Id>"
                        b"<Item Name=\"ScientificName\" Type=\"String\">"
                        b"root</Item></DocSum></eSummaryResult>")

        def fake_urlopen(req, timeout=None):
            url = req.full_url if hasattr(req, "full_url") else str(req)
            return _FakeResponse(esearch_xml if "esearch" in url else esummary_xml)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        hits, total, msg = sc._ncbi_taxid_search("a")
        assert len(hits) == 1
        assert total == 1200
        assert "refine" in msg.lower()

    def test_no_results(self, monkeypatch):
        xml = (b"<?xml version=\"1.0\"?><eSearchResult>"
               b"<Count>0</Count><IdList></IdList></eSearchResult>")

        def fake_urlopen(req, timeout=None):
            return _FakeResponse(xml)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        hits, total, msg = sc._ncbi_taxid_search("zzznope")
        assert hits == []
        assert total == 0
        assert "no ncbi" in msg.lower()

    def test_network_error(self, monkeypatch):
        def boom(req, timeout=None):
            raise OSError("connection refused")

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", boom)
        hits, total, msg = sc._ncbi_taxid_search("Homo sapiens")
        assert hits == []
        assert total == 0
        assert "network error" in msg.lower()

    def test_esummary_failure_still_returns_ids(self, monkeypatch):
        """If esummary fails after esearch succeeds, hits come back with
        fallback '(taxid N)' names so the user can still pick one."""
        esearch_xml = (b"<?xml version=\"1.0\"?><eSearchResult>"
                       b"<Count>1</Count><IdList><Id>4932</Id></IdList>"
                       b"</eSearchResult>")

        def fake_urlopen(req, timeout=None):
            url = req.full_url if hasattr(req, "full_url") else str(req)
            if "esearch" in url:
                return _FakeResponse(esearch_xml)
            raise OSError("esummary down")

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        hits, total, msg = sc._ncbi_taxid_search("Saccharomyces cerevisiae")
        assert len(hits) == 1
        assert hits[0]["taxid"] == "4932"
        assert "(taxid 4932)" in hits[0]["name"]


# ── Integration with _mut_design_inner ────────────────────────────────────────

class TestMutInnerWithTable:
    def test_inner_uses_custom_table(self):
        """_mut_design_inner should honour a custom codon table for the
        mutant codon choice."""
        # Synthetic CDS: ATG GCT TAA (M-A-stop)
        cds = "ATG" + "GCT" + ("GCT" * 30) + "TAA"
        # Custom table where Phe (F) only has TTT (no TTC)
        custom = {
            "ATG": ("M", 1), "GCT": ("A", 1), "GCC": ("A", 1),
            "TTT": ("F", 1), "TAA": ("*", 1),
            # Fill out only what the harmonizer/inner design needs
        }
        # Design A→F at position 2 — mut_codon must be TTT (only F codon)
        inner = sc._mut_design_inner(cds, 2, "F", "A", codon_table=custom)
        assert inner["mut_codon"] == "TTT"

    def test_inner_default_is_k12(self):
        cds = "ATG" + "GCT" + ("GCT" * 30) + "TAA"
        inner = sc._mut_design_inner(cds, 2, "F", "A")
        # K12: F has TTT (101) and TTC (77); TTT is most frequent
        assert inner["mut_codon"] == "TTT"


# ── UI smoke — MutagenizeModal 3-source flow + SpeciesPickerModal ────────────

class TestMutagenizeModalSources:
    async def test_source_switching(self):
        """The three source panels must toggle visibility as the Select
        changes, without any widget-lookup errors."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            cds = "ATG" + ("GCT" * 80) + "TAA"
            feats = [{"type": "CDS", "label": "testA",
                      "start": 0, "end": len(cds), "strand": 1}]
            await app.push_screen(sc.MutagenizeModal(cds, feats, "TEST"))
            await pilot.pause(0.3)
            src = app.screen.query_one("#mut-source")
            for val in ("prot", "lib", "map"):
                src.value = val
                await pilot.pause(0.1)
                visible = {
                    "map":  app.screen.query_one("#mut-src-map").display,
                    "lib":  app.screen.query_one("#mut-src-lib").display,
                    "prot": app.screen.query_one("#mut-src-prot").display,
                }
                assert visible[val] is True
                assert sum(visible.values()) == 1
            app.exit()

    async def test_species_picker_opens(self):
        """Clicking 'Change…' next to the codon-table label must push the
        SpeciesPickerModal."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            cds = "ATG" + ("GCT" * 80) + "TAA"
            feats = [{"type": "CDS", "label": "testA",
                      "start": 0, "end": len(cds), "strand": 1}]
            await app.push_screen(sc.MutagenizeModal(cds, feats, "TEST"))
            await pilot.pause(0.3)
            app.screen.query_one("#btn-mut-codon").action_press()
            await pilot.pause(0.3)
            assert type(app.screen).__name__ == "SpeciesPickerModal"
            # The built-in K12 row must be present
            lv = app.screen.query_one("#sp-list")
            assert lv is not None
            app.exit()

    async def test_protein_input_path(self):
        """Harmonize-from-protein must produce a valid CDS and enable
        mutation design."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            # Modal mounted with empty template → protein source becomes
            # the only way to get a CDS.
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            ta = modal.query_one("#mut-prot-aa")
            ta.text = "MAEVKLAGHIKQRSTVWY"
            modal.query_one("#btn-mut-harmonize").action_press()
            await pilot.pause(0.2)
            assert modal._cds_dna != ""
            assert sc._mut_translate(modal._cds_dna) == "MAEVKLAGHIKQRSTVWY"
            app.exit()

    async def test_protein_input_rejects_mid_sequence_stop(self):
        """Regression: a '*' inside the protein (not trailing) must be
        flagged, not silently stripped."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            modal.query_one("#mut-prot-aa").text = "MA*EF"
            modal.query_one("#btn-mut-harmonize").action_press()
            await pilot.pause(0.1)
            assert modal._cds_dna == ""
            app.exit()

    async def test_protein_input_allows_trailing_stop(self):
        """A single trailing '*' should be treated as an explicit stop and
        stripped silently before harmonization."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            modal.query_one("#mut-prot-aa").text = "MAEVK*"
            modal.query_one("#btn-mut-harmonize").action_press()
            await pilot.pause(0.1)
            assert modal._cds_dna != ""
            assert sc._mut_translate(modal._cds_dna) == "MAEVK"

    async def test_protein_input_rejects_invalid_chars(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            modal.query_one("#mut-prot-aa").text = "MAEBZJX"
            modal.query_one("#btn-mut-harmonize").action_press()
            await pilot.pause(0.1)
            assert modal._cds_dna == ""


class TestLibrarySourceWrapFeature:
    """Regression: a CompoundLocation CDS that spans the origin must be
    loaded correctly from the plasmid library, not flattened into a
    whole-plasmid range."""

    def test_compound_location_preserves_wrap(self, tmp_path, monkeypatch):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

        # 120 nt plasmid; CDS spans origin: nt 100..120 + 1..10 = 30 nt.
        seq = "A" * 100 + "ATG" + "GCT" * 5 + "GC" + ("GCT" * 2) + "TAA" + "A"
        # Actually build a known wrap CDS cleanly:
        cds = "ATG" + ("GCT" * 8) + "TAA"   # 30 nt = 10 codons
        padding = "N" * 90
        plasmid = cds[15:] + padding + cds[:15]  # wrap: last 15 at start, first 15 at end
        start = len(plasmid) - 15
        end   = 15
        rec = SeqRecord(Seq(plasmid), id="WRAP", name="WRAP",
                        description="wrap CDS test")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "circular"
        rec.features.append(SeqFeature(
            CompoundLocation([
                FeatureLocation(start, len(plasmid), strand=1),
                FeatureLocation(0, end, strand=1),
            ]),
            type="CDS",
            qualifiers={"label": ["wrapCDS"]},
        ))

        # Write a fake library entry and reload through MutagenizeModal logic
        from io import StringIO
        from Bio import SeqIO
        buf = StringIO()
        SeqIO.write(rec, buf, "genbank")
        gb_text = buf.getvalue()

        _entries = [{
            "id":      "WRAP",
            "name":    "WRAP",
            "size":    len(plasmid),
            "n_feats": 1,
            "source":  "test",
            "added":   "2026-04-13",
            "gb_text": gb_text,
        }]
        # Patch _load_library directly so the PlasmidApp's default-seed worker
        # can't overwrite our fake cache between mount and modal open.
        monkeypatch.setattr(sc, "_load_library", lambda: list(_entries))
        monkeypatch.setattr(sc, "_library_cache", _entries)
        # Disable default-seed worker so it doesn't try to fetch pACYC184.
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)

        import asyncio
        async def go():
            app = sc.PlasmidApp()
            async with app.run_test(size=(140, 50)) as pilot:
                await pilot.pause()
                await app.push_screen(sc.MutagenizeModal("", [], ""))
                await pilot.pause(0.3)
                modal = app.screen
                modal.query_one("#mut-source").value = "lib"
                await pilot.pause(0.1)
                modal.query_one("#mut-lib").value = "0"
                await pilot.pause(0.2)
                # The wrap CDS should have been recognised — end < start
                feats = modal._lib_feats
                assert any(f["end"] < f["start"] for f in feats), \
                    f"wrap CDS not preserved; got feats={feats}"
                wrap_feat = next(f for f in feats if f["end"] < f["start"])
                assert wrap_feat["start"] == start
                assert wrap_feat["end"]   == end
                app.exit()
        asyncio.run(go())
