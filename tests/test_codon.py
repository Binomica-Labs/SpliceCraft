"""
test_codon — codon usage registry + optimization pipeline.

Covers the persistent _codon_tables_* registry, the pure-function
_codon_optimize / _codon_fix_sites / _codon_cai / _codon_gc helpers, and
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


# ── Optimization ─────────────────────────────────────────────────────────────

class TestOptimize:
    def test_translate_of_optimize_is_original(self):
        aa = "MAEVKLAGHIKQRSTVWYFND"
        dna = sc._codon_optimize(aa, sc._CODON_BUILTIN_K12)
        assert sc._mut_translate(dna) == aa
        assert dna.endswith("TAA")

    def test_optimize_met_trp_use_only_codon(self):
        dna = sc._codon_optimize("MWMW", sc._CODON_BUILTIN_K12)
        # Met → ATG, Trp → TGG
        assert dna[:3]  == "ATG"
        assert dna[3:6] == "TGG"
        assert dna[6:9] == "ATG"
        assert dna[9:12] == "TGG"

    def test_optimize_rejects_unknown_aa(self):
        with pytest.raises(ValueError, match="No codons"):
            sc._codon_optimize("MAXA", sc._CODON_BUILTIN_K12)

    def test_distribution_matches_target(self):
        """For a leucine-heavy protein, CTG (K12's dominant Leu codon, ~49%)
        should be used most often."""
        aa = "L" * 100
        dna = sc._codon_optimize(aa, sc._CODON_BUILTIN_K12)
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
        # GGTCTC forward or GAGACC (rc) inside the CDS.
        # Synthesize a CDS that deliberately contains a BsaI site by hand
        # (optimize won't produce one on K12 typically). Easier: seed manually.
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
        dna = sc._codon_optimize(aa, sc._CODON_BUILTIN_K12)
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
    """Minimal urlopen() context-manager stand-in that yields canned bytes.

    Mirrors `http.client.HTTPResponse.read(amt=None)`: when the caller
    passes a max-byte hint (the response-size cap pattern), return at
    most that many bytes.
    """
    def __init__(self, body: bytes):
        self._body = body
    def __enter__(self):
        return self
    def __exit__(self, *exc):
        return False
    def read(self, amt=None):
        if amt is None:
            return self._body
        return self._body[:amt]


# Real NCBI eutils responses open with an EXTERNAL-DTD DOCTYPE (a PUBLIC id +
# .dtd URL, no internal `[...]` subset). `_safe_xml_parse` refuses ANY DTD
# unless the caller passes allow_dtd=True, so mocks that OMITTED the DOCTYPE
# silently diverged from production — which is exactly how the taxon-search
# "XML parser error" (2026-05-30) shipped despite a green suite. Every NCBI
# mock below now carries the real preamble so the tests exercise the same
# parse path the live app hits. See `_ncbi_taxid_search` / [PIT-19].
_ESEARCH_PREAMBLE = (
    b'<?xml version="1.0" encoding="UTF-8" ?>\n'
    b'<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" '
    b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">\n')
_ESUMMARY_PREAMBLE = (
    b'<?xml version="1.0" encoding="UTF-8" ?>\n'
    b'<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD esummary v1 20041029//EN" '
    b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20041029/esummary-v1.dtd">\n')


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
        esearch_xml = (_ESEARCH_PREAMBLE + b"<eSearchResult>"
                       b"<Count>3</Count><IdList>"
                       b"<Id>561</Id><Id>562</Id><Id>564</Id>"
                       b"</IdList></eSearchResult>")
        esummary_xml = (_ESUMMARY_PREAMBLE + b"<eSummaryResult>"
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
        esearch_xml = (_ESEARCH_PREAMBLE + b"<eSearchResult>"
                       b"<Count>1200</Count><IdList>"
                       b"<Id>1</Id></IdList></eSearchResult>")
        esummary_xml = (_ESUMMARY_PREAMBLE + b"<eSummaryResult>"
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
        xml = (_ESEARCH_PREAMBLE + b"<eSearchResult>"
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
        esearch_xml = (_ESEARCH_PREAMBLE + b"<eSearchResult>"
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

    def test_oversized_response_rejected(self, monkeypatch):
        """Regression guard for 2026-05-06 fix: a compromised / MITM'd
        upstream that streams gigabytes at us must not OOM the worker.
        Response is capped at `_NCBI_MAX_RESPONSE_BYTES`."""
        monkeypatch.setattr(sc, "_NCBI_MAX_RESPONSE_BYTES", 100)
        # Build a response well over the cap.
        big_xml = b"<?xml?><eSearchResult>" + (b"<Id>1</Id>" * 1000) + b"</eSearchResult>"

        def fake_urlopen(req, timeout=None):
            return _FakeResponse(big_xml)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        hits, total, msg = sc._ncbi_taxid_search("a")
        assert hits == []
        assert total == 0
        assert "oversized" in msg.lower()

    def test_external_dtd_doctype_response_parses(self, monkeypatch):
        """Regression (2026-05-30): real NCBI esearch/esummary responses open
        with an external-DTD DOCTYPE. `_ncbi_taxid_search` must parse them
        (allow_dtd=True), else every live taxon lookup fails with
        'Could not parse NCBI response: XML contains DTD/ENTITY — refusing to
        parse' — the codon-table picker's "XML parser error". The shared
        mocks above now all carry the DOCTYPE, but assert the happy path here
        explicitly so the guard can't be lost in a future mock refactor."""
        esearch = _ESEARCH_PREAMBLE + (
            b"<eSearchResult><Count>1</Count>"
            b"<IdList><Id>9606</Id></IdList></eSearchResult>")
        esummary = _ESUMMARY_PREAMBLE + (
            b"<eSummaryResult><DocSum><Id>9606</Id>"
            b"<Item Name=\"ScientificName\" Type=\"String\">"
            b"Homo sapiens</Item></DocSum></eSummaryResult>")

        def fake_urlopen(req, timeout=None):
            url = req.full_url if hasattr(req, "full_url") else str(req)
            return _FakeResponse(esearch if "esearch" in url else esummary)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        hits, total, msg = sc._ncbi_taxid_search("Homo sapiens")
        assert [h["taxid"] for h in hits] == ["9606"]
        assert hits[0]["name"] == "Homo sapiens"
        assert total == 1
        # The defining symptom of the bug was a parse-error status string.
        assert "parse" not in msg.lower()
        assert "dtd" not in msg.lower()


class TestNcbiTaxidSearchTerms:
    """Cascading term builder mirrored from ScriptoScope's genbank_search."""

    def test_empty(self):
        assert sc._ncbi_taxid_search_terms("") == []
        assert sc._ncbi_taxid_search_terms("   ") == []

    def test_single_token_one_strategy(self):
        # Single token has no broader form beyond `_ncbi_prep_term`'s
        # `OR {t}*` clause — so exactly one strategy, no cascade.
        terms = sc._ncbi_taxid_search_terms("Escherichia")
        assert terms == [sc._ncbi_prep_term("Escherichia")]

    def test_multi_word_cascade_strict_then_and_then_or(self):
        terms = sc._ncbi_taxid_search_terms("Homo sapiens")
        assert terms == [
            "Homo sapiens*",          # strict (== _ncbi_prep_term)
            "Homo* AND sapiens*",     # broaden 1: partial genus + species
            "Homo* OR sapiens*",      # broaden 2: any token → related taxa
        ]

    def test_strict_strategy_always_equals_prep_term(self):
        for q in ("Homo sapiens", "Saccharomyces cerevisiae", "Escher", "E"):
            assert sc._ncbi_taxid_search_terms(q)[0] == sc._ncbi_prep_term(q)

    def test_user_wildcard_or_field_passthrough_no_cascade(self):
        assert sc._ncbi_taxid_search_terms("Escher*") == ["Escher*"]
        assert sc._ncbi_taxid_search_terms("coli[Scientific Name]") == \
            ["coli[Scientific Name]"]


class TestNcbiSearchCascade:
    """Lax/broadening search behaviour: a strict query that returns
    nothing falls through to broader rounds (AND-wildcard, then
    OR-wildcard) so a related / imprecise name still surfaces taxa —
    while a network error aborts immediately."""

    _EMPTY = (_ESEARCH_PREAMBLE + b'<eSearchResult>'
              b'<Count>0</Count><IdList></IdList></eSearchResult>')

    @staticmethod
    def _hit(taxid: str) -> bytes:
        return _ESEARCH_PREAMBLE + (
            f'<eSearchResult><Count>1</Count>'
            f'<IdList><Id>{taxid}</Id></IdList></eSearchResult>'
            ).encode()

    @staticmethod
    def _summ(taxid: str, name: str) -> bytes:
        return _ESUMMARY_PREAMBLE + (
            f'<eSummaryResult><DocSum>'
            f'<Id>{taxid}</Id><Item Name="ScientificName" '
            f'Type="String">{name}</Item></DocSum></eSummaryResult>'
            ).encode()

    @staticmethod
    def _url(req):
        return req.full_url if hasattr(req, "full_url") else str(req)

    def test_broadens_to_or_when_stricter_rounds_empty(self, monkeypatch):
        from urllib.parse import unquote_plus
        calls: list[str] = []

        def fake(req, timeout=None):
            url = self._url(req)
            calls.append(url)
            if "esummary" in url:
                return _FakeResponse(self._summ("9606", "Homo sapiens"))
            # Only the OR-broadened round finds anything.
            if " OR " in unquote_plus(url):
                return _FakeResponse(self._hit("9606"))
            return _FakeResponse(self._EMPTY)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        # Genus typo — rescued by the species epithet in the OR round.
        hits, total, msg = sc._ncbi_taxid_search("Homon sapiens")
        assert [h["taxid"] for h in hits] == ["9606"]
        assert "broadened" in msg.lower()
        assert sum("esearch" in u for u in calls) == 3   # strict + AND + OR

    def test_broadens_to_and(self, monkeypatch):
        from urllib.parse import unquote_plus
        calls: list[str] = []

        def fake(req, timeout=None):
            url = self._url(req)
            calls.append(url)
            if "esummary" in url:
                return _FakeResponse(self._summ("4932", "Saccharomyces cerevisiae"))
            if " AND " in unquote_plus(url):
                return _FakeResponse(self._hit("4932"))
            return _FakeResponse(self._EMPTY)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        hits, total, msg = sc._ncbi_taxid_search("Sacchar cerev")
        assert [h["taxid"] for h in hits] == ["4932"]
        assert "broadened" in msg.lower()
        assert sum("esearch" in u for u in calls) == 2   # strict + AND only

    def test_strict_hit_stops_cascade_and_omits_broadened_hint(self, monkeypatch):
        calls: list[str] = []

        def fake(req, timeout=None):
            url = self._url(req)
            calls.append(url)
            if "esummary" in url:
                return _FakeResponse(self._summ("9606", "Homo sapiens"))
            return _FakeResponse(self._hit("9606"))   # every round would hit

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        hits, total, msg = sc._ncbi_taxid_search("Homo sapiens")
        assert [h["taxid"] for h in hits] == ["9606"]
        assert "broadened" not in msg.lower()
        assert sum("esearch" in u for u in calls) == 1   # stopped at strict

    def test_multi_word_all_empty_exhausts_cascade(self, monkeypatch):
        calls: list[str] = []

        def fake(req, timeout=None):
            calls.append(self._url(req))
            return _FakeResponse(self._EMPTY)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake)
        hits, total, msg = sc._ncbi_taxid_search("zzz nope")
        assert hits == []
        assert total == 0
        assert "no ncbi" in msg.lower()
        assert sum("esearch" in u for u in calls) == 3   # tried all 3

    def test_network_error_aborts_cascade_on_first_failure(self, monkeypatch):
        calls: list[str] = []

        def boom(req, timeout=None):
            calls.append(self._url(req))
            raise OSError("connection refused")

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", boom)
        hits, total, msg = sc._ncbi_taxid_search("Homon sapiens")
        assert hits == []
        assert "network error" in msg.lower()
        assert len(calls) == 1   # did NOT retry broader strategies


class TestKazusaSizeCap:
    def test_kazusa_oversized_response_rejected(self, monkeypatch):
        """Regression guard for 2026-05-06 fix: cap Kazusa's HTML
        response size to bound worker memory."""
        monkeypatch.setattr(sc, "_KAZUSA_MAX_RESPONSE_BYTES", 100)
        big_html = b"<html>" + (b"X" * 10_000) + b"</html>"

        def fake_urlopen(req, timeout=None):
            return _FakeResponse(big_html)

        import urllib.request
        monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
        raw, msg = sc._codon_fetch_kazusa("83333")
        assert raw is None
        assert "oversized" in msg.lower()


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
            # Fill out only what the optimizer/inner design needs
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
        """The four source panels must toggle visibility as the Select
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
            for val in ("prot", "lib", "parts", "map"):
                src.value = val
                await pilot.pause(0.1)
                visible = {
                    "map":   app.screen.query_one("#mut-src-map").display,
                    "lib":   app.screen.query_one("#mut-src-lib").display,
                    "parts": app.screen.query_one("#mut-src-parts").display,
                    "prot":  app.screen.query_one("#mut-src-prot").display,
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
        """Optimize-from-protein must produce a valid CDS and enable
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
            modal.query_one("#btn-mut-optimize").action_press()
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
            modal.query_one("#btn-mut-optimize").action_press()
            await pilot.pause(0.1)
            assert modal._cds_dna == ""
            app.exit()

    async def test_protein_input_allows_trailing_stop(self):
        """A single trailing '*' should be treated as an explicit stop and
        stripped silently before optimization."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            modal.query_one("#mut-prot-aa").text = "MAEVK*"
            modal.query_one("#btn-mut-optimize").action_press()
            await pilot.pause(0.1)
            assert modal._cds_dna != ""
            assert sc._mut_translate(modal._cds_dna) == "MAEVK"

    async def test_protein_stops_selector_appends_triple(self):
        """The stops selector controls how many stop codons are appended
        when the protein carries no trailing '*'."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            modal.query_one("#mut-prot-aa").text = "MAEV"
            modal.query_one("#mut-stops").value = "3"
            await pilot.pause(0.1)
            modal.query_one("#btn-mut-optimize").action_press()
            await pilot.pause(0.2)
            assert len(modal._cds_dna) == 4 * 3 + 3 * 3      # body + 3 stops
            assert sc._mut_translate(modal._cds_dna) == "MAEV"

    async def test_protein_double_trailing_stop_overrides_selector(self):
        """A trailing '**' run is honored as two stops even when the
        selector says 1."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            modal.query_one("#mut-prot-aa").text = "MAEV**"
            modal.query_one("#mut-stops").value = "1"
            await pilot.pause(0.1)
            modal.query_one("#btn-mut-optimize").action_press()
            await pilot.pause(0.2)
            assert len(modal._cds_dna) == 4 * 3 + 3 * 2      # body + 2 stops
            assert sc._mut_translate(modal._cds_dna) == "MAEV"

    async def test_protein_rejects_more_than_three_stops(self):
        """More than three trailing stop codons is rejected, not optimized."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            modal.query_one("#mut-prot-aa").text = "MAEV****"
            modal.query_one("#btn-mut-optimize").action_press()
            await pilot.pause(0.1)
            assert modal._cds_dna == ""

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
            modal.query_one("#btn-mut-optimize").action_press()
            await pilot.pause(0.1)
            assert modal._cds_dna == ""


class TestMutagenizeNoPlasmid:
    """Mutagenize must launch with no loaded plasmid. The 'Current map
    features' option is excluded from the source dropdown in that case
    so the user only sees sources that can actually produce a CDS."""

    async def test_modal_opens_with_empty_template(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            src = modal.query_one("#mut-source")
            values = [v for (_label, v) in src._options]  # type: ignore[attr-defined]
            assert "map" not in values
            # The other three sources must all be present.
            assert "lib" in values and "parts" in values and "prot" in values
            # Default starts at the first remaining source — "lib" — and
            # the lib panel is the visible one.
            assert src.value == "lib"
            assert modal.query_one("#mut-src-lib").display is True
            assert modal.query_one("#mut-src-map").display is False
            app.exit()

    async def test_action_open_mutagenize_no_record(self, monkeypatch):
        """PlasmidApp.action_open_mutagenize must push the modal even when
        _current_record is None (regression: previously errored out with a
        notify-and-return). Suppress the default-seed worker so it can't
        load pACYC184 in the background and mask the no-record state."""
        monkeypatch.setattr(sc.PlasmidApp, "_seed_default_library",
                            lambda self: None)
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            assert app._current_record is None
            app.action_open_mutagenize()
            await pilot.pause(0.3)
            assert type(app.screen).__name__ == "MutagenizeModal"
            app.exit()

    async def test_modal_includes_map_when_template_present(self):
        cds = "ATG" + ("GCT" * 30) + "TAA"
        feats = [{"type": "CDS", "label": "x",
                  "start": 0, "end": len(cds), "strand": 1}]
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal(cds, feats, "x"))
            await pilot.pause(0.3)
            modal = app.screen
            src = modal.query_one("#mut-source")
            values = [v for (_label, v) in src._options]  # type: ignore[attr-defined]
            assert "map" in values
            assert src.value == "map"
            app.exit()


class TestPartsBinSource:
    """The Parts Bin source feeds a stored part's `sequence` (5'→3' insert,
    no tails) into Mutagenize as a single-CDS pseudo-plasmid.

    Eligibility filter mirrors the map/library sources: ≥ 30 bp and a
    multiple of 3."""

    def _ok_part(self, name: str = "myCds") -> dict:
        cds = "ATG" + ("GCT" * 30) + "TAA"   # 96 nt = 32 codons
        return {
            "name": name, "type": "CDS", "position": "Pos 3",
            "oh5": "AGGT", "oh3": "GCTT", "sequence": cds,
            "grammar": "gb_l0",
        }

    async def test_options_filter_short_and_offgrid(self, monkeypatch):
        bin_entries = [
            self._ok_part("good_cds"),
            {"name": "tooShort", "type": "CDS", "sequence": "ATGAAATAA"},  # 9 bp
            {"name": "offgrid",  "type": "CDS",
             "sequence": "ATG" + ("GCT" * 30) + "T"},                      # 94 bp, %3!=0
            {"name": "empty",    "type": "Promoter", "sequence": ""},
        ]
        monkeypatch.setattr(sc, "_load_parts_bin", lambda: list(bin_entries))
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            opts = modal._build_parts_options()
            # Only the one good part passes.
            assert len(opts) == 1
            label, val = opts[0]
            assert "good_cds" in label
            assert val == "0"
            app.exit()

    async def test_options_empty_bin(self, monkeypatch):
        monkeypatch.setattr(sc, "_load_parts_bin", lambda: [])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            opts = modal._build_parts_options()
            assert len(opts) == 1
            assert opts[0][1] == "_none"
            app.exit()

    async def test_select_loads_cds(self, monkeypatch):
        """Picking a part from the dropdown loads its `sequence` as the
        CDS, sets meta with origin='parts', and enables the preview."""
        part = self._ok_part("aeBlue")
        monkeypatch.setattr(sc, "_load_parts_bin", lambda: [part])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "parts"
            await pilot.pause(0.1)
            # Re-build options against the patched parts bin: compose ran
            # before our monkeypatch took effect, so the dropdown was seeded
            # with the (empty) real bin. set_options refreshes from the patch.
            modal.query_one("#mut-parts").set_options(modal._build_parts_options())
            await pilot.pause(0.1)
            modal.query_one("#mut-parts").value = "0"
            await pilot.pause(0.2)
            assert modal._cds_dna == part["sequence"].upper()
            assert modal._cds_meta is not None
            assert modal._cds_meta["origin"] == "parts"
            assert modal._cds_meta["name"] == "aeBlue"
            assert modal._cds_meta["start"] == 0
            assert modal._cds_meta["end"] == len(part["sequence"])
            assert modal._cds_meta["strand"] == 1
            # Translation should be valid (starts with M, ends with stop).
            protein = sc._mut_translate(modal._cds_dna)
            assert protein.startswith("M")
            app.exit()

    async def test_select_then_design_primers(self, monkeypatch):
        """End-to-end: pick part → enter mutation → design enables Save."""
        part = self._ok_part("partA")
        monkeypatch.setattr(sc, "_load_parts_bin", lambda: [part])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "parts"
            await pilot.pause(0.1)
            modal.query_one("#mut-parts").set_options(modal._build_parts_options())
            await pilot.pause(0.1)
            modal.query_one("#mut-parts").value = "0"
            await pilot.pause(0.2)
            # Protein is M + 30×A + (stop) → mutate position 5 (A) → V.
            modal.query_one("#mut-input").value = "A5V"
            modal.query_one("#btn-mut-design").action_press()
            await pilot.pause(0.2)
            assert modal._inner is not None
            assert modal._outer is not None
            assert modal.query_one("#btn-mut-save").disabled is False
            app.exit()

    async def test_source_switch_clears_parts_cds(self, monkeypatch):
        """Switching off 'parts' must clear the CDS so the user re-picks."""
        part = self._ok_part("partB")
        monkeypatch.setattr(sc, "_load_parts_bin", lambda: [part])
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal("", [], ""))
            await pilot.pause(0.3)
            modal = app.screen
            modal.query_one("#mut-source").value = "parts"
            await pilot.pause(0.1)
            modal.query_one("#mut-parts").set_options(modal._build_parts_options())
            await pilot.pause(0.1)
            modal.query_one("#mut-parts").value = "0"
            await pilot.pause(0.2)
            assert modal._cds_dna != ""
            modal.query_one("#mut-source").value = "prot"
            await pilot.pause(0.1)
            assert modal._cds_dna == ""
            assert modal._cds_meta is None
            app.exit()


class TestLibrarySourceWrapFeature:
    """Regression: a CompoundLocation CDS that spans the origin must be
    loaded correctly from the plasmid library, not flattened into a
    whole-plasmid range."""

    def test_compound_location_preserves_wrap(self, tmp_path, monkeypatch):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

        # 120 nt plasmid; CDS spans origin: nt 100..120 + 1..10 = 30 nt.
        # Build a known wrap CDS cleanly:
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


# ── AA click resolution under real CSS  (2026-05-07 regression guard) ────────

class TestMutagenizeClickAlignment:
    """End-to-end check that clicking an amino acid in the mutagenize
    preview lands on the codon the user actually clicked. Pre-fix
    ``on_click`` used ``event.screen_x - self.region.x`` which includes
    the widget's CSS ``border: solid`` + ``padding: 0 1`` overhead (4
    cols of horizontal chrome on the left edge); the resolved column
    was 2 cols right of the click target, so each click resolved to the
    codon ~one AA to the LEFT of where the user clicked. ``on_click``
    now uses ``event.x`` (already content-relative) so the click coord
    matches the rendered AA position exactly."""

    async def test_aa_click_lands_on_clicked_codon(self):
        # 33-bp CDS = 11 aa. The CDS feature filter requires >= 30 bp
        # and len % 3 == 0; 33 satisfies both.
        cds = "ATGGCCAGCAAATTCCATTGGGCAGAAGCCTAA"   # M A S K F H W A E A *
        assert len(cds) == 33 and len(cds) % 3 == 0
        feats = [{"type": "CDS", "label": "testCDS",
                  "start": 0, "end": len(cds), "strand": 1}]
        app = sc.PlasmidApp()
        async with app.run_test(size=(160, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.MutagenizeModal(cds, feats, "TEST"))
            await pilot.pause(0.3)
            modal = app.screen
            cds_select = modal.query_one("#mut-cds")
            cds_select.value = f"0:{len(cds)}:1"
            await pilot.pause(0.3)
            preview = modal.query_one("#mut-preview", sc._MutPreview)
            assert preview._dna_mode is True
            assert preview._protein.startswith("MASKFHW")
            # Drive synthetic Click events whose ``event.x`` matches
            # what real Textual delivers: ``screen_x - widget.region.x``,
            # i.e. measured from the widget's OUTER edge (gutter
            # included). The handler then subtracts ``content_offset``
            # to reach the content-area col. Reading the live
            # ``content_offset.x`` from the mounted widget keeps the
            # test honest under whatever border / padding the modal
            # CSS applies.
            from types import SimpleNamespace
            n = len(cds)
            num_w = len(str(n))
            pad = num_w + 2
            gutter_x = int(preview.content_offset.x)
            assert gutter_x > 0, (
                "test only meaningful when the widget actually has a "
                f"non-zero gutter; got content_offset.x={gutter_x}"
            )

            def _click_at_content_col(c: int) -> None:
                """Simulate a real Textual click at content col `c`,
                AA row of chunk 0 (= row 2). ``event.x`` includes the
                gutter offset because Textual delivers it that way."""
                preview._cursor_aa = -1
                evt = SimpleNamespace(
                    x=gutter_x + c, y=2,
                    screen_x=gutter_x + c, screen_y=2,
                    chain=1, button=1,
                )
                preview.on_click(evt)

            # AA letters land at content cols pad+1, pad+4, pad+7, …
            # corresponding to codons 0, 1, 2, …
            _click_at_content_col(pad + 1)
            await pilot.pause(0.05)
            assert preview._cursor_aa == 0, (
                f"clicked 'M' (codon 0); got cursor on "
                f"{preview._cursor_aa} "
                f"('{preview._protein[preview._cursor_aa] if preview._cursor_aa >= 0 else 'NONE'}')"
            )
            _click_at_content_col(pad + 4)
            await pilot.pause(0.05)
            assert preview._cursor_aa == 1   # 'A'
            _click_at_content_col(pad + 7)
            await pilot.pause(0.05)
            assert preview._cursor_aa == 2   # 'S'
            app.exit()


# ── TSV import parser ─────────────────────────────────────────────────────────

class TestParseCodonTsv:
    def test_codon_count_two_col(self):
        raw = sc._parse_codon_tsv("ATG\t127\nGCT\t120\n")
        assert raw["ATG"] == ("M", 127)
        assert raw["GCT"] == ("A", 120)

    def test_codon_aa_count_three_col(self):
        raw = sc._parse_codon_tsv("ATG\tM\t100\nTAA\t*\t5\n")
        assert raw["ATG"] == ("M", 100)
        assert raw["TAA"] == ("*", 5)

    def test_three_letter_aa(self):
        raw = sc._parse_codon_tsv("ATG Met 9\nTAA Stop 1\n")
        assert raw["ATG"] == ("M", 9)
        assert raw["TAA"] == ("*", 1)

    def test_skips_header_comment_blank(self):
        raw = sc._parse_codon_tsv(
            "# my table\ncodon\taa\tcount\n\nATG\tM\t12\n")
        assert set(raw) == {"ATG"}

    def test_u_folded_to_t(self):
        raw = sc._parse_codon_tsv("AUG 5\n")
        assert "ATG" in raw

    def test_fraction_only_scaled(self):
        raw = sc._parse_codon_tsv("GCT\t0.5\nATG\t0.25\n")
        assert raw["GCT"][1] == 500
        assert raw["ATG"][1] == 250

    def test_fraction_and_count_prefers_count(self):
        raw = sc._parse_codon_tsv("ATG\tM\t0.9\t127\n")
        assert raw["ATG"] == ("M", 127)

    def test_comma_delimited(self):
        raw = sc._parse_codon_tsv("ATG,M,12\nGCT,A,4\n")
        assert raw["ATG"] == ("M", 12)

    def test_aa_mismatch_raises(self):
        with pytest.raises(ValueError):
            sc._parse_codon_tsv("ATG\tA\t5\n")   # ATG is Met, not Ala

    def test_missing_count_raises(self):
        with pytest.raises(ValueError):
            sc._parse_codon_tsv("ATG\tM\n")

    def test_negative_count_raises(self):
        with pytest.raises(ValueError):
            sc._parse_codon_tsv("ATG\tM\t-5\n")

    def test_duplicate_codon_raises(self):
        with pytest.raises(ValueError):
            sc._parse_codon_tsv("ATG 5\nATG 6\n")

    def test_no_rows_raises(self):
        with pytest.raises(ValueError):
            sc._parse_codon_tsv("# just a comment\n\n")

    def test_non_str_raises(self):
        with pytest.raises(ValueError):
            sc._parse_codon_tsv(None)   # type: ignore[arg-type]

    def test_import_modal_and_handler_exist(self):
        # The paste-import feature is reachable: the modal class exists and
        # SpeciesPickerModal wires an Import-TSV handler.
        assert hasattr(sc, "CodonTsvImportModal")
        assert hasattr(sc.SpeciesPickerModal, "_import_tsv")


# ── Reachability — Settings entry point + Synthesis live dropdown refresh ────

class TestCodonTableReachability:
    """Codon-table collections are reachable from Settings (a new launcher)
    AND from Synthesis (a Manage button next to the dropdown), and the
    Synthesis dropdown rebuilds from the live registry so a freshly
    fetched / imported / deleted table shows up without reopening."""

    async def test_settings_codon_tables_opens_manager_and_sets_default(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            modal = sc.SettingsModal()
            app.push_screen(modal)
            await pilot.pause()
            await pilot.pause()
            from textual.widgets import Button
            # The launcher opens the codon-table manager (SpeciesPickerModal).
            modal.query_one("#set-codon-tables", Button).action_press()
            await pilot.pause()
            await pilot.pause()
            assert type(app.screen).__name__ == "SpeciesPickerModal"
            app.pop_screen()
            await pilot.pause()
            # "Use Selected" → callback persists the taxid as the launch
            # default that SynthesisScreen._init_codon_table honors on open.
            modal._codon_table_default_picked(
                {"name": "E. coli K12", "taxid": "83333"})
            assert sc._get_setting("active_codon_table", "") == "83333"
            # A taxid-less table can't be the default (no crash, unchanged).
            modal._codon_table_default_picked({"name": "custom", "taxid": ""})
            assert sc._get_setting("active_codon_table", "") == "83333"
            # Cancel (None) leaves the default untouched.
            modal._codon_table_default_picked(None)
            assert sc._get_setting("active_codon_table", "") == "83333"
            app.exit()

    async def test_synthesis_manage_button_and_live_refresh(self, monkeypatch):
        # Drive the codon registry through a controllable in-memory list so
        # the test never touches the real codon_tables.json.
        registry = {"v": [
            {"name": "E. coli K12", "taxid": "83333",
             "raw": {"GCT": ("A", 10)}, "source": "builtin"},
        ]}
        monkeypatch.setattr(
            sc, "_codon_tables_load",
            lambda: [dict(e) for e in registry["v"]])
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.SynthesisScreen())
            await pilot.pause(0.3)
            screen = app.screen
            from textual.widgets import Button
            # Manage button sits next to the codon dropdown.
            assert screen.query_one("#btn-syn-codon-manage", Button) is not None
            # Simulate the manager adding a user table, then the picker
            # closing with that table "Used".
            registry["v"].append(
                {"name": "My Lab Strain", "taxid": "999999",
                 "raw": {"GCT": ("A", 7)}, "source": "user"})
            screen._codon_manager_closed(
                {"name": "My Lab Strain", "taxid": "999999"})
            await pilot.pause(0.2)
            # The rebuilt dropdown now lists the new table …
            assert ("My Lab Strain|999999"
                    in [v for _l, v in screen._codon_table_options()])
            # … and it became the active selection (Select.Changed applied it).
            assert screen._codon_table_choice == "My Lab Strain|999999"
            app.exit()

    async def test_synthesis_refresh_falls_back_when_selection_deleted(
            self, monkeypatch):
        """If the currently-selected table is deleted in the manager, the
        refresh falls back to a still-present option instead of a dangling
        selection."""
        registry = {"v": [
            {"name": "E. coli K12", "taxid": "83333",
             "raw": {"GCT": ("A", 10)}, "source": "builtin"},
            {"name": "My Lab Strain", "taxid": "999999",
             "raw": {"GCT": ("A", 7)}, "source": "user"},
        ]}
        monkeypatch.setattr(
            sc, "_codon_tables_load",
            lambda: [dict(e) for e in registry["v"]])
        app = sc.PlasmidApp()
        async with app.run_test(size=(170, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.SynthesisScreen())
            await pilot.pause(0.3)
            screen = app.screen
            # Select the user strain, then "delete" it from the registry.
            screen._codon_manager_closed(
                {"name": "My Lab Strain", "taxid": "999999"})
            await pilot.pause(0.2)
            assert screen._codon_table_choice == "My Lab Strain|999999"
            registry["v"] = [registry["v"][0]]    # only K12 remains
            screen._codon_manager_closed(None)     # cancel after delete
            await pilot.pause(0.2)
            # Selection no longer dangles on the deleted table.
            assert screen._codon_table_choice == "E. coli K12|83333"
            app.exit()


# ── Codon optimization: bulletproofing (zero-tolerance correctness) ──────────
#
# The codon optimizer feeds real cloning/synthesis work, so these guard the
# invariants that "cannot fail": the optimized DNA encodes EXACTLY the input
# protein, contains only ACGT, is in-frame, carries no premature internal
# stop, and is deterministic. We assert them as a property sweep over many
# random proteins plus targeted edge cases. Translation here is done with the
# raw standard genetic code (`_CODON_GENETIC_CODE`) so the check is
# independent of the translator used inside the app.

_AA20 = "ACDEFGHIKLMNPQRSTVWY"
_STOP_CODONS = {"TAA", "TAG", "TGA"}


def _translate_std(dna: str) -> str:
    """3-frame-0 translation via the standard code; '?' for any non-codon.
    Does NOT stop at stops (so we can see the whole sequence including the
    appended stop run)."""
    gc = sc._CODON_GENETIC_CODE
    return "".join(gc.get(dna[i:i + 3], "?") for i in range(0, len(dna) - 2, 3))


def _tail_codons(dna: str, n: int) -> list:
    return [dna[i:i + 3] for i in range(len(dna) - 3 * n, len(dna), 3)]


class TestCodonAllocate:
    """The shared apportionment helper underpinning both residue and stop
    placement. Its length guarantee is what makes 'no rogue/empty codon'
    structurally true."""

    _ALA = [("GCT", 0.40), ("GCC", 0.35), ("GCA", 0.15), ("GCG", 0.10)]

    def test_length_is_exactly_n_for_every_n(self):
        for n in range(0, 64):
            out = sc._codon_allocate(self._ALA, n)
            assert len(out) == n
            assert all(c in {"GCT", "GCC", "GCA", "GCG"} for c in out)

    def test_zero_and_negative_give_empty(self):
        assert sc._codon_allocate(self._ALA, 0) == []
        assert sc._codon_allocate(self._ALA, -7) == []

    def test_single_codon_repeats(self):
        assert sc._codon_allocate([("ATG", 1.0)], 5) == ["ATG"] * 5

    def test_deterministic(self):
        assert sc._codon_allocate(self._ALA, 23) == sc._codon_allocate(self._ALA, 23)

    def test_dominant_codon_wins(self):
        from collections import Counter
        out = sc._codon_allocate([("GCT", 0.7), ("GCC", 0.2), ("GCA", 0.1)], 100)
        c = Counter(out)
        assert c["GCT"] == max(c.values())
        assert c["GCT"] >= 60

    def test_under_normalized_table_padded_to_n(self):
        # Fractions sum to 0.5 — apportionment falls short, must still pad.
        out = sc._codon_allocate([("GCT", 0.25), ("GCC", 0.25)], 8)
        assert len(out) == 8
        assert all(c in {"GCT", "GCC"} for c in out)

    def test_over_normalized_table_truncated_to_n(self):
        # Fractions sum to 2.0 — apportionment overshoots, must not exceed n.
        out = sc._codon_allocate([("GCT", 1.0), ("GCC", 1.0)], 5)
        assert len(out) == 5


class TestOptimizeProperties:
    """Core 'optimization actually works / no rogue bases' invariants."""

    def test_property_sweep_2000_proteins(self):
        import random
        rng = random.Random(20260530)
        K12 = sc._CODON_BUILTIN_K12
        for _ in range(2000):
            p = "".join(rng.choice(_AA20) for _ in range(rng.randint(0, 60)))
            dna = sc._codon_optimize(p, K12)
            assert all(b in "ACGT" for b in dna), (p, dna)     # no rogue base
            assert len(dna) == 3 * len(p) + 3                  # frame + length
            assert len(dna) % 3 == 0
            assert _translate_std(dna) == p + "*"              # exact protein
            assert dna[-3:] == "TAA"                           # single stop=TAA
            assert "*" not in _translate_std(dna[:-3])         # no internal stop
            for i, aa in enumerate(p):                         # every codon right
                assert K12[dna[3 * i:3 * i + 3]][0] == aa
            assert sc._codon_optimize(p, K12) == dna           # deterministic

    def test_all_twenty_homopolymers(self):
        for aa in _AA20:
            dna = sc._codon_optimize(aa * 17, sc._CODON_BUILTIN_K12)
            assert _translate_std(dna) == aa * 17 + "*"
            assert all(b in "ACGT" for b in dna)
            assert len(dna) == 17 * 3 + 3

    def test_lowercase_input_same_as_upper(self):
        up = sc._codon_optimize("MAKLEND", sc._CODON_BUILTIN_K12)
        lo = sc._codon_optimize("maklend", sc._CODON_BUILTIN_K12)
        assert lo == up
        assert sc._mut_translate(lo) == "MAKLEND"

    def test_empty_protein_is_lone_stop(self):
        assert sc._codon_optimize("", sc._CODON_BUILTIN_K12) == "TAA"

    def test_single_residue(self):
        assert sc._codon_optimize("M", sc._CODON_BUILTIN_K12) == "ATGTAA"

    def test_unknown_aa_raises(self):
        for bad in ("MAXA", "MZA", "MBK", "MUK", "MOK"):
            with pytest.raises(ValueError, match="No codons"):
                sc._codon_optimize(bad, sc._CODON_BUILTIN_K12)

    def test_optimize_then_fix_sites_preserves_protein(self):
        # Mirrors the Mutagenize worker: optimize → BsaI scrub. The protein
        # and the stop run must both survive the synonymous swaps.
        body = "MAKLEDGGRSTVWYFHIQNP" * 3
        dna = sc._codon_optimize(body, sc._CODON_BUILTIN_K12, stops=2)
        fixed, _fixes = sc._codon_fix_sites(
            dna, body, sc._CODON_BUILTIN_K12, {"BsaI": "GGTCTC"})
        assert len(fixed) == len(dna)
        assert all(b in "ACGT" for b in fixed)
        assert sc._mut_translate(fixed) == sc._mut_translate(dna) == body


class TestOptimizeStops:
    """Double / triple stop-codon handling (2026-05-30)."""

    def test_trailing_stops_honored_1_2_3(self):
        for k in (1, 2, 3):
            dna = sc._codon_optimize("MGK" + "*" * k, sc._CODON_BUILTIN_K12)
            assert len(dna) == 3 * 3 + 3 * k
            assert all(c in _STOP_CODONS for c in _tail_codons(dna, k))
            assert _translate_std(dna[:9]) == "MGK"          # body intact
            assert sc._mut_translate(dna) == "MGK"           # truncates at stop

    def test_stops_kwarg_counts(self):
        for k in (0, 1, 2, 3):
            dna = sc._codon_optimize("MGK", sc._CODON_BUILTIN_K12, stops=k)
            assert len(dna) == 9 + 3 * k
            assert all(c in _STOP_CODONS for c in _tail_codons(dna, k))
            assert sc._mut_translate(dna) == "MGK"

    def test_trailing_run_overrides_kwarg(self):
        dna = sc._codon_optimize("MGK***", sc._CODON_BUILTIN_K12, stops=1)
        assert len(dna) == 9 + 9                              # 3 stops, not 1
        assert all(c in _STOP_CODONS for c in _tail_codons(dna, 3))

    def test_mid_sequence_stop_raises(self):
        for bad in ("M*GK", "*MGK", "MG*K", "M**GK", "MG*K*"):
            with pytest.raises(ValueError, match="only allowed at the end"):
                sc._codon_optimize(bad, sc._CODON_BUILTIN_K12)

    def test_single_stop_is_always_TAA(self):
        assert sc._codon_optimize("MGK", sc._CODON_BUILTIN_K12).endswith("TAA")
        assert sc._codon_optimize("MGK*", sc._CODON_BUILTIN_K12).endswith("TAA")
        assert sc._codon_optimize(
            "MGK", sc._CODON_BUILTIN_K12, stops=1).endswith("TAA")

    def test_stops_zero_emits_no_stop(self):
        dna = sc._codon_optimize("MGK", sc._CODON_BUILTIN_K12, stops=0)
        assert len(dna) == 9
        assert sc._mut_translate(dna) == "MGK"

    def test_negative_stops_clamped_to_zero(self):
        assert (sc._codon_optimize("MGK", sc._CODON_BUILTIN_K12, stops=-5)
                == sc._codon_optimize("MGK", sc._CODON_BUILTIN_K12, stops=0))

    def test_multi_stops_are_frequency_matched(self):
        from collections import Counter
        # E. coli K12: TAA is the most-used stop → it leads the run.
        tail = _tail_codons(
            sc._codon_optimize("MGK", sc._CODON_BUILTIN_K12, stops=3), 3)
        assert all(c in _STOP_CODONS for c in tail)
        assert Counter(tail).most_common(1)[0][0] == "TAA"
        # A synthetic TGA-dominant organism → the run is all TGA.
        tga = dict(sc._CODON_BUILTIN_K12)
        tga["TGA"] = ("*", 900)
        tga["TAA"] = ("*", 50)
        tga["TAG"] = ("*", 50)
        assert _tail_codons(
            sc._codon_optimize("MGK", tga, stops=3), 3) == ["TGA", "TGA", "TGA"]

    def test_no_internal_stop_with_multi_stops(self):
        dna = sc._codon_optimize(
            "M" + "G" * 30 + "K", sc._CODON_BUILTIN_K12, stops=3)
        assert "*" not in _translate_std(dna[:-9])           # body has no stop

    def test_table_without_stops_falls_back_to_taa(self):
        # A custom table that declares no stop codons must still terminate.
        nostop = {c: v for c, v in sc._CODON_BUILTIN_K12.items()
                  if v[0] != "*"}
        dna = sc._codon_optimize("MGK", nostop, stops=2)
        assert _tail_codons(dna, 2) == ["TAA", "TAA"]
        assert sc._mut_translate(dna) == "MGK"


class TestCodonTableConsistency:
    """The optimizer trusts the loaded table's codon→AA labels, so every
    loader MUST pin those labels to the standard genetic code — otherwise an
    optimized sequence could encode the wrong residue. These cross-checks
    guard that 'no mismatch' contract for the builtin, Kazusa, and TSV paths.
    """

    def test_builtin_k12_is_complete_and_standard(self):
        K12 = sc._CODON_BUILTIN_K12
        assert len(K12) == 64
        for codon, (aa, _ct) in K12.items():
            assert aa == sc._CODON_GENETIC_CODE[codon]
        aamap, _ = sc._codon_build_aa_map(K12)
        assert all(aamap.get(a) for a in _AA20)              # all 20 present

    def test_kazusa_table_labels_match_genetic_code(self):
        raw = sc._codon_parse_kazusa_html(_FAKE_KAZUSA_HTML)
        assert raw is not None and len(raw) == 64
        # The parser assigns AAs from _CODON_GENETIC_CODE, never the file —
        # so no codon can be mislabeled.
        for codon, (aa, _ct) in raw.items():
            assert aa == sc._CODON_GENETIC_CODE[codon]

    def test_kazusa_table_optimizes_without_mismatch(self):
        import random
        raw = sc._codon_parse_kazusa_html(_FAKE_KAZUSA_HTML)
        rng = random.Random(11)
        for _ in range(300):
            p = "".join(rng.choice(_AA20) for _ in range(rng.randint(1, 40)))
            dna = sc._codon_optimize(p, raw)
            assert _translate_std(dna) == p + "*"
            assert all(b in "ACGT" for b in dna)

    def test_tsv_table_labels_match_and_optimize(self):
        # A complete TSV built straight from the genetic code.
        text = "\n".join(f"{c}\t{a}\t10"
                         for c, a in sc._CODON_GENETIC_CODE.items())
        raw = sc._parse_codon_tsv(text)
        assert len(raw) == 64
        for codon, (aa, _ct) in raw.items():
            assert aa == sc._CODON_GENETIC_CODE[codon]
        dna = sc._codon_optimize("MKLAEVNDPQRSTWYF", raw)
        assert sc._mut_translate(dna) == "MKLAEVNDPQRSTWYF"

    def test_tsv_codon_aa_mismatch_is_rejected(self):
        # ATG encodes M; a file claiming L must be refused, not silently
        # loaded (which would later mistranslate the optimized output).
        with pytest.raises(ValueError, match="encodes"):
            sc._parse_codon_tsv("ATG\tL\t5\n")


class TestKazusaIncompleteDownload:
    """A Kazusa fetch cut short by a dropped connection / interrupted Wi-Fi
    must be REJECTED (returns None), never silently accepted as a partial
    table that would later make the optimizer raise mid-design. The parser
    requires all 64 codons."""

    @staticmethod
    def _first_n_codon_rows(n: int) -> str:
        """The fake table truncated to its first `n` codon rows."""
        import re
        kept, out = 0, []
        for ln in _FAKE_KAZUSA_HTML.splitlines():
            is_codon = bool(re.search(r"\b[ACGTU]{3}\b\s+\d", ln)) \
                and "Codon" not in ln
            if is_codon:
                if kept >= n:
                    continue
                kept += 1
            out.append(ln)
        return "\n".join(out)

    def test_complete_table_parses_to_64(self):
        raw = sc._codon_parse_kazusa_html(_FAKE_KAZUSA_HTML)
        assert raw is not None and len(raw) == 64

    def test_first_40_codons_rejected(self):
        assert sc._codon_parse_kazusa_html(self._first_n_codon_rows(40)) is None

    def test_exactly_63_codons_rejected(self):
        # One missing codon (Gly/GGG row removed) → still rejected.
        partial = "\n".join(ln for ln in _FAKE_KAZUSA_HTML.splitlines()
                            if "GGG" not in ln)
        assert sc._codon_parse_kazusa_html(partial) is None

    def test_midflight_byte_cutoff_rejected(self):
        # Raw string sliced mid-stream (possibly mid-line), as a dropped
        # connection would leave it.
        for cut in (800, 1400, len(_FAKE_KAZUSA_HTML) // 2):
            assert sc._codon_parse_kazusa_html(_FAKE_KAZUSA_HTML[:cut]) is None

    def test_empty_and_garbage_rejected(self):
        assert sc._codon_parse_kazusa_html("") is None
        assert sc._codon_parse_kazusa_html("<html><body></body></html>") is None

    def test_fetch_complete_returns_full_table(self, monkeypatch):
        import urllib.request

        def fake(req, timeout=None):
            return _FakeResponse(_FAKE_KAZUSA_HTML.encode())

        monkeypatch.setattr(urllib.request, "urlopen", fake)
        raw, msg = sc._codon_fetch_kazusa("83333")
        assert raw is not None and len(raw) == 64
        assert "fetched" in msg.lower()

    def test_fetch_truncated_returns_parse_error(self, monkeypatch):
        import urllib.request

        def fake(req, timeout=None):
            return _FakeResponse(_FAKE_KAZUSA_HTML[:1400].encode())

        monkeypatch.setattr(urllib.request, "urlopen", fake)
        raw, msg = sc._codon_fetch_kazusa("83333")
        assert raw is None
        assert "parse" in msg.lower()

    def test_fetch_not_found_message(self, monkeypatch):
        import urllib.request

        def fake(req, timeout=None):
            return _FakeResponse(b"<html><body>Species not found</body></html>")

        monkeypatch.setattr(urllib.request, "urlopen", fake)
        raw, msg = sc._codon_fetch_kazusa("999999999")
        assert raw is None
        assert "not found" in msg.lower()


class TestOptimizeProteinEndpoint:
    """`optimize-protein` agent endpoint. The handler is pure (ignores
    `app`, reads the codon registry), so we call it directly. The `stops`
    param is OPTIONAL and backward compatible — omitting it keeps the
    historical single-TAA behaviour."""

    def test_default_is_single_taa(self):
        r = sc._h_optimize_protein(None, {"protein": "MGK"})
        assert r["ok"] is True
        assert r["dna"].endswith("TAA")
        assert r["length"] == 12          # 3 residues + 1 stop

    def test_stops_param_appends_that_many(self):
        r = sc._h_optimize_protein(None, {"protein": "MGK", "stops": 3})
        assert r["length"] == 9 + 9
        assert sc._mut_translate(r["dna"]) == "MGK"

    def test_trailing_stars_honored(self):
        r = sc._h_optimize_protein(None, {"protein": "MGK**"})
        assert r["length"] == 9 + 6       # two stops

    def test_out_of_range_stops_400(self):
        _body, code = sc._h_optimize_protein(None, {"protein": "MGK", "stops": 9})
        assert code == 400

    def test_non_integer_stops_400(self):
        _body, code = sc._h_optimize_protein(
            None, {"protein": "MGK", "stops": "lots"})
        assert code == 400

    def test_mid_sequence_stop_400(self):
        _body, code = sc._h_optimize_protein(None, {"protein": "M*GK"})
        assert code == 400

    def test_missing_protein_400(self):
        _body, code = sc._h_optimize_protein(None, {})
        assert code == 400


class TestNcbiTaxonPickerModalStyle:
    """The NCBI taxon picker is a centered modal dialog (like the species
    picker), not a full-screen panel. Empty initial query → no network."""

    async def test_renders_as_centered_dialog(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.NcbiTaxonPickerModal(""))
            await pilot.pause()
            box = app.screen.query_one("#ncbi-box")
            # Centered + bounded box, not the full 171×43 screen.
            assert box.region.width <= 92        # ~90 + border, not 171
            assert box.region.x > 0              # inset from the left edge
            assert box.region.height < 43        # not full height
            # The list it wraps is present and styled as a bordered list.
            assert app.screen.query_one("#ncbi-list") is not None


class TestPickerDataTables:
    """The codon-usage list and the NCBI results list are true DataTables
    (native zebra striping + columns), not ListViews — and selecting a row
    drives the Use/Fetch button + commit correctly."""

    async def test_codon_list_is_datatable_and_use_commits(self):
        from textual.widgets import DataTable, Button
        app = sc.PlasmidApp()
        async with app.run_test(size=(140, 50)) as pilot:
            await pilot.pause()
            picked = {}
            await app.push_screen(
                sc.SpeciesPickerModal(),
                callback=lambda e: picked.update(e or {}))
            await pilot.pause(0.3)
            modal = app.screen
            dt = modal.query_one("#sp-list", DataTable)
            assert dt.zebra_stripes is True
            assert len(dt.columns) == 3          # Species / Taxid / Source
            assert len(dt.rows) >= 1             # builtin K12 always present
            # A row is highlighted → Use is live; Use dismisses with it.
            dt.move_cursor(row=0)
            await pilot.pause(0.1)
            assert modal.query_one("#btn-sp-use", Button).disabled is False
            expect = dict(modal._entries[0])
            modal.query_one("#btn-sp-use", Button).action_press()
            await pilot.pause(0.2)
            assert picked.get("taxid") == expect["taxid"]

    async def test_ncbi_list_is_datatable_and_populates(self):
        from textual.widgets import DataTable, Button
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 43)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.NcbiTaxonPickerModal(""))
            await pilot.pause(0.2)
            modal = app.screen
            dt = modal.query_one("#ncbi-list", DataTable)
            assert dt.zebra_stripes is True
            assert len(dt.columns) == 2          # Species / Taxid
            # Feed canned hits straight to the result handler (no network).
            modal._search_done(
                [{"name": "Homo sapiens", "taxid": "9606"},
                 {"name": "Mus musculus", "taxid": "10090"}], 2, "2 hits")
            await pilot.pause(0.1)
            assert len(dt.rows) == 2
            assert modal.query_one("#btn-ncbi-use", Button).disabled is False
            assert modal._hits[0]["taxid"] == "9606"


class TestProteinOptimizeToDna:
    """Synthesis Protein tab → 'Optimize → DNA': codon-optimizes the protein
    and hands the CDS to the DNA tab as a fresh, editable fragment."""

    async def test_optimizes_and_loads_into_dna_tab(self):
        from textual.widgets import Select, TabbedContent
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.SynthesisScreen())
            await pilot.pause(0.3)
            screen = app.screen
            pe = screen.query_one("#syn-protein-editor", sc.ProteinEditor)
            pe.load("MAEVKLAGHIKQRSTVWY")
            await pilot.pause(0.1)
            screen.query_one("#syn-codon-stops", Select).value = "2"
            await pilot.pause(0.1)
            screen.query_one("#btn-syn-optimize-dna").action_press()
            await pilot.pause(0.6)            # threaded worker + apply
            # The DNA tab is now active and holds the optimized CDS.
            assert screen.query_one(
                "#syn-tabs", TabbedContent).active == "syn-tab-dna"
            ed = screen.query_one("#syn-editor", sc.SynthesisEditor)
            seq, feats = ed.get_state()
            assert sc._mut_translate(seq) == "MAEVKLAGHIKQRSTVWY"
            assert len(seq) == 18 * 3 + 2 * 3       # body + 2 stops
            assert all(b in "ACGT" for b in seq)
            assert screen._dirty is True            # fresh unsaved fragment
            assert screen._loaded_id is None
            assert any(f.get("type") == "CDS" for f in feats)
            app.exit()

    async def test_trailing_stops_override_selector(self):
        from textual.widgets import Select
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.SynthesisScreen())
            await pilot.pause(0.3)
            screen = app.screen
            screen.query_one("#syn-protein-editor", sc.ProteinEditor).load("MGK***")
            await pilot.pause(0.1)
            screen.query_one("#syn-codon-stops", Select).value = "1"
            await pilot.pause(0.1)
            screen.query_one("#btn-syn-optimize-dna").action_press()
            await pilot.pause(0.6)
            ed = screen.query_one("#syn-editor", sc.SynthesisEditor)
            seq = ed.get_state()[0]
            assert len(seq) == 3 * 3 + 3 * 3        # MGK + 3 trailing stops
            assert sc._mut_translate(seq) == "MGK"
            app.exit()

    async def test_prompts_before_clobbering_dirty_dna_tab(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.SynthesisScreen())
            await pilot.pause(0.3)
            screen = app.screen
            screen._dirty = True                    # DNA tab has unsaved edits
            screen.query_one("#syn-protein-editor", sc.ProteinEditor).load("MAEVK")
            await pilot.pause(0.1)
            screen.query_one("#btn-syn-optimize-dna").action_press()
            await pilot.pause(0.3)
            assert type(app.screen).__name__ == "SynthesisReplaceDnaConfirmModal"
            app.screen.query_one("#btn-srd-cancel").action_press()
            await pilot.pause(0.2)
            # Cancel returns to Synthesis with the DNA tab untouched.
            assert type(app.screen).__name__ == "SynthesisScreen"
            app.exit()

    async def test_staleguard_apply_bails_when_unmounted(self):
        """If the optimize worker's result lands after the Synthesis screen
        was torn down, the apply must no-op (is_mounted guard) rather than
        write into dead widgets. Force the guard deterministically."""
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.SynthesisScreen())
            await pilot.pause(0.3)
            screen = app.screen
            ed = screen.query_one("#syn-editor", sc.SynthesisEditor)
            before = ed.get_state()[0]
            cls = type(screen)
            cls.is_mounted = property(lambda self: False)   # simulate teardown
            try:
                screen._apply_optimize_to_dna(
                    "MGK", "ATGGGTAAATAA", [], sc._CODON_BUILTIN_K12, "x", 1)
                # Guard held — the DNA editor was left untouched.
                assert ed.get_state()[0] == before
            finally:
                del cls.is_mounted                           # restore Widget's
            app.exit()


class TestFixSitesHardening:
    """The cut-site remover is IUPAC-aware (degenerate sites match real DNA)
    and robust against odd input (2026-05-30 hardening)."""

    def test_hit_set_matches_degenerate_overlapping(self):
        # GGWCC = GG[AT]CC → matches GGACC@3 and GGTCC@10 (overlap-safe scan).
        hits = sc._forbidden_hit_set("AAAGGACCGGGGTCCAA", ("GGWCC",))
        assert sorted(p for _s, p in hits) == [3, 10]

    def test_hit_set_exact_site_unchanged(self):
        hits = sc._forbidden_hit_set("GGTCTCAAGGTCTC", ("GGTCTC",))
        assert sorted(p for _s, p in hits) == [0, 8]

    def test_property_random_forbidden_sets_preserve_protein(self):
        import random
        rng = random.Random(99)
        K12 = sc._CODON_BUILTIN_K12
        enz = sc._all_enzymes()
        names = list(enz)
        for _ in range(300):
            p = "".join(rng.choice(_AA20) for _ in range(rng.randint(20, 70)))
            cds = sc._codon_optimize(p, K12)
            sites = {n: enz[n][0] for n in rng.sample(names, rng.randint(1, 8))}
            fixed, _f = sc._codon_fix_sites(cds, p, K12, sites)
            assert sc._mut_translate(fixed) == sc._mut_translate(cds)   # synonymous
            assert all(b in "ACGT" for b in fixed)                      # no rogue
            assert len(fixed) == len(cds)                               # 3→3 swaps
            pats = tuple({s.upper() for s in sites.values()})
            # Never INTRODUCES a forbidden site (degenerate or exact).
            assert not (sc._forbidden_hit_set(fixed, pats)
                        - sc._forbidden_hit_set(cds, pats))

    def test_invalid_site_skipped_not_fatal(self):
        out, fixes = sc._codon_fix_sites(
            "ATGGGTAAA", "MG", sc._CODON_BUILTIN_K12, {"Weird": "GGXZ!"})
        assert out == "ATGGGTAAA" and fixes == []

    def test_empty_sites_is_noop(self):
        out, fixes = sc._codon_fix_sites(
            "ATGGGTAAA", "MG", sc._CODON_BUILTIN_K12, {})
        assert out == "ATGGGTAAA" and fixes == []

    def test_idempotent(self):
        body = "MAEVKLAGGRSTWND" * 3
        cds = sc._codon_optimize(body, sc._CODON_BUILTIN_K12)
        sites = {"BsaI": "GGTCTC", "EcoRI": "GAATTC"}
        f1, _ = sc._codon_fix_sites(cds, body, sc._CODON_BUILTIN_K12, sites)
        f2, fx2 = sc._codon_fix_sites(f1, body, sc._CODON_BUILTIN_K12, sites)
        assert f2 == f1 and fx2 == []


class TestForbiddenSitesSetting:
    """`_codon_forbidden_enzymes` setting → `_codon_fix_sites` site map."""

    @staticmethod
    def _patch(monkeypatch, names):
        monkeypatch.setattr(
            sc, "_get_setting",
            lambda k, d=None: names if k == "codon_forbidden_enzymes" else d)

    def test_resolves_names_to_sites(self, monkeypatch):
        self._patch(monkeypatch, ["BsaI", "EcoRI"])
        sites = sc._codon_forbidden_sites()
        assert sites.get("BsaI") == "GGTCTC"
        assert sites.get("EcoRI") == "GAATTC"

    def test_skips_unknown_names(self, monkeypatch):
        self._patch(monkeypatch, ["BsaI", "NotARealEnzyme"])
        sites = sc._codon_forbidden_sites()
        assert "BsaI" in sites and "NotARealEnzyme" not in sites

    def test_empty_list_means_no_scrub(self, monkeypatch):
        self._patch(monkeypatch, [])
        assert sc._codon_forbidden_sites() == {}

    def test_label_reflects_count(self, monkeypatch):
        self._patch(monkeypatch, ["BsaI", "EcoRI", "NdeI"])
        assert sc._forbidden_sites_label() == "Avoid sites (3)"


class TestForbiddenSitesModalAndButtons:
    """The picker modal + its 'Avoid sites' buttons in Mutato / Synthesis."""

    async def test_modal_toggles_and_returns_selection(self):
        from textual.widgets import DataTable, Button
        import types
        app = sc.PlasmidApp()
        async with app.run_test(size=(120, 50)) as pilot:
            await pilot.pause()
            result = {}
            await app.push_screen(
                sc.ForbiddenSitesModal(["BsaI"]),
                callback=lambda r: result.update(picked=r))
            await pilot.pause(0.2)
            modal = app.screen
            assert len(modal.query_one("#fsm-table", DataTable).columns) == 3
            assert "BsaI" in modal._selected            # pre-selected
            assert "EcoRI" in modal._rows               # common set listed
            # Toggle EcoRI on via the row handler.
            i = modal._rows.index("EcoRI")
            modal._toggle(types.SimpleNamespace(cursor_row=i))
            await pilot.pause(0.1)
            assert "EcoRI" in modal._selected
            modal.query_one("#btn-fsm-done", Button).action_press()
            await pilot.pause(0.2)
            assert set(result["picked"]) == {"BsaI", "EcoRI"}

    async def test_synthesis_has_avoid_sites_button(self):
        app = sc.PlasmidApp()
        async with app.run_test(size=(171, 50)) as pilot:
            await pilot.pause()
            await app.push_screen(sc.SynthesisScreen())
            await pilot.pause(0.3)
            btn = app.screen.query_one("#btn-syn-forbidden")
            assert str(btn.label).startswith("Avoid sites")
            app.exit()
