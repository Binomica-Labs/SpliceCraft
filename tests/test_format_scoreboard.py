"""The reader scoreboard (`scripts/format_scoreboard.py`) must not rot.

Its ground truth is the census in `tests/_format_oracles.py`: a reader that
disagrees with it is scored as having lost content. So the census is pinned
here against hand-counted files, and the harness's own judgements — what is a
loss, what is the same record, what may be saved about a PRIVATE file —
against constructed cases. No network: the corpus is not needed.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import splicecraft as sc
from tests._format_oracles import dna_census, embl_census, genbank_census

_ROOT = Path(__file__).resolve().parent.parent
_SCRIPT = _ROOT / "scripts" / "format_scoreboard.py"
if not _SCRIPT.is_file():          # the sdist ships tests/ but not scripts/
    pytest.skip("scripts/format_scoreboard.py not present (source distribution)",
                allow_module_level=True)
_spec = importlib.util.spec_from_file_location("format_scoreboard", _SCRIPT)
assert _spec is not None and _spec.loader is not None
board = importlib.util.module_from_spec(_spec)
sys.modules["format_scoreboard"] = board     # dataclasses resolve their module by name
_spec.loader.exec_module(board)

_CASES = {name: (kind, text, enc) for name, kind, text, enc in board.synthetic_cases()}


def _census_of(name: str) -> "list[dict]":
    kind, text, enc = _CASES[name]
    return genbank_census(text.encode(enc).decode(enc))


# ── the census ────────────────────────────────────────────────────────────────

class TestGenbankCensus:
    def test_counts_the_canonical_record_by_hand(self):
        (c,) = _census_of("canonical")
        assert (c["declared_len"], c["seq_len"], c["topology"]) == (1200, 1200, "circular")
        assert c["keys"] == ["source", "promoter", "CDS", "misc_feature",
                             "misc_feature", "CDS", "primer_bind", "terminator"]
        # source 2 + promoter 1 + CDS 3 (its note runs over two lines: one
        # qualifier) + 5 more labels.
        assert c["n_qualifiers"] == 11

    @pytest.mark.parametrize("name", [
        "crlf", "cr_only", "utf8_bom", "trailing_whitespace", "vendor_column_42",
        "features_header_spacing", "tab_indented_features", "over_indented_key",
        "slash_line_in_quoted_value", "blank_line_in_features",
    ])
    def test_layouts_that_move_columns_count_the_same(self, name):
        """What a file DECLARES does not depend on its spacing — which is
        exactly what the readers are being tested against."""
        (base,) = _census_of("canonical")
        (c,) = _census_of(name)
        for k in ("seq_len", "n_features", "keys", "n_qualifiers"):
            assert c[k] == base[k], (name, k)

    def test_an_unclosed_quote_does_not_swallow_the_next_feature(self):
        (c,) = _census_of("unterminated_quote")
        assert c["n_features"] == 8

    def test_every_base_counts_without_position_numbers(self):
        (c,) = _census_of("origin_without_numbers")
        assert c["seq_len"] == 1200

    def test_a_contig_record_has_no_known_length(self):
        (c,) = _census_of("contig_only")
        assert c["contig"] and c["seq_len"] is None

    def test_records_split_at_their_terminators(self):
        assert [c["seq_len"] for c in _census_of("multi_record")] == [1200, 1200]

    def test_no_locus_line_means_no_record(self):
        assert _census_of("no_locus_line") == []

    def test_a_junk_position_number_is_not_bases(self):
        """Calibrated on fuzz: a stray character in the number column is
        still the number column, which every reader skips."""
        from tests._format_oracles import _sequence_letters
        assert _sequence_letters("       61 acgtacgtac gtacg") == 15
        assert _sequence_letters("      X61 acgtacgtac gtacg") == 15
        assert _sequence_letters("acgtacgtac gtacg") == 15

    def test_a_qualifier_with_no_feature_above_it_belongs_to_nothing(self):
        text = _CASES["canonical"][1].replace("     source          1..1200\n", "")
        (c,) = genbank_census(text)
        assert c["n_features"] == 7 and c["n_qualifiers"] == 9

    def test_sequence_lines_are_never_feature_keys(self):
        """A record that lost its ORIGIN header must not have its sequence
        lines counted as features (the position number is a token)."""
        text = _CASES["canonical"][1].replace("\nORIGIN\n", "\n")
        (c,) = genbank_census(text)
        assert c["n_features"] == 8


class TestEmblCensus:
    _EMBL = (
        "ID   X1; SV 1; circular; genomic DNA; STD; SYN; 30 BP.\nXX\n"
        "FT   source          1..30\n"
        'FT                   /organism="synthetic DNA construct"\n'
        "FT   CDS             1..30\n"
        'FT                   /note="two line\n'
        'FT                   note"\n'
        "FT                   /pseudo\nXX\n"
        "SQ   Sequence 30 BP; 8 A; 7 C; 8 G; 7 T; 0 other;\n"
        "     acgtacgtac gtacgtacgt acgtacgtac                                    30\n//\n")

    def test_counts_by_hand(self):
        (c,) = embl_census(self._EMBL)
        assert (c["declared_len"], c["seq_len"], c["topology"]) == (30, 30, "circular")
        assert c["keys"] == ["source", "CDS"]
        assert c["qualifier_names"] == ["organism", "note", "pseudo"]

    def test_line_codes_are_read_case_blind(self):
        (c,) = embl_census(self._EMBL.replace("FT   CDS", "ft   CDS"))
        assert c["keys"] == ["source", "CDS"]

    def test_a_record_shipped_without_bases_has_no_known_length(self):
        text = self._EMBL.split("     acgt")[0] + "//\n"
        (c,) = embl_census(text)
        assert c["seq_len"] is None


def _packet(kind: int, payload: bytes) -> bytes:
    return bytes([kind]) + len(payload).to_bytes(4, "big") + payload


class TestDnaCensus:
    def test_counts_a_packet_stream_by_hand(self):
        feats = (b'<Features><Feature name="a&amp;b" type="CDS"><Segment range="1-5" '
                 b'color="#ff0000"/></Feature><Feature name="p" type="primer_bind">'
                 b'<Segment range="6-9"/></Feature></Features>')
        primers = (b'<Primers><Primer name="x"><BindingSite location="1-5"/>'
                   b'<BindingSite location="7-9"/></Primer></Primers>')
        data = (_packet(0x09, b"cookie") + _packet(0x00, b"\x01" + b"ACGT" * 5)
                + _packet(0x0A, feats) + _packet(0x05, primers) + _packet(0x07, b"xz"))
        c = dna_census(data)
        assert (c["seq_len"], c["topology"]) == (20, "circular")
        assert c["feature_names"] == ["a&b", "p"]
        assert c["feature_types"] == ["CDS", "primer_bind"]
        assert (c["n_colored"], c["n_primers"], c["n_binding_sites"]) == (1, 1, 2)
        assert c["has_history"] and not c["truncated"]

    def test_a_length_past_the_end_is_reported_not_trusted(self):
        data = _packet(0x00, b"\x00ACGT")[:-2]
        assert dna_census(data)["truncated"]

    @pytest.mark.parametrize("path", sorted((_ROOT / "tests").glob("*.dna"))[:3],
                             ids=lambda p: p.stem)
    def test_agrees_with_splicecraft_on_the_repo_fixtures(self, path):
        """Calibration on real files: the census and SpliceCraft's reader
        count the same bases and the same non-primer features."""
        c = dna_census(path.read_bytes())
        rec = sc.load_genbank(str(path))
        assert c["seq_len"] == len(rec.seq)
        n = sum(1 for t in c["feature_types"] if t != "primer_bind")
        assert n == sum(1 for f in rec.features if f.type not in ("primer_bind", "source"))


# ── the harness's judgements ─────────────────────────────────────────────────

def _rec(*feats, topology="circular", n=30) -> SeqRecord:
    r = SeqRecord(Seq("A" * n), id="x")
    r.annotations["topology"] = topology
    r.features = [SeqFeature(SimpleLocation(a, b, strand=s), type=t,
                             qualifiers=dict(q)) for t, a, b, s, q in feats]
    return r


def _cen(keys, quals, *, seq_len=30, topology="circular") -> dict:
    return {"declared_len": seq_len, "seq_len": seq_len, "topology": topology,
            "contig": False, "n_features": len(keys), "keys": list(keys),
            "n_qualifiers": len(quals), "qualifier_names": list(quals)}


class TestJudgements:
    F = ("misc_feature", 0, 10, 1, {"label": ["a"]})

    def test_a_faithful_read_has_no_issues(self):
        assert board.census_issues([_rec(self.F)], [_cen(["misc_feature"], ["label"])], "genbank") == []

    def test_each_kind_of_loss_is_named(self):
        cen = _cen(["misc_feature", "CDS"], ["label"])
        assert board.census_issues([_rec(self.F)], [cen], "genbank") == ["feature_count"]
        cen = _cen(["CDS"], ["label"])
        assert board.census_issues([_rec(self.F)], [cen], "genbank") == ["feature_keys"]
        cen = _cen(["misc_feature"], ["label", "note"])
        assert board.census_issues([_rec(self.F)], [cen], "genbank") == ["qualifiers"]
        cen = _cen(["misc_feature"], ["label"], seq_len=31)
        assert board.census_issues([_rec(self.F)], [cen], "genbank") == ["bases"]
        cen = _cen(["misc_feature"], ["label"])
        assert board.census_issues([_rec(self.F, topology=None)], [cen], "genbank") == ["topology"]
        assert board.census_issues([], [cen], "genbank") == ["record_count"]

    @pytest.mark.parametrize("key,garbage", [
        ("misc_feature", False), ("5'UTR", False), ("-10_signal", False),
        ("isc_feature    2925", True), ("/organism=\"x", True), ("6", True),
        ("terminator    =", True)])
    def test_a_feature_key_no_writer_could_mean_is_garbage(self, key, garbage):
        """Shapes seen in the fuzz run: Biopython slices the key out of fixed
        columns, so a shifted line becomes a feature called `/organism=…`."""
        assert board.garbage_keys([_rec((key, 0, 10, 1, {}))]) == int(garbage)

    def test_a_repaired_key_or_a_decoded_marker_is_not_a_loss(self):
        cen = _cen(["misc_feature_"], ["label", "SpliceCraft_strand"])
        assert board.census_issues([_rec(self.F)], [cen], "genbank") == []

    def test_no_arrow_has_two_spellings(self):
        a = _rec(("misc_feature", 0, 10, None, {}))
        b = _rec(("misc_feature", 0, 10, 0, {}))
        assert board._first_diff(board.signature([a]), board.signature([b])) == ""

    def test_a_moved_feature_is_a_location_change(self):
        a = _rec(("misc_feature", 0, 10, 1, {}))
        b = _rec(("misc_feature", 1, 10, 1, {}))
        assert board._first_diff(board.signature([a]), board.signature([b])) == "location"

    def test_ascii_folding_is_reported_apart_from_a_real_change(self):
        def note(v):
            return _rec(("misc_feature", 0, 10, 1, {"note": [v]}))
        sig = board.signature
        assert board._first_diff(sig([note("5 µM")]), sig([note("5 uM")])) == "non_ascii_transliterated"
        assert board._first_diff(sig([note("abc")]), sig([note("abd")])) == "qualifiers"
        # A qualifier ADDED or DROPPED is never transliteration, whatever its
        # text — `all()` over nothing used to say it was.
        both = _rec(("misc_feature", 0, 10, 1, {"note": ["5 µM"], "color": ["#fff"]}))
        assert board._first_diff(sig([note("5 µM")]), sig([both])) == "qualifiers"
        assert board._first_diff(sig([both]), sig([note("5 µM")])) == "qualifiers"

    def test_untyped_values_compare_as_text(self):
        a = _rec(("CDS", 0, 9, 1, {"codon_start": [1]}))
        b = _rec(("CDS", 0, 9, 1, {"codon_start": ["1"]}))
        assert board._first_diff(board.signature([a]), board.signature([b])) == ""

    def test_an_undefined_sequence_does_not_break_the_signature(self):
        r = SeqRecord(Seq(None, length=10), id="contig")
        assert board.signature([r])[0][0] == "<undefined 10>"

    def test_synthetic_grades(self):
        def row(label, kind, outcome, issues=(), census_records=1):
            x = {"outcome": outcome}
            if issues:
                x["issues"] = list(issues)
            return {"label": label, "kind": kind, "census_records": census_records,
                    "by_reader": {"r": x}}
        rows = [row("a", "valid", "accepted"), row("b", "valid", "refused"),
                row("c", "valid", "accepted", ["bases"]), row("d", "malformed", "refused"),
                row("e", "malformed", "crashed"), row("f", "malformed", "empty", census_records=0)]
        assert board._synthetic_grades(rows, "r") == {
            "a": "pass", "b": "REFUSED", "c": "SILENT LOSS bases", "d": "pass",
            "e": "CRASHED", "f": "pass"}

    def test_case_names_are_unique_and_kinds_known(self):
        names = [n for n, *_ in board.synthetic_cases()]
        assert len(names) == len(set(names))
        assert {k for k, *_ in _CASES.values()} <= {"valid", "malformed"}


# ── end to end through the real readers ──────────────────────────────────────

def _item(tmp_path, name: str, *, private: bool = False) -> Any:
    kind, text, enc = _CASES[name]
    p = tmp_path / f"{name}.gb"
    p.write_bytes(text.encode(enc))
    it = board._make_item("synthetic", p, private=private, label=name, kind=kind)
    assert it is not None
    return it


class TestScoring:
    @pytest.mark.parametrize("name", ["canonical", "features_header_spacing",
                                      "spaced_locus_name", "multi_record"])
    def test_splicecraft_reads_them_faithfully_and_writes_to_spec(self, tmp_path, name):
        it = _item(tmp_path, name)
        census = genbank_census(board.decode_text(it.data()))
        row = board.score(it, "splicecraft", sc, census, tmp_path, write=True)
        assert row["outcome"] == "accepted" and row["issues"] == []
        assert row["write"] == "ok" and row["spec_problems"] == 0
        assert row["round_trip"] == "same"

    def test_a_private_file_leaves_no_name_or_path_in_the_saved_row(self, tmp_path):
        p = tmp_path / "secretfile.gb"
        p.write_text(_CASES["canonical"][1].replace(
            "LOCUS       pSynth          ", "LOCUS       Secretname in pVec "))
        it = board._make_item("local", p, private=True)
        assert it is not None and it.label == ""
        census = genbank_census(board.decode_text(it.data()))
        row = {"segment": "local", "key": it.key, "label": it.label, "_path": str(p),
               "_private": True, "by_reader": {
                   "biopython": board.score(it, "biopython", sc, census, tmp_path, write=False)}}
        saved = json.dumps(board.public_rows([row]), default=str)
        assert row["by_reader"]["biopython"]["outcome"] == "refused"
        assert "Secretname" not in saved and "secretfile" not in saved and str(tmp_path) not in saved

    def test_scrub_keeps_the_parsers_sentence_only(self):
        s = board._scrub
        assert s("Did not recognise the LOCUS line layout:\nLOCUS       Secret 30 bp") == (
            "Did not recognise the LOCUS line layout:")
        assert "SecretId" not in s("Expected sequence length 30, found 29 (SecretId).")
        assert "Secret" not in s("Invalid whitespace in 'Secret name' for LOCUS line")
        # A path names its file, and `.dna` originals are named after library
        # entries: the path must go whole, wherever it sits in the sentence.
        for msg in ("Could not parse X file /home/u/data/dna_originals/Secret_Plasmid.dna: bad",
                    "Could not parse X file C:\\Users\\u\\Secret_Plasmid.dna: bad",
                    "no record in ~/Downloads/Secret Name.gb"):
            out = s(msg)
            assert "Secret" not in out and "home" not in out and "Users" not in out, out
