"""SpliceCraft's own GenBank reader, `splicecraft_insdc` [INV-219].

The contract this file holds it to:

* On every file Biopython reads correctly, the SAME record: name, id,
  description, dbxrefs, every annotation (references too), every feature
  (type, location object for location object, qualifiers in order) and the
  sequence.
* Where Biopython refuses, crashes or reads wrong, what the file means, with
  an `InsdcParserWarning` and a count on `rec._insdc_repairs`.
* Never a location of None, never an exception but ValueError on bad input,
  and time linear in the size of the file.

Biopython is the ORACLE here, not a dependency of the reader's output: every
"same record" test reads the text with both and compares.
"""
from __future__ import annotations

import io
import random
import time
import warnings

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    AfterPosition, BeforePosition, CompoundLocation, ExactPosition,
    OneOfPosition, Reference, SeqFeature, SimpleLocation, WithinPosition)
from Bio.SeqRecord import SeqRecord

import splicecraft as sc
import splicecraft_insdc as H

_SEQ = ("atgaaagatctggaattcaaagggccctagaaagcatgcaaaatcgatgtcgacaaagaa"
        "ttcaaatcctaggaaaaggatccaaaactcgagcccaaaaaatttgggcccaaaatcgat")


def _origin(seq: str) -> str:
    return "\n".join(
        f"{i + 1:>9} " + " ".join(seq[j:j + 10]
                                   for j in range(i, min(i + 60, len(seq)), 10))
        for i in range(0, len(seq), 60))


def _locus(name: str = "pDEMO1", n: int = 120, topo: str = "circular") -> str:
    return (f"LOCUS       {name:<16}{n:>12} bp    DNA     {topo:<8} SYN "
            f"12-OCT-2026")


def _record(features: str, *, header: str = "", seq: str = _SEQ,
            locus: "str | None" = None) -> str:
    """A minimal standard record around a FEATURES block."""
    return "\n".join(x for x in (
        locus or _locus(n=len(seq)), header,
        "FEATURES             Location/Qualifiers", features,
        "ORIGIN", _origin(seq), "//", "") if x is not None)


_FULL = f"""{_locus()}
DEFINITION  Synthetic demonstration vector carrying GeneX under Pdemo, with a
            terminator.
ACCESSION   XX000001 XX000002
VERSION     XX000001.2
DBLINK      BioProject: PRJNA000001
            BioSample: SAMN00000001
KEYWORDS    demo; synthetic construct.
SOURCE      synthetic DNA construct
  ORGANISM  synthetic DNA construct
            other sequences; artificial sequences; vectors.
REFERENCE   1  (bases 1 to 120)
  AUTHORS   Doe,J. and Roe,R.
  CONSRTM   Demo Consortium
  TITLE     A demonstration vector for parser tests, with a title long enough
            to wrap onto a second line
  JOURNAL   J. Demo. 1 (1), 1-2 (2026)
   PUBMED   12345678
  REMARK    Erratum: none.
REFERENCE   2  (bases 1 to 60; 70 to 120)
  AUTHORS   Doe,J.
  TITLE     Direct Submission
  JOURNAL   Submitted (12-OCT-2026) Nowhere Lab, Demo Street 1, Demo City
COMMENT     A plain comment line.
            ##Assembly-Data-START##
            Assembly Method       :: DemoAssembler v. 1.0
            Sequencing Technology :: Sanger dideoxy sequencing
            ##Assembly-Data-END##
FEATURES             Location/Qualifiers
     source          1..120
                     /organism="synthetic DNA construct"
                     /mol_type="other DNA"
     promoter        <1..20
                     /label="Pdemo"
     CDS             join(21..50,61..>90)
                     /gene="GeneX"
                     /codon_start=1
                     /note="a ""quoted"" word and a long note that wraps over
                     more than one line of the feature table"
                     /translation="MKVLAAGIVLLLAKKLLLKVMMMMMM"
     misc_feature    complement(join(91..100,105..110))
                     /label="rev part"
     misc_feature    one-of(3,5)..10
     misc_feature    (10.12)..20
     misc_feature    40^41
                     /pseudo
     misc_feature    order(1..5,7..9)
     terminator      complement(111..120)
                     /label="Tdemo"
     rep_origin      XX000002.1:1..10
BASE COUNT       30 a     30 c     30 g     30 t
ORIGIN
{_origin(_SEQ)}
//
"""


def _quiet(fn, *a):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fn(*a)


def _bio(text: str) -> list:
    return _quiet(lambda t: list(SeqIO.parse(io.StringIO(t), "genbank")), text)


def _ours(text: str) -> list:
    return _quiet(H.parse_genbank, text)


def _seq_or_len(rec):
    try:
        return str(rec.seq)
    except Exception:                       # noqa: BLE001 — undefined sequence
        return ("<undefined>", len(rec.seq))


def _diff(a, b) -> str:
    """'' when two records are the same in everything a reader decides."""
    for attr in ("name", "id", "description", "dbxrefs"):
        if getattr(a, attr) != getattr(b, attr):
            return attr
    if set(a.annotations) != set(b.annotations):
        return f"annotation keys {sorted(set(a.annotations) ^ set(b.annotations))}"
    for k, va in a.annotations.items():
        vb = b.annotations[k]
        if k == "references":
            if len(va) != len(vb):
                return "reference count"
            for ra, rb in zip(va, vb):
                for at in ("authors", "consrtm", "title", "journal",
                           "medline_id", "pubmed_id", "comment"):
                    if getattr(ra, at) != getattr(rb, at):
                        return f"reference {at}"
                if [repr(x) for x in ra.location] != [repr(x) for x in rb.location]:
                    return "reference location"
        elif va != vb:
            return f"annotation {k}"
    if _seq_or_len(a) != _seq_or_len(b):
        return "sequence"
    if len(a.features) != len(b.features):
        return "feature count"
    for i, (fa, fb) in enumerate(zip(a.features, b.features)):
        if fa.type != fb.type:
            return f"feature {i} type"
        if repr(fa.location) != repr(fb.location):
            return f"feature {i} location"
        if list(fa.qualifiers.items()) != list(fb.qualifiers.items()):
            return f"feature {i} qualifiers"
    return ""


def _assert_same(text: str) -> list:
    bio, ours = _bio(text), _ours(text)
    assert len(bio) == len(ours) >= 1
    for a, b in zip(bio, ours):
        assert _diff(a, b) == "", _diff(a, b)
    return ours


# ── the same record as Biopython ────────────────────────────────────────────

class TestSameRecordAsBiopython:

    def test_a_full_ncbi_style_record(self):
        (rec,) = _assert_same(_FULL)
        # The reading itself, so "same" is not two readers agreeing on junk.
        assert rec.id == "XX000001.2" and rec.name == "pDEMO1"
        assert rec.annotations["topology"] == "circular"
        assert rec.dbxrefs == ["BioProject:PRJNA000001", "BioSample:SAMN00000001"]
        assert [str(r.location[0]) for r in rec.annotations["references"]] == [
            "[0:120]", "[0:60]"]
        assert rec.annotations["structured_comment"]["Assembly-Data"][
            "Assembly Method"] == "DemoAssembler v. 1.0"
        cds = rec.features[2]
        assert cds.qualifiers["note"] == [
            'a "quoted" word and a long note that wraps over more than one '
            "line of the feature table"]
        assert rec.features[-1].location.ref == "XX000002.1"
        assert rec.features[6].qualifiers == {"pseudo": [""]}
        assert rec._insdc_repairs == {}

    def test_a_protein_record(self):
        text = """LOCUS       DEMO_PROT                 30 aa            linear   PRI 12-OCT-2026
DEFINITION  demo protein.
ACCESSION   XP_000001
VERSION     XP_000001.1
DBSOURCE    REFSEQ: accession XM_000001.1
KEYWORDS    RefSeq.
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            other sequences; artificial sequences.
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="synthetic construct"
     Protein         1..30
                     /product="demo protein"
     Region          5..20
                     /region_name="Demo domain"
     Site            order(7,9,11)
                     /site_type="active"
     Bond            bond(3,25)
                     /bond_type="disulfide"
ORIGIN
        1 mkvlaagivl llakkllllk vmmmmmmmmm
//
"""
        (rec,) = _assert_same(text)
        assert rec.annotations["molecule_type"] == "protein"
        assert all(f.location.strand is None for f in rec.features)

    def test_a_contig_record_has_no_bases_of_its_own(self):
        text = (f"{_locus('pCONTIG', 500, 'linear')}\n"
                "FEATURES             Location/Qualifiers\n"
                "     source          1..500\n"
                "CONTIG      join(XX000001.1:1..250,XX000002.1:1..250)\n//\n")
        (rec,) = _assert_same(text)
        assert not rec.seq.defined and len(rec.seq) == 500
        assert rec.annotations["contig"].startswith("join(")

    def test_several_records_in_one_text(self):
        second = _FULL.replace("pDEMO1", "pDEMO2").replace("XX000001", "XX000009")
        recs = _assert_same(_FULL + second)
        assert [r.name for r in recs] == ["pDEMO1", "pDEMO2"]

    @pytest.mark.parametrize("seed", range(3))
    def test_random_records_written_by_biopython(self, seed):
        """Locations of every kind (fuzzy, within, one-of, between, mixed
        strand, order), wrapped and quoted values, references, protein
        records — written by Biopython's own writer, read by both."""
        rng = random.Random(seed)
        for i in range(60):
            buf = io.StringIO()
            _quiet(SeqIO.write, _random_record(rng, i), buf, "genbank")
            _assert_same(buf.getvalue())

    def test_what_splicecraft_itself_writes(self, realistic_plasmid):
        """Every library entry is SpliceCraft's own writer's output: arrowless
        features, an origin-spanning feature, a primer with its sequence."""
        rec = realistic_plasmid
        n = len(rec.seq)
        rec.features.append(SeqFeature(CompoundLocation(
            [SimpleLocation(n - 30, n, 1), SimpleLocation(0, 40, 1)]),
            type="misc_feature", qualifiers={"label": ["across the origin"]}))
        rec.features.append(SeqFeature(SimpleLocation(100, 140, 0),
                                       type="misc_feature",
                                       qualifiers={"label": ["no arrow"]}))
        rec.features.append(SeqFeature(
            SimpleLocation(200, 224, -1), type="primer_bind",
            qualifiers={"label": ["oDemo-R"],
                        "primer_seq": ["GATTACAGATTACAGATTACAGAT"]}))
        text = sc._record_to_gb_text(rec)
        _assert_same(text)


def _random_record(rng: random.Random, i: int) -> SeqRecord:
    protein = rng.random() < 0.1
    n = rng.randint(80, 1500)
    seq = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY" if protein else "ACGT")
                  for _ in range(n))
    words = ("alpha beta gamma 5' 3' (partial) [x] a=b c/d e;f 10..20 "
             '"quoted" verylongwordwithoutanyspacesthatmustbewrappedsomewhere'
             ).split()

    def text(k: int) -> str:
        return " ".join(rng.choice(words) for _ in range(rng.randint(1, k)))

    rec = SeqRecord(Seq(seq), id=f"rec{i}.1", name=f"rec{i}",
                    description=text(10).rstrip("."))
    rec.annotations.update(
        molecule_type="protein" if protein else rng.choice(["DNA", "mRNA"]),
        data_file_division="SYN", date="12-OCT-2026",
        accessions=[f"rec{i}"], sequence_version=1,
        keywords=[text(3)], source=text(4), organism=text(3),
        taxonomy=["Bacteria", "Pseudomonadota"])
    if not protein:
        rec.annotations["topology"] = rng.choice(["linear", "circular"])
    ref = Reference()
    ref.location = [SimpleLocation(0, n)]
    ref.authors, ref.title, ref.journal = text(6), text(25), text(6)
    rec.annotations["references"] = [ref]
    rec.annotations["comment"] = "\n".join(text(12) for _ in range(2))
    for _ in range(rng.randint(0, 14)):
        strand = None if protein else rng.choice([1, -1])
        s = rng.randrange(0, n - 20)
        e = rng.randint(s + 1, min(n, s + 300))
        kind = rng.random()
        if kind < 0.5:
            loc = SimpleLocation(s, e, strand)
        elif kind < 0.6:
            loc = SimpleLocation(BeforePosition(s), AfterPosition(e), strand)
        elif kind < 0.67 and e - s > 4:
            loc = SimpleLocation(WithinPosition(s, left=s, right=s + 2),
                                 WithinPosition(e, left=e - 2, right=e), strand)
        elif kind < 0.74 and e - s > 4:
            loc = SimpleLocation(OneOfPosition(
                s, choices=[ExactPosition(s), ExactPosition(s + 2)]),
                ExactPosition(e), strand)
        elif kind < 0.8:
            loc = SimpleLocation(s, s, strand)
        else:
            cuts = sorted(rng.sample(range(1, n), 4))
            parts = [SimpleLocation(cuts[0], cuts[1], strand),
                     SimpleLocation(cuts[2], cuts[3], strand)]
            if strand == -1:
                parts.reverse()
            loc = CompoundLocation(parts, operator=rng.choice(["join", "order"]))
        f = SeqFeature(loc, type=rng.choice(
            ["Region", "Site"] if protein else
            ["CDS", "misc_feature", "promoter", "3'UTR", "-10_signal"]))
        for _ in range(rng.randint(0, 4)):
            k = rng.choice(["label", "note", "codon_start", "translation"])
            v = (str(rng.randint(1, 3)) if k == "codon_start" else
                 "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY")
                         for _ in range(rng.randint(5, 200)))
                 if k == "translation" else text(30))
            f.qualifiers.setdefault(k, []).append(v)
        if rng.random() < 0.1:
            f.qualifiers["pseudo"] = [""]
        rec.features.append(f)
    return rec


# ── where Biopython refuses, crashes or reads wrong ─────────────────────────

class TestWhereBiopythonFails:

    def test_spaced_locus_name(self):
        text = _record("     source          1..120", locus=(
            "LOCUS       GeneX in pShuttle Amp High Copy        120 bp DNA     "
            "circular   12-OCT-2026"))
        with pytest.raises(ValueError):
            _bio(text)
        (rec,) = _ours(text)
        assert rec.name == "GeneX in pShuttle Amp High Copy"
        assert rec.annotations["topology"] == "circular"
        assert rec.annotations["date"] == "12-OCT-2026"
        assert rec._insdc_repairs.get("locus_line") == 1

    @pytest.mark.parametrize("locus", [
        "LOCUS       pDemo 120 bp DNA circular 12-OCT-2026",       # 7 tokens
        "LOCUS       pDemo            120 bp    DNA     linear",    # no division
        "LOCUS       pDemo                     120 bp    DNA      circular SYN 12-OCT-2026",
    ])
    def test_topology_is_read_wherever_the_field_sits(self, locus):
        """Biopython loses it on a 7-token line, refuses it one column over."""
        (rec,) = _ours(_record("", locus=locus))
        assert rec.annotations["topology"] == ("linear" if "linear" in locus
                                               else "circular")
        assert len(rec.seq) == 120

    def test_origin_spanning_feature_on_a_record_not_marked_circular(self):
        """Biopython: location None, and the next save of the record crashed."""
        text = _record("     misc_feature    100..20\n"
                       '                     /label="wraps"',
                       locus=_locus(topo="linear"))
        (bio,) = _bio(text)
        assert bio.features[0].location is None
        (rec,) = _ours(text)
        loc = rec.features[0].location
        assert [(int(p.start), int(p.end)) for p in loc.parts] == [(99, 120), (0, 20)]
        assert rec._insdc_repairs == {"origin_spanning_not_circular": 1}

    def test_negative_position_runs_off_the_start(self):
        (rec,) = _ours(_record("     misc_feature    -5..20"))
        loc = rec.features[0].location
        assert isinstance(loc.start, BeforePosition) and int(loc.start) == 0
        assert int(loc.end) == 20
        assert rec._insdc_repairs == {"negative_position": 1}

    def test_unnumbered_sequence_lines_keep_every_base(self):
        """Biopython cuts ten COLUMNS off every line: 10 bases a line lost,
        with only a warning. (An unbroken 60-letter line it refuses.)"""
        seq = _SEQ
        text = (_locus() + "\nFEATURES             Location/Qualifiers\n"
                "ORIGIN\n" + "\n".join(
                    " ".join(seq[j:j + 10] for j in range(i, i + 60, 10))
                    for i in range(0, 120, 60)) + "\n//\n")
        assert len(_bio(text)[0].seq) == 100
        (rec,) = _ours(text)
        assert str(rec.seq) == seq.upper()
        assert rec._insdc_repairs == {"unnumbered_sequence": 1}

    def test_position_numbers_at_the_end_of_the_line(self):
        text = (_locus() + "\nFEATURES             Location/Qualifiers\n"
                "ORIGIN\n" + "\n".join(f"{_SEQ[i:i + 60]} {i + 60}"
                                       for i in range(0, 120, 60)) + "\n//\n")
        assert str(_ours(text)[0].seq) == _SEQ.upper()

    def test_tab_indented_feature_table(self):
        """Biopython reads this by column: a garbage key or an IndexError."""
        text = _record("\tsource\t\t1..120\n"
                       '\t\t\t/organism="synthetic DNA construct"\n'
                       "\tCDS\t\tcomplement(10..50)\n"
                       '\t\t\t/label="GeneX"\n'
                       '\t\t\t/note="wrapped over\n'
                       '\t\t\ttwo lines"')
        (rec,) = _ours(text)
        assert [(f.type, str(f.location)) for f in rec.features] == [
            ("source", "[0:120](+)"), ("CDS", "[9:50](-)")]
        assert rec.features[1].qualifiers == {
            "label": ["GeneX"], "note": ["wrapped over two lines"]}
        assert rec._insdc_repairs == {"feature_table_layout": 1}

    def test_a_whole_table_shifted_left(self):
        text = _record("source          1..120\n"
                       '                /organism="synthetic DNA construct"\n'
                       "promoter        <1..20\n"
                       '                /label="Pdemo"')
        (rec,) = _ours(text)
        assert [f.type for f in rec.features] == ["source", "promoter"]
        assert rec.features[1].qualifiers == {"label": ["Pdemo"]}

    def test_a_quote_that_never_closes_ends_at_the_next_feature(self):
        text = _record("     misc_feature    1..10\n"
                       '                     /note="never closed\n'
                       "     promoter        11..20\n"
                       '                     /label="Pdemo"')
        (rec,) = _ours(text)
        assert [f.type for f in rec.features] == ["misc_feature", "promoter"]
        assert rec.features[0].qualifiers == {"note": ["never closed"]}
        assert rec.features[1].qualifiers == {"label": ["Pdemo"]}
        assert rec._insdc_repairs == {"unclosed_quote": 1}

    @pytest.mark.parametrize("line", ['/product "a protein"',
                                      '/product ="a protein"',
                                      '/product  =  "a protein"'])
    def test_a_missing_or_spaced_equals_sign(self, line):
        """Biopython makes the whole line one garbage qualifier NAME."""
        (rec,) = _ours(_record("     CDS             1..30\n"
                               f"                     {line}"))
        assert rec.features[0].qualifiers == {"product": ["a protein"]}

    def test_a_nested_operator_is_read_as_one(self):
        (rec,) = _ours(_record(
            "     misc_feature    join(1..5,join(10..15,20..25))"))
        assert [(int(p.start), int(p.end))
                for p in rec.features[0].location.parts] == [(0, 5), (9, 15), (19, 25)]

    def test_both_kinds_of_fuzzy_end_at_once(self):
        """`one-of(..)..(a.b)`: Biopython gives up (location None)."""
        text = _record("     misc_feature    one-of(3,5)..(10.12)")
        assert _bio(text)[0].features[0].location is None
        loc = _ours(text)[0].features[0].location
        assert isinstance(loc.start, OneOfPosition) and int(loc.start) == 2
        assert isinstance(loc.end, WithinPosition) and int(loc.end) == 12

    def test_a_utf8_bom_in_front_of_the_text(self):
        """Biopython: zero records. File loads strip it on decode; text
        handed in directly (an agent `gb_text`) did not."""
        assert _bio("\ufeff" + _FULL) == []
        (rec,) = _ours("\ufeff" + _FULL)
        assert _diff(_bio(_FULL)[0], rec) == ""

    @pytest.mark.parametrize("eol", ["\r\n", "\r"])
    def test_other_line_endings(self, eol):
        (rec,) = _ours(_FULL.replace("\n", eol))
        assert _diff(_bio(_FULL)[0], rec) == ""

    def test_a_lower_case_or_tabbed_features_line(self):
        for header in ("features             Location/Qualifiers",
                       "FEATURES\tLocation/Qualifiers"):
            text = _record("     source          1..120").replace(
                "FEATURES             Location/Qualifiers", header)
            assert [f.type for f in _ours(text)[0].features] == ["source"]


# ── corrupt input: never worse than Biopython, never a crash ────────────────

class TestDamagedFiles:
    """Each case below came out of mutation fuzzing (2026-10-13) and was read
    WORSE than Biopython reads it, or crashed."""

    def test_a_misspelt_sub_keyword_does_not_pour_into_the_field_above(self):
        text = _FULL.replace("  AUTHORS   Doe,J. and Roe,R.",
                             "  AUTHO;S   Doe,J. and Roe,R.")
        (bio,) = _bio(text)
        (rec,) = _ours(text)
        assert _diff(bio, rec) == ""
        assert [str(x) for x in rec.annotations["references"][0].location] == [
            "[0:120]"]

    def test_a_missing_features_line_does_not_pour_into_the_header(self):
        """The table becomes unknown header lines — skipped, as Biopython
        skips them — not one enormous COMMENT."""
        text = _FULL.replace("FEATURES             Location/Qualifiers\n", "")
        (bio,) = _bio(text)
        (rec,) = _ours(text)
        assert _diff(bio, rec) == ""
        assert rec.annotations["comment"] == "A plain comment line."

    def test_a_lineage_line_two_columns_short_still_continues(self):
        text = _FULL.replace(
            "            other sequences; artificial sequences; vectors.",
            "          other sequences; artificial sequences; vectors.")
        (rec,) = _ours(text)
        assert rec.annotations["taxonomy"] == [
            "other sequences", "artificial sequences", "vectors"]

    def test_stray_qualifiers_above_the_first_feature(self):
        """The first key line lost: its wrapped qualifier used to set the KEY
        column, and every wrapped line after it became a feature of its own,
        taking its text with it."""
        text = _record(
            '                     /note="a long note whose key line was\n'
            '                     lost to a damaged file"\n'
            "     CDS             1..30\n"
            '                     /product="a protein whose name wraps\n'
            '                     onto a second line"\n'
            "                     /label=an unquoted label that a writer\n"
            "                     wrapped\n"
            "     promoter        40..60")
        (rec,) = _ours(text)
        assert [f.type for f in rec.features] == ["CDS", "promoter"]
        assert rec.features[0].qualifiers == {
            "product": ["a protein whose name wraps onto a second line"],
            "label": ["an unquoted label that a writer wrapped"]}
        assert _diff(_bio(text)[0], rec) == ""

    @pytest.mark.parametrize("junk", ['"', "x", "é"])
    def test_junk_where_the_indent_was(self, junk):
        text = _record("     CDS             1..30\n"
                       '                     /label="GeneX"\n'
                       f'{junk}                    /codon_start=1\n'
                       '                     /product="a protein"')
        (rec,) = _ours(text)
        assert [f.type for f in rec.features] == ["CDS"]
        assert rec.features[0].qualifiers == {
            "label": ["GeneX"], "codon_start": ["1"], "product": ["a protein"]}

    def test_junk_before_the_first_qualifier_does_not_move_the_column(self):
        """The qualifier column is learnt from the first qualifier line. With
        junk in front of THAT one, the column was taken from the junk (0) and
        every feature line after it read as text glued into its value."""
        text = _record("     CDS             1..30\n"
                       '"                    /label="GeneX"\n'
                       "     promoter        40..60\n"
                       '                     /label="Pdemo"\n'
                       "     terminator      70..90")
        (rec,) = _ours(text)
        assert [(f.type, f.qualifiers) for f in rec.features] == [
            ("CDS", {"label": ["GeneX"]}), ("promoter", {"label": ["Pdemo"]}),
            ("terminator", {})]

    def test_a_header_keyword_inside_the_feature_table(self):
        text = _record("     CDS             1..30\n"
                       '                     /label="GeneX"\n'
                       "COMMENT     written after the features\n"
                       "            and continued\n"
                       "     promoter        40..60")
        (rec,) = _ours(text)
        assert [f.type for f in rec.features] == ["CDS", "promoter"]
        assert rec.features[0].qualifiers == {"label": ["GeneX"]}
        # Skipped as text — not read as two features with unreadable
        # locations, which is what the user would have been told.
        assert rec._insdc_repairs == {}

    def test_a_custom_key_after_a_header_line_in_the_table_is_kept(self):
        text = _record("     CDS             1..30\n"
                       "COMMENT     written after the features\n"
                       "            and continued\n"
                       "     my_part         40..60\n"
                       '                     /label="kept"')
        (rec,) = _ours(text)
        assert [(f.type, f.qualifiers) for f in rec.features] == [
            ("CDS", {}), ("my_part", {"label": ["kept"]})]

    def test_a_duplicated_locus_line_is_still_one_record(self):
        (rec,) = _ours(_locus() + "\n" + _FULL)
        assert str(rec.seq) == _SEQ.upper()

    def test_a_repeated_origin_line_is_not_bases(self):
        text = _FULL.replace("ORIGIN\n", "ORIGIN\nORIGIN\n")
        assert str(_ours(text)[0].seq) == _SEQ.upper()

    def test_a_no_break_space_in_the_sequence_is_white_space(self):
        text = _FULL.replace("atgaaagatc tggaattcaa", "atgaaagatc\u00a0tggaattcaa")
        assert str(_ours(text)[0].seq) == _SEQ.upper()

    def test_a_carriage_return_inside_a_value_is_data(self):
        """CRLF and CR-only line ends are line ends; a lone CR inside a line of
        LF text is DATA (a `.dna` note typed with one), as Biopython reads it.
        Normalised into a line break, it changed the note and made an
        unedited old `.dna` import look edited."""
        text = _record("     misc_feature    1..10\n"
                       '                     /note="first\rsecond"')
        (rec,) = _assert_same(text)
        assert rec.features[0].qualifiers["note"] == ["first\rsecond"]

    def test_a_non_ascii_letter_in_the_sequence_is_named(self):
        text = _FULL.replace("atgaaagatc", "atgaaégatc")
        with pytest.raises(ValueError, match="'é', which is not a sequence letter"):
            H.parse_genbank(text)

    def test_a_record_cut_off_before_its_sequence_is_refused(self):
        cut = _FULL[:_FULL.index("     terminator")]
        with pytest.raises(ValueError, match="ends before its sequence"):
            H.parse_genbank(cut)
        # …but a WGS master record has no bases by design.
        wgs = (_locus("pWGS", 0, "linear") + "\nWGS         XXXX01000001-XXXX01000100\n")
        (rec,) = _ours(wgs)
        assert rec.annotations["wgs"] == ["XXXX01000001", "XXXX01000100"]

    def test_a_record_without_its_end_line_still_reads(self):
        """Biopython warns and reads it; so does this."""
        (rec,) = _ours(_FULL.replace("\n//\n", "\n"))
        assert str(rec.seq) == _SEQ.upper()

    def test_a_structured_comment_reopened_after_another_block(self):
        """Biopython raises KeyError here."""
        text = _FULL.replace(
            "            ##Assembly-Data-END##",
            "            ##Assembly-Data-END##\n"
            "            ##Other-START##\n"
            "            Key :: value\n"
            "            ##Other-END##\n"
            "            ##Assembly-Data-START##\n"
            "            a continuation line")
        (rec,) = _ours(text)
        assert rec.annotations["structured_comment"]["Other"] == {"Key": "value"}

    @pytest.mark.parametrize("locus", ["LOCUS       pX ² bp DNA linear",
                                       "LOCUS       pX 1²0 bp DNA linear"])
    def test_a_unicode_digit_in_the_locus_length(self, locus):
        """`'²'.isdigit()` is True and `int('²')` raises."""
        recs = _ours(locus + "\nORIGIN\n        1 acgt\n//\n")
        assert str(recs[0].seq) == "ACGT"

    def test_unreadable_locations_drop_their_feature_never_none(self):
        text = _record("     misc_feature    1..10\n"
                       "     misc_feature    10..wat\n"
                       "     misc_feature    join(1..5,\n"
                       "     promoter        20..30")
        (rec,) = _ours(text)
        assert [f.type for f in rec.features] == ["misc_feature", "promoter"]
        assert rec._insdc_repairs == {"unreadable_location": 2}

    @pytest.mark.parametrize("opener", ["join(", "complement(", "order("])
    def test_a_location_nested_thousands_deep(self, opener):
        """Recursed into RecursionError — not a ValueError, so it got past
        every caller (`_gb_text_to_record` is how every library entry and
        agent `gb_text` loads)."""
        deep = opener * 3000 + "1..2" + ")" * 3000
        (rec,) = _ours(_record(f"     misc_feature    {deep}\n"
                               "     promoter        20..30"))
        assert [f.type for f in rec.features] == ["promoter"]
        with pytest.raises(ValueError, match="nested too deeply"):
            H.parse_location(deep)

    def test_a_long_wrapped_value_reads_in_linear_time(self):
        """Values were built with `+=`: 200,000 lines took 45 s."""
        lines = "\n".join("                     word word word word"
                          for _ in range(100_000))
        text = _record('     misc_feature    1..10\n'
                       '                     /note="start\n' + lines + '"')
        t0 = time.perf_counter()
        (rec,) = _ours(text)
        assert time.perf_counter() - t0 < 5
        assert rec.features[0].qualifiers["note"][0].count("word") == 400_000

    def test_seeded_mutations_never_crash(self):
        """Damage a valid file at random: ValueError or records, never another
        exception, a None location or a key no writer could mean."""
        rng = random.Random(20261013)
        base = _FULL.split("\n")
        for _ in range(400):
            lines = list(base)
            for _ in range(rng.randint(1, 4)):
                if not lines:
                    lines = ["LOCUS"]
                i = rng.randrange(len(lines))
                op = rng.randrange(6)
                if op == 0:
                    del lines[i]
                elif op == 1:
                    lines.insert(i, lines[i])
                elif op == 2:
                    lines[i] = lines[i].replace("    ", "\t", 1)
                elif op == 3:
                    lines[i] = " " * rng.randint(1, 4) + lines[i]
                elif op == 4 and lines[i]:
                    j = rng.randrange(len(lines[i]))
                    lines[i] = (lines[i][:j] + rng.choice('"()<>.,:^/=\t?9é')
                                + lines[i][j + 1:])
                else:
                    del lines[i:]
            try:
                recs = _ours("\n".join(lines))
            except ValueError:
                continue
            for rec in recs:
                for f in rec.features:
                    assert f.location is not None
                    assert f.type and not any(c in f.type for c in ' /="\t')


# ── through the app's own doors ─────────────────────────────────────────────

class TestAppDoors:

    def test_origin_spanning_feature_on_a_linear_file_saves_and_exports(
            self, tmp_path):
        """Finding 2 of the format scoreboard [INV-218]: the record loaded
        with location None, then the next library save AND the export raised
        AttributeError."""
        p = tmp_path / "wraps.gb"
        p.write_text(_record("     enhancer        100..20\n"
                             '                     /label="wraps"',
                             locus=_locus(topo="linear")))
        rec = sc.load_genbank(str(p))
        (feat,) = rec.features
        assert feat.location is not None and len(feat.location.parts) == 2
        back = sc._gb_text_to_record(sc._record_to_gb_text(rec), cache=False)
        assert [(int(x.start), int(x.end)) for x in back.features[0].location.parts] \
            == [(99, 120), (0, 20)]
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(rec, out)
        again = sc.load_genbank(str(out))
        assert str(again.features[0].location) == str(feat.location)

    def test_unnumbered_sequence_file(self, tmp_path):
        p = tmp_path / "raw.gb"
        p.write_text(_locus() + "\nFEATURES             Location/Qualifiers\n"
                     "ORIGIN\n" + _SEQ[:60] + "\n" + _SEQ[60:] + "\n//\n")
        assert str(sc.load_genbank(str(p)).seq) == _SEQ.upper()

    def test_tab_indented_file(self, tmp_path):
        p = tmp_path / "tabs.gb"
        p.write_text(_record("\tCDS\t\t10..50\n\t\t\t/label=\"GeneX\""))
        rec = sc.load_genbank(str(p))
        assert [(f.type, f.qualifiers.get("label")) for f in rec.features] == [
            ("CDS", ["GeneX"])]

    def test_agent_text_with_a_bom_and_cr_line_ends(self):
        rec = sc._gb_text_to_record("\ufeff" + _FULL.replace("\n", "\r"),
                                    cache=False)
        assert str(rec.seq) == _SEQ.upper()
        assert len(rec.features) == 10

    @staticmethod
    def _open_through_the_app(rec) -> "tuple[list, list]":
        """Drive `PlasmidApp._on_open_file_done`, stubbing only its I/O: the
        toast and the library write."""
        app = sc.PlasmidApp.__new__(sc.PlasmidApp)
        notes: list = []
        imported: list = []
        app.notify = lambda msg, **kw: notes.append((msg, kw.get("severity")))
        app._import_and_persist = imported.append
        sc.PlasmidApp._on_open_file_done(app, rec)
        assert imported == [rec]
        return notes, imported

    def test_file_open_says_when_a_feature_was_left_out(self, tmp_path):
        """A dropped feature was only a log line: silent loss to the user."""
        p = tmp_path / "damaged.gb"
        p.write_text(_record("     misc_feature    1..10\n"
                             "     misc_feature    10..wat\n"
                             "     enhancer        100..20",
                             locus=_locus(topo="linear")))
        notes, _ = self._open_through_the_app(sc.load_genbank(str(p)))
        ((msg, severity),) = notes
        assert severity == "warning"
        assert ("1 feature had a location that could not be read and was "
                "left out") in msg
        assert "1 feature runs across the origin" in msg

    def test_file_open_is_quiet_about_layout_repairs(self, tmp_path):
        p = tmp_path / "tabs.gb"
        p.write_text(_record("\tCDS\t\t10..50\n\t\t\t/label=\"GeneX\""))
        rec = sc.load_genbank(str(p))
        assert rec._insdc_repairs == {"feature_table_layout": 1}
        notes, _ = self._open_through_the_app(rec)
        assert notes == []

    def test_repair_note_counts(self):
        rec = SeqRecord(Seq("A"))
        rec._insdc_repairs = {"unreadable_location": 2,
                              "origin_spanning_not_circular": 3}
        assert H.repair_note(rec) == (
            "2 features had a location that could not be read and were left "
            "out; 3 features run across the origin, but the file marks the "
            "molecule linear (details in the log).")
        assert H.repair_note(SeqRecord(Seq("A"))) is None

    def test_cr_only_text_with_a_length_in_a_note(self):
        """With no `\\n` in it the LOCUS repair saw the whole text as the
        LOCUS line, and a note saying "500 bp" made it rewrite the record
        into one LOCUS line — refused as cut short."""
        text = _record("     misc_feature    1..60\n"
                       '                     /note="a 500 bp insert"')
        rec = sc._gb_text_to_record(text.replace("\n", "\r"), cache=False)
        assert str(rec.seq) == _SEQ.upper()
        assert rec.features[0].qualifiers["note"] == ["a 500 bp insert"]
        assert rec.name == "pDEMO1"

    def test_a_cut_off_text_says_so(self):
        with pytest.raises(ValueError, match="ends before its sequence"):
            sc._gb_text_to_record(_FULL[:_FULL.index("ORIGIN")], cache=False)
