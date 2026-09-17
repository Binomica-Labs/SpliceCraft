"""Spec-conformance guards for every interchange format SpliceCraft writes.

These tests exist because "Biopython reads it back" is NOT the same claim as
"any tool that reads standard GenBank can read it". The whole 2026-09-17 audit
came out of that gap: every real NCBI record with a CDS exported with the
entire protein translation on ONE line (7,132 characters for SARS-CoV-2
ORF1ab), which Biopython round-trips perfectly and an 80-column reader does
not.

So the assertions here are written against the FORMAT SPECS — the NCBI GenBank
release-notes layout and the INSDC Feature Table definition, the ENA flat-file
ID line, GFF3 1.26, RFC 4180 — not against what SpliceCraft happens to emit.
`_gb_problems` / `_gff3_problems` below are deliberately small independent
re-implementations of the rules, not calls back into the code under test.
"""

from __future__ import annotations

import pathlib
import re

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import splicecraft as sc


# ── independent format checkers ───────────────────────────────────────────────

_GB_DIVISIONS = set(
    "PRI ROD MAM VRT INV PLN BCT VRL PHG SYN UNA EST PAT STS GSS HTG HTC ENV "
    "CON TSA UNK".split()
)


def _gb_problems(text: str) -> "list[str]":
    """Check GenBank flat-file text against the published layout."""
    out: list[str] = []
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines = lines[:-1]
    for i, ln in enumerate(lines, 1):
        if len(ln) > 80:
            out.append(f"line {i}: {len(ln)} columns (max 80): {ln[:48]!r}")
        if not ln.isascii():
            bad = [c for c in ln if not c.isascii()]
            out.append(f"line {i}: non-ASCII {bad[:4]!r}")
        if "\t" in ln:
            out.append(f"line {i}: TAB character")
    if not lines:
        return ["empty file"] + out
    loc = lines[0]
    if not loc.startswith("LOCUS       "):
        out.append("no LOCUS line")
        return out
    if not re.search(r"(?<= )\d+$", loc[12:40]):
        out.append(f"LOCUS: length not right-justified to col 40: {loc[12:40]!r}")
    if loc[41:43] not in ("bp", "aa"):
        out.append(f"LOCUS: cols 42-43 = {loc[41:43]!r}, expected bp/aa")
    if loc[55:63].strip() not in ("linear", "circular", ""):
        out.append(f"LOCUS: topology field = {loc[55:63]!r}")
    if loc[64:67].strip() and loc[64:67] not in _GB_DIVISIONS:
        out.append(f"LOCUS: division {loc[64:67]!r} unknown")
    if loc[68:79].strip() and not re.fullmatch(r"\d{2}-[A-Z]{3}-\d{4}",
                                               loc[68:79]):
        out.append(f"LOCUS: date {loc[68:79]!r} not DD-MMM-YYYY")
    if "FEATURES             Location/Qualifiers" not in text:
        out.append("no FEATURES header")
    if not text.rstrip().endswith("//"):
        out.append("no // terminator")
    in_ft = False
    for i, ln in enumerate(lines, 1):
        if ln.startswith("FEATURES"):
            in_ft = True
            continue
        if in_ft and ln[:6] in ("ORIGIN", "//"):
            break
        if not in_ft:
            continue
        if ln[:5] != "     ":
            out.append(f"line {i}: feature line not indented 5: {ln[:12]!r}")
            continue
        key = ln[5:20]
        if key.strip():
            if len(key.strip()) > 15:
                out.append(f"line {i}: feature key >15 chars: {key.strip()!r}")
            if ln[20:21] not in (" ", ""):
                out.append(f"line {i}: no space before the location column")
        elif ln[:21] != " " * 21:
            out.append(f"line {i}: qualifier not indented 21: {ln[:21]!r}")
    # ORIGIN block
    if "\nORIGIN" in text:
        body = text.split("\nORIGIN", 1)[1].split("\n", 1)[1]
        pos = 1
        for ln in body.split("\n"):
            if ln.strip() in ("//", ""):
                break
            if ln[:9] != str(pos).rjust(9):
                out.append(f"ORIGIN: expected base counter {pos} in cols 1-9, "
                           f"got {ln[:9]!r}")
                break
            blocks = [b for b in ln[10:].split(" ") if b]
            if any(not re.fullmatch(r"[a-z]+", b) for b in blocks):
                out.append(f"ORIGIN: non-lowercase block in {ln[:40]!r}")
                break
            pos += sum(len(b) for b in blocks)
    return out


_GFF3_RESERVED = {"ID", "Name", "Alias", "Parent", "Target", "Gap",
                  "Derives_from", "Note", "Dbxref", "Ontology_term",
                  "Is_circular"}


def _gff3_problems(text: str) -> "list[str]":
    """Check GFF3 text against specification 1.26."""
    out: list[str] = []
    lines = text.split("\n")
    if not lines[0].startswith("##gff-version 3"):
        out.append("first line is not '##gff-version 3'")
    ids: dict = {}
    for i, raw in enumerate(lines, 1):
        if not raw or raw.startswith("#"):
            continue
        cols = raw.split("\t")
        if len(cols) != 9:
            out.append(f"line {i}: {len(cols)} columns, GFF3 needs exactly 9")
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = cols
        if not re.fullmatch(r"[A-Za-z0-9.:^*$@!+_?|%-]+", seqid):
            out.append(f"line {i}: seqid {seqid!r} needs percent-encoding")
        if not source or not ftype:
            out.append(f"line {i}: empty source/type column")
        if re.search(r"\s", ftype):
            out.append(f"line {i}: whitespace in the type column {ftype!r}")
        if not (start.isdigit() and end.isdigit()):
            out.append(f"line {i}: non-integer coordinates {start!r} {end!r}")
        elif int(start) < 1 or int(end) < int(start):
            out.append(f"line {i}: bad span {start}..{end}")
        if strand not in ("+", "-", ".", "?"):
            out.append(f"line {i}: strand {strand!r}")
        if phase not in (".", "0", "1", "2"):
            out.append(f"line {i}: phase {phase!r}")
        if ftype == "CDS" and phase == ".":
            out.append(f"line {i}: CDS without a phase")
        seen: list = []
        for kv in attrs.split(";"):
            if not kv or kv == ".":
                continue
            if "=" not in kv:
                out.append(f"line {i}: attribute {kv!r} has no '='")
                continue
            k, _, v = kv.partition("=")
            if k in seen:
                out.append(f"line {i}: duplicate attribute tag {k!r}")
            seen.append(k)
            if any(c in v for c in "\t\n\r"):
                out.append(f"line {i}: raw control char in attribute {k!r}")
            for m in re.finditer("%", v):
                if not re.match(r"%[0-9A-Fa-f]{2}", v[m.start():m.start() + 3]):
                    out.append(f"line {i}: bare '%' in attribute {k!r}")
                    break
            if k == "ID":
                ids.setdefault(v, set()).add(ftype)
    for gid, types in ids.items():
        if len(types) > 1:
            out.append(f"ID={gid!r} shared by different types {sorted(types)}")
    return out


# ── fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def plasmid():
    """A circular record shaped like the real thing: a source feature, a CDS
    with a translation long enough to wrap, an origin-spanning wrap feature,
    a reverse feature and an arrowless one."""
    seq = "ATGCGTACGT" * 200            # 2000 bp
    rec = SeqRecord(Seq(seq), id="pTest", name="pTest",
                    description="conformance fixture")
    rec.annotations.update(molecule_type="DNA", topology="circular",
                           data_file_division="SYN", date="17-SEP-2026")
    rec.features = [
        SeqFeature(FeatureLocation(0, 2000, strand=1), type="source",
                   qualifiers={"organism": ["synthetic construct"],
                               "mol_type": ["other DNA"]}),
        SeqFeature(FeatureLocation(100, 700, strand=1), type="CDS",
                   qualifiers={"label": ["bla"], "codon_start": ["1"],
                               "transl_table": ["11"],
                               "translation": ["M" * 199]}),
        SeqFeature(FeatureLocation(800, 900, strand=-1), type="promoter",
                   qualifiers={"label": ["Prev"], "note": ["reverse"]}),
        SeqFeature(FeatureLocation(1000, 1100, strand=0),
                   type="misc_feature", qualifiers={"label": ["arrowless"]}),
        SeqFeature(CompoundLocation([FeatureLocation(1900, 2000, strand=1),
                                     FeatureLocation(0, 50, strand=1)]),
                   type="misc_feature", qualifiers={"label": ["wrap"]}),
    ]
    return rec


# ── GenBank ───────────────────────────────────────────────────────────────────

class TestGenBankConformance:

    def test_export_is_spec_clean(self, plasmid, tmp_path):
        out = tmp_path / "p.gb"
        sc._export_genbank_to_path(plasmid, out)
        assert _gb_problems(out.read_text()) == []

    def test_save_path_is_spec_clean(self, plasmid):
        """`_record_to_gb_text` is not just the exporter — it is how every
        library entry is STORED, so a violation here is baked into the
        user's own data, not only into a file they hand out."""
        assert _gb_problems(sc._record_to_gb_text(plasmid)) == []

    def test_translation_wraps_at_eighty_columns(self, plasmid):
        """THE regression this whole file exists for. `/translation` is the
        one qualifier every INSDC reader rejoins with no separator, which is
        why NCBI wraps it mid-token — so refusing to break it put every real
        CDS on a single multi-thousand-column line."""
        text = sc._record_to_gb_text(plasmid)
        assert max(len(ln) for ln in text.split("\n")) <= 80
        back = sc._gb_text_to_record(text, cache=False)
        cds = next(f for f in back.features if f.type == "CDS")
        assert cds.qualifiers["translation"] == ["M" * 199]

    def test_long_primer_seq_wraps_and_survives(self, tmp_path):
        """A 70-mer Gibson primer overflows the qualifier line. Whitespace in
        an oligo is never data, so this one wraps — and both load paths
        normalise it back."""
        oligo = "GCGCCGTCTCAAATGAATAAATGTATTCCAATGATAATTAATGGAATGATCCGGTAAGCTTACGTACGTA"
        rec = SeqRecord(Seq("ATGC" * 100), id="pP", name="pP", description="d")
        rec.annotations.update(molecule_type="DNA", topology="circular")
        rec.features = [SeqFeature(FeatureLocation(10, 80, strand=1),
                                   type="primer_bind",
                                   qualifiers={"label": ["p1"],
                                               "primer_seq": [oligo]})]
        out = tmp_path / "p.gb"
        sc._export_genbank_to_path(rec, out)
        assert _gb_problems(out.read_text()) == []
        for back in (sc.load_genbank(str(out)),
                     sc._gb_text_to_record(out.read_text(), cache=False)):
            assert back.features[0].qualifiers["primer_seq"] == [oligo]

    def test_non_ascii_is_transliterated(self, plasmid, tmp_path):
        """A `®` in a feature note is not hypothetical — the `.dna` files in
        `tests/` carry one. GenBank is 7-bit ASCII."""
        plasmid.features[2].qualifiers["note"] = ["grow at 37 °C ± 2, β-lactam"]
        plasmid.description = "Ünïcode — description"
        out = tmp_path / "p.gb"
        sc._export_genbank_to_path(plasmid, out)
        text = out.read_text()
        assert _gb_problems(text) == []
        back = sc.load_genbank(str(out))
        note = next(f for f in back.features
                    if f.qualifiers.get("label") == ["Prev"])
        assert note.qualifiers["note"] == ["grow at 37 deg C +/- 2, beta-lactam"]
        assert back.description.startswith("Unicode - description")

    def test_non_ascii_display_name_is_transliterated(self, plasmid,
                                                      tmp_path):
        """The display name reaches the file through the `SpliceCraft-name:`
        COMMENT marker, which is stamped AFTER the normalisation pass — so
        folding the annotations alone left this one route open."""
        plasmid._tui_display_name = "pCAMBIA beta-lactamase 37 \u00b0C \u00b5L"
        out = tmp_path / "p.gb"
        sc._export_genbank_to_path(plasmid, out)
        text = out.read_text()
        assert text.isascii()
        assert _gb_problems(text) == []
        # and the caller's record still holds what the user typed
        assert "\u00b0" in plasmid._tui_display_name

    def test_storage_keeps_unicode(self, plasmid):
        """...but the LIBRARY keeps what the user typed. Transliteration is
        an export-time concession to the file format, not a data edit."""
        plasmid.features[2].qualifiers["note"] = ["grow at 37 °C"]
        assert "°C" in sc._record_to_gb_text(plasmid)

    def test_tabs_folded_on_export(self, plasmid, tmp_path):
        plasmid.features[2].qualifiers["note"] = ["a\tb"]
        out = tmp_path / "p.gb"
        sc._export_genbank_to_path(plasmid, out)
        assert "\t" not in out.read_text()

    def test_over_long_feature_key_keeps_the_location_column(self, tmp_path):
        """An INSDC feature key is <= 15 characters. Biopython's 21-column
        template is truncated rather than extended, so a longer key eats the
        space before the location and the feature line stops being splittable."""
        rec = SeqRecord(Seq("ATGC" * 50), id="pK", name="pK", description="d")
        rec.annotations.update(molecule_type="DNA", topology="linear")
        rec.features = [SeqFeature(FeatureLocation(10, 50, strand=1),
                                   type="an_extremely_long_feature_key",
                                   qualifiers={"label": ["k"]})]
        out = tmp_path / "k.gb"
        sc._export_genbank_to_path(rec, out)
        assert _gb_problems(out.read_text()) == []

    def test_topology_survives_a_non_column_aligned_locus(self, tmp_path):
        """ApE and older SnapGene put `circular` in a slightly different
        column, and Biopython — which reads the field POSITIONALLY — then
        reports no topology at all. The export path defaults a topology-less
        record to linear, so a circular plasmid used to come out the far side
        as a linear fragment."""
        src = tmp_path / "ape.gb"
        src.write_text(
            "LOCUS       pApe       400 bp    ds-DNA   circular     "
            "15-APR-2024\n"
            "FEATURES             Location/Qualifiers\n"
            "     misc_feature    11..100\n"
            '                     /label="x"\n'
            "ORIGIN\n"
            + "".join(
                f"{i * 60 + 1:>9} " + " ".join(
                    ("ATGCGTACGT" * 6)[j:j + 10] for j in range(0, 60, 10)
                ) + "\n" for i in range(400 // 60)
            )
            + f"{361:>9} " + " ".join(("ATGCGTACGT" * 4)[j:j + 10]
                                      for j in range(0, 40, 10)) + "\n"
            + "//\n"
        )
        rec = sc.load_genbank(str(src))
        assert rec.annotations.get("topology") == "circular"
        out = tmp_path / "out.gb"
        sc._export_genbank_to_path(rec, out)
        assert sc._topology_from_gb_text(out.read_text()) == "circular"

    @pytest.mark.parametrize("prefix,encoding", [
        ("﻿", "utf-8"),        # BOM — every Windows "save as UTF-8"
        ("", "cp1252"),             # a Latin-1 author name in a REFERENCE
    ])
    def test_import_tolerates_real_world_encodings(self, plasmid, tmp_path,
                                                   prefix, encoding):
        text = prefix + sc._record_to_gb_text(plasmid).replace(
            "conformance fixture", "conformance fixture (Müller)")
        p = tmp_path / "enc.gb"
        p.write_bytes(text.encode(encoding))
        rec = sc.load_genbank(str(p))
        assert len(rec.seq) == len(plasmid.seq)

    def test_import_tolerates_cr_only_line_endings(self, plasmid, tmp_path):
        p = tmp_path / "cr.gb"
        p.write_bytes(sc._record_to_gb_text(plasmid)
                      .replace("\n", "\r").encode("utf-8"))
        assert len(sc.load_genbank(str(p)).seq) == len(plasmid.seq)

    def test_file_round_trip_matches_library_round_trip(self, plasmid,
                                                        tmp_path):
        """Opening a `.gb` FILE used to skip the arrowless decode and the
        primer-seq repair that the library path applies, so the same record
        came back different depending on which door it walked through."""
        out = tmp_path / "p.gb"
        sc._export_genbank_to_path(plasmid, out)
        text = out.read_text()
        by_file = sc.load_genbank(str(out))
        by_text = sc._gb_text_to_record(text, cache=False)
        def sig(rec):
            return [(f.type, str(f.location), sorted(f.qualifiers.items()))
                    for f in rec.features]
        assert sig(by_file) == sig(by_text)
        arrow = next(f for f in by_file.features
                     if f.qualifiers.get("label") == ["arrowless"])
        assert arrow.location.strand == 0
        assert "SpliceCraft_strand" not in arrow.qualifiers

    def test_line_endings_are_lf_only(self, plasmid, tmp_path):
        out = tmp_path / "p.gb"
        sc._export_genbank_to_path(plasmid, out)
        assert b"\r" not in out.read_bytes()


# ── EMBL ──────────────────────────────────────────────────────────────────────

class TestEmblConformance:

    def test_id_line_has_every_field(self, plasmid, tmp_path):
        """ENA's ID line is seven semicolon-separated fields. Biopython leaves
        the data class empty always, and the sequence version empty unless the
        record id happens to carry a `.N` suffix."""
        out = tmp_path / "p.embl"
        sc._export_embl_to_path(plasmid, out)
        first = out.read_text().split("\n")[0]
        assert first.startswith("ID   ")
        fields = [f.strip() for f in first[5:].split(";")]
        assert len(fields) == 7, first
        assert all(fields), f"empty field in {first!r}"
        assert fields[1].startswith("SV ")
        assert fields[6].endswith("BP.")

    def test_multiline_note_stays_parseable(self, plasmid, tmp_path):
        """A newline inside a qualifier value made the EMBL writer emit a
        continuation flush against column 0, which Biopython's own parser then
        REJECTED — so the export produced a file nothing could open."""
        from Bio import SeqIO
        plasmid.features[2].qualifiers["note"] = ["line one\nline two"]
        out = tmp_path / "p.embl"
        sc._export_embl_to_path(plasmid, out)
        back = SeqIO.read(str(out), "embl")
        assert len(back.features) >= len(plasmid.features)

    def test_arrowless_survives(self, plasmid, tmp_path):
        out = tmp_path / "p.embl"
        sc._export_embl_to_path(plasmid, out)
        back = sc.load_genbank(str(out))
        arrow = next(f for f in back.features
                     if f.qualifiers.get("label") == ["arrowless"])
        assert arrow.location.strand == 0

    def test_accession_with_a_space_does_not_raise(self, plasmid, tmp_path):
        """EMBL forbids a space in the accession and Biopython raises from
        deep inside the writer. The LOCUS-safe name is always legal."""
        plasmid.id = "my plasmid"
        out = tmp_path / "p.embl"
        sc._export_embl_to_path(plasmid, out)
        assert out.exists()

    def test_lines_fit_eighty_columns(self, plasmid, tmp_path):
        out = tmp_path / "p.embl"
        sc._export_embl_to_path(plasmid, out)
        text = out.read_text()
        assert max(len(ln) for ln in text.split("\n")) <= 80
        assert text.isascii()


# ── GFF3 ──────────────────────────────────────────────────────────────────────

class TestGff3Conformance:

    def test_export_is_spec_clean(self, plasmid):
        assert _gff3_problems(sc._record_to_gff3(plasmid)) == []

    def test_no_duplicate_note_tags(self, plasmid):
        """GFF3 attribute tags are unique per line; multiple values are
        comma-separated. A feature with two unmapped qualifiers — i.e. nearly
        every real CDS — used to emit `Note=` twice, and parsers respond by
        keeping one and dropping the rest."""
        text = sc._record_to_gff3(plasmid)
        for ln in text.split("\n"):
            if not ln or ln.startswith("#"):
                continue
            tags = [kv.split("=", 1)[0] for kv in ln.split("\t")[8].split(";")]
            assert len(tags) == len(set(tags)), ln

    def test_wrap_feature_keeps_its_arc_order(self, plasmid):
        """A GenBank `join(1901..2000,1..50)` parses to [tail, head]. Exporting
        it head-first and re-importing SWAPPED the two arcs, which for a CDS is
        a different protein."""
        text = sc._record_to_gff3(plasmid)
        rows = [ln.split("\t") for ln in text.split("\n")
                if ln and not ln.startswith("#")]
        wrap = [r for r in rows if "Name=wrap" in r[8]]
        assert len(wrap) == 2
        assert (int(wrap[0][3]), int(wrap[0][4])) == (1901, 2000)
        assert (int(wrap[1][3]), int(wrap[1][4])) == (1, 50)

    def test_split_cds_phases_accumulate(self):
        """GFF3 phase is per-row: how many bases to drop from THIS row to hit
        a codon boundary. Repeating the first row's phase told every
        downstream translator to start each exon in frame 0."""
        rec = SeqRecord(Seq("ATGC" * 100), id="pC", name="pC", description="d")
        rec.annotations.update(molecule_type="DNA", topology="linear")
        rec.features = [SeqFeature(
            CompoundLocation([FeatureLocation(0, 100, strand=1),
                              FeatureLocation(200, 300, strand=1)]),
            type="CDS", qualifiers={"label": ["split"]})]
        rows = [ln.split("\t") for ln in sc._record_to_gff3(rec).split("\n")
                if ln and not ln.startswith("#")]
        cds = [r for r in rows if r[2] == "CDS"]
        assert [r[7] for r in cds] == ["0", "2"]

    def test_linear_split_feature_is_not_read_as_a_wrap(self):
        """The same two-part shape means different things on a circle and on
        a line: on a circle it wraps the origin and belongs tail-first, on a
        line it is exon 1 then exon 2 and reversing it hands every reader a
        cyclically rotated protein. `_feat_bounds` has always drawn that
        distinction; the GFF3 writer did not."""
        rec = SeqRecord(Seq("ATGC" * 150), id="pL", name="pL", description="d")
        rec.annotations.update(molecule_type="DNA", topology="linear")
        rec.features = [SeqFeature(
            CompoundLocation([FeatureLocation(0, 50, strand=1),
                              FeatureLocation(500, 600, strand=1)]),
            type="CDS", qualifiers={"label": ["spliced"]})]
        rows = [ln.split("\t") for ln in sc._record_to_gff3(rec).split("\n")
                if ln and not ln.startswith("#") and "\tCDS\t" in ln]
        assert [(int(r[3]), int(r[4])) for r in rows] == [(1, 50), (501, 600)]

    def test_circular_wrap_is_tail_first_whichever_way_it_arrives(self):
        """A GenBank `join(501..600,1..50)` parses to [tail, head]; a
        programmatically-built one is usually [head, tail]. Both are the same
        feature and must export the same way."""
        for parts in ([FeatureLocation(500, 600, strand=1),
                       FeatureLocation(0, 50, strand=1)],
                      [FeatureLocation(0, 50, strand=1),
                       FeatureLocation(500, 600, strand=1)]):
            rec = SeqRecord(Seq("ATGC" * 150), id="pW", name="pW",
                            description="d")
            rec.annotations.update(molecule_type="DNA", topology="circular")
            rec.features = [SeqFeature(CompoundLocation(parts), type="CDS",
                                       qualifiers={"label": ["w"]})]
            rows = [ln.split("\t")
                    for ln in sc._record_to_gff3(rec).split("\n")
                    if ln and not ln.startswith("#") and "\tCDS\t" in ln]
            assert [(int(r[3]), int(r[4])) for r in rows] == [(501, 600),
                                                              (1, 50)]
            assert [r[7] for r in rows] == ["0", "2"]

    def test_source_qualifiers_reach_the_region_row(self, plasmid):
        """The `source` feature is folded into the synthetic `region` row, so
        without carrying its qualifiers across they are simply dropped."""
        region = next(ln for ln in sc._record_to_gff3(plasmid).split("\n")
                      if "\tregion\t" in ln)
        assert "organism=synthetic%20construct" in region
        assert "Is_circular=true" in region

    @pytest.mark.parametrize("value", ["alpha,beta", "k=v", "a;b", "100%",
                                       "line\nbreak", "tab\there"])
    def test_qualifier_values_round_trip(self, value):
        rec = SeqRecord(Seq("ATGC" * 100), id="pQ", name="pQ", description="d")
        rec.annotations.update(molecule_type="DNA", topology="circular")
        rec.features = [SeqFeature(FeatureLocation(10, 100, strand=1),
                                   type="misc_feature",
                                   qualifiers={"label": ["x"],
                                               "note": [value]})]
        text = sc._record_to_gff3(rec)
        assert _gff3_problems(text) == []
        parsed = sc._parse_gff3_text(text)
        back = sc._gff3_features_to_biopython(parsed, len(rec.seq))
        assert back[0].qualifiers["note"] == [value]


# ── FASTA ─────────────────────────────────────────────────────────────────────

class TestFastaConformance:

    def test_sequence_is_wrapped(self, tmp_path):
        out = tmp_path / "s.fa"
        sc._export_fasta_to_path("pTest", "ACGT" * 500, out)
        lines = out.read_text().rstrip("\n").split("\n")
        assert lines[0] == ">pTest"
        assert all(len(ln) <= sc._FASTA_LINE_WIDTH for ln in lines[1:])
        assert "".join(lines[1:]) == "ACGT" * 500

    def test_newline_in_the_name_cannot_split_the_record(self, tmp_path):
        """Everything up to the first whitespace is the record ID and a
        newline ENDS the header — so a name carrying one used to prepend its
        own tail to the exported bases."""
        out = tmp_path / "s.fa"
        sc._export_fasta_to_path("my\nplasmid", "ACGTACGT", out)
        text = out.read_text()
        assert text.count(">") == 1
        assert text.split("\n")[0] == ">my plasmid"
        assert text.split("\n")[1] == "ACGTACGT"

    def test_whitespace_in_the_sequence_is_stripped(self, tmp_path):
        out = tmp_path / "s.fa"
        sc._export_fasta_to_path("p", "ACGT ACGT\nACGT", out)
        assert out.read_text() == ">p\nACGTACGTACGT\n"

    def test_leading_chevron_is_not_doubled(self, tmp_path):
        out = tmp_path / "s.fa"
        sc._export_fasta_to_path(">already", "ACGT", out)
        assert out.read_text().startswith(">already\n")

    def test_round_trips_through_the_reader(self, tmp_path):
        out = tmp_path / "s.fa"
        seq = "ACGTTGCA" * 137
        sc._export_fasta_to_path("pRT", seq, out)
        rid, back = sc._parse_fasta_single(str(out))
        assert (rid, back) == ("pRT", seq)


# ── CSV (RFC 4180) ────────────────────────────────────────────────────────────

class TestPrimerCsvConformance:

    def test_row_terminators_are_not_doubled(self, tmp_path):
        """`csv.writer` already emits the RFC 4180 CRLF. Letting Python's
        text layer translate the LF half again produced `\\r\\r\\n`, which
        spreadsheet importers read as a blank row between every primer."""
        out = tmp_path / "p.csv"
        sc._export_primers_to_csv(
            [{"name": "p1", "sequence": "ACGTACGTACGT", "tm": 58.0},
             {"name": "p2, comma", "sequence": "TTTTGGGGCCCC"}],
            out,
        )
        raw = out.read_bytes()
        assert b"\r\r\n" not in raw
        assert raw.count(b"\r\n") == raw.count(b"\n")

    def test_round_trips_through_the_importer(self, tmp_path):
        out = tmp_path / "p.csv"
        sc._export_primers_to_csv(
            [{"name": 'quote " and, comma', "sequence": "ACGTACGTACGT"}], out)
        res = sc._import_primers_from_csv(str(out))
        assert res["primers"][0]["name"] == 'quote " and, comma'
        assert res["primers"][0]["sequence"] == "ACGTACGTACGT"


# ── 2026-09-17 (second pass): real third-party files as the oracle ───────────
#
# The first pass wrote strict validators from the specs and pointed them at
# records SpliceCraft itself had built. This pass instead pushed 27 real files
# — NCBI nucleotide + protein, ENA EMBL, Biopython's own awkward GenBank
# fixtures, our five `.dna` fixtures — through every writer, and judged the
# output with parsers SpliceCraft had never been tested against (`gffutils`,
# `snapgene_reader`, a newer Biopython). Ten of the 31 records could not be
# exported to GenBank AT ALL. Every finding below is a real file, not a
# hypothetical. [INV-196]

_FIXTURE_DNA = sorted(pathlib.Path(__file__).parent.glob("*.dna"))


def _protein_record(seq="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESF"):
    """A protein record shaped the way NCBI serves one: `aa` units, a
    whole-length `source` feature that carries NO strand (a protein flat
    file has no strand column), and a `DBSOURCE` line."""
    rec = SeqRecord(Seq(seq), id="NP_TEST", name="NP_TEST",
                    description="protein conformance fixture")
    rec.annotations.update(molecule_type="protein", topology="linear",
                           data_file_division="PRI", date="17-SEP-2026")
    rec.features = [
        SeqFeature(FeatureLocation(0, len(seq)), type="source",
                   qualifiers={"organism": ["Homo sapiens"]}),
        SeqFeature(FeatureLocation(0, 10), type="Site",
                   qualifiers={"site_type": ["acetylation"]}),
    ]
    return rec


class TestNumericQualifierValues:
    """Regression guard for the 2026-09-17 fix: a qualifier value that is an
    `int` blocked GenBank export outright."""

    def test_int_codon_start_still_exports(self, tmp_path):
        """BioPython's `.dna` parser hands back `codon_start=1` as an INT.
        GenBank text has no numeric qualifier value, so the re-parse always
        yields `'1'` — and the exporter's round-trip guard compared `(1,)`
        with `('1',)` and refused the file. Every `.dna` a user opened was
        therefore un-exportable as GenBank."""
        rec = SeqRecord(Seq("ATGGCC" * 40), id="T", name="T",
                        description="d",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular",
                                     "data_file_division": "SYN",
                                     "date": "17-SEP-2026"})
        rec.features = [SeqFeature(FeatureLocation(0, 240, strand=1),
                                   type="CDS",
                                   qualifiers={"codon_start": [1],
                                               "transl_table": [11]})]
        out = tmp_path / "int.gb"
        sc._export_genbank_to_path(rec, out)        # used to raise ValueError
        text = out.read_text()
        assert "/codon_start=1" in text             # INSDC: unquoted integer
        assert "/transl_table=11" in text
        assert _gb_problems(text) == []

    def test_int_on_a_free_text_qualifier_is_quoted(self):
        """The quieter half of the same bug: BioPython writes a bare int
        UNQUOTED whatever the key, so an int on `/note` emitted `/note=1`.
        INSDC requires free text to be quoted."""
        rec = SeqRecord(Seq("ATGGCC" * 10), id="T", name="T", description="d",
                        annotations={"molecule_type": "DNA",
                                     "topology": "linear",
                                     "data_file_division": "SYN",
                                     "date": "17-SEP-2026"})
        rec.features = [SeqFeature(FeatureLocation(0, 60, strand=1),
                                   type="misc_feature",
                                   qualifiers={"note": [7]})]
        text = sc._record_to_gb_text(sc._normalize_for_genbank(rec))
        assert '/note="7"' in text
        assert "/note=7" not in text

    @pytest.mark.parametrize("path", _FIXTURE_DNA,
                             ids=[p.stem for p in _FIXTURE_DNA])
    def test_real_dna_fixture_round_trips_to_genbank(self, path, tmp_path):
        """End-to-end on the real files: open each committed `.dna` and
        export it as GenBank. All five raised before the fix."""
        rec = sc.load_genbank(str(path))
        out = tmp_path / "from_dna.gb"
        sc._export_genbank_to_path(rec, out)
        text = out.read_text()
        assert _gb_problems(text) == []
        back = sc._gb_text_to_record(text, cache=False)
        assert str(back.seq).upper() == str(rec.seq).upper()
        assert len(back.features) == len(rec.features)


class TestProteinRecords:
    """Regression guard for the 2026-09-17 fix: no protein record could be
    exported as GenBank."""

    def test_protein_record_exports(self, tmp_path):
        """A protein flat file has no strand column: BioPython writes
        `source 1..N` and reads it back STRANDLESS. `_normalize_for_genbank`
        stamped the nucleotide +1 on it regardless, so the round-trip guard
        reported a strand divergence and refused every `.gp` NCBI serves."""
        out = tmp_path / "p.gp"
        sc._export_genbank_to_path(_protein_record(), out)
        text = out.read_text()
        assert _gb_problems(text) == []
        assert " aa " in text.split("\n")[0]

    def test_protein_source_strand_matches_what_the_reader_returns(self):
        """The measured rule this fix encodes, pinned so a refactor cannot
        silently reintroduce the nucleotide assumption: on a protein record
        the normaliser must leave `source` strandless, because that is what
        the re-parse produces."""
        norm = sc._normalize_for_genbank(_protein_record())
        src = next(f for f in norm.features if f.type == "source")
        assert src.location.strand is None

        back = sc._gb_text_to_record(sc._record_to_gb_text(norm), cache=False)
        src_back = next(f for f in back.features if f.type == "source")
        assert src_back.location.strand == src.location.strand

    def test_nucleotide_source_still_gets_a_strand(self, plasmid):
        """The other half of the same branch — unchanged for DNA, where the
        reader DOES return +1."""
        norm = sc._normalize_for_genbank(plasmid)
        src = next(f for f in norm.features if f.type == "source")
        assert src.location.strand == 1


class TestDbSourceWrapping:
    """Regression guard for the 2026-09-17 fix: `DBSOURCE` was emitted on one
    unwrapped line."""

    def test_dbsource_wraps_at_eighty_columns(self):
        """Every protein record NCBI serves puts its whole xref list in
        DBSOURCE — 15,577 characters and 252 wrapped lines for UniProt
        P69905. BioPython writes that tag with `_write_single_line`, which
        only warns, so we emitted all of it on a single line. It is the same
        shape as the `/translation` bug the first pass found, one field over."""
        rec = _protein_record()
        rec.annotations["db_source"] = "xrefs: " + ", ".join(
            f"ACC{i:05d}.1" for i in range(400))
        text = sc._record_to_gb_text(rec)
        assert max(len(ln) for ln in text.split("\n")) <= 80
        assert _gb_problems(text) == []
        assert text.count("\nDBSOURCE") == 1        # one tag, many continuations

    def test_dbsource_survives_the_round_trip(self):
        """Wrapping must not mutate the value — the library persists every
        record as GenBank text, so a lossy fold would corrupt stored data."""
        rec = _protein_record()
        original = "xrefs: " + ", ".join(f"ACC{i:05d}.1" for i in range(400))
        rec.annotations["db_source"] = original
        back = sc._gb_text_to_record(sc._record_to_gb_text(rec), cache=False)
        got = back.annotations.get("db_source") or ""
        if isinstance(got, list):
            got = got[0] if got else ""
        assert " ".join(got.split()) == " ".join(original.split())

    def test_locus_line_is_never_folded(self):
        """The reason this is a narrow tag allow-list and not "wrap every
        header": LOCUS is a FIXED-COLUMN line."""
        rec = _protein_record()
        rec.annotations["db_source"] = " ".join(f"seg{i:04d}" for i in range(90))
        first = sc._record_to_gb_text(rec).split("\n")[0]
        assert first.startswith("LOCUS       ")
        assert re.search(r"(?<= )\d+$", first[12:40])


class TestDnaFeaturesPacket:
    """Regression guard for the 2026-09-17 fix: a `.dna` whose only feature
    was `source` carried an EMPTY features packet."""

    @staticmethod
    def _features_payloads(data):
        return [payload for tb, _len, payload
                in sc._iter_commercialsaas_packets(data)
                if tb == sc._COMMERCIALSAAS_PACKET_FEATURES]

    def test_a_features_packet_is_never_empty(self, plasmid):
        """The general invariant. `<Features nextValidID="0"/>` with no
        children is legal XML and unreadable in practice: the widely used
        third-party `snapgene_reader` raises `KeyError: 'Feature'` and never
        reaches the sequence."""
        for rec in (plasmid, _protein_record()):
            for payload in self._features_payloads(
                    sc._write_commercialsaas_dna_bytes(rec)):
                assert b"<Feature " in payload or b"<Feature>" in payload

    def test_source_only_record_omits_the_packet(self):
        """NCBI pUC19 (L09137 — the accession in our own README) has exactly
        one feature, `source`, which the `.dna` schema does not store. It
        exported a file no other tool could open."""
        rec = SeqRecord(Seq("ATGC" * 500), id="pUC", name="pUC",
                        description="d",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec.features = [SeqFeature(FeatureLocation(0, 2000, strand=1),
                                   type="source",
                                   qualifiers={"organism": ["synthetic"]})]
        data = sc._write_commercialsaas_dna_bytes(rec)
        assert self._features_payloads(data) == []
        # ...and the file is still a complete, readable `.dna`.
        from Bio import SeqIO
        import io as _io
        parsed = SeqIO.read(_io.BytesIO(data), sc._BIOPYTHON_DNA_FMT)
        assert str(parsed.seq).upper() == "ATGC" * 500
        assert parsed.annotations.get("topology") == "circular"


class TestGff3SourceRow:
    """Regression guard for the 2026-09-17 fix: the `source` feature was lost
    on every GenBank → GFF3 → GenBank cycle."""

    def test_region_row_rebuilds_the_source_feature(self, plasmid, tmp_path):
        """`_record_to_gff3` deliberately folds `source` onto the `region`
        row (NCBI's own convention) so its /organism and /mol_type survive —
        and `_parse_gff3_text` then dropped that row outright. Two halves of
        one format disagreeing, which is the shape that keeps recurring."""
        text = sc._record_to_gff3(plasmid)
        text += f"##FASTA\n>{plasmid.id}\n{str(plasmid.seq)}\n"
        p = tmp_path / "rt.gff3"
        p.write_text(text)
        back = sc._gff3_path_to_record(str(p))
        src = next((f for f in back.features if f.type == "source"), None)
        assert src is not None
        assert src.qualifiers.get("organism") == ["synthetic construct"]
        assert src.qualifiers.get("mol_type") == ["other DNA"]
        assert (int(src.location.start), int(src.location.end)) == (0, 2000)
        assert len(back.features) == len(plasmid.features)

    def test_overlay_path_adds_no_source(self, plasmid, tmp_path):
        """The behaviour that must NOT change: grafting a full-length
        `source` onto an already-loaded plasmid would be a spurious
        annotation, so the canvas-overlay path still ignores the row."""
        p = tmp_path / "overlay.gff3"
        p.write_text(sc._record_to_gff3(plasmid))
        target = SeqRecord(Seq(str(plasmid.seq)), id="t", name="t",
                           annotations={"molecule_type": "DNA",
                                        "topology": "circular"})
        added = sc._gff3_apply_to_loaded_record(target, str(p))
        assert added == len(plasmid.features) - 1        # every one but source
        assert not any(f.type == "source" for f in target.features)

    def test_a_bare_region_row_does_not_invent_a_source(self, tmp_path):
        """A hand-written GFF3 whose region row carries only bookkeeping
        must not gain an empty `source`."""
        p = tmp_path / "bare.gff3"
        p.write_text(
            "##gff-version 3\n"
            "##sequence-region chr1 1 12\n"
            "chr1\t.\tregion\t1\t12\t.\t.\t.\tID=chr1\n"
            "chr1\t.\tgene\t2\t8\t.\t+\t.\tID=g1;Name=g\n"
            "##FASTA\n>chr1\nATGCATGCATGC\n")
        back = sc._gff3_path_to_record(str(p))
        assert [f.type for f in back.features] == ["gene"]


class TestNucleotideFastaMessage:

    def test_a_protein_fasta_says_it_is_a_protein(self, tmp_path):
        """Regression guard for 2026-09-17: handing this nucleotide-only
        path an NCBI `.faa` reported `Non-IUPAC characters in sequence:
        EFLPQ`, which is true and useless."""
        p = tmp_path / "prot.fasta"
        p.write_text(">sp|P69905|HBA_HUMAN\nMVHLTPEEKSAVTALWGKVNVDEVGGEALGRL\n")
        with pytest.raises(ValueError, match="(?i)protein sequence"):
            sc._parse_fasta_single(str(p))

    def test_a_nucleotide_fasta_still_loads(self, tmp_path):
        p = tmp_path / "dna.fasta"
        p.write_text(">plasmid\nATGCATGCATGCRYKMSWN\n")
        name, seq = sc._parse_fasta_single(str(p))
        assert name == "plasmid"
        assert seq == "ATGCATGCATGCRYKMSWN"
