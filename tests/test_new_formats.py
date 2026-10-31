"""Minimum-viable coverage for the format-import / bulk-export
helpers landed this session. Targets the highest-leverage failure
classes (round-trip integrity, size-cap enforcement, error messages)
without exhaustive matrix coverage."""

from __future__ import annotations

from pathlib import Path

import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    CompoundLocation, FeatureLocation, SeqFeature,
)
from Bio.SeqRecord import SeqRecord

import splicecraft as sc
import splicecraft_fileio as _fileio  # the format-I/O caps moved here (Phase D);
# the parsers read them in this sibling's namespace, so patch it, not `sc`.


def _make_circular_record(seq: str = "ACGTACGTACGTACGTACGTACGTACGT") -> SeqRecord:
    """Tiny circular record with one normal feature + one wrap feature."""
    rec = SeqRecord(Seq(seq), id="plasmid", name="plasmid")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    # Normal feature.
    rec.features.append(SeqFeature(
        FeatureLocation(4, 12, strand=1),
        type="CDS", qualifiers={"label": ["gene1"]},
    ))
    # Wrap feature: origin-spanning.
    rec.features.append(SeqFeature(
        CompoundLocation([
            FeatureLocation(24, len(seq), strand=1),
            FeatureLocation(0, 3,  strand=1),
        ]),
        type="misc_feature", qualifiers={"label": ["wrap_feat"]},
    ))
    return rec


# ── GFF3 round-trip ────────────────────────────────────────────────────────────

def test_gff3_roundtrip_preserves_wrap_feature(tmp_path):
    """Export a record with a wrap feature to GFF3, re-import via the
    `##FASTA` directive path, and assert the wrap feature comes back
    as a CompoundLocation. Closes pitfall #26's inverse round-trip
    guarantee.
    """
    rec = _make_circular_record()
    # Emit GFF3 text and append a ##FASTA section with the sequence.
    gff_text = sc._record_to_gff3(rec)
    gff_text += "\n##FASTA\n"
    gff_text += f">{rec.id}\n{str(rec.seq)}\n"
    p = tmp_path / "wrap.gff3"
    p.write_text(gff_text)

    parsed = sc._gff3_path_to_record(str(p))
    assert str(parsed.seq) == str(rec.seq)
    # Find the wrap feature by label.
    wrap = next(
        f for f in parsed.features
        if "wrap_feat" in (f.qualifiers.get("label") or [])
    )
    assert isinstance(wrap.location, CompoundLocation)
    assert len(wrap.location.parts) == 2


def test_gff3_no_fasta_rejects_standalone_import(tmp_path):
    gff = "##gff-version 3\nplasmid\tSpliceCraft\tregion\t1\t100\t.\t+\t.\tID=plasmid\n"
    p = tmp_path / "no-fasta.gff3"
    p.write_text(gff)
    with pytest.raises(ValueError, match="##FASTA"):
        sc._gff3_path_to_record(str(p))


def test_gff3_apply_to_loaded_rejects_length_mismatch(tmp_path):
    rec = _make_circular_record()
    other_text = (
        "##gff-version 3\n"
        "##sequence-region plasmid 1 99\n"
        "plasmid\tSpliceCraft\tregion\t1\t99\t.\t+\t.\tID=plasmid\n"
    )
    p = tmp_path / "mismatch.gff3"
    p.write_text(other_text)
    with pytest.raises(ValueError, match="length"):
        sc._gff3_apply_to_loaded_record(rec, str(p))


def test_gff3_apply_to_loaded_rejects_fasta_payload(tmp_path):
    """If the GFF3 file IS a full record, the apply-overlay path must
    refuse it — otherwise the imported sequence gets silently
    discarded."""
    rec = _make_circular_record()
    gff = sc._record_to_gff3(rec) + f"\n##FASTA\n>plasmid\n{str(rec.seq)}\n"
    p = tmp_path / "full.gff3"
    p.write_text(gff)
    with pytest.raises(ValueError, match="##FASTA"):
        sc._gff3_apply_to_loaded_record(rec, str(p))


# ── EMBL export ────────────────────────────────────────────────────────────────

def test_embl_export_roundtrip(tmp_path):
    rec = _make_circular_record()
    p = tmp_path / "out.embl"
    summary = sc._export_embl_to_path(rec, str(p))
    assert summary["bp"] == len(rec.seq)
    # SeqIO.read should accept the emitted file.
    back = SeqIO.read(str(p), "embl")
    assert str(back.seq) == str(rec.seq)


def test_embl_export_rejects_empty():
    with pytest.raises(ValueError):
        sc._export_embl_to_path(
            SeqRecord(Seq("")), "/tmp/should-not-be-written.embl",
        )


# ── AB1 / FASTQ caps ───────────────────────────────────────────────────────────

def test_ab1_path_rejects_oversize(tmp_path, monkeypatch):
    """`_ab1_path_to_record` rejects files over `_BULK_IMPORT_MAX_BYTES`
    BEFORE any SeqIO parse. Pre-fix a 1 GB hostile .ab1 OOM'd the
    worker because SeqIO.parse streams unbounded."""
    p = tmp_path / "huge.ab1"
    p.write_bytes(b"\x00" * 100)
    monkeypatch.setattr(_fileio, "_BULK_IMPORT_MAX_BYTES", 50)
    with pytest.raises(ValueError):
        sc._ab1_path_to_record(str(p))


def _abif(entries, head_stride: int = 28, head_count: "int | None" = None) -> bytes:
    """A minimal ABIF trace: ``entries`` is ``[(name, number, elem_code,
    elem_num, data), ...]`` laid out the way an instrument writes them — a
    26-byte header after the marker, the data blocks, then the directory.

    ``head_stride`` / ``head_count`` override the two directory fields the
    header declares, which real instruments always write honestly and a
    hostile file need not.
    """
    import struct
    body = b""
    offsets = []
    base = 4 + 26
    for _n, _i, _c, _k, data in entries:
        offsets.append(base + len(body))
        body += data
    dir_off = base + len(body)
    directory = b"".join(
        struct.pack(">4sI2H4I", name.encode(), num, code, 1, k, len(data),
                    off if len(data) > 4 else
                    int.from_bytes(data.ljust(4, b"\0"), "big"), 0)
        for (name, num, code, k, data), off in zip(entries, offsets))
    head = struct.pack(">H4sI2H3I", 101, b"tdir", 1, 1023, head_stride,
                       len(entries) if head_count is None else head_count,
                       len(directory), dir_off)
    return b"ABIF" + head + body + directory


_BASES = b"ACGTACGTTGCA" * 4


def test_a_well_formed_synthetic_trace_loads(tmp_path):
    """The builder's premise: an honest trace parses."""
    p = tmp_path / "ok.ab1"
    p.write_bytes(_abif([("PBAS", 2, 2, len(_BASES), _BASES),
                         ("PCON", 2, 2, len(_BASES), bytes([30]) * len(_BASES))]))
    rec = sc._ab1_path_to_record(str(p))
    assert str(rec.seq) == _BASES.decode()
    assert rec.letter_annotations["phred_quality"] == [30] * len(_BASES)


@pytest.mark.parametrize("entries", [
    # an element count that disagrees with the data size: Biopython asserts
    [("PBAS", 2, 2, len(_BASES) + 7, _BASES)],
    # qualities stored as words, not characters: `.decode()` on a tuple
    [("PBAS", 2, 2, len(_BASES), _BASES),
     ("PCON", 2, 3, len(_BASES) // 2, bytes(len(_BASES)))],
], ids=["count-vs-size", "wrong-element-type"])
def test_a_corrupt_trace_is_a_value_error_not_a_crash(tmp_path, entries):
    """The Sanger picker, `load-file` and `read-trace` all catch ValueError;
    an AssertionError or AttributeError out of the parser exited the app
    (audit 2026-10-23)."""
    p = tmp_path / "bad.ab1"
    p.write_bytes(_abif(entries))
    with pytest.raises(ValueError):
        sc._ab1_path_to_record(str(p))


class TestAbifHeaderDirectoryBounds:
    """Biopython walks `range(head_elem_num)` seeking to
    `head_offset + index * head_elem_size`, trusting both header fields.
    Neither is bounded by the file, so a tiny trace can declare a directory
    no file could hold — and with a stride of ZERO every step re-reads the
    same entry, so the parse SUCCEEDS after doing the work. 217 bytes bought
    ~2-3 hours of it and then returned an ordinary 50 bp record: the size cap
    can't see a small file and no exception handler can catch a parse that
    never raises (audit 2026-10-24)."""

    def _entries(self):
        return [("PBAS", 2, 2, len(_BASES), _BASES),
                ("PCON", 2, 2, len(_BASES), bytes([30]) * len(_BASES))]

    def test_a_stride_shorter_than_one_entry_is_refused(self, tmp_path):
        p = tmp_path / "stride0.ab1"
        p.write_bytes(_abif(self._entries(), head_stride=0,
                            head_count=200_000))
        with pytest.raises(ValueError, match="(?i)corrupt ab1 header"):
            sc._ab1_path_to_record(str(p))

    def test_a_directory_that_cannot_fit_in_the_file_is_refused(self, tmp_path):
        p = tmp_path / "toomany.ab1"
        p.write_bytes(_abif(self._entries(), head_count=200_000))
        with pytest.raises(ValueError, match="(?i)corrupt ab1 header"):
            sc._ab1_path_to_record(str(p))

    def test_the_refusal_does_not_walk_the_claimed_directory(self, tmp_path):
        """The point of checking the header instead of catching the fallout:
        the cost must not scale with the number the file asked for. Pre-fix
        this spun for minutes on the way to answering."""
        import time
        p = tmp_path / "dos.ab1"
        p.write_bytes(_abif(self._entries(), head_stride=0,
                            head_count=100_000_000))
        t0 = time.process_time()
        with pytest.raises(ValueError, match="(?i)corrupt ab1 header"):
            sc._ab1_path_to_record(str(p))
        spent = time.process_time() - t0
        # CPU time: the pathology is a 100-million-entry header loop, not
        # I/O, and a wall-clock budget under `-n auto` fails on starvation.
        assert spent < 1.0, f"refusal took {spent:.3f} s CPU (budget 1.0)"

    def test_an_honest_header_still_parses(self, tmp_path):
        """The guard is arithmetic on the real layout, so the builder's own
        well-formed trace — and every real instrument file — is untouched."""
        p = tmp_path / "ok.ab1"
        p.write_bytes(_abif(self._entries()))
        rec = sc._ab1_path_to_record(str(p))
        assert str(rec.seq) == _BASES.decode()


def _fq(path, reads, offset):
    """Write ``reads`` — ``[(bases, [phred, ...]), ...]`` — with this offset."""
    with open(path, "w") as fh:
        for i, (bases, quals) in enumerate(reads):
            fh.write(f"@r{i}\n{bases}\n+\n"
                     + "".join(chr(q + offset) for q in quals) + "\n")
    return str(path)


class TestFastqQualityEncoding:
    """Every quality character at or above `@` with the highest past `J` is
    BOTH an old Illumina Phred+64 file and a modern Phred+33 one whose every
    base is >= Q31 (HiFi, duplex, consensus reads). Read as +64 that second
    kind lost 31 points a base — Q44 became Q13 — and the heterogeneity
    filter then threw away the real evidence (audit 2026-10-23). Phred+64 was
    only ever written by short-read Illumina pipelines, so a long read
    settles it."""

    def test_long_high_accuracy_reads_stay_phred33(self, tmp_path):
        import random
        rng = random.Random(5)
        reads = [("".join(rng.choice("ACGT") for _ in range(1200)),
                  [rng.randint(31, 60) for _ in range(1200)])
                 for _ in range(20)]
        p = _fq(tmp_path / "hifi.fastq", reads, 33)
        assert _fileio._fastq_quality_format(p) == "fastq"
        got = _fileio._fastq_path_to_records(p)
        assert [r.letter_annotations["phred_quality"] for r in got] \
            == [q for _b, q in reads]

    def test_old_illumina_short_reads_are_still_phred64(self, tmp_path):
        import random
        rng = random.Random(6)
        reads = [("".join(rng.choice("ACGT") for _ in range(76)),
                  [rng.choice([2] * 3 + list(range(10, 41))) for _ in range(76)])
                 for _ in range(50)]
        p = _fq(tmp_path / "ga2.fastq", reads, 64)
        assert _fileio._fastq_quality_format(p) == "fastq-illumina"
        got = _fileio._fastq_path_to_records(p)
        assert [r.letter_annotations["phred_quality"] for r in got] \
            == [q for _b, q in reads]

    def test_a_character_no_phred64_file_holds_means_phred33(self, tmp_path):
        # Q80 + 33 is `q`, past the `i` (Q41) no Phred+64 file goes beyond.
        p = _fq(tmp_path / "q80.fastq", [("ACGT" * 20, [80] * 80)], 33)
        assert _fileio._fastq_quality_format(p) == "fastq"

    def test_low_qualities_are_phred33_as_before(self, tmp_path):
        p = _fq(tmp_path / "ont.fastq", [("ACGT" * 20, [5] * 40 + [35] * 40)], 33)
        assert _fileio._fastq_quality_format(p) == "fastq"


def test_fastq_path_rejects_oversize(tmp_path, monkeypatch):
    p = tmp_path / "huge.fastq"
    p.write_text("@r\nACGT\n+\n!!!!\n")
    monkeypatch.setattr(_fileio, "_BULK_IMPORT_MAX_BYTES", 5)
    with pytest.raises(ValueError):
        sc._fastq_path_to_records(str(p))


def test_fastq_path_rejects_too_many_reads(tmp_path, monkeypatch):
    """Multi-read FASTQ above `_FASTQ_MAX_READS` is refused so a
    100k-read import doesn't blow the library JSON size budget."""
    body = "\n".join(
        f"@r{i}\nACGT\n+\n!!!!"
        for i in range(20)
    ) + "\n"
    p = tmp_path / "many.fastq"
    p.write_text(body)
    monkeypatch.setattr(_fileio, "_FASTQ_MAX_READS", 5)
    with pytest.raises(ValueError, match="reads"):
        sc._fastq_path_to_records(str(p))


# ── _safe_export_filename ──────────────────────────────────────────────────────

def test_safe_export_filename_windows_reserved():
    """Names matching CON/PRN/AUX/NUL/COM1-9/LPT1-9 get an underscore
    prefix so NTFS doesn't refuse the open."""
    assert sc._safe_export_filename("CON", "gb").startswith("_")
    assert sc._safe_export_filename("nul", "gb").startswith("_")
    assert sc._safe_export_filename("com3", "gb").startswith("_")


def test_safe_export_filename_strips_separators():
    """Path-traversal characters become `_` so we can't escape the
    target directory."""
    out = sc._safe_export_filename("../../../etc/passwd", "gb")
    assert "/" not in out
    assert "\\" not in out
    assert ".." not in out.split(".")[0]


def test_safe_export_filename_strips_control_chars():
    out = sc._safe_export_filename("foo\x00bar\x01baz", "gb")
    assert "\x00" not in out
    assert "\x01" not in out


# ── _bulk_export_collection ────────────────────────────────────────────────────

def _seed_collection(name: str, plasmid_names: list[str]) -> None:
    """Populate the collections cache with one collection containing
    plasmids whose names appear in `plasmid_names`. Each plasmid has
    valid gb_text so non-`.dna` export paths can parse it."""
    rec = _make_circular_record()
    gb_text = sc._record_to_gb_text(rec)
    entries = [{
        "id":      f"id_{n}",
        "name":    n,
        "size":    len(rec.seq),
        "n_feats": len(rec.features),
        "gb_text": gb_text,
    } for n in plasmid_names]
    sc._save_collections([{"name": name, "plasmids": entries}])


def test_bulk_export_collection_genbank_writes_gb_text_directly(tmp_path):
    """Fast path: GenBank export should write the entry's cached
    `gb_text` directly, skipping the BioPython round-trip."""
    _seed_collection("Backbones", ["pUC19", "pACYC184"])
    result = sc._bulk_export_collection("Backbones", "genbank", tmp_path)
    assert result["total"] == 2
    assert len(result["written"]) == 2
    assert not result["failures"]
    # Files exist + match the cached gb_text.
    for w in result["written"]:
        assert Path(w["path"]).exists()


def test_bulk_export_collection_unknown_format_raises():
    with pytest.raises(ValueError, match="unknown export format"):
        sc._bulk_export_collection("Nope", "xyz", "/tmp/wontwrite")


def test_bulk_export_collection_case_insensitive_collision(tmp_path):
    """Two plasmids whose names differ only in case must NOT
    overwrite each other on the target FS. Mirrors `_dna_sidecar_path`
    case-fold defence."""
    _seed_collection("Backbones", ["pUC19", "puc19"])
    result = sc._bulk_export_collection("Backbones", "genbank", tmp_path)
    assert len(result["written"]) == 2
    # Filenames should differ even on case-folding.
    casefolded = {Path(w["path"]).name.casefold() for w in result["written"]}
    assert len(casefolded) == 2


def test_bulk_export_collection_unknown_collection(tmp_path):
    sc._save_collections([])  # empty
    with pytest.raises(ValueError, match="collection not found"):
        sc._bulk_export_collection("Missing", "genbank", tmp_path)


# ── _backup_info damaged-row schema ───────────────────────────────────────────

def test_backup_info_damaged_file_returns_error_row(tmp_path):
    """Corrupt JSON should surface as `{n_entries: None, error: ...}`
    instead of being silently dropped — closes the "user can't see
    what's broken" gap in the Restore UI."""
    p = tmp_path / "library.json"
    p.write_text("{this is not valid json}")
    info = sc._backup_info(p)
    assert info is not None
    assert info["n_entries"] is None
    assert info["error"]


class TestGff3WholeMoleculeFeatureKeepsItsReadingOrder:
    """A feature whose parts tile the whole circle edge to edge has two legal
    reading orders — start at bp 0 and run forward, or start at the internal
    boundary and wrap through the origin — and they cover identical bases, so
    only the ROW ORDER tells them apart. `_record_to_gff3` writes rows in
    reading order, but the reader's minus-strand branch re-sorted them by
    coordinate: `complement(join(51..200,1..50))` came back as
    `complement(join(1..50,51..200))` and extracted a different sequence,
    while the GenBank path kept it (audit 2026-10-24). A spliced gene listed
    ascending the NCBI/Ensembl way always leaves intron gaps, so gating on
    the gapless shape keeps that reading untouched."""

    N = 200

    def _rec(self, parts, strand, ftype="misc_feature"):
        import random
        rng = random.Random(7)
        seq = "".join(rng.choice("ACGT") for _ in range(self.N))
        rec = SeqRecord(Seq(seq), id="pW", name="pW",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"})
        rec.features.append(SeqFeature(
            CompoundLocation([FeatureLocation(a, b, strand) for a, b in parts]),
            type=ftype, qualifiers={"label": ["W"]}))
        return rec

    @pytest.mark.parametrize("strand", [1, -1], ids=["plus", "minus"])
    @pytest.mark.parametrize("parts,label", [
        ([(0, 50), (50, 200)], "crossing"),
        ([(50, 200), (0, 50)], "whole-circle-from-zero"),
    ])
    def test_a_gapless_whole_circle_feature_round_trips(self, strand, parts,
                                                        label):
        rec = self._rec(parts, strand)
        want = str(rec.features[0].extract(rec.seq))
        back = sc._gff3_features_to_biopython(
            sc._parse_gff3_text(sc._record_to_gff3(rec)), self.N, circular=True)
        got = next(f for f in back if f.type == "misc_feature")
        assert str(got.extract(rec.seq)) == want, (label, got.location)

    def test_an_ascending_minus_strand_spliced_gene_is_still_read_backwards(self):
        """The reading this gate must not disturb: NCBI and Ensembl list a
        minus-strand gene's exons ASCENDING, and joining them in file order
        read every one of them exon-last (audit 2026-09-22). Those exons have
        intron gaps, so the whole-molecule shape never matches."""
        text = ("##gff-version 3\n"
                "s\tx\tregion\t1\t200\t.\t+\t.\tID=s;Is_circular=true\n"
                "s\tx\tCDS\t11\t40\t.\t-\t0\tID=g1\n"
                "s\tx\tCDS\t61\t90\t.\t-\t0\tID=g1\n"
                f"##FASTA\n>s\n{'ACGT' * 50}\n")
        feats = sc._gff3_features_to_biopython(
            sc._parse_gff3_text(text), 200, circular=True)
        loc = next(f for f in feats if f.type == "CDS").location
        assert [(int(p.start), int(p.end)) for p in loc.parts] == [
            (60, 90), (10, 40)]


class TestGff3OriginCrossingCdsKeepsItsFrame:
    """The writer gives each row its phase in READING order: the row the
    ribosome meets first has phase 0, and the one after the origin carries
    whatever the first left over. The reader picked the "5' row" by lowest
    start (+) / highest end (-) — the wrong row for a CDS that crosses the
    origin — and stamped its phase as `/codon_start`, so `MKPPKGPFPA*` came
    back as `ETAERPVPGL` (audit 2026-10-23)."""

    @pytest.mark.parametrize("strand", [1, -1])
    @pytest.mark.parametrize("leftover", [1, 2])
    def test_a_round_trip_keeps_the_reading_frame(self, strand, leftover):
        import random
        from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
        n = 300
        rng = random.Random(9 + leftover)
        seq = "".join(rng.choice("ACGT") for _ in range(n))
        first_len = 30 + leftover             # not a whole number of codons
        if strand == 1:
            parts = [FeatureLocation(n - first_len, n, 1),
                     FeatureLocation(0, 90 - first_len, 1)]
        else:
            parts = [FeatureLocation(0, first_len, -1),
                     FeatureLocation(n - (90 - first_len), n, -1)]
        rec = SeqRecord(Seq(seq), id="s", name="s",
                        annotations={"molecule_type": "DNA",
                                     "topology": "circular"},
                        features=[SeqFeature(CompoundLocation(parts),
                                             type="CDS")])
        want = str(rec.features[0].extract(rec.seq).translate())
        back = sc._gff3_features_to_biopython(
            sc._parse_gff3_text(sc._record_to_gff3(rec)), n, circular=True)
        cds = next(f for f in back if f.type == "CDS")
        assert "codon_start" not in cds.qualifiers, cds.qualifiers
        assert str(cds.extract(rec.seq).translate()) == want


class TestGff3NoteFreeText:
    """`Note` carries SpliceCraft's `key%3Dvalue` encoding for qualifiers
    with no GFF3 equivalent, AND free text from every other writer. A free
    note holding an escaped `=` became a qualifier named after the text
    before it — `induced with IPTG at a final OD600 ` — and the next library
    save wrote that out as a GenBank qualifier (audit 2026-10-23). The
    writer percent-encodes the whole key, so a key with a raw space is
    free text."""

    def _feats(self, attrs):
        n = 120
        text = ("##gff-version 3\n"
                f"s\tx\tmisc_feature\t11\t40\t.\t+\t.\tID=f1;{attrs}\n"
                f"##FASTA\n>s\n{'ACGT' * (n // 4)}\n")
        return sc._gff3_features_to_biopython(sc._parse_gff3_text(text), n)

    def test_an_escaped_equals_in_free_text_stays_in_the_note(self):
        f = self._feats("Note=induced with IPTG at a final OD600 %3D 0.6, "
                        "harvested at 4 h %3D good")[0]
        assert f.qualifiers.get("note") == [
            "induced with IPTG at a final OD600 = 0.6",
            " harvested at 4 h = good"], f.qualifiers
        assert set(f.qualifiers) <= {"note", "label"}, f.qualifiers

    def test_the_writers_own_encoding_still_round_trips(self):
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        rec = SeqRecord(Seq("ACGT" * 30), id="s", name="s",
                        annotations={"molecule_type": "DNA",
                                     "topology": "linear"})
        quals = {"Primer Tm": ["62.5 C"], "note": ["a=b, c", "plain"],
                 "ApEinfo_fwdcolor": ["#ff0000"]}
        rec.features.append(SeqFeature(FeatureLocation(10, 40, 1),
                                       type="misc_feature",
                                       qualifiers=dict(quals)))
        back = sc._gff3_features_to_biopython(
            sc._parse_gff3_text(sc._record_to_gff3(rec)), 120)
        got = next(f for f in back if f.type == "misc_feature").qualifiers
        for k, v in quals.items():
            assert got.get(k) == v, (k, got)

    def test_an_identifier_key_is_still_a_qualifier(self):
        f = self._feats("Note=codon_start%3D2")[0]
        assert f.qualifiers.get("codon_start") == ["2"]


class TestGff3FastaSectionIsValidated:
    """The `##FASTA` bases went straight into the record: spaces, tabs and
    position numbers became sequence letters (a 40 bp molecule imported as
    42 or 47, every feature misplaced), and a non-ASCII byte surfaced later
    as an opaque codec error (audit 2026-10-23)."""

    def _gff(self, fasta_lines):
        return ("##gff-version 3\n##sequence-region s 1 40\n"
                "s\tx\tmisc_feature\t5\t12\t.\t+\t.\tID=f1\n"
                "##FASTA\n>s\n" + "\n".join(fasta_lines) + "\n")

    def test_whitespace_and_position_numbers_are_not_bases(self):
        bases = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA"
        lines = ["1 acgttgcaac gttgcaacgt", "21\tTGCAACGTTG CAACGTTGCA "]
        parsed = sc._parse_gff3_text(self._gff(lines))
        assert parsed["fasta_seq"] == bases

    @pytest.mark.parametrize("junk", ["ACGTJACGT", "ACGTéACGT", "ACG#T"])
    def test_a_non_nucleotide_character_is_refused_by_name(self, junk):
        with pytest.raises(ValueError, match="(?i)iupac|nucleotide"):
            sc._parse_gff3_text(self._gff([junk]))


class TestGff3RoundTripFuzz:
    """Every feature shape, both strands, CDS and non-CDS — does it come back
    extracting the SAME bases? The GENBANK path is the co-oracle: a case
    GenBank keeps and GFF3 loses is a GFF3 fault, and a case both normalise is
    a convention the whole app shares. That is how the minus-strand
    whole-circle loss was found (audit 2026-10-24).

    GOTCHA that cost 400 phantom failures: build the cases with CANONICAL
    part orders. Biopython parses `complement(join(100..200,300..400))` into
    DESCENDING parts — reading order, 5' exon first — so an ascending-minus
    CompoundLocation is a shape no real file produces, and feeding them made
    the reader's correct NCBI-convention reading look like data loss."""

    N = 300

    def _shapes(self):
        n = self.N
        a, b, c = 63, 180, 234
        return [
            ("plain",              [(a, b)]),
            ("wrap gapped",        [(c, n), (0, a)]),
            ("whole-circle @0",    [(0, b), (b, n)]),
            ("whole-circle @mid",  [(b, n), (0, b)]),
            ("spliced 3 exon",     [(a, a + 9), (b, b + 9), (c, c + 9)]),
            ("spliced wrap",       [(c, n), (0, a), (a + 5, b)]),
            ("whole-circle 3part", [(0, a), (a, b), (b, n)]),
        ]

    def test_every_shape_round_trips(self):
        import random
        rng = random.Random(4242)
        seq = "".join(rng.choice("ACGT") for _ in range(self.N))
        losses = []
        for label, parts in self._shapes():
            for strand in (1, -1):
                for ftype in ("misc_feature", "CDS"):
                    # Reading order: the plus traversal, reversed on minus.
                    use = parts if strand == 1 else list(reversed(parts))
                    loc = (CompoundLocation(
                        [FeatureLocation(x, y, strand) for x, y in use])
                        if len(use) > 1
                        else FeatureLocation(use[0][0], use[0][1], strand))
                    rec = SeqRecord(Seq(seq), id="pF", name="pF",
                                    annotations={"molecule_type": "DNA",
                                                 "topology": "circular"})
                    rec.features.append(SeqFeature(
                        loc, type=ftype, qualifiers={"label": ["W"]}))
                    want = str(loc.extract(Seq(seq)))
                    back = sc._gff3_features_to_biopython(
                        sc._parse_gff3_text(sc._record_to_gff3(rec)),
                        self.N, circular=True)
                    got_f = next(f for f in back if f.type == ftype)
                    if str(got_f.extract(Seq(seq))) != want:
                        # Only a GFF3-ONLY loss counts: check GenBank too.
                        lib = sc._gb_text_to_record(
                            sc._record_to_gb_text(rec),
                            cache=False).features[0]
                        if str(lib.extract(Seq(seq))) == want:
                            losses.append((label, strand, ftype,
                                           str(loc), str(got_f.location)))
        assert losses == [], losses
