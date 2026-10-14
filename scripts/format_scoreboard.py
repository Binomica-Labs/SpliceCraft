#!/usr/bin/env python3
"""format_scoreboard — score sequence-file READERS on a real-world corpus.

Step 1 of owning the file-format layer: before any reader replaces
Biopython's, measure what the readers do to real files, so a replacement is
held to numbers rather than impressions. Two readers are scored, each on the
same bytes:

  biopython    Bio.SeqIO exactly as any caller gets it
  splicecraft  what SpliceCraft does with those bytes today: File > Open for
               a file, `_gb_text_to_record` for a library entry

A third reader is one function added to `READERS`.

    python scripts/format_scoreboard.py fetch             # public corpus
    python scripts/format_scoreboard.py run               # everything
    python scripts/format_scoreboard.py run --segments ncbi,synthetic
    python scripts/format_scoreboard.py run --local ~/Downloads --show 20

THE CORPUS lives outside the repo, in $SPLICECRAFT_CORPUS_DIR or
~/.cache/splicecraft-corpus. `fetch` puts the PUBLIC part there: NCBI records
(curated + a sample of vector records + one 4.6 Mb chromosome for speed) and
Biopython's own GenBank / EMBL / `.dna` test fixtures, which exist to be
awkward. PRIVATE files are read where they are and never copied: your library
(the content-addressed blobs), your `.dna` originals, and any folder passed
with --local or listed one per line in <corpus>/local_dirs.txt. Generated
segments need no files: `synthetic` (one canonical record re-written the ways
real tools write GenBank) and `fuzz` (seeded mutations of real files).

WHAT IS MEASURED, per segment and reader:
  accepted / refused / crashed  a refusal is a ValueError, the clean way to say
                                no; anything else escaping is a crash
  silent loss                   accepted, but the records disagree with the
                                file's own census (tests/_format_oracles.py):
                                bases, features, feature keys, qualifiers,
                                topology, record count — counted without a
                                parser, so a reader cannot agree with itself
  warnings                      what the reader admitted it had to guess
  speed                         MB/s over the segment, and per-file p50/p95
  write side                    each accepted record written back out by the
                                same tool: refusals, spec errors against the
                                independent validator, and round-trip drift
                                when the output is read back
  agreement                     both readers accepted and produced the same
                                sequence + features

SAFETY (CLAUDE.md data-dir rules): SpliceCraft is imported only after
$XDG_DATA_HOME, $SPLICECRAFT_DATA_DIR and $SPLICECRAFT_LOG point into a temp
dir, and the run refuses to start if that did not take. The real data dir is
located BEFORE that, and only ever read. Private items are reported by content
hash, never by name or path, unless --paths is given — and even then only on
the terminal; the saved JSON never carries them. Nothing here is imported by
the app; tests/test_format_scoreboard.py keeps the harness honest.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import random
import re
import signal
import statistics
import sys
import tempfile
import time
import urllib.parse
import urllib.request
import warnings
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tests._format_oracles import (  # noqa: E402  (stdlib-only, no app import)
    _gb_problems, dna_census, embl_census, genbank_census,
)

# The commercial editor's name and Biopython's format key for its `.dna`
# files, kept out of the source the way splicecraft_fileio keeps them.
_CS_NAME = bytes.fromhex("536e617047656e65").decode()
_DNA_FMT = _CS_NAME.lower()

GENBANK_EXTS = {".gb", ".gbk", ".genbank", ".gbff", ".gp", ".gpff", ".seq"}
EMBL_EXTS = {".embl"}
DNA_EXTS = {".dna"}
SEGMENTS = ("ncbi", "biopython", "local", "library", "dna", "synthetic", "fuzz")
PRIVATE_SEGMENTS = {"local", "library", "dna"}
PARSE_TIMEOUT_S = 30.0
_UA = "splicecraft-format-scoreboard (+https://github.com/Binomica-Labs/SpliceCraft)"


def corpus_dir() -> Path:
    raw = os.environ.get("SPLICECRAFT_CORPUS_DIR", "").strip()
    return Path(raw).expanduser() if raw else Path.home() / ".cache" / "splicecraft-corpus"


def real_data_dir() -> Path:
    """Where SpliceCraft keeps the user's data, resolved the way the app
    resolves it — BEFORE the sandbox moves the variables it reads."""
    override = os.environ.get("SPLICECRAFT_DATA_DIR")
    if override:
        return Path(override).expanduser()
    try:
        from platformdirs import user_data_dir
        return Path(user_data_dir("splicecraft", appauthor=False, roaming=False))
    except ImportError:
        return Path.home() / ".local" / "share" / "splicecraft"


# ── corpus items ──────────────────────────────────────────────────────────────

@dataclass
class Item:
    segment: str
    fmt: str                    # "genbank" | "embl" | "dna"
    path: Path                  # generated items are materialised to a temp file
    key: str                    # sha256[:12] of the bytes
    size: int
    private: bool
    label: str = ""             # public items only: file or case name
    producer: str = "unknown"
    kind: str = ""              # synthetic: "valid" | "malformed"

    def data(self) -> bytes:
        return self.path.read_bytes()


def fmt_for(path: Path) -> "str | None":
    ext = path.suffix.lower()
    if ext in GENBANK_EXTS:
        return "genbank"
    if ext in EMBL_EXTS:
        return "embl"
    if ext in DNA_EXTS:
        return "dna"
    return None


def decode_text(data: bytes) -> str:
    for enc in ("utf-8-sig", "cp1252"):
        try:
            return data.decode(enc)
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", "replace")


def detect_producer(text: str, fmt: str) -> str:
    """Best-effort guess at which tool wrote a file, from what it leaves in
    the text. Only labels the breakdown; nothing is scored on it."""
    if fmt == "dna":
        return "commercialsaas"
    head = text[:200_000]
    low = head.lower()
    if "splicecraft" in low:
        return "splicecraft"
    if re.search(r"^VERSION\s+\S+\.\d+", head, re.M) and "\nREFERENCE" in head:
        return "ncbi"
    if "twist bioscience" in low:
        return "twist"
    if _CS_NAME.lower() in low:
        return "commercialsaas"
    if "apeinfo_" in low:
        return "ape"
    if "benchling" in low:
        return "benchling"
    if "geneious" in low:
        return "geneious"
    if "plannotate" in low or "plasmidsaurus" in low:
        return "plasmidsaurus"
    if fmt == "embl" and re.search(r"^ID\s+\S+;\s*SV", head, re.M):
        return "ena"
    return "unknown"


def _make_item(segment: str, path: Path, *, private: bool, label: str = "",
               kind: str = "") -> "Item | None":
    fmt = fmt_for(path)
    if fmt is None:
        return None
    try:
        data = path.read_bytes()
    except OSError:
        return None
    if not data:
        return None
    return Item(segment=segment, fmt=fmt, path=path,
                key=hashlib.sha256(data).hexdigest()[:12], size=len(data),
                private=private, label="" if private else (label or path.name),
                producer=detect_producer(decode_text(data[:200_000]), fmt),
                kind=kind)


def _walk(root: Path) -> "list[Path]":
    if root.is_file():
        return [root]
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if not d.startswith(".") and d not in ("node_modules", "__pycache__")]
        out.extend(Path(dirpath) / f for f in filenames)
    return sorted(out)


def discover(segments, *, corpus: Path, data_dir: Path, local: "list[Path]",
             workdir: Path, sample: int, seed: int, fuzz_bases: int,
             fuzz_per_base: int) -> "list[Item]":
    items: "list[Item]" = []

    def add_tree(segment: str, root: Path, *, private: bool, only=None) -> None:
        if not root.exists():
            return
        seen = {i.key for i in items if i.segment == segment}
        for p in _walk(root):
            if only is not None:
                try:
                    with open(p, "rb") as fh:
                        if not only(fh.read(64)):
                            continue
                except OSError:
                    continue
            it = _make_item(segment, p, private=private)
            if it is not None and it.key not in seen:
                seen.add(it.key)
                items.append(it)

    if "ncbi" in segments:
        add_tree("ncbi", corpus / "ncbi", private=False)
    if "biopython" in segments:
        add_tree("biopython", corpus / "biopython", private=False)
    if "local" in segments:
        dirs = list(local)
        listing = corpus / "local_dirs.txt"
        if listing.is_file():
            dirs += [Path(ln.strip()).expanduser()
                     for ln in listing.read_text(encoding="utf-8").splitlines()
                     if ln.strip() and not ln.lstrip().startswith("#")]
        for d in dirs:
            add_tree("local", d, private=True)
    if "library" in segments:
        add_tree("library", data_dir / "plasmid_blobs", private=True,
                 only=lambda head: head.lstrip().startswith(b"LOCUS"))
    if "dna" in segments:
        add_tree("dna", data_dir / "dna_originals", private=True)
    if "synthetic" in segments:
        gen = workdir / "synthetic"
        gen.mkdir(parents=True, exist_ok=True)
        for name, kind, text, enc in synthetic_cases():
            p = gen / f"{name}.gb"
            p.write_bytes(text.encode(enc))
            it = _make_item("synthetic", p, private=False, label=name, kind=kind)
            if it is not None:
                items.append(it)
    if sample:
        rng = random.Random(seed)
        by_seg: "dict[str, list[Item]]" = {}
        for it in items:
            by_seg.setdefault(it.segment, []).append(it)
        items = []
        for seg_items in by_seg.values():
            items.extend(seg_items if len(seg_items) <= sample
                         else sorted(rng.sample(seg_items, sample), key=lambda i: str(i.path)))
    if "fuzz" in segments:
        items.extend(fuzz_items(items, workdir / "fuzz", seed=seed,
                                n_bases=fuzz_bases, per_base=fuzz_per_base))
    return items


# ── synthetic: one record, written the ways real tools write GenBank ────────

def _synthetic_seq(n: int = 1200, seed: int = 7) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("acgt") for _ in range(n))


def _origin_block(seq: str, *, width: int = 9, chunk: int = 10, per_line: int = 60,
                  numbers: bool = True) -> "list[str]":
    lines = []
    for i in range(0, len(seq), per_line):
        part = " ".join(seq[j:j + chunk] for j in range(i, min(i + per_line, len(seq)), chunk))
        lines.append(f"{i + 1:>{width}} {part}" if numbers else part)
    return lines


_SYN_FEATURES = [
    ("source", "1..1200", ['/organism="synthetic DNA construct"', '/mol_type="other DNA"']),
    ("promoter", "10..60", ['/label="Pdemo"']),
    ("CDS", "100..402", ['/label="GeneX"', "/codon_start=1",
                         '/note="a note long enough that the writer has to wrap',
                         'it onto a second line of the qualifier column"']),
    ("misc_feature", "complement(500..560)", ['/label="rev"']),
    ("misc_feature", "join(1150..1200,1..20)", ['/label="wrap"']),
    ("CDS", "join(600..650,700..760)", ['/label="spliced"']),
    ("primer_bind", "complement(900..920)", ['/label="P1"']),
    ("terminator", "1000..1040", ['/label="Tdemo"']),
]


def _syn_record(*, name: str = "pSynth", seq: "str | None" = None,
                locus: "str | None" = None, features=None,
                features_header: str = "FEATURES             Location/Qualifiers",
                key_col: int = 5, value_col: int = 21, origin=None,
                header_extra: "list[str] | None" = None, terminator: bool = True) -> str:
    seq = seq if seq is not None else _synthetic_seq()
    features = _SYN_FEATURES if features is None else features
    locus = locus if locus is not None else (
        f"LOCUS       {name:<16}{len(seq):>12} bp    DNA     circular SYN 01-JAN-2026")
    lines = [locus,
             "DEFINITION  Synthetic scoreboard fixture.",
             f"ACCESSION   {name}",
             f"VERSION     {name}",
             "KEYWORDS    .",
             "SOURCE      synthetic DNA construct",
             "  ORGANISM  synthetic DNA construct",
             *(header_extra or []),
             features_header]
    for key, loc, quals in features:
        lines.append(" " * key_col + key.ljust(value_col - key_col) + loc)
        lines.extend(" " * value_col + q for q in quals)
    lines.append("ORIGIN")
    lines.extend(origin if origin is not None else _origin_block(seq))
    if terminator:
        lines.append("//")
    return "\n".join(lines) + "\n"


def synthetic_cases() -> "list[tuple[str, str, str, str]]":
    """(name, kind, text, encoding). `valid` = well-formed, or a layout real
    tools write that a forgiving reader should take; `malformed` = genuinely
    broken, where a clean refusal is fine and a silent mis-read is not."""
    base = _syn_record()
    seq = _synthetic_seq()
    cases: "list[tuple[str, str, str, str]]" = []

    def add(name, kind, text, enc="utf-8"):
        cases.append((name, kind, text, enc))

    add("canonical", "valid", base)
    add("crlf", "valid", base.replace("\n", "\r\n"))
    add("cr_only", "valid", base.replace("\n", "\r"))
    add("utf8_bom", "valid", "﻿" + base)
    add("cp1252_note", "valid", base.replace("Pdemo", "Pµdemo"), "cp1252")
    add("utf8_note", "valid", base.replace("Pdemo", "Pµdemo °C"))
    add("trailing_whitespace", "valid", "\n".join(ln + "   " for ln in base.split("\n")))
    add("uppercase_sequence", "valid", _syn_record(origin=_origin_block(seq.upper())))
    add("spaced_locus_name", "valid", _syn_record(
        locus=f"LOCUS       GeneX in Demo Vector   {len(seq)} bp DNA     circular   01-JAN-2026"))
    add("features_header_spacing", "valid", _syn_record(
        features_header="FEATURES" + " " * 33 + "Location/Qualifiers"))
    add("vendor_column_42", "valid", _syn_record(
        locus=f"LOCUS       GeneX in Demo Vector                {len(seq)} bp DNA     circular   01-JAN-2026",
        features_header="FEATURES" + " " * 33 + "Location/Qualifiers", value_col=41))
    add("tab_indented_features", "valid", re.sub(
        r"(?m)^ {21}", "\t\t\t", re.sub(r"(?m)^ {5}(\S+) +", lambda m: "\t" + m.group(1) + "\t\t", base)))
    add("locus_seven_tokens", "valid", _syn_record(
        locus=f"LOCUS       pSynth {len(seq)} bp DNA circular 01-JAN-2026"))
    add("locus_name_only", "valid", _syn_record(locus="LOCUS       pSynth"))
    add("locus_no_division_or_date", "valid", _syn_record(
        locus=f"LOCUS       pSynth{len(seq):>22} bp    DNA     circular"))
    add("locus_topology_one_column_over", "valid", _syn_record(
        locus=f"LOCUS       pSynth{len(seq):>22} bp    ds-DNA   circular     01-JAN-2026"))
    add("missing_terminator", "valid", _syn_record(terminator=False))
    add("blank_line_in_features", "valid", base.replace(
        '     terminator', '\n     terminator', 1))
    add("base_count_line", "valid", base.replace(
        "\nORIGIN", "\nBASE COUNT      300 a    300 c    300 g    300 t\nORIGIN", 1))
    add("extra_header_keywords", "valid", _syn_record(header_extra=[
        "DBLINK      BioProject: PRJNA000000", "PROJECT     demo", "XYZFIELD    unknown keyword"]))
    add("trailing_underscore_key", "valid", base.replace("     misc_feature    ", "     misc_feature_   "))
    add("over_indented_key", "valid", base.replace("     promoter        ", "      promoter       ", 1))
    add("slash_line_in_quoted_value", "valid", base.replace(
        'it onto a second line of the qualifier column"',
        '/and a path-like continuation of the same note"'))
    add("qualifier_without_value", "valid", base.replace(
        '/label="rev"', '/label="rev"\n                     /pseudo'))
    add("empty_quoted_value", "valid", base.replace('/label="Tdemo"', '/label="Tdemo"\n                     /note=""'))
    add("escaped_quotes", "valid", base.replace('/label="P1"', '/label="P1"\n                     /note="he said ""hi"""'))
    add("duplicate_qualifiers", "valid", base.replace('/label="GeneX"', '/label="GeneX"\n                     /note="one"\n                     /note="two"'))
    add("long_qualifier_value", "valid", base.replace(
        '/label="Tdemo"', '/label="Tdemo"\n' + "\n".join(
            " " * 21 + (('/note="' if i == 0 else "") + ("word " * 11).strip()
                        + ('"' if i == 59 else ""))
            for i in range(60))))
    add("location_variants", "valid", _syn_record(features=_SYN_FEATURES + [
        ("misc_feature", "order(200..210,220..230)", ['/label="ordered"']),
        ("misc_feature", "<1..>50", ['/label="fuzzy"']),
        ("misc_feature", "300^301", ['/label="between"']),
        ("misc_feature", "450", ['/label="single"']),
        ("misc_feature", "join(complement(800..810),complement(780..790))", ['/label="revjoin"']),
    ]))
    add("wrapped_location", "valid", _syn_record(features=_SYN_FEATURES + [
        ("CDS", "join(100..110,120..130,140..150,160..170,180..190,200..210,",
         ["220..230,240..250)", '/label="many_exons"'])]))
    add("multi_record", "valid", base + _syn_record(name="pSynth2", seq=_synthetic_seq(seed=8)))
    add("origin_without_numbers", "malformed", _syn_record(
        origin=_origin_block(seq, numbers=False)))
    add("origin_one_long_line", "malformed", _syn_record(origin=[f"{1:>9} {seq}"]))
    add("origin_narrow_numbers", "malformed", _syn_record(origin=_origin_block(seq, width=3)))
    add("sequence_with_invalid_chars", "malformed", base.replace(seq[:10], "acgtxxacgt", 1))
    add("unterminated_quote", "malformed", base.replace('/label="Pdemo"', '/label="Pdemo'))
    add("truncated_mid_sequence", "malformed", base[:len(base) - 400])
    add("no_locus_line", "malformed", "\n".join(base.split("\n")[1:]))
    add("contig_only", "malformed", _syn_record(origin=[]).replace(
        "ORIGIN\n//", "CONTIG      join(AB000001.1:1..600,AB000002.1:1..600)\n//"))
    return cases


# ── fuzz: seeded mutations of real files ─────────────────────────────────────

_FUZZ_OPS = ("delete_line", "duplicate_line", "swap_lines", "shift_indent",
             "replace_char", "truncate", "tabify", "drop_terminator",
             "insert_blank", "delete_char", "lower_keyword")


def _mutate(text: str, rng: random.Random) -> "tuple[str, list[str]]":
    ops = []
    lines = text.split("\n")
    for _ in range(rng.randint(1, 3)):
        op = rng.choice(_FUZZ_OPS)
        ops.append(op)
        if not lines:
            break
        i = rng.randrange(len(lines))
        if op == "delete_line":
            del lines[i]
        elif op == "duplicate_line":
            lines.insert(i, lines[i])
        elif op == "swap_lines" and i + 1 < len(lines):
            lines[i], lines[i + 1] = lines[i + 1], lines[i]
        elif op == "shift_indent":
            n = rng.choice((-3, -2, -1, 1, 2, 3))
            lines[i] = (" " * n + lines[i]) if n > 0 else lines[i][min(-n, len(lines[i]) - len(lines[i].lstrip(" "))):]
        elif op == "replace_char" and lines[i]:
            j = rng.randrange(len(lines[i]))
            lines[i] = lines[i][:j] + chr(rng.randrange(33, 127)) + lines[i][j + 1:]
        elif op == "truncate":
            joined = "\n".join(lines)
            lines = joined[:rng.randrange(len(joined) + 1)].split("\n")
        elif op == "tabify":
            lines[i] = re.sub(r" {2,8}", "\t", lines[i], count=1)
        elif op == "drop_terminator":
            lines = [ln for ln in lines if ln.strip() != "//"]
        elif op == "insert_blank":
            lines.insert(i, "")
        elif op == "delete_char" and lines[i]:
            j = rng.randrange(len(lines[i]))
            lines[i] = lines[i][:j] + lines[i][j + 1:]
        elif op == "lower_keyword" and lines[i][:1].isalpha():
            word = lines[i].split()[0]
            lines[i] = word.lower() + lines[i][len(word):]
    return "\n".join(lines), ops


def fuzz_items(pool: "list[Item]", out: Path, *, seed: int, n_bases: int,
               per_base: int) -> "list[Item]":
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    bases = [i for i in pool if not i.private and i.fmt in ("genbank", "embl")
             and i.size < 400_000]
    bases = sorted(bases, key=lambda i: i.key)
    if len(bases) > n_bases:
        bases = rng.sample(bases, n_bases)
    items = []
    for b in bases:
        text = decode_text(b.data())
        for k in range(per_base):
            mutated, ops = _mutate(text, rng)
            p = out / f"{b.key}-{k}{b.path.suffix.lower() if b.fmt == 'genbank' else '.embl'}"
            p.write_bytes(mutated.encode("utf-8"))
            it = _make_item("fuzz", p, private=False, label=f"{b.label or b.key}:{'+'.join(ops)}")
            if it is not None:
                items.append(it)
    return items


# ── readers ───────────────────────────────────────────────────────────────────

class ParseTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise ParseTimeout(f"no answer within {PARSE_TIMEOUT_S:.0f} s")


def read_biopython(item: Item, sc) -> list:
    from Bio import SeqIO
    if item.fmt == "dna":
        return list(SeqIO.parse(str(item.path), _DNA_FMT))
    return list(SeqIO.parse(str(item.path), item.fmt))


def open_like_the_app(sc, path: Path, fmt: str) -> list:
    """File > Open: a file goes through `load_genbank`; a GenBank file that
    turns out to hold several records is parsed again as records and becomes
    a collection."""
    if fmt == "genbank":
        try:
            return [sc.load_genbank(str(path))]
        except sc._MultipleRecordsError:
            return sc._parse_genbank_file_records(str(path))
    return [sc.load_genbank(str(path))]


def read_splicecraft(item: Item, sc) -> list:
    if item.segment == "library":
        return [sc._gb_text_to_record(item.path.read_text(encoding="utf-8"), cache=False)]
    return open_like_the_app(sc, item.path, item.fmt)


READERS = {"biopython": read_biopython, "splicecraft": read_splicecraft}


def write_biopython(records, sc) -> "list[str]":
    from Bio import SeqIO
    out = []
    for rec in records:
        buf = io.StringIO()
        SeqIO.write([rec], buf, "genbank")
        out.append(buf.getvalue())
    return out


def write_splicecraft(records, sc, tmp: Path) -> "list[str]":
    out = []
    for rec in records:
        p = tmp / "export.gb"
        sc._export_genbank_to_path(rec, str(p))
        out.append(p.read_text(encoding="utf-8"))
    return out


# ── comparisons ───────────────────────────────────────────────────────────────

# Qualifiers a reader is DESIGNED to consume (SpliceCraft's arrowless marker
# becomes a strand); counting them would score a decode as a loss.
_CONSUMED_QUALIFIERS = {"SpliceCraft_strand"}


def _n_qualifiers(rec) -> int:
    n = 0
    for f in rec.features:
        for k, v in (f.qualifiers or {}).items():
            if k in _CONSUMED_QUALIFIERS:
                continue
            n += len(v) if isinstance(v, (list, tuple)) else 1
    return n


def _key_norm(k: str) -> str:
    return str(k).rstrip("_")


_VALID_KEY_RE = re.compile(r"[A-Za-z0-9_'*-]+")


def garbage_keys(records) -> int:
    """Features whose TYPE cannot be a feature key — spaces, `/`, `=` or a
    bare number in it. Needs no census: a reader that emits one has read a
    location, a qualifier or a sequence line as a feature."""
    return sum(1 for rec in records for f in rec.features
               if not _VALID_KEY_RE.fullmatch(str(f.type)) or str(f.type).isdigit())


def census_issues(records, census, fmt: str) -> "list[str]":
    """Where the parsed records disagree with what the file declares."""
    if fmt == "dna":
        return _dna_issues(records, census)
    if garbage_keys(records):
        return ["garbage_feature_key"]
    if len(records) != len(census):
        return ["record_count"]
    issues = []
    for rec, c in zip(records, census):
        n = len(rec.seq) if rec.seq is not None else 0
        if c["seq_len"] is not None and n != c["seq_len"]:
            issues.append("bases")
        if len(rec.features) != c["n_features"]:
            issues.append("feature_count")
        elif sorted(map(_key_norm, (f.type for f in rec.features))) != sorted(map(_key_norm, c["keys"])):
            issues.append("feature_keys")
        expected_q = sum(1 for q in c["qualifier_names"] if q not in _CONSUMED_QUALIFIERS)
        if _n_qualifiers(rec) != expected_q:
            issues.append("qualifiers")
        if c["topology"] and (rec.annotations or {}).get("topology") != c["topology"]:
            issues.append("topology")
    return sorted(set(issues))


def _strip_tags(s: str) -> str:
    return re.sub(r"<[^>]+>", "", s or "")


def _dna_issues(records, c) -> "list[str]":
    if len(records) != 1:
        return ["record_count"]
    rec = records[0]
    issues = []
    if c["seq_len"] is not None and len(rec.seq) != c["seq_len"]:
        issues.append("bases")
    # The feature packet may itself hold `primer_bind` features, and a
    # reader turns the PRIMER packet's binding sites into `primer_bind` too,
    # filtering some by design — so primer_bind is compared on neither side.
    # (`source` IS compared: a `.dna` may hold one, and no reader invents it.)
    feats = [f for f in rec.features if f.type != "primer_bind"]
    keep = [i for i, t in enumerate(c.get("feature_types") or []) if t != "primer_bind"]
    if len(feats) != len(keep):
        issues.append("feature_count")
    else:
        # A reader may keep the name as the label or, when the feature has a
        # label of its own, beside it — either way the name is still there.
        kept = [{str(v) for k in ("label", "name") for v in (f.qualifiers.get(k) or [])}
                for f in feats]
        names = [_strip_tags(c["feature_names"][i]) for i in keep]
        if any(n and n not in k for n, k in zip(names, kept)):
            issues.append("feature_names")
    if c["topology"] and (rec.annotations or {}).get("topology") != c["topology"]:
        issues.append("topology")
    return sorted(set(issues))


def dna_metadata_dropped(records, c) -> "list[str]":
    """What the format stores that a reader did not bring back. Not scored as
    loss (a reader may not model it) but reported, because a viewer needs it."""
    if len(records) != 1:
        return []
    rec = records[0]
    out = []
    if c["n_colored"]:
        colored = sum(1 for f in rec.features
                      if any("color" in k.lower() for k in (f.qualifiers or {})))
        if colored == 0 and not getattr(rec, "_dna_original_bytes", None):
            out.append("colours")
    if c["n_primers"]:
        seqs = sum(1 for f in rec.features if f.type == "primer_bind"
                   and any("seq" in k.lower() for k in (f.qualifiers or {})))
        if seqs == 0 and not getattr(rec, "_dna_primer_entries", None):
            out.append("primer_sequences")
    if c["has_history"] and not getattr(rec, "_dna_history_xml", None):
        out.append("history")
    return out


def _loc_sig(loc):
    """A location as data: operator, and per part the position KINDS (exact,
    before, after, between, …), coordinates, strand and remote reference.
    Strand None and 0 are one value — "no arrow" — as SpliceCraft defines it."""
    if loc is None:
        return None
    parts = []
    for p in getattr(loc, "parts", [loc]):
        ends = []
        for pos in (p.start, p.end):
            try:
                ends.append((type(pos).__name__, int(pos)))
            except (TypeError, ValueError):
                ends.append((type(pos).__name__, None))
        parts.append((tuple(ends), p.strand or 0, getattr(p, "ref", None)))
    return (getattr(loc, "operator", None), tuple(parts))


def signature(records) -> list:
    """What two readings must share to count as the same record. Qualifier
    values compare as text: GenBank has no types, so a `.dna` reader's int
    `1` and the `"1"` read back from GenBank are the same value."""
    sig = []
    for rec in records:
        feats = []
        for f in rec.features:
            quals = tuple(sorted((k, tuple(str(x) for x in (v if isinstance(v, (list, tuple)) else [v])))
                                 for k, v in (f.qualifiers or {}).items()
                                 if k not in _CONSUMED_QUALIFIERS))
            feats.append((_key_norm(f.type), _loc_sig(f.location), quals))
        try:
            bases = str(rec.seq).upper()
        except ValueError:                  # a CONTIG record: length, no bases
            bases = f"<undefined {len(rec.seq)}>"
        sig.append((bases, feats))
    return sig


def _first_diff(a: list, b: list) -> str:
    if len(a) != len(b):
        return "record_count"
    for (sa, fa), (sb, fb) in zip(a, b):
        if sa != sb:
            return "bases"
        if len(fa) != len(fb):
            return "feature_count"
        for x, y in zip(fa, fb):
            if x[0] != y[0]:
                return "feature_key"
            if x[1] != y[1]:
                return "location"
            if x[2] != y[2]:
                # Text the writer turned into 7-bit ASCII on purpose (INSDC
                # files are ASCII) is reported apart from a real change: the
                # SAME qualifier names on both sides, and every value that
                # changed had non-ASCII in it. An added or dropped qualifier
                # is never that.
                left, right = set(x[2]) - set(y[2]), set(y[2]) - set(x[2])
                if (left and {k for k, _ in left} == {k for k, _ in right}
                        and {k for k, _ in x[2]} == {k for k, _ in y[2]}
                        and all(any(not v.isascii() for v in vals) for _, vals in left)):
                    return "non_ascii_transliterated"
                return "qualifiers"
    return ""


def split_records(text: str) -> "list[str]":
    """Records of a flat file, each ending at its `//` line (blank lines
    between records dropped, a trailing unterminated record kept)."""
    parts, cur = [], []
    for ln in text.split("\n"):
        if not cur and not ln.strip():
            continue
        cur.append(ln)
        if ln.strip() == "//":
            parts.append("\n".join(cur) + "\n")
            cur = []
    if any(x.strip() for x in cur):
        parts.append("\n".join(cur))
    return parts


# ── one item, one reader ──────────────────────────────────────────────────────

def _scrub(msg: str) -> str:
    """A message about a PRIVATE file, cut down to the parser's own opening
    clause. Everything from the first path separator, quote, bracket or `~`
    on is the FILE talking — a quoted LOCUS line, a record id, a path whose
    last part names the file (`.dna` originals are named after library
    entries; a path may hold spaces, so it cannot be matched as one token)."""
    msg = msg.split("\n")[0]
    msg = re.split(r"[\\/'\"(\[~]", msg, maxsplit=1)[0]
    msg = re.split(r":\s", msg, maxsplit=1)[0]
    return re.sub(r"\b[A-Za-z0-9_.-]*\d[A-Za-z0-9_.-]*\b", "#", msg).strip()[:100]


def _error_text(exc: BaseException, private: bool = False) -> str:
    msg = " ".join(str(exc).split()) if not private else _scrub(str(exc))
    return f"{type(exc).__name__}: {msg[:160]}"


def score_item(item: Item, readers, sc, tmp: Path, *, write: bool) -> dict:
    """One file through every reader, against one census of it."""
    data = item.data()
    text = decode_text(data) if item.fmt != "dna" else ""
    census = (dna_census(data) if item.fmt == "dna"
              else embl_census(text) if item.fmt == "embl" else genbank_census(text))
    row = {"segment": item.segment, "key": item.key, "fmt": item.fmt,
           "size": item.size, "producer": item.producer, "kind": item.kind,
           "label": item.label, "_path": str(item.path), "_private": item.private,
           "census_records": 1 if item.fmt == "dna" else len(census),
           "by_reader": {}}
    for r in readers:
        row["by_reader"][r] = score(item, r, sc, census, tmp, write=write)
    sigs = [rr.get("_sig") for rr in row["by_reader"].values()]
    if len(sigs) == 2 and all(s is not None for s in sigs):
        row["agreement"] = _first_diff(sigs[0], sigs[1]) or "same"
    for rr in row["by_reader"].values():
        rr.pop("_sig", None)
    return row


def score(item: Item, reader: str, sc, census, tmp: Path, *, write: bool) -> dict:
    row: dict = {"reader": reader}
    fn = READERS[reader]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        signal.signal(signal.SIGALRM, _alarm)
        signal.setitimer(signal.ITIMER_REAL, PARSE_TIMEOUT_S)
        t0 = time.perf_counter()
        try:
            records = fn(item, sc)
            row["outcome"] = "accepted" if records else "empty"
        except ParseTimeout as exc:
            records, row["outcome"], row["error"] = None, "timeout", _error_text(exc, item.private)
        except ValueError as exc:
            records, row["outcome"], row["error"] = None, "refused", _error_text(exc, item.private)
        except Exception as exc:  # noqa: BLE001 — classifying it IS the point
            records, row["outcome"], row["error"] = None, "crashed", _error_text(exc, item.private)
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        row["seconds"] = time.perf_counter() - t0
    leaked = [w for w in caught if issubclass(w.category, ResourceWarning)]
    caught = [w for w in caught if not issubclass(w.category, ResourceWarning)]
    if leaked:
        row["leaked_handles"] = len(leaked)
    row["warnings"] = len(caught)
    if caught:
        row["warning_kinds"] = sorted({
            _scrub(str(w.message)) if item.private else " ".join(str(w.message).split())[:60]
            for w in caught})[:5]
    if not records:
        return row
    row["n_records"] = len(records)
    row["issues"] = census_issues(records, census, item.fmt)
    if item.fmt == "dna":
        row["metadata_dropped"] = dna_metadata_dropped(records, census)
    row["_sig"] = signature(records)
    if write:
        row.update(write_side(records, reader, sc, tmp, private=item.private))
    return row


def write_side(records, reader: str, sc, tmp: Path, *, private: bool) -> dict:
    out: dict = {}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            texts = (write_biopython(records, sc) if reader == "biopython"
                     else write_splicecraft(records, sc, tmp))
    except Exception as exc:  # noqa: BLE001
        out["write"] = "refused" if isinstance(exc, ValueError) else "crashed"
        out["write_error"] = _error_text(exc, private)
        return out
    out["write"] = "ok"
    problems = [p for t in texts for rt in split_records(t) for p in _gb_problems(rt)]
    out["spec_problems"] = len(problems)
    if problems:
        out["spec_kinds"] = dict(Counter(map(_problem_kind, problems)).most_common(4))
    try:
        p = tmp / "reread.gb"
        p.write_text("".join(texts), encoding="utf-8")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if reader == "biopython":
                from Bio import SeqIO
                back = list(SeqIO.parse(str(p), "genbank"))
            else:
                back = open_like_the_app(sc, p, "genbank")
        diff = _first_diff(signature(records), signature(back))
        out["round_trip"] = diff or "same"
    except Exception as exc:  # noqa: BLE001
        out["round_trip"] = f"unreadable ({type(exc).__name__})"
    return out


def _problem_kind(p: str) -> str:
    """'line 12: 93 columns (max 80): ...' -> 'over 80 columns'."""
    body = p.split(":", 1)[1].strip() if p.startswith("line ") else p
    if "columns (max 80)" in body:
        return "over 80 columns"
    if body.startswith("non-ASCII"):
        return "non-ASCII"
    return re.sub(r"[0-9]+|'[^']*'|\"[^\"]*\"", "#", body)[:50]


# ── the run ───────────────────────────────────────────────────────────────────

def sandbox() -> Path:
    """Point SpliceCraft's data dir and log into a temp dir BEFORE it is
    imported (its data dir is fixed at import time)."""
    root = Path(tempfile.mkdtemp(prefix="sc-verify-scoreboard-"))
    os.environ["XDG_DATA_HOME"] = str(root / "data")
    os.environ.pop("SPLICECRAFT_DATA_DIR", None)
    os.environ["SPLICECRAFT_LOG"] = str(root / "splicecraft.log")
    os.environ.setdefault("SPLICECRAFT_SKIP_LOCK", "1")
    return root


def import_splicecraft(root: Path):
    import splicecraft as sc
    if not str(sc._state._DATA_DIR).startswith(str(root)):
        raise SystemExit(f"refusing to run: SpliceCraft's data dir is "
                         f"{sc._state._DATA_DIR}, not the sandbox")
    return sc


def run(args) -> int:
    data_dir = real_data_dir()                # read-only, resolved pre-sandbox
    root = sandbox()
    sc = import_splicecraft(root)
    import Bio
    segments = [s.strip() for s in args.segments.split(",") if s.strip()]
    unknown = sorted(set(segments) - set(SEGMENTS))
    if unknown:
        raise SystemExit(f"unknown segment(s) {unknown}; choose from {', '.join(SEGMENTS)}")
    corpus = corpus_dir()
    work = root / "work"
    tmp = root / "tmp"
    tmp.mkdir(parents=True)
    items = discover(segments, corpus=corpus, data_dir=data_dir,
                     local=[Path(p).expanduser() for p in args.local], workdir=work,
                     sample=args.sample, seed=args.seed,
                     fuzz_bases=args.fuzz_bases, fuzz_per_base=args.fuzz_per_base)
    if not items:
        print("No corpus items found. Run `fetch` first, or pass --local.")
        return 1
    readers = [r.strip() for r in args.readers.split(",") if r.strip()]
    print(f"{len(items)} files · {sum(i.size for i in items) / 1e6:.1f} MB · "
          f"readers: {', '.join(readers)} · Biopython {Bio.__version__} · "
          f"SpliceCraft {sc.__version__}")
    rows = []
    t_start = time.perf_counter()
    harness_errors: "list[str]" = []
    for n, item in enumerate(items, 1):
        try:
            rows.append(score_item(item, readers, sc, tmp,
                                   write=not args.no_write and item.segment != "fuzz"))
        except Exception as exc:  # noqa: BLE001 — one bad item must not cost the run
            harness_errors.append(f"{item.segment}/{item.key}: {_error_text(exc, item.private)}")
        if n % 100 == 0 or n == len(items):
            el = time.perf_counter() - t_start
            print(f"  {n}/{len(items)} files · {el:.0f} s", flush=True)
    report = build_report(rows, readers)
    print_report(report, rows, readers, show=args.show, paths=args.paths)
    if harness_errors:
        print(f"\n{len(harness_errors)} item(s) could not be scored (harness, not reader):")
        for e in harness_errors[:10]:
            print(f"   {e}")
    meta = {"when": time.strftime("%Y-%m-%dT%H:%M:%S"), "biopython": Bio.__version__,
            "splicecraft": sc.__version__, "python": sys.version.split()[0],
            "segments": segments, "readers": readers, "sample": args.sample,
            "seed": args.seed}
    out_dir = corpus / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {"meta": meta, "summary": report, "items": public_rows(rows),
               "harness_errors": harness_errors}
    stamp = time.strftime("%Y%m%d-%H%M%S")
    for name in (f"scoreboard-{stamp}.json", "latest.json"):
        (out_dir / name).write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")
    print(f"\nSaved {out_dir / f'scoreboard-{stamp}.json'}")
    return 0


def public_rows(rows) -> list:
    """Rows as saved: no path, and no label for a private file. (Its error
    and warning texts were already scrubbed when they were recorded.)"""
    return [{k: v for k, v in r.items()
             if not k.startswith("_") and not (k == "label" and r["_private"])}
            for r in rows]


# ── the board ─────────────────────────────────────────────────────────────────

def _pct(a: int, b: int) -> str:
    return f"{100 * a / b:.1f}%" if b else "—"


def build_report(rows, readers) -> dict:
    report: dict = {}
    for seg in [s for s in SEGMENTS if any(r["segment"] == s for r in rows)]:
        srows = [r for r in rows if r["segment"] == seg]
        entry: dict = {"files": len(srows), "mb": round(sum(r["size"] for r in srows) / 1e6, 2),
                       "producers": dict(Counter(r["producer"] for r in srows).most_common()),
                       "readers": {}}
        for rd in readers:
            rr = [r["by_reader"][rd] for r in srows]
            acc = [x for x in rr if x["outcome"] == "accepted"]
            lossy = [x for x in acc if x.get("issues")]
            secs = [x["seconds"] for x in rr]
            mb_acc = sum(r["size"] for r, x in zip(srows, rr) if x["outcome"] == "accepted") / 1e6
            t_acc = sum(x["seconds"] for x in acc)
            e = {
                "accepted": len(acc),
                "refused": sum(1 for x in rr if x["outcome"] == "refused"),
                "crashed": sum(1 for x in rr if x["outcome"] == "crashed"),
                "timeout": sum(1 for x in rr if x["outcome"] == "timeout"),
                "empty": sum(1 for x in rr if x["outcome"] == "empty"),
                "silent_loss": len(lossy),
                "loss_kinds": dict(Counter(k for x in lossy for k in x["issues"]).most_common()),
                "with_warnings": sum(1 for x in rr if x.get("warnings")),
                "leaked_handles": sum(1 for x in rr if x.get("leaked_handles")),
                "crash_kinds": dict(Counter(x["error"].split(":")[0] for x in rr
                                            if x["outcome"] == "crashed").most_common()),
                "mb_per_s": round(mb_acc / t_acc, 2) if t_acc else None,
                "p50_ms": round(1000 * statistics.median(secs), 2) if secs else None,
                "p95_ms": round(1000 * sorted(secs)[int(0.95 * (len(secs) - 1))], 2) if secs else None,
                "write_refused": sum(1 for x in acc if x.get("write") in ("refused", "crashed")),
                "write_spec_files": sum(1 for x in acc if x.get("spec_problems")),
                "write_spec_kinds": dict(sum((Counter(x.get("spec_kinds") or {}) for x in acc),
                                             Counter()).most_common(6)),
                "round_trip_drift": sum(1 for x in acc if x.get("round_trip") not in (None, "same")),
                "round_trip_kinds": dict(Counter(x["round_trip"] for x in acc
                                                 if x.get("round_trip") not in (None, "same")).most_common()),
            }
            if any(r["fmt"] == "dna" for r in srows):
                e["dna_metadata_dropped"] = dict(Counter(
                    k for x in acc for k in x.get("metadata_dropped") or []).most_common())
            if seg == "synthetic":
                e["synthetic"] = _synthetic_grades(srows, rd)
            entry["readers"][rd] = e
        both = [r for r in srows if "agreement" in r]
        entry["both_accepted"] = len(both)
        entry["agree"] = sum(1 for r in both if r["agreement"] == "same")
        entry["disagree_kinds"] = dict(Counter(r["agreement"] for r in both
                                               if r["agreement"] != "same").most_common())
        report[seg] = entry
    return report


def _synthetic_grades(srows, rd) -> dict:
    """valid case: read faithfully, or it failed. malformed case: a clean
    refusal or a faithful read both pass; a mis-read or a crash fails."""
    grades = {}
    for r in srows:
        x = r["by_reader"][rd]
        ok_read = ((x["outcome"] == "accepted" and not x.get("issues"))
                   or (x["outcome"] == "empty" and r["census_records"] == 0))
        if r["kind"] == "valid":
            g = "pass" if ok_read else (
                "SILENT LOSS " + ",".join(x["issues"]) if x["outcome"] == "accepted"
                else x["outcome"].upper())
        else:
            g = "pass" if ok_read or x["outcome"] == "refused" else (
                "SILENT LOSS " + ",".join(x.get("issues") or []) if x["outcome"] == "accepted"
                else x["outcome"].upper())
        grades[r["label"]] = g
    return grades


def print_report(report, rows, readers, *, show: int, paths: bool) -> None:
    for seg, e in report.items():
        prod = ", ".join(f"{k} {v}" for k, v in e["producers"].items())
        print(f"\n━━ {seg}: {e['files']} files, {e['mb']} MB  [{prod}]")
        hdr = (f"   {'reader':<12}{'accepted':>9}{'refused':>8}{'crashed':>8}{'silent':>8}"
               f"{'warned':>8}{'MB/s':>8}{'p50 ms':>8}{'p95 ms':>8}{'w-refuse':>9}{'w-spec':>7}{'w-drift':>8}")
        print(hdr)
        for rd in readers:
            x = e["readers"][rd]
            print(f"   {rd:<12}{x['accepted']:>9}{x['refused']:>8}{x['crashed'] + x['timeout']:>8}"
                  f"{x['silent_loss']:>8}{x['with_warnings']:>8}{x['mb_per_s'] or 0:>8}"
                  f"{x['p50_ms'] or 0:>8}{x['p95_ms'] or 0:>8}{x['write_refused']:>9}"
                  f"{x['write_spec_files']:>7}{x['round_trip_drift']:>8}")
        for rd in readers:
            x = e["readers"][rd]
            notes = []
            if x["loss_kinds"]:
                notes.append("silent loss: " + ", ".join(f"{k} {v}" for k, v in x["loss_kinds"].items()))
            if x["crash_kinds"]:
                notes.append("crashes: " + ", ".join(f"{k} {v}" for k, v in x["crash_kinds"].items()))
            if x.get("dna_metadata_dropped"):
                notes.append("not brought back: " + ", ".join(
                    f"{k} {v}" for k, v in x["dna_metadata_dropped"].items()))
            if x["write_spec_kinds"]:
                notes.append("spec errors on write: " + ", ".join(
                    f"{k} ×{v}" for k, v in x["write_spec_kinds"].items()))
            if x["round_trip_kinds"]:
                notes.append("round-trip drift: " + ", ".join(
                    f"{k} {v}" for k, v in x["round_trip_kinds"].items()))
            for nline in notes:
                print(f"   {rd}: {nline}")
            if x.get("synthetic"):
                bad = {k: v for k, v in x["synthetic"].items() if v != "pass"}
                print(f"   {rd}: synthetic {len(x['synthetic']) - len(bad)}/{len(x['synthetic'])} pass"
                      + ("" if not bad else " — " + "; ".join(f"{k}: {v}" for k, v in bad.items())))
        if e.get("both_accepted"):
            print(f"   agreement: {e['agree']}/{e['both_accepted']} identical"
                  + ("" if not e["disagree_kinds"] else " — differ in " + ", ".join(
                      f"{k} {v}" for k, v in e["disagree_kinds"].items())))
    if show:
        print("\n━━ failures (first per reader)")
        shown = 0
        for r in rows:
            for rd, x in r["by_reader"].items():
                if x["outcome"] in ("crashed", "timeout") or x.get("issues"):
                    who = r["_path"] if (paths or not r["_private"]) else f"<private {r['key']}>"
                    why = x.get("error") or ("silent loss: " + ",".join(x["issues"]))
                    print(f"   [{r['segment']}/{rd}] {who}: {why}")
                    shown += 1
                    if shown >= show:
                        return


# ── fetch: the public corpus ─────────────────────────────────────────────────

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
# Hand-picked for structure, not biology: the repo's reference plasmid,
# small and large genomes with introns, trans-splicing, origin-adjacent
# features, very long /translation values, and a protein record.
_NCBI_CURATED = (
    ("nuccore", "L09137"), ("nuccore", "J01749"),
    ("nuccore", "J02459"), ("nuccore", "NC_045512"), ("nuccore", "NC_012920"),
    ("nuccore", "NC_001133"), ("nuccore", "U49845"), ("nuccore", "NC_005816"),
    ("nuccore", "NC_000932"), ("protein", "NP_000508"), ("protein", "P69905"),
)
_NCBI_SPEED = ("nuccore", "NC_000913.3")           # 4.6 Mb, ~4,500 genes
_NCBI_VECTOR_QUERIES = (
    '"cloning vector"[Title] AND 2000:20000[SLEN]',
    '"expression vector"[Title] AND 2000:20000[SLEN]',
    '"shuttle vector"[Title] AND 2000:30000[SLEN]',
)


def _get(url: str, *, timeout: float = 120.0) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _eutils(endpoint: str, **params) -> bytes:
    params.setdefault("tool", "splicecraft-format-scoreboard")
    time.sleep(0.4)                                  # NCBI: <= 3 requests/s without a key
    return _get(_EUTILS + endpoint + "?" + urllib.parse.urlencode(params))


def _save_records(text: str, out: Path) -> int:
    n = 0
    for rec in split_records(text):
        m = re.match(r"LOCUS\s+(\S+)", rec)
        if not m:
            continue
        name = re.sub(r"[^A-Za-z0-9_.-]", "_", m.group(1))
        (out / f"{name}.{'gp' if ' aa ' in rec.split(chr(10), 1)[0] else 'gb'}").write_text(rec, encoding="utf-8")
        n += 1
    return n


def fetch_ncbi(out: Path, *, n_sample: int, speed: bool) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for db in ("nuccore", "protein"):
        ids = [i for d, i in _NCBI_CURATED if d == db]
        text = _eutils("efetch.fcgi", db=db, id=",".join(ids),
                       rettype="gbwithparts" if db == "nuccore" else "gb",
                       retmode="text").decode("utf-8", "replace")
        print(f"  curated {db}: {_save_records(text, out)} records")
    ids: "list[str]" = []
    per = max(1, n_sample // len(_NCBI_VECTOR_QUERIES))
    for term in _NCBI_VECTOR_QUERIES:
        body = json.loads(_eutils("esearch.fcgi", db="nuccore", term=term, retmax=per, retmode="json"))
        ids += [i for i in body.get("esearchresult", {}).get("idlist", []) if i not in ids]
    for k in range(0, len(ids), 50):
        text = _eutils("efetch.fcgi", db="nuccore", id=",".join(ids[k:k + 50]),
                       rettype="gb", retmode="text").decode("utf-8", "replace")
        print(f"  vectors {k + 1}-{min(k + 50, len(ids))}: {_save_records(text, out)} records")
    if speed:
        db, acc = _NCBI_SPEED
        text = _eutils("efetch.fcgi", db=db, id=acc, rettype="gbwithparts", retmode="text").decode("utf-8", "replace")
        print(f"  speed record: {_save_records(text, out)} record, {len(text) / 1e6:.1f} MB")


def fetch_biopython(out: Path, ref: str) -> None:
    api = "https://api.github.com/repos/biopython/biopython/contents/Tests/{}?ref={}"
    for sub in ("GenBank", "EMBL", _CS_NAME):
        listing = json.loads(_get(api.format(sub, ref)))
        dest = out / sub.lower()
        dest.mkdir(parents=True, exist_ok=True)
        n = 0
        for ent in listing:
            name = ent.get("name", "")
            if ent.get("type") != "file" or fmt_for(Path(name)) is None:
                continue
            (dest / name).write_bytes(_get(ent["download_url"]))
            n += 1
        print(f"  biopython {sub.lower()}: {n} files")


def cmd_fetch(args) -> int:
    import Bio
    corpus = corpus_dir()
    print(f"corpus: {corpus}")
    if not args.skip_ncbi:
        fetch_ncbi(corpus / "ncbi", n_sample=args.ncbi, speed=not args.no_speed)
    if not args.skip_biopython:
        ref = args.biopython_ref or "biopython-" + Bio.__version__.replace(".", "")
        fetch_biopython(corpus / "biopython", ref)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="download the public corpus")
    f.add_argument("--ncbi", type=int, default=300, help="vector records to sample (default 300)")
    f.add_argument("--no-speed", action="store_true", help="skip the 4.6 Mb speed record")
    f.add_argument("--biopython-ref", default="", help="git ref for the fixtures (default: installed version)")
    f.add_argument("--skip-ncbi", action="store_true")
    f.add_argument("--skip-biopython", action="store_true")
    r = sub.add_parser("run", help="score the readers")
    r.add_argument("--segments", default=",".join(SEGMENTS))
    r.add_argument("--readers", default="biopython,splicecraft")
    r.add_argument("--local", action="append", default=[], help="a private folder or file (repeatable)")
    r.add_argument("--sample", type=int, default=0, help="at most N files per segment")
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--fuzz-bases", type=int, default=40)
    r.add_argument("--fuzz-per-base", type=int, default=12)
    r.add_argument("--no-write", action="store_true", help="skip the write side")
    r.add_argument("--show", type=int, default=0, help="list N failures")
    r.add_argument("--paths", action="store_true", help="name private files in --show output")
    args = ap.parse_args(argv)
    return cmd_fetch(args) if args.cmd == "fetch" else run(args)


if __name__ == "__main__":
    sys.exit(main())
