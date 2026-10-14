"""Independent oracles for sequence-file formats: spec validators and a
structural census.

Shared by `test_format_conformance.py` (are SpliceCraft's EXPORTS to spec?)
and `scripts/format_scoreboard.py` (what do READERS do with real files?). The
whole value of this module is that it is NOT the code under test: it never
imports SpliceCraft or Biopython, and each rule is re-derived from the format
documentation — the NCBI GenBank release notes + INSDC feature table, the ENA
flat-file layout, GFF3 1.26 — and from the `.dna` packet layout.
"""

from __future__ import annotations

import html
import re


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


# ── structural census: what a file DECLARES, counted without a parser ────────
#
# A reader that drops or invents content cannot be caught by comparing it with
# itself, and comparing it with Biopython only finds where the two DIFFER. The
# census counts what the file text plainly holds — records, bases, features,
# qualifiers — with the forgiving rules a reader should follow, so "accepted
# the file but lost its features" shows up as a number that disagrees. Same
# discipline as the validators above: no call into the code under test, and
# no Biopython.

_CENSUS_UNITS = ("bp", "aa", "rc")
_FEATURE_KEY_RE = re.compile(r"[A-Za-z0-9_'*-]+")
_QUALIFIER_START_RE = re.compile(r"/([A-Za-z0-9_][A-Za-z0-9_'*-]*)(=|$)")
_EMBL_ID_LEN_RE = re.compile(r"(\d+)\s+(BP|AA)\.?\s*$", re.IGNORECASE)
# What may follow a feature table in column 1 (Biopython's own list).
_GB_TABLE_ENDERS = frozenset({"BASE", "WGS", "WGS_SCAFLD", "TSA", "TLS"})


def _census_lines(text: str) -> "list[str]":
    if text.startswith("﻿"):
        text = text[1:]
    return [ln.expandtabs(8).rstrip()
            for ln in text.replace("\r\n", "\n").replace("\r", "\n").split("\n")]


def _census_record(declared, topology) -> dict:
    return {"declared_len": declared, "topology": topology, "seq_len": None,
            "contig": False, "n_features": 0, "n_qualifiers": 0,
            "keys": [], "qualifier_names": []}


def _locus_census(line: str) -> "tuple[int | None, str | None]":
    """Length and topology from a LOCUS line, read by TOKEN: the length is the
    last ``<digits> <unit>`` pair, the topology a token after the unit."""
    toks = line.split()
    for i in range(len(toks) - 1, 1, -1):
        if (toks[i].lower() in _CENSUS_UNITS
                and toks[i - 1].isascii() and toks[i - 1].isdigit()):
            tail = [t.lower() for t in toks[i + 1:]]
            topo = ("circular" if "circular" in tail
                    else "linear" if "linear" in tail else None)
            return int(toks[i - 1]), topo
    return None, None


class _FeatureTally:
    """Counts feature keys and qualifier starts in a feature table, one body
    line at a time. A qualifier's quoted value may run over several lines,
    and one of those lines may itself start with `/` — so quote state is
    tracked and continuation lines are never counted as qualifiers."""

    def __init__(self, rec: dict, key_indent_max: int = 19):
        self.rec = rec
        self.key_indent_max = key_indent_max
        self.open_quote = False

    def _is_key_line(self, body: str) -> bool:
        s = body.strip()
        if not s or s.startswith("/"):
            return False
        first = s.split()[0]
        indent = len(body) - len(body.lstrip(" "))
        return (indent <= self.key_indent_max and not first.isdigit()
                and bool(_FEATURE_KEY_RE.fullmatch(first)))

    def feed(self, body: str) -> None:
        if self.open_quote:
            if not self._is_key_line(body):
                if body.count('"') % 2 == 1:
                    self.open_quote = False
                return
            # A value never continues at the KEY indent: the quote was never
            # closed, and this line starts the next feature.
            self.open_quote = False
        s = body.strip()
        if not s:
            return
        if s.startswith("/"):
            m = _QUALIFIER_START_RE.match(s)
            if m and self.rec["n_features"] == 0:
                return                          # belongs to no feature
            if m:
                self.rec["n_qualifiers"] += 1
                self.rec["qualifier_names"].append(m.group(1))
                value = s.split("=", 1)[1] if "=" in s else ""
                if value.startswith('"') and value.count('"') % 2 == 1:
                    self.open_quote = True
            return
        if self._is_key_line(body):
            self.rec["n_features"] += 1
            self.rec["keys"].append(s.split()[0])
        # anything else is a location continuing from the line above


def _sequence_letters(line: str) -> int:
    """Bases on one sequence line: everything that is not whitespace or a
    digit — except a leading token that contains a digit, which is the
    position number (GenBank puts it first) even when a stray character has
    landed in it. EMBL's trailing count is digits, so it drops out anyway."""
    toks = line.split()
    if len(toks) > 1 and any(c.isdigit() for c in toks[0]):
        toks = toks[1:]
    return sum(1 for tok in toks for c in tok if not c.isdigit())


def _close_sequence(rec: "dict | None") -> None:
    """A sequence block holding no letters says nothing about the length
    (EMBL patent records declare "Sequence 358 AA" and ship no bases), so
    report it as unknown rather than as zero."""
    if rec is not None and rec["seq_len"] == 0:
        rec["seq_len"] = None


def genbank_census(text: str) -> "list[dict]":
    """Per record of GenBank `text`, what the file declares:

    ``declared_len`` / ``topology`` (LOCUS, by token), ``seq_len`` (characters
    in the ORIGIN block that are neither whitespace nor the position numbers;
    None for a record without one, e.g. CONTIG-only), ``contig``,
    ``n_features`` + ``keys``, ``n_qualifiers`` + ``qualifier_names``.

    Read forgivingly — tabs are expanded, the FEATURES and ORIGIN headers are
    recognised by keyword whatever their spacing, a key is any token indented
    less than the qualifier column, a missing `//` ends the record at EOF —
    because the point is to know what a tolerant reader SHOULD have seen."""
    recs: "list[dict]" = []
    cur = None
    state = ""
    tally = None
    for line in _census_lines(text):
        if line.startswith("LOCUS") and (len(line) == 5 or line[5] == " "):
            cur = _census_record(*_locus_census(line))
            recs.append(cur)
            state, tally = "header", _FeatureTally(cur)
            continue
        if cur is None or tally is None:
            continue
        if line.strip() == "//":
            _close_sequence(cur)
            cur, state = None, ""
            continue
        if state == "origin":
            cur["seq_len"] += _sequence_letters(line)
            continue
        if line[:1] not in ("", " "):           # something in column 1
            word = line.split()[0].upper()
            if word == "FEATURES":
                state = "features"
            elif word == "ORIGIN":
                state = "origin"
                cur["seq_len"] = 0
            elif word == "CONTIG":
                state = "contig"
                cur["contig"] = True
            elif state == "features" and word not in _GB_TABLE_ENDERS:
                tally.feed(line)                # a key pushed into column 1
            else:
                state = "header"
            continue
        if state == "features":
            tally.feed(line)
    _close_sequence(cur)
    return recs


def embl_census(text: str) -> "list[dict]":
    """`genbank_census` for EMBL: records open on `ID`, feature lines carry an
    `FT` prefix, the sequence follows `SQ` and carries its counts at the END
    of each line."""
    recs: "list[dict]" = []
    cur = None
    state = ""
    tally = None
    for line in _census_lines(text):
        tag = line[:2].upper()
        if tag == "ID" and line[2:3] in ("", " "):
            m = _EMBL_ID_LEN_RE.search(line)
            low = line.lower()
            topo = ("circular" if "circular" in low
                    else "linear" if "linear" in low else None)
            cur = _census_record(int(m.group(1)) if m else None, topo)
            recs.append(cur)
            state, tally = "header", _FeatureTally(cur, key_indent_max=3)
            continue
        if cur is None or tally is None:
            continue
        if line.strip() == "//":
            _close_sequence(cur)
            cur, state = None, ""
            continue
        if state == "sequence":
            cur["seq_len"] += _sequence_letters(line)
            continue
        if tag == "SQ":
            state = "sequence"
            cur["seq_len"] = 0
        elif tag == "CO":
            cur["contig"] = True
        elif tag == "FT":
            tally.feed(line[2:])
    _close_sequence(cur)
    return recs


def dna_census(data: bytes) -> dict:
    """What a `.dna` file holds, from its packet stream (1-byte type, 4-byte
    big-endian length, payload) rather than from any reader of the format.

    ``seq_len`` / ``topology`` (DNA packet 0x00), ``n_features`` +
    ``feature_names`` + ``feature_types`` (``<Feature`` elements of packet
    0x0A; a missing type is ``misc_feature``, as the format defines), ``n_colored``
    (features carrying a colour), ``n_primers`` / ``n_binding_sites``
    (packet 0x05), ``has_history`` (packet 0x07), ``packet_types``, and
    ``truncated`` when a declared length runs past the end of the file."""
    out: dict = {"seq_len": None, "topology": None, "n_features": 0,
                 "feature_names": [], "feature_types": [], "n_colored": 0,
                 "n_primers": 0,
                 "n_binding_sites": 0, "has_history": False,
                 "packet_types": [], "truncated": False}
    pos, n = 0, len(data)
    while pos + 5 <= n:
        kind = data[pos]
        size = int.from_bytes(data[pos + 1:pos + 5], "big")
        payload = data[pos + 5:pos + 5 + size]
        if len(payload) < size:
            out["truncated"] = True
            break
        out["packet_types"].append(kind)
        if kind == 0x00 and size >= 1:
            out["seq_len"] = size - 1
            out["topology"] = "circular" if payload[0] & 0x01 else "linear"
        elif kind == 0x0A:
            xml = payload.decode("utf-8", "replace")
            heads = re.findall(r"<Feature\b[^>]*>", xml)
            out["n_features"] = len(heads)
            out["feature_names"] = [
                html.unescape(m.group(1)) if (m := re.search(r'\bname="([^"]*)"', h))
                else "" for h in heads]
            out["feature_types"] = [
                html.unescape(m.group(1)) if (m := re.search(r'\btype="([^"]*)"', h))
                else "misc_feature" for h in heads]
            out["n_colored"] = sum(
                1 for block in re.findall(r"<Feature\b.*?</Feature>", xml, re.S)
                if 'color="#' in block)
        elif kind == 0x05:
            xml = payload.decode("utf-8", "replace")
            out["n_primers"] = len(re.findall(r"<Primer\b", xml))
            out["n_binding_sites"] = len(re.findall(r"<BindingSite\b", xml))
        elif kind == 0x07:
            out["has_history"] = True
        pos += 5 + size
    return out
