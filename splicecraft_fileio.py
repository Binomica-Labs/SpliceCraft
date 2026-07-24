"""splicecraft_fileio — single-file sequence-format I/O (Phase D, layer L2).

The plain file-format read/write core, extracted from the hub: FASTA / AB1 /
FASTQ / GFF3 ingest (path -> SeqRecord) + GenBank / GFF / FASTA / EMBL export
(record -> path), with their codec helpers. Biopython is imported lazily inside
the entry points exactly as in the hub. fileio now owns the WHOLE file + network
I/O layer: the .dna CommercialSaaS codec, the _save/_load/_delete_dna_original
blob store, the NCBI fetch_genbank / fetch_protein fetchers, and the zip /
Plasmidsaurus handling were all migrated here from the hub (data-safety
preserved — blob writes go through the persistence chokepoint, and NCBI fetch
reaches the egress guard via _state._demo_block_network_hook). Re-exported by the
hub so `sc.<name>` + every call site (modals, agent endpoints, the loader)
resolves unchanged.
"""
from __future__ import annotations

import os
import re
from copy import copy as _shallow_copy
from datetime import date as _date
from pathlib import Path

import splicecraft_state as _state
from splicecraft_logging import _log, _log_event, _timed
from splicecraft_persistence import (
    _atomic_write_text, _fsync_parent_dir, _refuse_unauthorized_delete,
    _refuse_unauthorized_write, _safe_file_size_check,
)
from splicecraft_biology import _rc
from splicecraft_util import (
    _CONTROL_CHARS_RE, _DEFAULT_TYPE_COLORS, _feat_bounds, _is_windows_reserved_stem,
    _natural_sort_key, _pick_single_record, _safe_xml_parse, _sanitize_label,
)
from splicecraft_record import (
    _GB_LOCUS_NAME_MAX, _gb_text_to_record, _normalize_primer_seq, _record_to_gb_text,
    _restore_display_name_from_comment,
)
from splicecraft_history import _CommercialSaaSHistoryNode, _history_human_dt
from splicecraft_net import _NCBI_TIMEOUT_S, _sanitize_accession


# Foreign-file ingest cap (GenBank / .dna / GFF3 the user *opens*, NOT the
# user's own data-dir JSON — that keeps the 1 GB cap above). A file passes
# the on-disk byte check and is then eagerly materialised by Biopython into
# a SeqRecord several× its size in RAM, so a ~1 GB record could OOM the
# process. No real plasmid / construct / chromosome an editor opens is this
# big; 256 MB is a generous ceiling that blocks the OOM without refusing any
# plausible file. (Chosen 2026-06-14; see attack-surface sweep.)
_GB_INGEST_MAX_BYTES = 256 * 1024 * 1024   # 256 MB

# Per-file size cap. Plasmids are typically <100 KB; anything larger is
# probably a chromosome dump or an unrelated file. Refuse rather than
# OOM the parser. 50 MB leaves headroom for huge cosmids / BACs.
_BULK_IMPORT_MAX_BYTES = 50 * 1024 * 1024

# AB1 Sanger trace files (single record, base-called via BioPython's
# `abi` SeqIO key). Plasmidsaurus zips include these and we now surface
# them as loadable entries.
_AB1_EXTS: frozenset[str] = frozenset({".ab1", ".abi"})

# FASTQ multi-read files (sequencing output). Imported as a new
# collection per file, same flow as multi-record FASTA. Quality scores
# are read but not surfaced in the plasmid map (the editor's primary
# axis is annotation, not basecall quality).
_FASTQ_EXTS: frozenset[str] = frozenset({".fastq", ".fq"})

# FASTQ multi-read import caps. Without them a 1 GB Illumina FASTQ
# (≈8 M reads) blows RAM before any guard fires, and a 100k-read FASTQ
# would create 100k library entries blowing `_BULK_IMPORT_MAX_BYTES`
# on the serialised library JSON. Plasmid editor != read aligner.
_FASTQ_MAX_READS = 1000


def _normalize_for_genbank(record):
    """Return a shallow copy of `record` with NCBI-required fields filled in.

    Idempotent — existing values are preserved. Only fills gaps. Caller's
    record is never mutated.
    """
    from datetime import datetime as _dt

    rec = _shallow_copy(record)
    anns = dict(getattr(record, "annotations", None) or {})

    anns.setdefault("molecule_type", "DNA")
    # 2026-05-27 (audit-3 H1): default to "linear" not "circular" so a
    # record imported from GFF3 / FASTA without an explicit topology
    # is NOT silently re-labelled circular on first GenBank save.
    # Topology is biologically load-bearing — getting it wrong changes
    # how every downstream tool (PCR sim, primer design, restriction
    # scan wrap detection) treats the sequence. Pre-fix all topology-
    # less imports flipped to circular on the first save.
    anns.setdefault("topology", "linear")
    anns.setdefault("data_file_division", "SYN")

    if not anns.get("date"):
        anns["date"] = _dt.now().strftime("%d-%b-%Y").upper()

    if not anns.get("accessions"):
        acc = rec.id if rec.id and rec.id != "<unknown id>" else ""
        anns["accessions"] = [acc]

    if not anns.get("organism"):
        anns["organism"] = "synthetic construct"
    if not anns.get("source"):
        anns["source"] = anns["organism"]
    if not anns.get("taxonomy"):
        anns["taxonomy"] = ["other sequences", "artificial sequences"]

    rec.annotations = anns

    # LOCUS name: spec is 16 chars; NCBI accepts up to 28 in practice.
    # Biopython itself warns if >16 but does not fail.
    # 2026-05-27 (audit-3 M2): surface a log warning + notify when
    # truncation actually happens so the user knows the LOCUS no
    # longer matches their display name.
    if rec.name and len(rec.name) > _GB_LOCUS_NAME_MAX:
        original_name = rec.name
        rec.name = rec.name[:_GB_LOCUS_NAME_MAX]
        _log.warning(
            "GenBank export: LOCUS name truncated from %d to %d chars: "
            "%r → %r",
            len(original_name), _GB_LOCUS_NAME_MAX,
            original_name, rec.name,
        )
    if not rec.name or rec.name == "<unknown name>":
        rec.name = (rec.id or "PLASMID")[:_GB_LOCUS_NAME_MAX] or "PLASMID"

    if not rec.description or rec.description == "<unknown description>":
        rec.description = rec.name

    if not rec.id or rec.id == "<unknown id>":
        rec.id = rec.name

    return rec


def _export_genbank_to_path(record, path) -> dict:
    """Write `record` to `path` as a GenBank file. Atomic + round-trip verified.

    Returns a small summary dict `{"path", "bp", "features"}` for UI reporting.

    Raises:
      OSError on filesystem failures (write, replace, fsync).
      ValueError if the round-trip parse fails or the parsed record
        disagrees with the source on sequence length, sequence content,
        or feature count — meaning the export is not byte-safe.

    The round-trip happens BEFORE the target file is touched, so a failed
    export never leaves a half-written / corrupt .gb at `path`.
    """
    from pathlib import Path as _Path

    p = _Path(path).expanduser()
    normalized = _normalize_for_genbank(record)
    text = _record_to_gb_text(normalized)

    # Round-trip verify before touching the filesystem
    try:
        parsed = _gb_text_to_record(text)
    except Exception as exc:
        raise ValueError(f"export round-trip parse failed: {exc}") from exc
    if len(parsed.seq) != len(normalized.seq):
        raise ValueError(
            f"export round-trip sequence length mismatch "
            f"({len(parsed.seq)} vs {len(normalized.seq)})"
        )
    if str(parsed.seq).upper() != str(normalized.seq).upper():
        raise ValueError("export round-trip sequence content mismatch")
    if len(parsed.features) != len(normalized.features):
        raise ValueError(
            f"export round-trip feature count mismatch "
            f"({len(parsed.features)} vs {len(normalized.features)})"
        )
    # 2026-05-27 (audit-3 H3): feature-COUNT check is necessary but
    # not sufficient. Biopython has historic edge cases where a
    # CompoundLocation flattens on write — the count survives but
    # the wrap structure is lost. Compare the SIGNATURE of every
    # feature (type + location-string + sorted qualifier items)
    # before touching disk so a flatten-on-write surfaces as a
    # round-trip failure rather than silent corruption.
    def _feature_signature(feat) -> "tuple":
        loc_str = str(getattr(feat, "location", "") or "")
        ftype   = getattr(feat, "type", "") or ""
        quals   = getattr(feat, "qualifiers", None) or {}
        # Qualifier values are typically list[str]; sort by key
        # then stringify each value list so the signature is
        # order-stable across writer round-trips.
        qual_sig = tuple(sorted(
            (k, tuple(v) if isinstance(v, (list, tuple)) else (str(v),))
            for k, v in quals.items()
        ))
        return (ftype, loc_str, qual_sig)

    src_sigs = sorted(_feature_signature(f) for f in normalized.features)
    dst_sigs = sorted(_feature_signature(f) for f in parsed.features)
    if src_sigs != dst_sigs:
        # Compare as MULTISETS rather than zipping the two sorted lists:
        # one changed signature re-sorts, so the zip flags every feature
        # after it as divergent too and buries the one that actually moved.
        from collections import Counter as _Counter
        src_only = sorted((_Counter(src_sigs) - _Counter(dst_sigs)).elements())
        dst_only = sorted((_Counter(dst_sigs) - _Counter(src_sigs)).elements())

        def _handle(sig) -> str:
            """`CDS 'GeneX' at [10:40]` — something findable in the map."""
            ftype, loc_str, qual_sig = sig
            quals = dict(qual_sig)
            for key in ("label", "gene", "product", "note"):
                if quals.get(key):
                    return f"{ftype or 'feature'} {quals[key][0]!r} at {loc_str}"
            return f"{ftype or 'feature'} at {loc_str}"

        def _no_strand(loc_str: str) -> str:
            for suffix in ("(+)", "(-)", "(?)"):
                loc_str = loc_str.replace(suffix, "")
            return loc_str

        # Say what diverged, in words. This used to print the raw signature
        # tuples, which is unreadable in the common case: a feature built
        # with `strand=None` differs from its round-tripped self ONLY by a
        # `(?)` the reader has to spot inside two near-identical tuples.
        detail = ""
        if src_only:
            bad = src_only[0]
            # Same type + qualifiers ⇒ this is the same feature, come back
            # changed; that pairing is what makes a specific diagnosis safe.
            mate = next((d for d in dst_only
                          if d[0] == bad[0] and d[2] == bad[2]), None)
            if mate is not None and _no_strand(bad[1]) == _no_strand(mate[1]):
                if "(" not in bad[1]:
                    detail = (
                        f": {_handle(bad)} has no strand. GenBank cannot "
                        f"represent a strandless feature, so it parsed back "
                        f"as {mate[1]}. Set strand=1 (forward) or strand=-1 "
                        f"(reverse) on the feature, then export again."
                    )
                else:
                    detail = (f": {_handle(bad)} changed strand "
                              f"({bad[1]} → {mate[1]})")
            elif mate is not None:
                detail = (f": {_handle(bad)} moved "
                          f"({bad[1]} → {mate[1]})")
            else:
                detail = (f": {_handle(bad)} did not survive the round-trip "
                          f"— its type or qualifiers changed")
        raise ValueError(
            f"export round-trip feature signature mismatch — "
            f"{len(src_only)} of {len(src_sigs)} features diverged{detail}"
        )

    _atomic_write_text(p, text)

    _log.info(
        "Exported GenBank to %s (%d bp, %d features)",
        p, len(normalized.seq), len(normalized.features),
    )
    return {"path": str(p), "bp": len(normalized.seq),
            "features": len(normalized.features)}


def _record_to_gff3(record) -> str:
    """Serialise `record` to GFF3 text (specification 1.26).

    GFF3 columns: seqid, source, type, start, end, score, strand, phase,
    attributes. Coordinates are 1-based inclusive — note the off-by-one
    versus SpliceCraft's internal 0-based half-open `[start, end)`. For
    each feature we emit one line per `FeatureLocation` part so wrap
    features (origin-spanning `CompoundLocation`) become two
    same-ID rows joined by a shared `ID=...` attribute, the standard
    GFF3 convention for split features. Circular records carry
    `Is_circular=true` on a synthesised top-level region row.
    """
    from urllib.parse import quote as _q

    seqid = (record.id or record.name or "plasmid").strip() or "plasmid"
    # GFF3 reserves a tighter set of seqid characters than GenBank LOCUS;
    # percent-encode anything outside `[A-Za-z0-9._:^*$@!+_?\-|]`.
    safe_seqid = _q(seqid, safe=".:_-")
    n = len(record.seq)
    is_circular = (
        (record.annotations or {}).get("topology", "").lower() == "circular"
    )

    out: list[str] = []
    out.append("##gff-version 3")
    if n:
        out.append(f"##sequence-region {safe_seqid} 1 {n}")
    # Synthesise a top-level region row so downstream consumers can see
    # the topology flag — Bio.SeqRecord doesn't surface it through the
    # features list otherwise. `region` is the GFF3 SO term for this.
    region_attrs = [f"ID={safe_seqid}"]
    if is_circular:
        region_attrs.append("Is_circular=true")
    if n:
        out.append("\t".join((
            safe_seqid, "SpliceCraft", "region",
            "1", str(n), ".", "+", ".",
            ";".join(region_attrs),
        )))

    auto_id = 0
    for feat in record.features:
        ftype = (feat.type or "misc_feature").strip() or "misc_feature"
        if ftype == "source":
            # Source features map to the synthetic `region` row above;
            # emitting both would double-list the whole-record span.
            continue
        try:
            strand_int = int(feat.location.strand or 0)
        except (AttributeError, TypeError, ValueError):
            strand_int = 0
        # NOTE: per-part strand is computed below in the parts loop
        # (see `part_gff_strand`); the feature-level strand glyph was
        # an earlier-revision artefact and is unused after the
        # mixed-strand `CompoundLocation` fix.
        # Pull a display name + extra qualifiers. Same precedence other
        # parts of the codebase use.
        name = ""
        try:
            for q in ("label", "gene", "product"):
                v = feat.qualifiers.get(q)
                if isinstance(v, list) and v:
                    name = str(v[0])
                    break
        except AttributeError:
            pass
        auto_id += 1
        feat_id = f"feat{auto_id}"
        attr_parts: list[str] = [f"ID={feat_id}"]
        if name:
            attr_parts.append(f"Name={_q(name, safe='')}")
        # All other qualifiers as Note-prefixed key/value, except the
        # ones we already mapped (label/gene/product) — except gene and
        # product can convey extra info, so emit them too as proper
        # GFF3 attributes when present.
        try:
            quals = feat.qualifiers
        except AttributeError:
            quals = {}
        for k, vlist in (quals or {}).items():
            if k in ("label",):
                continue
            if not isinstance(vlist, list):
                vlist = [str(vlist)]
            joined = ",".join(_q(str(v), safe="") for v in vlist)
            if k in ("gene", "product"):
                attr_parts.append(f"{k}={joined}")
            else:
                # GFF3 spec is strict — treat any unknown qualifier as
                # a Note (free-text). Multiple Notes get comma-joined.
                attr_parts.append(f"Note={_q(k + '=' + ','.join(str(v) for v in vlist), safe='')}")
        # Phase: CDS features default to 0 unless the qualifier
        # supplies a codon_start (1-based 1/2/3 → 0/1/2 phase).
        phase = "."
        if ftype.upper() == "CDS":
            try:
                cs = int((quals.get("codon_start", ["1"]) or ["1"])[0])
            except (TypeError, ValueError, IndexError):
                cs = 1
            phase = str(max(0, min(2, cs - 1)))
        # Iterate location parts so wrap features get one row per arc,
        # sharing the same ID — the GFF3 split-feature convention.
        # Order parts in biological 5'→3' direction so the split-feature
        # rows read in the order the ribosome (or polymerase) traverses
        # them. For forward-strand: ascending genomic start. For
        # reverse-strand: descending. A canonical-wrap forward feature
        # (`join(tail..total, 1..head)`) is already declared in 5'→3'
        # order (tail first), so sorting by start would REVERSE it; we
        # detect canonical wrap and keep it in declared order. Same for
        # reverse-strand wrap (head first, descending). Pre-2026-05-11
        # the loop iterated `parts` in declared order — fine for
        # GenBank-parsed wraps (which declare in biological order) but
        # wrong for programmatically constructed CompoundLocations.
        try:
            parts_seq = list(getattr(feat.location, "parts", None)
                              or [feat.location])
        except AttributeError:
            parts_seq = []
        try:
            is_wrap_canonical = (
                len(parts_seq) == 2
                and int(parts_seq[0].start) == 0
                and int(parts_seq[-1].end) == len(record.seq)
                and int(parts_seq[0].end) < int(parts_seq[-1].start)
            )
        except (AttributeError, TypeError, ValueError):
            is_wrap_canonical = False
        if is_wrap_canonical:
            # Tail first for + strand (biological 5'→3'); head first for - strand.
            parts = ([parts_seq[-1], parts_seq[0]] if strand_int != -1
                     else [parts_seq[0], parts_seq[-1]])
        else:
            try:
                parts_sorted = sorted(parts_seq, key=lambda p: int(p.start))
            except (AttributeError, TypeError, ValueError):
                parts_sorted = parts_seq
            parts = (parts_sorted if strand_int != -1
                     else list(reversed(parts_sorted)))
        for part in parts:
            try:
                p_s = int(part.start)
                p_e = int(part.end)
            except (AttributeError, TypeError, ValueError):
                continue
            if p_e <= p_s:
                continue
            # Per-part strand: a mixed-strand `CompoundLocation` has
            # `feat.location.strand == None` (Biopython returns None when
            # parts disagree), so the feature-level `gff_strand`
            # computation above flattened to "." and lost the per-part
            # info. Probe each `part.strand` so mixed-strand joins
            # (rare but legal in biology) emit the right +/- per arc.
            # Falls back to the feature-level strand when the part
            # itself has no strand attribute.
            part_strand = getattr(part, "strand", None)
            if part_strand is None:
                part_strand = strand_int or 0
            part_gff_strand = ("+" if part_strand == 1
                                else "-" if part_strand == -1 else ".")
            out.append("\t".join((
                safe_seqid,
                "SpliceCraft",
                ftype,
                # 1-based inclusive: start = p_s + 1; end = p_e.
                str(p_s + 1),
                str(p_e),
                ".",
                part_gff_strand,
                phase,
                ";".join(attr_parts),
            )))

    out.append("")  # trailing newline so cat-friendly tools don't
                    # complain about a missing final EOL.
    return "\n".join(out)


def _parse_gff3_text(text: str) -> dict:
    """Parse GFF3 text into a structured dict for the loader.

    Returns ``{seqid, length, is_circular, features, fasta_seq}``:
      * ``seqid``        — first seqid encountered (from sequence-region
                            directive, region row, or the first feature)
      * ``length``       — value from ``##sequence-region`` or the region
                            row, else None.
      * ``is_circular``  — True if the synthesised region row carries
                            ``Is_circular=true``.
      * ``features``     — list of dicts: ``{type, start_0, end, strand,
                            qualifiers, gff_id, phase}``. Wrap split-feat
                            rows share the same ``gff_id``; callers merge
                            them into one CompoundLocation.
      * ``fasta_seq``    — sequence string from the inline ``##FASTA``
                            directive (if present), else None.

    Coordinates are converted from GFF3 1-based inclusive to SpliceCraft
    0-based half-open. Raises ValueError on malformed input.
    """
    from urllib.parse import unquote as _u
    seqid = None
    length = None
    is_circular = False
    features: list[dict] = []
    fasta_seq: "str | None" = None
    in_fasta = False
    fasta_buf: list[str] = []

    for raw in text.splitlines():
        line = raw.rstrip("\n").rstrip("\r")
        if in_fasta:
            if line.startswith(">"):
                # 2026-05-27 (audit-3 M8): multi-record `##FASTA`
                # blocks are silently parsed as "first wins" pre-fix
                # — but the second `>` header losing its record is a
                # data-loss bug for the rare GFF3 file that bundles
                # multiple sequences. Raise loud so the user knows
                # the file isn't single-record and can split it.
                if fasta_buf:
                    raise ValueError(
                        "GFF3 ##FASTA block contains multiple records; "
                        "SpliceCraft requires a single-record ##FASTA. "
                        "Split the file and import each record "
                        "separately."
                    )
                continue
            if line:
                fasta_buf.append(line.strip())
            continue
        if not line:
            continue
        if line.startswith("##FASTA"):
            in_fasta = True
            continue
        if line.startswith("##sequence-region"):
            parts = line.split()
            if len(parts) >= 4:
                if seqid is None:
                    seqid = _u(parts[1])
                try:
                    length = int(parts[3])
                except ValueError:
                    pass
            continue
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        seq_col = _u(cols[0])
        if seqid is None:
            seqid = seq_col
        ftype = cols[2].strip()
        try:
            start_1 = int(cols[3])
            end_1   = int(cols[4])
        except ValueError:
            continue
        if end_1 < start_1:
            continue
        strand = cols[6].strip()
        # Sweep #27: validate strand against the GFF3 spec
        # ({"+", "-", ".", "?"}). Pre-sweep ANY string silently mapped
        # to strand 0; a malformed GFF3 with embedded shell metas /
        # terminal escapes / HTML in the strand column would parse and
        # potentially surface in user-facing toasts unescaped. We now
        # SKIP the row (with a debug log) when the strand column isn't
        # one of the canonical four; matches the silent-skip semantics
        # GFF3 readers conventionally use for malformed records.
        if strand not in ("+", "-", ".", "?"):
            _log.debug(
                "GFF3 row skipped: invalid strand %r (must be one of "
                "'+', '-', '.', '?')", strand[:32],
            )
            continue
        strand_int = 1 if strand == "+" else (-1 if strand == "-" else 0)
        phase_s = cols[7].strip()
        try:
            phase = int(phase_s) if phase_s != "." else None
        except ValueError:
            phase = None
        attrs_raw = cols[8]
        quals: dict[str, list[str]] = {}
        gff_id = ""
        for kv in attrs_raw.split(";"):
            kv = kv.strip()
            if not kv or "=" not in kv:
                continue
            k, _, v = kv.partition("=")
            k = _u(k)
            vals = [_u(x) for x in v.split(",")]
            if k == "ID":
                gff_id = vals[0] if vals else ""
                continue
            if k == "Name":
                quals.setdefault("label", []).extend(vals)
                continue
            if k == "Is_circular":
                if ftype == "region" and vals and vals[0].lower() == "true":
                    is_circular = True
                continue
            if k == "Note":
                # `Note=key=value` round-trip (mirrors `_record_to_gff3`'s
                # encoding for unknown qualifiers).
                for v_one in vals:
                    if "=" in v_one:
                        nk, _, nv = v_one.partition("=")
                        quals.setdefault(nk, []).extend(
                            nv.split(",") if nv else [""]
                        )
                    else:
                        quals.setdefault("note", []).append(v_one)
                continue
            quals.setdefault(k, []).extend(vals)
        if ftype == "region":
            if length is None:
                length = end_1
            continue
        if ftype == "source":
            continue
        features.append({
            "type":       ftype,
            "start_0":    start_1 - 1,
            "end":        end_1,
            "strand":     strand_int,
            "qualifiers": quals,
            "gff_id":     gff_id,
            "phase":      phase,
        })
    if fasta_buf:
        fasta_seq = "".join(fasta_buf).upper()

    return {
        "seqid":       seqid or "plasmid",
        "length":      length,
        "is_circular": is_circular,
        "features":    features,
        "fasta_seq":   fasta_seq,
    }


def _gff3_features_to_biopython(
    parsed: dict, total: int,
) -> list:
    """Convert parsed GFF3 feature rows into BioPython SeqFeature
    objects. Same-`gff_id` rows are rejoined as a CompoundLocation
    (canonical wrap inverse of `_record_to_gff3`). Rows with coords
    outside ``[0, total)`` are dropped with a per-row warning log
    (sweep #25 2026-05-23 — pre-fix the drop was silent, so a GFF3
    with one bad row in a multi-part feature destroyed the entire
    feature with no user-visible signal).
    """
    from Bio.SeqFeature import (
        SeqFeature, FeatureLocation, CompoundLocation,
    )
    by_id: dict[str, list[dict]] = {}
    no_id: list[dict] = []
    for f in parsed["features"]:
        if f["gff_id"]:
            by_id.setdefault(f["gff_id"], []).append(f)
        else:
            no_id.append(f)

    out: list = []

    def _make_loc(parts: list[dict], gid: str = ""):
        # Build FeatureLocation per part, then merge into a
        # CompoundLocation when there are 2+ parts (wrap features).
        locs = []
        for p in parts:
            s, e = p["start_0"], p["end"]
            if s < 0 or e > total or e <= s:
                # Sweep #25: log instead of silent drop. Use INFO not
                # WARNING so a noisy GFF3 (e.g. coordinates relative
                # to a different reference) doesn't spam the log.
                _log.info(
                    "GFF3: dropping out-of-range feature %r "
                    "(type=%s start=%d end=%d total=%d)",
                    gid or "(no-id)", p.get("type") or "?", s, e, total,
                )
                return None
            locs.append(FeatureLocation(s, e, strand=p["strand"] or 0))
        if not locs:
            return None
        if len(locs) == 1:
            return locs[0]
        return CompoundLocation(locs)

    for gid, parts in by_id.items():
        loc = _make_loc(parts, gid=gid)
        if loc is None:
            continue
        first = parts[0]
        feat = SeqFeature(loc, type=first["type"],
                            qualifiers=dict(first["qualifiers"]))
        out.append(feat)
    for f in no_id:
        loc = _make_loc([f], gid=f.get("type") or "")
        if loc is None:
            continue
        feat = SeqFeature(loc, type=f["type"],
                            qualifiers=dict(f["qualifiers"]))
        out.append(feat)
    return out


def _gff3_path_to_record(path: str):
    """Load a GFF3 file as a SeqRecord.

    Requires an inline ``##FASTA`` directive — GFF3 alone carries no
    sequence and SpliceCraft is a sequence editor. Topology is set from
    the synthesised region row's ``Is_circular=true`` attribute (mirrors
    `_record_to_gff3`'s export convention).

    For sequence-less GFF3 files, use ``_gff3_apply_to_loaded_record``
    instead — that path treats the file as a feature-transfer overlay
    on the currently-loaded plasmid.

    Raises ValueError on malformed GFF3 or missing ``##FASTA`` section.
    """
    from pathlib import Path as _P
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    ok, reason = _safe_file_size_check(
        _P(path), _GB_INGEST_MAX_BYTES, "GFF3",
    )
    if not ok:
        raise ValueError(reason or "GFF3 file rejected")
    text = _P(path).read_text(encoding="utf-8", errors="replace")
    parsed = _parse_gff3_text(text)
    if parsed["fasta_seq"] is None:
        raise ValueError(
            "GFF3 file has no ##FASTA section. Standalone GFF3 import "
            "requires inline sequence; use the canvas-overlay path "
            "instead to apply features to the currently-loaded plasmid."
        )
    seq = parsed["fasta_seq"]
    total = len(seq)
    if total == 0:
        raise ValueError("GFF3 ##FASTA section is empty.")
    features = _gff3_features_to_biopython(parsed, total)
    seqid = parsed["seqid"]
    rec = SeqRecord(Seq(seq), id=seqid, name=seqid,
                     description=f"Imported from GFF3: {_P(path).name}",
                     features=features)
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = (
        "circular" if parsed["is_circular"] else "linear"
    )
    return rec


def _gff3_apply_to_loaded_record(record, path: str) -> int:
    """Feature-transfer mode for sequence-less GFF3 files.

    Reads features from `path` and appends them to ``record.features``
    in place. Returns the number of features added. Rejects rows whose
    coords don't fit within ``len(record.seq)`` (which the caller's
    sequence guards against tampering).

    Raises ValueError when the GFF3 file IS a complete record (has a
    ##FASTA section) — that case should use `_gff3_path_to_record`
    instead.
    """
    from pathlib import Path as _P
    ok, reason = _safe_file_size_check(
        _P(path), _GB_INGEST_MAX_BYTES, "GFF3",
    )
    if not ok:
        raise ValueError(reason or "GFF3 file rejected")
    text = _P(path).read_text(encoding="utf-8", errors="replace")
    parsed = _parse_gff3_text(text)
    if parsed["fasta_seq"] is not None:
        raise ValueError(
            "GFF3 file carries its own sequence (##FASTA). Use the "
            "standalone-import path so the imported plasmid isn't "
            "silently grafted onto the currently-loaded one."
        )
    total = len(record.seq)
    parsed_len = parsed.get("length")
    if parsed_len and parsed_len != total:
        raise ValueError(
            f"GFF3 declares length {parsed_len:,} but loaded plasmid "
            f"is {total:,} bp — coordinate frames don't match. Refusing "
            f"to apply features that would land at wrong positions."
        )
    new_feats = _gff3_features_to_biopython(parsed, total)
    record.features = list(record.features) + new_feats
    return len(new_feats)


def _export_gff_to_path(record, path) -> dict:
    """Write `record` to `path` as GFF3. Atomic write. Returns
    ``{"path", "bp", "features"}``.

    No round-trip verify (GFF3 has no canonical 1:1 reader in the
    standard library and Biopython's GFF support lives in BCBio.GFF
    which is an optional dep). The serialiser is deterministic and
    test-covered, so the same source record always produces the
    same output."""
    from pathlib import Path as _Path

    p = _Path(path).expanduser()
    text = _record_to_gff3(record)
    _atomic_write_text(p, text)
    _log.info(
        "Exported GFF3 to %s (%d bp, %d features)",
        p, len(record.seq),
        len([f for f in record.features if f.type != "source"]),
    )
    return {
        "path":     str(p),
        "bp":       len(record.seq),
        "features": len([f for f in record.features
                          if f.type != "source"]),
    }


def _export_fasta_to_path(name: str, sequence: str, path) -> dict:
    """Write `sequence` to `path` as a single-record FASTA. Atomic write.

    Returns `{"path", "bp", "name"}` on success. Raises:
      ValueError  — empty name or empty sequence.
      OSError     — filesystem failures (write, replace, fsync).

    The sequence is written on a single line (no hard-wrap at 80 chars);
    that matches what Biopython's default SeqIO writer emits for us
    elsewhere and keeps downstream `grep`/`awk` one-liners simple.
    """
    from pathlib import Path as _Path

    header = (name or "").strip()
    seq = (sequence or "").strip().upper()
    if not header:
        raise ValueError("FASTA export needs a non-empty record name.")
    if not seq:
        raise ValueError("FASTA export needs a non-empty sequence.")

    p = _Path(path).expanduser()
    _atomic_write_text(p, f">{header}\n{seq}\n")

    _log.info("Exported FASTA to %s (%s, %d bp)", p, header, len(seq))
    return {"path": str(p), "bp": len(seq), "name": header}


def _export_embl_to_path(record, path) -> dict:
    """Write `record` to `path` as EMBL flatfile via BioPython's SeqIO.

    EMBL is the European Nucleotide Archive's flatfile format — same
    feature-table model as GenBank, different text layout. Round-trip
    via SeqIO is straightforward; the writer preserves features,
    qualifiers, and circular topology. Atomic write.

    Returns ``{path, bp, features}`` on success. Raises:
      ValueError — record has no sequence.
      OSError    — filesystem failures (write, replace, fsync).
    """
    from pathlib import Path as _Path
    from io import StringIO
    from Bio import SeqIO

    if record is None or not getattr(record, "seq", None):
        raise ValueError("EMBL export needs a record with a sequence.")
    bp = len(record.seq)
    if bp == 0:
        raise ValueError("EMBL export needs a non-empty sequence.")

    # 2026-05-27 (audit-3 H6): route through `_normalize_for_genbank`
    # before the EMBL writer. Pre-fix EMBL export bypassed the
    # normalisation step that fills in molecule_type / topology /
    # accessions, so records constructed without these fields (GFF3
    # import, programmatically-built records) would raise mid-write
    # from inside Biopython's EMBL writer. EMBL's required-fields
    # set is a superset of GenBank's so the same normaliser covers it.
    normalized = _normalize_for_genbank(record)
    buf = StringIO()
    SeqIO.write([normalized], buf, "embl")
    text = buf.getvalue()

    p = _Path(path).expanduser()
    _atomic_write_text(p, text)

    n_feats = len([f for f in (record.features or [])
                   if f.type != "source"])
    _log.info("Exported EMBL to %s (%d bp, %d features)", p, bp, n_feats)
    return {"path": str(p), "bp": bp, "features": n_feats}


def _is_ab1_path(path) -> bool:
    try:
        suffix = getattr(path, "suffix", None)
        if suffix is None:
            suffix = Path(str(path)).suffix
    except Exception:
        return False
    return suffix.lower() in _AB1_EXTS


def _is_fastq_path(path) -> bool:
    try:
        suffix = getattr(path, "suffix", None)
        if suffix is None:
            suffix = Path(str(path)).suffix
    except Exception:
        return False
    return suffix.lower() in _FASTQ_EXTS


def _parse_fasta_single(path: str) -> tuple[str, str]:
    """Parse a FASTA file that must contain **exactly one** record and
    return ``(record_id, sequence)``.

    Multi-record FASTA files are rejected: the domesticator / parts bin
    flow only makes sense for a single part, so we surface a helpful
    error rather than silently picking the first record.

    Raises ``ValueError`` with a user-friendly message on any failure
    (read errors, zero records, multiple records, empty or non-IUPAC
    sequence). The sequence is upper-cased on success and validated
    against the IUPAC alphabet plus ``-``/``*``/``X`` for gap / stop /
    unknown.

    Size + symlink guard via ``_safe_file_size_check`` matches the
    `OpenFileModal` ingest pattern: rejects symlinks outright and
    refuses files larger than ``_BULK_IMPORT_MAX_BYTES`` (50 MB) so a
    multi-GB FASTA piped into the Domesticator's picker doesn't OOM
    the worker. Single-record FASTAs at the 50 MB ceiling already
    represent a 50 Mb sequence — bigger than this app supports.
    """
    from Bio import SeqIO
    ok, reason = _safe_file_size_check(
        Path(path), _BULK_IMPORT_MAX_BYTES, "FASTA",
    )
    if not ok:
        raise ValueError(reason or "FASTA file rejected by size check.")
    try:
        records = list(SeqIO.parse(path, "fasta"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Failed to read FASTA: {exc}") from exc
    if not records:
        raise ValueError("No FASTA records found in file.")
    if len(records) > 1:
        raise ValueError(
            f"Multi-sequence FASTA not supported ({len(records)} records "
            "found). Please provide a single-record FASTA."
        )
    rec = records[0]
    seq = str(rec.seq).upper()
    if not seq:
        raise ValueError("FASTA record has empty sequence.")
    valid = set("ACGTURYMKSWBDHVN-X*")
    bad = sorted(set(seq) - valid)
    if bad:
        raise ValueError(
            f"Non-IUPAC characters in sequence: {''.join(bad[:8])}"
        )
    return (rec.id or "fasta", seq)


def _fasta_path_to_record(path: str):
    """Parse a single-record FASTA at `path` into a `SeqRecord` ready
    to feed `_apply_record`.

    Defaults to `topology="linear"` because FASTA carries no topology
    annotation and assuming circular would mis-orient the user when
    they imported a chromosome chunk; the map's view-mode toggle
    flips to circular if the user knows otherwise. `molecule_type`
    defaults to ``"DNA"`` — protein FASTA is rejected upstream by
    `_parse_fasta_single`'s IUPAC check (which accepts ``X``/``*``
    but matches the codebase's plasmid-centric assumption that DNA
    is the right molecule).

    Raises whatever `_parse_fasta_single` raises for malformed input;
    callers are expected to surface the exception text via the modal
    status line.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    rec_id, seq = _parse_fasta_single(path)
    rec = SeqRecord(Seq(seq), id=rec_id, name=rec_id,
                     description=f"Imported from FASTA: {Path(path).name}")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = "linear"
    return rec


def _ab1_path_to_record(path: str):
    """Parse a Sanger trace `.ab1` / `.abi` file into a `SeqRecord`.

    BioPython's `SeqIO.parse(path, "abi")` base-calls the trace
    automatically and produces a single SeqRecord with quality scores
    in `letter_annotations["phred_quality"]`. We force `topology=linear`
    because Sanger reads are linear sequencing fragments — circular
    interpretation would mis-orient the user.

    Size + symlink check via `_safe_file_size_check` at entry — typical
    AB1 traces are 200 KB to a few MB, but a hostile / corrupted file
    could be arbitrarily large. The cap mirrors `_BULK_IMPORT_MAX_BYTES`
    (50 MB) used by every other import path.

    Raises ValueError on malformed traces / missing base-call channels
    or files exceeding the size cap.
    """
    import struct
    from Bio import SeqIO
    ok, reason = _safe_file_size_check(
        Path(path), _BULK_IMPORT_MAX_BYTES, "AB1",
    )
    if not ok:
        raise ValueError(reason or "AB1 trace rejected")
    try:
        records = list(SeqIO.parse(path, "abi"))
    except (OSError, ValueError, struct.error, KeyError) as exc:
        # Biopython's ABI reader struct.unpack's attacker-controlled directory
        # entries; a truncated/garbage trace can raise struct.error or KeyError.
        # Normalise all of them to the friendly ValueError (mirrors load_genbank).
        raise ValueError(f"could not parse AB1 trace: {exc}") from exc
    if not records:
        raise ValueError(f"no base-called sequence in {Path(path).name}")
    if len(records) > 1:
        raise ValueError(
            f"AB1 file contains {len(records)} records; expected one"
        )
    rec = records[0]
    # AB1 records leave `name` / `id` as the sample name, which is fine.
    # Force molecule_type/topology so downstream library/save paths see
    # the same fields as every other linear import.
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"]      = "linear"
    if not rec.description or rec.description == rec.id:
        rec.description = f"Imported from AB1 trace: {Path(path).name}"
    return rec


def _fastq_path_to_records(path: str) -> "list":
    """Parse a multi-read `.fastq` / `.fq` file into a list of
    `SeqRecord` objects. Each read becomes one record. Topology is
    forced to linear (reads are by definition linear fragments).

    Designed to feed `_handle_multi_record_fasta`'s collection-import
    flow — FASTQ becomes a new collection with one entry per read.
    Quality scores survive on `letter_annotations` but are not
    surfaced in the plasmid map (which is annotation-axis, not
    quality-axis).

    Two caps protect against pathological inputs:
      * File size capped at `_BULK_IMPORT_MAX_BYTES` (50 MB) via
        `_safe_file_size_check` — a multi-GB Illumina FASTQ would
        OOM `SeqIO.parse(...)`'s `list(...)` wrapper otherwise.
      * Read-count capped at `_FASTQ_MAX_READS` (1000) — even at the
        50 MB file cap a typical Illumina lane delivers ≫1000 short
        reads; the plasmid library wasn't designed to hold raw
        sequencing reads. Users with large datasets belong on the
        Plasmidsaurus tab (consensus-only) or in a dedicated aligner.

    Raises ValueError on malformed input, zero records, oversize, or
    read-count over cap.
    """
    from Bio import SeqIO
    ok, reason = _safe_file_size_check(
        Path(path), _BULK_IMPORT_MAX_BYTES, "FASTQ",
    )
    if not ok:
        raise ValueError(reason or "FASTQ file rejected")
    # Iterate + stop one past the cap instead of `list(SeqIO.parse(...))`: the
    # 50 MB file cap does NOT bound RAM, because 50 MB of short reads — each
    # carrying a Biopython phred-quality int list — materialises to ~1 GB
    # before a post-hoc len() check could fire. Breaking at the cap keeps peak
    # memory to _FASTQ_MAX_READS+1 records.
    records = []
    try:
        for rec in SeqIO.parse(path, "fastq"):
            records.append(rec)
            if len(records) > _FASTQ_MAX_READS:
                break
    except (OSError, ValueError) as exc:
        raise ValueError(f"could not parse FASTQ: {exc}") from exc
    if len(records) > _FASTQ_MAX_READS:
        raise ValueError(
            f"FASTQ contains more than {_FASTQ_MAX_READS:,} reads; cap is "
            f"{_FASTQ_MAX_READS:,} per file. Split the file or use a "
            f"dedicated read aligner — the plasmid library isn't "
            f"designed for raw sequencing reads."
        )
    if not records:
        raise ValueError(f"no reads in {Path(path).name}")
    for rec in records:
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"]      = "linear"
    return records



# ── CommercialSaaS .dna codec (Phase D, moved from hub) ─────────────────────
# Pure byte<->record codec for the SnapGene-style .dna binary format: the
# packet TLV reader/writer, history-XML extract/pack/inject, the from-scratch
# .dna emit (cookie/dna/features/notes packets), and the augment-from-packets
# recovery path. No data-dir writes, no network (the blob store, the loader,
# and fetch_genbank stay hub-side). Lazy struct/lzma/xml.etree/primer3 inside
# the fns exactly as in the hub.
# ── .dna packet I/O (low-level) ──────────────────────────────────────────────
#
# Binary format: each packet is a TLV — 1-byte type + 4-byte
# big-endian length + N bytes payload. The cookie packet (0x09)
# MUST come first; everything after is in implementation order.
# BioPython parses ~5 packet types (cookie, dna, primers, notes,
# features); the rest are silently dropped on read and impossible
# to write at all (BioPython has no .dna writer).
#
# These helpers go below the BioPython layer: they walk the raw
# byte stream, surface every packet (known or not), and let us
# round-trip files by splicing modified packets back into the
# original byte stream — preserving every packet we don't yet
# understand. Packet-type catalog (in progress) lives in
# `scripts/commercialsaas_inspect.py`'s ``KNOWN_PACKETS``.

_COMMERCIALSAAS_PACKET_HISTORY  = 0x07   # xz-compressed <HistoryTree> XML
_COMMERCIALSAAS_PACKET_COOKIE   = 0x09   # 8-byte format magic + 3 shorts
_COMMERCIALSAAS_HISTORY_MAX_XML = 32 * 1024 * 1024   # 32 MB hard cap on
                                                # decompressed payload —
                                                # protects against
                                                # decompression-bomb
                                                # crafted .dna files.

# Per-packet cap for the OTHER (uncompressed) XML packets — features (0x0A),
# primers (0x05), notes (0x06). ElementTree materialises a tree ~30x the raw
# XML byte size, so an uncapped ~100 MB flat-`<Feature>` packet (which fits
# inside the whole-file cap) would balloon to multiple GB and OOM the worker
# (the hosted demo droplet / small VMs fall over). Real feature/primer/notes
# XML is a few KB; 8 MB is wildly generous but bounds the amplification.
_COMMERCIALSAAS_PACKET_MAX_XML = 8 * 1024 * 1024


def _iter_commercialsaas_packets(data: bytes):
    """Yield ``(type_byte, length, payload_bytes)`` for every packet
    in a CommercialSaaS .dna byte stream. Raises ``ValueError`` on a
    declared-length overrun (audit-3 H2) OR a missing/wrong cookie
    packet at offset 0 (audit-3 M9, 2026-05-27).

    Empty input yields nothing (back-compat with test fixtures that
    pass `b""`). Non-empty input MUST start with a valid cookie
    packet — type byte 0x09, length 14, payload prefixed with the
    8-byte format magic. Pre-fix a malformed file lacking the cookie
    still parsed past the iterator and downstream consumers happily
    extracted features from junk — now the load refuses upfront so
    the user sees a clear "not a .dna file" error.

    Stops cleanly on EOF where a header doesn't fit (5-byte header;
    a partial trailer is treated as benign end-of-stream rather
    than a malformed-file error)."""
    import struct as _struct
    offset = 0
    n = len(data)
    # 2026-05-27 (audit-3 M9): cookie-packet validation gate.
    # Refuse files that don't start with the 0x09 cookie carrying the
    # 8-byte format magic. Empty input is allowed (yields nothing) so
    # back-compat tests keep working.
    if n > 0:
        if n < 5 + 8:
            raise ValueError(
                "CommercialSaaS .dna file too short for cookie packet "
                f"({n} bytes; need at least 13)"
            )
        if data[0] != _COMMERCIALSAAS_PACKET_COOKIE:
            raise ValueError(
                f"CommercialSaaS .dna file does not start with the "
                f"cookie packet (first byte 0x{data[0]:02X}, expected "
                f"0x{_COMMERCIALSAAS_PACKET_COOKIE:02X}). Not a valid "
                f".dna file."
            )
        if not data[5:5 + len(_COMMERCIALSAAS_COOKIE_MAGIC)].startswith(
                _COMMERCIALSAAS_COOKIE_MAGIC):
            raise ValueError(
                "CommercialSaaS .dna cookie packet payload doesn't "
                "carry the expected format magic. File is corrupt "
                "or not a valid .dna file."
            )
    while offset < n:
        if offset + 5 > n:
            # Trailing < 5 bytes: ambiguous whether benign EOF or
            # silent truncation. Treat as benign (the cookie + length
            # framing is designed to tolerate this — pre-2026-05-27
            # behaviour preserved here on purpose).
            return
        type_byte = data[offset]
        length = _struct.unpack(">I", data[offset + 1:offset + 5])[0]
        payload_start = offset + 5
        payload_end   = payload_start + length
        if payload_end > n:
            raise ValueError(
                f"CommercialSaaS packet length overrun at offset "
                f"{offset}: type=0x{type_byte:02X} declared {length} "
                f"bytes but only {n - payload_start} bytes remain. "
                f"File is truncated or corrupted; refusing to load "
                f"to avoid round-trip data loss."
            )
        yield (type_byte, length, data[payload_start:payload_end])
        offset = payload_end


def _build_commercialsaas_packet(type_byte: int, payload: bytes) -> bytes:
    """Serialise a single packet to bytes — type + 4-byte BE length +
    payload. Used by the writer + history-replace paths."""
    import struct as _struct
    if not (0 <= type_byte <= 0xFF):
        raise ValueError(f"packet type byte out of range: {type_byte!r}")
    if len(payload) > 0xFFFFFFFF:
        raise ValueError(f"payload too large for 32-bit length: "
                         f"{len(payload)} bytes")
    return bytes([type_byte]) + _struct.pack(">I", len(payload)) + payload


def _extract_commercialsaas_history_xml(data: bytes, *, packets=None) -> "str | None":
    """Find the 0x07 history packet in a .dna byte stream and return
    its decompressed XML (UTF-8 text) — or ``None`` if no history
    packet exists. Decompression is xz (LZMA); a malformed payload
    raises ``ValueError`` rather than returning a truncated string.

    Streaming decompression with a per-call cap defeats decompression
    bombs: a 10 MB compressed payload that expands to gigabytes is
    aborted at ``_COMMERCIALSAAS_HISTORY_MAX_XML`` rather than
    OOM-ing the worker. The decoder reads up to ``cap + 1`` bytes; if
    it fills, we know the input exceeded the cap.

    ``packets``: an already-materialised ``_iter_commercialsaas_packets(data)``
    list, shared by the file-load path so the buffer is walked once for all
    three .dna readers. ``None`` iterates ``data`` directly (every other caller).
    """
    import lzma as _lzma
    cap = _COMMERCIALSAAS_HISTORY_MAX_XML
    for type_byte, length, payload in (
            packets if packets is not None
            else _iter_commercialsaas_packets(data)):
        if type_byte != _COMMERCIALSAAS_PACKET_HISTORY:
            continue
        try:
            decoder = _lzma.LZMADecompressor()
            # Read up to cap + 1; if the decoder's not exhausted, we know
            # the input was bigger than the cap.
            decompressed = decoder.decompress(payload, max_length=cap + 1)
        except _lzma.LZMAError as exc:
            raise ValueError(
                f".dna history packet (0x07) is not valid xz: {exc}"
            ) from exc
        if len(decompressed) > cap or not decoder.eof:
            raise ValueError(
                f".dna history XML too large after decompression: "
                f">{cap:,} bytes (cap {cap:,})"
            )
        try:
            return decompressed.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(
                f".dna history XML is not valid UTF-8: {exc}"
            ) from exc
    return None


def _pack_commercialsaas_history_payload(xml_text: str) -> bytes:
    """Inverse of `_extract_commercialsaas_history_xml`'s decompression
    step: xz-compress the XML so it can be written as a 0x07
    payload. Uses the default LZMA preset (matches what real
    CommercialSaaS files appear to use). Caller is responsible for
    wrapping in `_build_commercialsaas_packet(0x07, …)`."""
    import lzma as _lzma
    encoded = xml_text.encode("utf-8")
    return _lzma.compress(encoded)


def _inject_commercialsaas_history(data: bytes, new_xml: "str | None") -> bytes:
    """Replace (or insert / remove) the 0x07 history packet in a
    .dna byte stream while preserving every other packet verbatim.
    Returns a new bytes object — input is not mutated.

    - ``new_xml`` non-empty string  → replace existing 0x07 (or
      insert a fresh one immediately after the cookie packet).
    - ``new_xml`` ``None`` or empty → strip the history packet
      entirely. Result is a valid .dna without history.

    Insertion position: immediately after the cookie packet. Real
    CommercialSaaS files position 0x07 in implementation-defined order
    among the metadata packets; we follow the convention of "as
    early as possible after cookie" since that's where CommercialSaaS
    seems to write it on round-trip."""
    new_packet: "bytes | None" = None
    if new_xml:
        payload = _pack_commercialsaas_history_payload(new_xml)
        new_packet = _build_commercialsaas_packet(_COMMERCIALSAAS_PACKET_HISTORY,
                                              payload)
    # Collect once so we can decide where to insert based on whether
    # an existing history packet is present. Two-pass keeps the
    # in-place-replace and the after-cookie-insert paths from
    # competing — the prior single-pass version emitted a new
    # packet at the cookie position AND skipped the existing 0x07
    # later, effectively moving the history slot forward in the
    # file (regression caught 2026-05-06).
    packets = list(_iter_commercialsaas_packets(data))
    has_history = any(t == _COMMERCIALSAAS_PACKET_HISTORY for t, _, _ in packets)
    out: list[bytes] = []
    emitted = False
    for type_byte, _length, payload in packets:
        if type_byte == _COMMERCIALSAAS_PACKET_HISTORY:
            # Replace in place (or skip if removing).
            if new_packet is not None and not emitted:
                out.append(new_packet)
                emitted = True
            continue
        out.append(_build_commercialsaas_packet(type_byte, payload))
        if (not has_history
                and type_byte == _COMMERCIALSAAS_PACKET_COOKIE
                and new_packet is not None
                and not emitted):
            # No existing history — insert the new one immediately
            # after the cookie.
            out.append(new_packet)
            emitted = True
    # Edge case: file had no cookie AND no history (malformed?);
    # append the new history at the end so it's not lost.
    if new_packet is not None and not emitted:
        out.append(new_packet)
    return b"".join(out)


# ── CommercialSaaS .dna packet writers (Phase 3 — from-scratch .dna emit) ──────────
#
# Builds a minimum-viable `.dna` file from a SeqRecord: cookie +
# DNA + features + notes + (optional) history. Other packet types
# real CommercialSaaS files carry — primers, alignments, custom enzymes,
# enzyme visibility, additional sequence properties, etc. — are
# omitted; CommercialSaaS Viewer treats most of them as optional metadata
# and fills sensible defaults on first open. Round-trip fidelity
# for those cosmetic packets is a Phase 5 hardening item.
#
# Bytes read by `_iter_commercialsaas_packets` should round-trip back to
# bytes written here (with re-derived DNA / features XML), so the
# writer is testable against the existing reader without needing
# CommercialSaaS Viewer in the loop.

_COMMERCIALSAAS_PACKET_DNA      = 0x00
_COMMERCIALSAAS_PACKET_PRIMERS  = 0x05
_COMMERCIALSAAS_PACKET_NOTES    = 0x06
_COMMERCIALSAAS_PACKET_ADDPROPS = 0x08
_COMMERCIALSAAS_PACKET_FEATURES = 0x0A

# Cookie payload: 8-byte format magic + 3 unsigned shorts
# (seqType, exportVersion, importVersion). Values cribbed from
# real outputs of the commercial editor (`00 01 00 0f 00 13` =
# seqType 1, exp 15, imp 19) so files we emit advertise the same
# versions the current commercial release writes — gives us the best
# chance of being accepted without "old format" warnings. The 8-byte
# magic is stored hex-encoded so the trademarked string doesn't appear
# verbatim in source.
_COMMERCIALSAAS_COOKIE_MAGIC      = bytes.fromhex("536e617047656e65")  # 8 bytes
_COMMERCIALSAAS_COOKIE_SEQ_TYPE   = 1
_COMMERCIALSAAS_COOKIE_EXP_VER    = 15
_COMMERCIALSAAS_COOKIE_IMP_VER    = 19


def _build_commercialsaas_cookie_packet() -> bytes:
    """Build the 0x09 cookie packet — 14 bytes: 8-byte format magic + 3
    unsigned shorts."""
    import struct as _struct
    payload = _COMMERCIALSAAS_COOKIE_MAGIC + _struct.pack(
        ">HHH", _COMMERCIALSAAS_COOKIE_SEQ_TYPE,
        _COMMERCIALSAAS_COOKIE_EXP_VER, _COMMERCIALSAAS_COOKIE_IMP_VER,
    )
    return _build_commercialsaas_packet(_COMMERCIALSAAS_PACKET_COOKIE, payload)


def _build_commercialsaas_dna_packet(seq: str, *, circular: bool) -> bytes:
    """Build the 0x00 DNA packet — 1-byte flags + N-byte ASCII
    sequence. Real CommercialSaaS files appear to use lowercase bases
    in the payload; we lowercase to match the convention. Flag
    bit 0x01 = circular; other bits cleared (their meaning is
    not fully documented and CommercialSaaS defaults them on read).

    2026-05-27 (audit-3 H5): validate the sequence up front instead
    of relying on a bare ``.encode("ASCII")`` strict crash. A stray
    unicode char (BOM, en-dash, full-width letter) used to raise
    ``UnicodeEncodeError`` mid-write after the cookie packet was
    already serialised — the user saw an opaque traceback and lost
    the save. Now raises a clear ValueError BEFORE any output is
    emitted, naming the offending character so the user can find it.
    """
    if not isinstance(seq, str):
        raise ValueError(
            f"DNA packet sequence must be str (got {type(seq).__name__})"
        )
    bad: list[str] = []
    for ch in seq:
        if ord(ch) > 127:
            bad.append(ch)
            if len(bad) >= 3:
                break
    if bad:
        raise ValueError(
            f"DNA packet sequence contains non-ASCII character(s) "
            f"{', '.join(repr(c) for c in bad)}"
            f"{' (truncated)' if len(bad) >= 3 else ''}"
            f" — strip whitespace/unicode artefacts before export."
        )
    flags = 0x01 if circular else 0x00
    payload = bytes([flags]) + seq.lower().encode("ascii", "strict")
    return _build_commercialsaas_packet(_COMMERCIALSAAS_PACKET_DNA, payload)


def _build_commercialsaas_features_packet_from_record(record) -> bytes:
    """Build the 0x0A features packet by serialising every non-source
    feature from the record into CommercialSaaS's XML schema. Features
    with `CompoundLocation` parts emit one `<Segment>` per part
    (e.g., wrap features get 2 segments; spliced CDSes get one per
    exon)."""
    import xml.etree.ElementTree as _ET
    real_feats = [f for f in record.features if f.type != "source"]
    root = _ET.Element("Features",
                          nextValidID=str(len(real_feats)))
    for i, feat in enumerate(real_feats):
        attrs = {
            "recentID": str(i),
            "name":     _commercialsaas_feat_name(feat),
            "type":     feat.type or "misc_feature",
            "allowSegmentOverlaps": "0",
            "consecutiveTranslationNumbering": "1",
        }
        # Strand → CommercialSaaS's `directionality` attribute.
        # Forward = "1", reverse = "2", omit for unknown / unstranded.
        strand = feat.location.strand
        if strand == 1:
            attrs["directionality"] = "1"
        elif strand == -1:
            attrs["directionality"] = "2"
        feat_el = _ET.SubElement(root, "Feature", attrs)
        # Color: derive from `_DEFAULT_TYPE_COLORS` so newly-written
        # features get a sensible default that matches what SpliceCraft
        # renders. CommercialSaaS's library-wide colour map differs slightly,
        # but it gracefully accepts any 6-digit hex.
        color = _DEFAULT_TYPE_COLORS.get(feat.type or "", "#a6acb3")
        # Segments: one per CompoundLocation part; one for simple.
        for part in _commercialsaas_iter_location_parts(feat.location):
            start_1based = int(part.start) + 1
            end_1based   = int(part.end)
            _ET.SubElement(feat_el, "Segment", {
                "range": f"{start_1based}-{end_1based}",
                "color": color,
                "type":  "standard",
            })
        # Qualifiers: skip `label` (already in the `name` attribute).
        # Skip `translation` for non-CDS features (it's CDS-only and
        # CommercialSaaS re-derives it). Otherwise emit each value as a
        # `<V text=>` (or `<V int=>` when the value is an integer).
        for qname, qvals in (feat.qualifiers or {}).items():
            if qname == "label":
                continue
            q_el = _ET.SubElement(feat_el, "Q", {"name": qname})
            for v in (qvals if isinstance(qvals, list) else [qvals]):
                if isinstance(v, bool):
                    # Treat as int for the few CommercialSaaS attrs that use
                    # 0/1; bool subclasses int so isinstance comes
                    # first.
                    _ET.SubElement(q_el, "V", {"int": str(int(v))})
                elif isinstance(v, int):
                    _ET.SubElement(q_el, "V", {"int": str(v)})
                else:
                    # Sweep #11 (2026-05-20): strip control chars from
                    # the text attribute before serialising. Pre-fix
                    # hostile qualifier values containing `\x01`,
                    # `\x02`, etc. produced XML that BioPython
                    # accepted on write but `_safe_xml_parse` would
                    # choke on during round-trip — the file's
                    # `_augment_dna_record_from_packets` then silently
                    # dropped the colour + label overrides, losing
                    # the user's feature-styling work. Mirrors
                    # `_commercialsaas_feat_name` (line ~8230) which
                    # already strips controls from the feature label.
                    sanitised = "".join(
                        c if (c >= " " or c in "\t\n")
                        else " "
                        for c in str(v)
                    )
                    _ET.SubElement(q_el, "V", {"text": sanitised})
    body = _ET.tostring(root, encoding="unicode")
    xml = '<?xml version="1.0"?>' + body
    return _build_commercialsaas_packet(_COMMERCIALSAAS_PACKET_FEATURES,
                                     xml.encode("utf-8"))


def _build_commercialsaas_notes_packet(record) -> bytes:
    """Build the 0x06 notes packet from the SeqRecord's metadata.
    Always sets <Type>Synthetic</Type> (matches what CommercialSaaS writes
    for newly-built constructions); fills <Created> + <LastModified>
    from today's date, and <CreatedBy> as 'SpliceCraft' so users can
    tell at a glance which file came from where."""
    import xml.etree.ElementTree as _ET
    from datetime import datetime as _dt
    now = _dt.now()
    root = _ET.Element("Notes")
    _ET.SubElement(root, "Type").text = "Synthetic"
    _ET.SubElement(root, "ConfirmedExperimentally").text = "0"
    created = _ET.SubElement(root, "Created", {
        "UTC": now.strftime("%H:%M:%S"),
    })
    created.text = now.strftime("%Y.%m.%d")
    modified = _ET.SubElement(root, "LastModified", {
        "UTC": now.strftime("%H:%M:%S"),
    })
    modified.text = now.strftime("%Y.%m.%d")
    _ET.SubElement(root, "CreatedBy").text = "SpliceCraft"
    # Description from the record (mapped to `<Comments>` per BioPython
    # parser convention).
    desc = getattr(record, "description", "") or ""
    if desc and desc != "<unknown description>":
        _ET.SubElement(root, "Comments").text = str(desc)
    body = _ET.tostring(root, encoding="unicode")
    return _build_commercialsaas_packet(_COMMERCIALSAAS_PACKET_NOTES,
                                     body.encode("utf-8"))


def _extract_commercialsaas_file_date(data: bytes, *, packets=None) -> "str | None":
    """Best-effort: pull the source file's own date (normalised to ISO
    ``YYYY-MM-DD``) from the CommercialSaaS Notes packet's ``<Created>`` (or
    ``<LastModified>``) element — the inverse of the ``<Created>`` that
    `_build_commercialsaas_notes_packet` writes. Lets an imported plasmid's
    top History entry show WHEN the source file was made (user request
    2026-06-10), since the construction-history nodes themselves carry no
    per-step dates. Returns None when absent / unparseable — never raises.

    ``packets``: shared materialised packet list (see
    `_extract_commercialsaas_history_xml`); ``None`` iterates ``data``."""
    try:
        import xml.etree.ElementTree as _ET
        for type_byte, _length, payload in (
                packets if packets is not None
                else _iter_commercialsaas_packets(data)):
            if type_byte != _COMMERCIALSAAS_PACKET_NOTES:
                continue
            if len(payload) > _COMMERCIALSAAS_PACKET_MAX_XML:
                _log.warning("commercialsaas file-date: notes packet too "
                             "large (%d bytes) — skipping", len(payload))
                return None
            try:
                # Route the untrusted Notes XML through `_safe_xml_parse`
                # (defangs billion-laughs / DOCTYPE) — the same hardening the
                # features / primers / history packet readers already use. A
                # raw `_ET.fromstring` here was the LONE .dna XML reader that a
                # crafted Notes packet could detonate as an entity-expansion
                # bomb. `_safe_xml_parse` takes str, and `ET.fromstring` rejects
                # a str carrying an `<?xml … encoding=…?>` declaration, so
                # transcode and strip a leading declaration first (its declared
                # encoding is moot once the bytes are decoded to unicode).
                _text = re.sub(r"^\s*<\?xml[^>]*\?>\s*", "",
                                payload.decode("utf-8", "replace"), count=1)
                root = _safe_xml_parse(_text)
            except (_ET.ParseError, ValueError):
                return None
            for tag in ("Created", "LastModified"):
                el = root.find(tag)
                text = (el.text or "").strip() if el is not None else ""
                if text:
                    iso = text.replace(".", "-")      # "2026.06.09" → ISO
                    if _history_human_dt(iso):         # validates it parses
                        return iso
            return None      # notes present but no usable date
    except Exception:
        _log.debug("commercialsaas file-date extract failed", exc_info=True)
    return None


def _commercialsaas_feat_name(feat) -> str:
    """Return the display name for a feature in the CommercialSaaS
    convention: prefer the `/label` qualifier, fall back to the
    feature's type. Strips control bytes for safety."""
    quals = feat.qualifiers or {}
    label = (quals.get("label") or quals.get("product") or [feat.type or "?"])
    name = str(label[0] if isinstance(label, list) else label)
    return _CONTROL_CHARS_RE.sub("", name)[:200] or "feature"


def _build_commercialsaas_addprops_packet_default() -> bytes:
    """Build a default 0x08 AdditionalSequenceProperties packet —
    Upstream/DownstreamStickiness="0" (= blunt) and Upstream/Downstream
    Modification=FivePrimePhosphorylated. Real CommercialSaaS files
    emit this even on circular plasmids where end-stickiness has no
    biological meaning; it's part of the editor's standard packet
    inventory and Viewer reads it for the "Sequence Properties"
    inspector. Pinning the 289-byte default matches all three FFE_*
    fixtures byte-for-byte.

    SpliceCraft doesn't currently model strand stickiness or end
    modifications, so a constant default is the right output —
    extending this when we gain a richer linear-fragment model is a
    future hardening item."""
    xml = (
        "<AdditionalSequenceProperties>"
        "<UpstreamStickiness>0</UpstreamStickiness>"
        "<DownstreamStickiness>0</DownstreamStickiness>"
        "<UpstreamModification>FivePrimePhosphorylated</UpstreamModification>"
        "<DownstreamModification>FivePrimePhosphorylated</DownstreamModification>"
        "</AdditionalSequenceProperties>"
    )
    return _build_commercialsaas_packet(_COMMERCIALSAAS_PACKET_ADDPROPS,
                                       xml.encode("utf-8"))


def _build_commercialsaas_primers_packet_default() -> bytes:
    """Build a minimum 0x05 Primers packet — root ``<Primers>`` with a
    single ``<HybridizationParams>`` child carrying the defaults the
    commercial editor writes on save (10 bp min continuous match, 1
    mismatch allowed, 40°C min Tm, 5'-end matching with 15 bp). No
    actual ``<Primer>`` entries — primer features (``primer_bind``)
    are still emitted via the 0x0A features packet.

    Real CommercialSaaS files always carry this packet even when no
    user-tracked primers exist — emitting it lets our from-scratch
    files mirror the editor's expected packet inventory and stops the
    Primers panel from defaulting to "(empty)" instead of the user's
    configured hybridization defaults. Three of the FFE_* test
    fixtures all carry exactly this 217-byte default; pinning it byte-
    for-byte gives us a regression target for ``CommercialSaaS Viewer``
    acceptance.
    """
    xml = (
        '<?xml version="1.0"?>'
        '<Primers nextValidID="0">'
        '<HybridizationParams '
        'minContinuousMatchLen="10" '
        'allowMismatch="1" '
        'minMeltingTemperature="40" '
        'showAdditionalFivePrimeMatches="1" '
        'minimumFivePrimeAnnealing="15"'
        '/>'
        '</Primers>'
    )
    return _build_commercialsaas_packet(_COMMERCIALSAAS_PACKET_PRIMERS,
                                       xml.encode("utf-8"))


def _commercialsaas_iter_location_parts(location):
    """Yield `(start, end)` simple parts for a feature location.
    Handles both `SimpleLocation` (single part) and `CompoundLocation`
    (multi-part; emits one part per sub-location)."""
    parts = getattr(location, "parts", None)
    if parts:
        for p in parts:
            yield p
    else:
        yield location


@_timed("op.write_commercialsaas_dna")
def _write_commercialsaas_dna_bytes(record, *,
                                 history_xml: "str | None" = None) -> bytes:
    """Return a complete `.dna` byte stream from a SeqRecord.

    Packet order: cookie → DNA → features → notes → (optional)
    history. CommercialSaaS tolerates implementation-defined order for
    non-cookie packets; this order matches what real CommercialSaaS
    output writes and what BioPython's parser handles smoothly.

    The result is round-trippable through `_iter_commercialsaas_packets`,
    re-parseable by BioPython's commercialsaas reader (sequence + features
    + topology + notes), and — if Phase 5 validation succeeds —
    accepted by CommercialSaaS Viewer. Until that validation is done,
    treat the writer as "expected to work; please test against
    CommercialSaaS Viewer and report rejections".
    """
    if record is None:
        raise ValueError("record is None")
    seq = str(getattr(record, "seq", "") or "")
    if not seq:
        raise ValueError("record has empty sequence")
    annotations = getattr(record, "annotations", None) or {}
    is_circ = (annotations.get("topology", "") or "").lower() == "circular"
    parts: list[bytes] = []
    parts.append(_build_commercialsaas_cookie_packet())
    parts.append(_build_commercialsaas_dna_packet(seq, circular=is_circ))
    parts.append(_build_commercialsaas_features_packet_from_record(record))
    parts.append(_build_commercialsaas_notes_packet(record))
    # Default 0x08 AdditionalSequenceProperties — strand stickiness +
    # end modifications. Real files carry this even on circular
    # plasmids where end-stickiness is meaningless; emitting the same
    # 289-byte default keeps the Sequence Properties inspector from
    # falling back to "(empty)" when Viewer reads our output.
    parts.append(_build_commercialsaas_addprops_packet_default())
    # Default 0x05 Primers packet — every real CommercialSaaS file
    # carries this with the same `HybridizationParams` defaults even
    # when no user primers are tracked. Emitting it ourselves keeps
    # the from-scratch output's packet inventory aligned with what the
    # editor produces, so Viewer can render the Primers panel cleanly
    # and never has to fall back to "(empty)" defaults.
    parts.append(_build_commercialsaas_primers_packet_default())
    if history_xml:
        payload = _pack_commercialsaas_history_payload(history_xml)
        parts.append(_build_commercialsaas_packet(
            _COMMERCIALSAAS_PACKET_HISTORY, payload))
    return b"".join(parts)


def _augment_dna_record_from_packets(
    rec, data: bytes, *, packets=None,
) -> list[dict]:
    """Recover info BioPython's ``.dna`` parser drops:
      * **per-feature colours** from the 0x0A Features packet
        (``<Segment color="#RRGGBB"/>`` attributes). BioPython parses
        ``<Feature>`` name/type/location but throws the colour away —
        SpliceCraft then falls back to its rotating ``_FEATURE_PALETTE``,
        which gives correct-but-unfamiliar colours that don't match
        what the user saw in the original editor. This helper stamps
        ``ApEinfo_revcolor`` + ``ApEinfo_fwdcolor`` qualifiers on every
        non-source feature so the colour-read path in
        ``PlasmidMap._parse`` picks them up.
      * **primer sequence stamps** on every ``primer_bind`` feature so
        the seq panel renders them with the full primer machinery
        (flap detection, weak-primer arrow, partial-binding tooltip).
        Derived from the bound region's bases — forward primers take
        the top-strand sequence directly; reverse primers take the
        reverse-complement. Skipped if BioPython already provided a
        ``primer_seq`` (defensive — future BioPython versions may
        decode the 0x05 packet themselves).
      * **standalone <Primer> entries** from the 0x05 Primers packet,
        when present. Most user-saved ``.dna`` files keep the 0x05
        packet at its empty default (just ``HybridizationParams``), but
        files that the user has run primer-design on inside the editor
        carry real entries here; we surface them into ``primers.json``
        so the user's primer library mirrors what they had in the
        source file.

    Mutates ``rec`` in place. Returns a list of primer dicts (the
    ``primers.json`` shape) — one per primer_bind feature plus one per
    standalone 0x05 ``<Primer>`` entry, with duplicates by sequence
    already collapsed within this call. The caller (``_apply_record``)
    dedupes against the existing primer DB before persisting.
    """
    import xml.etree.ElementTree as _ET

    # Local Tm calculator — primer3 if available, 2+4 fallback otherwise.
    # Captured once at the top of the augment so we don't pay the import
    # cost for every primer entry we build below. Imported primers get a
    # computed Tm so the Primer Library table renders the same way it
    # does for designed primers (it does `f"{tm:.1f}°C"` on the value).
    try:
        import primer3 as _primer3
        def _calc_tm(s: str) -> float:
            try:
                return float(_primer3.calc_tm(s))
            except Exception:
                # primer3 occasionally barfs on weird sequences (very
                # short, contains N, etc.); fall back to the 2+4 rule.
                # Log so a wave of import-time degenerate primers
                # surfaces as a diagnosable bundle entry instead of
                # silent mis-Tm on every imported primer.
                _log.exception(
                    "import _calc_tm: primer3.calc_tm fell back to "
                    "GC approximation for %d-mer", len(s))
                gc = sum(1 for c in s.upper() if c in "GC")
                at = sum(1 for c in s.upper() if c in "AT")
                return float(2 * at + 4 * gc)
    except ImportError:
        def _calc_tm(s: str) -> float:
            gc = sum(1 for c in s.upper() if c in "GC")
            at = sum(1 for c in s.upper() if c in "AT")
            return float(2 * at + 4 * gc)

    feature_colors: list[str] = []
    feature_names_from_xml: list[str] = []
    standalone_primers: list[dict] = []

    # Iterate the packet stream ONCE: the file-load path materialises
    # `_iter_commercialsaas_packets(data)` into `packets` and shares it across
    # this augment + the history-XML + file-date readers, collapsing three full
    # passes over the (capped) .dna buffer into one. `packets=None` (every other
    # caller) iterates `data` directly, exactly as before.
    for type_byte, _length, payload in (
            packets if packets is not None
            else _iter_commercialsaas_packets(data)):
        if type_byte == _COMMERCIALSAAS_PACKET_FEATURES:
            if len(payload) > _COMMERCIALSAAS_PACKET_MAX_XML:
                _log.warning("dna augment: 0x0A features packet too large "
                             "(%d bytes) — skipping to avoid XML-expansion OOM",
                             len(payload))
                continue
            try:
                root = _safe_xml_parse(payload.decode("utf-8"))
            except (_ET.ParseError, UnicodeDecodeError, ValueError):
                _log.warning("dna augment: 0x0A features packet parse failed")
                continue
            if root is None:
                continue
            for feat_el in root.findall(".//Feature"):
                seg = feat_el.find("Segment")
                feature_colors.append(
                    seg.get("color", "") if seg is not None else ""
                )
                # Capture the raw XML `name` attribute too. BioPython's
                # `.dna` parser has been observed to mangle whitespace
                # in feature names (GH #17, a user 2026-05-13:
                # spaces replaced with backslashes after import). We
                # pin the label to whatever the XML actually contains
                # so the user sees what their authoring tool wrote,
                # not whatever BioPython did on the way in. Strip the
                # control-char set we'd refuse to write anyway —
                # NUL / CR / LF would break a single-row sidebar
                # render — but SPACES + every other printable char
                # survive verbatim.
                xml_name = feat_el.get("name", "") or ""
                xml_name = _CONTROL_CHARS_RE.sub("", xml_name)[:200]
                feature_names_from_xml.append(xml_name)
        elif type_byte == _COMMERCIALSAAS_PACKET_PRIMERS:
            if len(payload) > _COMMERCIALSAAS_PACKET_MAX_XML:
                _log.warning("dna augment: 0x05 primers packet too large "
                             "(%d bytes) — skipping to avoid XML-expansion OOM",
                             len(payload))
                continue
            try:
                root = _safe_xml_parse(payload.decode("utf-8"))
            except (_ET.ParseError, UnicodeDecodeError, ValueError):
                _log.warning("dna augment: 0x05 primers packet parse failed")
                continue
            if root is None:
                continue
            today = _date.today().isoformat()
            for prim_el in root.findall(".//Primer"):
                pseq = (prim_el.get("sequence") or "").upper().replace("U", "T")
                if not pseq:
                    continue
                pname = (prim_el.get("name") or "").strip()
                if not pname:
                    pname = f"primer_{len(standalone_primers) + 1}"
                standalone_primers.append({
                    "name":        pname,
                    "sequence":    pseq,
                    "tm":          round(_calc_tm(pseq), 1),
                    "primer_type": "imported",
                    "source":      ".dna import",
                    "pos_start":   None,
                    "pos_end":     None,
                    "strand":      None,
                    "date":        today,
                    "status":      "Imported",
                })

    # Stamp colour qualifiers + override feature labels on features
    # by enumeration order. The 0x0A packet only carries the features
    # the editor itself created; any `source` row in the SeqRecord
    # comes from BioPython's LOCUS parsing (not the 0x0A packet) so
    # it doesn't consume a slot. Out-of-order or extra features are
    # tolerated — we stop when we run off the end of either list.
    #
    # Label override: BioPython's `.dna` parser has been observed to
    # mangle whitespace in feature names (GH #17 — user-typed
    # "Integration Seq" or "Lambda T0 Terminator" comes out with
    # backslashes inserted around the spaces). We pin
    # `qualifiers["label"]` to the raw XML name attribute so the
    # user sees what their authoring tool wrote, not whatever
    # BioPython produced. Skipped when the XML name is empty (some
    # third-party .dna writers omit the attribute), in which case
    # BioPython's parsed label survives untouched.
    feat_idx = 0
    for f in rec.features:
        if f.type == "source":
            continue
        # Color stamp (existing behaviour).
        if feat_idx < len(feature_colors):
            c = feature_colors[feat_idx]
            if c and isinstance(c, str):
                c = c.strip()
                # Defensive: only accept plausible CSS hex colours
                # so a malformed packet can't sneak arbitrary strings
                # into our qualifiers.
                if c.startswith("#") and len(c) in (4, 7):
                    f.qualifiers["ApEinfo_revcolor"] = [c]
                    f.qualifiers["ApEinfo_fwdcolor"] = [c]
        # Label override from raw XML.
        if feat_idx < len(feature_names_from_xml):
            xml_name = feature_names_from_xml[feat_idx]
            if xml_name:
                f.qualifiers["label"] = [xml_name]
        feat_idx += 1
        if (feat_idx >= len(feature_colors)
                and feat_idx >= len(feature_names_from_xml)):
            break

    # Build primer DB entries from the primer_bind features. Two
    # sources for the primer sequence:
    #   * `primer_seq` qualifier already on the feature (the round-trip
    #     case — splicecraft stamps this when it writes a .dna, and ApE
    #     uses the same convention). USE this verbatim so a primer
    #     carrying a 5' flap keeps its full length in the DB; deriving
    #     from the bound region would drop the flap.
    #   * Otherwise derive from the bound region (forward primers take
    #     the top strand directly; reverse primers take the RC).
    # Pre-2026-05-10 the code `continue`d on the qualifier-present
    # branch, which skipped the DB-entry append — so any .dna file
    # round-tripped through splicecraft (or exported from a tool that
    # stamps primer_seq) lost its primers from the imported DB.
    seq_str = str(rec.seq).upper() if getattr(rec, "seq", None) else ""
    n = len(seq_str)
    today = _date.today().isoformat()
    primer_bind_entries: list[dict] = []
    for f in rec.features:
        if f.type != "primer_bind":
            continue
        try:
            bounds = _feat_bounds(f, n)
        except (TypeError, ValueError, AttributeError):
            bounds = None
        if bounds is None:
            continue
        start, end, strand = bounds
        strand = strand or 1
        # Wrap-aware slice: origin-spanning primer_bind features have
        # `end < start` after `_feat_bounds` normalisation. Pre-fix the
        # raw `int(f.location.start)/.end` flattened the `CompoundLocation`
        # to (min, max) and `seq_str[0:n]` stamped the WHOLE plasmid as
        # primer_seq. Rare in practice (most primers don't cross the
        # origin) but biologically wrong when it happens.
        if end < start:
            if not (0 <= start < n and 0 <= end <= n):
                continue
            sliced_top = seq_str[start:] + seq_str[:end]
        else:
            if not (0 <= start < end <= n):
                continue
            sliced_top = seq_str[start:end]
        existing_pseq = f.qualifiers.get("primer_seq", [])
        if existing_pseq and isinstance(existing_pseq, list):
            bound_seq = _normalize_primer_seq(existing_pseq[0])
            if not bound_seq:
                bound_seq = sliced_top
                if strand < 0:
                    bound_seq = _rc(bound_seq)
        else:
            bound_seq = sliced_top
            if strand < 0:
                bound_seq = _rc(bound_seq)
            f.qualifiers["primer_seq"] = [bound_seq]
        labels = f.qualifiers.get("label", [])
        pname = str(labels[0]).strip() if labels else f"primer_{start}_{end}"
        primer_bind_entries.append({
            "name":        pname,
            "sequence":    bound_seq,
            "tm":          round(_calc_tm(bound_seq), 1),
            "primer_type": "imported",
            "source":      ".dna import",
            "pos_start":   start,
            "pos_end":     end,
            "strand":      strand,
            "date":        today,
            "status":      "Imported",
        })

    # Merge + dedupe by sequence (case-insensitive). Standalone 0x05
    # entries come first so their explicit ``<Primer name="...">`` wins
    # over an auto-generated ``primer_<start>_<end>`` name when the
    # same sequence appears in both places.
    merged: list[dict] = []
    seen: set[str] = set()
    for entry in standalone_primers + primer_bind_entries:
        key = (entry["sequence"] or "").upper()
        if not key or key in seen:
            continue
        merged.append(entry)
        seen.add(key)
    return merged


@_timed("op.parse_commercialsaas_history")
def _parse_commercialsaas_history(xml_text: str) -> "_CommercialSaaSHistoryNode | None":
    """Parse `<HistoryTree>` XML into a node tree. Returns the root
    `<Node>` (the result plasmid) or ``None`` if the XML is empty /
    has no nodes. Raises ``ValueError`` on malformed XML.

    Routes through `_safe_xml_parse` to defang billion-laughs / DOCTYPE
    entity-expansion attacks: .dna files come from external sources
    (collaborators, online repositories, scraped archives) and the
    history XML packet is the most attacker-controlled payload in the
    binary stream.
    """
    import xml.etree.ElementTree as _ET
    if not xml_text or not xml_text.strip():
        return None
    try:
        root = _safe_xml_parse(xml_text)
    except _ET.ParseError as exc:
        raise ValueError(f"Invalid .dna history XML: {exc}") from exc
    if root.tag != "HistoryTree":
        raise ValueError(
            f"Expected root <HistoryTree>, got <{root.tag}>"
        )
    nodes = root.findall("Node")
    if not nodes:
        return None
    if len(nodes) > 1:
        # CommercialSaaS's convention is one top-level Node per file. If we
        # see more, take the first and warn. 2026-05-27 (audit-3 M5):
        # stash the sibling raw elements on the wrapper so the serialise
        # path can re-emit them — pre-fix the siblings were dropped on
        # round-trip even though we held a reference, costing the user
        # any provenance trees the source file carried alongside the
        # primary tree.
        _log.warning(
            "CommercialSaaS history has %d top-level <Node> elements; "
            "expected 1. Using the first; preserving %d sibling(s) "
            "for round-trip.", len(nodes), len(nodes) - 1,
        )
    wrapper = _CommercialSaaSHistoryNode(nodes[0])
    wrapper._sibling_elements = list(nodes[1:])
    return wrapper



# ── zip / Plasmidsaurus member handling (Phase D, moved from hub) ───────────
# [SUB-plasmidsaurus] zip member safety (path-traversal + size caps) + the
# Plasmidsaurus run-zip / per-base-TSV parsers + gbk display-label helper. The
# members are returned as text/metadata; the CALLER parses records (no
# load_genbank dependency here). Lazy zipfile/stat/codecs/re inside the fns.
# ── Sequencing-data ingestion (Plasmidsaurus etc.) ────────────────────────────
#
# Plasmidsaurus delivers each run as a .zip with per-sample directories
# containing a consensus assembly (`*.gbk` / `*.gb`), raw reads, and
# QC plots. SpliceCraft's first-pass workflow:
#
#   1. User picks the .zip via `PlasmidsaurusZipModal` (file browser).
#   2. We list every `.gbk` / `.gb` member inside (size-capped to keep
#      the picker responsive against absurd zips).
#   3. User picks a member; we `_extract_gbk_member` it into memory and
#      parse via `_gb_text_to_record`.
#   4. User picks a target plasmid from the library — we run a
#      `_pairwise_align` between the two sequences and surface the
#      result on `AlignmentScreen`.
#
# Future: a Plasmidsaurus API key + per-account download lands here as
# a network sibling of the zip path. Same downstream alignment +
# visualisation pipeline.

# Per-zip + per-member size caps. Plasmidsaurus zips for a single run
# are typically <50 MB and the consensus assembly is <100 KB. Anything
# larger is suspicious; refuse rather than OOM the picker.
_PLASMIDSAURUS_ZIP_MAX_BYTES    = 500 * 1024 * 1024   # 500 MB whole zip
_PLASMIDSAURUS_MEMBER_MAX_BYTES = 50  * 1024 * 1024   # 50 MB per .gbk
_PLASMIDSAURUS_MAX_MEMBERS      = 2000                # listing cap
_GBK_EXTS = (".gbk", ".gb", ".genbank")


def _is_safe_zip_member_name(name: str) -> bool:
    """Return True if `name` is safe to surface in the picker AND safe
    to feed back into `zf.getinfo(name)`. Rejects:

      * absolute paths (`/etc/passwd.gbk`, `C:\\Users\\…`)
      * `..` segments (`../../escape.gbk`)
      * NUL bytes (terminal-display ANSI smuggling)
      * embedded ANSI escape codes (likewise)

    Stdlib `zipfile.open(info)` doesn't extract to a path, so a hostile
    name is currently bounded — but a future `extract()` / `extractall()`
    switch would silently flow into path traversal. This rejection
    fences off the whole class of bugs at the listing stage.
    """
    if not name or "\x00" in name:
        return False
    # ANSI escape (\x1b) — used in terminal-injection attacks
    if "\x1b" in name:
        return False
    # Any other C0 control byte / DEL / C1 block (CR, bell, backspace, the
    # 8-bit CSI 0x9b, …) mangles a DataTable row render and is never a
    # legitimate member name; plus lone surrogates (a non-UTF-8 filename
    # decoded with surrogateescape) that can't even UTF-8-encode for
    # display / logging. Closes the same class as commit e573159 did for logs.
    if any(ord(c) < 0x20 or 0x7f <= ord(c) <= 0x9f for c in name):
        return False
    try:
        name.encode("utf-8")
    except UnicodeEncodeError:
        return False
    # Windows absolute (C:, D:, etc.)
    if len(name) >= 2 and name[1] == ":" and name[0].isalpha():
        return False
    # Backslash separator (Windows paths in malicious zips)
    norm = name.replace("\\", "/")
    if norm.startswith("/"):
        return False
    parts = norm.split("/")
    if any(p in ("..", ".") for p in parts):
        return False
    return True


def _list_gbk_members_in_zip(zip_path: Path) -> "list[dict]":
    """Walk a .zip and return a list of dicts describing every .gbk /
    .gb / .genbank member. Each dict carries:

        ``{"name": <member name>, "size": <uncompressed bytes>}``

    Raises ValueError on unreadable / oversized / non-zip files.
    Cap-protected: refuses zips above `_PLASMIDSAURUS_ZIP_MAX_BYTES`,
    individual members above `_PLASMIDSAURUS_MEMBER_MAX_BYTES`, and
    listings beyond `_PLASMIDSAURUS_MAX_MEMBERS` to keep the picker
    snappy and resistant to malformed archives.

    Sweep #26 (2026-05-23): closes the TOCTOU window between the
    size check and the zip open by opening via ``os.open(path,
    O_RDONLY)`` (which dereferences the symlink ONCE) then
    ``fstat`` on the fd (immune to a concurrent path swap) and
    ``zipfile.ZipFile(fileobj=os.fdopen(...))``. Pre-sweep a hostile
    local process could swap the file between the path-based
    ``stat()`` and the path-based ``ZipFile(str(p))`` call to
    bypass the size cap.
    """
    import zipfile
    import stat as _stat
    p = Path(zip_path)
    try:
        fd = os.open(str(p), os.O_RDONLY)
    except FileNotFoundError as exc:
        raise ValueError(f"zip not found: {p}") from exc
    except OSError as exc:
        raise ValueError(f"could not open zip: {exc}") from exc
    try:
        st = os.fstat(fd)
        if not _stat.S_ISREG(st.st_mode):
            raise ValueError(f"not a regular file: {p}")
        size = st.st_size
        if size > _PLASMIDSAURUS_ZIP_MAX_BYTES:
            raise ValueError(
                f"zip too large ({size:,} bytes; cap "
                f"{_PLASMIDSAURUS_ZIP_MAX_BYTES:,})"
            )
        try:
            fobj = os.fdopen(fd, "rb")
            fd = -1   # ownership transferred to fobj
        except OSError as exc:
            raise ValueError(f"could not open zip: {exc}") from exc
        try:
            zf = zipfile.ZipFile(fobj, "r")
        except (zipfile.BadZipFile, OSError) as exc:
            fobj.close()
            raise ValueError(f"could not open zip: {exc}") from exc
    finally:
        if fd != -1:
            try:
                os.close(fd)
            except OSError:
                pass
    members: list[dict] = []
    try:
        for info in zf.infolist():
            if len(members) >= _PLASMIDSAURUS_MAX_MEMBERS:
                _log.warning(
                    "plasmidsaurus: zip listing truncated at %d members",
                    _PLASMIDSAURUS_MAX_MEMBERS,
                )
                break
            if info.is_dir():
                continue
            name = info.filename
            if not _is_safe_zip_member_name(name):
                _log.warning(
                    "plasmidsaurus: skipping unsafe zip member name %r",
                    name,
                )
                continue
            base = name.rsplit("/", 1)[-1]
            if not base or base.startswith("."):
                continue
            ext = ""
            for e in _GBK_EXTS:
                if name.lower().endswith(e):
                    ext = e
                    break
            if not ext:
                continue
            if info.file_size > _PLASMIDSAURUS_MEMBER_MAX_BYTES:
                _log.warning(
                    "plasmidsaurus: skipping oversized member %s (%d bytes)",
                    name, info.file_size,
                )
                continue
            members.append({
                "name": name,
                "size": int(info.file_size),
            })
    finally:
        # Sweep #26: explicit close of both the zipfile AND the
        # underlying fileobj. `ZipFile.close()` does NOT close a
        # caller-supplied fileobj (per stdlib docs), so without
        # this `fobj.close()` we'd leak the underlying fd.
        try:
            zf.close()
        except OSError:
            pass
        try:
            fobj.close()
        except OSError:
            pass
    members.sort(key=lambda m: _natural_sort_key(m["name"]))
    return members


def _extract_gbk_member(zip_path: Path, member_name: str) -> str:
    """Read a single .gbk member out of a zip and return its decoded
    text. Caller passes the result through `_gb_text_to_record` to get
    a SeqRecord. Raises ValueError on missing / oversized / unreadable
    member, on UTF-8 / latin-1 decode failure, or on an unsafe member
    name (path traversal / NUL / ANSI smuggling).

    Sweep #26 (2026-05-23): TOCTOU-safe open via ``os.open`` +
    ``fileobj`` — see ``_list_gbk_members_in_zip`` rationale.
    """
    import zipfile
    import stat as _stat
    if not _is_safe_zip_member_name(member_name):
        raise ValueError(f"unsafe zip member name: {member_name!r}")
    try:
        fd = os.open(str(zip_path), os.O_RDONLY)
    except OSError as exc:
        raise ValueError(f"could not open zip: {exc}") from exc
    try:
        if not _stat.S_ISREG(os.fstat(fd).st_mode):
            raise ValueError(f"not a regular file: {zip_path}")
        try:
            fobj = os.fdopen(fd, "rb")
            fd = -1
        except OSError as exc:
            raise ValueError(f"could not open zip: {exc}") from exc
    finally:
        if fd != -1:
            try:
                os.close(fd)
            except OSError:
                pass
    try:
        with zipfile.ZipFile(fobj, "r") as zf:
            try:
                # Sweep #9 (2026-05-19): route through the resolver
                # so Windows-built zips with backslash separators in
                # their central directory still match a member name
                # we stored with forward slashes.
                info = _zf_get_member_info(zf, member_name)
            except KeyError as exc:
                raise ValueError(
                    f"member not in zip: {member_name!r}"
                ) from exc
            if info.file_size > _PLASMIDSAURUS_MEMBER_MAX_BYTES:
                raise ValueError(
                    f"member too large ({info.file_size:,} bytes; cap "
                    f"{_PLASMIDSAURUS_MEMBER_MAX_BYTES:,})"
                )
            # Bounded streaming read: the central-directory `file_size`
            # is attacker-controlled. A malicious zip can claim 1 KB
            # but stream 5 GB of compressed-to-zero data and OOM the
            # process. Cap the actual read at MAX+1 and abort if we
            # exceeded — pre-fix `fh.read()` was unbounded and only
            # the (untrusted) header was checked.
            cap = _PLASMIDSAURUS_MEMBER_MAX_BYTES
            with zf.open(info, "r") as fh:
                raw = fh.read(cap + 1)
            if len(raw) > cap:
                raise ValueError(
                    f"member exceeded cap during decompression "
                    f"(claimed {info.file_size:,} bytes; cap {cap:,}) "
                    f"— possible zip-bomb"
                )
    except (zipfile.BadZipFile, OSError) as exc:
        raise ValueError(f"could not read zip: {exc}") from exc
    finally:
        try:
            fobj.close()
        except OSError:
            pass
    # GenBank is ASCII per the spec; fall back to latin-1 only if a
    # stray high-bit byte slipped in (some sequencer pipelines do).
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1", errors="replace")


# ── Plasmidsaurus run-zip structured parser ───────────────────────────────────
#
# A standard Plasmidsaurus results zip groups files by category and
# names each by `<run>_<sample-idx>_<sample-name>.<ext>`. Categories
# seen in the wild:
#
#   <run>_genbank-files/<base>.gbk          — annotated consensus
#   <run>_fasta-files/<base>.fasta          — bare consensus FASTA
#   <run>_summary-files/<base>.txt          — k-mer % + contamination %
#   <run>_per-base-data/<base>.tsv          — per-bp coverage / quality
#   <run>_histograms/<base>.png             — read-length histogram
#   <run>_coverage-plots/<base>.png         — coverage plot
#   <run>_interactive-map/<base>.html       — Bokeh interactive map
#   <run>_ab1-files/<base>.<n-of-m>.ab1     — Sanger trace
#   <run>_gel.png                            — run-level virtual gel
#
# The `_parse_plasmidsaurus_zip` walk groups every member under its
# sample's canonical base name, surfaces run-level files separately,
# and returns small dicts the Sequencing toolbar's sub-tabs can render
# without re-reading the zip per tab. The summary-file body and the
# per-base TSV header are streamed inline (both tiny — <500 KB) so
# the QC tab can render contamination / coverage stats without a
# second pass.

# Cap on the inline-cached summary-file text we keep on the SequencingScreen
# instance. Per the Plasmidsaurus summary files inspected (≤200 bytes
# typical), 4 KB is generous; refuse anything larger to keep the
# in-memory cache bounded even on malformed runs.
_PLASMIDSAURUS_SUMMARY_MAX_BYTES = 4 * 1024

# Cap on per-base TSV reads. Real Plasmidsaurus per-base data is
# ~50 bytes/row × plasmid bp (so a 10 kbp plasmid ships ~500 KB).
# 100 MB is generous for the largest plasmids we'd ever see (200 kbp
# ≈ 10 MB) while still bounding a hostile zip that claims a tiny
# `file_size` in the central directory but streams unbounded data on
# decompress. Sacred-invariant defence-in-depth: also enforced inside
# `_summarize_perbase_tsv` via a chunked read so a single bombing
# line without newlines can't OOM `io.TextIOWrapper`'s buffer.
_PLASMIDSAURUS_PERBASE_MAX_BYTES = 100 * 1024 * 1024


def _normalize_zip_member(name: str) -> str:
    """Return `name` with backslashes folded to forward slashes. Some
    Windows-built zips ship `category\\sample.ext` member names; the
    rest of this module assumes forward slashes."""
    return (name or "").replace("\\", "/")


def _zf_get_member_info(zf, name: str):
    """Resolve a zip member by name, transparently handling
    Windows-built zips whose central directory stores `\\`
    separators while our internal stored names use `/`.

    Tries `getinfo(name)` first (the common case where the name
    we stored is the same as the name in the zip). On KeyError,
    falls back to looking the member up by its normalised form —
    matches against `_normalize_zip_member(info.filename)` so we
    catch backslash-stored entries whose normalised form equals
    the requested name.

    Raises KeyError if no match. Sacred robustness invariant added
    in Sweep #9 (2026-05-19) — without this, a Windows-built
    Plasmidsaurus zip on a POSIX host silently degraded the
    Sequencing tab to error rows because every `zf.getinfo` call
    KeyError'd against a stored-forward-slash name."""
    try:
        return zf.getinfo(name)
    except KeyError:
        pass
    norm_target = _normalize_zip_member(name)
    for info in zf.infolist():
        if _normalize_zip_member(info.filename) == norm_target:
            return info
    raise KeyError(name)


def _parse_plasmidsaurus_zip(zip_path: Path) -> dict:
    """Inspect a Plasmidsaurus results zip and return a structured
    summary keyed by sample.

    Returns ``{"run_id", "samples": [{"base", "name", "gbk", "fasta",
    "summary", "summary_text", "perbase", "perbase_coverage",
    "histogram", "coverage_plot", "interactive_map", "ab1_files"}, ...],
    "run_files": [{"name", "size", "category"}, ...], "total_files",
    "total_size"}``.

    Each ``gbk``/``fasta``/``summary``/``perbase`` value is the zip
    member name (or None when absent). ``summary_text`` is a verbatim
    decode of the summary file body when it is below
    ``_PLASMIDSAURUS_SUMMARY_MAX_BYTES``; consumers parse contamination
    + k-mer percentages from this. ``perbase_coverage`` carries
    coverage stats {"mean", "min", "max", "n_pos", "above_20x"} derived
    from the per-base TSV's `reads_all` column. Both summary_text and
    perbase_coverage are best-effort — missing / unreadable members
    leave them as ``""`` and ``{}``.

    All size caps from `_list_gbk_members_in_zip` apply. Members with
    paths failing `_is_safe_zip_member_name` are silently skipped.

    Raises ValueError on the zip itself being missing / oversized /
    corrupt — same surface as `_list_gbk_members_in_zip`.
    """
    import zipfile
    import stat as _stat
    p = Path(zip_path)
    # TOCTOU-safe open (sweep #26 pattern): open the fd ONCE, fstat it (immune
    # to a concurrent path swap / symlink), then hand the fileobj to ZipFile —
    # mirrors `_list_gbk_members_in_zip` / `_extract_gbk_member` instead of the
    # path-based exists()/is_file()/stat() + ZipFile(str(p)) it used to use.
    try:
        fd = os.open(str(p), os.O_RDONLY)
    except FileNotFoundError as exc:
        raise ValueError(f"zip not found: {p}") from exc
    except OSError as exc:
        raise ValueError(f"could not open zip: {exc}") from exc
    fobj = None
    try:
        st = os.fstat(fd)
        if not _stat.S_ISREG(st.st_mode):
            raise ValueError(f"not a regular file: {p}")
        if st.st_size > _PLASMIDSAURUS_ZIP_MAX_BYTES:
            raise ValueError(
                f"zip too large ({st.st_size:,} bytes; cap "
                f"{_PLASMIDSAURUS_ZIP_MAX_BYTES:,})"
            )
        try:
            fobj = os.fdopen(fd, "rb")
            fd = -1   # ownership transferred to fobj
        except OSError as exc:
            raise ValueError(f"could not open zip: {exc}") from exc
        try:
            zf = zipfile.ZipFile(fobj, "r")
        except (zipfile.BadZipFile, OSError) as exc:
            fobj.close()
            raise ValueError(f"could not open zip: {exc}") from exc
    finally:
        if fd != -1:
            try:
                os.close(fd)
            except OSError:
                pass

    # Sweep #35 (2026-05-26): everything from here through the loop
    # is inside the try/finally that closes `zf`. Pre-fix, the
    # dict/tuple setup below sat in an unguarded window — currently
    # exception-free, but a future maintainer wiring in any call
    # that can raise (per-zip lookups, custom suffix tables, etc.)
    # would silently leak the ZipFile handle until GC ran. Moving
    # the setup under the same try keeps the closure guarantee
    # tight regardless of future changes.
    try:
        # Per-sample dict keyed by canonical base name (filename without
        # the category-folder prefix, minus extension).
        samples: dict[str, dict] = {}
        run_files: list[dict] = []
        total_files = 0
        total_size = 0
        # Run ID: parent folder name prefix shared by category dirs
        # (e.g., `34XK5N_genbank-files/...` → run_id=`34XK5N`).
        # Inferred by majority vote so a malformed zip with mixed
        # prefixes doesn't pick a stale one.
        prefix_votes: dict[str, int] = {}

        # Category folder → field-name on the sample dict. Anchored on
        # the underscore-separator pattern Plasmidsaurus uses so a
        # sample whose own filename happens to contain `_genbank-files`
        # doesn't false-positive. We probe per-suffix rather than full
        # prefix because the run-ID prefix varies per zip.
        _CATEGORY_SUFFIXES = (
            ("_genbank-files",      "gbk"),
            ("_fasta-files",        "fasta"),
            ("_summary-files",      "summary"),
            ("_per-base-data",      "perbase"),
            ("_histograms",         "histogram"),
            ("_coverage-plots",     "coverage_plot"),
            ("_interactive-map",    "interactive_map"),
            ("_ab1-files",          "ab1"),   # multiple per sample
        )
        for idx, info in enumerate(zf.infolist()):
            if idx >= _PLASMIDSAURUS_MAX_MEMBERS:
                # Parity with the listing walk's cap (line ~15627): a zip of
                # millions of tiny members would otherwise balloon the per-
                # sample dicts. Log the truncation rather than silently drop.
                _log.warning(
                    "plasmidsaurus: zip exceeds %d members — parsing only "
                    "the first %d", _PLASMIDSAURUS_MAX_MEMBERS,
                    _PLASMIDSAURUS_MAX_MEMBERS,
                )
                break
            if info.is_dir():
                continue
            raw_name = _normalize_zip_member(info.filename)
            if not _is_safe_zip_member_name(raw_name):
                _log.warning(
                    "plasmidsaurus: skipping unsafe zip member name %r",
                    raw_name,
                )
                continue
            total_files += 1
            total_size += int(info.file_size or 0)
            parts = raw_name.split("/")
            if len(parts) < 2:
                # Top-level file (e.g., `<run>_gel.png`, `README.txt`).
                base = parts[0]
                # Run-level: anchor the prefix vote on a typical
                # top-level naming, e.g. `34XK5N_gel.png` → `34XK5N`.
                stem = base.rsplit(".", 1)[0]
                if "_" in stem:
                    prefix_votes[stem.split("_", 1)[0]] = (
                        prefix_votes.get(stem.split("_", 1)[0], 0) + 1
                    )
                run_files.append({
                    "name": base,
                    "size": int(info.file_size or 0),
                    "category": "run",
                })
                continue
            folder = parts[0]
            # Match a category suffix anchored on `_<suffix>` at end.
            matched_field: "str | None" = None
            for suffix, field in _CATEGORY_SUFFIXES:
                if folder.endswith(suffix):
                    matched_field = field
                    # Vote: everything before the matched suffix is
                    # the run-ID candidate (e.g. `34XK5N` from
                    # `34XK5N_genbank-files`). A zip whose folders
                    # share one prefix wins by majority.
                    prefix_candidate = folder[: -len(suffix)]
                    if prefix_candidate:
                        prefix_votes[prefix_candidate] = (
                            prefix_votes.get(prefix_candidate, 0) + 1
                        )
                    break
            if matched_field is None:
                # Unknown category folder. If the file itself is a
                # .gbk / .gb / .genbank, treat it as a standalone
                # sample so test-fixture-style zips (no
                # `_genbank-files/` subfolder, just `sample_A/x.gbk`)
                # still surface every plasmid. Otherwise surface
                # under `run_files` so a custom Plasmidsaurus build
                # with an extra folder is still inspectable.
                leaf_raw = parts[-1]
                low = leaf_raw.lower()
                if any(low.endswith(ext) for ext in _GBK_EXTS):
                    matched_field = "gbk"
                else:
                    run_files.append({
                        "name": raw_name,
                        "size": int(info.file_size or 0),
                        "category": folder,
                    })
                    continue
            leaf = parts[-1]
            # Canonical sample base: filename minus dot-extensions.
            # `.1-of-2.ab1` collapses to the base+strip-of-tail logic
            # below so all ab1s for a sample group under the same key.
            base = leaf
            for _ in range(4):  # at most a couple of dot segments
                if "." not in base:
                    break
                stem, _, ext = base.rpartition(".")
                if not stem:
                    break
                # Stop on the canonical biology extensions; everything
                # past `.ab1`, `.gbk`, `.fasta`, `.tsv`, `.txt`,
                # `.png`, `.html` is part of the sample name.
                if ext.lower() in (
                    "gbk", "gb", "genbank", "fasta", "fa", "tsv", "txt",
                    "png", "jpg", "jpeg", "html", "ab1",
                ):
                    base = stem
                    continue
                # `.1-of-2` style sub-suffix on ab1 — strip too so all
                # trace files for a sample collapse to one entry.
                if matched_field == "ab1" and "-of-" in ext:
                    base = stem
                    continue
                break
            entry = samples.setdefault(base, {
                "base": base,
                "name": base,
                "gbk": None, "fasta": None, "summary": None,
                "perbase": None, "histogram": None,
                "coverage_plot": None, "interactive_map": None,
                "ab1_files": [],
                "summary_text": "", "perbase_coverage": {},
            })
            if matched_field == "ab1":
                entry["ab1_files"].append(raw_name)
            else:
                entry[matched_field] = raw_name

        # Stream the small text bodies inline so the QC sub-tab can
        # render contamination / coverage without re-opening the zip.
        for sample in samples.values():
            sm = sample.get("summary")
            if sm:
                try:
                    info = _zf_get_member_info(zf, sm)
                    if info.file_size <= _PLASMIDSAURUS_SUMMARY_MAX_BYTES:
                        with zf.open(info, "r") as fh:
                            data = fh.read(
                                _PLASMIDSAURUS_SUMMARY_MAX_BYTES + 1,
                            )
                        if len(data) <= _PLASMIDSAURUS_SUMMARY_MAX_BYTES:
                            try:
                                sample["summary_text"] = data.decode(
                                    "utf-8",
                                )
                            except UnicodeDecodeError:
                                sample["summary_text"] = data.decode(
                                    "latin-1", errors="replace",
                                )
                except (KeyError, OSError, zipfile.BadZipFile):
                    _log.exception(
                        "plasmidsaurus: failed to read summary for %s",
                        sample.get("base", "?"),
                    )
            pb = sample.get("perbase")
            if pb:
                try:
                    info = _zf_get_member_info(zf, pb)
                    # Per-base TSV can run multi-MB on large plasmids;
                    # stream-decode rather than load fully into memory.
                    # Two-layer cap: refuse outright on a central-directory
                    # `file_size` that already overshoots the limit, then
                    # pass the cap into the streamer so a hostile zip with
                    # an under-claimed size that decompresses into GB still
                    # bails before OOM.
                    if info.file_size > _PLASMIDSAURUS_PERBASE_MAX_BYTES:
                        _log.warning(
                            "plasmidsaurus: skipping oversize perbase %s "
                            "(%d bytes; cap %d)",
                            pb, info.file_size,
                            _PLASMIDSAURUS_PERBASE_MAX_BYTES,
                        )
                    else:
                        with zf.open(info, "r") as fh:
                            sample["perbase_coverage"] = (
                                _summarize_perbase_tsv(
                                    fh,
                                    max_bytes=(
                                        _PLASMIDSAURUS_PERBASE_MAX_BYTES
                                    ),
                                )
                            )
                except (KeyError, OSError, zipfile.BadZipFile):
                    _log.exception(
                        "plasmidsaurus: failed to read perbase for %s",
                        sample.get("base", "?"),
                    )
    finally:
        zf.close()
        # ZipFile(fileobj) does NOT close the caller-supplied fileobj (stdlib
        # docs) — close it explicitly so the underlying fd isn't leaked.
        try:
            if fobj is not None:
                fobj.close()
        except OSError:
            pass

    # Resolve the run ID by majority vote on the prefixes we saw.
    run_id = ""
    if prefix_votes:
        run_id = max(prefix_votes.items(), key=lambda kv: kv[1])[0]

    # Natural-sort samples by name so DEMO32 < DEMO34 < ... < DEMO310.
    sample_list = sorted(
        samples.values(),
        key=lambda s: _natural_sort_key(s.get("name") or ""),
    )
    return {
        "run_id": run_id,
        "samples": sample_list,
        "run_files": run_files,
        "total_files": total_files,
        "total_size": total_size,
    }


# Display-name cap for API-imported samples (mirrors the hub-side
# `_BULK_IMPORT_MAX_NAME_LEN`; redeclared here because that const lives in the
# hub L7 and fileio can't import upward).
_PLASMIDSAURUS_SAMPLE_NAME_MAX = 256


def _plasmidsaurus_zip_to_entries(
        zip_path: Path, *,
        run_id: str = "") -> "tuple[list[dict], list[str]]":
    """Turn a downloaded Plasmidsaurus results zip into ready-to-save library
    entry dicts — the HEADLESS import behind the ``download-plasmidsaurus``
    agent endpoint. (The GUI fetch delegates to the interactive Sequencing
    importer instead, so it keeps alignment + the collision modal.)

    Lists every .gbk/.gb/.genbank member, extracts it, parses to a record,
    and builds the SAME entry shape ``_bulk_align_worker``'s add-path produces
    (``name/id/size/n_feats/source/added/gb_text/status``), tagging ``source``
    as ``plasmidsaurus:<run_id>:<sample>``. Per-member failures become warning
    strings rather than aborting the whole import. Returns
    ``(entries, warnings)``. Propagates ValueError only when the zip itself is
    unreadable / oversized / not a zip (from ``_list_gbk_members_in_zip``)."""
    entries: "list[dict]" = []
    warnings: "list[str]" = []
    members = _list_gbk_members_in_zip(zip_path)   # ValueError on a bad zip
    for m in members:
        name = m.get("name", "")
        try:
            gb_text = _extract_gbk_member(zip_path, name)
            rec = _gb_text_to_record(gb_text)
            sample = Path(name).stem
            display = (sample or rec.name or rec.id or "plasmid").strip()
            display = "".join(
                c if c.isprintable() else "_" for c in display
            )[:_PLASMIDSAURUS_SAMPLE_NAME_MAX]
            entries.append({
                "name":    display,
                "id":      str(rec.id or display).strip(),
                "size":    len(rec.seq),
                "n_feats": len([f for f in rec.features
                                if f.type != "source"]),
                "source":  (f"plasmidsaurus:{run_id}:{sample}"
                            if run_id else f"plasmidsaurus:{sample}"),
                "added":   _date.today().isoformat(),
                "gb_text": gb_text,
                "status":  "",
            })
        except (OSError, ValueError, RuntimeError) as exc:
            warnings.append(f"{name}: {exc}")
            _log.warning(
                "Plasmidsaurus import: skipped %r (%s)", name, exc)
    return entries, warnings


# A single per-base TSV row is ~tens of chars; a run this long with no newline
# is malformed / a zip-bomb line. Bounds `pending` so the streaming split stays
# O(n) instead of O(n^2) re-splitting an ever-growing buffer every chunk.
_PERBASE_MAX_LINE_BYTES = 1 * 1024 * 1024


def _summarize_perbase_tsv(fh, *, max_bytes: "int | None" = None) -> dict:
    """Stream the per-base TSV file handle and return summary coverage
    stats: ``{"mean", "min", "max", "n_pos", "above_20x"}``.

    Assumes the Plasmidsaurus column layout (col 0=pos, col 1=ref,
    col 2=reads_all, ...). Defensive on header detection — header row
    is skipped if col 2 isn't an integer. Returns an empty dict on
    parse failure so the caller can degrade gracefully.

    Streams in 64 KB chunks with a hard byte cap (``max_bytes``). A
    zip-bomb that decompresses into a single multi-GB line without
    newlines would otherwise let ``io.TextIOWrapper`` buffer the whole
    line in memory before yielding it. ``max_bytes=None`` skips the
    cap for trusted callers (no current caller hits this path — the
    Plasmidsaurus parser always passes its cap).
    """
    total = 0
    n = 0
    lo = 10 ** 9
    hi = -1
    above_20x = 0
    import codecs
    chunk_size = 64 * 1024
    try:
        utf8_decoder = codecs.getincrementaldecoder("utf-8")(
            errors="replace",
        )
    except LookupError:
        # Should never happen on a stdlib build, but degrade rather
        # than raising — the caller treats {} as "couldn't summarise".
        return {}
    consumed = 0
    pending = ""
    truncated = False

    def _consume_line(line: str, header_seen_box: list[bool]) -> None:
        nonlocal n, total, lo, hi, above_20x
        cols = line.rstrip("\r").split("\t")
        if len(cols) < 3:
            return
        if not header_seen_box[0]:
            header_seen_box[0] = True
            # Header probe: skip iff col 2 isn't numeric. Sweep #9
            # (2026-05-19): use `float()` rather than `int()` for
            # the numeric test so a TSV that ships fractional
            # `reads_all` (e.g. a sub-sampled assay) doesn't
            # silently discard its first data row.
            try:
                float(cols[2])
            except (ValueError, OverflowError):
                return
        try:
            v = int(float(cols[2]))
        except (ValueError, OverflowError):
            # `int(float("1e400"))` / `int(float("inf"))` raise OverflowError,
            # not ValueError — without this a single poisoned coverage row would
            # escape to the outer `except Exception` and discard ALL QC stats.
            return
        n += 1
        total += v
        if v < lo:
            lo = v
        if v > hi:
            hi = v
        if v >= 20:
            above_20x += 1

    header_seen = [False]
    try:
        while True:
            try:
                chunk = fh.read(chunk_size)
            except (OSError, ValueError):
                _log.exception(
                    "plasmidsaurus: perbase TSV read failed mid-stream"
                )
                return {}
            if not chunk:
                # Tail flush — last line without trailing newline.
                if pending:
                    _consume_line(pending, header_seen)
                break
            consumed += len(chunk)
            if max_bytes is not None and consumed > max_bytes:
                truncated = True
            text = utf8_decoder.decode(chunk, final=truncated)
            pending += text
            if "\n" in text:
                lines = pending.split("\n")
                pending = lines.pop()  # last fragment may be partial
                for line in lines:
                    _consume_line(line, header_seen)
            elif len(pending) > _PERBASE_MAX_LINE_BYTES:
                # No newline yet and `pending` has grown past a sane line
                # length: bail rather than keep O(n^2)-growing + re-splitting it.
                _log.warning(
                    "plasmidsaurus: perbase TSV line exceeds %d bytes with no "
                    "newline - aborting summary", _PERBASE_MAX_LINE_BYTES,
                )
                return {}
            if truncated:
                # Flush any partial tail we have, then stop.
                if pending:
                    _consume_line(pending, header_seen)
                _log.warning(
                    "plasmidsaurus: perbase TSV exceeded %d-byte cap; "
                    "summary truncated", max_bytes,
                )
                break
    except Exception:
        _log.exception("plasmidsaurus: perbase TSV summarize failed")
        return {}
    if n == 0:
        return {}
    return {
        "mean":      total / n,
        "min":       lo if lo != 10 ** 9 else 0,
        "max":       hi if hi >= 0 else 0,
        "n_pos":     n,
        "above_20x": above_20x,
    }


def _parse_plasmidsaurus_summary(text: str) -> dict:
    """Parse the body of a Plasmidsaurus `<run>_<n>_<name>.txt` summary
    file. Returns ``{"kmer_moles_pct": float, "kmer_mass_pct": float,
    "contamination_pct": float, "contamination_source": str,
    "raw_text": str}``. Missing fields default to ``None``.

    Format observed in the wild::

            1-mer (%)  2-mer (%)
        moles       97.5        2.5
        mass        95.1        4.9


        *************************


        E. coli genomic contamination: 18.0%

    Both the column count (`<n>-mer`) and the contamination source
    line vary across runs; we extract the first numeric token after
    `moles`/`mass` (the canonical 1-mer percentages) and grep the
    contamination line for the % value + source organism.
    """
    out: dict = {
        "kmer_moles_pct":      None,
        "kmer_mass_pct":       None,
        "contamination_pct":   None,
        "contamination_source": "",
        "raw_text":             text,
    }
    if not text:
        return out
    import re as _re
    moles_match = _re.search(
        r"^\s*moles\s+([0-9.]+)", text, _re.MULTILINE,
    )
    mass_match = _re.search(
        r"^\s*mass\s+([0-9.]+)",  text, _re.MULTILINE,
    )
    if moles_match:
        try:
            out["kmer_moles_pct"] = float(moles_match.group(1))
        except (TypeError, ValueError):
            pass
    if mass_match:
        try:
            out["kmer_mass_pct"] = float(mass_match.group(1))
        except (TypeError, ValueError):
            pass
    contam_match = _re.search(
        r"([^\n]+?)\s*contamination:\s*([0-9.]+)%", text, _re.IGNORECASE,
    )
    if contam_match:
        out["contamination_source"] = contam_match.group(1).strip()
        try:
            out["contamination_pct"] = float(contam_match.group(2))
        except (TypeError, ValueError):
            pass
    return out


def _display_label_for_gbk(name: str) -> str:
    """INV-73 (2026-05-25): convert a Plasmidsaurus gbk basename
    (e.g. ``RUN42_1_DEMO34.gbk``) into a display-ready label
    (``DEMO34``) by stripping the run-id prefix, the file extension,
    and replacing any remaining underscores with spaces. User
    feedback: "no underscores in names in the TUI".

    Falls back to the basename minus extension if no run prefix
    is found, so non-Plasmidsaurus filenames still get clean
    display labels.
    """
    if not name:
        return ""
    leaf = name.rsplit("/", 1)[-1]
    for ext in (".gbk", ".gb", ".genbank", ".dna"):
        if leaf.lower().endswith(ext):
            leaf = leaf[: -len(ext)]
            break
    m = re.match(
        r"^[A-Z0-9]+_\d+_(.+)$", leaf, flags=re.IGNORECASE,
    )
    if m:
        leaf = m.group(1)
    return leaf.replace("_", " ").strip()



# ── .dna blob store + export (Phase D, moved from hub) ──────────────────────
# The per-entry .dna sidecar blob store (save/load/delete) under
# _state._DNA_ORIGINALS_DIR + the from-scratch .dna export. DATA-SAFETY: every
# write/delete goes through the persistence chokepoint
# (_refuse_unauthorized_write / _refuse_unauthorized_delete). Lazy os/tempfile/
# hashlib/Path inside the fns, as in the hub.
_DNA_SIDECAR_MAX_BYTES = 50 * 1024 * 1024   # 50 MB — same as bulk import
# Cap the sidecar basename so the full path stays under NTFS's 260-char
# default limit on reasonable installs. Raised to 200 from the implicit
# OS limit so a hostile / accidentally-pasted multi-KB entry_id can't
# trip ENAMETOOLONG at write time.
_DNA_SIDECAR_BASENAME_MAX = 200


def _dna_sidecar_path(entry_id: str) -> "Path":
    """Path for the sidecar file for `entry_id`. The id is sanitised
    via ``Path(...).name`` so a user-controlled id like ``../../etc/passwd``
    or ``foo/bar`` can't break out of the originals dir. NUL bytes are
    rejected (POSIX would raise ``ValueError`` on the join, but normalising
    here gives a stable sentinel). Empty / non-string ids fall back to a
    sentinel so the path is always under ``_DNA_ORIGINALS_DIR``.

    The basename is case-folded AND suffixed with an 8-char SHA-1 prefix
    of the original (un-folded) entry_id. Without the case-fold + hash:
    on case-insensitive filesystems (macOS APFS default, NTFS) two
    library entries ``pUC19`` and ``puc19`` collided on the same on-disk
    sidecar — silently overwriting each other's round-trip bytes, so
    exporting the older entry to `.dna` emitted the wrong molecule.
    The hash is computed on the raw id (pre-sanitisation) so two ids
    that collapse to the same `safe` after separator scrubbing
    (``a/b`` vs ``a_b``) also stay distinct. Total basename is capped
    so the resulting path stays under NTFS's 260-char total-path limit
    on reasonable installs (suffix `-<8>.dna` reserves 13 chars).
    """
    raw = str(entry_id) if entry_id else "_unknown_"
    cleaned = raw.replace("/", "_").replace("\\", "_").replace("\x00", "_")
    safe = Path(cleaned).name or "_unknown_"
    if safe in (".", "..") or not safe.strip("."):
        safe = "_unknown_"
    import hashlib as _hashlib
    digest = _hashlib.sha1(raw.encode("utf-8", errors="replace")).hexdigest()[:8]
    safe_lower = safe.lower()
    if _is_windows_reserved_stem(safe_lower):
        safe_lower = f"_{safe_lower}"
    max_safe = _DNA_SIDECAR_BASENAME_MAX - len(digest) - len(".dna") - 1
    if len(safe_lower) > max_safe:
        safe_lower = safe_lower[:max_safe]
    return _state._DNA_ORIGINALS_DIR / f"{safe_lower}-{digest}.dna"


def _dna_sidecar_legacy_path(entry_id: str) -> "Path":
    """Pre-0.8.9 sidecar path (case-sensitive, no hash discriminator).
    Migration scaffolding — `_load_dna_original` falls back to this
    when the canonical path is missing, and `_save_dna_original` /
    `_delete_dna_original` clean it up alongside the canonical write.
    Safe to remove a few releases after all-users-migrated."""
    raw = str(entry_id) if entry_id else "_unknown_"
    cleaned = raw.replace("/", "_").replace("\\", "_").replace("\x00", "_")
    safe = Path(cleaned).name or "_unknown_"
    if safe in (".", "..") or not safe.strip("."):
        safe = "_unknown_"
    return _state._DNA_ORIGINALS_DIR / f"{safe}.dna"


def _save_dna_original(entry_id: str, data: bytes) -> bool:
    """Write the original .dna bytes for `entry_id` into the sidecar
    dir. Returns True on success, False on any failure (logged but
    never raised — the sidecar is a nice-to-have, not a blocker
    for import). Atomic write via tempfile + os.replace, mirroring
    the JSON helpers' safety convention.

    Refuses to write if `data` exceeds the per-sidecar cap. Also
    refuses to overwrite an existing sidecar with EMPTY bytes (a
    bug in the caller would silently nuke the round-trip data
    otherwise)."""
    if not entry_id:
        return False
    if not data:
        _log.warning("dna sidecar: refusing to save empty bytes for %r",
                       entry_id)
        return False
    if len(data) > _DNA_SIDECAR_MAX_BYTES:
        _log.warning("dna sidecar: refusing to save %d bytes for %r "
                       "(cap %d)", len(data), entry_id,
                       _DNA_SIDECAR_MAX_BYTES)
        return False
    import tempfile as _tempfile
    import os as _os
    try:
        _state._DNA_ORIGINALS_DIR.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        _log.warning("dna sidecar: mkdir failed for %s: %s",
                       _state._DNA_ORIGINALS_DIR, exc)
        return False
    target = _dna_sidecar_path(entry_id)
    legacy = _dna_sidecar_legacy_path(entry_id)
    # [INV-66 extension] L2 chokepoint covers the .dna sidecar too. The
    # 2026-05-22 caught failure was a probe nuking JSON files; sidecars
    # are equally "the data is the product" since they carry the
    # CommercialSaaS round-trip bytes. An unsandboxed
    # `sc._save_dna_original("some_id", b"...")` would otherwise write
    # straight into the user's real ~/.local/share/splicecraft/
    # dna_originals/ and overwrite a valid sidecar by id collision.
    _refuse_unauthorized_write(target, "dna sidecar")
    try:
        fd, tmp_path = _tempfile.mkstemp(prefix=".tmp_",
                                            dir=str(_state._DNA_ORIGINALS_DIR))
        try:
            with _os.fdopen(fd, "wb") as f:
                f.write(data)
                f.flush()
                _os.fsync(f.fileno())
            _os.replace(tmp_path, target)
            _fsync_parent_dir(target)
        except Exception:
            try:
                _os.unlink(tmp_path)
            except OSError:
                pass
            raise
        if legacy != target:
            try:
                legacy.unlink(missing_ok=True)
            except OSError:
                pass
    except OSError as exc:
        _log.warning("dna sidecar: write failed for %r: %s",
                       entry_id, exc)
        return False
    return True


def _load_dna_original(entry_id: str) -> "bytes | None":
    """Return the sidecar bytes for `entry_id`, or None if missing /
    unreadable. Used by the export path to splice fresh history into
    the original .dna byte stream. Read size-capped via
    `_safe_file_size_check` so a hand-edited or filesystem-corrupted
    sidecar can't OOM the export path."""
    if not entry_id:
        return None
    target = _dna_sidecar_path(entry_id)
    legacy = _dna_sidecar_legacy_path(entry_id)
    try:
        if target.exists():
            ok, reason = _safe_file_size_check(
                target, _DNA_SIDECAR_MAX_BYTES, "dna sidecar",
            )
            if not ok:
                _log.warning("dna sidecar: %s", reason)
                return None
            return target.read_bytes()
        if legacy != target and legacy.exists():
            ok, reason = _safe_file_size_check(
                legacy, _DNA_SIDECAR_MAX_BYTES, "dna sidecar (legacy)",
            )
            if not ok:
                _log.warning("dna sidecar: %s", reason)
                return None
            return legacy.read_bytes()
        return None
    except OSError as exc:
        _log.warning("dna sidecar: read failed for %r: %s", entry_id, exc)
        return None


def _delete_dna_original(entry_id: str) -> bool:
    """Remove the sidecar for `entry_id` (both canonical and any
    pre-0.8.9 legacy path). Returns True if at least one file was
    removed; False on any failure (logged) or if neither existed.
    Idempotent — safe to call when removing a library entry whether or
    not it had a sidecar. Uses `unlink(missing_ok=True)` so a race
    where another process deletes the file between our `exists()` and
    `unlink()` doesn't false-negative on the otherwise-successful op."""
    if not entry_id:
        return False
    target = _dna_sidecar_path(entry_id)
    legacy = _dna_sidecar_legacy_path(entry_id)
    # [INV-66 extension + INV-75 sweep #27] L2 chokepoint: refuse to
    # unlink a user sidecar from an unsandboxed import. Same threat
    # model as `_save_dna_original` — a probe calling
    # `_delete_dna_original` would otherwise wipe the user's real
    # sidecar. Routed through `_refuse_unauthorized_delete` (sweep
    # #27); semantically equivalent to the prior write-gate call but
    # uses the delete-flavoured helper for grep-discoverability.
    _refuse_unauthorized_delete(target, "dna sidecar")
    paths = (target,) if target == legacy else (target, legacy)
    removed = False
    for p in paths:
        try:
            had = p.exists()
            p.unlink(missing_ok=True)
            if had:
                removed = True
        except OSError as exc:
            _log.warning("dna sidecar: delete failed for %r at %s: %s",
                           entry_id, p, exc)
    return removed


def _export_commercialsaas_dna(entry: dict, out_path: "Path | str") -> "str":
    """Write a `.dna` file for `entry`.

    Two paths:
      1. **Splice mode** — preferred when the entry has a sidecar
         (i.e., was originally imported from a CommercialSaaS `.dna`).
         Reads the sidecar bytes, splices in the entry's current
         `history_xml` via `_inject_commercialsaas_history`, writes.
         Preserves every CommercialSaaS-specific packet we don't yet
         understand (alignments, custom enzymes, etc.).
      2. **From-scratch mode** — when no sidecar exists. Parses the
         entry's `gb_text` into a SeqRecord and runs
         `_write_commercialsaas_dna_bytes` to emit a minimum-viable
         `.dna` (cookie + DNA + features + notes + optional
         history). Sequence + features are byte-correct; cosmetic
         packets (primers, alignments, enzyme-visibility) are
         omitted. CommercialSaaS Viewer is expected to fill defaults on
         first open.

    Returns the absolute output path on success. Raises
    ``ValueError`` when the entry has no usable id or no parseable
    GenBank text. Raises ``OSError`` on file-write failure."""
    from pathlib import Path as _Path
    eid = str(entry.get("id") or "")
    if not eid:
        raise ValueError("entry has no id; cannot resolve a path")
    original = _load_dna_original(eid)
    history_xml = entry.get("history_xml") or None
    if original is not None:
        # Splice mode: replace history packet, leave everything else
        # byte-identical.
        out_bytes = _inject_commercialsaas_history(original, history_xml)
    else:
        # From-scratch mode: rebuild from the GenBank text.
        gb_text = entry.get("gb_text") or ""
        if not gb_text:
            raise ValueError(
                f"entry {entry.get('name')!r} has no `gb_text` to "
                f"build from; cannot export"
            )
        try:
            record = _gb_text_to_record(gb_text)
        except Exception as exc:
            raise ValueError(
                f"entry {entry.get('name')!r}: GenBank text is not "
                f"parseable ({exc}); cannot build .dna"
            ) from exc
        out_bytes = _write_commercialsaas_dna_bytes(
            record, history_xml=history_xml,
        )
    out = _Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Atomic write via tempfile + replace (matches sidecar convention).
    import tempfile as _tempfile, os as _os
    fd, tmp_path = _tempfile.mkstemp(prefix=".tmp_", dir=str(out.parent))
    try:
        with _os.fdopen(fd, "wb") as f:
            f.write(out_bytes)
            f.flush()
            _os.fsync(f.fileno())
        _os.replace(tmp_path, out)
        _fsync_parent_dir(out)
    except Exception:
        try:
            _os.unlink(tmp_path)
        except OSError:
            pass
        raise
    return str(out.resolve())



# ── .dna / GenBank file loader (Phase D, moved from hub) ────────────────────
# The format-detecting single-file loader: sniff GenBank/EMBL/.dna by extension
# + content, parse via Biopython, then (for .dna) augment from the raw packets
# (colours/primers/history). _stamp_history_root_date stays hub-side (cloning
# L3) and is reached via _state._stamp_history_root_date_hook. fetch_genbank
# (network) stays hub-side pending the _DEMO_MODE migration.
# Local `.dna` augment-read cap. Mirrors the network caps (pitfall #20)
# for an untrusted-blob class that arrives via download/share, not just
# the network path. `_augment_dna_record_from_packets` walks raw bytes;
# without a cap a hostile `.dna` could OOM before any packet validation.
_DNA_AUGMENT_MAX_BYTES = 100 * 1024 * 1024  # 100 MB


# Third-party API contract: BioPython's SeqIO format identifier for the
# popular commercial plasmid editor's binary `.dna` format. Stored hex-encoded
# so the trademarked string never appears verbatim in our source.
_BIOPYTHON_DNA_FMT = bytes.fromhex("736e617067656e65").decode("ascii")


def _detect_plasmid_format(path: str) -> str:
    """Pick a Biopython SeqIO format key from a file path's extension.

    Supported:
      - GenBank        (.gb, .gbk, .genbank)       → "genbank"
      - EMBL           (.embl)                     → "embl"
      - popular commercial plasmid editor binary (.dna)
                                                  → ``_BIOPYTHON_DNA_FMT``

    Extensions are matched case-insensitively. Unknown extensions get
    a cheap content sniff: if the first non-empty line starts with
    "LOCUS " (GenBank) or "ID " (EMBL), we route accordingly; the
    8-byte CommercialSaaS cookie at offset 0 routes to `.dna`.
    Otherwise we still default to GenBank — that's the historical
    behaviour and the most common plasmid format; the parser will
    then raise a clear error.
    """
    from pathlib import Path
    suffix = Path(path).suffix.lower()
    if suffix == ".dna":
        return _BIOPYTHON_DNA_FMT
    if suffix in (".gb", ".gbk", ".genbank"):
        return "genbank"
    if suffix == ".embl":
        return "embl"
    # Unknown extension — peek at the first few bytes to give a
    # better error class than "GenBank parser said no". Reads are
    # capped so a 10 GB symlink can't OOM us before `load_genbank`
    # gets its own size guard.
    try:
        with open(path, "rb") as fh:
            head = fh.read(256)
    except OSError:
        return "genbank"
    if head.startswith(_COMMERCIALSAAS_COOKIE_MAGIC):
        return _BIOPYTHON_DNA_FMT
    try:
        text_head = head.decode("ascii", errors="replace").lstrip()
    except (UnicodeDecodeError, ValueError):
        text_head = ""
    if text_head.startswith("LOCUS "):
        return "genbank"
    if text_head.startswith("ID "):
        return "embl"
    return "genbank"


@_timed("op.load_genbank")
def load_genbank(path: str):
    """Load a plasmid file (GenBank .gb/.gbk or .dna). Returns
    SeqRecord.

    Despite the name (kept for backward compatibility), this also
    handles native `.dna` files (the popular commercial plasmid editor's
    binary format) via the BioPython binary parser keyed by
    ``_BIOPYTHON_DNA_FMT``. Dispatch is based on file extension.

    For `.dna` files, BioPython leaves `record.id` / `record.name`
    as `<unknown id>` / `<unknown name>` sentinels (the editor's own
    name/title metadata is not exposed through SeqIO). Backfill both
    from the file stem so the library and map title show something
    human-readable instead of `<unknown name>`.

    Size + symlink check via `_safe_file_size_check`: pre-fix a CLI
    invocation `splicecraft /path/to/huge.gb` (e.g. a 1 GB chromosome
    dump) reached `SeqIO.parse` directly and SeqIO eagerly parsed it
    into RAM. The cap matches the modal-flow cap.

    Raises ValueError with a user-friendly message if the file has no
    records, multiple records, is oversized, or is a symlink to a
    character device.
    """
    import struct
    from Bio import SeqIO
    from pathlib import Path as _P
    ok, reason = _safe_file_size_check(
        _P(path), _GB_INGEST_MAX_BYTES, "Plasmid",
    )
    if not ok:
        raise ValueError(reason or "Plasmid file is unsafe to load")
    fmt = _detect_plasmid_format(path)
    try:
        records = list(SeqIO.parse(path, fmt))
    except (ValueError, struct.error) as exc:
        # CommercialSaaS parser raises ValueError on most malformed files, but
        # Biopython's binary unpacking can also leak struct.error when a
        # .dna file has a truncated header / nonsense packet length.
        # Rewrap both as user-friendly ValueError so callers don't have
        # to special-case struct.error.
        if fmt == _BIOPYTHON_DNA_FMT:
            raise ValueError(
                f"Could not parse popular commercial plasmid editor "
                f"file {path}: {exc}. If this file was exported from an "
                f"older release, try re-exporting as .dna from a current "
                f"version of the editor."
            ) from exc
        raise
    rec = _pick_single_record(records, path)

    # CommercialSaaS (and occasionally minimally-annotated GenBank) records
    # leave id/name as Biopython sentinels. Fall back to the filename
    # so the UI has something meaningful to display.
    stem = _P(path).stem or "plasmid"
    # Sanitize: GenBank LOCUS names can't contain spaces (INSDC spec)
    # so the rec.id / rec.name fall back to an underscored, truncated
    # form. SACRED (2026-05-21): preserve the ORIGINAL stem with
    # spaces as `_tui_display_name` so panels show "my plasmid" while
    # GenBank serialisation still gets the LOCUS-safe "my_plasmid".
    # See feedback_no_underscores_in_names: file-being-imported wins.
    # `_sanitize_label` also strips control chars / lone surrogates a hostile
    # filename stem could otherwise bake into the LOCUS (rec.name/id) or the
    # display name; it preserves internal spaces, so the no-underscore rule
    # still holds.
    #
    # GenBank LOCUS names allow only letters/digits/underscore, and
    # Biopython's writer maps BOTH spaces AND hyphens to "_" (so an exported
    # `FRAG-CYP76AD6` round-trips with LOCUS `FRAG_CYP76AD6`). `_sanitize_label`
    # KEEPS hyphens, so replacing only spaces left hyphenated filenames
    # (`FRAG-…`/`PCR-…`, the petunia running-log) with `stem == safe_stem` —
    # the display-name override below never fired and the underscored LOCUS
    # won. Folding hyphens in too makes `safe_stem` a faithful LOCUS form, so
    # the override preserves the hyphenated filename as the display name (and a
    # backfilled rec.name/id can never carry a LOCUS-invalid hyphen). [INV-98]
    safe_stem = _sanitize_label(
        stem.replace(" ", "_").replace("-", "_"), max_len=16) or "plasmid"
    if not rec.id or rec.id.startswith("<unknown"):
        rec.id = safe_stem
    if not rec.name or rec.name.startswith("<unknown"):
        rec.name = safe_stem
    # SC-E: a `.gb` SpliceCraft itself exported carries the real display name
    # in a COMMENT marker — the MOST authoritative source (it's the construct's
    # true identity, regardless of how the file got named). Restore it FIRST so
    # it wins over the filename stem below. `SeqIO.parse` was used here (not
    # `_gb_text_to_record`), so this path needs its own restore call. [INV-98]
    _restore_display_name_from_comment(rec)
    # Display-name override: whenever the LOCUS-safe form differs from
    # the real filename stem — because the stem carried spaces or hyphens,
    # OR was longer than the 16-char LOCUS cap, OR held chars the sanitiser
    # dropped — keep the FULL stem visible in the UI. Pre-2026-06-22 this
    # only fired when the stem contained a space, so an underscored long
    # filename like `My_Long_Plasmid_v2` silently saved under the
    # 16-char truncation `My_Long_Plasmid_` (agent-API snag #15). The
    # full stem can itself carry underscores (the user's real filename) —
    # that's fine; the no-underscore rule only forbids LETTING the LOCUS
    # sanitisation mangle a spaced name. _apply_record / library_load
    # both read this attr. Gated on the COMMENT restore: a SpliceCraft-stamped
    # name already set above wins over the filename stem. [INV-98 / INV-15]
    if stem != safe_stem and not getattr(rec, "_tui_display_name", None):
        try:
            rec._tui_display_name = _sanitize_label(stem, max_len=200)
        except Exception:
            pass

    # `.dna` augmentation: re-read the file bytes to recover per-feature
    # colours + primer-bind sequence stamps + 0x05 primer-DB entries that
    # BioPython's `.dna` parser silently drops. Stashed on the record
    # as a transient attribute so `_apply_record` can flush the primer
    # entries into `primers.json` after the load completes. Failures here
    # are non-fatal — the record loads with palette colours and bare
    # primer features (the historical behaviour).
    if fmt == _BIOPYTHON_DNA_FMT:
        try:
            with open(path, "rb") as fh:
                _dna_bytes = fh.read(_DNA_AUGMENT_MAX_BYTES + 1)
            if len(_dna_bytes) > _DNA_AUGMENT_MAX_BYTES:
                _log.warning(
                    "dna augment skipped for %s: exceeds %d MB cap; "
                    "record loads with palette colours instead of source colours",
                    path, _DNA_AUGMENT_MAX_BYTES // (1024 * 1024),
                )
            else:
                # Walk the packet TLV stream ONCE and share it across all three
                # .dna readers (augment + history-XML + file-date) instead of
                # re-parsing the (capped) buffer three times. A malformed stream
                # raises ValueError here, caught by the outer handler exactly as
                # the augment's own raise was before.
                _dna_packets = list(_iter_commercialsaas_packets(_dna_bytes))
                extra_primers = _augment_dna_record_from_packets(
                    rec, _dna_bytes, packets=_dna_packets)
                if extra_primers:
                    rec._dna_primer_entries = extra_primers
                # Stash the construction-history packet + the raw bytes
                # on the record so a SINGLE-file open (→ `add_entry`)
                # round-trips history the way the bulk-import path
                # already does. Pre-fix only bulk import extracted the
                # 0x07 history + saved the `.dna` sidecar; opening one
                # `.dna` and re-exporting it silently dropped both the
                # history AND every CommercialSaaS-specific packet
                # (alignments, enzymes, …). `add_entry` consumes these.
                try:
                    hist_xml = _extract_commercialsaas_history_xml(
                        _dna_bytes, packets=_dna_packets)
                    if hist_xml:
                        # Stamp the file's OWN date (from its Notes packet)
                        # onto the dateless history root so the import's top
                        # History entry shows when the source file was made.
                        rec._dna_history_xml = _state._stamp_history_root_date_hook(
                            hist_xml,
                            _extract_commercialsaas_file_date(
                                _dna_bytes, packets=_dna_packets))
                except ValueError as exc:
                    _log.warning("dna history extract failed for %s: %s",
                                 path, exc)
                # Original bytes enable splice-mode export (byte-exact
                # round-trip preserving unknown packets). Only retained
                # when within the augment cap — see the skip branch above.
                rec._dna_original_bytes = _dna_bytes
        except (OSError, ValueError) as exc:
            _log.warning(
                "dna augment failed for %s (%s); record loads with "
                "palette colours instead of source colours",
                path, exc,
            )
    return rec



# ── NCBI Entrez fetch (Phase D, moved from hub) ─────────────────────────────
# Network record fetch by accession (GenBank + protein). The web-demo egress
# guard `_demo_block_network` STAYS hub-side (it reads _DEMO_MODE + has 5 other
# callers); reached here via _state._demo_block_network_hook (fail-closed).
# Entrez.tool stamps _state._sc_version (= hub __version__). Lazy Bio.Entrez/
# SeqIO/socket/io inside the fns. _NCBI_TIMEOUT_S now lives in splicecraft_net
# (shared with online BLAST); imported in the header.

# Cap on the GenBank-text response size from `Entrez.efetch`. The record-list
# API (`_NCBI_MAX_RESPONSE_BYTES = 4 MB`) is too tight for a real chromosome
# accession (e.g. an E. coli genome serialises to ~10 MB of GB text).
# 64 MB lets every legitimate plasmid / cosmid / BAC / small chromosome
# through while refusing a multi-GB pathological response (compromised /
# misconfigured server, MITM). Mirrors the `resp.read(MAX + 1)` + bail
# pattern from invariant #20.
_NCBI_GB_MAX_RESPONSE_BYTES = 64 * 1024 * 1024


@_timed("op.fetch_genbank")
def fetch_genbank(accession: str, email: str = "splicecraft@local"):
    """Fetch a GenBank record by accession from NCBI Entrez. Returns SeqRecord.

    Raises ValueError with a user-friendly message if NCBI returns no
    records (obsolete accession), multiple records, or a response larger
    than `_NCBI_GB_MAX_RESPONSE_BYTES`. A 30 s socket timeout is applied
    so a silent network stall surfaces as an error instead of pinning
    the worker thread forever; the size cap defends against a server
    that is reachable but misbehaving (a fast 30 s of GB-text could
    OOM the worker without it).
    """
    _state._demo_block_network_hook("NCBI fetch")
    # Defence-in-depth: clamp the accession to the safe charset here, not only
    # at the agent endpoint, so the 5 local UI/CLI callers are covered too. A
    # smuggled `;`/`/`/`..`/space can't reach Entrez's URL builder.
    _acc = _sanitize_accession(accession)
    if _acc is None:
        raise ValueError(f"invalid NCBI accession: {accession!r}")
    accession = _acc
    import io
    import socket
    import time as _time
    from Bio import Entrez, SeqIO
    Entrez.email = email
    # Sweep #27: NCBI's E-utilities policy asks every caller to identify
    # itself via `Entrez.tool` + `Entrez.email` so they can rate-limit
    # and contact problem clients. Pre-sweep we set only `email`; NCBI
    # could (and occasionally does) throttle or block unidentified
    # traffic. `tool` is a free-form ASCII identifier — using
    # `SpliceCraft/<version>` mirrors the User-Agent set on PyPI /
    # Kazusa fetches and makes our traffic identifiable in NCBI logs.
    Entrez.tool = f"SpliceCraft/{_state._sc_version}"
    # Biopython's Entrez._open adds its OWN retry loop (max_tries=3,
    # sleep_between_tries=15s) UNDER our 2-attempt loop below; left at the
    # default a stalling server could burn ~2×(3×30s socket + 2×15s) ≈ 240s,
    # not the ~60s our retry intends. We do our own retry, so cap Biopython
    # at one try per efetch call.
    Entrez.max_tries = 1
    # Serialise the process-global socket-timeout window (see _NCBI_FETCH_LOCK):
    # a concurrent fetch must not interleave set/restore and leave this socket
    # timeout-less → a pinned worker thread. Released in the `finally` below.
    _state._NCBI_FETCH_LOCK.acquire()
    prev_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(_NCBI_TIMEOUT_S)
    # One-retry pattern with 250 ms backoff to absorb transient
    # NCBI flakes (rate-limited 429, brief gateway 502, etc.). Two
    # attempts total — bounded by 2 × 30 s socket timeout + backoff.
    raw: "bytes | str | None" = None
    last_exc: "BaseException | None" = None
    try:
        for attempt in range(2):
            try:
                with Entrez.efetch(
                    db="nucleotide", id=accession,
                    rettype="gb", retmode="text"
                ) as handle:
                    raw = handle.read(_NCBI_GB_MAX_RESPONSE_BYTES + 1)
                break
            except (OSError, socket.timeout) as exc:
                last_exc = exc
                if attempt == 0:
                    _log.debug(
                        "NCBI fetch attempt 1 failed for %r (%s); retrying",
                        accession, exc,
                    )
                    _log_event(
                        "net.retry", endpoint="ncbi",
                        accession=accession, attempt=1,
                        exc_type=type(exc).__name__,
                    )
                    _time.sleep(0.25)
                    continue
                # Both attempts failed — re-raise so the caller's
                # outer handler turns this into a user-facing error.
                raise
    finally:
        socket.setdefaulttimeout(prev_timeout)
        _state._NCBI_FETCH_LOCK.release()
    if raw is None:
        # Defensive — if we got here without raw, something went
        # wrong that the loop didn't catch. Fall back to last_exc.
        raise last_exc if last_exc is not None else OSError(
            f"NCBI fetch failed for {accession!r}"
        )
    # Entrez handles return either bytes or str depending on Biopython
    # version; normalise to text for SeqIO.
    if isinstance(raw, bytes):
        raw_len = len(raw)
        text = raw.decode("utf-8", errors="replace")
    else:
        raw_len = len(raw.encode("utf-8", errors="replace"))
        text = raw
    if raw_len > _NCBI_GB_MAX_RESPONSE_BYTES:
        raise ValueError(
            f"NCBI accession {accession!r} returned more than "
            f"{_NCBI_GB_MAX_RESPONSE_BYTES // (1024 * 1024)} MB; refusing "
            f"to load. If this is a legitimate large record, fetch via "
            f"NCBI directly and use File > Open."
        )
    # An NCBI outage / rate-limit can return an HTML error page (or an
    # eutils <ERROR> document) with a 200 status instead of GenBank text;
    # SeqIO would then yield zero records and surface a cryptic "no
    # records" message. Detect the markup up front and say what's wrong.
    if text.lstrip()[:1] == "<":
        raise ValueError(
            f"NCBI returned an error page instead of a record for "
            f"{accession!r} — the service may be down or rate-limiting. "
            f"Wait a moment and try again."
        )
    records = list(SeqIO.parse(io.StringIO(text), "genbank"))
    rec = _pick_single_record(records, f"NCBI accession {accession!r}")
    # 2026-05-27 (audit-3 H4): verify the returned record actually
    # matches the requested accession. NCBI occasionally redirects an
    # obsolete accession to a different record; pre-fix we accepted
    # whatever came back and stamped it with the requested id, so a
    # user asking for L09137 could silently get a different plasmid.
    # Match against the record's id + accessions list (versioned or
    # not); fail loud on mismatch.
    requested = (accession or "").strip().upper()
    # Strip a `.N` version suffix from the user's request — accession
    # comparison is version-agnostic by NCBI convention.
    requested_base = requested.split(".", 1)[0] if requested else ""
    returned_ids: "list[str]" = []
    if rec.id and rec.id != "<unknown id>":
        returned_ids.append(rec.id.upper())
    for acc in (rec.annotations or {}).get("accessions") or []:
        if isinstance(acc, str) and acc:
            returned_ids.append(acc.upper())
    returned_bases = {rid.split(".", 1)[0] for rid in returned_ids}
    if (requested_base
            and requested_base not in returned_bases
            and requested not in {rid.upper() for rid in returned_ids}):
        raise ValueError(
            f"NCBI returned a record that doesn't match the requested "
            f"accession {accession!r}: got id={rec.id!r}, "
            f"accessions={(rec.annotations or {}).get('accessions') or []!r}. "
            f"The accession may be obsolete or redirected — fetch the "
            f"replacement directly from NCBI if intended."
        )
    return rec


def fetch_protein(accession: str, email: str = "splicecraft@local"):
    """Fetch a protein record by accession from NCBI Entrez (``db=protein``).
    Returns ``(description, aa_sequence)`` — the FASTA header text and the
    uppercase amino-acid sequence. Raises ValueError if NCBI returns no /
    multiple records, an oversized or empty response, or a record that
    doesn't match the requested accession. Mirrors `fetch_genbank`'s
    30 s timeout, one-retry, size-cap, AND demo-mode egress guard."""
    _state._demo_block_network_hook("NCBI protein fetch")
    # Defence-in-depth: clamp the accession to the safe charset here too (see
    # fetch_genbank) so every local caller is covered, not just the agent path.
    _acc = _sanitize_accession(accession)
    if _acc is None:
        raise ValueError(f"invalid NCBI accession: {accession!r}")
    accession = _acc
    import io
    import socket
    import time as _time
    from Bio import Entrez, SeqIO
    Entrez.email = email
    Entrez.tool = f"SpliceCraft/{_state._sc_version}"
    # Biopython's Entrez._open adds its OWN retry loop (max_tries=3,
    # sleep_between_tries=15s) UNDER our 2-attempt loop below; left at the
    # default a stalling server could burn ~2×(3×30s socket + 2×15s) ≈ 240s,
    # not the ~60s our retry intends. We do our own retry, so cap Biopython
    # at one try per efetch call.
    Entrez.max_tries = 1
    # Serialise the process-global socket-timeout window (see _NCBI_FETCH_LOCK).
    _state._NCBI_FETCH_LOCK.acquire()
    prev_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(_NCBI_TIMEOUT_S)
    raw: "bytes | str | None" = None
    last_exc: "BaseException | None" = None
    try:
        for attempt in range(2):
            try:
                with Entrez.efetch(
                    db="protein", id=accession,
                    rettype="fasta", retmode="text",
                ) as handle:
                    raw = handle.read(_NCBI_GB_MAX_RESPONSE_BYTES + 1)
                break
            except (OSError, socket.timeout) as exc:
                last_exc = exc
                if attempt == 0:
                    _log.debug("NCBI protein fetch attempt 1 failed for %r "
                               "(%s); retrying", accession, exc)
                    _time.sleep(0.25)
                    continue
                raise
    finally:
        socket.setdefaulttimeout(prev_timeout)
        _state._NCBI_FETCH_LOCK.release()
    if raw is None:
        raise last_exc if last_exc is not None else OSError(
            f"NCBI protein fetch failed for {accession!r}")
    if isinstance(raw, bytes):
        raw_len = len(raw)
        text = raw.decode("utf-8", errors="replace")
    else:
        raw_len = len(raw.encode("utf-8", errors="replace"))
        text = raw
    if raw_len > _NCBI_GB_MAX_RESPONSE_BYTES:
        raise ValueError(
            f"NCBI protein {accession!r} returned more than "
            f"{_NCBI_GB_MAX_RESPONSE_BYTES // (1024 * 1024)} MB; refusing "
            f"to load.")
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    rec = _pick_single_record(records, f"NCBI protein {accession!r}")
    seq = str(rec.seq).strip().upper()
    if not seq:
        raise ValueError(
            f"NCBI protein {accession!r} returned an empty sequence.")
    # Accession-match guard (mirrors fetch_genbank's H4 check): the FASTA
    # id should contain the requested accession (version-agnostic).
    requested = (accession or "").strip().upper()
    req_base = requested.split(".", 1)[0]
    rid = (rec.id or "").upper()
    if req_base and req_base not in rid and requested not in rid:
        raise ValueError(
            f"NCBI returned a protein that doesn't match the requested "
            f"accession {accession!r}: got id={rec.id!r}. The accession may "
            f"be obsolete or redirected.")
    return (rec.description or rec.id or accession).strip(), seq
