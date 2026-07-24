"""splicecraft_codon — the codon optimizer (Phase D, layer L2).

The mission-critical codon-optimization core, extracted from the hub: the
genetic-code table + the reverse-translation optimizer (_codon_optimize /
_codon_allocate / _codon_build_aa_map), CAI / GC metrics, forbidden-site
scrubbing (_codon_fix_sites / _codon_forbidden_sites), and codon-table TSV
parsing / search. The table ACCESSORS (_codon_tables_load/save/get) +
_CODON_GENETIC_CODE live in dataaccess; the enzyme catalog is reached via
_state._all_enzymes_hook (the restriction scanner's hook). The Kazusa /
genome-datasets network + the Rich usage-chart renderer stay hub-side
(separate concerns) for a later step. _MUT_* (mutagenesis codon usage) is
primer code, not here. Re-exported by the hub so sc.<name> + every call site
resolves unchanged.

ZERO-TOLERANCE subsystem (project_codon_optimizer_mission_critical) — guarded
by test_codon's property tests: round-trip, ACGT-only, in-frame, no internal
stop, no codon->AA mismatch.
"""
from __future__ import annotations

import re
from datetime import date as _date
from pathlib import Path

import splicecraft_state as _state
from splicecraft_logging import _log
from splicecraft_biology import _forbidden_hit_set, _iupac_pattern, _mut_revcomp
from splicecraft_util import _natural_sort_key, _sanitize_path
from splicecraft_net import (
    _HMM_DB_RETRY_BACKOFF_S, _NCBI_MAX_RESPONSE_BYTES, _build_hardened_url_opener,
    _sanitize_accession,
)
from splicecraft_fileio import _is_safe_zip_member_name
from splicecraft_dataaccess import (
    _CODON_GENETIC_CODE, _codon_tables_load, _codon_tables_save, _get_setting,
)


# Standard genetic code for CDS translation (no biopython dependency)
_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# NCBI genetic-code tables other than the standard (table 1) are resolved
# lazily from Biopython and cached. Table 1 returns the hand-rolled
# `_CODON_TABLE` above (fast path + identical behaviour for the
# overwhelmingly-common case).
_CODON_TABLE_BY_ID: "dict[int, dict[str, str]]" = {1: _CODON_TABLE}


def _codon_table_for(table_id: "int | None") -> "dict[str, str]":
    """Return the codon→AA map for an NCBI genetic-code id (the GenBank
    ``/transl_table`` qualifier).

    Table 1 (or None / falsy) is the standard code. Other ids — 2
    (vertebrate mito), 4 (Mycoplasma / Spiroplasma + mold/protozoan
    mito), 5 (invertebrate mito), 11 (bacterial/plastid), … — are built
    once from Biopython's ``CodonTable`` and cached. Reassigned codons
    (e.g. ``TGA`` = Trp in tables 2/4/5, ``AGR`` = stop in table 2) then
    translate CORRECTLY instead of rendering a wrong residue + a false
    premature-stop ⚠ on the map. An unknown / invalid id falls back to
    the standard code with a warning rather than crashing — a hand-edited
    ``/transl_table=99`` must not break the map. Stops map to ``"*"`` to
    match ``_CODON_TABLE`` exactly."""
    try:
        tid = int(table_id) if table_id else 1
    except (TypeError, ValueError):
        tid = 1
    cached = _CODON_TABLE_BY_ID.get(tid)
    if cached is not None:
        return cached
    try:
        from Bio.Data import CodonTable as _BioCodonTable
        ct = _BioCodonTable.unambiguous_dna_by_id[tid]
        m = dict(ct.forward_table)          # sense codons (ACGT only)
        for stop in ct.stop_codons:
            m[stop] = "*"
    except Exception:  # noqa: BLE001  (KeyError for bad id, ImportError, …)
        _log.warning(
            "Unknown genetic-code table %r; falling back to the standard "
            "code (table 1).", table_id,
        )
        m = _CODON_TABLE
    _CODON_TABLE_BY_ID[tid] = m
    return m


_STOP_CODONS = frozenset(("TAA", "TAG", "TGA"))

# Codon-selection strategies for `_codon_optimize`. Deliberately NOT offering
# "harmonize": true Angov-style harmonization needs the SOURCE organism's usage
# table to preserve the relative rare-codon positions that pace cotranslational
# folding, and we only consume the target table. Calling this "harmonize" is the
# exact mislabel that was removed in 2026-05-01; naming it back would undo that.
_CODON_MODE_FREQUENCY = "frequency"
_CODON_MODE_MAX_CAI   = "max_cai"
_CODON_MODES = (_CODON_MODE_FREQUENCY, _CODON_MODE_MAX_CAI)


_CODON_FIX_POS_RE = re.compile(r"codon (\d+)")


def _codon_fix_mutation_positions(mutations: list[str]) -> list[int]:
    """Given the string list returned by :func:`_codon_fix_sites`, return
    each mutation's 0-based codon-start nucleotide position in the insert.

    The mutation format is fixed by ``_codon_fix_sites`` — ``(codon N …)``
    where ``N`` is 1-based. A missing / malformed entry gets ``-1`` so
    callers can filter without raising.
    """
    out: list[int] = []
    for m in mutations:
        match = _CODON_FIX_POS_RE.search(m) if isinstance(m, str) else None
        if match:
            out.append((int(match.group(1)) - 1) * 3)
        else:
            out.append(-1)
    return out


# Forbidden sites for the optimizer's restriction-site fixer. Keys are the
# forward site only; the fixer adds the reverse complement automatically if
# the site is non-palindromic.
_CODON_DEFAULT_FORBIDDEN: dict[str, str] = {
    "BsaI":    "GGTCTC",
    "BsmBI":   "CGTCTC",
    "BbsI":    "GAAGAC",
    "EcoRI":   "GAATTC",
    "NdeI":    "CATATG",
    "XhoI":    "CTCGAG",
    "BamHI":   "GGATCC",
    "HindIII": "AAGCTT",
    "NcoI":    "CCATGG",
    "SalI":    "GTCGAC",
    "KpnI":    "GGTACC",
    "SacI":    "GAGCTC",
}


# Host-expression hazards that are NOT restriction sites. `_codon_fix_sites`
# has always accepted any IUPAC pattern, so these could be passed by hand —
# but `_codon_forbidden_sites` only resolves ENZYME NAMES out of the enzyme
# catalog, so there was no built-in route to a polyA signal or a splice donor
# and every caller had to hand-write the dict (agent-API field report,
# 2026-07-23). Degenerate codes are used where the consensus is degenerate;
# each entry is a sequence motif whose presence INSIDE a CDS causes trouble in
# that host, not a cut site.
_CODON_HOST_HAZARD_MOTIFS: "dict[str, dict[str, str]]" = {
    "bacterial": {
        # Internal ribosome-binding site — drives spurious internal
        # translation initiation and truncated product.
        "Shine-Dalgarno": "AGGAGG",
    },
    "plant": {
        # Canonical plant polyadenylation signal: premature transcript
        # cleavage part-way through the CDS.
        "plant polyA signal": "AATAAA",
        # 5' splice-donor core, GT(A/G)AG — cryptic intron excision.
        "cryptic splice donor": "GTRAG",
    },
    "mammalian": {
        "polyA signal":         "AATAAA",
        "polyA signal variant": "ATTAAA",
        "cryptic splice donor": "GTRAG",
    },
}


def _codon_hazard_motifs(hosts) -> "dict[str, str]":
    """Merge the built-in hazard motifs for each host in ``hosts`` into the
    ``{name: motif}`` shape `_codon_fix_sites` consumes.

    Raises ValueError naming the unknown host rather than silently returning
    a smaller set — a typo'd host that quietly scrubbed nothing would read as
    "my construct is clean" when nothing was ever checked.
    """
    if isinstance(hosts, str):
        hosts = [hosts]
    out: dict[str, str] = {}
    seen_motifs: set = set()
    for host in hosts or []:
        key = str(host or "").strip().lower()
        group = _CODON_HOST_HAZARD_MOTIFS.get(key)
        if group is None:
            raise ValueError(
                f"unknown host {host!r} — expected one of "
                f"{', '.join(repr(h) for h in sorted(_CODON_HOST_HAZARD_MOTIFS))}")
        # Dedup by the PATTERN, not the name: AATAAA is both the plant and the
        # mammalian polyA signal (under different names), and carrying it twice
        # would scrub the same site twice and list it twice in remaining_motifs.
        # First host named wins the label — the merge is host-order-stable.
        for name, motif in group.items():
            if motif in seen_motifs:
                continue
            seen_motifs.add(motif)
            out[name] = motif
    return out


def _codon_tables_add(name: str, taxid: str, raw: dict,
                      source: str = "user") -> dict:
    """Insert or replace a table in the registry. Dedup key is taxid when
    non-empty, else name. Returns the stored entry."""
    entries = _codon_tables_load()
    taxid = str(taxid or "").strip()
    name  = (name or "?").strip() or "?"
    entry = {
        "name":   name,
        "taxid":  taxid,
        "source": source,
        "added":  _date.today().isoformat(),
        "raw":    dict(raw),
    }
    def _same(e):
        if taxid and e.get("taxid") == taxid:
            return True
        if not taxid and e.get("name") == name:
            return True
        return False
    kept = [e for e in entries if not _same(e)]
    kept.append(entry)
    _codon_tables_save(kept)
    return entry


_CODON_TSV_MAX_CHARS = 1_000_000   # 64 codons; even a verbose export is KB


# Three-letter (and stop-alias) → one-letter for TSV amino-acid columns.
_CODON_TSV_AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "STOP": "*", "TER": "*", "END": "*",
}


def _parse_codon_tsv(text: str) -> dict:
    """Parse a tab/whitespace/comma-delimited codon-usage table into the
    in-memory ``{codon: (aa, count)}`` registry shape. Pure — no I/O.

    Each data row needs a 3-base codon and a usage count. An optional
    amino-acid column (1-letter, 3-letter, or ``*`` / ``Stop``) is
    validated against the standard code and otherwise derived from
    ``_CODON_GENETIC_CODE``. Rows whose first token isn't an ACGT/U codon
    (headers), blank lines, and ``#`` comments are skipped. Counts may be
    ints or floats (rounded); a fraction-only row (0..1) is scaled by
    1000 so relative preference survives. ``U`` is accepted and folded to
    ``T``.

    Raises ``ValueError`` with a readable message on a bad codon, a
    non-numeric / negative count, an AA/codon mismatch, a duplicate
    codon, or a file with zero usable rows — callers surface it in a
    status line rather than crashing.
    """
    if not isinstance(text, str):
        raise ValueError("codon table must be text")
    if len(text) > _CODON_TSV_MAX_CHARS:
        raise ValueError("codon table is too large (1 MB cap)")
    raw: "dict[str, tuple[str, int]]" = {}
    for lineno, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        toks = [t for t in re.split(r"[\s,]+", s) if t]
        if not toks:
            continue
        codon = toks[0].upper().replace("U", "T")
        if len(codon) != 3 or any(b not in "ACGT" for b in codon):
            continue   # header / non-codon row — skip silently
        if codon not in _CODON_GENETIC_CODE:
            raise ValueError(f"line {lineno}: {codon!r} is not a valid codon")
        if codon in raw:
            raise ValueError(f"line {lineno}: duplicate codon {codon!r}")
        expected_aa = _CODON_GENETIC_CODE[codon]
        aa_given: "str | None" = None
        numeric: list[str] = []
        for t in toks[1:]:
            tu = t.upper()
            if len(tu) == 1 and (tu in "ACDEFGHIKLMNPQRSTVWY" or tu in ("*", ".")):
                # Only `*` / `.` denote stop. Ambiguity codes — X (any AA),
                # B/Z/J, Sec `U`, Pyl `O` — must fall through and be ignored
                # (same as any other non-standard letter). The old substring
                # test `tu in "*.X"` matched `X` too, mis-mapping an "X" AA
                # column onto stop and spuriously rejecting the whole table
                # ("codon GCT encodes A but the file says *").
                aa_given = "*" if tu in ("*", ".") else tu
            elif tu in _CODON_TSV_AA3:
                aa_given = _CODON_TSV_AA3[tu]
            else:
                try:
                    float(t)
                    numeric.append(t)
                except ValueError:
                    pass   # stray non-numeric, non-AA token — ignore
        if aa_given is not None and aa_given != expected_aa:
            raise ValueError(
                f"line {lineno}: codon {codon!r} encodes {expected_aa!r} "
                f"but the file says {aa_given!r}"
            )
        if not numeric:
            raise ValueError(
                f"line {lineno}: no count/frequency column for {codon!r}"
            )
        ints = [t for t in numeric if float(t) == int(float(t))]
        count = (int(round(float(ints[-1]))) if ints
                 else int(round(float(numeric[-1]) * 1000)))
        if count < 0:
            raise ValueError(f"line {lineno}: negative count for {codon!r}")
        raw[codon] = (expected_aa, count)
    if not raw:
        raise ValueError(
            "no codon rows found — expected lines like 'GCT A 120' or "
            "'GCT 120' (codon then count)"
        )
    return raw


def _codon_name_parts(name: str) -> tuple[str, str]:
    """Return (genus, species) as lowercased tokens from an entry name.

    Genus = first whitespace-delimited token; species = second token (or "").
    Names like "E. coli K12" yield ("e.", "coli"); "Escherichia coli" yields
    ("escherichia", "coli"). No normalization between abbreviated and
    unabbreviated genera — users search what they see.
    """
    parts = (name or "").strip().split()
    genus   = parts[0].lower() if parts else ""
    species = parts[1].lower() if len(parts) > 1 else ""
    return genus, species


def _codon_search(query: str, entries: "list | None" = None) -> list[dict]:
    """Ranked search over taxid, genus, species, and full name.

    Rank 0: taxid exact match
    Rank 1: taxid prefix match
    Rank 2: genus prefix (first whitespace token of name)
    Rank 3: species prefix (second whitespace token)
    Rank 4: substring anywhere in the full name

    Results are sorted by (rank, name) so same-genus entries cluster and
    the strongest match wins. An empty/whitespace query returns the
    registry naturally sorted by display name (`Escherichia coli K12`
    before `Escherichia coli K-12 MG1655`) so the unfiltered table
    reads in human-friendly order rather than disk-insertion order.
    Same-rank ties under a non-empty query also use the natural-sort
    secondary key (e.g. `E. coli K12` before `K-12 MG1655`).
    """
    if entries is None:
        entries = _codon_tables_load()
    q = (query or "").strip().lower()
    if not q:
        # Empty query → natural-sort the registry by display name so
        # `Escherichia coli K12` lands near the top in stable order
        # regardless of the JSON-file insertion order.
        return sorted(
            entries,
            key=lambda e: _natural_sort_key(str(e.get("name") or "")),
        )
    ranked: list[tuple[int, tuple, dict]] = []
    for e in entries:
        name_lc  = str(e.get("name", "")).lower()
        taxid_lc = str(e.get("taxid", "")).lower()
        genus, species = _codon_name_parts(e.get("name", ""))
        if taxid_lc and taxid_lc == q:
            rank = 0
        elif taxid_lc and taxid_lc.startswith(q):
            rank = 1
        elif genus and genus.startswith(q):
            rank = 2
        elif species and species.startswith(q):
            rank = 3
        elif q in name_lc:
            rank = 4
        else:
            continue
        # Secondary key uses natural sort so same-rank matches come
        # back in human-friendly order (`E. coli K12` before `K-12 MG1655`).
        ranked.append((rank, _natural_sort_key(name_lc), e))
    ranked.sort(key=lambda t: (t[0], t[1]))
    return [t[2] for t in ranked]


def _codon_build_aa_map(raw: dict, *, genetic_code: "dict | None" = None) -> tuple[dict, dict]:
    """Given {codon: (aa, count)}, return (aa_codons, codon_frac) where
    aa_codons[aa] = [(codon, frac), ...] sorted by fraction descending, and
    codon_frac[codon] = fractional usage for its amino acid.

    Codon-integrity defense-in-depth: the amino-acid label is taken from
    `_CODON_GENETIC_CODE`, NOT the `aa` stored in `raw`. Every table loader
    already forces the canonical label, so for any real table this is a no-op
    (byte-identical output). But a hand-built in-memory table (a test, a
    `codon_table=` override) could mislabel a codon — and a mislabel here
    would silently emit a wrong-residue codon that `_mut_translate` reads back
    as a different protein, with no error raised, in a subsystem where a
    silent wrong-protein round-trip is catastrophic. Deriving the AA from the
    standard code closes that hole; codons outside the 64 ACGT keys (N / gaps)
    are skipped."""
    from collections import defaultdict
    # `genetic_code` defaults to the standard NCBI table 1; pass an alternate
    # (via `_codon_table_for`) so a reassigned codon (TGA=Trp, TAA=Gln, …) is
    # LABELLED with the residue the target host actually reads. Default keeps
    # every existing call byte-identical.
    gc = genetic_code if genetic_code is not None else _CODON_GENETIC_CODE
    aa_total: dict = defaultdict(int)
    codon_aa: dict = {}
    for codon, (_aa, count) in raw.items():
        canon = gc.get(codon)
        if canon is None:               # non-ACGT codon (N / gap) — not optimizable
            continue
        codon_aa[codon] = canon
        aa_total[canon] += int(count)
    codon_frac: dict = {}
    for codon, (_aa, count) in raw.items():
        canon = codon_aa.get(codon)
        if canon is None:
            continue
        total = aa_total.get(canon, 0) or 1
        codon_frac[codon] = count / total
    aa_codons: dict = defaultdict(list)
    for codon, canon in codon_aa.items():
        if canon == "*":
            continue
        aa_codons[canon].append((codon, codon_frac[codon]))
    for aa in aa_codons:
        aa_codons[aa].sort(key=lambda x: -x[1])
    return dict(aa_codons), codon_frac


def _codon_allocate(codons: list, n: int) -> list:
    """Pick ``n`` codons from ``[(codon, frac), ...]`` so the chosen multiset
    matches the frequency distribution as closely as possible (largest-
    remainder apportionment), then interleave them so no single codon
    clusters at the front. Deterministic: equal remainders break by input
    order, so the same ``(codons, n)`` always yields the same list.

    Shared by amino-acid positions AND frequency-matched stop codons. The
    return is GUARANTEED to be exactly length ``n`` (empty when ``n <= 0``):
    for a normalized table the apportionment already sums to ``n``, but a
    hand-rolled / truncated table whose fractions don't sum to 1 is padded
    with (or truncated to) the most-frequent codon so the caller can never
    receive a short list that would leave a position unfilled — i.e. NO
    rogue/empty codon ever reaches the output. Every codon in ``codons`` is
    synonymous (same residue, or all stops), so padding never changes the
    encoded amino acid."""
    if n <= 0:
        return []
    if len(codons) == 1:
        return [codons[0][0]] * n
    targets: list = []
    remainders: list = []
    allocated = 0
    for codon, frac in codons:
        exact = n * frac
        floored = int(exact)
        targets.append(floored)
        remainders.append((exact - floored, len(targets) - 1))
        allocated += floored
    shortage = n - allocated
    remainders.sort(key=lambda x: -x[0])
    for k in range(max(0, shortage)):
        # `% len` guards a malformed table whose rounding leaves a shortage
        # larger than the codon count; for a normalized table
        # shortage <= len(remainders), so this is a plain index.
        targets[remainders[k % len(remainders)][1]] += 1
    queues = [[codon] * cnt
              for (codon, _frac), cnt in zip(codons, targets) if cnt > 0]
    interleaved: list = []
    i = 0
    while any(queues):
        q = queues[i % len(queues)]
        if q:
            interleaved.append(q.pop(0))
        i += 1
    # Length guarantee (defensive — a no-op for any table whose per-AA
    # fractions sum to 1): pad short with the most-frequent synonym,
    # truncate long.
    if len(interleaved) < n:
        interleaved += [codons[0][0]] * (n - len(interleaved))
    return interleaved[:n]


def _codon_optimize(protein: str, raw: dict, *, stops: int = 1,
                    transl_table: "int | None" = None,
                    mode: str = _CODON_MODE_FREQUENCY) -> str:
    """Codon-optimize a protein against a target organism's usage table.
    Raises ValueError on unknown amino acids or an unknown ``mode``.

    ``mode`` picks the strategy (default ``"frequency"`` — byte-identical to
    every pre-2026-07-23 call):

      * ``"frequency"`` — distribute synonymous codons so each amino acid's
        codon distribution MATCHES the target's overall usage frequencies.
        Deliberately uses rare codons at their natural rate, which is why the
        result can score a LOWER CAI than the wild-type gene, and why an
        individual residue can be assigned a worse synonym than a native gene
        used at that position (the reported CGT→CGA Arg case). That is the
        strategy working as designed, not a defect — but it is a surprise if
        you read "optimize" as "maximise CAI", so:
      * ``"max_cai"`` — give EVERY position its single most-frequent synonym.
        Maximises CAI by construction and therefore can never downgrade a
        residue that already had the best codon. The trade-off is real:
        one codon per residue removes the natural rare-codon pacing and
        produces a far more repetitive sequence, so pair it with
        ``_codon_fix_sites`` / repeat-diversification if you care about
        synthesis complexity or recombination between similar constructs.

    Both modes share the stop-codon handling below, and stops don't affect CAI
    (``_codon_cai`` skips them), so the choice is confined to the coding body.

    Stop codons (2026-05-30): a run of trailing ``*`` in ``protein`` is
    honored verbatim — ``"MGK*"`` → one stop, ``"MGK**"`` → two, ``"MGK***"``
    → three — and ``stops`` is then ignored. A ``*`` anywhere but that
    trailing run raises ValueError (a premature stop is never silently
    encoded). When the protein carries no trailing ``*``, ``stops`` (default
    1, negatives clamped to 0) stop codons are appended. A SINGLE stop is
    always ``TAA`` (the strongest terminator / lowest readthrough in E. coli,
    and the historical default that downstream site-scrubbing assumes); 2+
    stops are frequency-matched to the table's OWN stop-codon usage so the
    emitted run resists readthrough with organism-appropriate diversity,
    falling back to ``TAA`` only when the table declares no stop codons.

    The amino-acid body is apportioned by `_codon_allocate`, which excludes
    stop codons from every residue's synonym pool — so a body position can
    NEVER be a stop (no premature internal stop) and is always a real codon
    (no rogue/empty base). Output length is exactly ``3*len(body) +
    3*n_stops``.

    Distinct from Angov-style codon HARMONIZATION (Angov 2011), which
    requires a SOURCE organism's codon usage to preserve relative
    rare-codon positions in the target (those positions encode
    translation pauses important for cotranslational folding). We only
    consume the target table, so this is pure optimization, not
    harmonization. Renamed 2026-05-01 to stop misleading users.

    ``transl_table`` (default: standard NCBI table 1) selects the target host's
    genetic code, so an organism with a reassigned codon (Mycoplasma table 4
    TGA=Trp, ciliate table 6 TAA/TAG=Gln, a mito code, …) optimises — and
    terminates — against the residues/stops it ACTUALLY reads instead of
    silently emitting a standard-code gene that mistranslates in the host. The
    default (``None`` / 1) path is byte-identical to before. An empty protein
    still yields the requested stop run (``stops`` default 1 ⇒ ``TAA``) so the
    ``len == 3·len(body) + 3·n_stops`` invariant holds — see
    ``test_empty_protein_is_lone_stop``."""
    if mode not in _CODON_MODES:
        raise ValueError(
            f"unknown codon mode {mode!r} — expected one of "
            f"{', '.join(repr(m) for m in _CODON_MODES)}")
    gc = (_CODON_GENETIC_CODE if transl_table in (None, 1)
          else _codon_table_for(transl_table))
    aa_codons, codon_frac = _codon_build_aa_map(raw, genetic_code=gc)
    # Peel a trailing run of stop requests; a '*' anywhere else is an error.
    body = protein
    n_trailing = 0
    while body and body[-1] == "*":
        body = body[:-1]
        n_trailing += 1
    if "*" in body:
        raise ValueError(
            "stop codon '*' is only allowed at the end of the protein")
    n_stops = n_trailing if n_trailing else max(0, int(stops))

    aa_positions: dict = {}
    for i, aa in enumerate(body):
        aa_positions.setdefault(aa.upper(), []).append(i)
    codon_at = [""] * len(body)
    for aa, positions in aa_positions.items():
        codons_for_aa = aa_codons.get(aa, [])
        if not codons_for_aa:
            raise ValueError(f"No codons for amino acid '{aa}' in this table")
        if mode == _CODON_MODE_MAX_CAI:
            # `aa_codons` is sorted by fraction descending, so [0] is the peak
            # synonym — the same codon `_codon_cai` divides by. Assigning it
            # everywhere maximises CAI and makes a downgrade impossible.
            chosen = [codons_for_aa[0][0]] * len(positions)
        else:
            chosen = _codon_allocate(codons_for_aa, len(positions))
        for pos, codon in zip(positions, chosen):
            codon_at[pos] = codon

    # Stop codons come from the TARGET genetic code, not a hardcoded standard
    # stop: under an alt code (e.g. table 6, TAA/TAG = Gln) a literal "TAA"
    # would read THROUGH. "TAA" stays the single-stop default whenever it IS a
    # stop of `gc` (every standard-ish code), so the default path is unchanged.
    code_stops = [c for c, a in gc.items() if a == "*"]
    single_stop = ("TAA" if "TAA" in code_stops
                   else (code_stops[0] if code_stops else "TAA"))
    if n_stops <= 1:
        tail = single_stop * n_stops
    else:
        stop_codons = sorted(
            ((c, codon_frac.get(c, 0.0))
             for c, (a, _ct) in raw.items() if gc.get(c) == "*"),
            key=lambda x: -x[1])
        tail = "".join(_codon_allocate(stop_codons or [(single_stop, 1.0)], n_stops))
    return "".join(codon_at) + tail


def _codon_swap_ok(seq: str, codon_start: int, alt: str, site: str, idx: int,
                   all_forbidden, before_hits, maxlen: int,
                   *, window: bool = True) -> bool:
    """Decide whether replacing ``seq[codon_start:codon_start+3]`` with the
    synonymous codon ``alt`` is an acceptable forbidden-site fix:

      (1) it must REMOVE the target forbidden occurrence ``(site, idx)``, and
      (2) it must introduce NO new forbidden hit anywhere.

    ``window=True`` (production) rescans only the slice the 3-base swap could
    affect; ``window=False`` rescans the whole candidate (the original logic),
    kept as the equivalence ORACLE that `test_codon` compares against on
    adversarial inputs. Both must return identical verdicts — this is a
    zero-tolerance subsystem (a missed forbidden site corrupts a synthetic
    gene), so the windowing is proven, not assumed.

    Why the window is exact: `_forbidden_hit_set` is a LINEAR scan, so a hit
    spanning ``[p, p+len)`` differs between before/after ONLY if it overlaps the
    changed bases ``[codon_start, codon_start+3)`` — i.e. its start ``p`` lies in
    ``[codon_start-maxlen+1, codon_start+2]`` (maxlen = longest forbidden site).
    Outside that band before≡after, so the full diff equals the windowed diff,
    and the target at ``idx`` (if removable at all) necessarily falls in it.
    """
    if window:
        a_lo = max(0, codon_start - maxlen + 1)
        a_hi = codon_start + 3                       # exclusive; affected starts
        s_hi = min(len(seq), a_hi + maxlen - 1)      # cover a maxlen site's tail
        # Candidate substring with the swap applied — no full-sequence rebuild.
        sub = seq[a_lo:codon_start] + alt + seq[codon_start + 3:s_hi]
        sub_hits = _forbidden_hit_set(sub, all_forbidden)
        after_win = {(p, a_lo + off) for (p, off) in sub_hits if a_lo + off < a_hi}
        # (1) The swap must reach idx AND the site must not survive there.
        if not (a_lo <= idx < a_hi) or (site, idx) in after_win:
            return False
        # (2) New forbidden hits can only appear inside the window, so compare
        #     `before` restricted to the same band — outside it is unchanged.
        before_win = {(p, q) for (p, q) in before_hits if a_lo <= q < a_hi}
        return not (after_win - before_win)
    # ── Reference: full-sequence rescan (pre-optimization logic / test oracle).
    cand = seq[:codon_start] + alt + seq[codon_start + 3:]
    after_hits = _forbidden_hit_set(cand, all_forbidden)
    if (site, idx) in after_hits:
        return False
    return not (after_hits - before_hits)


def _codon_fix_sites(dna: str, protein: str, raw: dict,
                     sites: "dict | None" = None,
                     *, has_appended_stop: bool = True,
                     transl_table: "int | None" = None) -> tuple:
    """Substitute synonymous codons to remove internal restriction sites.

    ``sites`` is a forward-strand ``{name: site}`` dict; reverse complements
    are added automatically for non-palindromic sites. Returns
    ``(new_dna, fixes)``.

    ``has_appended_stop`` (2026-05-27 audit-5 H2): set to True when ``dna``
    ends in a SYNTHETIC stop codon that the caller appended (e.g.
    ``_codon_optimize`` always appends TAA). The boundary check then
    skips the last codon to avoid silently substituting an appended
    stop. Set to False when the caller passes a raw CDS region without
    an appended stop — without this kwarg pre-2026-05-27 the last 1-2
    codons were silently skipped for any forbidden site overlapping
    them, and the caller's "remaining sites" check then aborted with
    "no synonymous alternative" when the code had actually refused to
    try the swap.

    Hardening (2026-04-21) — a candidate swap is accepted only if it:
      1. Actually removes the target site at the current position, AND
      2. Introduces **no new** forbidden site (forward or RC) anywhere
         in the full sequence — counted against the full input site set,
         not just the enzyme currently being iterated. This guards
         against the classic failure mode of fixing BsaI by accidentally
         spawning an Esp3I (or the RC of either) a few bases away.

    Multiple occurrences of the same site are processed left-to-right;
    each swap only needs to remove its own position, so repeated sites
    of the same enzyme are handled correctly (pre-2026-04-21 the check
    was ``site not in test`` which failed when two copies were present).
    """
    if sites is None:
        sites = _CODON_DEFAULT_FORBIDDEN
    expanded: dict = {}
    for name, site in sites.items():
        site = str(site or "").upper()
        if not site:
            continue
        # Validate as an IUPAC recognition site; skip (rather than crash the
        # whole optimize) an enzyme whose site carries a stray non-IUPAC char.
        try:
            _iupac_pattern(site)
        except ValueError:
            _log.warning("Codon site-scrub: skipping %r (%s) — not a valid "
                          "IUPAC recognition site", site, name)
            continue
        expanded[name] = site
        rc = _mut_revcomp(site)
        if rc != site:
            expanded[f"{name}_rc"] = rc
    # Flat tuple of every forbidden pattern (forward + RC). Used by the
    # per-swap cross-check to veto swaps that would introduce a NEW
    # pattern anywhere (different enzyme, different strand, different
    # position — the check is global).
    all_forbidden = tuple(expanded.values())
    # Longest forbidden site = the half-width of the window a 3-base codon swap
    # can affect (see `_codon_swap_ok`). Computed once; the catalog is fixed.
    maxlen = max((len(s) for s in all_forbidden), default=0)
    # Precompiled (cached) IUPAC matchers for the per-enzyme outer scan, so
    # a degenerate site is actually located (not literal-substring searched).
    site_pats = {nm: _iupac_pattern(st) for nm, st in expanded.items()}
    _gc = (_CODON_GENETIC_CODE if transl_table in (None, 1)
           else _codon_table_for(transl_table))
    aa_codons, _ = _codon_build_aa_map(raw, genetic_code=_gc)
    dna_list = list(dna)
    fixes: list[str] = []
    for enzyme, site in expanded.items():
        pat = site_pats[enzyme]
        pos = 0
        while True:
            seq = "".join(dna_list)
            m = pat.search(seq, pos)
            if m is None:
                break
            idx = m.start()
            fixed = False
            lo_codon = max(0, (idx // 3) - 1)
            hi_codon = (idx + len(site)) // 3 + 2
            before_hits = _forbidden_hit_set(seq, all_forbidden)
            for codon_idx in range(lo_codon, hi_codon):
                codon_start = codon_idx * 3
                # 2026-05-27 (audit-5 H2): the `- 3` skip is only
                # valid when the caller appended a synthetic stop
                # codon to `dna_list`. Raw CDS regions (passed by
                # `_design_gb_primers`) don't end in an appended stop,
                # and skipping their last codon silently leaves a
                # forbidden site overlapping the C-terminus unfixed.
                last_safe = (
                    len(dna_list) - 3 if has_appended_stop
                    else len(dna_list)
                )
                if codon_start + 3 > last_safe:
                    break
                if codon_idx >= len(protein):
                    break
                aa = protein[codon_idx].upper()
                current = "".join(dna_list[codon_start:codon_start + 3])
                for alt, frac in aa_codons.get(aa, []):
                    if alt == current:
                        continue
                    # Windowed forbidden-site check (see `_codon_swap_ok`):
                    # rescans only the slice the 3-base swap can affect, not the
                    # whole candidate — provably identical to the full rescan,
                    # and the equivalence is pinned by test_codon's oracle test.
                    if not _codon_swap_ok(seq, codon_start, alt, site, idx,
                                          all_forbidden, before_hits, maxlen):
                        continue
                    dna_list[codon_start:codon_start + 3] = list(alt)
                    strand = " (rc)" if enzyme.endswith("_rc") else ""
                    fixes.append(
                        f"{enzyme.replace('_rc', '')}{strand} at nt {idx+1}: "
                        f"{current}→{alt} (codon {codon_idx+1} {aa}, "
                        f"freq={frac:.3f})"
                    )
                    fixed = True
                    break
                if fixed:
                    break
            if not fixed:
                pos = idx + 1
    return "".join(dna_list), fixes


_CODON_GC_WINDOW_DEFAULT = 50


def _codon_gc_window_range(dna: str,
                           window: int = _CODON_GC_WINDOW_DEFAULT) -> tuple:
    """``(min_pct, max_pct)`` GC over every ``window``-base window.

    This is what a caller re-measures AFTER `_codon_fix_gc_window` to learn
    whether the requested band was actually reachable: a Lys/Ile/Phe/Asn
    stretch has no synonym richer than one G/C per codon, so it tops out near
    33% and a 35% floor cannot be met there by any number of swaps. Reporting
    the achieved range beats implying success.
    """
    n = len(dna)
    if n == 0:
        return (0.0, 0.0)
    w = min(window, n) if window > 0 else n
    run = sum(1 for b in dna[:w] if b in "GCgc")
    lo = hi = run / w * 100.0
    for i in range(1, n - w + 1):
        if dna[i - 1] in "GCgc":
            run -= 1
        if dna[i + w - 1] in "GCgc":
            run += 1
        pct = run / w * 100.0
        lo = min(lo, pct)
        hi = max(hi, pct)
    return (lo, hi)


def _codon_fix_gc_window(dna: str, protein: str, raw: dict, *,
                         window: int = _CODON_GC_WINDOW_DEFAULT,
                         min_gc: "float | None" = None,
                         max_gc: "float | None" = None,
                         sites: "dict | None" = None,
                         has_appended_stop: bool = True,
                         transl_table: "int | None" = None,
                         avoid_kmers: "set | None" = None,
                         kmer_len: int = 0) -> tuple:
    """Pull every ``window``-base stretch of ``dna`` into the
    ``[min_gc, max_gc]`` band by swapping synonymous codons. Both bounds are
    PERCENT (matching `_codon_gc`); pass either alone for a one-sided floor or
    ceiling. Returns ``(new_dna, fixes)``; a no-op returning ``(dna, [])`` when
    neither bound is given.

    Whole-sequence GC can sit at a comfortable 50% while some 50-base stretch
    is 15% or 85% — and it is the LOCAL window that stalls synthesis and drives
    secondary structure, which a global GC number cannot see.

    A swap is applied only if it keeps the residue AND introduces no new
    ``sites`` hit anywhere, so running this after `_codon_fix_sites` can't
    re-create a forbidden site the scrub just removed. Pass ``avoid_kmers``
    (with ``kmer_len``) — the shared-run set `_codon_diversify` worked
    against — to extend the same guard to repeats: without it, running the GC
    pass after diversification silently re-creates shared runs the
    diversification just broke, because a GC swap has no reason to prefer a
    codon the neighbouring construct doesn't use. Windows no synonymous
    swap can rescue (a stretch of Met/Trp has no choices at all) are left alone
    and reported via the returned fixes being silent about them — the caller
    re-measures rather than trusting the pass to have succeeded.
    """
    if min_gc is None and max_gc is None:
        return dna, []
    if window <= 0:
        raise ValueError("'window' must be a positive number of bases")
    if (min_gc is not None and max_gc is not None and min_gc > max_gc):
        raise ValueError(
            f"min_gc ({min_gc}) must not exceed max_gc ({max_gc})")

    _gc_code = (_CODON_GENETIC_CODE if transl_table in (None, 1)
                else _codon_table_for(transl_table))
    aa_codons, _frac = _codon_build_aa_map(raw, genetic_code=_gc_code)

    # Same forward+RC expansion as `_codon_fix_sites`, so "don't create a new
    # site" means the same thing in both passes.
    all_forbidden: list[str] = []
    for _name, site in (sites or {}).items():
        s = str(site or "").upper()
        if not s:
            continue
        try:
            _iupac_pattern(s)
        except ValueError:
            continue
        all_forbidden.append(s)
        rc = _mut_revcomp(s)
        if rc != s:
            all_forbidden.append(rc)
    forbidden = tuple(all_forbidden)

    dna_list = list(dna)
    last_safe = len(dna_list) - 3 if has_appended_stop else len(dna_list)

    def _worst_window(seq: str, skip: set):
        """Largest band violation as ``(lo, hi, need_more_gc)``, else None.
        Rolling count keeps this O(len(seq)) so the repair loop stays linear
        per swap instead of rescanning every window from scratch. Windows in
        ``skip`` are passed over, so one unrescuable stretch doesn't stop the
        pass from fixing every other window."""
        n = len(seq)
        if n == 0:
            return None
        w = min(window, n)
        run = sum(1 for b in seq[:w] if b in "GCgc")
        worst = None            # (violation, key)
        for lo in range(0, n - w + 1):
            if lo:
                if seq[lo - 1] in "GCgc":
                    run -= 1
                if seq[lo + w - 1] in "GCgc":
                    run += 1
            pct = run / w * 100.0
            if min_gc is not None and pct < min_gc:
                key, gap = (lo, lo + w, True), min_gc - pct
            elif max_gc is not None and pct > max_gc:
                key, gap = (lo, lo + w, False), pct - max_gc
            else:
                continue
            if key in skip:
                continue
            if worst is None or gap > worst[0]:
                worst = (gap, key)
        return worst[1] if worst else None

    fixes: list[str] = []
    stuck: set = set()
    # One swap can only ever help one codon, so a codon-count budget bounds the
    # loop even if a pathological table keeps offering zero-gain swaps.
    for _ in range(max(1, len(dna_list) // 3)):
        seq = "".join(dna_list)
        target = _worst_window(seq, stuck)
        if target is None:
            break
        lo, hi, need_more = target
        before_hits = _forbidden_hit_set(seq, forbidden) if forbidden else set()
        # Candidate codons: every codon overlapping the offending window.
        cands = []
        for codon_idx in range(lo // 3, min(hi // 3 + 1, len(protein))):
            codon_start = codon_idx * 3
            if codon_start + 3 > last_safe:
                continue
            aa = protein[codon_idx].upper()
            current = "".join(dna_list[codon_start:codon_start + 3])
            cur_gc = sum(1 for b in current if b in "GCgc")
            for alt, afrac in aa_codons.get(aa, []):
                if alt == current:
                    continue
                delta = sum(1 for b in alt if b in "GCgc") - cur_gc
                if (need_more and delta > 0) or (not need_more and delta < 0):
                    # Best GC movement first; among equals prefer the commoner
                    # codon so the repair costs as little CAI as possible.
                    cands.append((-abs(delta), -afrac, codon_idx,
                                  codon_start, aa, current, alt))
        cands.sort()
        applied = False
        # Bind a definitely-set local so the membership test below doesn't
        # trip the "in <set | None>" type check inside the generator.
        kmer_guard: set = avoid_kmers or set()
        check_kmers = bool(kmer_guard) and kmer_len > 0
        for _d, _f, codon_idx, codon_start, aa, current, alt in cands:
            dna_list[codon_start:codon_start + 3] = list(alt)
            # `after` is O(n) to build — only materialise it when a guard
            # actually reads it. In the common GC-only case (no forbidden
            # sites, no repeat k-mers) neither branch runs, so we skip the
            # rebuild entirely; building it unconditionally per candidate was
            # the O(n²) that made a 6 kb CDS take seconds.
            if forbidden and not (_forbidden_hit_set("".join(dna_list),
                                                     forbidden) <= before_hits):
                dna_list[codon_start:codon_start + 3] = list(current)
                continue
            if check_kmers:
                # Only windows overlapping the 3 bases we changed can be new,
                # so test a local slice, not the whole sequence.
                lo_k = max(0, codon_start - kmer_len + 1)
                hi_k = min(len(dna_list) - kmer_len, codon_start + 2)
                if any("".join(dna_list[i:i + kmer_len]) in kmer_guard
                       for i in range(lo_k, hi_k + 1)):
                    dna_list[codon_start:codon_start + 3] = list(current)
                    continue
            direction = "raise" if need_more else "lower"
            fixes.append(
                f"GC {direction} at nt {lo + 1}-{hi}: {current}→{alt} "
                f"(codon {codon_idx + 1} {aa})"
            )
            applied = True
            break
        if not applied:
            # No synonymous move can rescue this window; remember it so the
            # next iteration attacks a different one instead of spinning.
            stuck.add((lo, hi, need_more))
    return "".join(dna_list), fixes


_CODON_REPEAT_RUN_DEFAULT = 25

# The GC-window and repeat-diversification passes are O(n²) — each makes up to
# one synonymous swap per codon, and every swap rescans the whole sequence to
# find the next-worst window / next shared run. That's fine for real proteins
# (a 1500-codon CDS finishes in ~1.5 s) but a titin-scale input would peg a
# worker thread for a minute. Callers that request those passes on a longer CDS
# are refused with a clear limit rather than silently hung; the plain optimize
# and the (linear) motif scrub stay uncapped. Chosen well above any realistic
# single-ORF codon-optimization target.
_CODON_SCRUB_MAX_CODONS = 1500


def _codon_kmer_set(references, k: int) -> set:
    """Every ``k``-mer across ``references`` — the shared-run lookup that both
    the repeat breaker and the GC pass's repeat guard test against."""
    out: set = set()
    if k <= 0:
        return out
    for ref in references or []:
        text = str(ref or "").upper()
        for i in range(len(text) - k + 1):
            out.add(text[i:i + k])
    return out


def _codon_shared_runs(dna: str, references, min_run: int) -> "list[int]":
    """Start offsets in ``dna`` of every ``min_run``-base window that occurs
    verbatim in at least one reference. Overlapping hits are all reported; a
    caller breaking a long run re-scans after each swap."""
    if min_run <= 0 or not dna:
        return []
    kmers = _codon_kmer_set(references, min_run)
    if not kmers:
        return []
    up = dna.upper()
    return [i for i in range(len(up) - min_run + 1)
            if up[i:i + min_run] in kmers]


def _codon_diversify(dna: str, protein: str, raw: dict, references, *,
                     min_run: int = _CODON_REPEAT_RUN_DEFAULT,
                     sites: "dict | None" = None,
                     has_appended_stop: bool = True,
                     transl_table: "int | None" = None) -> tuple:
    """Break every perfect-identity run of ``min_run``+ bases that ``dna``
    shares with any sequence in ``references``, using synonymous codon swaps.
    Returns ``(new_dna, fixes)``.

    Two constructs carrying a long identical stretch can recombine when they
    meet in the same cell, so a panel of same-protein variants optimized one
    at a time is a real hazard: each is optimized against the same table, so
    they converge on the same codons and end up sharing hundreds of bases (the
    reported case was a 363 bp perfect repeat across a co-transformed panel).
    Optimizing each variant against the ones already built removes that.

    Direct repeats only — the shared-sequence recombination substrate. A swap
    is applied only if it keeps the residue, actually breaks the window it was
    chosen for, and introduces no new ``sites`` hit, so this composes with the
    motif scrub in either order. Runs no synonymous swap can break (a
    Met/Trp-only stretch is invariant) are left alone; re-scan with
    `_codon_shared_runs` to see what survived rather than assuming none did.
    """
    refs = [str(r or "").upper() for r in (references or [])]
    refs = [r for r in refs if len(r) >= min_run]
    if not refs or min_run <= 0:
        return dna, []

    _gc_code = (_CODON_GENETIC_CODE if transl_table in (None, 1)
                else _codon_table_for(transl_table))
    aa_codons, _frac = _codon_build_aa_map(raw, genetic_code=_gc_code)

    all_forbidden: list[str] = []
    for _name, site in (sites or {}).items():
        s = str(site or "").upper()
        if not s:
            continue
        try:
            _iupac_pattern(s)
        except ValueError:
            continue
        all_forbidden.append(s)
        rc = _mut_revcomp(s)
        if rc != s:
            all_forbidden.append(rc)
    forbidden = tuple(all_forbidden)

    kmers = _codon_kmer_set(refs, min_run)

    dna_list = list(dna)
    last_safe = len(dna_list) - 3 if has_appended_stop else len(dna_list)
    fixes: list[str] = []
    stuck: set = set()

    for _ in range(max(1, len(dna_list) // 3)):
        seq = "".join(dna_list)
        hit = next((i for i in range(len(seq) - min_run + 1)
                    if i not in stuck and seq[i:i + min_run] in kmers), None)
        if hit is None:
            break
        before_hits = (_forbidden_hit_set(seq, forbidden) if forbidden
                       else set())
        # Break from the MIDDLE outward: a change near the centre splits the
        # run into two pieces that are each shorter than min_run, where an
        # edge change just slides the window along by one base.
        mid = hit + min_run // 2
        order = sorted(range(hit // 3, min((hit + min_run) // 3 + 1,
                                           len(protein))),
                       key=lambda ci: abs(ci * 3 + 1 - mid))
        applied = False
        for codon_idx in order:
            codon_start = codon_idx * 3
            if codon_start + 3 > last_safe:
                continue
            aa = protein[codon_idx].upper()
            current = "".join(dna_list[codon_start:codon_start + 3])
            # Commonest synonym first so diversification costs least CAI.
            for alt, afrac in aa_codons.get(aa, []):
                if alt == current:
                    continue
                dna_list[codon_start:codon_start + 3] = list(alt)
                # Did the swap actually break THIS window? Check the local
                # slice, not a freshly-joined whole sequence (that rebuild per
                # candidate was O(n²) on a long CDS).
                if "".join(dna_list[hit:hit + min_run]) in kmers:
                    dna_list[codon_start:codon_start + 3] = list(current)
                    continue        # didn't actually break this window
                if forbidden and not (_forbidden_hit_set("".join(dna_list),
                                                         forbidden)
                                      <= before_hits):
                    dna_list[codon_start:codon_start + 3] = list(current)
                    continue
                fixes.append(
                    f"repeat at nt {hit + 1}-{hit + min_run}: "
                    f"{current}→{alt} (codon {codon_idx + 1} {aa}, "
                    f"freq={afrac:.3f})"
                )
                applied = True
                break
            if applied:
                break
        if not applied:
            stuck.add(hit)
    return "".join(dna_list), fixes


def _codon_forbidden_sites() -> "dict[str, str]":
    """Resolve the persisted ``codon_forbidden_enzymes`` name list into the
    ``{name: recognition_site}`` map that `_codon_fix_sites` consumes, drawing
    sites from the merged built-in + custom enzyme set. A name that's unknown
    (e.g. a custom enzyme since deleted) or whose site isn't valid IUPAC is
    skipped. An empty / missing list means 'scrub nothing'. Default ['BsaI'].
    """
    names = _get_setting("codon_forbidden_enzymes", ["BsaI"])
    if not isinstance(names, list):
        names = ["BsaI"]
    enz = _state._all_enzymes_hook()
    out: dict[str, str] = {}
    for n in names:
        info = enz.get(str(n))
        if not info:
            continue
        site = str(info[0] or "").upper()
        if not site:
            continue
        try:
            _iupac_pattern(site)
        except ValueError:
            continue
        out[str(n)] = site
    return out


def _codon_cai(dna: str, raw: dict, *, transl_table: "int | None" = None) -> float:
    """Codon Adaptation Index (geometric mean of per-codon freq ÷ peak freq
    of its amino-acid synonymy group). Skips stops and unknown codons.

    ``transl_table`` labels codons by the target host's genetic code (default
    standard), so a reassigned codon is scored against the RIGHT synonymy group
    — pre-fix CAI keyed on the table's stored label, mis-scoring alt-code
    tables. Display-only; never affects an emitted sequence."""
    import math
    gc = (_CODON_GENETIC_CODE if transl_table in (None, 1)
          else _codon_table_for(transl_table))
    aa_codons, codon_frac = _codon_build_aa_map(raw, genetic_code=gc)
    w: list[float] = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        entry = raw.get(codon)
        aa = gc.get(codon)
        if not entry or aa is None or aa == "*":
            continue
        peak = aa_codons[aa][0][1] if aa in aa_codons else 0.0
        if peak > 0:
            w.append(codon_frac.get(codon, 0.0) / peak)
    if not w:
        return 0.0
    return math.exp(sum(math.log(max(v, 1e-10)) for v in w) / len(w))


def _codon_gc(dna: str) -> float:
    """GC%. Empty string → 0."""
    if not dna:
        return 0.0
    gc = sum(1 for c in dna.upper() if c in "GC")
    return gc / len(dna) * 100.0


# ── codon usage chart (Phase D, moved from hub) ─────────────────────────────
# Renders a per-amino-acid synonymous-codon usage chart as a Rich-markup
# STRING (no Rich object dependency). _AA_NAME_3 (AA 1->3 letter) is used only
# here.
# Single-letter → three-letter amino-acid names for the genetic-code chart
# (`_render_codon_chart`). Stop is spelled "Stop" to match the canonical
# textbook wall chart. The full-name catalog lives on `AminoAcidPickerModal`.
_AA_NAME_3: dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Stop",
}


def _render_codon_chart(raw: dict, *, rna: bool = True) -> str:
    """Render a codon-usage table as the classic textbook genetic-code grid,
    returned as a Rich-markup string.

    The layout matches the canonical wall chart: the four 1st-base blocks
    (U/C/A/G) stack vertically, the 2nd base runs across the four columns,
    and the four lines inside every cell are the 3rd base (U/C/A/G). Each
    codon is annotated with its usage *within its amino-acid family* (the
    relative synonymous usage, as a percentage), and each amino-acid family's
    single most-used codon is highlighted bold green for easy visual
    identification. The choice is FAMILY-WIDE (matching the % shown), so a
    family split across two cells — Leu, Ser, Arg, and the three stops
    (UAA/UAG/UGA, treated as one family) — yields exactly one highlight, not
    one per cell. A codon with no usage in the table is never highlighted.
    Synonymous residues are bracketed and named (3-letter; stops as "Stop").
    Codons missing from the table render a dim "·" placeholder.

    Bases display as RNA (U) by default to match the iconic chart; lookups
    are always by the stored DNA (T) key. Pure — display only, no I/O."""
    _aa_codons, codon_frac = _codon_build_aa_map(raw)

    def _count(codon: str) -> int:
        v = raw.get(codon)
        return int(v[1]) if v else 0

    # The single most-used codon in each amino-acid family is highlighted green.
    # The choice is FAMILY-WIDE (matching the relative-synonymous-usage % shown),
    # so a family split across two cells — Leu (UU+CU), Ser (UC+AG), Arg (CG+AG),
    # and the stops UAA/UAG/UGA (UA+UG) — yields exactly ONE highlight, never one
    # per cell. Stops are one family. Deterministic tie-break by codon; a family
    # with no usage in the table gets no highlight.
    _fam: dict = {}
    for _codon, _aa in _CODON_GENETIC_CODE.items():
        _fam.setdefault(_aa, []).append(_codon)
    dominant: set = set()
    for _members in _fam.values():
        champ = sorted(_members, key=lambda c: (-_count(c), c))[0]
        if _count(champ) > 0:
            dominant.add(champ)

    order = "TCAG"                       # U, C, A, G in DNA (T) coordinates
    CW, AAW, LGUT, RGUT = 17, 4, 2, 2    # cell / AA-label / gutter widths
    label_row = {1: 0, 2: 0, 3: 1, 4: 1}  # which row of a run carries the name

    def _disp(b: str) -> str:
        return "U" if (rna and b == "T") else b

    def _cell(first: str, second: str) -> list:
        """Four markup lines (3rd base = U/C/A/G) for one grid cell, each
        exactly CW visible columns wide (markup tags are zero-width)."""
        codons = [first + second + third for third in order]
        aas = [_CODON_GENETIC_CODE[c] for c in codons]
        # Group the four residues into runs of equal AA (for the brackets). The
        # green highlight is the family-wide champion computed above (`dominant`).
        run_of: dict = {}
        i = 0
        while i < 4:
            j = i
            while j < 4 and aas[j] == aas[i]:
                j += 1
            for r in range(i, j):
                run_of[r] = (i, j - i)
            i = j
        out: list = []
        for r in range(4):
            cdn = codons[r]
            show = "".join(_disp(b) for b in cdn)
            frac = codon_frac.get(cdn)
            pct = ("[dim]   ·[/dim]" if frac is None
                   else f"{round(frac * 100):>3d}%")     # 4 visible cols
            field = f"{show} {pct}"                       # 8 visible cols
            if cdn in dominant:
                field = f"[b green]{field}[/]"            # stark + easy to spot
            start, length = run_of[r]
            pos = r - start
            if length == 1:
                g = "─"
            elif pos == 0:
                g = "╮"
            elif pos == length - 1:
                g = "╯"
            elif pos == label_row[length]:
                g = "┤"
            else:
                g = "│"
            if pos == label_row[length]:
                nm = _AA_NAME_3.get(aas[r], aas[r])
                lab = (f"[red]{nm:<{AAW}}[/red]" if aas[r] == "*"
                       else f"{nm:<{AAW}}")
            else:
                lab = " " * AAW
            out.append(f" {field} [dim]{g}[/dim] {lab} ")   # CW visible
        return out

    def _rule(left: str, mid: str, right: str) -> str:
        return (" " * LGUT) + left + mid.join(["─" * CW] * 4) + right + \
               (" " * RGUT)

    width = LGUT + 1 + CW * 4 + 3 + 1 + RGUT
    lines: list = [f"{'second base':^{width}}"]
    # Column header: 2nd-base letter centred over each column.
    lines.append((" " * (LGUT + 1))
                 + " ".join(f"[b]{_disp(s):^{CW}}[/b]" for s in order)
                 + " " + (" " * RGUT))
    lines.append(_rule("┌", "┬", "┐"))
    for bi, first in enumerate(order):
        cells = [_cell(first, second) for second in order]
        for r in range(4):
            row = "│".join(cells[c][r] for c in range(4))
            left = f" [b]{_disp(first)}[/b]" if r == 1 else "  "
            right = f" [b]{_disp(order[r])}[/b]"
            lines.append(f"{left}│{row}│{right}")
        lines.append(_rule("├", "┼", "┤") if bi < 3 else _rule("└", "┴", "┘"))
    return "\n".join(lines)



# ── codon-table NETWORK builders (Phase D, moved from hub) ──────────────────
# Build a codon-usage table from a live source: Kazusa (HTML scrape) or an
# NCBI genome (datasets CDS zip -> highly-expressed-gene table). Egress is
# gated via _state._demo_block_network_hook + the splicecraft_net hardened
# opener; lazy urllib/json/zipfile/io/socket/itertools inside the fns.
def _build_heg_table_from_cds(cds_fasta_text: str,
                              mode: str = "heg") -> "tuple[dict, dict]":
    """Build a codon-usage table from a genome's CDS FASTA (the
    ``cds_from_genomic.fna`` shape NCBI Datasets ships). Pure — no I/O, no
    network — so it unit-tests against an inline FASTA string.

    ``mode``:
      * ``"heg"``    — count only highly-expressed genes (ribosomal proteins),
        the right signal for heterologous-expression optimization. Amino acids
        absent from the r-protein set (ribosomal proteins can be Cys/Trp-poor)
        are BACKFILLED from the whole-genome counts so every protein still
        optimizes — backfill only supplies a synonym pool for an otherwise-
        missing residue; it never overrides the HEG bias where it exists.
      * ``"genome"`` — count every CDS (whole-genome average; what Kazusa
        approximates).

    Returns ``(raw, stats)`` where ``raw`` is the ``{codon: (aa, count)}``
    registry shape consumed by `_codon_optimize` / `_codon_tables_add`, and
    ``stats`` carries mode / n_cds_total / n_cds_heg / n_codons / aa_coverage /
    backfilled / gc3 for the status line. Each record is read in frame 0 (CDS
    FASTA is in-frame); codons containing non-ACGT (N, gaps) are skipped — a
    codon is counted iff it is one of the 64 keys of `_CODON_GENETIC_CODE`.

    Raises ``ValueError`` on an unknown ``mode``, when no usable CDS were found,
    or (heg mode) when the genome carries no ribosomal-protein CDS — callers
    surface it in a status line.
    """
    if mode not in ("heg", "genome"):
        raise ValueError(f"unknown mode {mode!r} (expected 'heg' or 'genome')")
    import io
    from collections import Counter

    def _is_rprotein(header: str) -> bool:
        h = header.lower()
        if "transferase" in h:           # rimI / prmA modification enzymes
            return False                 # carry "ribosomal protein" but aren't one
        if "ribosomal protein" in h:
            return True
        # NCBI cds_from_genomic deflines carry [gene=rplB] / [gene=rpsL] /
        # [gene=rpmA] for the 50S/30S/large-subunit r-proteins.
        return bool(re.search(r"\[gene=rp[lsm]", h))

    genome_counts: "Counter" = Counter()
    heg_counts: "Counter" = Counter()
    n_cds_total = 0
    n_cds_heg = 0

    def _consume(hdr: "str | None", parts: list) -> None:
        """Count one CDS record's in-frame ACGT codons into the two Counters."""
        nonlocal n_cds_total, n_cds_heg
        if hdr is None:
            return
        s = "".join(parts).upper().replace("U", "T")
        usable = (len(s) // 3) * 3
        if usable <= 0:
            return
        n_cds_total += 1
        is_heg = _is_rprotein(hdr)
        if is_heg:
            n_cds_heg += 1
        for i in range(0, usable, 3):
            codon = s[i:i + 3]
            if codon in _CODON_GENETIC_CODE:   # 64 ACGT keys → skips N/gaps
                genome_counts[codon] += 1
                if is_heg:
                    heg_counts[codon] += 1

    # Stream line-by-line (lazy — avoids materialising a whole-file line list
    # AND a records list; at the 256 MB download cap that's the difference
    # between ~1× and ~2-3× the input in peak memory). Each record is consumed
    # the instant its sequence is complete; only one record is held at a time.
    header: "str | None" = None
    seq_parts: list = []
    for raw_line in io.StringIO(cds_fasta_text):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            _consume(header, seq_parts)
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    _consume(header, seq_parts)

    if not genome_counts:
        raise ValueError("no usable CDS found — expected an in-frame CDS "
                         "FASTA (e.g. NCBI cds_from_genomic.fna)")
    base = heg_counts if mode == "heg" else genome_counts
    if not base:
        raise ValueError("no highly-expressed (ribosomal-protein) CDS found "
                         "in this genome — try whole-genome mode")

    raw: dict = {c: (_CODON_GENETIC_CODE[c], int(n)) for c, n in base.items()}
    backfilled: list = []
    if mode == "heg":
        have_aa = {_CODON_GENETIC_CODE[c] for c in raw}
        for c, n in genome_counts.items():
            aa = _CODON_GENETIC_CODE[c]
            if aa not in have_aa:
                raw[c] = (aa, int(n))
                if aa not in backfilled:
                    backfilled.append(aa)

    n_codons = sum(n for _aa, n in raw.values())
    gc3 = ((sum(n for c, (_aa, n) in raw.items() if c[2] in "GC")
            / n_codons * 100) if n_codons else 0.0)
    stats = {
        "mode":        mode,
        "n_cds_total": n_cds_total,
        "n_cds_heg":   n_cds_heg,
        "n_codons":    n_codons,
        "aa_coverage": len({_aa for _aa, _n in raw.values()} - {"*"}),
        "backfilled":  sorted(backfilled),
        "gc3":         round(gc3, 1),
    }
    return raw, stats


def _codon_parse_kazusa_html(html: str) -> "dict | None":
    """Parse Kazusa showcodon.cgi GCG-format HTML. Returns {codon: (aa, count)}
    or None on failure."""
    pre = re.search(r"<[Pp][Rr][Ee]>(.*?)</[Pp][Rr][Ee]>", html, re.DOTALL)
    text = pre.group(1) if pre else html
    pat = re.compile(r"\b([ACGTU]{3})\b\s+(\d+(?:\.\d+)?)", re.IGNORECASE)
    raw: dict = {}
    for m in pat.finditer(text):
        rna = m.group(1).upper()
        dna = rna.replace("U", "T")
        if dna not in _CODON_GENETIC_CODE or dna in raw:
            continue
        try:
            count = round(float(m.group(2)))
        except ValueError:
            continue
        raw[dna] = (_CODON_GENETIC_CODE[dna], count)
    # 2026-05-27 (audit-5 domesticator M1): require all 64 codons.
    # Pre-fix the threshold was 60, silently accepting up to 4
    # missing codons; if a rare codon (e.g. ATA-Ile) was missing
    # for an organism, `_codon_optimize` raised `ValueError("No
    # codons for amino acid …")` only when the protein happened to
    # contain that AA — surface unpredictably at design time. Now
    # malformed tables are rejected upfront.
    if len([c for c in raw if raw[c][0] != "?"]) < 64:
        return None
    return raw


_KAZUSA_MAX_RESPONSE_BYTES = 1 * 1024 * 1024


# NCBI Datasets CDS-FASTA zip download (codon-table genome builder). A
# bacterial genome's CDS set zips to ~1 MB; 256 MB caps a hostile / oversized
# response. The builder targets prokaryotic hosts — an accidentally-eukaryotic
# CDS set trips this cap and surfaces a "use a smaller / prokaryotic genome"
# error rather than ballooning memory.
_NCBI_CDS_ZIP_MAX_BYTES    = 256 * 1024 * 1024
# Local CDS FASTA cap (the offline `_file_build_codon_table` source). Matches the
# NCBI zip cap — it's the same CDS data, just read from disk instead of fetched.
_CDS_FILE_MAX_BYTES        = 256 * 1024 * 1024


# A real Datasets CDS zip carries a handful of members; cap the walk so a
# hostile / MITM'd response can't make us iterate millions of entries.
_NCBI_CDS_ZIP_MAX_MEMBERS  = 10_000


_NCBI_DATASETS_TIMEOUT_S   = 60.0


def _codon_fetch_kazusa(taxid: str, timeout: float = 15.0) -> tuple:
    """Fetch codon usage from Kazusa for an NCBI taxid. Returns
    (raw_dict_or_None, status_message). Pure network call — callers should
    invoke from a worker thread."""
    _state._demo_block_network_hook("Kazusa codon fetch")
    import urllib.request
    import urllib.error as _urllib_error  # sweep #25 — narrow excepts
    taxid = str(taxid).strip()
    if not taxid.isdigit():
        return None, f"Invalid taxid '{taxid}' (must be numeric)"
    url = (f"https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi"
           f"?species={taxid}&aa=1&style=GCG")
    # Route through the shared hardened opener (verifying SSL, bounded
    # redirects, https-downgrade refusal, private-IP/SSRF filter) like every
    # other external fetch — the egress guard above already fired.
    opener = _build_hardened_url_opener()
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with opener.open(req, timeout=timeout) as r:
            raw = r.read(_KAZUSA_MAX_RESPONSE_BYTES + 1)
        if len(raw) > _KAZUSA_MAX_RESPONSE_BYTES:
            _log.warning("Kazusa response exceeded %d bytes; aborting",
                          _KAZUSA_MAX_RESPONSE_BYTES)
            return None, "Kazusa returned an oversized response"
        html = raw.decode("utf-8", errors="replace")
    except (OSError, _urllib_error.URLError) as exc:
        # Sweep #25 (2026-05-23): narrowed from bare `Exception`.
        _log.exception("Kazusa fetch failed for taxid %s", taxid)
        return None, f"Network error: {exc}"
    low = html.lower()
    if "not found" in low or "no data" in low:
        return None, f"Taxid {taxid} not found in Kazusa database"
    raw = _codon_parse_kazusa_html(html)
    if raw is None:
        return None, f"Could not parse Kazusa table for taxid {taxid}"
    return raw, f"Fetched from Kazusa: {len(raw)} codons (taxid {taxid})"


# ── Genome → codon table (NCBI Datasets HEG builder) ──────────────────────────
# Build a codon-usage table from a genome's CDS instead of Kazusa: resolve a
# taxid → RefSeq reference accession, download the CDS_FASTA zip from the NCBI
# Datasets v2 API, extract the in-frame `cds_from_genomic.fna`, hand it to
# `_build_heg_table_from_cds`. Network reads are size-capped + retried once
# (250 ms) through the shared hardened opener ([PIT-20] / [RECIPE] "New HTTP
# fetch"); the zip is opened with the existing zip-safety guards. Worker-thread
# helpers — never call from the UI thread.
_NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"


def _genome_datasets_request(url: str, max_bytes: int, timeout: float) -> bytes:
    """GET `url` through the shared hardened opener, one 250 ms-backoff retry on
    transient errors, body bounded at `max_bytes` (ValueError if exceeded —
    [PIT-20]). HTTP/URL errors propagate for the caller to message."""
    _state._demo_block_network_hook("Genome download")
    import socket
    import time as _time
    import urllib.error
    import urllib.request
    opener = _build_hardened_url_opener()
    req = urllib.request.Request(
        url, headers={"User-Agent": f"SpliceCraft/{_state._sc_version}"})
    last_exc: "BaseException | None" = None
    for attempt in range(2):                       # 1 try + 1 retry
        try:
            with opener.open(req, timeout=timeout) as resp:
                body = resp.read(max_bytes + 1)
            if len(body) > max_bytes:
                raise ValueError(
                    f"response exceeded {max_bytes:,}-byte cap "
                    f"(genome too large — try a prokaryote / smaller assembly)")
            return body
        except urllib.error.HTTPError:
            raise                                  # 4xx/5xx is permanent
        except (urllib.error.URLError, socket.timeout) as exc:
            last_exc = exc
            if attempt == 0:
                _time.sleep(_HMM_DB_RETRY_BACKOFF_S)
                continue
            raise
    assert last_exc is not None                    # unreachable
    raise last_exc


def _genome_resolve_reference_accession(
        taxid: str, timeout: float = _NCBI_DATASETS_TIMEOUT_S) -> tuple:
    """Resolve an NCBI taxid → its RefSeq reference-genome assembly accession.
    Returns (accession, organism_name) on success or (None, error_message)."""
    import json as _json
    taxid = str(taxid).strip()
    if not taxid.isdigit():
        return None, f"Invalid taxid '{taxid}' (must be numeric)"
    url = (f"{_NCBI_DATASETS_BASE}/genome/taxon/{taxid}/dataset_report"
           f"?filters.reference_only=true&filters.assembly_source=refseq"
           f"&page_size=1")
    try:
        body = _genome_datasets_request(url, _NCBI_MAX_RESPONSE_BYTES, timeout)
        report = _json.loads(body.decode("utf-8", errors="replace"))
    except ValueError as exc:                      # size cap / malformed JSON
        return None, f"Could not read assembly report for taxid {taxid}: {exc}"
    except Exception as exc:                        # network / HTTP
        _log.exception("Datasets taxon resolve failed for %s", taxid)
        return None, f"Could not resolve taxid {taxid}: {exc}"
    reports = report.get("reports") if isinstance(report, dict) else None
    if not reports:
        return None, (f"No RefSeq reference genome for taxid {taxid} — supply a "
                      f"specific assembly accession (GCF_…) instead")
    rec = reports[0] or {}
    acc = str(rec.get("accession", "") or "")
    org = rec.get("organism") if isinstance(rec.get("organism"), dict) else {}
    name = str((org or {}).get("organism_name", "") or "")
    if not acc:
        return None, f"Reference assembly for taxid {taxid} carried no accession"
    return acc, name


def _genome_extract_cds_fasta(zip_bytes: bytes) -> str:
    """Extract the `cds_from_genomic.fna` member from an in-memory NCBI Datasets
    CDS zip, reusing the zip-safety guards (`_is_safe_zip_member_name`, member-
    size cap, bounded read — [SUB-plasmidsaurus]). ValueError when the archive
    is not a zip, carries no CDS FASTA (unannotated assembly), or trips a cap."""
    import io
    import zipfile
    try:
        with zipfile.ZipFile(io.BytesIO(zip_bytes), "r") as zf:
            # Cap the member walk and require an EXACT basename so a decoy
            # like `evilcds_from_genomic.fna` can't be selected by `endswith`.
            import itertools
            member = next(
                (info for info in itertools.islice(
                     zf.infolist(), _NCBI_CDS_ZIP_MAX_MEMBERS)
                 if Path(info.filename).name == "cds_from_genomic.fna"
                 and _is_safe_zip_member_name(info.filename)), None)
            if member is None:
                raise ValueError(
                    "assembly has no annotated CDS (no cds_from_genomic.fna) — "
                    "try a RefSeq (GCF_) accession that includes annotation")
            if member.file_size > _NCBI_CDS_ZIP_MAX_BYTES:
                raise ValueError(
                    f"CDS member too large ({member.file_size:,} bytes; cap "
                    f"{_NCBI_CDS_ZIP_MAX_BYTES:,})")
            with zf.open(member, "r") as fh:
                raw = fh.read(_NCBI_CDS_ZIP_MAX_BYTES + 1)
    except zipfile.BadZipFile:
        raise ValueError("downloaded file is not a valid zip — the assembly "
                         "accession may not exist or has no CDS download")
    if len(raw) > _NCBI_CDS_ZIP_MAX_BYTES:
        raise ValueError("CDS member exceeded cap during decompression "
                         "— possible zip-bomb")
    return raw.decode("utf-8", errors="replace")


def _genome_build_codon_table(
        query: str, mode: str = "heg",
        timeout: float = _NCBI_DATASETS_TIMEOUT_S) -> tuple:
    """Build a codon table from an NCBI genome. `query` is an assembly accession
    (`GCF_…`/`GCA_…`) or an NCBI taxid (digits → resolve to the RefSeq reference
    assembly). `mode` is 'heg' or 'genome'. Returns (raw_or_None, message,
    meta_or_None); meta carries accession / taxid / organism / stats for the
    caller's display name + saved entry. Worker-thread only (network + parse)."""
    query = str(query or "").strip()
    if not query:
        return None, "Enter a genome assembly accession (GCF_…) or NCBI taxid", None
    if mode not in ("heg", "genome"):
        return None, f"Invalid mode '{mode}' (expected 'heg' or 'genome')", None
    taxid = ""
    organism = ""
    accession = query
    if query.isdigit():
        taxid = query
        resolved, info = _genome_resolve_reference_accession(query, timeout)
        if resolved is None:
            return None, info, None                # info = error message
        accession, organism = resolved, info       # info = organism name
    # Sanitize before the accession reaches the request URL — blocks path /
    # query injection (e.g. "GCF_1/../x?y=z") from the un-sanitized modal
    # input; the taxid path is already digit-only. `_sanitize_accession`
    # permits the GCF_…/GCA_… charset and rejects URL/shell metacharacters.
    accession = _sanitize_accession(accession) or ""
    if not accession.upper().startswith(("GCF_", "GCA_")):
        return None, (f"'{query}' is not a valid assembly accession "
                      f"(expected GCF_… / GCA_…) or a numeric taxid"), None
    url = (f"{_NCBI_DATASETS_BASE}/genome/accession/{accession}/download"
           f"?include_annotation_type=CDS_FASTA")
    try:
        zip_bytes = _genome_datasets_request(url, _NCBI_CDS_ZIP_MAX_BYTES, timeout)
        fasta = _genome_extract_cds_fasta(zip_bytes)
        raw, stats = _build_heg_table_from_cds(fasta, mode)
    except ValueError as exc:                       # cap / no-CDS / build error
        return None, str(exc), None
    except Exception as exc:                         # network / HTTP
        _log.exception("Genome codon-table build failed for %s", accession)
        return None, f"Download failed for {accession}: {exc}", None
    label = "ribosomal-protein (HEG)" if mode == "heg" else "whole-genome"
    n_cds = stats["n_cds_heg"] if mode == "heg" else stats["n_cds_total"]
    bf = (f", backfilled {', '.join(stats['backfilled'])}"
          if stats.get("backfilled") else "")
    msg = (f"Built {label} table from {accession}"
           + (f" ({organism})" if organism else "")
           + f": {stats['n_codons']:,} codons from {n_cds} CDS, "
           + f"{stats['aa_coverage']}/20 AAs, GC3 {stats['gc3']}%{bf}")
    meta = {"accession": accession, "taxid": taxid,
            "organism": organism, "stats": stats}
    return raw, msg, meta


def _file_build_codon_table(path: "str | Path", mode: str = "heg") -> tuple:
    """Build a codon table from a LOCAL CDS FASTA file — the in-frame
    `cds_from_genomic.fna` shape (nucleotide CDS records), optionally gzipped.
    No network; the offline analogue of `_genome_build_codon_table`. `mode` is
    'heg' (ribosomal-protein bias, for heterologous expression) or 'genome'
    (whole-file average). Returns (raw_or_None, message, meta_or_None); meta
    carries organism (the filename) + stats for the caller's display name + saved
    entry.

    A PROTEOME (amino-acid) FASTA can't be used — codon usage needs the nucleotide
    CDS — so a non-ACGT file surfaces the shared builder's 'no usable CDS' error."""
    if mode not in ("heg", "genome"):
        return None, f"Invalid mode '{mode}' (expected 'heg' or 'genome')", None
    import gzip
    p = _sanitize_path(str(path) if path is not None else None)
    if p is None:
        return None, ("Enter a path to a CDS FASTA file "
                      "(.fna / .fasta, optionally .gz)"), None
    try:
        if not p.is_file():
            return None, f"No such file: {p}", None
        size = p.stat().st_size
        if size > _CDS_FILE_MAX_BYTES:
            return None, (f"File too large ({size:,} bytes > "
                          f"{_CDS_FILE_MAX_BYTES:,})"), None
        data = p.read_bytes()
    except OSError as exc:
        return None, f"Could not read {p.name}: {exc}", None
    if data[:2] == b"\x1f\x8b":                       # gzip magic — NCBI ships CDS .gz
        try:
            data = gzip.decompress(data)
        except (OSError, EOFError) as exc:
            return None, f"Could not gunzip {p.name}: {exc}", None
        if len(data) > _CDS_FILE_MAX_BYTES:
            return None, (f"Decompressed file too large "
                          f"(> {_CDS_FILE_MAX_BYTES:,} bytes)"), None
    try:
        raw, stats = _build_heg_table_from_cds(
            data.decode("utf-8", "replace"), mode)
    except ValueError as exc:                         # no-CDS / no-rprotein / build
        return None, str(exc), None
    except Exception as exc:                          # pragma: no cover — defensive
        _log.exception("Local CDS codon-table build failed for %s", p.name)
        return None, f"Build failed: {exc}", None
    label = "ribosomal-protein (HEG)" if mode == "heg" else "whole-genome"
    n_cds = stats["n_cds_heg"] if mode == "heg" else stats["n_cds_total"]
    bf = (f", backfilled {', '.join(stats['backfilled'])}"
          if stats.get("backfilled") else "")
    msg = (f"Built {label} table from {p.name}: "
           f"{stats['n_codons']:,} codons from {n_cds} CDS, "
           f"{stats['aa_coverage']}/20 AAs, GC3 {stats['gc3']}%{bf}")
    meta = {"accession": "", "taxid": "", "organism": p.stem,
            "source_file": p.name, "stats": stats}
    return raw, msg, meta
