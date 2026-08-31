"""splicecraft_dataaccess — the user-data access layer (Phase D, layer 1).

First piece: `_typed_clone`, the deep-clone-on-read/write helper that enforces
the data layer's SACRED invariant #17 — a caller mutating a returned entry
can't poison the in-memory cache, and a caller mutating its argument after a
save can't leak post-save edits into the next reader. Every `_load_X` returns
`_typed_clone(cache)`; every `_save_X` re-seats the cache via
`_typed_clone(entries)`. It belongs here because that contract IS the data
layer.

This module grows to hold the domain `_load_X` / `_save_X` accessors as Phase D
moves them off the hub, so the modal / screen / agent siblings can import them
here instead of importing the hub (which would be a cycle). It imports ONLY
stdlib + lower/same-layer siblings (splicecraft_state, splicecraft_persistence,
splicecraft_logging); the hub re-exports every name (`from splicecraft_dataaccess
import _typed_clone as _typed_clone` ...) so `sc.<name>` + every existing call
site resolve unchanged.
"""
from __future__ import annotations

import re

from copy import deepcopy
from datetime import date as _date
from typing import Any as _Any, TypeVar as _TypeVar

import splicecraft_state as _state
from splicecraft_logging import _log, _log_event, _repr_for_log
from splicecraft_persistence import _safe_load_json, _safe_save_json
from splicecraft_record import _gb_text_to_record
from splicecraft_biology import _rc
from splicecraft_util import _feat_label, _fuzzy_match, _natural_sort_key, _sanitize_gel_id, _sanitize_plasmid_status

# Immutable types `_typed_clone` returns as-is (no copy needed). Everything else
# recurses (dict / list / tuple) or falls through to `deepcopy` for any
# unexpected type, so the sacred-#17 contract ("caller-side mutations can't
# poison the cache") is preserved end to end.
_IMMUTABLE_CLONE_TYPES = (str, int, float, bool, bytes, type(None))
_TC = _TypeVar("_TC")


def _typed_clone(obj: _TC) -> _TC:
    t = type(obj)
    if t is dict:
        return {k: _typed_clone(v) for k, v in obj.items()}  # type: ignore[return-value]
    if t is list:
        return [_typed_clone(v) for v in obj]  # type: ignore[return-value]
    if t is tuple:
        return tuple(_typed_clone(v) for v in obj)  # type: ignore[return-value]
    if isinstance(obj, _IMMUTABLE_CLONE_TYPES):
        return obj
    return deepcopy(obj)


# ── Feature colours (type → colour map) ─────────────────────────────────────
# First domain accessor cluster moved off the hub (Phase D). Clean closure:
# splicecraft_state (cache + file path + the save lock), splicecraft_persistence
# (safe load/save), the logger, and `_typed_clone` above. No migration / mirror
# / schema validation, so it's the proving cut for the accessor-extraction
# pattern. The save lock is reached as `_state._cache_lock` (the hub keeps a
# same-object `_cache_lock` alias for its own sites).


def _load_feature_colors() -> dict[str, str]:
    """Return the user's customised type → colour map. Missing file / empty
    entries → empty dict. Callers should combine this with
    ``_DEFAULT_TYPE_COLORS`` — that precedence is handled by
    ``_resolve_feature_color`` (hub-side).

    Uses ``_typed_clone`` on read so the cache contract matches every
    other library (invariant #17). Values are ``str`` (immutable) today,
    so the practical risk is nil — but using the canonical helper keeps
    the pattern honest in case a future schema bump adds a nested
    structure to the value side.
    """
    # Sweep #26: double-checked locking — cache-hit fast path is lock-free.
    cached = _state._feature_colors_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._feature_colors_cache is None:
            entries, warning = _safe_load_json(_state._FEATURE_COLORS_FILE, "Feature colors")
            if warning:
                _log.warning(warning)
            result: dict[str, str] = {}
            for e in entries:
                if not isinstance(e, dict):
                    continue
                ft  = e.get("feature_type")
                col = e.get("color")
                if isinstance(ft, str) and ft and isinstance(col, str) and col:
                    result[ft] = col
            _state._feature_colors_cache = _typed_clone(result)
        return _typed_clone(_state._feature_colors_cache)


def _save_feature_colors(mapping: dict[str, str]) -> None:
    """Persist the type → colour map. Written as a list of ``{"feature_type":
    ..., "color": ...}`` dicts so it shares the schema-envelope shape with
    the other libraries (sacred invariant #7).

    Re-seats the cache via ``_typed_clone`` so a caller that keeps
    mutating its mapping after the save doesn't leak post-save edits
    into the next reader (invariant #17, full deepcopy-on-save side).
    """
    entries = [{"feature_type": ft, "color": col}
               for ft, col in mapping.items()]
    with _state._cache_lock:
        _safe_save_json(_state._FEATURE_COLORS_FILE, entries, "Feature colors")
        _state._feature_colors_cache = _typed_clone(mapping)


# ── Enzyme collections ──────────────────────────────────────────────────────


def _load_enzyme_collections() -> list[dict]:
    """Deep-copy on read so callers can mutate freely; pitfall #17.

    Sweep #26: double-checked locking — cache-hit fast path is lock-free."""
    cached = _state._enzyme_collections_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._enzyme_collections_cache is None:
            entries, warning = _safe_load_json(
                _state._ENZYME_COLLECTIONS_FILE, "Enzyme collections",
            )
            if warning:
                _log.warning(warning)
            _state._enzyme_collections_cache = [
                e for e in entries if isinstance(e, dict)
            ]
        return _typed_clone(_state._enzyme_collections_cache)


def _save_enzyme_collections(entries: list[dict]) -> None:
    with _state._cache_lock:
        _safe_save_json(
            _state._ENZYME_COLLECTIONS_FILE, entries, "Enzyme collections",
        )
        _state._enzyme_collections_cache = _typed_clone(entries)


# ── Experiment projects ─────────────────────────────────────────────────────


def _load_experiment_projects() -> "list[dict]":
    """Return a deepcopy of the experiment-projects list so callers can
    mutate entries (rename, edit experiments list) without poisoning
    the in-memory cache (invariant #17).

    Sweep #26: double-checked locking — cache-hit fast path is lock-free.
    """
    cached = _state._experiment_projects_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._experiment_projects_cache is None:
            entries, warning = _safe_load_json(
                _state._EXPERIMENT_PROJECTS_FILE, "Experiment projects",
            )
            if warning:
                _log.warning(warning)
            _state._experiment_projects_cache = [
                e for e in entries if isinstance(e, dict)
            ]
    return _typed_clone(_state._experiment_projects_cache)


def _save_experiment_projects(entries: "list[dict]") -> None:
    """Persist the full experiment-projects list through the four-layer
    data-safety net (invariant #31). Deepcopies into the cache so caller
    mutations after save can't poison subsequent loaders (invariant #17),
    under the save lock (invariant #41 — concurrency)."""
    with _state._cache_lock:
        _safe_save_json(
            _state._EXPERIMENT_PROJECTS_FILE, entries, "Experiment projects",
        )
        _state._experiment_projects_cache = _typed_clone(entries)


# ── Gels ────────────────────────────────────────────────────────────────────


def _load_gels() -> "list[dict]":
    """Cached + deepcopy-on-read load (invariant #17). Filters non-dict
    entries defensively (hand-edited JSON / schema drift). Sweep #26
    double-checked locking: the cache-hit fast path runs without the lock
    so a worker holding it can't freeze every UI-thread `_load_gels`; only
    the cold populate-from-disk path acquires it, double-checking after."""
    cached = _state._gels_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._gels_cache is None:
            entries, warning = _safe_load_json(_state._GELS_FILE, "Gels")
            if warning:
                _log.warning(warning)
            _state._gels_cache = [e for e in entries if isinstance(e, dict)]
        return _typed_clone(_state._gels_cache)


def _save_gels(entries: "list[dict]") -> None:
    """Persist through `_safe_save_json` (atomic + four-layer data-safety).
    Takes the save lock so concurrent saves can't desync cache vs disk
    (invariant #41)."""
    with _state._cache_lock:
        _safe_save_json(_state._GELS_FILE, entries, "Gels")
        _state._gels_cache = _typed_clone(entries)


# ── Protein motifs (builtin defaults merged with user edits on read) ──────────
# `_PROTEIN_MOTIFS` (the builtin motif library) is used only by
# `_load_protein_motifs` (merge with the user's protein_motifs.json), so it
# travels here with the accessors.


# Builtin protein-motif / domain library — common tags + linkers +
# protease sites a synthetic-biology user reaches for when composing
# a recombinant protein. Stored as a module-level list so it's
# always available without needing a user-managed file. Entry shape
# mirrors the DNA feature library so the side-panel renderer can
# reuse the same row-build code.
# Sweep #16 (2026-05-21): each built-in motif carries its own unique
# ``color`` hex so the user can tell them apart at a glance in the
# motif library AND in the dithered feature bars above the AA row.
# The palette spreads visually within each family (e.g. Tags shift
# from deep blue → cyan → indigo → violet so a side-by-side stack of
# His6 + FLAG + Myc + V5 reads as four obviously different bars).
# User edits via the Edit button persist to `protein_motifs.json`
# and override these defaults.
_PROTEIN_MOTIFS: list[dict] = [
    # ── Affinity tags (cool spectrum: blue → cyan → indigo → violet) ──
    {"name": "His6",        "feature_type": "Tag",
     "sequence": "HHHHHH", "color": "#1E40AF",
     "description": "Hexahistidine affinity tag (Ni-NTA / IMAC purification)."},
    {"name": "His8",        "feature_type": "Tag",
     "sequence": "HHHHHHHH", "color": "#3B82F6",
     "description": "Octahistidine — tighter Ni-NTA binding than 6xHis."},
    {"name": "His10",       "feature_type": "Tag",
     "sequence": "HHHHHHHHHH", "color": "#60A5FA",
     "description": "Decahistidine — even tighter, for harsh-wash purification."},
    {"name": "FLAG",        "feature_type": "Tag",
     "sequence": "DYKDDDDK", "color": "#0E7490",
     "description": "FLAG tag (anti-FLAG M2 affinity purification)."},
    {"name": "3xFLAG",      "feature_type": "Tag",
     "sequence": "DYKDHDGDYKDHDIDYKDDDDK", "color": "#06B6D4",
     "description": "Triple FLAG — higher sensitivity for low-expression targets."},
    {"name": "HA",          "feature_type": "Tag",
     "sequence": "YPYDVPDYA", "color": "#67E8F9",
     "description": "Influenza hemagglutinin epitope (anti-HA antibodies)."},
    {"name": "Myc",         "feature_type": "Tag",
     "sequence": "EQKLISEEDL", "color": "#4338CA",
     "description": "c-Myc epitope tag (9E10 antibody)."},
    {"name": "V5",          "feature_type": "Tag",
     "sequence": "GKPIPNPLLGLDST", "color": "#6366F1",
     "description": "V5 epitope tag (paramyxovirus)."},
    {"name": "Strep-II",    "feature_type": "Tag",
     "sequence": "WSHPQFEK", "color": "#A78BFA",
     "description": "Strep-Tactin affinity tag (mild biotin elution)."},
    {"name": "T7",          "feature_type": "Tag",
     "sequence": "MASMTGGQQMG", "color": "#7C3AED",
     "description": "T7 leader peptide (anti-T7 monoclonal)."},
    # ── Localisation signals (warm yellows / golds) ────────────────
    {"name": "NLS (SV40)",  "feature_type": "Signal",
     "sequence": "PKKKRKV", "color": "#EAB308",
     "description": "Classical SV40 large T-antigen nuclear localisation signal."},
    {"name": "NLS (bipartite)", "feature_type": "Signal",
     "sequence": "KRPAATKKAGQAKKKK", "color": "#FBBF24",
     "description": "Nucleoplasmin bipartite NLS."},
    {"name": "NES",         "feature_type": "Signal",
     "sequence": "LPPLERLTL", "color": "#F59E0B",
     "description": "HIV-Rev nuclear export signal (CRM1-dependent)."},
    # ── Linkers (neutral greys, light → dark, warm at the end) ─────
    {"name": "GSG",         "feature_type": "Linker",
     "sequence": "GSG", "color": "#94A3B8",
     "description": "Minimal flexible linker."},
    {"name": "GGS",         "feature_type": "Linker",
     "sequence": "GGS", "color": "#64748B",
     "description": "Short flexible linker."},
    {"name": "(GGGGS)x3",   "feature_type": "Linker",
     "sequence": "GGGGSGGGGSGGGGS", "color": "#475569",
     "description": "Classic flexible linker for scFv / fusion proteins."},
    {"name": "(GGGGS)x4",   "feature_type": "Linker",
     "sequence": "GGGGSGGGGSGGGGSGGGGS", "color": "#57534E",
     "description": "Longer flexible linker for domain separation."},
    {"name": "EAAAK x3",    "feature_type": "Linker",
     "sequence": "EAAAKEAAAKEAAAK", "color": "#78716C",
     "description": "Rigid α-helical linker."},
    # ── Protease cleavage sites (reds → pinks) ──────────────────────
    {"name": "TEV",         "feature_type": "Cleavage",
     "sequence": "ENLYFQG", "color": "#B91C1C",
     "description": "TEV protease site (cuts between Q and G)."},
    {"name": "PreScission", "feature_type": "Cleavage",
     "sequence": "LEVLFQGP", "color": "#DC2626",
     "description": "HRV 3C / PreScission protease site (cuts between Q and G)."},
    {"name": "Thrombin",    "feature_type": "Cleavage",
     "sequence": "LVPRGS", "color": "#F87171",
     "description": "Thrombin cleavage site."},
    {"name": "Factor Xa",   "feature_type": "Cleavage",
     "sequence": "IEGR", "color": "#EC4899",
     "description": "Factor Xa protease site."},
    {"name": "Furin",       "feature_type": "Cleavage",
     "sequence": "RRRR", "color": "#DB2777",
     "description": "Furin recognition site (R-X-K/R-R minimal)."},
    # ── Self-cleaving 2A peptides (greens) ─────────────────────────
    {"name": "P2A",         "feature_type": "2A",
     "sequence": "GSGATNFSLLKQAGDVEENPGP", "color": "#15803D",
     "description": "Porcine teschovirus 2A self-cleaving peptide."},
    {"name": "T2A",         "feature_type": "2A",
     "sequence": "GSGEGRGSLLTCGDVEENPGP", "color": "#16A34A",
     "description": "Thosea asigna 2A peptide."},
    {"name": "E2A",         "feature_type": "2A",
     "sequence": "GSGQCTNYALLKLAGDVESNPGP", "color": "#22C55E",
     "description": "Equine rhinitis A 2A peptide."},
    {"name": "F2A",         "feature_type": "2A",
     "sequence": "GSGVKQTLNFDLLKLAGDVESNPGP", "color": "#10B981",
     "description": "Foot-and-mouth-disease virus 2A peptide."},
    # ── Common functional motifs (fuchsias / magentas) ─────────────
    {"name": "Kozak start", "feature_type": "Motif",
     "sequence": "M", "color": "#C026D3",
     "description": "Start codon (methionine) — required N-terminal."},
    {"name": "Stop",        "feature_type": "Motif",
     "sequence": "*", "color": "#A21CAF",
     "description": "Translation stop codon."},
    {"name": "FLAG+Stop",   "feature_type": "Motif",
     "sequence": "DYKDDDDK*", "color": "#D946EF",
     "description": "FLAG tag followed by stop — quick C-terminal tagging."},
]


def _load_protein_motifs() -> list[dict]:
    """Return the merged protein-motif library: built-in entries
    overridden by user edits stored in `protein_motifs.json`. Deep-
    copied so callers can mutate freely.

    Sweep #26: double-checked locking — cache-hit fast path is lock-free."""
    cached = _state._protein_motifs_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._protein_motifs_cache is None:
            user_entries, warning = _safe_load_json(
                _state._PROTEIN_MOTIFS_FILE, "Protein motifs",
            )
            if warning:
                _log.warning(warning)
            user_entries = [e for e in user_entries if isinstance(e, dict)]
            user_by_name: dict[str, dict] = {
                str(e.get("name") or ""): e for e in user_entries if e.get("name")
            }
            merged: list[dict] = []
            builtin_names: set[str] = set()
            for builtin in _PROTEIN_MOTIFS:
                name = str(builtin.get("name") or "")
                builtin_names.add(name)
                merged.append(user_by_name.get(name, dict(builtin)))
            # User-added novel motifs (name NOT in builtins) append in
            # insertion order so user-defined items land predictably.
            for e in user_entries:
                name = str(e.get("name") or "")
                if name and name not in builtin_names:
                    merged.append(dict(e))
            _state._protein_motifs_cache = merged
        return _typed_clone(_state._protein_motifs_cache)


def _save_protein_motifs(entries: list[dict]) -> None:
    """Persist `entries` (user-modified motifs only — not the full
    merged list) to `protein_motifs.json`. Caller passes the list of
    entries to STORE; the merge with built-ins happens on read."""
    with _state._cache_lock:
        _safe_save_json(_state._PROTEIN_MOTIFS_FILE, entries, "Protein motifs")
        # Invalidate the cache; next _load_protein_motifs rebuilds
        # the merged list. Reseating with the user-only list would
        # leave the merged form stale.
        _state._protein_motifs_cache = None


# ── Codon-table registry (MISSION-CRITICAL: _codon_optimize depends on it) ──
# The genetic-code map + builtin K12 table + the raw<->json converters +
# the load/save accessors. `_codon_raw_from_json` forces each codon's AA to
# the standard genetic code (never the stored label) so a corrupt table
# can't feed a wrong AA->codon grouping into the optimizer. `_CODON_GENETIC_CODE`
# is also used hub-side (optimizer + translation, 17 refs) via the re-export.


_CODON_GENETIC_CODE: dict[str, str] = {
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


_CODON_BUILTIN_K12: dict[str, tuple[str, int]] = {
    "GGG": ("G",  44), "GGA": ("G",  47), "GGT": ("G", 109), "GGC": ("G", 171),
    "GAG": ("E",  94), "GAA": ("E", 224), "GAT": ("D", 194), "GAC": ("D", 105),
    "GTG": ("V", 135), "GTA": ("V",  59), "GTT": ("V",  86), "GTC": ("V",  60),
    "GCG": ("A", 197), "GCA": ("A", 108), "GCT": ("A",  55), "GCC": ("A", 162),
    "AGG": ("R",   8), "AGA": ("R",   7), "AGT": ("S",  37), "AGC": ("S",  85),
    "AAG": ("K",  62), "AAA": ("K", 170), "AAT": ("N", 112), "AAC": ("N", 125),
    "ATG": ("M", 127), "ATA": ("I",  19), "ATT": ("I", 156), "ATC": ("I",  93),
    "ACG": ("T",  59), "ACA": ("T",  33), "ACT": ("T",  41), "ACC": ("T", 117),
    "TGG": ("W",  55), "TGT": ("C",  30), "TGC": ("C",  41),
    "TAT": ("Y",  86), "TAC": ("Y",  75),
    "TTG": ("L",  61), "TTA": ("L",  78), "TTT": ("F", 101), "TTC": ("F",  77),
    "TCG": ("S",  41), "TCA": ("S",  40), "TCT": ("S",  29), "TCC": ("S",  28),
    "CGG": ("R",  21), "CGA": ("R",  22), "CGT": ("R", 108), "CGC": ("R", 133),
    "CAG": ("Q", 142), "CAA": ("Q",  62), "CAT": ("H",  81), "CAC": ("H",  67),
    "CTG": ("L", 240), "CTA": ("L",  27), "CTT": ("L",  61), "CTC": ("L",  54),
    "CCG": ("P", 137), "CCA": ("P",  34), "CCT": ("P",  43), "CCC": ("P",  33),
    "TAA": ("*",   9), "TAG": ("*",   0), "TGA": ("*",   5),
}


def _codon_raw_to_json(raw: dict) -> dict:
    """Convert in-memory {codon: (aa, count)} → JSON-safe {codon: [aa, count]}."""
    return {c: [aa, int(n)] for c, (aa, n) in raw.items()}


def _codon_raw_from_json(blob: dict) -> dict:
    """Inverse of `_codon_raw_to_json`. Accepts tuples or 2-item lists.

    HARDENED (2026-06-12): the codon KEY must be a real codon and its
    amino acid is forced to the standard genetic code — never trusted from
    the stored label. A corrupted or hand-edited `codon_tables.json` (e.g.
    `ATG` mislabelled `L`, or a 2-nt key) would otherwise feed a wrong
    AA→codon grouping straight into `_codon_optimize`, silently breaking
    the round-trip / no-codon→AA-mismatch invariant the optimizer
    guarantees (a mission-critical, zero-tolerance subsystem). This mirrors
    the validation the interactive TSV importer already enforces
    (`_parse_codon_tsv`). Unparseable rows are dropped (not raised) so a
    partly-corrupt table still loads its good rows — and a dropped codon
    surfaces loudly later as a "missing amino acid" optimizer error rather
    than as silent wrong DNA."""
    out: dict = {}
    if not isinstance(blob, dict):
        return out
    for c, v in blob.items():
        if not (isinstance(v, (list, tuple)) and len(v) == 2):
            continue
        codon = str(c).upper().replace("U", "T")
        aa = _CODON_GENETIC_CODE.get(codon)
        if aa is None:
            continue                      # bad length / non-ACGT / unknown codon
        try:
            count = int(v[1])
        except (TypeError, ValueError):
            continue
        if count < 0:
            continue
        out[codon] = (aa, count)
    return out


def _codon_tables_load() -> list[dict]:
    """Load codon-table registry. Returns a deep copy so caller-side
    mutation of the `raw` dict (e.g. an in-progress Kazusa import edit)
    can't poison the cache — matches invariant #17. Seeds built-in
    E. coli K12 on first run so the library is never empty."""
    # Double-checked locking (match `_load_features`): lock-free cache hit; the
    # populate (incl. the K12 seed + its save) holds `_cache_lock` so it can't
    # race a concurrent `_codon_tables_save` (agent thread) and lose the write.
    cached = _state._codon_tables_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._codon_tables_cache is None:
            entries, warning = _safe_load_json(_state._CODON_TABLES_FILE, "Codon table library")
            if warning:
                _log.warning("Codon table library: %s", warning)
            fixed: list = []
            for e in entries:
                if not isinstance(e, dict):
                    continue
                raw = _codon_raw_from_json(e.get("raw", {}))
                if not raw:
                    continue
                fixed.append({
                    "name":   e.get("name", "?"),
                    "taxid":  str(e.get("taxid", "")),
                    "source": e.get("source", "user"),
                    "added":  e.get("added", ""),
                    "raw":    raw,
                })
            # Seed built-in K12 if not present
            if not any(e.get("taxid") == "83333" for e in fixed):
                fixed.insert(0, {
                    "name":   "E. coli K12",
                    "taxid":  "83333",
                    "source": "builtin",
                    "added":  _date.today().isoformat(),
                    "raw":    dict(_CODON_BUILTIN_K12),
                })
                _state._codon_tables_cache = fixed
                # Module-level seed path — no `app` context to notify. Log-only
                # so a disk-full first launch surfaces in the log bundle and
                # the cache (in memory) still has the K12 seed.
                try:
                    _codon_tables_save(fixed)
                except (OSError, RuntimeError):
                    _log.exception("Codon tables: K12 seed save failed (in-memory "
                                   "cache populated; disk write deferred)")
            else:
                _state._codon_tables_cache = fixed
        return _typed_clone(_state._codon_tables_cache)


def _codon_tables_save(entries: list[dict]) -> None:
    """Persist registry to disk via _safe_save_json (atomic, .bak).
    Re-seats the cache with a deep copy so callers that retain a reference
    to the input list (and continue mutating it) can't leak post-save
    edits into the next reader. Matches invariant #17."""
    serializable = [{
        "name":   e.get("name", "?"),
        "taxid":  str(e.get("taxid", "")),
        "source": e.get("source", "user"),
        "added":  e.get("added", ""),
        "raw":    _codon_raw_to_json(e.get("raw", {})),
    } for e in entries]
    with _state._cache_lock:
        _safe_save_json(_state._CODON_TABLES_FILE, serializable, "Codon table library")
        _state._codon_tables_cache = _typed_clone(entries)


# ── HMM database catalog (pyhmmer profile DBs) ──────────────────────────────
# Builtin catalog + id/url sanitisers + entry normaliser + load/save. The
# sanitisers + normaliser are also used by the hub HMM-management handlers
# (via the re-export).


# Custom URL accepted only if it parses as http/https with a host;
# the protocol allowlist matches `[INV-36]`'s `$SPLICECRAFT_PYPI_URL`
# override pattern.
_HMM_DB_URL_MAX_LEN = 2048


_BUILTIN_HMM_DB_CATALOG: "tuple[dict, ...]" = (
    {
        "id":          "pfam-a",
        "name":        "Pfam-A",
        "url":         "https://ftp.ebi.ac.uk/pub/databases/Pfam/"
                       "current_release/Pfam-A.hmm.gz",
        "version_url": "https://ftp.ebi.ac.uk/pub/databases/Pfam/"
                       "current_release/Pfam.version.gz",
        "format":      "hmm-gz",
        "builtin":     True,
        "description": "Pfam-A: curated protein family HMMs from EBI "
                       "(~300 MB download, ~3 GB on disk after hmmpress).",
    },
    {
        "id":          "ncbifam",
        "name":        "NCBIfam",
        "url":         "https://ftp.ncbi.nlm.nih.gov/hmm/current/"
                       "NCBIfam.HMM.gz",
        # NCBIfam doesn't publish a one-line version file; we use
        # the source file's HTTP Last-Modified header as the version
        # signature via the version-check helper's fallback path.
        "version_url": "",
        "format":      "hmm-gz",
        "builtin":     True,
        "description": "NCBIfam: HMMs for prokaryotic protein families "
                       "(~600 MB download, ~4 GB on disk).",
    },
)


def _typed_clone_hmm_catalog(entries: "list | None") -> list:
    """Deepcopy-on-read for the catalog cache. Mirrors the pattern
    every other persisted-list helper uses (cf. `[PIT-17]`)."""
    if entries is None:
        return []
    return [dict(e) for e in entries if isinstance(e, dict)]


def _sanitize_hmm_db_id(raw: object) -> "str | None":
    """Return a sanitised id, or None on reject. Same rule set as
    `_sanitize_experiment_id` / `_sanitize_gel_id`: non-empty, no NUL,
    no path separators / shell metas, no `..`, ≤64 chars. The id
    doubles as the on-disk directory name."""
    if not isinstance(raw, str):
        return None
    s = raw.strip()
    if not s or len(s) > 64:
        return None
    if "\x00" in s or ".." in s:
        return None
    if "/" in s or "\\" in s:
        return None
    for ch in s:
        # Allow ASCII letters, digits, underscore, hyphen ONLY.
        # Python's `str.isalnum()` is Unicode-permissive (accepts
        # `é`, `ñ`, CJK, etc.) — we explicitly check the codepoint
        # range so the id stays safe as a filesystem path component
        # on every OS. Sweep #28 hardening: a unicode-looking id
        # could pass cross-platform path checks on one OS but break
        # on another, or differ between the bytes-on-disk and what
        # the user typed (NFC vs NFD normalisation).
        cp = ord(ch)
        is_lower = (0x61 <= cp <= 0x7A)   # a-z
        is_upper = (0x41 <= cp <= 0x5A)   # A-Z
        is_digit = (0x30 <= cp <= 0x39)   # 0-9
        is_safe  = (ch == "_" or ch == "-")
        if not (is_lower or is_upper or is_digit or is_safe):
            return None
    return s


def _sanitize_hmm_db_url(raw: object) -> "str | None":
    """http(s) URL ≤ `_HMM_DB_URL_MAX_LEN`, otherwise None. Mirrors
    `[INV-36]`'s `$SPLICECRAFT_PYPI_URL` validator.

    Sweep #28: reject if the RAW input carries any whitespace /
    control char ANYWHERE — even at the trailing edge. A leading /
    trailing newline on a copy-pasted URL could trick a downstream
    parser into header-smuggling. The check fires before `.strip()`
    so we surface "you pasted a URL with hidden whitespace" rather
    than silently scrubbing it."""
    if not isinstance(raw, str):
        return None
    if any(ch.isspace() or ord(ch) < 0x20 for ch in raw):
        return None
    if not raw or len(raw) > _HMM_DB_URL_MAX_LEN:
        return None
    if not (raw.startswith("http://") or raw.startswith("https://")):
        return None
    return raw


def _normalise_hmm_db_entry(entry: dict, *, builtin_default: bool = False
                              ) -> "dict | None":
    """Return a fresh dict shaped to the catalog schema, or None if
    the entry is unsalvageable. Idempotent — safe to apply to entries
    loaded from disk and to entries built in-process."""
    if not isinstance(entry, dict):
        return None
    entry_id = _sanitize_hmm_db_id(entry.get("id"))
    if entry_id is None:
        return None
    url = _sanitize_hmm_db_url(entry.get("url"))
    if url is None:
        return None
    name_raw = entry.get("name")
    if not isinstance(name_raw, str) or not name_raw.strip():
        name = entry_id
    else:
        name = name_raw.strip()[:200]
    version_url_raw = entry.get("version_url")
    if isinstance(version_url_raw, str) and version_url_raw.strip():
        version_url = _sanitize_hmm_db_url(version_url_raw) or ""
    else:
        version_url = ""
    fmt = entry.get("format") or "hmm-gz"
    if fmt not in ("hmm-gz", "hmm"):
        fmt = "hmm-gz"
    desc_raw = entry.get("description")
    description = (desc_raw.strip()[:500]
                   if isinstance(desc_raw, str) else "")
    return {
        "id":          entry_id,
        "name":        name,
        "url":         url,
        "version_url": version_url,
        "format":      fmt,
        "builtin":     bool(entry.get("builtin", builtin_default)),
        "description": description,
    }


def _load_hmm_db_catalog() -> list:
    """Return the registered HMM-DB catalog (builtins + user-added).
    Deep-clone on read so callers can freely mutate without affecting
    the cache (cf. `[PIT-17]`).

    First-load merges the built-in defaults into whatever is on disk;
    subsequent loads return whatever was saved. Built-ins are still
    re-injected if the user has removed them (defensive — the user
    UI doesn't let you delete a builtin, but a hand-edited
    catalog.json shouldn't permanently nuke them either).
    """
    with _state._cache_lock:
        if _state._hmm_db_catalog_cache is None:
            entries, warning = _safe_load_json(
                _state._HMM_DB_CATALOG_FILE, "HMM database catalog",
            )
            if warning:
                _log.warning("HMM database catalog: %s", warning)
            normalised: list[dict] = []
            seen_ids: set[str] = set()
            for entry in entries:
                norm = _normalise_hmm_db_entry(entry)
                if norm is None:
                    continue
                if norm["id"] in seen_ids:
                    continue
                seen_ids.add(norm["id"])
                normalised.append(norm)
            # Re-inject any built-in the user removed (or that wasn't
            # there on first launch). User overrides for builtin URL
            # are preserved.
            for builtin in _BUILTIN_HMM_DB_CATALOG:
                if builtin["id"] not in seen_ids:
                    normalised.append(
                        _normalise_hmm_db_entry(
                            dict(builtin), builtin_default=True,
                        ) or {}
                    )
                    seen_ids.add(builtin["id"])
            _state._hmm_db_catalog_cache = normalised
        return _typed_clone_hmm_catalog(_state._hmm_db_catalog_cache)


def _save_hmm_db_catalog(entries: list) -> None:
    """Persist the catalog. Routed through `_safe_save_json` for the
    L2 chokepoint + atomic-write + .bak chain. Cache reseat under
    `_cache_lock` per `[INV-50]`."""
    cleaned: list[dict] = []
    seen: set[str] = set()
    for entry in entries:
        n = _normalise_hmm_db_entry(entry)
        if n is None or n["id"] in seen:
            continue
        seen.add(n["id"])
        cleaned.append(n)
    with _state._cache_lock:
        _safe_save_json(
            _state._HMM_DB_CATALOG_FILE, cleaned, "HMM database catalog",
        )
        _state._hmm_db_catalog_cache = _typed_clone_hmm_catalog(cleaned)


# ── Custom restriction enzymes ──────────────────────────────────────────────
# `_save_custom_enzymes` triggers a DOMAIN side-effect — rebuild the restriction
# `_SCAN_CATALOG` + bust the enzyme caches so a new / redefined custom enzyme
# shows up in scans (a silent-biology hazard if skipped). That side-effect lives
# hub-side (it reaches the scanner) and runs via the registered hook below.


def _load_custom_enzymes() -> list[dict]:
    """Deep-copy on read (pitfall #17).

    Sweep #26: double-checked locking — cache-hit fast path is
    lock-free so an unrelated lock-holder can't freeze UI reads."""
    cached = _state._custom_enzymes_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._custom_enzymes_cache is None:
            entries, warning = _safe_load_json(
                _state._CUSTOM_ENZYMES_FILE, "Custom enzymes",
            )
            if warning:
                _log.warning(warning)
            _state._custom_enzymes_cache = [
                e for e in entries if isinstance(e, dict)
            ]
        return _typed_clone(_state._custom_enzymes_cache)


def _save_custom_enzymes(entries: list[dict]) -> None:
    with _state._cache_lock:
        _safe_save_json(
            _state._CUSTOM_ENZYMES_FILE, entries, "Custom enzymes",
        )
        _state._custom_enzymes_cache = _typed_clone(entries)
    # Rebuild `_SCAN_CATALOG` + bust the enzyme caches so the new custom enzyme
    # shows up in the next restriction scan — hub-side effect, fired via the
    # registered hook. None only during the import window (no save runs then).
    hook = getattr(_state, "_after_custom_enzyme_save_hook", None)
    if hook is not None:
        hook()


# ── Entry vectors (Golden-Braid / MoClo acceptor plasmids) ──────────────────
# `_save_entry_vectors`' post-save side-effect (bust the EV digest + acceptor-TU
# caches) stays hub-side via `_state._after_entry_vectors_save_hook`. The
# name-trim backfill is a pure helper that travels with the cluster.


def _backfill_entry_vector_names(
        entries: list[dict]) -> "tuple[list[dict], int]":
    """Trim leading / trailing whitespace from each entry vector's
    `name` field. Returns `(entries, n_changed)`. Idempotent.

    **Hardening:** only mutates entries that already HAVE a `name`
    key — never adds the field where missing. Skips non-dict entries
    + non-string names (defensive against hand-edited JSON)."""
    n_changed = 0
    for e in entries:
        if not isinstance(e, dict) or "name" not in e:
            continue
        raw = e["name"]
        if not isinstance(raw, str):
            continue
        trimmed = raw.strip()
        if trimmed != raw:
            e["name"] = trimmed
            n_changed += 1
    return entries, n_changed


def _load_entry_vectors() -> list[dict]:
    with _state._cache_lock:
        if _state._entry_vectors_cache is None:
            entries, warning = _safe_load_json(
                _state._ENTRY_VECTORS_FILE, "Entry vectors",
            )
            if warning:
                _log.warning(warning)
            _state._entry_vectors_cache = [
                e for e in entries
                if isinstance(e, dict)
                and isinstance(e.get("grammar_id"), str)
            ]
        if not _state._entry_vectors_name_trim_done:
            _state._entry_vectors_cache, n_changed = (
                _backfill_entry_vector_names(_state._entry_vectors_cache)
            )
            _state._entry_vectors_name_trim_done = True
            if n_changed > 0:
                _log.info(
                    "entry-vectors: name trim backfill cleaned %d "
                    "entries", n_changed,
                )
                _log_event(
                    "entry_vectors.name_trim_backfill",
                    n_changed=n_changed,
                )
                save_target = _state._entry_vectors_cache
            else:
                save_target = None
        else:
            save_target = None
    if save_target is not None:
        try:
            _save_entry_vectors(save_target)
        except (OSError, RuntimeError):
            _log.exception(
                "entry-vectors: name trim backfill save failed",
            )
    assert _state._entry_vectors_cache is not None
    return _typed_clone(_state._entry_vectors_cache)


def _save_entry_vectors(entries: list[dict]) -> None:
    with _state._cache_lock:
        _safe_save_json(_state._ENTRY_VECTORS_FILE, entries, "Entry vectors")
        _state._entry_vectors_cache = _typed_clone(entries)
    # Bust the EV digest + acceptor-TU caches (hub-side: a reconfigured vector
    # set changes which overhangs/stuffers match) via the registered hook.
    hook = getattr(_state, "_after_entry_vectors_save_hook", None)
    if hook is not None:
        hook()


# ── Cloning grammars (Golden-Braid L0 / MoClo level rules) ──────────────────
# `_save_custom_grammars`' post-save side-effect (bust the assembly-fragment +
# EV-role-detect caches) stays hub-side via `_state._after_custom_grammars_save_hook`.


def _load_custom_grammars() -> list[dict]:
    """Sweep #26: double-checked locking — cache-hit fast path is lock-free."""
    cached = _state._grammars_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._grammars_cache is None:
            entries, warning = _safe_load_json(_state._GRAMMARS_FILE, "Cloning grammars")
            if warning:
                _log.warning(warning)
            entries = [e for e in entries if isinstance(e, dict) and isinstance(e.get("id"), str)]
            _state._grammars_cache = entries
        return _typed_clone(_state._grammars_cache)


def _save_custom_grammars(entries: list[dict]) -> None:
    with _state._cache_lock:
        _safe_save_json(_state._GRAMMARS_FILE, entries, "Cloning grammars")
        _state._grammars_cache = _typed_clone(entries)
    # Bust the assembly-fragment digest + EV-role-detect caches (hub-side: a
    # grammar enzyme / level_up_enzyme change shifts fragment overhangs + role
    # detection) via the registered hook.
    hook = getattr(_state, "_after_custom_grammars_save_hook", None)
    if hook is not None:
        hook()


# ── Primers + primer collections (re-entrant mirror, INV-50) ────────────────
# The primer library + its collections layer + the active-collection mirror +
# the dedup/name-trim helpers move as ONE group: `_save_primers` mirrors into
# the active collection INSIDE the `_state._cache_lock` block (so the mirror
# can't drift, sweep #10 / INV-50) and the mirror re-enters `_save_primer_
# collections` via the re-entrant RLock — which only works with the whole
# group sibling-internal. The active-collection NAME (a setting) is resolved
# hub-side via `_state._active_primer_collection_name_hook` to avoid pulling
# the settings layer in.


def _backfill_primer_names(
        entries: list[dict]) -> "tuple[list[dict], int]":
    """Trim leading / trailing whitespace from each primer entry's
    `name` field. Returns `(entries, n_changed)`. Idempotent.
    Sequence + tm fields are left untouched — only `name` is the
    user-facing label that participates in `==` matching.

    **Hardening:** only mutates entries that already HAVE a `name`
    key — never adds the field where missing. Skips non-dict entries
    + non-string names (defensive against hand-edited JSON)."""
    n_changed = 0
    for e in entries:
        if not isinstance(e, dict) or "name" not in e:
            continue
        raw = e["name"]
        if not isinstance(raw, str):
            continue
        trimmed = raw.strip()
        if trimmed != raw:
            e["name"] = trimmed
            n_changed += 1
    return entries, n_changed


def _load_primers() -> list[dict]:
    """Load + cache primers, with one-shot name-trim backfill on
    first read (PIT-36 whitespace-trim invariant applied to the
    primer subsystem too). Deep-copy on read per invariant #17."""
    with _state._cache_lock:
        if _state._primers_cache is None:
            entries, warning = _safe_load_json(
                _state._PRIMERS_FILE, "Primer library",
            )
            if warning:
                _log.warning(warning)
            _state._primers_cache = [e for e in entries if isinstance(e, dict)]
        if not _state._primers_name_trim_done:
            _state._primers_cache, n_changed = _backfill_primer_names(
                _state._primers_cache,
            )
            _state._primers_name_trim_done = True
            if n_changed > 0:
                _log.info(
                    "primers: name trim backfill cleaned %d entries",
                    n_changed,
                )
                _log_event(
                    "primers.name_trim_backfill", n_changed=n_changed,
                )
                save_target = _state._primers_cache
            else:
                save_target = None
        else:
            save_target = None
    if save_target is not None:
        try:
            _save_primers(save_target)
        except (OSError, RuntimeError):
            _log.exception(
                "primers: name trim backfill save failed",
            )
    assert _state._primers_cache is not None
    return _typed_clone(_state._primers_cache)


def _dedupe_primers_by_sequence(entries: list[dict]) -> list[dict]:
    """Collapse primer entries whose ``sequence`` field duplicates an
    earlier entry's sequence (case-insensitive). The first entry of
    each unique sequence wins — callers that prepend MRU at index 0
    therefore keep their newest copy, while older duplicates are
    silently dropped.

    The primer-library policy across the app has always been "one entry
    per unique sequence":

      * ``PrimerDesignScreen._save_primers_btn`` rejects a designed
        primer when its sequence already exists.
      * ``_apply_record`` and ``_bulk_import_folder`` dedupe
        .dna-imported primers against the existing library before
        appending.

    Pre-2026-05-10 those dedupe paths only filtered NEW additions —
    duplicates that snuck into ``primers.json`` from earlier sessions
    (manual edits, imports before the dedupe landed, etc.) stayed
    forever. Calling this helper on every save means the next write
    of ``primers.json`` cleans up legacy duplicates without the user
    having to do anything.

    Entries lacking a usable ``sequence`` (None / empty / non-string)
    are kept verbatim — losing them silently would be worse than
    leaving the user with a one-off "unknown" row to investigate.
    """
    seen: set[str] = set()
    out: list[dict] = []
    for e in entries:
        if not isinstance(e, dict):
            out.append(e)
            continue
        raw_seq = e.get("sequence")
        if not isinstance(raw_seq, str):
            out.append(e)
            continue
        key = raw_seq.strip().upper()
        if not key:
            out.append(e)
            continue
        if key in seen:
            continue
        seen.add(key)
        out.append(e)
    return out


def _save_primers(entries: list[dict]) -> None:
    # Dedupe by sequence on every save (sacred policy: one entry per
    # unique sequence). See `_dedupe_primers_by_sequence`. Note that
    # NAME collisions are NOT auto-collapsed here — `.dna` imports
    # legitimately stash same-name different-sequence variants, and
    # silently dropping any would lose primer data. Name-collision
    # cleanup is opt-in via `PrimerDuplicatesModal` at startup.
    # The primer-usage cache stays valid here because adding/removing
    # primer-library entries doesn't change which plasmids carry
    # `primer_bind` features — that's `_save_library`/`_save_collections`'
    # job to invalidate.
    deduped = _dedupe_primers_by_sequence(entries)
    with _state._cache_lock:
        _safe_save_json(_state._PRIMERS_FILE, deduped, "Primer library")
        # Deep-copy on save so the cache is independent of the caller's
        # reference — pre-fix `list(entries)` shared dict refs with the
        # caller, so subsequent caller mutations (e.g. an aborted modal)
        # leaked back into the next `_load_primers`.
        _state._primers_cache = _typed_clone(deduped)
        # Mirror into the active primer collection (collections layer
        # for primers, same pattern as plasmid library →
        # collections.json mirror, pitfall #10). Inside the lock so
        # the mirror file can't drift from the live file (sweep #10,
        # [INV-50]). _save_primer_collections re-acquires the RLock.
        _sync_active_primer_collection_primers(deduped)


def _load_primer_collections() -> list[dict]:
    """Deep-copy on read so callers can mutate freely; pitfall #17.
    One-shot name-trim backfill on first load: walks every collection
    and applies `_backfill_primer_names` to its embedded `primers`
    list (parallel to `_ensure_collections_cache_populated_and_migrated`
    for the library, and to the same backfill applied to the live
    `primers.json` in `_load_primers`).
    """
    with _state._cache_lock:
        if _state._primer_collections_cache is None:
            entries, warning = _safe_load_json(
                _state._PRIMER_COLLECTIONS_FILE, "Primer collections",
            )
            if warning:
                _log.warning(warning)
            _state._primer_collections_cache = [
                e for e in entries if isinstance(e, dict)
            ]
        if not _state._primer_collections_backfill_done:
            n_changed_total = 0
            for pc in _state._primer_collections_cache:
                if not isinstance(pc, dict):
                    continue
                primers = pc.get("primers")
                if not isinstance(primers, list):
                    continue
                _, n = _backfill_primer_names(primers)
                n_changed_total += n
            _state._primer_collections_backfill_done = True
            if n_changed_total > 0:
                _log.info(
                    "primer-collections: name trim backfill cleaned "
                    "%d embedded primer(s)", n_changed_total,
                )
                _log_event(
                    "primer_collections.name_trim_backfill",
                    n_changed=n_changed_total,
                )
                save_target = _state._primer_collections_cache
            else:
                save_target = None
        else:
            save_target = None
    if save_target is not None:
        try:
            _save_primer_collections(save_target)
        except (OSError, RuntimeError):
            _log.exception(
                "primer-collections: name trim backfill save failed",
            )
    assert _state._primer_collections_cache is not None
    return _typed_clone(_state._primer_collections_cache)


def _save_primer_collections(entries: list[dict]) -> None:
    with _state._cache_lock:
        _safe_save_json(
            _state._PRIMER_COLLECTIONS_FILE, entries, "Primer collections",
        )
        _state._primer_collections_cache = _typed_clone(entries)


def _sync_active_primer_collection_primers(entries: list[dict]) -> None:
    """Mirror the live primer library into the active collection so the
    on-disk record never drifts from `primers.json`. Silent no-op if no
    collection is active or the active name has been deleted. Caller
    MUST already hold `_cache_lock` (this re-acquires via RLock for
    `_save_primer_collections`). See `[INV-50]` for the save-chain
    lock-release gap that this design avoids.
    """
    _get_name = getattr(_state, "_active_primer_collection_name_hook", None)
    name = _get_name() if _get_name is not None else None
    if not name:
        return
    snapshot = [_typed_clone(e) for e in entries if isinstance(e, dict)]
    colls = _load_primer_collections()
    for c in colls:
        if c.get("name") == name:
            c["primers"] = snapshot
            _save_primer_collections(colls)
            return


# ── Plasmid library + collections (the 160 MB SACRED path) ──────────────────
# Thin accessors only. The DANGEROUS logic stays hub-side and is reached via
# _state hooks: the `_ensure_*` migration (its backfills pull GenBank parse/
# serialise — not movable), the async active-collection mirror, and the post-
# save cache-busts. `_load_*` fire the ensure hook then return a clone of the
# `_state` cache the hub populated; `_save_*` write + reseat + fire the mirror
# (in-lock, re-enters `_save_collections` via the re-entrant RLock) + the
# after-save hook (post-lock). The data-safety chokepoint is unchanged — it
# lives in `_safe_save_json`. See `_iter_*_readonly` (hub-side) for hot reads.


def _load_library() -> list[dict]:
    # Deep-copy on read (sacred invariant #17): library entries carry nested
    # dicts (qualifiers, history, primer-pair info) — caller-side mutation would
    # otherwise poison the cache for every subsequent reader.
    hook = _state._ensure_library_hook
    if hook is not None:
        hook()
    assert _state._library_cache is not None, "library cache unpopulated (ensure hook not registered?)"
    return _typed_clone(_state._library_cache)


def _save_library(entries: list[dict], *, async_sync: bool = False) -> None:
    """Persist the live library + mirror it into the active collection.
    ``async_sync=True`` (LibraryPanel delete path) defers the mirror to a
    background worker so a 100+ MB collections.json doesn't hang the UI. The
    mirror runs INSIDE `_state._cache_lock` (sweep #10 / INV-50: the RLock lets
    the mirror's own `_save_collections` chain re-enter, and stops a concurrent
    save from drifting library.json vs collections.json). The migration + mirror
    + cache-busts stay hub-side via `_state` hooks."""
    with _state._cache_lock:
        _safe_save_json(_state._LIBRARY_FILE, entries, "Plasmid library")
        # Deep-copy into the cache so it's independent of the caller's reference.
        _state._library_cache = _typed_clone(entries)
        # Mirror into the active collection — hub-side (async worker subsystem),
        # INSIDE the lock so the mirror file can't drift from library.json.
        _mirror = getattr(_state, "_sync_active_collection_plasmids_hook", None)
        if _mirror is not None:
            _mirror(entries, async_write=async_sync)
    # Post-lock cache invalidations (primer-usage + bulk-align k-mer) — hub-side.
    _after = getattr(_state, "_after_library_save_hook", None)
    if _after is not None:
        _after()


def _load_collections() -> list[dict]:
    """Return a deepcopy of the collections list so callers can mutate entries
    (rename, edit plasmids list) without poisoning the in-memory cache (#17)."""
    hook = _state._ensure_collections_hook
    if hook is not None:
        hook()
    assert _state._collections_cache is not None, "collections cache unpopulated (ensure hook not registered?)"
    return _typed_clone(_state._collections_cache)


def _save_collections(entries: list[dict]) -> None:
    with _state._cache_lock:
        _safe_save_json(_state._COLLECTIONS_FILE, entries, "Plasmid collections")
        _state._collections_cache = _typed_clone(entries)
    # Post-save: invalidate the primer-usage cache (the plasmid set changed) —
    # hub-side via the after-save hook.
    _after = getattr(_state, "_after_collections_save_hook", None)
    if _after is not None:
        _after()


# ── Model collections (BABS model picker — Ollama model refs, see [INV-139]) ──
# Mirrors the library/collections architecture but MUCH lighter: a collection
# groups model REFERENCES (Ollama tags / `hf.co/<repo>` pull strings) under a
# name, never embedded model data (the models live in Ollama's own store). A
# collection dict is {name, description, models:[{ref, note, added}], saved};
# install status is resolved LIVE at display time, so a ref may sit in a
# collection before (or after) the model is pulled. Same safe-save chokepoint
# (#7) + deep-copy-on-read/save cache discipline (#17) as every accessor here.
# No migration hook (schema v1, new file) — self-populating, double-checked.
def _load_model_collections() -> "list[dict]":
    """Deep-copy on read so callers can freely mutate the collections + their
    model lists (pitfall #17). Lock-free cache-hit fast path."""
    cached = _state._model_collections_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._model_collections_cache is None:
            entries, warning = _safe_load_json(_state._MODEL_COLLECTIONS_FILE, "Model collections")
            if warning:
                _log.warning(warning)
            _state._model_collections_cache = _typed_clone(
                [e for e in entries if isinstance(e, dict)])
        return _typed_clone(_state._model_collections_cache)


def _save_model_collections(entries: "list[dict]") -> None:
    """Persist the model-collection list (atomic write + .bak via the
    chokepoint, invariant #7) and re-seat the cache via `_typed_clone`
    (invariant #17 — a caller mutating its list post-save can't leak edits
    into the next reader)."""
    with _state._cache_lock:
        _safe_save_json(_state._MODEL_COLLECTIONS_FILE, entries, "Model collections")
        _state._model_collections_cache = _typed_clone(entries)


def _iter_model_collections_readonly() -> "list[dict]":
    """Read-only view of the model-collection cache (NO clone). Callers MUST
    NOT mutate any returned dict/list — that poisons the cache (#17). For hot
    reads (name-collision checks, listing). Mutators use `_load_model_collections`."""
    cached = _state._model_collections_cache
    if cached is not None:
        return cached
    with _state._cache_lock:
        if _state._model_collections_cache is None:
            entries, warning = _safe_load_json(_state._MODEL_COLLECTIONS_FILE, "Model collections")
            if warning:
                _log.warning(warning)
            _state._model_collections_cache = _typed_clone(
                [e for e in entries if isinstance(e, dict)])
        return _state._model_collections_cache


def _find_model_collection(name: str) -> "dict | None":
    """O(N) walk + clone of the collection matching `name` (trimmed, exact),
    or None. Clone so the caller can mutate safely (#17)."""
    name = (name or "").strip()
    if not name:
        return None
    for c in _iter_model_collections_readonly():
        if isinstance(c, dict) and (c.get("name") or "").strip() == name:
            return _typed_clone(c)
    return None


def _model_collection_name_taken(name: str) -> bool:
    """Case-insensitive name-collision check over the readonly view."""
    key = (name or "").strip().lower()
    if not key:
        return False
    return any(isinstance(c, dict) and (c.get("name") or "").strip().lower() == key
               for c in _iter_model_collections_readonly())


# ── Parts bin / parts-bin collections (multi-bin part storage) ───────────────
# Same architecture as library/collections: `parts_bin.json` holds the ACTIVE
# bin's parts; `parts_bin_collections.json` is the canonical multi-bin record.
# These 4 accessors are THIN — the `_ensure_*` migration + sequence backfill,
# the active-parts-bin mirror (`_sync_active_parts_bin_parts`, synchronous since
# bins are small), and the assembly-fragment-cache bust all STAY hub-side and
# are reached via `_state` hooks (registered by the hub after every target is
# defined). `_iter_parts_bin_readonly`, `_find_parts_bin`, the active-name
# getter/setter, and `_switch_active_parts_bin` stay hub-side.
def _load_parts_bin() -> list[dict]:
    # Deep-copy on read (sacred invariant #17): parts entries carry nested dicts
    # (qualifiers, primer pairs, mutation lists) — caller-side mutation would
    # otherwise poison the cache for every subsequent reader.
    hook = _state._ensure_parts_bin_hook
    if hook is not None:
        hook()
    assert _state._parts_bin_cache is not None, "parts-bin cache unpopulated (ensure hook not registered?)"
    return _typed_clone(_state._parts_bin_cache)


def _save_parts_bin(entries: list[dict]) -> None:
    """Persist the live parts bin + mirror it into the active parts-bin
    collection. The mirror runs INSIDE `_state._cache_lock` (sweep #10: the RLock
    lets the mirror's own `_save_parts_bin_collections` chain re-enter, and stops
    a concurrent save from drifting parts_bin.json vs parts_bin_collections.json).
    The migration + mirror + assembly-fragment-cache bust stay hub-side via the
    `_state` hooks."""
    with _state._cache_lock:
        _safe_save_json(_state._PARTS_BIN_FILE, entries, "Parts bin")
        # Deep-copy into the cache so it's independent of the caller's reference.
        _state._parts_bin_cache = _typed_clone(entries)
        # Mirror into the active parts-bin collection — hub-side, INSIDE the lock
        # so the mirror file can't drift from parts_bin.json.
        _mirror = getattr(_state, "_sync_active_parts_bin_parts_hook", None)
        if _mirror is not None:
            _mirror(entries)
    # Post-lock: a part edit shifts what `_assembly_fragment_from_source` yields,
    # so invalidate the assembly-fragment digest cache — hub-side.
    _after = getattr(_state, "_after_parts_bin_save_hook", None)
    if _after is not None:
        _after()


def _load_parts_bin_collections() -> list[dict]:
    """Return a deepcopy of the parts-bin collections list so callers can mutate
    entries (rename, edit parts list) without poisoning the in-memory cache (#17)."""
    hook = _state._ensure_parts_bin_collections_hook
    if hook is not None:
        hook()
    assert _state._parts_bin_collections_cache is not None, "parts-bin collections cache unpopulated (ensure hook not registered?)"
    return _typed_clone(_state._parts_bin_collections_cache)


def _save_parts_bin_collections(entries: list[dict]) -> None:
    """Persist the full parts-bin collections list. Deepcopies into the cache so
    caller mutations after save can't poison subsequent loaders (#17)."""
    with _state._cache_lock:
        _safe_save_json(_state._PARTS_BIN_COLLECTIONS_FILE, entries, "Parts-bin collections")
        _state._parts_bin_collections_cache = _typed_clone(entries)
    # Post-save: a bin add/delete/rename/parts-edit shifts the active bin's parts
    # identity → invalidate the assembly-fragment digest cache — hub-side.
    _after = getattr(_state, "_after_parts_bin_collections_save_hook", None)
    if _after is not None:
        _after()


# ── Feature library ──────────────────────────────────────────────────────────
# Clean cluster (no migration / mirror / cache-bust): just the cache + the
# `_state._features_generation` counter consumers watch to rebuild derived
# indices. The upsert / scan helpers stay hub-side.
def _load_features() -> list[dict]:
    """Return an independent (deep-copied) list of feature library
    entries. Callers can mutate the returned dicts freely without
    poisoning the cache — important for ``FeatureLibraryScreen``
    which buffers in-place edits (rename / color / strand / etc.) and
    then either persists or abandons. A shallow ``list(_state._features_cache)``
    used to share dict refs with the cache, so an abandoned mutation
    would survive in the cache and leak into the next ``_load_features``
    consumer (a freshly opened FeatureLibraryScreen, the
    DomesticatorModal feature picker, etc.) as if it had been saved.
    """
    # Sweep #26: double-checked locking — cache-hit fast path is lock-free.
    cached = _state._features_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._features_cache is None:
            entries, warning = _safe_load_json(_state._FEATURES_FILE, "Feature library")
            if warning:
                _log.warning(warning)
            entries = [e for e in entries if isinstance(e, dict)]
            _state._features_cache = entries
            # A fresh disk read is the result of either first-load or an
            # external invalidation (test harness setting `_state._features_cache =
            # None`, or a hand-edit of features.json). Either way the contents
            # may have changed since the last write, so bump the generation so
            # consumers know to rebuild any derived indices.
            _state._features_generation += 1
        return _typed_clone(_state._features_cache)


def _save_features(entries: list[dict]) -> None:
    """Persist `entries` and seed the in-memory cache with a deepcopy
    so subsequent caller-side mutations of `entries` (or any dict
    inside it) cannot leak into the cache after the save returns.
    Without the deepcopy, the dicts in `_state._features_cache` would alias
    the dicts in the caller's list — so e.g. a FeatureLibraryScreen
    that saved, then made another change, then abandoned, would leave
    the post-save mutations stuck in the cache.
    """
    with _state._cache_lock:
        _safe_save_json(_state._FEATURES_FILE, entries, "Feature library")
        _state._features_cache = _typed_clone(entries)
        _state._features_generation += 1


# ── Protein collections ──────────────────────────────────────────────────────
# Clean cluster: `_ensure_*` is a pure disk-load+filter (no GenBank backfill /
# migration), its only caller is `_load_protein_collections`, so it moves
# wholesale — no hook. Save is a plain write+reseat. The add/scan helpers stay
# hub-side.
def _ensure_protein_collections_cache() -> None:
    if _state._protein_collections_cache is not None:
        return
    with _state._cache_lock:
        if _state._protein_collections_cache is None:
            entries, warning = _safe_load_json(
                _state._PROTEIN_COLLECTIONS_FILE, "Protein collections")
            if warning:
                _log.warning(warning)
            _state._protein_collections_cache = [
                e for e in (entries or []) if isinstance(e, dict)]


def _load_protein_collections() -> list[dict]:
    """Return a deepcopy of the protein-collections list so callers can
    mutate freely without poisoning the in-memory cache ([PIT-17])."""
    _ensure_protein_collections_cache()
    assert _state._protein_collections_cache is not None
    return _typed_clone(_state._protein_collections_cache)


def _save_protein_collections(entries: list[dict]) -> None:
    """Persist the full protein-collections list (envelope via
    `_safe_save_json` — `.bak` + atomic replace, re-raises on failure)
    and reseat the cache, under `_cache_lock`."""
    with _state._cache_lock:
        _safe_save_json(_state._PROTEIN_COLLECTIONS_FILE, entries, "Protein collections")
        _state._protein_collections_cache = _typed_clone(entries)


# ── Experiments notebook ─────────────────────────────────────────────────────
# Mirror cluster (like library/collections, synchronous mirror). Load applies
# the legacy tag-format migration per body via `_state._migrate_experiment_body_hook`
# (the migrator stays hub-side — it's shared with the editor body-readers). Save
# mirrors into the active experiment project via `_sync_active_project_experiments`
# (hub-side, fired INSIDE the lock per sweep #9). The editor / attachment helpers
# stay hub-side.
def _load_experiments() -> "list[dict]":
    """Load the experiments-notebook entries. Cached + deepcopy-on-read
    per sacred invariant #17 so a caller mutating returned dicts can't
    poison the cache for the next reader. Filters non-dict entries
    defensively (hand-edited JSON / schema drift). Applies the
    legacy-tag-format migration (hub-side via the body hook) to every
    entry's body so the editor only ever sees the new single-sigil tokens."""
    # Double-checked locking (match `_load_features`): the cache-hit fast path
    # is lock-free; the populate holds `_cache_lock` so it can't race a
    # concurrent `_save_experiments` (agent thread) and clobber its cache reseat
    # (lost-update → the UI shows stale entries / the agent's write is dropped).
    cached = _state._experiments_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._experiments_cache is None:
            entries, warning = _safe_load_json(_state._EXPERIMENTS_FILE, "Experiments")
            if warning:
                _log.warning(warning)
            _migrate = getattr(_state, "_migrate_experiment_body_hook", None)
            cleaned: "list[dict]" = []
            n_migrated = 0
            for e in entries:
                if not isinstance(e, dict):
                    continue
                body = e.get("body_md")
                if isinstance(body, str) and _migrate is not None:
                    migrated = _migrate(body)
                    if migrated != body:
                        e["body_md"] = migrated
                        n_migrated += 1
                cleaned.append(e)
            if n_migrated:
                _log_event(
                    "experiments.tag.migrated",
                    n_entries=n_migrated, n_loaded=len(cleaned),
                )
            _state._experiments_cache = cleaned
        return _typed_clone(_state._experiments_cache)


def _save_experiments(entries: "list[dict]") -> None:
    """Persist the experiments-notebook entries through the full four-
    layer data-safety net (invariant #31). Deep-clones into the cache
    on save (invariant #17), under `_state._cache_lock` so concurrent saves
    don't desync disk/cache ordering (invariant #41 — concurrency).

    After the primary save, mirrors the entries into the active
    experiment project's `experiments` list — sacred contract: every
    writeable experiments path goes through this helper. The mirror is
    hub-side (`_sync_active_project_experiments`), fired via the hook
    INSIDE the lock (sweep #9: RLock re-entry; a concurrent reader can't
    observe experiments.json updated while experiment_projects.json lags)."""
    with _state._cache_lock:
        _safe_save_json(_state._EXPERIMENTS_FILE, entries, "Experiments")
        _state._experiments_cache = _typed_clone(entries)
        # Mirror inside the lock — hub-side via the hook.
        _mirror = getattr(_state, "_sync_active_project_experiments_hook", None)
        if _mirror is not None:
            _mirror(entries)


# ── Settings ─────────────────────────────────────────────────────────────────
# The hot path (`_get_setting`/`_set_setting` fire on every keystroke / render).
# Only these 4 accessors move; the type-validation web (`_validate_settings` +
# `_SETTINGS_SCHEMA` + safe-identifier check) and the coalesced background disk-
# flush subsystem (`_settings_flush_worker` daemon + `_settings_flush_sync` +
# UI-failure notify) stay hub-side, reached via `_state._validate_settings_hook`
# / `_state._settings_schedule_flush_hook`.
def _load_settings() -> dict:
    """Return the persistent settings dict. Stored on disk as a list of
    ``{"key": ..., "value": ...}`` envelope entries so it shares the
    schema layout (sacred invariant #7) with every other JSON file.

    Type-validates against `_SETTINGS_SCHEMA` (hub-side, via the validate
    hook) so a hand-edited settings.json (or a partial restore from a
    future-version snapshot) can't propagate wrong-type values into the
    reactive UI. Unknown keys are preserved for forward-compat.
    """
    # Sweep #26: double-checked locking — the cache-hit fast path is lock-free.
    cached = _state._settings_cache
    if cached is not None:
        return _typed_clone(cached)
    with _state._cache_lock:
        if _state._settings_cache is None:
            entries, warning = _safe_load_json(_state._SETTINGS_FILE, "Settings")
            if warning:
                _log.warning(warning)
            settings: dict = {}
            for e in entries:
                if not isinstance(e, dict):
                    continue
                k, v = e.get("key"), e.get("value")
                if isinstance(k, str):
                    settings[k] = v
            _validate = getattr(_state, "_validate_settings_hook", None)
            if _validate is not None:
                cleaned, warns = _validate(settings)
                for w in warns:
                    _log.warning(w)
            else:
                cleaned = settings
            _state._settings_cache = cleaned
        return _typed_clone(_state._settings_cache)


def _save_settings(settings: dict) -> None:
    """Synchronous write — kept for the few callsites (tests, migrations)
    that want the disk state visible immediately on return. UI toggles
    use `_set_setting`, which mirrors to `_state._settings_cache` synchronously
    and defers the disk write to a background thread."""
    entries = [{"key": k, "value": v} for k, v in settings.items()]
    with _state._cache_lock:
        _safe_save_json(_state._SETTINGS_FILE, entries, "Settings")
        _state._settings_cache = _typed_clone(settings)


def _get_setting(key: str, default: "_Any" = None) -> "_Any":
    """Return the persisted value for ``key``, or ``default`` if absent.
    Typed as ``Any`` so callers like ``_active_grammar_id`` can declare
    a narrower return type (`-> str`) without a cast — pyright treats
    ``Any`` as universally assignable, sidestepping the
    ``Unknown | None`` propagation from the loose `_load_settings()`
    return.

    Sweep #25 (2026-05-23): reads the cache directly instead of paying
    a full ``_load_settings()`` deepcopy on every call. With 59 call-
    sites — including 4 per ``EditCustomizeModal`` open and 24 per
    ``_get_active_collection_name()`` — the per-call full-dict clone
    burned ~3–5 ms across the app per second of typical use. Settings
    values are JSON primitives (str/int/float/bool/None) for ~99% of
    keys; ``_typed_clone`` of an immutable returns the same object,
    so the cost is negligible. The handful of container values
    (e.g. ``experiments_custom_dict``) still clone defensively so a
    caller mutation can't poison the cache.
    """
    if _state._settings_cache is None:
        # Trigger cache populate (one-time per session).
        _load_settings()
    with _state._cache_lock:
        cache = _state._settings_cache
        if cache is None:
            return default
        if key not in cache:
            return default
        return _typed_clone(cache[key])


# Setting keys whose VALUE is a secret — redacted in the change-log + event
# stream by `_set_setting`, and deliberately kept out of the agent settings
# allowlist so a remote agent can't read them back.
_SENSITIVE_SETTING_KEYS: "frozenset[str]" = frozenset({
    "plasmidsaurus_client_secret",
    # Optional online-lookup provider keys (web + patent search). Secret, and
    # kept out of `_AGENT_SETTINGS_ALLOWLIST` so an agent can't read them back.
    "brave_search_api_key",
    "patentsview_api_key",
})


def _set_setting(key: str, value) -> None:
    """Update one setting key. Cache is updated synchronously so a
    subsequent `_load_settings()` (in this process) sees the new value
    immediately. The disk write is scheduled on a daemon thread (hub-side
    via the schedule-flush hook) so a keypress that toggles a setting
    doesn't block the UI on fsync; bursts of toggles coalesce into fewer
    disk writes.

    Set ``SPLICECRAFT_SKIP_SETTINGS_FLUSH=1`` to bypass the daemon
    thread and write synchronously — used by tests so they get
    deterministic disk state without a trailing daemon thread on exit.
    """
    # Sweep #26 (2026-05-25): hold `_state._cache_lock` across the RMW so
    # two concurrent `_set_setting(...)` calls from agent threads can't both
    # load the same pre-mutation cache, each mutate a local copy with a
    # different key, and the second reseat the global cache — silently
    # losing the first writer's key.
    with _state._cache_lock:
        settings = _load_settings()
        # Log only when the value actually changes — avoids spamming the log
        # when callers re-set the same value on settings hydration. Values
        # run through `_repr_for_log` to truncate over-long lists / cap any
        # accidental sequence content.
        prev = settings.get(key, "<unset>")
        if prev != value:
            # Never write a secret-valued setting (e.g. the Plasmidsaurus
            # client secret) to the log or the event stream — redact both the
            # new value AND the previous one to a fixed marker.
            if key in _SENSITIVE_SETTING_KEYS:
                v_log = p_log = "<redacted>"
            else:
                v_log = _repr_for_log(value)
                p_log = _repr_for_log(prev)
            _log.info("setting changed: %s = %r (was %r)", key, v_log, p_log)
            _log_event(
                "settings.changed", key=key, value=v_log, prev=p_log,
            )
        settings[key] = value
        # Sweep #9: `_typed_clone` (deepcopy), not shallow `dict(...)` — a
        # nested-list value would otherwise share the caller's reference and a
        # later caller-side mutation would leak into the cache (invariant #17).
        _state._settings_cache = _typed_clone(settings)
    # Schedule the coalesced background disk flush — hub-side via the hook
    # (keeps the daemon worker + SKIP_FLUSH env handling + UI-failure notify
    # off the data layer).
    _schedule = getattr(_state, "_settings_schedule_flush_hook", None)
    if _schedule is not None:
        _schedule(settings)


# ── Library / collections read-only views + finders + active-collection pointer
# Hot-path accessors moved from the hub (Phase D) so the modal/screen siblings
# that call them bare can import them. The readonly iters return the cache list
# WITHOUT cloning (callers MUST NOT mutate — #17); finders clone only the matched
# entry. The `_ensure_*` migration stays hub-side, reached via the ensure hooks.
def _iter_library_readonly() -> list[dict]:
    """Read-only view of the library cache — returns the cached list
    directly without cloning. Callers MUST NOT mutate any returned
    dict or its nested lists; doing so poisons the in-memory cache
    (pitfall #17). Use for hot read paths (listing endpoints, search
    filters, primer-usage scans, name-collision precheck) where the
    multi-second `_typed_clone` cost on a 100+ MB library is
    meaningful. For any path that mutates entries, use
    ``_load_library()`` instead. Mirrors ``_iter_collections_readonly``
    (sweep #23) and added as part of sweep #25 (2026-05-23) to close
    the perf gap on ~40 library-read callsites.
    """
    _ensure = _state._ensure_library_hook
    if _ensure is not None:
        _ensure()
    assert _state._library_cache is not None
    return _state._library_cache


def _find_library_entry_by_id(entry_id: "str | None") -> "dict | None":
    """Return a deep-clone of the library entry whose ``id`` matches,
    or ``None`` if no match. Accepts ``None`` / non-string / empty
    string and returns ``None`` for those (callsite-friendly).

    Sweep #11 (2026-05-20): pre-fix ~80 callsites did
    ``next((e for e in _load_library() if e.get("id") == X), None)``,
    which clones the ENTIRE library on every lookup. On a 100+ MB
    library this is multi-second per click for status changes,
    history-viewer opens, picker dismisses. This helper takes
    ``_cache_lock``, walks ``_state._library_cache`` directly, and clones
    only the matched entry — O(N) walk with O(1) cloned bytes
    instead of O(N) cloned bytes.

    Returns a deep copy of the entry so callers can mutate freely
    without poisoning the cache (same contract as ``_load_library``
    per invariant #17). Returns ``None`` for empty/non-string
    ``entry_id``; the empty-string skip prevents accidental
    aliasing against id-less library rows (defensive against
    hand-edited JSON).
    """
    if not isinstance(entry_id, str) or not entry_id:
        return None
    _ensure = _state._ensure_library_hook
    if _ensure is not None:
        _ensure()
    assert _state._library_cache is not None
    with _state._cache_lock:
        for e in _state._library_cache:
            if e.get("id") == entry_id:
                return _typed_clone(e)
    return None


def _find_library_entry_by_name(entry_name: str) -> "dict | None":
    """Return a deep-clone of the first library entry whose ``name``
    matches, or ``None`` if no match. Sibling of
    ``_find_library_entry_by_id`` for the rare endpoints (currently
    ``_h_align_plasmidsaurus_zip``) that accept either an id or a
    name. Added by sweep #25 (2026-05-23).

    Note: ``name`` is NOT guaranteed unique in the library
    (collision modals can save a deliberate COPY); returns the first
    hit, walking in insertion order. Callers needing the canonical
    one should prefer ``_find_library_entry_by_id``.
    """
    if not isinstance(entry_name, str) or not entry_name:
        return None
    _ensure = _state._ensure_library_hook
    if _ensure is not None:
        _ensure()
    assert _state._library_cache is not None
    with _state._cache_lock:
        for e in _state._library_cache:
            if e.get("name") == entry_name:
                return _typed_clone(e)
    return None


def _iter_collections_readonly() -> list[dict]:
    """Read-only view of the collections cache — returns the cached list
    directly without cloning. Callers MUST NOT mutate any returned dict
    or its nested lists; doing so poisons the in-memory cache (pitfall
    #17). Use for hot read paths (search modals, listing endpoints) where
    the ~30–50 ms `_typed_clone` cost is meaningful on big libraries.
    For any path that mutates entries, use `_load_collections()` instead.
    """
    _ensure = _state._ensure_collections_hook
    if _ensure is not None:
        _ensure()
    assert _state._collections_cache is not None
    return _state._collections_cache


def _get_active_collection_name() -> "str | None":
    val = _get_setting("active_collection", None)
    return val if isinstance(val, str) and val else None


def _set_active_collection_name(name: "str | None") -> None:
    """Persist (or clear) the active-collection pointer."""
    prev = _get_active_collection_name()
    target = name or ""
    if prev != target:
        _log_event("collection.switched", prev=prev or "", new=target)
    _set_setting("active_collection", target)


def _find_collection(name: str) -> "dict | None":
    """Return a deep-clone of the collection whose ``name`` matches, or
    ``None`` if not found. Mirrors ``_find_library_entry_by_id`` —
    walks the readonly cache view and only ``_typed_clone``s the
    matched entry instead of the entire collections list.

    Sweep #26 (2026-05-25): pre-fix did ``for c in _load_collections()``
    which deep-clones the entire ~160 MB collections.json file every
    call. 13 call sites — each pays 400 ms – 1.6 s. Agent
    ``bulk-import-folder`` paid 2× (early-out + atomic re-check).
    """
    if not isinstance(name, str) or not name:
        return None
    with _state._cache_lock:
        for c in _iter_collections_readonly():
            if isinstance(c, dict) and c.get("name") == name:
                return _typed_clone(c)
    return None


def _collection_name_taken(name: str) -> bool:
    """Dup-name guard for create / rename. Pure check, no side effects.
    Sweep #26: routed through the readonly walk so the check is O(N)
    in memory rather than O(N) in deep-clone bytes."""
    if not isinstance(name, str) or not name:
        return False
    with _state._cache_lock:
        for c in _iter_collections_readonly():
            if isinstance(c, dict) and c.get("name") == name:
                return True
    return False


# ── Golden Braid / MoClo grammar + enzyme data + derived accessors (Phase D) ─
# Pure GB grammar/enzyme definition constants + `_all_grammars` (merges builtin
# with `_load_custom_grammars`) + `_get_entry_vector` (walks `_load_entry_vectors`).
# Heavily referenced across the constructor/domestication hub code (re-exported).
_NEB_ENZYMES: dict[str, tuple[str, int, int]] = {

    # ── Common Type IIP — 6-bp palindromic cutters ─────────────────────────────
    "EcoRI":     ("GAATTC",       1,  5),  # G^AATTC   / CTTAA^G     5' overhang
    "EcoRV":     ("GATATC",       3,  3),  # GAT^ATC   / GAT^ATC     blunt
    "BamHI":     ("GGATCC",       1,  5),  # G^GATCC   / CCTAG^G     5' overhang
    "HindIII":   ("AAGCTT",       1,  5),  # A^AGCTT   / TTCGA^A     5' overhang
    "NcoI":      ("CCATGG",       1,  5),  # C^CATGG   / GGTAC^C     5' overhang
    "NdeI":      ("CATATG",       2,  4),  # CA^TATG   / GTAT^AC     5' overhang
    "XhoI":      ("CTCGAG",       1,  5),  # C^TCGAG   / GAGCT^C     5' overhang
    "SalI":      ("GTCGAC",       1,  5),  # G^TCGAC   / CAGCT^G     5' overhang
    "KpnI":      ("GGTACC",       5,  1),  # GGTAC^C   / G^GTACC     3' overhang
    # Acc65I is KpnI's NEOschizomer, not a synonym: same GGTACC, cut on the
    # OTHER side of the tetramer, so it leaves a 4-nt 5' overhang where KpnI
    # leaves a 3' one. It gets a real catalog entry (an alias to KpnI would
    # hand back the wrong overhang) — Asp718I is its true synonym and lives in
    # `_ENZYME_ALIASES`.
    "Acc65I":    ("GGTACC",       1,  5),  # G^GTACC   / CCATG^G     5' overhang
    "SacI":      ("GAGCTC",       5,  1),  # GAGCT^C   / G^AGCTC     3' overhang
    "SacII":     ("CCGCGG",       4,  2),  # CCGC^GG   / CC^GCGG     3' overhang
    "SpeI":      ("ACTAGT",       1,  5),  # A^CTAGT   / TGATC^A     5' overhang
    "XbaI":      ("TCTAGA",       1,  5),  # T^CTAGA   / AGATC^T     5' overhang
    "NotI":      ("GCGGCCGC",     2,  6),  # GC^GGCCGC / CGCCGG^CG   5' overhang (8-cutter)
    "PstI":      ("CTGCAG",       5,  1),  # CTGCA^G   / G^CTGCA     3' overhang
    "SphI":      ("GCATGC",       5,  1),  # GCATG^C   / C^GCATG     3' overhang
    "ClaI":      ("ATCGAT",       2,  4),  # AT^CGAT   / CGAT^AT     5' overhang
    "NheI":      ("GCTAGC",       1,  5),  # G^CTAGC   / CGATC^G     5' overhang
    "AvaI":      ("CYCGRG",       1,  5),  # C^YCGRG                 5' overhang (degenerate)
    "AvaII":     ("GGWCC",        1,  4),  # G^GWCC    / CCWG^G      5' overhang
    "AvrII":     ("CCTAGG",       1,  5),  # C^CTAGG   / GGATC^C     5' overhang
    "BclI":      ("TGATCA",       1,  5),  # T^GATCA   / AGTCA^T     5' overhang (dam-sensitive)
    "BglII":     ("AGATCT",       1,  5),  # A^GATCT   / TCTAG^A     5' overhang
    "BsiWI":     ("CGTACG",       1,  5),  # C^GTACG   / GCATG^C     5' overhang
    "BspEI":     ("TCCGGA",       1,  5),  # T^CCGGA   / AGGCC^T     5' overhang
    "BsrGI":     ("TGTACA",       1,  5),  # T^GTACA   / ACATG^T     5' overhang
    "BssHII":    ("GCGCGC",       1,  5),  # G^CGCGC   / CGCGC^G     5' overhang
    "BstBI":     ("TTCGAA",       2,  4),  # TT^CGAA   / AAGC^TT     5' overhang
    "BstEII":    ("GGTNACC",      1,  6),  # G^GTNACC / CCANTG^G  5-nt 5' overhang GTNAC
    "BstXI":     ("CCANNNNNNTGG", 8,  4),  # CCANNNNN^NTGG  4-nt 3' overhang (12-bp recog)
    "BstYI":     ("RGATCY",       1,  5),  # R^GATCY   / YCTAG^R     5' overhang
    "CpoI":      ("CGGWCCG",      2,  5),  # CG^GWCCG  / CGCC^WGG    5' overhang
    "DraI":      ("TTTAAA",       3,  3),  # TTT^AAA   / TTT^AAA     blunt
    "DraIII":    ("CACNNNGTG",    6,  3),  # CACNNN^GTG/ GTG^NNNGTG  3' overhang
    "EagI":      ("CGGCCG",       1,  5),  # C^GGCCG   / GCCGG^C     5' overhang (NotI subset)
    "Eco47III":  ("AGCGCT",       3,  3),  # AGC^GCT                 blunt
    "Eco53kI":   ("GAGCTC",       3,  3),  # GAG^CTC                 blunt (SacI neoschizomer)
    "EcoNI":     ("CCTNNNNNAGG",  5,  6),  # CCTNN^NNNAGG            5' overhang
    "FseI":      ("GGCCGGCC",     6,  2),  # GGCCGG^CC / CC^GGCCGG   3' overhang (8-cutter)
    "FspI":      ("TGCGCA",       3,  3),  # TGC^GCA                 blunt
    "HaeII":     ("RGCGCY",       5,  1),  # RGCGC^Y   / R^GCGCY     3' overhang
    "HaeIII":    ("GGCC",         2,  2),  # GG^CC                   blunt (4-cutter)
    "HincII":    ("GTYRAC",       3,  3),  # GTY^RAC                 blunt
    "HindII":    ("GTYRAC",       3,  3),  # GTY^RAC                 blunt (HincII isoschizomer)
    "HpaI":      ("GTTAAC",       3,  3),  # GTT^AAC                 blunt
    "HpaII":     ("CCGG",         1,  3),  # C^CGG     / CGG^C        5' overhang (4-cutter)
    "MfeI":      ("CAATTG",       1,  5),  # C^AATTG   / GTTAA^C     EcoRI-compatible ends
    "MluI":      ("ACGCGT",       1,  5),  # A^CGCGT   / TGCGC^A     5' overhang
    "MscI":      ("TGGCCA",       3,  3),  # TGG^CCA                 blunt
    "MspI":      ("CCGG",         1,  3),  # C^CGG     / CGG^C        5' overhang (HpaII isoschizomer)
    "MunI":      ("CAATTG",       1,  5),  # C^AATTG                 MfeI isoschizomer
    "NarI":      ("GGCGCC",       2,  4),  # GG^CGCC   / CGCG^CC     3' overhang
    "NruI":      ("TCGCGA",       3,  3),  # TCG^CGA                 blunt
    "NsiI":      ("ATGCAT",       5,  1),  # ATGCA^T   / T^ATGCA     PstI-compatible ends
    "NspI":      ("RCATGY",       5,  1),  # RCATG^Y   / R^CATGY     3' overhang
    "PacI":      ("TTAATTAA",     5,  3),  # TTAAT^TAA / TTA^ATTAA   3' overhang (8-cutter)
    "PaeR7I":    ("CTCGAG",       1,  5),  # C^TCGAG                 XhoI isoschizomer
    "PciI":      ("ACATGT",       1,  5),  # A^CATGT   / TGTAC^A     5' overhang
    "PmeI":      ("GTTTAAAC",     4,  4),  # GTTT^AAAC               blunt (8-cutter)
    "PmlI":      ("CACGTG",       3,  3),  # CAC^GTG                 blunt
    "PscI":      ("ACATGT",       1,  5),  # A^CATGT                 PciI isoschizomer
    "PvuI":      ("CGATCG",       4,  2),  # CGATC^G   / G^CGATC     3' overhang
    "PvuII":     ("CAGCTG",       3,  3),  # CAG^CTG                 blunt
    "RsrII":     ("CGGWCCG",      2,  5),  # CG^GWCCG                CpoI isoschizomer
    "SbfI":      ("CCTGCAGG",     6,  2),  # CCTGCA^GG / CC^TGCAGG   PstI-compatible (8-cutter)
    "ScaI":      ("AGTACT",       3,  3),  # AGT^ACT                 blunt
    "SfiI":      ("GGCCNNNNNGGCC",8,  5),  # GGCCN^NNN^NGGCC         3-nt 3' overhang (13-bp)
    "SgrAI":     ("CRCCGGYG",     2,  6),  # CR^CCGGYG / GCCGGR^C    5' overhang
    "SmaI":      ("CCCGGG",       3,  3),  # CCC^GGG                 blunt
    "SnaBI":     ("TACGTA",       3,  3),  # TAC^GTA                 blunt
    "SrfI":      ("GCCCGGGC",     4,  4),  # GCCC^GGGC               blunt (8-cutter)
    "StuI":      ("AGGCCT",       3,  3),  # AGG^CCT                 blunt
    "SwaI":      ("ATTTAAAT",     4,  4),  # ATTT^AAAT               blunt (8-cutter)
    "Tth111I":   ("GACNNNGTC",    4,  5),  # GACN^NNGTC              1-base 3' overhang
    "XmaI":      ("CCCGGG",       1,  5),  # C^CCGGG   / GGGCC^C     5' overhang (SmaI isoschizomer)
    "XmnI":      ("GAANNNNTTC",   5,  5),  # GAANN^NNTTC             blunt

    # ── Rare 8-cutters ─────────────────────────────────────────────────────────
    "AscI":      ("GGCGCGCC",     2,  6),  # GG^CGCGCC / CGCGCC^GG   5' overhang
    "AsiSI":     ("GCGATCGC",     5,  3),  # GCGAT^CGC / GCG^ATCGC   3' overhang

    # ── Degenerate / IUPAC recognition sequences ───────────────────────────────
    "AccI":      ("GTMKAC",       2,  4),  # GT^MKAC / CAKM^TG       2-nt 5' overhang
    "AclI":      ("AACGTT",       2,  4),  # AA^CGTT / TTGC^AA       2-nt 5' overhang CG
    "AfeI":      ("AGCGCT",       3,  3),  # AGC^GCT                 blunt (Eco47III isoschizomer)
    "AflII":     ("CTTAAG",       1,  5),  # C^TTAAG                 MfeI-compatible ends
    "AflIII":    ("ACRYGT",       1,  5),  # A^CRYGT                 MluI-compatible ends
    "AgeI":      ("ACCGGT",       1,  5),  # A^CCGGT   / TGGCC^A     5' overhang
    "AhdI":      ("GACNNNNNGTC",  6,  5),  # GACNNNN^NGTC            1-base 3' overhang
    "AluI":      ("AGCT",         2,  2),  # AG^CT                   blunt (4-cutter)
    "ApaI":      ("GGGCCC",       5,  1),  # GGGCC^C   / G^GGCCC     3' overhang
    "ApaLI":     ("GTGCAC",       1,  5),  # G^TGCAC                 SphI-compatible ends
    "ApoI":      ("RAATTY",       1,  5),  # R^AATTY                 EcoRI isoschizomer (degenerate)
    "AatII":     ("GACGTC",       5,  1),  # GACGT^C   / G^ACGTC     3' overhang
    "BaeGI":     ("GKGCMC",       5,  1),  # GKGCM^C   / G^KGCMC     3' overhang
    "BglI":      ("GCCNNNNNGGC",  7,  4),  # GCCNNNN^NGGC            3' overhang
    "BmgBI":     ("CACGTC",       3,  3),  # CAC^GTC                 blunt
    "BsaAI":     ("YACGTR",       3,  3),  # YAC^GTR                 blunt
    "BsaBI":     ("GATNNNNATC",   5,  5),  # GATN4^ATC               blunt
    "BsaHI":     ("GRCGYC",       2,  4),  # GR^CGYC                 3' overhang
    "BsaWI":     ("WCCGGW",       1,  5),  # W^CCGGW                 5' overhang
    "BseYI":     ("CCCAGC",       1,  5),  # C^CCAGC   / GCTGG^G     5' overhang
    "BsiEI":     ("CGRYCG",       4,  2),  # CGRY^CG                 3' overhang
    "BsiHKAI":   ("GWGCWC",       5,  1),  # GWGCW^C                 3' overhang
    "BsrFI":     ("RCCGGY",       1,  5),  # R^CCGGY                 5' overhang
    "Bsp1286I":  ("GDGCHC",       5,  1),  # GDGCH^C                 3' overhang
    "BspHI":     ("TCATGA",       1,  5),  # T^CATGA                 NcoI-compatible ends
    "BsrI":      ("ACTGG",        6,  4),  # ACTGGN^/^NCCAGT 2-nt 3' overhang Type IIS
    "BstAPI":    ("GCANNNNNTGC",  7,  4),  # GCAN^NNN^NTGC           3-nt 3' overhang
    "BstNI":     ("CCWGG",        2,  3),  # CC^WGG    / WGG^CC      3' overhang
    "BstUI":     ("CGCG",         2,  2),  # CG^CG                   blunt (4-cutter; methylation-sensitive)
    "BstZ17I":   ("GTATAC",       3,  3),  # GTA^TAC                 blunt
    "BtgI":      ("CCRYGG",       1,  5),  # C^CRYGG                 5' overhang
    "Cac8I":     ("GCNNGC",       3,  3),  # GCN^NGC                 blunt
    "CviAII":    ("CATG",         1,  3),  # C^ATG                   NcoI-compatible ends (4-cutter)
    "CviQI":     ("GTAC",         1,  3),  # G^TAC                   KpnI subset (4-cutter)
    "DpnI":      ("GATC",         2,  2),  # GA^TC                   blunt; cuts only methylated
    "DpnII":     ("GATC",         0,  4),  # ^GATC     / GATC^       5' overhang (4-cutter)
    "DrdI":      ("GACNNNNNNGTC", 7,  5),  # GACNNNNN^NGTC           3' overhang
    "EcoO109I":  ("RGGNCCY",      2,  5),  # RG^GNCCY                5' overhang
    "HphI":      ("GGTGA",       13, 12),  # GGTGA(8/7) downstream   Type IIS
    "KasI":      ("GGCGCC",       1,  5),  # G^GCGCC                 5' overhang (NarI isoschizomer)
    "MboI":      ("GATC",         0,  4),  # ^GATC                   DpnII isoschizomer
    "MboII":     ("GAAGA",       13, 12),  # GAAGA(8/7) downstream   Type IIS
    "MlyI":      ("GAGTC",       10, 10),  # GAGTC(5/5) downstream   blunt, Type IIS
    "MmeI":      ("TCCRAC",      26, 24),  # TCCRAC(20/18)           Type IIS far-cutter
    "MspA1I":    ("CMGCKG",       3,  3),  # CMG^CKG                 blunt
    "NgoMIV":    ("GCCGGC",       1,  5),  # G^CCGGC                 5' overhang (EagI-compatible)
    "NmeAIII":   ("GCCGAG",      27, 25),  # GCCGAG(21/19)           Type IIS far-cutter
    "PflMI":     ("CCANNNNNTGG",  7,  4),  # CCANN4^NTGG             3' overhang (BstXI isoschizomer)
    "PspOMI":    ("GGGCCC",       1,  5),  # G^GGCCC                 ApaI isoschizomer (5' overhang)
    "Sau3AI":    ("GATC",         0,  4),  # ^GATC                   BamHI-compatible ends (4-cutter)
    "SfcI":      ("CTRYAG",       1,  5),  # C^TRYAG                 5' overhang
    "SspI":      ("AATATT",       3,  3),  # AAT^ATT                 blunt
    "TaqI":      ("TCGA",         1,  3),  # T^CGA     / CGT^A       5' overhang (heat-stable)
    "Van91I":    ("CCANNNNNTGG",  7,  4),  # PflMI isoschizomer
    "ZraI":      ("GACGTC",       3,  3),  # GAC^GTC                 blunt (AatII-related)

    # ── Type IIS — cut outside recognition sequence ────────────────────────────
    # fwd/rev positions are still offsets from start of recognition seq.
    # For an n-bp recognition sequence cutting d1/d2 downstream:
    #   fwd = n + d1,  rev = n + d2
    "BaeI":      ("ACNNNNGTAYC", -10,-15), # (10/15)…(12/7): upstream pair only — downstream cut not represented (11-bp recog)
    "BbsI":      ("GAAGAC",       8, 12),  # GAAGAC(2/6)  BpiI isoschizomer
    "BcoDI":     ("GTCTC",        6, 10),  # GTCTC(1/5)   BsaI 5-bp variant
    "BceAI":     ("ACGGC",       17, 19),  # ACGGC(12/14)            Type IIS far-cutter
    "BciVI":     ("GTATCC",      12, 11),  # GTATCC(6/5)             1-nt 3' overhang
    "BfuAI":     ("ACCTGC",      10, 14),  # ACCTGC(4/8)  BspMI isoschizomer
    "BmrI":      ("ACTGGG",      11, 10),  # ACTGGG(5/4)             1-nt 3' overhang Type IIS
    "BpiI":      ("GAAGAC",       8, 12),  # BbsI isoschizomer
    "BsaI":      ("GGTCTC",       7, 11),  # GGTCTC(1/5)  Golden Gate workhorse
    "BsaXI":     ("ACNNNNNCTCC", -9,-12),  # (9/12)…(10/7): upstream pair only — downstream cut not represented
    # BsbI removed 2026-05-11 (issue #14, a user): real REBASE id 329
    # but no commercial supplier — users can't actually buy or order this
    # enzyme. Showing it on the map suggests a digest option that doesn't
    # exist in any wet lab. Keep removed unless a vendor starts producing.
    # BseJI removed 2026-05-11: previous entry claimed BbsI-isoschizomer
    # behaviour with site GAAGAC, but real BseJI recognises GATNNNNATC
    # (blunt) per REBASE/NEB. Users who want BbsI behaviour should use
    # `BbsI` / `BpiI` (already in this catalog). If a real BseJI entry
    # is requested, add it back as: ("GATNNNNATC", 5, 5) — blunt 5-bp
    # 5'-OH cut after the 5th base of recognition.
    "BseLI":     ("CCNNNNNNNGG",  7,  4),  # 3' overhang
    "BseMII":    ("CTCAG",       15, 13),  # CTCAG(10/8)             Type IIS far-cutter
    "BseRI":     ("GAGGAG",      16, 14),  # GAGGAG(10/8)            2-nt 3' overhang
    "BsgI":      ("GTGCAG",      22, 20),  # GTGCAG(16/14) far-cutter
    "BslI":      ("CCNNNNNNNGG",  7,  4),  # 3' overhang (BseLI variant)
    "BsmAI":     ("GTCTC",        6, 10),  # GTCTC(1/5)   BsaI isoschizomer (5-bp)
    "BsmBI":     ("CGTCTC",       7, 11),  # CGTCTC(1/5)  Esp3I isoschizomer
    "BsmFI":     ("GGGAC",       15, 19),  # GGGAC(10/14)
    "BsmI":      ("GAATGC",       7,  5),  # GAATGC(1/-1)            2-nt 3' overhang Type IIS
    # BspLU11III removed 2026-05-11: not a real enzyme. Closest REBASE
    # match is BspLU11I (site ACATGT) — PciI isoschizomer — but BspLU11I
    # also has no commercial supplier, and PciI itself is already in this
    # catalog with the correct (1, 5) tuple.
    "BspMI":     ("ACCTGC",      10, 14),  # ACCTGC(4/8)
    "BspQI":     ("GCTCTTC",      8, 11),  # SapI isoschizomer
    "BspTNI":    ("GGTCTC",       7, 11),  # BsaI isoschizomer
    "BsrBI":     ("CCGCTC",       3,  3),  # cuts within recog (special case)
    "BsrDI":     ("GCAATG",       8,  6),  # GCAATG(2/0)
    "BssSI":     ("CACGAG",       1,  5),  # C^ACGAG / GTGCT^C       4-nt 5' overhang
    "BtgZI":     ("GCGATG",      16, 20),  # GCGATG(10/14)
    "BtsCI":     ("GGATG",        7,  5),  # GGATG(2/0)              2-nt 3' overhang
    "BtsI":      ("GCAGTG",       8,  6),  # GCAGTG(2/0)
    "BtsIMutI":  ("CAGTG",        7,  5),  # CAGTG(2/0)              2-nt 3' overhang (canonical REBASE/NEB capitalisation)
    "EarI":      ("CTCTTC",       7, 10),  # CTCTTC(1/4)  SapI-related 3-nt 5' overhang
    "Esp3I":     ("CGTCTC",       7, 11),  # BsmBI isoschizomer
    "PaqCI":     ("CACCTGC",     11, 15),  # CACCTGC(4/8)
    "SapI":      ("GCTCTTC",      8, 11),  # GCTCTTC(1/4)
    "BsmBI-v2":  ("CGTCTC",       7, 11),  # v2/HF variant

    # ── High-Fidelity (HF) and v2 variants — same recognition/cut as canonical ─
    "AgeI-HF":   ("ACCGGT",       1,  5),
    "BamHI-HF":  ("GGATCC",       1,  5),
    "BclI-HF":   ("TGATCA",       1,  5),
    "BmtI":      ("GCTAGC",       5,  1),  # GCTAG^C / G^CTAGC       3' overhang CTAG (NheI neoschizomer)
    "BsiWI-HF":  ("CGTACG",       1,  5),
    "BsrFI-v2":  ("RCCGGY",       1,  5),
    "BsrGI-HF":  ("TGTACA",       1,  5),
    "BssSI-v2":  ("CACGAG",       1,  5),
    "BstEII-HF": ("GGTNACC",      1,  6),
    "BstZ17I-HF":("GTATAC",       3,  3),
    "DraIII-HF": ("CACNNNGTG",    6,  3),
    "EcoRI-HF":  ("GAATTC",       1,  5),
    "EcoRV-HF":  ("GATATC",       3,  3),
    "HindIII-HF":("AAGCTT",       1,  5),
    "KpnI-HF":   ("GGTACC",       5,  1),
    "MfeI-HF":   ("CAATTG",       1,  5),
    "MluI-HF":   ("ACGCGT",       1,  5),
    "MunI-HF":   ("CAATTG",       1,  5),
    "NcoI-HF":   ("CCATGG",       1,  5),
    "NheI-HF":   ("GCTAGC",       1,  5),
    "NotI-HF":   ("GCGGCCGC",     2,  6),
    "NruI-HF":   ("TCGCGA",       3,  3),
    "NsiI-HF":   ("ATGCAT",       5,  1),
    "PstI-HF":   ("CTGCAG",       5,  1),
    "PvuI-HF":   ("CGATCG",       4,  2),
    "PvuII-HF":  ("CAGCTG",       3,  3),
    "SacI-HF":   ("GAGCTC",       5,  1),
    "SalI-HF":   ("GTCGAC",       1,  5),
    "SbfI-HF":   ("CCTGCAGG",     6,  2),
    "ScaI-HF":   ("AGTACT",       3,  3),
    "SpeI-HF":   ("ACTAGT",       1,  5),
    "SphI-HF":   ("GCATGC",       5,  1),
    "TaqI-v2":   ("TCGA",         1,  3),
    "XhoI-HF":   ("CTCGAG",       1,  5),
}


# ── Commercial synonyms: alias → the `_NEB_ENZYMES` name that IS this enzyme ──
#
# The same specificity is sold under different names by different suppliers —
# Thermo/Fermentas ship BsaI as **Eco31I**, SapI as **LguI**, PaqCI as
# **AarI** — and a scan that only knows one spelling answers "no sites" for the
# other. That is the [INV-187] failure shape: a safety check that CANNOT fail.
# Rather than duplicate every synonym into `_NEB_ENZYMES` (which would make the
# per-enzyme scan ~20% slower for zero new biology), the resolver consults this
# table AFTER the catalog, so a user-added custom enzyme of the same name still
# wins.
#
# HARD RULE: an entry here must cut IDENTICALLY to its target — same
# recognition sequence AND same offsets. A NEOschizomer (same site, different
# cut: Acc65I/KpnI, XmaI/SmaI, ApaI/PspOMI, NheI/BmtI, KasI/NarI, AatII/ZraI,
# SacI/Eco53kI) is a DIFFERENT enzyme leaving a DIFFERENT overhang and gets its
# own `_NEB_ENZYMES` entry instead. Aliasing one would be the confidently-wrong
# answer this table exists to prevent — pinned by
# `test_dna_sanity.py::TestEnzymeNameResolution`.
_ENZYME_ALIASES: dict[str, str] = {
    # ── Type IIS — the modular-cloning names ──────────────────────────────────
    "AarI":      "PaqCI",     # CACCTGC(4/8)   Thermo
    "BveI":      "BspMI",     # ACCTGC(4/8)    Thermo (= BfuAI)
    "Eam1104I":  "EarI",      # CTCTTC(1/4)    Thermo (= Ksp632I)
    "Eco31I":    "BsaI",      # GGTCTC(1/5)    Thermo
    "LguI":      "SapI",      # GCTCTTC(1/4)   Thermo (= BspQI)
    # ── Type IIP — true isoschizomers of a catalog entry ──────────────────────
    "Alw44I":    "ApaLI",     # G^TGCAC        Thermo
    "Asp718I":   "Acc65I",    # G^GTACC        Roche
    "BcuI":      "SpeI",      # A^CTAGT        Thermo
    "BshTI":     "AgeI",      # A^CCGGT        Thermo
    "Bsp119I":   "BstBI",     # TT^CGAA        Thermo (= AsuII)
    "Bsp120I":   "PspOMI",    # G^GGCCC        Thermo
    "Bsp1407I":  "BsrGI",     # T^GTACA        Thermo
    "BspDI":     "ClaI",      # AT^CGAT        NEB
    "BspTI":     "AflII",     # C^TTAAG        Thermo
    "Bst1107I":  "BstZ17I",   # GTA^TAC        Thermo
    "Bsu15I":    "ClaI",      # AT^CGAT        Thermo
    "Cfr9I":     "XmaI",      # C^CCGGG        Thermo
    "Cfr42I":    "SacII",     # CCGC^GG        Thermo (= KspI)
    "Ecl136II":  "Eco53kI",   # GAG^CTC        Thermo (= EcoICRI)
    "Eco32I":    "EcoRV",     # GAT^ATC        Thermo
    "Eco72I":    "PmlI",      # CAC^GTG        Thermo
    "Eco91I":    "BstEII",    # G^GTNACC       Thermo (= BstPI)
    "Eco105I":   "SnaBI",     # TAC^GTA        Thermo
    "Eco147I":   "StuI",      # AGG^CCT        Thermo (= AatI)
    "Kpn2I":     "BspEI",     # T^CCGGA        Thermo
    "MlsI":      "MscI",      # TGG^CCA        Thermo
    "Mph1103I":  "NsiI",      # ATGCA^T        Thermo
    "MssI":      "PmeI",      # GTTT^AAAC      Thermo
    "NdeII":     "DpnII",     # ^GATC          (= MboI / Sau3AI)
    "PaeI":      "SphI",      # GCATG^C        Thermo
    "PauI":      "BssHII",    # G^CGCGC        Thermo
    "PdmI":      "XmnI",      # GAANN^NNTTC    Thermo
    "Pfl23II":   "BsiWI",     # C^GTACG        Thermo (= SplI)
    "Psp1406I":  "AclI",      # AA^CGTT        Thermo
    "SgsI":      "AscI",      # GG^CGCGCC      Thermo
    "SmiI":      "SwaI",      # ATTT^AAAT      Thermo
    "SstI":      "SacI",      # GAGCT^C        classic
    "SstII":     "SacII",     # CCGC^GG        classic
    "TspMI":     "XmaI",      # C^CCGGG        NEB
}


_GB_L0_PARTS: list[tuple] = [
    # ── Promoters (Pos 1, combined Promoter+5'UTR: GGAG → AATG) ───────
    ("CaMV 35S",          "Promoter",   "Pos 1",   "GGAG", "AATG", "pUPD2", "Spectinomycin"),
    ("Nos",               "Promoter",   "Pos 1",   "GGAG", "AATG", "pUPD2", "Spectinomycin"),
    ("AtUBQ10",           "Promoter",   "Pos 1",   "GGAG", "AATG", "pUPD2", "Spectinomycin"),
    ("ZmUBI1",            "Promoter",   "Pos 1",   "GGAG", "AATG", "pUPD2", "Spectinomycin"),
    ("AtRPS5a",           "Promoter",   "Pos 1",   "GGAG", "AATG", "pUPD2", "Spectinomycin"),
    # ── CDS with stop (Positions 3-4: AATG → GCTT) ─────────────────────
    ("eGFP",              "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("mCherry",           "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("mVenus",            "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("mTurquoise2",       "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("GUS (uidA)",        "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("Luciferase (LUC+)", "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("NptII (KanR)",      "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("hptII (HygR)",      "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("Bar (BastaR)",      "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    ("Cas9 (SpCas9)",     "CDS",        "Pos 3-4", "AATG", "GCTT", "pUPD2", "Spectinomycin"),
    # ── CDS without stop (Position 3: AATG → TTCG) ─────────────────────
    ("eGFP (no stop)",    "CDS-NS",     "Pos 3",   "AATG", "TTCG", "pUPD2", "Spectinomycin"),
    ("mCherry (no stop)", "CDS-NS",     "Pos 3",   "AATG", "TTCG", "pUPD2", "Spectinomycin"),
    # ── C-terminal tags (Position 4: TTCG → GCTT) ──────────────────────
    ("GFP C-tag",         "C-tag",      "Pos 4",   "TTCG", "GCTT", "pUPD2", "Spectinomycin"),
    ("HA tag",            "C-tag",      "Pos 4",   "TTCG", "GCTT", "pUPD2", "Spectinomycin"),
    ("6xHis tag",         "C-tag",      "Pos 4",   "TTCG", "GCTT", "pUPD2", "Spectinomycin"),
    # ── Terminators (Position 5: GCTT → CGCT) ──────────────────────────
    ("Nos terminator",    "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
    ("CaMV 35S term",     "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
    ("OCS terminator",    "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
    ("rbcS terminator",   "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
    ("HSP18.2 term",      "Terminator", "Pos 5",   "GCTT", "CGCT", "pUPD2", "Spectinomycin"),
]


_GB_TYPE_COLORS: dict[str, str] = {
    # Legacy slots:
    "Promoter":         "green",
    "Promoter-only":    "green",
    "5' UTR":           "cyan",
    "CDS":              "yellow",
    "OPERON":           "yellow",      # whole operon as a CDS-equivalent part
    "CDS-NS":           "dark_orange",
    "C-tag":            "magenta",
    "Terminator":       "blue",
    # GB 2.0 canonical additions (2026-05-10) — pick colors that group
    # with their nearest legacy equivalent so the palette stays
    # readable (5'NT shades of green, 5'UTR shades of cyan, translated
    # shades of yellow/orange, 3'NT shades of blue).
    "Operator-A":       "dark_green",   # OP-PROM 5'NT variant
    "Operator-B":       "dark_green",
    "Min Promoter":     "green",
    "Distal 5' UTR":    "cyan",
    "Signal peptide":   "dark_orange",  # N-terminal coding extension
    "CDS-NS (CT)":      "dark_orange",  # canonical no-stop CDS
    "CT-tag":           "magenta",      # canonical C-tag
    "CDS-after-SP":     "yellow",       # full CDS after SP cleavage
    "3' UTR":           "blue",
    "Terminator-only":  "blue",
}


_GB_POSITIONS: dict[str, tuple[str, str, str]] = {
    # ── 5' Non-Transcribed (combined Promoter + 5'UTR + Link) ─────────
    # Combined Promoter+5'UTR+ATG (GB 2.0 PromUTR; the common BASIC form).
    # Legacy position label "Pos 1" preserved for back-compat with
    # existing user parts in `parts_bin.json` that were stored under
    # this label string.
    "Promoter":         ("Pos 1",     "GGAG", "AATG"),
    # Separate Promoter (no LINK/+ATG) — pairs with a `5' UTR` part.
    "Promoter-only":    ("Pos 1a",    "GGAG", "CCAT"),
    # ── Operator/Promoter variants (OP-PROM-A/B workflow) ─────────────
    "Operator-A":       ("Pos 01-02", "GGAG", "TCCC"),  # OP-PROM-A operator
    "Operator-B":       ("Pos 02",    "TGAC", "TCCC"),  # OP-PROM-B operator
    "Min Promoter":     ("Pos 03-12", "TCCC", "AATG"),  # pairs with either Operator
    # ── 5' UTR ────────────────────────────────────────────────────────
    # Historical "5' UTR" slot — technically the LINK position (Pos 12)
    # in canonical GB 2.0. Kept under this name for backward compat.
    "5' UTR":           ("Pos 1b",    "CCAT", "AATG"),
    # Canonical GB 2.0 5'UTR (Pos 03-11), distinct from the LINK above.
    "Distal 5' UTR":    ("Pos 03-11", "TCCC", "CCAT"),
    # ── Translated region ─────────────────────────────────────────────
    # Signal peptide (SECRETED workflow N-terminal extension; coding).
    "Signal peptide":   ("Pos 13",    "AATG", "AGCC"),
    # Full CDS with stop codon (BASIC workflow).
    "CDS":              ("Pos 3-4",   "AATG", "GCTT"),
    # A whole NATIVE OPERON as one CDS-equivalent L0 part — same fusion-site
    # slot + overhangs as CDS (Pos 3-4, AATG→GCTT): the AATG overhang carries
    # the FIRST gene's start codon (ATG-fusion, exactly like a CDS) and GCTT
    # meets a downstream terminator after the LAST gene's stop. Used by the
    # Native Operon Domestication workbench so a cured operon clones straight
    # into Golden Braid between a promoter and a terminator.
    "OPERON":           ("Pos 3-4",   "AATG", "GCTT"),
    # Legacy 2-part CDS split (splits at TL2/TL3 boundary, non-canonical
    # but kept for backward compat with existing user parts).
    "CDS-NS":           ("Pos 3",     "AATG", "TTCG"),
    "C-tag":            ("Pos 4",     "TTCG", "GCTT"),
    # Canonical GB 2.0 CDS variants (added 2026-05-10):
    # CDS-no-stop for CT-FUSION (positions 13-15; pairs with `CT-tag`).
    "CDS-NS (CT)":      ("Pos 13-15", "AATG", "GCAG"),
    "CT-tag":           ("Pos 16",    "GCAG", "GCTT"),  # canonical C-term tag
    # CDS body for SECRETED workflow (after Signal peptide, full to stop).
    "CDS-after-SP":     ("Pos 14-16", "AGCC", "GCTT"),
    # ── 3' Non-Translated ─────────────────────────────────────────────
    # Combined 3'UTR+Terminator (BASIC workflow).
    "Terminator":       ("Pos 5",     "GCTT", "CGCT"),
    # Canonical split variants (added 2026-05-10):
    "3' UTR":           ("Pos 17",    "GCTT", "GGTA"),
    "Terminator-only":  ("Pos 21",    "GGTA", "CGCT"),
}


_GB_CODING_PART_TYPES: frozenset[str] = frozenset({
    # Legacy splits:
    "CDS", "CDS-NS", "C-tag",
    # GB 2.0 canonical translational parts (added 2026-05-10):
    "Signal peptide", "CDS-NS (CT)", "CT-tag", "CDS-after-SP",
    # A whole native operon cloned as one CDS-equivalent L0 part — "coding" so
    # the AATG overhang carries the first gene's ATG (the ATG-fusion skip in
    # `_atg_offset_for_part`); curing itself is per-CDS in the SOE designer.
    "OPERON",
})


_GB_PART_TYPE_TO_INSDC: dict[str, str] = {
    # Legacy slots:
    "Promoter":         "promoter",
    "Promoter-only":    "promoter",
    "5' UTR":           "5'UTR",
    "CDS":              "CDS",
    "OPERON":           "operon",
    "CDS-NS":           "CDS",
    "C-tag":            "CDS",
    "Terminator":       "terminator",
    # GB 2.0 canonical additions (2026-05-10):
    "Operator-A":       "promoter",   # operator+promoter assembled together
    "Operator-B":       "promoter",
    "Min Promoter":     "promoter",
    "Distal 5' UTR":    "5'UTR",
    "Signal peptide":   "sig_peptide",
    "CDS-NS (CT)":      "CDS",
    "CT-tag":           "CDS",
    "CDS-after-SP":     "CDS",
    "3' UTR":           "3'UTR",
    "Terminator-only":  "terminator",
}


_GB_L0_ENZYME_NAME = "Esp3I"       # Esp3I is the isoschizomer of BsmBI


_GB_L0_ENZYME_SITE = "CGTCTC"      # recognition; rc = "GAGACG"


_GB_SPACER         = "A"           # 1 nt between recognition and the overhang


_GB_PAD            = "GCGC"        # 4 nt of extra bases for efficient end-cutting


_GB_DOMESTICATION_FORBIDDEN: dict[str, str] = {
    # Esp3I self-cuts during L0 domestication; BsaI would re-cut during any
    # downstream L1 assembly — both must be absent from the final part.
    "BsaI":  "GGTCTC",
    "Esp3I": "CGTCTC",
}


_BUILTIN_GRAMMARS: dict[str, dict] = {
    "gb_l0": {
        "id":              "gb_l0",
        "name":            "Golden Braid L0",
        "enzyme":          _GB_L0_ENZYME_NAME,
        # Iterative GB cycle: L0 → L1 cuts with Esp3I (`enzyme`); L1 →
        # L2 cuts with BsaI (`level_up_enzyme`); L2 → L3 wraps around
        # to Esp3I again (parity on source level — see
        # `_enzyme_for_level_up`). The two enzymes alternate each
        # level so the assembled product survives the next cut.
        "level_up_enzyme": "BsaI",
        "site":            _GB_L0_ENZYME_SITE,
        "spacer":          _GB_SPACER,
        "pad":             _GB_PAD,
        "forbidden_sites": dict(_GB_DOMESTICATION_FORBIDDEN),
        "positions": [
            {"name": pos, "type": ptype, "oh5": oh5, "oh3": oh3,
             "color": _GB_TYPE_COLORS.get(ptype, "white")}
            for ptype, (pos, oh5, oh3) in _GB_POSITIONS.items()
        ],
        "coding_types":    sorted(_GB_CODING_PART_TYPES),
        "type_to_insdc":   dict(_GB_PART_TYPE_TO_INSDC),
        "catalog":         list(_GB_L0_PARTS),
        "editable":        False,
    },
    # Plant MoClo (Weber et al. 2011, Engler et al. 2014). BsaI at L0,
    # BpiI/BbsI at L1 — both scrubbed during domestication. Ships
    # without a built-in catalog because Plant MoClo's reference parts
    # depend heavily on the user's host system; users seed via "New
    # Part" or by duplicating into a custom grammar.
    "moclo_plant": {
        "id":              "moclo_plant",
        "name":            "MoClo Plant (Weber 2011)",
        "enzyme":          "BsaI",
        "level_up_enzyme": "BpiI",
        "site":            "GGTCTC",
        "spacer":          "A",
        "pad":             "GCGC",
        # BsaI for the current L0 cut; BpiI (= BbsI) for the next-level
        # MoClo assembly, which uses a different Type IIS site so the
        # L0 part survives the L1 reaction without re-cutting.
        "forbidden_sites": {"BsaI": "GGTCTC", "BpiI": "GAAGAC"},
        "positions": [
            {"name": "Pos 1", "type": "Promoter",   "oh5": "GGAG", "oh3": "AATG", "color": "green"},
            {"name": "Pos 2", "type": "5' UTR",     "oh5": "AATG", "oh3": "AGGT", "color": "cyan"},
            {"name": "Pos 3", "type": "CDS",        "oh5": "AGGT", "oh3": "GCTT", "color": "yellow"},
            {"name": "Pos 4", "type": "C-tag",      "oh5": "GCTT", "oh3": "GGTA", "color": "magenta"},
            {"name": "Pos 5", "type": "Terminator", "oh5": "GGTA", "oh3": "CGCT", "color": "blue"},
        ],
        "coding_types":    ["CDS", "C-tag"],
        "type_to_insdc": {
            "Promoter":   "promoter",
            "5' UTR":     "5'UTR",
            "CDS":        "CDS",
            "C-tag":      "CDS",
            "Terminator": "terminator",
        },
        "catalog":         [],
        "editable":        False,
    },
}


def _all_grammars() -> dict[str, dict]:
    """Return all grammars (built-in + user-defined) keyed by id.

    Built-ins come first; user-defined grammars override builtin IDs
    if they ever collide (defensive — UI prevents this on save). The
    returned dicts are independent copies, so callers may mutate them
    without poisoning the cache.
    """
    out: dict[str, dict] = {gid: deepcopy(g) for gid, g in _BUILTIN_GRAMMARS.items()}
    for g in _load_custom_grammars():
        gid = g.get("id")
        if isinstance(gid, str):
            # Custom grammars are always editable regardless of what
            # the JSON file says — stops a mis-flagged file from
            # locking the user out of their own definitions.
            g = dict(g)
            g["editable"] = True
            out[gid] = g
    return out


def _get_entry_vector(
    grammar_id: str, role: str = "",
) -> "dict | None":
    """Return the entry-vector dict for ``(grammar_id, role)``, or
    None if none has been assigned yet.

    `role` is the per-grammar slot — for Golden Braid the Constructor
    has four roles (``Alpha1``, ``Alpha2``, ``Omega1``, ``Omega2``)
    so a single grammar carries multiple L1 acceptors. The empty
    role (default) is the singleton L0 entry vector used by the
    Domesticator / `_clone_part_into_entry_vector` workflow — kept
    backward-compat with pre-2026-05-07 entry_vectors.json files
    where every entry has no `role` field.
    """
    if not isinstance(grammar_id, str) or not grammar_id:
        return None
    role = role or ""
    for e in _load_entry_vectors():
        if e.get("grammar_id") == grammar_id and (e.get("role") or "") == role:
            return e
    return None


# ── Selection-marker detection + grammar dropdown options (moved, Phase D) ──
# `_detect_selection_marker` parses gb_text via `_gb_text_to_record` (record sibling,
# same layer); `_grammar_dropdown_options` reads `_all_grammars`/`_BUILTIN_GRAMMARS`.
_SELECTION_MARKER_KEYWORDS: "tuple[tuple[str, str], ...]" = (
    ("ampicillin",     "Ampicillin"),
    ("ampr",           "Ampicillin"),
    ("kanamycin",      "Kanamycin"),
    ("kanr",           "Kanamycin"),
    ("neomycin",       "Kanamycin"),
    ("neor",           "Kanamycin"),
    ("spectinomycin",  "Spectinomycin"),
    ("specr",          "Spectinomycin"),
    ("aada",           "Spectinomycin"),
    ("smr",            "Spectinomycin"),
    ("chloramphenicol", "Chloramphenicol"),
    ("cmr",            "Chloramphenicol"),
    ("cat",            "Chloramphenicol"),
    ("tetracycline",   "Tetracycline"),
    ("tetr",           "Tetracycline"),
    ("hygromycin",     "Hygromycin"),
    ("hygr",           "Hygromycin"),
    ("zeocin",         "Zeocin"),
    ("zeor",           "Zeocin"),
    ("gentamicin",     "Gentamicin"),
    ("gmr",            "Gentamicin"),
    ("erythromycin",   "Erythromycin"),
    ("ermr",           "Erythromycin"),
    # ermB/ermC are the standard Gram-positive MLS markers and the actual
    # selection on common Lactobacillus/Lactococcus shuttle vectors — but
    # only "ermr" was listed, and matching is whole-token, so a backbone
    # labelled "ermB" reported NO marker at all (found 2026-08-05).
    ("erm",            "Erythromycin"),
    ("ermb",           "Erythromycin"),
    ("ermc",           "Erythromycin"),
    ("ermam",          "Erythromycin"),
    ("bla",            "Ampicillin"),   # last — matches "bla" alone
)


def _detect_selection_marker(gb_text: str) -> "str | None":
    """Scan a plasmid's GenBank text for a feature whose
    label / gene / product / note contains one of
    ``_SELECTION_MARKER_KEYWORDS`` and return the matching display
    name (e.g. ``"Kanamycin"``). Returns ``None`` when no recognised
    marker is found. Used by the Domesticator's part-save path so
    the saved `marker` field reflects the user's configured entry
    vector instead of the historical pUPD2 / Spectinomycin defaults.

    Conservative on parse failure: malformed or empty `gb_text`
    yields ``None`` so the caller can fall back to a placeholder
    rather than asserting a marker we can't actually verify.
    """
    if not isinstance(gb_text, str) or not gb_text.strip():
        return None
    try:
        rec = _gb_text_to_record(gb_text)
    except Exception:
        return None
    qual_keys = ("label", "gene", "product", "note", "standard_name")
    # Sweep #34 (2026-05-26): word-boundary match instead of bare
    # substring. Pre-fix a CDS labelled `"category"` would hit
    # `"cat" in "category"` → returned "Chloramphenicol" silently.
    # Same false-positive vector with `"bla"` inside `"blast"`,
    # `"smr"` inside `"smrt-seq"`, etc. Tokenise on non-alnum
    # separators so the match needs an isolated word — false-
    # positives shrink to near-zero while real labels (`cat`,
    # `bla`, `kanR-cat-tet`) still resolve.
    for feat in getattr(rec, "features", []) or []:
        bag: list[str] = []
        for k in qual_keys:
            vals = feat.qualifiers.get(k) if hasattr(feat, "qualifiers") \
                else None
            if isinstance(vals, list):
                bag.extend(str(v) for v in vals)
            elif isinstance(vals, str):
                bag.append(vals)
        for s in bag:
            tokens = {
                t.lower()
                for t in re.split(r"[^A-Za-z0-9]+", s)
                if t
            }
            if not tokens:
                continue
            for kw, display in _SELECTION_MARKER_KEYWORDS:
                if kw in tokens:
                    return display
    return None


def _grammar_dropdown_options() -> list[tuple[str, str]]:
    """Return ``[(display_name, id), …]`` for every grammar, in the
    canonical order used wherever a Select dropdown lists grammars
    (DomesticatorModal "Grammar" picker today; future menus likely):

      1. **Golden Braid L0 first** — the default reference grammar.
         Pinned at position 1 regardless of any other ordering
         shenanigans (e.g., a custom grammar id-sorted before
         ``gb_l0``).
      2. Other built-in grammars (MoClo Plant, etc.) in
         ``_BUILTIN_GRAMMARS`` insertion order.
      3. Custom grammars from ``cloning_grammars.json`` last, tagged
         ``(custom)`` for visual disambiguation.
    """
    grammars = _all_grammars()
    out: list[tuple[str, str]] = []
    if "gb_l0" in grammars:
        g = grammars["gb_l0"]
        out.append((f"{g.get('name', 'Golden Braid L0')}", "gb_l0"))
    for gid in _BUILTIN_GRAMMARS:
        if gid == "gb_l0" or gid not in grammars:
            continue
        g = grammars[gid]
        out.append((f"{g.get('name', gid)}", gid))
    for gid, g in grammars.items():
        if gid in _BUILTIN_GRAMMARS:
            continue
        out.append((f"{g.get('name', gid)}  (custom)", gid))
    return out


# ── Feature-entry extraction from a record (moved from hub, Phase D) ────────
def _extract_feature_entries_from_record(record) -> list[dict]:
    """Return one feature-library entry dict per non-source feature.

    Wrap features (origin-spanning CompoundLocations) are flattened into the
    forward-strand genomic sequence before export, so the entry's ``sequence``
    is always the 5'→3' DNA that would be re-inserted. Reverse-strand
    features store the revcomp (i.e. the 5'→3' of the feature as read), which
    matches how the Add Feature modal expects input.
    """
    try:
        from Bio.SeqFeature import CompoundLocation
    except ImportError:
        CompoundLocation = tuple()  # type: ignore[assignment]
    seq = str(getattr(record, "seq", "") or "").upper()
    total = len(seq)
    entries: list[dict] = []
    for feat in getattr(record, "features", []) or []:
        if feat.type == "source":
            continue
        loc = feat.location
        strand = getattr(loc, "strand", 1) or 1
        # Assemble the forward-strand genomic sequence under the feature,
        # respecting wrap/compound locations.
        if isinstance(loc, CompoundLocation):
            parts_seq = []
            for part in loc.parts:
                s = int(part.start) % total if total else 0
                e = int(part.end)   % (total or 1) if total else 0
                if total and e <= s:
                    parts_seq.append(seq[s:] + seq[:e])
                else:
                    parts_seq.append(seq[s:e])
            fwd = "".join(parts_seq)
        else:
            s = int(loc.start)
            e = int(loc.end)
            fwd = seq[s:e]
        # Store 5'→3' of the feature as read. For reverse-strand CDS that is
        # the revcomp of the genomic slice.
        if strand == -1 and fwd:
            feat_seq = _rc(fwd)
        else:
            feat_seq = fwd
        entries.append({
            "name":         _feat_label(feat),
            "feature_type": feat.type,
            "sequence":     feat_seq,
            "strand":       1 if strand != -1 else -1,
            "qualifiers":   {k: list(v) if isinstance(v, (list, tuple)) else [v]
                             for k, v in (feat.qualifiers or {}).items()},
            "description":  "",
        })
    return entries


# ── Gel/HMM/codon/entry-vector derived accessors (moved from hub, Phase D) ──
def _clear_entry_vectors_for_grammar(grammar_id: str) -> int:
    """Drop every entry-vector binding for ``grammar_id``, regardless
    of role. Used by grammar-delete flows so a grammar with multiple
    role bindings (Constructor's Alpha1/Alpha2/Omega1/Omega2) doesn't
    leave orphan rows in entry_vectors.json that would silently come
    back to life if a future custom grammar slugs to the same id.

    Returns the number of bindings removed; 0 when nothing matched.
    No-op for empty grammar_id (defensive — callers can pass a
    user-typed value without pre-checking).

    Sweep #10 (2026-05-20): RMW under `_cache_lock`.
    """
    if not isinstance(grammar_id, str) or not grammar_id:
        return 0
    with _state._cache_lock:
        entries = _load_entry_vectors()
        kept = [e for e in entries if e.get("grammar_id") != grammar_id]
        n_dropped = len(entries) - len(kept)
        if n_dropped > 0:
            _save_entry_vectors(kept)
    return n_dropped


def _codon_tables_get(key: str) -> "dict | None":
    """Look up a table by taxid or name (case-insensitive). Returns the
    in-memory entry (with 'raw' as tuples) or None."""
    key = (key or "").strip().lower()
    if not key:
        return None
    for e in _codon_tables_load():
        if str(e.get("taxid", "")).lower() == key:
            return e
        if str(e.get("name", "")).lower() == key:
            return e
    return None


def _find_hmm_db_entry(entry_id: str) -> "dict | None":
    """Return a deep-cloned catalog entry by id, or None on miss."""
    safe = _sanitize_hmm_db_id(entry_id)
    if safe is None:
        return None
    for e in _load_hmm_db_catalog():
        if e.get("id") == safe:
            return e
    return None


def _hmm_db_name_taken(name: str,
                         exclude_id: "str | None" = None) -> bool:
    """Case-insensitive name lookup. The UI uses display name as the
    "is this taken?" check (matches collection / gel modal conventions)
    so two entries can't share a display name even if they have
    different ids."""
    if not isinstance(name, str):
        return False
    target = name.strip().casefold()
    if not target:
        return False
    for e in _load_hmm_db_catalog():
        if exclude_id is not None and e.get("id") == exclude_id:
            continue
        if (e.get("name") or "").strip().casefold() == target:
            return True
    return False


def _gel_name_taken(name: str) -> bool:
    """Dup-name guard. Trims + case-sensitive compare so 'Friday
    digest' and 'Friday digest ' coexist gracefully but two identical
    names won't (matches the user's intuition that names are visible
    identifiers)."""
    if not isinstance(name, str):
        return False
    needle = name.strip()
    if not needle:
        return False
    for g in _load_gels():
        if (g.get("name") or "").strip() == needle:
            return True
    return False


# ── Gel finder (moved from hub, Phase D) ────────────────────────────────────
def _find_gel(gel_id: str) -> "dict | None":
    safe = _sanitize_gel_id(gel_id)
    if safe is None:
        return None
    for g in _load_gels():
        if g.get("id") == safe:
            return g
    return None


# ── Entry-vector setter + parts-bin readonly view (moved from hub, Phase D) ─
def _set_entry_vector(
    grammar_id: str, vector: "dict | None", role: str = "",
) -> None:
    """Replace (or remove, with ``vector=None``) the entry vector
    bound to ``(grammar_id, role)``. Persists via `_safe_save_json`.

    Only the matching ``(grammar_id, role)`` entry is replaced —
    other roles for the same grammar (e.g. Alpha2 / Omega1 / Omega2
    when setting Alpha1) survive untouched.

    Sweep #10 (2026-05-20): full load-modify-save sequence runs under
    `_cache_lock` so concurrent rapid clicks across grammar tabs (or
    `_auto_bind_entry_vectors_from_entries` looping per-role) can't
    lose updates via classic RMW racing. RLock allows re-entry from
    the nested `_save_entry_vectors` lock acquisition.
    """
    if not isinstance(grammar_id, str) or not grammar_id:
        return
    role = role or ""
    with _state._cache_lock:
        entries = [
            e for e in _load_entry_vectors()
            if not (
                e.get("grammar_id") == grammar_id
                and (e.get("role") or "") == role
            )
        ]
        if vector is not None:
            v = dict(vector)
            v["grammar_id"] = grammar_id
            v["role"]       = role
            entries.append(v)
        _save_entry_vectors(entries)


def _iter_parts_bin_readonly() -> list[dict]:
    """Read-only view of the parts-bin cache — returns the cached list
    directly without cloning. Callers MUST NOT mutate (pitfall #17).
    Use for hot read paths (Constructor palette repop, classifier
    probes) where the per-row deepcopy of qualifier dicts is wasted
    when the caller only reads fields. For any path that mutates
    entries, use ``_load_parts_bin()``. Added by sweep #25
    (2026-05-23) mirroring ``_iter_library_readonly`` /
    ``_iter_collections_readonly``.
    """
    _ensure = _state._ensure_parts_bin_hook
    if _ensure is not None:
        _ensure()
    assert _state._parts_bin_cache is not None
    return _state._parts_bin_cache


# ── Active-primer-collection pointer + library/collection fuzzy search (Phase D)
def _search_collections_library(query: str, *,
                                  limit: int = 200) -> list[dict]:
    """Cross-collection plasmid search. Returns up to `limit` matches
    across every collection on disk. Each match is
    ``{collection, name, id, size, status, n_feats}``. The query is
    fuzzy-matched against the plasmid name (same `_fuzzy_match`
    semantics as the library panel's in-collection search) so a user
    typing `puc` finds `pUC19`. An empty `query` returns the first
    `limit` plasmids across all collections (useful as a "list
    everything" affordance for agents).

    Used by the cross-collection search modal (`LibrarySearchModal`,
    File → Find plasmid (all collections)…) and the `search-library`
    agent endpoint.
    Results are natural-sorted by ``(collection_name, plasmid_name)``
    so ``Backbones 2`` lands before ``Backbones 10`` and within a
    collection ``pBin2`` lands before ``pBin10``.
    """
    out: list[dict] = []
    q = (query or "").strip()
    # Read-only walk — we copy each match's scalar fields into a fresh
    # dict below, so the cached entries themselves are never exposed
    # to the caller. Lets us skip the per-call `_typed_clone` of the
    # entire collections cache (~30–50 ms on a 5k-plasmid library).
    for c in _iter_collections_readonly():
        cname = c.get("name") or "?"
        for entry in c.get("plasmids", []) or []:
            if not isinstance(entry, dict):
                continue
            entry_id = entry.get("id") or ""
            # Skip entries without an `id` — without one we can't
            # round-trip the dismiss payload to `_apply_record`, and
            # the loader would alias every untagged entry to the
            # first one it finds.
            if not entry_id:
                continue
            name = entry.get("name") or entry_id
            if q and not _fuzzy_match(q, name):
                continue
            out.append({
                "collection": cname,
                "name":       name,
                "id":         entry_id,
                "size":       entry.get("size", 0),
                "status":     _sanitize_plasmid_status(entry.get("status")),
                "n_feats":    entry.get("n_feats", 0),
            })
    # Natural-sort by (collection, name) — group rows from the same
    # collection together, with plasmids inside each group in
    # human-friendly numeric order. `heapq.nsmallest` keeps this O(N
    # + limit·log(limit)) instead of the full O(N·log(N)) `out.sort`
    # would cost, which matters once the user crosses ~1k plasmids
    # (typical limit=200). For small libraries the two paths are
    # indistinguishable; for a 5k-plasmid catalog the heap form saves
    # ~30 ms per search keystroke.
    import heapq as _heapq
    return _heapq.nsmallest(
        limit, out,
        key=lambda r: (
            _natural_sort_key(r.get("collection") or ""),
            _natural_sort_key(r.get("name") or ""),
        ),
    )


def _get_active_primer_collection_name() -> "str | None":
    val = _get_setting("active_primer_collection", None)
    return val if isinstance(val, str) and val else None


def _set_active_primer_collection_name(name: "str | None") -> None:
    prev = _get_active_primer_collection_name()
    target = name or ""
    if prev != target:
        _log_event(
            "primer_collection.switched",
            prev=prev or "", new=target,
        )
    _set_setting("active_primer_collection", target)


# ── Active-name getters/setters + finders (Phase D: project / enzyme-collection /
# parts-bin). Settings-backed data ops, same shape as the active-primer-collection
# pair above; relocated from the hub so the agent endpoints + hub callers resolve
# them here. The enzyme catalog is reached via `_state._all_enzymes_hook`.
def _find_enzyme_collection(name: str) -> "dict | None":
    """Return the collection dict matching ``name`` (case-sensitive), or
    None. Operates on the readonly cache view to avoid the deepcopy cost
    on hot read paths (modal grammar selector, restriction-scan dispatch)."""
    if not isinstance(name, str) or not name:
        return None
    if _state._enzyme_collections_cache is None:
        _load_enzyme_collections()  # populates the cache as a side effect
    for entry in (_state._enzyme_collections_cache or []):
        if isinstance(entry, dict) and entry.get("name") == name:
            return _typed_clone(entry)
    return None


def _get_active_enzyme_collection_name() -> "str | None":
    """Return the active enzyme-collection name, or None when the user
    hasn't picked one (= scan against the full NEB master catalog)."""
    val = _get_setting("active_enzyme_collection", None)
    return val if isinstance(val, str) and val else None


def _set_active_enzyme_collection_name(name: "str | None") -> None:
    """Persist (or clear) the active enzyme-collection pointer."""
    prev = _get_active_enzyme_collection_name()
    target = name or ""
    if prev != target:
        _log_event(
            "enzyme_collection.switched",
            prev=prev or "", new=target,
        )
    _set_setting("active_enzyme_collection", target)


def _active_enzyme_allowed_set() -> "frozenset[str] | None":
    """Resolve the enzyme allow-list for the active enzyme collection.

    Returns ``None`` when no collection is active or the active one is
    empty/missing — caller treats None as "scan the full master
    catalog" (built-in NEB ∪ user-added custom enzymes). Filters
    against `_state._all_enzymes_hook()` so a stale name in a collection (manually
    edited JSON, deleted built-in, or a custom enzyme since removed)
    never reaches the scanner."""
    name = _get_active_enzyme_collection_name()
    if not name:
        return None
    coll = _find_enzyme_collection(name)
    if not coll:
        return None
    enzymes = coll.get("enzymes") or []
    if not isinstance(enzymes, list):
        return None
    known = _state._all_enzymes_hook()
    parsed = frozenset(
        n for n in enzymes
        if isinstance(n, str) and n in known
    )
    return parsed or None


def _get_active_parts_bin_name() -> "str | None":
    val = _get_setting("active_parts_bin", None)
    return val if isinstance(val, str) and val else None


def _set_active_parts_bin_name(name: "str | None") -> None:
    _set_setting("active_parts_bin", name or "")


def _find_parts_bin(name: str) -> "dict | None":
    for b in _load_parts_bin_collections():
        if b.get("name") == name:
            return b
    return None


def _get_active_project_name() -> "str | None":
    val = _get_setting("active_project", None)
    return val if isinstance(val, str) and val else None


def _set_active_project_name(name: "str | None") -> None:
    """Persist (or clear) the active-project pointer. Emits a
    `project.switched` event on actual change so debug bundles
    capture the switch context (mirrors `collection.switched`)."""
    prev = _get_active_project_name()
    target = name or ""
    if prev != target:
        _log_event("project.switched", prev=prev or "", new=target)
    _set_setting("active_project", target)


def _find_project(name: str) -> "dict | None":
    for p in _load_experiment_projects():
        if p.get("name") == name:
            return p
    return None


# ── OT-2 protocol library + custom-labware library ──────────────────────────────
# Two single-file stores, each a list of named collections that EMBED their items
# (saved protocol designs / custom Opentrons labware definitions) — the plasmid-
# collections shape. Saved through the sacred `_safe_save_json` engine; wired into
# the backup / master-delete / migrate registries + the pytest sandbox (see
# `splicecraft_backup._PERSISTENT_CACHE_REGISTRY`). The active-collection pointer
# is a persisted setting.

def _load_protocol_collections() -> "list[dict]":
    """Every saved-protocol collection: ``{name, description?, protocols: [...]}``."""
    if _state._protocol_collections_cache is not None:
        return _typed_clone(_state._protocol_collections_cache)
    with _state._cache_lock:
        if _state._protocol_collections_cache is None:
            entries, warning = _safe_load_json(
                _state._PROTOCOL_COLLECTIONS_FILE, "Protocol collections")
            if warning:
                _log.warning(warning)
            _state._protocol_collections_cache = [e for e in entries if isinstance(e, dict)]
        return _typed_clone(_state._protocol_collections_cache)


def _save_protocol_collections(entries: "list[dict]") -> None:
    with _state._cache_lock:
        _safe_save_json(_state._PROTOCOL_COLLECTIONS_FILE, entries, "Protocol collections")
        _state._protocol_collections_cache = _typed_clone(entries)


def _find_protocol_collection(name: str) -> "dict | None":
    key = (name or "").strip().casefold()
    for c in _load_protocol_collections():
        if (c.get("name") or "").strip().casefold() == key:
            return c
    return None


def _get_active_protocol_collection_name() -> str:
    return str(_get_setting("active_protocol_collection", "") or "")


def _set_active_protocol_collection_name(name: str) -> None:
    _set_setting("active_protocol_collection", str(name or ""))


def _load_custom_labware() -> "list[dict]":
    """Every custom-labware collection: ``{name, description?, labware: [...]}``,
    each labware ``{name, definition, ...}``."""
    if _state._custom_labware_cache is not None:
        return _typed_clone(_state._custom_labware_cache)
    with _state._cache_lock:
        if _state._custom_labware_cache is None:
            entries, warning = _safe_load_json(_state._CUSTOM_LABWARE_FILE, "Custom labware")
            if warning:
                _log.warning(warning)
            _state._custom_labware_cache = [e for e in entries if isinstance(e, dict)]
        return _typed_clone(_state._custom_labware_cache)


def _save_custom_labware(entries: "list[dict]") -> None:
    with _state._cache_lock:
        _safe_save_json(_state._CUSTOM_LABWARE_FILE, entries, "Custom labware")
        _state._custom_labware_cache = _typed_clone(entries)


def _find_labware_collection(name: str) -> "dict | None":
    key = (name or "").strip().casefold()
    for c in _load_custom_labware():
        if (c.get("name") or "").strip().casefold() == key:
            return c
    return None


def _get_active_labware_collection_name() -> str:
    return str(_get_setting("active_labware_collection", "") or "")


def _set_active_labware_collection_name(name: str) -> None:
    _set_setting("active_labware_collection", str(name or ""))
