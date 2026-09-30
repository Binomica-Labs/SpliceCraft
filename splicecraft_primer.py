"""splicecraft_primer — primer / mutagenesis design (Phase D, layer L2).

The primer subsystem, extracted from the hub. Holds the PURE, app-free core:

  * Site-directed mutagenesis: the _MUT_* E. coli codon-usage data, the
    biophysics math (_mut_tm / _mut_hairpin_dg / _mut_homodimer_dg / _mut_gc_pct /
    _mut_ends_gc, via primer3 thermodynamics), mutation parse/translate, the
    constant-outer annealing design (_mut_design_outer/_fwd_anneal/_rev_anneal),
    and the record-region INNER design (_mut_design_inner / _mut_design_modified_outer
    / _mut_extract_cds) that picks the mutant codon from the usage map.
  * Binding / display rederivation — the catastrophic core: _rederive_primer_binding
    (longest 3'-anchored match wins; origin-rotation invariant) + _primer_binding_sites
    (origin-wrap-aware site scan) + _primer_check_confidence.
  * Generic primer designers: _primer_tm / _pick_binding_region / _binding_max_len /
    _design_detection_primers / _design_cloning_primers(_raw) / _design_generic_primers.
  * QuikChange-scrub engine [INV-97] (clone-free restriction-site removal): the
    planner (_scrub_design + its _scrub_* helpers) that finds the minimal SILENT
    substitutions destroying each site, and the improved-QuikChange primer designer
    (_scrub_qc_primers / _scrub_qc_verify). Reaches its two hub-pinned deps via
    _state hooks — _all_enzymes_hook (custom-enzyme catalog) and _translate_cds_hook
    (the canonical translator; stays hub-side for its genetic-code-table dep).

Deps are siblings ≤ L2 only: biology L0 (_rc/_circ_slice/_iupac_compatible/_iupac_pattern/
_scan_restriction_sites/_forbidden_hit_set/_mut_revcomp/_slice_circular), util L0
(_normalize_dna_for_align), logging L0, and codon L2 (_codon_build_aa_map — the
mutant-codon pick, the single cross-subsystem dep and why this sibling is L2; codon
does not import primer, so no cycle). NO record / dataaccess / UI coupling — the two
hub-pinned biology deps come through _state hooks, not imports.

CATASTROPHIC-CLASS subsystem (project_primer_design_catastrophic): primers must BIND
+ DISPLAY exactly where they anneal, and a scrub must never silently change a protein
or leave a site behind. Guarded by test_mutagenize + test_scrub +
test_primer_binding_rotation (origin-rotation golden on a real plasmid). What STAYS
hub-side: the grammar-coupled orchestrators (_design_gb_primers, _design_operon_soe_primers),
the Golden-Braid scrub (_scrub_gb_* — needs real digest+ligate assembly = cloning L3),
the primer-usage index + collection cache (_primer_usage_*, _primer_collection_name_taken),
the UI cursor/preview (_mut_next_cursor, _paint_primer_*) — all anchored on the
PlasmidApp God-class, the dataaccess save layer, or an above-L2 sibling. Re-exported by
the hub so sc.<name> + every call site resolves unchanged.
"""
from __future__ import annotations

import re
import threading
from typing import Callable as _Callable

import splicecraft_state as _state
from splicecraft_logging import _log, _timed
from splicecraft_util import (
    _IUPAC_DNA_CHARS, _IUPAC_WEAKEST_BASE, _PRIMER3_LOCK,
    _degenerate_tm_bracket, _normalize_dna_for_align, _p3_thermo, _p3_tm,
)
from splicecraft_biology import (
    _circ_slice, _forbidden_hit_set, _iupac_compatible, _iupac_pattern,
    _mut_revcomp, _rc, _scan_restriction_sites, _search_subsequence,
    _slice_circular,
)
# The mutagenesis-inner design picks the mutant codon via the codon-usage map —
# the one cross-subsystem dep, and why splicecraft_primer sits at L2 (codon is L2;
# codon does not import primer, so no cycle).
from splicecraft_codon import _codon_build_aa_map


_MUT_CODON_USAGE = {
    "GGG": ("G", 44),  "GGA": ("G", 47),  "GGT": ("G", 109), "GGC": ("G", 171),
    "GAG": ("E", 94),  "GAA": ("E", 224), "GAT": ("D", 194), "GAC": ("D", 105),
    "GTG": ("V", 135), "GTA": ("V", 59),  "GTT": ("V", 86),  "GTC": ("V", 60),
    "GCG": ("A", 197), "GCA": ("A", 108), "GCT": ("A", 55),  "GCC": ("A", 162),
    "AGG": ("R", 8),   "AGA": ("R", 7),   "AGT": ("S", 37),  "AGC": ("S", 85),
    "AAG": ("K", 62),  "AAA": ("K", 170), "AAT": ("N", 112), "AAC": ("N", 125),
    "ATG": ("M", 127), "ATA": ("I", 19),  "ATT": ("I", 156), "ATC": ("I", 93),
    "ACG": ("T", 59),  "ACA": ("T", 33),  "ACT": ("T", 41),  "ACC": ("T", 117),
    "TGG": ("W", 55),  "TGT": ("C", 30),  "TGC": ("C", 41),
    "TAT": ("Y", 86),  "TAC": ("Y", 75),
    "TTG": ("L", 61),  "TTA": ("L", 78),  "TTT": ("F", 101), "TTC": ("F", 77),
    "TCG": ("S", 41),  "TCA": ("S", 40),  "TCT": ("S", 29),  "TCC": ("S", 28),
    "CGG": ("R", 21),  "CGA": ("R", 22),  "CGT": ("R", 108), "CGC": ("R", 133),
    "CAG": ("Q", 142), "CAA": ("Q", 62),  "CAT": ("H", 81),  "CAC": ("H", 67),
    "CTG": ("L", 240), "CTA": ("L", 27),  "CTT": ("L", 61),  "CTC": ("L", 54),
    "CCG": ("P", 137), "CCA": ("P", 34),  "CCT": ("P", 43),  "CCC": ("P", 33),
    "TAA": ("*", 9),   "TAG": ("*", 0),   "TGA": ("*", 5),
}


_MUT_CODON_TO_AA = {c: aa for c, (aa, _) in _MUT_CODON_USAGE.items()}


_MUT_STOPS       = {"TAA", "TAG", "TGA"}


def _mut_aa_to_codons() -> dict:
    from collections import defaultdict
    totals = defaultdict(int)
    for c, (aa, n) in _MUT_CODON_USAGE.items():
        totals[aa] += n
    result: dict = defaultdict(list)
    for c, (aa, n) in _MUT_CODON_USAGE.items():
        if aa == "*":
            continue
        result[aa].append((c, n / totals[aa] if totals[aa] else 0.0))
    for aa in result:
        result[aa].sort(key=lambda x: -x[1])
    return dict(result)


_MUT_AA_TO_CODONS = _mut_aa_to_codons()


_MUT_BSAI_FWD_TAIL = "CCCC" + "GGTCTCA" + "AATG"   # 15 nt; AATG = A(extra)+ATG ovhg


_MUT_BSAI_REV_TAIL = "CCCC" + "GGTCTCA" + "AACG"   # 15 nt; AACG = revcomp(CGTT)


# ONE set of reaction conditions for every nearest-neighbour number the app
# reports: primer3's own defaults (50 mM monovalent, 1.5 mM Mg2+, 0.6 mM dNTP,
# 50 nM oligo) — the same conditions `_primer_tm`, `_pick_binding_region` and
# the `.dna` import compute under. Until 2026-09-22 the mutagenesis designers
# used 0.2 mM dNTP / 250 nM oligo instead, so the SAME oligo read 2-3 °C
# higher on the mutagenesis screen than on Check Primer, and a PCR program
# built from a stored mutagenesis Tm annealed that much too hot.
# Typed as a plain `dict` so `**_MUT_P3` type-checks into the primer3 helpers.
_MUT_P3: dict = dict(mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0)


# ── Tm profiles: which reaction conditions a reported Tm is computed under ──
# A Tm means nothing without its conditions, and tools pick different ones: the
# same 25-mer reads 62.4 °C here and 58.2 °C in Benchling. A profile names one
# set; every Tm in a response that took a `tm_profile` is computed under it,
# and the response names the profile back. Both use primer3's SantaLucia 1998
# nearest-neighbour model with the SantaLucia salt correction (primer3's
# defaults) — only the concentrations differ. `primer3_default` IS `_MUT_P3`,
# the app-wide conditions above, so a caller that never names a profile gets
# exactly the numbers it always got. `benchling_compatible` is Benchling's
# default SantaLucia setting — no Mg2+ or dNTP, 5x the oligo — which reads an
# oligo ~4-5 °C LOWER (the missing Mg2+ outweighs the extra oligo); it
# reproduces Tms read off Benchling to 0.1 °C (tests/test_primer_pair.py).
#
# A profile moves the Tm ONLY. Secondary structure (`_primer_structure`) is
# always judged in a real PCR buffer, `_MUT_P3`: leaving Mg2+ out is a way of
# REPORTING a Tm, but every PCR tube has Mg2+ in it, and scoring a hairpin
# without it would understate exactly the risk being checked.
_TM_PROFILES: "dict[str, dict]" = {
    "primer3_default": {
        "label": "Primer3 defaults",
        "conditions": dict(_MUT_P3),
        "about": ("SpliceCraft's own conditions — what every screen, the "
                  "mutagenesis designers and the .dna import show: primer3's "
                  "defaults, a typical PCR buffer."),
    },
    "benchling_compatible": {
        "label": "Benchling-compatible",
        "conditions": dict(mv_conc=50.0, dv_conc=0.0, dntp_conc=0.0,
                           dna_conc=250.0),
        "about": ("Benchling's default SantaLucia melting-temperature settings "
                  "(no Mg2+ or dNTP, 250 nM oligo), for comparing a Tm with "
                  "one Benchling shows. If you changed Benchling's Melting "
                  "Temp settings, compare the conditions instead."),
    },
}
_TM_PROFILE_DEFAULT = "primer3_default"
_TM_PROFILE_METHOD = ("SantaLucia 1998 nearest-neighbour, SantaLucia salt "
                      "correction (primer3 calc_tm)")
# Spellings a caller might reasonably send for one of the profiles above
# (after case-folding and `-` / space → `_`).
_TM_PROFILE_ALIASES = {
    "default": "primer3_default",
    "primer3": "primer3_default",
    "benchling": "benchling_compatible",
}


def _tm_profile_resolve(name: object) -> "str | None":
    """The canonical profile key for `name`, or None when it names none.
    None / blank → the default profile; case, `-` and spaces are forgiven and
    `_TM_PROFILE_ALIASES` apply, so `"Benchling"` finds `benchling_compatible`.
    Anything that is not a string is unknown (None) — never the default."""
    if name is None:
        return _TM_PROFILE_DEFAULT
    if not isinstance(name, str):
        return None
    key = "_".join(name.strip().lower().replace("-", " ").split())
    if not key:
        return _TM_PROFILE_DEFAULT
    if key in _TM_PROFILES:
        return key
    return _TM_PROFILE_ALIASES.get(key)


def _tm_profile_kwargs(name: str = _TM_PROFILE_DEFAULT) -> dict:
    """The primer3 `calc_tm` keywords (salt / dNTP / oligo concentrations) of
    a Tm profile — a copy, so a caller can't edit the registry. Raises
    ValueError for a name `_tm_profile_resolve` doesn't recognise."""
    key = _tm_profile_resolve(name)
    if key is None:
        raise ValueError(
            f"unknown tm_profile {name!r}; expected one of "
            f"{', '.join(_TM_PROFILES)}")
    return dict(_TM_PROFILES[key]["conditions"])


def _tm_profile_catalog() -> "list[dict]":
    """Every Tm profile with its conditions spelled out in lab units — what
    `list-tm-profiles` returns, so a caller can see what a number came from."""
    out = []
    for key, prof in _TM_PROFILES.items():
        c = prof["conditions"]
        out.append({
            "key": key,
            "label": prof["label"],
            "default": key == _TM_PROFILE_DEFAULT,
            "monovalent_mm": c["mv_conc"],
            "mg_mm": c["dv_conc"],
            "dntp_mm": c["dntp_conc"],
            "oligo_nm": c["dna_conc"],
            "method": _TM_PROFILE_METHOD,
            "about": prof["about"],
        })
    return out


def _mut_parse(s: str) -> tuple:
    """Parse a mutation string like 'W140F'. Returns (wt_aa, pos_1based, mut_aa)."""
    m = re.fullmatch(r"([A-Za-z\*])(\d+)([A-Za-z\*])", s.strip())
    if not m:
        raise ValueError(f"Cannot parse '{s}'. Use format: [WT][pos][MUT], e.g. W140F")
    return m.group(1).upper(), int(m.group(2)), m.group(3).upper()


def _mut_translate(dna: str) -> str:
    aa: list = []
    for i in range(0, len(dna) - 2, 3):
        c = dna[i:i+3].upper()
        if c in _MUT_STOPS:
            break
        aa.append(_MUT_CODON_TO_AA.get(c, "?"))
    return "".join(aa)


# Thermodynamic results are pure functions of the sequence (the only other
# input, `_MUT_P3`, is a frozen module constant), but the inner-design grid
# (up to 23×23 = 529 windows) and the QuikChange scrub (10×9×9 = 810 pairs)
# call them ~1,500× per design over heavily-overlapping sequences. Memoize the
# primer3 SUCCESS path only — NEVER the approximation fallback — so (a) a cached
# real value can never shadow the fallback the property tests force by making
# primer3 unimportable, and (b) the cache stays a pure primer3 mirror. Bounded
# (str→float entries are tiny) with oldest-insertion eviction; clearable for
# test isolation. NOT re-exported, so the public surface is unchanged.
_MUT_TM_CACHE: "dict[str, float]" = {}
_MUT_HAIRPIN_CACHE: "dict[str, float]" = {}
_MUT_HOMODIMER_CACHE: "dict[str, float]" = {}
_MUT_THERMO_CACHE_MAX = 8192
# The three caches above are read/filled from BOTH the agent-API design
# endpoints (design-primers / design-mutagenesis / simulate-pcr / check-primer,
# on ThreadingMixIn worker threads) AND the Textual design workers. Guard the
# evict+insert so two concurrent designs can't `pop(next(iter()))` the same key
# (KeyError) or iterate a dict mid-insert (RuntimeError) — 2026-07 finding.
_MUT_THERMO_CACHE_LOCK = threading.Lock()


def _mut_thermo_cache_put(cache: "dict[str, float]", seq: str, val: float) -> None:
    with _MUT_THERMO_CACHE_LOCK:
        if len(cache) >= _MUT_THERMO_CACHE_MAX:
            cache.pop(next(iter(cache)))       # evict oldest insertion (FIFO)
        cache[seq] = val


def _mut_thermo_cache_clear() -> None:
    """Drop all memoized primer3 thermodynamic results. The fallback property
    tests call this so a real value cached by an earlier primer-design test
    can't mask the primer3-unavailable path they deliberately force."""
    with _MUT_THERMO_CACHE_LOCK:
        _MUT_TM_CACHE.clear()
        _MUT_HAIRPIN_CACHE.clear()
        _MUT_HOMODIMER_CACHE.clear()


def _mut_tm(seq: str) -> float:
    hit = _MUT_TM_CACHE.get(seq)
    if hit is not None:
        return hit
    try:
        val = _p3_tm(seq, **_MUT_P3)
    except Exception:
        # A degenerate window (a consensus N, a typed IUPAC code) is a MIX:
        # rate it by its weakest member's nearest-neighbour Tm — the same
        # scale as every clean window — before the crude rule below.
        bracket = _degenerate_tm_bracket(seq, **_MUT_P3)
        if bracket is not None:
            _mut_thermo_cache_put(_MUT_TM_CACHE, seq, bracket[0])
            return bracket[0]
        # Fall back to the crude 2×AT + 4×GC approximation when
        # primer3 is missing or raises (degenerate input, NaN config).
        # Log so a wave of failures shows up as one diagnosable
        # symptom in the bug-report bundle instead of silent
        # mis-temperature on every primer. NOT cached (see cache note).
        _log.exception(
            "_mut_tm: primer3.calc_tm fell back to GC approximation "
            "for %d-mer", len(seq))
        su = seq.upper()
        gc = sum(1 for c in su if c in "GC")
        at = sum(1 for c in su if c in "AT")
        # Degenerate/IUPAC bases (N, R, Y, …) contribute to neither term; count
        # them as the AT/GC midpoint (3) so a user-pasted degenerate oligo isn't
        # Tm-UNDERestimated (which would over-extend a binding-region search).
        other = len(su) - gc - at
        return 2 * at + 4 * gc + 3 * other
    _mut_thermo_cache_put(_MUT_TM_CACHE, seq, val)
    return val


# primer3-py reports every ΔG in CALORIES per mole. Oligo literature, every
# supplier's spec sheet and every QC rule of thumb ("reject below −9 kcal/mol")
# are in KILOcalories, so the raw number is 1000× what a reader expects: a
# perfectly ordinary hairpin came back as `-738.3` and an ordinary self-dimer
# as `-4993.7`, which is unusable as a gate and reads as catastrophic if you
# take the unit at face value. Convert once, here, so every caller and every
# reported field is in kcal/mol (2026-08-26 agent field report #3).
_CAL_PER_KCAL = 1000.0


def _mut_hairpin_dg(seq: str) -> "float | None":
    """Hairpin ΔG in **kcal/mol** (negative = more stable), or ``None`` when
    it could not be computed.

    ``None`` rather than ``0.0``: zero is the BEST possible value, so a
    missing primer3 would make every QC gate pass on a primer nobody
    measured. Callers that only rank treat it as no penalty; callers that
    REPORT it must pass the None through."""
    hit = _MUT_HAIRPIN_CACHE.get(seq)
    if hit is not None:
        return hit
    try:
        val = _p3_thermo("hairpin", seq, **_MUT_P3).dg / _CAL_PER_KCAL
    except Exception:
        _log.exception(
            "_mut_hairpin_dg: primer3.calc_hairpin raised on %d-mer; "
            "reporting None (not measured)", len(seq))
        return None
    _mut_thermo_cache_put(_MUT_HAIRPIN_CACHE, seq, val)
    return val


def _mut_homodimer_dg(seq: str) -> "float | None":
    """Self-dimer ΔG in **kcal/mol**, or ``None`` when it could not be
    computed — see `_mut_hairpin_dg` for why not 0.0."""
    hit = _MUT_HOMODIMER_CACHE.get(seq)
    if hit is not None:
        return hit
    try:
        val = _p3_thermo("homodimer", seq, **_MUT_P3).dg / _CAL_PER_KCAL
    except Exception:
        _log.exception(
            "_mut_homodimer_dg: primer3.calc_homodimer raised on "
            "%d-mer; reporting None (not measured)", len(seq))
        return None
    _mut_thermo_cache_put(_MUT_HOMODIMER_CACHE, seq, val)
    return val


def _mut_gc_pct(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) * 100 if seq else 0.0


def _mut_ends_gc(seq: str) -> bool:
    return bool(seq) and seq[-1].upper() in "GC"


def _mut_score_outer(anneal: str, target_tm: float = 60.0,
                     tm_fn: "_Callable[[str], float] | None" = None) -> float:
    # `tm_fn` lets a caller rank under another Tm profile; left out, it is
    # `_mut_tm` (looked up at CALL time, so a test patching it still counts).
    t  = (tm_fn or _mut_tm)(anneal)
    gc = _mut_gc_pct(anneal)
    # ΔG is kcal/mol; an unmeasurable one scores as no penalty (ranking only —
    # the REPORTED field keeps the None so a QC gate can't pass vacuously).
    hp = _mut_hairpin_dg(anneal) or 0.0
    return (
        abs(t - target_tm) * 2.0
        + (0 if _mut_ends_gc(anneal) else 4.0)
        + max(0.0, -hp - 1.0) / 0.4
        + abs(gc - 50) * 0.1
    )


def _mut_design_fwd_anneal(dna: str) -> "dict | None":
    body = dna[3:]
    best = None
    for length in range(18, 28):
        anneal = body[:length]
        if len(anneal) < 18:
            continue
        s = _mut_score_outer(anneal)
        if best is None or s < best["score"]:
            best = {
                "anneal": anneal,
                "full":   _MUT_BSAI_FWD_TAIL + anneal,
                "tm_anneal": _mut_tm(anneal),
                "gc":     _mut_gc_pct(anneal),
                "score":  s,
            }
    return best


def _mut_design_rev_anneal(dna: str) -> "dict | None":
    end_rc = _mut_revcomp(dna)
    best = None
    for length in range(18, 28):
        anneal = end_rc[:length]
        if len(anneal) < 18:
            continue
        s = _mut_score_outer(anneal)
        if best is None or s < best["score"]:
            best = {
                "anneal": anneal,
                "full":   _MUT_BSAI_REV_TAIL + anneal,
                "tm_anneal": _mut_tm(anneal),
                "gc":     _mut_gc_pct(anneal),
                "score":  s,
            }
    return best


def _mut_design_outer(dna: str) -> dict:
    """Constant FWD/REV outer primers with BsaI-AATG / BsaI-AACG tails."""
    fwd = _mut_design_fwd_anneal(dna)
    rev = _mut_design_rev_anneal(dna)
    if fwd is None or rev is None:
        raise RuntimeError("CDS is too short to design outer primers (need ≥ 21 nt).")
    return {
        "fwd": fwd, "rev": rev,
        "b3_overhang": "AATG",
        "b5_overhang": "CGTT",
        "fwd_anneal_start": 3,
    }


# ── circular primer-binding re-derivation (Phase D, moved from hub) ─────────
# THE catastrophic 'where does a primer ACTUALLY land on the (circular) map'
# core: longest 3'-anchored template match, flap excluded, wrap-aware,
# hint-tie-broken. Verified by the real-plasmid golden (site-containment +
# origin-rotation invariance + multi-site hint). pos_end==total => ends at origin.
# Shortest contiguous match we'll trust as a genuine primer binding when
# re-deriving from the template (Domesticator binding is ≥18 bp; 12 is a
# safe floor that still rejects spurious short coincidences).
_PRIMER_REBIND_MIN: int = 12


def _rederive_primer_binding(primer_seq: str, strand: int, template: str,
                             total: int, hint_start: int = 0,
                             *, circular: bool = True,
                             ) -> "tuple[int, int] | None":
    """Find where a primer's annealing region ACTUALLY binds the
    (circular) ``template``, so a stale / mis-saved ``pos_start`` /
    ``pos_end`` can't park the primer off its true site on the map.

    Returns ``(pos_start, pos_end)`` in top-strand coordinates —
    half-open, with ``pos_end < pos_start`` when the binding wraps the
    origin (and ``pos_end == total`` when it ends exactly at the
    origin). Returns ``None`` when no clean binding is found, so the
    caller keeps the stored positions.

    The binding region is the LONGEST contiguous stretch at the primer's
    3' end that matches the template — forward: the primer's own 3'
    suffix on the top strand; reverse: that suffix's reverse-complement
    on the top strand (which equals a prefix of ``rc(primer)``). The 5'
    flap (enzyme site / spacer / fusion overhang) is whatever doesn't
    match. When a primer legitimately binds more than one site, the
    occurrence nearest ``hint_start`` (the stored position) wins, so the
    primer still lands where the user designed it."""
    seq = (primer_seq or "").upper()
    if not seq or not template or total <= 0:
        return None
    # Circular search space: append the wrap-around head so a binding
    # that spans the origin is found as one contiguous slice. Cap the
    # tail at the primer length so we never scan more than necessary.
    # A LINEAR template's ends don't join, so a primer can't anneal
    # across them (real-world affinity) — search the bare template.
    if circular:
        tail = template[:max(0, min(len(seq), total) - 1)]
        aug = template + tail
    else:
        aug = template
    rc = _rc(seq) if strand < 0 else ""
    max_L = min(len(seq), total)
    hint = (int(hint_start) % total) if total else 0
    for L in range(max_L, _PRIMER_REBIND_MIN - 1, -1):
        target = (seq[len(seq) - L:] if strand >= 0 else rc[:L])
        starts: list[int] = []
        i = aug.find(target)
        while i != -1 and i < total:
            starts.append(i)
            i = aug.find(target, i + 1)
        if not starts:
            continue
        # Closest occurrence to the stored hint (circular distance).
        def _cdist(p: int) -> int:
            d = abs(p - hint)
            return min(d, total - d)
        m = min(starts, key=_cdist)
        pos_end = m + L
        if pos_end > total:          # wraps the origin
            pos_end -= total
        return (m, pos_end)
    # No clean contiguous suffix — keep the stored positions (caller decides).
    # NB: we deliberately do NOT slide the whole primer to a best-offset here.
    # `_attach_pcr_primers_to_record` calls with hint_start=0, so a best-offset
    # window would sit at the origin and could anchor a non-binding primer
    # there — a mis-placed cloning primer is catastrophic. The length-short
    # *stored*-feature repair (a fragment built 1 bp short of its primer) is
    # handled in `PlasmidMap._parse`, where start/end/length are all known.
    return None


# ── primer-CHECK binding finder (Phase D, moved from hub) ──────────────────
# Mismatch-tolerant 3'-anchored binding-site list (which sites a primer binds, with
# identity %) + the confidence glyph. Verified by the real-plasmid golden
# (site lists + origin-rotation invariance). Uses the biology align primitives.
_PRIMER_CHECK_SEED_LEN  = 12     # exact 3'-anchor required for a binding call


_PRIMER_CHECK_MAX_SITES = 200    # per-template binding-site cap (repeat guard)


# Shortest primer the Primer Check tab and `check-primer-pair` accept (cf the
# PCR simulator's `_PCR_MIN_PRIMER_LEN`): below it a 3'-anchored "site" is
# noise and primer3's Tm is meaningless (a 1-mer reads −1,000,000 °C).
_PRIMER_CHECK_MIN_LEN   = 10


def _primer_binding_sites(
    primer: str, top: str, total: int, *,
    circular: bool = True,
    seed_len: int = _PRIMER_CHECK_SEED_LEN,
    min_identity_pct: float = 0.0,
    max_sites: int = _PRIMER_CHECK_MAX_SITES,
) -> "list[dict]":
    """3'-anchored binding sites of `primer` on the top strand `top`
    (length `total`, ASSUMED pre-normalised — uppercase IUPAC, no whitespace).

    A site requires an EXACT match over the primer's 3'-terminal `seed_len`
    bases (clamped to the primer length); identity is then computed over the
    FULL primer. Returns sites sorted best-first (identity desc)::

        {"strand": +1 | -1,   # +1 forward (primer == top-strand sense),
                              # -1 reverse (primer anneals TO the top strand)
         "foot_start": int,   # 0-based footprint start on the top strand,
                              # canonical [0, total); footprint spans
                              # [foot_start, foot_start+length) around the circle
         "length": int,       # primer length
         "ident_pct": float,  # full-primer identity 0..100
         "mismatches": int}

    The primer's 3' end is the HIGH-coord edge of a forward footprint and the
    foot_start (LOW-coord) edge of a reverse footprint — so an amplicon runs
    from a forward site's foot_start to a reverse site's foot_start+length.

    Raises ValueError (via `_normalize_dna_for_align`) on a foreign character
    in `primer`.
    """
    P = _normalize_dna_for_align(primer or "")
    L = len(P)
    if L == 0 or total <= 0 or L > total:
        return []
    seed = max(1, min(int(seed_len), L))
    anchor = P[-seed:]
    try:
        # Exact 3'-anchor hits, wrap-aware, via the tested matcher — as TWO
        # single-strand searches, one per primer role: `anchor` on the top
        # strand is the 3' end of a FORWARD-role primer, rc(anchor) on the top
        # strand the 3' end of a REVERSE-role primer. One both-strand search
        # is NOT equivalent: it skips the reverse pass for a palindromic query
        # and collapses a '+' and '-' hit on one span into the '+', so a primer
        # whose 3' seed is palindromic (a restriction site at the 3' end — the
        # normal shape of a site-tailed primer with a short seed) never had its
        # reverse-role site evaluated: a perfect reverse primer came back as a
        # forward site 13 bp upstream at ~56 % identity.
        fwd_hits = _search_subsequence(
            top, anchor, max_mismatches=0,
            circular=circular, both_strands=False,
        )
        rev_hits = _search_subsequence(
            top, _rc(anchor), max_mismatches=0,
            circular=circular, both_strands=False,
        )
    except ValueError:
        return []
    iupac_ok = _iupac_compatible
    sites: "list[dict]" = []
    seen: "set[tuple[int, int]]" = set()
    roles = [(h, 1) for h in fwd_hits] + [(h, -1) for h in rev_hits]
    for h, s_strand in roles:
        hs, he = h["start"], h["end"]
        if s_strand == 1:
            foot_start = he - L           # 5' edge of the forward footprint
        else:
            foot_start = hs               # 3'/left edge of the reverse footprint
        if circular:
            window = _circ_slice(top, foot_start, L, total)
        else:
            if foot_start < 0 or foot_start + L > total:
                continue                  # primer hangs off a linear end
            window = top[foot_start:foot_start + L]
        if len(window) != L:
            continue
        oriented = window if s_strand == 1 else _rc(window)
        mm = 0
        for i in range(L):
            if not iupac_ok(P[i], oriented[i]):
                mm += 1
        ident = 100.0 * (L - mm) / L
        if ident < min_identity_pct:
            continue
        canon = foot_start % total
        key = (canon, s_strand)
        if key in seen:
            continue
        seen.add(key)
        sites.append({
            "strand":     s_strand,
            "foot_start": canon,
            "length":     L,
            "ident_pct":  ident,
            "mismatches": mm,
        })
    # Rank THEN cap. Capping first kept the first `max_sites` sites in scan
    # order and sorted those, so on a repeat-rich template the best site could
    # be dropped and `best_identity` / "binds as designed" answered from the
    # leftovers. A caller that needs to know about truncation asks for one
    # more than it shows.
    sites.sort(key=lambda s: (-s["ident_pct"], s["foot_start"], s["strand"]))
    return sites[:max(0, int(max_sites))]


def _primer_check_confidence(pct: "float | int | None") -> "tuple[str, str]":
    """Map a primer-binding / amplicon identity to a (glyph, Rich-colour)
    confidence badge for the Primer Check results table. Mirrors the
    alignment-status tiers so ✓/⚠/~/✗ read consistently app-wide."""
    if not isinstance(pct, (int, float)):
        return ("?", "white")
    v = float(pct)
    if v >= 99.999:
        return ("✓", "bright_cyan")
    if v >= 90.0:
        return ("✓", "green")
    if v >= 75.0:
        return ("⚠", "yellow")
    if v >= 60.0:
        return ("~", "dark_orange")
    return ("✗", "red")


# ── generic primer design (Phase D, moved from hub) ────────────────────────
# Tm + binding-region selection + the cloning / detection / generic primer
# designers (primer3 thermodynamics; enzyme catalog via _state._all_enzymes_hook).
# Verified by the real-plasmid design golden (byte-identical output). The GB /
# domestication-scrub designers (_design_gb_primers / _scrub_*) stay hub-side.
def _primer_tm(seq: str,
               tm_profile: str = _TM_PROFILE_DEFAULT) -> "float | None":
    """Melting temperature (°C, 1 dp) of an oligo under a Tm profile
    (`_TM_PROFILES`; the default is the app-wide conditions) — primer3's
    nearest-neighbour model when available, else the 2(A+T)+4(G+C) rule.
    Module-level so the CSV import (and any caller) can compute a Tm without
    the local ``_calc_tm`` closures the design / .dna paths use. Returns None
    for empty input; raises ValueError for an unknown profile."""
    s = (seq or "").strip().upper()
    if not s:
        return None
    conditions = _tm_profile_kwargs(tm_profile)
    try:
        return round(_p3_tm(s, **conditions), 1)
    except Exception:
        # Degenerate oligo → its weakest member's NN Tm (what the anneal must
        # respect), not the 2+4 rule, which runs ~8-12 °C away from it.
        bracket = _degenerate_tm_bracket(s, **conditions)
        if bracket is not None:
            return round(bracket[0], 1)
        gc = sum(1 for c in s if c in "GC")
        at = sum(1 for c in s if c in "AT")
        # Degenerate/IUPAC bases → AT/GC midpoint (3) so a degenerate oligo
        # isn't Tm-underestimated (mirrors `_mut_tm`).
        other = len(s) - gc - at
        return float(2 * at + 4 * gc + 3 * other)


# Hard cap on the TOTAL synthesised oligo length: the 5' tail (pad + enzyme
# site + spacer + overhang) PLUS the 3' binding region. Oligo-synthesis cost
# is driven by the whole oligo and standard (cheap) synthesis tops out around
# here, so each design grows its binding region to reach `target_tm` but
# never past `_PRIMER_MAX_OLIGO_LEN − len(tail)` (2026-06-09, user spec).
# For a low-GC (AT-rich — e.g. codon-optimised for a low-GC host like
# E. faecium) part this lets the binding extend well past the old fixed
# 25 nt so it can actually reach ~60 °C; for a high-GC part the "closest to
# target" pick stays short. If even the capped binding can't reach
# `target_tm`, the closest is returned + flagged (`_GB_TM_OFFTARGET_MARGIN`)
# rather than bloating the oligo past the synthesis budget.
_PRIMER_MAX_OLIGO_LEN = 50


def _binding_max_len(tail_len: int, min_len: int = 18) -> int:
    """Largest binding-region length that keeps the TOTAL oligo
    (`tail_len` + binding) within `_PRIMER_MAX_OLIGO_LEN`. Never returns
    below `min_len`: an over-long tail can't shrink the binding below the
    minimum usable annealing length — that (rare) design is surfaced by the
    downstream low-Tm advisory instead of an unusably short arm."""
    return max(min_len, _PRIMER_MAX_OLIGO_LEN - max(0, int(tail_len)))


def _pick_binding_region(seq: str, target_tm: float = 60.0,
                         min_len: int = 18, max_len: int = 25,
                         tm_profile: str = _TM_PROFILE_DEFAULT,
                         ) -> tuple[str, float]:
    """Return the prefix of `seq` (length min_len..max_len) whose Tm is
    closest to `target_tm`. Uses primer3-py's SantaLucia Tm calculation under
    the Tm profile's conditions (`_TM_PROFILES`), so `target_tm` means the same
    thing the profile's reported Tm does.

    Returns (binding_sequence, tm). If primer3-py is not installed, falls
    back to a crude 2+4 rule estimate.
    """
    # Type the dispatcher as a `(str) -> float` Callable so pyright can
    # accept both primer3.calc_tm (which has a richly-typed signature
    # with extra defaulted kwargs) and the fallback approximation
    # below. Using two separate names avoids the param-name mismatch
    # pyright flags when a `def _tm(s)` re-defines the same binding.
    def _tm_fallback(s: str) -> float:
        u = s.upper()
        gc = sum(1 for c in u if c in "GC")
        at = sum(1 for c in u if c in "AT")
        # Degenerate / ambiguous bases score the AT–GC midpoint so an
        # N-containing arm isn't Tm-UNDERestimated into looking unusable
        # (mirrors `_primer_tm` / `_mut_tm`).
        return float(2 * at + 4 * gc + 3 * (len(u) - gc - at))
    _tm: "_Callable[..., float]"
    conditions = _tm_profile_kwargs(tm_profile)
    try:
        import primer3  # noqa: F401 — availability probe; `_p3_tm` computes
    except ImportError:
        _tm = _tm_fallback
    else:
        # primer3's nearest-neighbour model REFUSES any non-ACGT base
        # (`ValueError: Sequence contains non-ACGT base 'N' at position …`),
        # and real templates carry them: a Plasmidsaurus consensus has N at
        # low-coverage positions, NCBI records carry IUPAC ambiguity codes,
        # and the Synthesis composer lets you TYPE IUPAC ([INV-53]). Pre-fix
        # only `ImportError` was caught, so designing primers over such a
        # region raised straight out of the worker instead of degrading to
        # the approximation — every sibling Tm helper (`_primer_tm`,
        # `_mut_tm`, `_tm_cached`, the `.dna` import `_calc_tm`) already
        # falls back per-call; this one didn't.
        #
        # Decide ONCE for the whole candidate window rather than per length:
        # mixing a nearest-neighbour Tm for the clean prefixes with a 2+4
        # estimate for the ones that reach the ambiguous base would compare
        # two different scales and pick the length whose ESTIMATOR happened
        # to land nearest `target_tm`.
        # NB distinct name (see the two-separate-names note above): a
        # `def _tm(...)` here would re-declare the annotated binding and
        # pyright's `reportRedeclaration` fails the release lint gate.
        # A degenerate prefix is rated by its WEAKEST member's NN Tm
        # (ambiguity codes resolved to A/T where they allow it) — still the
        # nearest-neighbour scale, so a window that reaches a consensus N no
        # longer drops the whole design onto the 2+4 rule.
        def _tm_nn(s: str) -> float:
            try:
                return _p3_tm(s.upper().translate(_IUPAC_WEAKEST_BASE),
                              **conditions)
            except (ImportError, ValueError, OSError, RuntimeError,
                    TypeError):
                return _tm_fallback(s)
        _window = seq[:max_len].upper()
        _tm = (_tm_nn
               if (_window and set(_window) <= _IUPAC_DNA_CHARS)
               else _tm_fallback)

    # Defensive init: if the caller forgot the len(seq) >= min_len guard,
    # the loop below won't execute and we'd otherwise return Tm=0 with a
    # too-short binding. Compute Tm for whatever is there so downstream
    # validation (low Tm, short primer) still trips honestly.
    best_seq = seq[:max(min_len, 1)]
    best_tm  = _tm(best_seq) if best_seq else 0.0
    best_diff = float("inf")
    # Pick the length whose Tm is CLOSEST to target — minimising
    # |tm - target|. When `target_tm` is unreachable within the
    # [min_len, max_len] window (too high for an AT-rich arm even at
    # max_len, too low for a GC-rich arm even at min_len), the nearest
    # achievable Tm wins — the longest candidate when the target sits
    # above the whole window, the shortest when below — so the caller
    # always gets the next-best binding instead of a failure.
    for n in range(min_len, min(max_len + 1, len(seq) + 1)):
        candidate = seq[:n]
        tm = _tm(candidate)
        diff = abs(tm - target_tm)
        if diff < best_diff:
            best_seq, best_tm, best_diff = candidate, tm, diff
    return best_seq, best_tm


@_timed("op.primer3.detection_design")
def _design_detection_primers(
    template_seq: str,
    target_start: int,
    target_end: int,
    product_min: int = 450,
    product_max: int = 550,
    target_tm: float = 60.0,
    primer_len: int = 25,
    tm_profile: str = _TM_PROFILE_DEFAULT,
) -> dict:
    """Design diagnostic PCR primers WITHIN a selected region using Primer3.

    Both primers bind INSIDE the region (target_start..target_end) and the
    amplicon is product_min..product_max bp. This is the standard approach
    for detection/screening primers: you pick a gene or feature and want a
    ~500 bp diagnostic band from within it.

    Uses SEQUENCE_INCLUDED_REGION (not SEQUENCE_TARGET) so Primer3 places
    both primers inside the selected region rather than trying to flank it.

    `tm_profile` (`_TM_PROFILES`) sets Primer3's salt / dNTP / oligo
    concentrations, so `target_tm` and the returned Tms are on that profile's
    scale — the SAME numbers `_primer_tm(seq, tm_profile)` gives. Primer3 has
    one set of conditions per run, so its own hairpin / dimer screening runs
    under the profile too; `_primer_pair_qc` re-judges those in PCR buffer.

    Returns a dict with keys: fwd_seq, rev_seq, fwd_tm, rev_tm, fwd_pos,
    rev_pos, product_size, or an 'error' key on failure.
    """
    import primer3
    cond = _tm_profile_kwargs(tm_profile)
    seq   = template_seq.upper()
    total = len(seq)
    wraps = target_end < target_start

    # Primer3 is linear-only. For a wrap region we rotate the template
    # so the region becomes contiguous at [0, region_len), run Primer3,
    # then unrotate the returned positions via (coord + rotation) % total.
    if wraps:
        rotation    = target_start
        p3_seq      = seq[target_start:] + seq[:target_start]
        region_len  = (total - target_start) + target_end
        p3_start    = 0
    else:
        rotation    = 0
        p3_seq      = seq
        region_len  = target_end - target_start
        p3_start    = target_start

    if region_len < 1:
        return {"error": "Target region is empty."}
    if region_len < product_min:
        return {
            "error": f"Region ({region_len} bp) is shorter than minimum "
                     f"product size ({product_min} bp). Select a larger "
                     f"region or reduce the product size."
        }

    seq_args = {
        "SEQUENCE_TEMPLATE": p3_seq,
        # INCLUDED_REGION: primers must bind WITHIN this region.
        # This is the key difference from SEQUENCE_TARGET (which
        # would require primers to sit OUTSIDE the target).
        "SEQUENCE_INCLUDED_REGION": [p3_start, region_len],
    }
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        # primer_len is the OPTIMAL length — Primer3 will expand
        # or contract within the min/max range to find the best Tm.
        "PRIMER_OPT_SIZE": primer_len,
        "PRIMER_MIN_SIZE": max(15, primer_len - 8),
        "PRIMER_MAX_SIZE": min(36, primer_len + 8),
        "PRIMER_OPT_TM": target_tm,
        "PRIMER_MIN_TM": target_tm - 3,
        "PRIMER_MAX_TM": target_tm + 3,
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
        "PRIMER_NUM_RETURN": 1,
        # The Tm profile's reaction conditions. The default profile's ARE
        # Primer3's defaults — passed anyway so the two can't drift apart.
        "PRIMER_SALT_MONOVALENT": cond["mv_conc"],
        "PRIMER_SALT_DIVALENT": cond["dv_conc"],
        "PRIMER_DNTP_CONC": cond["dntp_conc"],
        "PRIMER_DNA_CONC": cond["dna_conc"],
    }
    try:
        # Primer3's design runs its thermodynamic alignment on the same C
        # statics as `_p3_thermo` — serialise it with every other alignment.
        with _PRIMER3_LOCK:
            result = primer3.design_primers(seq_args=seq_args,
                                            global_args=global_args)
    except (OSError, RuntimeError, ValueError, KeyError, TypeError) as exc:
        # Sweep #25 (2026-05-23): `(OSError, Exception)` is `Exception`
        # since `Exception` subsumes `OSError` — that tuple was a
        # bare-except in disguise (the INV-65 grep for `(AttributeError,
        # Exception)` missed the `(OSError, Exception)` shape).
        # Narrowed to the actual Primer3 failure modes: missing lib
        # (OSError), bad params (ValueError), and the misc errors the
        # C wrapper raises (RuntimeError / KeyError / TypeError when
        # the result dict shape doesn't match).
        return {"error": f"Primer3 rejected parameters: {exc}"}

    n_found = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    if n_found == 0:
        explain = result.get("PRIMER_LEFT_EXPLAIN", "")
        return {"error": f"Primer3 found no valid pair. {explain}"}

    fwd_pos = result["PRIMER_LEFT_0"]     # (start, length) on p3_seq
    rev_pos = result["PRIMER_RIGHT_0"]    # (start, length) — start is 3' end on p3_seq

    # Unrotate positions back to original-template coordinates. Apply
    # `% total` ONLY when the target region wraps origin (Primer3 ran
    # on a rotated template). On the linear path, `rotation == 0` and
    # the modulo is a no-op for most values BUT silently flips a primer
    # 3'-ending at exact bp `total - 1` into `rev_end = total % total = 0`
    # — which `_add_selected_to_map` then reads as a wrap-encoded primer
    # (`p_end < p_start`) and stamps a wrap CompoundLocation with a 0-bp
    # head part. Mirrors `_design_cloning_primers_raw` (gate-on-wraps).
    if wraps:
        fwd_start = (fwd_pos[0] + rotation) % total
        fwd_end   = (fwd_pos[0] + fwd_pos[1] + rotation) % total
        rev_start = (rev_pos[0] - rev_pos[1] + 1 + rotation) % total
        rev_end   = (rev_pos[0] + 1 + rotation) % total
    else:
        fwd_start = fwd_pos[0]
        fwd_end   = fwd_pos[0] + fwd_pos[1]
        rev_start = rev_pos[0] - rev_pos[1] + 1
        rev_end   = rev_pos[0] + 1

    return {
        "fwd_seq":      result["PRIMER_LEFT_0_SEQUENCE"],
        "rev_seq":      result["PRIMER_RIGHT_0_SEQUENCE"],
        "fwd_tm":       round(result["PRIMER_LEFT_0_TM"], 1),
        "rev_tm":       round(result["PRIMER_RIGHT_0_TM"], 1),
        "fwd_pos":      (fwd_start, fwd_end),
        "rev_pos":      (rev_start, rev_end),
        "product_size": result["PRIMER_PAIR_0_PRODUCT_SIZE"],
    }


@_timed("op.primer3.cloning_design")
def _design_cloning_primers_raw(
    template_seq: str,
    start: int,
    end: int,
    site_5: str,
    site_3: str,
    name_5: str = "5'site",
    name_3: str = "3'site",
    target_tm: float = 60.0,
    padding: str = "GCGC",
    tm_profile: str = _TM_PROFILE_DEFAULT,
) -> dict:
    """Design cloning primers with arbitrary recognition-site tails + padding.

    Accepts raw site sequences (not just NEB enzyme names) so users can
    enter custom cutter sequences.

    Structure (5'→3'):
        Forward: [padding] [5' site]    [binding region →]
        Reverse: [padding] [RC 3' site] [← binding region RC]

    Returns dict with keys: fwd_full, rev_full, fwd_binding, rev_binding,
    fwd_tm, rev_tm, re_5prime, re_3prime, site_5, site_3, insert_seq,
    fwd_pos, rev_pos, or 'error'.
    """
    site_5 = site_5.upper()
    site_3 = site_3.upper()
    if not site_5 or not set(site_5) <= set("ACGTRYWSMKBDHVN"):
        return {"error": f"Invalid 5' site sequence: {site_5!r}"}
    if not site_3 or not set(site_3) <= set("ACGTRYWSMKBDHVN"):
        return {"error": f"Invalid 3' site sequence: {site_3!r}"}

    total  = len(template_seq)
    insert = _slice_circular(template_seq.upper(), start, end)
    wraps  = end < start
    if len(insert) < 18:
        return {"error": "Region too short (< 18 bp)."}

    # Cap each binding so the TOTAL oligo (tail + binding) stays within
    # `_PRIMER_MAX_OLIGO_LEN`; tail = padding + RE site.
    fwd_bind, fwd_tm = _pick_binding_region(
        insert, target_tm, max_len=_binding_max_len(len(padding) + len(site_5)),
        tm_profile=tm_profile)
    rev_bind, rev_tm = _pick_binding_region(
        _rc(insert), target_tm,
        max_len=_binding_max_len(len(padding) + len(site_3)),
        tm_profile=tm_profile)

    fwd_full = padding + site_5 + fwd_bind
    rev_full = padding + _rc(site_3) + rev_bind

    if wraps:
        fwd_pos = (start, (start + len(fwd_bind)) % total)
        rev_pos = ((end - len(rev_bind)) % total, end)
    else:
        fwd_pos = (start, start + len(fwd_bind))
        rev_pos = (end - len(rev_bind), end)

    return {
        "fwd_full":    fwd_full,
        "rev_full":    rev_full,
        "fwd_binding": fwd_bind,
        "rev_binding": rev_bind,
        "fwd_tm":      round(fwd_tm, 1),
        "rev_tm":      round(rev_tm, 1),
        "re_5prime":   name_5,
        "re_3prime":   name_3,
        "site_5":      site_5,
        "site_3":      site_3,
        "insert_seq":  insert,
        "fwd_pos":     fwd_pos,
        "rev_pos":     rev_pos,
    }


def _design_cloning_primers(
    template_seq: str,
    start: int,
    end: int,
    re_5prime: str,
    re_3prime: str,
    target_tm: float = 60.0,
    padding: str = "GCGC",
    tm_profile: str = _TM_PROFILE_DEFAULT,
) -> dict:
    """Design cloning primers using enzyme names from the combined
    catalog (built-in NEB ∪ user-added custom). Delegates to
    _design_cloning_primers_raw after looking up recognition sites."""
    catalog = _state._all_enzymes_hook()
    if re_5prime not in catalog:
        return {"error": f"Unknown enzyme: {re_5prime}"}
    if re_3prime not in catalog:
        return {"error": f"Unknown enzyme: {re_3prime}"}
    site_5, _, _ = catalog[re_5prime]
    site_3, _, _ = catalog[re_3prime]
    return _design_cloning_primers_raw(
        template_seq, start, end, site_5, site_3,
        name_5=re_5prime, name_3=re_3prime,
        target_tm=target_tm, padding=padding, tm_profile=tm_profile,
    )


@_timed("op.primer3.generic_design")
def _design_generic_primers(
    template_seq: str,
    start: int,
    end: int,
    target_tm: float = 60.0,
    tm_profile: str = _TM_PROFILE_DEFAULT,
) -> dict:
    """Design simple binding primers (no tails, no RE sites, no overhangs).

    Forward primer: optimal binding region at the start of the region.
    Reverse primer: optimal binding region at the end (reverse-complement).
    """
    total  = len(template_seq)
    insert = _slice_circular(template_seq.upper(), start, end)
    wraps  = end < start
    if len(insert) < 18:
        return {"error": "Region too short (< 18 bp)."}
    # No tail (binding-only primers) → the whole oligo IS the binding, so it
    # may grow up to the full `_PRIMER_MAX_OLIGO_LEN`.
    fwd_bind, fwd_tm = _pick_binding_region(
        insert, target_tm, max_len=_binding_max_len(0), tm_profile=tm_profile)
    rev_bind, rev_tm = _pick_binding_region(
        _rc(insert), target_tm, max_len=_binding_max_len(0),
        tm_profile=tm_profile)
    if wraps:
        fwd_pos = (start, (start + len(fwd_bind)) % total)
        rev_pos = ((end - len(rev_bind)) % total, end)
    else:
        fwd_pos = (start, start + len(fwd_bind))
        rev_pos = (end - len(rev_bind), end)
    return {
        "fwd_seq":  fwd_bind,
        "rev_seq":  rev_bind,
        "fwd_tm":   round(fwd_tm, 1),
        "rev_tm":   round(rev_tm, 1),
        "fwd_pos":  fwd_pos,
        "rev_pos":  rev_pos,
    }


# ── primer / primer-pair QC (check-primer-pair, design results, Primer Check) ─
# Tm under a profile, plus the secondary structures Primer3 screens for when it
# DESIGNS — each reported as ΔG at 37 °C (kcal/mol, the unit supplier sheets
# use) AND as the structure's MELTING TEMPERATURE, the number Primer3 itself
# tests (they match its design output, PRIMER_LEFT_0_HAIRPIN_TH and friends,
# digit for digit, floored at 0 °C as it reports them). The melting temperature
# is what decides: a hairpin at a harmless-looking −2.6 kcal/mol can melt at
# 62 °C, i.e. stay folded through the whole anneal.
#
# What gets FLAGGED was measured, not assumed (2026-09-30, 200 random designs a
# mode). Primer3's 47 °C limit on every structure flags 0 % of its own detection
# designs — it filters on exactly that — but 55 % of generic and 82 % of
# cloning designs: a GCGC + restriction-site tail folds onto the arm (median
# hairpin 48 °C) in primers labs use every day, and a warning on four cloning
# pairs in five is a warning nobody reads. Such a hairpin is REAL but rarely
# decisive — primers sit in vast excess, refold in microseconds, and from cycle
# 2 the template carries the tail's complement. What fails a PCR is a 3' end
# that cannot anneal or extend. So:
#   * a HAIRPIN is flagged when its stem takes in the primer's 3'-terminal
#     `_PRIMER_3P_END` bases (Primer3's own "3' end", PRIMER_MAX_END_STABILITY)
#     AND it is still folded at the anneal — Tm (concentration-independent for a
#     hairpin) at or above `_PRIMER_ANNEAL_BELOW_TM` under the lower primer's Tm
#     in PCR buffer, the conventional annealing temperature. 0 % / 12 % / 12 %
#     of detection / generic / cloning designs. Every hairpin's ΔG and Tm is
#     reported regardless, with `hairpin_3p` saying whether it holds the 3' end.
#   * a DIMER (self or pair; anywhere or on a 3' end) keeps Primer3's 47 °C
#     limit (PRIMER_MAX_SELF_ANY_TH / _SELF_END_TH / PAIR_MAX_COMPL_ANY_TH /
#     _COMPL_END_TH), below any anneal on purpose: a paired 3' end can be
#     extended while the reaction is set up and ramped, which is how a
#     primer-dimer band is made. ~1 % of designs; every real primer-dimer.
_PRIMER_DIMER_TM_LIMIT = 47.0
_PRIMER_ANNEAL_BELOW_TM = 5.0
_PRIMER_3P_END = 5
# primer3's thermodynamic alignment refuses a hairpin / self-dimer of a strand
# longer than this, and any pairing in which BOTH strands are longer.
_THAL_MAX_LEN = 60
# The 5 °C amplify-feature already warns at.
_PRIMER_PAIR_MAX_TM_DELTA = 5.0
_IUPAC_AMBIGUOUS = frozenset("RYSWKMBDHVN")


def _same_structure(a: "dict | None", b: "dict | None") -> bool:
    """Two `_primer_structure` results that are the SAME structure — the most
    stable dimer anywhere happening to be 3'-anchored — so it is flagged once."""
    return (a is not None and b is not None and a["dg"] == b["dg"]
            and a["tm"] == b["tm"])


def _primer_structure(kind: str, seq1: str, seq2: str = "") -> "dict | None":
    """One secondary structure in PCR buffer (`_MUT_P3` — never a Tm
    profile's conditions, see `_TM_PROFILES`): ``{"dg": kcal/mol at 37 °C
    (2 dp), "tm": melting temperature °C (1 dp, floored at 0)}``, or None when
    primer3 couldn't measure it (not installed, or it refused the input —
    callers check `_THAL_MAX_LEN` first, so they can say why). `kind` as in
    `_p3_thermo`. No structure at all comes back as ΔG 0 / Tm 0. A hairpin
    also carries ``"end_paired"``: does its stem take in any of the oligo's
    3'-terminal `_PRIMER_3P_END` bases."""
    hairpin = kind == "hairpin"
    try:
        res = _p3_thermo(kind, seq1, seq2, structure=hairpin, **_MUT_P3)
    except (ImportError, OSError, RuntimeError, ValueError, TypeError) as exc:
        _log.warning("_primer_structure: primer3 %s on %d/%d-mer not "
                     "measured: %s", kind, len(seq1), len(seq2), exc)
        return None
    out: dict = {"dg": round(res.dg / _CAL_PER_KCAL, 2),
                 "tm": round(max(0.0, float(res.tm)), 1)}
    if hairpin:
        # primer3 draws a hairpin as ONE line under the sequence: `/` and `\`
        # mark the two sides of the stem, `-` a free base.
        drawn = ""
        if res.structure_found:
            for line in (res.ascii_structure_lines or []):
                if line.startswith("SEQ\t"):
                    drawn = line[4:]
                    break
        out["end_paired"] = any(ch in "/\\" for ch in drawn[-_PRIMER_3P_END:])
    return out


def _anneal_ref(*arms: str) -> "float | None":
    """The conventional annealing temperature for these arms — the lower Tm IN
    PCR BUFFER (the default profile, whatever profile the Tm is reported
    under: the anneal happens in the tube) minus `_PRIMER_ANNEAL_BELOW_TM`."""
    tms = [t for t in (_primer_tm(a) for a in arms) if t is not None]
    return round(min(tms) - _PRIMER_ANNEAL_BELOW_TM, 1) if tms else None


def _primer_oligo_qc(primer: str, *, binding: "str | None" = None,
                     tm_profile: str = _TM_PROFILE_DEFAULT,
                     role: str = "primer",
                     anneal_ref: "float | None" = None) -> dict:
    """Tm + secondary-structure QC of ONE oligo — no template needed.

    `primer` is the whole oligo as ordered, 5'→3' (IUPAC; any 5' tail
    included); `binding` its 3' annealing arm — a 3' suffix of `primer` — or
    None when the whole oligo anneals. Tm and GC% are the ARM's (it alone
    anneals in the first cycles, so it sets the annealing temperature); the
    hairpin and self-dimers are the WHOLE oligo's (all of it is in the tube).
    `role` only words the warnings ("the forward primer …"). A hairpin is
    judged against `anneal_ref` — by default this oligo's own
    (`_anneal_ref`); a pair passes the lower of its two.

    Returns ``{primer, length, binding, binding_length, tail, tm, gc_pct,
    hairpin_dg, hairpin_tm, hairpin_3p, homodimer_dg, homodimer_tm,
    homodimer_3p_dg, homodimer_3p_tm, anneal_ref_c, warnings}`` (+
    ``tm_range`` for a degenerate arm) — `hairpin_3p` is whether the hairpin
    holds the 3' end, which is what gets it flagged (with a Tm at or above
    `anneal_ref_c`). A ``None`` ΔG/Tm means NOT measured, and a warning says
    why — never "no structure", which reads 0. Raises ValueError for an empty
    oligo, a `binding` that isn't a 3' suffix, or an unknown profile."""
    p = (primer or "").strip().upper()
    b = p if binding is None else (binding or "").strip().upper()
    if not p or not b:
        raise ValueError("empty primer or binding arm")
    if not p.endswith(b):
        raise ValueError("the binding arm must be a 3' suffix of the primer")
    cond = _tm_profile_kwargs(tm_profile)
    if anneal_ref is None:
        anneal_ref = _anneal_ref(b)
    warnings: "list[str]" = []
    out: dict = {
        "primer": p, "length": len(p),
        "binding": b, "binding_length": len(b),
        "tail": p[:len(p) - len(b)],
        "tm": _primer_tm(b, tm_profile),
        "gc_pct": round(100.0 * sum(1 for c in b if c in "GCS") / len(b), 1),
        "anneal_ref_c": anneal_ref,
    }
    bracket = _degenerate_tm_bracket(b, **cond)
    if bracket is not None:
        out["tm_range"] = [round(bracket[0], 1), round(bracket[1], 1)]
    names = (("hairpin", "hairpin"), ("homodimer", "homodimer"),
             ("homodimer_3p", "end_stability"))
    if len(p) > _THAL_MAX_LEN:
        for name, _k in names:
            out[f"{name}_dg"] = out[f"{name}_tm"] = None
        out["hairpin_3p"] = None
        warnings.append(
            f"the {role} primer is {len(p)} nt — primer3 folds oligos only "
            f"up to {_THAL_MAX_LEN} nt, so its hairpin and self-dimers were "
            f"not checked")
        out["warnings"] = warnings
        return out
    found = {name: _primer_structure(kind, p,
                                     p if kind == "end_stability" else "")
             for name, kind in names}
    for name, _k in names:
        st = found[name]
        out[f"{name}_dg"] = None if st is None else st["dg"]
        out[f"{name}_tm"] = None if st is None else st["tm"]
    hp = found["hairpin"]
    out["hairpin_3p"] = None if hp is None else hp["end_paired"]
    if (hp is not None and hp["end_paired"] and anneal_ref is not None
            and hp["tm"] >= anneal_ref):
        warnings.append(
            f"the {role} primer's 3' end is held in a hairpin that melts at "
            f"{hp['tm']:.1f} °C (ΔG {hp['dg']:.2f} kcal/mol) — still folded "
            f"at a {anneal_ref:.1f} °C anneal ({_PRIMER_ANNEAL_BELOW_TM:.0f} "
            f"°C under the lower primer Tm), so it can't prime")
    homo, homo3 = found["homodimer"], found["homodimer_3p"]
    for st, words in ((homo, "a self-dimer"),
                      (homo3, "a self-dimer on its 3' end (extendable)")):
        if st is None or st["tm"] <= _PRIMER_DIMER_TM_LIMIT:
            continue
        if st is homo and _same_structure(homo, homo3):
            continue            # one structure: the 3' warning says more
        warnings.append(
            f"the {role} primer forms {words} stable up to {st['tm']:.1f} °C "
            f"(ΔG {st['dg']:.2f} kcal/mol) — Primer3 rejects one stable "
            f"above {_PRIMER_DIMER_TM_LIMIT:.0f} °C")
    if any(v is None for v in found.values()):
        warnings.append(
            f"the {role} primer's hairpin / self-dimer could not be computed "
            f"(primer3 unavailable) — it is NOT secondary-structure checked")
    elif set(p) & _IUPAC_AMBIGUOUS:
        warnings.append(
            f"the {role} primer carries ambiguity codes — its structures "
            f"count them as non-pairing, so the worst member of the mix can "
            f"fold more than shown")
    out["warnings"] = warnings
    return out


def _primer_pair_qc(fwd: str, rev: str, *,
                    fwd_binding: "str | None" = None,
                    rev_binding: "str | None" = None,
                    tm_profile: str = _TM_PROFILE_DEFAULT) -> dict:
    """Thermodynamic QC of one PCR primer pair — no template needed.

    Each oligo gets `_primer_oligo_qc` (its Tm from the annealing arm, its
    structures from the whole oligo), with ONE annealing reference for both —
    the reaction has one anneal, set by the weaker arm. The pair adds the Tm
    difference of the two arms and the primer-dimers between the two WHOLE
    oligos: anywhere, and with either 3' end anchored (the worse direction).

    Returns ``{forward, reverse, pair: {tm_delta, heterodimer_dg,
    heterodimer_tm, heterodimer_3p_dg, heterodimer_3p_tm}, tm_profile,
    dg_units, anneal_ref_c, dimer_tm_limit, warnings}`` with every oligo's
    warnings folded into the one list (empty = nothing to flag). Raises
    ValueError as `_primer_oligo_qc` does."""
    key = _tm_profile_resolve(tm_profile)
    if key is None:
        raise ValueError(f"unknown tm_profile {tm_profile!r}")
    fb = (fwd_binding if fwd_binding is not None else fwd or "").strip()
    rb = (rev_binding if rev_binding is not None else rev or "").strip()
    ref = _anneal_ref(fb.upper(), rb.upper()) if fb and rb else None
    f = _primer_oligo_qc(fwd, binding=fwd_binding, tm_profile=key,
                         role="forward", anneal_ref=ref)
    r = _primer_oligo_qc(rev, binding=rev_binding, tm_profile=key,
                         role="reverse", anneal_ref=ref)
    warnings = f.pop("warnings") + r.pop("warnings")
    pair: dict = {"tm_delta": None}
    if f["tm"] is not None and r["tm"] is not None:
        pair["tm_delta"] = round(abs(f["tm"] - r["tm"]), 1)
        if pair["tm_delta"] > _PRIMER_PAIR_MAX_TM_DELTA:
            warnings.append(f"primer Tms differ by {pair['tm_delta']:.1f} °C")
    fp, rp = f["primer"], r["primer"]
    if len(fp) > _THAL_MAX_LEN and len(rp) > _THAL_MAX_LEN:
        pair.update(heterodimer_dg=None, heterodimer_tm=None,
                    heterodimer_3p_dg=None, heterodimer_3p_tm=None)
        warnings.append(
            f"both primers are longer than {_THAL_MAX_LEN} nt — primer3 "
            f"cannot align them against each other, so primer-dimers were "
            f"not checked")
    else:
        het = _primer_structure("heterodimer", fp, rp)
        ends = [st for st in (_primer_structure("end_stability", fp, rp),
                              _primer_structure("end_stability", rp, fp))
                if st is not None]
        end = max(ends, key=lambda st: st["tm"]) if len(ends) == 2 else None
        pair["heterodimer_dg"] = None if het is None else het["dg"]
        pair["heterodimer_tm"] = None if het is None else het["tm"]
        pair["heterodimer_3p_dg"] = None if end is None else end["dg"]
        pair["heterodimer_3p_tm"] = None if end is None else end["tm"]
        if het is None or end is None:
            warnings.append(
                "the primer-dimer between the two primers could not be "
                "computed (primer3 unavailable) — the pair is NOT "
                "secondary-structure checked")
        for st, words in ((het, "pair with each other"),
                          (end, "pair at a 3' end (extendable)")):
            if st is None or st["tm"] <= _PRIMER_DIMER_TM_LIMIT:
                continue
            if st is het and _same_structure(het, end):
                continue        # one structure: the 3' warning says more
            warnings.append(
                f"the two primers {words} up to {st['tm']:.1f} °C "
                f"(ΔG {st['dg']:.2f} kcal/mol) — Primer3 rejects a pair "
                f"stable above {_PRIMER_DIMER_TM_LIMIT:.0f} °C")
    return {"forward": f, "reverse": r, "pair": pair, "tm_profile": key,
            "dg_units": "kcal/mol", "anneal_ref_c": ref,
            "dimer_tm_limit": _PRIMER_DIMER_TM_LIMIT, "warnings": warnings}


# ── mutagenesis-INNER design (Phase D, moved from hub) ─────────────────────
# The QuikChange inner primers that INTRODUCE the codon change (mutant codon
# chosen via the codon-usage map _codon_build_aa_map, codon L2), the folded
# modified-outer, and CDS extraction (origin-wrap + reverse-strand aware).
# Verified by the real-CDS mutagenesis-inner golden (byte-identical output).
_MUT_MIN_SOE_FRAG  = 60                             # nt; below this → edge case


# A folded "modified outer" primer (2-primer direct PCR, no SOE) MUST span
# the mutant codon, otherwise it amplifies wild-type — a silent wrong product.
# Cap its length and require a matched extension anchor past the mutation; when
# no spanning primer fits, the caller falls back to the regular SOE design.
_MUT_OUTER_MAX_LEN = 45    # nt; mutation farther than this from the end → SOE


_MUT_OUTER_ANCHOR  = 8     # matched template bases required past the mutation


def _mut_design_modified_outer(dna_mut: str, near_start: bool,
                               nt_start: int) -> "dict | None":
    """Edge-case: mutation < _MUT_MIN_SOE_FRAG nt from a CDS end → fold the
    mutant codon into a single outer primer so the PCR is a 2-primer direct
    reaction (no SOE).

    The folded primer MUST span the mutant codon at `nt_start` (with a matched
    anchor for clean extension), else it would carry wild-type sequence and
    silently amplify a WT product. Returns None when no spanning primer fits
    within `_MUT_OUTER_MAX_LEN`, so the caller falls back to the regular
    (always-correct) SOE inner+outer design."""
    mut_end = nt_start + 3                      # exclusive end of mutant codon
    if nt_start < 0 or mut_end > len(dna_mut):
        return None
    if near_start:
        # FWD outer anneals from base 3 (after the BsaI tail) toward 3'; it
        # must start at/before the codon and reach mut_end + anchor.
        anneal_start = 3
        if nt_start < anneal_start:
            return None
        end_hi = min(len(dna_mut), anneal_start + _MUT_OUTER_MAX_LEN)
        best = None
        for end in range(mut_end + _MUT_OUTER_ANCHOR, end_hi + 1):
            anneal = dna_mut[anneal_start:end]
            if len(anneal) < 18:
                continue
            s = _mut_score_outer(anneal)
            if best is None or s < best["score"]:
                best = {
                    "anneal":    anneal,
                    "full":      _MUT_BSAI_FWD_TAIL + anneal,
                    "tm_anneal": _mut_tm(anneal),
                    "gc":        _mut_gc_pct(anneal),
                    "score":     s,
                }
        if best is None:
            return None
        best["label"]    = "modified_FWD_outer"
        best["partner"]  = "REV_outer (unchanged)"
        best["replaces"] = "FWD_outer"
        return best
    # near_end: REV outer anneals to the bottom strand at the 3' end; the
    # top-strand window [start, len) must reach back past the mutant codon.
    seq_len  = len(dna_mut)
    lo_start = max(0, seq_len - _MUT_OUTER_MAX_LEN)
    hi_start = nt_start - _MUT_OUTER_ANCHOR
    best = None
    for start in range(lo_start, hi_start + 1):
        tail = dna_mut[start:]
        if len(tail) < 18:
            continue
        anneal = _mut_revcomp(tail)
        s = _mut_score_outer(anneal)
        if best is None or s < best["score"]:
            best = {
                "anneal":    anneal,
                "full":      _MUT_BSAI_REV_TAIL + anneal,
                "tm_anneal": _mut_tm(anneal),
                "gc":        _mut_gc_pct(anneal),
                "score":     s,
            }
    if best is None:
        return None
    best["label"]    = "modified_REV_outer"
    best["partner"]  = "FWD_outer (unchanged)"
    best["replaces"] = "REV_outer"
    return best


def _mut_design_inner(dna: str, mut_pos_1: int, mut_aa: str, wt_aa: str,
                      codon_table: "dict | None" = None) -> dict:
    """Inner mutagenic pair (FWD carries mutant codon; REV = revcomp(FWD)).

    `codon_table` is an optional {codon: (aa, count)} map used to pick the
    mutant codon. Defaults to E. coli K12 (_MUT_AA_TO_CODONS)."""
    idx      = mut_pos_1 - 1
    nt_start = idx * 3

    wt_codon  = dna[nt_start:nt_start + 3]
    if len(wt_codon) < 3:
        raise ValueError(
            f"Position {mut_pos_1} is past the end of the CDS."
        )
    wt_actual = _MUT_CODON_TO_AA.get(wt_codon, "?")
    if wt_actual != wt_aa:
        raise ValueError(
            f"Position {mut_pos_1}: mutation says WT='{wt_aa}' but DNA codon "
            f"'{wt_codon}' encodes '{wt_actual}'."
        )

    if codon_table:
        aa_map, _ = _codon_build_aa_map(codon_table)
    else:
        aa_map = _MUT_AA_TO_CODONS

    if mut_aa == "*":
        # Prefer the selected table's most-frequent stop codon (TAA fallback)
        # so the suggestion respects the organism instead of always TAA.
        # _codon_build_aa_map excludes stops, so read them from the raw table.
        stops = [
            (str(c).upper(), v)
            for c, v in (codon_table or {}).items()
            if isinstance(v, (list, tuple)) and len(v) >= 2 and v[0] == "*"
        ]

        def _stop_count(item) -> float:
            try:
                return float(item[1][1])
            except (TypeError, ValueError, IndexError):
                return 0.0

        mut_codon = max(stops, key=_stop_count)[0] if stops else "TAA"
    else:
        mut_codon = next(
            (c for c, _f in aa_map.get(mut_aa, []) if c != wt_codon),
            None,
        )
        if mut_codon is None:
            raise ValueError(f"No alternative codon available for '{mut_aa}' "
                             "in the selected codon table")

    mut_dna = dna[:nt_start] + mut_codon + dna[nt_start + 3:]

    TM_TARGET      = 60.0
    TM_MIN, TM_MAX = 55.0, 75.0
    GC_MIN, GC_MAX = 35.0, 68.0
    seq_len = len(mut_dna)

    candidates: list = []
    # Why each window was turned down, so the error names the real dead end
    # (the old message blamed the sequence ends for Tm / GC failures too).
    rejected = {"length": 0, "tm": 0, "gc": 0}
    for left_ext in range(5, 28):
        for right_ext in range(5, 28):
            lo  = max(0, nt_start - left_ext)
            hi  = min(seq_len, nt_start + 3 + right_ext)
            fwd = mut_dna[lo:hi]
            if len(fwd) < 15 or len(fwd) > 58:
                # A LENGTH dead end only when a sequence end cut the window
                # short. The smallest extensions are under 15 nt by
                # construction, so counting them put "too close to a sequence
                # end" into every failure, mid-sequence ones included.
                if lo > nt_start - left_ext or hi < nt_start + 3 + right_ext:
                    rejected["length"] += 1
                continue
            t  = _mut_tm(fwd)
            gc = _mut_gc_pct(fwd)
            if not (TM_MIN <= t <= TM_MAX):
                rejected["tm"] += 1
                continue
            if not (GC_MIN <= gc <= GC_MAX):
                rejected["gc"] += 1
                continue
            hp = _mut_hairpin_dg(fwd)
            hd = _mut_homodimer_dg(fwd)
            # kcal/mol; None (not measured) scores as no penalty but is
            # reported verbatim below.
            score = (
                abs(t - TM_TARGET) * 2.0
                + (0 if _mut_ends_gc(fwd) else 4.0)
                + max(0.0, -(hp or 0.0) - 1.0) / 0.4
                + max(0.0, -(hd or 0.0) - 2.0) / 0.4
                + abs(gc - 50) * 0.1
                - (len(fwd) * 0.15 if abs(t - TM_TARGET) <= 1.0 else 0)
            )
            candidates.append({
                "fwd": fwd, "rev": _mut_revcomp(fwd),
                "tm": t, "gc": gc, "length": len(fwd),
                # kcal/mol (None = not measured); see `_mut_hairpin_dg`.
                "hairpin_dg": None if hp is None else round(hp, 2),
                "homodimer_dg": None if hd is None else round(hd, 2),
                "dg_units": "kcal/mol",
                "score": score, "lo": lo,
            })

    if not candidates:
        why_of = {
            "length": "the mutation is too close to a sequence end to fit a "
                      "15 nt primer around it",
            "tm":     f"no window around it reaches a Tm of "
                      f"{TM_MIN:.0f}-{TM_MAX:.0f} °C (too AT- or GC-rich)",
            "gc":     f"no window around it has {GC_MIN:.0f}-{GC_MAX:.0f} % "
                      f"GC",
        }
        # EVERY reason that rejected a window, commonest first — not just
        # `max(rejected, ...)`. That returned the FIRST key on a tie, and the
        # dict is ordered `length, tm, gc`, so a pure Tm or GC dead end was
        # reported as "too close to a sequence end" and sent the user to move a
        # mutation that was nowhere near one (audit 2026-09-22). The counts are
        # included because they say which constraint to relax.
        ranked_reasons = [(k, n) for k, n in
                          sorted(rejected.items(), key=lambda kv: -kv[1])
                          if n > 0 and k in why_of]
        if ranked_reasons:
            why = "; ".join(f"{why_of[k]} ({n} window"
                            f"{'' if n == 1 else 's'})"
                            for k, n in ranked_reasons)
        else:
            why = "no candidate window was generated at all"
        raise RuntimeError(
            f"No valid inner primers found for {wt_aa}{mut_pos_1}{mut_aa}: "
            f"{why}."
        )

    seen: dict = {}
    for c in sorted(candidates, key=lambda x: x["score"]):
        if c["fwd"] not in seen:
            seen[c["fwd"]] = c
    ranked = sorted(seen.values(), key=lambda x: x["score"])[:5]
    for i, c in enumerate(ranked):
        c["rank"] = i + 1

    best_lo = ranked[0]["lo"]
    best_hi = best_lo + ranked[0]["length"]
    fwd_anneal_start = 3
    frag_a = best_hi - fwd_anneal_start
    frag_b = seq_len - best_lo

    near_start = frag_a < _MUT_MIN_SOE_FRAG
    near_end   = frag_b < _MUT_MIN_SOE_FRAG

    edge_case = None
    if near_start or near_end:
        modified_outer = _mut_design_modified_outer(
            mut_dna, near_start=near_start, nt_start=nt_start)
        # Only offer the 2-primer direct shortcut when the folded primer
        # actually spans the mutation; otherwise edge_case stays None and the
        # regular SOE inner+outer design (always correct) is shown/saved.
        # Pre-fix the modified outer used a fixed start-anchored window and
        # could omit the mutation for codons ~11-19 → silent WT product.
        if modified_outer is not None:
            edge_case = {
                "near_start":     near_start,
                "near_end":       near_end,
                "frag_a":         frag_a,
                "frag_b":         frag_b,
                "modified_outer": modified_outer,
            }

    return {
        "mutation":    f"{wt_aa}{mut_pos_1}{mut_aa}",
        "nt_position": nt_start + 1,
        "wt_codon":    wt_codon,
        "mut_codon":   mut_codon,
        "nt_changes":  sum(a != b for a, b in zip(wt_codon, mut_codon)),
        "candidates":  ranked,
        "edge_case":   edge_case,
    }


def _mut_extract_cds(full_seq: str, start: int, end: int, strand: int) -> str:
    """Return the CDS DNA in its biological 5'→3' orientation, handling
    origin-wrap (end < start) and reverse-strand features."""
    if end < start:
        sub = full_seq[start:] + full_seq[:end]
    else:
        sub = full_seq[start:end]
    sub = sub.upper()
    if strand == -1:
        sub = _mut_revcomp(sub)
    return sub


# ═══ QuikChange-scrub engine (moved from the hub) ══════════════════════════
# Clone-free restriction-site removal: the planner (_scrub_design + helpers)
# and the improved-QuikChange primer designer (_scrub_qc_*). Reaches the two
# hub-pinned deps via _state hooks (_all_enzymes_hook / _translate_cds_hook).

# ── Scrub: clone-free restriction-site removal (improved QuikChange) ────────
#
# "Scrub" cures a plasmid of chosen restriction recognition sites by
# introducing the MINIMAL set of point substitutions that destroy each
# site WITHOUT (a) changing any overlapping CDS's protein (synonymous /
# silent), (b) creating a NEW copy of any forbidden site anywhere, or
# (c) — softly — touching annotated non-coding elements when avoidable.
# The result is a single contiguous circular sequence the user makes in
# the lab by improved-QuikChange whole-plasmid PCR + self-recircularization
# (DpnI, transform — no ligase, no fragment assembly). See docs/invariants
# [INV-97] and the MutagenizeModal "Scrub" tab.
#
# Why it leans on the hardened primitives instead of rolling its own:
#   * Site enumeration AND the final "did I introduce a new site?" guard go
#     through `_scan_restriction_sites(circular=True)`, so origin-spanning
#     sites (sacred inv #6) and reverse-strand coords (sacred inv #2) are
#     handled by code that already has tests.
#   * Synonymy is verified by RE-TRANSLATING every overlapping CDS with
#     `_translate_cds` (the canonical translator) rather than hand-mapping
#     codons — wrap, reverse strand, /codon_start and /transl_table all come
#     out right for free, including a site straddling TWO CDSes on opposite
#     strands (both frames must stay synonymous or the change is rejected).
#   * Every cure is a 1-base→1-base SUBSTITUTION, so the cured sequence is
#     the SAME LENGTH and feature coordinates never shift — the caller can
#     swap `record.seq` and keep features verbatim (no `_rebuild_record_
#     with_edit`, no coordinate migration).

# Marquee use case is clone-free MoClo / Golden Braid domestication, so the
# default target set is the Type IIS cutters assembly can't tolerate inside a
# part. Esp3I covers its isoschizomer BsmBI (same CGTCTC site).
_SCRUB_DEFAULT_ENZYMES: "tuple[str, ...]" = ("BsaI", "Esp3I", "BbsI")
# Edits within this many bp share a single QuikChange primer pair (one PCR
# round); edits farther apart need separate sequential rounds.
_SCRUB_PRIMER_FOOTPRINT = 30
# Cap on substitutions tried to kill one site. A 6-cutter almost always dies
# to a single change; 3 covers a constrained CDS where the obvious wobble
# base isn't free. Beyond this we report the site as un-scrubbable rather
# than mangling a long stretch.
_SCRUB_MAX_CHANGES = 3

# Homopolymer hygiene for a cure (2026-08-26 field report). A single-base
# swap that happens to butt onto an existing run can LENGTHEN it — curing
# `GGAAAAA GAAGAC` by G→A spells `GGAAAAAAAAGAC`, an 8×A run, where the
# equally valid A→G at the third base spells `GGAAAAA GAGGAC` and leaves the
# longest run at 5. Both are single transitions with the same GC delta, so
# the old score separated them only by the lexicographic determinism
# tie-break — and `A` sorts before `G`, so the run-extending cure won every
# time. Long runs are a synthesis and sequencing hazard (slippage), and the
# report's case doubled it: the site sat in a TANDEM-duplicated 35S enhancer,
# so the same bad choice landed twice.
#
# `_SCRUB_HOMOPOLYMER_FREE` is the run length a cure may reach for free.
# Scoring is RELATIVE to what was already there, so a cure is never blamed
# for a run it didn't make — only for making one longer.
_SCRUB_HOMOPOLYMER_FREE = 5
# How far either side of the site to look. A run is only interesting where it
# meets the edit, and this bounds the per-candidate cost.
_SCRUB_HOMOPOLYMER_FLANK = 15


def _scrub_longest_run(window: str) -> int:
    """Longest single-base homopolymer run in `window` (0 for empty).

    Local to the scrub engine (L2) rather than reusing seqanalysis'
    `_lint_homopolymer_runs` (L3) — L2 cannot import upward, and this needs
    one integer per candidate, not a list of annotated findings."""
    best = run = 0
    prev = ""
    for ch in window:
        run = run + 1 if ch == prev else 1
        prev = ch
        if run > best:
            best = run
    return best


def _circ_window(seq: str, start: int, length: int, n: int) -> str:
    """Return `length` bases of `seq` starting at forward coord `start`,
    wrapping the origin when `start + length > n` (`n == len(seq)`)."""
    end = start + length
    if end <= n:
        return seq[start:end]
    return seq[start:] + seq[:end - n]


def _scrub_pos_in_feat(g: int, s: int, e: int) -> bool:
    """Wrap-aware membership: is genome coord `g` inside span `[s, e)`
    (where `e < s` signals an origin wrap)? Module-level mirror of
    `PlasmidMap._bp_in` for the Scrub path."""
    if e >= s:
        return s <= g < e
    return g >= s or g < e


def _scrub_is_transition(a: str, b: str) -> bool:
    """A↔G / C↔T are transitions; the other four swaps are transversions.
    Used only as a tie-break — transitions are the milder substitution."""
    return {a, b} in ({"A", "G"}, {"C", "T"})


def _scrub_resolve_sites(enzymes) -> "dict[str, str]":
    """Resolve enzyme names → forward recognition sites from the merged
    built-in + custom catalog. Unknown names and non-IUPAC sites are
    skipped (one bad custom enzyme never aborts the whole scrub)."""
    enz = _state._all_enzymes_hook()
    out: "dict[str, str]" = {}
    for nm in enzymes or ():
        info = enz.get(str(nm))
        if not info:
            continue
        site = str(info[0] or "").upper()
        if not site:
            continue
        try:
            _iupac_pattern(site)
        except ValueError:
            continue
        out[str(nm)] = site
    return out


def _scrub_expand_forbidden(forward: "dict[str, str]") -> "tuple[str, ...]":
    """Forward sites + their reverse complements (deduped) as the flat
    tuple `_forbidden_hit_set` consumes — so the 'no new site' guard vetoes
    a swap that would spawn a forbidden site on EITHER strand."""
    out: "list[str]" = []
    for site in forward.values():
        if site not in out:
            out.append(site)
        rc = _rc(site)
        if rc not in out:
            out.append(rc)
    return tuple(out)


def _scrub_scan_targets(seq: str, allowed: "frozenset[str]",
                        circular: bool) -> "list[dict]":
    """One entry per recognition-site INSTANCE of the allowed enzymes:
    `{enzyme, strand, rec_start, rec_end, positions}` where `positions` is
    the wrap-aware list of genome coords the recognition covers. Counts
    only LABELED resite pieces (sacred inv #6 — a wrap hit's unlabeled
    head piece is a continuation, not a second site)."""
    n = len(seq)
    hits = _scan_restriction_sites(
        seq, min_recognition_len=1, unique_only=False,
        circular=circular, allowed_enzymes=allowed)
    out: "list[dict]" = []
    for h in hits:
        if h.get("type") != "resite" or not h.get("label"):
            continue
        rs = h.get("rec_start", h["start"])
        re_ = h.get("rec_end", h["end"])
        if re_ < rs:                      # origin wrap
            positions = list(range(rs, n)) + list(range(0, re_))
        else:
            positions = list(range(rs, re_))
        if not positions:
            continue
        out.append({
            "enzyme": h["label"], "strand": h.get("strand", 1),
            "rec_start": rs, "rec_end": re_, "positions": positions,
        })
    out.sort(key=lambda t: (t["rec_start"], t["enzyme"], t["strand"]))
    return out


def _scrub_overlapping_feats(target: dict, feats: list) -> tuple:
    """Split `feats` into (cds_feats, other_feats) that overlap the target
    recognition window. Only 'CDS' constrains synonymy; everything else
    (gene/promoter/RBS/rep_origin/…) is 'annotated non-coding' and only
    earns a soft penalty for being touched."""
    posset = set(target["positions"])
    cds: list = []
    other: list = []
    for f in feats:
        s = f.get("start")
        e = f.get("end")
        if s is None or e is None:
            continue
        if not any(_scrub_pos_in_feat(g, s, e) for g in posset):
            continue
        if f.get("type") == "CDS":
            cds.append(f)
        else:
            other.append(f)
    return cds, other


def _scrub_region_label(cds_feats: list, other_feats: list) -> str:
    """Human label for where a site sits, for the preview table."""
    if cds_feats:
        names = ", ".join(str(f.get("label") or "CDS") for f in cds_feats)
        return f"CDS: {names}"
    if other_feats:
        names = ", ".join(str(f.get("label") or f.get("type") or "?")
                          for f in other_feats)
        return f"non-coding (in {names})"
    return "non-coding"


def _scrub_cds_protein(seq: str, f: dict) -> str:
    """Translate one CDS feature dict on `seq` via the canonical translator,
    honouring wrap / strand / _exons / codon_start / transl_table."""
    return _state._translate_cds_hook(
        seq,
        f.get("_orig_start", f["start"]),
        f.get("_orig_end", f["end"]),
        f.get("strand", 1),
        f.get("_exons"),
        int(f.get("codon_start", 1) or 1),
        int(f.get("transl_table", 1) or 1),
    )


def _scrub_introduces_site(orig: str, test: str, target: dict,
                           all_forbidden: "tuple[str, ...]", n: int) -> bool:
    """True if `test` has a forbidden site near the edit that `orig` didn't.

    A substitution can only create a new site within (L-1) bp of a changed
    base (L = longest forbidden site), so a padded circular window around
    the recognition covers every site this change could spawn — comparing
    before/after hit sets on that SAME window is wrap-safe (a site present
    at the window's own boundary in both is not 'new'). The caller's final
    full re-scan is the authoritative whole-plasmid guard; this is the fast
    per-candidate pre-filter so we PICK a clean candidate."""
    L = max((len(s) for s in all_forbidden), default=8)
    pad = max(L - 1, 1)
    win_len = len(target["positions"]) + 2 * pad
    if win_len >= n:                       # tiny plasmid — scan it all
        a = orig + orig[:L - 1]
        b = test + test[:L - 1]
    else:
        start = (target["rec_start"] - pad) % n
        a = _circ_window(orig, start, win_len, n)
        b = _circ_window(test, start, win_len, n)
    return bool(_forbidden_hit_set(b, all_forbidden)
                - _forbidden_hit_set(a, all_forbidden))


def _scrub_cds_reading_positions(f: dict, n: int) -> tuple:
    """Genome positions of a CDS feature in 5'→3' CODING order (after the
    ``/codon_start`` offset), plus its strand. Codon ``i`` is positions
    ``[3i:3i+3]``; on the reverse strand the coding base at each position is
    the complement of the top strand. Mirrors the slice logic in
    `_cds_aa_list` / `_translate_cds` so frequency look-ups land on exactly
    the codons the synonymy check translates."""
    s = f.get("_orig_start", f["start"])
    e = f.get("_orig_end", f["end"])
    strand = f.get("strand", 1)
    exons = f.get("_exons")
    cs = max(0, min(2, int(f.get("codon_start", 1) or 1) - 1))
    if exons:
        gpos: list = []
        for xs, xe in exons:
            gpos.extend(range(xs, xe))
    elif e < s:                              # origin wrap
        gpos = list(range(s, n)) + list(range(0, e))
    else:
        gpos = list(range(s, e))
    if strand == -1:
        gpos = gpos[::-1]
    return gpos[cs:], strand


def _scrub_one_site(seq: str, target: dict, feats: list,
                    forward: "dict[str, str]",
                    all_forbidden: "tuple[str, ...]", n: int, *,
                    codon_frac: "dict | None" = None,
                    max_changes: int = _SCRUB_MAX_CHANGES) -> tuple:
    """Find the minimal substitution set that destroys ONE site instance.

    Returns `(changes, region, None)` where `changes = {genome_pos:
    new_base}`, or `(None, region, reason)` when the site can't be cured
    silently. A candidate is accepted only if it (1) destroys this
    instance, (3) leaves every overlapping CDS's protein unchanged, and
    (2) introduces no new forbidden site near the edit. Among the accepted,
    the score prefers: fewest changes → not touching annotated bases → not
    LENGTHENING a homopolymer run past `_SCRUB_HOMOPOLYMER_FREE` → (when
    `codon_frac` is given) the host-frequent synonymous codon → transitions
    over transversions → GC-neutral → lexicographic (deterministic, so the
    same plasmid always cures the same way).

    The homopolymer term ranks above codon usage on purpose: a run long
    enough to stall synthesis is a hard failure, a slightly rarer synonymous
    codon is a soft cost. It is measured RELATIVE to the sequence that was
    already there, so a cure is never penalised for a run it inherited.

    `codon_frac` (``{codon: usage_fraction}``) is a TIE-BREAK ONLY: synonymy
    is already guaranteed by check (3), so even a wrong frame map could only
    pick a less-preferred — never a non-silent — cure."""
    from itertools import combinations, product
    positions = target["positions"]
    site_len = len(positions)
    cds_feats, other_feats = _scrub_overlapping_feats(target, feats)
    region = _scrub_region_label(cds_feats, other_feats)
    orig_aa = [_scrub_cds_protein(seq, f) for f in cds_feats]
    annotated: set = set()
    for f in other_feats:
        s = f.get("start")
        e = f.get("end")
        for g in positions:
            if _scrub_pos_in_feat(g, s, e):
                annotated.add(g)
    # Reading-frame maps for overlapping CDSes, so a coding cure can prefer
    # the host-frequent synonymous codon. Built once per site (not per
    # candidate). Skipped entirely when no codon table was supplied.
    cds_frames: list = []
    if codon_frac:
        for f in cds_feats:
            gpos, st = _scrub_cds_reading_positions(f, n)
            cds_frames.append((gpos, st, {g: j for j, g in enumerate(gpos)}))
    fwd_site = forward[target["enzyme"]]
    pat_f = _iupac_pattern(fwd_site)
    pat_r = _iupac_pattern(_rc(fwd_site))
    # Homopolymer baseline: the same window, measured on the UNEDITED
    # sequence. Scoring the delta (not the absolute run) is what keeps a cure
    # inside an existing poly-A tract from being scored as if it made it.
    hp_start = (target["rec_start"] - _SCRUB_HOMOPOLYMER_FLANK) % n if n else 0
    hp_len = min(site_len + 2 * _SCRUB_HOMOPOLYMER_FLANK, n) if n else 0
    hp_before = max(_scrub_longest_run(_circ_window(seq, hp_start, hp_len, n))
                    if hp_len else 0, _SCRUB_HOMOPOLYMER_FREE)
    best: "tuple | None" = None
    for k in range(1, max_changes + 1):
        for combo in combinations(positions, k):
            # Each chosen position is forced to a DIFFERENT base, so every
            # candidate makes exactly `k` real changes (no silent no-ops).
            alts = [[b for b in "ACGT" if b != seq[g]] for g in combo]
            for repl in product(*alts):
                tl = list(seq)
                for j, g in enumerate(combo):
                    tl[g] = repl[j]
                test = "".join(tl)
                # (1) this instance must be gone (check the exact window,
                #     both strands so a reverse-strand hit is covered).
                win = _circ_window(test, target["rec_start"], site_len, n)
                if pat_f.fullmatch(win) or pat_r.fullmatch(win):
                    continue
                # (3) synonymous in every overlapping CDS frame.
                if any(_scrub_cds_protein(test, f) != orig_aa[i]
                       for i, f in enumerate(cds_feats)):
                    continue
                # (2) no NEW forbidden site near the edit.
                if _scrub_introduces_site(seq, test, target,
                                          all_forbidden, n):
                    continue
                # Codon-usage preference: total fraction of the resulting
                # codons at each changed coding position (higher = better).
                freq = 0.0
                for gpos, st, idx in cds_frames:
                    for g in combo:
                        jj = idx.get(g)
                        if jj is None:
                            continue
                        trip = gpos[(jj // 3) * 3:(jj // 3) * 3 + 3]
                        if len(trip) < 3:
                            continue
                        codon = ("".join(_rc(test[p]) for p in trip)
                                 if st == -1 else
                                 "".join(test[p] for p in trip))
                        freq += (codon_frac or {}).get(codon, 0.0)
                # How much longer this cure makes the longest local run than
                # it already was. 0 for almost every cure — it only bites
                # where a swap butts onto a run, which is exactly the case
                # the plain lexicographic tie-break used to get wrong.
                hp_cost = 0
                if hp_len:
                    hp_cost = max(0, _scrub_longest_run(
                        _circ_window(test, hp_start, hp_len, n)) - hp_before)
                score = (
                    k,
                    sum(1 for g in combo if g in annotated),
                    hp_cost,
                    -round(freq, 6),
                    sum(1 for j in range(k)
                        if not _scrub_is_transition(seq[combo[j]], repl[j])),
                    sum(1 for j in range(k)
                        if (seq[combo[j]] in "GC") != (repl[j] in "GC")),
                    repl,
                )
                if best is None or score < best[0]:
                    best = (score, {combo[j]: repl[j] for j in range(k)})
        if best is not None:        # found the minimal-k solution; stop
            break
    if best is None:
        return None, region, ("no silent substitution removes this site "
                              "without creating another")
    return best[1], region, None


def _scrub_cluster_edits(positions: list, n: int,
                         footprint: int = _SCRUB_PRIMER_FOOTPRINT) -> list:
    """Group edit positions into QuikChange rounds: consecutive positions
    within `footprint` bp share one primer pair. Merges the first and last
    clusters when they sit within `footprint` across the origin."""
    if not positions:
        return []
    ps = sorted(set(positions))
    clusters: "list[list[int]]" = [[ps[0]]]
    for p in ps[1:]:
        if p - clusters[-1][-1] <= footprint:
            clusters[-1].append(p)
        else:
            clusters.append([p])
    if len(clusters) > 1:
        # Distance from the last cluster's tail, across the origin, to the
        # first cluster's head.
        if (clusters[0][0] + n) - clusters[-1][-1] <= footprint:
            clusters[0] = clusters[-1] + clusters[0]
            clusters.pop()
    return clusters


def _scrub_design(seq: str, feats: "list | None" = None,
                  enzymes=None, *, circular: bool = True,
                  codon_raw: "dict | None" = None) -> dict:
    """Plan a clone-free scrub of `seq`. Returns a dict with the cured
    sequence, the per-base edits, which sites were removed vs skipped, the
    QuikChange round clustering, and any warnings. Pure biology — designs
    no primers (see `_scrub_qc_primers`) and mutates nothing on disk.

    `codon_raw` (an organism's ``{codon: (aa, count)}`` usage table) is an
    optional preference: when a coding site has more than one silent cure,
    the host-frequent synonymous codon wins. Omit it and ties fall back to
    the deterministic transition/GC/lexicographic order."""
    seq = (seq or "").upper()
    n = len(seq)
    forward = _scrub_resolve_sites(
        enzymes if enzymes is not None else _SCRUB_DEFAULT_ENZYMES)
    codon_frac: "dict | None" = None
    if codon_raw:
        try:
            _, codon_frac = _codon_build_aa_map(codon_raw)
        except Exception:
            codon_frac = None
    result: dict = {
        "ok": True, "orig_seq": seq, "cured_seq": seq,
        "enzymes": sorted(forward.keys()),
        "edits": [], "sites_removed": [], "sites_skipped": [],
        "clusters": [], "n_rounds": 0, "warnings": [],
    }
    if not seq:
        result["warnings"].append("No sequence loaded.")
        return result
    if not forward:
        result["warnings"].append("No valid enzymes selected to scrub.")
        return result
    allowed = frozenset(forward.keys())
    all_forbidden = _scrub_expand_forbidden(forward)
    feats = [f for f in (feats or [])
             if f.get("start") is not None and f.get("end") is not None]

    working = list(seq)
    failed: set = set()
    # Each successful cure removes >=1 target and (guard #2) adds none, so the
    # target count strictly drops; failed targets are parked. The cap is a
    # belt-and-braces backstop against a pathological cascade.
    initial = _scrub_scan_targets(seq, allowed, circular)
    max_iter = 2 * len(initial) + 8
    it = 0
    # Rescan the whole plasmid only when the sequence actually CHANGED. A
    # skipped (un-curable) site leaves `working` untouched, so the next round's
    # scan would return the identical target set minus the just-failed one —
    # re-running the full `_scan_restriction_sites` there was pure O(n)-per-skip
    # waste. `targets is None` means "sequence changed, rescan"; after a skip we
    # reuse the current scan and just drop the failed head.
    targets: "list | None" = None
    cur = "".join(working)   # bound before the loop; recomputed on each rescan
    while it < max_iter:
        it += 1
        if targets is None:
            cur = "".join(working)
            targets = [t for t in _scrub_scan_targets(cur, allowed, circular)
                       if (t["enzyme"], t["rec_start"], t["strand"]) not in failed]
        if not targets:
            break
        t = targets[0]
        ident = (t["enzyme"], t["rec_start"], t["strand"])
        changes, region, reason = _scrub_one_site(
            cur, t, feats, forward, all_forbidden, n, codon_frac=codon_frac)
        if changes is None:
            result["sites_skipped"].append({
                "enzyme": t["enzyme"], "pos": t["rec_start"],
                "strand": t["strand"], "region": region, "reason": reason})
            failed.add(ident)
            targets = targets[1:]   # seq unchanged → reuse scan, drop failed head
            continue
        for g, nb in sorted(changes.items()):
            result["edits"].append({
                "pos": g, "frm": cur[g], "to": nb,
                "region": region, "enzyme": t["enzyme"]})
            working[g] = nb
        result["sites_removed"].append({
            "enzyme": t["enzyme"], "pos": t["rec_start"],
            "strand": t["strand"], "region": region})
        targets = None   # sequence changed → rescan on the next round

    cured = "".join(working)

    # Substitution-only invariant — should be impossible to violate, but a
    # length drift would corrupt every downstream feature coordinate, so
    # abort loudly rather than ship it.
    if len(cured) != n:
        _log.error("Scrub: length drift %d->%d — aborting cure", n, len(cured))
        return {**result, "ok": False, "cured_seq": seq, "edits": [],
                "sites_removed": [],
                "warnings": result["warnings"] + ["Internal error: aborted "
                "to protect sequence integrity."]}
    result["cured_seq"] = cured

    # Authoritative wrap-aware guard: by the loop's exit condition the only
    # residual selected-enzyme sites are the ones we explicitly parked as
    # failed. Anything else means a cure cascade re-formed or spawned a site
    # — never observed (guard #2 forbids it), but if it ever happens we must
    # NOT claim the plasmid is clean.
    for t in _scrub_scan_targets(cured, allowed, circular):
        ident = (t["enzyme"], t["rec_start"], t["strand"])
        if ident in failed:
            continue
        _log.error("Scrub: unexpected residual %s at %d after cure pass",
                   t["enzyme"], t["rec_start"])
        result["sites_removed"] = [
            r for r in result["sites_removed"]
            if (r["enzyme"], r["pos"], r["strand"]) != ident]
        result["sites_skipped"].append({
            "enzyme": t["enzyme"], "pos": t["rec_start"],
            "strand": t["strand"], "region": "?",
            "reason": "could not be removed without side effects"})

    # Homopolymer report. The scorer PREFERS a cure that doesn't lengthen a
    # run, but on a site where every silent option extends one it still has to
    # pick something — and that is precisely the case a reader must be told
    # about, because a long run is a synthesis/sequencing hazard and it can
    # land twice when the site sits in a tandem repeat (2026-08-26 field
    # report). Measured on the whole molecule so the answer is the same one a
    # vendor's pre-flight would give, not a window-local approximation.
    def _whole_run(s: str) -> int:
        # A run can cross the origin on a circular molecule, so measure on a
        # doubled head — capped at `n`, or an all-one-base plasmid would
        # report a run longer than itself.
        probe = s + s[:_SCRUB_HOMOPOLYMER_FLANK] if circular else s
        return min(_scrub_longest_run(probe), n)

    hp_before = _whole_run(seq)
    hp_after = _whole_run(cured)
    result["max_homopolymer_run"] = hp_after
    if result["edits"] and hp_after > max(hp_before, _SCRUB_HOMOPOLYMER_FREE):
        result["warnings"].append(
            f"The cure lengthens a homopolymer run to {hp_after}×"
            f" (was {hp_before}×) — no silent alternative avoided it. Long "
            "runs are synthesis- and sequencing-hard; if this site sits in a "
            "repeated element the run is introduced once per copy.")

    clusters = _scrub_cluster_edits([e["pos"] for e in result["edits"]], n)
    result["clusters"] = [{"positions": c} for c in clusters]
    result["n_rounds"] = len(clusters)
    if n > 8000:
        result["warnings"].append(
            f"Plasmid is {n:,} bp — whole-plasmid QuikChange amplifies "
            "linearly and less efficiently above ~8 kb; use a long extension "
            "time (~1 min/kb) and more template.")
    return result


# ── Scrub: improved-QuikChange whole-plasmid primer design ──────────────────
#
# Per scrubbed locus (cluster of nearby cures), design ONE partial-overlap
# QuikChange primer pair (Liu & Naismith 2008) that amplifies the whole
# cured plasmid: a shared central region carries the cure(s) with ~12 bp
# correct flank each side, and each primer has a unique 3' extension (a
# perfectly-matched anchor) so the pair doesn't fully self-prime (far less
# primer-dimer than classic full-overlap QuikChange) while the products
# still self-anneal into a nicked circle the cell seals after transform —
# no ligase, no assembly. Primers are sliced from the CURED template, so
# what we display IS what anneals on the product (primer design is
# catastrophic-class — binding == display by construction). We compute Tm
# on circular slices directly, so an origin-adjacent locus needs no Primer3
# rotation ([PIT-06] is only for Primer3's linear placement engine).

_SCRUB_QC_FLANK_RANGE = range(10, 19)    # correct bases each side of the cure
_SCRUB_QC_EXT_RANGE   = range(8, 17)     # unique 3' anchor length (improved)
_SCRUB_QC_LEN_MIN, _SCRUB_QC_LEN_MAX = 25, 48
_SCRUB_QC_MIN_TEMPLATE = 60              # whole-plasmid PCR needs room
_SCRUB_QC_TARGET_TM = 72.0              # NN target for the matched flanks


def _circ_extract(seq: str, start: int, length: int, n: int) -> str:
    """Extract `length` bases of `seq` from `start` (mod n), wrapping the
    origin. `start` may be negative or >= n — it's normalized. Lets a
    primer straddle the origin so an origin-adjacent cure is still primed
    correctly without a Primer3 rotation."""
    if n <= 0 or length <= 0:
        return ""
    start %= n
    if start + length <= n:
        return seq[start:start + length]
    return "".join(seq[(start + i) % n] for i in range(length))


def _scrub_cluster_span(positions: list, n: int) -> tuple:
    """Smallest circular arc containing all `positions`: returns
    `(start, end)` where `end < start` signals the arc wraps the origin.
    Found as the complement of the largest inter-position gap."""
    ps = sorted(set(positions))
    if len(ps) == 1:
        return ps[0], ps[0]
    best_gap = -1
    best_i = 0
    for i in range(len(ps)):
        gap = (ps[(i + 1) % len(ps)] - ps[i]) % n
        if gap > best_gap:
            best_gap = gap
            best_i = i
    return ps[(best_i + 1) % len(ps)], ps[best_i]


def _scrub_qc_tm(primer: str, n_mismatch: int) -> float:
    """Stratagene QuikChange Tm: 81.5 + 0.41·%GC − 675/N − %mismatch. Unlike
    the nearest-neighbour `_mut_tm` (which assumes a perfectly-matched
    primer), this folds in the internal cure mismatch that lowers the
    first-cycle annealing temperature — so the reported number is the one
    the QuikChange ≥78 °C guideline actually refers to."""
    N = len(primer)
    if N == 0:
        return 0.0
    return 81.5 + 0.41 * _mut_gc_pct(primer) - 675.0 / N \
        - (n_mismatch / N) * 100.0


def _scrub_qc_primers(cured_seq: str, positions: list, *,
                      circular: bool = True, overlap: str = "improved",
                      round_no: int = 1) -> dict:
    """Design the improved-QuikChange primer pair for ONE cluster of cure
    positions on `cured_seq`. `overlap` is "improved" (partial overlap,
    default) or "classic" (full overlap). Returns a result dict with
    fwd_seq/rev_seq + their template coords, Tm/GC, overlap length, mismatch
    count and any warnings — or one carrying an `error` string when no
    acceptable pair fits.

    `hairpin_dg` / `homodimer_dg` are in **kcal/mol** (`dg_units` says so),
    reporting the WORSE of the two primers so the number can be used as a
    gate; ``None`` means primer3 couldn't measure it, which is warned about
    rather than passed off as a clean result."""
    n = len(cured_seq)
    res: dict = {"round": round_no, "positions": sorted(positions),
                 "warnings": []}
    if n < _SCRUB_QC_MIN_TEMPLATE:
        res["error"] = (f"Plasmid is only {n} bp — too small for "
                        "whole-plasmid QuikChange.")
        return res
    if not positions:
        res["error"] = "No cure positions in this round."
        return res

    start, end = _scrub_cluster_span(positions, n)
    width = (end - start) % n          # 0 for a single-base cluster
    ext_choices = [0] if overlap == "classic" else list(_SCRUB_QC_EXT_RANGE)
    posset = set(positions)

    def _count_mm(fp_start: int, fp_len: int) -> int:
        return sum(1 for i in range(fp_len) if (fp_start + i) % n in posset)

    # Winner stored as a positional tuple (not a dict) so each field keeps
    # its own type — a heterogeneous dict would widen the int coords to
    # `str | float` and trip the downstream int-typed helpers.
    best_score: "float | None" = None
    best: "tuple | None" = None
    for flank in _SCRUB_QC_FLANK_RANGE:
        ov_start = (start - flank) % n
        ov_len = width + 1 + 2 * flank
        for ext_f in ext_choices:
            fwd_len = ov_len + ext_f
            if not (_SCRUB_QC_LEN_MIN <= fwd_len <= _SCRUB_QC_LEN_MAX) \
                    or fwd_len >= n:
                continue
            fwd = _circ_extract(cured_seq, ov_start, fwd_len, n)
            gc_f = _mut_gc_pct(fwd)
            if not (35 <= gc_f <= 68):
                continue
            tm_f = _mut_tm(fwd)
            for ext_r in ext_choices:
                rev_start = (start - flank - ext_r) % n
                rev_len = ov_len + ext_r
                if not (_SCRUB_QC_LEN_MIN <= rev_len <= _SCRUB_QC_LEN_MAX) \
                        or rev_len >= n:
                    continue
                rev = _mut_revcomp(
                    _circ_extract(cured_seq, rev_start, rev_len, n))
                gc_r = _mut_gc_pct(rev)
                if not (35 <= gc_r <= 68):
                    continue
                tm_r = _mut_tm(rev)
                score = (
                    abs(tm_f - _SCRUB_QC_TARGET_TM)
                    + abs(tm_r - _SCRUB_QC_TARGET_TM)
                    + (0 if _mut_ends_gc(fwd) else 6)
                    + (0 if _mut_ends_gc(rev) else 6)
                    + abs(gc_f - 50) * 0.1 + abs(gc_r - 50) * 0.1
                    + abs(ov_len - 21) * 0.5
                    + (fwd_len + rev_len) * 0.02
                    + abs(tm_f - tm_r) * 0.5
                )
                if best_score is None or score < best_score:
                    best_score = score
                    best = (fwd, rev, ov_start, fwd_len, rev_start, rev_len,
                            tm_f, tm_r, gc_f, gc_r, ov_len)
    if best is None:
        res["error"] = ("No QuikChange primer pair met the length/GC "
                        "constraints for this locus.")
        return res

    (fwd, rev, fwd_start, fwd_len, rev_start, rev_len,
     tm_f, tm_r, gc_f, gc_r, ov_len) = best
    mm_f = _count_mm(fwd_start, fwd_len)
    mm_r = _count_mm(rev_start, rev_len)
    # kcal/mol. Worst of the pair: the most stable hairpin (most negative)
    # and the least stable self-dimer is NOT what we want — take the most
    # negative of each. `None` (not measured) drops out of the comparison
    # rather than winning it as a 0.0.
    hp_vals = [v for v in (_mut_hairpin_dg(fwd), _mut_hairpin_dg(rev))
               if v is not None]
    dim_vals = [v for v in (_mut_homodimer_dg(fwd), _mut_homodimer_dg(rev))
                if v is not None]
    hp = min(hp_vals) if hp_vals else None
    dim = min(dim_vals) if dim_vals else None
    res.update({
        "fwd_seq": fwd, "rev_seq": rev,
        "fwd_start": fwd_start, "fwd_len": fwd_len, "fwd_strand": 1,
        "rev_start": rev_start, "rev_len": rev_len, "rev_strand": -1,
        "fwd_tm": round(tm_f, 1), "rev_tm": round(tm_r, 1),
        "fwd_tm_qc": round(_scrub_qc_tm(fwd, mm_f), 1),
        "rev_tm_qc": round(_scrub_qc_tm(rev, mm_r), 1),
        "fwd_gc": round(gc_f, 1), "rev_gc": round(gc_r, 1),
        "overlap_len": ov_len, "n_mismatch": mm_f,
        # ΔG in kcal/mol — the WORST (most negative) of the two primers, so
        # the field can be used as a gate. `None` means primer3 couldn't
        # measure it, never "no structure".
        "hairpin_dg": None if hp is None else round(hp, 2),
        "homodimer_dg": None if dim is None else round(dim, 2),
        "dg_units": "kcal/mol",
        "overlap_style": "classic" if overlap == "classic" else "improved",
    })
    if not _mut_ends_gc(fwd) or not _mut_ends_gc(rev):
        res["warnings"].append("No 3' G/C clamp on a primer.")
    if min(res["fwd_tm_qc"], res["rev_tm_qc"]) < 78:
        res["warnings"].append(
            f"QuikChange Tm below the 78 °C guideline "
            f"(min {min(res['fwd_tm_qc'], res['rev_tm_qc'])} °C).")
    if (hp is not None and hp < -9.0) or (dim is not None and dim < -9.0):
        res["warnings"].append("Strong predicted hairpin/dimer.")
    if hp is None or dim is None:
        # Say so rather than letting a missing measurement read as a clean
        # one — an absent hairpin/dimer number is not a good hairpin/dimer
        # number.
        res["warnings"].append(
            "Hairpin/self-dimer ΔG could not be computed (primer3 "
            "unavailable) — these primers are NOT secondary-structure "
            "checked.")
    return res


def _scrub_qc_verify(orig: str, cured: str, rounds: list, n: int) -> tuple:
    """Prove the improved-QuikChange primers reconstitute the CURED plasmid
    EXACTLY — the seamless-cure guarantee for the QuikChange route (the Golden
    Braid route has `_scrub_gb_verify`).

    The lab product: each primer anneals to the ORIGINAL template and is
    extended around the whole plasmid; the product strand is the primer's own
    bases over its footprint (the cure rides in the primer) + a faithful copy
    of the template everywhere else; the two strands self-anneal through the
    shared overlap into the cured nicked circle. So simulate it: start from
    `orig`, lay every round's forward primer AND the top-strand bases of its
    reverse primer onto their footprints, and the result MUST equal `cured`.
    A cure that fell outside every primer's footprint would leave the original
    base there → product != cured → caught. Returns `(ok, errors)`."""
    if not orig or len(orig) != len(cured) or n <= 0:
        return False, ["sequence / length mismatch"]
    product = list(orig)
    any_primer = False
    for r in rounds:
        if r.get("error"):
            continue
        any_primer = True
        fwd = str(r.get("fwd_seq", ""))
        fstart = int(r.get("fwd_start", 0))
        for j, base in enumerate(fwd):
            product[(fstart + j) % n] = base
        rev_top = _mut_revcomp(str(r.get("rev_seq", "")))
        rstart = int(r.get("rev_start", 0))
        for j, base in enumerate(rev_top):
            product[(rstart + j) % n] = base
    if not any_primer:
        return False, ["no usable primer rounds to verify"]
    if "".join(product) != cured:
        return False, ["QuikChange primers would not reconstitute the cured "
                       "plasmid exactly — a cure fell outside primer reach"]
    return True, []


# ══════════════════════════════════════════════════════════════════════
#  Allele-specific PCR (ARMS) primer design
# ══════════════════════════════════════════════════════════════════════
#
# Allele-specific PCR IS the act of pinning the primer's 3'-terminal base on
# the variant: a matched 3' end extends, a mismatched one does not. None of
# the existing designers can express that. `_design_detection_primers` and
# `_design_generic_primers` pick whatever window scores best and slide freely;
# `_mut_design_fwd_anneal` / `_mut_design_rev_anneal` anchor the FIVE-prime end
# (they build Golden-Gate SDM primers whose 5' tail carries the site), so they
# are the wrong instrument entirely — sliding a 5'-anchored primer never fixes
# which base lands last.
#
# So the 3' terminus is fixed first and the length swept around it, which is
# the opposite of every other designer here.

_AS_MIN_LEN = 18
_AS_MAX_LEN = 27
# ARMS convention: a deliberate second mismatch near the 3' end further
# destabilises the mismatched duplex, so the discriminating primer fails to
# extend more cleanly. Position 3 from the 3' end is the usual choice — close
# enough to matter to the polymerase, far enough not to abolish the matched
# duplex too. OPT-IN: it changes the primer you order, and a caller who did
# not ask for a mismatch should not silently receive one.
_AS_DESTAB_OFFSET = 3
_AS_STRONG_CLAMP = frozenset("GC")


def _as_complement(base: str) -> str:
    return {"A": "T", "T": "A", "G": "C", "C": "G"}.get(base.upper(), "N")


def _design_allele_specific_primer(
        template: str, variant_pos: int, *,
        alt_base: "str | None" = None,
        orientation: str = "auto",
        target_tm: float = 60.0,
        destabilize: bool = False,
        circular: bool = False,
        tm_profile: str = _TM_PROFILE_DEFAULT) -> "dict | None":
    """One allele-specific primer whose 3'-terminal base sits ON `variant_pos`.

    `alt_base` is the allele the primer must be specific FOR; omit it to
    target whatever the template already carries at that position.
    `orientation` is ``"fwd"``, ``"rev"`` or ``"auto"`` (try both, keep the
    better-scoring one). Lengths `_AS_MIN_LEN`..`_AS_MAX_LEN` are swept with
    the 3' end PINNED — the length grows 5'-wards, never past the variant.

    Returns ``{seq, orientation, length, tm, gc, score, three_prime_base,
    gc_clamp, clamp_note, destabilized, destabilize_pos, binding_start,
    binding_end, matches_template}`` or None when the template is too short on
    the needed side.

    ``gc_clamp`` is False when the discriminating base is A or T, and
    ``clamp_note`` SAYS SO in words. That is not a scoring detail: an A/T 3'
    terminus cannot form a strong clamp no matter how the rest of the primer
    is chosen, and quietly ranking it lower would leave the caller wondering
    why every candidate scored badly instead of telling them the allele itself
    is the constraint.

    `tm_profile` (`_TM_PROFILES`) puts `target_tm` and the reported `tm` on
    that profile's scale; the default keeps the mutagenesis engine's own
    `_mut_tm`, so a default design is unchanged to the last digit.
    """
    t = (template or "").upper()
    n = len(t)
    if n == 0 or not (0 <= variant_pos < n):
        return None
    _cond = _tm_profile_kwargs(tm_profile)           # ValueError if unknown

    def _tm_profiled(s: str) -> float:
        try:
            return _p3_tm(s, **_cond)
        except Exception:
            # Same fallback ladder as `_mut_tm`: the weakest member of a
            # degenerate window, else the 2+4 rule.
            br = _degenerate_tm_bracket(s, **_cond)
            return br[0] if br is not None else float(_mut_tm(s))
    # Distinct names, not a redefinition: pyright's reportRedeclaration.
    tm_fn: "_Callable[[str], float]" = (
        _mut_tm if _tm_profile_resolve(tm_profile) == _TM_PROFILE_DEFAULT
        else _tm_profiled)
    ref_base = t[variant_pos]
    allele = (alt_base or ref_base).upper()
    if allele not in "ACGT":
        return None
    want = ("fwd", "rev") if orientation == "auto" else (orientation,)
    best = None
    for orient in want:
        for length in range(_AS_MIN_LEN, _AS_MAX_LEN + 1):
            if orient == "fwd":
                # 3' end AT the variant; the primer reads 5'->3' up to it.
                start = variant_pos - length + 1
                if start < 0:
                    if not circular:
                        continue
                    seq = "".join(t[(start + i) % n] for i in range(length))
                else:
                    seq = t[start:variant_pos + 1]
                seq = seq[:-1] + allele
                b_start, b_end = (start % n), (variant_pos + 1)
            else:
                # Reverse primer: its 3' end anneals to the variant, so it is
                # the reverse complement of the window STARTING at the variant.
                end = variant_pos + length
                if end > n:
                    if not circular:
                        continue
                    win = "".join(t[(variant_pos + i) % n]
                                  for i in range(length))
                else:
                    win = t[variant_pos:end]
                win = allele + win[1:]
                seq = _mut_revcomp(win)
                b_start, b_end = variant_pos, (end % n if circular else end)
            if len(seq) < _AS_MIN_LEN:
                continue
            destab_pos = None
            if destabilize and len(seq) > _AS_DESTAB_OFFSET:
                i = len(seq) - _AS_DESTAB_OFFSET
                cur = seq[i]
                # Any base other than the template's is a mismatch; pick a
                # transversion, which destabilises more than a transition.
                swap = {"A": "C", "C": "A", "G": "T", "T": "G"}.get(cur, cur)
                seq = seq[:i] + swap + seq[i + 1:]
                destab_pos = i
            score = _mut_score_outer(seq, target_tm, tm_fn=tm_fn)
            if best is None or score < best["score"]:
                three = seq[-1]
                best = {
                    "seq": seq,
                    "orientation": orient,
                    "length": len(seq),
                    "tm": tm_fn(seq),
                    "gc": _mut_gc_pct(seq),
                    "score": score,
                    "three_prime_base": three,
                    "gc_clamp": three in _AS_STRONG_CLAMP,
                    "clamp_note": (
                        None if three in _AS_STRONG_CLAMP else
                        f"the discriminating base is {three} — an A/T 3' "
                        f"terminus cannot form a GC clamp, so this primer "
                        f"relies on the mismatch alone for specificity. "
                        f"Consider a destabilizing second mismatch "
                        f"(destabilize: true) or a touchdown/higher "
                        f"annealing temperature."),
                    "destabilized": destab_pos is not None,
                    "destabilize_pos": destab_pos,
                    "binding_start": b_start,
                    "binding_end": b_end,
                    "matches_template": allele == ref_base,
                    "ref_base": ref_base,
                    "allele": allele,
                    "variant_pos": variant_pos,
                }
    return best
