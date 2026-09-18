"""CRISPR guide design — PAM scanning, guide triage, off-target search within
sequences you actually hold, and the cloning oligos.

Layer L1: imports biology L0 (`_rc` / `_iupac_pattern` / `_iter_match_starts`)
+ util/logging L0. No app, no I/O, no network.

WHAT THIS DOES NOT DO — read before trusting it
-----------------------------------------------
There is **no genome-wide off-target prediction here**, and none is implied.
Doing that honestly needs an index of the target genome, which this program does
not have and will not silently fake. `_guide_offtargets` searches the sequences
the CALLER hands it — the plasmid, a host genome the user actually loaded — and
says so in its own output (`searched_bp`). A guide reported with "no
off-targets" means "none in what you gave me", which is a different and much
weaker claim than "none in the genome".

Likewise `_score_guide` does not emit a trained efficiency score. Published
on-target models (Doench/Rule Set 2 and successors) are regressions over large
screen datasets; reproducing one from memory would give a number that LOOKS like
an efficiency prediction and is not. Instead this reports NAMED, checkable
properties — GC window, Pol III terminator, homopolymer runs, self-complementarity
— each traceable to a concrete mechanism, and a coarse tier derived from them.
A named flag a user can reason about beats a fabricated 0-100.

Same discipline as `splicecraft_splice`'s module docstring: say what the model
is, and refuse to dress a heuristic up as a measurement.

Geometry
--------
Everything is reported in FORWARD-strand coordinates with an explicit `strand`,
matching sacred invariant #2 (a reverse-strand hit stores the forward
coordinate). Reverse-strand guides are found by scanning the forward sequence
for the reverse complement of the PAM rather than by building an RC of the whole
sequence — one less coordinate transform to get wrong.

Circular plasmids are scanned across the origin the way `_scan_restriction_sites`
does (sacred #6): the search runs over `seq + seq[:pad]`, and a guide whose span
crosses bp 0 is reported with `wraps=True`.
"""

from __future__ import annotations

from splicecraft_biology import _rc, _iupac_pattern, _iter_match_starts
from splicecraft_logging import _log

__all__ = [
    "_CAS_VARIANTS", "_cas_variant", "_cas_variant_names",
    "_find_guides", "_score_guide", "_guide_offtargets",
    "_guide_cloning_oligos", "_design_guides",
]


# ── Cas variants ──────────────────────────────────────────────────────────────
# Three DISTINCT PAM geometries rather than a long catalogue of near-duplicates:
# a short 3' PAM, a long degenerate 3' PAM, and a 5' PAM with a staggered cut.
# `cut_offsets` are measured from the PROTOSPACER 5' end on the guide's own
# strand, as (top_strand_nick, bottom_strand_nick); equal values = blunt.
_CAS_VARIANTS: "dict[str, dict]" = {
    "spcas9": {
        "label":       "SpCas9",
        "pam":         "NGG",
        "pam_side":    "3prime",
        "guide_len":   20,
        "cut_offsets": (17, 17),          # blunt, 3 bp 5' of the PAM
        "notes":       "S. pyogenes Cas9. The default for most editing and "
                       "the PAM most vectors are built around.",
    },
    "sacas9": {
        "label":       "SaCas9",
        "pam":         "NNGRRT",
        "pam_side":    "3prime",
        "guide_len":   21,
        "cut_offsets": (18, 18),          # blunt
        "notes":       "S. aureus Cas9. Small enough for single-AAV "
                       "delivery; the longer PAM means far fewer sites.",
    },
    "lbcas12a": {
        "label":       "LbCas12a (Cpf1)",
        "pam":         "TTTV",
        "pam_side":    "5prime",
        "guide_len":   23,
        "cut_offsets": (18, 23),          # staggered, 5-nt 5' overhang
        "notes":       "PAM sits 5' of the spacer and the cut is staggered, "
                       "leaving a 5' overhang rather than a blunt end.",
    },
}

# Pol III (U6 / H1) terminates on a run of T, so a guide carrying one is
# truncated at transcription — the single most common reason a cloned guide
# does nothing. This is a mechanism, not a correlation, which is why it is a
# hard flag rather than a score penalty.
_POLIII_TERMINATOR = "TTTT"
_MAX_HOMOPOLYMER = 4          # runs LONGER than this are flagged
_GC_LOW, _GC_HIGH = 40.0, 70.0
_SEED_LEN = 12                # PAM-proximal bases where a mismatch hurts most
_OFFTARGET_MAX_MISMATCH = 4
_MAX_GUIDES = 5000            # bound on a scan of a large record


def _cas_variant_names() -> "list[str]":
    """Variant keys, stable order, for a UI picker or an agent enum."""
    return list(_CAS_VARIANTS)


def _cas_variant(name: object) -> dict:
    """Look up a Cas variant. Raises ValueError naming the alternatives, so a
    typo in an agent payload gets an actionable 400 rather than a silent
    fallback to SpCas9 — quietly designing for the wrong nuclease would be a
    wasted cloning week."""
    key = str(name or "").strip().lower().replace("-", "").replace(" ", "")
    if key in _CAS_VARIANTS:
        return _CAS_VARIANTS[key]
    raise ValueError(
        f"unknown Cas variant {name!r}; choose one of "
        f"{', '.join(sorted(_CAS_VARIANTS))}"
    )


def _take_circular(seq: str, start: int, length: int) -> str:
    """`length` bases from `start`, wrapping the origin when needed.

    NOT `splicecraft_biology._slice_circular`, which takes an END boundary and
    cannot repeat: this one takes a LENGTH and tolerates `length > len(seq)`,
    which a PAM-padded haystack on a very small plasmid needs. Two different
    contracts, so two names — sharing one would make every call site ambiguous
    about whether the third argument is a position or a count.

    Returns "" for a non-positive length or an empty sequence rather than
    raising: callers are scan loops, and one degenerate candidate should not
    abort a whole plasmid's worth of good ones.
    """
    n = len(seq)
    if n == 0 or length <= 0:
        return ""
    start %= n
    if start + length <= n:
        return seq[start:start + length]
    # Wrap. Repeat only as far as needed rather than doubling a 5 Mb string.
    out = [seq[start:]]
    remaining = length - (n - start)
    while remaining > 0:
        take = min(remaining, n)
        out.append(seq[:take])
        remaining -= take
    return "".join(out)


def _gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for c in seq.upper() if c in "GCS")
    return gc * 100.0 / len(seq)


def _longest_homopolymer(seq: str) -> int:
    best = run = 0
    prev = ""
    for c in seq.upper():
        run = run + 1 if c == prev else 1
        prev = c
        best = max(best, run)
    return best


def _self_complement_len(seq: str) -> int:
    """Longest suffix of `seq` that is the reverse complement of a prefix —
    i.e. how far the guide can fold back on itself.

    A hairpin in the spacer competes with the scaffold fold and is a real
    (mechanistic) cause of a dead guide. Reported in bases so the caller can
    set its own tolerance instead of trusting a hidden threshold.
    """
    s = seq.upper()
    n = len(s)
    if n < 6:
        return 0
    rc = _rc(s)
    best = 0
    # A stem longer than half the guide cannot form without overlapping itself.
    for k in range(3, n // 2 + 1):
        if s[:k] == rc[:k]:
            best = k
    return best


def _find_guides(sequence: str, *, variant: str = "spcas9",
                 circular: bool = False,
                 max_guides: int = _MAX_GUIDES) -> "list[dict]":
    """Every protospacer+PAM site in `sequence`, both strands.

    Each hit: ``{guide, pam, strand, start, end, pam_start, pam_end, cut_site,
    wraps, variant}`` — all coordinates 0-based, half-open, in the FORWARD
    frame whatever the strand (sacred #2). `cut_site` is the top-strand nick
    position in the forward frame.

    Reverse-strand sites are found by searching for the reverse complement of
    the PAM on the forward strand, so no coordinate is ever mapped through an
    RC of the whole sequence. Circular records are scanned across the origin
    (sacred #6) and a site spanning bp 0 is flagged `wraps`.
    """
    seq = (sequence or "").upper()
    n = len(seq)
    v = _cas_variant(variant)
    glen, pam = int(v["guide_len"]), str(v["pam"])
    plen = len(pam)
    need = glen + plen
    if n < need:
        return []
    cap = max(1, int(max_guides or _MAX_GUIDES))

    # Circular: extend by one full site minus a base so a PAM straddling the
    # origin is still found. Linear: no extension.
    pad = (need - 1) if circular else 0
    haystack = seq + _take_circular(seq, 0, pad) if pad else seq

    fwd_pat = _iupac_pattern(pam)
    rev_pat = _iupac_pattern(_rc(pam))
    three_prime = v["pam_side"] == "3prime"
    top_off, _bot_off = v["cut_offsets"]

    out: list[dict] = []
    seen: set = set()

    def _emit(gstart: int, pstart: int, strand: int) -> None:
        if len(out) >= cap:
            return
        guide = _take_circular(seq, gstart, glen)
        if strand == -1:
            guide = _rc(guide)
        if len(guide) != glen or any(c not in "ACGT" for c in guide):
            # An N inside the protospacer means the template is ambiguous
            # there; a guide ordered against it would be wrong. Skip rather
            # than emit a guide with an N in it.
            return
        pam_seq = _take_circular(seq, pstart, plen)
        if strand == -1:
            pam_seq = _rc(pam_seq)
        # Cut site, forward frame. On the minus strand the protospacer runs
        # right-to-left, so the offset counts back from its forward END.
        if strand == 1:
            cut = (gstart + top_off) % n
        else:
            cut = (gstart + glen - top_off) % n
        key = (gstart % n, strand)
        if key in seen:
            return
        seen.add(key)
        out.append({
            "guide":     guide,
            "pam":       pam_seq,
            "strand":    strand,
            "start":     gstart % n,
            "end":       (gstart + glen) % n,
            "pam_start": pstart % n,
            "pam_end":   (pstart + plen) % n,
            "cut_site":  cut,
            "wraps":     (gstart % n) + glen > n or (pstart % n) + plen > n,
            "variant":   str(variant).lower(),
        })

    for pat, strand in ((fwd_pat, 1), (rev_pat, -1)):
        for hit in _iter_match_starts(pat, haystack):
            if hit >= n:                      # only origin-anchored duplicates
                continue
            if strand == 1:
                # Plus strand: 3' PAM sits after the spacer, 5' PAM before it.
                gstart = hit - glen if three_prime else hit + plen
            else:
                # Minus strand: the forward view is reverse-complemented, so a
                # 3' PAM appears BEFORE the spacer in forward coordinates.
                gstart = hit + plen if three_prime else hit - glen
            if not circular and not (0 <= gstart and gstart + glen <= n):
                continue
            _emit(gstart, hit, strand)
        if len(out) >= cap:
            _log.info("crispr: guide scan hit the %d-site cap on a %d bp "
                      "record; results truncated", cap, n)
            break
    out.sort(key=lambda g: (g["start"], -g["strand"]))
    return out


def _score_guide(guide: str, *, u6_driven: bool = True) -> dict:
    """Triage one spacer on named, mechanistic properties.

    Returns ``{gc_pct, homopolymer, self_complement, starts_with_g, flags,
    tier}``. `flags` is a list of short human strings; `tier` is
    ``"good" | "ok" | "poor"``.

    NOT an efficiency prediction — see the module docstring. Each flag names a
    mechanism the user can check:

      * ``"Pol III terminator (TTTT)"`` — U6/H1 stops transcribing, so the
        guide is truncated before it is ever loaded. Fatal, not cosmetic.
      * ``"GC x%"`` outside 40-70 — duplex stability; too low binds poorly,
        too high tolerates mismatches (worse specificity).
      * ``"homopolymer run of N"`` — synthesis and sequencing both slip.
      * ``"self-complementary N bp"`` — the spacer folds instead of pairing.
      * ``"does not start with G"`` — U6 initiates best on G; advisory only,
        and suppressed when `u6_driven` is False (a T7 or tRNA-processed
        construct has no such preference).
    """
    g = (guide or "").upper()
    flags: list[str] = []
    gc = _gc_pct(g)
    homo = _longest_homopolymer(g)
    selfc = _self_complement_len(g)
    fatal = False

    if not g:
        return {"gc_pct": 0.0, "homopolymer": 0, "self_complement": 0,
                "starts_with_g": False, "flags": ["empty guide"],
                "tier": "poor"}
    if _POLIII_TERMINATOR in g:
        flags.append("Pol III terminator (TTTT)")
        fatal = True
    if gc < _GC_LOW:
        flags.append(f"GC {gc:.0f}% (low)")
    elif gc > _GC_HIGH:
        flags.append(f"GC {gc:.0f}% (high)")
    if homo > _MAX_HOMOPOLYMER:
        flags.append(f"homopolymer run of {homo}")
    if selfc >= 6:
        flags.append(f"self-complementary {selfc} bp")
    starts_g = g.startswith("G")
    if u6_driven and not starts_g:
        flags.append("does not start with G")

    if fatal:
        tier = "poor"
    elif not flags:
        tier = "good"
    elif len(flags) == 1 and flags[0] == "does not start with G":
        # A missing 5' G is routinely fixed by just adding one, so it alone
        # should not demote a guide that is otherwise clean.
        tier = "good"
    else:
        tier = "ok" if len(flags) == 1 else "poor"
    return {"gc_pct": round(gc, 1), "homopolymer": homo,
            "self_complement": selfc, "starts_with_g": starts_g,
            "flags": flags, "tier": tier}


def _guide_offtargets(guide: str, sequence: str, *,
                      variant: str = "spcas9", circular: bool = False,
                      max_mismatch: int = _OFFTARGET_MAX_MISMATCH,
                      require_pam: bool = True,
                      max_hits: int = 200) -> dict:
    """Near-matches of `guide` in `sequence`, both strands.

    Returns ``{searched_bp, max_mismatch, require_pam, hits, truncated}`` where
    each hit is ``{start, strand, mismatches, seed_mismatches, sequence, pam,
    pam_ok}``.

    `searched_bp` is in the result ON PURPOSE: this is a search of what the
    caller supplied, never a genome-wide claim, and the number is the honest
    scope of the answer. A caller that reports "no off-targets" without also
    reporting `searched_bp` is overstating it.

    `seed_mismatches` counts mismatches in the PAM-proximal `_SEED_LEN` bases,
    where Cas9 tolerates them far less — a hit with 3 mismatches all in the
    distal end is a real risk, the same 3 in the seed usually is not. Hits are
    ordered most-dangerous first (fewest seed mismatches, then fewest total).
    """
    g = (guide or "").upper()
    seq = (sequence or "").upper()
    n, glen = len(seq), len(g)
    if not g or n < glen:
        return {"searched_bp": n, "max_mismatch": int(max_mismatch),
                "require_pam": bool(require_pam), "hits": [], "truncated": False}
    v = _cas_variant(variant)
    pam, plen = str(v["pam"]), len(str(v["pam"]))
    three_prime = v["pam_side"] == "3prime"
    pam_re = _iupac_pattern(pam)
    limit = max(0, int(max_mismatch))
    cap = max(1, int(max_hits))
    # Seed = the PAM-proximal end of the spacer.
    seed_slice = (slice(glen - _SEED_LEN, glen) if three_prime
                  else slice(0, _SEED_LEN))

    hits: list[dict] = []
    truncated = False
    span = n if circular else n - glen + 1
    for strand in (1, -1):
        probe = g if strand == 1 else _rc(g)
        for start in range(span):
            window = _take_circular(seq, start, glen) if circular \
                else seq[start:start + glen]
            if len(window) != glen:
                continue
            mism = 0
            for a, b in zip(probe, window):
                if a != b:
                    mism += 1
                    if mism > limit:
                        break
            if mism > limit:
                continue
            # PAM check, on the same side the variant uses.
            if three_prime:
                pstart = start + glen if strand == 1 else start - plen
            else:
                pstart = start - plen if strand == 1 else start + glen
            pam_seq = ""
            if circular:
                pam_seq = _take_circular(seq, pstart, plen)
            elif 0 <= pstart and pstart + plen <= n:
                pam_seq = seq[pstart:pstart + plen]
            oriented = pam_seq if strand == 1 else _rc(pam_seq)
            pam_ok = bool(oriented) and bool(pam_re.fullmatch(oriented))
            if require_pam and not pam_ok:
                continue
            seed_probe = probe[seed_slice]
            seed_window = window[seed_slice]
            seed_mm = sum(1 for a, b in zip(seed_probe, seed_window) if a != b)
            hits.append({
                "start": start % n, "strand": strand, "mismatches": mism,
                "seed_mismatches": seed_mm,
                "sequence": window if strand == 1 else _rc(window),
                "pam": oriented, "pam_ok": pam_ok,
            })
            if len(hits) >= cap:
                truncated = True
                break
        if truncated:
            break
    hits.sort(key=lambda h: (h["seed_mismatches"], h["mismatches"], h["start"]))
    return {"searched_bp": n, "max_mismatch": limit,
            "require_pam": bool(require_pam), "hits": hits,
            "truncated": truncated}


# Type IIS cloning overhangs. Defaults are the BbsI pair from the
# pX330 / pSpCas9(BB) family that most lab guide vectors descend from;
# a BsmBI lentiCRISPR backbone wants different ones. Documented rather than
# hidden, because a wrong overhang is a cloning failure the user only discovers
# a week later.
_GUIDE_OVERHANG_TOP = "CACC"
_GUIDE_OVERHANG_BOTTOM = "AAAC"


def _guide_cloning_oligos(guide: str, *,
                          overhang_top: str = _GUIDE_OVERHANG_TOP,
                          overhang_bottom: str = _GUIDE_OVERHANG_BOTTOM,
                          add_g: bool = True) -> dict:
    """The annealed oligo pair that clones `guide` into a Type IIS guide vector.

    Returns ``{top, bottom, guide, duplex, added_g, overhang_top,
    overhang_bottom}``.

    The two oligos anneal over the DUPLEX CORE (`duplex`) and each contributes a
    single-stranded overhang at one end, which is what ligates into the cut
    vector::

        top     5'-CACC G[spacer]        -3'
        bottom  3'-     C[rc spacer]CAAA -5'

    So the bottom oligo is `overhang_bottom + rc(duplex)`, NOT
    `overhang_bottom + rc(guide)`: it has to carry the complement of the
    transcription-start G as well, or the two oligos come out different lengths
    and leave a one-base gap at that end. That missing C is the classic way this
    pair is written down wrong.

    `add_g`: U6 initiates best on G, so the duplex needs exactly ONE at its 5'
    end. A guide that already starts with G supplies its own and gets no extra
    (which is why the standard protocol says use `CACC` for a G-initiated
    spacer and `CACCG` otherwise) — `added_g` reports which happened instead of
    leaving the user to count bases.
    """
    g = (guide or "").upper()
    if not g or any(c not in "ACGT" for c in g):
        raise ValueError(
            "guide must be a non-empty unambiguous ACGT spacer; got "
            f"{guide!r}"
        )
    top_oh = (overhang_top or "").upper()
    bot_oh = (overhang_bottom or "").upper()
    for label, oh in (("overhang_top", top_oh), ("overhang_bottom", bot_oh)):
        if oh and any(c not in "ACGT" for c in oh):
            raise ValueError(f"{label} must be unambiguous ACGT; got {oh!r}")
    added = bool(add_g) and not g.startswith("G")
    duplex = ("G" + g) if added else g
    return {
        "top":             top_oh + duplex,
        "bottom":          bot_oh + _rc(duplex),
        "duplex":          duplex,
        "guide":           g,
        "added_g":         added,
        "overhang_top":    top_oh,
        "overhang_bottom": bot_oh,
    }


def _design_guides(sequence: str, *, variant: str = "spcas9",
                   circular: bool = False, u6_driven: bool = True,
                   region: "tuple[int, int] | None" = None,
                   offtarget_in: "str | None" = None,
                   max_mismatch: int = _OFFTARGET_MAX_MISMATCH,
                   limit: int = 50) -> dict:
    """Scan, triage, and (optionally) off-target-check in one call — the shape
    the UI and the agent endpoint both want.

    `region` restricts results to guides whose CUT falls inside a
    ``(start, end)`` half-open window, which is how you aim at an exon or a
    specific codon rather than the whole plasmid.

    `offtarget_in` is any sequence to search for near-matches — usually the
    plasmid itself, or a host genome the user has loaded. When omitted NO
    off-target claim is made at all, and `offtarget_searched_bp` is 0 so the
    caller cannot mistake silence for a clean result.
    """
    seq = (sequence or "").upper()
    guides = _find_guides(seq, variant=variant, circular=circular)
    if region:
        lo, hi = int(region[0]), int(region[1])
        n = len(seq) or 1
        if circular and lo > hi:
            # Region wraps the origin.
            guides = [g for g in guides
                      if g["cut_site"] >= lo % n or g["cut_site"] < hi % n]
        else:
            guides = [g for g in guides if lo <= g["cut_site"] < hi]
    for g in guides:
        g["score"] = _score_guide(g["guide"], u6_driven=u6_driven)
    searched = 0
    if offtarget_in:
        searched = len(offtarget_in)
        for g in guides:
            ot = _guide_offtargets(
                g["guide"], offtarget_in, variant=variant,
                circular=circular and offtarget_in == seq,
                max_mismatch=max_mismatch)
            # The guide's own site is a match against itself; a caller wants
            # OTHER sites, so drop the exact on-target hit.
            others = [h for h in ot["hits"] if h["mismatches"] > 0]
            g["offtargets"] = others
            g["n_offtargets"] = len(others)
    tier_rank = {"good": 0, "ok": 1, "poor": 2}
    guides.sort(key=lambda g: (
        tier_rank.get(g["score"]["tier"], 3),
        g.get("n_offtargets", 0),
        g["start"],
    ))
    lim = max(1, int(limit))
    return {
        "variant":  str(variant).lower(),
        "label":    _cas_variant(variant)["label"],
        "circular": bool(circular),
        "n_found":  len(guides),
        "guides":   guides[:lim],
        "truncated": len(guides) > lim,
        "offtarget_searched_bp": searched,
    }
