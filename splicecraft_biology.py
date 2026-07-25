"""splicecraft_biology — pure biology primitives extracted from
splicecraft.py as a controlled test of the single-file seam.

The single-file rule (entire app in splicecraft.py) is intentional;
see CLAUDE.md + docs/architecture.md. This module is the first
deliberate extraction, scoped to entities that:

1. Have NO `PlasmidApp` coupling (no `self.notify`, no Textual,
   no reactive attrs).
2. Are pure functions or top-level constants.
3. Are imported back into `splicecraft.py` so external callers
   keep `splicecraft._rc(...)` etc. unchanged.

If this stays clean (no cross-imports back into splicecraft.py, no
test churn), it's the precedent for future extractions. See
CONTRIBUTING.md's "three-test rule" for the criteria.

Sacred invariants this module owns:
  #3  — `_rc` handles full IUPAC via `_IUPAC_COMP` (not just ACGT).
  #4  — `_iupac_pattern` is bounded-LRU cached in `_PATTERN_CACHE`.
  #8  — `_feat_len` returns `(total - start) + end` when `end < start`
       so wrap features have the right length.

These are the same invariants documented in `CLAUDE.md`; the
extraction did not change them, only their physical home.
"""
from __future__ import annotations

import base64 as _base64
import functools
import gzip as _gzip
import math as _math
import re
import threading as _threading
from collections import OrderedDict

# Same-layer (L0) sibling deps for the restriction-scan engine
# (Phase D): the LRU caches + catalog/enzyme getters live in _state;
# the @_timed scan decorator + _log come from splicecraft_logging.
import splicecraft_state as _state
from splicecraft_logging import _log, _timed
from splicecraft_util import _normalize_dna_for_align  # same-layer L0 (util doesn't import biology); used by _search_subsequence


# ── IUPAC + reverse complement ────────────────────────────────────────────


_IUPAC_RE: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "W": "[AT]", "S": "[CG]",
    "M": "[AC]", "K": "[GT]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


# Pattern cache (sacred invariant #4). Bounded LRU so a long-lived
# process scanning many recognition sites can't grow the cache
# indefinitely; the catalog is ~120 enzymes (palindromic + RC variants
# < 256), so 256 is comfortably above steady state. Public dict so
# tests can `.clear()` and inspect membership.
_PATTERN_CACHE: "OrderedDict[str, re.Pattern[str]]" = OrderedDict()
_PATTERN_CACHE_MAX = 256

# Sweep #26 (2026-05-23): cut-position validation bounds for custom enzymes.
# Type IIS enzymes cut OUTSIDE the recognition site (offset can be larger than
# `len(site)`); the ±30 window is generous enough to accept every NEB Type IIS
# in production while still catching typos. Shared so the UI
# (`AddCustomEnzymeModal`) and the agent custom-enzyme validator read one
# definition rather than re-hardcoding ±30. Re-exported into the hub.
_ENZYME_CUT_RANGE = 30


def _iupac_pattern(site: str) -> "re.Pattern[str]":
    # Sweep #22: case-fold cache key. Patterns ARE always built from
    # `site.upper()` internally, so two calls with `"gaattc"` and
    # `"GAATTC"` produce identical regex objects but pre-fix occupied
    # two separate slots — wasting one cap unit per mixed-case
    # variant. Normalize to uppercase before lookup AND store.
    #
    # 2026-05-27 (audit-5 restriction M3): reject unknown characters
    # instead of silently letting them into the regex. Pre-fix any
    # char not in `_IUPAC_RE` fell through to `c` itself — a custom
    # enzyme site like ``"GAATU"`` (RNA U typo) compiled to a literal
    # ``U`` pattern that never matched DNA, and a stray regex
    # metacharacter (``*``, ``(``, ``?``) became part of the pattern.
    # User-defined enzyme sites are an attack surface; validate here.
    key = site.upper()
    # An empty / whitespace-only site would compile to `re.compile("")`, which
    # matches at EVERY position → a phantom cut at every base of any plasmid.
    # Both `_rebuild_scan_catalog` and `_enzyme_cuts_impl` wrap this call in an
    # `except ValueError … skip`, so raising here neutralises a corrupt
    # custom_enzymes.json entry (e.g. `["", 1, 1]`) at a single chokepoint.
    if not key.strip():
        raise ValueError(f"recognition site {site!r} is empty / whitespace-only")
    with _state._COMPUTE_CACHE_LOCK:
        pat = _PATTERN_CACHE.get(key)
        if pat is not None:
            _PATTERN_CACHE.move_to_end(key)
            return pat
    bad = [c for c in key if c not in _IUPAC_RE]
    if bad:
        raise ValueError(
            f"recognition site {site!r} contains non-IUPAC "
            f"character(s) {', '.join(repr(c) for c in bad[:6])}"
            f"{' (truncated)' if len(bad) > 6 else ''}"
        )
    pat = re.compile("".join(_IUPAC_RE[c] for c in key))
    with _state._COMPUTE_CACHE_LOCK:
        _PATTERN_CACHE[key] = pat
        if len(_PATTERN_CACHE) > _PATTERN_CACHE_MAX:
            _PATTERN_CACHE.popitem(last=False)
    return pat


_IUPAC_COMP = str.maketrans(
    # U→A so `_rc` and minus-strand translation are correct for RNA-form
    # bases the loaders accept (`_IUPAC_NUC_CHARS` includes U). Pre-fix a
    # `U` passed through `.translate()` unchanged, so a minus-strand CDS
    # over it mistranslated (the stray U survived into the codon → "?").
    "ACGTURYWSMKBDHVN",
    "TGCAAYRWSKMVHDBN",
)

# Case-preserving ACGT complement used by the sequence-panel renderer.
_DNA_COMP_PRESERVE_CASE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


# Tiny LRU on `_rc` — the cached value of a 200 kb plasmid's RC is
# itself ~200 kb, so a 4-entry cap is enough to cover the working set
# (current sequence + a couple of recent rotations / undo snapshots)
# without ballooning RAM. Benches show 35–113× speedup on cache hit
# for cosmid-size sequences (scripts/perf_probe.py).
@functools.lru_cache(maxsize=4)
def _rc(seq: str) -> str:
    return seq.upper().translate(_IUPAC_COMP)[::-1]


# ── Wrap-aware coordinate helpers ─────────────────────────────────────────


def _feat_len(start: int, end: int, total: int) -> int:
    """Circular-aware feature length. A wrap feature (end < start) is
    (total - start) + end bp long; a linear feature is end - start."""
    return (total - start) + end if end < start else end - start


def _seq_len(record) -> int:
    """Length of ``record.seq`` in bp, or 0 if the record has no
    sequence attached. BioPython's ``SeqRecord.seq`` is typed as
    ``Seq | MutableSeq | None`` because the dataclass allows records
    without sequences (rare — e.g. annotation-only GenBank views).
    SpliceCraft always loads records with sequence content, but
    routing length lookups through this helper sidesteps the
    ``"None" is not assignable to "Sized"`` pyright noise at every
    ``len(rec.seq)`` call site without an inline None guard."""
    seq = getattr(record, "seq", None)
    return len(seq) if seq is not None else 0


def _slice_circular(seq: str, start: int, end: int) -> str:
    """Circular-aware slice. If end > start this is a normal slice; if
    end < start the slice wraps the origin and returns seq[start:] + seq[:end].
    end == start is treated as empty (not "wrap whole plasmid") — callers
    that want the latter should pass explicit boundaries. Used by the
    primer-design helpers so a region straddling the origin can be
    primer-designed without special casing at every call site.
    """
    if end >= start:
        return seq[start:end]
    return seq[start:] + seq[:end]


# ── RNA secondary-structure free energy (Turner-2004, pure-Python) ──────────
#
# A dependency-free minimum-free-energy RNA folder + structure evaluator,
# reproducing the ViennaRNA dangles=2 model. Parameters are the standard
# Turner-2004 nearest-neighbor set, embedded as a gzip+base64 constant so
# nothing is read from disk and no compiled RNA library is required. Every
# energy term + the MFE folder are validated to match ViennaRNA exactly
# (see tests/test_rna_fold.py, which asserts against a frozen reference).
#
# Energies are centi-kcal/mol (int) internally; kcal/mol at the API edge.


_RNA_INF = 1 << 30
_RNA_MAXLOOP = 30                       # ViennaRNA default max interior-loop size
_RNA_FOLD_MAX_LEN = 600                 # O(n^3) DP — cap to bound worst-case time

# ViennaRNA pair-type order: CG=0 GC=1 GU=2 UG=3 AU=4 UA=5, no-pair=6.
_RNA_PT = {('C', 'G'): 0, ('G', 'C'): 1, ('G', 'U'): 2,
           ('U', 'G'): 3, ('A', 'U'): 4, ('U', 'A'): 5}
_RNA_BI = {'@': 0, 'A': 1, 'C': 2, 'G': 3, 'U': 4}    # 5-wide (mm/dangle)
_RNA_BI4 = {'A': 0, 'C': 1, 'G': 2, 'U': 3}           # 4-wide (int22)

# Turner-2004 parameters (ViennaRNA `params_save` dump, gzip+base64). The
# numeric values are the published experimental free energies; only the
# DEF placeholders are resolved (to a uniform -50, measured against the
# reference engine). Decoded + parsed once, lazily, on first fold/eval.
_RNA_TURNER_PARAMS_GZ_B64 = (
    "H4sIAAzFJGoC/+29X48kuZEn+Kz8FA7oYfe06xqS/v/hABE5M7ECdnXAYQncW6NWqhkVrrvU6G4dtM"
    "B++ItwM9KMThrdPTsjI7goAfJmBdOcNJI/N9L4o/G3v23+7z/Zf/vb939pfvz006cfPv/y+afm3758"
    "/7n5/8zv1cvLb5uff/n05//35Z9+1zSvl+b6v8vr+nS3p1t/sZC2t+cfmuZ3//TSNK3p1fXZdben0b"
    "en7ikdfnkJf9RB9rBmrE9jVHjR9ZcXkr5lN/om1vzzv/xreDf8cn2+0E9D+KOmWbPbcf1dKfiFvdWE"
    "Wvo/0mt6wV/Yn1Jx+Kb1j7AC1/S2ApDBKoBvvdXVt/F3n7/+8tdP3//45fPPZ5v7WonrS6/F97cn1H"
    "FYf1HD+uzxl5fwR/2y/umq59zd0hqEe/zlBd8Ef3QVgBLWnpnn9fcZ0y9YnH8T/NEC3bq26qQw/RLq"
    "BM0zz/RHy/rLOPM/DXXyxcGb4I+mCdNxBeK3wpPeem3uH778/MOnX/781+/++unLTz9++fpy6485dC"
    "b0u+/YeduZ/snG5a40jOeOjW0+qlFOsdcqRMkLDXEYjBMNzOvvL1SeZlXrWNX2pIc87OKq6US9FQ5N"
    "Y1R4tvDUmN5mAwbbjsr2f6sFaay5UhstctJdSN8E1rKbzGvh9xdeEZIe8U2HpfHL0LFvBUhnVQpfhj"
    "298eXdr2o1c7bVzNtaLaSPSuMXcZM+0GrRWGNa4Muj7C56sqppQfpYh6bPJftB4Z/y61/Bp2m4gc+n"
    "16eZ2Jdr0ZQRni/wV14ajSSkQbqB73yv6OWLYdJDRxKUvgr00Cwr5Adou4HK7lcQtetrIdunWTZKTE"
    "x6gj81Qfo6sBaSHg2ThnaeqMfW9LVh149ms6ztBfaYpwfUeyKTimnDsiXpridpaNqc9NyRBE8jiAz8"
    "tLY5pNePZqPR8rMvMlhYaEEDH029sCkGTTdQ2r9KBenGcGnDPlY83U7HW601v6bVBnOu1dpu5jXvtj"
    "VHvRloz7eal1gE6dOtZsz7jjUcLUP68tJY49+WL1+vc/Ovn76nr1n2mc9eJyhydirtZ8hxuiTdjurZ"
    "qwZdl3vms9FETMelYRBv0kVp/a1q36r2cVXLfVE20xUzb/FjihCFWclBBOMqbjz6ciximNWTV21a7d"
    "Og6Bl+edlkwFNPxezwDC/H2cChlw9s+fmtat+q9jFVy35c9Nc3TQu+Zb/FinzL/pb9v1m28FHZTFq+"
    "fUDekr1jCr5lf8v+3z47+30x3ZvAhRuNZ3wZmrkppiO+jO6Dvhq4h/CGqt3fUrQ9y2iOSPdvLNs8b9"
    "W+tdrztprwUfmVkxbcQzoD4pncGe1ywFfajv3HfF9arZa3Ve3+luLm0wnp/oA0UlnOl92OZ6vWjR9U"
    "tW+t9jGt9qaq8e/LD3///pcvL+l2P22UTcm2ZGBxQPqFyCGDitPIp/ASHZPumTSWN6hNGjfYkCDTqz"
    "gdXq6i8pJsTm7REQfnJcuc0SNjgmB5g4rTnvKADIaZ5m/huWWhwCwIBVCxkfawNumXDeMGpDF7EsoO"
    "6S0TBF8+8pfvZE8qrnlEStqt+ShIj4dqnpat1dtbDeeiB1vtgGK7Zb+11XakE2pPNNbSDkX1RvXrFS"
    "u1WvI52U5UBvo68TTM+HGrHZ7AAJk0po9IM84LSI+GpJE5gk/ORFGcHjOqTXrAl0+MXLOmZ4PpFyLU"
    "jGqTHtKq4Zu4tFT2DC3eraX28MlZ23r9pUHbwrP7TtEvnMA6AmdHU9p/9ICJuQRp/AWkGy69lg1pX/"
    "YI3JX1ydM9Bz9Schlt1nQHpJHo2DHSY0i/hL81irFBYJNiZDWHPoaadzPVHP7hh0On6JdDrQYSvNVC"
    "+le3mtd1YXrre7daqyfh5TDWej4UeYeeqRqMrC5lse53KP+yfP7HdSH05W8/fZurfJurfJurfJur/P"
    "q5iv+ifJuufJuufJuufJuuvMd05S+fvv7795+H9QzWH1bf43rUqlkPYjXrEazmdgRrPeeGtL3oHIk/"
    "oZCciOgSR2ibHlfo+UmQ7TENL2DeQTo9rEFVo2bYnkrbbxE88wcMaFB5Zh61ZgqTHkz3I1b9hT45xH"
    "tHp3oH7Tku1NyQbkZ8/npp+AmkMdt/BqlFuqMDo2fzAnY80RvAGdEQp/1Er+WzsDhdso/+5XeTpmY4"
    "PTBa5FcCqGd06IMDdP3W6bARgOkZbOzAZit4fmEK7ubbcbWXcACzX+gw5sAx/x7SmD2xbDhP+OXrL1"
    "rfGuH18vvfv15wADQwgpYwjaDnkexBSc9sNq5X9qSvIMN6Xl736jkm9Wz7tCIjHY3ISQ/s80LHC7P1"
    "bBQXCPV0WE+YVmeflI3VGYTsHemWMrA9c9KKpamezvc7KIineKakpCi7SyvCJaZ9NQr1NNl62lPt+W"
    "7Zx+p5/d6G9rRPXM9bNtbz+uF74n6Hel5epe/SWP4udRF6GaB7laI3/12K8T5GX4gt3m/1pO/STBY4"
    "9wXcyT77/TTdAembhfT1DDhaFD2h1Qd8FrPfZXwOgvQQ6hm+SzCvwWempDRb3R/vU6hnDd+ltT0r+C"
    "7d6hm+S888Ph37Lp2w7xtrzbLT8tQ72M1bPcN3KdNgptyeZ/p9UHGznW7P/Hcpft4jG9YIe9Im1JO+"
    "S89YzzHU09bSnraO9sx/l56rno6v487M64w0r+vuMp93fB133r6ft0fqbfbdXZ79u0T9/tzfpSXU09"
    "bSnraO9syv49KSdrLvWk97er70mHmyjeZLT17PCuy7rWS+ZCuZL9mnny8NoZ5VzJfsk3+XsN9v9Xzq"
    "7xKvp6uj36uYL9ln/y4N1J5P/V2aQj2f+7sE/d784dnXcSbUswb/0lrPCr5Lt3rm9+Oeaz5/q6etpT"
    "1tHXh/+nWcJwxsCBSMOxDCNuee98ie+4PSMXfgOevZLQl3AOka2ef7Z0OI7V3peUm4A09Zz2FIuANP"
    "3562jvb036pnrifjDjw13hl34DnrOXUJd+AJ+73Vpku4A8+Md84dePr2tHW051N/l3w9mS/82evp13"
    "ZPPT65L3xz3id7zuz9suHWjh3pdupS7sAz1rNbUu7As7enraM9w9ruievJuQPPPK/j3IFn/n5G3IEn"
    "7/dn/i6t2VvuwLO3p62jPZ/5u+TraSuZL9lK5ku2kvmSrWS+ZCuZL9lK5ku2lvmSrWS+ZCuZL9lK5k"
    "u2kvmSrWS+ZCuZL9k65kucO/DM3yXOHXjm71LEHXhiHEXcgScfn7aW9rR1tOfTf5eAO2B8sIGI3sTu"
    "Vdw8H5oNNbXFmkKUDPFlZ7LhH9nnnjTU9LWaNr0kNQ26vKlV7iENNXUH2hTpbXfNjthv8ZMQFYh4cM"
    "189Iw6KJeRZHeS9E7vF14eEJWvaQYyO9m8oSah9+Gg/VsRdbhND7XK3dqUH6fj7SE2Gs+GFoqzRUTF"
    "Z/nehKhNTd8FMuY0oiaGIp7uOaL8NKWDSF+550OzCVFyTYHAJ77sTPYOt7AkTYiqo00vSU2DLvlGe2"
    "v2oTbNSxOi9toUo0XdNZvzUjdPQlRYoDw9ogo1fTJEVdKml6Sm98DEOyFqt02fBlG2Ghtlq7FRthob"
    "ZauxUbYeG2WrsVG2Ghtlq7FRthobZWtBVHBPPj2iCjV9MkRV0qaXpKZPbKN22/QJEHV5jb3nh5x5wy"
    "M8vdea2l9VU328Kof8UqWaPqP3HEbotqaVeM+vNa3Eew6I2vVKPsHeCSDqY2p6ej8qg6hK2jTjPddl"
    "7/kbs98DUffwnt8JURWsohFRH1LT07O+bU1fq2nTSrzniChXy6yvBk8vIOpcTWHWh2l9GlHm1yDq6d"
    "qUj9YIUTV4zwFRNXjPAVG2Ghtlq7FRthobZauxUbYeG2WrsVG2FkRV0qaVeM/RRtlaEFWDpxcQ9TE1"
    "fQ9EVdKmlXjPAVFVeM9d5D1/5t53kff8mddRLvKeP/M6ykXe82e2US7ynj89ooJX8ukR9TE1/ZU2yk"
    "Xe8+du00tS0+e0US7ynj8/osKKDy4Czj0fmk2I+pCaGrgbOfc8VtPXatr0ktQUVN2kf332oTbNSxOi"
    "dttUqftn736lIu/5syPqXE1x1jexWd8JRMGs762Iero2xW9+iqhNTe+BiXdC1H6bPguibDU2ylZjo2"
    "w1NspWY6NsPTbKVmOjbC2IqqRNL0lNn9lG2VoQRV7JZ0fUx9T0PRBVSZtekpo+sY3ab9PHI8pdavGe"
    "u8uO9/xpTnO4y1N6z/NtWon33F1q8Z67Sy3ec3fZ8Z4/GaIqadNKvOfuUov33F1q8Z67y473PDMFeW"
    "v26VlfNENxl1q85+5Si/fcXWrxnrtLLd5zd9nxnj8Zoipp00q85+5Si/fcXWrxnrvLjvf8yRBlq7FR"
    "thobZeuxUbYaG2WrsVG2Ghtlq7FRthZE1eDpdZcd7/mTIaqSNq3Ee+4utXjPbTXcc/sG7vljmLLW1e"
    "I9t9Vwz2013HNbDffcVsM9t9Vwz2013HNbDffcVsM9t9Vwz2013HNbDffcVsM9t9Vwz+0buOePRFQl"
    "bVqJ99xWwz231XDPbTXcc1sN99xWwz231XDPbTXcc1sN99xWwz231XDPbTXcc1sN99xWwz231XDPbT"
    "Xcc1sP99xWwz231XDPbTXcc1sN99xWwz231XDPbTXcc1sN99xWwz231XDPbTXcc1sN99xWwz231XDP"
    "bTXcc1sN99xWwz231XDPbTXcc1sN99xWwz231XDPbTXcc1sN99xWwz231XDPbTXcc1sN99xWwz231X"
    "DPbTXcc1sN99xWwz231XDPbS3e8+YPtXjPrzWthHt+relrNW1aiff8WtNKvOeAqBo8vYCoGrjngKhK"
    "2rQS7zkgqgbvOSKqglU0IqoCXh8iqo42rcR7johytcz6avD0AqJq4J4Doipp00q854CoGrzngChbjY"
    "2y1dgoW42NstXYKFuPjbLV2ChbjY2y1dgoW42NsrUgqgZPLyCqBu45IKqSNq3Eew6Ienrv+ctvmy9f"
    "fzH6u89ff/nrp+9//PL551vtXzex0OESptzzodlQU1uqaWu60svOZMM/8s89aajpazVtetnWlHR5S6"
    "vcRRpq6vbbtB0nde/s23+kJyEqeCqn+faCsVchfab/uMR5aSgVniKiyKeavMakZeGnT8rGnrsPos61"
    "qR4f1KaXbZu2rCqZz1CU3eskewdRnf41iAo1HVhLsHR7TR+AzK+SjlG0UHpSHFF+FTivL8g+H5pNiJ"
    "Jraubiy85kwz+yzz1pQlQdbXpJahp0yTfaW7MPtWlemhC106Y3uN87Gz4p+SchKngqnx5RhZo+GaIq"
    "adNLUtN7YOKdELXXps+DKFuNjbLV2ChbjY2y1dgoW4+NstXYKFuNjbLV2ChbjY2ytSAqeCqfHlGFmj"
    "4Zoipp00tS0ye2UXtt+gyIurxG3vPVd4HPnOuonH1XT++1pvbX1NT0h6uCXqhBHfX1JTV9Ru+57nI1"
    "rcR7fq1pJd5zQFTBJ/00eyeAqA+p6WnveQ5RlbRp6j2XMPHrst8FUXtt+kSIqmAVjYj6kJqenvVta/"
    "paTZtW4j1HRLlaZn01eHoBUedqCrO+mc36TiAKZn1vRtTTtSnM+jKIqsF7DoiqwXsOiLLV2ChbjY2y"
    "1dgoW42NsvXYKFuNjbK1IKqSNq3Ee442ytaCqBo8vYCoj6npeyCqkjatxHsOiKrCex7fJPrMvR/fJP"
    "rM66j4JtFnXkfFN4k+s42KbxJ9ekQFr+TTI+pjavorbVR8k+hzt+klqelz2qj4JtHnRxSu+PTqeM8/"
    "H5pNiPqImuKxiezzWE1fq2nTy7amoOo2/auzD7WpIE2I2mnTpuvVvbNv/5GehCjvlXx+RJ2qKcz6NO"
    "NMnEEUzvreiqhna1OY9eUQFdf0Lph4J0TttenzIMpWY6NsNTbKVmOjbDU2ytZjo2w1NsrWgqhK2vSy"
    "relT2yhbC6K8V/L5EfUhNX0XRFXSppdtTZ/ZRu216TMgyl1q8Z67y473/GlOc7jLU3rP821aiffcXW"
    "rxnrtLLd5zd9nxnj8Zoipp00q85+5Si/fcXWrxnrtL2XuemYK8Nfv0rC+eobhLLd5zd6nFe+4utXjP"
    "3aUW77m7lL3nz4aoStq0Eu+5u9TiPXeXWrzn7lL2nj8bomw1NspWY6NsPTbKVmOjbDU2ylZjo2w1Ns"
    "rWgqgaPL3uUvaePxuiKmnTSrzn7lKL99xWwz23b+CeP4Ypa10t3nNbDffcVsM9t9Vwz2013HNbDffc"
    "VsM9t9Vwz2013HNbDffcVsM9t9Vwz2013HNbDffcnueePxRRlbRpJd5zWw333FbDPbfVcM9tNdxzWw"
    "333FbDPbfVcM9tNdxzWw333FbDPbfVcM9tNdxzWw333FbDPbfVcM9tPdxzWw333FbDPbfVcM9tNdxz"
    "Ww333FbDPbfVcM9tNdxzWw333FbDPbfVcM9tNdxzWw333FbDPbfVcM9tNdxzWw333FbDPbfVcM9tNd"
    "xzWw333FbDPbfVcM9tNdxzWw333FbDPbfVcM9tNdxzWw333FbDPbfVcM9tNdxzWw333FbDPbe1eM+b"
    "P9TiPb/WtBLu+bWmr9W0aSXe82tNK/GeA6Jq8PQComrgngOiKmnTSrzngKgavOeIqApW0YioCnh9iK"
    "g62rQS7zkiytUy66vB0wuIqoF7DoiqpE0r8Z4DomrwngOibDU2ylZjo2w1NspWY6NsPTbKVmOjbDU2"
    "ylZjo2w1NsrWgqgaPL2AqBq454CoStq0Eu85IOrpvecvv22+fP3FmFuVXy+/t9w7rc1659DNi9E0Pr"
    "1m6PUfA8uAJoN/NKtEO2FalsAyuZ8ZSjOK3nNNrxlrlfXaA/g28BDJEusVmRqqsWDal8k9xidrffur"
    "9dXNmt00mN7V0x3SM1vrogT8tLbMP//Lv0Iay3yN+hMq17Mm7JSQYVh/wk/YSktBAsuM+nNg3TZwdX"
    "rSE9PQtrKEYT8ZlPNlpv15uNa3/lx/atefMG329XSH9MzWuiTRdNTFmA5j6GJTPWHcdgyfDSuz5W0b"
    "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"
    "cucP4c8fmCQ5kU+e7z11/++un7H798/jmvE3AXesYxQ/8KOX5uJGTiIr/h+R7CN53+x9+///fPpEZH"
    "4YAb2I8P96gRnWzcxqNr+oW9X6stf7tnz4E9R76LD885eS70hEERnuug0KRGrmO8DdEnni9nBd5b+K"
    "bRl6+/fP7p66fvt0MM+STwBJ6kYk/NgoX2xANrzKi2wTcXtQnuFgUaNzyMJw/mCc8heY7sGSjPXBMZ"
    "N3johT9HxuraPF+kjCPP9xC+KfXf/ut3P3766dMPP9++qv/a/J/Nn//+u6/f/f3rj9dvxOe/NP+p+f"
    "Ofb48vv/v+b3/78bu/fP73nz5/bv7jf/rvn3/6wbr/4/b1/aff/eZmpP78999c///dX/4L/OvP13/9"
    "Ofzry/VfX67/am4Cv4E1IfvP9QP/m7Xz1n+1y+0/NwfmrYJ/+uOf/vh/3Sr3py9fv/ztWsH/9sc//c"
    "cfPv3jPzc//O5/fdXtV/O/WC1+uD5/WAu6/ePTP6jAcX33dTj8Bg+O3pT/8vOfb6/+9P33zdoKn689"
    "/HPz6afPzU3/n5u//VvzHz5//fzTv//PBjv9f/4HKO32yn/++4/ff/7HH79++WUt7tYmX67j47qoW/"
    "/3X/+fwM7xtKfgDEMfym3XV6vp9/P1EwLegNva+1qz//L5H59uLf7zy2/sq704+7q+FDmCt1vA4fcL"
    "FOYPuvbhd/x7pGPT37v1d4yLe/PoXAv7759/+cmX9mrt6+XShG/buufxm9frz/ArbBOsIRFuv+Lfwr"
    "b2um9w/fXV/y1MR2d4w8X/Ctu7txOqt1/9GzAQ4LD+6mz8XvwV/xYZod36Xhdqhm8Y4Vf/t6Dpjap1"
    "/fU1fi/Uwb36v8XAw/DeS/YNbqOFgV95O0BcrNuvLi5thtb+6Qtv6ybYf3P7zvzm2kP2taELURRI/c"
    "uf/vnl/weHW7qo0dMDAA=="
)


def _rna_pairtype(a, b):
    return _RNA_PT.get((a, b), 6)


def _rna_pairtable(db):
    """Dot-bracket -> partner array (pt[k]=partner or -1). Raises on
    imbalance or a stray glyph so a malformed structure fails loud."""
    pt = [-1] * len(db)
    st = []
    for k, c in enumerate(db):
        if c == '(':
            st.append(k)
        elif c == ')':
            if not st:
                raise ValueError(f"unbalanced ')' at {k}")
            o = st.pop()
            pt[o] = k
            pt[k] = o
        elif c != '.':
            raise ValueError(f"bad dot-bracket char {c!r} at {k}")
    if st:
        raise ValueError(f"unbalanced '(' at {st}")
    return pt


def _rna_reshape(flat, shape):
    n = 1
    for d in shape:
        n *= d
    if len(flat) != n:
        raise ValueError(f"reshape {shape}: expected {n}, got {len(flat)}")

    def build(fl, sh):
        if len(sh) == 1:
            return list(fl)
        step = 1
        for d in sh[1:]:
            step *= d
        return [build(fl[k * step:(k + 1) * step], sh[1:]) for k in range(sh[0])]
    return build(flat, shape)


def _rna_resolve_def(tbl, val):
    for i, x in enumerate(tbl):
        if x is None:
            tbl[i] = val
        elif isinstance(x, list):
            _rna_resolve_def(x, val)


def _rna_parse_params(text):
    """Parse a ViennaRNA params file (as text) into {section: flat list}
    plus {special-section: {seq: energy}} for the tri/tetra/hexaloops."""
    secs, special, cur = {}, {}, None
    SPECIAL = {'Hexaloops', 'Tetraloops', 'Triloops'}
    for raw in text.splitlines():
        s = raw.strip()
        if not s or s.startswith('##') or s.startswith('/*'):
            continue
        if s.startswith('#'):
            cur = s[1:].strip()
            special[cur] = {} if cur in SPECIAL else None
            if cur not in SPECIAL:
                secs[cur] = []
            continue
        if cur is None:
            continue
        if cur in SPECIAL:
            parts = s.split()
            if len(parts) >= 2 and parts[0].isalpha():
                special[cur][parts[0].replace('T', 'U').upper()] = int(parts[1])
            continue
        s = re.sub(r'/\*.*?\*/', ' ', s)
        for t in s.split():
            if t == 'INF':
                secs[cur].append(_RNA_INF)
            elif t == 'DEF':
                secs[cur].append(None)
            elif '.' in t:
                secs[cur].append(float(t))
            else:
                secs[cur].append(int(t))
    return secs, special


class _RNAModel:
    """Turner-2004 energy model + Zuker MFE folder. Construct once via
    `_rna_model()`; all methods are pure (thread-safe to share)."""

    def __init__(self, par_text):
        secs, special = _rna_parse_params(par_text)
        self.stack = _rna_reshape(secs['stack'], (7, 7))
        self.hairpin = secs['hairpin']
        self.bulge = secs['bulge']
        self.internal = secs['internal']
        self.mm_hairpin = _rna_reshape(secs['mismatch_hairpin'], (7, 5, 5))
        self.mm_internal = _rna_reshape(secs['mismatch_internal'], (7, 5, 5))
        self.mm_internal_1n = _rna_reshape(secs['mismatch_internal_1n'], (7, 5, 5))
        self.mm_internal_23 = _rna_reshape(secs['mismatch_internal_23'], (7, 5, 5))
        self.mm_multi = _rna_reshape(secs['mismatch_multi'], (7, 5, 5))
        self.mm_ext = _rna_reshape(secs['mismatch_exterior'], (7, 5, 5))
        self.dangle5 = _rna_reshape(secs['dangle5'], (7, 5))
        self.dangle3 = _rna_reshape(secs['dangle3'], (7, 5))
        self.int11 = _rna_reshape(secs['int11'], (7, 7, 5, 5))
        self.int21 = _rna_reshape(secs['int21'], (7, 7, 5, 5, 5))
        self.int22 = _rna_reshape(secs['int22'], (6, 6, 4, 4, 4, 4))
        misc = secs['Misc']
        self.terminalAU = misc[2]
        self.lxc = misc[4]
        ml = secs['ML_params']
        self.ml_base, self.ml_close, self.ml_branch = ml[0], ml[2], ml[4]
        self.tetra = special['Tetraloops']
        self.tri = special['Triloops']
        self.hexa = special['Hexaloops']
        # DEF -> -50 (a uniform default, measured against the reference)
        for tbl in (self.stack, self.mm_hairpin, self.mm_internal,
                    self.mm_internal_1n, self.mm_internal_23, self.mm_multi,
                    self.mm_ext, self.dangle5, self.dangle3,
                    self.int11, self.int21, self.int22):
            _rna_resolve_def(tbl, -50)

    # ---- loop energies (centi-kcal) ----
    def _loop_init(self, table, size):
        if size <= 30:
            return table[size]
        # ViennaRNA's C `(int)` cast TRUNCATES toward zero (no half-up
        # rounding); matching it keeps large-loop (>30 nt) energies exact
        # to the cent. A prior `+ 0.5` here diverged by +0.01 kcal on big
        # loops (the validation suite's max loop was 22 nt, so it never
        # caught it — found in the 2026-06-07 adversarial audit).
        return table[30] + int(self.lxc * _math.log(size / 30.0))

    def _stem_d2(self, t, s5i, s3i, mm):
        if s5i and s3i:
            e = mm[t][s5i][s3i]
        elif s5i:
            e = self.dangle5[t][s5i]
        elif s3i:
            e = self.dangle3[t][s3i]
        else:
            e = 0
        if t not in (0, 1):
            e += self.terminalAU
        return e

    def energy_stack(self, s, i, j):
        return self.stack[_rna_pairtype(s[i], s[j])][_rna_pairtype(s[j - 1], s[i + 1])]

    def energy_hairpin(self, s, i, j):
        size = j - i - 1
        if size < 3:
            return _RNA_INF
        sub = s[i:j + 1]
        if size == 3 and sub in self.tri:
            return self.tri[sub]
        if size == 4 and sub in self.tetra:
            return self.tetra[sub]
        if size == 6 and sub in self.hexa:
            return self.hexa[sub]
        t = _rna_pairtype(s[i], s[j])
        e = self._loop_init(self.hairpin, size)
        if size == 3:
            if t not in (0, 1):
                e += self.terminalAU
        else:
            e += self.mm_hairpin[t][_RNA_BI[s[i + 1]]][_RNA_BI[s[j - 1]]]
        return e

    def energy_bulge(self, s, i, j, a, b):
        size = (a - i - 1) + (j - b - 1)
        e = self._loop_init(self.bulge, size)
        if size == 1:
            e += self.stack[_rna_pairtype(s[i], s[j])][_rna_pairtype(s[b], s[a])]
        else:
            for (x, y) in ((s[i], s[j]), (s[a], s[b])):
                if _rna_pairtype(x, y) not in (0, 1):
                    e += self.terminalAU
        return e

    def energy_internal(self, s, i, j, a, b):
        n1 = a - i - 1
        n2 = j - b - 1
        t1 = _rna_pairtype(s[i], s[j])
        t2 = _rna_pairtype(s[b], s[a])
        v, is_special = None, True
        if n1 == 1 and n2 == 1:
            v = self.int11[t1][t2][_RNA_BI[s[i + 1]]][_RNA_BI[s[j - 1]]]
        elif n1 == 1 and n2 == 2:
            v = self.int21[t1][t2][_RNA_BI[s[i + 1]]][_RNA_BI[s[b + 1]]][_RNA_BI[s[j - 1]]]
        elif n1 == 2 and n2 == 1:
            v = self.int21[t2][t1][_RNA_BI[s[b + 1]]][_RNA_BI[s[i + 1]]][_RNA_BI[s[a - 1]]]
        elif n1 == 2 and n2 == 2:
            v = self.int22[t1][t2][_RNA_BI4[s[i + 1]]][_RNA_BI4[s[a - 1]]][
                _RNA_BI4[s[b + 1]]][_RNA_BI4[s[j - 1]]]
        else:
            is_special = False
        if is_special and v is not None:
            return v
        size = n1 + n2
        e = self._loop_init(self.internal, size)
        e += min(300, abs(n1 - n2) * 60)        # NINIO asymmetry (m=60, max=300)
        if n1 == 1 or n2 == 1:
            mm = self.mm_internal_1n
        elif (n1, n2) in ((2, 3), (3, 2)):
            mm = self.mm_internal_23
        else:
            mm = self.mm_internal
        e += mm[t1][_RNA_BI[s[i + 1]]][_RNA_BI[s[j - 1]]]
        e += mm[t2][_RNA_BI[s[b + 1]]][_RNA_BI[s[a - 1]]]
        return e

    def energy_multiloop(self, s, i, j, inner):
        e = self.ml_close + self.ml_branch * (len(inner) + 1)
        e += self.ml_base * ((j - i - 1) - sum(b - a + 1 for (a, b) in inner))
        for (a, b) in inner:
            t = _rna_pairtype(s[a], s[b])
            e += self._stem_d2(t, _RNA_BI[s[a - 1]], _RNA_BI[s[b + 1]], self.mm_multi)
        tc = _rna_pairtype(s[j], s[i])
        e += self._stem_d2(tc, _RNA_BI[s[j - 1]], _RNA_BI[s[i + 1]], self.mm_multi)
        return e

    def _inner_pairs(self, pt, i, j):
        res, k = [], i + 1
        while k < j:
            if pt[k] > k:
                res.append((k, pt[k]))
                k = pt[k] + 1
            else:
                k += 1
        return res

    def _energy_enclosed(self, s, pt, i, j):
        inner = self._inner_pairs(pt, i, j)
        if not inner:
            return self.energy_hairpin(s, i, j)
        if len(inner) == 1:
            (a, b) = inner[0]
            lg, rg = a - i - 1, j - b - 1
            if lg == 0 and rg == 0:
                e = self.energy_stack(s, i, j)
            elif lg == 0 or rg == 0:
                e = self.energy_bulge(s, i, j, a, b)
            else:
                e = self.energy_internal(s, i, j, a, b)
            return e + self._energy_enclosed(s, pt, a, b)
        e = self.energy_multiloop(s, i, j, inner)
        for (a, b) in inner:
            e += self._energy_enclosed(s, pt, a, b)
        return e

    def eval_structure(self, s, db):
        pt = _rna_pairtable(db)
        n = len(s)
        # Fail loud on a structurally-valid bracketing that's chemically
        # impossible: a non-canonical pair (the energy tables have no
        # type-6 entry — would IndexError) or a hairpin loop < 3 nt (would
        # leak the _RNA_INF sentinel as a finite ~1e7 energy). The MFE
        # folder never emits either; this guards the eval-only path.
        for a in range(n):
            b = pt[a]
            if b > a:
                if _rna_pairtype(s[a], s[b]) == 6:
                    raise ValueError(
                        f"structure has a non-canonical pair at {a},{b} "
                        f"({s[a]}·{s[b]})")
                if (b - a - 1 < 3
                        and not any(pt[x] > x for x in range(a + 1, b))):
                    raise ValueError(
                        f"structure has an infeasible hairpin loop "
                        f"(<3 nt) closed at {a},{b}")
        total, k = 0, 0
        while k < n:
            if pt[k] > k:
                i, j = k, pt[k]
                t = _rna_pairtype(s[i], s[j])
                s5 = _RNA_BI[s[i - 1]] if i - 1 >= 0 else 0
                s3 = _RNA_BI[s[j + 1]] if j + 1 < n else 0
                total += self._stem_d2(t, s5, s3, self.mm_ext)
                total += self._energy_enclosed(s, pt, i, j)
                k = j + 1
            else:
                k += 1
        if total >= _RNA_INF:
            raise ValueError("structure contains an infeasible loop")
        return total / 100.0

    # ---- d2 helix-end contributions for the folder ----
    def _d2_ext(self, s, i, j, n):
        t = _rna_pairtype(s[i], s[j])
        s5 = _RNA_BI[s[i - 1]] if i > 0 else 0
        s3 = _RNA_BI[s[j + 1]] if j + 1 < n else 0
        return self._stem_d2(t, s5, s3, self.mm_ext)

    def _d2_ml(self, s, i, j, n):
        t = _rna_pairtype(s[i], s[j])
        s5 = _RNA_BI[s[i - 1]] if i > 0 else 0
        s3 = _RNA_BI[s[j + 1]] if j + 1 < n else 0
        return self._stem_d2(t, s5, s3, self.mm_multi)

    def _d2_ml_close(self, s, i, j):
        tc = _rna_pairtype(s[j], s[i])
        return self._stem_d2(tc, _RNA_BI[s[j - 1]], _RNA_BI[s[i + 1]], self.mm_multi)

    def _loop_e(self, s, i, j, p, q):
        lg, rg = p - i - 1, j - q - 1
        if lg == 0 and rg == 0:
            return self.energy_stack(s, i, j)
        if lg == 0 or rg == 0:
            return self.energy_bulge(s, i, j, p, q)
        return self.energy_internal(s, i, j, p, q)

    def fold(self, s):
        """Minimum-free-energy fold -> (dot_bracket, dg_kcal)."""
        n = len(s)
        if n < 5:
            return '.' * n, 0.0
        INF, ML = _RNA_INF, _RNA_MAXLOOP
        a, b = self.ml_close, self.ml_branch
        V = [[INF] * n for _ in range(n)]
        M = [[INF] * n for _ in range(n)]
        M1 = [[INF] * n for _ in range(n)]
        for d in range(3, n):
            for i in range(0, n - d):
                j = i + d
                if _rna_pairtype(s[i], s[j]) != 6:
                    best = self.energy_hairpin(s, i, j)
                    pmax = min(j - 1, i + ML + 1)
                    for p in range(i + 1, pmax + 1):
                        lg = p - i - 1
                        qmin = max(p + 1, j - 1 - (ML - lg))
                        for q in range(qmin, j):
                            if V[p][q] >= INF or _rna_pairtype(s[p], s[q]) == 6:
                                continue
                            cand = self._loop_e(s, i, j, p, q) + V[p][q]
                            if cand < best:
                                best = cand
                    base = a + b + self._d2_ml_close(s, i, j)
                    for u in range(i + 2, j - 1):
                        if M[i + 1][u] < INF and M1[u + 1][j - 1] < INF:
                            cand = base + M[i + 1][u] + M1[u + 1][j - 1]
                            if cand < best:
                                best = cand
                    V[i][j] = best
                m1 = M1[i][j - 1] if M1[i][j - 1] < INF else INF
                if V[i][j] < INF:
                    cand = V[i][j] + b + self._d2_ml(s, i, j, n)
                    if cand < m1:
                        m1 = cand
                M1[i][j] = m1
                m = M[i][j - 1] if M[i][j - 1] < INF else INF
                for k in range(i, j):
                    if V[k][j] < INF:
                        cand = V[k][j] + b + self._d2_ml(s, k, j, n)
                        if cand < m:
                            m = cand
                for u in range(i, j):
                    if M[i][u] < INF and M1[u + 1][j] < INF:
                        cand = M[i][u] + M1[u + 1][j]
                        if cand < m:
                            m = cand
                M[i][j] = m
        F = [0] * n
        for j in range(1, n):
            best = F[j - 1]
            for i in range(0, j):
                if V[i][j] < INF:
                    prev = F[i - 1] if i > 0 else 0
                    cand = prev + V[i][j] + self._d2_ext(s, i, j, n)
                    if cand < best:
                        best = cand
            F[j] = best
        pairs = []
        self._tb_ext(s, n, V, M, M1, F, n - 1, pairs)
        db = ['.'] * n
        for (x, y) in pairs:
            db[x], db[y] = '(', ')'
        return ''.join(db), F[n - 1] / 100.0

    def _tb_ext(self, s, n, V, M, M1, F, j, pairs):
        while j > 0:
            if F[j] == F[j - 1]:
                j -= 1
                continue
            for i in range(0, j):
                if V[i][j] >= _RNA_INF:
                    continue
                prev = F[i - 1] if i > 0 else 0
                if F[j] == prev + V[i][j] + self._d2_ext(s, i, j, n):
                    pairs.append((i, j))
                    self._tb_V(s, n, V, M, M1, i, j, pairs)
                    j = i - 1
                    break
            else:
                break

    def _tb_V(self, s, n, V, M, M1, i, j, pairs):
        if V[i][j] == self.energy_hairpin(s, i, j):
            return
        pmax = min(j - 1, i + _RNA_MAXLOOP + 1)
        for p in range(i + 1, pmax + 1):
            lg = p - i - 1
            qmin = max(p + 1, j - 1 - (_RNA_MAXLOOP - lg))
            for q in range(qmin, j):
                if V[p][q] >= _RNA_INF or _rna_pairtype(s[p], s[q]) == 6:
                    continue
                if V[i][j] == self._loop_e(s, i, j, p, q) + V[p][q]:
                    pairs.append((p, q))
                    self._tb_V(s, n, V, M, M1, p, q, pairs)
                    return
        base = self.ml_close + self.ml_branch + self._d2_ml_close(s, i, j)
        for u in range(i + 2, j - 1):
            if M[i + 1][u] < _RNA_INF and M1[u + 1][j - 1] < _RNA_INF and \
               V[i][j] == base + M[i + 1][u] + M1[u + 1][j - 1]:
                self._tb_M(s, n, V, M, M1, i + 1, u, pairs)
                self._tb_M1(s, n, V, M, M1, u + 1, j - 1, pairs)
                return
        raise RuntimeError(f"V traceback failed at {i},{j}")

    def _tb_M(self, s, n, V, M, M1, i, j, pairs):
        if j > i and M[i][j] == M[i][j - 1]:
            self._tb_M(s, n, V, M, M1, i, j - 1, pairs)
            return
        for k in range(i, j):
            if V[k][j] < _RNA_INF and \
               M[i][j] == V[k][j] + self.ml_branch + self._d2_ml(s, k, j, n):
                pairs.append((k, j))
                self._tb_V(s, n, V, M, M1, k, j, pairs)
                return
        for u in range(i, j):
            if M[i][u] < _RNA_INF and M1[u + 1][j] < _RNA_INF and \
               M[i][j] == M[i][u] + M1[u + 1][j]:
                self._tb_M(s, n, V, M, M1, i, u, pairs)
                self._tb_M1(s, n, V, M, M1, u + 1, j, pairs)
                return
        raise RuntimeError(f"M traceback failed at {i},{j}")

    def _tb_M1(self, s, n, V, M, M1, i, j, pairs):
        if j > i and M1[i][j] == M1[i][j - 1]:
            self._tb_M1(s, n, V, M, M1, i, j - 1, pairs)
            return
        if V[i][j] < _RNA_INF and \
           M1[i][j] == V[i][j] + self.ml_branch + self._d2_ml(s, i, j, n):
            pairs.append((i, j))
            self._tb_V(s, n, V, M, M1, i, j, pairs)
            return
        raise RuntimeError(f"M1 traceback failed at {i},{j}")

    # ---- bound-state heterodimer (cofold) ----
    def _junction(self, s, i, j, cut):
        """Energy of a cut-spanning pair's duplex junction — replaces the
        hairpin for inter-strand pairs. The two cut-facing bases dangle
        across the backbone break, scored via the exterior mismatch table
        in the inward (reversed-pair) orientation."""
        s5 = _RNA_BI[s[j - 1]] if j - 1 >= cut else 0
        s3 = _RNA_BI[s[i + 1]] if i + 1 < cut else 0
        return self._stem_d2(_rna_pairtype(s[j], s[i]), s5, s3, self.mm_ext)

    def cofold(self, a_seq, b_seq):
        """Bound-state heterodimer free energy of strands A & B (kcal/mol).
        Concatenates A+B with an inter-strand cut and computes the BOUND
        complex (DuplexInit always paid — the ribosome-bound state).
        Matches ViennaRNA RNAcofold on binding duplexes; intra-strand
        structure inside the duplex (the footprint) is forbidden, which is
        the constrained bound state the translation-initiation model
        uses. Energy only (no traceback)."""
        s = a_seq + b_seq
        cut = len(a_seq)
        n = len(s)
        if n < 2 or cut == 0 or cut == n:
            return 0.0
        INF, MLP = _RNA_INF, _RNA_MAXLOOP
        a, b = self.ml_close, self.ml_branch
        V = [[INF] * n for _ in range(n)]
        M = [[INF] * n for _ in range(n)]
        M1 = [[INF] * n for _ in range(n)]
        for d in range(1, n):
            for i in range(0, n - d):
                j = i + d
                if _rna_pairtype(s[i], s[j]) != 6:
                    if i < cut <= j:
                        best = self._junction(s, i, j, cut)
                    elif j - i - 1 >= 3:
                        best = self.energy_hairpin(s, i, j)
                    else:
                        best = INF
                    pmax = min(j - 1, i + MLP + 1)
                    for p in range(i + 1, pmax + 1):
                        lg = p - i - 1
                        qmin = max(p + 1, j - 1 - (MLP - lg))
                        for q in range(qmin, j):
                            if V[p][q] >= INF or _rna_pairtype(s[p], s[q]) == 6:
                                continue
                            if (i < cut <= p) or (q < cut <= j):
                                continue        # loop's unpaired stretch spans the cut
                            c = self._loop_e(s, i, j, p, q) + V[p][q]
                            if c < best:
                                best = c
                    base = a + b + self._d2_ml_close(s, i, j)
                    for u in range(i + 2, j - 1):
                        if M[i + 1][u] < INF and M1[u + 1][j - 1] < INF:
                            c = base + M[i + 1][u] + M1[u + 1][j - 1]
                            if c < best:
                                best = c
                    V[i][j] = best
                m1 = M1[i][j - 1] if M1[i][j - 1] < INF else INF
                if V[i][j] < INF:
                    c = V[i][j] + b + self._d2_ml(s, i, j, n)
                    if c < m1:
                        m1 = c
                M1[i][j] = m1
                m = M[i][j - 1] if M[i][j - 1] < INF else INF
                for k in range(i, j):
                    if V[k][j] < INF:
                        c = V[k][j] + b + self._d2_ml(s, k, j, n)
                        if c < m:
                            m = c
                for u in range(i, j):
                    if M[i][u] < INF and M1[u + 1][j] < INF:
                        c = M[i][u] + M1[u + 1][j]
                        if c < m:
                            m = c
                M[i][j] = m
        F = [0] * n
        for j in range(1, n):
            best = F[j - 1]
            for i in range(0, j):
                if V[i][j] < INF:
                    prev = F[i - 1] if i > 0 else 0
                    c = prev + V[i][j] + self._d2_ext(s, i, j, n)
                    if c < best:
                        best = c
            F[j] = best
        # + DuplexInit (Turner-2004 Misc, +4.10 kcal) for the bound state.
        return (F[n - 1] + 410) / 100.0


_RNA_MODEL_SINGLETON = None
_RNA_MODEL_LOCK = _threading.Lock()


def _rna_model():
    # Double-checked locking: the lock-free fast path serves the common
    # case (singleton already built); the lock serialises the one-time
    # build so concurrent first-requests (the threaded agent server) can't
    # each parse a model. Construction is pure + idempotent, so the lock is
    # belt-and-braces, but it makes the "thread-safe singleton" real.
    global _RNA_MODEL_SINGLETON
    if _RNA_MODEL_SINGLETON is None:
        with _RNA_MODEL_LOCK:
            if _RNA_MODEL_SINGLETON is None:
                text = _gzip.decompress(
                    _base64.b64decode(_RNA_TURNER_PARAMS_GZ_B64)).decode('ascii')
                _RNA_MODEL_SINGLETON = _RNAModel(text)
    return _RNA_MODEL_SINGLETON


def _rna_normalize(seq):
    if not isinstance(seq, str):
        raise ValueError("sequence must be a string")
    s = seq.strip().upper().replace('T', 'U')
    if not s:
        raise ValueError("empty sequence")
    bad = set(s) - {'A', 'C', 'G', 'U'}
    if bad:
        raise ValueError(f"RNA folding needs unambiguous A/C/G/U; got {sorted(bad)}")
    return s


def _rna_fold(seq, *, max_len=_RNA_FOLD_MAX_LEN):
    """Fold an RNA/DNA sequence to its minimum-free-energy secondary
    structure. Returns (dot_bracket, dg_kcal_per_mol). DNA T is read as
    U. Raises ValueError on empty / ambiguous / over-length input."""
    s = _rna_normalize(seq)
    if len(s) > max_len:
        raise ValueError(f"sequence too long to fold ({len(s)} > {max_len} nt cap)")
    return _rna_model().fold(s)


def _rna_mfe(seq, *, max_len=_RNA_FOLD_MAX_LEN):
    """Minimum free energy (kcal/mol) only — the structure is discarded."""
    return _rna_fold(seq, max_len=max_len)[1]


def _rna_eval_structure(seq, dot_bracket):
    """Free energy (kcal/mol) of a GIVEN secondary structure on `seq`."""
    s = _rna_normalize(seq)
    if len(s) != len(dot_bracket):
        raise ValueError("sequence / structure length mismatch")
    return _rna_model().eval_structure(s, dot_bracket)


_RNA_COFOLD_MAX_LEN = 400               # combined A+B length cap (O(n^3) DP)


def _rna_cofold(seq_a, seq_b, *, max_len=_RNA_COFOLD_MAX_LEN):
    """Bound-state heterodimer free energy (kcal/mol) of two strands — the
    ΔG of strand B bound to strand A (e.g. the 16S anti-SD tail hybridized
    to an mRNA window). DNA `T` is read as `U`. The bound state is forced
    (DuplexInit always paid), matching ViennaRNA RNAcofold on binding
    duplexes; a weak / non-complementary pair returns a high (unfavorable)
    ΔG rather than reporting 'unbound'. Raises ValueError on empty /
    ambiguous / over-length input."""
    a = _rna_normalize(seq_a)
    b = _rna_normalize(seq_b)
    if len(a) + len(b) > max_len:
        raise ValueError(
            f"combined length too long to cofold "
            f"({len(a) + len(b)} > {max_len} nt cap)")
    return _rna_model().cofold(a, b)


# ── Ribosome binding site strength (E. coli translation initiation) ─────────
#
# A biophysically-grounded RELATIVE estimate of translation-initiation
# strength, built on the validated RNA folder + cofold. The STRUCTURAL
# energies are exact (`_rna_fold` / `_rna_cofold`, validated to the cent vs
# ViennaRNA); the constants below — the Boltzmann factor β, the
# spacing-penalty curve, and the start-codon ΔG — are literature-standard
# empirical CALIBRATION values, not first-principles. So only RATIOS
# between RBSs are meaningful: this is a tuning / ranking score, NOT an
# absolute expression rate. Validated by relative ranking on the canonical
# determinants (SD strength, 5'UTR occlusion, spacing, start codon), not
# against an absolute thermodynamic oracle.
#
#   ΔG_total = ΔG_hybrid(best SD register) + ΔG_start + ΔG_spacing − ΔG_mRNA
#   strength ∝ exp(−β · ΔG_total)

_RBS_ANTI_SD = 'ACCUCCUUA'         # E. coli 16S rRNA 3' tail (anti-SD), 5'->3'
_RBS_BETA = 0.45                   # mol/kcal — apparent Boltzmann factor (calibration)
_RBS_OPT_SPACING = 5               # optimal SD-to-start aligned spacing (nt)
_RBS_WINDOW = 35                   # nt up/downstream of the start folded for ΔG_mRNA
_RBS_SPACING_SCAN = range(3, 13)   # aligned-spacing registers scanned for the SD
# start-codon : initiator-tRNA(fMet) hybridisation ΔG (kcal/mol, favourable);
# a non-canonical start gets 0 (no favourable initiation). Calibration.
_RBS_START_DG = {'AUG': -1.19, 'GUG': -0.075, 'UUG': -0.075, 'CUG': -0.03,
                 'AUU': -0.03, 'AUC': -0.03, 'AUA': -0.03}


def _rbs_spacing_penalty(d):
    """ΔG penalty (kcal/mol) for an SD-to-start spacing of `d` nt deviating
    from the ~5-nt optimum. Asymmetric: too-short (steric clash with the
    ribosome) is penalised far harder than too-long (entropic). Calibration."""
    if d == _RBS_OPT_SPACING:
        return 0.0
    if d < _RBS_OPT_SPACING:
        return 0.20 * (_RBS_OPT_SPACING - d) ** 2
    return 0.05 * (d - _RBS_OPT_SPACING) ** 2


_RBS_STRENGTH_MAX_LEN = 100000          # mRNA length cap. The fold window is
#                                          bounded (≤2·_RBS_WINDOW) so compute is
#                                          ~constant in length, but cap the O(n)
#                                          normalize/scan so a pathological input
#                                          can't balloon memory. Matches the cap
#                                          discipline of _rna_fold / _rna_cofold /
#                                          _rbs_design.


def _rbs_strength(mrna, start_pos):
    """Relative E. coli translation-initiation strength of the ribosome
    binding site preceding the start codon at `start_pos` (0-based) in
    `mrna` (RNA or DNA; T read as U). Returns a dict::

        {dg_total, dg_mrna, dg_hybrid, spacing, rel_strength}

    `rel_strength` ∝ exp(−β·dg_total): only RATIOS between RBSs are
    meaningful (a ranking score, not an absolute rate). Captures
    SD:anti-SD complementarity, the 5'UTR structure that occludes the site
    (incl. the upstream standby region, via the folded window), the
    SD-to-start spacing, and the start codon. Raises ValueError on bad
    input; returns rel_strength 0.0 when the start is too close to the 5'
    end for an SD to fit."""
    s = mrna.strip().upper().replace('T', 'U') if isinstance(mrna, str) else None
    if not s:
        raise ValueError("mRNA must be a non-empty string")
    bad = set(s) - {'A', 'C', 'G', 'U'}
    if bad:
        raise ValueError(f"mRNA needs unambiguous A/C/G/U; got {sorted(bad)}")
    if len(s) > _RBS_STRENGTH_MAX_LEN:
        raise ValueError(f"mRNA too long ({len(s)} > {_RBS_STRENGTH_MAX_LEN} nt cap)")
    if (not isinstance(start_pos, int) or isinstance(start_pos, bool)
            or not (0 <= start_pos <= len(s) - 3)):
        raise ValueError(f"start_pos {start_pos!r} out of range for length {len(s)}")
    dg_start = _RBS_START_DG.get(s[start_pos:start_pos + 3], 0.0)
    w0 = max(0, start_pos - _RBS_WINDOW)
    w1 = min(len(s), start_pos + _RBS_WINDOW)
    dg_mrna = _rna_mfe(s[w0:w1])
    best = None
    best_d, best_hybrid = 0, 0.0
    for d in _RBS_SPACING_SCAN:
        end = start_pos - d
        begin = end - len(_RBS_ANTI_SD)
        if begin < 0:
            continue
        dg_h = _rna_cofold(s[begin:end], _RBS_ANTI_SD)
        dg_f = dg_h + dg_start + _rbs_spacing_penalty(d)
        if best is None or dg_f < best:
            best, best_d, best_hybrid = dg_f, d, dg_h
    if best is None:                       # start too close to the 5' end
        # dg_total is None (not float('inf')) so the dict is JSON-valid —
        # `Infinity` is not legal JSON and broke the agent endpoint.
        return {'dg_total': None, 'dg_mrna': round(dg_mrna, 2),
                'dg_hybrid': None, 'spacing': None, 'rel_strength': 0.0}
    dg_total = best - dg_mrna
    return {'dg_total': round(dg_total, 2), 'dg_mrna': round(dg_mrna, 2),
            'dg_hybrid': round(best_hybrid, 2), 'spacing': best_d,
            'rel_strength': round(_math.exp(-_RBS_BETA * dg_total), 3)}


# Graded Shine-Dalgarno library (complementarity to the anti-SD, strong → none)
# + spacer lengths, for reverse RBS design. The forward model ranks them;
# this just spans the strength range so a target can be matched.
_RBS_DESIGN_SD_LADDER = ['UAAGGAGGU', 'AAGGAGGU', 'AGGAGGA', 'AGGAGG', 'UAAGGAG',
                         'GGAGGU', 'AGGAG', 'GGAGG', 'AGGA', 'GAGGA', 'AGAGA',
                         'AAGAA', 'ACAUA', '']
_RBS_DESIGN_SPACERS = {4: 'AAUA', 5: 'AAUAA', 6: 'AACAAU', 7: 'AACAAUA',
                       8: 'AACAAUAA', 9: 'AACAAUAAU'}
_RBS_DESIGN_UPSTREAM = 'UUAAUUAAUU'      # low-structure 5' context (standby region)
_RBS_DESIGN_MAX_CDS = 50000              # cds length cap (the design does
#                                          O(len) passes — bounds worst-case cost)
_RBS_DESIGN_COARSE_SPACER = 6            # spacer length used to rank SDs in the
#                                          coarse pass of the coarse-to-fine search


def _rbs_design(cds, target_strength, *, upstream=_RBS_DESIGN_UPSTREAM):
    """Design a 5'UTR (Shine-Dalgarno + spacer) preceding `cds` (which must
    begin with the start codon) to achieve a target RELATIVE RBS strength.
    Searches a graded SD × spacer library, scores each construct with
    `_rbs_strength`, and returns the design closest to `target_strength`::

        {utr, full, sd, spacing, rel_strength, dg_total,
         achievable_min, achievable_max, on_target}

    Strength is relative (see `_rbs_strength`). A target outside the
    CDS-achievable range yields the nearest achievable design and
    `on_target=False`. Raises ValueError on bad input."""
    c = cds.strip().upper().replace('T', 'U') if isinstance(cds, str) else ''
    if not c or (set(c) - {'A', 'C', 'G', 'U'}):
        raise ValueError("cds must be a non-empty A/C/G/U(T) string")
    if len(c) < 3:
        raise ValueError("cds must include the start codon (>= 3 nt)")
    if len(c) > _RBS_DESIGN_MAX_CDS:
        raise ValueError(
            f"cds too long ({len(c)} > {_RBS_DESIGN_MAX_CDS} nt cap)")
    if (isinstance(target_strength, bool)
            or not isinstance(target_strength, (int, float))
            or not _math.isfinite(target_strength)
            or target_strength < 0):
        raise ValueError(
            "target_strength must be a finite non-negative number")
    up = (upstream or '').strip().upper().replace('T', 'U')
    if set(up) - {'A', 'C', 'G', 'U'}:
        raise ValueError("upstream must be A/C/G/U(T)")
    # Coarse-to-fine search. A full 14 SD × 6 spacer = 84-eval scan costs
    # ~5-7 s. Instead: (coarse) rank every SD at one representative spacer —
    # 14 evals — then (fine) scan every spacer for the SD that brackets the
    # target plus its two ladder neighbours — ~15 evals. The fine pass
    # covers the target's neighbourhood densely, so the achievable extreme
    # NEAR the target (the one that drives `on_target`) is exact; the FAR
    # extreme is an estimate. ~29 evals total (~2.5-3x faster).
    best = None
    lo, hi = float('inf'), -1.0
    spacers = _RBS_DESIGN_SPACERS
    coarse_len = _RBS_DESIGN_COARSE_SPACER

    def _consider(sd, slen, sp):
        nonlocal best, lo, hi
        utr = up + sd + sp
        r = _rbs_strength(utr + c, len(utr))
        v = r['rel_strength']
        lo, hi = min(lo, v), max(hi, v)
        if best is None or abs(v - target_strength) < abs(
                best['rel_strength'] - target_strength):
            best = {'utr': utr, 'sd': sd, 'spacing': slen,
                    'rel_strength': v, 'dg_total': r['dg_total']}
        return v

    coarse = []
    for i, sd in enumerate(_RBS_DESIGN_SD_LADDER):
        coarse.append((_consider(sd, coarse_len, spacers[coarse_len]), i))
    coarse.sort(key=lambda x: abs(x[0] - target_strength))
    best_i = coarse[0][1]
    last = len(_RBS_DESIGN_SD_LADDER) - 1
    for i in {best_i, max(0, best_i - 1), min(last, best_i + 1)}:
        sd = _RBS_DESIGN_SD_LADDER[i]
        for slen, sp in spacers.items():
            if slen == coarse_len:
                continue                     # already scored in the coarse pass
            _consider(sd, slen, sp)
    assert best is not None              # the SD ladder is non-empty -> always set
    best['full'] = best['utr'] + c
    best['achievable_min'] = round(lo, 3)
    best['achievable_max'] = round(hi, 3)
    best['on_target'] = lo <= target_strength <= hi
    return best


def _assemble_operon(genes, *, promoter='', terminator='',
                     leader=_RBS_DESIGN_UPSTREAM):
    """Assemble a contiguous bacterial operon — promoter + (RBS + CDS) per
    gene + terminator — CONTEXT-AWARE: each RBS is reverse-designed against
    the REAL upstream sequence (the promoter, or the preceding gene's 3'
    end), so the achieved in-context strength tracks the target. (Designing
    each RBS in isolation then concatenating does NOT — the upstream can
    occlude it.) When a gene's target is unreachable in its context (e.g.
    the previous CDS's 3' end sequesters the SD), the nearest achievable is
    used and that gene's `on_target` is False — a real, useful signal.

    `genes`: a non-empty list of dicts {cds, target_strength, name?} — each
    `cds` begins with the start codon (DNA `T` read as `U`). `leader` is
    the low-structure 5' standby used when there is no promoter. Returns::

        {sequence, layout, genes}

    `sequence` is DNA (T). `layout` is the ordered element map
    [{kind, name, start, end}] (kind ∈ promoter/rbs/cds/terminator) where
    `sequence[start:end]` is EXACTLY that element — contiguous, no gaps or
    overlaps. `genes` is the per-gene report [{name, target, cds_len, rbs,
    spacing, rel_strength, on_target}]. Raises ValueError on bad input."""
    if not isinstance(genes, (list, tuple)) or not genes:
        raise ValueError("genes must be a non-empty list")

    def _norm(x, label):
        x = (x or '').strip().upper().replace('T', 'U')
        if set(x) - {'A', 'C', 'G', 'U'}:
            raise ValueError(f"{label} must be A/C/G/U(T)")
        return x
    prom = _norm(promoter, 'promoter')
    term = _norm(terminator, 'terminator')
    lead = _norm(leader, 'leader')

    assembled = prom
    layout, anchors = [], []
    if prom:
        layout.append({'kind': 'promoter', 'name': 'promoter',
                       'start': 0, 'end': len(prom)})
    for i, g in enumerate(genes):
        if not isinstance(g, dict):
            raise ValueError(f"gene {i} must be a dict")
        cds = g.get('cds')
        if not isinstance(cds, str):
            raise ValueError(f"gene {i}: missing string 'cds'")
        cds = cds.strip().upper().replace('T', 'U')
        target = g.get('target_strength', g.get('target'))
        name = str(g.get('name') or f"gene{i + 1}")
        if assembled:
            ctx = assembled[-_RBS_WINDOW:]
            d = _rbs_design(cds, target, upstream=ctx)    # validates cds + target
            rbs = d['utr'][len(ctx):]                      # strip the already-present upstream
        else:
            d = _rbs_design(cds, target, upstream=lead)
            rbs = d['utr']                                 # the leader is the operon's 5' start
        rbs_start = len(assembled)
        layout.append({'kind': 'rbs', 'name': f"{name} RBS",
                       'start': rbs_start, 'end': rbs_start + len(rbs)})
        assembled += rbs
        cds_start = len(assembled)
        layout.append({'kind': 'cds', 'name': name,
                       'start': cds_start, 'end': cds_start + len(cds)})
        assembled += cds
        anchors.append({'name': name, 'target': target, 'cds_len': len(cds),
                        'rbs': rbs.replace('U', 'T'), 'spacing': d['spacing'],
                        'cds_start': cds_start})
    if term:
        layout.append({'kind': 'terminator', 'name': 'terminator',
                       'start': len(assembled), 'end': len(assembled) + len(term)})
        assembled += term

    report = []
    for a in anchors:
        try:
            rel = _rbs_strength(assembled, a['cds_start'])['rel_strength']
        except ValueError:
            rel = None
        tgt = a['target']
        on = (rel is not None and isinstance(tgt, (int, float))
              and abs(rel - tgt) <= 0.25 * max(tgt, 1e-9))
        report.append({'name': a['name'], 'target': tgt, 'cds_len': a['cds_len'],
                       'rbs': a['rbs'], 'spacing': a['spacing'],
                       'rel_strength': rel, 'on_target': on})
    return {'sequence': assembled.replace('U', 'T'), 'layout': layout,
            'genes': report}




# ── Restriction-site scanner + enzyme digest (Phase D, extracted from the
# hub). Sacred invariants #1/#2/#6 live here. Pure given _state (caches +
# `_scan_catalog_hook`/`_all_enzymes_hook` getters) + `_rc`/`_iupac_pattern`
# above. `_rebuild_scan_catalog` (writes the catalog) + `_all_enzymes` (reads
# dataaccess) stay hub-side and feed in through the _state getters.
# ──────────────────────────────────────────────────────────────────────────

def _iter_match_starts(pat, s: str):
    """Yield the start index of EVERY match of `pat` in `s`, INCLUDING
    self-overlapping tandem matches.

    `re.finditer` is non-overlapping: it resumes each search past the end of
    the previous match, so it silently drops the second of two overlapping
    recognition sites — a site whose string has a border (a proper prefix that
    equals a proper suffix), e.g. NotI ``GCGGCCGC`` (border "GC"), BstUI
    ``CGCG``, HhaI ``GCGC``. That undercounts cutters: a real 2-cutter is
    badged ``unique`` (the user then linearises on a site that actually cuts
    twice) and `_enzyme_cuts` hands too few fragments into the cloning
    simulation. Advancing the search ONE base past each match's start recovers
    every occurrence — the same overlap-preserving intent as the ``(?=...)``
    lookahead already used by `_forbidden_hit_set` / `_search_subsequence`.
    Yields in ascending start order (like `finditer`), so downstream `seen`
    dedup and result ordering are unchanged bar the recovered sites."""
    pos = 0
    slen = len(s)
    while pos <= slen:
        m = pat.search(s, pos)
        if m is None:
            break
        start = m.start()
        yield start
        pos = start + 1   # +1 (not m.end()) keeps overlapping tandem hits


def _scan_restriction_sites(
    seq: str,
    min_recognition_len: int = 6,
    unique_only: bool = True,
    circular: bool = True,
    allowed_enzymes: "frozenset[str] | None" = None,
) -> list[dict]:
    """Cached entry point for restriction-site scans. Identical
    signature + return shape to the inner `_scan_restriction_sites_impl`;
    consults `_RESTR_SCAN_CACHE` first so a `r`-toggle on a 5 Mb record
    drops from ~3 s to ~5 ms after the first scan.

    `allowed_enzymes` (GH #13, 2026-05-14): when supplied, the scan
    restricts to ONLY these names — overrides the `min_recognition_len`
    / `unique_only` filters since the user has hand-picked the set.
    Pass `None` to use the standard filter rules."""
    # Cache key uses `hash(seq)` rather than `id(seq)`: id-based keys
    # break for transient strings (e.g. unit tests build a fresh seq
    # per test, GC'd between calls — CPython's allocator can hand the
    # same address to a later string and produce a stale cache hit).
    # CPython interns the string hash on the first call, so subsequent
    # scans of the same seq are still O(1). Same fix applied to
    # `_ENZYME_CUTS_CACHE`.
    allowed_key = (
        tuple(sorted(allowed_enzymes)) if allowed_enzymes else None
    )
    # `len(seq)` is folded into the key so two different sequences can only
    # ever share a slot if they collide on BOTH `hash()` AND length — a strictly
    # narrower (still statistically-negligible) collision cross-section at zero
    # memory cost (unlike full value-keying, which would pin a multi-Mb string
    # per slot; see `_ENZYME_CUTS_CACHE`).
    key = (hash(seq), len(seq), int(min_recognition_len),
           bool(unique_only), bool(circular), allowed_key)
    with _state._COMPUTE_CACHE_LOCK:
        hit = _state._RESTR_SCAN_CACHE.get(key)
        if hit is not None:
            # Move to end (LRU touch) so a steady-state "I'm scanning the
            # same record" stays hot.
            _state._RESTR_SCAN_CACHE.move_to_end(key)
            # Defensive copy (match `_enzyme_cuts`): a consumer that mutates
            # the returned resite/recut dicts must not poison the shared cache
            # entry for the next scan of this sequence.
            return [dict(d) for d in hit]
    result = _scan_restriction_sites_impl(
        seq, min_recognition_len, unique_only, circular,
        allowed_enzymes=allowed_enzymes,
    )
    with _state._COMPUTE_CACHE_LOCK:
        if len(_state._RESTR_SCAN_CACHE) >= _state._RESTR_SCAN_CACHE_MAX:
            _state._RESTR_SCAN_CACHE.popitem(last=False)
        _state._RESTR_SCAN_CACHE[key] = result
    return [dict(d) for d in result]


@_timed("op.scan_restriction", threshold_ms=25)
def _scan_restriction_sites_impl(
    seq: str,
    min_recognition_len: int = 6,
    unique_only: bool = True,
    circular: bool = True,
    *,
    allowed_enzymes: "frozenset[str] | None" = None,
) -> list[dict]:
    """Scan both strands; return resite + recut dicts for every hit.

    resite — the recognition sequence span (colored bar)
    recut  — the cut position (single-bp marker: ↓ above or ↑ below DNA)

    min_recognition_len — skip enzymes whose recognition site is shorter than this
                          (default 6 to reduce noise from 4-cutters)
    unique_only         — if True, only include enzymes that cut exactly once
                          (forward + reverse strand combined; default True)
    circular            — if True (default), recognition sequences that span
                          the origin (bp n-1 → bp 0) are also found. SpliceCraft
                          is a plasmid viewer, so circularity is on by default.

    Wrap-around resites are emitted as TWO pieces so the existing linear-span
    rendering in the map / sequence panel stays correct: one piece on the
    "tail" (start..n) with the enzyme label, one unlabeled piece on the
    "head" (0..tail_len). The single-bp recut marker is placed at its real
    absolute position modulo n.
    """
    seq_u = seq.upper()
    n = len(seq_u)
    # Per-enzyme results collected first so we can filter to unique cutters
    by_enzyme: dict[str, list[dict]] = {}
    seen: set[tuple[str, int, int]] = set()   # deduplicate palindromes

    # For circular sequences, scan an augmented copy that includes up to
    # site_len-1 bp re-attached from the beginning so matches starting near
    # the end (that would otherwise be truncated) are found too.
    max_site_len = max((e[2] for e in _state._scan_catalog_hook()), default=0)
    # `max_site_len > 1` guard mirrors `_enzyme_cuts_impl`: at 0 (empty catalog)
    # the slice would be `seq_u[:-1]`, appending the whole sequence bar its last
    # base — a pointless O(n) copy of a ~2n string. Harmless today (an empty
    # catalog scans no enzymes) but exactly the asymmetry that bites a later
    # refactor or a catalog-monkeypatching test.
    if circular and n > 0 and max_site_len > 1:
        scan_seq = seq_u + seq_u[: max_site_len - 1]
    else:
        scan_seq = seq_u

    def _emit_resite(hits, p, site_len, strand, color, name,
                     cut_col, ext_cut_bp,
                     top_cut_bp=-1, bottom_cut_bp=-1):
        """Emit one or two resite dicts depending on wrap. Labels only on the
        first piece so the map doesn't double-print. For wrapped sites, the
        cut_col / ext_cut_bp fields are only meaningful on the piece that
        actually contains the cut; we attach them to the tail piece by default
        and clear them on the head piece.

        `top_cut_bp` / `bottom_cut_bp` are absolute top-strand-coordinate
        positions where the enzyme cleaves each strand. They're stored on
        every piece (including wrap continuations) so a click anywhere on
        the bar can render the per-strand cut split.
        """
        common = {
            "top_cut_bp":    top_cut_bp,
            "bottom_cut_bp": bottom_cut_bp,
        }
        if p + site_len <= n:
            hits.append({
                "type":       "resite",
                "start":      p,
                "end":        p + site_len,
                "strand":     strand,
                "color":      color,
                "label":      name,
                "cut_col":    cut_col,
                "ext_cut_bp": ext_cut_bp,
                # Full recognition span (== start/end when it doesn't
                # wrap). Carried explicitly so the click-highlight knows
                # the TRUE recognition bounds even on a wrapped piece.
                "rec_start":  p,
                "rec_end":    p + site_len,
                **common,
            })
            return
        # Wraps origin: tail [p, n) + head [0, (p + site_len) - n).
        tail_len = n - p
        head_len = (p + site_len) - n
        # cut_col (bar-relative) maps to whichever piece actually contains the
        # cut. ext_cut_bp (absolute) is unrelated to the tail/head split — it's
        # only meaningful when cut_col is None (Type IIS cuts outside the
        # recognition sequence). Attach it to both pieces so the cut arrow is
        # drawn regardless of which chunk contains the external cut position;
        # the chunk-level `chunk_start <= ext_cut_bp < chunk_end` test makes
        # the render idempotent. Regression guard added 2026-04-13.
        tail_cut_col = cut_col if (cut_col is not None and cut_col < tail_len) else None
        head_cut_col = ((cut_col - tail_len) if (cut_col is not None and cut_col >= tail_len)
                        else None)
        # Both pieces carry the SAME wrap-encoded recognition span
        # `[p, head_len)` (rec_end < rec_start signals the origin wrap)
        # so clicking the labeled tail can still highlight the full
        # recognition — tail bases [p, n) AND head bases [0, head_len).
        # Pre-2026-05-30 the tail piece reported `rec_end = n`, so the
        # click-highlight coloured only the pre-origin bases ("too few
        # purple bases" on a near-origin BsaI / Esp3I site).
        hits.append({
            "type":       "resite",
            "start":      p,
            "end":        n,
            "strand":     strand,
            "color":      color,
            "label":      name,
            "cut_col":    tail_cut_col,
            "ext_cut_bp": ext_cut_bp,
            "rec_start":  p,
            "rec_end":    head_len,
            **common,
        })
        hits.append({
            "type":       "resite",
            "start":      0,
            "end":        head_len,
            "strand":     strand,
            "color":      color,
            "label":      "",     # unlabeled continuation
            "cut_col":    head_cut_col,
            "ext_cut_bp": ext_cut_bp,
            "rec_start":  p,
            "rec_end":    head_len,
            **common,
        })

    for entry in _state._scan_catalog_hook():
        name, site, site_len, fwd_cut, rev_cut, color, pat, is_palindrome, rc_pat = entry
        # `allowed_enzymes` (GH #13): when set, the user has hand-
        # picked a working list — bypass the length filter so a 4-cutter
        # like Sau3AI shows up even when the global `restr_min_len`
        # is 6. We DO still honour `unique_only` below, since "unique
        # cutters of MY hand-picked list" is the common intent.
        if allowed_enzymes is not None:
            if name not in allowed_enzymes:
                continue
        else:
            if site_len < min_recognition_len:
                continue
        hits: list[dict] = []

        # Forward strand scan (over augmented sequence if circular)
        for p in _iter_match_starts(pat, scan_seq):
            if p >= n:
                continue   # duplicate of match already found at p - n
            key = (name, p, 1)
            if key in seen:
                continue
            seen.add(key)
            # On LINEAR molecules a negative-cut enzyme (e.g. BaeI fwd_cut=-10)
            # matching within |cut| bp of the 5' end makes (p+fwd_cut)%n wrap
            # to a phantom cut near the 3' end that doesn't biologically exist.
            # Mirror the reverse-strand guard below: drop the hit when a raw cut
            # falls outside the molecule. Circular molecules wrap correctly.
            if not circular and (
                    (p + fwd_cut) < 0 or (p + fwd_cut) > n
                    or (p + rev_cut) < 0 or (p + rev_cut) > n):
                continue
            # ext_cut_bp: absolute cut position when cut falls outside recognition
            _ext = ((p + fwd_cut) % n) if (fwd_cut <= 0 or fwd_cut >= site_len) else None
            _cc  = fwd_cut if 0 < fwd_cut < site_len else None
            # For forward-strand binding (whether palindromic or Type IIS),
            # fwd_cut counts from the recognition's 5' end on the top strand
            # and rev_cut counts from the recognition's 3' end on the bottom
            # strand (= 5' end of bottom strand reading right-to-left). Both
            # measured in top-strand coordinates: top cut at p+fwd_cut,
            # bottom cut at p+rev_cut. For palindromes these are mirror
            # images (rev_cut == site_len - fwd_cut); for Type IIS like BsaI
            # both fall outside the recognition.
            _top_cut = (p + fwd_cut) % n if n > 0 else 0
            _bot_cut = (p + rev_cut) % n if n > 0 else 0
            _emit_resite(hits, p, site_len, 1, color, name, _cc, _ext,
                         top_cut_bp=_top_cut, bottom_cut_bp=_bot_cut)
            hits.append({
                "type":   "recut",
                "start":  _top_cut,
                "end":    _top_cut + 1,
                "strand": 1,
                "color":  color,
                "label":  name,
            })

        # Reverse strand handling — uses precomputed is_palindrome and rc_pat
        # from _SCAN_CATALOG (no per-call _rc / _iupac_pattern work).
        if not is_palindrome:
            # Non-palindromic: scan for RC on forward strand to find
            # reverse-strand binding sites at their correct positions.
            for p in _iter_match_starts(rc_pat, scan_seq):
                if p >= n:
                    continue   # duplicate of match already found at p - n
                key = (name, p, -1)
                if key in seen:
                    continue
                seen.add(key)
                # Cut column within the bar: enzyme's fwd_cut mapped to
                # the reversed orientation displayed on the forward strand.
                # Symmetry: a forward-strand cut at offset `c` from the
                # recognition's 5' end appears on a reverse-bound site at
                # offset `site_len - c` from the bar's left edge.
                rev_cut_col = site_len - fwd_cut
                # 2026-05-27 (audit-5 restriction M1): on LINEAR
                # molecules, a reverse-strand non-palindromic
                # enzyme whose Type IIS cut would land BEFORE the
                # 5' end (e.g. BsaI with rev_cut=11 matching at p=2
                # on a linear plasmid → cut at -3 wraps to n-3)
                # can produce a cut position that doesn't biologically
                # exist. Drop the hit on linear molecules when the
                # cut would wrap. Circular molecules wrap correctly.
                _top_cut_raw = p + site_len - rev_cut
                _bot_cut_raw = p + site_len - fwd_cut
                if not circular and (
                        _top_cut_raw < 0 or _top_cut_raw > n
                        or _bot_cut_raw < 0 or _bot_cut_raw > n):
                    continue
                _top_cut_bp = _top_cut_raw % n   # top-strand cut in fwd coords
                _bot_cut_bp = _bot_cut_raw % n if n > 0 else 0
                _top_cut_outside = ((_top_cut_bp - p) % n) >= site_len
                _cc  = rev_cut_col if 0 <= rev_cut_col < site_len else None
                _ext = _top_cut_bp if _top_cut_outside else None
                _emit_resite(hits, p, site_len, -1, color, name, _cc, _ext,
                             top_cut_bp=_top_cut_bp,
                             bottom_cut_bp=_bot_cut_bp)
                hits.append({
                    "type":   "recut",
                    "start":  _bot_cut_bp,
                    "end":    _bot_cut_bp + 1,
                    "strand": -1,
                    "color":  color,
                    "label":  name,
                })

        if hits:
            by_enzyme[name] = hits

    feats: list[dict] = []
    # `placed` tracks (start, end, recognition_site) — keying on the
    # recognition string in addition to the span means HF / iso
    # variants of the SAME enzyme (e.g., EcoRI vs EcoRI-HF, both with
    # site "GAATTC") still collapse, but two enzymes with DIFFERENT
    # recognition patterns whose hits happen to land on the same
    # bp range (e.g., AccI/GTMKAC and BstZ17I/GTATAC both matching
    # GTATAC at one position) stay independent. Pre-2026-05-18 this
    # was `set[tuple[int, int]]` which over-collapsed and caused
    # `unique_only=True`/`False` to disagree on which enzyme to
    # surface when the catalog-order winner was a multi-cutter
    # filtered out by `unique_only=True` but kept by `unique_only=False`.
    placed: set[tuple[int, int, str]] = set()
    site_of: dict[str, str] = {
        entry[0]: entry[1] for entry in _state._scan_catalog_hook()
    }
    # When a custom enzyme list is active, the user has hand-picked
    # the set — surface every hit of those enzymes regardless of
    # cut-count. The `unique_only` filter is a discovery aid for the
    # default "show me 6+ bp unique cutters" workflow; it actively
    # hides multi-cutters (BsaI in a Golden Gate plasmid, EcoRI in a
    # repeat-laden synthetic construct) which is the opposite of what
    # the user wants after they typed those enzymes in by hand.
    effective_unique_only = (
        unique_only and allowed_enzymes is None
    )
    for name, hits in by_enzyme.items():
        # Count LABELED resites only — a wrap-around hit is emitted as one
        # labeled piece + one unlabeled continuation, but counts as 1 site.
        n_sites = sum(
            1 for h in hits if h["type"] == "resite" and h.get("label")
        )
        if effective_unique_only and n_sites != 1:
            continue
        # Skip isoschizomers / HF-variants of the SAME recognition that
        # land on an already-placed site. The recognition string is
        # part of the key so genuinely-different enzymes (different
        # recognition patterns) with accidental position overlap stay
        # independent.
        site_key = site_of.get(name, "")
        positions = {
            (h["start"], h["end"], site_key) for h in hits
            if h["type"] == "resite" and h.get("label")
        }
        if positions & placed:
            continue
        placed |= positions
        # Tag the per-enzyme cut-count on each labeled resite (NOT the label —
        # the label is a lookup key in Scrub / digest code). Renderers append a
        # superscript badge ("EcoRI²") from this when it's > 1, so the user
        # spots a multi-cutter at a glance and watches the count fall as they
        # edit out a site. Reactive for free: the scan re-runs on every record
        # change.
        if n_sites > 1:
            for h in hits:
                if h["type"] == "resite" and h.get("label"):
                    h["cut_count"] = n_sites
        feats.extend(hits)
    return feats


def _enzyme_cuts(seq: str, enzyme_names: list[str], *,
                  circular: bool = True) -> list[dict]:
    """Return all cuts on `seq` from the given enzymes, sorted by top
    cut position. Each entry is
    ``{top, bot, kind, overhang_seq, enzyme}`` where ``top`` and
    ``bot`` are absolute 0-based top-strand coords of the top-strand
    and bottom-strand cuts respectively.

    Unknown enzyme names are silently dropped (caller validates).
    Empty `enzyme_names` returns ``[]``.

    Results are LRU-cached on `(hash(seq), tuple(sorted enzymes), circular)`
    so that the Constructor's "Traditional" tab (which calls
    ``str(rec.seq).upper()`` afresh on every Simulate, allocating a
    new string object each time) still hits the cache on repeat
    clicks. Hash collisions are statistically irrelevant at the
    16-entry LRU size (~2^-32 per pair).

    Why not `id(seq)`: the old key only worked when the caller held
    onto the exact same string object across calls. The traditional-
    cloning flow doesn't, so the cache used to be permanently cold."""
    key = (hash(seq), len(seq), tuple(sorted(set(enzyme_names))), bool(circular))
    with _state._COMPUTE_CACHE_LOCK:
        hit = _state._ENZYME_CUTS_CACHE.get(key)
        if hit is not None:
            _state._ENZYME_CUTS_CACHE.move_to_end(key)
            # Defensive copy — callers occasionally mutate fragment dicts
            # (e.g., when filtering features), and a shared mutable list
            # would poison subsequent hits.
            return [dict(c) for c in hit]
    result = _enzyme_cuts_impl(seq, enzyme_names, circular=circular)
    with _state._COMPUTE_CACHE_LOCK:
        if len(_state._ENZYME_CUTS_CACHE) >= _state._ENZYME_CUTS_CACHE_MAX:
            _state._ENZYME_CUTS_CACHE.popitem(last=False)
        _state._ENZYME_CUTS_CACHE[key] = [dict(c) for c in result]
    return result


def _enzyme_cuts_impl(seq: str, enzyme_names: list[str], *,
                        circular: bool = True) -> list[dict]:
    """Underlying scanner — see `_enzyme_cuts` for the cached entry
    point. Same signature + return shape; this one is the actual
    work."""
    n = len(seq)
    if n == 0 or not enzyme_names:
        return []
    seq_u = seq.upper()
    out: dict[tuple[int, int], dict] = {}
    # Combined catalog so user-added custom enzymes participate in
    # every cloning flow that funnels through `_enzyme_cuts`. Via the
    # `_state` getter so this works once the fn moves to the L0 sibling.
    catalog = _state._all_enzymes_hook()
    # Hoisted out of the per-enzyme loop: build the wrap-augmented scan string
    # ONCE, sized to the LONGEST recognition site among the named enzymes,
    # instead of rebuilding `seq_u + seq_u[:site_len-1]` (an O(n) string copy)
    # for every enzyme. A longer augmentation than a short-site enzyme needs is
    # harmless — the only extra regex matches land at p >= n and are dropped by
    # the `if p >= n: continue` guard below (each duplicates a match already
    # found at p-n), so the cut set stays byte-identical. Mirrors
    # `_scan_restriction_sites_impl`, which already computes max_site_len once.
    max_site_len = max(
        (len(catalog[e][0]) for e in enzyme_names if e in catalog),
        default=0,
    )
    if circular and n > 0 and max_site_len > 1:
        scan_seq = seq_u + seq_u[: max_site_len - 1]
    else:
        scan_seq = seq_u
    for ename in enzyme_names:
        if ename not in catalog:
            continue
        site, fwd_cut, rev_cut = catalog[ename]
        site_u   = site.upper()
        site_len = len(site_u)
        try:
            pat      = _iupac_pattern(site_u)
            rc_site  = _rc(site_u)
        except ValueError:
            # Sweep #30 (2026-05-28): a corrupt custom-enzyme site (e.g. a
            # hand-edited custom_enzymes.json with a non-IUPAC char) makes
            # _iupac_pattern raise. Skip the bad enzyme rather than
            # crashing the whole digest — mirrors the per-enzyme guard
            # INV-85 added to _rebuild_scan_catalog. Pre-fix one bad site
            # silently killed the trad-cloning Simulate / MoClo workers.
            # (A valid site's reverse complement is always valid IUPAC, so
            # the rc_pat compile below can't raise once this passed.)
            _log.warning(
                "enzyme_cuts: skipping %r — bad recognition site %r",
                ename, site,
            )
            continue
        is_pal   = (rc_site == site_u)

        def _emit(top_bp_raw: int, bot_bp_raw: int):
            # Use raw (pre-modulo) values for kind detection AND to find
            # the overhang's earlier-cut anchor — post-modulo, a cut
            # that crosses the origin can flip the top<bot ordering
            # (e.g., raw top=99, raw bot=103, n=100 → post-mod top=99,
            # bot=3, which would falsely suggest a 3' overhang). The
            # raw difference equals |fwd_cut - rev_cut| which is fixed
            # by the enzyme.
            top_bp = top_bp_raw % n
            bot_bp = bot_bp_raw % n
            overhang_len = abs(top_bp_raw - bot_bp_raw)
            oh_start = (top_bp if top_bp_raw <= bot_bp_raw else bot_bp)
            oh_end   = (oh_start + overhang_len) % n if n else 0
            if overhang_len == 0:
                overhang = ""
            elif oh_end > oh_start:
                overhang = seq_u[oh_start:oh_end]
            else:
                # Origin-wrap: overhang region crosses bp 0.
                overhang = seq_u[oh_start:] + seq_u[:oh_end]
            kind = ("blunt" if top_bp_raw == bot_bp_raw
                    else "5'" if top_bp_raw < bot_bp_raw
                    else "3'")
            # Collapse COINCIDENT cuts: two enzymes that sever the identical
            # phosphodiester bond — isoschizomers (EcoRI / EcoRI-HF) or the GATC
            # trio (DpnII / Sau3AI / MboI) — make ONE physical cut, not two.
            # Keying on (top, bot, enzyme) emitted a duplicate cut per enzyme,
            # and `_fragments_from_cuts` then sliced the zero-gap span into a
            # phantom whole-plasmid fragment (a single-cut digest reported two
            # fragments of [n, n]). Key on the bond alone and merge the enzyme
            # attribution so fragment topology is correct while every enzyme
            # that makes the cut is still named.
            key = (top_bp, bot_bp)
            prev = out.get(key)
            if prev is not None:
                names = prev["enzyme"].split("/")
                if ename not in names:
                    prev["enzyme"] = "/".join(names + [ename])
                return
            out[key] = {
                "top":          top_bp,
                "bot":          bot_bp,
                "kind":         kind,
                "overhang_seq": overhang,
                "enzyme":       ename,
            }

        for p in _iter_match_starts(pat, scan_seq):
            if p >= n:
                continue
            # Linear molecules: skip a negative-cut enzyme whose cut would wrap
            # past the 5' end into a phantom 3'-end fragment boundary (mirrors
            # the restriction-overlay scan guard). Circular wraps correctly.
            if not circular and (
                    (p + fwd_cut) < 0 or (p + fwd_cut) > n
                    or (p + rev_cut) < 0 or (p + rev_cut) > n):
                continue
            _emit(p + fwd_cut, p + rev_cut)
        if not is_pal:
            rc_pat = _iupac_pattern(rc_site)
            for p in _iter_match_starts(rc_pat, scan_seq):
                if p >= n:
                    continue
                # On a reverse-strand binding, the cut positions mirror
                # around the recognition midpoint. Top cut on the bound
                # site (= bottom strand of unbound) is `site_len - rev_cut`
                # bases from the recognition's 5' end on the unbound
                # forward strand; bottom cut is `site_len - fwd_cut`.
                _rev_top_raw = p + site_len - rev_cut
                _rev_bot_raw = p + site_len - fwd_cut
                # Linear molecules: a reverse-strand Type IIS cut that falls
                # past the 5' end would wrap via _emit's `% n` into a phantom
                # 3'-end fragment boundary. Mirror the forward-path guard above
                # (and the restriction-overlay reverse-strand guard) and drop
                # it; circular molecules wrap correctly.
                if not circular and (
                        _rev_top_raw < 0 or _rev_top_raw > n
                        or _rev_bot_raw < 0 or _rev_bot_raw > n):
                    continue
                _emit(_rev_top_raw, _rev_bot_raw)
    return sorted(out.values(), key=lambda c: (c["top"], c["enzyme"]))


def _split_features_at_cuts(features: list[dict], n: int,
                              cut_top_positions: list[int],
                              circular: bool) -> dict[int, list[dict]]:
    """Slot features into fragments based on cut positions. Returns
    ``{fragment_index: [features]}``. Features that span a cut are
    split into two halves (one in each adjacent fragment); wrap features
    in a circular input get the same treatment around the origin too.

    Fragment indexing matches ``_fragments_from_cuts`` ordering:
      - circular: `i` ranges over `[0, len(cuts))`, fragment `i` runs from
        `cuts[i].top` to `cuts[(i+1) % n_cuts].top`.
      - linear: `i` ranges over `[0, len(cuts)+1)`, fragment 0 starts at 0
        and fragment `len(cuts)` ends at `n`.
    """
    if not features:
        return {}
    n_cuts = len(cut_top_positions)
    if n_cuts == 0:
        # All features in fragment 0.
        return {0: list(features)}

    # Wrap features (end < start on a circular input) need to be
    # split into a tail half [start, n) + a head half [0, end)
    # before the slotting algorithm runs. The latter assumes
    # start ≤ end, so a wrap feature would otherwise route through
    # `_slot_for(end-1)` for a position BEFORE the wrap, mis-
    # slotting the head half into a non-adjacent fragment. Each
    # half is tagged with `_wrap_origin_split` so a downstream
    # caller can rejoin them if it cares about origin-spanning
    # annotations on the result.
    expanded: list[dict] = []
    for f in features:
        try:
            fs = int(f.get("start", 0))
            fe = int(f.get("end",   0))
        except (TypeError, ValueError):
            continue
        if circular and fe < fs and 0 <= fe and fs <= n:
            expanded.append({**f, "start": fs, "end": n,
                              "_wrap_origin_split": "tail"})
            expanded.append({**f, "start": 0, "end": fe,
                              "_wrap_origin_split": "head"})
        else:
            expanded.append(f)
    features = expanded

    def _slot_for(bp: int) -> int:
        """Return the fragment index containing `bp` (0-based)."""
        if circular:
            # Find the cut k such that bp ∈ (cuts[k], cuts[k+1]] (mod n).
            for i in range(n_cuts):
                a = cut_top_positions[i]
                b = cut_top_positions[(i + 1) % n_cuts]
                if a < b:
                    if a <= bp < b:
                        return i
                else:
                    # Wrap fragment crosses origin
                    if bp >= a or bp < b:
                        return i
            return 0
        # Linear: cut at position c → bases [c, next_c) belong to fragment k+1.
        for i, c in enumerate(cut_top_positions):
            if bp < c:
                return i
        return n_cuts

    # Absolute span of each fragment, as a list of linear intervals (the
    # origin-crossing fragment needs two; a single circular cut makes ONE
    # fragment that covers the whole molecule).
    spans: list[list[tuple[int, int]]] = []
    if circular:
        for i in range(n_cuts):
            a = cut_top_positions[i]
            b = cut_top_positions[(i + 1) % n_cuts]
            spans.append([(a, b)] if a < b else [(a, n), (0, b)])
    else:
        bounds = [0] + list(cut_top_positions) + [n]
        spans = [[(bounds[i], bounds[i + 1])]
                 for i in range(len(bounds) - 1)]

    out: dict[int, list[dict]] = {}
    for f in features:
        s = int(f.get("start", 0))
        e = int(f.get("end",   0))
        if e == s:
            # Zero-length marker — no span to intersect; slot by position.
            slot = _slot_for(s)
            out.setdefault(slot, []).append(dict(f))
            continue
        # Slot by OVERLAP, not by endpoint. Endpoint-slotting emitted a piece
        # into the fragment holding `start` and the one holding `end-1` and
        # NOTHING else, which lost annotation in two ways:
        #   (1) a feature LONGER than one fragment covers the fragments
        #       BETWEEN its ends entirely — a 230 bp gene over three 100 bp
        #       fragments annotated only the two end pieces and left the
        #       middle fragment (100% inside the gene) bare;
        #   (2) on a circular molecule a feature long enough to run past every
        #       cut lands BOTH endpoints back in the SAME fragment, which
        #       endpoint-slotting read as "contained in one fragment" — so the
        #       gene vanished from every other fragment it actually covers.
        # Intersecting the feature's span with each fragment's span reports
        # both correctly and is what "which fragments carry this gene?" means.
        cut_inside = any(s < c < e for c in cut_top_positions)
        for i, ivs in enumerate(spans):
            if not any(max(s, lo) < min(e, hi) for lo, hi in ivs):
                continue
            starts_here = any(lo <= s < hi for lo, hi in ivs)
            ends_here   = any(lo <= e - 1 < hi for lo, hi in ivs)
            piece = dict(f)
            if starts_here and ends_here:
                # Both ends in this fragment — BUT if a cut falls STRICTLY
                # inside the feature, an excised piece was removed from its
                # middle (the classic case: cloning into lacZα's MCS excises
                # the stuffer between two cuts that both sit inside lacZα).
                # The remnant rides whole here, yet the gene is DISRUPTED —
                # tag it ``_split="whole"`` so the labeller flags it, even
                # though it wasn't split across two fragments.
                if cut_inside:
                    piece["_split"] = "whole"
            elif starts_here:
                piece["_split"] = "head"
            elif ends_here:
                piece["_split"] = "tail"
            else:
                # Fragment lies wholly INSIDE the feature. `_fragments_from_cuts`
                # clamps the piece to the fragment's full extent.
                piece["_split"] = "mid"
            out.setdefault(i, []).append(piece)
    return out


def _fragments_from_cuts(seq: str, cuts: list[dict], *,
                          circular: bool,
                          features: "list[dict] | None" = None,
                          source_label: str = "") -> list[dict]:
    """Slice `seq` at the given cut positions into a list of Fragment dicts.

    For circular input with ≥1 cut: `len(cuts)` fragments arranged around
    the origin. For circular input with 0 cuts: 1 fragment (the whole
    plasmid linearised at position 0; both ends marked "linear" — the
    caller should usually error out before reaching this case).

    For linear input: `len(cuts) + 1` fragments; the leftmost fragment
    starts at 0, the rightmost ends at `n`. End edges of these are
    marked `kind="linear"` since there's no enzyme there.

    Each fragment's `top_seq` is the contiguous top-strand slice from
    one cut's `top` to the next cut's `top` (or origin / endpoint).
    Features are slotted via `_split_features_at_cuts` and shifted into
    fragment-local 0-based coords."""
    n = len(seq)
    if n == 0:
        return []
    if not cuts:
        if circular:
            return [{
                "top_seq": seq,
                "left":  {"overhang_seq": "", "kind": "linear", "enzyme": ""},
                "right": {"overhang_seq": "", "kind": "linear", "enzyme": ""},
                "features": [dict(f) for f in (features or [])],
                "source_label": source_label,
            }]
        return [{
            "top_seq": seq,
            "left":  {"overhang_seq": "", "kind": "linear", "enzyme": ""},
            "right": {"overhang_seq": "", "kind": "linear", "enzyme": ""},
            "features": [dict(f) for f in (features or [])],
            "source_label": source_label,
        }]
    # De-dup COINCIDENT top-strand cuts. Two enzymes severing the same top bond
    # with different overhangs stay as separate cuts (`_enzyme_cuts` merges only
    # on the full (top,bot) bond). Sliced by top position, a coincident pair
    # drives the circular `a==b` branch below to emit a phantom FULL-LENGTH
    # fragment for BOTH cuts, and `_slot_for` then mis-slots every feature into
    # one — a wrong digest/clone. Collapse to one boundary per top position.
    # Guarded so the normal all-distinct case is byte-for-byte unchanged (no
    # re-sort): only the already-degenerate coincident case is normalised.
    if len({c["top"] for c in cuts}) != len(cuts):
        seen_tops = set()
        deduped = []
        for c in sorted(cuts, key=lambda c: c["top"]):
            if c["top"] not in seen_tops:
                seen_tops.add(c["top"])
                deduped.append(c)
        cuts = deduped
    cut_tops = [c["top"] for c in cuts]
    feat_slots = _split_features_at_cuts(features or [], n, cut_tops,
                                           circular=circular)
    fragments: list[dict] = []
    if circular:
        for i, c in enumerate(cuts):
            nxt = cuts[(i + 1) % len(cuts)]
            a, b = c["top"], nxt["top"]
            if a < b:
                top_seq = seq[a:b]
                offset  = a
            else:
                top_seq = seq[a:] + seq[:b]
                offset  = a   # fragment-local coord = (abs - a) % n
            frag_len = len(top_seq)
            local_feats: list[dict] = []
            for f in feat_slots.get(i, []):
                fs = int(f.get("start", 0))
                fe = int(f.get("end",   0))
                # Shift into fragment-local.
                if a < b:
                    new_s = fs - offset
                    new_e = fe - offset
                else:
                    # Origin-crossing fragment — the local map is MODULAR, and
                    # two things go wrong with a naive `(x - offset) % n`:
                    #
                    # (1) A slotted piece may overhang this fragment (the head/
                    #     tail halves emitted by `_split_features_at_cuts` keep
                    #     the feature's FULL span; the clamp below is what trims
                    #     them). Modularly, any absolute bp outside
                    #     `[a, n) ∪ [0, b)` maps ABOVE `frag_len` — so the plain
                    #     `min(...)` clamp pinned an overhanging START to the
                    #     fragment's END and emitted an INVERTED feature
                    #     (`start > end`). Downstream that renders as a phantom
                    #     near-full-length annotation on the cloning product
                    #     (a 49 bp AmpR remnant shipped as a 489 bp "AmpR").
                    #     Overhang is only ever "began before this fragment"
                    #     (start → 0) or "ran past it" (end → frag_len).
                    #
                    # (2) The END needs the half-open form `((e-a-1) % n) + 1`:
                    #     the naive `(e-a) % n` sends an end landing exactly on
                    #     the fragment's 3' edge to 0 instead of `frag_len`
                    #     (reachable on the single-cut `a == b` linearisation).
                    new_s = (fs - offset) % n
                    if new_s >= frag_len:
                        new_s = 0
                    if fe == fs:
                        new_e = new_s          # zero-length marker feature
                    else:
                        new_e = ((fe - offset - 1) % n) + 1
                        if new_e > frag_len:
                            new_e = frag_len
                # Clamp to fragment bounds.
                new_s = max(0, min(new_s, frag_len))
                new_e = max(0, min(new_e, frag_len))
                if new_s > new_e:
                    # The piece straddles this fragment's OWN 5'/3' ends. Only
                    # reachable on the single-cut circular digest (`a == b`),
                    # where the lone cut falls inside the feature. The fragment
                    # is LINEAR, so a `start > end` "wrap" is meaningless here —
                    # emit the two real pieces so no annotated base is lost.
                    local_feats.append({**f, "start": new_s, "end": frag_len})
                    local_feats.append({**f, "start": 0,     "end": new_e})
                    continue
                local_feats.append({**f, "start": new_s, "end": new_e})
            fragments.append({
                "top_seq": top_seq,
                "left":  {"overhang_seq": c["overhang_seq"],
                           "kind": c["kind"],
                           "enzyme": c["enzyme"]},
                "right": {"overhang_seq": nxt["overhang_seq"],
                           "kind": nxt["kind"],
                           "enzyme": nxt["enzyme"]},
                "features": local_feats,
                "source_label": source_label,
            })
        return fragments
    # Linear: walk left → right
    boundaries = [0] + cut_tops + [n]
    for i in range(len(boundaries) - 1):
        a, b = boundaries[i], boundaries[i + 1]
        top_seq = seq[a:b]
        local_feats: list[dict] = []
        for f in feat_slots.get(i, []):
            fs = int(f.get("start", 0))
            fe = int(f.get("end",   0))
            local_feats.append({
                **f,
                "start": max(0, min(fs - a, len(top_seq))),
                "end":   max(0, min(fe - a, len(top_seq))),
            })
        left  = ({"overhang_seq": "", "kind": "linear", "enzyme": ""}
                 if i == 0 else
                 {"overhang_seq": cuts[i - 1]["overhang_seq"],
                  "kind":         cuts[i - 1]["kind"],
                  "enzyme":       cuts[i - 1]["enzyme"]})
        right = ({"overhang_seq": "", "kind": "linear", "enzyme": ""}
                 if i == len(boundaries) - 2 else
                 {"overhang_seq": cuts[i]["overhang_seq"],
                  "kind":         cuts[i]["kind"],
                  "enzyme":       cuts[i]["enzyme"]})
        fragments.append({
            "top_seq":      top_seq,
            "left":         left,
            "right":        right,
            "features":     local_feats,
            "source_label": source_label,
        })
    return fragments


def _digest_with_enzymes(seq: str, enzyme_names: list[str], *,
                          circular: bool = True,
                          features: "list[dict] | None" = None,
                          source_label: str = "") -> list[dict]:
    """One-call digest: cut `seq` with `enzyme_names`, return Fragments.

    Fragments are sorted in cut order around the molecule (or 5'→3' for
    linear). Caller passes the input's features in absolute 0-based
    coords; they're slotted + shifted onto the appropriate fragments.

    Empty `enzyme_names` (or all-unknown) returns the input as a single
    uncut fragment."""
    cuts = _enzyme_cuts(seq, enzyme_names, circular=circular)
    return _fragments_from_cuts(seq, cuts, circular=circular,
                                  features=features,
                                  source_label=source_label)

# What is intentionally NOT extracted (yet):
#
#   _translate_cds — depends on `_GENETIC_CODE` and Biopython's translate.
#   Could move but increases the import surface meaningfully.
#
#   _feat_bounds — touches Biopython's `SeqFeature` / `CompoundLocation`
#   types; extracting would force splicecraft_biology to import Biopython
#   eagerly at module load, costing ~250 ms of startup time on every
#   `splicecraft-cli` call (which only imports for type hints today).
#
#   _bp_in — a METHOD on PlasmidMap, not a module-level function.
#
# Future extractions are welcome but must pass the three-test rule
# in CONTRIBUTING.md: no PlasmidApp coupling, reduces complexity at
# the call site, every existing test passes unchanged.


# ── Shared pure scan/revcomp helpers (relocated from hub for codon/primer) ──
# IUPAC-aware forbidden-site scan + IUPAC reverse-complement. Pure (only
# _iupac_pattern / _IUPAC_COMP); used by the codon optimizer, primer design,
# and the scrub paths — an L0 biology home keeps them reachable from every
# sibling without a codon/primer dependency.
def _forbidden_hit_set(seq: str, patterns) -> set[tuple[str, int]]:
    """Return ``{(pattern, position)}`` for every occurrence of every
    pattern in *seq*.

    IUPAC-aware (2026-05-30 hardening): a degenerate recognition site
    (``GGWCC``, ``CCANNNNNNTGG``, ``RGATCY`` …) is matched against concrete
    A/C/G/T via `_iupac_pattern`, not as a literal substring — so forbidding
    a degenerate-site enzyme actually finds (and the caller removes) its real
    sites instead of silently matching nothing. For an exact A/C/G/T site
    this is identical to the previous `seq.find` scan. Overlapping
    occurrences are reported via a zero-width lookahead (codon swaps are 3→3
    bases, so positions stay stable, letting the caller compare before/after
    hit sets directly). A pattern that isn't valid IUPAC is skipped."""
    out: set[tuple[str, int]] = set()
    for p in patterns:
        try:
            pat = _iupac_pattern(p)
        except ValueError:
            continue
        for m in re.finditer(f"(?={pat.pattern})", seq):
            out.add((p, m.start()))
    return out


def _mut_revcomp(seq: str) -> str:
    return seq.upper().translate(_IUPAC_COMP)[::-1]


# ── General sequence-alignment primitives (moved from hub) ─────────────────
# IUPAC degenerate-base matching + fuzzy 3'-anchored subsequence search +
# circular slice. Shared by the Myers/pairwise aligner, annotation-transfer,
# and primer-check — an L0 home below all of them.
# IUPAC ambiguity → set of unambiguous bases the code can stand for.
# Used by alignment match/mismatch counting + segment-state derivation
# so a primer with `N` against an `A` target counts as a MATCH (the
# `N` is compatible with `A`) rather than a mismatch. Pre-2026-05-27
# the pairwise aligner used scalar match/mismatch scoring with strict
# character equality — consensus reads and ambiguity-aware primers
# scored artificially low. `U` maps to {T} so RNA-pasted-as-DNA still
# aligns (also normalised to T in `_normalize_dna_for_align`).
_IUPAC_BASE_SET: "dict[str, frozenset[str]]" = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


def _iupac_compatible(a: str, b: str) -> bool:
    """Return True when IUPAC codes ``a`` and ``b`` share at least one
    unambiguous base. ``"A" vs "A"`` → True; ``"A" vs "N"`` → True;
    ``"A" vs "C"`` → False; ``"R" vs "Y"`` → False (R={A,G}, Y={C,T});
    ``"R" vs "M"`` → True (both include `A`).

    Case-insensitive. Returns False for any character not in the IUPAC
    nucleotide alphabet — including the gap glyph ``"-"`` (gaps are
    handled BEFORE this check at every call site).
    """
    if not a or not b:
        return False
    sa = _IUPAC_BASE_SET.get(a.upper())
    sb = _IUPAC_BASE_SET.get(b.upper())
    if sa is None or sb is None:
        return False
    return bool(sa & sb)


# Hard cap on hits returned by `_search_subsequence` — a degenerate
# query (all-`N`, or a high mismatch budget) on a large plasmid could
# otherwise match nearly every position. The scan stops once the cap is
# reached; the Ctrl+F handler surfaces "first N" when it bites.
_SEQ_SEARCH_MAX_HITS = 10_000


# Matches the FIRST non-concrete base in a template; `None` ⇒ the template is
# pure A/C/G/T, the precondition for `_search_subsequence`'s regex fast path.
_NON_ACGT_RE = re.compile(r"[^ACGT]")


@functools.lru_cache(maxsize=64)
def _iupac_overlap_pattern(site: str) -> "re.Pattern[str]":
    """`_iupac_pattern` wrapped in a zero-width lookahead so `finditer` yields
    EVERY start — including overlapping repeats (the three starts of ``"ATA"``
    in ``"ATATA"``) that a plain non-overlapping `finditer` drops. Used only by
    `_search_subsequence`'s exact fast path; LRU-cached because find-next
    re-queries the same needle on repeat."""
    return re.compile("(?=(" + _iupac_pattern(site).pattern + "))")


def _search_subsequence(
    seq: str,
    query: str,
    *,
    max_mismatches: int = 0,
    circular: bool = False,
    both_strands: bool = True,
    max_hits: int = _SEQ_SEARCH_MAX_HITS,
) -> "list[dict]":
    """Find every position where ``query`` occurs in ``seq`` within
    ``max_mismatches`` substitutions (Hamming distance — NO indels).

    IUPAC-aware: each base pair is compared with `_iupac_compatible`, so
    an ``N`` (or any ambiguity code) in either the query or the template
    matches every base it admits. Both inputs are scrubbed + validated
    through `_normalize_dna_for_align` (uppercase, ``U``→``T``, FASTA /
    whitespace / digits stripped); a foreign character propagates its
    ``ValueError`` so the caller can show the offending char.

    Returns hits sorted by start (then ``+`` before ``-``, then fewest
    mismatches), each::

        {"start": int, "end": int, "strand": "+"|"-", "mismatches": int}

    Coordinates are FORWARD (top-strand) for BOTH strands: a ``-`` hit at
    ``[start, end)`` means the query's reverse-complement sits on the top
    strand there (the query anneals to the BOTTOM strand across that
    span) — the same forward-coordinate convention the restriction scan
    uses for bottom-strand cuts (sacred invariant #2). ``start`` is always
    in ``[0, len(seq))``; ``end`` can exceed ``len(seq)`` for a hit that
    wraps the origin (circular only).

    circular: also find matches spanning the origin (scans
              ``seq + seq[:m-1]``). Leave False for a linear record.
    both_strands: also search the reverse-complement of the query, so a
                  motif is found on either strand (handy for annotating
                  inverted / direct repeats). A hit that lands on the same
                  span on both strands collapses to one ``+`` entry.
    max_hits: hard cap; the scan stops once reached. The caller should
              surface truncation when ``len(result) == max_hits``.

    Empty query, query longer than the sequence, or
    ``max_mismatches < 0`` → ``[]``.
    """
    if max_mismatches < 0:
        return []
    # Scrub + validate (may raise ValueError on a foreign char — let it
    # propagate; the modal turns it into a status-line message).
    s = _normalize_dna_for_align(seq or "")
    q = _normalize_dna_for_align(query or "")
    # Coordinate-safety guard (2026-06-02): this matcher returns indices
    # INTO the normalised haystack `s`, which the caller maps straight
    # back to a bp position on the loaded record (focus_span / select /
    # annotate). `_normalize_dna_for_align` only strips length-CHANGING
    # junk (whitespace / digits / FASTA headers); a loaded record's
    # sequence never carries those (U→T and upper-casing preserve
    # length), so `len(s) == len(seq)`. If a future caller ever feeds an
    # un-cleaned haystack, fail loud here rather than silently jump to /
    # select the WRONG base — coordinate drift in a search hit is the
    # same misleading-position class as a mis-drawn primer.
    if len(s) != len(seq or ""):
        raise ValueError(
            "search haystack changed length under normalisation "
            "(stray whitespace / digits / FASTA header?) — match "
            "coordinates would not line up with the sequence"
        )
    n = len(s)
    m = len(q)
    if m == 0 or n == 0 or m > n:
        return []
    # k >= m makes every window a "match"; clamp so the budget is
    # meaningful and the early-exit still bounds work.
    k = min(int(max_mismatches), m)
    iupac_ok = _iupac_compatible
    # Exact searches (no mismatch budget) over a template of only concrete
    # bases are equivalent to a compiled-regex scan and run ~50-100× faster
    # (C-level finditer vs a per-char Python double loop). When the TEMPLATE
    # carries ambiguity codes the equivalence breaks — a concrete query base
    # must still match a template `N`/`R` via `_iupac_compatible`, which a
    # query-derived regex can't express — so that case (and any k > 0) keeps
    # the Python loop. The match SET is identical either way (proven in
    # test_seq_search_nav.py::test_regex_fastpath_matches_python_loop).
    fastpath = (k == 0) and (_NON_ACGT_RE.search(s) is None)

    def _scan(pattern: str, strand: str, out: "list[dict]") -> None:
        # Origin-wrapping windows for circular records: a start at p in
        # [0, n) covers [p, p+m), reading across the origin via the
        # appended head. Linear records stop at n-m.
        scan_text = (s + s[: m - 1]) if circular else s
        last_start = (n - 1) if circular else (n - m)
        if fastpath:
            # Lookahead keeps overlapping hits; finditer yields increasing
            # positions, so `out` is built in the same order as the loop below
            # and the downstream dedup/sort sees an identical hit set.
            for mo in _iupac_overlap_pattern(pattern).finditer(scan_text):
                p = mo.start()
                if p > last_start:        # defensive — lookahead can't exceed
                    break
                out.append({"start": p, "end": p + m,
                            "strand": strand, "mismatches": 0})
                if len(out) >= max_hits:
                    return
            return
        for p in range(0, last_start + 1):
            mm = 0
            ok = True
            for i in range(m):
                if not iupac_ok(pattern[i], scan_text[p + i]):
                    mm += 1
                    if mm > k:
                        ok = False
                        break
            if ok:
                out.append({
                    "start": p,
                    "end": p + m,
                    "strand": strand,
                    "mismatches": mm,
                })
                if len(out) >= max_hits:
                    return

    hits: "list[dict]" = []
    _scan(q, "+", hits)
    if both_strands and len(hits) < max_hits:
        rc = _rc(q)
        # A fully palindromic query (rc == q) would duplicate every
        # forward hit on the reverse pass — skip it.
        if rc != q:
            _scan(rc, "-", hits)

    # De-dup: a '+' and '-' hit on the SAME span (semi-palindromic
    # window) collapse to the '+' entry so find-next doesn't visit the
    # same bp twice.
    seen: "dict[tuple, dict]" = {}
    for h in hits:
        key = (h["start"], h["end"])
        prev = seen.get(key)
        if prev is None:
            seen[key] = h
        elif h["mismatches"] < prev["mismatches"]:
            # Fewer mismatches = the strictly better match — keep it regardless
            # of strand (pre-fix a '+' always won even when the '-' hit annealed
            # more tightly, mis-ranking mismatch-tolerant Ctrl-F results).
            seen[key] = h
        elif (h["mismatches"] == prev["mismatches"]
              and prev["strand"] == "-" and h["strand"] == "+"):
            # Tie on mismatches → prefer the '+' (forward) representation, so a
            # semi-palindromic span still collapses to one '+' entry. Exact
            # (k=0) search is unaffected: all mismatches are 0 ⇒ this tie branch.
            seen[key] = h
    deduped = list(seen.values())
    deduped.sort(
        key=lambda h: (h["start"], 0 if h["strand"] == "+" else 1,
                       h["mismatches"]),
    )
    return deduped[:max_hits]


def _circ_slice(seq: str, start: int, length: int, total: int) -> str:
    """`length` bases of `seq` starting at `start` (taken mod `total`),
    wrapping the origin as needed. `start` may be negative or >= total.
    Returns "" on degenerate input."""
    if total <= 0 or length <= 0 or not seq:
        return ""
    start %= total
    if start + length <= total:
        return seq[start:start + length]
    out = seq[start:]
    need = length - len(out)
    # Primers are short (<= _PCR_MAX_PRIMER_LEN) and the template is >= the
    # primer length at every call site, so this wraps at most once for real
    # inputs; the loop is belt-and-suspenders for a length > total.
    while need > 0:
        chunk = seq[:need]
        if not chunk:
            break
        out += chunk
        need -= len(chunk)
    return out[:length]
