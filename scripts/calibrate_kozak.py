#!/usr/bin/env python3
"""Train and re-measure the vendored plant translation-initiation matrix.

Emits `splicecraft_kozak_model.py` — the log-odds PWM that
`splicecraft_kozak` scores an ATG context with — and prints every number
recorded in its metadata. Nothing here writes to the SpliceCraft data dir and
nothing here imports the hub.

    python3 scripts/calibrate_kozak.py                     # fetch + measure
    python3 scripts/calibrate_kozak.py --cache DIR         # reuse a download
    python3 scripts/calibrate_kozak.py --emit splicecraft_kozak_model.py

**Why this script exists rather than a citation.** The shipped -3 purine /
+4 G rule is Kozak's 1987 mammalian work, and the first attempt to measure a
plant replacement for it OVERSTATED the gain by 70% (reported +0.09 AUC; the
honest number is +0.054). Three separate things did that, and each is now
closed off here:

  1. the train/test halves were cut from a list ordered BY SPECIES, so the
     model was trained on two species and tested on a third
  2. splice variants of one gene were counted as independent transcripts,
     and they share a start context
  3. the scoring window was chosen on the same split the AUC was reported
     from

So: one transcript per GENE, a 75/25 split BY GENE stratified by species, the
window chosen on an inner split of the TRAIN half only, and a gene-resampled
bootstrap CI on the difference that is the actual claim.

Two confounds are measured rather than assumed away:

  * **Coding composition.** A real start's downstream half is CODING while a
    5'UTR uATG's is UTR, so a window reaching past +4 discriminates by base
    composition alone -- worth +0.03 AUC of signal that says nothing about
    initiation. The window therefore STOPS at +4, and the grid below prints
    what +6 and +9 would have bought so the refusal is visible.
  * **Gene-predictor circularity.** Gnomon-predicted (XM_) start codons are
    partly a gene predictor's own Kozak-like prior. Curated (NM_) and
    predicted starts are scored apart; if the model only worked on XM_ it
    would be learning the predictor, not the plant.
"""
from __future__ import annotations

import argparse
import bisect
import json
import random
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import date
from io import StringIO
from math import log2
from pathlib import Path

UA = {"User-Agent": "SpliceCraft-kozak-calibration/1.0"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SEED = 20261031
BASES = "ACGT"
PER_SPECIES, POOL = 420, 4000
# Proportional smoothing, exactly as `splicecraft_splice` -- a flat Laplace +1
# over two populations of very different size makes an IMPOSSIBLE base outscore
# the right one at an invariant position. Here that matters at the three ATG
# columns, which must read 0.00 bits.
PSEUDO_FRAC = 0.001
# Grid searched on the inner dev split. `down` past 4 is reported but refused.
UP_GRID, DOWN_GRID = (6, 8, 10, 12, 15, 20), (4, 6, 9)

# 4 dicots + 3 monocots. Clade spread is the point: the splice matrices are
# clade-matched because intron recognition is NOT uniform across plants, so
# whether initiation context is uniform has to be measured. Leave-one-species-
# out below is that measurement.
SPECIES = {
    "ath": ("Arabidopsis thaliana", "dicot"),
    "gma": ("Glycine max", "dicot"),
    "sly": ("Solanum lycopersicum", "dicot"),
    "nta": ("Nicotiana tabacum", "dicot"),
    "osa": ("Oryza sativa Japonica Group", "monocot"),
    "zma": ("Zea mays", "monocot"),
    "hvu": ("Hordeum vulgare", "monocot"),
}


# ── fetch ─────────────────────────────────────────────────────────────────
def _get(path: str, **params) -> str:
    url = EUTILS + path + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers=UA)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read(300_000_000).decode("utf-8", "replace")
        except Exception as exc:                      # noqa: BLE001 - retry
            if attempt == 3:
                raise
            print(f"    retry after {exc}", file=sys.stderr)
            time.sleep(4)
    raise RuntimeError("unreachable")


def fetch(cache: Path) -> dict:
    """RefSeq plant mRNAs whose CDS start is a verified complete-ORF ATG.

    Every filter here is a QC check that fired at least once, and the counts
    are printed so a silent change in RefSeq is visible rather than absorbed."""
    from Bio import SeqIO
    from Bio.SeqFeature import ExactPosition

    drops: "Counter[str]" = Counter()
    out: dict = {}
    for tag, (organism, clade) in SPECIES.items():
        term = (f'"{organism}"[ORGN] AND biomol_mrna[PROP] AND '
                'srcdb_refseq[PROP] AND 600:4000[SLEN]')
        res = json.loads(_get("esearch.fcgi", db="nucleotide", term=term,
                              retmax=POOL, retmode="json"))["esearchresult"]
        ids = res["idlist"]
        # A SEEDED SAMPLE of the whole matching set, not the first page
        # esearch happens to order first.
        rng = random.Random(f"{SEED}:{tag}")
        if len(ids) > PER_SPECIES:
            ids = rng.sample(ids, PER_SPECIES)
        print(f"  {tag}: {res.get('count')} matching, sampled {len(ids)}",
              file=sys.stderr)
        for i in range(0, len(ids), 100):
            txt = _get("efetch.fcgi", db="nucleotide",
                       id=",".join(ids[i:i + 100]), rettype="gb",
                       retmode="text")
            for rec in SeqIO.parse(StringIO(txt), "genbank"):
                s = str(rec.seq).upper()
                if set(s) - set(BASES):
                    drops["non_acgt"] += 1; continue
                cdss = [f for f in rec.features if f.type == "CDS"]
                if not cdss:
                    drops["no_cds"] += 1; continue
                f = cdss[0]
                if f.location.strand != 1:
                    drops["minus_strand"] += 1; continue
                if {"pseudo", "pseudogene"} & set(f.qualifiers):
                    drops["pseudo"] += 1; continue
                if {"exception", "transl_except"} & set(f.qualifiers):
                    drops["exception"] += 1; continue
                # A `<1..` partial start is a truncated record, not a start.
                if not (isinstance(f.location.start, ExactPosition)
                        and isinstance(f.location.end, ExactPosition)):
                    drops["fuzzy_bounds"] += 1; continue
                c, e = int(f.location.start), int(f.location.end)
                if s[c:c + 3] != "ATG":
                    drops["not_atg"] += 1; continue
                if (e - c) % 3:
                    drops["not_in_frame"] += 1; continue
                if s[e - 3:e] not in ("TAA", "TAG", "TGA"):
                    drops["no_stop_codon"] += 1; continue
                if c < 1:
                    drops["no_five_utr"] += 1; continue
                q = f.qualifiers
                out[f"{tag}:{rec.id}"] = {
                    "species": tag, "clade": clade, "acc": rec.id,
                    "prefix": rec.id.split("_")[0],
                    "gene": (q.get("gene") or q.get("locus_tag")
                             or q.get("protein_id") or [rec.id])[0],
                    "seq": s, "cds": c, "cds_end": e,
                }
            time.sleep(0.4)
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "kozak_mrna.json").write_text(json.dumps(out))
    print(f"  -> {len(out)} transcripts; QC drops: {dict(drops)}",
          file=sys.stderr)
    return out


# ── model ─────────────────────────────────────────────────────────────────
def count_matrix(windows, width: int) -> "list[dict[str, float]]":
    cols = [{b: 0.0 for b in BASES} for _ in range(width)]
    used = 0
    for w in windows:
        if len(w) != width or set(w) - set(BASES):
            continue
        used += 1
        for i, b in enumerate(w):
            cols[i][b] += 1.0
    if not used:
        raise ValueError(f"no usable training windows of width {width}")
    pseudo = max(PSEUDO_FRAC * used, 1e-9) / len(BASES)
    for col in cols:
        for b in col:
            col[b] += pseudo
        total = sum(col.values())
        for b in col:
            col[b] /= total
    return cols


def train(real, decoy, width: int) -> "list[dict[str, float]]":
    p_real, p_dec = count_matrix(real, width), count_matrix(decoy, width)
    return [{b: log2(p_real[i][b] / p_dec[i][b]) for b in BASES}
            for i in range(width)]


def score(window: str, pwm, mask=None) -> float:
    return sum(pwm[i][b] for i, b in enumerate(window)
               if mask is None or mask[i])


def rule_score(window: str, up: int) -> int:
    """What the shipped -3 purine / +4 G rule computes, on the same window."""
    return int(window[up - 3] in "AG") + int(window[up + 3] == "G")


def auc(a, b) -> float:
    """P(score(positive) > score(decoy)), ties at 0.5 -- Mann-Whitney. The
    rule takes only 3 values, so tie handling is not optional here."""
    if not a or not b:
        return float("nan")
    sb = sorted(b); n = len(sb); total = 0.0
    for x in a:
        lo, hi = bisect.bisect_left(sb, x), bisect.bisect_right(sb, x)
        total += (lo + (hi - lo) / 2) / n
    return total / len(a)


def contexts(rec: dict, up: int, down: int):
    """(start context, 5'UTR uATG contexts, in-CDS ATG contexts).

    The 5'UTR set is the population that matters: it is exactly what a uORF
    screen has to rank. The in-CDS set is coding on BOTH sides, so it cannot
    be discriminated by the UTR/CDS composition shift."""
    s, c, e = rec["seq"], rec["cds"], rec["cds_end"]

    def ctx(p):
        return s[p - up:p + down] if p >= up and p + down <= len(s) else None

    utr, cds = [], []
    i = s.find("ATG")
    while i != -1:
        if i != c:
            w = ctx(i)
            if w is not None:
                (utr if i < c else cds if i < e else []).append(w)
        i = s.find("ATG", i + 1)
    return ctx(c), utr, cds


def one_per_gene(recs: dict) -> dict:
    """Splice variants of one gene share a start context; counting them twice
    inflates every number downstream. Deterministic pick per gene."""
    by_gene = defaultdict(list)
    for key, r in recs.items():
        by_gene[(r["species"], r["gene"])].append((key, r))
    return {gk: random.Random(f"{SEED}:{gk}").choice(sorted(lst))[1]
            for gk, lst in by_gene.items()}


def bundles(genes: dict, keys, up: int, down: int) -> list:
    """Per-gene bundles, so a bootstrap resamples GENES and carries each
    gene's own decoys with it. Duplicate start contexts are dropped: a
    paralogue pair sharing a window would otherwise straddle the split."""
    out, seen = [], set()
    for gk in keys:
        r = genes[gk]
        pos, utr, cds = contexts(r, up, down)
        if pos is None or pos in seen:
            continue
        seen.add(pos)
        out.append({"gene": gk, "species": r["species"], "clade": r["clade"],
                    "prefix": r["prefix"], "pos": pos, "utr": utr, "cds": cds})
    return out


def flat(bs, key):
    return [w for b in bs for w in b[key]]


# ── report + emit ─────────────────────────────────────────────────────────
MODULE_TEMPLATE = '''"""Vendored plant translation-initiation matrix. GENERATED -- do not hand-edit.

Regenerate with `python3 scripts/calibrate_kozak.py --emit
splicecraft_kozak_model.py`. Every number in KOZAK_MODEL_META is held out;
that script's docstring records what each one controls for.
"""

KOZAK_MODEL_META = {meta}

# Per-position log2-odds, real plant start codons vs the 5'UTR upstream ATGs
# they compete with. Window {window} with the A of the ATG at offset {up}.
#
# **This matrix reads GC-rich upstream and that is CORRECT.** It is a log-odds
# ratio, not a frequency table: real starts ARE A-rich upstream
# (`{consensus}` is the most frequent base per position), but the uATGs they
# are being discriminated FROM are even more AT-rich ({at_real:.3f} vs
# {at_decoy:.3f} A+T over the upstream window), so the ratio tilts the other
# way. Read `real_start_consensus` in the metadata for the biology; read this
# matrix only as a discriminator.
INITIATION_PWM = [
{pwm}
]

# Score quantiles of {n_quant} HELD-OUT real plant start codons, 0th to 100th
# percentile in steps of 5. Turns a score into "stronger than N% of real plant
# starts", which is the only calibrated statement this model can make.
REAL_START_QUANTILES = [
{quantiles}
]
'''


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=None,
                    help="directory holding/receiving kozak_mrna.json")
    ap.add_argument("--emit", default=None,
                    help="write the vendored model module to this path")
    ap.add_argument("--no-fetch", action="store_true",
                    help="fail rather than hit NCBI when the cache is cold")
    args = ap.parse_args()

    cache = Path(args.cache) if args.cache else Path("kozak-calibration")
    blob = cache / "kozak_mrna.json"
    if blob.exists():
        recs = json.loads(blob.read_text())
        print(f"reference set: {len(recs)} transcripts (cached {blob})")
    elif args.no_fetch:
        print(f"no cache at {blob} and --no-fetch given", file=sys.stderr)
        return 2
    else:
        print("fetching RefSeq plant mRNAs ...")
        recs = fetch(cache)

    genes = one_per_gene(recs)
    print(f"{len(recs)} transcripts -> {len(genes)} genes "
          f"({len(recs) - len(genes)} same-gene variants dropped)")
    print(f"  species: {dict(Counter(r['species'] for r in genes.values()))}")
    print(f"  prefix : {dict(Counter(r['prefix'] for r in genes.values()))}")

    # 75/25 by gene, stratified by species.
    rng = random.Random(SEED)
    train_g, test_g = [], []
    for sp in sorted({r["species"] for r in genes.values()}):
        g = sorted(gk for gk, r in genes.items() if r["species"] == sp)
        rng.shuffle(g)
        cut = int(round(0.75 * len(g)))
        train_g += g[:cut]; test_g += g[cut:]
    inner = list(train_g)
    random.Random(SEED + 1).shuffle(inner)
    icut = int(round(0.75 * len(inner)))
    fit_g, dev_g = inner[:icut], inner[icut:]
    print(f"split by gene: train {len(train_g)} / test {len(test_g)}  "
          f"(inner fit {len(fit_g)} / dev {len(dev_g)})")

    # ── window grid, decided on the inner dev split and NOWHERE else ──────
    print("\nwindow grid -- AUC real starts vs 5'UTR uATGs, inner dev split:")
    print(f"  {'window':>10s} {'full':>7s} {'upstream':>9s} {'+4 only':>8s}")
    grid = {}
    for up in UP_GRID:
        for down in DOWN_GRID:
            fb, db = (bundles(genes, fit_g, up, down),
                      bundles(genes, dev_g, up, down))
            pwm = train([b["pos"] for b in fb], flat(fb, "utr"), up + down)
            m_up = [i < up for i in range(up + down)]
            m_p4 = [i == up + 3 for i in range(up + down)]
            dp, du = [b["pos"] for b in db], flat(db, "utr")
            grid[(up, down)] = a = auc([score(w, pwm) for w in dp],
                                       [score(w, pwm) for w in du])
            print(f"  -{up:>2d}..+{down:<2d}   {a:>7.3f} "
                  f"{auc([score(w,pwm,m_up) for w in dp],[score(w,pwm,m_up) for w in du]):>9.3f} "
                  f"{auc([score(w,pwm,m_p4) for w in dp],[score(w,pwm,m_p4) for w in du]):>8.3f}")
    UP, DOWN = max((k for k in grid if k[1] == 4), key=lambda k: grid[k])
    leak = max(grid, key=lambda k: grid[k])
    print(f"\n  chosen (down=+4, no coding leak): -{UP}..+{DOWN} "
          f"dev AUC {grid[(UP, DOWN)]:.3f}")
    print(f"  REFUSED -{leak[0]}..+{leak[1]} at dev AUC {grid[leak]:.3f}: "
          f"+{grid[leak]-grid[(UP,DOWN)]:.3f} of that is the CODING/UTR "
          f"composition shift, not initiation context")
    W = UP + DOWN

    # ── the shipped matrix: fit on TRAIN, every number below from TEST ────
    tr = bundles(genes, train_g, UP, DOWN)
    te = bundles(genes, test_g, UP, DOWN)
    PWM = train([b["pos"] for b in tr], flat(tr, "utr"), W)
    m_up = [i < UP for i in range(W)]

    s_pos = [score(b["pos"], PWM) for b in te]
    s_utr = [score(w, PWM) for w in flat(te, "utr")]
    s_cds = [score(w, PWM) for w in flat(te, "cds")]
    r_pos = [rule_score(b["pos"], UP) for b in te]
    r_utr = [rule_score(w, UP) for w in flat(te, "utr")]
    r_cds = [rule_score(w, UP) for w in flat(te, "cds")]
    a_utr, a_cds = auc(s_pos, s_utr), auc(s_pos, s_cds)
    ar_utr, ar_cds = auc(r_pos, r_utr), auc(r_pos, r_cds)
    a_up = auc([score(b["pos"], PWM, m_up) for b in te],
               [score(w, PWM, m_up) for w in flat(te, "utr")])

    print(f"\n{'='*72}\nHELD OUT: {len(te)} genes, {len(s_utr)} 5'UTR uATGs, "
          f"{len(s_cds)} in-CDS ATGs\n{'='*72}")
    print(f"{'':30s} {'vs 5UTR uATG':>13s} {'vs in-CDS ATG':>14s}")
    print(f"{'-3/+4 rule (shipped today)':30s} {ar_utr:>13.3f} {ar_cds:>14.3f}")
    print(f"{'PWM, upstream positions only':30s} {a_up:>13.3f} {'':>14s}")
    print(f"{'PWM, full window':30s} {a_utr:>13.3f} {a_cds:>14.3f}")

    boot, brng = [], random.Random(SEED + 7)
    for _ in range(2000):
        samp = [te[brng.randrange(len(te))] for _ in range(len(te))]
        u = flat(samp, "utr")
        if not u:
            continue
        boot.append(
            auc([score(b["pos"], PWM) for b in samp], [score(w, PWM) for w in u])
            - auc([rule_score(b["pos"], UP) for b in samp],
                  [rule_score(w, UP) for w in u]))
    boot.sort()
    ci = [round(boot[int(0.025 * len(boot))], 3),
          round(boot[int(0.975 * len(boot))], 3)]
    print(f"\ndAUC over the rule: {a_utr - ar_utr:+.3f}  95% CI "
          f"[{ci[0]:+.3f}, {ci[1]:+.3f}]  (2000 gene-resamples)")

    # The capability claim: the rule has 3 buckets and cannot order WITHIN one.
    print("\nwhat the rule cannot do -- separation inside each rule bucket:")
    buckets = {}
    for b2 in (0, 1, 2):
        rp = [score(b["pos"], PWM) for b in te if rule_score(b["pos"], UP) == b2]
        rd = [score(w, PWM) for w in flat(te, "utr") if rule_score(w, UP) == b2]
        if len(rp) > 10 and len(rd) > 10:
            buckets[b2] = round(auc(rp, rd), 3)
            print(f"  rule score {b2}: n_real={len(rp):4d} n_uATG={len(rd):4d}"
                  f"   PWM AUC {buckets[b2]:.3f}")

    print("\nper-species held-out AUC:")
    per_species = {}
    for sp in sorted({b["species"] for b in te}):
        g = [b for b in te if b["species"] == sp]
        u = flat(g, "utr")
        if not u:
            continue
        per_species[sp] = round(auc([score(b["pos"], PWM) for b in g],
                                    [score(w, PWM) for w in u]), 3)
        print(f"  {sp} ({g[0]['clade']:7s}) n={len(g):3d}  PWM "
              f"{per_species[sp]:.3f}  rule "
              f"{auc([rule_score(b['pos'],UP) for b in g],[rule_score(w,UP) for w in u]):.3f}")

    print("\nleave-one-species-out (trained on the other six):")
    loso = {}
    allb = bundles(genes, train_g + test_g, UP, DOWN)
    for sp in sorted({b["species"] for b in allb}):
        out_b = [b for b in allb if b["species"] == sp]
        u = flat(out_b, "utr")
        if not u:
            continue
        P = train([b["pos"] for b in allb if b["species"] != sp],
                  flat([b for b in allb if b["species"] != sp], "utr"), W)
        loso[sp] = round(auc([score(b["pos"], P) for b in out_b],
                             [score(w, P) for w in u]), 3)
        print(f"  {sp} ({out_b[0]['clade']:7s}) n={len(out_b):3d}  PWM "
              f"{loso[sp]:.3f}  rule "
              f"{auc([rule_score(b['pos'],UP) for b in out_b],[rule_score(w,UP) for w in u]):.3f}")

    print("\ncurated vs gene-predictor starts (held out):")
    circ = {}
    for pref in sorted({b["prefix"] for b in te}):
        g = [b for b in te if b["prefix"] == pref]
        u = flat(g, "utr")
        if len(g) < 20 or not u:
            continue
        circ[pref] = round(auc([score(b["pos"], PWM) for b in g],
                               [score(w, PWM) for w in u]), 3)
        print(f"  {pref}_ n={len(g):3d}  PWM {circ[pref]:.3f}")

    # ── calibration: percentile table + the two tier cuts ─────────────────
    srt = sorted(s_pos)
    quant = [round(srt[min(len(srt) - 1, int(p / 100 * len(srt)))], 4)
             for p in range(0, 101, 5)]
    tiers = {}
    for name, pct in (("adequate", 5), ("strong", 50)):
        cut = srt[int(pct / 100 * len(srt))]
        tiers[name] = {
            "cut": round(cut, 4),
            "real_starts_at_or_above": round(
                sum(1 for s in s_pos if s >= cut) / len(s_pos), 3),
            "utr_uatgs_at_or_above": round(
                sum(1 for s in s_utr if s >= cut) / len(s_utr), 3),
            "cds_atgs_at_or_above": round(
                sum(1 for s in s_cds if s >= cut) / len(s_cds), 3),
        }
    print("\ntier cuts (anchored on the real-start distribution):")
    for name, t in tiers.items():
        print(f"  {name:9s} score >= {t['cut']:>7.3f}   keeps "
              f"{t['real_starts_at_or_above']:.1%} of real starts, admits "
              f"{t['utr_uatgs_at_or_above']:.1%} of 5'UTR uATGs")

    # Coverage: a uATG too close to the 5' end has no full window, and a
    # partial window's log-odds is NOT comparable to this distribution.
    tot_utr = scored_utr = 0
    for gk in train_g + test_g:
        r = genes[gk]
        s, c = r["seq"], r["cds"]
        i = s.find("ATG")
        while i != -1 and i < c:
            tot_utr += 1
            if i >= UP and i + DOWN <= len(s):
                scored_utr += 1
            i = s.find("ATG", i + 1)
    cover = round(scored_utr / max(1, tot_utr), 3)
    print(f"\ncoverage: {cover:.1%} of real 5'UTR uATGs sit far enough from "
          f"the 5' end to score ({scored_utr}/{tot_utr})")

    real_freq = count_matrix([b["pos"] for b in tr], W)
    dec_freq = count_matrix(flat(tr, "utr"), W)
    consensus = "".join(max(c, key=lambda b: c[b]) for c in real_freq)
    at_real = sum(real_freq[i]["A"] + real_freq[i]["T"] for i in range(UP)) / UP
    at_dec = sum(dec_freq[i]["A"] + dec_freq[i]["T"] for i in range(UP)) / UP
    print(f"\nmost frequent base per position, real starts: {consensus}")
    print(f"  upstream A+T  real {at_real:.3f}  vs uATG decoys {at_dec:.3f}")
    bits = [round(max(c.values()) - min(c.values()), 2) for c in PWM]
    print(f"  per-position bits: {bits}")
    if any(bits[UP + i] > 0.001 for i in range(3)):
        print("  !! the three ATG columns are NOT ~0 bits -- the window "
              "geometry or the pseudocounting is wrong", file=sys.stderr)
        return 3

    meta = {
        "window": f"-{UP}..+{DOWN}",
        "up": UP, "down": DOWN, "atg_offset": UP,
        "source": (f"{len(genes)} genes from {len(SPECIES)} plant species "
                   "(4 dicot, 3 monocot), RefSeq mRNAs with an exact, "
                   "complete-ORF, ATG-verified CDS start and an annotated "
                   "5'UTR; one transcript per gene"),
        "decoy": ("the 5'UTR upstream ATGs of the same transcripts -- the "
                  "population a uORF screen actually has to rank"),
        "fetched": date.today().isoformat(),
        "n_genes_train": len(tr), "n_genes_test": len(te),
        "n_utr_uatg_test": len(s_utr), "n_cds_atg_test": len(s_cds),
        "auc_vs_utr_uatg": round(a_utr, 3),
        "auc_vs_cds_atg": round(a_cds, 3),
        "auc_upstream_only": round(a_up, 3),
        "rule_auc_vs_utr_uatg": round(ar_utr, 3),
        "rule_auc_vs_cds_atg": round(ar_cds, 3),
        "delta_auc_over_rule": round(a_utr - ar_utr, 3),
        "delta_auc_ci95": ci,
        "auc_within_rule_bucket": buckets,
        "auc_per_species": per_species,
        "auc_leave_one_species_out": loso,
        "auc_by_record_class": circ,
        "tiers": tiers,
        "uatg_coverage": cover,
        "real_start_consensus": consensus,
        "upstream_at_fraction": {"real_starts": round(at_real, 3),
                                 "utr_uatgs": round(at_dec, 3)},
        "position_bits": bits,
        "refused": {
            "window_past_plus4": (
                f"-{leak[0]}..+{leak[1]} scored {grid[leak]:.3f} on the dev "
                f"split vs {grid[(UP, DOWN)]:.3f} for the shipped window, but "
                "the gain is the CODING/UTR base-composition shift: a real "
                "start's downstream half is coding and a 5'UTR uATG's is not. "
                "That discriminates starts without saying anything about "
                "initiation, so the window stops at +4."),
            "efficiency_rate": (
                "no initiation RATE or protein-yield number is produced. "
                "This ranks ATG contexts against each other; it is not a "
                "trained predictor of how much protein an ATG makes."),
            "partial_windows": (
                f"an ATG with fewer than {UP} nt upstream is NOT scored. A "
                "shorter window sums fewer log-odds terms and is not "
                "comparable to the calibrated distribution."),
        },
    }

    if args.emit:
        pwm_txt = ",\n".join(
            "    {" + ", ".join(f'"{b}": {col[b]:+.6f}' for b in BASES) + "}"
            for col in PWM)
        quant_txt = "\n".join(
            f"    {q:+.4f},  # {p:>3d}th percentile"
            for p, q in zip(range(0, 101, 5), quant))
        Path(args.emit).write_text(MODULE_TEMPLATE.format(
            meta=json.dumps(meta, indent=4, sort_keys=True)
                     .replace(": true", ": True").replace(": false", ": False"),
            window=f"-{UP}..+{DOWN}", up=UP, consensus=consensus,
            at_real=at_real, at_decoy=at_dec, pwm=pwm_txt,
            n_quant=len(s_pos), quantiles=quant_txt))
        print(f"\n-> wrote {args.emit}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
