#!/usr/bin/env python3
"""Reproduce every measured number in `splicecraft_regulettuce`.

Fetches its own reference sets from public NCBI records, re-measures the
cis-element enrichments, the TATA positional peak and the 3'-end operating
curve, and prints them in the form the module's constants record. Nothing here
writes to the SpliceCraft data dir; nothing here imports the hub.

    python3 scripts/calibrate_regulettuce.py            # full run (~3 min)
    python3 scripts/calibrate_regulettuce.py --cache DIR  # reuse a download

Why this script exists rather than a prose citation: every threshold in that
module has to be defensible against a re-measurement, and the first attempt at
the promoter set was WRONG in a way only a QC check caught (92% of GenBank
single-gene `gene`/`mRNA` features are annotated at the start codon, not the
transcription start). The ATG check below is that QC, and it runs every time.
"""
from __future__ import annotations

import argparse
import bisect
import json
import random
import re
import statistics
import sys
import time
import urllib.parse
import urllib.request
from io import StringIO
from pathlib import Path

UA = {"User-Agent": "SpliceCraft-regulettuce-calibration/1.0"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UP, DOWN, W3P = 400, 100, 160
# Chromosome slices carrying TAIR/RefSeq annotation. Four 2 Mb windows is
# enough for ~3,000 transcripts with an annotated 5'UTR and keeps the fetch
# under a minute; widening it moves no number by more than sampling error.
SLICES = [("NC_003070.9", 1, 2_000_000), ("NC_003070.9", 10_000_000, 12_000_000),
          ("NC_003074.8", 1, 2_000_000), ("NC_003076.8", 5_000_000, 7_000_000)]
MRNA_TERMS = {
    "ath": ("Arabidopsis thaliana[ORGN] AND biomol_mrna[PROP] AND "
            "srcdb_refseq[PROP] AND 800:3000[SLEN]", 400),
    "osa": ("Oryza sativa[ORGN] AND biomol_mrna[PROP] AND "
            "srcdb_refseq[PROP] AND 800:3000[SLEN]", 200),
}


def _get(path: str, **params) -> str:
    url = EUTILS + path + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers=UA)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read(200_000_000).decode("utf-8", "replace")
        except Exception as exc:                      # noqa: BLE001 - retry
            if attempt == 2:
                raise
            print(f"  retry after {exc}", file=sys.stderr)
            time.sleep(3)
    raise RuntimeError("unreachable")


def _cached(cache: "Path | None", name: str, build):
    if cache is not None:
        f = cache / f"{name}.json"
        if f.exists():
            print(f"  (cached) {name}", file=sys.stderr)
            return json.loads(f.read_text())
    data = build()
    if cache is not None:
        cache.mkdir(parents=True, exist_ok=True)
        (cache / f"{name}.json").write_text(json.dumps(data))
    return data


def fetch_tss(cache) -> dict:
    """Verified plant TSS windows: mRNAs with an annotated 5'UTR only."""
    def build():
        from Bio import SeqIO
        out, seen_mrna = {}, 0
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        for acc, a, b in SLICES:
            print(f"  {acc}:{a}-{b}", file=sys.stderr)
            rec = next(SeqIO.parse(StringIO(_get(
                "efetch.fcgi", db="nucleotide", id=acc,
                rettype="gbwithparts", retmode="text",
                seq_start=a, seq_stop=b)), "genbank"))
            s = str(rec.seq).upper()
            cds: "dict[str, list[int]]" = {}
            for f in rec.features:
                if f.type != "CDS":
                    continue
                g = (f.qualifiers.get("gene")
                     or f.qualifiers.get("locus_tag") or [None])[0]
                if g:
                    cds.setdefault(g, []).append(
                        int(f.location.start) if f.location.strand == 1
                        else int(f.location.end) - 1)
            for f in rec.features:
                if f.type != "mRNA":
                    continue
                g = (f.qualifiers.get("gene")
                     or f.qualifiers.get("locus_tag") or [None])[0]
                if not g or g not in cds:
                    continue
                seen_mrna += 1
                st = f.location.strand
                if st == 1:
                    tss = int(f.location.start)
                    # A 5'UTR is what makes the 5' end a TRANSCRIPTION start
                    # rather than a start codon. Without this filter 92% of
                    # the set is ATGs.
                    if not any(c > tss for c in cds[g]):
                        continue
                    if tss < UP or tss + DOWN > len(s):
                        continue
                    win = s[tss - UP:tss + DOWN]
                elif st == -1:
                    tss = int(f.location.end) - 1
                    if not any(c < tss for c in cds[g]):
                        continue
                    if tss + UP >= len(s) or tss - DOWN + 1 < 0:
                        continue
                    win = "".join(comp.get(c, "N") for c in
                                  reversed(s[tss - DOWN + 1:tss + UP + 1]))
                else:
                    continue
                if len(win) == UP + DOWN and not (set(win) - set("ACGT")):
                    out[f"{acc}:{g}:{tss}"] = win
            time.sleep(0.4)
        print(f"  {seen_mrna} mRNA features seen -> {len(out)} usable windows",
              file=sys.stderr)
        return out
    return _cached(cache, "tss", build)


def fetch_genomic(cache) -> dict:
    def build():
        out = {}
        for acc, a, b in SLICES:
            txt = _get("efetch.fcgi", db="nucleotide", id=acc,
                       rettype="fasta", retmode="text",
                       seq_start=a, seq_stop=b)
            out[f"{acc}:{a}"] = "".join(
                ln.strip() for ln in txt.splitlines() if not ln.startswith(">"))
            time.sleep(0.4)
        return out
    return _cached(cache, "genomic", build)


def fetch_mrna(cache) -> dict:
    def build():
        out = {}
        for tag, (term, n) in MRNA_TERMS.items():
            ids = json.loads(_get("esearch.fcgi", db="nucleotide", term=term,
                                  retmax=n, retmode="json"
                                  ))["esearchresult"]["idlist"]
            for i in range(0, len(ids), 80):
                txt = _get("efetch.fcgi", db="nucleotide",
                           id=",".join(ids[i:i + 80]), rettype="fasta",
                           retmode="text")
                name, buf = None, []
                for ln in txt.splitlines():
                    if ln.startswith(">"):
                        if name:
                            out[name] = "".join(buf)
                        name, buf = f"{tag}:{ln[1:].split()[0]}", []
                    elif ln.strip():
                        buf.append(ln.strip())
                if name:
                    out[name] = "".join(buf)
                time.sleep(0.4)
        return out
    return _cached(cache, "mrna", build)


def _trim_tail(s: str) -> str:
    """Drop a trailing poly(A) tail. 3'UTRs are A-rich, so a few A's are real
    sequence and only a long run is a tail."""
    t = s.rstrip("Nn").upper()
    k = 0
    while k < len(t) and t[len(t) - 1 - k] == "A":
        k += 1
    return t[:len(t) - k] if k >= 8 else t


def _auc(a, b) -> float:
    sb = sorted(b)
    n = len(sb)
    total = 0.0
    for x in a:
        lo = bisect.bisect_left(sb, x)
        hi = bisect.bisect_right(sb, x)
        total += (lo + (hi - lo) / 2) / n
    return total / len(a)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", type=Path, default=None,
                    help="directory to cache downloads in")
    args = ap.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import splicecraft_regulettuce as rl

    random.seed(1)
    print("fetching reference sets from NCBI", file=sys.stderr)
    tss_all = list(fetch_tss(args.cache).values())
    genomic = fetch_genomic(args.cache)
    mrna = fetch_mrna(args.cache)

    # ---- QC: the check that caught the start-codon mistake ----------------
    atg = sum(1 for w in tss_all if w[UP:UP + 3] == "ATG")
    print("\n== QC ==")
    print(f"ATG exactly at the annotated start: {atg}/{len(tss_all)} = "
          f"{atg / len(tss_all):.3f}  (chance ~0.016)")
    if atg / len(tss_all) > 0.10:
        print("FAIL: these are start codons, not transcription starts. The "
              "5'UTR filter is not working -- every promoter number below "
              "would be measured on the wrong 400 bp.", file=sys.stderr)
        return 1
    tss = [w for w in tss_all if w[UP:UP + 3] != "ATG"]
    prom = [w[:UP] for w in tss]

    bg = []
    for s in genomic.values():
        s = s.upper()
        for _ in range(3000):
            i = random.randint(0, len(s) - UP)
            w = s[i:i + UP]
            if not (set(w) - set("ACGT")):
                bg.append(w)
    print(f"promoters {len(prom)}   null windows {len(bg)}")

    # ---- cis-element enrichment -------------------------------------------
    print("\n== cis-element catalogue (shipped vs re-measured) ==")
    print(f"{'element':24s} {'ship p':>7s} {'meas p':>7s} "
          f"{'ship null':>10s} {'meas null':>10s} {'drift':>6s}")
    worst = 0.0
    for row in rl._cis_element_catalogue():
        name = row["name"]
        p = sum(1 for w in prom
                if rl._scan_cis_elements(w, circular=False,
                                         names=[name])) / len(prom)
        b = sum(1 for w in bg
                if rl._scan_cis_elements(w, circular=False,
                                         names=[name])) / len(bg)
        drift = max(abs(p - row["promoter_rate"]), abs(b - row["null_rate"]))
        worst = max(worst, drift)
        print(f"{name:24s} {row['promoter_rate']:7.3f} {p:7.3f} "
              f"{row['null_rate']:10.3f} {b:10.3f} {drift:6.3f}")
    print(f"worst drift: {worst:.3f}")

    # ---- TATA --------------------------------------------------------------
    pat = re.compile("(?=TATA[AT]A)")
    gen = "".join(s.upper() for s in genomic.values()).replace("N", "")
    hits = len(pat.findall(gen))
    print("\n== TATA ==")
    print(f"genomic density: 1 per {len(gen) // hits} bp one strand "
          f"(module records {rl._TATA_GENOMIC_BP_PER_HIT})")
    best, bestn = None, -1
    for lo in range(-200, 40, 10):
        c = sum(1 for w in tss for m in pat.finditer(w)
                if lo <= m.start() - UP < lo + 10)
        if c > bestn:
            best, bestn = lo, c
    print(f"peak 10 bp bin: {best:+d}..{best + 9:+d} ({bestn} hits); "
          f"module records {rl._TATA_PEAK}")
    pin = sum(1 for w in tss if pat.search(w[UP - 40:UP - 20 + 6])) / len(tss)
    bin_ = sum(1 for w in bg if pat.search(w[UP - 40:UP - 20 + 6])) / len(bg)
    print(f"[-40,-20]: promoters {pin:.3f}  null {bin_:.3f}  "
          f"= {pin / bin_:.2f}x  (module records "
          f"{rl._TATA_PEAK_ENRICHMENT}x)")

    # ---- 3'-end operating curve -------------------------------------------
    pos, cod = [], []
    for s in mrna.values():
        t = _trim_tail(s)
        if len(t) < 2 * W3P + 200 or (set(t) - set("ACGT")):
            continue
        pos.append(t[-W3P:])
        cod.append(t[random.randint(100, len(t) - W3P - 200):][:W3P])
    gw = []
    for s in genomic.values():
        s = s.upper()
        for _ in range(1500):
            i = random.randint(0, len(s) - W3P)
            w = s[i:i + W3P]
            if not (set(w) - set("ACGT")):
                gw.append(w)

    def score(w):
        return max(rl._3prime_score(w)[0], rl._3prime_score(rl._rc(w))[0])

    sp, sc_, sg = [score(w) for w in pos], [score(w) for w in cod], \
        [score(w) for w in gw]
    print("\n== 3'-end-like regions ==")
    print(f"n: positives {len(pos)}  coding {len(cod)}  genomic {len(gw)}")
    print(f"AUC vs coding {_auc(sp, sc_):.3f} (module records "
          f"{rl._3P_AUC_VS_CODING})   vs genomic {_auc(sp, sg):.3f} "
          f"(module records {rl._3P_AUC_VS_GENOMIC})")
    print(f"means: pos {statistics.mean(sp):.3f}  coding "
          f"{statistics.mean(sc_):.3f}  genomic {statistics.mean(sg):.3f}")
    for thr in (0.60, 0.62, 0.64, 0.66, 0.68, 0.70):
        mark = "  <- default" if abs(thr - rl._3P_MIN_SCORE) < 1e-9 else (
            "  <- strong" if abs(thr - rl._3P_STRONG_SCORE) < 1e-9 else "")
        print(f"  >= {thr:.2f}: sensitivity "
              f"{sum(1 for x in sp if x >= thr) / len(sp):.3f}   coding FPR "
              f"{sum(1 for x in sc_ if x >= thr) / len(sc_):.4f}   genomic FPR "
              f"{sum(1 for x in sg if x >= thr) / len(sg):.4f}{mark}")

    # ---- the refusals ------------------------------------------------------
    print("\n== refused models (re-deriving the numbers that killed them) ==")
    print("Y-patch, >=8 nt C/T run in [-120,0]:")
    def yrun(w: str, k: int) -> bool:
        """Does the 120 bp upstream of this window's midpoint hold a C/T run
        of at least `k`? Both sets are 400/500 bp windows whose offset UP is
        the start site, so the same slice applies to each."""
        seg = w[UP - 120:UP]
        return max((len(m) for m in re.findall(r"[CT]+", seg)),
                   default=0) >= k
    for k in (8, 10):
        p = sum(1 for w in tss if yrun(w, k)) / len(tss)
        b = sum(1 for w in bg if yrun(w, k)) / len(bg)
        print(f"  >={k}: promoters {p:.3f}  null {b:.3f}  = {p / b:.2f}x")
    print("  (the module ships no Y-patch scanner: 1.05x at >=8 is no signal)")
    print("\nSee the module docstring for the poly(A) positional result "
          "(argmax within +-20 nt: 27.2% vs 29.1% chance) and the promoter "
          "composite (AUC 0.646); both are reproduced by the probes recorded "
          "there and neither is shipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
