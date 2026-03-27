#!/usr/bin/env python3
"""
splicecraft.py
==============
SpliceCraft — terminal circular plasmid map viewer.

Features:
  - Fetch any GenBank record by accession (pUC19 = L09137)
  - Load local .gb / .gbk files
  - Circular map with per-strand feature rings and arrowheads
  - Rotate origin freely with ← → keys or mouse scroll
  - Click map to select feature; click sidebar row to highlight on map
  - Feature detail panel
  - Plasmid library panel (left, SnapGene-style collection, persistent JSON)
  - DNA sequence viewer / editor (bottom, press e to edit, Ctrl+S to save)
  - CDS amino acid translation shown on feature selection

Run standalone:
    python3 splicecraft.py
    python3 splicecraft.py L09137          # fetch pUC19 on launch
    python3 splicecraft.py myplasmid.gb    # open local file
"""

import json
import math
import re
import sys
from io import StringIO
from pathlib import Path

# ── Dependency check ───────────────────────────────────────────────────────────

_REQUIRED = {
    "textual":  "textual",
    "Bio":      "biopython",
}

def _check_deps():
    missing = []
    for module, package in _REQUIRED.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(package)
    if missing:
        print("Missing dependencies — install with:")
        print(f"  pip install {' '.join(missing)}")
        sys.exit(1)

_check_deps()

from textual import on, work
from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical, ScrollableContainer
from textual.coordinate import Coordinate
from textual.events import Click, MouseDown, MouseMove, MouseUp, MouseScrollDown, MouseScrollUp
from textual.message import Message
from textual.reactive import reactive
from textual.screen import ModalScreen
from textual.widget import Widget
from textual.widgets import (
    Button, DataTable, Footer, Header, Input, Label,
    Select, Static, TabbedContent, TabPane,
)
from rich.text import Text

# ── Feature appearance ─────────────────────────────────────────────────────────

# Visually distinct per-feature palette (assigned by feature index, not type)
_FEATURE_PALETTE: list[str] = [
    "color(39)",   "color(118)",  "color(208)",  "color(213)",  "color(51)",
    "color(220)",  "color(196)",  "color(46)",   "color(201)",  "color(129)",
    "color(166)",  "color(33)",   "color(226)",  "color(160)",  "color(87)",
    "color(105)",  "color(154)",  "color(203)",  "color(81)",   "color(185)",
]

# ── Library persistence ────────────────────────────────────────────────────────

_LIBRARY_FILE = Path(__file__).parent / "plasmid_library.json"

def _load_library() -> list[dict]:
    if _LIBRARY_FILE.exists():
        try:
            return json.loads(_LIBRARY_FILE.read_text())
        except Exception:
            pass
    return []

def _save_library(entries: list[dict]) -> None:
    try:
        _LIBRARY_FILE.write_text(json.dumps(entries, indent=2))
    except Exception:
        pass

# ── Restriction sites ──────────────────────────────────────────────────────────

# Comprehensive NEB enzyme catalog.
# Each entry maps enzyme name → (recognition_sequence, fwd_cut, rev_cut).
#
# Cut position convention (0-based, relative to start of recognition sequence):
#   0          = before first base of recognition seq
#   len(seq)   = after last base of recognition seq
#   >len(seq)  = downstream of recognition seq (Type IIS)
#   negative   = upstream of recognition seq
#
# IUPAC ambiguity codes used in recognition sequences:
#   R=A/G  Y=C/T  W=A/T  S=G/C  M=A/C  K=G/T
#   B=not-A  D=not-C  H=not-G  V=not-T  N=any
#
# Examples:
#   EcoRI  G^AATTC / CTTAA^G  → fwd=1, rev=5
#   BamHI  G^GATCC / CCTAG^G  → fwd=1, rev=5
#   SmaI   CCC^GGG / GGG^CCC  → fwd=3, rev=3  (blunt)
#   BsaI   GGTCTC(1/5)        → recognition=6bp, fwd=6+1=7, rev=6+5=11
#
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
    "BstEII":    ("GGTNACC",      1,  5),  # G^GTNACC  / CCANT^GG    5' overhang
    "BstXI":     ("CCANNNNNTGG",  8,  4),  # CCANN4^TGG/ CCA^N5TGG   3' overhang
    "BstYI":     ("RGATCY",       1,  5),  # R^GATCY   / YCTAG^R     5' overhang
    "CpoI":      ("CGGWCCG",      2,  5),  # CG^GWCCG  / CGCC^WGG    5' overhang
    "DraI":      ("TTTAAA",       3,  3),  # TTT^AAA   / TTT^AAA     blunt
    "DraIII":    ("CACNNNGTG",    6,  3),  # CACNNN^GTG/ GTG^NNNGTG  3' overhang
    "EagI":      ("CGGCCG",       1,  5),  # C^GGCCG   / GCCGG^C     5' overhang (NotI subset)
    "Eco47III":  ("AGCGCT",       3,  3),  # AGC^GCT                 blunt
    "Eco53kI":   ("GAGCTC",       5,  1),  # GAGCT^C                 3' overhang (SacI isoschizomer)
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
    "SfiI":      ("GGCCNNNNNGGCC",8,  4),  # GGCCN4^NGGCC            3' overhang (13-bp)
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
    "AccI":      ("GTYRAC",       2,  4),  # GT^YRAC                 3' overhang
    "AclI":      ("AACGTT",       2,  4),  # AA^CGTT                 3' overhang
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
    "BseYI":     ("CCCAGC",       5,  1),  # CCCAG^C   / C^CCAGC     3' overhang
    "BsiEI":     ("CGRYCG",       4,  2),  # CGRY^CG                 3' overhang
    "BsiHKAI":   ("GWGCWC",       5,  1),  # GWGCW^C                 3' overhang
    "BsrFI":     ("RCCGGY",       1,  5),  # R^CCGGY                 5' overhang
    "Bsp1286I":  ("GDGCHC",       5,  1),  # GDGCH^C                 3' overhang
    "BspHI":     ("TCATGA",       1,  5),  # T^CATGA                 NcoI-compatible ends
    "BsrI":      ("ACTGG",        1,  1),  # ACT^GG                  degenerate Type IIS-like
    "BstAPI":    ("GCANNNNNTGC",  7,  5),  # GCANNNN^NTGC            3' overhang
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
    "HphI":      ("GGTGA",        8, 12),  # GGTGA(8/7) downstream   Type IIS
    "KasI":      ("GGCGCC",       1,  5),  # G^GCGCC                 5' overhang (NarI isoschizomer)
    "MboI":      ("GATC",         0,  4),  # ^GATC                   DpnII isoschizomer
    "MboII":     ("GAAGA",        8, 12),  # GAAGA(8/7) downstream   Type IIS
    "MlyI":      ("GAGTC",       10, 10),  # GAGTC(5/5) downstream   blunt, Type IIS
    "MmeI":      ("TCCRAC",      20, 18),  # TCCRAC(20/18)           Type IIS far-cutter
    "MspA1I":    ("CMGCKG",       3,  3),  # CMG^CKG                 blunt
    "NgoMIV":    ("GCCGGC",       1,  5),  # G^CCGGC                 5' overhang (EagI-compatible)
    "NmeAIII":   ("GCCGAG",      21, 19),  # GCCGAG(15/13)           Type IIS far-cutter
    "NspI":      ("RCATGY",       5,  1),  # RCATG^Y                 3' overhang
    "PflMI":     ("CCANNNNNTGG",  7,  4),  # CCANN4^NTGG             3' overhang (BstXI isoschizomer)
    "PspOMI":    ("GGGCCC",       1,  5),  # G^GGCCC                 ApaI isoschizomer (5' overhang)
    "Sau3AI":    ("GATC",         0,  4),  # ^GATC                   BamHI-compatible ends (4-cutter)
    "SbfI":      ("CCTGCAGG",     6,  2),  # CCTGCA^GG               already PstI-compatible (see above)
    "SfcI":      ("CTRYAG",       1,  5),  # C^TRYAG                 5' overhang
    "SspI":      ("AATATT",       3,  3),  # AAT^ATT                 blunt
    "TaqI":      ("TCGA",         1,  3),  # T^CGA     / CGT^A       5' overhang (heat-stable)
    "Van91I":    ("CCANNNNNTGG",  7,  4),  # PflMI isoschizomer
    "ZraI":      ("GACGTC",       3,  3),  # GAC^GTC                 blunt (AatII-related)

    # ── Type IIS — cut outside recognition sequence ────────────────────────────
    # fwd/rev positions are still offsets from start of recognition seq.
    # For an n-bp recognition sequence cutting d1/d2 downstream:
    #   fwd = n + d1,  rev = n + d2
    "BaeI":      ("ACNNNNGTAYYC",-10, -7), # cuts 10/15 upstream (negative = upstream)
    "BbsI":      ("GAAGAC",       8, 12),  # GAAGAC(2/6)  BpiI isoschizomer
    "BcoDI":     ("GTCTC",        6, 10),  # GTCTC(1/5)   BsaI 5-bp variant
    "BceAI":     ("ACGGC",        9, 11),  # ACGGC(4/6)
    "BciVI":     ("GTATCC",      12,  6),  # GTATCC(6/5)  asymmetric
    "BfuAI":     ("ACCTGC",      10, 14),  # ACCTGC(4/8)  BspMI isoschizomer
    "BmrI":      ("ACTGGN",       6,  5),  # cuts 1 past end of recog
    "BpiI":      ("GAAGAC",       8, 12),  # BbsI isoschizomer
    "BsaI":      ("GGTCTC",       7, 11),  # GGTCTC(1/5)  Golden Gate workhorse
    "BsaXI":     ("ACNNNNNCTCC",  3, 12),  # cuts upstream and downstream (unusual)
    "BsbI":      ("CAACAC",      17, 15),  # far-cutter
    "BseJI":     ("GAAGAC",       8, 12),  # BbsI isoschizomer
    "BseLI":     ("CCNNNNNNNGG",  7,  4),  # 3' overhang
    "BseMII":    ("CTCAG",       10,  8),  # CTCAG(5/3)
    "BseRI":     ("GAGGAG",      10, 14),  # GAGGAG(10/8)
    "BsgI":      ("GTGCAG",      22, 20),  # GTGCAG(16/14) far-cutter
    "BslI":      ("CCNNNNNNNGG",  7,  4),  # 3' overhang (BseLI variant)
    "BsmAI":     ("GTCTC",        6, 10),  # GTCTC(1/5)   BsaI isoschizomer (5-bp)
    "BsmBI":     ("CGTCTC",       7, 11),  # CGTCTC(1/5)  Esp3I isoschizomer
    "BsmFI":     ("GGGAC",       15, 19),  # GGGAC(10/14)
    "BsmI":      ("GAATGC",       7,  1),  # GAATGC(1/-1) asymmetric Type IIS
    "BspLU11III":("ACATGT",       1,  5),  # PciI isoschizomer (not strictly IIS)
    "BspMI":     ("ACCTGC",      10, 14),  # ACCTGC(4/8)
    "BspQI":     ("GCTCTTC",      8, 11),  # SapI isoschizomer
    "BspTNI":    ("GGTCTC",       7, 11),  # BsaI isoschizomer
    "BsrBI":     ("CCGCTC",       3,  3),  # cuts within recog (special case)
    "BsrDI":     ("GCAATG",       8,  6),  # GCAATG(2/0)
    "BssSI":     ("CACGAG",      -5, -1),  # cuts upstream of recognition
    "BtgZI":     ("GCGATG",      16, 20),  # GCGATG(10/14)
    "BtsCI":     ("GGATG",        5,  3),  # GGATG(2/0)
    "BtsI":      ("GCAGTG",       8,  6),  # GCAGTG(2/0)
    "BtsImutI":  ("CAGTG",        5,  3),  # BtsCI variant
    "EarI":      ("CTCTTC",      10,  6),  # CTCTTC(4/1)  SapI-related
    "Esp3I":     ("CGTCTC",       7, 11),  # BsmBI isoschizomer
    "PaqCI":     ("CACCTGC",     11, 15),  # CACCTGC(4/8)
    "SapI":      ("GCTCTTC",      8, 11),  # GCTCTTC(1/4)
    "BsmBI-v2":  ("CGTCTC",       7, 11),  # v2/HF variant

    # ── High-Fidelity (HF) and v2 variants — same recognition/cut as canonical ─
    "AgeI-HF":   ("ACCGGT",       1,  5),
    "BamHI-HF":  ("GGATCC",       1,  5),
    "BclI-HF":   ("TGATCA",       1,  5),
    "BmtI":      ("GCTAGC",       1,  5),  # NheI isoschizomer / HF variant
    "BsiWI-HF":  ("CGTACG",       1,  5),
    "BsrFI-v2":  ("RCCGGY",       1,  5),
    "BsrGI-HF":  ("TGTACA",       1,  5),
    "BssSI-v2":  ("CACGAG",      -5, -1),
    "BstEII-HF": ("GGTNACC",      1,  5),
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

# Flat recognition-sequence-only dict derived from _NEB_ENZYMES (used by
# _scan_restriction_sites and _RESTR_COLOR below).
_RESTRICTION_SITES: dict[str, str] = {
    name: seq for name, (seq, _, _) in _NEB_ENZYMES.items()
}


_RESTR_PALETTE: list[str] = [
    "color(220)", "color(208)", "color(196)", "color(160)",
    "color(105)", "color(129)", "color(57)",  "color(21)",
    "color(33)",  "color(39)",  "color(51)",  "color(87)",
    "color(118)", "color(154)", "color(190)", "color(226)",
    "color(185)", "color(180)",
]
_RESTR_COLOR: dict[str, str] = {
    name: _RESTR_PALETTE[i % len(_RESTR_PALETTE)]
    for i, name in enumerate(_RESTRICTION_SITES)
}


_IUPAC_RE: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "W": "[AT]", "S": "[CG]",
    "M": "[AC]", "K": "[GT]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def _iupac_pattern(site: str) -> "re.Pattern[str]":
    return re.compile("".join(_IUPAC_RE.get(c, c) for c in site.upper()))


def _rc(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _scan_restriction_sites(
    seq: str,
    min_recognition_len: int = 6,
    unique_only: bool = True,
) -> list[dict]:
    """Scan both strands; return resite + recut dicts for every hit.

    resite — the recognition sequence span (colored bar)
    recut  — the cut position (single-bp marker: ↓ above or ↑ below DNA)

    min_recognition_len — skip enzymes whose recognition site is shorter than this
                          (default 6 to reduce noise from 4-cutters)
    unique_only         — if True, only include enzymes that cut exactly once
                          (forward + reverse strand combined; default True)
    """
    seq_u = seq.upper()
    n = len(seq_u)
    # Per-enzyme results collected first so we can filter to unique cutters
    by_enzyme: dict[str, list[dict]] = {}
    seen: set[tuple[str, int, int]] = set()   # deduplicate palindromes

    for name, (site, fwd_cut, rev_cut) in _NEB_ENZYMES.items():
        if len(site) < min_recognition_len:
            continue
        color    = _RESTR_COLOR[name]
        pat      = _iupac_pattern(site)
        site_len = len(site)
        hits: list[dict] = []

        # Forward strand scan
        for m in pat.finditer(seq_u):
            p = m.start()
            key = (name, p, 1)
            if key in seen:
                continue
            seen.add(key)
            hits.append({
                "type":    "resite",
                "start":   p,
                "end":     p + site_len,
                "strand":  1,
                "color":   color,
                "label":   name,
                "cut_col": fwd_cut if 0 < fwd_cut < site_len else None,
            })
            cut_bp = min(p + fwd_cut, n - 1)
            hits.append({
                "type":   "recut",
                "start":  cut_bp,
                "end":    cut_bp + 1,
                "strand": 1,
                "color":  color,
                "label":  name,
            })

        # Reverse strand scan
        rc_site = _rc(site)
        rc_pat  = _iupac_pattern(rc_site) if rc_site != site.upper() else pat
        for m in rc_pat.finditer(seq_u):
            p = m.start()
            orig_start = n - p - site_len
            key = (name, orig_start, -1)
            if key in seen:
                continue
            seen.add(key)
            hits.append({
                "type":    "resite",
                "start":   orig_start,
                "end":     orig_start + site_len,
                "strand":  -1,
                "color":   color,
                "label":   name,
                "cut_col": rev_cut if 0 < rev_cut < site_len else None,
            })
            cut_bp = min(orig_start + rev_cut, n - 1)
            hits.append({
                "type":   "recut",
                "start":  cut_bp,
                "end":    cut_bp + 1,
                "strand": -1,
                "color":  color,
                "label":  name,
            })

        if hits:
            by_enzyme[name] = hits

    feats: list[dict] = []
    for name, hits in by_enzyme.items():
        # Count recognition-sequence hits (resite only) across both strands
        if unique_only:
            n_sites = sum(1 for h in hits if h["type"] == "resite")
            if n_sites != 1:
                continue
        feats.extend(hits)
    return feats


def _assign_chunk_features(
    chunk_feats: list[dict], chunk_start: int, chunk_end: int
) -> tuple[list[list[dict]], list[list[dict]]]:
    """
    Forward-strand features always go above DNA; reverse-strand always below.
    Within each group, overlapping features are stacked into greedy lanes
    (each lane is a list of non-overlapping features).  Capped at 3 lanes/side.
    """
    def _greedy_lanes(feats: list[dict]) -> list[list[dict]]:
        sorted_f = sorted(feats, key=lambda f: max(f["start"], chunk_start))
        lanes: list[list[dict]] = []
        lane_ends: list[int]    = []
        for f in sorted_f:
            bar_s = max(f["start"], chunk_start)
            bar_e = min(f["end"],   chunk_end)
            if bar_e - bar_s <= 0:
                continue
            placed = False
            for i, end in enumerate(lane_ends):
                if bar_s >= end:
                    lanes[i].append(f)
                    lane_ends[i] = bar_e
                    placed = True
                    break
            if not placed:
                lanes.append([f])
                lane_ends.append(bar_e)
        return lanes[:3]

    fwd = [f for f in chunk_feats if f["strand"] >= 0]
    rev = [f for f in chunk_feats if f["strand"] <  0]
    return _greedy_lanes(fwd), _greedy_lanes(rev)


def _render_feature_row_pair(
    result: "Text",
    feats: list[dict],
    chunk_start: int,
    chunk_end: int,
    prefix_w: int,
    is_below_dna: bool,
    show_connectors: bool,
) -> None:
    """
    Append one label row + optional connector row + one braille-bar row to result.
    For above-DNA: label / [connector] / bar.
    For below-DNA: bar / [connector] / label.
    Multiple non-overlapping features share the same pair of rows horizontally.
    """
    content_w = chunk_end - chunk_start
    label_arr: list[tuple[str, str]] = [(" ", "")] * content_w
    bar_arr:   list[tuple[str, str]] = [(" ", "")] * content_w
    conn_arr:  list[tuple[str, str]] = [(" ", "")] * content_w

    for f in feats:
        bar_s = max(f["start"], chunk_start) - chunk_start
        bar_e = min(f["end"],   chunk_end)   - chunk_start
        bar_len = bar_e - bar_s
        if bar_len <= 0:
            continue
        starts_here = f["start"] >= chunk_start
        ends_here   = f["end"]   <= chunk_end
        strand      = f["strand"]
        color       = f["color"]
        label       = f.get("label", f.get("type", ""))

        feat_type = f.get("type", "")

        if feat_type == "resite":
            # ── Parenthesis-style RE site: ( EnzymeName ) ─────────────────
            # Bar row:   ( EnzymeName )  — bold white name, colored parens
            # Label row: ↓ or ↑ at the intra-site cut column (if cut is
            #            within the recognition site; Type IIS cuts elsewhere)
            cut_col = f.get("cut_col")   # 0-based offset from f["start"], or None

            # Opening / closing parens
            if starts_here and bar_len >= 1:
                bar_arr[bar_s] = ("(", color)
            if ends_here and bar_len >= 1:
                bar_arr[bar_s + bar_len - 1] = (")", color)

            # Enzyme name: bold white, centered in the interior
            interior_start = (1 if starts_here else 0)
            interior_end   = (bar_len - 1 if ends_here else bar_len)
            interior_len   = interior_end - interior_start
            if interior_len > 0:
                name_str  = label[:interior_len]
                name_pad  = interior_len - len(name_str)
                name_lpad = name_pad // 2
                for j, ch in enumerate(name_str):
                    pos = bar_s + interior_start + name_lpad + j
                    if 0 <= pos < content_w:
                        bar_arr[pos] = (ch, "bold white")

            # Cut marker in the label row so it doesn't obscure the name
            if cut_col is not None:
                visible_offset = cut_col - max(0, chunk_start - f["start"])
                cut_pos = bar_s + visible_offset
                if 0 <= cut_pos < content_w:
                    cut_ch = "↑" if is_below_dna else "↓"
                    label_arr[cut_pos] = (cut_ch, "bold " + color)

            # Connector tick at midpoint
            mid = bar_s + bar_len // 2
            if 0 <= mid < content_w:
                conn_arr[mid] = ("┊", color)
            continue   # skip the regular label/bar/conn logic below

        elif feat_type == "recut":
            continue   # cut position is rendered inside the resite bar; skip here

        # ── Standard feature (non-RE) ──────────────────────────────────────

        # Label (centered in feature span)
        lbl = label[:bar_len]
        pad = bar_len - len(lbl)
        pl  = pad // 2
        lbl_str = " " * pl + lbl + " " * (pad - pl)
        for i, ch in enumerate(lbl_str):
            if 0 <= bar_s + i < content_w:
                label_arr[bar_s + i] = (ch, color)

        # Braille bar
        if bar_len == 1:
            bar_str = "▲" if is_below_dna else "▼"
        elif strand >= 0:
            bar_str = "⣿" * (bar_len - (1 if ends_here   else 0)) + ("▶" if ends_here   else "")
        else:
            bar_str = ("◀" if starts_here else "") + "⣿" * (bar_len - (1 if starts_here else 0))
        for i, ch in enumerate(bar_str):
            if 0 <= bar_s + i < content_w:
                bar_arr[bar_s + i] = (ch, color)

        # Connector tick at midpoint of feature span
        mid = bar_s + bar_len // 2
        if 0 <= mid < content_w:
            conn_arr[mid] = ("┊", color)

    def _write_arr(arr: list[tuple[str, str]]) -> None:
        result.append(" " * prefix_w)
        run: list[str] = []
        sty = ""
        for ch, s in arr:
            if s == sty:
                run.append(ch)
            else:
                if run:
                    result.append("".join(run), style=sty)
                run = [ch]
                sty = s
        if run:
            result.append("".join(run), style=sty)
        result.append("\n")

    if not is_below_dna:
        _write_arr(label_arr)
        if show_connectors:
            _write_arr(conn_arr)
        _write_arr(bar_arr)
    else:
        _write_arr(bar_arr)
        if show_connectors:
            _write_arr(conn_arr)
        _write_arr(label_arr)


def _build_seq_text(seq: str, feats: list[dict], line_width: int = 60,
                    sel_range: "tuple[int,int] | None" = None,
                    user_sel:  "tuple[int,int] | None" = None,
                    cursor_pos: int = -1,
                    show_connectors: bool = False) -> Text:
    """Rich Text of the sequence with per-position feature coloring.

    sel_range  — feature highlight: bold + underline on feature bases
    user_sel   — shift-click selection: subtle background, used by edit dialog
    cursor_pos — click cursor: │ inserted before cursor_pos
    Annotation bars appear below each DNA line (one bar per overlapping
    non-site feature, capped at 4, largest first).
    Each feature renders as:  label row  /  [connector row]  /  braille bar row.
    """
    n = len(seq)
    styles = ["color(252)"] * n
    for f in reversed(feats):          # reversed so first feature wins
        col = f["color"]
        for i in range(f["start"], min(f["end"], n)):
            styles[i] = col

    sel_s  = sel_range[0] if sel_range else -1
    sel_e  = sel_range[1] if sel_range else -1
    usr_s  = user_sel[0]  if user_sel  else -1
    usr_e  = user_sel[1]  if user_sel  else -1

    seq_upper = seq.upper()
    _COMP     = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")   # base complement
    result    = Text(no_wrap=False)

    # Annotation-bar features: exclude old "site" and "recut" (cut pos is
    # embedded inside the resite bar; recut only used by the map overlays).
    annot_feats = sorted(
        [f for f in feats if f.get("type") not in ("site", "recut")],
        key=lambda f: -(f["end"] - f["start"]),
    )

    for chunk_start in range(0, n, line_width):
        chunk_end = min(chunk_start + line_width, n)

        # ── Assign features to above / below lanes ──
        chunk_feats = [
            f for f in annot_feats
            if f["start"] < chunk_end and f["end"] > chunk_start
        ]
        above_lanes, below_lanes = _assign_chunk_features(chunk_feats, chunk_start, chunk_end)

        # ── Feature rows ABOVE DNA (forward strand) ──
        for lane in above_lanes:
            _render_feature_row_pair(result, lane, chunk_start, chunk_end,
                                     10, False, show_connectors)

        # ── Double-stranded DNA block ─────────────────────────────────────
        # ── Rows 1+2: double-stranded DNA block ──────────────────────────
        # Forward strand:  "   1  5'─ATGC…─3'"
        # RC strand:       "      3'─TACG…─5'"
        #
        # ALIGNMENT GUARANTEE: no extra characters are ever inserted into
        # either strand — both always emit exactly (chunk_end-chunk_start)
        # base characters between their fixed-width prefixes and suffixes.
        # The cursor is shown as a reverse-video highlight on the base IN
        # PLACE rather than as an inserted │ glyph, so column counts match.

        def _strand_chars(bases: "list[str]") -> None:
            """Append base chars with RLE styling into result."""
            run: list[str] = []
            sty = ""
            for ch, s in bases:
                if s == sty:
                    run.append(ch)
                else:
                    if run:
                        result.append("".join(run), style=sty)
                    run, sty = [ch], s
            if run:
                result.append("".join(run), style=sty)

        fwd_bases: list[tuple[str, str]] = []
        rc_bases:  list[tuple[str, str]] = []
        for i in range(chunk_start, chunk_end):
            base   = styles[i]
            in_usr = (usr_s <= i < usr_e)
            in_sel = (sel_s <= i < sel_e)
            is_cur = (cursor_pos == i)
            if is_cur:
                sty = "reverse bold white"
            elif in_usr:
                sty = base + " on color(237)"
            elif in_sel:
                sty = "bold underline " + base
            else:
                sty = base
            fwd_bases.append((seq_upper[i],               sty))
            rc_bases.append( (seq_upper[i].translate(_COMP), sty))

        # Forward strand
        result.append(f"{chunk_start + 1:>8}  5'─", style="color(245)")
        _strand_chars(fwd_bases)
        result.append("─3'\n", style="color(245)")

        # Reverse-complement strand (aligned column-for-column)
        result.append(" " * 10 + "3'─", style="color(245)")
        _strand_chars(rc_bases)
        result.append("─5'\n", style="color(245)")

        # ── Feature rows BELOW DNA (reverse strand) ──
        for lane in below_lanes:
            _render_feature_row_pair(result, lane, chunk_start, chunk_end,
                                     10, True, show_connectors)

    return result


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

def _copy_to_clipboard(text: str) -> bool:
    """Copy text to system clipboard via xclip / xsel / wl-copy. Returns True on success."""
    import subprocess
    for cmd in (
        ["xclip", "-selection", "clipboard"],
        ["xsel", "--clipboard", "--input"],
        ["wl-copy"],
    ):
        try:
            subprocess.run(cmd, input=text.encode(), check=True, timeout=3,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return True
        except Exception:
            continue
    try:
        import pyperclip
        pyperclip.copy(text)
        return True
    except Exception:
        return False


def _translate_cds(full_seq: str, start: int, end: int, strand: int) -> str:
    """Translate a CDS region to single-letter AA string (stop codon → *)."""
    sub = full_seq[start:end].upper()
    if strand == -1:
        sub = sub.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    aa = [_CODON_TABLE.get(sub[i:i+3], "?") for i in range(0, len(sub) - 2, 3)]
    result = "".join(aa)
    if result and not result.endswith("*"):
        result += "*"
    return result


# ── Char aspect detection ──────────────────────────────────────────────────────

def _detect_char_aspect() -> float:
    """
    Return char_height / char_width for the current terminal by reading
    pixel dimensions via TIOCGWINSZ.  Falls back to 2.0 if unavailable.
    """
    try:
        import fcntl, os, struct, termios
        fds: list[tuple[int, bool]] = []
        try:
            fds.append((os.open("/dev/tty", os.O_RDWR | os.O_NOCTTY), True))
        except OSError:
            pass
        fds += [(0, False), (1, False), (2, False)]
        try:
            for fd, _ in fds:
                try:
                    buf = struct.pack("HHHH", 0, 0, 0, 0)
                    res = fcntl.ioctl(fd, termios.TIOCGWINSZ, buf)
                    rows, cols, px_w, px_h = struct.unpack("HHHH", res)
                    if px_w > 0 and px_h > 0 and rows > 0 and cols > 0:
                        ratio = (px_h / rows) / (px_w / cols)
                        if 0.8 <= ratio <= 5.0:
                            return round(ratio, 3)
                except OSError:
                    continue
        finally:
            for fd, should_close in fds:
                if should_close:
                    try:
                        os.close(fd)
                    except OSError:
                        pass
    except Exception:
        pass
    return 2.0


def _format_bp(bp: int) -> str:
    """Human-readable bp label with up to 2 decimal places in 'k' units."""
    if bp < 1000:
        return str(bp)
    if bp % 1000 == 0:
        return f"{bp // 1000}k"
    if bp % 100 == 0:
        return f"{bp / 1000:.1f}k"
    return f"{bp / 1000:.2f}k"


def _feat_label(feat) -> str:
    for q in ("label", "gene", "product", "standard_name", "note", "bound_moiety"):
        if q in feat.qualifiers:
            v = feat.qualifiers[q]
            s = v[0] if isinstance(v, list) else v
            return s[:28]
    return feat.type

def _nice_tick(total: int) -> int:
    """A tick interval that gives ~6-10 ticks for this plasmid size."""
    for t in [50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 25000, 50000]:
        if 4 <= total // t <= 14:
            return t
    return max(1, total // 8)

# ── GenBank I/O ────────────────────────────────────────────────────────────────

def fetch_genbank(accession: str, email: str = "splicecraft@local"):
    """Fetch a GenBank record by accession from NCBI Entrez. Returns SeqRecord."""
    from Bio import Entrez, SeqIO
    Entrez.email = email
    with Entrez.efetch(
        db="nucleotide", id=accession, rettype="gb", retmode="text"
    ) as handle:
        record = SeqIO.read(handle, "genbank")
    return record

def load_genbank(path: str):
    """Load a GenBank (.gb/.gbk) file. Returns SeqRecord."""
    from Bio import SeqIO
    return SeqIO.read(path, "genbank")

def _record_to_gb_text(record) -> str:
    """Serialize a SeqRecord to GenBank format text."""
    from Bio import SeqIO
    buf = StringIO()
    SeqIO.write(record, buf, "genbank")
    return buf.getvalue()

def _gb_text_to_record(text: str):
    """Parse GenBank format text back to a SeqRecord."""
    from Bio import SeqIO
    return SeqIO.read(StringIO(text), "genbank")

# ── Core drawing ───────────────────────────────────────────────────────────────

class _Canvas:
    """A mutable 2-D character grid that renders to a Rich Text object."""

    def __init__(self, width: int, height: int):
        self.w = width
        self.h = height
        self._chars:  list[list[str]] = [[" "] * width for _ in range(height)]
        self._styles: list[list[str]] = [[""]  * width for _ in range(height)]

    def put(self, col: int, row: int, ch: str, style: str = ""):
        if 0 <= col < self.w and 0 <= row < self.h:
            self._chars[row][col]  = ch
            self._styles[row][col] = style

    def put_text(self, col: int, row: int, text: str, style: str = ""):
        for j, ch in enumerate(text):
            self.put(col + j, row, ch, style)

    def to_rich_text(self) -> Text:
        result = Text(no_wrap=True, overflow="crop")
        for r in range(self.h):
            for c in range(self.w):
                ch = self._chars[r][c]
                st = self._styles[r][c]
                result.append(ch, style=st) if st else result.append(ch)
            if r < self.h - 1:
                result.append("\n")
        return result


class _BrailleCanvas:
    """
    Sub-character resolution canvas using Unicode braille (U+2800–U+28FF).

    Each terminal cell (col, row) encodes a 2-wide × 4-tall dot grid —
    8 pixels per character cell.  Braille dot layout:

        px%2=0  px%2=1
        dot1    dot4    ← py%4=0   (bits 0, 3)
        dot2    dot5    ← py%4=1   (bits 1, 4)
        dot3    dot6    ← py%4=2   (bits 2, 5)
        dot7    dot8    ← py%4=3   (bits 6, 7)

    Codepoint = 0x2800 + bitmask of active dots.
    Colors: higher-priority write wins per cell.
    """

    _DOT_BITS: list[list[int]] = [
        [0, 3],
        [1, 4],
        [2, 5],
        [6, 7],
    ]

    def __init__(self, cols: int, rows: int):
        self.cols = cols
        self.rows = rows
        self._bits:   list[list[int]] = [[0]  * cols for _ in range(rows)]
        self._colors: list[list[str]] = [[" "] * cols for _ in range(rows)]
        self._prio:   list[list[int]] = [[0]  * cols for _ in range(rows)]

    def set_pixel(self, px: int, py: int,
                  color: str = "", priority: int = 1) -> None:
        col, row = px // 2, py // 4
        if not (0 <= col < self.cols and 0 <= row < self.rows):
            return
        self._bits[row][col] |= 1 << self._DOT_BITS[py % 4][px % 2]
        if color and priority >= self._prio[row][col]:
            self._colors[row][col] = color
            self._prio[row][col]   = priority

    def combine(self, text_canvas: "_Canvas") -> Text:
        """
        Return a Rich Text object.
        Non-space cells from *text_canvas* are drawn on top;
        braille pixels fill the rest.
        Consecutive blank cells are batched into a single append call.
        """
        result = Text(no_wrap=True, overflow="crop")
        rows   = min(self.rows, text_canvas.h)
        cols   = min(self.cols, text_canvas.w)
        tc_chars  = text_canvas._chars
        tc_styles = text_canvas._styles
        bc_bits   = self._bits
        bc_colors = self._colors
        for row in range(rows):
            blank_run = 0
            tc_row  = tc_chars[row]
            tcs_row = tc_styles[row]
            bc_bits_row   = bc_bits[row]
            bc_colors_row = bc_colors[row]
            for col in range(cols):
                tc_ch = tc_row[col]
                if tc_ch == " " and not bc_bits_row[col]:
                    blank_run += 1
                else:
                    if blank_run:
                        result.append(" " * blank_run)
                        blank_run = 0
                    if tc_ch != " ":
                        st = tcs_row[col]
                        result.append(tc_ch, style=st) if st else result.append(tc_ch)
                    else:
                        ch = chr(0x2800 + bc_bits_row[col])
                        c  = bc_colors_row[col]
                        result.append(ch, style=c) if c != " " else result.append(ch)
            if blank_run:
                result.append(" " * blank_run)
            if row < rows - 1:
                result.append("\n")
        return result


def _arrow_char(tangent_angle: float) -> str:
    """Pick a directional arrow char for an arrowhead."""
    t = tangent_angle % (2 * math.pi)
    sector = int((t + math.pi / 4) / (math.pi / 2)) % 4
    return ["▶", "▼", "◀", "▲"][sector]


# ── PlasmidMap widget ──────────────────────────────────────────────────────────

class PlasmidMap(Widget):
    """
    Circular plasmid map.

    Keyboard: ← → to rotate when focused (click map first); [ ] always rotate.
    Mouse:    scroll to rotate; click to select a feature.
    """

    DEFAULT_CSS = """
    PlasmidMap {
        width: 1fr;
        height: 1fr;
        background: $background;
    }
    PlasmidMap:focus { border: solid $accent; }
    """

    can_focus = True

    BINDINGS = [
        Binding("left",        "rotate_cw",        "Rotate ←",      show=True),
        Binding("right",       "rotate_ccw",        "Rotate →",      show=True),
        Binding("shift+left",  "rotate_cw_lg",     "Rotate ←←",     show=False),
        Binding("shift+right", "rotate_ccw_lg",    "Rotate →→",     show=False),
        Binding("home",        "reset_origin",     "Reset",         show=False),
        Binding("comma",       "aspect_dec",       "Circle wider",   show=False),
        Binding("full_stop",   "aspect_inc",       "Circle taller",  show=False),
        Binding("v",           "toggle_map_view",  "Toggle view",    show=False),
    ]

    origin_bp:    reactive[int]   = reactive(0)
    selected_idx: reactive[int]   = reactive(-1)
    _aspect:      reactive[float] = reactive(2.0)
    _map_mode:    reactive[str]   = reactive("circular")

    # ── Messages ───────────────────────────────────────────────────────────────

    class FeatureSelected(Message):
        def __init__(self, idx: int, feat_dict: dict | None):
            self.idx       = idx
            self.feat_dict = feat_dict
            super().__init__()

    # ── Lifecycle ──────────────────────────────────────────────────────────────

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.record  = None
        self._feats:          list[dict] = []
        self._restr_feats:    list[dict] = []   # restriction site overlay
        self._total:          int  = 0
        self._show_connectors: bool = False

    def on_mount(self) -> None:
        detected = _detect_char_aspect()
        if detected != self._aspect:
            self._aspect = detected

    def load_record(self, record) -> None:
        self.record       = record
        self._total       = len(record.seq)
        self.origin_bp    = 0
        self.selected_idx = -1
        self._feats       = self._parse(record)
        self._restr_feats = []
        self.refresh()

    def _parse(self, record) -> list[dict]:
        feats = []
        for feat in record.features:
            if feat.type in ("source",):
                continue
            start  = int(feat.location.start)
            end    = int(feat.location.end)
            strand = getattr(feat.location, "strand", 1) or 1
            idx    = len(feats)
            feats.append({
                "type":   feat.type,
                "start":  start,
                "end":    end,
                "strand": strand,
                "color":  _FEATURE_PALETTE[idx % len(_FEATURE_PALETTE)],
                "label":  _feat_label(feat),
            })
        return feats

    # ── Actions ────────────────────────────────────────────────────────────────

    def _step(self, coarse: bool = False) -> int:
        if not self._total:
            return 1
        return max(1, self._total // (10 if coarse else 200))

    def action_rotate_cw(self):
        if self._total:
            self.origin_bp = (self.origin_bp - self._step()) % self._total

    def action_rotate_ccw(self):
        if self._total:
            self.origin_bp = (self.origin_bp + self._step()) % self._total

    def action_rotate_cw_lg(self):
        if self._total:
            self.origin_bp = (self.origin_bp - self._step(True)) % self._total

    def action_rotate_ccw_lg(self):
        if self._total:
            self.origin_bp = (self.origin_bp + self._step(True)) % self._total

    def action_reset_origin(self):
        self.origin_bp = 0

    def action_aspect_inc(self):
        self._aspect = round(min(5.0, self._aspect + 0.05), 3)
        self.notify(f"Aspect {self._aspect:.2f}  (press , to widen)", timeout=1.5)

    def action_aspect_dec(self):
        self._aspect = round(max(0.5, self._aspect - 0.05), 3)
        self.notify(f"Aspect {self._aspect:.2f}  (press . to heighten)", timeout=1.5)

    def action_toggle_map_view(self):
        self._map_mode = "linear" if self._map_mode == "circular" else "circular"
        self.refresh()

    def select_feature(self, idx: int) -> None:
        self.selected_idx = idx
        self.refresh()

    # ── Mouse ──────────────────────────────────────────────────────────────────

    def on_mouse_scroll_up(self, _: MouseScrollUp):
        self.action_rotate_cw()

    def on_mouse_scroll_down(self, _: MouseScrollDown):
        self.action_rotate_ccw()

    def _feat_at(self, x: int, y: int) -> int:
        if not self.record or not self._total:
            return -1
        w, h = self.size.width, self.size.height
        cx, cy, rx, ry = self._geometry(w, h)
        if rx == 0 or ry == 0:
            return -1
        dc_n = (x - cx) / rx
        dr_n = (y - cy) / ry
        r_norm = math.sqrt(dc_n ** 2 + dr_n ** 2)
        if r_norm < 0.75 or r_norm > 1.35:
            return -1
        angle = math.atan2(dr_n, dc_n)
        angle_norm = (angle + math.pi / 2) % (2 * math.pi)
        bp = int(self.origin_bp + self._total * angle_norm / (2 * math.pi)) % self._total
        for i, f in enumerate(self._feats):
            if self._bp_in(bp, f):
                return i
        return -1

    def on_click(self, event: Click):
        if not self.record:
            return
        if self._map_mode == "linear":
            idx = self._feat_at_linear(event.x, event.y)
        else:
            idx = self._feat_at(event.x, event.y)
        self.selected_idx = idx
        f = self._feats[idx] if idx >= 0 else None
        self.post_message(self.FeatureSelected(idx, f))

    def _bp_in(self, bp: int, f: dict) -> bool:
        s, e = f["start"], f["end"]
        return (s <= bp < e) if e >= s else (bp >= s or bp < e)

    # ── Geometry ───────────────────────────────────────────────────────────────

    def _geometry(self, w: int, h: int) -> tuple[int, int, int, int]:
        cx     = w // 2
        cy     = h // 2
        aspect = self._aspect
        rx_from_w = cx - 16
        rx_from_h = int((cy - 3) * aspect)
        rx = max(8, min(rx_from_w, rx_from_h))
        ry = max(4, round(rx / aspect))
        return cx, cy, rx, ry

    def _angle_to_xy(self, angle: float, cx: int, cy: int,
                     rx: int, ry: int, dr: int = 0) -> tuple[int, int]:
        scale = ry / rx if rx else 0.5
        x = round(cx + (rx + dr) * math.cos(angle))
        y = round(cy + (ry + dr * scale) * math.sin(angle))
        return x, y

    def _bp_to_angle(self, bp: int) -> float:
        return 2 * math.pi * ((bp - self.origin_bp) % self._total) / self._total - math.pi / 2

    # ── Render ─────────────────────────────────────────────────────────────────

    # cache: (w, h, origin_bp, selected_idx, _aspect, n_feats, n_restr) → Text
    _render_cache: "tuple | None" = None

    def render(self) -> Text:
        if not self.record:
            return Text(
                "\n\n   No record loaded.\n"
                "   Press  f  to fetch from GenBank (e.g. L09137 for pUC19)\n"
                "   Press  o  to open a local .gb file\n",
                style="dim italic",
            )
        w, h = self.size.width, self.size.height
        if w < 30 or h < 14:
            return Text(f"  Window too small ({w}×{h})", style="dim red")
        key = (w, h, self.origin_bp, self.selected_idx, self._aspect,
               len(self._feats), len(self._restr_feats), self._map_mode,
               self._show_connectors)
        if self._render_cache and self._render_cache[0] == key:
            return self._render_cache[1]
        result = self._draw_linear(w, h) if self._map_mode == "linear" else self._draw(w, h)
        self._render_cache = (key, result)
        return result

    def _draw(self, w: int, h: int) -> Text:
        canvas = _Canvas(w, h)
        bc     = _BrailleCanvas(w, h)

        total  = self._total
        cx, cy, rx, ry = self._geometry(w, h)

        _ry_rx   = ry / rx if rx else 0.5
        TWO_PI   = 2 * math.pi
        HALF_PI  = math.pi / 2
        # Braille-space offsets for the circle centre
        cx2, cy4 = cx * 2, cy * 4
        rx2, ry4 = rx * 2, ry * 4

        # Enough steps to fill every braille column (2 dots/cell → circumference
        # in braille-px ≈ 2π*rx*2; ×2 safety margin → 4π*rx ≈ 13*rx)
        n_steps = max(360, int(13 * rx))

        # Local refs avoid repeated attribute lookups in hot loops
        bc_bits   = bc._bits
        bc_colors = bc._colors
        bc_prio   = bc._prio
        bc_cols   = bc.cols
        bc_rows   = bc.rows
        # Flat braille dot-bit lookup: index = (py & 3) << 1 | (px & 1)
        _BIT = (0, 3, 1, 4, 2, 5, 6, 7)

        def bp2a(bp: int) -> float:
            return self._bp_to_angle(bp)

        def a2xy(angle: float, dr: int = 0) -> tuple[int, int]:
            return self._angle_to_xy(angle, cx, cy, rx, ry, dr)

        # ── Backbone (inlined — no function call overhead) ─────────────────────
        _step_angle = TWO_PI / n_steps
        for step in range(n_steps):
            a  = step * _step_angle - HALF_PI
            px = round(cx2 + rx2 * math.cos(a))
            py = round(cy4 + ry4 * math.sin(a))
            pcol, prow = px >> 1, py >> 2
            if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                bc_bits[prow][pcol]   |= 1 << _BIT[(py & 3) << 1 | (px & 1)]
                if bc_prio[prow][pcol] == 0:
                    bc_colors[prow][pcol] = "color(238)"
                    bc_prio[prow][pcol]   = 1

        # ── Position ticks (inside the circle) ───────────────────────────────
        # TICK_DR_MARK  — radial inset for the ┼ graduation mark
        # TICK_DR_LABEL — radial inset for the bp number label
        # Both are negative so they sit inside the backbone ring.
        # These values scale automatically with the aspect ratio (, / . keys).
        TICK_DR_MARK  = -2
        TICK_DR_LABEL = -5

        tick_int = _nice_tick(total)
        bp = 0
        while bp < total:
            angle  = bp2a(bp)
            tx, ty = a2xy(angle, dr=TICK_DR_MARK)
            canvas.put(tx, ty, "┼", "color(250)")
            label  = _format_bp(bp)
            lx, ly = a2xy(angle, dr=TICK_DR_LABEL)
            # Inside the circle: text points inward, so alignment is flipped
            # vs outside placement (right-side labels hang left, left-side hang right).
            if math.cos(angle) >= 0:
                canvas.put_text(lx - len(label) + 1, ly, label, "color(245)")
            else:
                canvas.put_text(lx, ly, label, "color(245)")
            bp += tick_int

        # ── Restriction site marks ─────────────────────────────────────────────
        # resite: draw a thin arc just outside (+) or inside (-) the backbone
        # recut:  draw a radial tick crossing the backbone at the cut position
        restr_labels: list[tuple[float, str, str]] = []
        for rf in self._restr_feats:
            color = rf["color"]
            if rf["type"] == "resite":
                dr_lo = 4 if rf["strand"] >= 0 else -5
                dr_hi = dr_lo + 1
                start_a = bp2a(rf["start"])
                end_a   = bp2a(rf["end"])
                span_a  = (end_a - start_a) % TWO_PI or TWO_PI
                arc_steps = max(4, int(span_a / TWO_PI * n_steps))
                inv_arc   = 1.0 / arc_steps if arc_steps else 1.0
                for dr_i in range(2):
                    dr     = dr_lo + (dr_hi - dr_lo) * dr_i
                    _rx2dr = (rx + dr) * 2
                    _ry4dr = (ry + dr * _ry_rx) * 4
                    for s in range(arc_steps + 1):
                        a  = start_a + s * inv_arc * span_a
                        px = round(cx2 + _rx2dr * math.cos(a))
                        py = round(cy4 + _ry4dr * math.sin(a))
                        pcol, prow = px >> 1, py >> 2
                        if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                            bc_bits[prow][pcol]   |= 1 << _BIT[(py & 3) << 1 | (px & 1)]
                            if 1 >= bc_prio[prow][pcol]:
                                bc_colors[prow][pcol] = color
                                bc_prio[prow][pcol]   = 1
                # Collect label at midpoint
                mid_bp = (rf["start"] + rf["end"]) // 2
                restr_labels.append((bp2a(mid_bp), rf["label"], color))
            elif rf["type"] == "recut":
                # Radial tick crossing the backbone at cut position
                cut_a = bp2a(rf["start"])
                for dr in range(-1, 3):
                    tx, ty = a2xy(cut_a, dr=dr)
                    if 0 <= tx < w and 0 <= ty < h:
                        canvas.put(tx, ty, "┼" if dr == 0 else "·", color)

        # ── Features (large → small) ──────────────────────────────────────────
        feats_sorted = sorted(
            enumerate(self._feats),
            key=lambda iv: -(iv[1]["end"] - iv[1]["start"]),
        )
        label_slots: list[tuple[float, str, str]] = []
        N_DR = 4   # radial samples per arc (was 6; 4 is visually identical)

        for orig_idx, f in feats_sorted:
            is_sel   = (orig_idx == self.selected_idx)
            color    = f["color"]
            style    = ("reverse " + color) if is_sel else color
            f_prio   = 3 if is_sel else 2

            start_bp, end_bp = f["start"], f["end"]
            strand           = f["strand"]
            dr_lo = 1 if strand >= 0 else -3
            dr_hi = dr_lo + 2

            start_a = bp2a(start_bp)
            end_a   = bp2a(end_bp)
            span_a  = (end_a - start_a) % TWO_PI or TWO_PI

            arc_steps = max(8, int(span_a / TWO_PI * n_steps))
            inv_arc   = 1.0 / arc_steps if arc_steps else 1.0

            # Precompute cos/sin once; reuse across all N_DR radial passes
            arc_cos = [math.cos(start_a + s * inv_arc * span_a) for s in range(arc_steps + 1)]
            arc_sin = [math.sin(start_a + s * inv_arc * span_a) for s in range(arc_steps + 1)]

            inv_ndr1 = 1.0 / max(1, N_DR - 1)
            dr_span  = dr_hi - dr_lo
            for dr_i in range(N_DR):
                dr     = dr_lo + dr_span * dr_i * inv_ndr1
                _rx2dr = (rx + dr) * 2
                _ry4dr = (ry + dr * _ry_rx) * 4
                for ca, sa in zip(arc_cos, arc_sin):
                    px   = round(cx2 + _rx2dr * ca)
                    py   = round(cy4 + _ry4dr * sa)
                    pcol = px >> 1
                    prow = py >> 2
                    if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                        bc_bits[prow][pcol]   |= 1 << _BIT[(py & 3) << 1 | (px & 1)]
                        if f_prio >= bc_prio[prow][pcol]:
                            bc_colors[prow][pcol] = style
                            bc_prio[prow][pcol]   = f_prio

            if strand >= 0:
                tip_a, tip_tan = end_a, end_a + HALF_PI
            else:
                tip_a, tip_tan = start_a, start_a - HALF_PI
            ax, ay = a2xy(tip_a, dr=dr_lo + 1)
            canvas.put(ax, ay, _arrow_char(tip_tan), "bold " + style)

            mid_bp = (start_bp + (end_bp - start_bp) // 2) % total
            label_slots.append((bp2a(mid_bp), f["label"], color))

        # Add restriction site labels (from resite entries only)
        for angle, lbl, color in restr_labels:
            label_slots.append((angle, lbl, color))

        # ── Labels: place each as close to the arc as possible ───────────────
        # Greedily try increasing dr until the label's bounding box doesn't
        # overlap any already-placed label.
        # ↓ Tune this to control how far labels sit from the arc.
        LABEL_DR_MIN = 9          # minimum radial clearance from arc edge
        dr_min = LABEL_DR_MIN
        dr_max = max(rx // 2 + 6, LABEL_DR_MIN + 10)

        # placed: list of (x0, x1, row) bounding boxes of accepted labels
        placed_boxes: list[tuple[int, int, int]] = []
        final_labels: list[tuple[float, str, str, int, int, int]] = []
        # angle, lbl, color, chosen_dr, lx, ly

        for angle, lbl, color in label_slots:
            on_right = math.cos(angle) >= 0
            chosen = None
            for dr in range(dr_min, dr_max + 1):
                lx, ly = a2xy(angle, dr=dr)
                if not (0 <= ly < h):
                    continue
                lbl_x0 = lx if on_right else max(0, lx - len(lbl) + 1)
                lbl_x1 = lx + len(lbl) - 1 if on_right else lx
                # Check against every already-placed box
                ok = True
                for bx0, bx1, by in placed_boxes:
                    if by == ly and not (lbl_x1 < bx0 or lbl_x0 > bx1):
                        ok = False
                        break
                if ok:
                    chosen = (dr, lx, ly, lbl_x0, lbl_x1)
                    break
            if chosen is None:
                # Couldn't fit without overlap — place at max dr anyway
                lx, ly = a2xy(angle, dr=dr_max)
                lbl_x0 = lx if on_right else max(0, lx - len(lbl) + 1)
                lbl_x1 = lx + len(lbl) - 1 if on_right else lx
                chosen = (dr_max, lx, ly, lbl_x0, lbl_x1)
            dr_c, lx, ly, lbl_x0, lbl_x1 = chosen
            placed_boxes.append((lbl_x0, lbl_x1, ly))
            final_labels.append((angle, lbl, color, dr_c, lx, ly))

        # Render
        for angle, lbl, color, dr_c, lx, ly in final_labels:
            on_right = math.cos(angle) >= 0

            # Dot just outside the arc
            dot_x, dot_y = a2xy(angle, dr=3)
            canvas.put(dot_x, dot_y, "·", color)

            # Optional connector line from dot to label
            if self._show_connectors:
                lbl_mid_x = lx + len(lbl) // 2 if on_right else lx - len(lbl) // 2
                steps = max(1, abs(lbl_mid_x - dot_x) + abs(ly - dot_y))
                for t in range(1, steps):
                    px = dot_x + (lbl_mid_x - dot_x) * t // steps
                    py = dot_y + (ly        - dot_y) * t // steps
                    if 0 <= px < w and 0 <= py < h and canvas._chars[py][px] == " ":
                        canvas.put(px, py, "·", color)

            if on_right:
                canvas.put_text(lx, ly, lbl, color)
            else:
                canvas.put_text(max(0, lx - len(lbl) + 1), ly, lbl, color)

        # ── Center info ───────────────────────────────────────────────────────
        name     = (self.record.name or self.record.id or "?")[:w // 3]
        size_txt = f"{total:,} bp"
        orig_txt = f"▲ {self.origin_bp:,}"
        for i, (txt, sty) in enumerate([
            (name,     "bold white"),
            (size_txt, "color(245)"),
            (orig_txt, "dim cyan"),
        ]):
            canvas.put_text(cx - len(txt) // 2, cy - 1 + i, txt, sty)

        return bc.combine(canvas)

    # ── Linear map ─────────────────────────────────────────────────────────────

    def _feat_at_linear(self, x: int, y: int) -> int:
        """Return feature index at terminal cell (x, y) in linear view, or -1."""
        if not self._total:
            return -1
        w, h      = self.size.width, self.size.height
        margin_l  = 5
        margin_r  = 2
        usable_w  = w - margin_l - margin_r
        backbone_row = max(4, h // 2)
        if x < margin_l or x >= w - margin_r or usable_w <= 0:
            return -1
        bp = int((x - margin_l) / usable_w * self._total)
        above = y < backbone_row
        below = y > backbone_row
        for i, f in enumerate(self._feats):
            s, e = f["start"], f["end"]
            in_range = (s <= bp <= e) if e > s else (bp >= s or bp <= e)
            if not in_range:
                continue
            if above and f["strand"] >= 0:
                return i
            if below and f["strand"] < 0:
                return i
        return -1

    def _draw_linear(self, w: int, h: int) -> Text:
        """Render a horizontal linear plasmid map."""
        canvas = _Canvas(w, h)
        bc     = _BrailleCanvas(w, h)
        total  = self._total

        if not total:
            canvas.put_text(w // 2 - 9, h // 2, "No record loaded", "dim")
            return bc.combine(canvas)

        # ── Layout ──
        margin_l     = 5
        margin_r     = 2
        usable_w     = w - margin_l - margin_r
        px_w         = usable_w * 2
        px_start     = margin_l * 2
        backbone_row = max(4, h // 2)
        backbone_py  = backbone_row * 4 + 1   # braille pixel row

        _BIT      = (0, 3, 1, 4, 2, 5, 6, 7)
        bc_bits   = bc._bits
        bc_colors = bc._colors
        bc_prio   = bc._prio
        bc_cols   = bc.cols
        bc_rows   = bc.rows

        def bp_to_px(bp: int) -> int:
            return px_start + int(bp / total * px_w)

        def _set(px: int, py: int, color: str, prio: int) -> None:
            pcol, prow = px >> 1, py >> 2
            if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                bc_bits[prow][pcol]   |= 1 << _BIT[(py & 3) << 1 | (px & 1)]
                if prio >= bc_prio[prow][pcol]:
                    bc_colors[prow][pcol] = color
                    bc_prio[prow][pcol]   = prio

        # ── Backbone ──
        for px in range(px_start, px_start + px_w + 1):
            _set(px, backbone_py, "color(238)", 1)

        # ── Ticks + bp labels ──
        tick_int  = _nice_tick(total)
        label_row = backbone_row + 1
        bp = 0
        while bp <= total:
            tx = margin_l + bp * usable_w // total
            canvas.put(tx, backbone_row, "┼", "color(250)")
            lbl = _format_bp(bp)
            if label_row < h:
                canvas.put_text(tx - len(lbl) // 2, label_row, lbl, "color(245)")
            bp += tick_int
        # Right cap
        end_tx = min(margin_l + usable_w, w - 1)
        canvas.put(end_tx, backbone_row, "┤", "color(250)")

        # ── Restriction site marks ──
        # resite: thin braille arc above (fwd) or below (rev) backbone
        # recut:  radial tick at cut position
        for rf in self._restr_feats:
            color = rf["color"]
            if rf["type"] == "resite":
                x0 = margin_l + rf["start"] * usable_w // total
                x1 = margin_l + rf["end"]   * usable_w // total
                x0 = max(margin_l, min(x0, w - margin_r - 1))
                x1 = max(margin_l, min(x1, w - margin_r - 1))
                dy = -2 if rf["strand"] >= 0 else 2
                bar_py = backbone_py + dy * 4
                for bx in range(x0 * 2, x1 * 2 + 1):
                    pcol, prow = bx >> 1, bar_py >> 2
                    if 0 <= pcol < bc_cols and 0 <= prow < bc_rows:
                        bc_bits[prow][pcol]   |= 1 << _BIT[(bar_py & 3) << 1 | (bx & 1)]
                        if 1 >= bc_prio[prow][pcol]:
                            bc_colors[prow][pcol] = color
                            bc_prio[prow][pcol]   = 1
            elif rf["type"] == "recut":
                cut_x = margin_l + rf["start"] * usable_w // total
                if margin_l <= cut_x < w - margin_r:
                    row_above = backbone_row - 1
                    row_below = backbone_row + 1
                    canvas.put(cut_x, backbone_row, "┼", color)
                    if 0 <= row_above < h:
                        canvas.put(cut_x, row_above, "↓" if rf["strand"] >= 0 else " ", color)
                    if 0 <= row_below < h:
                        canvas.put(cut_x, row_below, "↑" if rf["strand"] < 0 else " ", color)

        # ── Lane assignment (greedy interval scheduling) ──
        fwd_ends: list[int] = []   # rightmost braille-x used per fwd lane
        rev_ends: list[int] = []
        feat_meta: list[tuple[bool, int]] = []   # (is_fwd, lane_idx)

        for feat in self._feats:
            is_fwd = feat["strand"] >= 0
            x0, x1 = bp_to_px(feat["start"]), bp_to_px(feat["end"])
            ends   = fwd_ends if is_fwd else rev_ends
            lane   = len(ends)
            for li, ex in enumerate(ends):
                if x0 > ex + 2:
                    lane = li
                    ends[li] = x1
                    break
            else:
                ends.append(x1)
            feat_meta.append((is_fwd, lane))

        # ── Draw features ──
        lane_gap    = 3    # braille rows between backbone and first lane
        bar_half    = 2    # half-height of feature bar in braille rows
        lane_stride = 7    # braille rows between lane centres

        for i, (feat, (is_fwd, lane)) in enumerate(zip(self._feats, feat_meta)):
            start_bp = feat["start"]
            end_bp   = feat["end"]
            strand   = feat["strand"]
            color    = feat["color"]
            label    = feat.get("label", feat.get("type", ""))

            # Handle wrap-around features
            if end_bp > start_bp:
                segments = [(bp_to_px(start_bp), bp_to_px(end_bp))]
            else:
                segments = [(bp_to_px(start_bp), px_start + px_w),
                            (px_start, bp_to_px(end_bp))]

            center_py = (backbone_py - lane_gap - bar_half - lane * lane_stride
                         if is_fwd else
                         backbone_py + lane_gap + bar_half + lane * lane_stride)
            bar_top    = center_py - bar_half
            bar_bottom = center_py + bar_half

            if bar_top < 0 or bar_bottom >= h * 4:
                continue

            is_sel = (i == self.selected_idx)
            prio   = 3 if is_sel else 2
            style  = ("reverse " + color) if is_sel else color
            feat_ty = center_py >> 2

            # Label row: above bar for forward, below for reverse
            label_ty = feat_ty - 1 if is_fwd else feat_ty + 1

            for sx0, sx1 in segments:
                sx0 = max(px_start, min(sx0, px_start + px_w))
                sx1 = max(px_start, min(sx1, px_start + px_w))
                if sx1 <= sx0:
                    continue

                # Fill bar in braille
                for py in range(bar_top, bar_bottom + 1):
                    for px in range(sx0, sx1):
                        _set(px, py, style, prio)

                # Arrowhead character
                if strand >= 0:
                    canvas.put(min(sx1 >> 1, w - margin_r - 1), feat_ty, "▶", style)
                else:
                    canvas.put(max(sx0 >> 1, margin_l), feat_ty, "◀", style)

                # Label above (fwd) or below (rev) the bar
                x0c, x1c = sx0 >> 1, sx1 >> 1
                max_lbl_w = x1c - x0c
                if max_lbl_w >= 1 and 0 <= label_ty < h:
                    lbl = label[:max_lbl_w]
                    lx  = x0c + (max_lbl_w - len(lbl)) // 2
                    canvas.put_text(lx, label_ty, lbl, style)

        # ── Header ──
        name = (self.record.name or self.record.id or "?")[:w // 3]
        canvas.put_text(margin_l, 0, f"{name}  {total:,} bp", "bold white")
        hint = "[ linear  ·  v = circular ]"
        canvas.put_text(w - len(hint) - 1, 0, hint, "dim")

        return bc.combine(canvas)


# ── Sidebar feature table ──────────────────────────────────────────────────────

class FeatureSidebar(Widget):
    DEFAULT_CSS = """
    FeatureSidebar {
        width: 32;
        border-left: solid $primary;
    }
    #feat-table  { height: 1fr; }
    #detail-box  { height: 8; border-top: solid $accent; padding: 0 1; }
    #sidebar-hdr { background: $primary; padding: 0 1; }
    """

    class RowActivated(Message):
        def __init__(self, idx: int):
            self.idx = idx
            super().__init__()

    def compose(self) -> ComposeResult:
        yield Static(" Features", id="sidebar-hdr")
        yield DataTable(id="feat-table", cursor_type="row", zebra_stripes=True)
        yield Static("", id="detail-box")

    def on_mount(self):
        t = self.query_one("#feat-table", DataTable)
        t.add_columns("Type", "Label", "bp", "±")

    def populate(self, feats: list[dict]) -> None:
        t = self.query_one("#feat-table", DataTable)
        t.clear()
        for f in feats:
            strand_sym = "+" if f["strand"] == 1 else ("−" if f["strand"] == -1 else "·")
            bp_str     = f"{f['start']+1}‥{f['end']}"
            t.add_row(
                Text(f["type"][:12],  style=f["color"]),
                Text(f["label"][:14], style=f["color"]),
                bp_str,
                strand_sym,
            )

    def show_detail(self, f: dict | None) -> None:
        box = self.query_one("#detail-box", Static)
        if f is None:
            box.update(Text(""))
            return
        strand_sym = "+" if f["strand"] == 1 else ("−" if f["strand"] == -1 else "·")
        span = f["end"] - f["start"]
        t = Text()
        t.append(f["type"],  style=f"bold {f['color']}")
        t.append("\n")
        t.append(f["label"], style="white")
        t.append("\n")
        t.append(f"{f['start']+1}‥{f['end']} ({span:,} bp)", style="dim")
        t.append("\n")
        t.append(f"Strand: {strand_sym}", style="dim")
        box.update(t)

    _prog_row: int = -1   # cursor moves driven by highlight_row, not the user

    def highlight_row(self, idx: int) -> None:
        """Move cursor to row; suppresses the resulting RowActivated echo."""
        if idx < 0:
            return
        t = self.query_one("#feat-table", DataTable)
        self._prog_row = idx
        try:
            t.move_cursor(row=idx)
        except Exception:
            self._prog_row = -1

    @on(DataTable.RowHighlighted, "#feat-table")
    def _row_highlighted(self, event: DataTable.RowHighlighted):
        if event.cursor_row == self._prog_row:
            self._prog_row = -1
            return
        self.post_message(self.RowActivated(event.cursor_row))


# ── Library panel ──────────────────────────────────────────────────────────────

class LibraryPanel(Widget):
    """Left-hand plasmid library — persistent SnapGene-style collection."""

    DEFAULT_CSS = """
    LibraryPanel {
        width: 24;
        border-right: solid $primary;
    }
    #lib-hdr   { background: $primary; padding: 0 1; }
    #lib-table { height: 1fr; }
    #lib-btns  { height: 3; }
    #lib-btns Button { min-width: 5; margin: 0 0 0 1; }
    """

    class PlasmidLoad(Message):
        """User selected a library entry to load."""
        def __init__(self, entry: dict):
            self.entry = entry
            super().__init__()

    class AddCurrentRequested(Message):
        """User pressed '+' to add the currently-loaded record."""
        pass

    def compose(self) -> ComposeResult:
        yield Static(" Library", id="lib-hdr")
        yield DataTable(id="lib-table", cursor_type="row", zebra_stripes=True)
        with Horizontal(id="lib-btns"):
            yield Button("+", id="btn-lib-add", variant="success",
                         tooltip="Add current plasmid")
            yield Button("−", id="btn-lib-del", variant="error",
                         tooltip="Remove selected")

    def on_mount(self):
        t = self.query_one("#lib-table", DataTable)
        t.add_columns("Name", "bp")
        self._repopulate()

    def _repopulate(self):
        t = self.query_one("#lib-table", DataTable)
        t.clear()
        for entry in _load_library():
            t.add_row(
                entry["name"][:14],
                f"{entry['size']:,}",
                key=entry["id"],
            )

    def add_entry(self, record) -> None:
        """Serialize record and persist to library JSON."""
        gb_text = _record_to_gb_text(record)
        entries = _load_library()
        entries = [e for e in entries if e.get("id") != record.id]
        import datetime
        entries.insert(0, {
            "name":    record.name or record.id,
            "id":      record.id,
            "size":    len(record.seq),
            "n_feats": len([f for f in record.features if f.type != "source"]),
            "source":  getattr(record, "_tui_source", f"id:{record.id}"),
            "added":   datetime.date.today().isoformat(),
            "gb_text": gb_text,
        })
        _save_library(entries)
        self._repopulate()

    @on(DataTable.RowSelected, "#lib-table")
    def _row_selected(self, event: DataTable.RowSelected):
        key = event.row_key.value if event.row_key else None
        if key is None:
            return
        for entry in _load_library():
            if entry.get("id") == key:
                self.post_message(self.PlasmidLoad(entry))
                return

    @on(Button.Pressed, "#btn-lib-add")
    def _btn_add(self):
        self.post_message(self.AddCurrentRequested())

    def set_dirty(self, dirty: bool) -> None:
        """Show/hide unsaved-changes marker in the panel header."""
        self.query_one("#lib-hdr", Static).update(
            " * Library" if dirty else " Library"
        )

    @on(Button.Pressed, "#btn-lib-del")
    def _btn_del(self):
        t = self.query_one("#lib-table", DataTable)
        if t.row_count == 0:
            return
        row_keys = list(t.rows.keys())
        if not (0 <= t.cursor_row < len(row_keys)):
            return
        rk = row_keys[t.cursor_row]
        entries = [e for e in _load_library() if e.get("id") != rk.value]
        _save_library(entries)
        self._repopulate()


# ── Sequence panel ─────────────────────────────────────────────────────────────

class SequencePanel(Widget):
    """
    Bottom DNA sequence viewer.

    Click  on the sequence → place cursor + select feature at that position.
    Shift+click           → extend selection for editing.
    Shift+E               → open insert/replace dialog at cursor / selection.
    """

    DEFAULT_CSS = """
    SequencePanel {
        height: 14;
        border-top: solid $primary;
    }
    #seq-hdr    { background: $primary; padding: 0 1; height: 1; }
    #seq-scroll { height: 1fr; }
    #seq-trans  {
        height: 3; border-top: solid $primary-darken-2;
        padding: 0 1; display: none;
    }
    SequencePanel.has-trans #seq-trans { display: block; }
    """

    # ── Messages ───────────────────────────────────────────────────────────────

    class SequenceChanged(Message):
        """Emitted when the sequence is modified (commit=True = full rebuild)."""
        def __init__(self, seq: str, commit: bool = False):
            self.seq    = seq
            self.commit = commit
            super().__init__()

    class SequenceClick(Message):
        """User clicked on a base; app should select feature there."""
        def __init__(self, bp: int, double: bool = False):
            self.bp     = bp
            self.double = double
            super().__init__()

    # ── Lifecycle ──────────────────────────────────────────────────────────────

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._seq:          str                     = ""
        self._feats:        list[dict]              = []
        self._sel_range:    "tuple[int,int] | None" = None  # feature highlight
        self._user_sel:     "tuple[int,int] | None" = None  # drag/shift selection
        self._cursor_pos:   int                     = -1    # -1 = no cursor
        self._view_cache_key: "tuple | None"        = None
        self._view_cache_txt: "Text | None"         = None
        self._show_connectors:  bool = False
        # Drag-to-select state
        self._drag_start_bp:    int  = -1
        self._has_dragged:      bool = False
        self._mouse_button_held: bool = False
        self._drag_was_shift:   bool = False
        self._last_was_drag:    bool = False
        # Set by _click_to_bp when the click lands on a resite bar row
        self._last_resite_click: "dict | None" = None

    def compose(self) -> ComposeResult:
        yield Static(
            " Sequence  (click: select · Shift+click: select region · Shift+E: edit)",
            id="seq-hdr",
        )
        with ScrollableContainer(id="seq-scroll"):
            yield Static("", id="seq-view")
        yield Static("", id="seq-trans")

    # ── Public API ─────────────────────────────────────────────────────────────

    def update_seq(self, seq: str, feats: list[dict]) -> None:
        """Called after loading a record or committing an edit."""
        self._seq        = seq
        self._feats      = feats
        self._sel_range  = None
        self._user_sel   = None
        self._cursor_pos = -1
        self.remove_class("has-trans")
        self._refresh_view()

    def highlight_feature(self, feat: "dict | None") -> None:
        """Highlight a feature's region in the sequence; show CDS translation."""
        if feat is None or not self._seq:
            self._sel_range = None
            self.remove_class("has-trans")
            self._refresh_view()
            return

        start, end = feat["start"], min(feat["end"], len(self._seq))
        self._sel_range = (start, end)
        self._user_sel  = None          # clear shift-selection on programmatic highlight
        self._refresh_view()

        # Scroll to feature start (accounting for annotation bar rows)
        self._scroll_to_row(self._bp_to_content_row(start))

        trans_box = self.query_one("#seq-trans", Static)
        if feat.get("type") == "CDS":
            aa   = _translate_cds(self._seq, start, end, feat.get("strand", 1))
            n_aa = len(aa.rstrip("*"))
            t    = Text()
            t.append(f" {feat['label']}  ", style=f"bold {feat['color']}")
            t.append(f"({n_aa} aa)  ", style="dim")
            t.append(aa[:self.size.width - 4], style=feat["color"])
            trans_box.update(t)
            self.add_class("has-trans")
        else:
            self.remove_class("has-trans")

    def select_feature_range(self, feat: dict) -> None:
        """Double-click highlight: set user_sel to the entire feature span."""
        if not self._seq or feat is None:
            return
        start = feat["start"]
        end   = min(feat["end"], len(self._seq))
        self._user_sel   = (start, end)
        self._sel_range  = None
        self._cursor_pos = start
        self._refresh_view()
        self._scroll_to_row(self._bp_to_content_row(start))

    # ── Mouse / click ──────────────────────────────────────────────────────────

    def on_mouse_down(self, event: MouseDown) -> None:
        if event.button != 1:
            return
        bp = self._click_to_bp(event.screen_x, event.screen_y)
        if bp < 0:
            return
        self._mouse_button_held = True
        self._drag_start_bp     = bp
        self._has_dragged       = False
        self._drag_was_shift    = event.shift
        if event.shift and self._cursor_pos >= 0:
            # Shift+click: extend selection from cursor to here
            s = min(self._cursor_pos, bp)
            e = max(self._cursor_pos, bp) + 1
            self._user_sel   = (s, e)
            self._cursor_pos = bp
            self._sel_range  = None
        else:
            # Plain click: place cursor, clear selection
            self._cursor_pos = bp
            self._user_sel   = None
        self._refresh_view()
        self._scroll_to_row(self._bp_to_content_row(bp))

    def on_mouse_move(self, event: MouseMove) -> None:
        if not self._mouse_button_held or self._drag_start_bp < 0:
            return
        bp = self._click_to_bp(event.screen_x, event.screen_y)
        if bp < 0 or bp == self._drag_start_bp and not self._has_dragged:
            return
        self._has_dragged = True
        s = min(self._drag_start_bp, bp)
        e = max(self._drag_start_bp, bp) + 1
        self._user_sel   = (s, e)
        self._cursor_pos = bp
        self._sel_range  = None
        self._refresh_view()

    def on_mouse_up(self, event: MouseUp) -> None:
        if event.button != 1:
            return
        self._last_was_drag     = self._has_dragged
        self._mouse_button_held = False
        self._drag_start_bp     = -1
        self._has_dragged       = False
        self._drag_was_shift    = False

    def on_click(self, event: Click) -> None:
        """Handle single/double click — fires after mouse_up at same position."""
        if self._last_was_drag:
            self._last_was_drag = False
            return
        self._last_was_drag = False
        self._last_resite_click = None
        bp = self._click_to_bp(event.screen_x, event.screen_y)
        if bp < 0:
            return

        # If the click landed on a restriction site bar, highlight that span
        resite = self._last_resite_click
        self._last_resite_click = None
        if resite is not None:
            self._sel_range  = (resite["start"], min(resite["end"], len(self._seq)))
            self._user_sel   = None
            self._cursor_pos = resite["start"]
            self._refresh_view()
            self._scroll_to_row(self._bp_to_content_row(resite["start"]))
            return

        double = event.chain >= 2
        self.post_message(self.SequenceClick(bp, double=double))

    def _click_to_bp(self, screen_x: int, screen_y: int) -> int:
        """Map absolute screen coords to a bp index, or -1 if not on a base."""
        if not self._seq:
            return -1
        try:
            scroll      = self.query_one("#seq-scroll", ScrollableContainer)
            reg         = scroll.region          # absolute screen position of viewport
            vp_x        = screen_x - reg.x       # column within the viewport
            vp_y        = screen_y - reg.y       # visible row within the viewport
            if vp_x < 0 or vp_y < 0 or vp_x >= reg.width or vp_y >= reg.height:
                return -1
            content_row = vp_y + int(scroll.scroll_y)
        except Exception:
            return -1

        line_width  = max(20, self.size.width - 16)
        annot_feats = sorted(
            [f for f in self._feats if f.get("type") not in ("site", "recut")],
            key=lambda f: -(f["end"] - f["start"]),
        )
        rpg = 2 + (1 if self._show_connectors else 0)  # rows per feature group
        n   = len(self._seq)
        row = 0
        seq_col     = vp_x - 10   # offset for feature bar rows (prefix_w = 10)
        seq_col_dna = vp_x - 13   # offset for DNA rows (10 + 3 for '5'─' / '3'─')

        for chunk_start in range(0, n, line_width):
            chunk_end   = min(chunk_start + line_width, n)
            chunk_feats = [f for f in annot_feats
                           if f["start"] < chunk_end and f["end"] > chunk_start]
            above_lanes, below_lanes = _assign_chunk_features(chunk_feats, chunk_start, chunk_end)

            # Above feature rows (forward strand)
            for lane in above_lanes:
                for _ in range(rpg):
                    if row == content_row:
                        for f in lane:
                            bar_s = max(f["start"], chunk_start) - chunk_start
                            bar_e = min(f["end"],   chunk_end)   - chunk_start
                            if bar_s <= seq_col < bar_e:
                                if f.get("type") == "resite":
                                    self._last_resite_click = f
                                return (f["start"] + f["end"]) // 2
                        return lane[0]["start"]
                    row += 1

            # DNA rows: fwd strand + RC strand (2 rows)
            for _ in range(2):
                if row == content_row:
                    if 0 <= seq_col_dna < (chunk_end - chunk_start):
                        return chunk_start + seq_col_dna
                    return -1
                row += 1

            # Below feature rows (reverse strand)
            for lane in below_lanes:
                for _ in range(rpg):
                    if row == content_row:
                        for f in lane:
                            bar_s = max(f["start"], chunk_start) - chunk_start
                            bar_e = min(f["end"],   chunk_end)   - chunk_start
                            if bar_s <= seq_col < bar_e:
                                if f.get("type") == "resite":
                                    self._last_resite_click = f
                                return (f["start"] + f["end"]) // 2
                        return lane[0]["start"]
                    row += 1

            if row > content_row:
                break
        return -1

    # ── Helpers ────────────────────────────────────────────────────────────────

    def _annot_feats_sorted(self) -> list:
        return sorted(
            [f for f in self._feats if f.get("type") not in ("site", "recut")],
            key=lambda f: -(f["end"] - f["start"]),
        )

    def _bp_to_content_row(self, bp: int) -> int:
        """Return the content row index (0-based) of the DNA line containing bp."""
        line_width  = max(20, self.size.width - 16)
        annot_feats = self._annot_feats_sorted()
        rpg = 2 + (1 if self._show_connectors else 0)
        n   = len(self._seq)
        row = 0
        for chunk_start in range(0, n, line_width):
            chunk_end   = min(chunk_start + line_width, n)
            chunk_feats = [f for f in annot_feats
                           if f["start"] < chunk_end and f["end"] > chunk_start]
            above_lanes, below_lanes = _assign_chunk_features(chunk_feats, chunk_start, chunk_end)
            above_rows = len(above_lanes) * rpg
            if bp < chunk_end:
                return row + above_rows   # forward-strand DNA row within this chunk
            row += above_rows + 2 + len(below_lanes) * rpg
        return row

    def _scroll_to_row(self, row: int) -> None:
        try:
            self.query_one("#seq-scroll", ScrollableContainer).scroll_to(
                0, row, animate=False
            )
        except Exception:
            pass

    def _refresh_view(self) -> None:
        view = self.query_one("#seq-view", Static)
        if not self._seq:
            view.update(Text("  No sequence loaded.", style="dim italic"))
            return
        # 10-char line-number prefix + "5'─" (3) + seq + "─3'" (3) = prefix 16
        line_width = max(20, self.size.width - 16)
        key = (id(self._seq), id(self._feats), line_width,
               self._sel_range, self._user_sel, self._cursor_pos,
               self._show_connectors)
        if key != self._view_cache_key:
            self._view_cache_txt = _build_seq_text(
                self._seq, self._feats,
                line_width      = line_width,
                sel_range       = self._sel_range,
                user_sel        = self._user_sel,
                cursor_pos      = self._cursor_pos,
                show_connectors = self._show_connectors,
            )
            self._view_cache_key = key
        view.update(self._view_cache_txt)

    def on_resize(self, _) -> None:
        self._refresh_view()


# ── Sequence edit dialog ───────────────────────────────────────────────────────

class EditSeqDialog(ModalScreen):
    """Insert or replace a sequence region via a small dialog.

    mode="insert"  — inserts new_seq at position start
    mode="replace" — replaces seq[start:end] with new_seq
    Dismisses with (new_seq, mode, start, end) on OK, or None on cancel.
    """

    _VALID = frozenset("ATCGNRYSWKMBDHV")   # IUPAC DNA codes

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def __init__(self, mode: str, existing: str = "",
                 start: int = 0, end: int = 0):
        super().__init__()
        self._mode     = mode
        self._existing = existing
        self._start    = start
        self._end      = end

    def compose(self) -> ComposeResult:
        with Vertical(id="edit-dlg"):
            if self._mode == "insert":
                yield Static(
                    f" Insert sequence at position {self._start + 1} ",
                    id="edit-title",
                )
            else:
                excerpt = self._existing[:35] + ("…" if len(self._existing) > 35 else "")
                yield Static(
                    f" Replace {self._end - self._start} bp  "
                    f"({self._start + 1}‥{self._end}) ",
                    id="edit-title",
                )
                yield Static(f"  Current: {excerpt}", id="edit-current")
            yield Label("New sequence  (A T C G N and IUPAC codes only):")
            yield Input(id="edit-input", placeholder="ATCG…")
            yield Static("", id="edit-err")
            with Horizontal(id="edit-btns"):
                yield Button("OK", id="btn-ok", variant="primary")
                yield Button("Cancel", id="btn-cancel")

    def on_mount(self) -> None:
        self.query_one("#edit-input", Input).focus()

    @on(Input.Changed, "#edit-input")
    def _live_validate(self, event: Input.Changed) -> None:
        val = event.value.strip().upper()
        err = self.query_one("#edit-err", Static)
        if not val:
            err.update("")
            return
        bad = sorted(set(c for c in val if c not in self._VALID))
        if bad:
            chars = "  ".join(repr(c) for c in bad)
            err.update(Text(
                f"Invalid: {chars} — only A T C G N (IUPAC) allowed",
                style="bold red",
            ))
        else:
            err.update(Text(f"{len(val)} bp", style="dim green"))

    def _try_submit(self) -> None:
        val = self.query_one("#edit-input", Input).value.strip().upper()
        err = self.query_one("#edit-err", Static)
        if not val:
            err.update(Text("Please enter a sequence.", style="bold red"))
            return
        bad = [c for c in val if c not in self._VALID]
        if bad:
            return   # live validation already shows the error
        self.dismiss((val, self._mode, self._start, self._end))

    @on(Button.Pressed, "#btn-ok")
    def _ok(self, _) -> None:
        self._try_submit()

    @on(Button.Pressed, "#btn-cancel")
    def _cancel_btn(self, _) -> None:
        self.dismiss(None)

    @on(Input.Submitted, "#edit-input")
    def _submitted(self, _) -> None:
        self._try_submit()

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Fetch modal ────────────────────────────────────────────────────────────────

class FetchModal(ModalScreen):
    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def compose(self) -> ComposeResult:
        with Vertical(id="fetch-box"):
            yield Static(" Fetch GenBank Record ", id="fetch-title")
            yield Label("NCBI Accession  (e.g.  M77789  for pUC19,  Y14837  for pUC57):")
            yield Input(placeholder="MW463917.1", id="fetch-acc", value="MW463917.1")
            yield Label("Email  (required by NCBI Entrez):")
            yield Input(placeholder="you@example.com", id="fetch-email")
            with Horizontal(id="fetch-btns"):
                yield Button("Fetch", id="btn-fetch", variant="primary")
                yield Button("Cancel", id="btn-cancel-fetch")
            yield Static("", id="fetch-status", markup=True)

    @on(Button.Pressed, "#btn-fetch")
    def _fetch(self):
        acc   = self.query_one("#fetch-acc",   Input).value.strip()
        email = self.query_one("#fetch-email", Input).value.strip() or "splicecraft@local"
        if not acc:
            self.query_one("#fetch-status", Static).update("[red]Enter an accession.[/red]")
            return
        self.query_one("#fetch-status", Static).update(
            f"[dim]Fetching {acc!r} from NCBI…[/dim]"
        )
        self._do_fetch(acc, email)

    @work(thread=True)
    def _do_fetch(self, acc: str, email: str):
        try:
            record = fetch_genbank(acc, email)
            self.app.call_from_thread(self.dismiss, record)
        except Exception as exc:
            def _err():
                self.query_one("#fetch-status", Static).update(
                    f"[red]Error: {exc}[/red]"
                )
            self.app.call_from_thread(_err)

    @on(Input.Submitted)
    def _submitted(self):
        self.query_one("#btn-fetch", Button).press()

    @on(Button.Pressed, "#btn-cancel-fetch")
    def _cancel_btn(self):
        self.dismiss(None)

    def action_cancel(self):
        self.dismiss(None)


# ── Open-file modal ────────────────────────────────────────────────────────────

class OpenFileModal(ModalScreen):
    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def compose(self) -> ComposeResult:
        with Vertical(id="open-box"):
            yield Static(" Open GenBank File ", id="open-title")
            yield Label("File path  (.gb / .gbk):")
            yield Input(placeholder="/path/to/plasmid.gb", id="open-path")
            with Horizontal(id="open-btns"):
                yield Button("Open", id="btn-open", variant="primary")
                yield Button("Cancel", id="btn-cancel-open")
            yield Static("", id="open-status", markup=True)

    @on(Button.Pressed, "#btn-open")
    def _open(self):
        path = self.query_one("#open-path", Input).value.strip()
        if not path:
            self.query_one("#open-status", Static).update("[red]Enter a file path.[/red]")
            return
        try:
            record = load_genbank(path)
            record._tui_source = path   # remember where it came from
            self.dismiss(record)
        except Exception as exc:
            self.query_one("#open-status", Static).update(f"[red]{exc}[/red]")

    @on(Input.Submitted)
    def _submitted(self):
        self.query_one("#btn-open", Button).press()

    @on(Button.Pressed, "#btn-cancel-open")
    def _cancel_btn(self):
        self.dismiss(None)

    def action_cancel(self):
        self.dismiss(None)


# ── Dropdown menu modal ────────────────────────────────────────────────────────

class DropdownScreen(ModalScreen):
    """Transparent overlay showing a positioned dropdown menu."""

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    DEFAULT_CSS = """
    DropdownScreen {
        background: transparent;
    }
    """

    def __init__(self, items: list, x: int, y: int) -> None:
        super().__init__()
        self._items = items   # (label, action_str | None)
        self._x = x
        self._y = y
        self._cursor = next(
            (i for i, (_, a) in enumerate(items) if a is not None), 0
        )

    def compose(self) -> ComposeResult:
        yield Static(
            self._render_content(),
            id="dropdown-box",
        )

    def on_mount(self) -> None:
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        box_h   = len(self._items) + 2
        box = self.query_one("#dropdown-box", Static)
        box.styles.offset = (self._x, self._y)
        box.styles.width  = inner_w
        box.styles.height = box_h
        box.styles.border = ("solid", "white")
        box.styles.background = "black"

    def _render_content(self) -> Text:
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        sep_line = "\u2500" * (inner_w - 2)
        result = Text()
        for i, (label, action) in enumerate(self._items):
            is_sep      = (label == "---")
            is_selected = (i == self._cursor and not is_sep and action is not None)
            is_disabled = (action is None and not is_sep)

            if is_sep:
                line = Text(sep_line + "\n", style="white")
            else:
                padded = f" {label:<{inner_w - 3}}"
                if is_selected:
                    line = Text(padded + "\n", style="reverse white")
                elif is_disabled:
                    line = Text(padded + "\n", style="dim white")
                else:
                    line = Text(padded + "\n", style="white")
            result.append_text(line)
        return result

    def _refresh_box(self) -> None:
        box = self.query_one("#dropdown-box", Static)
        box.update(self._render_content())

    def on_key(self, event) -> None:
        items = self._items
        if event.key == "up":
            pos = self._cursor - 1
            while pos >= 0 and (items[pos][0] == "---" or items[pos][1] is None):
                pos -= 1
            if pos >= 0:
                self._cursor = pos
                self._refresh_box()
            event.stop()
        elif event.key == "down":
            pos = self._cursor + 1
            while pos < len(items) and (items[pos][0] == "---" or items[pos][1] is None):
                pos += 1
            if pos < len(items):
                self._cursor = pos
                self._refresh_box()
            event.stop()
        elif event.key == "enter":
            label, action = items[self._cursor]
            if action is not None:
                self.dismiss(action)
            event.stop()

    def on_click(self, event: Click) -> None:
        bx = self._x
        by = self._y
        inner_w = max((len(lbl) for lbl, _ in self._items), default=10) + 4
        bh = len(self._items) + 2
        cx, cy = event.screen_x, event.screen_y
        if bx <= cx < bx + inner_w and by <= cy < by + bh:
            row_in_box = cy - by - 1  # -1 for top border
            if 0 <= row_in_box < len(self._items):
                label, action = self._items[row_in_box]
                if label == "---" or action is None:
                    event.stop()
                    return
                self.dismiss(action)
        else:
            self.dismiss(None)
        event.stop()

    def action_cancel(self) -> None:
        self.dismiss(None)


# ── Menu bar widget ────────────────────────────────────────────────────────────

class MenuBar(Widget):
    DEFAULT_CSS = """
    MenuBar {
        height: 1;
        background: $primary-darken-3;
        layout: horizontal;
    }
    MenuBar Static {
        width: auto;
        padding: 0 2;
        color: $text;
    }
    MenuBar Static:hover {
        background: $primary;
    }
    """

    MENUS = ["File", "Edit", "Enzymes", "Features", "Primers", "Genes"]

    def compose(self) -> ComposeResult:
        for name in self.MENUS:
            yield Static(name, classes="menu-item", id=f"menu-{name.lower()}")

    def on_click(self, event: Click) -> None:
        for name in self.MENUS:
            widget_id = f"menu-{name.lower()}"
            try:
                item = self.query_one(f"#{widget_id}", Static)
                region = item.region
                if (region.x <= event.screen_x < region.x + region.width and
                        region.y <= event.screen_y < region.y + region.height):
                    x = region.x
                    y = region.y + 1
                    self.app.open_menu(name, x, y)
                    break
            except Exception:
                pass


# ── Unsaved-changes quit dialog ────────────────────────────────────────────────

class UnsavedQuitModal(ModalScreen):
    """Shown when the user tries to quit with unsaved edits."""

    BINDINGS = [Binding("escape", "cancel", "Cancel")]

    def compose(self) -> ComposeResult:
        with Vertical(id="quit-dlg"):
            yield Static(" Unsaved Changes ", id="quit-title")
            yield Static(
                "  You have unsaved edits. What would you like to do?",
                id="quit-msg",
            )
            with Horizontal(id="quit-btns"):
                yield Button("Save & Quit",      id="btn-save-quit", variant="primary")
                yield Button("Abandon Changes",  id="btn-abandon",   variant="error")
                yield Button("Cancel",           id="btn-cancel-quit")

    @on(Button.Pressed, "#btn-save-quit")
    def _save_quit(self, _): self.dismiss("save")

    @on(Button.Pressed, "#btn-abandon")
    def _abandon(self, _):   self.dismiss("abandon")

    @on(Button.Pressed, "#btn-cancel-quit")
    def _cancel_btn(self, _): self.dismiss(None)

    def action_cancel(self): self.dismiss(None)


# ── Main app ───────────────────────────────────────────────────────────────────

class PlasmidApp(App):
    TITLE       = "SpliceCraft"
    TRANSITIONS = {}          # instant screen open/close — no slide animations
    _preload_record = None
    _current_record = None   # last-loaded SeqRecord
    _source_path:   "str | None" = None   # file the current record was loaded from
    _unsaved:        bool         = False  # True when there are unsaved edits
    _MAX_UNDO = 50
    _restr_unique_only: bool = True
    _restr_min_len: int = 6

    CSS = """
Screen { background: $background; }

/* ── Layout ─────────────────────────────────────────────── */
MenuBar { height: 1; dock: top; }
#main-row   { height: 1fr; }
#center-col { width: 1fr; height: 1fr; }
#map-row    { height: 1fr; }

#status-bar {
    height: 1;
    background: $primary-darken-2;
    color: $text-muted;
    padding: 0 1;
}

/* ── Fetch modal ─────────────────────────────────────────── */
FetchModal { align: center middle; }
#fetch-box {
    width: 70; height: auto; max-height: 90%;
    background: $surface; border: solid $accent; padding: 1 2;
}
#fetch-title { background: $accent; padding: 0 1; margin-bottom: 1; }
#fetch-box Label { color: $text-muted; margin-top: 1; }
#fetch-btns { height: 3; margin-top: 1; }
#fetch-btns Button { margin-right: 1; }
#fetch-status { height: 1; margin-top: 1; }

/* ── Open-file modal ─────────────────────────────────────── */
OpenFileModal { align: center middle; }
#open-box {
    width: 70; height: auto;
    background: $surface; border: solid $primary; padding: 1 2;
}
#open-title { background: $primary; padding: 0 1; margin-bottom: 1; }
#open-box Label { color: $text-muted; margin-top: 1; }
#open-btns { height: 3; margin-top: 1; }
#open-btns Button { margin-right: 1; }
#open-status { height: 1; margin-top: 1; }

/* ── Edit-sequence dialog ────────────────────────────────── */
EditSeqDialog { align: center middle; }
#edit-dlg {
    width: 72; height: auto;
    background: $surface; border: solid $warning; padding: 1 2;
}
#edit-title   { background: $warning-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#edit-current { color: $text-muted; margin-bottom: 1; }
#edit-dlg Label { color: $text-muted; margin-top: 1; }
#edit-err  { height: 2; margin-top: 1; }
#edit-btns { height: 3; margin-top: 1; }
#edit-btns Button { margin-right: 1; }

/* ── Unsaved-quit dialog ─────────────────────────────────── */
UnsavedQuitModal { align: center middle; }
#quit-dlg {
    width: 60; height: auto;
    background: $surface; border: solid $error; padding: 1 2;
}
#quit-title { background: $error-darken-2; color: $text; padding: 0 1; margin-bottom: 1; }
#quit-msg   { color: $text-muted; margin-bottom: 1; }
#quit-btns  { height: 3; margin-top: 1; }
#quit-btns Button { margin-right: 1; }
"""

    BINDINGS = [
        Binding("f",           "fetch",            "Fetch GenBank", show=True),
        Binding("o",           "open_file",        "Open .gb file", show=True),
        Binding("a",           "add_to_library",   "Add to lib",    show=True),
        Binding("E",           "edit_seq",         "Edit seq",      show=True,  key_display="E"),
        Binding("S",           "save",             "Save",          show=True,  key_display="S"),
        Binding("[",           "rotate_cw",        "Rotate ←",      show=True,  priority=True),
        Binding("]",           "rotate_ccw",       "Rotate →",      show=True,  priority=True),
        Binding("shift+[",     "rotate_cw_lg",     "Rotate ←←",     show=False, priority=True),
        Binding("shift+]",     "rotate_ccw_lg",    "Rotate →→",     show=False, priority=True),
        Binding("home",        "reset_origin",     "Reset origin",  show=True,  priority=True),
        Binding("v",           "toggle_map_view",  "⊙/─ View",      show=True,  priority=True),
        Binding("l",           "toggle_connectors","Connectors",    show=True,  priority=True),
        Binding("delete",      "delete_feature",   "Del feature",   show=True,  priority=True),
        Binding("q",           "quit",             "Quit",          show=True),
    ]

    def compose(self) -> ComposeResult:
        yield Header()
        yield MenuBar()
        with Horizontal(id="main-row"):
            yield LibraryPanel(id="library")
            with Vertical(id="center-col"):
                with Horizontal(id="map-row"):
                    yield PlasmidMap(id="plasmid-map")
                    yield FeatureSidebar(id="sidebar")
                yield SequencePanel(id="seq-panel")
        yield Static(
            Text(
                "  [ ] rotate   ← → cursor/map   Shift coarse   Home reset"
                "   f fetch   o open   a add-to-lib   E edit seq   S save   , / . circle",
                style="color(245)",
                no_wrap=True,
            ),
            id="status-bar",
        )
        yield Footer()

    # ── Delegate rotation keys to PlasmidMap ───────────────────────────────────

    def action_rotate_cw(self):
        self.query_one("#plasmid-map", PlasmidMap).action_rotate_cw()

    def action_rotate_ccw(self):
        self.query_one("#plasmid-map", PlasmidMap).action_rotate_ccw()

    def action_rotate_cw_lg(self):
        self.query_one("#plasmid-map", PlasmidMap).action_rotate_cw_lg()

    def action_rotate_ccw_lg(self):
        self.query_one("#plasmid-map", PlasmidMap).action_rotate_ccw_lg()

    def action_reset_origin(self):
        self.query_one("#plasmid-map", PlasmidMap).action_reset_origin()

    def action_toggle_map_view(self):
        self.query_one("#plasmid-map", PlasmidMap).action_toggle_map_view()

    def action_toggle_connectors(self):
        sp = self.query_one("#seq-panel", SequencePanel)
        pm = self.query_one("#plasmid-map", PlasmidMap)
        sp._show_connectors  = not sp._show_connectors
        pm._show_connectors  = sp._show_connectors
        sp._view_cache_key   = None   # invalidate seq panel cache
        pm._render_cache     = None   # invalidate map cache
        sp._refresh_view()
        pm.refresh()
        state = "on" if sp._show_connectors else "off"
        self.notify(f"Label connectors {state}")

    def action_edit_seq(self) -> None:
        sp = self.query_one("#seq-panel", SequencePanel)
        if not sp._seq:
            self.notify("No sequence loaded.", severity="warning")
            return
        if sp._user_sel is not None:
            # Replace the shift-selected region
            s, e     = sp._user_sel
            existing = sp._seq[s:e]
            self.push_screen(
                EditSeqDialog("replace", existing, s, e),
                callback=self._edit_dialog_result,
            )
        elif sp._cursor_pos >= 0:
            # Insert at cursor position
            pos = sp._cursor_pos
            self.push_screen(
                EditSeqDialog("insert", start=pos, end=pos),
                callback=self._edit_dialog_result,
            )
        else:
            self.notify(
                "Click on the sequence to place a cursor, "
                "or Shift+click to select a region.",
                severity="information",
            )

    def _edit_dialog_result(self, result) -> None:
        if result is None:
            return
        self._push_undo()
        new_bases, mode, s, e = result
        sp      = self.query_one("#seq-panel", SequencePanel)
        pm      = self.query_one("#plasmid-map", PlasmidMap)
        old_seq = sp._seq

        if mode == "insert":
            new_seq = old_seq[:s] + new_bases + old_seq[s:]
        else:
            new_seq = old_seq[:s] + new_bases + old_seq[e:]

        new_cursor = s + len(new_bases)

        if self._current_record is not None:
            new_record = self._rebuild_record_with_edit(new_seq, mode, s, e, new_bases)
            self._current_record = new_record
            pm.load_record(new_record)
            self.query_one("#sidebar", FeatureSidebar).populate(pm._feats)
            restr = _scan_restriction_sites(new_seq, min_recognition_len=self._restr_min_len, unique_only=self._restr_unique_only)
            pm._restr_feats = restr
            pm.refresh()
            sp.update_seq(new_seq, pm._feats + restr)
            self.notify(f"Sequence updated  ({len(new_seq):,} bp)")
        else:
            restr = _scan_restriction_sites(new_seq, min_recognition_len=self._restr_min_len, unique_only=self._restr_unique_only)
            pm._restr_feats = restr
            pm.refresh()
            sp.update_seq(new_seq, pm._feats + restr)

        self._mark_dirty()

        # Restore cursor after update_seq resets it
        sp._cursor_pos = new_cursor
        sp._user_sel   = None
        sp._refresh_view()

    def _rebuild_record_with_edit(self, new_seq: str, mode: str,
                                   s: int, e: int, new_bases: str):
        """Rebuild SeqRecord after an insert/replace, shifting feature coords precisely."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        ins_len  = len(new_bases)
        del_len  = 0 if mode == "insert" else (e - s)
        delta    = ins_len - del_len
        new_len  = len(new_seq)
        # For insert: the edit region is [s, s); e == s
        # For replace: the edit region is [s, e)
        edit_end = s if mode == "insert" else e

        new_record = SeqRecord(
            Seq(new_seq),
            id=self._current_record.id,
            name=self._current_record.name,
            description=self._current_record.description,
            annotations=dict(self._current_record.annotations),
        )

        for feat in self._current_record.features:
            fs = int(feat.location.start)
            fe = int(feat.location.end)

            if mode == "insert":
                if fe <= s:
                    new_fs, new_fe = fs, fe           # entirely before insert
                elif fs >= s:
                    new_fs, new_fe = fs + ins_len, fe + ins_len  # entirely after
                else:
                    new_fs, new_fe = fs, fe + ins_len  # spans insert → expand
            else:  # replace [s, e)
                if fe <= s:
                    new_fs, new_fe = fs, fe             # entirely before
                elif fs >= e:
                    new_fs, new_fe = fs + delta, fe + delta  # entirely after
                elif fs <= s and fe >= e:
                    new_fs, new_fe = fs, fe + delta     # spans replaced region
                elif fs < s:
                    new_fs, new_fe = fs, s + ins_len    # overlaps start of region
                else:
                    new_fs, new_fe = s, fe + delta      # overlaps end of region

            new_fe = max(new_fs + 1, min(new_fe, new_len))
            new_record.features.append(SeqFeature(
                FeatureLocation(new_fs, new_fe,
                                strand=getattr(feat.location, "strand", 1)),
                type=feat.type,
                qualifiers=dict(feat.qualifiers),
            ))

        return new_record

    def _rebuild_record_without_feature(self, feat_idx: int):
        """Create a new SeqRecord with the feat_idx-th non-source feature removed."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation

        new_record = SeqRecord(
            Seq(str(self._current_record.seq)),
            id=self._current_record.id,
            name=self._current_record.name,
            description=self._current_record.description,
            annotations=dict(self._current_record.annotations),
        )
        non_source_idx = 0
        for feat in self._current_record.features:
            if feat.type == "source":
                new_record.features.append(SeqFeature(
                    feat.location, type=feat.type,
                    qualifiers=dict(feat.qualifiers),
                ))
                continue
            if non_source_idx != feat_idx:
                new_record.features.append(SeqFeature(
                    feat.location, type=feat.type,
                    qualifiers=dict(feat.qualifiers),
                ))
            non_source_idx += 1
        return new_record

    def action_delete_feature(self):
        pm = self.query_one("#plasmid-map", PlasmidMap)
        if pm.selected_idx < 0 or not pm._feats:
            self.notify("No feature selected.", severity="warning")
            return
        feat = pm._feats[pm.selected_idx]
        label = feat.get("label") or feat.get("type", "feature")
        self._push_undo()
        new_record = self._rebuild_record_without_feature(pm.selected_idx)
        sp = self.query_one("#seq-panel", SequencePanel)
        self._apply_snapshot(str(new_record.seq), sp._cursor_pos, new_record)
        self.notify(f"Deleted '{label}'  (Ctrl+Z to undo)")

    def action_add_to_library(self):
        if self._current_record is None:
            self.notify("No record loaded to add.", severity="warning")
            return
        lib = self.query_one("#library", LibraryPanel)
        lib.add_entry(self._current_record)
        self.notify(f"Added {self._current_record.name} to library.")

    # ── Mount: auto-load preloaded record ──────────────────────────────────────

    def on_mount(self) -> None:
        self._undo_stack: list = []
        self._redo_stack: list = []
        if self._preload_record is not None:
            def _load_preload():
                record = self._preload_record
                path   = getattr(record, "_tui_source", None)
                self._apply_record(record)
                self._source_path = path
            self.call_after_refresh(_load_preload)
        elif not _load_library():
            self._seed_default_library()

    @work(thread=True)
    def _seed_default_library(self) -> None:
        """Fetch MW463917.1 and pre-populate the library on first run."""
        try:
            record = fetch_genbank("MW463917.1")
            def _add():
                lib = self.query_one("#library", LibraryPanel)
                lib.add_entry(record)
                self._apply_record(record)
            self.call_from_thread(_add)
        except Exception:
            pass

    # ── Keyboard: cursor movement, copy, undo/redo ─────────────────────────────

    def on_key(self, event) -> None:
        sp = self.query_one("#seq-panel", SequencePanel)

        # ── Ctrl+C: copy selection ────────────────────────────────────────────
        if event.key == "ctrl+c":
            seq = sp._seq
            if seq:
                sel = sp._user_sel or sp._sel_range
                if sel:
                    text = seq[sel[0]:sel[1]].upper()
                    if _copy_to_clipboard(text):
                        self.notify(f"Copied {len(text)} bp to clipboard")
                    else:
                        self.notify("Clipboard unavailable (install xclip / xsel / wl-copy)",
                                    severity="warning")
                else:
                    self.notify("No selection — Shift+click or double-click a feature first",
                                severity="information")
            event.stop()
            return

        # ── Ctrl+Z: undo ──────────────────────────────────────────────────────
        if event.key == "ctrl+z":
            self._action_undo()
            event.stop()
            return

        # ── Ctrl+Shift+Z / Ctrl+Y: redo ───────────────────────────────────────
        if event.key in ("ctrl+shift+z", "ctrl+Z", "ctrl+y"):
            self._action_redo()
            event.stop()
            return

        # ── Arrow keys: move sequence cursor ──────────────────────────────────
        if sp._cursor_pos < 0 or not sp._seq:
            return
        if event.key == "left":
            new_pos = max(0, sp._cursor_pos - 1)
        elif event.key == "right":
            new_pos = min(len(sp._seq), sp._cursor_pos + 1)
        else:
            return
        event.stop()
        sp._cursor_pos = new_pos
        sp._user_sel   = None
        sp._refresh_view()
        sp._scroll_to_row(sp._bp_to_content_row(new_pos))

    # ── Undo / redo ────────────────────────────────────────────────────────────

    def _push_undo(self) -> None:
        sp = self.query_one("#seq-panel", SequencePanel)
        if not sp._seq:
            return
        snapshot = (sp._seq, sp._cursor_pos, self._current_record)
        self._undo_stack.append(snapshot)
        if len(self._undo_stack) > self._MAX_UNDO:
            self._undo_stack.pop(0)
        self._redo_stack.clear()

    def _apply_snapshot(self, seq: str, cursor_pos: int, record) -> None:
        pm      = self.query_one("#plasmid-map",  PlasmidMap)
        sidebar = self.query_one("#sidebar",      FeatureSidebar)
        sp      = self.query_one("#seq-panel",    SequencePanel)
        self._current_record = record
        if record is not None:
            pm.load_record(record)
            sidebar.populate(pm._feats)
            restr = _scan_restriction_sites(seq, min_recognition_len=self._restr_min_len, unique_only=self._restr_unique_only)
            pm._restr_feats = restr
            pm.refresh()
            sp.update_seq(seq, pm._feats + restr)
        else:
            sp.update_seq(seq, [])
        sp._cursor_pos = cursor_pos
        sp._refresh_view()
        self._mark_dirty()

    def _action_undo(self) -> None:
        if not self._undo_stack:
            self.notify("Nothing to undo", severity="information")
            return
        sp = self.query_one("#seq-panel", SequencePanel)
        # Save current state to redo stack
        self._redo_stack.append((sp._seq, sp._cursor_pos, self._current_record))
        seq, cursor_pos, record = self._undo_stack.pop()
        self._apply_snapshot(seq, cursor_pos, record)
        remaining = len(self._undo_stack)
        self.notify(f"Undo  ({remaining} left)")

    def _action_redo(self) -> None:
        if not self._redo_stack:
            self.notify("Nothing to redo", severity="information")
            return
        sp = self.query_one("#seq-panel", SequencePanel)
        self._undo_stack.append((sp._seq, sp._cursor_pos, self._current_record))
        seq, cursor_pos, record = self._redo_stack.pop()
        self._apply_snapshot(seq, cursor_pos, record)
        remaining = len(self._redo_stack)
        self.notify(f"Redo  ({remaining} left)")

    # ── Dirty-state helpers ────────────────────────────────────────────────────

    def _mark_dirty(self) -> None:
        self._unsaved = True
        if self._current_record:
            n = len(self._current_record.seq)
            self.title = f"SpliceCraft — *{self._current_record.name}  ({n:,} bp)"
        try:
            self.query_one("#library", LibraryPanel).set_dirty(True)
        except Exception:
            pass

    def _mark_clean(self) -> None:
        self._unsaved = False
        if self._current_record:
            n = len(self._current_record.seq)
            self.title = f"SpliceCraft — {self._current_record.name}  ({n:,} bp)"
        try:
            self.query_one("#library", LibraryPanel).set_dirty(False)
        except Exception:
            pass

    def _do_save(self) -> bool:
        """Write current record to its source file. Returns True on success."""
        if self._current_record is None:
            self.notify("Nothing to save.", severity="warning")
            return False
        if not self._source_path:
            self.notify(
                "No source file — open a .gb file first, or use 'a' to add to library.",
                severity="warning",
            )
            return False
        try:
            Path(self._source_path).write_text(_record_to_gb_text(self._current_record))
            self._mark_clean()
            self.notify(f"Saved → {self._source_path}")
            return True
        except Exception as exc:
            self.notify(f"Save failed: {exc}", severity="error")
            return False

    def action_save(self) -> None:
        self._do_save()

    def action_quit(self) -> None:
        if self._unsaved:
            self.push_screen(UnsavedQuitModal(), callback=self._on_quit_response)
        else:
            self.exit()

    def _on_quit_response(self, result) -> None:
        if result == "save":
            if self._do_save():
                self.exit()
        elif result == "abandon":
            self.exit()
        # None = cancel → stay

    # ── Fetch / open ───────────────────────────────────────────────────────────

    def action_fetch(self):
        self.push_screen(FetchModal(), callback=self._apply_record)

    def action_open_file(self):
        def _cb(record):
            if record is None:
                return
            path = getattr(record, "_tui_source", None)
            self._apply_record(record)
            self._source_path = path   # restore after _apply_record clears it
        self.push_screen(OpenFileModal(), callback=_cb)

    # ── Central record loader ──────────────────────────────────────────────────

    def _apply_record(self, record) -> None:
        """Load a SeqRecord into all panels."""
        if record is None:
            return
        self._current_record = record
        self._source_path    = None   # caller sets this if it came from a file

        pm      = self.query_one("#plasmid-map", PlasmidMap)
        sidebar = self.query_one("#sidebar",     FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel",   SequencePanel)

        pm.load_record(record)
        sidebar.populate(pm._feats)

        seq_str = str(record.seq)
        restr   = _scan_restriction_sites(seq_str, min_recognition_len=self._restr_min_len, unique_only=self._restr_unique_only)

        # Store restriction sites on the map for visual overlay
        pm._restr_feats = restr
        pm.refresh()

        # Sequence panel: feature coloring = record feats + restriction sites
        seq_pnl.update_seq(seq_str, pm._feats + restr)

        self._mark_clean()
        self.notify(
            f"Loaded {record.name}  ({len(record.seq):,} bp, "
            f"{len(pm._feats)} features, {len(restr)} restriction sites)"
        )

    # ── Feature selection: map ↔ sidebar ↔ sequence panel ─────────────────────

    @on(SequencePanel.SequenceClick)
    def _seq_click(self, event: SequencePanel.SequenceClick) -> None:
        """Single or double click on the sequence — select the smallest feature at bp."""
        pm      = self.query_one("#plasmid-map", PlasmidMap)
        sidebar = self.query_one("#sidebar",     FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel",   SequencePanel)
        bp      = event.bp
        best_idx  = -1
        best_span = float("inf")
        for i, f in enumerate(pm._feats):
            s, e = f["start"], f["end"]
            if s <= bp < e and (e - s) < best_span:
                best_span = e - s
                best_idx  = i
        if best_idx >= 0:
            f = pm._feats[best_idx]
            pm.select_feature(best_idx)
            sidebar.show_detail(f)
            sidebar.highlight_row(best_idx)
            if event.double:
                # Double-click: highlight entire feature span as editable selection
                seq_pnl.select_feature_range(f)
            else:
                seq_pnl.highlight_feature(f)

    @on(PlasmidMap.FeatureSelected)
    def _map_feat_selected(self, event: PlasmidMap.FeatureSelected):
        sidebar = self.query_one("#sidebar",   FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel", SequencePanel)
        sidebar.show_detail(event.feat_dict)
        if event.idx >= 0:
            sidebar.highlight_row(event.idx)
        seq_pnl.highlight_feature(event.feat_dict)

    @on(FeatureSidebar.RowActivated)
    def _sidebar_row_activated(self, event: FeatureSidebar.RowActivated):
        pm      = self.query_one("#plasmid-map", PlasmidMap)
        sidebar = self.query_one("#sidebar",     FeatureSidebar)
        seq_pnl = self.query_one("#seq-panel",   SequencePanel)
        pm.select_feature(event.idx)
        f = pm._feats[event.idx] if 0 <= event.idx < len(pm._feats) else None
        sidebar.show_detail(f)
        seq_pnl.highlight_feature(f)

    # ── Library events ─────────────────────────────────────────────────────────

    @on(LibraryPanel.PlasmidLoad)
    def _library_load(self, event: LibraryPanel.PlasmidLoad):
        gb_text = event.entry.get("gb_text", "")
        if not gb_text:
            self.notify(f"Library entry has no stored sequence.", severity="warning")
            return
        try:
            record = _gb_text_to_record(gb_text)
            self._apply_record(record)
        except Exception as exc:
            self.notify(f"Failed to load from library: {exc}", severity="error")

    @on(LibraryPanel.AddCurrentRequested)
    def _library_add_current(self, _):
        self.action_add_to_library()

    # ── Menu bar ───────────────────────────────────────────────────────────────

    def _rescan_restrictions(self) -> None:
        """Re-scan restriction sites with current settings and update UI."""
        sp = self.query_one("#seq-panel", SequencePanel)
        pm = self.query_one("#plasmid-map", PlasmidMap)
        if not sp._seq:
            return
        restr = _scan_restriction_sites(
            sp._seq,
            min_recognition_len=self._restr_min_len,
            unique_only=self._restr_unique_only,
        )
        pm._restr_feats = restr
        pm.refresh()
        sp.update_seq(sp._seq, pm._feats + restr)

    def open_menu(self, name: str, x: int, y: int) -> None:
        """Build menu item list for name and push DropdownScreen."""
        ck = "\u2713"  # checkmark
        nc = " "
        u  = ck if self._restr_unique_only else nc
        m6 = ck if self._restr_min_len == 6  else nc
        m4 = ck if self._restr_min_len == 4  else nc

        menus = {
            "File": [
                ("Open .gb file",   "open_file"),
                ("Fetch from NCBI", "fetch"),
                ("---",             None),
                ("Add to Library",  "add_to_library"),
                ("Save",            "save"),
                ("---",             None),
                ("Quit",            "quit"),
            ],
            "Edit": [
                ("Edit Sequence",   "edit_seq"),
                ("---",             None),
                ("Undo",            "undo"),
                ("Redo",            "redo"),
                ("---",             None),
                ("Delete Feature",  "delete_feature"),
            ],
            "Enzymes": [
                (f"[{u}] Unique cutters",   "toggle_restr_unique"),
                (f"[{m6}] 6+ bp sites",     "toggle_restr_min6"),
                (f"[{m4}] 4+ bp sites",     "toggle_restr_min4"),
                ("---",                      None),
                ("Toggle connectors",        "toggle_connectors"),
            ],
            "Features": [
                ("Add Feature...",   "add_feature"),
                ("Delete Feature",   "delete_feature"),
                ("---",              None),
                ("Toggle connectors","toggle_connectors"),
            ],
            "Primers": [
                ("Design Primer... (coming soon)", None),
            ],
            "Genes": [
                ("Annotate from NCBI... (coming soon)", None),
            ],
        }
        items = menus.get(name, [])
        self.push_screen(
            DropdownScreen(items, x, y),
            callback=self._menu_action,
        )

    def _menu_action(self, action: "str | None") -> None:
        if action is None:
            return
        # Handle toggle actions directly since they need state updates
        if action in ("toggle_restr_unique", "toggle_restr_min6", "toggle_restr_min4"):
            getattr(self, f"action_{action}")()
        else:
            self.call_action(action)

    def action_toggle_restr_unique(self) -> None:
        self._restr_unique_only = not self._restr_unique_only
        self._rescan_restrictions()
        state = "on" if self._restr_unique_only else "off"
        self.notify(f"Unique cutters {state}")

    def action_toggle_restr_min6(self) -> None:
        self._restr_min_len = 6
        self._rescan_restrictions()
        self.notify("Showing 6+ bp recognition sites")

    def action_toggle_restr_min4(self) -> None:
        self._restr_min_len = 4
        self._rescan_restrictions()
        self.notify("Showing 4+ bp recognition sites")

    def action_add_feature(self) -> None:
        self.notify("Add feature: coming soon", severity="information")

    def action_undo(self) -> None:
        self._action_undo()

    def action_redo(self) -> None:
        self._action_redo()

    # ── Sequence edits ─────────────────────────────────────────────────────────

    @on(SequencePanel.SequenceChanged)
    def _seq_changed(self, event: SequencePanel.SequenceChanged):
        # Update restriction site overlay whenever sequence changes
        pm    = self.query_one("#plasmid-map", PlasmidMap)
        restr = _scan_restriction_sites(event.seq, min_recognition_len=self._restr_min_len, unique_only=self._restr_unique_only)
        pm._restr_feats = restr
        pm.refresh()


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    app = PlasmidApp()

    if arg:
        if Path(arg).exists():
            try:
                record = load_genbank(arg)
                record._tui_source = str(Path(arg).resolve())
            except Exception as exc:
                print(f"Could not load {arg!r}: {exc}", file=sys.stderr)
                sys.exit(1)
        else:
            print(f"Fetching {arg!r} from NCBI…", flush=True)
            try:
                record = fetch_genbank(arg)
                print(f"  Got: {record.name}  ({len(record.seq)} bp)")
            except Exception as exc:
                print(f"Fetch failed: {exc}", file=sys.stderr)
                sys.exit(1)
        app._preload_record = record

    app.run()


if __name__ == "__main__":
    main()
