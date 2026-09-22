"""splicecraft_presets — the built-in feature preset catalogue (L0, pure data).

A shipped, read-only library of the DNA elements that show up in almost every
plasmid: resistance markers, origins, promoters, terminators, polyadenylation
signals, reporters, affinity tags, protease sites, recombination sites. The
point is that a user who has never annotated a plasmid by hand can still say
"what is this region?" or "put a T7 promoter here" without hunting a sequence
down and retyping it.

WHERE THE SEQUENCES COME FROM
=============================
Not from memory, and not from a commercial feature database. Every sequence in
this file was extracted from a public GenBank record, by this route:

1. ~400 annotated vector records were fetched from NCBI nuccore (cloning,
   expression, shuttle, binary, lentiviral, yeast, plant and reporter vectors,
   plus the canonical primary records: pUC19 L09137, pBR322 J01749, pACYC184
   X06403, pEGFP-N1 U55762, pGEX-4T-1 U13853, SV40 J02400, T7 V01146,
   f1 J02448, CaMV V00141, Tn5 U00004, and others).
2. Every annotated feature in that corpus was indexed by its label, and the
   sequence was extracted with ``SeqFeature.extract`` so a feature annotated on
   the minus strand is stored in its own functional orientation.
3. For each element, the sequence kept is the variant carried by the MOST
   INDEPENDENT ACCESSIONS under that label. ``source`` records the reference
   accession and how many independent records agree byte-for-byte. Most entries
   are corroborated by 3 or more unrelated submissions; the 15 that are backed
   by a single authoritative record say so.
4. Coding entries were checked mechanically against Biopython's translator:
   length divisible by three, a recognised start codon, a terminal stop, no
   internal stop. Each translation was then read against the protein its label
   claims (N-terminal residues and length) — which is how the two contaminated
   candidates below were caught. `tests/test_feature_presets.py` re-runs the
   mechanical half on every test run and pins the whole catalogue with an
   aggregate digest, so no base can change silently.

Two entries are a documented SUBSEQUENCE of an annotated feature rather than
the whole feature, because the source record folds vector context into the same
annotation. Both say so in ``source`` and neither is edited — only trimmed:

* ``GST tag`` — the record's CDS runs on through the vector thrombin site, a
  FLAG tag, the polylinker and a stop codon, which would silently terminate an
  N-terminal fusion. Only the GST open reading frame is kept.
* ``TEV protease site`` — sliced out of a His6+TEV tag block.

Curation caught a contaminated candidate that a from-memory dataset would
have shipped: a ``tdTomato`` whose annotation ran through a C-terminal V5 tag
in every corpus record (dropped — no clean copy existed).

Label-count consensus has one blind spot, and it shipped once: several public
records label an APH(3')-IIIa kanamycin kinase "aadA (kanamycin resistance)".
Counting labels made that the "aadA" winner, and the first catalogue named a
kanamycin gene as the spectinomycin marker. The entry is now
``KanR (aph(3')-IIIa)``, and ``SmR/SpecR (aadA)`` is a separately curated
AadA1. The lesson is the guard: a resistance marker's identity is checked by
TRANSLATING it and comparing the protein, never by its most common label.

WHAT THIS MODULE IS NOT
=======================
It is not user data. Presets live in code, never in ``features.json``, and
nothing here is ever written to the data directory. ``_load_features`` remains
the user's own library, unchanged. The two meet only in
``_merge_presets_with_library``, which builds a read-only view for browsing and
for annotation scans, and in ``_preset_to_library_entry``, which converts one
preset into an ordinary library entry when the user explicitly imports it.
A user entry always wins over a preset of the same (name, feature_type).

Layer 0 — pure stdlib, no sibling imports, no state.
"""

from __future__ import annotations

import copy

# Schema version for the preset entry shape. Bump when a field is added or
# renamed so a consumer can branch on it.
_PRESET_SCHEMA_VERSION = 1

# Display order for the category filter.
_PRESET_CATEGORIES: tuple[str, ...] = (
    "Promoter (bacterial)",
    "Promoter (mammalian)",
    "Promoter (yeast)",
    "Promoter (plant)",
    "Operator",
    "Regulator",
    "Translation",
    "Terminator",
    "polyA signal",
    "Origin",
    "Resistance",
    "Reporter",
    "Tag",
    "Protease site",
    "Linker / 2A",
    "Localisation",
    "Recombination",
    "Viral element",
    "CRISPR",
    "Cloning site",
)


# The catalogue. Entry shape mirrors a feature-library entry
# (name / feature_type / strand / color / sequence / description) so the
# existing row-build + annotate code paths work on a preset unchanged, plus
# three preset-only fields: `category` (browser filter), `aliases` (search
# synonyms) and `source` (provenance).
_FEATURE_PRESETS: list[dict] = [

    # ── Promoter (bacterial) ───────────────────────────────────────────
    {
        "name":         "AmpR promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Native promoter driving bla in pUC/pBR322-family "
                        "vectors. Constitutive and weak; sits immediately "
                        "upstream of the AmpR coding sequence.",
        "source":       "GenBank MT891329.1; identical sequence annotated 'AmpR "
                        "promoter' in 37 independent records",
        "aliases":      ["bla promoter"],
        "sequence":     "CGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGA"
                        "CAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGT",
    },
    {
        "name":         "araBAD promoter (pBAD)",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "L-arabinose-inducible promoter of the E. coli araBAD "
                        "operon. Needs AraC in the same cell. Titratable over a "
                        "wide range, and repressed by glucose, which makes it "
                        "the usual choice for toxic proteins.",
        "source":       "GenBank OR900359.1, annotated '@arabad'",
        "aliases":      ["pBAD", "araBAD", "arabinose"],
        "sequence":     "AAGAAACCAATTGTCCATATTGCATCAGACATTGCCGTCACTGCGTCTTTTACTGGCTCT"
                        "TCTCGCTAACCAAACCGGTAACCCCGCTTATTAAAAGCATTCTGTAACAAAGCGGGACCA"
                        "AAGCCATGACAAAAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTG"
                        "ATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGG"
                        "ATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCAT",
    },
    {
        "name":         "CmR promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Native promoter driving cat on pACYC-family vectors. "
                        "Constitutive.",
        "source":       "GenBank PZ357609.1; identical sequence annotated 'cat "
                        "promoter' in 8 independent records",
        "aliases":      ["cat promoter"],
        "sequence":     "TGATCGGCACGTAAGAGGTTCCAACTTTCACCATAATGAAATAAGATCACTACCGGGCGT"
                        "ATTTTTTGAGTTATCGAGATTTTCAGGAGCTAAGGAAGCTAAA",
    },
    {
        "name":         "EM7 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Synthetic constitutive bacterial promoter. Usually "
                        "placed in front of a resistance gene so one cassette "
                        "selects in both E. coli and mammalian cells.",
        "source":       "GenBank LT009451.1; identical sequence annotated "
                        "'synthetic EM7 promoter' in 4 independent records",
        "aliases":      ["EM7"],
        "sequence":     "GTTGACAATTAATCATCGGCATAGTATATCGGCATAGTATAATACGACAAGGTGAGGAAC"
                        "TAAACC",
    },
    {
        "name":         "J23119 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Strongest member of the Anderson constitutive promoter "
                        "family (BBa_J23119). No inducer, no repressor; the "
                        "reference promoter for iGEM-style part "
                        "characterisation.",
        "source":       "GenBank PZ357609.1; identical sequence annotated "
                        "'J23119 promoter' in 9 independent records",
        "aliases":      ["Anderson", "BBa_J23119", "constitutive"],
        "sequence":     "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC",
    },
    {
        "name":         "lac promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "E. coli lac promoter including the CAP site. Weak and "
                        "IPTG-inducible when LacI is present; drives lacZ-alpha "
                        "in pUC-family vectors for blue/white screening.",
        "source":       "GenBank LT009450.1; identical sequence annotated 'lac "
                        "promoter' in 10 independent records",
        "aliases":      ["lac", "Plac"],
        "sequence":     "GCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTAT"
                        "GCTTCCGGCTCGTATGTTGTGTGG",
    },
    {
        "name":         "pTetL promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Tetracycline-repressed promoter from the Tn10 tet "
                        "regulon. Repressed by TetR and de-repressed by "
                        "anhydrotetracycline.",
        "source":       "GenBank PZ366119.1; identical sequence annotated "
                        "'pTetL' in 4 independent records",
        "aliases":      ["tet promoter", "pTet"],
        "sequence":     "CGTTCAACAAACGGGCCATATTGTTGTATAAGTGATGAAATACTGAATTTAAAACTTAGT"
                        "TTATATGTGGTAAAATGTTTTAATCAAGTTTAGGAGGAATTAATTATGAAGTGTAATGAA"
                        "TAATGAATGTAACAGGGTTCAATTAAAAGAGGGAAGCGTATCATTAACCCTATAAACTAC"
                        "GTCTGCCCTCATTATTGGAGGGTGAAAT",
    },
    {
        "name":         "SP6 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Bacteriophage SP6 promoter. Used with SP6 RNA "
                        "polymerase for in-vitro transcription and probe "
                        "synthesis.",
        "source":       "GenBank MT891329.1; identical sequence annotated 'SP6 "
                        "promoter' in 6 independent records",
        "aliases":      ["SP6"],
        "sequence":     "ATTTAGGTGACACTATAGA",
    },
    {
        "name":         "T3 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Bacteriophage T3 promoter. Paired with a T7 or SP6 "
                        "promoter on the other side of a polylinker for "
                        "strand-specific in-vitro transcription.",
        "source":       "GenBank LT009450.1; identical sequence annotated 'T3 "
                        "promoter' in 10 independent records",
        "aliases":      ["T3"],
        "sequence":     "AATTAACCCTCACTAAA",
    },
    {
        "name":         "T7 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Bacteriophage T7 class III (phi10) promoter. "
                        "Transcribed only by T7 RNA polymerase, so expression "
                        "needs a DE3 lysogen or a T7 polymerase plasmid. "
                        "Transcription starts at the G immediately 3' of this "
                        "element.",
        "source":       "GenBank LT009450.1; identical sequence annotated 'T7 "
                        "promoter' in 19 independent records",
        "aliases":      ["T7", "phi10"],
        "sequence":     "TAATACGACTCACTATA",
    },
    {
        "name":         "tac promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Hybrid of the trp -35 and lac -10 elements. Much "
                        "stronger than lac and still IPTG-inducible. The "
                        "promoter in pGEX-family expression vectors.",
        "source":       "GenBank U13853.1, annotated 'tac'",
        "aliases":      ["tac", "Ptac"],
        "sequence":     "TTGACAATTAATCATCGGCTCGTATAATG",
    },
    {
        "name":         "trc promoter",
        "feature_type": "promoter",
        "category":     "Promoter (bacterial)",
        "strand":       1,
        "color":        "#00CED1",
        "description":  "Hybrid trp/lac promoter, a one-base spacer variant of "
                        "tac. Strong and IPTG-inducible.",
        "source":       "GenBank PZ357631.1; identical sequence annotated 'trc "
                        "promoter' in 8 independent records",
        "aliases":      ["trc", "Ptrc"],
        "sequence":     "TTGACAATTAATCATCCGGCTCGTATAATG",
    },

    # ── Promoter (mammalian) ───────────────────────────────────────────
    {
        "name":         "CMV enhancer",
        "feature_type": "regulatory",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "Human cytomegalovirus immediate-early enhancer. Placed "
                        "directly upstream of the CMV promoter; together they "
                        "give the strongest common mammalian expression.",
        "source":       "GenBank LC897328.1; identical sequence annotated 'CMV "
                        "enhancer' in 13 independent records",
        "aliases":      ["CMV", "hCMV", "enhancer"],
        "sequence":     "GACATTGATTATTGACTAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCC"
                        "CATATATGGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCA"
                        "ACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGA"
                        "CTTTCCATTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCACTTGGCAGTACATC"
                        "AAGTGTATCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCT"
                        "GGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTAT"
                        "TAGTCATCGCTATTACCATG",
    },
    {
        "name":         "CMV enhancer+promoter",
        "feature_type": "promoter",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "The full hCMV immediate-early enhancer and promoter as "
                        "one block. Use this when you want the complete element "
                        "rather than the two pieces separately.",
        "source":       "GenBank LT009450.1; identical sequence annotated "
                        "'hCMV-IE promoter and enhancer' in 10 independent "
                        "records",
        "aliases":      ["CMV", "hCMV-IE"],
        "sequence":     "GGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCC"
                        "CCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGACTTTCCA"
                        "TTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCACTTGGCAGTACATCAAGTGTA"
                        "TCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTGGCATTA"
                        "TGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCAT"
                        "CGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGA"
                        "CTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCA"
                        "AAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGG"
                        "TAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCT",
    },
    {
        "name":         "CMV promoter",
        "feature_type": "promoter",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "Human cytomegalovirus immediate-early promoter. Very "
                        "strong in most cell lines, but silences over time in "
                        "stable lines and works poorly in primary "
                        "haematopoietic and stem cells.",
        "source":       "GenBank LC897328.1; identical sequence annotated 'CMV "
                        "promoter' in 13 independent records",
        "aliases":      ["CMV", "pCMV"],
        "sequence":     "GTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGGATTT"
                        "CCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGAC"
                        "TTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGG"
                        "TGGGAGGTCTATATAAGCAGAGCT",
    },
    {
        "name":         "CMV/TetO2 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "CMV promoter with two tet operators inserted "
                        "downstream. Silent when TetR is bound, active on "
                        "tetracycline or doxycycline.",
        "source":       "GenBank MW987523.1; identical sequence annotated "
                        "'CMV-TetO2 promoter' in 3 independent records",
        "aliases":      ["TetO2", "T-REx", "dox"],
        "sequence":     "ACATTGATTATTGACTAGTTATTAATAGTAATCAATTACGGGGTCATTAGTTCATAGCCC"
                        "ATATATGGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAA"
                        "CGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGAC"
                        "TTTCCATTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCACTTGGCAGTACATCA"
                        "AGTGTATCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTG"
                        "GCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATT"
                        "AGTCATCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCG"
                        "GTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTG"
                        "GCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAAT"
                        "GGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCATCCCTATCAGTGATAGAGATCA"
                        "GATCTCCCTATCAGTGATAGAGA",
    },
    {
        "name":         "EF-1alpha promoter",
        "feature_type": "promoter",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "Human elongation factor 1-alpha promoter with intron "
                        "A. Weaker than CMV at peak but far more resistant to "
                        "silencing, so it is the usual choice for stable lines "
                        "and primary cells.",
        "source":       "GenBank MT891329.1; identical sequence annotated "
                        "'EF-1alpha promoter' in 9 independent records",
        "aliases":      ["EF1a", "EF-1a", "EEF1A1"],
        "sequence":     "GCTCCGGTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGG"
                        "GAGGGGTCGGCAATTGAACCGGTGCCTAGAGAAGGTGGCGCGGGGTAAACTGGGAAAGTG"
                        "ATGTCGTGTACTGGCTCCGCCTTTTTCCCGAGGGTGGGGGAGAACCGTATATAAGTGCAG"
                        "TAGTCGCCGTGAACGTTCTTTTTCGCAACGGGTTTGCCGCCAGAACACAGGTAAGTGCCG"
                        "TGTGTGGTTCCCGCGGGCCTGGCCTCTTTACGGGTTATGGCCCTTGCGTGCCTTGAATTA"
                        "CTTCCACGCCCCTGGCTGCAGTACGTGATTCTTGATCCCGAGCTTCGGGTTGGAAGTGGG"
                        "TGGGAGAGTTCGAGGCCTTGCGCTTAAGGAGCCCCTTCGCCTCGTGCTTGAGTTGAGGCC"
                        "TGGCCTGGGCGCTGGGGCCGCCGCGTGCGAATCTGGTGGCACCTTCGCGCCTGTCTCGCT"
                        "GCTTTCGATAAGTCTCTAGCCATTTAAAATTTTTGATGACCTGCTGCGACGCTTTTTTTC"
                        "TGGCAAGATAGTCTTGTAAATGCGGGCCAAGATCTGCACACTGGTATTTCGGTTTTTGGG"
                        "GCCGCGGGCGGCGACGGGGCCCGTGCGTCCCAGCGCACATGTTCGGCGAGGCGGGGCCTG"
                        "CGAGCGCGGCCACCGAGAATCGGACGGGGGTAGTCTCAAGCTGGCCGGCCTGCTCTGGTG"
                        "CCTGGCCTCGCGCCGCCGTGTATCGCCCCGCCCTGGGCGGCAAGGCTGGCCCGGTCGGCA"
                        "CCAGTTGCGTGAGCGGAAAGATGGCCGCTTCCCGGCCCTGCTGCAGGGAGCTCAAAATGG"
                        "AGGACGCGGCGCTCGGGAGAGCGGGCGGGTGAGTCACCCACACAAAGGAAAAGGGCCTTT"
                        "CCGTCCTCAGCCGTCGCTTCATGTGACTCCACGGAGTACCGGGCGCCGTCCAGGCACCTC"
                        "GATTAGTTCTCGAGCTTTTGGAGTACGTCGTCTTTAGGTTGGGGGGAGGGGTTTTATGCG"
                        "ATGGAGTTTCCCCACACTGAGTGGGTGGAGACTGAAGTTAGGCCAGCTTGGCACTTGATG"
                        "TAATTCTCCTTGGAATTTGCCCTTTTTGAGTTTGGATCTTGGTTCATTCTCAAGCCTCAG"
                        "ACAGTGGTTCAAAGTTTTTTTCTTCCATTTCAGGTGTCGTGA",
    },
    {
        "name":         "PGK promoter",
        "feature_type": "promoter",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "Mouse phosphoglycerate kinase 1 promoter. Modest, very "
                        "reliable, ubiquitous expression; standard for "
                        "selection cassettes in targeting constructs.",
        "source":       "GenBank MW503936.1; identical sequence annotated 'PGK "
                        "promoter' in 4 independent records",
        "aliases":      ["PGK", "hPGK"],
        "sequence":     "GGGTAGGGGAGGCGCTTTTCCCAAGGCAGTCTGGAGCATGCGCTTTAGCAGCCCCGCTGG"
                        "GCACTTGGCGCTACACAAGTGGCCTCTGGCCTCGCACACATTCCACATCCACCGGTAGGC"
                        "GCCAACCGGCTCCGTTCTTTGGTGGCCCCTTCGCGCCACCTTCTACTCCTCCCCTAGTCA"
                        "GGAAGTTCCCCCCCGCCCCGCAGCTCGCGTCGTGCAGGACGTGACAAATGGAAGTAGCAT"
                        "GTCTCACTAGGCTCGTGCAGATGGACAGCACCGCTGAGCAATGGAAGCGGGTAGGCCTTT"
                        "GGGGCAGCGGCCAATAGCAGCTTTGCTCCTTCGCTTTCTGGGCTCAGAGGCTGGGAAGGG"
                        "GTGGGTCCGGGGGCGGGCTCAGGGGCGGGCTCAGGGGCGGGGCGGGCGCCCGAAGGTCCT"
                        "CCGGAGGCCCGGCATTCTGCACGCTTCAAAAGCGCACGTCTGCCGCGCTGTTCTCCTCTT"
                        "CCTCATCTCCGGGCCTTTCG",
    },
    {
        "name":         "RSV promoter",
        "feature_type": "promoter",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "Rous sarcoma virus LTR enhancer/promoter. In "
                        "third-generation lentiviral transfer vectors it "
                        "replaces the 5' LTR U3 so transcription of the genome "
                        "no longer needs Tat.",
        "source":       "GenBank PZ267684.1; identical sequence annotated 'rous "
                        "sarcoma virus enhancer/promoter; rsv promoter' in 5 "
                        "independent records",
        "aliases":      ["RSV", "LTR"],
        "sequence":     "TGTAGTCTTATGCAATACTCTTGTAGTCTTGCAACATGGTAACGATGAGTTAGCAACATG"
                        "CCTTACAAGGAGAGAAAAAGCACCGTGCATGCCGATTGGTGGAAGTAAGGTGGTACGATC"
                        "GTGCCTTATTAGGAAGGCAACAGACGGGTCTGACATGGATTGGACGAACCACTGAATTGC"
                        "CGCATTGCAGAGATATTGTATTTAAGTGCCTAGCTCGATACATAAAC",
    },
    {
        "name":         "SV40 early promoter",
        "feature_type": "promoter",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "SV40 early promoter with its enhancer and origin. "
                        "Moderate strength; the usual driver for a mammalian "
                        "resistance marker rather than the payload.",
        "source":       "GenBank LT009450.1; identical sequence annotated 'SV40 "
                        "early promoter' in 10 independent records",
        "aliases":      ["SV40", "pSV40"],
        "sequence":     "CTGTGGAATGTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCAGGCAGGCAGAAG"
                        "TATGCAAAGCATGCATCTCAATTAGTCAGCAACCAGGTGTGGAAAGTCCCCAGGCTCCCC"
                        "AGCAGGCAGAAGTATGCAAAGCATGCATCTCAATTAGTCAGCAACCATAGTCCCGCCCCT"
                        "AACTCCGCCCATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATGGCTG"
                        "ACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCTGCC",
    },
    {
        "name":         "U6 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (mammalian)",
        "strand":       1,
        "color":        "#00B4D8",
        "description":  "Human U6 snRNA RNA-polymerase-III promoter. Drives "
                        "shRNA and sgRNA. Transcription starts at a fixed "
                        "position, so the guide should begin with G for best "
                        "initiation.",
        "source":       "GenBank PZ036137.1; identical sequence annotated 'U6 "
                        "promoter' in 5 independent records",
        "aliases":      ["U6", "pol III", "shRNA", "sgRNA"],
        "sequence":     "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAG"
                        "ATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGA"
                        "AAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCAT"
                        "ATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGA"
                        "C",
    },

    # ── Promoter (yeast) ───────────────────────────────────────────────
    {
        "name":         "GAL1 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (yeast)",
        "strand":       1,
        "color":        "#48CAE4",
        "description":  "Galactose-inducible GAL1 promoter. Off in glucose, "
                        "strongly on in galactose; the standard inducible "
                        "promoter in S. cerevisiae.",
        "source":       "GenBank MW820852.1; identical sequence annotated 'GAL1 "
                        "promoter' in 4 independent records",
        "aliases":      ["GAL1", "galactose"],
        "sequence":     "CGGATTAGAAGCCGCCGAGCGGGTGACAGCCCTCCGAAGGAAGACTCTCCTCCGTGCGTC"
                        "CTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATGTGCCTCGCGCCGCACTGCTCCGAA"
                        "CAATAAAGATTCTACAATACTAGCTTTTATGGTTATGAAGAGGAAAAATTGGCAGTAACC"
                        "TGGCCCCACAAACCTTCAAATGAACGAATCAAATTAACAACCATAGGATGATAATGCGAT"
                        "TAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAGCGAAGCGATGATTTTTGATCTATT"
                        "AACAGATATATAAATGCAAAAACTGCATAACCACTTTAACTAATACTTTCAACATTTTCG"
                        "GTTTGTATTACTTCTTATTCAAATGTAATAAAAGTATCAACAAAAAATTGTTAATATACC"
                        "TCTATACTTTAACGTCAAGGAG",
    },
    {
        "name":         "GAL1/GAL10 bidirectional promoter",
        "feature_type": "promoter",
        "category":     "Promoter (yeast)",
        "strand":       1,
        "color":        "#48CAE4",
        "description":  "The intergenic region driving GAL1 in one direction "
                        "and GAL10 in the other. Co-induces two genes from one "
                        "element.",
        "source":       "GenBank MT338523.1; identical sequence annotated "
                        "'GAL1/10' in 3 independent records",
        "aliases":      ["GAL1-10", "bidirectional"],
        "sequence":     "TTTCAAAAATTCTTACTTTTTTTTTGGATGGACGCAAAGAAGTTTAATAATCATATTACA"
                        "TGGCATTACCACCATATACATATCCATATCTAATCTTACTTATATGTTGTGGAAATGTAA"
                        "AGAGCCCCATTATCTTAGCCTAAAAAAACCTTCTCTTTGGAACTTTCAGTAATACGCTTA"
                        "ACTGCTCATTGCTATATTGAAGTACGGATTAGAAGCCGCCGAGCGGGTGACAGCCCTCCG"
                        "AAGGAAGACTCTCCTCCGTGCGTCCTCGTCTTCACCGGTCGCGTTCCTGAAACGCAGATG"
                        "TGCCTCGCGCCGCACTGCTCCGAACAATAAAGATTCTACAATACTAGCTTTTATGGTTAT"
                        "GAAGAGGAAAAATTGGCAGTAACCTGGCCCCACAAACCTTCAAATGAACGAATCAAATTA"
                        "ACAACCATAGGATGATAATGCGATTAGTTTTTTAGCCTTATTTCTGGGGTAATTAATCAG"
                        "CGAAGCGATGATTTTTGATCTATTAACAGATATATAAATGCAAAAACTGCATAACCACTT"
                        "TAACTAATACTTTCAACATTTTCGGTTTGTATTACTTCTTATTCAAATGTAATAAAAGTA"
                        "TCAACAAAAAATTGTTAATATACCTCTATACTTTAACGTCAAGGAGAAAAAACCC",
    },
    {
        "name":         "TEF1 promoter",
        "feature_type": "promoter",
        "category":     "Promoter (yeast)",
        "strand":       1,
        "color":        "#48CAE4",
        "description":  "Ashbya gossypii TEF promoter. Strong and constitutive "
                        "in yeast; the promoter half of the standard yeast "
                        "marker cassettes.",
        "source":       "GenBank PP239280.1; identical sequence annotated 'TEF "
                        "promoter; derived from A. gossypii' in 4 independent "
                        "records",
        "aliases":      ["TEF", "TEF1"],
        "sequence":     "GACATGAGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATG"
                        "TGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACA"
                        "TTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCA"
                        "GGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAA"
                        "TATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTT"
                        "GCTAGGATACAGTTCTCACATCACATCCGAACATAAACAACC",
    },

    # ── Promoter (plant) ───────────────────────────────────────────────
    {
        "name":         "CaMV 2x35S promoter (enhanced)",
        "feature_type": "promoter",
        "category":     "Promoter (plant)",
        "strand":       1,
        "color":        "#90E0EF",
        "description":  "35S promoter with a duplicated enhancer. Several-fold "
                        "stronger than the single-enhancer version.",
        "source":       "GenBank PX048278.1; identical sequence annotated 'CaMV "
                        "35S promoter (enhanced)' in 3 independent records",
        "aliases":      ["2x35S", "enhanced 35S"],
        "sequence":     "AACATGGTGGAGCACGACACTCTCGTCTACTCCAAGAATATCAAAGATACAGTCTCAGAA"
                        "GACCAAAGGGCTATTGAGACTTTTCAACAAAGGGTAATATCGGGAAACCTCCTCGGATTC"
                        "CATTGCCCAGCTATCTGTCACTTCATCAAAAGGACAGTAGAAAAGGAAGGTGGCACCTAC"
                        "AAATGCCATCATTGCGATAAAGGAAAGGCTATCGTTCAAGATGCCTCTGCCGACAGTGGT"
                        "CCCAAAGATGGACCCCCACCCACGAGGAGCATCGTGGAAAAAGAAGACGTTCCAACCACG"
                        "TCTTCAAAGCAAGTGGATTGATGTGAACATGGTGGAGCACGACACTCTCGTCTACTCCAA"
                        "GAATATCAAAGATACAGTCTCAGAAGACCAAAGGGCTATTGAGACTTTTCAACAAAGGGT"
                        "AATATCGGGAAACCTCCTCGGATTCCATTGCCCAGCTATCTGTCACTTCATCAAAAGGAC"
                        "AGTAGAAAAGGAAGGTGGCACCTACAAATGCCATCATTGCGATAAAGGAAAGGCTATCGT"
                        "TCAAGATGCCTCTGCCGACAGTGGTCCCAAAGATGGACCCCCACCCACGAGGAGCATCGT"
                        "GGAAAAAGAAGACGTTCCAACCACGTCTTCAAAGCAAGTGGATTGATGTGATATCTCCAC"
                        "TGACGTAAGGGATGACGCACAATCCCACTATCCTTCGCAAGACCCTTCCTCTATATAAGG"
                        "AAGTTCATTTCATTTGGAGAGGACACGCTGA",
    },
    {
        "name":         "CaMV 35S promoter",
        "feature_type": "promoter",
        "category":     "Promoter (plant)",
        "strand":       1,
        "color":        "#90E0EF",
        "description":  "Cauliflower mosaic virus 35S promoter. Strong and "
                        "near-constitutive in most dicots; the default promoter "
                        "of plant transformation.",
        "source":       "GenBank PZ562718.1, annotated 'camv 35s promoter'",
        "aliases":      ["35S", "CaMV"],
        "sequence":     "TGAGACTTTTCAACAAAGGGTAATTTCGGGAAACCTCCTCGGATTCCATTGCCCAGCTAT"
                        "CTGTCACTTCATCGAAAGGACAGTAGAAAAGGAAGGTGGCTCCTACAAATGCCATCATTG"
                        "CGATAAAGGAAAGGCTATCATTCAAGATGCCTCTGCCGACAGTGGTCCCAAAGATGGACC"
                        "CCCACCCACGAGGAGCATCGTGGAAAAAGAAGACGTTCCAACCACGTCTTCAAAGCAAGT"
                        "GGATTGATGTGACATCTCCACTGACGTAAGGGATGACGCACAATCCCACTATCCTTCGCA"
                        "AGACCCTTCCTCTATATAAGGAAGTTCATTTCATTTGGAGAGGACA",
    },
    {
        "name":         "NOS promoter",
        "feature_type": "promoter",
        "category":     "Promoter (plant)",
        "strand":       1,
        "color":        "#90E0EF",
        "description":  "Agrobacterium nopaline synthase promoter. Weak and "
                        "constitutive; usually drives the plant selection "
                        "marker rather than the payload.",
        "source":       "GenBank OQ753826.1; identical sequence annotated 'nos "
                        "promoter' in 10 independent records",
        "aliases":      ["nos", "pNOS"],
        "sequence":     "GATCATGAGCGGAGAATTAAGGGAGTCACGTTATGACCCCCGCCGATGACGCGGGACAAG"
                        "CCGTTTTACGTTTGGAACTGACAGAACCGCAACGTTGAAGGAGCCACTCAGCCGCGGGTT"
                        "TCTGGAGTTTAATGAGCTAAGCACATACGTCAGAAACCATTATTGCGCGTTCAAAAGTCG"
                        "CCTAAGGTCACTATCAGCTAGCAAATATTTCTTGTCAAAAATGCTCCACTGACGTTCCAT"
                        "AAATTCCCCTCGGTATCCAATTAGAGTCTCATATTCACTCTCAATCC",
    },

    # ── Operator ───────────────────────────────────────────────────────
    {
        "name":         "lac operator (lacO, symmetric)",
        "feature_type": "protein_bind",
        "category":     "Operator",
        "strand":       1,
        "color":        "#F08080",
        "description":  "Perfectly palindromic lac operator used in T7lac-style "
                        "vectors. Binds LacI more tightly than the natural O1 "
                        "operator; placing it just downstream of a T7 promoter "
                        "is what suppresses uninduced expression.",
        "source":       "GenBank PZ357631.1; identical sequence annotated 'lac "
                        "operator' in 8 independent records",
        "aliases":      ["lacO", "lac operator", "T7lac"],
        "sequence":     "AATTGTGAGCGCTCACAATT",
    },
    {
        "name":         "tet operator (tetO2)",
        "feature_type": "protein_bind",
        "category":     "Operator",
        "strand":       1,
        "color":        "#F08080",
        "description":  "TetR binding site. Two copies flanking a minimal "
                        "promoter give the tetracycline-regulated "
                        "(Tet-On/Tet-Off) systems their switch.",
        "source":       "GenBank PP125207.1 and PP480513.1; identical sequence "
                        "annotated 'tet operator' in both",
        "aliases":      ["tetO", "tetO2", "TetR site"],
        "sequence":     "TCCCTATCAGTGATAGAGA",
    },

    # ── Regulator ──────────────────────────────────────────────────────
    {
        "name":         "araC",
        "feature_type": "CDS",
        "category":     "Regulator",
        "strand":       1,
        "color":        "#DAA520",
        "description":  "AraC regulator of the arabinose operon. Represses pBAD "
                        "without arabinose and activates it with arabinose, so "
                        "a pBAD vector must carry it.",
        "source":       "GenBank PX640807.1; identical sequence annotated "
                        "'araC' in 10 independent records",
        "aliases":      ["araC", "arabinose regulator"],
        "sequence":     "ATGGCTGAAGCGCAAAATGATCCCCTGCTGCCGGGATACTCGTTTAATGCCCATCTGGTG"
                        "GCGGGTTTAACGCCGATTGAGGCCAACGGTTATCTCGATTTTTTCATCGATCGTCCGTTA"
                        "GGCATGAAAGGCTACATTCTCAACCTGACCATTCGTGGGCAAGGCGTGGTTAAGAACCAA"
                        "GGGCGTGAGTTTGTCTGTCGGCCTGGCGACATTTTACTGTTTCCTCCGGGTGAGATCCAC"
                        "CACTATGGGCGTCATCCGGAAGCGCGTGAATGGTATCACCAGTGGGTCTATTTTCGCCCG"
                        "CGTGCCTATTGGCACGAATGGTTGAACTGGCCCAGTATCTTTGCCAATACGGGCTTCTTC"
                        "CGTCCCGATGAGGCACATCAGCCGCACTTTTCAGATCTGTTTGGTCAGATCATTAATGCC"
                        "GGTCAGGGTGAAGGTCGGTATTCGGAACTGTTAGCGATTAACCTTCTGGAACAGTTGCTG"
                        "CTTCGCCGCATGGAGGCGATTAACGAGTCTCTGCATCCACCGATGGATAATCGCGTGCGT"
                        "GAAGCCTGCCAGTACATTTCGGATCATCTGGCAGACTCTAATTTCGACATCGCTTCTGTA"
                        "GCGCAACATGTGTGCTTGAGTCCGAGTCGCCTGAGCCATCTGTTTCGCCAGCAACTGGGC"
                        "ATCAGCGTGTTATCATGGCGCGAAGATCAGCGCATTTCCCAGGCGAAACTGCTCTTGTCC"
                        "ACCACACGCATGCCGATTGCTACGGTTGGTCGCAATGTGGGTTTCGACGATCAGCTCTAC"
                        "TTTAGCCGCGTATTCAAGAAATGTACCGGAGCTAGCCCATCGGAATTTCGCGCAGGATGC"
                        "GAAGAGAAAGTCAACGATGTTGCGGTGAAACTGAGCTAA",
    },
    {
        "name":         "ccdB (counter-selection)",
        "feature_type": "CDS",
        "category":     "Regulator",
        "strand":       1,
        "color":        "#DAA520",
        "description":  "CcdB DNA gyrase poison. Lethal in ordinary E. coli, so "
                        "an uncut or self-ligated vector carrying it yields no "
                        "colonies. Propagate only in a ccdB-resistant strain "
                        "such as DB3.1.",
        "source":       "GenBank KJ541667.1; identical sequence annotated "
                        "'ccdB' in 15 independent records",
        "aliases":      ["ccdB", "counter selection", "death gene"],
        "sequence":     "ATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTA"
                        "CAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGT"
                        "CTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGC"
                        "TGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTG"
                        "GCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGA"
                        "ATATAA",
    },
    {
        "name":         "lacI",
        "feature_type": "CDS",
        "category":     "Regulator",
        "strand":       1,
        "color":        "#DAA520",
        "description":  "Lac repressor. Binds lacO and blocks transcription "
                        "until IPTG or lactose relieves it. Note the natural "
                        "GTG start codon.",
        "source":       "GenBank PQ628036.1; identical sequence annotated "
                        "'lacI' in 10 independent records",
        "aliases":      ["lac repressor", "lacI"],
        "sequence":     "GTGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTT"
                        "TCCCGCGTGGTGAACCAGGCCAGCCACGTTTCTGCGAAAACGCGGGAAAAAGTGGAAGCG"
                        "GCGATGGCGGAGCTGAATTACATTCCCAACCGCGTGGCACAACAACTGGCGGGCAAACAG"
                        "TCGTTGCTGATTGGCGTTGCCACCTCCAGTCTGGCCCTGCACGCGCCGTCGCAAATTGTC"
                        "GCGGCGATTAAATCTCGCGCCGATCAACTGGGTGCCAGCGTGGTGGTGTCGATGGTAGAA"
                        "CGAAGCGGCGTCGAAGCCTGTAAAGCGGCGGTGCACAATCTTCTCGCGCAACGCGTCAGT"
                        "GGGCTGATCATTAACTATCCGCTGGATGACCAGGATGCCATTGCTGTGGAAGCTGCCTGC"
                        "ACTAATGTTCCGGCGTTATTTCTTGATGTCTCTGACCAGACACCCATCAACAGTATTATT"
                        "TTCTCCCATGAAGACGGTACGCGACTGGGCGTGGAGCATCTGGTCGCATTGGGTCACCAG"
                        "CAAATCGCGCTGTTAGCGGGCCCATTAAGTTCTGTCTCGGCGCGTCTGCGTCTGGCTGGC"
                        "TGGCATAAATATCTCACTCGCAATCAAATTCAGCCGATAGCGGAACGGGAAGGCGACTGG"
                        "AGTGCCATGTCCGGTTTTCAACAAACCATGCAAATGCTGAATGAGGGCATCGTTCCCACT"
                        "GCGATGCTGGTTGCCAACGATCAGATGGCGCTGGGCGCAATGCGCGCCATTACCGAGTCC"
                        "GGGCTGCGCGTTGGTGCGGATATCTCGGTAGTGGGATACGACGATACCGAAGACAGCTCA"
                        "TGTTATATCCCGCCGTTAACCACCATCAAACAGGATTTTCGCCTGCTGGGGCAAACCAGC"
                        "GTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCC"
                        "GTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGC"
                        "GCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAG"
                        "TGA",
    },
    {
        "name":         "rop/ROM",
        "feature_type": "CDS",
        "category":     "Regulator",
        "strand":       1,
        "color":        "#DAA520",
        "description":  "Rop (ROM) protein. Stabilises the RNA I / RNA II "
                        "kissing complex and so LOWERS ColE1 plasmid copy "
                        "number. Deleting it is half the reason pUC vectors are "
                        "high-copy.",
        "source":       "GenBank PZ357623.1; identical sequence annotated 'rop' "
                        "in 13 independent records",
        "aliases":      ["rom", "rop"],
        "sequence":     "GTGACCAAACAGGAAAAAACCGCCCTTAACATGGCCCGCTTTATCAGAAGCCAGACATTA"
                        "ACGCTTCTGGAGAAACTCAACGAGCTGGACGCGGATGAACAGGCAGACATCTGTGAATCG"
                        "CTTCACGACCACGCTGATGAGCTTTACCGCAGCTGCCTCGCGCGTTTCGGTGATGACGGT"
                        "GAAAACCTCTGA",
    },

    # ── Translation ────────────────────────────────────────────────────
    {
        "name":         "Kozak sequence",
        "feature_type": "regulatory",
        "category":     "Translation",
        "strand":       1,
        "color":        "#00FF7F",
        "description":  "Consensus Kozak context GCCACC immediately before the "
                        "ATG. The purine at -3 and the G at +4 matter most for "
                        "efficient initiation in mammalian cells. Includes the "
                        "ATG.",
        "source":       "GenBank MT891329.1; identical sequence annotated "
                        "'Kozak sequence' in 6 independent records",
        "aliases":      ["Kozak", "GCCACC"],
        "sequence":     "GCCACCATGG",
    },
    {
        "name":         "T7 gene 10 RBS",
        "feature_type": "RBS",
        "category":     "Translation",
        "strand":       1,
        "color":        "#00FF7F",
        "description":  "Ribosome binding site and leader from bacteriophage T7 "
                        "gene 10. Very strong initiation in E. coli; keep 6-8 "
                        "bases between the AAGGAGA Shine-Dalgarno and the start "
                        "codon.",
        "source":       "GenBank OR900359.1, annotated '@t7g10rbs'",
        "aliases":      ["RBS", "g10", "Shine-Dalgarno"],
        "sequence":     "TTTGTTTAACTTTAAGAAGGAGA",
    },

    # ── Terminator ─────────────────────────────────────────────────────
    {
        "name":         "ADH1 terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "S. cerevisiae ADH1 transcription terminator. Compact "
                        "and reliable at the 3' end of a yeast expression "
                        "cassette.",
        "source":       "GenBank PP239280.1; identical sequence annotated 'ADH1 "
                        "terminator; derived from S. cerevisiae' in 3 "
                        "independent records",
        "aliases":      ["ADH1", "ADH1t"],
        "sequence":     "GCGAATTTCTTATGATTTATGATTTTTATTATTAAATAAGTTATAAAAAAAATAAGTGTA"
                        "TACAAATTTTAAAGTGACTCTTAGGTTTTAAAACGAAAATTCTTATTCTTGAGTAACTCT"
                        "TTCCTGTAGGTCAGGTTGCTTTCTCAGGTATAGCATGAGGTCGCTCTTATTGACCACACC"
                        "TCTACCGGCATGCCGAGCAAATGCCTGCAAATCGCTCCCCATTTCACCCAATTGTAGATA"
                        "TGCTAACTCCAGCAATGAGTTGATGAATCTCGGTGTGTATTTTATGTCCTCAGAGGACAA"
                        "CACCTGT",
    },
    {
        "name":         "BBa_B1008 terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "Registry bidirectional terminator part, used to close "
                        "a BioBrick transcription unit. Bidirectional, so it "
                        "also blocks read-through arriving from the far side.",
        "source":       "GenBank PX640807.1; identical sequence annotated "
                        "'bidirectional-terminator_bba_b1008' in 10 independent "
                        "records",
        "aliases":      ["BioBrick", "B1008", "iGEM"],
        "sequence":     "CGCCAAAAACCCCGCCCCTGACAGGGCGGGGTTTTTCCGC",
    },
    {
        "name":         "CYC1 terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "S. cerevisiae CYC1 terminator. The usual partner for "
                        "GAL1- and TEF1-driven cassettes.",
        "source":       "GenBank LT727450.1; identical sequence annotated 'CYC1 "
                        "terminator' in 3 independent records",
        "aliases":      ["CYC1"],
        "sequence":     "CCCCTTTTCCTTTGTCGATATCATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCC"
                        "CCCCACATCCGCTCTAACCGAAAAGGAAGGAGTTAGACAACCTGAAGTCTAGGTCCCTAT"
                        "TTATTTTTTTATAGTTATGTTAGTATTAAGAACGTTATTTATATTTCAAATTTTTCTTTT"
                        "TTTTCTGTACAGACGCGTGTACGCATGTAACATTATACTGAAAACCTTGCTTGAGAAGGT"
                        "TTTGGGACGCTCGAAGGCTTTAATTTGCAAG",
    },
    {
        "name":         "lambda t0 terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "Bacteriophage lambda t0 terminator. Compact, and "
                        "orthogonal in sequence to the rrnB terminators, which "
                        "avoids repeat-driven recombination when you need two "
                        "terminators on one plasmid.",
        "source":       "GenBank OK148689.1; identical sequence annotated "
                        "'bacteriophage lambda T0' in 7 independent records",
        "aliases":      ["lambda t0", "T0"],
        "sequence":     "GACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGC"
                        "TCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAAT",
    },
    {
        "name":         "NOS terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "Nopaline synthase terminator. The standard 3' end for "
                        "a plant expression cassette.",
        "source":       "GenBank OQ753826.1; identical sequence annotated 'nos "
                        "terminator' in 10 independent records",
        "aliases":      ["nos terminator", "tNOS"],
        "sequence":     "CCCCCGAATTTCCCCGATCGTTCAAACATTTGGCAATAAAGTTTCTTAAGATTGAATCCT"
                        "GTTGCCGGTCTTGCGATGATTATCATATAATTTCTGTTGAATTACGTTAAGCATGTAATA"
                        "ATTAACATGTAATGCATGACGTTATTTATGAGATGGGTTTTTATGATTAGAGTCCCGCAA"
                        "TTATACATTTAATACGCGATAGAAAACAAAATATAGCGCGCAAACTAGGATAAATTATCG"
                        "CGCGCGGTGTCATCTATGTTACTAGA",
    },
    {
        "name":         "rrnB T1 terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "E. coli rrnB T1 rho-independent terminator. Strong and "
                        "bidirectional in practice; the standard insulator "
                        "between cassettes.",
        "source":       "GenBank OK148689.1; identical sequence annotated "
                        "'Escherichia coli rrnB T1' in 7 independent records",
        "aliases":      ["rrnB", "T1"],
        "sequence":     "AGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGT"
                        "TTGTCGGTGAACGCTCTCCTGAGTAGGACAAATCCGCC",
    },
    {
        "name":         "rrnB T1+T2 terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "Tandem rrnB T1 and T2 terminators. The strongest "
                        "common transcriptional stop in E. coli vectors; use it "
                        "to insulate a sensitive cassette.",
        "source":       "GenBank PX895169.1; identical sequence annotated 'T1T2 "
                        "bidirectional transcription terminator' in 4 "
                        "independent records",
        "aliases":      ["T1T2", "rrnB"],
        "sequence":     "GAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTG"
                        "TTTGTCGGTGAACGCTCTCCTGAGTAGGACAAATCCGCCGCTCTAGCTAAGCAGAAGGCC"
                        "ATCCTGACGGATGGCCTTTTTGCGTTTCTACAAACTCTTGTTAACTCTAGAGCTGCCTGC"
                        "CGCGTTTCGGTGATGAAGATCTTCCCGATGATTAATTAATTCAGAACGCTCGGTTGCCGC"
                        "CGGGCGTTTTTTATGCAGCAATGGCAAGAACGTTGCTCTAGA",
    },
    {
        "name":         "T7 terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "Bacteriophage T7 T-phi terminator from gene 10. Placed "
                        "downstream of the coding sequence in T7 expression "
                        "vectors to stop read-through.",
        "source":       "GenBank PV764404.1; identical sequence annotated 'T7 "
                        "terminator' in 9 independent records",
        "aliases":      ["T7 term", "Tphi"],
        "sequence":     "TAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
    },
    {
        "name":         "T7Te terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "Short T7 early terminator. 28 bp, so it fits where "
                        "rrnB will not.",
        "source":       "GenBank PZ357631.1; identical sequence annotated 'T7TE "
                        "terminator' in 8 independent records",
        "aliases":      ["T7Te", "TE"],
        "sequence":     "GGCTCACCTTCGGGTGGGCCTTTCTGCG",
    },
    {
        "name":         "U6 terminator",
        "feature_type": "terminator",
        "category":     "Terminator",
        "strand":       1,
        "color":        "#DC143C",
        "description":  "Pol III termination signal: a run of T residues "
                        "closing a U6-driven shRNA or sgRNA cassette.",
        "source":       "GenBank PZ036137.1; identical sequence annotated 'u6 "
                        "terminator' in 5 independent records",
        "aliases":      ["pol III terminator", "TTTTTT"],
        "sequence":     "GCGGACTTCGGTCCGCTTTTT",
    },

    # ── polyA signal ───────────────────────────────────────────────────
    {
        "name":         "BGH poly(A) signal",
        "feature_type": "polyA_signal",
        "category":     "polyA signal",
        "strand":       1,
        "color":        "#FF6347",
        "description":  "Bovine growth hormone polyadenylation signal. Compact "
                        "and efficient; the default 3' end in pcDNA-family "
                        "vectors.",
        "source":       "GenBank MW503936.1; identical sequence annotated 'BGH "
                        "poly(A) signal' in 6 independent records",
        "aliases":      ["BGH", "bGH polyA"],
        "sequence":     "CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCC"
                        "TGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTC"
                        "TGAGTAGGTGTCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATT"
                        "GGGAAGACAATAGCAGGCATGCTGGGGA",
    },
    {
        "name":         "CaMV poly(A) signal",
        "feature_type": "polyA_signal",
        "category":     "polyA signal",
        "strand":       1,
        "color":        "#FF6347",
        "description":  "Cauliflower mosaic virus polyadenylation signal. An "
                        "alternative plant 3' end when you already used the NOS "
                        "terminator elsewhere.",
        "source":       "GenBank PX048278.1; identical sequence annotated 'CaMV "
                        "poly(A) signal' in 11 independent records",
        "aliases":      ["35S polyA"],
        "sequence":     "TTTCTCCATAATAATGTGTGAGTAGTTCCCAGATAAGGGAATTAGGGTTCCTATAGGGTT"
                        "TCGCTCATGTGTTGAGCATATAAGAAACCCTTAGTATGTATTTGTATTTGTAAAATACTT"
                        "CTATCAATAAAATTTCTAATTCCTAAAACCAAAATCCAGTACTAAAATCCAGATC",
    },
    {
        "name":         "Rabbit beta-globin poly(A)",
        "feature_type": "polyA_signal",
        "category":     "polyA signal",
        "strand":       1,
        "color":        "#FF6347",
        "description":  "Rabbit beta-globin 3' UTR and polyadenylation signal, "
                        "the 3' end of the CAG expression cassette. Well "
                        "tolerated across cell types.",
        "source":       "GenBank LT727518.1, annotated '3' UTR of rabbit "
                        "beta-globin, incl. polyA'",
        "aliases":      ["beta-globin polyA", "rbGlob"],
        "sequence":     "GATCTTTTTCCCTCTGCCAAAAATTATGGGGACATCATGAAGCCCCTTGAGCATCTGACT"
                        "TCTGGCTAATAAAGGAAATTTATTTTCATTGC",
    },
    {
        "name":         "SV40 late poly(A) signal",
        "feature_type": "polyA_signal",
        "category":     "polyA signal",
        "strand":       1,
        "color":        "#FF6347",
        "description":  "SV40 late polyadenylation signal. Stronger than the "
                        "early signal and commonly used on a marker cassette.",
        "source":       "GenBank MZ090948.1; identical sequence annotated 'SV40 "
                        "late poly(A) signal' in 10 independent records",
        "aliases":      ["SV40 late polyA"],
        "sequence":     "CAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAA"
                        "AATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCA"
                        "ATAAACAAGTTAACAACAACAATTGCATTCATTTTATGTTTCAGGTTCAGGGGGAGGTGT"
                        "GGGAGGTTTTTTAAAGCAAGTAAAACCTCTACAAATGTGGTA",
    },
    {
        "name":         "SV40 poly(A) signal",
        "feature_type": "polyA_signal",
        "category":     "polyA signal",
        "strand":       1,
        "color":        "#FF6347",
        "description":  "SV40 early polyadenylation signal. Short, which "
                        "matters inside a size-limited viral genome.",
        "source":       "GenBank MG356850.1; identical sequence annotated 'SV40 "
                        "poly(A) signal' in 8 independent records",
        "aliases":      ["SV40 polyA"],
        "sequence":     "AACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACA"
                        "AATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCT"
                        "TA",
    },
    {
        "name":         "Synthetic poly(A) / pause site",
        "feature_type": "polyA_signal",
        "category":     "polyA signal",
        "strand":       1,
        "color":        "#FF6347",
        "description":  "Synthetic polyadenylation and transcriptional pause "
                        "element. Used to insulate a downstream cassette from "
                        "read-through transcription.",
        "source":       "GenBank MZ090948.1; identical sequence annotated "
                        "'synthetic poly(A) signal/transcriptional pause site' "
                        "in 10 independent records",
        "aliases":      ["synthetic polyA", "pause site"],
        "sequence":     "AATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTA"
                        "CTAACATACGCTCTCCATCAAAACAAAACGAAACAAAACAAACTAGCAAAATAGGCTGTC"
                        "CCCAGTGCAAGTGCAGGTGCCAGAACATTTCTCT",
    },

    # ── Origin ─────────────────────────────────────────────────────────
    {
        "name":         "2-micron ori",
        "feature_type": "rep_origin",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "Saccharomyces cerevisiae 2-micron plasmid origin. "
                        "High-copy episomal maintenance in yeast (30-100 "
                        "copies), at the price of unequal segregation.",
        "source":       "GenBank LT727380.1; identical sequence annotated '2mu "
                        "ori' in 6 independent records",
        "aliases":      ["2u", "2 micron", "yeast episomal"],
        "sequence":     "AGAAAGTATAGGAACTTCAGAGCGCTTTTGAAAACCAAAAGCGCTCTGAAGACGCACTTT"
                        "CAAAAAACCAAAAACGCACCGGACTGTAACGAGCTACTAAAATATTGCGAATACCGCTTC"
                        "CACAAACATTGCTCAAAAGTATCTCTTTGCTATATATCTCTGTGCTATATCCCTATATAA"
                        "CCTACCCATCCACCTTTCGCTCCTTGAACTTGCATCTAAACTCGACCTCTACATTTTTTA"
                        "TGTTTATCTCTAGTATTACTCTTTAGACAAAAAAATTGTAGTAAGAACTATTCATAGAGT"
                        "GAATCGAAAACAATACGAAAATGTAAACATTTCCTATACGTAGTATATAGAGACAAAATA"
                        "GAAGAAACCGTTCATAATTTTCTGACCAATGAAGAATCATCAACGCTATCACTTTCTGTT"
                        "CACAAAGTATGCGCAATCCACATCGGTATAGAATATAATCGGGGATGCCTTTATCTTGAA"
                        "AAAATGCACCCGCAGCTTCGCTAGTAATCAGTAAACGCGGGAAGTGGAGTCAGGCTTTTT"
                        "TTATGGAAGAGAAAATAGACACCAAAGTAGCCTTCTTCTAACCTTAACGGACCTACAGTG"
                        "CAAAAAGTTATCAAGAGACTGCATTATAGAGCGCACAAAGGAGAAAAAAAGTAATCTAAG"
                        "ATGCTTTGTTAGAAAAATAGCGCTCTCGGGATGCATTTTTGTAGAACAAAAAAGAAGTAT"
                        "AGATTCTTTGTTGGTAAAATAGCGCTCTCGCGTTGCATTTCTGTTCTGTAAAAATGCAGC"
                        "TCAGATTCTTTGTTTGAAAAATTAGCGCTCTCGCGTTGCATTTTTGTTTTACAAAAATGA"
                        "AGCACAGATTCTTCGTTGGTAAAATAGCGCTTTCGCGTTGCATTTCTGTTCTGTAAAAAT"
                        "GCAGCTCAGATTCTTTGTTTGAAAAATTAGCGCTCTCGCGTTGCATTTTTGTTCTACAAA"
                        "ATGAAGCACAGATGCTTCGTT",
    },
    {
        "name":         "bom (mobility region)",
        "feature_type": "misc_feature",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "Basis-of-mobility region from pBR322. A relic in most "
                        "modern vectors; it allows mobilisation when the nic "
                        "site and helper functions are supplied in trans.",
        "source":       "GenBank PZ357623.1; identical sequence annotated 'bom' "
                        "in 12 independent records",
        "aliases":      ["mobility", "nic"],
        "sequence":     "CGCAGCCATGACCCAGTCACGTAGCGATAGCGGAGTGTATACTGGCTTAACTATGCGGCA"
                        "TCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTA"
                        "AGGAGAAAATACCGCATCAGG",
    },
    {
        "name":         "CEN6/ARS4",
        "feature_type": "rep_origin",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "Yeast centromere plus autonomously replicating "
                        "sequence. Holds the plasmid at 1-2 copies per cell and "
                        "segregates it faithfully, the low-copy counterpart to "
                        "the 2-micron origin.",
        "source":       "GenBank PP239280.1; identical sequence annotated "
                        "'cen/ars derived from s. cerevisiae s288c' in 3 "
                        "independent records",
        "aliases":      ["CEN", "ARS", "centromere"],
        "sequence":     "ATCACGTGCTATAAAAATAATTATAATTTAAATTTTTTAATATAAATATATAAATTAAAA"
                        "ATAGAAAGTAAAAAAAGAAATTAAAGAAAAAATAGTTTTTGTTTTCCGAAGATGTAAAAG"
                        "ACTCTAGGGGGATCGCCAACAAATACTACCTTTTATCTTGCTCTTCCTGCTCTCAGGTAT"
                        "TAATGCCGAATTGTTTCATCTTGTCTGTGTAGAAGACCACACACGAAAATCCTGTGATTT"
                        "TACATTTTACTTATCGTTAATCGAATGTATATCTATTTAATCTGCTTTTCTTGTCTAATA"
                        "AATATATATGTAAAGTACGCTTTTTGTTGAAATTTTTTAAACCTTTGTTTATTTTTTTTT"
                        "CTTCATTCCGTAACTCTTCTACCTTCTTTATTTACTTTCTAAAATCCAAATACAAAACAT"
                        "AAAAATAAATAAACACAGAGTAAATTCCCAAATTATTCCATCATTAAAAGATACGAGGCG"
                        "CGTGTAAGTTACAGGCAAGCGATC",
    },
    {
        "name":         "f1 ori",
        "feature_type": "rep_origin",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "Filamentous phage f1 origin. With a helper phage the "
                        "vector packages single-stranded DNA of one strand, "
                        "which is what makes a phagemid a phagemid.",
        "source":       "GenBank LT009450.1; identical sequence annotated 'f1 "
                        "origin' in 10 independent records",
        "aliases":      ["f1", "phagemid", "M13 ori"],
        "sequence":     "ACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCG"
                        "CTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCA"
                        "CGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTA"
                        "GTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTAGGGTGATGGTTCACGTAGTGGGC"
                        "CATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTG"
                        "GACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTAT"
                        "AAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTA"
                        "ACGCGAATTTTAACAAAATATTAACGCTTACAATTT",
    },
    {
        "name":         "ori (pUC/pMB1, high-copy)",
        "feature_type": "rep_origin",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "pMB1/ColE1-derived origin carrying the pUC point "
                        "mutation that removes copy-number control. 500-700 "
                        "copies per cell. Incompatible with any other "
                        "ColE1-family plasmid in the same strain.",
        "source":       "GenBank MN232102.1; identical sequence annotated "
                        "'high-copy-number ColE1/pMB1/pBR322/pUC origin of "
                        "replication' in 10 independent records",
        "aliases":      ["pUC ori", "ColE1", "pMB1"],
        "sequence":     "TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACC"
                        "AGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTT"
                        "CAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTT"
                        "CAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGC"
                        "TGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAA"
                        "GGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGAC"
                        "CTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGG"
                        "GAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGA"
                        "GCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACT"
                        "TGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAA",
    },
    {
        "name":         "oriT (RP4/incP)",
        "feature_type": "oriT",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "RP4/incP origin of transfer. Lets a conjugative helper "
                        "mobilise the plasmid from donor to recipient, which is "
                        "how you move DNA into strains that will not take "
                        "electroporation.",
        "source":       "GenBank PZ357609.1; identical sequence annotated "
                        "'oriT' in 10 independent records",
        "aliases":      ["conjugation", "RP4", "incP", "mob"],
        "sequence":     "CCGGCCAGCCTCGCAGAGCAGGATTCCCGTTGAGCACCGCCAGGTGCGAATAAGGGACAG"
                        "TGAAGAAGGAACACCCGCTCGCGGGTGGGCCTACTTCACCTATCCTGCCC",
    },
    {
        "name":         "p15A ori",
        "feature_type": "rep_origin",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "p15A replication origin from pACYC184. About 10-12 "
                        "copies per cell and compatible with ColE1/pUC "
                        "plasmids, so it is the standard second-plasmid origin "
                        "for co-expression.",
        "source":       "GenBank PZ366119.1; identical sequence annotated 'p15A "
                        "origin of replication' in 4 independent records",
        "aliases":      ["pACYC", "p15A"],
        "sequence":     "TTGAGATCGTTTTGGTCTGCGCGTAATCTCTTGCTCTGAAAACGAAAAAACCGCCTTGCA"
                        "GGGCGGTTTTTCGAAGGTTCTCTGAGCTACCAACTCTTTGAACCGAGGTAACTGGCTTGG"
                        "AGGAGCGCAGTCACCAAAACTTGTCCTTTCAGTTTAGCCTTAACCGGCGCATGACTTCAA"
                        "GACTAACTCCTCTAAATCAATTACCAGTGGCTGCTGCCAGTGGTGCTTTTGCATGTCTTT"
                        "CCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGACTGAACGGGGG"
                        "GTTCGTGCATACAGTCCAGCTTGGAGCGAACTGCCTACCCGGAACTGAGTGTCAGGCGTG"
                        "GAATGAGACAAACGCGGCCATAACAGCGGAATGACACCGGTAAACCGAAAGGCAGGAACA"
                        "GGAGAGCGCACGAGGGAGCCGCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGG"
                        "TTTCGCCACCACTGATTTGAGCGTCAGATTTCGTGATGCTTGTCAGGGGGGCGGAGCCTA"
                        "TGGAAA",
    },
    {
        "name":         "R6K gamma ori",
        "feature_type": "rep_origin",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "R6K gamma origin. Replicates only in strains "
                        "expressing the pir gene, so the plasmid is a suicide "
                        "vector everywhere else. Standard for allelic exchange "
                        "and transposon delivery.",
        "source":       "GenBank PZ357609.1; identical sequence annotated 'R6K "
                        "gamma ori' in 10 independent records",
        "aliases":      ["R6K", "pir", "suicide"],
        "sequence":     "TGTCAGCCGTTAAGTGTTCCTGTGTCACTCAAAATTGCTTTGAGAGGCTCTAAGGGCTTC"
                        "TCAGTGCGTTACATCCCTGGCTTGTTGTCCACAACCGTTAAACCTTAAAAGCTTTAAAAG"
                        "CCTTATATATTCTTTTTTTTCTTATAAAACTTAAAACCTTAGAGGCTATTTAAGTTGCTG"
                        "ATTTATATTAATTTTATTGTTCAAACATGAGAGCTTAGTACGTGAAACATGAGAGCTTAG"
                        "TACGTTAGCCATGAGAGCTTAGTACGTTAGCCATGAGGGTTTAGTTCGTTAAACATGAGA"
                        "GCTTAGTACGTTAAACATGAGAGCTTAGTACGTGAAACATGAGAGCTTAGTACGTACTAT"
                        "CAACAGGTTGAACTGCTGATCTTCAGATC",
    },
    {
        "name":         "SV40 ori",
        "feature_type": "rep_origin",
        "category":     "Origin",
        "strand":       1,
        "color":        "#9370DB",
        "description":  "SV40 replication origin. Drives episomal amplification "
                        "only in cells expressing large T antigen (COS-1/COS-7, "
                        "HEK293T).",
        "source":       "GenBank LT009450.1; identical sequence annotated 'SV40 "
                        "ori' in 11 independent records",
        "aliases":      ["SV40"],
        "sequence":     "GACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGA"
                        "AGTAGTGAGGAGGCTT",
    },

    # ── Resistance ─────────────────────────────────────────────────────
    {
        "name":         "AmpR (bla)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "TEM-1 beta-lactamase. Confers ampicillin and "
                        "carbenicillin resistance in E. coli. The workhorse "
                        "marker of pUC/pBR322-family vectors. Secreted to the "
                        "periplasm, so satellite colonies appear on old plates.",
        "source":       "GenBank AF027126.1; identical sequence annotated 'bla' "
                        "in 39 independent records",
        "aliases":      ["ampicillin", "beta-lactamase", "AmpR"],
        "sequence":     "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCT"
                        "GTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCA"
                        "CGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCC"
                        "GAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCC"
                        "CGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTG"
                        "GTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTA"
                        "TGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATC"
                        "GGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTT"
                        "GATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATG"
                        "CCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCT"
                        "TCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGC"
                        "TCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCT"
                        "CGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTAC"
                        "ACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCC"
                        "TCACTGATTAAGCATTGGTAA",
    },
    {
        "name":         "BlastR (bsr)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Blasticidin S deaminase from Bacillus cereus. Confers "
                        "blasticidin resistance in mammalian cells. Short "
                        "coding sequence, useful when cargo space is tight.",
        "source":       "GenBank LT009451.1; identical sequence annotated 'bsr' "
                        "in 4 independent records",
        "aliases":      ["blasticidin", "bsr", "bsd"],
        "sequence":     "ATGGCCAAGCCTTTGTCTCAAGAAGAATCCACCCTCATTGAAAGAGCAACGGCTACAATC"
                        "AACAGCATCCCCATCTCTGAAGACTACAGCGTCGCCAGCGCAGCTCTCTCTAGCGACGGC"
                        "CGCATCTTCACTGGTGTCAATGTATATCATTTTACTGGGGGACCTTGTGCAGAACTCGTG"
                        "GTGCTGGGCACTGCTGCTGCTGCGGCAGCTGGCAACCTGACTTGTATCGTCGCGATCGGA"
                        "AATGAGAACAGGGGCATCTTGAGCCCCTGCGGACGGTGCCGACAGGTGCTTCTCGATCTG"
                        "CATCCTGGGATCAAAGCCATAGTGAAGGACAGTGATGGACAGCCGACGGCAGTTGGGATT"
                        "CGTGAATTGCTGCCCTCTGGTTATGTGTGGGAGGGCTAA",
    },
    {
        "name":         "CmR (cat)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Chloramphenicol acetyltransferase. Confers "
                        "chloramphenicol resistance. The marker on "
                        "pACYC/p15A-family vectors, so it pairs with an AmpR or "
                        "KanR plasmid in the same cell.",
        "source":       "GenBank KJ541667.1; identical sequence annotated 'cat' "
                        "in 13 independent records",
        "aliases":      ["chloramphenicol", "cat"],
        "sequence":     "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCCCAATGGCATCGTAAAGAA"
                        "CATTTTGAGGCATTTCAGTCAGTTGCTCAATGTACCTATAACCAGACCGTTCAGCTGGAT"
                        "ATTACGGCCTTTTTAAAGACCGTAAAGAAAAATAAGCACAAGTTTTATCCGGCCTTTATT"
                        "CACATTCTTGCCCGCCTGATGAATGCTCATCCGGAATTCCGTATGGCAATGAAAGACGGT"
                        "GAGCTGGTGATATGGGATAGTGTTCACCCTTGTTACACCGTTTTCCATGAGCAAACTGAA"
                        "ACGTTTTCATCGCTCTGGAGTGAATACCACGACGATTTCCGGCAGTTTCTACACATATAT"
                        "TCGCAAGATGTGGCGTGTTACGGTGAAAACCTGGCCTATTTCCCTAAAGGGTTTATTGAG"
                        "AATATGTTTTTCGTCTCAGCCAATCCCTGGGTGAGTTTCACCAGTTTTGATTTAAACGTG"
                        "GCCAATATGGACAACTTCTTCGCCCCCGTTTTCACCATGGGCAAATATTATACGCAAGGC"
                        "GACAAGGTGCTGATGCCGCTGGCGATTCAGGTTCATCATGCCGTCTGTGATGGCTTCCAT"
                        "GTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGCGGGGCGTAA",
    },
    {
        "name":         "HygR (hph)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Hygromycin B phosphotransferase. Confers hygromycin B "
                        "resistance in bacteria, plants, yeast and mammalian "
                        "cells, which makes it the common marker for plant "
                        "transformation.",
        "source":       "GenBank OQ753831.1; identical sequence annotated 'hpt' "
                        "in 5 independent records",
        "aliases":      ["hygromycin", "hph", "hpt"],
        "sequence":     "ATGAAAAAGCCTGAACTCACCGCGACGTCTGTCGAGAAGTTTCTGATCGAAAAGTTCGAC"
                        "AGCGTCTCCGACCTGATGCAGCTCTCGGAGGGCGAAGAATCTCGTGCTTTCAGCTTCGAT"
                        "GTAGGAGGGCGTGGATATGTCCTGCGGGTAAATAGCTGCGCCGATGGTTTCTACAAAGAT"
                        "CGTTATGTTTATCGGCACTTTGCATCGGCCGCGCTCCCGATTCCGGAAGTGCTTGACATT"
                        "GGGGCATTCAGCGAGAGCCTGACCTATTGCATCTCCCGCCGTGCACAGGGTGTCACGTTG"
                        "CAAGACCTGCCTGAAACCGAACTGCCCGCTGTTCTGCAGCCGGTCGCGGAGGCCATGGAT"
                        "GCGATCGCTGCGGCCGATCTTAGCCAGACGAGCGGGTTCGGCCCATTCGGACCGCAAGGA"
                        "ATCGGTCAATACACTACATGGCGTGATTTCATATGCGCGATTGCTGATCCCCATGTGTAT"
                        "CACTGGCAAACTGTGATGGACGACACCGTCAGTGCGTCCGTCGCGCAGGCTCTCGATGAG"
                        "CTGATGCTTTGGGCCGAGGACTGCCCCGAAGTCCGGCACCTCGTGCACGCGGATTTCGGC"
                        "TCCAACAATGTCCTGACGGACAATGGCCGCATAACAGCGGTCATTGACTGGAGCGAGGCG"
                        "ATGTTCGGGGATTCCCAATACGAGGTCGCCAACATCTTCTTCTGGAGGCCGTGGTTGGCT"
                        "TGTATGGAGCAGCAGACGCGCTACTTCGAGCGGAGGCATCCGGAGCTTGCAGGATCGCCG"
                        "CGGCTCCGGGCGTATATGCTCCGCATTGGTCTTGACCAACTCTATCAGAGCTTGGTTGAC"
                        "GGCAATTTCGATGATGCAGCTTGGGCGCAGGGTCGATGCGACGCAATCGTCCGATCCGGA"
                        "GCCGGGACTGTCGGGCGTACACAAATCGCCCGCAGAAGCGCGGCCGTCTGGACCGATGGC"
                        "TGTGTAGAAGTACTCGCCGATAGTGGAAACCGACGCCCCAGCACTCGTCCGAGGGCAAAG"
                        "GAATAG",
    },
    {
        "name":         "KanR (aph(3')-Ia)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Aminoglycoside 3'-phosphotransferase from Tn903. "
                        "Confers kanamycin resistance in E. coli. The kanamycin "
                        "marker of pET and pACYC-family vectors.",
        "source":       "GenBank PZ357627.1; identical sequence annotated 'aph' "
                        "in 9 independent records",
        "aliases":      ["kanamycin", "aph3", "neo"],
        "sequence":     "ATGAGCCATATTCAACGGGAAACGTCTTGCTCCAGGCCGCGATTAAATTCCAACATGGAT"
                        "GCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCAGGTGCGACAATC"
                        "TATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGC"
                        "GTTGCCAATGATGTTACAGATGAGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCT"
                        "CTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCG"
                        "ATCCCCGGGAAAACAGCATTCCAGGTATTAGAAGAATATCCTGATTCAGGTGAAAATATT"
                        "GTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAATTGTCCT"
                        "TTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTG"
                        "GTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGCCTGTTGAACAAGTCTGGAAA"
                        "GAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCA"
                        "CTTGATAACCTTATTTTTGACGAGGGGAAATTAATAGGTTGTATTGATGTTGGACGAGTC"
                        "GGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCT"
                        "CCTTCATTACAGAAACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAA"
                        "TTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAA",
    },
    {
        "name":         "NeoR/KanR (nptII)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Neomycin phosphotransferase II from Tn5. Confers "
                        "kanamycin resistance in bacteria and G418 (geneticin) "
                        "resistance in mammalian cells, which is why it doubles "
                        "as a mammalian selection marker.",
        "source":       "GenBank PV031341.1; identical sequence annotated "
                        "'aph(3')-II' in 6 independent records",
        "aliases":      ["neomycin", "G418", "nptII", "aph3-II"],
        "sequence":     "ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTC"
                        "GGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCA"
                        "GCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTG"
                        "CAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTG"
                        "CTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAG"
                        "GATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATG"
                        "CGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGC"
                        "ATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAA"
                        "GAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGAC"
                        "GGCGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAAT"
                        "GGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGAC"
                        "ATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTC"
                        "CTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTT"
                        "GACGAGTTCTTCTGA",
    },
    {
        "name":         "PuroR (pac)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Puromycin N-acetyltransferase from Streptomyces "
                        "alboniger. Confers puromycin resistance in mammalian "
                        "cells. Kills fast (2-4 days), so it is the usual "
                        "marker in lentiviral knockdown vectors.",
        "source":       "GenBank LT009450.1; identical sequence annotated 'pac' "
                        "in 6 independent records",
        "aliases":      ["puromycin", "pac"],
        "sequence":     "ATGGCCACCGAGTACAAGCCCACGGTGCGCCTCGCCACCCGCGACGACGTCCCCCGGGCC"
                        "GTACGCACCCTCGCCGCCGCGTTCGCCGACTACCCCGCCACGCGCCACACCGTCGACCCG"
                        "GACCGCCACATCGAGCGGGTCACCGAGCTGCAAGAACTCTTCCTCACGCGCGTCGGGCTC"
                        "GACATCGGCAAGGTGTGGGTCGCGGACGACGGCGCCGCGGTGGCGGTCTGGACCACGCCG"
                        "GAGAGCGTCGAAGCGGGGGCGGTGTTCGCCGAGATCGGCTCGCGCATGGCCGAGTTGAGC"
                        "GGTTCCCGGCTGGCCGCGCAGCAACAGATGGAAGGCCTCCTGGCGCCGCACCGGCCCAAG"
                        "GAGCCCGCGTGGTTCCTGGCCACCGTCGGCGTCTCGCCCGACCACCAGGGCAAGGGTCTG"
                        "GGCAGCGCCGTCGTGCTCCCCGGAGTGGAGGCGGCCGAGCGCGCTGGGGTGCCCGCCTTC"
                        "CTGGAGACCTCCGCGCCCCGCAACCTCCCCTTCTACGAGCGGCTCGGCTTCACCGTCACC"
                        "GCCGACGTCGAGGTGCCCGAAGGACCGCGCACCTGGTGCATGACCCGCAAGCCCGGTGCC"
                        "TGA",
    },
    {
        # Re-identified 2026-09-22: this sequence was catalogued as aadA
        # because public records LABEL it "aadA (kanamycin resistance)". The
        # protein is 264/264 identical to APH(3')-IIIa (aphA-3, Tn1545) — a
        # kinase with the APH catalytic loop (HGDxxxxN) and no
        # nucleotidyltransferase motif — so the "kanamycin" half of that
        # label was right and "aadA" was the error. Annotating from the old
        # entry told users a kanamycin plasmid was spectinomycin-resistant.
        "name":         "KanR (aph(3')-IIIa)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Aminoglycoside 3'-phosphotransferase type IIIa "
                        "(aphA-3, carried on the Gram-positive transposon "
                        "Tn1545). Confers kanamycin and neomycin "
                        "resistance. Public "
                        "records frequently label this gene 'aadA'; it is a "
                        "phosphotransferase, not the spectinomycin/"
                        "streptomycin adenylyltransferase.",
        "source":       "GenBank PX048277.1; identical sequence in 10+ "
                        "independent records (mislabelled 'aadA' in several); "
                        "protein 264/264 identical to APH(3')-IIIa",
        "aliases":      ["aphA-3", "aph(3')-IIIa", "aph3-III"],
        "sequence":     "ATGGCTAAAATGAGAATATCACCGGAATTGAAAAAACTGATCGAAAAATACCGCTGCGTA"
                        "AAAGATACGGAAGGAATGTCTCCTGCTAAGGTATATAAGCTGGTGGGAGAAAATGAAAAC"
                        "CTATATTTAAAAATGACGGACAGCCGGTATAAAGGGACCACCTATGATGTGGAACGGGAA"
                        "AAGGACATGATGCTATGGCTGGAAGGAAAGCTGCCTGTTCCAAAGGTCCTGCACTTTGAA"
                        "CGGCATGATGGCTGGAGCAATCTGCTCATGAGTGAGGCCGATGGCGTCCTTTGCTCGGAA"
                        "GAGTATGAAGATGAACAAAGCCCTGAAAAGATTATCGAGCTGTATGCGGAGTGCATCAGG"
                        "CTCTTTCACTCCATCGACATATCGGATTGTCCCTATACGAATAGCTTAGACAGCCGCTTA"
                        "GCCGAATTGGATTACTTACTGAATAACGATCTGGCCGATGTGGATTGCGAAAACTGGGAA"
                        "GAAGACACTCCATTTAAAGATCCGCGCGAGCTGTATGATTTTTTAAAGACGGAAAAGCCC"
                        "GAAGAGGAACTTGTCTTTTCCCACGGCGACCTGGGAGACAGCAACATCTTTGTGAAAGAT"
                        "GGCAAAGTAAGTGGCTTTATTGATCTTGGGAGAAGCGGCAGGGCGGACAAGTGGTATGAC"
                        "ATTGCCTTCTGCGTCCGGTCGATCAGGGAGGATATCGGGGAAGAACAGTATGTCGAGCTA"
                        "TTTTTTGACTTACTGGGGATCAAGCCTGATTGGGAGAAAATAAAATATTATATTTTACTG"
                        "GATGAATTGTTTTAG",
    },
    {
        # Curated 2026-09-22 by the catalogue's own method (multi-accession
        # consensus over NCBI vector records), with the protein checked
        # against AadA1 — the label-count winner for "aadA" in the same
        # corpus was the APH(3')-IIIa above, which is how the old entry went
        # wrong.
        "name":         "SmR/SpecR (aadA)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Aminoglycoside 3''-adenylyltransferase (aadA1, "
                        "ANT(3'')-Ia). Confers streptomycin AND spectinomycin "
                        "resistance; does not confer kanamycin resistance.",
        "source":       "GenBank M60473.1; identical sequence in 6 independent "
                        "records; protein 263/263 identical to AadA1",
        "aliases":      ["spectinomycin", "streptomycin", "aadA", "aadA1"],
        "sequence":     "ATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATC"
                        "GAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGC"
                        "GGCCTGAAGCCACACAGTGATATTGATTTGCTGGTTACGGTGACCGTAAGGCTTGATGAA"
                        "ACAACGCGGCGAGCTTTGATCAACGACCTTTTGGAAACTTCGGCTTCCCCTGGAGAGAGC"
                        "GAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGT"
                        "TATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGT"
                        "ATCTTCGAGCCAGCCACGATCGACATTGATCTGGCTATCTTGCTGACAAAAGCAAGAGAA"
                        "CATAGCGTTGCCTTGGTAGGTCCAGCGGCGGAGGAACTCTTTGATCCGGTTCCTGAACAG"
                        "GATCTATTTGAGGCGCTAAATGAAACCTTAACGCTATGGAACTCGCCGCCCGACTGGGCT"
                        "GGCGATGAGCGAAATGTAGTGCTTACGTTGTCCCGCATTTGGTACAGCGCAGTAACCGGC"
                        "AAAATCGCGCCGAAGGATGTCGCTGCCGACTGGGCAATGGAGCGCCTGCCGGCCCAGTAT"
                        "CAGCCCGTCATACTTGAAGCTAGACAGGCTTATCTTGGACAAGAAGAAGATCGCTTGGCC"
                        "TCGCGCGCAGATCAGTTGGAAGAATTTGTCCACTACGTGAAAGGCGAGATCACCAAGGTA"
                        "GTCGGCAAATAA",
    },
    {
        "name":         "TetR/TetA (tetA)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Tetracycline efflux pump from pBR322. Confers "
                        "tetracycline resistance. Selection is leaky at low "
                        "drug concentration and the pump is toxic when "
                        "overexpressed.",
        "source":       "GenBank PZ357623.1; identical sequence annotated 'tcr' "
                        "in 6 independent records",
        "aliases":      ["tetracycline", "tetA"],
        "sequence":     "ATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTGTAGGC"
                        "ATAGGCTTGGTTATGCCGGTACTGCCGGGCCTCTTGCGGGATATCGTCCATTCCGACAGC"
                        "ATCGCCAGTCACTATGGCGTGCTGCTAGCGCTATATGCGTTGATGCAATTTCTATGCGCA"
                        "CCCGTTCTCGGAGCACTGTCCGACCGCTTTGGCCGCCGCCCAGTCCTGCTCGCTTCGCTA"
                        "CTTGGAGCCACTATCGACTACGCGATCATGGCGACCACACCCGTCCTGTGGATCCTCTAC"
                        "GCCGGACGCATCGTGGCCGGCATCACCGGCGCCACAGGTGCGGTTGCTGGCGCCTATATC"
                        "GCCGACATCACCGATGGGGAAGATCGGGCTCGCCACTTCGGGCTCATGAGCGCTTGTTTC"
                        "GGCGTGGGTATGGTGGCAGGCCCCGTGGCCGGGGGACTGTTGGGCGCCATCTCCTTGCAT"
                        "GCACCATTCCTTGCGGCGGCGGTGCTCAACGGCCTCAACCTACTACTGGGCTGCTTCCTA"
                        "ATGCAGGAGTCGCATAAGGGAGAGCGTCGACCGATGCCCTTGAGAGCCTTCAACCCAGTC"
                        "AGCTCCTTCCGGTGGGCGCGGGGCATGACTATCGTCGCCGCACTTATGACTGTCTTCTTT"
                        "ATCATGCAACTCGTAGGACAGGTGCCGGCAGCGCTCTGGGTCATTTTCGGCGAGGACCGC"
                        "TTTCGCTGGAGCGCGACGATGATCGGCCTGTCGCTTGCGGTATTCGGAATCTTGCACGCC"
                        "CTCGCTCAAGCCTTCGTCACTGGTCCCGCCACCAAACGTTTCGGCGAGAAGCAGGCCATT"
                        "ATCGCCGGCATGGCGGCCGACGCGCTGGGCTACGTCTTGCTGGCGTTCGCGACGCGAGGC"
                        "TGGATGGCCTTCCCCATTATGATTCTTCTCGCTTCCGGCGGCATCGGGATGCCCGCGTTG"
                        "CAGGCCATGCTGTCCAGGCAGGTAGATGACGACCATCAGGGACAGCTTCAAGGATCGCTC"
                        "GCGGCTCTTACCAGCCTAACTTCGATCATTGGACCGCTGATCGTCACGGCGATTTATGCC"
                        "GCCTCGGCGAGCACATGGAACGGGTTGGCATGGATTGTAGGCGCCGCCCTATACCTTGTC"
                        "TGCCTCCCCGCGTTGCGTCGCGGTGCATGGAGCCGGGCCACCTCGACCTGA",
    },
    {
        "name":         "TRP1",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Phosphoribosylanthranilate isomerase. Selects for "
                        "tryptophan prototrophy in a trp1 yeast strain.",
        "source":       "GenBank PQ010743.1; identical sequence annotated "
                        "'trp1' in 3 independent records",
        "aliases":      ["TRP1", "auxotrophy"],
        "sequence":     "ATGTCTGTTATTAATTTCACAGGTAGTTCTGGTCCATTGGTGAAAGTTTGCGGCTTGCAG"
                        "AGCACAGAGGCCGCAGAATGTGCTCTAGATTCCGATGCTGACTTGCTGGGTATTATATGT"
                        "GTGCCCAATAGAAAGAGAACAATTGACCCGGTTATTGCAAGGAAAATTTCAAGTCTTGTA"
                        "AAAGCATATAAAAATAGTTCAGGCACTCCGAAATACTTGGTTGGCGTGTTTCGTAATCAA"
                        "CCTAAGGAGGATGTTTTGGCTCTGGTCAATGATTACGGCATTGATATCGTCCAACTGCAT"
                        "GGAGATGAGTCGTGGCAAGAATACCAAGAGTTCCTCGGTTTGCCAGTTATTAAAAGACTC"
                        "GTATTTCCAAAAGACTGCAACATACTACTCAGTGCAGCTTCACAGAAACCTCATTCGTTT"
                        "ATTCCCTTGTTTGATTCAGAAGCAGGTGGGACAGGTGAACTTTTGGATTGGAACTCGATT"
                        "TCTGACTGGGTTGGAAGGCAAGAGAGCCCCGAAAGCTTACATTTTATGTTAGCTGGTGGA"
                        "CTGACGCCAGAAAATGTTGGTGATGCGCTTAGATTAAATGGCGTTATTGGTGTTGATGTA"
                        "AGCGGAGGTGTGGAGACAAATGGTGTAAAAGACTCTAACAAAATAGCAAATTTCGTCAAA"
                        "AATGCTAAGAAATAG",
    },
    {
        "name":         "URA3",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Orotidine-5'-phosphate decarboxylase. Selects for "
                        "uracil prototrophy, and counter-selects on 5-FOA, "
                        "which is what makes it recyclable.",
        "source":       "GenBank KJ922019.1; identical sequence annotated "
                        "'ura3' in 5 independent records",
        "aliases":      ["URA3", "5-FOA", "auxotrophy"],
        "sequence":     "ATGTCGAAAGCTACATATAAGGAACGTGCTGCTACTCATCCTAGTCCTGTTGCTGCCAAG"
                        "CTATTTAATATCATGCACGAAAAGCAAACAAACTTGTGTGCTTCATTGGATGTTCGTACC"
                        "ACCAAGGAATTACTGGAGTTAGTTGAAGCATTAGGTCCCAAAATTTGTTTACTAAAAACA"
                        "CATGTGGATATCTTGACTGATTTTTCCATGGAGGGCACAGTTAAGCCGCTAAAGGCATTA"
                        "TCCGCCAAGTACAATTTTTTACTCTTCGAAGACAGAAAATTTGCTGACATTGGTAATACA"
                        "GTCAAATTGCAGTACTCTGCGGGTGTATACAGAATAGCAGAATGGGCAGACATTACGAAT"
                        "GCACACGGTGTGGTGGGCCCAGGTATTGTTAGCGGTTTGAAGCAGGCGGCAGAAGAAGTA"
                        "ACAAAGGAACCTAGAGGCCTTTTGATGTTAGCAGAATTGTCATGCAAGGGCTCCCTATCT"
                        "ACTGGAGAATATACTAAGGGTACTGTTGACATTGCGAAGAGCGACAAAGATTTTGTTATC"
                        "GGCTTTATTGCTCAAAGAGACATGGGTGGAAGAGATGAAGGTTACGATTGGTTGATTATG"
                        "ACACCCGGTGTGGGTTTAGATGACAAGGGAGACGCATTGGGTCAACAGTATAGAACCGTG"
                        "GATGATGTGGTCTCTACAGGATCTGACATTATTATTGTTGGAAGAGGACTATTTGCAAAG"
                        "GGAAGGGATGCTAAGGTAGAGGGTGAACGTTACAGAAAAGCAGGCTGGGAAGCATATTTG"
                        "AGAAGATGCGGCCAGCAAAACTAA",
    },
    {
        "name":         "ZeoR (Sh ble)",
        "feature_type": "CDS",
        "category":     "Resistance",
        "strand":       1,
        "color":        "#FF8C42",
        "description":  "Bleomycin-binding protein from Streptoalloteichus "
                        "hindustanus. Confers zeocin resistance in bacteria, "
                        "yeast, plants and mammalian cells. Stoichiometric, not "
                        "enzymatic, so expression level must stay high.",
        "source":       "GenBank X52869.1, annotated 'unnamed protein product; "
                        "ble protein (aa 1-124)'",
        "aliases":      ["zeocin", "bleomycin", "ble", "Sh ble"],
        "sequence":     "ATGGCCAAGTTGACCAGTGCCGTTCCGGTGCTCACCGCGCGCGACGTCGCCGGAGCGGTC"
                        "GAGTTCTGGACCGACCGGCTCGGGTTCTCCCGGGACTTCGTGGAGGACGACTTCGCCGGT"
                        "GTGGTCCGGGACGACGTGACCCTGTTCATCAGCGCGGTCCAGGACCAGGTGGTGCCGGAC"
                        "AACACCCTGGCCTGGGTGTGGGTGCGCGGCCTGGACGAGCTGTACGCCGAGTGGTCGGAG"
                        "GTCGTGTCCACGAACTTCCGGGACGCCTCCGGGCCGGCCATGACCGAGATCGGCGAGCAG"
                        "CCGTGGGGGCGGGAGTTCGCCCTGCGCGACCCGGCCGGCAACTGCGTGCACTTCGTGGCC"
                        "GAGGAGCAGGACTGA",
    },

    # ── Reporter ───────────────────────────────────────────────────────
    {
        "name":         "EGFP",
        "feature_type": "CDS",
        "category":     "Reporter",
        "strand":       1,
        "color":        "#FFD700",
        "description":  "Enhanced green fluorescent protein. Excitation 488 nm, "
                        "emission 507 nm. Human-codon-optimised with the "
                        "F64L/S65T folding mutations.",
        "source":       "GenBank LC884826.1; identical sequence annotated "
                        "'EGFP' in 7 independent records",
        "aliases":      ["GFP", "EGFP", "green"],
        "sequence":     "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC"
                        "GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTAC"
                        "GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACC"
                        "CTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG"
                        "CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTC"
                        "TTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTG"
                        "GTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCAC"
                        "AAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC"
                        "GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCC"
                        "GACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCAC"
                        "TACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTC"
                        "CTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA",
    },
    {
        "name":         "Firefly luciferase (luc+)",
        "feature_type": "CDS",
        "category":     "Reporter",
        "strand":       1,
        "color":        "#FFD700",
        "description":  "Photinus pyralis luciferase, codon-improved. Wide "
                        "dynamic range, so it is the usual reporter for "
                        "promoter-activity assays. Needs luciferin.",
        "source":       "GenBank AF027126.1; identical sequence annotated "
                        "'luc+' in 4 independent records",
        "aliases":      ["luciferase", "luc", "firefly"],
        "sequence":     "ATGGAAGACGCCAAAAACATAAAGAAAGGCCCGGCGCCATTCTATCCGCTGGAAGATGGA"
                        "ACCGCTGGAGAGCAACTGCATAAGGCTATGAAGAGATACGCCCTGGTTCCTGGAACAATT"
                        "GCTTTTACAGATGCACATATCGAGGTGGACATCACTTACGCTGAGTACTTCGAAATGTCC"
                        "GTTCGGTTGGCAGAAGCTATGAAACGATATGGGCTGAATACAAATCACAGAATCGTCGTA"
                        "TGCAGTGAAAACTCTCTTCAATTCTTTATGCCGGTGTTGGGCGCGTTATTTATCGGAGTT"
                        "GCAGTTGCGCCCGCGAACGACATTTATAATGAACGTGAATTGCTCAACAGTATGGGCATT"
                        "TCGCAGCCTACCGTGGTGTTCGTTTCCAAAAAGGGGTTGCAAAAAATTTTGAACGTGCAA"
                        "AAAAAGCTCCCAATCATCCAAAAAATTATTATCATGGATTCTAAAACGGATTACCAGGGA"
                        "TTTCAGTCGATGTACACGTTCGTCACATCTCATCTACCTCCCGGTTTTAATGAATACGAT"
                        "TTTGTGCCAGAGTCCTTCGATAGGGACAAGACAATTGCACTGATCATGAACTCCTCTGGA"
                        "TCTACTGGTCTGCCTAAAGGTGTCGCTCTGCCTCATAGAACTGCCTGCGTGAGATTCTCG"
                        "CATGCCAGAGATCCTATTTTTGGCAATCAAATCATTCCGGATACTGCGATTTTAAGTGTT"
                        "GTTCCATTCCATCACGGTTTTGGAATGTTTACTACACTCGGATATTTGATATGTGGATTT"
                        "CGAGTCGTCTTAATGTATAGATTTGAAGAAGAGCTGTTTCTGAGGAGCCTTCAGGATTAC"
                        "AAGATTCAAAGTGCGCTGCTGGTGCCAACCCTATTCTCCTTCTTCGCCAAAAGCACTCTG"
                        "ATTGACAAATACGATTTATCTAATTTACACGAAATTGCTTCTGGTGGCGCTCCCCTCTCT"
                        "AAGGAAGTCGGGGAAGCGGTTGCCAAGAGGTTCCATCTGCCAGGTATCAGGCAAGGATAT"
                        "GGGCTCACTGAGACTACATCAGCTATTCTGATTACACCCGAGGGGGATGATAAACCGGGC"
                        "GCGGTCGGTAAAGTTGTTCCATTTTTTGAAGCGAAGGTTGTGGATCTGGATACCGGGAAA"
                        "ACGCTGGGCGTTAATCAAAGAGGCGAACTGTGTGTGAGAGGTCCTATGATTATGTCCGGT"
                        "TATGTAAACAATCCGGAAGCGACCAACGCCTTGATTGACAAGGATGGATGGCTACATTCT"
                        "GGAGACATAGCTTACTGGGACGAAGACGAACACTTCTTCATCGTTGACCGCCTGAAGTCT"
                        "CTGATTAAGTACAAAGGCTATCAGGTGGCTCCCGCTGAATTGGAATCCATCTTGCTCCAA"
                        "CACCCCAACATCTTCGACGCAGGTGTCGCAGGTCTTCCCGACGATGACGCCGGTGAACTT"
                        "CCCGCCGCCGTTGTTGTTTTGGAGCACGGAAAGACGATGACGGAAAAAGAGATCGTGGAT"
                        "TACGTCGCCAGTCAAGTAACAACCGCGAAAAAGTTGCGCGGAGGAGTTGTGTTTGTGGAC"
                        "GAAGTACCGAAAGGTCTTACCGGAAAACTCGACGCAAGAAAAATCAGAGAGATCCTCATA"
                        "AAGGCCAAGAAGGGCGGAAAGATCGCCGTGTAA",
    },
    {
        "name":         "lacZ-alpha",
        "feature_type": "CDS",
        "category":     "Reporter",
        "strand":       1,
        "color":        "#FFD700",
        "description":  "LacZ alpha peptide with the pUC polylinker inside it. "
                        "Complements a lacZ-delta-M15 host to give blue "
                        "colonies; an insert in the polylinker breaks it, so "
                        "recombinants stay white.",
        "source":       "GenBank OK148689.1; identical sequence annotated "
                        "'lacZalpha' in 6 independent records",
        "aliases":      ["lacZ", "blue white", "alpha complementation"],
        "sequence":     "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTA"
                        "CCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTT"
                        "ACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAG"
                        "GCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATG"
                        "CGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGGTGCACTCTCAGT"
                        "ACAATCTGCTCTGATGCCGCATAG",
    },
    {
        "name":         "mCherry",
        "feature_type": "CDS",
        "category":     "Reporter",
        "strand":       1,
        "color":        "#FFD700",
        "description":  "Monomeric red fluorescent protein. Excitation 587 nm, "
                        "emission 610 nm. Photostable and genuinely monomeric, "
                        "so it behaves in fusions.",
        "source":       "GenBank MW987534.1; identical sequence annotated "
                        "'mCherry' in 7 independent records",
        "aliases":      ["mCherry", "red", "RFP"],
        "sequence":     "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAG"
                        "GTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGC"
                        "CGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCC"
                        "TTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCAC"
                        "CCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGC"
                        "GTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGAC"
                        "GGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTA"
                        "ATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGC"
                        "GCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCT"
                        "GAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTC"
                        "AACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAA"
                        "CGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA",
    },
    {
        "name":         "mRFP1",
        "feature_type": "CDS",
        "category":     "Reporter",
        "strand":       1,
        "color":        "#FFD700",
        "description":  "First monomeric red fluorescent protein, the ancestor "
                        "of the mFruit series. Dimmer and less photostable than "
                        "mCherry.",
        "source":       "GenBank OQ753831.1; identical sequence annotated "
                        "'mRFP' in 5 independent records",
        "aliases":      ["mRFP", "RFP"],
        "sequence":     "ATGGCTTCCTCCGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGT"
                        "TCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGT"
                        "ACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATC"
                        "CTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCGTACGTTAAACACCCGGCTGACATCCCG"
                        "GACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAA"
                        "GACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTAC"
                        "AAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACC"
                        "ATGGGTTGGGAAGCGTCCACCGAACGTATGTACCCGGAAGACGGTGCTCTGAAAGGTGAA"
                        "ATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACC"
                        "TACATGGCTAAAAAACCGGTTCAGCTGCCGGGTGCTTACAAAACCGACATCAAACTGGAC"
                        "ATCACCTCCCACAACGAAGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGT"
                        "CACTCCACCGGTGCTTAA",
    },
    {
        "name":         "sfGFP",
        "feature_type": "CDS",
        "category":     "Reporter",
        "strand":       1,
        "color":        "#FFD700",
        "description":  "Superfolder GFP. Folds and fluoresces even when fused "
                        "to poorly folding partners, which makes it the safer "
                        "choice for a fusion tag.",
        "source":       "GenBank PZ366119.1; identical sequence annotated "
                        "'sfGFP' in 4 independent records",
        "aliases":      ["superfolder", "sfGFP"],
        "sequence":     "ATGCGTAAAGGCGAAGAGCTGTTCACTGGTGTCGTCCCTATTCTGGTGGAACTGGATGGT"
                        "GATGTCAACGGTCATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGT"
                        "AAACTGACGCTGAAGTTCATCTGTACTACTGGTAAACTGCCGGTACCTTGGCCGACTCTG"
                        "GTAACGACGCTGACTTATGGTGTTCAGTGCTTTGCTCGTTATCCGGACCATATGAAGCAG"
                        "CATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCACGATTTCCTTT"
                        "AAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTA"
                        "AACCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAATATCCTGGGCCATAAG"
                        "CTGGAATACAATTTTAACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGC"
                        "ATTAAAGCGAATTTTAAAATTCGCCACAACGTGGAGGATGGCAGCGTGCAGCTGGCTGAT"
                        "CACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCTGCTGCCAGACAATCACTAT"
                        "CTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCATATGGTTCTG"
                        "CTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAATAA",
    },

    # ── Tag ────────────────────────────────────────────────────────────
    {
        "name":         "3xFLAG tag",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "Three tandem FLAG epitopes. Far more sensitive by "
                        "Western blot than a single copy.",
        "source":       "GenBank OQ427855.1, annotated '3flag'",
        "aliases":      ["3xFLAG", "triple FLAG"],
        "sequence":     "GACTACAAAGACCATGACGGTGATTATAAAGATCATGACATCGATTACAAGGATGACGAT"
                        "GACAAG",
    },
    {
        "name":         "6xHis tag",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "Hexahistidine tag for IMAC (Ni-NTA or Co) "
                        "purification. No stop codon, so it can sit at either "
                        "end of a coding sequence. Keep it in frame.",
        "source":       "GenBank LC823121.1 and LC853333.1; identical sequence "
                        "annotated '6xHis affinity tag' in both",
        "aliases":      ["His6", "His tag", "IMAC", "Ni-NTA"],
        "sequence":     "CATCACCATCACCATCAC",
    },
    {
        "name":         "8xHis tag",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "Octahistidine tag. Binds the resin more tightly than "
                        "6xHis, which helps when the protein elutes early or "
                        "the lysate is dirty.",
        "source":       "GenBank KJ541667.1; identical sequence annotated "
                        "'8xHis tag' in 3 independent records",
        "aliases":      ["His8"],
        "sequence":     "CACCACCATCACCACCATCACCAC",
    },
    {
        "name":         "FLAG tag",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "DYKDDDDK epitope. Immunoprecipitates with anti-FLAG "
                        "resin and elutes under native conditions with FLAG "
                        "peptide.",
        "source":       "GenBank LT009451.1 and LT009452.1; identical sequence "
                        "annotated 'FLAG epitope' in both",
        "aliases":      ["FLAG", "DYKDDDDK"],
        "sequence":     "GACTACAAGGATGACGATGACAAA",
    },
    {
        "name":         "GST tag",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "Schistosoma japonicum glutathione S-transferase, a 26 "
                        "kDa solubility and affinity tag eluted with reduced "
                        "glutathione. Dimeric, so cleave it off before any "
                        "oligomeric-state work. This is the GST open reading "
                        "frame only: the vector thrombin site, FLAG tag, "
                        "polylinker and stop codon that the source record folds "
                        "into the same CDS are trimmed away so it fuses in "
                        "frame.",
        "source":       "GenBank U67875.1, bases 1-660 of the annotated "
                        "feature; annotated 'glutathione S-transferase'",
        "aliases":      ["GST", "glutathione"],
        "sequence":     "ATGTCCCCTATACTAGGTTATTGGAAAATTAAGGGCCTTGTGCAACCCACTCGACTTCTT"
                        "TTGGAATATCTTGAAGAAAAATATGAAGAGCATTTGTATGAGCGCGATGAAGGTGATAAA"
                        "TGGCGAAACAAAAAGTTTGAATTGGGTTTGGAGTTTCCCAATCTTCCTTATTATATTGAT"
                        "GGTGATGTTAAATTAACACAGTCTATGGCCATCATACGTTATATAGCTGACAAGCACAAC"
                        "ATGTTGGGTGGTTGTCCAAAAGAGCGTGCAGAGATTTCAATGCTTGAAGGAGCGGTTTTG"
                        "GATATTAGATACGGTGTTTCGAGAATTGCATATAGTAAAGACTTTGAAACTCTCAAAGTT"
                        "GATTTTCTTAGCAAGCTACCTGAAATGCTGAAAATGTTCGAAGATCGTTTATGTCATAAA"
                        "ACATATTTAAATGGTGATCATGTAACCCATCCTGACTTCATGTTGTATGACGCTCTTGAT"
                        "GTTGTTTTATACATGGACCCAATGTGCCTGGATGCGTTCCCAAAATTAGTTTGTTTTAAA"
                        "AAACGTATTGAAGCTATCCCACAAATTGATAAGTACTTGAAATCCAGCAAGTATATAGCA"
                        "TGGCCTTTGCAGGGCTGGCAAGCCACGTTTGGTGGTGGCGACCATCCTCCAAAATCGGAT",
    },
    {
        "name":         "HA tag",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "YPYDVPDYA epitope from influenza haemagglutinin. "
                        "Small, well tolerated and backed by excellent "
                        "antibodies.",
        "source":       "GenBank LT009453.1 and LT009454.1; identical sequence "
                        "annotated 'HA epitope' in both",
        "aliases":      ["HA", "YPYDVPDYA"],
        "sequence":     "TACCCCTACGACGTGCCCGACTACGCC",
    },
    {
        "name":         "Myc tag",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "EQKLISEEDL epitope from human c-Myc. Commonly paired "
                        "with a His tag at the opposite terminus for two-step "
                        "detection.",
        "source":       "GenBank MT338523.1; identical sequence annotated "
                        "'myc-tag' in 3 independent records",
        "aliases":      ["myc", "c-myc", "EQKLISEEDL"],
        "sequence":     "GAACAAAAGTTAATTTCTGAAGAGGACTTG",
    },
    {
        "name":         "Strep-tag II",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "WSHPQFEK tag. Binds Strep-Tactin and elutes gently "
                        "with desthiobiotin, which suits intact complexes.",
        "source":       "GenBank PV031344.1, annotated 'region: streptag II'",
        "aliases":      ["Strep", "StrepII", "WSHPQFEK"],
        "sequence":     "TGGAGCCACCCGCAGTTCGAAAAG",
    },
    {
        "name":         "V5 tag",
        "feature_type": "misc_feature",
        "category":     "Tag",
        "strand":       1,
        "color":        "#DA70D6",
        "description":  "GKPIPNPLLGLDST epitope from simian virus 5. Low "
                        "background in mammalian cells.",
        "source":       "GenBank PZ267685.1 and PZ267686.1; identical sequence "
                        "annotated 'epitope tag from simian virus 5; V5 tag' in "
                        "both",
        "aliases":      ["V5", "GKPIPNPLLGLDST"],
        "sequence":     "GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG",
    },

    # ── Protease site ──────────────────────────────────────────────────
    {
        "name":         "PreScission/3C site",
        "feature_type": "misc_feature",
        "category":     "Protease site",
        "strand":       1,
        "color":        "#FF69B4",
        "description":  "LEVLFQ/GP site for HRV 3C (PreScission) protease. Cuts "
                        "efficiently at 4 C, which suits unstable proteins.",
        "source":       "GenBank KJ541667.1 and KJ541669.2; identical sequence "
                        "annotated 'prescission cleavage site' in both",
        "aliases":      ["PreScission", "3C", "LEVLFQGP"],
        "sequence":     "CTGGAGGTGCTCTTCCAGGGTCCG",
    },
    {
        "name":         "TEV protease site",
        "feature_type": "misc_feature",
        "category":     "Protease site",
        "strand":       1,
        "color":        "#FF69B4",
        "description":  "ENLYFQ/S recognition site for tobacco etch virus "
                        "protease. TEV cuts between the Q and the following "
                        "residue, so one extra residue stays on the product. "
                        "Keep it in frame between the tag and the protein.",
        "source":       "GenBank JN792439.2, bases 46-66 of the annotated "
                        "feature; annotated 'n-ternimal tag; contains his6 and "
                        "tev-protease cleavage site'",
        "aliases":      ["TEV", "ENLYFQ"],
        "sequence":     "GAGAACCTGTACTTCCAATCC",
    },
    {
        "name":         "Thrombin site",
        "feature_type": "misc_feature",
        "category":     "Protease site",
        "strand":       1,
        "color":        "#FF69B4",
        "description":  "LVPR/GS thrombin cleavage site. The tag-removal site "
                        "in pGEX-family vectors. Thrombin can nick at secondary "
                        "sites, so watch the digest.",
        "source":       "GenBank U13853.1, annotated 'encodes thrombin "
                        "recognition site'",
        "aliases":      ["thrombin", "LVPRGS"],
        "sequence":     "CTGGTTCCGCGTGGATCC",
    },

    # ── Linker / 2A ────────────────────────────────────────────────────
    {
        "name":         "(GGGGS)x3 linker",
        "feature_type": "misc_feature",
        "category":     "Linker / 2A",
        "strand":       1,
        "color":        "#C792EA",
        "description":  "Three GGGGS repeats. Flexible, protease-resistant "
                        "spacer that keeps two fused domains from interfering "
                        "with each other.",
        "source":       "GenBank PP297162.1, annotated 'flexible linker "
                        "(ggggs)x3 zf-opt'",
        "aliases":      ["GS linker", "GGGGS"],
        "sequence":     "GGAGGAGGAGGATCTGGAGGAGGAGGATCTGGAGGAGGAGGATCT",
    },
    {
        "name":         "P2A peptide",
        "feature_type": "misc_feature",
        "category":     "Linker / 2A",
        "strand":       1,
        "color":        "#C792EA",
        "description":  "Porcine teschovirus-1 2A peptide. Usually the most "
                        "efficient 2A; add a GSG spacer in front of it to raise "
                        "skipping further.",
        "source":       "GenBank MW079274.1, annotated 'p2a site'",
        "aliases":      ["P2A", "2A"],
        "sequence":     "GCGACCAACTTTAGCCTGCTGAAACAGGCGGGCGATGTGGAAGAAAACCCAGGACCG",
    },
    {
        "name":         "T2A peptide",
        "feature_type": "misc_feature",
        "category":     "Linker / 2A",
        "strand":       1,
        "color":        "#C792EA",
        "description":  "Thosea asigna virus 2A self-cleaving peptide. Ribosome "
                        "skipping gives two separate proteins from one open "
                        "reading frame. Leaves ~18 residues on the upstream "
                        "protein and one proline on the downstream one.",
        "source":       "GenBank MW503936.1; identical sequence annotated "
                        "'similar to T2A' in 4 independent records",
        "aliases":      ["T2A", "2A", "self-cleaving"],
        "sequence":     "GAGGGCAGAGGAAGTCTGCTAACATGCGGTGACGTCGAGGAGAATCCTGGACCT",
    },

    # ── Localisation ───────────────────────────────────────────────────
    {
        "name":         "SV40 NLS",
        "feature_type": "misc_feature",
        "category":     "Localisation",
        "strand":       1,
        "color":        "#ADFF2F",
        "description":  "PKKKRKV nuclear localisation signal from SV40 large T "
                        "antigen. The default NLS for driving a fusion protein "
                        "into the nucleus.",
        "source":       "GenBank MW503936.1; identical sequence annotated "
                        "'similar to SV40 NLS' in 4 independent records",
        "aliases":      ["NLS", "PKKKRKV"],
        "sequence":     "CCCAAGAAGAAGAGGAAGGTG",
    },

    # ── Recombination ──────────────────────────────────────────────────
    {
        "name":         "attB (BxbI)",
        "feature_type": "misc_recomb",
        "category":     "Recombination",
        "strand":       1,
        "color":        "#48D1CC",
        "description":  "Bacteriophage BxbI attB site. Recombines with attP "
                        "irreversibly, which is what makes BxbI landing pads "
                        "stable.",
        "source":       "GenBank MW503936.1; identical sequence annotated 'attB "
                        "(BxbI)' in 4 independent records",
        "aliases":      ["attB", "BxbI", "integrase"],
        "sequence":     "GGCCGGCTTGTCGACGACGGCGGTCTCCGTCGTCAGGATCATCCGG",
    },
    {
        "name":         "attR1 (Gateway)",
        "feature_type": "misc_recomb",
        "category":     "Recombination",
        "strand":       1,
        "color":        "#48D1CC",
        "description":  "attR1 recombination site. With attR2 it flanks the "
                        "ccdB cassette of a Gateway destination vector; LR "
                        "clonase swaps the insert in.",
        "source":       "GenBank LT009450.1; identical sequence annotated "
                        "'Gateway recombination site attR1' in 10 independent "
                        "records",
        "aliases":      ["attR1", "Gateway", "LR"],
        "sequence":     "CATAGTGACTGGATATGTTGTGTTTTACAGTATTATGTAGTCTGTTTTTTATGCAAAATC"
                        "TAATTTAATATATTGATATTTATATCATTTTACGTTTCTCGTTCAGCTTTTTTGTACAAA"
                        "CTTGT",
    },
    {
        "name":         "attR2 (Gateway)",
        "feature_type": "misc_recomb",
        "category":     "Recombination",
        "strand":       1,
        "color":        "#48D1CC",
        "description":  "attR2 recombination site, the downstream partner of "
                        "attR1 in a Gateway destination vector.",
        "source":       "GenBank LT009450.1; identical sequence annotated "
                        "'Gateway recombination site attR2' in 10 independent "
                        "records",
        "aliases":      ["attR2", "Gateway"],
        "sequence":     "CATAGTGACTGGATATGTTGTGTTTTACAGTATTATGTAGTCTGTTTTTTATGCAAAATC"
                        "TAATTTAATATATTGATATTTATATCATTTTACGTTTCTCGTTCAGCTTTCTTGTACAAA"
                        "GTGGT",
    },
    {
        "name":         "FRT",
        "feature_type": "misc_recomb",
        "category":     "Recombination",
        "strand":       1,
        "color":        "#48D1CC",
        "description":  "Flp recombinase target. Same excise-or-invert logic as "
                        "loxP, and orthogonal to it, so Cre and Flp can be used "
                        "in the same genome.",
        "source":       "GenBank LC884821.1; identical sequence annotated 'FRT' "
                        "in 5 independent records",
        "aliases":      ["FRT", "Flp", "FLP"],
        "sequence":     "GAAGTTCCTATTCCGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC",
    },
    {
        "name":         "lox2272",
        "feature_type": "misc_recomb",
        "category":     "Recombination",
        "strand":       1,
        "color":        "#48D1CC",
        "description":  "loxP variant with a mutated spacer. It recombines with "
                        "itself but not with loxP, which is what lets two "
                        "independent Cre events share one construct.",
        "source":       "GenBank MN044709.1 and MN044710.1; identical sequence "
                        "annotated 'lox2272 site' in both",
        "aliases":      ["lox2272", "Cre"],
        "sequence":     "ATAACTTCGTATAAAGTATCCTATACGAAGTTAT",
    },
    {
        "name":         "loxP",
        "feature_type": "misc_recomb",
        "category":     "Recombination",
        "strand":       1,
        "color":        "#48D1CC",
        "description":  "Cre recombinase target: two 13 bp inverted repeats "
                        "around an 8 bp asymmetric spacer. Two loxP sites in "
                        "the same orientation excise the intervening DNA; in "
                        "opposite orientation they invert it.",
        "source":       "GenBank JQ394985.1; identical sequence annotated "
                        "'loxP' in 6 independent records",
        "aliases":      ["loxP", "Cre", "lox"],
        "sequence":     "ATAACTTCGTATAGCATACATTATACGAAGTTAT",
    },
    {
        "name":         "T-DNA left border",
        "feature_type": "repeat_region",
        "category":     "Recombination",
        "strand":       1,
        "color":        "#48D1CC",
        "description":  "Left border repeat of the Agrobacterium T-DNA. Marks "
                        "the end of the transferred segment.",
        "source":       "GenBank PX927304.1, annotated 'left border repeat from "
                        "nopaline C58 T-DNA'",
        "aliases":      ["LB", "T-DNA"],
        "sequence":     "TGGCAGGATATATTGTGGTGTAAAC",
    },
    {
        "name":         "T-DNA right border",
        "feature_type": "repeat_region",
        "category":     "Recombination",
        "strand":       1,
        "color":        "#48D1CC",
        "description":  "Right border repeat of the Agrobacterium T-DNA. "
                        "Transfer starts here, so everything you want in the "
                        "plant must sit between the right and left borders.",
        "source":       "GenBank PX048277.1; identical sequence annotated 'RB "
                        "T-DNA repeat' in 12 independent records",
        "aliases":      ["RB", "T-DNA"],
        "sequence":     "TGACAGGATATATTGGCGGGTAAAC",
    },

    # ── Viral element ──────────────────────────────────────────────────
    {
        "name":         "Chimeric intron",
        "feature_type": "intron",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "Synthetic intron built from a human beta-globin donor "
                        "and an immunoglobulin acceptor. An intron in the 5' "
                        "UTR raises steady-state mRNA for many transgenes.",
        "source":       "GenBank PZ285976.1, annotated 'chimeric intron'",
        "aliases":      ["intron", "hybrid intron"],
        "sequence":     "GTAAGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCTTGTCGAGA"
                        "CAGAGAAGACTCTTGCGTTTCTGATAGGCACCTATTGGTCTTACTGACATCCACTTTGCC"
                        "TTTCTCTCCACAG",
    },
    {
        "name":         "HIV-1 cPPT/CTS",
        "feature_type": "misc_feature",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "Central polypurine tract and central termination "
                        "sequence. Improves nuclear import and raises "
                        "lentiviral titre and transduction.",
        "source":       "GenBank LT009450.1; identical sequence annotated "
                        "'HIV-1 cPPT' in 10 independent records",
        "aliases":      ["cPPT", "CTS"],
        "sequence":     "AAAAGAAAAGGGGGGA",
    },
    {
        "name":         "HIV-1 psi packaging signal",
        "feature_type": "misc_feature",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "HIV-1 packaging signal. Required in cis on the "
                        "transfer vector for the genome to be packaged into a "
                        "lentiviral particle.",
        "source":       "GenBank LT009450.1; identical sequence annotated "
                        "'HIV-1 psi packaging signal' in 10 independent records",
        "aliases":      ["psi", "packaging"],
        "sequence":     "TGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAG",
    },
    {
        "name":         "HIV-1 RRE",
        "feature_type": "misc_feature",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "Rev response element. Rev binds it to export the "
                        "unspliced viral RNA from the nucleus, which a "
                        "lentiviral transfer vector needs.",
        "source":       "GenBank LT009450.1; identical sequence annotated "
                        "'HIV-1 Rev response element' in 10 independent records",
        "aliases":      ["RRE", "Rev"],
        "sequence":     "AGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCGTCAAT"
                        "GACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTT"
                        "GCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCA"
                        "GCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCT",
    },
    {
        "name":         "IRES (EMCV)",
        "feature_type": "misc_feature",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "Encephalomyocarditis virus internal ribosome entry "
                        "site. Lets a second cistron translate from the same "
                        "mRNA, though typically at lower level than the first.",
        "source":       "GenBank JQ394985.1 and JQ394986.1; identical sequence "
                        "annotated 'IRES' in both",
        "aliases":      ["IRES", "EMCV"],
        "sequence":     "AACGTTACTGGCCGAAGCCGCTTGGAATAAGGCCGGTGTGCGTTTGTCTATATGTTATTT"
                        "TCCACCATATTGCCGTCTTTTGGCAATGTGAGGGCCCGGAAACCTGGCCCTGTCTTCTTG"
                        "ACGAGCATTCCTAGGGGTCTTTCCCCTCTCGCCAAAGGAATGCAAGGTCTGTTGAATGTC"
                        "GTGAAGGAAGCAGTTCCTCTGGAAGCTTCTTGAAGACAAACAACGTCTGTAGCGACCCTT"
                        "TGCAGGCAGCGGAACCCCCCACCTGGCGACAGGTGCCTCTGCGGCCAAAAGCCACGTGTA"
                        "TAAGATACACCTGCAAAGGCGGCACAACCCCAGTGCCACGTTGTGAGTTGGATAGTTGTG"
                        "GAAAGAGTCAAATGGCTCTCCTCAAGCGTATTCAACAAGGGGCTGAAGGATGCCCAGAAG"
                        "GTACCCCATTGTATGGGATCTGATCTGGGGCCTCGGTGCACATGCTTTACATGTGTTTAG"
                        "TCGAGGTTAAAAAAACGTCTAGGCCCCCCGAACCACGGGGACGTGGTTTTCCTTTGAAAA"
                        "ACACGATGATAA",
    },
    {
        "name":         "Lentiviral 3' LTR (dU3)",
        "feature_type": "LTR",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "Self-inactivating 3' LTR with the U3 enhancer deleted. "
                        "After reverse transcription the deletion is copied to "
                        "the 5' LTR, which is what makes the provirus "
                        "self-inactivating.",
        "source":       "GenBank LT009450.1; identical sequence annotated "
                        "'RSV/HIV-1 dU3/3' LTR' in 10 independent records",
        "aliases":      ["3' LTR", "SIN", "dU3"],
        "sequence":     "CTGGAAGGGCTAATTCACTCCCAACGAAGACAAGATCTGCTTTTTGCTTGTACTGGGTCT"
                        "CTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTT"
                        "AAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGAC"
                        "TCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA",
    },
    {
        "name":         "Lentiviral 5' LTR (RSV/HIV-1)",
        "feature_type": "LTR",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "Chimeric RSV/HIV-1 5' long terminal repeat from a "
                        "third-generation transfer vector. Tat-independent "
                        "transcription of the vector genome.",
        "source":       "GenBank LT009450.1; identical sequence annotated "
                        "'RSV/HIV-1 5' LTR' in 10 independent records",
        "aliases":      ["5' LTR", "LTR"],
        "sequence":     "GCACCGTGCATGCCGATTGGTGGAAGTAAGGTGGTACGATCGTGCCTTATTAGGAAGGCA"
                        "ACAGACGGGTCTGACATGGATTGGACGAACCACTGAATTGCCGCATTGCAGAGATATTGT"
                        "ATTTAAGTGCCTAGCTCGATACATAAACGGGTCTCTCTGGTTAGACCAGATCTGAGCCTG"
                        "GGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGT"
                        "GCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACC"
                        "CTTTTAGTCAGTGTGGAAAATCTCTAGCA",
    },
    {
        "name":         "WPRE",
        "feature_type": "misc_feature",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "Woodchuck hepatitis virus post-transcriptional "
                        "regulatory element. Placed between the stop codon and "
                        "the polyA signal, it raises transgene expression "
                        "several-fold.",
        "source":       "GenBank PZ285978.1 and PZ285979.1; identical sequence "
                        "annotated 'WPRE' in both",
        "aliases":      ["WPRE", "woodchuck"],
        "sequence":     "AATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCT"
                        "CCTTTTACGCTATGTGGATACGCTGCTTTAATGCCTTTGTATCATGCTATTGCTTCCCGT"
                        "ATGGCTTTCATTTTCTCCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTATGAGGAGTTG"
                        "TGGCCCGTTGTCAGGCAACGTGGCGTGGTGTGCACTGTGTTTGCTGACGCAACCCCCACT"
                        "GGTTGGGGCATTGCCACCACCTGTCAGCTCCTTTCCGGGACTTTCGCTTTCCCCCTCCCT"
                        "ATTGCCACGGCGGAACTCATCGCCGCCTGCCTTGCCCGCTGCTGGACAGGGGCTCGGCTG"
                        "TTGGGCACTGACAATTCCGTGGTGTTGTCGGGGAAATCATCGTCCTTTCCTTGGCTGCTC"
                        "GCCTGTGTTGCCACCTGGATTCTGCGCGGGACGCCCTTCTGCTACGTCCCTTCGGCCCTC"
                        "AATCCAGCGGACCTTCCTTCCCGCGGCCTGCTGCCGGCTCTGCGGCCTCTTCCGCGTCTT"
                        "CGCCTTCGCCCTCAGACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGC",
    },
    {
        "name":         "WPRE3 (minimal)",
        "feature_type": "misc_feature",
        "category":     "Viral element",
        "strand":       1,
        "color":        "#8B4513",
        "description":  "Minimised WPRE. Keeps most of the expression boost in "
                        "247 bp, which matters when the viral genome is near "
                        "its packaging limit.",
        "source":       "GenBank PZ036137.1; identical sequence annotated "
                        "'WPRE3' in 5 independent records",
        "aliases":      ["WPRE3"],
        "sequence":     "ATAATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTG"
                        "CTCCTTTTACGCTATGTGGATACGCTGCTTTAATGCCTTTGTATCATGCTATTGCTTCCC"
                        "GTATGGCTTTCATTTTCTCCTCCTTGTATAAATCCTGGTTAGTTCTTGCCACGGCGGAAC"
                        "TCATCGCCGCCTGCCTTGCCCGCTGCTGGACAGGGGCTCGGCTGTTGGGCACTGACAATT"
                        "CCGTGG",
    },

    # ── CRISPR ─────────────────────────────────────────────────────────
    {
        "name":         "sgRNA scaffold (SpCas9)",
        "feature_type": "misc_RNA",
        "category":     "CRISPR",
        "strand":       1,
        "color":        "#7B68EE",
        "description":  "Single-guide RNA constant region for S. pyogenes Cas9. "
                        "Place the 20 nt spacer immediately upstream; this part "
                        "carries the tracrRNA-derived structure Cas9 binds.",
        "source":       "GenBank LC906470.1 and LC906471.1; identical sequence "
                        "annotated 'sgRNA_scaffold' in both",
        "aliases":      ["sgRNA", "scaffold", "tracrRNA", "Cas9"],
        "sequence":     "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGT"
                        "GGCACCGAGTCGGTGC",
    },

    # ── Cloning site ───────────────────────────────────────────────────
    {
        "name":         "I-SceI site",
        "feature_type": "misc_feature",
        "category":     "Cloning site",
        "strand":       1,
        "color":        "#20B2AA",
        "description":  "18 bp homing endonuclease recognition site. Long "
                        "enough to be absent from most genomes, so cutting it "
                        "makes a single defined double-strand break.",
        "source":       "GenBank OQ753826.1; identical sequence annotated "
                        "'I-SceI' in 10 independent records",
        "aliases":      ["I-SceI", "meganuclease"],
        "sequence":     "TAGGGATAACAGGGTAAT",
    },
    {
        "name":         "M13 forward primer site",
        "feature_type": "primer_bind",
        "category":     "Cloning site",
        "strand":       1,
        "color":        "#20B2AA",
        "description":  "Binding site for the M13/pUC forward sequencing "
                        "primer, just upstream of the polylinker in pUC- and "
                        "pBluescript-family vectors.",
        "source":       "GenBank AY497507.1 and AY745747.1; identical sequence "
                        "annotated 'puc/m13 forward sequencing primer' in both",
        "aliases":      ["M13F", "M13 forward", "sequencing"],
        "sequence":     "CGCCAGGGTTTTCCCAGTCACGAC",
    },
    {
        "name":         "M13 reverse primer site",
        "feature_type": "primer_bind",
        "category":     "Cloning site",
        "strand":       1,
        "color":        "#20B2AA",
        "description":  "Binding site for the M13/pUC reverse sequencing "
                        "primer, just downstream of the polylinker.",
        "source":       "GenBank AY497507.1 and AY745747.1; identical sequence "
                        "annotated 'puc/m13 reverse sequencing primer' in both",
        "aliases":      ["M13R", "M13 reverse", "sequencing"],
        "sequence":     "TCACACAGGAAACAGCTATGAC",
    },
    {
        "name":         "pUC19 MCS (EcoRI-HindIII)",
        "feature_type": "misc_feature",
        "category":     "Cloning site",
        "strand":       1,
        "color":        "#20B2AA",
        "description":  "The pUC19 polylinker, EcoRI through HindIII. Thirteen "
                        "unique sites inside lacZ-alpha, which is what couples "
                        "cloning to blue/white screening.",
        "source":       "GenBank OK148689.1; identical sequence annotated "
                        "'pUC19 multiple cloning site (EcoRI to HindIII)' in 6 "
                        "independent records",
        "aliases":      ["MCS", "polylinker", "pUC19"],
        "sequence":     "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT",
    },
]


def _shadow_key(entry: dict) -> tuple[str, str]:
    """The ``(name, feature_type)`` identity a user entry shadows a preset by.

    Both halves are forced to ``str``. `features.json` is a user-editable file
    and `_load_features` only checks that a row is a dict — so a hand-edit (or
    a partially-written file recovered from a backup) can leave a LIST or DICT
    in ``name``. Building the set key from those raw values raised
    ``TypeError: unhashable type`` and took out every merged view at once: the
    Synthesis side panel, the preset browser's "already in your library"
    marks, and the whole annotate-from-presets scan. One bad row must cost the
    user that row, not the feature.

    A stringified corrupt key simply won't match any preset, which is the
    right outcome — the row stays visible so the user can see and fix it.
    """
    if not isinstance(entry, dict):
        return ("", "")
    return (str(entry.get("name", "")), str(entry.get("feature_type", "")))


def _coerce_preset_strand(value) -> int:
    """Clamp ``value`` to a strand in ``{-1, 0, 1, 2}``, defaulting to +1.

    Never raises. The shipped catalogue is all +1, but this helper also runs
    over dicts that reached it from a UI prefill or an agent payload, and a
    bare ``int("oops")`` there would surface as a traceback in a modal
    callback instead of a usable entry.
    """
    try:
        n = int(value)
    except (TypeError, ValueError):
        return 1
    return n if n in (-1, 0, 1, 2) else 1


def _preset_features() -> list[dict]:
    """Return a deep copy of the preset catalogue.

    Deep-copied on every read for the same reason ``_load_features`` is
    ([PIT-17]): callers buffer and mutate entries (the browser marks rows,
    the import path rewrites qualifiers), and a shallow copy would let those
    edits leak back into the module-level constant for the life of the
    process.
    """
    return copy.deepcopy(_FEATURE_PRESETS)


def _preset_categories() -> tuple[str, ...]:
    """Category names in display order (only those actually used)."""
    used = {e.get("category") for e in _FEATURE_PRESETS}
    return tuple(c for c in _PRESET_CATEGORIES if c in used)


def _find_preset(name: str, feature_type: "str | None" = None) -> "dict | None":
    """Look a preset up by ``name`` (exact, case-insensitive), optionally
    pinned to a ``feature_type``. Returns a deep copy, or None."""
    if not isinstance(name, str) or not name.strip():
        return None
    want = name.strip().lower()
    for e in _FEATURE_PRESETS:
        if str(e.get("name", "")).lower() != want:
            continue
        if feature_type and e.get("feature_type") != feature_type:
            continue
        return copy.deepcopy(e)
    return None


def _preset_matches(preset: dict, needle: str) -> bool:
    """Case-insensitive substring search across name, feature type,
    category, description and aliases. Empty needle matches everything."""
    q = (needle or "").strip().lower()
    if not q:
        return True
    if not isinstance(preset, dict):
        return False
    hay = [
        str(preset.get("name", "")),
        str(preset.get("feature_type", "")),
        str(preset.get("category", "")),
        str(preset.get("description", "")),
    ]
    aliases = preset.get("aliases")
    if isinstance(aliases, list):
        hay.extend(str(a) for a in aliases)
    return any(q in h.lower() for h in hay)


def _preset_to_library_entry(preset: dict) -> dict:
    """Convert a preset into an ordinary feature-library entry.

    Drops the preset-only fields (``category`` / ``aliases`` / ``source``)
    and folds the provenance into ``qualifiers['note']`` so it survives into
    a GenBank export and the user can still see where the sequence came
    from. The result is shaped exactly like an ``AddFeatureModal`` save, so
    it can go straight through ``_save_features``.
    """
    if not isinstance(preset, dict):
        raise TypeError("preset must be a dict")
    name = str(preset.get("name", "") or "").strip()
    seq = str(preset.get("sequence", "") or "").upper()
    source = str(preset.get("source", "") or "").strip()
    note = "SpliceCraft preset"
    if source:
        note += " — " + source
    quals: dict[str, list[str]] = {"label": [name]}
    if note:
        quals["note"] = [note]
    return {
        "name":         name,
        "feature_type": str(preset.get("feature_type", "") or "misc_feature"),
        "strand":       _coerce_preset_strand(preset.get("strand", 1)),
        "color":        str(preset.get("color", "") or ""),
        "sequence":     seq,
        "qualifiers":   quals,
        "description":  str(preset.get("description", "") or ""),
    }


def _merge_presets_with_library(user_entries: "list[dict] | None") -> list[dict]:
    """Return a READ-ONLY browse/scan view: the user's own feature library
    followed by every preset the user has not already shadowed.

    A user entry always wins: a preset whose ``(name, feature_type)`` matches
    a user entry is dropped from the view entirely, so editing a preset into
    your own library replaces it rather than doubling it. User entries keep
    their own order and come first; presets follow in catalogue order, each
    flagged ``preset=True`` so the UI can label the row and refuse in-place
    edits.

    The returned list is safe to mutate (both halves are deep copies) and is
    NEVER passed to ``_save_features`` — presets live in code, not in the
    user's data directory.
    """
    out: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for e in (user_entries or []):
        if not isinstance(e, dict):
            continue
        out.append(copy.deepcopy(e))
        seen.add(_shadow_key(e))
    for p in _FEATURE_PRESETS:
        if _shadow_key(p) in seen:
            continue
        entry = copy.deepcopy(p)
        entry["preset"] = True
        out.append(entry)
    return out
