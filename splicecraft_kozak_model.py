"""Vendored plant translation-initiation matrix. GENERATED -- do not hand-edit.

Regenerate with `python3 scripts/calibrate_kozak.py --emit
splicecraft_kozak_model.py`. Every number in KOZAK_MODEL_META is held out;
that script's docstring records what each one controls for.
"""

KOZAK_MODEL_META = {
    "atg_offset": 12,
    "auc_by_record_class": {
        "NM": 0.786,
        "XM": 0.751
    },
    "auc_leave_one_species_out": {
        "ath": 0.758,
        "gma": 0.747,
        "hvu": 0.805,
        "nta": 0.744,
        "osa": 0.785,
        "sly": 0.776,
        "zma": 0.746
    },
    "auc_per_species": {
        "ath": 0.76,
        "gma": 0.706,
        "hvu": 0.83,
        "nta": 0.721,
        "osa": 0.775,
        "sly": 0.804,
        "zma": 0.753
    },
    "auc_upstream_only": 0.715,
    "auc_vs_cds_atg": 0.718,
    "auc_vs_utr_uatg": 0.766,
    "auc_within_rule_bucket": {
        "0": 0.624,
        "1": 0.683,
        "2": 0.658
    },
    "decoy": "the 5'UTR upstream ATGs of the same transcripts -- the population a uORF screen actually has to rank",
    "delta_auc_ci95": [
        0.035,
        0.074
    ],
    "delta_auc_over_rule": 0.054,
    "down": 4,
    "fetched": "2026-10-31",
    "n_cds_atg_test": 17319,
    "n_genes_test": 623,
    "n_genes_train": 1860,
    "n_utr_uatg_test": 1391,
    "position_bits": [
        0.6,
        0.38,
        0.48,
        0.74,
        0.66,
        0.56,
        0.83,
        0.95,
        1.0,
        1.86,
        1.81,
        1.26,
        0.0,
        0.0,
        0.0,
        2.48
    ],
    "real_start_consensus": "AAAAAAGCAACAATGG",
    "refused": {
        "efficiency_rate": "no initiation RATE or protein-yield number is produced. This ranks ATG contexts against each other; it is not a trained predictor of how much protein an ATG makes.",
        "partial_windows": "an ATG with fewer than 12 nt upstream is NOT scored. A shorter window sums fewer log-odds terms and is not comparable to the calibrated distribution.",
        "window_past_plus4": "-12..+9 scored 0.791 on the dev split vs 0.760 for the shipped window, but the gain is the CODING/UTR base-composition shift: a real start's downstream half is coding and a 5'UTR uATG's is not. That discriminates starts without saying anything about initiation, so the window stops at +4."
    },
    "rule_auc_vs_cds_atg": 0.663,
    "rule_auc_vs_utr_uatg": 0.712,
    "source": "2526 genes from 7 plant species (4 dicot, 3 monocot), RefSeq mRNAs with an exact, complete-ORF, ATG-verified CDS start and an annotated 5'UTR; one transcript per gene",
    "tiers": {
        "adequate": {
            "cds_atgs_at_or_above": 0.876,
            "cut": -2.6207,
            "real_starts_at_or_above": 0.95,
            "utr_uatgs_at_or_above": 0.764
        },
        "strong": {
            "cds_atgs_at_or_above": 0.193,
            "cut": 1.0858,
            "real_starts_at_or_above": 0.501,
            "utr_uatgs_at_or_above": 0.151
        }
    },
    "uatg_coverage": 0.958,
    "up": 12,
    "upstream_at_fraction": {
        "real_starts": 0.518,
        "utr_uatgs": 0.593
    },
    "window": "-12..+4"
}

# Per-position log2-odds, real plant start codons vs the 5'UTR upstream ATGs
# they compete with. Window -12..+4 with the A of the ATG at offset 12.
#
# **This matrix reads GC-rich upstream and that is CORRECT.** It is a log-odds
# ratio, not a frequency table: real starts ARE A-rich upstream
# (`AAAAAAGCAACAATGG` is the most frequent base per position), but the uATGs they
# are being discriminated FROM are even more AT-rich (0.518 vs
# 0.593 A+T over the upstream window), so the ratio tilts the other
# way. Read `real_start_consensus` in the metadata for the biology; read this
# matrix only as a discriminator.
INITIATION_PWM = [
    {"A": -0.248227, "C": +0.345055, "G": +0.309767, "T": -0.250773},
    {"A": +0.077853, "C": +0.137740, "G": +0.084110, "T": -0.246460},
    {"A": -0.011718, "C": +0.174961, "G": +0.195215, "T": -0.289228},
    {"A": -0.172164, "C": +0.276043, "G": +0.382317, "T": -0.361674},
    {"A": -0.078878, "C": +0.329561, "G": +0.181585, "T": -0.334685},
    {"A": -0.046991, "C": +0.138154, "G": +0.268129, "T": -0.287719},
    {"A": -0.272027, "C": +0.109622, "G": +0.532519, "T": -0.298910},
    {"A": -0.105939, "C": +0.543463, "G": +0.034472, "T": -0.407739},
    {"A": +0.321469, "C": +0.222200, "G": +0.022922, "T": -0.680547},
    {"A": +0.456813, "C": -0.502626, "G": +0.703983, "T": -1.154819},
    {"A": +0.140327, "C": +0.898801, "G": -0.537359, "T": -0.908907},
    {"A": +0.277177, "C": +0.189133, "G": +0.132209, "T": -0.980779},
    {"A": +0.000000, "C": +0.000000, "G": +0.000000, "T": +0.000000},
    {"A": +0.000000, "C": +0.000000, "G": +0.000000, "T": +0.000000},
    {"A": +0.000000, "C": +0.000000, "G": +0.000000, "T": +0.000000},
    {"A": -0.633980, "C": -0.984149, "G": +1.308848, "T": -1.175553}
]

# Score quantiles of 623 HELD-OUT real plant start codons, 0th to 100th
# percentile in steps of 5. Turns a score into "stronger than N% of real plant
# starts", which is the only calibrated statement this model can make.
REAL_START_QUANTILES = [
    -4.7444,  #   0th percentile
    -2.6207,  #   5th percentile
    -1.5822,  #  10th percentile
    -1.0157,  #  15th percentile
    -0.6362,  #  20th percentile
    -0.2610,  #  25th percentile
    -0.0972,  #  30th percentile
    +0.2790,  #  35th percentile
    +0.6603,  #  40th percentile
    +0.8666,  #  45th percentile
    +1.0858,  #  50th percentile
    +1.2818,  #  55th percentile
    +1.5801,  #  60th percentile
    +1.8665,  #  65th percentile
    +2.1247,  #  70th percentile
    +2.4413,  #  75th percentile
    +2.7999,  #  80th percentile
    +3.0740,  #  85th percentile
    +3.5828,  #  90th percentile
    +4.1681,  #  95th percentile
    +5.9707,  # 100th percentile
]
