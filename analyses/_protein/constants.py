"""Amino-acid alphabets, residue masses, scales, pKa sets and the DIWV table."""

from __future__ import annotations



# ============================================================
# PROTEIN ANALYSIS
# ============================================================


# ============================================================
# CONSTANTS
# ============================================================

STANDARD_AMINO_ACIDS = set(
    "ACDEFGHIKLMNPQRSTVWY"
)

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

# IUPAC protein codes beyond the 20 standard residues (core.alphabets).
# U (selenocysteine) and O (pyrrolysine) are defined residues; B, Z, J and
# X name a set of residues, and "?" marks a position that translation could
# not resolve. All keep their position in position-dependent analyses.
NONSTANDARD_RESIDUES = set("UO")
AMBIGUOUS_RESIDUE_OPTIONS = {"B": "DN", "Z": "EQ", "J": "IL"}
UNRESOLVED_RESIDUES = set("BZJX?")
STOP_SYMBOL = "*"


RNA_IUPAC = {
    "A": ("A",),
    "U": ("U",),
    "C": ("C",),
    "G": ("G",),

    "R": ("A", "G"),
    "Y": ("C", "U"),
    "S": ("G", "C"),
    "W": ("A", "U"),
    "K": ("G", "U"),
    "M": ("A", "C"),

    "B": ("C", "G", "U"),
    "D": ("A", "G", "U"),
    "H": ("A", "C", "U"),
    "V": ("A", "C", "G"),
    "N": ("A", "U", "C", "G"),
}


AMINO_ACID_RESIDUE_MASSES = {
    # Average residue masses in Da (free amino-acid mass - H2O).
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
    # Selenocysteine and pyrrolysine: Biopython IUPACData.protein_weights
    # (free amino acid, average) minus WATER_MASS.
    "U": 150.0380,
    "O": 237.2982,
}

# One H2O is retained for the complete peptide chain.
WATER_MASS = 18.01524


HYDROPHOBIC_RESIDUES = set("AVILMFWY")
POLAR_RESIDUES = set("STNQCY")
POSITIVE_RESIDUES = set("KRH")
NEGATIVE_RESIDUES = set("DE")
AROMATIC_RESIDUES = set("FWY")
ALIPHATIC_RESIDUES = set("AVIL")


# ============================================================
# KYTE-DOOLITTLE HYDROPATHY
# ============================================================

KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}


# ============================================================
# FLEXIBILITY SCALE
# ============================================================

FLEXIBILITY_SCALE = {
    "A": 0.984,
    "C": 0.906,
    "E": 1.094,
    "D": 1.068,
    "G": 1.031,
    "F": 0.915,
    "I": 0.927,
    "H": 0.950,
    "K": 1.102,
    "M": 0.952,
    "L": 0.935,
    "N": 1.048,
    "Q": 1.037,
    "P": 1.049,
    "S": 1.046,
    "R": 1.008,
    "T": 0.997,
    "W": 0.904,
    "V": 0.931,
    "Y": 0.929,
}


# ============================================================
# pKa VALUES
# ============================================================

# Bjellqvist/ProtParam-style pKa values used for sequence-level
# charge and theoretical pI calculations.
POSITIVE_PKA = {
    # Bjellqvist pKa scale for theoretical pI / charge.
    "Nterm": 7.50,
    "K": 10.00,
    "R": 12.00,
    "H": 5.98,
}

NEGATIVE_PKA = {
    # Bjellqvist pKa scale for theoretical pI / charge.
    "Cterm": 3.55,
    "D": 4.05,
    "E": 4.45,
    "C": 9.00,
    "Y": 10.00,
}

# Terminal pKa values depend on the identity of the terminal residue.
# Keep the existing terminal-residue adjustment mechanism, but use
# values consistent with the same pKa scale.
NTERM_PKA = {
    # Bjellqvist terminal-residue adjustments.
    "A": 7.59,
    "E": 7.70,
    "M": 7.00,
    "P": 8.36,
    "S": 6.93,
    "T": 6.82,
    "V": 7.44,
}

CTERM_PKA = {
    # Bjellqvist terminal-residue adjustments.
    "D": 4.55,
    "E": 4.75,
}


# ============================================================
# INSTABILITY INDEX
# ============================================================

INSTABILITY_ORDER = (
    "ACEDGFIHKMLNQPSRTWVY"
)

INSTABILITY_ROWS = {

    "A": [
        1.0, 44.94, 1.0, -7.49, 1.0,
        1.0, 1.0, -7.49, 1.0, 1.0,
        1.0, 1.0, 1.0, 20.26, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0
    ],

    "C": [
        1.0, 1.0, 1.0, 20.26, 1.0,
        1.0, 1.0, 33.60, 1.0, 33.60,
        20.26, 1.0, -6.54, 20.26, 1.0,
        1.0, 33.60, 24.68, -6.54, 1.0
    ],

    "E": [
        1.0, 44.94, 33.60, 20.26, 1.0,
        1.0, 20.26, -6.54, 1.0, 1.0,
        1.0, 1.0, 20.26, 20.26, 20.26,
        1.0, 1.0, -14.03, 1.0, 1.0
    ],

    "D": [
        1.0, 1.0, 1.0, 1.0, 1.0,
        -6.54, 1.0, 1.0, -7.49, 1.0,
        1.0, 1.0, 1.0, 1.0, 20.26,
        -6.54, -14.03, 1.0, 1.0, 1.0
    ],

    "G": [
        -7.49, 1.0, -6.54, 1.0, 13.34,
        1.0, -7.49, 1.0, -7.49, 1.0,
        1.0, -7.49, 1.0, 1.0, 1.0,
        1.0, -7.49, 13.34, 1.0, -7.49
    ],

    "F": [
        1.0, 1.0, 1.0, 13.34, 1.0,
        1.0, 1.0, 1.0, -14.03, 1.0,
        1.0, 1.0, 1.0, 20.26, 1.0,
        1.0, 1.0, 1.0, 1.0, 33.60
    ],

    "I": [
        1.0, 1.0, 44.94, 1.0, 1.0,
        1.0, 1.0, 13.34, -7.49, 1.0,
        20.26, 1.0, 1.0, -1.88, 1.0,
        1.0, 1.0, 1.0, -7.49, 1.0
    ],

    "H": [
        1.0, 1.0, 1.0, 1.0, -9.37,
        -9.37, 44.94, 1.0, 24.68, 1.0,
        1.0, 24.68, 1.0, -1.88, 1.0,
        1.0, -6.54, -1.88, 1.0, 44.94
    ],

    "K": [
        1.0, 1.0, 1.0, 1.0, -7.49,
        1.0, -7.49, 1.0, 1.0, 33.60,
        -7.49, 1.0, 24.64, -6.54, 1.0,
        33.60, 1.0, 1.0, -7.49, 1.0
    ],

    "M": [
        13.34, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 58.28, 1.0, -1.88,
        1.0, 1.0, -6.54, 44.94, 44.94,
        -6.54, -1.88, 1.0, 1.0, 24.68
    ],

    "L": [
        1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, -7.49, 1.0,
        1.0, 1.0, 33.60, 20.26, 1.0,
        20.26, 1.0, 24.68, 1.0, 1.0
    ],

    "N": [
        1.0, -1.88, 1.0, 1.0, -14.03,
        -14.03, 44.94, 1.0, 24.68, 1.0,
        1.0, 1.0, -6.54, -1.88, 1.0,
        1.0, -7.49, -9.37, 1.0, 1.0
    ],

    "Q": [
        1.0, -6.54, 20.26, 20.26, 1.0,
        -6.54, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 20.26, 20.26, 44.94,
        1.0, 1.0, 1.0, -6.54, -6.54
    ],

    "P": [
        20.26, -6.54, 18.38, -6.54, 1.0,
        20.26, 1.0, 1.0, 1.0, -6.54,
        1.0, 1.0, 20.26, 20.26, 20.26,
        -6.54, 1.0, -1.88, 20.26, 1.0
    ],

    "S": [
        1.0, 33.60, 20.26, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 20.26, 44.94, 20.26,
        20.26, 1.0, 1.0, 1.0, 1.0
    ],

    "R": [
        1.0, 1.0, 1.0, 1.0, -7.49,
        1.0, 1.0, 20.26, 1.0, 1.0,
        1.0, 13.34, 20.26, 20.26, 44.94,
        58.28, 1.0, 58.28, 1.0, -6.54
    ],

    "T": [
        1.0, 1.0, 20.26, 1.0, -7.49,
        13.34, 1.0, 1.0, 1.0, 1.0,
        1.0, -14.03, -6.54, 1.0, 1.0,
        1.0, 1.0, -14.03, 1.0, 1.0
    ],

    "W": [
        -14.03, 1.0, 1.0, 1.0, -9.37,
        1.0, 1.0, 24.68, 1.0, 24.68,
        13.34, 13.34, 1.0, 1.0, 1.0,
        1.0, -14.03, 1.0, -7.49, 1.0
    ],

    "V": [
        1.0, 1.0, 1.0, -14.03, -7.49,
        1.0, 1.0, 1.0, -1.88, 1.0,
        1.0, 1.0, 1.0, 20.26, 1.0,
        1.0, -7.49, 1.0, 1.0, -6.54
    ],

    "Y": [
        24.68, 1.0, -6.54, 24.68, -7.49,
        1.0, 1.0, 13.34, 1.0, 44.94,
        1.0, 1.0, 1.0, 13.34, 1.0,
        -15.91, -7.49, -9.37, 1.0, 13.34
    ],
}


def build_instability_weights():

    weights = {}

    for first, values in INSTABILITY_ROWS.items():

        weights[first] = {
            amino_acid: value
            for amino_acid, value
            in zip(
                INSTABILITY_ORDER,
                values
            )
        }

    return weights


INSTABILITY_WEIGHTS = (
    build_instability_weights()
)
