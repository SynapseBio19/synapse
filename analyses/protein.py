import math
from itertools import product

from data.genetic_code import CODON_TABLE


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


# ============================================================
# BASIC VALIDATION
# ============================================================


def _strip_terminal_stop(sequence):
    """Drop one terminal "*" (a stop, not a residue)."""
    return sequence[:-1] if sequence.endswith(STOP_SYMBOL) else sequence


def validate_protein_sequence(protein):
    """True for IUPAC protein letters, "?" and at most a terminal "*"."""

    if not isinstance(protein, str):
        return False

    sequence = _strip_terminal_stop(
        "".join(protein.upper().split())
    )

    if not sequence:
        return False

    return all(
        residue in STANDARD_AMINO_ACIDS
        or residue in NONSTANDARD_RESIDUES
        or residue in UNRESOLVED_RESIDUES
        for residue in sequence
    )


def get_resolved_protein(protein):

    if not isinstance(protein, str):
        return ""

    return "".join(
        residue
        for residue in protein.upper()
        if residue in STANDARD_AMINO_ACIDS
    )


def get_positional_protein(protein):
    """Keep every position; non-standard residues become "?".

    Position-dependent analyses (sliding windows, dipeptide scales, signal
    peptide) must not delete residues: doing so shifts every later position.
    Standard residues are kept, every other symbol (B, Z, J, X, U, O, an
    internal "*") becomes "?", and the windows or dipeptides containing a
    "?" are skipped (audit P2-5, D14-F6). Whitespace and one terminal "*"
    (a stop) are removed; they are not residues.
    """
    if not isinstance(protein, str):
        return ""
    sequence = _strip_terminal_stop("".join(protein.upper().split()))
    return "".join(
        residue if residue in STANDARD_AMINO_ACIDS else "?"
        for residue in sequence
    )


def analyze_protein_quality(protein):

    if not isinstance(protein, str):
        protein = ""

    raw = "".join(
        protein.upper().split()
    )

    terminal_stop = raw.endswith(STOP_SYMBOL)
    sequence = _strip_terminal_stop(raw)

    standard = 0
    nonstandard = 0
    ambiguous = 0
    invalid = 0

    invalid_positions = []

    for position, residue in enumerate(
        sequence
    ):

        if residue in STANDARD_AMINO_ACIDS:
            standard += 1

        elif residue in NONSTANDARD_RESIDUES:
            nonstandard += 1

        elif residue in UNRESOLVED_RESIDUES:
            ambiguous += 1

        else:
            invalid += 1

            invalid_positions.append({
                "position": position,
                "position_1_based": position + 1,
                "residue": residue,
            })

    total = len(sequence)

    return {
        "length":
            total,

        "standard_residues":
            standard,

        "nonstandard_residues":
            nonstandard,

        "ambiguous_residues":
            ambiguous,

        "terminal_stop":
            terminal_stop,

        "invalid_residues":
            invalid,

        "standard_percentage":
            standard / total * 100
            if total else 0.0,

        "ambiguous_percentage":
            ambiguous / total * 100
            if total else 0.0,

        "invalid_percentage":
            invalid / total * 100
            if total else 0.0,

        "invalid_positions":
            invalid_positions,

        "valid":
            invalid == 0,

        "fully_resolved":
            (
                ambiguous == 0
                and invalid == 0
            ),
    }


def expand_iupac_codon(codon):

    codon = (
        codon
        .upper()
        .replace("T", "U")
    )

    if len(codon) != 3:
        return []

    possibilities = []

    for base in codon:

        if base not in RNA_IUPAC:
            return []

        possibilities.append(
            RNA_IUPAC[base]
        )

    return [
        "".join(combination)
        for combination in product(
            *possibilities
        )
    ]


def rna_codon_table(genetic_code=1):
    """Codon table (RNA alphabet, "STOP" for stops) for an NCBI genetic code."""
    if genetic_code in (None, 1):
        return CODON_TABLE
    from core.seqtools import codon_dict
    return codon_dict(int(genetic_code), rna=True, stop_symbol="STOP")


def resolve_codon(codon, codon_table=None):
    codon_table = CODON_TABLE if codon_table is None else codon_table

    codon = (
        codon
        .upper()
        .replace("T", "U")
    )

    expanded = expand_iupac_codon(
        codon
    )

    if not expanded:

        return {
            "codon": codon,
            "status": "invalid",
            "amino_acid": None,
            "possibilities": [],
            "expanded_codons": [],
            "contains_stop": False,
        }

    translations = set()

    for concrete_codon in expanded:

        amino_acid = codon_table.get(
            concrete_codon
        )

        if amino_acid is not None:

            translations.add(
                amino_acid
            )

    if not translations:

        return {
            "codon": codon,
            "status": "invalid",
            "amino_acid": None,
            "possibilities": [],
            "expanded_codons": expanded,
            "contains_stop": False,
        }

    contains_stop = (
        "STOP" in translations
    )

    amino_acids = sorted(
        amino_acid
        for amino_acid in translations
        if amino_acid != "STOP"
    )

    if len(translations) == 1:

        result = next(
            iter(translations)
        )

        if result == "STOP":

            return {
                "codon": codon,
                "status": "stop",
                "amino_acid": None,
                "possibilities": ["STOP"],
                "expanded_codons": expanded,
                "contains_stop": True,
            }

        return {
            "codon": codon,
            "status": "resolved",
            "amino_acid": result,
            "possibilities": [result],
            "expanded_codons": expanded,
            "contains_stop": False,
        }

    if not contains_stop:

        return {
            "codon": codon,
            "status":
                "ambiguous_amino_acid",
            "amino_acid": None,
            "possibilities":
                amino_acids,
            "expanded_codons": expanded,
            "contains_stop": False,
        }

    return {
        "codon": codon,
        "status": "ambiguous_stop",
        "amino_acid": None,
        "possibilities":
            amino_acids + ["STOP"],
        "expanded_codons": expanded,
        "contains_stop": True,
    }


# ============================================================
# TRANSLATION
# ============================================================


def translate_protein(rna):

    if not isinstance(rna, str):
        return ""

    rna = (
        "".join(rna.upper().split())
        .replace("T", "U")
    )

    protein = ""

    for i in range(
        0,
        len(rna) - 2,
        3
    ):

        codon = rna[
            i:i + 3
        ]

        result = resolve_codon(
            codon
        )

        if result["status"] == "stop":
            break

        if result["status"] == "resolved":

            protein += (
                result["amino_acid"]
            )

        else:

            protein += "?"

    return protein


def analyze_translation(rna, initiator_start=False, genetic_code=1):
    """Translate RNA in frame +1 with the standard genetic code.

    initiator_start=True states that the first codon is the selected start
    codon of an ORF. If it is an initiation codon of the genetic code
    (analyses.orf.is_initiation_codon; table 1: AUG, CUG, GUG, UUG) it is
    decoded as methionine, because initiator tRNA inserts (formyl-)Met at
    the start codon (audits P0-8, D1); any other first codon is read with
    the elongation code. For arbitrary sequence windows leave it False.
    genetic_code selects the NCBI translation table (default 1, standard).
    """
    from analyses.orf import is_initiation_codon

    if not isinstance(rna, str):
        return {
            "rna": "",
            "reading_frame": "+1",
            "protein": "",
            "protein_length": 0,
            "codon_count_available": 0,
            "codons_processed": 0,
            "codon_records": [],
            "start_codon": None,
            "starts_with_aug": False,
            "stop_found": False,
            "stop_codon": None,
            "stop_position": None,
            "trailing_nucleotides": 0,
            "ambiguous_codon_count": 0,
            "deterministic_ambiguous_count": 0,
            "uncertain_stop_count": 0,
            "invalid_codon_count": 0,
        }

    rna = (
        "".join(rna.upper().split())
        .replace("T", "U")
    )

    available_codons = (
        len(rna) // 3
    )

    trailing_nucleotides = (
        len(rna) % 3
    )

    codon_table = rna_codon_table(genetic_code)
    records = []

    protein = ""

    stop_found = False
    stop_codon = None
    stop_position = None

    ambiguous_codons = 0
    deterministic_ambiguous = 0
    uncertain_stop_codons = 0
    invalid_codons = 0

    for codon_index in range(
        available_codons
    ):

        start = codon_index * 3
        end = start + 3

        codon = rna[
            start:end
        ]

        result = resolve_codon(codon, codon_table)

        initiator = False
        if (
            initiator_start
            and codon_index == 0
            and is_initiation_codon(codon, genetic_code)
        ):
            result = dict(result)
            result["status"] = "resolved"
            result["amino_acid"] = "M"
            result["possibilities"] = ["M"]
            initiator = True

        iupac_ambiguous = any(
            base not in "AUCG"
            for base in codon
        )

        if iupac_ambiguous:
            ambiguous_codons += 1

        record = {
            "number":
                codon_index + 1,
            "codon":
                codon,
            "start":
                start + 1,
            "end":
                end,
            "status":
                result["status"],
            "amino_acid":
                result["amino_acid"],
            "possibilities":
                result["possibilities"],
            "contains_stop":
                result["contains_stop"],
            "iupac_ambiguous":
                iupac_ambiguous,
            "initiator":
                initiator,
        }

        records.append(record)

        if result["status"] == "stop":

            stop_found = True
            stop_codon = codon
            stop_position = start + 1
            break

        if result["status"] == "resolved":

            protein += (
                result["amino_acid"]
            )

            if iupac_ambiguous:
                deterministic_ambiguous += 1

        elif result["status"] == "ambiguous_amino_acid":

            protein += "?"

        elif result["status"] == "ambiguous_stop":

            protein += "?"
            uncertain_stop_codons += 1

        else:

            protein += "?"
            invalid_codons += 1

    starts_with_aug = (
        len(rna) >= 3
        and rna[:3] == "AUG"
    )

    return {
        "rna":
            rna,

        "reading_frame":
            "+1",

        "genetic_code":

            int(genetic_code or 1),

        "protein":
            protein,

        "protein_length":
            len(protein),

        "codon_count_available":
            available_codons,

        "codons_processed":
            len(records),

        "codon_records":
            records,

        "start_codon":
            rna[:3]
            if len(rna) >= 3
            else None,

        "starts_with_aug":
            starts_with_aug,

        "stop_found":
            stop_found,

        "stop_codon":
            stop_codon,

        "stop_position":
            stop_position,

        "trailing_nucleotides":
            trailing_nucleotides,

        "ambiguous_codon_count":
            ambiguous_codons,

        "deterministic_ambiguous_count":
            deterministic_ambiguous,

        "uncertain_stop_count":
            uncertain_stop_codons,

        "invalid_codon_count":
            invalid_codons,

        "initiator_start":
            bool(initiator_start),
    }


def analyze_protein(protein):

    if not isinstance(protein, str):
        protein = ""

    sequence = "".join(
        protein.upper().split()
    )

    counts = {
        amino_acid: 0
        for amino_acid in AA_ORDER
    }

    ambiguous = 0
    nonstandard_counts = {}

    for residue in sequence:

        if residue in STANDARD_AMINO_ACIDS:
            counts[residue] += 1

        elif residue in NONSTANDARD_RESIDUES:
            nonstandard_counts[residue] = nonstandard_counts.get(residue, 0) + 1

        elif residue in UNRESOLVED_RESIDUES:
            ambiguous += 1

    counts = {
        amino_acid: count
        for amino_acid, count
        in counts.items()
        if count > 0
    }

    return {
        "counts":
            counts,

        "nonstandard_counts":
            nonstandard_counts,

        "ambiguous_count":
            ambiguous,

        "resolved_length":
            sum(counts.values()),

        "total_length":
            len(sequence),
    }


def calculate_amino_acid_percentages(
    protein
):

    analysis = analyze_protein(
        protein
    )

    resolved_length = (
        analysis["resolved_length"]
    )

    if resolved_length == 0:
        return {}

    return {
        amino_acid:
            count
            / resolved_length
            * 100
        for amino_acid, count
        in analysis["counts"].items()
    }


# ============================================================
# MOLECULAR WEIGHT
# ============================================================


def calculate_protein_molecular_weight(
    protein,
    translation_analysis=None
):

    if not isinstance(protein, str):
        return {
            "min": None,
            "max": None,
            "ambiguous": False,
            "undetermined": True,
            "known_mass": 0.0,
        }

    sequence = protein.upper()

    if not sequence:
        return {
            "min": 0.0,
            "max": 0.0,
            "ambiguous": False,
            "undetermined": False,
            "known_mass": 0.0,
        }

    if translation_analysis is None:

        # B, Z and J give a mass range over their residues; X and "?" leave
        # the mass undetermined; a terminal "*" is a stop, not a residue.
        sequence = _strip_terminal_stop("".join(sequence.split()))
        total = WATER_MASS
        total_min = WATER_MASS
        total_max = WATER_MASS
        unknown = False
        ranged = False

        for residue in sequence:

            if residue in AMINO_ACID_RESIDUE_MASSES:
                mass = AMINO_ACID_RESIDUE_MASSES[residue]
                total += mass
                total_min += mass
                total_max += mass
            elif residue in AMBIGUOUS_RESIDUE_OPTIONS:
                masses = [AMINO_ACID_RESIDUE_MASSES[option]
                          for option in AMBIGUOUS_RESIDUE_OPTIONS[residue]]
                total_min += min(masses)
                total_max += max(masses)
                ranged = True
            else:
                unknown = True

        if unknown:
            return {
                "min": None,
                "max": None,
                "ambiguous": True,
                "undetermined": True,
                "known_mass": total,
            }

        return {
            "min": total_min,
            "max": total_max,
            "ambiguous": ranged,
            "undetermined": False,
            "known_mass": total,
        }

    # Translation-aware calculation.
    #
    # The protein represented by an ambiguous codon can have
    # multiple possible masses. An ambiguous-stop codon can also
    # terminate translation, so a simple per-residue min/max is
    # insufficient. We therefore keep separate mass ranges for
    # translation paths that are still active and paths that have
    # already terminated.
    active_min = 0.0
    active_max = 0.0

    completed_min = None
    completed_max = None

    ambiguous = False
    known_mass = WATER_MASS
    protein_index = 0

    for record in translation_analysis.get(
        "codon_records",
        []
    ):

        status = record["status"]

        if status == "stop":

            if active_min is not None:
                completed_min = (
                    active_min
                    if completed_min is None
                    else min(
                        completed_min,
                        active_min
                    )
                )

                completed_max = (
                    active_max
                    if completed_max is None
                    else max(
                        completed_max,
                        active_max
                    )
                )

            active_min = None
            active_max = None
            break

        residue = (
            sequence[protein_index]
            if protein_index < len(sequence)
            else "?"
        )

        if residue in AMINO_ACID_RESIDUE_MASSES:

            mass = AMINO_ACID_RESIDUE_MASSES[
                residue
            ]

            active_min += mass
            active_max += mass
            known_mass += mass

        elif status in {
            "ambiguous_amino_acid",
            "ambiguous_stop",
        }:

            ambiguous = True

            possible_amino_acids = [
                amino_acid
                for amino_acid in record.get(
                    "possibilities",
                    []
                )
                if amino_acid
                in AMINO_ACID_RESIDUE_MASSES
            ]

            masses = [
                AMINO_ACID_RESIDUE_MASSES[
                    amino_acid
                ]
                for amino_acid
                in possible_amino_acids
            ]

            if not masses:

                return {
                    "min": None,
                    "max": None,
                    "ambiguous": True,
                    "undetermined": True,
                    "known_mass": known_mass,
                }

            if status == "ambiguous_stop":

                # STOP is one possible outcome. It terminates
                # the current active path without adding a residue.
                if "STOP" in record.get(
                    "possibilities",
                    []
                ):

                    stop_min = active_min
                    stop_max = active_max

                    completed_min = (
                        stop_min
                        if completed_min is None
                        else min(
                            completed_min,
                            stop_min
                        )
                    )

                    completed_max = (
                        stop_max
                        if completed_max is None
                        else max(
                            completed_max,
                            stop_max
                        )
                    )

                # Non-stop amino-acid outcomes continue.
                active_min += min(masses)
                active_max += max(masses)

            else:

                active_min += min(masses)
                active_max += max(masses)

        else:

            return {
                "min": None,
                "max": None,
                "ambiguous": True,
                "undetermined": True,
                "known_mass": known_mass,
            }

        protein_index += 1

    # If translation never terminated, the active paths are valid
    # final products as well.
    final_min = completed_min
    final_max = completed_max

    # The two bounds become None together (undetermined residue).
    if active_min is not None and active_max is not None:

        final_min = (
            active_min
            if final_min is None
            else min(
                final_min,
                active_min
            )
        )

        final_max = (
            active_max
            if final_max is None
            else max(
                final_max,
                active_max
            )
        )

    if final_min is None or final_max is None:

        return {
            "min": None,
            "max": None,
            "ambiguous": ambiguous,
            "undetermined": True,
            "known_mass": known_mass,
        }

    final_min += WATER_MASS
    final_max += WATER_MASS

    return {
        "min": final_min,
        "max": final_max,
        "ambiguous": ambiguous,
        "undetermined": False,
        "known_mass": known_mass,
    }


def _get_pka_tables(protein):

    positive = POSITIVE_PKA.copy()
    negative = NEGATIVE_PKA.copy()

    if protein:

        first = protein[0]
        last = protein[-1]

        if first in NTERM_PKA:

            positive["Nterm"] = (
                NTERM_PKA[first]
            )

        if last in CTERM_PKA:

            negative["Cterm"] = (
                CTERM_PKA[last]
            )

    return positive, negative


def calculate_charge_at_ph(
    protein,
    pH
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return 0.0

    try:
        pH = float(pH)
    except (TypeError, ValueError):
        raise ValueError("pH must be a numeric value.")

    if not math.isfinite(pH):
        raise ValueError("pH must be finite.")

    if not 0.0 <= pH <= 14.0:
        raise ValueError("pH must be between 0 and 14.")

    positive_pka, negative_pka = (
        _get_pka_tables(
            protein
        )
    )

    positive_charge = 0.0
    negative_charge = 0.0

    # Henderson-Hasselbalch form for basic groups:
    # fraction protonated = 1 / (1 + 10^(pH-pKa))
    for group, pka in positive_pka.items():

        count = (
            1
            if group == "Nterm"
            else protein.count(group)
        )

        if count:
            positive_charge += (
                count
                / (
                    1.0
                    + 10.0 ** (pH - pka)
                )
            )

    # Henderson-Hasselbalch form for acidic groups:
    # fraction deprotonated = 1 / (1 + 10^(pKa-pH))
    for group, pka in negative_pka.items():

        count = (
            1
            if group == "Cterm"
            else protein.count(group)
        )

        if count:
            negative_charge += (
                count
                / (
                    1.0
                    + 10.0 ** (pka - pH)
                )
            )

    return (
        positive_charge
        - negative_charge
    )


def calculate_isoelectric_point(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return None

    # The Bjellqvist pKa model is monotonic in net charge over
    # the pH range 0-14 for the sequence-level Henderson-
    # Hasselbalch model used here. Bisection therefore gives
    # the theoretical pI where calculated net charge is zero.
    low = 0.0
    high = 14.0

    low_charge = calculate_charge_at_ph(
        protein,
        low
    )
    high_charge = calculate_charge_at_ph(
        protein,
        high
    )

    if low_charge < 0 or high_charge > 0:
        return None

    for _ in range(200):

        middle = (
            low + high
        ) / 2.0

        charge = calculate_charge_at_ph(
            protein,
            middle
        )

        if charge > 0:
            low = middle
        else:
            high = middle

    return (
        low + high
    ) / 2.0


def calculate_gravy(protein):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return None

    return sum(
        KYTE_DOOLITTLE[
            residue
        ]
        for residue in protein
    ) / len(protein)


def calculate_hydropathy_profile(
    protein,
    window=9
):

    protein = get_positional_protein(
        protein
    )

    try:
        window = int(window)
    except (TypeError, ValueError):
        raise ValueError("window must be an integer.")

    if window < 1:
        raise ValueError("window must be at least 1.")

    if len(protein) < window:
        return []

    profile = []

    for start in range(
        len(protein) - window + 1
    ):

        segment = protein[
            start:
            start + window
        ]

        if "?" in segment:

            # Unresolved position: no score instead of a shifted window.

            continue

        score = (
            sum(
                KYTE_DOOLITTLE[residue]
                for residue in segment
            )
            / window
        )

        profile.append({
            "start":
                start + 1,

            "end":
                start + window,

            "score":
                score,

            "hydropathy":
                score,

            "sequence":
                segment
        })

    return profile


def find_hydrophobic_regions(
    protein,
    window=19,
    threshold=1.6
):

    # ExPASy ProtScale notes that windows of 19 or 21 residues
    # make membrane-spanning hydrophobic segments stand out;
    # 1.6 is a commonly used Kyte-Doolittle threshold.
    profile = calculate_hydropathy_profile(
        protein,
        window
    )

    return [
        region
        for region in profile
        if region["score"] >= threshold
    ]


def find_transmembrane_domains(
    protein,
    window=19,
    hydropathy_threshold=1.6,
    max_charged=2
):

    protein = get_positional_protein(
        protein
    )

    try:
        window = int(window)
    except (TypeError, ValueError):
        raise ValueError("window must be an integer.")

    if window < 3:
        raise ValueError("window must be at least 3.")

    if window % 2 == 0:
        raise ValueError(
            "Transmembrane window must be odd so that it has a central residue."
        )

    if len(protein) < window:
        return []

    candidates = []

    for start in range(
        len(protein) - window + 1
    ):

        segment = protein[
            start:
            start + window
        ]

        if "?" in segment:

            # Unresolved position: no score instead of a shifted window.

            continue

        hydropathy = (
            sum(
                KYTE_DOOLITTLE[residue]
                for residue in segment
            )
            / window
        )

        charged = sum(
            residue in POSITIVE_RESIDUES
            or residue in NEGATIVE_RESIDUES
            for residue in segment
        )

        if (
            hydropathy >= hydropathy_threshold
            and charged <= max_charged
        ):

            candidates.append({
                "start": start + 1,
                "end": start + window,
            })

    if not candidates:
        return []

    # Consecutive sliding windows overlap heavily. Merge only
    # overlapping candidate windows, then recompute all reported
    # properties on the merged sequence instead of carrying the
    # first window's score into the merged region.
    merged_intervals = []

    for candidate in candidates:

        if not merged_intervals:
            merged_intervals.append(
                candidate.copy()
            )
            continue

        previous = merged_intervals[-1]

        if candidate["start"] <= previous["end"]:

            previous["end"] = max(
                previous["end"],
                candidate["end"]
            )

        else:

            merged_intervals.append(
                candidate.copy()
            )

    domains = []

    for interval in merged_intervals:

        start = interval["start"]
        end = interval["end"]

        sequence = protein[
            start - 1:end
        ]

        hydropathy = (
            sum(
                KYTE_DOOLITTLE[residue]
                for residue in sequence
            )
            / len(sequence)
        )

        charged = sum(
            residue in POSITIVE_RESIDUES
            or residue in NEGATIVE_RESIDUES
            for residue in sequence
        )

        domains.append({
            "start": start,
            "end": end,
            "length": len(sequence),
            "hydropathy": hydropathy,
            "charged_residues": charged,
            "sequence": sequence,
            "method": (
                "Kyte-Doolittle sliding-window "
                "sequence heuristic"
            ),
        })

    return domains


def calculate_aliphatic_index(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return None

    length = len(protein)

    a = protein.count("A") * 100 / length
    v = protein.count("V") * 100 / length
    i = protein.count("I") * 100 / length
    leu = protein.count("L") * 100 / length

    return (
        a
        + 2.9 * v
        + 3.9 * (i + leu)
    )


# ============================================================
# AROMATICITY
# ============================================================

def calculate_aromaticity(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return None

    return sum(
        protein.count(residue)
        for residue in AROMATIC_RESIDUES
    ) / len(protein)


# ============================================================
# EXTINCTION COEFFICIENT
# ============================================================


def calculate_extinction_coefficient(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    tryptophan = protein.count("W")
    tyrosine = protein.count("Y")
    cysteine = protein.count("C")

    base = (
        tryptophan * 5500
        + tyrosine * 1490
    )

    # Pace et al. values at 280 nm in water.
    # All cysteines reduced: no disulfide contribution.
    reduced = base

    # Maximum possible number of disulfide bonds from the
    # available cysteines. If the Cys count is odd, one Cys
    # necessarily remains unpaired.
    disulfide_pairs = cysteine // 2
    unpaired_cysteine = cysteine % 2

    oxidized = (
        base
        + disulfide_pairs * 125
    )

    return {
        "reduced_cysteine": reduced,
        "oxidized_cysteine": oxidized,
        "tryptophan": tryptophan,
        "tyrosine": tyrosine,
        "cysteine": cysteine,
        "disulfide_pairs": disulfide_pairs,
        "unpaired_cysteine": unpaired_cysteine,
        "wavelength_nm": 280,
        "units": "M^-1 cm^-1",
    }


def calculate_instability_index(
    protein
):

    protein = get_positional_protein(
        protein
    )

    if len(protein) < 2:
        return None

    score = 0.0
    scored = 0

    for i in range(
        len(protein) - 1
    ):

        first = protein[i]
        second = protein[i + 1]

        if first == "?" or second == "?":
            # Unresolved residue: this dipeptide has no weight. Skipping it
            # keeps every other dipeptide at its true position.
            continue

        scored += 1

        score += (
            INSTABILITY_WEIGHTS[
                first
            ][
                second
            ]
        )

    if scored == 0:
        return None

    # 10/L * sum over dipeptides (Guruprasad et al. 1990); L is the number of
    # scored residues, which equals len(protein) for a fully resolved sequence.
    return (
        10.0
        / (scored + 1)
        * score
    )


# ============================================================
# SECONDARY STRUCTURE TENDENCIES
# ============================================================


def calculate_secondary_structure_propensity(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:

        return {
            "helix": 0.0,
            "turn": 0.0,
            "sheet": 0.0,
            "method": "Chou-Fasman",
            "units": "propensity",
        }

    # Classic Chou-Fasman conformational parameters:
    # P(alpha), P(beta), P(turn).
    # Values above 1 indicate a relative preference for that
    # conformational class; an average propensity is NOT a
    # percentage of residues adopting that structure.
    chou_fasman = {
        "A": (1.42, 0.83, 0.66),
        "R": (0.98, 0.93, 0.95),
        "N": (0.67, 0.89, 1.56),
        "D": (1.01, 0.54, 1.46),
        "C": (0.70, 1.19, 1.19),
        "Q": (1.11, 1.10, 0.98),
        "E": (1.51, 0.37, 0.74),
        "G": (0.57, 0.75, 1.56),
        "H": (1.00, 0.87, 0.95),
        "I": (1.08, 1.60, 0.47),
        "L": (1.21, 1.30, 0.59),
        "K": (1.16, 0.74, 1.01),
        "M": (1.45, 1.05, 0.60),
        "F": (1.13, 1.38, 0.60),
        "P": (0.57, 0.55, 1.52),
        "S": (0.77, 0.75, 1.43),
        "T": (0.83, 1.19, 0.96),
        "W": (1.08, 1.37, 0.96),
        "Y": (0.69, 1.47, 1.14),
        "V": (1.06, 1.70, 0.50),
    }

    length = len(protein)

    return {
        "helix": (
            sum(
                chou_fasman[residue][0]
                for residue in protein
            )
            / length
        ),
        "turn": (
            sum(
                chou_fasman[residue][2]
                for residue in protein
            )
            / length
        ),
        "sheet": (
            sum(
                chou_fasman[residue][1]
                for residue in protein
            )
            / length
        ),
        "method": "Chou-Fasman",
        "units": "propensity",
    }


def calculate_secondary_structure_fraction(
    protein
):
    """
    Backward-compatible wrapper.

    The historical function name is retained so existing callers
    do not break, but the returned values are propensities, not
    structural fractions.
    """
    return calculate_secondary_structure_propensity(
        protein
    )


def calculate_flexibility_profile(
    protein,
    window=9
):

    protein = get_positional_protein(
        protein
    )

    try:
        window = int(window)
    except (TypeError, ValueError):
        raise ValueError("window must be an integer.")

    if window < 1:
        raise ValueError("window must be at least 1.")

    if len(protein) < window:
        return []

    profile = []

    for start in range(
        len(protein) - window + 1
    ):

        segment = protein[
            start:
            start + window
        ]

        if "?" in segment:

            # Unresolved position: no score instead of a shifted window.

            continue

        score = (
            sum(
                FLEXIBILITY_SCALE[residue]
                for residue in segment
            )
            / window
        )

        profile.append({
            "start": start + 1,
            "end": start + window,
            "score": score,
            "sequence": segment,
            "scale": "Vihinen normalized flexibility",
        })

    return profile


def predict_signal_peptide(
    protein
):

    protein = get_positional_protein(
        protein
    )

    if len(protein) < 15:

        return {
            "predicted": False,
            "evidence": "low",
            "reason": "Sequence too short for a meaningful signal-peptide screen.",
            "hydrophobic_core": None,
            "cleavage_site": None,
            "method": "N-terminal sequence heuristic",
        }

    # This is deliberately a screen, not a SignalP-equivalent
    # predictor. A classical secretory signal peptide usually has
    # an N-region, a hydrophobic h-region and a C-region containing
    # a plausible cleavage site.
    n_terminal = protein[
        :min(35, len(protein))
    ]

    hydrophobic_core = None

    for window in range(8, 16):

        if len(n_terminal) < window:
            continue

        for start in range(
            0,
            min(
                len(n_terminal) - window + 1,
                13
            )
        ):

            segment = n_terminal[
                start:
                start + window
            ]

            if "?" in segment:
                # Unresolved position: no score instead of a shifted window.
                continue

            score = (
                sum(
                    KYTE_DOOLITTLE[residue]
                    for residue in segment
                )
                / window
            )

            charged = sum(
                residue in POSITIVE_RESIDUES
                or residue in NEGATIVE_RESIDUES
                for residue in segment
            )

            if (
                score >= 1.4
                and charged <= 1
            ):

                candidate = {
                    "start": start + 1,
                    "end": start + window,
                    "score": score,
                    "sequence": segment,
                    "charged_residues": charged,
                }

                if (
                    hydrophobic_core is None
                    or candidate["score"]
                    > hydrophobic_core["score"]
                ):
                    hydrophobic_core = candidate

    if hydrophobic_core is None:

        return {
            "predicted": False,
            "evidence": "low",
            "reason": (
                "No convincing N-terminal hydrophobic core "
                "was detected."
            ),
            "hydrophobic_core": None,
            "cleavage_site": None,
            "method": "N-terminal sequence heuristic",
        }

    # Simplified -3/-1 cleavage-site screen. Small residues are
    # favored at these positions in classical signal peptidase
    # cleavage motifs. This is only a plausibility filter.
    cleavage_residues = set(
        "AVGSTC"
    )

    cleavage_candidates = []

    minimum_cleavage = max(
        hydrophobic_core["end"] + 2,
        15
    )

    maximum_cleavage = min(
        30,
        len(n_terminal) - 1
    )

    for cleavage in range(
        minimum_cleavage,
        maximum_cleavage + 1
    ):

        minus_1 = n_terminal[
            cleavage - 1
        ]

        minus_3 = n_terminal[
            cleavage - 3
        ]

        if (
            minus_1 in cleavage_residues
            and minus_3 in cleavage_residues
            and minus_1 != "P"
            and minus_3 != "P"
        ):

            cleavage_candidates.append(
                cleavage
            )

    if not cleavage_candidates:

        return {
            "predicted": False,
            "evidence": "low",
            "reason": (
                "N-terminal hydrophobic core detected, "
                "but no plausible cleavage-site pattern "
                "was found."
            ),
            "hydrophobic_core": hydrophobic_core,
            "cleavage_site": None,
            "method": "N-terminal sequence heuristic",
        }

    cleavage_site = min(
        cleavage_candidates,
        key=lambda position:
            abs(
                position
                - (
                    hydrophobic_core["end"]
                    + 4
                )
            )
    )

    evidence = (
        "moderate"
        if hydrophobic_core["score"] >= 1.8
        else "low"
    )

    return {
        "predicted": True,
        "evidence": evidence,
        "reason": (
            "N-terminal hydrophobic core and a plausible "
            "-3/-1 cleavage-site pattern were detected."
        ),
        "hydrophobic_core": hydrophobic_core,
        "cleavage_site": cleavage_site,
        "method": "N-terminal sequence heuristic",
    }


def analyze_protein_properties(
    protein
):

    resolved = get_resolved_protein(
        protein
    )

    length = len(resolved)

    if length == 0:

        return {
            "resolved_residues":
                0,

            "hydrophobic_percentage":
                0.0,

            "polar_percentage":
                0.0,

            "positive_percentage":
                0.0,

            "negative_percentage":
                0.0,

            "aromatic_percentage":
                0.0
        }

    return {
        "resolved_residues":
            length,

        "hydrophobic_percentage":
            sum(
                residue
                in HYDROPHOBIC_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100,

        "polar_percentage":
            sum(
                residue
                in POLAR_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100,

        "positive_percentage":
            sum(
                residue
                in POSITIVE_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100,

        "negative_percentage":
            sum(
                residue
                in NEGATIVE_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100,

        "aromatic_percentage":
            sum(
                residue
                in AROMATIC_RESIDUES
                for residue
                in resolved
            )
            / length
            * 100
    }


# ============================================================
# COMPLETE PROTEIN ANALYSIS
# ============================================================


def analyze_complete_protein(
    protein,
    rna=None,
    initiator_start=False,
    genetic_code=1,
):

    if not isinstance(protein, str):
        protein = ""

    protein = "".join(
        protein.upper().split()
    )

    translation = None

    if rna is not None:

        translation = analyze_translation(
            rna,
            initiator_start=initiator_start,
            genetic_code=genetic_code,
        )

        protein = translation[
            "protein"
        ]

    quality = (
        analyze_protein_quality(
            protein
        )
    )

    composition = (
        analyze_protein(
            protein
        )
    )

    percentages = (
        calculate_amino_acid_percentages(
            protein
        )
    )

    properties = (
        analyze_protein_properties(
            protein
        )
    )

    molecular_weight = (
        calculate_protein_molecular_weight(
            protein,
            translation
        )
    )

    pI = (
        calculate_isoelectric_point(
            protein
        )
    )

    charges = {}

    for ph in (
        5.0,
        7.0,
        7.4,
        9.0
    ):

        charges[str(ph)] = (
            calculate_charge_at_ph(
                protein,
                ph
            )
        )

    gravy = calculate_gravy(
        protein
    )

    hydropathy = (
        calculate_hydropathy_profile(
            protein,
            window=9
        )
    )

    hydrophobic_regions = (
        find_hydrophobic_regions(
            protein,
            window=19,
            threshold=1.6
        )
    )

    tm_domains = (
        find_transmembrane_domains(
            protein,
            window=19,
            hydropathy_threshold=1.6,
            max_charged=2
        )
    )

    instability = (
        calculate_instability_index(
            protein
        )
    )

    aliphatic = (
        calculate_aliphatic_index(
            protein
        )
    )

    aromaticity = (
        calculate_aromaticity(
            protein
        )
    )

    extinction = (
        calculate_extinction_coefficient(
            protein
        )
    )

    secondary_structure = (
        calculate_secondary_structure_propensity(
            protein
        )
    )

    flexibility = (
        calculate_flexibility_profile(
            protein,
            window=9
        )
    )

    signal_peptide = (
        predict_signal_peptide(
            protein
        )
    )

    return {
        "protein":
            protein,

        "length":
            len(protein),

        "quality":
            quality,

        "composition":
            composition,

        "percentages":
            percentages,

        "properties":
            properties,

        "molecular_weight":
            molecular_weight,

        "pI":
            pI,

        "charges":
            charges,

        "gravy":
            gravy,

        "hydropathy":
            hydropathy,

        "hydrophobic_regions":
            hydrophobic_regions,

        "transmembrane_domains":
            tm_domains,

        "instability_index":
            instability,

        "aliphatic_index":
            aliphatic,

        "aromaticity":
            aromaticity,

        "extinction_coefficient":
            extinction,

        "secondary_structure":
            secondary_structure,

        "flexibility":
            flexibility,

        "signal_peptide":
            signal_peptide,

        "translation":
            translation
    }


def create_protein_preview(
    protein,
    max_length=60
):

    if len(protein) <= max_length:
        return protein

    half = max_length // 2

    return (
        protein[:half]
        + "..."
        + protein[-half:]
    )


# ============================================================
# REPORT
# ============================================================

# The terminal dialogs of this module are in ui/protein_menu.py (P2-2).


# ============================================================
# PROTEIN ANALYSIS UI - SYNAPSE PROFESSIONAL V5
# ============================================================
#
# The analysis engine above remains unchanged.
#
# Menu:
#   [1] Protein Overview
#   [2] Translation & Sequence
#   [3] Amino Acid Composition
#   [4] Physicochemical & Structural Properties
#   [5] Complete Protein Report
#   [0] Back
#
# ============================================================


# ============================================================
# 2 - TRANSLATION & SEQUENCE
# ============================================================


# ============================================================
# 3 - AMINO ACID COMPOSITION
# ============================================================


# ============================================================
# 4 - PHYSICOCHEMICAL & STRUCTURAL PROPERTIES
# ============================================================


# ============================================================
# 5 - COMPLETE PROTEIN REPORT
# ============================================================


# ============================================================
# MENU
# ============================================================


# ============================================================
# MAIN PROTEIN UI
# ============================================================


# ============================================================
# PIPELINE ADAPTER
# ============================================================
