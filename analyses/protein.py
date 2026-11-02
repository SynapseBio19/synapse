import math
from itertools import product

from data.genetic_code import CODON_TABLE
from data.amino_acids import AMINO_ACIDS
from ui.style import major_rule, minor_rule
from ui import style


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

    if active_min is not None:

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
    l = protein.count("L") * 100 / length

    return (
        a
        + 2.9 * v
        + 3.9 * (i + l)
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

def print_protein_report(
    protein,
    amino_acid_counts=None,
    rna=None
):

    analysis = analyze_complete_protein(
        protein,
        rna
    )

    protein = analysis["protein"]

    style.header('PROTEIN ANALYSIS')

    # ========================================================
    # TRANSLATION
    # ========================================================

    translation = analysis[
        "translation"
    ]

    if translation is not None:

        style.section('Translation Information')

        style.field('Reading Frame', f"{translation['reading_frame']}")

        style.field('First Codon', f"{translation['start_codon'] or 'N/A'}")

        style.field('AUG at Start', f"{'Yes' if translation['starts_with_aug'] else 'No'}")

        style.field('Codons Processed', f"{translation['codons_processed']}")

        if translation["stop_found"]:

            style.field('Stop Codon', f"{translation['stop_codon']}")

            style.field('Stop Position', f"{translation['stop_position']}")

        else:

            style.field('Stop Codon', "Not detected")

        style.field('Trailing Bases', f"{translation['trailing_nucleotides']}")

    # ========================================================
    # BASIC INFORMATION
    # ========================================================

    style.section('Protein Sequence')

    style.field('Sequence', f"{create_protein_preview(protein)}")

    style.field('Protein Length', f"{analysis['length']} aa")

    style.field('Resolved Residues', f"{analysis['quality']['standard_residues']}")

    style.field('Ambiguous Residues', f"{analysis['quality']['ambiguous_residues']}")

    # ========================================================
    # QUALITY
    # ========================================================

    quality = analysis[
        "quality"
    ]

    style.section('Protein Quality')

    style.field('Standard Residues', f"{quality['standard_residues']}")

    if quality.get("nonstandard_residues"):
        print(
            f"Sec/Pyl (U, O):      "
            f"{quality['nonstandard_residues']}"
        )

    style.field('Ambiguous Residues', f"{quality['ambiguous_residues']}")

    style.field('Invalid Residues', f"{quality['invalid_residues']}")

    style.field('Ambiguous Content', f"{quality['ambiguous_percentage']:.2f}%")

    if not quality["valid"]:

        style.field('Status', "INVALID")

    elif not quality["fully_resolved"]:

        style.field('Status', "VALID - "
            "contains unresolved residues")

    else:

        style.field('Status', "VALID - "
            "fully resolved")

    # ========================================================
    # MOLECULAR WEIGHT
    # ========================================================

    mw = analysis[
        "molecular_weight"
    ]

    style.section('Molecular Weight')

    if mw.get("undetermined"):

        if mw.get("known_mass") is not None:
            style.field('Estimated MW', "Undetermined")

            style.field('Known-residue mass', f"{mw['known_mass']:,.0f} Da")
        else:
            style.field('Estimated MW', "Undetermined")

        style.field('Status', "Unresolved sequence prevents an exact mass")

    elif mw["ambiguous"]:

        style.field('Estimated MW', f"{mw['min']:,.0f} - "
            f"{mw['max']:,.0f} Da")

        style.field('Status', "Range due to ambiguous translation")

    else:

        style.field('Estimated MW', f"{mw['min']:,.0f} Da")

    # ========================================================
    # AMINO ACID COMPOSITION
    # ========================================================

    composition = analysis[
        "composition"
    ]

    percentages = analysis[
        "percentages"
    ]

    style.section('Amino Acid Composition')

    widths = (3, 4, 16, 11, 6)
    print(style.columns(("AA", "3L", "Name", "Type", "Count", "Percent"), widths))

    print(minor_rule())

    for code in AA_ORDER:

        count = composition[
            "counts"
        ].get(
            code,
            0
        )

        if count == 0:
            continue

        amino_acid = (
            AMINO_ACIDS[code]
        )

        print(style.columns((
            code, amino_acid["three_letter"], amino_acid["name"],
            amino_acid["type"], count, f"{percentages[code]:.2f}%",
        ), widths))

    if composition[
        "ambiguous_count"
    ] > 0:

        print()

        print(
            f"?  | Unresolved residues | "
            f"Count: "
            f"{composition['ambiguous_count']}"
        )

        print(
            "Percentages exclude unresolved residues."
        )

    # ========================================================
    # COMPOSITION PROPERTIES
    # ========================================================

    properties = analysis[
        "properties"
    ]

    style.section('Composition Properties')

    style.field('Hydrophobic', f"{properties['hydrophobic_percentage']:.2f}%")

    style.field('Polar', f"{properties['polar_percentage']:.2f}%")

    style.field('Positive', f"{properties['positive_percentage']:.2f}%")

    style.field('Negative', f"{properties['negative_percentage']:.2f}%")

    style.field('Aromatic', f"{properties['aromatic_percentage']:.2f}%")

    # ========================================================
    # CHARGE / pI
    # ========================================================

    style.section('Charge Properties')

    if analysis["pI"] is not None:

        style.field('Theoretical pI', f"{analysis['pI']:.2f}")

        for ph, charge in (
            analysis["charges"].items()
        ):

            print(
                f"Net Charge at pH "
                f"{ph:>3}:        "
                f"{charge:+.2f}"
            )

    else:

        print(
            "Charge analysis unavailable."
        )

    # ========================================================
    # HYDROPATHY
    # ========================================================

    style.section('Hydropathy')

    if analysis["gravy"] is not None:

        style.field('GRAVY', f"{analysis['gravy']:+.3f}")

        if analysis["gravy"] > 0.5:

            style.field('Interpretation', "Overall hydrophobic tendency")

        elif analysis["gravy"] < -0.5:

            style.field('Interpretation', "Overall hydrophilic tendency")

        else:

            style.field('Interpretation', "Mixed / intermediate")

    else:

        print(
            "GRAVY unavailable."
        )

    # ========================================================
    # MOLECULAR PROPERTIES
    # ========================================================

    style.section('Protein Molecular Properties')

    instability = (
        analysis[
            "instability_index"
        ]
    )

    if instability is not None:

        style.field('Instability Index', f"{instability:.2f}")

        if analysis["length"] >= 15:
            classification = (
                "Predicted stable"
                if instability < 40
                else "Potentially unstable"
            )
        else:
            classification = (
                "Not classified - sequence too short"
            )

        style.field('Stability Class', f"{classification}")

    else:

        style.field('Instability Index', "Not available")

    if analysis[
        "aliphatic_index"
    ] is not None:

        style.field('Aliphatic Index', f"{analysis['aliphatic_index']:.2f}")

    if analysis[
        "aromaticity"
    ] is not None:

        style.field('Aromaticity', f"{analysis['aromaticity']:.3f}")

    # ========================================================
    # EXTINCTION
    # ========================================================

    extinction = analysis[
        "extinction_coefficient"
    ]

    style.section('Extinction Coefficient')

    style.field('Wavelength', "280 nm")

    style.field('Tryptophan', f"{extinction['tryptophan']}")

    style.field('Tyrosine', f"{extinction['tyrosine']}")

    style.field('Cysteine', f"{extinction['cysteine']}")

    style.field('All Cys reduced', f"{extinction['reduced_cysteine']} "
        f"M^-1 cm^-1")

    style.field('Max. disulfide form', f"{extinction['oxidized_cysteine']} "
        f"M^-1 cm^-1")

    style.field('Disulfide pairs', f"{extinction['disulfide_pairs']}")

    if extinction["unpaired_cysteine"]:
        style.field('Unpaired cysteine', "1")

    # ========================================================
    # SECONDARY STRUCTURE
    # ========================================================

    secondary = analysis[
        "secondary_structure"
    ]

    style.section('Secondary Structure Tendencies')

    style.field('Helix Propensity', f"{secondary['helix']:.3f}")

    style.field('Turn Propensity', f"{secondary['turn']:.3f}")

    style.field('Sheet Propensity', f"{secondary['sheet']:.3f}")

    style.field('Method', "Chou-Fasman residue propensities")

    print(
        "Note: Propensity values are not "
        "percentages and do not constitute "
        "3D structure prediction."
    )

    # ========================================================
    # HYDROPHOBIC REGIONS
    # ========================================================

    regions = analysis[
        "hydrophobic_regions"
    ]

    style.section('Hydrophobic Regions')

    if regions:

        for region in regions[:10]:

            print(
                f"{region['start']}-"
                f"{region['end']} | "
                f"{region['hydropathy']:+.2f} | "
                f"{region['sequence']}"
            )

        if len(regions) > 10:

            print(
                f"... {len(regions) - 10} "
                f"additional windows"
            )

    else:

        print(
            "No strongly hydrophobic windows detected."
        )

    # ========================================================
    # TRANSMEMBRANE
    # ========================================================

    tm_domains = analysis[
        "transmembrane_domains"
    ]

    style.section('Transmembrane Domain Heuristic')

    if tm_domains:

        for index, domain in enumerate(
            tm_domains,
            start=1
        ):

            print(
                f"TM {index}:             "
                f"{domain['start']}-"
                f"{domain['end']} | "
                f"{domain['hydropathy']:+.2f}"
            )

    else:

        print(
            "No strong TM-like region detected."
        )

    print(
        "Note: Sequence-only heuristic."
    )

    # ========================================================
    # SIGNAL PEPTIDE
    # ========================================================

    signal = analysis[
        "signal_peptide"
    ]

    style.section('Signal Peptide Heuristic')

    if signal["predicted"]:

        style.field('Prediction', "Possible N-terminal signal peptide")

        style.field('Heuristic evidence', f"{signal['evidence']}")

        core = signal.get(
            "hydrophobic_core"
        )

        if core:

            style.field('Hydrophobic Core', f"{core['start']}-"
                f"{core['end']}")

            style.field('Core Hydropathy', f"{core['score']:+.2f}")

        if signal.get("cleavage_site") is not None:

            print(
                f"Cleavage-site screen:"
                f" {signal['cleavage_site']}"
            )

    else:

        style.field('Prediction', "No convincing pattern")

    print(
        "Note: Sequence-only heuristic."
    )

    # ========================================================
    # FLEXIBILITY
    # ========================================================

    flexibility = analysis[
        "flexibility"
    ]

    style.section('Flexibility Analysis')

    if flexibility:

        scores = [
            entry["score"]
            for entry in flexibility
        ]

        style.field('Average', f"{sum(scores) / len(scores):.3f}")

        style.field('Minimum', f"{min(scores):.3f}")

        style.field('Maximum', f"{max(scores):.3f}")

    else:

        print(
            "Sequence too short for "
            "window-based analysis."
        )

    # ========================================================
    # AMBIGUOUS CODONS
    # ========================================================

    if translation is not None:

        ambiguous_records = [
            record
            for record
            in translation[
                "codon_records"
            ]
            if (
                record["iupac_ambiguous"]
                or record["status"]
                in {
                    "ambiguous_amino_acid",
                    "ambiguous_stop",
                    "invalid"
                }
            )
        ]

        if ambiguous_records:

            style.section('Ambiguous Codon Analysis')

            for record in (
                ambiguous_records
            ):

                possibilities = (
                    record["possibilities"]
                )

                result = (
                    "/".join(
                        possibilities
                    )
                    if possibilities
                    else
                    "Unresolved"
                )

                print(
                    f"Codon "
                    f"{record['number']:>4}: "
                    f"{record['codon']} | "
                    f"{result}"
                )

    # ========================================================
    # SCIENTIFIC INTERPRETATION
    # ========================================================

    style.section('Scientific Interpretation')

    print(
        "Protein properties are calculated from the "
        "resolved amino-acid sequence."
    )

    print(
        "Molecular weight is an estimated sequence-based "
        "value and does not account for post-translational "
        "modifications."
    )

    print(
        "Theoretical pI and charge are estimated from "
        "ionizable groups and sequence-derived pKa values."
    )

    print(
        "GRAVY describes the overall hydropathy of the "
        "protein sequence."
    )

    print(
        "Instability Index is a statistical sequence-based "
        "indicator and is not an experimental stability test."
    )

    print(
        "Hydrophobic regions, transmembrane domains and "
        "signal peptides are heuristic sequence predictions."
    )

    print(
        "Secondary-structure values are Chou-Fasman "
        "propensities and are not a 3D structure prediction."
    )

    print(
        "Ambiguous IUPAC codons are explicitly tracked "
        "instead of silently producing potentially incorrect "
        "amino-acid assignments."
    )

    print(major_rule())


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


def _protein_ui_pause():
    style.pause()


def _frame_display(frame):
    if frame is None:
        return "N/A"

    value = str(frame)

    if not value.startswith(("+", "-")):
        value = "+" + value

    return value


def _sequence_applicability(length):
    """
    Qualitative applicability label only.
    It is NOT a statistical confidence probability.
    """
    if length < 15:
        return (
            "LOW",
            "Very short sequence; individual residues can "
            "strongly affect sequence-derived descriptors."
        )

    if length < 30:
        return (
            "MODERATE",
            "Short sequence; sequence-derived descriptors "
            "should be interpreted cautiously."
        )

    return (
        "STANDARD",
        "Sequence length is more suitable for the "
        "displayed sequence-level descriptors."
    )


def _protein_limitations_block(analysis):
    level, reason = _sequence_applicability(
        analysis["length"]
    )

    style.section('ANALYSIS APPLICABILITY')

    style.field('Level', f"{level}")

    style.field('Reason', f"{reason}")

    print()
    print(
        "Applicability level is qualitative; "
        "it is not a statistical confidence probability."
    )

    print(
        "All predictions are sequence-level descriptors "
        "or heuristics."
    )

    print(
        "They do not establish biological function, "
        "experimental stability, 3D structure, "
        "cellular localization or biological activity."
    )


def _protein_input_context(
    analysis,
    source="Selected ORF",
    orf_id=None,
    orf_length=None,
    strand=None,
    reading_frame=None,
    orf_status=None,
):
    style.section('INPUT')

    style.field('Source', f"{source}")

    if orf_id is not None:
        style.field('ORF', f"{orf_id}")

    if orf_length is not None:
        style.field('ORF length', f"{orf_length} bp")

    if strand is not None:
        style.field('Strand', f"{strand}")

    if reading_frame is not None:
        style.field('Reading frame', f"{_frame_display(reading_frame)}")

    if orf_status is not None:
        style.field('ORF status', f"{orf_status}")

    style.field('Protein length', f"{analysis['length']} aa")


def _protein_status_line(analysis):
    quality = analysis["quality"]

    if not quality["valid"]:
        return "INVALID"

    if not quality["fully_resolved"]:
        return "VALID - unresolved residues present"

    return "VALID - fully resolved"


def _protein_quality_block(analysis):
    quality = analysis["quality"]

    style.section('SEQUENCE QUALITY')

    style.field('Standard residues', f"{quality['standard_residues']}")

    if quality.get("nonstandard_residues"):
        print(
            f"Sec/Pyl (U, O):      "
            f"{quality['nonstandard_residues']}"
        )

    style.field('Ambiguous residues', f"{quality['ambiguous_residues']}")

    style.field('Invalid residues', f"{quality['invalid_residues']}")

    style.field('Status', f"{_protein_status_line(analysis)}")


def _protein_composition_block(analysis):
    properties = analysis["properties"]

    style.section('COMPOSITION')

    style.field('Hydrophobic', f"{properties['hydrophobic_percentage']:.2f}%")

    style.field('Polar', f"{properties['polar_percentage']:.2f}%")

    style.field('Positive', f"{properties['positive_percentage']:.2f}%")

    style.field('Negative', f"{properties['negative_percentage']:.2f}%")

    style.field('Aromatic', f"{properties['aromatic_percentage']:.2f}%")


def _protein_key_properties(analysis):
    style.section('KEY PROPERTIES')

    mw = analysis["molecular_weight"]

    if mw["ambiguous"]:
        style.field('Molecular weight', f"{mw['min']:,.0f} - "
            f"{mw['max']:,.0f} Da")
    else:
        style.field('Molecular weight', f"{mw['min']:,.0f} Da")

    if analysis["pI"] is not None:
        style.field('Theoretical pI', f"{analysis['pI']:.2f}")

    if analysis["gravy"] is not None:
        style.field('GRAVY', f"{analysis['gravy']:+.3f}")

    if analysis["instability_index"] is not None:
        style.field('Instability index', f"{analysis['instability_index']:.2f}")

    style.field('Aliphatic index', f"{analysis['aliphatic_index']:.2f}")

    style.field('Aromaticity', f"{analysis['aromaticity']:.3f}")


def _translation_quality_block(translation):
    style.section('TRANSLATION QUALITY')

    style.field('Start codon', f"{translation['start_codon'] or 'Not detected'}")

    style.field('Start detected', f"{'Yes' if translation['start_codon'] else 'No'}")

    style.field('Stop codon', f"{translation['stop_codon'] if translation['stop_found'] else 'Not detected'}")

    style.field('Complete endpoint', f"{'Yes' if translation['stop_found'] else 'No'}")

    print(
        f"Trailing nucleotides:"
        f" {translation['trailing_nucleotides']}"
    )

    style.field('Translation status', f"{'Complete' if translation['stop_found'] else 'Partial'}")


def _safe_extinction_values(analysis):
    extinction = analysis["extinction_coefficient"]

    reduced = extinction.get(
        "reduced",
        extinction.get(
            "reduced_coefficient",
            extinction.get(
                "reduced_cysteine",
                0
            )
        )
    )

    oxidized = extinction.get(
        "oxidized",
        extinction.get(
            "oxidized_coefficient",
            extinction.get(
                "oxidized_cysteine",
                0
            )
        )
    )

    return (
        extinction,
        reduced,
        oxidized
    )


def _hydropathy_interpretation(gravy):
    if gravy is None:
        return "Unavailable"

    if gravy >= 1.0:
        return "Strong hydrophobic tendency"

    if gravy >= 0.3:
        return "Moderate hydrophobic tendency"

    if gravy <= -1.0:
        return "Strong hydrophilic tendency"

    if gravy < -0.3:
        return "Hydrophilic tendency"

    return "Mixed hydrophobic / hydrophilic character"


def _charge_interpretation(analysis):
    charges = analysis.get(
        "charges",
        {}
    )

    charge = charges.get(
        "7.4"
    )

    if charge is None:
        charge = charges.get(
            "7.0"
        )

    if charge is None:
        return "Unavailable"

    if charge > 1:
        return "Net positive tendency near neutral pH"

    if charge < -1:
        return "Net negative tendency near neutral pH"

    return "Approximately neutral near neutral pH"



def _stability_interpretation(analysis):
    index = analysis.get(
        "instability_index"
    )

    if index is None:
        return "Unavailable - fewer than 2 resolved residues"

    if index < 40:
        return (
            "Sequence-based stable tendency "
            "(II < 40)"
        )

    return (
        "Sequence-based unstable tendency "
        "(II >= 40)"
    )


def _print_protein_overview(
    analysis,
    source="Selected ORF",
    orf_id=None,
    orf_length=None,
    strand=None,
    reading_frame=None,
    orf_status=None,
):
    style.header(
        "PROTEIN OVERVIEW"
    )

    _protein_input_context(
        analysis,
        source=source,
        orf_id=orf_id,
        orf_length=orf_length,
        strand=strand,
        reading_frame=reading_frame,
        orf_status=orf_status,
    )

    style.section('PROTEIN')

    style.field('Length', f"{analysis['length']} aa")

    style.field('Sequence', f"{analysis['protein']}")

    translation = analysis["translation"]

    if translation is not None:
        style.field('Translation', f"{'Complete' if translation['stop_found'] else 'Partial'}")

    print()
    _protein_quality_block(
        analysis
    )

    print()
    _protein_key_properties(
        analysis
    )

    print()
    _protein_composition_block(
        analysis
    )

    style.section('ANALYSIS QUALITY')

    level, reason = _sequence_applicability(
        analysis["length"]
    )

    style.field('Applicability level', f"{level}")

    style.field('Reason', f"{reason}")

    if analysis["length"] < 30:
        print()
        print(
            "This is a qualitative applicability label, "
            "not a statistical confidence probability."
        )

    print()
    print(
        "This overview contains sequence-level "
        "descriptors only."
    )


# ============================================================
# 2 - TRANSLATION & SEQUENCE
# ============================================================


def _print_translation_and_sequence(analysis):
    style.header(
        "TRANSLATION & SEQUENCE"
    )

    translation = analysis["translation"]

    if translation is None:
        print(
            "Translation context: Not supplied."
        )

        print(
            "Protein sequence was supplied directly."
        )

    else:
        style.section('TRANSLATION')

        style.field('Reading frame', f"{_frame_display(translation['reading_frame'])}")

        style.field('First codon', f"{translation['start_codon'] or 'N/A'}")

        style.field('AUG at start', f"{'Yes' if translation['starts_with_aug'] else 'No'}")

        style.field('Codons processed', f"{translation['codons_processed']}")

        _translation_quality_block(
            translation
        )

        style.section('CODON TRANSLATION')

        for record in translation[
            "codon_records"
        ]:
            aa = record[
                "amino_acid"
            ]

            if record[
                "status"
            ] == "stop":
                result = "STOP"

            elif aa is not None:
                result = aa

            elif record[
                "possibilities"
            ]:
                result = "/".join(
                    record[
                        "possibilities"
                    ]
                )

            else:
                result = "?"

            print(
                f"{record['number']:>3}: "
                f"{record['codon']} -> "
                f"{result}"
            )

    style.section('PROTEIN SEQUENCE')

    sequence = analysis[
        "protein"
    ]

    width = 60

    for start in range(
        0,
        len(sequence),
        width
    ):
        end = min(
            start + width,
            len(sequence)
        )

        print(
            f"{start + 1:>5}-{end:<5} "
            f"{sequence[start:end]}"
        )

    print()
    style.field('Length', f"{len(sequence)} aa")

    print()
    _protein_limitations_block(
        analysis
    )


# ============================================================
# 3 - AMINO ACID COMPOSITION
# ============================================================


def _print_amino_acid_composition(analysis):
    style.header(
        "AMINO ACID COMPOSITION"
    )

    composition = analysis[
        "composition"
    ]

    percentages = analysis[
        "percentages"
    ]

    counts = composition[
        "counts"
    ]

    widths = (3, 4, 16, 11, 6)
    print(style.columns(("AA", "3L", "Name", "Type", "Count", "Percent"), widths))

    print(minor_rule())

    ordered = sorted(
        counts.items(),
        key=lambda item: (
            -item[1],
            item[0]
        )
    )

    for code, count in ordered:
        amino_acid = AMINO_ACIDS[
            code
        ]

        print(style.columns((
            code, amino_acid["three_letter"], amino_acid["name"],
            amino_acid["type"], count, f"{percentages[code]:.2f}%",
        ), widths)
        )

    print(minor_rule())

    style.field('Resolved residues', f"{composition['resolved_length']}")

    style.field('Unresolved residues', f"{composition['ambiguous_count']}")

    style.section('COMPOSITION GROUPS')

    _protein_composition_block(
        analysis
    )

    if counts:
        residue, count = max(
            counts.items(),
            key=lambda item: item[1]
        )

        name = AMINO_ACIDS[
            residue
        ][
            "name"
        ]

        percent = percentages[
            residue
        ]

        print()
        style.field('Most abundant residue', f"{residue} ({name})")

        style.field('Occurrences', f"{count} / "
            f"{composition['resolved_length']}")

        style.field('Relative abundance', f"{percent:.2f}%")

    print()
    _protein_limitations_block(
        analysis
    )


# ============================================================
# 4 - PHYSICOCHEMICAL & STRUCTURAL PROPERTIES
# ============================================================


def _print_physicochemical_structural(analysis):
    style.header(
        "PHYSICOCHEMICAL & STRUCTURAL PROPERTIES"
    )

    print()
    _protein_limitations_block(
        analysis
    )

    style.section('MOLECULAR PROPERTIES')

    _protein_key_properties(
        analysis
    )

    print()
    style.field('Hydropathy', f"{_hydropathy_interpretation(analysis['gravy'])}")

    style.field('Charge', f"{_charge_interpretation(analysis)}")

    style.field('Stability', f"{_stability_interpretation(analysis)}")

    style.section('CHARGE PROFILE')

    for ph, charge in analysis[
        "charges"
    ].items():
        print(
            f"Net charge pH {ph:>4}: "
            f"{charge:+.2f}"
        )

    style.section('EXTINCTION COEFFICIENT')

    extinction, reduced, oxidized = (
        _safe_extinction_values(
            analysis
        )
    )

    style.field('Wavelength', "280 nm")

    style.field('Tryptophan', f"{extinction['tryptophan']}")

    style.field('Tyrosine', f"{extinction['tyrosine']}")

    style.field('Cysteine', f"{extinction['cysteine']}")

    style.field('All Cys reduced', f"{reduced} M^-1 cm^-1")

    style.field('Max. disulfide form', f"{oxidized} M^-1 cm^-1")

    style.field('Disulfide pairs', f"{extinction['disulfide_pairs']}")

    if extinction["unpaired_cysteine"]:
        style.field('Unpaired cysteine', "1")

    style.section('SECONDARY STRUCTURE SCREENING')

    secondary = analysis[
        "secondary_structure"
    ]

    style.field('Helix propensity', f"{secondary['helix']:.3f}")

    style.field('Turn propensity', f"{secondary['turn']:.3f}")

    style.field('Sheet propensity', f"{secondary['sheet']:.3f}")

    style.field('Interpretation', "Chou-Fasman sequence propensity")

    style.field('Limitation', "Not a 3D structure prediction")

    style.section('HYDROPHOBIC REGION SCREENING')

    regions = analysis[
        "hydrophobic_regions"
    ]

    if regions:
        for index, region in enumerate(
            regions,
            start=1
        ):
            print(
                f"Region {index}:        "
                f"{region['start']}-"
                f"{region['end']} | "
                f"Hydropathy "
                f"{region['hydropathy']:+.2f}"
            )
    else:
        print(
            "No strongly hydrophobic windows detected."
        )

    style.section('TRANSMEMBRANE SCREENING')

    tm_domains = analysis[
        "transmembrane_domains"
    ]

    if tm_domains:
        for index, domain in enumerate(
            tm_domains,
            start=1
        ):
            print(
                f"TM {index}:             "
                f"{domain['start']}-"
                f"{domain['end']} | "
                f"Hydropathy "
                f"{domain['hydropathy']:+.2f}"
            )
    else:
        style.field('Prediction', "No strong TM-like region detected")

    style.field('Method', "Sequence-only heuristic")

    style.section('SIGNAL PEPTIDE SCREENING')

    signal = analysis[
        "signal_peptide"
    ]

    if signal[
        "predicted"
    ]:
        style.field('Prediction', "Possible N-terminal signal peptide")

        style.field('Heuristic evidence', f"{signal['evidence']}")
    else:
        style.field('Prediction', "No convincing pattern")

    style.field('Method', "Sequence-only heuristic")

    style.section('FLEXIBILITY SCREENING')

    flexibility = analysis[
        "flexibility"
    ]

    if flexibility:
        scores = [
            entry["score"]
            for entry in flexibility
        ]

        style.field('Average', f"{sum(scores) / len(scores):.3f}")

        style.field('Minimum', f"{min(scores):.3f}")

        style.field('Maximum', f"{max(scores):.3f}")
    else:
        print(
            "Sequence too short for "
            "window-based analysis."
        )

    print()
    print(
        "Structural results are sequence-level "
        "screening outputs, not experimental "
        "structure determinations."
    )


# ============================================================
# 5 - COMPLETE PROTEIN REPORT
# ============================================================


def _print_complete_protein_report(
    analysis,
    source="Selected ORF",
    orf_id=None,
    orf_length=None,
    strand=None,
    reading_frame=None,
    orf_status=None,
):
    style.header(
        "COMPLETE PROTEIN REPORT"
    )

    style.section('1. INPUT')

    style.field('Source', f"{source}")

    if orf_id is not None:
        style.field('ORF', f"{orf_id}")

    if orf_length is not None:
        style.field('ORF length', f"{orf_length} bp")

    if strand is not None:
        style.field('Strand', f"{strand}")

    if reading_frame is not None:
        style.field('Reading frame', f"{_frame_display(reading_frame)}")

    if orf_status is not None:
        style.field('ORF status', f"{orf_status}")

    style.field('Protein length', f"{analysis['length']} aa")

    style.section('2. TRANSLATION')

    translation = analysis[
        "translation"
    ]

    if translation is None:
        print(
            "Translation context: Not supplied"
        )
    else:
        style.field('Reading frame', f"{_frame_display(translation['reading_frame'])}")

        style.field('Start codon', f"{translation['start_codon'] or 'Not detected'}")

        style.field('Stop codon', f"{translation['stop_codon'] if translation['stop_found'] else 'Not detected'}")

        style.field('Complete endpoint', f"{'Yes' if translation['stop_found'] else 'No'}")

        style.field('Translation status', f"{'Complete' if translation['stop_found'] else 'Partial'}")

        print(
            f"Trailing nucleotides:"
            f" {translation['trailing_nucleotides']}"
        )

    style.field('Protein sequence', f"{analysis['protein']}")

    style.section('3. SEQUENCE QUALITY')

    _protein_quality_block(
        analysis
    )

    style.section('4. COMPOSITION')

    _protein_composition_block(
        analysis
    )

    counts = analysis[
        "composition"
    ][
        "counts"
    ]

    if counts:
        residue, count = max(
            counts.items(),
            key=lambda item: item[1]
        )

        style.field('Most abundant residue', f"{residue} "
            f"({AMINO_ACIDS[residue]['name']})")

        style.field('Occurrences', f"{count}")

    style.section('5. PHYSICOCHEMISTRY')

    _protein_key_properties(
        analysis
    )

    style.field('Hydropathy', f"{_hydropathy_interpretation(analysis['gravy'])}")

    style.field('Charge', f"{_charge_interpretation(analysis)}")

    style.field('Stability', f"{_stability_interpretation(analysis)}")

    style.section('6. STRUCTURAL SCREENING')

    secondary = analysis[
        "secondary_structure"
    ]

    style.field('Helix propensity', f"{secondary['helix']:.3f}")

    style.field('Turn propensity', f"{secondary['turn']:.3f}")

    style.field('Sheet propensity', f"{secondary['sheet']:.3f}")

    style.field('TM-like regions', f"{len(analysis['transmembrane_domains'])}")

    signal = analysis[
        "signal_peptide"
    ]

    style.field('Signal peptide', f"{'Possible' if signal['predicted'] else 'Not detected'}")

    style.field('Structural method', "Sequence-level heuristics")

    style.section('7. EXTINCTION COEFFICIENT')

    extinction, reduced, oxidized = (
        _safe_extinction_values(
            analysis
        )
    )

    style.field('Wavelength', "280 nm")

    style.field('Tryptophan', f"{extinction['tryptophan']}")

    style.field('Tyrosine', f"{extinction['tyrosine']}")

    style.field('Cysteine', f"{extinction['cysteine']}")

    style.field('All Cys reduced', f"{reduced} M^-1 cm^-1")

    style.field('Max. disulfide form', f"{oxidized} M^-1 cm^-1")

    style.field('Disulfide pairs', f"{extinction['disulfide_pairs']}")

    if extinction["unpaired_cysteine"]:
        style.field('Unpaired cysteine', "1")

    style.section('8. INTERPRETATION')

    level, reason = _sequence_applicability(
        analysis["length"]
    )

    style.field('Applicability level', f"{level}")

    style.field('Reason', f"{reason}")

    print()
    print(
        "All predictions are sequence-level "
        "descriptors or heuristics."
    )

    print(
        "They do not establish biological function, "
        "experimental stability, 3D structure, "
        "cellular localization or biological activity."
    )

    print()
    print(major_rule())


# ============================================================
# MENU
# ============================================================


def _protein_ui_menu():

    style.header('PROTEIN ANALYSIS')

    style.menu([
        ("Overview", [('1', 'Protein Overview'), ('2', 'Translation & Sequence')]),
        ("Analysis", [('3', 'Amino Acid Composition'),
                      ('4', 'Physicochemical & Structural Properties')]),
        ("Report", [('5', 'Complete Protein Report')]),
    ])


# ============================================================
# MAIN PROTEIN UI
# ============================================================


def protein_analysis_menu(
    protein,
    rna=None,
    source="Selected ORF",
    orf_id=None,
    orf_length=None,
    strand=None,
    reading_frame=None,
    orf_status=None,
    initiator_start=False,
    genetic_code=1,
):
    """
    Main Protein Analysis UI.

    If RNA is supplied, it is used for translation.
    There is no complete-DNA fallback.
    """

    if rna is not None:

        translation = analyze_translation(
            rna,
            initiator_start=initiator_start,
            genetic_code=genetic_code,
        )

        # Preserve the biological frame from the selected ORF.
        if reading_frame is not None:
            translation[
                "reading_frame"
            ] = reading_frame

        active_protein = translation[
            "protein"
        ]

    else:

        translation = None
        active_protein = protein

    if not isinstance(
        active_protein,
        str
    ):
        active_protein = ""

    active_protein = (
        active_protein
        .upper()
    )

    if not active_protein:

        style.header(
            "PROTEIN ANALYSIS ERROR"
        )

        style.note("No resolvable protein sequence is available.")

        _protein_ui_pause()

        return

    analysis = analyze_complete_protein(
        active_protein,
        rna=None,
    )

    # Recalculate molecular weight with the actual translation
    # records so ambiguous codons can contribute a mass range.
    if translation is not None:
        analysis[
            "molecular_weight"
        ] = calculate_protein_molecular_weight(
            active_protein,
            translation
        )

    analysis[
        "translation"
    ] = translation

    while True:

        _protein_ui_menu()

        selection = style.ask_choice()

        if selection == "0":
            return

        elif selection == "1":

            _print_protein_overview(
                analysis,
                source=source,
                orf_id=orf_id,
                orf_length=orf_length,
                strand=strand,
                reading_frame=reading_frame,
                orf_status=orf_status,
            )

            _protein_ui_pause()

        elif selection == "2":

            _print_translation_and_sequence(
                analysis
            )

            _protein_ui_pause()

        elif selection == "3":

            _print_amino_acid_composition(
                analysis
            )

            _protein_ui_pause()

        elif selection == "4":

            _print_physicochemical_structural(
                analysis
            )

            _protein_ui_pause()

        elif selection == "5":

            _print_complete_protein_report(
                analysis,
                source=source,
                orf_id=orf_id,
                orf_length=orf_length,
                strand=strand,
                reading_frame=reading_frame,
                orf_status=orf_status,
            )

            _protein_ui_pause()

        else:

            style.invalid_selection()


# ============================================================
# PIPELINE ADAPTER
# ============================================================


def run_protein_analysis_from_record(
    protein,
    name="Protein record",
    pipeline=None,
):
    """
    Start Protein Analysis from a protein sequence (a protein FASTA or
    GenPept record; audit D14).

    The sequence must be in the IUPAC protein alphabet. When a pipeline is
    given, the complete protein analysis is stored before the interactive
    menu opens.
    """

    from core.alphabets import PROTEIN as _PROTEIN, require_sequence_type

    require_sequence_type(
        protein,
        (_PROTEIN,),
        declared=_PROTEIN,
        context="Protein analysis",
    )

    clean_protein = "".join(
        protein.upper().split()
    )

    if pipeline is not None:
        context = getattr(pipeline, "context", None)
        analysis = analyze_complete_protein(clean_protein)
        if context is not None:
            context.protein_source = "Protein record"
            context.store_protein_analysis(analysis, protein_sequence=clean_protein)

    protein_analysis_menu(
        protein=clean_protein,
        source="Protein record",
        orf_id=name,
    )


def run_protein_analysis_from_orf(
    orf_sequence,
    orf_id="Selected ORF",
    strand=None,
    reading_frame=None,
    orf_status=None,
    pipeline=None,
    genetic_code=1,
):
    """
    Start Protein Analysis from exactly one selected ORF.

    When a pipeline is given, the complete protein analysis (computed with the
    initiator-Met rule) is stored before the interactive menu opens.

    The complete DNA sequence is never used as a fallback.
    """

    if not isinstance(
        orf_sequence,
        str
    ):

        print(
            "Protein analysis error: "
            "selected ORF sequence must be a string."
        )

        return

    clean_orf = "".join(
        orf_sequence
        .upper()
        .split()
    )

    allowed = set(
        RNA_IUPAC
    ) | {"T"}

    if (
        not clean_orf
        or any(
            base not in allowed
            for base in clean_orf
        )
    ):

        print(
            "Protein analysis error: "
            "selected ORF contains invalid "
            "nucleotide symbols."
        )

        return

    rna = clean_orf.replace(
        "T",
        "U"
    )

    # Every Synapse ORF begins at its start codon, so the first codon is
    # decoded as initiator Met (also for GTG/TTG starts).
    translation = analyze_translation(
        rna,
        initiator_start=True,
        genetic_code=genetic_code,
    )

    if reading_frame is not None:

        translation[
            "reading_frame"
        ] = reading_frame

    protein = translation[
        "protein"
    ]

    if not protein:

        print(
            "Protein analysis error: "
            "selected ORF produced no "
            "resolvable protein sequence."
        )

        return

    if pipeline is not None:
        context = getattr(pipeline, "context", None)
        try:
            analysis = analyze_complete_protein("", rna=rna, initiator_start=True, genetic_code=genetic_code)
            if context is not None:
                context.protein_source = "Selected ORF"
                context.store_protein_analysis(analysis, protein_sequence=analysis.get("protein", ""))
        except Exception:
            import logging
            logging.getLogger("synapse.protein").exception(
                "Could not store protein result in the pipeline."
            )

    protein_analysis_menu(
        protein=protein,
        rna=rna,
        source="Selected ORF",
        orf_id=orf_id,
        orf_length=len(clean_orf),
        strand=strand,
        reading_frame=(
            reading_frame
            if reading_frame is not None
            else "+1"
        ),
        orf_status=orf_status,
        initiator_start=True,
        genetic_code=genetic_code,
    )