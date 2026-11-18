"""Alphabets, IUPAC handling, normalization and translation."""

from __future__ import annotations

import re
from core.seqtools import codon_dict, IUPAC_COMPLEMENT as _CORE_IUPAC_COMPLEMENT
from typing import Dict, Iterable, List, Optional, Set

from analyses._mutation.models import MutationAnalysisError


# ============================================================================
# CONSTANTS
# ============================================================================

DNA_ALPHABET = set("ACGTN")

IUPAC_CODES: Dict[str, Set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

IUPAC_FROM_BASES = {
    frozenset(value): key
    for key, value in IUPAC_CODES.items()
}

# IUPAC complement table shared via core.seqtools.
IUPAC_COMPLEMENT = dict(_CORE_IUPAC_COMPLEMENT)


# Standard code (NCBI table 1) from core.seqtools.
CODON_TABLE = codon_dict(1)


AA_NAMES = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "Stop",
    "X": "Unknown",
}


PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}


# ============================================================================
# BASIC SEQUENCE FUNCTIONS
# ============================================================================

def normalize_dna(
    sequence: str,
    *,
    allow_n: bool = True,
) -> str:

    if sequence is None:
        raise MutationAnalysisError(
            "DNA sequence cannot be None."
        )

    dna = re.sub(
        r"\s+",
        "",
        str(sequence),
    ).upper()

    if not dna:
        raise MutationAnalysisError(
            "DNA sequence cannot be empty."
        )

    alphabet = (
        DNA_ALPHABET
        if allow_n
        else set("ACGT")
    )

    invalid = sorted(
        set(dna) - alphabet
    )

    if invalid:
        raise MutationAnalysisError(
            "Invalid DNA character(s): "
            + ", ".join(invalid)
        )

    return dna


def normalize_iupac_sequence(
    sequence: str,
) -> str:

    if sequence is None:
        raise MutationAnalysisError(
            "Sequence cannot be None."
        )

    normalized = re.sub(
        r"\s+",
        "",
        str(sequence),
    ).upper()

    if not normalized:
        raise MutationAnalysisError(
            "Sequence cannot be empty."
        )

    invalid = sorted(
        set(normalized) - set(IUPAC_CODES)
    )

    if invalid:
        raise MutationAnalysisError(
            "Invalid IUPAC nucleotide symbol(s): "
            + ", ".join(invalid)
        )

    return normalized


def validate_dna_sequence(
    sequence: str,
    *,
    allow_n: bool = True,
) -> bool:

    normalize_dna(
        sequence,
        allow_n=allow_n,
    )

    return True


def gc_content(sequence: str) -> float:

    dna = normalize_dna(sequence)

    return (
        100.0
        * (
            dna.count("G")
            + dna.count("C")
        )
        / len(dna)
    )


def reverse_complement(
    sequence: str,
) -> str:

    dna = normalize_dna(sequence)

    table = str.maketrans(
        "ACGTN",
        "TGCAN",
    )

    return dna.translate(table)[::-1]


def reverse_complement_iupac(
    sequence: str,
) -> str:

    dna = normalize_iupac_sequence(
        sequence
    )

    return "".join(
        IUPAC_COMPLEMENT[b]
        for b in reversed(dna)
    )


def iupac_matches(
    code_a: str,
    code_b: str,
) -> bool:

    a = IUPAC_CODES.get(
        str(code_a).upper()
    )

    b = IUPAC_CODES.get(
        str(code_b).upper()
    )

    if not a or not b:
        return False

    return bool(a & b)


def iupac_code(
    bases: Iterable[str],
) -> str:

    normalized = frozenset(
        str(base).upper()
        for base in bases
    )

    if not normalized:
        raise MutationAnalysisError(
            "Cannot create an IUPAC code from no bases."
        )

    try:
        return IUPAC_FROM_BASES[
            normalized
        ]

    except KeyError as exc:

        raise MutationAnalysisError(
            "Unsupported base set for IUPAC encoding: "
            + str(sorted(normalized))
        ) from exc


# ============================================================================
# TRANSLATION
# ============================================================================

def translate_dna(
    sequence: str,
    frame: int = 0,
    stop_at_stop: bool = False,
) -> str:

    dna = normalize_iupac_sequence(
        sequence
    )

    if frame not in (0, 1, 2):
        raise MutationAnalysisError(
            "Frame must be 0, 1, or 2."
        )

    protein: List[str] = []

    for i in range(
        frame,
        len(dna) - 2,
        3,
    ):

        codon = dna[i:i + 3]

        if any(
            base not in "ACGT"
            for base in codon
        ):
            aa = "X"

        else:
            aa = CODON_TABLE.get(
                codon,
                "X",
            )

        protein.append(aa)

        if stop_at_stop and aa == "*":
            break

    return "".join(protein)


def amino_acid_name(
    symbol: Optional[str],
) -> Optional[str]:

    if symbol is None:
        return None

    return AA_NAMES.get(
        symbol,
        "Unknown",
    )
