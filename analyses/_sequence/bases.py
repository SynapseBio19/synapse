"""DNA base tables, normalization, nucleotide counts and complements."""

from __future__ import annotations



# ============================================================
# SEQUENCE ANALYSIS
# ============================================================


# ============================================================
# CONSTANTS
# ============================================================

DNA_STANDARD_BASES = {
    "A",
    "T",
    "C",
    "G"
}


DNA_AMBIGUOUS_CODES = {
    "R": "A/G",
    "Y": "C/T",
    "S": "G/C",
    "W": "A/T",
    "K": "G/T",
    "M": "A/C",
    "B": "C/G/T",
    "D": "A/G/T",
    "H": "A/C/T",
    "V": "A/C/G",
    "N": "Any"
}


DNA_IUPAC_BASES = {
    "A": {"A"},
    "T": {"T"},
    "C": {"C"},
    "G": {"G"},
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
    "N": {"A", "T", "C", "G"}
}


DNA_COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N"
}


# Average nucleotide residue masses for ssDNA.
# Used compositionally together with complementary bases
# for an estimated dsDNA molecular weight.

DNA_NUCLEOTIDE_RESIDUE_MASSES = {
    "A": 313.21,
    "T": 304.20,
    "C": 289.18,
    "G": 329.21
}


# Approximate average molecular weight of one base pair.

AVERAGE_DSDNA_BP_MASS = 660.0


# ============================================================
# INTERNAL HELPERS
# ============================================================


# ============================================================
# INPUT NORMALIZATION
# ============================================================

def normalize_dna_sequence(dna):
    """
    Normalize sequence input for sequence-level analysis.

    Whitespace is ignored and letters are converted to uppercase.
    Symbols that are not DNA/IUPAC symbols are preserved so that
    the quality analysis can report them as invalid rather than
    silently deleting them.
    """

    if dna is None:
        return ""

    return "".join(
        str(dna).split()
    ).upper()




def _is_standard_base(base):

    return base in DNA_STANDARD_BASES


def _is_iupac_base(base):

    return base in DNA_IUPAC_BASES


def _get_standard_bases(dna):

    return [
        base
        for base in dna
        if base in DNA_STANDARD_BASES
    ]


def _get_valid_iupac_bases(dna):

    return [
        base
        for base in dna
        if base in DNA_IUPAC_BASES
    ]


# ============================================================
# NUCLEOTIDE COUNTS
# ============================================================


def count_nucleotides(dna):
    """Counts of A, T, C and G; U (the RNA form of T) is counted as T.

    Counting U keeps GC/AT content and skews correct for an RNA sequence
    passed to this DNA module directly (audit D14-F10).
    """

    nucleotide_counts = {
        "A": 0,
        "T": 0,
        "C": 0,
        "G": 0
    }

    for nucleotide in dna:

        if nucleotide == "U":
            nucleotide = "T"

        if nucleotide in nucleotide_counts:

            nucleotide_counts[
                nucleotide
            ] += 1

    return nucleotide_counts


# ============================================================
# AMBIGUOUS BASE COUNTS
# ============================================================


def count_ambiguous_bases(dna):

    counts = {}

    for code in DNA_AMBIGUOUS_CODES:

        count = dna.count(code)

        if count > 0:

            counts[code] = count

    return counts


# ============================================================
# COMPLEMENT
# ============================================================


def get_complement(dna):

    complement = []

    for base in dna:

        complement.append(
            DNA_COMPLEMENT.get(
                base,
                base
            )
        )

    return "".join(complement)


def get_reverse_complement(dna):

    return get_complement(dna)[::-1]
