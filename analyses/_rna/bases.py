"""RNA base tables, codon table, transcription, base counts and complements."""

from __future__ import annotations

from core.seqtools import codon_dict


# ============================================================
# RNA ANALYSIS
# ============================================================


# ============================================================
# CONSTANTS
# ============================================================

RNA_STANDARD_BASES = {
    "A",
    "U",
    "C",
    "G"
}


RNA_IUPAC_BASES = {
    "A": {"A"},
    "U": {"U"},
    "C": {"C"},
    "G": {"G"},

    "R": {"A", "G"},
    "Y": {"C", "U"},
    "S": {"G", "C"},
    "W": {"A", "U"},
    "K": {"G", "U"},
    "M": {"A", "C"},

    "B": {"C", "G", "U"},
    "D": {"A", "G", "U"},
    "H": {"A", "C", "U"},
    "V": {"A", "C", "G"},
    "N": {"A", "U", "C", "G"}
}


RNA_AMBIGUOUS_CODES = {
    "R": "A/G",
    "Y": "C/U",
    "S": "G/C",
    "W": "A/U",
    "K": "G/U",
    "M": "A/C",
    "B": "C/G/U",
    "D": "A/G/U",
    "H": "A/C/U",
    "V": "A/C/G",
    "N": "Any"
}


RNA_COMPLEMENT = {
    "A": "U",
    "U": "A",
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


RNA_NUCLEOTIDE_RESIDUE_MASSES = {
    "A": 329.21,
    "U": 306.17,
    "C": 305.18,
    "G": 345.21
}


# Standard code (NCBI table 1) from core.seqtools, RNA alphabet.
RNA_CODON_TABLE = codon_dict(1, rna=True, stop_symbol="STOP")


def rna_codon_table(genetic_code=1):
    """Codon table (RNA alphabet, "STOP" for stops) for an NCBI genetic code.

    Audit P1-9: the three-frame codon view used to read table 1 even when a
    selected ORF had been found with another table, so the preview could
    disagree with the ORF's own translation.
    """
    if genetic_code in (None, 1):
        return RNA_CODON_TABLE
    return codon_dict(int(genetic_code), rna=True, stop_symbol="STOP")


def orf_genetic_code(selected_orf, default=1):
    """The genetic code of a selected ORF, or *default* when it carries none."""
    if selected_orf is None:
        return default
    value = _orf_value(selected_orf, "genetic_code", None)
    try:
        return int(value) if value else default
    except (TypeError, ValueError):
        return default


# ============================================================
# TRANSCRIPTION
# ============================================================

def transcribe_rna(dna):
    """Transcribe a coding-strand DNA sequence to RNA (T -> U).

    Whitespace is removed and letters are uppercased; a sequence that is
    already RNA is returned unchanged. Symbols outside the IUPAC nucleotide
    alphabet raise ValueError (audit D14-F9).
    """

    from core.alphabets import DNA_IUPAC, RNA_IUPAC, normalize_sequence, symbols_outside

    sequence = normalize_sequence(dna)
    invalid = symbols_outside(sequence, DNA_IUPAC | RNA_IUPAC)
    if invalid:
        shown = ", ".join(f"{symbol!r} at {position}" for position, symbol in invalid[:5])
        raise ValueError(
            f"Transcription needs a nucleotide sequence; found {shown}"
            + (" ..." if len(invalid) > 5 else "") + "."
        )

    return sequence.replace(
        "T",
        "U"
    )


# ============================================================
# BASIC BASE ANALYSIS
# ============================================================

def count_rna_bases(rna):

    counts = {
        "A": 0,
        "U": 0,
        "C": 0,
        "G": 0
    }

    for base in rna:

        if base in counts:
            counts[base] += 1

    return counts


# ============================================================
# RNA COMPLEMENT
# ============================================================

def get_rna_complement(rna: str) -> str:

    return "".join(
        RNA_COMPLEMENT.get(
            base,
            base
        )
        for base in rna
    )


def get_rna_reverse_complement(rna):

    return get_rna_complement(
        rna
    )[::-1]


# ============================================================
# RNA REPORT
# ============================================================


# ============================================================
# SYNAPSE PIPELINE INTEGRATION
# ============================================================
#
# This section is the ONLY pipeline/UI layer in this module.
#
# Input resolution:
#   1. Explicitly supplied selected ORF
#   2. ORF stored in the pipeline context
#   3. Complete DNA sequence
#
# The original RNA analysis engine above is kept intact.
# The functions below only connect that engine to Synapse.
# ============================================================


def _orf_value(orf, key, default=None):
    """Read a value from either a dict-like ORF or an ORF object."""

    if orf is None:
        return default

    if isinstance(orf, dict):
        return orf.get(key, default)

    return getattr(orf, key, default)
