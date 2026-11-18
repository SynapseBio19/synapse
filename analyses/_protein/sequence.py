"""Protein sequence validation, resolved and positional views, sequence quality."""

from __future__ import annotations

from analyses._protein.constants import (
    NONSTANDARD_RESIDUES,
    STANDARD_AMINO_ACIDS,
    STOP_SYMBOL,
    UNRESOLVED_RESIDUES,
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
