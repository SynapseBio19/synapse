"""DNA input checks used by the main menu.

Thin layer over ``core.alphabets`` (audit D13): ``normalize_dna`` has the
same contract as ``core.seqtools.normalize_dna`` and
``core.alphabets.normalize_sequence`` - whitespace removed, uppercased,
nothing converted or dropped.
"""

from core.alphabets import DNA_IUPAC, normalize_sequence, symbols_outside


def normalize_dna(dna):

    return normalize_sequence(dna)

def validate_dna(dna):
    """True for a non-empty sequence of IUPAC DNA symbols."""

    return bool(dna) and not symbols_outside(dna, DNA_IUPAC)

def get_dna_errors(dna):

    return [
        {"character": nucleotide, "position": position}
        for position, nucleotide in symbols_outside(dna, DNA_IUPAC)
    ]
