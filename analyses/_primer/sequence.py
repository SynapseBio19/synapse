"""DNA sequence helpers shared by every part of the primer module."""

from __future__ import annotations

from core.seqtools import IUPAC_COMPLEMENT as _CORE_IUPAC_COMPLEMENT
from functools import lru_cache
from operator import eq as _equal


DNA_ALPHABET = set("ACGT")
IUPAC_ALPHABET = set("ACGTRYSWKMBDHVN")

# IUPAC complement table shared via core.seqtools.
IUPAC_COMPLEMENT = dict(_CORE_IUPAC_COMPLEMENT)


# ============================================================================
# VALIDATION
# ============================================================================

def normalize_dna_sequence(sequence: str, allow_iupac: bool = False) -> str:
    """Normalize DNA by removing whitespace, uppercasing and converting U to T."""
    if not isinstance(sequence, str):
        raise TypeError("DNA sequence must be a string.")
    return _normalize_dna_cached(sequence, allow_iupac)


# Primer design normalizes the same short strings about two million times
# for a 3 kb template (P2-6). The result depends only on the string and the
# alphabet flag, so it is cached; an invalid input is never cached and
# raises the same error every time.
@lru_cache(maxsize=65536)
def _normalize_dna_cached(sequence: str, allow_iupac: bool) -> str:
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    if not normalized:
        raise ValueError("DNA sequence is empty.")
    alphabet = IUPAC_ALPHABET if allow_iupac else DNA_ALPHABET
    invalid = sorted(set(normalized) - alphabet)
    if invalid:
        raise ValueError("Invalid DNA character(s): " + ", ".join(invalid))
    return normalized


def validate_target_region(sequence: str, target_start: int, target_end: int) -> None:
    """Validate a 0-based, half-open target interval."""
    if target_start < 0 or target_end < 0:
        raise ValueError("Target coordinates must be non-negative.")
    if target_start >= target_end:
        raise ValueError("target_start must be smaller than target_end.")
    if target_end > len(sequence):
        raise ValueError("Target region extends beyond the sequence.")


# ============================================================================
# BASIC SEQUENCE ANALYSIS
# ============================================================================

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement."""
    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    # Every IUPAC symbol has a complement in the table, so translating and
    # reversing gives the same string as a per-base lookup.
    return sequence.translate(_IUPAC_COMPLEMENT_TABLE)[::-1]


_IUPAC_COMPLEMENT_TABLE = str.maketrans(IUPAC_COMPLEMENT)


def gc_percent(sequence: str) -> float:
    """Return GC content as a percentage."""
    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    gc = sequence.count("G") + sequence.count("C")
    return 100.0 * gc / len(sequence)


def count_gc_clamp(sequence: str) -> int:
    """Count consecutive G/C bases at the 3' terminus."""
    sequence = normalize_dna_sequence(sequence)
    count = 0
    for base in reversed(sequence):
        if base in "GC":
            count += 1
        else:
            break
    return count


def longest_homopolymer(sequence: str) -> int:
    """Return the longest consecutive identical-base run."""
    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    best = 1
    current = 1
    for previous, base in zip(sequence, sequence[1:]):
        current = current + 1 if previous == base else 1
        best = max(best, current)
    return best


def longest_tandem_repeat(sequence: str, max_unit: int = 3) -> int:
    """Return the largest number of consecutive copies of a short repeat unit."""
    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    if len(sequence) < 2:
        return 1
    best = 1
    max_unit = min(max_unit, len(sequence) // 2)
    for unit_length in range(1, max_unit + 1):
        # Copy k + 1 of a unit starting at s repeats copy k exactly when
        # sequence[p] == sequence[p + unit_length] for the unit_length
        # positions of copy k. A run of m such matches from s therefore
        # holds 1 + m // unit_length copies, and the longest run gives the
        # most copies for this unit length (P2-6: the same result as
        # comparing the copies one by one).
        matches = bytes(map(_equal, sequence, sequence[unit_length:]))
        longest_run = max(map(len, matches.split(b"\x00")))
        best = max(best, 1 + longest_run // unit_length)
    return best
