"""Prefix counts for sliding-window scans over long sequences.

A window scan that slices and recounts every window costs
O(length x window). With prefix counts, the number of a symbol in any
window is one subtraction, so a scan costs O(length). Used by the DNA and RNA window scans (analyses/_sequence,
analyses/_rna). The counts are exact
integers, and the callers apply the same formulas to them as before, so
the results are identical; only the time changes (a 4.6 Mb genome:
minutes before, seconds now; docs/BENCHMARKS.md).
"""

from __future__ import annotations

from array import array
from itertools import accumulate
from operator import add
from typing import Collection, Dict, Iterable

DNA_STANDARD = frozenset("ACGT")


def _prefix(flags: Iterable[int]) -> array:
    """prefix[i] = number of true flags among the first i positions."""
    counts = array("q", [0])
    counts.extend(accumulate(flags))
    return counts


def base_prefix_counts(sequence: str, bases: Iterable[str]) -> Dict[str, array]:
    """Prefix count of each base in ``bases`` (case-sensitive, like str.count)."""
    return {base: _prefix(map(base.__eq__, sequence)) for base in bases}


def nonstandard_prefix_count(sequence: str, standard: Collection[str] = DNA_STANDARD) -> array:
    """Prefix count of positions whose symbol is not in ``standard`` (default A, C, G, T)."""
    member = frozenset(standard).__contains__
    return _prefix(map(lambda base: not member(base), sequence))


def dinucleotide_prefix_count(dna: str, dinucleotide: str) -> array:
    """prefix[i] = occurrences of ``dinucleotide`` starting before position i.

    A window [s, e) holds ``prefix[e - 1] - prefix[s]`` occurrences, the
    same number as ``dna[s:e].count(dinucleotide)`` for a dinucleotide
    that cannot overlap itself (such as "CG").
    """
    return _prefix(map(dinucleotide.__eq__, map(add, dna, dna[1:])))


# A profile with more windows than this is stored thinned by analyze_sequence
# and analyze_rna (the DNA statistics still use every window). 200,000 windows of a 100-nt
# profile cover a 200 kb sequence at step 1; a 4.6 Mb genome at step 1
# needed 3.5 GB for the two full profiles (docs/BENCHMARKS.md).
DEFAULT_MAX_PROFILE_POINTS = 200_000


def window_count(length, window_size, step):
    """Number of windows a profile of a sequence of ``length`` has."""
    if window_size <= 0 or step <= 0 or length < window_size:
        return 0
    return (length - window_size) // step + 1
