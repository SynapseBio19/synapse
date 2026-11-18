"""Sequence-level hairpin, self-dimer and heterodimer screens (heuristics,
not free-energy calculations)."""

from __future__ import annotations

from typing import Optional

from analyses._primer.sequence import normalize_dna_sequence, reverse_complement
from analyses._primer.models import StructureResult


# ============================================================================
# SECONDARY STRUCTURE HEURISTICS
# ============================================================================

def _alignment_stem(left: str, right: str) -> tuple[int, Optional[int]]:
    """Find the longest contiguous complementary alignment."""
    left = normalize_dna_sequence(left)
    right = normalize_dna_sequence(right)
    right_rc = reverse_complement(right)
    # The stem is the longest common substring of left and the reverse
    # complement of right. The position is the one the original diagonal
    # scan reported: it walked the diagonals shift = i - j in ascending
    # order, positions within one in ascending order, and kept the first
    # run to reach a new maximum. Among all windows of the maximal length
    # that is the one with the smallest shift, then the smallest start i.
    # A common window of length k implies one of every shorter length, so
    # the maximal length is found by bisection (P2-6).
    left_length = len(left)
    right_length = len(right_rc)

    def has_common_window(size: int) -> bool:
        windows = {right_rc[j:j + size] for j in range(right_length - size + 1)}
        return any(left[i:i + size] in windows for i in range(left_length - size + 1))

    low, high = 0, min(left_length, right_length)
    while low < high:
        middle = (low + high + 1) // 2
        if has_common_window(middle):
            low = middle
        else:
            high = middle - 1
    if low == 0:
        return 0, None

    last_start: dict[str, int] = {}
    for j in range(right_length - low + 1):
        last_start[right_rc[j:j + low]] = j
    best_key: Optional[tuple[int, int]] = None
    for i in range(left_length - low + 1):
        j = last_start.get(left[i:i + low])
        if j is not None and (best_key is None or (i - j, i) < best_key):
            best_key = (i - j, i)
    assert best_key is not None
    return low, best_key[1]


def _three_prime_complementarity(
    first: str,
    second: Optional[str] = None,
    window: int = 8,
) -> int:
    """Estimate contiguous terminal complementarity as a screening heuristic."""
    first = normalize_dna_sequence(first)
    second = first if second is None else normalize_dna_sequence(second)
    first_terminal = first[-window:]
    second_terminal = second[-window:]
    # reverse_complement(second_terminal[:length]) is the last `length`
    # bases of the reverse complement of second_terminal, so it is
    # computed once.
    second_reverse = reverse_complement(second_terminal)
    best = 0
    for length in range(1, min(len(first_terminal), len(second_terminal)) + 1):
        first_part = first_terminal[-length:]
        if first_part == second_reverse[len(second_reverse) - length:]:
            best = length
    return best


# ============================================================================
# V8 STRUCTURE CLASSIFICATION + HARD PAIR SAFETY GATES
# ============================================================================
#
# V7 exposed a real classification bug:
# the legacy structure functions used max_stem=4 as the severity threshold,
# which made a 20 bp or 24 bp stem appear as LOW in some reports.
#
# V8 fixes the structure severity itself and fixes the pair-ranking bug where
# the heterodimer length was read from a non-existent "stem_length" attribute.
# ============================================================================

def _v8_structure_severity(max_stem: int, three_prime_stem: int = 0) -> str:
    """Classify structural risk from the measured stem, not the old default."""
    max_stem = int(max_stem or 0)
    three_prime_stem = int(three_prime_stem or 0)

    if max_stem >= 10 or three_prime_stem >= 4:
        return "VERY HIGH"
    if max_stem >= 7 or three_prime_stem >= 3:
        return "HIGH"
    if max_stem >= 4 or three_prime_stem >= 2:
        return "MODERATE"
    return "LOW"


def analyze_hairpin(sequence: str, max_stem: int = 4) -> StructureResult:
    """Analyze hairpins using the existing detection algorithm and V8 severity."""
    sequence = normalize_dna_sequence(sequence)

    best = 0
    best_position = None
    length = len(sequence)
    max_test_stem = min(8, length // 2)
    # The reverse complement of sequence[r:r + stem] is the slice
    # [length - r - stem:length - r] of the whole reverse complement, so
    # it is computed once (P2-6: this loop dominated primer design).
    reverse = reverse_complement(sequence) if length else ""

    # Stems are tried from the longest down, and only the first match at
    # the longest matching stem sets best and best_position, so the
    # search ends there.
    #
    # A right arm starting at right_start pairs with the left arm when
    # left == reverse[k:k + stem] with k = length - right_start - stem.
    # right_start runs from left_start + stem + 3 to length - stem, so k
    # runs from 0 to length - left_start - 2 * stem - 3, and any such
    # occurrence ends before length - left_start - stem - 3. One bounded
    # str.find replaces the loop over right_start (P2-6).
    for stem in range(max_test_stem, 2, -1):
        for left_start in range(0, length - stem - 3):
            left = sequence[left_start:left_start + stem]
            if reverse.find(left, 0, length - left_start - stem - 3) != -1:
                best = stem
                best_position = left_start
                break
        if best:
            break

    result = StructureResult(
        max_stem=best,
        three_prime_stem=0,
        best_position=best_position,
        severity=_v8_structure_severity(best),
        explanation=(
            f"Longest detected intramolecular stem: {best} bp; "
            "3-nt minimum-loop heuristic."
        ),
    )
    return result


def analyze_self_dimer(sequence: str, max_stem: int = 4) -> StructureResult:
    """Analyze self-dimer formation with V8 severity thresholds."""
    sequence = normalize_dna_sequence(sequence)
    stem, position = _alignment_stem(sequence, sequence)
    terminal = _three_prime_complementarity(sequence)

    return StructureResult(
        max_stem=stem,
        three_prime_stem=terminal,
        best_position=position,
        severity=_v8_structure_severity(stem, terminal),
        explanation=(
            f"Longest self-complementary stem: {stem} bp; "
            f"estimated 3'-complementarity: {terminal} bp."
        ),
    )


def analyze_heterodimer(
    forward: str,
    reverse: str,
    max_stem: int = 4,
) -> StructureResult:
    """Analyze heterodimer formation with V8 severity thresholds."""
    forward = normalize_dna_sequence(forward)
    reverse = normalize_dna_sequence(reverse)

    stem, position = _alignment_stem(forward, reverse)
    terminal = max(
        _three_prime_complementarity(forward, reverse),
        _three_prime_complementarity(reverse, forward),
    )

    return StructureResult(
        max_stem=stem,
        three_prime_stem=terminal,
        best_position=position,
        severity=_v8_structure_severity(stem, terminal),
        explanation=(
            f"Longest heterodimer stem: {stem} bp; "
            f"estimated 3'-complementarity: {terminal} bp."
        ),
    )
