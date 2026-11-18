"""Needleman-Wunsch global alignment (banded when exact), identity, edit
distance and alignment blocks."""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from analyses._mutation.sequence import iupac_matches, normalize_iupac_sequence
from analyses._mutation.models import AlignmentBlock, MutationAnalysisError


# ============================================================================
# ALIGNMENT
# ============================================================================

def _nucleotide_match(
    a: str,
    b: str,
) -> bool:

    return iupac_matches(
        a,
        b,
    )


# Alignments with more cells than this try an exact banded computation first.
BANDED_ALIGNMENT_MIN_CELLS = 1_000_000


def _banded_align(ref, var, half_width, substitution, gap_score):
    """Global alignment restricted to |j - i| <= half_width; None if (n, m) is outside."""
    rows, cols = len(ref), len(var)
    if abs(cols - rows) > half_width:
        return None
    width = 2 * half_width + 1
    negative = float("-inf")
    diagonal_step, up_step, left_step = 1, 2, 3
    trace = []
    previous = [negative] * width
    row = bytearray(width)
    for k in range(width):                     # row 0: j = k - half_width
        j = k - half_width
        if 0 <= j <= cols:
            previous[k] = j * gap_score
            row[k] = left_step if j > 0 else 0
    trace.append(row)
    for i in range(1, rows + 1):
        scores = substitution[ref[i - 1]]
        current = [negative] * width
        row = bytearray(width)
        for k in range(width):
            j = i - half_width + k
            if j < 0 or j > cols:
                continue
            if j == 0:
                current[k] = i * gap_score
                row[k] = up_step
                continue
            diagonal = previous[k] + scores[j - 1]
            up = previous[k + 1] + gap_score if k + 1 < width else negative
            left = current[k - 1] + gap_score if k > 0 else negative
            if diagonal >= up and diagonal >= left:
                current[k] = diagonal
                row[k] = diagonal_step
            elif up >= left:
                current[k] = up
                row[k] = up_step
            else:
                current[k] = left
                row[k] = left_step
        trace.append(row)
        previous = current
    score = previous[cols - rows + half_width]
    aligned_ref, aligned_var = [], []
    i, j = rows, cols
    while i > 0 or j > 0:
        direction = trace[i][j - i + half_width]
        if direction == diagonal_step:
            aligned_ref.append(ref[i - 1])
            aligned_var.append(var[j - 1])
            i -= 1
            j -= 1
        elif direction == up_step:
            aligned_ref.append(ref[i - 1])
            aligned_var.append("-")
            i -= 1
        elif direction == left_step:
            aligned_ref.append("-")
            aligned_var.append(var[j - 1])
            j -= 1
        else:
            raise MutationAnalysisError("Alignment traceback failed.")
    return "".join(reversed(aligned_ref)), "".join(reversed(aligned_var)), int(score)


def _banded_is_exact(score, rows, cols, half_width, match_score, gap_score):
    """True when no alignment leaving the band can reach the banded score.

    Leaving |j - i| <= w and ending at offset cols - rows takes at least
    G = 2w + 2 - |cols - rows| gap columns, and an alignment with G gap
    columns scores at most match * (rows + cols - G) / 2 + gap * G (every other
    column a match). That bound falls with G when gap < match / 2. A banded
    score strictly above it means every optimal path lies in the band; the
    traceback then makes the same tie-breaking decisions as the full matrix.
    """
    gaps = 2 * half_width + 2 - abs(cols - rows)
    bound = match_score * (rows + cols - gaps) / 2 + gap_score * gaps
    return score > bound


def global_align(
    reference: str,
    variant: str,
    *,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_score: int = -2,
) -> Tuple[str, str, int]:
    """Needleman-Wunsch global alignment with a linear gap cost.

    Ties are broken deterministically: diagonal > up (gap in the variant) >
    left (gap in the reference). The scores are kept for two rows only and the
    traceback as one byte per cell (audit P2-6), about 1 byte per cell instead
    of the ~45 of a matrix of Python objects; the recurrence and the
    tie-breaking are unchanged, so is every alignment.

    Above BANDED_ALIGNMENT_MIN_CELLS cells an exact banded computation is
    tried first (half-width 32 or more, doubled until _banded_is_exact
    proves that no alignment outside the band can compete); it returns the
    same alignment as the full matrix in O(n*w) time. Dissimilar sequences
    fall back to the full O(n*m) computation.
    """

    ref = normalize_iupac_sequence(
        reference
    )

    var = normalize_iupac_sequence(
        variant
    )

    rows = len(ref)
    cols = len(var)
    diagonal_step, up_step, left_step = 1, 2, 3

    # IUPAC-aware substitution score of each reference base against the
    # whole variant, computed once per distinct base.
    substitution: Dict[str, List[int]] = {}
    for base in set(ref):
        substitution[base] = [
            match_score if _nucleotide_match(base, other) else mismatch_score
            for other in var
        ]

    if ((rows + 1) * (cols + 1) > BANDED_ALIGNMENT_MIN_CELLS
            and 2 * gap_score < match_score and mismatch_score <= match_score):
        half_width = max(32, abs(cols - rows) + 16)
        while 2 * half_width < max(rows, cols):
            banded = _banded_align(ref, var, half_width, substitution, gap_score)
            if banded is not None and _banded_is_exact(
                    banded[2], rows, cols, half_width, match_score, gap_score):
                return banded
            half_width *= 2

    trace = [bytearray(cols + 1) for _ in range(rows + 1)]
    trace[0][1:] = bytes([left_step]) * cols
    previous = [j * gap_score for j in range(cols + 1)]

    for i in range(1, rows + 1):
        scores = substitution[ref[i - 1]]
        current = [0] * (cols + 1)
        current[0] = i * gap_score
        row = trace[i]
        row[0] = up_step

        for j in range(1, cols + 1):
            diagonal = previous[j - 1] + scores[j - 1]
            up = previous[j] + gap_score
            left = current[j - 1] + gap_score

            if diagonal >= up and diagonal >= left:
                current[j] = diagonal
                row[j] = diagonal_step
            elif up >= left:
                current[j] = up
                row[j] = up_step
            else:
                current[j] = left
                row[j] = left_step

        previous = current

    aligned_ref: List[str] = []
    aligned_var: List[str] = []

    i = rows
    j = cols

    while i > 0 or j > 0:

        direction = trace[i][j]

        if direction == diagonal_step:
            aligned_ref.append(ref[i - 1])
            aligned_var.append(var[j - 1])
            i -= 1
            j -= 1

        elif direction == up_step:
            aligned_ref.append(ref[i - 1])
            aligned_var.append("-")
            i -= 1

        elif direction == left_step:
            aligned_ref.append("-")
            aligned_var.append(var[j - 1])
            j -= 1

        else:
            raise MutationAnalysisError(
                "Alignment traceback failed."
            )

    return (
        "".join(reversed(aligned_ref)),
        "".join(reversed(aligned_var)),
        previous[-1],
    )


def alignment_identity(
    aligned_reference: str,
    aligned_variant: str,
) -> float:
    """
    Calculate Alignment Identity.

    Exact nucleotide matches among aligned positions containing two
    nucleotides. Gap positions are excluded from the denominator.
    """

    if len(aligned_reference) != len(
        aligned_variant
    ):
        raise MutationAnalysisError(
            "Aligned sequences must have equal length."
        )

    if not aligned_reference:
        return 100.0

    comparable = 0
    identical = 0

    for a, b in zip(
        aligned_reference,
        aligned_variant,
    ):

        if a == "-" or b == "-":
            continue

        comparable += 1

        if a == b:
            identical += 1

    if comparable == 0:
        return 0.0

    return (
        100.0
        * identical
        / comparable
    )


def _full_edit_distance(ref: str, var: str) -> int:
    previous = list(range(len(var) + 1))
    for i, ref_base in enumerate(ref, 1):
        current = [i]
        for j, var_base in enumerate(var, 1):
            current.append(min(
                current[-1] + 1,
                previous[j] + 1,
                previous[j - 1] + (0 if ref_base == var_base else 1),
            ))
        previous = current
    return previous[-1]


def _banded_edit_distance(ref: str, var: str, band: int) -> int:
    """Edit distance over paths with |i - j| <= band (larger outside it)."""
    rows, cols = len(ref), len(var)
    width = 2 * band + 1
    infinite = rows + cols + 1
    previous = [infinite] * width
    for k in range(width):                 # row 0: j = k - band
        j = k - band
        if 0 <= j <= cols:
            previous[k] = j
    for i in range(1, rows + 1):
        ref_base = ref[i - 1]
        current = [infinite] * width
        for k in range(width):
            j = i - band + k
            if j < 0 or j > cols:
                continue
            if j == 0:
                current[k] = i
                continue
            best = previous[k] + (0 if ref_base == var[j - 1] else 1)
            if k + 1 < width and previous[k + 1] + 1 < best:
                best = previous[k + 1] + 1
            if k > 0 and current[k - 1] + 1 < best:
                best = current[k - 1] + 1
            current[k] = best
        previous = current
    return previous[cols - rows + band]


def _edit_distance(ref: str, var: str) -> int:
    """Levenshtein distance; exact banded computation for large inputs (P2-6).

    A path of cost d never leaves |i - j| <= d, because every step off the
    diagonal is an insertion or deletion of cost 1. A banded result <= band
    is therefore the true distance; otherwise the band is doubled (Ukkonen's
    observation). Small inputs use the full matrix directly.
    """
    rows, cols = len(ref), len(var)
    if (rows + 1) * (cols + 1) > BANDED_ALIGNMENT_MIN_CELLS:
        band = max(abs(rows - cols), 32)
        while band < max(rows, cols):
            distance = _banded_edit_distance(ref, var, band)
            if distance <= band:
                return distance
            band *= 2
    return _full_edit_distance(ref, var)


def sequence_similarity(
    reference: str,
    variant: str,
) -> float:
    """
    Calculate Sequence Similarity using normalized edit distance.

    Unlike Alignment Identity, this metric does not ignore insertions or
    deletions. A substitution, insertion, or deletion contributes one edit.

    Formula:
        100 * (1 - edit_distance / max(sequence lengths))
    """

    ref = normalize_iupac_sequence(reference)
    var = normalize_iupac_sequence(variant)

    max_length = max(len(ref), len(var))

    if max_length == 0:
        return 100.0

    distance = _edit_distance(ref, var)

    return max(
        0.0,
        100.0 * (
            1.0
            - distance / max_length
        ),
    )


# ============================================================================
# ALIGNMENT BLOCKS
# ============================================================================

def build_alignment_blocks(
    aligned_reference: str,
    aligned_variant: str,
) -> List[AlignmentBlock]:

    if len(aligned_reference) != len(
        aligned_variant
    ):
        raise MutationAnalysisError(
            "Aligned sequences must have equal length."
        )

    blocks: List[AlignmentBlock] = []

    ref_pos = 1
    var_pos = 1

    current_operation: Optional[str] = None
    current_ref: List[str] = []
    current_var: List[str] = []

    current_ref_start = 1
    current_var_start = 1

    def flush() -> None:

        nonlocal current_operation
        nonlocal current_ref
        nonlocal current_var

        if current_operation is None:
            return

        ref_seq = "".join(
            current_ref
        )

        var_seq = "".join(
            current_var
        )

        ref_end = (
            current_ref_start
            + len(ref_seq.replace("-", ""))
            - 1
        )

        var_end = (
            current_var_start
            + len(var_seq.replace("-", ""))
            - 1
        )

        if current_operation == "insertion":
            ref_end = (
                current_ref_start - 1
            )

        if current_operation == "deletion":
            var_end = (
                current_var_start - 1
            )

        blocks.append(
            AlignmentBlock(
                operation=current_operation,
                reference_start=current_ref_start,
                reference_end=ref_end,
                variant_start=current_var_start,
                variant_end=var_end,
                reference_sequence=ref_seq,
                variant_sequence=var_seq,
            )
        )

        current_operation = None
        current_ref = []
        current_var = []

    for a, b in zip(
        aligned_reference,
        aligned_variant,
    ):

        if a != "-" and b != "-":

            # Exact equality is a true match. IUPAC compatibility is useful
            # for alignment scoring, but different symbols must remain
            # observable as sequence-level events.
            operation = (
                "match"
                if a == b
                else "substitution"
            )

        elif a == "-" and b != "-":
            operation = "insertion"

        elif a != "-" and b == "-":
            operation = "deletion"

        else:
            raise MutationAnalysisError(
                "Invalid alignment column '--'."
            )

        if operation != current_operation:

            flush()

            current_operation = operation
            current_ref_start = ref_pos
            current_var_start = var_pos

        current_ref.append(a)
        current_var.append(b)

        if a != "-":
            ref_pos += 1

        if b != "-":
            var_pos += 1

    flush()

    return blocks
