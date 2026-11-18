"""Mutation calls from an alignment, their classification and notation, and
left normalization of indels."""

from __future__ import annotations

from typing import List, Optional, Tuple

from analyses._mutation.sequence import normalize_iupac_sequence
from analyses._mutation.models import Mutation, MutationAnalysisError
from analyses._mutation.alignment import build_alignment_blocks


# ============================================================================
# MUTATION CLASSIFICATION
# ============================================================================

def is_transition(
    reference_base: str,
    alternative_base: str,
) -> bool:

    ref = reference_base.upper()
    alt = alternative_base.upper()

    if ref == alt:
        return False

    return (
        {ref, alt}
        in (
            {"A", "G"},
            {"C", "T"},
        )
    )


def classify_substitution(
    reference_base: str,
    alternative_base: str,
) -> str:

    ref = reference_base.upper()
    alt = alternative_base.upper()

    if ref == alt:
        return "Match"

    if (
        ref not in "ACGT"
        or alt not in "ACGT"
    ):
        return "Ambiguous"

    if is_transition(ref, alt):
        return "Transition"

    return "Transversion"


def _mutation_notation(
    mutation_type: str,
    position: int,
    end_position: int,
    reference: str,
    alternative: str,
) -> str:

    if mutation_type == "SNP":

        return (
            f"{position}"
            f"{reference}>{alternative}"
        )

    if mutation_type == "Insertion":

        return (
            f"{position}ins"
            f"{alternative}"
        )

    if mutation_type == "Deletion":

        if position == end_position:
            return (
                f"{position}del"
                f"{reference}"
            )

        return (
            f"{position}_{end_position}"
            f"del{reference}"
        )

    return (
        f"{position}_{end_position}"
        f"{reference}>{alternative}"
    )


def _make_mutation(
    mutation_id: int,
    mutation_type: str,
    position: int,
    end_position: int,
    reference: str,
    alternative: str,
) -> Mutation:

    ref = reference.replace(
        "-",
        "",
    )

    alt = alternative.replace(
        "-",
        "",
    )

    notation = _mutation_notation(
        mutation_type,
        position,
        end_position,
        ref,
        alt,
    )

    transition = None

    if mutation_type == "SNP":

        transition = is_transition(
            ref,
            alt,
        )

    return Mutation(
        mutation_id=mutation_id,
        mutation_type=mutation_type,
        position=position,
        end_position=end_position,
        reference=ref,
        alternative=alt,
        length_change=(
            len(alt) - len(ref)
        ),
        notation=notation,
        transition=transition,
    )


# ============================================================================
# MUTATION DETECTION
# ============================================================================

def detect_mutations_from_alignment(
    aligned_reference: str,
    aligned_variant: str,
) -> List[Mutation]:

    blocks = build_alignment_blocks(
        aligned_reference,
        aligned_variant,
    )

    mutations: List[Mutation] = []

    mutation_id = 1

    for block in blocks:

        if block.operation == "match":
            continue

        if block.operation == "substitution":

            ref_seq = block.reference_sequence
            var_seq = block.variant_sequence

            for offset, (
                ref_base,
                alt_base,
            ) in enumerate(
                zip(
                    ref_seq,
                    var_seq,
                )
            ):

                # Exact character equality means no mutation. IUPAC
                # compatibility must not suppress a real ambiguity event
                # such as G -> R or T -> Y.
                if ref_base == alt_base:
                    continue

                position = (
                    block.reference_start
                    + offset
                )

                mutation_type = (
                    "SNP"
                    if (
                        ref_base in "ACGT"
                        and alt_base in "ACGT"
                    )
                    else "Substitution"
                )

                mutations.append(
                    _make_mutation(
                        mutation_id,
                        mutation_type,
                        position,
                        position,
                        ref_base,
                        alt_base,
                    )
                )

                mutation_id += 1

        elif block.operation == "insertion":

            mutations.append(
                _make_mutation(
                    mutation_id,
                    "Insertion",
                    block.reference_start,
                    block.reference_start - 1,
                    "",
                    block.variant_sequence,
                )
            )

            mutation_id += 1

        elif block.operation == "deletion":

            mutations.append(
                _make_mutation(
                    mutation_id,
                    "Deletion",
                    block.reference_start,
                    block.reference_end,
                    block.reference_sequence,
                    "",
                )
            )

            mutation_id += 1

    return mutations


def _detect_simple_indel(
    reference: str,
    variant: str,
) -> Optional[List[Mutation]]:
    """
    Detect a single contiguous insertion/deletion using the maximal common
    prefix/suffix representation. This avoids ambiguous global-alignment
    placement inside short repeats for pure indel test cases.
    """
    if len(reference) == len(variant):
        return None

    prefix = 0
    max_prefix = min(len(reference), len(variant))
    while prefix < max_prefix and reference[prefix] == variant[prefix]:
        prefix += 1

    suffix = 0
    max_suffix = min(
        len(reference) - prefix,
        len(variant) - prefix,
    )
    while (
        suffix < max_suffix
        and reference[len(reference) - 1 - suffix]
        == variant[len(variant) - 1 - suffix]
    ):
        suffix += 1

    ref_end = len(reference) - suffix if suffix else len(reference)
    var_end = len(variant) - suffix if suffix else len(variant)

    ref_middle = reference[prefix:ref_end]
    var_middle = variant[prefix:var_end]

    if ref_middle and not var_middle:
        # Deletion starts at the first changed reference base.
        return [
            _make_mutation(
                1,
                "Deletion",
                prefix + 1,
                prefix + len(ref_middle),
                ref_middle,
                "",
            )
        ]

    if var_middle and not ref_middle:
        # Insertion is reported between the bases on either side.
        return [
            _make_mutation(
                1,
                "Insertion",
                prefix + 1,
                prefix,
                "",
                var_middle,
            )
        ]

    return None


# ============================================================================
# INDEL NORMALIZATION
# ============================================================================

def normalize_indel_left(
    reference: str,
    variant: str,
    position: int,
    reference_allele: str,
    alternative_allele: str,
) -> Tuple[
    int,
    str,
    str,
]:

    ref = normalize_iupac_sequence(
        reference
    )

    if position < 1:
        raise MutationAnalysisError(
            "Position must be >= 1."
        )

    if (
        reference_allele
        and alternative_allele
    ):

        return (
            position,
            reference_allele,
            alternative_allele,
        )

    # Conservative left-normalization for simple indels.
    pos = position
    ref_allele = reference_allele
    alt_allele = alternative_allele

    while pos > 1:

        if not ref_allele and alt_allele:

            previous = ref[pos - 2]

            if alt_allele[-1] == previous:
                alt_allele = (
                    previous
                    + alt_allele[:-1]
                )
                pos -= 1
            else:
                break

        elif ref_allele and not alt_allele:

            previous = ref[pos - 2]

            if ref_allele[-1] == previous:
                ref_allele = (
                    previous
                    + ref_allele[:-1]
                )
                pos -= 1
            else:
                break

        else:
            break

    return (
        pos,
        ref_allele,
        alt_allele,
    )
