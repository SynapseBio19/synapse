"""Comparing a reference with one or several variants: alignment, mutation
calls, coding consequences and statistics in one result."""

from __future__ import annotations

import math
from typing import Dict, List, Mapping, Sequence, Tuple

from analyses._mutation.sequence import normalize_iupac_sequence, reverse_complement_iupac
from analyses._mutation.models import (
    Mutation,
    MutationAnalysisError,
    MutationAnalysisResult,
    MutationStatistics,
)
from analyses._mutation.alignment import alignment_identity, global_align, sequence_similarity
from analyses._mutation.detection import detect_mutations_from_alignment, _detect_simple_indel
from analyses._mutation.coding import annotate_mutations_with_coding_consequences


def detect_mutations(
    reference: str,
    variant: str,
    *,
    alignment_match_score: int = 2,
    alignment_mismatch_score: int = -1,
    alignment_gap_score: int = -2,
) -> MutationAnalysisResult:
    ref = normalize_iupac_sequence(reference)
    var = normalize_iupac_sequence(variant)

    # A reverse-complement sequence represents the same nucleotide sequence
    # in the opposite orientation. It should not be reported as dozens of
    # substitutions/indels merely because the coordinate orientation changed.
    reverse_complement_match = (
        var == reverse_complement_iupac(ref)
        and var != ref
    )

    aligned_ref, aligned_var, score = global_align(
        ref,
        var,
        match_score=alignment_match_score,
        mismatch_score=alignment_mismatch_score,
        gap_score=alignment_gap_score,
    )

    if reverse_complement_match:
        return MutationAnalysisResult(
            reference_name="Reference",
            variant_name="Variant",
            reference_length=len(ref),
            variant_length=len(var),
            mutations=[],
            statistics=calculate_mutation_statistics(
                [],
                reference_length=len(ref),
            ),
            aligned_reference=aligned_ref,
            aligned_variant=aligned_var,
            alignment_score=score,
            identity_percent=100.0,
            sequence_similarity_percent=100.0,
            notes=[
                "Variant is the reverse complement of the reference sequence.",
                "No sequence-level mutations were reported because the difference is orientation-only.",
            ],
        )

    # For a pure length change, use the maximal common prefix/suffix
    # representation. Needleman-Wunsch can otherwise split one biological
    # indel into multiple events when short repeats are present.
    simple_indel = _detect_simple_indel(ref, var)

    if simple_indel is not None:
        mutations = simple_indel
    else:
        mutations = detect_mutations_from_alignment(
            aligned_ref,
            aligned_var,
        )

    annotate_mutations_with_coding_consequences(
        mutations,
        ref,
        var,
    )

    statistics = calculate_mutation_statistics(
        mutations,
        reference_length=len(ref),
    )

    return MutationAnalysisResult(
        reference_name="Reference",
        variant_name="Variant",
        reference_length=len(ref),
        variant_length=len(var),
        mutations=mutations,
        statistics=statistics,
        aligned_reference=aligned_ref,
        aligned_variant=aligned_var,
        alignment_score=score,
        identity_percent=alignment_identity(
            aligned_ref,
            aligned_var,
        ),
        sequence_similarity_percent=sequence_similarity(
            ref,
            var,
        ),
    )


def calculate_mutation_statistics(
    mutations: Sequence[Mutation],
    *,
    reference_length: int,
) -> MutationStatistics:

    if reference_length < 1:
        raise MutationAnalysisError(
            "Reference length must be positive."
        )

    stats = MutationStatistics()

    stats.total_mutations = len(
        mutations
    )

    for mutation in mutations:

        if mutation.mutation_type in {
            "SNP",
            "Substitution",
        }:

            stats.substitutions += 1

        if mutation.mutation_type == "SNP":

            stats.snps += 1

            if mutation.transition is True:
                stats.transitions += 1

            elif mutation.transition is False:
                stats.transversions += 1

        elif mutation.mutation_type == "Insertion":

            stats.insertions += 1

            stats.inserted_bases += len(
                mutation.alternative
            )

            if (
                len(mutation.alternative)
                % 3
                == 0
            ):
                stats.in_frame_indels += 1

            else:
                stats.frameshifts += 1

        elif mutation.mutation_type == "Deletion":

            stats.deletions += 1

            stats.deleted_bases += len(
                mutation.reference
            )

            if (
                len(mutation.reference)
                % 3
                == 0
            ):
                stats.in_frame_indels += 1

            else:
                stats.frameshifts += 1

        elif mutation.mutation_type == "Complex":

            stats.complex_changes += 1

        effect = mutation.effect

        if effect == "Synonymous":
            stats.synonymous += 1

        elif effect == "Missense":
            stats.missense += 1

        elif effect == "Nonsense":
            stats.nonsense += 1

        elif effect == "Start-loss":
            stats.start_loss += 1

        elif effect == "Stop-loss":
            stats.stop_loss += 1

        elif effect in {
            None,
            "Unknown",
        }:
            stats.unknown_effects += 1

    stats.net_length_change = (
        stats.inserted_bases
        - stats.deleted_bases
    )

    if stats.transversions == 0:

        if stats.transitions > 0:
            stats.transition_transversion_ratio = (
                math.inf
            )
        else:
            stats.transition_transversion_ratio = (
                None
            )

    else:

        stats.transition_transversion_ratio = (
            stats.transitions
            / stats.transversions
        )

    stats.affected_codons = len({
        mutation.codon_index
        for mutation in mutations
        if mutation.codon_index is not None
    })

    stats.mutation_density_per_kb = (
        stats.total_mutations
        / reference_length
        * 1000.0
    )

    return stats


# ============================================================================
# COMPARISON
# ============================================================================

def compare_sequences(
    reference: str,
    variant: str,
    *,
    reference_name: str = "Reference",
    variant_name: str = "Variant",
    alignment_match_score: int = 2,
    alignment_mismatch_score: int = -1,
    alignment_gap_score: int = -2,
) -> MutationAnalysisResult:

    result = detect_mutations(
        reference,
        variant,
        alignment_match_score=alignment_match_score,
        alignment_mismatch_score=alignment_mismatch_score,
        alignment_gap_score=alignment_gap_score,
    )

    result.reference_name = (
        reference_name
    )

    result.variant_name = (
        variant_name
    )

    longest = max(result.reference_length, result.variant_length)
    if longest > MAX_RECOMMENDED_ALIGNMENT_LENGTH:
        result.notes.append(
            f"Sequences up to {longest} nt were aligned globally "
            f"(recommended maximum {MAX_RECOMMENDED_ALIGNMENT_LENGTH} nt): the "
            "global alignment visits (n+1)x(m+1) cells and keeps one traceback "
            "byte per cell (about 1 byte per cell: 4 MB at 2,000 nt, 100 MB at "
            "10,000 nt). Similar sequences are aligned in an exact band (about "
            "0.2 s at 10,000 nt); unrelated ones take time growing with n x m "
            "(about 23 s at 10,000 nt). "
            "For genome-scale comparisons use VCF-based variant analysis."
        )

    return result


def compare_multiple_variants(
    reference: str,
    variants: Sequence[str],
    *,
    reference_name: str = "Reference",
) -> Dict[str, MutationAnalysisResult]:

    results: Dict[
        str,
        MutationAnalysisResult,
    ] = {}

    for index, variant in enumerate(
        variants,
        1,
    ):

        name = (
            f"Variant {index}"
        )

        results[name] = compare_sequences(
            reference,
            variant,
            reference_name=reference_name,
            variant_name=name,
        )

    return results


def _mutation_key(
    mutation: Mutation,
) -> Tuple:

    return (
        mutation.mutation_type,
        mutation.position,
        mutation.end_position,
        mutation.reference,
        mutation.alternative,
    )


def shared_mutations(
    results: Mapping[
        str,
        MutationAnalysisResult,
    ],
) -> List[Mutation]:

    if not results:
        return []

    mutation_sets = [
        {
            _mutation_key(m)
            for m in result.mutations
        }
        for result in results.values()
    ]

    shared_keys = set.intersection(
        *mutation_sets
    )

    first_result = next(
        iter(results.values())
    )

    return [
        mutation
        for mutation in first_result.mutations
        if _mutation_key(mutation)
        in shared_keys
    ]


def unique_mutations(
    result: MutationAnalysisResult,
    other_results: Sequence[
        MutationAnalysisResult
    ],
) -> List[Mutation]:

    other_keys = {
        _mutation_key(mutation)
        for other in other_results
        for mutation in other.mutations
    }

    return [
        mutation
        for mutation in result.mutations
        if _mutation_key(mutation)
        not in other_keys
    ]


# ============================================================================
# CODING-FRAME ASSUMPTIONS (audit P1-10)
# ============================================================================

# Global alignment needs O(n*m) time and memory. Above this length a warning
# is attached; the analysis still runs. The figure is a practical limit for an
# interactive pure-Python tool, not a scientific threshold.
MAX_RECOMMENDED_ALIGNMENT_LENGTH = 2000


# ============================================================================
# PUBLIC API
# ============================================================================

analyze_mutations = detect_mutations
