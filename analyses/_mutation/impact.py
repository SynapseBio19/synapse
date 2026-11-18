"""Overlap of mutations with restriction sites, primers and guide positions."""

from __future__ import annotations

from typing import List, Mapping, Sequence, Tuple

from analyses._mutation.sequence import iupac_matches, normalize_iupac_sequence
from analyses._mutation.models import Mutation, MutationImpact, RestrictionSiteImpact


# ============================================================================
# RESTRICTION / PRIMER / CRISPR IMPACT
# ============================================================================

def _find_all_occurrences(
    sequence: str,
    motif: str,
) -> List[int]:

    seq = normalize_iupac_sequence(
        sequence
    )

    motif = normalize_iupac_sequence(
        motif
    )

    hits = []

    if len(motif) > len(seq):
        return hits

    for i in range(
        0,
        len(seq) - len(motif) + 1,
    ):

        if all(
            iupac_matches(
                seq[i + j],
                motif[j],
            )
            for j in range(len(motif))
        ):
            hits.append(i + 1)

    return hits


def restriction_site_impact(
    reference: str,
    variant: str,
    enzymes: Mapping[str, str],
) -> List[RestrictionSiteImpact]:

    ref = normalize_iupac_sequence(
        reference
    )

    var = normalize_iupac_sequence(
        variant
    )

    results = []

    for enzyme, motif in enzymes.items():

        ref_sites = _find_all_occurrences(
            ref,
            motif,
        )

        var_sites = _find_all_occurrences(
            var,
            motif,
        )

        lost = [
            position
            for position in ref_sites
            if position not in var_sites
        ]

        gained = [
            position
            for position in var_sites
            if position not in ref_sites
        ]

        unchanged = len(
            set(ref_sites)
            & set(var_sites)
        )

        results.append(
            RestrictionSiteImpact(
                enzyme=str(enzyme),
                recognition_sequence=str(
                    motif
                ),
                lost_sites=lost,
                gained_sites=gained,
                unchanged_sites=unchanged,
            )
        )

    return results


def _interval_overlap(
    mutation: Mutation,
    start: int,
    end: int,
) -> bool:

    mutation_start = mutation.position

    mutation_end = max(
        mutation.position,
        mutation.end_position,
    )

    return not (
        mutation_end < start
        or mutation_start > end
    )


def primer_mutation_overlap(
    mutations: Sequence[Mutation],
    primers: Sequence[
        Tuple[str, int, int]
    ],
) -> MutationImpact:

    impact = MutationImpact()

    for mutation in mutations:

        for primer_id, start, end in primers:

            if _interval_overlap(
                mutation,
                start,
                end,
            ):

                if primer_id not in (
                    impact.primer_overlap_ids
                ):

                    impact.primer_overlap_ids.append(
                        str(primer_id)
                    )

    impact.primer_overlap = bool(
        impact.primer_overlap_ids
    )

    return impact


def crispr_mutation_overlap(
    mutations: Sequence[Mutation],
    targets: Sequence[
        Tuple[str, int, int]
    ],
) -> MutationImpact:

    impact = MutationImpact()

    for mutation in mutations:

        for target_id, start, end in targets:

            if _interval_overlap(
                mutation,
                start,
                end,
            ):

                if target_id not in (
                    impact.crispr_overlap_ids
                ):

                    impact.crispr_overlap_ids.append(
                        str(target_id)
                    )

    impact.crispr_target_overlap = bool(
        impact.crispr_overlap_ids
    )

    return impact
