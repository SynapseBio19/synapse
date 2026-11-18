"""Validation of a result, the text and JSON reports, the advanced analysis
and the mutation map."""

from __future__ import annotations

import json
import math
from core.textformat import field_line, major_rule, minor_rule
from dataclasses import asdict
from typing import Any, Dict, Optional, Sequence

from analyses._mutation.models import (
    AdvancedMutationAnalysis,
    Mutation,
    MutationAnalysisError,
    MutationAnalysisResult,
)
from analyses._mutation.hgvs import hgvs_c_allele, mutation_hgvs_c, mutation_hgvs_p
from analyses._mutation.context import (
    build_mutation_spectrum,
    extract_mutation_context,
    find_mutation_hotspots,
    sequence_quality,
)


# ============================================================================
# RESULT VALIDATION
# ============================================================================

def validate_mutation_result(
    result: MutationAnalysisResult,
) -> bool:

    if result.reference_length < 1:
        raise MutationAnalysisError(
            "Reference length is invalid."
        )

    if result.variant_length < 1:
        raise MutationAnalysisError(
            "Variant length is invalid."
        )

    if len(result.mutations) != (
        result.statistics.total_mutations
    ):
        raise MutationAnalysisError(
            "Mutation count and statistics disagree."
        )

    if not (
        0.0 <= result.identity_percent <= 100.0
        and 0.0 <= result.sequence_similarity_percent <= 100.0
    ):
        raise MutationAnalysisError(
            "Identity/similarity percentages are out of range."
        )

    expected_net = (
        result.variant_length
        - result.reference_length
    )

    if (
        result.statistics.net_length_change
        != expected_net
    ):
        raise MutationAnalysisError(
            "Net length change is inconsistent with sequence lengths."
        )

    ids = [
        mutation.mutation_id
        for mutation in result.mutations
    ]

    if ids != list(
        range(
            1,
            len(ids) + 1,
        )
    ):
        raise MutationAnalysisError(
            "Mutation IDs are not sequential."
        )

    return True


# ============================================================================
# FORMATTING / REPORTING
# ============================================================================

def format_mutation(
    mutation: Mutation,
    reference_sequence: Optional[str] = None,
) -> str:
    """One report line; the DNA change in HGVS syntax (see mutation_hgvs_c)."""

    extra = ""

    if mutation.effect:
        protein = mutation_hgvs_p(mutation)
        extra = (
            f"  [{mutation.effect}"
            f"{' ' + protein if protein else ''}]"
        )

    return (
        f"#{mutation.mutation_id:<3} "
        f"{mutation.mutation_type:<12} "
        f"position {mutation.position:<6} "
        f"{mutation_hgvs_c(mutation, reference_sequence)}"
        f"{extra}"
    )


def build_mutation_report(
    result: MutationAnalysisResult,
    *,
    include_alignment: bool = False,
) -> str:

    lines = [
        major_rule(),
        "MUTATION ANALYSIS REPORT",
        major_rule(),
        field_line('Reference', f"{result.reference_name}"),
        field_line('Variant', f"{result.variant_name}"),
        field_line('Reference length', f"{result.reference_length} bp"),
        field_line('Variant length', f"{result.variant_length} bp"),
        field_line('Alignment identity', f"{result.identity_percent:.2f}%"),
        field_line('Sequence similarity', f"{result.sequence_similarity_percent:.2f}%"),
        field_line('Mutations', f"{len(result.mutations)}"),
        "",
        "MUTATION SUMMARY",
        minor_rule(),
        "Alignment identity = exact matches among aligned non-gap bases.",
        "Sequence similarity = normalized edit-distance similarity.",
    ]

    s = result.statistics

    ratio = (
        "inf"
        if s.transition_transversion_ratio
        == math.inf
        else (
            f"{s.transition_transversion_ratio:.3f}"
            if s.transition_transversion_ratio
            is not None
            else "N/A"
        )
    )

    lines.extend(
        [
            field_line('SNPs', f"{s.snps}"),
            field_line('Substitutions', f"{s.substitutions}"),
            field_line('Insertions', f"{s.insertions}"),
            field_line('Deletions', f"{s.deletions}"),
            field_line('Complex changes', f"{s.complex_changes}"),
            field_line('Transitions', f"{s.transitions}"),
            field_line('Transversions', f"{s.transversions}"),
            field_line('Ti/Tv ratio', f"{ratio}"),
            field_line('Inserted bases', f"{s.inserted_bases}"),
            field_line('Deleted bases', f"{s.deleted_bases}"),
            field_line('Net length', f"{s.net_length_change:+d} bp"),
            field_line('Mutation density', f"{s.mutation_density_per_kb:.2f}/kb"),
            "",
            "CODING CONSEQUENCES",
            minor_rule(),
            field_line('Synonymous', f"{s.synonymous}"),
            field_line('Missense', f"{s.missense}"),
            field_line('Nonsense', f"{s.nonsense}"),
            field_line('Start-loss', f"{s.start_loss}"),
            field_line('Stop-loss', f"{s.stop_loss}"),
            field_line('Frameshifts', f"{s.frameshifts}"),
            field_line('In-frame indels', f"{s.in_frame_indels}"),
            "",
            "MUTATIONS",
            minor_rule(),
        ]
    )

    if result.mutations:

        reference_sequence = result.aligned_reference.replace("-", "") or None

        if reference_sequence:
            lines.append(field_line(
                'HGVS allele',
                f"{hgvs_c_allele(result.mutations, reference_sequence)}"))

        for mutation in result.mutations:

            lines.append(
                format_mutation(
                    mutation,
                    reference_sequence,
                )
            )

    else:

        lines.append(
            "No mutations detected."
        )

    if result.coding_analysis:

        c = result.coding_analysis

        lines.extend(
            [
                "",
                "CODING ANALYSIS",
                minor_rule(),
                field_line('Reference protein', f"{c.reference_protein}"),
                field_line('Variant protein', f"{c.alternative_protein}"),
                field_line('Affected codons', f"{len(c.affected_codons)}"),
                field_line('Frameshift', f"{'YES' if c.frameshift else 'NO'}"),
            ]
        )

        if c.frameshift_position is not None:

            lines.append(
                field_line('Frameshift position', f"{c.frameshift_position}")
            )

        if c.premature_stop_position is not None:

            lines.append(
                field_line('Premature stop', f"codon "
                f"{c.premature_stop_position}")
            )

    if include_alignment:

        lines.extend(
            [
                "",
                "GLOBAL ALIGNMENT",
                minor_rule(),
                result.aligned_reference,
                result.aligned_variant,
            ]
        )

    return "\n".join(lines)


# The terminal dialogs of this module are in ui/mutation_menu.py (P2-2).


# ============================================================================
# MUTATION SIGNATURE
# ============================================================================

def mutation_signature(
    mutations: Sequence[Mutation],
) -> Dict[str, Any]:

    spectrum = build_mutation_spectrum(
        mutations
    )

    return {
        "total": len(mutations),
        "snps": sum(
            1
            for mutation in mutations
            if mutation.mutation_type == "SNP"
        ),
        "insertions": sum(
            1
            for mutation in mutations
            if mutation.mutation_type
            == "Insertion"
        ),
        "deletions": sum(
            1
            for mutation in mutations
            if mutation.mutation_type
            == "Deletion"
        ),
        "transitions": spectrum.transitions,
        "transversions": spectrum.transversions,
        "spectrum_12": spectrum.substitutions_12,
        "spectrum_6": spectrum.substitutions_6,
    }


# ============================================================================
# ADVANCED ANALYSIS
# ============================================================================

def create_advanced_analysis(
    result: MutationAnalysisResult,
) -> AdvancedMutationAnalysis:

    advanced = AdvancedMutationAnalysis()

    advanced.reference_quality = sequence_quality(
        result.aligned_reference.replace(
            "-",
            "",
        )
    )

    advanced.variant_quality = sequence_quality(
        result.aligned_variant.replace(
            "-",
            "",
        )
    )

    advanced.spectra = build_mutation_spectrum(
        result.mutations
    )

    for mutation in result.mutations:

        advanced.contexts.append(
            extract_mutation_context(
                mutation,
                result.aligned_reference.replace(
                    "-",
                    "",
                ),
                result.aligned_variant.replace(
                    "-",
                    "",
                ),
            )
        )

    advanced.hotspots = find_mutation_hotspots(
        result.mutations,
        result.reference_length,
    )

    if (
        advanced.reference_quality.ambiguous_count
        or advanced.variant_quality.ambiguous_count
    ):

        advanced.ambiguity_warnings.append(
            "Ambiguous/IUPAC bases may limit exact mutation interpretation."
        )

    if (
        len(result.aligned_reference)
        != len(result.aligned_variant)
    ):

        advanced.normalization_notes.append(
            "Alignment lengths were inconsistent."
        )

    return advanced


def advanced_mutation_report(
    result: MutationAnalysisResult,
    advanced: Optional[
        AdvancedMutationAnalysis
    ] = None,
) -> str:

    if advanced is None:
        advanced = create_advanced_analysis(
            result
        )

    lines = [
        build_mutation_report(
            result
        ),
        "",
        "ADVANCED SEQUENCE ANALYSIS",
        minor_rule(),
    ]

    if (
        advanced.reference_quality
        and advanced.variant_quality
    ):

        lines.extend(
            [
                field_line('Reference GC', f"{advanced.reference_quality.gc_percent:.2f}%"),
                field_line('Variant GC', f"{advanced.variant_quality.gc_percent:.2f}%"),
                field_line('Reference complexity', f"{advanced.reference_quality.complexity_percent:.2f}%"),
                field_line('Variant complexity', f"{advanced.variant_quality.complexity_percent:.2f}%"),
            ]
        )

    lines.extend(
        [
            "",
            "MUTATION SPECTRUM",
            minor_rule(),
        ]
    )

    if advanced.spectra:

        for key, count in (
            advanced.spectra
            .substitutions_12
            .items()
        ):

            lines.append(
                f"{key:<8} {count}"
            )

    lines.extend(
        [
            "",
            "HOTSPOTS",
            minor_rule(),
        ]
    )

    if advanced.hotspots:

        for hotspot in advanced.hotspots:

            lines.append(
                f"{hotspot.start}-{hotspot.end}: "
                f"{hotspot.mutation_count} mutations "
                f"({hotspot.density_per_kb:.2f}/kb)"
            )

    else:

        lines.append(
            "No hotspots detected under current thresholds."
        )

    if advanced.ambiguity_warnings:

        lines.extend(
            [
                "",
                "AMBIGUITY / QUALITY WARNINGS",
                minor_rule(),
            ]
        )

        lines.extend(
            f"- {warning}"
            for warning in advanced.ambiguity_warnings
        )

    return "\n".join(lines)


def mutation_map(
    mutations: Sequence[Mutation],
    sequence_length: int,
    *,
    width: int = 70,
) -> str:

    if sequence_length < 1:
        raise MutationAnalysisError(
            "sequence_length must be positive."
        )

    width = max(
        20,
        width,
    )

    line = [" "] * width

    for mutation in mutations:

        index = min(
            width - 1,
            max(
                0,
                int(
                    (
                        mutation.position - 1
                    )
                    / sequence_length
                    * width
                ),
            ),
        )

        line[index] = "^"

    return "".join(line)


def result_to_json(
    result: MutationAnalysisResult,
    *,
    include_advanced: bool = True,
) -> str:

    payload = result.to_dict()

    if include_advanced:

        payload["advanced"] = asdict(
            create_advanced_analysis(
                result
            )
        )

    return json.dumps(
        payload,
        indent=2,
    )
mutation_report = build_mutation_report
