"""``analyze_variants`` and ``analyze_vcf_text``."""

from __future__ import annotations

from typing import Optional, Sequence

from analyses._variant.models import Variant, VariantAnalysisResult, VariantFilter
from analyses._variant.vcf import parse_vcf, validate_variant
from analyses._variant.normalization import normalize_variants
from analyses._variant.annotation import annotate_from_info
from analyses._variant.statistics import (
    calculate_cohort_summary,
    calculate_variant_statistics,
    filter_variants,
)


# ============================================================================
# High-level analysis
# ============================================================================

def analyze_variants(
    variants: Sequence[Variant],
    *,
    name: str = "Variant Analysis",
    filter_criteria: Optional[VariantFilter] = None,
    annotate_info: bool = True,
    normalize: bool = False,
    reference_sequence: Optional[str] = None,
) -> VariantAnalysisResult:
    """Run the main variant-analysis workflow."""
    working = list(variants)

    for variant in working:
        issues = validate_variant(variant)
        if issues:
            variant.annotation.setdefault("validation_warnings", []).extend(issues)

    if annotate_info:
        working = annotate_from_info(working)

    if normalize:
        working = normalize_variants(
            working,
            reference_sequence=reference_sequence,
        )

    result = VariantAnalysisResult(
        name=name,
        variants=working,
        statistics=calculate_variant_statistics(working),
    )

    if filter_criteria is not None:
        result.filter_result = filter_variants(
            working,
            filter_criteria,
        )

    result.cohort_summary = calculate_cohort_summary(working)
    return result


def analyze_vcf_text(
    text: str,
    *,
    name: str = "VCF Analysis",
    filter_criteria: Optional[VariantFilter] = None,
    normalize: bool = False,
    reference_sequence: Optional[str] = None,
) -> VariantAnalysisResult:
    variants, _, _ = parse_vcf(text)
    return analyze_variants(
        variants,
        name=name,
        filter_criteria=filter_criteria,
        normalize=normalize,
        reference_sequence=reference_sequence,
    )
