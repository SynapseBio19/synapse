"""Filters, statistics, cohort summary, set comparison, grouping and intervals."""

from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

from analyses._variant.models import (
    CohortVariantSummary,
    SampleVariantSummary,
    Variant,
    VariantAnalysisError,
    VariantComparison,
    VariantFilter,
    VariantFilterResult,
    VariantStatistics,
)
from analyses._variant.vcf import genotype_category, parse_genotype
from analyses._variant.annotation import synchronize_variant_fields


# ============================================================================
# Filtering
# ============================================================================

def variant_passes_filter(
    variant: Variant,
    criteria: VariantFilter,
) -> Tuple[bool, List[str]]:
    """Evaluate one variant against all configured criteria."""
    reasons: List[str] = []

    if criteria.minimum_quality is not None:
        if variant.quality is None:
            if not criteria.include_unknown_quality:
                reasons.append("Missing quality.")
        elif variant.quality < criteria.minimum_quality:
            reasons.append("Quality below minimum.")

    if criteria.maximum_quality is not None:
        if variant.quality is not None and variant.quality > criteria.maximum_quality:
            reasons.append("Quality above maximum.")

    if criteria.require_pass and variant.filter_status != "PASS":
        reasons.append("FILTER is not PASS.")

    if criteria.allowed_filter_statuses is not None:
        allowed = {s.upper() for s in criteria.allowed_filter_statuses}
        if variant.filter_status.upper() not in allowed:
            reasons.append("FILTER status not allowed.")

    if criteria.variant_types is not None:
        allowed = {s.upper() for s in criteria.variant_types}
        if variant.variant_type.upper() not in allowed:
            reasons.append("Variant type not allowed.")

    if criteria.chromosomes is not None:
        if variant.chromosome not in criteria.chromosomes:
            reasons.append("Chromosome not allowed.")

    if criteria.genes is not None:
        if variant.gene is None:
            reasons.append("Gene annotation missing.")
        elif variant.gene not in criteria.genes:
            reasons.append("Gene not allowed.")

    if criteria.consequences is not None:
        if variant.consequence is None:
            reasons.append("Consequence annotation missing.")
        elif variant.consequence not in criteria.consequences:
            reasons.append("Consequence not allowed.")

    if criteria.minimum_allele_frequency is not None:
        if variant.allele_frequency is None:
            if not criteria.include_unknown_frequency:
                reasons.append("Missing allele frequency.")
        elif variant.allele_frequency < criteria.minimum_allele_frequency:
            reasons.append("Allele frequency below minimum.")

    if criteria.maximum_allele_frequency is not None:
        if variant.allele_frequency is None:
            if not criteria.include_unknown_frequency:
                reasons.append("Missing allele frequency.")
        elif variant.allele_frequency > criteria.maximum_allele_frequency:
            reasons.append("Allele frequency above maximum.")

    if criteria.minimum_allele_count is not None:
        if variant.allele_count is None:
            reasons.append("Missing allele count.")
        elif variant.allele_count < criteria.minimum_allele_count:
            reasons.append("Allele count below minimum.")

    if criteria.maximum_allele_count is not None:
        if variant.allele_count is not None and variant.allele_count > criteria.maximum_allele_count:
            reasons.append("Allele count above maximum.")

    return not reasons, reasons


def filter_variants(
    variants: Sequence[Variant],
    criteria: VariantFilter,
) -> VariantFilterResult:
    variants = synchronize_variant_fields(variants)
    passed, failed, reasons = [], [], {}

    for variant in variants:
        ok, why = variant_passes_filter(variant, criteria)
        if ok:
            passed.append(variant)
        else:
            failed.append(variant)
            reasons[variant.key] = why

    return VariantFilterResult(passed=passed, failed=failed, reasons=reasons)


# ============================================================================
# Statistics
# ============================================================================

def is_transition_variant(variant: Variant) -> Optional[bool]:
    if not variant.is_snv:
        return None
    pair = {variant.reference, variant.alternative}
    return pair in ({"A", "G"}, {"C", "T"})


def calculate_variant_statistics(
    variants: Sequence[Variant],
    *,
    rare_frequency_threshold: float = 0.01,
) -> VariantStatistics:
    variants = synchronize_variant_fields(variants)
    stats = VariantStatistics()
    stats.total_variants = len(variants)
    stats.unique_variants = len({v.key for v in variants})

    qualities = []
    frequencies = []
    genes = set()
    chromosomes = set()
    consequences = Counter()

    for variant in variants:
        variant_type = variant.variant_type
        if variant_type == "SNV":
            stats.snvs += 1
            ti_tv = is_transition_variant(variant)
            if ti_tv is True:
                stats.transition_count += 1
            elif ti_tv is False:
                stats.transversion_count += 1
        elif variant_type == "MNV":
            stats.mnvs += 1
        elif variant_type == "INSERTION":
            stats.insertions += 1
            stats.indels += 1
        elif variant_type == "DELETION":
            stats.deletions += 1
            stats.indels += 1
        elif variant_type == "INDEL":
            stats.indels += 1
        elif variant_type == "COMPLEX":
            stats.complex += 1
        elif variant_type == "SYMBOLIC":
            stats.symbolic += 1
        else:
            stats.unknown += 1

        # VCF FILTER: "PASS" = passed all filters; "." = filters not
        # applied (not a failure); anything else names failed filters.
        if variant.filter_status == "PASS":
            stats.passed_filters += 1
        elif variant.filter_status in {"UNFILTERED", "UNKNOWN", "."}:
            stats.unfiltered += 1
        else:
            stats.failed_filters += 1

        if variant.quality is None:
            stats.missing_quality += 1
        else:
            qualities.append(variant.quality)

        if variant.allele_frequency is not None:
            frequencies.append(variant.allele_frequency)
            if variant.allele_frequency <= rare_frequency_threshold:
                stats.rare_variants += 1

        if variant.gene:
            genes.add(variant.gene)
        chromosomes.add(variant.chromosome)

        if variant.consequence:
            consequences[variant.consequence] += 1

    stats.genes_affected = len(genes)
    stats.chromosomes_affected = len(chromosomes)
    stats.consequence_counts = dict(consequences)

    if qualities:
        stats.mean_quality = statistics.mean(qualities)
        stats.median_quality = statistics.median(qualities)
        stats.min_quality = min(qualities)
        stats.max_quality = max(qualities)

    if frequencies:
        stats.mean_allele_frequency = statistics.mean(frequencies)
        stats.median_allele_frequency = statistics.median(frequencies)

    if stats.transversion_count:
        stats.transition_transversion_ratio = (
            stats.transition_count / stats.transversion_count
        )
    elif stats.transition_count:
        stats.transition_transversion_ratio = math.inf

    return stats


# ============================================================================
# Cohort / sample analysis
# ============================================================================

def sample_variant_summary(
    variants: Sequence[Variant],
    sample: str,
) -> SampleVariantSummary:
    summary = SampleVariantSummary(sample=sample)

    for variant in variants:
        summary.total_variant_records += 1
        gt = variant.sample_data.get(sample, {}).get("GT")
        category = genotype_category(gt)

        if category == "missing":
            summary.missing_calls += 1
        elif category == "homozygous_ref":
            summary.homozygous_ref_calls += 1
        elif category == "heterozygous":
            summary.heterozygous_calls += 1
            summary.called_variants += 1
        elif category == "homozygous_alt":
            summary.homozygous_alt_calls += 1
            summary.called_variants += 1
        elif category == "partial":
            summary.missing_calls += 1

        alleles, _ = parse_genotype(gt)
        if alleles:
            summary.allele_count += sum(
                allele is not None and allele > 0 for allele in alleles
            )

    return summary


def calculate_cohort_summary(
    variants: Sequence[Variant],
) -> CohortVariantSummary:
    samples = sorted({
        sample
        for variant in variants
        for sample in variant.sample_data
    })

    prevalence = Counter()

    for variant in variants:
        for sample, fields in variant.sample_data.items():
            gt = fields.get("GT")
            alleles, _ = parse_genotype(gt)
            if alleles and any(
                allele is not None and allele > 0
                for allele in alleles
            ):
                prevalence[variant.key] += 1

    unique_keys = {v.key for v in variants}
    shared = sum(count >= 2 for count in prevalence.values())
    private = sum(count == 1 for count in prevalence.values())

    return CohortVariantSummary(
        sample_count=len(samples),
        variant_count=len(variants),
        unique_variant_count=len(unique_keys),
        shared_variant_count=shared,
        private_variant_count=private,
        variant_prevalence=dict(prevalence),
        sample_summaries=[
            sample_variant_summary(variants, sample)
            for sample in samples
        ],
    )


# ============================================================================
# Comparison / overlap
# ============================================================================

def variant_key_set(variants: Sequence[Variant]) -> Set[Tuple[str, int, str, str]]:
    return {variant.key for variant in variants}


def compare_variant_sets(
    left: Sequence[Variant],
    right: Sequence[Variant],
) -> VariantComparison:
    left_map = {v.key: v for v in left}
    right_map = {v.key: v for v in right}

    left_keys = set(left_map)
    right_keys = set(right_map)
    shared_keys = left_keys & right_keys
    union = left_keys | right_keys

    return VariantComparison(
        only_in_left=[left_map[k] for k in sorted(left_keys - right_keys)],
        only_in_right=[right_map[k] for k in sorted(right_keys - left_keys)],
        shared=[left_map[k] for k in sorted(shared_keys)],
        left_count=len(left_keys),
        right_count=len(right_keys),
        shared_count=len(shared_keys),
        jaccard_similarity=(
            len(shared_keys) / len(union) if union else 1.0
        ),
    )


def variants_by_chromosome(
    variants: Sequence[Variant],
) -> Dict[str, List[Variant]]:
    grouped = defaultdict(list)
    for variant in variants:
        grouped[variant.chromosome].append(variant)
    return dict(grouped)


def variants_by_gene(
    variants: Sequence[Variant],
) -> Dict[str, List[Variant]]:
    grouped = defaultdict(list)
    for variant in variants:
        grouped[variant.gene or "Unknown"].append(variant)
    return dict(grouped)


def variants_by_consequence(
    variants: Sequence[Variant],
) -> Dict[str, List[Variant]]:
    grouped = defaultdict(list)
    for variant in variants:
        grouped[variant.consequence or "unknown"].append(variant)
    return dict(grouped)


def variants_in_interval(
    variants: Sequence[Variant],
    chromosome: str,
    start: int,
    end: int,
) -> List[Variant]:
    if start < 1 or end < start:
        raise VariantAnalysisError("Invalid interval.")
    return [
        v for v in variants
        if v.chromosome == chromosome
        and v.position <= end
        and (v.end_position or v.position) >= start
    ]


def interval_variant_density(
    variants: Sequence[Variant],
    chromosome: str,
    start: int,
    end: int,
) -> float:
    interval_length = end - start + 1
    if interval_length <= 0:
        raise VariantAnalysisError("Invalid interval.")
    count = len(variants_in_interval(variants, chromosome, start, end))
    return count / interval_length * 1000.0
