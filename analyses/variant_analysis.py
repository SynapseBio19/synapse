"""
Synapse - Variant Analysis
==========================

Professional variant-level analysis layer for Project Synapse.

Scope
-----
This module is deliberately separated from mutation.py.

mutation.py:
    sequence-level mutation detection and molecular consequence analysis.

variant_analysis.py:
    variant records, VCF-like input, filtering, normalization, annotation
    containers, cohort analysis, variant comparison, frequencies, quality
    metrics, summaries, reporting, JSON/TSV export, and CLI.

The module is self-contained and does not require third-party packages.
When mutation.py is available, sequence-derived consequence analysis can be
connected through the optional adapter functions near the end of the file.

Important scientific boundary
-----------------------------
This module does not pretend to provide clinical interpretation or database
annotation without evidence. Fields such as clinical significance,
population frequency, gene, transcript, and pathogenicity are represented
explicitly and are only populated when supplied by the user or an external
annotation source.
"""

# Phase 9: the implementation moved to the internal package
# ``analyses._variant`` (one module per responsibility, see its __init__).
# This module remains the public interface: it re-exports every name it
# defined before and keeps ``__all__``.

from __future__ import annotations

from core.textformat import minor_rule
from analyses._variant.models import (
    classify_variant_type,
    CohortVariantSummary,
    CONSEQUENCE_PRIORITY,
    DNA_BASES,
    GT_ALLELE_SEPARATOR_RE,
    SampleVariantSummary,
    VALID_FILTER_STATUS,
    VALID_VARIANT_TYPES,
    Variant,
    VariantAnalysisError,
    VariantAnalysisResult,
    VariantComparison,
    VariantEvidence,
    VariantFilter,
    VariantFilterResult,
    VariantStatistics,
)
from analyses._variant.vcf import (
    _coerce_info_value,
    _field_number,
    format_info_field,
    genotype_category,
    _genotype_index,
    _HEADER_FIELD_RE,
    _normalize_filter_status,
    parse_float,
    parse_genotype,
    parse_info_field,
    parse_int,
    parse_vcf,
    parse_vcf_file,
    parse_vcf_header_numbers,
    _remap_sample_data_for_alt,
    RESERVED_FORMAT_NUMBER,
    RESERVED_INFO_NUMBER,
    _select_allele_specific_info,
    validate_allele,
    validate_variant,
    VariantRecord,
)
from analyses._variant.normalization import (
    _add_validation_flag,
    _annotate_contig_resolution,
    check_reference_alleles,
    _chr_prefix_partner,
    CONTIG_MATCHING_MODES,
    ContigResolution,
    left_normalize_variant,
    MITOCHONDRIAL_CONTIG_NAMES,
    normalize_variants,
    _reference_allele_matches,
    resolve_contig_names,
    trim_common_bases,
    validate_contig_matching,
    validate_contig_references,
    validate_reference_sequence,
)
from analyses._variant.annotation import (
    annotate_from_info,
    annotate_from_mutation_analysis,
    annotate_variant,
    attach_evidence,
    check_features_against_variants,
    consequence_priority,
    highest_priority_consequence,
    prioritize_variants,
    synchronize_variant_fields,
)
from analyses._variant.statistics import (
    calculate_cohort_summary,
    calculate_variant_statistics,
    compare_variant_sets,
    filter_variants,
    interval_variant_density,
    is_transition_variant,
    sample_variant_summary,
    variant_key_set,
    variant_passes_filter,
    variants_by_chromosome,
    variants_by_consequence,
    variants_by_gene,
    variants_in_interval,
)
from analyses._variant.export import (
    build_variant_report,
    format_variant,
    variants_to_tsv,
    variants_to_vcf,
    write_vcf,
)
from analyses._variant.workflow import analyze_variants, analyze_vcf_text
from analyses._variant.selftest import (
    assert_variant_analysis_self_test,
    run_variant_analysis_self_test,
)


# ============================================================================
# Interactive UI
# ============================================================================


# ============================================================================
# Public API
# ============================================================================

__all__ = [
    "VariantAnalysisError",
    "Variant",
    "VariantRecord",
    "VariantFilter",
    "VariantFilterResult",
    "VariantStatistics",
    "SampleVariantSummary",
    "CohortVariantSummary",
    "VariantComparison",
    "VariantAnalysisResult",
    "VariantEvidence",
    "parse_float",
    "parse_int",
    "validate_allele",
    "validate_variant",
    "classify_variant_type",
    "parse_info_field",
    "format_info_field",
    "parse_genotype",
    "genotype_category",
    "parse_vcf",
    "parse_vcf_file",
    "trim_common_bases",
    "left_normalize_variant",
    "normalize_variants",
    "CONTIG_MATCHING_MODES",
    "ContigResolution",
    "validate_contig_matching",
    "validate_contig_references",
    "validate_reference_sequence",
    "resolve_contig_names",
    "check_reference_alleles",
    "annotate_variant",
    "annotate_from_info",
    "synchronize_variant_fields",
    "variant_passes_filter",
    "filter_variants",
    "is_transition_variant",
    "calculate_variant_statistics",
    "sample_variant_summary",
    "calculate_cohort_summary",
    "variant_key_set",
    "compare_variant_sets",
    "variants_by_chromosome",
    "variants_by_gene",
    "variants_by_consequence",
    "consequence_priority",
    "prioritize_variants",
    "variants_in_interval",
    "interval_variant_density",
    "variants_to_vcf",
    "write_vcf",
    "variants_to_tsv",
    "analyze_variants",
    "analyze_vcf_text",
    "format_variant",
    "build_variant_report",
    "annotate_from_mutation_analysis",
    "attach_evidence",
    "highest_priority_consequence",
    "run_variant_analysis_self_test",
    "assert_variant_analysis_self_test",
]


if __name__ == "__main__":
    from ui import style  # only when run as a script
    tests = run_variant_analysis_self_test()
    passed = sum(tests.values())
    style.header('Variant analysis self-test')
    for name, success in tests.items():
        print(f"[{'PASS' if success else 'FAIL'}] {name}")
    print(minor_rule())
    print(f"{passed}/{len(tests)} tests passed.")
    if passed != len(tests):
        raise SystemExit(1)
