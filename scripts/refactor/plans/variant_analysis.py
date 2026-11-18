"""Split plan for analyses/variant_analysis.py (Phase 9)."""

SOURCE = "analyses/variant_analysis.py"
PACKAGE = "analyses._variant"
PACKAGE_DIR = "analyses/_variant"
PACKAGE_DOC = """Implementation of :mod:`analyses.variant_analysis` (split out in Phase 9).

Import from ``analyses.variant_analysis``; it re-exports every name defined
here, keeps ``__all__`` and is the only public interface. The modules, in
dependency order:

- ``models``: constants, the dataclasses, ``VariantAnalysisError`` and
  ``classify_variant_type`` (used by ``Variant``)
- ``vcf``: VCF parsing (header numbers, INFO, genotypes, ``VariantRecord``
  and the multi-allelic split)
- ``normalization``: trimming, left normalization, contig matching,
  reference checks
- ``annotation``: annotation from INFO and from the mutation analysis,
  evidence, prioritization, feature overlap
- ``statistics``: filters, statistics, cohort summary, set comparison,
  grouping and intervals
- ``export``: VCF and TSV export, the text report
- ``workflow``: ``analyze_variants`` and ``analyze_vcf_text``
- ``selftest``: the built-in self-test
"""

MODULES = {
    "models": ("Constants, result dataclasses and the module error.", [
        "DNA_BASES", "VALID_FILTER_STATUS", "VALID_VARIANT_TYPES", "CONSEQUENCE_PRIORITY",
        "GT_ALLELE_SEPARATOR_RE", "VariantAnalysisError", "classify_variant_type", "Variant", "VariantFilter",
        "VariantFilterResult", "VariantStatistics", "SampleVariantSummary", "CohortVariantSummary",
        "VariantComparison", "VariantAnalysisResult", "VariantEvidence",
    ]),
    "vcf": ("VCF parsing: header numbers, INFO and FORMAT fields, genotypes, multi-allelic split.", [
        "RESERVED_INFO_NUMBER", "RESERVED_FORMAT_NUMBER", "_HEADER_FIELD_RE", "parse_vcf_header_numbers",
        "_field_number", "_genotype_index", "_select_allele_specific_info", "_remap_sample_data_for_alt",
        "VariantRecord", "parse_float", "parse_int", "validate_allele", "validate_variant",
        "parse_info_field", "_coerce_info_value", "format_info_field", "parse_genotype", "genotype_category",
        "_normalize_filter_status", "parse_vcf", "parse_vcf_file",
    ]),
    "normalization": ("Allele trimming, left normalization, contig matching and reference checks.", [
        "trim_common_bases", "_reference_allele_matches", "left_normalize_variant", "CONTIG_MATCHING_MODES",
        "MITOCHONDRIAL_CONTIG_NAMES", "ContigResolution", "validate_contig_matching",
        "validate_reference_sequence", "validate_contig_references", "_chr_prefix_partner",
        "resolve_contig_names", "_add_validation_flag", "_annotate_contig_resolution", "normalize_variants",
        "check_reference_alleles",
    ]),
    "annotation": ("Annotation, evidence, consequence priority and feature overlap.", [
        "annotate_variant", "annotate_from_info", "synchronize_variant_fields", "consequence_priority",
        "prioritize_variants", "check_features_against_variants", "annotate_from_mutation_analysis",
        "attach_evidence", "highest_priority_consequence",
    ]),
    "statistics": ("Filters, statistics, cohort summary, set comparison, grouping and intervals.", [
        "variant_passes_filter", "filter_variants", "is_transition_variant", "calculate_variant_statistics",
        "sample_variant_summary", "calculate_cohort_summary", "variant_key_set", "compare_variant_sets",
        "variants_by_chromosome", "variants_by_gene", "variants_by_consequence", "variants_in_interval",
        "interval_variant_density",
    ]),
    "export": ("VCF and TSV export and the text report.", [
        "variants_to_vcf", "write_vcf", "variants_to_tsv", "format_variant", "build_variant_report",
    ]),
    "workflow": ("``analyze_variants`` and ``analyze_vcf_text``.", [
        "analyze_variants", "analyze_vcf_text",
    ]),
    "selftest": ("The built-in self-test of the variant analysis (menu entry and test suite).", [
        "run_variant_analysis_self_test", "assert_variant_analysis_self_test",
    ]),
}

FACADE_NAMES = {"__all__"}
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 9: the implementation moved to the internal package\n"
    "# ``analyses._variant`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface: it re-exports every name it\n"
    "# defined before and keeps ``__all__``.\n"
)
