"""The built-in self-test of the variant analysis (menu entry and test suite)."""

from __future__ import annotations

from typing import Dict

from analyses._variant.models import classify_variant_type, Variant, VariantEvidence, VariantFilter
from analyses._variant.vcf import (
    genotype_category,
    parse_info_field,
    parse_vcf,
    validate_variant,
    VariantRecord,
)
from analyses._variant.normalization import left_normalize_variant
from analyses._variant.annotation import (
    annotate_from_info,
    annotate_variant,
    attach_evidence,
    prioritize_variants,
)
from analyses._variant.statistics import (
    calculate_cohort_summary,
    calculate_variant_statistics,
    compare_variant_sets,
    filter_variants,
    variants_in_interval,
)
from analyses._variant.export import variants_to_tsv, variants_to_vcf
from analyses._variant.workflow import analyze_variants


# ============================================================================
# Self tests
# ============================================================================

def run_variant_analysis_self_test() -> Dict[str, bool]:
    tests: Dict[str, bool] = {}

    # Basic classification.
    tests["snv_classification"] = (
        classify_variant_type("A", "G") == "SNV"
    )
    tests["mnv_classification"] = (
        classify_variant_type("AT", "GC") == "MNV"
    )
    tests["insertion_classification"] = (
        classify_variant_type("A", "ATG") == "INSERTION"
    )
    tests["deletion_classification"] = (
        classify_variant_type("ATG", "A") == "DELETION"
    )

    # INFO.
    info = parse_info_field("DP=20;AF=0.25;DB")
    tests["info_parser"] = (
        info["DP"] == 20 and
        info["AF"] == 0.25 and
        info["DB"] is True
    )

    # Genotypes.
    tests["genotype_parser"] = (
        genotype_category("0/1") == "heterozygous"
        and genotype_category("1/1") == "homozygous_alt"
        and genotype_category("0/0") == "homozygous_ref"
        and genotype_category("./.") == "missing"
    )

    # VCF parsing.
    vcf = """##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\trs1\tA\tG\t50\tPASS\tAF=0.1\tGT:DP\t0/1:20\t1/1:30
chr1\t20\t.\tC\tT,A\t40\tPASS\tDP=10\tGT\t0/0\t0/1
"""
    parsed, _, samples = parse_vcf(vcf)
    tests["vcf_parser"] = (
        len(parsed) == 3
        and samples == ["S1", "S2"]
        and parsed[0].variant_type == "SNV"
    )

    # Multi-allelic splitting.
    record = VariantRecord(
        "chr1", 100, "A", ["C", "G"], variant_id="multi"
    )
    tests["multiallelic_split"] = len(record.split()) == 2

    # Normalization.
    v = Variant("chr1", 5, "AA", "AT")
    normalized = left_normalize_variant(v)
    tests["normalization"] = normalized.normalized and normalized.position >= 1

    # Filtering.
    variants = [
        Variant("chr1", 1, "A", "G", quality=50, filter_status="PASS"),
        Variant("chr1", 2, "C", "T", quality=5, filter_status="PASS"),
    ]
    filtered = filter_variants(
        variants,
        VariantFilter(minimum_quality=20),
    )
    tests["filtering"] = (
        len(filtered.passed) == 1
        and len(filtered.failed) == 1
    )

    # Statistics.
    stats = calculate_variant_statistics(variants)
    tests["statistics"] = (
        stats.total_variants == 2
        and stats.snvs == 2
        and stats.transition_count == 2
    )

    # Set comparison.
    left = [Variant("chr1", 1, "A", "G"), Variant("chr1", 2, "C", "T")]
    right = [Variant("chr1", 2, "C", "T"), Variant("chr1", 3, "G", "A")]
    comparison = compare_variant_sets(left, right)
    tests["set_comparison"] = (
        comparison.shared_count == 1
        and comparison.jaccard_similarity == 1 / 3
    )

    # Annotation.
    annotated = annotate_variant(
        Variant("chr1", 10, "A", "G"),
        gene="GENE1",
        consequence="missense",
        allele_frequency=0.005,
    )
    tests["annotation"] = (
        annotated.gene == "GENE1"
        and annotated.consequence == "missense"
        and annotated.allele_frequency == 0.005
    )

    # Cohort.
    cohort_variants = [
        Variant(
            "chr1", 1, "A", "G",
            sample_data={
                "S1": {"GT": "0/1"},
                "S2": {"GT": "0/0"},
            },
        ),
        Variant(
            "chr1", 2, "C", "T",
            sample_data={
                "S1": {"GT": "1/1"},
                "S2": {"GT": "0/1"},
            },
        ),
    ]
    cohort = calculate_cohort_summary(cohort_variants)
    tests["cohort_analysis"] = (
        cohort.sample_count == 2
        and cohort.unique_variant_count == 2
        and cohort.shared_variant_count == 1
    )

    # Prioritization.
    prioritized = prioritize_variants([
        Variant("chr1", 1, "A", "G", consequence="synonymous"),
        Variant("chr1", 2, "A", "T", consequence="frameshift"),
    ])
    tests["prioritization"] = prioritized[0].consequence == "frameshift"

    # Interval.
    interval = variants_in_interval(
        [Variant("chr1", 10, "A", "G"), Variant("chr1", 50, "C", "T")],
        "chr1", 1, 20,
    )
    tests["interval_analysis"] = len(interval) == 1

    # Serialization.
    exported = variants_to_vcf([Variant("chr1", 1, "A", "G")])
    reparsed, _, _ = parse_vcf(exported)
    tests["vcf_roundtrip"] = (
        len(reparsed) == 1
        and reparsed[0].key == ("chr1", 1, "A", "G")
    )

    tsv = variants_to_tsv([Variant("chr1", 1, "A", "G")])
    tests["tsv_export"] = "CHROM" in tsv and "chr1" in tsv

    # Evidence.
    evidence_variant = attach_evidence(
        Variant("chr1", 1, "A", "G"),
        VariantEvidence(
            source="example",
            field="frequency",
            value=0.01,
            evidence_level="observed",
        ),
    )
    tests["evidence_container"] = (
        len(evidence_variant.annotation["evidence"]) == 1
    )

    # High-level workflow.
    workflow = analyze_variants(
        variants,
        filter_criteria=VariantFilter(minimum_quality=10),
    )
    tests["high_level_workflow"] = (
        workflow.count == 2
        and workflow.filter_result is not None
        and len(workflow.filter_result.passed) == 1
    )

    # JSON.
    tests["json_export"] = '"variants"' in workflow.to_json()

    # Validation.
    invalid = Variant("chr1", 1, "A", "G", allele_frequency=1.5)
    tests["validation"] = bool(validate_variant(invalid))

    # Integration: VCF INFO AF must feed filtering and statistics directly.
    af_vcf = """##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t100\tPASS\tAF=0.10
chr1\t200\t.\tC\tT\t100\tPASS\tAF=0.01
"""
    af_variants, _, _ = parse_vcf(af_vcf)
    af_filter = filter_variants(
        af_variants,
        VariantFilter(maximum_allele_frequency=0.05),
    )
    af_stats = calculate_variant_statistics(af_variants)
    tests["af_integration"] = (
        len(af_filter.passed) == 1
        and len(af_filter.failed) == 1
        and af_stats.mean_allele_frequency == 0.055
    )

    # Missing values are excluded when a threshold is explicitly requested.
    missing_filter = filter_variants(
        [Variant("chr1", 1, "A", "G")],
        VariantFilter(minimum_quality=10),
    )
    tests["missing_value_filtering"] = (
        len(missing_filter.passed) == 0
        and len(missing_filter.failed) == 1
    )

    # Multiallelic AF and genotype splitting must be allele-specific.
    multi_vcf = """##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t300\t.\tA\tC,G\t50\tPASS\tAF=0.1,0.2\tGT:AD\t1/2:10,5,7
"""
    multi_variants, _, _ = parse_vcf(multi_vcf)
    tests["multiallelic_allele_specific"] = (
        len(multi_variants) == 2
        and multi_variants[0].info.get("AF") == 0.1
        and multi_variants[1].info.get("AF") == 0.2
        and multi_variants[0].sample_data["S1"]["GT"] == "1/0"
        and multi_variants[1].sample_data["S1"]["GT"] == "0/1"
        and multi_variants[0].sample_data["S1"]["AD"] == "10,5"
        and multi_variants[1].sample_data["S1"]["AD"] == "10,7"
    )

    # Normalization must preserve evidence/annotation state.
    annotated_original = annotate_variant(
        Variant("chr1", 5, "AA", "AT"),
        gene="GENE1",
        consequence="missense",
    )
    attach_evidence(
        annotated_original,
        VariantEvidence(source="test", field="note", value="preserve"),
    )
    normalized_annotated = left_normalize_variant(annotated_original)
    tests["normalization_preserves_annotation"] = (
        normalized_annotated.gene == "GENE1"
        and normalized_annotated.consequence == "missense"
        and len(normalized_annotated.annotation.get("evidence", [])) == 1
    )

    # Round-trip must retain INFO AF.
    roundtrip_source = [Variant("chr1", 400, "G", "A", info={"AF": 0.25})]
    roundtrip_text = variants_to_vcf(roundtrip_source)
    roundtrip_variants, _, _ = parse_vcf(roundtrip_text)
    annotate_from_info(roundtrip_variants)
    tests["af_roundtrip"] = (
        len(roundtrip_variants) == 1
        and roundtrip_variants[0].allele_frequency == 0.25
    )

    return tests


def assert_variant_analysis_self_test() -> None:
    tests = run_variant_analysis_self_test()
    failed = [name for name, passed in tests.items() if not passed]
    if failed:
        raise AssertionError(
            "Variant analysis self-test failed: " + ", ".join(failed)
        )
