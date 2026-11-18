"""VCF and TSV export and the text report."""

from __future__ import annotations

import csv
import io
import math
from core.textformat import field_line, major_rule, minor_rule
from typing import Optional, Sequence

from analyses._variant.models import Variant, VariantAnalysisResult
from analyses._variant.vcf import format_info_field


# ============================================================================
# VCF serialization
# ============================================================================

def variants_to_vcf(
    variants: Sequence[Variant],
    *,
    metadata: Optional[Sequence[str]] = None,
    sample_names: Optional[Sequence[str]] = None,
) -> str:
    """Serialize variants into a valid basic VCF representation."""
    lines = [
        "##fileformat=VCFv4.3",
    ]

    if metadata:
        lines.extend(
            line if line.startswith("##") else f"##{line}"
            for line in metadata
            if line
            and not line.startswith("#CHROM")
            and line != "##fileformat=VCFv4.3"
        )

    samples = list(sample_names or sorted({
        sample
        for variant in variants
        for sample in variant.sample_data
    }))

    header = [
        "#CHROM", "POS", "ID", "REF", "ALT",
        "QUAL", "FILTER", "INFO",
    ]

    if samples:
        header.extend(["FORMAT", *samples])

    lines.append("\t".join(header))

    for variant in sorted(
        variants,
        key=lambda v: (v.chromosome, v.position, v.reference, v.alternative),
    ):
        row = [
            variant.chromosome,
            str(variant.position),
            variant.variant_id or ".",
            variant.reference,
            variant.alternative,
            "." if variant.quality is None else str(variant.quality),
            "." if variant.filter_status in (None, "", "UNFILTERED") else variant.filter_status,
            format_info_field(variant.info),
        ]

        if samples:
            format_keys = variant.format_keys or ["GT"]
            row.append(":".join(format_keys))
            for sample in samples:
                values = variant.sample_data.get(sample, {})
                row.append(
                    ":".join(str(values.get(key, ".")) for key in format_keys)
                )

        lines.append("\t".join(row))

    return "\n".join(lines) + "\n"


def write_vcf(
    variants: Sequence[Variant],
    path: str,
    *,
    metadata: Optional[Sequence[str]] = None,
    sample_names: Optional[Sequence[str]] = None,
) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(
            variants_to_vcf(
                variants,
                metadata=metadata,
                sample_names=sample_names,
            )
        )


def variants_to_tsv(variants: Sequence[Variant]) -> str:
    """Export a normalized analysis table."""
    columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "TYPE",
        "QUAL", "FILTER", "GENE", "TRANSCRIPT", "REGION",
        "CONSEQUENCE", "CODING_CHANGE", "PROTEIN_CHANGE",
        "AF", "AC", "AN",
    ]
    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)

    for v in variants:
        writer.writerow([
            v.chromosome,
            v.position,
            v.variant_id or ".",
            v.reference,
            v.alternative,
            v.variant_type,
            "." if v.quality is None else v.quality,
            v.filter_status,
            v.gene or ".",
            v.transcript or ".",
            v.region or ".",
            v.consequence or ".",
            v.coding_change or ".",
            v.protein_change or ".",
            "." if v.allele_frequency is None else v.allele_frequency,
            "." if v.allele_count is None else v.allele_count,
            "." if v.allele_number is None else v.allele_number,
        ])

    return output.getvalue()


# ============================================================================
# Reporting
# ============================================================================

def format_variant(variant: Variant) -> str:
    quality = "." if variant.quality is None else f"{variant.quality:.2f}"
    return (
        f"{variant.chromosome}:{variant.position} "
        f"{variant.reference}>{variant.alternative} "
        f"[{variant.variant_type}] "
        f"QUAL={quality}"
    )


def build_variant_report(
    result: VariantAnalysisResult,
    *,
    max_variants: Optional[int] = 100,
) -> str:
    s = result.statistics

    lines = [
        major_rule(),
        "SYNAPSE VARIANT ANALYSIS",
        major_rule(),
        field_line('Analysis', f"{result.name}"),
        field_line('Total variants', f"{s.total_variants}"),
        field_line('Unique variants', f"{s.unique_variants}"),
        "",
        "VARIANT TYPES",
        minor_rule(),
        field_line('SNVs', f"{s.snvs}"),
        field_line('MNVs', f"{s.mnvs}"),
        field_line('Insertions', f"{s.insertions}"),
        field_line('Deletions', f"{s.deletions}"),
        field_line('Indels', f"{s.indels}"),
        field_line('Complex', f"{s.complex}"),
        field_line('Symbolic', f"{s.symbolic}"),
        "",
        "QUALITY",
        minor_rule(),
        field_line('PASS', f"{s.passed_filters}"),
        field_line('Failed filters', f"{s.failed_filters}"),
        field_line('Not filtered (.)', f"{s.unfiltered}"),
        field_line('Missing QUAL', f"{s.missing_quality}"),
        field_line('Mean QUAL', f"{'.' if s.mean_quality is None else f'{s.mean_quality:.2f}'}"),
        field_line('Median QUAL', f"{'.' if s.median_quality is None else f'{s.median_quality:.2f}'}"),
        field_line('Min QUAL', f"{'.' if s.min_quality is None else f'{s.min_quality:.2f}'}"),
        field_line('Max QUAL', f"{'.' if s.max_quality is None else f'{s.max_quality:.2f}'}"),
        "",
        "ALLELE FREQUENCY",
        minor_rule(),
        field_line('Mean AF', f"{'.' if s.mean_allele_frequency is None else f'{s.mean_allele_frequency:.5f}'}"),
        field_line('Median AF', f"{'.' if s.median_allele_frequency is None else f'{s.median_allele_frequency:.5f}'}"),
        field_line('Rare variants', f"{s.rare_variants}"),
        "",
        "ANNOTATION",
        minor_rule(),
        field_line('Genes affected', f"{s.genes_affected}"),
        field_line('Chromosomes', f"{s.chromosomes_affected}"),
    ]

    if s.consequence_counts:
        lines.extend(["", "CONSEQUENCES", minor_rule()])
        for consequence, count in sorted(
            s.consequence_counts.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            lines.append(f"{consequence:<25} {count}")

    ratio = (
        "inf" if s.transition_transversion_ratio == math.inf
        else "N/A" if s.transition_transversion_ratio is None
        else f"{s.transition_transversion_ratio:.3f}"
    )
    lines.extend([
        "",
        "SUBSTITUTION SPECTRUM",
        minor_rule(),
        field_line('Transitions', f"{s.transition_count}"),
        field_line('Transversions', f"{s.transversion_count}"),
        field_line('Ti/Tv ratio', f"{ratio}"),
        "",
        "VARIANTS",
        minor_rule(),
    ])

    listed = result.variants if max_variants is None else result.variants[:max_variants]
    for index, variant in enumerate(listed, start=1):
        lines.append(f"{index:>4}. {format_variant(variant)}")
        if variant.gene:
            lines.append(f"       Gene: {variant.gene}")
        if variant.consequence:
            lines.append(f"       Consequence: {variant.consequence}")
        if variant.protein_change:
            lines.append(f"       Protein: {variant.protein_change}")

    if max_variants is not None and len(result.variants) > max_variants:
        lines.append(
            f"... {len(result.variants) - max_variants} additional variants omitted."
        )

    if result.filter_result is not None:
        lines.extend([
            "",
            "FILTER RESULT",
            minor_rule(),
            field_line('Passed', f"{len(result.filter_result.passed)}"),
            field_line('Failed', f"{len(result.filter_result.failed)}"),
        ])

    if result.cohort_summary:
        c = result.cohort_summary
        lines.extend([
            "",
            "COHORT",
            minor_rule(),
            field_line('Samples', f"{c.sample_count}"),
            field_line('Shared variants', f"{c.shared_variant_count}"),
            field_line('Private variants', f"{c.private_variant_count}"),
        ])

    if result.notes:
        lines.extend(["", "NOTES", minor_rule()])
        lines.extend(f"- {note}" for note in result.notes)

    lines.append(major_rule())
    return "\n".join(lines)
