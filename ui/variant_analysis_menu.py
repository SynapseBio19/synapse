"""Terminal dialogs of analyses/variant_analysis.py (moved here in P2-2).

The scientific functions stay in ``analyses.variant_analysis``; this module draws the
screens, reads the answers and calls them.
"""

from __future__ import annotations

from typing import Dict
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple
import math
from analyses.variant_analysis import (  # noqa: F401
    Variant,
    VariantAnalysisError,
    VariantAnalysisResult,
    VariantFilter,
    VariantFilterResult,
    analyze_variants,
    build_variant_report,
    calculate_cohort_summary,
    calculate_variant_statistics,
    check_reference_alleles,
    compare_variant_sets,
    filter_variants,
    format_variant,
    normalize_variants,
    parse_float,
    parse_int,
    parse_vcf,
    parse_vcf_file,
    run_variant_analysis_self_test,
    synchronize_variant_fields,
    variants_to_tsv,
    write_vcf,
)
from ui import style
from ui.style import minor_rule


def print_variant_report(
    result: VariantAnalysisResult,
    *,
    max_variants: Optional[int] = 100,
) -> None:
    print(build_variant_report(result, max_variants=max_variants))


def _prompt_variant() -> Variant:
    print()
    chromosome = style.ask('Chromosome/contig')
    position = parse_int(style.ask('Position'))
    if position is None:
        raise VariantAnalysisError("Position must be an integer.")
    reference = style.ask('Reference allele')
    alternative = style.ask('Alternative allele')
    quality = parse_float(style.ask('Quality (Enter for unknown)'))

    variant = Variant(
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternative=alternative,
        quality=quality,
        filter_status="UNFILTERED",
    )
    return variant


def _prompt_filter() -> VariantFilter:
    style.header('VARIANT FILTER')

    quality_text = style.ask('Minimum QUAL (Enter for none)')
    minimum_quality = parse_float(quality_text)

    max_af_text = style.ask('Maximum AF (Enter for none)')
    maximum_af = parse_float(max_af_text)

    require_pass = style.confirm("Require FILTER=PASS?")
    include_unknown_quality = style.confirm("Include variants with unknown QUAL?")
    include_unknown_frequency = style.confirm("Include variants with unknown AF?")

    return VariantFilter(
        minimum_quality=minimum_quality,
        maximum_allele_frequency=maximum_af,
        require_pass=require_pass,
        include_unknown_quality=include_unknown_quality,
        include_unknown_frequency=include_unknown_frequency,
    )


def _prompt_filename(default_name: str, extension: str) -> str:
    name = style.ask(f"Output filename [{default_name}.{extension}]")
    if not name:
        name = f"{default_name}.{extension}"
    elif not name.lower().endswith(f".{extension}"):
        name += f".{extension}"
    return name


def _print_filter_result(result: VariantFilterResult) -> None:
    print()
    style.field('Passed', f"{len(result.passed)}")
    style.field('Failed', f"{len(result.failed)}")
    if result.failed:
        style.section('FILTER FAILURES')
        for variant in result.failed:
            reasons = result.reasons.get(variant.key, [])
            print(f"{format_variant(variant)}")
            for reason in reasons:
                print(f"  - {reason}")


# Reference-only records (ALT ".") of the last load, reported after loading.
_REFERENCE_ONLY: List[str] = []


def _load_vcf_text_interactive() -> Tuple[List[Variant], List[str], List[str]]:
    print()
    path = style.ask('VCF file path (.vcf or .vcf.gz; Enter to paste text instead)').strip('"')
    if path:
        try:
            return parse_vcf_file(path, skipped=_REFERENCE_ONLY)
        except OSError as error:
            raise VariantAnalysisError(f"Could not read VCF file: {error}") from error
    print("Paste VCF text. Enter an empty line to finish.")
    lines: List[str] = []
    while True:
        line = input()
        if line == "":
            break
        lines.append(line)
    if not lines:
        raise VariantAnalysisError("No VCF text supplied.")
    return parse_vcf("\n".join(lines), skipped=_REFERENCE_ONLY)


def _prompt_contig_reference(pipeline, variants: Sequence[Variant]) -> Dict[str, str]:
    """Optionally map one VCF contig to the sequence loaded in Synapse."""
    context = getattr(pipeline, "context", None)
    sequence = getattr(context, "sequence", "") if context is not None else ""
    if not sequence:
        style.note("No Synapse sequence loaded: variants are trimmed only (no REF check, no repeat shift).")
        return {}
    contigs = sorted({v.chromosome for v in variants})
    header = getattr(context, "header", "") or "loaded sequence"
    print()
    style.field('Loaded Synapse sequence', f"{header} ({len(sequence)} bp)")
    style.field('Contigs in the VCF', f"{', '.join(contigs)}")
    style.note("VCF positions must be 1-based positions on this exact sequence.")
    choice = style.ask('Contig to map onto the loaded sequence (Enter = none)')
    if not choice:
        return {}
    if choice not in contigs:
        print(f"Contig '{choice}' is not in the VCF; no mapping used.")
        return {}
    return {choice: sequence}


def run_variant_analysis_ui(
    variants: Optional[Sequence[Variant]] = None,
    *,
    pipeline=None,
    name: str = "Variant Analysis",
) -> Optional[VariantAnalysisResult]:
    """Interactive CLI for the Variant Analysis module.

    The UI keeps the original dataset intact during filtering and normalization,
    supports real file exports, provides set comparison through VCF text input,
    and refreshes the analysis result after state-changing operations.
    """
    working = synchronize_variant_fields(list(variants or []))
    last_result: Optional[VariantAnalysisResult] = None
    last_filter_result: Optional[VariantFilterResult] = None
    metadata: List[str] = []
    sample_names: List[str] = sorted({s for v in working for s in v.sample_data})
    filtered_view: Optional[List[Variant]] = None

    def refresh_result() -> Optional[VariantAnalysisResult]:
        nonlocal last_result
        if not working:
            last_result = None
            return None
        last_result = analyze_variants(working, name=name)
        if pipeline is not None:
            _store_variant_pipeline_result(pipeline, "variant", last_result.to_dict())
        return last_result

    while True:
        style.header('VARIANT ANALYSIS')
        style.field('Loaded variants', f"{len(working)}")
        if filtered_view is not None:
            style.field('Current filter view', f"{len(filtered_view)} passed")
        style.menu([
            ('Input', [('1', 'Load VCF'), ('2', 'Add Variant Manually')]),
            ('Processing', [('3', 'Normalize Variants'), ('4', 'Filter Variants')]),
            ('Analysis', [
                ('5', 'Variant Statistics'),
                ('6', 'Analyze Variants'),
                ('7', 'Compare Variant Sets'),
                ('8', 'Cohort / Sample Analysis'),
            ]),
            ('Output', [('9', 'Export VCF'), ('10', 'Export TSV'), ('11', 'Full Variant Report')]),
            ('System', [('12', 'Run Variant Self-Test')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return last_result

        if choice == "1":
            try:
                _REFERENCE_ONLY.clear()
                parsed, parsed_metadata, parsed_samples = _load_vcf_text_interactive()
                working = synchronize_variant_fields(parsed)
                metadata = parsed_metadata
                sample_names = parsed_samples
                filtered_view = None
                last_filter_result = None
                refresh_result()
                print(f"Loaded {len(working)} variant(s).")
                if _REFERENCE_ONLY:
                    style.note(f"{len(_REFERENCE_ONLY)} record(s) without an ALT allele (ALT '.', "
                               "reference-only sites) are not variants and were not loaded.")
            except (VariantAnalysisError, ValueError) as error:
                style.error(f"{error}")

        elif choice == "2":
            try:
                variant = _prompt_variant()
                working.append(variant)
                working = synchronize_variant_fields(working)
                filtered_view = None
                last_filter_result = None
                refresh_result()
                style.field('Added', f"{format_variant(variant)}")
            except VariantAnalysisError as error:
                style.error(f"{error}")

        elif choice == "6":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            try:
                result = refresh_result()
                if result is not None:
                    print()
                    print_variant_report(result)
            except VariantAnalysisError as error:
                style.error(f"{error}")

        elif choice == "4":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            try:
                criteria = _prompt_filter()
                last_filter_result = filter_variants(working, criteria)
                filtered_view = list(last_filter_result.passed)
                _print_filter_result(last_filter_result)
            except VariantAnalysisError as error:
                style.error(f"{error}")

        elif choice == "7":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            try:
                print()
                print("Paste the SECOND VCF set. Enter an empty line to finish.")
                lines: List[str] = []
                while True:
                    line = input()
                    if line == "":
                        break
                    lines.append(line)
                if not lines:
                    raise VariantAnalysisError("No second VCF set supplied.")
                right, _, _ = parse_vcf("\n".join(lines))
                comparison = compare_variant_sets(working, right)
                style.header('VARIANT SET COMPARISON')
                style.field('Left variants', f"{comparison.left_count}")
                style.field('Right variants', f"{comparison.right_count}")
                style.field('Shared variants', f"{comparison.shared_count}")
                style.field('Jaccard similarity', f"{comparison.jaccard_similarity:.4f}")
                print()
                print("Only in current set:")
                for v in comparison.only_in_left:
                    print(f"  - {format_variant(v)}")
                print("Only in second set:")
                for v in comparison.only_in_right:
                    print(f"  - {format_variant(v)}")
            except (VariantAnalysisError, ValueError) as error:
                style.error(f"{error}")

        elif choice == "5":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            stats = calculate_variant_statistics(working)
            style.header('VARIANT STATISTICS')
            style.field('Total', f"{stats.total_variants}")
            style.field('Unique', f"{stats.unique_variants}")
            style.field('SNVs', f"{stats.snvs}")
            style.field('MNVs', f"{stats.mnvs}")
            style.field('Insertions', f"{stats.insertions}")
            style.field('Deletions', f"{stats.deletions}")
            style.field('Indels', f"{stats.indels}")
            style.field('Missing QUAL', f"{stats.missing_quality}")
            style.field('Mean QUAL', f"{'.' if stats.mean_quality is None else f'{stats.mean_quality:.2f}'}")
            style.field('Mean AF', f"{'.' if stats.mean_allele_frequency is None else f'{stats.mean_allele_frequency:.5f}'}")
            style.field('Median AF', f"{'.' if stats.median_allele_frequency is None else f'{stats.median_allele_frequency:.5f}'}")
            style.field('Rare variants', f"{stats.rare_variants}")
            style.field('Genes affected', f"{stats.genes_affected}")
            style.field('Chromosomes', f"{stats.chromosomes_affected}")
            style.field('Transitions', f"{stats.transition_count}")
            style.field('Transversions', f"{stats.transversion_count}")
            ratio = stats.transition_transversion_ratio
            style.field('Ti/Tv ratio', f"{'inf' if ratio == math.inf else '.' if ratio is None else f'{ratio:.3f}'}")

        elif choice == "8":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            cohort = calculate_cohort_summary(working)
            style.header('COHORT / SAMPLE ANALYSIS')
            style.field('Samples', f"{cohort.sample_count}")
            style.field('Shared variants', f"{cohort.shared_variant_count}")
            style.field('Private variants', f"{cohort.private_variant_count}")
            for sample in cohort.sample_summaries:
                print(
                    f"{sample.sample}: called={sample.called_variants}, "
                    f"het={sample.heterozygous_calls}, "
                    f"hom_alt={sample.homozygous_alt_calls}, "
                    f"missing={sample.missing_calls}, "
                    f"alt_alleles={sample.allele_count}"
                )

        elif choice == "3":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            try:
                contig_map = _prompt_contig_reference(pipeline, working)
                if contig_map:
                    report = check_reference_alleles(working, contig_map)
                    print(
                        f"REF check against the loaded sequence: {report['checked']} checked, "
                        f"{report['mismatched']} mismatch(es), {report['not_checked']} not checked."
                    )
                working = normalize_variants(working, contig_references=contig_map or {})
                working = synchronize_variant_fields(working)
                context = getattr(pipeline, "context", None)
                if contig_map and context is not None and hasattr(context, "set_mapped_variants"):
                    contig = next(iter(contig_map))
                    context.set_mapped_variants(contig, working)
                    style.note(f"The variants of {contig} are kept for this sequence; later primer "
                               "and guide designs on it are checked against them.")
                filtered_view = None
                last_filter_result = None
                refresh_result()
                without_reference = sum(
                    "REFERENCE_CONTIG_NOT_FOUND" in v.annotation.get("validation_flags", [])
                    for v in working
                )
                print(
                    f"Normalized {len(working)} variant(s): "
                    f"{len(working) - without_reference} with a reference sequence "
                    "(REF checked, left-aligned unless flagged), "
                    f"{without_reference} without one (common bases trimmed only)."
                )
            except VariantAnalysisError as error:
                style.error(f"{error}")

        elif choice == "9":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            try:
                filename = _prompt_filename("synapse_variants", "vcf")
                write_vcf(working, filename, metadata=metadata, sample_names=sample_names)
                style.success(f"VCF exported: {filename}")
            except OSError as error:
                style.error(f"could not write the VCF file: {error}")

        elif choice == "10":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            try:
                filename = _prompt_filename("synapse_variants", "tsv")
                with open(filename, "w", encoding="utf-8", newline="") as handle:
                    handle.write(variants_to_tsv(working))
                style.success(f"TSV exported: {filename}")
            except OSError as error:
                style.error(f"could not write the TSV file: {error}")

        elif choice == "11":
            if not working:
                style.note("No variants loaded.")
                style.pause()
                continue
            result = refresh_result()
            if result is not None:
                print()
                print_variant_report(result)

        elif choice == "12":
            tests = run_variant_analysis_self_test()
            style.header('VARIANT ANALYSIS SELF-TEST')
            for test_name, passed in tests.items():
                print(f"[{'PASS' if passed else 'FAIL'}] {test_name}")
            print(minor_rule())
            print(f"{sum(tests.values())}/{len(tests)} tests passed.")

        else:
            style.invalid_selection()


def _store_variant_pipeline_result(pipeline, key: str, value: object) -> None:
    """Store a result without making the pipeline a hard dependency."""
    import logging
    logger = logging.getLogger("synapse.variant")
    if pipeline is None:
        return
    try:
        if hasattr(pipeline, "store") and callable(pipeline.store):
            pipeline.store(key, value)
        elif hasattr(pipeline, "results") and isinstance(pipeline.results, dict):
            pipeline.results[key] = value
        elif isinstance(pipeline, dict):
            pipeline[key] = value
        else:
            logger.warning("Pipeline object has no result store; %s result not stored.", key)
    except Exception:
        # Analysis success must not be lost because pipeline integration
        # differs, but the failure is logged (audit Phase V).
        logger.exception("Could not store %s result in the pipeline.", key)
