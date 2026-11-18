"""
Synapse - Mutation Analysis
===========================

Professional sequence-level mutation analysis module.

Features
--------
- DNA validation and normalization
- IUPAC ambiguity support
- Reverse-complement utilities
- Needleman-Wunsch global alignment
- Mutation/event detection
- SNP / transition / transversion classification
- Insertions / deletions / complex substitutions
- Mutation statistics
- Coding consequence analysis
- Synonymous / missense / nonsense
- Start-loss / stop-loss
- Frameshift / in-frame indel detection
- Multiple-variant comparison
- Shared / unique mutation analysis
- Sequence quality metrics
- Local mutation context
- Mutation spectra
- Mutation hotspots
- Restriction-site impact
- Primer overlap
- CRISPR target overlap
- Mutation maps
- JSON serialization
- Interactive CLI
- Deterministic regression self-test

This module intentionally performs sequence-level analysis.
Biological pathogenicity, clinical significance, experimental success,
or functional effect must not be inferred from sequence-level descriptors
alone.
"""

# Phase 8: the implementation moved to the internal package
# ``analyses._mutation`` (one module per responsibility, see its __init__).
# This module remains the public interface: it re-exports every name it
# defined before, private helpers included, and keeps ``__all__``.

from __future__ import annotations

from core.textformat import minor_rule
from analyses._mutation.sequence import (
    AA_NAMES,
    amino_acid_name,
    CODON_TABLE,
    DNA_ALPHABET,
    gc_content,
    iupac_code,
    IUPAC_CODES,
    IUPAC_COMPLEMENT,
    IUPAC_FROM_BASES,
    iupac_matches,
    normalize_dna,
    normalize_iupac_sequence,
    PURINES,
    PYRIMIDINES,
    reverse_complement,
    reverse_complement_iupac,
    translate_dna,
    validate_dna_sequence,
)
from analyses._mutation.models import (
    AdvancedMutationAnalysis,
    AlignmentBlock,
    CodingAnalysis,
    Mutation,
    MutationAnalysisError,
    MutationAnalysisResult,
    MutationContext,
    MutationHotspot,
    MutationImpact,
    MutationSpectrum,
    MutationStatistics,
    RestrictionSiteImpact,
    SequenceCoordinate,
    SequencePromptCancelled,
    SequenceQuality,
)
from analyses._mutation.alignment import (
    alignment_identity,
    _banded_align,
    BANDED_ALIGNMENT_MIN_CELLS,
    _banded_edit_distance,
    _banded_is_exact,
    build_alignment_blocks,
    _edit_distance,
    _full_edit_distance,
    global_align,
    _nucleotide_match,
    sequence_similarity,
)
from analyses._mutation.detection import (
    classify_substitution,
    detect_mutations_from_alignment,
    _detect_simple_indel,
    is_transition,
    _make_mutation,
    _mutation_notation,
    normalize_indel_left,
)
from analyses._mutation.coding import (
    analyze_coding_consequences,
    annotate_mutations_with_coding_consequences,
    assess_coding_frame,
    classify_coding_substitution,
    CODING_FRAME_ASSUMPTION,
    compare_proteins,
    _first_difference_position,
    orf_reference_cds,
    _translate_canonical,
)
from analyses._mutation.hgvs import (
    hgvs_c_allele,
    mutation_hgvs_c,
    mutation_hgvs_p,
    _reference_span,
)
from analyses._mutation.comparison import (
    analyze_mutations,
    calculate_mutation_statistics,
    compare_multiple_variants,
    compare_sequences,
    detect_mutations,
    MAX_RECOMMENDED_ALIGNMENT_LENGTH,
    _mutation_key,
    shared_mutations,
    unique_mutations,
)
from analyses._mutation.context import (
    annotate_regions,
    annotate_sequence_regions,
    build_mutation_spectrum,
    extract_mutation_context,
    find_mutation_hotspots,
    sequence_complexity,
    sequence_quality,
)
from analyses._mutation.impact import (
    crispr_mutation_overlap,
    _find_all_occurrences,
    _interval_overlap,
    primer_mutation_overlap,
    restriction_site_impact,
)
from analyses._mutation.reporting import (
    advanced_mutation_report,
    build_mutation_report,
    create_advanced_analysis,
    format_mutation,
    mutation_map,
    mutation_report,
    mutation_signature,
    result_to_json,
    validate_mutation_result,
)
from analyses._mutation.selftest import assert_mutation_self_test, run_mutation_self_test


__all__ = [
    # Core
    "AlignmentBlock",
    "Mutation",
    "MutationStatistics",
    "CodingAnalysis",
    "MutationAnalysisResult",
    "MutationAnalysisError",

    # Sequence
    "normalize_dna",
    "validate_dna_sequence",
    "normalize_iupac_sequence",
    "gc_content",
    "reverse_complement",
    "reverse_complement_iupac",

    # Translation
    "translate_dna",
    "amino_acid_name",

    # Alignment
    "global_align",
    "alignment_identity",
    "sequence_similarity",
    "build_alignment_blocks",

    # Mutation classification
    "is_transition",
    "classify_substitution",
    "detect_mutations_from_alignment",
    "detect_mutations",
    "analyze_mutations",

    # Coding
    "classify_coding_substitution",
    "analyze_coding_consequences",
    "annotate_mutations_with_coding_consequences",

    # Statistics
    "calculate_mutation_statistics",

    # Comparison
    "compare_sequences",
    "compare_multiple_variants",
    "shared_mutations",
    "unique_mutations",

    # Regions
    "annotate_regions",
    "annotate_sequence_regions",

    # Quality
    "sequence_quality",
    "sequence_complexity",

    # Context / spectrum
    "MutationContext",
    "MutationSpectrum",
    "MutationHotspot",
    "extract_mutation_context",
    "build_mutation_spectrum",
    "find_mutation_hotspots",

    # IUPAC
    "IUPAC_CODES",
    "iupac_matches",
    "iupac_code",

    # Indels / HGVS
    "normalize_indel_left",
    "mutation_hgvs_c",
    "hgvs_c_allele",
    "mutation_hgvs_p",

    # Protein
    "compare_proteins",

    # Cross-module impact
    "RestrictionSiteImpact",
    "MutationImpact",
    "restriction_site_impact",
    "primer_mutation_overlap",
    "crispr_mutation_overlap",

    # Advanced
    "AdvancedMutationAnalysis",
    "create_advanced_analysis",
    "advanced_mutation_report",

    # Validation / output
    "validate_mutation_result",
    "format_mutation",
    "build_mutation_report",
    "mutation_report",
    "mutation_signature",
    "mutation_map",
    "result_to_json",

    # UI

    # Tests
    "run_mutation_self_test",
    "assert_mutation_self_test",
]


# ============================================================================
# STANDALONE EXECUTION
# ============================================================================

if __name__ == "__main__":
    from ui import style  # only when run as a script

    results = run_mutation_self_test()

    passed = sum(
        results.values()
    )

    style.header('Mutation analysis self-test')

    for name, success in (
        results.items()
    ):

        print(
            f"[{'PASS' if success else 'FAIL'}] "
            f"{name}"
        )

    print(
        minor_rule()
    )

    print(
        f"{passed}/{len(results)} tests passed."
    )

    if passed != len(results):

        raise SystemExit(1)
