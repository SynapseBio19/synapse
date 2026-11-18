"""Split plan for analyses/mutation.py (Phase 8, package C)."""

SOURCE = "analyses/mutation.py"
PACKAGE = "analyses._mutation"
PACKAGE_DIR = "analyses/_mutation"
PACKAGE_DOC = """Implementation of :mod:`analyses.mutation` (split out in Phase 8).

Import from ``analyses.mutation``; it re-exports every name defined here,
keeps ``__all__`` and is the only public interface. The modules, in
dependency order:

- ``sequence``: alphabets, IUPAC handling, normalization, translation
- ``models``: the result dataclasses and errors
- ``alignment``: global (banded) alignment, identity, edit distance, blocks
- ``detection``: mutation calls from an alignment, notation, left alignment
- ``coding``: codon and protein consequences, the coding-frame check
- ``hgvs``: HGVS c. and p. descriptions
- ``comparison``: comparing sequences and sets of variants, statistics
- ``context``: regions, sequence quality, context, spectrum, hotspots
- ``impact``: restriction-site, primer and guide overlap of mutations
- ``reporting``: validation, text and JSON reports, the advanced analysis
- ``selftest``: the built-in self-test
"""

MODULES = {
    "sequence": ("Alphabets, IUPAC handling, normalization and translation.", [
        "DNA_ALPHABET", "IUPAC_CODES", "IUPAC_FROM_BASES", "IUPAC_COMPLEMENT", "CODON_TABLE", "AA_NAMES",
        "PURINES", "PYRIMIDINES", "normalize_dna", "normalize_iupac_sequence", "validate_dna_sequence",
        "gc_content", "reverse_complement", "reverse_complement_iupac", "iupac_matches", "iupac_code",
        "translate_dna", "amino_acid_name",
    ]),
    "models": ("Result dataclasses and errors of the mutation analysis.", [
        "SequencePromptCancelled", "MutationAnalysisError", "AlignmentBlock", "Mutation", "MutationStatistics",
        "CodingAnalysis", "MutationAnalysisResult", "SequenceCoordinate", "SequenceQuality", "MutationContext",
        "MutationSpectrum", "MutationHotspot", "RestrictionSiteImpact", "MutationImpact",
        "AdvancedMutationAnalysis",
    ]),
    "alignment": ("Needleman-Wunsch global alignment (banded when exact), identity, edit\n"
                  "distance and alignment blocks.", [
        "_nucleotide_match", "BANDED_ALIGNMENT_MIN_CELLS", "_banded_align", "_banded_is_exact", "global_align",
        "alignment_identity", "_full_edit_distance", "_banded_edit_distance", "_edit_distance",
        "sequence_similarity", "build_alignment_blocks",
    ]),
    "detection": ("Mutation calls from an alignment, their classification and notation, and\n"
                  "left normalization of indels.", [
        "is_transition", "classify_substitution", "_mutation_notation", "_make_mutation",
        "detect_mutations_from_alignment", "_detect_simple_indel", "normalize_indel_left",
    ]),
    "coding": ("Codon and protein consequences of mutations and the coding-frame check.", [
        "classify_coding_substitution", "_translate_canonical", "_first_difference_position",
        "analyze_coding_consequences", "annotate_mutations_with_coding_consequences",
        "CODING_FRAME_ASSUMPTION", "assess_coding_frame", "orf_reference_cds", "compare_proteins",
    ]),
    "hgvs": ("HGVS c. and p. descriptions (3' rule, dup, delins, alleles).", [
        "mutation_hgvs_c", "_reference_span", "hgvs_c_allele", "mutation_hgvs_p",
    ]),
    "comparison": ("Comparing a reference with one or several variants: alignment, mutation\n"
                   "calls, coding consequences and statistics in one result.", [
        "detect_mutations", "analyze_mutations",
        "calculate_mutation_statistics", "MAX_RECOMMENDED_ALIGNMENT_LENGTH", "compare_sequences",
        "compare_multiple_variants", "_mutation_key", "shared_mutations", "unique_mutations",
    ]),
    "context": ("Region annotation, sequence quality, mutation context, spectrum and hotspots.", [
        "annotate_regions", "annotate_sequence_regions", "sequence_complexity", "sequence_quality",
        "extract_mutation_context", "build_mutation_spectrum", "find_mutation_hotspots",
    ]),
    "impact": ("Overlap of mutations with restriction sites, primers and guide positions.", [
        "_find_all_occurrences", "restriction_site_impact", "_interval_overlap", "primer_mutation_overlap",
        "crispr_mutation_overlap",
    ]),
    "reporting": ("Validation of a result, the text and JSON reports, the advanced analysis\n"
                  "and the mutation map.", [
        "validate_mutation_result", "format_mutation", "build_mutation_report", "mutation_signature",
        "create_advanced_analysis", "advanced_mutation_report", "mutation_map", "result_to_json",
        "mutation_report",
    ]),
    "selftest": ("The built-in self-test of the mutation module (menu entry and test suite).", [
        "run_mutation_self_test", "assert_mutation_self_test",
    ]),
}

FACADE_NAMES = {"__all__"}
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 8: the implementation moved to the internal package\n"
    "# ``analyses._mutation`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface: it re-exports every name it\n"
    "# defined before, private helpers included, and keeps ``__all__``.\n"
)
