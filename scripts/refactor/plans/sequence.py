"""Split plan for analyses/sequence.py (Phase 9)."""

SOURCE = "analyses/sequence.py"
PACKAGE = "analyses._sequence"
PACKAGE_DIR = "analyses/_sequence"
PACKAGE_DOC = """Implementation of :mod:`analyses.sequence` (split out in Phase 9).

Import from ``analyses.sequence``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``bases``: base tables, normalization, counts, complement
- ``composition``: percentages, GC/AT content and ranges, skews, quality,
  molecular weight
- ``profiles``: sliding-window GC and complexity profiles, homopolymers
- ``cpg``: CpG statistics and CpG islands
- ``summary``: sequence forms, interpretation, the summary and
  ``analyze_sequence``
- ``selftest``: the built-in self-test
"""

MODULES = {
    "bases": ("DNA base tables, normalization, nucleotide counts and complements.", [
        "DNA_STANDARD_BASES", "DNA_AMBIGUOUS_CODES", "DNA_IUPAC_BASES", "DNA_COMPLEMENT",
        "DNA_NUCLEOTIDE_RESIDUE_MASSES", "AVERAGE_DSDNA_BP_MASS", "normalize_dna_sequence", "_is_standard_base",
        "_is_iupac_base", "_get_standard_bases", "_get_valid_iupac_bases", "count_nucleotides",
        "count_ambiguous_bases", "get_complement", "get_reverse_complement",
    ]),
    "composition": ("Composition, GC/AT content, skews, sequence quality and molecular weight.", [
        "calculate_sequence_statistics", "_calculate_gc_for_window", "analyze_sequence_quality",
        "calculate_nucleotide_percentages", "calculate_gc_content", "calculate_at_content",
        "calculate_dna_content_ranges", "calculate_gc_skew", "calculate_at_skew",
        "calculate_dna_molecular_weight", "calculate_average_dna_molecular_weight",
    ]),
    "profiles": ("Sliding-window GC and complexity profiles, their statistics, homopolymers.", [
        "describe_profile_availability", "calculate_gc_content_profile", "calculate_sequence_complexity",
        "calculate_complexity_profile", "calculate_profile_statistics", "analyze_gc_profile",
        "analyze_complexity_profile", "find_homopolymers",
    ]),
    "cpg": ("CpG observed/expected statistics and CpG island search.", [
        "DEFAULT_CPG_ISLAND_CRITERIA", "analyze_cpg", "find_cpg_islands",
    ]),
    "summary": ("Sequence forms, interpretation, the summary and analyze_sequence.", [
        "get_sequence_forms", "create_sequence_interpretation", "create_sequence_summary", "analyze_sequence",
    ]),
    "selftest": ("The built-in self-test of the sequence overview (menu entry and test suite).", [
        "run_sequence_self_test",
    ]),
}

FACADE_NAMES = set()
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 9: the implementation moved to the internal package\n"
    "# ``analyses._sequence`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface and re-exports every name it\n"
    "# defined before.\n"
)
