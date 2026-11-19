"""Split plan for analyses/orf.py (Phase 9)."""

SOURCE = "analyses/orf.py"
PACKAGE = "analyses._orf"
PACKAGE_DIR = "analyses/_orf"
PACKAGE_DOC = """Implementation of :mod:`analyses.orf` (split out in Phase 9).

Import from ``analyses.orf``; it re-exports every name defined here and is
the only public interface. The modules, in dependency order:

- ``codes``: start/stop codons, policies and defaults, genetic-code tables,
  IUPAC codon handling, start-codon classification
- ``features``: per-ORF evaluation (ambiguity, codon usage, Kozak and
  Shine-Dalgarno context, CDS quality, heuristic score)
- ``finder``: ``find_orfs`` with partial and origin-spanning ORFs
- ``summary``: ``analyze_orfs``, density, frame statistics, relationships,
  ranking and comparison
- ``selection``: translation and analysis of a selected ORF,
  ``analyze_orfs_enhanced``
"""

MODULES = {
    "codes": ("Start and stop codons, search policies, genetic-code tables and codon classification.", [
        "STANDARD_START_CODONS", "ALTERNATIVE_START_CODONS", "ALL_START_CODONS", "STANDARD_STOP_CODONS",
        "DEFAULT_MIN_ORF_LENGTH", "START_CODON_POLICIES", "DEFAULT_START_CODON_POLICY", "DEFAULT_REPORT_NESTED",
        "DEFAULT_RANKING", "ORF_RANKING_MODES", "LEGACY_ORF_SETTINGS", "DNA_IUPAC", "GENETIC_CODE",
        "expand_iupac_codon", "_stop_codons_of", "codon_table_for", "start_codons_for", "initiation_codons",
        "is_initiation_codon", "_genetic_code_name", "analyze_codon_identity", "_is_alternative_start",
        "classify_start_codon", "resolve_start_codons",
    ]),
    "features": ("Per-ORF evaluation: ambiguity, codon usage, start context, CDS quality, heuristic score.", [
        "analyze_orf_ambiguity", "analyze_codon_usage", "calculate_codon_diversity", "analyze_kozak_context",
        "analyze_shine_dalgarno", "analyze_start_context", "calculate_cds_quality",
        "calculate_orf_heuristic_score", "calculate_orf_confidence",
    ]),
    "finder": ("The ORF search: complete, 5'-partial and origin-spanning ORFs.", [
        "_orf_sort_key", "_keep_longest_per_stop", "MISSING_START_STATUSES", "orf_has_start_codon",
        "_five_prime_partial_orfs", "_origin_spanning_orfs", "orf_end_on_sequence", "find_orfs",
    ]),
    "summary": ("analyze_orfs, density, frame statistics, relationships, ranking and comparison.", [
        "calculate_orf_density", "calculate_frame_statistics", "summarize_orfs", "analyze_orfs",
        "create_orf_preview", "_orf_genomic_interval", "analyze_orf_relationships", "rank_orf_candidates",
        "compare_orf_candidates",
    ]),
    "selection": ("Translation and analysis of a selected ORF, and the enhanced analysis.", [
        "translate_orf_sequence", "analyze_selected_orf", "calculate_orf_feature_summary", "analyze_orfs_enhanced",
    ]),
}

FACADE_NAMES = set()
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 9: the implementation moved to the internal package\n"
    "# ``analyses._orf`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface and re-exports every name it\n"
    "# defined before.\n"
)
