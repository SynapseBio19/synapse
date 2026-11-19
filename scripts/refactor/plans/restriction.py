"""Split plan for analyses/restriction.py (Phase 9)."""

SOURCE = "analyses/restriction.py"
PACKAGE = "analyses._restriction"
PACKAGE_DIR = "analyses/_restriction"
PACKAGE_DOC = """Implementation of :mod:`analyses.restriction` (split out in Phase 9).

Import from ``analyses.restriction``; it re-exports every name defined here
and is the only public interface. The modules, in dependency order:

- ``sites``: sequence validation, IUPAC site matching (also Type IIS),
  cut positions, overhangs, per-site and per-enzyme analysis
- ``digest``: restriction maps, cut positions, linear and circular
  fragments, digests, gel bands, cutter lists
- ``compatibility``: sticky ends, cloning compatibility, enzyme-pair
  recommendations, double digests
- ``api``: the combined analysis, summaries and the ``get_*`` / ``run_*``
  entry points, the enzyme-database check
- ``selftest``: the built-in self-test
"""

MODULES = {
    "sites": ("Sequence validation, IUPAC recognition-site matching, cut positions and overhangs.", [
        "IUPAC_BASES", "_normalize_restriction_sequence", "_restriction_sequence_status",
        "validate_restriction_sequence", "_require_valid_restriction_sequence", "bases_match",
        "recognition_site_matches", "find_restriction_site_matches", "find_recognition_site_matches",
        "_scan_recognition_sites", "find_type_iis_sites", "find_restriction_sites", "calculate_cut_positions",
        "reverse_complement", "get_overhang_sequence", "_site_match_for_position", "analyze_restriction_site",
        "analyze_enzyme",
    ]),
    "digest": ("Restriction maps, fragments, digests, gel bands and cutter lists.", [
        "create_restriction_map", "get_cut_positions", "calculate_linear_fragments",
        "calculate_circular_fragments", "perform_digest", "calculate_gel_bands", "summarize_restriction_map",
        "find_unique_cutters", "find_non_cutters", "find_internal_sites", "get_physical_cut_boundaries",
    ]),
    "compatibility": ("Sticky-end compatibility, enzyme-pair recommendations and double digests.", [
        "get_sticky_end_sequence", "_site_overhang_polarity", "sticky_ends_are_compatible",
        "analyze_cloning_compatibility", "recommend_enzyme_pairs", "analyze_double_digest",
    ]),
    "api": ("The combined analysis, summaries and the get_*/run_* entry points.", [
        "analyze_restriction", "summarize_enzyme_sites", "validate_enzyme_database",
        "compare_single_enzyme_digest", "run_restriction_analysis", "get_restriction_summary",
        "get_restriction_sites_for_enzyme", "get_restriction_map", "get_digest_prediction",
        "get_gel_band_prediction", "get_unique_cutters", "get_non_cutters", "get_enzyme_pair_recommendations",
        "run_double_digest", "get_cloning_compatible_sites",
    ]),
    "selftest": ("The built-in self-test of the restriction analysis (menu entry and test suite).", [
        "run_restriction_self_test",
    ]),
}

FACADE_NAMES = set()
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 9: the implementation moved to the internal package\n"
    "# ``analyses._restriction`` (one module per responsibility, see its\n"
    "# __init__). This module remains the public interface and re-exports\n"
    "# every name it defined before, and RESTRICTION_ENZYMES, which tests read.\n"
)
