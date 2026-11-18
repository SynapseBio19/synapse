"""Split plan for analyses/cloning.py (Phase 8, package D)."""

SOURCE = "analyses/cloning.py"
PACKAGE = "analyses._cloning"
PACKAGE_DIR = "analyses/_cloning"
PACKAGE_DOC = """Implementation of :mod:`analyses.cloning` (split out in Phase 8).

Import from ``analyses.cloning``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``cache``: the per-search result cache (``_strategy_search_cache``)
- ``sites``: sequences, enzyme data, overhangs, site finding, restriction maps
- ``ends``: cut ends, their compatibility and the junctions
- ``strategy``: evaluating and scoring one cloning strategy
- ``search``: searching all enzyme combinations, the automatic planner
- ``summary``: entry points, summaries, ranking and comparison
- ``selftest``: the built-in self-test

Tests that patch a helper must patch the module that calls it, for example
``analyses._cloning.sites._analyze_enzyme``; the cache state is
``analyses._cloning.cache._STRATEGY_SEARCH_CACHE``.
"""

MODULES = {
    "cache": ("Result cache for one strategy search (P2-6).", [
        "_STRATEGY_SEARCH_CACHE", "_strategy_search_cache", "_IMMUTABLE", "_T", "_copy_result", "_cached",
    ]),
    "sites": ("Sequences, enzyme data, overhangs, recognition sites and restriction maps.", [
        "IUPAC_BASES", "normalize_sequence", "_normalize_topology", "reverse_complement", "get_enzyme_data",
        "get_recognition_sequence", "get_cut_type", "get_overhang", "get_overhang_polarity", "is_blunt_enzyme",
        "find_restriction_sites", "analyze_enzyme", "_analyze_enzyme", "build_restriction_map",
        "_validate_sequence", "_validate_sequence_uncached", "validate_dna_sequence",
    ]),
    "ends": ("Cut ends, their compatibility and the resulting junctions.", [
        "make_site", "make_missing_site", "_end_signature", "_site_polarity", "ends_compatible",
        "analyze_end_compatibility", "analyze_junction", "_ambiguity_reason", "_end_definition",
    ]),
    "strategy": ("Evaluation and scoring of one cloning strategy.", [
        "_internal_site_check", "analyze_directional_strategy", "_digest_fragments", "_selected_digest",
        "_score_strategy", "score_cloning_strategy", "evaluate_cloning_strategy", "analyze_enzyme_pair",
        "analyze_cloning", "assess_enzyme_suitability", "find_compatible_enzymes",
    ]),
    "search": ("Search over enzyme combinations and the automatic cloning planner.", [
        "_candidate_enzyme_names", "_candidate_key", "_strategy_identity", "find_best_cloning_strategies",
        "_find_all_cloning_strategies", "_site_diagnostics", "_no_strategy_diagnostics",
        "automatic_cloning_planner",
    ]),
    "summary": ("Entry points, summaries, ranking and comparison of cloning results.", [
        "run_advanced_cloning_analysis", "run_cloning_analysis", "get_best_cloning_strategy",
        "get_cloning_junctions", "get_cloning_summary", "rank_cloning_strategy", "compare_cloning_strategies",
        "get_design_reasons", "_cloning_ui_normalize_choice",
    ]),
    "selftest": ("The built-in self-test of the cloning module (menu entry and test suite).", [
        "run_cloning_self_test",
    ]),
}

FACADE_NAMES = set()
REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 8: the implementation moved to the internal package\n"
    "# ``analyses._cloning`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface and re-exports every name it\n"
    "# defined before, private helpers included. ``_STRATEGY_SEARCH_CACHE``\n"
    "# here is only the value at import time; the live cache state is\n"
    "# ``analyses._cloning.cache._STRATEGY_SEARCH_CACHE``.\n"
)
