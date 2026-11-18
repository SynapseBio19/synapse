
# Phase 8: the implementation moved to the internal package
# ``analyses._cloning`` (one module per responsibility, see its __init__).
# This module remains the public interface and re-exports every name it
# defined before, private helpers included. ``_STRATEGY_SEARCH_CACHE``
# here is only the value at import time; the live cache state is
# ``analyses._cloning.cache._STRATEGY_SEARCH_CACHE``.

from __future__ import annotations

# Names analyses.cloning imported and callers use as its attributes (the
# golden master and tests read cloning.RESTRICTION_ENZYMES and check that
# cloning.bases_match is the restriction module's matcher).
from data.restriction_enzymes import RESTRICTION_ENZYMES
from analyses.restriction import bases_match, find_recognition_site_matches
from core.textformat import major_rule, minor_rule
from analyses._cloning.cache import (
    _cached,
    _copy_result,
    _IMMUTABLE,
    _STRATEGY_SEARCH_CACHE,
    _strategy_search_cache,
    _T,
)
from analyses._cloning.sites import (
    _analyze_enzyme,
    analyze_enzyme,
    build_restriction_map,
    find_restriction_sites,
    get_cut_type,
    get_enzyme_data,
    get_overhang,
    get_overhang_polarity,
    get_recognition_sequence,
    is_blunt_enzyme,
    IUPAC_BASES,
    normalize_sequence,
    _normalize_topology,
    reverse_complement,
    validate_dna_sequence,
    _validate_sequence,
    _validate_sequence_uncached,
)
from analyses._cloning.ends import (
    _ambiguity_reason,
    analyze_end_compatibility,
    analyze_junction,
    _end_definition,
    _end_signature,
    ends_compatible,
    make_missing_site,
    make_site,
    _site_polarity,
)
from analyses._cloning.strategy import (
    analyze_cloning,
    analyze_directional_strategy,
    analyze_enzyme_pair,
    assess_enzyme_suitability,
    _digest_fragments,
    evaluate_cloning_strategy,
    find_compatible_enzymes,
    _internal_site_check,
    score_cloning_strategy,
    _score_strategy,
    _selected_digest,
)
from analyses._cloning.search import (
    automatic_cloning_planner,
    _candidate_enzyme_names,
    _candidate_key,
    _find_all_cloning_strategies,
    find_best_cloning_strategies,
    _no_strategy_diagnostics,
    _site_diagnostics,
    _strategy_identity,
)
from analyses._cloning.summary import (
    _cloning_ui_normalize_choice,
    compare_cloning_strategies,
    get_best_cloning_strategy,
    get_cloning_junctions,
    get_cloning_summary,
    get_design_reasons,
    rank_cloning_strategy,
    run_advanced_cloning_analysis,
    run_cloning_analysis,
)
from analyses._cloning.selftest import run_cloning_self_test


# ============================================================
# MODULE SELF TEST
# ============================================================

if __name__ == "__main__":
    from ui import style  # only when run as a script
    tests = run_cloning_self_test()

    style.header('Cloning engine self-test')

    for name, passed in tests.items():
        print(
            f"{name:<48}"
            f"{'PASS' if passed else 'FAIL'}"
        )

    print(
        minor_rule()
    )

    print(
        "RESULT: "
        + (
            "ALL CLONING ENGINE TESTS PASSED"
            if tests[
                "all_passed"
            ]
            else
            "ONE OR MORE CLONING ENGINE TESTS FAILED"
        )
    )

    print(
        major_rule()
    )


# ============================================================
# SELF-TEST ALIAS
# ============================================================

# run_cloning_ui_self_test and the report printers are in ui/cloning_menu.py
# (P2-2).
