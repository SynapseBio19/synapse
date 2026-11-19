
# Phase 9: the implementation moved to the internal package
# ``analyses._restriction`` (one module per responsibility, see its
# __init__). This module remains the public interface and re-exports
# every name it defined before, and the enzyme-data names it imported.

from __future__ import annotations

# Imported by analyses.restriction before the split; the restriction menu and
# tests read them here.
from data.restriction_enzymes import RESTRICTION_ENZYMES, get_restriction_enzyme, get_type_iis_enzyme
from core.textformat import minor_rule
from analyses._restriction.sites import (
    analyze_enzyme,
    analyze_restriction_site,
    bases_match,
    calculate_cut_positions,
    find_recognition_site_matches,
    find_restriction_site_matches,
    find_restriction_sites,
    find_type_iis_sites,
    get_overhang_sequence,
    IUPAC_BASES,
    _normalize_restriction_sequence,
    recognition_site_matches,
    _require_valid_restriction_sequence,
    _restriction_sequence_status,
    reverse_complement,
    _scan_recognition_sites,
    _site_match_for_position,
    validate_restriction_sequence,
)
from analyses._restriction.digest import (
    calculate_circular_fragments,
    calculate_gel_bands,
    calculate_linear_fragments,
    create_restriction_map,
    find_internal_sites,
    find_non_cutters,
    find_unique_cutters,
    get_cut_positions,
    get_physical_cut_boundaries,
    perform_digest,
    summarize_restriction_map,
)
from analyses._restriction.compatibility import (
    analyze_cloning_compatibility,
    analyze_double_digest,
    get_sticky_end_sequence,
    recommend_enzyme_pairs,
    _site_overhang_polarity,
    sticky_ends_are_compatible,
)
from analyses._restriction.api import (
    analyze_restriction,
    compare_single_enzyme_digest,
    get_cloning_compatible_sites,
    get_digest_prediction,
    get_enzyme_pair_recommendations,
    get_gel_band_prediction,
    get_non_cutters,
    get_restriction_map,
    get_restriction_sites_for_enzyme,
    get_restriction_summary,
    get_unique_cutters,
    run_double_digest,
    run_restriction_analysis,
    summarize_enzyme_sites,
    validate_enzyme_database,
)
from analyses._restriction.selftest import run_restriction_self_test


# The terminal dialogs of this module are in ui/restriction_menu.py (P2-2).


if __name__ == "__main__":
    from ui import style  # only when run as a script
    self_test = run_restriction_self_test()

    style.header('Restriction engine self-test')

    for name, passed in self_test.items():
        if name == "all_passed":
            continue
        print(f"{name:<40}{'PASS' if passed else 'FAIL'}")

    print(minor_rule())
    print(
        "RESULT: "
        + (
            "ALL TESTS PASSED"
            if self_test["all_passed"]
            else "ONE OR MORE TESTS FAILED"
        )
    )


# ============================================================
# RESTRICTION & CLONING UI
# ============================================================
#
# UI layer only. The restriction engine above is unchanged.
#
# Top level:
#   [1] Restriction Analysis
#   [2] Cloning Analysis
#
# Restriction:
#   [1] Restriction Overview
#   [2] Restriction Map
#   [3] Digest & Gel Analysis
#   [4] Complete Restriction Report
#
# ============================================================
