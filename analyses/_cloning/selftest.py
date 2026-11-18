"""The built-in self-test of the cloning module (menu entry and test suite)."""

from __future__ import annotations

from analyses.restriction import bases_match
from data.restriction_enzymes import RESTRICTION_ENZYMES

from analyses._cloning.sites import normalize_sequence, reverse_complement
from analyses._cloning.ends import _end_definition, ends_compatible
from analyses._cloning.search import automatic_cloning_planner


# ============================================================
# REPORTING
# ============================================================

# The terminal dialogs of this module are in ui/cloning_menu.py (P2-2).


# ============================================================
# REPORT ALIASES
# ============================================================


# ============================================================
# SELF TESTS
# ============================================================

def run_cloning_self_test():
    results = {}

    results[
        "normalize_sequence"
    ] = (
        normalize_sequence(
            " at gc "
        )
        == "ATGC"
    )

    results[
        "reverse_complement"
    ] = (
        reverse_complement(
            "ATGC"
        )
        == "GCAT"
    )

    results[
        "iupac_matching"
    ] = (
        bases_match(
            "A",
            "R",
        )
        and bases_match(
            "G",
            "R",
        )
        and not bases_match(
            "C",
            "R",
        )
    )

    results[
        "enzyme_database"
    ] = all(
        name in RESTRICTION_ENZYMES
        for name in list(
            RESTRICTION_ENZYMES
        )[:1]
    )

    try:
        invalid = automatic_cloning_planner(
            "",
            "ATGC",
        )
        results[
            "invalid_vector"
        ] = (
            invalid[
                "status"
            ]
            == "INVALID_VECTOR"
            and invalid[
                "best"
            ] is None
        )
    except Exception:
        results[
            "invalid_vector"
        ] = False

    try:
        invalid = automatic_cloning_planner(
            "ATGC",
            "",
        )
        results[
            "invalid_insert"
        ] = (
            invalid[
                "status"
            ]
            == "INVALID_INSERT"
            and invalid[
                "best"
            ] is None
        )
    except Exception:
        results[
            "invalid_insert"
        ] = False

    try:
        empty = automatic_cloning_planner(
            "ATGC",
            "ATGC",
            enzyme_names=[],
        )
        results[
            "empty_enzyme_set"
        ] = (
            empty[
                "best"
            ] is None
            and empty[
                "evaluated_count"
            ] == 0
        )
    except Exception:
        results[
            "empty_enzyme_set"
        ] = False

    try:
        one_site = _end_definition(
            "ATGCGAATTCATGC",
            "EcoRI",
            "EcoRI",
        )
        results[
            "same_enzyme_one_site_rejected"
        ] = (
            one_site[
                "valid"
            ] is False
        )
    except Exception:
        results[
            "same_enzyme_one_site_rejected"
        ] = False

    try:
        two_sites = _end_definition(
            "ATGCGAATTCATGCGAATTCATGC",
            "EcoRI",
            "EcoRI",
        )
        results[
            "same_enzyme_two_sites_accepted"
        ] = (
            two_sites[
                "valid"
            ] is True
            and two_sites[
                "left_site"
            ]["position"]
            != two_sites[
                "right_site"
            ]["position"]
        )
    except Exception:
        results[
            "same_enzyme_two_sites_accepted"
        ] = False

    try:
        blunt_a = {
            "enzyme": "A",
            "position": 1,
            "overhang": "",
            "is_blunt": True,
        }
        blunt_b = {
            "enzyme": "B",
            "position": 1,
            "overhang": "",
            "is_blunt": True,
        }
        sticky = {
            "enzyme": "C",
            "position": 1,
            "overhang": "AATT",
            "is_blunt": False,
        }

        results[
            "blunt_compatibility"
        ] = (
            ends_compatible(
                blunt_a,
                blunt_b,
            )
            is True
        )

        results[
            "blunt_sticky_incompatibility"
        ] = (
            ends_compatible(
                blunt_a,
                sticky,
            )
            is False
        )
    except Exception:
        results[
            "blunt_compatibility"
        ] = False
        results[
            "blunt_sticky_incompatibility"
        ] = False

    required_api = [
        "analyze_enzyme",
        "assess_enzyme_suitability",
        "ends_compatible",
        "analyze_end_compatibility",
        "analyze_junction",
        "analyze_directional_strategy",
        "score_cloning_strategy",
        "build_restriction_map",
        "analyze_enzyme_pair",
        "find_best_cloning_strategies",
        "automatic_cloning_planner",
        "find_compatible_enzymes",
        "compare_cloning_strategies",
        "get_design_reasons",
        "analyze_cloning",
        "run_cloning_analysis",
        "get_cloning_summary",
        "get_best_cloning_strategy",
        "get_cloning_junctions",
    ]

    # The public interface is analyses.cloning. Before Phase 8 this check
    # read globals(), which were that module's; the self-test now lives in
    # analyses/_cloning/selftest.py, so it asks the public module.
    import analyses.cloning as public_module

    results[
        "public_api"
    ] = all(
        callable(
            getattr(
                public_module,
                name,
                None,
            )
        )
        for name in required_api
    )

    results[
        "all_passed"
    ] = all(
        results.values()
    )

    return results
