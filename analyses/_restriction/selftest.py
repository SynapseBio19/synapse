"""The built-in self-test of the restriction analysis (menu entry and test suite)."""

from __future__ import annotations

from data.restriction_enzymes import get_restriction_enzyme, RESTRICTION_ENZYMES

from analyses._restriction.sites import (
    bases_match,
    find_restriction_site_matches,
    reverse_complement,
    validate_restriction_sequence,
)
from analyses._restriction.digest import (
    calculate_circular_fragments,
    calculate_gel_bands,
    calculate_linear_fragments,
)
from analyses._restriction.compatibility import sticky_ends_are_compatible
from analyses._restriction.api import validate_enzyme_database


def run_restriction_self_test():
    """Deterministic structural tests for the final restriction engine."""
    results = {}

    results["validation"] = (
        validate_restriction_sequence("ATGC")
        and not validate_restriction_sequence("ATGX")
        and validate_restriction_sequence("ATGN")
        and not validate_restriction_sequence("")
    )

    results["iupac"] = (
        bases_match("A", "R")
        and bases_match("G", "R")
        and not bases_match("C", "R")
    )

    results["reverse_complement"] = (
        reverse_complement("ATGCN") == "NGCAT"
        and reverse_complement("ARYN") == "NRYT"
    )

    results["linear_fragments"] = (
        calculate_linear_fragments(
            100,
            [25, 75],
        )
        == [
            {
                "fragment": 1,
                "start": 1,
                "end": 25,
                "start_0based": 0,
                "end_0based_exclusive": 25,
                "length": 25,
                "wraps_origin": False,
            },
            {
                "fragment": 2,
                "start": 26,
                "end": 75,
                "start_0based": 25,
                "end_0based_exclusive": 75,
                "length": 50,
                "wraps_origin": False,
            },
            {
                "fragment": 3,
                "start": 76,
                "end": 100,
                "start_0based": 75,
                "end_0based_exclusive": 100,
                "length": 25,
                "wraps_origin": False,
            },
        ]
    )

    circular = calculate_circular_fragments(
        100,
        [25, 75],
    )
    results["circular_fragment_total"] = (
        sum(
            fragment["length"]
            for fragment in circular
        ) == 100
        and len(circular) == 2
    )

    single = calculate_circular_fragments(
        100,
        [25],
    )
    results["single_circular_cut"] = (
        len(single) == 1
        and single[0]["length"] == 100
    )

    results["gel_band_collapse"] = (
        len(
            calculate_gel_bands(
                [
                    {"fragment": 1, "length": 100},
                    {"fragment": 2, "length": 100},
                    {"fragment": 3, "length": 50},
                ]
            )
        ) == 2
        and calculate_gel_bands(
            [
                {"fragment": 1, "length": 100},
                {"fragment": 2, "length": 100},
                {"fragment": 3, "length": 50},
            ]
        )[0]["copies"] == 2
    )

    # Same displayed overhang and reverse-complement representation
    # are both accepted as compatible database encodings.
    site_a = {
        "is_blunt": False,
        "overhang_sequence": "AATT",
    }
    site_b = {
        "is_blunt": False,
        "overhang_sequence": "AATT",
    }
    results["sticky_compatibility"] = (
        sticky_ends_are_compatible(site_a, site_b)
        and sticky_ends_are_compatible(
            site_a,
            {
                "is_blunt": False,
                "overhang_sequence": "AATT",
            },
        )
        and not sticky_ends_are_compatible(
            site_a,
            {
                "is_blunt": False,
                "overhang_sequence": "GGCC",
            },
        )
    )

    results["database_validation_api"] = isinstance(
        validate_enzyme_database(),
        dict,
    )

    # Regression test for linear digest coordinate convention.
    # A 75-bp molecule cut at 0-based boundaries 9 and 28
    # produces 9, 19, and 47 bp fragments.
    # This is intentionally NOT 10, 19, 46: the cut coordinate
    # represents the boundary after base 9 / after base 28.
    regression_fragments = calculate_linear_fragments(
        75,
        [9, 28],
    )
    results["linear_digest_coordinate_regression"] = (
        [fragment["length"] for fragment in regression_fragments]
        == [9, 19, 47]
        and sum(
            fragment["length"]
            for fragment in regression_fragments
        ) == 75
        and len(regression_fragments) == 3
    )

    # The public interface is analyses.restriction. Before Phase 9 this check
    # read globals(), which were that module's; the self-test now lives in
    # analyses/_restriction/selftest.py, so it asks the public module.
    import analyses.restriction as public_module

    results["public_api"] = all(
        hasattr(public_module, name)
        for name in (
            "analyze_enzyme",
            "create_restriction_map",
            "perform_digest",
            "calculate_gel_bands",
            "find_unique_cutters",
            "find_non_cutters",
            "find_internal_sites",
            "recommend_enzyme_pairs",
            "analyze_restriction",
            "run_restriction_analysis",
            "get_restriction_summary",
            "get_restriction_sites_for_enzyme",
            "get_restriction_map",
            "get_digest_prediction",
            "get_gel_band_prediction",
            "get_unique_cutters",
            "get_non_cutters",
            "get_enzyme_pair_recommendations",
            "run_double_digest",
            "get_cloning_compatible_sites",
        )
    )

    try:
        first_enzyme = next(
            iter(RESTRICTION_ENZYMES)
        )

        results["enzyme_database"] = isinstance(
            get_restriction_enzyme(first_enzyme),
            dict,
        )

        results["site_match_api"] = isinstance(
            find_restriction_site_matches(
                "ATGC" * 100,
                first_enzyme,
            ),
            list,
        )
    except Exception:
        results["enzyme_database"] = False
        results["site_match_api"] = False

    results["all_passed"] = all(
        results.values()
    )

    return results
