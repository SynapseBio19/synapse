"""Search over enzyme combinations and the automatic cloning planner."""

from __future__ import annotations

from data.restriction_enzymes import RESTRICTION_ENZYMES
from itertools import combinations

from analyses._cloning.cache import _strategy_search_cache
from analyses._cloning.sites import analyze_enzyme, _normalize_topology, _validate_sequence
from analyses._cloning.strategy import evaluate_cloning_strategy


# ============================================================
# AUTOMATIC STRATEGY SEARCH
# ============================================================

def _candidate_enzyme_names(
    enzyme_names=None,
):
    if enzyme_names is None:
        return list(
            RESTRICTION_ENZYMES.keys()
        )

    return [
        name
        for name in dict.fromkeys(
            enzyme_names
        )
        if name in RESTRICTION_ENZYMES
    ]


def _candidate_key(
    result,
):
    strategy = result.get(
        "strategy",
        {},
    )
    validity = result.get(
        "validity",
        {},
    )
    score = result.get(
        "score",
        {},
    )

    return (
        bool(
            validity.get(
                "valid",
                False,
            )
        ),
        bool(
            strategy.get(
                "directional",
                False,
            )
        ),
        bool(
            validity.get(
                "junctions_compatible",
                False,
            )
        ),
        not bool(
            validity.get(
                "reverse_orientation_compatible",
                True,
            )
        ),
        not bool(
            strategy.get(
                "self_ligation_risk",
                True,
            )
        ),
        bool(
            validity.get(
                "no_internal_sites",
                False,
            )
        ),
        int(
            score.get(
                "score",
                0,
            )
        ),
        str(
            result.get(
                "vector",
                {}
            ).get(
                "left",
                {}
            ).get(
                "enzyme",
                "",
            )
        ),
        str(
            result.get(
                "vector",
                {}
            ).get(
                "right",
                {}
            ).get(
                "enzyme",
                "",
            )
        ),
        str(
            result.get(
                "insert",
                {}
            ).get(
                "left",
                {}
            ).get(
                "enzyme",
                "",
            )
        ),
        str(
            result.get(
                "insert",
                {}
            ).get(
                "right",
                {}
            ).get(
                "enzyme",
                "",
            )
        ),
    )


def _strategy_identity(
    result,
):
    vector = result.get(
        "vector",
        {},
    )
    insert = result.get(
        "insert",
        {},
    )

    return (
        vector.get(
            "left",
            {}
        ).get(
            "enzyme"
        ),
        vector.get(
            "right",
            {}
        ).get(
            "enzyme"
        ),
        insert.get(
            "left",
            {}
        ).get(
            "enzyme"
        ),
        insert.get(
            "right",
            {}
        ).get(
            "enzyme"
        ),
        vector.get(
            "left",
            {}
        ).get(
            "position"
        ),
        vector.get(
            "right",
            {}
        ).get(
            "position"
        ),
        insert.get(
            "left",
            {}
        ).get(
            "position"
        ),
        insert.get(
            "right",
            {}
        ).get(
            "position"
        ),
    )


@_strategy_search_cache()
def find_best_cloning_strategies(
    vector_sequence,
    insert_sequence,
    enzyme_names=None,
    top_n=10,
    vector_topology="linear",
    insert_topology="linear",
):
    """
    Exhaustively evaluate physical two-end architectures.

    Crucially, the complete candidate set is generated and ranked
    BEFORE top_n is applied.
    """
    names = _candidate_enzyme_names(
        enzyme_names
    )

    if len(names) < 1:
        return []

    all_candidates = []

    # Different enzymes on both ends of both molecules.
    for vector_left, vector_right in combinations(
        names,
        2,
    ):
        for insert_left, insert_right in combinations(
            names,
            2,
        ):
            result = evaluate_cloning_strategy(
                vector_sequence,
                insert_sequence,
                vector_left,
                vector_right,
                insert_left,
                insert_right,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )

            if result.get(
                "status"
            ) in {
                "VALID",
                "INVALID",
            }:
                all_candidates.append(
                    result
                )

    # Same enzyme at both vector ends and both insert ends.
    for enzyme_name in names:
        result = evaluate_cloning_strategy(
            vector_sequence,
            insert_sequence,
            enzyme_name,
            enzyme_name,
            enzyme_name,
            enzyme_name,
            vector_topology=vector_topology,
            insert_topology=insert_topology,
        )

        if result.get(
            "status"
        ) in {
            "VALID",
            "INVALID",
        }:
            all_candidates.append(
                result
            )

    # Vector: two enzymes; insert: same enzyme.
    for vector_left, vector_right in combinations(
        names,
        2,
    ):
        for insert_enzyme in names:
            result = evaluate_cloning_strategy(
                vector_sequence,
                insert_sequence,
                vector_left,
                vector_right,
                insert_enzyme,
                insert_enzyme,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )

            if result.get(
                "status"
            ) in {
                "VALID",
                "INVALID",
            }:
                all_candidates.append(
                    result
                )

    # Vector: same enzyme; insert: two enzymes.
    for vector_enzyme in names:
        for insert_left, insert_right in combinations(
            names,
            2,
        ):
            result = evaluate_cloning_strategy(
                vector_sequence,
                insert_sequence,
                vector_enzyme,
                vector_enzyme,
                insert_left,
                insert_right,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )

            if result.get(
                "status"
            ) in {
                "VALID",
                "INVALID",
            }:
                all_candidates.append(
                    result
                )

    # Deduplicate by physical architecture.
    unique = {}

    for candidate in all_candidates:
        unique[
            _strategy_identity(
                candidate
            )
        ] = candidate

    candidates = list(
        unique.values()
    )

    candidates.sort(
        key=_candidate_key,
        reverse=True,
    )

    return candidates[
        :max(
            0,
            int(top_n),
        )
    ]


# ============================================================
# COMPLETE SEARCH INTERNAL
# ============================================================

@_strategy_search_cache()
def _find_all_cloning_strategies(
    vector_sequence,
    insert_sequence,
    enzyme_names=None,
    vector_topology="linear",
    insert_topology="linear",
):
    names = _candidate_enzyme_names(
        enzyme_names
    )

    if not names:
        return []

    all_candidates = []

    for vector_left, vector_right in combinations(
        names,
        2,
    ):
        for insert_left, insert_right in combinations(
            names,
            2,
        ):
            all_candidates.append(
                evaluate_cloning_strategy(
                    vector_sequence,
                    insert_sequence,
                    vector_left,
                    vector_right,
                    insert_left,
                    insert_right,
                    vector_topology=vector_topology,
                    insert_topology=insert_topology,
                )
            )

    for enzyme_name in names:
        all_candidates.append(
            evaluate_cloning_strategy(
                vector_sequence,
                insert_sequence,
                enzyme_name,
                enzyme_name,
                enzyme_name,
                enzyme_name,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )
        )

    for vector_left, vector_right in combinations(
        names,
        2,
    ):
        for insert_enzyme in names:
            all_candidates.append(
                evaluate_cloning_strategy(
                    vector_sequence,
                    insert_sequence,
                    vector_left,
                    vector_right,
                    insert_enzyme,
                    insert_enzyme,
                    vector_topology=vector_topology,
                    insert_topology=insert_topology,
                )
            )

    for vector_enzyme in names:
        for insert_left, insert_right in combinations(
            names,
            2,
        ):
            all_candidates.append(
                evaluate_cloning_strategy(
                    vector_sequence,
                    insert_sequence,
                    vector_enzyme,
                    vector_enzyme,
                    insert_left,
                    insert_right,
                    vector_topology=vector_topology,
                    insert_topology=insert_topology,
                )
            )

    unique = {}

    for candidate in all_candidates:
        if candidate.get(
            "status"
        ) not in {
            "VALID",
            "INVALID",
        }:
            continue

        unique[
            _strategy_identity(
                candidate
            )
        ] = candidate

    candidates = list(
        unique.values()
    )

    candidates.sort(
        key=_candidate_key,
        reverse=True,
    )

    return candidates


# ============================================================
# AUTOMATIC PLANNER
# ============================================================

def _site_diagnostics(
    sequence,
    names,
    role,
    topology="linear",
):
    rows = []

    for name in names:
        analysis = analyze_enzyme(
            sequence,
            name,
            topology=topology,
        )

        count = analysis[
            "site_count"
        ]

        if analysis["possible_site_count"]:
            status = "AMBIGUITY_DEPENDENT"
        elif count == 0:
            status = "NO_SITE"
        elif count == 1:
            status = "UNIQUE"
        elif count == 2:
            status = "TWO_SITES"
        else:
            status = "MULTIPLE"

        rows.append({
            "role": role,
            "enzyme": name,
            "site_count": count,
            "sites": list(
                analysis[
                    "sites"
                ]
            ),
            "possible_site_count": analysis[
                "possible_site_count"
            ],
            "possible_sites": list(
                analysis[
                    "possible_sites"
                ]
            ),
            "status": status,
        })

    return rows


def _no_strategy_diagnostics(
    vector_sequence,
    insert_sequence,
    names,
    vector_topology="linear",
    insert_topology="linear",
):
    vector_rows = _site_diagnostics(
        vector_sequence,
        names,
        "vector",
        topology=vector_topology,
    )
    insert_rows = _site_diagnostics(
        insert_sequence,
        names,
        "insert",
        topology=insert_topology,
    )

    # Rows with ambiguity-dependent sites cannot define a clean end.
    vector_end_possible = any(
        row["status"] in {"UNIQUE", "TWO_SITES"}
        for row in vector_rows
    )

    insert_end_possible = any(
        row["status"] in {"UNIQUE", "TWO_SITES"}
        for row in insert_rows
    )

    if not vector_end_possible:
        primary_reason = (
            "The vector contains no enzyme with a usable "
            "unique or two-site restriction architecture."
        )
    elif not insert_end_possible:
        primary_reason = (
            "The insert contains no enzyme with a usable "
            "unique or two-site restriction architecture."
        )
    else:
        primary_reason = (
            "Restriction sites exist, but no complete two-end "
            "architecture passed the cloning checks."
        )

    recommendations = []

    if not vector_end_possible:
        recommendations.append(
            "Choose enzymes that cut the vector once, or exactly twice "
            "when the same enzyme is intentionally used at both ends."
        )

    if not insert_end_possible:
        recommendations.append(
            "Choose enzymes that cut the insert once, or exactly twice "
            "when the same enzyme is intentionally used at both ends."
        )

    if vector_end_possible and insert_end_possible:
        recommendations.append(
            "Inspect the ranked rejected architectures for incompatible "
            "junctions, internal sites, or non-directional orientation."
        )

    return {
        "primary_reason": primary_reason,
        "recommendations": recommendations,
        "vector": {
            "enzyme_diagnostics": vector_rows,
        },
        "insert": {
            "enzyme_diagnostics": insert_rows,
        },
    }


def automatic_cloning_planner(
    vector_sequence,
    insert_sequence,
    enzyme_names=None,
    top_n=5,
    vector_topology="linear",
    insert_topology="linear",
):
    """
    Canonical automatic planner.

    The complete search is always performed before top_n is applied.
    Therefore evaluated_count and valid_count describe the full search,
    not merely the displayed rows. vector_topology / insert_topology are
    passed to every evaluated strategy (see evaluate_cloning_strategy).
    """
    topology = {
        "vector": _normalize_topology(vector_topology),
        "insert": _normalize_topology(insert_topology),
    }

    vector_validation = _validate_sequence(
        vector_sequence,
        "vector",
    )
    insert_validation = _validate_sequence(
        insert_sequence,
        "insert",
    )

    if not vector_validation[
        "valid"
    ]:
        return {
            "status": "INVALID_VECTOR",
            "best": None,
            "strategies": [],
            "valid_strategies": [],
            "evaluated_count": 0,
            "valid_count": 0,
            "rejected_count": 0,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
            "topology": topology,
            "diagnostics": {
                "primary_reason": (
                    "The vector sequence is empty or contains "
                    "unsupported symbols."
                ),
                "recommendations": [
                    "Provide a valid DNA sequence."
                ],
            },
        }

    if not insert_validation[
        "valid"
    ]:
        return {
            "status": "INVALID_INSERT",
            "best": None,
            "strategies": [],
            "valid_strategies": [],
            "evaluated_count": 0,
            "valid_count": 0,
            "rejected_count": 0,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
            "topology": topology,
            "diagnostics": {
                "primary_reason": (
                    "The insert sequence is empty or contains "
                    "unsupported symbols."
                ),
                "recommendations": [
                    "Provide a valid DNA sequence."
                ],
            },
        }

    names = _candidate_enzyme_names(
        enzyme_names
    )

    if not names:
        return {
            "status": "NO_ENZYMES_AVAILABLE",
            "best": None,
            "strategies": [],
            "valid_strategies": [],
            "evaluated_count": 0,
            "valid_count": 0,
            "rejected_count": 0,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
            "topology": topology,
            "diagnostics": {
                "primary_reason": (
                    "No recognized restriction enzymes were supplied."
                ),
                "recommendations": [
                    "Provide at least one enzyme from the Synapse enzyme database."
                ],
            },
        }

    all_candidates = _find_all_cloning_strategies(
        vector_validation[
            "sequence"
        ],
        insert_validation[
            "sequence"
        ],
        names,
        vector_topology=topology["vector"],
        insert_topology=topology["insert"],
    )

    valid_candidates = [
        candidate
        for candidate in all_candidates
        if candidate.get(
            "valid",
            False,
        )
    ]

    rejected_count = (
        len(all_candidates)
        - len(valid_candidates)
    )

    try:
        limit = max(
            0,
            int(top_n),
        )
    except (
        TypeError,
        ValueError,
    ):
        limit = 5

    if valid_candidates:
        best = valid_candidates[0]
        status = "VALID STRATEGY FOUND"
        diagnostics = {
            "primary_reason": (
                "At least one strategy passed all implemented "
                "sequence-level cloning checks."
            ),
            "recommendations": [
                "Review the ranked alternatives before choosing the final design."
            ],
        }
    elif all_candidates:
        best = all_candidates[0]
        status = "NO FULLY VALID STRATEGY FOUND"
        diagnostics = {
            "primary_reason": (
                "Candidate architectures were evaluated, but none "
                "passed every hard validity check."
            ),
            "recommendations": [
                "Review the best rejected strategy and its warnings.",
                "Try alternative enzyme pairs or a different cloning architecture."
            ],
        }
    else:
        best = None
        status = "NO STRATEGY COULD BE EVALUATED"
        diagnostics = _no_strategy_diagnostics(
            vector_validation[
                "sequence"
            ],
            insert_validation[
                "sequence"
            ],
            names,
            vector_topology=topology["vector"],
            insert_topology=topology["insert"],
        )

    return {
        "status": status,
        "best": best,
        "strategies": all_candidates[
            :limit
        ],
        "valid_strategies": valid_candidates[
            :limit
        ],
        "evaluated_count": len(
            all_candidates
        ),
        "valid_count": len(
            valid_candidates
        ),
        "rejected_count": rejected_count,
        "vector_validation": vector_validation,
        "insert_validation": insert_validation,
        "topology": topology,
        "search": {
            "enzyme_count": len(names),
            "enzyme_names": names,
        },
        "diagnostics": diagnostics,
    }
