"""Entry points, summaries, ranking and comparison of cloning results."""

from __future__ import annotations

from analyses._cloning.search import automatic_cloning_planner, _candidate_key


# ============================================================
# PIPELINE API
# ============================================================

def run_advanced_cloning_analysis(
    vector_sequence,
    insert_sequence,
    enzyme_names=None,
    top_n=5,
    vector_topology="linear",
    insert_topology="linear",
):
    planner = automatic_cloning_planner(
        vector_sequence,
        insert_sequence,
        enzyme_names=enzyme_names,
        top_n=top_n,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )

    valid = (
        planner["best"] is not None
        and planner["best"].get(
            "valid",
            False,
        )
    )

    return {
        "status": planner["status"],
        "valid": valid,
        "vector_validation": planner[
            "vector_validation"
        ],
        "insert_validation": planner[
            "insert_validation"
        ],
        "planner": planner,
    }


def run_cloning_analysis(
    vector_sequence,
    insert_sequence,
    mode="auto",
    enzyme_names=None,
    top_n=5,
    vector_topology="linear",
    insert_topology="linear",
):
    result = run_advanced_cloning_analysis(
        vector_sequence,
        insert_sequence,
        enzyme_names=enzyme_names,
        top_n=top_n,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )

    result["mode"] = mode
    return result


def get_best_cloning_strategy(
    vector_sequence,
    insert_sequence,
    enzyme_names=None,
    vector_topology="linear",
    insert_topology="linear",
):
    return automatic_cloning_planner(
        vector_sequence,
        insert_sequence,
        enzyme_names=enzyme_names,
        top_n=1,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    ).get(
        "best"
    )


def get_cloning_junctions(
    analysis,
):
    if not isinstance(
        analysis,
        dict,
    ):
        return {}

    if "junctions" in analysis:
        junctions = (
            analysis.get(
                "junctions"
            )
            or {}
        )
    else:
        planner = analysis.get(
            "planner",
            {},
        )
        best = planner.get(
            "best"
        )

        if not isinstance(
            best,
            dict,
        ):
            return {}

        junctions = (
            best.get(
                "junctions"
            )
            or {}
        )

    return {
        "left": dict(
            junctions.get(
                "left",
                {},
            )
        ),
        "right": dict(
            junctions.get(
                "right",
                {},
            )
        ),
    }


def get_cloning_summary(
    analysis,
):
    if not analysis:
        return {
            "status": "NO_ANALYSIS",
            "valid": False,
        }

    planner = analysis.get(
        "planner",
        analysis,
    )

    if not isinstance(
        planner,
        dict,
    ):
        return {
            "status": "UNKNOWN",
            "valid": False,
        }

    best = planner.get(
        "best"
    )

    if best is None:
        return {
            "status": planner.get(
                "status",
                "UNKNOWN",
            ),
            "valid": False,
            "evaluated_strategies": planner.get(
                "evaluated_count",
                0,
            ),
            "valid_strategies": planner.get(
                "valid_count",
                0,
            ),
            "diagnostics": planner.get(
                "diagnostics",
                {},
            ),
        }

    vector = best[
        "vector"
    ]
    insert = best[
        "insert"
    ]
    strategy = best[
        "strategy"
    ]
    score = best[
        "score"
    ]

    return {
        "status": planner[
            "status"
        ],
        "valid": best[
            "validity"
        ][
            "valid"
        ],
        "score": score[
            "score"
        ],
        "score_maximum": score[
            "maximum"
        ],
        "category": score[
            "category"
        ],
        "directional": strategy[
            "directional"
        ],
        "orientation": strategy[
            "orientation"
        ],
        "vector_left_enzyme": vector[
            "left"
        ][
            "enzyme"
        ],
        "vector_right_enzyme": vector[
            "right"
        ][
            "enzyme"
        ],
        "insert_left_enzyme": insert[
            "left"
        ][
            "enzyme"
        ],
        "insert_right_enzyme": insert[
            "right"
        ][
            "enzyme"
        ],
        "evaluated_strategies": planner[
            "evaluated_count"
        ],
        "valid_strategies": planner[
            "valid_count"
        ],
        "warnings": best.get(
            "warnings",
            [],
        ),
        "digest": best.get(
            "digest",
            {},
        ),
    }


def rank_cloning_strategy(
    result,
):
    return _candidate_key(
        result
    )


# ============================================================
# STRATEGY COMPARISON
# ============================================================

def compare_cloning_strategies(
    first,
    second,
):
    first_score = first.get(
        "score",
        {},
    ).get(
        "score",
        0,
    )

    second_score = second.get(
        "score",
        {},
    ).get(
        "score",
        0,
    )

    return {
        "score_difference": (
            first_score
            - second_score
        ),
        "first_better": (
            first_score
            > second_score
        ),
        "second_better": (
            second_score
            > first_score
        ),
        "same_score": (
            first_score
            == second_score
        ),
        "directionality_difference": (
            first.get(
                "strategy",
                {},
            ).get(
                "directional",
                False,
            )
            !=
            second.get(
                "strategy",
                {},
            ).get(
                "directional",
                False,
            )
        ),
        "self_ligation_difference": (
            first.get(
                "strategy",
                {},
            ).get(
                "self_ligation_risk",
                True,
            )
            !=
            second.get(
                "strategy",
                {},
            ).get(
                "self_ligation_risk",
                True,
            )
        ),
    }


def get_design_reasons(
    analysis,
):
    if not analysis:
        return {
            "positive": [],
            "warnings": [],
        }

    score = analysis.get(
        "score",
        {},
    )

    reasons = list(
        score.get(
            "reasons",
            [],
        )
    )

    warnings = list(
        analysis.get(
            "warnings",
            [],
        )
    )

    return {
        "positive": reasons,
        "warnings": warnings,
    }

# ============================================================
# USER INTERFACE
# ============================================================
#
# Synapse module UI.
# The cloning engine above is intentionally unchanged.
# This section contains only the CLI interaction layer.
# ============================================================

def _cloning_ui_normalize_choice(value):
    return str(value).strip()
