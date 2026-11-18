"""Evaluation and scoring of one cloning strategy."""

from __future__ import annotations

from data.restriction_enzymes import RESTRICTION_ENZYMES

from analyses._cloning.sites import (
    analyze_enzyme,
    get_enzyme_data,
    get_overhang,
    is_blunt_enzyme,
    _normalize_topology,
    _validate_sequence,
)
from analyses._cloning.ends import (
    _ambiguity_reason,
    analyze_end_compatibility,
    analyze_junction,
    _end_definition,
    ends_compatible,
    make_missing_site,
    make_site,
)


def _internal_site_check(
    sequence,
    definition,
    topology="linear",
):
    if not definition["valid"]:
        return {
            "valid": False,
            "additional_sites": [],
            "additional_site_count": 0,
            "reason": definition[
                "reason"
            ],
        }

    expected_positions = {
        definition["left_site"][
            "position"
        ],
        definition["right_site"][
            "position"
        ],
    }

    selected_enzymes = {
        definition["left_site"][
            "enzyme"
        ],
        definition["right_site"][
            "enzyme"
        ],
    }

    additional = []

    for enzyme_name in selected_enzymes:
        analysis = analyze_enzyme(
            sequence,
            enzyme_name,
            topology=topology,
        )

        for position in analysis[
            "sites"
        ]:
            if position not in expected_positions:
                additional.append({
                    "enzyme": enzyme_name,
                    "position": position,
                })

    return {
        "valid": not additional,
        "additional_sites": additional,
        "additional_site_count": len(
            additional
        ),
        "reason": (
            "No additional selected-enzyme "
            "restriction sites detected."
            if not additional
            else
            "Additional selected-enzyme "
            "restriction sites detected."
        ),
    }


# ============================================================
# ORIENTATION
# ============================================================

def analyze_directional_strategy(
    vector_left,
    vector_right,
    insert_left,
    insert_right,
):
    left = analyze_end_compatibility(
        vector_left,
        insert_left,
    )
    right = analyze_end_compatibility(
        vector_right,
        insert_right,
    )

    reverse_left = analyze_end_compatibility(
        vector_left,
        insert_right,
    )
    reverse_right = analyze_end_compatibility(
        vector_right,
        insert_left,
    )

    intended_valid = (
        left["compatible"] is True
        and right["compatible"] is True
    )

    reverse_valid = (
        reverse_left["compatible"] is True
        and reverse_right["compatible"] is True
    )

    directional = (
        intended_valid
        and not reverse_valid
    )

    vector_end_distinction = (
        vector_left["enzyme"]
        != vector_right["enzyme"]
        or (
            vector_left.get(
                "position"
            )
            != vector_right.get(
                "position"
            )
        )
    )

    self_ligation_risk = (
        ends_compatible(
            vector_left,
            vector_right,
        )
        is True
    )

    if directional:
        orientation = "Directional"
        status = (
            "Both intended junctions are compatible "
            "and the reverse orientation is blocked."
        )
    elif intended_valid:
        orientation = "Non-directional"
        status = (
            "Both intended junctions are compatible, "
            "but the reverse orientation is also possible."
        )
    else:
        orientation = "Invalid"
        status = (
            "One or both intended junctions are incompatible."
        )

    return {
        "left": left,
        "right": right,
        "reverse_left": reverse_left,
        "reverse_right": reverse_right,
        "directional": directional,
        "orientation": orientation,
        "status": status,
        "reverse_orientation_compatible": (
            reverse_valid
        ),
        "self_ligation_risk": (
            self_ligation_risk
        ),
        "vector_end_distinction": (
            vector_end_distinction
        ),
    }


# ============================================================
# DIGEST / FRAGMENT ANALYSIS
# ============================================================

def _digest_fragments(
    sequence_length,
    cut_positions,
    topology="linear",
):
    positions = sorted(
        set(cut_positions)
    )

    if topology == "circular":
        if not positions:
            return [{
                "fragment": 1,
                "length": sequence_length,
            }]

        fragments = []

        for index, position in enumerate(
            positions
        ):
            next_position = (
                positions[
                    (index + 1)
                    % len(positions)
                ]
            )

            if index == len(positions) - 1:
                length = (
                    sequence_length
                    - position
                    + positions[0]
                )
            else:
                length = (
                    next_position
                    - position
                )

            fragments.append({
                "fragment": index + 1,
                "length": length,
            })

        return fragments

    if not positions:
        return [{
            "fragment": 1,
            "start": 1,
            "end": sequence_length,
            "length": sequence_length,
        }]

    fragments = []

    previous = 0

    for index, position in enumerate(
        positions,
        start=1,
    ):
        length = position - previous

        if length > 0:
            fragments.append({
                "fragment": len(
                    fragments
                ) + 1,
                "start": previous + 1,
                "end": position,
                "length": length,
            })

        previous = position

    if previous < sequence_length:
        fragments.append({
            "fragment": len(
                fragments
            ) + 1,
            "start": previous + 1,
            "end": sequence_length,
            "length": (
                sequence_length
                - previous
            ),
        })

    return fragments


def _selected_digest(
    sequence,
    enzymes,
    topology,
):
    cuts = []

    for enzyme_name in dict.fromkeys(
        enzymes
    ):
        analysis = analyze_enzyme(
            sequence,
            enzyme_name,
            topology=topology,
        )

        # A fragment boundary is the top-strand cut, i.e. the number of
        # bases 5' of the cut: recognition start (1-based) - 1 + cut_top,
        # taken modulo the length on a circular molecule (audit S17; the
        # recognition start was used before, wrong by cut_top - 1).
        cut_top = int(get_enzyme_data(enzyme_name)["cut_top"])
        for position in analysis[
            "sites"
        ]:
            cut = position - 1 + cut_top
            if topology == "circular" and sequence:
                cut %= len(sequence)
            cuts.append({
                "enzyme": enzyme_name,
                "position": cut,
                "recognition_start": position,
            })

    positions = sorted({
        cut["position"]
        for cut in cuts
    })

    fragments = _digest_fragments(
        len(sequence),
        positions,
        topology=topology,
    )

    return {
        "enzymes": list(
            dict.fromkeys(
                enzymes
            )
        ),
        "cuts": cuts,
        "unique_cut_positions": positions,
        "cut_count": len(cuts),
        "fragment_count": len(
            fragments
        ),
        "fragments": fragments,
        "topology": topology,
    }


# ============================================================
# SCORING
# ============================================================

def _score_strategy(
    validity,
    strategy,
    vector_definition,
    insert_definition,
    vector_internal,
    insert_internal,
):
    components = {}

    components["required_sites"] = (
        20
        if (
            vector_definition["valid"]
            and insert_definition["valid"]
        )
        else 0
    )

    components["unique_end_architecture"] = (
        20
        if (
            vector_internal["valid"]
            and insert_internal["valid"]
        )
        else 0
    )

    components["junction_compatibility"] = (
        25
        if validity[
            "junctions_compatible"
        ]
        else 0
    )

    components["directionality"] = (
        15
        if strategy[
            "directional"
        ]
        else 0
    )

    components["vector_end_orthogonality"] = (
        10
        if not strategy[
            "self_ligation_risk"
        ]
        else 0
    )

    components["compatibility_metadata"] = (
        10
        if validity[
            "compatibility_known"
        ]
        else 0
    )

    total = sum(
        components.values()
    )

    if total >= 90:
        category = "Excellent"
    elif total >= 75:
        category = "Very Good"
    elif total >= 60:
        category = "Good"
    elif total >= 40:
        category = "Moderate"
    else:
        category = "Poor"

    reasons = []

    if components[
        "required_sites"
    ]:
        reasons.append(
            "All selected enzymes define the required two-end architecture."
        )
    else:
        reasons.append(
            "The selected enzymes do not define complete two-end architectures."
        )

    if components[
        "unique_end_architecture"
    ]:
        reasons.append(
            "No additional selected-enzyme sites interfere with the intended ends."
        )
    else:
        reasons.append(
            "Additional selected-enzyme sites may interfere with the design."
        )

    if components[
        "junction_compatibility"
    ]:
        reasons.append(
            "Both intended vector-insert junctions are compatible."
        )
    else:
        reasons.append(
            "One or both intended vector-insert junctions are not compatible."
        )

    if components[
        "directionality"
    ]:
        reasons.append(
            "The intended orientation is compatible and the reverse orientation is blocked."
        )
    elif validity[
        "junctions_compatible"
    ]:
        reasons.append(
            "The intended junctions work, but reverse orientation remains possible."
        )

    if components[
        "vector_end_orthogonality"
    ]:
        reasons.append(
            "The two vector ends are not mutually compatible."
        )
    else:
        reasons.append(
            "The two vector ends are mutually compatible at the sequence level."
        )

    if not validity[
        "compatibility_known"
    ]:
        reasons.append(
            "At least one junction could not be fully resolved from enzyme metadata."
        )

    return {
        "score": total,
        "maximum": 100,
        "percentage": total,
        "category": category,
        "components": components,
        "reasons": reasons,
    }


def score_cloning_strategy(
    strategy,
    validity,
):
    # Backward-compatible 10-point score.
    score = 0
    reasons = []

    if validity.get(
        "all_required_sites",
        validity.get(
            "required_sites",
            False,
        ),
    ):
        score += 2
        reasons.append(
            "All required restriction sites are present."
        )
    else:
        reasons.append(
            "Required restriction sites are missing."
        )

    if validity.get(
        "vector_unique",
        False,
    ) and validity.get(
        "insert_unique",
        False,
    ):
        score += 2
        reasons.append(
            "Selected enzymes define clean two-end architectures."
        )
    else:
        reasons.append(
            "Selected enzymes do not define clean two-end architectures."
        )

    if validity.get(
        "junctions_compatible",
        False,
    ):
        score += 2
        reasons.append(
            "Both cloning junctions are compatible."
        )
    else:
        reasons.append(
            "One or both cloning junctions are incompatible."
        )

    if strategy.get(
        "directional",
        False,
    ):
        score += 3
        reasons.append(
            "The strategy provides directional cloning."
        )
    else:
        reasons.append(
            "The strategy does not provide strong directional control."
        )

    if not strategy.get(
        "self_ligation_risk",
        True,
    ):
        score += 1
        reasons.append(
            "Vector ends are not mutually compatible."
        )
    else:
        reasons.append(
            "Vector ends are mutually compatible."
        )

    category = (
        "Excellent"
        if score >= 8
        else "Good"
        if score >= 6
        else "Moderate"
        if score >= 4
        else "Poor"
    )

    return {
        "score": score,
        "maximum": 10,
        "category": category,
        "reasons": reasons,
    }


# ============================================================
# CANONICAL STRATEGY EVALUATOR
# ============================================================

def evaluate_cloning_strategy(
    vector_sequence,
    insert_sequence,
    vector_left_enzyme,
    vector_right_enzyme,
    insert_left_enzyme,
    insert_right_enzyme,
    vector_topology="linear",
    insert_topology="linear",
):
    """Evaluate one four-enzyme cloning architecture.

    vector_topology / insert_topology ("linear" or "circular") decide
    whether recognition sites spanning a molecule's origin are found. On a
    circular molecule the segment between its two selected sites, in
    coordinate order, is taken as the region replaced (vector) or excised
    (insert); this assumption is listed in the result.
    """
    vector_topology = _normalize_topology(
        vector_topology
    )
    insert_topology = _normalize_topology(
        insert_topology
    )

    vector_validation = _validate_sequence(
        vector_sequence,
        "vector",
    )
    insert_validation = _validate_sequence(
        insert_sequence,
        "insert",
    )

    if not vector_validation["valid"]:
        return {
            "status": "INVALID_VECTOR",
            "valid": False,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
        }

    if not insert_validation["valid"]:
        return {
            "status": "INVALID_INSERT",
            "valid": False,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
        }

    vector = vector_validation[
        "sequence"
    ]
    insert = insert_validation[
        "sequence"
    ]

    vector_definition = _end_definition(
        vector,
        vector_left_enzyme,
        vector_right_enzyme,
        role="vector",
        topology=vector_topology,
    )

    insert_definition = _end_definition(
        insert,
        insert_left_enzyme,
        insert_right_enzyme,
        role="insert",
        topology=insert_topology,
    )

    if not vector_definition[
        "valid"
    ]:
        return {
            "status": "INVALID_VECTOR_ENDS",
            "valid": False,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
            "vector_definition": vector_definition,
            "insert_definition": insert_definition,
            "reason": vector_definition[
                "reason"
            ],
        }

    if not insert_definition[
        "valid"
    ]:
        return {
            "status": "INVALID_INSERT_ENDS",
            "valid": False,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
            "vector_definition": vector_definition,
            "insert_definition": insert_definition,
            "reason": insert_definition[
                "reason"
            ],
        }

    vector_left = vector_definition[
        "left_site"
    ]
    vector_right = vector_definition[
        "right_site"
    ]
    insert_left = insert_definition[
        "left_site"
    ]
    insert_right = insert_definition[
        "right_site"
    ]

    vector_internal = _internal_site_check(
        vector,
        vector_definition,
        topology=vector_topology,
    )
    insert_internal = _internal_site_check(
        insert,
        insert_definition,
        topology=insert_topology,
    )

    left_junction = analyze_junction(
        vector_left,
        insert_left,
        "LEFT",
    )
    right_junction = analyze_junction(
        vector_right,
        insert_right,
        "RIGHT",
    )

    strategy = analyze_directional_strategy(
        vector_left,
        vector_right,
        insert_left,
        insert_right,
    )

    compatibility_known = (
        left_junction[
            "compatible"
        ] is not None
        and right_junction[
            "compatible"
        ] is not None
    )

    junctions_compatible = (
        left_junction[
            "compatible"
        ] is True
        and right_junction[
            "compatible"
        ] is True
    )

    reverse_orientation = (
        strategy[
            "reverse_orientation_compatible"
        ]
    )

    validity = {
        "all_required_sites": (
            vector_definition[
                "valid"
            ]
            and insert_definition[
                "valid"
            ]
        ),
        "required_sites": (
            vector_definition[
                "valid"
            ]
            and insert_definition[
                "valid"
            ]
        ),
        "vector_unique": (
            vector_internal[
                "valid"
            ]
        ),
        "insert_unique": (
            insert_internal[
                "valid"
            ]
        ),
        "no_internal_sites": (
            vector_internal[
                "valid"
            ]
            and insert_internal[
                "valid"
            ]
        ),
        "junctions_compatible": (
            junctions_compatible
        ),
        "compatibility_known": (
            compatibility_known
        ),
        "reverse_orientation_compatible": (
            reverse_orientation
        ),
        "vector_self_compatible": (
            strategy[
                "self_ligation_risk"
            ]
        ),
        "valid": (
            vector_definition[
                "valid"
            ]
            and insert_definition[
                "valid"
            ]
            and vector_internal[
                "valid"
            ]
            and insert_internal[
                "valid"
            ]
            and junctions_compatible
            and compatibility_known
        ),
    }

    score = _score_strategy(
        validity,
        strategy,
        vector_definition,
        insert_definition,
        vector_internal,
        insert_internal,
    )

    legacy_score = score_cloning_strategy(
        strategy,
        validity,
    )

    vector_digest = _selected_digest(
        vector,
        [
            vector_left["enzyme"],
            vector_right["enzyme"],
        ],
        vector_topology,
    )

    insert_digest = _selected_digest(
        insert,
        [
            insert_left["enzyme"],
            insert_right["enzyme"],
        ],
        insert_topology,
    )

    assumptions = [
        f"The {role} is circular: the segment between its two selected "
        "sites, in coordinate order, is taken as the region "
        f"{'replaced' if role == 'vector' else 'excised'}."
        for role, topology in (
            ("vector", vector_topology),
            ("insert", insert_topology),
        )
        if topology == "circular"
    ]

    warnings = []

    if not vector_internal[
        "valid"
    ]:
        warnings.append(
            "Vector contains additional selected-enzyme restriction sites."
        )

    if not insert_internal[
        "valid"
    ]:
        warnings.append(
            "Insert contains additional selected-enzyme restriction sites."
        )

    if not junctions_compatible:
        if left_junction[
            "compatible"
        ] is False:
            warnings.append(
                "Left cloning junction is incompatible."
            )
        if right_junction[
            "compatible"
        ] is False:
            warnings.append(
                "Right cloning junction is incompatible."
            )

    if reverse_orientation:
        warnings.append(
            "The reverse insert orientation is also compatible."
        )

    if strategy[
        "self_ligation_risk"
    ]:
        warnings.append(
            "Vector ends are mutually compatible; "
            "sequence-level self-ligation is possible."
        )

    if (
        vector_definition[
            "mode"
        ] == "same-enzyme"
        or insert_definition[
            "mode"
        ] == "same-enzyme"
    ):
        warnings.append(
            "At least one molecule uses the same enzyme at both ends; "
            "orientation is not enforced by enzyme identity alone."
        )

    if not warnings:
        warnings.append(
            "No major sequence-level design warnings detected."
        )

    return {
        "status": (
            "VALID"
            if validity[
                "valid"
            ]
            else "INVALID"
        ),
        "valid": validity[
            "valid"
        ],
        "mode": "two-enzyme",
        "vector": {
            "length": len(vector),
            "topology": vector_topology,
            "left": vector_left,
            "right": vector_right,
        },
        "insert": {
            "length": len(insert),
            "topology": insert_topology,
            "left": insert_left,
            "right": insert_right,
        },
        "assumptions": assumptions,
        "vector_definition": (
            vector_definition
        ),
        "insert_definition": (
            insert_definition
        ),
        "vector_internal_sites": (
            vector_internal
        ),
        "insert_internal_sites": (
            insert_internal
        ),
        "junctions": {
            "left": left_junction,
            "right": right_junction,
        },
        "reverse_orientation": {
            "left": strategy[
                "reverse_left"
            ],
            "right": strategy[
                "reverse_right"
            ],
            "compatible": (
                reverse_orientation
            ),
        },
        "strategy": {
            "directional": strategy[
                "directional"
            ],
            "orientation": strategy[
                "orientation"
            ],
            "status": strategy[
                "status"
            ],
            "self_ligation_risk": strategy[
                "self_ligation_risk"
            ],
            "vector_self_compatibility": strategy[
                "self_ligation_risk"
            ],
            "insert_self_compatibility": (
                ends_compatible(
                    insert_left,
                    insert_right,
                )
                is True
            ),
        },
        "digest": {
            "vector": vector_digest,
            "insert": insert_digest,
        },
        "validity": validity,
        "score": score,
        "legacy_score": legacy_score,
        "warnings": warnings,
    }


# ============================================================
# BACKWARD-COMPATIBLE TWO-ENZYME API
# ============================================================

def analyze_enzyme_pair(
    vector_sequence,
    insert_sequence,
    vector_left_enzyme,
    vector_right_enzyme,
    insert_left_enzyme,
    insert_right_enzyme,
    vector_topology="linear",
    insert_topology="linear",
):
    return evaluate_cloning_strategy(
        vector_sequence,
        insert_sequence,
        vector_left_enzyme,
        vector_right_enzyme,
        insert_left_enzyme,
        insert_right_enzyme,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )


# ============================================================
# SINGLE-ENZYME ANALYSIS
# ============================================================

def analyze_cloning(
    vector_sequence,
    insert_sequence,
    vector_enzyme,
    insert_enzyme,
    vector_topology="linear",
    insert_topology="linear",
):
    vector_topology = _normalize_topology(
        vector_topology
    )
    insert_topology = _normalize_topology(
        insert_topology
    )

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
            "valid": False,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
        }

    if not insert_validation[
        "valid"
    ]:
        return {
            "status": "INVALID_INSERT",
            "valid": False,
            "vector_validation": vector_validation,
            "insert_validation": insert_validation,
        }

    vector_analysis = analyze_enzyme(
        vector_validation[
            "sequence"
        ],
        vector_enzyme,
        topology=vector_topology,
    )
    insert_analysis = analyze_enzyme(
        insert_validation[
            "sequence"
        ],
        insert_enzyme,
        topology=insert_topology,
    )

    vector_site = (
        make_site(
            vector_analysis,
            vector_analysis[
                "sites"
            ][0],
        )
        if vector_analysis[
            "sites"
        ]
        else make_missing_site(
            vector_analysis
        )
    )

    insert_site = (
        make_site(
            insert_analysis,
            insert_analysis[
                "sites"
            ][0],
        )
        if insert_analysis[
            "sites"
        ]
        else make_missing_site(
            insert_analysis
        )
    )

    compatibility = (
        analyze_end_compatibility(
            vector_site,
            insert_site,
        )
    )

    vector_ambiguity = _ambiguity_reason(
        vector_analysis,
        "vector",
    )
    insert_ambiguity = _ambiguity_reason(
        insert_analysis,
        "insert",
    )
    vector_unique = (
        vector_analysis[
            "site_count"
        ] == 1
        and vector_ambiguity is None
    )
    insert_unique = (
        insert_analysis[
            "site_count"
        ] == 1
        and insert_ambiguity is None
    )

    validity = {
        "all_required_sites": (
            vector_analysis[
                "site_count"
            ] >= 1
            and insert_analysis[
                "site_count"
            ] >= 1
        ),
        "vector_unique": vector_unique,
        "insert_unique": insert_unique,
        "no_internal_sites": (
            vector_unique
            and insert_unique
        ),
        "junctions_compatible": (
            compatibility[
                "compatible"
            ] is True
        ),
        "compatibility_known": (
            compatibility[
                "compatible"
            ] is not None
        ),
        "valid": (
            vector_unique
            and insert_unique
            and compatibility[
                "compatible"
            ] is True
        ),
    }

    strategy = {
        "strategy": "Single-enzyme cloning",
        "orientation": (
            "Non-directional"
            if validity["valid"]
            else "Invalid"
        ),
        "directional": False,
        "self_ligation_risk": True,
        "status": (
            "Compatible ends."
            if compatibility[
                "compatible"
            ]
            else "Ends are incompatible."
        ),
        "left": compatibility,
        "right": compatibility,
    }

    score = score_cloning_strategy(
        strategy,
        validity,
    )

    warnings = []

    if vector_analysis[
        "site_count"
    ] != 1:
        warnings.append(
            "Vector does not contain exactly one selected site."
        )

    if insert_analysis[
        "site_count"
    ] != 1:
        warnings.append(
            "Insert does not contain exactly one selected site."
        )

    for reason in (vector_ambiguity, insert_ambiguity):
        if reason:
            warnings.append(reason)

    if compatibility[
        "compatible"
    ] is not True:
        warnings.append(
            "Vector and insert ends are not confirmed compatible."
        )

    return {
        "mode": "single-enzyme",
        "vector_enzyme": vector_enzyme,
        "insert_enzyme": insert_enzyme,
        "vector": vector_analysis,
        "insert": insert_analysis,
        "vector_site": vector_site,
        "insert_site": insert_site,
        "compatibility": compatibility,
        "strategy": strategy,
        "score": score,
        "validity": validity,
        "warnings": warnings,
    }


# ============================================================
# ENZYME SUITABILITY / COMPATIBILITY HELPERS
# ============================================================

def assess_enzyme_suitability(
    vector_sequence,
    insert_sequence,
    enzyme_name,
    vector_topology="linear",
    insert_topology="linear",
):
    vector = analyze_enzyme(
        vector_sequence,
        enzyme_name,
        topology=vector_topology,
    )
    insert = analyze_enzyme(
        insert_sequence,
        enzyme_name,
        topology=insert_topology,
    )

    ambiguity = [
        reason
        for reason in (
            _ambiguity_reason(vector, "vector"),
            _ambiguity_reason(insert, "insert"),
        )
        if reason
    ]

    suitable = (
        vector["site_count"] == 1
        and insert["site_count"] == 1
        and not ambiguity
        and (
            get_overhang(
                enzyme_name
            )
            or is_blunt_enzyme(
                enzyme_name
            )
        )
    )

    reasons = []

    if vector[
        "site_count"
    ] == 1:
        reasons.append(
            "Vector contains one selected site."
        )
    elif vector[
        "site_count"
    ] == 0:
        reasons.append(
            "Vector does not contain the selected site."
        )
    else:
        reasons.append(
            "Vector contains multiple selected sites."
        )

    if insert[
        "site_count"
    ] == 1:
        reasons.append(
            "Insert contains one selected site."
        )
    elif insert[
        "site_count"
    ] == 0:
        reasons.append(
            "Insert does not contain the selected site."
        )
    else:
        reasons.append(
            "Insert contains multiple selected sites."
        )

    reasons.extend(ambiguity)

    return {
        "enzyme": enzyme_name,
        "vector": vector,
        "insert": insert,
        "suitable": suitable,
        "reasons": reasons,
    }


def find_compatible_enzymes(
    vector_sequence,
    insert_sequence,
    enzyme_names=None,
    vector_topology="linear",
    insert_topology="linear",
):
    names = (
        list(
            RESTRICTION_ENZYMES.keys()
        )
        if enzyme_names is None
        else [
            name
            for name in enzyme_names
            if name in RESTRICTION_ENZYMES
        ]
    )

    return [
        assess_enzyme_suitability(
            vector_sequence,
            insert_sequence,
            name,
            vector_topology=vector_topology,
            insert_topology=insert_topology,
        )
        for name in names
    ]
