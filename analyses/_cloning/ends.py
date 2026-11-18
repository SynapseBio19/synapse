"""Cut ends, their compatibility and the resulting junctions."""

from __future__ import annotations

from data.restriction_enzymes import RESTRICTION_ENZYMES

from analyses._cloning.sites import analyze_enzyme, get_overhang_polarity, reverse_complement


# ============================================================
# SITE OBJECTS
# ============================================================

def make_site(
    analysis,
    position,
):
    return {
        "enzyme": analysis["enzyme"],
        "position": position,
        "recognition_sequence": (
            analysis[
                "recognition_sequence"
            ]
        ),
        "cut_type": analysis[
            "cut_type"
        ],
        "overhang": (
            analysis.get(
                "overhang",
                "",
            )
            or ""
        ).upper(),
        "overhang_type": analysis.get(
            "overhang_type",
            "Unknown",
        ),
        "overhang_length": analysis.get(
            "overhang_length",
            0,
        ),
        "overhang_polarity": analysis.get(
            "overhang_polarity",
            "unknown",
        ),
        "is_blunt": bool(
            analysis.get(
                "is_blunt",
                False,
            )
        ),
    }


def make_missing_site(
    analysis,
):
    return {
        "enzyme": analysis["enzyme"],
        "position": None,
        "recognition_sequence": (
            analysis[
                "recognition_sequence"
            ]
        ),
        "cut_type": analysis[
            "cut_type"
        ],
        "overhang": (
            analysis.get(
                "overhang",
                "",
            )
            or ""
        ),
        "overhang_type": analysis.get(
            "overhang_type",
            "Unknown",
        ),
        "overhang_length": analysis.get(
            "overhang_length",
            0,
        ),
        "overhang_polarity": analysis.get(
            "overhang_polarity",
            "unknown",
        ),
        "is_blunt": bool(
            analysis.get(
                "is_blunt",
                False,
            )
        ),
    }


# ============================================================
# END COMPATIBILITY
# ============================================================

def _end_signature(site):
    if site is None:
        return {
            "known": False,
            "kind": "INVALID",
            "overhang": "",
        }

    if site.get(
        "is_blunt",
        False,
    ):
        return {
            "known": True,
            "kind": "BLUNT",
            "overhang": "",
        }

    overhang = (
        site.get(
            "overhang",
            "",
        )
        or ""
    ).upper()

    if not overhang:
        return {
            "known": False,
            "kind": "UNKNOWN",
            "overhang": "",
        }

    return {
        "known": True,
        "kind": "COHESIVE",
        "overhang": overhang,
        "polarity": _site_polarity(site),
    }


def _site_polarity(site):
    """Overhang polarity of a site dict ('5_prime', '3_prime', 'blunt', 'unknown').

    Sites created before polarity was recorded fall back to the enzyme
    database entry.
    """
    polarity = str(site.get("overhang_polarity") or "").lower()
    if polarity in {"5_prime", "3_prime", "blunt"}:
        return polarity

    enzyme_name = site.get("enzyme")
    if enzyme_name in RESTRICTION_ENZYMES:
        return get_overhang_polarity(enzyme_name)

    return "unknown"


def ends_compatible(
    first,
    second,
):
    if first is None or second is None:
        return False

    first_signature = _end_signature(
        first
    )
    second_signature = _end_signature(
        second
    )

    if not first_signature["known"]:
        return None

    if not second_signature["known"]:
        return None

    if (
        first_signature["kind"]
        == "BLUNT"
        and second_signature["kind"]
        == "BLUNT"
    ):
        return True

    if (
        first_signature["kind"]
        == "BLUNT"
        or second_signature["kind"]
        == "BLUNT"
    ):
        return False

    first_overhang = (
        first_signature["overhang"]
    )
    second_overhang = (
        second_signature["overhang"]
    )

    if not first_overhang or not second_overhang:
        return None

    # End chemistry: a 5' overhang can only anneal to a 5' overhang and a
    # 3' overhang only to a 3' overhang (e.g. HindIII 5'-AGCT is NOT
    # compatible with SacI 3'-AGCT). Unknown polarity -> undetermined.
    first_polarity = first_signature.get("polarity", "unknown")
    second_polarity = second_signature.get("polarity", "unknown")

    if "unknown" in (first_polarity, second_polarity):
        return None

    if first_polarity != second_polarity:
        return False

    return (
        first_overhang == second_overhang
        or first_overhang
        == reverse_complement(
            second_overhang
        )
    )


def analyze_end_compatibility(
    vector_site,
    insert_site,
):
    compatible = ends_compatible(
        vector_site,
        insert_site,
    )

    if compatible is True:
        status = "COMPATIBLE"
        mechanism = (
            "Compatible cohesive/blunt DNA ends."
        )
    elif compatible is False:
        status = "INCOMPATIBLE"
        mechanism = (
            "The selected DNA ends are incompatible."
        )
    else:
        status = "UNKNOWN"
        mechanism = (
            "Compatibility cannot be established "
            "from the available enzyme metadata."
        )

    return {
        "compatible": compatible,
        "status": status,
        "type": (
            "Blunt"
            if (
                vector_site
                and insert_site
                and vector_site.get(
                    "is_blunt",
                    False,
                )
                and insert_site.get(
                    "is_blunt",
                    False,
                )
            )
            else "Cohesive"
            if compatible is not None
            else "Unknown"
        ),
        "mechanism": mechanism,
        "vector_end": (
            vector_site.get(
                "overhang",
                "",
            )
            or "-"
            if vector_site
            else "-"
        ),
        "insert_end": (
            insert_site.get(
                "overhang",
                "",
            )
            or "-"
            if insert_site
            else "-"
        ),
    }


# ============================================================
# JUNCTIONS
# ============================================================

def analyze_junction(
    vector_site,
    insert_site,
    side,
):
    compatibility = (
        analyze_end_compatibility(
            vector_site,
            insert_site,
        )
    )

    return {
        "side": side,
        "compatible": compatibility[
            "compatible"
        ],
        "status": compatibility[
            "status"
        ],
        "mechanism": compatibility[
            "mechanism"
        ],
        "vector_enzyme": vector_site[
            "enzyme"
        ],
        "insert_enzyme": insert_site[
            "enzyme"
        ],
        "vector_position": vector_site[
            "position"
        ],
        "insert_position": insert_site[
            "position"
        ],
        "vector_overhang": (
            vector_site.get(
                "overhang",
                "",
            )
            or "-"
        ),
        "insert_overhang": (
            insert_site.get(
                "overhang",
                "",
            )
            or "-"
        ),
        "vector_end_type": (
            "BLUNT"
            if vector_site.get(
                "is_blunt",
                False,
            )
            else "COHESIVE"
            if vector_site.get(
                "overhang",
                "",
            )
            else "UNKNOWN"
        ),
        "insert_end_type": (
            "BLUNT"
            if insert_site.get(
                "is_blunt",
                False,
            )
            else "COHESIVE"
            if insert_site.get(
                "overhang",
                "",
            )
            else "UNKNOWN"
        ),
    }


# ============================================================
# TWO-END DEFINITIONS
# ============================================================

def _ambiguity_reason(analysis, role):
    """Explain why POSSIBLE sites leave an enzyme's site count open, else None.

    An ambiguity-dependent site may or may not exist, so an enzyme with one
    is neither a confirmed unique cutter nor a confirmed non-cutter -- the
    rule restriction.py applies to unique cutters and enzyme pairs.
    """
    count = analysis.get("possible_site_count", 0)
    if not count:
        return None
    positions = ", ".join(
        str(position)
        for position in analysis.get("possible_sites", [])
    )
    return (
        f"{analysis['enzyme']} has {count} ambiguity-dependent (POSSIBLE) "
        f"site(s) in the {role} at position(s) {positions}: IUPAC symbols "
        "inside the recognition sequence leave the number of real sites "
        "undetermined."
    )


def _end_definition(
    sequence,
    left_enzyme,
    right_enzyme,
    role="molecule",
    topology="linear",
):
    left_analysis = analyze_enzyme(
        sequence,
        left_enzyme,
        topology=topology,
    )
    right_analysis = analyze_enzyme(
        sequence,
        right_enzyme,
        topology=topology,
    )

    ambiguity = [
        reason
        for reason in (
            _ambiguity_reason(left_analysis, role),
            _ambiguity_reason(right_analysis, role)
            if right_enzyme != left_enzyme
            else None,
        )
        if reason
    ]

    if ambiguity:
        return {
            "valid": False,
            "mode": (
                "same-enzyme"
                if left_enzyme == right_enzyme
                else "two-enzyme"
            ),
            "reason": " ".join(ambiguity),
            "ambiguity_dependent": True,
            "left_analysis": left_analysis,
            "right_analysis": right_analysis,
            "left_site": None,
            "right_site": None,
        }

    if (
        left_enzyme
        == right_enzyme
    ):
        sites = left_analysis[
            "sites"
        ]

        if len(sites) != 2:
            return {
                "valid": False,
                "mode": "same-enzyme",
                "reason": (
                    "The same enzyme must occur "
                    "exactly twice to define two "
                    "physical cloning ends."
                ),
                "left_analysis": left_analysis,
                "right_analysis": right_analysis,
                "left_site": None,
                "right_site": None,
            }

        return {
            "valid": True,
            "mode": "same-enzyme",
            "reason": (
                "Exactly two recognition sites "
                "define the two physical ends."
            ),
            "left_analysis": left_analysis,
            "right_analysis": right_analysis,
            "left_site": make_site(
                left_analysis,
                sites[0],
            ),
            "right_site": make_site(
                left_analysis,
                sites[1],
            ),
        }

    left_sites = left_analysis[
        "sites"
    ]
    right_sites = right_analysis[
        "sites"
    ]

    if (
        len(left_sites) != 1
        or len(right_sites) != 1
    ):
        return {
            "valid": False,
            "mode": "two-enzyme",
            "reason": (
                "Different enzymes must each "
                "occur exactly once to define "
                "two clean physical ends."
            ),
            "left_analysis": left_analysis,
            "right_analysis": right_analysis,
            "left_site": (
                make_site(
                    left_analysis,
                    left_sites[0],
                )
                if len(left_sites) == 1
                else None
            ),
            "right_site": (
                make_site(
                    right_analysis,
                    right_sites[0],
                )
                if len(right_sites) == 1
                else None
            ),
        }

    left_site = make_site(
        left_analysis,
        left_sites[0],
    )
    right_site = make_site(
        right_analysis,
        right_sites[0],
    )

    if left_site["position"] > right_site[
        "position"
    ]:
        left_site, right_site = (
            right_site,
            left_site,
        )

    return {
        "valid": True,
        "mode": "two-enzyme",
        "reason": (
            "Each selected enzyme occurs "
            "exactly once."
        ),
        "left_analysis": left_analysis,
        "right_analysis": right_analysis,
        "left_site": left_site,
        "right_site": right_site,
    }
