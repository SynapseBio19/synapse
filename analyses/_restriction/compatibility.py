"""Sticky-end compatibility, enzyme-pair recommendations and double digests."""

from __future__ import annotations

from data.restriction_enzymes import get_restriction_enzyme, RESTRICTION_ENZYMES
from itertools import combinations

from analyses._restriction.sites import (
    analyze_enzyme,
    _normalize_restriction_sequence,
    _require_valid_restriction_sequence,
    reverse_complement,
)
from analyses._restriction.digest import perform_digest


def get_sticky_end_sequence(site):
    if site.get("is_blunt", False):
        return None
    return site.get("overhang_sequence", "") or None


def _site_overhang_polarity(site):
    """Return '5_prime', '3_prime', 'blunt' or 'unknown' for a site dict."""
    if site.get("is_blunt", False):
        return "blunt"
    top = site.get("top_cut_0based")
    bottom = site.get("bottom_cut_0based")
    if isinstance(top, int) and isinstance(bottom, int) and top != bottom:
        return "5_prime" if top < bottom else "3_prime"
    declared = str(site.get("overhang_type", "")).lower()
    if declared in {"5_prime", "3_prime", "blunt"}:
        return declared
    enzyme = get_restriction_enzyme(site.get("enzyme")) if site.get("enzyme") else None
    if enzyme:
        declared = str(enzyme.get("overhang_type", "")).lower()
        if declared in {"5_prime", "3_prime", "blunt"}:
            return declared
    return "unknown"


def sticky_ends_are_compatible(first_site, second_site):
    """Evaluate sticky-end compatibility with explicit orientation handling."""
    first_blunt = bool(
        first_site.get("is_blunt", False)
    )
    second_blunt = bool(
        second_site.get("is_blunt", False)
    )

    if first_blunt and second_blunt:
        return True

    if first_blunt or second_blunt:
        return False

    first = str(
        first_site.get(
            "overhang_sequence",
            "",
        )
        or ""
    ).upper()

    second = str(
        second_site.get(
            "overhang_sequence",
            "",
        )
        or ""
    ).upper()

    if not first or not second:
        return False

    if len(first) != len(second):
        return False

    # 5' and 3' overhangs never anneal to each other, even with identical
    # letters (HindIII 5'-AGCT vs SacI 3'-AGCT). Audit P0-3.
    first_polarity = _site_overhang_polarity(first_site)
    second_polarity = _site_overhang_polarity(second_site)
    if (
        first_polarity in {"5_prime", "3_prime"}
        and second_polarity in {"5_prime", "3_prime"}
        and first_polarity != second_polarity
    ):
        return False

    # Depending on how a database stores the forward-strand overhang,
    # compatible mating ends may appear as the same sequence or as
    # reverse complements.
    return (
        first == second
        or first == reverse_complement(second)
    )


def analyze_cloning_compatibility(vector_site, insert_site):
    compatible = sticky_ends_are_compatible(
        vector_site,
        insert_site,
    )

    vector_blunt = bool(
        vector_site.get("is_blunt", False)
    )
    insert_blunt = bool(
        insert_site.get("is_blunt", False)
    )

    if vector_blunt and insert_blunt:
        mechanism = "Blunt-end ligation"
        compatibility_type = "BLUNT"
    elif compatible:
        mechanism = "Compatible sticky-end ligation"
        compatibility_type = "STICKY"
    else:
        mechanism = "Incompatible ends"
        compatibility_type = "INCOMPATIBLE"

    return {
        "compatible": compatible,
        "type": compatibility_type,
        "mechanism": mechanism,
        "vector_overhang": vector_site.get(
            "overhang_sequence",
            "",
        ),
        "insert_overhang": insert_site.get(
            "overhang_sequence",
            "",
        ),
        "vector_enzyme": vector_site.get("enzyme"),
        "insert_enzyme": insert_site.get("enzyme"),
    }


def recommend_enzyme_pairs(sequence, enzyme_names=None, top_n=5, topology="linear"):
    """Rank distinct enzyme pairs by transparent sequence-level criteria."""
    sequence = _require_valid_restriction_sequence(sequence)

    if enzyme_names is None:
        enzyme_names = list(RESTRICTION_ENZYMES.keys())

    candidates = []

    for enzyme_name in dict.fromkeys(enzyme_names):
        enzyme = get_restriction_enzyme(enzyme_name)
        if enzyme is None:
            continue

        analysis = analyze_enzyme(
            sequence,
            enzyme_name,
            topology=topology,
        )

        candidates.append({
            "name": enzyme_name,
            "enzyme": enzyme,
            "analysis": analysis,
            "sites": analysis["sites"],
            "site_count": analysis["site_count"],
        })

    recommendations = []

    for first, second in combinations(
        candidates,
        2,
    ):
        # A "recommended enzyme pair" is a cloning/double-digest
        # recommendation, so both enzymes must actually cut the
        # sequence exactly once at confirmed sites. Non-cutters,
        # multi-cutters and ambiguity-dependent cutters are not
        # valid recommendations.
        if not all(
            candidate["site_count"] == 1
            and candidate["analysis"]["confirmed_site_count"] == 1
            and candidate["analysis"]["possible_site_count"] == 0
            for candidate in (first, second)
        ):
            continue

        score = 10
        reasons = []
        warnings = []

        for candidate in (first, second):
            analysis = candidate["analysis"]

            if (
                analysis["site_count"] == 1
                and analysis["confirmed_site_count"] == 1
            ):
                reasons.append(
                    f"{candidate['name']} cuts exactly once at a confirmed site"
                )
            elif analysis["site_count"] == 0:
                warnings.append(
                    f"{candidate['name']} does not cut"
                )
            elif analysis["possible_site_count"] > 0:
                warnings.append(
                    f"{candidate['name']} has ambiguity-dependent site(s)"
                )
            else:
                warnings.append(
                    f"{candidate['name']} cuts "
                    f"{analysis['site_count']} times"
                )

        first_cut = (
            first["sites"][0]["top_cut_0based"]
            if first["site_count"] == 1
            else None
        )
        second_cut = (
            second["sites"][0]["top_cut_0based"]
            if second["site_count"] == 1
            else None
        )

        distinct_physical_cuts = (
            first_cut is not None
            and second_cut is not None
            and first_cut != second_cut
        )

        if distinct_physical_cuts:
            score += 3
            reasons.append(
                "Two distinct physical cut positions"
            )
        elif (
            first_cut is not None
            and second_cut is not None
        ):
            warnings.append(
                "Both enzymes cut at the same physical position"
            )

        first_sticky = not bool(
            first["enzyme"].get("is_blunt", False)
        )
        second_sticky = not bool(
            second["enzyme"].get("is_blunt", False)
        )

        if first_sticky and second_sticky:
            score += 2
            reasons.append(
                "Both enzymes generate sticky ends"
            )
        elif not first_sticky and not second_sticky:
            reasons.append(
                "Both enzymes generate blunt ends"
            )
        else:
            warnings.append(
                "Mixed blunt/sticky cut types"
            )

        recognition_a = _normalize_restriction_sequence(
            first["enzyme"].get(
                "recognition_sequence",
                "",
            )
        )
        recognition_b = _normalize_restriction_sequence(
            second["enzyme"].get(
                "recognition_sequence",
                "",
            )
        )

        if recognition_a != recognition_b:
            score += 1
            reasons.append(
                "Distinct recognition sequences"
            )

        if first_sticky and second_sticky:
            overhang_a = str(
                first["enzyme"].get(
                    "overhang_length",
                    0,
                )
            )
            overhang_b = str(
                second["enzyme"].get(
                    "overhang_length",
                    0,
                )
            )
            if overhang_a == overhang_b:
                reasons.append(
                    "Both enzymes generate sticky ends of the same length"
                )
            else:
                reasons.append(
                    "Sticky-end lengths differ; compatibility must be evaluated by sequence"
                )

        recommendations.append({
            "enzyme_a": first["name"],
            "enzyme_b": second["name"],
            "score": score,
            "maximum": 16,
            "reasons": reasons,
            "warnings": warnings,
            "site_count_a": first["site_count"],
            "site_count_b": second["site_count"],
            "distinct_cut_positions": distinct_physical_cuts,
            "cut_position_a_0based": first_cut,
            "cut_position_b_0based": second_cut,
            "sticky_ends": first_sticky and second_sticky,
        })

    recommendations.sort(
        key=lambda item: (
            item["score"],
            item["distinct_cut_positions"],
            -len(item["warnings"]),
        ),
        reverse=True,
    )

    return recommendations[:max(0, int(top_n))]


def analyze_double_digest(
    sequence,
    enzyme_a,
    enzyme_b,
    topology="linear",
):
    if enzyme_a == enzyme_b:
        raise ValueError(
            "For a double digest, select two different enzymes."
        )

    return perform_digest(
        sequence,
        [enzyme_a, enzyme_b],
        topology,
    )
