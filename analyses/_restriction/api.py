"""The combined analysis, summaries and the get_*/run_* entry points."""

from __future__ import annotations

from data.restriction_enzymes import get_restriction_enzyme, RESTRICTION_ENZYMES

from analyses._restriction.sites import (
    analyze_enzyme,
    calculate_cut_positions,
    _normalize_restriction_sequence,
    _require_valid_restriction_sequence,
    _restriction_sequence_status,
)
from analyses._restriction.digest import (
    calculate_gel_bands,
    create_restriction_map,
    find_internal_sites,
    find_non_cutters,
    find_unique_cutters,
    perform_digest,
)
from analyses._restriction.compatibility import analyze_double_digest, recommend_enzyme_pairs


def analyze_restriction(
    sequence,
    enzyme_names=None,
    topology="linear",
):
    """Complete restriction analysis used by the Synapse pipeline."""
    sequence = _require_valid_restriction_sequence(sequence)

    topology = str(topology).lower()
    if topology not in {"linear", "circular"}:
        raise ValueError(
            "Topology must be 'linear' or 'circular'."
        )

    if enzyme_names is None:
        enzyme_names = list(RESTRICTION_ENZYMES.keys())

    enzyme_names = list(dict.fromkeys(enzyme_names))

    unknown = [
        name
        for name in enzyme_names
        if get_restriction_enzyme(name) is None
    ]
    if unknown:
        raise ValueError(
            "Unknown restriction enzyme(s): "
            + ", ".join(unknown)
        )

    enzyme_results = {
        enzyme_name: analyze_enzyme(
            sequence,
            enzyme_name,
            topology=topology,
        )
        for enzyme_name in enzyme_names
    }

    restriction_map = create_restriction_map(
        sequence,
        enzyme_names,
        topology=topology,
    )

    digest = perform_digest(
        sequence,
        enzyme_names,
        topology,
    )

    gel_bands = calculate_gel_bands(
        digest["fragments"]
    )

    unique_cutters = find_unique_cutters(
        sequence,
        enzyme_names,
        topology=topology,
    )

    non_cutters = find_non_cutters(
        sequence,
        enzyme_names,
        topology=topology,
    )

    internal_sites = find_internal_sites(
        sequence,
        enzyme_names,
        topology=topology,
    )

    recommended_pairs = recommend_enzyme_pairs(
        sequence,
        enzyme_names,
        top_n=5,
        topology=topology,
    )

    confirmed_site_total = sum(
        result["confirmed_site_count"]
        for result in enzyme_results.values()
    )
    possible_site_total = sum(
        result["possible_site_count"]
        for result in enzyme_results.values()
    )

    return {
        "sequence_length": len(sequence),
        "topology": topology,
        "enzymes": enzyme_names,
        "enzyme_results": enzyme_results,
        "restriction_map": restriction_map,
        "digest": digest,
        "gel_bands": gel_bands,
        "unique_cutters": unique_cutters,
        "non_cutters": non_cutters,
        "internal_sites": internal_sites,
        "recommended_pairs": recommended_pairs,
        "site_statistics": {
            "total_site_events": len(restriction_map),
            "confirmed_site_events": confirmed_site_total,
            "possible_site_events": possible_site_total,
            "physical_cut_boundaries": digest[
                "unique_cut_count"
            ],
        },
        "sequence_status": _restriction_sequence_status(
            sequence
        ),
    }


def summarize_enzyme_sites(sequence, enzyme_name, topology="linear"):
    analysis = analyze_enzyme(
        sequence,
        enzyme_name,
        topology=topology,
    )

    return {
        "enzyme": enzyme_name,
        "site_count": analysis["site_count"],
        "confirmed_site_count": analysis[
            "confirmed_site_count"
        ],
        "possible_site_count": analysis[
            "possible_site_count"
        ],
        "confirmed_positions": [
            site["site_position"]
            for site in analysis["sites"]
            if site["match_status"] == "CONFIRMED"
        ],
        "possible_positions": [
            site["site_position"]
            for site in analysis["sites"]
            if site["match_status"] == "POSSIBLE"
        ],
    }


def validate_enzyme_database(enzyme_names=None):
    """Validate all metadata required by the restriction engine."""
    if enzyme_names is None:
        enzyme_names = list(RESTRICTION_ENZYMES.keys())

    required = {
        "recognition_sequence",
        "cut_top",
        "cut_bottom",
        "is_blunt",
        "overhang_type",
        "overhang_length",
    }

    valid = []
    invalid = []

    for name in enzyme_names:
        enzyme = get_restriction_enzyme(name)

        if enzyme is None:
            invalid.append({
                "enzyme": name,
                "missing": ["DATABASE_ENTRY"],
            })
            continue

        missing = sorted(
            field
            for field in required
            if field not in enzyme
        )

        problems = list(missing)

        if not missing:
            recognition = _normalize_restriction_sequence(
                enzyme.get("recognition_sequence", "")
            )

            try:
                calculate_cut_positions(
                    0,
                    enzyme,
                )
            except Exception as exc:
                problems.append(
                    f"INVALID_CUT_DATA: {exc}"
                )

            try:
                declared_overhang = int(
                    enzyme.get(
                        "overhang_length",
                        0,
                    )
                )
                calculated_overhang = calculate_cut_positions(
                    0,
                    enzyme,
                )["overhang_length"]
                if declared_overhang != calculated_overhang:
                    problems.append(
                        "OVERHANG_LENGTH_MISMATCH"
                    )
            except (TypeError, ValueError):
                problems.append(
                    "INVALID_OVERHANG_LENGTH"
                )

            if bool(enzyme.get("is_blunt", False)) != (
                calculate_cut_positions(0, enzyme)["overhang_length"] == 0
            ):
                problems.append(
                    "BLUNT_FLAG_MISMATCH"
                )

            if not recognition:
                problems.append(
                    "EMPTY_RECOGNITION_SEQUENCE"
                )

        if problems:
            invalid.append({
                "enzyme": name,
                "missing": problems,
            })
        else:
            valid.append(name)

    return {
        "valid": valid,
        "invalid": invalid,
        "valid_count": len(valid),
        "invalid_count": len(invalid),
    }


def compare_single_enzyme_digest(
    sequence,
    enzyme_name,
    topology="linear",
):
    result = perform_digest(
        sequence,
        enzyme_name,
        topology=topology,
    )

    if result["cut_count"] == 0:
        interpretation = "No cut detected."
    elif result["unique_cut_count"] == 1:
        interpretation = "One physical cut boundary detected."
    else:
        interpretation = (
            f"{result['unique_cut_count']} physical cut boundaries detected."
        )

    return {
        "digest": result,
        "interpretation": interpretation,
        "gel_bands": calculate_gel_bands(
            result["fragments"]
        ),
    }


def run_restriction_analysis(
    sequence,
    enzyme_names=None,
    topology="linear",
):
    """Canonical non-printing pipeline entry point."""
    if not isinstance(sequence, str):
        return {
            "status": "INVALID_SEQUENCE",
            "valid": False,
            "analysis": None,
        }

    normalized = _normalize_restriction_sequence(sequence)

    status = _restriction_sequence_status(
        normalized
    )

    if not status["valid"]:
        return {
            "status": "INVALID_SEQUENCE",
            "valid": False,
            "sequence_length": len(normalized),
            "sequence_status": status,
            "analysis": None,
        }

    topology = str(topology).lower()

    if topology not in {"linear", "circular"}:
        return {
            "status": "INVALID_TOPOLOGY",
            "valid": False,
            "sequence_length": len(normalized),
            "sequence_status": status,
            "analysis": None,
        }

    try:
        analysis = analyze_restriction(
            normalized,
            enzyme_names=enzyme_names,
            topology=topology,
        )
    except ValueError as exc:
        return {
            "status": "ANALYSIS_ERROR",
            "valid": False,
            "sequence_length": len(normalized),
            "sequence_status": status,
            "error": str(exc),
            "analysis": None,
        }

    return {
        "status": "OK",
        "valid": True,
        "sequence_length": len(normalized),
        "topology": topology,
        "sequence_status": status,
        "analysis": analysis,
    }


def get_restriction_summary(result):
    if not result or not result.get("valid"):
        return {
            "status": (
                result.get("status", "NO_ANALYSIS")
                if result
                else "NO_ANALYSIS"
            ),
            "valid": False,
        }

    analysis = result["analysis"]

    return {
        "status": "OK",
        "valid": True,
        "sequence_length": analysis["sequence_length"],
        "topology": analysis["topology"],
        "enzymes_analyzed": len(analysis["enzymes"]),
        "total_sites": len(
            analysis["restriction_map"]
        ),
        "confirmed_sites": analysis[
            "site_statistics"
        ]["confirmed_site_events"],
        "possible_sites": analysis[
            "site_statistics"
        ]["possible_site_events"],
        "cut_count": analysis["digest"]["cut_count"],
        "unique_cut_count": analysis["digest"][
            "unique_cut_count"
        ],
        "fragment_count": analysis["digest"][
            "fragment_count"
        ],
        "unique_cutters": len(
            analysis["unique_cutters"]
        ),
        "non_cutters": len(
            analysis["non_cutters"]
        ),
        "multiple_site_enzymes": len(
            analysis["internal_sites"]
        ),
        "recommended_pairs": len(
            analysis["recommended_pairs"]
        ),
    }


def get_restriction_sites_for_enzyme(
    sequence,
    enzyme_name,
    topology="linear",
):
    result = analyze_enzyme(
        sequence,
        enzyme_name,
        topology=topology,
    )

    return {
        "enzyme": result["enzyme"],
        "recognition_sequence": result[
            "recognition_sequence"
        ],
        "site_count": result["site_count"],
        "confirmed_site_count": result[
            "confirmed_site_count"
        ],
        "possible_site_count": result[
            "possible_site_count"
        ],
        "sites": result["sites"],
        "cut_type": result["cut_type"],
        "overhang_type": result["overhang_type"],
        "overhang_length": result["overhang_length"],
    }


def get_restriction_map(sequence, enzyme_names=None, topology="linear"):
    return create_restriction_map(
        sequence,
        enzyme_names=enzyme_names,
        topology=topology,
    )


def get_digest_prediction(
    sequence,
    enzyme_names,
    topology="linear",
):
    return perform_digest(
        sequence,
        enzyme_names,
        topology=topology,
    )


def get_gel_band_prediction(digest_result):
    if not digest_result or "fragments" not in digest_result:
        return []
    return calculate_gel_bands(
        digest_result["fragments"]
    )


def get_unique_cutters(sequence, enzyme_names=None, topology="linear"):
    return find_unique_cutters(
        sequence,
        enzyme_names=enzyme_names,
        topology=topology,
    )


def get_non_cutters(sequence, enzyme_names=None, topology="linear"):
    return find_non_cutters(
        sequence,
        enzyme_names=enzyme_names,
        topology=topology,
    )


def get_enzyme_pair_recommendations(
    sequence,
    enzyme_names=None,
    top_n=5,
    topology="linear",
):
    return recommend_enzyme_pairs(
        sequence,
        enzyme_names=enzyme_names,
        top_n=top_n,
        topology=topology,
    )


def run_double_digest(
    sequence,
    enzyme_a,
    enzyme_b,
    topology="linear",
):
    try:
        digest = analyze_double_digest(
            sequence,
            enzyme_a,
            enzyme_b,
            topology=topology,
        )
    except ValueError as exc:
        return {
            "status": "INVALID",
            "valid": False,
            "error": str(exc),
            "digest": None,
            "gel_bands": [],
        }

    return {
        "status": "OK",
        "valid": True,
        "digest": digest,
        "gel_bands": calculate_gel_bands(
            digest["fragments"]
        ),
    }


def get_cloning_compatible_sites(
    vector_sequence,
    insert_sequence,
    enzyme_name,
    vector_topology="linear",
    insert_topology="linear",
):
    """Restriction-to-Cloning handoff with explicit confirmed-site counts."""
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

    vector_unique = (
        vector["site_count"] == 1
        and vector["confirmed_site_count"] == 1
    )
    insert_unique = (
        insert["site_count"] == 1
        and insert["confirmed_site_count"] == 1
    )

    return {
        "enzyme": enzyme_name,
        "vector": vector,
        "insert": insert,
        "vector_unique": vector_unique,
        "insert_unique": insert_unique,
        "ready_for_cloning": (
            vector_unique
            and insert_unique
        ),
    }
