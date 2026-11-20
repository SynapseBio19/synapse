"""Restriction maps, fragments, digests, gel bands and cutter lists."""

from __future__ import annotations

from data.restriction_enzymes import (
    get_restriction_enzyme,
    get_type_iis_enzyme,
    RESTRICTION_ENZYMES,
)

from analyses._restriction.sites import (
    analyze_enzyme,
    analyze_restriction_site,
    calculate_cut_positions,
    find_restriction_site_matches,
    find_restriction_sites,
    find_type_iis_sites,
    _require_valid_restriction_sequence,
)


def create_restriction_map(sequence, enzyme_names=None, topology="linear"):
    """Create a complete site map with deterministic ordering.

    Each entry is an analyze_restriction_site() result plus "wraps_origin"
    (True for a site spanning the origin of a circular sequence).
    """
    sequence = _require_valid_restriction_sequence(sequence)

    if enzyme_names is None:
        enzyme_names = list(RESTRICTION_ENZYMES.keys())

    enzyme_names = list(dict.fromkeys(enzyme_names))
    sites = []

    for enzyme_name in enzyme_names:
        if get_restriction_enzyme(enzyme_name) is None:
            continue

        for match in find_restriction_site_matches(
            sequence,
            enzyme_name,
            topology=topology,
        ):
            site = analyze_restriction_site(
                sequence,
                enzyme_name,
                match["position"],
                topology=topology,
            )
            site["wraps_origin"] = match["wraps_origin"]
            sites.append(site)

    sites.sort(
        key=lambda site: (
            site["site_position_0based"],
            site["enzyme"],
        )
    )

    return sites


def get_cut_positions(sequence, enzyme_names, include_possible=False, topology="linear"):
    """Collect enzyme cut events while preserving staggered-cut metadata.

    Only CONFIRMED sites (no IUPAC ambiguity inside the recognition site)
    produce cut events by default. A POSSIBLE site, e.g. a run of N, is not
    evidence of a physical cut and would otherwise create phantom fragments
    (audit P0-5). Pass include_possible=True for a worst-case view; each
    event carries its "status".
    """
    sequence = _require_valid_restriction_sequence(sequence)

    if isinstance(enzyme_names, str):
        enzyme_names = [enzyme_names]

    cuts = []

    for enzyme_name in dict.fromkeys(enzyme_names):
        if get_type_iis_enzyme(enzyme_name) is not None:
            for site in find_type_iis_sites(sequence, enzyme_name, topology=topology):
                if site["status"] != "CONFIRMED" and not include_possible:
                    continue
                if not site["cuts_inside"]:
                    continue
                cuts.append({
                    "position": site["top_cut_0based"],
                    "position_0based": site["top_cut_0based"],
                    "position_1based": site["top_cut"],
                    "enzyme": enzyme_name,
                    "site_position": site["site_position"],
                    "site_position_0based": site["site_position_0based"],
                    "site_position_1based": site["site_position"],
                    "top_cut": site["top_cut"],
                    "bottom_cut": site["bottom_cut"],
                    "top_cut_0based": site["top_cut_0based"],
                    "bottom_cut_0based": site["bottom_cut_0based"],
                    "top_cut_1based": site["top_cut"],
                    "bottom_cut_1based": site["bottom_cut"],
                    "cut_info": {"strand": site["strand"], "overhang": site["overhang"],
                                 "overhang_length": site["overhang_length"]},
                    "status": site["status"],
                    "wraps_origin": site["wraps_origin"],
                    "strand": site["strand"],
                })
            continue

        enzyme = get_restriction_enzyme(enzyme_name)
        if enzyme is None:
            continue

        for match in find_restriction_site_matches(
            sequence,
            enzyme_name,
            topology=topology,
        ):
            if match["status"] != "CONFIRMED" and not include_possible:
                continue

            position = match["position"]
            cut_info = calculate_cut_positions(
                position,
                enzyme,
            )

            # On a circular molecule a cut beyond the last base belongs at
            # the start (origin-spanning site).
            physical_cut = cut_info["top_cut"]
            if str(topology).lower() == "circular":
                physical_cut %= len(sequence)

            cuts.append({
                "position": physical_cut,
                "position_0based": physical_cut,
                "position_1based": physical_cut + 1,
                "enzyme": enzyme_name,
                "site_position": position + 1,
                "site_position_0based": position,
                "site_position_1based": position + 1,
                "top_cut": cut_info["top_cut_1based"],
                "bottom_cut": cut_info["bottom_cut_1based"],
                "top_cut_0based": cut_info["top_cut"],
                "bottom_cut_0based": cut_info["bottom_cut"],
                "top_cut_1based": cut_info["top_cut_1based"],
                "bottom_cut_1based": cut_info["bottom_cut_1based"],
                "cut_info": cut_info,
                "status": match["status"],
                "wraps_origin": match.get("wraps_origin", False),
            })

    cuts.sort(
        key=lambda item: (
            item["position_0based"],
            item["enzyme"],
        )
    )
    return cuts


def calculate_linear_fragments(sequence_length, cut_positions):
    """Return fragments from unique physical cut coordinates."""
    sequence_length = int(sequence_length)

    if sequence_length <= 0:
        return []

    positions = sorted({
        int(position)
        for position in cut_positions
        if 0 < int(position) < sequence_length
    })

    boundaries = [0] + positions + [sequence_length]
    fragments = []

    for index, (start, end) in enumerate(
        zip(boundaries, boundaries[1:]),
        start=1,
    ):
        if end <= start:
            continue

        fragments.append({
            "fragment": index,
            "start": start + 1,
            "end": end,
            "start_0based": start,
            "end_0based_exclusive": end,
            "length": end - start,
            "wraps_origin": False,
        })

    return fragments


def calculate_circular_fragments(sequence_length, cut_positions):
    """Return circular fragments from unique physical cut coordinates."""
    sequence_length = int(sequence_length)

    if sequence_length <= 0:
        return []

    positions = sorted({
        int(position) % sequence_length
        for position in cut_positions
    })

    if not positions:
        return [{
            "fragment": 1,
            "start": 1,
            "end": sequence_length,
            "start_0based": 0,
            "end_0based_exclusive": sequence_length,
            "length": sequence_length,
            "wraps_origin": False,
        }]

    if len(positions) == 1:
        cut = positions[0]
        return [{
            "fragment": 1,
            "start": cut + 1,
            "end": cut if cut > 0 else sequence_length,
            "start_0based": cut,
            "end_0based_exclusive": cut,
            "length": sequence_length,
            "wraps_origin": True,
        }]

    fragments = []

    for index, start in enumerate(positions, start=1):
        next_cut = (
            positions[index]
            if index < len(positions)
            else positions[0] + sequence_length
        )

        fragments.append({
            "fragment": index,
            "start": start + 1,
            "end": (
                next_cut
                if next_cut <= sequence_length
                else next_cut - sequence_length
            ),
            "start_0based": start,
            "end_0based_exclusive": next_cut % sequence_length,
            "length": next_cut - start,
            "wraps_origin": next_cut > sequence_length,
        })

    return fragments


def perform_digest(sequence, enzyme_names, topology="linear", include_possible=False):
    """Perform an in-silico digest with physical-boundary accounting.

    Ambiguity-dependent (POSSIBLE) sites are reported in "possible_cuts" but
    do not define fragments unless include_possible=True.
    """
    sequence = _require_valid_restriction_sequence(sequence)

    topology = str(topology).lower()
    if topology not in {"linear", "circular"}:
        raise ValueError(
            "Topology must be 'linear' or 'circular'."
        )

    if isinstance(enzyme_names, str):
        enzyme_names = [enzyme_names]

    enzyme_names = list(dict.fromkeys(enzyme_names))

    unknown = [
        name
        for name in enzyme_names
        if get_restriction_enzyme(name) is None
        and get_type_iis_enzyme(name) is None
    ]
    if unknown:
        raise ValueError(
            "Unknown restriction enzyme(s): "
            + ", ".join(unknown)
        )

    cuts = get_cut_positions(
        sequence,
        enzyme_names,
        include_possible=include_possible,
        topology=topology,
    )

    possible_cuts = [
        cut
        for cut in get_cut_positions(
            sequence,
            enzyme_names,
            include_possible=True,
            topology=topology,
        )
        if cut["status"] != "CONFIRMED"
    ]

    # The top-strand cut coordinate represents the physical fragment
    # boundary for standard bp-length reporting. Multiple enzymes at the
    # same coordinate are one physical boundary.
    if topology == "circular":
        unique_cut_positions = sorted({
            cut["position_0based"] % len(sequence)
            for cut in cuts
        })
        fragments = calculate_circular_fragments(
            len(sequence),
            unique_cut_positions,
        )
    else:
        unique_cut_positions = sorted({
            cut["position_0based"]
            for cut in cuts
            if 0 < cut["position_0based"] < len(sequence)
        })
        fragments = calculate_linear_fragments(
            len(sequence),
            unique_cut_positions,
        )

    for fragment in fragments:
        fragment["enzymes"] = []
        fragment["cut_events"] = []

    # Each fragment is bounded by cut coordinates. Store the enzyme events
    # touching either boundary for UI/reporting purposes. The cuts are
    # indexed by coordinate, so each fragment looks up its two boundaries
    # instead of scanning every cut (fragments x cuts: the full enzyme panel
    # on a 4.6 Mb genome did not finish in 15 minutes; docs/BENCHMARKS.md).
    # Events keep the order of the cut list, as before. (In a circular
    # digest a cut at coordinate 0 that bounds a fragment across the origin
    # is a boundary of that fragment like any other.)
    cuts_at = {}
    for index, cut in enumerate(cuts):
        cuts_at.setdefault(cut["position_0based"], []).append(index)

    for fragment in fragments:
        boundaries = {
            fragment["start_0based"],
            fragment["end_0based_exclusive"],
        }
        event_indices = sorted(
            index
            for boundary in boundaries
            for index in cuts_at.get(boundary, ())
        )
        events = [cuts[index] for index in event_indices]
        fragment["cut_events"] = [
            {
                "enzyme": cut["enzyme"],
                "position_1based": cut["position_1based"],
            }
            for cut in events
        ]
        fragment["enzymes"] = sorted({
            cut["enzyme"]
            for cut in events
        })

    return {
        "sequence_length": len(sequence),
        "topology": topology,
        "enzymes": enzyme_names,
        "cut_count": len(cuts),
        "unique_cut_count": len(unique_cut_positions),
        "cuts": cuts,
        "cut_positions": unique_cut_positions,
        "fragment_count": len(fragments),
        "fragments": fragments,
        "has_cut": bool(cuts),
        "possible_cuts": possible_cuts,
        "possible_cuts_included": bool(include_possible),
    }


def calculate_gel_bands(fragments):
    """Group equal-sized fragments into predicted bands."""
    bands = {}

    for fragment in fragments or []:
        length = int(fragment.get("length", 0))
        if length <= 0:
            continue

        band = bands.setdefault(
            length,
            {
                "size": length,
                "copies": 0,
                "fragments": [],
            },
        )

        band["copies"] += 1
        band["fragments"].append(
            fragment["fragment"]
        )

    result = sorted(
        bands.values(),
        key=lambda band: band["size"],
        reverse=True,
    )

    for band in result:
        band["relative_abundance"] = band["copies"]

    return result


def summarize_restriction_map(sequence, enzyme_names=None, topology="linear"):
    sites = create_restriction_map(
        sequence,
        enzyme_names,
        topology=topology,
    )

    summary = {}

    for site in sites:
        enzyme = site["enzyme"]
        entry = summary.setdefault(
            enzyme,
            {
                "site_count": 0,
                "confirmed_site_count": 0,
                "possible_site_count": 0,
                "positions": [],
            },
        )

        entry["site_count"] += 1
        entry["confirmed_site_count"] += (
            site["match_status"] == "CONFIRMED"
        )
        entry["possible_site_count"] += (
            site["match_status"] == "POSSIBLE"
        )
        entry["positions"].append(
            site["site_position"]
        )

    return summary


def find_unique_cutters(sequence, enzyme_names=None, topology="linear"):
    sequence = _require_valid_restriction_sequence(sequence)

    if enzyme_names is None:
        enzyme_names = list(RESTRICTION_ENZYMES.keys())

    return [
        enzyme_name
        for enzyme_name in enzyme_names
        if get_restriction_enzyme(enzyme_name) is not None
        and len(find_restriction_sites(sequence, enzyme_name, topology=topology)) == 1
        and analyze_enzyme(
            sequence,
            enzyme_name,
            topology=topology,
        )["confirmed_site_count"] == 1
    ]


def find_non_cutters(sequence, enzyme_names=None, topology="linear"):
    sequence = _require_valid_restriction_sequence(sequence)

    if enzyme_names is None:
        enzyme_names = list(RESTRICTION_ENZYMES.keys())

    return [
        enzyme_name
        for enzyme_name in enzyme_names
        if get_restriction_enzyme(enzyme_name) is not None
        and len(find_restriction_sites(sequence, enzyme_name, topology=topology)) == 0
    ]


def find_internal_sites(sequence, enzyme_names=None, topology="linear"):
    sequence = _require_valid_restriction_sequence(sequence)

    if enzyme_names is None:
        enzyme_names = list(RESTRICTION_ENZYMES.keys())

    warnings = []

    for enzyme_name in enzyme_names:
        if get_restriction_enzyme(enzyme_name) is None:
            continue

        analysis = analyze_enzyme(
            sequence,
            enzyme_name,
            topology=topology,
        )

        if analysis["site_count"] > 1:
            warnings.append({
                "enzyme": enzyme_name,
                "site_count": analysis["site_count"],
                "confirmed_site_count": analysis[
                    "confirmed_site_count"
                ],
                "possible_site_count": analysis[
                    "possible_site_count"
                ],
                "positions": [
                    site["site_position"]
                    for site in analysis["sites"]
                ],
                "status": "WARNING",
            })

    return warnings


# ============================================================
# REPORT
# ============================================================

def get_physical_cut_boundaries(digest_result):
    """Return sorted unique physical cut coordinates from a digest."""
    if not digest_result:
        return []

    if "cut_positions" in digest_result:
        return list(digest_result["cut_positions"])

    return sorted({
        int(cut["position_0based"])
        for cut in digest_result.get("cuts", [])
    })
