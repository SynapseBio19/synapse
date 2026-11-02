from functools import lru_cache
from itertools import combinations

from data.restriction_enzymes import (
    RESTRICTION_ENZYMES,
    get_restriction_enzyme,
)
from ui.style import major_rule, minor_rule


# ============================================================
# IUPAC DNA
# ============================================================

IUPAC_BASES = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

# ============================================================
# SYNAPSE RESTRICTION ENGINE - FINAL CONSISTENCY / QUALITY LAYER v3
# ============================================================
#
# Production restriction engine.
#
# Design rules:
# - Coordinates are 0-based internally.
# - Human-readable coordinates are 1-based.
# - IUPAC-dependent sites are separated from confirmed sites.
# - A physical digest boundary is counted once even if multiple enzymes
#   cut at the same coordinate.
# - Staggered cuts retain both strand coordinates and the overhang.
# - Public functions preserve the Synapse pipeline API.
# ============================================================


def _normalize_restriction_sequence(sequence):
    """Normalize DNA without silently deleting non-whitespace characters."""
    if not isinstance(sequence, str):
        return ""
    return "".join(sequence.upper().split())


def _restriction_sequence_status(sequence):
    normalized = _normalize_restriction_sequence(sequence)

    invalid_positions = [
        index + 1
        for index, base in enumerate(normalized)
        if base not in IUPAC_BASES
    ]

    return {
        "sequence": normalized,
        "valid": bool(normalized) and not invalid_positions,
        "length": len(normalized),
        "invalid_positions": invalid_positions,
        "invalid_symbols": sorted({
            normalized[position - 1]
            for position in invalid_positions
        }),
        "ambiguous_positions": [
            index + 1
            for index, base in enumerate(normalized)
            if base not in {"A", "C", "G", "T"}
            and base in IUPAC_BASES
        ],
    }


def validate_restriction_sequence(sequence):
    """Return True only for non-empty DNA using the supported IUPAC alphabet."""
    return _restriction_sequence_status(sequence)["valid"]


def _require_valid_restriction_sequence(sequence):
    status = _restriction_sequence_status(sequence)
    if not status["valid"]:
        if not status["sequence"]:
            raise ValueError("DNA sequence is empty.")
        invalid = ", ".join(
            f"{symbol}@{position}"
            for position, symbol in (
                (p, status["sequence"][p - 1])
                for p in status["invalid_positions"]
            )
        )
        raise ValueError(
            f"Invalid DNA sequence. Unsupported symbols: {invalid}."
        )
    return status["sequence"]


def bases_match(sequence_base, recognition_base):
    """Return whether two IUPAC symbols have at least one possible overlap."""
    sequence_options = IUPAC_BASES.get(
        str(sequence_base).upper()
    )
    recognition_options = IUPAC_BASES.get(
        str(recognition_base).upper()
    )

    if sequence_options is None or recognition_options is None:
        return False

    return bool(sequence_options & recognition_options)


def recognition_site_matches(sequence, position, recognition_sequence):
    """Check one recognition site using IUPAC-compatible matching."""
    sequence = _normalize_restriction_sequence(sequence)
    recognition_sequence = _normalize_restriction_sequence(
        recognition_sequence
    )

    if position < 0 or not recognition_sequence:
        return False

    end = position + len(recognition_sequence)
    if end > len(sequence):
        return False

    return all(
        bases_match(a, b)
        for a, b in zip(
            sequence[position:end],
            recognition_sequence,
        )
    )


def find_restriction_site_matches(sequence, enzyme_name, topology="linear"):
    """Return all sequence-compatible sites with CONFIRMED/POSSIBLE status.

    topology="circular" also finds recognition sites that span the origin
    (end -> start junction); such sites have "wraps_origin": True and a
    0-based "position" near the end of the sequence. Sites are searched on
    the forward strand, which is complete for the palindromic enzymes in the
    Synapse database.
    """
    enzyme = get_restriction_enzyme(enzyme_name)

    if enzyme is None:
        raise ValueError(f"Unknown restriction enzyme: {enzyme_name}")

    return find_recognition_site_matches(
        sequence,
        enzyme.get("recognition_sequence", ""),
        topology=topology,
    )


def find_recognition_site_matches(sequence, recognition_sequence, topology="linear"):
    """The single Synapse recognition-site matcher (restriction and cloning).

    Same result format as find_restriction_site_matches(). A window counts
    as a site when every base is IUPAC-compatible with the recognition
    sequence; it is CONFIRMED only when the window has no ambiguous symbol,
    otherwise POSSIBLE. Raises ValueError for an empty or non-IUPAC sequence.
    """
    sequence = _require_valid_restriction_sequence(sequence)
    circular = str(topology).lower() == "circular"

    recognition = _normalize_restriction_sequence(recognition_sequence)

    if not recognition:
        return []

    # Fresh dicts per call: the scan is cached, and a caller that edited a
    # returned match must not reach into the cache.
    return [dict(match) for match in _scan_recognition_sites(sequence, recognition, circular)]


@lru_cache(maxsize=256)
def _scan_recognition_sites(sequence, recognition, circular):
    """Scan one sequence for one recognition sequence; cached (audit S18).

    The scan depends only on these three values, but the cloning planner
    re-ran it for every candidate strategy: a 1.3 kb vector with 8 enzymes
    produced 5,002 calls for 16 distinct combinations. Caching removes the
    repetition without changing any result. The bound keeps a planner run's
    worth of scans and nothing more.
    """
    matches = []

    search_space = sequence
    last_start = len(sequence) - len(recognition)
    if circular and len(sequence) >= len(recognition):
        search_space = sequence + sequence[:len(recognition) - 1]
        last_start = len(sequence) - 1

    for position in range(
        last_start + 1
    ):
        fragment = search_space[
            position:position + len(recognition)
        ]

        if not all(
            bases_match(a, b)
            for a, b in zip(fragment, recognition)
        ):
            continue

        ambiguous_offsets = [
            index + 1
            for index, base in enumerate(fragment)
            if base not in {"A", "C", "G", "T"}
        ]

        matches.append({
            "position": position,
            "position_1based": position + 1,
            "status": (
                "POSSIBLE"
                if ambiguous_offsets
                else "CONFIRMED"
            ),
            "ambiguous_offsets": ambiguous_offsets,
            "sequence": fragment,
            "recognition_sequence": recognition,
            "wraps_origin": position + len(recognition) > len(sequence),
        })

    return tuple(matches)


def find_restriction_sites(sequence, enzyme_name, topology="linear"):
    """Legacy-compatible site list: all confirmed or ambiguity-dependent matches."""
    return [
        match["position"]
        for match in find_restriction_site_matches(
            sequence,
            enzyme_name,
            topology=topology,
        )
    ]


def calculate_cut_positions(site_position, enzyme):
    """Calculate validated 0-based strand cut coordinates."""
    try:
        site_position = int(site_position)
    except (TypeError, ValueError) as exc:
        raise ValueError("Site position must be an integer.") from exc

    if site_position < 0:
        raise ValueError("Site position cannot be negative.")

    recognition = _normalize_restriction_sequence(
        enzyme.get("recognition_sequence", "")
    )

    if not recognition:
        raise ValueError("Enzyme has no recognition sequence.")

    try:
        top_offset = int(enzyme["cut_top"])
        bottom_offset = int(enzyme["cut_bottom"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            "Enzyme database entry is missing valid cut_top/cut_bottom values."
        ) from exc

    if not 0 <= top_offset <= len(recognition):
        raise ValueError(
            f"Invalid cut_top={top_offset} for {enzyme.get('recognition_sequence')}."
        )

    if not 0 <= bottom_offset <= len(recognition):
        raise ValueError(
            f"Invalid cut_bottom={bottom_offset} for {enzyme.get('recognition_sequence')}."
        )

    top_cut = site_position + top_offset
    bottom_cut = site_position + bottom_offset

    return {
        "top_cut": top_cut,
        "bottom_cut": bottom_cut,
        "top_cut_0based": top_cut,
        "bottom_cut_0based": bottom_cut,
        "top_cut_1based": top_cut + 1,
        "bottom_cut_1based": bottom_cut + 1,
        "recognition_start": site_position,
        "recognition_end": site_position + len(recognition) - 1,
        "recognition_start_1based": site_position + 1,
        "recognition_end_1based": site_position + len(recognition),
        "overhang_length": abs(bottom_cut - top_cut),
    }


def reverse_complement(sequence):
    """Reverse-complement standard and IUPAC DNA symbols."""
    complement = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B",
        "D": "H", "H": "D", "N": "N",
    }

    normalized = _normalize_restriction_sequence(sequence)

    invalid = [
        base for base in normalized
        if base not in complement
    ]
    if invalid:
        raise ValueError(
            f"Unsupported DNA symbols in reverse complement: "
            f"{', '.join(sorted(set(invalid)))}"
        )

    return "".join(
        complement[base]
        for base in reversed(normalized)
    )


def get_overhang_sequence(sequence, site_position, enzyme):
    """Extract the staggered region in forward-strand coordinates."""
    sequence = _require_valid_restriction_sequence(sequence)
    cuts = calculate_cut_positions(site_position, enzyme)

    top_cut = cuts["top_cut"]
    bottom_cut = cuts["bottom_cut"]

    if (
        enzyme.get("is_blunt", False)
        or top_cut == bottom_cut
    ):
        return ""

    start = min(top_cut, bottom_cut)
    end = max(top_cut, bottom_cut)

    if start < 0 or end > len(sequence):
        raise ValueError(
            "Restriction cut coordinates fall outside the supplied sequence."
        )

    return sequence[start:end]


def _site_match_for_position(sequence, site_position, recognition):
    return recognition_site_matches(
        sequence,
        site_position,
        recognition,
    )


def analyze_restriction_site(sequence, enzyme_name, site_position, topology="linear"):
    """Complete coordinate-safe analysis of one restriction site.

    topology="circular" also accepts a site that spans the origin; its end
    and cut coordinates are reported modulo the sequence length (the
    convention of perform_digest) and its overhang is read across the
    origin. Linear results are unchanged.
    """
    sequence = _require_valid_restriction_sequence(sequence)
    enzyme = get_restriction_enzyme(enzyme_name)

    if enzyme is None:
        raise ValueError(
            f"Unknown restriction enzyme: {enzyme_name}"
        )

    recognition = _normalize_restriction_sequence(
        enzyme.get("recognition_sequence", "")
    )

    circular = str(topology).lower() == "circular"
    length = len(sequence)
    search_space = sequence
    if circular and length >= len(recognition):
        search_space = sequence + sequence[:len(recognition) - 1]
        if site_position >= length:
            raise ValueError(
                "Site positions on a circular sequence must lie within the sequence."
            )

    if not _site_match_for_position(
        search_space,
        site_position,
        recognition,
    ):
        raise ValueError(
            "The specified position does not contain a compatible recognition site."
        )

    matches = find_restriction_site_matches(
        sequence,
        enzyme_name,
        topology=topology,
    )
    match = next(
        (
            item
            for item in matches
            if item["position"] == site_position
        ),
        None,
    )

    cuts = calculate_cut_positions(
        site_position,
        enzyme,
    )

    overhang = get_overhang_sequence(
        search_space,
        site_position,
        enzyme,
    )

    is_blunt = bool(
        enzyme.get("is_blunt", False)
        or cuts["top_cut"] == cuts["bottom_cut"]
    )

    overhang_type = str(
        enzyme.get(
            "overhang_type",
            "Blunt" if is_blunt else "Unknown",
        )
    )

    site_end = site_position + len(recognition)
    top_cut = cuts["top_cut"]
    bottom_cut = cuts["bottom_cut"]
    if circular:
        site_end = (site_end - 1) % length + 1
        top_cut %= length
        bottom_cut %= length

    return {
        "enzyme": enzyme_name,
        "recognition_sequence": recognition,
        "site_position": site_position + 1,
        "site_position_0based": site_position,
        "site_end": site_end,
        "site_end_1based": site_end,
        "top_cut": top_cut + 1,
        "top_cut_0based": top_cut,
        "bottom_cut": bottom_cut + 1,
        "bottom_cut_0based": bottom_cut,
        "overhang_type": overhang_type,
        "overhang_length": cuts["overhang_length"],
        "overhang_sequence": overhang,
        "is_blunt": is_blunt,
        "cut_type": "Blunt" if is_blunt else "Sticky",
        "match_status": (
            match["status"] if match else "UNKNOWN"
        ),
        "ambiguity_dependent": bool(
            match and match["status"] == "POSSIBLE"
        ),
    }


def analyze_enzyme(sequence, enzyme_name, topology="linear"):
    """Analyze every compatible site of one enzyme."""
    sequence = _require_valid_restriction_sequence(sequence)
    enzyme = get_restriction_enzyme(enzyme_name)

    if enzyme is None:
        raise ValueError(
            f"Unknown restriction enzyme: {enzyme_name}"
        )

    matches = find_restriction_site_matches(
        sequence,
        enzyme_name,
        topology=topology,
    )

    sites = [
        analyze_restriction_site(
            sequence,
            enzyme_name,
            match["position"],
            topology=topology,
        )
        for match in matches
    ]

    confirmed_count = sum(
        site["match_status"] == "CONFIRMED"
        for site in sites
    )
    possible_count = sum(
        site["match_status"] == "POSSIBLE"
        for site in sites
    )

    return {
        "enzyme": enzyme_name,
        "recognition_sequence": _normalize_restriction_sequence(
            enzyme.get("recognition_sequence", "")
        ),
        "site_count": len(sites),
        "confirmed_site_count": confirmed_count,
        "possible_site_count": possible_count,
        "sites": sites,
        "cut_type": (
            "Blunt"
            if enzyme.get("is_blunt", False)
            else "Sticky"
        ),
        "overhang_type": str(
            enzyme.get("overhang_type", "Unknown")
        ),
        "overhang_length": int(
            enzyme.get(
                "overhang_length",
                0,
            )
        ),
        "database_valid": True,
    }


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
                "site_position": position,
                "site_position_1based": position + 1,
                "top_cut": cut_info["top_cut"],
                "bottom_cut": cut_info["bottom_cut"],
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
    # touching either boundary for UI/reporting purposes.
    if topology == "linear":
        for fragment in fragments:
            boundaries = {
                fragment["start_0based"],
                fragment["end_0based_exclusive"],
            }
            events = [
                cut
                for cut in cuts
                if cut["position_0based"] in boundaries
            ]
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
    else:
        for fragment in fragments:
            boundaries = {
                fragment["start_0based"],
                fragment["end_0based_exclusive"],
            }
            events = [
                cut
                for cut in cuts
                if (
                    cut["position_0based"] in boundaries
                    or (
                        fragment["wraps_origin"]
                        and cut["position_0based"] == 0
                        and 0 in boundaries
                    )
                )
            ]
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

    results["public_api"] = all(
        name in globals()
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


def print_restriction_report(analysis):
    """Print a human-readable report without changing analysis data."""
    if not analysis:
        print("No restriction analysis available.")
        return

    print()
    print(major_rule())
    print("SYNAPSE — RESTRICTION ANALYSIS")
    print(major_rule())

    print("\nSequence")
    print(minor_rule())
    print(f"Length:              {analysis['sequence_length']} bp")
    print(f"Topology:            {analysis['topology'].capitalize()}")
    print(f"Enzymes analyzed:    {len(analysis['enzymes'])}")

    status = analysis.get("sequence_status", {})
    ambiguous = status.get("ambiguous_positions", [])
    if ambiguous:
        print(
            f"IUPAC ambiguous bases: {len(ambiguous)} "
            f"(positions: {', '.join(map(str, ambiguous[:20]))}"
            f"{'…' if len(ambiguous) > 20 else ''})"
        )

    print("\nRestriction Enzymes")
    print(minor_rule())
    print(
        f"{'Enzyme':<12} {'Recognition':<14} "
        f"{'Sites':>6} {'Confirmed':>10} {'Possible':>9} {'Type':>8}"
    )

    for enzyme_name in analysis["enzymes"]:
        result = analysis["enzyme_results"][enzyme_name]
        print(
            f"{enzyme_name:<12} "
            f"{result['recognition_sequence']:<14} "
            f"{result['site_count']:>6} "
            f"{result['confirmed_site_count']:>10} "
            f"{result['possible_site_count']:>9} "
            f"{result['cut_type']:>8}"
        )

    print("\nRestriction Map")
    print(minor_rule())
    sites = analysis.get("restriction_map", [])

    if not sites:
        print("No restriction sites detected.")
    else:
        for site in sites:
            overhang = site.get("overhang_sequence") or "-"
            status_text = site.get("match_status", "UNKNOWN")
            print(
                f"{site['enzyme']:<12} "
                f"site {site['site_position']:>7} "
                f"| cut {site['top_cut']:>7}/{site['bottom_cut']:<7} "
                f"| {site['cut_type']:<6} "
                f"| {site['overhang_type']:<12} "
                f"| {overhang:<8} "
                f"| {status_text}"
            )

    digest = analysis["digest"]

    print("\nDigest")
    print(minor_rule())
    print(f"Enzyme cut events:       {digest['cut_count']}")
    print(f"Physical cut boundaries: {digest['unique_cut_count']}")
    print(f"Expected fragments:      {digest['fragment_count']}")

    if digest["fragments"]:
        for fragment in digest["fragments"]:
            enzymes = ", ".join(fragment.get("enzymes", [])) or "-"
            print(
                f"Fragment {fragment['fragment']:>3}: "
                f"{fragment['length']:>8} bp"
                f" | Cut by: {enzymes}"
            )

    print("\nPredicted Gel Bands")
    print(minor_rule())
    bands = analysis.get("gel_bands", [])
    if not bands:
        print("No bands available.")
    else:
        for band in bands:
            print(
                f"{band['size']:>8} bp | "
                f"{band['copies']} fragment"
                f"{'s' if band['copies'] != 1 else ''}"
            )

    print("\nUnique Cutters")
    print(minor_rule())
    print(", ".join(analysis["unique_cutters"]) or "None detected.")

    print("\nNon-Cutting Enzymes")
    print(minor_rule())
    print(", ".join(analysis["non_cutters"]) or "None detected.")

    print("\nMultiple-Site Warnings")
    print(minor_rule())
    warnings = analysis.get("internal_sites", [])
    if not warnings:
        print("No multiple-site warnings detected.")
    else:
        for warning in warnings:
            positions = ", ".join(map(str, warning["positions"]))
            print(
                f"WARNING: {warning['enzyme']} has "
                f"{warning['site_count']} sites at {positions}"
            )

    print("\nRecommended Enzyme Pairs")
    print(minor_rule())
    recommendations = analysis.get("recommended_pairs", [])
    if not recommendations:
        print("No suitable enzyme pairs detected.")
    else:
        for index, recommendation in enumerate(recommendations, 1):
            print(
                f"{index}. {recommendation['enzyme_a']} + "
                f"{recommendation['enzyme_b']} "
                f"({recommendation['score']}/{recommendation['maximum']})"
            )
            for reason in recommendation["reasons"]:
                print(f"   + {reason}")
            for warning in recommendation["warnings"]:
                print(f"   ! {warning}")

    print("\nInterpretation")
    print(minor_rule())
    if digest["cut_count"] == 0:
        print("No selected restriction enzyme cuts the sequence.")
    elif analysis["topology"] == "circular":
        print(
            "Circular DNA: the number of fragments equals the number "
            "of distinct physical cut boundaries."
        )
    else:
        print(
            "Linear DNA: n distinct internal cut boundaries produce "
            "n + 1 fragments."
        )

    if analysis.get("site_statistics"):
        stats = analysis["site_statistics"]
        print(
            f"Site events: {stats['total_site_events']} | "
            f"confirmed: {stats['confirmed_site_events']} | "
            f"possible: {stats['possible_site_events']}"
        )

    print(
        "\nNote: This is sequence-level in-silico analysis. "
        "Experimental digestion depends on enzyme-specific conditions, "
        "DNA context, methylation and other laboratory factors."
    )
    print(major_rule())


if __name__ == "__main__":
    self_test = run_restriction_self_test()

    print(major_rule())
    print("SYNAPSE RESTRICTION ENGINE SELF-TEST")
    print(major_rule())

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


def _synapse_ui_pause():
    input("\nPress Enter to continue...")


def _synapse_ui_select_enzymes(
    allow_all=True,
    minimum=1,
    maximum=None,
):
    enzyme_names = list(RESTRICTION_ENZYMES.keys())

    print("\nAvailable restriction enzymes:")
    print(minor_rule())

    for number, enzyme_name in enumerate(enzyme_names, start=1):
        enzyme = RESTRICTION_ENZYMES[enzyme_name]
        print(
            f"{number:>2}. {enzyme_name:<10} "
            f"{enzyme['recognition_sequence']}"
        )

    print("\nEnter enzyme names separated by commas.")
    print("Example: EcoRI,HindIII")
    if allow_all:
        print("Enter ALL to analyze all enzymes.")
    print("Enter 0 to go back.")

    while True:
        raw = input("\nEnzymes: ").strip()

        if raw == "0":
            return None

        if allow_all and raw.upper() == "ALL":
            return enzyme_names

        selected = list(dict.fromkeys(
            name.strip()
            for name in raw.split(",")
            if name.strip()
        ))

        invalid = [
            name for name in selected
            if name not in RESTRICTION_ENZYMES
        ]

        if invalid:
            print("Unknown enzyme(s): " + ", ".join(invalid))
            continue

        if len(selected) < minimum:
            print(f"Please select at least {minimum} enzyme(s).")
            continue

        if maximum is not None and len(selected) > maximum:
            print(f"Please select no more than {maximum} enzyme(s).")
            continue

        return selected


def _synapse_ui_select_topology(default="linear"):
    print("\nDNA Topology")
    print(minor_rule())
    print("[1] Linear")
    print("[2] Circular")
    print(f"[Enter] Keep current ({str(default).capitalize()})")
    print("[0] Back")

    while True:
        choice = input("\nSelect topology: ").strip()

        if choice == "0":
            return None
        if choice == "":
            return str(default).lower()
        if choice == "1":
            return "linear"
        if choice == "2":
            return "circular"

        print("Invalid topology selection.")


def _restriction_overview_ui(dna, pipeline=None, topology="linear"):
    print()
    print(major_rule())
    print("RESTRICTION OVERVIEW")
    print(major_rule())

    enzymes = _synapse_ui_select_enzymes()
    if enzymes is None:
        return None

    selected_topology = _synapse_ui_select_topology(topology)
    if selected_topology is None:
        return None

    try:
        analysis = analyze_restriction(
            dna,
            enzymes,
            selected_topology,
        )

        if pipeline is not None:
            try:
                pipeline.context.store_restriction_analysis(analysis)
            except AttributeError:
                pass

        print_restriction_report(analysis)
        _synapse_ui_pause()
        return analysis

    except (ValueError, TypeError, KeyError) as error:
        print(f"\nRestriction analysis error: {error}")
        _synapse_ui_pause()
        return None


def _restriction_map_ui(dna, topology="linear"):
    print()
    print(major_rule())
    print("RESTRICTION MAP")
    print(major_rule())

    enzymes = _synapse_ui_select_enzymes()
    if enzymes is None:
        return None

    selected_topology = _synapse_ui_select_topology(topology)
    if selected_topology is None:
        return None

    try:
        sites = get_restriction_map(
            dna,
            enzyme_names=enzymes,
            topology=selected_topology,
        )

        print()
        print(minor_rule())

        if not sites:
            print("No restriction sites detected.")
        else:
            for site in sites:
                position = site.get(
                    "site_position",
                    site.get("position", "-"),
                )
                print(
                    f"{site.get('enzyme', '-'):<10} "
                    f"site {position:>6} | "
                    f"cut {site.get('top_cut', '-')}/"
                    f"{site.get('bottom_cut', '-')} | "
                    f"{site.get('cut_type', '-'):<6} | "
                    f"{site.get('overhang_sequence') or '-'}"
                )

        _synapse_ui_pause()
        return sites

    except (ValueError, TypeError, KeyError) as error:
        print(f"\nRestriction map error: {error}")
        _synapse_ui_pause()
        return None


def _restriction_digest_gel_ui(dna, topology="linear"):
    print()
    print(major_rule())
    print("DIGEST & GEL ANALYSIS")
    print(major_rule())

    enzymes = _synapse_ui_select_enzymes(
        allow_all=False,
        minimum=1,
    )
    if enzymes is None:
        return None

    selected_topology = _synapse_ui_select_topology(topology)
    if selected_topology is None:
        return None

    try:
        digest = get_digest_prediction(
            dna,
            enzymes,
            topology=selected_topology,
        )
        bands = get_gel_band_prediction(digest)

        print()
        print(f"Topology:       {selected_topology.capitalize()}")
        print(f"Enzymes:        {', '.join(enzymes)}")
        print(f"Cut events:     {digest['cut_count']}")
        print(f"Cut boundaries: {digest['unique_cut_count']}")
        print(f"Fragments:      {digest['fragment_count']}")

        possible = digest.get("possible_cuts", [])
        if possible:
            print(
                f"\nNote: {len(possible)} ambiguity-dependent site(s) "
                "(IUPAC symbols such as N inside the recognition "
                "sequence) were NOT counted as cuts."
            )

        for enzyme_name in enzymes:
            enzyme = get_restriction_enzyme(enzyme_name) or {}
            note = enzyme.get("methylation_note")
            if note:
                print(f"Methylation note ({enzyme_name}): {note}")

        print("\nFragments")
        print(minor_rule())
        for fragment in digest.get("fragments", []):
            print(
                f"Fragment {fragment['fragment']:>2}: "
                f"{fragment['length']:>8} bp"
            )

        print("\nPredicted Gel Bands")
        print(minor_rule())
        if not bands:
            print("No bands available.")
        else:
            for band in bands:
                copies = band["copies"]
                print(
                    f"{band['size']:>8} bp | "
                    f"{copies} fragment"
                    f"{'s' if copies != 1 else ''}"
                )

        _synapse_ui_pause()
        return {"digest": digest, "gel_bands": bands}

    except (ValueError, TypeError, KeyError) as error:
        print(f"\nDigest analysis error: {error}")
        _synapse_ui_pause()
        return None


def _restriction_complete_report_ui(
    dna,
    pipeline=None,
    topology="linear",
):
    print()
    print(major_rule())
    print("COMPLETE RESTRICTION REPORT")
    print(major_rule())

    enzymes = _synapse_ui_select_enzymes()
    if enzymes is None:
        return None

    selected_topology = _synapse_ui_select_topology(topology)
    if selected_topology is None:
        return None

    try:
        analysis = analyze_restriction(
            dna,
            enzymes,
            selected_topology,
        )

        if pipeline is not None:
            try:
                pipeline.context.store_restriction_analysis(analysis)
            except AttributeError:
                pass

        print_restriction_report(analysis)
        _synapse_ui_pause()
        return analysis

    except (ValueError, TypeError, KeyError) as error:
        print(f"\nRestriction analysis error: {error}")
        _synapse_ui_pause()
        return None


def run_restriction_analysis_ui(
    dna,
    pipeline=None,
    header="Sequence",
    topology="linear",
):
    """Small user-facing Restriction menu."""

    while True:
        print()
        print(major_rule())
        print("RESTRICTION ANALYSIS")
        print(major_rule())

        print()
        print("[1] Restriction Overview")
        print("[2] Restriction Map")
        print("[3] Digest & Gel Analysis")
        print("[4] Complete Restriction Report")
        print()
        print("[0] Back")
        print(minor_rule())

        choice = input("\nSelection: ").strip()

        if choice == "0":
            return None

        if choice == "1":
            _restriction_overview_ui(
                dna,
                pipeline=pipeline,
                topology=topology,
            )
        elif choice == "2":
            _restriction_map_ui(dna, topology=topology)
        elif choice == "3":
            _restriction_digest_gel_ui(
                dna,
                topology=topology,
            )
        elif choice == "4":
            _restriction_complete_report_ui(
                dna,
                pipeline=pipeline,
                topology=topology,
            )
        else:
            print("\nInvalid selection.")


def run_restriction_cloning_ui(
    dna,
    header,
    sequences,
    pipeline=None,
    topology="linear",
):
    """
    Top-level module 5 UI.

    First choose Restriction or Cloning.
    Cloning remains owned by the cloning module.
    """

    try:
        from analyses.cloning import run_cloning_analysis_ui
    except ImportError:
        run_cloning_analysis_ui = None

    while True:
        print()
        print(major_rule())
        print("RESTRICTION & CLONING")
        print(major_rule())

        print()
        print("[1] Restriction Analysis")
        print("[2] Cloning Analysis")
        print()
        print("[0] Back")
        print(minor_rule())

        choice = input("\nSelection: ").strip()

        if choice == "0":
            return None

        if choice == "1":
            run_restriction_analysis_ui(
                dna=dna,
                pipeline=pipeline,
                header=header,
                topology=topology,
            )
            continue

        if choice == "2":
            if run_cloning_analysis_ui is None:
                print("\nCloning module could not be imported.")
                _synapse_ui_pause()
                continue

            run_cloning_analysis_ui(
                vector_sequence=dna,
                vector_header=header,
                sequences=sequences,
                pipeline=pipeline,
                vector_topology=topology,
            )
            continue

        print("\nInvalid selection.")
