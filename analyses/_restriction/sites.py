"""Sequence validation, IUPAC recognition-site matching, cut positions and overhangs."""

from __future__ import annotations

from data.restriction_enzymes import get_restriction_enzyme, get_type_iis_enzyme
from functools import lru_cache


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


def find_type_iis_sites(sequence, enzyme_name, topology="linear"):
    """Sites of a Type IIS enzyme on both strands, with their cut coordinates (P2-7).

    The recognition sequence is asymmetric, so it is searched as given (strand
    "+") and as its reverse complement (strand "-"). The arithmetic uses cut
    boundaries (bases before the cut, top-strand coordinates). For a "+" site at p,
    top = p + cut_top and bottom = p + cut_bottom; for a "-" site the cuts are
    mirrored to the other side: top = p + L - cut_bottom, bottom = p + L - cut_top.
    The overhang is the top-strand sequence between the two cuts.

    Records follow the restriction-map convention: site_position, top_cut
    and bottom_cut are 1-based (for a cut, the first base after it); the
    *_0based keys hold the 0-based site start and the cut boundaries.

    On a linear molecule a site whose cuts fall outside the sequence does not
    cut ("cuts_inside" False); on a circular molecule the cuts wrap around the
    origin. As for every enzyme, a site with an IUPAC ambiguity inside the
    recognition sequence is POSSIBLE, not CONFIRMED.
    """
    enzyme = get_type_iis_enzyme(enzyme_name)
    if enzyme is None:
        raise ValueError(f"Not a Type IIS enzyme in the Synapse database: {enzyme_name}")
    sequence = _require_valid_restriction_sequence(sequence)
    circular = str(topology).lower() == "circular"
    length = len(sequence)
    site = _normalize_restriction_sequence(enzyme["recognition_sequence"])
    size = len(site)
    cut_top, cut_bottom = int(enzyme["cut_top"]), int(enzyme["cut_bottom"])

    sites = []
    orientations = [("+", site)]
    if reverse_complement(site) != site:
        orientations.append(("-", reverse_complement(site)))
    for strand, pattern in orientations:
        for match in _scan_recognition_sites(sequence, pattern, circular):
            start = match["position"]
            if strand == "+":
                top, bottom = start + cut_top, start + cut_bottom
            else:
                top, bottom = start + size - cut_bottom, start + size - cut_top
            if circular:
                inside = length >= size + max(cut_top, cut_bottom)
                left = min(top, bottom)
                overhang = "".join(sequence[(left + i) % length] for i in range(abs(bottom - top)))
                top_physical, bottom_physical = top % length, bottom % length
            else:
                inside = 0 < top < length and 0 < bottom < length
                overhang = sequence[min(top, bottom):max(top, bottom)] if inside else ""
                top_physical, bottom_physical = top, bottom
            sites.append({
                "enzyme": enzyme_name,
                "strand": strand,
                "site_position": start + 1,
                "site_position_0based": start,
                "site_position_1based": start + 1,
                "site_sequence": match["sequence"],
                "status": match["status"],
                "wraps_origin": match["wraps_origin"] or (circular and not 0 <= min(top, bottom) <= max(top, bottom) <= length),
                "cuts_inside": inside,
                "top_cut": top_physical + 1,
                "bottom_cut": bottom_physical + 1,
                "top_cut_0based": top_physical,
                "bottom_cut_0based": bottom_physical,
                "top_cut_1based": top_physical + 1,
                "bottom_cut_1based": bottom_physical + 1,
                "overhang": overhang,
                "overhang_length": abs(bottom - top),
            })
    sites.sort(key=lambda item: (item["top_cut_0based"], item["strand"]))
    return sites


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
