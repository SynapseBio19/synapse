"""Sequences, enzyme data, overhangs, recognition sites and restriction maps."""

from __future__ import annotations

from analyses.restriction import find_recognition_site_matches
from data.restriction_enzymes import RESTRICTION_ENZYMES

from analyses._cloning.cache import _cached


# ============================================================
# IUPAC / DNA HELPERS
# ============================================================

IUPAC_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def normalize_sequence(sequence):
    if sequence is None:
        return ""

    return "".join(
        str(sequence).upper().split()
    )


def _normalize_topology(topology):
    value = str(topology).lower()
    if value not in {"linear", "circular"}:
        raise ValueError(
            "Topology must be 'linear' or 'circular'."
        )
    return value


def reverse_complement(sequence):
    sequence = normalize_sequence(sequence)

    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
        "N": "N",
    }

    return "".join(
        complement.get(base, "N")
        for base in reversed(sequence)
    )


# ============================================================
# ENZYME DATABASE ACCESS
# ============================================================

def get_enzyme_data(enzyme_name):
    if enzyme_name not in RESTRICTION_ENZYMES:
        raise ValueError(
            f"Unknown restriction enzyme: {enzyme_name}"
        )

    return RESTRICTION_ENZYMES[
        enzyme_name
    ]


def get_recognition_sequence(enzyme_name):
    enzyme = get_enzyme_data(
        enzyme_name
    )

    return str(
        enzyme.get(
            "recognition_sequence",
            "",
        )
    ).upper()


def get_cut_type(enzyme_name):
    enzyme = get_enzyme_data(
        enzyme_name
    )

    if enzyme.get(
        "is_blunt",
        False,
    ):
        return "Blunt"

    return str(
        enzyme.get(
            "cut_type",
            enzyme.get(
                "type",
                enzyme.get(
                    "overhang_type",
                    "Unknown",
                ),
            ),
        )
    )


def get_overhang(enzyme_name):
    enzyme = get_enzyme_data(
        enzyme_name
    )

    value = enzyme.get(
        "overhang",
        enzyme.get(
            "overhang_sequence",
            None,
        ),
    )

    if value is not None:
        return str(value).upper()

    # Derive the single-stranded overhang from the database cut positions.
    # For the palindromic Type II enzymes in the Synapse database the
    # overhang, read 5'->3' on the top strand, is the part of the
    # recognition site between the two strand cuts (G^AATT_C -> AATT,
    # CA^TA_TG -> TA, CTGCA^G with bottom cut at 1 -> TGCA).
    # This replaces a hand-written fallback table that contained errors
    # (NdeI "CATATG", ClaI "TCGA"); see docs/CHANGES.md, audit P0-4.
    recognition = get_recognition_sequence(enzyme_name)
    try:
        top = int(enzyme["cut_top"])
        bottom = int(enzyme["cut_bottom"])
    except (KeyError, TypeError, ValueError):
        return ""

    if top == bottom:
        return ""

    start, end = min(top, bottom), max(top, bottom)
    if start < 0 or end > len(recognition):
        return ""

    return recognition[start:end]


def get_overhang_polarity(enzyme_name):
    """Return '5_prime', '3_prime' or 'blunt' for an enzyme's DNA ends.

    The polarity is part of end chemistry: a 5' AGCT overhang (HindIII)
    cannot anneal to a 3' AGCT overhang (SacI) although the letters agree.
    """
    enzyme = get_enzyme_data(
        enzyme_name
    )

    if enzyme.get("is_blunt", False):
        return "blunt"

    try:
        top = int(enzyme["cut_top"])
        bottom = int(enzyme["cut_bottom"])
    except (KeyError, TypeError, ValueError):
        declared = str(enzyme.get("overhang_type", "")).lower()
        if declared in {"5_prime", "3_prime", "blunt"}:
            return declared
        return "unknown"

    if top == bottom:
        return "blunt"

    return "5_prime" if top < bottom else "3_prime"


def is_blunt_enzyme(enzyme_name):
    enzyme = get_enzyme_data(
        enzyme_name
    )

    if enzyme.get(
        "is_blunt",
        False,
    ):
        return True

    return "blunt" in get_cut_type(
        enzyme_name
    ).lower()


# ============================================================
# RESTRICTION SITE ANALYSIS USED BY CLONING
# ============================================================

def find_restriction_sites(
    sequence,
    recognition_sequence,
    include_possible=False,
    topology="linear",
):
    """1-based start positions of recognition sites.

    Only CONFIRMED sites (no IUPAC ambiguity inside the site) are returned
    unless include_possible=True; a run of N is not a site. On a circular
    sequence, sites spanning the origin are included (start near the end).
    Raises ValueError for an empty or non-IUPAC sequence.
    """
    return [
        match["position_1based"]
        for match in find_recognition_site_matches(
            normalize_sequence(sequence),
            recognition_sequence,
            topology=_normalize_topology(topology),
        )
        if include_possible or match["status"] == "CONFIRMED"
    ]


def analyze_enzyme(
    sequence,
    enzyme_name,
    topology="linear",
):
    return _cached(
        ("enzyme", sequence, enzyme_name, topology),
        lambda: _analyze_enzyme(sequence, enzyme_name, topology),
    )


def _analyze_enzyme(
    sequence,
    enzyme_name,
    topology,
):
    sequence = normalize_sequence(
        sequence
    )
    topology = _normalize_topology(
        topology
    )

    recognition_sequence = (
        get_recognition_sequence(
            enzyme_name
        )
    )

    matches = find_recognition_site_matches(
        sequence,
        recognition_sequence,
        topology=topology,
    )
    # "sites" holds CONFIRMED sites only. Ambiguity-dependent (POSSIBLE)
    # sites are never used as cut positions, but they prevent an enzyme
    # from counting as a unique or absent cutter (see _ambiguity_reason).
    sites = [
        match["position_1based"]
        for match in matches
        if match["status"] == "CONFIRMED"
    ]
    possible_sites = [
        match["position_1based"]
        for match in matches
        if match["status"] == "POSSIBLE"
    ]

    overhang = get_overhang(
        enzyme_name
    )

    return {
        "enzyme": enzyme_name,
        "recognition_sequence": (
            recognition_sequence
        ),
        "site_count": len(sites),
        "sites": sites,
        "possible_site_count": len(possible_sites),
        "possible_sites": possible_sites,
        "site_matches": matches,
        "topology": topology,
        "cut_type": get_cut_type(
            enzyme_name
        ),
        "overhang": overhang,
        "overhang_type": (
            "Blunt"
            if is_blunt_enzyme(
                enzyme_name
            )
            else (
                "Cohesive"
                if overhang
                else "Unknown"
            )
        ),
        "overhang_length": len(
            overhang
        ),
        "overhang_polarity": get_overhang_polarity(
            enzyme_name
        ),
        "is_blunt": is_blunt_enzyme(
            enzyme_name
        ),
    }


def build_restriction_map(
    sequence,
    enzyme_names=None,
    topology="linear",
):
    sequence = normalize_sequence(
        sequence
    )

    if enzyme_names is None:
        enzyme_names = list(
            RESTRICTION_ENZYMES.keys()
        )

    restriction_map = []

    for enzyme_name in enzyme_names:
        if enzyme_name not in RESTRICTION_ENZYMES:
            continue

        analysis = analyze_enzyme(
            sequence,
            enzyme_name,
            topology=topology,
        )

        for position in analysis[
            "sites"
        ]:
            restriction_map.append({
                "enzyme": enzyme_name,
                "position": position,
                "recognition_sequence": (
                    analysis[
                        "recognition_sequence"
                    ]
                ),
                "cut_type": analysis[
                    "cut_type"
                ],
                "overhang": analysis[
                    "overhang"
                ],
            })

    restriction_map.sort(
        key=lambda site: (
            site["position"],
            site["enzyme"],
        )
    )

    return restriction_map


# ============================================================
# VALIDATION
# ============================================================

def _validate_sequence(
    sequence,
    role,
):
    return _cached(
        ("validate", sequence, role),
        lambda: _validate_sequence_uncached(sequence, role),
    )


def _validate_sequence_uncached(
    sequence,
    role,
):
    normalized = normalize_sequence(
        sequence
    )

    invalid_positions = [
        index + 1
        for index, base
        in enumerate(normalized)
        if base not in IUPAC_BASES
    ]

    return {
        "role": role,
        "valid": (
            bool(normalized)
            and not invalid_positions
        ),
        "sequence": normalized,
        "length": len(normalized),
        "invalid_positions": (
            invalid_positions
        ),
        "invalid_symbols": sorted({
            normalized[index - 1]
            for index in invalid_positions
        }),
        "ambiguous_symbols": sorted(
            set(normalized) - set("ACGT")
        ),
    }


def validate_dna_sequence(
    sequence,
    role="sequence",
):
    return _validate_sequence(
        sequence,
        role,
    )
