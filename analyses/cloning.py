# ============================================================
# SYNAPSE - CLONING ENGINE
# ============================================================
#
# Canonical cloning module.
#
# This file is intentionally separate from restriction.py.
# It consumes the restriction-enzyme database but owns:
#   - restriction-site lookup needed for cloning
#   - selected two-end architecture
#   - junction compatibility
#   - directional cloning analysis
#   - reverse-orientation analysis
#   - vector self-ligation assessment
#   - internal/off-target site checks
#   - digest-fragment reasoning
#   - strategy scoring/ranking
#   - automatic cloning planning
#   - terminal reports
#   - pipeline compatibility
#
# No wet-lab success is inferred.
# ============================================================

import copy
from contextlib import contextmanager
from itertools import combinations
from typing import Any, Callable, TypeVar

from data.restriction_enzymes import RESTRICTION_ENZYMES

# Restriction sites are recognized by the restriction engine's matcher, so
# restriction and cloning analyses apply the same CONFIRMED/POSSIBLE rules
# (audit D4). bases_match is re-exported for existing callers.
from analyses.restriction import bases_match, find_recognition_site_matches
from ui.style import major_rule, minor_rule
from ui import style


# ============================================================
# STRATEGY SEARCH CACHE
# ============================================================

# The planner evaluates tens of thousands of strategies over the same two
# molecules, and every evaluation used to re-validate both sequences and
# re-scan them for each enzyme (P2-6). Both steps are pure functions of
# their arguments. Inside _strategy_search_cache() each distinct result is
# computed once and every caller receives its own deep copy, so no result
# shares objects with another. Outside the context nothing is cached.
_STRATEGY_SEARCH_CACHE = None


@contextmanager
def _strategy_search_cache():
    global _STRATEGY_SEARCH_CACHE
    if _STRATEGY_SEARCH_CACHE is not None:   # already inside a search
        yield
        return
    _STRATEGY_SEARCH_CACHE = {}
    try:
        yield
    finally:
        _STRATEGY_SEARCH_CACHE = None


_IMMUTABLE = (str, int, float, bool, type(None))
_T = TypeVar("_T")


def _copy_result(value: Any) -> Any:
    """An independent copy of a result tree (dicts, lists, plain values)."""
    if isinstance(value, _IMMUTABLE):
        return value
    if isinstance(value, dict):
        return {key: _copy_result(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_copy_result(item) for item in value]
    return copy.deepcopy(value)   # anything else: no shortcut


def _cached(key: tuple, compute: Callable[[], _T]) -> _T:
    cache = _STRATEGY_SEARCH_CACHE
    if cache is None:
        return compute()
    if key not in cache:
        cache[key] = compute()
    return _copy_result(cache[key])


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
# REPORTING
# ============================================================

def _yes_no(value):
    return "YES" if value else "NO"


def print_design_reasons(
    analysis,
):
    evaluation = get_design_reasons(
        analysis
    )

    style.section('DESIGN EVALUATION')

    for reason in evaluation[
        "positive"
    ]:
        print(
            f"+ {reason}"
        )

    for warning in evaluation[
        "warnings"
    ]:
        print(
            f"! {warning}"
        )


def ascii_orientation_visualization(
    analysis,
):
    if not analysis:
        return

    vector = analysis[
        "vector"
    ]
    insert = analysis[
        "insert"
    ]
    strategy = analysis[
        "strategy"
    ]

    style.section('CLONING ARCHITECTURE')
    print(
        "VECTOR:"
    )
    print(
        f"{vector['left']['enzyme']}"
        f" ========================== "
        f"{vector['right']['enzyme']}"
    )
    print(
        "                 | INSERT |"
    )
    print(
        "INSERT:"
    )
    print(
        f"{insert['left']['enzyme']}"
        f" >>>>>>>>>>>>>>>>>>>>>>>>>>> "
        f"{insert['right']['enzyme']}"
    )
    print()
    style.field('Orientation', f"{strategy['orientation']}")


def print_cloning_report(
    analysis,
):
    style.header('CLONING ANALYSIS')

    if not analysis:
        print()
        style.note("No analysis available.")
        print(
            major_rule()
        )
        return

    if analysis.get(
        "planner"
    ) is not None:
        planner = analysis[
            "planner"
        ]
        best = planner.get(
            "best"
        )

        if best is not None:
            _print_single_strategy_report(
                best
            )
        else:
            print_planner_report(
                planner
            )

        return

    if "vector" not in analysis:
        print()
        style.field('Status', f"{analysis.get('status', 'UNKNOWN')}")
        print(
            major_rule()
        )
        return

    if analysis.get("mode") == "single-enzyme":
        _print_single_enzyme_report(
            analysis
        )
        return

    _print_single_strategy_report(
        analysis
    )


def _pass_fail(value):
    return "PASS" if value else "FAIL"


def _print_single_enzyme_report(
    analysis,
):
    """Report for analyze_cloning(): one enzyme opens the vector and cuts the insert.

    The two-enzyme report reads left and right ends, which this result does
    not have; it shows the fields analyze_cloning() returns.
    """
    validity = analysis["validity"]
    strategy = analysis["strategy"]
    compatibility = analysis["compatibility"]
    score = analysis["score"]

    for title, site, molecule in (
        ("VECTOR", analysis["vector_site"], analysis["vector"]),
        ("INSERT", analysis["insert_site"], analysis["insert"]),
    ):
        style.section(title)
        style.field('Enzyme', f"{site['enzyme']}")
        style.field('Sites found', f"{molecule['site_count']}")
        style.field('First site position', f"{site['position'] if site['position'] is not None else '-'}")
        style.field('Recognition', f"{site['recognition_sequence']}")
        style.field('Overhang', f"{site['overhang'] or '-'} ({site['overhang_type']})")

    style.section('CLONING JUNCTION')
    style.field('Status', f"{compatibility['status']}")
    style.field('Ends', f"{compatibility['vector_end'] or '-'} -> {compatibility['insert_end'] or '-'}")
    style.field('Mechanism', f"{compatibility['mechanism']}")

    style.section('ORIENTATION')
    style.field('Strategy', f"{strategy['strategy']}")
    style.field('Orientation', f"{strategy['orientation']}")
    style.field('Directional', f"{_yes_no(strategy['directional'])}")
    style.field('Self-ligation risk', f"{'POSSIBLE' if strategy['self_ligation_risk'] else 'REDUCED'}")

    style.section('DESIGN SCORE')
    style.field('Score', f"{score['score']}/{score['maximum']}")
    style.field('Category', f"{score['category']}")
    for reason in score["reasons"]:
        print(f"- {reason}")

    style.section('VALIDATION')
    style.field('Required sites', _pass_fail(validity['all_required_sites']))
    style.field('Vector site unique', _pass_fail(validity['vector_unique']))
    style.field('Insert site unique', _pass_fail(validity['insert_unique']))
    style.field('Junction', _pass_fail(validity['junctions_compatible']))
    style.field('Overall valid', f"{'YES' if validity['valid'] else 'NO'}")

    if analysis.get("warnings"):
        style.section('WARNINGS')
        for warning in analysis["warnings"]:
            print(f"- {warning}")

    style.section('INTERPRETATION')
    if validity["valid"]:
        print("The selected design passes the implemented sequence-level cloning checks.")
    else:
        print("The selected design does not pass all implemented sequence-level cloning checks.")
    print(
        "\nThis is a sequence-level design assessment. "
        "It does not guarantee experimental cleavage, ligation, "
        "enzyme activity, methylation behavior, buffer compatibility, "
        "or wet-lab success."
    )
    print(major_rule())


def _print_single_strategy_report(
    analysis,
):
    vector = analysis[
        "vector"
    ]
    insert = analysis[
        "insert"
    ]
    strategy = analysis[
        "strategy"
    ]
    score = analysis[
        "score"
    ]
    validity = analysis[
        "validity"
    ]

    style.section('VECTOR')
    style.field('Left enzyme', f"{vector['left']['enzyme']}")
    style.field('Left position', f"{vector['left']['position']}")
    style.field('Left recognition', f"{vector['left']['recognition_sequence']}")
    style.field('Left overhang', f"{vector['left']['overhang'] or '-'}")
    style.field('Right enzyme', f"{vector['right']['enzyme']}")
    style.field('Right position', f"{vector['right']['position']}")
    style.field('Right recognition', f"{vector['right']['recognition_sequence']}")
    style.field('Right overhang', f"{vector['right']['overhang'] or '-'}")

    style.section('INSERT')
    style.field('Left enzyme', f"{insert['left']['enzyme']}")
    style.field('Left position', f"{insert['left']['position']}")
    style.field('Left recognition', f"{insert['left']['recognition_sequence']}")
    style.field('Left overhang', f"{insert['left']['overhang'] or '-'}")
    style.field('Right enzyme', f"{insert['right']['enzyme']}")
    style.field('Right position', f"{insert['right']['position']}")
    style.field('Right recognition', f"{insert['right']['recognition_sequence']}")
    style.field('Right overhang', f"{insert['right']['overhang'] or '-'}")

    style.section('CLONING JUNCTIONS')

    for side in (
        "left",
        "right",
    ):
        junction = analysis[
            "junctions"
        ][side]

        print(
            f"{side.upper():<7}"
            f"{junction['vector_enzyme']}"
            f" ({junction['vector_overhang']}) -> "
            f"{junction['insert_enzyme']}"
            f" ({junction['insert_overhang']}) "
            f"{junction['status']}"
        )
        print(
            f"        {junction['mechanism']}"
        )

    style.section('ORIENTATION')
    style.field('Intended orientation', f"{'COMPATIBLE' if validity['junctions_compatible'] else 'NOT VALID'}")
    style.field('Reverse orientation', f"{'COMPATIBLE' if validity['reverse_orientation_compatible'] else 'BLOCKED'}")
    style.field('Directional', f"{_yes_no(strategy['directional'])}")
    style.field('Self-ligation risk', f"{'POSSIBLE' if strategy['self_ligation_risk'] else 'REDUCED'}")

    style.section('DESIGN SCORE')
    style.field('Score', f"{score['score']}/{score['maximum']}")
    style.field('Category', f"{score['category']}")

    for component, value in score[
        "components"
    ].items():
        print(
            f"{component.replace('_', ' ').title():30}"
            f"{value:>3} pts"
        )

    style.section('VALIDATION')
    style.field('Required sites', f"{'PASS' if validity['required_sites'] else 'FAIL'}")
    style.field('Vector architecture', f"{'PASS' if validity['vector_unique'] else 'FAIL'}")
    style.field('Insert architecture', f"{'PASS' if validity['insert_unique'] else 'FAIL'}")
    style.field('Junctions', f"{'PASS' if validity['junctions_compatible'] else 'FAIL'}")
    style.field('No internal sites', f"{'PASS' if validity['no_internal_sites'] else 'FAIL'}")
    style.field('Overall valid', f"{'YES' if validity['valid'] else 'NO'}")

    print_design_reasons(
        analysis
    )

    ascii_orientation_visualization(
        analysis
    )

    style.section('INTERPRETATION')

    if validity[
        "valid"
    ]:
        print(
            "The selected design passes the implemented "
            "sequence-level cloning checks."
        )
    else:
        print(
            "The selected design does not pass all implemented "
            "sequence-level cloning checks."
        )

    print(
        "\nThis is a sequence-level design assessment. "
        "It does not guarantee experimental cleavage, ligation, "
        "enzyme activity, methylation behavior, buffer compatibility, "
        "or wet-lab success."
    )

    print(
        major_rule()
    )


def print_planner_report(
    planner,
):
    style.header('AUTOMATIC CLONING PLANNER')

    if not planner:
        print()
        style.note("No planner result available.")
        print(
            major_rule()
        )
        return

    print()
    style.field('Status', f"{planner.get('status', 'UNKNOWN')}")
    style.field('Strategies evaluated', f"{planner.get('evaluated_count', 0)}")
    style.field('Valid strategies', f"{planner.get('valid_count', 0)}")
    style.field('Rejected strategies', f"{planner.get('rejected_count', 0)}")

    best = planner.get(
        "best"
    )

    if best is None:
        diagnostics = planner.get(
            "diagnostics",
            {},
        )

        style.section('NO STRATEGY AVAILABLE')
        print(
            diagnostics.get(
                "primary_reason",
                "No cloning strategy could be evaluated.",
            )
        )

        recommendations = diagnostics.get(
            "recommendations",
            [],
        )

        if recommendations:
            print(
                "\nRecommendations:"
            )
            for recommendation in recommendations:
                print(
                    f"- {recommendation}"
                )

        for role in (
            "vector",
            "insert",
        ):
            section = diagnostics.get(
                role,
                {},
            )
            rows = section.get(
                "enzyme_diagnostics",
                [],
            )

            if not rows:
                continue

            print(
                f"\n{role.upper()} SITE DIAGNOSTICS"
            )
            print(
                minor_rule()
            )

            for row in rows:
                print(
                    f"{row['enzyme']:<12}"
                    f"{row['site_count']:>3} site(s)  "
                    f"{row['status']:<10}"
                    f"{row['sites']}"
                )

        print(
            major_rule()
        )
        return

    style.section('BEST STRATEGY')

    style.field('Vector', f"{best['vector']['left']['enzyme']} + "
        f"{best['vector']['right']['enzyme']}")
    style.field('Insert', f"{best['insert']['left']['enzyme']} + "
        f"{best['insert']['right']['enzyme']}")
    style.field('Score', f"{best['score']['score']}/"
        f"{best['score']['maximum']}")
    style.field('Category', f"{best['score']['category']}")
    style.field('Orientation', f"{best['strategy']['orientation']}")
    style.field('Directional', f"{_yes_no(best['strategy']['directional'])}")
    style.field('Valid', f"{_yes_no(best['validity']['valid'])}")

    style.section('JUNCTIONS')

    for side in (
        "left",
        "right",
    ):
        junction = best[
            "junctions"
        ][side]

        print(
            f"{side.upper():<7}"
            f"{junction['vector_enzyme']}"
            f" ({junction['vector_overhang']}) -> "
            f"{junction['insert_enzyme']}"
            f" ({junction['insert_overhang']}) "
            f"{junction['status']}"
        )

    style.section('ORIENTATION CHECK')

    style.field('Intended orientation', f"{'COMPATIBLE' if best['validity']['junctions_compatible'] else 'NOT VALID'}")
    style.field('Reverse orientation', f"{'COMPATIBLE' if best['validity']['reverse_orientation_compatible'] else 'BLOCKED'}")

    style.section('VALIDATION')

    validity = best[
        "validity"
    ]

    style.field('Required sites', f"{'PASS' if validity['required_sites'] else 'FAIL'}")
    style.field('Vector architecture', f"{'PASS' if validity['vector_unique'] else 'FAIL'}")
    style.field('Insert architecture', f"{'PASS' if validity['insert_unique'] else 'FAIL'}")
    style.field('Junctions', f"{'PASS' if validity['junctions_compatible'] else 'FAIL'}")
    style.field('No internal sites', f"{'PASS' if validity['no_internal_sites'] else 'FAIL'}")
    style.field('Overall valid', f"{'YES' if validity['valid'] else 'NO'}")

    print_design_reasons(
        best
    )

    style.section('RANKED STRATEGIES')

    for index, candidate in enumerate(
        planner.get(
            "strategies",
            [],
        ),
        start=1,
    ):
        vector = candidate[
            "vector"
        ]
        insert = candidate[
            "insert"
        ]
        strategy = candidate[
            "strategy"
        ]
        score = candidate[
            "score"
        ]

        vector_pair = (
            f"{vector['left']['enzyme']} + "
            f"{vector['right']['enzyme']}"
        )
        insert_pair = (
            f"{insert['left']['enzyme']} + "
            f"{insert['right']['enzyme']}"
        )

        print(
            f"{index:>2}. "
            f"V: {vector_pair:<22} "
            f"I: {insert_pair:<22} "
            f"{score['score']:>3}/100 "
            f"{'DIR' if strategy['directional'] else 'NON-DIR':<7} "
            f"{'VALID' if candidate['validity']['valid'] else 'REJECTED'}"
        )

    print(
        "\nImportant:"
    )
    print(
        "Ranking is based on sequence-level design properties only."
    )

    print(
        major_rule()
    )


# ============================================================
# REPORT ALIASES
# ============================================================

print_advanced_cloning_report = (
    print_cloning_report
)


# ============================================================
# SELF TESTS
# ============================================================

def run_cloning_self_test():
    results = {}

    results[
        "normalize_sequence"
    ] = (
        normalize_sequence(
            " at gc "
        )
        == "ATGC"
    )

    results[
        "reverse_complement"
    ] = (
        reverse_complement(
            "ATGC"
        )
        == "GCAT"
    )

    results[
        "iupac_matching"
    ] = (
        bases_match(
            "A",
            "R",
        )
        and bases_match(
            "G",
            "R",
        )
        and not bases_match(
            "C",
            "R",
        )
    )

    results[
        "enzyme_database"
    ] = all(
        name in RESTRICTION_ENZYMES
        for name in list(
            RESTRICTION_ENZYMES
        )[:1]
    )

    try:
        invalid = automatic_cloning_planner(
            "",
            "ATGC",
        )
        results[
            "invalid_vector"
        ] = (
            invalid[
                "status"
            ]
            == "INVALID_VECTOR"
            and invalid[
                "best"
            ] is None
        )
    except Exception:
        results[
            "invalid_vector"
        ] = False

    try:
        invalid = automatic_cloning_planner(
            "ATGC",
            "",
        )
        results[
            "invalid_insert"
        ] = (
            invalid[
                "status"
            ]
            == "INVALID_INSERT"
            and invalid[
                "best"
            ] is None
        )
    except Exception:
        results[
            "invalid_insert"
        ] = False

    try:
        empty = automatic_cloning_planner(
            "ATGC",
            "ATGC",
            enzyme_names=[],
        )
        results[
            "empty_enzyme_set"
        ] = (
            empty[
                "best"
            ] is None
            and empty[
                "evaluated_count"
            ] == 0
        )
    except Exception:
        results[
            "empty_enzyme_set"
        ] = False

    try:
        one_site = _end_definition(
            "ATGCGAATTCATGC",
            "EcoRI",
            "EcoRI",
        )
        results[
            "same_enzyme_one_site_rejected"
        ] = (
            one_site[
                "valid"
            ] is False
        )
    except Exception:
        results[
            "same_enzyme_one_site_rejected"
        ] = False

    try:
        two_sites = _end_definition(
            "ATGCGAATTCATGCGAATTCATGC",
            "EcoRI",
            "EcoRI",
        )
        results[
            "same_enzyme_two_sites_accepted"
        ] = (
            two_sites[
                "valid"
            ] is True
            and two_sites[
                "left_site"
            ]["position"]
            != two_sites[
                "right_site"
            ]["position"]
        )
    except Exception:
        results[
            "same_enzyme_two_sites_accepted"
        ] = False

    try:
        blunt_a = {
            "enzyme": "A",
            "position": 1,
            "overhang": "",
            "is_blunt": True,
        }
        blunt_b = {
            "enzyme": "B",
            "position": 1,
            "overhang": "",
            "is_blunt": True,
        }
        sticky = {
            "enzyme": "C",
            "position": 1,
            "overhang": "AATT",
            "is_blunt": False,
        }

        results[
            "blunt_compatibility"
        ] = (
            ends_compatible(
                blunt_a,
                blunt_b,
            )
            is True
        )

        results[
            "blunt_sticky_incompatibility"
        ] = (
            ends_compatible(
                blunt_a,
                sticky,
            )
            is False
        )
    except Exception:
        results[
            "blunt_compatibility"
        ] = False
        results[
            "blunt_sticky_incompatibility"
        ] = False

    required_api = [
        "analyze_enzyme",
        "assess_enzyme_suitability",
        "ends_compatible",
        "analyze_end_compatibility",
        "analyze_junction",
        "analyze_directional_strategy",
        "score_cloning_strategy",
        "build_restriction_map",
        "analyze_enzyme_pair",
        "find_best_cloning_strategies",
        "automatic_cloning_planner",
        "find_compatible_enzymes",
        "compare_cloning_strategies",
        "get_design_reasons",
        "print_cloning_report",
        "print_planner_report",
        "analyze_cloning",
        "run_cloning_analysis",
        "get_cloning_summary",
        "get_best_cloning_strategy",
        "get_cloning_junctions",
    ]

    results[
        "public_api"
    ] = all(
        callable(
            globals().get(
                name
            )
        )
        for name in required_api
    )

    results[
        "all_passed"
    ] = all(
        results.values()
    )

    return results


# ============================================================
# MODULE SELF TEST
# ============================================================

if __name__ == "__main__":
    tests = run_cloning_self_test()

    style.header('Cloning engine self-test')

    for name, passed in tests.items():
        print(
            f"{name:<48}"
            f"{'PASS' if passed else 'FAIL'}"
        )

    print(
        minor_rule()
    )

    print(
        "RESULT: "
        + (
            "ALL CLONING ENGINE TESTS PASSED"
            if tests[
                "all_passed"
            ]
            else
            "ONE OR MORE CLONING ENGINE TESTS FAILED"
        )
    )

    print(
        major_rule()
    )

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




def _cloning_ui_select_insert(
    sequences,
    vector_index=None,
):
    if not sequences:
        print()
        style.note("No additional sequences are available.")
        return None, None, None

    records = []

    if isinstance(sequences, dict):
        records = list(sequences.items())
    else:
        for index, record in enumerate(sequences):
            if isinstance(record, dict):
                records.append((index, record))

    available = []

    style.header("SELECT INSERT")

    for key, record in records:
        sequence = record.get("sequence", "")
        header = record.get("header", str(key))

        if not sequence:
            continue

        if vector_index is not None and key == vector_index:
            continue

        available.append((key, header, sequence))

    if not available:
        style.note("No suitable insert sequence is available.")
        return None, None, None

    style.menu([("Insert", [(index, header) for index, (_, header, _) in
                            enumerate(available, start=1)])])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None, None, None

        for display_index, (key, header, sequence) in enumerate(
            available,
            start=1,
        ):
            if choice == str(display_index):
                return key, header, sequence

        style.invalid_selection()


def _cloning_ui_record_topology(
    sequences,
    key,
):
    """Topology of the record chosen as insert ("linear" if not recorded)."""
    try:
        record = sequences[key]
    except (KeyError, IndexError, TypeError):
        return "linear"

    if not isinstance(record, dict):
        return "linear"

    return record.get("topology", "linear")


def _cloning_ui_select_enzyme(
    prompt,
    allow_back=True,
):
    enzyme_names = list(RESTRICTION_ENZYMES.keys())

    style.header("SELECT RESTRICTION ENZYME")

    style.menu([("Enzymes", list(enumerate(enzyme_names, start=1)))],
               zero=style.BACK_LABEL if allow_back else None)

    while True:
        choice = style.ask_choice(prompt)

        if allow_back and choice == "0":
            return None

        try:
            index = int(choice) - 1
        except ValueError:
            style.invalid_selection()
            continue

        if 0 <= index < len(enzyme_names):
            return enzyme_names[index]

        style.invalid_selection()


def _cloning_ui_select_two_distinct_enzymes():
    first = _cloning_ui_select_enzyme(
        "Left enzyme"
    )

    if first is None:
        return None, None

    while True:
        second = _cloning_ui_select_enzyme(
            "Right enzyme"
        )

        if second is None:
            return None, None

        if second == first:
            print()
            style.error("The two enzymes must be different.")
            continue

        return first, second


def _cloning_ui_store_result(
    pipeline,
    result,
):
    if pipeline is None:
        return

    try:
        pipeline.context.store_cloning_analysis(
            result
        )
    except AttributeError:
        pass


def _cloning_ui_print_sequence_context(
    vector_header,
    vector_sequence,
    insert_header,
    insert_sequence,
    vector_topology="linear",
    insert_topology="linear",
):
    style.field('Vector', f"{vector_header}")
    style.field('Vector length', f"{len(normalize_sequence(vector_sequence))} bp ({vector_topology})")
    style.field('Insert', f"{insert_header}")
    style.field('Insert length', f"{len(normalize_sequence(insert_sequence))} bp ({insert_topology})")


def _cloning_ui_automatic_planner(
    vector_sequence,
    insert_sequence,
    pipeline=None,
    vector_topology="linear",
    insert_topology="linear",
):
    style.header(
        "AUTOMATIC CLONING PLANNER"
    )

    print(
        "\nSynapse evaluates the available enzyme architectures "
        "and ranks the resulting sequence-level strategies."
    )

    style.menu([
        (None, [('1', 'Find best strategy'), ('2', 'Show top strategies')]),
    ])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice in {"1", "2"}:
            break

        style.invalid_selection()

    if choice == "1":
        top_n = 1
    else:
        while True:
            raw = style.ask('Number of strategies to display [default 5]')

            if not raw:
                top_n = 5
                break

            try:
                top_n = int(raw)
                if top_n < 1:
                    raise ValueError
                break
            except ValueError:
                style.error("Please enter a positive integer.")

    result = automatic_cloning_planner(
        vector_sequence,
        insert_sequence,
        top_n=top_n,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )

    _cloning_ui_store_result(
        pipeline,
        result,
    )

    print_planner_report(result)
    style.pause()
    return result


def _cloning_ui_single_enzyme(
    vector_sequence,
    insert_sequence,
    pipeline=None,
    vector_topology="linear",
    insert_topology="linear",
):
    style.header(
        "SINGLE-ENZYME ANALYSIS"
    )

    print(
        "\nThe same restriction enzyme is evaluated "
        "for vector and insert."
    )

    style.menu([
        (None, [('1', 'Select enzyme')]),
    ])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            break

        style.invalid_selection()

    enzyme = _cloning_ui_select_enzyme(
        "Selection"
    )

    if enzyme is None:
        return None

    result = analyze_cloning(
        vector_sequence,
        insert_sequence,
        enzyme,
        enzyme,
        vector_topology=vector_topology,
        insert_topology=insert_topology,
    )

    _cloning_ui_store_result(
        pipeline,
        result,
    )

    print_cloning_report(result)
    style.pause()
    return result


def _cloning_ui_select_two_enzyme_architecture(
    molecule_name,
):
    style.header(
        f"SELECT {molecule_name.upper()} ENZYMES"
    )

    print(
        f"\nSelect the two restriction enzymes for the "
        f"{molecule_name.lower()} ends."
    )
    print(
        "The two enzymes must be different."
    )
    print()

    left_enzyme = _cloning_ui_select_enzyme(
        f"{molecule_name} left enzyme"
    )

    if left_enzyme is None:
        return None, None

    right_enzyme = _cloning_ui_select_enzyme(
        f"{molecule_name} right enzyme"
    )

    if right_enzyme is None:
        return None, None

    if right_enzyme == left_enzyme:
        print()
        style.error("The two enzymes must be different.")
        style.pause()
        return _cloning_ui_select_two_enzyme_architecture(
            molecule_name
        )

    return left_enzyme, right_enzyme


def _cloning_ui_two_enzyme_strategy(
    vector_sequence,
    insert_sequence,
    pipeline=None,
    vector_topology="linear",
    insert_topology="linear",
):
    style.header(
        "TWO-ENZYME STRATEGY ANALYSIS"
    )

    print(
        "\nVector and insert are configured independently."
    )
    print(
        "This allows Synapse to evaluate four enzyme positions:"
    )
    print(
        "  Vector left / Vector right / Insert left / Insert right"
    )

    style.menu([
        (None, [('1', 'Configure enzyme architecture')]),
    ])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            break

        style.invalid_selection()

    vector_left, vector_right = (
        _cloning_ui_select_two_enzyme_architecture(
            "Vector"
        )
    )

    if vector_left is None:
        return None

    insert_left, insert_right = (
        _cloning_ui_select_two_enzyme_architecture(
            "Insert"
        )
    )

    if insert_left is None:
        return None

    style.header(
        "TWO-ENZYME STRATEGY"
    )

    print()
    style.field('Vector', f"{vector_left} + {vector_right}")
    style.field('Insert', f"{insert_left} + {insert_right}")

    style.menu([
        (None, [('1', 'Run analysis')]),
    ])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            break

        style.invalid_selection()

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

    _cloning_ui_store_result(
        pipeline,
        result,
    )

    print_cloning_report(result)
    style.pause()
    return result


def run_cloning_analysis_ui(
    vector_sequence,
    insert_sequence=None,
    vector_header="Vector",
    insert_header="Insert",
    pipeline=None,
    sequences=None,
    vector_index=None,
    vector_topology="linear",
    insert_topology=None,
):
    """
    Interactive CLI for the canonical cloning engine.

    UI structure:
        CLONING ANALYSIS
            [1] Automatic Cloning Planner
            [2] Single-Enzyme Analysis
            [3] Two-Enzyme Strategy Analysis
            [0] Back

    Sequence-specific input and enzyme selection are requested only
    after the corresponding analysis is selected.

    The cloning engine itself is not modified. An insert chosen from
    ``sequences`` keeps the topology of its record; otherwise
    insert_topology (default linear) applies.
    """
    vector_validation = validate_dna_sequence(
        vector_sequence,
        "vector",
    )

    if not vector_validation["valid"]:
        style.header(
            "CLONING ANALYSIS"
        )
        print()
        style.error("Invalid vector sequence.")

        if vector_validation["invalid_positions"]:
            style.field('Invalid positions', f"{vector_validation['invalid_positions']}")

        print(minor_rule())
        return None

    if insert_sequence is None:
        insert_key, insert_header, insert_sequence = (
            _cloning_ui_select_insert(
                sequences,
                vector_index=vector_index,
            )
        )

        if insert_sequence is None:
            return None

        if insert_topology is None:
            insert_topology = _cloning_ui_record_topology(
                sequences,
                insert_key,
            )

    vector_topology = _normalize_topology(
        vector_topology
    )
    insert_topology = _normalize_topology(
        insert_topology or "linear"
    )

    insert_validation = validate_dna_sequence(
        insert_sequence,
        "insert",
    )

    if not insert_validation["valid"]:
        style.header(
            "CLONING ANALYSIS"
        )
        print()
        style.error("Invalid insert sequence.")

        if insert_validation["invalid_positions"]:
            style.field('Invalid positions', f"{insert_validation['invalid_positions']}")

        print(minor_rule())
        return None

    while True:
        style.header(
            "CLONING ANALYSIS"
        )

        _cloning_ui_print_sequence_context(
            vector_header,
            vector_sequence,
            insert_header,
            insert_sequence,
            vector_topology=vector_topology,
            insert_topology=insert_topology,
        )

        style.menu([
            ("Planning", [('1', 'Automatic Cloning Planner')]),
            ("Analysis", [('2', 'Single-Enzyme Analysis'), ('3', 'Two-Enzyme Strategy Analysis')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            result = _cloning_ui_automatic_planner(
                vector_sequence,
                insert_sequence,
                pipeline=pipeline,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )
            if result is not None:
                continue
            return None

        if choice == "2":
            result = _cloning_ui_single_enzyme(
                vector_sequence,
                insert_sequence,
                pipeline=pipeline,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )
            if result is not None:
                continue
            return None

        if choice == "3":
            result = _cloning_ui_two_enzyme_strategy(
                vector_sequence,
                insert_sequence,
                pipeline=pipeline,
                vector_topology=vector_topology,
                insert_topology=insert_topology,
            )
            if result is not None:
                continue
            return None

        style.invalid_selection()


# ============================================================
# SELF-TEST ALIAS
# ============================================================

def run_cloning_ui_self_test():
    return callable(
        globals().get(
            "run_cloning_analysis_ui"
        )
    )
