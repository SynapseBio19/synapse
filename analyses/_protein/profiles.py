"""Sequence profiles and region predictions (heuristics, see METHODS)."""

from __future__ import annotations

from analyses._protein.constants import (
    FLEXIBILITY_SCALE,
    KYTE_DOOLITTLE,
    NEGATIVE_RESIDUES,
    POSITIVE_RESIDUES,
)
from analyses._protein.sequence import get_positional_protein, get_resolved_protein


def calculate_gravy(protein):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return None

    return sum(
        KYTE_DOOLITTLE[
            residue
        ]
        for residue in protein
    ) / len(protein)


def calculate_hydropathy_profile(
    protein,
    window=9
):

    protein = get_positional_protein(
        protein
    )

    try:
        window = int(window)
    except (TypeError, ValueError):
        raise ValueError("window must be an integer.")

    if window < 1:
        raise ValueError("window must be at least 1.")

    if len(protein) < window:
        return []

    profile = []

    for start in range(
        len(protein) - window + 1
    ):

        segment = protein[
            start:
            start + window
        ]

        if "?" in segment:

            # Unresolved position: no score instead of a shifted window.

            continue

        score = (
            sum(
                KYTE_DOOLITTLE[residue]
                for residue in segment
            )
            / window
        )

        profile.append({
            "start":
                start + 1,

            "end":
                start + window,

            "score":
                score,

            "hydropathy":
                score,

            "sequence":
                segment
        })

    return profile


def find_hydrophobic_regions(
    protein,
    window=19,
    threshold=1.6
):

    # ExPASy ProtScale notes that windows of 19 or 21 residues
    # make membrane-spanning hydrophobic segments stand out;
    # 1.6 is a commonly used Kyte-Doolittle threshold.
    profile = calculate_hydropathy_profile(
        protein,
        window
    )

    return [
        region
        for region in profile
        if region["score"] >= threshold
    ]


def find_transmembrane_domains(
    protein,
    window=19,
    hydropathy_threshold=1.6,
    max_charged=2
):

    protein = get_positional_protein(
        protein
    )

    try:
        window = int(window)
    except (TypeError, ValueError):
        raise ValueError("window must be an integer.")

    if window < 3:
        raise ValueError("window must be at least 3.")

    if window % 2 == 0:
        raise ValueError(
            "Transmembrane window must be odd so that it has a central residue."
        )

    if len(protein) < window:
        return []

    candidates = []

    for start in range(
        len(protein) - window + 1
    ):

        segment = protein[
            start:
            start + window
        ]

        if "?" in segment:

            # Unresolved position: no score instead of a shifted window.

            continue

        hydropathy = (
            sum(
                KYTE_DOOLITTLE[residue]
                for residue in segment
            )
            / window
        )

        charged = sum(
            residue in POSITIVE_RESIDUES
            or residue in NEGATIVE_RESIDUES
            for residue in segment
        )

        if (
            hydropathy >= hydropathy_threshold
            and charged <= max_charged
        ):

            candidates.append({
                "start": start + 1,
                "end": start + window,
            })

    if not candidates:
        return []

    # Consecutive sliding windows overlap heavily. Merge only
    # overlapping candidate windows, then recompute all reported
    # properties on the merged sequence instead of carrying the
    # first window's score into the merged region.
    merged_intervals = []

    for candidate in candidates:

        if not merged_intervals:
            merged_intervals.append(
                candidate.copy()
            )
            continue

        previous = merged_intervals[-1]

        if candidate["start"] <= previous["end"]:

            previous["end"] = max(
                previous["end"],
                candidate["end"]
            )

        else:

            merged_intervals.append(
                candidate.copy()
            )

    domains = []

    for interval in merged_intervals:

        start = interval["start"]
        end = interval["end"]

        sequence = protein[
            start - 1:end
        ]

        hydropathy = (
            sum(
                KYTE_DOOLITTLE[residue]
                for residue in sequence
            )
            / len(sequence)
        )

        charged = sum(
            residue in POSITIVE_RESIDUES
            or residue in NEGATIVE_RESIDUES
            for residue in sequence
        )

        domains.append({
            "start": start,
            "end": end,
            "length": len(sequence),
            "hydropathy": hydropathy,
            "charged_residues": charged,
            "sequence": sequence,
            "method": (
                "Kyte-Doolittle sliding-window "
                "sequence heuristic"
            ),
        })

    return domains


# ============================================================
# SECONDARY STRUCTURE TENDENCIES
# ============================================================


def calculate_secondary_structure_propensity(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:

        return {
            "helix": 0.0,
            "turn": 0.0,
            "sheet": 0.0,
            "method": "Chou-Fasman",
            "units": "propensity",
        }

    # Classic Chou-Fasman conformational parameters:
    # P(alpha), P(beta), P(turn).
    # Values above 1 indicate a relative preference for that
    # conformational class; an average propensity is NOT a
    # percentage of residues adopting that structure.
    chou_fasman = {
        "A": (1.42, 0.83, 0.66),
        "R": (0.98, 0.93, 0.95),
        "N": (0.67, 0.89, 1.56),
        "D": (1.01, 0.54, 1.46),
        "C": (0.70, 1.19, 1.19),
        "Q": (1.11, 1.10, 0.98),
        "E": (1.51, 0.37, 0.74),
        "G": (0.57, 0.75, 1.56),
        "H": (1.00, 0.87, 0.95),
        "I": (1.08, 1.60, 0.47),
        "L": (1.21, 1.30, 0.59),
        "K": (1.16, 0.74, 1.01),
        "M": (1.45, 1.05, 0.60),
        "F": (1.13, 1.38, 0.60),
        "P": (0.57, 0.55, 1.52),
        "S": (0.77, 0.75, 1.43),
        "T": (0.83, 1.19, 0.96),
        "W": (1.08, 1.37, 0.96),
        "Y": (0.69, 1.47, 1.14),
        "V": (1.06, 1.70, 0.50),
    }

    length = len(protein)

    return {
        "helix": (
            sum(
                chou_fasman[residue][0]
                for residue in protein
            )
            / length
        ),
        "turn": (
            sum(
                chou_fasman[residue][2]
                for residue in protein
            )
            / length
        ),
        "sheet": (
            sum(
                chou_fasman[residue][1]
                for residue in protein
            )
            / length
        ),
        "method": "Chou-Fasman",
        "units": "propensity",
    }


def calculate_secondary_structure_fraction(
    protein
):
    """
    Backward-compatible wrapper.

    The historical function name is retained so existing callers
    do not break, but the returned values are propensities, not
    structural fractions.
    """
    return calculate_secondary_structure_propensity(
        protein
    )


def calculate_flexibility_profile(
    protein,
    window=9
):

    protein = get_positional_protein(
        protein
    )

    try:
        window = int(window)
    except (TypeError, ValueError):
        raise ValueError("window must be an integer.")

    if window < 1:
        raise ValueError("window must be at least 1.")

    if len(protein) < window:
        return []

    profile = []

    for start in range(
        len(protein) - window + 1
    ):

        segment = protein[
            start:
            start + window
        ]

        if "?" in segment:

            # Unresolved position: no score instead of a shifted window.

            continue

        score = (
            sum(
                FLEXIBILITY_SCALE[residue]
                for residue in segment
            )
            / window
        )

        profile.append({
            "start": start + 1,
            "end": start + window,
            "score": score,
            "sequence": segment,
            "scale": "Vihinen normalized flexibility",
        })

    return profile


def predict_signal_peptide(
    protein
):

    protein = get_positional_protein(
        protein
    )

    if len(protein) < 15:

        return {
            "predicted": False,
            "evidence": "low",
            "reason": "Sequence too short for a meaningful signal-peptide screen.",
            "hydrophobic_core": None,
            "cleavage_site": None,
            "method": "N-terminal sequence heuristic",
        }

    # This is deliberately a screen, not a SignalP-equivalent
    # predictor. A classical secretory signal peptide usually has
    # an N-region, a hydrophobic h-region and a C-region containing
    # a plausible cleavage site.
    n_terminal = protein[
        :min(35, len(protein))
    ]

    hydrophobic_core = None

    for window in range(8, 16):

        if len(n_terminal) < window:
            continue

        for start in range(
            0,
            min(
                len(n_terminal) - window + 1,
                13
            )
        ):

            segment = n_terminal[
                start:
                start + window
            ]

            if "?" in segment:
                # Unresolved position: no score instead of a shifted window.
                continue

            score = (
                sum(
                    KYTE_DOOLITTLE[residue]
                    for residue in segment
                )
                / window
            )

            charged = sum(
                residue in POSITIVE_RESIDUES
                or residue in NEGATIVE_RESIDUES
                for residue in segment
            )

            if (
                score >= 1.4
                and charged <= 1
            ):

                candidate = {
                    "start": start + 1,
                    "end": start + window,
                    "score": score,
                    "sequence": segment,
                    "charged_residues": charged,
                }

                if (
                    hydrophobic_core is None
                    or candidate["score"]
                    > hydrophobic_core["score"]
                ):
                    hydrophobic_core = candidate

    if hydrophobic_core is None:

        return {
            "predicted": False,
            "evidence": "low",
            "reason": (
                "No convincing N-terminal hydrophobic core "
                "was detected."
            ),
            "hydrophobic_core": None,
            "cleavage_site": None,
            "method": "N-terminal sequence heuristic",
        }

    # Simplified -3/-1 cleavage-site screen. Small residues are
    # favored at these positions in classical signal peptidase
    # cleavage motifs. This is only a plausibility filter.
    cleavage_residues = set(
        "AVGSTC"
    )

    cleavage_candidates = []

    minimum_cleavage = max(
        hydrophobic_core["end"] + 2,
        15
    )

    maximum_cleavage = min(
        30,
        len(n_terminal) - 1
    )

    for cleavage in range(
        minimum_cleavage,
        maximum_cleavage + 1
    ):

        minus_1 = n_terminal[
            cleavage - 1
        ]

        minus_3 = n_terminal[
            cleavage - 3
        ]

        if (
            minus_1 in cleavage_residues
            and minus_3 in cleavage_residues
            and minus_1 != "P"
            and minus_3 != "P"
        ):

            cleavage_candidates.append(
                cleavage
            )

    if not cleavage_candidates:

        return {
            "predicted": False,
            "evidence": "low",
            "reason": (
                "N-terminal hydrophobic core detected, "
                "but no plausible cleavage-site pattern "
                "was found."
            ),
            "hydrophobic_core": hydrophobic_core,
            "cleavage_site": None,
            "method": "N-terminal sequence heuristic",
        }

    cleavage_site = min(
        cleavage_candidates,
        key=lambda position:
            abs(
                position
                - (
                    hydrophobic_core["end"]
                    + 4
                )
            )
    )

    evidence = (
        "moderate"
        if hydrophobic_core["score"] >= 1.8
        else "low"
    )

    return {
        "predicted": True,
        "evidence": evidence,
        "reason": (
            "N-terminal hydrophobic core and a plausible "
            "-3/-1 cleavage-site pattern were detected."
        ),
        "hydrophobic_core": hydrophobic_core,
        "cleavage_site": cleavage_site,
        "method": "N-terminal sequence heuristic",
    }
