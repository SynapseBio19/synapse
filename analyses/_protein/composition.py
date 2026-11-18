"""Composition, molecular weight, extinction coefficient and composition indices."""

from __future__ import annotations

from analyses._protein.constants import (
    AA_ORDER,
    AMBIGUOUS_RESIDUE_OPTIONS,
    AMINO_ACID_RESIDUE_MASSES,
    AROMATIC_RESIDUES,
    INSTABILITY_WEIGHTS,
    NONSTANDARD_RESIDUES,
    STANDARD_AMINO_ACIDS,
    UNRESOLVED_RESIDUES,
    WATER_MASS,
)
from analyses._protein.sequence import (
    get_positional_protein,
    get_resolved_protein,
    _strip_terminal_stop,
)


def analyze_protein(protein):

    if not isinstance(protein, str):
        protein = ""

    sequence = "".join(
        protein.upper().split()
    )

    counts = {
        amino_acid: 0
        for amino_acid in AA_ORDER
    }

    ambiguous = 0
    nonstandard_counts = {}

    for residue in sequence:

        if residue in STANDARD_AMINO_ACIDS:
            counts[residue] += 1

        elif residue in NONSTANDARD_RESIDUES:
            nonstandard_counts[residue] = nonstandard_counts.get(residue, 0) + 1

        elif residue in UNRESOLVED_RESIDUES:
            ambiguous += 1

    counts = {
        amino_acid: count
        for amino_acid, count
        in counts.items()
        if count > 0
    }

    return {
        "counts":
            counts,

        "nonstandard_counts":
            nonstandard_counts,

        "ambiguous_count":
            ambiguous,

        "resolved_length":
            sum(counts.values()),

        "total_length":
            len(sequence),
    }


def calculate_amino_acid_percentages(
    protein
):

    analysis = analyze_protein(
        protein
    )

    resolved_length = (
        analysis["resolved_length"]
    )

    if resolved_length == 0:
        return {}

    return {
        amino_acid:
            count
            / resolved_length
            * 100
        for amino_acid, count
        in analysis["counts"].items()
    }


# ============================================================
# MOLECULAR WEIGHT
# ============================================================


def calculate_protein_molecular_weight(
    protein,
    translation_analysis=None
):

    if not isinstance(protein, str):
        return {
            "min": None,
            "max": None,
            "ambiguous": False,
            "undetermined": True,
            "known_mass": 0.0,
        }

    sequence = protein.upper()

    if not sequence:
        return {
            "min": 0.0,
            "max": 0.0,
            "ambiguous": False,
            "undetermined": False,
            "known_mass": 0.0,
        }

    if translation_analysis is None:

        # B, Z and J give a mass range over their residues; X and "?" leave
        # the mass undetermined; a terminal "*" is a stop, not a residue.
        sequence = _strip_terminal_stop("".join(sequence.split()))
        total = WATER_MASS
        total_min = WATER_MASS
        total_max = WATER_MASS
        unknown = False
        ranged = False

        for residue in sequence:

            if residue in AMINO_ACID_RESIDUE_MASSES:
                mass = AMINO_ACID_RESIDUE_MASSES[residue]
                total += mass
                total_min += mass
                total_max += mass
            elif residue in AMBIGUOUS_RESIDUE_OPTIONS:
                masses = [AMINO_ACID_RESIDUE_MASSES[option]
                          for option in AMBIGUOUS_RESIDUE_OPTIONS[residue]]
                total_min += min(masses)
                total_max += max(masses)
                ranged = True
            else:
                unknown = True

        if unknown:
            return {
                "min": None,
                "max": None,
                "ambiguous": True,
                "undetermined": True,
                "known_mass": total,
            }

        return {
            "min": total_min,
            "max": total_max,
            "ambiguous": ranged,
            "undetermined": False,
            "known_mass": total,
        }

    # Translation-aware calculation.
    #
    # The protein represented by an ambiguous codon can have
    # multiple possible masses. An ambiguous-stop codon can also
    # terminate translation, so a simple per-residue min/max is
    # insufficient. We therefore keep separate mass ranges for
    # translation paths that are still active and paths that have
    # already terminated.
    active_min = 0.0
    active_max = 0.0

    completed_min = None
    completed_max = None

    ambiguous = False
    known_mass = WATER_MASS
    protein_index = 0

    for record in translation_analysis.get(
        "codon_records",
        []
    ):

        status = record["status"]

        if status == "stop":

            if active_min is not None:
                completed_min = (
                    active_min
                    if completed_min is None
                    else min(
                        completed_min,
                        active_min
                    )
                )

                completed_max = (
                    active_max
                    if completed_max is None
                    else max(
                        completed_max,
                        active_max
                    )
                )

            active_min = None
            active_max = None
            break

        residue = (
            sequence[protein_index]
            if protein_index < len(sequence)
            else "?"
        )

        if residue in AMINO_ACID_RESIDUE_MASSES:

            mass = AMINO_ACID_RESIDUE_MASSES[
                residue
            ]

            active_min += mass
            active_max += mass
            known_mass += mass

        elif status in {
            "ambiguous_amino_acid",
            "ambiguous_stop",
        }:

            ambiguous = True

            possible_amino_acids = [
                amino_acid
                for amino_acid in record.get(
                    "possibilities",
                    []
                )
                if amino_acid
                in AMINO_ACID_RESIDUE_MASSES
            ]

            masses = [
                AMINO_ACID_RESIDUE_MASSES[
                    amino_acid
                ]
                for amino_acid
                in possible_amino_acids
            ]

            if not masses:

                return {
                    "min": None,
                    "max": None,
                    "ambiguous": True,
                    "undetermined": True,
                    "known_mass": known_mass,
                }

            if status == "ambiguous_stop":

                # STOP is one possible outcome. It terminates
                # the current active path without adding a residue.
                if "STOP" in record.get(
                    "possibilities",
                    []
                ):

                    stop_min = active_min
                    stop_max = active_max

                    completed_min = (
                        stop_min
                        if completed_min is None
                        else min(
                            completed_min,
                            stop_min
                        )
                    )

                    completed_max = (
                        stop_max
                        if completed_max is None
                        else max(
                            completed_max,
                            stop_max
                        )
                    )

                # Non-stop amino-acid outcomes continue.
                active_min += min(masses)
                active_max += max(masses)

            else:

                active_min += min(masses)
                active_max += max(masses)

        else:

            return {
                "min": None,
                "max": None,
                "ambiguous": True,
                "undetermined": True,
                "known_mass": known_mass,
            }

        protein_index += 1

    # If translation never terminated, the active paths are valid
    # final products as well.
    final_min = completed_min
    final_max = completed_max

    # The two bounds become None together (undetermined residue).
    if active_min is not None and active_max is not None:

        final_min = (
            active_min
            if final_min is None
            else min(
                final_min,
                active_min
            )
        )

        final_max = (
            active_max
            if final_max is None
            else max(
                final_max,
                active_max
            )
        )

    if final_min is None or final_max is None:

        return {
            "min": None,
            "max": None,
            "ambiguous": ambiguous,
            "undetermined": True,
            "known_mass": known_mass,
        }

    final_min += WATER_MASS
    final_max += WATER_MASS

    return {
        "min": final_min,
        "max": final_max,
        "ambiguous": ambiguous,
        "undetermined": False,
        "known_mass": known_mass,
    }


def calculate_aliphatic_index(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return None

    length = len(protein)

    a = protein.count("A") * 100 / length
    v = protein.count("V") * 100 / length
    i = protein.count("I") * 100 / length
    leu = protein.count("L") * 100 / length

    return (
        a
        + 2.9 * v
        + 3.9 * (i + leu)
    )


# ============================================================
# AROMATICITY
# ============================================================

def calculate_aromaticity(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    if not protein:
        return None

    return sum(
        protein.count(residue)
        for residue in AROMATIC_RESIDUES
    ) / len(protein)


# ============================================================
# EXTINCTION COEFFICIENT
# ============================================================


def calculate_extinction_coefficient(
    protein
):

    protein = get_resolved_protein(
        protein
    )

    tryptophan = protein.count("W")
    tyrosine = protein.count("Y")
    cysteine = protein.count("C")

    base = (
        tryptophan * 5500
        + tyrosine * 1490
    )

    # Pace et al. values at 280 nm in water.
    # All cysteines reduced: no disulfide contribution.
    reduced = base

    # Maximum possible number of disulfide bonds from the
    # available cysteines. If the Cys count is odd, one Cys
    # necessarily remains unpaired.
    disulfide_pairs = cysteine // 2
    unpaired_cysteine = cysteine % 2

    oxidized = (
        base
        + disulfide_pairs * 125
    )

    return {
        "reduced_cysteine": reduced,
        "oxidized_cysteine": oxidized,
        "tryptophan": tryptophan,
        "tyrosine": tyrosine,
        "cysteine": cysteine,
        "disulfide_pairs": disulfide_pairs,
        "unpaired_cysteine": unpaired_cysteine,
        "wavelength_nm": 280,
        "units": "M^-1 cm^-1",
    }


def calculate_instability_index(
    protein
):

    protein = get_positional_protein(
        protein
    )

    if len(protein) < 2:
        return None

    score = 0.0
    scored = 0

    for i in range(
        len(protein) - 1
    ):

        first = protein[i]
        second = protein[i + 1]

        if first == "?" or second == "?":
            # Unresolved residue: this dipeptide has no weight. Skipping it
            # keeps every other dipeptide at its true position.
            continue

        scored += 1

        score += (
            INSTABILITY_WEIGHTS[
                first
            ][
                second
            ]
        )

    if scored == 0:
        return None

    # 10/L * sum over dipeptides (Guruprasad et al. 1990); L is the number of
    # scored residues, which equals len(protein) for a fully resolved sequence.
    return (
        10.0
        / (scored + 1)
        * score
    )
