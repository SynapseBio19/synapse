"""Composition, GC/AT content, skews, sequence quality and molecular weight."""

from __future__ import annotations

from analyses._sequence.bases import (
    AVERAGE_DSDNA_BP_MASS,
    count_nucleotides,
    DNA_AMBIGUOUS_CODES,
    DNA_COMPLEMENT,
    DNA_IUPAC_BASES,
    DNA_NUCLEOTIDE_RESIDUE_MASSES,
    DNA_STANDARD_BASES,
    _get_valid_iupac_bases,
)


# ============================================================
# SEQUENCE STATISTICS
# ============================================================

def calculate_sequence_statistics(
    dna,
    nucleotide_counts=None,
    nucleotide_percentages=None
):
    """
    Compact top-level statistics for the Sequence Overview.

    All counts are based on standard A/T/C/G bases.
    Ambiguous and invalid symbols are reported separately.
    """

    if nucleotide_counts is None:
        nucleotide_counts = count_nucleotides(dna)

    if nucleotide_percentages is None:
        nucleotide_percentages = (
            calculate_nucleotide_percentages(
                nucleotide_counts
            )
        )

    standard_length = sum(
        nucleotide_counts.values()
    )

    gc_count = (
        nucleotide_counts["G"]
        + nucleotide_counts["C"]
    )

    at_count = (
        nucleotide_counts["A"]
        + nucleotide_counts["T"]
    )

    return {
        "length": len(dna),
        "standard_length": standard_length,
        "gc_count": gc_count,
        "at_count": at_count,
        "gc_content": (
            gc_count / standard_length * 100
            if standard_length
            else 0.0
        ),
        "at_content": (
            at_count / standard_length * 100
            if standard_length
            else 0.0
        ),
        "nucleotide_counts": dict(
            nucleotide_counts
        ),
        "nucleotide_percentages": dict(
            nucleotide_percentages
        )
    }


def _calculate_gc_for_window(window):

    if not window:
        return 0.0

    gc_count = (
        window.count("G")
        + window.count("C")
    )

    return (
        gc_count
        / len(window)
    ) * 100


# ============================================================
# SEQUENCE QUALITY
# ============================================================


def analyze_sequence_quality(dna):

    standard_count = 0
    ambiguous_count = 0
    invalid_count = 0

    ambiguous_positions = []
    invalid_positions = []

    for position, nucleotide in enumerate(dna):

        if nucleotide in DNA_STANDARD_BASES:

            standard_count += 1

        elif nucleotide in DNA_AMBIGUOUS_CODES:

            ambiguous_count += 1

            ambiguous_positions.append({

                "position":
                    position,

                "position_1_based":
                    position + 1,

                "base":
                    nucleotide,

                "possibilities":
                    sorted(
                        DNA_IUPAC_BASES[
                            nucleotide
                        ]
                    )
            })

        else:

            invalid_count += 1

            invalid_positions.append({

                "position":
                    position,

                "position_1_based":
                    position + 1,

                "symbol":
                    nucleotide
            })

    sequence_length = len(dna)

    if sequence_length > 0:

        standard_percentage = (
            standard_count
            / sequence_length
        ) * 100

        ambiguous_percentage = (
            ambiguous_count
            / sequence_length
        ) * 100

        invalid_percentage = (
            invalid_count
            / sequence_length
        ) * 100

    else:

        standard_percentage = 0.0
        ambiguous_percentage = 0.0
        invalid_percentage = 0.0

    return {

        "standard_bases":
            standard_count,

        "ambiguous_bases":
            ambiguous_count,

        "invalid_symbols":
            invalid_count,

        "standard_percentage":
            standard_percentage,

        "ambiguous_percentage":
            ambiguous_percentage,

        "invalid_percentage":
            invalid_percentage,

        "ambiguous_positions":
            ambiguous_positions,

        "invalid_positions":
            invalid_positions,

        "valid":
            invalid_count == 0,

        "unambiguous":
            (
                ambiguous_count == 0
                and invalid_count == 0
            )
    }


# ============================================================
# NUCLEOTIDE PERCENTAGES
# ============================================================


def calculate_nucleotide_percentages(
    nucleotide_counts
):

    total = sum(
        nucleotide_counts.values()
    )

    if total == 0:

        return {
            "A": 0.0,
            "T": 0.0,
            "C": 0.0,
            "G": 0.0
        }

    percentages = {}

    for nucleotide, count in (
        nucleotide_counts.items()
    ):

        percentages[nucleotide] = (
            count / total
        ) * 100

    return percentages


# ============================================================
# GC / AT CONTENT
# ============================================================


def calculate_gc_content(dna):

    counts = count_nucleotides(dna)

    known_bases = sum(
        counts.values()
    )

    if known_bases == 0:

        return 0.0

    return (

        (
            counts["G"]
            + counts["C"]
        )

        / known_bases

    ) * 100


def calculate_at_content(dna):

    counts = count_nucleotides(dna)

    known_bases = sum(
        counts.values()
    )

    if known_bases == 0:

        return 0.0

    return (

        (
            counts["A"]
            + counts["T"]
        )

        / known_bases

    ) * 100


# ============================================================
# IUPAC-AWARE GC / AT RANGE
# ============================================================


def calculate_dna_content_ranges(dna):

    valid_symbols = _get_valid_iupac_bases(
        dna
    )

    evaluated_length = len(
        valid_symbols
    )

    if evaluated_length == 0:

        return {

            "gc_min": 0.0,
            "gc_max": 0.0,

            "at_min": 0.0,
            "at_max": 0.0,

            "ambiguous": False,

            "evaluated_length": 0
        }

    gc_min_count = 0
    gc_max_count = 0

    at_min_count = 0
    at_max_count = 0

    ambiguous = False

    for symbol in valid_symbols:

        possible_bases = (
            DNA_IUPAC_BASES[symbol]
        )

        if len(possible_bases) > 1:

            ambiguous = True

        possible_gc = [

            base in {"G", "C"}

            for base in possible_bases
        ]

        possible_at = [

            base in {"A", "T"}

            for base in possible_bases
        ]

        gc_min_count += int(
            min(possible_gc)
        )

        gc_max_count += int(
            max(possible_gc)
        )

        at_min_count += int(
            min(possible_at)
        )

        at_max_count += int(
            max(possible_at)
        )

    return {

        "gc_min":
            gc_min_count
            / evaluated_length
            * 100,

        "gc_max":
            gc_max_count
            / evaluated_length
            * 100,

        "at_min":
            at_min_count
            / evaluated_length
            * 100,

        "at_max":
            at_max_count
            / evaluated_length
            * 100,

        "ambiguous":
            ambiguous,

        "evaluated_length":
            evaluated_length
    }


# ============================================================
# NUCLEOTIDE SKEW
# ============================================================


def calculate_gc_skew(dna):

    g_count = dna.count("G")
    c_count = dna.count("C")

    total_gc = (
        g_count + c_count
    )

    if total_gc == 0:

        return 0.0

    return (
        g_count - c_count
    ) / total_gc


def calculate_at_skew(dna):

    a_count = dna.count("A")
    t_count = dna.count("T")

    total_at = (
        a_count + t_count
    )

    if total_at == 0:

        return 0.0

    return (
        a_count - t_count
    ) / total_at



# ============================================================
# DNA MOLECULAR WEIGHT
# ============================================================


def calculate_dna_molecular_weight(dna):

    """
    Composition-based residue estimate.

    The exact dsDNA estimate is obtained by pairing each possible
    base with its Watson-Crick complement. For ambiguous IUPAC
    symbols, a minimum and maximum are reported.

    This is deliberately labeled an estimate: nucleotide residue
    mass conventions and terminal-group assumptions differ between
    molecular-weight calculators.
    """

    ss_minimum_weight = 0.0
    ss_maximum_weight = 0.0

    ds_minimum_weight = 0.0
    ds_maximum_weight = 0.0

    ambiguous = False
    evaluated_positions = 0

    for symbol in dna:

        possible_bases = DNA_IUPAC_BASES.get(
            symbol
        )

        if possible_bases is None:
            continue

        evaluated_positions += 1

        if len(possible_bases) > 1:
            ambiguous = True

        possible_ss_weights = []
        possible_ds_weights = []

        for base in possible_bases:

            complement = DNA_COMPLEMENT[
                base
            ]

            ss_weight = (
                DNA_NUCLEOTIDE_RESIDUE_MASSES[
                    base
                ]
            )

            ds_weight = (
                DNA_NUCLEOTIDE_RESIDUE_MASSES[
                    base
                ]
                +
                DNA_NUCLEOTIDE_RESIDUE_MASSES[
                    complement
                ]
            )

            possible_ss_weights.append(
                ss_weight
            )

            possible_ds_weights.append(
                ds_weight
            )

        ss_minimum_weight += min(
            possible_ss_weights
        )

        ss_maximum_weight += max(
            possible_ss_weights
        )

        ds_minimum_weight += min(
            possible_ds_weights
        )

        ds_maximum_weight += max(
            possible_ds_weights
        )

    return {

        "ss_min":
            ss_minimum_weight,

        "ss_max":
            ss_maximum_weight,

        "ds_min":
            ds_minimum_weight,

        "ds_max":
            ds_maximum_weight,

        # Backward-compatible keys.
        "min":
            ds_minimum_weight,

        "max":
            ds_maximum_weight,

        "ambiguous":
            ambiguous,

        "evaluated_positions":
            evaluated_positions,

        "method":
            "Composition-based nucleotide residue estimate"
    }



# ============================================================
# AVERAGE DNA MOLECULAR WEIGHT
# ============================================================


def calculate_average_dna_molecular_weight(
    dna
):

    """
    Conventional approximate dsDNA mass using 660 Da per bp.

    This is intentionally separate from the composition-based
    residue estimate above.
    """

    valid_length = len(
        _get_valid_iupac_bases(dna)
    )

    return (
        valid_length
        * AVERAGE_DSDNA_BP_MASS
    )
