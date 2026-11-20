"""Composition, content profiles, quality, molecular weight, complexity and homopolymers."""

from __future__ import annotations

import math
from functools import lru_cache

from analyses._rna.bases import (
    count_rna_bases,
    RNA_AMBIGUOUS_CODES,
    RNA_IUPAC_BASES,
    RNA_NUCLEOTIDE_RESIDUE_MASSES,
    RNA_STANDARD_BASES,
)
from core.windows import base_prefix_counts, nonstandard_prefix_count


def calculate_rna_base_percentages(
    base_counts
):

    total = sum(
        base_counts.values()
    )

    if total == 0:

        return {
            "A": 0.0,
            "U": 0.0,
            "C": 0.0,
            "G": 0.0
        }

    return {
        base:
            (
                count / total
            ) * 100
        for base, count
        in base_counts.items()
    }


# ============================================================
# GC / AU CONTENT
# ============================================================

def calculate_rna_gc_content(rna):

    counts = count_rna_bases(rna)

    total = sum(
        counts.values()
    )

    if total == 0:
        return 0.0

    return (
        (
            counts["G"]
            + counts["C"]
        )
        / total
    ) * 100


def calculate_rna_au_content(rna):

    counts = count_rna_bases(rna)

    total = sum(
        counts.values()
    )

    if total == 0:
        return 0.0

    return (
        (
            counts["A"]
            + counts["U"]
        )
        / total
    ) * 100


def calculate_rna_content_ranges(rna):

    valid_positions = 0

    gc_min = 0
    gc_max = 0

    au_min = 0
    au_max = 0

    ambiguous = False

    for symbol in rna:

        possible = RNA_IUPAC_BASES.get(
            symbol
        )

        if possible is None:
            continue

        valid_positions += 1

        if len(possible) > 1:
            ambiguous = True

        gc_values = [
            base in {"G", "C"}
            for base in possible
        ]

        au_values = [
            base in {"A", "U"}
            for base in possible
        ]

        gc_min += int(
            min(gc_values)
        )

        gc_max += int(
            max(gc_values)
        )

        au_min += int(
            min(au_values)
        )

        au_max += int(
            max(au_values)
        )

    if valid_positions == 0:

        return {
            "gc_min": 0.0,
            "gc_max": 0.0,
            "au_min": 0.0,
            "au_max": 0.0,
            "ambiguous": False,
            "evaluated_length": 0
        }

    return {
        "gc_min":
            gc_min
            / valid_positions
            * 100,

        "gc_max":
            gc_max
            / valid_positions
            * 100,

        "au_min":
            au_min
            / valid_positions
            * 100,

        "au_max":
            au_max
            / valid_positions
            * 100,

        "ambiguous":
            ambiguous,

        "evaluated_length":
            valid_positions
    }


# ============================================================
# SKEW
# ============================================================

def calculate_rna_gc_skew(rna):

    counts = count_rna_bases(rna)

    denominator = (
        counts["G"]
        + counts["C"]
    )

    if denominator == 0:
        return 0.0

    return (
        counts["G"]
        - counts["C"]
    ) / denominator


def calculate_rna_au_skew(rna):

    counts = count_rna_bases(rna)

    denominator = (
        counts["A"]
        + counts["U"]
    )

    if denominator == 0:
        return 0.0

    return (
        counts["A"]
        - counts["U"]
    ) / denominator


# ============================================================
# SLIDING WINDOW GC PROFILE
# ============================================================

def calculate_rna_gc_profile(
    rna,
    window_size=100,
    step=1
):

    profile = []

    if (
        window_size <= 0
        or step <= 0
        or len(rna) < window_size
    ):
        return profile

    # Prefix counts instead of slicing every window: the same integer
    # counts and formula, in O(length) (core/windows.py).
    counts = base_prefix_counts(rna, ("G", "C"))
    g_prefix, c_prefix = counts["G"], counts["C"]
    invalid_prefix = nonstandard_prefix_count(rna, RNA_STANDARD_BASES)

    for start in range(
        0,
        len(rna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        if invalid_prefix[end] == invalid_prefix[start]:

            gc = (
                (
                    (g_prefix[end] - g_prefix[start])
                    + (c_prefix[end] - c_prefix[start])
                )
                / window_size
            ) * 100

            profile.append({
                "start":
                    start,

                "end":
                    end - 1,

                "start_1_based":
                    start + 1,

                "end_1_based":
                    end,

                "gc_content":
                    gc,

                "resolved":
                    True
            })

        else:

            profile.append({
                "start":
                    start,

                "end":
                    end - 1,

                "start_1_based":
                    start + 1,

                "end_1_based":
                    end,

                "gc_content":
                    None,

                "resolved":
                    False
            })

    return profile


# ============================================================
# AMBIGUITY / QUALITY
# ============================================================

def analyze_rna_ambiguity(rna):

    ambiguous_counts = {}
    ambiguous_positions = []

    for position, base in enumerate(rna):

        if base in RNA_AMBIGUOUS_CODES:

            ambiguous_positions.append({
                "position":
                    position,

                "position_1_based":
                    position + 1,

                "base":
                    base,

                "possibilities":
                    sorted(
                        RNA_IUPAC_BASES[base]
                    )
            })

            ambiguous_counts[base] = (
                ambiguous_counts.get(
                    base,
                    0
                ) + 1
            )

    count = len(
        ambiguous_positions
    )

    percentage = (
        count / len(rna) * 100
        if rna
        else 0.0
    )

    return {
        "ambiguous_count":
            count,

        "ambiguous_percentage":
            percentage,

        "ambiguous_counts":
            ambiguous_counts,

        "ambiguous_positions":
            ambiguous_positions,

        "unambiguous":
            count == 0
    }


def analyze_rna_quality(rna):

    standard = 0
    ambiguous = 0
    invalid = 0

    invalid_positions = []

    for position, symbol in enumerate(
        rna
    ):

        if symbol in RNA_STANDARD_BASES:

            standard += 1

        elif symbol in RNA_AMBIGUOUS_CODES:

            ambiguous += 1

        else:

            invalid += 1

            invalid_positions.append({
                "position":
                    position,

                "position_1_based":
                    position + 1,

                "symbol":
                    symbol
            })

    length = len(rna)

    return {
        "standard_count":
            standard,

        "ambiguous_count":
            ambiguous,

        "invalid_count":
            invalid,

        "standard_percentage":
            standard / length * 100
            if length else 0.0,

        "ambiguous_percentage":
            ambiguous / length * 100
            if length else 0.0,

        "invalid_percentage":
            invalid / length * 100
            if length else 0.0,

        "invalid_positions":
            invalid_positions,

        "valid":
            invalid == 0,

        "unambiguous":
            (
                ambiguous == 0
                and invalid == 0
            )
    }


# ============================================================
# MOLECULAR WEIGHT
# ============================================================

def calculate_rna_molecular_weight(rna):

    minimum_weight = 0.0
    maximum_weight = 0.0

    ambiguous = False
    evaluated_positions = 0

    for symbol in rna:

        possible = RNA_IUPAC_BASES.get(
            symbol
        )

        if possible is None:
            continue

        evaluated_positions += 1

        if len(possible) > 1:
            ambiguous = True

        weights = [
            RNA_NUCLEOTIDE_RESIDUE_MASSES[
                base
            ]
            for base in possible
        ]

        minimum_weight += min(
            weights
        )

        maximum_weight += max(
            weights
        )

    return {
        "min":
            minimum_weight,

        "max":
            maximum_weight,

        "ambiguous":
            ambiguous,

        "evaluated_positions":
            evaluated_positions
    }


# ============================================================
# COMPLEXITY
# ============================================================

def calculate_rna_complexity(rna):

    counts = count_rna_bases(rna)

    # Counts in the order A, U, C, G (count_rna_bases), as before.
    return _rna_complexity_from_counts(tuple(counts.values()))


@lru_cache(maxsize=65536)
def _rna_complexity_from_counts(counts):
    """Normalized Shannon entropy from the A, U, C, G counts (cached; same arithmetic)."""

    total = sum(
        counts
    )

    if total == 0:
        return 0.0

    entropy = 0.0

    for count in counts:

        if count == 0:
            continue

        probability = (
            count / total
        )

        entropy -= (
            probability
            * math.log2(probability)
        )

    return entropy / 2.0


# ============================================================
# COMPLEXITY PROFILE
# ============================================================

def calculate_rna_complexity_profile(
    rna,
    window_size=100,
    step=1
):

    profile = []

    if (
        window_size <= 0
        or step <= 0
        or len(rna) < window_size
    ):
        return profile

    # Per-window counts of A, U, C, G by prefix subtraction: the counts
    # count_rna_bases takes from the slice, in the same order.
    counts = base_prefix_counts(rna, ("A", "U", "C", "G"))
    prefixes = [counts[base] for base in ("A", "U", "C", "G")]

    for start in range(
        0,
        len(rna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        profile.append({
            "start":
                start,

            "end":
                end - 1,

            "start_1_based":
                start + 1,

            "end_1_based":
                end,

            "complexity":
                _rna_complexity_from_counts(tuple(
                    prefix[end] - prefix[start]
                    for prefix in prefixes
                ))
        })

    return profile


# ============================================================
# HOMOPOLYMERS
# ============================================================

def find_rna_homopolymers(
    rna,
    minimum_length=4
):

    regions = []

    if not rna:
        return regions

    start = 0

    for index in range(
        1,
        len(rna) + 1
    ):

        if (
            index < len(rna)
            and rna[index] == rna[start]
            and rna[start]
            in RNA_STANDARD_BASES
        ):
            continue

        length = index - start

        if (
            length >= minimum_length
            and rna[start]
            in RNA_STANDARD_BASES
        ):

            regions.append({
                "base":
                    rna[start],

                "length":
                    length,

                "start":
                    start + 1,

                "end":
                    index
            })

        start = index

    return regions


# ============================================================
# POLY-A / POLY-U
# ============================================================

def find_rna_homopolymer_regions(
    rna,
    base,
    minimum_length=4
):

    regions = []

    start = None

    for index, symbol in enumerate(
        rna
    ):

        if symbol == base:

            if start is None:
                start = index

        else:

            if start is not None:

                length = index - start

                if length >= minimum_length:

                    regions.append({
                        "base":
                            base,

                        "start":
                            start + 1,

                        "end":
                            index,

                        "length":
                            length
                    })

                start = None

    if start is not None:

        length = len(rna) - start

        if length >= minimum_length:

            regions.append({
                "base":
                    base,

                "start":
                    start + 1,

                "end":
                    len(rna),

                "length":
                    length
            })

    return regions
