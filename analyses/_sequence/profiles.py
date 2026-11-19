"""Sliding-window GC and complexity profiles, their statistics, homopolymers."""

from __future__ import annotations

import math

from analyses._sequence.bases import DNA_STANDARD_BASES, _get_standard_bases
from analyses._sequence.composition import _calculate_gc_for_window


# ============================================================
# PROFILE AVAILABILITY
# ============================================================

def describe_profile_availability(
    dna,
    window_size,
    profile_name
):
    """
    Explain why a sliding-window profile is unavailable for a
    sequence shorter than the configured window.
    """

    if len(dna) >= window_size:
        return {
            "available": True,
            "reason": None,
            "message": None
        }

    return {
        "available": False,
        "reason": "SEQUENCE_SHORTER_THAN_WINDOW",
        "message": (
            f"{profile_name} is not available for this sequence. "
            f"Sequence length: {len(dna)} bp. "
            f"Minimum configured window: {window_size} bp."
        )
    }


# ============================================================
# SLIDING-WINDOW GC CONTENT
# ============================================================


def calculate_gc_content_profile(
    dna,
    window_size=100,
    step=1
):

    profile = []

    if window_size <= 0:
        return profile

    if step <= 0:
        return profile

    if len(dna) < window_size:
        return profile

    for start in range(
        0,
        len(dna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        window = dna[
            start:end
        ]

        if all(
            base in DNA_STANDARD_BASES
            for base in window
        ):

            gc_content = (
                _calculate_gc_for_window(
                    window
                )
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

                "gc_content":
                    gc_content,

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
# SEQUENCE COMPLEXITY
# ============================================================


def calculate_sequence_complexity(dna):

    standard_bases = _get_standard_bases(
        dna
    )

    total = len(
        standard_bases
    )

    if total == 0:

        return 0.0

    counts = {

        base:
            standard_bases.count(base)

        # Sorted: a set iterates in a different order in every process, which
        # changed the last digit of the entropy sum and, through ties, the
        # reported position of the complexity maximum (fixed 2026-11-19).
        for base in sorted(DNA_STANDARD_BASES)
    }

    entropy = 0.0

    for count in counts.values():

        if count == 0:

            continue

        probability = (
            count / total
        )

        entropy -= (

            probability
            * math.log2(
                probability
            )
        )

    return entropy / 2.0


# ============================================================
# SEQUENCE COMPLEXITY PROFILE
# ============================================================


def calculate_complexity_profile(
    dna,
    window_size=100,
    step=1
):

    profile = []

    if window_size <= 0:

        return profile

    if step <= 0:

        return profile

    if len(dna) < window_size:

        return profile

    for start in range(
        0,
        len(dna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        window = dna[
            start:end
        ]

        complexity = (
            calculate_sequence_complexity(
                window
            )
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
                complexity
        })

    return profile


# ============================================================
# PROFILE STATISTICS
# ============================================================


def calculate_profile_statistics(
    profile,
    value_key
):

    values = []

    for item in profile:

        value = item.get(
            value_key
        )

        if value is not None:

            values.append(value)

    if not values:

        return {

            "count": 0,

            "minimum": None,

            "maximum": None,

            "average": None,

            "minimum_position": None,

            "maximum_position": None
        }

    minimum = min(values)
    maximum = max(values)

    average = (
        sum(values)
        / len(values)
    )

    minimum_item = min(
        (
            item
            for item in profile
            if item.get(value_key) is not None
        ),
        key=lambda item:
            item[value_key]
    )

    maximum_item = max(
        (
            item
            for item in profile
            if item.get(value_key) is not None
        ),
        key=lambda item:
            item[value_key]
    )

    return {

        "count":
            len(values),

        "minimum":
            minimum,

        "maximum":
            maximum,

        "average":
            average,

        "minimum_position": {

            "start":
                minimum_item[
                    "start_1_based"
                ],

            "end":
                minimum_item[
                    "end_1_based"
                ]
        },

        "maximum_position": {

            "start":
                maximum_item[
                    "start_1_based"
                ],

            "end":
                maximum_item[
                    "end_1_based"
                ]
        }
    }


# ============================================================
# GC PROFILE ANALYSIS
# ============================================================


def analyze_gc_profile(
    dna,
    window_size=100,
    step=1
):

    profile = calculate_gc_content_profile(
        dna,
        window_size,
        step
    )

    statistics = calculate_profile_statistics(
        profile,
        "gc_content"
    )

    return {

        "window_size":
            window_size,

        "step":
            step,

        "profile":
            profile,

        "statistics":
            statistics
    }


# ============================================================
# COMPLEXITY PROFILE ANALYSIS
# ============================================================


def analyze_complexity_profile(
    dna,
    window_size=100,
    step=1
):

    profile = calculate_complexity_profile(
        dna,
        window_size,
        step
    )

    statistics = calculate_profile_statistics(
        profile,
        "complexity"
    )

    return {

        "window_size":
            window_size,

        "step":
            step,

        "profile":
            profile,

        "statistics":
            statistics
    }


# ============================================================
# HOMOPOLYMERS
# ============================================================


def find_homopolymers(
    dna,
    min_length=4
):

    homopolymers = []

    if not dna:

        return homopolymers

    start = 0

    for i in range(
        1,
        len(dna) + 1
    ):

        if (

            i == len(dna)

            or dna[i] != dna[start]
        ):

            length = (
                i - start
            )

            base = dna[start]

            if (

                base in DNA_STANDARD_BASES

                and length >= min_length
            ):

                homopolymers.append({

                    "base":
                        base,

                    "start":
                        start,

                    "end":
                        i - 1,

                    "start_1_based":
                        start + 1,

                    "end_1_based":
                        i,

                    "length":
                        length
                })

            start = i

    return homopolymers
