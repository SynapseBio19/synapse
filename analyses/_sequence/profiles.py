"""Sliding-window GC and complexity profiles, their statistics, homopolymers."""

from __future__ import annotations

import math
from array import array
from functools import lru_cache

from analyses._sequence.bases import DNA_STANDARD_BASES, _get_standard_bases
from core.windows import (
    base_prefix_counts,
    DEFAULT_MAX_PROFILE_POINTS,  # noqa: F401  (re-exported to analyses.sequence)
    nonstandard_prefix_count,
    window_count,
)


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


def gc_window_values(
    dna,
    window_size=100,
    step=1
):
    """Yield (0-based start, GC % or None) for every window of the GC profile.

    None marks a window with a non-standard base. Prefix counts instead of
    slicing every window: the same integer counts and the same formula as
    _calculate_gc_for_window, in O(length) (core/windows.py).
    """

    if window_size <= 0 or step <= 0 or len(dna) < window_size:
        return

    counts = base_prefix_counts(dna, ("G", "C"))
    g_prefix, c_prefix = counts["G"], counts["C"]
    invalid_prefix = nonstandard_prefix_count(dna)

    for start in range(
        0,
        len(dna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        if invalid_prefix[end] == invalid_prefix[start]:

            gc_count = (
                (g_prefix[end] - g_prefix[start])
                + (c_prefix[end] - c_prefix[start])
            )

            yield start, (
                gc_count
                / window_size
            ) * 100

        else:

            yield start, None


def calculate_gc_content_profile(
    dna,
    window_size=100,
    step=1
):

    return [
        {

            "start":
                start,

            "end":
                start + window_size - 1,

            "start_1_based":
                start + 1,

            "end_1_based":
                start + window_size,

            "gc_content":
                gc_content,

            "resolved":
                gc_content is not None
        }
        for start, gc_content in gc_window_values(dna, window_size, step)
    ]



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

    # Sorted: a set iterates in a different order in every process, which
    # changed the last digit of the entropy sum and, through ties, the
    # reported position of the complexity maximum (fixed 2026-11-19).
    return _complexity_from_counts(tuple(
        standard_bases.count(base)
        for base in _SORTED_STANDARD_BASES
    ))


_SORTED_STANDARD_BASES = tuple(sorted(DNA_STANDARD_BASES))


@lru_cache(maxsize=65536)
def _complexity_from_counts(counts):
    """Normalized Shannon entropy (0-1) from the A, C, G, T counts, in that order.

    Many windows share a composition, so the result is cached; the
    arithmetic and its order are unchanged, so cached and computed values
    are the same float.
    """

    total = sum(counts)

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
            * math.log2(
                probability
            )
        )

    return entropy / 2.0


# ============================================================
# SEQUENCE COMPLEXITY PROFILE
# ============================================================


def complexity_window_values(
    dna,
    window_size=100,
    step=1
):
    """Yield (0-based start, normalized entropy) for every profile window.

    Per-window base counts by prefix subtraction: the same counts that
    calculate_sequence_complexity takes from the slice.
    """

    if window_size <= 0 or step <= 0 or len(dna) < window_size:
        return

    counts = base_prefix_counts(dna, _SORTED_STANDARD_BASES)
    prefixes = [counts[base] for base in _SORTED_STANDARD_BASES]

    for start in range(
        0,
        len(dna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        yield start, _complexity_from_counts(tuple(
            prefix[end] - prefix[start]
            for prefix in prefixes
        ))


def calculate_complexity_profile(
    dna,
    window_size=100,
    step=1
):

    return [
        {

            "start":
                start,

            "end":
                start + window_size - 1,

            "start_1_based":
                start + 1,

            "end_1_based":
                start + window_size,

            "complexity":
                complexity
        }
        for start, complexity in complexity_window_values(dna, window_size, step)
    ]


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
# STREAMING STATISTICS AND LONG SEQUENCES
# ============================================================

def profile_statistics_from_values(window_values, window_size):
    """calculate_profile_statistics without building the profile.

    The same result, value for value: the minimum and maximum keep the
    first window that reaches them (as min() and max() do), and the
    average is sum() over the values in window order, so the float is
    identical (sum() of floats is compensated on Python >= 3.12 and plain
    on 3.10; either way it is the same function on the same sequence).
    """

    values = array("d")
    minimum = maximum = None
    minimum_start = maximum_start = 0

    for start, value in window_values:

        if value is None:
            continue

        values.append(value)

        if minimum is None or value < minimum:
            minimum, minimum_start = value, start

        if maximum is None or value > maximum:
            maximum, maximum_start = value, start

    if not values:

        return {

            "count": 0,

            "minimum": None,

            "maximum": None,

            "average": None,

            "minimum_position": None,

            "maximum_position": None
        }

    return {

        "count":
            len(values),

        "minimum":
            minimum,

        "maximum":
            maximum,

        "average":
            sum(values) / len(values),

        "minimum_position": {
            "start": minimum_start + 1,
            "end": minimum_start + window_size
        },

        "maximum_position": {
            "start": maximum_start + 1,
            "end": maximum_start + window_size
        }
    }


def _profile_stride(length, window_size, step, max_points):
    total = window_count(length, window_size, step)
    if max_points is None or max_points <= 0 or total <= max_points:
        return 1, total
    return -(-total // max_points), total


def _analyze_profile(dna, window_size, step, max_points, window_values, build, value_key):
    stride, total = _profile_stride(len(dna), window_size, step, max_points)

    if stride == 1:
        profile = build(dna, window_size, step)
        statistics = calculate_profile_statistics(profile, value_key)
    else:
        statistics = profile_statistics_from_values(
            window_values(dna, window_size, step),
            window_size
        )
        # Every stride-th window of the full profile is the profile at
        # step x stride: the same windows, starting at 0.
        profile = build(dna, window_size, step * stride)

    result = {

        "window_size":
            window_size,

        "step":
            step,

        "profile":
            profile,

        "statistics":
            statistics
    }

    if stride > 1:
        result["profile_stride"] = stride
        result["profile_windows_total"] = total

    return result


# ============================================================
# GC PROFILE ANALYSIS
# ============================================================


def analyze_gc_profile(
    dna,
    window_size=100,
    step=1,
    max_points=None
):
    """GC profile and its statistics.

    With ``max_points`` and more windows than that, the statistics are
    still computed over every window, but the stored profile keeps every
    k-th window only: k is "profile_stride" and "profile_windows_total"
    the full window count. Without ``max_points`` (the default) the full
    profile is stored, as before.
    """

    return _analyze_profile(dna, window_size, step, max_points,
                            gc_window_values, calculate_gc_content_profile, "gc_content")


# ============================================================
# COMPLEXITY PROFILE ANALYSIS
# ============================================================


def analyze_complexity_profile(
    dna,
    window_size=100,
    step=1,
    max_points=None
):
    """Complexity profile and its statistics; ``max_points`` as in analyze_gc_profile."""

    return _analyze_profile(dna, window_size, step, max_points,
                            complexity_window_values, calculate_complexity_profile, "complexity")


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
