"""CpG observed/expected statistics and CpG island search."""

from __future__ import annotations

from analyses._sequence.bases import DNA_IUPAC_BASES, DNA_STANDARD_BASES, _get_valid_iupac_bases


# ============================================================
# CpG CRITERIA
# ============================================================

DEFAULT_CPG_ISLAND_CRITERIA = {
    "window_size": 200,
    "min_gc_content": 50.0,
    "min_oe_ratio": 0.60,
    "definition": (
        "Classical Gardiner-Garden/Frommer-style screening "
        "criteria"
    )
}


# ============================================================
# CpG ANALYSIS
# ============================================================


def analyze_cpg(dna):

    evaluated_length = len(
        _get_valid_iupac_bases(dna)
    )

    if evaluated_length == 0:

        return {

            "count": 0,

            "density": 0.0,

            "oe_ratio": 0.0,

            "ambiguous_contexts": 0
        }

    definite_cpg_count = 0
    ambiguous_contexts = 0

    for i in range(
        len(dna) - 1
    ):

        first = dna[i]
        second = dna[i + 1]

        if (

            first not in DNA_IUPAC_BASES

            or second not in DNA_IUPAC_BASES
        ):

            continue

        if (

            first == "C"

            and second == "G"
        ):

            definite_cpg_count += 1

            continue

        first_possibilities = (
            DNA_IUPAC_BASES[first]
        )

        second_possibilities = (
            DNA_IUPAC_BASES[second]
        )

        if (

            "C" in first_possibilities

            and

            "G" in second_possibilities
        ):

            ambiguous_contexts += 1

    c_count = dna.count("C")
    g_count = dna.count("G")

    density = (

        definite_cpg_count
        / evaluated_length

    ) * 100

    if (

        c_count == 0

        or g_count == 0
    ):

        oe_ratio = 0.0

    else:

        expected_cpg = (

            c_count
            * g_count

        ) / evaluated_length

        if expected_cpg == 0:

            oe_ratio = 0.0

        else:

            oe_ratio = (

                definite_cpg_count
                / expected_cpg
            )

    return {

        "count":
            definite_cpg_count,

        "density":
            density,

        "oe_ratio":
            oe_ratio,

        "ambiguous_contexts":
            ambiguous_contexts
    }


# ============================================================
# CpG ISLANDS
# ============================================================


def find_cpg_islands(
    dna,
    window_size=200,
    min_gc=50.0,
    min_oe=0.60
):

    qualifying_windows = []

    if window_size <= 0:
        return []

    if len(dna) < window_size:
        return []

    for start in range(
        len(dna) - window_size + 1
    ):

        end = (
            start + window_size
        )

        window = dna[
            start:end
        ]

        if any(
            base not in DNA_STANDARD_BASES
            for base in window
        ):
            continue

        c_count = window.count("C")
        g_count = window.count("G")
        cpg_count = window.count("CG")

        gc_content = (
            (c_count + g_count)
            / window_size
        ) * 100

        if c_count == 0 or g_count == 0:
            oe_ratio = 0.0
        else:
            expected_cpg = (
                c_count * g_count
            ) / window_size

            oe_ratio = (
                cpg_count
                / expected_cpg
            )

        if (
            gc_content >= min_gc
            and oe_ratio >= min_oe
        ):
            qualifying_windows.append(
                (start, end - 1)
            )

    if not qualifying_windows:
        return []

    # Merge overlapping qualifying windows.
    merged_regions = []

    current_start, current_end = (
        qualifying_windows[0]
    )

    for start, end in qualifying_windows[1:]:

        if start <= current_end + 1:
            current_end = max(
                current_end,
                end
            )
        else:
            merged_regions.append(
                (current_start, current_end)
            )

            current_start = start
            current_end = end

    merged_regions.append(
        (current_start, current_end)
    )

    islands = []

    for start, end in merged_regions:

        region = dna[
            start:end + 1
        ]

        length = len(region)

        if length < window_size:
            continue

        c_count = region.count("C")
        g_count = region.count("G")
        cpg_count = region.count("CG")

        gc_content = (
            (c_count + g_count)
            / length
        ) * 100

        if c_count == 0 or g_count == 0:
            oe_ratio = 0.0
        else:
            expected_cpg = (
                c_count * g_count
            ) / length

            oe_ratio = (
                cpg_count
                / expected_cpg
            )

        # Do not report a merged region unless the FINAL
        # reported region itself still satisfies the criteria.
        if (
            gc_content < min_gc
            or oe_ratio < min_oe
        ):
            continue

        islands.append({

            "start":
                start,

            "end":
                end,

            "start_1_based":
                start + 1,

            "end_1_based":
                end + 1,

            "length":
                length,

            "gc_content":
                gc_content,

            "cpg_count":
                cpg_count,

            "oe_ratio":
                oe_ratio,

            "criteria":
                {
                    "min_length":
                        window_size,

                    "min_gc_content":
                        min_gc,

                    "min_oe_ratio":
                        min_oe
                }
        })

    return islands
