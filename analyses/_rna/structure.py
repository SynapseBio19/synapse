"""The RNA hairpin screen (heuristic) and its summary."""

from __future__ import annotations

from analyses._rna.bases import RNA_STANDARD_BASES
from analyses._rna.composition import find_rna_homopolymers


# ============================================================
# RNA HAIRPIN SCREENING
# ============================================================

def find_rna_hairpins(
    rna,
    min_stem=4,
    min_loop=3,
    max_loop=12,
    allow_gu=False,
):
    """Candidate stem-loops: stems of min_stem..min_stem+3 nt around a loop.

    A candidate is reported when at least 75 % of the stem positions with
    standard bases pair. Watson-Crick pairs (A-U, G-C) count; with
    allow_gu=True the G-U wobble pair, a regular but weaker pair in RNA
    helices, counts as well and "wobble_pairs" says how many a stem uses.
    The default stays Watson-Crick only: without free energies, counting
    G-U like G-C multiplies the candidates (58 -> 260 on the golden-master
    reference) instead of finding stronger stems, so the structure summary
    reports the G-U count separately. A sequence screen, not a prediction.
    """

    structures = []

    sequence = rna.upper()

    if (
        min_stem < 1
        or min_loop < 1
        or max_loop < min_loop
    ):
        return structures

    for stem_start in range(
        len(sequence)
    ):

        for stem_length in range(
            min_stem,
            min_stem + 4
        ):

            left_end = (
                stem_start
                + stem_length
            )

            if left_end >= len(sequence):
                continue

            for loop_length in range(
                min_loop,
                max_loop + 1
            ):

                right_start = (
                    left_end
                    + loop_length
                )

                right_end = (
                    right_start
                    + stem_length
                )

                if right_end > len(sequence):
                    continue

                left = sequence[
                    stem_start:left_end
                ]

                right = sequence[
                    right_start:right_end
                ]

                paired = 0
                wobble = 0
                possible_pairs = 0

                for left_base, right_base in zip(
                    left,
                    reversed(right)
                ):

                    if (
                        left_base
                        not in RNA_STANDARD_BASES
                        or right_base
                        not in RNA_STANDARD_BASES
                    ):
                        continue

                    possible_pairs += 1

                    if (
                        (
                            left_base == "A"
                            and right_base == "U"
                        )
                        or
                        (
                            left_base == "U"
                            and right_base == "A"
                        )
                        or
                        (
                            left_base == "C"
                            and right_base == "G"
                        )
                        or
                        (
                            left_base == "G"
                            and right_base == "C"
                        )
                    ):

                        paired += 1

                    elif allow_gu and {left_base, right_base} == {"G", "U"}:

                        paired += 1
                        wobble += 1

                if possible_pairs == 0:
                    continue

                pairing_fraction = (
                    paired
                    / possible_pairs
                )

                if pairing_fraction >= 0.75:

                    structures.append({
                        "start":
                            stem_start + 1,

                        "end":
                            right_end,

                        "stem_length":
                            stem_length,

                        "loop_length":
                            loop_length,

                        "pairing_fraction":
                            pairing_fraction,

                        "wobble_pairs":
                            wobble,

                        "left_stem":
                            left,

                        "right_stem":
                            right
                    })

    unique = []
    seen = set()

    for structure in structures:

        key = (
            structure["start"],
            structure["end"],
            structure["stem_length"],
            structure["loop_length"]
        )

        if key in seen:
            continue

        seen.add(key)
        unique.append(
            structure
        )

    unique.sort(
        key=lambda item: (
            item["pairing_fraction"],
            item["stem_length"]
        ),
        reverse=True
    )

    return unique


# ============================================================
# RNA STRUCTURAL SUMMARY
# ============================================================

def summarize_rna_structure(rna):

    hairpins = find_rna_hairpins(
        rna
    )

    homopolymers = (
        find_rna_homopolymers(
            rna
        )
    )

    return {
        "hairpins":
            hairpins,

        "hairpin_count":
            len(hairpins),

        # Candidates when G-U wobble pairs count as pairs (P2-5); reported
        # beside the Watson-Crick screen, not instead of it.
        "hairpin_count_with_gu":
            len(find_rna_hairpins(rna, allow_gu=True)),

        "homopolymers":
            homopolymers,

        "homopolymer_count":
            len(homopolymers)
    }
