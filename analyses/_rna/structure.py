"""The RNA hairpin screen (heuristic) and its summary."""

from __future__ import annotations

from analyses._rna.bases import RNA_STANDARD_BASES
from analyses._rna.composition import find_rna_homopolymers


# ============================================================
# RNA HAIRPIN SCREENING
# ============================================================

_WATSON_CRICK, _WOBBLE, _MISMATCH = 1, 2, 3


def _pair_classes(allow_gu):
    """Class of every pair of standard bases (left + right); others are absent."""
    classes = {}
    for left in sorted(RNA_STANDARD_BASES):
        for right in sorted(RNA_STANDARD_BASES):
            if left + right in ("AU", "UA", "CG", "GC"):
                classes[left + right] = _WATSON_CRICK
            elif allow_gu and {left, right} == {"G", "U"}:
                classes[left + right] = _WOBBLE
            else:
                classes[left + right] = _MISMATCH
    return classes


def _hairpin_candidates(sequence, min_stem, min_loop, max_loop, pair_classes):
    """Yield (stem start, stem length, loop length, paired, possible, wobble) per candidate.

    Unordered; find_rna_hairpins sorts them. See find_rna_hairpins for the
    extension from the innermost pair.
    """

    length = len(sequence)
    max_stem = min_stem + 3
    global_limit = max_stem // 4

    for inner in range(length):

        for loop_length in range(
            min_loop,
            max_loop + 1
        ):

            partner = inner + loop_length + 1

            if partner >= length:
                break

            paired = 0
            wobble = 0
            possible_pairs = 0
            mismatches = 0

            for offset in range(max_stem):

                left_index = inner - offset
                right_index = partner + offset

                if left_index < 0 or right_index >= length:
                    break

                pair_class = pair_classes.get(
                    sequence[left_index]
                    + sequence[right_index]
                )

                if pair_class is not None:

                    possible_pairs += 1

                    if pair_class == _WATSON_CRICK:
                        paired += 1
                    elif pair_class == _WOBBLE:
                        paired += 1
                        wobble += 1
                    else:
                        mismatches += 1
                        if mismatches > global_limit:
                            break

                stem_length = offset + 1

                if (
                    stem_length >= min_stem
                    and possible_pairs
                    and 4 * mismatches <= possible_pairs
                ):

                    yield ((
                        left_index,
                        stem_length,
                        loop_length,
                        paired,
                        possible_pairs,
                        wobble,
                    ))


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

    pair_classes = _pair_classes(allow_gu)

    if (
        min_stem < 1
        or min_loop < 1
        or max_loop < min_loop
    ):
        return structures

    # Every candidate is a stem (stem_start, stem_length) around a loop. All
    # candidates with the same innermost pair (left base i, loop length)
    # share their inner pairs, so the stem is extended outward from that
    # pair once and each stem length min_stem..min_stem+3 is read off the
    # running counts. A candidate needs paired / possible >= 0.75, i.e. at
    # most stem_length / 4 mismatches, so the extension stops at the first
    # mismatch beyond (min_stem + 3) // 4. The counts per candidate are
    # those of pairing it in full, so the candidates are the same; they are
    # then put into the order of the former loops (stem start, stem length,
    # loop length) before the final sort. This makes a 300 kb sequence about
    # five times faster (docs/BENCHMARKS.md).
    found = list(_hairpin_candidates(sequence, min_stem, min_loop, max_loop, pair_classes))

    # The order of the former nested loops: stem start, stem length, loop.
    found.sort()

    for stem_start, stem_length, loop_length, paired, possible_pairs, wobble in found:

        left_end = stem_start + stem_length
        right_start = left_end + loop_length
        right_end = right_start + stem_length

        pairing_fraction = (
            paired
            / possible_pairs
        )

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
                sequence[stem_start:left_end],

            "right_stem":
                sequence[right_start:right_end]
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


def count_rna_hairpins(
    rna,
    min_stem=4,
    min_loop=3,
    max_loop=12,
    allow_gu=False,
):
    """len(find_rna_hairpins(...)) without building the candidates.

    Every candidate is a distinct (start, stem, loop), so the count equals
    the length of the de-duplicated list.
    """

    if (
        min_stem < 1
        or min_loop < 1
        or max_loop < min_loop
    ):
        return 0

    return sum(1 for _ in _hairpin_candidates(
        rna.upper(), min_stem, min_loop, max_loop, _pair_classes(allow_gu)))


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
            count_rna_hairpins(rna, allow_gu=True),

        "homopolymers":
            homopolymers,

        "homopolymer_count":
            len(homopolymers)
    }
