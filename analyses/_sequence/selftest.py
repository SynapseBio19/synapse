"""The built-in self-test of the sequence overview (menu entry and test suite)."""

from __future__ import annotations

from analyses._sequence.bases import count_nucleotides, get_complement, get_reverse_complement
from analyses._sequence.composition import (
    calculate_at_content,
    calculate_at_skew,
    calculate_average_dna_molecular_weight,
    calculate_dna_content_ranges,
    calculate_gc_content,
    calculate_gc_skew,
)
from analyses._sequence.profiles import (
    calculate_sequence_complexity,
    describe_profile_availability,
    find_homopolymers,
)
from analyses._sequence.cpg import analyze_cpg
from analyses._sequence.summary import analyze_sequence


# ============================================================
# SEQUENCE ENGINE SELF-TEST
# ============================================================

def run_sequence_self_test():

    results = {}

    test = "ATGC"

    # Composition.
    counts = count_nucleotides(test)
    results["nucleotide_counts"] = (
        counts == {
            "A": 1,
            "T": 1,
            "C": 1,
            "G": 1
        }
    )

    # GC / AT.
    results["gc_content"] = (
        abs(
            calculate_gc_content(test)
            - 50.0
        )
        < 1e-12
    )

    results["at_content"] = (
        abs(
            calculate_at_content(test)
            - 50.0
        )
        < 1e-12
    )

    # Skew.
    results["gc_skew"] = (
        abs(
            calculate_gc_skew(test)
        )
        < 1e-12
    )

    results["at_skew"] = (
        abs(
            calculate_at_skew(test)
        )
        < 1e-12
    )

    # Complexity: all four bases equally represented => 1.0.
    results["complexity"] = (
        abs(
            calculate_sequence_complexity(test)
            - 1.0
        )
        < 1e-12
    )

    # Complement.
    results["complement"] = (
        get_complement(test)
        == "TACG"
    )

    results["reverse_complement"] = (
        get_reverse_complement(test)
        == "GCAT"
    )

    # CpG.
    cpg_test = "CG"
    cpg = analyze_cpg(cpg_test)

    results["cpg_count"] = (
        cpg["count"] == 1
    )

    # Homopolymer.
    homopolymer = find_homopolymers(
        "AAAA"
    )

    results["homopolymer"] = (
        len(homopolymer) == 1
        and homopolymer[0]["length"] == 4
        and homopolymer[0]["start_1_based"] == 1
        and homopolymer[0]["end_1_based"] == 4
    )

    # IUPAC range.
    ranges = calculate_dna_content_ranges(
        "AN"
    )

    results["iupac_range"] = (
        ranges["gc_min"] == 0.0
        and ranges["gc_max"] == 50.0
        and ranges["at_min"] == 50.0
        and ranges["at_max"] == 100.0
    )

    # Profile availability.
    profile_info = describe_profile_availability(
        "ATGC",
        100,
        "GC content profile"
    )

    results["short_profile"] = (
        profile_info["available"] is False
        and profile_info["reason"]
        == "SEQUENCE_SHORTER_THAN_WINDOW"
    )

    # Conventional dsDNA estimate.
    results["average_ds_weight"] = (
        calculate_average_dna_molecular_weight(
            "ATGC"
        )
        == 2640.0
    )

    # Full pipeline.
    analysis = analyze_sequence(
        "ATGCGCGT"
    )

    results["full_pipeline"] = (
        analysis["sequence_length"] == 8
        and "interpretation" in analysis
        and "sequence_statistics" in analysis
        and "sequence_forms" in analysis
        and "cpg_island_criteria" in analysis
    )

    return results
