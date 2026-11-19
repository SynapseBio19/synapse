
# Phase 9: the implementation moved to the internal package
# ``analyses._sequence`` (one module per responsibility, see its __init__).
# This module remains the public interface and re-exports every name it
# defined before.

from __future__ import annotations

from core.textformat import minor_rule
from analyses._sequence.bases import (
    AVERAGE_DSDNA_BP_MASS,
    count_ambiguous_bases,
    count_nucleotides,
    DNA_AMBIGUOUS_CODES,
    DNA_COMPLEMENT,
    DNA_IUPAC_BASES,
    DNA_NUCLEOTIDE_RESIDUE_MASSES,
    DNA_STANDARD_BASES,
    get_complement,
    get_reverse_complement,
    _get_standard_bases,
    _get_valid_iupac_bases,
    _is_iupac_base,
    _is_standard_base,
    normalize_dna_sequence,
)
from analyses._sequence.composition import (
    analyze_sequence_quality,
    calculate_at_content,
    calculate_at_skew,
    calculate_average_dna_molecular_weight,
    calculate_dna_content_ranges,
    calculate_dna_molecular_weight,
    calculate_gc_content,
    _calculate_gc_for_window,
    calculate_gc_skew,
    calculate_nucleotide_percentages,
    calculate_sequence_statistics,
)
from analyses._sequence.profiles import (
    analyze_complexity_profile,
    analyze_gc_profile,
    calculate_complexity_profile,
    calculate_gc_content_profile,
    calculate_profile_statistics,
    calculate_sequence_complexity,
    describe_profile_availability,
    find_homopolymers,
)
from analyses._sequence.cpg import analyze_cpg, DEFAULT_CPG_ISLAND_CRITERIA, find_cpg_islands
from analyses._sequence.summary import (
    analyze_sequence,
    create_sequence_interpretation,
    create_sequence_summary,
    get_sequence_forms,
)
from analyses._sequence.selftest import run_sequence_self_test


if __name__ == "__main__":
    from ui import style  # only when run as a script

    self_test = run_sequence_self_test()

    style.header('Sequence engine self-test')

    for name, passed in self_test.items():

        print(
            f"{name:<30}"
            f"{'PASS' if passed else 'FAIL'}"
        )

    print(minor_rule())

    if all(
        self_test.values()
    ):

        print(
            "RESULT: ALL TESTS PASSED"
        )

    else:

        print(
            "RESULT: ONE OR MORE TESTS FAILED"
        )