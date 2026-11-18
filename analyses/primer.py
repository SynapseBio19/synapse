"""
Synapse — Professional Primer Design Engine v2
================================================

A self-contained primer-design module for PCR-style primer design.

v2 adds a documented nearest-neighbor thermodynamic model with configurable:
- monovalent salt (Na+) concentration
- Mg2+ concentration
- dNTP concentration
- total primer concentration
- primer self-complementarity symmetry correction

The default Tm model remains configurable. The recommended default is the
nearest-neighbor model. Wallace and empirical salt-adjusted screening remain
available for comparison.

Scientific scope and limitations:
- Nearest-neighbor Tm uses SantaLucia-style DNA/DNA thermodynamic parameters.
- Mixed-ion handling uses an effective monovalent ionic-strength approximation;
  it is an approximation, not a full unified-ion Owczarzy/Primer3 calculation.
- Mg2+ is adjusted for dNTP binding through a simple free-Mg approximation.
- Hairpin, self-dimer and heterodimer detection remain sequence-level heuristics;
  they do not calculate full secondary-structure partition functions or ΔG.
- Specificity screening is limited to reference sequences supplied by the caller.
  No genome-wide specificity claim is made without an actual database search.
- Experimental primer performance must be validated under the actual reaction
  conditions.

Coordinates:
- All reported sequence coordinates are 0-based, half-open [start, end).
"""

# Phase 7.1: the implementation moved to the internal package
# ``analyses._primer`` (one module per responsibility, see its __init__).
# This module remains the public interface and re-exports every name it
# defined before, private helpers included.

from __future__ import annotations

from core.textformat import major_rule
from analyses._primer.sequence import (
    count_gc_clamp,
    DNA_ALPHABET,
    gc_percent,
    IUPAC_ALPHABET,
    IUPAC_COMPLEMENT,
    _IUPAC_COMPLEMENT_TABLE,
    longest_homopolymer,
    longest_tandem_repeat,
    _normalize_dna_cached,
    normalize_dna_sequence,
    reverse_complement,
    validate_target_region,
)
from analyses._primer.models import (
    PrimerCandidate,
    PrimerDesignConfig,
    PrimerDesignReport,
    PrimerPair,
    SpecificityHit,
    StructureResult,
    ThermodynamicResult,
)
from analyses._primer.thermodynamics import (
    ABSOLUTE_ZERO_C,
    calculate_duplex_concentration_term,
    calculate_duplex_thermodynamics,
    calculate_effective_ionic_strength,
    calculate_primer_thermodynamics,
    calculate_salt_correction,
    calculate_terminal_stability,
    calculate_tm,
    calculate_tm_advanced,
    calculate_tm_nearest_neighbor,
    _calculate_tm_nearest_neighbor_v3_raw,
    calculate_tm_salt_adjusted,
    calculate_tm_wallace,
    compare_tm_methods,
    DEFAULT_NN_PARAMETER_SET,
    DEFAULT_PROFESSIONAL_THERMODYNAMIC_CONFIG,
    _duplex_thermodynamics,
    _duplex_thermodynamics_cached,
    DuplexThermodynamics,
    GAS_CONSTANT,
    get_thermodynamic_conditions,
    get_thermodynamic_summary,
    HIGH_SALT_SCREENING_CONFIG,
    INITIATION_ENTHALPY,
    INITIATION_ENTROPY,
    _is_self_complementary,
    LOW_SALT_SCREENING_CONFIG,
    NEAREST_NEIGHBOR_PARAMETERS,
    _nearest_neighbor_sum,
    _nearest_neighbor_sum_cached,
    _nearest_neighbor_sum_uncached,
    NEAREST_NEIGHBOR_TM_METHODS,
    NearestNeighborTmResult,
    NN_AT_TERMINAL,
    NN_GC_TERMINAL,
    NN_INITIATION,
    _nn_parameter_set,
    NN_PARAMETER_SETS,
    NN_PARAMETERS,
    NN_SYMMETRY_DS,
    PrimerThermodynamicProfile,
    SALT_MODELS,
    _salt_phosphate_count,
    _SANTALUCIA_2004_STACKS,
    SYMMETRY_ENTROPY,
    _terminal_at_count,
    TERMINAL_AT_ENTHALPY,
    TERMINAL_AT_ENTROPY,
    ThermodynamicConfig,
    validate_thermodynamic_conditions,
)
from analyses._primer.structure import (
    _alignment_stem,
    analyze_hairpin,
    analyze_heterodimer,
    analyze_self_dimer,
    _three_prime_complementarity,
    _v8_structure_severity,
)
from analyses._primer.specificity import primer_pair_binding_sites, screen_primer_specificity
from analyses._primer.design import (
    analyze_primer,
    build_primer_pair,
    _candidate_rejections,
    _candidate_thermodynamics,
    _candidate_windows,
    _category,
    design_primers,
    _distance_score,
    generate_primer_candidates,
    _generation_sort_key,
    get_amplicon,
    get_primer_summary,
    _make_candidate,
    _make_candidate_measurements,
    rank_primer_pairs,
    score_primer_candidate,
    summarize_primer_pair,
)
from analyses._primer.professional import (
    AdvancedPrimerCandidate,
    AdvancedPrimerPair,
    analyze_primer_candidate_thermodynamics,
    analyze_primer_pair_thermodynamics,
    analyze_primer_professional,
    diagnose_thermodynamic_failure,
    get_professional_primer_summary,
    ProfessionalPrimerConfig,
    ProfessionalPrimerDesignReport,
    score_thermodynamic_pair,
    score_thermodynamic_profile,
)
from analyses._primer.screening import (
    generate_professional_primer_candidates,
    _v5_candidate_quality,
    _v5_structure_penalty,
    _v6_relaxed_generation_config,
    _v7_apply_selected_tm_model,
    _v7_generate_professional_primer_candidates,
    _v7_measure_all_candidates,
    _v7_measure_candidates,
    _V7_NEAREST_NEIGHBOR_METHODS,
    _v7_screening_score,
    _v7_screening_tm,
    _v7_screening_upper_bound,
    _v7_tm_tier,
)
from analyses._primer.ranking import (
    rank_professional_primer_pairs,
    summarize_professional_primer_pair,
    _v11_is_diverse,
    _v11_model_tm,
    _v11_pair_quality,
    _v11_piecewise_penalty,
    _v11_rank_professional_primer_pairs,
    _v11_select_pairs,
    _v11_sequence_similarity,
    _v11_tm_optimum,
)
from analyses._primer.workflow import (
    design_primers_for_restriction_site,
    design_primers_professional,
    _professional_primer_report,
    restriction_site_choices,
    RestrictionSitePrimerDesign,
)
from analyses._primer.selftest import (
    PRIMER_V3_PUBLIC_API,
    run_primer_self_test,
    run_primer_v3_self_test,
    run_primer_v3_self_test_final,
    run_primer_v4_screening_self_test,
    run_thermodynamic_self_test,
)


# ============================================================================
# MODULE ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    from ui import style  # only when run as a script
    test_results = run_primer_self_test()
    style.header('PRIMER MODULE SELF-TEST')
    for test_name, passed in test_results.items():
        print(f"[{'PASS' if passed else 'FAIL'}] {test_name}")
    print(major_rule())
    print("ALL TESTS PASSED:", all(test_results.values()))


# ============================================================================
# V3 MODULE ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    from ui import style  # only when run as a script
    results = run_primer_v3_self_test()

    style.header('PROFESSIONAL PRIMER ENGINE V3 SELF-TEST')

    for name, passed in results.items():
        print(
            f"[{'PASS' if passed else 'FAIL'}] "
            f"{name}"
        )

    print(major_rule())

    print(
        "TOTAL TESTS:",
        len(results),
    )

    print(
        "ALL TESTS PASSED:",
        all(results.values()),
    )

    print(major_rule())
