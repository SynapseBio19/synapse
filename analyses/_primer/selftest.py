"""The built-in self-tests of the primer module and the list of the v3 API."""

from __future__ import annotations

from analyses._primer.sequence import (
    count_gc_clamp,
    gc_percent,
    longest_homopolymer,
    longest_tandem_repeat,
    normalize_dna_sequence,
    reverse_complement,
)
from analyses._primer.models import PrimerDesignConfig
from analyses._primer.thermodynamics import (
    calculate_duplex_thermodynamics,
    calculate_effective_ionic_strength,
    calculate_primer_thermodynamics,
    calculate_terminal_stability,
    calculate_tm_nearest_neighbor,
    calculate_tm_wallace,
    compare_tm_methods,
    NEAREST_NEIGHBOR_PARAMETERS,
    ThermodynamicConfig,
)
from analyses._primer.structure import analyze_hairpin, analyze_heterodimer, analyze_self_dimer
from analyses._primer.specificity import screen_primer_specificity
from analyses._primer.design import design_primers, get_amplicon, get_primer_summary
from analyses._primer.professional import (
    analyze_primer_professional,
    get_professional_primer_summary,
    ProfessionalPrimerConfig,
)
from analyses._primer.workflow import design_primers_professional


# ============================================================================
# REGRESSION TESTS
# ============================================================================

def run_primer_self_test() -> dict[str, bool]:
    """Run deterministic implementation and regression tests."""
    results = {}

    assert normalize_dna_sequence(" atgc\n") == "ATGC"
    results["sequence_normalization"] = True

    assert reverse_complement("ATGC") == "GCAT"
    results["reverse_complement"] = True

    assert gc_percent("ATGC") == 50.0
    results["gc_percent"] = True

    assert calculate_tm_wallace("ATGC") == 12.0
    results["wallace_tm"] = True

    assert count_gc_clamp("ATGCGCT") == 0
    assert count_gc_clamp("ATGCGC") == 4
    assert count_gc_clamp("ATGCGCTGC") == 2
    results["gc_clamp"] = True

    assert longest_homopolymer("AATTTGC") == 3
    results["homopolymer"] = True

    assert longest_tandem_repeat("ATATATGC") >= 3
    results["tandem_repeat"] = True

    thermo = calculate_tm_nearest_neighbor(
        "GCGATCGATCGC",
        salt_mM=50,
        mg_mM=1.5,
        dntp_mM=0.2,
        primer_concentration_nM=250,
        return_details=True,
    )
    assert thermo.model == "nearest_neighbor"
    assert thermo.tm > 0
    results["nearest_neighbor_tm"] = True

    thermo_high_salt = calculate_tm_nearest_neighbor(
        "GCGATCGATCGC",
        salt_mM=100,
        mg_mM=1.5,
        dntp_mM=0.2,
        primer_concentration_nM=250,
    )
    assert thermo_high_salt > thermo.tm
    results["salt_response"] = True

    thermo_low_conc = calculate_tm_nearest_neighbor(
        "GCGATCGATCGC",
        salt_mM=50,
        mg_mM=1.5,
        dntp_mM=0.2,
        primer_concentration_nM=25,
    )
    assert thermo_low_conc < thermo.tm
    results["concentration_response"] = True

    hairpin = analyze_hairpin("GCGAAAAACGC")
    assert hairpin.max_stem >= 3
    results["hairpin_detection"] = True

    dimer = analyze_self_dimer("GCGATCGC")
    assert dimer.max_stem >= 3
    results["self_dimer_detection"] = True

    heterodimer = analyze_heterodimer("GCGATCGC", "GCGATCGC")
    assert heterodimer.max_stem >= 3
    results["heterodimer_detection"] = True

    test_sequence = (
        "TTATCGGATCTGCCGTGCGTACGACCGACACTACGCGGTTGGCTAAATAATGCTGGCGCTCGTACTCCTCCAGCGTCGATCTTCTGGTAACTAAGACTGGGACGTAGCACGGGTTGAGAAGAGTGTCTCTACGCCCCAGAAAGCGCCGTACCCTGGCTCCTTCGGATAGAGAAACGCATTGCGTTTAGGTCCATTTTATGGTGCCCGCACTGGATCCTGCGCACAGTTCATTCACATAACTGAGTCACCCAGGCTACGGATCATTCGTAAATCTCAGCTTGGGACTACTGCTTTCCTAGA"
    )

    report = design_primers(test_sequence, tm_method="nearest_neighbor")
    assert report.forward_candidates_evaluated > 0
    assert report.reverse_candidates_evaluated > 0
    results["design_workflow"] = True

    assert report.tm_model == "nearest_neighbor"
    assert report.thermodynamic_conditions["mg_mM"] == 1.5
    results["thermodynamic_reporting"] = True

    if report.ranked_pairs:
        pair = report.ranked_pairs[0]
        assert pair.amplicon_length == pair.amplicon_end - pair.amplicon_start
        assert len(pair.amplicon_sequence) == pair.amplicon_length
        assert get_amplicon(test_sequence, pair.forward, pair.reverse) == pair.amplicon_sequence
        assert 0.0 <= pair.amplicon_gc_percent <= 100.0
    results["amplicon_invariants"] = True

    hits = screen_primer_specificity("ATGC", {"reference": "CCCATGCGGG"}, max_mismatches=0)
    assert any(hit.exact for hit in hits)
    results["specificity_screen"] = True

    summary = get_primer_summary(report)
    assert "pairs_evaluated" in summary
    assert "tm_model" in summary
    assert "thermodynamic_conditions" in summary
    results["structured_api"] = True

    return results


# ============================================================================
# THERMODYNAMIC SELF-TESTS
# ============================================================================

def run_thermodynamic_self_test() -> dict[str, bool]:
    """
    Deterministic tests for the V3 thermodynamic layer.
    """

    results = {}

    # ------------------------------------------------------------------
    # Condition validation
    # ------------------------------------------------------------------

    conditions = ThermodynamicConfig()

    conditions.validate()

    results[
        "thermodynamic_conditions"
    ] = True

    # ------------------------------------------------------------------
    # Effective salt
    # ------------------------------------------------------------------

    effective = (
        calculate_effective_ionic_strength(
            sodium_mM=50.0,
            magnesium_mM=1.5,
            dntp_mM=0.2,
        )
    )

    assert effective > 50.0

    results[
        "effective_ionic_strength"
    ] = True

    # ------------------------------------------------------------------
    # Nearest-neighbor parameter availability
    # ------------------------------------------------------------------

    for dinucleotide in (
        "AA",
        "AT",
        "TA",
        "CG",
        "GC",
        "GG",
    ):
        assert (
            dinucleotide
            in NEAREST_NEIGHBOR_PARAMETERS
        )

    results[
        "nearest_neighbor_parameters"
    ] = True

    # ------------------------------------------------------------------
    # Nearest-neighbor Tm
    # ------------------------------------------------------------------

    test_primer = (
        "CTGGTAACTAAGACTGGGAC"
    )

    nn_tm = (
        calculate_tm_nearest_neighbor(
            test_primer
        )
    )

    assert (
        0.0
        < nn_tm
        < 100.0
    )

    results[
        "nearest_neighbor_tm"
    ] = True

    # ------------------------------------------------------------------
    # Concentration response
    # ------------------------------------------------------------------

    low_concentration = (
        calculate_tm_nearest_neighbor(
            test_primer,
            primer_nM=25.0,
        )
    )

    high_concentration = (
        calculate_tm_nearest_neighbor(
            test_primer,
            primer_nM=1000.0,
        )
    )

    assert (
        high_concentration
        > low_concentration
    )

    results[
        "primer_concentration_response"
    ] = True

    # ------------------------------------------------------------------
    # Salt response
    # ------------------------------------------------------------------

    low_salt = (
        calculate_tm_nearest_neighbor(
            test_primer,
            sodium_mM=20.0,
        )
    )

    high_salt = (
        calculate_tm_nearest_neighbor(
            test_primer,
            sodium_mM=100.0,
        )
    )

    assert (
        high_salt
        > low_salt
    )

    results[
        "salt_response"
    ] = True

    # ------------------------------------------------------------------
    # Mg2+ / dNTP response
    # ------------------------------------------------------------------

    low_mg = (
        calculate_tm_nearest_neighbor(
            test_primer,
            magnesium_mM=0.5,
            dntp_mM=0.2,
        )
    )

    high_mg = (
        calculate_tm_nearest_neighbor(
            test_primer,
            magnesium_mM=2.5,
            dntp_mM=0.2,
        )
    )

    assert (
        high_mg
        > low_mg
    )

    results[
        "magnesium_response"
    ] = True

    # ------------------------------------------------------------------
    # dNTP effect
    # ------------------------------------------------------------------

    low_dntp = (
        calculate_effective_ionic_strength(
            sodium_mM=50.0,
            magnesium_mM=2.0,
            dntp_mM=0.1,
        )
    )

    high_dntp = (
        calculate_effective_ionic_strength(
            sodium_mM=50.0,
            magnesium_mM=2.0,
            dntp_mM=1.0,
        )
    )

    assert (
        high_dntp
        < low_dntp
    )

    results[
        "dntp_response"
    ] = True

    # ------------------------------------------------------------------
    # Full duplex thermodynamics
    # ------------------------------------------------------------------

    duplex = (
        calculate_duplex_thermodynamics(
            test_primer
        )
    )

    assert (
        duplex.sequence
        == test_primer
    )

    assert (
        duplex.delta_h_kcal_mol
        < 0.0
    )

    assert (
        duplex.delta_s_cal_mol_k
        < 0.0
    )

    results[
        "duplex_thermodynamics"
    ] = True

    # ------------------------------------------------------------------
    # Primer profile
    # ------------------------------------------------------------------

    profile = (
        calculate_primer_thermodynamics(
            test_primer
        )
    )

    assert (
        profile.method
        == "nearest_neighbor"
    )

    assert (
        profile.tm_C
        == duplex.tm_C
    )

    results[
        "primer_thermodynamic_profile"
    ] = True

    # ------------------------------------------------------------------
    # Terminal stability
    # ------------------------------------------------------------------

    terminal = (
        calculate_terminal_stability(
            test_primer
        )
    )

    assert (
        terminal["window"]
        > 0
    )

    results[
        "terminal_stability"
    ] = True

    # ------------------------------------------------------------------
    # Method comparison
    # ------------------------------------------------------------------

    comparison = (
        compare_tm_methods(
            test_primer
        )
    )

    assert (
        "nearest_neighbor_tm_C"
        in comparison
    )

    results[
        "tm_method_comparison"
    ] = True

    # ------------------------------------------------------------------
    # Professional configuration
    # ------------------------------------------------------------------

    professional_config = (
        ProfessionalPrimerConfig()
    )

    professional_config.validate()

    results[
        "professional_configuration"
    ] = True

    # ------------------------------------------------------------------
    # Professional workflow
    # ------------------------------------------------------------------

    test_sequence = (
        "TTATCGGATCTGCCGTGCGTACGACCGACACTACGCGGTT"
        "GGCTAAATAATGCTGGCGCTCGTACTCCTCCAGCGTCGATC"
        "TTCTGGTAACTAAGACTGGGACGTAGCACGGGTTGAGAAGA"
        "GTGTCTCTACGCCCCAGAAAGCGCCGTACCCTGGCTCCTTC"
        "GGATAGAGAAACGCATTGCGTTTAGGTCCATTTTATGGTGC"
        "CCGCACTGGATCCTGCGCACAGTTCATTCACATAACTGAGT"
        "CACCCAGGCTACGGATCATTCGTAAATCTCAGCTTGGGACTA"
        "CTGCTTTCCTAGA"
    )

    report = design_primers_professional(
        test_sequence,
        config=ProfessionalPrimerConfig(
            primer=PrimerDesignConfig(
                min_amplicon_length=80,
                max_amplicon_length=300,
                optimal_amplicon_length=150,
                min_tm=45.0,
                max_tm=80.0,
            ),
            thermodynamics=ThermodynamicConfig(
                sodium_mM=50.0,
                magnesium_mM=1.5,
                dntp_mM=0.2,
                primer_nM=250.0,
            ),
            tm_method="nearest_neighbor",
            top_pairs=5,
        ),
    )

    assert (
        report.sequence_length
        == len(test_sequence)
    )

    assert (
        report.forward_candidates_evaluated
        > 0
    )

    assert (
        report.reverse_candidates_evaluated
        > 0
    )

    results[
        "professional_workflow"
    ] = True

    # ------------------------------------------------------------------
    # Structured professional API
    # ------------------------------------------------------------------

    summary = (
        get_professional_primer_summary(
            report
        )
    )

    assert (
        "thermodynamic_conditions"
        in summary
    )

    assert (
        "tm_method"
        in summary
    )

    results[
        "professional_structured_api"
    ] = True

    return results


# ============================================================================
# V3 TEST ENTRY POINT
# ============================================================================

def run_primer_v3_self_test() -> dict[str, bool]:
    """
    Run both the original V2 regression suite and the V3 thermodynamic suite.

    The original suite is intentionally executed first so that V3 cannot
    silently break backward compatibility.
    """

    original = run_primer_self_test()

    advanced = run_thermodynamic_self_test()

    combined = {}

    for name, passed in original.items():
        combined[
            f"legacy_{name}"
        ] = passed

    for name, passed in advanced.items():
        combined[
            f"v3_{name}"
        ] = passed

    return combined


# ============================================================================
# PUBLIC API INVENTORY
# ============================================================================

PRIMER_V3_PUBLIC_API = (
    # Existing V2 API.
    "normalize_dna_sequence",
    "validate_target_region",
    "reverse_complement",
    "gc_percent",
    "count_gc_clamp",
    "longest_homopolymer",
    "longest_tandem_repeat",
    "calculate_tm_wallace",
    "calculate_tm_salt_adjusted",
    "calculate_tm",
    "analyze_hairpin",
    "analyze_self_dimer",
    "analyze_heterodimer",
    "generate_primer_candidates",
    "screen_primer_specificity",
    "build_primer_pair",
    "rank_primer_pairs",
    "design_primers",
    "summarize_primer_pair",
    "print_primer_report",
    "get_primer_summary",
    "analyze_primer",
    "get_amplicon",
    "run_primer_self_test",

    # V3 thermodynamic API.
    "ThermodynamicConfig",
    "DuplexThermodynamics",
    "PrimerThermodynamicProfile",
    "calculate_effective_ionic_strength",
    "calculate_salt_correction",
    "calculate_tm_nearest_neighbor",
    "calculate_duplex_thermodynamics",
    "calculate_duplex_concentration_term",
    "calculate_terminal_stability",
    "calculate_primer_thermodynamics",
    "calculate_tm_advanced",
    "ProfessionalPrimerConfig",
    "AdvancedPrimerCandidate",
    "AdvancedPrimerPair",
    "ProfessionalPrimerDesignReport",
    "generate_professional_primer_candidates",
    "rank_professional_primer_pairs",
    "design_primers_professional",
    "summarize_professional_primer_pair",
    "print_professional_primer_report",
    "get_professional_primer_summary",
    "analyze_primer_professional",
    "diagnose_thermodynamic_failure",
    "compare_tm_methods",
    "get_thermodynamic_summary",
    "run_thermodynamic_self_test",
    "run_primer_v3_self_test",
)


# ============================================================================
# FINAL V3 SELF-TEST ENTRY POINT
# ============================================================================

def run_primer_v3_self_test_final() -> dict[str, bool]:
    """
    Final compatibility suite.

    Runs the legacy suite after the compatibility adapter has been installed,
    then runs all V3 thermodynamic tests.
    """

    legacy = run_primer_self_test()
    v3 = run_thermodynamic_self_test()

    combined = {}

    for name, passed in legacy.items():
        combined[
            f"legacy_{name}"
        ] = passed

    for name, passed in v3.items():
        combined[
            f"v3_{name}"
        ] = passed

    return combined


# Superseded definition of _primer_ui_design() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# Superseded definition of run_primer_analysis_ui() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# ============================================================================
# SYNAPSE — PROFESSIONAL PRIMER ENGINE V4 SCREENING OVERRIDES
# ============================================================================
#
# The V2/V3 engine above is intentionally preserved.  The functions below
# replace only the professional candidate-selection layer and the CLI entry
# point.  The thermodynamic equations, legacy API and existing analysis
# routines remain available.
#
# Main improvements:
# - The professional NN workflow no longer pre-filters candidates with the
#   Wallace Tm model.
# - 3' GC clamp is treated as a scoring preference, not an automatic failure.
# - Hairpin/self-dimer heuristics are used as graded penalties; only clearly
#   severe structures are hard-rejected in the professional workflow.
# - Candidate search uses the complete configured length range before ranking.
# - If the strict 57–63 °C window contains no candidates, a transparent
#   fallback tier searches 55–65 °C and then reports the relaxation.
# - The report explains why candidates were rejected instead of simply saying
#   "0 candidates".
# - Existing legacy generate_primer_candidates() behavior is untouched.
# ============================================================================



# Superseded definition of generate_professional_primer_candidates() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# Superseded definition of design_primers_professional() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# ============================================================================
# V4 REGRESSION CHECKS
# ============================================================================

def run_primer_v4_screening_self_test() -> dict[str, bool]:
    """Test the improved professional screening layer."""
    results = {}

    sequence = (
        "ATGCGTACGTTAGCTAGCGATCGTACGATCGGCTAGCTAGG"
        "ATCGCTAGCTGATCGATCGGATCGCTAGCTAGCTACGATCG"
        "ATCGGCTAGCTAGCTGATCGATCGGCTAGCTAGCTAGCTAA"
        "GATCGCTAGCTAGCTGATCGATCGGATCGCTAGCTAGCTGA"
    )

    config = ProfessionalPrimerConfig(
        primer=PrimerDesignConfig(
            min_amplicon_length=40,
            max_amplicon_length=200,
            optimal_amplicon_length=120,
        ),
        thermodynamics=ThermodynamicConfig(),
        tm_method="nearest_neighbor",
        top_pairs=5,
    )

    report = design_primers_professional(
        sequence,
        config=config,
    )

    assert report.forward_candidates_evaluated > 0
    assert report.reverse_candidates_evaluated > 0
    results["professional_candidate_generation"] = True

    if report.ranked_pairs:
        best = report.ranked_pairs[0]
        assert best.pair.forward.orientation == "FORWARD"
        assert best.pair.reverse.orientation == "REVERSE"
        assert best.pair.amplicon_length > 0
    results["pair_invariants"] = True

    primer = "ATGCGTACGTTAGCTAGCGA"
    analysis = analyze_primer_professional(
        primer,
        conditions=ThermodynamicConfig(),
    )
    assert analysis["length"] == 20
    assert 0.0 <= analysis["gc_percent"] <= 100.0
    results["single_primer_analysis"] = True

    return results
