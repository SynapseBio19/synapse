"""Split plan for analyses/primer.py (Phase 7.1, package C)."""

SOURCE = "analyses/primer.py"
PACKAGE = "analyses._primer"
PACKAGE_DIR = "analyses/_primer"
PACKAGE_DOC = """Implementation of :mod:`analyses.primer` (split out in Phase 7.1).

Import from ``analyses.primer``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``sequence``: DNA normalization, reverse complement, GC, repeats
- ``models``: the legacy design configuration and result dataclasses
- ``thermodynamics``: nearest-neighbor model, salt correction, Tm methods,
  condition presets
- ``structure``: hairpin, self-dimer and heterodimer screens
- ``specificity``: specificity screen against supplied references, binding sites
- ``design``: the legacy design engine (``design_primers``, ``analyze_primer``)
- ``professional``: the professional data model, thermodynamic scores,
  single-primer analysis and diagnostics
- ``screening``: the professional candidate pool and its Tm-tier screen
- ``ranking``: professional pair scoring, diversity selection, ranking
- ``workflow``: ``design_primers_professional`` and restriction-site design
- ``selftest``: the built-in self-tests and ``PRIMER_V3_PUBLIC_API``

The ``_v5``-``_v11`` prefixes are history, not versions: every such
function is live code (see the section comments).
"""

MODULES = {
    "sequence": ("DNA sequence helpers shared by every part of the primer module.", [
        "DNA_ALPHABET", "IUPAC_ALPHABET", "IUPAC_COMPLEMENT", "normalize_dna_sequence",
        "_normalize_dna_cached", "validate_target_region", "reverse_complement",
        "_IUPAC_COMPLEMENT_TABLE", "gc_percent", "count_gc_clamp", "longest_homopolymer",
        "longest_tandem_repeat",
    ]),
    "models": ("Configuration and result dataclasses of the legacy primer-design engine.", [
        "PrimerDesignConfig", "ThermodynamicResult", "StructureResult", "SpecificityHit",
        "PrimerCandidate", "PrimerPair", "PrimerDesignReport",
    ]),
    "thermodynamics": ("Melting temperature and duplex thermodynamics.\n\n"
                       "Wallace and salt-adjusted screening formulas, the SantaLucia\n"
                       "nearest-neighbor model with its parameter sets, salt and concentration\n"
                       "terms, the backward-compatible ``calculate_tm_nearest_neighbor`` and the\n"
                       "condition presets.", [
        "NEAREST_NEIGHBOR_TM_METHODS", "NN_PARAMETERS", "NN_INITIATION", "NN_AT_TERMINAL",
        "NN_GC_TERMINAL", "NN_SYMMETRY_DS", "GAS_CONSTANT", "calculate_tm_wallace",
        "calculate_tm_salt_adjusted", "calculate_tm", "get_thermodynamic_conditions",
        "ABSOLUTE_ZERO_C", "NEAREST_NEIGHBOR_PARAMETERS", "INITIATION_ENTHALPY",
        "INITIATION_ENTROPY", "TERMINAL_AT_ENTHALPY", "TERMINAL_AT_ENTROPY", "SYMMETRY_ENTROPY",
        "_SANTALUCIA_2004_STACKS", "NN_PARAMETER_SETS", "DEFAULT_NN_PARAMETER_SET", "SALT_MODELS",
        "_nn_parameter_set", "_salt_phosphate_count", "ThermodynamicConfig",
        "DuplexThermodynamics", "PrimerThermodynamicProfile", "validate_thermodynamic_conditions",
        "calculate_effective_ionic_strength", "calculate_salt_correction", "_is_self_complementary",
        "_terminal_at_count", "_nearest_neighbor_sum", "_nearest_neighbor_sum_uncached",
        "_nearest_neighbor_sum_cached", "calculate_duplex_concentration_term",
        "_calculate_tm_nearest_neighbor_v3_raw", "calculate_duplex_thermodynamics",
        "_duplex_thermodynamics", "calculate_terminal_stability", "calculate_primer_thermodynamics",
        "_duplex_thermodynamics_cached", "calculate_tm_advanced", "compare_tm_methods",
        "DEFAULT_PROFESSIONAL_THERMODYNAMIC_CONFIG", "HIGH_SALT_SCREENING_CONFIG",
        "LOW_SALT_SCREENING_CONFIG", "get_thermodynamic_summary", "NearestNeighborTmResult",
        "calculate_tm_nearest_neighbor",
    ]),
    "structure": ("Sequence-level hairpin, self-dimer and heterodimer screens (heuristics,\n"
                  "not free-energy calculations).", [
        "_alignment_stem", "_three_prime_complementarity", "_v8_structure_severity",
        "analyze_hairpin", "analyze_self_dimer", "analyze_heterodimer",
    ]),
    "specificity": ("Specificity screen against caller-supplied reference sequences and the\n"
                    "binding sites of a pair. No genome-wide claim is made.", [
        "screen_primer_specificity", "primer_pair_binding_sites",
    ]),
    "design": ("The legacy primer-design engine: candidate scoring and generation, pair\n"
               "building and ranking, ``design_primers`` and ``analyze_primer``.", [
        "_distance_score", "_category", "_candidate_rejections", "score_primer_candidate",
        "_make_candidate", "_candidate_thermodynamics", "_make_candidate_measurements",
        "generate_primer_candidates", "_candidate_windows", "_generation_sort_key",
        "build_primer_pair", "rank_primer_pairs", "design_primers", "summarize_primer_pair",
        "get_primer_summary", "analyze_primer", "get_amplicon",
    ]),
    "professional": ("The professional primer model: advanced candidates and pairs, their\n"
                     "thermodynamic scores, the configuration and report, single-primer\n"
                     "analysis and the failure diagnosis.", [
        "AdvancedPrimerCandidate", "AdvancedPrimerPair", "score_thermodynamic_profile",
        "score_thermodynamic_pair", "analyze_primer_candidate_thermodynamics",
        "analyze_primer_pair_thermodynamics", "ProfessionalPrimerConfig",
        "ProfessionalPrimerDesignReport", "get_professional_primer_summary",
        "analyze_primer_professional", "diagnose_thermodynamic_failure",
    ]),
    "screening": ("The professional candidate pool: measurement of every window, the\n"
                  "screening score and the Tm-tier screen with its upper bound.", [
        "_v5_structure_penalty", "_v5_candidate_quality", "_v6_relaxed_generation_config",
        "_v7_apply_selected_tm_model", "_v7_generate_professional_primer_candidates",
        "_V7_NEAREST_NEIGHBOR_METHODS", "_v7_screening_tm", "_v7_screening_score",
        "_v7_screening_upper_bound", "_v7_tm_tier", "_v7_measure_candidates",
        "_v7_measure_all_candidates", "generate_professional_primer_candidates",
    ]),
    "ranking": ("Professional pair scoring, diversity selection and ranking.", [
        "_v11_model_tm", "_v11_tm_optimum", "_v11_piecewise_penalty", "_v11_pair_quality",
        "_v11_sequence_similarity", "_v11_is_diverse", "_v11_select_pairs",
        "_v11_rank_professional_primer_pairs", "rank_professional_primer_pairs",
        "summarize_professional_primer_pair",
    ]),
    "workflow": ("Professional primer design for a target region and for a restriction site.", [
        "design_primers_professional", "_professional_primer_report",
        "RestrictionSitePrimerDesign", "restriction_site_choices",
        "design_primers_for_restriction_site",
    ]),
    "selftest": ("The built-in self-tests of the primer module and the list of the v3 API.", [
        "run_primer_self_test", "run_thermodynamic_self_test", "run_primer_v3_self_test",
        "PRIMER_V3_PUBLIC_API", "run_primer_v3_self_test_final", "run_primer_v4_screening_self_test",
    ]),
}

REWRITES = []
EXTRA = {}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 7.1: the implementation moved to the internal package\n"
    "# ``analyses._primer`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface and re-exports every name it\n"
    "# defined before, private helpers included.\n"
)
