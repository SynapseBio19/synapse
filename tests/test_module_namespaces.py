"""The namespaces of analyses.blast and analyses.primer stay importable (Phase 7.1).

Phase 7.1 moves the implementation of both modules into the internal
packages ``analyses/_blast`` and ``analyses/_primer``; ``analyses/blast.py``
and ``analyses/primer.py`` remain as the public modules. Every name either
module defined at top level at commit 9077f1e, private helpers included
(the menus and tests use some of them), must still be reachable through
the public module.

One name is removed on purpose: ``analyses.primer._V5_ACTIVE_CONFIG``, the
hidden global configuration of the candidate screen. It is now passed
explicitly (``docs/CHANGES.md``, Phase 7.1).
"""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import textwrap
import unittest

PROJECT_ROOT = Path(__file__).resolve().parent.parent

BLAST_NAMES = (
    'NCBI_BLAST_URL', 'DEFAULT_TOOL', 'DEFAULT_EMAIL', 'API_TIMEOUT_SECONDS',
    'NCBI_MIN_REQUEST_INTERVAL_SECONDS', 'NCBI_MIN_POLL_INTERVAL_SECONDS',
    'DEFAULT_POLL_INTERVAL_SECONDS', 'DEFAULT_MAX_WAIT_SECONDS', 'DEFAULT_MAX_HITS',
    'DEFAULT_EVIDENCE_EVALUE_THRESHOLD', 'LOCAL_DB_ROOT_ENV', 'LOCAL_DB_MANIFEST_NAME',
    '_resolve_default_local_db_root', 'DEFAULT_LOCAL_DB_ROOT', 'LOCAL_DATABASE_CATALOG',
    'BLAST_EXECUTABLES', 'BLAST_UPDATE_SCRIPT', 'DNA_BASES', 'PROTEIN_AMINO_ACIDS',
    'AMBIGUOUS_DNA', 'PROGRAMS', 'QUERY_TYPES', 'programs_for_query_type',
    'check_program_query_type', 'logger', 'STATUS_RETRIES', '_is_transient_ncbi_error',
    'DATABASES', 'BLAST_TAG_ALIASES', 'BlastParameters', 'BlastHSP', 'BlastHit',
    'BlastQuery', 'BlastRunMetadata', 'BlastRun', 'HitEvidence', 'BlastComparison',
    'SelfTestResult', 'BlastAnalysisError', '_utc_now', '_evalue_text', '_local_tag',
    '_child_text', '_descendants', '_to_int', '_to_float', 'normalize_sequence',
    'infer_sequence_type', 'build_query', '_resolve_synapse_user_path',
    'load_fasta_queries', 'NCBIBlastClient', 'LocalBlastDatabase', 'BlastExecutableManager',
    '_parse_blastdbcmd_info', '_extract_blast_volume_index',
    '_parse_blast_alias_volume_names', 'LocalBlastDatabaseManager',
    '_build_local_database_command_reference', '_volume_prefix_from_stored_path',
    '_create_runtime_blast_alias', '_needs_runtime_blast_alias',
    '_build_local_blast_environment', '_should_run_short_query_rescue', 'LocalBlastEngine',
    '_parse_hsp', '_infer_scientific_name_from_definition', '_parse_hit',
    '_parse_parameters', '_parse_database_statistics', 'describe_run_from_xml',
    'parse_blast_xml', 'calculate_evidence_score', 'rank_hits', 'build_evidence_profile',
    'group_hits_by_taxonomy', 'alignment_text', 'extract_orf_context', 'build_orf_query',
    'run_blast_analysis', 'hsp_to_dict', 'hit_to_dict', 'run_to_dict', 'run_to_json',
    'export_json', 'export_csv', 'build_text_report', 'export_text', 'compare_blast_runs',
    'SAMPLE_XML', '_assert', '_test_sequence_validation', '_test_type_inference',
    '_test_program_database_compatibility', '_test_query_builder',
    '_test_fasta_query_loader', '_test_xml_parser', '_test_hit_ranking',
    '_test_evidence_profile', '_test_alignment', '_test_serialization', '_test_comparison',
    '_test_orf_adapter', '_test_export_schema', '_test_ncbi_client_configuration',
    '_test_ncbi_request_parameters', '_test_local_database_command_reference',
    '_test_local_blast_environment', '_test_local_blast_environment_preserves_existing',
    '_test_default_database_root_resolution', '_test_local_database_structure',
    '_test_nested_runtime_alias', '_test_blastdbcmd_info_parser',
    '_test_nested_local_database_discovery', '_test_nested_multivolume_database_discovery',
    '_test_alias_inside_first_volume_directory', '_test_local_database_manifest',
    '_test_local_parameter_compatibility', '_test_synapse_user_path_resolution',
    '_test_executable_manager', '_test_short_query_rescue_policy',
    '_test_short_query_task_configuration', '_test_taxonomy_definition_fallback',
    '_test_taxonomy_grouping_with_fallback', 'run_self_test',
)

PRIMER_NAMES = (
    'DNA_ALPHABET', 'IUPAC_ALPHABET', 'IUPAC_COMPLEMENT', 'NEAREST_NEIGHBOR_TM_METHODS',
    'NN_PARAMETERS', 'NN_INITIATION', 'NN_AT_TERMINAL', 'NN_GC_TERMINAL', 'NN_SYMMETRY_DS',
    'GAS_CONSTANT', 'PrimerDesignConfig', 'ThermodynamicResult', 'StructureResult',
    'SpecificityHit', 'PrimerCandidate', 'PrimerPair', 'PrimerDesignReport',
    'normalize_dna_sequence', '_normalize_dna_cached', 'validate_target_region',
    'reverse_complement', '_IUPAC_COMPLEMENT_TABLE', 'gc_percent', 'count_gc_clamp',
    'longest_homopolymer', 'longest_tandem_repeat', 'calculate_tm_wallace',
    'calculate_tm_salt_adjusted', 'calculate_tm', 'get_thermodynamic_conditions',
    '_alignment_stem', '_three_prime_complementarity', '_distance_score', '_category',
    '_candidate_rejections', 'score_primer_candidate', '_make_candidate',
    '_candidate_thermodynamics', '_make_candidate_measurements',
    'generate_primer_candidates', '_candidate_windows', '_generation_sort_key',
    'screen_primer_specificity', 'build_primer_pair', 'rank_primer_pairs', 'design_primers',
    'summarize_primer_pair', 'get_primer_summary', 'analyze_primer', 'get_amplicon',
    'run_primer_self_test', 'ABSOLUTE_ZERO_C', 'NEAREST_NEIGHBOR_PARAMETERS',
    'INITIATION_ENTHALPY', 'INITIATION_ENTROPY', 'TERMINAL_AT_ENTHALPY',
    'TERMINAL_AT_ENTROPY', 'SYMMETRY_ENTROPY', '_SANTALUCIA_2004_STACKS',
    'NN_PARAMETER_SETS', 'DEFAULT_NN_PARAMETER_SET', 'SALT_MODELS', '_nn_parameter_set',
    '_salt_phosphate_count', 'ThermodynamicConfig', 'DuplexThermodynamics',
    'PrimerThermodynamicProfile', 'validate_thermodynamic_conditions',
    'calculate_effective_ionic_strength', 'calculate_salt_correction',
    '_is_self_complementary', '_terminal_at_count', '_nearest_neighbor_sum',
    '_nearest_neighbor_sum_uncached', '_nearest_neighbor_sum_cached',
    'calculate_duplex_concentration_term', '_calculate_tm_nearest_neighbor_raw',
    'calculate_duplex_thermodynamics', '_duplex_thermodynamics',
    'calculate_terminal_stability', 'calculate_primer_thermodynamics',
    '_duplex_thermodynamics_cached', 'calculate_tm_advanced', 'AdvancedPrimerCandidate',
    'AdvancedPrimerPair', 'score_thermodynamic_profile', 'score_thermodynamic_pair',
    'analyze_primer_candidate_thermodynamics', 'analyze_primer_pair_thermodynamics',
    'ProfessionalPrimerConfig', 'ProfessionalPrimerDesignReport',
    'get_professional_primer_summary', 'analyze_primer_professional',
    'diagnose_thermodynamic_failure', 'compare_tm_methods', 'run_thermodynamic_self_test',
    'run_primer_v3_self_test', 'DEFAULT_PROFESSIONAL_THERMODYNAMIC_CONFIG',
    'HIGH_SALT_SCREENING_CONFIG', 'LOW_SALT_SCREENING_CONFIG', 'get_thermodynamic_summary',
    'PRIMER_V3_PUBLIC_API', 'NearestNeighborTmResult', 'calculate_tm_nearest_neighbor',
    'run_primer_v3_self_test_final', 'primer_pair_binding_sites',
    'run_primer_v4_screening_self_test', '_structure_penalty', '_screening_quality',
    '_measurement_config', '_apply_selected_tm_model',
    'generate_professional_primer_candidates', '_screening_tm', '_screening_score',
    '_screening_upper_bound', '_tm_tier', '_measure_candidates', '_measure_all_candidates',
    '_structure_severity', 'analyze_hairpin', 'analyze_self_dimer', 'analyze_heterodimer',
    '_pair_model_tm', '_model_tm_optimum', '_piecewise_penalty', '_pair_quality',
    '_sequence_similarity', '_is_diverse', '_select_diverse_pairs',
    'rank_professional_primer_pairs', 'summarize_professional_primer_pair',
    'design_primers_professional', '_professional_primer_report',
    'RestrictionSitePrimerDesign', 'restriction_site_choices',
    'design_primers_for_restriction_site',
)

# Phase 8: analyses.mutation, recorded from 727208a before its split.
MUTATION_NAMES = (
    'DNA_ALPHABET', 'IUPAC_CODES', 'IUPAC_FROM_BASES', 'IUPAC_COMPLEMENT', 'CODON_TABLE',
    'AA_NAMES', 'PURINES', 'PYRIMIDINES', 'SequencePromptCancelled',
    'MutationAnalysisError', 'AlignmentBlock', 'Mutation', 'MutationStatistics',
    'CodingAnalysis', 'MutationAnalysisResult', 'SequenceCoordinate', 'SequenceQuality',
    'MutationContext', 'MutationSpectrum', 'MutationHotspot', 'RestrictionSiteImpact',
    'MutationImpact', 'AdvancedMutationAnalysis', 'normalize_dna',
    'normalize_iupac_sequence', 'validate_dna_sequence', 'gc_content', 'reverse_complement',
    'reverse_complement_iupac', 'iupac_matches', 'iupac_code', 'translate_dna',
    'amino_acid_name', '_nucleotide_match', 'BANDED_ALIGNMENT_MIN_CELLS', '_banded_align',
    '_banded_is_exact', 'global_align', 'alignment_identity', '_full_edit_distance',
    '_banded_edit_distance', '_edit_distance', 'sequence_similarity',
    'build_alignment_blocks', 'is_transition', 'classify_substitution',
    '_mutation_notation', '_make_mutation', 'detect_mutations_from_alignment',
    '_detect_simple_indel', 'detect_mutations', 'classify_coding_substitution',
    '_translate_canonical', '_first_difference_position', 'analyze_coding_consequences',
    'annotate_mutations_with_coding_consequences', 'calculate_mutation_statistics',
    'compare_sequences', 'compare_multiple_variants', '_mutation_key', 'shared_mutations',
    'unique_mutations', 'annotate_regions', 'annotate_sequence_regions',
    'sequence_complexity', 'sequence_quality', 'extract_mutation_context',
    'build_mutation_spectrum', 'find_mutation_hotspots', 'normalize_indel_left',
    'MAX_RECOMMENDED_ALIGNMENT_LENGTH', 'CODING_FRAME_ASSUMPTION', 'assess_coding_frame',
    'orf_reference_cds', 'mutation_hgvs_c', '_reference_span', 'hgvs_c_allele',
    'mutation_hgvs_p', 'compare_proteins', '_find_all_occurrences',
    'restriction_site_impact', '_interval_overlap', 'primer_mutation_overlap',
    'crispr_mutation_overlap', 'validate_mutation_result', 'format_mutation',
    'build_mutation_report', 'mutation_signature', 'create_advanced_analysis',
    'advanced_mutation_report', 'mutation_map', 'result_to_json', 'run_mutation_self_test',
    'assert_mutation_self_test', 'analyze_mutations', 'mutation_report', '__all__',
)

REMOVED = {"analyses.primer": ("_V5_ACTIVE_CONFIG",)}

# Phase 8 (2026-11-18): the historical _vN_ helpers were renamed; the old
# names are gone on purpose (docs/CHANGES.md, Phase 8, package B).
RENAMED_PRIMER_NAMES = {
    '_v5_structure_penalty': '_structure_penalty',
    '_v5_candidate_quality': '_screening_quality',
    '_v6_relaxed_generation_config': '_measurement_config',
    '_v7_apply_selected_tm_model': '_apply_selected_tm_model',
    '_v7_generate_professional_primer_candidates': 'generate_professional_primer_candidates',
    '_V7_NEAREST_NEIGHBOR_METHODS': 'NEAREST_NEIGHBOR_TM_METHODS',
    '_v7_screening_tm': '_screening_tm',
    '_v7_screening_score': '_screening_score',
    '_v7_screening_upper_bound': '_screening_upper_bound',
    '_v7_tm_tier': '_tm_tier',
    '_v7_measure_candidates': '_measure_candidates',
    '_v7_measure_all_candidates': '_measure_all_candidates',
    '_v8_structure_severity': '_structure_severity',
    '_v11_model_tm': '_pair_model_tm',
    '_v11_tm_optimum': '_model_tm_optimum',
    '_v11_piecewise_penalty': '_piecewise_penalty',
    '_v11_pair_quality': '_pair_quality',
    '_v11_sequence_similarity': '_sequence_similarity',
    '_v11_is_diverse': '_is_diverse',
    '_v11_select_pairs': '_select_diverse_pairs',
    '_v11_rank_professional_primer_pairs': 'rank_professional_primer_pairs',
    '_calculate_tm_nearest_neighbor_v3_raw': '_calculate_tm_nearest_neighbor_raw',
}


class PublicModuleNamespaceTests(unittest.TestCase):
    def _check(self, module, names):
        missing = [name for name in names if not hasattr(module, name)]
        self.assertEqual(missing, [], f"{module.__name__} no longer provides these names")

    def test_blast_namespace(self):
        from analyses import blast
        self._check(blast, BLAST_NAMES)

    def test_primer_namespace(self):
        from analyses import primer
        self._check(primer, PRIMER_NAMES)

    def test_mutation_namespace(self):
        from analyses import mutation
        self._check(mutation, MUTATION_NAMES)

    def test_mutation_all_names_exist(self):
        from analyses import mutation
        self._check(mutation, mutation.__all__)

    def test_renamed_primer_helpers(self):
        from analyses import primer
        self._check(primer, tuple(RENAMED_PRIMER_NAMES.values()))
        stale = [old for old in RENAMED_PRIMER_NAMES if hasattr(primer, old)]
        self.assertEqual(stale, [])

    def test_the_declared_primer_api_exists(self):
        # PRIMER_V3_PUBLIC_API listed two report printers of ui/primer_menu.py
        # until 2026-11-18; the list now names only analyses.primer itself.
        from analyses import primer
        self._check(primer, primer.PRIMER_V3_PUBLIC_API)


class InternalModuleIndependenceTests(unittest.TestCase):
    """Each module of the internal packages is complete on its own.

    The public modules import everything, so a statement left behind in
    them (or a module that only works after another one ran) would go
    unnoticed in the rest of the suite. Each internal module is imported
    here in a clean interpreter state, without the public module.
    """

    PROBE = textwrap.dedent('''
        import importlib, json, pkgutil, sys
        import analyses._blast, analyses._mutation, analyses._primer
        failures = {}
        for package in (analyses._blast, analyses._primer, analyses._mutation):
            for info in pkgutil.iter_modules(package.__path__):
                name = package.__name__ + "." + info.name
                for loaded in [m for m in sys.modules if m.startswith("analyses")]:
                    del sys.modules[loaded]
                try:
                    importlib.import_module(name)
                except Exception as error:
                    failures[name] = repr(error)
                if {"analyses.blast", "analyses.primer", "analyses.mutation"} & set(sys.modules):
                    failures[name] = "imports the public module"
        for loaded in [m for m in sys.modules if m.startswith("analyses")]:
            del sys.modules[loaded]
        from analyses._primer import thermodynamics
        stacks = thermodynamics.NN_PARAMETER_SETS["santalucia_2004"]["stacks"]
        print(json.dumps({"failures": failures, "aa_2004": stacks["AA"], "aa_1998": thermodynamics.NEAREST_NEIGHBOR_PARAMETERS["AA"]}))
    ''')

    def test_every_internal_module_imports_alone(self):
        result = subprocess.run([sys.executable, "-c", self.PROBE], cwd=PROJECT_ROOT,
                                capture_output=True, text=True, check=True)
        report = json.loads(result.stdout.strip().splitlines()[-1])
        self.assertEqual(report["failures"], {})
        # The SantaLucia & Hicks (2004) AA/TT stack is set by a statement that
        # defines no name; it must run with the thermodynamics module itself.
        self.assertEqual(report["aa_2004"], [-7.6, -21.3])
        self.assertEqual(report["aa_1998"], [-7.9, -22.2])

    def test_public_names_are_the_implementation_objects(self):
        from analyses import blast, primer
        from analyses._blast import models, ncbi
        from analyses._primer import screening, thermodynamics
        self.assertIs(blast.BlastRun, models.BlastRun)
        self.assertIs(blast.NCBIBlastClient, ncbi.NCBIBlastClient)
        self.assertIs(primer.calculate_tm_nearest_neighbor, thermodynamics.calculate_tm_nearest_neighbor)
        self.assertIs(primer._nearest_neighbor_sum_cached, thermodynamics._nearest_neighbor_sum_cached)
        self.assertIs(primer._measure_candidates, screening._measure_candidates)


if __name__ == "__main__":
    unittest.main()
