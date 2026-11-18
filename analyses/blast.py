"""
Synapse - BLAST Search Engine
=============================

Professional BLAST integration and result-analysis module for Project Synapse.

Scientific design
-----------------
This module does NOT reimplement the BLAST algorithm and does not fabricate
database hits. It uses the NCBI BLAST URL API for remote searches and the official
NCBI BLAST+ command-line applications for local searches. Synapse performs
transparent post-processing and database management around those engines.

Supported BLAST programs
------------------------
- blastn   : nucleotide query -> nucleotide database
- blastp   : protein query   -> protein database
- blastx   : translated nucleotide query -> protein database
- tblastn  : protein query -> translated nucleotide database
- tblastx  : translated nucleotide query -> translated nucleotide database

Core capabilities
-----------------
- NCBI remote BLAST submission and polling.
- Local NCBI BLAST+ execution against managed databases.
- Local BLAST database status, verification, download/update and removal.
- Custom FASTA -> local BLAST database creation through makeblastdb.
- XML result parsing into typed Synapse records.
- Identity, coverage, E-value, bit score and alignment metrics.
- HSP-level alignment inspection.
- Transparent hit ranking / evidence profiling.
- Taxonomic grouping when taxonomy information is returned by BLAST.
- Reproducible run metadata and parameter recording.
- ORF-aware query preparation.
- Comparison of independent BLAST runs.
- JSON, CSV and text report export.
- Deterministic self-test suite that never calls the network.

Scientific boundaries
---------------------
- A BLAST hit is a sequence-similarity result, not proof of biological identity,
  function, pathogenicity, clinical relevance or evolutionary relationship.
- Interpretation is deliberately evidence-based and transparent.
- NCBI database contents and versions change over time; run metadata therefore
  records the requested database, program, parameters and timestamp.
- Remote BLAST requires internet access and is subject to NCBI service limits.
- Local BLAST requires an installed NCBI BLAST+ package and an installed local
  database.
- Official NCBI database downloads are explicit user actions; Synapse never
  downloads multi-gigabyte databases automatically.
- This module does not claim genome-wide biological conclusions beyond the
  searched database and returned alignments.

No third-party Python packages are required.
"""

# Phase 7.1: the implementation moved to the internal package
# ``analyses._blast`` (one module per responsibility, see its __init__).
# This module remains the public interface and re-exports every name it
# defined before, private helpers included.

from __future__ import annotations

from analyses._blast.paths import _resolve_synapse_user_path
from analyses._blast.config import (
    AMBIGUOUS_DNA,
    API_TIMEOUT_SECONDS,
    BLAST_EXECUTABLES,
    BLAST_TAG_ALIASES,
    BLAST_UPDATE_SCRIPT,
    DATABASES,
    DEFAULT_EMAIL,
    DEFAULT_EVIDENCE_EVALUE_THRESHOLD,
    DEFAULT_LOCAL_DB_ROOT,
    DEFAULT_MAX_HITS,
    DEFAULT_MAX_WAIT_SECONDS,
    DEFAULT_POLL_INTERVAL_SECONDS,
    DEFAULT_TOOL,
    DNA_BASES,
    LOCAL_DATABASE_CATALOG,
    LOCAL_DB_MANIFEST_NAME,
    LOCAL_DB_ROOT_ENV,
    NCBI_BLAST_URL,
    NCBI_MIN_POLL_INTERVAL_SECONDS,
    NCBI_MIN_REQUEST_INTERVAL_SECONDS,
    PROGRAMS,
    PROTEIN_AMINO_ACIDS,
    QUERY_TYPES,
    _resolve_default_local_db_root,
)
from analyses._blast.models import (
    BlastAnalysisError,
    BlastComparison,
    BlastHit,
    BlastHSP,
    BlastParameters,
    BlastQuery,
    BlastRun,
    BlastRunMetadata,
    _evalue_text,
    HitEvidence,
    LocalBlastDatabase,
    rank_hits,
    SelfTestResult,
    _utc_now,
)
from analyses._blast.query import (
    build_orf_query,
    build_query,
    check_program_query_type,
    extract_orf_context,
    infer_sequence_type,
    load_fasta_queries,
    normalize_sequence,
    programs_for_query_type,
)
from analyses._blast.parsing import (
    _child_text,
    _descendants,
    describe_run_from_xml,
    _infer_scientific_name_from_definition,
    _local_tag,
    parse_blast_xml,
    _parse_database_statistics,
    _parse_hit,
    _parse_hsp,
    _parse_parameters,
    _to_float,
    _to_int,
)
from analyses._blast.ncbi import _is_transient_ncbi_error, logger, NCBIBlastClient, STATUS_RETRIES
from analyses._blast.executables import BlastExecutableManager
from analyses._blast.ncbi_metadata import (
    check_download_space,
    DownloadSpaceCheck,
    fetch_database_metadata,
    format_bytes,
    free_disk_bytes,
    local_release_date,
    metadata_url,
    NCBI_DB_METADATA_URL,
    NcbiDatabaseMetadata,
    newer_release_available,
    parse_database_metadata,
)
from analyses._blast.runtime import (
    _build_local_blast_environment,
    _build_local_database_command_reference,
    _create_runtime_blast_alias,
    _needs_runtime_blast_alias,
    _volume_prefix_from_stored_path,
)
from analyses._blast.databases import (
    _extract_blast_volume_index,
    LocalBlastDatabaseManager,
    _parse_blast_alias_volume_names,
    _parse_blastdbcmd_info,
)
from analyses._blast.local_engine import LocalBlastEngine, _should_run_short_query_rescue
from analyses._blast.evidence import (
    alignment_text,
    build_evidence_profile,
    calculate_evidence_score,
    group_hits_by_taxonomy,
)
from analyses._blast.reporting import (
    build_text_report,
    compare_blast_runs,
    export_csv,
    export_json,
    export_text,
    hit_to_dict,
    hsp_to_dict,
    run_to_dict,
    run_to_json,
)
from analyses._blast.workflow import run_blast_analysis
from analyses._blast.selftest import (
    _assert,
    run_self_test,
    SAMPLE_XML,
    _test_alias_inside_first_volume_directory,
    _test_alignment,
    _test_blastdbcmd_info_parser,
    _test_comparison,
    _test_default_database_root_resolution,
    _test_evidence_profile,
    _test_executable_manager,
    _test_export_schema,
    _test_fasta_query_loader,
    _test_hit_ranking,
    _test_local_blast_environment,
    _test_local_blast_environment_preserves_existing,
    _test_local_database_command_reference,
    _test_local_database_manifest,
    _test_local_database_structure,
    _test_local_parameter_compatibility,
    _test_ncbi_client_configuration,
    _test_ncbi_request_parameters,
    _test_nested_local_database_discovery,
    _test_nested_multivolume_database_discovery,
    _test_nested_runtime_alias,
    _test_orf_adapter,
    _test_program_database_compatibility,
    _test_query_builder,
    _test_sequence_validation,
    _test_serialization,
    _test_short_query_rescue_policy,
    _test_short_query_task_configuration,
    _test_synapse_user_path_resolution,
    _test_taxonomy_definition_fallback,
    _test_taxonomy_grouping_with_fallback,
    _test_type_inference,
    _test_xml_parser,
)


# ============================================================================
# MODULE ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    from ui import style  # only when run as a script
    if not run_self_test(verbose=True):
        raise SystemExit(1)

    print()
    style.success("Synapse BLAST module is ready.")
    style.note("No network search was executed during self-test.")
    style.note("Use blast_menu(...) from main.py for interactive BLAST Search.")
