"""Split plan for analyses/blast.py (Phase 7.1, package B)."""
import ast

SOURCE = "analyses/blast.py"
PACKAGE = "analyses._blast"
PACKAGE_DIR = "analyses/_blast"
PACKAGE_DOC = """Implementation of :mod:`analyses.blast` (split out in Phase 7.1).

Import from ``analyses.blast``; it re-exports every name defined here and
is the only public interface. The modules, in dependency order:

- ``paths``: the ``analyses`` directory and user-path resolution
- ``config``: endpoints, rate limits, programs, databases, local catalog
- ``models``: the dataclasses, ``BlastAnalysisError`` and hit ranking
- ``query``: query validation, FASTA loading, the ORF adapter
- ``parsing``: BLAST XML parsing
- ``ncbi``: the NCBI BLAST URL API client
- ``executables``: locating the BLAST+ programs
- ``runtime``: runtime alias views and the BLAST environment
- ``databases``: the local database manager
- ``local_engine``: local BLAST+ execution
- ``evidence``: evidence score, taxonomy grouping, alignment text
- ``reporting``: serialization, exports, text report, run comparison
- ``workflow``: ``run_blast_analysis``
- ``selftest``: the built-in self-test
"""

_tree = ast.parse(open(SOURCE, encoding="utf-8").read())
_TESTS = [n.name for n in _tree.body if isinstance(n, ast.FunctionDef) and n.name.startswith("_test_")]

MODULES = {
    "paths": ("Locations Synapse resolves relative to its installation.", [
        "_resolve_synapse_user_path",
    ]),
    "config": ("Fixed configuration of the BLAST module: NCBI endpoint and rate limits,\n"
               "programs, remote databases and the local database catalog.", [
        "NCBI_BLAST_URL", "DEFAULT_TOOL", "DEFAULT_EMAIL", "API_TIMEOUT_SECONDS",
        "NCBI_MIN_REQUEST_INTERVAL_SECONDS", "NCBI_MIN_POLL_INTERVAL_SECONDS",
        "DEFAULT_POLL_INTERVAL_SECONDS", "DEFAULT_MAX_WAIT_SECONDS", "DEFAULT_MAX_HITS",
        "DEFAULT_EVIDENCE_EVALUE_THRESHOLD", "LOCAL_DB_ROOT_ENV", "LOCAL_DB_MANIFEST_NAME",
        "_resolve_default_local_db_root", "DEFAULT_LOCAL_DB_ROOT", "LOCAL_DATABASE_CATALOG",
        "BLAST_EXECUTABLES", "BLAST_UPDATE_SCRIPT", "DNA_BASES", "PROTEIN_AMINO_ACIDS",
        "AMBIGUOUS_DNA", "PROGRAMS", "QUERY_TYPES", "DATABASES", "BLAST_TAG_ALIASES",
    ]),
    "models": ("Data model of BLAST queries, parameters, runs, hits and databases.\n\n"
               "``rank_hits`` lives here because ``BlastRun.ranked_hits`` uses it.", [
        "BlastParameters", "BlastHSP", "BlastHit", "BlastQuery", "BlastRunMetadata", "BlastRun",
        "HitEvidence", "BlastComparison", "SelfTestResult", "BlastAnalysisError", "_utc_now",
        "_evalue_text", "rank_hits", "LocalBlastDatabase",
    ]),
    "query": ("Query construction and validation: sequence normalization, type inference,\n"
              "program compatibility, FASTA queries and the ORF adapter.", [
        "programs_for_query_type", "check_program_query_type", "normalize_sequence",
        "infer_sequence_type", "build_query", "load_fasta_queries", "extract_orf_context",
        "build_orf_query",
    ]),
    "parsing": "BLAST XML (outfmt 5) parsing into the data model.",
    "ncbi": ("Client of the NCBI BLAST URL API: submission, polling with the NCBI rate\n"
             "limits, result retrieval and resuming a saved RID.", [
        "logger", "STATUS_RETRIES", "_is_transient_ncbi_error", "NCBIBlastClient",
    ]),
    "executables": ("Locating the NCBI BLAST+ programs and the database update script.", [
        "BlastExecutableManager",
    ]),
    "runtime": ("Runtime views of local databases: command references, temporary alias\n"
                "files for nested volumes, and the BLAST process environment.", [
        "_build_local_database_command_reference", "_volume_prefix_from_stored_path",
        "_create_runtime_blast_alias", "_needs_runtime_blast_alias", "_build_local_blast_environment",
    ]),
    "databases": ("The local BLAST database manager: discovery, verification, manifest,\n"
                  "NCBI download, custom databases and removal.", [
        "_parse_blastdbcmd_info", "_extract_blast_volume_index", "_parse_blast_alias_volume_names",
        "LocalBlastDatabaseManager",
    ]),
    "local_engine": "Running a search with local BLAST+ programs.",
    "evidence": "Evidence score, evidence profile, taxonomy grouping and alignment text.",
    "reporting": "Serialization, CSV/JSON/text exports, the text report and run comparison.",
    "workflow": "``run_blast_analysis``: one search, remote or local, from a sequence or an ORF.",
    "selftest": "The built-in self-test of the BLAST module (menu entry and test suite).",
}

MODULES["parsing"] = (MODULES["parsing"], [
    "_local_tag", "_child_text", "_descendants", "_to_int", "_to_float", "_parse_hsp",
    "_infer_scientific_name_from_definition", "_parse_hit", "_parse_parameters",
    "_parse_database_statistics", "describe_run_from_xml", "parse_blast_xml",
])
MODULES["local_engine"] = (MODULES["local_engine"], ["_should_run_short_query_rescue", "LocalBlastEngine"])
MODULES["evidence"] = (MODULES["evidence"], [
    "calculate_evidence_score", "build_evidence_profile", "group_hits_by_taxonomy", "alignment_text",
])
MODULES["reporting"] = (MODULES["reporting"], [
    "hsp_to_dict", "hit_to_dict", "run_to_dict", "run_to_json", "export_json", "export_csv",
    "build_text_report", "export_text", "compare_blast_runs",
])
MODULES["workflow"] = (MODULES["workflow"], ["run_blast_analysis"])
MODULES["selftest"] = (MODULES["selftest"], ["SAMPLE_XML", "_assert", *_TESTS, "run_self_test"])

# Before Phase 7.1 these lookups ran in analyses/blast.py. The code now lives
# one directory deeper, so it uses the analyses directory explicitly.
REWRITES = [
    ("Path(__file__).resolve().parent.parent", "ANALYSES_DIR.parent"),
    ("Path(__file__).resolve().parent", "ANALYSES_DIR"),
]
EXTRA = {
    "ANALYSES_DIR": ("paths", (
        "# The ``analyses`` directory. Before Phase 7.1 the BLAST code lived in\n"
        "# analyses/blast.py and used ``Path(__file__).resolve().parent`` for it.\n"
        "ANALYSES_DIR = Path(__file__).resolve().parent.parent\n")),
}
EXPORT_EXTRA = set()

FACADE_NOTE = (
    "# Phase 7.1: the implementation moved to the internal package\n"
    "# ``analyses._blast`` (one module per responsibility, see its __init__).\n"
    "# This module remains the public interface and re-exports every name it\n"
    "# defined before, private helpers included.\n"
)
