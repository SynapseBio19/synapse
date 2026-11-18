"""The built-in self-test of the BLAST module (menu entry and test suite)."""

from __future__ import annotations

import json
import os
import tempfile
from core.textformat import header_text, minor_rule
from pathlib import Path
from typing import List

from analyses._blast.paths import ANALYSES_DIR, _resolve_synapse_user_path
from analyses._blast.config import (
    BLAST_EXECUTABLES,
    DEFAULT_LOCAL_DB_ROOT,
    NCBI_MIN_POLL_INTERVAL_SECONDS,
    NCBI_MIN_REQUEST_INTERVAL_SECONDS,
    PROGRAMS,
)
from analyses._blast.models import (
    BlastAnalysisError,
    BlastHit,
    BlastParameters,
    rank_hits,
    SelfTestResult,
)
from analyses._blast.query import (
    build_orf_query,
    build_query,
    infer_sequence_type,
    load_fasta_queries,
    normalize_sequence,
)
from analyses._blast.parsing import _infer_scientific_name_from_definition, parse_blast_xml
from analyses._blast.ncbi import NCBIBlastClient
from analyses._blast.executables import BlastExecutableManager
from analyses._blast.runtime import (
    _build_local_blast_environment,
    _build_local_database_command_reference,
    _create_runtime_blast_alias,
    _needs_runtime_blast_alias,
)
from analyses._blast.databases import LocalBlastDatabaseManager, _parse_blastdbcmd_info
from analyses._blast.local_engine import _should_run_short_query_rescue
from analyses._blast.evidence import alignment_text, build_evidence_profile, group_hits_by_taxonomy
from analyses._blast.reporting import compare_blast_runs, run_to_dict, run_to_json


# ============================================================================
# DISPLAY
# ============================================================================


# ============================================================================
# INTERACTIVE MENU
# ============================================================================


# ============================================================================
# SELF-TEST SUITE
# ============================================================================

SAMPLE_XML = """<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.17.0+</BlastOutput_version>
  <BlastOutput_db>nt</BlastOutput_db>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>ref|TEST001|</Hit_id>
          <Hit_def>Example reference sequence</Hit_def>
          <Hit_accession>TEST001</Hit_accession>
          <Hit_len>120</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>100</Hsp_bit-score>
              <Hsp_score>110</Hsp_score>
              <Hsp_evalue>1e-25</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>20</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>29</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>19</Hsp_identity>
              <Hsp_positive>19</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>20</Hsp_align-len>
              <Hsp_qseq>ATGCGTACGTACGTACGTAC</Hsp_qseq>
              <Hsp_hseq>ATGCGTACGTACGTACGTAT</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||| </Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_statistics>
    <Statistics>
      <Statistics_db-num>1000</Statistics_db-num>
      <Statistics_db-len>1000000</Statistics_db-len>
    </Statistics>
  </BlastOutput_statistics>
</BlastOutput>
"""


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _test_sequence_validation() -> None:
    _assert(normalize_sequence(" a c g t n ") == "ACGTN", "DNA normalization failed")
    try:
        normalize_sequence("ACGTZ", "DNA")
    except BlastAnalysisError:
        pass
    else:
        raise AssertionError("Invalid DNA was accepted")


def _test_type_inference() -> None:
    _assert(infer_sequence_type("ACGTNN") == "DNA", "DNA inference failed")
    _assert(infer_sequence_type("MKTLLV") == "PROTEIN", "Protein inference failed")


def _test_program_database_compatibility() -> None:
    BlastParameters(program="blastn", database="nt").normalized()
    BlastParameters(program="blastp", database="nr").normalized()

    try:
        BlastParameters(program="blastn", database="nr").normalized()
    except BlastAnalysisError:
        pass
    else:
        raise AssertionError("Incompatible program/database accepted")


def _test_query_builder() -> None:
    query = build_query(
        "ACGT",
        name="Test",
        sequence_type="DNA",
        source="Unit Test",
        region_start=5,
    )
    _assert(query.length == 4, "Query length incorrect")
    _assert(query.region_start == 5, "Query start incorrect")
    _assert(query.region_end == 8, "Query end incorrect")


def _test_fasta_query_loader() -> None:
    with tempfile.TemporaryDirectory(prefix="synapse_fasta_query_test_") as temp_dir:
        fasta = Path(temp_dir) / "query.fasta"
        fasta.write_text(
            ">Query_A\r\nACGTACGT\r\n>Query_B\r\nGGCCGGCC\r\n",
            encoding="utf-8",
        )

        queries = load_fasta_queries(fasta)
        _assert(len(queries) == 2, "FASTA loader did not recover both records")
        _assert(queries[0].name == "Query_A", "First FASTA header was parsed incorrectly")
        _assert(queries[0].sequence == "ACGTACGT", "First FASTA sequence was parsed incorrectly")
        _assert(queries[1].name == "Query_B", "Second FASTA header was parsed incorrectly")
        _assert(queries[1].source == "FASTA File", "FASTA query source was not recorded")


def _test_xml_parser() -> None:
    query = build_query(
        "ATGCGTACGTACGTACGTAC",
        name="Test_Query",
        sequence_type="DNA",
    )
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="TESTRID",
        parameters=params,
        submitted_at="2026-01-01T00:00:00+00:00",
        completed_at="2026-01-01T00:01:00+00:00",
    )

    _assert(run.hit_count == 1, "XML parser did not recover one hit")
    _assert(run.hits[0].accession == "TEST001", "Accession parsing failed")
    _assert(run.hits[0].best_identity_percent == 95.0, "Identity calculation failed")
    _assert(run.hits[0].best_coverage_percent == 100.0, "Coverage calculation failed")


def _test_hit_ranking() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )

    ranked = rank_hits(run.hits)
    _assert(ranked[0].accession == "TEST001", "Ranking failed")


def _test_evidence_profile() -> None:
    # SAMPLE_XML aligns query positions 1-20; a 20-nt query keeps the
    # fixture consistent with the query-coverage definition.
    query = build_query("ATGCGTACGTACGTACGTAC", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )

    profile = build_evidence_profile(run)[0]
    _assert(profile.identity_percent == 95.0, "Evidence identity failed")
    _assert(profile.query_coverage_percent == 100.0, "Evidence coverage failed")
    _assert(0 <= profile.evidence_score <= 100, "Evidence score out of range")
    _assert("experimental" not in profile.evidence_band.lower(), "Invalid evidence label")


def _test_alignment() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )
    text = alignment_text(run.hits[0])
    _assert("Query" in text and "Subject" in text, "Alignment rendering failed")


def _test_serialization() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )
    payload = json.loads(run_to_json(run))
    _assert(payload["metadata"]["rid"] == "RID", "JSON metadata failed")
    _assert(payload["hits"][0]["accession"] == "TEST001", "JSON hit failed")


def _test_comparison() -> None:
    query = build_query("ACGT", name="A", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run_a = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="A",
        parameters=params,
        submitted_at="now",
    )

    query_b = build_query("ACGT", name="B", sequence_type="DNA")
    run_b = parse_blast_xml(
        SAMPLE_XML,
        query=query_b,
        rid="B",
        parameters=params,
        submitted_at="now",
    )

    comparison = compare_blast_runs(run_a, run_b)
    _assert(comparison.common_accessions == ["TEST001"], "Comparison failed")


def _test_orf_adapter() -> None:
    orf = {
        "id": "ORF-001",
        "sequence": "ATGCGT",
        "start": 101,
    }
    query = build_orf_query(orf)
    _assert(query.name == "ORF-001", "ORF name extraction failed")
    _assert(query.region_start == 101, "ORF coordinate extraction failed")
    _assert(query.region_end == 106, "ORF end coordinate failed")
    _assert(query.source == "Selected ORF", "ORF source metadata failed")


def _test_export_schema() -> None:
    query = build_query("ACGT", name="Q", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    run = parse_blast_xml(
        SAMPLE_XML,
        query=query,
        rid="RID",
        parameters=params,
        submitted_at="now",
    )
    payload = run_to_dict(run)
    _assert("evidence_profile" in payload, "Evidence profile missing")
    _assert("database_statistics" in payload, "Database statistics missing")
    _assert("interpretation_notes" in payload, "Interpretation notes missing")


def _test_ncbi_client_configuration() -> None:
    client = NCBIBlastClient(tool="Synapse Test Tool", email="test@example.org")
    _assert(client.tool == "Synapse_Test_Tool", "NCBI tool identifier was not normalized")
    _assert(client.email == "test@example.org", "NCBI email configuration was lost")
    _assert(
        client.min_request_interval >= NCBI_MIN_REQUEST_INTERVAL_SECONDS,
        "NCBI request interval is below the documented minimum",
    )
    _assert(
        client.min_poll_interval >= NCBI_MIN_POLL_INTERVAL_SECONDS,
        "NCBI polling interval is below the documented minimum",
    )


def _test_ncbi_request_parameters() -> None:
    client = NCBIBlastClient(tool="Synapse", email="dev@example.org")
    query = build_query("ACGT", name="NCBI_Query", sequence_type="DNA")
    params = BlastParameters(program="blastn", database="nt")
    normalized = params.normalized()
    payload = {
        "CMD": "Put",
        "PROGRAM": normalized.program,
        "DATABASE": normalized.database,
        "QUERY": f">{query.name}\r\n{query.sequence}",
        "EXPECT": normalized.evalue,
        "HITLIST_SIZE": normalized.max_target_seqs,
        "TOOL": client.tool,
        "EMAIL": client.email,
    }
    _assert(payload["CMD"] == "Put", "NCBI command parameter missing")
    _assert(payload["TOOL"] == "Synapse", "NCBI tool parameter missing")
    _assert(payload["EMAIL"] == "dev@example.org", "NCBI email parameter missing")
    _assert(payload["DATABASE"] == "nt", "NCBI database parameter incorrect")


def _test_local_database_command_reference() -> None:
    """Ensure local BLAST uses a basename/cwd pair for robust Windows paths."""
    with tempfile.TemporaryDirectory(prefix="Synapse Project With Spaces ") as temp_dir:
        database = Path(temp_dir) / "nested folder" / "synapse_test"
        database.parent.mkdir(parents=True, exist_ok=True)
        reference, working_directory = _build_local_database_command_reference(database)
        _assert(reference == "synapse_test", "Database reference should be the prefix name")
        _assert(
            working_directory == database.parent.resolve(),
            "Database working directory was resolved incorrectly",
        )


def _test_local_blast_environment() -> None:
    """Ensure BLASTDB includes the database directory without losing existing paths."""
    with tempfile.TemporaryDirectory(prefix="synapse_env_test_") as temp_dir:
        database = Path(temp_dir) / "db" / "testdb"
        env = _build_local_blast_environment(database)
        _assert(
            env["BLASTDB"] == str(database.parent.resolve()),
            "BLASTDB was not initialized with the database directory",
        )


def _test_local_blast_environment_preserves_existing() -> None:
    """Ensure an existing BLASTDB search path is preserved."""
    with tempfile.TemporaryDirectory(prefix="synapse_env_test_") as temp_dir:
        database = Path(temp_dir) / "db" / "testdb"
        old_value = os.environ.get("BLASTDB")
        try:
            os.environ["BLASTDB"] = "existing_database_root"
            env = _build_local_blast_environment(database)
        finally:
            if old_value is None:
                os.environ.pop("BLASTDB", None)
            else:
                os.environ["BLASTDB"] = old_value

        parts = env["BLASTDB"].split(os.pathsep)
        _assert(str(database.parent.resolve()) in parts, "Database root was not prepended")
        _assert("existing_database_root" in parts, "Existing BLASTDB path was lost")


def _test_default_database_root_resolution() -> None:
    """Ensure the default root is a project-level data/blast_db path."""
    expected_name = Path("data") / "blast_db"
    resolved = DEFAULT_LOCAL_DB_ROOT.resolve()
    _assert(
        resolved.name == expected_name.name,
        "Default BLAST root must end in data/blast_db",
    )
    _assert(
        resolved.parent.name == expected_name.parent.name,
        "Default BLAST root parent must be data",
    )


def _test_local_database_structure() -> None:
    with tempfile.TemporaryDirectory(prefix="synapse_db_test_") as temp_dir:
        root = Path(temp_dir)
        (root / "testdb.nhr").write_bytes(b"nhr")
        (root / "testdb.nin").write_bytes(b"nin")
        (root / "testdb.nsq").write_bytes(b"nsq")

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("testdb", "nucleotide")

        _assert(database.installed, "Local database was not detected")
        _assert(database.structurally_valid, "Local database structure was not valid")
        _assert(database.file_count == 3, "Unexpected local database file count")


def _test_nested_runtime_alias() -> None:
    """Verify runtime aliases point to sibling nested volume prefixes."""
    with tempfile.TemporaryDirectory(prefix="synapse_runtime_alias_test_") as temp_dir:
        root = Path(temp_dir)
        container = root / "refseq_rna"
        first_volume = container / "refseq_rna.00"
        second_volume = container / "refseq_rna.01"
        first_volume.mkdir(parents=True)
        second_volume.mkdir(parents=True)

        for volume in (first_volume, second_volume):
            prefix = volume / volume.name
            for extension in (".nhr", ".nin", ".nsq"):
                (Path(str(prefix) + extension)).write_bytes(b"test")

        alias = first_volume / "refseq_rna.nal"
        alias.write_text(
            "DBLIST refseq_rna.00 refseq_rna.01\r\n"
            "TITLE Synapse nested runtime alias test\r\n",
            encoding="utf-8",
        )

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("refseq_rna", "nucleotide")
        _assert(database.structurally_valid, "Nested runtime-alias test database was not valid")
        _assert(_needs_runtime_blast_alias(database), "Nested alias did not require runtime alias")

        with tempfile.TemporaryDirectory(prefix="synapse_runtime_alias_out_", dir=str(container)) as runtime_dir:
            reference, working_directory = _create_runtime_blast_alias(
                database, runtime_dir
            )
            runtime_alias = Path(working_directory) / (reference + ".nal")
            payload = runtime_alias.read_text(encoding="utf-8")
            _assert(runtime_alias.suffix == ".nal", "Runtime alias extension was incorrect")
            _assert(reference == "refseq_rna", "Runtime BLAST reference must use the canonical database root")
            _assert("DBLIST refseq_rna.00 refseq_rna.01" in payload, "Runtime alias did not point to volume 00")
            _assert("DBLIST refseq_rna.00 refseq_rna.01" in payload, "Runtime alias did not point to volume 01")
            _assert("DBLIST refseq_rna.00 refseq_rna.01" in payload, "Runtime DBLIST was malformed")
            _assert((Path(working_directory) / "refseq_rna.00.nsq").exists(), "Runtime flat view is missing volume 00")
            _assert((Path(working_directory) / "refseq_rna.01.nsq").exists(), "Runtime flat view is missing volume 01")
            _assert((Path(working_directory) / "refseq_rna.00.nsq").stat().st_ino == (first_volume / "refseq_rna.00.nsq").stat().st_ino, "Runtime volume 00 was copied instead of hard-linked")


def _test_blastdbcmd_info_parser() -> None:
    sample = """Database: Synapse Local Test Database
        3 sequences; 153 total bases

BLASTDB Version: 5

Volumes:
        C:\\Synapse\\synapse_test
"""
    info = _parse_blastdbcmd_info(sample)
    _assert(info["title"] == "Synapse Local Test Database", "Database title parsing failed")
    _assert(info["sequence_count"] == 3, "Sequence count parsing failed")
    _assert(info["total_length"] == 153, "Total length parsing failed")
    _assert(info["blastdb_version"] == "5", "BLASTDB version parsing failed")
    _assert(len(info["volumes"]) == 1, "Volume parsing failed")


def _test_nested_local_database_discovery() -> None:
    """Ensure Synapse detects BLAST DB v5-style files in subdirectories."""
    with tempfile.TemporaryDirectory(prefix="synapse_nested_db_test_") as temp_dir:
        root = Path(temp_dir)
        nested = root / "local_test"
        nested.mkdir(parents=True)

        # BLAST DB v5 commonly exposes a .ndb core file. Additional files make
        # the fixture resemble a real nucleotide database layout.
        (nested / "synapse_test.ndb").write_bytes(b"ndb")
        (nested / "synapse_test.nhr").write_bytes(b"nhr")
        (nested / "synapse_test.nin").write_bytes(b"nin")
        (nested / "synapse_test.nsq").write_bytes(b"nsq")

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("synapse_test")

        _assert(database.installed, "Nested BLAST database was not detected")
        _assert(database.structurally_valid, "BLAST DB v5 structure was not recognized")
        _assert(database.name == "local_test/synapse_test", "Nested database name was not resolved")
        _assert(
            Path(database.path) == (nested / "synapse_test").resolve(),
            "Nested database prefix was resolved incorrectly",
        )

        listed = manager.list_databases()
        _assert(
            any(item.name == "local_test/synapse_test" for item in listed),
            "Nested database was not returned by list_databases",
        )


def _test_nested_multivolume_database_discovery() -> None:
    """Verify the extracted NCBI-style folder-per-volume layout."""
    with tempfile.TemporaryDirectory(prefix="synapse_multivolume_db_test_") as temp_dir:
        root = Path(temp_dir)
        container = root / "refseq_rna"
        container.mkdir(parents=True)
        alias = container / "refseq_rna.nal"
        alias.write_text(
            "DBLIST refseq_rna.00 refseq_rna.01\r\n"
            "TITLE Synapse RefSeq RNA test\r\n",
            encoding="utf-8",
        )

        for volume in ("00", "01"):
            folder = container / f"refseq_rna.{volume}"
            folder.mkdir()
            for extension in (".nhr", ".nin", ".nsq"):
                (folder / f"refseq_rna.{volume}{extension}").write_bytes(b"test")

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("refseq_rna", "nucleotide")

        _assert(database.installed, "Nested multi-volume database was not detected")
        _assert(database.structurally_valid, "Complete multi-volume database was not recognized")
        _assert(database.name == "refseq_rna", "Logical database name changed unexpectedly")
        _assert(database.volume_count == 2, "Expected two BLAST volumes")
        _assert(database.expected_volume_count == 2, "Alias DBLIST volume count was not parsed")
        _assert(not database.missing_volumes, "Complete database reported missing volumes")
        _assert(not database.incomplete_volumes, "Complete volumes were marked incomplete")
        _assert(
            all(Path(path).is_dir() for path in database.volume_paths),
            "Volume directory paths were not discovered correctly",
        )

        # Alias numbering can be unpadded while the extracted filesystem
        # uses the canonical two-digit NCBI volume directories. The resolver
        # must match .0 to .00 (and .1 to .01) by numeric index.
        alias.write_text(
            "DBLIST refseq_rna.0 refseq_rna.1\r\n"
            "TITLE Synapse RefSeq RNA test\r\n",
            encoding="utf-8",
        )
        unpadded = manager.scan_database("refseq_rna", "nucleotide")
        _assert(unpadded.structurally_valid, "Unpadded .0/.1 alias did not match .00/.01 volumes")
        _assert(unpadded.volume_count == 2, "Unpadded alias volume count was incorrect")
        _assert(not unpadded.missing_volumes, "Unpadded alias incorrectly reported missing volumes")
        _assert(not unpadded.incomplete_volumes, "Unpadded alias marked complete volumes incomplete")

        # Restore the canonical padded alias for the remaining damage/missing
        # volume assertions.
        alias.write_text(
            "DBLIST refseq_rna.00 refseq_rna.01\r\n"
            "TITLE Synapse RefSeq RNA test\r\n",
            encoding="utf-8",
        )

        # A misleading prefix inside a volume directory must never replace
        # the parent database alias.
        decoy = container / "refseq_rna.00" / "refseq_rna.nhr"
        decoy.write_bytes(b"decoy")
        resolved = manager.scan_database("refseq_rna", "nucleotide")
        _assert(
            Path(resolved.path) == (container / "refseq_rna").resolve(),
            "A volume-internal prefix incorrectly became the database prefix",
        )

        # Remove one required file to ensure the manager refuses to call the
        # database structurally complete when an extracted volume is damaged.
        (container / "refseq_rna.00" / "refseq_rna.00.nsq").unlink()
        damaged = manager.scan_database("refseq_rna", "nucleotide")
        _assert(not damaged.structurally_valid, "Damaged volume was accepted as valid")
        _assert(
            "refseq_rna.00" in damaged.incomplete_volumes,
            "Damaged volume was not reported",
        )
        _assert(
            damaged.missing_volume_files.get("refseq_rna.00") == [".nsq"],
            "Missing volume file was not reported precisely",
        )

        # Remove the second volume completely to verify missing-volume reporting.
        for path in (container / "refseq_rna.01").iterdir():
            path.unlink()
        (container / "refseq_rna.01").rmdir()
        missing = manager.scan_database("refseq_rna", "nucleotide")
        _assert(
            "refseq_rna.01" in missing.missing_volumes,
            "Missing BLAST volume was not reported",
        )


def _test_alias_inside_first_volume_directory() -> None:
    """Verify an alias stored inside .00 can still expose all sibling volumes."""
    with tempfile.TemporaryDirectory(prefix="synapse_alias_in_first_volume_") as temp_dir:
        root = Path(temp_dir)
        container = root / "refseq_rna"
        container.mkdir(parents=True)
        first = container / "refseq_rna.00"
        first.mkdir()
        (first / "refseq_rna.nal").write_text(
            "DBLIST refseq_rna.00 refseq_rna.01\r\nTITLE Test\r\n",
            encoding="utf-8",
        )
        # Shared alias metadata may also be stored beside the alias.
        (first / "refseq_rna.ndb").write_bytes(b"ndb")
        for volume in ("00", "01"):
            folder = container / f"refseq_rna.{volume}"
            if not folder.exists():
                folder.mkdir()
            for extension in (".nhr", ".nin", ".nsq"):
                (folder / f"refseq_rna.{volume}{extension}").write_bytes(b"test")

        manager = LocalBlastDatabaseManager(root=root)
        database = manager.scan_database("refseq_rna", "nucleotide")
        _assert(database.structurally_valid, "Alias inside .00 was not recognized as complete")
        _assert(database.volume_count == 2, "Sibling volumes were not discovered")
        _assert(database.expected_volume_count == 2, "Alias DBLIST count was not preserved")
        _assert(not database.missing_volumes, "Sibling volume was incorrectly reported missing")
        _assert(not database.incomplete_volumes, "Complete sibling volumes were marked incomplete")


def _test_local_database_manifest() -> None:
    with tempfile.TemporaryDirectory(prefix="synapse_manifest_test_") as temp_dir:
        root = Path(temp_dir)
        manager = LocalBlastDatabaseManager(root=root)
        manager.ensure_root()

        manifest = {
            "custom_db": {
                "database_type": "protein",
                "source": "Custom FASTA",
            }
        }
        manager._save_manifest(manifest)
        loaded = manager._load_manifest()

        _assert(loaded["custom_db"]["database_type"] == "protein", "Manifest roundtrip failed")


def _test_local_parameter_compatibility() -> None:
    manager = LocalBlastDatabaseManager()
    with tempfile.TemporaryDirectory(prefix="synapse_db_test_") as temp_dir:
        root = Path(temp_dir)
        (root / "protein_db.phr").write_bytes(b"phr")
        (root / "protein_db.pin").write_bytes(b"pin")
        (root / "protein_db.psq").write_bytes(b"psq")
        database = manager.__class__(root=root).scan_database("protein_db", "protein")

        query = build_query("MKTLLV", name="Protein_Query", sequence_type="PROTEIN")
        params = BlastParameters(program="blastp", database="protein_db")

        _assert(
            PROGRAMS[params.program]["database_type"] == database.database_type,
            "Program/database type compatibility failed",
        )
        _assert(query.sequence_type == "PROTEIN", "Protein query type failed")


def _test_synapse_user_path_resolution() -> None:
    application_root = ANALYSES_DIR.parent
    relative_input = Path("relative_test") / "result.json"
    relative = _resolve_synapse_user_path(relative_input)
    system_root = Path(os.environ.get("SystemRoot", r"C:\Windows")).resolve()
    expected_base = application_root if (
        os.name == "nt" and Path.cwd().resolve() == system_root / "System32"
    ) else Path.cwd().resolve()
    _assert(
        relative == (expected_base / relative_input).resolve(),
        "Relative user path was resolved incorrectly",
    )

    absolute = _resolve_synapse_user_path(application_root / "data")
    _assert(absolute == (application_root / "data").resolve(),
            "Absolute path was modified incorrectly")


def _test_executable_manager() -> None:
    manager = BlastExecutableManager()
    available = manager.available()

    _assert(set(BLAST_EXECUTABLES) == set(available), "Executable catalog mismatch")
    _assert(manager.resolve_perl() is None or isinstance(manager.resolve_perl(), str),
            "Perl discovery returned an invalid value")


def _test_short_query_rescue_policy() -> None:
    query = build_query(
        "ATGAAACGTGCTGACTTTGGTGACGCTTACCGTAAAGGTC" * 2,
        name="Short_DNA_Query",
        sequence_type="DNA",
    )
    default_parameters = BlastParameters(
        program="blastn",
        database="nt",
        evalue=10.0,
    )
    _assert(
        _should_run_short_query_rescue(query, default_parameters),
        "Default short-query rescue policy should be eligible",
    )

    explicit_word_size = BlastParameters(
        program="blastn",
        database="nt",
        word_size=7,
    )
    _assert(
        not _should_run_short_query_rescue(query, explicit_word_size),
        "Explicit word size must disable automatic rescue",
    )


def _test_short_query_task_configuration() -> None:
    parameters = BlastParameters(
        program="blastn",
        database="nt",
        task="blastn-short",
        filter_query="no",
    ).normalized()
    _assert(parameters.task == "blastn-short", "Short-query task was not preserved")
    _assert(parameters.filter_query == "no", "Short-query filtering mode was not preserved")


def _test_taxonomy_definition_fallback() -> None:
    name = _infer_scientific_name_from_definition(
        "PREDICTED: Nicrophorus vespilloides protein UBASH3A homolog"
    )
    _assert(name == "Nicrophorus vespilloides", "Failed to infer binomial organism name")

    bracketed = _infer_scientific_name_from_definition(
        "Reference transcript [Homo sapiens]"
    )
    _assert(bracketed == "Homo sapiens", "Failed to parse bracketed organism name")

    _assert(
        _infer_scientific_name_from_definition("PREDICTED: hypothetical protein") is None,
        "Should not infer a non-taxonomic phrase",
    )


def _test_taxonomy_grouping_with_fallback() -> None:
    hit_a = BlastHit(
        hit_id="A",
        accession="A",
        definition="PREDICTED: Nicrophorus vespilloides protein UBASH3A homolog",
        length=100,
    )
    hit_b = BlastHit(
        hit_id="B",
        accession="B",
        definition="PREDICTED: Nicrophorus vespilloides transcript variant",
        length=100,
    )
    hit_c = BlastHit(
        hit_id="C",
        accession="C",
        definition="Example reference sequence",
        length=100,
    )
    for hit in (hit_a, hit_b, hit_c):
        inferred = _infer_scientific_name_from_definition(hit.definition)
        if inferred:
            hit.scientific_names.append(inferred)

    groups = group_hits_by_taxonomy((hit_a, hit_b, hit_c))
    _assert(groups.get("Nicrophorus vespilloides") == 2, "Taxonomy grouping failed")
    _assert(groups.get("Unknown / not provided") == 1, "Unknown taxonomy count failed")


def run_self_test(verbose: bool = True) -> bool:
    tests = [
        ("sequence_validation", _test_sequence_validation),
        ("type_inference", _test_type_inference),
        ("program_database_compatibility", _test_program_database_compatibility),
        ("query_builder", _test_query_builder),
        ("fasta_query_loader", _test_fasta_query_loader),
        ("xml_parser", _test_xml_parser),
        ("hit_ranking", _test_hit_ranking),
        ("evidence_profile", _test_evidence_profile),
        ("alignment_renderer", _test_alignment),
        ("json_serialization", _test_serialization),
        ("run_comparison", _test_comparison),
        ("orf_adapter", _test_orf_adapter),
        ("export_schema", _test_export_schema),
        ("ncbi_client_configuration", _test_ncbi_client_configuration),
        ("ncbi_request_parameters", _test_ncbi_request_parameters),
        ("local_database_command_reference", _test_local_database_command_reference),
        ("local_blast_environment", _test_local_blast_environment),
        ("local_blast_environment_preserves_existing", _test_local_blast_environment_preserves_existing),
        ("default_database_root_resolution", _test_default_database_root_resolution),
        ("local_database_structure", _test_local_database_structure),
        ("nested_local_database_discovery", _test_nested_local_database_discovery),
        ("nested_multivolume_database_discovery", _test_nested_multivolume_database_discovery),
        ("nested_runtime_alias", _test_nested_runtime_alias),
        ("alias_inside_first_volume_directory", _test_alias_inside_first_volume_directory),
        ("blastdbcmd_info_parser", _test_blastdbcmd_info_parser),
        ("local_database_manifest", _test_local_database_manifest),
        ("local_parameter_compatibility", _test_local_parameter_compatibility),
        ("synapse_user_path_resolution", _test_synapse_user_path_resolution),
        ("executable_manager", _test_executable_manager),
        ("short_query_rescue_policy", _test_short_query_rescue_policy),
        ("short_query_task_configuration", _test_short_query_task_configuration),
        ("taxonomy_definition_fallback", _test_taxonomy_definition_fallback),
        ("taxonomy_grouping_with_fallback", _test_taxonomy_grouping_with_fallback),
    ]

    results: List[SelfTestResult] = []

    for name, test in tests:
        try:
            test()
            results.append(SelfTestResult(name, True))
        except Exception as exc:
            results.append(SelfTestResult(name, False, str(exc)))

    if verbose:
        print(header_text('BLAST self-test'))
        for result in results:
            marker = "PASS" if result.passed else "FAIL"
            suffix = f" - {result.detail}" if result.detail else ""
            print(f"[{marker}] {result.name}{suffix}")
        passed = sum(result.passed for result in results)
        print(minor_rule())
        print(f"{passed}/{len(results)} tests passed.")

    return all(result.passed for result in results)
