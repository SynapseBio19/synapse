"""``run_blast_analysis``: one search, remote or local, from a sequence or an ORF."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from analyses._blast.config import DEFAULT_EMAIL, DEFAULT_MAX_HITS, DEFAULT_TOOL, PROGRAMS
from analyses._blast.models import BlastAnalysisError, BlastParameters, BlastRun
from analyses._blast.query import build_orf_query, build_query
from analyses._blast.ncbi import NCBIBlastClient
from analyses._blast.databases import LocalBlastDatabaseManager
from analyses._blast.local_engine import LocalBlastEngine


def run_blast_analysis(
    sequence: str,
    *,
    sequence_name: str = "Selected Sequence",
    sequence_type: Optional[str] = None,
    source: str = "Selected Sequence",
    selected_orf: Optional[Any] = None,
    program: str = "blastn",
    database: str = "nt",
    evalue: float = 10.0,
    max_target_seqs: int = DEFAULT_MAX_HITS,
    email: str = DEFAULT_EMAIL,
    tool: str = DEFAULT_TOOL,
    engine: str = "remote",
    local_database_path: Optional[str | Path] = None,
    local_database_root: Optional[str | Path] = None,
    client: Optional[NCBIBlastClient] = None,
) -> BlastRun:
    """
    Programmatic Synapse entry point.

    If selected_orf is supplied, the ORF becomes the BLAST query while global
    coordinates are retained in the query metadata.

    engine:
        "remote" -> NCBI BLAST URL API
        "local"  -> installed NCBI BLAST+ against a local database
    """
    if selected_orf is not None:
        query = build_orf_query(
            selected_orf,
            sequence_name=sequence_name,
            parent_sequence=sequence,
        )
    else:
        query = build_query(
            sequence,
            name=sequence_name,
            sequence_type=sequence_type,
            source=source,
        )

    parameters = BlastParameters(
        program=program,
        database=database,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        email=email,
        tool=tool,
    )

    engine_normalized = engine.lower().strip()
    if engine_normalized == "remote":
        parameters = parameters.normalized()
        api_client = client or NCBIBlastClient(
            tool=parameters.tool,
            email=parameters.email,
        )
        return api_client.run(query, parameters)

    if engine_normalized != "local":
        raise BlastAnalysisError("engine must be 'local' or 'remote'.")

    manager = LocalBlastDatabaseManager(root=local_database_root)
    if local_database_path:
        # A BLAST database is addressed by a path prefix (".../name" for
        # name.nin, name.nsq, ...), not by a file (audit D6). The prefix is
        # scanned and verified like a managed database, which also records
        # its type, title, date and size (S7).
        prefix = Path(local_database_path).expanduser().resolve()
        manager = LocalBlastDatabaseManager(root=prefix.parent)
        database_record = manager.verify_database(prefix.name)
        if not database_record.exists:
            raise BlastAnalysisError(
                f"No BLAST database found at prefix {prefix} "
                f"({database_record.verification_message or database_record.verification_status})."
            )
        database_record.source = "User-supplied path"
    else:
        database_record = manager.scan_database(database)

    if not database_record.exists:
        raise BlastAnalysisError(
            f"Local BLAST database '{database}' is not available under "
            f"{manager.root}."
        )

    if database_record.database_type not in {"nucleotide", "protein"}:
        raise BlastAnalysisError("Local database type is invalid.")

    expected_db_type = PROGRAMS.get(program.lower().strip(), {}).get("database_type")
    if expected_db_type != database_record.database_type:
        raise BlastAnalysisError(
            f"{program} requires a {expected_db_type} database, "
            f"but {database_record.name} is {database_record.database_type}."
        )

    # Local validation is intentionally independent from the remote catalog.
    if parameters.evalue <= 0 or parameters.max_target_seqs < 1:
        raise BlastAnalysisError("Invalid BLAST parameters.")

    engine_client = LocalBlastEngine()
    return engine_client.run(query, parameters, database_record)
