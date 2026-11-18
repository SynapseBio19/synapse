"""Running a search with local BLAST+ programs."""

from __future__ import annotations

from contextlib import ExitStack
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from analyses._blast.config import DEFAULT_MAX_WAIT_SECONDS, PROGRAMS
from analyses._blast.models import (
    BlastAnalysisError,
    BlastParameters,
    BlastQuery,
    BlastRun,
    LocalBlastDatabase,
    _utc_now,
)
from analyses._blast.query import check_program_query_type
from analyses._blast.parsing import parse_blast_xml
from analyses._blast.executables import BlastExecutableManager
from analyses._blast.runtime import (
    _build_local_blast_environment,
    _build_local_database_command_reference,
    _create_runtime_blast_alias,
    _needs_runtime_blast_alias,
)


def _should_run_short_query_rescue(
    query: BlastQuery,
    parameters: BlastParameters,
) -> bool:
    """Return whether Synapse may add one automatic short-query rescue run."""
    return (
        parameters.program.lower().strip() == "blastn"
        and query.sequence_type == "DNA"
        and query.length <= 200
        and parameters.word_size is None
        and parameters.task is None
        and parameters.filter_query is None
    )


class LocalBlastEngine:
    """Execute BLAST+ locally against an installed database."""

    def __init__(
        self,
        executable_manager: Optional[BlastExecutableManager] = None,
    ) -> None:
        self.executable_manager = executable_manager or BlastExecutableManager()

    def run(
        self,
        query: BlastQuery,
        parameters: BlastParameters,
        database: LocalBlastDatabase,
        *,
        timeout: float = DEFAULT_MAX_WAIT_SECONDS,
        allow_sensitive_retry: bool = True,
    ) -> BlastRun:
        if not database.exists:
            raise BlastAnalysisError(
                f"Local BLAST database '{database.name}' is not installed or incomplete."
            )

        expected_type = PROGRAMS.get(parameters.program.lower().strip(), {}).get(
            "database_type"
        )
        if expected_type and expected_type != database.database_type:
            raise BlastAnalysisError(
                f"{parameters.program} requires a {expected_type} database, "
                f"but '{database.name}' is {database.database_type}."
            )

        program = parameters.program.lower().strip()
        if program not in PROGRAMS:
            raise BlastAnalysisError(f"Unsupported BLAST program: {program}")

        blast_executable = self.executable_manager.require(program)

        # Validate query type against the selected BLAST program.
        check_program_query_type(program, query)

        local_parameters = BlastParameters(
            program=program,
            database=database.name,
            evalue=parameters.evalue,
            max_target_seqs=parameters.max_target_seqs,
            word_size=parameters.word_size,
            gap_open=parameters.gap_open,
            gap_extend=parameters.gap_extend,
            matrix=parameters.matrix,
            filter_query=parameters.filter_query,
            num_threads=parameters.num_threads,
            task=parameters.task,
            tool=parameters.tool,
            email=parameters.email,
        )

        # The ExitStack removes the runtime alias view (below) when the search
        # ends, also after an error. Before 2026-11-18 it was never removed.
        with tempfile.TemporaryDirectory(prefix="synapse_blast_") as temp_dir, ExitStack() as cleanup:
            temp = Path(temp_dir)
            query_file = temp / "query.fasta"
            output_file = temp / "result.xml"

            # CRLF as written, not CR CR LF (D24).
            query_file.write_text(
                f">{query.name}\r\n{query.sequence}\r\n",
                encoding="utf-8",
                newline="",
            )

            # Set only when a runtime alias is needed; the environment below uses it.
            runtime_parent: Optional[Path] = None
            if _needs_runtime_blast_alias(database):
                runtime_parent = Path(database.path).expanduser().resolve().parent
                if re.fullmatch(r".+\.\d+", runtime_parent.name):
                    runtime_parent = runtime_parent.parent
                runtime_directory = Path(
                    tempfile.mkdtemp(
                        prefix="synapse_blast_runtime_",
                        dir=str(runtime_parent),
                    )
                )
                # Hard links only: removing them leaves the database untouched.
                cleanup.callback(shutil.rmtree, runtime_directory, ignore_errors=True)
                database_reference, database_working_directory = _create_runtime_blast_alias(
                    database, runtime_directory
                )
            else:
                database_reference, database_working_directory = (
                    _build_local_database_command_reference(database.path)
                )

            command = [
                blast_executable,
                "-query",
                str(query_file),
                "-db",
                database_reference,
                "-out",
                str(output_file),
                "-outfmt",
                "5",
                "-evalue",
                str(local_parameters.evalue),
                "-max_target_seqs",
                str(local_parameters.max_target_seqs),
            ]

            if local_parameters.num_threads is not None:
                if local_parameters.num_threads < 1:
                    raise BlastAnalysisError("num_threads must be at least 1.")
                command.extend(["-num_threads", str(local_parameters.num_threads)])

            if local_parameters.task:
                command.extend(["-task", local_parameters.task])

            if local_parameters.word_size is not None:
                command.extend(["-word_size", str(local_parameters.word_size)])
            if local_parameters.gap_open is not None:
                command.extend(["-gapopen", str(local_parameters.gap_open)])
            if local_parameters.gap_extend is not None:
                command.extend(["-gapextend", str(local_parameters.gap_extend)])
            if local_parameters.matrix:
                command.extend(["-matrix", local_parameters.matrix])
            if local_parameters.filter_query:
                command.extend(["-seg" if program in {"blastp", "blastx", "tblastn", "tblastx"} else "-dust",
                                local_parameters.filter_query])

            if runtime_parent is not None:
                env = _build_local_blast_environment(
                    runtime_parent
                )
            else:
                env = _build_local_blast_environment(database.path)

            started = _utc_now()
            try:
                completed = subprocess.run(
                    command,
                    cwd=str(database_working_directory),
                    env=env,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    check=False,
                )
            except subprocess.TimeoutExpired as exc:
                raise BlastAnalysisError(
                    f"Local BLAST search exceeded the timeout of {timeout:.0f} seconds."
                ) from exc
            except OSError as exc:
                raise BlastAnalysisError(f"Could not start local BLAST: {exc}") from exc

            stderr = (completed.stderr or "").strip()
            stdout = (completed.stdout or "").strip()
            if completed.returncode != 0:
                detail = stderr or stdout or "No diagnostic output was returned."
                raise BlastAnalysisError(
                    f"Local {program} failed with exit code {completed.returncode}.\r\n"
                    f"{detail[-4000:]}"
                )

            if not output_file.is_file():
                raise BlastAnalysisError(
                    "Local BLAST completed without producing an XML result file."
                )

            xml_text = output_file.read_text(encoding="utf-8", errors="replace")

        rid = f"LOCAL-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
        completed_at = _utc_now()
        run = parse_blast_xml(
            xml_text,
            query=query,
            rid=rid,
            parameters=local_parameters,
            submitted_at=started,
            completed_at=completed_at,
            database_type_override=database.database_type,
        )
        run.metadata.engine = "Local BLAST+"
        run.metadata.local_database_path = database.path
        run.metadata.ncbi_status = "LOCAL_COMPLETE"
        run.metadata.database_details = {
            "name": database.name,
            "title": database.title,
            "database_date": database.database_date,
            "blastdb_version": database.blastdb_version,
            "sequence_count": database.sequence_count,
            "total_length": database.total_length,
            "files_modified_at": database.modified_at,
            "verification_status": database.verification_status,
        }
        if stderr:
            run.metadata.warnings.append(stderr[-2000:])
        run.interpretation_notes.append(
            "This result was generated by the locally installed NCBI BLAST+ executable."
        )
        run.interpretation_notes.append(
            "The local database is a specific installed snapshot; database contents "
            "may differ from the current NCBI remote database."
        )

        # Short DNA queries can miss local similarities with the standard BLASTN
        # seed model. For an otherwise default local BLASTN search, Synapse performs
        # one conservative sensitivity rescue using NCBI's blastn-short task. The
        # reporting E-value is kept unchanged, so the rescue does not silently turn
        # weak matches into significant hits. Explicit user parameters are never
        # overridden.
        short_query_eligible = _should_run_short_query_rescue(query, parameters)

        if allow_sensitive_retry and short_query_eligible and not run.hits:
            rescue_parameters = BlastParameters(
                program=program,
                database=database.name,
                evalue=parameters.evalue,
                max_target_seqs=parameters.max_target_seqs,
                word_size=None,
                gap_open=parameters.gap_open,
                gap_extend=parameters.gap_extend,
                matrix=parameters.matrix,
                filter_query="no",
                num_threads=parameters.num_threads,
                task="blastn-short",
                tool=parameters.tool,
                email=parameters.email,
            )

            rescue_run = self.run(
                query,
                rescue_parameters,
                database,
                timeout=timeout,
                allow_sensitive_retry=False,
            )
            rescue_run.interpretation_notes.insert(
                0,
                "Standard BLASTN returned no hits; Synapse automatically performed "
                "one short-query sensitivity rescue using the blastn-short task.",
            )
            rescue_run.interpretation_notes.insert(
                1,
                "The sensitivity rescue retained the original E-value threshold and "
                "did not override explicit user search parameters.",
            )
            if rescue_run.hits:
                rescue_run.interpretation_notes.insert(
                    2,
                    "Displayed hits are from the sensitive short-query search and "
                    "should be interpreted using identity, coverage and E-value together.",
                )
            else:
                rescue_run.interpretation_notes.insert(
                    2,
                    "No hits were found in either the standard or short-query search.",
                )
            return rescue_run

        return run
