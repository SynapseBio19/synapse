"""Terminal dialogs of analyses/blast.py (moved here in P2-2).

The scientific functions stay in ``analyses.blast``; this module draws the
screens, reads the answers and calls them.
"""

from __future__ import annotations

from core.alphabets import AMBIGUOUS
from core.alphabets import classify_sequence
from typing import Any
from typing import Optional
import json
from analyses.blast import (  # noqa: F401
    BLAST_UPDATE_SCRIPT,
    BlastAnalysisError,
    BlastHit,
    BlastParameters,
    BlastQuery,
    BlastRun,
    DATABASES,
    DEFAULT_LOCAL_DB_ROOT,
    DEFAULT_MAX_HITS,
    LOCAL_DATABASE_CATALOG,
    LOCAL_DB_ROOT_ENV,
    LocalBlastDatabase,
    LocalBlastDatabaseManager,
    LocalBlastEngine,
    NCBIBlastClient,
    PROGRAMS,
    _evalue_text,
    _infer_scientific_name_from_definition,
    alignment_text,
    build_evidence_profile,
    build_orf_query,
    build_query,
    build_text_report,
    compare_blast_runs,
    export_csv,
    export_json,
    export_text,
    group_hits_by_taxonomy,
    infer_sequence_type,
    load_fasta_queries,
    programs_for_query_type,
    rank_hits,
    run_self_test,
)
from ui import style
from ui.style import minor_rule


def _taxonomy_source(hit: BlastHit) -> str:
    """Return how Synapse obtained the displayed taxonomy label."""
    if hit.taxonomy_ids:
        return "NCBI BLAST taxonomy fields"
    if hit.scientific_names:
        return "NCBI BLAST scientific-name field"
    if _infer_scientific_name_from_definition(hit.definition):
        return "inferred from BLAST hit definition"
    return "not available"


def print_hit_summary(run: BlastRun, *, limit: int = 15) -> None:
    style.header('BLAST RESULT SUMMARY')
    style.field('Query', f"{run.query.name}")
    style.field('Type', f"{run.query.sequence_type}")
    style.field('Length', f"{run.query.length}")
    style.field('Program', f"{run.metadata.program}")
    style.field('Database', f"{run.metadata.database}")
    style.field("RID", run.metadata.rid)
    style.field('Hits', f"{run.hit_count}")
    style.field('Significant hits', f"{run.significant_hit_count}")
    print(minor_rule())

    widths = (5, 17, 9, 9, 10)
    print(style.columns(("Rank", "Accession", "Identity", "Coverage", "E-value", "Bit score"), widths))
    print(minor_rule())

    for rank, hit in enumerate(rank_hits(run.hits)[:limit], start=1):
        print(style.columns((
            rank,
            (hit.accession or hit.hit_id)[:17],
            f"{hit.best_identity_percent:.2f}%",
            f"{hit.best_coverage_percent:.2f}%",
            _evalue_text(hit.best_evalue),
            f"{hit.best_bit_score:.2f}",
        ), widths))


def print_hit_detail(hit: BlastHit) -> None:
    style.header('BLAST HIT DETAIL')
    style.field('Accession', f"{hit.accession}")
    style.field('Hit ID', f"{hit.hit_id}")
    style.field('Definition', f"{hit.definition}")
    style.field('Subject length', f"{hit.length}")
    style.field('HSPs', f"{len(hit.hsps)}")
    print()

    for hsp in hit.hsps:
        print(
            f"HSP {hsp.number}: identity={hsp.identity_percent:.2f}% | "
            f"coverage={hsp.query_coverage:.2f}% | "
            f"E-value={hsp.evalue:.3g} | "
            f"bit score={hsp.bit_score:.2f}"
        )


def _print_local_database_status(manager: LocalBlastDatabaseManager) -> None:
    style.header("Local BLAST database status")
    style.field('Root', f"{manager.root}")
    print(minor_rule())
    widths = (16, 11, 22, 8)
    print(style.columns(("Database", "Type", "Status", "Volumes", "Size"), widths))
    print(minor_rule())

    for database in manager.list_databases():
        if database.exists:
            status = database.verification_status
            size = (
                f"{database.size_gb:.2f} GB"
                if database.size_gb >= 1
                else f"{database.size_mb:.1f} MB"
            )
        elif database.installed:
            if database.missing_volume_files:
                first_volume, first_missing = next(iter(database.missing_volume_files.items()))
                status = f"INCOMPLETE ({first_volume}: {','.join(first_missing)})"
            elif database.missing_volumes:
                status = f"INCOMPLETE (missing {len(database.missing_volumes)})"
            else:
                status = "INCOMPLETE"
            size = (
                f"{database.size_gb:.2f} GB"
                if database.size_gb >= 1
                else f"{database.size_mb:.1f} MB"
            )
        else:
            status = "NOT INSTALLED"
            size = "-"

        if database.expected_volume_count is not None:
            volume_text = (
                f"{database.volume_count}/{database.expected_volume_count}"
            )
        elif database.volume_count:
            volume_text = str(database.volume_count)
        else:
            volume_text = "1" if database.installed else "-"

        print(style.columns(
            (database.name, database.database_type, status, volume_text, size), widths))


def _ensure_local_database_verified(
    manager: LocalBlastDatabaseManager,
    database: LocalBlastDatabase,
) -> LocalBlastDatabase:
    """Ensure a selected local database passes a real BLASTDB verification.

    Structural discovery alone is not enough for a multi-volume NCBI database.
    If the database has not yet been verified, run ``blastdbcmd -info`` through
    the same runtime-view mechanism used by the database manager. This prevents
    a search from starting against a database that Synapse can see on disk but
    BLAST+ cannot actually open.
    """
    if database.verification_status == "VERIFIED":
        return database

    verified = manager.verify_database(database.name, run_blastdbcmd=True)
    if verified.verification_status != "VERIFIED":
        message = verified.verification_message or "No verification diagnostic was returned."
        raise BlastAnalysisError(
            f"Local BLAST database '{database.name}' is not verified. {message}"
        )
    return verified


def _select_local_database(
    manager: LocalBlastDatabaseManager,
    database_type: str,
) -> Optional[LocalBlastDatabase]:
    available = [
        database
        for database in manager.list_databases()
        if database.database_type == database_type and database.exists
    ]

    if not available:
        print()
        style.note("No compatible local BLAST database is installed.")
        style.note("Use [4] BLAST Database Manager to download or create one.")
        return None

    # Prefer a fully verified database and, among equal states, prefer the
    # canonical NCBI RefSeq database over small test/custom databases.
    status_rank = {
        "VERIFIED": 0,
        "STRUCTURALLY_VALID": 1,
        "USER_SUPPLIED": 2,
    }
    available.sort(
        key=lambda database: (
            status_rank.get(database.verification_status, 3),
            0 if database.name.lower() == "refseq_rna" else 1,
            -database.size_bytes,
            database.name.lower(),
        )
    )

    entries = []
    for index, database in enumerate(available, start=1):
        if database.verification_status == "VERIFIED":
            status_label = "VERIFIED"
        else:
            status_label = database.verification_status
        sequence_info = (
            f" | {database.sequence_count:,} seq"
            if database.sequence_count is not None
            else ""
        )
        entries.append((
            index,
            f"{database.name} | "
            f"{database.size_gb:.2f} GB | {status_label}{sequence_info}",
        ))
    style.menu([("Local database", entries)], zero=style.CANCEL_LABEL)

    choice = style.ask('Database')
    if choice == "0":
        return None

    try:
        return available[int(choice) - 1]
    except (ValueError, IndexError):
        style.invalid_selection()
        return None


def blast_database_manager_menu(
    manager: Optional[LocalBlastDatabaseManager] = None,
) -> None:
    """Interactive manager for local BLAST+ databases."""
    manager = manager or LocalBlastDatabaseManager()

    while True:
        style.header('BLAST DATABASE MANAGER')
        style.field('Database root', f"{manager.root}")
        style.menu([
            ('Database', [
                ('1', 'Database Status'),
                ('2', 'Download / Update NCBI Database'),
                ('3', 'Verify Database'),
                ('4', 'Create Custom Database from FASTA'),
                ('5', 'Remove Database'),
                ('6', 'Database Information'),
            ]),
            ('Tools', [('7', 'BLAST+ Tool Status'), ('8', 'Change Database Root')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return

        if choice == "1":
            _print_local_database_status(manager)
            continue

        if choice == "2":
            entries = []
            for index, (name, info) in enumerate(
                LOCAL_DATABASE_CATALOG.items(), start=1
            ):
                recommendation = "RECOMMENDED START" if info["recommended"] else "ADVANCED"
                entries.append((
                    index,
                    f"{name:<18} "
                    f"{info['type']:<12} "
                    f"{recommendation:<16} "
                    f"{info['description']}",
                ))
            style.menu([("NCBI databases available through update_blastdb.pl", entries)],
                       zero=style.CANCEL_LABEL)

            raw = style.ask('Database')
            if raw == "0":
                continue

            names = list(LOCAL_DATABASE_CATALOG)
            try:
                name = names[int(raw) - 1]
            except (ValueError, IndexError):
                style.invalid_selection()
                continue

            print()
            style.field('Selected database', f"{name}")
            style.warning(
                "The database will be downloaded from NCBI into the "
                f"following directory:\r\n{manager.root}"
            )
            print(
                "Large databases can require substantial disk space and network "
                "bandwidth. Synapse does not start this automatically."
            )
            confirm = style.ask('Type DOWNLOAD to continue')
            if confirm != "DOWNLOAD":
                style.note("Download cancelled.")
                continue

            try:
                database = manager.update_ncbi_database(name)
                print()
                style.success(f"Update completed: {database.name}")
                style.field('Status', f"{database.verification_status}")
                style.field('Files', f"{database.file_count}")
                style.field('Size', f"{database.size_gb:.3f} GB")
            except BlastAnalysisError as exc:
                style.error(f"{exc}")
            continue

        if choice == "3":
            _print_local_database_status(manager)
            print()
            name = style.ask('Database name to verify (0 cancel)')
            if name == "0":
                continue
            try:
                database = manager.verify_database(name, run_blastdbcmd=True)
                style.section('VERIFICATION RESULT')
                style.field('Database', f"{database.name}")
                style.field('Status', f"{database.verification_status}")
                style.field('Message', f"{database.verification_message}")
            except BlastAnalysisError as exc:
                style.error(f"{exc}")
            continue

        if choice == "4":
            fasta = style.ask('FASTA file path')
            if not fasta:
                style.error("No FASTA path supplied.")
                continue

            name = style.ask('New database name')
            if not name:
                style.error("No database name supplied.")
                continue

            style.menu([
                ('Database Type', [('1', 'Nucleotide'), ('2', 'Protein')]),
            ], zero="Cancel")
            db_choice = style.ask('Type')
            if db_choice == "0":
                continue
            database_type = {"1": "nucleotide", "2": "protein"}.get(db_choice)
            if database_type is None:
                style.invalid_selection()
                continue

            try:
                database = manager.create_custom_database(
                    fasta,
                    name,
                    database_type,
                )
                style.section('CUSTOM DATABASE CREATED')
                style.field('Name', f"{database.name}")
                style.field('Type', f"{database.database_type}")
                style.field('Path', f"{database.path}")
                style.field('Files', f"{database.file_count}")
                style.field('Size', f"{database.size_mb:.2f} MB")
            except BlastAnalysisError as exc:
                style.error(f"{exc}")
            continue

        if choice == "5":
            _print_local_database_status(manager)
            print()
            name = style.ask('Database name to remove (0 cancel)')
            if name == "0":
                continue

            database = manager.scan_database(name)
            if not database.installed:
                style.error("Database is not installed.")
                continue

            confirm = style.ask(f"Type REMOVE {database.name} to permanently remove its local files")
            if confirm != f"REMOVE {database.name}":
                style.note("Removal cancelled.")
                continue

            try:
                removed = manager.remove_database(name)
                style.success(f"Removed {len(removed)} files.")
            except BlastAnalysisError as exc:
                style.error(f"{exc}")
            continue

        if choice == "6":
            _print_local_database_status(manager)
            print()
            name = style.ask('Database name (0 cancel)')
            if name == "0":
                continue
            try:
                info = manager.database_info(name)
                style.section('DATABASE INFORMATION')
                print(json.dumps(info, indent=2, ensure_ascii=False))
            except BlastAnalysisError as exc:
                style.error(f"{exc}")
            continue

        if choice == "7":
            style.section('BLAST+ TOOL STATUS')
            tools = manager.executable_manager.available()
            for name, path in tools.items():
                style.field(name, path or "NOT FOUND")
            style.field("perl", manager.executable_manager.resolve_perl() or "NOT FOUND")
            print(
                f"{BLAST_UPDATE_SCRIPT:<14}: "
                f"{manager.executable_manager.resolve_update_script() or 'NOT FOUND'}"
            )
            continue

        if choice == "8":
            raw_root = style.ask('New database root (blank = default project data/blast_db, 0 cancel)')
            if raw_root == "0":
                continue
            try:
                manager = LocalBlastDatabaseManager(root=raw_root or DEFAULT_LOCAL_DB_ROOT)
                manager.ensure_root()
                style.field('Database root changed to', f"{manager.root}")
                print(
                    f"To make this persistent across sessions, set the "
                    f"{LOCAL_DB_ROOT_ENV} environment variable."
                )
            except OSError as exc:
                style.error(f"{exc}")
            continue

        style.invalid_selection()


def _query_type_for_menu(sequence: str, declared: Optional[str]) -> Optional[str]:
    """Query type for the menu: declared, inferred, or asked when undetermined.

    Returns None when the user cancels. Invalid sequences raise.
    """
    try:
        return infer_sequence_type(sequence, declared)
    except BlastAnalysisError as exc:
        if classify_sequence(sequence, declared).evidence != AMBIGUOUS:
            raise
        print(f"\r\n{exc}")
    style.menu([
        (None, [('1', 'DNA'), ('2', 'Protein')]),
    ], zero="Cancel")
    # An unknown key is answered and asked again, as in every other menu;
    # until 2026-11-12 it cancelled like 0 (TERMINAL_UI_STANDARD 9).
    while True:
        choice = style.ask('Query type')
        if choice == "0":
            return None
        if choice in ("1", "2"):
            return {"1": "DNA", "2": "PROTEIN"}[choice]
        style.invalid_selection()


def blast_menu(
    sequence: Optional[str] = None,
    sequence_name: str = "Selected Sequence",
    selected_orf: Optional[Any] = None,
    sequence_type: Optional[str] = None,
) -> Optional[BlastRun]:
    """
    Interactive Synapse BLAST Search menu.

    Synapse supports two explicit search engines:
        1. Local BLAST+ against a managed local database.
        2. NCBI remote BLAST through the NCBI URL API.

    If an ORF is selected, it is used as the BLAST query context and its
    global coordinates are retained. A user can also load an external FASTA
    file as an explicit query source; no implicit first-record-only behavior
    is used for multi-record FASTA files.

    *sequence_type* is the selected sequence's declared type, if known; a
    sequence whose type cannot be told from its letters is asked for.
    """
    base_sequence = sequence
    current_name = sequence_name
    current_orf = selected_orf
    active_query: Optional[BlastQuery] = None
    last_run: Optional[BlastRun] = None
    comparison_run: Optional[BlastRun] = None
    manager = LocalBlastDatabaseManager()

    while True:
        style.header('BLAST SEARCH')

        if active_query is not None:
            style.field('Query context', f"{active_query.name}")
            style.field('Query source', f"{active_query.source}")
            style.field('Query length', f"{active_query.length} {active_query.sequence_type.lower()}")
        elif current_orf is not None:
            try:
                orf_query = build_orf_query(current_orf, sequence_name=current_name, parent_sequence=base_sequence)
                style.field('Query context', f"Selected ORF ({orf_query.name})")
                style.field('Query length', f"{orf_query.length} bp")
            except BlastAnalysisError:
                style.field('Query context', "Selected ORF (unavailable)")
        elif base_sequence:
            style.field('Query context', f"{current_name}")
            style.field('Query source', "Selected Synapse Sequence")
            style.field('Query length', f"{len(''.join(base_sequence.split()))} nt")
        else:
            style.field('Query context', "No sequence loaded")

        if last_run:
            style.field('Last engine', f"{last_run.metadata.engine}")
            style.field('Last database', f"{last_run.metadata.database}")

        style.menu([
            ('Search', [
                ('1', 'Run BLAST Search'),
                ('2', 'Resume Search by Request ID (RID)'),
                ('3', 'Search Configuration'),
                ('4', 'BLAST Database Manager'),
            ]),
            ('Result Analysis', [
                ('5', 'Result Summary'),
                ('6', 'Hit Detail & Alignment'),
                ('7', 'Evidence Profile'),
                ('8', 'Taxonomic Overview'),
                ('9', 'Compare BLAST Runs'),
            ]),
            ('Output', [('10', 'Export Results'), ('11', 'Full BLAST Report')]),
            ('System', [
                ('12', 'BLAST Self-Test'),
                ('13', 'Supported Programs & Databases'),
                ('14', 'Database Status'),
            ]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            # Returned to main.py so the run can be stored in the pipeline.
            return last_run

        if choice == "1":
            style.menu([
                ('Query Source', [
                    ('1', 'Selected Synapse Sequence'),
                    ('2', 'Load FASTA File'),
                ]),
            ], zero="Cancel")
            query_source_choice = style.ask('Source')

            if query_source_choice == "0":
                continue

            search_query: Optional[BlastQuery] = None

            if query_source_choice == "1":
                active_query = None
                if not base_sequence and current_orf is None:
                    style.error("No DNA/protein sequence or ORF is available.")
                    style.pause()
                    continue

                try:
                    if current_orf is not None:
                        search_query = build_orf_query(
                            current_orf,
                            sequence_name=current_name,
                            parent_sequence=base_sequence,
                        )
                    else:
                        selected_type = _query_type_for_menu(base_sequence or "", sequence_type)
                        if selected_type is None:
                            continue
                        search_query = build_query(
                            base_sequence or "",
                            name=current_name,
                            sequence_type=selected_type,
                            source="Selected Synapse Sequence",
                        )
                except BlastAnalysisError as exc:
                    style.error(f"{exc}")
                    style.pause()
                    continue

            elif query_source_choice == "2":
                raw_path = style.ask('FASTA file path').strip('\"')
                if not raw_path:
                    style.error("No FASTA path supplied.")
                    continue

                try:
                    fasta_queries = load_fasta_queries(raw_path)
                except BlastAnalysisError as exc:
                    style.error(f"{exc}")
                    style.pause()
                    continue

                if len(fasta_queries) == 1:
                    search_query = fasta_queries[0]
                else:
                    style.menu([("FASTA records", [
                        (index,
                         f"{fasta_query.name} | "
                         f"{fasta_query.length} {fasta_query.sequence_type.lower()}")
                        for index, fasta_query in enumerate(fasta_queries, start=1)
                    ])], zero=style.CANCEL_LABEL)
                    record_choice = style.ask('Record')
                    if record_choice == "0":
                        continue
                    try:
                        search_query = fasta_queries[int(record_choice) - 1]
                    except (ValueError, IndexError):
                        style.invalid_selection()
                        continue

                active_query = search_query
                style.field('Loaded FASTA query', f"{search_query.name} "
                    f"({search_query.length} {search_query.sequence_type.lower()})")

            else:
                style.invalid_selection()
                continue

            if search_query is None:
                style.error("No BLAST query was selected.")
                style.pause()
                continue

            available_query_type = search_query.sequence_type

            compatible_programs = programs_for_query_type(available_query_type)

            style.menu([
                ('Search Engine', [('1', 'Local BLAST+'), ('2', 'NCBI BLAST')]),
            ], zero="Cancel")
            engine_choice = style.ask('Engine')
            if engine_choice == "0":
                continue

            engine = {"1": "local", "2": "remote"}.get(engine_choice)
            if engine is None:
                style.invalid_selection()
                continue

            style.menu([("BLAST program", [
                (index, f"{name.upper()} - {PROGRAMS[name]['description']}")
                for index, name in enumerate(compatible_programs, start=1)
            ])], zero=style.CANCEL_LABEL)

            program_choice = style.ask('Program')
            if program_choice == "0":
                continue

            try:
                program_index = int(program_choice)
                program = compatible_programs[program_index - 1]
            except (ValueError, IndexError):
                style.invalid_selection()
                continue

            database: Optional[str] = None
            local_database: Optional[LocalBlastDatabase] = None

            if engine == "local":
                database_type = PROGRAMS[program]["database_type"]
                local_database = _select_local_database(manager, database_type)
                if local_database is None:
                    continue
                try:
                    local_database = _ensure_local_database_verified(
                        manager, local_database
                    )
                except BlastAnalysisError as exc:
                    style.error(f"{exc}")
                    style.pause()
                    continue
                database = local_database.name
                print(
                    f"Verified local database: {local_database.name} | "
                    f"{local_database.sequence_count:,} sequences"
                    if local_database.sequence_count is not None
                    else f"Verified local database: {local_database.name}"
                )
            else:
                compatible_databases = [
                    name
                    for name, info in DATABASES.items()
                    if info["type"] == PROGRAMS[program]["database_type"]
                ]

                style.menu([("NCBI remote database", [
                    (index, f"{name} - {DATABASES[name]['description']}")
                    for index, name in enumerate(compatible_databases, start=1)
                ])], zero=style.CANCEL_LABEL)

                db_choice = style.ask('Database')
                if db_choice == "0":
                    continue

                try:
                    db_index = int(db_choice)
                    database = compatible_databases[db_index - 1]
                except (ValueError, IndexError):
                    style.invalid_selection()
                    continue

            raw_evalue = style.ask('E-value [10]')
            raw_max = style.ask(f"Maximum target sequences [{DEFAULT_MAX_HITS}]")

            raw_threads = ""
            raw_word_size = ""
            raw_filter = ""
            if engine == "local":
                raw_threads = style.ask('CPU threads [4]')
                raw_word_size = style.ask('Word size [default]')
                raw_filter = style.ask('Filtering [default/none]').lower()

            try:
                evalue = float(raw_evalue) if raw_evalue else 10.0
                max_hits = int(raw_max) if raw_max else DEFAULT_MAX_HITS
                threads = int(raw_threads) if raw_threads else 4
                word_size = int(raw_word_size) if raw_word_size else None
                filter_query = None
                if raw_filter in {"none", "no", "off"}:
                    filter_query = "no"
                elif raw_filter in {"default", "", "yes", "on"}:
                    filter_query = None
                else:
                    raise BlastAnalysisError(
                        "Filtering must be 'default' or 'none'."
                    )
                if threads < 1:
                    raise BlastAnalysisError("CPU threads must be at least 1.")
                if word_size is not None and word_size < 1:
                    raise BlastAnalysisError("Word size must be at least 1.")

                query = search_query
                if query is None:
                    raise BlastAnalysisError("No BLAST query was selected.")

                parameters = BlastParameters(
                    program=program,
                    database=database or "",
                    evalue=evalue,
                    max_target_seqs=max_hits,
                    word_size=word_size,
                    filter_query=filter_query,
                    num_threads=threads if engine == "local" else None,
                )

                if engine == "local":
                    assert local_database is not None
                    print("\nRunning local BLAST+ search...")
                    new_run = LocalBlastEngine().run(
                        query,
                        parameters,
                        local_database,
                    )
                else:
                    normalized = parameters.normalized()
                    style.section("PRIVACY NOTICE")
                    print(
                        "The query sequence will be sent to NCBI (blast.ncbi.nlm.nih.gov) "
                        "and processed on NCBI servers. Do not submit confidential or "
                        "unpublished sequences unless this is acceptable."
                    )
                    if not style.confirm("Send this query to NCBI?"):
                        style.note("Remote search cancelled.")
                        continue
                    print("\nSubmitting BLAST search to NCBI...")
                    print("The search may take several seconds or longer.")
                    client = NCBIBlastClient(
                        tool=normalized.tool,
                        email=normalized.email,
                    )
                    new_run = client.run(query, normalized)

                comparison_run = last_run
                last_run = new_run
                print()
                style.success("BLAST search completed.")
                print_hit_summary(last_run)
                if not last_run.hits:
                    print()
                    style.note("No hits were returned. This is a valid BLAST result; it does not indicate that the query is biologically irrelevant.")

            except (BlastAnalysisError, AssertionError) as exc:
                style.error(f"{exc}")

        elif choice == "3":
            style.section('SEARCH CONFIGURATION')
            print("Supported programs:")
            for name, info in PROGRAMS.items():
                print(
                    f"- {name.upper():<8} "
                    f"query={info['query_type']:<8} "
                    f"database={info['database_type']:<11}"
                )
            print()
            style.field('Default E-value', f"{10.0}")
            style.field('Default max target seqs', f"{DEFAULT_MAX_HITS}")
            print("Default local CPU threads: 4")
            style.field('Local filtering', "BLAST+ default")
            style.field('Local database root', f"{manager.root}")
            print(
                "NCBI remote searches use the NCBI URL API. "
                "Local searches use installed BLAST+ executables."
            )

        elif choice == "4":
            blast_database_manager_menu(manager)
            manager = LocalBlastDatabaseManager(root=manager.root)

        elif choice == "5":
            if not last_run:
                style.note("Run [1] BLAST Search first.")
                style.pause()
                continue
            print_hit_summary(last_run)

        elif choice == "6":
            if not last_run or not last_run.hits:
                style.note("No BLAST hits available.")
                style.pause()
                continue

            ranked = rank_hits(last_run.hits)
            style.menu([("Ranked hits", [
                (index,
                 f"{hit.accession or hit.hit_id} | "
                 f"{hit.best_identity_percent:.2f}% identity | "
                 f"{hit.best_coverage_percent:.2f}% coverage")
                for index, hit in enumerate(ranked[:20], start=1)
            ])], zero=style.CANCEL_LABEL)

            raw = style.ask('Hit number')
            if raw == "0":
                continue

            # Only the listed numbers select a hit. Before 2026-11-12 a
            # larger number gave "list index out of range" and a negative
            # one silently showed a hit from the end of the list.
            shown = min(len(ranked), 20)
            if not raw.isdigit() or not 1 <= int(raw) <= shown:
                style.invalid_selection()
                continue

            try:
                index = int(raw)
                hit = ranked[index - 1]
                print_hit_detail(hit)
                print()
                print(alignment_text(hit))
            except (ValueError, IndexError, BlastAnalysisError) as exc:
                style.error(f"{exc}")

        elif choice == "7":
            if not last_run:
                style.note("Run [1] BLAST Search first.")
                style.pause()
                continue

            style.section('EVIDENCE PROFILE')
            for index, profile in enumerate(build_evidence_profile(last_run), start=1):
                print(
                    f"{index}. {profile.accession or profile.hit_id}: "
                    f"{profile.evidence_score:.2f}/100 | "
                    f"{profile.evidence_band}"
                )
                for note in profile.notes:
                    print(f"    - {note}")

        elif choice == "8":
            if not last_run:
                style.note("Run [1] BLAST Search first.")
                style.pause()
                continue

            groups = group_hits_by_taxonomy(last_run.hits)
            style.section('TAXONOMIC OVERVIEW')

            identified = sum(1 for hit in last_run.hits if hit.scientific_names)
            unknown = len(last_run.hits) - identified
            inferred = sum(
                1
                for hit in last_run.hits
                if hit.scientific_names
                and not hit.taxonomy_ids
                and _taxonomy_source(hit) == "inferred from BLAST hit definition"
            )

            style.field('Hits analyzed', f"{len(last_run.hits)}")
            style.field('Taxa identified', f"{identified}")
            style.field('Unknown', f"{unknown}")
            if inferred:
                style.field('Definition-inferred', f"{inferred} "
                    "(not an authoritative taxonomy assignment)")
            print()

            if not groups:
                style.note("No taxonomy labels were available from this result set.")
            else:
                print(style.columns(("Taxon / organism", "Hits"), (40,)))
                print(minor_rule())
                for name, count in groups.items():
                    print(style.columns((name, count), (40,)))

            if inferred:
                print()
                style.note("Some organism names were inferred conservatively "
                    "from BLAST hit definitions because explicit taxonomy "
                    "fields were not present in the returned result.")

        elif choice == "9":
            if not last_run:
                style.note("Run [1] BLAST Search first.")
                style.pause()
                continue

            if comparison_run is None:
                print(
                    "A second BLAST run is required. "
                    "Run another search first."
                )
                continue

            comparison = compare_blast_runs(comparison_run, last_run)
            style.section('BLAST RUN COMPARISON')
            style.field('Query A', f"{comparison.query_a}")
            style.field('Query B', f"{comparison.query_b}")
            style.field('Common hits', f"{len(comparison.common_accessions)}")
            style.field('Unique to A', f"{len(comparison.unique_to_a)}")
            style.field('Unique to B', f"{len(comparison.unique_to_b)}")

        elif choice == "10":
            if not last_run:
                style.note("Run [1] BLAST Search first.")
                style.pause()
                continue

            path = style.ask('Export path (.json/.csv/.txt)')
            if not path:
                style.error("No path supplied.")
                continue

            try:
                lower = path.lower()
                if lower.endswith(".json"):
                    export_json(last_run, path)
                elif lower.endswith(".csv"):
                    export_csv(last_run, path)
                elif lower.endswith(".txt"):
                    export_text(last_run, path)
                else:
                    style.error("Use .json, .csv or .txt.")
                    continue
                style.success(f"Exported: {path}")
            except OSError as exc:
                style.error(f"Export failed: {exc}")

        elif choice == "11":
            if not last_run:
                style.note("Run [1] BLAST Search first.")
            else:
                print(build_text_report(last_run))

        elif choice == "12":
            run_self_test(verbose=True)

        elif choice == "13":
            style.section('BLAST PROGRAMS')
            for name, info in PROGRAMS.items():
                print(
                    f"{name.upper():<9} "
                    f"query={info['query_type']:<8} "
                    f"database={info['database_type']:<11} "
                    f"{info['description']}"
                )

            style.section('REMOTE DATABASES')
            for name, info in DATABASES.items():
                print(f"{name:<18} {info['type']:<12} {info['description']}")

            style.section('LOCAL NCBI DATABASE CATALOG')
            for name, info in LOCAL_DATABASE_CATALOG.items():
                print(
                    f"{name:<18} {info['type']:<12} "
                    f"{info['description']}"
                )

        elif choice == "14":
            _print_local_database_status(manager)

        elif choice == "2":
            # P1-7: collect a search submitted earlier instead of repeating it.
            style.section("RESUME A PREVIOUS SEARCH")
            print(
                "Enter the request ID (RID) of a search submitted earlier. Only "
                "the RID is sent to NCBI; the query sequence is not transmitted "
                "again. NCBI keeps a request ID for a limited time, so an older "
                "one may have expired."
            )
            rid = style.ask('Request ID (RID)')
            if not rid:
                style.error("No request ID supplied.")
                continue

            if not style.confirm(f"Contact NCBI for request {rid}?"):
                style.note("Resume cancelled.")
                continue

            try:
                resumed = NCBIBlastClient().resume(rid)
            except BlastAnalysisError as exc:
                style.error(f"{exc}")
                style.pause()
                continue

            last_run = resumed
            active_query = resumed.query
            style.success(
                f"Resumed request {rid}: {len(resumed.hits)} hit(s) for "
                f"{resumed.query.name} ({resumed.query.length} "
                f"{resumed.query.sequence_type.lower()})."
            )
            for note in resumed.metadata.warnings:
                style.note(f"{note}")

        else:
            style.invalid_selection()
