import logging
import os
from datetime import datetime

from core.pipeline import SynapsePipeline
from core.models import SequenceRecord
from core.report import export_analysis_report
from core.logging_setup import configure_logging
from core.paths import example_sequence_file, user_data_dir

from fileio.loader import load_sequences

from core.alphabets import (
    AMBIGUOUS,
    INVALID,
    PROTEIN,
    RNA,
    classify_sequence,
)

from analyses.sequence import (
    analyze_sequence_quality,
    count_nucleotides,
    count_ambiguous_bases,
    calculate_gc_content,
    calculate_at_content,
    calculate_dna_content_ranges,
    calculate_gc_skew,
    calculate_at_skew,
    calculate_nucleotide_percentages,
    get_complement,
    get_reverse_complement,
    calculate_dna_molecular_weight,
    calculate_sequence_complexity,
    find_homopolymers,
    analyze_cpg,
    find_cpg_islands,
)

from analyses.orf import (
    DEFAULT_MIN_ORF_LENGTH,
)

from analyses.literature import (
    build_sequence_literature_context,
)

# The module dialogs live in ui/ (P2-2); the analyses/ modules hold the
# science and do not import the interface.
from ui.rna_menu import run_rna_analysis_ui
from ui.protein_menu import (
    run_protein_analysis_from_orf,
    run_protein_analysis_from_record,
)
from ui.orf_menu import run_orf_analysis_ui
from ui.restriction_menu import run_restriction_cloning_ui
from ui.mutation_menu import run_mutation_variant_ui
from ui.primer_menu import run_primer_analysis_ui
from ui.crispr_menu import crispr_menu
from ui.blast_menu import blast_menu
from ui.literature_menu import main as literature_main

from ui.report import print_pipeline_status, print_results
from ui.menu import select_gene, show_menu
from ui.style import major_rule
from ui import style


LOG_PATH = configure_logging()
logger = logging.getLogger("synapse.main")


def report_error(title, error):
    """Show an actionable error and keep the full traceback in the log file."""
    logger.exception("%s", title)
    style.header(title)
    style.error(f"{type(error).__name__}: {error}")
    style.note(f"Details (traceback) were written to: {LOG_PATH}")
    style.pause()


# ============================================================
# MAIN PROGRAM
# ============================================================


# Synapse's folder for exports and relative paths: the project folder when
# Synapse runs from source, a per-user folder when it is installed
# (core.paths; Phase 9). The example FASTA is installed with the data package.
current_folder = str(user_data_dir())

default_sequence_file = str(example_sequence_file())


# ============================================================
# LOAD SEQUENCE FILE (FASTA or GenBank)
# ============================================================


def choose_sequence_file():
    """Ask for a FASTA/GenBank file; Enter loads the bundled example."""

    while True:

        style.header('LOAD SEQUENCE FILE')
        style.field("Formats", "FASTA (.fasta/.fa) or GenBank (.gb/.gbk)")
        style.field("Topology", "FASTA headers may contain topology=circular")
        style.field("Example", default_sequence_file)
        style.note("Press Enter to load the example; type Q to quit.")
        print()

        raw = style.ask("File path").strip('"')

        if raw.upper() == "Q":
            raise SystemExit

        path = raw or default_sequence_file

        if not os.path.isabs(path):
            candidate = os.path.join(current_folder, path)
            if not os.path.exists(path) and os.path.exists(candidate):
                path = candidate

        try:
            loaded = load_sequences(path)
            logger.info("Loaded %d record(s) from %s", len(loaded), path)
            return loaded

        except (OSError, ValueError) as error:
            logger.warning("Could not load %s: %s", path, error)
            style.header('SEQUENCE FILE ERROR')
            style.error(str(error))


# ============================================================
# SHARED MENU ACTIONS
# ============================================================


def store_blast_run(pipeline, blast_run):
    """Store a finished BLAST run in the pipeline."""
    if blast_run is None:
        return
    from dataclasses import asdict as _asdict

    pipeline.store(
        "blast",
        blast_run,
        parameters=_asdict(blast_run.metadata),
        limitations=[
            "BLAST hits show sequence similarity; they do "
            "not prove function, orthology or identity.",
        ],
    )


def _overview_parameters(topology):
    """The settings of the sequence overview, read from the functions it calls."""
    import inspect

    def default(function, name):
        return inspect.signature(function).parameters[name].default

    return {
        "topology": topology,
        "homopolymer_min_length": default(find_homopolymers, "min_length"),
        "cpg_island_window": default(find_cpg_islands, "window_size"),
        "cpg_island_min_gc_percent": default(find_cpg_islands, "min_gc"),
        "cpg_island_min_observed_expected": default(find_cpg_islands, "min_oe"),
    }


def export_report_interactive(pipeline):
    """Ask for a path and export every stored result."""
    default_report = os.path.join(
        current_folder,
        "exports",
        "synapse_report_"
        + datetime.now().strftime("%Y%m%d_%H%M%S")
        + ".json",
    )
    style.section("Export analysis report")
    style.note("Report of all results stored in this session; use .json (complete) or .md (summary).")
    raw_path = style.ask("Path", default=default_report).strip('"')
    try:
        saved = export_analysis_report(
            pipeline, raw_path or default_report
        )
        style.success(f"Report written: {saved.resolve()}")
    except (OSError, ValueError) as error:
        report_error("REPORT EXPORT ERROR", error)
    style.pause()


# ============================================================
# SEQUENCE TYPE (audit D14)
# ============================================================


def print_sequence_type_error(record):
    """Explain why a record cannot be analyzed."""
    result = classify_sequence(record.sequence, record.declared_type)
    reason = result.reasons[0]
    style.header('SEQUENCE VALIDATION ERROR')
    style.error(f"{reason[0].upper()}{reason[1:]}.")
    print()
    for position, symbol in result.invalid_symbols[:20]:
        print(f"Invalid character '{symbol}' at position {position}")
    if len(result.invalid_symbols) > 20:
        print(f"... and {len(result.invalid_symbols) - 20} more")
    print(major_rule())
    style.pause("Press Enter to select another sequence...")


def resolve_sequence_type(record):
    """Return the record with a determined type, or None to choose another.

    A sequence whose type cannot be told from its letters is not guessed:
    the user declares DNA, RNA or protein.
    """
    while record.type_evidence == AMBIGUOUS:
        style.header('SEQUENCE TYPE')
        style.note("The type of this sequence cannot be determined from its letters: "
                   f"{record.type_reasons[0]}.")
        style.menu([("Declare the type", [(1, "DNA"), (2, "RNA"), (3, "Protein")])],
                   zero="Cancel (choose another sequence)")
        choice = style.ask_choice()
        if choice == "0":
            return None
        declared = {"1": "DNA", "2": "RNA", "3": "PROTEIN"}.get(choice)
        if declared is None:
            print()
            style.invalid_selection()
            continue
        candidate = record.with_declared_type(declared)
        if candidate.type_evidence == INVALID:
            style.error(f"Not possible: {candidate.type_reasons[0]}.")
            continue
        logger.info("Sequence type declared by the user: %s", declared)
        record = candidate
    if record.type_evidence == INVALID:
        print_sequence_type_error(record)
        return None
    return record


def run_protein_record(record, selected_index):
    """Menu for a protein record: protein analysis, BLAST, report."""
    pipeline = SynapsePipeline(
        sequence=record.sequence,
        header=record.description,
        sequence_id=selected_index,
        topology=record.topology,
    )
    pipeline.attach_record(record)

    while True:
        style.header('PROTEIN RECORD')
        style.field("Record", record.description)
        style.field('Length', f"{record.length} aa")
        style.field('Type', f"protein ({record.type_evidence.lower()})")
        style.menu([
            ("Analysis", [(1, "Protein Analysis"), (2, "BLAST Search")]),
            ("Session", [("R", "Export Analysis Report (JSON / Markdown)"),
                         ("S", "Select Another Sequence")]),
        ], zero=style.EXIT_LABEL)

        choice = style.ask_choice().upper()

        if choice == "1":
            try:
                run_protein_analysis_from_record(
                    record.sequence,
                    name=record.identifier,
                    pipeline=pipeline,
                )
            except Exception as error:
                report_error("PROTEIN ANALYSIS ERROR", error)
        elif choice == "2":
            try:
                store_blast_run(
                    pipeline,
                    blast_menu(
                        sequence=record.sequence,
                        sequence_name=record.description,
                        sequence_type="PROTEIN",
                    ),
                )
            except Exception as error:
                report_error("BLAST SEARCH ERROR", error)
        elif choice == "R":
            export_report_interactive(pipeline)
        elif choice == "S":
            style.note("Returning to sequence selection.")
            return
        elif choice == "0":
            print("\nThank you for using Synapse.")
            raise SystemExit
        else:
            print()
            style.invalid_selection()
            style.pause()


# ============================================================
# SESSION STEPS (extracted from main() in Phase 7)
# ============================================================


def compute_sequence_overview(dna):
    """Every value the sequence overview shows and stores, by name.

    The keys are the parameter names of ui.report.print_results; the
    functions run in the order main() used to call them.
    """
    nucleotide_counts = count_nucleotides(dna)
    return {
        "nucleotide_counts": nucleotide_counts,
        "nucleotide_percentages": calculate_nucleotide_percentages(nucleotide_counts),
        "gc_content": calculate_gc_content(dna),
        "at_content": calculate_at_content(dna),
        "content_ranges": calculate_dna_content_ranges(dna),
        "gc_skew": calculate_gc_skew(dna),
        "at_skew": calculate_at_skew(dna),
        "sequence_complexity": calculate_sequence_complexity(dna),
        "homopolymers": find_homopolymers(dna),
        "cpg_analysis": analyze_cpg(dna),
        "cpg_islands": find_cpg_islands(dna),
        "molecular_weight": calculate_dna_molecular_weight(dna),
        "sequence_quality": analyze_sequence_quality(dna),
        "ambiguous_counts": count_ambiguous_bases(dna),
        "complement": get_complement(dna),
        "reverse_complement": get_reverse_complement(dna),
    }


def show_sequence_overview(selected_index, header, dna, topology, overview):
    """Draw the sequence overview from compute_sequence_overview()'s values."""
    print_results(selected_index, header, dna, topology, **overview)


def run_orf_step(pipeline, dna, topology):
    """ORF analysis; a selected ORF is stored in the pipeline and returned."""
    selected_orf = run_orf_analysis_ui(
        dna,
        include_partial=True,
        min_length=DEFAULT_MIN_ORF_LENGTH,
        topology=topology,
    )
    if selected_orf is None:
        return None
    pipeline.set_selected_orf(selected_orf)
    pipeline.store("orf", selected_orf, parameters=selected_orf.get("search_parameters", {}))
    style.header('ORF CONNECTED TO PIPELINE')
    print()
    style.field('Selected', f"ORF-{selected_orf['rank']:03d}")
    style.field('Length', f"{selected_orf['length']} bp")
    style.field('Protein', f"{selected_orf['protein_length']} aa")
    style.field('Location', f"{selected_orf['start'] + 1}-"
                f"{selected_orf.get('stop_on_circle', selected_orf['stop']) + 1} "
                f"({selected_orf['strand']} strand, 1-based)")
    print()
    style.success("ORF stored in the pipeline; RNA, protein and CRISPR guide design can use it.")
    if selected_orf.get("genetic_code", 1) != 1:
        style.note(
            "ORF found with NCBI table "
            f"{selected_orf['genetic_code']}; protein analysis and the RNA "
            "codon view use this table."
        )
    style.pause()
    return selected_orf


def show_protein_requires_orf():
    """Protein analysis needs a selected ORF; there is no whole-sequence fallback."""
    style.header('PROTEIN ANALYSIS')
    print()
    style.note("No ORF has been selected.")
    print()
    style.note("Please run [2] ORF Analysis first and select an ORF.")
    print()
    style.note("Protein Analysis will not use the full DNA sequence as a fallback.")
    style.pause()


def run_protein_step(pipeline, selected_orf):
    """Protein analysis of the selected ORF, translated with the ORF's genetic code."""
    orf_rank = selected_orf.get("rank")
    run_protein_analysis_from_orf(
        orf_sequence=selected_orf.get("sequence", ""),
        orf_id=f"ORF-{orf_rank:03d}" if orf_rank is not None else "Selected ORF",
        strand=selected_orf.get("strand"),
        reading_frame=selected_orf.get("reading_frame", selected_orf.get("frame")),
        orf_status=selected_orf.get("status"),
        pipeline=pipeline,
        genetic_code=selected_orf.get("genetic_code", 1),
    )


def show_planner_status():
    """AI Experiment Planner: COMING SOON (owner decision of 2026-11-13)."""
    print(
        "\nAI Experiment Planner: COMING SOON. This feature is "
        "planned but not yet available. Nothing is run or stored."
    )
    print(
        "It is intentionally not implemented in this phase "
        "(P3-1, docs/AUDIT_STATUS.md)."
    )
    style.pause()



def main():
    """Run the interactive Synapse session.

    Audit S13: this flow used to run at module level, so importing
    ``main`` started the application and tests could only drive it in a
    subprocess. It is a function now, called only from the guard below.
    """

    sequences = choose_sequence_file()


    # ============================================================
    # SEQUENCE SELECTION LOOP
    # ============================================================


    while True:

        selected_index, record = (
            select_gene(sequences)
        )

        if record is None:
            sequences = choose_sequence_file()
            continue

        header = record["header"]

        for note in record.get("notes", []):
            style.note(note)

        # ========================================================
        # SEQUENCE TYPE (DNA, RNA or protein; audit D14)
        # ========================================================

        sequence_record = resolve_sequence_type(
            SequenceRecord.from_legacy_dict(record, selected_index)
        )

        if sequence_record is None:
            continue

        if sequence_record.alphabet == PROTEIN:
            run_protein_record(sequence_record, selected_index)
            continue

        topology = sequence_record.topology
        dna = sequence_record.sequence

        if sequence_record.alphabet == RNA:
            # The nucleotide analyses work in the DNA alphabet.
            dna = dna.replace("U", "T")
            style.note("RNA sequence, analyzed in the DNA alphabet (U read as T).")

        # ========================================================
        # SYNAPSE PIPELINE
        # ========================================================

        pipeline = SynapsePipeline(
            sequence=dna,
            header=header,
            sequence_id=selected_index,
            topology=topology,
        )

        try:
            pipeline.attach_record(sequence_record)
        except ValueError as error:
            logger.warning("Provenance not attached: %s", error)

        # Selected ORF is shared by ORF -> RNA -> Protein.
        selected_orf = None


        # ========================================================
        # SEQUENCE OVERVIEW
        # ========================================================

        overview = compute_sequence_overview(dna)
        show_sequence_overview(selected_index, header, dna, topology, overview)
        pipeline.context.store_sequence_analysis(
            parameters=_overview_parameters(topology), result=dict(overview))


        # ========================================================
        # ANALYSIS MENU
        # ========================================================


        while True:

            choice = (
                show_menu()
                .strip()
                .upper()
            )


            # ====================================================
            # 1 - SEQUENCE OVERVIEW
            # ====================================================


            if choice == "1":

                show_sequence_overview(selected_index, header, dna, topology, overview)

                style.pause()


            # ====================================================
            # 2 - ORF ANALYSIS
            # ====================================================


            elif choice == "2":

                try:
                    orf_result = run_orf_step(pipeline, dna, topology)
                    if orf_result is not None:
                        selected_orf = orf_result
                except Exception as error:
                    report_error("ORF ANALYSIS ERROR", error)


            # ====================================================
            # 3 - RNA ANALYSIS
            # ====================================================


            elif choice == "3":

                try:

                    rna_analysis = run_rna_analysis_ui(
                        dna,
                        pipeline=pipeline,
                        selected_orf=selected_orf,
                        header=header,
                    )

                    # The RNA module prints the pipeline connection itself and
                    # stores the result (core.pipeline.store_rna_analysis).
                    if rna_analysis is not None:
                        logger.info(
                            "RNA analysis stored (source: %s)",
                            rna_analysis.get("pipeline_input", {}).get("source"),
                        )

                except Exception as error:

                    report_error("RNA ANALYSIS ERROR", error)


            # ====================================================
            # 4 - PROTEIN ANALYSIS
            # ====================================================


            elif choice == "4":

                if selected_orf is None:
                    show_protein_requires_orf()
                    continue

                try:
                    run_protein_step(pipeline, selected_orf)
                except Exception as error:
                    report_error("PROTEIN ANALYSIS ERROR", error)


            # ====================================================
            # 5 - RESTRICTION & CLONING
            # ====================================================


            elif choice == "5":

                try:

                    run_restriction_cloning_ui(
                        dna=dna,
                        header=header,
                        sequences=sequences,
                        pipeline=pipeline,
                        topology=topology,
                    )

                except Exception as error:

                    report_error("RESTRICTION & CLONING ERROR", error)


            # ====================================================
            # 6 - PRIMER DESIGN
            # ====================================================


            elif choice == "6":

                try:

                    run_primer_analysis_ui(
                        dna=dna,
                        header=header,
                        pipeline=pipeline,
                    )

                except Exception as error:

                    report_error("PRIMER DESIGN ERROR", error)


            # ====================================================
            # 7 - MUTATION & VARIANT ANALYSIS
            # ====================================================


            elif choice == "7":

                try:

                    run_mutation_variant_ui(
                        dna=dna,
                        header=header,
                        sequences=sequences,
                        pipeline=pipeline,
                    )

                except Exception as error:

                    report_error("MUTATION & VARIANT ANALYSIS ERROR", error)


            # ====================================================
            # 8 - CRISPR GUIDE DESIGN
            # ====================================================


            elif choice == "8":

                try:

                    crispr_result = crispr_menu(
                        sequence=dna,
                        sequence_name=header,
                        selected_orf=selected_orf,
                        pipeline=pipeline,
                    )

                    if crispr_result is not None:
                        pipeline.store(
                            "crispr",
                            crispr_result,
                            parameters={
                                "system": crispr_result.system.name,
                                "region": (
                                    f"{crispr_result.region.start}-"
                                    f"{crispr_result.region.end}"
                                ),
                                "feature_type": crispr_result.region.feature_type,
                            },
                            limitations=[
                                "Guide ranking is a sequence-level heuristic, "
                                "not an editing-efficiency prediction.",
                                "Off-target screening covers only supplied "
                                "sequences, not a genome.",
                            ],
                        )

                except Exception as error:

                    report_error("CRISPR GUIDE DESIGN ERROR", error)


            # ====================================================
            # 9 - BLAST SEARCH
            # ====================================================


            elif choice == "9":

                try:

                    store_blast_run(
                        pipeline,
                        blast_menu(
                            sequence=dna,
                            sequence_name=header,
                            selected_orf=selected_orf,
                            sequence_type="DNA",
                        ),
                    )

                except Exception as error:

                    report_error("BLAST SEARCH ERROR", error)


            # ====================================================
            # 10 - LITERATURE SEARCH
            # ====================================================


            elif choice == "10":

                try:

                    literature_main(build_sequence_literature_context(
                        pipeline.context.header,
                        pipeline.context.selected_orf if isinstance(pipeline.context.selected_orf, dict) else None))

                except Exception as error:

                    report_error("LITERATURE SEARCH ERROR", error)


            # ====================================================
            # 11 - AI EXPERIMENT PLANNER
            # ====================================================


            elif choice == "11":

                show_planner_status()


            # ====================================================
            # P - PIPELINE STATUS
            # ====================================================


            elif choice == "R":
                export_report_interactive(pipeline)

            elif choice == "P":

                pipeline.print_status(print_pipeline_status)

                style.pause()


            # ====================================================
            # S - SELECT ANOTHER SEQUENCE
            # ====================================================


            elif choice == "S":

                print()
                style.note(
                    "Returning to "
                    "sequence selection..."
                )

                break


            # ====================================================
            # 0 - EXIT
            # ====================================================


            elif choice == "0":

                print(
                    "\nThank you for using Synapse."
                )

                raise SystemExit


            # ====================================================
            # INVALID INPUT
            # ====================================================


            else:

                print()
                style.invalid_selection()

                style.pause()


def cli():
    """Entry point of the ``synapse`` command and ``python main.py`` (Phase 9).

    Ends the session cleanly when the input ends (Ctrl+Z / Ctrl+D, or the
    end of piped input) or on Ctrl+C, instead of printing a Python traceback.
    Returns the process exit code.
    """
    try:
        main()
    except (EOFError, KeyboardInterrupt) as reason:
        print()
        style.note("Synapse closed (input ended)." if isinstance(reason, EOFError) else "Synapse closed.")
        return 0 if isinstance(reason, EOFError) else 130
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
