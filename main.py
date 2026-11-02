import logging
import os
from datetime import datetime

from core.pipeline import SynapsePipeline
from core.models import SequenceRecord
from core.report import export_analysis_report
from core.logging_setup import configure_logging

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

from analyses.rna import (
    run_rna_analysis_ui,
)

from analyses.protein import (
    run_protein_analysis_from_orf,
    run_protein_analysis_from_record,
)

from analyses.orf import (
    run_orf_analysis_ui,
    DEFAULT_MIN_ORF_LENGTH,
)

from analyses.restriction import (
    run_restriction_cloning_ui,
)

from analyses.mutation import (
    run_mutation_variant_ui,
)

from analyses.primer import (
    run_primer_analysis_ui,
)

from analyses.crispr import (
    crispr_menu,
)

from analyses.blast import (
    blast_menu,
)

from analyses.literature import (
    main as literature_main,
    build_sequence_literature_context,
)

from ui.report import print_results
from ui.menu import select_gene, show_menu
from ui.style import major_rule, minor_rule


LOG_PATH = configure_logging()
logger = logging.getLogger("synapse.main")


def report_error(title, error):
    """Show an actionable error and keep the full traceback in the log file."""
    logger.exception("%s", title)
    print()
    print(major_rule())
    print(title)
    print(major_rule())
    print()
    print(f"{type(error).__name__}: {error}")
    print()
    print(f"Details (traceback) were written to: {LOG_PATH}")
    input("\nPress Enter to continue...")


# ============================================================
# MAIN PROGRAM
# ============================================================


current_folder = os.path.dirname(os.path.abspath(__file__))

default_sequence_file = os.path.join(
    current_folder,
    "examples",
    "dna.fasta",
)


# ============================================================
# LOAD SEQUENCE FILE (FASTA or GenBank)
# ============================================================


def choose_sequence_file():
    """Ask for a FASTA/GenBank file; Enter loads the bundled example."""

    while True:

        print()
        print(major_rule())
        print("SYNAPSE - LOAD SEQUENCE FILE")
        print(major_rule())
        print("FASTA (.fasta/.fa) or GenBank (.gb/.gbk).")
        print("FASTA headers may contain topology=circular.")
        print(f"Press Enter for the example: {default_sequence_file}")
        print("Type Q to quit.")

        raw = input("\nFile path: ").strip().strip('"')

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
            print()
            print(major_rule())
            print("SEQUENCE FILE ERROR")
            print(major_rule())
            print(f"\n{error}")


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


def export_report_interactive(pipeline):
    """Ask for a path and export every stored result."""
    default_report = os.path.join(
        current_folder,
        "exports",
        "synapse_report_"
        + datetime.now().strftime("%Y%m%d_%H%M%S")
        + ".json",
    )
    print("\nReport of all results stored in this session.")
    print("Use .json (complete) or .md (summary).")
    raw_path = input(
        f"Path [Enter = {default_report}]: "
    ).strip().strip('"')
    try:
        saved = export_analysis_report(
            pipeline, raw_path or default_report
        )
        print(f"Report written: {saved.resolve()}")
    except (OSError, ValueError) as error:
        report_error("REPORT EXPORT ERROR", error)
    input(
        "\nPress Enter to continue..."
    )


# ============================================================
# SEQUENCE TYPE (audit D14)
# ============================================================


def print_sequence_type_error(record):
    """Explain why a record cannot be analyzed."""
    result = classify_sequence(record.sequence, record.declared_type)
    reason = result.reasons[0]
    print()
    print(major_rule())
    print("SEQUENCE VALIDATION ERROR")
    print(major_rule())
    print(f"\n{reason[0].upper()}{reason[1:]}.\n")
    for position, symbol in result.invalid_symbols[:20]:
        print(f"Invalid character '{symbol}' at position {position}")
    if len(result.invalid_symbols) > 20:
        print(f"... and {len(result.invalid_symbols) - 20} more")
    print(major_rule())
    input("\nPress Enter to select another sequence...")


def resolve_sequence_type(record):
    """Return the record with a determined type, or None to choose another.

    A sequence whose type cannot be told from its letters is not guessed:
    the user declares DNA, RNA or protein.
    """
    while record.type_evidence == AMBIGUOUS:
        print()
        print(major_rule())
        print("SEQUENCE TYPE")
        print(major_rule())
        print("\nThe type of this sequence cannot be determined from its letters:")
        print(f"{record.type_reasons[0]}.")
        print("\n[1] DNA")
        print("[2] RNA")
        print("[3] Protein")
        print("[0] Cancel (choose another sequence)")
        choice = input("\nSequence type: ").strip()
        if choice == "0":
            return None
        declared = {"1": "DNA", "2": "RNA", "3": "PROTEIN"}.get(choice)
        if declared is None:
            print("\nInvalid selection.")
            continue
        candidate = record.with_declared_type(declared)
        if candidate.type_evidence == INVALID:
            print(f"\nNot possible: {candidate.type_reasons[0]}.")
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
        print()
        print(major_rule())
        print("SYNAPSE - PROTEIN RECORD")
        print(major_rule())
        print(f"\n{record.description}")
        print(
            f"Length: {record.length} aa | "
            f"Type: protein ({record.type_evidence.lower()})"
        )
        print("\n[1] Protein Analysis")
        print("[2] BLAST Search")
        print("[R] Export Analysis Report (JSON / Markdown)")
        print("[S] Select Another Sequence")
        print("[0] Exit")
        choice = input("\nSelect an option: ").strip().upper()

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
            print("\nReturning to sequence selection...")
            return
        elif choice == "0":
            print("\nThank you for using Synapse!")
            raise SystemExit
        else:
            print("\nInvalid selection.")
            input("\nPress Enter to continue...")


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
            print(f"\nNote: {note}")

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
            print("\nNote: RNA sequence, analyzed in the DNA alphabet (U read as T).")

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
        # SEQUENCE ANALYSIS
        # ========================================================


        nucleotide_counts = (
            count_nucleotides(dna)
        )

        nucleotide_percentages = (
            calculate_nucleotide_percentages(
                nucleotide_counts
            )
        )

        gc_content = (
            calculate_gc_content(dna)
        )

        at_content = (
            calculate_at_content(dna)
        )

        content_ranges = (
            calculate_dna_content_ranges(
                dna
            )
        )

        gc_skew = (
            calculate_gc_skew(dna)
        )

        at_skew = (
            calculate_at_skew(dna)
        )

        sequence_complexity = (
            calculate_sequence_complexity(
                dna
            )
        )

        homopolymers = (
            find_homopolymers(dna)
        )

        cpg_analysis = (
            analyze_cpg(dna)
        )

        cpg_islands = (
            find_cpg_islands(dna)
        )

        molecular_weight = (
            calculate_dna_molecular_weight(
                dna
            )
        )

        sequence_quality = (
            analyze_sequence_quality(
                dna
            )
        )

        ambiguous_counts = (
            count_ambiguous_bases(
                dna
            )
        )

        complement = (
            get_complement(dna)
        )

        reverse_complement = (
            get_reverse_complement(dna)
        )


        # ========================================================
        # SEQUENCE OVERVIEW
        # ========================================================


        print_results(
            selected_index,
            header,
            dna,
            topology,
            nucleotide_counts,
            nucleotide_percentages,
            gc_content,
            at_content,
            content_ranges,
            gc_skew,
            at_skew,
            sequence_complexity,
            homopolymers,
            cpg_analysis,
            cpg_islands,
            molecular_weight,
            sequence_quality,
            ambiguous_counts,
            complement,
            reverse_complement,
        )

        pipeline.context.store_sequence_analysis({
            "nucleotide_counts": nucleotide_counts,
            "nucleotide_percentages": nucleotide_percentages,
            "gc_content": gc_content,
            "at_content": at_content,
            "content_ranges": content_ranges,
            "gc_skew": gc_skew,
            "at_skew": at_skew,
            "sequence_complexity": sequence_complexity,
            "homopolymers": homopolymers,
            "cpg_analysis": cpg_analysis,
            "cpg_islands": cpg_islands,
            "molecular_weight": molecular_weight,
            "sequence_quality": sequence_quality,
            "ambiguous_counts": ambiguous_counts,
            "complement": complement,
            "reverse_complement": reverse_complement,
        })


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

                print_results(
                    selected_index,
                    header,
                    dna,
                    topology,
                    nucleotide_counts,
                    nucleotide_percentages,
                    gc_content,
                    at_content,
                    content_ranges,
                    gc_skew,
                    at_skew,
                    sequence_complexity,
                    homopolymers,
                    cpg_analysis,
                    cpg_islands,
                    molecular_weight,
                    sequence_quality,
                    ambiguous_counts,
                    complement,
                    reverse_complement,
                )

                input(
                    "\nPress Enter to continue..."
                )


            # ====================================================
            # 2 - ORF ANALYSIS
            # ====================================================


            elif choice == "2":

                try:

                    selected_orf_result = run_orf_analysis_ui(
                        dna,
                        include_partial=True,
                        min_length=DEFAULT_MIN_ORF_LENGTH,
                    )

                    if selected_orf_result is not None:

                        selected_orf = selected_orf_result

                        pipeline.set_selected_orf(
                            selected_orf
                        )

                        pipeline.store(
                            "orf",
                            selected_orf,
                            parameters=selected_orf.get(
                                "search_parameters", {}
                            ),
                        )

                        print()
                        print(major_rule())
                        print("ORF CONNECTED TO PIPELINE")
                        print(major_rule())
                        print()
                        print(
                            f"Selected: ORF-{selected_orf['rank']:03d}"
                        )
                        print(
                            f"Length:   {selected_orf['length']} bp"
                        )
                        print(
                            f"Protein:  {selected_orf['protein_length']} aa"
                        )
                        print(
                            f"Location: {selected_orf['start'] + 1}-"
                            f"{selected_orf['stop'] + 1} "
                            f"({selected_orf['strand']} strand, 1-based)"
                        )
                        print()
                        print("✓ ORF stored in pipeline")
                        print("✓ RNA can use this ORF")
                        print("✓ Protein can use this ORF")
                        print("✓ CRISPR Guide Design can use this ORF")
                        if selected_orf.get("genetic_code", 1) != 1:
                            print(
                                "! ORF found with NCBI table "
                                f"{selected_orf['genetic_code']}; Protein analysis uses this "
                                "table, the RNA module's codon view uses the standard code."
                            )

                        input(
                            "\nPress Enter to continue..."
                        )

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

                    print()
                    print(major_rule())
                    print("PROTEIN ANALYSIS")
                    print(major_rule())
                    print()
                    print("No ORF has been selected.")
                    print()
                    print(
                        "Please run [2] ORF Analysis first and "
                        "select an ORF."
                    )
                    print()
                    print(
                        "Protein Analysis will not use the full DNA "
                        "sequence as a fallback."
                    )

                    input(
                        "\nPress Enter to continue..."
                    )

                    continue

                try:

                    orf_sequence = selected_orf.get(
                        "sequence",
                        ""
                    )

                    orf_rank = selected_orf.get(
                        "rank"
                    )

                    if orf_rank is not None:
                        orf_id = f"ORF-{orf_rank:03d}"
                    else:
                        orf_id = "Selected ORF"

                    strand = selected_orf.get(
                        "strand"
                    )

                    reading_frame = selected_orf.get(
                        "reading_frame",
                        selected_orf.get(
                            "frame"
                        )
                    )

                    orf_status = selected_orf.get(
                        "status"
                    )

                    run_protein_analysis_from_orf(
                        orf_sequence=orf_sequence,
                        orf_id=orf_id,
                        strand=strand,
                        reading_frame=reading_frame,
                        orf_status=orf_status,
                        pipeline=pipeline,
                        genetic_code=selected_orf.get("genetic_code", 1),
                    )

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

                # Audit D15: advertised but not implemented (P3-1, deferred).
                print(
                    "\nAI Experiment Planner: not implemented. "
                    "Nothing is run or stored."
                )
                print(
                    "It is planned as an evidence-bound planner (P3-1, "
                    "docs/AUDIT_STATUS.md) and is deferred."
                )

                input(
                    "\nPress Enter to continue..."
                )


            # ====================================================
            # P - PIPELINE STATUS
            # ====================================================


            elif choice == "R":
                export_report_interactive(pipeline)

            elif choice == "P":

                pipeline.print_status()

                input(
                    "\nPress Enter to continue..."
                )


            # ====================================================
            # S - SELECT ANOTHER SEQUENCE
            # ====================================================


            elif choice == "S":

                print(
                    "\nReturning to "
                    "sequence selection..."
                )

                break


            # ====================================================
            # 0 - EXIT
            # ====================================================


            elif choice == "0":

                print(
                    "\nThank you for using Synapse!"
                )

                raise SystemExit


            # ====================================================
            # INVALID INPUT
            # ====================================================


            else:

                print(
                    "\nInvalid selection."
                )

                input(
                    "\nPress Enter to continue..."
                )


if __name__ == "__main__":
    main()
