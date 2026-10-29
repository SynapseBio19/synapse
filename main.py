import os

from core.pipeline import SynapsePipeline

from fileio.fasta import load_fasta

from core.validation import (
    validate_dna,
    normalize_dna,
    get_dna_errors,
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
)

from ui.report import print_results
from ui.menu import select_gene, show_menu


# ============================================================
# MAIN PROGRAM
# ============================================================


current_folder = os.path.dirname(__file__)

fasta_file = os.path.join(
    current_folder,
    "examples",
    "dna.fasta",
)


# ============================================================
# LOAD FASTA FILE
# ============================================================


try:

    sequences = load_fasta(
        fasta_file
    )

except ValueError as error:

    print()
    print("=" * 60)
    print("FASTA FILE ERROR")
    print("=" * 60)

    print(f"\n{error}")

    print()
    print("=" * 60)

    input(
        "\nPress Enter to exit..."
    )

    raise SystemExit


# ============================================================
# SEQUENCE SELECTION LOOP
# ============================================================


while True:

    selected_index, record = (
        select_gene(sequences)
    )

    header = record["header"]
    dna = record["sequence"]

    topology = record.get(
        "topology",
        "linear",
    )

    dna = normalize_dna(dna)

    # ========================================================
    # SYNAPSE PIPELINE
    # ========================================================

    pipeline = SynapsePipeline(
        sequence=dna,
        header=header,
        sequence_id=selected_index,
        topology=topology,
    )

    # Selected ORF is shared by ORF -> RNA -> Protein.
    selected_orf = None


    # ========================================================
    # DNA VALIDATION
    # ========================================================


    if not validate_dna(dna):

        errors = get_dna_errors(dna)

        print()
        print("=" * 60)
        print("SEQUENCE VALIDATION ERROR")
        print("=" * 60)

        print(
            "\nInvalid DNA sequence.\n"
        )

        for error in errors:

            print(
                f"Invalid character "
                f"'{error['character']}' "
                f"at position "
                f"{error['position']}"
            )

        print("=" * 60)

        input(
            "\nPress Enter to select "
            "another sequence..."
        )

        continue


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

                    print()
                    print("=" * 68)
                    print("ORF CONNECTED TO PIPELINE")
                    print("=" * 68)
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
                    print()
                    print("✓ ORF stored in pipeline")
                    print("✓ RNA can use this ORF")
                    print("✓ Protein can use this ORF")
                    print("✓ CRISPR Guide Design can use this ORF")

                    input(
                        "\nPress Enter to continue..."
                    )

            except Exception as error:

                print()
                print("=" * 68)
                print("ORF ANALYSIS ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


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

                if rna_analysis is not None:

                    print()
                    print("=" * 68)
                    print("RNA CONNECTED TO PIPELINE")
                    print("=" * 68)
                    print()

                    source = (
                        rna_analysis
                        .get("pipeline_input", {})
                        .get("source", "Unknown")
                    )

                    print(
                        f"Source: {source}"
                    )

                    print(
                        f"RNA Length: "
                        f"{len(rna_analysis.get('rna', ''))} nt"
                    )

                    print(
                        "✓ RNA analysis stored in pipeline"
                    )

                    if selected_orf is not None:
                        print(
                            "✓ Selected ORF was used as RNA input"
                        )
                    else:
                        print(
                            "✓ Full DNA sequence was used as RNA input"
                        )

                    input(
                        "\nPress Enter to continue..."
                    )

            except Exception as error:

                print()
                print("=" * 68)
                print("RNA ANALYSIS ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


        # ====================================================
        # 4 - PROTEIN ANALYSIS
        # ====================================================


        elif choice == "4":

            if selected_orf is None:

                print()
                print("=" * 68)
                print("PROTEIN ANALYSIS")
                print("=" * 68)
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
                )

            except Exception as error:

                print()
                print("=" * 68)
                print("PROTEIN ANALYSIS ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


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

                print()
                print("=" * 68)
                print("RESTRICTION & CLONING ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


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

                print()
                print("=" * 68)
                print("PRIMER DESIGN ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


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

                print()
                print("=" * 68)
                print("MUTATION & VARIANT ANALYSIS ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


        # ====================================================
        # 8 - CRISPR GUIDE DESIGN
        # ====================================================


        elif choice == "8":

            try:

                crispr_menu(
                    sequence=dna,
                    sequence_name=header,
                    selected_orf=selected_orf,
                )

            except Exception as error:

                print()
                print("=" * 68)
                print("CRISPR GUIDE DESIGN ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


        # ====================================================
        # 9 - BLAST SEARCH
        # ====================================================


        elif choice == "9":

            try:

                blast_menu(
                    sequence=dna,
                    sequence_name=header,
                    selected_orf=selected_orf,
                )

            except Exception as error:

                print()
                print("=" * 68)
                print("BLAST SEARCH ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


        # ====================================================
        # 10 - LITERATURE SEARCH
        # ====================================================


        elif choice == "10":

            try:

                literature_main()

            except Exception as error:

                print()
                print("=" * 68)
                print("LITERATURE SEARCH ERROR")
                print("=" * 68)
                print()
                print(
                    f"{type(error).__name__}: {error}"
                )

                input(
                    "\nPress Enter to continue..."
                )


        # ====================================================
        # 11 - AI EXPERIMENT PLANNER
        # ====================================================


        elif choice == "11":

            print(
                "\nAI Experiment Planner "
                "is coming soon."
            )

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
