"""Terminal dialogs of analyses/rna.py (moved here in P2-2).

The scientific functions stay in ``analyses.rna``; this module draws the
screens, reads the answers and calls them.
"""

from analyses.rna import (  # noqa: F401
    RNA_AMBIGUOUS_CODES,
    _orf_value,
    analyze_rna,
    analyze_rna_from_pipeline,
    analyze_rna_three_frames,
    create_rna_preview,
)
from ui import style
from ui.style import major_rule, minor_rule


def print_rna_translation_report(rna, genetic_code=1):
    """Print the three forward RNA reading-frame previews."""

    translation = analyze_rna_three_frames(rna, genetic_code=genetic_code)

    style.header('RNA TRANSLATION & CODON ANALYSIS')

    for frame in translation["frames"]:

        print()
        print(
            f"Frame +{frame['frame']}"
        )

        print(minor_rule())

        style.field('Codons', f"{frame['codon_count']}")

        style.field('Protein length', f"{frame['protein_length']} aa")

        style.field('Stop detected', f"{'YES' if frame['complete'] else 'NO'}")

        style.field('Stop position', f"{frame['stop_position'] or '-'}")

        style.field('Protein preview', f"{frame['protein'][:60] or '-'}")

    best = translation["best_frame"]

    if best is not None:

        print()
        print(
            f"Longest/complete frame by heuristic: "
            f"+{best['frame']}"
        )

    print()
    print(
        "Interpretation: frame comparison is sequence-based "
        "and does not establish biological translation."
    )

    print(major_rule())


def print_rna_feature_report(analysis):
    """Print feature-level RNA screening results."""

    score = analysis["feature_score"]

    style.header('RNA FEATURES')

    style.section('Sequence Feature Screening')

    style.field('Screening score', f"{score['score']}/100")

    style.field('Screening band', f"{score['classification']}")

    style.section('Sequence Features')

    style.field('Complexity', f"{analysis['complexity']:.3f}")

    style.field('Homopolymers', f"{len(analysis['homopolymers'])}")

    style.field('Hairpin candidates', f"{analysis['structure']['hairpin_count']}")

    style.field('Poly-A regions', f"{len(analysis['poly_a'])}")

    style.field('Poly-U regions', f"{len(analysis['poly_u'])}")

    style.section('Structure Screening')
    style.field('Method', "Heuristic stem-loop screening")
    style.field('Thermodynamics', "Not calculated")

    style.section('Start / Stop Landscape')

    style.field('AUG starts', f"{analysis['start_stop']['start_count']}")

    style.field('Stop codons', f"{analysis['start_stop']['stop_count']}")

    style.section('Interpretation')

    for reason in score["reasons"]:
        print(f"- {reason}")

    print()
    print(
        "The screening score describes sequence-level feature detection. "
        "It is not a measure of biological quality, expression, "
        "stability or biological function."
    )

    print(major_rule())


def print_rna_report(
    dna,
    analysis=None,
    source_name=None
):

    if analysis is None:
        analysis = analyze_rna(
            dna
        )

    rna = analysis["rna"]

    counts = analysis["base_counts"]

    percentages = (
        analysis["base_percentages"]
    )

    quality = analysis["quality"]

    ambiguity = analysis["ambiguity"]

    ranges = analysis["content_ranges"]

    molecular_weight = (
        analysis["molecular_weight"]
    )

    pipeline_info = analysis.get(
        "pipeline_input",
        {}
    )

    source_name = (
        source_name
        or pipeline_info.get(
            "source",
            "Direct DNA Input"
        )
    )

    style.header('COMPLETE RNA REPORT')

    # ========================================================
    # INPUT CONTEXT
    # ========================================================

    style.section('Input Context')

    style.field('Source', f"{source_name}")

    if (
        pipeline_info.get("mode") == "orf"
        and pipeline_info.get("selected_orf") is not None
    ):
        orf = pipeline_info["selected_orf"]
        rank = _orf_value(orf, "rank", None)

        try:
            orf_name = (
                f"ORF-{int(rank):03d}"
                if rank is not None
                else "Selected ORF"
            )
        except (TypeError, ValueError):
            orf_name = str(rank)

        style.field('ORF', f"{orf_name}")
        style.field('ORF length', f"{pipeline_info.get('dna_length', len(dna))} bp")
        style.field('Strand', f"{_orf_value(orf, 'strand', '-')}")
        style.field('Reading frame', f"{_orf_value(orf, 'frame', '-')}")
        style.field('ORF status', f"{_orf_value(orf, 'status', '-')}")
        style.field('Analysis scope', "Selected ORF")

    else:
        style.field('Sequence length', f"{pipeline_info.get('dna_length', len(dna))} bp")
        style.field('Analysis scope', "Entire DNA sequence")

    # ========================================================
    # TRANSCRIPTION
    # ========================================================

    style.section('Transcription Information')

    if pipeline_info.get("mode") == "orf":
        orf = pipeline_info.get("selected_orf")

        rank = _orf_value(
            orf,
            "rank",
            None
        )

        try:
            input_name = (
                f"ORF-{int(rank):03d}"
                if rank is not None
                else "Selected ORF"
            )
        except (TypeError, ValueError):
            input_name = str(rank)

        style.field('Input source', "Selected ORF")

        style.field('Input', f"{input_name}")

    else:
        style.field('Input source', "Full DNA Sequence")

        style.field('Input', "Complete DNA sequence")

    style.field('Interpretation', "Coding / sense strand")

    style.field('DNA Direction', "5' -> 3'")

    style.field('RNA Direction', "5' -> 3'")

    style.field('Transcription Rule', "T -> U")

    # ========================================================
    # BASIC INFORMATION
    # ========================================================

    style.section('RNA Sequence Information')

    style.field('RNA Length', f"{analysis['length']} nt")

    if ranges["ambiguous"]:

        style.field('GC Content Range', f"{ranges['gc_min']:.2f}% - "
            f"{ranges['gc_max']:.2f}%")

        style.field('AU Content Range', f"{ranges['au_min']:.2f}% - "
            f"{ranges['au_max']:.2f}%")

    else:

        style.field('GC Content', f"{analysis['gc_content']:.2f}%")

        style.field('AU Content', f"{analysis['au_content']:.2f}%")

    style.field('GC Skew', f"{analysis['gc_skew']:+.3f}")

    style.field('AU Skew', f"{analysis['au_skew']:+.3f}")

    # ========================================================
    # BASE COMPOSITION
    # ========================================================

    style.section('RNA Base Composition')

    for base in (
        "A",
        "U",
        "C",
        "G"
    ):

        print(
            f"{base} : "
            f"{counts[base]:>6}   "
            f"{percentages[base]:>6.2f}%"
        )

    # ========================================================
    # MOLECULAR WEIGHT
    # ========================================================

    style.section('Molecular Weight')

    if molecular_weight["ambiguous"]:

        style.field('Estimated ssRNA MW', f"{molecular_weight['min']:,.0f} - "
            f"{molecular_weight['max']:,.0f} Da")

    else:

        style.field('Estimated ssRNA MW', f"{molecular_weight['min']:,.0f} Da")

    style.field('Method', "Composition-based estimate")

    # ========================================================
    # COMPLEXITY
    # ========================================================

    style.section('RNA Sequence Complexity')

    style.field('Normalized Entropy', f"{analysis['complexity']:.3f}")

    # ========================================================
    # CODON STATISTICS
    # ========================================================

    codons = analysis["codon_usage"]

    style.section('Codon Statistics')

    style.field('Reading Frame', "Frame 0")

    style.field('Resolved Codons', f"{codons['total_codons']}")

    style.field('Unique Codons', f"{codons['unique_codons']}")

    style.field('Codon Diversity', f"{analysis['codon_diversity']:.3f}")

    if codons["gc3_content"] is not None:

        style.field('GC3 Content', f"{codons['gc3_content']:.2f}%")

    else:

        style.field('GC3 Content', "Not available")

    style.field('Interpretation', "Frame-dependent; coding relevance "
        "requires a valid CDS context")

    # ========================================================
    # INITIATION CONTEXT
    # ========================================================

    style.section('Translation Initiation Context')

    kozak = analysis["kozak"]

    if kozak["available"]:

        style.field('Kozak Context', f"{kozak['context']}")

        style.field('Kozak Score', f"{kozak['score']}/3")

        print(
            f"Kozak Classification:"
            f" {kozak['classification']}"
        )

    else:

        style.field('Kozak Context', "Unavailable")

    shine = analysis["shine_dalgarno"]

    if shine["detected"]:

        style.field('Shine-Dalgarno', f"{shine['motif']}")

        style.field('SD Distance', f"{shine['distance_to_start']} nt")

    else:

        style.field('Shine-Dalgarno', "Not detected")

    # ========================================================
    # HOMOPOLYMERS
    # ========================================================

    style.section('Homopolymer Regions')

    if analysis["homopolymers"]:

        for region in analysis["homopolymers"]:

            print(
                f"{region['base']} | "
                f"Length {region['length']} | "
                f"Position "
                f"{region['start']}-"
                f"{region['end']}"
            )

    else:

        print(
            "None detected "
            "(minimum length: 4)"
        )

    # ========================================================
    # POLY-A / POLY-U
    # ========================================================

    style.section('Poly-A / Poly-U Regions')

    style.field('Poly-A Regions', f"{len(analysis['poly_a'])}")

    style.field('Poly-U Regions', f"{len(analysis['poly_u'])}")

    # ========================================================
    # SECONDARY STRUCTURE
    # ========================================================

    structure = analysis["structure"]

    style.section('RNA Secondary Structure Screening')

    style.field('Potential Hairpins', f"{structure['hairpin_count']}")
    style.field('Including G-U pairs', f"{structure['hairpin_count_with_gu']}")

    if structure["hairpins"]:

        for index, hairpin in enumerate(
            structure["hairpins"][:5],
            start=1
        ):

            print(
                f"\nHairpin {index}"
            )

            style.field('Position', f"{hairpin['start']}-"
                f"{hairpin['end']}")

            style.field('Stem Length', f"{hairpin['stem_length']} bp")

            style.field('Loop Length', f"{hairpin['loop_length']} nt")

            style.field('Pairing Fraction', f"{hairpin['pairing_fraction']:.2f}")

    print()
    style.field('Method', "Heuristic stem-loop screening")

    style.field('Thermodynamics', "Not calculated")

    # ========================================================
    # QUALITY
    # ========================================================

    style.section('RNA Sequence Quality')

    style.field('Standard Bases', f"{quality['standard_count']}")

    style.field('Ambiguous Bases', f"{quality['ambiguous_count']}")

    style.field('Invalid Symbols', f"{quality['invalid_count']}")

    style.field('Ambiguous Content', f"{quality['ambiguous_percentage']:.2f}%")

    if not quality["valid"]:

        style.field('Status', "INVALID")

    elif quality["ambiguous_count"] > 0:

        style.field('Status', "VALID - IUPAC ambiguity present")

    else:

        style.field('Status', "VALID - Unambiguous RNA")

    # ========================================================
    # IUPAC
    # ========================================================

    if ambiguity["ambiguous_counts"]:

        style.section('IUPAC Ambiguous Base Composition')

        for code, count in (
            ambiguity["ambiguous_counts"].items()
        ):

            print(
                f"{code} : "
                f"{count:>6}   "
                f"({RNA_AMBIGUOUS_CODES[code]})"
            )

    # ========================================================
    # AMBIGUOUS POSITIONS
    # ========================================================

    if ambiguity["ambiguous_positions"]:

        style.section('Ambiguous Positions')

        for item in (
            ambiguity["ambiguous_positions"]
        ):

            possibilities = "/".join(
                item["possibilities"]
            )

            print(
                f"Position "
                f"{item['position_1_based']}: "
                f"{item['base']} -> "
                f"{possibilities}"
            )

    # ========================================================
    # INVALID POSITIONS
    # ========================================================

    if quality["invalid_positions"]:

        style.section('Invalid Symbols')

        for item in quality["invalid_positions"]:

            print(
                f"Position "
                f"{item['position_1_based']}: "
                f"{item['symbol']}"
            )

    # ========================================================
    # COMPLEMENT
    # ========================================================

    style.section('RNA Strand Transformations')

    style.field('Complement', f"{create_rna_preview(analysis['complement'])}")

    style.field('Reverse Complement', f"{create_rna_preview(analysis['reverse_complement'])}")

    # ========================================================
    # SEQUENCE PREVIEW
    # ========================================================

    style.section('Sequence Preview')

    print(
        f"DNA 5' -> 3':       "
        f"{create_rna_preview(dna)}"
    )

    print(
        f"RNA 5' -> 3':       "
        f"{create_rna_preview(rna)}"
    )

    # ========================================================
    # START / STOP LANDSCAPE
    # ========================================================

    landscape = analysis["start_stop"]

    style.section('Start / Stop Landscape')

    style.field('AUG Start Codons', f"{landscape['start_count']}")

    style.field('Stop Codons', f"{landscape['stop_count']}")

    if landscape["starts"]:

        print(
            "AUG Positions:       "
            + ", ".join(
                str(item["position_1_based"])
                for item in landscape["starts"][:20]
            )
        )

    if landscape["stops"]:

        print(
            "Stop Positions:      "
            + ", ".join(
                str(item["position_1_based"])
                for item in landscape["stops"][:20]
            )
        )

    # ========================================================
    # THREE-FRAME TRANSLATION PREVIEW
    # ========================================================

    translation_frames = (
        analysis["translation_frames"]
    )

    style.section('Three-Frame Translation Preview')

    for frame in translation_frames["frames"]:

        print(
            f"Frame +{frame['frame']}: "
            f"{frame['protein'][:60] or '-'}"
        )

    # ========================================================
    # FEATURE SCORE
    # ========================================================

    feature_score = (
        analysis["feature_score"]
    )

    style.section('Sequence Feature Screening')

    style.field('Screening Score', f"{feature_score['score']}/100")

    style.field('Screening Band', f"{feature_score['classification']}")

    # ========================================================
    # SUMMARY
    # ========================================================

    summary = analysis["summary"]

    style.section('RNA Summary')

    style.field('Length', f"{summary['length']} nt")

    style.field('GC Content', f"{summary['gc_content']:.2f}%")

    style.field('AU Content', f"{summary['au_content']:.2f}%")

    style.field('Ambiguous Bases', f"{summary['ambiguous_count']}")

    style.field('Invalid Symbols', f"{summary['invalid_count']}")

    style.field('AUG Starts', f"{summary['start_count']}")

    style.field('Stop Codons', f"{summary['stop_count']}")

    # ========================================================
    # SCIENTIFIC INTERPRETATION
    # ========================================================

    style.section('Scientific Interpretation')

    print(
        "The supplied DNA is interpreted as the coding "
        "(sense) strand in 5' -> 3' orientation."
    )

    print(
        "RNA is generated by replacing T with U while "
        "preserving IUPAC ambiguity."
    )

    print(
        "GC/AU content, skew and complexity are "
        "sequence-level descriptors."
    )

    print(
        "Codon statistics are frame-dependent and should "
        "only be interpreted as coding information when "
        "the RNA represents a CDS."
    )

    print(
        "Kozak and Shine-Dalgarno detection are "
        "sequence-based heuristics and do not prove "
        "translation initiation."
    )

    print(
        "Hairpin detection identifies candidate stem-loop "
        "motifs from Watson-Crick pairing (the G-U count also "
        "accepts wobble pairs) but does not calculate "
        "thermodynamic stability or minimum free energy."
    )

    print(
        "Poly-A/Poly-U regions are sequence motifs and "
        "do not by themselves prove RNA processing, "
        "termination or transcript stability."
    )

    print()

    print(major_rule())


def print_rna_pipeline_input(analysis):
    """Display the resolved input exactly once."""

    info = analysis.get(
        "pipeline_input",
        {}
    )

    style.section('Pipeline Input')

    style.field('Source', f"{info.get('source', 'Unknown')}")

    if info.get("mode") == "orf":
        style.field('ORF length', f"{info.get('dna_length', 0)} bp")
    else:
        style.field('Sequence length', f"{info.get('dna_length', 0)} bp")

    if (
        info.get("mode") == "orf"
        and info.get("selected_orf") is not None
    ):
        orf = info["selected_orf"]

        rank = _orf_value(
            orf,
            "rank",
            None
        )

        if rank is not None:
            try:
                orf_name = f"ORF-{int(rank):03d}"
            except (TypeError, ValueError):
                orf_name = str(rank)
        else:
            orf_name = "Selected ORF"

        style.field('ORF', f"{orf_name}")

        style.field('Strand', f"{_orf_value(orf, 'strand', '-')}")

        style.field('Reading frame', f"{_orf_value(orf, 'frame', '-')}")

        style.field('ORF status', f"{_orf_value(orf, 'status', '-')}")

        style.field('Analysis scope', "Selected ORF")

    else:
        style.field('Selected ORF', "None")

        style.field('Analysis scope', "Entire DNA sequence")


def print_rna_pipeline_connection(analysis):
    """Display the single pipeline connection confirmation."""

    info = analysis.get(
        "pipeline_input",
        {}
    )

    style.header('RNA CONNECTED TO PIPELINE')

    print()
    style.field('Source', f"{info.get('source', 'Unknown')}")

    style.field('RNA Length', f"{len(analysis.get('rna', ''))} nt")

    if info.get("mode") == "orf":
        style.success("Selected ORF was used as RNA input")
    else:
        style.success("Full DNA sequence was used as RNA input")

    style.success("RNA analysis stored in pipeline")


def print_rna_overview(analysis):
    """Compact professional overview for the RNA module."""

    summary = analysis["summary"]
    quality = analysis["quality"]
    score = analysis["feature_score"]
    info = analysis.get("pipeline_input", {})

    style.header('RNA OVERVIEW')

    style.section('Input Context')

    style.field('Source', f"{info.get('source', 'Unknown')}")

    if (
        info.get("mode") == "orf"
        and info.get("selected_orf") is not None
    ):
        orf = info["selected_orf"]

        rank = _orf_value(
            orf,
            "rank",
            None
        )

        try:
            orf_name = (
                f"ORF-{int(rank):03d}"
                if rank is not None
                else "Selected ORF"
            )
        except (TypeError, ValueError):
            orf_name = str(rank)

        style.field('ORF', f"{orf_name}")

        style.field('ORF length', f"{info.get('dna_length', summary['length'])} bp")

        style.field('Strand', f"{_orf_value(orf, 'strand', '-')}")

        style.field('Reading frame', f"{_orf_value(orf, 'frame', '-')}")

        style.field('ORF status', f"{_orf_value(orf, 'status', '-')}")

        style.field('Analysis scope', "Selected ORF")

    else:
        style.field('Sequence length', f"{info.get('dna_length', summary['length'])} bp")

        style.field('Selected ORF', "None")

        style.field('Analysis scope', "Entire DNA sequence")

    style.field('RNA length', f"{summary['length']} nt")

    style.section('Composition')
    style.field('GC content', f"{summary['gc_content']:.2f}%")

    style.field('AU content', f"{summary['au_content']:.2f}%")

    style.field('Complexity', f"{summary['complexity']:.3f}")

    style.section('Sequence Quality')
    style.field('Ambiguous bases', f"{summary['ambiguous_count']}")

    style.field('Invalid symbols', f"{summary['invalid_count']}")

    if quality["invalid_count"] > 0:
        status = "INVALID"

    elif quality["ambiguous_count"] > 0:
        status = "VALID - IUPAC ambiguity"

    else:
        status = "VALID - Unambiguous"

    style.field('Status', f"{status}")

    style.section('Coding Landscape')
    style.field('AUG starts', f"{summary['start_count']}")

    style.field('Stop codons', f"{summary['stop_count']}")

    style.section('Feature Screening')
    style.field('Feature score', f"{score['score']}/100")

    style.field('Classification', f"{score['classification']}")

    print()
    print(
        "This overview is sequence-level analysis. "
        "It does not establish biological function."
    )


def print_rna_sequence_view(analysis):
    """Display the complete resolved RNA sequence in a readable format."""

    rna = analysis.get(
        "rna",
        ""
    )

    style.header('RNA SEQUENCE')

    print()
    style.field('Length', f"{len(rna)} nt")

    print()
    print("5' -> 3'")
    print(minor_rule())

    if not rna:
        print("[empty]")
        return

    line_width = 60

    for start in range(
        0,
        len(rna),
        line_width
    ):
        end = min(
            start + line_width,
            len(rna)
        )

        print(
            f"{start + 1:>6}-{end:<6} "
            f"{rna[start:end]}"
        )

    print(minor_rule())


def _pause_rna_ui():
    style.pause()


def run_rna_analysis_ui(
    dna,
    pipeline=None,
    selected_orf=None,
    header=None,
    **analysis_kwargs
):
    """
    Main RNA module UI.

    This is the single public UI entry point used by main.py.

    Behavior:
        - selected ORF exists -> analyze selected ORF
        - no ORF exists -> analyze complete DNA
        - RNA analysis is performed once
        - pipeline connection is displayed once
        - user receives a dedicated RNA analysis menu

    'header' is retained for compatibility with older main.py code.
    """

    del header

    analysis = analyze_rna_from_pipeline(
        dna,
        pipeline=pipeline,
        selected_orf=selected_orf,
        **analysis_kwargs
    )

    style.header('RNA ANALYSIS')

    print_rna_pipeline_input(
        analysis
    )

    print_rna_pipeline_connection(
        analysis
    )

    _pause_rna_ui()

    while True:

        style.header('RNA ANALYSIS')

        style.menu([
            ("Overview", [('1', 'RNA Overview')]),
            ("Analysis", [('2', 'Translation & Codon Analysis'), ('3', 'RNA Features')]),
            ("Report", [('4', 'Complete RNA Report')]),
            ("Sequence", [('5', 'RNA Sequence')]),
        ])

        choice = style.ask_choice().lower()

        if choice == "1":

            print_rna_overview(
                analysis
            )

            _pause_rna_ui()

        elif choice == "2":

            print_rna_translation_report(
                analysis["rna"],
                genetic_code=analysis.get("genetic_code", 1),
            )

            _pause_rna_ui()

        elif choice == "3":

            print_rna_feature_report(
                analysis
            )

            _pause_rna_ui()

        elif choice == "4":

            print_rna_report(
                analysis.get("dna", ""),
                analysis=analysis,
                source_name=analysis.get("pipeline_input", {}).get("source")
            )

            _pause_rna_ui()

        elif choice == "5":

            print_rna_sequence_view(
                analysis
            )

            _pause_rna_ui()

        elif choice == "0":

            break

        else:

            style.invalid_selection()

    return analysis


run_rna_analysis_pipeline_ui = run_rna_analysis_ui


run_rna_ui_pipeline = run_rna_analysis_ui
