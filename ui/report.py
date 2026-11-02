from analyses.sequence import (
    calculate_gc_content_profile,
    calculate_complexity_profile,
    calculate_profile_statistics,
    calculate_average_dna_molecular_weight,
)
from ui.style import major_rule, minor_rule
from ui import style


# ============================================================
# REPORT HELPERS
# ============================================================


# ============================================================
# REPORT DATA HELPERS
# ============================================================

def _sequence_statistics(
    dna,
    nucleotide_counts,
    nucleotide_percentages,
):
    gc_count = (
        nucleotide_counts["G"]
        + nucleotide_counts["C"]
    )

    at_count = (
        nucleotide_counts["A"]
        + nucleotide_counts["T"]
    )

    standard_count = sum(
        nucleotide_counts.values()
    )

    return {
        "length": len(dna),
        "standard_length": standard_count,
        "gc_count": gc_count,
        "at_count": at_count,
        "gc_content": (
            gc_count / standard_count * 100
            if standard_count
            else 0.0
        ),
        "at_content": (
            at_count / standard_count * 100
            if standard_count
            else 0.0
        ),
        "nucleotide_counts": nucleotide_counts,
        "nucleotide_percentages": nucleotide_percentages,
    }


def _print_profile_unavailable(
    profile_name,
    sequence_length,
    window_size,
):
    style.note(f"{profile_name} is not available for this sequence.")
    style.field('Sequence length', f"{sequence_length} bp")
    style.field('Minimum window', f"{window_size} bp")
    style.field('Reason', "Sequence is shorter than the configured window.")


def _print_cpg_criteria():
    style.field('Criteria', "200 bp window | >=50% GC | O/E >=0.60")
    style.field('Definition', "Classical CpG island screening criteria")



def _preview(sequence, max_length=80):

    if len(sequence) <= max_length:
        return sequence

    half = max_length // 2

    return (
        sequence[:half]
        + "..."
        + sequence[-half:]
    )


def _print_profile_statistics(
    statistics,
    value_label,
    unit=""
):

    if statistics["count"] == 0:

        print(
            f"No {value_label.lower()} "
            "profile values available."
        )

        return

    style.field('Average', f"{statistics['average']:.2f}{unit}")

    style.field('Minimum', f"{statistics['minimum']:.2f}{unit}")

    style.field('Maximum', f"{statistics['maximum']:.2f}{unit}")

    minimum_position = (
        statistics["minimum_position"]
    )

    maximum_position = (
        statistics["maximum_position"]
    )

    if minimum_position is not None:

        style.field('Minimum position', f"{minimum_position['start']}-"
            f"{minimum_position['end']} bp")

    if maximum_position is not None:

        style.field('Maximum position', f"{maximum_position['start']}-"
            f"{maximum_position['end']} bp")


# ============================================================
# SEQUENCE REPORT
# ============================================================

def print_results(
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
    reverse_complement
):

    # ========================================================
    # NEW PROFILES
    # ========================================================

    profile_window = 100
    profile_step = 1

    gc_profile = calculate_gc_content_profile(
        dna,
        window_size=profile_window,
        step=profile_step
    )

    complexity_profile = calculate_complexity_profile(
        dna,
        window_size=profile_window,
        step=profile_step
    )

    gc_profile_statistics = (
        calculate_profile_statistics(
            gc_profile,
            "gc_content"
        )
    )

    complexity_profile_statistics = (
        calculate_profile_statistics(
            complexity_profile,
            "complexity"
        )
    )

    # ========================================================
    # HEADER
    # ========================================================

    style.header('SEQUENCE OVERVIEW')

    # ========================================================
    # SEQUENCE INFORMATION
    # ========================================================

    style.section('Sequence Information')

    style.field('Sequence', f"{selected_index}")

    style.field('Header', f"{header}")

    style.field('Topology', f"{topology}")

    style.field('Sequence Length', f"{len(dna)} bp")

    # ========================================================
    # SEQUENCE STATISTICS
    # ========================================================

    statistics = _sequence_statistics(
        dna,
        nucleotide_counts,
        nucleotide_percentages
    )

    style.section('Sequence Statistics')

    style.field('Length', f"{statistics['length']} bp")

    style.field('Standard Bases', f"{statistics['standard_length']}")

    style.field('GC Count', f"{statistics['gc_count']}")

    style.field('AT Count', f"{statistics['at_count']}")

    style.field('GC Content', f"{statistics['gc_content']:.2f}%")

    style.field('AT Content', f"{statistics['at_content']:.2f}%")

    # ========================================================
    # SEQUENCE SUMMARY
    # ========================================================

    style.section('Sequence Summary')

    if sequence_quality["invalid_symbols"] > 0:

        sequence_status = "INVALID"

    elif sequence_quality["ambiguous_bases"] > 0:

        sequence_status = "VALID - AMBIGUOUS"

    else:

        sequence_status = "VALID"

    style.field('Status', f"{sequence_status}")

    style.field('GC Content', f"{gc_content:.2f}%")

    style.field('AT Content', f"{at_content:.2f}%")

    style.field('Complexity', f"{sequence_complexity:.3f}")

    style.field('CpG Sites', f"{cpg_analysis['count']}")

    style.field('CpG Islands', f"{len(cpg_islands)}")

    style.field('Homopolymer Regions', f"{len(homopolymers)}")

    style.field('Ambiguous Bases', f"{sequence_quality['ambiguous_bases']}")

    style.field('Invalid Symbols', f"{sequence_quality['invalid_symbols']}")

    # ========================================================
    # GC / AT CONTENT
    # ========================================================

    style.section('GC / AT Content')

    if content_ranges["ambiguous"]:

        style.field('GC Content Range', f"{content_ranges['gc_min']:.2f}% - "
            f"{content_ranges['gc_max']:.2f}%")

        style.field('AT Content Range', f"{content_ranges['at_min']:.2f}% - "
            f"{content_ranges['at_max']:.2f}%")

        style.field('Content Status', "Range due to ambiguous IUPAC bases")

    else:

        style.field('GC Content', f"{gc_content:.2f}%")

        style.field('AT Content', f"{at_content:.2f}%")

    # ========================================================
    # NUCLEOTIDE COMPOSITION
    # ========================================================

    style.section('Nucleotide Composition')

    for base in (
        "A",
        "T",
        "C",
        "G"
    ):

        style.field(base, f"{nucleotide_counts[base]:>6}  "
                    f"{nucleotide_percentages[base]:>6.2f}%")

    # ========================================================
    # NUCLEOTIDE SKEW
    # ========================================================

    style.section('Nucleotide Skew')

    style.field('GC Skew', f"{gc_skew:+.3f}")

    style.field('AT Skew', f"{at_skew:+.3f}")

    style.field('Interpretation', "Positive values indicate relative enrichment "
        "of the first nucleotide in the pair.")

    # ========================================================
    # MOLECULAR WEIGHT
    # ========================================================

    style.section('Molecular Weight')

    if molecular_weight["ambiguous"]:

        style.field('Estimated dsDNA MW', f"{molecular_weight['min'] / 1000:,.2f} - "
            f"{molecular_weight['max'] / 1000:,.2f} kDa")

        style.field('MW Status', "Range due to ambiguous IUPAC bases")

    else:

        style.field('Estimated dsDNA MW', f"{molecular_weight['min'] / 1000:,.2f} kDa")

    style.field('Method', "Composition-based nucleotide residue estimate")

    # ========================================================
    # GC CONTENT PROFILE
    # ========================================================

    style.section('GC Content Profile')

    if len(dna) < profile_window:

        _print_profile_unavailable(
            "GC Content Profile",
            len(dna),
            profile_window
        )

    else:

        style.field('Window', f"{profile_window} bp")

        style.field('Step', f"{profile_step} bp")

        _print_profile_statistics(
            gc_profile_statistics,
            "GC Content",
            "%"
        )

        style.section('Profile Preview')

        preview_points = (
            gc_profile[::max(
                1,
                len(gc_profile) // 10
            )][:10]
        )

        for item in preview_points:

            if item["gc_content"] is None:

                value = "N/A"

            else:

                value = (
                    f"{item['gc_content']:.2f}%"
                )

            print(
                f"{item['start_1_based']:>6}-"
                f"{item['end_1_based']:<6} "
                f"{value}"
            )

    # ========================================================
    # COMPLEXITY PROFILE
    # ========================================================

    style.section('Sequence Complexity Profile')

    if len(dna) < profile_window:

        _print_profile_unavailable(
            "Sequence Complexity Profile",
            len(dna),
            profile_window
        )

    else:

        style.field('Window', f"{profile_window} bp")

        style.field('Step', f"{profile_step} bp")

        _print_profile_statistics(
            complexity_profile_statistics,
            "Complexity"
        )

        style.section('Profile Preview')

        preview_points = (
            complexity_profile[::max(
                1,
                len(complexity_profile) // 10
            )][:10]
        )

        for item in preview_points:

            print(
                f"{item['start_1_based']:>6}-"
                f"{item['end_1_based']:<6} "
                f"{item['complexity']:.3f}"
            )

    # ========================================================
    # HOMOPOLYMERS
    # ========================================================

    style.section('Homopolymers')

    if homopolymers:

        for region in homopolymers[:20]:

            print(
                f"{region['base']} × "
                f"{region['length']} | "
                f"Position "
                f"{region['start_1_based']}-"
                f"{region['end_1_based']}"
            )

        if len(homopolymers) > 20:

            print(
                f"... and "
                f"{len(homopolymers) - 20} more"
            )

    else:

        print(
            "No homopolymer regions detected."
        )

    # ========================================================
    # CpG ANALYSIS
    # ========================================================

    style.section('CpG Analysis')

    style.field('CpG Count', f"{cpg_analysis['count']}")

    style.field('CpG Density', f"{cpg_analysis['density']:.3f}%")

    style.field('Observed/Expected', f"{cpg_analysis['oe_ratio']:.3f}")

    style.field('Ambiguous Contexts', f"{cpg_analysis['ambiguous_contexts']}")

    # ========================================================
    # CpG ISLANDS
    # ========================================================

    style.section('CpG Islands')

    _print_cpg_criteria()

    if cpg_islands:

        for index, island in enumerate(
            cpg_islands[:20],
            start=1
        ):

            print(
                f"Island {index}:       "
                f"{island['start_1_based']}-"
                f"{island['end_1_based']} | "
                f"{island['length']} bp | "
                f"GC "
                f"{island['gc_content']:.2f}% | "
                f"O/E "
                f"{island['oe_ratio']:.2f}"
            )

        if len(cpg_islands) > 20:

            print(
                f"... and "
                f"{len(cpg_islands) - 20} more"
            )

    else:

        print(
            "No CpG islands detected "
            "under the configured criteria."
        )

    # ========================================================
    # SEQUENCE QUALITY
    # ========================================================

    style.section('Sequence Quality')

    style.field('Standard Bases', f"{sequence_quality['standard_bases']}")

    style.field('Ambiguous Bases', f"{sequence_quality['ambiguous_bases']}")

    style.field('Invalid Symbols', f"{sequence_quality['invalid_symbols']}")

    style.field('Standard', f"{sequence_quality['standard_percentage']:.2f}%")

    style.field('Ambiguous', f"{sequence_quality['ambiguous_percentage']:.2f}%")

    style.field('Invalid', f"{sequence_quality['invalid_percentage']:.2f}%")

    if sequence_quality["invalid_positions"]:

        print()
        style.field(
            "Invalid positions",
            ", ".join(
                str(item["position_1_based"])
                for item in
                sequence_quality[
                    "invalid_positions"
                ][:20]
            )
        )

    if sequence_quality["ambiguous_positions"]:

        style.field(
            "Ambiguous positions",
            ", ".join(
                str(item["position_1_based"])
                for item in
                sequence_quality[
                    "ambiguous_positions"
                ][:20]
            )
        )

    # ========================================================
    # IUPAC ANALYSIS
    # ========================================================

    style.section('IUPAC / Ambiguous Bases')

    if ambiguous_counts:

        for code, count in (
            ambiguous_counts.items()
        ):

            print(
                f"{code}:                  "
                f"{count:>6} | "
                f"{code}"
            )

    else:

        print(
            "No ambiguous IUPAC bases detected."
        )

    # ========================================================
    # COMPLEMENTARY SEQUENCES
    # ========================================================

    style.section('Complementary Sequences')

    style.field('Orientation', "complement 3' -> 5'; reverse complement 5' -> 3'")

    style.field('Complement', _preview(complement))

    style.field('Reverse complement', _preview(reverse_complement))

    # ========================================================
    # SEQUENCE PREVIEW
    # ========================================================

    style.section('Sequence Preview')

    style.field("Sequence (5' -> 3')", _preview(dna))

    style.field('Length shown', f"{min(len(dna), 80)} / "
        f"{len(dna)} bp")

    # ========================================================
    # INTERPRETATION
    # ========================================================

    style.section('Interpretation')

    if sequence_quality["invalid_symbols"] > 0:

        style.field('Sequence status', "INVALID")

        style.warning(
            "The sequence contains invalid "
            "symbols and should be reviewed "
            "before downstream analysis."
        )

    elif sequence_quality["ambiguous_bases"] > 0:

        style.field('Sequence status', "VALID - AMBIGUOUS")

        style.note(
            "The sequence is valid but contains "
            "IUPAC ambiguity codes. Some values "
            "are therefore represented as ranges."
        )

    else:

        style.field('Sequence status', "VALID")

        print(
            "The sequence contains only standard "
            "DNA bases and passed sequence-level "
            "quality screening."
        )

    if sequence_complexity >= 0.75:

        style.field('Complexity class', "HIGH")

    elif sequence_complexity >= 0.50:

        style.field('Complexity class', "MODERATE")

    else:

        style.field('Complexity class', "LOW")

    if gc_content >= 60:

        style.field('GC composition', "Relatively high")

    elif gc_content <= 40:

        style.field('GC composition', "Relatively low")

    else:

        style.field('GC composition', "Relatively balanced")

    style.field('CpG enrichment', f"O/E {cpg_analysis['oe_ratio']:.3f}")

    style.field('CpG island screen', f"{'Detected' if cpg_islands else 'Not detected'}")

    style.field('Homopolymers', f"{'Detected' if homopolymers else 'None detected'}")

    print()
    style.note(
        "Profile analyses are sequence-level "
        "descriptors and do not establish biological "
        "function or experimental behavior."
    )

    print()
    print(major_rule())