from analyses.sequence import (
    calculate_gc_content_profile,
    calculate_complexity_profile,
    calculate_profile_statistics,
    calculate_average_dna_molecular_weight,
)


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
    print(
        f"{profile_name} is not available for this sequence."
    )
    print(
        f"Sequence length:     {sequence_length} bp"
    )
    print(
        f"Minimum window:      {window_size} bp"
    )
    print(
        "Reason:               "
        "Sequence is shorter than the configured window."
    )


def _print_cpg_criteria():
    print(
        "Criteria:             "
        "200 bp window | >=50% GC | O/E >=0.60"
    )
    print(
        "Definition:           "
        "Classical CpG island screening criteria"
    )



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

    print(
        f"Average:             "
        f"{statistics['average']:.2f}{unit}"
    )

    print(
        f"Minimum:             "
        f"{statistics['minimum']:.2f}{unit}"
    )

    print(
        f"Maximum:             "
        f"{statistics['maximum']:.2f}{unit}"
    )

    minimum_position = (
        statistics["minimum_position"]
    )

    maximum_position = (
        statistics["maximum_position"]
    )

    if minimum_position is not None:

        print(
            f"Minimum position:    "
            f"{minimum_position['start']}-"
            f"{minimum_position['end']} bp"
        )

    if maximum_position is not None:

        print(
            f"Maximum position:    "
            f"{maximum_position['start']}-"
            f"{maximum_position['end']} bp"
        )


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

    print()
    print("=" * 70)
    print("SEQUENCE OVERVIEW")
    print("=" * 70)

    # ========================================================
    # SEQUENCE INFORMATION
    # ========================================================

    print()
    print("Sequence Information")
    print("-" * 45)

    print(
        f"Sequence:            "
        f"{selected_index}"
    )

    print(
        f"Header:              "
        f"{header}"
    )

    print(
        f"Topology:            "
        f"{topology}"
    )

    print(
        f"Sequence Length:     "
        f"{len(dna)} bp"
    )

    # ========================================================
    # SEQUENCE STATISTICS
    # ========================================================

    statistics = _sequence_statistics(
        dna,
        nucleotide_counts,
        nucleotide_percentages
    )

    print()
    print("Sequence Statistics")
    print("-" * 45)

    print(
        f"Length:               "
        f"{statistics['length']} bp"
    )

    print(
        f"Standard Bases:       "
        f"{statistics['standard_length']}"
    )

    print(
        f"GC Count:             "
        f"{statistics['gc_count']}"
    )

    print(
        f"AT Count:             "
        f"{statistics['at_count']}"
    )

    print(
        f"GC Content:           "
        f"{statistics['gc_content']:.2f}%"
    )

    print(
        f"AT Content:           "
        f"{statistics['at_content']:.2f}%"
    )

    # ========================================================
    # SEQUENCE SUMMARY
    # ========================================================

    print()
    print("Sequence Summary")
    print("-" * 45)

    if sequence_quality["invalid_symbols"] > 0:

        sequence_status = "INVALID"

    elif sequence_quality["ambiguous_bases"] > 0:

        sequence_status = "VALID - AMBIGUOUS"

    else:

        sequence_status = "VALID"

    print(
        f"Status:              "
        f"{sequence_status}"
    )

    print(
        f"GC Content:          "
        f"{gc_content:.2f}%"
    )

    print(
        f"AT Content:          "
        f"{at_content:.2f}%"
    )

    print(
        f"Complexity:          "
        f"{sequence_complexity:.3f}"
    )

    print(
        f"CpG Sites:           "
        f"{cpg_analysis['count']}"
    )

    print(
        f"CpG Islands:         "
        f"{len(cpg_islands)}"
    )

    print(
        f"Homopolymer Regions: "
        f"{len(homopolymers)}"
    )

    print(
        f"Ambiguous Bases:     "
        f"{sequence_quality['ambiguous_bases']}"
    )

    print(
        f"Invalid Symbols:     "
        f"{sequence_quality['invalid_symbols']}"
    )

    # ========================================================
    # GC / AT CONTENT
    # ========================================================

    print()
    print("GC / AT Content")
    print("-" * 45)

    if content_ranges["ambiguous"]:

        print(
            f"GC Content Range:    "
            f"{content_ranges['gc_min']:.2f}% - "
            f"{content_ranges['gc_max']:.2f}%"
        )

        print(
            f"AT Content Range:    "
            f"{content_ranges['at_min']:.2f}% - "
            f"{content_ranges['at_max']:.2f}%"
        )

        print(
            "Content Status:      "
            "Range due to ambiguous IUPAC bases"
        )

    else:

        print(
            f"GC Content:          "
            f"{gc_content:.2f}%"
        )

        print(
            f"AT Content:          "
            f"{at_content:.2f}%"
        )

    # ========================================================
    # NUCLEOTIDE COMPOSITION
    # ========================================================

    print()
    print("Nucleotide Composition")
    print("-" * 45)

    for base in (
        "A",
        "T",
        "C",
        "G"
    ):

        print(
            f"{base}:                  "
            f"{nucleotide_counts[base]:>6}  "
            f"{nucleotide_percentages[base]:>6.2f}%"
        )

    # ========================================================
    # NUCLEOTIDE SKEW
    # ========================================================

    print()
    print("Nucleotide Skew")
    print("-" * 45)

    print(
        f"GC Skew:             "
        f"{gc_skew:+.3f}"
    )

    print(
        f"AT Skew:             "
        f"{at_skew:+.3f}"
    )

    print(
        "Interpretation:      "
        "Positive values indicate relative enrichment "
        "of the first nucleotide in the pair."
    )

    # ========================================================
    # MOLECULAR WEIGHT
    # ========================================================

    print()
    print("Molecular Weight")
    print("-" * 45)

    if molecular_weight["ambiguous"]:

        print(
            f"Estimated dsDNA MW:  "
            f"{molecular_weight['min'] / 1000:,.2f} - "
            f"{molecular_weight['max'] / 1000:,.2f} kDa"
        )

        print(
            "MW Status:           "
            "Range due to ambiguous IUPAC bases"
        )

    else:

        print(
            f"Estimated dsDNA MW:  "
            f"{molecular_weight['min'] / 1000:,.2f} kDa"
        )

    print(
        "Method:              "
        "Composition-based nucleotide residue estimate"
    )

    # ========================================================
    # GC CONTENT PROFILE
    # ========================================================

    print()
    print("GC Content Profile")
    print("-" * 45)

    if len(dna) < profile_window:

        _print_profile_unavailable(
            "GC Content Profile",
            len(dna),
            profile_window
        )

    else:

        print(
            f"Window:              "
            f"{profile_window} bp"
        )

        print(
            f"Step:                "
            f"{profile_step} bp"
        )

        _print_profile_statistics(
            gc_profile_statistics,
            "GC Content",
            "%"
        )

        print()
        print("Profile Preview")
        print("-" * 45)

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

    print()
    print("Sequence Complexity Profile")
    print("-" * 45)

    if len(dna) < profile_window:

        _print_profile_unavailable(
            "Sequence Complexity Profile",
            len(dna),
            profile_window
        )

    else:

        print(
            f"Window:              "
            f"{profile_window} bp"
        )

        print(
            f"Step:                "
            f"{profile_step} bp"
        )

        _print_profile_statistics(
            complexity_profile_statistics,
            "Complexity"
        )

        print()
        print("Profile Preview")
        print("-" * 45)

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

    print()
    print("Homopolymers")
    print("-" * 45)

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

    print()
    print("CpG Analysis")
    print("-" * 45)

    print(
        f"CpG Count:           "
        f"{cpg_analysis['count']}"
    )

    print(
        f"CpG Density:         "
        f"{cpg_analysis['density']:.3f}%"
    )

    print(
        f"Observed/Expected:   "
        f"{cpg_analysis['oe_ratio']:.3f}"
    )

    print(
        f"Ambiguous Contexts:  "
        f"{cpg_analysis['ambiguous_contexts']}"
    )

    # ========================================================
    # CpG ISLANDS
    # ========================================================

    print()
    print("CpG Islands")
    print("-" * 45)

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

    print()
    print("Sequence Quality")
    print("-" * 45)

    print(
        f"Standard Bases:      "
        f"{sequence_quality['standard_bases']}"
    )

    print(
        f"Ambiguous Bases:     "
        f"{sequence_quality['ambiguous_bases']}"
    )

    print(
        f"Invalid Symbols:     "
        f"{sequence_quality['invalid_symbols']}"
    )

    print(
        f"Standard:            "
        f"{sequence_quality['standard_percentage']:.2f}%"
    )

    print(
        f"Ambiguous:           "
        f"{sequence_quality['ambiguous_percentage']:.2f}%"
    )

    print(
        f"Invalid:             "
        f"{sequence_quality['invalid_percentage']:.2f}%"
    )

    if sequence_quality["invalid_positions"]:

        print()
        print(
            "Invalid Positions:   "
            + ", ".join(
                str(item["position_1_based"])
                for item in
                sequence_quality[
                    "invalid_positions"
                ][:20]
            )
        )

    if sequence_quality["ambiguous_positions"]:

        print(
            "Ambiguous Positions: "
            + ", ".join(
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

    print()
    print("IUPAC / Ambiguous Bases")
    print("-" * 45)

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

    print()
    print("Complementary Sequences")
    print("-" * 45)

    print(
        "Original:             5' -> 3'"
    )

    print(
        "Complement:           3' -> 5'"
    )

    print(
        f"                     {_preview(complement)}"
    )

    print(
        "Reverse Complement:   5' -> 3'"
    )

    print(
        f"                     {_preview(reverse_complement)}"
    )

    # ========================================================
    # SEQUENCE PREVIEW
    # ========================================================

    print()
    print("Sequence Preview")
    print("-" * 45)

    print(
        f"5' -> 3':            "
        f"{_preview(dna)}"
    )

    print(
        f"Length shown:        "
        f"{min(len(dna), 80)} / "
        f"{len(dna)} bp"
    )

    # ========================================================
    # INTERPRETATION
    # ========================================================

    print()
    print("Interpretation")
    print("-" * 45)

    if sequence_quality["invalid_symbols"] > 0:

        print(
            "Sequence status:     INVALID"
        )

        print(
            "The sequence contains invalid "
            "symbols and should be reviewed "
            "before downstream analysis."
        )

    elif sequence_quality["ambiguous_bases"] > 0:

        print(
            "Sequence status:     VALID - AMBIGUOUS"
        )

        print(
            "The sequence is valid but contains "
            "IUPAC ambiguity codes. Some values "
            "are therefore represented as ranges."
        )

    else:

        print(
            "Sequence status:     VALID"
        )

        print(
            "The sequence contains only standard "
            "DNA bases and passed sequence-level "
            "quality screening."
        )

    if sequence_complexity >= 0.75:

        print(
            "Complexity class:    HIGH"
        )

    elif sequence_complexity >= 0.50:

        print(
            "Complexity class:    MODERATE"
        )

    else:

        print(
            "Complexity class:    LOW"
        )

    if gc_content >= 60:

        print(
            "GC composition:      Relatively high"
        )

    elif gc_content <= 40:

        print(
            "GC composition:      Relatively low"
        )

    else:

        print(
            "GC composition:      Relatively balanced"
        )

    print(
        f"CpG enrichment:      "
        f"O/E {cpg_analysis['oe_ratio']:.3f}"
    )

    print(
        f"CpG island screen:   "
        f"{'Detected' if cpg_islands else 'Not detected'}"
    )

    print(
        f"Homopolymers:        "
        f"{'Detected' if homopolymers else 'None detected'}"
    )

    print()
    print(
        "Profile analyses are sequence-level "
        "descriptors and do not establish biological "
        "function or experimental behavior."
    )

    print()
    print("=" * 70)