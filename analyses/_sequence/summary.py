"""Sequence forms, interpretation, the summary and analyze_sequence."""

from __future__ import annotations

from analyses._sequence.bases import (
    count_ambiguous_bases,
    count_nucleotides,
    get_complement,
    normalize_dna_sequence,
)
from analyses._sequence.composition import (
    analyze_sequence_quality,
    calculate_at_content,
    calculate_at_skew,
    calculate_average_dna_molecular_weight,
    calculate_dna_content_ranges,
    calculate_dna_molecular_weight,
    calculate_gc_content,
    calculate_gc_skew,
    calculate_nucleotide_percentages,
    calculate_sequence_statistics,
)
from analyses._sequence.profiles import (
    analyze_complexity_profile,
    analyze_gc_profile,
    calculate_sequence_complexity,
    describe_profile_availability,
    find_homopolymers,
)
from analyses._sequence.cpg import analyze_cpg, DEFAULT_CPG_ISLAND_CRITERIA, find_cpg_islands


# ============================================================
# COMPLEMENT ORIENTATION
# ============================================================

def get_sequence_forms(dna):
    """
    Return sequence forms with explicit strand orientation.

    Original and reverse complement are represented 5'->3'.
    The direct complement is represented 3'->5' relative to the
    original strand.
    """

    complement = get_complement(dna)
    reverse_complement = complement[::-1]

    return {
        "original": {
            "sequence": dna,
            "orientation": "5' -> 3'"
        },
        "complement": {
            "sequence": complement,
            "orientation": "3' -> 5'"
        },
        "reverse_complement": {
            "sequence": reverse_complement,
            "orientation": "5' -> 3'"
        }
    }


# ============================================================
# INTERPRETATION
# ============================================================

def create_sequence_interpretation(
    dna,
    quality,
    sequence_complexity,
    cpg,
    cpg_islands,
    homopolymers,
    content_ranges
):
    """
    Generate a concise, sequence-level interpretation.

    This describes sequence properties only. It does not infer
    biological function, expression, pathogenicity, or phenotype.
    """

    statements = []

    if quality["invalid_symbols"] == 0:
        statements.append(
            "The sequence contains only recognized DNA/IUPAC symbols "
            "and passed symbol-level quality screening."
        )
    else:
        statements.append(
            f"The sequence contains {quality['invalid_symbols']} "
            "invalid symbol(s) and requires cleanup before "
            "unambiguous downstream analysis."
        )

    if quality["ambiguous_bases"] > 0:
        statements.append(
            f"{quality['ambiguous_bases']} IUPAC ambiguous base(s) "
            "are present; analyses that require an exact base "
            "sequence may therefore have unresolved positions."
        )

    gc_min = content_ranges["gc_min"]
    gc_max = content_ranges["gc_max"]

    if abs(gc_max - gc_min) < 1e-12:
        gc_value = gc_min
        if gc_value < 40:
            gc_class = "relatively low"
        elif gc_value <= 60:
            gc_class = "relatively balanced"
        else:
            gc_class = "relatively high"

        statements.append(
            f"GC content is {gc_value:.2f}% and is {gc_class} "
            "for a simple sequence-level composition descriptor."
        )
    else:
        statements.append(
            f"GC content is constrained to {gc_min:.2f}%–"
            f"{gc_max:.2f}% because ambiguous bases are present."
        )

    if sequence_complexity >= 0.75:
        complexity_class = "high"
    elif sequence_complexity >= 0.50:
        complexity_class = "moderate"
    else:
        complexity_class = "low"

    statements.append(
        f"Sequence complexity is {complexity_class} "
        f"({sequence_complexity:.3f}): the normalized Shannon "
        "entropy of A/C/G/T composition. It measures base "
        "composition only and does not detect repeats; a tandem "
        "repeat of all four bases scores 1.000."
    )

    if cpg["count"] > 0:
        statements.append(
            f"{cpg['count']} definite CpG site(s) were detected; "
            f"CpG density is {cpg['density']:.3f}% of evaluated "
            "positions."
        )
    else:
        statements.append(
            "No definite CpG dinucleotides were detected."
        )

    if cpg_islands:
        statements.append(
            f"{len(cpg_islands)} CpG island-like region(s) passed "
            "the configured classical screening criteria."
        )
    else:
        statements.append(
            "No CpG island-like region passed the configured "
            "criteria."
        )

    if homopolymers:
        longest = max(
            homopolymers,
            key=lambda item: item["length"]
        )
        statements.append(
            f"{len(homopolymers)} homopolymer region(s) were "
            f"detected; the longest is {longest['length']} bp."
        )
    else:
        statements.append(
            "No homopolymer regions met the configured minimum "
            "length."
        )

    statements.append(
        "These results are sequence-level descriptors and do not "
        "establish biological function or experimental behavior."
    )

    return statements



# ============================================================
# SEQUENCE SUMMARY
# ============================================================


def create_sequence_summary(
    dna,
    nucleotide_counts=None,
    nucleotide_percentages=None,
    quality=None,
    content_ranges=None,
    cpg=None,
    cpg_islands=None,
    homopolymers=None,
    molecular_weight=None,
    sequence_complexity=None
):

    if nucleotide_counts is None:

        nucleotide_counts = (
            count_nucleotides(dna)
        )

    if nucleotide_percentages is None:

        nucleotide_percentages = (
            calculate_nucleotide_percentages(
                nucleotide_counts
            )
        )

    if quality is None:

        quality = (
            analyze_sequence_quality(
                dna
            )
        )

    if content_ranges is None:

        content_ranges = (
            calculate_dna_content_ranges(
                dna
            )
        )

    if cpg is None:

        cpg = analyze_cpg(dna)

    if cpg_islands is None:

        cpg_islands = find_cpg_islands(
            dna
        )

    if homopolymers is None:

        homopolymers = find_homopolymers(
            dna
        )

    if molecular_weight is None:

        molecular_weight = (
            calculate_dna_molecular_weight(
                dna
            )
        )

    if sequence_complexity is None:

        sequence_complexity = (
            calculate_sequence_complexity(
                dna
            )
        )

    if not quality["valid"]:

        status = "INVALID"

    elif quality["ambiguous_bases"] > 0:

        status = "VALID - AMBIGUOUS"

    else:

        status = "VALID"

    if content_ranges["ambiguous"]:

        gc_min = content_ranges["gc_min"]
        gc_max = content_ranges["gc_max"]

        at_min = content_ranges["at_min"]
        at_max = content_ranges["at_max"]

    else:

        gc_value = calculate_gc_content(dna)
        at_value = calculate_at_content(dna)

        gc_min = gc_value
        gc_max = gc_value

        at_min = at_value
        at_max = at_value

    return {

        "length":
            len(dna),

        "status":
            status,

        "gc_content":
            calculate_gc_content(dna),

        "gc_min":
            gc_min,

        "gc_max":
            gc_max,

        "at_content":
            calculate_at_content(dna),

        "at_min":
            at_min,

        "at_max":
            at_max,

        "sequence_complexity":
            sequence_complexity,

        "cpg_count":
            cpg["count"],

        "cpg_density":
            cpg["density"],

        "cpg_oe_ratio":
            cpg["oe_ratio"],

        "cpg_islands":
            len(cpg_islands),

        "homopolymer_regions":
            len(homopolymers),

        "ambiguous_bases":
            quality["ambiguous_bases"],

        "invalid_symbols":
            quality["invalid_symbols"],

        "molecular_weight_min":
            molecular_weight["min"],

        "molecular_weight_max":
            molecular_weight["max"],

        "nucleotide_counts":
            nucleotide_counts,

        "nucleotide_percentages":
            nucleotide_percentages
    }


# ============================================================
# UNIFIED SEQUENCE ANALYSIS
# ============================================================


def analyze_sequence(
    dna,
    profile_window=100,
    profile_step=1
):

    # --------------------------------------------------------
    # INPUT NORMALIZATION
    # --------------------------------------------------------

    dna = normalize_dna_sequence(
        dna
    )

    # --------------------------------------------------------
    # BASIC ANALYSIS
    # --------------------------------------------------------

    nucleotide_counts = (
        count_nucleotides(dna)
    )

    nucleotide_percentages = (
        calculate_nucleotide_percentages(
            nucleotide_counts
        )
    )

    quality = (
        analyze_sequence_quality(
            dna
        )
    )

    content_ranges = (
        calculate_dna_content_ranges(
            dna
        )
    )

    cpg = analyze_cpg(dna)

    homopolymers = find_homopolymers(
        dna
    )

    cpg_islands = find_cpg_islands(
        dna
    )

    molecular_weight = (
        calculate_dna_molecular_weight(
            dna
        )
    )

    sequence_statistics = (
        calculate_sequence_statistics(
            dna,
            nucleotide_counts,
            nucleotide_percentages
        )
    )

    sequence_complexity = (
        calculate_sequence_complexity(
            dna
        )
    )

    # --------------------------------------------------------
    # PROFILES
    # --------------------------------------------------------

    gc_profile = analyze_gc_profile(
        dna,
        profile_window,
        profile_step
    )

    complexity_profile = (
        analyze_complexity_profile(
            dna,
            profile_window,
            profile_step
        )
    )

    gc_profile_availability = (
        describe_profile_availability(
            dna,
            profile_window,
            "GC content profile"
        )
    )

    complexity_profile_availability = (
        describe_profile_availability(
            dna,
            profile_window,
            "Sequence complexity profile"
        )
    )

    # --------------------------------------------------------
    # SEQUENCE FORMS
    # --------------------------------------------------------

    sequence_forms = get_sequence_forms(
        dna
    )

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    summary = create_sequence_summary(

        dna,

        nucleotide_counts=
            nucleotide_counts,

        nucleotide_percentages=
            nucleotide_percentages,

        quality=
            quality,

        content_ranges=
            content_ranges,

        cpg=
            cpg,

        cpg_islands=
            cpg_islands,

        homopolymers=
            homopolymers,

        molecular_weight=
            molecular_weight,

        sequence_complexity=
            sequence_complexity
    )

    # --------------------------------------------------------
    # INTERPRETATION
    # --------------------------------------------------------

    interpretation = (
        create_sequence_interpretation(
            dna,
            quality,
            sequence_complexity,
            cpg,
            cpg_islands,
            homopolymers,
            content_ranges
        )
    )

    # --------------------------------------------------------
    # FINAL RESULT
    # --------------------------------------------------------

    return {

        # ====================================================
        # SEQUENCE INFORMATION
        # ====================================================

        "length":
            len(dna),

        "sequence_length":
            len(dna),

        "sequence_statistics":
            sequence_statistics,

        # ====================================================
        # NUCLEOTIDE COMPOSITION
        # ====================================================

        "nucleotide_counts":
            nucleotide_counts,

        "nucleotide_percentages":
            nucleotide_percentages,

        # ====================================================
        # GC / AT
        # ====================================================

        "gc_content":
            calculate_gc_content(dna),

        "at_content":
            calculate_at_content(dna),

        "gc_content_range":
            content_ranges,

        "at_content_range":
            content_ranges,

        # ====================================================
        # SKEW
        # ====================================================

        "gc_skew":
            calculate_gc_skew(dna),

        "at_skew":
            calculate_at_skew(dna),

        # ====================================================
        # COMPLEXITY
        # ====================================================

        "sequence_complexity":
            sequence_complexity,

        "complexity_profile":
            complexity_profile["profile"],

        "complexity_profile_statistics":
            complexity_profile["statistics"],

        "complexity_profile_availability":
            complexity_profile_availability,

        # ====================================================
        # GC PROFILE
        # ====================================================

        "gc_profile":
            gc_profile["profile"],

        "gc_profile_statistics":
            gc_profile["statistics"],

        "gc_profile_availability":
            gc_profile_availability,

        # ====================================================
        # HOMOPOLYMERS
        # ====================================================

        "homopolymers":
            homopolymers,

        # ====================================================
        # CpG
        # ====================================================

        "cpg":
            cpg,

        "cpg_islands":
            cpg_islands,

        "cpg_island_criteria":
            {
                "window_size":
                    profile_window
                    if False
                    else DEFAULT_CPG_ISLAND_CRITERIA[
                        "window_size"
                    ],

                "min_gc_content":
                    DEFAULT_CPG_ISLAND_CRITERIA[
                        "min_gc_content"
                    ],

                "min_oe_ratio":
                    DEFAULT_CPG_ISLAND_CRITERIA[
                        "min_oe_ratio"
                    ],

                "definition":
                    DEFAULT_CPG_ISLAND_CRITERIA[
                        "definition"
                    ]
            },

        # ====================================================
        # MOLECULAR WEIGHT
        # ====================================================

        "molecular_weight":
            molecular_weight,

        "average_molecular_weight":
            calculate_average_dna_molecular_weight(
                dna
            ),

        # ====================================================
        # QUALITY
        # ====================================================

        "quality":
            quality,

        "ambiguous_bases":
            count_ambiguous_bases(dna),

        # ====================================================
        # SEQUENCE FORMS / ORIENTATION
        # ====================================================

        "complement":
            sequence_forms[
                "complement"
            ][
                "sequence"
            ],

        "reverse_complement":
            sequence_forms[
                "reverse_complement"
            ][
                "sequence"
            ],

        "sequence_forms":
            sequence_forms,

        # ====================================================
        # INTERPRETATION
        # ====================================================

        "interpretation":
            interpretation,

        # ====================================================
        # SUMMARY
        # ====================================================

        "summary":
            summary
    }
