"""The feature score, the analysis summary, analyze_rna and the preview."""

from __future__ import annotations

from analyses._rna.bases import (
    count_rna_bases,
    get_rna_complement,
    get_rna_reverse_complement,
    transcribe_rna,
)
from analyses._rna.composition import (
    analyze_rna_ambiguity,
    analyze_rna_quality,
    calculate_rna_au_content,
    calculate_rna_au_skew,
    calculate_rna_base_percentages,
    calculate_rna_complexity,
    calculate_rna_complexity_profile,
    calculate_rna_content_ranges,
    calculate_rna_gc_content,
    calculate_rna_gc_profile,
    calculate_rna_gc_skew,
    calculate_rna_molecular_weight,
    find_rna_homopolymer_regions,
    find_rna_homopolymers,
)
from analyses._rna.coding import (
    analyze_rna_codon_usage,
    analyze_rna_kozak_context,
    analyze_rna_shine_dalgarno,
    analyze_rna_three_frames,
    calculate_rna_codon_diversity,
    find_rna_start_stop_landscape,
)
from analyses._rna.structure import summarize_rna_structure


# ============================================================
# RNA FEATURE SCORE
# ============================================================

def calculate_rna_feature_score(analysis):
    """
    Heuristic sequence-level feature score.

    This is NOT a biological function or expression predictor.
    """

    score = 100
    reasons = []

    quality = analysis["quality"]
    ambiguity = analysis["ambiguity"]

    if quality["invalid_count"] > 0:
        score -= 40
        reasons.append("Invalid RNA symbols detected.")

    if ambiguity["ambiguous_percentage"] == 0:
        reasons.append("No IUPAC ambiguity detected.")

    elif ambiguity["ambiguous_percentage"] <= 5:
        score -= 5
        reasons.append("Low IUPAC ambiguity.")

    elif ambiguity["ambiguous_percentage"] <= 15:
        score -= 15
        reasons.append("Moderate IUPAC ambiguity.")

    else:
        score -= 25
        reasons.append("High IUPAC ambiguity.")

    complexity = analysis["complexity"]

    if complexity >= 0.75:
        reasons.append("High sequence complexity.")

    elif complexity >= 0.50:
        score -= 5
        reasons.append("Moderate sequence complexity.")

    else:
        score -= 15
        reasons.append("Low sequence complexity.")

    homopolymer_count = len(
        analysis["homopolymers"]
    )

    if homopolymer_count == 0:
        reasons.append("No long homopolymers detected.")

    elif homopolymer_count <= 2:
        score -= 5
        reasons.append("Few homopolymer regions detected.")

    else:
        score -= 10
        reasons.append("Multiple homopolymer regions detected.")

    hairpin_count = analysis["structure"]["hairpin_count"]

    if hairpin_count > 0:
        reasons.append(
            "Candidate stem-loop motifs detected; "
            "thermodynamic stability not calculated."
        )

    landscape = analysis["start_stop"]

    if landscape["start_count"] > 0:
        reasons.append(
            f"{landscape['start_count']} AUG start codon(s) detected."
        )

    if landscape["stop_count"] > 0:
        reasons.append(
            f"{landscape['stop_count']} stop codon(s) detected."
        )

    score = max(
        0,
        min(100, score)
    )

    if score >= 80:
        classification = "High"

    elif score >= 60:
        classification = "Moderate"

    elif score >= 40:
        classification = "Low"

    else:
        classification = "Very Low"

    return {
        "score": score,
        "classification": classification,
        "reasons": reasons
    }


# ============================================================
# RNA SUMMARY
# ============================================================

def summarize_rna_analysis(analysis):
    """Create a compact sequence-level RNA summary."""

    quality = analysis["quality"]
    ambiguity = analysis["ambiguity"]
    landscape = analysis["start_stop"]

    return {
        "length": analysis["length"],
        "gc_content": analysis["gc_content"],
        "au_content": analysis["au_content"],
        "ambiguous_count": ambiguity["ambiguous_count"],
        "ambiguous_percentage": ambiguity["ambiguous_percentage"],
        "invalid_count": quality["invalid_count"],
        "start_count": landscape["start_count"],
        "stop_count": landscape["stop_count"],
        "complexity": analysis["complexity"],
        "homopolymer_count": len(
            analysis["homopolymers"]
        ),
        "hairpin_count": analysis["structure"]["hairpin_count"],
        "feature_score": analysis["feature_score"]["score"],
        "feature_classification":
            analysis["feature_score"]["classification"]
    }


# ============================================================
# THREE-FRAME TRANSLATION REPORT
# ============================================================

# The terminal dialogs of this module are in ui/rna_menu.py (P2-2).


# ============================================================
# RNA FEATURE REPORT
# ============================================================


# ============================================================
# COMPLETE RNA ANALYSIS
# ============================================================

def analyze_rna(
    dna,
    profile_window=100,
    profile_step=1,
    genetic_code=1
):

    rna = transcribe_rna(
        dna
    )

    base_counts = (
        count_rna_bases(
            rna
        )
    )

    base_percentages = (
        calculate_rna_base_percentages(
            base_counts
        )
    )

    ambiguity = (
        analyze_rna_ambiguity(
            rna
        )
    )

    quality = (
        analyze_rna_quality(
            rna
        )
    )

    content_ranges = (
        calculate_rna_content_ranges(
            rna
        )
    )

    molecular_weight = (
        calculate_rna_molecular_weight(
            rna
        )
    )

    codon_usage = (
        analyze_rna_codon_usage(
            rna
        )
    )

    structure = (
        summarize_rna_structure(
            rna
        )
    )

    start_stop = (
        find_rna_start_stop_landscape(
            rna
        )
    )

    # Build the core analysis first so the heuristic feature
    # score can use all previously calculated sequence features.
    analysis = {
        "rna":
            rna,

        "length":
            len(rna),

        "base_counts":
            base_counts,

        "base_percentages":
            base_percentages,

        "gc_content":
            calculate_rna_gc_content(
                rna
            ),

        "au_content":
            calculate_rna_au_content(
                rna
            ),

        "content_ranges":
            content_ranges,

        "gc_skew":
            calculate_rna_gc_skew(
                rna
            ),

        "au_skew":
            calculate_rna_au_skew(
                rna
            ),

        "gc_profile":
            calculate_rna_gc_profile(
                rna,
                profile_window,
                profile_step
            ),

        "complexity":
            calculate_rna_complexity(
                rna
            ),

        "complexity_profile":
            calculate_rna_complexity_profile(
                rna,
                profile_window,
                profile_step
            ),

        "ambiguity":
            ambiguity,

        "molecular_weight":
            molecular_weight,

        "quality":
            quality,

        "homopolymers":
            find_rna_homopolymers(
                rna
            ),

        "codon_usage":
            codon_usage,

        "codon_diversity":
            calculate_rna_codon_diversity(
                rna
            ),

        "kozak":
            analyze_rna_kozak_context(
                rna
            ),

        "shine_dalgarno":
            analyze_rna_shine_dalgarno(
                rna
            ),

        "poly_a":
            find_rna_homopolymer_regions(
                rna,
                "A"
            ),

        "poly_u":
            find_rna_homopolymer_regions(
                rna,
                "U"
            ),

        "complement":
            get_rna_complement(
                rna
            ),

        "reverse_complement":
            get_rna_reverse_complement(
                rna
            ),

        "structure":
            structure,

        "start_stop":
            start_stop
    }

    analysis["feature_score"] = (
        calculate_rna_feature_score(
            analysis
        )
    )

    analysis["summary"] = (
        summarize_rna_analysis(
            analysis
        )
    )

    analysis["genetic_code"] = (
        1 if genetic_code in (None, 1) else int(genetic_code)
    )

    analysis["translation_frames"] = (
        analyze_rna_three_frames(
            rna,
            genetic_code=genetic_code,
        )
    )

    return analysis


# ============================================================
# SEQUENCE PREVIEW
# ============================================================

def create_rna_preview(
    rna,
    max_length=60
):

    if len(rna) <= max_length:
        return rna

    half = max_length // 2

    return (
        rna[:half]
        + "..."
        + rna[-half:]
    )
