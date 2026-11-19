"""analyze_orfs, density, frame statistics, relationships, ranking and comparison."""

from __future__ import annotations

from analyses._orf.codes import (
    codon_table_for,
    DEFAULT_MIN_ORF_LENGTH,
    DEFAULT_RANKING,
    DEFAULT_REPORT_NESTED,
    DEFAULT_START_CODON_POLICY,
    _genetic_code_name,
    ORF_RANKING_MODES,
    start_codons_for,
    _stop_codons_of,
)
from analyses._orf.finder import find_orfs, _orf_sort_key


# ============================================================
# ORF DENSITY
# ============================================================

def calculate_orf_density(
    orfs,
    sequence_length
):

    if sequence_length <= 0:
        return 0.0

    if not orfs:
        return 0.0

    total_orf_bases = sum(
        orf["length"]
        for orf in orfs
    )

    return (
        total_orf_bases
        / sequence_length
    ) * 100


# ============================================================
# FRAME STATISTICS
# ============================================================

def calculate_frame_statistics(
    orfs
):

    frames = {
        "+1": [],
        "+2": [],
        "+3": [],
        "-1": [],
        "-2": [],
        "-3": []
    }

    for orf in orfs:

        key = (
            f"{orf['frame']:+d}"
        )

        if key in frames:
            frames[key].append(
                orf
            )

    statistics = {}

    for frame, frame_orfs in (
        frames.items()
    ):

        if not frame_orfs:

            statistics[frame] = {
                "count": 0,
                "complete": 0,
                "partial": 0,
                "longest": None,
                "average_length": 0.0
            }

            continue

        complete = sum(
            1
            for orf in frame_orfs
            if orf["status"] == "Complete"
        )

        longest = max(
            frame_orfs,
            key=lambda orf: orf["length"]
        )

        average_length = (
            sum(
                orf["length"]
                for orf in frame_orfs
            )
            / len(frame_orfs)
        )

        statistics[frame] = {
            "count":
                len(frame_orfs),

            "complete":
                complete,

            "partial":
                len(frame_orfs) - complete,

            "longest":
                longest,

            "average_length":
                average_length
        }

    return statistics


# ============================================================
# ORF SUMMARY
# ============================================================

def summarize_orfs(
    orfs,
    sequence_length=None
):

    if not orfs:

        return {
            "total": 0,
            "forward": 0,
            "reverse": 0,
            "complete": 0,
            "partial": 0,
            "canonical_start": 0,
            "alternative_start": 0,
            "ambiguous": 0,
            "translation_safe": 0,
            "kozak_available": 0,
            "shine_dalgarno_detected": 0,
            "average_length": 0.0,
            "longest": None,
            "orf_density": 0.0,
            "frame_statistics":
                calculate_frame_statistics([])
        }

    forward = sum(
        1
        for orf in orfs
        if orf["strand"] == "+"
    )

    reverse = sum(
        1
        for orf in orfs
        if orf["strand"] == "-"
    )

    complete = sum(
        1
        for orf in orfs
        if orf["status"] == "Complete"
    )

    partial = (
        len(orfs) - complete
    )

    canonical_start = sum(
        1
        for orf in orfs
        if orf["start_type"]
        == "Canonical"
    )

    alternative_start = sum(
        1
        for orf in orfs
        if orf["start_type"]
        == "Alternative"
    )

    ambiguous = sum(
        1
        for orf in orfs
        if orf["ambiguous_codon_count"] > 0
    )

    translation_safe = sum(
        1
        for orf in orfs
        if orf["translation_safe"]
    )

    kozak_available = sum(
        1
        for orf in orfs
        if (
            orf["start_context"]
            ["kozak"]["available"]
        )
    )

    shine_dalgarno_detected = sum(
        1
        for orf in orfs
        if (
            orf["start_context"]
            ["shine_dalgarno"]["detected"]
        )
    )

    average_length = (
        sum(
            orf["length"]
            for orf in orfs
        )
        / len(orfs)
    )

    longest = max(
        orfs,
        key=lambda orf: orf["length"]
    )

    if sequence_length is None:
        sequence_length = 0

    return {
        "total":
            len(orfs),

        "forward":
            forward,

        "reverse":
            reverse,

        "complete":
            complete,

        "partial":
            partial,

        "canonical_start":
            canonical_start,

        "alternative_start":
            alternative_start,

        "ambiguous":
            ambiguous,

        "translation_safe":
            translation_safe,

        "kozak_available":
            kozak_available,

        "shine_dalgarno_detected":
            shine_dalgarno_detected,

        "average_length":
            average_length,

        "longest":
            longest,

        "orf_density":
            calculate_orf_density(
                orfs,
                sequence_length
            ),

        "frame_statistics":
            calculate_frame_statistics(
                orfs
            )
    }


# ============================================================
# UNIFIED ORF ANALYSIS
# ============================================================

def analyze_orfs(
    dna,
    include_partial=True,
    min_length=DEFAULT_MIN_ORF_LENGTH,
    start_codon_policy=None,
    nested=DEFAULT_REPORT_NESTED,
    ranking=DEFAULT_RANKING,
    genetic_code=1,
    topology="linear",
    include_5prime_partial=False,
):

    """
    Unified high-level ORF analysis.

    Returns all ORFs together with summary, frame statistics,
    density, ranked candidates and the search parameters used.
    """

    orfs = find_orfs(
        dna,
        include_partial=include_partial,
        min_length=min_length,
        start_codon_policy=start_codon_policy,
        nested=nested,
        ranking=ranking,
        genetic_code=genetic_code,
        topology=topology,
        include_5prime_partial=include_5prime_partial,
    )

    summary = summarize_orfs(
        orfs,
        len(dna)
    )

    return {
        "sequence_length":
            len(dna),

        "minimum_length":
            min_length,

        "include_partial":
            include_partial,

        "parameters": {
            "min_length_nt": min_length,
            "include_partial": include_partial,
            "include_5prime_partial": include_5prime_partial,
            "topology": str(topology).lower(),
            "start_codon_policy": start_codon_policy or DEFAULT_START_CODON_POLICY,
            "start_codons": sorted(start_codons_for(genetic_code, start_codon_policy)),
            "stop_codons": sorted(_stop_codons_of(codon_table_for(genetic_code))),
            "nested": nested,
            "ranking": ranking,
            "genetic_code": _genetic_code_name(genetic_code),
            "genetic_code_table": int(genetic_code or 1),
        },

        "orfs":
            orfs,

        "summary":
            summary,

        "frame_statistics":
            summary["frame_statistics"],

        "orf_density":
            summary["orf_density"],

        "longest_orf":
            summary["longest"],

        "top_candidate":
            orfs[0]
            if orfs
            else None
    }


# ============================================================
# SEQUENCE PREVIEW
# ============================================================

def create_orf_preview(
    sequence,
    max_length=60
):

    if len(sequence) <= max_length:
        return sequence

    half = max_length // 2

    return (
        sequence[:half]
        + "..."
        + sequence[-half:]
    )


# ============================================================
# ORF REPORT
# ============================================================

# The terminal dialogs of this module are in ui/orf_menu.py (P2-2).


# ============================================================
# ORF OVERLAP / NESTED ANALYSIS
# ============================================================

def _orf_genomic_interval(orf):
    """Return a normalized 0-based inclusive genomic interval."""
    return min(orf["start"], orf["stop"]), max(orf["start"], orf["stop"])


def analyze_orf_relationships(orfs):
    """
    Detect overlapping and nested ORFs.

    Relationships are sequence-coordinate relationships only. They do not
    imply that the ORFs are biologically functional or simultaneously used.
    """
    relationships = []
    overlap_pairs = []
    nested_pairs = []

    for i in range(len(orfs)):
        a = orfs[i]
        a_start, a_end = _orf_genomic_interval(a)

        for j in range(i + 1, len(orfs)):
            b = orfs[j]
            b_start, b_end = _orf_genomic_interval(b)

            overlap_start = max(a_start, b_start)
            overlap_end = min(a_end, b_end)

            if overlap_start > overlap_end:
                continue

            overlap_length = overlap_end - overlap_start + 1

            a_contains_b = (
                a_start <= b_start and b_end <= a_end
            )
            b_contains_a = (
                b_start <= a_start and a_end <= b_end
            )

            # (outer, inner) when one ORF lies inside the other.
            nesting = None
            if a_contains_b and not b_contains_a:
                relation = "Nested"
                nesting = (a, b)
                nested_pairs.append((a, b))
            elif b_contains_a and not a_contains_b:
                relation = "Nested"
                nesting = (b, a)
                nested_pairs.append((b, a))
            else:
                relation = "Overlapping"

            overlap_fraction_a = overlap_length / a["length"]
            overlap_fraction_b = overlap_length / b["length"]

            record = {
                "orf_a_rank": a["rank"],
                "orf_b_rank": b["rank"],
                "relation": relation,
                "overlap_start": overlap_start + 1,
                "overlap_end": overlap_end + 1,
                "overlap_length": overlap_length,
                "overlap_fraction_a": overlap_fraction_a,
                "overlap_fraction_b": overlap_fraction_b,
                "same_strand": a["strand"] == b["strand"],
                "same_frame": a["frame"] == b["frame"],
            }

            if nesting is not None:
                record["outer_orf_rank"] = nesting[0]["rank"]
                record["inner_orf_rank"] = nesting[1]["rank"]

            relationships.append(record)

            if relation == "Overlapping":
                overlap_pairs.append(record)

    return {
        "relationship_count": len(relationships),
        "overlap_count": len(overlap_pairs),
        "nested_count": len(nested_pairs),
        "relationships": relationships,
        "overlaps": overlap_pairs,
        "nested": [
            relationship
            for relationship in relationships
            if relationship["relation"] == "Nested"
        ],
    }


# ============================================================
# ORF RANKING / COMPARISON
# ============================================================

def rank_orf_candidates(orfs, ranking=DEFAULT_RANKING):
    """Return ORFs ranked by length (default) or the heuristic score."""
    if ranking not in ORF_RANKING_MODES:
        raise ValueError(f"ranking must be one of {ORF_RANKING_MODES}.")
    ranked = sorted(
        orfs,
        key=lambda orf: _orf_sort_key(orf, ranking),
        reverse=True,
    )

    for rank, orf in enumerate(ranked, start=1):
        orf["rank"] = rank

    return ranked


def compare_orf_candidates(orfs, top_n=5, ranking=None):
    """Create a compact comparison table for the best ORF candidates.

    ranking=None keeps the existing ranks (no re-ranking); pass "length" or
    "heuristic" to re-rank.
    """
    if ranking is None and all("rank" in orf for orf in orfs):
        ranked = sorted(orfs, key=lambda orf: orf["rank"])
    else:
        ranked = rank_orf_candidates(orfs, ranking=ranking or DEFAULT_RANKING)
    return [
        {
            "rank": orf["rank"],
            "strand": orf["strand"],
            "frame": orf["frame"],
            "length": orf["length"],
            "protein_length": orf["protein_length"],
            "start_codon": orf["start_codon"],
            "status": orf["status"],
            "heuristic_score": orf["heuristic_score"]["score"],
            "heuristic_classification": orf["heuristic_score"]["classification"],
            "cds_quality_score": orf["cds_quality"]["score"],
            "cds_quality_classification": orf["cds_quality"]["classification"],
            "translation_safe": orf["translation_safe"],
            "ambiguous_codons": orf["ambiguous_codon_count"],
        }
        for orf in ranked[:max(0, top_n)]
    ]
