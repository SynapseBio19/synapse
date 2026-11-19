"""Translation and analysis of a selected ORF, and the enhanced analysis."""

from __future__ import annotations

from analyses._orf.codes import (
    codon_table_for,
    DEFAULT_MIN_ORF_LENGTH,
    DEFAULT_RANKING,
    DEFAULT_REPORT_NESTED,
    DNA_IUPAC,
    expand_iupac_codon,
    is_initiation_codon,
)
from analyses._orf.finder import orf_has_start_codon
from analyses._orf.summary import (
    analyze_orf_relationships,
    analyze_orfs,
    compare_orf_candidates,
    rank_orf_candidates,
)


# ============================================================
# ORF TRANSLATION PREVIEW
# ============================================================

def translate_orf_sequence(orf_sequence, initiator_start=True, genetic_code=None):
    """
    Translate an ORF using the local genetic code.

    orf_sequence is an ORF dict from find_orfs() or a DNA string. An ORF is
    translated with its own "genetic_code" (audit D2); an explicit,
    different genetic_code raises ValueError instead of being used
    silently. A string carries no table, so genetic_code defaults to 1.

    Ambiguous codons are translated as '?' when they do not resolve to a
    single amino acid. STOP terminates the returned protein and is not added
    as an amino-acid character.

    initiator_start=True (default, because every Synapse ORF begins at its
    start codon): a first codon that is an initiation codon of the genetic
    code (is_initiation_codon; table 1: ATG, CTG, GTG, TTG) is decoded as
    methionine, as initiator tRNA does (audits P0-8, D1). Any other first
    codon, and every later codon, is read with the elongation code.
    """
    if isinstance(orf_sequence, dict):
        orf_code = orf_sequence.get("genetic_code") or 1
        if genetic_code is not None and int(genetic_code) != int(orf_code):
            raise ValueError(
                f"The ORF was found with NCBI table {orf_code}; "
                f"translating it with table {genetic_code} would be inconsistent."
            )
        genetic_code = orf_code
        if not orf_has_start_codon(orf_sequence):
            initiator_start = False
        orf_sequence = orf_sequence["sequence"]

    if genetic_code is None:
        genetic_code = 1

    codon_table = codon_table_for(genetic_code)
    protein = []

    for position in range(0, len(orf_sequence) - 2, 3):
        codon = orf_sequence[position:position + 3].upper()

        if (
            position == 0
            and initiator_start
            and is_initiation_codon(codon, genetic_code)
        ):
            protein.append("M")
            continue

        if any(base not in DNA_IUPAC for base in codon):
            protein.append("?")
            continue

        expanded = expand_iupac_codon(codon)
        translations = {
            codon_table.get(expanded_codon)
            for expanded_codon in expanded
            if codon_table.get(expanded_codon) is not None
        }

        if translations == {"STOP"}:
            break

        if len(translations) == 1 and "STOP" not in translations:
            protein.append(next(iter(translations)))
        else:
            protein.append("?")

    return "".join(protein)


def analyze_selected_orf(orf):
    """Return a focused downstream-ready analysis for one selected ORF."""
    protein = translate_orf_sequence(orf)

    return {
        "rank": orf["rank"],
        "strand": orf["strand"],
        "frame": orf["frame"],
        "start": orf["start"],
        "stop": orf["stop"],
        "length": orf["length"],
        "status": orf["status"],
        "sequence": orf["sequence"],
        "protein_sequence": protein,
        "protein_length": len(protein),
        "heuristic_score": orf["heuristic_score"],
        "cds_quality": orf["cds_quality"],
        "start_context": orf["start_context"],
        "codon_usage": orf["codon_usage"],
        "translation_safe": orf["translation_safe"],
    }


# ============================================================
# ORF FEATURE SUMMARY
# ============================================================

def calculate_orf_feature_summary(orfs):
    """Summarize high-value candidate features without claiming function."""
    if not orfs:
        return {
            "high_heuristic_score": 0,
            "moderate_heuristic_score": 0,
            "low_heuristic_score": 0,
            "canonical_complete": 0,
            "alternative_complete": 0,
            "translation_safe_fraction": 0.0,
        }

    high = sum(
        orf["heuristic_score"]["classification"] == "High"
        for orf in orfs
    )
    moderate = sum(
        orf["heuristic_score"]["classification"] == "Moderate"
        for orf in orfs
    )
    low = len(orfs) - high - moderate

    canonical_complete = sum(
        orf["status"] == "Complete"
        and orf["start_type"] == "Canonical"
        for orf in orfs
    )
    alternative_complete = sum(
        orf["status"] == "Complete"
        and orf["start_type"] == "Alternative"
        for orf in orfs
    )

    translation_safe_fraction = (
        sum(orf["translation_safe"] for orf in orfs) / len(orfs)
    )

    return {
        "high_heuristic_score": high,
        "moderate_heuristic_score": moderate,
        "low_heuristic_score": low,
        "canonical_complete": canonical_complete,
        "alternative_complete": alternative_complete,
        "translation_safe_fraction": translation_safe_fraction,
    }


# ============================================================
# ENHANCED UNIFIED ORF ANALYSIS
# ============================================================

def analyze_orfs_enhanced(
    dna,
    include_partial=True,
    min_length=DEFAULT_MIN_ORF_LENGTH,
    top_n=5,
    start_codon_policy=None,
    nested=DEFAULT_REPORT_NESTED,
    ranking=DEFAULT_RANKING,
    genetic_code=1,
    topology="linear",
    include_5prime_partial=False,
):
    """
    Full ORF analysis with ranking, candidate comparison, overlap/nesting
    detection and downstream translation of a selected candidate.
    """
    result = analyze_orfs(
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

    ranked = rank_orf_candidates(result["orfs"], ranking=ranking)
    relationships = analyze_orf_relationships(ranked)

    result["orfs"] = ranked
    result["top_candidates"] = compare_orf_candidates(
        ranked,
        top_n=top_n,
        ranking=ranking,
    )
    result["relationships"] = relationships
    result["feature_summary"] = calculate_orf_feature_summary(ranked)
    result["candidate_count"] = len(ranked)

    if ranked:
        result["top_candidate"] = ranked[0]

    return result
