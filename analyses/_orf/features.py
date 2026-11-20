"""Per-ORF evaluation: ambiguity, codon usage, start context, CDS quality, heuristic score."""

from __future__ import annotations

from collections import Counter

from analyses._orf.codes import (
    expand_iupac_codon,
    GENETIC_CODE,
    _is_alternative_start,
    _stop_codons_of,
)


# ============================================================
# ORF AMBIGUITY ANALYSIS
# ============================================================

_UNAMBIGUOUS = frozenset("ATCG")


def _empty_ambiguity():
    return {
        "ambiguous_codon_count": 0,
        "ambiguous_codons": [],
        "deterministic_ambiguous_count": 0,
        "deterministic_ambiguous_codons": [],
        "uncertain_codon_count": 0,
        "uncertain_codons": [],
        "possible_stop_codon_count": 0,
        "possible_stop_codons": [],
        "translation_safe": True,
    }


def analyze_orf_ambiguity(orf_sequence, codon_table=None):
    codon_table = GENETIC_CODE if codon_table is None else codon_table

    # A sequence of A, C, G and T only has no ambiguous codon: the loop
    # below would find none, so the empty result is returned directly (a
    # genome-wide search evaluates about 100,000 ORF candidates per Mb).
    if not set(orf_sequence) - _UNAMBIGUOUS:
        return _empty_ambiguity()

    ambiguous_codons = []
    deterministic_ambiguous_codons = []
    uncertain_codons = []
    possible_stop_codons = []

    for position in range(
        0,
        len(orf_sequence),
        3
    ):

        codon = orf_sequence[
            position:
            position + 3
        ]

        if len(codon) < 3:
            continue

        is_ambiguous = any(
            base not in "ATCG"
            for base in codon
        )

        if not is_ambiguous:
            continue

        expanded_codons = (
            expand_iupac_codon(
                codon
            )
        )

        translations = set()

        for expanded_codon in (
            expanded_codons
        ):

            amino_acid = codon_table.get(
                expanded_codon
            )

            if amino_acid is not None:
                translations.add(
                    amino_acid
                )

        amino_acids = sorted(
            amino_acid
            for amino_acid in translations
            if amino_acid != "STOP"
        )

        contains_stop = (
            "STOP" in translations
        )

        if len(translations) == 1:

            translation_status = (
                "Deterministic"
            )

            resolved_translation = (
                next(iter(translations))
            )

            deterministic_ambiguous_codons.append({
                "codon": codon,
                "position": position + 1,
                "translation":
                    resolved_translation
            })

        else:

            translation_status = (
                "Uncertain"
            )

            resolved_translation = None

            uncertain_codons.append({
                "codon": codon,
                "position": position + 1,
                "possibilities":
                    amino_acids,
                "contains_stop":
                    contains_stop
            })

        if contains_stop:

            possible_stop_codons.append({
                "codon": codon,
                "position": position + 1
            })

        ambiguous_codons.append({
            "codon": codon,
            "position": position + 1,
            "translation_status":
                translation_status,
            "translation":
                resolved_translation,
            "possibilities":
                amino_acids,
            "contains_stop":
                contains_stop
        })

    return {
        "ambiguous_codon_count":
            len(ambiguous_codons),

        "ambiguous_codons":
            ambiguous_codons,

        "deterministic_ambiguous_count":
            len(
                deterministic_ambiguous_codons
            ),

        "deterministic_ambiguous_codons":
            deterministic_ambiguous_codons,

        "uncertain_codon_count":
            len(uncertain_codons),

        "uncertain_codons":
            uncertain_codons,

        "possible_stop_codon_count":
            len(possible_stop_codons),

        "possible_stop_codons":
            possible_stop_codons,

        "translation_safe":
            len(uncertain_codons) == 0
    }


# ============================================================
# CODON STATISTICS
# ============================================================

def analyze_codon_usage(coding_sequence, codon_table=None):
    """Codon counts and GC3 of the codons before the first stop codon.

    Rewritten on 2026-11-20 for speed with the same result: the codons
    are the full triplets from position 0, cut at the first codon that is
    a stop of the table; GC3 counts the codons of A, C, G and T only (the
    integer counts and the division are unchanged). Checked against the
    former loop in tests/test_orf_feature_speed.py.
    """
    stop_codons = _stop_codons_of(codon_table)

    codons = [
        coding_sequence[position:position + 3]
        for position in range(0, len(coding_sequence) - 2, 3)
    ]

    for index, codon in enumerate(codons):
        if codon in stop_codons:
            del codons[index:]
            break

    counts = Counter(codons)

    standard_codon_count = 0
    standard_total = 0
    gc3_total = 0

    for codon, count in counts.items():
        if not set(codon) - _UNAMBIGUOUS:
            standard_codon_count += 1
            standard_total += count
            if codon[2] in "GC":
                gc3_total += count

    if standard_total:

        gc3_content = (
            gc3_total
            / standard_total
            * 100
        )

    else:

        gc3_content = None

    return {
        "total_codons":
            len(codons),

        "unique_codons":
            len(counts),

        "codon_counts":
            dict(counts),

        "standard_codon_count":
            standard_codon_count,

        "gc3_content":
            gc3_content
    }


def codon_diversity_from_usage(usage):
    """calculate_codon_diversity from an analyze_codon_usage result."""
    total = usage["total_codons"]
    if total == 0:
        return 0.0
    return usage["unique_codons"] / total


# ============================================================
# CODON DIVERSITY
# ============================================================

def calculate_codon_diversity(coding_sequence, codon_table=None):

    usage = analyze_codon_usage(coding_sequence, codon_table)

    total = usage["total_codons"]
    unique = usage["unique_codons"]

    if total == 0:
        return 0.0

    return unique / total


# ============================================================
# START CONTEXT — KOZAK
# ============================================================

def analyze_kozak_context(
    sequence,
    start_position
):

    """
    Simplified Kozak-context analysis.

    Consensus:

        GCCRCCATGG

    Evaluated positions:
        -6 = G
        -3 = A/G
        +4 = G

    This is a sequence-context heuristic and does not prove
    translation initiation.
    """

    if (
        start_position < 6
        or start_position + 3 >= len(sequence)
    ):

        return {
            "available": False,
            "context": None,
            "score": None,
            "classification":
                "Unavailable"
        }

    context_start = (
        start_position - 6
    )

    context_end = (
        start_position + 4
    )

    context = sequence[
        context_start:
        context_end
    ]

    score = 0

    if sequence[
        start_position - 6
    ] == "G":

        score += 1

    if sequence[
        start_position - 3
    ] in "AG":

        score += 1

    if sequence[
        start_position + 3
    ] == "G":

        score += 1

    if score == 3:
        classification = "Strong"

    elif score == 2:
        classification = "Moderate"

    elif score == 1:
        classification = "Weak"

    else:
        classification = "Poor"

    return {
        "available": True,
        "context": context,
        "score": score,
        "classification":
            classification
    }


# ============================================================
# START CONTEXT — SHINE-DALGARNO
# ============================================================

def analyze_shine_dalgarno(
    sequence,
    start_position
):

    upstream_start = max(
        0,
        start_position - 20
    )

    upstream_end = max(
        0,
        start_position - 4
    )

    upstream = sequence[
        upstream_start:
        upstream_end
    ]

    motifs = [
        "AGGAGG",
        "GGAGG",
        "AAGGAG",
        "AGGAG",
        "GAGG"
    ]

    best_motif = None
    best_position = None

    for motif in motifs:

        position = upstream.find(
            motif
        )

        if position != -1:

            if (
                best_motif is None
                or len(motif)
                > len(best_motif)
            ):

                best_motif = motif

                best_position = (
                    upstream_start
                    + position
                )

    # Motif and position are set together, so both are None or neither is.
    if best_motif is None or best_position is None:

        return {
            "available": True,
            "detected": False,
            "motif": None,
            "position": None,
            "distance_to_start": None
        }

    distance = (
        start_position
        - (
            best_position
            + len(best_motif)
        )
    )

    return {
        "available": True,
        "detected": True,
        "motif": best_motif,
        "position": best_position,
        "distance_to_start": distance
    }


# ============================================================
# START CONTEXT
# ============================================================

def analyze_start_context(
    sequence,
    start_position
):

    return {
        "kozak":
            analyze_kozak_context(
                sequence,
                start_position
            ),

        "shine_dalgarno":
            analyze_shine_dalgarno(
                sequence,
                start_position
            )
    }


# ============================================================
# CDS QUALITY
# ============================================================

def calculate_cds_quality(
    orf
):

    score = 0
    reasons = []

    if orf["status"] == "Complete":

        score += 30

        reasons.append(
            "Complete start-to-stop ORF"
        )

    else:

        score += 10

        reasons.append(
            "Partial ORF"
        )

    if orf["length"] >= 300:

        score += 20

        reasons.append(
            "Substantial ORF length"
        )

    elif orf["length"] >= 150:

        score += 15

    else:

        score += 5

    if orf["start_codon"] == "ATG":

        score += 20

        reasons.append(
            "Canonical ATG start"
        )

    elif _is_alternative_start(orf):

        score += 10

        reasons.append(
            "Alternative start codon"
        )

    if orf["translation_safe"]:

        score += 15

        reasons.append(
            "Translation is unambiguous"
        )

    else:

        score -= 10

        reasons.append(
            "Translation contains uncertainty"
        )

    if orf["ambiguous_codon_count"] == 0:

        score += 10

    elif (
        orf["ambiguous_codon_count"] <= 2
    ):

        score += 5

    else:

        score -= 5

    if (
        orf["possible_stop_codon_count"]
        > 0
    ):

        score -= 5

        reasons.append(
            "Possible ambiguous stop codons"
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
        "classification":
            classification,
        "reasons": reasons
    }


# ============================================================
# ORF HEURISTIC SCORE
# ============================================================

def calculate_orf_heuristic_score(
    orf
):

    score = 0

    if orf["status"] == "Complete":
        score += 25
    else:
        score += 10

    if orf["start_codon"] == "ATG":
        score += 20

    elif _is_alternative_start(orf):
        score += 10

    if orf["length"] >= 600:
        score += 20

    elif orf["length"] >= 300:
        score += 15

    elif orf["length"] >= 150:
        score += 10

    else:
        score += 5

    if orf["translation_safe"]:
        score += 15
    else:
        score -= 5

    if orf["ambiguous_codon_count"] == 0:
        score += 10

    elif (
        orf["ambiguous_codon_count"] <= 2
    ):
        score += 5

    else:
        score -= 5

    kozak = (
        orf["start_context"]
        ["kozak"]
    )

    if kozak["available"]:

        if kozak["classification"] == "Strong":
            score += 10

        elif kozak["classification"] == "Moderate":
            score += 6

        elif kozak["classification"] == "Weak":
            score += 3

    shine = (
        orf["start_context"]
        ["shine_dalgarno"]
    )

    if shine["detected"]:
        score += 5

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
        "classification":
            classification
    }


# Former name, kept so existing callers keep working (D18). The value is a
# heuristic ranking score, not a confidence or a probability.
calculate_orf_confidence = calculate_orf_heuristic_score
