"""Codon usage, translation-initiation context, start/stop landscape and frame translation."""

from __future__ import annotations

from analyses._rna.bases import rna_codon_table, RNA_STANDARD_BASES


# ============================================================
# CODON STATISTICS
# ============================================================

def analyze_rna_codon_usage(rna):

    codons = []
    ambiguous_codons = []

    for position in range(
        0,
        len(rna) - 2,
        3
    ):

        codon = rna[
            position:
            position + 3
        ]

        if any(
            base not in RNA_STANDARD_BASES
            for base in codon
        ):

            ambiguous_codons.append({
                "codon":
                    codon,

                "position":
                    position + 1
            })

            continue

        codons.append(codon)

    counts = {}

    for codon in codons:

        counts[codon] = (
            counts.get(
                codon,
                0
            ) + 1
        )

    gc3_values = [
        1
        if codon[2] in {"G", "C"}
        else 0
        for codon in codons
    ]

    gc3_content = (
        sum(gc3_values)
        / len(gc3_values)
        * 100
        if gc3_values
        else None
    )

    return {
        "frame":
            0,

        "total_codons":
            len(codons),

        "unique_codons":
            len(counts),

        "codon_counts":
            counts,

        "ambiguous_codons":
            ambiguous_codons,

        "ambiguous_codon_count":
            len(ambiguous_codons),

        "gc3_content":
            gc3_content
    }


def calculate_rna_codon_diversity(rna):

    usage = analyze_rna_codon_usage(
        rna
    )

    if usage["total_codons"] == 0:
        return 0.0

    return (
        usage["unique_codons"]
        / usage["total_codons"]
    )


# ============================================================
# START CODON / INITIATION CONTEXT
# ============================================================

def analyze_rna_kozak_context(
    rna,
    start_position=0
):

    if (
        start_position < 6
        or start_position + 3 >= len(rna)
    ):

        return {
            "available": False,
            "context": None,
            "score": None,
            "classification":
                "Unavailable"
        }

    context = rna[
        start_position - 6:
        start_position + 4
    ]

    score = 0

    if rna[
        start_position - 6
    ] == "G":

        score += 1

    if rna[
        start_position - 3
    ] in {"A", "G"}:

        score += 1

    if rna[
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


def analyze_rna_shine_dalgarno(
    rna,
    start_position=0
):

    upstream_start = max(
        0,
        start_position - 20
    )

    upstream_end = max(
        0,
        start_position - 4
    )

    upstream = rna[
        upstream_start:
        upstream_end
    ]

    motifs = [
        "AGGAGG",
        "AAGGAG",
        "GGAGG",
        "AGGAG",
        "GAGG"
    ]

    best_motif = None
    best_position = None

    for motif in motifs:

        position = upstream.find(
            motif
        )

        if position == -1:
            continue

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
            "distance_to_start":
                None
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
        "position":
            best_position + 1,
        "distance_to_start":
            distance
    }


# ============================================================
# START / STOP CODON LANDSCAPE
# ============================================================

def find_rna_start_stop_landscape(rna):
    """
    Locate all definite AUG start codons and UAA/UAG/UGA
    stop codons in the RNA sequence.

    This is a sequence-level landscape and does not infer
    biological translation.
    """

    starts = []
    stops = []

    for position in range(len(rna) - 2):

        codon = rna[position:position + 3]

        if codon == "AUG":
            starts.append({
                "codon": codon,
                "position": position,
                "position_1_based": position + 1
            })

        elif codon in {"UAA", "UAG", "UGA"}:
            stops.append({
                "codon": codon,
                "position": position,
                "position_1_based": position + 1
            })

    return {
        "starts": starts,
        "start_count": len(starts),
        "stops": stops,
        "stop_count": len(stops)
    }


# ============================================================
# THREE-FRAME TRANSLATION
# ============================================================

def translate_rna_frame(rna, frame=0, genetic_code=1):
    """
    Translate one RNA reading frame.

    frame:
        0 = +1
        1 = +2
        2 = +3

    Ambiguous codons are represented as '?'.
    Translation stops at the first definite stop codon.

    genetic_code selects the NCBI table (audit P1-9).
    """

    table = rna_codon_table(genetic_code)

    if frame not in (0, 1, 2):
        raise ValueError(
            "RNA frame must be 0, 1 or 2."
        )

    protein = []
    codons = []
    stop_position = None

    for position in range(
        frame,
        len(rna) - 2,
        3
    ):

        codon = rna[
            position:position + 3
        ]

        codons.append(codon)

        amino_acid = table.get(
            codon
        )

        if amino_acid == "STOP":

            stop_position = position + 1
            break

        if amino_acid is None:
            protein.append("?")
        else:
            protein.append(amino_acid)

    return {
        "frame": frame + 1,
        "frame_offset": frame,
        "genetic_code": 1 if genetic_code in (None, 1) else int(genetic_code),
        "codons": codons,
        "codon_count": len(codons),
        "protein": "".join(protein),
        "protein_length": len(protein),
        "stop_position": stop_position,
        "complete": stop_position is not None
    }


def analyze_rna_three_frames(rna, genetic_code=1):
    """Return translation previews for all three forward RNA frames.

    All three frames use *genetic_code* (audit P1-9).
    """

    frames = [
        translate_rna_frame(rna, frame, genetic_code=genetic_code)
        for frame in range(3)
    ]

    return {
        "genetic_code": 1 if genetic_code in (None, 1) else int(genetic_code),
        "frames": frames,
        "best_frame": max(
            frames,
            key=lambda item: (
                item["complete"],
                item["protein_length"]
            )
        ) if frames else None
    }
