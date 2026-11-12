import math
from core.seqtools import codon_dict


# ============================================================
# RNA ANALYSIS
# ============================================================


# ============================================================
# CONSTANTS
# ============================================================

RNA_STANDARD_BASES = {
    "A",
    "U",
    "C",
    "G"
}


RNA_IUPAC_BASES = {
    "A": {"A"},
    "U": {"U"},
    "C": {"C"},
    "G": {"G"},

    "R": {"A", "G"},
    "Y": {"C", "U"},
    "S": {"G", "C"},
    "W": {"A", "U"},
    "K": {"G", "U"},
    "M": {"A", "C"},

    "B": {"C", "G", "U"},
    "D": {"A", "G", "U"},
    "H": {"A", "C", "U"},
    "V": {"A", "C", "G"},
    "N": {"A", "U", "C", "G"}
}


RNA_AMBIGUOUS_CODES = {
    "R": "A/G",
    "Y": "C/U",
    "S": "G/C",
    "W": "A/U",
    "K": "G/U",
    "M": "A/C",
    "B": "C/G/U",
    "D": "A/G/U",
    "H": "A/C/U",
    "V": "A/C/G",
    "N": "Any"
}


RNA_COMPLEMENT = {
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C",

    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",

    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N"
}


RNA_NUCLEOTIDE_RESIDUE_MASSES = {
    "A": 329.21,
    "U": 306.17,
    "C": 305.18,
    "G": 345.21
}


# Standard code (NCBI table 1) from core.seqtools, RNA alphabet.
RNA_CODON_TABLE = codon_dict(1, rna=True, stop_symbol="STOP")


def rna_codon_table(genetic_code=1):
    """Codon table (RNA alphabet, "STOP" for stops) for an NCBI genetic code.

    Audit P1-9: the three-frame codon view used to read table 1 even when a
    selected ORF had been found with another table, so the preview could
    disagree with the ORF's own translation.
    """
    if genetic_code in (None, 1):
        return RNA_CODON_TABLE
    return codon_dict(int(genetic_code), rna=True, stop_symbol="STOP")


def orf_genetic_code(selected_orf, default=1):
    """The genetic code of a selected ORF, or *default* when it carries none."""
    if selected_orf is None:
        return default
    value = _orf_value(selected_orf, "genetic_code", None)
    try:
        return int(value) if value else default
    except (TypeError, ValueError):
        return default


# ============================================================
# TRANSCRIPTION
# ============================================================

def transcribe_rna(dna):
    """Transcribe a coding-strand DNA sequence to RNA (T -> U).

    Whitespace is removed and letters are uppercased; a sequence that is
    already RNA is returned unchanged. Symbols outside the IUPAC nucleotide
    alphabet raise ValueError (audit D14-F9).
    """

    from core.alphabets import DNA_IUPAC, RNA_IUPAC, normalize_sequence, symbols_outside

    sequence = normalize_sequence(dna)
    invalid = symbols_outside(sequence, DNA_IUPAC | RNA_IUPAC)
    if invalid:
        shown = ", ".join(f"{symbol!r} at {position}" for position, symbol in invalid[:5])
        raise ValueError(
            f"Transcription needs a nucleotide sequence; found {shown}"
            + (" ..." if len(invalid) > 5 else "") + "."
        )

    return sequence.replace(
        "T",
        "U"
    )


# ============================================================
# BASIC BASE ANALYSIS
# ============================================================

def count_rna_bases(rna):

    counts = {
        "A": 0,
        "U": 0,
        "C": 0,
        "G": 0
    }

    for base in rna:

        if base in counts:
            counts[base] += 1

    return counts


def calculate_rna_base_percentages(
    base_counts
):

    total = sum(
        base_counts.values()
    )

    if total == 0:

        return {
            "A": 0.0,
            "U": 0.0,
            "C": 0.0,
            "G": 0.0
        }

    return {
        base:
            (
                count / total
            ) * 100
        for base, count
        in base_counts.items()
    }


# ============================================================
# GC / AU CONTENT
# ============================================================

def calculate_rna_gc_content(rna):

    counts = count_rna_bases(rna)

    total = sum(
        counts.values()
    )

    if total == 0:
        return 0.0

    return (
        (
            counts["G"]
            + counts["C"]
        )
        / total
    ) * 100


def calculate_rna_au_content(rna):

    counts = count_rna_bases(rna)

    total = sum(
        counts.values()
    )

    if total == 0:
        return 0.0

    return (
        (
            counts["A"]
            + counts["U"]
        )
        / total
    ) * 100


def calculate_rna_content_ranges(rna):

    valid_positions = 0

    gc_min = 0
    gc_max = 0

    au_min = 0
    au_max = 0

    ambiguous = False

    for symbol in rna:

        possible = RNA_IUPAC_BASES.get(
            symbol
        )

        if possible is None:
            continue

        valid_positions += 1

        if len(possible) > 1:
            ambiguous = True

        gc_values = [
            base in {"G", "C"}
            for base in possible
        ]

        au_values = [
            base in {"A", "U"}
            for base in possible
        ]

        gc_min += int(
            min(gc_values)
        )

        gc_max += int(
            max(gc_values)
        )

        au_min += int(
            min(au_values)
        )

        au_max += int(
            max(au_values)
        )

    if valid_positions == 0:

        return {
            "gc_min": 0.0,
            "gc_max": 0.0,
            "au_min": 0.0,
            "au_max": 0.0,
            "ambiguous": False,
            "evaluated_length": 0
        }

    return {
        "gc_min":
            gc_min
            / valid_positions
            * 100,

        "gc_max":
            gc_max
            / valid_positions
            * 100,

        "au_min":
            au_min
            / valid_positions
            * 100,

        "au_max":
            au_max
            / valid_positions
            * 100,

        "ambiguous":
            ambiguous,

        "evaluated_length":
            valid_positions
    }


# ============================================================
# SKEW
# ============================================================

def calculate_rna_gc_skew(rna):

    counts = count_rna_bases(rna)

    denominator = (
        counts["G"]
        + counts["C"]
    )

    if denominator == 0:
        return 0.0

    return (
        counts["G"]
        - counts["C"]
    ) / denominator


def calculate_rna_au_skew(rna):

    counts = count_rna_bases(rna)

    denominator = (
        counts["A"]
        + counts["U"]
    )

    if denominator == 0:
        return 0.0

    return (
        counts["A"]
        - counts["U"]
    ) / denominator


# ============================================================
# SLIDING WINDOW GC PROFILE
# ============================================================

def calculate_rna_gc_profile(
    rna,
    window_size=100,
    step=1
):

    profile = []

    if (
        window_size <= 0
        or step <= 0
        or len(rna) < window_size
    ):
        return profile

    for start in range(
        0,
        len(rna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        window = rna[
            start:end
        ]

        if all(
            base in RNA_STANDARD_BASES
            for base in window
        ):

            gc = (
                (
                    window.count("G")
                    + window.count("C")
                )
                / window_size
            ) * 100

            profile.append({
                "start":
                    start,

                "end":
                    end - 1,

                "start_1_based":
                    start + 1,

                "end_1_based":
                    end,

                "gc_content":
                    gc,

                "resolved":
                    True
            })

        else:

            profile.append({
                "start":
                    start,

                "end":
                    end - 1,

                "start_1_based":
                    start + 1,

                "end_1_based":
                    end,

                "gc_content":
                    None,

                "resolved":
                    False
            })

    return profile


# ============================================================
# AMBIGUITY / QUALITY
# ============================================================

def analyze_rna_ambiguity(rna):

    ambiguous_counts = {}
    ambiguous_positions = []

    for position, base in enumerate(rna):

        if base in RNA_AMBIGUOUS_CODES:

            ambiguous_positions.append({
                "position":
                    position,

                "position_1_based":
                    position + 1,

                "base":
                    base,

                "possibilities":
                    sorted(
                        RNA_IUPAC_BASES[base]
                    )
            })

            ambiguous_counts[base] = (
                ambiguous_counts.get(
                    base,
                    0
                ) + 1
            )

    count = len(
        ambiguous_positions
    )

    percentage = (
        count / len(rna) * 100
        if rna
        else 0.0
    )

    return {
        "ambiguous_count":
            count,

        "ambiguous_percentage":
            percentage,

        "ambiguous_counts":
            ambiguous_counts,

        "ambiguous_positions":
            ambiguous_positions,

        "unambiguous":
            count == 0
    }


def analyze_rna_quality(rna):

    standard = 0
    ambiguous = 0
    invalid = 0

    invalid_positions = []

    for position, symbol in enumerate(
        rna
    ):

        if symbol in RNA_STANDARD_BASES:

            standard += 1

        elif symbol in RNA_AMBIGUOUS_CODES:

            ambiguous += 1

        else:

            invalid += 1

            invalid_positions.append({
                "position":
                    position,

                "position_1_based":
                    position + 1,

                "symbol":
                    symbol
            })

    length = len(rna)

    return {
        "standard_count":
            standard,

        "ambiguous_count":
            ambiguous,

        "invalid_count":
            invalid,

        "standard_percentage":
            standard / length * 100
            if length else 0.0,

        "ambiguous_percentage":
            ambiguous / length * 100
            if length else 0.0,

        "invalid_percentage":
            invalid / length * 100
            if length else 0.0,

        "invalid_positions":
            invalid_positions,

        "valid":
            invalid == 0,

        "unambiguous":
            (
                ambiguous == 0
                and invalid == 0
            )
    }


# ============================================================
# MOLECULAR WEIGHT
# ============================================================

def calculate_rna_molecular_weight(rna):

    minimum_weight = 0.0
    maximum_weight = 0.0

    ambiguous = False
    evaluated_positions = 0

    for symbol in rna:

        possible = RNA_IUPAC_BASES.get(
            symbol
        )

        if possible is None:
            continue

        evaluated_positions += 1

        if len(possible) > 1:
            ambiguous = True

        weights = [
            RNA_NUCLEOTIDE_RESIDUE_MASSES[
                base
            ]
            for base in possible
        ]

        minimum_weight += min(
            weights
        )

        maximum_weight += max(
            weights
        )

    return {
        "min":
            minimum_weight,

        "max":
            maximum_weight,

        "ambiguous":
            ambiguous,

        "evaluated_positions":
            evaluated_positions
    }


# ============================================================
# COMPLEXITY
# ============================================================

def calculate_rna_complexity(rna):

    counts = count_rna_bases(rna)

    total = sum(
        counts.values()
    )

    if total == 0:
        return 0.0

    entropy = 0.0

    for count in counts.values():

        if count == 0:
            continue

        probability = (
            count / total
        )

        entropy -= (
            probability
            * math.log2(probability)
        )

    return entropy / 2.0


# ============================================================
# COMPLEXITY PROFILE
# ============================================================

def calculate_rna_complexity_profile(
    rna,
    window_size=100,
    step=1
):

    profile = []

    if (
        window_size <= 0
        or step <= 0
        or len(rna) < window_size
    ):
        return profile

    for start in range(
        0,
        len(rna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        window = rna[
            start:end
        ]

        profile.append({
            "start":
                start,

            "end":
                end - 1,

            "start_1_based":
                start + 1,

            "end_1_based":
                end,

            "complexity":
                calculate_rna_complexity(
                    window
                )
        })

    return profile


# ============================================================
# HOMOPOLYMERS
# ============================================================

def find_rna_homopolymers(
    rna,
    minimum_length=4
):

    regions = []

    if not rna:
        return regions

    start = 0

    for index in range(
        1,
        len(rna) + 1
    ):

        if (
            index < len(rna)
            and rna[index] == rna[start]
            and rna[start]
            in RNA_STANDARD_BASES
        ):
            continue

        length = index - start

        if (
            length >= minimum_length
            and rna[start]
            in RNA_STANDARD_BASES
        ):

            regions.append({
                "base":
                    rna[start],

                "length":
                    length,

                "start":
                    start + 1,

                "end":
                    index
            })

        start = index

    return regions


# ============================================================
# POLY-A / POLY-U
# ============================================================

def find_rna_homopolymer_regions(
    rna,
    base,
    minimum_length=4
):

    regions = []

    start = None

    for index, symbol in enumerate(
        rna
    ):

        if symbol == base:

            if start is None:
                start = index

        else:

            if start is not None:

                length = index - start

                if length >= minimum_length:

                    regions.append({
                        "base":
                            base,

                        "start":
                            start + 1,

                        "end":
                            index,

                        "length":
                            length
                    })

                start = None

    if start is not None:

        length = len(rna) - start

        if length >= minimum_length:

            regions.append({
                "base":
                    base,

                "start":
                    start + 1,

                "end":
                    len(rna),

                "length":
                    length
            })

    return regions


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

    if best_motif is None:

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
# RNA COMPLEMENT
# ============================================================

def get_rna_complement(rna):

    return "".join(
        RNA_COMPLEMENT.get(
            base,
            base
        )
        for base in rna
    )


def get_rna_reverse_complement(rna):

    return get_rna_complement(
        rna
    )[::-1]


# ============================================================
# RNA HAIRPIN SCREENING
# ============================================================

def find_rna_hairpins(
    rna,
    min_stem=4,
    min_loop=3,
    max_loop=12,
    allow_gu=False,
):
    """Candidate stem-loops: stems of min_stem..min_stem+3 nt around a loop.

    A candidate is reported when at least 75 % of the stem positions with
    standard bases pair. Watson-Crick pairs (A-U, G-C) count; with
    allow_gu=True the G-U wobble pair, a regular but weaker pair in RNA
    helices, counts as well and "wobble_pairs" says how many a stem uses.
    The default stays Watson-Crick only: without free energies, counting
    G-U like G-C multiplies the candidates (58 -> 260 on the golden-master
    reference) instead of finding stronger stems, so the structure summary
    reports the G-U count separately. A sequence screen, not a prediction.
    """

    structures = []

    sequence = rna.upper()

    if (
        min_stem < 1
        or min_loop < 1
        or max_loop < min_loop
    ):
        return structures

    for stem_start in range(
        len(sequence)
    ):

        for stem_length in range(
            min_stem,
            min_stem + 4
        ):

            left_end = (
                stem_start
                + stem_length
            )

            if left_end >= len(sequence):
                continue

            for loop_length in range(
                min_loop,
                max_loop + 1
            ):

                right_start = (
                    left_end
                    + loop_length
                )

                right_end = (
                    right_start
                    + stem_length
                )

                if right_end > len(sequence):
                    continue

                left = sequence[
                    stem_start:left_end
                ]

                right = sequence[
                    right_start:right_end
                ]

                paired = 0
                wobble = 0
                possible_pairs = 0

                for left_base, right_base in zip(
                    left,
                    reversed(right)
                ):

                    if (
                        left_base
                        not in RNA_STANDARD_BASES
                        or right_base
                        not in RNA_STANDARD_BASES
                    ):
                        continue

                    possible_pairs += 1

                    if (
                        (
                            left_base == "A"
                            and right_base == "U"
                        )
                        or
                        (
                            left_base == "U"
                            and right_base == "A"
                        )
                        or
                        (
                            left_base == "C"
                            and right_base == "G"
                        )
                        or
                        (
                            left_base == "G"
                            and right_base == "C"
                        )
                    ):

                        paired += 1

                    elif allow_gu and {left_base, right_base} == {"G", "U"}:

                        paired += 1
                        wobble += 1

                if possible_pairs == 0:
                    continue

                pairing_fraction = (
                    paired
                    / possible_pairs
                )

                if pairing_fraction >= 0.75:

                    structures.append({
                        "start":
                            stem_start + 1,

                        "end":
                            right_end,

                        "stem_length":
                            stem_length,

                        "loop_length":
                            loop_length,

                        "pairing_fraction":
                            pairing_fraction,

                        "wobble_pairs":
                            wobble,

                        "left_stem":
                            left,

                        "right_stem":
                            right
                    })

    unique = []
    seen = set()

    for structure in structures:

        key = (
            structure["start"],
            structure["end"],
            structure["stem_length"],
            structure["loop_length"]
        )

        if key in seen:
            continue

        seen.add(key)
        unique.append(
            structure
        )

    unique.sort(
        key=lambda item: (
            item["pairing_fraction"],
            item["stem_length"]
        ),
        reverse=True
    )

    return unique


# ============================================================
# RNA STRUCTURAL SUMMARY
# ============================================================

def summarize_rna_structure(rna):

    hairpins = find_rna_hairpins(
        rna
    )

    homopolymers = (
        find_rna_homopolymers(
            rna
        )
    )

    return {
        "hairpins":
            hairpins,

        "hairpin_count":
            len(hairpins),

        # Candidates when G-U wobble pairs count as pairs (P2-5); reported
        # beside the Watson-Crick screen, not instead of it.
        "hairpin_count_with_gu":
            len(find_rna_hairpins(rna, allow_gu=True)),

        "homopolymers":
            homopolymers,

        "homopolymer_count":
            len(homopolymers)
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


# ============================================================
# RNA REPORT
# ============================================================


# ============================================================
# SYNAPSE PIPELINE INTEGRATION
# ============================================================
#
# This section is the ONLY pipeline/UI layer in this module.
#
# Input resolution:
#   1. Explicitly supplied selected ORF
#   2. ORF stored in the pipeline context
#   3. Complete DNA sequence
#
# The original RNA analysis engine above is kept intact.
# The functions below only connect that engine to Synapse.
# ============================================================


def _orf_value(orf, key, default=None):
    """Read a value from either a dict-like ORF or an ORF object."""

    if orf is None:
        return default

    if isinstance(orf, dict):
        return orf.get(key, default)

    return getattr(orf, key, default)


def _clean_sequence(sequence):
    """Normalize a sequence without changing IUPAC symbols."""

    return (
        str(sequence or "")
        .upper()
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .replace("\t", "")
    )


def _get_pipeline_context(pipeline):
    """Return the pipeline context regardless of pipeline implementation."""

    if pipeline is None:
        return None

    return getattr(pipeline, "context", pipeline)


def _get_pipeline_selected_orf(pipeline):
    """Get the selected ORF from the pipeline if one exists."""

    context = _get_pipeline_context(pipeline)

    if context is None:
        return None

    selected = getattr(context, "selected_orf", None)

    if selected is None:
        selected = getattr(pipeline, "selected_orf", None)

    return selected


def resolve_rna_source(dna, selected_orf=None, pipeline=None):
    """
    Resolve the biological input for RNA analysis.

    Priority:
        1. Explicit selected_orf
        2. selected_orf stored in pipeline
        3. Complete DNA sequence

    This is important because users are allowed to enter RNA Analysis
    directly without selecting an ORF first.
    """

    full_dna = _clean_sequence(dna)

    if selected_orf is None:
        selected_orf = _get_pipeline_selected_orf(pipeline)

    if selected_orf is not None:
        orf_sequence = _clean_sequence(
            _orf_value(selected_orf, "sequence", "")
        )

        if orf_sequence:
            return {
                "dna": orf_sequence,
                "source": "Selected ORF",
                "mode": "orf",
                "orf": selected_orf,
            }

    return {
        "dna": full_dna,
        "source": "Full DNA Sequence",
        "mode": "full_sequence",
        "orf": None,
    }


def _store_rna_in_pipeline(analysis, pipeline):
    """
    Store RNA analysis in whichever pipeline storage API is available.

    This function is deliberately defensive so that RNA works with the
    current pipeline implementation without forcing a specific class design.
    """

    if pipeline is None:
        return

    context = _get_pipeline_context(pipeline)

    if context is None:
        return

    store_method = getattr(
        context,
        "store_rna_analysis",
        None
    )

    if callable(store_method):
        try:
            store_method(
                analysis,
                rna_sequence=analysis.get("rna")
            )
        except TypeError:
            store_method(analysis)

    # Always keep a direct reference when the context allows attributes.
    try:
        context.rna_analysis = analysis
    except Exception:
        pass


def analyze_rna_from_pipeline(
    dna,
    pipeline=None,
    selected_orf=None,
    **analysis_kwargs
):
    """
    Run the original analyze_rna() engine with pipeline-aware input.

    'header' is accepted and ignored for compatibility with main.py.
    It must NOT be forwarded into analyze_rna(), because analyze_rna()
    does not define a header parameter.
    """

    # Compatibility with callers that still pass header=...
    analysis_kwargs.pop("header", None)

    source = resolve_rna_source(
        dna,
        selected_orf=selected_orf,
        pipeline=pipeline,
    )

    if not source["dna"]:
        raise ValueError(
            "No DNA sequence is available for RNA analysis."
        )

    # Audit P1-9: a selected ORF carries the table it was found with, so the
    # codon view must use it unless the caller states another one.
    analysis_kwargs.setdefault(
        "genetic_code",
        orf_genetic_code(source["orf"]),
    )

    analysis = analyze_rna(
        source["dna"],
        **analysis_kwargs
    )

    analysis["pipeline_input"] = {
        "source": source["source"],
        "mode": source["mode"],
        "selected_orf": source["orf"],
        "dna": source["dna"],
        "dna_length": len(source["dna"]),
    }

    _store_rna_in_pipeline(
        analysis,
        pipeline
    )

    return analysis


# ============================================================
# PUBLIC PIPELINE ALIASES
# ============================================================
#
# Keep these aliases so existing main.py / pipeline.py imports
# continue to work. There is only ONE actual UI implementation.
# ============================================================

run_rna_pipeline = analyze_rna_from_pipeline