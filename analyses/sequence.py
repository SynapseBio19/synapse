import math
from core.textformat import minor_rule


# ============================================================
# SEQUENCE ANALYSIS
# ============================================================


# ============================================================
# CONSTANTS
# ============================================================

DNA_STANDARD_BASES = {
    "A",
    "T",
    "C",
    "G"
}


DNA_AMBIGUOUS_CODES = {
    "R": "A/G",
    "Y": "C/T",
    "S": "G/C",
    "W": "A/T",
    "K": "G/T",
    "M": "A/C",
    "B": "C/G/T",
    "D": "A/G/T",
    "H": "A/C/T",
    "V": "A/C/G",
    "N": "Any"
}


DNA_IUPAC_BASES = {
    "A": {"A"},
    "T": {"T"},
    "C": {"C"},
    "G": {"G"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "T", "C", "G"}
}


DNA_COMPLEMENT = {
    "A": "T",
    "T": "A",
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


# Average nucleotide residue masses for ssDNA.
# Used compositionally together with complementary bases
# for an estimated dsDNA molecular weight.

DNA_NUCLEOTIDE_RESIDUE_MASSES = {
    "A": 313.21,
    "T": 304.20,
    "C": 289.18,
    "G": 329.21
}


# Approximate average molecular weight of one base pair.

AVERAGE_DSDNA_BP_MASS = 660.0


# ============================================================
# INTERNAL HELPERS
# ============================================================


# ============================================================
# INPUT NORMALIZATION
# ============================================================

def normalize_dna_sequence(dna):
    """
    Normalize sequence input for sequence-level analysis.

    Whitespace is ignored and letters are converted to uppercase.
    Symbols that are not DNA/IUPAC symbols are preserved so that
    the quality analysis can report them as invalid rather than
    silently deleting them.
    """

    if dna is None:
        return ""

    return "".join(
        str(dna).split()
    ).upper()


# ============================================================
# SEQUENCE STATISTICS
# ============================================================

def calculate_sequence_statistics(
    dna,
    nucleotide_counts=None,
    nucleotide_percentages=None
):
    """
    Compact top-level statistics for the Sequence Overview.

    All counts are based on standard A/T/C/G bases.
    Ambiguous and invalid symbols are reported separately.
    """

    if nucleotide_counts is None:
        nucleotide_counts = count_nucleotides(dna)

    if nucleotide_percentages is None:
        nucleotide_percentages = (
            calculate_nucleotide_percentages(
                nucleotide_counts
            )
        )

    standard_length = sum(
        nucleotide_counts.values()
    )

    gc_count = (
        nucleotide_counts["G"]
        + nucleotide_counts["C"]
    )

    at_count = (
        nucleotide_counts["A"]
        + nucleotide_counts["T"]
    )

    return {
        "length": len(dna),
        "standard_length": standard_length,
        "gc_count": gc_count,
        "at_count": at_count,
        "gc_content": (
            gc_count / standard_length * 100
            if standard_length
            else 0.0
        ),
        "at_content": (
            at_count / standard_length * 100
            if standard_length
            else 0.0
        ),
        "nucleotide_counts": dict(
            nucleotide_counts
        ),
        "nucleotide_percentages": dict(
            nucleotide_percentages
        )
    }


# ============================================================
# PROFILE AVAILABILITY
# ============================================================

def describe_profile_availability(
    dna,
    window_size,
    profile_name
):
    """
    Explain why a sliding-window profile is unavailable for a
    sequence shorter than the configured window.
    """

    if len(dna) >= window_size:
        return {
            "available": True,
            "reason": None,
            "message": None
        }

    return {
        "available": False,
        "reason": "SEQUENCE_SHORTER_THAN_WINDOW",
        "message": (
            f"{profile_name} is not available for this sequence. "
            f"Sequence length: {len(dna)} bp. "
            f"Minimum configured window: {window_size} bp."
        )
    }


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
# CpG CRITERIA
# ============================================================

DEFAULT_CPG_ISLAND_CRITERIA = {
    "window_size": 200,
    "min_gc_content": 50.0,
    "min_oe_ratio": 0.60,
    "definition": (
        "Classical Gardiner-Garden/Frommer-style screening "
        "criteria"
    )
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




def _is_standard_base(base):

    return base in DNA_STANDARD_BASES


def _is_iupac_base(base):

    return base in DNA_IUPAC_BASES


def _get_standard_bases(dna):

    return [
        base
        for base in dna
        if base in DNA_STANDARD_BASES
    ]


def _get_valid_iupac_bases(dna):

    return [
        base
        for base in dna
        if base in DNA_IUPAC_BASES
    ]


def _calculate_gc_for_window(window):

    if not window:
        return 0.0

    gc_count = (
        window.count("G")
        + window.count("C")
    )

    return (
        gc_count
        / len(window)
    ) * 100


# ============================================================
# NUCLEOTIDE COUNTS
# ============================================================


def count_nucleotides(dna):
    """Counts of A, T, C and G; U (the RNA form of T) is counted as T.

    Counting U keeps GC/AT content and skews correct for an RNA sequence
    passed to this DNA module directly (audit D14-F10).
    """

    nucleotide_counts = {
        "A": 0,
        "T": 0,
        "C": 0,
        "G": 0
    }

    for nucleotide in dna:

        if nucleotide == "U":
            nucleotide = "T"

        if nucleotide in nucleotide_counts:

            nucleotide_counts[
                nucleotide
            ] += 1

    return nucleotide_counts


# ============================================================
# AMBIGUOUS BASE COUNTS
# ============================================================


def count_ambiguous_bases(dna):

    counts = {}

    for code in DNA_AMBIGUOUS_CODES:

        count = dna.count(code)

        if count > 0:

            counts[code] = count

    return counts


# ============================================================
# COMPLEMENT
# ============================================================


def get_complement(dna):

    complement = []

    for base in dna:

        complement.append(
            DNA_COMPLEMENT.get(
                base,
                base
            )
        )

    return "".join(complement)


def get_reverse_complement(dna):

    return get_complement(dna)[::-1]


# ============================================================
# SEQUENCE QUALITY
# ============================================================


def analyze_sequence_quality(dna):

    standard_count = 0
    ambiguous_count = 0
    invalid_count = 0

    ambiguous_positions = []
    invalid_positions = []

    for position, nucleotide in enumerate(dna):

        if nucleotide in DNA_STANDARD_BASES:

            standard_count += 1

        elif nucleotide in DNA_AMBIGUOUS_CODES:

            ambiguous_count += 1

            ambiguous_positions.append({

                "position":
                    position,

                "position_1_based":
                    position + 1,

                "base":
                    nucleotide,

                "possibilities":
                    sorted(
                        DNA_IUPAC_BASES[
                            nucleotide
                        ]
                    )
            })

        else:

            invalid_count += 1

            invalid_positions.append({

                "position":
                    position,

                "position_1_based":
                    position + 1,

                "symbol":
                    nucleotide
            })

    sequence_length = len(dna)

    if sequence_length > 0:

        standard_percentage = (
            standard_count
            / sequence_length
        ) * 100

        ambiguous_percentage = (
            ambiguous_count
            / sequence_length
        ) * 100

        invalid_percentage = (
            invalid_count
            / sequence_length
        ) * 100

    else:

        standard_percentage = 0.0
        ambiguous_percentage = 0.0
        invalid_percentage = 0.0

    return {

        "standard_bases":
            standard_count,

        "ambiguous_bases":
            ambiguous_count,

        "invalid_symbols":
            invalid_count,

        "standard_percentage":
            standard_percentage,

        "ambiguous_percentage":
            ambiguous_percentage,

        "invalid_percentage":
            invalid_percentage,

        "ambiguous_positions":
            ambiguous_positions,

        "invalid_positions":
            invalid_positions,

        "valid":
            invalid_count == 0,

        "unambiguous":
            (
                ambiguous_count == 0
                and invalid_count == 0
            )
    }


# ============================================================
# NUCLEOTIDE PERCENTAGES
# ============================================================


def calculate_nucleotide_percentages(
    nucleotide_counts
):

    total = sum(
        nucleotide_counts.values()
    )

    if total == 0:

        return {
            "A": 0.0,
            "T": 0.0,
            "C": 0.0,
            "G": 0.0
        }

    percentages = {}

    for nucleotide, count in (
        nucleotide_counts.items()
    ):

        percentages[nucleotide] = (
            count / total
        ) * 100

    return percentages


# ============================================================
# GC / AT CONTENT
# ============================================================


def calculate_gc_content(dna):

    counts = count_nucleotides(dna)

    known_bases = sum(
        counts.values()
    )

    if known_bases == 0:

        return 0.0

    return (

        (
            counts["G"]
            + counts["C"]
        )

        / known_bases

    ) * 100


def calculate_at_content(dna):

    counts = count_nucleotides(dna)

    known_bases = sum(
        counts.values()
    )

    if known_bases == 0:

        return 0.0

    return (

        (
            counts["A"]
            + counts["T"]
        )

        / known_bases

    ) * 100


# ============================================================
# IUPAC-AWARE GC / AT RANGE
# ============================================================


def calculate_dna_content_ranges(dna):

    valid_symbols = _get_valid_iupac_bases(
        dna
    )

    evaluated_length = len(
        valid_symbols
    )

    if evaluated_length == 0:

        return {

            "gc_min": 0.0,
            "gc_max": 0.0,

            "at_min": 0.0,
            "at_max": 0.0,

            "ambiguous": False,

            "evaluated_length": 0
        }

    gc_min_count = 0
    gc_max_count = 0

    at_min_count = 0
    at_max_count = 0

    ambiguous = False

    for symbol in valid_symbols:

        possible_bases = (
            DNA_IUPAC_BASES[symbol]
        )

        if len(possible_bases) > 1:

            ambiguous = True

        possible_gc = [

            base in {"G", "C"}

            for base in possible_bases
        ]

        possible_at = [

            base in {"A", "T"}

            for base in possible_bases
        ]

        gc_min_count += int(
            min(possible_gc)
        )

        gc_max_count += int(
            max(possible_gc)
        )

        at_min_count += int(
            min(possible_at)
        )

        at_max_count += int(
            max(possible_at)
        )

    return {

        "gc_min":
            gc_min_count
            / evaluated_length
            * 100,

        "gc_max":
            gc_max_count
            / evaluated_length
            * 100,

        "at_min":
            at_min_count
            / evaluated_length
            * 100,

        "at_max":
            at_max_count
            / evaluated_length
            * 100,

        "ambiguous":
            ambiguous,

        "evaluated_length":
            evaluated_length
    }


# ============================================================
# NUCLEOTIDE SKEW
# ============================================================


def calculate_gc_skew(dna):

    g_count = dna.count("G")
    c_count = dna.count("C")

    total_gc = (
        g_count + c_count
    )

    if total_gc == 0:

        return 0.0

    return (
        g_count - c_count
    ) / total_gc


def calculate_at_skew(dna):

    a_count = dna.count("A")
    t_count = dna.count("T")

    total_at = (
        a_count + t_count
    )

    if total_at == 0:

        return 0.0

    return (
        a_count - t_count
    ) / total_at


# ============================================================
# SLIDING-WINDOW GC CONTENT
# ============================================================


def calculate_gc_content_profile(
    dna,
    window_size=100,
    step=1
):

    profile = []

    if window_size <= 0:
        return profile

    if step <= 0:
        return profile

    if len(dna) < window_size:
        return profile

    for start in range(
        0,
        len(dna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        window = dna[
            start:end
        ]

        if all(
            base in DNA_STANDARD_BASES
            for base in window
        ):

            gc_content = (
                _calculate_gc_for_window(
                    window
                )
            )

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
                    gc_content,

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
# SEQUENCE COMPLEXITY
# ============================================================


def calculate_sequence_complexity(dna):

    standard_bases = _get_standard_bases(
        dna
    )

    total = len(
        standard_bases
    )

    if total == 0:

        return 0.0

    counts = {

        base:
            standard_bases.count(base)

        # Sorted: a set iterates in a different order in every process, which
        # changed the last digit of the entropy sum and, through ties, the
        # reported position of the complexity maximum (fixed 2026-11-19).
        for base in sorted(DNA_STANDARD_BASES)
    }

    entropy = 0.0

    for count in counts.values():

        if count == 0:

            continue

        probability = (
            count / total
        )

        entropy -= (

            probability
            * math.log2(
                probability
            )
        )

    return entropy / 2.0


# ============================================================
# SEQUENCE COMPLEXITY PROFILE
# ============================================================


def calculate_complexity_profile(
    dna,
    window_size=100,
    step=1
):

    profile = []

    if window_size <= 0:

        return profile

    if step <= 0:

        return profile

    if len(dna) < window_size:

        return profile

    for start in range(
        0,
        len(dna) - window_size + 1,
        step
    ):

        end = (
            start + window_size
        )

        window = dna[
            start:end
        ]

        complexity = (
            calculate_sequence_complexity(
                window
            )
        )

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
                complexity
        })

    return profile


# ============================================================
# PROFILE STATISTICS
# ============================================================


def calculate_profile_statistics(
    profile,
    value_key
):

    values = []

    for item in profile:

        value = item.get(
            value_key
        )

        if value is not None:

            values.append(value)

    if not values:

        return {

            "count": 0,

            "minimum": None,

            "maximum": None,

            "average": None,

            "minimum_position": None,

            "maximum_position": None
        }

    minimum = min(values)
    maximum = max(values)

    average = (
        sum(values)
        / len(values)
    )

    minimum_item = min(
        (
            item
            for item in profile
            if item.get(value_key) is not None
        ),
        key=lambda item:
            item[value_key]
    )

    maximum_item = max(
        (
            item
            for item in profile
            if item.get(value_key) is not None
        ),
        key=lambda item:
            item[value_key]
    )

    return {

        "count":
            len(values),

        "minimum":
            minimum,

        "maximum":
            maximum,

        "average":
            average,

        "minimum_position": {

            "start":
                minimum_item[
                    "start_1_based"
                ],

            "end":
                minimum_item[
                    "end_1_based"
                ]
        },

        "maximum_position": {

            "start":
                maximum_item[
                    "start_1_based"
                ],

            "end":
                maximum_item[
                    "end_1_based"
                ]
        }
    }


# ============================================================
# GC PROFILE ANALYSIS
# ============================================================


def analyze_gc_profile(
    dna,
    window_size=100,
    step=1
):

    profile = calculate_gc_content_profile(
        dna,
        window_size,
        step
    )

    statistics = calculate_profile_statistics(
        profile,
        "gc_content"
    )

    return {

        "window_size":
            window_size,

        "step":
            step,

        "profile":
            profile,

        "statistics":
            statistics
    }


# ============================================================
# COMPLEXITY PROFILE ANALYSIS
# ============================================================


def analyze_complexity_profile(
    dna,
    window_size=100,
    step=1
):

    profile = calculate_complexity_profile(
        dna,
        window_size,
        step
    )

    statistics = calculate_profile_statistics(
        profile,
        "complexity"
    )

    return {

        "window_size":
            window_size,

        "step":
            step,

        "profile":
            profile,

        "statistics":
            statistics
    }


# ============================================================
# HOMOPOLYMERS
# ============================================================


def find_homopolymers(
    dna,
    min_length=4
):

    homopolymers = []

    if not dna:

        return homopolymers

    start = 0

    for i in range(
        1,
        len(dna) + 1
    ):

        if (

            i == len(dna)

            or dna[i] != dna[start]
        ):

            length = (
                i - start
            )

            base = dna[start]

            if (

                base in DNA_STANDARD_BASES

                and length >= min_length
            ):

                homopolymers.append({

                    "base":
                        base,

                    "start":
                        start,

                    "end":
                        i - 1,

                    "start_1_based":
                        start + 1,

                    "end_1_based":
                        i,

                    "length":
                        length
                })

            start = i

    return homopolymers


# ============================================================
# CpG ANALYSIS
# ============================================================


def analyze_cpg(dna):

    evaluated_length = len(
        _get_valid_iupac_bases(dna)
    )

    if evaluated_length == 0:

        return {

            "count": 0,

            "density": 0.0,

            "oe_ratio": 0.0,

            "ambiguous_contexts": 0
        }

    definite_cpg_count = 0
    ambiguous_contexts = 0

    for i in range(
        len(dna) - 1
    ):

        first = dna[i]
        second = dna[i + 1]

        if (

            first not in DNA_IUPAC_BASES

            or second not in DNA_IUPAC_BASES
        ):

            continue

        if (

            first == "C"

            and second == "G"
        ):

            definite_cpg_count += 1

            continue

        first_possibilities = (
            DNA_IUPAC_BASES[first]
        )

        second_possibilities = (
            DNA_IUPAC_BASES[second]
        )

        if (

            "C" in first_possibilities

            and

            "G" in second_possibilities
        ):

            ambiguous_contexts += 1

    c_count = dna.count("C")
    g_count = dna.count("G")

    density = (

        definite_cpg_count
        / evaluated_length

    ) * 100

    if (

        c_count == 0

        or g_count == 0
    ):

        oe_ratio = 0.0

    else:

        expected_cpg = (

            c_count
            * g_count

        ) / evaluated_length

        if expected_cpg == 0:

            oe_ratio = 0.0

        else:

            oe_ratio = (

                definite_cpg_count
                / expected_cpg
            )

    return {

        "count":
            definite_cpg_count,

        "density":
            density,

        "oe_ratio":
            oe_ratio,

        "ambiguous_contexts":
            ambiguous_contexts
    }


# ============================================================
# CpG ISLANDS
# ============================================================


def find_cpg_islands(
    dna,
    window_size=200,
    min_gc=50.0,
    min_oe=0.60
):

    qualifying_windows = []

    if window_size <= 0:
        return []

    if len(dna) < window_size:
        return []

    for start in range(
        len(dna) - window_size + 1
    ):

        end = (
            start + window_size
        )

        window = dna[
            start:end
        ]

        if any(
            base not in DNA_STANDARD_BASES
            for base in window
        ):
            continue

        c_count = window.count("C")
        g_count = window.count("G")
        cpg_count = window.count("CG")

        gc_content = (
            (c_count + g_count)
            / window_size
        ) * 100

        if c_count == 0 or g_count == 0:
            oe_ratio = 0.0
        else:
            expected_cpg = (
                c_count * g_count
            ) / window_size

            oe_ratio = (
                cpg_count
                / expected_cpg
            )

        if (
            gc_content >= min_gc
            and oe_ratio >= min_oe
        ):
            qualifying_windows.append(
                (start, end - 1)
            )

    if not qualifying_windows:
        return []

    # Merge overlapping qualifying windows.
    merged_regions = []

    current_start, current_end = (
        qualifying_windows[0]
    )

    for start, end in qualifying_windows[1:]:

        if start <= current_end + 1:
            current_end = max(
                current_end,
                end
            )
        else:
            merged_regions.append(
                (current_start, current_end)
            )

            current_start = start
            current_end = end

    merged_regions.append(
        (current_start, current_end)
    )

    islands = []

    for start, end in merged_regions:

        region = dna[
            start:end + 1
        ]

        length = len(region)

        if length < window_size:
            continue

        c_count = region.count("C")
        g_count = region.count("G")
        cpg_count = region.count("CG")

        gc_content = (
            (c_count + g_count)
            / length
        ) * 100

        if c_count == 0 or g_count == 0:
            oe_ratio = 0.0
        else:
            expected_cpg = (
                c_count * g_count
            ) / length

            oe_ratio = (
                cpg_count
                / expected_cpg
            )

        # Do not report a merged region unless the FINAL
        # reported region itself still satisfies the criteria.
        if (
            gc_content < min_gc
            or oe_ratio < min_oe
        ):
            continue

        islands.append({

            "start":
                start,

            "end":
                end,

            "start_1_based":
                start + 1,

            "end_1_based":
                end + 1,

            "length":
                length,

            "gc_content":
                gc_content,

            "cpg_count":
                cpg_count,

            "oe_ratio":
                oe_ratio,

            "criteria":
                {
                    "min_length":
                        window_size,

                    "min_gc_content":
                        min_gc,

                    "min_oe_ratio":
                        min_oe
                }
        })

    return islands



# ============================================================
# DNA MOLECULAR WEIGHT
# ============================================================


def calculate_dna_molecular_weight(dna):

    """
    Composition-based residue estimate.

    The exact dsDNA estimate is obtained by pairing each possible
    base with its Watson-Crick complement. For ambiguous IUPAC
    symbols, a minimum and maximum are reported.

    This is deliberately labeled an estimate: nucleotide residue
    mass conventions and terminal-group assumptions differ between
    molecular-weight calculators.
    """

    ss_minimum_weight = 0.0
    ss_maximum_weight = 0.0

    ds_minimum_weight = 0.0
    ds_maximum_weight = 0.0

    ambiguous = False
    evaluated_positions = 0

    for symbol in dna:

        possible_bases = DNA_IUPAC_BASES.get(
            symbol
        )

        if possible_bases is None:
            continue

        evaluated_positions += 1

        if len(possible_bases) > 1:
            ambiguous = True

        possible_ss_weights = []
        possible_ds_weights = []

        for base in possible_bases:

            complement = DNA_COMPLEMENT[
                base
            ]

            ss_weight = (
                DNA_NUCLEOTIDE_RESIDUE_MASSES[
                    base
                ]
            )

            ds_weight = (
                DNA_NUCLEOTIDE_RESIDUE_MASSES[
                    base
                ]
                +
                DNA_NUCLEOTIDE_RESIDUE_MASSES[
                    complement
                ]
            )

            possible_ss_weights.append(
                ss_weight
            )

            possible_ds_weights.append(
                ds_weight
            )

        ss_minimum_weight += min(
            possible_ss_weights
        )

        ss_maximum_weight += max(
            possible_ss_weights
        )

        ds_minimum_weight += min(
            possible_ds_weights
        )

        ds_maximum_weight += max(
            possible_ds_weights
        )

    return {

        "ss_min":
            ss_minimum_weight,

        "ss_max":
            ss_maximum_weight,

        "ds_min":
            ds_minimum_weight,

        "ds_max":
            ds_maximum_weight,

        # Backward-compatible keys.
        "min":
            ds_minimum_weight,

        "max":
            ds_maximum_weight,

        "ambiguous":
            ambiguous,

        "evaluated_positions":
            evaluated_positions,

        "method":
            "Composition-based nucleotide residue estimate"
    }



# ============================================================
# AVERAGE DNA MOLECULAR WEIGHT
# ============================================================


def calculate_average_dna_molecular_weight(
    dna
):

    """
    Conventional approximate dsDNA mass using 660 Da per bp.

    This is intentionally separate from the composition-based
    residue estimate above.
    """

    valid_length = len(
        _get_valid_iupac_bases(dna)
    )

    return (
        valid_length
        * AVERAGE_DSDNA_BP_MASS
    )



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



# ============================================================
# SEQUENCE ENGINE SELF-TEST
# ============================================================

def run_sequence_self_test():

    results = {}

    test = "ATGC"

    # Composition.
    counts = count_nucleotides(test)
    results["nucleotide_counts"] = (
        counts == {
            "A": 1,
            "T": 1,
            "C": 1,
            "G": 1
        }
    )

    # GC / AT.
    results["gc_content"] = (
        abs(
            calculate_gc_content(test)
            - 50.0
        )
        < 1e-12
    )

    results["at_content"] = (
        abs(
            calculate_at_content(test)
            - 50.0
        )
        < 1e-12
    )

    # Skew.
    results["gc_skew"] = (
        abs(
            calculate_gc_skew(test)
        )
        < 1e-12
    )

    results["at_skew"] = (
        abs(
            calculate_at_skew(test)
        )
        < 1e-12
    )

    # Complexity: all four bases equally represented => 1.0.
    results["complexity"] = (
        abs(
            calculate_sequence_complexity(test)
            - 1.0
        )
        < 1e-12
    )

    # Complement.
    results["complement"] = (
        get_complement(test)
        == "TACG"
    )

    results["reverse_complement"] = (
        get_reverse_complement(test)
        == "GCAT"
    )

    # CpG.
    cpg_test = "CG"
    cpg = analyze_cpg(cpg_test)

    results["cpg_count"] = (
        cpg["count"] == 1
    )

    # Homopolymer.
    homopolymer = find_homopolymers(
        "AAAA"
    )

    results["homopolymer"] = (
        len(homopolymer) == 1
        and homopolymer[0]["length"] == 4
        and homopolymer[0]["start_1_based"] == 1
        and homopolymer[0]["end_1_based"] == 4
    )

    # IUPAC range.
    ranges = calculate_dna_content_ranges(
        "AN"
    )

    results["iupac_range"] = (
        ranges["gc_min"] == 0.0
        and ranges["gc_max"] == 50.0
        and ranges["at_min"] == 50.0
        and ranges["at_max"] == 100.0
    )

    # Profile availability.
    profile_info = describe_profile_availability(
        "ATGC",
        100,
        "GC content profile"
    )

    results["short_profile"] = (
        profile_info["available"] is False
        and profile_info["reason"]
        == "SEQUENCE_SHORTER_THAN_WINDOW"
    )

    # Conventional dsDNA estimate.
    results["average_ds_weight"] = (
        calculate_average_dna_molecular_weight(
            "ATGC"
        )
        == 2640.0
    )

    # Full pipeline.
    analysis = analyze_sequence(
        "ATGCGCGT"
    )

    results["full_pipeline"] = (
        analysis["sequence_length"] == 8
        and "interpretation" in analysis
        and "sequence_statistics" in analysis
        and "sequence_forms" in analysis
        and "cpg_island_criteria" in analysis
    )

    return results


if __name__ == "__main__":
    from ui import style  # only when run as a script

    self_test = run_sequence_self_test()

    style.header('Sequence engine self-test')

    for name, passed in self_test.items():

        print(
            f"{name:<30}"
            f"{'PASS' if passed else 'FAIL'}"
        )

    print(minor_rule())

    if all(
        self_test.values()
    ):

        print(
            "RESULT: ALL TESTS PASSED"
        )

    else:

        print(
            "RESULT: ONE OR MORE TESTS FAILED"
        )