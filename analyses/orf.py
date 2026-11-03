from itertools import product
from collections import Counter

from analyses.sequence import get_reverse_complement
from core.seqtools import codon_dict
from ui.style import RULE_WIDTH, major_rule, minor_rule
from ui import style


# ============================================================
# ORF ANALYSIS
# ============================================================


# ============================================================
# CONSTANTS
# ============================================================

STANDARD_START_CODONS = {"ATG"}

ALTERNATIVE_START_CODONS = {
    "GTG",
    "TTG"
}

ALL_START_CODONS = (
    STANDARD_START_CODONS
    | ALTERNATIVE_START_CODONS
)

STANDARD_STOP_CODONS = {
    "TAA",
    "TAG",
    "TGA"
}

# ------------------------------------------------------------
# Default search policy (audit P0-9, 2026-10-29)
# ------------------------------------------------------------
# Defaults follow NCBI ORFfinder conventions: minimal ORF length 75 nt
# (stop codon included), ATG-only start codons, and one ORF per stop codon
# (the longest, i.e. outermost in-frame start). These are screening
# conventions, not biological thresholds: short ORFs, non-ATG starts and
# nested ORFs exist and can be enabled explicitly.
#
# Previous defaults (still selectable): min length 10 nt, starts
# ATG/GTG/TTG, every in-frame start reported, heuristic score ranking.
DEFAULT_MIN_ORF_LENGTH = 75

START_CODON_POLICIES = {
    "atg_only": frozenset(STANDARD_START_CODONS),
    "atg_and_alternative": frozenset(ALL_START_CODONS),
}

DEFAULT_START_CODON_POLICY = "atg_only"
DEFAULT_REPORT_NESTED = False
DEFAULT_RANKING = "length"
ORF_RANKING_MODES = ("length", "heuristic")

LEGACY_ORF_SETTINGS = {
    "min_length": 10,
    "start_codon_policy": "atg_and_alternative",
    "nested": True,
    "ranking": "heuristic",
}


DNA_IUPAC = {
    "A": ("A",),
    "T": ("T",),
    "C": ("C",),
    "G": ("G",),

    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),

    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "T", "C", "G")
}


# ============================================================
# GENETIC CODE
# ============================================================

# Standard code (NCBI table 1) from core.seqtools; stops are "STOP" here.
GENETIC_CODE = codon_dict(1, stop_symbol="STOP")


# ============================================================
# IUPAC CODON HANDLING
# ============================================================

def expand_iupac_codon(codon):

    codon = codon.upper()

    if len(codon) != 3:
        return []

    possibilities = []

    for base in codon:

        if base not in DNA_IUPAC:
            return []

        possibilities.append(
            DNA_IUPAC[base]
        )

    return [
        "".join(combination)
        for combination in product(
            *possibilities
        )
    ]


def _stop_codons_of(codon_table):
    if codon_table is None:
        return STANDARD_STOP_CODONS
    return {codon for codon, aa in codon_table.items() if aa == "STOP"}


def codon_table_for(genetic_code=1):
    """Return the local codon->amino-acid dict ('STOP' for stops) for an NCBI table.

    Table 1 returns the module's own GENETIC_CODE (cross-checked against
    core.seqtools in the test suite); other tables come from core.seqtools.
    """
    if genetic_code in (None, 1):
        return GENETIC_CODE
    from core.seqtools import get_genetic_code
    table = get_genetic_code(int(genetic_code))
    return {codon: ("STOP" if aa == "*" else aa) for codon, aa in table.codons.items()}


def start_codons_for(genetic_code=1, start_codon_policy=None):
    """Start codons for a genetic code and a start-codon policy.

    Table 1 keeps the historical Synapse policies (ATG, or ATG+GTG+TTG).
    For any table, "atg_only" means ATG. For tables other than 1,
    "atg_and_alternative" uses the start codons NCBI lists for that table;
    "table_starts" does so for every table (table 1: TTG, CTG, ATG).
    """
    policy = start_codon_policy or DEFAULT_START_CODON_POLICY
    if policy not in START_CODON_POLICIES and policy != "table_starts":
        raise ValueError(
            f"Unknown start codon policy '{policy}'. "
            f"Use one of: {', '.join(list(START_CODON_POLICIES) + ['table_starts'])}."
        )
    if policy == "atg_only":
        return set(STANDARD_START_CODONS)
    if policy == "table_starts" or genetic_code not in (None, 1):
        from core.seqtools import get_genetic_code
        return set(get_genetic_code(int(genetic_code or 1)).start_codons)
    return set(START_CODON_POLICIES[policy])


def initiation_codons(genetic_code=1):
    """Codons that can initiate an ORF in a genetic code (audit D1).

    The union of the start codons of every start-codon policy for the
    table: the NCBI start codons of that table (core.seqtools, checked
    against Biopython) plus, for table 1, Synapse's documented
    "atg_and_alternative" policy (ATG, GTG, TTG). As the first codon of an
    ORF such a codon is the initiation codon and is decoded as methionine;
    anywhere else it is read with the elongation code.
    """
    policies = list(START_CODON_POLICIES) + ["table_starts"]
    return frozenset().union(
        *(start_codons_for(genetic_code, policy) for policy in policies)
    )


def is_initiation_codon(codon, genetic_code=1):
    """True when every IUPAC reading of the codon is an initiation codon.

    This is the rule the ORF finder applies to a "definite" start codon,
    so the codon that starts a Synapse ORF is always decoded as initiator.
    Accepts DNA or RNA letters.
    """
    codon = str(codon).upper().replace("U", "T")
    return analyze_codon_identity(codon, initiation_codons(genetic_code)) == "definite"


def _genetic_code_name(genetic_code=1):
    if genetic_code in (None, 1):
        return "NCBI translation table 1 (standard)"
    from core.seqtools import get_genetic_code
    table = get_genetic_code(int(genetic_code))
    return f"NCBI translation table {table.table_id} ({table.name})"


def analyze_codon_identity(
    codon,
    target_codons
):

    expanded = expand_iupac_codon(
        codon
    )

    if not expanded:
        return "invalid"

    matches = [
        expanded_codon in target_codons
        for expanded_codon in expanded
    ]

    if all(matches):
        return "definite"

    if any(matches):
        return "possible"

    return "no"


# ============================================================
# START CODON CLASSIFICATION
# ============================================================

def _is_alternative_start(orf):
    """A non-ATG start that initiates in the ORF's genetic code (audit S19)."""
    codon = orf["start_codon"]
    return codon != "ATG" and is_initiation_codon(codon, orf.get("genetic_code") or 1)


def classify_start_codon(codon, genetic_code=1):
    """Canonical (ATG), Alternative, Possible or None for a genetic code.

    The start set is ``initiation_codons(genetic_code)`` (audit S19; it was
    the fixed ATG/GTG/TTG set, so CTG in table 1 or ATA/ATT in table 2
    were "None").
    """

    status = analyze_codon_identity(
        codon,
        initiation_codons(genetic_code)
    )

    if status == "definite":

        if codon == "ATG":
            return "Canonical"

        if codon in initiation_codons(genetic_code):
            return "Alternative"

        return "Ambiguous-Definite"

    if status == "possible":
        return "Possible"

    return "None"


# ============================================================
# ORF AMBIGUITY ANALYSIS
# ============================================================

def analyze_orf_ambiguity(orf_sequence, codon_table=None):
    codon_table = GENETIC_CODE if codon_table is None else codon_table

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
    stop_codons = _stop_codons_of(codon_table)

    codons = []

    for position in range(
        0,
        len(coding_sequence) - 2,
        3
    ):

        codon = coding_sequence[
            position:
            position + 3
        ]

        if len(codon) < 3:
            continue

        if codon in stop_codons:
            break

        codons.append(codon)

    counts = Counter(codons)

    standard_codons = [
        codon
        for codon in counts
        if all(
            base in "ATCG"
            for base in codon
        )
    ]

    gc3_values = []

    for codon in codons:

        if all(
            base in "ATCG"
            for base in codon
        ):

            gc3_values.append(
                1
                if codon[2] in "GC"
                else 0
            )

    if gc3_values:

        gc3_content = (
            sum(gc3_values)
            / len(gc3_values)
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
            len(standard_codons),

        "gc3_content":
            gc3_content
    }


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

    if best_motif is None:

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


# ============================================================
# ORF FINDER
# ============================================================

def resolve_start_codons(start_codons=None, start_codon_policy=None):
    """Return the start-codon set from an explicit set or a named policy."""
    if start_codons is not None:
        return set(start_codons)
    policy = start_codon_policy or DEFAULT_START_CODON_POLICY
    if policy not in START_CODON_POLICIES:
        raise ValueError(
            f"Unknown start codon policy '{policy}'. "
            f"Use one of: {', '.join(START_CODON_POLICIES)}."
        )
    return set(START_CODON_POLICIES[policy])


def _orf_sort_key(orf, ranking):
    if ranking == "heuristic":
        return (
            orf["heuristic_score"]["score"],
            orf["cds_quality"]["score"],
            orf["length"],
            orf["status"] == "Complete",
            orf["start_type"] == "Canonical",
        )
    # "length": longest ORF first (ORFfinder-style); ties broken by
    # completeness, canonical start and then the heuristic score.
    return (
        orf["length"],
        orf["status"] == "Complete",
        orf["start_type"] == "Canonical",
        orf["heuristic_score"]["score"],
    )


def _keep_longest_per_stop(orfs):
    """Keep only the outermost ORF for each (strand, frame, 3' end)."""
    best = {}
    for orf in orfs:
        three_prime = orf["stop"] if orf["strand"] == "+" else orf["start"]
        key = (orf["strand"], orf["frame"], three_prime, orf["status"])
        current = best.get(key)
        if current is None or orf["length"] > current["length"]:
            best[key] = orf
    return list(best.values())


def find_orfs(
    dna,
    include_partial=True,
    min_length=DEFAULT_MIN_ORF_LENGTH,
    start_codons=None,
    stop_codons=None,
    start_codon_policy=None,
    nested=DEFAULT_REPORT_NESTED,
    ranking=DEFAULT_RANKING,
    genetic_code=1,
):
    """Find open reading frames on both strands.

    Parameters
    ----------
    min_length:
        Minimum ORF length in nucleotides, stop codon included.
    start_codon_policy:
        "atg_only" (default) or "atg_and_alternative" (ATG/GTG/TTG). Ignored
        when ``start_codons`` is given explicitly.
    nested:
        False (default): report only the longest ORF per stop codon.
        True: report every in-frame start as its own (nested) ORF.
    ranking:
        "length" (default, longest first) or "heuristic" (legacy
        heuristic score).

    Coordinates in the result are 0-based and inclusive on the forward
    strand ("start" <= "stop"), for both strands.

    The sequence is uppercased with whitespace removed (soft-masked input
    is searched like any other); symbols outside the IUPAC DNA alphabet,
    RNA's U included, raise ValueError (audit D14-F8). ``core.api``
    converts RNA explicitly before calling this.
    """

    from core.alphabets import DNA_IUPAC, normalize_sequence, symbols_outside

    dna = normalize_sequence(dna)
    invalid = symbols_outside(dna, DNA_IUPAC)
    if invalid:
        shown = ", ".join(f"{symbol!r} at {position}" for position, symbol in invalid[:5])
        raise ValueError(
            f"ORF search needs a DNA sequence; found {shown}"
            + (" ..." if len(invalid) > 5 else "") + "."
        )

    if ranking not in ORF_RANKING_MODES:
        raise ValueError(
            f"ranking must be one of {ORF_RANKING_MODES}."
        )

    codon_table = codon_table_for(genetic_code)
    if start_codons is None and (genetic_code not in (None, 1) or start_codon_policy == "table_starts"):
        start_codons = start_codons_for(genetic_code, start_codon_policy)

    start_codons = resolve_start_codons(
        start_codons,
        start_codon_policy,
    )

    if stop_codons is None:
        stop_codons = _stop_codons_of(codon_table)

    start_codons = set(
        start_codons
    )

    stop_codons = set(
        stop_codons
    )

    if min_length < 0:

        raise ValueError(
            "min_length must be "
            "greater than or equal to 0."
        )

    reverse_dna = (
        get_reverse_complement(dna)
    )

    strands = [
        ("+", dna),
        ("-", reverse_dna)
    ]

    sequence_length = len(dna)

    orfs = []

    for strand, sequence in strands:

        for frame_offset in range(3):

            frame_number = (
                frame_offset + 1
                if strand == "+"
                else -(frame_offset + 1)
            )

            position = frame_offset

            while (
                position
                <= len(sequence) - 3
            ):

                codon = sequence[
                    position:
                    position + 3
                ]

                start_status = (
                    analyze_codon_identity(
                        codon,
                        start_codons
                    )
                )

                if start_status == "definite":

                    stop_position = (
                        position + 3
                    )

                    stop_codon = None

                    possible_stop_sites = []

                    while (
                        stop_position
                        <= len(sequence) - 3
                    ):

                        current_codon = (
                            sequence[
                                stop_position:
                                stop_position + 3
                            ]
                        )

                        stop_status = (
                            analyze_codon_identity(
                                current_codon,
                                stop_codons
                            )
                        )

                        if (
                            stop_status
                            == "definite"
                        ):

                            stop_codon = (
                                current_codon
                            )

                            break

                        if (
                            stop_status
                            == "possible"
                        ):

                            possible_stop_sites.append({
                                "codon":
                                    current_codon,

                                "position":
                                    stop_position + 1
                            })

                        stop_position += 3

                    # ====================================================
                    # COMPLETE ORF
                    # ====================================================

                    if stop_codon is not None:

                        orf_sequence = sequence[
                            position:
                            stop_position + 3
                        ]

                        if strand == "+":

                            genomic_start = (
                                position
                            )

                            genomic_stop = (
                                stop_position + 2
                            )

                        else:

                            genomic_start = (
                                sequence_length
                                - (
                                    stop_position
                                    + 3
                                )
                            )

                            genomic_stop = (
                                sequence_length
                                - position
                                - 1
                            )

                        if (
                            len(orf_sequence)
                            >= min_length
                        ):

                            ambiguity = (
                                analyze_orf_ambiguity(
                                    orf_sequence, codon_table=codon_table
                                )
                            )

                            codon_count = (
                                len(orf_sequence)
                                // 3
                            )

                            protein_length = (
                                codon_count - 1
                            )

                            start_type = (
                                "Canonical"
                                if codon == "ATG"
                                else "Alternative"
                            )

                            start_context = (
                                analyze_start_context(
                                    sequence,
                                    position
                                )
                            )

                            codon_usage = (
                                analyze_codon_usage(
                                    orf_sequence, codon_table=codon_table
                                )
                            )

                            orf = {
                                "genetic_code":
                                    genetic_code,
                                "strand":
                                    strand,

                                "frame":
                                    frame_number,

                                "start":
                                    genomic_start,

                                "stop":
                                    genomic_stop,

                                "length":
                                    len(orf_sequence),

                                "codon_count":
                                    codon_count,

                                "protein_length":
                                    protein_length,

                                "start_codon":
                                    codon,

                                "start_type":
                                    start_type,

                                "start_status":
                                    start_status,

                                "stop_codon":
                                    stop_codon,

                                "status":
                                    "Complete",

                                "sequence":
                                    orf_sequence,

                                "possible_stop_sites":
                                    possible_stop_sites,

                                "ambiguous_codon_count":
                                    ambiguity[
                                        "ambiguous_codon_count"
                                    ],

                                "ambiguous_codons":
                                    ambiguity[
                                        "ambiguous_codons"
                                    ],

                                "deterministic_ambiguous_count":
                                    ambiguity[
                                        "deterministic_ambiguous_count"
                                    ],

                                "uncertain_codon_count":
                                    ambiguity[
                                        "uncertain_codon_count"
                                    ],

                                "possible_stop_codon_count":
                                    ambiguity[
                                        "possible_stop_codon_count"
                                    ],

                                "translation_safe":
                                    ambiguity[
                                        "translation_safe"
                                    ],

                                "start_context":
                                    start_context,

                                "codon_usage":
                                    codon_usage,

                                "codon_diversity":
                                    calculate_codon_diversity(
                                        orf_sequence, codon_table=codon_table
                                    )
                            }

                            orf["cds_quality"] = (
                                calculate_cds_quality(
                                    orf
                                )
                            )

                            orf["heuristic_score"] = (
                                calculate_orf_heuristic_score(
                                    orf
                                )
                            )

                            orfs.append(orf)

                    # ====================================================
                    # PARTIAL ORF
                    # ====================================================

                    elif include_partial:

                        remaining_length = (
                            len(sequence)
                            - position
                        )

                        usable_length = (
                            remaining_length
                            - (
                                remaining_length
                                % 3
                            )
                        )

                        if usable_length >= 3:

                            orf_sequence = sequence[
                                position:
                                position + usable_length
                            ]

                            if (
                                len(orf_sequence)
                                >= min_length
                            ):

                                last_position = (
                                    position
                                    + usable_length
                                    - 1
                                )

                                if strand == "+":

                                    genomic_start = (
                                        position
                                    )

                                    genomic_stop = (
                                        last_position
                                    )

                                else:

                                    genomic_start = (
                                        sequence_length
                                        - (
                                            position
                                            + usable_length
                                        )
                                    )

                                    genomic_stop = (
                                        sequence_length
                                        - position
                                        - 1
                                    )

                                ambiguity = (
                                    analyze_orf_ambiguity(
                                        orf_sequence, codon_table=codon_table
                                    )
                                )

                                codon_count = (
                                    len(orf_sequence)
                                    // 3
                                )

                                protein_length = (
                                    codon_count
                                )

                                start_type = (
                                    "Canonical"
                                    if codon == "ATG"
                                    else "Alternative"
                                )

                                start_context = (
                                    analyze_start_context(
                                        sequence,
                                        position
                                    )
                                )

                                codon_usage = (
                                    analyze_codon_usage(
                                        orf_sequence, codon_table=codon_table
                                    )
                                )

                                orf = {
                                    "genetic_code":
                                        genetic_code,
                                    "strand":
                                        strand,

                                    "frame":
                                        frame_number,

                                    "start":
                                        genomic_start,

                                    "stop":
                                        genomic_stop,

                                    "length":
                                        len(orf_sequence),

                                    "codon_count":
                                        codon_count,

                                    "protein_length":
                                        protein_length,

                                    "start_codon":
                                        codon,

                                    "start_type":
                                        start_type,

                                    "start_status":
                                        start_status,

                                    "stop_codon":
                                        None,

                                    "status":
                                        "Partial - Missing Stop",

                                    "sequence":
                                        orf_sequence,

                                    "possible_stop_sites":
                                        possible_stop_sites,

                                    "ambiguous_codon_count":
                                        ambiguity[
                                            "ambiguous_codon_count"
                                        ],

                                    "ambiguous_codons":
                                        ambiguity[
                                            "ambiguous_codons"
                                        ],

                                    "deterministic_ambiguous_count":
                                        ambiguity[
                                            "deterministic_ambiguous_count"
                                        ],

                                    "uncertain_codon_count":
                                        ambiguity[
                                            "uncertain_codon_count"
                                        ],

                                    "possible_stop_codon_count":
                                        ambiguity[
                                            "possible_stop_codon_count"
                                        ],

                                    "translation_safe":
                                        ambiguity[
                                            "translation_safe"
                                        ],

                                    "start_context":
                                        start_context,

                                    "codon_usage":
                                        codon_usage,

                                    "codon_diversity":
                                        calculate_codon_diversity(
                                            orf_sequence, codon_table=codon_table
                                        )
                                }

                                orf["cds_quality"] = (
                                    calculate_cds_quality(
                                        orf
                                    )
                                )

                                orf["heuristic_score"] = (
                                    calculate_orf_heuristic_score(
                                        orf
                                    )
                                )

                                orfs.append(orf)

                position += 3

    if not nested:
        orfs = _keep_longest_per_stop(orfs)

    # ========================================================
    # RANK ORFS
    # ========================================================

    orfs.sort(
        key=lambda orf: _orf_sort_key(orf, ranking),
        reverse=True
    )

    for rank, orf in enumerate(
        orfs,
        start=1
    ):

        orf["rank"] = rank

    return orfs


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

def print_orf_report(
    orfs,
    min_length=DEFAULT_MIN_ORF_LENGTH
):

    style.header('ORF ANALYSIS')

    print()
    style.field('Minimum ORF Length', f"{min_length} bp")

    if not orfs:

        print("\nNo ORFs found.")

        style.section('ORF Definition')

        style.field('Search Space', "Both strands / six reading frames")

        style.field('Start Codons', "ATG, GTG, TTG")

        style.field('Stop Codons', "TAA, TAG, TGA")

        style.field('Minimum Length', f"{min_length} bp")

        print()
        print(major_rule())

        return

    # Infer sequence length from genomic coordinates.
    sequence_length = max(
        orf["stop"]
        for orf in orfs
    ) + 1

    summary = summarize_orfs(
        orfs,
        sequence_length
    )

    # ========================================================
    # SUMMARY
    # ========================================================

    style.section('ORF Summary')

    style.field('ORFs Found', f"{summary['total']}")

    style.field('Forward Strand', f"{summary['forward']}")

    style.field('Reverse Strand', f"{summary['reverse']}")

    style.field('Complete ORFs', f"{summary['complete']}")

    style.field('Partial ORFs', f"{summary['partial']}")

    style.field('Canonical Start (ATG)', f"{summary['canonical_start']}")

    style.field('Alternative Start', f"{summary['alternative_start']}")

    style.field('ORFs with Ambiguity', f"{summary['ambiguous']}")

    style.field('Translation-Safe', f"{summary['translation_safe']}")

    style.field('Kozak Context Available', f"{summary['kozak_available']}")

    style.field('Shine-Dalgarno Detected', f"{summary['shine_dalgarno_detected']}")

    style.field('Average ORF Length', f"{summary['average_length']:.1f} bp")

    style.field('Longest ORF', f"{summary['longest']['length']} bp")

    style.field('ORF Density', f"{summary['orf_density']:.2f}%")

    # ========================================================
    # FRAME STATISTICS
    # ========================================================

    style.section('Reading Frame Statistics')

    for frame in (
        "+1",
        "+2",
        "+3",
        "-1",
        "-2",
        "-3"
    ):

        statistics = (
            summary["frame_statistics"]
            [frame]
        )

        print(
            f"Frame {frame}: "
            f"{statistics['count']} ORFs | "
            f"{statistics['complete']} complete | "
            f"{statistics['partial']} partial"
        )

        if statistics["longest"]:

            print(
                f"           Longest: "
                f"{statistics['longest']['length']} bp"
            )

    # ========================================================
    # SEARCH PARAMETERS
    # ========================================================

    style.section('Search Parameters')

    style.field('Search Space', "Both strands / six reading frames")

    style.field('Start Codons', "ATG, GTG, TTG")

    style.field('Canonical Start', "ATG")

    style.field('Alternative Starts', "GTG, TTG")

    style.field('Stop Codons', "TAA, TAG, TGA")

    style.field('Minimum Length', f"{min_length} bp")

    style.field('Stop Rule', "First definite downstream "
        "in-frame stop")

    style.field('Ambiguous Starts', "Only definite starts initiate ORFs")

    style.field('Ambiguous Stops', "Possible-only stops do not terminate ORFs")

    # ========================================================
    # RANKED ORFS
    # ========================================================

    style.section('Ranked ORF Candidates')

    for orf in orfs:

        print()

        print(
            f"ORF #{orf['rank']}"
        )

        print(minor_rule())

        style.field('Heuristic score', f"{orf['heuristic_score']['score']}/100 "
            f"({orf['heuristic_score']['classification']})")

        style.field('CDS Quality', f"{orf['cds_quality']['score']}/100 "
            f"({orf['cds_quality']['classification']})")

        style.field('Strand', f"{orf['strand']}")

        style.field('Reading Frame', f"{orf['frame']:+d}")

        style.field('Start Position', f"{orf['start'] + 1}")

        style.field('Stop Position', f"{orf['stop'] + 1}")

        style.field('Length', f"{orf['length']} bp")

        style.field('Codons', f"{orf['codon_count']}")

        style.field('Potential Protein', f"{orf['protein_length']} aa")

        style.field('Start Codon', f"{orf['start_codon']}")

        style.field('Start Type', f"{orf['start_type']}")

        style.field('Status', f"{orf['status']}")

        if orf["stop_codon"] is None:

            style.field('Stop Codon', "Not detected")

        else:

            style.field('Stop Codon', f"{orf['stop_codon']}")

        # ====================================================
        # START CONTEXT
        # ====================================================

        style.section('Start Context')

        kozak = (
            orf["start_context"]
            ["kozak"]
        )

        if kozak["available"]:

            style.field('Kozak Context', f"{kozak['context']}")

            style.field('Kozak Score', f"{kozak['score']}/3")

            print(
                f"Kozak Classification:"
                f" {kozak['classification']}"
            )

        else:

            style.field('Kozak Context', "Unavailable")

        shine = (
            orf["start_context"]
            ["shine_dalgarno"]
        )

        if shine["detected"]:

            style.field('Shine-Dalgarno', f"{shine['motif']}")

            style.field('SD Distance', f"{shine['distance_to_start']} nt")

        else:

            style.field('Shine-Dalgarno', "Not detected")

        # ====================================================
        # CODON STATISTICS
        # ====================================================

        codon_usage = (
            orf["codon_usage"]
        )

        style.section('Codon Statistics')

        style.field('Total Codons', f"{codon_usage['total_codons']}")

        style.field('Unique Codons', f"{codon_usage['unique_codons']}")

        style.field('Codon Diversity', f"{orf['codon_diversity']:.3f}")

        if (
            codon_usage["gc3_content"]
            is not None
        ):

            style.field('GC3 Content', f"{codon_usage['gc3_content']:.2f}%")

        else:

            style.field('GC3 Content', "Not available")

        # ====================================================
        # AMBIGUITY
        # ====================================================

        style.section('Sequence Ambiguity')

        style.field('Ambiguous Codons', f"{orf['ambiguous_codon_count']}")

        style.field('Deterministic', f"{orf['deterministic_ambiguous_count']}")

        print(
            f"Translation-Uncertain:"
            f" {orf['uncertain_codon_count']}"
        )

        print(
            f"Possible Stop Codons:"
            f" {orf['possible_stop_codon_count']}"
        )

        if orf["translation_safe"]:

            if (
                orf["ambiguous_codon_count"]
                == 0
            ):

                style.field('Translation Status', " Unambiguous")

            else:

                style.field('Translation Status', " Deterministically resolvable")

        else:

            style.field('Translation Status', " Contains uncertainty")

        # ====================================================
        # AMBIGUOUS CODON DETAILS
        # ====================================================

        if orf["ambiguous_codons"]:

            style.section('Ambiguous Codon Details')

            for item in (
                orf["ambiguous_codons"]
            ):

                codon_number = (
                    (
                        item["position"] - 1
                    )
                    // 3
                ) + 1

                if (
                    item[
                        "translation_status"
                    ]
                    == "Deterministic"
                ):

                    print(
                        f"Codon {codon_number}: "
                        f"{item['codon']} -> "
                        f"{item['translation']}"
                    )

                else:

                    possibilities = (
                        "/".join(
                            item[
                                "possibilities"
                            ]
                        )
                    )

                    if item[
                        "contains_stop"
                    ]:

                        possibilities += (
                            "/STOP"
                        )

                    print(
                        f"Codon {codon_number}: "
                        f"{item['codon']} -> "
                        f"{possibilities}"
                    )

        # ====================================================
        # SEQUENCE
        # ====================================================

        style.section('Sequence')

        print(
            create_orf_preview(
                orf["sequence"]
            )
        )

    # ========================================================
    # INTERPRETATION
    # ========================================================

    style.section('Interpretation')

    print(
        "ORFs are sequence-defined candidate "
        "coding regions."
    )

    print(
        "ORF detection alone does not prove "
        "biological expression or protein function."
    )

    print(
        "ATG is treated as the canonical start codon, "
        "while GTG and TTG are treated as alternative "
        "start codons."
    )

    print(
        "The analysis evaluates all six reading frames "
        "on both DNA strands."
    )

    print(
        "IUPAC ambiguity is explicitly evaluated rather "
        "than silently ignored."
    )

    print(
        "Kozak and Shine-Dalgarno results are sequence-based "
        "context indicators and do not prove translation."
    )

    print(
        "Heuristic and CDS quality scores are "
        "ranking systems, not experimental evidence."
    )

    print()

    print(major_rule())



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

            if a_contains_b and not b_contains_a:
                relation = "Nested"
                outer = a
                inner = b
                nested_pairs.append((a, b))
            elif b_contains_a and not a_contains_b:
                relation = "Nested"
                outer = b
                inner = a
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

            if relation == "Nested":
                record["outer_orf_rank"] = outer["rank"]
                record["inner_orf_rank"] = inner["rank"]

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


# ============================================================
# ORF COMPARISON REPORT
# ============================================================

def print_orf_candidate_comparison(orfs, top_n=5):
    """Print a compact comparison of the highest-ranked candidates."""
    candidates = compare_orf_candidates(orfs, top_n=top_n)

    style.section('Top ORF Candidate Comparison')

    if not candidates:
        style.note("No ORF candidates available.")
        return

    for candidate in candidates:
        print(
            f"#{candidate['rank']} | "
            f"Frame {candidate['frame']:+d} | "
            f"{candidate['length']} bp | "
            f"{candidate['protein_length']} aa | "
            f"{candidate['start_codon']} | "
            f"{candidate['status']} | "
            f"Heuristic {candidate['heuristic_score']}/100 | "
            f"CDS {candidate['cds_quality_score']}/100"
        )


def print_orf_relationship_report(orfs):
    """Print detected overlap/nesting relationships."""
    relationships = analyze_orf_relationships(orfs)

    style.section('ORF Relationships')

    style.field('Overlapping pairs', f"{relationships['overlap_count']}")
    style.field('Nested pairs', f"{relationships['nested_count']}")

    if not relationships["relationships"]:
        print("No overlapping or nested ORF relationships detected.")
        return

    for item in relationships["relationships"]:
        print(
            f"ORF #{item['orf_a_rank']} ↔ ORF #{item['orf_b_rank']} | "
            f"{item['relation']} | "
            f"overlap {item['overlap_length']} bp | "
            f"same strand: {'YES' if item['same_strand'] else 'NO'} | "
            f"same frame: {'YES' if item['same_frame'] else 'NO'}"
        )

# ============================================================
# ORF UI / SELECTION
# ============================================================

def _orf_status_label(orf):
    return "COMPLETE" if orf["status"] == "Complete" else "PARTIAL"


def print_orf_summary_card(analysis):
    """Compact executive summary for the ORF UI."""
    style.header("ORF ANALYSIS SUMMARY")

    summary = analysis.get("summary", {})
    top = analysis.get("top_candidate")

    style.field('Sequence length', f"{analysis.get('sequence_length', 0)} bp")
    style.field('ORF candidates', f"{summary.get('total', 0)}")

    if top:
        style.field('Best candidate', f"ORF #{top['rank']}")
        style.field('Heuristic score', f"{top['heuristic_score']['score']}/100 "
            f"({top['heuristic_score']['classification']})")
        style.field('CDS quality', f"{top['cds_quality']['score']}/100 "
            f"({top['cds_quality']['classification']})")
        style.field('Status', f"{top['status']}")
    else:
        style.field('Best candidate', "None")

    print()
    style.field('Candidates', f"{summary.get('complete', 0)} complete, "
        f"{summary.get('partial', 0)} partial; "
        f"{summary.get('forward', 0)} forward, "
        f"{summary.get('reverse', 0)} reverse")


def print_orf_candidates(orfs, top_n=None):
    """Display a clean candidate list without asking for input."""
    style.header("ORF CANDIDATES")

    if not orfs:
        style.note("No ORF candidates available.")
        return

    candidates = orfs if top_n is None else orfs[:max(0, top_n)]

    widths = (8, 6, 9, 8, 6, 10)
    print(style.columns(
        ("ID", "Frame", "Length", "Protein", "Start", "Status", "Heuristic"), widths))
    print(minor_rule())

    for orf in candidates:
        print(style.columns((
            f"ORF-{orf['rank']:03d}",
            f"{orf['frame']:+d}",
            f"{orf['length']} bp",
            f"{orf['protein_length']} aa",
            orf["start_codon"],
            _orf_status_label(orf),
            f"{orf['heuristic_score']['score']}/100",
        ), widths))


def print_orf_candidate_details(orf):
    """Detailed view for one candidate."""
    style.header(f"ORF-{orf['rank']:03d} — CANDIDATE DETAILS")

    style.field('Identifier', f"ORF-{orf['rank']:03d}")
    style.field('Rank', f"#{orf['rank']}")
    style.field('Strand', f"{orf['strand']}")
    style.field('Reading frame', f"{orf['frame']:+d}")
    style.field('Genomic start', f"{orf['start'] + 1}")
    style.field('Genomic stop', f"{orf['stop'] + 1}")
    style.field('Length', f"{orf['length']} bp")
    style.field('Potential protein', f"{orf['protein_length']} aa")
    style.field('Start codon', f"{orf['start_codon']}")
    style.field('Start type', f"{orf['start_type']}")
    style.field('Status', f"{orf['status']}")

    style.section('Quality')
    style.field('Heuristic score', f"{orf['heuristic_score']['score']}/100 "
        f"({orf['heuristic_score']['classification']})")
    style.field('CDS quality', f"{orf['cds_quality']['score']}/100 "
        f"({orf['cds_quality']['classification']})")
    style.field('Translation safe', f"{'YES' if orf['translation_safe'] else 'NO'}")

    style.section('Start Context')

    kozak = orf["start_context"]["kozak"]
    if kozak["available"]:
        style.field('Kozak context', f"{kozak['context']}")
        style.field('Kozak score', f"{kozak['score']}/3")
        style.field('Kozak classification', f"{kozak['classification']}")
    else:
        style.field('Kozak context', "Unavailable")

    shine = orf["start_context"]["shine_dalgarno"]
    if shine["detected"]:
        style.field('Shine-Dalgarno', f"{shine['motif']}")
        style.field('SD distance', f"{shine['distance_to_start']} nt")
    else:
        style.field('Shine-Dalgarno', "Not detected")

    style.section('Sequence')
    print(create_orf_preview(orf["sequence"]))


def print_orf_comparison_table(orfs, top_n=5):
    """Display the highest-ranked candidates side by side."""
    candidates = compare_orf_candidates(orfs, top_n=top_n)

    style.header("COMPARE ORFs")

    if not candidates:
        style.note("No ORF candidates available.")
        return

    print(
        f"{'Feature':<24}"
        + "".join(f"ORF #{c['rank']:<14}" for c in candidates)
    )
    # A table rule spans its own table; the standard width is the floor.
    print(minor_rule() if 24 + 14 * len(candidates) <= RULE_WIDTH
          else "-" * (24 + 14 * len(candidates)))

    rows = [
        ("Strand", lambda c: c["strand"]),
        ("Reading frame", lambda c: f"{c['frame']:+d}"),
        ("Length", lambda c: f"{c['length']} bp"),
        ("Protein", lambda c: f"{c['protein_length']} aa"),
        ("Start codon", lambda c: c["start_codon"]),
        ("Status", lambda c: c["status"]),
        ("Heuristic score", lambda c: f"{c['heuristic_score']}/100"),
        ("CDS quality", lambda c: f"{c['cds_quality_score']}/100"),
        (
            "Translation safe",
            lambda c: "YES" if c["translation_safe"] else "NO",
        ),
        ("Ambiguous codons", lambda c: str(c["ambiguous_codons"])),
    ]

    for label, getter in rows:
        values = "".join(f"{getter(c):<14}" for c in candidates)
        print(f"{label:<24}{values}")


def print_orf_relationships(orfs):
    """Display sequence-coordinate relationships between ORFs."""
    data = analyze_orf_relationships(orfs)

    style.header("ORF RELATIONSHIPS")

    style.field('Relationships found', f"{data['relationship_count']}")
    style.field('Overlapping pairs', f"{data['overlap_count']}")
    style.field('Nested pairs', f"{data['nested_count']}")

    if not data["relationships"]:
        print("\nNo overlapping or nested ORF relationships detected.")
        return

    print()

    for item in data["relationships"]:
        print(
            f"ORF #{item['orf_a_rank']} <-> ORF #{item['orf_b_rank']}"
        )
        style.field('Relationship', f"{item['relation']}")
        style.field('Overlap', f"{item['overlap_length']} bp")
        style.field('ORF A coverage', f"{item['overlap_fraction_a'] * 100:.1f}%")
        style.field('ORF B coverage', f"{item['overlap_fraction_b'] * 100:.1f}%")
        style.field('Same strand', f"{'YES' if item['same_strand'] else 'NO'}")
        style.field('Same frame', f"{'YES' if item['same_frame'] else 'NO'}")
        print(minor_rule())


def select_orf(orfs):
    """
    Single source of truth for ORF selection.

    Important:
    The caller must not call another selection prompt after this function
    returns. The selected ORF is already the object used by downstream RNA
    and protein analysis.
    """
    if not orfs:
        print()
        style.note("No ORF candidates available for selection.")
        return None

    style.header("SELECT ORF")

    style.menu([("Available candidates", [
        (index,
         f"ORF-{orf['rank']:03d} | "
         f"Frame {orf['frame']:+d} | "
         f"{orf['length']} bp | "
         f"{orf['protein_length']} aa | "
         f"{orf['start_codon']} | "
         f"{orf['status']} | "
         f"Heuristic {orf['heuristic_score']['score']}/100")
        for index, orf in enumerate(orfs, start=1)
    ])])

    while True:
        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice.isdigit():
            index = int(choice)

            if 1 <= index <= len(orfs):
                selected = orfs[index - 1]

                style.header("ORF SELECTED")

                style.field('Identifier', f"ORF-{selected['rank']:03d}")
                style.field('Rank', f"#{selected['rank']}")
                style.field('Strand', f"{selected['strand']}")
                style.field('Reading frame', f"{selected['frame']:+d}")
                style.field('Length', f"{selected['length']} bp")
                style.field('Protein length', f"{selected['protein_length']} aa")
                style.field('Start codon', f"{selected['start_codon']}")
                style.field('Status', f"{selected['status']}")
                style.field('Heuristic score', f"{selected['heuristic_score']['score']}/100")
                style.field('CDS quality', f"{selected['cds_quality']['score']}/100")

                style.section('Pipeline Status')
                style.success("ORF Selected")
                style.success("Ready for RNA Analysis")
                style.success("Ready for Protein Analysis")
                style.success("Ready for CRISPR guide design.")

                return selected

        style.invalid_selection()
        style.note("Enter a candidate number or 0 to go back.")


def _orf_settings_menu(settings):
    """Let the user change the ORF search policy. Returns updated settings."""
    settings = dict(settings)
    while True:
        style.header("ORF SEARCH SETTINGS")
        style.menu([
            ("Search policy", [
                (1, f"Minimum ORF length:  {settings['min_length']} nt (stop codon included)"),
                (2, f"Start codons:        {settings['start_codon_policy']}"),
                (3, f"Nested ORFs:         {'reported' if settings['nested'] else 'longest per stop only'}"),
                (4, f"Ranking:             {settings['ranking']}"),
                (5, f"Partial ORFs:        {'included' if settings['include_partial'] else 'excluded'}"),
                (6, f"Genetic code:        {_genetic_code_name(settings.get('genetic_code', 1))}"),
            ]),
            ("Defaults", [
                (7, "Restore Synapse defaults (ORFfinder-style)"),
                (8, "Restore legacy settings (before 2026-10-29)"),
            ]),
        ])
        style.note("Defaults follow NCBI ORFfinder conventions (75 nt, ATG only, one ORF per "
                   "stop). They are screening conventions, not biological thresholds; adjust "
                   "them to your organism and question.")

        choice = style.ask_choice()

        if choice == "0":
            return settings
        if choice == "1":
            raw = style.ask('Minimum length in nt [75]')
            try:
                value = int(raw) if raw else 75
                if value < 3:
                    raise ValueError
                settings["min_length"] = value
            except ValueError:
                style.error("Please enter an integer >= 3.")
        elif choice == "2":
            settings["start_codon_policy"] = (
                "atg_and_alternative"
                if settings["start_codon_policy"] == "atg_only"
                else "atg_only"
            )
        elif choice == "3":
            settings["nested"] = not settings["nested"]
        elif choice == "4":
            settings["ranking"] = (
                "heuristic" if settings["ranking"] == "length" else "length"
            )
        elif choice == "5":
            settings["include_partial"] = not settings["include_partial"]
        elif choice == "7":
            settings.update(
                min_length=DEFAULT_MIN_ORF_LENGTH,
                start_codon_policy=DEFAULT_START_CODON_POLICY,
                nested=DEFAULT_REPORT_NESTED,
                ranking=DEFAULT_RANKING,
                genetic_code=1,
            )
        elif choice == "8":
            settings.update(LEGACY_ORF_SETTINGS)
            settings["genetic_code"] = 1
        elif choice == "6":
            from core.seqtools import available_genetic_codes
            for table in available_genetic_codes():
                print(f"  {table.table_id:>2}  {table.name}")
            raw = style.ask('NCBI table number [1]')
            try:
                value = int(raw) if raw else 1
                _genetic_code_name(value)
                settings["genetic_code"] = value
            except ValueError as exc:
                style.error(f"Genetic code not available: {exc}")
        else:
            style.invalid_selection()


def run_orf_analysis_ui(
    dna,
    include_partial=True,
    min_length=DEFAULT_MIN_ORF_LENGTH,
    start_codon_policy=None,
    nested=DEFAULT_REPORT_NESTED,
    ranking=DEFAULT_RANKING,
    genetic_code=1,
):
    """
    Professional ORF module.

    Workflow:
        Settings -> Discovery -> Candidates -> Compare -> Relationships -> Select

    Returns the selected ORF or None. The selected ORF carries the search
    parameters that produced it under "search_parameters".
    """
    settings = {
        "include_partial": include_partial,
        "min_length": min_length,
        "start_codon_policy": start_codon_policy or DEFAULT_START_CODON_POLICY,
        "nested": nested,
        "ranking": ranking,
        "genetic_code": genetic_code,
    }

    def run_search():
        return analyze_orfs_enhanced(
            dna,
            include_partial=settings["include_partial"],
            min_length=settings["min_length"],
            start_codon_policy=settings["start_codon_policy"],
            nested=settings["nested"],
            ranking=settings["ranking"],
            genetic_code=settings["genetic_code"],
        )

    analysis = run_search()

    while True:
        style.header("ORF ANALYSIS")

        style.field('Sequence length', f"{analysis['sequence_length']} bp")
        style.field('ORF candidates', f"{analysis['candidate_count']}")
        style.field('Search', f">= {settings['min_length']} nt, "
            f"{settings['start_codon_policy']}, "
            f"{'nested' if settings['nested'] else 'longest per stop'}, "
            f"ranked by {settings['ranking']}")

        if settings["genetic_code"] != 1:
            style.field('Genetic code', f"{_genetic_code_name(settings['genetic_code'])}")

        top = analysis.get("top_candidate")

        if top:
            style.field('Best candidate', f"ORF-{top['rank']:03d} "
                f"({top['length']} bp, {top['start_codon']}, {top['status']})")
        else:
            style.field('Best candidate', "None")

        style.menu([
            ("Discovery", [(1, "ORF Discovery"), (2, "ORF Candidates")]),
            ("Evaluation", [(3, "Compare ORFs"), (4, "ORF Relationships")]),
            ("Selection", [(5, "Select ORF")]),
            ("Settings", [(6, "Search Settings")]),
        ])

        choice = style.ask_choice()

        if choice == "1":
            print_orf_report(
                analysis["orfs"],
                min_length=settings["min_length"],
            )
            style.pause()

        elif choice == "2":
            print_orf_summary_card(analysis)
            print_orf_candidates(analysis["orfs"])
            style.pause()

        elif choice == "3":
            print_orf_comparison_table(analysis["orfs"])
            style.pause()

        elif choice == "4":
            print_orf_relationships(analysis["orfs"])
            style.pause()

        elif choice == "5":
            selected = select_orf(analysis["orfs"])
            if selected is not None:
                selected["search_parameters"] = dict(analysis.get("parameters", {}))
                return selected

        elif choice == "6":
            settings = _orf_settings_menu(settings)
            analysis = run_search()

        elif choice == "0":
            return None

        else:
            style.invalid_selection()


# Backward-compatible UI aliases
show_orf_candidates = print_orf_candidates
show_orf_comparison = print_orf_comparison_table
show_orf_relationships = print_orf_relationships