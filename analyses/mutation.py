"""
Synapse - Mutation Analysis
===========================

Professional sequence-level mutation analysis module.

Features
--------
- DNA validation and normalization
- IUPAC ambiguity support
- Reverse-complement utilities
- Needleman-Wunsch global alignment
- Mutation/event detection
- SNP / transition / transversion classification
- Insertions / deletions / complex substitutions
- Mutation statistics
- Coding consequence analysis
- Synonymous / missense / nonsense
- Start-loss / stop-loss
- Frameshift / in-frame indel detection
- Multiple-variant comparison
- Shared / unique mutation analysis
- Sequence quality metrics
- Local mutation context
- Mutation spectra
- Mutation hotspots
- Restriction-site impact
- Primer overlap
- CRISPR target overlap
- Mutation maps
- JSON serialization
- Interactive CLI
- Deterministic regression self-test

This module intentionally performs sequence-level analysis.
Biological pathogenicity, clinical significance, experimental success,
or functional effect must not be inferred from sequence-level descriptors
alone.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
import json
import math
import re
from core.seqtools import IUPAC_COMPLEMENT as _CORE_IUPAC_COMPLEMENT, codon_dict
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)
from ui.style import major_rule, minor_rule
from ui import style


# ============================================================================
# CONSTANTS
# ============================================================================

DNA_ALPHABET = set("ACGTN")

IUPAC_CODES: Dict[str, Set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
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
    "N": {"A", "C", "G", "T"},
}

IUPAC_FROM_BASES = {
    frozenset(value): key
    for key, value in IUPAC_CODES.items()
}

# IUPAC complement table shared via core.seqtools.
IUPAC_COMPLEMENT = dict(_CORE_IUPAC_COMPLEMENT)


# Standard code (NCBI table 1) from core.seqtools.
CODON_TABLE = codon_dict(1)


AA_NAMES = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "Stop",
    "X": "Unknown",
}


PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}


# ============================================================================
# EXCEPTIONS
# ============================================================================

class SequencePromptCancelled(Exception):
    """The user left a sequence prompt with [0]; the caller returns to its menu."""


class MutationAnalysisError(ValueError):
    """Raised when mutation-analysis input or state is invalid."""


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass(frozen=True)
class AlignmentBlock:
    operation: str
    reference_start: int
    reference_end: int
    variant_start: int
    variant_end: int
    reference_sequence: str
    variant_sequence: str


@dataclass
class Mutation:
    mutation_id: int
    mutation_type: str
    position: int
    end_position: int
    reference: str
    alternative: str
    length_change: int
    notation: str

    transition: Optional[bool] = None

    region: str = "Unknown"

    codon_index: Optional[int] = None
    codon_position: Optional[int] = None

    reference_codon: Optional[str] = None
    alternative_codon: Optional[str] = None

    reference_amino_acid: Optional[str] = None
    alternative_amino_acid: Optional[str] = None

    effect: Optional[str] = None
    consequence: Optional[str] = None

    frameshift: bool = False
    premature_stop: bool = False

    notes: List[str] = field(default_factory=list)


@dataclass
class MutationStatistics:
    total_mutations: int = 0
    substitutions: int = 0
    snps: int = 0
    insertions: int = 0
    deletions: int = 0
    complex_changes: int = 0

    transitions: int = 0
    transversions: int = 0

    transition_transversion_ratio: Optional[float] = None

    inserted_bases: int = 0
    deleted_bases: int = 0

    net_length_change: int = 0

    synonymous: int = 0
    missense: int = 0
    nonsense: int = 0
    start_loss: int = 0
    stop_loss: int = 0

    frameshifts: int = 0
    in_frame_indels: int = 0

    unknown_effects: int = 0
    affected_codons: int = 0

    mutation_density_per_kb: float = 0.0


@dataclass
class CodingAnalysis:
    reference_protein: str
    alternative_protein: str

    reference_coding_length: int
    alternative_coding_length: int

    protein_length_change: int

    affected_codons: List[int] = field(default_factory=list)

    synonymous_count: int = 0
    missense_count: int = 0
    nonsense_count: int = 0
    start_loss_count: int = 0
    stop_loss_count: int = 0

    frameshift: bool = False
    frameshift_position: Optional[int] = None

    premature_stop_position: Optional[int] = None

    notes: List[str] = field(default_factory=list)


@dataclass
class MutationAnalysisResult:
    reference_name: str
    variant_name: str

    reference_length: int
    variant_length: int

    mutations: List[Mutation] = field(default_factory=list)

    statistics: MutationStatistics = field(
        default_factory=MutationStatistics
    )

    coding_analysis: Optional[CodingAnalysis] = None

    aligned_reference: str = ""
    aligned_variant: str = ""

    alignment_score: Optional[int] = None

    # Alignment Identity:
    # Exact nucleotide matches among aligned, non-gap positions.
    identity_percent: float = 0.0

    # Sequence Similarity:
    # Edit-distance similarity across the complete ungapped sequences.
    sequence_similarity_percent: float = 0.0

    notes: List[str] = field(default_factory=list)

    @property
    def mutation_count(self) -> int:
        return len(self.mutations)

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(
            self.to_dict(),
            indent=indent,
            sort_keys=False,
        )


@dataclass(frozen=True)
class SequenceCoordinate:
    zero_based: int
    one_based: int

    @classmethod
    def from_one_based(
        cls,
        position: int,
    ) -> "SequenceCoordinate":

        if position < 1:
            raise MutationAnalysisError(
                "One-based coordinates must be >= 1."
            )

        return cls(
            zero_based=position - 1,
            one_based=position,
        )


@dataclass
class SequenceQuality:
    length: int
    gc_percent: float
    n_count: int
    ambiguous_count: int

    valid: bool

    longest_homopolymer: int
    longest_homopolymer_base: Optional[str]

    low_complexity: bool
    complexity_percent: float

    warnings: List[str] = field(default_factory=list)


@dataclass
class MutationContext:
    position: int

    upstream: str
    downstream: str

    reference_context: str
    variant_context: str

    kmer_size: int

    trinucleotide_context: Optional[str] = None


@dataclass
class MutationSpectrum:
    substitutions_12: Dict[str, int] = field(
        default_factory=dict
    )

    substitutions_6: Dict[str, int] = field(
        default_factory=dict
    )

    transitions: int = 0
    transversions: int = 0


@dataclass
class MutationHotspot:
    start: int
    end: int

    mutation_count: int
    density_per_kb: float


@dataclass
class RestrictionSiteImpact:
    enzyme: str
    recognition_sequence: str

    lost_sites: List[int] = field(default_factory=list)
    gained_sites: List[int] = field(default_factory=list)

    unchanged_sites: int = 0


@dataclass
class MutationImpact:
    restriction_sites_lost: int = 0
    restriction_sites_gained: int = 0

    primer_overlap: bool = False
    primer_overlap_ids: List[str] = field(default_factory=list)

    crispr_target_overlap: bool = False
    crispr_overlap_ids: List[str] = field(default_factory=list)

    orf_disrupted: bool = False

    rna_region: Optional[str] = None

    notes: List[str] = field(default_factory=list)


@dataclass
class AdvancedMutationAnalysis:
    reference_quality: Optional[SequenceQuality] = None
    variant_quality: Optional[SequenceQuality] = None

    spectra: Optional[MutationSpectrum] = None

    contexts: List[MutationContext] = field(
        default_factory=list
    )

    hotspots: List[MutationHotspot] = field(
        default_factory=list
    )

    impacts: MutationImpact = field(
        default_factory=MutationImpact
    )

    ambiguity_warnings: List[str] = field(
        default_factory=list
    )

    normalization_notes: List[str] = field(
        default_factory=list
    )


# ============================================================================
# BASIC SEQUENCE FUNCTIONS
# ============================================================================

def normalize_dna(
    sequence: str,
    *,
    allow_n: bool = True,
) -> str:

    if sequence is None:
        raise MutationAnalysisError(
            "DNA sequence cannot be None."
        )

    dna = re.sub(
        r"\s+",
        "",
        str(sequence),
    ).upper()

    if not dna:
        raise MutationAnalysisError(
            "DNA sequence cannot be empty."
        )

    alphabet = (
        DNA_ALPHABET
        if allow_n
        else set("ACGT")
    )

    invalid = sorted(
        set(dna) - alphabet
    )

    if invalid:
        raise MutationAnalysisError(
            "Invalid DNA character(s): "
            + ", ".join(invalid)
        )

    return dna


def normalize_iupac_sequence(
    sequence: str,
) -> str:

    if sequence is None:
        raise MutationAnalysisError(
            "Sequence cannot be None."
        )

    normalized = re.sub(
        r"\s+",
        "",
        str(sequence),
    ).upper()

    if not normalized:
        raise MutationAnalysisError(
            "Sequence cannot be empty."
        )

    invalid = sorted(
        set(normalized) - set(IUPAC_CODES)
    )

    if invalid:
        raise MutationAnalysisError(
            "Invalid IUPAC nucleotide symbol(s): "
            + ", ".join(invalid)
        )

    return normalized


def validate_dna_sequence(
    sequence: str,
    *,
    allow_n: bool = True,
) -> bool:

    normalize_dna(
        sequence,
        allow_n=allow_n,
    )

    return True


def gc_content(sequence: str) -> float:

    dna = normalize_dna(sequence)

    return (
        100.0
        * (
            dna.count("G")
            + dna.count("C")
        )
        / len(dna)
    )


def reverse_complement(
    sequence: str,
) -> str:

    dna = normalize_dna(sequence)

    table = str.maketrans(
        "ACGTN",
        "TGCAN",
    )

    return dna.translate(table)[::-1]


def reverse_complement_iupac(
    sequence: str,
) -> str:

    dna = normalize_iupac_sequence(
        sequence
    )

    return "".join(
        IUPAC_COMPLEMENT[b]
        for b in reversed(dna)
    )


def iupac_matches(
    code_a: str,
    code_b: str,
) -> bool:

    a = IUPAC_CODES.get(
        str(code_a).upper()
    )

    b = IUPAC_CODES.get(
        str(code_b).upper()
    )

    if not a or not b:
        return False

    return bool(a & b)


def iupac_code(
    bases: Iterable[str],
) -> str:

    normalized = frozenset(
        str(base).upper()
        for base in bases
    )

    if not normalized:
        raise MutationAnalysisError(
            "Cannot create an IUPAC code from no bases."
        )

    try:
        return IUPAC_FROM_BASES[
            normalized
        ]

    except KeyError as exc:

        raise MutationAnalysisError(
            "Unsupported base set for IUPAC encoding: "
            + str(sorted(normalized))
        ) from exc


# ============================================================================
# TRANSLATION
# ============================================================================

def translate_dna(
    sequence: str,
    frame: int = 0,
    stop_at_stop: bool = False,
) -> str:

    dna = normalize_iupac_sequence(
        sequence
    )

    if frame not in (0, 1, 2):
        raise MutationAnalysisError(
            "Frame must be 0, 1, or 2."
        )

    protein: List[str] = []

    for i in range(
        frame,
        len(dna) - 2,
        3,
    ):

        codon = dna[i:i + 3]

        if any(
            base not in "ACGT"
            for base in codon
        ):
            aa = "X"

        else:
            aa = CODON_TABLE.get(
                codon,
                "X",
            )

        protein.append(aa)

        if stop_at_stop and aa == "*":
            break

    return "".join(protein)


def amino_acid_name(
    symbol: Optional[str],
) -> Optional[str]:

    if symbol is None:
        return None

    return AA_NAMES.get(
        symbol,
        "Unknown",
    )


# ============================================================================
# ALIGNMENT
# ============================================================================

def _nucleotide_match(
    a: str,
    b: str,
) -> bool:

    return iupac_matches(
        a,
        b,
    )


def global_align(
    reference: str,
    variant: str,
    *,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_score: int = -2,
) -> Tuple[str, str, int]:

    ref = normalize_iupac_sequence(
        reference
    )

    var = normalize_iupac_sequence(
        variant
    )

    rows = len(ref) + 1
    cols = len(var) + 1

    score = [
        [0] * cols
        for _ in range(rows)
    ]

    trace = [
        [""] * cols
        for _ in range(rows)
    ]

    for i in range(1, rows):

        score[i][0] = (
            i * gap_score
        )

        trace[i][0] = "U"

    for j in range(1, cols):

        score[0][j] = (
            j * gap_score
        )

        trace[0][j] = "L"

    for i in range(1, rows):

        for j in range(1, cols):

            diagonal = (
                score[i - 1][j - 1]
                + (
                    match_score
                    if _nucleotide_match(
                        ref[i - 1],
                        var[j - 1],
                    )
                    else mismatch_score
                )
            )

            up = (
                score[i - 1][j]
                + gap_score
            )

            left = (
                score[i][j - 1]
                + gap_score
            )

            # Deterministic tie-breaking:
            # diagonal > up > left
            best = max(
                diagonal,
                up,
                left,
            )

            score[i][j] = best

            if diagonal == best:
                trace[i][j] = "D"
            elif up == best:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

    aligned_ref: List[str] = []
    aligned_var: List[str] = []

    i = len(ref)
    j = len(var)

    while i > 0 or j > 0:

        direction = trace[i][j]

        if direction == "D":

            aligned_ref.append(
                ref[i - 1]
            )

            aligned_var.append(
                var[j - 1]
            )

            i -= 1
            j -= 1

        elif direction == "U":

            aligned_ref.append(
                ref[i - 1]
            )

            aligned_var.append("-")

            i -= 1

        elif direction == "L":

            aligned_ref.append("-")

            aligned_var.append(
                var[j - 1]
            )

            j -= 1

        else:
            raise MutationAnalysisError(
                "Alignment traceback failed."
            )

    return (
        "".join(reversed(aligned_ref)),
        "".join(reversed(aligned_var)),
        score[-1][-1],
    )


def alignment_identity(
    aligned_reference: str,
    aligned_variant: str,
) -> float:
    """
    Calculate Alignment Identity.

    Exact nucleotide matches among aligned positions containing two
    nucleotides. Gap positions are excluded from the denominator.
    """

    if len(aligned_reference) != len(
        aligned_variant
    ):
        raise MutationAnalysisError(
            "Aligned sequences must have equal length."
        )

    if not aligned_reference:
        return 100.0

    comparable = 0
    identical = 0

    for a, b in zip(
        aligned_reference,
        aligned_variant,
    ):

        if a == "-" or b == "-":
            continue

        comparable += 1

        if a == b:
            identical += 1

    if comparable == 0:
        return 0.0

    return (
        100.0
        * identical
        / comparable
    )


def sequence_similarity(
    reference: str,
    variant: str,
) -> float:
    """
    Calculate Sequence Similarity using normalized edit distance.

    Unlike Alignment Identity, this metric does not ignore insertions or
    deletions. A substitution, insertion, or deletion contributes one edit.

    Formula:
        100 * (1 - edit_distance / max(sequence lengths))
    """

    ref = normalize_iupac_sequence(reference)
    var = normalize_iupac_sequence(variant)

    max_length = max(len(ref), len(var))

    if max_length == 0:
        return 100.0

    previous = list(range(len(var) + 1))

    for i, ref_base in enumerate(ref, 1):

        current = [i]

        for j, var_base in enumerate(var, 1):

            substitution_cost = (
                0
                if ref_base == var_base
                else 1
            )

            current.append(
                min(
                    current[-1] + 1,
                    previous[j] + 1,
                    previous[j - 1] + substitution_cost,
                )
            )

        previous = current

    distance = previous[-1]

    return max(
        0.0,
        100.0 * (
            1.0
            - distance / max_length
        ),
    )




# ============================================================================
# ALIGNMENT BLOCKS
# ============================================================================

def build_alignment_blocks(
    aligned_reference: str,
    aligned_variant: str,
) -> List[AlignmentBlock]:

    if len(aligned_reference) != len(
        aligned_variant
    ):
        raise MutationAnalysisError(
            "Aligned sequences must have equal length."
        )

    blocks: List[AlignmentBlock] = []

    ref_pos = 1
    var_pos = 1

    current_operation: Optional[str] = None
    current_ref: List[str] = []
    current_var: List[str] = []

    current_ref_start = 1
    current_var_start = 1

    def flush() -> None:

        nonlocal current_operation
        nonlocal current_ref
        nonlocal current_var

        if current_operation is None:
            return

        ref_seq = "".join(
            current_ref
        )

        var_seq = "".join(
            current_var
        )

        ref_end = (
            current_ref_start
            + len(ref_seq.replace("-", ""))
            - 1
        )

        var_end = (
            current_var_start
            + len(var_seq.replace("-", ""))
            - 1
        )

        if current_operation == "insertion":
            ref_end = (
                current_ref_start - 1
            )

        if current_operation == "deletion":
            var_end = (
                current_var_start - 1
            )

        blocks.append(
            AlignmentBlock(
                operation=current_operation,
                reference_start=current_ref_start,
                reference_end=ref_end,
                variant_start=current_var_start,
                variant_end=var_end,
                reference_sequence=ref_seq,
                variant_sequence=var_seq,
            )
        )

        current_operation = None
        current_ref = []
        current_var = []

    for a, b in zip(
        aligned_reference,
        aligned_variant,
    ):

        if a != "-" and b != "-":

            # Exact equality is a true match. IUPAC compatibility is useful
            # for alignment scoring, but different symbols must remain
            # observable as sequence-level events.
            operation = (
                "match"
                if a == b
                else "substitution"
            )

        elif a == "-" and b != "-":
            operation = "insertion"

        elif a != "-" and b == "-":
            operation = "deletion"

        else:
            raise MutationAnalysisError(
                "Invalid alignment column '--'."
            )

        if operation != current_operation:

            flush()

            current_operation = operation
            current_ref_start = ref_pos
            current_var_start = var_pos

        current_ref.append(a)
        current_var.append(b)

        if a != "-":
            ref_pos += 1

        if b != "-":
            var_pos += 1

    flush()

    return blocks


# ============================================================================
# MUTATION CLASSIFICATION
# ============================================================================

def is_transition(
    reference_base: str,
    alternative_base: str,
) -> bool:

    ref = reference_base.upper()
    alt = alternative_base.upper()

    if ref == alt:
        return False

    return (
        {ref, alt}
        in (
            {"A", "G"},
            {"C", "T"},
        )
    )


def classify_substitution(
    reference_base: str,
    alternative_base: str,
) -> str:

    ref = reference_base.upper()
    alt = alternative_base.upper()

    if ref == alt:
        return "Match"

    if (
        ref not in "ACGT"
        or alt not in "ACGT"
    ):
        return "Ambiguous"

    if is_transition(ref, alt):
        return "Transition"

    return "Transversion"


def _mutation_notation(
    mutation_type: str,
    position: int,
    end_position: int,
    reference: str,
    alternative: str,
) -> str:

    if mutation_type == "SNP":

        return (
            f"{position}"
            f"{reference}>{alternative}"
        )

    if mutation_type == "Insertion":

        return (
            f"{position}ins"
            f"{alternative}"
        )

    if mutation_type == "Deletion":

        if position == end_position:
            return (
                f"{position}del"
                f"{reference}"
            )

        return (
            f"{position}_{end_position}"
            f"del{reference}"
        )

    return (
        f"{position}_{end_position}"
        f"{reference}>{alternative}"
    )


def _make_mutation(
    mutation_id: int,
    mutation_type: str,
    position: int,
    end_position: int,
    reference: str,
    alternative: str,
) -> Mutation:

    ref = reference.replace(
        "-",
        "",
    )

    alt = alternative.replace(
        "-",
        "",
    )

    notation = _mutation_notation(
        mutation_type,
        position,
        end_position,
        ref,
        alt,
    )

    transition = None

    if mutation_type == "SNP":

        transition = is_transition(
            ref,
            alt,
        )

    return Mutation(
        mutation_id=mutation_id,
        mutation_type=mutation_type,
        position=position,
        end_position=end_position,
        reference=ref,
        alternative=alt,
        length_change=(
            len(alt) - len(ref)
        ),
        notation=notation,
        transition=transition,
    )


# ============================================================================
# MUTATION DETECTION
# ============================================================================

def detect_mutations_from_alignment(
    aligned_reference: str,
    aligned_variant: str,
) -> List[Mutation]:

    blocks = build_alignment_blocks(
        aligned_reference,
        aligned_variant,
    )

    mutations: List[Mutation] = []

    mutation_id = 1

    for block in blocks:

        if block.operation == "match":
            continue

        if block.operation == "substitution":

            ref_seq = block.reference_sequence
            var_seq = block.variant_sequence

            for offset, (
                ref_base,
                alt_base,
            ) in enumerate(
                zip(
                    ref_seq,
                    var_seq,
                )
            ):

                # Exact character equality means no mutation. IUPAC
                # compatibility must not suppress a real ambiguity event
                # such as G -> R or T -> Y.
                if ref_base == alt_base:
                    continue

                position = (
                    block.reference_start
                    + offset
                )

                mutation_type = (
                    "SNP"
                    if (
                        ref_base in "ACGT"
                        and alt_base in "ACGT"
                    )
                    else "Substitution"
                )

                mutations.append(
                    _make_mutation(
                        mutation_id,
                        mutation_type,
                        position,
                        position,
                        ref_base,
                        alt_base,
                    )
                )

                mutation_id += 1

        elif block.operation == "insertion":

            mutations.append(
                _make_mutation(
                    mutation_id,
                    "Insertion",
                    block.reference_start,
                    block.reference_start - 1,
                    "",
                    block.variant_sequence,
                )
            )

            mutation_id += 1

        elif block.operation == "deletion":

            mutations.append(
                _make_mutation(
                    mutation_id,
                    "Deletion",
                    block.reference_start,
                    block.reference_end,
                    block.reference_sequence,
                    "",
                )
            )

            mutation_id += 1

    return mutations




def _detect_simple_indel(
    reference: str,
    variant: str,
) -> Optional[List[Mutation]]:
    """
    Detect a single contiguous insertion/deletion using the maximal common
    prefix/suffix representation. This avoids ambiguous global-alignment
    placement inside short repeats for pure indel test cases.
    """
    if len(reference) == len(variant):
        return None

    prefix = 0
    max_prefix = min(len(reference), len(variant))
    while prefix < max_prefix and reference[prefix] == variant[prefix]:
        prefix += 1

    suffix = 0
    max_suffix = min(
        len(reference) - prefix,
        len(variant) - prefix,
    )
    while (
        suffix < max_suffix
        and reference[len(reference) - 1 - suffix]
        == variant[len(variant) - 1 - suffix]
    ):
        suffix += 1

    ref_end = len(reference) - suffix if suffix else len(reference)
    var_end = len(variant) - suffix if suffix else len(variant)

    ref_middle = reference[prefix:ref_end]
    var_middle = variant[prefix:var_end]

    if ref_middle and not var_middle:
        # Deletion starts at the first changed reference base.
        return [
            _make_mutation(
                1,
                "Deletion",
                prefix + 1,
                prefix + len(ref_middle),
                ref_middle,
                "",
            )
        ]

    if var_middle and not ref_middle:
        # Insertion is reported between the bases on either side.
        return [
            _make_mutation(
                1,
                "Insertion",
                prefix + 1,
                prefix,
                "",
                var_middle,
            )
        ]

    return None



def detect_mutations(
    reference: str,
    variant: str,
    *,
    alignment_match_score: int = 2,
    alignment_mismatch_score: int = -1,
    alignment_gap_score: int = -2,
) -> MutationAnalysisResult:
    ref = normalize_iupac_sequence(reference)
    var = normalize_iupac_sequence(variant)

    # A reverse-complement sequence represents the same nucleotide sequence
    # in the opposite orientation. It should not be reported as dozens of
    # substitutions/indels merely because the coordinate orientation changed.
    reverse_complement_match = (
        var == reverse_complement_iupac(ref)
        and var != ref
    )

    aligned_ref, aligned_var, score = global_align(
        ref,
        var,
        match_score=alignment_match_score,
        mismatch_score=alignment_mismatch_score,
        gap_score=alignment_gap_score,
    )

    if reverse_complement_match:
        return MutationAnalysisResult(
            reference_name="Reference",
            variant_name="Variant",
            reference_length=len(ref),
            variant_length=len(var),
            mutations=[],
            statistics=calculate_mutation_statistics(
                [],
                reference_length=len(ref),
            ),
            aligned_reference=aligned_ref,
            aligned_variant=aligned_var,
            alignment_score=score,
            identity_percent=100.0,
            sequence_similarity_percent=100.0,
            notes=[
                "Variant is the reverse complement of the reference sequence.",
                "No sequence-level mutations were reported because the difference is orientation-only.",
            ],
        )

    # For a pure length change, use the maximal common prefix/suffix
    # representation. Needleman-Wunsch can otherwise split one biological
    # indel into multiple events when short repeats are present.
    simple_indel = _detect_simple_indel(ref, var)

    if simple_indel is not None:
        mutations = simple_indel
    else:
        mutations = detect_mutations_from_alignment(
            aligned_ref,
            aligned_var,
        )

    annotate_mutations_with_coding_consequences(
        mutations,
        ref,
        var,
    )

    statistics = calculate_mutation_statistics(
        mutations,
        reference_length=len(ref),
    )

    return MutationAnalysisResult(
        reference_name="Reference",
        variant_name="Variant",
        reference_length=len(ref),
        variant_length=len(var),
        mutations=mutations,
        statistics=statistics,
        aligned_reference=aligned_ref,
        aligned_variant=aligned_var,
        alignment_score=score,
        identity_percent=alignment_identity(
            aligned_ref,
            aligned_var,
        ),
        sequence_similarity_percent=sequence_similarity(
            ref,
            var,
        ),
    )


def classify_coding_substitution(
    reference_codon: str,
    alternative_codon: str,
    codon_index: Optional[int] = None,
) -> str:
    """Classify a codon substitution (standard genetic code).

    ``codon_index`` is the 1-based codon number within the CDS. "Start-loss"
    applies only to codon 1 (the initiator); a Met->X change at any other
    codon is an ordinary missense change (audit P0-6). When ``codon_index``
    is None the position is unknown and Met->X is reported as "Missense".
    """

    ref = normalize_dna(
        reference_codon,
        allow_n=False,
    )

    alt = normalize_dna(
        alternative_codon,
        allow_n=False,
    )

    if len(ref) != 3 or len(alt) != 3:
        raise MutationAnalysisError(
            "Codons must contain exactly 3 bases."
        )

    ref_aa = CODON_TABLE.get(
        ref,
        "X",
    )

    alt_aa = CODON_TABLE.get(
        alt,
        "X",
    )

    if ref_aa == alt_aa:
        return "Synonymous"

    if ref_aa != "*" and alt_aa == "*":
        return "Nonsense"

    if ref_aa == "*" and alt_aa != "*":
        return "Stop-loss"

    if codon_index == 1 and ref_aa == "M" and alt_aa != "M":
        return "Start-loss"

    return "Missense"


def _translate_canonical(
    sequence: str,
) -> str:

    dna = normalize_dna(
        sequence,
        allow_n=False,
    )

    protein = []

    for i in range(
        0,
        len(dna) - 2,
        3,
    ):

        protein.append(
            CODON_TABLE.get(
                dna[i:i + 3],
                "X",
            )
        )

    return "".join(protein)



def _first_difference_position(reference: str, alternative: str) -> Optional[int]:
    """Return the first one-based coordinate where two sequences differ."""
    for index, (ref_base, alt_base) in enumerate(zip(reference, alternative), 1):
        if ref_base != alt_base:
            return index
    if len(reference) != len(alternative):
        return min(len(reference), len(alternative)) + 1
    return None


def analyze_coding_consequences(
    reference_cds: str,
    alternative_cds: str,
) -> CodingAnalysis:
    """
    Analyze coding consequences without treating an indel as a cascade of
    independent downstream amino-acid substitutions.

    Rules:
    - Equal coding length: compare complete codons one-to-one.
    - Length difference divisible by 3: classify the event as an in-frame
      indel; do not count shifted downstream codons as missense mutations.
    - Length difference not divisible by 3: classify the event as a frameshift;
      downstream amino-acid differences are not counted as individual
      missense mutations.
    """
    ref = normalize_iupac_sequence(reference_cds)
    alt = normalize_iupac_sequence(alternative_cds)

    if len(ref) < 3:
        raise MutationAnalysisError(
            "Reference CDS must contain at least one complete codon."
        )

    # Ambiguous IUPAC symbols are valid nucleotide input, but they do not
    # define a unique amino acid. Preserve the sequence-level information
    # without inventing a protein consequence.
    canonical_ref = set(ref) <= set("ACGT")
    canonical_alt = set(alt) <= set("ACGT")

    # A reverse-complement sequence represents the same nucleotide sequence
    # in the opposite orientation. For coding-consequence analysis, normalize
    # that orientation before comparing codons. Otherwise a pure orientation
    # change would be incorrectly reported as dozens of missense changes.
    ref_protein = (
        _translate_canonical(ref)
        if canonical_ref
        else ""
    )
    alt_protein = (
        _translate_canonical(alt)
        if canonical_alt
        else ""
    )

    reverse_complement_match = (
        len(ref) == len(alt)
        and alt == reverse_complement_iupac(ref)
        and alt != ref
    )

    if reverse_complement_match and canonical_ref and canonical_alt:
        coding = CodingAnalysis(
            reference_protein=ref_protein,
            alternative_protein=ref_protein,
            reference_coding_length=len(ref),
            alternative_coding_length=len(alt),
            protein_length_change=0,
            notes=[
                "Alternative CDS is the reverse complement of the reference CDS.",
                "No coding consequences were reported because the difference is orientation-only.",
            ],
        )
        return coding

    if not (canonical_ref and canonical_alt):
        return CodingAnalysis(
            reference_protein="",
            alternative_protein="",
            reference_coding_length=len(ref),
            alternative_coding_length=len(alt),
            protein_length_change=0,
            notes=[
                "Coding consequence: Unknown.",
                "Exact amino-acid consequences were not inferred because "
                "ambiguous IUPAC bases are present.",
                "A codon containing an ambiguous nucleotide cannot be assigned "
                "a unique amino acid without resolving the ambiguity.",
            ],
        )

    coding = CodingAnalysis(
        reference_protein=ref_protein,
        alternative_protein=alt_protein,
        reference_coding_length=len(ref),
        alternative_coding_length=len(alt),
        protein_length_change=len(alt_protein) - len(ref_protein),
    )

    delta = len(alt) - len(ref)

    # No length change: true codon-by-codon substitution analysis.
    if delta == 0:
        complete_codons = min(len(ref), len(alt)) // 3

        for codon_zero in range(complete_codons):
            start = codon_zero * 3
            ref_codon = ref[start:start + 3]
            alt_codon = alt[start:start + 3]

            if ref_codon == alt_codon:
                continue

            effect = classify_coding_substitution(
                ref_codon,
                alt_codon,
                codon_index=codon_zero + 1,
            )

            coding.affected_codons.append(codon_zero + 1)

            if effect == "Synonymous":
                coding.synonymous_count += 1
            elif effect == "Missense":
                coding.missense_count += 1
            elif effect == "Nonsense":
                coding.nonsense_count += 1
                if coding.premature_stop_position is None:
                    coding.premature_stop_position = codon_zero + 1
            elif effect == "Start-loss":
                coding.start_loss_count += 1
            elif effect == "Stop-loss":
                coding.stop_loss_count += 1

        return coding

    # Any non-zero length change is structurally an indel event here.
    first_diff = _first_difference_position(ref, alt)
    if first_diff is None:
        first_diff = min(len(ref), len(alt)) + 1

    codon_index = (first_diff - 1) // 3 + 1
    coding.affected_codons.append(codon_index)

    if delta % 3 == 0:
        coding.notes.append(
            f"In-frame indel: coding length changed by {delta:+d} bp."
        )
        return coding

    coding.frameshift = True
    coding.frameshift_position = first_diff
    coding.notes.append(
        f"Frameshift: coding length changed by {delta:+d} bp."
    )

    # A stop after a frameshift is descriptive only; it is not counted as a
    # separate nonsense substitution.
    frame_aa_index = max(0, (first_diff - 1) // 3)
    for aa_index in range(frame_aa_index, len(alt_protein)):
        if alt_protein[aa_index] == "*":
            coding.premature_stop_position = aa_index + 1
            break

    return coding



def annotate_mutations_with_coding_consequences(
    mutations: List[Mutation],
    reference_cds: str,
    alternative_cds: str,
) -> CodingAnalysis:
    """
    Add coding annotations to detected events.

    Important:
    - SNPs/substitutions are classified codon-by-codon only when both
      sequences have identical length.
    - Insertions/deletions are classified structurally as either in-frame
      or frameshift events.
    - IUPAC ambiguity is accepted at sequence level, but an exact amino-acid
      consequence is not invented for ambiguous codons.
    """
    ref = normalize_iupac_sequence(reference_cds)
    alt = normalize_iupac_sequence(alternative_cds)

    canonical_ref = set(ref) <= set("ACGT")
    canonical_alt = set(alt) <= set("ACGT")

    # Structural annotations are always safe and do not require canonical CDS.
    for mutation in mutations:
        mutation.effect = None
        mutation.consequence = None
        mutation.frameshift = False
        mutation.premature_stop = False

        mutation.codon_index = (
            (max(1, mutation.position) - 1) // 3 + 1
        )
        mutation.codon_position = (
            (max(1, mutation.position) - 1) % 3 + 1
        )

        if mutation.mutation_type in {"Insertion", "Deletion"}:
            if mutation.length_change % 3 == 0:
                mutation.effect = "In-frame indel"
                mutation.consequence = "In-frame indel"
            else:
                mutation.effect = "Frameshift"
                mutation.consequence = "Frameshift"
                mutation.frameshift = True
            continue

        if mutation.mutation_type == "Complex":
            if mutation.length_change % 3 == 0:
                mutation.effect = "In-frame indel"
                mutation.consequence = "In-frame indel"
            else:
                mutation.effect = "Frameshift"
                mutation.consequence = "Frameshift"
                mutation.frameshift = True

    # Ambiguous bases: preserve sequence-level detection but do not fabricate
    # exact protein consequences.
    if not (canonical_ref and canonical_alt):
        coding = CodingAnalysis(
            reference_protein="",
            alternative_protein="",
            reference_coding_length=len(ref),
            alternative_coding_length=len(alt),
            protein_length_change=0,
            notes=[
                "Coding consequence: Unknown.",
                "Exact amino-acid consequences were not inferred because "
                "ambiguous IUPAC bases are present.",
                "A codon containing an ambiguous nucleotide cannot be assigned "
                "a unique amino acid without resolving the ambiguity.",
            ],
        )
        for mutation in mutations:
            if mutation.mutation_type in {"SNP", "Substitution"}:
                mutation.effect = "Unknown"
                mutation.consequence = "Unknown"
                mutation.notes.append(
                    "Ambiguous IUPAC base prevents exact codon consequence."
                )
        return coding

    coding = analyze_coding_consequences(
        ref,
        alt,
    )

    # Only equal-length sequences support direct SNP codon comparison.
    # For indel-containing comparisons, this deliberately avoids the old
    # downstream-missense cascade.
    if len(ref) != len(alt):
        return coding

    complete_codons = len(ref) // 3

    for mutation in mutations:
        if mutation.mutation_type not in {"SNP", "Substitution"}:
            continue

        codon_zero = (mutation.position - 1) // 3
        codon_start = codon_zero * 3

        if codon_zero >= complete_codons:
            mutation.effect = "Unknown"
            mutation.consequence = "Unknown"
            mutation.notes.append(
                "Mutation lies in an incomplete terminal codon."
            )
            continue

        ref_codon = ref[codon_start:codon_start + 3]
        alt_codon = alt[codon_start:codon_start + 3]

        if len(ref_codon) != 3 or len(alt_codon) != 3:
            mutation.effect = "Unknown"
            mutation.consequence = "Unknown"
            continue

        mutation.codon_index = codon_zero + 1
        mutation.codon_position = (mutation.position - 1) % 3 + 1
        mutation.reference_codon = ref_codon
        mutation.alternative_codon = alt_codon
        mutation.reference_amino_acid = CODON_TABLE[ref_codon]
        mutation.alternative_amino_acid = CODON_TABLE[alt_codon]

        effect = classify_coding_substitution(
            ref_codon,
            alt_codon,
            codon_index=codon_zero + 1,
        )
        mutation.effect = effect
        mutation.consequence = effect

        if effect == "Nonsense":
            mutation.premature_stop = True

    return coding


def calculate_mutation_statistics(
    mutations: Sequence[Mutation],
    *,
    reference_length: int,
) -> MutationStatistics:

    if reference_length < 1:
        raise MutationAnalysisError(
            "Reference length must be positive."
        )

    stats = MutationStatistics()

    stats.total_mutations = len(
        mutations
    )

    for mutation in mutations:

        if mutation.mutation_type in {
            "SNP",
            "Substitution",
        }:

            stats.substitutions += 1

        if mutation.mutation_type == "SNP":

            stats.snps += 1

            if mutation.transition is True:
                stats.transitions += 1

            elif mutation.transition is False:
                stats.transversions += 1

        elif mutation.mutation_type == "Insertion":

            stats.insertions += 1

            stats.inserted_bases += len(
                mutation.alternative
            )

            if (
                len(mutation.alternative)
                % 3
                == 0
            ):
                stats.in_frame_indels += 1

            else:
                stats.frameshifts += 1

        elif mutation.mutation_type == "Deletion":

            stats.deletions += 1

            stats.deleted_bases += len(
                mutation.reference
            )

            if (
                len(mutation.reference)
                % 3
                == 0
            ):
                stats.in_frame_indels += 1

            else:
                stats.frameshifts += 1

        elif mutation.mutation_type == "Complex":

            stats.complex_changes += 1

        effect = mutation.effect

        if effect == "Synonymous":
            stats.synonymous += 1

        elif effect == "Missense":
            stats.missense += 1

        elif effect == "Nonsense":
            stats.nonsense += 1

        elif effect == "Start-loss":
            stats.start_loss += 1

        elif effect == "Stop-loss":
            stats.stop_loss += 1

        elif effect in {
            None,
            "Unknown",
        }:
            stats.unknown_effects += 1

    stats.net_length_change = (
        stats.inserted_bases
        - stats.deleted_bases
    )

    if stats.transversions == 0:

        if stats.transitions > 0:
            stats.transition_transversion_ratio = (
                math.inf
            )
        else:
            stats.transition_transversion_ratio = (
                None
            )

    else:

        stats.transition_transversion_ratio = (
            stats.transitions
            / stats.transversions
        )

    stats.affected_codons = len({
        mutation.codon_index
        for mutation in mutations
        if mutation.codon_index is not None
    })

    stats.mutation_density_per_kb = (
        stats.total_mutations
        / reference_length
        * 1000.0
    )

    return stats


# ============================================================================
# COMPARISON
# ============================================================================

def compare_sequences(
    reference: str,
    variant: str,
    *,
    reference_name: str = "Reference",
    variant_name: str = "Variant",
    alignment_match_score: int = 2,
    alignment_mismatch_score: int = -1,
    alignment_gap_score: int = -2,
) -> MutationAnalysisResult:

    result = detect_mutations(
        reference,
        variant,
        alignment_match_score=alignment_match_score,
        alignment_mismatch_score=alignment_mismatch_score,
        alignment_gap_score=alignment_gap_score,
    )

    result.reference_name = (
        reference_name
    )

    result.variant_name = (
        variant_name
    )

    longest = max(result.reference_length, result.variant_length)
    if longest > MAX_RECOMMENDED_ALIGNMENT_LENGTH:
        result.notes.append(
            f"Sequences up to {longest} nt were aligned globally "
            f"(recommended maximum {MAX_RECOMMENDED_ALIGNMENT_LENGTH} nt): the "
            "global alignment fills a dynamic-programming matrix of "
            "(n+1)x(m+1) cells, measured at roughly 45 bytes per cell, so "
            "about 180 MB is needed at 2,000 nt and about 4.5 GB at 10,000 nt. "
            "For genome-scale comparisons use VCF-based variant analysis."
        )

    return result


def compare_multiple_variants(
    reference: str,
    variants: Sequence[str],
    *,
    reference_name: str = "Reference",
) -> Dict[str, MutationAnalysisResult]:

    results: Dict[
        str,
        MutationAnalysisResult,
    ] = {}

    for index, variant in enumerate(
        variants,
        1,
    ):

        name = (
            f"Variant {index}"
        )

        results[name] = compare_sequences(
            reference,
            variant,
            reference_name=reference_name,
            variant_name=name,
        )

    return results


def _mutation_key(
    mutation: Mutation,
) -> Tuple:

    return (
        mutation.mutation_type,
        mutation.position,
        mutation.end_position,
        mutation.reference,
        mutation.alternative,
    )


def shared_mutations(
    results: Mapping[
        str,
        MutationAnalysisResult,
    ],
) -> List[Mutation]:

    if not results:
        return []

    mutation_sets = [
        {
            _mutation_key(m)
            for m in result.mutations
        }
        for result in results.values()
    ]

    shared_keys = set.intersection(
        *mutation_sets
    )

    first_result = next(
        iter(results.values())
    )

    return [
        mutation
        for mutation in first_result.mutations
        if _mutation_key(mutation)
        in shared_keys
    ]


def unique_mutations(
    result: MutationAnalysisResult,
    other_results: Sequence[
        MutationAnalysisResult
    ],
) -> List[Mutation]:

    other_keys = {
        _mutation_key(mutation)
        for other in other_results
        for mutation in other.mutations
    }

    return [
        mutation
        for mutation in result.mutations
        if _mutation_key(mutation)
        not in other_keys
    ]


# ============================================================================
# REGIONS
# ============================================================================

def annotate_regions(
    mutations: Sequence[Mutation],
    regions: Mapping[str, Tuple[int, int]],
) -> List[Mutation]:

    for mutation in mutations:

        for name, (
            start,
            end,
        ) in regions.items():

            if (
                mutation.position >= start
                and mutation.position <= end
            ):

                mutation.region = name
                break

    return list(mutations)


def annotate_sequence_regions(
    mutations: Sequence[Mutation],
    regions: Mapping[str, Tuple[int, int]],
) -> List[Mutation]:

    return annotate_regions(
        mutations,
        regions,
    )


# ============================================================================
# QUALITY
# ============================================================================

def sequence_complexity(
    sequence: str,
) -> float:

    dna = normalize_iupac_sequence(
        sequence
    )

    canonical = [
        base
        for base in dna
        if base in "ACGT"
    ]

    if not canonical:
        return 0.0

    counts = Counter(
        canonical
    )

    total = len(canonical)

    entropy = 0.0

    for count in counts.values():

        p = count / total

        entropy -= (
            p * math.log2(p)
        )

    # Maximum Shannon entropy for four bases = 2 bits.
    return (
        entropy / 2.0
        * 100.0
    )


def sequence_quality(
    sequence: str,
    *,
    low_complexity_threshold: float = 75.0,
) -> SequenceQuality:

    dna = normalize_iupac_sequence(
        sequence
    )

    counts = Counter(dna)

    gc = (
        100.0
        * (
            counts["G"]
            + counts["C"]
        )
        / len(dna)
    )

    n_count = counts["N"]

    ambiguous_count = sum(
        counts[base]
        for base in dna
        if base not in "ACGT"
    )

    longest = 1
    longest_base = dna[0]

    current = 1

    for i in range(
        1,
        len(dna),
    ):

        if dna[i] == dna[i - 1]:
            current += 1
        else:
            current = 1

        if current > longest:

            longest = current
            longest_base = dna[i]

    complexity = sequence_complexity(
        dna
    )

    warnings: List[str] = []

    if ambiguous_count:
        warnings.append(
            f"{ambiguous_count} ambiguous/IUPAC base(s) detected."
        )

    if longest >= 6:
        warnings.append(
            f"Long homopolymer detected: "
            f"{longest}{longest_base}."
        )

    if complexity < low_complexity_threshold:
        warnings.append(
            "Sequence shows reduced nucleotide complexity."
        )

    return SequenceQuality(
        length=len(dna),
        gc_percent=gc,
        n_count=n_count,
        ambiguous_count=ambiguous_count,
        valid=True,
        longest_homopolymer=longest,
        longest_homopolymer_base=longest_base,
        low_complexity=(
            complexity
            < low_complexity_threshold
        ),
        complexity_percent=complexity,
        warnings=warnings,
    )


# ============================================================================
# CONTEXT / SPECTRUM
# ============================================================================

def extract_mutation_context(
    mutation: Mutation,
    reference: str,
    variant: str,
    *,
    kmer_size: int = 5,
) -> MutationContext:

    if kmer_size < 1 or kmer_size % 2 == 0:
        raise MutationAnalysisError(
            "kmer_size must be a positive odd integer."
        )

    ref = normalize_iupac_sequence(
        reference
    )

    var = normalize_iupac_sequence(
        variant
    )

    radius = kmer_size // 2

    center = (
        mutation.position - 1
    )

    start = max(
        0,
        center - radius,
    )

    end = min(
        len(ref),
        center + radius + 1,
    )

    ref_context = ref[
        start:end
    ]

    if mutation.mutation_type == "Insertion":

        variant_center = min(
            len(var),
            mutation.position - 1,
        )

    else:

        variant_center = min(
            len(var) - 1,
            center,
        )

    variant_start = max(
        0,
        variant_center - radius,
    )

    variant_end = min(
        len(var),
        variant_center + radius + 1,
    )

    var_context = var[
        variant_start:variant_end
    ]

    upstream = ref[
        max(
            0,
            center - radius,
        ):
        center
    ]

    downstream = ref[
        center + 1:
        min(
            len(ref),
            center + radius + 1,
        )
    ]

    trinucleotide = None

    if (
        mutation.mutation_type == "SNP"
        and len(ref_context) >= 3
    ):

        local = center - start

        if (
            local >= 1
            and local + 1 < len(ref_context)
        ):

            trinucleotide = (
                ref_context[
                    local - 1:
                    local + 2
                ]
            )

    return MutationContext(
        position=mutation.position,
        upstream=upstream,
        downstream=downstream,
        reference_context=ref_context,
        variant_context=var_context,
        kmer_size=kmer_size,
        trinucleotide_context=trinucleotide,
    )


def build_mutation_spectrum(
    mutations: Sequence[Mutation],
) -> MutationSpectrum:

    spectrum = MutationSpectrum()

    substitution_types = (
        "A>G", "G>A",
        "C>T", "T>C",
        "A>C", "C>A",
        "A>T", "T>A",
        "C>G", "G>C",
        "G>T", "T>G",
    )

    spectrum.substitutions_12 = {
        key: 0
        for key in substitution_types
    }

    for mutation in mutations:

        if mutation.mutation_type != "SNP":
            continue

        ref = mutation.reference
        alt = mutation.alternative

        key = f"{ref}>{alt}"

        if key in spectrum.substitutions_12:
            spectrum.substitutions_12[key] += 1

        if mutation.transition is True:
            spectrum.transitions += 1

        elif mutation.transition is False:
            spectrum.transversions += 1

    # Six-class pyrimidine-oriented spectrum.
    spectrum.substitutions_6 = {
        "C>T": 0,
        "C>A": 0,
        "C>G": 0,
        "T>C": 0,
        "T>A": 0,
        "T>G": 0,
    }

    for mutation in mutations:

        if mutation.mutation_type != "SNP":
            continue

        ref = mutation.reference
        alt = mutation.alternative

        # Orient purine substitutions through complement.
        if ref in "AG":

            ref = IUPAC_COMPLEMENT[
                ref
            ]

            alt = IUPAC_COMPLEMENT[
                alt
            ]

        key = f"{ref}>{alt}"

        if key in spectrum.substitutions_6:
            spectrum.substitutions_6[key] += 1

    return spectrum


def find_mutation_hotspots(
    mutations: Sequence[Mutation],
    sequence_length: int,
    *,
    window_size: int = 20,
    minimum_mutations: int = 2,
) -> List[MutationHotspot]:

    if sequence_length < 1:
        raise MutationAnalysisError(
            "sequence_length must be positive."
        )

    if window_size < 1:
        raise MutationAnalysisError(
            "window_size must be positive."
        )

    positions = sorted(
        mutation.position
        for mutation in mutations
    )

    hotspots: List[
        MutationHotspot
    ] = []

    if not positions:
        return hotspots

    for start in range(
        1,
        sequence_length + 1,
        max(1, window_size // 2),
    ):

        end = min(
            sequence_length,
            start + window_size - 1,
        )

        count = sum(
            start <= position <= end
            for position in positions
        )

        if count >= minimum_mutations:

            density = (
                count
                / (end - start + 1)
                * 1000.0
            )

            hotspots.append(
                MutationHotspot(
                    start=start,
                    end=end,
                    mutation_count=count,
                    density_per_kb=density,
                )
            )

    # Merge overlapping windows.
    merged: List[
        MutationHotspot
    ] = []

    for hotspot in hotspots:

        if not merged:

            merged.append(hotspot)
            continue

        previous = merged[-1]

        if hotspot.start <= previous.end:

            new_start = previous.start
            new_end = max(
                previous.end,
                hotspot.end,
            )

            new_count = sum(
                start <= position <= new_end
                for position in positions
                for start in [new_start]
            )

            merged[-1] = MutationHotspot(
                start=new_start,
                end=new_end,
                mutation_count=new_count,
                density_per_kb=(
                    new_count
                    / (
                        new_end
                        - new_start
                        + 1
                    )
                    * 1000.0
                ),
            )

        else:

            merged.append(hotspot)

    return merged


# ============================================================================
# INDEL NORMALIZATION
# ============================================================================

def normalize_indel_left(
    reference: str,
    variant: str,
    position: int,
    reference_allele: str,
    alternative_allele: str,
) -> Tuple[
    int,
    str,
    str,
]:

    ref = normalize_iupac_sequence(
        reference
    )

    if position < 1:
        raise MutationAnalysisError(
            "Position must be >= 1."
        )

    if (
        reference_allele
        and alternative_allele
    ):

        return (
            position,
            reference_allele,
            alternative_allele,
        )

    # Conservative left-normalization for simple indels.
    pos = position
    ref_allele = reference_allele
    alt_allele = alternative_allele

    while pos > 1:

        if not ref_allele and alt_allele:

            previous = ref[pos - 2]

            if alt_allele[-1] == previous:
                alt_allele = (
                    previous
                    + alt_allele[:-1]
                )
                pos -= 1
            else:
                break

        elif ref_allele and not alt_allele:

            previous = ref[pos - 2]

            if ref_allele[-1] == previous:
                ref_allele = (
                    previous
                    + ref_allele[:-1]
                )
                pos -= 1
            else:
                break

        else:
            break

    return (
        pos,
        ref_allele,
        alt_allele,
    )


# ============================================================================
# CODING-FRAME ASSUMPTIONS (audit P1-10)
# ============================================================================

# Global alignment needs O(n*m) time and memory. Above this length a warning
# is attached; the analysis still runs. The figure is a practical limit for an
# interactive pure-Python tool, not a scientific threshold.
MAX_RECOMMENDED_ALIGNMENT_LENGTH = 2000

CODING_FRAME_ASSUMPTION = (
    "Coding consequences assume that position 1 of the reference is the first "
    "base of codon 1 (CDS coordinates). For a genomic or non-coding sequence "
    "this assumption is wrong and amino-acid effects are meaningless."
)


def assess_coding_frame(cds: str) -> List[str]:
    """Return warnings when a sequence does not look like a CDS in frame +1.

    These checks are plausibility checks, not CDS validation: a sequence can
    pass all of them and still not be the real coding sequence.
    """
    seq = normalize_iupac_sequence(cds)
    warnings: List[str] = []
    if len(seq) % 3:
        warnings.append(
            f"Length {len(seq)} is not a multiple of 3; the last "
            f"{len(seq) % 3} base(s) do not form a codon."
        )
    if len(seq) >= 3 and seq[:3] not in {"ATG", "GTG", "TTG", "CTG"}:
        warnings.append(
            f"Does not begin with a start codon ({seq[:3]}); "
            "frame +1 may not be the coding frame."
        )
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    internal_stops = [
        index + 1
        for index, codon in enumerate(codons[:-1])
        if codon in {"TAA", "TAG", "TGA"}
    ]
    if internal_stops:
        shown = ", ".join(str(i) for i in internal_stops[:5])
        warnings.append(
            f"{len(internal_stops)} in-frame stop codon(s) before the end "
            f"(codon {shown}{', ...' if len(internal_stops) > 5 else ''}); "
            "reading frame +1 is unlikely to be the coding frame."
        )
    if codons and codons[-1] not in {"TAA", "TAG", "TGA"}:
        warnings.append("Does not end with a stop codon (possibly partial CDS).")
    return warnings


def orf_reference_cds(orf: Any, reference: str, variant: str) -> Tuple[str, Optional[str], str]:
    """Derive reference/variant CDS from a pipeline ORF.

    Returns (reference_cds, variant_cds or None, description). The variant
    CDS is only derived when the variant has the reference length (pure
    substitutions), because only then do ORF coordinates carry over directly.
    """
    from core.orf_context import resolve_orf_region

    region = resolve_orf_region(orf, reference)
    if not region.validated or region.start_1based is None:
        raise MutationAnalysisError(
            "The selected ORF could not be located on the reference sequence."
        )
    ref_cds = region.coding_sequence
    description = (
        f"{region.name}: {region.start_1based}-{region.end_1based} "
        f"({region.strand} strand)"
    )
    var_cds = None
    if len(variant) == len(reference):
        forward = variant[region.start_1based - 1:region.end_1based]
        var_cds = reverse_complement_iupac(forward) if region.strand == "-" else forward
    return ref_cds, var_cds, description


# ============================================================================
# HGVS-LIKE NOTATION
# ============================================================================

def mutation_hgvs_c(
    mutation: Mutation,
) -> str:
    """HGVS-like coding notation ("c." + Synapse notation).

    Not validated against the HGVS specification or a tool such as
    Mutalyzer. Position 1 is the first base of the supplied sequence, which
    equals HGVS c. numbering only if that sequence is the CDS starting at
    the A of the start codon.
    """

    return (
        "c."
        + mutation.notation
    )


def mutation_hgvs_p(
    mutation: Mutation,
) -> Optional[str]:
    """Predicted protein consequence in HGVS style, one-letter codes.

    The consequence is inferred from the DNA, not observed, so it is
    written in parentheses (HGVS recommendations, den Dunnen et al. 2016):
    p.(R12C) missense, p.(R12=) synonymous, p.(W24*) nonsense, p.(M1?)
    for a change of the initiator codon (the effect on initiation is
    unknown) and p.(*110Qext*?) for a lost stop codon whose new stop is not
    determined. Positions count codons of the supplied sequence, which is
    HGVS p. numbering only if that sequence is the CDS. Not checked
    against an HGVS validator. None when no codon change is known.
    """

    if (
        mutation.reference_amino_acid
        is None
        or mutation.alternative_amino_acid
        is None
        or mutation.codon_index is None
    ):
        return None

    ref = mutation.reference_amino_acid
    alt = mutation.alternative_amino_acid
    index = mutation.codon_index

    if ref == alt:
        return f"p.({ref}{index}=)"

    if index == 1 and ref == "M":
        return "p.(M1?)"

    if ref == "*":
        return f"p.(*{index}{alt}ext*?)"

    return f"p.({ref}{index}{alt})"


def compare_proteins(
    reference_protein: str,
    alternative_protein: str,
) -> Dict[str, Any]:

    ref = str(
        reference_protein
    )

    alt = str(
        alternative_protein
    )

    length = min(
        len(ref),
        len(alt),
    )

    mismatches = []

    for i in range(length):

        if ref[i] != alt[i]:

            mismatches.append(
                {
                    "position": i + 1,
                    "reference": ref[i],
                    "alternative": alt[i],
                }
            )

    return {
        "reference_length": len(ref),
        "alternative_length": len(alt),
        "length_change": len(alt) - len(ref),
        "identity_percent": (
            100.0
            if max(len(ref), len(alt)) == 0
            else (
                100.0
                * (
                    len(ref)
                    - len(mismatches)
                )
                / max(
                    len(ref),
                    len(alt),
                )
            )
        ),
        "mismatches": mismatches,
    }


# ============================================================================
# RESTRICTION / PRIMER / CRISPR IMPACT
# ============================================================================

def _find_all_occurrences(
    sequence: str,
    motif: str,
) -> List[int]:

    seq = normalize_iupac_sequence(
        sequence
    )

    motif = normalize_iupac_sequence(
        motif
    )

    hits = []

    if len(motif) > len(seq):
        return hits

    for i in range(
        0,
        len(seq) - len(motif) + 1,
    ):

        if all(
            iupac_matches(
                seq[i + j],
                motif[j],
            )
            for j in range(len(motif))
        ):
            hits.append(i + 1)

    return hits


def restriction_site_impact(
    reference: str,
    variant: str,
    enzymes: Mapping[str, str],
) -> List[RestrictionSiteImpact]:

    ref = normalize_iupac_sequence(
        reference
    )

    var = normalize_iupac_sequence(
        variant
    )

    results = []

    for enzyme, motif in enzymes.items():

        ref_sites = _find_all_occurrences(
            ref,
            motif,
        )

        var_sites = _find_all_occurrences(
            var,
            motif,
        )

        lost = [
            position
            for position in ref_sites
            if position not in var_sites
        ]

        gained = [
            position
            for position in var_sites
            if position not in ref_sites
        ]

        unchanged = len(
            set(ref_sites)
            & set(var_sites)
        )

        results.append(
            RestrictionSiteImpact(
                enzyme=str(enzyme),
                recognition_sequence=str(
                    motif
                ),
                lost_sites=lost,
                gained_sites=gained,
                unchanged_sites=unchanged,
            )
        )

    return results


def _interval_overlap(
    mutation: Mutation,
    start: int,
    end: int,
) -> bool:

    mutation_start = mutation.position

    mutation_end = max(
        mutation.position,
        mutation.end_position,
    )

    return not (
        mutation_end < start
        or mutation_start > end
    )


def primer_mutation_overlap(
    mutations: Sequence[Mutation],
    primers: Sequence[
        Tuple[str, int, int]
    ],
) -> MutationImpact:

    impact = MutationImpact()

    for mutation in mutations:

        for primer_id, start, end in primers:

            if _interval_overlap(
                mutation,
                start,
                end,
            ):

                if primer_id not in (
                    impact.primer_overlap_ids
                ):

                    impact.primer_overlap_ids.append(
                        str(primer_id)
                    )

    impact.primer_overlap = bool(
        impact.primer_overlap_ids
    )

    return impact


def crispr_mutation_overlap(
    mutations: Sequence[Mutation],
    targets: Sequence[
        Tuple[str, int, int]
    ],
) -> MutationImpact:

    impact = MutationImpact()

    for mutation in mutations:

        for target_id, start, end in targets:

            if _interval_overlap(
                mutation,
                start,
                end,
            ):

                if target_id not in (
                    impact.crispr_overlap_ids
                ):

                    impact.crispr_overlap_ids.append(
                        str(target_id)
                    )

    impact.crispr_target_overlap = bool(
        impact.crispr_overlap_ids
    )

    return impact


# ============================================================================
# RESULT VALIDATION
# ============================================================================

def validate_mutation_result(
    result: MutationAnalysisResult,
) -> bool:

    if result.reference_length < 1:
        raise MutationAnalysisError(
            "Reference length is invalid."
        )

    if result.variant_length < 1:
        raise MutationAnalysisError(
            "Variant length is invalid."
        )

    if len(result.mutations) != (
        result.statistics.total_mutations
    ):
        raise MutationAnalysisError(
            "Mutation count and statistics disagree."
        )

    if not (
        0.0 <= result.identity_percent <= 100.0
        and 0.0 <= result.sequence_similarity_percent <= 100.0
    ):
        raise MutationAnalysisError(
            "Identity/similarity percentages are out of range."
        )

    expected_net = (
        result.variant_length
        - result.reference_length
    )

    if (
        result.statistics.net_length_change
        != expected_net
    ):
        raise MutationAnalysisError(
            "Net length change is inconsistent with sequence lengths."
        )

    ids = [
        mutation.mutation_id
        for mutation in result.mutations
    ]

    if ids != list(
        range(
            1,
            len(ids) + 1,
        )
    ):
        raise MutationAnalysisError(
            "Mutation IDs are not sequential."
        )

    return True


# ============================================================================
# FORMATTING / REPORTING
# ============================================================================

def format_mutation(
    mutation: Mutation,
) -> str:

    extra = ""

    if mutation.effect:
        protein = mutation_hgvs_p(mutation)
        extra = (
            f"  [{mutation.effect}"
            f"{' ' + protein if protein else ''}]"
        )

    return (
        f"#{mutation.mutation_id:<3} "
        f"{mutation.mutation_type:<12} "
        f"position {mutation.position:<6} "
        f"{mutation.notation}"
        f"{extra}"
    )


def build_mutation_report(
    result: MutationAnalysisResult,
    *,
    include_alignment: bool = False,
) -> str:

    lines = [
        major_rule(),
        "MUTATION ANALYSIS REPORT",
        major_rule(),
        style.field_line('Reference', f"{result.reference_name}"),
        style.field_line('Variant', f"{result.variant_name}"),
        style.field_line('Reference length', f"{result.reference_length} bp"),
        style.field_line('Variant length', f"{result.variant_length} bp"),
        style.field_line('Alignment identity', f"{result.identity_percent:.2f}%"),
        style.field_line('Sequence similarity', f"{result.sequence_similarity_percent:.2f}%"),
        style.field_line('Mutations', f"{len(result.mutations)}"),
        "",
        "MUTATION SUMMARY",
        minor_rule(),
        "Alignment identity = exact matches among aligned non-gap bases.",
        "Sequence similarity = normalized edit-distance similarity.",
    ]

    s = result.statistics

    ratio = (
        "inf"
        if s.transition_transversion_ratio
        == math.inf
        else (
            f"{s.transition_transversion_ratio:.3f}"
            if s.transition_transversion_ratio
            is not None
            else "N/A"
        )
    )

    lines.extend(
        [
            style.field_line('SNPs', f"{s.snps}"),
            style.field_line('Substitutions', f"{s.substitutions}"),
            style.field_line('Insertions', f"{s.insertions}"),
            style.field_line('Deletions', f"{s.deletions}"),
            style.field_line('Complex changes', f"{s.complex_changes}"),
            style.field_line('Transitions', f"{s.transitions}"),
            style.field_line('Transversions', f"{s.transversions}"),
            style.field_line('Ti/Tv ratio', f"{ratio}"),
            style.field_line('Inserted bases', f"{s.inserted_bases}"),
            style.field_line('Deleted bases', f"{s.deleted_bases}"),
            style.field_line('Net length', f"{s.net_length_change:+d} bp"),
            style.field_line('Mutation density', f"{s.mutation_density_per_kb:.2f}/kb"),
            "",
            "CODING CONSEQUENCES",
            minor_rule(),
            style.field_line('Synonymous', f"{s.synonymous}"),
            style.field_line('Missense', f"{s.missense}"),
            style.field_line('Nonsense', f"{s.nonsense}"),
            style.field_line('Start-loss', f"{s.start_loss}"),
            style.field_line('Stop-loss', f"{s.stop_loss}"),
            style.field_line('Frameshifts', f"{s.frameshifts}"),
            style.field_line('In-frame indels', f"{s.in_frame_indels}"),
            "",
            "MUTATIONS",
            minor_rule(),
        ]
    )

    if result.mutations:

        for mutation in result.mutations:

            lines.append(
                format_mutation(
                    mutation
                )
            )

    else:

        lines.append(
            "No mutations detected."
        )

    if result.coding_analysis:

        c = result.coding_analysis

        lines.extend(
            [
                "",
                "CODING ANALYSIS",
                minor_rule(),
                style.field_line('Reference protein', f"{c.reference_protein}"),
                style.field_line('Variant protein', f"{c.alternative_protein}"),
                style.field_line('Affected codons', f"{len(c.affected_codons)}"),
                style.field_line('Frameshift', f"{'YES' if c.frameshift else 'NO'}"),
            ]
        )

        if c.frameshift_position is not None:

            lines.append(
                style.field_line('Frameshift position', f"{c.frameshift_position}")
            )

        if c.premature_stop_position is not None:

            lines.append(
                style.field_line('Premature stop', f"codon "
                f"{c.premature_stop_position}")
            )

    if include_alignment:

        lines.extend(
            [
                "",
                "GLOBAL ALIGNMENT",
                minor_rule(),
                result.aligned_reference,
                result.aligned_variant,
            ]
        )

    return "\n".join(lines)


def print_mutation_report(
    result: MutationAnalysisResult,
    *,
    include_alignment: bool = False,
) -> None:

    print(
        build_mutation_report(
            result,
            include_alignment=include_alignment,
        )
    )


# ============================================================================
# MUTATION SIGNATURE
# ============================================================================

def mutation_signature(
    mutations: Sequence[Mutation],
) -> Dict[str, Any]:

    spectrum = build_mutation_spectrum(
        mutations
    )

    return {
        "total": len(mutations),
        "snps": sum(
            1
            for mutation in mutations
            if mutation.mutation_type == "SNP"
        ),
        "insertions": sum(
            1
            for mutation in mutations
            if mutation.mutation_type
            == "Insertion"
        ),
        "deletions": sum(
            1
            for mutation in mutations
            if mutation.mutation_type
            == "Deletion"
        ),
        "transitions": spectrum.transitions,
        "transversions": spectrum.transversions,
        "spectrum_12": spectrum.substitutions_12,
        "spectrum_6": spectrum.substitutions_6,
    }


# ============================================================================
# ADVANCED ANALYSIS
# ============================================================================

def create_advanced_analysis(
    result: MutationAnalysisResult,
) -> AdvancedMutationAnalysis:

    advanced = AdvancedMutationAnalysis()

    advanced.reference_quality = sequence_quality(
        result.aligned_reference.replace(
            "-",
            "",
        )
    )

    advanced.variant_quality = sequence_quality(
        result.aligned_variant.replace(
            "-",
            "",
        )
    )

    advanced.spectra = build_mutation_spectrum(
        result.mutations
    )

    for mutation in result.mutations:

        advanced.contexts.append(
            extract_mutation_context(
                mutation,
                result.aligned_reference.replace(
                    "-",
                    "",
                ),
                result.aligned_variant.replace(
                    "-",
                    "",
                ),
            )
        )

    advanced.hotspots = find_mutation_hotspots(
        result.mutations,
        result.reference_length,
    )

    if (
        advanced.reference_quality.ambiguous_count
        or advanced.variant_quality.ambiguous_count
    ):

        advanced.ambiguity_warnings.append(
            "Ambiguous/IUPAC bases may limit exact mutation interpretation."
        )

    if (
        len(result.aligned_reference)
        != len(result.aligned_variant)
    ):

        advanced.normalization_notes.append(
            "Alignment lengths were inconsistent."
        )

    return advanced


def advanced_mutation_report(
    result: MutationAnalysisResult,
    advanced: Optional[
        AdvancedMutationAnalysis
    ] = None,
) -> str:

    if advanced is None:
        advanced = create_advanced_analysis(
            result
        )

    lines = [
        build_mutation_report(
            result
        ),
        "",
        "ADVANCED SEQUENCE ANALYSIS",
        minor_rule(),
    ]

    if (
        advanced.reference_quality
        and advanced.variant_quality
    ):

        lines.extend(
            [
                style.field_line('Reference GC', f"{advanced.reference_quality.gc_percent:.2f}%"),
                style.field_line('Variant GC', f"{advanced.variant_quality.gc_percent:.2f}%"),
                style.field_line('Reference complexity', f"{advanced.reference_quality.complexity_percent:.2f}%"),
                style.field_line('Variant complexity', f"{advanced.variant_quality.complexity_percent:.2f}%"),
            ]
        )

    lines.extend(
        [
            "",
            "MUTATION SPECTRUM",
            minor_rule(),
        ]
    )

    if advanced.spectra:

        for key, count in (
            advanced.spectra
            .substitutions_12
            .items()
        ):

            lines.append(
                f"{key:<8} {count}"
            )

    lines.extend(
        [
            "",
            "HOTSPOTS",
            minor_rule(),
        ]
    )

    if advanced.hotspots:

        for hotspot in advanced.hotspots:

            lines.append(
                f"{hotspot.start}-{hotspot.end}: "
                f"{hotspot.mutation_count} mutations "
                f"({hotspot.density_per_kb:.2f}/kb)"
            )

    else:

        lines.append(
            "No hotspots detected under current thresholds."
        )

    if advanced.ambiguity_warnings:

        lines.extend(
            [
                "",
                "AMBIGUITY / QUALITY WARNINGS",
                minor_rule(),
            ]
        )

        lines.extend(
            f"- {warning}"
            for warning in advanced.ambiguity_warnings
        )

    return "\n".join(lines)


def mutation_map(
    mutations: Sequence[Mutation],
    sequence_length: int,
    *,
    width: int = 70,
) -> str:

    if sequence_length < 1:
        raise MutationAnalysisError(
            "sequence_length must be positive."
        )

    width = max(
        20,
        width,
    )

    line = [" "] * width

    for mutation in mutations:

        index = min(
            width - 1,
            max(
                0,
                int(
                    (
                        mutation.position - 1
                    )
                    / sequence_length
                    * width
                ),
            ),
        )

        line[index] = "^"

    return "".join(line)


def result_to_json(
    result: MutationAnalysisResult,
    *,
    include_advanced: bool = True,
) -> str:

    payload = result.to_dict()

    if include_advanced:

        payload["advanced"] = asdict(
            create_advanced_analysis(
                result
            )
        )

    return json.dumps(
        payload,
        indent=2,
    )


# ============================================================================
# UI HELPERS
# ============================================================================

def _extract_sequence_record(
    record_key: Any,
    record: Any,
) -> Tuple[str, Any]:

    if isinstance(record, dict):

        header = str(
            record.get("header")
            or record.get("name")
            or record.get("id")
            or record_key
        )

        sequence = (
            record.get("sequence")
            or record.get("dna")
        )

        return header, sequence

    if (
        isinstance(record, (tuple, list))
        and len(record) >= 2
    ):

        return (
            str(record[0]),
            record[1],
        )

    return (
        str(record_key),
        record,
    )


def _iter_sequence_records(
    sequences: Any,
) -> List[Tuple[str, str]]:

    available = []

    if isinstance(
        sequences,
        dict,
    ):

        iterator = sequences.items()

    elif isinstance(
        sequences,
        (list, tuple),
    ):

        iterator = enumerate(
            sequences,
            1,
        )

    else:

        return available

    for key, record in iterator:

        try:

            header, sequence = (
                _extract_sequence_record(
                    key,
                    record,
                )
            )

            if sequence is None:
                continue

            sequence = str(
                sequence
            ).strip()

            if sequence:

                available.append(
                    (
                        str(header),
                        sequence,
                    )
                )

        except Exception:
            continue

    return available


def _prompt_sequence(
    prompt: str,
    sequences=None,
) -> Tuple[str, str]:

    available = _iter_sequence_records(
        sequences
    )

    if available:

        # [0] raises SequencePromptCancelled; both menu branches that ask
        # for sequences catch it and return to the mutation menu.
        style.menu([
            (f"Available sequences for {prompt}",
             [(i, header) for i, (header, _) in enumerate(available, 1)]),
            (None, [("M", "Enter sequence manually")]),
        ], zero=style.CANCEL_LABEL)

        while True:

            choice = style.ask_choice(f"{prompt} selection")

            if choice == "0":
                raise SequencePromptCancelled(prompt)

            if choice.upper() == "M":
                break

            if (
                choice.isdigit()
                and 1 <= int(choice)
                <= len(available)
            ):

                header, sequence = (
                    available[
                        int(choice) - 1
                    ]
                )

                return (
                    header,
                    normalize_iupac_sequence(
                        sequence
                    ),
                )

            style.invalid_selection()

    header = style.ask(f"{prompt} name")

    if not header:
        header = prompt

    sequence = style.ask(f"{prompt} DNA sequence")

    return (
        header,
        normalize_iupac_sequence(
            sequence
        ),
    )


def _store_pipeline_result(
    pipeline,
    key: str,
    value: Any,
) -> None:

    if pipeline is None:
        return

    try:

        if (
            hasattr(
                pipeline,
                "store",
            )
            and callable(
                pipeline.store
            )
        ):

            pipeline.store(
                key,
                value,
            )

        elif (
            hasattr(
                pipeline,
                "results",
            )
            and isinstance(
                pipeline.results,
                dict,
            )
        ):

            pipeline.results[
                key
            ] = value

        elif isinstance(
            pipeline,
            dict,
        ):

            pipeline[key] = value

        else:
            import logging
            logging.getLogger("synapse.mutation").warning(
                "Pipeline object has no result store; %s result not stored.", key
            )

    except Exception:
        # Storing must never destroy a finished analysis, but the failure
        # is logged instead of being hidden (audit Phase V).
        import logging
        logging.getLogger("synapse.mutation").exception(
            "Could not store %s result in the pipeline.", key
        )


# ============================================================================
# MUTATION UI
# ============================================================================

def run_mutation_analysis_ui(
    dna=None,
    pipeline=None,
    header="Sequence",
    sequences=None,
):

    last_result = None

    while True:

        style.header('MUTATION ANALYSIS')
        style.field("Reference", header if dna else "none loaded")
        if dna:
            style.field("Length", f"{len(dna)} bp")

        style.menu([
            ("Comparison", [('1', 'Compare Sequences')]),
            ("Analysis", [('2', 'Mutation Statistics'), ('3', 'Coding Consequence Analysis')]),
            ("Output", [('4', 'Mutation Report')]),
            ("System", [('5', 'Run Mutation Self-Test')]),
        ])

        choice = style.ask_choice()

        # ------------------------------------------------------------------
        # BACK
        # ------------------------------------------------------------------

        if choice == "0":
            return last_result

        # ------------------------------------------------------------------
        # COMPARE
        # ------------------------------------------------------------------

        if choice == "1":

            try:

                if dna:

                    reference_name = (
                        header
                    )

                    reference_sequence = (
                        normalize_iupac_sequence(
                            dna
                        )
                    )

                else:

                    (
                        reference_name,
                        reference_sequence,
                    ) = _prompt_sequence(
                        "Reference",
                        sequences,
                    )

                (
                    variant_name,
                    variant_sequence,
                ) = _prompt_sequence(
                    "Variant",
                    sequences,
                )

                last_result = compare_sequences(
                    reference_sequence,
                    variant_sequence,
                    reference_name=reference_name,
                    variant_name=variant_name,
                )

                # Coding annotations are only meaningful when the user
                # supplies/uses CDS sequences. Do not silently claim that
                # an arbitrary DNA sequence is coding.
                last_result.statistics = (
                    calculate_mutation_statistics(
                        last_result.mutations,
                        reference_length=(
                            last_result.reference_length
                        ),
                    )
                )

                print()

                print_mutation_report(
                    last_result
                )

                _store_pipeline_result(
                    pipeline,
                    "mutation",
                    last_result.to_dict(),
                )

            except SequencePromptCancelled:

                continue

            except MutationAnalysisError as error:

                style.error(f"{error}")

            except Exception as error:

                style.error(f"Mutation analysis failed: {error}")

        # ------------------------------------------------------------------
        # STATISTICS
        # ------------------------------------------------------------------

        elif choice == "2":

            if last_result is None:

                print(
                    "Run a sequence comparison first."
                )

                continue

            s = last_result.statistics

            if (
                s.transition_transversion_ratio
                == math.inf
            ):

                ratio = "inf"

            elif (
                s.transition_transversion_ratio
                is None
            ):

                ratio = "N/A"

            else:

                ratio = (
                    f"{s.transition_transversion_ratio:.3f}"
                )

            style.header('MUTATION STATISTICS')

            style.field('Total mutations', f"{s.total_mutations}")

            style.field('Alignment identity', f"{last_result.identity_percent:.2f}%")

            style.field('Sequence similarity', f"{last_result.sequence_similarity_percent:.2f}%")

            style.field('SNPs', f"{s.snps}")

            style.field('Insertions', f"{s.insertions}")

            style.field('Deletions', f"{s.deletions}")

            style.field('Transitions', f"{s.transitions}")

            style.field('Transversions', f"{s.transversions}")

            style.field('Ti/Tv ratio', f"{ratio}")

            style.field('Net length change', f"{s.net_length_change:+d} bp")

            style.field('Mutation density', f"{s.mutation_density_per_kb:.2f}/kb")

            style.field('Synonymous', f"{s.synonymous}")

            style.field('Missense', f"{s.missense}")

            style.field('Nonsense', f"{s.nonsense}")

            style.field('Frameshifts', f"{s.frameshifts}")

            style.field('In-frame indels', f"{s.in_frame_indels}")

        # ------------------------------------------------------------------
        # CODING CONSEQUENCE
        # ------------------------------------------------------------------

        elif choice == "3":

            try:

                if last_result is None:

                    if dna:

                        ref_name = header

                        ref = normalize_dna(
                            dna,
                            allow_n=False,
                        )

                    else:

                        (
                            ref_name,
                            ref,
                        ) = _prompt_sequence(
                            "Reference",
                            sequences,
                        )

                    (
                        var_name,
                        var,
                    ) = _prompt_sequence(
                        "Variant",
                        sequences,
                    )

                    last_result = compare_sequences(
                        ref,
                        var,
                        reference_name=ref_name,
                        variant_name=var_name,
                    )

                orf_ref_cds = orf_var_cds = None
                selected_orf = getattr(
                    getattr(pipeline, "context", None), "selected_orf", None
                )
                if selected_orf and dna:
                    if style.confirm("Use the selected pipeline ORF as CDS?", default_yes=True):
                        reference_plain = last_result.aligned_reference.replace("-", "")
                        variant_plain = last_result.aligned_variant.replace("-", "")
                        orf_ref_cds, orf_var_cds, orf_text = orf_reference_cds(
                            selected_orf, reference_plain, variant_plain
                        )
                        print(f"CDS from {orf_text}")
                        if orf_var_cds is None:
                            print(
                                "Variant length differs from the reference: enter the "
                                "variant CDS explicitly (ORF coordinates do not carry over)."
                            )
                ref_cds = orf_ref_cds or style.ask('Reference CDS (press Enter to use reference sequence)')

                if not ref_cds:

                    ref_cds = (
                        last_result
                        .aligned_reference
                        .replace(
                            "-",
                            "",
                        )
                    )

                var_cds = orf_var_cds or style.ask('Variant CDS (press Enter to use variant sequence)')

                if not var_cds:

                    var_cds = (
                        last_result
                        .aligned_variant
                        .replace(
                            "-",
                            "",
                        )
                    )

                coding = (
                    annotate_mutations_with_coding_consequences(
                        last_result.mutations,
                        ref_cds,
                        var_cds,
                    )
                )

                frame_warnings = assess_coding_frame(ref_cds)
                coding.notes.append(CODING_FRAME_ASSUMPTION)
                coding.notes.extend(
                    f"Reference CDS: {warning}" for warning in frame_warnings
                )
                if frame_warnings:
                    print()
                    print("WARNING - coding frame assumption:")
                    for warning in frame_warnings:
                        print(f"  - {warning}")
                last_result.coding_analysis = (
                    coding
                )

                last_result.statistics = (
                    calculate_mutation_statistics(
                        last_result.mutations,
                        reference_length=(
                            last_result.reference_length
                        ),
                    )
                )

                style.header('CODING CONSEQUENCE ANALYSIS')

                style.field('Reference protein', f"{coding.reference_protein}")

                style.field('Variant protein', f"{coding.alternative_protein}")

                style.field('Synonymous', f"{coding.synonymous_count}")

                style.field('Missense', f"{coding.missense_count}")

                style.field('Nonsense', f"{coding.nonsense_count}")

                style.field('Start-loss', f"{coding.start_loss_count}")

                style.field('Stop-loss', f"{coding.stop_loss_count}")

                style.field('Frameshift', f"{'YES' if coding.frameshift else 'NO'}")

                if (
                    coding.frameshift_position
                    is not None
                ):

                    style.field('Frameshift position', f"{coding.frameshift_position}")

                if (
                    coding.premature_stop_position
                    is not None
                ):

                    style.field('Premature stop', f"codon "
                        f"{coding.premature_stop_position}")

            except SequencePromptCancelled:

                continue

            except MutationAnalysisError as error:

                style.error(f"{error}")

            except Exception as error:

                style.error(f"Coding analysis failed: {error}")

        # ------------------------------------------------------------------
        # REPORT
        # ------------------------------------------------------------------

        elif choice == "4":

            if last_result is None:

                print(
                    "Run a sequence comparison first."
                )

            else:

                print()

                print_mutation_report(
                    last_result,
                    include_alignment=True,
                )

        # ------------------------------------------------------------------
        # SELF TEST
        # ------------------------------------------------------------------

        elif choice == "5":

            tests = (
                run_mutation_self_test()
            )

            style.header('MUTATION SELF-TEST')

            for name, passed in (
                tests.items()
            ):

                print(
                    f"[{'PASS' if passed else 'FAIL'}] "
                    f"{name}"
                )

            print(
                minor_rule()
            )

            print(
                f"{sum(tests.values())}/"
                f"{len(tests)} tests passed."
            )

        else:

            style.invalid_selection()


# ============================================================================
# TOP-LEVEL MUTATION + VARIANT UI
# ============================================================================

def run_mutation_variant_ui(
    dna=None,
    header="Sequence",
    sequences=None,
    pipeline=None,
):

    while True:

        style.header('MUTATION & VARIANT ANALYSIS')
        style.field("Sequence", header if dna else "none loaded")
        if dna:
            style.field("Length", f"{len(dna)} bp")

        style.menu([
            ("Analysis", [('1', 'Mutation Analysis'), ('2', 'Variant Analysis')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":

            return run_mutation_analysis_ui(
                dna=dna,
                pipeline=pipeline,
                header=header,
                sequences=sequences,
            )

        if choice == "2":

            try:

                from analyses.variant_analysis import (
                    run_variant_analysis_ui
                )

                return run_variant_analysis_ui(
                    pipeline=pipeline,
                    name=(
                        f"Variant Analysis - "
                        f"{header}"
                    ),
                )

            except ImportError as error:

                style.header('VARIANT ANALYSIS IMPORT ERROR')

                print()

                print(
                    f"{type(error).__name__}: "
                    f"{error}"
                )

                print()

                print(
                    "Check that "
                    "analyses/variant_analysis.py "
                    "exists and can be imported."
                )

                style.pause()

            except Exception as error:

                style.header('VARIANT ANALYSIS ERROR')

                print()

                print(
                    f"{type(error).__name__}: "
                    f"{error}"
                )

                style.pause()

        else:

            style.invalid_selection()


# ============================================================================
# SELF TEST
# ============================================================================

def run_mutation_self_test() -> Dict[str, bool]:

    tests: Dict[str, bool] = {}

    # ------------------------------------------------------------------
    # SNP
    # ------------------------------------------------------------------

    r = detect_mutations(
        "ATGCCG",
        "ATGTCG",
    )

    tests["snp_detection"] = (
        len(r.mutations) == 1
        and r.mutations[0].mutation_type
        == "SNP"
        and r.mutations[0].position
        == 4
        and r.mutations[0].transition
        is True
    )

    # ------------------------------------------------------------------
    # Multiple substitutions
    # ------------------------------------------------------------------

    r = detect_mutations(
        "ATGCCG",
        "GTGCTG",
    )

    tests["multiple_substitutions"] = (
        len(r.mutations) == 2
    )

    # ------------------------------------------------------------------
    # Insertion
    # ------------------------------------------------------------------

    r = detect_mutations(
        "ATGCCG",
        "ATGGCCG",
    )

    tests["insertion_detection"] = (
        r.statistics.insertions == 1
        and r.statistics.inserted_bases == 1
        and r.statistics.net_length_change == 1
    )

    # ------------------------------------------------------------------
    # Deletion
    # ------------------------------------------------------------------

    r = detect_mutations(
        "ATGGCCG",
        "ATGCCG",
    )

    tests["deletion_detection"] = (
        r.statistics.deletions == 1
        and r.statistics.deleted_bases == 1
        and r.statistics.net_length_change == -1
    )

    # ------------------------------------------------------------------
    # Synonymous
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "GCT",
        "GCC",
    )

    tests["synonymous"] = (
        coding.synonymous_count == 1
        and coding.missense_count == 0
    )

    # ------------------------------------------------------------------
    # Missense
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "GCT",
        "GTT",
    )

    tests["missense"] = (
        coding.missense_count == 1
        and coding.synonymous_count == 0
    )

    # ------------------------------------------------------------------
    # Nonsense
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "CAG",
        "TAG",
    )

    tests["nonsense"] = (
        coding.nonsense_count == 1
    )

    # ------------------------------------------------------------------
    # Start loss
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "ATG",
        "ACG",
    )

    tests["start_loss"] = (
        coding.start_loss_count == 1
    )

    # ------------------------------------------------------------------
    # Stop loss
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "TAA",
        "CAA",
    )

    tests["stop_loss"] = (
        coding.stop_loss_count == 1
    )

    # ------------------------------------------------------------------
    # Frameshift
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "ATGCCCGAA",
        "ATGCCGAA",
    )

    tests["frameshift"] = (
        coding.frameshift is True
    )

    # ------------------------------------------------------------------
    # In-frame indel
    # ------------------------------------------------------------------

    coding = analyze_coding_consequences(
        "ATGAAACCC",
        "ATGCCC",
    )

    tests["in_frame_indel"] = (
        coding.frameshift is False
    )

    # ------------------------------------------------------------------
    # IUPAC
    # ------------------------------------------------------------------

    tests["iupac_matching"] = (
        iupac_matches(
            "R",
            "A",
        )
        and not iupac_matches(
            "R",
            "C",
        )
    )

    # ------------------------------------------------------------------
    # IUPAC reverse complement
    # ------------------------------------------------------------------

    tests["iupac_reverse_complement"] = (
        reverse_complement_iupac(
            "ARYN"
        )
        == "NRYT"
    )

    # ------------------------------------------------------------------
    # Translation
    # ------------------------------------------------------------------

    tests["translation"] = (
        translate_dna(
            "ATGGCC"
        )
        == "MA"
    )

    # ------------------------------------------------------------------
    # Reverse complement
    # ------------------------------------------------------------------

    tests["reverse_complement"] = (
        reverse_complement(
            "ATGC"
        )
        == "GCAT"
    )

    # ------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------

    r = detect_mutations(
        "AAGC",
        "GAGG",
    )

    tests["statistics"] = (
        r.statistics.snps == 2
        and r.statistics.transitions == 1
        and r.statistics.transversions == 1
    )

    # ------------------------------------------------------------------
    # Multiple variants
    # ------------------------------------------------------------------

    multi = compare_multiple_variants(
        "ATGCCG",
        [
            "ATGTCG",
            "ATGCCG",
            "ATGGCCG",
        ],
    )

    tests["multiple_variants"] = (
        len(multi) == 3
        and multi[
            "Variant 2"
        ].mutation_count == 0
        and multi[
            "Variant 3"
        ].statistics.insertions == 1
    )

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    tests["serialization"] = (
        '"mutations"'
        in detect_mutations(
            "ATGC",
            "ATGT",
        ).to_json()
    )

    # ------------------------------------------------------------------
    # Advanced analysis
    # ------------------------------------------------------------------

    advanced_result = detect_mutations(
        "ATGCCG",
        "ATGTCG",
    )

    advanced = create_advanced_analysis(
        advanced_result
    )

    tests["advanced_analysis"] = (
        advanced.reference_quality
        is not None
        and advanced.variant_quality
        is not None
        and advanced.spectra
        is not None
        and len(advanced.contexts) == 1
    )

    # ------------------------------------------------------------------
    # JSON advanced
    # ------------------------------------------------------------------

    tests["json_export"] = (
        '"advanced"'
        in result_to_json(
            advanced_result
        )
    )

    # ------------------------------------------------------------------
    # Restriction impact
    # ------------------------------------------------------------------

    restriction = restriction_site_impact(
        "GAATTC",
        "GACTTC",
        {
            "EcoRI": "GAATTC"
        },
    )

    tests["restriction_impact"] = (
        len(restriction) == 1
        and len(
            restriction[0].lost_sites
        ) == 1
    )

    # ------------------------------------------------------------------
    # Primer overlap
    # ------------------------------------------------------------------

    primers = primer_mutation_overlap(
        detect_mutations(
            "ATGCCG",
            "ATGTCG",
        ).mutations,
        [
            ("P1", 1, 5)
        ],
    )

    tests["primer_overlap"] = (
        primers.primer_overlap is True
    )

    # ------------------------------------------------------------------
    # CRISPR overlap
    # ------------------------------------------------------------------

    targets = crispr_mutation_overlap(
        detect_mutations(
            "ATGCCG",
            "ATGTCG",
        ).mutations,
        [
            ("G1", 3, 6)
        ],
    )

    tests["crispr_overlap"] = (
        targets.crispr_target_overlap is True
    )

    # ------------------------------------------------------------------
    # Quality
    # ------------------------------------------------------------------

    quality = sequence_quality(
        "ATGCCG"
    )

    tests["sequence_quality"] = (
        quality.valid
        and quality.length == 6
        and quality.gc_percent
        > 0
    )

    # ------------------------------------------------------------------
    # Invalid sequence
    # ------------------------------------------------------------------

    try:

        normalize_dna(
            "ATGX"
        )

        tests[
            "invalid_dna_rejected"
        ] = False

    except MutationAnalysisError:

        tests[
            "invalid_dna_rejected"
        ] = True

    # ------------------------------------------------------------------
    # Empty sequence
    # ------------------------------------------------------------------

    try:

        normalize_dna(
            ""
        )

        tests[
            "empty_dna_rejected"
        ] = False

    except MutationAnalysisError:

        tests[
            "empty_dna_rejected"
        ] = True

    # ------------------------------------------------------------------
    # Result consistency
    # ------------------------------------------------------------------

    result = detect_mutations(
        "ATGCCG",
        "ATGTCG",
    )

    try:

        tests[
            "result_validation"
        ] = validate_mutation_result(
            result
        )

    except MutationAnalysisError:

        tests[
            "result_validation"
        ] = False

    # ------------------------------------------------------------------
    # Mutation map
    # ------------------------------------------------------------------

    mutation_map_text = mutation_map(
        result.mutations,
        result.reference_length,
    )

    tests["mutation_map"] = (
        "^" in mutation_map_text
    )


    # ------------------------------------------------------------------
    # Professional identity vs. sequence similarity
    # ------------------------------------------------------------------

    deletion_result = detect_mutations(
        "ATGAAACCC",
        "ATGCCC",
    )

    tests["identity_similarity_separation"] = (
        deletion_result.identity_percent == 100.0
        and deletion_result.sequence_similarity_percent
        < deletion_result.identity_percent
        and abs(
            deletion_result.sequence_similarity_percent
            - (100.0 * (1.0 - 3.0 / 9.0))
        ) < 1e-9
    )

    # ------------------------------------------------------------------
    # Ambiguous substitution must not be a complex change
    # ------------------------------------------------------------------

    ambiguous_result = detect_mutations(
        "ATGCCG",
        "ATGRCG",
    )

    tests["ambiguous_substitution_statistics"] = (
        len(ambiguous_result.mutations) == 1
        and ambiguous_result.mutations[0].mutation_type
        == "Substitution"
        and ambiguous_result.statistics.substitutions == 1
        and ambiguous_result.statistics.complex_changes == 0
        and ambiguous_result.statistics.snps == 0
        and ambiguous_result.statistics.unknown_effects == 1
    )

    # ------------------------------------------------------------------
    # TEST_8 regression suite
    # ------------------------------------------------------------------
    # This suite targets the previously observed failures:
    # - nonsense test must contain exactly one SNP
    # - in-frame indels must not create downstream missense cascades
    # - IUPAC R must be accepted
    # - terminal substitutions must remain valid sequence-level events
    test8_ref = (
        "ATGAAACGTGCTGACTTTGGTGACGCTTACCGTAAAGGTCCTTAAACTGCTGACCGTTGACGAA"
        "CTGAAAGCTGGTGCTTACGGTGACCTGATCGCTGAA"
    )

    # Build every test variant from the same reference so that each intended
    # regression contains only the mutation it claims to contain.
    test8_synonymous = test8_ref[:5] + "G" + test8_ref[6:]
    test8_missense = test8_ref[:10] + "T" + test8_ref[11:]

    # Codon 23 is GAA (Glu). GAA -> TAA is a single SNP at position 67.
    test8_nonsense = test8_ref[:66] + "T" + test8_ref[67:]

    # Three independent SNPs: one transition, one transition, one transversion.
    test8_multi = (
        test8_ref[:5] + "G" + test8_ref[6:11]
        + "C" + test8_ref[12:66]
        + "T" + test8_ref[67:]
    )

    # Two transitions: A>G and T>C.
    test8_transitions = (
        test8_ref[:5] + "G" + test8_ref[6:11]
        + "C" + test8_ref[12:]
    )

    # Two transversions: C>A and T>G.
    test8_transversions = (
        test8_ref[:10] + "A" + "G" + test8_ref[12:]
    )

    test8_inframe_insertion = (
        test8_ref[:30] + "GCT" + test8_ref[30:]
    )
    test8_inframe_deletion = (
        test8_ref[:30] + test8_ref[33:]
    )
    test8_frameshift_insertion = (
        test8_ref[:30] + "G" + test8_ref[30:]
    )
    test8_frameshift_deletion = (
        test8_ref[:30] + test8_ref[31:]
    )

    test8_ambiguous = (
        test8_ref[:20] + "R" + test8_ref[21:]
    )

    test8_cases = {
        "TEST_8_SNP_SYNONYMOUS": test8_synonymous,
        "TEST_8_SNP_MISSENSE": test8_missense,
        "TEST_8_SNP_NONSENSE": test8_nonsense,
        "TEST_8_MULTI_SNP": test8_multi,
        "TEST_8_TRANSITIONS": test8_transitions,
        "TEST_8_TRANSVERSIONS": test8_transversions,
        "TEST_8_INFRAME_INSERTION": test8_inframe_insertion,
        "TEST_8_INFRAME_DELETION": test8_inframe_deletion,
        "TEST_8_FRAMESHIFT_INSERTION": test8_frameshift_insertion,
        "TEST_8_FRAMESHIFT_DELETION": test8_frameshift_deletion,
        "TEST_8_AMBIGUOUS": test8_ambiguous,
    }

    try:
        r8 = {
            name: detect_mutations(test8_ref, sequence)
            for name, sequence in test8_cases.items()
        }

        tests["TEST_8_SNP_SYNONYMOUS"] = (
            r8["TEST_8_SNP_SYNONYMOUS"].statistics.snps == 1
            and r8["TEST_8_SNP_SYNONYMOUS"].statistics.synonymous == 1
        )

        tests["TEST_8_SNP_MISSENSE"] = (
            r8["TEST_8_SNP_MISSENSE"].statistics.snps == 1
            and r8["TEST_8_SNP_MISSENSE"].statistics.missense == 1
        )

        tests["TEST_8_SNP_NONSENSE"] = (
            len(r8["TEST_8_SNP_NONSENSE"].mutations) == 1
            and r8["TEST_8_SNP_NONSENSE"].statistics.snps == 1
            and r8["TEST_8_SNP_NONSENSE"].statistics.nonsense == 1
            and r8["TEST_8_SNP_NONSENSE"].statistics.insertions == 0
            and r8["TEST_8_SNP_NONSENSE"].statistics.deletions == 0
            and r8["TEST_8_SNP_NONSENSE"].statistics.frameshifts == 0
        )

        tests["TEST_8_MULTI_SNP"] = (
            r8["TEST_8_MULTI_SNP"].statistics.snps == 3
            and (
                r8["TEST_8_MULTI_SNP"].statistics.transitions
                + r8["TEST_8_MULTI_SNP"].statistics.transversions
            ) == 3
        )

        tests["TEST_8_TRANSITIONS"] = (
            r8["TEST_8_TRANSITIONS"].statistics.snps == 2
            and r8["TEST_8_TRANSITIONS"].statistics.transitions == 2
            and r8["TEST_8_TRANSITIONS"].statistics.transversions == 0
        )

        tests["TEST_8_TRANSVERSIONS"] = (
            r8["TEST_8_TRANSVERSIONS"].statistics.snps == 2
            and r8["TEST_8_TRANSVERSIONS"].statistics.transversions == 2
            and r8["TEST_8_TRANSVERSIONS"].statistics.transitions == 0
        )

        tests["TEST_8_INFRAME_INSERTION"] = (
            r8["TEST_8_INFRAME_INSERTION"].statistics.insertions == 1
            and r8["TEST_8_INFRAME_INSERTION"].statistics.in_frame_indels == 1
            and r8["TEST_8_INFRAME_INSERTION"].statistics.missense == 0
            and r8["TEST_8_INFRAME_INSERTION"].statistics.frameshifts == 0
        )

        tests["TEST_8_INFRAME_DELETION"] = (
            r8["TEST_8_INFRAME_DELETION"].statistics.deletions == 1
            and r8["TEST_8_INFRAME_DELETION"].statistics.in_frame_indels == 1
            and r8["TEST_8_INFRAME_DELETION"].statistics.missense == 0
            and r8["TEST_8_INFRAME_DELETION"].statistics.frameshifts == 0
        )

        tests["TEST_8_FRAMESHIFT_INSERTION"] = (
            r8["TEST_8_FRAMESHIFT_INSERTION"].statistics.insertions == 1
            and r8["TEST_8_FRAMESHIFT_INSERTION"].statistics.frameshifts == 1
        )

        tests["TEST_8_FRAMESHIFT_DELETION"] = (
            r8["TEST_8_FRAMESHIFT_DELETION"].statistics.deletions == 1
            and r8["TEST_8_FRAMESHIFT_DELETION"].statistics.frameshifts == 1
        )

        tests["TEST_8_AMBIGUOUS"] = (
            len(r8["TEST_8_AMBIGUOUS"].mutations) == 1
            and r8["TEST_8_AMBIGUOUS"].mutations[0].effect == "Unknown"
        )

    except Exception:
        for key in test8_cases:
            tests[key] = False

    # ------------------------------------------------------------------
    # Reverse-complement orientation
    # ------------------------------------------------------------------

    reverse_result = detect_mutations(
        test8_ref,
        reverse_complement_iupac(test8_ref),
    )

    tests["TEST_8_REVERSE"] = (
        len(reverse_result.mutations) == 0
        and reverse_result.reference_length
        == reverse_result.variant_length
        and any(
            "reverse complement"
            in note.lower()
            for note in reverse_result.notes
        )
    )

    # ------------------------------------------------------------------
    # Direct coding consequence: reverse-complement orientation
    # ------------------------------------------------------------------

    reverse_coding = analyze_coding_consequences(
        test8_ref,
        reverse_complement_iupac(test8_ref),
    )

    tests["TEST_8_REVERSE_CODING"] = (
        reverse_coding.missense_count == 0
        and reverse_coding.nonsense_count == 0
        and reverse_coding.synonymous_count == 0
        and reverse_coding.start_loss_count == 0
        and reverse_coding.stop_loss_count == 0
        and reverse_coding.frameshift is False
        and reverse_coding.reference_protein
        == reverse_coding.alternative_protein
        and any(
            "orientation-only"
            in note.lower()
            for note in reverse_coding.notes
        )
    )

    # ------------------------------------------------------------------
    # Direct coding consequence: ambiguous IUPAC base
    # ------------------------------------------------------------------

    ambiguous_coding = analyze_coding_consequences(
        test8_ref,
        test8_ambiguous,
    )

    tests["TEST_8_AMBIGUOUS_CODING"] = (
        ambiguous_coding.reference_protein == ""
        and ambiguous_coding.alternative_protein == ""
        and any(
            "ambiguous IUPAC bases"
            in note
            for note in ambiguous_coding.notes
        )
    )

    return tests


def assert_mutation_self_test() -> None:

    results = run_mutation_self_test()

    failed = [
        name
        for name, passed
        in results.items()
        if not passed
    ]

    if failed:

        raise AssertionError(
            "Mutation self-test failed: "
            + ", ".join(failed)
        )


# ============================================================================
# PUBLIC API
# ============================================================================

analyze_mutations = detect_mutations
mutation_report = build_mutation_report


__all__ = [
    # Core
    "AlignmentBlock",
    "Mutation",
    "MutationStatistics",
    "CodingAnalysis",
    "MutationAnalysisResult",
    "MutationAnalysisError",

    # Sequence
    "normalize_dna",
    "validate_dna_sequence",
    "normalize_iupac_sequence",
    "gc_content",
    "reverse_complement",
    "reverse_complement_iupac",

    # Translation
    "translate_dna",
    "amino_acid_name",

    # Alignment
    "global_align",
    "alignment_identity",
    "sequence_similarity",
    "build_alignment_blocks",

    # Mutation classification
    "is_transition",
    "classify_substitution",
    "detect_mutations_from_alignment",
    "detect_mutations",
    "analyze_mutations",

    # Coding
    "classify_coding_substitution",
    "analyze_coding_consequences",
    "annotate_mutations_with_coding_consequences",

    # Statistics
    "calculate_mutation_statistics",

    # Comparison
    "compare_sequences",
    "compare_multiple_variants",
    "shared_mutations",
    "unique_mutations",

    # Regions
    "annotate_regions",
    "annotate_sequence_regions",

    # Quality
    "sequence_quality",
    "sequence_complexity",

    # Context / spectrum
    "MutationContext",
    "MutationSpectrum",
    "MutationHotspot",
    "extract_mutation_context",
    "build_mutation_spectrum",
    "find_mutation_hotspots",

    # IUPAC
    "IUPAC_CODES",
    "iupac_matches",
    "iupac_code",

    # Indels / HGVS
    "normalize_indel_left",
    "mutation_hgvs_c",
    "mutation_hgvs_p",

    # Protein
    "compare_proteins",

    # Cross-module impact
    "RestrictionSiteImpact",
    "MutationImpact",
    "restriction_site_impact",
    "primer_mutation_overlap",
    "crispr_mutation_overlap",

    # Advanced
    "AdvancedMutationAnalysis",
    "create_advanced_analysis",
    "advanced_mutation_report",

    # Validation / output
    "validate_mutation_result",
    "format_mutation",
    "build_mutation_report",
    "mutation_report",
    "print_mutation_report",
    "mutation_signature",
    "mutation_map",
    "result_to_json",

    # UI
    "run_mutation_analysis_ui",
    "run_mutation_variant_ui",

    # Tests
    "run_mutation_self_test",
    "assert_mutation_self_test",
]


# ============================================================================
# STANDALONE EXECUTION
# ============================================================================

if __name__ == "__main__":

    results = run_mutation_self_test()

    passed = sum(
        results.values()
    )

    style.header('Mutation analysis self-test')

    for name, success in (
        results.items()
    ):

        print(
            f"[{'PASS' if success else 'FAIL'}] "
            f"{name}"
        )

    print(
        minor_rule()
    )

    print(
        f"{passed}/{len(results)} tests passed."
    )

    if passed != len(results):

        raise SystemExit(1)
