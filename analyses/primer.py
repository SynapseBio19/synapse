"""
Synapse — Professional Primer Design Engine v2
================================================

A self-contained primer-design module for PCR-style primer design.

v2 adds a documented nearest-neighbor thermodynamic model with configurable:
- monovalent salt (Na+) concentration
- Mg2+ concentration
- dNTP concentration
- total primer concentration
- primer self-complementarity symmetry correction

The default Tm model remains configurable. The recommended default is the
nearest-neighbor model. Wallace and empirical salt-adjusted screening remain
available for comparison.

Scientific scope and limitations:
- Nearest-neighbor Tm uses SantaLucia-style DNA/DNA thermodynamic parameters.
- Mixed-ion handling uses an effective monovalent ionic-strength approximation;
  it is an approximation, not a full unified-ion Owczarzy/Primer3 calculation.
- Mg2+ is adjusted for dNTP binding through a simple free-Mg approximation.
- Hairpin, self-dimer and heterodimer detection remain sequence-level heuristics;
  they do not calculate full secondary-structure partition functions or ΔG.
- Specificity screening is limited to reference sequences supplied by the caller.
  No genome-wide specificity claim is made without an actual database search.
- Experimental primer performance must be validated under the actual reaction
  conditions.

Coordinates:
- All reported sequence coordinates are 0-based, half-open [start, end).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import lru_cache
from operator import eq as _equal
from math import log, log10, sqrt
from typing import Optional, Sequence, cast
from core.seqtools import IUPAC_COMPLEMENT as _CORE_IUPAC_COMPLEMENT
from ui.style import major_rule, minor_rule
from ui import style

DNA_ALPHABET = set("ACGT")
IUPAC_ALPHABET = set("ACGTRYSWKMBDHVN")

# IUPAC complement table shared via core.seqtools.
IUPAC_COMPLEMENT = dict(_CORE_IUPAC_COMPLEMENT)

# Names of the nearest-neighbor Tm model. ProfessionalPrimerConfig accepts
# "nearest-neighbor", so every dispatch accepts it too (it raised
# "Unknown Tm method" in candidate generation before 2026-11-12).
NEAREST_NEIGHBOR_TM_METHODS = frozenset(
    {"nearest_neighbor", "nearest-neighbor", "nn", "santalucia", "santa_lucia"}
)

# SantaLucia 1998 DNA/DNA nearest-neighbor parameters.
# Values: (delta_H kcal/mol, delta_S cal/(mol*K)).
# Keys are written as 5'-XY-3' / 3'-X'Y'-5' dinucleotide stacks.
NN_PARAMETERS = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}

# SantaLucia initiation terms.
NN_INITIATION = (0.2, -5.7)
NN_AT_TERMINAL = (2.2, 6.9)
NN_GC_TERMINAL = (0.0, -0.1)
NN_SYMMETRY_DS = -1.4
GAS_CONSTANT = 1.98720425864083  # cal/(K*mol)


# ============================================================================
# DATA MODELS
# ============================================================================

@dataclass(frozen=True)
class PrimerDesignConfig:
    """Central configuration for primer-design thresholds and thermodynamics."""

    min_length: int = 18
    max_length: int = 25
    optimal_length: int = 20

    min_gc: float = 40.0
    max_gc: float = 60.0
    optimal_gc: float = 50.0

    min_tm: float = 57.0
    max_tm: float = 63.0
    optimal_tm: float = 60.0
    max_tm_delta: float = 2.0

    min_amplicon_length: int = 80
    max_amplicon_length: int = 1000
    optimal_amplicon_length: int = 250

    max_homopolymer: int = 4
    max_repeat_unit: int = 3
    max_repeat_copies: int = 3

    gc_clamp_min: int = 1
    gc_clamp_max: int = 2

    max_hairpin_stem: int = 4
    max_self_dimer_stem: int = 4
    max_heterodimer_stem: int = 4
    max_three_prime_complementarity: int = 3

    # Thermodynamic reaction conditions.
    salt_mM: float = 50.0
    mg_mM: float = 1.5
    dntp_mM: float = 0.2
    primer_concentration_nM: float = 250.0
    include_symmetry_correction: bool = True

    max_candidates_per_orientation: int = 3000
    top_pairs: int = 10


@dataclass
class ThermodynamicResult:
    tm: float
    delta_h_kcal: float
    delta_s_cal: float
    concentration_nM: float
    salt_mM: float
    mg_mM: float
    dntp_mM: float
    effective_monovalent_mM: float
    model: str
    symmetry_correction_applied: bool
    explanation: str


@dataclass
class StructureResult:
    max_stem: int
    three_prime_stem: int
    best_position: Optional[int]
    severity: str
    explanation: str


@dataclass
class SpecificityHit:
    sequence_name: str
    start: int
    end: int
    mismatches: int
    three_prime_mismatches: int
    exact: bool
    orientation: str


@dataclass
class PrimerCandidate:
    sequence: str
    orientation: str
    start: int
    end: int
    length: int
    gc_percent: float
    tm: float
    tm_model: str
    thermodynamics: ThermodynamicResult
    gc_clamp: int
    max_homopolymer: int
    max_repeat_copies: int
    hairpin: StructureResult
    self_dimer: StructureResult
    score: float
    category: str
    rejection_reasons: list[str] = field(default_factory=list)


@dataclass
class PrimerPair:
    forward: PrimerCandidate
    reverse: PrimerCandidate
    amplicon_start: int
    amplicon_end: int
    amplicon_length: int
    amplicon_sequence: str
    amplicon_gc_percent: float
    tm_delta: float
    heterodimer: StructureResult
    score: float
    category: str
    specificity_hits: list[SpecificityHit] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass
class PrimerDesignReport:
    sequence_length: int
    target_start: int
    target_end: int
    target_length: int
    tm_model: str
    thermodynamic_conditions: dict
    forward_candidates_evaluated: int
    reverse_candidates_evaluated: int
    forward_candidates_accepted: int
    reverse_candidates_accepted: int
    pairs_evaluated: int
    pairs_accepted: int
    ranked_pairs: list[PrimerPair]
    failure_reasons: list[str]
    warnings: list[str]


# ============================================================================
# VALIDATION
# ============================================================================

def normalize_dna_sequence(sequence: str, allow_iupac: bool = False) -> str:
    """Normalize DNA by removing whitespace, uppercasing and converting U to T."""
    if not isinstance(sequence, str):
        raise TypeError("DNA sequence must be a string.")
    return _normalize_dna_cached(sequence, allow_iupac)


# Primer design normalizes the same short strings about two million times
# for a 3 kb template (P2-6). The result depends only on the string and the
# alphabet flag, so it is cached; an invalid input is never cached and
# raises the same error every time.
@lru_cache(maxsize=65536)
def _normalize_dna_cached(sequence: str, allow_iupac: bool) -> str:
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    if not normalized:
        raise ValueError("DNA sequence is empty.")
    alphabet = IUPAC_ALPHABET if allow_iupac else DNA_ALPHABET
    invalid = sorted(set(normalized) - alphabet)
    if invalid:
        raise ValueError("Invalid DNA character(s): " + ", ".join(invalid))
    return normalized


def validate_target_region(sequence: str, target_start: int, target_end: int) -> None:
    """Validate a 0-based, half-open target interval."""
    if target_start < 0 or target_end < 0:
        raise ValueError("Target coordinates must be non-negative.")
    if target_start >= target_end:
        raise ValueError("target_start must be smaller than target_end.")
    if target_end > len(sequence):
        raise ValueError("Target region extends beyond the sequence.")


# ============================================================================
# BASIC SEQUENCE ANALYSIS
# ============================================================================

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement."""
    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    # Every IUPAC symbol has a complement in the table, so translating and
    # reversing gives the same string as a per-base lookup.
    return sequence.translate(_IUPAC_COMPLEMENT_TABLE)[::-1]


_IUPAC_COMPLEMENT_TABLE = str.maketrans(IUPAC_COMPLEMENT)


def gc_percent(sequence: str) -> float:
    """Return GC content as a percentage."""
    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    gc = sequence.count("G") + sequence.count("C")
    return 100.0 * gc / len(sequence)


def count_gc_clamp(sequence: str) -> int:
    """Count consecutive G/C bases at the 3' terminus."""
    sequence = normalize_dna_sequence(sequence)
    count = 0
    for base in reversed(sequence):
        if base in "GC":
            count += 1
        else:
            break
    return count


def longest_homopolymer(sequence: str) -> int:
    """Return the longest consecutive identical-base run."""
    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    best = 1
    current = 1
    for previous, base in zip(sequence, sequence[1:]):
        current = current + 1 if previous == base else 1
        best = max(best, current)
    return best


def longest_tandem_repeat(sequence: str, max_unit: int = 3) -> int:
    """Return the largest number of consecutive copies of a short repeat unit."""
    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    if len(sequence) < 2:
        return 1
    best = 1
    max_unit = min(max_unit, len(sequence) // 2)
    for unit_length in range(1, max_unit + 1):
        # Copy k + 1 of a unit starting at s repeats copy k exactly when
        # sequence[p] == sequence[p + unit_length] for the unit_length
        # positions of copy k. A run of m such matches from s therefore
        # holds 1 + m // unit_length copies, and the longest run gives the
        # most copies for this unit length (P2-6: the same result as
        # comparing the copies one by one).
        matches = bytes(map(_equal, sequence, sequence[unit_length:]))
        longest_run = max(map(len, matches.split(b"\x00")))
        best = max(best, 1 + longest_run // unit_length)
    return best


# ============================================================================
# THERMODYNAMICS
# ============================================================================

# Superseded definition of calculate_tm_nearest_neighbor() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


def calculate_tm_wallace(sequence: str) -> float:
    """Wallace-rule screening estimate: Tm = 2(A+T) + 4(G+C)."""
    sequence = normalize_dna_sequence(sequence)
    at = sequence.count("A") + sequence.count("T")
    gc = sequence.count("G") + sequence.count("C")
    return float(2 * at + 4 * gc)


def calculate_tm_salt_adjusted(sequence: str, salt_mM: float = 50.0) -> float:
    """Empirical monovalent-salt-adjusted Tm estimate."""
    sequence = normalize_dna_sequence(sequence)
    if salt_mM <= 0:
        raise ValueError("salt_mM must be greater than zero.")
    n = len(sequence)
    gc = sequence.count("G") + sequence.count("C")
    salt_molar = salt_mM / 1000.0
    return 64.9 + 41.0 * (gc - 16.4) / n + 16.6 * log10(salt_molar)


def calculate_tm(
    sequence: str,
    method: str = "nearest_neighbor",
    salt_mM: float = 50.0,
    mg_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_concentration_nM: float = 250.0,
    include_symmetry_correction: bool = True,
) -> float:
    """Calculate Tm using Wallace, salt-adjusted, or nearest-neighbor model."""
    method = method.lower().strip()
    if method == "wallace":
        return calculate_tm_wallace(sequence)
    if method in {"salt", "salt_adjusted"}:
        return calculate_tm_salt_adjusted(sequence, salt_mM)
    if method in NEAREST_NEIGHBOR_TM_METHODS:
        return float(calculate_tm_nearest_neighbor(
            sequence,
            salt_mM=salt_mM,
            mg_mM=mg_mM,
            dntp_mM=dntp_mM,
            primer_concentration_nM=primer_concentration_nM,
            include_symmetry_correction=include_symmetry_correction,
        ))
    raise ValueError(
        "Unknown Tm method. Use 'nearest_neighbor', 'wallace', or 'salt_adjusted'."
    )


def get_thermodynamic_conditions(config: PrimerDesignConfig) -> dict:
    """Return the configured reaction conditions in a serializable dictionary."""
    return {
        "salt_mM": config.salt_mM,
        "mg_mM": config.mg_mM,
        "dntp_mM": config.dntp_mM,
        "primer_concentration_nM": config.primer_concentration_nM,
        "include_symmetry_correction": config.include_symmetry_correction,
    }


# ============================================================================
# SECONDARY STRUCTURE HEURISTICS
# ============================================================================

def _alignment_stem(left: str, right: str) -> tuple[int, Optional[int]]:
    """Find the longest contiguous complementary alignment."""
    left = normalize_dna_sequence(left)
    right = normalize_dna_sequence(right)
    right_rc = reverse_complement(right)
    # The stem is the longest common substring of left and the reverse
    # complement of right. The position is the one the original diagonal
    # scan reported: it walked the diagonals shift = i - j in ascending
    # order, positions within one in ascending order, and kept the first
    # run to reach a new maximum. Among all windows of the maximal length
    # that is the one with the smallest shift, then the smallest start i.
    # A common window of length k implies one of every shorter length, so
    # the maximal length is found by bisection (P2-6).
    left_length = len(left)
    right_length = len(right_rc)

    def has_common_window(size: int) -> bool:
        windows = {right_rc[j:j + size] for j in range(right_length - size + 1)}
        return any(left[i:i + size] in windows for i in range(left_length - size + 1))

    low, high = 0, min(left_length, right_length)
    while low < high:
        middle = (low + high + 1) // 2
        if has_common_window(middle):
            low = middle
        else:
            high = middle - 1
    if low == 0:
        return 0, None

    last_start: dict[str, int] = {}
    for j in range(right_length - low + 1):
        last_start[right_rc[j:j + low]] = j
    best_key: Optional[tuple[int, int]] = None
    for i in range(left_length - low + 1):
        j = last_start.get(left[i:i + low])
        if j is not None and (best_key is None or (i - j, i) < best_key):
            best_key = (i - j, i)
    assert best_key is not None
    return low, best_key[1]


def _three_prime_complementarity(
    first: str,
    second: Optional[str] = None,
    window: int = 8,
) -> int:
    """Estimate contiguous terminal complementarity as a screening heuristic."""
    first = normalize_dna_sequence(first)
    second = first if second is None else normalize_dna_sequence(second)
    first_terminal = first[-window:]
    second_terminal = second[-window:]
    # reverse_complement(second_terminal[:length]) is the last `length`
    # bases of the reverse complement of second_terminal, so it is
    # computed once.
    second_reverse = reverse_complement(second_terminal)
    best = 0
    for length in range(1, min(len(first_terminal), len(second_terminal)) + 1):
        first_part = first_terminal[-length:]
        if first_part == second_reverse[len(second_reverse) - length:]:
            best = length
    return best


# Superseded definition of analyze_hairpin() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# Superseded definition of analyze_self_dimer() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# Superseded definition of analyze_heterodimer() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# ============================================================================
# SCORING
# ============================================================================

def _distance_score(value: float, optimum: float, scale: float) -> float:
    if scale <= 0:
        return 1.0 if value == optimum else 0.0
    return max(0.0, 1.0 - abs(value - optimum) / scale)


def _category(score: float) -> str:
    if score >= 90:
        return "EXCELLENT"
    if score >= 80:
        return "GOOD"
    if score >= 70:
        return "ACCEPTABLE"
    if score >= 60:
        return "MARGINAL"
    return "POOR"


def _candidate_rejections(
    *, length: int, gc: float, tm: float, homopolymer: int,
    repeats: int, gc_clamp: int, hairpin: StructureResult,
    self_dimer: StructureResult, config: PrimerDesignConfig,
) -> list[str]:
    reasons = []
    if not config.min_length <= length <= config.max_length:
        reasons.append("primer length outside allowed range")
    if not config.min_gc <= gc <= config.max_gc:
        reasons.append("GC content outside allowed range")
    if not config.min_tm <= tm <= config.max_tm:
        reasons.append("Tm outside allowed range")
    if homopolymer > config.max_homopolymer:
        reasons.append("excessive homopolymer")
    if repeats > config.max_repeat_copies:
        reasons.append("excessive tandem repeat")
    if not config.gc_clamp_min <= gc_clamp <= config.gc_clamp_max:
        reasons.append("3' GC clamp outside preferred range")
    if hairpin.max_stem >= config.max_hairpin_stem:
        reasons.append("strong hairpin detected")
    if self_dimer.max_stem >= config.max_self_dimer_stem or self_dimer.three_prime_stem >= config.max_three_prime_complementarity:
        reasons.append("strong self-dimer / 3' complementarity detected")
    return reasons


def score_primer_candidate(
    *, length: int, gc: float, tm: float, gc_clamp: int,
    homopolymer: int, repeats: int, hairpin: StructureResult,
    self_dimer: StructureResult, config: PrimerDesignConfig,
) -> float:
    """Calculate a 100-point individual-primer score."""
    tm_score = 25.0 * _distance_score(tm, config.optimal_tm, 5.0)
    gc_score = 12.5 * _distance_score(gc, config.optimal_gc, 15.0)
    length_score = 12.5 * _distance_score(float(length), float(config.optimal_length), 5.0)
    if config.gc_clamp_min <= gc_clamp <= config.gc_clamp_max:
        clamp_score = 20.0
    elif gc_clamp == 0:
        clamp_score = 8.0
    else:
        clamp_score = 12.0
    structure_penalty = hairpin.max_stem * 2.5 + self_dimer.max_stem * 2.5 + self_dimer.three_prime_stem * 2.0
    structure_score = max(0.0, 20.0 - min(20.0, structure_penalty))
    quality_penalty = max(0, homopolymer - config.max_homopolymer) * 3.0 + max(0, repeats - config.max_repeat_copies) * 2.0
    quality_score = max(0.0, 10.0 - min(10.0, quality_penalty))
    total = tm_score + gc_score + length_score + clamp_score + structure_score + quality_score
    return round(max(0.0, min(100.0, total)), 2)


# ============================================================================
# CANDIDATE GENERATION
# ============================================================================

def _make_candidate(
    sequence: str,
    orientation: str,
    start: int,
    end: int,
    config: PrimerDesignConfig,
    tm_method: str,
) -> PrimerCandidate:
    template = sequence[start:end]
    primer = template if orientation == "FORWARD" else reverse_complement(template)
    thermodynamics = _candidate_thermodynamics(primer, config, tm_method)
    return _make_candidate_measurements(
        sequence, orientation, start, end, config, tm_method, primer, thermodynamics,
    )


def _candidate_thermodynamics(
    primer: str, config: PrimerDesignConfig, tm_method: str,
) -> "NearestNeighborTmResult | ThermodynamicResult":
    """Return the Tm record of one candidate under the selected Tm model."""
    # return_details=True always returns the detailed record.
    return cast("NearestNeighborTmResult", calculate_tm_nearest_neighbor(
        primer,
        salt_mM=config.salt_mM,
        mg_mM=config.mg_mM,
        dntp_mM=config.dntp_mM,
        primer_concentration_nM=config.primer_concentration_nM,
        include_symmetry_correction=config.include_symmetry_correction,
        return_details=True,
    )) if tm_method.lower() in NEAREST_NEIGHBOR_TM_METHODS else ThermodynamicResult(
        tm=calculate_tm(
            primer,
            tm_method,
            config.salt_mM,
            config.mg_mM,
            config.dntp_mM,
            config.primer_concentration_nM,
            config.include_symmetry_correction,
        ),
        delta_h_kcal=0.0,
        delta_s_cal=0.0,
        concentration_nM=config.primer_concentration_nM,
        salt_mM=config.salt_mM,
        mg_mM=config.mg_mM,
        dntp_mM=config.dntp_mM,
        effective_monovalent_mM=config.salt_mM,
        model=tm_method.lower(),
        symmetry_correction_applied=False,
        explanation="Empirical screening model; nearest-neighbor thermodynamic details are not applicable.",
    )


def _make_candidate_measurements(
    sequence: str,
    orientation: str,
    start: int,
    end: int,
    config: PrimerDesignConfig,
    tm_method: str,
    primer: str,
    thermodynamics,
) -> PrimerCandidate:
    tm = float(thermodynamics.tm)
    gc = gc_percent(primer)
    gc_clamp = count_gc_clamp(primer)
    homopolymer = longest_homopolymer(primer)
    repeats = longest_tandem_repeat(primer, config.max_repeat_unit)
    hairpin = analyze_hairpin(primer, config.max_hairpin_stem)
    self_dimer = analyze_self_dimer(primer, config.max_self_dimer_stem)
    rejection_reasons = _candidate_rejections(
        length=len(primer), gc=gc, tm=tm, homopolymer=homopolymer,
        repeats=repeats, gc_clamp=gc_clamp, hairpin=hairpin,
        self_dimer=self_dimer, config=config,
    )
    score = score_primer_candidate(
        length=len(primer), gc=gc, tm=tm, gc_clamp=gc_clamp,
        homopolymer=homopolymer, repeats=repeats, hairpin=hairpin,
        self_dimer=self_dimer, config=config,
    )
    return PrimerCandidate(
        sequence=primer,
        orientation=orientation,
        start=start,
        end=end,
        length=len(primer),
        gc_percent=round(gc, 2),
        tm=round(tm, 2),
        tm_model=thermodynamics.model,
        thermodynamics=thermodynamics,
        gc_clamp=gc_clamp,
        max_homopolymer=homopolymer,
        max_repeat_copies=repeats,
        hairpin=hairpin,
        self_dimer=self_dimer,
        score=score,
        category=_category(score),
        rejection_reasons=rejection_reasons,
    )


def generate_primer_candidates(
    sequence: str,
    target_start: int,
    target_end: int,
    orientation: str,
    config: PrimerDesignConfig = PrimerDesignConfig(),
    tm_method: str = "nearest_neighbor",
) -> tuple[list[PrimerCandidate], int]:
    """Generate accepted primer candidates using 0-based half-open coordinates."""
    sequence = normalize_dna_sequence(sequence)
    validate_target_region(sequence, target_start, target_end)
    orientation, windows = _candidate_windows(target_start, target_end, orientation, config)
    candidates = []
    evaluated = 0
    for start, end in windows:
        evaluated += 1
        candidate = _make_candidate(sequence, orientation, start, end, config, tm_method)
        if not candidate.rejection_reasons:
            candidates.append(candidate)
    candidates.sort(key=_generation_sort_key(config))
    return candidates[:config.max_candidates_per_orientation], evaluated


def _candidate_windows(
    target_start: int,
    target_end: int,
    orientation: str,
    config: PrimerDesignConfig,
) -> tuple[str, list[tuple[int, int]]]:
    """Return the orientation and the (start, end) windows in search order."""
    orientation = orientation.upper()
    if orientation not in {"FORWARD", "REVERSE"}:
        raise ValueError("orientation must be FORWARD or REVERSE.")
    windows = []
    if orientation == "FORWARD":
        for start in range(target_start, target_end):
            for length in range(config.min_length, config.max_length + 1):
                end = start + length
                if end > target_end:
                    break
                windows.append((start, end))
    else:
        for end in range(target_end, target_start, -1):
            for length in range(config.min_length, config.max_length + 1):
                start = end - length
                if start < target_start:
                    break
                windows.append((start, end))
    return orientation, windows


def _generation_sort_key(config: PrimerDesignConfig):
    """Order accepted candidates by score, then closeness to the optima."""
    return lambda candidate: (
        -candidate.score,
        abs(candidate.tm - config.optimal_tm),
        abs(candidate.gc_percent - config.optimal_gc),
        abs(candidate.length - config.optimal_length),
    )


# ============================================================================
# SPECIFICITY SCREENING
# ============================================================================

def screen_primer_specificity(
    primer: str,
    reference_sequences: Optional[dict[str, str] | Sequence[tuple[str, str]]] = None,
    max_mismatches: int = 2,
    three_prime_window: int = 5,
) -> list[SpecificityHit]:
    """Screen a primer against supplied references; this is not genome-wide BLAST."""
    primer = normalize_dna_sequence(primer)
    if max_mismatches < 0:
        raise ValueError("max_mismatches must be non-negative.")
    if three_prime_window <= 0:
        raise ValueError("three_prime_window must be positive.")
    if not reference_sequences:
        return []
    references = list(reference_sequences.items()) if isinstance(reference_sequences, dict) else list(reference_sequences)
    hits = []
    for name, raw_reference in references:
        reference = normalize_dna_sequence(raw_reference)
        orientations = (
            ("FORWARD_MATCH", primer),
            ("REVERSE_COMPLEMENT_MATCH", reverse_complement(primer)),
        )
        for orientation, query in orientations:
            if len(query) > len(reference):
                continue
            for start in range(len(reference) - len(query) + 1):
                window = reference[start:start + len(query)]
                mismatches = sum(a != b for a, b in zip(query, window))
                if mismatches > max_mismatches:
                    continue
                terminal_length = min(three_prime_window, len(query))
                three_prime_mismatches = sum(
                    a != b for a, b in zip(query[-terminal_length:], window[-terminal_length:])
                )
                hits.append(SpecificityHit(
                    sequence_name=name,
                    start=start,
                    end=start + len(query),
                    mismatches=mismatches,
                    three_prime_mismatches=three_prime_mismatches,
                    exact=(mismatches == 0),
                    orientation=orientation,
                ))
    hits.sort(key=lambda hit: (hit.mismatches, hit.three_prime_mismatches))
    return hits


# ============================================================================
# PRIMER PAIRS
# ============================================================================

def build_primer_pair(
    sequence: str,
    forward: PrimerCandidate,
    reverse: PrimerCandidate,
    config: PrimerDesignConfig = PrimerDesignConfig(),
    reference_sequences: Optional[dict[str, str] | Sequence[tuple[str, str]]] = None,
) -> Optional[PrimerPair]:
    """Construct and score one compatible primer pair."""
    sequence = normalize_dna_sequence(sequence)
    if forward.orientation != "FORWARD":
        raise ValueError("forward must be a FORWARD candidate.")
    if reverse.orientation != "REVERSE":
        raise ValueError("reverse must be a REVERSE candidate.")
    if forward.end > reverse.start:
        return None
    amplicon_start = forward.start
    amplicon_end = reverse.end
    amplicon_length = amplicon_end - amplicon_start
    if not config.min_amplicon_length <= amplicon_length <= config.max_amplicon_length:
        return None
    amplicon_sequence = sequence[amplicon_start:amplicon_end]
    tm_delta = abs(forward.tm - reverse.tm)
    heterodimer = analyze_heterodimer(forward.sequence, reverse.sequence, config.max_heterodimer_stem)
    specificity_hits = []
    warnings = []
    if reference_sequences:
        specificity_hits.extend(screen_primer_specificity(forward.sequence, reference_sequences))
        specificity_hits.extend(screen_primer_specificity(reverse.sequence, reference_sequences))
        if any(hit.exact for hit in specificity_hits):
            warnings.append("Exact match detected in supplied specificity references.")
        elif specificity_hits:
            warnings.append("Near-match detected in supplied specificity references.")
    primer_quality = ((forward.score + reverse.score) / 2.0) * 0.60
    tm_score = 15.0 * _distance_score(tm_delta, 0.0, config.max_tm_delta)
    amplicon_scale = max(50.0, float(config.max_amplicon_length - config.min_amplicon_length))
    amplicon_score = 10.0 * _distance_score(float(amplicon_length), float(config.optimal_amplicon_length), amplicon_scale)
    heterodimer_penalty = heterodimer.max_stem * 2.0 + heterodimer.three_prime_stem * 2.0
    heterodimer_score = max(0.0, 10.0 - min(10.0, heterodimer_penalty))
    specificity_score = 5.0
    if specificity_hits:
        specificity_score = max(0.0, 5.0 - min(5.0, sum(2.0 if hit.exact else 1.0 for hit in specificity_hits)))
    score = round(max(0.0, min(100.0, primer_quality + tm_score + amplicon_score + heterodimer_score + specificity_score)), 2)
    return PrimerPair(
        forward=forward,
        reverse=reverse,
        amplicon_start=amplicon_start,
        amplicon_end=amplicon_end,
        amplicon_length=amplicon_length,
        amplicon_sequence=amplicon_sequence,
        amplicon_gc_percent=round(gc_percent(amplicon_sequence), 2),
        tm_delta=round(tm_delta, 2),
        heterodimer=heterodimer,
        score=score,
        category=_category(score),
        specificity_hits=specificity_hits,
        warnings=warnings,
    )


def rank_primer_pairs(
    sequence: str,
    forward_candidates: Sequence[PrimerCandidate],
    reverse_candidates: Sequence[PrimerCandidate],
    config: PrimerDesignConfig = PrimerDesignConfig(),
    reference_sequences: Optional[dict[str, str] | Sequence[tuple[str, str]]] = None,
) -> tuple[list[PrimerPair], int]:
    """Evaluate compatible primer pairs and return the highest-ranked pairs."""
    pairs = []
    evaluated = 0
    for forward in forward_candidates:
        for reverse in reverse_candidates:
            if forward.end > reverse.start:
                continue
            evaluated += 1
            pair = build_primer_pair(sequence, forward, reverse, config, reference_sequences)
            if pair is not None:
                pairs.append(pair)
    pairs.sort(key=lambda pair: (
        -pair.score,
        pair.tm_delta,
        pair.heterodimer.max_stem,
        abs(pair.amplicon_length - config.optimal_amplicon_length),
    ))
    return pairs[:config.top_pairs], evaluated


# ============================================================================
# COMPLETE WORKFLOW
# ============================================================================

def design_primers(
    sequence: str,
    target_start: int = 0,
    target_end: Optional[int] = None,
    *,
    config: PrimerDesignConfig = PrimerDesignConfig(),
    tm_method: str = "nearest_neighbor",
    reference_sequences: Optional[dict[str, str] | Sequence[tuple[str, str]]] = None,
) -> PrimerDesignReport:
    """Run the complete primer-design workflow."""
    sequence = normalize_dna_sequence(sequence)
    if target_end is None:
        target_end = len(sequence)
    validate_target_region(sequence, target_start, target_end)
    method = tm_method.lower().strip()
    if method in NEAREST_NEIGHBOR_TM_METHODS:
        tm_model = "nearest_neighbor"
    elif method in {"salt", "salt_adjusted"}:
        tm_model = "salt_adjusted"
    elif method == "wallace":
        tm_model = "wallace"
    else:
        raise ValueError("Unknown Tm method.")
    forward_candidates, forward_evaluated = generate_primer_candidates(
        sequence, target_start, target_end, "FORWARD", config, tm_method
    )
    reverse_candidates, reverse_evaluated = generate_primer_candidates(
        sequence, target_start, target_end, "REVERSE", config, tm_method
    )
    ranked_pairs, pairs_evaluated = rank_primer_pairs(
        sequence, forward_candidates, reverse_candidates, config, reference_sequences
    )
    failure_reasons = []
    if not forward_candidates:
        failure_reasons.append("No acceptable forward primer candidates found.")
    if not reverse_candidates:
        failure_reasons.append("No acceptable reverse primer candidates found.")
    if forward_candidates and reverse_candidates and not ranked_pairs:
        failure_reasons.append("No acceptable primer pair satisfies amplicon and compatibility constraints.")
    warnings = []
    if not reference_sequences:
        warnings.append("Specificity was not genome-wide screened; no reference database was supplied.")
    if tm_model == "nearest_neighbor":
        warnings.append("Tm uses a SantaLucia-style nearest-neighbor model with an approximate mixed-ion correction; validate against the actual reaction chemistry.")
    elif tm_model == "wallace":
        warnings.append("Tm uses the Wallace screening rule; experimental conditions should be validated with a thermodynamic model.")
    else:
        warnings.append("Tm uses an empirical salt-adjusted model; Mg2+, dNTP and sequence-specific thermodynamics are not fully modeled.")
    return PrimerDesignReport(
        sequence_length=len(sequence),
        target_start=target_start,
        target_end=target_end,
        target_length=target_end - target_start,
        tm_model=tm_model,
        thermodynamic_conditions=get_thermodynamic_conditions(config),
        forward_candidates_evaluated=forward_evaluated,
        reverse_candidates_evaluated=reverse_evaluated,
        forward_candidates_accepted=len(forward_candidates),
        reverse_candidates_accepted=len(reverse_candidates),
        pairs_evaluated=pairs_evaluated,
        pairs_accepted=len(ranked_pairs),
        ranked_pairs=ranked_pairs,
        failure_reasons=failure_reasons,
        warnings=warnings,
    )


# ============================================================================
# REPORTING
# ============================================================================

def summarize_primer_pair(pair: PrimerPair, rank: int = 1) -> str:
    """Return a detailed human-readable primer-pair summary."""
    forward = pair.forward
    reverse = pair.reverse
    specificity = "NOT SCREENED" if not pair.specificity_hits else f"{len(pair.specificity_hits)} HIT(S)"
    return "\n".join([
        f"PRIMER PAIR #{rank}",
        f"Score:                 {pair.score:.2f}/100",
        f"Category:              {pair.category}",
        "",
        "FORWARD",
        f"Sequence:              {forward.sequence}",
        f"Position:              {forward.start + 1}-{forward.end} (1-based)",
        f"Length:                {forward.length} bp",
        f"GC:                    {forward.gc_percent:.2f}%",
        f"Tm:                    {forward.tm:.2f} °C",
        f"Tm model:              {forward.tm_model}",
        f"ΔH:                    {forward.thermodynamics.delta_h_kcal:.2f} kcal/mol",
        f"ΔS:                    {forward.thermodynamics.delta_s_cal:.2f} cal/(mol·K)",
        f"3' GC clamp:           {forward.gc_clamp} bp",
        f"Hairpin stem:          {forward.hairpin.max_stem} bp ({forward.hairpin.severity})",
        f"Self-dimer stem:       {forward.self_dimer.max_stem} bp ({forward.self_dimer.severity})",
        "",
        "REVERSE",
        f"Sequence:              {reverse.sequence}",
        f"Position:              {reverse.start + 1}-{reverse.end} (1-based)",
        f"Length:                {reverse.length} bp",
        f"GC:                    {reverse.gc_percent:.2f}%",
        f"Tm:                    {reverse.tm:.2f} °C",
        f"Tm model:              {reverse.tm_model}",
        f"ΔH:                    {reverse.thermodynamics.delta_h_kcal:.2f} kcal/mol",
        f"ΔS:                    {reverse.thermodynamics.delta_s_cal:.2f} cal/(mol·K)",
        f"3' GC clamp:           {reverse.gc_clamp} bp",
        f"Hairpin stem:          {reverse.hairpin.max_stem} bp ({reverse.hairpin.severity})",
        f"Self-dimer stem:       {reverse.self_dimer.max_stem} bp ({reverse.self_dimer.severity})",
        "",
        "PAIR",
        f"Amplicon:              {pair.amplicon_start}-{pair.amplicon_end}",
        f"Amplicon length:       {pair.amplicon_length} bp",
        f"Amplicon GC:           {pair.amplicon_gc_percent:.2f}%",
        f"ΔTm:                   {pair.tm_delta:.2f} °C",
        f"Heterodimer stem:      {pair.heterodimer.max_stem} bp ({pair.heterodimer.severity})",
        f"Specificity screening: {specificity}",
    ])


def print_primer_report(report: PrimerDesignReport) -> None:
    """Print the complete Synapse primer report."""
    style.header('PRIMER DESIGN')
    style.field('Sequence length', f"{report.sequence_length} bp")
    style.field('Target region', f"{report.target_start + 1}-{report.target_end} (1-based; {report.target_length} bp)")
    style.section("THERMODYNAMIC SETTINGS")
    style.field('Tm model', f"{report.tm_model}")
    conditions = report.thermodynamic_conditions
    style.field('Na+ / salt', f"{conditions['salt_mM']:.2f} mM")
    style.field('Mg2+', f"{conditions['mg_mM']:.2f} mM")
    style.field('dNTPs', f"{conditions['dntp_mM']:.2f} mM")
    style.field('Primer concentration', f"{conditions['primer_concentration_nM']:.1f} nM")
    style.field('Symmetry correction', f"{'YES' if conditions['include_symmetry_correction'] else 'NO'}")
    style.section("CANDIDATE GENERATION")
    style.field('Forward evaluated', f"{report.forward_candidates_evaluated}")
    style.field('Forward accepted', f"{report.forward_candidates_accepted}")
    style.field('Reverse evaluated', f"{report.reverse_candidates_evaluated}")
    style.field('Reverse accepted', f"{report.reverse_candidates_accepted}")
    style.section("PAIR EVALUATION")
    style.field('Pairs evaluated', f"{report.pairs_evaluated}")
    style.field('Pairs accepted/ranked', f"{report.pairs_accepted}")
    print()
    if report.ranked_pairs:
        for index, pair in enumerate(report.ranked_pairs, start=1):
            print(summarize_primer_pair(pair, index))
            print()
    else:
        style.section("NO ACCEPTABLE PRIMER PAIR FOUND", blank_before=False)
        print()
        if report.failure_reasons:
            style.section("MAIN LIMITING FACTORS", blank_before=False)
            for reason in report.failure_reasons:
                print(f"- {reason}")
            print()
    if report.warnings:
        style.section("WARNINGS", blank_before=False)
        for warning in report.warnings:
            print(f"- {warning}")
    print(major_rule())


# ============================================================================
# SYNAPSE API
# ============================================================================

def get_primer_summary(report: PrimerDesignReport) -> dict:
    """Convert the report into a pipeline-friendly dictionary."""
    best_pair = None
    if report.ranked_pairs:
        best = report.ranked_pairs[0]
        best_pair = {
            "forward": best.forward.sequence,
            "reverse": best.reverse.sequence,
            "forward_start": best.forward.start,
            "forward_end": best.forward.end,
            "reverse_start": best.reverse.start,
            "reverse_end": best.reverse.end,
            "amplicon_length": best.amplicon_length,
            "amplicon_gc_percent": best.amplicon_gc_percent,
            "tm_delta": best.tm_delta,
            "score": best.score,
            "category": best.category,
        }
    return {
        "sequence_length": report.sequence_length,
        "target_start": report.target_start,
        "target_end": report.target_end,
        "target_length": report.target_length,
        "tm_model": report.tm_model,
        "thermodynamic_conditions": report.thermodynamic_conditions,
        "forward_candidates_evaluated": report.forward_candidates_evaluated,
        "reverse_candidates_evaluated": report.reverse_candidates_evaluated,
        "forward_candidates_accepted": report.forward_candidates_accepted,
        "reverse_candidates_accepted": report.reverse_candidates_accepted,
        "pairs_evaluated": report.pairs_evaluated,
        "pairs_accepted": report.pairs_accepted,
        "best_score": report.ranked_pairs[0].score if report.ranked_pairs else None,
        "best_pair": best_pair,
        "failure_reasons": report.failure_reasons,
        "warnings": report.warnings,
    }


def analyze_primer(
    primer: str,
    *,
    tm_method: str = "nearest_neighbor",
    salt_mM: float = 50.0,
    mg_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_concentration_nM: float = 250.0,
    include_symmetry_correction: bool = True,
) -> dict:
    """Analyze one primer independently."""
    primer = normalize_dna_sequence(primer)
    thermo = calculate_tm_nearest_neighbor(
        primer,
        salt_mM=salt_mM,
        mg_mM=mg_mM,
        dntp_mM=dntp_mM,
        primer_concentration_nM=primer_concentration_nM,
        include_symmetry_correction=include_symmetry_correction,
        return_details=True,
    ) if tm_method.lower() in NEAREST_NEIGHBOR_TM_METHODS else None
    tm = calculate_tm(
        primer,
        tm_method,
        salt_mM,
        mg_mM,
        dntp_mM,
        primer_concentration_nM,
        include_symmetry_correction,
    )
    return {
        "sequence": primer,
        "length": len(primer),
        "gc_percent": round(gc_percent(primer), 2),
        "tm": round(tm, 2),
        "tm_model": tm_method.lower(),
        "thermodynamics": thermo,
        "gc_clamp": count_gc_clamp(primer),
        "max_homopolymer": longest_homopolymer(primer),
        "max_tandem_repeat_copies": longest_tandem_repeat(primer),
        "hairpin": analyze_hairpin(primer),
        "self_dimer": analyze_self_dimer(primer),
    }


def get_amplicon(sequence: str, forward: PrimerCandidate, reverse: PrimerCandidate) -> str:
    """Return the amplicon sequence defined by a primer pair."""
    sequence = normalize_dna_sequence(sequence)
    if forward.end > reverse.start:
        raise ValueError("Primer pair does not define a valid amplicon.")
    return sequence[forward.start:reverse.end]


# ============================================================================
# REGRESSION TESTS
# ============================================================================

def run_primer_self_test() -> dict[str, bool]:
    """Run deterministic implementation and regression tests."""
    results = {}

    assert normalize_dna_sequence(" atgc\n") == "ATGC"
    results["sequence_normalization"] = True

    assert reverse_complement("ATGC") == "GCAT"
    results["reverse_complement"] = True

    assert gc_percent("ATGC") == 50.0
    results["gc_percent"] = True

    assert calculate_tm_wallace("ATGC") == 12.0
    results["wallace_tm"] = True

    assert count_gc_clamp("ATGCGCT") == 0
    assert count_gc_clamp("ATGCGC") == 4
    assert count_gc_clamp("ATGCGCTGC") == 2
    results["gc_clamp"] = True

    assert longest_homopolymer("AATTTGC") == 3
    results["homopolymer"] = True

    assert longest_tandem_repeat("ATATATGC") >= 3
    results["tandem_repeat"] = True

    thermo = calculate_tm_nearest_neighbor(
        "GCGATCGATCGC",
        salt_mM=50,
        mg_mM=1.5,
        dntp_mM=0.2,
        primer_concentration_nM=250,
        return_details=True,
    )
    assert thermo.model == "nearest_neighbor"
    assert thermo.tm > 0
    results["nearest_neighbor_tm"] = True

    thermo_high_salt = calculate_tm_nearest_neighbor(
        "GCGATCGATCGC",
        salt_mM=100,
        mg_mM=1.5,
        dntp_mM=0.2,
        primer_concentration_nM=250,
    )
    assert thermo_high_salt > thermo.tm
    results["salt_response"] = True

    thermo_low_conc = calculate_tm_nearest_neighbor(
        "GCGATCGATCGC",
        salt_mM=50,
        mg_mM=1.5,
        dntp_mM=0.2,
        primer_concentration_nM=25,
    )
    assert thermo_low_conc < thermo.tm
    results["concentration_response"] = True

    hairpin = analyze_hairpin("GCGAAAAACGC")
    assert hairpin.max_stem >= 3
    results["hairpin_detection"] = True

    dimer = analyze_self_dimer("GCGATCGC")
    assert dimer.max_stem >= 3
    results["self_dimer_detection"] = True

    heterodimer = analyze_heterodimer("GCGATCGC", "GCGATCGC")
    assert heterodimer.max_stem >= 3
    results["heterodimer_detection"] = True

    test_sequence = (
        "TTATCGGATCTGCCGTGCGTACGACCGACACTACGCGGTTGGCTAAATAATGCTGGCGCTCGTACTCCTCCAGCGTCGATCTTCTGGTAACTAAGACTGGGACGTAGCACGGGTTGAGAAGAGTGTCTCTACGCCCCAGAAAGCGCCGTACCCTGGCTCCTTCGGATAGAGAAACGCATTGCGTTTAGGTCCATTTTATGGTGCCCGCACTGGATCCTGCGCACAGTTCATTCACATAACTGAGTCACCCAGGCTACGGATCATTCGTAAATCTCAGCTTGGGACTACTGCTTTCCTAGA"
    )

    report = design_primers(test_sequence, tm_method="nearest_neighbor")
    assert report.forward_candidates_evaluated > 0
    assert report.reverse_candidates_evaluated > 0
    results["design_workflow"] = True

    assert report.tm_model == "nearest_neighbor"
    assert report.thermodynamic_conditions["mg_mM"] == 1.5
    results["thermodynamic_reporting"] = True

    if report.ranked_pairs:
        pair = report.ranked_pairs[0]
        assert pair.amplicon_length == pair.amplicon_end - pair.amplicon_start
        assert len(pair.amplicon_sequence) == pair.amplicon_length
        assert get_amplicon(test_sequence, pair.forward, pair.reverse) == pair.amplicon_sequence
        assert 0.0 <= pair.amplicon_gc_percent <= 100.0
    results["amplicon_invariants"] = True

    hits = screen_primer_specificity("ATGC", {"reference": "CCCATGCGGG"}, max_mismatches=0)
    assert any(hit.exact for hit in hits)
    results["specificity_screen"] = True

    summary = get_primer_summary(report)
    assert "pairs_evaluated" in summary
    assert "tm_model" in summary
    assert "thermodynamic_conditions" in summary
    results["structured_api"] = True

    return results


# ============================================================================
# MODULE ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    test_results = run_primer_self_test()
    style.header('PRIMER MODULE SELF-TEST')
    for test_name, passed in test_results.items():
        print(f"[{'PASS' if passed else 'FAIL'}] {test_name}")
    print(major_rule())
    print("ALL TESTS PASSED:", all(test_results.values()))
# ============================================================================
# SYNAPSE — PRIMER DESIGN ENGINE V3 THERMODYNAMIC EXTENSION
# ============================================================================
#
# This extension deliberately preserves the complete V2 engine above.
#
# The original candidate generation, scoring, structure heuristics,
# specificity screening, pair ranking, reporting, structured API and
# regression tests remain available.
#
# V3 adds a real nearest-neighbor DNA/DNA thermodynamic screening layer
# with explicit experimental-condition parameters:
#
#     - Na+ / monovalent salt
#     - Mg2+
#     - dNTP concentration
#     - primer concentration
#
# The nearest-neighbor model is based on published DNA/DNA thermodynamic
# parameter sets commonly associated with SantaLucia-style calculations.
#
# IMPORTANT SCIENTIFIC SCOPE
# --------------------------
# This remains a primer-design screening engine, not a replacement for
# experimentally validated commercial thermodynamic software.
#
# The model accounts for duplex thermodynamics and ionic conditions at
# screening level. It does not model every laboratory effect, including
# DMSO, betaine, cosolvents, mismatches, modified bases, non-standard
# buffer systems, polymerase-specific behavior, or sequence-specific
# concentration effects.
#
# The public API remains backward compatible:
#
#     design_primers(...)
#     analyze_primer(...)
#     get_primer_summary(...)
#     get_amplicon(...)
#     print_primer_report(...)
#
# New public functions:
#
#     calculate_tm_nearest_neighbor(...)
#     calculate_duplex_thermodynamics(...)
#     calculate_effective_ionic_strength(...)
#     calculate_terminal_stability(...)
#     calculate_primer_thermodynamics(...)
#     validate_thermodynamic_conditions(...)
#     get_thermodynamic_summary(...)
#
# New configuration:
#
#     ThermodynamicConfig
#
# Existing calls that do not provide a thermodynamic configuration continue
# to use the existing V2 Wallace default.
# ============================================================================

from dataclasses import dataclass, field
from math import exp, log, log10
from typing import Optional, Sequence


# ============================================================================
# THERMODYNAMIC CONSTANTS
# ============================================================================

GAS_CONSTANT = 1.98720425864083
ABSOLUTE_ZERO_C = 273.15

# SantaLucia-style DNA/DNA nearest-neighbor parameters.
#
# ΔH values are kcal/mol.
# ΔS values are cal/(mol*K).
#
# Keys are written as dinucleotide steps on the strand being evaluated.
# Reverse-complementary steps use the corresponding symmetry-equivalent
# parameterization.
#
# The values are a standard screening set used for DNA duplex Tm work.
NEAREST_NEIGHBOR_PARAMETERS = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}


# Initiation terms.
# These are included explicitly rather than hidden inside the Tm formula.
INITIATION_ENTHALPY = 0.2
INITIATION_ENTROPY = -5.7

# Terminal AT correction.
TERMINAL_AT_ENTHALPY = 2.2
TERMINAL_AT_ENTROPY = 6.9

# Symmetry correction for self-complementary duplexes.
SYMMETRY_ENTROPY = -1.4


# ----------------------------------------------------------------------------
# Consistent nearest-neighbor parameter sets (2026-10-30, Stage 3)
# ----------------------------------------------------------------------------
# Until 2026-10-29 Synapse combined the SantaLucia (1998) stack parameters
# with the SantaLucia & Hicks (2004) initiation terms, i.e. two different
# published parameter sets. The default is now the internally consistent
# SantaLucia (1998) unified set with terminal-specific initiation, the set
# used by Primer3's recommended default (PRIMER_TM_FORMULA=1). The previous
# combination remains available as "synapse_legacy" for reproducibility.
#
# Values: (delta_H kcal/mol, delta_S cal/(mol*K)).
_SANTALUCIA_2004_STACKS = dict(NEAREST_NEIGHBOR_PARAMETERS)
_SANTALUCIA_2004_STACKS.update({
    "AA": (-7.6, -21.3),
    "TT": (-7.6, -21.3),
})

NN_PARAMETER_SETS = {
    "santalucia_1998": {
        "reference": "SantaLucia JR (1998) PNAS 95:1460-1465 (unified parameters)",
        "stacks": dict(NEAREST_NEIGHBOR_PARAMETERS),
        "initiation": (0.0, 0.0),
        "terminal": {"GC": (0.1, -2.8), "AT": (2.3, 4.1)},
        "salt_phosphates": "n_minus_1",
        "self_complementary_concentration_divisor": 1.0,
    },
    "santalucia_2004": {
        "reference": "SantaLucia JR, Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415-440",
        "stacks": _SANTALUCIA_2004_STACKS,
        "initiation": (0.2, -5.7),
        "terminal": {"GC": (0.0, 0.0), "AT": (2.2, 6.9)},
        "salt_phosphates": "n_minus_1",
        "self_complementary_concentration_divisor": 1.0,
    },
    "synapse_legacy": {
        "reference": "Synapse before 2026-10-30: 1998 stacks + 2004 initiation (mixed)",
        "stacks": dict(NEAREST_NEIGHBOR_PARAMETERS),
        "initiation": (INITIATION_ENTHALPY, INITIATION_ENTROPY),
        "terminal": {"GC": (0.0, 0.0), "AT": (TERMINAL_AT_ENTHALPY, TERMINAL_AT_ENTROPY)},
        "salt_phosphates": "n",
        "self_complementary_concentration_divisor": 2.0,
    },
}

DEFAULT_NN_PARAMETER_SET = "santalucia_1998"

# Divalent-cation handling. "monovalent_equivalent" now uses the published
# conversion of von Ahsen et al. (2001), Clin Chem 47:1956-1961, as used by
# Primer3 (oligotm.c divalent_to_monovalent) and Biopython (salt_correction):
#     [Na_eq] (mM) = [Na+] + 120 * sqrt([Mg2+] - [dNTP]),  0 if dNTP >= Mg2+.
# "synapse_legacy_4x" is the previous ad-hoc rule [Na+] + 4 * ([Mg2+] - [dNTP]).
SALT_MODELS = ("monovalent_equivalent", "na_only", "synapse_legacy_4x")


def _nn_parameter_set(name: str) -> dict:
    try:
        return NN_PARAMETER_SETS[name]
    except KeyError:
        raise ValueError(
            f"Unknown nearest-neighbor parameter set '{name}'. "
            f"Use one of: {', '.join(NN_PARAMETER_SETS)}."
        ) from None


def _salt_phosphate_count(sequence: str, parameter_set: str) -> int:
    """Number of phosphates per strand used in the SantaLucia salt term (N - 1)."""
    if _nn_parameter_set(parameter_set)["salt_phosphates"] == "n":
        return len(sequence)
    return len(sequence) - 1


# ============================================================================
# THERMODYNAMIC DATA MODELS
# ============================================================================

@dataclass(frozen=True)
class ThermodynamicConfig:
    """
    Experimental-condition configuration for primer Tm screening.

    Concentrations are expressed in common laboratory units:
        Na+      mM
        Mg2+     mM
        dNTP     mM
        primer   nM

    The defaults are deliberately conservative and are intended as a
    documented screening condition rather than a universal PCR recipe.
    """

    sodium_mM: float = 50.0
    magnesium_mM: float = 1.5
    dntp_mM: float = 0.2
    primer_nM: float = 250.0

    dna_concentration_mode: str = "duplex"
    temperature_reference_C: float = 25.0

    salt_model: str = "monovalent_equivalent"
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET

    def validate(self) -> None:
        validate_thermodynamic_conditions(self)


@dataclass(frozen=True)
class DuplexThermodynamics:
    """
    Thermodynamic result for one DNA/DNA duplex.

    ΔH is kcal/mol.
    ΔS is cal/(mol*K).
    ΔG is kcal/mol at the reported reference temperature.
    """

    sequence: str
    delta_h_kcal_mol: float
    delta_s_cal_mol_k: float
    delta_g_kcal_mol: float
    tm_C: float

    sodium_mM: float
    magnesium_mM: float
    dntp_mM: float
    primer_nM: float

    effective_salt_mM: float
    concentration_term: float

    self_complementary: bool
    terminal_at_count: int


@dataclass(frozen=True)
class PrimerThermodynamicProfile:
    """
    Thermodynamic profile used by candidate scoring.
    """

    method: str
    tm_C: float
    delta_h_kcal_mol: float
    delta_s_cal_mol_k: float
    delta_g_kcal_mol: float

    effective_salt_mM: float

    sodium_mM: float
    magnesium_mM: float
    dntp_mM: float
    primer_nM: float

    concentration_term: float
    terminal_at_count: int
    self_complementary: bool


# ============================================================================
# THERMODYNAMIC VALIDATION
# ============================================================================

def validate_thermodynamic_conditions(
    conditions: ThermodynamicConfig,
) -> None:
    """Validate all thermodynamic input parameters."""

    if conditions.sodium_mM <= 0:
        raise ValueError(
            "sodium_mM must be greater than zero."
        )

    if conditions.magnesium_mM < 0:
        raise ValueError(
            "magnesium_mM cannot be negative."
        )

    if conditions.dntp_mM < 0:
        raise ValueError(
            "dntp_mM cannot be negative."
        )

    if conditions.primer_nM <= 0:
        raise ValueError(
            "primer_nM must be greater than zero."
        )

    if conditions.temperature_reference_C <= -ABSOLUTE_ZERO_C:
        raise ValueError(
            "temperature_reference_C is below absolute zero."
        )

    if conditions.dna_concentration_mode.lower() not in {
        "duplex",
        "strand",
    }:
        raise ValueError(
            "dna_concentration_mode must be 'duplex' or 'strand'."
        )

    if conditions.salt_model.lower() not in SALT_MODELS:
        raise ValueError(
            f"salt_model must be one of: {', '.join(SALT_MODELS)}."
        )

    _nn_parameter_set(conditions.nn_parameter_set)


# ============================================================================
# EFFECTIVE IONIC CONDITIONS
# ============================================================================

def calculate_effective_ionic_strength(
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    *,
    salt_model: str = "monovalent_equivalent",
) -> float:
    """
    Estimate an effective monovalent-equivalent salt concentration.

    The purpose is to provide a transparent screening approximation when
    both Na+ and Mg2+ are supplied.

    dNTPs reduce free Mg2+ availability. A simple conservative correction is
    applied before converting Mg2+ into a monovalent-equivalent term.

    This is intentionally documented as an approximation; a complete PCR
    buffer speciation model is outside the scope of this engine.
    """

    if sodium_mM <= 0:
        raise ValueError(
            "sodium_mM must be greater than zero."
        )

    if magnesium_mM < 0:
        raise ValueError(
            "magnesium_mM cannot be negative."
        )

    if dntp_mM < 0:
        raise ValueError(
            "dntp_mM cannot be negative."
        )

    model = salt_model.lower()

    if model == "na_only":
        return float(sodium_mM)

    if model not in SALT_MODELS:
        raise ValueError(
            "Unknown salt model."
        )

    # Free Mg2+ after complexation by dNTPs (1:1 approximation).
    free_mg = max(
        0.0,
        magnesium_mM - dntp_mM,
    )

    if model == "synapse_legacy_4x":
        # Pre-2026-10-30 Synapse rule, kept only for reproducibility; it has
        # no published basis and underestimates the Mg2+ effect.
        magnesium_equivalent = 4.0 * free_mg
    else:
        # von Ahsen et al. (2001); identical to Primer3 and Biopython.
        magnesium_equivalent = 120.0 * sqrt(free_mg)

    effective = (
        float(sodium_mM)
        + magnesium_equivalent
    )

    return max(
        0.001,
        effective,
    )


def calculate_salt_correction(
    effective_salt_mM: float,
) -> float:
    """
    Calculate the logarithmic monovalent salt correction term.

    This correction is applied to duplex entropy in the nearest-neighbor
    screening model.
    """

    if effective_salt_mM <= 0:
        raise ValueError(
            "effective_salt_mM must be greater than zero."
        )

    salt_M = effective_salt_mM / 1000.0

    return 0.368 * log(
        salt_M
    )


# ============================================================================
# SEQUENCE-LEVEL THERMODYNAMIC HELPERS
# ============================================================================

def _is_self_complementary(
    sequence: str,
) -> bool:
    """Return True when a sequence equals its reverse complement."""

    sequence = normalize_dna_sequence(
        sequence
    )

    return (
        sequence
        == reverse_complement(sequence)
    )


def _terminal_at_count(
    sequence: str,
) -> int:
    """Count A/T residues at the two duplex termini."""

    sequence = normalize_dna_sequence(
        sequence
    )

    count = 0

    if sequence[0] in "AT":
        count += 1

    if sequence[-1] in "AT":
        count += 1

    return count


def _nearest_neighbor_sum(
    sequence: str,
    parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> tuple[float, float]:
    """
    Sum nearest-neighbor ΔH and ΔS values (1 M NaCl) for one parameter set,
    including initiation, terminal and symmetry terms.
    """
    # A primer's sum is needed for its Tm and again for its duplex profile;
    # it depends only on the two strings, so it is cached (P2-6).
    if isinstance(sequence, str) and isinstance(parameter_set, str):
        return _nearest_neighbor_sum_cached(sequence, parameter_set)
    return _nearest_neighbor_sum_uncached(sequence, parameter_set)


def _nearest_neighbor_sum_uncached(
    sequence: str,
    parameter_set: str,
) -> tuple[float, float]:
    sequence = normalize_dna_sequence(
        sequence
    )

    if len(sequence) < 2:
        raise ValueError(
            "Nearest-neighbor thermodynamics requires at least 2 bases."
        )

    parameters = _nn_parameter_set(parameter_set)
    stacks = parameters["stacks"]

    delta_h, delta_s = parameters["initiation"]

    for left, right in zip(
        sequence,
        sequence[1:],
    ):
        key = left + right

        if key not in stacks:
            raise ValueError(
                f"No nearest-neighbor parameter for dinucleotide {key}."
            )

        step_h, step_s = stacks[key]

        delta_h += step_h
        delta_s += step_s

    for terminal_base in (sequence[0], sequence[-1]):
        kind = "AT" if terminal_base in "AT" else "GC"
        term_h, term_s = parameters["terminal"][kind]
        delta_h += term_h
        delta_s += term_s

    if _is_self_complementary(
        sequence
    ):
        delta_s += SYMMETRY_ENTROPY

    return (
        delta_h,
        delta_s,
    )


_nearest_neighbor_sum_cached = lru_cache(maxsize=65536)(_nearest_neighbor_sum_uncached)


def calculate_duplex_concentration_term(
    primer_nM: float,
    *,
    self_complementary: bool = False,
    mode: str = "duplex",
    parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> float:
    """
    Calculate the concentration-dependent logarithmic entropy term.

    The distinction between duplex and strand concentration conventions is
    retained explicitly so the calculation is inspectable.
    """

    if primer_nM <= 0:
        raise ValueError(
            "primer_nM must be greater than zero."
        )

    mode = mode.lower()

    concentration_M = (
        primer_nM
        * 1e-9
    )

    if mode == "duplex":
        concentration_term = log(
            concentration_M / 4.0
        )
    elif mode == "strand":
        concentration_term = log(
            concentration_M
            / 2.0
        )
    else:
        raise ValueError(
            "mode must be 'duplex' or 'strand'."
        )

    if self_complementary:
        # SantaLucia (1998): CT for self-complementary duplexes (Primer3
        # "Equation A"); the legacy set used CT/2.
        divisor = _nn_parameter_set(parameter_set)[
            "self_complementary_concentration_divisor"
        ]
        concentration_term = log(
            concentration_M / divisor
        )

    return concentration_term


# ============================================================================
# NEAREST-NEIGHBOR Tm
# ============================================================================

# This is the V3 model itself. The public entry point is the
# backward-compatible adapter calculate_tm_nearest_neighbor() further down,
# which accepts the older parameter names and delegates here.
#
# Audit D16: this implementation used to carry the public name and was then
# rebound to the private name and shadowed by the adapter, so the definition a
# reader found first was not the one that ran. Naming it privately here leaves
# exactly one definition per name without changing which code executes.
def _calculate_tm_nearest_neighbor_v3_raw(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> float:
    """
    Calculate a nearest-neighbor DNA/DNA duplex Tm.

    Parameters:
        sodium_mM:
            Monovalent sodium concentration in mM.

        magnesium_mM:
            Mg2+ concentration in mM.

        dntp_mM:
            Total dNTP concentration in mM.

        primer_nM:
            Primer concentration in nM.

        dna_concentration_mode:
            'duplex' or 'strand'.

        salt_model:
            'monovalent_equivalent' or 'na_only'.

    Returns:
        Melting temperature in degrees Celsius.

    Scientific note:
        This is a sequence-level thermodynamic screening model. It is not
        intended to replace experimental calibration or a complete PCR
        buffer model.
    """

    conditions = ThermodynamicConfig(
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        dna_concentration_mode=dna_concentration_mode,
        salt_model=salt_model,
        nn_parameter_set=nn_parameter_set,
    )

    conditions.validate()

    sequence = normalize_dna_sequence(
        sequence
    )

    delta_h, delta_s = (
        _nearest_neighbor_sum(
            sequence,
            nn_parameter_set,
        )
    )

    effective_salt = (
        calculate_effective_ionic_strength(
            sodium_mM=sodium_mM,
            magnesium_mM=magnesium_mM,
            dntp_mM=dntp_mM,
            salt_model=salt_model,
        )
    )

    salt_correction = (
        calculate_salt_correction(
            effective_salt
        )
    )

    corrected_delta_s = (
        delta_s
        + _salt_phosphate_count(sequence, nn_parameter_set)
        * salt_correction
    )

    self_complementary = (
        _is_self_complementary(
            sequence
        )
    )

    concentration_term = (
        calculate_duplex_concentration_term(
            primer_nM,
            self_complementary=self_complementary,
            mode=dna_concentration_mode,
            parameter_set=nn_parameter_set,
        )
    )

    denominator = (
        corrected_delta_s
        + GAS_CONSTANT
        * concentration_term
    )

    if denominator == 0:
        raise ZeroDivisionError(
            "Thermodynamic denominator is zero."
        )

    tm_K = (
        (delta_h * 1000.0)
        / denominator
    )

    tm_C = (
        tm_K
        - ABSOLUTE_ZERO_C
    )

    return float(
        tm_C
    )


# ============================================================================
# FULL DUPLEX THERMODYNAMICS
# ============================================================================

def calculate_duplex_thermodynamics(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    reference_temperature_C: float = 25.0,
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> DuplexThermodynamics:
    """
    Return the full nearest-neighbor thermodynamic profile.
    """
    arguments = (
        sequence, sodium_mM, magnesium_mM, dntp_mM, primer_nM,
        dna_concentration_mode, salt_model, reference_temperature_C,
        nn_parameter_set,
    )
    try:
        hash(arguments)
    except TypeError:
        return _duplex_thermodynamics(*arguments)
    return _duplex_thermodynamics_cached(*arguments)


def _duplex_thermodynamics(
    sequence, sodium_mM, magnesium_mM, dntp_mM, primer_nM,
    dna_concentration_mode, salt_model, reference_temperature_C,
    nn_parameter_set,
) -> DuplexThermodynamics:
    conditions = ThermodynamicConfig(
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        dna_concentration_mode=dna_concentration_mode,
        temperature_reference_C=reference_temperature_C,
        salt_model=salt_model,
        nn_parameter_set=nn_parameter_set,
    )

    conditions.validate()

    sequence = normalize_dna_sequence(
        sequence
    )

    delta_h, delta_s = (
        _nearest_neighbor_sum(
            sequence,
            nn_parameter_set,
        )
    )

    effective_salt = (
        calculate_effective_ionic_strength(
            sodium_mM=sodium_mM,
            magnesium_mM=magnesium_mM,
            dntp_mM=dntp_mM,
            salt_model=salt_model,
        )
    )

    salt_correction = (
        calculate_salt_correction(
            effective_salt
        )
    )

    corrected_delta_s = (
        delta_s
        + _salt_phosphate_count(sequence, nn_parameter_set)
        * salt_correction
    )

    self_complementary = (
        _is_self_complementary(
            sequence
        )
    )

    concentration_term = (
        calculate_duplex_concentration_term(
            primer_nM,
            self_complementary=self_complementary,
            mode=dna_concentration_mode,
            parameter_set=nn_parameter_set,
        )
    )

    denominator = (
        corrected_delta_s
        + GAS_CONSTANT
        * concentration_term
    )

    tm_K = (
        (delta_h * 1000.0)
        / denominator
    )

    tm_C = (
        tm_K
        - ABSOLUTE_ZERO_C
    )

    reference_K = (
        reference_temperature_C
        + ABSOLUTE_ZERO_C
    )

    delta_g = (
        delta_h
        - (
            reference_K
            * corrected_delta_s
            / 1000.0
        )
    )

    return DuplexThermodynamics(
        sequence=sequence,
        delta_h_kcal_mol=round(
            delta_h,
            4,
        ),
        delta_s_cal_mol_k=round(
            corrected_delta_s,
            4,
        ),
        delta_g_kcal_mol=round(
            delta_g,
            4,
        ),
        tm_C=round(
            tm_C,
            4,
        ),
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        effective_salt_mM=round(
            effective_salt,
            4,
        ),
        concentration_term=round(
            concentration_term,
            6,
        ),
        self_complementary=self_complementary,
        terminal_at_count=_terminal_at_count(
            sequence
        ),
    )


# ============================================================================
# TERMINAL STABILITY
# ============================================================================

def calculate_terminal_stability(
    sequence: str,
    window: int = 5,
) -> dict:
    """
    Analyze the thermodynamic composition of the primer 3' terminus.

    This is complementary to the existing GC-clamp calculation. The result
    describes the terminal window rather than merely counting G/C residues.
    """

    sequence = normalize_dna_sequence(
        sequence
    )

    if window <= 0:
        raise ValueError(
            "window must be greater than zero."
        )

    terminal = sequence[
        -min(window, len(sequence)):
    ]

    gc = (
        terminal.count("G")
        + terminal.count("C")
    )

    at = (
        terminal.count("A")
        + terminal.count("T")
    )

    return {
        "window": len(terminal),
        "sequence": terminal,
        "gc_count": gc,
        "at_count": at,
        "gc_percent": (
            round(
                100.0 * gc / len(terminal),
                2,
            )
            if terminal
            else 0.0
        ),
        "gc_clamp": count_gc_clamp(
            sequence
        ),
    }


# ============================================================================
# PRIMER THERMODYNAMIC PROFILE
# ============================================================================

def calculate_primer_thermodynamics(
    sequence: str,
    *,
    conditions: ThermodynamicConfig = ThermodynamicConfig(),
) -> PrimerThermodynamicProfile:
    """
    Calculate the complete thermodynamic profile for one primer.
    """

    conditions.validate()

    duplex = (
        calculate_duplex_thermodynamics(
            sequence,
            sodium_mM=conditions.sodium_mM,
            magnesium_mM=conditions.magnesium_mM,
            dntp_mM=conditions.dntp_mM,
            primer_nM=conditions.primer_nM,
            dna_concentration_mode=(
                conditions.dna_concentration_mode
            ),
            salt_model=conditions.salt_model,
            reference_temperature_C=(
                conditions.temperature_reference_C
            ),
            nn_parameter_set=conditions.nn_parameter_set,
        )
    )

    return PrimerThermodynamicProfile(
        method="nearest_neighbor",
        tm_C=duplex.tm_C,
        delta_h_kcal_mol=(
            duplex.delta_h_kcal_mol
        ),
        delta_s_cal_mol_k=(
            duplex.delta_s_cal_mol_k
        ),
        delta_g_kcal_mol=(
            duplex.delta_g_kcal_mol
        ),
        effective_salt_mM=(
            duplex.effective_salt_mM
        ),
        sodium_mM=(
            conditions.sodium_mM
        ),
        magnesium_mM=(
            conditions.magnesium_mM
        ),
        dntp_mM=(
            conditions.dntp_mM
        ),
        primer_nM=(
            conditions.primer_nM
        ),
        concentration_term=(
            duplex.concentration_term
        ),
        terminal_at_count=(
            duplex.terminal_at_count
        ),
        self_complementary=(
            duplex.self_complementary
        ),
    )


# Primer design asks for the same duplex profile twice per candidate (once
# for the candidate's Tm, once for its thermodynamic profile). The result is
# an immutable value that depends only on the arguments, so it is cached;
# typed=True keeps 50 and 50.0 apart because the record stores the value as
# given (P2-6). An invalid input raises every time and is never cached.
_duplex_thermodynamics_cached = lru_cache(maxsize=65536, typed=True)(_duplex_thermodynamics)


# ============================================================================
# THERMODYNAMIC METHOD DISPATCH
# ============================================================================

def calculate_tm_advanced(
    sequence: str,
    *,
    method: str = "nearest_neighbor",
    thermodynamic_config: ThermodynamicConfig = (
        ThermodynamicConfig()
    ),
    salt_mM: float = 50.0,
) -> float:
    """
    Unified Tm dispatcher for the professional engine.
    """

    method = method.lower().strip()

    if method == "wallace":
        return calculate_tm_wallace(
            sequence
        )

    if method in {
        "salt",
        "salt_adjusted",
    }:
        return calculate_tm_salt_adjusted(
            sequence,
            salt_mM,
        )

    if method in {
        "nearest_neighbor",
        "nearest-neighbor",
        "nn",
        "santalucia",
    }:
        return calculate_tm_nearest_neighbor(
            sequence,
            sodium_mM=(
                thermodynamic_config.sodium_mM
            ),
            magnesium_mM=(
                thermodynamic_config.magnesium_mM
            ),
            dntp_mM=(
                thermodynamic_config.dntp_mM
            ),
            primer_nM=(
                thermodynamic_config.primer_nM
            ),
            dna_concentration_mode=(
                thermodynamic_config.dna_concentration_mode
            ),
            salt_model=(
                thermodynamic_config.salt_model
            ),
            nn_parameter_set=(
                thermodynamic_config.nn_parameter_set
            ),
        )

    raise ValueError(
        "Unknown Tm method. "
        "Use 'wallace', 'salt_adjusted', "
        "or 'nearest_neighbor'."
    )


# ============================================================================
# EXTENDED CANDIDATE DATA MODEL
# ============================================================================

@dataclass
class AdvancedPrimerCandidate:
    """
    Extended candidate model containing the original candidate plus
    thermodynamic information.
    """

    candidate: PrimerCandidate
    thermodynamics: PrimerThermodynamicProfile
    terminal_stability: dict
    thermodynamic_score: float


@dataclass
class AdvancedPrimerPair:
    """
    Extended primer-pair model containing the original pair plus
    thermodynamic compatibility.
    """

    pair: PrimerPair
    forward_thermodynamics: PrimerThermodynamicProfile
    reverse_thermodynamics: PrimerThermodynamicProfile

    thermodynamic_tm_delta: float
    thermodynamic_pair_score: float


# ============================================================================
# THERMODYNAMIC SCORING
# ============================================================================

def score_thermodynamic_profile(
    profile: PrimerThermodynamicProfile,
    *,
    optimal_tm: float = 60.0,
    min_tm: float = 57.0,
    max_tm: float = 63.0,
) -> float:
    """
    Score nearest-neighbor thermodynamics on a 25-point scale.
    """

    tm = profile.tm_C

    if tm < min_tm or tm > max_tm:
        return 0.0

    distance = abs(
        tm - optimal_tm
    )

    score = (
        25.0
        * max(
            0.0,
            1.0
            - distance / 5.0,
        )
    )

    return round(
        score,
        2,
    )


def score_thermodynamic_pair(
    forward: PrimerThermodynamicProfile,
    reverse: PrimerThermodynamicProfile,
    *,
    max_tm_delta: float = 2.0,
) -> tuple[float, float]:
    """
    Return pair thermodynamic score and ΔTm.
    """

    delta = abs(
        forward.tm_C
        - reverse.tm_C
    )

    if delta >= max_tm_delta:
        score = 0.0
    else:
        score = (
            15.0
            * max(
                0.0,
                1.0
                - delta / max_tm_delta,
            )
        )

    return (
        round(score, 2),
        round(delta, 2),
    )


# ============================================================================
# ADVANCED CANDIDATE ANALYSIS
# ============================================================================

def analyze_primer_candidate_thermodynamics(
    candidate: PrimerCandidate,
    *,
    conditions: ThermodynamicConfig = ThermodynamicConfig(),
) -> AdvancedPrimerCandidate:
    """Add thermodynamic analysis to an existing candidate."""

    profile = (
        calculate_primer_thermodynamics(
            candidate.sequence,
            conditions=conditions,
        )
    )

    terminal = calculate_terminal_stability(
        candidate.sequence
    )

    score = score_thermodynamic_profile(
        profile
    )

    return AdvancedPrimerCandidate(
        candidate=candidate,
        thermodynamics=profile,
        terminal_stability=terminal,
        thermodynamic_score=score,
    )


# ============================================================================
# ADVANCED PAIR ANALYSIS
# ============================================================================

def analyze_primer_pair_thermodynamics(
    pair: PrimerPair,
    *,
    conditions: ThermodynamicConfig = ThermodynamicConfig(),
    max_tm_delta: float = 2.0,
) -> AdvancedPrimerPair:
    """
    Add thermodynamic compatibility analysis to a primer pair.
    """

    forward = (
        calculate_primer_thermodynamics(
            pair.forward.sequence,
            conditions=conditions,
        )
    )

    reverse = (
        calculate_primer_thermodynamics(
            pair.reverse.sequence,
            conditions=conditions,
        )
    )

    pair_score, delta = (
        score_thermodynamic_pair(
            forward,
            reverse,
            max_tm_delta=max_tm_delta,
        )
    )

    return AdvancedPrimerPair(
        pair=pair,
        forward_thermodynamics=forward,
        reverse_thermodynamics=reverse,
        thermodynamic_tm_delta=delta,
        thermodynamic_pair_score=pair_score,
    )


# ============================================================================
# ADVANCED PRIMER DESIGN CONFIGURATION
# ============================================================================

@dataclass(frozen=True)
class ProfessionalPrimerConfig:
    """
    Complete professional primer-design configuration.

    This wraps the existing primer thresholds with explicit thermodynamic
    conditions without removing any existing configuration capability.
    """

    primer: PrimerDesignConfig = field(
        default_factory=PrimerDesignConfig
    )

    thermodynamics: ThermodynamicConfig = field(
        default_factory=ThermodynamicConfig
    )

    tm_method: str = "nearest_neighbor"

    terminal_window: int = 5

    require_thermodynamic_tm: bool = True

    top_pairs: int = 10
    max_pair_candidates: int = 80

    def validate(self) -> None:
        self.thermodynamics.validate()

        if self.tm_method.lower() not in {
            "wallace",
            "salt",
            "salt_adjusted",
            "nearest_neighbor",
            "nearest-neighbor",
            "nn",
            "santalucia",
        }:
            raise ValueError(
                "Unsupported Tm method."
            )

        if self.terminal_window <= 0:
            raise ValueError(
                "terminal_window must be greater than zero."
            )

        if self.top_pairs <= 0:
            raise ValueError(
                "top_pairs must be greater than zero."
            )

        if self.max_pair_candidates <= 0:
            raise ValueError(
                "max_pair_candidates must be greater than zero."
            )


# ============================================================================
# PROFESSIONAL CANDIDATE GENERATION
# ============================================================================

# Superseded definition of generate_professional_primer_candidates() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# ============================================================================
# PROFESSIONAL PAIR RANKING
# ============================================================================

# Superseded definition of rank_professional_primer_pairs() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# ============================================================================
# PROFESSIONAL REPORT MODEL
# ============================================================================

@dataclass
class ProfessionalPrimerDesignReport:
    """
    Complete report for the thermodynamic V3 workflow.
    """

    sequence_length: int

    target_start: int
    target_end: int
    target_length: int

    tm_method: str

    thermodynamic_conditions: ThermodynamicConfig

    forward_candidates_evaluated: int
    reverse_candidates_evaluated: int

    forward_candidates_accepted: int
    reverse_candidates_accepted: int

    pairs_evaluated: int
    pairs_accepted: int

    ranked_pairs: list[AdvancedPrimerPair]

    failure_reasons: list[str]
    warnings: list[str]


# ============================================================================
# PROFESSIONAL WORKFLOW
# ============================================================================

# Superseded definition of design_primers_professional() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# ============================================================================
# PROFESSIONAL REPORTING
# ============================================================================

def _print_thermodynamic_conditions(
    conditions: ThermodynamicConfig,
) -> None:
    """Print the exact thermodynamic assumptions."""

    style.section('THERMODYNAMIC CONDITIONS')

    style.field('Na+', f"{conditions.sodium_mM:.3f} mM")

    style.field('Mg2+', f"{conditions.magnesium_mM:.3f} mM")

    style.field('dNTPs', f"{conditions.dntp_mM:.3f} mM")

    style.field('Primer concentration', f"{conditions.primer_nM:.3f} nM")

    effective = (
        calculate_effective_ionic_strength(
            sodium_mM=conditions.sodium_mM,
            magnesium_mM=conditions.magnesium_mM,
            dntp_mM=conditions.dntp_mM,
            salt_model=conditions.salt_model,
        )
    )

    style.field('Effective salt', f"{effective:.3f} mM")

    style.field('Concentration mode', f"{conditions.dna_concentration_mode}")

    print(
        f"Salt model:             "
        f"{conditions.salt_model}"
        + (
            " (von Ahsen et al. 2001: Na + 120*sqrt(Mg - dNTP))"
            if conditions.salt_model == "monovalent_equivalent"
            else ""
        )
    )

    style.field('NN parameters', f"{conditions.nn_parameter_set} - "
        f"{NN_PARAMETER_SETS[conditions.nn_parameter_set]['reference']}")


# Superseded definition of summarize_professional_primer_pair() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# Superseded definition of print_professional_primer_report() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# ============================================================================
# PROFESSIONAL STRUCTURED API
# ============================================================================

def get_professional_primer_summary(
    report: ProfessionalPrimerDesignReport,
) -> dict:
    """Convert V3 results into a pipeline-friendly dictionary."""

    best_pair = None

    if report.ranked_pairs:
        best = report.ranked_pairs[0]

        best_pair = {
            "forward": best.pair.forward.sequence,
            "reverse": best.pair.reverse.sequence,
            "forward_start": (
                best.pair.forward.start
            ),
            "forward_end": (
                best.pair.forward.end
            ),
            "reverse_start": (
                best.pair.reverse.start
            ),
            "reverse_end": (
                best.pair.reverse.end
            ),
            "amplicon_length": (
                best.pair.amplicon_length
            ),
            "tm_delta": (
                best.thermodynamic_tm_delta
            ),
            "forward_tm": (
                best.forward_thermodynamics.tm_C
            ),
            "reverse_tm": (
                best.reverse_thermodynamics.tm_C
            ),
            "forward_delta_h": (
                best.forward_thermodynamics.delta_h_kcal_mol
            ),
            "reverse_delta_h": (
                best.reverse_thermodynamics.delta_h_kcal_mol
            ),
            "forward_delta_s": (
                best.forward_thermodynamics.delta_s_cal_mol_k
            ),
            "reverse_delta_s": (
                best.reverse_thermodynamics.delta_s_cal_mol_k
            ),
            "forward_delta_g": (
                best.forward_thermodynamics.delta_g_kcal_mol
            ),
            "reverse_delta_g": (
                best.reverse_thermodynamics.delta_g_kcal_mol
            ),
            "effective_salt_mM": (
                best.forward_thermodynamics.effective_salt_mM
            ),
            "thermodynamic_pair_score": (
                best.thermodynamic_pair_score
            ),
            "score": best.pair.score,
            "category": best.pair.category,
        }

    conditions = (
        report.thermodynamic_conditions
    )

    return {
        "sequence_length": (
            report.sequence_length
        ),
        "target_start": (
            report.target_start
        ),
        "target_end": (
            report.target_end
        ),
        "target_length": (
            report.target_length
        ),
        "tm_method": (
            report.tm_method
        ),
        "thermodynamic_conditions": {
            "sodium_mM": (
                conditions.sodium_mM
            ),
            "magnesium_mM": (
                conditions.magnesium_mM
            ),
            "dntp_mM": (
                conditions.dntp_mM
            ),
            "primer_nM": (
                conditions.primer_nM
            ),
            "salt_model": (
                conditions.salt_model
            ),
            "dna_concentration_mode": (
                conditions.dna_concentration_mode
            ),
            "effective_salt_mM": (
                calculate_effective_ionic_strength(
                    sodium_mM=(
                        conditions.sodium_mM
                    ),
                    magnesium_mM=(
                        conditions.magnesium_mM
                    ),
                    dntp_mM=(
                        conditions.dntp_mM
                    ),
                    salt_model=(
                        conditions.salt_model
                    ),
                )
            ),
        },
        "forward_candidates_evaluated": (
            report.forward_candidates_evaluated
        ),
        "reverse_candidates_evaluated": (
            report.reverse_candidates_evaluated
        ),
        "forward_candidates_accepted": (
            report.forward_candidates_accepted
        ),
        "reverse_candidates_accepted": (
            report.reverse_candidates_accepted
        ),
        "pairs_evaluated": (
            report.pairs_evaluated
        ),
        "pairs_accepted": (
            report.pairs_accepted
        ),
        "best_score": (
            report.ranked_pairs[0].pair.score
            if report.ranked_pairs
            else None
        ),
        "best_pair": best_pair,
        "failure_reasons": (
            report.failure_reasons
        ),
        "warnings": (
            report.warnings
        ),
    }


# ============================================================================
# PROFESSIONAL PRIMER ANALYSIS API
# ============================================================================

def analyze_primer_professional(
    primer: str,
    *,
    conditions: ThermodynamicConfig = (
        ThermodynamicConfig()
    ),
) -> dict:
    """
    Analyze a single primer with the full V3 thermodynamic layer.
    """

    primer = normalize_dna_sequence(
        primer
    )

    profile = (
        calculate_primer_thermodynamics(
            primer,
            conditions=conditions,
        )
    )

    terminal = (
        calculate_terminal_stability(
            primer
        )
    )

    return {
        "sequence": primer,
        "length": len(primer),
        "gc_percent": round(
            gc_percent(primer),
            2,
        ),
        "wallace_tm": round(
            calculate_tm_wallace(primer),
            2,
        ),
        "salt_adjusted_tm": round(
            calculate_tm_salt_adjusted(
                primer,
                conditions.sodium_mM,
            ),
            2,
        ),
        "nearest_neighbor_tm": (
            profile.tm_C
        ),
        "delta_h_kcal_mol": (
            profile.delta_h_kcal_mol
        ),
        "delta_s_cal_mol_k": (
            profile.delta_s_cal_mol_k
        ),
        "delta_g_kcal_mol": (
            profile.delta_g_kcal_mol
        ),
        "effective_salt_mM": (
            profile.effective_salt_mM
        ),
        "gc_clamp": (
            count_gc_clamp(primer)
        ),
        "terminal_stability": terminal,
        "max_homopolymer": (
            longest_homopolymer(primer)
        ),
        "max_tandem_repeat_copies": (
            longest_tandem_repeat(primer)
        ),
        "hairpin": analyze_hairpin(
            primer
        ),
        "self_dimer": analyze_self_dimer(
            primer
        ),
    }


# ============================================================================
# ADVANCED FAILURE DIAGNOSTICS
# ============================================================================

def diagnose_thermodynamic_failure(
    *,
    forward_candidates: Sequence[
        AdvancedPrimerCandidate
    ],
    reverse_candidates: Sequence[
        AdvancedPrimerCandidate
    ],
    pairs: Sequence[
        AdvancedPrimerPair
    ],
    config: ProfessionalPrimerConfig,
) -> list[str]:
    """
    Explain why thermodynamic screening rejected candidates.
    """

    reasons = []

    if not forward_candidates:
        reasons.append(
            "Nearest-neighbor Tm rejected all forward "
            "candidates."
        )

    if not reverse_candidates:
        reasons.append(
            "Nearest-neighbor Tm rejected all reverse "
            "candidates."
        )

    if (
        forward_candidates
        and reverse_candidates
        and not pairs
    ):
        reasons.append(
            "No forward/reverse combination satisfies "
            "the thermodynamic pair constraints."
        )

    if forward_candidates:
        forward_tms = [
            item.thermodynamics.tm_C
            for item in forward_candidates
        ]

        if max(forward_tms) < config.primer.min_tm:
            reasons.append(
                "Forward Tm values are consistently below "
                "the configured minimum."
            )

        if min(forward_tms) > config.primer.max_tm:
            reasons.append(
                "Forward Tm values are consistently above "
                "the configured maximum."
            )

    if reverse_candidates:
        reverse_tms = [
            item.thermodynamics.tm_C
            for item in reverse_candidates
        ]

        if max(reverse_tms) < config.primer.min_tm:
            reasons.append(
                "Reverse Tm values are consistently below "
                "the configured minimum."
            )

        if min(reverse_tms) > config.primer.max_tm:
            reasons.append(
                "Reverse Tm values are consistently above "
                "the configured maximum."
            )

    return reasons


# ============================================================================
# THERMODYNAMIC COMPARISON UTILITIES
# ============================================================================

def compare_tm_methods(
    sequence: str,
    *,
    conditions: ThermodynamicConfig = (
        ThermodynamicConfig()
    ),
) -> dict:
    """
    Compare Wallace, salt-adjusted and nearest-neighbor Tm estimates.
    """

    sequence = normalize_dna_sequence(
        sequence
    )

    wallace = (
        calculate_tm_wallace(
            sequence
        )
    )

    salt_adjusted = (
        calculate_tm_salt_adjusted(
            sequence,
            conditions.sodium_mM,
        )
    )

    nearest_neighbor = (
        calculate_tm_nearest_neighbor(
            sequence,
            sodium_mM=conditions.sodium_mM,
            magnesium_mM=conditions.magnesium_mM,
            dntp_mM=conditions.dntp_mM,
            primer_nM=conditions.primer_nM,
            dna_concentration_mode=(
                conditions.dna_concentration_mode
            ),
            salt_model=conditions.salt_model,
        )
    )

    return {
        "sequence": sequence,
        "wallace_tm_C": round(
            wallace,
            2,
        ),
        "salt_adjusted_tm_C": round(
            salt_adjusted,
            2,
        ),
        "nearest_neighbor_tm_C": round(
            nearest_neighbor,
            2,
        ),
        "wallace_minus_nn_C": round(
            wallace - nearest_neighbor,
            2,
        ),
        "salt_minus_nn_C": round(
            salt_adjusted - nearest_neighbor,
            2,
        ),
    }


# ============================================================================
# THERMODYNAMIC SELF-TESTS
# ============================================================================

def run_thermodynamic_self_test() -> dict[str, bool]:
    """
    Deterministic tests for the V3 thermodynamic layer.
    """

    results = {}

    # ------------------------------------------------------------------
    # Condition validation
    # ------------------------------------------------------------------

    conditions = ThermodynamicConfig()

    conditions.validate()

    results[
        "thermodynamic_conditions"
    ] = True

    # ------------------------------------------------------------------
    # Effective salt
    # ------------------------------------------------------------------

    effective = (
        calculate_effective_ionic_strength(
            sodium_mM=50.0,
            magnesium_mM=1.5,
            dntp_mM=0.2,
        )
    )

    assert effective > 50.0

    results[
        "effective_ionic_strength"
    ] = True

    # ------------------------------------------------------------------
    # Nearest-neighbor parameter availability
    # ------------------------------------------------------------------

    for dinucleotide in (
        "AA",
        "AT",
        "TA",
        "CG",
        "GC",
        "GG",
    ):
        assert (
            dinucleotide
            in NEAREST_NEIGHBOR_PARAMETERS
        )

    results[
        "nearest_neighbor_parameters"
    ] = True

    # ------------------------------------------------------------------
    # Nearest-neighbor Tm
    # ------------------------------------------------------------------

    test_primer = (
        "CTGGTAACTAAGACTGGGAC"
    )

    nn_tm = (
        calculate_tm_nearest_neighbor(
            test_primer
        )
    )

    assert (
        0.0
        < nn_tm
        < 100.0
    )

    results[
        "nearest_neighbor_tm"
    ] = True

    # ------------------------------------------------------------------
    # Concentration response
    # ------------------------------------------------------------------

    low_concentration = (
        calculate_tm_nearest_neighbor(
            test_primer,
            primer_nM=25.0,
        )
    )

    high_concentration = (
        calculate_tm_nearest_neighbor(
            test_primer,
            primer_nM=1000.0,
        )
    )

    assert (
        high_concentration
        > low_concentration
    )

    results[
        "primer_concentration_response"
    ] = True

    # ------------------------------------------------------------------
    # Salt response
    # ------------------------------------------------------------------

    low_salt = (
        calculate_tm_nearest_neighbor(
            test_primer,
            sodium_mM=20.0,
        )
    )

    high_salt = (
        calculate_tm_nearest_neighbor(
            test_primer,
            sodium_mM=100.0,
        )
    )

    assert (
        high_salt
        > low_salt
    )

    results[
        "salt_response"
    ] = True

    # ------------------------------------------------------------------
    # Mg2+ / dNTP response
    # ------------------------------------------------------------------

    low_mg = (
        calculate_tm_nearest_neighbor(
            test_primer,
            magnesium_mM=0.5,
            dntp_mM=0.2,
        )
    )

    high_mg = (
        calculate_tm_nearest_neighbor(
            test_primer,
            magnesium_mM=2.5,
            dntp_mM=0.2,
        )
    )

    assert (
        high_mg
        > low_mg
    )

    results[
        "magnesium_response"
    ] = True

    # ------------------------------------------------------------------
    # dNTP effect
    # ------------------------------------------------------------------

    low_dntp = (
        calculate_effective_ionic_strength(
            sodium_mM=50.0,
            magnesium_mM=2.0,
            dntp_mM=0.1,
        )
    )

    high_dntp = (
        calculate_effective_ionic_strength(
            sodium_mM=50.0,
            magnesium_mM=2.0,
            dntp_mM=1.0,
        )
    )

    assert (
        high_dntp
        < low_dntp
    )

    results[
        "dntp_response"
    ] = True

    # ------------------------------------------------------------------
    # Full duplex thermodynamics
    # ------------------------------------------------------------------

    duplex = (
        calculate_duplex_thermodynamics(
            test_primer
        )
    )

    assert (
        duplex.sequence
        == test_primer
    )

    assert (
        duplex.delta_h_kcal_mol
        < 0.0
    )

    assert (
        duplex.delta_s_cal_mol_k
        < 0.0
    )

    results[
        "duplex_thermodynamics"
    ] = True

    # ------------------------------------------------------------------
    # Primer profile
    # ------------------------------------------------------------------

    profile = (
        calculate_primer_thermodynamics(
            test_primer
        )
    )

    assert (
        profile.method
        == "nearest_neighbor"
    )

    assert (
        profile.tm_C
        == duplex.tm_C
    )

    results[
        "primer_thermodynamic_profile"
    ] = True

    # ------------------------------------------------------------------
    # Terminal stability
    # ------------------------------------------------------------------

    terminal = (
        calculate_terminal_stability(
            test_primer
        )
    )

    assert (
        terminal["window"]
        > 0
    )

    results[
        "terminal_stability"
    ] = True

    # ------------------------------------------------------------------
    # Method comparison
    # ------------------------------------------------------------------

    comparison = (
        compare_tm_methods(
            test_primer
        )
    )

    assert (
        "nearest_neighbor_tm_C"
        in comparison
    )

    results[
        "tm_method_comparison"
    ] = True

    # ------------------------------------------------------------------
    # Professional configuration
    # ------------------------------------------------------------------

    professional_config = (
        ProfessionalPrimerConfig()
    )

    professional_config.validate()

    results[
        "professional_configuration"
    ] = True

    # ------------------------------------------------------------------
    # Professional workflow
    # ------------------------------------------------------------------

    test_sequence = (
        "TTATCGGATCTGCCGTGCGTACGACCGACACTACGCGGTT"
        "GGCTAAATAATGCTGGCGCTCGTACTCCTCCAGCGTCGATC"
        "TTCTGGTAACTAAGACTGGGACGTAGCACGGGTTGAGAAGA"
        "GTGTCTCTACGCCCCAGAAAGCGCCGTACCCTGGCTCCTTC"
        "GGATAGAGAAACGCATTGCGTTTAGGTCCATTTTATGGTGC"
        "CCGCACTGGATCCTGCGCACAGTTCATTCACATAACTGAGT"
        "CACCCAGGCTACGGATCATTCGTAAATCTCAGCTTGGGACTA"
        "CTGCTTTCCTAGA"
    )

    report = design_primers_professional(
        test_sequence,
        config=ProfessionalPrimerConfig(
            primer=PrimerDesignConfig(
                min_amplicon_length=80,
                max_amplicon_length=300,
                optimal_amplicon_length=150,
                min_tm=45.0,
                max_tm=80.0,
            ),
            thermodynamics=ThermodynamicConfig(
                sodium_mM=50.0,
                magnesium_mM=1.5,
                dntp_mM=0.2,
                primer_nM=250.0,
            ),
            tm_method="nearest_neighbor",
            top_pairs=5,
        ),
    )

    assert (
        report.sequence_length
        == len(test_sequence)
    )

    assert (
        report.forward_candidates_evaluated
        > 0
    )

    assert (
        report.reverse_candidates_evaluated
        > 0
    )

    results[
        "professional_workflow"
    ] = True

    # ------------------------------------------------------------------
    # Structured professional API
    # ------------------------------------------------------------------

    summary = (
        get_professional_primer_summary(
            report
        )
    )

    assert (
        "thermodynamic_conditions"
        in summary
    )

    assert (
        "tm_method"
        in summary
    )

    results[
        "professional_structured_api"
    ] = True

    return results


# ============================================================================
# V3 TEST ENTRY POINT
# ============================================================================

def run_primer_v3_self_test() -> dict[str, bool]:
    """
    Run both the original V2 regression suite and the V3 thermodynamic suite.

    The original suite is intentionally executed first so that V3 cannot
    silently break backward compatibility.
    """

    original = run_primer_self_test()

    advanced = run_thermodynamic_self_test()

    combined = {}

    for name, passed in original.items():
        combined[
            f"legacy_{name}"
        ] = passed

    for name, passed in advanced.items():
        combined[
            f"v3_{name}"
        ] = passed

    return combined


# ============================================================================
# DOCUMENTED EXAMPLE CONFIGURATIONS
# ============================================================================

DEFAULT_PROFESSIONAL_THERMODYNAMIC_CONFIG = (
    ThermodynamicConfig(
        sodium_mM=50.0,
        magnesium_mM=1.5,
        dntp_mM=0.2,
        primer_nM=250.0,
        dna_concentration_mode="duplex",
        salt_model="monovalent_equivalent",
    )
)


HIGH_SALT_SCREENING_CONFIG = (
    ThermodynamicConfig(
        sodium_mM=100.0,
        magnesium_mM=1.5,
        dntp_mM=0.2,
        primer_nM=250.0,
        dna_concentration_mode="duplex",
        salt_model="monovalent_equivalent",
    )
)


LOW_SALT_SCREENING_CONFIG = (
    ThermodynamicConfig(
        sodium_mM=20.0,
        magnesium_mM=1.0,
        dntp_mM=0.2,
        primer_nM=250.0,
        dna_concentration_mode="duplex",
        salt_model="monovalent_equivalent",
    )
)


# ============================================================================
# EXPERIMENTAL-CONDITION SUMMARY
# ============================================================================

def get_thermodynamic_summary(
    conditions: ThermodynamicConfig,
) -> dict:
    """
    Return a structured summary of all thermodynamic assumptions.
    """

    conditions.validate()

    effective_salt = (
        calculate_effective_ionic_strength(
            sodium_mM=conditions.sodium_mM,
            magnesium_mM=conditions.magnesium_mM,
            dntp_mM=conditions.dntp_mM,
            salt_model=conditions.salt_model,
        )
    )

    return {
        "model": (
            "DNA/DNA nearest-neighbor"
        ),
        "sodium_mM": (
            conditions.sodium_mM
        ),
        "magnesium_mM": (
            conditions.magnesium_mM
        ),
        "dntp_mM": (
            conditions.dntp_mM
        ),
        "primer_nM": (
            conditions.primer_nM
        ),
        "effective_salt_mM": (
            round(
                effective_salt,
                4,
            )
        ),
        "salt_model": (
            conditions.salt_model
        ),
        "dna_concentration_mode": (
            conditions.dna_concentration_mode
        ),
        "reference_temperature_C": (
            conditions.temperature_reference_C
        ),
        "scientific_scope": (
            "Sequence-level thermodynamic screening; "
            "not a complete PCR buffer-speciation model."
        ),
    }


# ============================================================================
# PUBLIC API INVENTORY
# ============================================================================

PRIMER_V3_PUBLIC_API = (
    # Existing V2 API.
    "normalize_dna_sequence",
    "validate_target_region",
    "reverse_complement",
    "gc_percent",
    "count_gc_clamp",
    "longest_homopolymer",
    "longest_tandem_repeat",
    "calculate_tm_wallace",
    "calculate_tm_salt_adjusted",
    "calculate_tm",
    "analyze_hairpin",
    "analyze_self_dimer",
    "analyze_heterodimer",
    "generate_primer_candidates",
    "screen_primer_specificity",
    "build_primer_pair",
    "rank_primer_pairs",
    "design_primers",
    "summarize_primer_pair",
    "print_primer_report",
    "get_primer_summary",
    "analyze_primer",
    "get_amplicon",
    "run_primer_self_test",

    # V3 thermodynamic API.
    "ThermodynamicConfig",
    "DuplexThermodynamics",
    "PrimerThermodynamicProfile",
    "calculate_effective_ionic_strength",
    "calculate_salt_correction",
    "calculate_tm_nearest_neighbor",
    "calculate_duplex_thermodynamics",
    "calculate_duplex_concentration_term",
    "calculate_terminal_stability",
    "calculate_primer_thermodynamics",
    "calculate_tm_advanced",
    "ProfessionalPrimerConfig",
    "AdvancedPrimerCandidate",
    "AdvancedPrimerPair",
    "ProfessionalPrimerDesignReport",
    "generate_professional_primer_candidates",
    "rank_professional_primer_pairs",
    "design_primers_professional",
    "summarize_professional_primer_pair",
    "print_professional_primer_report",
    "get_professional_primer_summary",
    "analyze_primer_professional",
    "diagnose_thermodynamic_failure",
    "compare_tm_methods",
    "get_thermodynamic_summary",
    "run_thermodynamic_self_test",
    "run_primer_v3_self_test",
)


# ============================================================================
# V3 MODULE ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    results = run_primer_v3_self_test()

    style.header('PROFESSIONAL PRIMER ENGINE V3 SELF-TEST')

    for name, passed in results.items():
        print(
            f"[{'PASS' if passed else 'FAIL'}] "
            f"{name}"
        )

    print(major_rule())

    print(
        "TOTAL TESTS:",
        len(results),
    )

    print(
        "ALL TESTS PASSED:",
        all(results.values()),
    )

    print(major_rule())
# ============================================================================
# BACKWARD-COMPATIBILITY ADAPTER
# ============================================================================
#
# The original Synapse primer engine already exposed a nearest-neighbor API
# with legacy keyword names. The adapter below preserves those names while
# keeping the new V3 thermodynamic implementation underneath.
# ============================================================================

@dataclass(frozen=True)
class NearestNeighborTmResult:
    """Backward-compatible detailed nearest-neighbor Tm result."""

    tm: float
    model: str
    sequence: str
    delta_h_kcal_mol: float
    delta_s_cal_mol_k: float
    effective_salt_mM: float
    sodium_mM: float
    magnesium_mM: float
    dntp_mM: float
    primer_concentration_nM: float
    symmetry_correction_applied: bool = False
    explanation: str = (
        "Nearest-neighbor DNA/DNA thermodynamic screening model."
    )

    @property
    def delta_h_kcal(self) -> float:
        return self.delta_h_kcal_mol

    @property
    def delta_s_cal(self) -> float:
        return self.delta_s_cal_mol_k

    @property
    def concentration_nM(self) -> float:
        return self.primer_concentration_nM

    @property
    def salt_mM(self) -> float:
        return self.sodium_mM

    @property
    def mg_mM(self) -> float:
        return self.magnesium_mM


def calculate_tm_nearest_neighbor(
    sequence: str,
    *,
    sodium_mM: float = 50.0,
    magnesium_mM: float = 1.5,
    dntp_mM: float = 0.2,
    primer_nM: float = 250.0,
    dna_concentration_mode: str = "duplex",
    salt_model: str = "monovalent_equivalent",
    salt_mM: Optional[float] = None,
    mg_mM: Optional[float] = None,
    primer_concentration_nM: Optional[float] = None,
    include_symmetry_correction: bool = True,
    return_details: bool = False,
    nn_parameter_set: str = DEFAULT_NN_PARAMETER_SET,
) -> float | NearestNeighborTmResult:
    """
    Backward-compatible nearest-neighbor Tm API.

    Supported parameter names include both the V3 names and the names used
    by earlier Synapse primer-engine tests:

        sodium_mM
        magnesium_mM / mg_mM
        dntp_mM
        primer_nM / primer_concentration_nM
        salt_mM
        return_details
    """

    if salt_mM is not None:
        sodium_mM = float(salt_mM)

    if mg_mM is not None:
        magnesium_mM = float(mg_mM)

    if primer_concentration_nM is not None:
        primer_nM = float(
            primer_concentration_nM
        )

    if not include_symmetry_correction:
        # The symmetry term only matters for self-complementary oligos; it
        # cannot be switched off in the consistent parameter sets.
        import logging
        logging.getLogger("synapse.primer").info(
            "include_symmetry_correction=False is ignored; the symmetry term "
            "is part of the nearest-neighbor model."
        )

    tm = _calculate_tm_nearest_neighbor_v3_raw(
        sequence,
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        dna_concentration_mode=dna_concentration_mode,
        salt_model=salt_model,
        nn_parameter_set=nn_parameter_set,
    )

    if not return_details:
        return tm

    duplex = calculate_duplex_thermodynamics(
        sequence,
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_nM=primer_nM,
        dna_concentration_mode=dna_concentration_mode,
        salt_model=salt_model,
        nn_parameter_set=nn_parameter_set,
    )

    return NearestNeighborTmResult(
        tm=duplex.tm_C,
        model="nearest_neighbor",
        sequence=duplex.sequence,
        delta_h_kcal_mol=(
            duplex.delta_h_kcal_mol
        ),
        delta_s_cal_mol_k=(
            duplex.delta_s_cal_mol_k
        ),
        effective_salt_mM=(
            duplex.effective_salt_mM
        ),
        sodium_mM=sodium_mM,
        magnesium_mM=magnesium_mM,
        dntp_mM=dntp_mM,
        primer_concentration_nM=primer_nM,
        symmetry_correction_applied=(
            duplex.self_complementary
            and include_symmetry_correction
        ),
        explanation=(
            "Nearest-neighbor DNA/DNA thermodynamic screening "
            "with ionic and concentration corrections."
        ),
    )


# ============================================================================
# FINAL V3 SELF-TEST ENTRY POINT
# ============================================================================

def run_primer_v3_self_test_final() -> dict[str, bool]:
    """
    Final compatibility suite.

    Runs the legacy suite after the compatibility adapter has been installed,
    then runs all V3 thermodynamic tests.
    """

    legacy = run_primer_self_test()
    v3 = run_thermodynamic_self_test()

    combined = {}

    for name, passed in legacy.items():
        combined[
            f"legacy_{name}"
        ] = passed

    for name, passed in v3.items():
        combined[
            f"v3_{name}"
        ] = passed

    return combined


# ============================================================================
# SYNAPSE — PRIMER DESIGN UI
# ============================================================================
#
# Module-owned UI.
# The primer engine above is intentionally preserved.
# ============================================================================

def _primer_ui_store_result(pipeline, result):
    """Store a primer result when the active pipeline supports it."""
    if pipeline is None:
        return

    import logging
    logger = logging.getLogger("synapse.primer")
    try:
        context = getattr(pipeline, "context", None)
        store = getattr(context, "store_primer_analysis", None)
        if callable(store):
            store(result)
        else:
            logger.warning("Pipeline has no store_primer_analysis(); primer result not stored.")
    except Exception:
        logger.exception("Could not store primer result in the pipeline.")


def _primer_ui_read_target(sequence):
    """Read an optional target region."""
    style.section('TARGET REGION')
    print("Enter 1-based, inclusive positions (the first base is 1).")
    print("Press Enter to use the complete sequence.")
    print()

    start_input = style.ask('Target start [1]')
    end_input = style.ask(f"Target end [{len(sequence)}]")

    try:
        first = 1 if not start_input else int(start_input)
        last = len(sequence) if not end_input else int(end_input)
        if not 1 <= first <= last <= len(sequence):
            raise ValueError(
                f"Positions must satisfy 1 <= start <= end <= {len(sequence)}."
            )
        # The design functions take a 0-based, half-open interval.
        target_start, target_end = first - 1, last
        validate_target_region(sequence, target_start, target_end)
        return target_start, target_end
    except (TypeError, ValueError) as exc:
        style.error(f"{exc}")
        return None


def _primer_ui_read_conditions():
    """Read thermodynamic conditions using documented defaults."""
    style.section('THERMODYNAMIC CONDITIONS')
    print("Press Enter to keep the default value.")
    print()

    try:
        sodium_input = style.ask('Na+ / salt [50.0 mM]')
        mg_input = style.ask('Mg2+ [1.5 mM]')
        dntp_input = style.ask('dNTPs [0.2 mM]')
        primer_input = style.ask('Primer concentration [250 nM]')

        sodium_mM = 50.0 if not sodium_input else float(sodium_input)
        magnesium_mM = 1.5 if not mg_input else float(mg_input)
        dntp_mM = 0.2 if not dntp_input else float(dntp_input)
        primer_nM = 250.0 if not primer_input else float(primer_input)

        conditions = ThermodynamicConfig(
            sodium_mM=sodium_mM,
            magnesium_mM=magnesium_mM,
            dntp_mM=dntp_mM,
            primer_nM=primer_nM,
        )
        conditions.validate()
        return conditions
    except (TypeError, ValueError) as exc:
        style.error(f"{exc}")
        return None


def primer_pair_binding_sites(pair, three_prime_window: int = 5) -> list[dict]:
    """The two binding sites of a primer pair, 1-based inclusive (audit P2-4).

    Accepts a PrimerPair or an AdvancedPrimerPair. A primer's start/end are
    0-based half-open on the template; the forward primer's 3' end is the
    right end of its site, the reverse primer's 3' end the left end (it binds
    the top strand reverse-complemented). three_prime_window bases at the 3'
    end form the critical region, where a mismatch most affects extension.
    """
    pair = getattr(pair, "pair", pair)
    window = max(int(three_prime_window), 1)
    forward, reverse = pair.forward, pair.reverse
    f_start, f_end = forward.start + 1, forward.end
    r_start, r_end = reverse.start + 1, reverse.end
    return [
        {"name": "forward primer", "kind": "primer", "start": f_start, "end": f_end,
         "critical": (max(f_start, f_end - window + 1), f_end), "critical_label": "3' end"},
        {"name": "reverse primer", "kind": "primer", "start": r_start, "end": r_end,
         "critical": (r_start, min(r_end, r_start + window - 1)), "critical_label": "3' end"},
    ]


def _primer_ui_select_tm_model():
    """Select the Tm model."""
    methods = {
        "": "nearest_neighbor",  # Enter takes the [1] the prompt shows
        "1": "nearest_neighbor",
        "2": "wallace",
        "3": "salt_adjusted",
    }
    while True:
        style.menu([
            ("Tm model", [('1', 'Nearest-neighbor'), ('2', 'Wallace'), ('3', 'Salt-adjusted')]),
        ], zero="Cancel")
        choice = style.ask('Selection [1]')
        if choice == "0":
            return None
        if choice in methods:
            return methods[choice]
        style.invalid_selection()


# Superseded definition of _primer_ui_design() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


def _primer_ui_single_analysis(pipeline=None):
    """Analyze one existing primer sequence."""
    style.header("Single primer analysis")

    primer_input = style.ask('Primer sequence')

    if not primer_input:
        style.error("Primer sequence is empty.")
        return None

    try:
        primer = normalize_dna_sequence(primer_input)
    except (TypeError, ValueError) as exc:
        style.error(f"{exc}")
        return None

    method = _primer_ui_select_tm_model()
    if method is None:
        return None

    conditions = _primer_ui_read_conditions()

    if conditions is None:
        return None

    try:
        if method == "nearest_neighbor":
            result = analyze_primer_professional(
                primer,
                conditions=conditions,
            )

            style.header('PRIMER ANALYSIS')
            style.field('Sequence', f"{result['sequence']}")
            style.field('Length', f"{result['length']} bp")
            style.field('GC content', f"{result['gc_percent']:.2f}%")
            style.field('Nearest-neighbor Tm', f"{result['nearest_neighbor_tm']:.2f} °C")
            style.field('Wallace Tm', f"{result['wallace_tm']:.2f} °C")
            style.field('Salt-adjusted Tm', f"{result['salt_adjusted_tm']:.2f} °C")
            print(f"ΔH:                     {result['delta_h_kcal_mol']:.2f} kcal/mol")
            print(f"ΔS:                     {result['delta_s_cal_mol_k']:.2f} cal/(mol·K)")
            print(f"ΔG at 25 °C:            {result['delta_g_kcal_mol']:.2f} kcal/mol")
            style.field('Effective salt', f"{result['effective_salt_mM']:.2f} mM")
            print(f"3' GC clamp:            {result['gc_clamp']} bp")
            style.field('Longest homopolymer', f"{result['max_homopolymer']}")
            style.field('Longest tandem repeat', f"{result['max_tandem_repeat_copies']}")
            style.field('Hairpin stem', f"{result['hairpin'].max_stem} bp "
                f"({result['hairpin'].severity})")
            style.field('Self-dimer stem', f"{result['self_dimer'].max_stem} bp "
                f"({result['self_dimer'].severity})")
            print(major_rule())

        else:
            result = analyze_primer(
                primer,
                tm_method=method,
                salt_mM=conditions.sodium_mM,
                mg_mM=conditions.magnesium_mM,
                dntp_mM=conditions.dntp_mM,
                primer_concentration_nM=conditions.primer_nM,
            )

            style.header('PRIMER ANALYSIS')
            style.field('Sequence', f"{result['sequence']}")
            style.field('Length', f"{result['length']} bp")
            style.field('GC content', f"{result['gc_percent']:.2f}%")
            style.field('Tm', f"{result['tm']:.2f} °C")
            style.field('Tm model', f"{result['tm_model']}")
            print(f"3' GC clamp:            {result['gc_clamp']} bp")
            style.field('Longest homopolymer', f"{result['max_homopolymer']}")
            style.field('Longest tandem repeat', f"{result['max_tandem_repeat_copies']}")
            style.field('Hairpin stem', f"{result['hairpin'].max_stem} bp "
                f"({result['hairpin'].severity})")
            style.field('Self-dimer stem', f"{result['self_dimer'].max_stem} bp "
                f"({result['self_dimer'].severity})")
            print(major_rule())

        return result

    except (TypeError, ValueError, ZeroDivisionError) as exc:
        style.error(f"{exc}")
        return None


# Superseded definition of run_primer_analysis_ui() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# ============================================================================
# SYNAPSE — PROFESSIONAL PRIMER ENGINE V4 SCREENING OVERRIDES
# ============================================================================
#
# The V2/V3 engine above is intentionally preserved.  The functions below
# replace only the professional candidate-selection layer and the CLI entry
# point.  The thermodynamic equations, legacy API and existing analysis
# routines remain available.
#
# Main improvements:
# - The professional NN workflow no longer pre-filters candidates with the
#   Wallace Tm model.
# - 3' GC clamp is treated as a scoring preference, not an automatic failure.
# - Hairpin/self-dimer heuristics are used as graded penalties; only clearly
#   severe structures are hard-rejected in the professional workflow.
# - Candidate search uses the complete configured length range before ranking.
# - If the strict 57–63 °C window contains no candidates, a transparent
#   fallback tier searches 55–65 °C and then reports the relaxation.
# - The report explains why candidates were rejected instead of simply saying
#   "0 candidates".
# - Existing legacy generate_primer_candidates() behavior is untouched.
# ============================================================================

from dataclasses import replace as _dataclass_replace
from collections import Counter as _Counter


# Superseded definition of generate_professional_primer_candidates() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


# Superseded definition of design_primers_professional() removed in Stage 4 (never reachable;
# a later definition in this module replaces it). See docs/CHANGES.md.


def print_professional_primer_report(
    report: ProfessionalPrimerDesignReport,
) -> None:
    """Print the improved professional report with explicit screening status."""
    style.header('PROFESSIONAL PRIMER DESIGN')
    style.field('Sequence length', f"{report.sequence_length} bp")
    style.field('Target region', f"{report.target_start + 1}-{report.target_end} "
        f"(1-based; {report.target_length} bp)")
    style.field('Tm method', f"{report.tm_method}")
    _print_thermodynamic_conditions(report.thermodynamic_conditions)

    style.section('CANDIDATE GENERATION')
    style.field('Forward evaluated', f"{report.forward_candidates_evaluated}")
    style.field('Forward selected', f"{report.forward_candidates_accepted}")
    style.field('Reverse evaluated', f"{report.reverse_candidates_evaluated}")
    style.field('Reverse selected', f"{report.reverse_candidates_accepted}")

    style.section('SCREENING POLICY')
    print("Tier 1: strict preferred criteria")
    print("Tier 2: relaxed Tm window (+/- 2 °C)")
    print("Tier 3: best measured fallback candidates")
    print("3' GC clamp and secondary structure are scored, not blindly discarded.")

    style.section('PAIR EVALUATION')
    style.field('Pairs evaluated', f"{report.pairs_evaluated}")
    style.field('Pairs accepted/ranked', f"{report.pairs_accepted}")

    if report.ranked_pairs:
        print()
        for index, advanced_pair in enumerate(report.ranked_pairs, start=1):
            print(summarize_professional_primer_pair(advanced_pair, index))
            print()
    else:
        print()
        style.warning("No acceptable primer pair was found.")
        if report.failure_reasons:
            style.section("Main limiting factors")
            for reason in report.failure_reasons:
                print(f"- {reason}")

    if report.warnings:
        style.section("Warnings and diagnostics")
        for warning in report.warnings:
            style.warning(warning)

    print(major_rule())


def _primer_ui_design_menu_v4(dna, header, pipeline):
    """Improved compact Primer Design module menu."""
    sequence = normalize_dna_sequence(dna)

    while True:
        style.header("Primer design")
        style.field('Input', f"{header}")
        style.field('Length', f"{len(sequence)} bp")
        style.menu([
            ("Design", [('1', 'Design Primers'),
                        ('2', 'Design Primers Around a Restriction Site')]),
            ("Analysis", [('3', 'Single Primer Analysis')]),
        ])

        choice = style.ask_choice()

        if choice == "0":
            return None

        if choice == "1":
            target = _primer_ui_read_target(sequence)
            if target is None:
                continue
            result = _primer_ui_design(
                sequence,
                target[0],
                target[1],
                pipeline,
            )
            if result is not None:
                return result
            continue

        if choice == "3":
            result = _primer_ui_single_analysis(pipeline)
            if result is not None:
                return result
            continue

        if choice == "2":
            result = _primer_ui_restriction_site_design(sequence, pipeline)
            if result is not None:
                return result
            continue

        style.invalid_selection()


def run_primer_analysis_ui(
    dna,
    pipeline=None,
    header="Sequence",
):
    """Module-owned Primer Design UI entry point."""
    try:
        normalize_dna_sequence(dna)
    except (TypeError, ValueError) as exc:
        style.header("Primer design")
        style.error(f"{exc}")
        return None

    return _primer_ui_design_menu_v4(
        dna,
        header,
        pipeline,
    )


# ============================================================================
# V4 REGRESSION CHECKS
# ============================================================================

def run_primer_v4_screening_self_test() -> dict[str, bool]:
    """Test the improved professional screening layer."""
    results = {}

    sequence = (
        "ATGCGTACGTTAGCTAGCGATCGTACGATCGGCTAGCTAGG"
        "ATCGCTAGCTGATCGATCGGATCGCTAGCTAGCTACGATCG"
        "ATCGGCTAGCTAGCTGATCGATCGGCTAGCTAGCTAGCTAA"
        "GATCGCTAGCTAGCTGATCGATCGGATCGCTAGCTAGCTGA"
    )

    config = ProfessionalPrimerConfig(
        primer=PrimerDesignConfig(
            min_amplicon_length=40,
            max_amplicon_length=200,
            optimal_amplicon_length=120,
        ),
        thermodynamics=ThermodynamicConfig(),
        tm_method="nearest_neighbor",
        top_pairs=5,
    )

    report = design_primers_professional(
        sequence,
        config=config,
    )

    assert report.forward_candidates_evaluated > 0
    assert report.reverse_candidates_evaluated > 0
    results["professional_candidate_generation"] = True

    if report.ranked_pairs:
        best = report.ranked_pairs[0]
        assert best.pair.forward.orientation == "FORWARD"
        assert best.pair.reverse.orientation == "REVERSE"
        assert best.pair.amplicon_length > 0
    results["pair_invariants"] = True

    primer = "ATGCGTACGTTAGCTAGCGA"
    analysis = analyze_primer_professional(
        primer,
        conditions=ThermodynamicConfig(),
    )
    assert analysis["length"] == 20
    assert 0.0 <= analysis["gc_percent"] <= 100.0
    results["single_primer_analysis"] = True

    return results


# ============================================================================
# V5 PROFESSIONAL PRIMER SCREENING / RANKING OVERRIDES
# ============================================================================
#
# V5 keeps the existing Synapse primer engine intact and changes only the
# professional screening/ranking layer:
#   - strong heterodimers receive a substantially larger penalty
#   - strong 3' self-complementarity receives a larger penalty
#   - strong hairpins receive a larger penalty
#   - Wallace and salt-adjusted modes use the same tiered candidate strategy
#   - candidate ranking is quality-first rather than Tm-first
#   - "ACCEPTABLE" is reserved for pairs without high-risk structural flags
#
# The underlying thermodynamic calculations and legacy APIs remain intact.
#
# Reading guide for the _vN_ prefixes (audit D16).
#
# These names record the generation a helper was introduced in; they are not
# a chain of live alternatives. Each public name is bound exactly once, at the
# end of the section that defines its final implementation:
#
#   generate_professional_primer_candidates -> _v7_generate_professional_primer_candidates
#   rank_professional_primer_pairs          -> _v11_rank_professional_primer_pairs
#   _primer_ui_design                       -> _v7_primer_ui_design
#   design_primers_professional             -> the plain def near the end of the module
#
# Earlier generations used to rebind the same public names in sequence, so a
# reader had to scan the whole module to learn which implementation actually
# ran. Those superseded bindings and the 24 private definitions they alone
# kept reachable were removed in Phase 6; every surviving helper below is
# reached from one of the four entry points above. Helpers that still carry an
# older prefix (_v5_structure_penalty, _v6_relaxed_generation_config,
# _v8_structure_severity) are live and shared - the prefix is history, not a
# version switch.
# ============================================================================

def _v5_structure_penalty(stem: int, low: int, high: int) -> float:
    """Return a graded penalty for a structural stem length."""
    if stem <= low:
        return 0.0
    if stem >= high:
        return 30.0
    return 30.0 * (stem - low) / float(high - low)


def _v5_candidate_quality(item) -> float:
    """Score a candidate while strongly penalizing structural risk."""
    profile = item.thermodynamics
    candidate = item.candidate

    # The stems are measured on the candidate. Before 2026-11-12 they were
    # read from the thermodynamic profile, which has no such fields, so the
    # structural penalty was always zero.
    hairpin = int(candidate.hairpin.max_stem or 0)
    self_dimer = int(candidate.self_dimer.max_stem or 0)

    # Prefer primers close to the configured optimum, but structure dominates.
    tm_error = abs(
        float(profile.tm_C) - float(_V5_ACTIVE_CONFIG.primer.optimal_tm)
    )
    gc_error = abs(
        float(candidate.gc_percent) - float(_V5_ACTIVE_CONFIG.primer.optimal_gc)
    )

    penalty = (
        _v5_structure_penalty(hairpin, 3, 8)
        + _v5_structure_penalty(self_dimer, 4, 9)
    )

    return (
        float(item.thermodynamic_score) * 2.0
        + float(candidate.score)
        - penalty
        - tm_error * 2.0
        - gc_error
    )


# Bind the V5 implementations.
_V5_ACTIVE_CONFIG = ProfessionalPrimerConfig()


# ============================================================================
# V6 CANDIDATE POOL FIX
# ============================================================================
# V5 still delegated candidate generation to the legacy acceptance filter.
# That filter intentionally rejects hairpins, self-dimers and 3' GC-clamp
# deviations, so a relaxed Tm window alone could still produce zero
# candidates. V6 separates:
#   1) candidate enumeration/measurement
#   2) professional screening/ranking
#
# This is the correct architecture for the professional workflow. The legacy
# generate_primer_candidates() behavior remains unchanged for compatibility.
# ============================================================================

def _v6_relaxed_generation_config(base_config: PrimerDesignConfig) -> PrimerDesignConfig:
    """Create a measurement-only configuration with no hard quality filters."""
    return replace(
        base_config,
        min_gc=0.0,
        max_gc=100.0,
        min_tm=0.0,
        max_tm=100.0,
        max_homopolymer=999,
        max_repeat_copies=999,
        gc_clamp_min=0,
        gc_clamp_max=999,
        max_hairpin_stem=999,
        max_self_dimer_stem=999,
        max_three_prime_complementarity=999,
        max_candidates_per_orientation=100000,
    )




# ============================================================================
# V7 MULTI-MODE PROFESSIONAL WORKFLOW
# ============================================================================
# V6 fixed the candidate pool, but the CLI still routed Wallace and
# salt-adjusted modes through the legacy hard-filter workflow. V7 makes all
# three selectable Tm models use the same professional screening pipeline.
#
# For Wallace/salt-adjusted modes, the selected model supplies the displayed
# Tm and Tm-pair compatibility, while the existing nearest-neighbor
# thermodynamic profile remains available as supplemental structural
# thermodynamic information.
# ============================================================================

def _v7_apply_selected_tm_model(
    advanced: AdvancedPrimerCandidate,
    candidate: PrimerCandidate,
    config: ProfessionalPrimerConfig,
) -> AdvancedPrimerCandidate:
    """Make the selected Tm model authoritative for professional screening."""
    method = config.tm_method.lower().strip()
    if method in NEAREST_NEIGHBOR_TM_METHODS:
        return advanced

    selected_tm = float(candidate.tm)
    profile = advanced.thermodynamics

    advanced.thermodynamics = replace(
        profile,
        method=method,
        tm_C=round(selected_tm, 2),
    )

    # Recalculate the thermodynamic screening score using the selected model's
    # Tm while retaining the measured structural/NN information.
    advanced.thermodynamic_score = score_thermodynamic_profile(
        advanced.thermodynamics,
        optimal_tm=config.primer.optimal_tm,
        min_tm=config.primer.min_tm,
        max_tm=config.primer.max_tm,
    )
    return advanced


def _v7_generate_professional_primer_candidates(
    sequence: str,
    target_start: int,
    target_end: int,
    orientation: str,
    *,
    config=None,
):
    """Generate a complete professional candidate pool for any Tm model."""
    global _V5_ACTIVE_CONFIG

    config = config or ProfessionalPrimerConfig()
    config.validate()
    _V5_ACTIVE_CONFIG = config

    base_config = config.primer
    relaxed = _v6_relaxed_generation_config(base_config)

    limit = int(getattr(config, "max_pair_candidates", 80))

    # The candidates are ordered by their Tm tier first, and only the best
    # `limit` are kept. The tier needs only the candidate's Tm, so it is
    # computed for every window first; the full measurement (structures,
    # scores) then runs tier by tier and stops once the better tiers hold
    # at least `limit` accepted candidates, because every candidate of a
    # later tier sorts after them. The kept candidates go through the same
    # two stable sorts in the same search order as before, so the result
    # is identical to measuring every window (P2-6: at 10 kb most windows
    # lie outside the Tm window).
    sequence = normalize_dna_sequence(sequence)
    validate_target_region(sequence, target_start, target_end)
    orientation, windows = _candidate_windows(target_start, target_end, orientation, relaxed)
    if len(windows) > relaxed.max_candidates_per_orientation:
        # The generator would truncate the pool before the tier sort; keep
        # its exact behavior by measuring every window.
        return _v7_measure_all_candidates(
            sequence, target_start, target_end, orientation, config, relaxed, limit,
        )

    windows_by_tier: dict[int, list[int]] = {}
    profiles: dict[int, PrimerThermodynamicProfile] = {}
    for index, (start, end) in enumerate(windows):
        tm, profile = _v7_screening_tm(sequence, orientation, start, end, config, relaxed)
        if profile is not None:
            profiles[index] = profile
        windows_by_tier.setdefault(_v7_tm_tier(tm, base_config), []).append(index)

    # Within a tier the next sort key is the screening score. Without the
    # hairpin and self-dimer measurements it has an upper bound (a primer
    # without any stem, and no structural penalty). Windows are measured in
    # descending order of that bound, and the tier is left as soon as it
    # has filled the remaining places with candidates whose true score is
    # strictly greater than the bound of every window not yet measured:
    # those windows sort strictly after them.
    accepted: list[tuple[int, PrimerCandidate]] = []
    for tier in sorted(windows_by_tier):
        places = limit - len(accepted)
        bounded = sorted(
            (
                (-_v7_screening_upper_bound(
                    sequence, orientation, windows[index][0], windows[index][1],
                    config, relaxed, profiles.get(index),
                ), index)
                for index in windows_by_tier[tier]
            ),
        )
        tier_scores: list[float] = []
        for negative_bound, index in bounded:
            if len(tier_scores) >= places:
                tier_scores.sort(reverse=True)
                if tier_scores[places - 1] > -negative_bound:
                    break
            start, end = windows[index]
            candidate = _make_candidate(
                sequence, orientation, start, end, relaxed, config.tm_method,
            )
            if not candidate.rejection_reasons:
                accepted.append((index, candidate))
                tier_scores.append(_v7_screening_score(candidate, config))
        if len(accepted) >= limit:
            break

    accepted.sort(key=lambda item: item[0])
    candidates = [candidate for _, candidate in accepted]
    candidates.sort(key=_generation_sort_key(relaxed))
    measured = _v7_measure_candidates(candidates, config)
    return measured[:limit], len(windows)


_V7_NEAREST_NEIGHBOR_METHODS = NEAREST_NEIGHBOR_TM_METHODS


def _v7_screening_tm(sequence, orientation, start, end, config, relaxed):
    """Return the Tm that decides a window's screening tier.

    It is the value _v7_apply_selected_tm_model leaves in the candidate's
    thermodynamic profile: the nearest-neighbor profile Tm, or the selected
    model's Tm rounded as the candidate stores it. The nearest-neighbor
    profile is returned too (None for the other models) so that it is not
    computed again.
    """
    template = sequence[start:end]
    primer = template if orientation == "FORWARD" else reverse_complement(template)
    if config.tm_method.lower().strip() in _V7_NEAREST_NEIGHBOR_METHODS:
        profile = calculate_primer_thermodynamics(primer, conditions=config.thermodynamics)
        return float(profile.tm_C), profile
    candidate_tm = round(float(_candidate_thermodynamics(primer, relaxed, config.tm_method).tm), 2)
    return float(round(float(candidate_tm), 2)), None


def _v7_screening_score(candidate: PrimerCandidate, config) -> float:
    """Return the screening score _v7_measure_candidates gives a candidate."""
    advanced = analyze_primer_candidate_thermodynamics(
        candidate, conditions=config.thermodynamics,
    )
    advanced = _v7_apply_selected_tm_model(advanced, candidate, config)
    return float(_v5_candidate_quality(advanced))


def _v7_screening_upper_bound(sequence, orientation, start, end, config, relaxed, profile=None) -> float:
    """Return an upper bound of a window's screening score.

    It is the score of the same primer with no hairpin and no self-dimer:
    score_primer_candidate then awards the full structure score, and the
    structural penalty of _v5_candidate_quality is zero. Every other term
    is computed exactly as for the measured candidate, and each step is
    monotone, so the true score never exceeds the bound.
    """
    template = sequence[start:end]
    primer = template if orientation == "FORWARD" else reverse_complement(template)
    thermodynamics = _candidate_thermodynamics(primer, relaxed, config.tm_method)
    tm = float(thermodynamics.tm)
    gc = gc_percent(primer)
    no_structure = StructureResult(
        max_stem=0, three_prime_stem=0, best_position=None,
        severity="LOW", explanation="",
    )
    score = score_primer_candidate(
        length=len(primer), gc=gc, tm=tm, gc_clamp=count_gc_clamp(primer),
        homopolymer=longest_homopolymer(primer),
        repeats=longest_tandem_repeat(primer, relaxed.max_repeat_unit),
        hairpin=no_structure, self_dimer=no_structure, config=relaxed,
    )
    if profile is None:
        profile = calculate_primer_thermodynamics(primer, conditions=config.thermodynamics)
    thermodynamic_score = score_thermodynamic_profile(profile)
    method = config.tm_method.lower().strip()
    if method not in _V7_NEAREST_NEIGHBOR_METHODS:
        profile = replace(profile, method=method, tm_C=round(float(round(tm, 2)), 2))
        thermodynamic_score = score_thermodynamic_profile(
            profile,
            optimal_tm=config.primer.optimal_tm,
            min_tm=config.primer.min_tm,
            max_tm=config.primer.max_tm,
        )
    tm_error = abs(float(profile.tm_C) - float(config.primer.optimal_tm))
    gc_error = abs(float(round(gc, 2)) - float(config.primer.optimal_gc))
    return (
        float(thermodynamic_score) * 2.0
        + float(score)
        - 0.0
        - tm_error * 2.0
        - gc_error
    )


def _v7_tm_tier(tm: float, base_config: PrimerDesignConfig) -> int:
    if base_config.min_tm <= tm <= base_config.max_tm:
        return 1
    if (base_config.min_tm - 2.0) <= tm <= (base_config.max_tm + 2.0):
        return 2
    return 3


def _v7_measure_candidates(candidates, config):
    """Add the thermodynamic profile, tier and screening score, then order."""
    base_config = config.primer
    measured = []
    for candidate in candidates:
        advanced = analyze_primer_candidate_thermodynamics(
            candidate,
            conditions=config.thermodynamics,
        )

        advanced = _v7_apply_selected_tm_model(
            advanced,
            candidate,
            config,
        )

        advanced.screening_tier = _v7_tm_tier(
            float(advanced.thermodynamics.tm_C), base_config,
        )
        advanced.screening_score = _v5_candidate_quality(advanced)
        measured.append(advanced)

    measured.sort(
        key=lambda item: (
            item.screening_tier,
            -float(item.screening_score),
            abs(
                float(item.thermodynamics.tm_C)
                - float(base_config.optimal_tm)
            ),
            abs(
                float(item.candidate.gc_percent)
                - float(base_config.optimal_gc)
            ),
        )
    )
    return measured


def _v7_measure_all_candidates(sequence, target_start, target_end, orientation, config, relaxed, limit):
    """Measure every window (the behavior before the tier screen)."""
    candidates, evaluated = generate_primer_candidates(
        sequence,
        target_start,
        target_end,
        orientation,
        relaxed,
        config.tm_method,
    )
    measured = _v7_measure_candidates(candidates, config)
    return measured[:limit], evaluated


def _v7_primer_ui_design(sequence, target_start, target_end, pipeline):
    """Run all three Tm models through the professional workflow."""
    tm_method = _primer_ui_select_tm_model()

    if tm_method is None:
        return None

    conditions = _primer_ui_read_conditions()
    if conditions is None:
        return None

    print()
    print("RUNNING PRIMER DESIGN...")
    print()

    try:
        config = ProfessionalPrimerConfig(
            thermodynamics=conditions,
            tm_method=tm_method,
        )

        report = design_primers_professional(
            sequence,
            target_start=target_start,
            target_end=target_end,
            config=config,
        )

        print_professional_primer_report(report)
        _primer_ui_store_result(pipeline, report)
        return report

    except (TypeError, ValueError, ZeroDivisionError) as exc:
        style.error(f"{exc}")
        return None


# Bind V7 implementations.
generate_professional_primer_candidates = _v7_generate_professional_primer_candidates
_primer_ui_design = _v7_primer_ui_design


# ============================================================================
# V8 STRUCTURE CLASSIFICATION + HARD PAIR SAFETY GATES
# ============================================================================
#
# V7 exposed a real classification bug:
# the legacy structure functions used max_stem=4 as the severity threshold,
# which made a 20 bp or 24 bp stem appear as LOW in some reports.
#
# V8 fixes the structure severity itself and fixes the pair-ranking bug where
# the heterodimer length was read from a non-existent "stem_length" attribute.
# ============================================================================

def _v8_structure_severity(max_stem: int, three_prime_stem: int = 0) -> str:
    """Classify structural risk from the measured stem, not the old default."""
    max_stem = int(max_stem or 0)
    three_prime_stem = int(three_prime_stem or 0)

    if max_stem >= 10 or three_prime_stem >= 4:
        return "VERY HIGH"
    if max_stem >= 7 or three_prime_stem >= 3:
        return "HIGH"
    if max_stem >= 4 or three_prime_stem >= 2:
        return "MODERATE"
    return "LOW"


def analyze_hairpin(sequence: str, max_stem: int = 4) -> StructureResult:
    """Analyze hairpins using the existing detection algorithm and V8 severity."""
    sequence = normalize_dna_sequence(sequence)

    best = 0
    best_position = None
    length = len(sequence)
    max_test_stem = min(8, length // 2)
    # The reverse complement of sequence[r:r + stem] is the slice
    # [length - r - stem:length - r] of the whole reverse complement, so
    # it is computed once (P2-6: this loop dominated primer design).
    reverse = reverse_complement(sequence) if length else ""

    # Stems are tried from the longest down, and only the first match at
    # the longest matching stem sets best and best_position, so the
    # search ends there.
    #
    # A right arm starting at right_start pairs with the left arm when
    # left == reverse[k:k + stem] with k = length - right_start - stem.
    # right_start runs from left_start + stem + 3 to length - stem, so k
    # runs from 0 to length - left_start - 2 * stem - 3, and any such
    # occurrence ends before length - left_start - stem - 3. One bounded
    # str.find replaces the loop over right_start (P2-6).
    for stem in range(max_test_stem, 2, -1):
        for left_start in range(0, length - stem - 3):
            left = sequence[left_start:left_start + stem]
            if reverse.find(left, 0, length - left_start - stem - 3) != -1:
                best = stem
                best_position = left_start
                break
        if best:
            break

    result = StructureResult(
        max_stem=best,
        three_prime_stem=0,
        best_position=best_position,
        severity=_v8_structure_severity(best),
        explanation=(
            f"Longest detected intramolecular stem: {best} bp; "
            "3-nt minimum-loop heuristic."
        ),
    )
    return result


def analyze_self_dimer(sequence: str, max_stem: int = 4) -> StructureResult:
    """Analyze self-dimer formation with V8 severity thresholds."""
    sequence = normalize_dna_sequence(sequence)
    stem, position = _alignment_stem(sequence, sequence)
    terminal = _three_prime_complementarity(sequence)

    return StructureResult(
        max_stem=stem,
        three_prime_stem=terminal,
        best_position=position,
        severity=_v8_structure_severity(stem, terminal),
        explanation=(
            f"Longest self-complementary stem: {stem} bp; "
            f"estimated 3'-complementarity: {terminal} bp."
        ),
    )


def analyze_heterodimer(
    forward: str,
    reverse: str,
    max_stem: int = 4,
) -> StructureResult:
    """Analyze heterodimer formation with V8 severity thresholds."""
    forward = normalize_dna_sequence(forward)
    reverse = normalize_dna_sequence(reverse)

    stem, position = _alignment_stem(forward, reverse)
    terminal = max(
        _three_prime_complementarity(forward, reverse),
        _three_prime_complementarity(reverse, forward),
    )

    return StructureResult(
        max_stem=stem,
        three_prime_stem=terminal,
        best_position=position,
        severity=_v8_structure_severity(stem, terminal),
        explanation=(
            f"Longest heterodimer stem: {stem} bp; "
            f"estimated 3'-complementarity: {terminal} bp."
        ),
    )


# V8 bindings.


# ============================================================================
# V9 PRIMER RANKING / DEDUPLICATION
# ============================================================================
#
# V8 correctly classified structures, but the resulting score was too harsh:
# a pair with 6 bp moderate heterodimer/self-dimer and excellent Tm matching
# could fall near 40/100. V9 separates "quality score" from "risk category",
# removes duplicate-equivalent primer pairs, and uses explicit risk gates.
# ============================================================================

# ============================================================================
# V10 DIVERSE PAIR SELECTION
# ============================================================================
# Keep high-quality alternatives genuinely different instead of returning
# many one-base-shifted variants of the same primer pair.
# ============================================================================

# ============================================================================
# V11 FINAL PRIMER QUALITY ENGINE
# ============================================================================
# Final integration layer. The existing sequence generation and thermodynamic
# engine are intentionally preserved. V11 replaces only the final pair
# decision/ranking layer so that quality, risk, model choice, specificity and
# diversity are handled consistently.
# ============================================================================


def _v11_model_tm(pair, advanced, method):
    method = str(method).lower().replace("-", "_")
    if method in {"wallace"}:
        return float(pair.forward.tm), float(pair.reverse.tm)
    if method in {"salt", "salt_adjusted"}:
        return float(pair.forward.tm), float(pair.reverse.tm)
    return (
        float(advanced.forward_thermodynamics.tm_C),
        float(advanced.reverse_thermodynamics.tm_C),
    )


def _v11_tm_optimum(method):
    method = str(method).lower().replace("-", "_")
    if method == "wallace":
        return 55.0
    if method in {"salt", "salt_adjusted"}:
        return 58.0
    return 60.0


def _v11_piecewise_penalty(value, thresholds):
    """Return the penalty for a structural metric using (threshold, penalty)."""
    penalty = 0.0
    for threshold, amount in thresholds:
        if value >= threshold:
            penalty = amount
    return penalty


def _v11_pair_quality(pair, advanced, config):
    """Calculate a transparent 0-100 pair score plus explicit risk gates."""
    method = config.tm_method
    f_tm, r_tm = _v11_model_tm(pair, advanced, method)
    tm_opt = _v11_tm_optimum(method)
    tm_delta = abs(f_tm - r_tm)

    f_gc = float(pair.forward.gc_percent)
    r_gc = float(pair.reverse.gc_percent)
    f_len = int(pair.forward.length)
    r_len = int(pair.reverse.length)

    hetero = int(getattr(pair.heterodimer, "max_stem", 0) or 0)
    hetero_3p = int(getattr(pair.heterodimer, "three_prime_stem", 0) or 0)

    f_self = int(getattr(pair.forward.self_dimer, "max_stem", 0) or 0)
    r_self = int(getattr(pair.reverse.self_dimer, "max_stem", 0) or 0)
    f_self_3p = int(getattr(pair.forward.self_dimer, "three_prime_stem", 0) or 0)
    r_self_3p = int(getattr(pair.reverse.self_dimer, "three_prime_stem", 0) or 0)

    f_hairpin = int(getattr(pair.forward.hairpin, "max_stem", 0) or 0)
    r_hairpin = int(getattr(pair.reverse.hairpin, "max_stem", 0) or 0)

    # Start from a full score and subtract interpretable penalties.
    score = 100.0

    # Model-aware Tm quality. Wallace and salt-adjusted are ranked against
    # their own empirical scale; nearest-neighbor uses the configured target.
    score -= min(18.0, abs(f_tm - tm_opt) * 4.5)
    score -= min(18.0, abs(r_tm - tm_opt) * 4.5)
    score -= min(15.0, tm_delta * 7.5)

    # GC and length are soft preferences.
    for gc in (f_gc, r_gc):
        score -= min(6.0, abs(gc - 50.0) * 0.35)
    for length in (f_len, r_len):
        score -= min(5.0, abs(length - config.primer.optimal_length) * 0.9)

    # Structural penalties increase sharply once risk becomes meaningful.
    score -= _v11_piecewise_penalty(
        hetero,
        ((4, 2), (5, 4), (6, 7), (7, 13), (8, 22), (9, 32), (10, 45)),
    )
    score -= _v11_piecewise_penalty(
        hetero_3p,
        ((2, 2), (3, 10), (4, 22)),
    )

    for stem, terminal in (
        (f_self, f_self_3p),
        (r_self, r_self_3p),
    ):
        score -= _v11_piecewise_penalty(
            stem,
            ((4, 1), (5, 3), (6, 6), (7, 10), (8, 16), (9, 23), (10, 32)),
        )
        score -= _v11_piecewise_penalty(
            terminal,
            ((2, 2), (3, 9), (4, 20)),
        )

    for stem in (f_hairpin, r_hairpin):
        score -= _v11_piecewise_penalty(
            stem,
            ((4, 1), (5, 3), (6, 6), (7, 11), (8, 17), (9, 24), (10, 32)),
        )

    # Specificity hits are only meaningful when a reference set was supplied.
    hits = list(getattr(pair, "specificity_hits", []) or [])
    exact_hits = sum(1 for hit in hits if getattr(hit, "exact", False))
    terminal_hits = sum(
        1 for hit in hits if int(getattr(hit, "three_prime_mismatches", 99)) <= 1
    )
    if exact_hits:
        score -= min(25.0, exact_hits * 12.5)
    elif terminal_hits:
        score -= min(15.0, terminal_hits * 5.0)
    elif hits:
        score -= min(8.0, len(hits) * 1.5)

    # Explicit risk gates. These prevent a tiny Tm advantage from rescuing a
    # structurally dangerous pair.
    critical = bool(
        hetero >= 10
        or hetero_3p >= 4
        or f_self >= 10
        or r_self >= 10
        or f_self_3p >= 4
        or r_self_3p >= 4
        or f_hairpin >= 10
        or r_hairpin >= 10
        or exact_hits > 0
    )
    high_risk = bool(
        critical
        or hetero >= 8
        or hetero_3p >= 3
        or f_self >= 8
        or r_self >= 8
        or f_self_3p >= 3
        or r_self_3p >= 3
        or f_hairpin >= 7
        or r_hairpin >= 7
        or terminal_hits >= 2
    )

    reasons = []
    if hetero >= 7 or hetero_3p >= 3:
        reasons.append(
            f"heterodimer risk: {hetero} bp stem, {hetero_3p} bp 3'-complementarity"
        )
    if f_self >= 7 or f_self_3p >= 3:
        reasons.append(
            f"forward self-dimer risk: {f_self} bp stem, {f_self_3p} bp 3'-complementarity"
        )
    if r_self >= 7 or r_self_3p >= 3:
        reasons.append(
            f"reverse self-dimer risk: {r_self} bp stem, {r_self_3p} bp 3'-complementarity"
        )
    if f_hairpin >= 7:
        reasons.append(f"forward hairpin risk: {f_hairpin} bp stem")
    if r_hairpin >= 7:
        reasons.append(f"reverse hairpin risk: {r_hairpin} bp stem")
    if exact_hits:
        reasons.append(f"{exact_hits} exact specificity hit(s) in supplied references")
    elif terminal_hits:
        reasons.append(f"{terminal_hits} strong terminal specificity hit(s) in supplied references")
    elif hits:
        reasons.append(f"{len(hits)} near-match specificity hit(s) in supplied references")

    score = round(max(0.0, min(100.0, score)), 2)

    if critical:
        category = "POOR"
    elif high_risk:
        category = "MARGINAL"
    elif score >= 85.0:
        category = "EXCELLENT"
    elif score >= 72.0:
        category = "GOOD"
    elif score >= 60.0:
        category = "ACCEPTABLE"
    else:
        category = "MARGINAL"

    return score, high_risk, critical, reasons, category


def _v11_sequence_similarity(a, b):
    """Simple normalized identity used only for diversity selection."""
    a = str(a).upper()
    b = str(b).upper()
    if not a or not b:
        return 0.0
    overlap = min(len(a), len(b))
    same = sum(x == y for x, y in zip(a[:overlap], b[:overlap]))
    return same / max(len(a), len(b))


def _v11_is_diverse(candidate, selected):
    """Reject alternatives that are effectively the same design."""
    pair = candidate.pair
    for old in selected:
        old_pair = old.pair
        f_identity = _v11_sequence_similarity(pair.forward.sequence, old_pair.forward.sequence)
        r_identity = _v11_sequence_similarity(pair.reverse.sequence, old_pair.reverse.sequence)
        f_shift = abs(int(pair.forward.start) - int(old_pair.forward.start))
        r_shift = abs(int(pair.reverse.start) - int(old_pair.reverse.start))

        # Exact/reused primers are never useful as separate top alternatives.
        if pair.forward.sequence.upper() == old_pair.forward.sequence.upper():
            return False
        if pair.reverse.sequence.upper() == old_pair.reverse.sequence.upper():
            return False

        # Near-identical primers at nearly identical loci are duplicates in
        # practical terms even when the literal strings differ.
        if f_identity >= 0.90 and r_identity >= 0.90 and f_shift < 5 and r_shift < 5:
            return False
    return True


def _v11_select_pairs(pairs, top_n):
    """Select safe, high-quality and genuinely diverse alternatives."""
    if not pairs:
        return []

    safe = [x for x in pairs if not getattr(x, "critical_risk", False) and not getattr(x, "high_risk", False)]
    marginal = [x for x in pairs if getattr(x, "high_risk", False) and not getattr(x, "critical_risk", False)]
    poor = [x for x in pairs if getattr(x, "critical_risk", False)]

    selected = []

    # Prefer safe pairs. If the sequence genuinely has fewer than top_n safe
    # options, progressively relax to marginal and finally poor candidates.
    for pool in (safe, marginal, poor):
        for item in pool:
            if len(selected) >= top_n:
                break
            if _v11_is_diverse(item, selected):
                selected.append(item)
        if len(selected) >= top_n:
            break

    # Last-resort fill: never silently return fewer candidates merely because
    # the input sequence is highly repetitive. Risk labels remain intact.
    if len(selected) < top_n:
        used = {
            (x.pair.forward.sequence.upper(), x.pair.reverse.sequence.upper())
            for x in selected
        }
        for item in pairs:
            key = (item.pair.forward.sequence.upper(), item.pair.reverse.sequence.upper())
            if key in used:
                continue
            selected.append(item)
            used.add(key)
            if len(selected) >= top_n:
                break

    return selected


def _v11_rank_professional_primer_pairs(
    sequence,
    forward_candidates,
    reverse_candidates,
    *,
    config=None,
    reference_sequences=None,
):
    """Final pair evaluator with risk-first quality ranking and diversity."""
    config = config or ProfessionalPrimerConfig()
    config.validate()

    pairs = []
    evaluated = 0
    seen = set()

    # Candidate pools are already limited by max_pair_candidates upstream.
    for forward_item in forward_candidates:
        for reverse_item in reverse_candidates:
            forward = forward_item.candidate
            reverse = reverse_item.candidate

            if forward.end > reverse.start:
                continue
            evaluated += 1

            pair = build_primer_pair(
                sequence,
                forward,
                reverse,
                config.primer,
                reference_sequences,
            )
            if pair is None:
                continue

            key = (forward.sequence.upper(), reverse.sequence.upper())
            if key in seen:
                continue
            seen.add(key)

            advanced = analyze_primer_pair_thermodynamics(
                pair,
                conditions=config.thermodynamics,
                max_tm_delta=config.primer.max_tm_delta,
            )

            score, high_risk, critical, reasons, category = _v11_pair_quality(
                pair,
                advanced,
                config,
            )

            pair.score = score
            pair.category = category
            pair.warnings = list(getattr(pair, "warnings", []))
            pair.warnings.extend(f"V11 diagnostic: {reason}" for reason in reasons)

            advanced.screening_score = score
            advanced.high_risk = high_risk
            advanced.critical_risk = critical
            advanced.structural_warnings = reasons
            advanced.model_tm_delta = abs(_v11_model_tm(pair, advanced, config.tm_method)[0] - _v11_model_tm(pair, advanced, config.tm_method)[1])
            advanced.specificity_screened = bool(reference_sequences)

            pairs.append(advanced)

    # Safe quality first. Structural risk is a hard ordering dimension; Tm
    # matching is only used to break ties between similarly safe candidates.
    pairs.sort(
        key=lambda item: (
            bool(getattr(item, "critical_risk", False)),
            bool(getattr(item, "high_risk", False)),
            -float(getattr(item, "screening_score", 0.0)),
            float(getattr(item, "model_tm_delta", 999.0)),
            int(getattr(item.pair.heterodimer, "max_stem", 0)),
            abs(int(item.pair.amplicon_length) - int(config.primer.optimal_amplicon_length)),
        )
    )

    selected = _v11_select_pairs(pairs, int(config.top_pairs))
    return selected, evaluated


rank_professional_primer_pairs = _v11_rank_professional_primer_pairs


# ============================================================================
# V11 FINAL REPORTING FIXES
# ============================================================================


def summarize_professional_primer_pair(advanced_pair: AdvancedPrimerPair, rank: int = 1) -> str:
    """Final pair report with accurate specificity state and diagnostics."""
    pair = advanced_pair.pair
    forward = pair.forward
    reverse = pair.reverse
    f_thermo = advanced_pair.forward_thermodynamics
    r_thermo = advanced_pair.reverse_thermodynamics

    screened = bool(getattr(advanced_pair, "specificity_screened", False))
    if not screened:
        specificity = "NOT SCREENED"
    elif pair.specificity_hits:
        specificity = f"{len(pair.specificity_hits)} HIT(S)"
    else:
        specificity = "SCREENED — NO SIGNIFICANT HITS"

    lines = [
        f"PRIMER PAIR #{rank}",
        style.field_line('Score', f"{pair.score:.2f}/100"),
        style.field_line('Category', f"{pair.category}"),
        "",
        "FORWARD",
        style.field_line('Sequence', f"{forward.sequence}"),
        style.field_line('Position', f"{forward.start + 1}-{forward.end} (1-based)"),
        style.field_line('Length', f"{forward.length} bp"),
        style.field_line('GC', f"{forward.gc_percent:.2f}%"),
        style.field_line('Tm', f"{f_thermo.tm_C:.2f} °C"),
        f"ΔH:                    {f_thermo.delta_h_kcal_mol:.2f} kcal/mol",
        f"ΔS:                    {f_thermo.delta_s_cal_mol_k:.2f} cal/(mol·K)",
        f"ΔG at 25 °C:           {f_thermo.delta_g_kcal_mol:.2f} kcal/mol",
        style.field_line('Effective salt', f"{f_thermo.effective_salt_mM:.2f} mM"),
        f"3' GC clamp:           {forward.gc_clamp} bp",
        style.field_line('Hairpin stem', f"{forward.hairpin.max_stem} bp ({forward.hairpin.severity})"),
        style.field_line('Self-dimer stem', f"{forward.self_dimer.max_stem} bp ({forward.self_dimer.severity})"),
        "",
        "REVERSE",
        style.field_line('Sequence', f"{reverse.sequence}"),
        style.field_line('Position', f"{reverse.start + 1}-{reverse.end} (1-based)"),
        style.field_line('Length', f"{reverse.length} bp"),
        style.field_line('GC', f"{reverse.gc_percent:.2f}%"),
        style.field_line('Tm', f"{r_thermo.tm_C:.2f} °C"),
        f"ΔH:                    {r_thermo.delta_h_kcal_mol:.2f} kcal/mol",
        f"ΔS:                    {r_thermo.delta_s_cal_mol_k:.2f} cal/(mol·K)",
        f"ΔG at 25 °C:           {r_thermo.delta_g_kcal_mol:.2f} kcal/mol",
        style.field_line('Effective salt', f"{r_thermo.effective_salt_mM:.2f} mM"),
        f"3' GC clamp:           {reverse.gc_clamp} bp",
        style.field_line('Hairpin stem', f"{reverse.hairpin.max_stem} bp ({reverse.hairpin.severity})"),
        style.field_line('Self-dimer stem', f"{reverse.self_dimer.max_stem} bp ({reverse.self_dimer.severity})"),
        "",
        "PAIR",
        style.field_line('Amplicon', f"{pair.amplicon_start}-{pair.amplicon_end}"),
        style.field_line('Amplicon length', f"{pair.amplicon_length} bp"),
        f"ΔTm:                   {advanced_pair.thermodynamic_tm_delta:.2f} °C",
        style.field_line('Thermodynamic score', f"{advanced_pair.thermodynamic_pair_score:.2f}/15"),
        style.field_line('Heterodimer stem', f"{pair.heterodimer.max_stem} bp ({pair.heterodimer.severity})"),
        style.field_line('Specificity screening', f"{specificity}"),
    ]
    return "\n".join(lines)


def design_primers_professional(
    sequence: str,
    target_start: int = 0,
    target_end: Optional[int] = None,
    *,
    config: ProfessionalPrimerConfig = ProfessionalPrimerConfig(),
    reference_sequences: Optional[dict[str, str] | Sequence[tuple[str, str]]] = None,
) -> ProfessionalPrimerDesignReport:
    """Final integrated professional workflow."""
    config.validate()
    sequence = normalize_dna_sequence(sequence)
    if target_end is None:
        target_end = len(sequence)
    validate_target_region(sequence, target_start, target_end)

    forward_candidates, forward_evaluated = generate_professional_primer_candidates(
        sequence, target_start, target_end, "FORWARD", config=config
    )
    reverse_candidates, reverse_evaluated = generate_professional_primer_candidates(
        sequence, target_start, target_end, "REVERSE", config=config
    )
    return _professional_primer_report(
        sequence, target_start, target_end, config, reference_sequences,
        forward_candidates, forward_evaluated, reverse_candidates, reverse_evaluated,
    )


def _professional_primer_report(
    sequence: str,
    target_start: int,
    target_end: int,
    config: ProfessionalPrimerConfig,
    reference_sequences,
    forward_candidates,
    forward_evaluated: int,
    reverse_candidates,
    reverse_evaluated: int,
) -> ProfessionalPrimerDesignReport:
    """Pair, rank and report forward and reverse candidates."""
    ranked_pairs, pairs_evaluated = rank_professional_primer_pairs(
        sequence,
        forward_candidates,
        reverse_candidates,
        config=config,
        reference_sequences=reference_sequences,
    )

    failure_reasons = []
    if not forward_candidates:
        failure_reasons.append("No thermodynamically acceptable forward primer candidates found.")
    if not reverse_candidates:
        failure_reasons.append("No thermodynamically acceptable reverse primer candidates found.")
    if forward_candidates and reverse_candidates and not ranked_pairs:
        failure_reasons.append("No primer pair satisfies the configured amplicon and compatibility constraints.")

    warnings = []
    if not reference_sequences:
        warnings.append("Specificity was not genome-wide screened; no reference database was supplied.")
    else:
        total_hits = sum(len(x.pair.specificity_hits) for x in ranked_pairs)
        if total_hits:
            warnings.append(f"Specificity screening found {total_hits} hit(s) across the reported primer pairs in the supplied references.")
        else:
            warnings.append("Specificity screening was performed against the supplied reference sequences; no significant hits were detected in the reported pairs.")

    method = config.tm_method.lower().replace("-", "_")
    if method == "nearest_neighbor":
        warnings.append("Nearest-neighbor Tm is a sequence-level thermodynamic screening model; validate final primers under the actual PCR conditions.")
        warnings.append("Mg2+ and dNTP effects use an effective ionic-strength approximation rather than a complete buffer-speciation model.")
    elif method == "wallace":
        warnings.append("Wallace Tm is an empirical screening estimate; use sequence-specific thermodynamics for final experimental selection.")
    else:
        warnings.append("Salt-adjusted Tm is an empirical screening estimate; Mg2+, dNTP and sequence-specific thermodynamics are not fully represented by this model.")

    if ranked_pairs and all(x.pair.category in {"MARGINAL", "POOR"} for x in ranked_pairs):
        warnings.append("Every reported primer pair is rated MARGINAL or POOR under the current structural and thermodynamic criteria; these pairs are best-effort alternatives and should not be treated as validated designs.")

    return ProfessionalPrimerDesignReport(
        sequence_length=len(sequence),
        target_start=target_start,
        target_end=target_end,
        target_length=target_end - target_start,
        tm_method=config.tm_method,
        thermodynamic_conditions=config.thermodynamics,
        forward_candidates_evaluated=forward_evaluated,
        reverse_candidates_evaluated=reverse_evaluated,
        forward_candidates_accepted=len(forward_candidates),
        reverse_candidates_accepted=len(reverse_candidates),
        pairs_evaluated=pairs_evaluated,
        pairs_accepted=len(ranked_pairs),
        ranked_pairs=ranked_pairs,
        failure_reasons=failure_reasons,
        warnings=warnings,
    )


# ============================================================================
# PRIMERS AROUND A RESTRICTION SITE (Master Prompt 10.2)
# ============================================================================
# Workflow: restriction analysis -> choice of one confirmed site -> primer
# design on both sides of it -> thermodynamic evaluation and ranking ->
# specificity -> report. The forward primers are placed upstream of the
# site and the reverse primers downstream, so every amplicon contains the
# site and both cuts; each pair is reported with the predicted digest of
# its amplicon by the same enzyme (a diagnostic digest, for example to
# tell an intact site from a mutated one).
# ============================================================================


@dataclass
class RestrictionSitePrimerDesign:
    """Primer pairs that flank one confirmed restriction site."""

    enzyme: str
    site_start: int
    site_end: int
    protected_start: int
    protected_end: int
    forward_region: tuple[int, int]
    reverse_region: tuple[int, int]
    flank: int
    report: ProfessionalPrimerDesignReport
    amplicon_digests: list[dict]
    warnings: list[str]


def restriction_site_choices(sequence: str, enzyme_name: str) -> list[dict]:
    """Return the confirmed sites of one enzyme with their cut boundaries.

    Positions are 0-based: ``site_start``/``site_end`` delimit the
    recognition sequence (half-open), ``top_cut``/``bottom_cut`` are the cut
    boundaries on the two strands, and ``protected_start``/``protected_end``
    span the site and both cuts (for a Type IIS enzyme the cuts lie outside
    the site). Only CONFIRMED sites are listed: a site that matches only
    through an ambiguity code is not evidence of a cut.
    """
    from analyses import restriction as R

    sequence = normalize_dna_sequence(sequence, allow_iupac=True)
    type_iis = R.get_type_iis_enzyme(enzyme_name)
    if type_iis is None and R.get_restriction_enzyme(enzyme_name) is None:
        raise ValueError(f"Unknown restriction enzyme: {enzyme_name}.")
    recognition = (type_iis or R.get_restriction_enzyme(enzyme_name) or {})["recognition_sequence"]
    choices = []
    for cut in R.get_cut_positions(sequence, [enzyme_name], topology="linear"):
        if cut.get("status") != "CONFIRMED" or cut.get("wraps_origin"):
            continue
        site_start = int(cut["site_position_0based"])
        site_end = site_start + len(recognition)
        top_cut = int(cut["top_cut_0based"])
        bottom_cut = int(cut["bottom_cut_0based"])
        choices.append({
            "enzyme": enzyme_name,
            "site_start": site_start,
            "site_end": site_end,
            "site_sequence": sequence[site_start:site_end],
            "strand": cut.get("strand", "+"),
            "top_cut": top_cut,
            "bottom_cut": bottom_cut,
            "protected_start": min(site_start, top_cut, bottom_cut),
            "protected_end": max(site_end, top_cut, bottom_cut),
        })
    return choices


def design_primers_for_restriction_site(
    sequence: str,
    enzyme_name: str,
    site_start: int,
    *,
    flank: int = 300,
    config: ProfessionalPrimerConfig = ProfessionalPrimerConfig(),
    reference_sequences: Optional[dict[str, str] | Sequence[tuple[str, str]]] = None,
) -> RestrictionSitePrimerDesign:
    """Design primer pairs whose amplicon contains one restriction site.

    ``site_start`` is the 0-based start of a CONFIRMED recognition site of
    ``enzyme_name`` on the linear template. Forward primers are searched in
    the ``flank`` bases upstream of the site and its cuts, reverse primers
    in the ``flank`` bases downstream.
    """
    from analyses import restriction as R

    config.validate()
    sequence = normalize_dna_sequence(sequence)
    if int(flank) < config.primer.min_length:
        raise ValueError(
            f"The flank must be at least the minimum primer length ({config.primer.min_length} nt)."
        )
    flank = int(flank)
    choice = next(
        (item for item in restriction_site_choices(sequence, enzyme_name)
         if item["site_start"] == int(site_start)),
        None,
    )
    if choice is None:
        raise ValueError(
            f"No confirmed {enzyme_name} site starts at position {int(site_start) + 1} (1-based)."
        )
    protected_start = choice["protected_start"]
    protected_end = choice["protected_end"]
    forward_region = (max(0, protected_start - flank), protected_start)
    reverse_region = (protected_end, min(len(sequence), protected_end + flank))
    for name, (start, end) in (("upstream", forward_region), ("downstream", reverse_region)):
        if end - start < config.primer.min_length:
            raise ValueError(
                f"Only {end - start} bases lie {name} of the site; a primer needs "
                f"at least {config.primer.min_length}."
            )

    forward_candidates, forward_evaluated = generate_professional_primer_candidates(
        sequence, forward_region[0], forward_region[1], "FORWARD", config=config
    )
    reverse_candidates, reverse_evaluated = generate_professional_primer_candidates(
        sequence, reverse_region[0], reverse_region[1], "REVERSE", config=config
    )
    report = _professional_primer_report(
        sequence, forward_region[0], reverse_region[1], config, reference_sequences,
        forward_candidates, forward_evaluated, reverse_candidates, reverse_evaluated,
    )

    digests = []
    for rank, advanced in enumerate(report.ranked_pairs, 1):
        forward = advanced.pair.forward
        reverse = advanced.pair.reverse
        amplicon = sequence[forward.start:reverse.end]
        digest = R.perform_digest(amplicon, [enzyme_name], "linear")
        digests.append({
            "rank": rank,
            "amplicon_start": forward.start,
            "amplicon_end": reverse.end,
            "amplicon_length": len(amplicon),
            "site_offset": choice["site_start"] - forward.start,
            "sites_in_amplicon": int(digest["cut_count"]),
            "fragment_lengths": [int(fragment["length"]) for fragment in digest["fragments"]],
        })

    warnings = [
        "The digest is predicted from the sequence; incomplete digestion, "
        "methylation and star activity are not modeled.",
    ]
    extra = sorted({item["sites_in_amplicon"] for item in digests if item["sites_in_amplicon"] > 1})
    if extra:
        warnings.append(
            f"Some amplicons contain more than one {enzyme_name} site "
            f"({', '.join(str(n) for n in extra)}); their digest gives more than two fragments."
        )
    if not report.ranked_pairs:
        warnings.append("No primer pair flanks the site under the current constraints; "
                        "a larger flank or a wider Tm range may help.")

    return RestrictionSitePrimerDesign(
        enzyme=enzyme_name,
        site_start=choice["site_start"],
        site_end=choice["site_end"],
        protected_start=protected_start,
        protected_end=protected_end,
        forward_region=forward_region,
        reverse_region=reverse_region,
        flank=flank,
        report=report,
        amplicon_digests=digests,
        warnings=warnings,
    )


def _primer_ui_restriction_site_design(sequence, pipeline):
    """Design primers around one confirmed restriction site (workflow 10.2)."""
    from data.restriction_enzymes import RESTRICTION_ENZYMES, TYPE_IIS_ENZYMES

    style.section("RESTRICTION SITE")
    known = {name.upper(): name for name in list(RESTRICTION_ENZYMES) + list(TYPE_IIS_ENZYMES)}
    answer = style.ask("Enzyme (for example EcoRI or BsaI)").strip()
    enzyme = known.get(answer.upper())
    if enzyme is None:
        style.error(f"Unknown restriction enzyme: {answer or '(none)'}.")
        return None

    choices = restriction_site_choices(sequence, enzyme)
    if not choices:
        style.error(f"{enzyme} has no confirmed site in this sequence.")
        return None

    while True:
        style.section(f"CONFIRMED {enzyme.upper()} SITES")
        style.menu([(None, [
            (str(number), f"Position {item['site_start'] + 1}-{item['site_end']}  "
                          f"{item['site_sequence']}  (cut after {item['top_cut']})")
            for number, item in enumerate(choices, 1)
        ])], zero=style.CANCEL_LABEL)
        selection = style.ask_choice()
        if selection == "0":
            return None
        if selection.isdigit() and 1 <= int(selection) <= len(choices):
            chosen = choices[int(selection) - 1]
            break
        style.invalid_selection()

    flank_input = style.ask("Bases searched on each side [300]")
    try:
        flank = 300 if not flank_input else int(flank_input)
    except ValueError:
        style.error("The flank must be a whole number of bases.")
        return None

    tm_method = _primer_ui_select_tm_model()
    if tm_method is None:
        return None
    conditions = _primer_ui_read_conditions()
    if conditions is None:
        return None

    try:
        config = ProfessionalPrimerConfig(thermodynamics=conditions, tm_method=tm_method)
        design = design_primers_for_restriction_site(
            sequence, enzyme, chosen["site_start"], flank=flank, config=config,
        )
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        style.error(f"{exc}")
        return None

    print_professional_primer_report(design.report)
    print_restriction_site_design(design)
    _primer_ui_store_result(pipeline, design)
    return design


def print_restriction_site_design(design: RestrictionSitePrimerDesign) -> None:
    """Print the site, the search regions and each amplicon's digest."""
    style.section("RESTRICTION SITE AND DIAGNOSTIC DIGEST")
    style.field("Enzyme", design.enzyme)
    style.field("Site", f"{design.site_start + 1}-{design.site_end} (1-based)")
    style.field("Forward primers in", f"{design.forward_region[0] + 1}-{design.forward_region[1]}")
    style.field("Reverse primers in", f"{design.reverse_region[0] + 1}-{design.reverse_region[1]}")
    print()
    widths = [6, 18, 10, 7]
    print(style.columns(["Rank", "Amplicon", "Length", "Sites", "Fragments (bp)"], widths))
    for item in design.amplicon_digests:
        print(style.columns([
            item["rank"],
            f"{item['amplicon_start'] + 1}-{item['amplicon_end']}",
            f"{item['amplicon_length']} bp",
            item["sites_in_amplicon"],
            " + ".join(str(length) for length in item["fragment_lengths"]),
        ], widths))
    print()
    style.notes(design.warnings, prefix=style.WARNING_PREFIX)
