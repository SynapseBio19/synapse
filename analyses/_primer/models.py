"""Configuration and result dataclasses of the legacy primer-design engine."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


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
