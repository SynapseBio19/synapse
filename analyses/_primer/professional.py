"""The professional primer model: advanced candidates and pairs, their
thermodynamic scores, the configuration and report, single-primer
analysis and the failure diagnosis."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

from analyses._primer.sequence import (
    count_gc_clamp,
    gc_percent,
    longest_homopolymer,
    longest_tandem_repeat,
    normalize_dna_sequence,
)
from analyses._primer.models import PrimerCandidate, PrimerDesignConfig, PrimerPair
from analyses._primer.thermodynamics import (
    calculate_effective_ionic_strength,
    calculate_primer_thermodynamics,
    calculate_terminal_stability,
    calculate_tm_salt_adjusted,
    calculate_tm_wallace,
    PrimerThermodynamicProfile,
    ThermodynamicConfig,
)
from analyses._primer.structure import analyze_hairpin, analyze_self_dimer


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
    # Set by the candidate screen. They were attributes outside the
    # dataclass fields before 2026-11-12 and so never reached an export
    # (D23).
    screening_tier: Optional[int] = None
    screening_score: Optional[float] = None


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
    # The pair assessment of the ranking step (D23: attributes outside the
    # dataclass fields until 2026-11-12, so exports lost the risk flags and
    # the structural warnings). The defaults are the values the ranking
    # assumed for a pair it had not assessed.
    screening_score: float = 0.0
    high_risk: bool = False
    critical_risk: bool = False
    structural_warnings: list[str] = field(default_factory=list)
    model_tm_delta: float = 999.0
    specificity_screened: bool = False


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
