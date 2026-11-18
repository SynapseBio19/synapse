"""The legacy primer-design engine: candidate scoring and generation, pair
building and ranking, ``design_primers`` and ``analyze_primer``."""

from __future__ import annotations

from typing import cast, Optional, Sequence

from analyses._primer.sequence import (
    count_gc_clamp,
    gc_percent,
    longest_homopolymer,
    longest_tandem_repeat,
    normalize_dna_sequence,
    reverse_complement,
    validate_target_region,
)
from analyses._primer.models import (
    PrimerCandidate,
    PrimerDesignConfig,
    PrimerDesignReport,
    PrimerPair,
    StructureResult,
    ThermodynamicResult,
)
from analyses._primer.thermodynamics import (
    calculate_tm,
    calculate_tm_nearest_neighbor,
    get_thermodynamic_conditions,
    NEAREST_NEIGHBOR_TM_METHODS,
    NearestNeighborTmResult,
)
from analyses._primer.structure import analyze_hairpin, analyze_heterodimer, analyze_self_dimer
from analyses._primer.specificity import screen_primer_specificity


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


# The terminal dialogs of this module are in ui/primer_menu.py (P2-2).


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
