"""The professional candidate pool: measurement of every window, the
screening score and the Tm-tier screen with its upper bound."""

from __future__ import annotations

from dataclasses import replace

from analyses._primer.sequence import (
    count_gc_clamp,
    gc_percent,
    longest_homopolymer,
    longest_tandem_repeat,
    normalize_dna_sequence,
    reverse_complement,
    validate_target_region,
)
from analyses._primer.models import PrimerCandidate, PrimerDesignConfig, StructureResult
from analyses._primer.thermodynamics import (
    calculate_primer_thermodynamics,
    NEAREST_NEIGHBOR_TM_METHODS,
    PrimerThermodynamicProfile,
)
from analyses._primer.design import (
    _candidate_thermodynamics,
    _candidate_windows,
    generate_primer_candidates,
    _generation_sort_key,
    _make_candidate,
    score_primer_candidate,
)
from analyses._primer.professional import (
    AdvancedPrimerCandidate,
    analyze_primer_candidate_thermodynamics,
    ProfessionalPrimerConfig,
    score_thermodynamic_profile,
)


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


def _v5_candidate_quality(item, config) -> float:
    """Score a candidate while strongly penalizing structural risk.

    ``config`` is the ProfessionalPrimerConfig of the design. Before Phase 7.1
    it was read from a module global that the candidate generator set on
    every call (``_V5_ACTIVE_CONFIG``); both callers already held the same
    object as their ``config`` argument.
    """
    profile = item.thermodynamics
    candidate = item.candidate

    # The stems are measured on the candidate. Before 2026-11-12 they were
    # read from the thermodynamic profile, which has no such fields, so the
    # structural penalty was always zero.
    hairpin = int(candidate.hairpin.max_stem or 0)
    self_dimer = int(candidate.self_dimer.max_stem or 0)

    # Prefer primers close to the configured optimum, but structure dominates.
    tm_error = abs(
        float(profile.tm_C) - float(config.primer.optimal_tm)
    )
    gc_error = abs(
        float(candidate.gc_percent) - float(config.primer.optimal_gc)
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
    config = config or ProfessionalPrimerConfig()
    config.validate()

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
    return float(_v5_candidate_quality(advanced, config))


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
        advanced.screening_score = _v5_candidate_quality(advanced, config)
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


# Bind V7 implementations.
generate_professional_primer_candidates = _v7_generate_professional_primer_candidates
