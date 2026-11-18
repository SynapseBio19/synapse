"""Professional pair scoring, diversity selection and ranking."""

from __future__ import annotations

from core.textformat import field_line

from analyses._primer.design import build_primer_pair
from analyses._primer.professional import (
    AdvancedPrimerPair,
    analyze_primer_pair_thermodynamics,
    ProfessionalPrimerConfig,
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
        field_line('Score', f"{pair.score:.2f}/100"),
        field_line('Category', f"{pair.category}"),
        "",
        "FORWARD",
        field_line('Sequence', f"{forward.sequence}"),
        field_line('Position', f"{forward.start + 1}-{forward.end} (1-based)"),
        field_line('Length', f"{forward.length} bp"),
        field_line('GC', f"{forward.gc_percent:.2f}%"),
        field_line('Tm', f"{f_thermo.tm_C:.2f} °C"),
        f"ΔH:                    {f_thermo.delta_h_kcal_mol:.2f} kcal/mol",
        f"ΔS:                    {f_thermo.delta_s_cal_mol_k:.2f} cal/(mol·K)",
        f"ΔG at 25 °C:           {f_thermo.delta_g_kcal_mol:.2f} kcal/mol",
        field_line('Effective salt', f"{f_thermo.effective_salt_mM:.2f} mM"),
        f"3' GC clamp:           {forward.gc_clamp} bp",
        field_line('Hairpin stem', f"{forward.hairpin.max_stem} bp ({forward.hairpin.severity})"),
        field_line('Self-dimer stem', f"{forward.self_dimer.max_stem} bp ({forward.self_dimer.severity})"),
        "",
        "REVERSE",
        field_line('Sequence', f"{reverse.sequence}"),
        field_line('Position', f"{reverse.start + 1}-{reverse.end} (1-based)"),
        field_line('Length', f"{reverse.length} bp"),
        field_line('GC', f"{reverse.gc_percent:.2f}%"),
        field_line('Tm', f"{r_thermo.tm_C:.2f} °C"),
        f"ΔH:                    {r_thermo.delta_h_kcal_mol:.2f} kcal/mol",
        f"ΔS:                    {r_thermo.delta_s_cal_mol_k:.2f} cal/(mol·K)",
        f"ΔG at 25 °C:           {r_thermo.delta_g_kcal_mol:.2f} kcal/mol",
        field_line('Effective salt', f"{r_thermo.effective_salt_mM:.2f} mM"),
        f"3' GC clamp:           {reverse.gc_clamp} bp",
        field_line('Hairpin stem', f"{reverse.hairpin.max_stem} bp ({reverse.hairpin.severity})"),
        field_line('Self-dimer stem', f"{reverse.self_dimer.max_stem} bp ({reverse.self_dimer.severity})"),
        "",
        "PAIR",
        field_line('Amplicon', f"{pair.amplicon_start}-{pair.amplicon_end}"),
        field_line('Amplicon length', f"{pair.amplicon_length} bp"),
        f"ΔTm:                   {advanced_pair.thermodynamic_tm_delta:.2f} °C",
        field_line('Thermodynamic score', f"{advanced_pair.thermodynamic_pair_score:.2f}/15"),
        field_line('Heterodimer stem', f"{pair.heterodimer.max_stem} bp ({pair.heterodimer.severity})"),
        field_line('Specificity screening', f"{specificity}"),
    ]
    return "\n".join(lines)
