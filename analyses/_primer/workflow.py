"""Professional primer design for a target region and for a restriction site."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

from analyses._primer.sequence import normalize_dna_sequence, validate_target_region
from analyses._primer.professional import ProfessionalPrimerConfig, ProfessionalPrimerDesignReport
from analyses._primer.screening import generate_professional_primer_candidates
from analyses._primer.ranking import rank_professional_primer_pairs


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
