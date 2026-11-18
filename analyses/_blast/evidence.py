"""Evidence score, evidence profile, taxonomy grouping and alignment text."""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence

from analyses._blast.config import DEFAULT_EVIDENCE_EVALUE_THRESHOLD
from analyses._blast.models import BlastAnalysisError, BlastHit, BlastRun, HitEvidence, rank_hits


# ============================================================================
# HIT ANALYSIS
# ============================================================================

def calculate_evidence_score(
    hit: BlastHit,
    *,
    significance_evalue: float = DEFAULT_EVIDENCE_EVALUE_THRESHOLD,
) -> HitEvidence:
    """
    Transparent heuristic evidence profile.

    This is NOT a BLAST score and must not be interpreted as experimental
    confidence. It summarizes returned alignment evidence for Synapse display.

    Evidence bands are assigned only when the best E-value is <=
    ``significance_evalue``; otherwise the band is "NOT_SIGNIFICANT" whatever
    the identity and coverage (audit P0-1). Coverage is the NCBI-style query
    coverage of the union of all HSPs.
    """
    hsp = hit.best_hsp
    if hsp is None:
        return HitEvidence(
            hit_id=hit.hit_id,
            accession=hit.accession,
            identity_percent=0.0,
            query_coverage_percent=0.0,
            evalue=None,
            bit_score=0.0,
            alignment_length=0,
            evidence_score=0.0,
            evidence_band="NO_ALIGNMENT",
            notes=["No HSP was returned for this hit."],
        )

    coverage_value = hit.best_coverage_percent
    if coverage_value != coverage_value:  # NaN: query length unknown
        coverage_value = 0.0

    identity_component = min(100.0, max(0.0, hsp.identity_percent))
    coverage_component = min(100.0, max(0.0, coverage_value))

    if hsp.evalue <= 0:
        evalue_component = 100.0
    else:
        evalue_component = min(100.0, max(0.0, -math.log10(max(hsp.evalue, 1e-300)) * 5.0))

    score = (
        0.45 * identity_component
        + 0.35 * coverage_component
        + 0.20 * evalue_component
    )
    score = round(min(100.0, max(0.0, score)), 2)

    significant = hsp.evalue <= significance_evalue

    if not significant:
        band = "NOT_SIGNIFICANT"
    elif score >= 90:
        band = "VERY_STRONG_SEQUENCE_EVIDENCE"
    elif score >= 75:
        band = "STRONG_SEQUENCE_EVIDENCE"
    elif score >= 55:
        band = "MODERATE_SEQUENCE_EVIDENCE"
    else:
        band = "LIMITED_SEQUENCE_EVIDENCE"

    notes = [
        "Evidence profile summarizes sequence-level BLAST output.",
        "It is not a probability, experimental efficacy score or functional annotation.",
        f"Bands require E-value <= {significance_evalue:g} (configurable convention, "
        "not a universal threshold).",
    ]
    if not significant:
        notes.append(
            f"Best E-value {hsp.evalue:.3g} exceeds {significance_evalue:g}: a match of "
            "this quality is expected by chance in a database of this size."
        )

    return HitEvidence(
        hit_id=hit.hit_id,
        accession=hit.accession,
        identity_percent=round(hsp.identity_percent, 3),
        query_coverage_percent=round(coverage_value, 3),
        evalue=hsp.evalue,
        bit_score=hsp.bit_score,
        alignment_length=hsp.alignment_length,
        evidence_score=score,
        evidence_band=band,
        notes=notes,
    )


def build_evidence_profile(
    run: BlastRun,
    *,
    significance_evalue: float = DEFAULT_EVIDENCE_EVALUE_THRESHOLD,
) -> List[HitEvidence]:
    return [
        calculate_evidence_score(hit, significance_evalue=significance_evalue)
        for hit in rank_hits(run.hits)
    ]


def group_hits_by_taxonomy(hits: Sequence[BlastHit]) -> Dict[str, int]:
    """Group explicit or conservatively inferred organism labels."""
    groups: Dict[str, int] = {}
    for hit in hits:
        names = hit.scientific_names or ["Unknown / not provided"]
        for name in names:
            groups[name] = groups.get(name, 0) + 1
    return dict(sorted(groups.items(), key=lambda item: (-item[1], item[0])))

def alignment_text(
    hit: BlastHit,
    *,
    hsp_number: Optional[int] = None,
    width: int = 60,
) -> str:
    """Render an HSP in a compact academic alignment view."""
    if width < 20:
        raise BlastAnalysisError("Alignment width must be at least 20.")

    hsp = hit.best_hsp
    if hsp_number is not None:
        hsp = next((x for x in hit.hsps if x.number == hsp_number), None)

    if hsp is None:
        raise BlastAnalysisError("Requested HSP was not found.")

    query = hsp.query_sequence
    subject = hsp.subject_sequence
    middle = hsp.midline

    lines: List[str] = []
    for offset in range(0, len(query), width):
        q = query[offset:offset + width]
        s = subject[offset:offset + width]
        m = middle[offset:offset + width]

        q_start = hsp.query_start + offset
        s_start = hsp.subject_start + offset

        q_step = len(q.replace("-", ""))
        s_step = len(s.replace("-", ""))

        q_end = q_start + max(0, q_step - 1)
        s_end = s_start + max(0, s_step - 1)

        lines.append(f"Query    {q_start:>8}  {q}  {q_end:<8}")
        lines.append(f"                  {m}")
        lines.append(f"Subject  {s_start:>8}  {s}  {s_end:<8}")
        lines.append("")

    return "\r\n".join(lines).rstrip()
