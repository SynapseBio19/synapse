"""Annotation, evidence, consequence priority and feature overlap."""

from __future__ import annotations

import math
import re
from dataclasses import asdict
from typing import Any, Dict, List, Mapping, Optional, Sequence

from analyses._variant.models import (
    CONSEQUENCE_PRIORITY,
    Variant,
    VariantAnalysisError,
    VariantEvidence,
)


# ============================================================================
# Annotation
# ============================================================================

def annotate_variant(
    variant: Variant,
    *,
    gene: Optional[str] = None,
    transcript: Optional[str] = None,
    region: Optional[str] = None,
    consequence: Optional[str] = None,
    coding_change: Optional[str] = None,
    protein_change: Optional[str] = None,
    allele_frequency: Optional[float] = None,
    allele_count: Optional[int] = None,
    allele_number: Optional[int] = None,
    extra: Optional[Mapping[str, Any]] = None,
) -> Variant:
    """Attach supplied evidence/annotation to a variant."""
    if gene is not None:
        variant.gene = gene
    if transcript is not None:
        variant.transcript = transcript
    if region is not None:
        variant.region = region
    if consequence is not None:
        variant.consequence = consequence
    if coding_change is not None:
        variant.coding_change = coding_change
    if protein_change is not None:
        variant.protein_change = protein_change
    if allele_frequency is not None:
        if not 0 <= allele_frequency <= 1:
            raise VariantAnalysisError("Allele frequency must be between 0 and 1.")
        variant.allele_frequency = allele_frequency
    if allele_count is not None:
        variant.allele_count = allele_count
    if allele_number is not None:
        variant.allele_number = allele_number
    if extra:
        variant.annotation.update(dict(extra))
    return variant


def annotate_from_info(
    variants: Sequence[Variant],
    *,
    gene_keys: Sequence[str] = ("GENE", "SYMBOL"),
    transcript_keys: Sequence[str] = ("TRANSCRIPT", "TX"),
    consequence_keys: Sequence[str] = ("CSQ", "CONSEQUENCE"),
    frequency_keys: Sequence[str] = ("AF", "MAF"),
) -> List[Variant]:
    """Populate typed annotation fields from INFO without overwriting explicit data.

    Explicitly annotated ``Variant`` fields take precedence. INFO is used as the
    canonical fallback for VCF-derived records.
    """
    for variant in variants:
        if variant.gene is None:
            for key in gene_keys:
                if key in variant.info and variant.info[key] not in (None, ""):
                    value = variant.info[key]
                    variant.gene = value[0] if isinstance(value, list) else str(value)
                    break

        if variant.transcript is None:
            for key in transcript_keys:
                if key in variant.info and variant.info[key] not in (None, ""):
                    value = variant.info[key]
                    variant.transcript = value[0] if isinstance(value, list) else str(value)
                    break

        if variant.consequence is None:
            for key in consequence_keys:
                if key in variant.info and variant.info[key] not in (None, ""):
                    value = variant.info[key]
                    if isinstance(value, list):
                        value = value[0]
                    variant.consequence = str(value).split("|")[0]
                    break

        if variant.allele_frequency is None:
            for key in frequency_keys:
                if key not in variant.info:
                    continue
                value = variant.info[key]
                if isinstance(value, list):
                    value = value[0] if value else None
                try:
                    if value is not None:
                        frequency = float(value)
                        if 0.0 <= frequency <= 1.0:
                            variant.allele_frequency = frequency
                except (TypeError, ValueError):
                    pass
                break
    return list(variants)


# ============================================================================
# Canonical field synchronization
# ============================================================================

def synchronize_variant_fields(variants: Sequence[Variant]) -> List[Variant]:
    """Populate typed fields from INFO for every variant in a deterministic pass."""
    return annotate_from_info(list(variants))


# ============================================================================
# Consequence prioritization
# ============================================================================

def consequence_priority(consequence: Optional[str]) -> int:
    if consequence is None:
        return CONSEQUENCE_PRIORITY["unknown"]

    normalized = consequence.lower().replace("-", "_").replace(" ", "_")

    for key, score in CONSEQUENCE_PRIORITY.items():
        if key in normalized:
            return score
    return CONSEQUENCE_PRIORITY["unknown"]


def prioritize_variants(
    variants: Sequence[Variant],
) -> List[Variant]:
    """Rank variants by supplied consequence and evidence fields."""
    def sort_key(variant: Variant):
        quality = variant.quality if variant.quality is not None else -math.inf
        frequency = (
            variant.allele_frequency
            if variant.allele_frequency is not None
            else math.inf
        )
        return (
            -consequence_priority(variant.consequence),
            frequency,
            -quality,
            variant.chromosome,
            variant.position,
        )

    return sorted(variants, key=sort_key)


# ============================================================================
# Region / interval operations
# ============================================================================

def check_features_against_variants(
    features: Sequence[Dict[str, Any]],
    variants: Sequence[Variant],
    *,
    chromosome: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Variants inside primer or guide binding sites (audit P2-4).

    Each feature is a dict with "name", "kind", "start" and "end" (1-based
    inclusive on the sequence the variants are mapped to) and optionally
    "critical" (start, end) with a "critical_label", the part where a
    mismatch matters most: a primer's 3' end, a guide's seed and PAM. A
    variant overlaps when any base of its REF allele lies in the site.
    Overlaps in the critical part are HIGH, elsewhere in the site MODERATE.
    chromosome restricts the check to variants on that contig; the caller is
    responsible for the variants being on the analyzed sequence.
    """
    findings: List[Dict[str, Any]] = []
    for feature in features:
        start, end = int(feature["start"]), int(feature["end"])
        critical = feature.get("critical")
        for variant in variants:
            if chromosome is not None and variant.chromosome != chromosome:
                continue
            first = variant.position
            last = variant.position + max(len(variant.reference), 1) - 1
            if first > end or last < start:
                continue
            in_critical = bool(critical) and not (first > critical[1] or last < critical[0])
            findings.append({
                "feature": feature["name"],
                "kind": feature.get("kind", ""),
                "feature_start": start,
                "feature_end": end,
                "chromosome": variant.chromosome,
                "position": variant.position,
                "reference": variant.reference,
                "alternative": variant.alternative,
                "region": feature.get("critical_label", "critical region") if in_critical else "binding site",
                "severity": "HIGH" if in_critical else "MODERATE",
            })
    return findings


# The terminal dialogs of this module are in ui/variant_analysis_menu.py (P2-2).


# ============================================================================
# Mutation-module adapter
# ============================================================================

def annotate_from_mutation_analysis(
    variant: Variant,
    *,
    reference_sequence: str,
    alternative_sequence: str,
) -> Variant:
    """
    Optional bridge to mutation.py.

    It imports mutation.py lazily so variant_analysis.py remains independently
    usable. The adapter uses the sequence-level engine when available and
    attaches the most direct consequence information to this variant.
    """
    # Audit S15: the import named a top-level "mutation" module and the call
    # passed keyword arguments compare_sequences does not have, so this
    # function always failed.
    try:
        from analyses.mutation import compare_sequences, mutation_hgvs_p
    except ImportError as exc:
        raise VariantAnalysisError(
            "mutation.py could not be imported for sequence-derived annotation."
        ) from exc

    result = compare_sequences(
        reference_sequence,
        alternative_sequence,
        reference_name="Reference",
        variant_name="Variant",
    )

    if not result.mutations:
        return variant

    first = result.mutations[0]
    variant.variant_type = (
        "SNV" if first.mutation_type == "SNP"
        else first.mutation_type.upper()
    )
    variant.consequence = first.effect
    variant.coding_change = (
        f"{first.reference_codon}>{first.alternative_codon}"
        if first.reference_codon and first.alternative_codon
        else None
    )
    # The consequence class goes to `consequence`; `protein_change` holds the
    # predicted protein notation, e.g. p.(R12C), or nothing outside a codon.
    variant.protein_change = mutation_hgvs_p(first)
    variant.annotation["mutation_signature"] = {
        "position": first.position,
        "reference": first.reference,
        "alternative": first.alternative,
    }
    return variant


def attach_evidence(
    variant: Variant,
    evidence: VariantEvidence,
) -> Variant:
    variant.annotation.setdefault("evidence", []).append(asdict(evidence))
    return variant


def highest_priority_consequence(
    variant: Variant,
) -> Optional[str]:
    if not variant.consequence:
        return None
    consequences = re.split(r"[,&|;]", variant.consequence)
    consequences = [item.strip() for item in consequences if item.strip()]
    if not consequences:
        return None
    return max(consequences, key=consequence_priority)
