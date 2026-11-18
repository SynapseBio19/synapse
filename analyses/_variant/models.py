"""Constants, result dataclasses and the module error."""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple


# ============================================================================
# Constants
# ============================================================================

DNA_BASES = set("ACGTN")
VALID_FILTER_STATUS = {"PASS", "FAIL", "UNFILTERED", "UNKNOWN"}
VALID_VARIANT_TYPES = {
    "SNV",
    "MNV",
    "INSERTION",
    "DELETION",
    "INDEL",
    "COMPLEX",
    "SYMBOLIC",
    "UNKNOWN",
}

CONSEQUENCE_PRIORITY = {
    "frameshift": 100,
    "stop_gained": 95,
    "stop_lost": 90,
    "start_lost": 85,
    "splice_acceptor": 85,
    "splice_donor": 85,
    "missense": 70,
    "inframe_deletion": 55,
    "inframe_insertion": 55,
    "synonymous": 30,
    "utr": 15,
    "intronic": 10,
    "intergenic": 0,
    "unknown": -1,
}

GT_ALLELE_SEPARATOR_RE = re.compile(r"[\/|]")


# ============================================================================
# Exceptions
# ============================================================================

class VariantAnalysisError(ValueError):
    """Raised for invalid variant-analysis input."""


# ============================================================================
# Core dataclasses
# ============================================================================

@dataclass
class Variant:
    """Canonical internal representation of a genomic variant."""

    chromosome: str
    position: int
    reference: str
    alternative: str
    variant_id: Optional[str] = None
    quality: Optional[float] = None
    filter_status: str = "UNFILTERED"
    info: Dict[str, Any] = field(default_factory=dict)
    format_keys: List[str] = field(default_factory=list)
    sample_data: Dict[str, Dict[str, str]] = field(default_factory=dict)

    variant_type: str = "UNKNOWN"
    end_position: Optional[int] = None
    gene: Optional[str] = None
    transcript: Optional[str] = None
    region: Optional[str] = None
    consequence: Optional[str] = None
    protein_change: Optional[str] = None
    coding_change: Optional[str] = None
    allele_frequency: Optional[float] = None
    allele_count: Optional[int] = None
    allele_number: Optional[int] = None

    normalized: bool = False
    normalization_notes: List[str] = field(default_factory=list)
    annotation: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.chromosome = str(self.chromosome).strip()
        self.reference = str(self.reference).strip().upper()
        self.alternative = str(self.alternative).strip().upper()

        if not self.chromosome:
            raise VariantAnalysisError("Chromosome/contig cannot be empty.")
        if self.position < 1:
            raise VariantAnalysisError("Variant position must be >= 1.")
        if not self.reference:
            raise VariantAnalysisError("Reference allele cannot be empty.")
        if not self.alternative:
            raise VariantAnalysisError("Alternative allele cannot be empty.")

        self.filter_status = str(self.filter_status or "UNFILTERED").upper()

        if self.end_position is None:
            self.end_position = self.position + max(0, len(self.reference) - 1)

        if self.variant_type == "UNKNOWN":
            self.variant_type = classify_variant_type(
                self.reference, self.alternative
            )

    @property
    def key(self) -> Tuple[str, int, str, str]:
        return (
            self.chromosome,
            self.position,
            self.reference,
            self.alternative,
        )

    @property
    def allele_change(self) -> str:
        return f"{self.reference}>{self.alternative}"

    @property
    def is_snv(self) -> bool:
        return self.variant_type == "SNV"

    @property
    def is_indel(self) -> bool:
        return self.variant_type in {"INSERTION", "DELETION", "INDEL"}

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class VariantFilter:
    """Composable filtering criteria."""

    minimum_quality: Optional[float] = None
    maximum_quality: Optional[float] = None
    allowed_filter_statuses: Optional[Set[str]] = None
    variant_types: Optional[Set[str]] = None
    chromosomes: Optional[Set[str]] = None
    genes: Optional[Set[str]] = None
    consequences: Optional[Set[str]] = None
    minimum_allele_frequency: Optional[float] = None
    maximum_allele_frequency: Optional[float] = None
    minimum_allele_count: Optional[int] = None
    maximum_allele_count: Optional[int] = None
    require_pass: bool = False
    include_unknown_quality: bool = False
    include_unknown_frequency: bool = False


@dataclass
class VariantFilterResult:
    """Result of applying a VariantFilter."""

    passed: List[Variant] = field(default_factory=list)
    failed: List[Variant] = field(default_factory=list)
    reasons: Dict[Tuple[str, int, str, str], List[str]] = field(
        default_factory=dict
    )


@dataclass
class VariantStatistics:
    """Aggregate cohort/record statistics."""

    total_variants: int = 0
    unique_variants: int = 0
    snvs: int = 0
    mnvs: int = 0
    insertions: int = 0
    deletions: int = 0
    indels: int = 0
    complex: int = 0
    symbolic: int = 0
    unknown: int = 0
    passed_filters: int = 0
    failed_filters: int = 0
    unfiltered: int = 0
    missing_quality: int = 0
    mean_quality: Optional[float] = None
    median_quality: Optional[float] = None
    min_quality: Optional[float] = None
    max_quality: Optional[float] = None
    mean_allele_frequency: Optional[float] = None
    median_allele_frequency: Optional[float] = None
    rare_variants: int = 0
    genes_affected: int = 0
    chromosomes_affected: int = 0
    consequence_counts: Dict[str, int] = field(default_factory=dict)
    transition_count: int = 0
    transversion_count: int = 0
    transition_transversion_ratio: Optional[float] = None


@dataclass
class SampleVariantSummary:
    """Per-sample variant burden summary."""

    sample: str
    total_variant_records: int = 0
    called_variants: int = 0
    heterozygous_calls: int = 0
    homozygous_alt_calls: int = 0
    homozygous_ref_calls: int = 0
    missing_calls: int = 0
    allele_count: int = 0


@dataclass
class CohortVariantSummary:
    """Cohort-wide variant summary."""

    sample_count: int = 0
    variant_count: int = 0
    unique_variant_count: int = 0
    shared_variant_count: int = 0
    private_variant_count: int = 0
    variant_prevalence: Dict[Tuple[str, int, str, str], int] = field(
        default_factory=dict
    )
    sample_summaries: List[SampleVariantSummary] = field(default_factory=list)


@dataclass
class VariantComparison:
    """Comparison of two variant collections."""

    only_in_left: List[Variant] = field(default_factory=list)
    only_in_right: List[Variant] = field(default_factory=list)
    shared: List[Variant] = field(default_factory=list)
    left_count: int = 0
    right_count: int = 0
    shared_count: int = 0
    jaccard_similarity: float = 0.0


@dataclass
class VariantAnalysisResult:
    """Top-level analysis result."""

    name: str
    variants: List[Variant] = field(default_factory=list)
    statistics: VariantStatistics = field(default_factory=VariantStatistics)
    filter_result: Optional[VariantFilterResult] = None
    cohort_summary: Optional[CohortVariantSummary] = None
    notes: List[str] = field(default_factory=list)

    @property
    def count(self) -> int:
        return len(self.variants)

    def to_dict(self) -> dict:
        payload = asdict(self)
        if self.filter_result is not None:
            payload["filter_result"] = {
                "passed": [v.to_dict() for v in self.filter_result.passed],
                "failed": [v.to_dict() for v in self.filter_result.failed],
                "reasons": {
                    str(k): v
                    for k, v in self.filter_result.reasons.items()
                },
            }
        return payload

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)


def classify_variant_type(reference: str, alternative: str) -> str:
    """Classify a single biallelic allele pair."""
    ref = str(reference).upper()
    alt = str(alternative).upper()

    if alt.startswith("<") and alt.endswith(">"):
        return "SYMBOLIC"
    if alt == "*":
        return "DELETION"
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(ref) < len(alt):
        return "INSERTION"
    if len(ref) > len(alt):
        return "DELETION"
    return "COMPLEX"


# ============================================================================
# Clinical/evidence-safe evidence helpers
# ============================================================================

@dataclass
class VariantEvidence:
    """Evidence container; no interpretation is inferred from missing data."""

    source: str
    field: str
    value: Any
    evidence_level: Optional[str] = None
    accession: Optional[str] = None
    notes: Optional[str] = None
