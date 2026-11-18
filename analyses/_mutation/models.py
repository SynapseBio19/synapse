"""Result dataclasses and errors of the mutation analysis."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional


# ============================================================================
# EXCEPTIONS
# ============================================================================

class SequencePromptCancelled(Exception):
    """The user left a sequence prompt with [0]; the caller returns to its menu."""


class MutationAnalysisError(ValueError):
    """Raised when mutation-analysis input or state is invalid."""


# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass(frozen=True)
class AlignmentBlock:
    operation: str
    reference_start: int
    reference_end: int
    variant_start: int
    variant_end: int
    reference_sequence: str
    variant_sequence: str


@dataclass
class Mutation:
    mutation_id: int
    mutation_type: str
    position: int
    end_position: int
    reference: str
    alternative: str
    length_change: int
    notation: str

    transition: Optional[bool] = None

    region: str = "Unknown"

    codon_index: Optional[int] = None
    codon_position: Optional[int] = None

    reference_codon: Optional[str] = None
    alternative_codon: Optional[str] = None

    reference_amino_acid: Optional[str] = None
    alternative_amino_acid: Optional[str] = None

    effect: Optional[str] = None
    consequence: Optional[str] = None

    frameshift: bool = False
    premature_stop: bool = False

    notes: List[str] = field(default_factory=list)


@dataclass
class MutationStatistics:
    total_mutations: int = 0
    substitutions: int = 0
    snps: int = 0
    insertions: int = 0
    deletions: int = 0
    complex_changes: int = 0

    transitions: int = 0
    transversions: int = 0

    transition_transversion_ratio: Optional[float] = None

    inserted_bases: int = 0
    deleted_bases: int = 0

    net_length_change: int = 0

    synonymous: int = 0
    missense: int = 0
    nonsense: int = 0
    start_loss: int = 0
    stop_loss: int = 0

    frameshifts: int = 0
    in_frame_indels: int = 0

    unknown_effects: int = 0
    affected_codons: int = 0

    mutation_density_per_kb: float = 0.0


@dataclass
class CodingAnalysis:
    reference_protein: str
    alternative_protein: str

    reference_coding_length: int
    alternative_coding_length: int

    protein_length_change: int

    affected_codons: List[int] = field(default_factory=list)

    synonymous_count: int = 0
    missense_count: int = 0
    nonsense_count: int = 0
    start_loss_count: int = 0
    stop_loss_count: int = 0

    frameshift: bool = False
    frameshift_position: Optional[int] = None

    premature_stop_position: Optional[int] = None

    notes: List[str] = field(default_factory=list)


@dataclass
class MutationAnalysisResult:
    reference_name: str
    variant_name: str

    reference_length: int
    variant_length: int

    mutations: List[Mutation] = field(default_factory=list)

    statistics: MutationStatistics = field(
        default_factory=MutationStatistics
    )

    coding_analysis: Optional[CodingAnalysis] = None

    aligned_reference: str = ""
    aligned_variant: str = ""

    alignment_score: Optional[int] = None

    # Alignment Identity:
    # Exact nucleotide matches among aligned, non-gap positions.
    identity_percent: float = 0.0

    # Sequence Similarity:
    # Edit-distance similarity across the complete ungapped sequences.
    sequence_similarity_percent: float = 0.0

    notes: List[str] = field(default_factory=list)

    @property
    def mutation_count(self) -> int:
        return len(self.mutations)

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(
            self.to_dict(),
            indent=indent,
            sort_keys=False,
        )


@dataclass(frozen=True)
class SequenceCoordinate:
    zero_based: int
    one_based: int

    @classmethod
    def from_one_based(
        cls,
        position: int,
    ) -> "SequenceCoordinate":

        if position < 1:
            raise MutationAnalysisError(
                "One-based coordinates must be >= 1."
            )

        return cls(
            zero_based=position - 1,
            one_based=position,
        )


@dataclass
class SequenceQuality:
    length: int
    gc_percent: float
    n_count: int
    ambiguous_count: int

    valid: bool

    longest_homopolymer: int
    longest_homopolymer_base: Optional[str]

    low_complexity: bool
    complexity_percent: float

    warnings: List[str] = field(default_factory=list)


@dataclass
class MutationContext:
    position: int

    upstream: str
    downstream: str

    reference_context: str
    variant_context: str

    kmer_size: int

    trinucleotide_context: Optional[str] = None


@dataclass
class MutationSpectrum:
    substitutions_12: Dict[str, int] = field(
        default_factory=dict
    )

    substitutions_6: Dict[str, int] = field(
        default_factory=dict
    )

    transitions: int = 0
    transversions: int = 0


@dataclass
class MutationHotspot:
    start: int
    end: int

    mutation_count: int
    density_per_kb: float


@dataclass
class RestrictionSiteImpact:
    enzyme: str
    recognition_sequence: str

    lost_sites: List[int] = field(default_factory=list)
    gained_sites: List[int] = field(default_factory=list)

    unchanged_sites: int = 0


@dataclass
class MutationImpact:
    restriction_sites_lost: int = 0
    restriction_sites_gained: int = 0

    primer_overlap: bool = False
    primer_overlap_ids: List[str] = field(default_factory=list)

    crispr_target_overlap: bool = False
    crispr_overlap_ids: List[str] = field(default_factory=list)

    orf_disrupted: bool = False

    rna_region: Optional[str] = None

    notes: List[str] = field(default_factory=list)


@dataclass
class AdvancedMutationAnalysis:
    reference_quality: Optional[SequenceQuality] = None
    variant_quality: Optional[SequenceQuality] = None

    spectra: Optional[MutationSpectrum] = None

    contexts: List[MutationContext] = field(
        default_factory=list
    )

    hotspots: List[MutationHotspot] = field(
        default_factory=list
    )

    impacts: MutationImpact = field(
        default_factory=MutationImpact
    )

    ambiguity_warnings: List[str] = field(
        default_factory=list
    )

    normalization_notes: List[str] = field(
        default_factory=list
    )
